#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_alloc.h>

void
gk_species_damping_write_disabled(gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm, int frame)
{
}

void
gk_species_damping_write_enabled(gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  // DG metadata for damping rate.
  struct gkyl_msgpack_map_elem mpe_drate[] = {
    { .key = "poly_order", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 0 },
    { .key = "basis_type", .elem_type = GKYL_MP_STRING, .cval = "serendipity" },
    { .key = "Description", .elem_type = GKYL_MP_STRING, .cval = "Rate of the damping term." },
    { .key = "time", .elem_type = GKYL_MP_DOUBLE, .dval = tm },
    { .key = "frame", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = frame },
  };

  int mpe_drate_len = sizeof(mpe_drate)/sizeof(mpe_drate[0]);
  // Package metadata.
  int io_meta_len[] = {gks->io_meta_basic_len, mpe_drate_len, app->gk_geom->io_meta_basic_len};
  const struct gkyl_msgpack_map_elem* io_meta[] = {gks->io_meta_basic, mpe_drate, app->gk_geom->io_meta_basic};
  struct gkyl_msgpack_data *mt = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);

  // Write out the damping rate.
  const char *fmt = "%s-%s_damping_rate_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, frame);
  char fileNm[sz + 1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, frame);

  // Copy data from device to host before writing it out.
  if (app->use_gpu)
    gkyl_array_copy(gks->damping.rate_host, gks->damping.rate);

  gkyl_comm_array_write(gks->comm, &gks->grid, &gks->local, mt, gks->damping.rate_host, fileNm);
  app->stat.n_io += 1;

  gkyl_msgpack_data_release(mt);
  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_species_damping_write_init_only(gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm, int frame)
{
  gk_species_damping_write_enabled(app, gks, tm, frame);
  gks->damping.write_func = gk_species_damping_write_disabled;
}

static void
proj_on_basis_c2p_phase_func(const double *xcomp, double *xphys, void *ctx)
{
  struct gk_proj_on_basis_c2p_func_ctx *c2p_ctx = ctx;
  int cdim = c2p_ctx->cdim; // Assumes update range is a phase range.
  gkyl_velocity_map_eval_c2p(c2p_ctx->vel_map, &xcomp[cdim], &xphys[cdim]);
}

void
gk_species_damping_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_damping *damp)
{
  damp->type = gks->info.damping.type;
  damp->evolve = false; // Whether the rate is time dependent.

  int num_quad = gks->info.damping.num_quad? gks->info.damping.num_quad : 1; // Default is a p=0 mask.
  assert(num_quad == 1); // MF 2025/06/11: Limited to this for now.

  // Default function pointers.
  damp->write_func = gk_species_damping_write_disabled;

  if (damp->type) {
    // Allocate rate array.
    damp->rate = mkarr(app->use_gpu, num_quad == 1? 1 : gks->basis.num_basis, gks->local_ext.volume);
    damp->rate_host = damp->rate;
    if (app->use_gpu)
      damp->rate_host = mkarr(false, damp->rate->ncomp, damp->rate->size);

    if (damp->type == GKYL_GK_DAMPING_USER_INPUT) {
      struct gk_proj_on_basis_c2p_func_ctx proj_on_basis_c2p_ctx; // c2p function context.
      proj_on_basis_c2p_ctx.cdim = app->cdim;
      proj_on_basis_c2p_ctx.vdim = gks->local_vel.ndim;
      proj_on_basis_c2p_ctx.vel_map = gks->vel_map;
      gkyl_proj_on_basis *projup = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
          .grid = &gks->grid,
          .basis = &gks->basis,
          .num_quad = num_quad,
          .num_ret_vals = 1,
          .eval = gks->info.damping.rate_profile,
          .ctx = gks->info.damping.rate_profile_ctx,
          .c2p_func = proj_on_basis_c2p_phase_func,
          .c2p_func_ctx = &proj_on_basis_c2p_ctx,
        }
      );
      gkyl_proj_on_basis_advance(projup, 0.0, &gks->local, damp->rate_host);
      gkyl_proj_on_basis_release(projup);
      gkyl_array_copy(damp->rate, damp->rate_host);

      if (num_quad == 1)
        gkyl_array_scale_range(damp->rate, 1.0/pow(sqrt(2.0),gks->grid.ndim), &gks->local);
    }

    // Set function pointers chosen at runtime.
    if (damp->evolve) {
      damp->write_func = gk_species_damping_write_enabled;
    }
    else {
      damp->write_func = gk_species_damping_write_init_only;
    }
  }
}

void
gk_species_damping_advance(gkyl_gyrokinetic_app *app, const struct gk_species *gks, struct gk_damping *damp,
  const struct gkyl_array *phi, const struct gkyl_array *fin, struct gkyl_array *f_buffer,
  struct gkyl_array *rhs, struct gkyl_array *cflrate)
{
  if (damp->type) {
    struct timespec wst = gkyl_wall_clock();
    if (damp->type == GKYL_GK_DAMPING_USER_INPUT) {
      gkyl_array_set(f_buffer, 1.0, fin);
      gkyl_array_scale_by_cell(f_buffer, damp->rate);
      gkyl_array_accumulate(rhs, -1.0, f_buffer);
    }

    // Add the frequency to the CFL frequency.
    gkyl_array_accumulate(cflrate, 1.0, damp->rate);

    app->stat.species_damp_tm += gkyl_time_diff_now_sec(wst);
  }
}

void
gk_species_damping_write(gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm, int frame)
{
  gks->damping.write_func(app, gks, tm, frame);
}

void
gk_species_damping_release(const struct gkyl_gyrokinetic_app *app, const struct gk_damping *damp)
{
  if (damp->type) {
    gkyl_array_release(damp->rate);
    if (app->use_gpu) {
      gkyl_array_release(damp->rate_host);
    }

    if (damp->type == GKYL_GK_DAMPING_USER_INPUT) {
      // Nothing to release.
    }
  }
}
