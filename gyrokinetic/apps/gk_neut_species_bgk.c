#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_const.h>

static void
gkbgk_moms_disabled(gkyl_gyrokinetic_app *app, const struct gk_neut_species *gkns,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin)
{
  // Empty method.
}

static void
gkbgk_moms_enabled(gkyl_gyrokinetic_app *app, const struct gk_neut_species *gkns,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin)
{
  struct timespec wst = gkyl_wall_clock();

  // Compute Maxwellian moments (n, u_par, T/m).
  gk_neut_species_moment_calc(&gkns->lte.moms, gkns->local, app->local, fin);
  gkyl_dg_div_op_range(gkns->lte.moms.mem_geo, &app->basis, 0, gkns->lte.moms.marr,
    0, gkns->lte.moms.marr, 0, app->gk_geom->geo_int.jacobgeo, &app->local);
  
  // Calculate nu_ss.
  bgk->self_nu_func_neut(app, gkns, bgk);

  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);    
}

static void
gkbgk_self_nu_calc_constNu(gkyl_gyrokinetic_app *app, const struct gk_neut_species *gkns,
  struct gk_bgk_collisions *bgk)
{
  // Empty method.
}

static void
gkbgk_rhs_disabled(gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin, double dt, struct gkyl_array *rhs)
{
  // Empty method.
}

static void
gkbgk_rhs_enabled(gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin, double dt, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();
    
  // Compute the self-collisions Maxwellian.
  struct gk_lte *lte = &gkns->lte;
  gk_neut_species_lte_from_moms(app, gkns, lte, lte->moms.marr);

  // Multiply the Maxwellian by self-species collision frequency.
  gkyl_dg_mul_conf_phase_op_range(&app->basis, &gkns->basis, bgk->nu_fmax, 
    bgk->self_nu, lte->f_lte, &app->local, &gkns->local);

  // Multiply the Maxwellian by the configuration-space Jacobian.
  gkyl_dg_mul_conf_phase_op_range(&app->basis, &gkns->basis, bgk->nu_fmax,
    app->gk_geom->geo_int.jacobgeo, bgk->nu_fmax, &app->local, &gkns->local);

  // Apply BGK collisions.
  gkyl_bgk_collisions_advance(bgk->up_bgk, &app->local, &gkns->local, 
    bgk->nu_sum, bgk->nu_fmax, fin, bgk->implicit_step, dt, rhs, gkns->cflrate);
  
  app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);
}

static void
gkbgk_write_mom_disabled(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame)
{
  // Empty method.
}

static void
gkbgk_write_mom_enabled(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame)
{
  struct timespec wtm = gkyl_wall_clock();
  // Package metadata.
  gkyl_msgpack_map_elem_set_double(gkns->io_meta_conf_len, gkns->io_meta_conf, "time", tm);
  gkyl_msgpack_map_elem_set_uint(gkns->io_meta_conf_len, gkns->io_meta_conf, "frame", frame);
  struct gkyl_msgpack_map_elem desc[] = {
    { .key = "Description", .elem_type = GKYL_MP_STRING, .cval = "Sum of collision frequencies for BGK collisions." }
  };
  int io_meta_len[] = {gkns->io_meta_conf_len, app->gk_geom->io_meta_basic_len, 1};
  const struct gkyl_msgpack_map_elem* io_meta[] = {gkns->io_meta_conf, app->gk_geom->io_meta_basic, desc};
  struct gkyl_msgpack_data *mt = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);

  // Write out nu_sum.
  const char *fmt = "%s-%s_nu_sum_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, gkns->info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, gkns->info.name, frame);
  
  // Copy data from device to host before writing it out.
  if (app->use_gpu) {
    gkyl_array_copy(gkns->bgk.nu_sum_host, gkns->bgk.nu_sum);
  }
  
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, gkns->bgk.nu_sum_host, fileNm);
  app->stat.n_diag_io += 2;

  gkyl_msgpack_data_release(mt); 
  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wtm);
}

void 
gk_neut_species_bgk_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns, struct gk_bgk_collisions *bgk)
{
  bgk->collision_id = gkns->info.collisions.collision_id;
  bgk->write_diagnostics = gkns->info.collisions.write_diagnostics;

  // Empty methods.
  bgk->moms_func_neut = gkbgk_moms_disabled;
  bgk->rhs_func_neut = gkbgk_rhs_disabled;
  bgk->moms_func_implicit_neut = gkbgk_moms_disabled;
  bgk->rhs_func_implicit_neut = gkbgk_rhs_disabled;
  bgk->write_mom_func_neut = gkbgk_write_mom_disabled;

  if (bgk->collision_id == GKYL_BGK_COLLISIONS) {
    bgk->implicit_step = gkns->info.collisions.is_implicit;
    bgk->num_cross_collisions = gkns->info.collisions.num_cross_collisions;
    assert(bgk->num_cross_collisions == 0); // NYI.
    
    int cdim = app->cdim, vdim = gkns->info.vdim;
  
    // Allocate self-species collision frequency and sum of collision frequencies.
    bgk->self_nu = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    bgk->nu_sum = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  
    double nu_frac = gkns->info.collisions.nu_frac ? gkns->info.collisions.nu_frac : 1.0;
  
    if (gkns->info.collisions.self_nu) {
      // Project user's self-species collision frequency.
      bgk->norm_nu_self = false;
  
      struct gkyl_array *self_nu_ho = mkarr(false, app->basis.num_basis, app->local_ext.volume);
      gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&app->grid, &app->basis,
        app->poly_order+1, 1, gkns->info.collisions.self_nu, gkns->info.collisions.self_nu_ctx);
      gkyl_proj_on_basis_advance(proj, 0.0, &app->local, self_nu_ho);
      gkyl_proj_on_basis_release(proj);
      gkyl_array_copy(bgk->self_nu, self_nu_ho);
      gkyl_array_release(self_nu_ho);
  
      gkyl_array_scale(bgk->self_nu, nu_frac);
      gkyl_array_set(bgk->nu_sum, 1.0, bgk->self_nu);
  
      // Set pointers to functions chosen at runtime.
      bgk->self_nu_func_neut = gkbgk_self_nu_calc_constNu;
    }
    else {
      // Self-collision frequency computed in time.
      assert(false); // NYI
    }
  
    // Host-side copy for I/O.
    if (bgk->write_diagnostics) {
      if (app->use_gpu) {
        bgk->nu_sum_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
      }
      else {
        bgk->nu_sum_host = bgk->nu_sum;
      }
    }
  
    // Collision frequency times Maxwellian.
    bgk->nu_fmax = mkarr(app->use_gpu, gkns->basis.num_basis, gkns->local_ext.volume);

    // BGK updater.
    bgk->up_bgk = gkyl_bgk_collisions_new(&app->basis, &gkns->basis, app->use_gpu);

    // Methods chosen at runtime.
    if (bgk->implicit_step) {
      bgk->moms_func_neut = gkbgk_moms_enabled;
      bgk->rhs_func_neut = gkbgk_rhs_enabled;
      bgk->moms_func_implicit_neut = gkbgk_moms_disabled;
      bgk->rhs_func_implicit_neut = gkbgk_rhs_disabled;
    }
    else {
      bgk->moms_func_neut = gkbgk_moms_disabled;
      bgk->rhs_func_neut = gkbgk_rhs_disabled;
      bgk->moms_func_implicit_neut = gkbgk_moms_enabled;
      bgk->rhs_func_implicit_neut = gkbgk_rhs_enabled;
    }
    if (bgk->write_diagnostics)
      bgk->write_mom_func_neut = gkbgk_write_mom_enabled;
  }
}

void 
gk_neut_species_bgk_cross_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns, struct gk_bgk_collisions *bgk)
{  
}

void
gk_neut_species_bgk_moms(gkyl_gyrokinetic_app *app, const struct gk_neut_species *gkns,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin)
{
  bgk->moms_func_neut(app, gkns, bgk, fin);
}

void
gk_neut_species_bgk_moms_implicit(gkyl_gyrokinetic_app *app, const struct gk_neut_species *gkns,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin)
{
  bgk->moms_func_implicit_neut(app, gkns, bgk, fin);
}

void
gk_neut_species_bgk_cross_moms(gkyl_gyrokinetic_app *app, const struct gk_neut_species *gkns,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin)
{
  // Compute this in bgk_rhs
}

void
gk_neut_species_bgk_cross_moms_implicit(gkyl_gyrokinetic_app *app, const struct gk_neut_species *gkns,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin)
{
  // Compute this in bgk_rhs
}

void
gk_neut_species_bgk_rhs(gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  bgk->rhs_func_neut(app, gkns, bgk, fin, 0.0, rhs);
}

void
gk_neut_species_bgk_rhs_implicit(gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin, double dt, struct gkyl_array *rhs)
{
  bgk->rhs_func_implicit_neut(app, gkns, bgk, fin, dt, rhs);
}

void
gk_neut_species_bgk_write_mom(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame)
{
  gkns->bgk.write_mom_func_neut(app, gkns, tm, frame);
}

void 
gk_neut_species_bgk_release(const struct gkyl_gyrokinetic_app *app, const struct gk_bgk_collisions *bgk)
{
  if (bgk->collision_id == GKYL_BGK_COLLISIONS) {
    gkyl_array_release(bgk->nu_fmax);
    gkyl_bgk_collisions_release(bgk->up_bgk);

    if (bgk->write_diagnostics) {
      if (app->use_gpu) {
        gkyl_array_release(bgk->nu_sum_host);
      }
    }

    gkyl_array_release(bgk->nu_sum);
    gkyl_array_release(bgk->self_nu);
  }
}
