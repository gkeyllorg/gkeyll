#include <assert.h>
#include <gkyl_alloc.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_loss_cone_mask_gyrokinetic.h>

static enum gkyl_gk_trapped_passing_orbit_type
trapped_passing_orbit_from_species_bc(enum gkyl_gyrokinetic_bc_type bc)
{
  switch (bc) {
    case GKYL_BC_GK_SPECIES_SHEATH:
      return GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_SHEATH;
    case GKYL_BC_GK_SPECIES_REFLECT:
    case GKYL_BC_GK_SPECIES_ZERO_FLUX:
    case GKYL_BC_GK_SPECIES_PERIODIC:
    case GKYL_BC_GK_SPECIES_TWISTSHIFT:
      return GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_WALL;
    default:
      return GKYL_GK_TRAP_PASS_ORBIT_PASSING;
  }
}

static void
gk_species_fdot_multiplier_write_disabled(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_fdot_multiplier *fdot_mult, double tm, int frame)
{
}

static void
gk_species_fdot_multiplier_write_enabled(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_fdot_multiplier *fdot_mult, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  // DG metadata for multiplier.
  struct gkyl_msgpack_map_elem mpe_mult[] = {
    { .key = "poly_order", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 0 },
    { .key = "basis_type", .elem_type = GKYL_MP_STRING, .cval = "serendipity" },
    { .key = "Description", .elem_type = GKYL_MP_STRING, .cval = "Function multiplying the distribution time derivative." },
    { .key = "time", .elem_type = GKYL_MP_DOUBLE, .dval = tm },
    { .key = "frame", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = frame },
  };
  int mpe_mult_len = sizeof(mpe_mult) / sizeof(mpe_mult[0]);
  // Package metadata.
  int io_meta_len[] = { gks->io_meta_basic_len, mpe_mult_len, app->gk_geom->io_meta_basic_len };
  const struct gkyl_msgpack_map_elem *io_meta[] = { gks->io_meta_basic, mpe_mult,
                                                    app->gk_geom->io_meta_basic };
  struct gkyl_msgpack_data *mt = gkyl_msgpack_create_union(sizeof(io_meta_len) / sizeof(int),
    io_meta_len, io_meta);

  // Write out the combined multiplier.
  const char *fmt = "%s-%s_fdot_multiplier_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, frame);
  char fileNm[sz + 1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, frame);

  // Copy data from device to host before writing it out.
  if (app->use_gpu) {
    gkyl_array_copy(fdot_mult->multiplier_host, fdot_mult->multiplier);
  }

  gkyl_comm_array_write(gks->comm, &gks->grid, &gks->local, mt, fdot_mult->multiplier_host,
    fileNm);
  app->stat.n_io += 1;

  gkyl_msgpack_data_release(mt);
  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_fdot_multiplier_advance_times_cfl_disabled(gkyl_gyrokinetic_app *app,
  struct gk_species *gks, struct gk_fdot_multiplier *fdot_mult,
  const struct gkyl_array *phi, const struct gkyl_array *f, struct gkyl_array *cflrate)
{
}

static void
gk_fdot_multiplier_advance_times_cfl_enabled(gkyl_gyrokinetic_app *app,
  struct gk_species *gks, struct gk_fdot_multiplier *fdot_mult,
  const struct gkyl_array *phi, const struct gkyl_array *f, struct gkyl_array *cflrate)
{
  struct timespec wst = gkyl_wall_clock();

  // Reset combined to 1; each component multiplies its contribution in.
  gkyl_array_clear(fdot_mult->multiplier, 1.0);
  for (int i = 0; i < fdot_mult->num_multipliers; ++i) {
    struct gk_fdot_multiplier_comp *fdmul = &fdot_mult->comp[i];
    fdmul->advance_func(app, gks, fdmul, phi, f, cflrate, fdot_mult->multiplier);
  }
  gkyl_array_scale_by_cell(cflrate, fdot_mult->multiplier);

  app->stat.species_fdot_mult_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_fdot_multiplier_advance_times_rate_disabled(gkyl_gyrokinetic_app *app,
  struct gk_species *gks, struct gk_fdot_multiplier *fdot_mult,
  const struct gkyl_array *phi, const struct gkyl_array *f, struct gkyl_array *rhs)
{
}

static void
gk_fdot_multiplier_advance_times_rate_enabled(gkyl_gyrokinetic_app *app,
  struct gk_species *gks, struct gk_fdot_multiplier *fdot_mult,
  const struct gkyl_array *phi, const struct gkyl_array *f, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();

  gkyl_array_scale_by_cell(rhs, fdot_mult->multiplier);

  app->stat.species_fdot_mult_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_species_fdot_multiplier_advance_disabled(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_fdot_multiplier_comp *fdmul, const struct gkyl_array *phi, const struct gkyl_array *f,
  const struct gkyl_array *cflrate, struct gkyl_array *combined_multiplier)
{
}

// Advance: multiply combined_multiplier by this component's precomputed buffer array.
static void
gk_species_fdot_multiplier_advance_mult(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_fdot_multiplier_comp *fdmul, const struct gkyl_array *phi, const struct gkyl_array *f,
  const struct gkyl_array *cflrate, struct gkyl_array *combined_multiplier)
{
  gkyl_array_scale_by_cell(combined_multiplier, fdmul->buffer);
}

// Advance: apply constant scale to combined_multiplier.
static void
gk_species_fdot_multiplier_advance_const(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_fdot_multiplier_comp *fdmul, const struct gkyl_array *phi, const struct gkyl_array *f,
  const struct gkyl_array *cflrate, struct gkyl_array *combined_multiplier)
{
  gkyl_array_scale(combined_multiplier, fdmul->time_dilation_scale_const);
}

static void
gk_species_fdot_multiplier_advance_loss_cone_mult(gkyl_gyrokinetic_app *app,
  const struct gk_species *gks, struct gk_fdot_multiplier_comp *fdmul, const struct gkyl_array *phi,
  const struct gkyl_array *f, const struct gkyl_array *cflrate,
  struct gkyl_array *combined_multiplier)
{
  gkyl_comm_array_allgather(app->comm, &app->local, &app->global, phi, fdmul->phi_global);
  int zdim = app->cdim - 1;
  if (gks->lower_bc[zdim].type == GKYL_BC_GK_SPECIES_SHEATH) {
    gkyl_comm_array_allgather(app->comm, &app->local, &app->global,
      gks->phi_wall_lo.phi, fdmul->phi_wall_lo_global);
  }
  if (gks->upper_bc[zdim].type == GKYL_BC_GK_SPECIES_SHEATH) {
    gkyl_comm_array_allgather(app->comm, &app->local, &app->global,
      gks->phi_wall_up.phi, fdmul->phi_wall_up_global);
  }
  gkyl_loss_cone_mask_gyrokinetic_advance(fdmul->lcm_proj_op, &gks->local, &app->global,
    fdmul->bmag_global, fdmul->phi_global, fdmul->phi_wall_lo_global,
    fdmul->phi_wall_up_global, fdmul->buffer);
  gkyl_array_scale(fdmul->buffer, fdmul->time_dilation_scale_const);
  gkyl_array_scale_by_cell(combined_multiplier, fdmul->buffer);
}

// Compute cell-wise clamp factor min(1, omega_max/cflrate)*scale into combined_multiplier.
// Hopkins, P. F., & Most, E. R. (2025). Time-Dilation Methods for Extreme Multiscale Timestepping Problems. arXiv:2510.09756.
static double
compute_global_array_max(const gkyl_gyrokinetic_app *app,
  const struct gk_fdot_multiplier_comp *fdmul, const struct gkyl_array *array)
{
  double local_max;
  if (app->use_gpu) {
#ifdef GKYL_HAVE_CUDA
    gkyl_array_reduce(fdmul->omega_max_local_cu, array, GKYL_MAX);
    gkyl_cu_memcpy(&local_max, fdmul->omega_max_local_cu, sizeof(double), GKYL_CU_MEMCPY_D2H);
#endif
  }
  else {
    gkyl_array_reduce(&local_max, array, GKYL_MAX);
  }
  double global_max = DBL_MAX;
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &local_max, &global_max);
  return global_max;
}

static void
clamp_cflrate_by_omega_max(const struct gk_fdot_multiplier_comp *fdmul, const double omega_max,
  const struct gkyl_array *cflrate, struct gkyl_array *combined_multiplier)
{
  // buffer = min(1, omega_max/cflrate) * time_dilation_scale_const
  gkyl_array_invert_by_cell(fdmul->buffer, cflrate);
  gkyl_array_scale(fdmul->buffer, omega_max);
  gkyl_array_min_by_cell(fdmul->buffer, fdmul->buffer, 1.0);
  gkyl_array_scale(fdmul->buffer, fdmul->time_dilation_scale_const);
  gkyl_array_scale_by_cell(combined_multiplier, fdmul->buffer);
}

static void
gk_species_fdot_multiplier_advance_time_dilation_cfl_dt_omegaH(gkyl_gyrokinetic_app *app,
  const struct gk_species *gks, struct gk_fdot_multiplier_comp *fdmul, const struct gkyl_array *phi,
  const struct gkyl_array *f, const struct gkyl_array *cflrate,
  struct gkyl_array *combined_multiplier)
{
  double omega_max = 1.0 / gks->dt_omegaH;
  clamp_cflrate_by_omega_max(fdmul, omega_max, cflrate, combined_multiplier);
}

static void
gk_species_fdot_multiplier_advance_time_dilation_cfl_dt_user_specified(gkyl_gyrokinetic_app *app,
  const struct gk_species *gks, struct gk_fdot_multiplier_comp *fdmul, const struct gkyl_array *phi,
  const struct gkyl_array *f, const struct gkyl_array *cflrate,
  struct gkyl_array *combined_multiplier)
{
  double omega_max = 1.0 / fdmul->cfl_dt_min_value;
  clamp_cflrate_by_omega_max(fdmul, omega_max, cflrate, combined_multiplier);
}

static void
gk_species_fdot_multiplier_advance_time_dilation_cfl_dt_set_by_species(gkyl_gyrokinetic_app *app,
  const struct gk_species *gks, struct gk_fdot_multiplier_comp *fdmul, const struct gkyl_array *phi,
  const struct gkyl_array *f, const struct gkyl_array *cflrate,
  struct gkyl_array *combined_multiplier)
{
  double omega_max = 1.0 / fdmul->species_dt_is_set_from->dt_cfl_global_ho;
  clamp_cflrate_by_omega_max(fdmul, omega_max, cflrate, combined_multiplier);
}

static void
gk_species_fdot_multiplier_advance_time_dilation_cfl_factor_user_specified(
  gkyl_gyrokinetic_app *app, const struct gk_species *gks, struct gk_fdot_multiplier_comp *fdmul,
  const struct gkyl_array *phi, const struct gkyl_array *f, const struct gkyl_array *cflrate,
  struct gkyl_array *combined_multiplier)
{
  double omega_max = compute_global_array_max(app, fdmul, cflrate);
  omega_max = fdmul->cfl_factor_times_omega_max * omega_max;
  clamp_cflrate_by_omega_max(fdmul, omega_max, cflrate, combined_multiplier);
}

static void
gk_species_fdot_multiplier_advance_time_dilation_cfl_f_frac_global(gkyl_gyrokinetic_app *app,
  const struct gk_species *gks, struct gk_fdot_multiplier_comp *fdmul, const struct gkyl_array *phi,
  const struct gkyl_array *f, const struct gkyl_array *cflrate,
  struct gkyl_array *combined_multiplier)
{
  gkyl_array_reduce(fdmul->local_max_f, f, GKYL_MAX);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, app->basis.num_basis,
    fdmul->local_max_f, fdmul->global_max_f);
  gkyl_dg_array_mask_advance_threshold(fdmul->cfl_mask, fdmul->global_max_f[0]);

  gkyl_dg_array_mask_advance(fdmul->cfl_mask, f);
  gkyl_dg_array_mask_scale_by_cell(fdmul->cfl_mask, cflrate);
  const struct gkyl_array *mask_array = gkyl_dg_array_mask_get_mask(fdmul->cfl_mask);
  double omega_max = compute_global_array_max(app, fdmul, mask_array);
  clamp_cflrate_by_omega_max(fdmul, omega_max, cflrate, combined_multiplier);
}

static void
gk_species_fdot_multiplier_advance_time_dilation_cfl_f_frac_local(gkyl_gyrokinetic_app *app,
  const struct gk_species *gks, struct gk_fdot_multiplier_comp *fdmul, const struct gkyl_array *phi,
  const struct gkyl_array *f, const struct gkyl_array *cflrate,
  struct gkyl_array *combined_multiplier)
{
  gkyl_dg_array_mask_advance(fdmul->cfl_mask, f);
  gkyl_dg_array_mask_scale_by_cell(fdmul->cfl_mask, cflrate);
  const struct gkyl_array *mask_array = gkyl_dg_array_mask_get_mask(fdmul->cfl_mask);
  double omega_max = compute_global_array_max(app, fdmul, mask_array);
  clamp_cflrate_by_omega_max(fdmul, omega_max, cflrate, combined_multiplier);
}

static void
proj_on_basis_c2p_phase_func(const double *xcomp, double *xphys, void *ctx)
{
  struct gk_proj_on_basis_c2p_func_ctx *c2p_ctx = ctx;
  int cdim = c2p_ctx->cdim; // Assumes update range is a phase range.
  gkyl_position_map_eval_mc2nu(c2p_ctx->pos_map, xcomp, xphys);
  gkyl_velocity_map_eval_c2p(c2p_ctx->vel_map, &xcomp[cdim], &xphys[cdim]);
}

static void
proj_on_basis_c2p_position_func(const double *xcomp, double *xphys, void *ctx)
{
  struct gk_proj_on_basis_c2p_func_ctx *c2p_ctx = ctx;
  gkyl_position_map_eval_mc2nu(c2p_ctx->pos_map, xcomp, xphys);
}

static void
gk_species_fdot_multiplier_init_comp(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_fdot_multiplier_comp *fdmul,
  const struct gkyl_gyrokinetic_fdot_multiplier_comp *fdot_mult_inp)
{
  fdmul->type = fdot_mult_inp->type;

  // Default function pointer: no-op for GKYL_GK_FDOT_MULTIPLIER_NONE.
  fdmul->advance_func = gk_species_fdot_multiplier_advance_disabled;

  if (fdmul->type) {
    bool cellwise_const = fdot_mult_inp->cellwise_const;
    assert(cellwise_const); // MF 2025/06/11: Limited to this for now.

    struct gkyl_basis basis_mult;
    if (cellwise_const) {
      gkyl_cart_modal_serendip(&basis_mult, gks->basis.ndim, 0);
    }
    else {
      basis_mult = gks->basis;
    }

    fdmul->proj_on_basis_c2p_ctx.cdim = app->cdim;
    fdmul->proj_on_basis_c2p_ctx.vdim = gks->local_vel.ndim;
    fdmul->proj_on_basis_c2p_ctx.vel_map = gks->vel_map;
    fdmul->proj_on_basis_c2p_ctx.pos_map = app->position_map;

    fdmul->time_dilation_scale_const = fdot_mult_inp->time_dilation_scale_const;
    if (fdmul->time_dilation_scale_const <= 0.0 || fdmul->time_dilation_scale_const >= 1.0) {
      fdmul->time_dilation_scale_const = 1.0;
    }

    if (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_USER_INPUT) {
      fdmul->buffer = mkarr(app->use_gpu, basis_mult.num_basis, gks->local_ext.volume);
      struct gkyl_array *buffer_ho = mkarr(false, basis_mult.num_basis, gks->local_ext.volume);

      gkyl_proj_on_basis *projup = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp) {
        .grid = &gks->grid,
        .basis = &basis_mult,
        .num_quad = basis_mult.poly_order + 1,
        .num_ret_vals = 1,
        .eval = fdot_mult_inp->profile,
        .ctx = fdot_mult_inp->profile_ctx,
        .c2p_func = proj_on_basis_c2p_phase_func,
        .c2p_func_ctx = &fdmul->proj_on_basis_c2p_ctx,
      });
      gkyl_proj_on_basis_advance(projup, 0.0, &gks->local, buffer_ho);
      gkyl_proj_on_basis_release(projup);

      if (basis_mult.poly_order == 0) {
        gkyl_array_scale_range(buffer_ho, 1.0 / pow(sqrt(2.0), gks->grid.ndim), &gks->local);
      }

      gkyl_array_copy(fdmul->buffer, buffer_ho);
      gkyl_array_release(buffer_ho);

      fdmul->advance_func = gk_species_fdot_multiplier_advance_mult;
    }
    else if (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_LOSS_CONE) {
      // Operator that projects the loss cone mask.
      fdmul->phi_wall_lo_global = 0;
      fdmul->phi_wall_up_global = 0;
      int zdim = app->cdim - 1;
      struct gkyl_loss_cone_mask_gyrokinetic_inp inp_proj = {
        .conf_basis = &app->basis,
        .vel_map = gks->vel_map,
        .mass = gks->info.mass,
        .charge = gks->info.charge,
        .use_gpu = app->use_gpu,
        .lower_orbit = trapped_passing_orbit_from_species_bc(gks->lower_bc[zdim].type),
        .upper_orbit = trapped_passing_orbit_from_species_bc(gks->upper_bc[zdim].type),
      };
      fdmul->lcm_proj_op = gkyl_loss_cone_mask_gyrokinetic_inew(&inp_proj);

      fdmul->buffer = mkarr(app->use_gpu, basis_mult.num_basis, gks->local_ext.volume);
      fdmul->bmag_global = mkarr(app->use_gpu, app->gk_geom->geo_corn.bmag->ncomp,
        app->global_ext.volume);
      fdmul->phi_global = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);
      if (inp_proj.lower_orbit == GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_SHEATH)
        fdmul->phi_wall_lo_global = mkarr(app->use_gpu, app->basis.num_basis,
          app->global_ext.volume);
      if (inp_proj.upper_orbit == GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_SHEATH)
        fdmul->phi_wall_up_global = mkarr(app->use_gpu, app->basis.num_basis,
          app->global_ext.volume);

      gkyl_comm_array_allgather(app->comm, &app->local, &app->global, app->gk_geom->geo_corn.bmag,
        fdmul->bmag_global);

      fdmul->advance_func = gk_species_fdot_multiplier_advance_loss_cone_mult;
    }
    else if ((fdmul->type == GKYL_GK_FDOT_MULTIPLIER_FIXED_DT) ||
      (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_FIXED_FACTOR_TIMES_OMEGA_MAX) ||
      (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_FIXED_DT_OMEGAH) ||
      (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_DT_SET_BY_SPECIES) ||
      (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_MASK_F_THRESHOLD) ||
      (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_MASK_F_FRAC_LOCAL) ||
      (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_MASK_F_FRAC_GLOBAL) ||
      (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_CONSTANT)) {

      fdmul->cfl_dt_min_value = fdot_mult_inp->cfl_dt_min_value;
      fdmul->f_threshold = fdot_mult_inp->f_threshold;
      fdmul->cfl_factor_times_omega_max = fdot_mult_inp->cfl_factor_times_omega_max;

      if (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_CONSTANT) {
        fdmul->advance_func = gk_species_fdot_multiplier_advance_const;
      }
      else {
        fdmul->buffer = mkarr(app->use_gpu, basis_mult.num_basis, gks->local_ext.volume);

        if (app->use_gpu) {
#ifdef GKYL_HAVE_CUDA
          fdmul->omega_max_local_cu = (double *)gkyl_cu_malloc(sizeof(double));
          fdmul->local_max_f = (double *)gkyl_cu_malloc(sizeof(double) * gks->basis.num_basis);
          fdmul->global_max_f = (double *)gkyl_cu_malloc(sizeof(double) * gks->basis.num_basis);
#endif
        }
        else {
          fdmul->local_max_f = gkyl_malloc(sizeof(double) * gks->basis.num_basis);
          fdmul->global_max_f = gkyl_malloc(sizeof(double) * gks->basis.num_basis);
        }

        enum gkyl_dg_array_mask_types mask_type = GKYL_DG_ARRAY_MASK_NONE;
        switch (fdmul->type) {
          case GKYL_GK_FDOT_MULTIPLIER_FIXED_DT:
            mask_type = GKYL_DG_ARRAY_MASK_NONE;
            fdmul->advance_func =
              gk_species_fdot_multiplier_advance_time_dilation_cfl_dt_user_specified;
            break;
          case GKYL_GK_FDOT_MULTIPLIER_FIXED_FACTOR_TIMES_OMEGA_MAX:
            mask_type = GKYL_DG_ARRAY_MASK_NONE;
            fdmul->advance_func =
              gk_species_fdot_multiplier_advance_time_dilation_cfl_factor_user_specified;
            break;
          case GKYL_GK_FDOT_MULTIPLIER_FIXED_DT_OMEGAH:
            mask_type = GKYL_DG_ARRAY_MASK_NONE;
            fdmul->advance_func =
              gk_species_fdot_multiplier_advance_time_dilation_cfl_dt_omegaH;
            break;
          case GKYL_GK_FDOT_MULTIPLIER_DT_SET_BY_SPECIES:
            mask_type = GKYL_DG_ARRAY_MASK_NONE;
            fdmul->species_dt_is_set_from = gk_find_species(app, fdot_mult_inp->dt_set_by_species);
            fdmul->advance_func =
              gk_species_fdot_multiplier_advance_time_dilation_cfl_dt_set_by_species;
            break;
          case GKYL_GK_FDOT_MULTIPLIER_MASK_F_THRESHOLD:
            mask_type = GKYL_DG_ARRAY_MASK_C0_GREATER;
            fdmul->advance_func =
              gk_species_fdot_multiplier_advance_time_dilation_cfl_f_frac_local;
            break;
          case GKYL_GK_FDOT_MULTIPLIER_MASK_F_FRAC_LOCAL:
            mask_type = GKYL_DG_ARRAY_MASK_C0_GREATER_FRAC_CONF;
            fdmul->advance_func =
              gk_species_fdot_multiplier_advance_time_dilation_cfl_f_frac_local;
            break;
          case GKYL_GK_FDOT_MULTIPLIER_MASK_F_FRAC_GLOBAL:
            mask_type = GKYL_DG_ARRAY_MASK_C0_GREATER_FRAC;
            fdmul->advance_func =
              gk_species_fdot_multiplier_advance_time_dilation_cfl_f_frac_global;
            break;
          default:
            assert(false); // Unknown fdot_multiplier type.
            break;
        }

        struct gkyl_dg_array_mask_inp cfl_mask_inp = {
          .type = mask_type,
          .threshold = fdmul->f_threshold,
          .phase_rng = &gks->local,
          .phase_rng_ext = &gks->local_ext,
          .conf_rng = &app->local,
          .conf_rng_ext = &app->local_ext,
          .vel_rng = &gks->local_vel,
          .use_gpu = app->use_gpu,
        };
        fdmul->cfl_mask = gkyl_dg_array_mask_new(cfl_mask_inp);
      }
    }
  }
}

void
gk_species_fdot_multiplier_init(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_fdot_multiplier *fdot_mult)
{
  fdot_mult->num_multipliers = gks->info.time_rate_multiplier.num_multipliers;
  fdot_mult->write_func = gk_species_fdot_multiplier_write_disabled;
  fdot_mult->advance_times_cfl_func = gk_fdot_multiplier_advance_times_cfl_disabled;
  fdot_mult->advance_times_rate_func = gk_fdot_multiplier_advance_times_rate_disabled;

  if (fdot_mult->num_multipliers == 0) {
    return;
  }

  // Initialize multiplier components and determine if any component requests diagnostics.
  fdot_mult->write_diagnostics = false;
  for (int i = 0; i < fdot_mult->num_multipliers; ++i) {
    gk_species_fdot_multiplier_init_comp(app, gks, &fdot_mult->comp[i],
      &gks->info.time_rate_multiplier.multiplier[i]);
    if (gks->info.time_rate_multiplier.multiplier[i].write_diagnostics) {
      fdot_mult->write_diagnostics = true;
    }
  }

  // Allocate the single combined multiplier array (cellwise-constant, ncomp=1).
  fdot_mult->multiplier = mkarr(app->use_gpu, 1, gks->local_ext.volume);
  gkyl_array_clear(fdot_mult->multiplier, 1.0);
  if (fdot_mult->write_diagnostics) {
    fdot_mult->multiplier_host = app->use_gpu ?
      mkarr(false, fdot_mult->multiplier->ncomp, fdot_mult->multiplier->size) :
      gkyl_array_acquire(fdot_mult->multiplier);
    fdot_mult->write_func = gk_species_fdot_multiplier_write_enabled;
  }

  fdot_mult->advance_times_cfl_func = gk_fdot_multiplier_advance_times_cfl_enabled;
  fdot_mult->advance_times_rate_func = gk_fdot_multiplier_advance_times_rate_enabled;
}

void
gk_species_fdot_multiplier_advance_times_cfl(gkyl_gyrokinetic_app *app,
  struct gk_species *gks, struct gk_fdot_multiplier *fdot_mult,
  const struct gkyl_array *phi, const struct gkyl_array *f, struct gkyl_array *cflrate)
{
  fdot_mult->advance_times_cfl_func(app, gks, fdot_mult, phi, f, cflrate);
}

void
gk_species_fdot_multiplier_advance_times_rate(gkyl_gyrokinetic_app *app,
  struct gk_species *gks, struct gk_fdot_multiplier *fdot_mult,
  const struct gkyl_array *phi, const struct gkyl_array *f, struct gkyl_array *rhs)
{
  fdot_mult->advance_times_rate_func(app, gks, fdot_mult, phi, f, rhs);
}

void
gk_species_fdot_multiplier_write(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_fdot_multiplier *fdot_mult, double tm, int frame)
{
  fdot_mult->write_func(app, gks, fdot_mult, tm, frame);
}

double
gk_fdot_multiplier_get_time_dilation_scale_const(gkyl_gyrokinetic_app *app,
  const struct gk_fdot_multiplier *fdot_mult)
{
  double result = 1.0;
  for (int i = 0; i < fdot_mult->num_multipliers; ++i) {
    if (fdot_mult->comp[i].time_dilation_scale_const) {
      result *= fdot_mult->comp[i].time_dilation_scale_const;
    }
  }
  return result;
}

static void
gk_species_fdot_multiplier_release_comp(const struct gkyl_gyrokinetic_app *app,
  const struct gk_fdot_multiplier_comp *fdmul)
{
  if (!fdmul->type) return;

  if (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_USER_INPUT) {
    gkyl_array_release(fdmul->buffer);
  }
  else if (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_LOSS_CONE) {
    gkyl_array_release(fdmul->buffer);
    gkyl_array_release(fdmul->bmag_global);
    gkyl_array_release(fdmul->phi_global);
    if (fdmul->phi_wall_lo_global)
      gkyl_array_release(fdmul->phi_wall_lo_global);
    if (fdmul->phi_wall_up_global)
      gkyl_array_release(fdmul->phi_wall_up_global);
    gkyl_loss_cone_mask_gyrokinetic_release(fdmul->lcm_proj_op);
  }
  else if (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_CONSTANT) {
    // No allocations to free for CONSTANT type.
  }
  else if ((fdmul->type == GKYL_GK_FDOT_MULTIPLIER_FIXED_DT) ||
    (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_FIXED_FACTOR_TIMES_OMEGA_MAX) ||
    (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_FIXED_DT_OMEGAH) ||
    (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_DT_SET_BY_SPECIES) ||
    (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_MASK_F_THRESHOLD) ||
    (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_MASK_F_FRAC_LOCAL) ||
    (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_MASK_F_FRAC_GLOBAL)) {
    gkyl_array_release(fdmul->buffer);
    gkyl_dg_array_mask_release(fdmul->cfl_mask);
    if (app->use_gpu) {
#ifdef GKYL_HAVE_CUDA
      gkyl_cu_free(fdmul->omega_max_local_cu);
      gkyl_cu_free(fdmul->local_max_f);
      gkyl_cu_free(fdmul->global_max_f);
#endif
    }
    else {
      gkyl_free(fdmul->local_max_f);
      gkyl_free(fdmul->global_max_f);
    }
  }
}

void
gk_species_fdot_multiplier_release(const struct gkyl_gyrokinetic_app *app,
  const struct gk_fdot_multiplier *fdot_mult)
{
  if (fdot_mult->num_multipliers == 0) return;

  for (int i = 0; i < fdot_mult->num_multipliers; ++i) {
    gk_species_fdot_multiplier_release_comp(app, &fdot_mult->comp[i]);
  }
  gkyl_array_release(fdot_mult->multiplier);
  if (fdot_mult->write_diagnostics) {
    gkyl_array_release(fdot_mult->multiplier_host);
  }
}

void
gk_species_fdot_multiplier_reset(gkyl_gyrokinetic_app *app, double tm, struct gk_species *gks,
  struct gk_fdot_multiplier *fdot_mult,
  struct gkyl_gyrokinetic_fdot_multiplier fdot_mult_inp)
{
  gk_species_fdot_multiplier_release(app, fdot_mult);

  gks->info.time_rate_multiplier = fdot_mult_inp;

  gk_species_fdot_multiplier_init(app, gks, fdot_mult);
}
