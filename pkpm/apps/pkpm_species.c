#include <math.h>
#include <stdlib.h>
#include <gkyl_alloc.h>
#include <gkyl_app.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_bc_basic.h>
#include <gkyl_dynvec.h>
#include <gkyl_elem_type.h>
#include <gkyl_eqn_type.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_pkpm_priv.h>

#include <assert.h>
#include <time.h>

// initialize species object

// Reflecting wall of the distributions on the lower or upper edge of direction dir.
static struct gkyl_bc_pkpm_species_reflect *
pkpm_species_reflect_new(
  gkyl_pkpm_app *app, struct pkpm_species *s, int dir, enum gkyl_edge_loc edge
)
{
  return gkyl_bc_pkpm_species_reflect_inew(&(struct gkyl_bc_pkpm_species_reflect_inp){
    .dir = dir,
    .edge = edge,
    .cdim = app->cdim,
    .phase_basis = &app->basis,
    .phase_basis_on_dev = app->basis_on_dev.basis,
    .conf_num_basis = app->confBasis.num_basis,
    .num_blocks = PKPM_DIST_NUM_BLOCKS,
    .skin_r = edge == GKYL_LOWER_EDGE ? &s->lower_skin_dist[dir] : &s->upper_skin_dist[dir],
    .ghost_r = edge == GKYL_LOWER_EDGE ? &s->lower_ghost_dist[dir] : &s->upper_ghost_dist[dir],
    .phase_range_ext = &s->local_ext,
    .conf_range_ext = &app->local_ext,
    .use_gpu = app->use_gpu,
  });
}

void
pkpm_species_init(struct gkyl_pkpm *pkpm, struct gkyl_pkpm_app *app, struct pkpm_species *s)
{
  int cdim = app->cdim, vdim = app->vdim;
  int pdim = cdim + vdim;

  int cells[GKYL_MAX_DIM], ghost[GKYL_MAX_DIM];
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];

  int cells_vel[GKYL_MAX_DIM], ghost_vel[GKYL_MAX_DIM];
  double lower_vel[GKYL_MAX_DIM], upper_vel[GKYL_MAX_DIM];

  for (int d = 0; d < cdim; ++d) {
    cells[d] = pkpm->cells[d];
    lower[d] = pkpm->lower[d];
    upper[d] = pkpm->upper[d];
    ghost[d] = 1;
  }
  for (int d = 0; d < vdim; ++d) {
    // full phase space grid
    // full phase space grid (the velocity coordinate is w = v_par/xi, normalized)
    cells[cdim + d] = s->info.cells[d];
    lower[cdim + d] = s->info.lower[d];
    upper[cdim + d] = s->info.upper[d];
    ghost[cdim + d] = 0; // no ghost-cells in velocity space

    // only velocity space
    cells_vel[d] = s->info.cells[d];
    lower_vel[d] = s->info.lower[d];
    upper_vel[d] = s->info.upper[d];
    ghost_vel[d] = 0; // no ghost-cells in velocity space
  }
  // full phase space grid
  gkyl_rect_grid_init(&s->grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&s->grid, ghost, &s->global_ext, &s->global);

  // velocity space grid
  gkyl_rect_grid_init(&s->grid_vel, vdim, lower_vel, upper_vel, cells_vel);
  gkyl_create_grid_ranges(&s->grid_vel, ghost_vel, &s->local_ext_vel, &s->local_vel);

  // phase-space communicator
  s->comm = gkyl_comm_extend_comm(app->comm, &s->local_vel);

  // create local and local_ext from app local range: local = conf-local X local_vel.
  // The extended range has the configuration ghost cells and ONE velocity ghost cell
  // above w_max: the empty state beyond the characteristic velocity edge, whose face slot
  // holds the upper-edge flux (the lower edge is the lowest cell's own slot).
  // local_ext_conf is extended in configuration space only: the range of everything that
  // integrates, projects, syncs or reflects the distribution.
  struct gkyl_range local;
  gkyl_range_ten_prod(&local, &app->local, &s->local_vel);
  int elo[GKYL_MAX_DIM] = {0}, eup[GKYL_MAX_DIM] = {0};
  for (int d = 0; d < cdim; ++d) {
    elo[d] = eup[d] = 1;
  }
  eup[cdim] = 1;
  gkyl_range_extend(&s->local_ext, &local, elo, eup);
  gkyl_sub_range_init(&s->local, &s->local_ext, local.lower, local.upper);
  int lower_conf[GKYL_MAX_DIM], upper_conf[GKYL_MAX_DIM];
  for (int d = 0; d < pdim; ++d) {
    lower_conf[d] = d < cdim ? s->local_ext.lower[d] : local.lower[d];
    upper_conf[d] = d < cdim ? s->local_ext.upper[d] : local.upper[d];
  }
  gkyl_sub_range_init(&s->local_ext_conf, &s->local_ext, lower_conf, upper_conf);

  // allocate distribution function arrays: the normalized distributions
  // [F, G, M0_perp (x, y, z), N_perp (x, y, z)] (enum pkpm_dist_block)
  int nc_dist = PKPM_DIST_NUM_BLOCKS * app->basis.num_basis;
  s->f = mkarr(app->use_gpu, nc_dist, s->local_ext.volume);
  s->f1 = mkarr(app->use_gpu, nc_dist, s->local_ext.volume);
  // physical distributions (position-map Jacobian divided out) for the output and the L2 norm
  s->f_no_J = mkarr(app->use_gpu, nc_dist, s->local_ext.volume);
  s->fnew = mkarr(app->use_gpu, nc_dist, s->local_ext.volume);
  gkyl_array_clear(s->f, 0.0);
  gkyl_array_clear(s->f1, 0.0);
  gkyl_array_clear(s->f_no_J, 0.0);
  gkyl_array_clear(s->fnew, 0.0);
  s->blk_in = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
  s->blk_out = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
  // allocate fluid arrays [rho, rho ux, rho uy, rho uz, E]
  s->fluid = mkarr(app->use_gpu, 5 * app->confBasis.num_basis, app->local_ext.volume);
  s->fluid1 = mkarr(app->use_gpu, 5 * app->confBasis.num_basis, app->local_ext.volume);
  s->fluidnew = mkarr(app->use_gpu, 5 * app->confBasis.num_basis, app->local_ext.volume);

  // create host arrays if on GPUs so initialization occurs host-side
  s->f_host = s->f;
  s->f_no_J_host = s->f_no_J;
  s->fluid_host = s->fluid;
  if (app->use_gpu) {
    s->f_host = mkarr(false, nc_dist, s->local_ext.volume);
    s->f_no_J_host = mkarr(false, nc_dist, s->local_ext.volume);
    s->fluid_host = mkarr(false, 5 * app->confBasis.num_basis, app->local_ext.volume);
  }

  // Duplicate copy of momentum data in case time step fails.
  // Needed because of implicit source split which modifies solution and
  // is always successful, so if a time step fails due to the SSP RK3
  // we must restore the old solution before restarting the time step
  s->fluid_dup = mkarr(app->use_gpu, 5 * app->confBasis.num_basis, app->local_ext.volume);

  // Wave equation object for upwinding fluid equations
  // We use the 10 moment system since PKPM model generates a full pressure tensor
  // k0 = 0.0 and use_grad_closure = false because we do not need a closure
  // (the kinetic equations are the closure)

  // allocate cflrate for distribution function and momentum updates (scalar arrays)
  s->cflrate_f = mkarr(app->use_gpu, 1, s->local_ext.volume);
  s->cflrate_fluid = mkarr(app->use_gpu, 1, app->local_ext.volume);

  if (app->use_gpu) {
    s->omegaCfl_ptr_dist = gkyl_cu_malloc(sizeof(double));
    s->omegaCfl_ptr_fluid = gkyl_cu_malloc(sizeof(double));
  } else {
    s->omegaCfl_ptr_dist = gkyl_malloc(sizeof(double));
    s->omegaCfl_ptr_fluid = gkyl_malloc(sizeof(double));
  }

  // allocate array to store q/m*(E,B) if using a fully explicit update
  s->qmem = mkarr(app->use_gpu, 8 * app->confBasis.num_basis, app->local_ext.volume);

  int nc = app->confBasis.num_basis;
  // pkpm moments for diagnostics [rho, M1, p_par, p_perp, q_par, q_perp, r_parpar, r_parperp]
  // (the is_integrated flag selects the diagnostic moment set)
  pkpm_species_moment_init(app, s, &s->pkpm_moms_diag, true);
  // raw moments of the normalized distributions [F, G, M0_perp, N_perp] (mass 1), the kinetic
  // coupling moments [p_par, p_perp, q_par, q_perp, s (3), r (3)] formed from them (plan 17.14) and
  // the constraints
  pkpm_species_moment_init_mass(app, s, &s->pkpm_nmoms, false, 1.0);
  s->pkpm_moms = mkarr(app->use_gpu, 10 * nc, app->local_ext.volume);
  s->pkpm_moms_no_J = mkarr(app->use_gpu, 10 * nc, app->local_ext.volume);
  s->fluid_no_J = mkarr(app->use_gpu, 5 * nc, app->local_ext.volume);
  gkyl_array_clear(s->pkpm_moms, 0.0);
  s->pkpm_cmoms = mkarr(app->use_gpu, 3 * nc, app->local_ext.volume);
  s->pkpm_cmoms_host = s->pkpm_cmoms;
  if (app->use_gpu) {
    s->pkpm_cmoms_host = mkarr(false, 3 * nc, app->local_ext.volume);
  }

  // Refactored update variables (p=1 configuration space only).
  s->pkpm_u = mkarr(app->use_gpu, 3 * nc, app->local_ext.volume);
  s->pkpm_grad_u = mkarr(app->use_gpu, cdim * 3 * nc, app->local_ext.volume);
  gkyl_array_clear(s->pkpm_grad_u, 0.0);
  s->pkpm_theta = mkarr(app->use_gpu, nc, app->local_ext.volume);
  s->pkpm_rho = mkarr(app->use_gpu, nc, app->local_ext.volume);
  s->pkpm_rho3 = mkarr(app->use_gpu, 3 * nc, app->local_ext.volume);
  gkyl_array_clear(s->pkpm_rho3, 0.0);
  s->pkpm_mom3 = mkarr(app->use_gpu, 3 * nc, app->local_ext.volume);
  s->pkpm_grad_rho = mkarr(app->use_gpu, cdim * 3 * nc, app->local_ext.volume);
  gkyl_array_clear(s->pkpm_grad_rho, 0.0);
  s->pkpm_grad_ln_n = mkarr(app->use_gpu, cdim * nc, app->local_ext.volume);
  gkyl_array_clear(s->pkpm_grad_ln_n, 0.0);
  s->pkpm_c_0 = mkarr(app->use_gpu, nc, app->local_ext.volume);
  gkyl_array_clear(s->pkpm_c_0, 0.0);
  s->pkpm_alpha_0 = mkarr(app->use_gpu, nc, app->local_ext.volume);
  gkyl_array_clear(s->pkpm_alpha_0, 0.0);
  s->pkpm_c_h = mkarr(app->use_gpu, nc, app->local_ext.volume);
  gkyl_array_clear(s->pkpm_alpha_0, 0.0);
  gkyl_array_clear(s->pkpm_c_h, 0.0);
  gkyl_array_clear(s->pkpm_c_h, 0.0);
  s->pkpm_a_perp = mkarr(app->use_gpu, 3 * nc, app->local_ext.volume);
  gkyl_array_clear(s->pkpm_a_perp, 0.0);
  s->pkpm_xi = mkarr(app->use_gpu, 3 * nc, app->local_ext.volume);
  gkyl_array_clear(s->pkpm_xi, 0.0);
  s->pkpm_xi3 = mkarr(app->use_gpu, 3 * nc, app->local_ext.volume);
  gkyl_array_clear(s->pkpm_xi3, 0.0);
  s->pkpm_grad_xi = mkarr(app->use_gpu, cdim * 3 * nc, app->local_ext.volume);
  gkyl_array_clear(s->pkpm_grad_xi, 0.0);
  s->pkpm_scratch = mkarr(app->use_gpu, nc, app->local_ext.volume);
  s->pkpm_rhs_moms = mkarr(
    app->use_gpu, gkyl_dg_updater_moment_pkpm_num_mom(s->pkpm_nmoms.mcalc) * nc,
    app->local_ext.volume
  );
  s->pkpm_p_ij = mkarr(app->use_gpu, 6 * nc, app->local_ext.volume);
  gkyl_array_clear(s->pkpm_theta, 0.0);
  s->num_null_rho = 0;
  s->num_null_closure = 0;
  s->closure_cond_min = 1.0;
  s->closure_null_warned = false;

  // Updater for the configuration-space variables (nodal primitive variables, closures, weak
  // gradients, sources, diagnostics); every method takes the range it acts on.
  s->calc_pkpm_vars =
    gkyl_dg_calc_pkpm_vars_new(&app->grid, &app->confBasis, app->pos_map, app->use_gpu);

  // array for storing integrated fluid variables in each cell
  // integ_pkpm_mom : integral[rho, rhoux, rhouy, rhouz, rhoux^2, rhouy^2, rhouz^2, p_par, p_perp]
  s->integ_pkpm_mom = mkarr(app->use_gpu, 9, app->local_ext.volume);
  // arrays for I/O, fluid_io and pkpm_vars_io
  // fluid_io : [rho, rhoux, rhouy, rhouz, P_xx+rhoux^2, P_xy+rhouxuy, P_xz+rhouxuz, P_yy+rhouy^2, P_yz+rhouyuz, P_zz+rhouz^2, E]
  // pkpm_vars_io : [ux, uy, uz, T_perp/m, a_P, div(b), bb:grad(u), T_perp/m div(b), xi, c]
  s->fluid_io = mkarr(app->use_gpu, 11 * nc, app->local_ext.volume);
  s->pkpm_vars_io = mkarr(app->use_gpu, 10 * nc, app->local_ext.volume);
  s->fluid_io_host = s->fluid_io;
  s->pkpm_vars_io_host = s->pkpm_vars_io;
  if (app->use_gpu) {
    s->fluid_io_host = mkarr(false, 11 * nc, app->local_ext.volume);
    s->pkpm_vars_io_host = mkarr(false, 10 * nc, app->local_ext.volume);
  }
  // Record of the positivity limiter: per cell [density scalings, energy scalings,
  // inadmissible averages] accumulated over the run, and its p=1 image for IO.
  s->pkpm_pos_count = mkarr(app->use_gpu, 3, app->local_ext.volume);
  gkyl_array_clear(s->pkpm_pos_count, 0.0);
  s->pkpm_pos_io = mkarr(app->use_gpu, 3 * nc, app->local_ext.volume);
  s->pkpm_pos_io_host = s->pkpm_pos_io;
  if (app->use_gpu) {
    s->pkpm_pos_io_host = mkarr(false, 3 * nc, app->local_ext.volume);
  }

  // Surface-flux updaters and stored modal fluxes (each cell owns its lower faces).
  struct gkyl_dg_pkpm_conf_flux_surf_inp conf_flux_inp = {
    .phase_grid = &s->grid,
    .conf_basis = &app->confBasis,
    .phase_basis = &app->basis,
    .pos_map = app->pos_map,
    .use_gpu = app->use_gpu,
  };
  s->calc_conf_flux = gkyl_dg_pkpm_conf_flux_surf_inew(&conf_flux_inp);
  int nkc = gkyl_dg_pkpm_conf_flux_surf_num_surf_basis(s->calc_conf_flux);
  s->conf_flux_surf = mkarr(app->use_gpu, 8 * cdim * nkc, s->local_ext.volume);
  s->face_coef = mkarr(
    app->use_gpu, gkyl_dg_pkpm_conf_flux_surf_face_coef_ncomp(s->calc_conf_flux),
    app->local_ext.volume
  );

  struct gkyl_dg_pkpm_vel_flux_surf_inp vel_flux_inp = {
    .phase_grid = &s->grid,
    .conf_basis = &app->confBasis,
    .phase_basis = &app->basis,
    .use_gpu = app->use_gpu,
  };
  s->calc_vel_flux = gkyl_dg_pkpm_vel_flux_surf_inew(&vel_flux_inp);
  int nkv = gkyl_dg_pkpm_vel_flux_surf_num_surf_basis(s->calc_vel_flux);
  s->vel_flux_surf = mkarr(app->use_gpu, PKPM_DIST_NUM_BLOCKS * nkv, s->local_ext.volume);

  struct gkyl_dg_pkpm_fluid_flux_surf_inp fluid_flux_inp = {
    .conf_grid = &app->grid,
    .conf_basis = &app->confBasis,
    .phase_grid = &s->grid,
    .pos_map = app->pos_map,
    .use_gpu = app->use_gpu,
  };
  s->calc_fluid_flux = gkyl_dg_pkpm_fluid_flux_surf_inew(&fluid_flux_inp);
  int nkf = gkyl_dg_pkpm_fluid_flux_surf_num_surf_basis(s->calc_fluid_flux);
  s->fluid_flux_surf = mkarr(app->use_gpu, 5 * cdim * nkf, app->local_ext.volume);
  s->fluid_flux_surf_host = app->use_gpu ? mkarr(false, 5 * cdim * nkf, app->local_ext.volume) : 0;

  // Equation objects and hyper_dg solvers of the three passes.
  struct gkyl_dg_pkpm_kinetic_conf_inp kinetic_conf_inp = {
    .conf_basis = &app->confBasis,
    .phase_basis = &app->basis,
    .conf_range = &app->local,
    .phase_range = &s->local,
    .pos_map = app->pos_map,
    .auxfields =
      {
        .bvar = app->field->bvar,
        .u = s->pkpm_u,
        .xi = s->pkpm_xi,
        .conf_flux_surf = s->conf_flux_surf,
      },
    .use_gpu = app->use_gpu,
  };
  s->eqn_kinetic_conf = gkyl_dg_pkpm_kinetic_conf_inew(&kinetic_conf_inp);
  struct gkyl_dg_pkpm_harm_vol_inp harm_vol_inp = {
    .phase_grid = &s->grid,
    .conf_basis = &app->confBasis,
    .phase_basis = &app->basis,
    .qbym = s->info.charge / s->info.mass,
    .pos_map = app->pos_map,
    .use_gpu = app->use_gpu,
  };
  s->harm_vol = gkyl_dg_pkpm_harm_vol_inew(&harm_vol_inp);
  s->harm_coef =
    mkarr(app->use_gpu, gkyl_dg_pkpm_harm_vol_num_coef(s->harm_vol), app->local_ext.volume);
  s->harm_wnode =
    mkarr(app->use_gpu, gkyl_dg_pkpm_harm_vol_num_wnode(s->harm_vol), app->local_ext.volume);
  // gyrophase mixing: default rate 10 |Omega| (plan 17.31: the bulk overshoot is ~(drive/Omega)/c, and
  // problems driven at the ion gyroscale drive agyrotropy at ~Omega), switched off by a negative input
  double gyro_mix = s->info.gyro_mix == 0.0 ? 10.0 : s->info.gyro_mix;
  s->harm_exp = gkyl_dg_pkpm_harm_exp_new(&(struct gkyl_dg_pkpm_harm_exp_inp){
    .phase_grid = &s->grid,
    .conf_basis = &app->confBasis,
    .qbym = s->info.charge / s->info.mass,
    .mix = gyro_mix,
    .use_gpu = app->use_gpu,
  });
  s->harm_exp_coef =
    mkarr(app->use_gpu, gkyl_dg_pkpm_harm_exp_num_coef(s->harm_exp), app->local_ext.volume);
  s->harm_exp_nu =
    mkarr(app->use_gpu, gkyl_dg_pkpm_harm_exp_num_rate(s->harm_exp), s->local_ext.volume);
  s->harm_exp_act = mkarr(app->use_gpu, 1, app->local_ext.volume);
  s->harm_exp_col =
    mkarr(app->use_gpu, gkyl_dg_pkpm_harm_exp_num_col(s->harm_exp), app->local_ext.volume);
  struct gkyl_dg_pkpm_kinetic_vel_inp kinetic_vel_inp = {
    .conf_basis = &app->confBasis,
    .phase_basis = &app->basis,
    .conf_range = &app->local,
    .phase_range = &s->local,
    .auxfields =
      {
        .c_0 = s->pkpm_c_0,
        .alpha_0 = s->pkpm_alpha_0,
        .c_h = s->pkpm_c_h,
        .a_perp = s->pkpm_a_perp,
        .vel_flux_surf = s->vel_flux_surf,
      },
    .use_gpu = app->use_gpu,
  };
  s->eqn_kinetic_vel = gkyl_dg_pkpm_kinetic_vel_inew(&kinetic_vel_inp);
  struct gkyl_dg_pkpm_fluid_inp fluid_inp = {
    .conf_basis = &app->confBasis,
    .conf_range = &app->local,
    .pos_map = app->pos_map,
    .auxfields =
      {
        .vlasov_pkpm_moms = s->pkpm_moms,
        .bvar = app->field->bvar,
        .u = s->pkpm_u,
        .fluid_flux_surf = s->fluid_flux_surf,
      },
    .use_gpu = app->use_gpu,
  };
  s->eqn_fluid = gkyl_dg_pkpm_fluid_inew(&fluid_inp);

  int up_dirs_conf[GKYL_MAX_DIM] = {0, 1, 2}, zero_flux_conf[2 * GKYL_MAX_DIM] = {0};
  s->slvr_kinetic_conf = gkyl_hyper_dg_new(
    &s->grid, &app->basis, s->eqn_kinetic_conf, cdim, up_dirs_conf, zero_flux_conf, 1, app->use_gpu
  );
  int up_dirs_vel[GKYL_MAX_DIM] = {cdim}, zero_flux_vel[2 * GKYL_MAX_DIM] = {0};
  // The lower velocity edge has no ghost cell (boundary term with the lowest cell's own
  // flux slot); the upper edge is the regular face between the uppermost cell and the
  // velocity ghost cell, whose slot the closure fills.
  zero_flux_vel[cdim] = 1;
  s->slvr_kinetic_vel = gkyl_hyper_dg_new(
    &s->grid, &app->basis, s->eqn_kinetic_vel, 1, up_dirs_vel, zero_flux_vel, 1, app->use_gpu
  );
  s->slvr_fluid = gkyl_hyper_dg_new(
    &app->grid, &app->confBasis, s->eqn_fluid, cdim, up_dirs_conf, zero_flux_conf, 1, app->use_gpu
  );
  s->kinetic_tm = 0.0;
  s->fluid_tm = 0.0;

  // array for storing F_0^2 (0th Laguerre coefficient) in each cell
  s->L2_f = mkarr(app->use_gpu, 1, s->local_ext.volume);
  if (app->use_gpu) {
    s->red_L2_f = gkyl_cu_malloc(sizeof(double));
    s->red_integ_diag = gkyl_cu_malloc(sizeof(double[9]));
  }
  // allocate dynamic-vector to store all-reduced integrated moments and F_0^2
  s->integ_L2_f = gkyl_dynvec_new(GKYL_DOUBLE, 1);
  s->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, 9);
  s->pos_stat = gkyl_dynvec_new(GKYL_DOUBLE, 4);
  s->bflux_diag = gkyl_dynvec_new(GKYL_DOUBLE, 2 * cdim);
  s->is_first_bflux_write_call = true;
  for (int i = 0; i < 2 * GKYL_MAX_CDIM; ++i) {
    s->bflux[i] = 0.0;
  }
  s->is_first_pos_stat_write_call = true;
  for (int c = 0; c < 4; ++c) {
    s->pos_count_since[c] = 0;
  }
  s->is_first_integ_L2_write_call = true;
  s->is_first_integ_write_call = true;

  // Initialize applied acceleration for use in force update.
  s->app_accel = mkarr(app->use_gpu, 3 * app->confBasis.num_basis, app->local_ext.volume);
  gkyl_array_clear(s->app_accel, 0.0);
  s->has_app_accel = false;
  s->app_accel_evolve = false;
  // setup applied acceleration
  if (s->info.app_accel) {
    s->has_app_accel = true;
    if (s->info.app_accel_evolve) {
      s->app_accel_evolve = s->info.app_accel_evolve;
    }

    s->app_accel_host = s->app_accel;
    if (app->use_gpu) {
      s->app_accel_host = mkarr(false, 3 * app->confBasis.num_basis, app->local_ext.volume);
    }
    s->app_accel_proj = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
      .grid = &app->grid,
      .basis = &app->confBasis,
      .num_quad = app->confBasis.poly_order + 1,
      .num_ret_vals = 3,
      .eval = s->info.app_accel,
      .ctx = s->info.app_accel_ctx,
      .c2p_func = pkpm_pos_c2p,
      .c2p_func_ctx = app->pos_map,
    });
  }

  // determine collision type to use in PKPM update
  s->collision_id = s->info.collisions.collision_id;
  if (s->collision_id == GKYL_BGK_COLLISIONS) {
    pkpm_species_bgk_init(app, s, &s->bgk);
  }
  if (s->collision_id == GKYL_LBO_COLLISIONS) {
    pkpm_species_lbo_init(app, s, &s->lbo);
  }

  // momentum diffusion, if present, on the position map (orders 2 and 4)
  s->has_diffusion = false;
  if (s->info.diffusion.D) {
    s->has_diffusion = true;
    int order = s->info.diffusion.order < 2 ? 2 : s->info.diffusion.order;
    if (order != 2 && order != 4) {
      gkyl_exit("pkpm_species: momentum diffusion supports orders 2 and 4 only.");
    }
    struct gkyl_dg_pkpm_mom_diffusion_inp diff_inp = {
      .conf_grid = &app->grid,
      .conf_basis = &app->confBasis,
      .pos_map = app->pos_map,
      .order = order,
      .use_gpu = app->use_gpu,
    };
    for (int d = 0; d < cdim; ++d) {
      diff_inp.D[d] = s->info.diffusion.D;
    }
    s->mom_diff = gkyl_dg_pkpm_mom_diffusion_inew(&diff_inp);
  }

  // determine which directions are not periodic and get non-periodic BC info
  // also create local skin-ghost ranges for ease of handling non-periodic BCs
  int num_periodic_dir = app->num_periodic_dir, is_np[3] = {1, 1, 1};
  for (int d = 0; d < num_periodic_dir; ++d) {
    is_np[app->periodic_dirs[d]] = 0;
  }

  for (int dir = 0; dir < cdim; ++dir) {
    s->lower_bc[dir] = s->upper_bc[dir] = GKYL_SPECIES_COPY;
    if (is_np[dir]) {
      const enum gkyl_species_bc_type *bc;
      if (dir == 0) {
        bc = s->info.bcx;
      } else if (dir == 1) {
        bc = s->info.bcy;
      } else {
        bc = s->info.bcz;
      }

      s->lower_bc[dir] = bc[0];
      s->upper_bc[dir] = bc[1];
    }
    // Create local lower skin and ghost ranges for distribution function
    gkyl_skin_ghost_ranges(
      &s->lower_skin_dist[dir], &s->lower_ghost_dist[dir], dir, GKYL_LOWER_EDGE, &s->local_ext_conf,
      ghost
    );
    // Create local upper skin and ghost ranges for distribution function
    gkyl_skin_ghost_ranges(
      &s->upper_skin_dist[dir], &s->upper_ghost_dist[dir], dir, GKYL_UPPER_EDGE, &s->local_ext_conf,
      ghost
    );
  }

  // allocate buffer for applying BCs
  long buff_sz_dist = 0;
  long buff_sz_fluid = 0;
  // compute buffer size needed
  for (int dir = 0; dir < cdim; ++dir) {
    long vol_dist = GKYL_MAX2(s->lower_skin_dist[dir].volume, s->upper_skin_dist[dir].volume);
    buff_sz_dist = buff_sz_dist > vol_dist ? buff_sz_dist : vol_dist;
    long vol_fluid = GKYL_MAX2(app->lower_skin[dir].volume, app->upper_skin[dir].volume);
    buff_sz_fluid = buff_sz_fluid > vol_fluid ? buff_sz_fluid : vol_fluid;
  }

  s->bc_buffer_dist = mkarr(app->use_gpu, s->f->ncomp, buff_sz_dist);
  // buffer arrays for fixed function boundary conditions on distribution function
  s->bc_buffer_lo_fixed_dist = mkarr(app->use_gpu, s->f->ncomp, buff_sz_dist);
  s->bc_buffer_up_fixed_dist = mkarr(app->use_gpu, s->f->ncomp, buff_sz_dist);
  gkyl_array_clear(s->bc_buffer_lo_fixed_dist, 0.0);
  gkyl_array_clear(s->bc_buffer_up_fixed_dist, 0.0);

  s->bc_buffer_fluid = mkarr(app->use_gpu, 5 * app->confBasis.num_basis, buff_sz_fluid);
  // buffer arrays for fixed function boundary conditions on momentum
  s->bc_buffer_lo_fixed_fluid = mkarr(app->use_gpu, 5 * app->confBasis.num_basis, buff_sz_fluid);
  s->bc_buffer_up_fixed_fluid = mkarr(app->use_gpu, 5 * app->confBasis.num_basis, buff_sz_fluid);

  // Absorbing boundaries: the distribution functions vanish in the ghost cells; the
  // fluid state is the mirror image of the skin cell (even in the wall direction), so the
  // ghost trace at the wall is the skin trace and the primitive variables, their weak
  // gradients and the face averages of the flux coefficients are one-sided at the wall (a
  // zero ghost state would be a huge gradient of xi, n and u that drives spurious forces in
  // the wall cell; a copied one halves the wall-cell gradients). The fluid fluxes through
  // the wall are the moments of the kinetic flux against the zero ghost distributions
  // (pkpm_species_fluid_flux), not the fluid flux of the ghost state.
  s->bc_is_absorb = false;
  for (int d = 0; d < cdim; ++d) {
    s->bc_lo_reflect[d] = 0;
    s->bc_up_reflect[d] = 0;
    s->bc_lo_theta[d] = s->bc_up_theta[d] = 0;
    // Lower BC updater. Copy BCs by default.
    enum gkyl_bc_basic_type bctype_dist = GKYL_BC_COPY;
    enum gkyl_bc_basic_type bctype_fluid = GKYL_BC_COPY;
    if (s->lower_bc[d] == GKYL_SPECIES_COPY) {
      bctype_dist = GKYL_BC_COPY;
      bctype_fluid = GKYL_BC_COPY;
    } else if (s->lower_bc[d] == GKYL_SPECIES_ABSORB) {
      bctype_dist = GKYL_BC_ABSORB;
      bctype_fluid = GKYL_BC_REFLECT;
      s->bc_is_absorb = true;
      s->bc_lo_theta[d] = gkyl_bc_basic_new(
        d, GKYL_LOWER_EDGE, GKYL_BC_REFLECT, app->basis_on_dev.confBasis, &app->lower_skin[d],
        &app->lower_ghost[d], app->confBasis.num_basis, app->cdim, app->use_gpu
      );
    } else if (s->lower_bc[d] == GKYL_SPECIES_REFLECT) {
      // the distributions use bc_lo_reflect; the basic updater is not applied there
      bctype_dist = GKYL_BC_COPY;
      bctype_fluid = GKYL_BC_PKPM_MOM_REFLECT;
      s->bc_lo_reflect[d] = pkpm_species_reflect_new(app, s, d, GKYL_LOWER_EDGE);
    } else if (s->lower_bc[d] == GKYL_SPECIES_FIXED_FUNC) {
      bctype_dist = GKYL_BC_FIXED_FUNC;
      bctype_fluid = GKYL_BC_FIXED_FUNC;
    }

    // Distribution function non-periodic lower boundary conditions
    s->bc_lo_dist[d] = gkyl_bc_basic_new(
      d, GKYL_LOWER_EDGE, bctype_dist, app->basis_on_dev.basis, &s->lower_skin_dist[d],
      &s->lower_ghost_dist[d], s->f->ncomp, app->cdim, app->use_gpu
    );
    // Momentum non-periodic lower boundary conditions
    s->bc_lo_fluid[d] = gkyl_bc_basic_new(
      d, GKYL_LOWER_EDGE, bctype_fluid, app->basis_on_dev.confBasis, &app->lower_skin[d],
      &app->lower_ghost[d], s->fluid->ncomp, app->cdim, app->use_gpu
    );

    // Upper BC updater. Copy BCs by default.
    if (s->upper_bc[d] == GKYL_SPECIES_COPY) {
      bctype_dist = GKYL_BC_COPY;
      bctype_fluid = GKYL_BC_COPY;
    } else if (s->upper_bc[d] == GKYL_SPECIES_ABSORB) {
      bctype_dist = GKYL_BC_ABSORB;
      bctype_fluid = GKYL_BC_REFLECT;
      s->bc_is_absorb = true;
      s->bc_up_theta[d] = gkyl_bc_basic_new(
        d, GKYL_UPPER_EDGE, GKYL_BC_REFLECT, app->basis_on_dev.confBasis, &app->upper_skin[d],
        &app->upper_ghost[d], app->confBasis.num_basis, app->cdim, app->use_gpu
      );
    } else if (s->upper_bc[d] == GKYL_SPECIES_REFLECT) {
      bctype_dist = GKYL_BC_COPY;
      bctype_fluid = GKYL_BC_PKPM_MOM_REFLECT;
      s->bc_up_reflect[d] = pkpm_species_reflect_new(app, s, d, GKYL_UPPER_EDGE);
    } else if (s->upper_bc[d] == GKYL_SPECIES_FIXED_FUNC) {
      bctype_dist = GKYL_BC_FIXED_FUNC;
      bctype_fluid = GKYL_BC_FIXED_FUNC;
    }

    // Distribution function non-periodic upper boundary conditions
    s->bc_up_dist[d] = gkyl_bc_basic_new(
      d, GKYL_UPPER_EDGE, bctype_dist, app->basis_on_dev.basis, &s->upper_skin_dist[d],
      &s->upper_ghost_dist[d], s->f->ncomp, app->cdim, app->use_gpu
    );
    // Momentum non-periodic upper boundary conditions
    s->bc_up_fluid[d] = gkyl_bc_basic_new(
      d, GKYL_UPPER_EDGE, bctype_fluid, app->basis_on_dev.confBasis, &app->upper_skin[d],
      &app->upper_ghost[d], s->fluid->ncomp, app->cdim, app->use_gpu
    );
  }
}

void
pkpm_species_apply_ic(gkyl_pkpm_app *app, struct pkpm_species *species, double t0)
{
  int poly_order = app->poly_order;
  // The user provides the normalized distributions [F, G] on the normalized velocity grid
  // (w = v_par/xi; F a probability density in w, G its perpendicular-energy partner so
  // that M0(F) = 1, M1(F) = 0 and M2(F) + 2 M0(G) = 3) and the fluid state
  // [rho, rho ux, rho uy, rho uz, E].
  // Both are sampled at the physical coordinates of the position map (identity: unchanged).
  gkyl_proj_on_basis *proj_dist = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
    .grid = &species->grid,
    .basis = &app->basis,
    .num_quad = 8,
    .num_ret_vals = 2,
    .eval = species->info.init_dist,
    .ctx = species->info.ctx_dist,
    .c2p_func = pkpm_pos_c2p,
    .c2p_func_ctx = app->pos_map,
  });
  gkyl_proj_on_basis *proj_fluid = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
    .grid = &app->grid,
    .basis = &app->confBasis,
    .num_quad = poly_order + 1,
    .num_ret_vals = 5,
    .eval = species->info.init_fluid,
    .ctx = species->info.ctx_fluid,
    .c2p_func = pkpm_pos_c2p,
    .c2p_func_ctx = app->pos_map,
  });

  // run updaters; need to project onto extended range for ease of handling
  // subsequent operations over extended range such as primitive variable computations
  // This is needed to fill the corner cells as the corner cells may not be filled by
  // boundary conditions and we cannot divide by 0 anywhere or the weak divisions will fail
  // (the velocity ghost cell stays empty)
  gkyl_array_clear(species->f_host, 0.0);
  gkyl_proj_on_basis_advance(proj_dist, t0, &species->local_ext_conf, species->f_host);
  gkyl_proj_on_basis_release(proj_dist);
  gkyl_proj_on_basis_advance(proj_fluid, t0, &app->local_ext, species->fluid_host);
  gkyl_proj_on_basis_release(proj_fluid);

  // note: f_host and fluid_host are the same as f and fluid respectively when not on GPUs
  if (app->use_gpu) {
    gkyl_array_copy(species->f, species->f_host);
    gkyl_array_copy(species->fluid, species->fluid_host);
  }

  // Constraints M0(F) = 1 and M2(F) + 2 M0(G) = 3 exactly: nodal rescaling of the
  // projected distributions by their own moments (removes the projection error).
  pkpm_species_moment_calc(
    &species->pkpm_nmoms, species->local_ext_conf, app->local_ext, species->f
  );
  gkyl_array_copy(species->f1, species->f);
  gkyl_dg_calc_pkpm_vars_normalize_dist(
    species->calc_pkpm_vars, &app->local_ext, &species->local_ext, species->pkpm_nmoms.marr,
    species->f1, species->f
  );
  // The evolved quantities carry the position-map Jacobian: J [F, G] and J [rho, M, E]
  // (identity map: no change). The velocity ghost row of f stays empty.
  gkyl_vlasov_position_map_rescale_jacobpos(
    app->pos_map, &app->basis, &species->local_ext_conf, species->f, species->f
  );
  gkyl_vlasov_position_map_rescale_jacobpos_conf(
    app->pos_map, &app->local_ext, species->fluid, species->fluid
  );
  if (app->use_gpu) {
    gkyl_array_copy(species->f_host, species->f);
    gkyl_array_copy(species->fluid_host, species->fluid);
  }

  // we are pre-computing acceleration for now in case is time-independent
  pkpm_species_calc_app_accel(app, species, t0);

  // copy contents of initial conditions into buffer if specific BCs require them
  // *only works in x dimension for now*
  gkyl_bc_basic_buffer_fixed_func(
    species->bc_lo_dist[0], species->bc_buffer_lo_fixed_dist, species->f
  );
  gkyl_bc_basic_buffer_fixed_func(
    species->bc_up_dist[0], species->bc_buffer_up_fixed_dist, species->f
  );
  gkyl_bc_basic_buffer_fixed_func(
    species->bc_lo_fluid[0], species->bc_buffer_lo_fixed_fluid, species->fluid
  );
  gkyl_bc_basic_buffer_fixed_func(
    species->bc_up_fluid[0], species->bc_buffer_up_fixed_fluid, species->fluid
  );
}

void
pkpm_species_calc_app_accel(gkyl_pkpm_app *app, struct pkpm_species *species, double tm)
{
  if (species->has_app_accel) {
    gkyl_proj_on_basis_advance(
      species->app_accel_proj, tm, &app->local_ext, species->app_accel_host
    );
    if (app->use_gpu) {
      // note: app_accel_host is same as app_accel when not on GPUs
      gkyl_array_copy(species->app_accel, species->app_accel_host);
    }
  }
}

// Physical moments [rho, p_par, p_perp, M1] of the factored state over the extended
// range: rho is the trailing block of the fluid array, the others are rho times the
// normalized moments [1, T_par/m, theta, 0] of [F_0/n, G/n] (weak products). The
// normalized moments must have been computed for the same stage.
// Kinetic coupling moments [p_par, p_perp, q_par, q_perp, s, r] = [p M2, p M0_G, p xi M3, p xi M1_G,
// p m1, p xi (m2/2 + n0)] from the raw normalized moments, the nodal scalar pressure p = rho xi^2
// and the thermal speed xi (weak products projected onto p=1): the pressures, first-harmonic stress
// and heat fluxes the fluid system needs. The constraints
// [M0, M1, (M2 + 2 M0_G)/3] are copied out for the diagnostics.
static void
pkpm_species_coupling_moms(gkyl_pkpm_app *app, struct pkpm_species *species)
{
  int nc = app->confBasis.num_basis;
  // [p_par, p_perp, q_par, q_perp] = p [M2, M0_G, xi M3, xi M1_G] and the first-harmonic stress and
  // heat flux s = p m1, r = p xi (m2/2 + n0) (plan 17.14), p = rho xi^2, from the normalized moments
  // [M0, M2, M0_G, M1, M3, M1_G, mu (3), m1 (3), m2 (3), n0 (3)].
  static const int raw[7] = {1, 2, 4, 5, 9, 10, 11}; // M2, M0_G, M3, M1_G, m1 (3)
  for (int c = 0; c < 10; ++c) {
    if (c < 7) {
      gkyl_dg_mul_op_range(
        &app->confBasis, c, species->pkpm_moms, 2, species->pkpm_xi, raw[c],
        species->pkpm_nmoms.marr, &app->local_ext
      );
    } else {
      int a = c - 7;
      gkyl_array_set_offset(species->pkpm_scratch, 0.5, species->pkpm_nmoms.marr, (12 + a) * nc);
      gkyl_array_accumulate_offset(
        species->pkpm_scratch, 1.0, species->pkpm_nmoms.marr, (15 + a) * nc
      );
      gkyl_dg_mul_op_range(
        &app->confBasis, c, species->pkpm_moms, 2, species->pkpm_xi, 0, species->pkpm_scratch,
        &app->local_ext
      );
    }
    if (c == 2 || c == 3 || c >= 7) {
      gkyl_dg_mul_op_range(
        &app->confBasis, 0, species->pkpm_scratch, c, species->pkpm_moms, 0, species->pkpm_xi,
        &app->local_ext
      );
      gkyl_array_set_offset(species->pkpm_moms, 1.0, species->pkpm_scratch, c * nc);
    }
  }
  gkyl_dg_calc_pkpm_vars_constraints(
    species->calc_pkpm_vars, &app->local_ext, nc, species->pkpm_nmoms.marr, species->pkpm_cmoms
  );
  // the constraints of the physical distributions ([1, 0, 1]): the moments of J [F, G] carry J
  gkyl_vlasov_position_map_divide_jacobpos_conf(
    app->pos_map, &app->local_ext, species->pkpm_cmoms->ncomp, species->pkpm_cmoms,
    species->pkpm_cmoms
  );
}

void
pkpm_species_calc_pkpm_vars(
  gkyl_pkpm_app *app, struct pkpm_species *species, const struct gkyl_array *fin,
  const struct gkyl_array *fluidin
)
{
  struct timespec tm = gkyl_wall_clock();
  int nc = app->confBasis.num_basis;

  // Raw normalized moments of [F, G] over the extended range (ghost cells hold valid data
  // after the boundary conditions, the fluid state is mirrored into the ghost cells of
  // absorbing boundaries), the nodal primitive variables u and [xi, xi^2, p] of the fluid
  // state, then the coupling moments. theta = T_perp/T is the zeroth moment of G directly.
  // The evolved fluid and distributions carry the position-map Jacobian J. The primitive
  // variables, the flow and grad(ln n) need the PHYSICAL fluid (a weak gradient of J rho
  // across cells of different J is not J grad rho); the moments of J [F, G] carry J and
  // the coupling moments p (physical) x M_k(J F) are the J-weighted ones the fluid needs.
  gkyl_vlasov_position_map_divide_jacobpos_conf(
    app->pos_map, &app->local_ext, fluidin->ncomp, fluidin, species->fluid_no_J
  );
  pkpm_species_moment_calc(&species->pkpm_nmoms, species->local_ext_conf, app->local_ext, fin);
  gkyl_array_set_offset(species->pkpm_rho, 1.0, species->fluid_no_J, 0);
  // The flow u = M/rho and xi by nodal quotients: for p=1 the weak division, rho u = M weakly,
  // which the implicit field coupling needs to conserve momentum (it forms the new momentum
  // as rho (2 u_half - u)).
  species->num_null_rho += gkyl_dg_calc_pkpm_vars_prim(
    species->calc_pkpm_vars, &app->local_ext, species->fluid_no_J, species->pkpm_u, species->pkpm_xi
  );
  pkpm_species_coupling_moms(app, species);
  // theta = M0(G) of the physical G: M0(J G)/J.
  gkyl_array_set_offset(species->pkpm_theta, 1.0, species->pkpm_nmoms.marr, 2 * nc);
  gkyl_vlasov_position_map_divide_jacobpos_conf(
    app->pos_map, &app->local_ext, nc, species->pkpm_theta, species->pkpm_theta
  );
  // Weak gradient of u for the kinetic sources (bb:grad u); needs u in the ghost cells
  // (at absorbing boundaries that of the mirrored fluid state).
  gkyl_dg_calc_pkpm_vars_grad_u(
    species->calc_pkpm_vars, &app->local, species->pkpm_u, species->pkpm_grad_u
  );
  // grad(ln n) for the density-gradient source: weak gradient of [rho, 0, 0] divided by rho.
  gkyl_array_set_offset(species->pkpm_rho3, 1.0, species->pkpm_rho, 0);
  gkyl_dg_calc_pkpm_vars_grad_u(
    species->calc_pkpm_vars, &app->local, species->pkpm_rho3, species->pkpm_grad_rho
  );
  species->num_null_rho += gkyl_dg_calc_pkpm_vars_grad_ln_n(
    species->calc_pkpm_vars, &app->local, species->pkpm_rho, species->pkpm_grad_rho,
    species->pkpm_grad_ln_n
  );
  // Weak gradient of xi for the forces -w^2 b.grad xi and -2 w b.grad xi G (the second
  // block of pkpm_xi3 receives xi^2 and is never read).
  gkyl_array_set_offset(species->pkpm_xi3, 1.0, species->pkpm_xi, 0);
  gkyl_dg_calc_pkpm_vars_grad_u(
    species->calc_pkpm_vars, &app->local, species->pkpm_xi3, species->pkpm_grad_xi
  );
  // The face fluxes average u, xi and theta across every face. At absorbing walls u and xi of
  // the mirrored fluid ghost are the mirror image of the skin values already; theta, the
  // zeroth moment of the empty ghost G, is replaced by the mirror image of the skin cell, so
  // that the face averages are the interior wall traces.
  if (species->bc_is_absorb) {
    for (int d = 0; d < app->cdim; ++d) {
      if (species->bc_lo_theta[d]) {
        gkyl_bc_basic_advance(
          species->bc_lo_theta[d], species->bc_buffer_fluid, species->pkpm_theta
        );
      }
      if (species->bc_up_theta[d]) {
        gkyl_bc_basic_advance(
          species->bc_up_theta[d], species->bc_buffer_fluid, species->pkpm_theta
        );
      }
    }
    gkyl_comm_array_sync(app->comm, &app->local, &app->local_ext, species->pkpm_theta);
  }
  // The coefficients of the configuration-space flux, once per configuration face: the kinetic
  // flux of every velocity cell and the fluid flux read them.
  gkyl_dg_pkpm_conf_flux_surf_face_coef(
    species->calc_conf_flux, &app->local, app->field->bvar, species->pkpm_u, species->pkpm_xi,
    species->pkpm_theta, species->face_coef
  );

  app->stat.species_pkpm_vars_tm += gkyl_time_diff_now_sec(tm);
}

void
pkpm_species_apply_positivity(
  gkyl_pkpm_app *app, struct pkpm_species *species, struct gkyl_array *fluid
)
{
  if (!app->fluid_positivity) {
    return;
  }
  struct timespec tm = gkyl_wall_clock();
  int num_rho = 0, num_e = 0, num_avg = 0;
  gkyl_dg_calc_pkpm_vars_positivity(
    species->calc_pkpm_vars, &app->local, app->fluid_positivity_eps, fluid, species->pkpm_pos_count,
    &num_rho, &num_e, &num_avg
  );
  species->pos_count_since[0] += num_rho;
  species->pos_count_since[1] += num_e;
  species->pos_count_since[2] += num_avg;
  species->pos_count_since[3] += 1;
  app->stat.species_pkpm_vars_tm += gkyl_time_diff_now_sec(tm);
}

// Face fluxes of the fluid system [rho, M, E] of the stage: the moments of the kinetic
// configuration-space flux of the stage distributions (physical traces of rho, u, xi, b
// on both sides), so the fluid and the distributions see one discrete transport. The
// walls need no special case: an absorbing ghost (F = 0, fluid and theta mirrored) gives the
// interior side's contribution with the interior wall coefficients, and a reflecting ghost
// cancels the mass and energy fluxes.
void
pkpm_species_fluid_flux(
  gkyl_pkpm_app *app, struct pkpm_species *species, const struct gkyl_array *fin,
  const struct gkyl_array *fluidin
)
{
  gkyl_dg_pkpm_fluid_flux_surf_advance(
    species->calc_fluid_flux, &app->local, &app->local_ext, &species->local, &species->local_ext,
    app->field->bvar, species->pkpm_u, species->pkpm_xi, species->face_coef, species->fluid_no_J,
    fin, species->cflrate_fluid, species->fluid_flux_surf
  );
}

// Collisionless configuration-space pass of the kinetic update into rhs (cleared here):
// the HLL streaming fluxes, the exact volume terms and the explicit sources. The closures
// are computed from the moments of this partial right-hand side, after the sources and
// collision operators have been accumulated into it.
static void
pkpm_species_kinetic_conf_pass(
  gkyl_pkpm_app *app, struct pkpm_species *species, const struct gkyl_array *fin,
  struct gkyl_array *rhs
)
{
  gkyl_array_clear(rhs, 0.0);
  gkyl_dg_pkpm_conf_flux_surf_advance(
    species->calc_conf_flux, &app->local, &species->local, &species->local_ext, species->face_coef,
    fin, species->cflrate_f, species->conf_flux_surf
  );
  gkyl_hyper_dg_advance(species->slvr_kinetic_conf, &species->local, fin, species->cflrate_f, rhs);
  // the known first-harmonic w fluxes and sources of the eight blocks (modal, exact)
  gkyl_dg_pkpm_harm_vol_coef(
    species->harm_vol, &app->local, app->field->bvar, app->field->grad_b, species->pkpm_u,
    species->pkpm_grad_u, species->pkpm_xi, species->pkpm_grad_xi, species->pkpm_grad_ln_n,
    species->pkpm_theta, app->field->dtB_over_magB, app->field->tot_em, species->harm_coef,
    species->harm_wnode
  );
  gkyl_dg_pkpm_harm_vol_advance(
    species->harm_vol, &app->local, &species->local, species->harm_coef, fin, species->cflrate_f,
    rhs
  );
}

void
pkpm_species_accumulate_conf_phase(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct gkyl_array *out, double a,
  const struct gkyl_array *conf, const struct gkyl_array *phase
)
{
  int np = app->basis.num_basis;
  for (int b = 0; b < phase->ncomp / np; ++b) {
    gkyl_array_set_offset(s->blk_in, 1.0, phase, b * np);
    gkyl_array_clear(s->blk_out, 0.0);
    gkyl_dg_mul_conf_phase_op_accumulate_range(
      &app->confBasis, &app->basis, s->blk_out, a, conf, s->blk_in, &app->local, &s->local
    );
    gkyl_array_accumulate_offset_range(out, 1.0, s->blk_out, b * np, &s->local);
  }
}

void
pkpm_species_project_transverse(
  gkyl_pkpm_app *app, struct pkpm_species *species, struct gkyl_array *f
)
{
  // Extended ranges: b is formed in the ghost cells too, and the projection commutes with the
  // periodic copy and the reflecting-wall image, so ghost and skin stay consistent.
  gkyl_dg_calc_pkpm_vars_transverse(
    species->calc_pkpm_vars, &app->local_ext, &species->local_ext, app->field->bvar, f
  );
}

void
pkpm_species_harm_exp_advance(
  gkyl_pkpm_app *app, struct pkpm_species *species, int stage, double h, struct gkyl_array *fout
)
{
  // The mixing is frozen with the field and the state of stage 1 (the start of the step); a
  // restarted step freezes it again with the new step.
  if (stage == 1) {
    gkyl_dg_pkpm_harm_exp_coef(
      species->harm_exp, &app->local, app->field->tot_em, species->harm_exp_coef
    );
    gkyl_array_clear(species->harm_exp_act, 0.0);
  }
  gkyl_dg_pkpm_harm_exp_advance(
    species->harm_exp, stage, h, &app->local, &species->local, species->harm_exp_coef, species->f,
    species->harm_exp_nu, species->harm_exp_act, fout
  );
  gkyl_dg_pkpm_harm_exp_restore(
    species->harm_exp, &app->local, &species->local, species->harm_exp_act, species->harm_exp_col,
    fout
  );
}

// Closures from the moments of the partial right-hand side rhs_A (configuration pass):
// c_0 M0 = -M0(rhs_A), (alpha_0, c) from the first and scalar second moments and a_perp from the
// perpendicular rows I_perp int M0_perp dw (pkpm-vars-closure.mac).
static void
pkpm_species_kinetic_closure(
  gkyl_pkpm_app *app, struct pkpm_species *species, const struct gkyl_array *fin,
  const struct gkyl_array *rhs
)
{
  // Closures from the moments of the configuration pass, with the characteristic
  // velocity-edge fluxes inside the 3x3 solve (upwind switch from the closures of the
  // previous stage): the edge fluxes are written for the velocity pass, so that the
  // complete right-hand side keeps the constraints exactly.
  gkyl_array_clear(species->pkpm_rhs_moms, 0.0);
  pkpm_species_moment_calc_into(
    &species->pkpm_nmoms, species->local, app->local, rhs, species->pkpm_rhs_moms
  );
  int num_null = gkyl_dg_calc_pkpm_vars_closure(
    species->calc_pkpm_vars, &app->local, &species->local, &species->grid, species->pkpm_nmoms.marr,
    species->pkpm_rhs_moms, app->field->bvar, species->harm_wnode, species->pkpm_theta, fin,
    species->pkpm_c_0, species->pkpm_alpha_0, species->pkpm_c_h, species->pkpm_a_perp,
    species->vel_flux_surf, species->cflrate_f
  );
  species->num_null_rho += num_null;
  species->num_null_closure += num_null;
  species->closure_cond_min = fmin(
    species->closure_cond_min, gkyl_dg_calc_pkpm_vars_closure_cond_min(species->calc_pkpm_vars)
  );
  // A nulled node means the traces at the velocity-domain edges are large enough to make the
  // closure system singular: the constraints drift in that cell. Say so once, with the proxy.
  if (num_null > 0 && !species->closure_null_warned) {
    species->closure_null_warned = true;
    int rank = 0;
    gkyl_comm_get_rank(app->comm, &rank);
    fprintf(
      stderr,
      "WARNING (rank %d): species %s at t = %g: the closure solve was nulled at %d nodes (min det/det0 so far %g).\n"
      "  The normalized distribution is not negligible at the velocity-space boundary; the constraints drift there.\n"
      "  The velocity extents are likely too small for this species.\n",
      rank, species->info.name, app->tcurr, num_null, species->closure_cond_min
    );
  }
}

void
pkpm_species_calc_diag_vars(
  gkyl_pkpm_app *app, struct pkpm_species *species, const struct gkyl_array *fin,
  const struct gkyl_array *fluidin
)
{
  pkpm_species_calc_pkpm_vars(app, species, fin, fluidin);
  gkyl_array_clear(species->cflrate_f, 0.0);
  gkyl_array_clear(species->cflrate_fluid, 0.0);
  pkpm_species_fluid_flux(app, species, fin, fluidin);
  // fnew is free scratch between time steps.
  pkpm_species_kinetic_conf_pass(app, species, fin, species->fnew);
  pkpm_species_kinetic_closure(app, species, fin, species->fnew);
}

// Compute the RHS for species update, returning maximum stable
// time-step.
double
pkpm_species_rhs(
  gkyl_pkpm_app *app, struct pkpm_species *species, const struct gkyl_array *fin,
  const struct gkyl_array *fluidin, const struct gkyl_array *em, struct gkyl_array *rhs_f,
  struct gkyl_array *rhs_fluid
)
{
  gkyl_array_clear(species->cflrate_f, 0.0);
  gkyl_array_clear(rhs_fluid, 0.0);

  // The fluid face fluxes of the stage (and the wall losses) were computed for every
  // species by the forward-Euler step before any right-hand side.
  struct timespec tm = gkyl_wall_clock();
  // Kinetic equations: the collisionless configuration pass (clears rhs_f), the sources and
  // collisions, the closures from the moments of that partial right-hand side, then the
  // velocity pass accumulating into the same rhs_f.
  pkpm_species_kinetic_conf_pass(app, species, fin, rhs_f);
  // Sources and collisions enter before the closures, which keep the constraints exact
  // against them (the fluid carries the density and energy changes). The self-species LBO
  // conserves the constraints on its own (boundary corrections, zero-flux w edges); it sits
  // here so that a cross-species energy exchange reaches the thermal-speed closure.
  pkpm_species_source_kinetic_rhs(app, species, &species->src, fin, rhs_f);
  if (species->collision_id == GKYL_BGK_COLLISIONS) {
    pkpm_species_bgk_kinetic_rhs(app, species, &species->bgk, fin, rhs_f);
  }
  if (species->collision_id == GKYL_LBO_COLLISIONS) {
    pkpm_species_lbo_rhs(app, species, &species->lbo, fin, rhs_f);
  }
  pkpm_species_kinetic_closure(app, species, fin, rhs_f);
  gkyl_dg_pkpm_vel_flux_surf_advance(
    species->calc_vel_flux, &app->local, &species->local, species->harm_wnode, species->pkpm_theta,
    app->field->bvar, species->pkpm_alpha_0, species->pkpm_c_h, fin, species->vel_flux_surf
  );
  gkyl_hyper_dg_advance(species->slvr_kinetic_vel, &species->local, fin, species->cflrate_f, rhs_f);
  species->kinetic_tm += gkyl_time_diff_now_sec(tm);

  // Momentum equation: volume and surface terms from the stored face fluxes.
  tm = gkyl_wall_clock();
  gkyl_hyper_dg_advance(
    species->slvr_fluid, &app->local, fluidin, species->cflrate_fluid, rhs_fluid
  );
  pkpm_species_source_fluid_rhs(app, species, &species->src, rhs_fluid);
  if (species->collision_id == GKYL_BGK_COLLISIONS) {
    pkpm_species_bgk_fluid_rhs(app, species, &species->bgk, rhs_fluid);
  }
  species->fluid_tm += gkyl_time_diff_now_sec(tm);

  if (species->has_diffusion) {
    gkyl_dg_pkpm_mom_diffusion_advance(
      species->mom_diff, &app->local, fluidin, species->cflrate_fluid, rhs_fluid
    );
  }

  // If PKPM update is fully explicit, include the fluid-EM coupling in fluid update
  if (app->use_explicit_source) {
    double qbym = species->info.charge / species->info.mass;
    gkyl_array_set(species->qmem, qbym, em);

    // Accumulate applied acceleration and/or q/m*(external electromagnetic)
    // fields onto qmem to get the total acceleration
    if (species->has_app_accel) {
      gkyl_array_accumulate_range(species->qmem, 1.0, species->app_accel, &app->local);
    }
    if (app->field->has_ext_em) {
      gkyl_array_accumulate_range(species->qmem, qbym, app->field->ext_em, &app->local);
    }

    gkyl_dg_calc_pkpm_vars_source(
      species->calc_pkpm_vars, &app->local, species->qmem, fluidin, rhs_fluid
    );
  }

  app->stat.n_species_omega_cfl += 1;
  tm = gkyl_wall_clock();
  gkyl_array_reduce_range(species->omegaCfl_ptr_dist, species->cflrate_f, GKYL_MAX, &species->local);
  gkyl_array_reduce_range(
    species->omegaCfl_ptr_fluid, species->cflrate_fluid, GKYL_MAX, &app->local
  );

  double omegaCfl_ho_dist[1];
  double omegaCfl_ho_fluid[1];
  if (app->use_gpu) {
    gkyl_cu_memcpy(omegaCfl_ho_dist, species->omegaCfl_ptr_dist, sizeof(double), GKYL_CU_MEMCPY_D2H);
    gkyl_cu_memcpy(
      omegaCfl_ho_fluid, species->omegaCfl_ptr_fluid, sizeof(double), GKYL_CU_MEMCPY_D2H
    );
  } else {
    omegaCfl_ho_dist[0] = species->omegaCfl_ptr_dist[0];
    omegaCfl_ho_fluid[0] = species->omegaCfl_ptr_fluid[0];
  }
  double omegaCfl = fmax(omegaCfl_ho_dist[0], omegaCfl_ho_fluid[0]);

  app->stat.species_omega_cfl_tm += gkyl_time_diff_now_sec(tm);

  return app->cfl / omegaCfl;
}

// Determine which directions are periodic and which directions are not periodic,
// and then apply boundary conditions for distribution function
void
pkpm_species_apply_bc(gkyl_pkpm_app *app, const struct pkpm_species *species, struct gkyl_array *f)
{
  struct timespec wst = gkyl_wall_clock();

  int num_periodic_dir = app->num_periodic_dir, cdim = app->cdim;
  gkyl_comm_array_per_sync(
    species->comm, &species->local, &species->local_ext_conf, num_periodic_dir, app->periodic_dirs,
    f
  );

  int is_np_bc[3] = {1, 1, 1}; // flags to indicate if direction is periodic
  for (int d = 0; d < num_periodic_dir; ++d) {
    is_np_bc[app->periodic_dirs[d]] = 0;
  }

  for (int d = 0; d < cdim; ++d) {
    if (is_np_bc[d]) {
      switch (species->lower_bc[d]) {
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(species->bc_lo_dist[d], species->bc_buffer_dist, f);
          break;
        case GKYL_SPECIES_REFLECT:
          gkyl_bc_pkpm_species_reflect_advance(species->bc_lo_reflect[d], app->field->bvar, f);
          break;
        case GKYL_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(species->bc_lo_dist[d], species->bc_buffer_lo_fixed_dist, f);
          break;
        case GKYL_SPECIES_NO_SLIP:
        case GKYL_SPECIES_WEDGE:
          assert(false);
          break;
        default:
          break;
      }

      switch (species->upper_bc[d]) {
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(species->bc_up_dist[d], species->bc_buffer_dist, f);
          break;
        case GKYL_SPECIES_REFLECT:
          gkyl_bc_pkpm_species_reflect_advance(species->bc_up_reflect[d], app->field->bvar, f);
          break;
        case GKYL_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(species->bc_up_dist[d], species->bc_buffer_up_fixed_dist, f);
          break;
        case GKYL_SPECIES_NO_SLIP:
        case GKYL_SPECIES_WEDGE:
          assert(false);
          break;
        default:
          break;
      }
    }
  }

  gkyl_comm_array_sync(species->comm, &species->local, &species->local_ext_conf, f);

  app->stat.species_bc_tm += gkyl_time_diff_now_sec(wst);
}

// Determine which directions are periodic and which directions are not periodic,
// and then apply boundary conditions for distribution function
void
pkpm_fluid_species_apply_bc(
  gkyl_pkpm_app *app, const struct pkpm_species *species, struct gkyl_array *fluid
)
{
  struct timespec wst = gkyl_wall_clock();

  int num_periodic_dir = app->num_periodic_dir, cdim = app->cdim;
  gkyl_comm_array_per_sync(
    app->comm, &app->local, &app->local_ext, num_periodic_dir, app->periodic_dirs, fluid
  );

  int is_np_bc[3] = {1, 1, 1}; // flags to indicate if direction is periodic
  for (int d = 0; d < num_periodic_dir; ++d) {
    is_np_bc[app->periodic_dirs[d]] = 0;
  }

  for (int d = 0; d < cdim; ++d) {
    if (is_np_bc[d]) {
      switch (species->lower_bc[d]) {
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(species->bc_lo_fluid[d], species->bc_buffer_fluid, fluid);
          break;
        case GKYL_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(species->bc_lo_fluid[d], species->bc_buffer_lo_fixed_fluid, fluid);
          break;
        case GKYL_SPECIES_NO_SLIP:
        case GKYL_SPECIES_WEDGE:
          assert(false);
          break;
        default:
          break;
      }

      switch (species->upper_bc[d]) {
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(species->bc_up_fluid[d], species->bc_buffer_fluid, fluid);
          break;
        case GKYL_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(species->bc_up_fluid[d], species->bc_buffer_up_fixed_fluid, fluid);
          break;
        case GKYL_SPECIES_NO_SLIP:
        case GKYL_SPECIES_WEDGE:
          assert(false);
          break;
        default:
          break;
      }
    }
  }

  gkyl_comm_array_sync(app->comm, &app->local, &app->local_ext, fluid);

  app->stat.species_bc_tm += gkyl_time_diff_now_sec(wst);
}

void
pkpm_species_calc_L2(gkyl_pkpm_app *app, double tm, const struct pkpm_species *species)
{
  // The evolved f carries the position-map Jacobian: int f^2 dx = int (f/J) f dxi.
  gkyl_vlasov_position_map_divide_jacobpos(
    app->pos_map, &app->basis, &species->local, species->f, species->f_no_J
  );
  gkyl_dg_calc_prod_op_range(
    app->basis, 0, species->L2_f, 0, species->f_no_J, 0, species->f, species->local
  );
  gkyl_array_scale_range(species->L2_f, species->grid.cellVolume, &species->local);

  double L2[1] = {0.0};
  if (app->use_gpu) {
    gkyl_array_reduce_range(species->red_L2_f, species->L2_f, GKYL_SUM, &species->local);
    gkyl_cu_memcpy(L2, species->red_L2_f, sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    gkyl_array_reduce_range(L2, species->L2_f, GKYL_SUM, &species->local);
  }
  double L2_global[1] = {0.0};
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, L2, L2_global);

  gkyl_dynvec_append(species->integ_L2_f, tm, L2_global);
}

void
pkpm_species_coll_tm(gkyl_pkpm_app *app)
{
  for (int i = 0; i < app->num_species; ++i) {
    if (app->species[i].collision_id == GKYL_LBO_COLLISIONS) {
      struct gkyl_dg_updater_lbo_pkpm_tm tm =
        gkyl_dg_updater_lbo_pkpm_get_tm(app->species[i].lbo.coll_slvr);
      app->stat.species_lbo_coll_diff_tm[i] = tm.diff_tm;
      app->stat.species_lbo_coll_drag_tm[i] = tm.drag_tm;
    }
  }
}

void
pkpm_species_tm(gkyl_pkpm_app *app)
{
  app->stat.species_rhs_tm = 0.0;
  app->stat.fluid_species_rhs_tm = 0.0;
  for (int i = 0; i < app->num_species; ++i) {
    app->stat.species_rhs_tm += app->species[i].kinetic_tm;
    app->stat.fluid_species_rhs_tm += app->species[i].fluid_tm;
  }
}

// release resources for PKPM species
void
pkpm_species_release(const gkyl_pkpm_app *app, const struct pkpm_species *s)
{
  // release various arrays
  gkyl_array_release(s->f);
  gkyl_array_release(s->f_no_J);
  gkyl_array_release(s->f1);
  gkyl_array_release(s->fnew);
  gkyl_array_release(s->blk_in);
  gkyl_array_release(s->blk_out);
  gkyl_array_release(s->fluid);
  gkyl_array_release(s->fluid1);
  gkyl_array_release(s->fluidnew);
  gkyl_array_release(s->fluid_dup);

  gkyl_array_release(s->cflrate_f);
  gkyl_array_release(s->cflrate_fluid);

  gkyl_array_release(s->bc_buffer_dist);
  gkyl_array_release(s->bc_buffer_lo_fixed_dist);
  gkyl_array_release(s->bc_buffer_up_fixed_dist);
  gkyl_array_release(s->bc_buffer_fluid);
  gkyl_array_release(s->bc_buffer_lo_fixed_fluid);
  gkyl_array_release(s->bc_buffer_up_fixed_fluid);

  gkyl_comm_release(s->comm);

  if (app->use_gpu) {
    gkyl_array_release(s->f_host);
    gkyl_array_release(s->f_no_J_host);
    gkyl_array_release(s->fluid_host);
  }

  gkyl_array_release(s->qmem);

  // release moment data
  pkpm_species_moment_release(app, &s->pkpm_moms_diag);
  pkpm_species_moment_release(app, &s->pkpm_nmoms);

  // release pkpm arrays data and updaters
  gkyl_array_release(s->pkpm_u);
  gkyl_array_release(s->pkpm_grad_u);
  gkyl_array_release(s->pkpm_theta);
  gkyl_array_release(s->pkpm_alpha_0);
  gkyl_array_release(s->pkpm_c_h);
  gkyl_array_release(s->pkpm_xi);
  gkyl_array_release(s->pkpm_xi3);
  gkyl_array_release(s->pkpm_grad_xi);
  gkyl_array_release(s->pkpm_scratch);
  gkyl_array_release(s->pkpm_moms);
  gkyl_array_release(s->pkpm_moms_no_J);
  gkyl_array_release(s->fluid_no_J);
  gkyl_array_release(s->pkpm_cmoms);
  if (app->use_gpu) {
    gkyl_array_release(s->pkpm_cmoms_host);
  }
  gkyl_array_release(s->pkpm_rhs_moms);
  gkyl_array_release(s->pkpm_p_ij);
  gkyl_array_release(s->conf_flux_surf);
  gkyl_array_release(s->face_coef);
  gkyl_array_release(s->vel_flux_surf);
  gkyl_array_release(s->fluid_flux_surf);
  gkyl_array_release(s->pkpm_rho);
  gkyl_array_release(s->pkpm_rho3);
  gkyl_array_release(s->pkpm_mom3);
  if (s->has_diffusion) {
    gkyl_dg_pkpm_mom_diffusion_release(s->mom_diff);
  }
  gkyl_array_release(s->pkpm_grad_rho);
  gkyl_array_release(s->pkpm_grad_ln_n);
  gkyl_array_release(s->pkpm_c_0);
  gkyl_array_release(s->pkpm_a_perp);
  gkyl_array_release(s->integ_pkpm_mom);
  gkyl_array_release(s->fluid_io);
  gkyl_array_release(s->pkpm_vars_io);
  gkyl_array_release(s->pkpm_pos_count);
  gkyl_array_release(s->pkpm_pos_io);
  if (app->use_gpu) {
    gkyl_array_release(s->pkpm_pos_io_host);
  }
  gkyl_dg_calc_pkpm_vars_release(s->calc_pkpm_vars);
  gkyl_dg_pkpm_conf_flux_surf_release(s->calc_conf_flux);
  gkyl_dg_pkpm_vel_flux_surf_release(s->calc_vel_flux);
  gkyl_dg_pkpm_harm_vol_release(s->harm_vol);
  gkyl_array_release(s->harm_coef);
  gkyl_array_release(s->harm_wnode);
  gkyl_dg_pkpm_harm_exp_release(s->harm_exp);
  gkyl_array_release(s->harm_exp_coef);
  gkyl_array_release(s->harm_exp_nu);
  gkyl_array_release(s->harm_exp_act);
  gkyl_array_release(s->harm_exp_col);
  gkyl_dg_pkpm_fluid_flux_surf_release(s->calc_fluid_flux);
  gkyl_hyper_dg_release(s->slvr_kinetic_conf);
  gkyl_hyper_dg_release(s->slvr_kinetic_vel);
  gkyl_hyper_dg_release(s->slvr_fluid);
  gkyl_dg_eqn_release(s->eqn_kinetic_conf);
  gkyl_dg_eqn_release(s->eqn_kinetic_vel);
  gkyl_dg_eqn_release(s->eqn_fluid);

  gkyl_array_release(s->L2_f);
  gkyl_dynvec_release(s->integ_L2_f);
  gkyl_dynvec_release(s->integ_diag);
  gkyl_dynvec_release(s->pos_stat);
  gkyl_dynvec_release(s->bflux_diag);
  if (s->fluid_flux_surf_host) {
    gkyl_array_release(s->fluid_flux_surf_host);
  }
  pkpm_species_source_release(app, &s->src);

  gkyl_array_release(s->app_accel);
  if (s->has_app_accel) {
    if (app->use_gpu) {
      gkyl_array_release(s->app_accel_host);
    }
    gkyl_proj_on_basis_release(s->app_accel_proj);
  }

  if (s->collision_id == GKYL_BGK_COLLISIONS) {
    pkpm_species_bgk_release(app, &s->bgk);
  }
  if (s->collision_id == GKYL_LBO_COLLISIONS) {
    pkpm_species_lbo_release(app, &s->lbo);
  }

  // Copy BCs are allocated by default. Need to free.
  for (int d = 0; d < app->cdim; ++d) {
    gkyl_bc_basic_release(s->bc_lo_dist[d]);
    gkyl_bc_basic_release(s->bc_up_dist[d]);
    if (s->bc_lo_reflect[d]) {
      gkyl_bc_pkpm_species_reflect_release(s->bc_lo_reflect[d]);
    }
    if (s->bc_up_reflect[d]) {
      gkyl_bc_pkpm_species_reflect_release(s->bc_up_reflect[d]);
    }
    gkyl_bc_basic_release(s->bc_lo_fluid[d]);
    gkyl_bc_basic_release(s->bc_up_fluid[d]);
    if (s->bc_lo_theta[d]) {
      gkyl_bc_basic_release(s->bc_lo_theta[d]);
    }
    if (s->bc_up_theta[d]) {
      gkyl_bc_basic_release(s->bc_up_theta[d]);
    }
  }

  if (app->use_gpu) {
    gkyl_cu_free(s->omegaCfl_ptr_dist);
    gkyl_cu_free(s->omegaCfl_ptr_fluid);
    gkyl_cu_free(s->red_L2_f);
    gkyl_cu_free(s->red_integ_diag);
  } else {
    gkyl_free(s->omegaCfl_ptr_dist);
    gkyl_free(s->omegaCfl_ptr_fluid);
  }
}
