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
    .phase_grid = &s->kinetic->grid,
    .phase_basis = &app->basis,
    .phase_basis_on_dev = app->basis_on_dev.basis,
    .conf_num_basis = app->confBasis.num_basis,
    .num_blocks = PKPM_DIST_NUM_BLOCKS,
    .skin_r = edge == GKYL_LOWER_EDGE ? &s->kinetic->lower_skin[dir] : &s->kinetic->upper_skin[dir],
    .ghost_r = edge == GKYL_LOWER_EDGE ? &s->kinetic->lower_ghost[dir] :
                                         &s->kinetic->upper_ghost[dir],
    .phase_range_ext = &s->kinetic->local_ext,
    .conf_range_ext = &app->local_ext,
    .use_gpu = app->use_gpu,
  });
}

// The basic boundary condition of the distributions for a kinetic BC type. A reflecting wall
// uses the field-aligned reflection updater instead; the basic updater is not applied there.
static enum gkyl_bc_basic_type
pkpm_species_bc_basic_type_dist(enum gkyl_species_bc_type bc_type)
{
  switch (bc_type) {
    case GKYL_SPECIES_ABSORB:
      return GKYL_BC_ABSORB;
    case GKYL_SPECIES_FIXED_FUNC:
      return GKYL_BC_FIXED_FUNC;
    default:
      return GKYL_BC_COPY;
  }
}

// The basic boundary condition of the fluid state for a fluid BC type: at an absorbing wall the
// mirror image of the skin cell (even in the wall direction), at a reflecting wall the mirror
// image with the momentum reflected.
static enum gkyl_bc_basic_type
pkpm_species_bc_basic_type_fluid(enum gkyl_species_bc_type bc_type)
{
  switch (bc_type) {
    case GKYL_SPECIES_ABSORB:
      return GKYL_BC_REFLECT;
    case GKYL_SPECIES_REFLECT:
      return GKYL_BC_PKPM_MOM_REFLECT;
    case GKYL_SPECIES_FIXED_FUNC:
      return GKYL_BC_FIXED_FUNC;
    default:
      return GKYL_BC_COPY;
  }
}

void
pkpm_species_alloc(const struct gkyl_pkpm_species *inp, struct pkpm_species *s)
{
  s->type = inp->type;
  strcpy(s->name, inp->name);
  s->charge = inp->charge;
  s->mass = inp->mass;

  // Each aspect holds its input block and a copy of the identity.
  s->kinetic = gkyl_malloc(sizeof(struct pkpm_kinetic_species));
  s->kinetic->info = inp->kinetic;
  strcpy(s->kinetic->name, inp->name);
  s->kinetic->charge = inp->charge;
  s->kinetic->mass = inp->mass;

  s->fluid = gkyl_malloc(sizeof(struct pkpm_fluid_species));
  s->fluid->info = inp->fluid;
  strcpy(s->fluid->name, inp->name);
  s->fluid->charge = inp->charge;
  s->fluid->mass = inp->mass;
}

void
pkpm_species_init(struct gkyl_pkpm *pkpm, struct gkyl_pkpm_app *app, struct pkpm_species *s)
{
  int cdim = app->cdim, vdim = app->vdim;
  int pdim = cdim + vdim;

  // A PKPM species owns both blocks; the distributions start from one function of eight blocks.
  assert(s->type == GKYL_SPECIES_PKPM);
  assert(s->kinetic->info.num_init == 1);
  assert(s->kinetic->info.projection[0].proj_id == GKYL_PROJ_FUNC);
  assert(s->kinetic->info.projection[0].func && s->fluid->info.init);
  // defaults of the number of face quadrature nodes and of the fluid positivity floor
  s->kinetic->num_face_nodes =
    s->kinetic->info.num_face_nodes > 0 ? s->kinetic->info.num_face_nodes : 5;
  assert(s->kinetic->num_face_nodes >= 3);
  s->fluid->positivity_eps = s->fluid->info.positivity_eps > 0.0 ? s->fluid->info.positivity_eps :
                                                                   1.0e-2;

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
    // full phase space grid (the velocity coordinate is w = v_par/xi, normalized)
    cells[cdim + d] = s->kinetic->info.cells[d];
    lower[cdim + d] = s->kinetic->info.lower[d];
    upper[cdim + d] = s->kinetic->info.upper[d];
    ghost[cdim + d] = 0; // no ghost-cells in velocity space

    // only velocity space
    cells_vel[d] = s->kinetic->info.cells[d];
    lower_vel[d] = s->kinetic->info.lower[d];
    upper_vel[d] = s->kinetic->info.upper[d];
    ghost_vel[d] = 0; // no ghost-cells in velocity space
  }
  // full phase space grid
  gkyl_rect_grid_init(&s->kinetic->grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&s->kinetic->grid, ghost, &s->kinetic->global_ext, &s->kinetic->global);

  // velocity space grid
  gkyl_rect_grid_init(&s->kinetic->grid_vel, vdim, lower_vel, upper_vel, cells_vel);
  gkyl_create_grid_ranges(
    &s->kinetic->grid_vel, ghost_vel, &s->kinetic->local_ext_vel, &s->kinetic->local_vel
  );

  // phase-space communicator
  s->kinetic->comm = gkyl_comm_extend_comm(app->comm, &s->kinetic->local_vel);

  // create local and local_ext from app local range: local = conf-local X local_vel.
  // The extended range has the configuration ghost cells and ONE velocity ghost cell
  // above w_max: the empty state beyond the characteristic velocity edge, whose face slot
  // holds the upper-edge flux (the lower edge is the lowest cell's own slot).
  // local_ext_conf is extended in configuration space only: the range of everything that
  // integrates, projects, syncs or reflects the distribution.
  struct gkyl_range local;
  gkyl_range_ten_prod(&local, &app->local, &s->kinetic->local_vel);
  int elo[GKYL_MAX_DIM] = {0}, eup[GKYL_MAX_DIM] = {0};
  for (int d = 0; d < cdim; ++d) {
    elo[d] = eup[d] = 1;
  }
  eup[cdim] = 1;
  gkyl_range_extend(&s->kinetic->local_ext, &local, elo, eup);
  gkyl_sub_range_init(&s->kinetic->local, &s->kinetic->local_ext, local.lower, local.upper);
  int lower_conf[GKYL_MAX_DIM], upper_conf[GKYL_MAX_DIM];
  for (int d = 0; d < pdim; ++d) {
    lower_conf[d] = d < cdim ? s->kinetic->local_ext.lower[d] : local.lower[d];
    upper_conf[d] = d < cdim ? s->kinetic->local_ext.upper[d] : local.upper[d];
  }
  gkyl_sub_range_init(&s->kinetic->local_ext_conf, &s->kinetic->local_ext, lower_conf, upper_conf);

  // allocate distribution function arrays: the normalized distributions
  // [F, G, M0_perp (x, y, z), N_perp (x, y, z)] (enum gkyl_pkpm_dist_block)
  int nc_dist = PKPM_DIST_NUM_BLOCKS * app->basis.num_basis;
  s->kinetic->f = mkarr(app->use_gpu, nc_dist, s->kinetic->local_ext.volume);
  s->kinetic->f1 = mkarr(app->use_gpu, nc_dist, s->kinetic->local_ext.volume);
  // physical distributions (position-map Jacobian divided out) for the output and the L2 norm
  s->kinetic->f_no_J = mkarr(app->use_gpu, nc_dist, s->kinetic->local_ext.volume);
  s->kinetic->fnew = mkarr(app->use_gpu, nc_dist, s->kinetic->local_ext.volume);
  gkyl_array_clear(s->kinetic->f, 0.0);
  gkyl_array_clear(s->kinetic->f1, 0.0);
  gkyl_array_clear(s->kinetic->f_no_J, 0.0);
  gkyl_array_clear(s->kinetic->fnew, 0.0);
  s->kinetic->blk_in = mkarr(app->use_gpu, app->basis.num_basis, s->kinetic->local_ext.volume);
  s->kinetic->blk_out = mkarr(app->use_gpu, app->basis.num_basis, s->kinetic->local_ext.volume);
  // allocate fluid arrays [rho, rho ux, rho uy, rho uz, E]
  s->fluid->fluid = mkarr(app->use_gpu, 5 * app->confBasis.num_basis, app->local_ext.volume);
  s->fluid->fluid1 = mkarr(app->use_gpu, 5 * app->confBasis.num_basis, app->local_ext.volume);
  s->fluid->fluidnew = mkarr(app->use_gpu, 5 * app->confBasis.num_basis, app->local_ext.volume);

  // create host arrays if on GPUs so initialization occurs host-side
  s->kinetic->f_host = s->kinetic->f;
  s->kinetic->f_no_J_host = s->kinetic->f_no_J;
  s->fluid->fluid_host = s->fluid->fluid;
  if (app->use_gpu) {
    s->kinetic->f_host = mkarr(false, nc_dist, s->kinetic->local_ext.volume);
    s->kinetic->f_no_J_host = mkarr(false, nc_dist, s->kinetic->local_ext.volume);
    s->fluid->fluid_host = mkarr(false, 5 * app->confBasis.num_basis, app->local_ext.volume);
  }

  // allocate cflrate for distribution function and fluid updates (scalar arrays)
  s->kinetic->cflrate = mkarr(app->use_gpu, 1, s->kinetic->local_ext.volume);
  s->fluid->cflrate = mkarr(app->use_gpu, 1, app->local_ext.volume);

  if (app->use_gpu) {
    s->kinetic->omega_cfl = gkyl_cu_malloc(sizeof(double));
    s->fluid->omegaCfl_ptr = gkyl_cu_malloc(sizeof(double));
  } else {
    s->kinetic->omega_cfl = gkyl_malloc(sizeof(double));
    s->fluid->omegaCfl_ptr = gkyl_malloc(sizeof(double));
  }

  // q/m (E, B) and the applied acceleration for the explicit Lorentz force
  s->fluid->qmem = mkarr(app->use_gpu, 8 * app->confBasis.num_basis, app->local_ext.volume);

  int nc = app->confBasis.num_basis;
  // pkpm moments for diagnostics [rho, M1, p_par, p_perp, q_par, q_perp, r_parpar, r_parperp]
  pkpm_species_moment_init(app, s, &s->kinetic->pkpm_moms_diag, true);
  // raw moments of the normalized distributions [F, G, M0_perp, N_perp] (mass 1), the kinetic
  // coupling moments [p_par, p_perp, q_par, q_perp, s (3), r (3)] formed from them and the
  // constraints
  pkpm_species_moment_init_mass(app, s, &s->kinetic->pkpm_nmoms, false, 1.0);
  s->fluid->pkpm_moms = mkarr(app->use_gpu, 10 * nc, app->local_ext.volume);
  s->fluid->pkpm_moms_no_J = mkarr(app->use_gpu, 10 * nc, app->local_ext.volume);
  s->fluid->fluid_no_J = mkarr(app->use_gpu, 5 * nc, app->local_ext.volume);
  s->fluid->fluid_no_J_host = s->fluid->fluid_no_J;
  if (app->use_gpu) {
    s->fluid->fluid_no_J_host = mkarr(false, 5 * nc, app->local_ext.volume);
  }
  gkyl_array_clear(s->fluid->pkpm_moms, 0.0);
  s->kinetic->pkpm_cmoms = mkarr(app->use_gpu, 6 * nc, app->local_ext.volume);
  s->kinetic->pkpm_cmoms_host = s->kinetic->pkpm_cmoms;
  if (app->use_gpu) {
    s->kinetic->pkpm_cmoms_host = mkarr(false, 6 * nc, app->local_ext.volume);
  }

  // Configuration-space variables of the update (p=1).
  s->fluid->pkpm_u = mkarr(app->use_gpu, 3 * nc, app->local_ext.volume);
  s->kinetic->pkpm_u = mkarr(app->use_gpu, 3 * nc, app->local_ext.volume);
  s->kinetic->pkpm_grad_u = mkarr(app->use_gpu, cdim * 3 * nc, app->local_ext.volume);
  gkyl_array_clear(s->kinetic->pkpm_grad_u, 0.0);
  s->kinetic->pkpm_theta = mkarr(app->use_gpu, nc, app->local_ext.volume);
  s->fluid->pkpm_rho = mkarr(app->use_gpu, nc, app->local_ext.volume);
  s->fluid->pkpm_rho3 = mkarr(app->use_gpu, 3 * nc, app->local_ext.volume);
  gkyl_array_clear(s->fluid->pkpm_rho3, 0.0);
  s->fluid->pkpm_mom3 = mkarr(app->use_gpu, 3 * nc, app->local_ext.volume);
  s->fluid->pkpm_grad_rho = mkarr(app->use_gpu, cdim * 3 * nc, app->local_ext.volume);
  gkyl_array_clear(s->fluid->pkpm_grad_rho, 0.0);
  s->kinetic->pkpm_grad_ln_n = mkarr(app->use_gpu, cdim * nc, app->local_ext.volume);
  gkyl_array_clear(s->kinetic->pkpm_grad_ln_n, 0.0);
  s->kinetic->pkpm_c_0 = mkarr(app->use_gpu, nc, app->local_ext.volume);
  gkyl_array_clear(s->kinetic->pkpm_c_0, 0.0);
  s->kinetic->pkpm_alpha_0 = mkarr(app->use_gpu, nc, app->local_ext.volume);
  gkyl_array_clear(s->kinetic->pkpm_alpha_0, 0.0);
  s->kinetic->pkpm_c_h = mkarr(app->use_gpu, nc, app->local_ext.volume);
  gkyl_array_clear(s->kinetic->pkpm_c_h, 0.0);
  s->kinetic->pkpm_a_perp = mkarr(app->use_gpu, 3 * nc, app->local_ext.volume);
  gkyl_array_clear(s->kinetic->pkpm_a_perp, 0.0);
  s->kinetic->pkpm_closure_step = mkarr(app->use_gpu, 2 * nc, app->local_ext.volume);
  gkyl_array_clear(s->kinetic->pkpm_closure_step, 0.0);
  s->kinetic->pkpm_closure_step_host = s->kinetic->pkpm_closure_step;
  if (app->use_gpu) {
    s->kinetic->pkpm_closure_step_host = mkarr(false, 2 * nc, app->local_ext.volume);
  }
  s->fluid->pkpm_xi = mkarr(app->use_gpu, 3 * nc, app->local_ext.volume);
  gkyl_array_clear(s->fluid->pkpm_xi, 0.0);
  s->kinetic->pkpm_xi = mkarr(app->use_gpu, 3 * nc, app->local_ext.volume);
  gkyl_array_clear(s->kinetic->pkpm_xi, 0.0);
  s->fluid->pkpm_xi3 = mkarr(app->use_gpu, 3 * nc, app->local_ext.volume);
  gkyl_array_clear(s->fluid->pkpm_xi3, 0.0);
  s->kinetic->pkpm_grad_xi = mkarr(app->use_gpu, cdim * 3 * nc, app->local_ext.volume);
  gkyl_array_clear(s->kinetic->pkpm_grad_xi, 0.0);
  s->fluid->pkpm_scratch = mkarr(app->use_gpu, nc, app->local_ext.volume);
  s->kinetic->pkpm_rhs_moms = mkarr(
    app->use_gpu, gkyl_dg_updater_moment_pkpm_num_mom(s->kinetic->pkpm_nmoms.mcalc) * nc,
    app->local_ext.volume
  );
  gkyl_array_clear(s->kinetic->pkpm_theta, 0.0);
  s->kinetic->num_null_closure = 0;
  s->kinetic->closure_cond_min = 1.0;
  s->kinetic->closure_null_warned = false;

  // Updaters for the configuration-space variables, one per aspect: the fluid one forms the
  // primitive variables, their gradients, the limiter, the sources and the diagnostics, the
  // kinetic one the closures, the constraints and the transverse projection. Every method takes
  // the range it acts on.
  struct gkyl_dg_calc_pkpm_vars_inp calc_pkpm_vars_inp = {
    .conf_grid = &app->grid,
    .conf_basis = &app->confBasis,
    .pos_map = app->pos_map,
    .use_gpu = app->use_gpu,
  };
  s->fluid->calc_pkpm_vars = gkyl_dg_calc_pkpm_vars_inew(&calc_pkpm_vars_inp);
  s->kinetic->calc_pkpm_vars = gkyl_dg_calc_pkpm_vars_inew(&calc_pkpm_vars_inp);

  // array for storing integrated fluid variables in each cell
  // integ_mom : integral[rho, rhoux, rhouy, rhouz, rhoux^2, rhouy^2, rhouz^2, p_par, p_perp]
  s->fluid->integ_mom = mkarr(app->use_gpu, 9, app->local_ext.volume);
  // arrays for I/O, pkpm_prim_io and pkpm_vars_io
  // pkpm_prim_io : [P_xx, P_xy, P_xz, P_yy, P_yz, P_zz, q_x, q_y, q_z]
  // pkpm_vars_io : [ux, uy, uz, T_perp/m, a_P, div(b), bb:grad(u), T_perp/m div(b), xi, c]
  s->fluid->pkpm_prim_io = mkarr(app->use_gpu, 9 * nc, app->local_ext.volume);
  s->kinetic->pkpm_vars_io = mkarr(app->use_gpu, 10 * nc, app->local_ext.volume);
  s->fluid->pkpm_prim_io_host = s->fluid->pkpm_prim_io;
  s->kinetic->pkpm_vars_io_host = s->kinetic->pkpm_vars_io;
  if (app->use_gpu) {
    s->fluid->pkpm_prim_io_host = mkarr(false, 9 * nc, app->local_ext.volume);
    s->kinetic->pkpm_vars_io_host = mkarr(false, 10 * nc, app->local_ext.volume);
  }
  // Record of the positivity limiter: per cell [density scalings, energy scalings,
  // inadmissible averages] accumulated over the run, and its p=1 image for IO.
  s->fluid->pkpm_pos_count = mkarr(app->use_gpu, PKPM_POS_NUM_COUNT, app->local_ext.volume);
  gkyl_array_clear(s->fluid->pkpm_pos_count, 0.0);
  s->fluid->pkpm_pos_io = mkarr(app->use_gpu, PKPM_POS_NUM_COUNT * nc, app->local_ext.volume);
  s->fluid->pkpm_pos_io_host = s->fluid->pkpm_pos_io;
  if (app->use_gpu) {
    s->fluid->pkpm_pos_io_host = mkarr(false, PKPM_POS_NUM_COUNT * nc, app->local_ext.volume);
  }

  // The configuration-space face flux updater and the stored modal fluxes (each cell owns its
  // lower faces, and stores the flux each of the two sides of the face receives).
  struct gkyl_dg_pkpm_remap_flux_surf_inp face_flux_inp = {
    .phase_grid = &s->kinetic->grid,
    .conf_basis = &app->confBasis,
    .phase_basis = &app->basis,
    .pos_map = app->pos_map,
    .num_face_nodes = s->kinetic->num_face_nodes,
    .use_gpu = app->use_gpu,
  };
  s->kinetic->calc_face_flux = gkyl_dg_pkpm_remap_flux_surf_inew(&face_flux_inp);
  int nkc = gkyl_dg_pkpm_remap_flux_surf_num_surf_basis(s->kinetic->calc_face_flux);
  s->kinetic->conf_flux_surf =
    mkarr(app->use_gpu, 2 * PKPM_DIST_NUM_BLOCKS * cdim * nkc, s->kinetic->local_ext.volume);
  s->kinetic->flux_blend = mkarr(app->use_gpu, cdim, app->local_ext.volume);
  gkyl_array_clear(s->kinetic->flux_blend, 0.0);

  struct gkyl_dg_pkpm_vel_flux_surf_inp vel_flux_inp = {
    .phase_grid = &s->kinetic->grid,
    .conf_basis = &app->confBasis,
    .phase_basis = &app->basis,
    .use_gpu = app->use_gpu,
  };
  s->kinetic->calc_vel_flux = gkyl_dg_pkpm_vel_flux_surf_inew(&vel_flux_inp);
  int nkv = gkyl_dg_pkpm_vel_flux_surf_num_surf_basis(s->kinetic->calc_vel_flux);
  s->kinetic->vel_flux_surf =
    mkarr(app->use_gpu, PKPM_DIST_NUM_BLOCKS * nkv, s->kinetic->local_ext.volume);

  int nkf = gkyl_dg_pkpm_remap_flux_surf_num_surf_basis_fluid(s->kinetic->calc_face_flux);
  s->fluid->fluid_flux_surf = mkarr(app->use_gpu, 5 * cdim * nkf, app->local_ext.volume);
  s->fluid->fluid_flux_surf_host =
    app->use_gpu ? mkarr(false, 5 * cdim * nkf, app->local_ext.volume) : 0;

  // Equation objects and hyper_dg solvers of the three passes.
  struct gkyl_dg_pkpm_kinetic_conf_inp kinetic_conf_inp = {
    .conf_basis = &app->confBasis,
    .phase_basis = &app->basis,
    .conf_range = &app->local,
    .phase_range = &s->kinetic->local,
    .pos_map = app->pos_map,
    .bvar = app->field->bvar,
    .u = s->kinetic->pkpm_u,
    .xi = s->kinetic->pkpm_xi,
    .conf_flux_surf = s->kinetic->conf_flux_surf,
    .use_gpu = app->use_gpu,
  };
  s->kinetic->eqn_kinetic_conf = gkyl_dg_pkpm_kinetic_conf_inew(&kinetic_conf_inp);
  struct gkyl_dg_pkpm_harm_vol_inp harm_vol_inp = {
    .phase_grid = &s->kinetic->grid,
    .conf_basis = &app->confBasis,
    .phase_basis = &app->basis,
    .qbym = s->charge / s->mass,
    .pos_map = app->pos_map,
    .use_gpu = app->use_gpu,
  };
  s->kinetic->harm_vol = gkyl_dg_pkpm_harm_vol_inew(&harm_vol_inp);
  s->kinetic->harm_coef = mkarr(
    app->use_gpu, gkyl_dg_pkpm_harm_vol_num_coef(s->kinetic->harm_vol), app->local_ext.volume
  );
  s->kinetic->harm_wnode = mkarr(
    app->use_gpu, gkyl_dg_pkpm_harm_vol_num_wnode(s->kinetic->harm_vol), app->local_ext.volume
  );
  // gyrophase mixing: default rate 10 |Omega|, switched off by a negative input
  double gyro_mix = s->kinetic->info.gyro_mix == 0.0 ? 10.0 : s->kinetic->info.gyro_mix;
  s->kinetic->harm_exp = gkyl_dg_pkpm_harm_exp_inew(&(struct gkyl_dg_pkpm_harm_exp_inp){
    .phase_grid = &s->kinetic->grid,
    .conf_basis = &app->confBasis,
    .qbym = s->charge / s->mass,
    .mix = gyro_mix,
    .use_gpu = app->use_gpu,
  });
  s->kinetic->harm_exp_coef = mkarr(
    app->use_gpu, gkyl_dg_pkpm_harm_exp_num_coef(s->kinetic->harm_exp), app->local_ext.volume
  );
  s->kinetic->harm_exp_nu = mkarr(
    app->use_gpu, gkyl_dg_pkpm_harm_exp_num_rate(s->kinetic->harm_exp), s->kinetic->local_ext.volume
  );
  s->kinetic->harm_exp_act = mkarr(app->use_gpu, 1, app->local_ext.volume);
  s->kinetic->harm_exp_col =
    mkarr(app->use_gpu, gkyl_dg_pkpm_harm_exp_num_col(s->kinetic->harm_exp), app->local_ext.volume);
  struct gkyl_dg_pkpm_kinetic_vel_inp kinetic_vel_inp = {
    .conf_basis = &app->confBasis,
    .phase_basis = &app->basis,
    .conf_range = &app->local,
    .phase_range = &s->kinetic->local,
    .c_0 = s->kinetic->pkpm_c_0,
    .alpha_0 = s->kinetic->pkpm_alpha_0,
    .c_h = s->kinetic->pkpm_c_h,
    .a_perp = s->kinetic->pkpm_a_perp,
    .vel_flux_surf = s->kinetic->vel_flux_surf,
    .use_gpu = app->use_gpu,
  };
  s->kinetic->eqn_kinetic_vel = gkyl_dg_pkpm_kinetic_vel_inew(&kinetic_vel_inp);
  struct gkyl_dg_pkpm_fluid_inp fluid_inp = {
    .conf_basis = &app->confBasis,
    .conf_range = &app->local,
    .pos_map = app->pos_map,
    .vlasov_pkpm_moms = s->fluid->pkpm_moms,
    .bvar = app->field->bvar,
    .u = s->fluid->pkpm_u,
    .fluid_flux_surf = s->fluid->fluid_flux_surf,
    .use_gpu = app->use_gpu,
  };
  s->fluid->eqn_fluid = gkyl_dg_pkpm_fluid_inew(&fluid_inp);

  int up_dirs_conf[GKYL_MAX_DIM] = {0, 1, 2}, zero_flux_conf[2 * GKYL_MAX_DIM] = {0};
  s->kinetic->slvr_kinetic_conf = gkyl_hyper_dg_new(
    &s->kinetic->grid, &app->basis, s->kinetic->eqn_kinetic_conf, cdim, up_dirs_conf,
    zero_flux_conf, 1, app->use_gpu
  );
  int up_dirs_vel[GKYL_MAX_DIM] = {cdim}, zero_flux_vel[2 * GKYL_MAX_DIM] = {0};
  // The lower velocity edge has no ghost cell (boundary term with the lowest cell's own
  // flux slot); the upper edge is the regular face between the uppermost cell and the
  // velocity ghost cell, whose slot the closure fills.
  zero_flux_vel[cdim] = 1;
  s->kinetic->slvr_kinetic_vel = gkyl_hyper_dg_new(
    &s->kinetic->grid, &app->basis, s->kinetic->eqn_kinetic_vel, 1, up_dirs_vel, zero_flux_vel, 1,
    app->use_gpu
  );
  s->fluid->slvr_fluid = gkyl_hyper_dg_new(
    &app->grid, &app->confBasis, s->fluid->eqn_fluid, cdim, up_dirs_conf, zero_flux_conf, 1,
    app->use_gpu
  );
  // set by the source and the temperature source where the species has them
  s->fluid->src_scale = 0.0;
  s->fluid->src_fluid = 0;
  s->fluid->has_src_bgk_rhs = false;
  s->fluid->src_bgk_rhs = 0;

  s->kinetic->kinetic_tm = 0.0;
  s->fluid->fluid_tm = 0.0;

  // array for storing F^2 in each cell
  s->kinetic->L2_f = mkarr(app->use_gpu, 1, s->kinetic->local_ext.volume);
  if (app->use_gpu) {
    s->kinetic->red_L2_f = gkyl_cu_malloc(sizeof(double));
    s->fluid->red_integ_diag = gkyl_cu_malloc(sizeof(double[9]));
  }
  // allocate dynamic-vector to store all-reduced integrated moments and F^2
  s->kinetic->integ_L2_f = gkyl_dynvec_new(GKYL_DOUBLE, 1);
  s->fluid->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, 9);
  s->fluid->pos_stat = gkyl_dynvec_new(GKYL_DOUBLE, PKPM_POS_NUM_STAT);
  s->kinetic->bflux_diag = gkyl_dynvec_new(GKYL_DOUBLE, 2 * cdim);
  s->kinetic->is_first_bflux_write_call = true;
  for (int i = 0; i < 2 * GKYL_MAX_CDIM; ++i) {
    s->kinetic->bflux[i] = 0.0;
  }
  s->fluid->is_first_pos_stat_write_call = true;
  for (int c = 0; c < PKPM_POS_NUM_STAT; ++c) {
    s->fluid->pos_count_since[c] = 0;
  }
  s->kinetic->is_first_integ_L2_write_call = true;
  s->fluid->is_first_integ_write_call = true;

  // Initialize applied acceleration for use in force update.
  s->fluid->app_accel = mkarr(app->use_gpu, 3 * app->confBasis.num_basis, app->local_ext.volume);
  gkyl_array_clear(s->fluid->app_accel, 0.0);
  s->fluid->has_app_accel = false;
  s->fluid->app_accel_evolve = false;
  // setup applied acceleration
  if (s->fluid->info.app_accel) {
    s->fluid->has_app_accel = true;
    if (s->fluid->info.app_accel_evolve) {
      s->fluid->app_accel_evolve = s->fluid->info.app_accel_evolve;
    }

    s->fluid->app_accel_host = s->fluid->app_accel;
    if (app->use_gpu) {
      s->fluid->app_accel_host = mkarr(false, 3 * app->confBasis.num_basis, app->local_ext.volume);
    }
    s->fluid->app_accel_proj = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
      .grid = &app->grid,
      .basis = &app->confBasis,
      .num_quad = app->confBasis.poly_order + 1,
      .num_ret_vals = 3,
      .eval = s->fluid->info.app_accel,
      .ctx = s->fluid->info.app_accel_ctx,
      .c2p_func = pkpm_pos_c2p,
      .c2p_func_ctx = app->pos_map,
    });
  }

  // determine collision type to use in PKPM update
  s->kinetic->collision_id = s->kinetic->info.collisions.collision_id;
  // not available: collision diagnostics and BGK relaxation to a fixed temperature (source_bgk)
  assert(!s->kinetic->info.collisions.write_coll_diagnostics);
  assert(!s->kinetic->info.collisions.fixed_temp_relax);
  if (s->kinetic->collision_id == GKYL_BGK_COLLISIONS) {
    pkpm_species_bgk_init(app, s, &s->kinetic->bgk);
  }
  if (s->kinetic->collision_id == GKYL_LBO_COLLISIONS) {
    pkpm_species_lbo_init(app, s, &s->kinetic->lbo);
  }
  pkpm_species_source_bgk_init(app, s, &s->kinetic->src_bgk);

  // An implicit kinetic relaxation that carries a temperature change (the shape factor of
  // implicit BGK with partners, or of the implicit temperature source) needs the weighted stage
  // combination.
  bool has_implicit_cross = s->kinetic->collision_id == GKYL_BGK_COLLISIONS &&
                            s->kinetic->info.collisions.is_implicit &&
                            s->kinetic->info.collisions.num_cross_collisions > 0;
  s->kinetic->has_stage_weights = !app->use_gpu &&
                                  (has_implicit_cross || s->kinetic->src_bgk.is_implicit);
  if (s->kinetic->has_stage_weights) {
    s->kinetic->stage_lam = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    gkyl_array_clear(s->kinetic->stage_lam, 1.0);
    s->kinetic->up_stage_weights =
      gkyl_pkpm_stage_weights_inew(&(struct gkyl_pkpm_stage_weights_inp){
        .conf_basis = &app->confBasis,
        .phase_basis = &app->basis,
        .use_gpu = false,
      });
  }

  // The LBO acts at the configuration-space nodes on a nodal copy of the distributions. The BGK
  // operator and the temperature source are weak products and need none.
  s->kinetic->has_nodal_ops = !app->use_gpu && s->kinetic->collision_id == GKYL_LBO_COLLISIONS;
  if (s->kinetic->has_nodal_ops) {
    s->kinetic->f_nodal = mkarr(false, s->kinetic->f->ncomp, s->kinetic->local_ext.volume);
    s->kinetic->rhs_nodal = mkarr(false, s->kinetic->f->ncomp, s->kinetic->local_ext.volume);
  }

  // momentum diffusion, if present, on the position map (orders 2 and 4)
  s->fluid->has_diffusion = false;
  if (s->fluid->info.diffusion.D) {
    s->fluid->has_diffusion = true;
    int order = s->fluid->info.diffusion.order < 2 ? 2 : s->fluid->info.diffusion.order;
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
      diff_inp.D[d] = s->fluid->info.diffusion.D;
    }
    s->fluid->mom_diff = gkyl_dg_pkpm_mom_diffusion_inew(&diff_inp);
  }

  // determine which directions are not periodic and get non-periodic BC info
  // also create local skin-ghost ranges for ease of handling non-periodic BCs
  int num_periodic_dir = app->num_periodic_dir, is_np[3] = {1, 1, 1};
  for (int d = 0; d < num_periodic_dir; ++d) {
    is_np[app->periodic_dirs[d]] = 0;
  }

  for (int dir = 0; dir < cdim; ++dir) {
    s->kinetic->lower_bc[dir] = s->kinetic->upper_bc[dir] = GKYL_SPECIES_COPY;
    s->fluid->lower_bc[dir] = s->fluid->upper_bc[dir] = GKYL_SPECIES_COPY;
    if (is_np[dir]) {
      const struct gkyl_vlasov_bcs *bc;
      const enum gkyl_species_bc_type *bc_fluid;
      if (dir == 0) {
        bc = &s->kinetic->info.bcx;
        bc_fluid = s->fluid->info.bcx;
      } else if (dir == 1) {
        bc = &s->kinetic->info.bcy;
        bc_fluid = s->fluid->info.bcy;
      } else {
        bc = &s->kinetic->info.bcz;
        bc_fluid = s->fluid->info.bcz;
      }

      s->kinetic->lower_bc[dir] = bc->lower.type;
      s->kinetic->upper_bc[dir] = bc->upper.type;
      s->fluid->lower_bc[dir] = bc_fluid[0];
      s->fluid->upper_bc[dir] = bc_fluid[1];
    }
    // Create local lower skin and ghost ranges for distribution function
    gkyl_skin_ghost_ranges(
      &s->kinetic->lower_skin[dir], &s->kinetic->lower_ghost[dir], dir, GKYL_LOWER_EDGE,
      &s->kinetic->local_ext_conf, ghost
    );
    // Create local upper skin and ghost ranges for distribution function
    gkyl_skin_ghost_ranges(
      &s->kinetic->upper_skin[dir], &s->kinetic->upper_ghost[dir], dir, GKYL_UPPER_EDGE,
      &s->kinetic->local_ext_conf, ghost
    );
  }

  // allocate buffer for applying BCs
  long buff_sz_dist = 0;
  long buff_sz_fluid = 0;
  // compute buffer size needed
  for (int dir = 0; dir < cdim; ++dir) {
    long vol_dist =
      GKYL_MAX2(s->kinetic->lower_skin[dir].volume, s->kinetic->upper_skin[dir].volume);
    buff_sz_dist = buff_sz_dist > vol_dist ? buff_sz_dist : vol_dist;
    long vol_fluid = GKYL_MAX2(app->lower_skin[dir].volume, app->upper_skin[dir].volume);
    buff_sz_fluid = buff_sz_fluid > vol_fluid ? buff_sz_fluid : vol_fluid;
  }

  s->kinetic->bc_buffer = mkarr(app->use_gpu, s->kinetic->f->ncomp, buff_sz_dist);
  s->fluid->bc_buffer = mkarr(app->use_gpu, 5 * app->confBasis.num_basis, buff_sz_fluid);
  // buffer arrays for fixed function boundary conditions on the distribution functions and the
  // fluid state, for the edges that have one
  for (int d = 0; d < 3; ++d) {
    s->kinetic->bc_buffer_lo_fixed[d] = s->kinetic->bc_buffer_up_fixed[d] = 0;
    s->fluid->bc_buffer_lo_fixed[d] = s->fluid->bc_buffer_up_fixed[d] = 0;
  }
  for (int d = 0; d < cdim; ++d) {
    int nc_fluid = 5 * app->confBasis.num_basis;
    if (s->kinetic->lower_bc[d] == GKYL_SPECIES_FIXED_FUNC) {
      s->kinetic->bc_buffer_lo_fixed[d] =
        mkarr(app->use_gpu, s->kinetic->f->ncomp, s->kinetic->lower_skin[d].volume);
    }
    if (s->kinetic->upper_bc[d] == GKYL_SPECIES_FIXED_FUNC) {
      s->kinetic->bc_buffer_up_fixed[d] =
        mkarr(app->use_gpu, s->kinetic->f->ncomp, s->kinetic->upper_skin[d].volume);
    }
    if (s->fluid->lower_bc[d] == GKYL_SPECIES_FIXED_FUNC) {
      s->fluid->bc_buffer_lo_fixed[d] = mkarr(app->use_gpu, nc_fluid, app->lower_skin[d].volume);
    }
    if (s->fluid->upper_bc[d] == GKYL_SPECIES_FIXED_FUNC) {
      s->fluid->bc_buffer_up_fixed[d] = mkarr(app->use_gpu, nc_fluid, app->upper_skin[d].volume);
    }
  }

  // Absorbing boundaries: the distribution functions vanish in the ghost cells; the
  // fluid state is the mirror image of the skin cell (even in the wall direction), so the
  // ghost trace at the wall is the skin trace and the primitive variables, their weak
  // gradients and the face averages of the flux coefficients are one-sided at the wall. The
  // fluid fluxes through the wall are the moments of the kinetic flux against the zero ghost
  // distributions (pkpm_species_face_flux), not the fluid flux of the ghost state.
  s->kinetic->bc_is_absorb = false;
  for (int d = 0; d < cdim; ++d) {
    s->kinetic->bc_lo_reflect[d] = 0;
    s->kinetic->bc_up_reflect[d] = 0;
    s->kinetic->bc_lo_theta[d] = s->kinetic->bc_up_theta[d] = 0;
    // Lower BC updaters. Each block has its own type; copy BCs by default.
    enum gkyl_bc_basic_type bctype_dist = pkpm_species_bc_basic_type_dist(s->kinetic->lower_bc[d]);
    enum gkyl_bc_basic_type bctype_fluid = pkpm_species_bc_basic_type_fluid(s->fluid->lower_bc[d]);
    if (s->kinetic->lower_bc[d] == GKYL_SPECIES_ABSORB) {
      s->kinetic->bc_is_absorb = true;
      s->kinetic->bc_lo_theta[d] = gkyl_bc_basic_new(
        d, GKYL_LOWER_EDGE, GKYL_BC_REFLECT, app->basis_on_dev.confBasis, &app->lower_skin[d],
        &app->lower_ghost[d], app->confBasis.num_basis, app->cdim, app->use_gpu
      );
    } else if (s->kinetic->lower_bc[d] == GKYL_SPECIES_REFLECT) {
      s->kinetic->bc_lo_reflect[d] = pkpm_species_reflect_new(app, s, d, GKYL_LOWER_EDGE);
    }

    // Distribution function non-periodic lower boundary conditions
    s->kinetic->bc_lo[d] = gkyl_bc_basic_new(
      d, GKYL_LOWER_EDGE, bctype_dist, app->basis_on_dev.basis, &s->kinetic->lower_skin[d],
      &s->kinetic->lower_ghost[d], s->kinetic->f->ncomp, app->cdim, app->use_gpu
    );
    // Fluid state non-periodic lower boundary conditions
    s->fluid->bc_lo[d] = gkyl_bc_basic_new(
      d, GKYL_LOWER_EDGE, bctype_fluid, app->basis_on_dev.confBasis, &app->lower_skin[d],
      &app->lower_ghost[d], s->fluid->fluid->ncomp, app->cdim, app->use_gpu
    );

    // Upper BC updaters.
    bctype_dist = pkpm_species_bc_basic_type_dist(s->kinetic->upper_bc[d]);
    bctype_fluid = pkpm_species_bc_basic_type_fluid(s->fluid->upper_bc[d]);
    if (s->kinetic->upper_bc[d] == GKYL_SPECIES_ABSORB) {
      s->kinetic->bc_is_absorb = true;
      s->kinetic->bc_up_theta[d] = gkyl_bc_basic_new(
        d, GKYL_UPPER_EDGE, GKYL_BC_REFLECT, app->basis_on_dev.confBasis, &app->upper_skin[d],
        &app->upper_ghost[d], app->confBasis.num_basis, app->cdim, app->use_gpu
      );
    } else if (s->kinetic->upper_bc[d] == GKYL_SPECIES_REFLECT) {
      s->kinetic->bc_up_reflect[d] = pkpm_species_reflect_new(app, s, d, GKYL_UPPER_EDGE);
    }

    // Distribution function non-periodic upper boundary conditions
    s->kinetic->bc_up[d] = gkyl_bc_basic_new(
      d, GKYL_UPPER_EDGE, bctype_dist, app->basis_on_dev.basis, &s->kinetic->upper_skin[d],
      &s->kinetic->upper_ghost[d], s->kinetic->f->ncomp, app->cdim, app->use_gpu
    );
    // Fluid state non-periodic upper boundary conditions
    s->fluid->bc_up[d] = gkyl_bc_basic_new(
      d, GKYL_UPPER_EDGE, bctype_fluid, app->basis_on_dev.confBasis, &app->upper_skin[d],
      &app->upper_ghost[d], s->fluid->fluid->ncomp, app->cdim, app->use_gpu
    );
  }
}

void
pkpm_species_apply_ic(gkyl_pkpm_app *app, struct pkpm_species *species, double t0)
{
  int poly_order = app->poly_order;
  // The user provides the normalized distributions [F, G, M0_perp, N_perp] on the normalized
  // velocity grid (w = v_par/xi; F a probability density in w, G its perpendicular-energy
  // partner so that M0(F) = 1, M1(F) = 0 and M2(F) + 2 M0(G) = 3; the first-harmonic blocks
  // are zero for a gyrotropic start) and the fluid state [rho, rho ux, rho uy, rho uz, E].
  // Both are sampled at the physical coordinates of the position map (identity: unchanged).
  gkyl_proj_on_basis *proj_dist = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
    .grid = &species->kinetic->grid,
    .basis = &app->basis,
    .num_quad = 8,
    .num_ret_vals = PKPM_DIST_NUM_BLOCKS,
    .eval = species->kinetic->info.projection[0].func,
    .ctx = species->kinetic->info.projection[0].ctx_func,
    .c2p_func = pkpm_pos_c2p,
    .c2p_func_ctx = app->pos_map,
  });
  gkyl_proj_on_basis *proj_fluid = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
    .grid = &app->grid,
    .basis = &app->confBasis,
    .num_quad = poly_order + 1,
    .num_ret_vals = 5,
    .eval = species->fluid->info.init,
    .ctx = species->fluid->info.ctx,
    .c2p_func = pkpm_pos_c2p,
    .c2p_func_ctx = app->pos_map,
  });

  // run updaters; need to project onto extended range for ease of handling
  // subsequent operations over extended range such as primitive variable computations
  // This is needed to fill the corner cells as the corner cells may not be filled by
  // boundary conditions and we cannot divide by 0 anywhere or the weak divisions will fail
  // (the velocity ghost cell stays empty)
  gkyl_array_clear(species->kinetic->f_host, 0.0);
  gkyl_proj_on_basis_advance(
    proj_dist, t0, &species->kinetic->local_ext_conf, species->kinetic->f_host
  );
  gkyl_proj_on_basis_release(proj_dist);
  gkyl_proj_on_basis_advance(proj_fluid, t0, &app->local_ext, species->fluid->fluid_host);
  gkyl_proj_on_basis_release(proj_fluid);

  // note: f_host and fluid_host are the same as f and fluid respectively when not on GPUs
  if (app->use_gpu) {
    gkyl_array_copy(species->kinetic->f, species->kinetic->f_host);
    gkyl_array_copy(species->fluid->fluid, species->fluid->fluid_host);
  }

  // Constraints M0(F) = 1 and M2(F) + 2 M0(G) = 3 exactly: nodal rescaling of the
  // projected distributions by their own moments (removes the projection error).
  pkpm_species_moment_calc(
    &species->kinetic->pkpm_nmoms, species->kinetic->local_ext_conf, app->local_ext,
    species->kinetic->f
  );
  gkyl_array_copy(species->kinetic->f1, species->kinetic->f);
  gkyl_dg_calc_pkpm_vars_normalize_dist(
    species->kinetic->calc_pkpm_vars, &app->local_ext, &species->kinetic->local_ext,
    species->kinetic->pkpm_nmoms.marr, species->kinetic->f1, species->kinetic->f
  );
  // The evolved quantities carry the position-map Jacobian: J [F, G] and J [rho, M, E]
  // (identity map: no change). The velocity ghost row of f stays empty.
  gkyl_vlasov_position_map_rescale_jacobpos(
    app->pos_map, &app->basis, &species->kinetic->local_ext_conf, species->kinetic->f,
    species->kinetic->f
  );
  gkyl_vlasov_position_map_rescale_jacobpos_conf(
    app->pos_map, &app->local_ext, species->fluid->fluid, species->fluid->fluid
  );
  if (app->use_gpu) {
    gkyl_array_copy(species->kinetic->f_host, species->kinetic->f);
    gkyl_array_copy(species->fluid->fluid_host, species->fluid->fluid);
  }

  // pre-compute the applied acceleration (it is not recomputed when time-independent)
  pkpm_species_calc_app_accel(app, species, t0);

  // copy the initial conditions into the buffers of the fixed-function boundary conditions
  for (int d = 0; d < app->cdim; ++d) {
    if (species->kinetic->lower_bc[d] == GKYL_SPECIES_FIXED_FUNC) {
      gkyl_bc_basic_buffer_fixed_func(
        species->kinetic->bc_lo[d], species->kinetic->bc_buffer_lo_fixed[d], species->kinetic->f
      );
    }
    if (species->kinetic->upper_bc[d] == GKYL_SPECIES_FIXED_FUNC) {
      gkyl_bc_basic_buffer_fixed_func(
        species->kinetic->bc_up[d], species->kinetic->bc_buffer_up_fixed[d], species->kinetic->f
      );
    }
    if (species->fluid->lower_bc[d] == GKYL_SPECIES_FIXED_FUNC) {
      gkyl_bc_basic_buffer_fixed_func(
        species->fluid->bc_lo[d], species->fluid->bc_buffer_lo_fixed[d], species->fluid->fluid
      );
    }
    if (species->fluid->upper_bc[d] == GKYL_SPECIES_FIXED_FUNC) {
      gkyl_bc_basic_buffer_fixed_func(
        species->fluid->bc_up[d], species->fluid->bc_buffer_up_fixed[d], species->fluid->fluid
      );
    }
  }
}

void
pkpm_species_calc_app_accel(gkyl_pkpm_app *app, struct pkpm_species *species, double tm)
{
  if (species->fluid->has_app_accel) {
    gkyl_proj_on_basis_advance(
      species->fluid->app_accel_proj, tm, &app->local_ext, species->fluid->app_accel_host
    );
    if (app->use_gpu) {
      // note: app_accel_host is same as app_accel when not on GPUs
      gkyl_array_copy(species->fluid->app_accel, species->fluid->app_accel_host);
    }
  }
}

// Kinetic coupling moments [p_par, p_perp, q_par, q_perp, s, r] = [p M2, p M0_G, p xi M3, p xi M1_G,
// p m1, p xi (m2/2 + n0)] from the raw normalized moments, the nodal scalar pressure p = rho xi^2
// and the thermal speed xi (weak products projected onto p=1): the pressures, first-harmonic stress
// and heat fluxes the fluid system needs. The constraints
// [M0, M1, (M2 + 2 M0_G)/3] are copied out for the diagnostics.
static void
pkpm_species_coupling_moms(gkyl_pkpm_app *app, struct pkpm_species *species)
{
  int nc = app->confBasis.num_basis;
  // layout of the normalized moments: [M0, M2, M0_G, M1, M3, M1_G, mu (3), m1 (3), m2 (3), n0 (3)]
  static const int raw[7] = {1, 2, 4, 5, 9, 10, 11}; // M2, M0_G, M3, M1_G, m1 (3)
  for (int c = 0; c < 10; ++c) {
    if (c < 7) {
      gkyl_dg_mul_op_range(
        &app->confBasis, c, species->fluid->pkpm_moms, 2, species->fluid->pkpm_xi, raw[c],
        species->kinetic->pkpm_nmoms.marr, &app->local_ext
      );
    } else {
      int a = c - 7;
      gkyl_array_set_offset(
        species->fluid->pkpm_scratch, 0.5, species->kinetic->pkpm_nmoms.marr, (12 + a) * nc
      );
      gkyl_array_accumulate_offset(
        species->fluid->pkpm_scratch, 1.0, species->kinetic->pkpm_nmoms.marr, (15 + a) * nc
      );
      gkyl_dg_mul_op_range(
        &app->confBasis, c, species->fluid->pkpm_moms, 2, species->fluid->pkpm_xi, 0,
        species->fluid->pkpm_scratch, &app->local_ext
      );
    }
    if (c == 2 || c == 3 || c >= 7) {
      gkyl_dg_mul_op_range(
        &app->confBasis, 0, species->fluid->pkpm_scratch, c, species->fluid->pkpm_moms, 0,
        species->fluid->pkpm_xi, &app->local_ext
      );
      gkyl_array_set_offset(species->fluid->pkpm_moms, 1.0, species->fluid->pkpm_scratch, c * nc);
    }
  }
  gkyl_dg_calc_pkpm_vars_constraints(
    species->kinetic->calc_pkpm_vars, &app->local_ext, nc, species->kinetic->pkpm_nmoms.marr,
    species->kinetic->pkpm_cmoms
  );
  // the constraints of the physical distributions ([1, 0, 1]): the moments of J [F, G] carry J
  gkyl_vlasov_position_map_divide_jacobpos_conf(
    app->pos_map, &app->local_ext, species->kinetic->pkpm_cmoms->ncomp,
    species->kinetic->pkpm_cmoms, species->kinetic->pkpm_cmoms
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
    app->pos_map, &app->local_ext, fluidin->ncomp, fluidin, species->fluid->fluid_no_J
  );
  pkpm_species_moment_calc(
    &species->kinetic->pkpm_nmoms, species->kinetic->local_ext_conf, app->local_ext, fin
  );
  gkyl_array_set_offset(species->fluid->pkpm_rho, 1.0, species->fluid->fluid_no_J, 0);
  // The flow u = M/rho and xi by nodal quotients: for p=1 the weak division, rho u = M weakly,
  // which the implicit field coupling needs to conserve momentum (it forms the new momentum
  // as rho (2 u_half - u)).
  gkyl_dg_calc_pkpm_vars_prim(
    species->fluid->calc_pkpm_vars, &app->local_ext, species->fluid->fluid_no_J,
    species->fluid->pkpm_u, species->fluid->pkpm_xi
  );
  // the kinetic aspect reads its own copies of the flow and the thermal speed
  gkyl_array_copy(species->kinetic->pkpm_u, species->fluid->pkpm_u);
  gkyl_array_copy(species->kinetic->pkpm_xi, species->fluid->pkpm_xi);
  pkpm_species_coupling_moms(app, species);
  // theta = M0(G) of the physical G: M0(J G)/J.
  gkyl_array_set_offset(
    species->kinetic->pkpm_theta, 1.0, species->kinetic->pkpm_nmoms.marr, 2 * nc
  );
  gkyl_vlasov_position_map_divide_jacobpos_conf(
    app->pos_map, &app->local_ext, nc, species->kinetic->pkpm_theta, species->kinetic->pkpm_theta
  );
  // The face fluxes average u, xi and theta across every face. At absorbing walls u and xi of
  // the mirrored fluid ghost are the mirror image of the skin values already; theta, the
  // zeroth moment of the empty ghost G, is replaced by the mirror image of the skin cell, so
  // that the face averages are the interior wall traces.
  if (species->kinetic->bc_is_absorb) {
    for (int d = 0; d < app->cdim; ++d) {
      if (species->kinetic->bc_lo_theta[d]) {
        gkyl_bc_basic_advance(
          species->kinetic->bc_lo_theta[d], species->fluid->bc_buffer, species->kinetic->pkpm_theta
        );
      }
      if (species->kinetic->bc_up_theta[d]) {
        gkyl_bc_basic_advance(
          species->kinetic->bc_up_theta[d], species->fluid->bc_buffer, species->kinetic->pkpm_theta
        );
      }
    }
    gkyl_comm_array_sync(app->comm, &app->local, &app->local_ext, species->kinetic->pkpm_theta);
  }

  app->stat.species_pkpm_vars_tm += gkyl_time_diff_now_sec(tm);
}

void
pkpm_species_apply_positivity(
  gkyl_pkpm_app *app, struct pkpm_fluid_species *fl, struct gkyl_array *fluid
)
{
  struct timespec tm = gkyl_wall_clock();
  int num_rho = 0, num_e = 0, num_avg = 0;
  gkyl_dg_calc_pkpm_vars_positivity(
    fl->calc_pkpm_vars, &app->local, fl->positivity_eps, fluid, fl->pkpm_pos_count, &num_rho,
    &num_e, &num_avg
  );
  fl->pos_count_since[0] += num_rho;
  fl->pos_count_since[1] += num_e;
  fl->pos_count_since[2] += num_avg;
  fl->pos_count_since[3] += 1;
  app->stat.species_pkpm_vars_tm += gkyl_time_diff_now_sec(tm);
}

// Configuration-space face fluxes of the stage: the kinetic flux each side of a face receives,
// the fluid flux, as moments of the same face flux, and both rates. The walls need no special
// case: an absorbing ghost (F = 0, fluid and theta mirrored) gives the interior side's
// contribution with the interior wall coefficients, and a reflecting ghost cancels the mass and
// energy fluxes.
//
// The face maps the jumps of the frame exactly, so its gradients in the volume terms are in-cell
// except for the share of each jump the label-paired flux leaves unmapped: they follow the flux.
void
pkpm_species_face_flux(
  gkyl_pkpm_app *app, struct pkpm_species *species, const struct gkyl_array *fin
)
{
  gkyl_dg_pkpm_remap_flux_surf_advance(
    species->kinetic->calc_face_flux, &app->local, &app->local_ext, &species->kinetic->local,
    &species->kinetic->local_ext, app->field->bvar, species->fluid->fluid_no_J,
    species->fluid->pkpm_u, species->fluid->pkpm_xi, species->kinetic->pkpm_theta, fin,
    species->kinetic->conf_flux_surf, species->fluid->fluid_flux_surf, species->kinetic->cflrate,
    species->fluid->cflrate, species->kinetic->flux_blend
  );
  // The gradient of u for the kinetic sources (bb:grad u); needs u in the ghost cells (at
  // absorbing boundaries that of the mirrored fluid state).
  gkyl_dg_calc_pkpm_vars_in_cell_grad(
    species->fluid->calc_pkpm_vars, &app->local, species->kinetic->flux_blend,
    species->fluid->pkpm_u, species->kinetic->pkpm_grad_u
  );
  // grad(ln n) for the density-gradient source: the gradient of [rho, 0, 0] divided by rho.
  gkyl_array_set_offset(species->fluid->pkpm_rho3, 1.0, species->fluid->pkpm_rho, 0);
  gkyl_dg_calc_pkpm_vars_in_cell_grad(
    species->fluid->calc_pkpm_vars, &app->local, species->kinetic->flux_blend,
    species->fluid->pkpm_rho3, species->fluid->pkpm_grad_rho
  );
  gkyl_dg_calc_pkpm_vars_grad_ln_n(
    species->fluid->calc_pkpm_vars, &app->local, species->fluid->pkpm_rho,
    species->fluid->pkpm_grad_rho, species->kinetic->pkpm_grad_ln_n
  );
  // The gradient of xi for the forces -w^2 b.grad xi and -2 w b.grad xi G (the second block of
  // pkpm_xi3 receives xi^2 and is never read).
  gkyl_array_set_offset(species->fluid->pkpm_xi3, 1.0, species->fluid->pkpm_xi, 0);
  gkyl_dg_calc_pkpm_vars_in_cell_grad(
    species->fluid->calc_pkpm_vars, &app->local, species->kinetic->flux_blend,
    species->fluid->pkpm_xi3, species->kinetic->pkpm_grad_xi
  );
}

// Collisionless configuration-space pass of the kinetic update into rhs (cleared here):
// the HLL streaming fluxes, the exact volume terms and the explicit sources. The closures
// are computed from the moments of this partial right-hand side, after the sources and
// collision operators have been accumulated into it.
static void
pkpm_species_kinetic_conf_pass(
  gkyl_pkpm_app *app, struct pkpm_kinetic_species *kin, const struct gkyl_array *fin,
  struct gkyl_array *rhs
)
{
  gkyl_array_clear(rhs, 0.0);
  gkyl_hyper_dg_advance(kin->slvr_kinetic_conf, &kin->local, fin, kin->cflrate, rhs);
  // the known first-harmonic w fluxes and sources of the eight blocks (modal, exact)
  gkyl_dg_pkpm_harm_vol_coef(
    kin->harm_vol, &app->local, app->field->bvar, app->field->grad_b, kin->pkpm_u, kin->pkpm_grad_u,
    kin->pkpm_xi, kin->pkpm_grad_xi, kin->pkpm_grad_ln_n, kin->pkpm_theta,
    app->field->dtB_over_magB, app->field->tot_em, kin->harm_coef, kin->harm_wnode
  );
  gkyl_dg_pkpm_harm_vol_advance(
    kin->harm_vol, &app->local, &kin->local, kin->harm_coef, fin, kin->cflrate, rhs
  );
}

void
pkpm_species_accumulate_conf_phase(
  gkyl_pkpm_app *app, struct pkpm_kinetic_species *kin, struct gkyl_array *out, double a,
  const struct gkyl_array *conf, const struct gkyl_array *phase
)
{
  int np = app->basis.num_basis;
  for (int b = 0; b < phase->ncomp / np; ++b) {
    gkyl_array_set_offset(kin->blk_in, 1.0, phase, b * np);
    gkyl_array_clear(kin->blk_out, 0.0);
    gkyl_dg_mul_conf_phase_op_accumulate_range(
      &app->confBasis, &app->basis, kin->blk_out, a, conf, kin->blk_in, &app->local, &kin->local
    );
    gkyl_array_accumulate_offset_range(out, 1.0, kin->blk_out, b * np, &kin->local);
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
    species->kinetic->calc_pkpm_vars, &app->local_ext, &species->kinetic->local_ext,
    app->field->bvar, f
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
      species->kinetic->harm_exp, &app->local, app->field->tot_em, species->kinetic->harm_exp_coef
    );
    gkyl_array_clear(species->kinetic->harm_exp_act, 0.0);
  }
  // With stage weights the combination of the driver is weighted (pkpm_coll_stage.c), and the
  // stage value is built for those weights.
  gkyl_dg_pkpm_harm_exp_advance(
    species->kinetic->harm_exp, stage, h, &app->local, &species->kinetic->local,
    species->kinetic->harm_exp_coef, species->kinetic->f,
    species->kinetic->has_stage_weights ? species->kinetic->stage_lam : 0,
    species->kinetic->harm_exp_nu, species->kinetic->harm_exp_act, fout
  );
  gkyl_dg_pkpm_harm_exp_restore(
    species->kinetic->harm_exp, &app->local, &species->kinetic->local,
    species->kinetic->harm_exp_act, species->kinetic->harm_exp_col, fout
  );
}

// Closures from the moments of the partial right-hand side rhs (configuration pass, sources and
// collisions): c_0 M0 = -M0(rhs), (alpha_0, c) from the first and scalar second moments and
// a_perp from the perpendicular rows I_perp int M0_perp dw. The velocity-edge fluxes enter the
// solve and are stored for the velocity pass, so that the complete right-hand side keeps the
// constraints. Their edge speeds use the closures of the stage before; from_step takes those of
// the last stage of the last completed step (the first stage of a step, the diagnostics).
// Returns the number of nodes where the closures are undefined.
static int
pkpm_species_kinetic_closure(
  gkyl_pkpm_app *app, struct pkpm_kinetic_species *kin, const struct gkyl_array *fin,
  const struct gkyl_array *rhs, bool from_step
)
{
  if (from_step) {
    gkyl_array_set_offset(kin->pkpm_alpha_0, 1.0, kin->pkpm_closure_step, 0);
    gkyl_array_set_offset(kin->pkpm_c_h, 1.0, kin->pkpm_closure_step, app->confBasis.num_basis);
  }
  gkyl_array_clear(kin->pkpm_rhs_moms, 0.0);
  pkpm_species_moment_calc_into(&kin->pkpm_nmoms, kin->local, app->local, rhs, kin->pkpm_rhs_moms);
  return gkyl_dg_calc_pkpm_vars_closure(
    kin->calc_pkpm_vars, &app->local, &kin->local, &kin->grid, kin->pkpm_nmoms.marr,
    kin->pkpm_rhs_moms, app->field->bvar, kin->harm_wnode, kin->pkpm_theta, fin, kin->pkpm_c_0,
    kin->pkpm_alpha_0, kin->pkpm_c_h, kin->pkpm_a_perp, kin->vel_flux_surf, kin->cflrate
  );
}

void
pkpm_species_complete_step(gkyl_pkpm_app *app, struct pkpm_species *sp)
{
  struct pkpm_kinetic_species *kin = sp->kinetic;
  gkyl_array_set_offset(kin->pkpm_closure_step, 1.0, kin->pkpm_alpha_0, 0);
  gkyl_array_set_offset(kin->pkpm_closure_step, 1.0, kin->pkpm_c_h, app->confBasis.num_basis);
}

// Record the nulled closure nodes of a stage in the statistics of the species.
static void
pkpm_species_kinetic_closure_stat(gkyl_pkpm_app *app, struct pkpm_kinetic_species *kin, int num_null)
{
  kin->num_null_closure += num_null;
  kin->closure_cond_min =
    fmin(kin->closure_cond_min, gkyl_dg_calc_pkpm_vars_closure_cond_min(kin->calc_pkpm_vars));
  // A nulled node means the traces at the velocity-domain edges are large enough to make the
  // closure system singular: the constraints drift in that cell. Say so once, with the proxy.
  if (num_null > 0 && !kin->closure_null_warned) {
    kin->closure_null_warned = true;
    int rank = 0;
    gkyl_comm_get_rank(app->comm, &rank);
    fprintf(
      stderr,
      "WARNING (rank %d): species %s at t = %g: the closure solve was nulled at %d nodes (min det/det0 so far %g).\n"
      "  The normalized distribution is not negligible at the velocity-space boundary; the constraints drift there.\n"
      "  The velocity extents are likely too small for this species.\n",
      rank, kin->name, app->tcurr, num_null, kin->closure_cond_min
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
  gkyl_array_clear(species->kinetic->cflrate, 0.0);
  gkyl_array_clear(species->fluid->cflrate, 0.0);
  pkpm_species_face_flux(app, species, fin);
  // fnew is free scratch between time steps.
  pkpm_species_kinetic_conf_pass(app, species->kinetic, fin, species->kinetic->fnew);
  pkpm_species_kinetic_closure(app, species->kinetic, fin, species->kinetic->fnew, true);
}

// Right-hand side of the kinetic equations of a species. The face fluxes of the stage, the source
// distributions and the rates were staged for every species before any right-hand side.
double
pkpm_kinetic_species_rhs(
  gkyl_pkpm_app *app, struct pkpm_kinetic_species *kin, int stage, const struct gkyl_array *fin,
  const struct gkyl_array *em, struct gkyl_array *rhs
)
{
  struct timespec tm = gkyl_wall_clock();
  // Kinetic equations: the collisionless configuration pass (clears rhs), the sources and
  // collisions, the closures from the moments of that partial right-hand side, then the
  // velocity pass accumulating into the same rhs.
  pkpm_species_kinetic_conf_pass(app, kin, fin, rhs);
  // Sources and collisions enter before the closures, which keep the constraints exact
  // against them (the fluid carries the density and energy changes). The LBO conserves the
  // constraints on its own (boundary corrections, zero-flux w edges).
  pkpm_species_source_kinetic_rhs(app, kin, &kin->src, fin, rhs);
  pkpm_species_source_bgk_kinetic_rhs(app, kin, &kin->src_bgk, fin, rhs);
  if (kin->collision_id == GKYL_BGK_COLLISIONS) {
    pkpm_species_bgk_kinetic_rhs(app, kin, &kin->bgk, fin, rhs);
  }
  if (kin->collision_id == GKYL_LBO_COLLISIONS) {
    pkpm_species_lbo_rhs(app, kin, &kin->lbo, fin, rhs);
  }
  // The operators at the nodes accumulated into the nodal right-hand side.
  pkpm_species_nodal_accumulate(app, kin, rhs);
  int num_null = pkpm_species_kinetic_closure(app, kin, fin, rhs, stage == 1);
  pkpm_species_kinetic_closure_stat(app, kin, num_null);
  gkyl_dg_pkpm_vel_flux_surf_advance(
    kin->calc_vel_flux, &app->local, &kin->local, kin->harm_wnode, kin->pkpm_theta,
    app->field->bvar, kin->pkpm_alpha_0, kin->pkpm_c_h, fin, kin->vel_flux_surf
  );
  gkyl_hyper_dg_advance(kin->slvr_kinetic_vel, &kin->local, fin, kin->cflrate, rhs);
  kin->kinetic_tm += gkyl_time_diff_now_sec(tm);

  app->stat.n_species_omega_cfl += 1;
  tm = gkyl_wall_clock();
  gkyl_array_reduce_range(kin->omega_cfl, kin->cflrate, GKYL_MAX, &kin->local);

  double omega_cfl_ho[1];
  if (app->use_gpu) {
    gkyl_cu_memcpy(omega_cfl_ho, kin->omega_cfl, sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    omega_cfl_ho[0] = kin->omega_cfl[0];
  }
  double omega_cfl = omega_cfl_ho[0];

  app->stat.species_omega_cfl_tm += gkyl_time_diff_now_sec(tm);

  return app->cfl / omega_cfl;
}

// Right-hand side of the fluid equations of a species.
double
pkpm_fluid_species_rhs(
  gkyl_pkpm_app *app, struct pkpm_fluid_species *fl, const struct gkyl_array *fluidin,
  const struct gkyl_array *em, struct gkyl_array *rhs
)
{
  gkyl_array_clear(rhs, 0.0);

  // The energy relaxation of the explicit temperature source of the stage (staged with the
  // species), then the volume and surface terms from the stored face fluxes and the source.
  struct timespec tm = gkyl_wall_clock();
  if (fl->has_src_bgk_rhs) {
    gkyl_array_accumulate_range(rhs, 1.0, fl->src_bgk_rhs, &app->local);
  }
  gkyl_hyper_dg_advance(fl->slvr_fluid, &app->local, fluidin, fl->cflrate, rhs);
  if (fl->src_scale != 0.0) {
    gkyl_array_accumulate_range(rhs, fl->src_scale, fl->src_fluid, &app->local);
  }
  fl->fluid_tm += gkyl_time_diff_now_sec(tm);

  if (fl->has_diffusion) {
    gkyl_dg_pkpm_mom_diffusion_advance(fl->mom_diff, &app->local, fluidin, fl->cflrate, rhs);
  }

  // If PKPM update is fully explicit, include the fluid-EM coupling in fluid update
  if (app->use_explicit_sources) {
    double qbym = fl->charge / fl->mass;
    gkyl_array_set(fl->qmem, qbym, em);

    // Accumulate applied acceleration and/or q/m*(external electromagnetic)
    // fields onto qmem to get the total acceleration
    if (fl->has_app_accel) {
      gkyl_array_accumulate_range(fl->qmem, 1.0, fl->app_accel, &app->local);
    }
    if (app->field->has_ext_em) {
      gkyl_array_accumulate_range(fl->qmem, qbym, app->field->ext_em, &app->local);
    }

    gkyl_dg_calc_pkpm_vars_source(
      fl->calc_pkpm_vars, &app->local, fl->qmem, app->src_omega_p_sq, fluidin, fl->cflrate, rhs
    );
  }

  tm = gkyl_wall_clock();
  gkyl_array_reduce_range(fl->omegaCfl_ptr, fl->cflrate, GKYL_MAX, &app->local);

  double omega_cfl_ho[1];
  if (app->use_gpu) {
    gkyl_cu_memcpy(omega_cfl_ho, fl->omegaCfl_ptr, sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    omega_cfl_ho[0] = fl->omegaCfl_ptr[0];
  }
  double omega_cfl = omega_cfl_ho[0];

  app->stat.species_omega_cfl_tm += gkyl_time_diff_now_sec(tm);

  return app->cfl / omega_cfl;
}

void
pkpm_species_calc_self_moms(
  gkyl_pkpm_app *app, struct pkpm_species *sp, struct gkyl_array *fin,
  const struct gkyl_array *fluidin
)
{
  // Every kernel of the stage sees first-harmonic blocks perpendicular to the stage's b: project
  // M0_perp and N_perp of the stage state at the Gauss-Legendre nodes of b.
  pkpm_species_project_transverse(app, sp, fin);
  // Configuration-space variables of the stage: the moments of the normalized distributions,
  // the flow and thermal speed, the kinetic coupling moments and T_perp/T.
  pkpm_species_calc_pkpm_vars(app, sp, fin, fluidin);
  // The distributions of the stage at the configuration-space nodes, for the LBO.
  pkpm_species_nodal_dist(app, sp, fin);
  // compute necessary moments and boundary corrections for collisions
  if (sp->kinetic->collision_id == GKYL_LBO_COLLISIONS) {
    pkpm_species_lbo_moms(app, sp, &sp->kinetic->lbo, fin);
  }
  // node values of the fluid state for the computed frequencies and cross-species collisions
  pkpm_species_coll_cross_prim(app, sp, fluidin);

  // Configuration-space face fluxes of the stage (kinetic, fluid and their rates, which start
  // here) and the particle loss rates they imply through absorbing walls.
  struct timespec tm = gkyl_wall_clock();
  gkyl_array_clear(sp->kinetic->cflrate, 0.0);
  gkyl_array_clear(sp->fluid->cflrate, 0.0);
  pkpm_species_face_flux(app, sp, fin);
  pkpm_species_calc_bflux(app, sp);
  sp->fluid->fluid_tm += gkyl_time_diff_now_sec(tm);

  // Explicit temperature source: its energy relaxation for the fluid right-hand side and the
  // rate of the distributions for the kinetic one.
  pkpm_species_source_bgk_calc_fluid(app, sp, &sp->kinetic->src_bgk, fluidin);
}

void
pkpm_species_calc_coupled_vars(
  gkyl_pkpm_app *app, struct pkpm_species *sp, const struct gkyl_array *fin,
  const struct gkyl_array *fluidin
)
{
  // rate of the kinetic operator with the partner species: needs the node values (and the
  // self-collision frequencies) of every species
  pkpm_species_coll_cross_moms(app, sp);
  // rate of the boundary-flux source (which reads another species' losses) and the source
  // distributions of the stage
  pkpm_species_source_calc_rate(app, sp, &sp->kinetic->src);
  pkpm_species_source_calc_dist(app, sp, &sp->kinetic->src);
}

double
pkpm_species_rhs(
  gkyl_pkpm_app *app, struct pkpm_species *sp, int stage, const struct gkyl_array *fin,
  const struct gkyl_array *fluidin, const struct gkyl_array *emin, struct gkyl_array *fout,
  struct gkyl_array *fluidout
)
{
  double dt = DBL_MAX;
  dt = fmin(dt, pkpm_kinetic_species_rhs(app, sp->kinetic, stage, fin, emin, fout));
  dt = fmin(dt, pkpm_fluid_species_rhs(app, sp->fluid, fluidin, emin, fluidout));
  return dt;
}

void
pkpm_species_step_f(
  struct pkpm_species *sp, double dt, const struct gkyl_array *fin,
  const struct gkyl_array *fluidin, struct gkyl_array *fout, struct gkyl_array *fluidout
)
{
  gkyl_array_accumulate(gkyl_array_scale(fout, dt), 1.0, fin);
  gkyl_array_accumulate(gkyl_array_scale(fluidout, dt), 1.0, fluidin);
}

void
pkpm_species_finish_stage(
  gkyl_pkpm_app *app, struct pkpm_species *sp, int stage, double dt, struct gkyl_array *fout,
  struct gkyl_array *fluidout
)
{
  // at the nodes the gyrophase mixing damps, the first-harmonic blocks are replaced by the
  // integrating-factor stage; the fluid state is limited
  pkpm_species_harm_exp_advance(app, sp, stage, dt, fout);
  pkpm_species_apply_positivity(app, sp->fluid, fluidout);
}

// Determine which directions are periodic and which directions are not periodic,
// and then apply boundary conditions for distribution function
void
pkpm_species_apply_bc(
  gkyl_pkpm_app *app, const struct pkpm_kinetic_species *kin, struct gkyl_array *f
)
{
  struct timespec wst = gkyl_wall_clock();

  int num_periodic_dir = app->num_periodic_dir, cdim = app->cdim;
  gkyl_comm_array_per_sync(
    kin->comm, &kin->local, &kin->local_ext_conf, num_periodic_dir, app->periodic_dirs, f
  );

  int is_np_bc[3] = {1, 1, 1}; // flags to indicate if direction is periodic
  for (int d = 0; d < num_periodic_dir; ++d) {
    is_np_bc[app->periodic_dirs[d]] = 0;
  }

  for (int d = 0; d < cdim; ++d) {
    if (is_np_bc[d]) {
      switch (kin->lower_bc[d]) {
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(kin->bc_lo[d], kin->bc_buffer, f);
          break;
        case GKYL_SPECIES_REFLECT:
          gkyl_bc_pkpm_species_reflect_advance(kin->bc_lo_reflect[d], app->field->bvar, f);
          break;
        case GKYL_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(kin->bc_lo[d], kin->bc_buffer_lo_fixed[d], f);
          break;
        case GKYL_SPECIES_NO_SLIP:
        case GKYL_SPECIES_WEDGE:
          assert(false);
          break;
        default:
          break;
      }

      switch (kin->upper_bc[d]) {
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(kin->bc_up[d], kin->bc_buffer, f);
          break;
        case GKYL_SPECIES_REFLECT:
          gkyl_bc_pkpm_species_reflect_advance(kin->bc_up_reflect[d], app->field->bvar, f);
          break;
        case GKYL_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(kin->bc_up[d], kin->bc_buffer_up_fixed[d], f);
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

  gkyl_comm_array_sync(kin->comm, &kin->local, &kin->local_ext_conf, f);

  app->stat.species_bc_tm += gkyl_time_diff_now_sec(wst);
}

// Determine which directions are periodic and which directions are not periodic,
// and then apply boundary conditions for the fluid state
void
pkpm_fluid_species_apply_bc(
  gkyl_pkpm_app *app, const struct pkpm_fluid_species *fl, struct gkyl_array *fluid
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
      switch (fl->lower_bc[d]) {
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(fl->bc_lo[d], fl->bc_buffer, fluid);
          break;
        case GKYL_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(fl->bc_lo[d], fl->bc_buffer_lo_fixed[d], fluid);
          break;
        case GKYL_SPECIES_NO_SLIP:
        case GKYL_SPECIES_WEDGE:
          assert(false);
          break;
        default:
          break;
      }

      switch (fl->upper_bc[d]) {
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(fl->bc_up[d], fl->bc_buffer, fluid);
          break;
        case GKYL_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(fl->bc_up[d], fl->bc_buffer_up_fixed[d], fluid);
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
    app->pos_map, &app->basis, &species->kinetic->local, species->kinetic->f,
    species->kinetic->f_no_J
  );
  gkyl_dg_calc_prod_op_range(
    app->basis, 0, species->kinetic->L2_f, 0, species->kinetic->f_no_J, 0, species->kinetic->f,
    species->kinetic->local
  );
  gkyl_array_scale_range(
    species->kinetic->L2_f, species->kinetic->grid.cellVolume, &species->kinetic->local
  );

  double L2[1] = {0.0};
  if (app->use_gpu) {
    gkyl_array_reduce_range(
      species->kinetic->red_L2_f, species->kinetic->L2_f, GKYL_SUM, &species->kinetic->local
    );
    gkyl_cu_memcpy(L2, species->kinetic->red_L2_f, sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    gkyl_array_reduce_range(L2, species->kinetic->L2_f, GKYL_SUM, &species->kinetic->local);
  }
  double L2_global[1] = {0.0};
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, L2, L2_global);

  gkyl_dynvec_append(species->kinetic->integ_L2_f, tm, L2_global);
}

void
pkpm_species_nodal_dist(
  gkyl_pkpm_app *app, const struct pkpm_species *species, const struct gkyl_array *fin
)
{
  if (species->kinetic->has_nodal_ops) {
    struct timespec tm = gkyl_wall_clock();
    gkyl_dg_conf_nodal_phase_to_nodal(
      app->conf_nodal, &species->kinetic->local, fin, species->kinetic->f_nodal
    );
    gkyl_array_clear_range(species->kinetic->rhs_nodal, 0.0, &species->kinetic->local);
    app->stat.species_coll_tm += gkyl_time_diff_now_sec(tm);
  }
}

void
pkpm_species_nodal_accumulate(
  gkyl_pkpm_app *app, const struct pkpm_kinetic_species *kin, struct gkyl_array *out
)
{
  if (kin->has_nodal_ops) {
    struct timespec tm = gkyl_wall_clock();
    gkyl_dg_conf_nodal_phase_to_modal_accumulate(app->conf_nodal, &kin->local, kin->rhs_nodal, out);
    app->stat.species_coll_tm += gkyl_time_diff_now_sec(tm);
  }
}

void
pkpm_species_tm(gkyl_pkpm_app *app)
{
  app->stat.species_rhs_tm = 0.0;
  app->stat.fluid_species_rhs_tm = 0.0;
  for (int i = 0; i < app->num_species; ++i) {
    app->stat.species_rhs_tm += app->species[i].kinetic->kinetic_tm;
    app->stat.fluid_species_rhs_tm += app->species[i].fluid->fluid_tm;
  }
}

// release resources for PKPM species
void
pkpm_species_release(const gkyl_pkpm_app *app, const struct pkpm_species *s)
{
  // release various arrays
  gkyl_array_release(s->kinetic->f);
  gkyl_array_release(s->kinetic->f_no_J);
  gkyl_array_release(s->kinetic->f1);
  gkyl_array_release(s->kinetic->fnew);
  if (s->kinetic->has_nodal_ops) {
    gkyl_array_release(s->kinetic->f_nodal);
    gkyl_array_release(s->kinetic->rhs_nodal);
  }
  gkyl_array_release(s->kinetic->blk_in);
  gkyl_array_release(s->kinetic->blk_out);
  gkyl_array_release(s->fluid->fluid);
  gkyl_array_release(s->fluid->fluid1);
  gkyl_array_release(s->fluid->fluidnew);

  gkyl_array_release(s->kinetic->cflrate);
  gkyl_array_release(s->fluid->cflrate);

  gkyl_array_release(s->kinetic->bc_buffer);
  gkyl_array_release(s->fluid->bc_buffer);
  for (int d = 0; d < app->cdim; ++d) {
    if (s->kinetic->bc_buffer_lo_fixed[d]) {
      gkyl_array_release(s->kinetic->bc_buffer_lo_fixed[d]);
    }
    if (s->kinetic->bc_buffer_up_fixed[d]) {
      gkyl_array_release(s->kinetic->bc_buffer_up_fixed[d]);
    }
    if (s->fluid->bc_buffer_lo_fixed[d]) {
      gkyl_array_release(s->fluid->bc_buffer_lo_fixed[d]);
    }
    if (s->fluid->bc_buffer_up_fixed[d]) {
      gkyl_array_release(s->fluid->bc_buffer_up_fixed[d]);
    }
  }

  gkyl_comm_release(s->kinetic->comm);

  if (app->use_gpu) {
    gkyl_array_release(s->kinetic->f_host);
    gkyl_array_release(s->kinetic->f_no_J_host);
    gkyl_array_release(s->fluid->fluid_host);
    gkyl_array_release(s->fluid->fluid_no_J_host);
    gkyl_array_release(s->fluid->pkpm_prim_io_host);
    gkyl_array_release(s->kinetic->pkpm_vars_io_host);
  }

  gkyl_array_release(s->fluid->qmem);

  // release moment data
  pkpm_species_moment_release(app, &s->kinetic->pkpm_moms_diag);
  pkpm_species_moment_release(app, &s->kinetic->pkpm_nmoms);

  // release pkpm arrays data and updaters
  gkyl_array_release(s->fluid->pkpm_u);
  gkyl_array_release(s->kinetic->pkpm_u);
  gkyl_array_release(s->kinetic->pkpm_grad_u);
  gkyl_array_release(s->kinetic->pkpm_theta);
  gkyl_array_release(s->kinetic->pkpm_alpha_0);
  gkyl_array_release(s->kinetic->pkpm_c_h);
  gkyl_array_release(s->fluid->pkpm_xi);
  gkyl_array_release(s->kinetic->pkpm_xi);
  gkyl_array_release(s->fluid->pkpm_xi3);
  gkyl_array_release(s->kinetic->pkpm_grad_xi);
  gkyl_array_release(s->fluid->pkpm_scratch);
  gkyl_array_release(s->fluid->pkpm_moms);
  gkyl_array_release(s->fluid->pkpm_moms_no_J);
  gkyl_array_release(s->fluid->fluid_no_J);
  gkyl_array_release(s->kinetic->pkpm_cmoms);
  if (app->use_gpu) {
    gkyl_array_release(s->kinetic->pkpm_cmoms_host);
  }
  gkyl_array_release(s->kinetic->pkpm_rhs_moms);
  gkyl_array_release(s->kinetic->conf_flux_surf);
  gkyl_array_release(s->kinetic->flux_blend);
  gkyl_array_release(s->kinetic->vel_flux_surf);
  gkyl_array_release(s->fluid->fluid_flux_surf);
  gkyl_array_release(s->fluid->pkpm_rho);
  gkyl_array_release(s->fluid->pkpm_rho3);
  gkyl_array_release(s->fluid->pkpm_mom3);
  if (s->fluid->has_diffusion) {
    gkyl_dg_pkpm_mom_diffusion_release(s->fluid->mom_diff);
  }
  gkyl_array_release(s->fluid->pkpm_grad_rho);
  gkyl_array_release(s->kinetic->pkpm_grad_ln_n);
  gkyl_array_release(s->kinetic->pkpm_c_0);
  gkyl_array_release(s->kinetic->pkpm_a_perp);
  if (app->use_gpu) {
    gkyl_array_release(s->kinetic->pkpm_closure_step_host);
  }
  gkyl_array_release(s->kinetic->pkpm_closure_step);
  gkyl_array_release(s->fluid->integ_mom);
  gkyl_array_release(s->fluid->pkpm_prim_io);
  gkyl_array_release(s->kinetic->pkpm_vars_io);
  gkyl_array_release(s->fluid->pkpm_pos_count);
  gkyl_array_release(s->fluid->pkpm_pos_io);
  if (app->use_gpu) {
    gkyl_array_release(s->fluid->pkpm_pos_io_host);
  }
  gkyl_dg_calc_pkpm_vars_release(s->fluid->calc_pkpm_vars);
  gkyl_dg_calc_pkpm_vars_release(s->kinetic->calc_pkpm_vars);
  gkyl_dg_pkpm_remap_flux_surf_release(s->kinetic->calc_face_flux);
  gkyl_dg_pkpm_vel_flux_surf_release(s->kinetic->calc_vel_flux);
  gkyl_dg_pkpm_harm_vol_release(s->kinetic->harm_vol);
  gkyl_array_release(s->kinetic->harm_coef);
  gkyl_array_release(s->kinetic->harm_wnode);
  gkyl_dg_pkpm_harm_exp_release(s->kinetic->harm_exp);
  gkyl_array_release(s->kinetic->harm_exp_coef);
  gkyl_array_release(s->kinetic->harm_exp_nu);
  gkyl_array_release(s->kinetic->harm_exp_act);
  gkyl_array_release(s->kinetic->harm_exp_col);
  gkyl_hyper_dg_release(s->kinetic->slvr_kinetic_conf);
  gkyl_hyper_dg_release(s->kinetic->slvr_kinetic_vel);
  gkyl_hyper_dg_release(s->fluid->slvr_fluid);
  gkyl_dg_eqn_release(s->kinetic->eqn_kinetic_conf);
  gkyl_dg_eqn_release(s->kinetic->eqn_kinetic_vel);
  gkyl_dg_eqn_release(s->fluid->eqn_fluid);

  gkyl_array_release(s->kinetic->L2_f);
  gkyl_dynvec_release(s->kinetic->integ_L2_f);
  gkyl_dynvec_release(s->fluid->integ_diag);
  gkyl_dynvec_release(s->fluid->pos_stat);
  gkyl_dynvec_release(s->kinetic->bflux_diag);
  if (s->fluid->fluid_flux_surf_host) {
    gkyl_array_release(s->fluid->fluid_flux_surf_host);
  }
  pkpm_species_source_release(app, &s->kinetic->src);
  pkpm_species_source_bgk_release(app, &s->kinetic->src_bgk);
  pkpm_species_coll_cross_release(app, &s->kinetic->cross);
  if (s->kinetic->has_stage_weights) {
    gkyl_array_release(s->kinetic->stage_lam);
    gkyl_pkpm_stage_weights_release(s->kinetic->up_stage_weights);
  }

  gkyl_array_release(s->fluid->app_accel);
  if (s->fluid->has_app_accel) {
    if (app->use_gpu) {
      gkyl_array_release(s->fluid->app_accel_host);
    }
    gkyl_proj_on_basis_release(s->fluid->app_accel_proj);
  }

  if (s->kinetic->collision_id == GKYL_BGK_COLLISIONS) {
    pkpm_species_bgk_release(app, &s->kinetic->bgk);
  }
  if (s->kinetic->collision_id == GKYL_LBO_COLLISIONS) {
    pkpm_species_lbo_release(app, &s->kinetic->lbo);
  }

  // Copy BCs are allocated by default. Need to free.
  for (int d = 0; d < app->cdim; ++d) {
    gkyl_bc_basic_release(s->kinetic->bc_lo[d]);
    gkyl_bc_basic_release(s->kinetic->bc_up[d]);
    if (s->kinetic->bc_lo_reflect[d]) {
      gkyl_bc_pkpm_species_reflect_release(s->kinetic->bc_lo_reflect[d]);
    }
    if (s->kinetic->bc_up_reflect[d]) {
      gkyl_bc_pkpm_species_reflect_release(s->kinetic->bc_up_reflect[d]);
    }
    gkyl_bc_basic_release(s->fluid->bc_lo[d]);
    gkyl_bc_basic_release(s->fluid->bc_up[d]);
    if (s->kinetic->bc_lo_theta[d]) {
      gkyl_bc_basic_release(s->kinetic->bc_lo_theta[d]);
    }
    if (s->kinetic->bc_up_theta[d]) {
      gkyl_bc_basic_release(s->kinetic->bc_up_theta[d]);
    }
  }

  if (app->use_gpu) {
    gkyl_cu_free(s->kinetic->omega_cfl);
    gkyl_cu_free(s->fluid->omegaCfl_ptr);
    gkyl_cu_free(s->kinetic->red_L2_f);
    gkyl_cu_free(s->fluid->red_integ_diag);
  } else {
    gkyl_free(s->kinetic->omega_cfl);
    gkyl_free(s->fluid->omegaCfl_ptr);
  }

  if (s->fluid->src_bgk_rhs) {
    gkyl_array_release(s->fluid->src_bgk_rhs);
  }

  gkyl_free(s->kinetic);
  gkyl_free(s->fluid);
}
