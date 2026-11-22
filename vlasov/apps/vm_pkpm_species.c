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
#include <gkyl_vlasov_priv.h>

#include <assert.h>
#include <time.h>

// Construction route of a PKPM species (GKYL_SPECIES_PKPM): a species container with both
// aspects. This file allocates and initializes the two aspects together, sets their methods
// and those of the container, and holds the staging phases of a stage, in which the order of
// the couplings between the distributions, the fluid state and the field is written:
//   vm_pkpm_species_calc_self_moms     what reads this species and the field alone
//   vm_pkpm_species_calc_coupled_vars  what reads other species
//   vm_pkpm_kinetic_species_rhs, vm_pkpm_fluid_species_rhs   the two right-hand sides
//   vm_pkpm_species_finish_stage       integrating factor of the first harmonic, limiter
// The collision operators, the sources and the moments are in vm_pkpm_species_*.c; the implicit
// operators of a stage in vlasov_update_implicit.c.

// Reflecting wall of the distributions on the lower or upper edge of direction dir.
static struct gkyl_bc_pkpm_species_reflect *
vm_pkpm_species_reflect_new(
  gkyl_vlasov_app *app, struct vlasov_species *s, int dir, enum gkyl_edge_loc edge
)
{
  return gkyl_bc_pkpm_species_reflect_inew(&(struct gkyl_bc_pkpm_species_reflect_inp){
    .dir = dir,
    .edge = edge,
    .cdim = app->cdim,
    .phase_grid = &s->kinetic->grid,
    .phase_basis = &s->kinetic->basis,
    .phase_basis_on_dev = s->kinetic->basis_on_dev,
    .conf_num_basis = app->basis.num_basis,
    .num_blocks = PKPM_DIST_NUM_BLOCKS,
    .skin_r = edge == GKYL_LOWER_EDGE ? &s->kinetic->lower_skin[dir] : &s->kinetic->upper_skin[dir],
    .ghost_r = edge == GKYL_LOWER_EDGE ? &s->kinetic->lower_ghost[dir] :
                                         &s->kinetic->upper_ghost[dir],
    .phase_range_ext = &s->kinetic->local_ext,
    .conf_range_ext = &app->local_ext,
    .use_gpu = app->use_gpu,
  });
}

// The boundary conditions of a PKPM species, in either block.
static bool
vm_pkpm_species_bc_is_supported(enum gkyl_species_bc_type bc_type)
{
  return bc_type == GKYL_SPECIES_COPY || bc_type == GKYL_SPECIES_ABSORB ||
         bc_type == GKYL_SPECIES_REFLECT || bc_type == GKYL_SPECIES_FIXED_FUNC;
}

// The basic boundary condition of the distributions for a kinetic BC type. A reflecting wall
// uses the field-aligned reflection updater instead; the basic updater is not applied there.
static enum gkyl_bc_basic_type
vm_pkpm_species_bc_basic_type_dist(enum gkyl_species_bc_type bc_type)
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
vm_pkpm_species_bc_basic_type_fluid(enum gkyl_species_bc_type bc_type)
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

// Methods defined below.
static void vm_pkpm_kinetic_species_apply_ic(
  gkyl_vlasov_app *app, struct vm_species *kin, double t0
);
static void vm_pkpm_fluid_species_apply_ic(
  gkyl_vlasov_app *app, struct vm_fluid_species *fl, double t0
);
static void vm_pkpm_kinetic_species_calc_L2(gkyl_vlasov_app *app, struct vm_species *kin, double tm);
static void vm_pkpm_kinetic_species_write_L2(gkyl_vlasov_app *app, struct vm_species *kin);
static struct gkyl_app_restart_status vm_pkpm_kinetic_species_from_file(
  gkyl_vlasov_app *app, struct vm_species *kin, const char *fname
);
static struct gkyl_app_restart_status vm_pkpm_fluid_species_from_file(
  gkyl_vlasov_app *app, struct vm_fluid_species *fl, const char *fname
);
static void vm_pkpm_species_calc_integrated_mom(
  gkyl_vlasov_app *app, struct vlasov_species *s, double tm
);
static void vm_pkpm_species_write(
  gkyl_vlasov_app *app, struct vlasov_species *s, double tm, int frame
);
static void vm_pkpm_species_write_mom(
  gkyl_vlasov_app *app, struct vlasov_species *s, double tm, int frame
);
static void vm_pkpm_species_write_integrated_mom(gkyl_vlasov_app *app, struct vlasov_species *s);
static double vm_pkpm_species_wall_loss_rate(gkyl_vlasov_app *app, const struct vlasov_species *s);
static struct gkyl_app_restart_status vm_pkpm_species_read_from_frame(
  gkyl_vlasov_app *app, struct vlasov_species *s, int frame
);

// Methods of the two aspects that act on one array each.
static void
vm_pkpm_species_array_step_f(struct gkyl_array *out, double dt, const struct gkyl_array *inp)
{
  gkyl_array_accumulate(gkyl_array_scale(out, dt), 1.0, inp);
}

static void
vm_pkpm_species_array_combine(
  struct gkyl_array *out, double c1, const struct gkyl_array *arr1, double c2,
  const struct gkyl_array *arr2, const struct gkyl_range *rng
)
{
  array_combine(out, c1, arr1, c2, arr2, rng);
}

static void
vm_pkpm_species_array_copy(
  struct gkyl_array *out, const struct gkyl_array *inp, const struct gkyl_range *range
)
{
  gkyl_array_copy_range(out, inp, range);
}

// The sources of a PKPM species are inside its right-hand side, before the closures: nothing in
// the source phases of the stage.
static void
vm_pkpm_kinetic_species_no_source_moms(
  gkyl_vlasov_app *app, struct vm_species *kin, const struct gkyl_array *fin
)
{
}

static void
vm_pkpm_kinetic_species_no_rhs_source(
  gkyl_vlasov_app *app, struct vm_species *kin, double tcurr, struct gkyl_array *rhs
)
{
}

static void
vm_pkpm_fluid_species_no_rhs_source(
  gkyl_vlasov_app *app, struct vm_fluid_species *fl, struct gkyl_array *rhs
)
{
}

// The positivity limiter of the fluid state is applied in the finish-stage method, before the
// boundary conditions: nothing after them.
static void
vm_pkpm_fluid_species_no_limiter(
  gkyl_vlasov_app *app, struct vm_fluid_species *fl, struct gkyl_array *fluid
)
{
}

static void
vm_pkpm_kinetic_species_apply_bc(
  gkyl_vlasov_app *app, const struct vm_species *kin, struct gkyl_array *f, double tcurr
)
{
  vm_pkpm_species_apply_bc(app, kin, f);
}

// No explicit field coupling: the implicit coupling takes the current in the source step.
static void
vm_pkpm_species_no_field_coupling(
  gkyl_vlasov_app *app, struct vlasov_species *sp, const struct gkyl_array *fin,
  const struct gkyl_array *fluidin, struct gkyl_array *target
)
{
}

// Current of the fluid aspect coupled explicitly to the Maxwell field ('target' is emout): the
// fluid array is [rho, rho ux, rho uy, rho uz, E], and its momentum block over the mass is the
// number flux.
static void
vm_pkpm_species_accumulate_current(
  gkyl_vlasov_app *app, struct vlasov_species *sp, const struct gkyl_array *fin,
  const struct gkyl_array *fluidin, struct gkyl_array *target
)
{
  struct timespec wst = gkyl_wall_clock();
  double qbymeps = sp->charge / (sp->mass * app->field->info.epsilon0);
  gkyl_array_set_offset(sp->fluid->pkpm.mom3, 1.0, fluidin, app->basis.num_basis);
  gkyl_array_accumulate_range(target, -qbymeps, sp->fluid->pkpm.mom3, &app->local);
  vm_field_accumulate_ghost_current(app, app->field, qbymeps, sp->fluid->pkpm.mom3, target);
  app->stat.current_tm += gkyl_time_diff_now_sec(wst);
}

void
vlasov_pkpm_species_init(
  struct gkyl_vlasov_app *app, const struct gkyl_vlasov_species *inp, struct vlasov_species *sp
)
{
  *sp = (struct vlasov_species){};
  sp->type = GKYL_SPECIES_PKPM;

  strcpy(sp->name, inp->name);
  sp->charge = inp->charge;
  sp->mass = inp->mass;

  // Each aspect holds its input block and a copy of the identity.
  sp->kinetic = gkyl_malloc(sizeof(struct vm_species));
  *sp->kinetic = (struct vm_species){};
  sp->kinetic->info = inp->kinetic;
  strcpy(sp->kinetic->name, inp->name);
  sp->kinetic->charge = inp->charge;
  sp->kinetic->mass = inp->mass;

  sp->fluid = gkyl_malloc(sizeof(struct vm_fluid_species));
  *sp->fluid = (struct vm_fluid_species){};
  sp->fluid->info = inp->fluid;
  strcpy(sp->fluid->name, inp->name);
  sp->fluid->charge = inp->charge;
  sp->fluid->mass = inp->mass;

  sp->coll_cross = &sp->kinetic->pkpm.cross;
  sp->calc_self_moms_func = vm_pkpm_species_calc_self_moms;
  sp->calc_coupled_vars_func = vm_pkpm_species_calc_coupled_vars;
  sp->finish_stage_func = vm_pkpm_species_finish_stage;
  sp->complete_step_func = vm_pkpm_species_complete_step;
  // diagnostics, output and restart of the species as a whole
  sp->calc_integrated_mom_func = vm_pkpm_species_calc_integrated_mom;
  sp->write_func = vm_pkpm_species_write;
  sp->write_mom_func = vm_pkpm_species_write_mom;
  sp->write_integrated_mom_func = vm_pkpm_species_write_integrated_mom;
  sp->read_from_frame_func = vm_pkpm_species_read_from_frame;
  // Field coupling of the stage: the current of the fluid momentum when the Maxwell field is
  // coupled explicitly; the implicit coupling takes it in the source step.
  sp->accumulate_field_coupling_func = app->use_explicit_sources ?
                                         vm_pkpm_species_accumulate_current :
                                         vm_pkpm_species_no_field_coupling;
  sp->wall_loss_rate_func = vm_pkpm_species_wall_loss_rate;
  sp->release_func = vm_pkpm_species_release;
}

void
vm_pkpm_species_init(struct gkyl_vm *pkpm, struct gkyl_vlasov_app *app, struct vlasov_species *s)
{
  // One velocity dimension, the normalized parallel velocity w = v_par/xi (app->vdim is the
  // velocity dimension of the kinetic species of the run).
  int cdim = app->cdim, vdim = 1;
  int pdim = cdim + vdim;

  // Phase-space basis: the hybrid basis (p=1 in configuration space, p=2 in v_par).
  gkyl_cart_modal_hybrid(&s->kinetic->basis, cdim, vdim);
  gkyl_cart_modal_serendip(&s->kinetic->basis_vel, vdim, 2);
  if (app->use_gpu) {
    s->kinetic->basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
    gkyl_cart_modal_hybrid_cu_dev(s->kinetic->basis_on_dev, cdim, vdim);
  } else {
    s->kinetic->basis_on_dev = &s->kinetic->basis;
  }

  // A PKPM species owns both blocks; the distributions start from one function of eight blocks.
  assert(s->type == GKYL_SPECIES_PKPM);
  assert(s->kinetic->info.num_init == 1);
  assert(s->kinetic->info.projection[0].proj_id == GKYL_PROJ_FUNC);
  assert(s->kinetic->info.projection[0].func && s->fluid->info.init);
  // defaults of the number of face quadrature nodes and of the fluid positivity floor
  s->kinetic->pkpm.num_face_nodes =
    s->kinetic->info.num_face_nodes > 0 ? s->kinetic->info.num_face_nodes : 5;
  assert(s->kinetic->pkpm.num_face_nodes >= 3);
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
  gkyl_sub_range_init(
    &s->kinetic->pkpm.local_ext_conf, &s->kinetic->local_ext, lower_conf, upper_conf
  );

  // allocate distribution function arrays: the normalized distributions
  // [F, G, M0_perp (x, y, z), N_perp (x, y, z)] (enum gkyl_pkpm_dist_block)
  int nc_dist = PKPM_DIST_NUM_BLOCKS * s->kinetic->basis.num_basis;
  s->kinetic->f = mkarr(app->use_gpu, nc_dist, s->kinetic->local_ext.volume);
  s->kinetic->f1 = mkarr(app->use_gpu, nc_dist, s->kinetic->local_ext.volume);
  // physical distributions (position-map Jacobian divided out) for the output and the L2 norm
  s->kinetic->f_no_J = mkarr(app->use_gpu, nc_dist, s->kinetic->local_ext.volume);
  s->kinetic->fnew = mkarr(app->use_gpu, nc_dist, s->kinetic->local_ext.volume);
  gkyl_array_clear(s->kinetic->f, 0.0);
  gkyl_array_clear(s->kinetic->f1, 0.0);
  gkyl_array_clear(s->kinetic->f_no_J, 0.0);
  gkyl_array_clear(s->kinetic->fnew, 0.0);
  s->kinetic->pkpm.blk_in =
    mkarr(app->use_gpu, s->kinetic->basis.num_basis, s->kinetic->local_ext.volume);
  s->kinetic->pkpm.blk_out =
    mkarr(app->use_gpu, s->kinetic->basis.num_basis, s->kinetic->local_ext.volume);
  // allocate fluid arrays [rho, rho ux, rho uy, rho uz, E]
  s->fluid->fluid = mkarr(app->use_gpu, 5 * app->basis.num_basis, app->local_ext.volume);
  s->fluid->fluid1 = mkarr(app->use_gpu, 5 * app->basis.num_basis, app->local_ext.volume);
  s->fluid->fluidnew = mkarr(app->use_gpu, 5 * app->basis.num_basis, app->local_ext.volume);

  // create host arrays if on GPUs so initialization occurs host-side
  s->kinetic->f_host = s->kinetic->f;
  s->kinetic->pkpm.f_no_J_host = s->kinetic->f_no_J;
  s->fluid->fluid_host = s->fluid->fluid;
  if (app->use_gpu) {
    s->kinetic->f_host = mkarr(false, nc_dist, s->kinetic->local_ext.volume);
    s->kinetic->pkpm.f_no_J_host = mkarr(false, nc_dist, s->kinetic->local_ext.volume);
    s->fluid->fluid_host = mkarr(false, 5 * app->basis.num_basis, app->local_ext.volume);
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
  s->fluid->qmem = mkarr(app->use_gpu, 8 * app->basis.num_basis, app->local_ext.volume);

  int nc = app->basis.num_basis;
  // pkpm moments for diagnostics [rho, M1, p_par, p_perp, q_par, q_perp, r_parpar, r_parperp]
  vm_pkpm_species_moment_init(app, s, &s->kinetic->pkpm.moms_diag, true);
  // raw moments of the normalized distributions [F, G, M0_perp, N_perp] (mass 1), the kinetic
  // coupling moments [p_par, p_perp, q_par, q_perp, s (3), r (3)] formed from them and the
  // constraints
  vm_pkpm_species_moment_init_mass(app, s, &s->kinetic->pkpm.nmoms, false, 1.0);
  s->fluid->pkpm.moms = mkarr(app->use_gpu, 10 * nc, app->local_ext.volume);
  s->fluid->pkpm.moms_no_J = mkarr(app->use_gpu, 10 * nc, app->local_ext.volume);
  s->fluid->pkpm.fluid_no_J = mkarr(app->use_gpu, 5 * nc, app->local_ext.volume);
  s->fluid->pkpm.fluid_no_J_host = s->fluid->pkpm.fluid_no_J;
  if (app->use_gpu) {
    s->fluid->pkpm.fluid_no_J_host = mkarr(false, 5 * nc, app->local_ext.volume);
  }
  gkyl_array_clear(s->fluid->pkpm.moms, 0.0);
  s->kinetic->pkpm.cmoms = mkarr(app->use_gpu, 6 * nc, app->local_ext.volume);
  s->kinetic->pkpm.cmoms_host = s->kinetic->pkpm.cmoms;
  if (app->use_gpu) {
    s->kinetic->pkpm.cmoms_host = mkarr(false, 6 * nc, app->local_ext.volume);
  }

  // Configuration-space variables of the update (p=1).
  s->fluid->pkpm.u = mkarr(app->use_gpu, 3 * nc, app->local_ext.volume);
  s->kinetic->pkpm.u = mkarr(app->use_gpu, 3 * nc, app->local_ext.volume);
  s->kinetic->pkpm.grad_u = mkarr(app->use_gpu, cdim * 3 * nc, app->local_ext.volume);
  gkyl_array_clear(s->kinetic->pkpm.grad_u, 0.0);
  s->kinetic->pkpm.theta = mkarr(app->use_gpu, nc, app->local_ext.volume);
  s->fluid->pkpm.rho = mkarr(app->use_gpu, nc, app->local_ext.volume);
  s->fluid->pkpm.rho3 = mkarr(app->use_gpu, 3 * nc, app->local_ext.volume);
  gkyl_array_clear(s->fluid->pkpm.rho3, 0.0);
  s->fluid->pkpm.mom3 = mkarr(app->use_gpu, 3 * nc, app->local_ext.volume);
  s->fluid->pkpm.grad_rho = mkarr(app->use_gpu, cdim * 3 * nc, app->local_ext.volume);
  gkyl_array_clear(s->fluid->pkpm.grad_rho, 0.0);
  s->kinetic->pkpm.grad_ln_n = mkarr(app->use_gpu, cdim * nc, app->local_ext.volume);
  gkyl_array_clear(s->kinetic->pkpm.grad_ln_n, 0.0);
  s->kinetic->pkpm.c_0 = mkarr(app->use_gpu, nc, app->local_ext.volume);
  gkyl_array_clear(s->kinetic->pkpm.c_0, 0.0);
  s->kinetic->pkpm.alpha_0 = mkarr(app->use_gpu, nc, app->local_ext.volume);
  gkyl_array_clear(s->kinetic->pkpm.alpha_0, 0.0);
  s->kinetic->pkpm.c_h = mkarr(app->use_gpu, nc, app->local_ext.volume);
  gkyl_array_clear(s->kinetic->pkpm.c_h, 0.0);
  s->kinetic->pkpm.a_perp = mkarr(app->use_gpu, 3 * nc, app->local_ext.volume);
  gkyl_array_clear(s->kinetic->pkpm.a_perp, 0.0);
  s->kinetic->pkpm.closure_step = mkarr(app->use_gpu, 2 * nc, app->local_ext.volume);
  gkyl_array_clear(s->kinetic->pkpm.closure_step, 0.0);
  s->kinetic->pkpm.closure_step_host = s->kinetic->pkpm.closure_step;
  if (app->use_gpu) {
    s->kinetic->pkpm.closure_step_host = mkarr(false, 2 * nc, app->local_ext.volume);
  }
  s->fluid->pkpm.xi = mkarr(app->use_gpu, 3 * nc, app->local_ext.volume);
  gkyl_array_clear(s->fluid->pkpm.xi, 0.0);
  s->kinetic->pkpm.xi = mkarr(app->use_gpu, 3 * nc, app->local_ext.volume);
  gkyl_array_clear(s->kinetic->pkpm.xi, 0.0);
  s->fluid->pkpm.xi3 = mkarr(app->use_gpu, 3 * nc, app->local_ext.volume);
  gkyl_array_clear(s->fluid->pkpm.xi3, 0.0);
  s->kinetic->pkpm.grad_xi = mkarr(app->use_gpu, cdim * 3 * nc, app->local_ext.volume);
  gkyl_array_clear(s->kinetic->pkpm.grad_xi, 0.0);
  s->fluid->pkpm.scratch = mkarr(app->use_gpu, nc, app->local_ext.volume);
  s->kinetic->pkpm.rhs_moms = mkarr(
    app->use_gpu, gkyl_dg_updater_moment_pkpm_num_mom(s->kinetic->pkpm.nmoms.mcalc) * nc,
    app->local_ext.volume
  );
  gkyl_array_clear(s->kinetic->pkpm.theta, 0.0);
  s->kinetic->pkpm.num_null_closure = 0;
  s->kinetic->pkpm.closure_cond_min = 1.0;
  s->kinetic->pkpm.closure_null_warned = false;

  // Updaters for the configuration-space variables, one per aspect: the fluid one forms the
  // primitive variables, their gradients, the limiter, the sources and the diagnostics, the
  // kinetic one the closures, the constraints and the transverse projection. Every method takes
  // the range it acts on.
  struct gkyl_dg_calc_pkpm_vars_inp calc_pkpm_vars_inp = {
    .conf_grid = &app->grid,
    .conf_basis = &app->basis,
    .pos_map = app->pos_map,
    .use_gpu = app->use_gpu,
  };
  s->fluid->pkpm.calc_pkpm_vars = gkyl_dg_calc_pkpm_vars_inew(&calc_pkpm_vars_inp);
  s->kinetic->pkpm.calc_pkpm_vars = gkyl_dg_calc_pkpm_vars_inew(&calc_pkpm_vars_inp);

  // array for storing integrated fluid variables in each cell
  // integ_mom : integral[rho, rhoux, rhouy, rhouz, rhoux^2, rhouy^2, rhouz^2, p_par, p_perp]
  s->fluid->integ_mom = mkarr(app->use_gpu, 9, app->local_ext.volume);
  // arrays for I/O, prim_io and vars_io
  // prim_io : [P_xx, P_xy, P_xz, P_yy, P_yz, P_zz, q_x, q_y, q_z]
  // vars_io : [ux, uy, uz, T_perp/m, a_P, div(b), bb:grad(u), T_perp/m div(b), xi, c]
  s->fluid->pkpm.prim_io = mkarr(app->use_gpu, 9 * nc, app->local_ext.volume);
  s->kinetic->pkpm.vars_io = mkarr(app->use_gpu, 10 * nc, app->local_ext.volume);
  s->fluid->pkpm.prim_io_host = s->fluid->pkpm.prim_io;
  s->kinetic->pkpm.vars_io_host = s->kinetic->pkpm.vars_io;
  if (app->use_gpu) {
    s->fluid->pkpm.prim_io_host = mkarr(false, 9 * nc, app->local_ext.volume);
    s->kinetic->pkpm.vars_io_host = mkarr(false, 10 * nc, app->local_ext.volume);
  }
  // Record of the positivity limiter: per cell [density scalings, energy scalings,
  // inadmissible averages] accumulated over the run, and its p=1 image for IO.
  s->fluid->pkpm.pos_count = mkarr(app->use_gpu, PKPM_POS_NUM_COUNT, app->local_ext.volume);
  gkyl_array_clear(s->fluid->pkpm.pos_count, 0.0);
  s->fluid->pkpm.pos_count_saved = mkarr(app->use_gpu, PKPM_POS_NUM_COUNT, app->local_ext.volume);
  s->fluid->pkpm.pos_io = mkarr(app->use_gpu, PKPM_POS_NUM_COUNT * nc, app->local_ext.volume);
  s->fluid->pkpm.pos_io_host = s->fluid->pkpm.pos_io;
  if (app->use_gpu) {
    s->fluid->pkpm.pos_io_host = mkarr(false, PKPM_POS_NUM_COUNT * nc, app->local_ext.volume);
  }

  // The configuration-space face flux updater and the stored modal fluxes (each cell owns its
  // lower faces, and stores the flux each of the two sides of the face receives).
  struct gkyl_dg_pkpm_remap_flux_surf_inp face_flux_inp = {
    .phase_grid = &s->kinetic->grid,
    .conf_basis = &app->basis,
    .phase_basis = &s->kinetic->basis,
    .pos_map = app->pos_map,
    .num_face_nodes = s->kinetic->pkpm.num_face_nodes,
    .use_gpu = app->use_gpu,
  };
  s->kinetic->pkpm.calc_face_flux = gkyl_dg_pkpm_remap_flux_surf_inew(&face_flux_inp);
  int nkc = gkyl_dg_pkpm_remap_flux_surf_num_surf_basis(s->kinetic->pkpm.calc_face_flux);
  s->kinetic->pkpm.conf_flux_surf =
    mkarr(app->use_gpu, 2 * PKPM_DIST_NUM_BLOCKS * cdim * nkc, s->kinetic->local_ext.volume);
  s->kinetic->pkpm.flux_blend = mkarr(app->use_gpu, cdim, app->local_ext.volume);
  gkyl_array_clear(s->kinetic->pkpm.flux_blend, 0.0);

  struct gkyl_dg_pkpm_vel_flux_surf_inp vel_flux_inp = {
    .phase_grid = &s->kinetic->grid,
    .conf_basis = &app->basis,
    .phase_basis = &s->kinetic->basis,
    .use_gpu = app->use_gpu,
  };
  s->kinetic->pkpm.calc_vel_flux = gkyl_dg_pkpm_vel_flux_surf_inew(&vel_flux_inp);
  int nkv = gkyl_dg_pkpm_vel_flux_surf_num_surf_basis(s->kinetic->pkpm.calc_vel_flux);
  s->kinetic->pkpm.vel_flux_surf =
    mkarr(app->use_gpu, PKPM_DIST_NUM_BLOCKS * nkv, s->kinetic->local_ext.volume);

  int nkf = gkyl_dg_pkpm_remap_flux_surf_num_surf_basis_fluid(s->kinetic->pkpm.calc_face_flux);
  s->fluid->pkpm.fluid_flux_surf = mkarr(app->use_gpu, 5 * cdim * nkf, app->local_ext.volume);
  s->fluid->pkpm.fluid_flux_surf_host =
    app->use_gpu ? mkarr(false, 5 * cdim * nkf, app->local_ext.volume) : 0;

  // Equation objects and hyper_dg solvers of the three passes.
  struct gkyl_dg_pkpm_kinetic_conf_inp kinetic_conf_inp = {
    .conf_basis = &app->basis,
    .phase_basis = &s->kinetic->basis,
    .conf_range = &app->local,
    .phase_range = &s->kinetic->local,
    .pos_map = app->pos_map,
    .bvar = app->field->bvar,
    .u = s->kinetic->pkpm.u,
    .xi = s->kinetic->pkpm.xi,
    .conf_flux_surf = s->kinetic->pkpm.conf_flux_surf,
    .use_gpu = app->use_gpu,
  };
  s->kinetic->pkpm.eqn_kinetic_conf = gkyl_dg_pkpm_kinetic_conf_inew(&kinetic_conf_inp);
  struct gkyl_dg_pkpm_harm_vol_inp harm_vol_inp = {
    .phase_grid = &s->kinetic->grid,
    .conf_basis = &app->basis,
    .phase_basis = &s->kinetic->basis,
    .qbym = s->charge / s->mass,
    .pos_map = app->pos_map,
    .use_gpu = app->use_gpu,
  };
  s->kinetic->pkpm.harm_vol = gkyl_dg_pkpm_harm_vol_inew(&harm_vol_inp);
  s->kinetic->pkpm.harm_coef = mkarr(
    app->use_gpu, gkyl_dg_pkpm_harm_vol_num_coef(s->kinetic->pkpm.harm_vol), app->local_ext.volume
  );
  s->kinetic->pkpm.harm_wnode = mkarr(
    app->use_gpu, gkyl_dg_pkpm_harm_vol_num_wnode(s->kinetic->pkpm.harm_vol), app->local_ext.volume
  );
  // gyrophase mixing: default rate 10 |Omega|, switched off by a negative input
  double gyro_mix = s->kinetic->info.gyro_mix == 0.0 ? 10.0 : s->kinetic->info.gyro_mix;
  s->kinetic->pkpm.harm_exp = gkyl_dg_pkpm_harm_exp_inew(&(struct gkyl_dg_pkpm_harm_exp_inp){
    .phase_grid = &s->kinetic->grid,
    .conf_basis = &app->basis,
    .qbym = s->charge / s->mass,
    .mix = gyro_mix,
    .use_gpu = app->use_gpu,
  });
  s->kinetic->pkpm.harm_exp_coef = mkarr(
    app->use_gpu, gkyl_dg_pkpm_harm_exp_num_coef(s->kinetic->pkpm.harm_exp), app->local_ext.volume
  );
  s->kinetic->pkpm.harm_exp_nu = mkarr(
    app->use_gpu, gkyl_dg_pkpm_harm_exp_num_rate(s->kinetic->pkpm.harm_exp),
    s->kinetic->local_ext.volume
  );
  s->kinetic->pkpm.harm_exp_act = mkarr(app->use_gpu, 1, app->local_ext.volume);
  s->kinetic->pkpm.harm_exp_col = mkarr(
    app->use_gpu, gkyl_dg_pkpm_harm_exp_num_col(s->kinetic->pkpm.harm_exp), app->local_ext.volume
  );
  struct gkyl_dg_pkpm_kinetic_vel_inp kinetic_vel_inp = {
    .conf_basis = &app->basis,
    .phase_basis = &s->kinetic->basis,
    .conf_range = &app->local,
    .phase_range = &s->kinetic->local,
    .c_0 = s->kinetic->pkpm.c_0,
    .alpha_0 = s->kinetic->pkpm.alpha_0,
    .c_h = s->kinetic->pkpm.c_h,
    .a_perp = s->kinetic->pkpm.a_perp,
    .vel_flux_surf = s->kinetic->pkpm.vel_flux_surf,
    .use_gpu = app->use_gpu,
  };
  s->kinetic->pkpm.eqn_kinetic_vel = gkyl_dg_pkpm_kinetic_vel_inew(&kinetic_vel_inp);
  struct gkyl_dg_pkpm_fluid_inp fluid_inp = {
    .conf_basis = &app->basis,
    .conf_range = &app->local,
    .pos_map = app->pos_map,
    .vlasov_pkpm_moms = s->fluid->pkpm.moms,
    .bvar = app->field->bvar,
    .u = s->fluid->pkpm.u,
    .fluid_flux_surf = s->fluid->pkpm.fluid_flux_surf,
    .use_gpu = app->use_gpu,
  };
  s->fluid->pkpm.eqn_fluid = gkyl_dg_pkpm_fluid_inew(&fluid_inp);

  int up_dirs_conf[GKYL_MAX_DIM] = {0, 1, 2}, zero_flux_conf[2 * GKYL_MAX_DIM] = {0};
  s->kinetic->pkpm.slvr_kinetic_conf = gkyl_hyper_dg_new(
    &s->kinetic->grid, &s->kinetic->basis, s->kinetic->pkpm.eqn_kinetic_conf, cdim, up_dirs_conf,
    zero_flux_conf, 1, app->use_gpu
  );
  int up_dirs_vel[GKYL_MAX_DIM] = {cdim}, zero_flux_vel[2 * GKYL_MAX_DIM] = {0};
  // The lower velocity edge has no ghost cell (boundary term with the lowest cell's own
  // flux slot); the upper edge is the regular face between the uppermost cell and the
  // velocity ghost cell, whose slot the closure fills.
  zero_flux_vel[cdim] = 1;
  s->kinetic->pkpm.slvr_kinetic_vel = gkyl_hyper_dg_new(
    &s->kinetic->grid, &s->kinetic->basis, s->kinetic->pkpm.eqn_kinetic_vel, 1, up_dirs_vel,
    zero_flux_vel, 1, app->use_gpu
  );
  s->fluid->pkpm.slvr_fluid = gkyl_hyper_dg_new(
    &app->grid, &app->basis, s->fluid->pkpm.eqn_fluid, cdim, up_dirs_conf, zero_flux_conf, 1,
    app->use_gpu
  );
  // set by the source and the temperature source where the species has them
  s->fluid->pkpm.src_scale = 0.0;
  s->fluid->pkpm.src_fluid = 0;
  s->fluid->pkpm.has_src_bgk_rhs = false;
  s->fluid->pkpm.src_bgk_rhs = 0;

  // array for storing F^2 in each cell
  s->kinetic->L2_f = mkarr(app->use_gpu, 1, s->kinetic->local_ext.volume);
  if (app->use_gpu) {
    s->kinetic->red_L2_f = gkyl_cu_malloc(sizeof(double));
    s->fluid->red_integ_diag = gkyl_cu_malloc(sizeof(double[9]));
  }
  // allocate dynamic-vector to store all-reduced integrated moments and F^2
  s->kinetic->integ_L2_f = gkyl_dynvec_new(GKYL_DOUBLE, 1);
  s->fluid->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, 9);
  s->fluid->pkpm.pos_stat = gkyl_dynvec_new(GKYL_DOUBLE, PKPM_POS_NUM_STAT);
  s->kinetic->pkpm.bflux_diag = gkyl_dynvec_new(GKYL_DOUBLE, 2 * cdim);
  s->kinetic->pkpm.is_first_bflux_write_call = true;
  for (int i = 0; i < 2 * GKYL_MAX_CDIM; ++i) {
    s->kinetic->pkpm.bflux[i] = 0.0;
  }
  s->fluid->pkpm.is_first_pos_stat_write_call = true;
  for (int c = 0; c < PKPM_POS_NUM_STAT; ++c) {
    s->fluid->pkpm.pos_count_since[c] = 0;
  }
  s->kinetic->is_first_integ_L2_write_call = true;
  s->fluid->is_first_integ_write_call = true;

  // Initialize applied acceleration for use in force update.
  s->fluid->app_accel = mkarr(app->use_gpu, 3 * app->basis.num_basis, app->local_ext.volume);
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
      s->fluid->app_accel_host = mkarr(false, 3 * app->basis.num_basis, app->local_ext.volume);
    }
    s->fluid->app_accel_proj = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
      .grid = &app->grid,
      .basis = &app->basis,
      .num_quad = app->basis.poly_order + 1,
      .num_ret_vals = 3,
      .eval = s->fluid->info.app_accel,
      .ctx = s->fluid->info.app_accel_ctx,
      .c2p_func = vm_pos_c2p,
      .c2p_func_ctx = app->pos_map,
    });
  }

  // determine collision type to use in PKPM update
  s->kinetic->pkpm.collision_id = s->kinetic->info.collisions.collision_id;
  // not available: collision diagnostics and BGK relaxation to a fixed temperature (source_bgk)
  assert(!s->kinetic->info.collisions.write_coll_diagnostics);
  assert(!s->kinetic->info.collisions.fixed_temp_relax);
  if (s->kinetic->pkpm.collision_id == GKYL_BGK_COLLISIONS) {
    vm_pkpm_species_bgk_init(app, s, &s->kinetic->pkpm.bgk);
  }
  if (s->kinetic->pkpm.collision_id == GKYL_LBO_COLLISIONS) {
    vm_pkpm_species_lbo_init(app, s, &s->kinetic->pkpm.lbo);
  }
  vm_pkpm_species_source_bgk_init(app, s, &s->kinetic->pkpm.src_bgk);

  // An implicit kinetic relaxation that carries a temperature change (the shape factor of
  // implicit BGK with partners, or of the implicit temperature source) needs the weighted stage
  // combination.
  bool has_implicit_cross = s->kinetic->pkpm.collision_id == GKYL_BGK_COLLISIONS &&
                            s->kinetic->info.collisions.is_implicit &&
                            s->kinetic->info.collisions.num_cross_collisions > 0;
  s->kinetic->pkpm.has_stage_weights = !app->use_gpu &&
                                       (has_implicit_cross || s->kinetic->pkpm.src_bgk.is_implicit);
  if (s->kinetic->pkpm.has_stage_weights) {
    s->kinetic->pkpm.stage_lam = mkarr(false, app->basis.num_basis, app->local_ext.volume);
    gkyl_array_clear(s->kinetic->pkpm.stage_lam, 1.0);
    s->kinetic->pkpm.up_stage_weights =
      gkyl_pkpm_stage_weights_inew(&(struct gkyl_pkpm_stage_weights_inp){
        .conf_basis = &app->basis,
        .phase_basis = &s->kinetic->basis,
        .use_gpu = false,
      });
  }

  // The LBO acts at the configuration-space nodes on a nodal copy of the distributions. The BGK
  // operator and the temperature source are weak products and need none.
  s->kinetic->pkpm.has_nodal_ops = !app->use_gpu &&
                                   s->kinetic->pkpm.collision_id == GKYL_LBO_COLLISIONS;
  if (s->kinetic->pkpm.has_nodal_ops) {
    s->kinetic->pkpm.f_nodal = mkarr(false, s->kinetic->f->ncomp, s->kinetic->local_ext.volume);
    s->kinetic->pkpm.rhs_nodal = mkarr(false, s->kinetic->f->ncomp, s->kinetic->local_ext.volume);
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
      .conf_basis = &app->basis,
      .pos_map = app->pos_map,
      .order = order,
      .use_gpu = app->use_gpu,
    };
    for (int d = 0; d < cdim; ++d) {
      diff_inp.D[d] = s->fluid->info.diffusion.D;
    }
    s->fluid->pkpm.mom_diff = gkyl_dg_pkpm_mom_diffusion_inew(&diff_inp);
  }

  // determine which directions are not periodic and get non-periodic BC info
  // also create local skin-ghost ranges for ease of handling non-periodic BCs
  int num_periodic_dir = app->num_periodic_dir, is_np[3] = {1, 1, 1};
  for (int d = 0; d < num_periodic_dir; ++d) {
    is_np[app->periodic_dirs[d]] = 0;
  }

  for (int dir = 0; dir < cdim; ++dir) {
    s->kinetic->lower_bc[dir].type = s->kinetic->upper_bc[dir].type = GKYL_SPECIES_COPY;
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

      s->kinetic->lower_bc[dir].type = bc->lower.type;
      s->kinetic->upper_bc[dir].type = bc->upper.type;
      s->fluid->lower_bc[dir] = bc_fluid[0];
      s->fluid->upper_bc[dir] = bc_fluid[1];
    }
    // Create local lower skin and ghost ranges for distribution function
    gkyl_skin_ghost_ranges(
      &s->kinetic->lower_skin[dir], &s->kinetic->lower_ghost[dir], dir, GKYL_LOWER_EDGE,
      &s->kinetic->pkpm.local_ext_conf, ghost
    );
    // Create local upper skin and ghost ranges for distribution function
    gkyl_skin_ghost_ranges(
      &s->kinetic->upper_skin[dir], &s->kinetic->upper_ghost[dir], dir, GKYL_UPPER_EDGE,
      &s->kinetic->pkpm.local_ext_conf, ghost
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
  s->fluid->bc_buffer = mkarr(app->use_gpu, 5 * app->basis.num_basis, buff_sz_fluid);
  // buffer arrays for fixed function boundary conditions on the distribution functions and the
  // fluid state, for the edges that have one
  for (int d = 0; d < 3; ++d) {
    s->kinetic->bc_buffer_lo_fixed[d] = s->kinetic->bc_buffer_up_fixed[d] = 0;
    s->fluid->bc_buffer_lo_fixed[d] = s->fluid->bc_buffer_up_fixed[d] = 0;
  }
  for (int d = 0; d < cdim; ++d) {
    int nc_fluid = 5 * app->basis.num_basis;
    if (s->kinetic->lower_bc[d].type == GKYL_SPECIES_FIXED_FUNC) {
      s->kinetic->bc_buffer_lo_fixed[d] =
        mkarr(app->use_gpu, s->kinetic->f->ncomp, s->kinetic->lower_skin[d].volume);
    }
    if (s->kinetic->upper_bc[d].type == GKYL_SPECIES_FIXED_FUNC) {
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

  for (int d = 0; d < cdim; ++d) {
    bool is_periodic = false;
    for (int p = 0; p < app->num_periodic_dir; ++p) {
      is_periodic = is_periodic || app->periodic_dirs[p] == d;
    }
    if (!is_periodic && (!vm_pkpm_species_bc_is_supported(s->kinetic->lower_bc[d].type) ||
                         !vm_pkpm_species_bc_is_supported(s->kinetic->upper_bc[d].type) ||
                         !vm_pkpm_species_bc_is_supported(s->fluid->lower_bc[d]) ||
                         !vm_pkpm_species_bc_is_supported(s->fluid->upper_bc[d]))) {
      gkyl_exit(
        "vm_pkpm_species: a PKPM species supports copy, absorbing, reflecting and fixed-function boundary conditions."
      );
    }
  }

  // Absorbing boundaries: empty ghost distributions and a mirrored fluid ghost (see
  // gkyl_vlasov.h); the wall flux is the kinetic flux against the empty ghost.
  s->kinetic->pkpm.bc_is_absorb = false;
  for (int d = 0; d < cdim; ++d) {
    s->kinetic->pkpm.bc_lo_reflect[d] = 0;
    s->kinetic->pkpm.bc_up_reflect[d] = 0;
    s->kinetic->pkpm.bc_lo_theta[d] = s->kinetic->pkpm.bc_up_theta[d] = 0;
    // Lower BC updaters. Each block has its own type; copy BCs by default.
    enum gkyl_bc_basic_type bctype_dist =
      vm_pkpm_species_bc_basic_type_dist(s->kinetic->lower_bc[d].type);
    enum gkyl_bc_basic_type bctype_fluid =
      vm_pkpm_species_bc_basic_type_fluid(s->fluid->lower_bc[d]);
    if (s->kinetic->lower_bc[d].type == GKYL_SPECIES_ABSORB) {
      s->kinetic->pkpm.bc_is_absorb = true;
      s->kinetic->pkpm.bc_lo_theta[d] = gkyl_bc_basic_new(
        d, GKYL_LOWER_EDGE, GKYL_BC_REFLECT, app->basis_on_dev, &app->lower_skin[d],
        &app->lower_ghost[d], app->basis.num_basis, app->cdim, app->use_gpu
      );
    } else if (s->kinetic->lower_bc[d].type == GKYL_SPECIES_REFLECT) {
      s->kinetic->pkpm.bc_lo_reflect[d] = vm_pkpm_species_reflect_new(app, s, d, GKYL_LOWER_EDGE);
    }

    // Distribution function non-periodic lower boundary conditions
    s->kinetic->bc_lo[d] = gkyl_bc_basic_new(
      d, GKYL_LOWER_EDGE, bctype_dist, s->kinetic->basis_on_dev, &s->kinetic->lower_skin[d],
      &s->kinetic->lower_ghost[d], s->kinetic->f->ncomp, app->cdim, app->use_gpu
    );
    // Fluid state non-periodic lower boundary conditions
    s->fluid->bc_lo[d] = gkyl_bc_basic_new(
      d, GKYL_LOWER_EDGE, bctype_fluid, app->basis_on_dev, &app->lower_skin[d],
      &app->lower_ghost[d], s->fluid->fluid->ncomp, app->cdim, app->use_gpu
    );

    // Upper BC updaters.
    bctype_dist = vm_pkpm_species_bc_basic_type_dist(s->kinetic->upper_bc[d].type);
    bctype_fluid = vm_pkpm_species_bc_basic_type_fluid(s->fluid->upper_bc[d]);
    if (s->kinetic->upper_bc[d].type == GKYL_SPECIES_ABSORB) {
      s->kinetic->pkpm.bc_is_absorb = true;
      s->kinetic->pkpm.bc_up_theta[d] = gkyl_bc_basic_new(
        d, GKYL_UPPER_EDGE, GKYL_BC_REFLECT, app->basis_on_dev, &app->upper_skin[d],
        &app->upper_ghost[d], app->basis.num_basis, app->cdim, app->use_gpu
      );
    } else if (s->kinetic->upper_bc[d].type == GKYL_SPECIES_REFLECT) {
      s->kinetic->pkpm.bc_up_reflect[d] = vm_pkpm_species_reflect_new(app, s, d, GKYL_UPPER_EDGE);
    }

    // Distribution function non-periodic upper boundary conditions
    s->kinetic->bc_up[d] = gkyl_bc_basic_new(
      d, GKYL_UPPER_EDGE, bctype_dist, s->kinetic->basis_on_dev, &s->kinetic->upper_skin[d],
      &s->kinetic->upper_ghost[d], s->kinetic->f->ncomp, app->cdim, app->use_gpu
    );
    // Fluid state non-periodic upper boundary conditions
    s->fluid->bc_up[d] = gkyl_bc_basic_new(
      d, GKYL_UPPER_EDGE, bctype_fluid, app->basis_on_dev, &app->upper_skin[d],
      &app->upper_ghost[d], s->fluid->fluid->ncomp, app->cdim, app->use_gpu
    );
  }

  // Methods of the two aspects: initial conditions, the stages of a time step, the L2 norm of the
  // distributions and the read of a file.
  s->kinetic->apply_ic_func = vm_pkpm_kinetic_species_apply_ic;
  s->kinetic->calc_L2_func = vm_pkpm_kinetic_species_calc_L2;
  s->kinetic->write_L2_func = vm_pkpm_kinetic_species_write_L2;
  s->kinetic->from_file_func = vm_pkpm_kinetic_species_from_file;
  s->fluid->apply_ic_func = vm_pkpm_fluid_species_apply_ic;
  s->fluid->from_file_func = vm_pkpm_fluid_species_from_file;
  s->kinetic->rhs_func = vm_pkpm_kinetic_species_rhs;
  s->kinetic->calc_source_moms_func = vm_pkpm_kinetic_species_no_source_moms;
  s->kinetic->rhs_source_func = vm_pkpm_kinetic_species_no_rhs_source;
  s->kinetic->step_f_func = vm_pkpm_species_array_step_f;
  s->kinetic->combine_func = vm_pkpm_species_array_combine;
  s->kinetic->copy_func = vm_pkpm_species_array_copy;
  s->kinetic->bc_func = vm_pkpm_kinetic_species_apply_bc;

  s->fluid->rhs_func = vm_pkpm_fluid_species_rhs;
  s->fluid->rhs_source_func = vm_pkpm_fluid_species_no_rhs_source;
  s->fluid->step_f_func = vm_pkpm_species_array_step_f;
  s->fluid->combine_func = vm_pkpm_species_array_combine;
  s->fluid->copy_func = vm_pkpm_species_array_copy;
  s->fluid->bc_func = vm_pkpm_fluid_species_apply_bc;
  s->fluid->limiter_func = vm_pkpm_fluid_species_no_limiter;
  s->fluid->positivity_func = vm_pkpm_species_apply_positivity;
  // the fluid right-hand side applies its forces with the explicit coupling only
  s->fluid->has_explicit_force = app->use_explicit_sources;
}

// Initial distributions of a PKPM species. The user provides the normalized distributions [F, G,
// M0_perp, N_perp] on the normalized velocity grid (w = v_par/xi; F a probability density in w, G
// its perpendicular-energy partner so that M0(F) = 1, M1(F) = 0 and M2(F) + 2 M0(G) = 3; the
// first-harmonic blocks are zero for a gyrotropic start), sampled at the physical coordinates of
// the position map (identity: unchanged).
static void
vm_pkpm_kinetic_species_apply_ic(gkyl_vlasov_app *app, struct vm_species *kin, double t0)
{
  gkyl_proj_on_basis *proj_dist = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
    .grid = &kin->grid,
    .basis = &kin->basis,
    .num_quad = 8,
    .num_ret_vals = PKPM_DIST_NUM_BLOCKS,
    .eval = kin->info.projection[0].func,
    .ctx = kin->info.projection[0].ctx_func,
    .c2p_func = vm_pos_c2p,
    .c2p_func_ctx = app->pos_map,
  });

  // run updater; need to project onto extended range for ease of handling
  // subsequent operations over extended range such as primitive variable computations
  // This is needed to fill the corner cells as the corner cells may not be filled by
  // boundary conditions and we cannot divide by 0 anywhere or the weak divisions will fail
  // (the velocity ghost cell stays empty)
  gkyl_array_clear(kin->f_host, 0.0);
  gkyl_proj_on_basis_advance(proj_dist, t0, &kin->pkpm.local_ext_conf, kin->f_host);
  gkyl_proj_on_basis_release(proj_dist);

  // note: f_host is the same as f when not on GPUs
  if (app->use_gpu) {
    gkyl_array_copy(kin->f, kin->f_host);
  }

  // Constraints M0(F) = 1 and M2(F) + 2 M0(G) = 3 exactly: nodal rescaling of the
  // projected distributions by their own moments (removes the projection error).
  vm_pkpm_species_moment_calc(&kin->pkpm.nmoms, kin->pkpm.local_ext_conf, app->local_ext, kin->f);
  gkyl_array_copy(kin->f1, kin->f);
  gkyl_dg_calc_pkpm_vars_normalize_dist(
    kin->pkpm.calc_pkpm_vars, &app->local_ext, &kin->local_ext, kin->pkpm.nmoms.marr, kin->f1,
    kin->f
  );
  // The evolved distributions carry the position-map Jacobian, J [F, G] (identity map: no
  // change). The velocity ghost row of f stays empty.
  gkyl_vlasov_position_map_rescale_jacobpos(
    app->pos_map, &kin->basis, &kin->pkpm.local_ext_conf, kin->f, kin->f
  );
  if (app->use_gpu) {
    gkyl_array_copy(kin->f_host, kin->f);
  }

  // copy the initial conditions into the buffers of the fixed-function boundary conditions
  for (int d = 0; d < app->cdim; ++d) {
    if (kin->lower_bc[d].type == GKYL_SPECIES_FIXED_FUNC) {
      gkyl_bc_basic_buffer_fixed_func(kin->bc_lo[d], kin->bc_buffer_lo_fixed[d], kin->f);
    }
    if (kin->upper_bc[d].type == GKYL_SPECIES_FIXED_FUNC) {
      gkyl_bc_basic_buffer_fixed_func(kin->bc_up[d], kin->bc_buffer_up_fixed[d], kin->f);
    }
  }
}

// Initial fluid state [rho, rho ux, rho uy, rho uz, E] of a PKPM species, sampled at the physical
// coordinates of the position map, with the applied acceleration of the initial time.
static void
vm_pkpm_fluid_species_project_ic(gkyl_vlasov_app *app, struct vm_fluid_species *fl, double t0)
{
  gkyl_proj_on_basis *proj_fluid = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
    .grid = &app->grid,
    .basis = &app->basis,
    .num_quad = app->poly_order + 1,
    .num_ret_vals = 5,
    .eval = fl->info.init,
    .ctx = fl->info.ctx,
    .c2p_func = vm_pos_c2p,
    .c2p_func_ctx = app->pos_map,
  });

  // projected onto the extended range, like the distributions
  gkyl_proj_on_basis_advance(proj_fluid, t0, &app->local_ext, fl->fluid_host);
  gkyl_proj_on_basis_release(proj_fluid);

  // note: fluid_host is the same as fluid when not on GPUs
  if (app->use_gpu) {
    gkyl_array_copy(fl->fluid, fl->fluid_host);
  }
  // The evolved state carries the position-map Jacobian, J [rho, M, E] (identity map: no change).
  gkyl_vlasov_position_map_rescale_jacobpos_conf(
    app->pos_map, &app->local_ext, fl->fluid, fl->fluid
  );
  if (app->use_gpu) {
    gkyl_array_copy(fl->fluid_host, fl->fluid);
  }

  // pre-compute the applied acceleration (it is not recomputed when time-independent)
  vm_fluid_species_calc_app_accel(app, fl, t0);

  // copy the initial conditions into the buffers of the fixed-function boundary conditions
  for (int d = 0; d < app->cdim; ++d) {
    if (fl->lower_bc[d] == GKYL_SPECIES_FIXED_FUNC) {
      gkyl_bc_basic_buffer_fixed_func(fl->bc_lo[d], fl->bc_buffer_lo_fixed[d], fl->fluid);
    }
    if (fl->upper_bc[d] == GKYL_SPECIES_FIXED_FUNC) {
      gkyl_bc_basic_buffer_fixed_func(fl->bc_up[d], fl->bc_buffer_up_fixed[d], fl->fluid);
    }
  }
}

// Initial conditions of the fluid aspect: the projected state, limited.
static void
vm_pkpm_fluid_species_apply_ic(gkyl_vlasov_app *app, struct vm_fluid_species *fl, double t0)
{
  vm_pkpm_fluid_species_project_ic(app, fl, t0);
  vm_pkpm_species_apply_positivity(app, fl, fl->fluid);
}

// Kinetic coupling moments [p_par, p_perp, q_par, q_perp, s, r] = [p M2, p M0_G, p xi M3, p xi M1_G,
// p m1, p xi (m2/2 + n0)] from the raw normalized moments, the nodal scalar pressure p = rho xi^2
// and the thermal speed xi (weak products projected onto p=1): the pressures, first-harmonic stress
// and heat fluxes the fluid system needs. The six constraints [M0(F), M1(F),
// (M2(F) + 2 M0(G))/3, int M0_perp dw] are copied out for the diagnostics.
static void
vm_pkpm_species_coupling_moms(gkyl_vlasov_app *app, struct vlasov_species *species)
{
  int nc = app->basis.num_basis;
  // layout of the normalized moments: [M0, M2, M0_G, M1, M3, M1_G, mu (3), m1 (3), m2 (3), n0 (3)]
  static const int raw[7] = {1, 2, 4, 5, 9, 10, 11}; // M2, M0_G, M3, M1_G, m1 (3)
  for (int c = 0; c < 10; ++c) {
    if (c < 7) {
      gkyl_dg_mul_op_range(
        &app->basis, c, species->fluid->pkpm.moms, 2, species->fluid->pkpm.xi, raw[c],
        species->kinetic->pkpm.nmoms.marr, &app->local_ext
      );
    } else {
      int a = c - 7;
      gkyl_array_set_offset(
        species->fluid->pkpm.scratch, 0.5, species->kinetic->pkpm.nmoms.marr, (12 + a) * nc
      );
      gkyl_array_accumulate_offset(
        species->fluid->pkpm.scratch, 1.0, species->kinetic->pkpm.nmoms.marr, (15 + a) * nc
      );
      gkyl_dg_mul_op_range(
        &app->basis, c, species->fluid->pkpm.moms, 2, species->fluid->pkpm.xi, 0,
        species->fluid->pkpm.scratch, &app->local_ext
      );
    }
    if (c == 2 || c == 3 || c >= 7) {
      gkyl_dg_mul_op_range(
        &app->basis, 0, species->fluid->pkpm.scratch, c, species->fluid->pkpm.moms, 0,
        species->fluid->pkpm.xi, &app->local_ext
      );
      gkyl_array_set_offset(species->fluid->pkpm.moms, 1.0, species->fluid->pkpm.scratch, c * nc);
    }
  }
  gkyl_dg_calc_pkpm_vars_constraints(
    species->kinetic->pkpm.calc_pkpm_vars, &app->local_ext, nc, species->kinetic->pkpm.nmoms.marr,
    species->kinetic->pkpm.cmoms
  );
  // the constraints of the physical distributions ([1, 0, 1]): the moments of J [F, G] carry J
  gkyl_vlasov_position_map_divide_jacobpos_conf(
    app->pos_map, &app->local_ext, species->kinetic->pkpm.cmoms->ncomp,
    species->kinetic->pkpm.cmoms, species->kinetic->pkpm.cmoms
  );
}

void
vm_pkpm_species_calc_pkpm_vars(
  gkyl_vlasov_app *app, struct vlasov_species *species, const struct gkyl_array *fin,
  const struct gkyl_array *fluidin
)
{
  struct timespec tm = gkyl_wall_clock();
  int nc = app->basis.num_basis;

  // Raw normalized moments of the blocks over the extended range, the primitive variables of
  // the physical fluid state (position-map Jacobian divided out), then the coupling moments,
  // which carry J through the moments of J [F, G]. theta = T_perp/T is M0(G).
  gkyl_vlasov_position_map_divide_jacobpos_conf(
    app->pos_map, &app->local_ext, fluidin->ncomp, fluidin, species->fluid->pkpm.fluid_no_J
  );
  vm_pkpm_species_moment_calc(
    &species->kinetic->pkpm.nmoms, species->kinetic->pkpm.local_ext_conf, app->local_ext, fin
  );
  gkyl_array_set_offset(species->fluid->pkpm.rho, 1.0, species->fluid->pkpm.fluid_no_J, 0);
  // u = M/rho and [xi, xi^2, p] by nodal quotients (the weak division at p=1).
  gkyl_dg_calc_pkpm_vars_prim(
    species->fluid->pkpm.calc_pkpm_vars, &app->local_ext, species->fluid->pkpm.fluid_no_J,
    species->fluid->pkpm.u, species->fluid->pkpm.xi
  );
  // the kinetic aspect reads its own copies of the flow and the thermal speed
  gkyl_array_copy(species->kinetic->pkpm.u, species->fluid->pkpm.u);
  gkyl_array_copy(species->kinetic->pkpm.xi, species->fluid->pkpm.xi);
  vm_pkpm_species_coupling_moms(app, species);
  // theta = M0(G) of the physical G: M0(J G)/J.
  gkyl_array_set_offset(
    species->kinetic->pkpm.theta, 1.0, species->kinetic->pkpm.nmoms.marr, 2 * nc
  );
  gkyl_vlasov_position_map_divide_jacobpos_conf(
    app->pos_map, &app->local_ext, nc, species->kinetic->pkpm.theta, species->kinetic->pkpm.theta
  );
  // The face fluxes average u, xi and theta across every face. At absorbing walls u and xi of
  // the mirrored fluid ghost are the mirror image of the skin values already; theta, the
  // zeroth moment of the empty ghost G, is replaced by the mirror image of the skin cell, so
  // that the face averages are the interior wall traces.
  if (species->kinetic->pkpm.bc_is_absorb) {
    for (int d = 0; d < app->cdim; ++d) {
      if (species->kinetic->pkpm.bc_lo_theta[d]) {
        gkyl_bc_basic_advance(
          species->kinetic->pkpm.bc_lo_theta[d], species->fluid->bc_buffer,
          species->kinetic->pkpm.theta
        );
      }
      if (species->kinetic->pkpm.bc_up_theta[d]) {
        gkyl_bc_basic_advance(
          species->kinetic->pkpm.bc_up_theta[d], species->fluid->bc_buffer,
          species->kinetic->pkpm.theta
        );
      }
    }
    gkyl_comm_array_sync(app->comm, &app->local, &app->local_ext, species->kinetic->pkpm.theta);
  }

  app->stat.species_pkpm_vars_tm += gkyl_time_diff_now_sec(tm);
}

void
vm_pkpm_species_positivity_record_save(struct vm_fluid_species *fl)
{
  for (int c = 0; c < PKPM_POS_NUM_STAT; ++c) {
    fl->pkpm.pos_count_since_saved[c] = fl->pkpm.pos_count_since[c];
  }
  gkyl_array_copy(fl->pkpm.pos_count_saved, fl->pkpm.pos_count);
}

void
vm_pkpm_species_positivity_record_restore(struct vm_fluid_species *fl)
{
  for (int c = 0; c < PKPM_POS_NUM_STAT; ++c) {
    fl->pkpm.pos_count_since[c] = fl->pkpm.pos_count_since_saved[c];
  }
  gkyl_array_copy(fl->pkpm.pos_count, fl->pkpm.pos_count_saved);
}

void
vm_pkpm_species_apply_positivity(
  gkyl_vlasov_app *app, struct vm_fluid_species *fl, struct gkyl_array *fluid
)
{
  struct timespec tm = gkyl_wall_clock();
  int num_rho = 0, num_e = 0, num_avg = 0;
  gkyl_dg_calc_pkpm_vars_positivity(
    fl->pkpm.calc_pkpm_vars, &app->local, fl->positivity_eps, fluid, fl->pkpm.pos_count, &num_rho,
    &num_e, &num_avg
  );
  fl->pkpm.pos_count_since[0] += num_rho;
  fl->pkpm.pos_count_since[1] += num_e;
  fl->pkpm.pos_count_since[2] += num_avg;
  fl->pkpm.pos_count_since[3] += 1;
  app->stat.species_pkpm_vars_tm += gkyl_time_diff_now_sec(tm);
}

// Configuration-space face fluxes of the stage (the kinetic flux each side of a face receives,
// the fluid flux as moments of the same face flux, and both rates), then the in-cell gradients
// of u, rho and xi with the blend-weighted lift of the face jumps.
void
vm_pkpm_species_face_flux(
  gkyl_vlasov_app *app, struct vlasov_species *species, const struct gkyl_array *fin
)
{
  gkyl_dg_pkpm_remap_flux_surf_advance(
    species->kinetic->pkpm.calc_face_flux, &app->local, &app->local_ext, &species->kinetic->local,
    &species->kinetic->local_ext, app->field->bvar, species->fluid->pkpm.fluid_no_J,
    species->fluid->pkpm.u, species->fluid->pkpm.xi, species->kinetic->pkpm.theta, fin,
    species->kinetic->pkpm.conf_flux_surf, species->fluid->pkpm.fluid_flux_surf,
    species->kinetic->cflrate, species->fluid->cflrate, species->kinetic->pkpm.flux_blend
  );
  // The gradient of u for the kinetic sources (bb:grad u); needs u in the ghost cells (at
  // absorbing boundaries that of the mirrored fluid state).
  gkyl_dg_calc_pkpm_vars_in_cell_grad(
    species->fluid->pkpm.calc_pkpm_vars, &app->local, species->kinetic->pkpm.flux_blend,
    species->fluid->pkpm.u, species->kinetic->pkpm.grad_u
  );
  // grad(ln n) for the density-gradient source: the gradient of [rho, 0, 0] divided by rho.
  gkyl_array_set_offset(species->fluid->pkpm.rho3, 1.0, species->fluid->pkpm.rho, 0);
  gkyl_dg_calc_pkpm_vars_in_cell_grad(
    species->fluid->pkpm.calc_pkpm_vars, &app->local, species->kinetic->pkpm.flux_blend,
    species->fluid->pkpm.rho3, species->fluid->pkpm.grad_rho
  );
  gkyl_dg_calc_pkpm_vars_grad_ln_n(
    species->fluid->pkpm.calc_pkpm_vars, &app->local, species->fluid->pkpm.rho,
    species->fluid->pkpm.grad_rho, species->kinetic->pkpm.grad_ln_n
  );
  // The gradient of xi for the forces -w^2 b.grad xi and -2 w b.grad xi G (computed for all of
  // [xi, xi^2, p]; only the xi component is read).
  gkyl_array_set_offset(species->fluid->pkpm.xi3, 1.0, species->fluid->pkpm.xi, 0);
  gkyl_dg_calc_pkpm_vars_in_cell_grad(
    species->fluid->pkpm.calc_pkpm_vars, &app->local, species->kinetic->pkpm.flux_blend,
    species->fluid->pkpm.xi3, species->kinetic->pkpm.grad_xi
  );
}

// Collisionless configuration-space pass of the kinetic update into rhs (cleared here): the
// streaming terms with the stored face fluxes, then the first-harmonic w fluxes and sources.
// The closures are computed from the moments of this partial right-hand side, after the
// sources and collision operators have been accumulated into it.
static void
vm_pkpm_species_kinetic_conf_pass(
  gkyl_vlasov_app *app, struct vm_species *kin, const struct gkyl_array *fin, struct gkyl_array *rhs
)
{
  gkyl_array_clear(rhs, 0.0);
  gkyl_hyper_dg_advance(kin->pkpm.slvr_kinetic_conf, &kin->local, fin, kin->cflrate, rhs);
  // the known first-harmonic w fluxes and sources of the eight blocks (modal, exact)
  gkyl_dg_pkpm_harm_vol_coef(
    kin->pkpm.harm_vol, &app->local, app->field->bvar, app->field->grad_b, kin->pkpm.u,
    kin->pkpm.grad_u, kin->pkpm.xi, kin->pkpm.grad_xi, kin->pkpm.grad_ln_n, kin->pkpm.theta,
    app->field->dtB_over_magB, app->field->tot_em, kin->pkpm.harm_coef, kin->pkpm.harm_wnode
  );
  gkyl_dg_pkpm_harm_vol_advance(
    kin->pkpm.harm_vol, &app->local, &kin->local, kin->pkpm.harm_coef, fin, kin->cflrate, rhs
  );
}

void
vm_pkpm_species_accumulate_conf_phase(
  gkyl_vlasov_app *app, struct vm_species *kin, struct gkyl_array *out, double a,
  const struct gkyl_array *conf, const struct gkyl_array *phase
)
{
  int np = kin->basis.num_basis;
  for (int b = 0; b < phase->ncomp / np; ++b) {
    gkyl_array_set_offset(kin->pkpm.blk_in, 1.0, phase, b * np);
    gkyl_array_clear(kin->pkpm.blk_out, 0.0);
    gkyl_dg_mul_conf_phase_op_accumulate_range(
      &app->basis, &kin->basis, kin->pkpm.blk_out, a, conf, kin->pkpm.blk_in, &app->local,
      &kin->local
    );
    gkyl_array_accumulate_offset_range(out, 1.0, kin->pkpm.blk_out, b * np, &kin->local);
  }
}

void
vm_pkpm_species_project_transverse(
  gkyl_vlasov_app *app, struct vlasov_species *species, struct gkyl_array *f
)
{
  // Extended ranges: b is formed in the ghost cells too, and the projection commutes with the
  // periodic copy and the reflecting-wall image, so ghost and skin stay consistent.
  gkyl_dg_calc_pkpm_vars_transverse(
    species->kinetic->pkpm.calc_pkpm_vars, &app->local_ext, &species->kinetic->local_ext,
    app->field->bvar, f
  );
}

void
vm_pkpm_species_harm_exp_advance(
  gkyl_vlasov_app *app, struct vlasov_species *species, int stage, double h, struct gkyl_array *fout
)
{
  // The mixing is frozen with the field and the state of stage 1 (the start of the step); a
  // restarted step freezes it again with the new step.
  if (stage == 1) {
    gkyl_dg_pkpm_harm_exp_coef(
      species->kinetic->pkpm.harm_exp, &app->local, app->field->tot_em,
      species->kinetic->pkpm.harm_exp_coef
    );
    gkyl_array_clear(species->kinetic->pkpm.harm_exp_act, 0.0);
  }
  // With stage weights the combination of the driver is weighted (vlasov_update_implicit.c), and the
  // stage value is built for those weights.
  gkyl_dg_pkpm_harm_exp_advance(
    species->kinetic->pkpm.harm_exp, stage, h, &app->local, &species->kinetic->local,
    species->kinetic->pkpm.harm_exp_coef, species->kinetic->f,
    species->kinetic->pkpm.has_stage_weights ? species->kinetic->pkpm.stage_lam : 0,
    species->kinetic->pkpm.harm_exp_nu, species->kinetic->pkpm.harm_exp_act, fout
  );
  gkyl_dg_pkpm_harm_exp_restore(
    species->kinetic->pkpm.harm_exp, &app->local, &species->kinetic->local,
    species->kinetic->pkpm.harm_exp_act, species->kinetic->pkpm.harm_exp_col, fout
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
vm_pkpm_species_kinetic_closure(
  gkyl_vlasov_app *app, struct vm_species *kin, const struct gkyl_array *fin,
  const struct gkyl_array *rhs, bool from_step
)
{
  if (from_step) {
    gkyl_array_set_offset(kin->pkpm.alpha_0, 1.0, kin->pkpm.closure_step, 0);
    gkyl_array_set_offset(kin->pkpm.c_h, 1.0, kin->pkpm.closure_step, app->basis.num_basis);
  }
  gkyl_array_clear(kin->pkpm.rhs_moms, 0.0);
  vm_pkpm_species_moment_calc_into(
    &kin->pkpm.nmoms, kin->local, app->local, rhs, kin->pkpm.rhs_moms
  );
  return gkyl_dg_calc_pkpm_vars_closure(
    kin->pkpm.calc_pkpm_vars, &app->local, &kin->local, &kin->grid, kin->pkpm.nmoms.marr,
    kin->pkpm.rhs_moms, app->field->bvar, kin->pkpm.harm_wnode, kin->pkpm.theta, fin, kin->pkpm.c_0,
    kin->pkpm.alpha_0, kin->pkpm.c_h, kin->pkpm.a_perp, kin->pkpm.vel_flux_surf, kin->cflrate
  );
}

void
vm_pkpm_species_complete_step(gkyl_vlasov_app *app, struct vlasov_species *sp)
{
  struct vm_species *kin = sp->kinetic;
  gkyl_array_set_offset(kin->pkpm.closure_step, 1.0, kin->pkpm.alpha_0, 0);
  gkyl_array_set_offset(kin->pkpm.closure_step, 1.0, kin->pkpm.c_h, app->basis.num_basis);
}

// Record the nulled closure nodes of a stage in the statistics of the species.
static void
vm_pkpm_species_kinetic_closure_stat(gkyl_vlasov_app *app, struct vm_species *kin, int num_null)
{
  kin->pkpm.num_null_closure += num_null;
  kin->pkpm.closure_cond_min = fmin(
    kin->pkpm.closure_cond_min, gkyl_dg_calc_pkpm_vars_closure_cond_min(kin->pkpm.calc_pkpm_vars)
  );
  // A nulled node means the traces at the velocity-domain edges are large enough to make the
  // closure system singular: the constraints drift in that cell. Say so once, with the proxy.
  if (num_null > 0 && !kin->pkpm.closure_null_warned) {
    kin->pkpm.closure_null_warned = true;
    int rank = 0;
    gkyl_comm_get_rank(app->comm, &rank);
    fprintf(
      stderr,
      "WARNING (rank %d): species %s at t = %g: the closure solve was nulled at %d nodes (min det/det0 so far %g).\n"
      "  The normalized distribution is not negligible at the velocity-space boundary; the constraints drift there.\n"
      "  The velocity extents are likely too small for this species.\n",
      rank, kin->name, app->tcurr, num_null, kin->pkpm.closure_cond_min
    );
  }
}

void
vm_pkpm_species_calc_diag_vars(
  gkyl_vlasov_app *app, struct vlasov_species *species, const struct gkyl_array *fin,
  const struct gkyl_array *fluidin
)
{
  vm_pkpm_species_calc_pkpm_vars(app, species, fin, fluidin);
  gkyl_array_clear(species->kinetic->cflrate, 0.0);
  gkyl_array_clear(species->fluid->cflrate, 0.0);
  vm_pkpm_species_face_flux(app, species, fin);
  // fnew is free scratch between time steps.
  vm_pkpm_species_kinetic_conf_pass(app, species->kinetic, fin, species->kinetic->fnew);
  vm_pkpm_species_kinetic_closure(app, species->kinetic, fin, species->kinetic->fnew, true);
}

// Right-hand side of the kinetic equations of a species. The face fluxes of the stage, the source
// distributions and the rates were staged for every species before any right-hand side. The
// species reads the physical fields of the stage on the field object (em is not used).
double
vm_pkpm_kinetic_species_rhs(
  gkyl_vlasov_app *app, struct vm_species *kin, const struct gkyl_array *fin,
  const struct gkyl_array *em, struct gkyl_array *rhs
)
{
  struct timespec tm = gkyl_wall_clock();
  // Kinetic equations: the collisionless configuration pass (clears rhs), the sources and
  // collisions, the closures from the moments of that partial right-hand side, then the
  // velocity pass accumulating into the same rhs.
  vm_pkpm_species_kinetic_conf_pass(app, kin, fin, rhs);
  // Sources and collisions enter before the closures, which keep the constraints exact
  // against them (the fluid carries the density and energy changes). The LBO conserves the
  // constraints on its own (boundary corrections, zero-flux w edges).
  vm_pkpm_species_source_kinetic_rhs(app, kin, &kin->pkpm.src, fin, rhs);
  vm_pkpm_species_source_bgk_kinetic_rhs(app, kin, &kin->pkpm.src_bgk, fin, rhs);
  if (kin->pkpm.collision_id == GKYL_BGK_COLLISIONS) {
    vm_pkpm_species_bgk_kinetic_rhs(app, kin, &kin->pkpm.bgk, fin, rhs);
  }
  if (kin->pkpm.collision_id == GKYL_LBO_COLLISIONS) {
    vm_pkpm_species_lbo_rhs(app, kin, &kin->pkpm.lbo, fin, rhs);
  }
  // The nodal increment of the LBO, converted and added to rhs.
  vm_pkpm_species_nodal_accumulate(app, kin, rhs);
  int num_null = vm_pkpm_species_kinetic_closure(app, kin, fin, rhs, kin->pkpm.closure_from_step);
  vm_pkpm_species_kinetic_closure_stat(app, kin, num_null);
  gkyl_dg_pkpm_vel_flux_surf_advance(
    kin->pkpm.calc_vel_flux, &app->local, &kin->local, kin->pkpm.harm_wnode, kin->pkpm.theta,
    app->field->bvar, kin->pkpm.alpha_0, kin->pkpm.c_h, fin, kin->pkpm.vel_flux_surf
  );
  gkyl_hyper_dg_advance(kin->pkpm.slvr_kinetic_vel, &kin->local, fin, kin->cflrate, rhs);
  app->stat.species_rhs_tm += gkyl_time_diff_now_sec(tm);

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

// Right-hand side of the fluid equations of a species. The explicit force takes the physical
// fields of the stage (em, the evolved field, carries the position-map Jacobian).
double
vm_pkpm_fluid_species_rhs(
  gkyl_vlasov_app *app, struct vm_fluid_species *fl, const struct gkyl_array *fluidin,
  const struct gkyl_array *em, struct gkyl_array *rhs
)
{
  gkyl_array_clear(rhs, 0.0);

  // The energy relaxation of the explicit temperature source of the stage (staged with the
  // species), then the volume and surface terms from the stored face fluxes and the source.
  struct timespec tm = gkyl_wall_clock();
  if (fl->pkpm.has_src_bgk_rhs) {
    gkyl_array_accumulate_range(rhs, 1.0, fl->pkpm.src_bgk_rhs, &app->local);
  }
  gkyl_hyper_dg_advance(fl->pkpm.slvr_fluid, &app->local, fluidin, fl->cflrate, rhs);
  if (fl->pkpm.src_scale != 0.0) {
    gkyl_array_accumulate_range(rhs, fl->pkpm.src_scale, fl->pkpm.src_fluid, &app->local);
  }
  // Kinetic collision partners: the momentum and energy take the opposite of what the cross
  // operator of each partner receives from this species.
  for (int k = 0; k < fl->pkpm.num_kin_exch; ++k) {
    gkyl_array_accumulate_offset_range(
      rhs, -1.0, fl->pkpm.kin_exch[k], PKPM_FLUID_MOM * app->basis.num_basis, &app->local
    );
  }
  app->stat.fluid_species_rhs_tm += gkyl_time_diff_now_sec(tm);

  if (fl->has_diffusion) {
    gkyl_dg_pkpm_mom_diffusion_advance(fl->pkpm.mom_diff, &app->local, fluidin, fl->cflrate, rhs);
  }

  // If PKPM update is fully explicit, include the fluid-EM coupling in fluid update
  if (app->use_explicit_sources) {
    double qbym = fl->charge / fl->mass;
    gkyl_array_set(fl->qmem, qbym, app->field->em_no_J);

    // Accumulate applied acceleration and/or q/m*(external electromagnetic)
    // fields onto qmem to get the total acceleration
    if (fl->has_app_accel) {
      gkyl_array_accumulate_range(fl->qmem, 1.0, fl->app_accel, &app->local);
    }
    if (app->field->has_ext_em) {
      gkyl_array_accumulate_range(fl->qmem, qbym, app->field->ext_em, &app->local);
    }

    gkyl_dg_calc_pkpm_vars_source(
      fl->pkpm.calc_pkpm_vars, &app->local, fl->qmem, app->src_omega_p_sq, fluidin, fl->cflrate, rhs
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
vm_pkpm_species_calc_self_moms(
  gkyl_vlasov_app *app, struct vlasov_species *sp, int stage, struct gkyl_array *fin,
  const struct gkyl_array *fluidin
)
{
  // The closures of the first stage of a step take the speeds of the velocity-domain edges from
  // the last stage of the last completed step, those of the other stages from the stage before.
  sp->kinetic->pkpm.closure_from_step = stage == 1;
  // Every kernel of the stage sees first-harmonic blocks perpendicular to the stage's b: project
  // M0_perp and N_perp of the stage state at the Gauss-Legendre nodes of b.
  vm_pkpm_species_project_transverse(app, sp, fin);
  // Configuration-space variables of the stage: the moments of the normalized distributions,
  // the flow and thermal speed, the kinetic coupling moments and T_perp/T.
  vm_pkpm_species_calc_pkpm_vars(app, sp, fin, fluidin);
  // The distributions of the stage at the configuration-space nodes, for the LBO.
  vm_pkpm_species_nodal_dist(app, sp, fin);
  // compute necessary moments and boundary corrections for collisions
  if (sp->kinetic->pkpm.collision_id == GKYL_LBO_COLLISIONS) {
    vm_pkpm_species_lbo_moms(app, sp, &sp->kinetic->pkpm.lbo, fin);
  }
  // node values of the fluid state for the computed frequencies and cross-species collisions
  vm_pkpm_species_coll_cross_prim(app, sp, fluidin);

  // Configuration-space face fluxes of the stage (kinetic, fluid and their rates, which start
  // here) and the particle loss rates they imply through absorbing walls.
  struct timespec tm = gkyl_wall_clock();
  gkyl_array_clear(sp->kinetic->cflrate, 0.0);
  gkyl_array_clear(sp->fluid->cflrate, 0.0);
  vm_pkpm_species_face_flux(app, sp, fin);
  vm_pkpm_species_calc_bflux(app, sp);
  app->stat.fluid_species_rhs_tm += gkyl_time_diff_now_sec(tm);

  // Explicit temperature source: its energy relaxation for the fluid right-hand side and the
  // rate of the distributions for the kinetic one.
  vm_pkpm_species_source_bgk_calc_fluid(app, sp, &sp->kinetic->pkpm.src_bgk, fluidin);
}

void
vm_pkpm_species_calc_coupled_vars(
  gkyl_vlasov_app *app, struct vlasov_species *sp, const struct gkyl_array *fin,
  const struct gkyl_array *fluidin
)
{
  // rate of the kinetic operator with the partner species: needs the node values (and the
  // self-collision frequencies) of every species
  vm_pkpm_species_coll_cross_moms(app, sp);
  // rate of the boundary-flux source (which reads another species' losses) and the source
  // distributions of the stage
  vm_pkpm_species_source_calc_rate(app, sp, &sp->kinetic->pkpm.src);
  vm_pkpm_species_source_calc_dist(app, sp, &sp->kinetic->pkpm.src);
}

void
vm_pkpm_species_finish_stage(
  gkyl_vlasov_app *app, struct vlasov_species *sp, int stage, double dt, struct gkyl_array *fout,
  struct gkyl_array *fluidout
)
{
  // at the nodes the gyrophase mixing damps, the first-harmonic blocks are replaced by the
  // integrating-factor stage; the fluid state is limited
  vm_pkpm_species_harm_exp_advance(app, sp, stage, dt, fout);
  vm_pkpm_species_apply_positivity(app, sp->fluid, fluidout);
}

// Determine which directions are periodic and which directions are not periodic,
// and then apply boundary conditions for distribution function
void
vm_pkpm_species_apply_bc(gkyl_vlasov_app *app, const struct vm_species *kin, struct gkyl_array *f)
{
  struct timespec wst = gkyl_wall_clock();

  int num_periodic_dir = app->num_periodic_dir, cdim = app->cdim;
  gkyl_comm_array_per_sync(
    kin->comm, &kin->local, &kin->pkpm.local_ext_conf, num_periodic_dir, app->periodic_dirs, f
  );

  int is_np_bc[3] = {1, 1, 1}; // flags to indicate if direction is periodic
  for (int d = 0; d < num_periodic_dir; ++d) {
    is_np_bc[app->periodic_dirs[d]] = 0;
  }

  for (int d = 0; d < cdim; ++d) {
    if (is_np_bc[d]) {
      switch (kin->lower_bc[d].type) {
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(kin->bc_lo[d], kin->bc_buffer, f);
          break;
        case GKYL_SPECIES_REFLECT:
          gkyl_bc_pkpm_species_reflect_advance(kin->pkpm.bc_lo_reflect[d], app->field->bvar, f);
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

      switch (kin->upper_bc[d].type) {
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(kin->bc_up[d], kin->bc_buffer, f);
          break;
        case GKYL_SPECIES_REFLECT:
          gkyl_bc_pkpm_species_reflect_advance(kin->pkpm.bc_up_reflect[d], app->field->bvar, f);
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

  gkyl_comm_array_sync(kin->comm, &kin->local, &kin->pkpm.local_ext_conf, f);

  app->stat.species_bc_tm += gkyl_time_diff_now_sec(wst);
}

// Determine which directions are periodic and which directions are not periodic,
// and then apply boundary conditions for the fluid state
void
vm_pkpm_fluid_species_apply_bc(
  gkyl_vlasov_app *app, const struct vm_fluid_species *fl, struct gkyl_array *fluid
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

static void
vm_pkpm_kinetic_species_calc_L2(gkyl_vlasov_app *app, struct vm_species *kin, double tm)
{
  struct timespec wst = gkyl_wall_clock();
  // The evolved f carries the position-map Jacobian: int f^2 dx = int (f/J) f dxi.
  gkyl_vlasov_position_map_divide_jacobpos(
    app->pos_map, &kin->basis, &kin->local, kin->f, kin->f_no_J
  );
  gkyl_dg_calc_prod_op_range(kin->basis, 0, kin->L2_f, 0, kin->f_no_J, 0, kin->f, kin->local);
  gkyl_array_scale_range(kin->L2_f, kin->grid.cellVolume, &kin->local);

  double L2[1] = {0.0};
  if (app->use_gpu) {
    gkyl_array_reduce_range(kin->red_L2_f, kin->L2_f, GKYL_SUM, &kin->local);
    gkyl_cu_memcpy(L2, kin->red_L2_f, sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    gkyl_array_reduce_range(L2, kin->L2_f, GKYL_SUM, &kin->local);
  }
  double L2_global[1] = {0.0};
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, L2, L2_global);

  gkyl_dynvec_append(kin->integ_L2_f, tm, L2_global);

  app->stat.species_diag_calc_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_diag += 1;
}

void
vm_pkpm_species_nodal_dist(
  gkyl_vlasov_app *app, const struct vlasov_species *species, const struct gkyl_array *fin
)
{
  if (species->kinetic->pkpm.has_nodal_ops) {
    struct timespec tm = gkyl_wall_clock();
    gkyl_dg_conf_nodal_phase_to_nodal(
      app->conf_nodal, &species->kinetic->local, fin, species->kinetic->pkpm.f_nodal
    );
    gkyl_array_clear_range(species->kinetic->pkpm.rhs_nodal, 0.0, &species->kinetic->local);
    app->stat.species_coll_tm += gkyl_time_diff_now_sec(tm);
  }
}

void
vm_pkpm_species_nodal_accumulate(
  gkyl_vlasov_app *app, const struct vm_species *kin, struct gkyl_array *out
)
{
  if (kin->pkpm.has_nodal_ops) {
    struct timespec tm = gkyl_wall_clock();
    gkyl_dg_conf_nodal_phase_to_modal_accumulate(
      app->conf_nodal, &kin->local, kin->pkpm.rhs_nodal, out
    );
    app->stat.species_coll_tm += gkyl_time_diff_now_sec(tm);
  }
}
// Write the evolved state of a PKPM species for a frame: the distributions [F, G, M0_perp,
// N_perp] in <name>_<frame>, the fluid state [rho, rho ux, rho uy, rho uz, E] in
// <name>_fluid_<frame> and the closures of the last completed step.
static void
vm_pkpm_species_write(gkyl_vlasov_app *app, struct vlasov_species *s, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  struct gkyl_msgpack_data *mt = vlasov_array_meta_new((struct vlasov_output_meta){
    .frame = frame,
    .stime = tm,
    .poly_order = app->poly_order,
    .basis_type = s->kinetic->basis.id,
    .num_cdim = app->cdim,
    .num_vdim = 1,
  });

  const char *fmt = "%s-%s_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, s->name, frame);
  char fileNm[sz + 1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, s->name, frame);

  // The evolved distributions carry the position-map Jacobian; write the physical ones.
  gkyl_vlasov_position_map_divide_jacobpos(
    app->pos_map, &s->kinetic->basis, &s->kinetic->local, s->kinetic->f, s->kinetic->f_no_J
  );
  if (app->use_gpu) {
    // copy data from device to host before writing it out
    gkyl_array_copy(s->kinetic->pkpm.f_no_J_host, s->kinetic->f_no_J);
  }
  gkyl_comm_array_write(
    s->kinetic->comm, &s->kinetic->grid, &s->kinetic->local, mt, s->kinetic->pkpm.f_no_J_host,
    fileNm
  );

  vlasov_array_meta_release(mt);

  // The fluid state [rho, rho ux, rho uy, rho uz, E] of the species, likewise the physical one.
  struct gkyl_msgpack_data *mt_fluid = vlasov_array_meta_new((struct vlasov_output_meta){
    .frame = frame,
    .stime = tm,
    .poly_order = app->poly_order,
    .basis_type = app->basis.id,
  });

  const char *fmt_fluid = "%s-%s_fluid_%d.gkyl";
  int sz_fluid = gkyl_calc_strlen(fmt_fluid, app->name, s->name, frame);
  char fileNm_fluid[sz_fluid + 1]; // ensures no buffer overflow
  snprintf(fileNm_fluid, sizeof fileNm_fluid, fmt_fluid, app->name, s->name, frame);

  gkyl_vlasov_position_map_divide_jacobpos_conf(
    app->pos_map, &app->local, s->fluid->fluid->ncomp, s->fluid->fluid, s->fluid->pkpm.fluid_no_J
  );
  if (app->use_gpu) {
    // copy data from device to host before writing it out
    gkyl_array_copy(s->fluid->pkpm.fluid_no_J_host, s->fluid->pkpm.fluid_no_J);
  }
  gkyl_comm_array_write(
    app->comm, &app->grid, &app->local, mt_fluid, s->fluid->pkpm.fluid_no_J_host, fileNm_fluid
  );

  // What a restart needs besides the frame, as the stepper holds it: the closures [alpha_0, c]
  // of the last stage of the last completed step, from which the closures of the next step take
  // the speeds of the velocity-domain edges (the staggered fluid state of the operator split is
  // written with the split, vlasov_op_split_write_species).
  cstr fileNm_closure =
    cstr_from_fmt("%s-%s_pkpm_closure_restart_%d.gkyl", app->name, s->name, frame);
  if (app->use_gpu) {
    // copy data from device to host before writing it out
    gkyl_array_copy(s->kinetic->pkpm.closure_step_host, s->kinetic->pkpm.closure_step);
  }
  gkyl_comm_array_write(
    app->comm, &app->grid, &app->local, mt_fluid, s->kinetic->pkpm.closure_step_host,
    fileNm_closure.str
  );
  cstr_drop(&fileNm_closure);

  vlasov_array_meta_release(mt_fluid);

  // the source of the species, with the first frame
  if (frame == 0) {
    vm_pkpm_species_source_write(app, s, &s->kinetic->pkpm.src, tm, frame);
  }

  app->stat.species_io_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_io += 1;
}

// Write the diagnostics of a PKPM species for a frame: the moments, the pressure tensor and heat
// flux, the variables of the update, the constraints and the record of the positivity limiter.
static void
vm_pkpm_species_write_mom(gkyl_vlasov_app *app, struct vlasov_species *s, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  struct gkyl_msgpack_data *mt = vlasov_array_meta_new((struct vlasov_output_meta){
    .frame = frame,
    .stime = tm,
    .poly_order = app->poly_order,
    .basis_type = app->basis.id,
  });

  // Construct the file handles for the three quantities (PKPM moments, pressure tensor and heat
  // flux, PKPM update variables)
  const char *fmt = "%s-%s_pkpm_moms_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, s->name, frame);
  char fileNm[sz + 1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, s->name, frame);

  const char *fmt_prim = "%s-%s_pkpm_prim_%d.gkyl";
  int sz_prim = gkyl_calc_strlen(fmt_prim, app->name, s->name, frame);
  char fileNm_prim[sz_prim + 1]; // ensures no buffer overflow
  snprintf(fileNm_prim, sizeof fileNm_prim, fmt_prim, app->name, s->name, frame);

  const char *fmt_pkpm_vars = "%s-%s_pkpm_vars_%d.gkyl";
  int sz_pkpm_vars = gkyl_calc_strlen(fmt_pkpm_vars, app->name, s->name, frame);
  char fileNm_pkpm_vars[sz_pkpm_vars + 1]; // ensures no buffer overflow
  snprintf(fileNm_pkpm_vars, sizeof fileNm_pkpm_vars, fmt_pkpm_vars, app->name, s->name, frame);

  // Compute the PKPM variables including moments and primitive variables and construct the
  // arrays for writing out the pressure tensor, the heat flux and the other pkpm variables,
  // with b, its gradient and d_t B/|B| of the field being written.
  vm_field_calc_diag_vars(app, app->field);
  vm_pkpm_species_calc_diag_vars(app, s, s->kinetic->f, s->fluid->fluid);
  vm_pkpm_species_moment_diag_calc(
    app, s, s->kinetic->local, app->local, s->kinetic->f, s->fluid->fluid
  );
  // Physical coupling moments for the diagnostics (the evolved ones carry the position-map
  // Jacobian).
  gkyl_vlasov_position_map_divide_jacobpos_conf(
    app->pos_map, &app->local, s->fluid->pkpm.moms->ncomp, s->fluid->pkpm.moms,
    s->fluid->pkpm.moms_no_J
  );
  gkyl_dg_calc_pkpm_vars_diag(
    s->fluid->pkpm.calc_pkpm_vars, &app->local, s->fluid->pkpm.moms_no_J, app->field->bvar,
    app->field->grad_b, s->fluid->pkpm.u, s->kinetic->pkpm.grad_u, s->kinetic->pkpm.theta,
    s->kinetic->pkpm.alpha_0, s->fluid->pkpm.xi, s->kinetic->pkpm.c_h, s->fluid->pkpm.prim_io,
    s->kinetic->pkpm.vars_io
  );

  // copy data from device to host before writing it out
  if (app->use_gpu) {
    gkyl_array_copy(s->kinetic->pkpm.moms_diag.marr_host, s->kinetic->pkpm.moms_diag.marr);
    gkyl_array_copy(s->fluid->pkpm.prim_io_host, s->fluid->pkpm.prim_io);
    gkyl_array_copy(s->kinetic->pkpm.vars_io_host, s->kinetic->pkpm.vars_io);
    gkyl_array_copy(s->kinetic->pkpm.cmoms_host, s->kinetic->pkpm.cmoms);
  }

  gkyl_comm_array_write(
    app->comm, &app->grid, &app->local, mt, s->kinetic->pkpm.moms_diag.marr_host, fileNm
  );
  gkyl_comm_array_write(
    app->comm, &app->grid, &app->local, mt, s->fluid->pkpm.prim_io_host, fileNm_prim
  );
  gkyl_comm_array_write(
    app->comm, &app->grid, &app->local, mt, s->kinetic->pkpm.vars_io_host, fileNm_pkpm_vars
  );
  // The constraints [M0(F), M1(F), (M2(F) + 2 M0(G))/3, int M0_perp dw (x, y, z)] of the normalized
  // distributions, [1, 0, 1, 0, 0, 0] in every mode when the closures hold.
  const char *fmt_cons = "%s-%s_pkpm_constraints_%d.gkyl";
  int sz_cons = gkyl_calc_strlen(fmt_cons, app->name, s->name, frame);
  char fileNm_cons[sz_cons + 1];
  snprintf(fileNm_cons, sizeof fileNm_cons, fmt_cons, app->name, s->name, frame);
  gkyl_comm_array_write(
    app->comm, &app->grid, &app->local, mt, s->kinetic->pkpm.cmoms_host, fileNm_cons
  );
  // The cumulative record of the positivity limiter as the cell averages of a p=1 field:
  // [cells scaled for the density, for the internal energy, inadmissible averages].
  const char *fmt_pos = "%s-%s_pkpm_positivity_%d.gkyl";
  int sz_pos = gkyl_calc_strlen(fmt_pos, app->name, s->name, frame);
  char fileNm_pos[sz_pos + 1];
  snprintf(fileNm_pos, sizeof fileNm_pos, fmt_pos, app->name, s->name, frame);
  int nc = app->basis.num_basis;
  struct gkyl_array *count_host = s->fluid->pkpm.pos_count;
  if (app->use_gpu) {
    count_host = mkarr(false, PKPM_POS_NUM_COUNT, app->local_ext.volume);
    gkyl_array_copy(count_host, s->fluid->pkpm.pos_count);
  }
  gkyl_array_clear(s->fluid->pkpm.pos_io_host, 0.0);
  double c0 = sqrt(pow(2.0, app->cdim)); // constant-mode coefficient of a unit cell average
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local_ext, iter.idx);
    const double *cnt = gkyl_array_cfetch(count_host, loc);
    double *io = gkyl_array_fetch(s->fluid->pkpm.pos_io_host, loc);
    for (int c = 0; c < PKPM_POS_NUM_COUNT; ++c) {
      io[c * nc] = c0 * cnt[c];
    }
  }
  if (app->use_gpu) {
    gkyl_array_release(count_host);
  }
  gkyl_comm_array_write(
    app->comm, &app->grid, &app->local, mt, s->fluid->pkpm.pos_io_host, fileNm_pos
  );

  vlasov_array_meta_release(mt);

  app->stat.species_io_tm += gkyl_time_diff_now_sec(wst);
}

// Integrated diagnostics of a PKPM species: the integrated moments of the fluid state and the
// pressures, the particle loss rates through absorbing walls and the statistics of the
// positivity limiter since the previous call.
static void
vm_pkpm_species_calc_integrated_mom(gkyl_vlasov_app *app, struct vlasov_species *s, double tm)
{
  struct timespec wst = gkyl_wall_clock();
  double avals[9], avals_global[9];

  gkyl_array_clear(s->fluid->integ_mom, 0.0);

  vm_pkpm_species_calc_pkpm_vars(app, s, s->kinetic->f, s->fluid->fluid);
  gkyl_dg_calc_pkpm_vars_int(
    s->fluid->pkpm.calc_pkpm_vars, &app->local, s->fluid->fluid, s->fluid->pkpm.moms,
    s->fluid->pkpm.u, s->fluid->integ_mom
  );
  gkyl_array_scale_range(s->fluid->integ_mom, app->grid.cellVolume, &(app->local));
  if (app->use_gpu) {
    gkyl_array_reduce_range(s->fluid->red_integ_diag, s->fluid->integ_mom, GKYL_SUM, &(app->local));
    gkyl_cu_memcpy(avals, s->fluid->red_integ_diag, sizeof(double[9]), GKYL_CU_MEMCPY_D2H);
  } else {
    gkyl_array_reduce_range(avals, s->fluid->integ_mom, GKYL_SUM, &(app->local));
  }

  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 9, avals, avals_global);
  gkyl_dynvec_append(s->fluid->integ_diag, tm, avals_global);
  // particle loss rates through the absorbing walls of the last stage
  gkyl_dynvec_append(s->kinetic->pkpm.bflux_diag, tm, s->kinetic->pkpm.bflux);

  // Positivity limiter statistics since the previous call: [cells scaled for the
  // density, for the internal energy, cells with an inadmissible average, applications].
  double pos[PKPM_POS_NUM_STAT], pos_global[PKPM_POS_NUM_STAT];
  for (int c = 0; c < PKPM_POS_NUM_STAT; ++c) {
    pos[c] = s->fluid->pkpm.pos_count_since[c];
    s->fluid->pkpm.pos_count_since[c] = 0;
  }
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, PKPM_POS_NUM_STAT, pos, pos_global);
  gkyl_dynvec_append(s->fluid->pkpm.pos_stat, tm, pos_global);

  app->stat.species_diag_calc_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_diag += 1;
}

// Append a record of a species to its file <app>-<species>-<record>.gkyl (a new file at the first
// call of a run), and clear it.
static void
vm_pkpm_species_write_dynvec(
  gkyl_vlasov_app *app, const struct vlasov_species *s, const char *record, gkyl_dynvec dv,
  bool *is_first_write_call
)
{
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (rank == 0) {
    cstr fileNm = cstr_from_fmt("%s-%s-%s.gkyl", app->name, s->name, record);
    if (*is_first_write_call) {
      // write to a new file (this ensure previous output is removed)
      gkyl_dynvec_write(dv, fileNm.str);
      *is_first_write_call = false;
    } else {
      // append to existing file
      gkyl_dynvec_awrite(dv, fileNm.str);
    }
    cstr_drop(&fileNm);
  }
  gkyl_dynvec_clear(dv);
}

// Particle loss rate of a PKPM species through its absorbing walls, from the fluxes of the stage
// (vm_pkpm_species_calc_bflux, with the self moments).
static double
vm_pkpm_species_wall_loss_rate(gkyl_vlasov_app *app, const struct vlasov_species *s)
{
  double tot = 0.0;
  for (int i = 0; i < 2 * app->cdim; ++i) {
    tot += s->kinetic->pkpm.bflux[i];
  }
  return tot;
}

// Write the integrated diagnostics of a PKPM species: the integrated moments, and on the same
// cadence the wall loss rates and the statistics of the positivity limiter.
static void
vm_pkpm_species_write_integrated_mom(gkyl_vlasov_app *app, struct vlasov_species *s)
{
  struct timespec wst = gkyl_wall_clock();
  vm_pkpm_species_write_dynvec(
    app, s, "imom", s->fluid->integ_diag, &s->fluid->is_first_integ_write_call
  );
  vm_pkpm_species_write_dynvec(
    app, s, "bflux", s->kinetic->pkpm.bflux_diag, &s->kinetic->pkpm.is_first_bflux_write_call
  );
  vm_pkpm_species_write_dynvec(
    app, s, "positivity-stat", s->fluid->pkpm.pos_stat, &s->fluid->pkpm.is_first_pos_stat_write_call
  );
  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_diag_io += 1;
}

// Write the integrated L2 norm of the distributions of a PKPM species.
static void
vm_pkpm_kinetic_species_write_L2(gkyl_vlasov_app *app, struct vm_species *kin)
{
  struct timespec wst = gkyl_wall_clock();
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (rank == 0) {
    cstr fileNm = cstr_from_fmt("%s-%s-%s.gkyl", app->name, kin->name, "L2");
    if (kin->is_first_integ_L2_write_call) {
      // write to a new file (this ensure previous output is removed)
      gkyl_dynvec_write(kin->integ_L2_f, fileNm.str);
      kin->is_first_integ_L2_write_call = false;
    } else {
      // append to existing file
      gkyl_dynvec_awrite(kin->integ_L2_f, fileNm.str);
    }
    cstr_drop(&fileNm);
  }
  gkyl_dynvec_clear(kin->integ_L2_f);
  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_diag_io += 1;
}

// Read the distributions of a PKPM species from a file (the physical ones; the evolved ones carry
// the position-map Jacobian), with their boundary conditions.
static struct gkyl_app_restart_status
vm_pkpm_kinetic_species_from_file(gkyl_vlasov_app *app, struct vm_species *kin, const char *fname)
{
  struct gkyl_app_restart_status rstat = vlasov_header_from_file(app, fname);

  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    rstat.io_status = gkyl_comm_array_read(kin->comm, &kin->grid, &kin->local, kin->f_host, fname);
    if (app->use_gpu) {
      gkyl_array_copy(kin->f, kin->f_host);
    }
    if (GKYL_ARRAY_RIO_SUCCESS == rstat.io_status) {
      // The file holds the physical distributions; the evolved ones carry the Jacobian.
      gkyl_vlasov_position_map_rescale_jacobpos(
        app->pos_map, &kin->basis, &kin->local, kin->f, kin->f
      );
      vm_pkpm_species_apply_bc(app, kin, kin->f);
    }
  }

  return rstat;
}

// Read the fluid state of a PKPM species from a file (the physical one; the evolved one carries
// the position-map Jacobian), with its boundary conditions.
static struct gkyl_app_restart_status
vm_pkpm_fluid_species_from_file(gkyl_vlasov_app *app, struct vm_fluid_species *fl, const char *fname)
{
  struct gkyl_app_restart_status rstat = vlasov_header_from_file(app, fname);

  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    rstat.io_status =
      gkyl_comm_array_read(app->comm, &app->grid, &app->local, fl->fluid_host, fname);
    if (app->use_gpu) {
      gkyl_array_copy(fl->fluid, fl->fluid_host);
    }
    if (GKYL_ARRAY_RIO_SUCCESS == rstat.io_status) {
      // The file holds the physical state; the evolved one carries the Jacobian.
      gkyl_vlasov_position_map_rescale_jacobpos_conf(
        app->pos_map, &app->local, fl->fluid, fl->fluid
      );
      vm_pkpm_fluid_species_apply_bc(app, fl, fl->fluid);
    }
  }

  return rstat;
}

// Restart a PKPM species from a frame: its distributions, its fluid state, the closures of the
// last completed step and the record of the positivity limiter. The diagnostic records then
// append to their files and the app clock is set to the frame time.
static struct gkyl_app_restart_status
vm_pkpm_species_read_from_frame(gkyl_vlasov_app *app, struct vlasov_species *s, int frame)
{
  // Seed the species from the initial conditions before the read overwrites the interior:
  // this fills the buffers of the fixed-function boundary conditions and the cells no
  // boundary condition fills (the corners of the extended range). The limiter of the initial
  // fluid state is not applied: its records are those of the run.
  vm_pkpm_kinetic_species_apply_ic(app, s->kinetic, 0.0);
  vm_pkpm_fluid_species_project_ic(app, s->fluid, 0.0);

  cstr fileNm = cstr_from_fmt("%s-%s_%d.gkyl", app->name, s->name, frame);
  struct gkyl_app_restart_status rstat =
    vm_pkpm_kinetic_species_from_file(app, s->kinetic, fileNm.str);
  cstr_drop(&fileNm);

  // The applied acceleration is not part of the frame, and it is not recomputed in the
  // time-stepping loop when time-independent.
  vm_fluid_species_calc_app_accel(app, s->fluid, rstat.stime);

  cstr fileNm_fluid = cstr_from_fmt("%s-%s_fluid_%d.gkyl", app->name, s->name, frame);
  struct gkyl_app_restart_status rstat_fluid =
    vm_pkpm_fluid_species_from_file(app, s->fluid, fileNm_fluid.str);
  cstr_drop(&fileNm_fluid);

  // The closures the first step takes the speeds of its velocity-domain edges from (a frame
  // without the file starts them from zero, as a new run does).
  cstr fileNm_closure =
    cstr_from_fmt("%s-%s_pkpm_closure_restart_%d.gkyl", app->name, s->name, frame);
  struct gkyl_app_restart_status rstat_closure = vlasov_header_from_file(app, fileNm_closure.str);
  if (rstat_closure.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    rstat_closure.io_status = gkyl_comm_array_read(
      app->comm, &app->grid, &app->local, s->kinetic->pkpm.closure_step_host, fileNm_closure.str
    );
  }
  if (rstat_closure.io_status == GKYL_ARRAY_RIO_SUCCESS && app->use_gpu) {
    gkyl_array_copy(s->kinetic->pkpm.closure_step, s->kinetic->pkpm.closure_step_host);
  }
  cstr_drop(&fileNm_closure);

  // The cumulative record of the positivity limiter continues from the frame's (a frame
  // without the file restarts the record).
  cstr fileNm_pos = cstr_from_fmt("%s-%s_pkpm_positivity_%d.gkyl", app->name, s->name, frame);
  struct gkyl_app_restart_status rstat_pos = vlasov_header_from_file(app, fileNm_pos.str);
  if (rstat_pos.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    rstat_pos.io_status = gkyl_comm_array_read(
      app->comm, &app->grid, &app->local, s->fluid->pkpm.pos_io_host, fileNm_pos.str
    );
  }
  if (rstat_pos.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    int nc = app->basis.num_basis;
    struct gkyl_array *count_host = s->fluid->pkpm.pos_count;
    if (app->use_gpu) {
      count_host = mkarr(false, PKPM_POS_NUM_COUNT, app->local_ext.volume);
      gkyl_array_clear(count_host, 0.0);
    }
    double c0 = sqrt(pow(2.0, app->cdim)); // constant-mode coefficient of a unit cell average
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &app->local);
    while (gkyl_range_iter_next(&iter)) {
      long loc = gkyl_range_idx(&app->local_ext, iter.idx);
      const double *io = gkyl_array_cfetch(s->fluid->pkpm.pos_io_host, loc);
      double *cnt = gkyl_array_fetch(count_host, loc);
      for (int c = 0; c < PKPM_POS_NUM_COUNT; ++c) {
        cnt[c] = io[c * nc] / c0;
      }
    }
    if (app->use_gpu) {
      gkyl_array_copy(s->fluid->pkpm.pos_count, count_host);
      gkyl_array_release(count_host);
    }
  }
  cstr_drop(&fileNm_pos);

  // Diagnostic records append to their existing files after a restart.
  s->fluid->is_first_integ_write_call = false;
  s->kinetic->is_first_integ_L2_write_call = false;
  s->kinetic->pkpm.is_first_bflux_write_call = false;
  s->fluid->pkpm.is_first_pos_stat_write_call = false;

  // The first failure is reported (a frame without the closures or the positivity record is
  // not one).
  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    rstat.io_status = rstat_fluid.io_status;
  }
  // The app clock resumes at the frame time.
  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    app->tcurr = rstat.stime;
  }

  return rstat;
}

void
vlasov_pkpm_species_stat_write(gkyl_vlasov_app *app, struct vlasov_species *sp, FILE *fp)
{
  double null_loc = (double)sp->kinetic->pkpm.num_null_closure, null_glob = 0.0,
         cond_loc = sp->kinetic->pkpm.closure_cond_min, cond_glob = 1.0;
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &null_loc, &null_glob);
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MIN, 1, &cond_loc, &cond_glob);
  gkyl_vlasov_app_cout(app, fp, " %s_closure_null_nodes : %ld,\n", sp->name, (long)null_glob);
  gkyl_vlasov_app_cout(app, fp, " %s_closure_cond_min : %lg,\n", sp->name, cond_glob);
}

// release resources for PKPM species
void
vm_pkpm_species_release(const gkyl_vlasov_app *app, struct vlasov_species *s)
{
  // release various arrays
  gkyl_array_release(s->kinetic->f);
  gkyl_array_release(s->kinetic->f_no_J);
  gkyl_array_release(s->kinetic->f1);
  gkyl_array_release(s->kinetic->fnew);
  if (s->kinetic->pkpm.has_nodal_ops) {
    gkyl_array_release(s->kinetic->pkpm.f_nodal);
    gkyl_array_release(s->kinetic->pkpm.rhs_nodal);
  }
  gkyl_array_release(s->kinetic->pkpm.blk_in);
  gkyl_array_release(s->kinetic->pkpm.blk_out);
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
    gkyl_array_release(s->kinetic->pkpm.f_no_J_host);
    gkyl_array_release(s->fluid->fluid_host);
    gkyl_array_release(s->fluid->pkpm.fluid_no_J_host);
    gkyl_array_release(s->fluid->pkpm.prim_io_host);
    gkyl_array_release(s->kinetic->pkpm.vars_io_host);
  }

  gkyl_array_release(s->fluid->qmem);

  // release moment data
  vm_pkpm_species_moment_release(app, &s->kinetic->pkpm.moms_diag);
  vm_pkpm_species_moment_release(app, &s->kinetic->pkpm.nmoms);

  // release pkpm arrays data and updaters
  gkyl_array_release(s->fluid->pkpm.u);
  gkyl_array_release(s->kinetic->pkpm.u);
  gkyl_array_release(s->kinetic->pkpm.grad_u);
  gkyl_array_release(s->kinetic->pkpm.theta);
  gkyl_array_release(s->kinetic->pkpm.alpha_0);
  gkyl_array_release(s->kinetic->pkpm.c_h);
  gkyl_array_release(s->fluid->pkpm.xi);
  gkyl_array_release(s->kinetic->pkpm.xi);
  gkyl_array_release(s->fluid->pkpm.xi3);
  gkyl_array_release(s->kinetic->pkpm.grad_xi);
  gkyl_array_release(s->fluid->pkpm.scratch);
  gkyl_array_release(s->fluid->pkpm.moms);
  gkyl_array_release(s->fluid->pkpm.moms_no_J);
  gkyl_array_release(s->fluid->pkpm.fluid_no_J);
  gkyl_array_release(s->kinetic->pkpm.cmoms);
  if (app->use_gpu) {
    gkyl_array_release(s->kinetic->pkpm.cmoms_host);
  }
  gkyl_array_release(s->kinetic->pkpm.rhs_moms);
  gkyl_array_release(s->kinetic->pkpm.conf_flux_surf);
  gkyl_array_release(s->kinetic->pkpm.flux_blend);
  gkyl_array_release(s->kinetic->pkpm.vel_flux_surf);
  gkyl_array_release(s->fluid->pkpm.fluid_flux_surf);
  gkyl_array_release(s->fluid->pkpm.rho);
  gkyl_array_release(s->fluid->pkpm.rho3);
  gkyl_array_release(s->fluid->pkpm.mom3);
  if (s->fluid->has_diffusion) {
    gkyl_dg_pkpm_mom_diffusion_release(s->fluid->pkpm.mom_diff);
  }
  gkyl_array_release(s->fluid->pkpm.grad_rho);
  gkyl_array_release(s->kinetic->pkpm.grad_ln_n);
  gkyl_array_release(s->kinetic->pkpm.c_0);
  gkyl_array_release(s->kinetic->pkpm.a_perp);
  if (app->use_gpu) {
    gkyl_array_release(s->kinetic->pkpm.closure_step_host);
  }
  gkyl_array_release(s->kinetic->pkpm.closure_step);
  gkyl_array_release(s->fluid->integ_mom);
  gkyl_array_release(s->fluid->pkpm.prim_io);
  gkyl_array_release(s->kinetic->pkpm.vars_io);
  gkyl_array_release(s->fluid->pkpm.pos_count);
  gkyl_array_release(s->fluid->pkpm.pos_count_saved);
  gkyl_array_release(s->fluid->pkpm.pos_io);
  if (app->use_gpu) {
    gkyl_array_release(s->fluid->pkpm.pos_io_host);
  }
  gkyl_dg_calc_pkpm_vars_release(s->fluid->pkpm.calc_pkpm_vars);
  gkyl_dg_calc_pkpm_vars_release(s->kinetic->pkpm.calc_pkpm_vars);
  gkyl_dg_pkpm_remap_flux_surf_release(s->kinetic->pkpm.calc_face_flux);
  gkyl_dg_pkpm_vel_flux_surf_release(s->kinetic->pkpm.calc_vel_flux);
  gkyl_dg_pkpm_harm_vol_release(s->kinetic->pkpm.harm_vol);
  gkyl_array_release(s->kinetic->pkpm.harm_coef);
  gkyl_array_release(s->kinetic->pkpm.harm_wnode);
  gkyl_dg_pkpm_harm_exp_release(s->kinetic->pkpm.harm_exp);
  gkyl_array_release(s->kinetic->pkpm.harm_exp_coef);
  gkyl_array_release(s->kinetic->pkpm.harm_exp_nu);
  gkyl_array_release(s->kinetic->pkpm.harm_exp_act);
  gkyl_array_release(s->kinetic->pkpm.harm_exp_col);
  gkyl_hyper_dg_release(s->kinetic->pkpm.slvr_kinetic_conf);
  gkyl_hyper_dg_release(s->kinetic->pkpm.slvr_kinetic_vel);
  gkyl_hyper_dg_release(s->fluid->pkpm.slvr_fluid);
  gkyl_dg_eqn_release(s->kinetic->pkpm.eqn_kinetic_conf);
  gkyl_dg_eqn_release(s->kinetic->pkpm.eqn_kinetic_vel);
  gkyl_dg_eqn_release(s->fluid->pkpm.eqn_fluid);

  gkyl_array_release(s->kinetic->L2_f);
  gkyl_dynvec_release(s->kinetic->integ_L2_f);
  gkyl_dynvec_release(s->fluid->integ_diag);
  gkyl_dynvec_release(s->fluid->pkpm.pos_stat);
  gkyl_dynvec_release(s->kinetic->pkpm.bflux_diag);
  if (s->fluid->pkpm.fluid_flux_surf_host) {
    gkyl_array_release(s->fluid->pkpm.fluid_flux_surf_host);
  }
  vm_pkpm_species_source_release(app, &s->kinetic->pkpm.src);
  vm_pkpm_species_source_bgk_release(app, &s->kinetic->pkpm.src_bgk);
  vm_pkpm_species_coll_cross_release(app, &s->kinetic->pkpm.cross);
  if (s->kinetic->pkpm.has_stage_weights) {
    gkyl_array_release(s->kinetic->pkpm.stage_lam);
    gkyl_pkpm_stage_weights_release(s->kinetic->pkpm.up_stage_weights);
  }

  gkyl_array_release(s->fluid->app_accel);
  if (s->fluid->has_app_accel) {
    if (app->use_gpu) {
      gkyl_array_release(s->fluid->app_accel_host);
    }
    gkyl_proj_on_basis_release(s->fluid->app_accel_proj);
  }

  if (s->kinetic->pkpm.collision_id == GKYL_BGK_COLLISIONS) {
    vm_pkpm_species_bgk_release(app, &s->kinetic->pkpm.bgk);
  }
  if (s->kinetic->pkpm.collision_id == GKYL_LBO_COLLISIONS) {
    vm_pkpm_species_lbo_release(app, &s->kinetic->pkpm.lbo);
  }

  // Copy BCs are allocated by default. Need to free.
  for (int d = 0; d < app->cdim; ++d) {
    gkyl_bc_basic_release(s->kinetic->bc_lo[d]);
    gkyl_bc_basic_release(s->kinetic->bc_up[d]);
    if (s->kinetic->pkpm.bc_lo_reflect[d]) {
      gkyl_bc_pkpm_species_reflect_release(s->kinetic->pkpm.bc_lo_reflect[d]);
    }
    if (s->kinetic->pkpm.bc_up_reflect[d]) {
      gkyl_bc_pkpm_species_reflect_release(s->kinetic->pkpm.bc_up_reflect[d]);
    }
    gkyl_bc_basic_release(s->fluid->bc_lo[d]);
    gkyl_bc_basic_release(s->fluid->bc_up[d]);
    if (s->kinetic->pkpm.bc_lo_theta[d]) {
      gkyl_bc_basic_release(s->kinetic->pkpm.bc_lo_theta[d]);
    }
    if (s->kinetic->pkpm.bc_up_theta[d]) {
      gkyl_bc_basic_release(s->kinetic->pkpm.bc_up_theta[d]);
    }
  }

  if (app->use_gpu) {
    gkyl_cu_free(s->kinetic->basis_on_dev);
    gkyl_cu_free(s->kinetic->omega_cfl);
    gkyl_cu_free(s->fluid->omegaCfl_ptr);
    gkyl_cu_free(s->kinetic->red_L2_f);
    gkyl_cu_free(s->fluid->red_integ_diag);
  } else {
    gkyl_free(s->kinetic->omega_cfl);
    gkyl_free(s->fluid->omegaCfl_ptr);
  }

  if (s->fluid->pkpm.src_bgk_rhs) {
    gkyl_array_release(s->fluid->pkpm.src_bgk_rhs);
  }

  gkyl_free(s->kinetic);
  gkyl_free(s->fluid);
}
