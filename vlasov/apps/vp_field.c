#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_bc_basic.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_priv.h>

#include <assert.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <time.h>

// Configuration-space c2p for external-field/potential projection: map the
// projection's computational quadrature coordinates to physical ones via the
// position map, so the user-supplied function (defined in physical space) is
// sampled at the correct locations on a non-uniform conf mesh. For an identity
// map eval_mc2p is the identity, leaving uniform-grid behavior unchanged.
static void
vp_field_ext_c2p(const double *xcomp, double *xphys, void *ctx)
{
  struct vm_field_proj_c2p_ctx *c = ctx;
  gkyl_vlasov_position_map_eval_mc2p(c->pos_map, xcomp, xphys);
}

// --- Species coupling laws of the potentials ----------------------------------

// Electrostatic potential: source q*n, acceleration -(q/m)*grad(phi).
static double
vp_field_es_src_weight(const struct vp_potential *pot, double charge, double mass)
{
  return charge;
}

static double
vp_field_es_force_weight(const struct vp_potential *pot, double charge, double mass)
{
  return charge / mass;
}

// Gravitational potential: source alpha_g*m*n, acceleration -grad(phi_g) for
// every species.
static double
vp_field_grav_src_weight(const struct vp_potential *pot, double charge, double mass)
{
  return pot->coupling * mass;
}

static double
vp_field_grav_force_weight(const struct vp_potential *pot, double charge, double mass)
{
  return 1.0;
}

// --- Metric tensors on a mapped configuration mesh ---------------------------

// Number of components of a symmetric conf-space tensor: 1 when it is a
// constant scalar (identity position map), otherwise 1x->1, 2x->3, 3x->6.
static int
vp_field_tensor_num_comp(const struct gkyl_vlasov_app *app, bool is_const)
{
  int cdim = app->cdim;
  return is_const ? 1 : cdim + (int)ceil((pow(3.0, cdim - 1) - cdim) / 2.0);
}

// Set a per-cell constant diagonal tensor (symmetric upper-triangular storage,
// off-diagonals 0) on the local range of a host array: component (i,i) =
// scale * J^jpow / J_xi^2, with J_xi the position-map Jacobian in direction i
// and J = prod_i J_xi the total conf Jacobian.
static void
vp_field_set_metric_diag(
  const struct gkyl_vlasov_app *app, struct gkyl_array *arr_ho, double scale, int jpow
)
{
  int cdim = app->cdim, nb = app->basis.num_basis;
  double dg0 = pow(sqrt(2.0), cdim); // 0th DG coeff representing a constant value.
  int stride = app->basis.poly_order + 1; // per-direction block stride in jacob_pos.
  gkyl_array_clear(arr_ho, 0.0);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  while (gkyl_range_iter_next(&iter)) {
    long cidx = gkyl_range_idx(&app->local, iter.idx);
    const double *jacob_pos = gkyl_array_cfetch(app->pos_map->jacob_pos_host, cidx);
    const double *jacob_pos_gauss = gkyl_array_cfetch(app->pos_map->jacob_pos_gauss_host, cidx);
    double Jtot = jacob_pos_gauss[0];
    double *arr_d = gkyl_array_fetch(arr_ho, cidx);
    for (int i = 0; i < cdim; ++i) {
      double Jxi = jacob_pos[i * stride];
      int diag =
        i * cdim - (i * (i - 1)) / 2; // index of (i,i) in symmetric upper-triangular storage.
      arr_d[diag * nb] = scale * pow(Jtot, jpow) / (Jxi * Jxi) * dg0;
    }
  }
}

// Permittivity of a Poisson solve with scalar permittivity 'scale', on the
// global extended range the FEM solver reads. On a uniform grid (identity
// position map) this is the constant scalar. On a non-uniform conf mesh it is
// the diagonal metric permittivity tensor
//   eps^{ii} = scale * J / J_xi^2   (off-diagonals 0, diagonal position map)
// which makes the weak Poisson operator int eps^{ij} d_i(phi) d_j(psi) dxi equal
// the physical int scale * grad(phi).grad(psi) dx; the source already carries J.
static struct gkyl_array *
vp_field_new_epsilon(struct gkyl_vlasov_app *app, double scale)
{
  int nb = app->basis.num_basis;
  bool eps_const = app->pos_map->is_identity;
  int epsnum = vp_field_tensor_num_comp(app, eps_const);

  struct gkyl_array *epsilon = mkarr(app->use_gpu, epsnum * nb, app->global_ext.volume);
  gkyl_array_clear(epsilon, 0.0);
  if (eps_const) {
    gkyl_array_shiftc(epsilon, scale * pow(sqrt(2.0), app->cdim), 0);
  } else {
    // Build the local diagonal tensor from the position map, then allgather to
    // the global array the FEM solver reads (mirrors the source allgather).
    struct gkyl_array *eps_local = mkarr(app->use_gpu, epsnum * nb, app->local_ext.volume);
    struct gkyl_array *eps_local_ho = app->use_gpu ?
                                        mkarr(false, epsnum * nb, app->local_ext.volume) :
                                        gkyl_array_acquire(eps_local);
    vp_field_set_metric_diag(app, eps_local_ho, scale, 1);
    if (app->use_gpu) {
      gkyl_array_copy(eps_local, eps_local_ho);
    }
    gkyl_comm_array_allgather(app->comm, &app->local, &app->global, eps_local, epsilon);
    gkyl_array_release(eps_local);
    gkyl_array_release(eps_local_ho);
  }
  return epsilon;
}

// --- One potential: allocation, solve, diagnostics, release -------------------

// Allocate the source, potential, permittivity (scalar eps_scale), Poisson
// solver and energy diagnostic of a potential; its names and coupling laws are
// set by the caller. All potentials share the field's Poisson boundary
// conditions.
static void
vp_potential_init(
  struct gkyl_vlasov_app *app, struct vm_field *field, struct vp_potential *pot, double eps_scale
)
{
  int nb = app->basis.num_basis;

  pot->rho = mkarr(app->use_gpu, nb, app->local_ext.volume);
  pot->rho_global = mkarr(app->use_gpu, nb, app->global_ext.volume);

  pot->phi = mkarr(app->use_gpu, nb, app->local_ext.volume);
  pot->phi_global = mkarr(app->use_gpu, nb, app->global_ext.volume);
  pot->phi_host = app->use_gpu ? mkarr(false, pot->phi->ncomp, pot->phi->size) :
                                 gkyl_array_acquire(pot->phi);

  pot->epsilon = vp_field_new_epsilon(app, eps_scale);
  pot->fem_poisson = gkyl_fem_poisson_new(
    &app->global, &app->grid, app->basis, &field->info.poisson_bcs, NULL, pot->epsilon, NULL,
    app->pos_map->is_identity, app->use_gpu
  );

  if (app->use_gpu) {
    pot->energy_red = gkyl_cu_malloc(sizeof(double[1]));
    pot->energy_red_global = gkyl_cu_malloc(sizeof(double[1]));
  } else {
    pot->energy_red = gkyl_malloc(sizeof(double[1]));
    pot->energy_red_global = gkyl_malloc(sizeof(double[1]));
  }
  pot->integ_energy = gkyl_dynvec_new(GKYL_DOUBLE, 1);
}

static void
vp_potential_release(const struct gkyl_vlasov_app *app, struct vp_potential *pot)
{
  gkyl_dynvec_release(pot->integ_energy);
  if (app->use_gpu) {
    gkyl_cu_free(pot->energy_red);
    gkyl_cu_free(pot->energy_red_global);
  } else {
    gkyl_free(pot->energy_red);
    gkyl_free(pot->energy_red_global);
  }

  gkyl_fem_poisson_release(pot->fem_poisson);
  gkyl_array_release(pot->epsilon);

  gkyl_array_release(pot->phi_host);
  gkyl_array_release(pot->phi);
  gkyl_array_release(pot->phi_global);

  gkyl_array_release(pot->rho_global);
  gkyl_array_release(pot->rho);
}

// Solve one potential from its accumulated source: gather the source into the
// global array, solve the Poisson problem, and copy the portion of the global
// potential corresponding to this MPI process to the local potential.
static void
vp_potential_solve(gkyl_vlasov_app *app, struct vm_field *field, struct vp_potential *pot)
{
  gkyl_comm_array_allgather(app->comm, &app->local, &app->global, pot->rho, pot->rho_global);

  gkyl_fem_poisson_set_rhs(pot->fem_poisson, pot->rho_global, NULL);
  gkyl_fem_poisson_solve(pot->fem_poisson, pot->phi_global);

  gkyl_array_copy_range_to_range(pot->phi, pot->phi_global, &app->local, &field->global_sub_range);
}

// Append the integrated |grad phi|^2 (scaled by the potential's energy factor)
// at time tm to the potential's energy diagnostic.
static void
vp_potential_calc_energy(
  gkyl_vlasov_app *app, const struct vm_field *field, struct vp_potential *pot, double tm
)
{
  // The integrate updater already accounts for the cell volume, so the factor is 1.
  gkyl_array_integrate_advance(
    field->calc_grad_sq, pot->phi, 1.0, field->grad_sq_wgt, &app->local, &app->local,
    pot->energy_red
  );

  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, pot->energy_red, pot->energy_red_global);

  double energy_global[1] = {0.0};
  if (app->use_gpu) {
    gkyl_cu_memcpy(energy_global, pot->energy_red_global, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
  } else {
    energy_global[0] = pot->energy_red_global[0];
  }
  energy_global[0] *= pot->energy_fac;

  gkyl_dynvec_append(pot->integ_energy, tm, energy_global);
}

// Write the potential's frame, <app>-<frame_name>_<frame>.gkyl.
static void
vp_potential_write(
  gkyl_vlasov_app *app, struct vp_potential *pot, struct gkyl_msgpack_data *mt, int frame
)
{
  const char *fmt = "%s-%s_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, pot->frame_name, frame);
  char fileNm[sz + 1]; // Ensures no buffer overflow.
  snprintf(fileNm, sizeof fileNm, fmt, app->name, pot->frame_name, frame);

  // Copy data from device to host before writing it out.
  if (app->use_gpu) {
    gkyl_array_copy(pot->phi_host, pot->phi);
  }
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, pot->phi_host, fileNm);
}

// Write (first call) or append (later calls) the potential's energy diagnostic,
// <app>-<energy_name>.gkyl, and clear it.
static void
vp_potential_write_energy(gkyl_vlasov_app *app, struct vp_potential *pot, bool is_first_write)
{
  const char *fmt = "%s-%s.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, pot->energy_name);
  char fileNm[sz + 1]; // Ensures no buffer overflow.
  snprintf(fileNm, sizeof fileNm, fmt, app->name, pot->energy_name);

  int rank;
  gkyl_comm_get_rank(app->comm, &rank);

  if (rank == 0) {
    if (is_first_write) {
      // Write to a new file (this ensure previous output is removed).
      gkyl_dynvec_write(pot->integ_energy, fileNm);
    } else {
      // Append to existing file.
      gkyl_dynvec_awrite(pot->integ_energy, fileNm);
    }
  }
  gkyl_dynvec_clear(pot->integ_energy);
}

// --- Field ---------------------------------------------------------------------

struct vm_field *
vp_field_new(struct gkyl_vm *vm, struct gkyl_vlasov_app *app)
{
  // Initialize field object.
  struct vm_field *vpf = gkyl_malloc(sizeof(struct vm_field));
  *vpf = (struct vm_field){};

  vpf->info = vm->field;
  vpf->field_id = GKYL_FIELD_PHI;

  int cdim = app->cdim, nb = app->basis.num_basis;
  double dg0 = pow(sqrt(2.0), cdim); // 0th DG coeff representing a constant value.

  // Create global subrange we'll copy the field solver solution from (into local).
  int intersect = gkyl_sub_range_intersect(&vpf->global_sub_range, &app->global, &app->local);

  // Potentials to solve for: the electrostatic potential needs a positive
  // permittivity and a charged species; the gravitational potential needs a
  // positive gravitational coupling and a species with mass. Either, both, or
  // neither (external potentials only) may be present.
  bool has_charge = false, has_mass = false;
  for (int i = 0; i < vm->num_species; ++i) {
    has_charge = has_charge || vm->species[i].charge != 0.0;
    has_mass = has_mass || vm->species[i].mass != 0.0;
  }

  vpf->num_pots = 0;
  if (vpf->info.epsilon0 > 0.0 && has_charge) {
    // Electrostatics: - nabla . (epsilon0 * nabla phi) = sum_s q_s n_s.
    struct vp_potential *pot = &vpf->pots[vpf->num_pots++];
    strcpy(pot->frame_name, "field");
    strcpy(pot->energy_name, "field-energy");
    pot->coupling = vpf->info.epsilon0;
    pot->energy_fac = 1.0; // Diagnostic is int |grad phi|^2 (no epsilon0/2).
    pot->src_weight = vp_field_es_src_weight;
    pot->force_weight = vp_field_es_force_weight;
    vp_potential_init(app, vpf, pot, vpf->info.epsilon0);
  }
  if (vpf->info.alpha_g > 0.0 && has_mass) {
    // Self-gravity: nabla^2 phi_g = alpha_g * sum_s m_s n_s, i.e. a permittivity
    // of -1 with the gravitational coupling folded into the source.
    struct vp_potential *pot = &vpf->pots[vpf->num_pots++];
    strcpy(pot->frame_name, "field_grav");
    strcpy(pot->energy_name, "field-grav-energy");
    pot->coupling = vpf->info.alpha_g;
    // Diagnostic is int |grad phi_g|^2 / alpha_g; -1/2 of it is the
    // gravitational potential energy.
    pot->energy_fac = 1.0 / vpf->info.alpha_g;
    pot->src_weight = vp_field_grav_src_weight;
    pot->force_weight = vp_field_grav_force_weight;
    vp_potential_init(app, vpf, pot, -1.0);
  }

  // Coordinate map for projecting external fields/potentials: their user
  // functions are defined in physical space, so the projection must map its
  // computational quadrature points through the position map first. Identity map
  // -> identity c2p, so uniform grids are unaffected.
  vpf->ext_c2p_ctx = (struct vm_field_proj_c2p_ctx){.pos_map = app->pos_map};

  // Initialize external potentials.
  vpf->ext_pot = mkarr(app->use_gpu, 4 * app->basis.num_basis, app->local_ext.volume);
  gkyl_array_clear(vpf->ext_pot, 0.0);
  vpf->has_ext_pot = false;
  vpf->ext_pot_evolve = false;
  // setup external electromagnetic field
  if (vpf->info.external_potentials) {
    vpf->has_ext_pot = true;
    if (vpf->info.external_potentials_evolve) {
      vpf->ext_pot_evolve = vpf->info.external_potentials_evolve;
    }

    vpf->ext_pot_host = app->use_gpu ? mkarr(false, vpf->ext_pot->ncomp, vpf->ext_pot->size) :
                                       gkyl_array_acquire(vpf->ext_pot);
    // Project on the physical coordinates of the (possibly mapped) conf mesh.
    // 4 components: phi and the three components of A.
    vpf->ext_pot_proj = gkyl_eval_on_nodes_inew(&(struct gkyl_eval_on_nodes_inp){
      .grid = &app->grid,
      .basis = &app->basis,
      .num_ret_vals = 4,
      .eval = vpf->info.external_potentials,
      .ctx = vpf->info.external_potentials_ctx,
      .c2p_func = vp_field_ext_c2p,
      .c2p_func_ctx = &vpf->ext_c2p_ctx,
    });
  }

  // Initialize external EM fields.
  vpf->ext_em = mkarr(app->use_gpu, 6 * app->basis.num_basis, app->local_ext.volume);
  gkyl_array_clear(vpf->ext_em, 0.0);
  vpf->has_ext_em = false;
  vpf->ext_em_evolve = false;
  // setup external electromagnetic field
  if (vpf->info.ext_em) {
    vpf->has_ext_em = true;
    if (vpf->info.ext_em_evolve) {
      vpf->ext_em_evolve = vpf->info.ext_em_evolve;
    }

    vpf->ext_em_host = app->use_gpu ? mkarr(false, vpf->ext_em->ncomp, vpf->ext_em->size) :
                                      gkyl_array_acquire(vpf->ext_em);
    // Project on the physical coordinates of the (possibly mapped) conf mesh.
    vpf->ext_em_proj = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
      .grid = &app->grid,
      .basis = &app->basis,
      .qtype = GKYL_GAUSS_QUAD,
      .num_quad = app->basis.poly_order + 1,
      .num_ret_vals = 6,
      .eval = vpf->info.ext_em,
      .ctx = vpf->info.ext_em_ctx,
      .c2p_func = vp_field_ext_c2p,
      .c2p_func_ctx = &vpf->ext_c2p_ctx,
    });
  }

  // Vlasov-Poisson doesn't presently use external currents or limiters.
  vpf->has_app_current = vpf->app_current_evolve = false;
  vpf->limit_em = false;

  // Integrated |grad phi|^2 energy diagnostic shared by the potentials. For an
  // identity position map this is the plain |grad phi|^2 operator (GRAD_SQ) with
  // weight 1. For a mapped conf mesh the physical
  //   int |grad_x phi|^2 dx = int sum_i (J/J_xi^2)(d_xi phi)^2 dxi,
  // i.e. a diagonal metric (J = prod_i J_xi, per-cell constant) folded into the
  // full-gradient weighted operator (EPS_GRAD_SQ). The weight is J^2/J_xi^2: one
  // factor of J = sqrt(g) folds the physical volume into the metric
  // (g^ii = 1/J_xi^2), and a second factor of J corrects the integrate operator's
  // (computational) cell-volume factor to the physical cell volume J*dxi_comp.
  // With this weight EPS_GRAD_SQ matches the uniform-grid GRAD_SQ diagnostic
  // exactly for a constant-Jacobian map (verified for all cdim) and reduces to
  // GRAD_SQ for the identity map. The weight uses the same symmetric
  // upper-triangular layout as the Poisson permittivity tensor.
  bool wgt_const = app->pos_map->is_identity;
  int wgtnum = vp_field_tensor_num_comp(app, wgt_const);
  vpf->grad_sq_wgt = mkarr(app->use_gpu, wgtnum * nb, app->local_ext.volume);
  gkyl_array_clear(vpf->grad_sq_wgt, 0.0);
  if (wgt_const) {
    gkyl_array_shiftc(vpf->grad_sq_wgt, dg0, 0); // Sets the weight to 1 (unused by GRAD_SQ).
  } else {
    struct gkyl_array *wgt_ho = app->use_gpu ? mkarr(false, wgtnum * nb, app->local_ext.volume) :
                                               gkyl_array_acquire(vpf->grad_sq_wgt);
    vp_field_set_metric_diag(app, wgt_ho, 1.0, 2);
    if (app->use_gpu) {
      gkyl_array_copy(vpf->grad_sq_wgt, wgt_ho);
    }
    gkyl_array_release(wgt_ho);
  }

  vpf->calc_grad_sq = gkyl_array_integrate_new(
    &app->grid, &app->basis, 1,
    wgt_const ? GKYL_ARRAY_INTEGRATE_OP_GRAD_SQ : GKYL_ARRAY_INTEGRATE_OP_EPS_GRAD_SQ, app->use_gpu
  );
  vpf->is_first_energy_write_call = true;

  // Set the type-specific dispatch methods (Vlasov-Poisson). The potentials are
  // re-solved each stage (update_func) rather than carried in the RK state, so
  // the combine/copy/BC/current/limiter stage operations are no-ops.
  vpf->update_func = vp_field_update;
  vpf->combine_func = vp_field_combine;
  vpf->copy_range_func = vp_field_copy_range;
  vpf->apply_ic_func = vp_field_apply_ic;
  vpf->apply_bc_func = vp_field_apply_bc;
  vpf->limiter_func = vp_field_limiter;
  vpf->complete_update_func = vp_field_complete_update;
  vpf->calc_ext_em_func = vp_field_calc_ext_em;
  vpf->calc_app_current_func = vp_field_calc_app_current;
  vpf->calc_ext_pot_func = vp_field_calc_ext_pot;
  vpf->calc_energy_func = vp_field_calc_energy;
  vpf->write_func = vp_field_write;
  vpf->write_energy_func = vp_field_write_energy;
  vpf->from_file_func = vp_field_from_file;
  vpf->release_func = vp_field_release;

  return vpf;
}

// Vlasov-Poisson field update: solve for the potentials at the current time from
// the species densities. Elliptic solves (not part of the RK state vector);
// impose no CFL constraint of their own.
double
vp_field_update(
  gkyl_vlasov_app *app, double tcurr, const struct gkyl_array *fin[], const struct gkyl_array *emin,
  struct gkyl_array *emout
)
{
  vp_field_solve(app, app->field, fin);
  return DBL_MAX;
}

// Evaluate the external potentials and fields at the given time. Shared by the
// initial conditions and the restart: the potentials themselves are solved from
// the distribution wherever they are needed, so neither path solves them here.
static void
vp_field_calc_ext(gkyl_vlasov_app *app, struct vm_field *field, double tm)
{
  vp_field_calc_ext_pot(app, field, tm);
  vp_field_calc_ext_em(app, field, tm);
}

struct gkyl_app_restart_status
vp_field_from_file(gkyl_vlasov_app *app, struct vm_field *field, const char *fname)
{
  vp_field_calc_ext(app, field, app->tcurr);
  return (struct gkyl_app_restart_status){
    .io_status = GKYL_ARRAY_RIO_SUCCESS,
    .frame = 0,
    .stime = 0.0,
  };
}

// Vlasov-Poisson stage operations that are no-ops: the potentials are re-solved
// each stage and there are no EM RK stages, applied currents, EM BCs, or EM
// limiting for the potentials.
void
vp_field_combine(
  gkyl_vlasov_app *app, struct vm_field *field, struct gkyl_array *out, double c1,
  const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2
)
{
}

void
vp_field_copy_range(
  gkyl_vlasov_app *app, struct vm_field *field, struct gkyl_array *out, const struct gkyl_array *inp
)
{
}

void
vp_field_apply_bc(gkyl_vlasov_app *app, const struct vm_field *field, struct gkyl_array *em)
{
}

void
vp_field_limiter(gkyl_vlasov_app *app, struct vm_field *field, struct gkyl_array *em)
{
}

void
vp_field_complete_update(
  gkyl_vlasov_app *app, double dt, const struct gkyl_array *fin[],
  const struct gkyl_array *fluidin[], const struct gkyl_array *emin, struct gkyl_array *emout
)
{
}

void
vp_field_calc_ext_pot(gkyl_vlasov_app *app, struct vm_field *field, double tm)
{
  if (field->has_ext_pot) {
    gkyl_eval_on_nodes_advance(field->ext_pot_proj, tm, &app->local, field->ext_pot_host);
    if (app->use_gpu) {
      // Note: ext_pot_host is same as ext_pot when not on GPUs.
      gkyl_array_copy(field->ext_pot, field->ext_pot_host);
    }
  }
}

void
vp_field_calc_ext_em(gkyl_vlasov_app *app, struct vm_field *field, double tm)
{
  if (field->has_ext_em) {
    gkyl_proj_on_basis_advance(field->ext_em_proj, tm, &app->local_ext, field->ext_em_host);
    if (app->use_gpu) {
      // Note: ext_em_host is same as ext_em when not on GPUs.
      gkyl_array_copy(field->ext_em, field->ext_em_host);
    }
  }
}

void
vp_field_calc_app_current(gkyl_vlasov_app *app, struct vm_field *field, double tm)
{
  // No applied currents in Vlasov-Poisson.
}

void
vp_field_solve(gkyl_vlasov_app *app, struct vm_field *field, const struct gkyl_array *fin[])
{
  struct timespec wst = gkyl_wall_clock();

  // Accumulate the sources: each species owns its explicit contribution and
  // deposits its density, weighted by its coupling to each potential, onto that
  // potential's source (species without a Poisson coupling are a no-op). fin[]
  // is indexed over the overall species count.
  for (int p = 0; p < field->num_pots; ++p) {
    gkyl_array_clear(field->pots[p].rho, 0.0);
  }
  int num_species = app->num_species;
  for (int i = 0; i < num_species; ++i) {
    vlasov_species_accumulate_field_coupling(app, &app->species[i], fin[i], 0, 0);
  }

  for (int p = 0; p < field->num_pots; ++p) {
    vp_potential_solve(app, field, &field->pots[p]);
  }

  app->stat.field_rhs_tm += gkyl_time_diff_now_sec(wst);
}

void
vp_field_apply_ic(
  gkyl_vlasov_app *app, struct vm_field *field, const struct gkyl_array *fin[], double t0
)
{
  vp_field_calc_ext(app, field, t0);
}

void
vp_field_write(gkyl_vlasov_app *app, double tm, int frame, const struct gkyl_array *fin[])
{
  // The potentials are solved on demand so the written potentials are at time tm.
  vp_field_solve(app, app->field, fin);

  struct timespec wst = gkyl_wall_clock();

  struct gkyl_msgpack_data *mt = vlasov_array_meta_new((struct vlasov_output_meta){
    .frame = frame,
    .stime = tm,
    .poly_order = app->poly_order,
    .basis_type = app->basis.id,
  });

  for (int p = 0; p < app->field->num_pots; ++p) {
    vp_potential_write(app, &app->field->pots[p], mt, frame);
  }

  if (app->field->has_ext_em) {
    // Only write out external fields at t=0 or if they are time-dependent.
    if (frame == 0 || app->field->ext_em_evolve) {
      const char *fmt_ext_em = "%s-field_ext_em_%d.gkyl";
      int sz_ext_em = gkyl_calc_strlen(fmt_ext_em, app->name, frame);
      char fileNm_ext_em[sz_ext_em + 1]; // Ensures no buffer overflow.
      snprintf(fileNm_ext_em, sizeof fileNm_ext_em, fmt_ext_em, app->name, frame);

      // External EM field computed with project on basis, so just use host copy.
      vp_field_calc_ext_em(app, app->field, tm);

      gkyl_comm_array_write(
        app->comm, &app->grid, &app->local, mt, app->field->ext_em_host, fileNm_ext_em
      );
    }
  }
  if (app->field->has_ext_pot) {
    if (frame == 0 || app->field->ext_pot_evolve) {
      const char *fmt_ext_pot = "%s-field_ext_pot_%d.gkyl";
      int sz_ext_pot = gkyl_calc_strlen(fmt_ext_pot, app->name, frame);
      char fileNm_ext_pot[sz_ext_pot + 1]; // ensures no buffer overflow
      snprintf(fileNm_ext_pot, sizeof fileNm_ext_pot, fmt_ext_pot, app->name, frame);

      // External potentials computed with project on basis, so just use host copy.
      vp_field_calc_ext_pot(app, app->field, tm);
      gkyl_comm_array_write(
        app->comm, &app->grid, &app->local, mt, app->field->ext_pot_host, fileNm_ext_pot
      );
    }
  }

  vlasov_array_meta_release(mt);

  app->stat.field_io_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_field_io += 1;
}

void
vp_field_calc_energy(
  gkyl_vlasov_app *app, double tm, struct vm_field *field, const struct gkyl_array *fin[]
)
{
  // The potentials are solved on demand so the energies are at time tm.
  vp_field_solve(app, field, fin);

  struct timespec wst = gkyl_wall_clock();

  for (int p = 0; p < field->num_pots; ++p) {
    vp_potential_calc_energy(app, field, &field->pots[p], tm);
  }

  app->stat.field_diag_calc_tm += gkyl_time_diff_now_sec(wst);
}

void
vp_field_write_energy(gkyl_vlasov_app *app)
{
  struct timespec wst = gkyl_wall_clock();

  // The energy diagnostics of all potentials are written at the same cadence, so
  // one flag says whether their files are created or appended to.
  for (int p = 0; p < app->field->num_pots; ++p) {
    vp_potential_write_energy(app, &app->field->pots[p], app->field->is_first_energy_write_call);
  }
  app->field->is_first_energy_write_call = false;

  app->stat.n_field_diag_io += 1;
  app->stat.field_diag_io_tm += gkyl_time_diff_now_sec(wst);
}

void
vp_field_release(const gkyl_vlasov_app *app, struct vm_field *vpf)
{
  // Release resources for Vlasov-Poisson field.

  for (int p = 0; p < vpf->num_pots; ++p) {
    vp_potential_release(app, &vpf->pots[p]);
  }

  gkyl_array_integrate_release(vpf->calc_grad_sq);
  gkyl_array_release(vpf->grad_sq_wgt);

  gkyl_array_release(vpf->ext_pot);
  if (vpf->has_ext_pot) {
    gkyl_array_release(vpf->ext_pot_host);
    gkyl_eval_on_nodes_release(vpf->ext_pot_proj);
  }

  gkyl_array_release(vpf->ext_em);
  if (vpf->has_ext_em) {
    gkyl_array_release(vpf->ext_em_host);
    gkyl_proj_on_basis_release(vpf->ext_em_proj);
  }

  gkyl_free(vpf);
}
