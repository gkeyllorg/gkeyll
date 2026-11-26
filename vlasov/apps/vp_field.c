#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_bc_basic.h>
#include <gkyl_dg_differentiate.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_priv.h>

#include <assert.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <time.h>

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

// Square of the plasma frequency per unit number density, q^2/(epsilon0 m).
static double
vp_field_es_rate_sq_weight(const struct vp_potential *pot, double charge, double mass)
{
  return charge * charge / (pot->coupling * mass);
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

// Square of the Jeans rate per unit number density, alpha_g m.
static double
vp_field_grav_rate_sq_weight(const struct vp_potential *pot, double charge, double mass)
{
  return pot->coupling * mass;
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
// scale * J / J_xi^2, with J_xi the position-map Jacobian in direction i
// and J = prod_i J_xi the total conf Jacobian.
static void
vp_field_set_metric_diag(const struct gkyl_vlasov_app *app, struct gkyl_array *arr_ho, double scale)
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
      arr_d[diag * nb] = scale * Jtot / (Jxi * Jxi) * dg0;
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
    vp_field_set_metric_diag(app, eps_local_ho, scale);
    if (app->use_gpu) {
      gkyl_array_copy(eps_local, eps_local_ho);
    }
    gkyl_comm_array_allgather(app->comm, &app->local, &app->global, eps_local, epsilon);
    gkyl_array_release(eps_local);
    gkyl_array_release(eps_local_ho);
  }
  return epsilon;
}

// Per-cell constant DG scalar 'scale' times the conf Jacobian J (1 on an
// identity position map) on the local extended range, or allgathered to the
// global extended range the FEM solver reads. Used for the Helmholtz mass term
// and for the weight of the int phi^2 energy term, whose computational-space
// integrals need J to equal the physical ones.
static struct gkyl_array *
vp_field_new_jacobian_scalar(struct gkyl_vlasov_app *app, double scale, bool is_global)
{
  int nb = app->basis.num_basis;
  double dg0 = pow(sqrt(2.0), app->cdim); // 0th DG coeff representing a constant value.
  const struct gkyl_range *ext = is_global ? &app->global_ext : &app->local_ext;

  struct gkyl_array *arr = mkarr(app->use_gpu, nb, ext->volume);
  gkyl_array_clear(arr, 0.0);
  if (app->pos_map->is_identity) {
    gkyl_array_shiftc(arr, scale * dg0, 0);
    return arr;
  }

  struct gkyl_array *arr_local = is_global ? mkarr(app->use_gpu, nb, app->local_ext.volume) :
                                             gkyl_array_acquire(arr);
  struct gkyl_array *arr_local_ho = app->use_gpu ? mkarr(false, nb, app->local_ext.volume) :
                                                   gkyl_array_acquire(arr_local);
  gkyl_array_clear(arr_local_ho, 0.0);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  while (gkyl_range_iter_next(&iter)) {
    long cidx = gkyl_range_idx(&app->local, iter.idx);
    const double *jacob_pos_gauss = gkyl_array_cfetch(app->pos_map->jacob_pos_gauss_host, cidx);
    double *arr_d = gkyl_array_fetch(arr_local_ho, cidx);
    arr_d[0] = scale * jacob_pos_gauss[0] * dg0;
  }
  if (app->use_gpu) {
    gkyl_array_copy(arr_local, arr_local_ho);
  }
  if (is_global) {
    gkyl_comm_array_allgather(app->comm, &app->local, &app->global, arr_local, arr);
  }
  gkyl_array_release(arr_local);
  gkyl_array_release(arr_local_ho);
  return arr;
}

// --- One potential: allocation, solve, diagnostics, release -------------------

// Neutralizing background of the source (see neutralize_func): none, or on a position map J times
// the physical mean of the source.
static void
vp_potential_neutralize_none(gkyl_vlasov_app *app, struct vm_field *field, struct vp_potential *pot)
{
}

static void
vp_potential_neutralize_mapped(
  gkyl_vlasov_app *app, struct vm_field *field, struct vp_potential *pot
)
{
  gkyl_dg_calc_average_range(&app->basis, 0, pot->rho_avg, 0, pot->rho, app->local);
  gkyl_array_reduce_range(pot->rho_red, pot->rho_avg, GKYL_SUM, &app->local);
  double red_ho[1];
  if (app->use_gpu) {
    gkyl_cu_memcpy(red_ho, pot->rho_red, sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    red_ho[0] = pot->rho_red[0];
  }
  double red_global[1];
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, red_ho, red_global);
  // The cells have equal computational volumes: the physical mean is the sum of the cell
  // averages of J rho over the sum of J.
  double rho_mean = red_global[0] / pot->jac_sum;
  gkyl_array_accumulate_range(pot->rho, -rho_mean, pot->jac_dg, &app->local);
}

// Allocate the source, potential, permittivity (scalar eps_scale), screening
// (mu_sq > 0 adds the Helmholtz term eps_scale * mu_sq * phi), Poisson solver
// and energy diagnostic of a potential; its names and coupling laws are set by
// the caller. All potentials share the field's Poisson boundary conditions.
static void
vp_potential_init(
  struct gkyl_vlasov_app *app, struct vm_field *field, struct vp_potential *pot, double eps_scale,
  double mu_sq
)
{
  int nb = app->basis.num_basis;

  pot->rho = mkarr(app->use_gpu, nb, app->local_ext.volume);
  pot->rho_global = mkarr(app->use_gpu, nb, app->global_ext.volume);

  pot->phi = mkarr(app->use_gpu, nb, app->local_ext.volume);
  pot->phi_global = mkarr(app->use_gpu, nb, app->global_ext.volume);
  pot->phi_host = app->use_gpu ? mkarr(false, pot->phi->ncomp, pot->phi->size) :
                                 gkyl_array_acquire(pot->phi);

  // The FEM solver solves - nabla . (epsilon * nabla phi) - kSq * phi = rho, so
  // the screening term + eps_scale * mu_sq * phi is kSq = - eps_scale * mu_sq.
  pot->mu_sq = mu_sq;
  pot->ksq = NULL;
  pot->sq_wgt = NULL;
  pot->calc_sq = NULL;
  if (mu_sq > 0.0) {
    pot->ksq = vp_field_new_jacobian_scalar(app, -eps_scale * mu_sq, true);
    pot->sq_wgt = vp_field_new_jacobian_scalar(app, 1.0, false);
    pot->calc_sq = gkyl_array_integrate_new(
      &app->grid, &app->basis, 1, GKYL_ARRAY_INTEGRATE_OP_SQ_WEIGHTED, app->use_gpu
    );
  }

  // The neutralizing background of an unscreened fully periodic solve on a position map.
  bool is_periodic = true;
  for (int d = 0; d < app->cdim; ++d) {
    is_periodic = is_periodic && field->info.poisson_bcs.lo_type[d] == GKYL_POISSON_PERIODIC;
  }
  pot->neutralize_func = vp_potential_neutralize_none;
  pot->jac_dg = 0;
  pot->rho_avg = 0;
  pot->rho_red = 0;
  pot->jac_sum = 0.0;
  if (is_periodic && !(mu_sq > 0.0) && !app->pos_map->is_identity) {
    pot->jac_dg = vp_field_new_jacobian_scalar(app, 1.0, false);
    pot->rho_avg = mkarr(app->use_gpu, 1, app->local_ext.volume);
    pot->rho_red = app->use_gpu ? gkyl_cu_malloc(sizeof(double)) : gkyl_malloc(sizeof(double));
    double jac_sum_local = 0.0;
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &app->local);
    while (gkyl_range_iter_next(&iter)) {
      long cidx = gkyl_range_idx(&app->local, iter.idx);
      const double *jacob_pos_gauss = gkyl_array_cfetch(app->pos_map->jacob_pos_gauss_host, cidx);
      jac_sum_local += jacob_pos_gauss[0];
    }
    gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &jac_sum_local, &pot->jac_sum);
    pot->neutralize_func = vp_potential_neutralize_mapped;
  }

  pot->epsilon = vp_field_new_epsilon(app, eps_scale);
  pot->fem_poisson = gkyl_fem_poisson_new(
    &app->global, &app->grid, app->basis, &field->info.poisson_bcs, NULL, pot->epsilon, pot->ksq,
    app->pos_map->is_identity, app->use_gpu
  );

  // Reductions of [int |grad phi|^2, int phi^2]; the second stays 0 if unscreened.
  double zeros[2] = {0.0, 0.0};
  if (app->use_gpu) {
    pot->energy_red = gkyl_cu_malloc(sizeof(double[2]));
    pot->energy_red_global = gkyl_cu_malloc(sizeof(double[2]));
    gkyl_cu_memcpy(pot->energy_red, zeros, sizeof(double[2]), GKYL_CU_MEMCPY_H2D);
    gkyl_cu_memcpy(pot->energy_red_global, zeros, sizeof(double[2]), GKYL_CU_MEMCPY_H2D);
  } else {
    pot->energy_red = gkyl_malloc(sizeof(double[2]));
    pot->energy_red_global = gkyl_malloc(sizeof(double[2]));
    memcpy(pot->energy_red, zeros, sizeof(double[2]));
    memcpy(pot->energy_red_global, zeros, sizeof(double[2]));
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
  if (pot->jac_dg) {
    gkyl_array_release(pot->jac_dg);
    gkyl_array_release(pot->rho_avg);
    if (app->use_gpu) {
      gkyl_cu_free(pot->rho_red);
    } else {
      gkyl_free(pot->rho_red);
    }
  }
  if (pot->calc_sq) {
    gkyl_array_integrate_release(pot->calc_sq);
    gkyl_array_release(pot->sq_wgt);
    gkyl_array_release(pot->ksq);
  }

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
  pot->neutralize_func(app, field, pot);
  gkyl_comm_array_allgather(app->comm, &app->local, &app->global, pot->rho, pot->rho_global);

  gkyl_fem_poisson_set_rhs(pot->fem_poisson, pot->rho_global, NULL);
  gkyl_fem_poisson_solve(pot->fem_poisson, pot->phi_global);

  gkyl_array_copy_range_to_range(pot->phi, pot->phi_global, &app->local, &field->global_sub_range);
}

// Append the integrated |grad phi|^2 + mu_sq phi^2 (scaled by the potential's
// energy factor) at time tm to the potential's energy diagnostic.
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
  if (pot->calc_sq) {
    gkyl_array_integrate_advance(
      pot->calc_sq, pot->phi, 1.0, pot->sq_wgt, &app->local, &app->local, pot->energy_red + 1
    );
  }

  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 2, pot->energy_red, pot->energy_red_global);

  double red_global[2] = {0.0, 0.0};
  if (app->use_gpu) {
    gkyl_cu_memcpy(red_global, pot->energy_red_global, sizeof(double[2]), GKYL_CU_MEMCPY_D2H);
  } else {
    red_global[0] = pot->energy_red_global[0];
    red_global[1] = pot->energy_red_global[1];
  }
  double energy_global[1] = {pot->energy_fac * (red_global[0] + pot->mu_sq * red_global[1])};

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
    // Electrostatics: - nabla . (epsilon0 * nabla phi) + epsilon0 * mu_sq * phi =
    // sum_s q_s n_s, screened by a massive dark photon when mu_sq > 0.
    struct vp_potential *pot = &vpf->pots[vpf->num_pots++];
    strcpy(pot->frame_name, "field");
    strcpy(pot->energy_name, "field-energy");
    pot->coupling = vpf->info.epsilon0;
    pot->energy_fac = 1.0; // Diagnostic is int |grad phi|^2 + mu_sq phi^2 (no epsilon0/2).
    pot->src_weight = vp_field_es_src_weight;
    pot->force_weight = vp_field_es_force_weight;
    pot->rate_sq_weight = vp_field_es_rate_sq_weight;
    pot->is_electric = true;
    vp_potential_init(app, vpf, pot, vpf->info.epsilon0, vpf->info.mu_sq);
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
    pot->rate_sq_weight = vp_field_grav_rate_sq_weight;
    pot->is_electric = false;
    vp_potential_init(app, vpf, pot, -1.0, 0.0);
  }

  // Initialize external potentials.
  vpf->ext_pot = mkarr(app->use_gpu, 4 * app->basis.num_basis, app->local_ext.volume);
  gkyl_array_clear(vpf->ext_pot, 0.0);
  vpf->has_ext_pot = false;
  vpf->ext_pot_evolve = false;
  // setup external potentials
  if (vpf->info.external_potentials) {
    vpf->has_ext_pot = true;
    if (vpf->info.external_potentials_evolve) {
      vpf->ext_pot_evolve = vpf->info.external_potentials_evolve;
    }

    vpf->ext_pot_host = app->use_gpu ? mkarr(false, vpf->ext_pot->ncomp, vpf->ext_pot->size) :
                                       gkyl_array_acquire(vpf->ext_pot);
    // Evaluate at nodes, on the physical coordinates of the (possibly mapped) conf mesh.
    // 4 components: phi and the three components of A.
    vpf->ext_pot_proj = gkyl_eval_on_nodes_inew(&(struct gkyl_eval_on_nodes_inp){
      .grid = &app->grid,
      .basis = &app->basis,
      .num_ret_vals = 4,
      .eval = vpf->info.external_potentials,
      .ctx = vpf->info.external_potentials_ctx,
      .c2p_func = vm_pos_c2p,
      .c2p_func_ctx = app->pos_map,
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
      .c2p_func = vm_pos_c2p,
      .c2p_func_ctx = app->pos_map,
    });
  }

  // No applied currents in Vlasov-Poisson.
  vpf->has_app_current = vpf->app_current_evolve = false;

  // Integrated |grad phi|^2 energy diagnostic shared by the potentials. For an
  // identity position map this is the plain |grad phi|^2 operator (GRAD_SQ) with
  // weight 1. For a mapped conf mesh the physical
  //   int |grad_x phi|^2 dx = int sum_i (J/J_xi^2)(d_xi phi)^2 dxi,
  // i.e. the diagonal metric J/J_xi^2 (J = prod_i J_xi, per-cell constant; the
  // same metric as the Poisson permittivity with unit scale) folded into the
  // full-gradient weighted operator (EPS_GRAD_SQ), which integrates over the
  // computational cells. The weight uses the same symmetric upper-triangular
  // layout as the Poisson permittivity tensor.
  bool wgt_const = app->pos_map->is_identity;
  int wgtnum = vp_field_tensor_num_comp(app, wgt_const);
  vpf->grad_sq_wgt = mkarr(app->use_gpu, wgtnum * nb, app->local_ext.volume);
  gkyl_array_clear(vpf->grad_sq_wgt, 0.0);
  if (wgt_const) {
    gkyl_array_shiftc(vpf->grad_sq_wgt, dg0, 0); // Sets the weight to 1 (unused by GRAD_SQ).
  } else {
    struct gkyl_array *wgt_ho = app->use_gpu ? mkarr(false, wgtnum * nb, app->local_ext.volume) :
                                               gkyl_array_acquire(vpf->grad_sq_wgt);
    vp_field_set_metric_diag(app, wgt_ho, 1.0);
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

  // The forces of the kick of the staggered operator split are allocated by
  // vp_field_set_source_split.
  vpf->em_pot = 0;
  vpf->accel_grav = 0;
  vpf->phi_sum = 0;

  // A PKPM species takes b = B/|B|, its gradient, d_t B/|B| and the total field from the field. The
  // magnetic field of a Poisson field is the external one, static (d_t B/|B| = 0, and b and its
  // gradient are formed once); E in the total field is that of the kick.
  vpf->has_bvar = false;
  for (int i = 0; i < vm->num_species; ++i) {
    vpf->has_bvar = vpf->has_bvar || vm->species[i].type == GKYL_SPECIES_PKPM;
  }
  if (vpf->has_bvar) {
    if (!vpf->info.ext_em) {
      gkyl_exit(
        "vp_field: a PKPM species with the Poisson field needs an external magnetic field (ext_em), which sets b."
      );
    }
    if (vpf->info.ext_em_evolve) {
      gkyl_exit(
        "vp_field: a PKPM species with the Poisson field needs a static external field (ext_em_evolve = false)."
      );
    }
    long num_cells = app->local_ext.volume;
    vpf->tot_em = mkarr(app->use_gpu, 8 * nb, num_cells);
    gkyl_array_clear(vpf->tot_em, 0.0);
    vpf->bvar = mkarr(app->use_gpu, 3 * nb, num_cells);
    vpf->grad_b = mkarr(app->use_gpu, cdim * 3 * nb, num_cells);
    gkyl_array_clear(vpf->grad_b, 0.0);
    vpf->dt_em_no_J = mkarr(app->use_gpu, 8 * nb, num_cells);
    gkyl_array_clear(vpf->dt_em_no_J, 0.0);
    vpf->dtB_over_magB = mkarr(app->use_gpu, 3 * nb, num_cells);
    gkyl_array_clear(vpf->dtB_over_magB, 0.0);
    vpf->calc_bvar = gkyl_dg_calc_pkpm_bvar_inew(&(struct gkyl_dg_calc_pkpm_bvar_inp){
      .conf_grid = &app->grid,
      .conf_basis = &app->basis,
      .pos_map = app->pos_map,
      .use_gpu = app->use_gpu,
    });
  }

  // Set the type-specific dispatch methods (Vlasov-Poisson). The potentials are
  // solved each stage (update_func), so the combine/copy/BC/current/limiter
  // stage operations are no-ops.
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
  vpf->calc_diag_vars_func = vp_field_calc_diag_vars;
  vpf->from_file_func = vp_field_from_file;
  vpf->release_func = vp_field_release;

  return vpf;
}

// Vlasov-Poisson field update: solve for the potentials at the current time from
// the species densities. Elliptic solves (not part of the RK state vector);
// impose no CFL constraint of their own.
double
vp_field_update(
  gkyl_vlasov_app *app, int stage, double tcurr, const struct gkyl_array *fin[],
  const struct gkyl_array *fluidin[], const struct gkyl_array *emin, struct gkyl_array *emout
)
{
  vp_field_solve(app, app->field, fin, fluidin);
  return DBL_MAX;
}

// Field update of a stage of the transport step of the staggered operator split without kinetic
// species: no species reads the potentials there (the fluid aspects take their forces in the kick).
static double
vp_field_update_split_none(
  gkyl_vlasov_app *app, int stage, double tcurr, const struct gkyl_array *fin[],
  const struct gkyl_array *fluidin[], const struct gkyl_array *emin, struct gkyl_array *emout
)
{
  return DBL_MAX;
}

// b and its gradient of the static external magnetic field (PKPM species), over the extended
// range: the boundary conditions of the species read b in the ghost cells. tot_em holds the
// external field and the E of the last kick.
static void
vp_field_calc_bvar(gkyl_vlasov_app *app, struct vm_field *field)
{
  gkyl_dg_calc_pkpm_bvar_advance(field->calc_bvar, &app->local_ext, field->tot_em, field->bvar);
  gkyl_dg_calc_pkpm_bvar_grad(field->calc_bvar, &app->local, field->bvar, field->grad_b);
}

// The total field of PKPM species, [E + E_ext, B_ext], with the E of the kick (zero before the
// split has formed one).
static void
vp_field_calc_tot_em(gkyl_vlasov_app *app, struct vm_field *field)
{
  gkyl_array_clear(field->tot_em, 0.0);
  if (field->em_pot) {
    gkyl_array_accumulate_range(field->tot_em, 1.0, field->em_pot, &app->local_ext);
  }
  gkyl_array_accumulate_range(field->tot_em, 1.0, field->ext_em, &app->local_ext);
}

// The physical gradient out = -grad(phi) of a continuous p=1 potential, exact in the p=1 basis: on a
// uniform mesh the gradient of the expansion, on a position map that gradient divided in each
// direction by the per-cell constant J_d of the C^0 linear map (the J_d of the metric of the
// Poisson solve, so the work of the kick and the change of the field energy stay consistent).
static void
vp_field_grad_identity(
  gkyl_vlasov_app *app, struct vm_field *field, const struct gkyl_array *phi, struct gkyl_array *out
)
{
  for (int d = 0; d < app->cdim; ++d) {
    gkyl_dg_differentiate_op_local_range(
      &field->basis_grad, d, 1, app->grid.dx[d], d, out, 0, phi, &app->local
    );
  }
  gkyl_array_scale_range(out, -1.0, &app->local);
}

static void
vp_field_grad_mapped(
  gkyl_vlasov_app *app, struct vm_field *field, const struct gkyl_array *phi, struct gkyl_array *out
)
{
  for (int d = 0; d < app->cdim; ++d) {
    gkyl_dg_differentiate_op_local_range(
      &field->basis_grad, d, 1, app->grid.dx[d], d, field->grad_buf, 0, phi, &app->local
    );
    gkyl_dg_mul_op_range(&app->basis, d, out, d, field->rjacob_dir, d, field->grad_buf, &app->local);
  }
  gkyl_array_scale_range(out, -1.0, &app->local);
}

// The parts of the potentials in the forces of the kick: phi of the electrostatic potential adds to
// the potential of the electric field, the gravitational potential gives the acceleration of every
// species.
static void
vp_potential_kick_electric(
  gkyl_vlasov_app *app, struct vm_field *field, const struct vp_potential *pot
)
{
  gkyl_array_accumulate_range(field->phi_sum, 1.0, pot->phi, &app->local);
}

static void
vp_potential_kick_grav(gkyl_vlasov_app *app, struct vm_field *field, const struct vp_potential *pot)
{
  field->grad_func(app, field, pot->phi, field->accel_grav);
}

static void
vp_field_kick_tot_em_none(gkyl_vlasov_app *app, struct vm_field *field)
{
}

void
vp_field_set_source_split(gkyl_vlasov_app *app, struct vm_field *field)
{
  if (app->poly_order != 1) {
    gkyl_exit(
      "vp_field: the coupling of fluid aspects to the Poisson field requires polynomial order 1."
    );
  }
  int nb = app->basis.num_basis;
  long num_cells = app->local_ext.volume;
  field->em_pot = mkarr(app->use_gpu, 8 * nb, num_cells);
  gkyl_array_clear(field->em_pot, 0.0);
  field->accel_grav = mkarr(app->use_gpu, 3 * nb, num_cells);
  gkyl_array_clear(field->accel_grav, 0.0);
  field->phi_sum = mkarr(app->use_gpu, nb, num_cells);
  gkyl_array_clear(field->phi_sum, 0.0);
  gkyl_cart_modal_serendip(&field->basis_grad, app->cdim, 1);
  field->grad_func = vp_field_grad_identity;
  field->rjacob_dir = 0;
  field->grad_buf = 0;
  if (!app->pos_map->is_identity) {
    // 1/J_d of each cell as a constant DG field per direction (built on the host from the map).
    int cdim = app->cdim;
    double dg0 = pow(sqrt(2.0), cdim); // 0th DG coefficient of a constant.
    int stride = app->basis.poly_order + 1; // per-direction block stride in jacob_pos
    field->rjacob_dir = mkarr(app->use_gpu, cdim * nb, num_cells);
    struct gkyl_array *rjacob_ho = app->use_gpu ? mkarr(false, cdim * nb, num_cells) :
                                                  gkyl_array_acquire(field->rjacob_dir);
    gkyl_array_clear(rjacob_ho, 0.0);
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &app->local);
    while (gkyl_range_iter_next(&iter)) {
      long cidx = gkyl_range_idx(&app->local, iter.idx);
      const double *jacob_pos = gkyl_array_cfetch(app->pos_map->jacob_pos_host, cidx);
      double *r_d = gkyl_array_fetch(rjacob_ho, cidx);
      for (int d = 0; d < cdim; ++d) {
        r_d[d * nb] = dg0 / jacob_pos[d * stride];
      }
    }
    if (app->use_gpu) {
      gkyl_array_copy(field->rjacob_dir, rjacob_ho);
    }
    gkyl_array_release(rjacob_ho);
    field->grad_buf = mkarr(app->use_gpu, cdim * nb, num_cells);
    gkyl_array_clear(field->grad_buf, 0.0);
    field->grad_func = vp_field_grad_mapped;
  }
  for (int p = 0; p < field->num_pots; ++p) {
    struct vp_potential *pot = &field->pots[p];
    pot->kick_func = pot->is_electric ? vp_potential_kick_electric : vp_potential_kick_grav;
  }
  field->kick_tot_em_func = field->has_bvar ? vp_field_calc_tot_em : vp_field_kick_tot_em_none;

  // The kinetic species of a hybrid run take the potentials in each stage of the transport step,
  // solved from the stage densities of all species; without them the stages solve nothing.
  bool has_kinetic = false;
  for (int i = 0; i < app->num_species; ++i) {
    has_kinetic = has_kinetic || app->species[i].type == GKYL_SPECIES_VLASOV;
  }
  field->update_func = has_kinetic ? vp_field_update : vp_field_update_split_none;
}

void
vp_field_calc_kick(gkyl_vlasov_app *app, struct vm_field *field)
{
  int num_species = app->num_species;
  struct gkyl_array *distf[num_species];
  struct gkyl_array *fluid[num_species];
  vlasov_species_gather_state(app, distf, fluid);
  const struct gkyl_array *fin[num_species];
  const struct gkyl_array *fluidin[num_species];
  for (int i = 0; i < num_species; ++i) {
    fin[i] = distf[i];
    fluidin[i] = fluid[i];
  }
  vp_field_solve(app, field, fin, fluidin);

  struct timespec wst = gkyl_wall_clock();

  // E = -grad(phi + phi_ext) (the external potential is zero when there is none) and the
  // gravitational acceleration -grad(phi_g).
  gkyl_array_clear(field->phi_sum, 0.0);
  gkyl_array_accumulate_offset(field->phi_sum, 1.0, field->ext_pot, 0);
  for (int p = 0; p < field->num_pots; ++p) {
    field->pots[p].kick_func(app, field, &field->pots[p]);
  }
  field->grad_func(app, field, field->phi_sum, field->em_pot);
  field->kick_tot_em_func(app, field);

  app->stat.field_rhs_tm += gkyl_time_diff_now_sec(wst);
}

// The field quantities of PKPM species for the output: the E of the present state in the total
// field (b and its gradient are static).
void
vp_field_calc_diag_vars(gkyl_vlasov_app *app, struct vm_field *field)
{
  if (field->em_pot) {
    vp_field_calc_kick(app, field);
  }
}

// Evaluate the external potentials and fields at time tm for the initial conditions.
static void
vp_field_calc_ext(gkyl_vlasov_app *app, struct vm_field *field, double tm)
{
  vp_field_calc_ext_pot(app, field, tm);
  vp_field_calc_ext_em(app, field, tm);
}

struct gkyl_app_restart_status
vp_field_from_file(gkyl_vlasov_app *app, struct vm_field *field, const char *fname)
{
  // The time-stepping loop only re-evaluates the external potentials and fields
  // when they evolve. Since a quantity that is not evolving keeps its initial
  // value, pass 0.0 for it to guard against static inputs that were specified
  // as being time-dependent.
  vp_field_calc_ext_pot(app, field, field->ext_pot_evolve ? app->tcurr : 0.0);
  vp_field_calc_ext_em(app, field, field->ext_em_evolve ? app->tcurr : 0.0);
  if (field->has_bvar) {
    // b of the external field (the species boundary conditions read it); the E of the total
    // field is formed from the species once they are read (vlasov_op_split_restart).
    vp_field_calc_tot_em(app, field);
    vp_field_calc_bvar(app, field);
  }
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
vp_field_solve(
  gkyl_vlasov_app *app, struct vm_field *field, const struct gkyl_array *fin[],
  const struct gkyl_array *fluidin[]
)
{
  struct timespec wst = gkyl_wall_clock();

  // Accumulate the source of each potential: the number density of each species
  // times its coupling weight (fin[] and fluidin[] are indexed by species).
  for (int p = 0; p < field->num_pots; ++p) {
    gkyl_array_clear(field->pots[p].rho, 0.0);
  }
  int num_species = app->num_species;
  for (int i = 0; i < num_species; ++i) {
    vlasov_species_accumulate_field_coupling(app, &app->species[i], fin[i], fluidin[i], 0);
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
  // The forces of the kick of the initial state (the species are initialized first): with PKPM
  // species the first transport step of a run reads the E of the total field.
  if (field->em_pot) {
    vp_field_calc_kick(app, field);
  } else if (field->has_bvar) {
    vp_field_calc_tot_em(app, field);
  }
  if (field->has_bvar) {
    vp_field_calc_bvar(app, field);
  }
}

// Solve for the potentials from the present state: the distributions fin[] and the fluid states
// of the species.
static void
vp_field_solve_state(gkyl_vlasov_app *app, struct vm_field *field, const struct gkyl_array *fin[])
{
  int num_species = app->num_species;
  struct gkyl_array *fluid[num_species];
  vlasov_species_gather_state(app, 0, fluid);
  const struct gkyl_array *fluidin[num_species];
  for (int i = 0; i < num_species; ++i) {
    fluidin[i] = fluid[i];
  }
  vp_field_solve(app, field, fin, fluidin);
}

void
vp_field_write(gkyl_vlasov_app *app, double tm, int frame, const struct gkyl_array *fin[])
{
  // The potentials are solved on demand so the written potentials are at time tm.
  vp_field_solve_state(app, app->field, fin);

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

      // External potentials are evaluated on nodes, so just use host copy.
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
  vp_field_solve_state(app, field, fin);

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

  if (vpf->em_pot) {
    gkyl_array_release(vpf->em_pot);
    gkyl_array_release(vpf->accel_grav);
    gkyl_array_release(vpf->phi_sum);
    if (vpf->rjacob_dir) {
      gkyl_array_release(vpf->rjacob_dir);
      gkyl_array_release(vpf->grad_buf);
    }
  }
  if (vpf->has_bvar) {
    gkyl_array_release(vpf->tot_em);
    gkyl_array_release(vpf->bvar);
    gkyl_array_release(vpf->grad_b);
    gkyl_array_release(vpf->dt_em_no_J);
    gkyl_array_release(vpf->dtB_over_magB);
    gkyl_dg_calc_pkpm_bvar_release(vpf->calc_bvar);
  }

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
