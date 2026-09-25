#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_bc_basic.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_priv.h>

#include <assert.h>
#include <float.h>
#include <math.h>
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

struct vm_field *
vp_field_new(struct gkyl_vm *vm, struct gkyl_vlasov_app *app)
{
  // Initialize field object.
  struct vm_field *vpf = gkyl_malloc(sizeof(struct vm_field));

  vpf->info = vm->field;

  // Allocate arrays for charge density.
  vpf->rho_c = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  vpf->rho_c_global = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);

  // Allocate arrays for electrostatic potential.
  vpf->phi = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  vpf->phi_global = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);

  // Host potential for  I/O.
  vpf->phi_host = app->use_gpu ? mkarr(false, vpf->phi->ncomp, vpf->phi->size) :
                                 gkyl_array_acquire(vpf->phi);

  // Create global subrange we'll copy the field solver solution from (into local).
  int intersect = gkyl_sub_range_intersect(&vpf->global_sub_range, &app->global, &app->local);

  // Set the permittivity in the Poisson equation. On a uniform grid (identity
  // position map) this is the constant scalar epsilon0. On a non-uniform conf
  // mesh it is the diagonal metric permittivity tensor
  //   eps^{ii} = epsilon0 * J / J_xi^2   (off-diagonals 0, diagonal position map)
  // with J = prod_i J_xi the total conf Jacobian (per cell constant). This makes
  // the weak Poisson operator int eps^{ij} d_i(phi) d_j(psi) dxi equal the
  // physical int epsilon0 grad(phi).grad(psi) dx; the RHS rho_c already carries J.
  int cdim = app->cdim, nb = app->basis.num_basis;
  double dg0 = pow(sqrt(2.0), cdim); // 0th DG coeff representing a constant value.
  bool eps_const = app->pos_map->is_identity;
  // Symmetric permittivity tensor component count: 1x->1, 2x->3, 3x->6.
  int epsnum = eps_const ? 1 : cdim + (int)ceil((pow(3.0, cdim - 1) - cdim) / 2.0);

  vpf->epsilon = mkarr(app->use_gpu, epsnum * nb, app->global_ext.volume);
  gkyl_array_clear(vpf->epsilon, 0.0);
  if (eps_const) {
    gkyl_array_shiftc(vpf->epsilon, vpf->info.epsilon0 * dg0, 0);
  } else {
    // Build the local diagonal tensor from the position map, then allgather to
    // the global array the FEM solver reads (mirrors the rho_c allgather).
    struct gkyl_array *eps_local = mkarr(app->use_gpu, epsnum * nb, app->local_ext.volume);
    struct gkyl_array *eps_local_ho = app->use_gpu ?
                                        mkarr(false, epsnum * nb, app->local_ext.volume) :
                                        gkyl_array_acquire(eps_local);
    gkyl_array_clear(eps_local_ho, 0.0);
    int stride = app->basis.poly_order + 1; // per-direction block stride in jacob_pos.
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &app->local);
    while (gkyl_range_iter_next(&iter)) {
      long cidx = gkyl_range_idx(&app->local, iter.idx);
      const double *jacob_pos = gkyl_array_cfetch(app->pos_map->jacob_pos_host, cidx);
      const double *jacob_pos_gauss = gkyl_array_cfetch(app->pos_map->jacob_pos_gauss_host, cidx);
      double Jtot = jacob_pos_gauss[0];
      double *eps_d = gkyl_array_fetch(eps_local_ho, cidx);
      for (int i = 0; i < cdim; ++i) {
        double Jxi = jacob_pos[i * stride];
        int diag =
          i * cdim - (i * (i - 1)) / 2; // index of (i,i) in symmetric upper-triangular storage.
        eps_d[diag * nb] = vpf->info.epsilon0 * Jtot / (Jxi * Jxi) * dg0;
      }
    }
    if (app->use_gpu) {
      gkyl_array_copy(eps_local, eps_local_ho);
    }
    gkyl_comm_array_allgather(app->comm, &app->local, &app->global, eps_local, vpf->epsilon);
    gkyl_array_release(eps_local);
    gkyl_array_release(eps_local_ho);
  }

  // Create Poisson solver (variable permittivity when the conf mesh is mapped).
  vpf->fem_poisson = gkyl_fem_poisson_new(
    &app->global, &app->grid, app->basis, &vpf->info.poisson_bcs, NULL, vpf->epsilon, NULL,
    eps_const, app->use_gpu
  );

  vpf->field_id = GKYL_FIELD_PHI;

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

  if (app->use_gpu) {
    vpf->es_energy_red = gkyl_cu_malloc(sizeof(double[1]));
    vpf->es_energy_red_global = gkyl_cu_malloc(sizeof(double[1]));
  } else {
    vpf->es_energy_red = gkyl_malloc(sizeof(double[1]));
    vpf->es_energy_red_global = gkyl_malloc(sizeof(double[1]));
  }

  vpf->integ_energy = gkyl_dynvec_new(GKYL_DOUBLE, 1);

  // Electrostatic field-energy factor. For an identity position map this is the
  // scalar es_energy_fac=1 with the plain |grad phi|^2 operator (GRAD_SQ), exactly
  // as before. For a mapped conf mesh the physical field energy is
  //   int epsilon0/2 |grad_x phi|^2 dx = int sum_i (J/J_xi^2)(d_xi phi)^2 dxi,
  // i.e. a diagonal metric eps_ii = J/J_xi^2 (J = prod_i J_xi, per-cell constant;
  // no epsilon0 here, matching the existing |grad phi|^2 diagnostic convention)
  // folded into the full-gradient weighted operator (EPS_GRAD_SQ). The weight uses
  // the same symmetric upper-triangular layout as the Poisson permittivity tensor.
  bool es_fac_const = app->pos_map->is_identity;
  int es_epsnum = es_fac_const ? 1 : cdim + (int)ceil((pow(3.0, cdim - 1) - cdim) / 2.0);
  vpf->es_energy_fac = mkarr(app->use_gpu, es_epsnum * nb, app->local_ext.volume);
  gkyl_array_clear(vpf->es_energy_fac, 0.0);
  if (es_fac_const) {
    gkyl_array_shiftc(vpf->es_energy_fac, dg0, 0); // Sets es_energy_fac=1 (unused by GRAD_SQ).
  } else {
    struct gkyl_array *esfac_ho = app->use_gpu ?
                                    mkarr(false, es_epsnum * nb, app->local_ext.volume) :
                                    gkyl_array_acquire(vpf->es_energy_fac);
    gkyl_array_clear(esfac_ho, 0.0);
    int stride = app->basis.poly_order + 1; // per-direction block stride in jacob_pos.
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &app->local);
    while (gkyl_range_iter_next(&iter)) {
      long cidx = gkyl_range_idx(&app->local, iter.idx);
      const double *jacob_pos = gkyl_array_cfetch(app->pos_map->jacob_pos_host, cidx);
      const double *jacob_pos_gauss = gkyl_array_cfetch(app->pos_map->jacob_pos_gauss_host, cidx);
      double Jtot = jacob_pos_gauss[0];
      double *esfac_d = gkyl_array_fetch(esfac_ho, cidx);
      for (int i = 0; i < cdim; ++i) {
        double Jxi = jacob_pos[i * stride];
        int diag = i * cdim - (i * (i - 1)) / 2; // (i,i) in symmetric upper-triangular storage.
        // Diagonal metric weight = J^2/J_xi^2: one factor of J = sqrt(g) folds the
        // physical volume into the metric (g^ii = 1/J_xi^2), and a second factor of
        // J corrects the integrate operator's (computational) cell-volume factor to
        // the physical cell volume J*dxi_comp. With this weight EPS_GRAD_SQ matches
        // the uniform-grid GRAD_SQ diagnostic exactly for a constant-Jacobian map
        // (verified for all cdim) and reduces to GRAD_SQ for the identity map.
        esfac_d[diag * nb] = (Jtot * Jtot) / (Jxi * Jxi) * dg0;
      }
    }
    if (app->use_gpu) {
      gkyl_array_copy(vpf->es_energy_fac, esfac_ho);
    }
    gkyl_array_release(esfac_ho);
  }

  vpf->calc_es_energy = gkyl_array_integrate_new(
    &app->grid, &app->basis, 1,
    es_fac_const ? GKYL_ARRAY_INTEGRATE_OP_GRAD_SQ : GKYL_ARRAY_INTEGRATE_OP_EPS_GRAD_SQ,
    app->use_gpu
  );
  vpf->is_first_energy_write_call = true;

  // Set the type-specific dispatch methods (Vlasov-Poisson). The potential is
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

// Vlasov-Poisson field update: solve for the potential at the current time from
// the charge density. Elliptic solve (not part of the RK state vector); imposes
// no CFL constraint of its own.
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
// initial conditions and the restart: the potential itself is solved from the
// distribution wherever it is needed, so neither path solves it here.
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

// Vlasov-Poisson stage operations that are no-ops: the potential is re-solved
// each stage and there are no EM RK stages, applied currents, EM BCs, or EM
// limiting for the electrostatic potential.
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

  // Accumulate the charge density: each species owns its explicit contribution
  // (kinetic species accumulate q*m0; species without a Poisson coupling are a
  // no-op). fin[] is indexed over the overall species count.
  gkyl_array_clear(field->rho_c, 0.0);
  int num_species = app->num_species;
  for (int i = 0; i < num_species; ++i) {
    vlasov_species_accumulate_field_coupling(app, &app->species[i], fin[i], 0, field->rho_c);
  }

  // Gather charge density into global array.
  gkyl_comm_array_allgather(app->comm, &app->local, &app->global, field->rho_c, field->rho_c_global);

  // Solve the Poisson problem.
  gkyl_fem_poisson_set_rhs(field->fem_poisson, field->rho_c_global, NULL);
  gkyl_fem_poisson_solve(field->fem_poisson, field->phi_global);

  // Copy the portion of global potential corresponding to this MPI pcross to the local potential.
  gkyl_array_copy_range_to_range(
    field->phi, field->phi_global, &app->local, &field->global_sub_range
  );

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
  // The potential is solved on demand so the written phi is at time tm.
  vp_field_solve(app, app->field, fin);

  struct timespec wst = gkyl_wall_clock();

  struct gkyl_msgpack_data *mt = vlasov_array_meta_new((struct vlasov_output_meta){
    .frame = frame,
    .stime = tm,
    .poly_order = app->poly_order,
    .basis_type = app->basis.id,
  });

  const char *fmt = "%s-field_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, frame);
  char fileNm[sz + 1]; // Ensures no buffer overflow.
  snprintf(fileNm, sizeof fileNm, fmt, app->name, frame);

  // Copy data from device to host before writing it out.
  if (app->use_gpu) {
    gkyl_array_copy(app->field->phi_host, app->field->phi);
  }
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, app->field->phi_host, fileNm);

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
  // The potential is solved on demand so the energy is at time tm.
  vp_field_solve(app, field, fin);

  struct timespec wst = gkyl_wall_clock();

  gkyl_array_integrate_advance(
    field->calc_es_energy, field->phi, app->grid.cellVolume, field->es_energy_fac, &app->local,
    &app->local, field->es_energy_red
  );

  gkyl_comm_allreduce(
    app->comm, GKYL_DOUBLE, GKYL_SUM, 1, field->es_energy_red, field->es_energy_red_global
  );

  double energy_global[1] = {0.0};
  if (app->use_gpu) {
    gkyl_cu_memcpy(
      energy_global, field->es_energy_red_global, sizeof(double[1]), GKYL_CU_MEMCPY_D2H
    );
  } else {
    energy_global[0] = field->es_energy_red_global[0];
  }

  gkyl_dynvec_append(field->integ_energy, tm, energy_global);

  app->stat.field_diag_calc_tm += gkyl_time_diff_now_sec(wst);
}

void
vp_field_write_energy(gkyl_vlasov_app *app)
{
  struct timespec wst = gkyl_wall_clock();

  // Write out integrated field energy.
  const char *fmt = "%s-field-energy.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name);
  char fileNm[sz + 1]; // Ensures no buffer overflow.
  snprintf(fileNm, sizeof fileNm, fmt, app->name);

  int rank;
  gkyl_comm_get_rank(app->comm, &rank);

  if (rank == 0) {
    if (app->field->is_first_energy_write_call) {
      // Write to a new file (this ensure previous output is removed).
      gkyl_dynvec_write(app->field->integ_energy, fileNm);
      app->field->is_first_energy_write_call = false;
    } else {
      // Append to existing file.
      gkyl_dynvec_awrite(app->field->integ_energy, fileNm);
    }
  }
  gkyl_dynvec_clear(app->field->integ_energy);

  app->stat.n_field_diag_io += 1;
  app->stat.field_diag_io_tm += gkyl_time_diff_now_sec(wst);
}

void
vp_field_release(const gkyl_vlasov_app *app, struct vm_field *vpf)
{
  // Release resources for Vlasov-Poisson field.

  gkyl_dynvec_release(vpf->integ_energy);
  gkyl_array_integrate_release(vpf->calc_es_energy);
  gkyl_array_release(vpf->es_energy_fac);

  if (app->use_gpu) {
    gkyl_cu_free(vpf->es_energy_red);
    gkyl_cu_free(vpf->es_energy_red_global);
  } else {
    gkyl_free(vpf->es_energy_red);
    gkyl_free(vpf->es_energy_red_global);
  }

  gkyl_fem_poisson_release(vpf->fem_poisson);

  gkyl_array_release(vpf->epsilon);

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

  gkyl_array_release(vpf->phi_host);
  gkyl_array_release(vpf->phi);
  gkyl_array_release(vpf->phi_global);

  gkyl_array_release(vpf->rho_c_global);
  gkyl_array_release(vpf->rho_c);

  free(vpf);
}
