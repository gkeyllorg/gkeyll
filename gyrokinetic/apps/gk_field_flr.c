#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_gk_field_priv.h>
#include <gkyl_gk_species_priv.h>

#include <assert.h>
#include <math.h>

// Finite Larmor radius (FLR) model of the field equation, chosen in the field input:
// - Pade: the polarization is -div(eps grad_perp) A^{-1} with A = 1 - rho^2*nabla_perp^2, so
//   the field equation is solved for the modified potential Phi_0 without FLR effects and
//   phi = A Phi_0 (invert_flr).
// - gyromatrix: the polarization is that of the gyroaverages of the species (see
//   gk_species_flr_gyromatrix_accumulate_polarization). The field equation is solved
//   iteratively, with the Pade solve as preconditioner.

/** Pade approximant of the polarization. **/

static void
gk_field_flr_pade_es_energy(
  gkyl_gyrokinetic_app *app, const struct gk_field *field, double factor, double *out
)
{
  // The ES energy is (1/2) int eps grad_perp(phi).grad_perp(Phi_0),
  // computed with the polarization identity <a,b> = (<a+b,a+b> - <a-b,a-b>)/4.
  gkyl_array_set(field->flr_buff, 0.5, field->phi_smooth);
  gkyl_array_accumulate(field->flr_buff, 0.5, field->flr_phi0);
  gkyl_array_integrate_advance(
    field->calc_em_energy, field->flr_buff, factor, field->es_energy_fac, &app->local, &app->local,
    out
  );

  gkyl_array_set(field->flr_buff, 0.5, field->phi_smooth);
  gkyl_array_accumulate(field->flr_buff, -0.5, field->flr_phi0);
  gkyl_array_integrate_advance(
    field->calc_em_energy, field->flr_buff, factor, field->es_energy_fac, &app->local, &app->local,
    field->flr_energy_red
  );
  if (app->use_gpu) {
    double e[2];
    gkyl_cu_memcpy(&e[0], out, sizeof(double), GKYL_CU_MEMCPY_D2H);
    gkyl_cu_memcpy(&e[1], field->flr_energy_red, sizeof(double), GKYL_CU_MEMCPY_D2H);
    e[0] -= e[1];
    gkyl_cu_memcpy(out, &e[0], sizeof(double), GKYL_CU_MEMCPY_H2D);
  } else {
    out[0] -= field->flr_energy_red[0];
  }
}

static void
gk_field_flr_pade_new(struct gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  assert(app->cdim > 1);
  gkf->es_energy_func = gk_field_flr_pade_es_energy;
  gkf->invert_flr = gkf->info.flr.use_fem_operator ? gk_field_invert_flr_op :
                                                     gk_field_invert_flr_const;

  // Reference (squared) gyroradius in the operator A = 1 - rho^2*nabla_perp^2
  // used to retrieve phi from the modified potential Phi_0 (step 4 of the
  // algorithm in DR #797). Adiabatic species have no polarization, so the
  // sums run over the species in rho_c.
  double polarization_bmag = gkf->info.polarization_bmag ? gkf->info.polarization_bmag :
                                                           app->bmag_ref;
  double rhoSq_ref = 0.0, polarization_weight = 0.0;
  if (gkf->info.flr.avg_gyroradius) {
    // Polarization-weighted average of the species gyroradii,
    //   rho^2 = sum_s eps_s0*rho_s0^2 / sum_s eps_s0,  eps_s0 = n_s0*m_s/B^2,
    for (int i = 0; i < gkf->num_rhoc_species; ++i) {
      struct gk_species *gks = &app->species[gkf->rhoc_species_idx[i]];
      double eps_s0 = gks->info.polarization_density * gks->info.mass / pow(polarization_bmag, 2.0);
      rhoSq_ref += eps_s0 * gks->flr_rhoSq_ref;
      polarization_weight += eps_s0;
    }
    assert(polarization_weight > 0.0);
    rhoSq_ref /= polarization_weight;
  } else {
    // Gyroradius of the gyroaveraged species with the largest polarization
    // weight (the main ion), so steps 2-5 of DR #797 use the same rho.
    double eps_max = -1.0;
    for (int i = 0; i < gkf->num_rhoc_species; ++i) {
      struct gk_species *gks = &app->species[gkf->rhoc_species_idx[i]];
      double eps_s0 = gks->info.polarization_density * gks->info.mass / pow(polarization_bmag, 2.0);
      polarization_weight += eps_s0;
      if (gks->use_flr && eps_s0 > eps_max) {
        rhoSq_ref = gks->flr_rhoSq_ref;
        eps_max = eps_s0;
      }
    }
  }
  // At least one species must provide a reference gyroradius.
  assert(rhoSq_ref > 0.0);
  // rho_i^2/eps_pol, used in the local term and in the response of an adiabatic species.
  gkf->flr_local_fac = rhoSq_ref / polarization_weight;

  // Modified potential Phi_0 and a buffer, used in the local term and the
  // field energy diagnostic.
  gkf->flr_phi0 = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  gkf->flr_buff = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  gkf->flr_energy_red = app->use_gpu ? gkyl_cu_malloc(sizeof(double[1])) :
                                       gkyl_malloc(sizeof(double[1]));

  if (!gkf->info.flr.use_fem_operator) {
    return; // The local term phi = Phi_0 + (rho_i^2/eps_pol)*rho_c/J needs no operator.
  }

  // Spatially varying Dirichlet BCs are not supported with the FEM operator.
  for (int d = 0; d < app->cdim - 1; d++) {
    assert(gkf->poisson_bcs.lo_type[d] != GKYL_POISSON_DIRICHLET_VARYING);
    assert(gkf->poisson_bcs.up_type[d] != GKYL_POISSON_DIRICHLET_VARYING);
  }

  // Weight in the perpendicular Laplacian of A = 1 - rho^2*nabla_perp^2.
  gkf->flr_rhoSq =
    mkarr(app->use_gpu, (2 * (app->cdim / 3) + 1) * app->basis.num_basis, app->local_ext.volume);
  struct gkyl_array *Jgij[3] = {
    app->gk_geom->geo_int.gxxj, app->gk_geom->geo_int.gxyj, app->gk_geom->geo_int.gyyj
  };
  for (int i = 0; i < app->cdim - 2 / app->cdim; i++) {
    gkyl_array_set_offset(gkf->flr_rhoSq, rhoSq_ref, Jgij[i], i * app->basis.num_basis);
  }
  // The Laplacian weight carries the Jacobian (J*g^ij), so the identity term
  // must too: kSq=-J so the operator is J*(1 - rho^2*nabla_perp^2).
  gkf->flr_kSq = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  gkyl_array_set(gkf->flr_kSq, -1.0, app->gk_geom->geo_int.jacobgeo);

  gkf->flr_op = gkyl_fem_poisson_perp_new(
    &app->local, &app->grid, app->basis, &gkf->poisson_bcs, gkf->info.bias_line_list,
    gkf->flr_rhoSq, gkf->flr_kSq, app->use_gpu
  );
}

void
gk_field_invert_flr_const(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *phi)
{
  // phi = Phi_0 + (rho_i^2/eps_pol)*rho_c/J, rho_c the (J-weighted) charge
  // density used in the Poisson solve.
  gkyl_array_copy(field->flr_phi0, phi);
  gkyl_dg_mul_op_range(
    &app->basis, 0, field->flr_buff, 0, app->gk_geom->geo_int.jacobgeo_inv, 0, field->rho_c,
    &app->local
  );
  gkyl_array_accumulate_range(phi, field->flr_local_fac, field->flr_buff, &app->local);
}

void
gk_field_invert_flr_op(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *phi)
{
  // phi = (1 - rho^2*nabla_perp^2) Phi_0 with the FEM perpendicular operator.
  gkyl_array_copy(field->flr_phi0, phi);
  gkyl_fem_poisson_perp_lhs_apply(field->flr_op, phi, phi);
}

void
gk_field_invert_flr_none(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *phi)
{
  // Do nothing.
}

static void
gk_field_flr_pade_release(const struct gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  if (gkf->info.flr.use_fem_operator) {
    gkyl_array_release(gkf->flr_rhoSq);
    gkyl_array_release(gkf->flr_kSq);
    gkyl_fem_poisson_perp_release(gkf->flr_op);
  }
  gkyl_array_release(gkf->flr_phi0);
  gkyl_array_release(gkf->flr_buff);
  if (app->use_gpu) {
    gkyl_cu_free(gkf->flr_energy_red);
  } else {
    gkyl_free(gkf->flr_energy_red);
  }
}

/** Gyromatrix: polarization from the gyroaverages of the species. **/

// Dimension of the Krylov space, maximum number of restarts and tolerance (on the
// preconditioned residual relative to the solution) of the iterative solve of the field equation.
#define GK_FIELD_FLR_KRYLOV_DIM 10
#define GK_FIELD_FLR_MAX_RESTARTS 5
#define GK_FIELD_FLR_TOL 1.0e-6

static double
gk_field_flr_gyromatrix_inner(
  gkyl_gyrokinetic_app *app, const struct gk_field *field, const struct gkyl_array *lhs,
  const struct gkyl_array *rhs
)
{
  // Inner product int lhs*rhs dx over the global domain.
  const struct gk_field_flr_gyromatrix *gm = &field->flr_gyromatrix;
  gkyl_dg_mul_op_range(&app->basis, 0, gm->inner_buff, 0, lhs, 0, rhs, &app->local);
  gkyl_array_integrate_advance(
    gm->integ_op, gm->inner_buff, 1.0, NULL, &app->local, NULL, gm->integ_red
  );
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, gm->integ_red, gm->integ_red_global);

  double inner[1];
  if (app->use_gpu) {
    gkyl_cu_memcpy(inner, gm->integ_red_global, sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    inner[0] = gm->integ_red_global[0];
  }
  return inner[0];
}

static void
gk_field_flr_gyromatrix_adiabatic_none(
  gkyl_gyrokinetic_app *app, const struct gk_field *field, const struct gkyl_array *phi,
  struct gkyl_array *out
)
{
}

static void
gk_field_flr_gyromatrix_adiabatic(
  gkyl_gyrokinetic_app *app, const struct gk_field *field, const struct gkyl_array *phi,
  struct gkyl_array *out
)
{
  // out += K (phi - <phi>), the response of the adiabatic species.
  const struct gk_field_flr_gyromatrix *gm = &field->flr_gyromatrix;
  field->gk_adiab.dphi_func(app, field, phi, gm->inner_buff);
  gkyl_dg_mul_op_range(
    &app->basis, 0, gm->buff, 0, field->gk_adiab.kJ, 0, gm->inner_buff, &app->local
  );
  gkyl_array_accumulate_range(out, 1.0, gm->buff, &app->local);
}

static void
gk_field_flr_gyromatrix_polarization(
  gkyl_gyrokinetic_app *app, const struct gk_field *field, const struct gkyl_array *phi,
  struct gkyl_array *out
)
{
  // out = polarization charge density (times J) of the species.
  gkyl_array_clear(out, 0.0);
  for (int k = 0; k < field->num_rhoc_species; k++) {
    struct gk_species *gks = &app->species[field->rhoc_species_idx[k]];
    gk_species_flr_gyromatrix_accumulate_polarization(app, gks, phi, out);
  }
}

static void
gk_field_flr_gyromatrix_lhs(
  gkyl_gyrokinetic_app *app, const struct gk_field *field, const struct gkyl_array *phi,
  struct gkyl_array *out
)
{
  // out = A phi, with A the polarization of the species plus the response of the adiabatic
  // species.
  const struct gk_field_flr_gyromatrix *gm = &field->flr_gyromatrix;
  gk_field_flr_gyromatrix_polarization(app, field, phi, out);
  gm->accumulate_adiabatic_func(app, field, phi, out);
}

static void
gk_field_flr_gyromatrix_precond(
  gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *out
)
{
  // out = M^{-1} rho_c, with M^{-1} the Pade solve of the field equation, which takes the
  // charge density in rho_c and returns the potential in phi_smooth.
  const struct gk_field_flr_gyromatrix *gm = &field->flr_gyromatrix;
  gm->precond_func(app, field);
  gkyl_array_copy(out, field->phi_smooth);
}

static void
gk_field_flr_gyromatrix_rhs_phi(struct gkyl_gyrokinetic_app *app, struct gk_field *field)
{
  // Solve the field equation A phi = rho with GMRES, preconditioned from the left by the Pade
  // solve M^{-1} (M^{-1} A is close to the identity), starting from the previous solution.
  struct gk_field_flr_gyromatrix *gm = &field->flr_gyromatrix;
  int kdim = GK_FIELD_FLR_KRYLOV_DIM;
  // Hessenberg matrix (column-major with kdim+1 rows), Givens rotations that make it
  // triangular, right-hand side of the least squares problem and its solution.
  double hess[(GK_FIELD_FLR_KRYLOV_DIM + 1) * GK_FIELD_FLR_KRYLOV_DIM];
  double rot_cos[GK_FIELD_FLR_KRYLOV_DIM], rot_sin[GK_FIELD_FLR_KRYLOV_DIM];
  double lsq_rhs[GK_FIELD_FLR_KRYLOV_DIM + 1], coeff[GK_FIELD_FLR_KRYLOV_DIM];

  // Keep the charge density (rho_c is the input of the preconditioner).
  gkyl_array_copy(gm->charge, field->rho_c);

  // The residual is compared to the norm of the solution, known from the previous solve (so
  // that no additional application of the preconditioner is needed).
  double sol_norm = sqrt(gk_field_flr_gyromatrix_inner(app, field, gm->sol, gm->sol));
  double resid_norm = 0.0;
  for (int cycle = 0; cycle < GK_FIELD_FLR_MAX_RESTARTS; cycle++) {
    // Residual of the preconditioned system, M^{-1} (rho - A sol).
    gk_field_flr_gyromatrix_lhs(app, field, gm->sol, field->rho_c);
    gkyl_array_scale(field->rho_c, -1.0);
    gkyl_array_accumulate(field->rho_c, 1.0, gm->charge);
    gk_field_flr_gyromatrix_precond(app, field, gm->krylov[0]);
    resid_norm = sqrt(gk_field_flr_gyromatrix_inner(app, field, gm->krylov[0], gm->krylov[0]));
    if (sol_norm == 0.0) {
      // No previous solution: the residual is M^{-1} rho, close to the solution.
      sol_norm = resid_norm;
    }
    double tol = GK_FIELD_FLR_TOL * sol_norm;
    if (resid_norm <= tol) {
      break;
    }
    gkyl_array_scale(gm->krylov[0], 1.0 / resid_norm);
    lsq_rhs[0] = resid_norm;

    int num_vec = 0;
    for (int k = 0; k < kdim && resid_norm > tol; k++) {
      // Arnoldi step with modified Gram-Schmidt orthogonalization.
      double *hess_k = &hess[k * (kdim + 1)];
      gk_field_flr_gyromatrix_lhs(app, field, gm->krylov[k], field->rho_c);
      gk_field_flr_gyromatrix_precond(app, field, gm->krylov[k + 1]);
      for (int i = 0; i < k + 1; i++) {
        hess_k[i] = gk_field_flr_gyromatrix_inner(app, field, gm->krylov[k + 1], gm->krylov[i]);
        gkyl_array_accumulate(gm->krylov[k + 1], -hess_k[i], gm->krylov[i]);
      }
      hess_k[k + 1] =
        sqrt(gk_field_flr_gyromatrix_inner(app, field, gm->krylov[k + 1], gm->krylov[k + 1]));
      if (hess_k[k + 1] > 0.0) {
        gkyl_array_scale(gm->krylov[k + 1], 1.0 / hess_k[k + 1]);
      }

      // Apply the previous rotations to the new column, and a new rotation that eliminates
      // its last entry.
      for (int i = 0; i < k; i++) {
        double rotated = rot_cos[i] * hess_k[i] + rot_sin[i] * hess_k[i + 1];
        hess_k[i + 1] = -rot_sin[i] * hess_k[i] + rot_cos[i] * hess_k[i + 1];
        hess_k[i] = rotated;
      }
      double diag = hypot(hess_k[k], hess_k[k + 1]);
      rot_cos[k] = hess_k[k] / diag;
      rot_sin[k] = hess_k[k + 1] / diag;
      hess_k[k] = diag;
      lsq_rhs[k + 1] = -rot_sin[k] * lsq_rhs[k];
      lsq_rhs[k] = rot_cos[k] * lsq_rhs[k];

      resid_norm = fabs(lsq_rhs[k + 1]);
      num_vec = k + 1;
      gm->num_iter += 1;
    }

    // Solve the triangular least squares problem and update the solution.
    for (int i = num_vec - 1; i >= 0; i--) {
      coeff[i] = lsq_rhs[i];
      for (int j = i + 1; j < num_vec; j++) {
        coeff[i] -= hess[j * (kdim + 1) + i] * coeff[j];
      }
      coeff[i] /= hess[i * (kdim + 1) + i];
    }
    for (int i = 0; i < num_vec; i++) {
      gkyl_array_accumulate(gm->sol, coeff[i], gm->krylov[i]);
    }
  }

  gkyl_array_copy(field->phi_smooth, gm->sol);
  gkyl_array_copy(field->rho_c, gm->charge);

  if (gm->num_solve == 0) {
    gkyl_gyrokinetic_app_cout(
      app, stdout, "Gyromatrix field solve: %ld iterations, relative residual %.3e.\n",
      gm->num_iter, sol_norm > 0.0 ? resid_norm / sol_norm : 0.0
    );
  }
  gm->num_solve += 1;
}

static void
gk_field_flr_gyromatrix_es_energy(
  gkyl_gyrokinetic_app *app, const struct gk_field *field, double factor, double *out
)
{
  // The ES energy is (1/2) int phi*(P phi) dx, with P the polarization of the species.
  const struct gk_field_flr_gyromatrix *gm = &field->flr_gyromatrix;
  gk_field_flr_gyromatrix_polarization(app, field, field->phi_smooth, gm->buff);
  gkyl_dg_mul_op_range(
    &app->basis, 0, gm->inner_buff, 0, field->phi_smooth, 0, gm->buff, &app->local
  );
  gkyl_array_integrate_advance(
    gm->integ_op, gm->inner_buff, 0.5 * factor, NULL, &app->local, NULL, out
  );
}

static void
gk_field_flr_gyromatrix_new(struct gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  struct gk_field_flr_gyromatrix *gm = &gkf->flr_gyromatrix;

  // The polarization of every species must be that of its gyroaverages.
  for (int k = 0; k < gkf->num_rhoc_species; k++) {
    const struct gk_species *gks = &app->species[gkf->rhoc_species_idx[k]];
    if (!gks->use_flr) {
      gkyl_exit("Gyromatrix FLR model: every species in the field equation needs FLR parameters.");
    }
  }
  // The potential must be periodic or vanish at the perpendicular boundaries (the gyroaverages
  // extend it beyond them as an odd function), and have no biases (the preconditioner must be
  // linear).
  assert(!gkf->is_dirichletvar);
  assert(!gkf->info.bias_line_list || gkf->info.bias_line_list->num_bias_line == 0);
  for (int d = 0; d < app->cdim - 1; d++) {
    assert(
      gkf->poisson_bcs.lo_type[d] == GKYL_POISSON_PERIODIC ||
      (gkf->poisson_bcs.lo_type[d] == GKYL_POISSON_DIRICHLET &&
       gkf->poisson_bcs.lo_value[d].v[0] == 0.0)
    );
    assert(
      gkf->poisson_bcs.up_type[d] == GKYL_POISSON_PERIODIC ||
      (gkf->poisson_bcs.up_type[d] == GKYL_POISSON_DIRICHLET &&
       gkf->poisson_bcs.up_value[d].v[0] == 0.0)
    );
  }

  // The field equation is solved with the Pade solve as preconditioner.
  gm->precond_func = gkf->rhs_phi_func;
  gkf->rhs_phi_func = gk_field_flr_gyromatrix_rhs_phi;
  gkf->es_energy_func = gk_field_flr_gyromatrix_es_energy;
  gm->accumulate_adiabatic_func = gkf->has_adiabatic_species ?
                                    gk_field_flr_gyromatrix_adiabatic :
                                    gk_field_flr_gyromatrix_adiabatic_none;

  int num_basis = app->basis.num_basis;
  long size = app->local_ext.volume;
  gm->sol = mkarr(app->use_gpu, num_basis, size);
  gm->charge = mkarr(app->use_gpu, num_basis, size);
  gm->buff = mkarr(app->use_gpu, num_basis, size);
  gm->inner_buff = mkarr(app->use_gpu, num_basis, size);
  gm->krylov = gkyl_malloc(sizeof(struct gkyl_array * [GK_FIELD_FLR_KRYLOV_DIM + 1]));
  for (int i = 0; i < GK_FIELD_FLR_KRYLOV_DIM + 1; i++) {
    gm->krylov[i] = mkarr(app->use_gpu, num_basis, size);
  }
  gkyl_array_clear(gm->sol, 0.0);

  gm->integ_op = gkyl_array_integrate_new(
    &app->grid, &app->basis, 1, GKYL_ARRAY_INTEGRATE_OP_NONE, app->use_gpu
  );
  gm->integ_red = app->use_gpu ? gkyl_cu_malloc(sizeof(double)) : gkyl_malloc(sizeof(double));
  gm->integ_red_global = app->use_gpu ? gkyl_cu_malloc(sizeof(double)) :
                                        gkyl_malloc(sizeof(double));
  gm->num_solve = 0;
  gm->num_iter = 0;
}

static void
gk_field_flr_gyromatrix_release(const struct gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  struct gk_field_flr_gyromatrix *gm = &gkf->flr_gyromatrix;
  gkyl_array_release(gm->sol);
  gkyl_array_release(gm->charge);
  gkyl_array_release(gm->buff);
  gkyl_array_release(gm->inner_buff);
  for (int i = 0; i < GK_FIELD_FLR_KRYLOV_DIM + 1; i++) {
    gkyl_array_release(gm->krylov[i]);
  }
  gkyl_free(gm->krylov);
  gkyl_array_integrate_release(gm->integ_op);
  if (app->use_gpu) {
    gkyl_cu_free(gm->integ_red);
    gkyl_cu_free(gm->integ_red_global);
  } else {
    gkyl_free(gm->integ_red);
    gkyl_free(gm->integ_red_global);
  }
}

/** Selection of the model. **/

void
gk_field_flr_new(struct gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  // The Pade model is also the preconditioner of the gyromatrix model.
  gk_field_flr_pade_new(app, gkf);
  // The Pade model applies its operator after the Ampere and Ohm's law solves too. With the
  // gyromatrix model Ampere's law has no FLR operator: the current is gyroaveraged instead.
  gkf->invert_flr_apar = gkf->invert_flr;
  if (gkf->info.flr.type == GKYL_GK_FLR_GYROMATRIX) {
    gk_field_flr_gyromatrix_new(app, gkf);
    gkf->invert_flr_apar = gk_field_invert_flr_none;
  }
}

void
gk_field_flr_release(const struct gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  gk_field_flr_pade_release(app, gkf);
  if (gkf->info.flr.type == GKYL_GK_FLR_GYROMATRIX) {
    gk_field_flr_gyromatrix_release(app, gkf);
  }
}
