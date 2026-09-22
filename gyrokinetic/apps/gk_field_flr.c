#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_gk_field_priv.h>

#include <assert.h>
#include <math.h>

void
gk_field_flr_new(struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  assert(app->cdim > 1);
  f->invert_flr = f->info.flr.use_fem_operator? gk_field_invert_flr_op : gk_field_invert_flr_const;

  // Reference (squared) gyroradius in the operator A = 1 - rho^2*nabla_perp^2
  // used to retrieve phi from the modified potential Phi_0 (step 4 of the
  // algorithm in DR #797).
  double polarization_bmag = f->info.polarization_bmag ? f->info.polarization_bmag : app->bmag_ref;
  double rhoSq_ref = 0.0;
  if (f->info.flr.avg_gyroradius) {
    // Polarization-weighted average of the species gyroradii,
    //   rho^2 = sum_s eps_s0*rho_s0^2 / sum_s eps_s0,  eps_s0 = n_s0*m_s/B^2,
    double eps_sum = 0.0;
    for (int i = 0; i < app->num_species; ++i) {
      struct gk_species *s = &app->species[i];
      double eps_s0 = s->info.polarization_density * s->info.mass / pow(polarization_bmag, 2.0);
      rhoSq_ref += eps_s0 * s->flr_rhoSq_ref;
      eps_sum += eps_s0;
    }
    assert(eps_sum > 0.0);
    rhoSq_ref /= eps_sum;
  }
  else {
    // Gyroradius of the gyroaveraged species with the largest polarization
    // weight (the main ion), so steps 2-5 of DR #797 use the same rho.
    double eps_max = -1.0;
    for (int i = 0; i < app->num_species; ++i) {
      struct gk_species *s = &app->species[i];
      if (!s->use_flr) continue;
      double eps_s0 = s->info.polarization_density * s->info.mass / pow(polarization_bmag, 2.0);
      if (eps_s0 > eps_max) {
        rhoSq_ref = s->flr_rhoSq_ref;
        eps_max = eps_s0;
      }
    }
  }
  // At least one species must provide a reference gyroradius.
  assert(rhoSq_ref > 0.0);

  // Modified potential Phi_0 and a buffer, used in the local term and the
  // field energy diagnostic.
  f->flr_phi0 = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  f->flr_buff = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  f->flr_energy_red = app->use_gpu? gkyl_cu_malloc(sizeof(double[1])) : gkyl_malloc(sizeof(double[1]));

  if (!f->info.flr.use_fem_operator) {
    // Use the simplification phi = Phi_0 + (rho_i^2/eps_pol)*rho_c/J
    double polarization_weight = 0.0;
    for (int i=0; i<app->num_species; ++i) {
      struct gk_species *s = &app->species[i];
      polarization_weight += s->info.polarization_density*s->info.mass/pow(polarization_bmag,2);
    }
    f->flr_local_fac = rhoSq_ref/polarization_weight;
    return;
  }

  // Spatially varying Dirichlet BCs are not supported with the FEM operator.
  for (int d = 0; d < app->cdim - 1; d++) {
    assert(f->poisson_bcs.lo_type[d] != GKYL_POISSON_DIRICHLET_VARYING);
    assert(f->poisson_bcs.up_type[d] != GKYL_POISSON_DIRICHLET_VARYING);
  }

  // Weight in the perpendicular Laplacian of A = 1 - rho^2*nabla_perp^2.
  f->flr_rhoSq = mkarr(app->use_gpu, (2*(app->cdim/3)+1)*app->basis.num_basis, app->local_ext.volume);
  struct gkyl_array *Jgij[3] = {app->gk_geom->geo_int.gxxj, app->gk_geom->geo_int.gxyj, app->gk_geom->geo_int.gyyj};
  for (int i=0; i<app->cdim-2/app->cdim; i++) {
    gkyl_array_set_offset(f->flr_rhoSq, rhoSq_ref, Jgij[i], i*app->basis.num_basis);
  }
  // The Laplacian weight carries the Jacobian (J*g^ij), so the identity term
  // must too: kSq=-J so the operator is J*(1 - rho^2*nabla_perp^2).
  f->flr_kSq = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  gkyl_array_set(f->flr_kSq, -1.0, app->gk_geom->geo_int.jacobgeo);

  f->flr_op = gkyl_fem_poisson_perp_new(&app->local, &app->grid, app->basis, &f->poisson_bcs, f->info.bias_line_list, f->flr_rhoSq, f->flr_kSq, app->use_gpu);
}

void
gk_field_invert_flr_const(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *phi)
{
  // phi = Phi_0 + (rho_i^2/eps_pol)*rho_c/J, rho_c the (J-weighted) charge
  // density used in the Poisson solve.
  gkyl_array_copy(field->flr_phi0, phi);
  gkyl_dg_mul_op_range(&app->basis, 0, field->flr_buff, 0, app->gk_geom->geo_int.jacobgeo_inv,
    0, field->rho_c, &app->local);
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
}

void
gk_field_flr_release(const struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  if (f->info.flr.use_fem_operator) {
    gkyl_array_release(f->flr_rhoSq);
    gkyl_array_release(f->flr_kSq);
    gkyl_fem_poisson_perp_release(f->flr_op);
  }
  gkyl_array_release(f->flr_phi0);
  gkyl_array_release(f->flr_buff);
  if (app->use_gpu)
    gkyl_cu_free(f->flr_energy_red);
  else
    gkyl_free(f->flr_energy_red);
}
