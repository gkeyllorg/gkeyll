#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_gk_field_priv.h>

#include <assert.h>
#include <math.h>

void
gk_field_flr_new(struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  assert(app->cdim > 1);
  f->invert_flr = gk_field_invert_flr;

  // Reference (squared) gyroradius in the operator A = 1 - rho^2*nabla_perp^2
  // used to retrieve phi from the modified potential Phi_0 (step 4 of the
  // algorithm in flr_effects.tex).
  double rhoSq_ref = 0.0;
  if (f->info.flr.type == GKYL_GK_FLR_PADE_CONST) {
    // Single reference species (e.g. the main ion).
    assert(f->info.flr.gyroradius > 0.0);
    rhoSq_ref = pow(f->info.flr.gyroradius, 2.0);
  }
  else {
    // polarization-weighted average of the species gyroradii,
    //   rho^2 = sum_s eps_s0*rho_s0^2 / sum_s eps_s0,  eps_s0 = n_s0*m_s/B^2,
    double polarization_bmag = f->info.polarization_bmag ? f->info.polarization_bmag : app->bmag_ref;
    double eps_sum = 0.0;
    for (int i = 0; i < app->num_species; ++i) {
      struct gk_species *s = &app->species[i];
      double gyroradius_bmag = s->info.flr.bmag ? s->info.flr.bmag : app->bmag_ref;
      double rhoSq_s = s->info.flr.Tperp * s->info.mass / pow(s->info.charge * gyroradius_bmag, 2.0);
      double eps_s0 = s->info.polarization_density * s->info.mass / pow(polarization_bmag, 2.0);
      rhoSq_ref += eps_s0 * rhoSq_s;
      eps_sum += eps_s0;
    }
    rhoSq_ref /= eps_sum;
  }
  // The apply passes boundary values through Dirichlet rows, so spatially
  // varying Dirichlet BCs are not supported with FLR effects.
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
  f->flr_kSq = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  gkyl_array_shiftc(f->flr_kSq, -pow(sqrt(2.0),app->cdim), 0); // Sets kSq=-1.

  f->flr_op = gkyl_fem_poisson_perp_new(&app->local, &app->grid, app->basis, &f->poisson_bcs, f->info.bias_line_list, f->flr_rhoSq, f->flr_kSq, app->use_gpu);

}

void
gk_field_invert_flr(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *phi)
{
  // Retrieve phi from the potential obtained from the FLR charge.
  // phi = (1 - rho^2*nabla_perp^2) Phi_0
  gkyl_fem_poisson_perp_lhs_apply(field->flr_op, phi, phi);
}

void
gk_field_invert_flr_none(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *phi)
{
}

void
gk_field_flr_release(const struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  gkyl_array_release(f->flr_rhoSq);
  gkyl_array_release(f->flr_kSq);
  gkyl_fem_poisson_perp_release(f->flr_op);
}
