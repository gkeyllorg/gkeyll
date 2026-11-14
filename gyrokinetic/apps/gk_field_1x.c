#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_position_map.h>
#include <gkyl_util.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_gk_field_priv.h>
#include <gkyl_array_rio_priv.h>
#include <gkyl_comm_io.h>

#include <assert.h>
#include <float.h>
#include <time.h>

static void
gk_field_rhs_phi_1x(struct gkyl_gyrokinetic_app *app, struct gk_field *field)
{
  // Solve the Poisson equation in 1x with the parallel FEM projection.
  gk_field_fem_projection_par(app, field, field->rho_c, field->phi_smooth);
}

static void
gk_field_fem_release_1x(const gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  gkyl_array_release(gkf->rho_c);
  gkyl_array_release(gkf->rho_c_global_dg);
  gkyl_array_release(gkf->rho_c_global_smooth);
  gkyl_array_release(gkf->phi_smooth);
  gkyl_array_release(gkf->phi_fem);

  if (gkf->gkfield_id == GKYL_GK_FIELD_EM) {
    gkyl_array_release(gkf->apar_fem);
    gkyl_array_release(gkf->apardot_fem);
  }

  if (app->use_gpu) {
    gkyl_array_release(gkf->phi_host);
  }

  gkyl_array_release(gkf->epsilon);

  gkyl_fem_parproj_release(gkf->fem_parproj);

  gkyl_array_integrate_release(gkf->calc_em_energy);
}

void
gk_field_fem_new_1x(struct gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  // Create global subrange we'll copy the field solver solution from (into local).
  gkyl_sub_range_intersect(&gkf->global_sub_range, &app->global, &app->local);

  // Allocate arrays for charge density.
  gkf->rho_c = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  gkf->rho_c_global_dg = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);
  gkf->rho_c_global_smooth = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);

  // Allocate arrays for electrostatic potential.
  gkf->phi_fem = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);
  gkf->phi_smooth = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

  // Allocate electromagnetic arrays if needed.
  if (gkf->gkfield_id == GKYL_GK_FIELD_EM) {
    gkf->apar_fem = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    gkf->apardot_fem = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  }

  // Allocate phi_host for I/O.
  gkf->phi_host = gkf->phi_smooth;
  if (app->use_gpu) {
    gkf->phi_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
  }

  gkf->rhs_phi_func = gk_field_rhs_phi_1x;

  // Allocate array for the polarization weight times geometric coefficients.
  gkf->epsilon =
    mkarr(app->use_gpu, (2 * (app->cdim / 3) + 1) * app->basis.num_basis, app->local_ext.volume);

  double polarization_weight = 0.0;
  double polarization_bmag = gkf->info.polarization_bmag ? gkf->info.polarization_bmag :
                                                           app->bmag_ref;
  // Linearized polarization density
  for (int i = 0; i < app->num_species; ++i) {
    struct gk_species *gks = &app->species[i];
    polarization_weight +=
      gks->info.polarization_density * gks->info.mass / pow(polarization_bmag, 2);
  }
  // Need to set weight to kperpsq*polarizationWeight for use in potential smoothing.
  assert(
    gkf->info.kperpSq > 0.0 || gkf->has_adiabatic_species ||
    (!gkf->calc_init_field && !gkf->update_field)
  );
  gkyl_array_copy(gkf->epsilon, app->gk_geom->geo_int.jacobgeo);
  gkyl_array_scale(gkf->epsilon, polarization_weight);
  gkyl_array_scale(gkf->epsilon, gkf->info.kperpSq);

  gkf->accumulate_rhoc_func = gk_field_accumulate_rho_c_poisson;
  gkf->es_energy_fac_1d = 0.5 * polarization_weight * gkf->info.kperpSq;
  // The parallel projection solves the field equation with weight epsilon. An adiabatic species
  // adds K to epsilon, so the projection is then created once K is known.
  if (!gkf->has_adiabatic_species) {
    gk_field_fem_parproj_new(app, gkf, gkf->epsilon);
  }

  gkf->calc_em_energy =
    gkyl_array_integrate_new(&app->grid, &app->basis, 1, GKYL_ARRAY_INTEGRATE_OP_SQ, app->use_gpu);

  // Create operator needed for FLR effects.
  gkf->use_flr = false;
  gkf->invert_flr = gk_field_invert_flr_none;
  for (int i = 0; i < app->num_species; ++i) {
    struct gk_species *gks = &app->species[i];
    if (gks->info.flr.type) {
      gkf->use_flr = gkf->use_flr || gks->info.flr.type;
    }
  }

  gkf->release_func = gk_field_fem_release_1x;
}
