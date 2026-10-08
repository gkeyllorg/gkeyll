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

// Integral of a conf-space array over the global domain.
static double
gk_field_1x_integral(struct gkyl_gyrokinetic_app *app, const struct gkyl_array *arr)
{
  double *integ = app->use_gpu ? gkyl_cu_malloc(sizeof(double)) : gkyl_malloc(sizeof(double));
  double *integ_red = app->use_gpu ? gkyl_cu_malloc(sizeof(double)) : gkyl_malloc(sizeof(double));
  struct gkyl_array_integrate *int_op = gkyl_array_integrate_new(
    &app->grid, &app->basis, 1, GKYL_ARRAY_INTEGRATE_OP_NONE, app->use_gpu
  );
  gkyl_array_integrate_advance(int_op, arr, 1.0, NULL, &app->local, NULL, integ);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, integ, integ_red);

  double integ_ho[1];
  if (app->use_gpu) {
    gkyl_cu_memcpy(integ_ho, integ_red, sizeof(double), GKYL_CU_MEMCPY_D2H);
    gkyl_cu_free(integ);
    gkyl_cu_free(integ_red);
  } else {
    integ_ho[0] = integ_red[0];
    gkyl_free(integ);
    gkyl_free(integ_red);
  }
  gkyl_array_integrate_release(int_op);
  return integ_ho[0];
}

// Projection onto the continuous parallel basis, weighted by epsilon (solves the 1x field equation).
static void
gk_field_1x_parproj_new(struct gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  // Gather epsilon for (global) smoothing in z.
  struct gkyl_array *epsilon_global =
    mkarr(app->use_gpu, gkf->epsilon->ncomp, app->global_ext.volume);
  gkyl_comm_array_allgather(app->comm, &app->local, &app->global, gkf->epsilon, epsilon_global);

  enum gkyl_fem_parproj_bc_type fem_parproj_bc = GKYL_FEM_PARPROJ_NONE;
  for (int d = 0; d < app->num_periodic_dir; ++d) {
    if (app->periodic_dirs[d] == app->cdim - 1) {
      fem_parproj_bc = GKYL_FEM_PARPROJ_PERIODIC;
    }
  }

  gkf->fem_parproj = gkyl_fem_parproj_new(
    &app->global, &app->grid, &app->basis, fem_parproj_bc, 0, epsilon_global, 0, app->use_gpu
  );
  gkyl_array_release(epsilon_global);
}

// No flux-surface average in 1x: phi - <phi> = phi.
static void
gk_field_adiabatic_dphi_1x(
  gkyl_gyrokinetic_app *app, const struct gk_field *gkf, const struct gkyl_array *phi,
  struct gkyl_array *out
)
{
  gkyl_array_copy(out, phi);
}

// Add K = (q_s^2 n_s0/T_s)*J of the adiabatic species, from its t=0 moments, to the weight.
static void
gk_field_adiabatic_init_1x(gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  gk_field_adiabatic_profiles_calc(app, gkf);

  gkyl_array_accumulate(gkf->epsilon, 1.0, gkf->adiab.kJ);
  gk_field_1x_parproj_new(app, gkf);

  // The 1x energy uses a scalar factor, the J-weighted average of q_s^2 n_s0/T_s
  // (exact for uniform profiles).
  gkf->es_energy_fac_1d += 0.5 * gk_field_1x_integral(app, gkf->adiab.kJ) /
                           gk_field_1x_integral(app, app->gk_geom->geo_int.jacobgeo);
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
    if (i == gkf->adiab.species_idx) {
      continue; // An adiabatic species has no polarization.
    }
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
  if (gkf->has_adiabatic_species) {
    if (app->species[gkf->adiab.species_idx].info.charge < 0.0) {
      gkf->init_adiab_func = gk_field_adiabatic_init_1x;
      gkf->adiab.dphi_func = gk_field_adiabatic_dphi_1x;
    } else {
      assert(false); // Not implemented for ions.
    }
  } else {
    gk_field_1x_parproj_new(app, gkf);
  }

  gkf->calc_em_energy =
    gkyl_array_integrate_new(&app->grid, &app->basis, 1, GKYL_ARRAY_INTEGRATE_OP_SQ, app->use_gpu);

  // FLR effects are not implemented for this field type.
  gkf->use_flr = gkf->info.flr.type != GKYL_GK_FLR_NONE;
  gkf->invert_flr = gk_field_invert_flr_none;
  assert(!gkf->use_flr);

  gkf->release_func = gk_field_fem_release_1x;
}
