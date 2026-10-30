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
gk_field_ohm_solve_1x(
  struct gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *out
)
{
  struct timespec wst = gkyl_wall_clock();

  // Compute the LHS factor (= k_perp^2/mu_0 + sum_s q_s^2 n_s/m_s)
  // gkyl_array_set(field->dApartdtSlvr_lhs_factor, 1.0, field->dApartdtSlvr_kSq);
  gkyl_array_clear(field->dApartdtSlvr_lhs_factor, 0.0);
  gkyl_array_accumulate_range(
    field->dApartdtSlvr_lhs_factor, -1.0, field->dApartdtSlvr_kSq, &app->local
  );
  gkyl_array_accumulate_range(
    field->dApartdtSlvr_lhs_factor, 1.0, field->lapWeightAmpere, &app->local
  );

  // Weak division method dApar/dt = sum_s q_s int dv vpar d/dt(F_s) / ( (k_perp^2/mu_0 + sum_s q_s^2/m_s int dv F_s) )
  gkyl_dg_div_op_range(
    field->div_mem, &app->basis, 0, out, 0, field->currentDensdot, 0,
    field->dApartdtSlvr_lhs_factor, &app->local
  );

  app->stat.field_apar_solve_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_field_ampere_solve_1x_enabled(
  gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *out
)
{
  struct timespec wst = gkyl_wall_clock();

  // Weak division method Apar = sum_s q_s int dv vpar F_s / (k_perp^2/mu_0)
  gkyl_dg_div_op_range(
    field->div_mem, &app->basis, 0, out, 0, field->currentDens, 0, field->lapWeightAmpere,
    &app->local
  );

  // Smooth Apar after solving Ampere's law.
  gk_field_fem_projection_par(app, field, out, out);

  app->stat.field_apar_solve_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_field_ampere_solve_1x_none(
  gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *out
)
{
  // Do nothing.
}

static void
gk_field_em_rhs_enabled(
  gkyl_gyrokinetic_app *app, struct gk_field *field, const struct gkyl_array *f_in[],
  struct gkyl_array *rhs_in[]
)
{
  gk_field_accumulate_current_dens_dot(app, field, rhs_in);
  gk_field_accumulate_ohms_kSq(app, field, f_in);
  gk_field_ohm_solve_1x(app, field, field->apardot);
}

static void
gk_field_em_rhs_none(
  gkyl_gyrokinetic_app *app, struct gk_field *field, const struct gkyl_array *f_in[],
  struct gkyl_array *rhs_in[]
)
{
  // Do nothing.
}

static void
gk_field_fem_projection_par_none(
  gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *arr_dg,
  struct gkyl_array *arr_fem
)
{
  // Do nothing.
}

static void
gk_field_fem_release_1x(const gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  gkyl_array_release(gkf->rho_c);
  gkyl_array_release(gkf->rho_c_global_dg);
  gkyl_array_release(gkf->rho_c_global_smooth);
  gkyl_array_release(gkf->phi_smooth);
  gkyl_array_release(gkf->phi_fem);
  gkyl_array_release(gkf->apar);
  gkyl_array_release(gkf->apardot);

  if (gkf->is_em) {
    gkyl_array_release(gkf->apar_curr);
    gkyl_array_release(gkf->apar1);
    gkyl_array_release(gkf->aparnew);
    gkyl_array_release(gkf->amperesol);
    gkyl_array_release(gkf->currentDens);
    gkyl_array_release(gkf->currentDensdot);
    gkyl_array_release(gkf->currentDens_global);
    gkyl_array_release(gkf->dApartdtSlvr_lhs_factor);
    gkyl_array_release(gkf->lapWeightAmpere);
    gkyl_array_release(gkf->dApartdtSlvr_kSq);
    gkyl_dg_bin_op_mem_release(gkf->div_mem);
    if (app->use_gpu) {
      gkyl_array_release(gkf->apar_host);
      gkyl_array_release(gkf->apardot_host);
      gkyl_array_release(gkf->amperesol_host);
    }
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
  gkf->apar = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  gkf->apardot = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

  // Allocate additional arrays for EM case.
  if (gkf->is_em) {
    // Arrays for A_parallel RK stages and Ampere's law solution.
    gkf->apar_curr = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    gkf->apar1 = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    gkf->aparnew = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    gkf->amperesol = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

    gkf->apar_host = gkf->apar;
    gkf->apardot_host = gkf->apardot;
    gkf->amperesol_host = gkf->amperesol;
    if (app->use_gpu) {
      gkf->apar_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
      gkf->apardot_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
      gkf->amperesol_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
    }

    gkf->currentDens = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    gkf->currentDensdot = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    gkf->lapWeightAmpere =
      mkarr(app->use_gpu, (2 * (app->cdim / 3) + 1) * app->basis.num_basis, app->local_ext.volume);
    gkf->dApartdtSlvr_kSq = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

    gkf->currentDens_global = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);
    gkf->dApartdtSlvr_lhs_factor = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    gkf->div_mem = app->use_gpu ?
                   gkyl_dg_bin_op_mem_cu_dev_new(app->local.volume, app->basis.num_basis) :
                   gkyl_dg_bin_op_mem_new(app->local.volume, app->basis.num_basis);

    // Need to set weight to kperpsq*polarizationWeight for use in potential smoothing.
    gkyl_array_copy(gkf->lapWeightAmpere, app->gk_geom->geo_int.jacobgeo);
    gkyl_array_scale(gkf->lapWeightAmpere, gkf->info.kperpSq / gkf->info.mu0);

    gkf->fem_parproj_ampere_bc = GKYL_FEM_PARPROJ_NONE;
    for (int d = 0; d < app->num_periodic_dir; ++d) {
      if (app->periodic_dirs[d] == app->cdim - 1) {
        gkf->fem_parproj_ampere_bc = GKYL_FEM_PARPROJ_PERIODIC;
      }
    }

    assert(gkf->info.mu0 > 0.0);
    gkf->apar_energy_fac_1d = 0.5 / gkf->info.mu0 * gkf->info.kperpSq;
    gkf->apar_energy_fac_1d *= -1.0; // the -1 is to be consistent with ES energy.
  }

  // Allocate phi_host for I/O.
  gkf->phi_host = gkf->phi_smooth;
  if (app->use_gpu) {
    gkf->phi_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
  }

  // Set function pointers.
  gkf->rhs_phi_func = gk_field_rhs_phi_1x;
  if (gkf->is_em) {
    gkf->em_rhs_func = gkf->info.is_apar_static ? gk_field_em_rhs_none : gk_field_em_rhs_enabled;
    gkf->ampere_solve = gk_field_ampere_solve_1x_enabled;
    gkf->fem_projection_par_apar_func = gk_field_fem_projection_par_none;
  } else {
    gkf->em_rhs_func = gk_field_em_rhs_none;
    gkf->ampere_solve = gk_field_ampere_solve_1x_none;
  }

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
