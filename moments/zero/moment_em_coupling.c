#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_fv_proj.h>
#include <gkyl_moment_em_coupling.h>
#include <gkyl_moment_em_coupling_priv.h>
#include <gkyl_mat.h>

gkyl_moment_em_coupling*
gkyl_moment_em_coupling_new(struct gkyl_moment_em_coupling_inp inp)
{
  gkyl_moment_em_coupling *mom_em = gkyl_malloc(sizeof(gkyl_moment_em_coupling));

  mom_em->grid = *(inp.grid);
  mom_em->ndim = mom_em->grid.ndim;
  mom_em->nfluids = inp.nfluids;
  for (int i = 0; i < mom_em->nfluids; i++) {
    mom_em->param[i] = inp.param[i];
  }

  mom_em->epsilon0 = inp.epsilon0;
  mom_em->mu0 = inp.mu0;

  mom_em->is_charged_species = false;
  for (int i = 0; i < mom_em->nfluids; i++) {
    if (inp.param[i].charge != 0.0) {
      mom_em->is_charged_species = true;
    }
  }

  mom_em->static_field = inp.static_field; 
  mom_em->t_ramp_E = inp.t_ramp_E;
  if (mom_em->t_ramp_E != 0.0) {
    mom_em->ramp_app_E = true;
  }
  else {
    mom_em->ramp_app_E = false;
  }
  mom_em->t_ramp_curr = inp.t_ramp_curr;
  if (mom_em->t_ramp_curr != 0.0) {
    mom_em->ramp_app_curr = true;
  }
  else {
    mom_em->ramp_app_curr = false;
  }

  mom_em->has_collision = inp.has_collision;
  mom_em->use_rel = inp.use_rel;

  if (mom_em->has_collision) {
    for (int i = 0; i < mom_em->nfluids; i++) {
      for (int j = 0; j < mom_em->nfluids; j++) {
        mom_em->nu_base[i][j] = inp.nu_base[i][j];
      }
    }
  }
  else {
    for (int i = 0; i < mom_em->nfluids; i++) {
      for (int j = 0; j < mom_em->nfluids; j++) {
        mom_em->nu_base[i][j] = 0.0;
      }
    }
  }

  mom_em->use_explicit_em_coupling = inp.use_explicit_em_coupling;

  mom_em->has_nT_sources = inp.has_nT_sources;

  mom_em->has_frictional_sources = inp.has_frictional_sources;
  if (mom_em->has_frictional_sources) {
    mom_em->use_explicit_friction = inp.use_explicit_friction;
    mom_em->friction_Z = inp.friction_Z;
    mom_em->friction_T_elc = inp.friction_T_elc;
    mom_em->friction_Lambda_ee = inp.friction_Lambda_ee;
  }

  mom_em->has_volume_sources = inp.has_volume_sources;
  if (mom_em->has_volume_sources) {
    mom_em->volume_gas_gamma = inp.volume_gas_gamma;
    mom_em->volume_U0 = inp.volume_U0;
    mom_em->volume_R0 = inp.volume_R0;
  }

  mom_em->has_reactive_sources = inp.has_reactive_sources;
  if (mom_em->has_reactive_sources) {
    mom_em->reactivity_gas_gamma = inp.reactivity_gas_gamma;
    mom_em->reactivity_specific_heat_capacity = inp.reactivity_specific_heat_capacity;
    mom_em->reactivity_energy_of_formation = inp.reactivity_energy_of_formation;
    mom_em->reactivity_ignition_temperature = inp.reactivity_ignition_temperature;
    mom_em->reactivity_reaction_rate = inp.reactivity_reaction_rate;
  }

  mom_em->has_einstein_medium_sources = inp.has_einstein_medium_sources;
  if (mom_em->has_einstein_medium_sources) {
    mom_em->medium_gas_gamma = inp.medium_gas_gamma;
    mom_em->medium_kappa = inp.medium_kappa;
  }

  mom_em->has_gr_ultra_rel_sources = inp.has_gr_ultra_rel_sources;
  if (mom_em->has_gr_ultra_rel_sources) {
    mom_em->gr_ultra_rel_gas_gamma = inp.gr_ultra_rel_gas_gamma;
  }

  mom_em->has_gr_euler_sources = inp.has_gr_euler_sources;
  if (mom_em->has_gr_euler_sources) {
    mom_em->gr_euler_gas_gamma = inp.gr_euler_gas_gamma;
  }

  mom_em->has_gr_twofluid_sources = inp.has_gr_twofluid_sources;
  if (mom_em->has_gr_twofluid_sources) {
    mom_em->gr_twofluid_mass_elc = inp.gr_twofluid_mass_elc;
    mom_em->gr_twofluid_mass_ion = inp.gr_twofluid_mass_ion;
    mom_em->gr_twofluid_charge_elc = inp.gr_twofluid_charge_elc;
    mom_em->gr_twofluid_charge_ion = inp.gr_twofluid_charge_ion;
    mom_em->gr_twofluid_gas_gamma_elc = inp.gr_twofluid_gas_gamma_elc;
    mom_em->gr_twofluid_gas_gamma_ion = inp.gr_twofluid_gas_gamma_ion;
    mom_em->gr_twofluid_e_fact = inp.gr_twofluid_e_fact;
  }

  mom_em->has_vacuum_einstein_sources = inp.has_vacuum_einstein_sources;
  if (mom_em->has_vacuum_einstein_sources) {
    mom_em->vacuum_einstein_excision_threshold = inp.vacuum_einstein_excision_threshold;
    mom_em->vacuum_einstein_spacetime_slicing = inp.vacuum_einstein_spacetime_slicing;
    mom_em->vacuum_einstein_spacetime_evolution = inp.vacuum_einstein_spacetime_evolution;
  }

  mom_em->has_vacuum_einstein_conformal_sources = inp.has_vacuum_einstein_conformal_sources;
  if (mom_em->has_vacuum_einstein_conformal_sources) {
    mom_em->vacuum_einstein_conformal_excision_threshold = inp.vacuum_einstein_conformal_excision_threshold;
    mom_em->vacuum_einstein_conformal_spacetime_slicing = inp.vacuum_einstein_conformal_spacetime_slicing;
    mom_em->vacuum_einstein_conformal_spacetime_evolution = inp.vacuum_einstein_conformal_spacetime_evolution;
  }
  
  mom_em->has_gr_mhd_sources = inp.has_gr_mhd_sources;
  if (mom_em->has_gr_mhd_sources) {
    mom_em->gr_mhd_gas_gamma = inp.gr_mhd_gas_gamma;
  }

  return mom_em;
}

void
gkyl_moment_em_coupling_implicit_advance(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, const struct gkyl_range* update_range,
  struct gkyl_array* fluid[GKYL_MAX_SPECIES], const struct gkyl_array* app_accel[GKYL_MAX_SPECIES], const struct gkyl_array* p_rhs[GKYL_MAX_SPECIES],
  struct gkyl_array* em, const struct gkyl_array* app_current, const struct gkyl_array* ext_em, const struct gkyl_array* nT_sources[GKYL_MAX_SPECIES])
{
  int nfluids = mom_em->nfluids;
  double *fluid_s[GKYL_MAX_SPECIES];
  const double *app_accel_s[GKYL_MAX_SPECIES];
  const double *p_rhs_s[GKYL_MAX_SPECIES];
  const double *nT_sources_s[GKYL_MAX_SPECIES];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);

  while (gkyl_range_iter_next(&iter)) {
    long cell_idx = gkyl_range_idx(update_range, iter.idx);

    for (int i = 0; i < nfluids; i++) {
      fluid_s[i] = gkyl_array_fetch(fluid[i], cell_idx);
      app_accel_s[i] = gkyl_array_cfetch(app_accel[i], cell_idx);
      p_rhs_s[i] = gkyl_array_cfetch(p_rhs[i], cell_idx);
      nT_sources_s[i] = gkyl_array_cfetch(nT_sources[i], cell_idx);
    }

    double *em_arr = em ? gkyl_array_fetch(em, cell_idx) : 0;
    const double *app_current_arr = app_current ? gkyl_array_cfetch(app_current, cell_idx) : 0;
    const double *ext_em_arr = ext_em ? gkyl_array_cfetch(ext_em, cell_idx) : 0;

    implicit_source_coupling_update(mom_em, t_curr, dt, fluid_s, app_accel_s, p_rhs_s, em_arr, app_current_arr, ext_em_arr, nT_sources_s);
  }
}

void
gkyl_moment_em_coupling_explicit_advance(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, const struct gkyl_range* update_range,
  struct gkyl_array* fluid[GKYL_MAX_SPECIES], const struct gkyl_array* app_accel[GKYL_MAX_SPECIES], const struct gkyl_array* p_rhs[GKYL_MAX_SPECIES],
  struct gkyl_array* em, const struct gkyl_array* app_current, const struct gkyl_array* app_current1, const struct gkyl_array* app_current2,
  const struct gkyl_array* ext_em, const struct gkyl_array* nT_sources[GKYL_MAX_SPECIES], gkyl_fv_proj *proj_app_curr, int nstrang)
{
  int nfluids = mom_em->nfluids;
  double *fluid_s[GKYL_MAX_SPECIES];
  const double *app_accel_s[GKYL_MAX_SPECIES];
  const double *p_rhs_s[GKYL_MAX_SPECIES];
  const double *nT_sources_s[GKYL_MAX_SPECIES];

  double dt_local = 2.0 * dt;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);

  while (gkyl_range_iter_next(&iter)) {
    long cell_idx = gkyl_range_idx(update_range, iter.idx);

    for (int i = 0; i < nfluids; i++) {
      fluid_s[i] = gkyl_array_fetch(fluid[i], cell_idx);
      app_accel_s[i] = gkyl_array_cfetch(app_accel[i], cell_idx);
    }

    double *em_arr = em ? em_arr = gkyl_array_fetch(em, cell_idx) : 0;
    const double *app_current_arr = app_current ? gkyl_array_cfetch(app_current, cell_idx) : 0;
    const double *app_current1_arr = app_current1 ? gkyl_array_cfetch(app_current1, cell_idx) : 0;
    const double *app_current2_arr = app_current2 ? gkyl_array_cfetch(app_current2, cell_idx) : 0;
    const double *ext_em_arr = ext_em ? gkyl_array_cfetch(ext_em, cell_idx) : 0;

    if (mom_em->use_rel) {
      explicit_source_coupling_update(mom_em, t_curr, dt_local, fluid_s, app_accel_s, em_arr, app_current_arr, app_current1_arr, app_current2_arr,
        ext_em_arr, nstrang);
    }
  }
}

void
gkyl_moment_em_coupling_release(gkyl_moment_em_coupling* mom_em)
{
  gkyl_free(mom_em);
}