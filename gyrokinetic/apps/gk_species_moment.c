#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_gk_fluid_moments.h>

// Whether the moment is computed by the Maxwellian (or BiMaxwellian) moments updater.
static bool
gk_species_moment_uses_maxwellian_moms(const struct gk_species_moment *sm)
{
  return sm->type == GKYL_F_MOMENT_MAXWELLIAN || sm->type == GKYL_F_MOMENT_BIMAXWELLIAN;
}

static void
gk_species_moment_diag_jacobgeo_div_disabled(
  const struct gkyl_gyrokinetic_app *app, struct gk_species_moment *sm, struct gkyl_array *Jmom_in,
  struct gkyl_array *mom_out
)
{
  // Do nothing.
}

static void
gk_species_moment_diag_jacobgeo_div_enabled_1st_comp(
  const struct gkyl_gyrokinetic_app *app, struct gk_species_moment *sm, struct gkyl_array *Jmom_in,
  struct gkyl_array *mom_out
)
{
  // Only divide the first component.
  gkyl_dg_div_op_range(
    sm->mem_geo, &app->basis, 0, mom_out, 0, Jmom_in, 0, app->gk_geom->geo_int.jacobgeo, &app->local
  );
}

static void
gk_species_moment_diag_jacobgeo_div_enabled_all_comp(
  const struct gkyl_gyrokinetic_app *app, struct gk_species_moment *sm, struct gkyl_array *Jmom_in,
  struct gkyl_array *mom_out
)
{
  // Divide all components.
  for (int k = 0; k < sm->num_mom; k++) {
    gkyl_dg_div_op_range(
      sm->mem_geo, &app->basis, k, mom_out, k, Jmom_in, 0, app->gk_geom->geo_int.jacobgeo,
      &app->local
    );
  }
}

static void
gk_species_moment_calc_kinetic(
  const struct gkyl_gyrokinetic_app *app, const struct gk_species_moment *sm,
  const struct gkyl_range phase_rng, const struct gkyl_range conf_rng, const struct gkyl_array *fin
)
{
  if (sm->is_integrated) {
    gkyl_mom_gyrokinetic_advance(
      sm->mcalc, &phase_rng, &conf_rng, sm->species->gyro_phi, fin, sm->marr
    );
  } else {
    if (sm->type == GKYL_F_MOMENT_MAXWELLIAN) {
      gkyl_gk_maxwellian_moments_advance(
        sm->gyrokinetic_maxwellian_moms, &phase_rng, &conf_rng, fin, sm->marr
      );
    } else if (sm->type == GKYL_F_MOMENT_BIMAXWELLIAN) {
      gkyl_gk_bimaxwellian_moments_advance(
        sm->gyrokinetic_maxwellian_moms, &phase_rng, &conf_rng, fin, sm->marr
      );
    } else {
      gkyl_mom_gyrokinetic_advance(
        sm->mcalc, &phase_rng, &conf_rng, sm->species->gyro_phi, fin, sm->marr
      );
    }
  }
}

static void
gk_species_moment_release_kinetic(
  const struct gkyl_gyrokinetic_app *app, const struct gk_species_moment *sm
)
{
  gkyl_array_release(sm->marr);
  if (app->use_gpu) {
    gkyl_array_release(sm->marr_host);
  }

  if (sm->is_integrated) {
    gkyl_mom_gyrokinetic_release(sm->mcalc);
  } else {
    if (gk_species_moment_uses_maxwellian_moms(sm)) {
      gkyl_gk_maxwellian_moments_release(sm->gyrokinetic_maxwellian_moms);
    } else {
      gkyl_mom_gyrokinetic_release(sm->mcalc);
    }

    // Free the weak division memory.
    gkyl_dg_bin_op_mem_release(sm->mem_geo);
  }
}

static void
gk_species_moment_calc_fluid(
  const struct gkyl_gyrokinetic_app *app, const struct gk_species_moment *sm,
  const struct gkyl_range species_rng, const struct gkyl_range conf_rng,
  const struct gkyl_array *fin
)
{
  gkyl_gk_fluid_moments_advance(sm->gf_moms, &conf_rng, fin, app->field->phi_smooth, sm->marr);
}

static void
gk_species_moment_release_fluid(
  const struct gkyl_gyrokinetic_app *app, const struct gk_species_moment *sm
)
{
  gkyl_gk_fluid_moments_release(sm->gf_moms);
  gkyl_array_release(sm->marr);
  if (app->use_gpu) {
    gkyl_array_release(sm->marr_host);
  }
  if (!sm->is_integrated) {
    gkyl_dg_bin_op_mem_release(sm->mem_geo);
  }
}

/**
 * Compute a weighted kinetic moment with the current field.
 *
 * @param app Application providing the current potential.
 * @param sm Initialized species moment object.
 * @param phase_rng Phase-space selection.
 * @param conf_rng Configuration selection.
 * @param fin Input distribution.
 * @param weight_range Independently indexed weight storage range.
 * @param weight Current weight array.
 */
static void
gk_species_moment_calc_kinetic_weighted(
  const struct gkyl_gyrokinetic_app *app, const struct gk_species_moment *sm,
  const struct gkyl_range phase_rng, const struct gkyl_range conf_rng, const struct gkyl_array *fin,
  const struct gkyl_range *weight_range, const struct gkyl_array *weight
)
{
  gkyl_mom_weighted_gyrokinetic_advance(
    sm->mcalc_weighted, &phase_rng, &conf_rng, weight_range, weight, sm->species->gyro_phi, fin,
    sm->marr
  );
}

/**
 * Compute a configuration-weighted raw fluid moment.
 *
 * @param app Application providing the current potential.
 * @param sm Initialized species moment object.
 * @param phase_rng Species range (unused for fluid moments).
 * @param conf_rng Configuration selection.
 * @param fin Distribution or conservative fluid state.
 * @param weight_range Independently indexed weight storage range.
 * @param weight Current weight array.
 */
static void
gk_species_moment_calc_fluid_weighted(
  const struct gkyl_gyrokinetic_app *app, const struct gk_species_moment *sm,
  const struct gkyl_range phase_rng, const struct gkyl_range conf_rng, const struct gkyl_array *fin,
  const struct gkyl_range *weight_range, const struct gkyl_array *weight
)
{
  gkyl_gk_fluid_moments_advance_weighted(
    sm->gf_moms, &conf_rng, fin, app->field->phi_smooth, weight_range, weight, sm->marr
  );
}

/**
 * Release the weighted updater and its output/workspace arrays.
 *
 * @param app Application specifying the backend.
 * @param sm Initialized weighted species moment object.
 */
static void
gk_species_moment_release_weighted(
  const struct gkyl_gyrokinetic_app *app, const struct gk_species_moment *sm
)
{
  gkyl_mom_weighted_gyrokinetic_release(sm->mcalc_weighted);
  gkyl_array_release(sm->marr);
  if (app->use_gpu) {
    gkyl_array_release(sm->marr_host);
  }
  if (!sm->is_integrated) {
    gkyl_dg_bin_op_mem_release(sm->mem_geo);
  }
}

static void
gk_species_moment_init_fluid(
  struct gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gk_species_moment *sm,
  enum gkyl_distribution_moments mom_type, bool weighted
)
{
  sm->gf_moms = gkyl_gk_fluid_moments_inew(&(struct gkyl_gk_fluid_moments_inp){
    .basis = &app->basis,
    .grid = &app->grid,
    .range = &app->local,
    .range_ext = &app->local_ext,
    .bmag = app->gk_geom->geo_int.bmag,
    .mass = gks->info.mass,
    .charge = gks->info.charge,
    .num_moments = gks->num_moments,
    .moment = mom_type,
    .integrated = sm->is_integrated,
    .weighted = weighted,
    .use_gpu = app->use_gpu,
  });

  sm->num_mom = gkyl_gk_fluid_moments_num_mom(sm->gf_moms);

  sm->marr = mkarr(
    app->use_gpu, sm->num_mom * (sm->is_integrated ? 1 : app->basis.num_basis),
    app->local_ext.volume
  );
  sm->marr_host = app->use_gpu ? mkarr(false, sm->marr->ncomp, sm->marr->size) : sm->marr;

  if (!sm->is_integrated) {
    sm->mem_geo = app->use_gpu ?
                    gkyl_dg_bin_op_mem_cu_dev_new(app->local.volume, app->basis.num_basis) :
                    gkyl_dg_bin_op_mem_new(app->local.volume, app->basis.num_basis);
    sm->diag_jacobgeo_div_func = gk_species_moment_uses_maxwellian_moms(sm) ?
                                   gk_species_moment_diag_jacobgeo_div_enabled_1st_comp :
                                   gk_species_moment_diag_jacobgeo_div_enabled_all_comp;
  }

  sm->calc_func = weighted ? NULL : gk_species_moment_calc_fluid;
  sm->weighted_calc_func = weighted ? gk_species_moment_calc_fluid_weighted : NULL;
  sm->release_func = gk_species_moment_release_fluid;
}

static void
gk_species_moment_init_kinetic(
  struct gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gk_species_moment *sm,
  enum gkyl_distribution_moments mom_type, enum gkyl_mom_weight_type weight_type
)
{
  sm->calc_func = gk_species_moment_calc_kinetic;
  sm->release_func = gk_species_moment_release_kinetic;

  // The gyroaveraged potential of this species enters the Hamiltonian moments.
  sm->species = gks;

  if (weight_type != GKYL_F_MOMENT_WEIGHT_NONE) {
    sm->mcalc_weighted = gkyl_mom_weighted_gyrokinetic_new(
      gks->info.mass, gks->info.charge, &app->basis, &gks->basis, &gks->grid, gks->vel_map,
      app->gk_geom, mom_type, weight_type, sm->is_integrated, app->use_gpu
    );
    sm->mcalc = NULL;
    sm->num_mom = gkyl_mom_weighted_gyrokinetic_num_mom(sm->mcalc_weighted);
    sm->calc_func = NULL;
    sm->weighted_calc_func = gk_species_moment_calc_kinetic_weighted;
    sm->release_func = gk_species_moment_release_weighted;
  } else if (!sm->is_integrated && gk_species_moment_uses_maxwellian_moms(sm)) {
    struct gkyl_gk_maxwellian_moments_inp inp_mom = {
      .phase_grid = &gks->grid,
      .conf_basis = &app->basis,
      .phase_basis = &gks->basis,
      .conf_range = &app->local,
      .conf_range_ext = &app->local_ext,
      .gk_geom = app->gk_geom,
      .vel_map = gks->vel_map,
      .divide_jacobgeo = false,
      .mass = gks->info.mass,
      .use_gpu = app->use_gpu,
    };
    sm->gyrokinetic_maxwellian_moms = gkyl_gk_maxwellian_moments_inew(&inp_mom);
    if (sm->type == GKYL_F_MOMENT_MAXWELLIAN) {
      sm->num_mom = 3; // (n, u_par, T/m).
    } else {
      sm->num_mom = 4; // (n, u_par, T_par/m, T_perp/m).
    }
  } else {
    sm->mcalc = gkyl_mom_gyrokinetic_new(
      gks->info.mass, gks->info.charge, &app->basis, &gks->basis, &gks->grid, gks->vel_map,
      app->gk_geom, mom_type, sm->is_integrated, app->use_gpu
    );
    sm->num_mom = gkyl_mom_gyrokinetic_num_mom(sm->mcalc);
  }

  // Allocate output arrays for the selected updater.
  sm->marr = mkarr(
    app->use_gpu, sm->num_mom * (sm->is_integrated ? 1 : app->basis.num_basis),
    app->local_ext.volume
  );
  sm->marr_host = app->use_gpu ? mkarr(false, sm->marr->ncomp, sm->marr->size) : sm->marr;

  if (!sm->is_integrated) {
    sm->mem_geo = app->use_gpu ?
                    gkyl_dg_bin_op_mem_cu_dev_new(app->local.volume, app->basis.num_basis) :
                    gkyl_dg_bin_op_mem_new(app->local.volume, app->basis.num_basis);
    sm->diag_jacobgeo_div_func = gk_species_moment_uses_maxwellian_moms(sm) ?
                                   gk_species_moment_diag_jacobgeo_div_enabled_1st_comp :
                                   gk_species_moment_diag_jacobgeo_div_enabled_all_comp;
  }
}

void
gk_species_moment_init(
  struct gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gk_species_moment *sm,
  enum gkyl_distribution_moments mom_type, enum gkyl_mom_weight_type weight_type, bool is_integrated
)
{
  sm->is_integrated = is_integrated;
  sm->type = mom_type;
  sm->diag_jacobgeo_div_func = gk_species_moment_diag_jacobgeo_div_disabled;

  sm->weighted_calc_func = NULL;
  sm->mcalc_weighted = NULL;
  if (weight_type != GKYL_F_MOMENT_WEIGHT_NONE) {
    if (weight_type != GKYL_F_MOMENT_WEIGHT_CONF && weight_type != GKYL_F_MOMENT_WEIGHT_VEL) {
      gkyl_exit("Phase-space moment weights are not implemented.");
    }
    if (gks->is_fluid && weight_type != GKYL_F_MOMENT_WEIGHT_CONF) {
      gkyl_exit("Fluid moments only support configuration-space weights.");
    }
  }

  if (gks->is_fluid) {
    gk_species_moment_init_fluid(app, gks, sm, mom_type, weight_type != GKYL_F_MOMENT_WEIGHT_NONE);
  } else {
    gk_species_moment_init_kinetic(app, gks, sm, mom_type, weight_type);
  }
}

void
gk_species_moment_diag_jacobgeo_div(
  const struct gkyl_gyrokinetic_app *app, struct gk_species_moment *sm, struct gkyl_array *Jmom_in,
  struct gkyl_array *mom_out
)
{
  sm->diag_jacobgeo_div_func(app, sm, Jmom_in, mom_out);
}

void
gk_species_moment_calc(
  const struct gkyl_gyrokinetic_app *app, const struct gk_species_moment *sm,
  const struct gkyl_range species_rng, const struct gkyl_range conf_rng,
  const struct gkyl_array *fin
)
{
  assert(sm->calc_func);
  sm->calc_func(app, sm, species_rng, conf_rng, fin);
}

void
gk_species_moment_calc_weighted(
  const struct gkyl_gyrokinetic_app *app, const struct gk_species_moment *sm,
  const struct gkyl_range species_rng, const struct gkyl_range conf_rng,
  const struct gkyl_array *fin, const struct gkyl_range *weight_range,
  const struct gkyl_array *weight
)
{
  assert(sm->weighted_calc_func && weight_range && weight);
  sm->weighted_calc_func(app, sm, species_rng, conf_rng, fin, weight_range, weight);
}

void
gk_species_moment_release(const struct gkyl_gyrokinetic_app *app, const struct gk_species_moment *sm)
{
  sm->release_func(app, sm);
}
