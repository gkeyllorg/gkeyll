#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

static void
gk_neut_species_moment_diag_jacobgeo_div_disabled(const struct gkyl_gyrokinetic_app *app,
  struct gk_species_moment *sm, struct gkyl_array *Jmom_in, struct gkyl_array *mom_out)
{
  // Do nothing.
}

static void
gk_neut_species_moment_diag_jacobgeo_div_enabled_1st_comp(const struct gkyl_gyrokinetic_app *app,
  struct gk_species_moment *sm, struct gkyl_array *Jmom_in, struct gkyl_array *mom_out)
{
  // Only divide the first component.
  gkyl_dg_div_op_range(sm->mem_geo, &app->basis, 0, mom_out, 0, Jmom_in, 0, 
    app->gk_geom->geo_int.jacobgeo, &app->local);
}

static void
gk_neut_species_moment_diag_jacobgeo_div_enabled_all_comp(const struct gkyl_gyrokinetic_app *app,
  struct gk_species_moment *sm, struct gkyl_array *Jmom_in, struct gkyl_array *mom_out)
{
  // Divide all components.
  for (int k=0; k<sm->num_mom; k++)
    gkyl_dg_div_op_range(sm->mem_geo, &app->basis, k, mom_out, k, Jmom_in, 0,
      app->gk_geom->geo_int.jacobgeo, &app->local);  
}

static void
gk_neut_species_kinetic_moment_calc(const struct gk_species_moment *sm,
  const struct gkyl_range phase_rng, const struct gkyl_range conf_rng,
  const struct gkyl_array *fin)
{
  if (sm->is_maxwellian_moms) {
    gkyl_vlasov_lte_moments_advance(sm->vlasov_lte_moms, 
      &phase_rng, &conf_rng, fin, sm->marr);
  }
  else {
    gkyl_mom_calc_advance(sm->mom_calc,
      &phase_rng, &conf_rng, fin, sm->marr);
  }
}

static void
gk_neut_species_kinetic_moment_release(const struct gkyl_gyrokinetic_app *app, const struct gk_species_moment *sm)
{
  gkyl_array_release(sm->marr);
  if (app->use_gpu)
    gkyl_array_release(sm->marr_host);

  if (sm->is_integrated) {
    gkyl_mom_type_release(sm->mom_type);
    gkyl_mom_calc_release(sm->mom_calc);
  }
  else {
    if (sm->is_maxwellian_moms) {
      gkyl_vlasov_lte_moments_release(sm->vlasov_lte_moms);
    }
    else {
      gkyl_mom_type_release(sm->mom_type);
      gkyl_mom_calc_release(sm->mom_calc);
    }

    // Free the weak division memory.
    gkyl_dg_bin_op_mem_release(sm->mem_geo);
  }
}

static void
gk_neut_species_kinetic_moment_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s,
  struct gk_species_moment *sm, enum gkyl_distribution_moments mom_type, bool is_integrated)
{
  sm->diag_jacobgeo_div_func = gk_neut_species_moment_diag_jacobgeo_div_disabled;

  // Initialize kinetic neutral species moment object.
  if (sm->is_integrated) {
    // Create moment operator.
    struct gkyl_mom_vlasov_inp inp_mom = {
      .conf_basis = &app->basis,
      .phase_basis = &s->basis,
      .vel_range = &s->local_vel,
      .vel_map = s->vlasov_vel_map,
      .hamil_range = &s->hamil_range,
      .hamil = s->hamil,
      .model_id = s->model_id,
      .hamil_id = s->hamil_id,
      .mom_type = mom_type,
      .use_gpu = app->use_gpu,
    };
    sm->mom_type = gkyl_int_mom_vlasov_inew(&inp_mom);
    sm->mom_calc = gkyl_mom_calc_new(&s->grid, sm->mom_type, app->use_gpu);

    sm->num_mom = gkyl_mom_type_num_mom(sm->mom_type);

    // Allocate arrays to hold moments.
    sm->marr = mkarr(app->use_gpu, sm->num_mom, app->local_ext.volume);
    sm->marr_host = sm->marr;
    if (app->use_gpu) {
      sm->marr_host = mkarr(false, sm->num_mom, app->local_ext.volume); 
    }
  }
  else {
    // Create moment operator.
    if (sm->is_maxwellian_moms) {
      struct gkyl_vlasov_lte_moments_inp inp_mom = {
        .phase_grid = &s->grid,
        .vel_grid = &s->grid_vel,
        .conf_basis = &app->basis,
        .vel_basis = &s->basis_vel,
        .phase_basis = &s->basis,
        .conf_range =  &app->local,
        .conf_range_ext = &app->local_ext,
        .vel_range = &s->local_vel,
        .phase_range = &s->local,
        .hamil_range = &s->hamil_range,
        .hamil = s->hamil,
        .model_id = s->model_id,
        .hamil_id = s->hamil_id,
        .vel_map = s->vlasov_vel_map,
        .use_extended_hamil_def = false,
        .use_gpu = app->use_gpu,
      };
      sm->vlasov_lte_moms = gkyl_vlasov_lte_moments_inew(&inp_mom);
      sm->num_mom = 5; // (n, ux, uy, uz, T/m).
      sm->diag_jacobgeo_div_func = gk_neut_species_moment_diag_jacobgeo_div_enabled_1st_comp;
    }
    else {
      struct gkyl_mom_vlasov_inp inp_mom = {
        .conf_basis = &app->basis,
        .phase_basis = &s->basis,
        .vel_range = &s->local_vel,
        .vel_map = s->vlasov_vel_map,
        .hamil_range = &s->hamil_range,
        .hamil = s->hamil,
        .model_id = s->model_id,
        .hamil_id = s->hamil_id,
        .mom_type = mom_type,
        .use_gpu = app->use_gpu,
      };
      sm->mom_type = gkyl_mom_vlasov_inew(&inp_mom);
      sm->mom_calc = gkyl_mom_calc_new(&s->grid, sm->mom_type, app->use_gpu);

      sm->num_mom = gkyl_mom_type_num_mom(sm->mom_type);
      sm->diag_jacobgeo_div_func = gk_neut_species_moment_diag_jacobgeo_div_enabled_all_comp;
    }

    // Allocate arrays to hold moments.
    sm->marr = mkarr(app->use_gpu, sm->num_mom*app->basis.num_basis, app->local_ext.volume);
    sm->marr_host = sm->marr;
    if (app->use_gpu)
      sm->marr_host = mkarr(false, sm->num_mom*app->basis.num_basis, app->local_ext.volume);

    // Bin Op memory for rescaling moment by inverse of Jacobian
    if (app->use_gpu) {
      sm->mem_geo = gkyl_dg_bin_op_mem_cu_dev_new(app->local.volume, app->basis.num_basis);
    }
    else {
      sm->mem_geo = gkyl_dg_bin_op_mem_new(app->local.volume, app->basis.num_basis);
    }
  }

  sm->calc_func = gk_neut_species_kinetic_moment_calc;
  sm->release_func = gk_neut_species_kinetic_moment_release;
}

static void
gk_neut_species_fluid_moment_calc_m0(const struct gk_species_moment *sm,
  const struct gkyl_range phase_rng, const struct gkyl_range conf_rng,
  const struct gkyl_array *fin)
{
  gkyl_array_set_offset(sm->marr, 1.0/sm->mass, fin, 0);
}

static void
gk_neut_species_fluid_moment_calc_m1(const struct gk_species_moment *sm,
  const struct gkyl_range phase_rng, const struct gkyl_range conf_rng,
  const struct gkyl_array *fin)
{
  gkyl_array_set_offset(sm->marr, 1.0/sm->mass, fin, 1*sm->num_basis_conf);
}

static void
gk_neut_species_fluid_moment_calc_m2(const struct gk_species_moment *sm,
  const struct gkyl_range phase_rng, const struct gkyl_range conf_rng,
  const struct gkyl_array *fin)
{
  gkyl_array_set_offset(sm->marr, 2.0/sm->mass, fin, 4*sm->num_basis_conf);
}

static void
gk_neut_species_fluid_moment_calc(const struct gk_species_moment *sm,
  const struct gkyl_range phase_rng, const struct gkyl_range conf_rng,
  const struct gkyl_array *fin)
{
  if (sm->is_integrated) {
    gkyl_gk_neut_fluid_prim_vars_mass_momentum_flow_thermal_energy_advance(sm->nf_prim_vars, fin, sm->marr, 0);
  }
  else {
    if (sm->is_maxwellian_moms) {
      // Get LTE moments.
      gkyl_gk_neut_fluid_prim_vars_lte_advance(sm->nf_prim_vars, fin, sm->marr, 0);
    }
    else {
      sm->fluid_calc_M(sm, phase_rng, conf_rng, fin);
    }  
  }  
}

static void
gk_neut_species_fluid_moment_release(const struct gkyl_gyrokinetic_app *app, const struct gk_species_moment *sm)
{
  gkyl_array_release(sm->marr);
  if (app->use_gpu)
    gkyl_array_release(sm->marr_host);

  if (sm->is_integrated) {
    gkyl_gk_neut_fluid_prim_vars_release(sm->nf_prim_vars);
  }
  else {
    if (sm->is_maxwellian_moms) {
      gkyl_gk_neut_fluid_prim_vars_release(sm->nf_prim_vars);
    }
    else {
      // Nothing to release.
    }

    // Free the weak division memory.
    gkyl_dg_bin_op_mem_release(sm->mem_geo);
  }
}

static void
gk_neut_species_fluid_moment_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s,
  struct gk_species_moment *sm, enum gkyl_distribution_moments mom_type, bool is_integrated)
{
  sm->diag_jacobgeo_div_func = gk_neut_species_moment_diag_jacobgeo_div_disabled;

  // Initialize fluid neutral species moment object.
  if (sm->is_integrated) {
    sm->num_mom = 6; // rho, rho*ux, rho*uy, rho*uz, flowE, thermalE.

    // Allocate arrays to hold moments.
    sm->marr = mkarr(app->use_gpu, sm->num_mom, app->local_ext.volume);
    sm->marr_host = sm->marr;
    if (app->use_gpu) {
      sm->marr_host = mkarr(false, sm->num_mom, app->local_ext.volume); 
    }

    sm->nf_prim_vars = gkyl_gk_neut_fluid_prim_vars_new(s->info.gas_gamma, s->info.mass, &app->basis,
      &app->grid, &app->local_ext, GKYL_GK_NEUT_FLUID_PRIM_VARS_MASS_MOMENTUM_FLOW_THERMAL_ENERGY,
      true, app->use_gpu);
  }
  else {
    if (sm->is_maxwellian_moms) {
      // Compute (n, ux, uy, uz, T/m) moments.
      sm->num_mom = 5;
      sm->nf_prim_vars = gkyl_gk_neut_fluid_prim_vars_new(s->info.gas_gamma, s->info.mass,
        &app->basis, &app->grid, &app->local_ext, GKYL_GK_NEUT_FLUID_PRIM_VARS_LTE, false, app->use_gpu);
      sm->diag_jacobgeo_div_func = gk_neut_species_moment_diag_jacobgeo_div_enabled_1st_comp;
    }
    else {
      sm->mass = s->info.mass;
      sm->num_basis_conf = app->basis.num_basis;

      if (mom_type == GKYL_F_MOMENT_M0) {
        sm->num_mom = 1;
        sm->fluid_calc_M = gk_neut_species_fluid_moment_calc_m0;
      }
      else if (mom_type == GKYL_F_MOMENT_M1) {
        sm->num_mom = 3;
        sm->fluid_calc_M = gk_neut_species_fluid_moment_calc_m1;
      }
      else if (mom_type == GKYL_F_MOMENT_M2) {
        sm->num_mom = 1;
        sm->fluid_calc_M = gk_neut_species_fluid_moment_calc_m2;
      }
      else {
        // Not yet implemented.
        assert(false);
      }

      sm->diag_jacobgeo_div_func = gk_neut_species_moment_diag_jacobgeo_div_enabled_all_comp;
    }

    // Allocate arrays to hold moments.
    sm->marr = mkarr(app->use_gpu, sm->num_mom*app->basis.num_basis, app->local_ext.volume);
    sm->marr_host = sm->marr;
    if (app->use_gpu)
      sm->marr_host = mkarr(false, sm->marr->ncomp, sm->marr->size);

    // Bin Op memory for rescaling moment by inverse of Jacobian
    if (app->use_gpu) {
      sm->mem_geo = gkyl_dg_bin_op_mem_cu_dev_new(app->local.volume, app->basis.num_basis);
    }
    else {
      sm->mem_geo = gkyl_dg_bin_op_mem_new(app->local.volume, app->basis.num_basis);
    }
  }

  sm->calc_func = gk_neut_species_fluid_moment_calc;
  sm->release_func = gk_neut_species_fluid_moment_release;
}

void
gk_neut_species_moment_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s,
  struct gk_species_moment *sm, enum gkyl_distribution_moments mom_type, bool is_integrated)
{
  // Initialize neutral species moment object.
  sm->is_integrated = is_integrated;
  sm->is_maxwellian_moms = mom_type == GKYL_F_MOMENT_LTE;

  if (s->is_fluid) {
    gk_neut_species_fluid_moment_init(app, s, sm, mom_type, is_integrated);
  }
  else {
    gk_neut_species_kinetic_moment_init(app, s, sm, mom_type, is_integrated);
  }
}

void
gk_neut_species_moment_calc(const struct gk_species_moment *sm,
  const struct gkyl_range phase_rng, const struct gkyl_range conf_rng,
  const struct gkyl_array *fin)
{
  sm->calc_func(sm, phase_rng, conf_rng, fin);
}

void
gk_neut_species_moment_diag_jacobgeo_div(const struct gkyl_gyrokinetic_app *app,
  struct gk_species_moment *sm, struct gkyl_array *Jmom_in, struct gkyl_array *mom_out)
{
  sm->diag_jacobgeo_div_func(app, sm, Jmom_in, mom_out);
}

void
gk_neut_species_moment_release(const struct gkyl_gyrokinetic_app *app, const struct gk_species_moment *sm)
{
  sm->release_func(app, sm);
}
