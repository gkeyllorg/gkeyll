#include <gkyl_moment_priv.h>

// initialize source solver: this should be called after all species
// and fields are initialized
void
moment_coupling_init(const struct gkyl_moment_app *app, struct moment_coupling *src)
{
  struct gkyl_moment_em_coupling_inp src_inp = {
    .grid = &app->grid,
    .nfluids = app->num_species,
    // if there is a field, need to update electric field too, otherwise just updating fluid
    .epsilon0 = app->field.epsilon0 ? app->field.epsilon0 : 0.0, 
    .mu0 = app->field.mu0 ? app->field.mu0 : 0.0, 
    // is the field static?
    .static_field = app->field.is_static, 
    // linear ramping function for slowing turning on applied accelerations, E fields, or currents
    .t_ramp_E = app->field.t_ramp_E ? app->field.t_ramp_E : 0.0,
    .t_ramp_curr = app->field.t_ramp_curr ? app->field.t_ramp_curr : 0.0,
  };

  for (int i=0; i<app->num_species; ++i)
    src_inp.param[i] = (struct gkyl_moment_em_coupling_data) {
      .type = app->species[i].eqn_type,
      .charge = app->species[i].charge,
      .mass = app->species[i].mass,
      // The gradient-based closure defines its heat flux through k0, so k0=0.0 in the source solve to avoid double-applying it.
      // The neural-network closure supplies the heat flux directly, so k0 is retained here as the integrating-factor relaxation rate.
      .k0 = (app->species[i].has_grad_closure) ? 0.0 : app->species[i].k0,
    };

  src_inp.has_collision = app->has_collision;
  for (int s=0; s<app->num_species; ++s)
    for (int r=0; r<app->num_species; ++r)
      src_inp.nu_base[s][r] = app->nu_base[s][r];

  src_inp.has_nT_sources = false;
  for (int i=0; i<app->num_species; ++i)
    if (app->species[i].proj_nT_source)
      src_inp.has_nT_sources = true;

  // check for relativistic-species
  bool use_rel = 0;
  for (int n=0;  n<app->num_species; ++n) 
    use_rel = use_rel || (app->species[n].eqn_type == GKYL_EQN_COLDFLUID_SR) || app->field.use_explicit_em_coupling;

  src_inp.has_frictional_sources = false;
  src_inp.use_explicit_friction = false;
  for (int i = 0; i < app->num_species; i++) {
    if(app->species[i].has_friction) {
      src_inp.has_frictional_sources = true;

      if (app->species[i].friction_Z != 0.0 || app->species[i].friction_T_elc != 0.0 || app->species[i].friction_Lambda_ee != 0.0) {
        src_inp.friction_Z = app->species[i].friction_Z;
        src_inp.friction_T_elc = app->species[i].friction_T_elc;
        src_inp.friction_Lambda_ee = app->species[i].friction_Lambda_ee;
      }
    }

    if(app->species[i].use_explicit_friction) {
      src_inp.use_explicit_friction = true;
    }
  }

  src_inp.has_volume_sources = false;
  for (int i = 0; i < app->num_species; i++) {
    if (app->species[i].has_volume_sources) {
      src_inp.has_volume_sources = true;

      if(app->species[i].volume_gas_gamma != 0.0 || app->species[i].volume_U0 != 0.0 || app->species[i].volume_R0 != 0.0) {
        src_inp.volume_gas_gamma = app->species[i].volume_gas_gamma;
        src_inp.volume_U0 = app->species[i].volume_U0;
        src_inp.volume_R0 = app->species[i].volume_R0;
      }
    }
  }
  if(app->field.has_volume_sources) {
    src_inp.has_volume_sources = true;

    if(app->field.volume_gas_gamma != 0.0 || app->field.volume_U0 != 0.0 || app->field.volume_R0 != 0.0) {
      src_inp.volume_gas_gamma = app->field.volume_gas_gamma;
      src_inp.volume_U0 = app->field.volume_U0;
      src_inp.volume_R0 = app->field.volume_R0;
    }
  }

  src_inp.has_reactive_sources = false;
  for (int i = 0; i < app->num_species; i++) {
    if (app->species[i].has_reactivity) {
      src_inp.has_reactive_sources = true;

      if (app->species[i].reactivity_gas_gamma != 0.0 || app->species[i].reactivity_specific_heat_capacity != 0.0 ||
        app->species[i].reactivity_energy_of_formation != 0.0 || app->species[i].reactivity_ignition_temperature != 0.0 ||
        app->species[i].reactivity_reaction_rate != 0.0) {

        src_inp.reactivity_gas_gamma = app->species[i].reactivity_gas_gamma;
        src_inp.reactivity_specific_heat_capacity = app->species[i].reactivity_specific_heat_capacity;
        src_inp.reactivity_energy_of_formation = app->species[i].reactivity_energy_of_formation;
        src_inp.reactivity_ignition_temperature = app->species[i].reactivity_ignition_temperature;
        src_inp.reactivity_reaction_rate = app->species[i].reactivity_reaction_rate;
      }
    }
  }
  
  src_inp.has_einstein_medium_sources = false;
  for (int i = 0; i < app->num_species; i++) {
    if (app->species[i].has_einstein_medium) {
      src_inp.has_einstein_medium_sources = true;

      if (app->species[i].medium_gas_gamma != 0.0 || app->species[i].medium_kappa != 0.0) {
        src_inp.medium_gas_gamma = app->species[i].medium_gas_gamma;
        src_inp.medium_kappa = app->species[i].medium_kappa;
      }
    }
  }

  src_inp.has_gr_ultra_rel_sources = false;
  for (int i = 0; i < app->num_species; i++) {
    if (app->species[i].has_gr_ultra_rel) {
      src_inp.has_gr_ultra_rel_sources = true;

      if (app->species[i].gr_ultra_rel_gas_gamma != 0.0) {
        src_inp.gr_ultra_rel_gas_gamma = app->species[i].gr_ultra_rel_gas_gamma;
      }
    }
  }

  src_inp.has_gr_euler_sources = false;
  for (int i = 0; i < app->num_species; i++) {
    if (app->species[i].has_gr_euler) {
      src_inp.has_gr_euler_sources = true;

      if (app->species[i].gr_euler_gas_gamma != 0.0) {
        src_inp.gr_euler_gas_gamma = app->species[i].gr_euler_gas_gamma;
      }
    }
  }

  src_inp.has_gr_twofluid_sources = false;
  for (int i = 0; i < app->num_species; i++) {
    if (app->species[i].has_gr_twofluid) {
      src_inp.has_gr_twofluid_sources = true;

      if (app->species[i].gr_twofluid_mass_elc != 0.0 || app->species[i].gr_twofluid_mass_ion != 0.0 || app->species[i].gr_twofluid_charge_elc != 0.0 ||
        app->species[i].gr_twofluid_charge_ion != 0.0 || app->species[i].gr_twofluid_gas_gamma_elc != 0.0 || app->species[i].gr_twofluid_gas_gamma_ion != 0.0 ||
        app->species[i].gr_twofluid_e_fact != 0.0) {
        src_inp.gr_twofluid_mass_elc = app->species[i].gr_twofluid_mass_elc;
        src_inp.gr_twofluid_mass_ion = app->species[i].gr_twofluid_mass_ion;
        src_inp.gr_twofluid_charge_elc = app->species[i].gr_twofluid_charge_elc;
        src_inp.gr_twofluid_charge_ion = app->species[i].gr_twofluid_charge_ion;
        src_inp.gr_twofluid_gas_gamma_elc = app->species[i].gr_twofluid_gas_gamma_elc;
        src_inp.gr_twofluid_gas_gamma_ion = app->species[i].gr_twofluid_gas_gamma_ion;
        src_inp.gr_twofluid_e_fact = app->species[i].gr_twofluid_e_fact;
      }
    }
  }

  src_inp.has_gr_mhd_sources = false;
  for (int i = 0; i < app->num_species; i++) {
    if (app->species[i].has_gr_mhd) {
      src_inp.has_gr_mhd_sources = true;

      if (app->species[i].gr_mhd_gas_gamma != 0.0) {
        src_inp.gr_mhd_gas_gamma = app->species[i].gr_mhd_gas_gamma;
      }
    }
  }

  src_inp.has_vacuum_einstein_sources = false;
  for (int i = 0; i < app->num_species; i++) {
    if (app->species[i].has_vacuum_einstein) {
      src_inp.has_vacuum_einstein_sources = true;

      src_inp.vacuum_einstein_excision_threshold = app->species[i].vacuum_einstein_excision_threshold;
      src_inp.vacuum_einstein_spacetime_slicing = app->species[i].vacuum_einstein_spacetime_slicing;
      src_inp.vacuum_einstein_spacetime_evolution = app->species[i].vacuum_einstein_spacetime_evolution;
    }
  }

  src_inp.has_vacuum_einstein_conformal_sources = false;
  for (int i = 0; i < app->num_species; i++) {
    if (app->species[i].has_vacuum_einstein_conformal) {
      src_inp.has_vacuum_einstein_conformal_sources = true;

      src_inp.vacuum_einstein_conformal_excision_threshold = app->species[i].vacuum_einstein_conformal_excision_threshold;
      src_inp.vacuum_einstein_conformal_spacetime_slicing = app->species[i].vacuum_einstein_conformal_spacetime_slicing;
      src_inp.vacuum_einstein_conformal_spacetime_evolution = app->species[i].vacuum_einstein_conformal_spacetime_evolution;
    }
  }

  // save the use-rel bool
  src_inp.use_rel = use_rel;

  // check for explicit em-coupling
  src_inp.use_explicit_em_coupling = 0;
  if (app->field.use_explicit_em_coupling)
    src_inp.use_explicit_em_coupling = 1;

  // create updater to solve for sources
  src->slvr = gkyl_moment_em_coupling_new(src_inp);

  for (int n=0; n<app->num_species; ++n) {
    int meqn = app->species[n].num_equations;
    src->pr_rhs[n] = mkarr(false, meqn, app->local_ext.volume);
    src->non_ideal_cflrate[n] = mkarr(false, 1, app->local_ext.volume);
  }

  int ghost[3] = { 1, 1, 1 };
  // create non-ideal local extended range from local range
  // has one additional cell in each direction because non-ideal
  // variables are stored at cell vertices
  gkyl_create_vertex_ranges(&app->local, ghost, &src->non_ideal_local_ext,
    &src->non_ideal_local);

  // check if gradient-closure is present
  for (int i=0; i<app->num_species; ++i) {
    if (app->species[i].eqn_type == GKYL_EQN_TEN_MOMENT && app->species[i].has_grad_closure) {
      int nadj[3] = { 1, 4, 8 }; // cells adjacent to a vertex
      src->non_ideal_vars[i] = mkarr(false, nadj[app->ndim - 1]*10,
        src->non_ideal_local_ext.volume);
      struct gkyl_ten_moment_grad_closure_inp grad_closure_inp = {
        .grid = &app->grid,
        .k0 = app->species[i].k0,
        .cfl = app->cfl,
        .comm = app->comm,
        .update_range = &app->local,
        .heat_flux_range = &src->non_ideal_local,
      };
      src->grad_closure_slvr[i] = gkyl_ten_moment_grad_closure_new(&grad_closure_inp);
    }
  }

  // Check whether neural network-based closure is present.
  for (int i = 0; i < app->num_species; i++) {
    if (app->species[i].eqn_type == GKYL_EQN_TEN_MOMENT && app->species[i].has_nn_closure) {
      struct gkyl_ten_moment_nn_closure_inp nn_closure_inp = {
        .grid = &app->grid,
        .k0 = app->species[i].k0,
        .poly_order = app->species[i].poly_order,
        .ann = app->species[i].ann,
      };
      src->nn_closure_slvr[i] = gkyl_ten_moment_nn_closure_new(nn_closure_inp);
    }
  }

  // check if braginskii terms are present
  if (app->has_braginskii) {
    struct gkyl_moment_braginskii_inp brag_inp = {
      .grid = &app->grid,
      .nfluids = app->num_species,
      .epsilon0 = app->field.epsilon0,
      // Check for multiplicative collisionality factor, default is 1.0
      .coll_fac = app->coll_fac == 0 ? 1.0 : app->coll_fac,
    };
    for (int i=0; i<app->num_species; ++i) {
      // Braginskii coefficients depend on pressure and coefficient to obtain
      // pressure is different for different equation systems (gasGamma, vt, Tr(P))
      double p_fac = 1.0;
      if (app->species[i].eqn_type == GKYL_EQN_EULER) {
        p_fac =  gkyl_wv_euler_gas_gamma(app->species[i].equation);
      }
      else if (app->species[i].eqn_type == GKYL_EQN_ISO_EULER) {
        p_fac =  gkyl_wv_iso_euler_vt(app->species[i].equation);
      }
      brag_inp.param[i] = (struct gkyl_moment_braginskii_data) {
        .type_eqn = app->species[i].eqn_type,
        .type_brag = app->species[i].type_brag,
        .charge = app->species[i].charge,
        .mass = app->species[i].mass,
        .p_fac = p_fac,
      };
    }
    src->brag_slvr = gkyl_moment_braginskii_new(brag_inp);
  }
}

// update sources: 'nstrang' is 0 for the first Strang step and 1 for
// the second step
struct gkyl_update_status
moment_coupling_update(gkyl_moment_app *app, struct moment_coupling *src,
  int nstrang, double tcurr, double dt)
{
  int sidx[] = { 0, app->ndim };
  struct gkyl_array *fluids[GKYL_MAX_SPECIES];
  const struct gkyl_array *app_accels[GKYL_MAX_SPECIES];
  const struct gkyl_array *pr_rhs_const[GKYL_MAX_SPECIES];
  const struct gkyl_array *nT_sources[GKYL_MAX_SPECIES];

  double dt_suggested = DBL_MAX;
  struct gkyl_ten_moment_grad_closure_status stat;

  for (int i=0; i<app->num_species; ++i) {
    fluids[i] = app->species[i].f[sidx[nstrang]];

    if (app->species[i].app_accel_evolve) {
      gkyl_fv_proj_advance(app->species[i].app_accel_proj, tcurr, &app->local, app->species[i].app_accel);
    }
    app_accels[i] = app->species[i].app_accel;

    if (app->species[i].eqn_type == GKYL_EQN_TEN_MOMENT && app->species[i].has_grad_closure) {
      // Non-ideal variables are defined on an extended range with one additional "cell" in each direction.
      // This additional cell accounts for the fact that non-ideal variables are stored at cell vertices.
      stat = gkyl_ten_moment_grad_closure_advance(src->grad_closure_slvr[i],
        &src->non_ideal_local, &app->local,
        app->species[i].f[sidx[nstrang]], app->field.f[sidx[nstrang]],
        src->non_ideal_cflrate[i], dt, src->non_ideal_vars[i], src->pr_rhs[i]);

      if (!stat.success)
        return (struct gkyl_update_status) {
          .success = false,
          .dt_suggested = stat.dt_suggested
        };

      dt_suggested = fmin(dt_suggested, stat.dt_suggested);
    }

    if (app->species[i].eqn_type == GKYL_EQN_TEN_MOMENT && app->species[i].has_nn_closure) {
      // Non-ideal variables are defined on an extended range with one additional "cell" in each direction.
      // This additional cell accounts for the fact that non-ideal variables are stored at cell vertices.
      gkyl_ten_moment_nn_closure_advance(src->nn_closure_slvr[i], &src->non_ideal_local, &app->local, app->species[i].f[sidx[nstrang]],
        app->field.f[sidx[nstrang]], src->non_ideal_vars[i], src->pr_rhs[i]);
    }
  }

  if (app->has_braginskii) {
    gkyl_moment_braginskii_advance(src->brag_slvr,
      src->non_ideal_local, app->local,
      fluids, app->field.f[sidx[nstrang]],
      src->non_ideal_cflrate, src->non_ideal_vars, src->pr_rhs);
  }

  if (app->field.ext_em_evolve) {
    gkyl_fv_proj_advance(app->field.ext_em_proj, tcurr, &app->local, app->field.ext_em);
  }

  if (app->field.app_current_evolve) {
    if (app->field.use_explicit_em_coupling) {
      gkyl_fv_proj_advance(app->field.app_current_proj, tcurr, &app->local, app->field.app_current);
      gkyl_fv_proj_advance(app->field.app_current_proj, tcurr + dt*2, &app->local, app->field.app_current1);
      gkyl_fv_proj_advance(app->field.app_current_proj, tcurr + 2*dt/2.0, &app->local, app->field.app_current2);
    }
    else {
      gkyl_fv_proj_advance(app->field.app_current_proj, tcurr, &app->local, app->field.app_current);
    }
  }

  // Get the RHS pointer for accumulation during source update
  for (int i=0; i<app->num_species; ++i) {
    pr_rhs_const[i] = src->pr_rhs[i];
  }

  for (int i=0; i<app->num_species; ++i) {
    if (app->species[i].proj_nT_source
        && !(app->species[i].nT_source_set_only_once
             && app->species[i].nT_source_is_set))
    {
      gkyl_fv_proj_advance(app->species[i].proj_nT_source, tcurr,
                           &app->local, app->species[i].nT_source);
    }
    nT_sources[i] = app->species[i].nT_source;
    app->species[i].nT_source_is_set = true;
  }

  if (app->field.use_explicit_em_coupling) {
    gkyl_moment_em_coupling_explicit_advance(src->slvr, tcurr, dt, &app->local,
      fluids, app_accels, pr_rhs_const, 
      app->field.f[sidx[nstrang]], app->field.app_current, app->field.app_current1,
      app->field.app_current2, app->field.ext_em, 
      nT_sources, app->field.app_current_proj,nstrang);
  }
  else {
    gkyl_moment_em_coupling_implicit_advance(src->slvr, tcurr, dt, &app->local,
      fluids, app_accels, pr_rhs_const, 
      app->field.f[sidx[nstrang]], app->field.app_current, app->field.ext_em, 
      nT_sources);
  }

  for (int i=0; i<app->num_species; ++i) {
    moment_species_apply_bc(app, tcurr, &app->species[i], fluids[i]);
  }
  if (app->has_field) {
    moment_field_apply_bc(app, tcurr, &app->field, app->field.f[sidx[nstrang]]);
  }
  return (struct gkyl_update_status) {
    .success = true,
    .dt_suggested = dt_suggested
  };
}

// free sources
void
moment_coupling_release(const struct gkyl_moment_app *app, const struct moment_coupling *src)
{
  gkyl_moment_em_coupling_release(src->slvr);
  for (int i=0; i<app->num_species; ++i) {
    gkyl_array_release(src->pr_rhs[i]);
    gkyl_array_release(src->non_ideal_cflrate[i]);
    gkyl_array_release(src->non_ideal_vars[i]);
    if (app->species[i].eqn_type == GKYL_EQN_TEN_MOMENT && app->species[i].has_grad_closure) {
      gkyl_ten_moment_grad_closure_release(src->grad_closure_slvr[i]);
    }
    if (app->species[i].eqn_type == GKYL_EQN_TEN_MOMENT && app->species[i].has_nn_closure) {
      gkyl_ten_moment_nn_closure_release(src->nn_closure_slvr[i]);
    }
  }
  if (app->has_braginskii)
    gkyl_moment_braginskii_release(src->brag_slvr);
}

