#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

static void
gk_species_moment_diag_jacobgeo_div_disabled(const struct gkyl_gyrokinetic_app *app,
  struct gk_species_moment *sm, struct gkyl_array *Jmom_in, struct gkyl_array *mom_out)
{
  // Do nothing.
}

static void
gk_species_moment_diag_jacobgeo_div_enabled_1st_comp(const struct gkyl_gyrokinetic_app *app,
  struct gk_species_moment *sm, struct gkyl_array *Jmom_in, struct gkyl_array *mom_out)
{
  // Only divide the first component.
  gkyl_dg_div_op_range(sm->mem_geo, &app->basis, 0, mom_out, 0, Jmom_in, 0, 
    app->gk_geom->geo_int.jacobgeo, &app->local);  
}

static void
gk_species_moment_diag_jacobgeo_div_enabled_all_comp(const struct gkyl_gyrokinetic_app *app,
  struct gk_species_moment *sm, struct gkyl_array *Jmom_in, struct gkyl_array *mom_out)
{
  // Divide all components.
  for (int k=0; k<sm->num_mom; k++)
    gkyl_dg_div_op_range(sm->mem_geo, &app->basis, k, mom_out, k, Jmom_in, 0,
      app->gk_geom->geo_int.jacobgeo, &app->local);  
}

void
gk_species_moment_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s,
  struct gk_species_moment *sm, enum gkyl_distribution_moments mom_type, bool is_integrated)
{
  sm->is_integrated = is_integrated;
  sm->is_maxwellian_moms = mom_type == GKYL_F_MOMENT_MAXWELLIAN;
  sm->is_bimaxwellian_moms = mom_type == GKYL_F_MOMENT_BIMAXWELLIAN;
  sm->diag_jacobgeo_div_func = gk_species_moment_diag_jacobgeo_div_disabled;

  if (sm->is_integrated) {
    // Create moment operator.
    sm->mcalc = gkyl_dg_updater_moment_gyrokinetic_new(&s->grid, &app->basis, 
      &s->basis, &app->local_ext, s->info.mass, s->info.charge, s->vel_map, app->gk_geom,
      app->field->phi_smooth, mom_type, sm->is_integrated, app->use_gpu);    

    sm->num_mom = gkyl_dg_updater_moment_gyrokinetic_num_mom(sm->mcalc);

    // Allocate arrays to hold moments.
    sm->marr = mkarr(app->use_gpu, sm->num_mom, app->local_ext.volume);
    sm->marr_host = sm->marr;
    if (app->use_gpu) {
      sm->marr_host = mkarr(false, sm->num_mom, app->local_ext.volume);  
    }
  }
  else {
    // Create moment operator.
    if (sm->is_maxwellian_moms || sm->is_bimaxwellian_moms) {
      struct gkyl_gk_maxwellian_moments_inp inp_mom = {
        .phase_grid = &s->grid,
        .conf_basis = &app->basis,
        .phase_basis = &s->basis,
        .conf_range =  &app->local,
        .conf_range_ext = &app->local_ext,
        .gk_geom = app->gk_geom,
        .vel_map = s->vel_map,
        .divide_jacobgeo = false, 
        .mass = s->info.mass,
        .use_gpu = app->use_gpu,
      };
      sm->gyrokinetic_maxwellian_moms = gkyl_gk_maxwellian_moments_inew(  &inp_mom  );
      if (sm->is_maxwellian_moms) {
        sm->num_mom = 3; // (n, u_par, T/m)
      }
      else {
        sm->num_mom = 4; // (n, u_par, T_par/m, T_perp/m)
      }  
      sm->diag_jacobgeo_div_func = gk_species_moment_diag_jacobgeo_div_enabled_1st_comp;
    }
    else {
      sm->mcalc = gkyl_dg_updater_moment_gyrokinetic_new(&s->grid, &app->basis, 
        &s->basis, &app->local_ext, s->info.mass, s->info.charge, s->vel_map, app->gk_geom,
        app->field->phi_smooth, mom_type, sm->is_integrated, app->use_gpu);    

      sm->num_mom = gkyl_dg_updater_moment_gyrokinetic_num_mom(sm->mcalc);
      sm->diag_jacobgeo_div_func = gk_species_moment_diag_jacobgeo_div_enabled_all_comp;
    }

    // Allocate arrays to hold moments.
    sm->marr = mkarr(app->use_gpu, sm->num_mom*app->basis.num_basis, app->local_ext.volume);
    sm->marr_host = sm->marr;
    if (app->use_gpu) {
      sm->marr_host = mkarr(false, sm->num_mom*app->basis.num_basis, app->local_ext.volume);
    }
    // Bin Op memory for rescaling moment by inverse of Jacobian
    if (app->use_gpu) {
      sm->mem_geo = gkyl_dg_bin_op_mem_cu_dev_new(app->local.volume, app->basis.num_basis);
    }
    else {
      sm->mem_geo = gkyl_dg_bin_op_mem_new(app->local.volume, app->basis.num_basis);
    }
  }
}

void
gk_species_moment_calc(const struct gk_species_moment *sm,
  const struct gkyl_range phase_rng, const struct gkyl_range conf_rng,
  const struct gkyl_array *fin)
{
  if (sm->is_integrated) {
    gkyl_dg_updater_moment_gyrokinetic_advance(sm->mcalc, 
      &phase_rng, &conf_rng, fin, sm->marr);
  }
  else {
    if (sm->is_maxwellian_moms) {
      gkyl_gk_maxwellian_moments_advance(sm->gyrokinetic_maxwellian_moms, 
        &phase_rng, &conf_rng, fin, sm->marr);
    }
    else if (sm->is_bimaxwellian_moms) {
      gkyl_gk_bimaxwellian_moments_advance(sm->gyrokinetic_maxwellian_moms, 
        &phase_rng, &conf_rng, fin, sm->marr);
    } 
    else {
      gkyl_dg_updater_moment_gyrokinetic_advance(sm->mcalc, 
        &phase_rng, &conf_rng, fin, sm->marr);
    }
  }
}

void
gk_species_moment_diag_jacobgeo_div(const struct gkyl_gyrokinetic_app *app,
  struct gk_species_moment *sm, struct gkyl_array *Jmom_in, struct gkyl_array *mom_out)
{
  sm->diag_jacobgeo_div_func(app, sm, Jmom_in, mom_out);
}

void
gk_species_moment_release(const struct gkyl_gyrokinetic_app *app, const struct gk_species_moment *sm)
{
  gkyl_array_release(sm->marr);
  if (app->use_gpu)
    gkyl_array_release(sm->marr_host);

  if (sm->is_integrated) {
    gkyl_dg_updater_moment_gyrokinetic_release(sm->mcalc);
  }
  else {
    if (sm->is_maxwellian_moms || sm->is_bimaxwellian_moms) {
      gkyl_gk_maxwellian_moments_release(sm->gyrokinetic_maxwellian_moms);
    }
    else {
      gkyl_dg_updater_moment_gyrokinetic_release(sm->mcalc);
    }

    // Free the weak division memory.
    gkyl_dg_bin_op_mem_release(sm->mem_geo);
  }
}
