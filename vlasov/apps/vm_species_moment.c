#include <assert.h>
#include <gkyl_mom_canonical_pb.h>
#include <gkyl_vlasov_priv.h>

// Initialize species moment object.
void
vm_species_moment_init(struct gkyl_vlasov_app *app, struct vm_species *s,
  struct vm_species_moment *sm, enum gkyl_distribution_moments mom_type, bool is_integrated)
{
  sm->is_integrated = is_integrated;

  int num_mom;
  sm->is_vlasov_lte_moms = mom_type == GKYL_F_MOMENT_LTE;
  if (sm->is_vlasov_lte_moms) {
    struct gkyl_vlasov_lte_moments_inp inp_mom = {
      .phase_grid = &s->grid,
      .vel_grid = &s->grid_vel, 
      .conf_basis = &app->confBasis,
      .vel_basis = &app->velBasis, 
      .phase_basis = &app->basis,
      .conf_range =  &app->local,
      .conf_range_ext = &app->local_ext,
      .vel_range = &s->local_vel,
      .phase_range = &s->local,
      .gamma = s->gamma,
      .gamma_inv = s->gamma_inv,
      .h_ij = s->h_ij,
      .h_ij_inv = s->h_ij_inv,
      .det_h = s->det_h,
      .hamil = s->hamil,
      .model_id = s->model_id,
      .use_gpu = app->use_gpu,
    };
    // Compute (n, V_drift, T/m)
    sm->vlasov_lte_moms = gkyl_vlasov_lte_moments_inew(&inp_mom);
    num_mom = app->vdim + 2;
  }
  else {
    if (s->model_id == GKYL_MODEL_SR) {
      struct gkyl_mom_vlasov_sr_auxfields sr_inp = {.gamma = s->gamma};
      sm->mcalc = gkyl_dg_updater_moment_new(&s->grid, &app->confBasis, 
        &app->basis, &app->local, &s->local_vel, &s->local, s->model_id, &sr_inp, 
        mom_type, is_integrated, app->use_gpu);
      num_mom = gkyl_dg_updater_moment_num_mom(sm->mcalc);
    } else if ((s->model_id == GKYL_MODEL_CANONICAL_PB || s->model_id == GKYL_MODEL_CANONICAL_PB_GR)
      && (mom_type == GKYL_F_MOMENT_M1_FROM_H || mom_type == GKYL_F_MOMENT_ENERGY
      || (sm->is_integrated && mom_type == GKYL_F_MOMENT_M0M1M2))) {
      struct gkyl_mom_canonical_pb_auxfields can_pb_inp = {.hamil = s->hamil};
      sm->mcalc = gkyl_dg_updater_moment_new(&s->grid, &app->confBasis, 
        &app->basis, &app->local, &s->local_vel, &s->local, s->model_id, &can_pb_inp, 
        mom_type, is_integrated, app->use_gpu);
      num_mom = gkyl_dg_updater_moment_num_mom(sm->mcalc);
    }  
    else {
      // No auxiliary fields for moments if not SR 
      sm->mcalc = gkyl_dg_updater_moment_new(&s->grid, &app->confBasis, 
        &app->basis, &app->local, &s->local_vel, &s->local, s->model_id, 0, 
        mom_type, is_integrated, app->use_gpu);   
      num_mom = gkyl_dg_updater_moment_num_mom(sm->mcalc); 
    }
  }

  if (is_integrated) {
    sm->marr = mkarr(app->use_gpu, num_mom, app->local_ext.volume);
    sm->marr_host = sm->marr;
    if (app->use_gpu) {
      sm->marr_host = mkarr(false, num_mom, app->local_ext.volume);      
    }
  }
  else {
    sm->marr = mkarr(app->use_gpu, num_mom*app->confBasis.num_basis,
      app->local_ext.volume);
    sm->marr_host = sm->marr;
    if (app->use_gpu) {
      sm->marr_host = mkarr(false, num_mom*app->confBasis.num_basis,
        app->local_ext.volume);
    }
  }
}

void
vm_species_moment_calc(const struct vm_species_moment *sm,
  const struct gkyl_range phase_rng, const struct gkyl_range conf_rng,
  const struct gkyl_array *fin)
{
  if (sm->is_vlasov_lte_moms) {
    gkyl_vlasov_lte_moments_advance(sm->vlasov_lte_moms, 
      &phase_rng, &conf_rng, fin, sm->marr);
  }
  else {
    gkyl_dg_updater_moment_advance(sm->mcalc, 
      &phase_rng, &conf_rng, fin, sm->marr);
  }
}

// release memory for moment data object
void
vm_species_moment_release(const struct gkyl_vlasov_app *app, const struct vm_species_moment *sm)
{
  if (app->use_gpu) {
    gkyl_array_release(sm->marr_host);
  }
  gkyl_array_release(sm->marr);

  if(sm->is_vlasov_lte_moms) {
    gkyl_vlasov_lte_moments_release(sm->vlasov_lte_moms);
  }
  else {
    gkyl_dg_updater_moment_release(sm->mcalc);
  }
}
