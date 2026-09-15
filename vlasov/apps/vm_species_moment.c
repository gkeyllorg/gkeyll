#include <assert.h>
#include <gkyl_vlasov_priv.h>

// Initialize species moment object.
void
vm_species_moment_init(struct gkyl_vlasov_app *app, struct vm_species *vms,
  struct vm_species_moment *sm, enum gkyl_distribution_moments mom_type, bool is_integrated)
{
  sm->is_integrated = is_integrated;
  sm->is_vlasov_lte_moms = mom_type == GKYL_F_MOMENT_LTE;

  sm->mom_hamil = gkyl_array_acquire(vms->mom_hamil);
  sm->mom_hamil_range = &vms->mom_hamil_range;

  int num_mom;
  if (sm->is_vlasov_lte_moms) {
    struct gkyl_vlasov_lte_moments_inp inp_mom = {
      .phase_grid = &vms->grid,
      .vel_grid = &vms->grid_vel, 
      .conf_basis = &app->basis,
      .vel_basis = &vms->basis_vel, 
      .phase_basis = &vms->basis,
      .conf_range =  &app->local,
      .conf_range_ext = &app->local_ext,
      .vel_range = &vms->local_vel,
      .phase_range = &vms->local,
      .vel_map = vms->vel_map,
      .hamil_range = sm->mom_hamil_range,
      .hamil = sm->mom_hamil,
      .model_id = vms->model_id,
      .hamil_id = vms->mom_hamil_id,
      .gamma_inv = vms->gamma_inv,
      .h_ij = vms->h_ij,
      .h_ij_inv = vms->h_ij_inv,
      .det_h = vms->det_h,
      .use_extended_hamil_def = vms->info.use_extended_hamil_def,
      .effective_potential = vms->effective_potential,
      .use_gpu = app->use_gpu,
    };
    // Compute (n, V_drift, T/m)
    sm->vlasov_lte_moms = gkyl_vlasov_lte_moments_inew(&inp_mom);
    num_mom = app->vdim + 2;
    sm->marr = mkarr(app->use_gpu, num_mom*app->basis.num_basis,
      app->local_ext.volume); 
  }
  else {
    struct gkyl_mom_vlasov_inp inp_mom = {
      .conf_basis = &app->basis,
      .phase_basis = &vms->basis,
      .vel_range = &vms->local_vel,
      .vel_map = vms->vel_map,
      .hamil_range = &vms->mom_hamil_range,
      .hamil = vms->mom_hamil,
      .model_id = vms->model_id,
      .hamil_id = vms->mom_hamil_id,
      .mom_type = mom_type, 
      .use_gpu = app->use_gpu,
    };
    if (is_integrated) {
      sm->mom_type = gkyl_int_mom_vlasov_inew(&inp_mom);
      num_mom = gkyl_mom_type_num_mom(sm->mom_type); 
      sm->marr = mkarr(app->use_gpu, num_mom, app->local_ext.volume);       
    }
    else {
      sm->mom_type = gkyl_mom_vlasov_inew(&inp_mom);
      num_mom = gkyl_mom_type_num_mom(sm->mom_type); 
      sm->marr = mkarr(app->use_gpu, num_mom*app->basis.num_basis,
        app->local_ext.volume); 
    }
    sm->mom_calc = gkyl_mom_calc_new(&vms->grid, sm->mom_type, app->use_gpu);
  }  
  
  // Host moment for  I/O.
  sm->marr_host = app->use_gpu ? mkarr(false, sm->marr->ncomp, sm->marr->size)
                               : gkyl_array_acquire(sm->marr);
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
    gkyl_mom_calc_advance(sm->mom_calc, &phase_rng, &conf_rng, fin, sm->marr);
  }
}

// release memory for moment data object
void
vm_species_moment_release(const struct gkyl_vlasov_app *app, const struct vm_species_moment *sm)
{
  gkyl_array_release(sm->marr_host);
  gkyl_array_release(sm->marr);

  gkyl_array_release(sm->mom_hamil);

  if(sm->is_vlasov_lte_moms) {
    gkyl_vlasov_lte_moments_release(sm->vlasov_lte_moms);
  }
  else {
    gkyl_mom_type_release(sm->mom_type); 
    gkyl_mom_calc_release(sm->mom_calc); 
  }
}
