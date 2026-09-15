#include <assert.h>
#include <gkyl_vlasov_priv.h>

void 
vm_species_lte_init(struct gkyl_vlasov_app *app, struct vm_species *vms, struct vm_lte *lte, 
  struct correct_all_moms_inp corr_inp)
{
  int cdim = app->cdim, vdim = app->vdim;

  // allocate moments needed for lte update
  vm_species_moment_init(app, vms, &lte->moms, GKYL_F_MOMENT_LTE, false);

  struct gkyl_vlasov_lte_proj_on_basis_inp inp_proj = {
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
    .hamil_range = &vms->mom_hamil_range,
    .hamil = vms->mom_hamil,
    .model_id = vms->model_id,
    .hamil_id = vms->mom_hamil_id,
    .gamma_inv = vms->gamma_inv,
    .h_ij = vms->h_ij,
    .h_ij_inv = vms->h_ij_inv,
    .det_h = vms->det_h,
    .use_extended_hamil_def = vms->info.use_extended_hamil_def,
    .background_flows = vms->background_flows,
    .effective_potential = vms->effective_potential,
    .use_gpu = app->use_gpu,
  };
  lte->proj_lte = gkyl_vlasov_lte_proj_on_basis_inew( &inp_proj );

  lte->correct_all_moms = corr_inp.correct_all_moms;
  int max_iter = corr_inp.max_iter > 0 ? vms->info.correct.max_iter : 100;
  double iter_eps = corr_inp.iter_eps > 0 ? vms->info.correct.iter_eps  : 1e-12;
  bool use_last_converged = corr_inp.use_last_converged;
  
  if (lte->correct_all_moms) {
    struct gkyl_vlasov_lte_correct_inp inp_corr = {
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
      .hamil_range = &vms->mom_hamil_range,
      .hamil = vms->mom_hamil,
      .model_id = vms->model_id,
      .hamil_id = vms->mom_hamil_id,
      .gamma_inv = vms->gamma_inv,
      .h_ij = vms->h_ij,
      .h_ij_inv = vms->h_ij_inv,
      .det_h = vms->det_h,
      .use_extended_hamil_def = vms->info.use_extended_hamil_def,
      .background_flows = vms->background_flows,
      .effective_potential = vms->effective_potential,
      .use_gpu = app->use_gpu,
      .max_iter = max_iter,
      .eps = iter_eps,
      .use_last_converged = use_last_converged, 
    };
    lte->n_iter = 0;
    lte->corr_lte = gkyl_vlasov_lte_correct_inew( &inp_corr );

    lte->corr_stat = gkyl_dynvec_new(GKYL_DOUBLE, 7);
    lte->is_first_corr_status_write_call = true;
  }

  lte->f_lte = mkarr(app->use_gpu, vms->basis.num_basis, vms->local_ext.volume);
}

// Compute f_lte from input LTE moments
void
vm_species_lte_from_moms(gkyl_vlasov_app *app, const struct vm_species *vms,
  struct vm_lte *lte, const struct gkyl_array *moms_lte)
{
  struct timespec wst = gkyl_wall_clock();

  gkyl_array_clear(lte->f_lte, 0.0);

  // Project the LTE distribution function to obtain f_lte.
  // e.g., Maxwellian for non-relativistic and Maxwell-Juttner for relativistic.
  // Projection routine also corrects the density of the projected distribution function.
  gkyl_vlasov_lte_proj_on_basis_advance(lte->proj_lte, &vms->local, &app->local, 
    moms_lte, lte->f_lte);

  // Correct all the moments of the projected LTE distribution function.
  if (lte->correct_all_moms) {
    struct gkyl_vlasov_lte_correct_status status_corr;
    status_corr = gkyl_vlasov_lte_correct_all_moments(lte->corr_lte, lte->f_lte, moms_lte,
      &vms->local, &app->local);
    double corr_vec[7] = { 0.0 };
    corr_vec[0] = status_corr.num_iter;
    corr_vec[1] = status_corr.iter_converged;
    for (int i=0; i<app->vdim+2; ++i) {
      corr_vec[2+i] = status_corr.error[i];
    }
    double corr_vec_global[7] = { 0.0 };
    gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, 7, corr_vec, corr_vec_global);    
    gkyl_dynvec_append(lte->corr_stat, app->tcurr, corr_vec_global);

    lte->n_iter += status_corr.num_iter;
  } 

  app->stat.species_lte_tm += gkyl_time_diff_now_sec(wst);   
}

// Compute equivalent f_lte from fin
void
vm_species_lte(gkyl_vlasov_app *app, const struct vm_species *vms,
  struct vm_lte *lte, const struct gkyl_array *fin)
{
  vm_species_moment_calc(&lte->moms, vms->local, app->local, fin);

  vm_species_lte_from_moms(app, vms, lte, lte->moms.marr);
}

void
vm_species_lte_write_max_corr_status(gkyl_vlasov_app* app, struct vm_species *vms)
{
  if (vms->lte.correct_all_moms) {
    struct timespec wst = gkyl_wall_clock();

    int rank;
    gkyl_comm_get_rank(app->comm, &rank);
    if (rank == 0) {
      // Write out correction status.
      const char *fmt = "%s-%s-%s.gkyl";
      int sz = gkyl_calc_strlen(fmt, app->name, vms->info.name, "corr-max-stat");
      char fileNm[sz+1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, vms->info.name, "corr-max-stat");

      if (vms->lte.is_first_corr_status_write_call) {
        // Write to a new file (this ensure previous output is removed).
        gkyl_dynvec_write(vms->lte.corr_stat, fileNm);
        vms->lte.is_first_corr_status_write_call = false;
      }
      else {
        // Append to existing file.
        gkyl_dynvec_awrite(vms->lte.corr_stat, fileNm);
      }
    }
    gkyl_dynvec_clear(vms->lte.corr_stat);

    app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
    app->stat.n_diag_io += 1;    
  }
}

void 
vm_species_lte_release(const struct gkyl_vlasov_app *app, const struct vm_lte *lte)
{
  gkyl_array_release(lte->f_lte);

  vm_species_moment_release(app, &lte->moms);

  gkyl_vlasov_lte_proj_on_basis_release(lte->proj_lte);
  if (lte->correct_all_moms) {
    gkyl_vlasov_lte_correct_release(lte->corr_lte);
    gkyl_dynvec_release(lte->corr_stat);
  }
}
