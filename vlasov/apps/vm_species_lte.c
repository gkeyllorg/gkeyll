#include <assert.h>
#include <gkyl_vlasov_priv.h>

void 
vm_species_lte_init(struct gkyl_vlasov_app *app, struct vm_species *s, struct vm_lte *lte, 
  struct correct_all_moms_inp corr_inp)
{
  int cdim = app->cdim, vdim = app->vdim;

  // allocate moments needed for lte update
  vm_species_moment_init(app, s, &lte->moms, GKYL_F_MOMENT_LTE, false);

  struct gkyl_vlasov_lte_proj_on_basis_inp inp_proj = {
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
  lte->proj_lte = gkyl_vlasov_lte_proj_on_basis_inew( &inp_proj );

  lte->correct_all_moms = corr_inp.correct_all_moms;
  int max_iter = corr_inp.max_iter > 0 ? s->info.max_iter : 100;
  double iter_eps = corr_inp.iter_eps > 0 ? s->info.iter_eps  : 1e-12;
  bool use_last_converged = corr_inp.use_last_converged;
  
  if (lte->correct_all_moms) {
    struct gkyl_vlasov_lte_correct_inp inp_corr = {
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
      .max_iter = max_iter,
      .eps = iter_eps,
      .use_last_converged = use_last_converged, 
    };
    lte->niter = 0;
    lte->corr_lte = gkyl_vlasov_lte_correct_inew( &inp_corr );

    lte->corr_stat = gkyl_dynvec_new(GKYL_DOUBLE, 7);
    lte->is_first_corr_status_write_call = true;
  }

  lte->f_lte = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
}

// Compute f_lte from input LTE moments
void
vm_species_lte_from_moms(gkyl_vlasov_app *app, const struct vm_species *species,
  struct vm_lte *lte, const struct gkyl_array *moms_lte)
{
  struct timespec wst = gkyl_wall_clock();

  gkyl_array_clear(lte->f_lte, 0.0);

  // Project the LTE distribution function to obtain f_lte.
  // e.g., Maxwellian for non-relativistic and Maxwell-Juttner for relativistic.
  // Projection routine also corrects the density of the projected distribution function.
  gkyl_vlasov_lte_proj_on_basis_advance(lte->proj_lte, &species->local, &app->local, 
    moms_lte, lte->f_lte);

  // Correct all the moments of the projected LTE distribution function.
  if (lte->correct_all_moms) {
    struct gkyl_vlasov_lte_correct_status status_corr;
    status_corr = gkyl_vlasov_lte_correct_all_moments(lte->corr_lte, lte->f_lte, moms_lte,
      &species->local, &app->local);
    double corr_vec[7] = { 0.0 };
    corr_vec[0] = status_corr.num_iter;
    corr_vec[1] = status_corr.iter_converged;
    for (int i=0; i<app->vdim+2; ++i) {
      corr_vec[2+i] = status_corr.error[i];
    }
    double corr_vec_global[7] = { 0.0 };
    gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, 7, corr_vec, corr_vec_global);    
    gkyl_dynvec_append(lte->corr_stat, app->tcurr, corr_vec_global);

    lte->niter += status_corr.num_iter;
  } 

  app->stat.species_lte_tm += gkyl_time_diff_now_sec(wst);   
}

// Compute equivalent f_lte from fin
void
vm_species_lte(gkyl_vlasov_app *app, const struct vm_species *species,
  struct vm_lte *lte, const struct gkyl_array *fin)
{
  vm_species_moment_calc(&lte->moms, species->local, app->local, fin);

  vm_species_lte_from_moms(app, species, lte, lte->moms.marr);
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
