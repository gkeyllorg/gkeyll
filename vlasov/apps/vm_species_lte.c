#include <assert.h>
#include <gkyl_vlasov_priv.h>

void
vm_species_lte_init(
  struct gkyl_vlasov_app *app, struct vm_species *vms, struct vm_lte *lte,
  struct correct_all_moms_inp corr_inp
)
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
    .conf_range = &app->local,
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
  lte->proj_lte = gkyl_vlasov_lte_proj_on_basis_inew(&inp_proj);

  lte->correct_all_moms = corr_inp.correct_all_moms;
  // corr_inp arrives with defaults already applied by the caller (vm_species_init).
  int max_iter = corr_inp.max_iter;
  double iter_eps = corr_inp.iter_eps;
  bool use_last_converged = corr_inp.use_last_converged;

  if (lte->correct_all_moms) {
    struct gkyl_vlasov_lte_correct_inp inp_corr = {
      .phase_grid = &vms->grid,
      .vel_grid = &vms->grid_vel,
      .conf_basis = &app->basis,
      .vel_basis = &vms->basis_vel,
      .phase_basis = &vms->basis,
      .conf_range = &app->local,
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
    lte->corr_lte = gkyl_vlasov_lte_correct_inew(&inp_corr);
  }

  // Correction status of the LTE diagnostic; other users (e.g. BGK) add their own.
  lte->num_corr_channels = 0;
  lte->corr_channel_diag = vm_species_lte_add_corr_channel(lte, "lte");

  lte->f_lte = mkarr(app->use_gpu, vms->basis.num_basis, vms->local_ext.volume);
}

// Compute the LTE moments of fin, with the configuration-space Jacobian divided out of
// the density so they are physical
void
vm_species_lte_moms(
  gkyl_vlasov_app *app, const struct vm_species *vms, const struct vm_lte *lte,
  const struct gkyl_array *fin
)
{
  vm_species_moment_calc(&lte->moms, vms->local, app->local, fin);
  vm_species_moment_divide_Jx(app, vms, &lte->moms, lte->moms.marr);
}

// Compute f_lte from input LTE moments
int
vm_species_lte_add_corr_channel(struct vm_lte *lte, const char *name)
{
  assert(lte->num_corr_channels < VM_LTE_MAX_CORR_CHANNELS);
  int ch_idx = lte->num_corr_channels;
  lte->num_corr_channels += 1;

  struct vm_lte_corr_channel *ch = &lte->corr_channels[ch_idx];
  *ch = (struct vm_lte_corr_channel){.is_first_write_call = true};
  snprintf(ch->name, sizeof ch->name, "%s", name);
  ch->stats = gkyl_dynvec_new(GKYL_DOUBLE, VM_LTE_CORR_STAT_NCOMP);
  return ch_idx;
}

void
vm_species_lte_from_moms(
  gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_lte *lte,
  const struct gkyl_array *moms_lte, int corr_channel
)
{
  struct timespec wst = gkyl_wall_clock();

  gkyl_array_clear(lte->f_lte, 0.0);

  // Project the LTE distribution function to obtain f_lte.
  // e.g., Maxwellian for non-relativistic and Maxwell-Juttner for relativistic.
  // Projection routine also corrects the density of the projected distribution function.
  gkyl_vlasov_lte_proj_on_basis_advance(
    lte->proj_lte, &vms->local, &app->local, moms_lte, lte->f_lte
  );

  // Correct all the moments of the projected LTE distribution function.
  if (lte->correct_all_moms) {
    struct gkyl_vlasov_lte_correct_status status_corr;
    status_corr = gkyl_vlasov_lte_correct_all_moments(
      lte->corr_lte, lte->f_lte, moms_lte, &vms->local, &app->local
    );
    // Accumulate this correction into the caller's status channel (on this rank); the
    // reduction over ranks happens once per snapshot (vm_species_lte_calc_corr_status).
    struct vm_lte_corr_channel *ch = &lte->corr_channels[corr_channel];
    ch->num_corr += 1;
    ch->max_iter = fmax(ch->max_iter, status_corr.num_iter);
    ch->not_converged = ch->not_converged || status_corr.iter_converged; // 1 = not converged
    for (int i = 0; i < app->vdim + 2; ++i) {
      ch->max_error[i] = fmax(ch->max_error[i], status_corr.error[i]);
    }
    ch->num_dropped += status_corr.num_cells_dropped;

    lte->n_iter += status_corr.num_iter;
  }

  // The moments are physical, so weight f_lte by the configuration-space
  // Jacobian to match the stored J_x J_v f.
  if (!vms->pos_map->is_identity) {
    gkyl_vlasov_position_map_rescale_jacobpos(
      vms->pos_map, &vms->basis, &vms->local, lte->f_lte, lte->f_lte
    );
  }

  app->stat.species_lte_tm += gkyl_time_diff_now_sec(wst);
}

// Compute equivalent f_lte from fin
void
vm_species_lte(
  gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_lte *lte,
  const struct gkyl_array *fin, int corr_channel
)
{
  vm_species_lte_moms(app, vms, lte, fin);

  vm_species_lte_from_moms(app, vms, lte, lte->moms.marr, corr_channel);
}

void
vm_species_lte_calc_corr_status(gkyl_vlasov_app *app, struct vm_lte *lte, double tm)
{
  if (!lte->correct_all_moms) {
    return;
  }

  // Reduce every channel at once: the maxima (number of corrections, iterations, the
  // not-converged flag and the errors) and the sums of the dropped cells.
  int nch = lte->num_corr_channels;
  int nmax = VM_LTE_CORR_STAT_NCOMP - 1;
  double max_local[nch * nmax], max_global[nch * nmax];
  double sum_local[nch], sum_global[nch];
  for (int c = 0; c < nch; ++c) {
    const struct vm_lte_corr_channel *ch = &lte->corr_channels[c];
    double *mx = &max_local[c * nmax];
    mx[0] = ch->num_corr;
    mx[1] = ch->max_iter;
    mx[2] = ch->not_converged ? 1.0 : 0.0;
    for (int i = 0; i < 5; ++i) {
      mx[3 + i] = ch->max_error[i];
    }
    sum_local[c] = ch->num_dropped;
  }
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, nch * nmax, max_local, max_global);
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, nch, sum_local, sum_global);

  for (int c = 0; c < nch; ++c) {
    struct vm_lte_corr_channel *ch = &lte->corr_channels[c];
    const double *mx = &max_global[c * nmax];
    if (mx[0] > 0.0) {
      double row[VM_LTE_CORR_STAT_NCOMP];
      for (int i = 0; i < nmax; ++i) {
        row[i] = mx[i];
      }
      row[nmax] = sum_global[c];
      gkyl_dynvec_append(ch->stats, tm, row);
    }
    // Start the next interval.
    ch->num_corr = 0;
    ch->max_iter = 0.0;
    ch->not_converged = false;
    for (int i = 0; i < 5; ++i) {
      ch->max_error[i] = 0.0;
    }
    ch->num_dropped = 0.0;
  }
}

void
vm_species_lte_write_corr_status(gkyl_vlasov_app *app, struct vm_species *vms)
{
  if (!vms->lte.correct_all_moms) {
    return;
  }
  struct timespec wst = gkyl_wall_clock();

  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  for (int c = 0; c < vms->lte.num_corr_channels; ++c) {
    struct vm_lte_corr_channel *ch = &vms->lte.corr_channels[c];
    if (gkyl_dynvec_size(ch->stats) == 0) {
      continue;
    }
    if (rank == 0) {
      const char *fmt = "%s-%s-%s-corr-stat.gkyl";
      int sz = gkyl_calc_strlen(fmt, app->name, vms->name, ch->name);
      char fileNm[sz + 1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, vms->name, ch->name);
      if (ch->is_first_write_call) {
        // Write to a new file (this ensures previous output is removed).
        gkyl_dynvec_write(ch->stats, fileNm);
      } else {
        // Append to existing file.
        gkyl_dynvec_awrite(ch->stats, fileNm);
      }
    }
    ch->is_first_write_call = false;
    gkyl_dynvec_clear(ch->stats);
  }

  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_diag_io += 1;
}

void
vm_species_lte_release(const struct gkyl_vlasov_app *app, const struct vm_lte *lte)
{
  gkyl_array_release(lte->f_lte);

  vm_species_moment_release(app, &lte->moms);

  gkyl_vlasov_lte_proj_on_basis_release(lte->proj_lte);
  if (lte->correct_all_moms) {
    gkyl_vlasov_lte_correct_release(lte->corr_lte);
  }
  for (int c = 0; c < lte->num_corr_channels; ++c) {
    gkyl_dynvec_release(lte->corr_channels[c].stats);
  }
}
