#include <float.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_array_reduce.h>
#include <gkyl_vlasov_lte_correct.h>
#include <gkyl_vlasov_lte_correct_priv.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_sr_vars.h>
#include <gkyl_vlasov_lte_moments.h>
#include <gkyl_vlasov_lte_proj_on_basis.h>

#include <assert.h>

struct gkyl_vlasov_lte_correct *
gkyl_vlasov_lte_correct_inew(const struct gkyl_vlasov_lte_correct_inp *inp)
{
  gkyl_vlasov_lte_correct *up = gkyl_malloc(sizeof(*up));

  up->eps = inp->eps;
  up->max_iter = inp->max_iter;
  up->use_gpu = inp->use_gpu;
  up->use_last_converged = inp->use_last_converged;

  up->num_conf_basis = inp->conf_basis->num_basis;
  // (n, V_drift, T/m) being corrected
  // If the model is SR, V_drift is the spatial component of the four-velocity u_i = GammaV*V_drift
  int vdim = inp->phase_basis->ndim - inp->conf_basis->ndim;
  up->num_comp = vdim + 2;

  long conf_local_ncells = inp->conf_range->volume;
  long conf_local_ext_ncells = inp->conf_range_ext->volume;

  // Moments of the iterate, the accumulated correction (d) and the mismatch (dd) of the
  // target with the iterate; a per-cell mask of the cells still being corrected; and the
  // per-cell errors of the iterate and of the uncorrected projection, which the mask and
  // the thread-safe reduction of the maximum error act on.
  if (up->use_gpu) {
    up->moms_iter =
      gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_comp * up->num_conf_basis, conf_local_ext_ncells);
    up->d_moms =
      gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_comp * up->num_conf_basis, conf_local_ext_ncells);
    up->dd_moms =
      gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_comp * up->num_conf_basis, conf_local_ext_ncells);
    up->corr_mask = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, conf_local_ext_ncells);
    up->abs_diff_moms = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_comp, conf_local_ext_ncells);
    up->abs_diff_init = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_comp, conf_local_ext_ncells);
    up->error_cu = gkyl_cu_malloc(sizeof(double[up->num_comp]));
    up->mask_sum_cu = gkyl_cu_malloc(sizeof(double));
  } else {
    up->moms_iter =
      gkyl_array_new(GKYL_DOUBLE, up->num_comp * up->num_conf_basis, conf_local_ext_ncells);
    up->d_moms =
      gkyl_array_new(GKYL_DOUBLE, up->num_comp * up->num_conf_basis, conf_local_ext_ncells);
    up->dd_moms =
      gkyl_array_new(GKYL_DOUBLE, up->num_comp * up->num_conf_basis, conf_local_ext_ncells);
    up->corr_mask = gkyl_array_new(GKYL_DOUBLE, 1, conf_local_ext_ncells);
    up->abs_diff_moms = gkyl_array_new(GKYL_DOUBLE, up->num_comp, conf_local_ext_ncells);
    up->abs_diff_init = gkyl_array_new(GKYL_DOUBLE, up->num_comp, conf_local_ext_ncells);
  }
  // Allocate host-side error for checking convergence and returning in the status object
  up->error = gkyl_malloc(sizeof(double[up->num_comp]));

  // Moments structure
  struct gkyl_vlasov_lte_moments_inp inp_mom = {
    .phase_grid = inp->phase_grid,
    .vel_grid = inp->vel_grid,
    .conf_basis = inp->conf_basis,
    .vel_basis = inp->vel_basis,
    .phase_basis = inp->phase_basis,
    .conf_range = inp->conf_range,
    .conf_range_ext = inp->conf_range_ext,
    .vel_range = inp->vel_range,
    .phase_range = inp->phase_range,
    .vel_map = inp->vel_map,
    .hamil_range = inp->hamil_range,
    .hamil = inp->hamil,
    .model_id = inp->model_id,
    .hamil_id = inp->hamil_id,
    .gamma_inv = inp->gamma_inv,
    .h_ij = inp->h_ij,
    .h_ij_inv = inp->h_ij_inv,
    .det_h = inp->det_h,
    .use_extended_hamil_def = inp->use_extended_hamil_def,
    .effective_potential = inp->effective_potential,
    .use_gpu = inp->use_gpu,
  };
  up->moments_up = gkyl_vlasov_lte_moments_inew(&inp_mom);

  // Create a projection updater for projecting the LTE distribution function
  // Projection routine also corrects the density before returning
  // the LTE distribution function.
  struct gkyl_vlasov_lte_proj_on_basis_inp inp_proj = {
    .phase_grid = inp->phase_grid,
    .vel_grid = inp->vel_grid,
    .conf_basis = inp->conf_basis,
    .vel_basis = inp->vel_basis,
    .phase_basis = inp->phase_basis,
    .conf_range = inp->conf_range,
    .conf_range_ext = inp->conf_range_ext,
    .vel_range = inp->vel_range,
    .phase_range = inp->phase_range,
    .vel_map = inp->vel_map,
    .hamil_range = inp->hamil_range,
    .hamil = inp->hamil,
    .model_id = inp->model_id,
    .hamil_id = inp->hamil_id,
    .gamma_inv = inp->gamma_inv,
    .quad_type = inp->quad_type,
    .h_ij = inp->h_ij,
    .h_ij_inv = inp->h_ij_inv,
    .det_h = inp->det_h,
    .use_extended_hamil_def = inp->use_extended_hamil_def,
    .background_flows = inp->background_flows,
    .effective_potential = inp->effective_potential,
    .use_gpu = inp->use_gpu,
  };
  up->proj_lte = gkyl_vlasov_lte_proj_on_basis_inew(&inp_proj);

  return up;
}

// Per-cell errors of the moments in moms_iter against the target, dropping cells according
// to drop_mode, and the maximum error of each moment over the cells still being corrected.
static void
vlasov_lte_correct_errors(
  gkyl_vlasov_lte_correct *up, const struct gkyl_range *conf_local,
  const struct gkyl_array *moms_target, enum vlasov_lte_correct_drop_mode drop_mode
)
{
  if (up->use_gpu) {
    gkyl_vlasov_lte_correct_cell_errors_cu(
      conf_local, up->num_comp, up->num_conf_basis, up->eps, drop_mode, moms_target, up->moms_iter,
      up->abs_diff_init, up->corr_mask, up->abs_diff_moms
    );
    gkyl_array_reduce_range(up->error_cu, up->abs_diff_moms, GKYL_MAX, conf_local);
    gkyl_cu_memcpy(up->error, up->error_cu, sizeof(double[up->num_comp]), GKYL_CU_MEMCPY_D2H);
  } else {
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, conf_local);
    while (gkyl_range_iter_next(&iter)) {
      long loc = gkyl_range_idx(conf_local, iter.idx);
      vlasov_lte_correct_cell_errors(
        up->num_comp, up->num_conf_basis, up->eps, drop_mode, gkyl_array_cfetch(moms_target, loc),
        gkyl_array_cfetch(up->moms_iter, loc), gkyl_array_cfetch(up->abs_diff_init, loc),
        gkyl_array_fetch(up->corr_mask, loc), gkyl_array_fetch(up->abs_diff_moms, loc)
      );
    }
    gkyl_array_reduce_range(up->error, up->abs_diff_moms, GKYL_MAX, conf_local);
  }
}

// Number of cells dropped from the correction: the cells of the range whose mask is zero.
static int
vlasov_lte_correct_num_dropped(gkyl_vlasov_lte_correct *up, const struct gkyl_range *conf_local)
{
  double mask_sum = 0.0;
  if (up->use_gpu) {
    gkyl_array_reduce_range(up->mask_sum_cu, up->corr_mask, GKYL_SUM, conf_local);
    gkyl_cu_memcpy(&mask_sum, up->mask_sum_cu, sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    gkyl_array_reduce_range(&mask_sum, up->corr_mask, GKYL_SUM, conf_local);
  }
  return conf_local->volume - (long)round(mask_sum);
}

// Project the LTE distribution from the target moments plus the accumulated correction of
// the cells still being corrected (dropped cells have the correction removed).
static void
vlasov_lte_correct_project(
  gkyl_vlasov_lte_correct *up, struct gkyl_array *f_lte, const struct gkyl_array *moms_target,
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local
)
{
  gkyl_array_set(up->moms_iter, 1.0, moms_target);
  gkyl_array_accumulate(up->moms_iter, 1.0, up->d_moms);
  gkyl_vlasov_lte_proj_on_basis_advance(up->proj_lte, phase_local, conf_local, up->moms_iter, f_lte);
}

struct gkyl_vlasov_lte_correct_status
gkyl_vlasov_lte_correct_all_moments(
  gkyl_vlasov_lte_correct *up, struct gkyl_array *f_lte, const struct gkyl_array *moms_target,
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local
)
{
  int num_comp = up->num_comp;
  double tol = up->eps;
  int max_iter = up->max_iter;

  // Every cell starts in the correction with no accumulated correction.
  gkyl_array_clear(up->corr_mask, 1.0);
  gkyl_array_clear(up->d_moms, 0.0);

  // Fixed-point iteration on the moments handed to the projection: each pass adds the
  // mismatch of the target with the moments of the current iterate. Cells whose iterate
  // loses positivity are dropped: their correction is frozen and they leave the convergence
  // test, so one cell cannot hold the rest of the domain in the loop.
  int niter = 0;
  double max_error = DBL_MAX;
  while ((niter < max_iter) && (max_error > tol)) {
    gkyl_vlasov_lte_moments_advance(up->moments_up, phase_local, conf_local, f_lte, up->moms_iter);
    gkyl_array_set(up->dd_moms, -1.0, up->moms_iter);
    gkyl_array_accumulate(up->dd_moms, 1.0, moms_target);

    vlasov_lte_correct_errors(up, conf_local, moms_target, VLASOV_LTE_CORRECT_DROP_NONPOSITIVE);
    if (niter == 0) {
      // Error of the uncorrected projection, the fallback for cells that fail to converge.
      gkyl_array_copy(up->abs_diff_init, up->abs_diff_moms);
    }
    max_error = 0.0;
    for (int c = 0; c < num_comp; ++c) {
      max_error = fmax(max_error, up->error[c]);
    }

    gkyl_array_scale_by_cell(up->dd_moms, up->corr_mask);
    gkyl_array_accumulate(up->d_moms, 1.0, up->dd_moms);
    vlasov_lte_correct_project(up, f_lte, moms_target, phase_local, conf_local);

    niter += 1;
  }

  // Cells that did not converge within max_iter are dropped as well: all of them, or, when the
  // last iterate is to be kept, only those no closer to the target than the uncorrected
  // projection. Dropped cells are then reprojected from their target moments alone, which the
  // projection corrects for density, while the kept cells reproject to the same distribution.
  bool loop_converged = max_error <= tol;
  if (!loop_converged) {
    gkyl_vlasov_lte_moments_advance(up->moments_up, phase_local, conf_local, f_lte, up->moms_iter);
    vlasov_lte_correct_errors(
      up, conf_local, moms_target,
      up->use_last_converged ? VLASOV_LTE_CORRECT_DROP_NOT_IMPROVED :
                               VLASOV_LTE_CORRECT_DROP_ABOVE_TOL
    );
  }
  int num_dropped = vlasov_lte_correct_num_dropped(up, conf_local);
  if (num_dropped > 0) {
    gkyl_array_scale_by_cell(up->d_moms, up->corr_mask);
    vlasov_lte_correct_project(up, f_lte, moms_target, phase_local, conf_local);

    // Report the error of the returned distribution over every cell.
    gkyl_vlasov_lte_moments_advance(up->moments_up, phase_local, conf_local, f_lte, up->moms_iter);
    gkyl_array_clear(up->corr_mask, 1.0);
    vlasov_lte_correct_errors(up, conf_local, moms_target, VLASOV_LTE_CORRECT_DROP_NONE);
  }

  struct gkyl_vlasov_lte_correct_status status;
  status.iter_converged = !loop_converged || (num_dropped > 0);
  status.num_iter = niter;
  status.num_cells_dropped = num_dropped;
  for (int c = 0; c < num_comp; ++c) {
    status.error[c] = up->error[c];
  }
  return status;
}

void
gkyl_vlasov_lte_correct_release(gkyl_vlasov_lte_correct *up)
{
  gkyl_array_release(up->moms_iter);
  gkyl_array_release(up->d_moms);
  gkyl_array_release(up->dd_moms);
  gkyl_array_release(up->corr_mask);
  gkyl_array_release(up->abs_diff_moms);
  gkyl_array_release(up->abs_diff_init);
  if (up->use_gpu) {
    gkyl_cu_free(up->error_cu);
    gkyl_cu_free(up->mask_sum_cu);
  }
  gkyl_free(up->error);

  gkyl_vlasov_lte_moments_release(up->moments_up);
  gkyl_vlasov_lte_proj_on_basis_release(up->proj_lte);

  gkyl_free(up);
}

#ifndef GKYL_HAVE_CUDA

void
gkyl_vlasov_lte_correct_cell_errors_cu(
  const struct gkyl_range *conf_range, int num_comp, int nc, double tol, int drop_mode,
  const struct gkyl_array *moms_target, const struct gkyl_array *moms_iter,
  const struct gkyl_array *abs_diff_init, struct gkyl_array *corr_mask,
  struct gkyl_array *abs_diff_moms
)
{
  assert(false);
}

#endif
