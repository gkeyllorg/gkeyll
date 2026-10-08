#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_array_reduce.h>
#include <gkyl_dist_type.h>
#include <gkyl_vlasov_dist_correct.h>
#include <gkyl_vlasov_dist_correct_priv.h>
#include <gkyl_vlasov_dist_moments.h>
#include <gkyl_vlasov_dist_proj_on_basis.h>

#include <assert.h>

struct gkyl_vlasov_dist_correct *
gkyl_vlasov_dist_correct_inew(const struct gkyl_vlasov_dist_correct_inp *inp)
{
  gkyl_vlasov_dist_correct *up = gkyl_malloc(sizeof(*up));
  
  up->dist_proj = gkyl_dist_proj_type_acquire(inp->dist_proj);
  up->vel_map = 0;
  if (inp->vel_map != 0) {
    up->vel_map = gkyl_velocity_map_acquire(inp->vel_map);
  }

  up->eps = inp->eps;
  up->max_iter = inp->max_iter;
  up->use_gpu = inp->use_gpu;
  up->use_last_converged = inp->use_last_converged;

  up->num_conf_basis = inp->conf_basis->num_basis;
  up->vdim = inp->phase_basis->ndim - inp->conf_basis->ndim;
  up->num_mom = gkyl_dist_proj_type_num_mom(up->dist_proj);

  long conf_local_ncells = inp->conf_range->volume;
  long conf_local_ext_ncells = inp->conf_range_ext->volume;

  // Individual moment memory: the iteration of the moments, the differences (d) and differences of differences (dd)
  if (up->use_gpu) {
    up->moms_iter =
      gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_mom * up->num_conf_basis, conf_local_ext_ncells);
    up->d_moms =
      gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_mom * up->num_conf_basis, conf_local_ext_ncells);
    up->dd_moms =
      gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_mom * up->num_conf_basis, conf_local_ext_ncells);
    // Two additional GPU-specific allocations for iterating over the grid to find the absolute value of
    // the difference between the target and iterative moments, and the GPU-side array for performing the
    // thread-safe reduction to find the maximum error on the grid.
    up->abs_diff_moms = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_mom, conf_local_ext_ncells);
    up->error_cu = gkyl_cu_malloc(sizeof(double[up->num_mom]));
  } else {
    up->moms_iter =
      gkyl_array_new(GKYL_DOUBLE, up->num_mom * up->num_conf_basis, conf_local_ext_ncells);
    up->d_moms =
      gkyl_array_new(GKYL_DOUBLE, up->num_mom * up->num_conf_basis, conf_local_ext_ncells);
    up->dd_moms =
      gkyl_array_new(GKYL_DOUBLE, up->num_mom * up->num_conf_basis, conf_local_ext_ncells);
  }
  // Allocate host-side error for checking convergence and returning in the status object
  up->error = gkyl_malloc(sizeof(double[up->num_mom]));

  // Moments structure
  struct gkyl_vlasov_dist_moments_inp inp_mom = {
    .phase_grid = inp->phase_grid,
    .vel_grid = inp->vel_grid,
    .conf_basis = inp->conf_basis,
    .vel_basis = inp->vel_basis,
    .phase_basis = inp->phase_basis,
    .conf_range = inp->conf_range,
    .conf_range_ext = inp->conf_range_ext,
    .vel_range = inp->vel_range,
    .phase_range = inp->phase_range,
    .dist_id = up->dist_proj->dist_id,
    .model_id = inp->model_id,
    .mass = inp->mass,
    .use_gpu = inp->use_gpu,
  };
  up->moments_up = gkyl_vlasov_dist_moments_inew(&inp_mom);

  // Generalized distribution projection updater.
  struct gkyl_vlasov_dist_proj_on_basis_inp inp_proj = {
    .dist_proj = up->dist_proj,
    .phase_grid = inp->phase_grid,
    .vel_grid = inp->vel_grid,
    .conf_basis = inp->conf_basis,
    .vel_basis = inp->vel_basis,
    .phase_basis = inp->phase_basis,
    .conf_range = inp->conf_range,
    .conf_range_ext = inp->conf_range_ext,
    .vel_range = inp->vel_range,
    .phase_range = inp->phase_range,
    .model_id = inp->model_id,
    .use_gpu = inp->use_gpu,
  };
  up->proj_dist = gkyl_vlasov_dist_proj_on_basis_inew(&inp_proj);

  return up;
}

struct gkyl_vlasov_dist_correct_status
gkyl_vlasov_dist_correct_all_moments(
  gkyl_vlasov_dist_correct *up, struct gkyl_array *f_dist, const struct gkyl_array *moms_target,
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local
)
{
  int num_mom = up->num_mom;
  int vdim = up->vdim;
  int nc = up->num_conf_basis;
  double tol = up->eps; // tolerance of the iterative scheme
  int max_iter = up->max_iter;

  int niter = 0;
  bool corr_status = true;
  int ispositive_f_dist = true;

  // Set initial max error to start the iteration.
  double max_error = 1.0;
  for (int i = 0; i < num_mom; ++i) {
    up->error[i] = 1.0;
  }
  // Copy the initial max error to GPU so initial error is set correctly (no uninitialized values).
  if (up->use_gpu) {
    gkyl_cu_memcpy(up->error_cu, up->error, sizeof(double[num_mom]), GKYL_CU_MEMCPY_H2D);
  }

  // Clear the differences prior to iteration
  gkyl_array_clear(up->d_moms, 0.0);
  gkyl_array_clear(up->dd_moms, 0.0);

  // Iteration loop, max_iter iterations is usually sufficient for machine precision moments
  while ((ispositive_f_dist) && ((niter < max_iter) && (max_error > tol))) {
    // 1. Calculate the moments from the projected distribution.
    gkyl_vlasov_dist_moments_advance(up->moments_up, phase_local, conf_local, f_dist, up->moms_iter);

    // a. Calculate  ddMi^(k+1) =  Mi_corr - Mi_new
    // ddn = n_target - n;
    // Compute out = out + a*inp. Returns out.
    gkyl_array_set(up->dd_moms, -1.0, up->moms_iter);
    gkyl_array_accumulate(up->dd_moms, 1.0, moms_target);

    // b. Calculate  dMi^(k+1) = dn^k + ddMi^(k+1) | where dn^0 = 0
    // dm_new = dm_old + ddn;
    gkyl_array_accumulate(up->d_moms, 1.0, up->dd_moms);

    // End the iteration early if all moments converge
    if ((niter % 1) == 0) {
      if (up->use_gpu) {
        // We insure the reduction to find the maximum error is thread-safe on GPUs
        // by first calling a specialized kernel for computing the absolute value
        // of the difference of the cell averages, then calling reduce_range.
        gkyl_vlasov_dist_correct_all_moments_abs_diff_cu(
          conf_local, num_mom, nc, vdim,  moms_target, up->moms_iter, up->abs_diff_moms
        );
        gkyl_array_reduce_range(up->error_cu, up->abs_diff_moms, GKYL_MAX, conf_local);
        gkyl_cu_memcpy(up->error, up->error_cu, sizeof(double[num_mom]), GKYL_CU_MEMCPY_D2H);
      } else {
        struct gkyl_range_iter biter;

        // Reset the maximum error
        for (int i = 0; i < num_mom; ++i) {
          up->error[i] = 0.0;
        }
        // Iterate over the input configuration-space range to find the maximum error
        gkyl_range_iter_init(&biter, conf_local);
        while (gkyl_range_iter_next(&biter)) {
          long midx = gkyl_range_idx(conf_local, biter.idx);
          const double *moms_local = gkyl_array_cfetch(up->moms_iter, midx);
          const double *moms_target_local = gkyl_array_cfetch(moms_target, midx);
          // Check the error in the absolute value of the cell average
          // Note: for density and temperature, this error is a relative error compared to the target moment value
          // so that we can converge to the correct target moments in SI units and minimize finite precision issues.
          up->error[0] = fmax(
            fabs(moms_local[0 * nc] - moms_target_local[0 * nc]) / moms_target_local[0 * nc],
            fabs(up->error[0])
          );

          // However, V_drift may be ~ 0 and if it is, we need to use absolute error. We can converge safely using
          // absolute error if V_drift ~ O(1). Otherwise, we use relative error for V_drift.
          for (int d = 1; d < vdim+1; ++d) {
            if (fabs(moms_target_local[d * nc]) < 1.0) {
              up->error[d] =
                fmax(fabs(moms_local[d * nc] - moms_target_local[d * nc]), fabs(up->error[d]));
            } else {
              up->error[d] = fmax(
                fabs(moms_local[d * nc] - moms_target_local[d * nc]) / fabs(moms_target_local[d * nc]),
                fabs(up->error[d])
              );
            }
          }
          
          // Higher moments use relative error.
          for (int d=vdim+1; d<num_mom; ++d) {
            up->error[d] = fmax(
              fabs(moms_local[d*nc] - moms_target_local[d*nc]) /
                fabs(moms_target_local[d*nc]),
              fabs(up->error[d])
            );
          }
          // Check if density is positive, if not we will break out of the iteration
          ispositive_f_dist = (moms_local[0 * nc] > 0.0) && ispositive_f_dist;
        }
      }
    }
    // Find the maximum error looping over the error in each component
    max_error = 0.0; // reset maximum error
    for (int d = 0; d < num_mom; ++d) {
      max_error = fmax(max_error, up->error[d]);
    }

    // c. Calculate  M_i^(k+1) = M_i^k + dM_i^(k+1)
    // n = n_target + dm_new;
    gkyl_array_set(up->moms_iter, 1.0, moms_target);
    gkyl_array_accumulate(up->moms_iter, 1.0, up->d_moms);

    // 2. Update the distribution function using the corrected moments.
    // Projection routine also corrects the density before the next iteration.
    gkyl_vlasov_dist_proj_on_basis_moments_advance(
    up->proj_dist, conf_local, up->moms_iter
    );

    gkyl_vlasov_dist_proj_on_basis_advance(
    up->proj_dist, phase_local, conf_local, f_dist
    );


    niter += 1;
  }

  if ((niter < max_iter) && (ispositive_f_dist) && (max_error < tol)) {
    corr_status = 0;
  } else {
    corr_status = 1;
  }

  // If the algorithm fails to converge and we are *not* using the results of the failed convergence,
  // we project the distribution function with the target moments.
  // We correct the density and then recompute moments/errors for this new projection.
  if (corr_status == 1 && !up->use_last_converged) {
    gkyl_vlasov_dist_proj_on_basis_moments_advance(up->proj_dist, conf_local, moms_target);

    gkyl_vlasov_dist_proj_on_basis_advance(up->proj_dist, phase_local, conf_local, f_dist);
    if (up->use_gpu) {
      // We insure the reduction to find the maximum error is thread-safe on GPUs
      // by first calling a specialized kernel for computing the absolute value
      // of the difference of the cell averages, then calling reduce_range.
      gkyl_vlasov_dist_correct_all_moments_abs_diff_cu(
        conf_local, num_mom, nc, vdim, moms_target, up->moms_iter, up->abs_diff_moms
      );
      gkyl_array_reduce_range(up->error_cu, up->abs_diff_moms, GKYL_MAX, conf_local);
      gkyl_cu_memcpy(up->error, up->error_cu, sizeof(double[num_mom]), GKYL_CU_MEMCPY_D2H);
    } else {
      struct gkyl_range_iter biter;

      // Reset the maximum error
      for (int i = 0; i < num_mom; ++i) {
        up->error[i] = 0.0;
      }
      // Iterate over the input configuration-space range to find the maximum error
      gkyl_range_iter_init(&biter, conf_local);
      while (gkyl_range_iter_next(&biter)) {
        long midx = gkyl_range_idx(conf_local, biter.idx);
        const double *moms_local = gkyl_array_cfetch(up->moms_iter, midx);
        const double *moms_target_local = gkyl_array_cfetch(moms_target, midx);
        // Check the error in the absolute value of the cell average
        // Note: for density and temperature, this error is a relative error compared to the target moment value.
        up->error[0] = fmax(
          fabs(moms_local[0 * nc] - moms_target_local[0 * nc]) / moms_target_local[0 * nc],
          fabs(up->error[0])
        );

        // However, V_drift may be ~ 0 and if it is, we need to use absolute error.
        // Otherwise, we use relative error for V_drift.
        for (int d = 1; d < vdim + 1; ++d) {
          if (fabs(moms_target_local[d * nc]) < 1.0) {
            up->error[d] =
              fmax(fabs(moms_local[d * nc] - moms_target_local[d * nc]), fabs(up->error[d]));
          } else {
            up->error[d] = fmax(
              fabs(moms_local[d * nc] - moms_target_local[d * nc]) / fabs(moms_target_local[d * nc]),
              fabs(up->error[d])
            );
          }      
        }

        for (int d=vdim+1; d<num_mom; ++d) {
          up->error[d] = fmax(
            fabs(moms_local[d*nc] - moms_target_local[d*nc]) /
              fabs(moms_target_local[d*nc]),
            fabs(up->error[d])
          );
        }
      }
    }
  }

  struct gkyl_vlasov_dist_correct_status status;
  status.iter_converged = corr_status;
  status.num_iter = niter;
  status.error = up->error;

  return status;
}

void
gkyl_vlasov_dist_correct_release(gkyl_vlasov_dist_correct *up)
{
  if (up->vel_map != 0) {
    gkyl_velocity_map_release(up->vel_map);
  }

  gkyl_array_release(up->moms_iter);
  gkyl_array_release(up->d_moms);
  gkyl_array_release(up->dd_moms);
  if (up->use_gpu) {
    gkyl_array_release(up->abs_diff_moms);
    gkyl_cu_free(up->error_cu);
  }
  gkyl_free(up->error);

  gkyl_vlasov_dist_moments_release(up->moments_up);
  gkyl_vlasov_dist_proj_on_basis_release(up->proj_dist);
  gkyl_dist_proj_type_release(up->dist_proj);
  gkyl_free(up);
}

#ifndef GKYL_HAVE_CUDA

void
gkyl_vlasov_dist_correct_all_moments_abs_diff_cu(
  const struct gkyl_range *conf_range, int num_mom, int nc, int vdim, const struct gkyl_array *moms_target,
  const struct gkyl_array *moms_iter, struct gkyl_array *moms_abs_diff
)
{
  assert(false);
}

#endif