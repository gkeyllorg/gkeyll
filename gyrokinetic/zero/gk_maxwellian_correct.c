#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_array_reduce.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_gk_maxwellian_correct.h>
#include <gkyl_gk_maxwellian_correct_priv.h>
#include <gkyl_gk_maxwellian_moments.h>
#include <gkyl_gk_maxwellian_proj_on_basis.h>

#include <assert.h>

struct gkyl_gk_maxwellian_correct*
gkyl_gk_maxwellian_correct_inew(const struct gkyl_gk_maxwellian_correct_inp *inp)
{
  gkyl_gk_maxwellian_correct *up = gkyl_malloc(sizeof(*up));
  up->eps = inp->eps;
  up->max_iter = inp->max_iter;
  up->use_last_converged = inp->use_last_converged;
  up->use_gpu = inp->use_gpu;

  up->num_conf_basis = inp->conf_basis->num_basis;
  
  long conf_range_ncells = inp->conf_range->volume;
  long conf_range_ext_ncells = inp->conf_range_ext->volume;

  up->num_comp = 3; // correct n, u_par, T/m
  up->bimaxwellian = inp->bimaxwellian;
  if (up->bimaxwellian) {
    up->num_comp = 4; // correct n, u_par, T_par/m, T_perp/m
  }

  // Individual moment memory: the iteration of the moments, the differences (d) and differences of differences (dd)
  // In the gyrokinetic system, we are correcting three moments for Maxwellian corrections: n, u_par, T/m 
  // or if we are correcting a BiMaxwellian, we correct four moments: n, u_par, T_par/m, T_perp/m
  if (up->use_gpu) {
    up->moms_iter = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_comp*inp->conf_basis->num_basis, conf_range_ext_ncells);
    up->d_moms = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_comp*inp->conf_basis->num_basis, conf_range_ext_ncells);
    up->dd_moms = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_comp*inp->conf_basis->num_basis, conf_range_ext_ncells);
    // Two additional GPU-specific allocations for iterating over the grid to find the absolute value of 
    // the difference between the target and iterative moments, and the GPU-side array for performing the
    // thread-safe reduction to find the maximum error on the grid.
    up->abs_diff_moms = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_comp, conf_range_ext_ncells);
    up->error_cu = gkyl_cu_malloc(sizeof(double[up->num_comp]));
  }
  else {
    up->moms_iter = gkyl_array_new(GKYL_DOUBLE, up->num_comp*inp->conf_basis->num_basis, conf_range_ext_ncells);
    up->d_moms = gkyl_array_new(GKYL_DOUBLE, up->num_comp*inp->conf_basis->num_basis, conf_range_ext_ncells);
    up->dd_moms = gkyl_array_new(GKYL_DOUBLE, up->num_comp*inp->conf_basis->num_basis, conf_range_ext_ncells);
  }
  up->error = gkyl_malloc(sizeof(double[up->num_comp]));

  // Moments structure 
  struct gkyl_gk_maxwellian_moments_inp inp_mom = {
    .phase_grid = inp->phase_grid,
    .conf_basis = inp->conf_basis,
    .phase_basis = inp->phase_basis,
    .conf_range =  inp->conf_range,
    .conf_range_ext = inp->conf_range_ext,
    .mass = inp->mass, 
    .gk_geom = inp->gk_geom, 
    .vel_map = inp->vel_map,
    .divide_jacobgeo = inp->divide_jacobgeo, 
    .use_gpu = inp->use_gpu,
  };
  up->moments_up = gkyl_gk_maxwellian_moments_inew( &inp_mom );

  // Create a projection updater for projecting the gyrokinetic Maxwellian or bi-Maxwellian
  struct gkyl_gk_maxwellian_proj_on_basis_inp inp_proj = {
    .phase_grid = inp->phase_grid,
    .conf_basis = inp->conf_basis,
    .phase_basis = inp->phase_basis,
    .conf_range = inp->conf_range,
    .conf_range_ext = inp->conf_range_ext,
    .vel_range = inp->vel_range,
    .gk_geom = inp->gk_geom, 
    .vel_map = inp->vel_map,
    .mass = inp->mass, 
    .bimaxwellian = inp->bimaxwellian, 
    .divide_jacobgeo = inp->divide_jacobgeo, 
    .use_gpu = inp->use_gpu,
  };
  up->proj_max = gkyl_gk_maxwellian_proj_on_basis_inew(&inp_proj);

  return up;
}

struct gkyl_gk_maxwellian_correct_status
gkyl_gk_maxwellian_correct_all_moments(gkyl_gk_maxwellian_correct *up,
  struct gkyl_array *f_max, const struct gkyl_array *moms_target, 
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range)
{
  int num_comp = up->num_comp;
  int nc = up->num_conf_basis;
  double tol = up->eps;  // tolerance of the iterative scheme
  int max_iter = up->max_iter;

  int niter = 0;
  bool corr_status = true;
  int ispositive_f_lte = true;

  // Set initial max error to start the iteration.
  double max_error = 1.0;
  for (int i=0; i<num_comp; ++i) {
    up->error[i] = 1.0;
  }
  // Copy the initial max error to GPU so initial error is set correctly (no uninitialized values).
  if (up->use_gpu) {
    gkyl_cu_memcpy(up->error_cu, up->error, sizeof(double[num_comp]), GKYL_CU_MEMCPY_H2D);
  }

  // Clear the differences prior to iteration
  gkyl_array_clear(up->d_moms, 0.0);
  gkyl_array_clear(up->dd_moms, 0.0);

  // Iteration loop, max_iter iterations is usually sufficient for machine precision moments
  while ((ispositive_f_lte) && ((niter < max_iter) && (max_error > tol))) {
    // 1. Calculate the needed moments from the projected distribution
    // either Maxwellian (n, u_par, T/m) or bi-Maxwellian (n, u_par, Tpar/m, Tperp/m) moments
    if (up->bimaxwellian) {
      gkyl_gk_bimaxwellian_moments_advance(up->moments_up, phase_range, conf_range, 
        f_max, up->moms_iter);
    }
    else {
      gkyl_gk_maxwellian_moments_advance(up->moments_up, phase_range, conf_range, 
        f_max, up->moms_iter);
    }

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
        gkyl_gk_maxwellian_correct_all_moments_abs_diff_cu(conf_range, 
          num_comp, nc, moms_target, up->moms_iter, up->abs_diff_moms);
        gkyl_array_reduce_range(up->error_cu, up->abs_diff_moms, GKYL_MAX, conf_range);
        gkyl_cu_memcpy(up->error, up->error_cu, sizeof(double[num_comp]), GKYL_CU_MEMCPY_D2H);
      }
      else {
        struct gkyl_range_iter biter;

        // Reset the maximum error
        for (int i=0; i<num_comp; ++i) {
          up->error[i] = 0.0;
        }
        // Iterate over the input configuration-space range to find the maximum error
        gkyl_range_iter_init(&biter, conf_range);
        while (gkyl_range_iter_next(&biter)) {
          long midx = gkyl_range_idx(conf_range, biter.idx);
          const double *moms_local = gkyl_array_cfetch(up->moms_iter, midx);
          const double *moms_target_local = gkyl_array_cfetch(moms_target, midx);
          // Check the error in the absolute value of the cell average
          // Note: for density and temperature(s), this error is a relative error compared to the target moment value
          // so that we can converge to the correct target moments in SI units and minimize finite precision issues.
          up->error[0] = fmax(fabs(moms_local[0*nc] - moms_target_local[0*nc])/moms_target_local[0*nc],fabs(up->error[0]));
          up->error[2] = fmax(fabs(moms_local[2*nc] - moms_target_local[2*nc])/moms_target_local[2*nc],fabs(up->error[2]));
          // However, u_par may be ~ 0 and if it is, we normalize it with the target thermal veocity instead.
          if (fabs(moms_target_local[1*nc]) < sqrt(moms_target_local[2*nc])) {
            up->error[1] = fmax(fabs(moms_local[1*nc] - moms_target_local[1*nc])/sqrt(moms_target_local[2*nc]),fabs(up->error[1]));
          }
          else {
            up->error[1] = fmax(fabs(moms_local[1*nc] - moms_target_local[1*nc])/fabs(moms_target_local[1*nc]),fabs(up->error[1]));
          }
          // Check if density and temperature (or parallel temperature) are positive, 
          // if they aren't we will break out of the iteration
          ispositive_f_lte = (moms_local[0*nc]>0.0) && ispositive_f_lte;
          ispositive_f_lte = (moms_local[2*nc]>0.0) && ispositive_f_lte;

          // Also compute the error in Tperp/m if we are correcting a bi-Maxwellian
          if (up->bimaxwellian) {
            up->error[3] = fmax(fabs(moms_local[3*nc] - moms_target_local[3*nc])/moms_target_local[3*nc],fabs(up->error[3]));  
            ispositive_f_lte = (moms_local[3*nc]>0.0) && ispositive_f_lte;          
          }
        }
      }
    }
    // Find the maximum error looping over the error in each component
    max_error = 0.0; // reset maximum error 
    for (int d=0; d<num_comp; ++d) {
      max_error = fmax(max_error, up->error[d]);
    }

    // c. Calculate  n^(k+1) = M^k + dM^(k+1)
    // n = n_target + dm_new;
    gkyl_array_set(up->moms_iter, 1.0, moms_target);
    gkyl_array_accumulate(up->moms_iter, 1.0, up->d_moms);

    // 2. Update the gyrokinetic Maxwellian distribution function using the corrected moments.
    gkyl_gk_maxwellian_proj_on_basis_advance(up->proj_max,
      phase_range, conf_range, up->moms_iter, false, f_max);

    niter += 1;
  }

  if ((niter < max_iter) && (ispositive_f_lte) && (max_error < tol)) {
    corr_status = 0;
  } 
  else {
    corr_status = 1;
  }

  // If the algorithm fails to converge and we are *not* using the results of the failed convergence,
  // we project the distribution function with the target moments.
  // We correct the density and then recompute moments/errors for this new projection.
  if (corr_status == 1 && !up->use_last_converged) {
    gkyl_gk_maxwellian_proj_on_basis_advance(up->proj_max,
      phase_range, conf_range, moms_target, false, f_max);

    if (up->bimaxwellian) {
      gkyl_gk_bimaxwellian_moments_advance(up->moments_up, phase_range, conf_range, 
        f_max, up->moms_iter);
    }
    else {
      gkyl_gk_maxwellian_moments_advance(up->moments_up, phase_range, conf_range, 
        f_max, up->moms_iter);
    }

    if (up->use_gpu) {
      // We insure the reduction to find the maximum error is thread-safe on GPUs
      // by first calling a specialized kernel for computing the absolute value 
      // of the difference of the cell averages, then calling reduce_range.
      gkyl_gk_maxwellian_correct_all_moments_abs_diff_cu(conf_range, 
        num_comp, nc, moms_target, up->moms_iter, up->abs_diff_moms);
      gkyl_array_reduce_range(up->error_cu, up->abs_diff_moms, GKYL_MAX, conf_range);
      gkyl_cu_memcpy(up->error, up->error_cu, sizeof(double[num_comp]), GKYL_CU_MEMCPY_D2H);
    }
    else {
      struct gkyl_range_iter biter;

      // Reset the maximum error
      for (int i=0; i<num_comp; ++i) {
        up->error[i] = 0.0;
      }
      // Iterate over the input configuration-space range to find the maximum error
      gkyl_range_iter_init(&biter, conf_range);
      while (gkyl_range_iter_next(&biter)) {
        long midx = gkyl_range_idx(conf_range, biter.idx);
        const double *moms_local = gkyl_array_cfetch(up->moms_iter, midx);
        const double *moms_target_local = gkyl_array_cfetch(moms_target, midx);
        // Check the error in the absolute value of the cell average
        // Note: for density and temperature, this error is a relative error compared to the target moment value
        // so that we can converge to the correct target moments in SI units and minimize finite precision issues.
        up->error[0] = fmax(fabs(moms_local[0*nc] - moms_target_local[0*nc])/moms_target_local[0*nc],fabs(up->error[0]));
        up->error[2] = fmax(fabs(moms_local[2*nc] - moms_target_local[2*nc])/moms_target_local[2*nc],fabs(up->error[2]));
        // However, u_par may be ~ 0 and if it is, we need to use absolute error. We can converge safely using
        // absolute error if u_par ~ O(1). Otherwise, we use relative error for u_par. 
        if (fabs(moms_target_local[1*nc]) < 1.0) {
          up->error[1] = fmax(fabs(moms_local[1*nc] - moms_target_local[1*nc]),fabs(up->error[1]));
        }
        else {
          up->error[1] = fmax(fabs(moms_local[1*nc] - moms_target_local[1*nc])/moms_target_local[1*nc],fabs(up->error[1]));
        }

        // Also compute the error in Tperp/m if we are correcting a bi-Maxwellian
        if (up->bimaxwellian) {
          up->error[3] = fmax(fabs(moms_local[3*nc] - moms_target_local[3*nc])/moms_target_local[3*nc],fabs(up->error[3]));  
          ispositive_f_lte = (moms_local[3*nc]>0.0) && ispositive_f_lte;          
        }
      }
    }
  }

  struct gkyl_gk_maxwellian_correct_status status;
  status.iter_converged = corr_status;
  status.num_iter = niter;
  for (int i=0; i<num_comp; ++i) {
    status.error[i] = up->error[i];
  }
  return status;
}

void 
gkyl_gk_maxwellian_correct_release(gkyl_gk_maxwellian_correct *up)
{
  gkyl_array_release(up->moms_iter);
  gkyl_array_release(up->d_moms);
  gkyl_array_release(up->dd_moms);
  if (up->use_gpu) {
    gkyl_array_release(up->abs_diff_moms);
    gkyl_cu_free(up->error_cu);
  }
  gkyl_free(up->error);

  gkyl_gk_maxwellian_moments_release(up->moments_up);
  gkyl_gk_maxwellian_proj_on_basis_release(up->proj_max);

  gkyl_free(up);
}

#ifndef GKYL_HAVE_CUDA

void 
gkyl_gk_maxwellian_correct_all_moments_abs_diff_cu(const struct gkyl_range *conf_range, 
  int num_comp, int nc, 
  const struct gkyl_array *moms_target, const struct gkyl_array *moms_iter, 
  struct gkyl_array *moms_abs_diff)
{
  assert(false);
}

#endif
