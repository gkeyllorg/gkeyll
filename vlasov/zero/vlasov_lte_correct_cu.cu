/* -*- c++ -*- */

extern "C" {
#include <gkyl_vlasov_lte_correct.h>
#include <gkyl_vlasov_lte_correct_priv.h>
#include <gkyl_range.h>
}

static void
gkyl_parallelize_components_kernel_launch_dims(dim3* dimGrid, dim3* dimBlock, gkyl_range range, int ncomp)
{
  // Create a 2D thread grid so we launch ncomp*range.volume number of threads and can parallelize over components too
  dimBlock->y = ncomp;
  dimGrid->y = 1;
  dimBlock->x = gkyl_int_div_up(240, ncomp);
  dimGrid->x = gkyl_int_div_up(range.volume, dimBlock->x);
}

__global__ static void
gkyl_vlasov_lte_correct_all_moments_abs_diff_cu_ker(struct gkyl_range conf_range, 
  int num_comp, int nc, 
  const struct gkyl_array *moms_target, const struct gkyl_array *moms_iter, 
  struct gkyl_array *abs_diff_moms)
{
  int idx[GKYL_MAX_DIM];

  // 2D thread grid
  // linc2 = c where c is the component index 
  long linc2 = threadIdx.y + blockIdx.y*blockDim.y;
  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);

    const double *moms_local = (const double*) gkyl_array_cfetch(moms_iter, loc);
    const double *moms_target_local = (const double*) gkyl_array_cfetch(moms_target, loc);
    double *abs_diff_moms_local = (double*) gkyl_array_fetch(abs_diff_moms, loc);

    // Compute the absolute value of the difference of cell averages 
    // Note: max error found by follow-up thread-safe reduction operation
    // Also: for density and temperature, this error is a relative error compared to the target moment value
    // so that we can converge to the correct target moments in SI units and minimize finite precision issues.
    // However, V_drift may be ~ 0 and if it is, we need to use absolute error. We can converge safely using
    // absolute error if V_drift ~ O(1). Otherwise, we use relative error for V_drift. 
    if (linc2 != 0 && linc2 != num_comp-1 && moms_target_local[linc2*nc] < 1.0) {
      abs_diff_moms_local[linc2] = fabs(moms_local[linc2*nc] - moms_target_local[linc2*nc]);
    }
    else {
      abs_diff_moms_local[linc2] = fabs(moms_local[linc2*nc] - moms_target_local[linc2*nc])/moms_target_local[linc2*nc];
    }
  }
}

void
gkyl_vlasov_lte_correct_all_moments_abs_diff_cu(const struct gkyl_range *conf_range, 
  int num_comp, int nc, 
  const struct gkyl_array *moms_target, const struct gkyl_array *moms_iter, 
  struct gkyl_array *moms_abs_diff)
{
  dim3 dimGrid, dimBlock;
  gkyl_parallelize_components_kernel_launch_dims(&dimGrid, &dimBlock, *conf_range, num_comp);
  gkyl_vlasov_lte_correct_all_moments_abs_diff_cu_ker<<<dimGrid, dimBlock>>>(*conf_range, 
    num_comp, nc, moms_target->on_dev, moms_iter->on_dev, moms_abs_diff->on_dev);
}