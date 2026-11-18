/* -*- c++ -*- */

extern "C" {
#include <gkyl_vlasov_lte_correct.h>
#include <gkyl_vlasov_lte_correct_priv.h>
#include <gkyl_range.h>
}

// One thread per configuration-space cell: the per-cell error of every moment, with the
// cell dropped from the correction according to drop_mode.
__global__ static void
gkyl_vlasov_lte_correct_cell_errors_cu_ker(
  struct gkyl_range conf_range, int num_comp, int nc, double tol, int drop_mode,
  const struct gkyl_array *moms_target, const struct gkyl_array *moms_iter,
  const struct gkyl_array *abs_diff_init, struct gkyl_array *corr_mask,
  struct gkyl_array *abs_diff_moms
)
{
  int idx[GKYL_MAX_DIM];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < conf_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);

    vlasov_lte_correct_cell_errors(
      num_comp, nc, tol, drop_mode, (const double *)gkyl_array_cfetch(moms_target, loc),
      (const double *)gkyl_array_cfetch(moms_iter, loc),
      (const double *)gkyl_array_cfetch(abs_diff_init, loc),
      (double *)gkyl_array_fetch(corr_mask, loc), (double *)gkyl_array_fetch(abs_diff_moms, loc)
    );
  }
}

void
gkyl_vlasov_lte_correct_cell_errors_cu(
  const struct gkyl_range *conf_range, int num_comp, int nc, double tol, int drop_mode,
  const struct gkyl_array *moms_target, const struct gkyl_array *moms_iter,
  const struct gkyl_array *abs_diff_init, struct gkyl_array *corr_mask,
  struct gkyl_array *abs_diff_moms
)
{
  gkyl_vlasov_lte_correct_cell_errors_cu_ker<<<conf_range->nblocks, conf_range->nthreads>>>(
    *conf_range, num_comp, nc, tol, drop_mode, moms_target->on_dev, moms_iter->on_dev,
    abs_diff_init->on_dev, corr_mask->on_dev, abs_diff_moms->on_dev
  );
}
