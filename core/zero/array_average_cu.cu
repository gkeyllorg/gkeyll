/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_average.h>
#include <gkyl_array_average_priv.h>
#include <gkyl_array_ops.h>
}

__global__ static void
gkyl_array_average_set_ker_cu(struct gkyl_array_average *up)
{
  array_average_choose_kernel(up);
}

struct gkyl_array_average *
gkyl_array_average_cu_dev_new(struct gkyl_array_average *up)
{
  struct gkyl_array *weight_ho;
  if (up->isweighted) {
    weight_ho = gkyl_array_acquire(up->weight);
  } else {
    weight_ho = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->weight->ncomp, up->weight->size);
    gkyl_array_copy(weight_ho, up->weight);
  }
  gkyl_array_release(up->weight);
  up->weight = weight_ho->on_dev;

  // Copy struct to device and set the kernel there.
  struct gkyl_array_average *up_cu =
    (struct gkyl_array_average *)gkyl_cu_malloc(sizeof(struct gkyl_array_average));
  gkyl_cu_memcpy(up_cu, up, sizeof(struct gkyl_array_average), GKYL_CU_MEMCPY_H2D);
  gkyl_array_average_set_ker_cu<<<1, 1>>>(up_cu);

  up->weight = weight_ho;
  up->on_dev = up_cu;

  return up;
}

__global__ static void
gkyl_array_average_advance_cu_ker(
  const struct gkyl_array_average *up, const struct gkyl_array *fin, struct gkyl_array *avgout
)
{
  // One thread per input cell, accumulating into its output cell.
  int idx[GKYL_MAX_DIM];
  for (unsigned long tid = threadIdx.x + blockIdx.x * blockDim.x; tid < up->local.volume;
       tid += blockDim.x * gridDim.x) {
    gkyl_sub_range_inv_idx(&up->local, tid, idx);
    long lidx = gkyl_range_idx(&up->local, idx);
    long lidx_out = array_average_out_lidx(up, idx);

    const double *fin_c = (const double *)gkyl_array_cfetch(fin, lidx);
    const double *win_c = (const double *)gkyl_array_cfetch(up->weight, up->isweighted ? lidx : 0);
    double *out_c = (double *)gkyl_array_fetch(avgout, lidx_out);

    up->kernel(up->subvol, win_c, fin_c, out_c);
  }
}

__global__ static void
gkyl_array_average_broadcast_cu_ker(const struct gkyl_array_average *up, struct gkyl_array *avgout)
{
  // Copy the average into every other cell along the averaged dimensions.
  int idx[GKYL_MAX_DIM];
  for (unsigned long tid = threadIdx.x + blockIdx.x * blockDim.x; tid < up->local_avg.volume;
       tid += blockDim.x * gridDim.x) {
    gkyl_sub_range_inv_idx(&up->local_avg, tid, idx);
    long lidx = gkyl_range_idx(&up->local_avg, idx);
    long lidx_rep = array_average_out_lidx(up, idx);
    if (lidx != lidx_rep) {
      const double *rep_c = (const double *)gkyl_array_cfetch(avgout, lidx_rep);
      double *out_c = (double *)gkyl_array_fetch(avgout, lidx);
      for (int k = 0; k < avgout->ncomp; k++) {
        out_c[k] = rep_c[k];
      }
    }
  }
}

void
gkyl_array_average_advance_cu(
  const struct gkyl_array_average *up, const struct gkyl_array *fin, struct gkyl_array *avgout
)
{
  gkyl_array_clear_range(avgout, 0.0, &up->local_avg);

  gkyl_array_average_advance_cu_ker<<<up->local.nblocks, up->local.nthreads>>>(
    up->on_dev, fin->on_dev, avgout->on_dev
  );

  if (up->full_dim_out) {
    gkyl_array_average_broadcast_cu_ker<<<up->local_avg.nblocks, up->local_avg.nthreads>>>(
      up->on_dev, avgout->on_dev
    );
  }

  if (up->isweighted) {
    gkyl_dg_div_op_range(
      up->div_mem, &up->basis_avg, 0, avgout, 0, avgout, 0, up->weight_avg, &up->local_avg
    );
  }
}
