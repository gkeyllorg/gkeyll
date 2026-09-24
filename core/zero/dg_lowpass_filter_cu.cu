/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_dg_lowpass_filter.h>
#include <gkyl_dg_lowpass_filter_priv.h>
#include <gkyl_util.h>
}

__global__ static void
gkyl_dg_lowpass_filter_advance_cu_ker(int dir, int M, int num_basis,
  bool reflect_lo, bool reflect_up,
  const double *GKYL_RESTRICT weights, const double *GKYL_RESTRICT sign_plain,
  const double *GKYL_RESTRICT sign_mirror, struct gkyl_range range,
  const struct gkyl_array *GKYL_RESTRICT fdo, struct gkyl_array *GKYL_RESTRICT ftar)
{
  int idx_tar[GKYL_MAX_DIM];
  int idx_do[GKYL_MAX_DIM];

  for (unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < range.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&range, tid, idx_tar);

    long linidx_tar = gkyl_range_idx(&range, idx_tar);
    double *ftar_c = (double *) gkyl_array_fetch(ftar, linidx_tar);

    for (int c=0; c<num_basis; c++)
      ftar_c[c] = 0.0;

    for (int d=0; d<range.ndim; d++)
      idx_do[d] = idx_tar[d];

    // Loop over the donor cells contributing to this target cell.
    double wsum = 0.0;
    for (int k=-M; k<M+1; k++) {
      bool mirrored = false;
      int idx_k = idx_tar[dir]+k;
      if (idx_k < range.lower[dir]) {
        if (!reflect_lo) continue; // Interior edge: drop donor, renormalize.
        idx_do[dir] = dg_lpf_mirror_idx(idx_k, range.lower[dir], range.upper[dir], &mirrored);
      }
      else if (idx_k > range.upper[dir]) {
        if (!reflect_up) continue; // Interior edge: drop donor, renormalize.
        idx_do[dir] = dg_lpf_mirror_idx(idx_k, range.lower[dir], range.upper[dir], &mirrored);
      }
      else {
        idx_do[dir] = idx_k;
      }

      long linidx_do = gkyl_range_idx(&range, idx_do);
      const double *fdo_c = (const double *) gkyl_array_cfetch(fdo, linidx_do);

      const double *sgn = mirrored? sign_mirror : sign_plain;

      double w = weights[k+M];
      wsum += w;
      for (int c=0; c<num_basis; c++)
        ftar_c[c] += w*sgn[c]*fdo_c[c];
    }

    // Renormalize to preserve the zero mode when donors were dropped.
    for (int c=0; c<num_basis; c++)
      ftar_c[c] /= wsum;
  }
}

void
gkyl_dg_lowpass_filter_advance_cu(gkyl_dg_lowpass_filter *up,
  struct gkyl_array *GKYL_RESTRICT fdo, struct gkyl_array *GKYL_RESTRICT ftar)
{
  int nblocks = up->range.nblocks, nthreads = up->range.nthreads;

  gkyl_dg_lowpass_filter_advance_cu_ker<<<nblocks, nthreads>>>
    (up->dir, up->half_width, up->num_basis, up->reflect_lo, up->reflect_up,
     up->weights_cu, up->sign_plain_cu, up->sign_mirror_cu, up->range,
     fdo->on_dev, ftar->on_dev);
}
