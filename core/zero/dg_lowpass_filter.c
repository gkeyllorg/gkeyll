#include <gkyl_dg_lowpass_filter.h>
#include <gkyl_dg_lowpass_filter_priv.h>
#include <gkyl_alloc.h>

struct gkyl_dg_lowpass_filter*
gkyl_dg_lowpass_filter_new(int dir, int half_width, double cutoff_wavelength,
  const struct gkyl_basis *basis, const struct gkyl_rect_grid *grid,
  const struct gkyl_range *range, bool use_gpu)
{
  // Allocate space for new updater.
  struct gkyl_dg_lowpass_filter *up = gkyl_malloc(sizeof(*up));

  up->use_gpu = use_gpu;
  up->ndim = basis->ndim;
  up->dir = dir;
  up->half_width = half_width;
  up->num_basis = basis->num_basis;
  up->grid = *grid;
  up->range = *range;

  // Reflect the stencil only at faces that coincide with the domain boundary.
  up->reflect_lo = range->lower[dir] == 1;
  up->reflect_up = range->upper[dir] == grid->cells[dir];

  // Perform some basic checks:
  assert(grid->ndim == range->ndim);
  assert(0 <= dir && dir < grid->ndim);
  assert(half_width > 0);

  // Normalized cutoff frequency in cycles per cell.
  double fc = grid->dx[dir]/cutoff_wavelength;
  assert(0.0 < fc && fc <= 0.5);

  up->weights = gkyl_malloc((2*half_width+1)*sizeof(double));
  dg_lpf_calc_weights(half_width, fc, up->weights);

  // Sign each coefficient picks up when its donor cell is reached by mirroring.
  up->sign_plain = gkyl_malloc(up->num_basis*sizeof(double));
  up->sign_mirror = gkyl_malloc(up->num_basis*sizeof(double));
  for (int c=0; c<up->num_basis; c++) up->sign_plain[c] = 1.0;
  basis->flip_odd_sign(dir, up->sign_plain, up->sign_mirror);

  up->weights_cu = up->sign_plain_cu = up->sign_mirror_cu = NULL;
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    up->weights_cu = gkyl_cu_malloc((2*half_width+1)*sizeof(double));
    up->sign_plain_cu = gkyl_cu_malloc(up->num_basis*sizeof(double));
    up->sign_mirror_cu = gkyl_cu_malloc(up->num_basis*sizeof(double));
    gkyl_cu_memcpy(up->weights_cu, up->weights, (2*half_width+1)*sizeof(double), GKYL_CU_MEMCPY_H2D);
    gkyl_cu_memcpy(up->sign_plain_cu, up->sign_plain, up->num_basis*sizeof(double), GKYL_CU_MEMCPY_H2D);
    gkyl_cu_memcpy(up->sign_mirror_cu, up->sign_mirror, up->num_basis*sizeof(double), GKYL_CU_MEMCPY_H2D);
  }
#endif

  return up;
}

void
gkyl_dg_lowpass_filter_advance(gkyl_dg_lowpass_filter *up,
  struct gkyl_array *GKYL_RESTRICT fdo, struct gkyl_array *GKYL_RESTRICT ftar)
{
  assert(fdo != ftar); // Stencil operation, can't be done in-place.
  assert(fdo->ncomp == ftar->ncomp);
  assert(fdo->size == ftar->size);
  assert(fdo->ncomp == up->num_basis);

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_dg_lowpass_filter_advance_cu(up, fdo, ftar);
    return;
  }
#endif

  int dir = up->dir;
  int M = up->half_width;
  int num_basis = up->num_basis;
  int idx_do[GKYL_MAX_DIM];

  // Loop over the target range.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->range);
  while (gkyl_range_iter_next(&iter)) {

    long linidx_tar = gkyl_range_idx(&up->range, iter.idx);
    double *ftar_c = gkyl_array_fetch(ftar, linidx_tar);

    gkyl_copy_int_arr(up->range.ndim, iter.idx, idx_do);

    for (int c=0; c<num_basis; c++)
      ftar_c[c] = 0.0;

    // Loop over the donor cells contributing to this target cell.
    double wsum = 0.0;
    for (int k=-M; k<M+1; k++) {
      bool mirrored = false;
      int idx_k = iter.idx[dir]+k;
      if (idx_k < up->range.lower[dir]) {
        if (!up->reflect_lo) continue; // Interior edge: drop donor, renormalize.
        idx_do[dir] = dg_lpf_mirror_idx(idx_k, up->range.lower[dir], up->range.upper[dir], &mirrored);
      }
      else if (idx_k > up->range.upper[dir]) {
        if (!up->reflect_up) continue; // Interior edge: drop donor, renormalize.
        idx_do[dir] = dg_lpf_mirror_idx(idx_k, up->range.lower[dir], up->range.upper[dir], &mirrored);
      }
      else {
        idx_do[dir] = idx_k;
      }

      long linidx_do = gkyl_range_idx(&up->range, idx_do);
      const double *fdo_c = gkyl_array_cfetch(fdo, linidx_do);

      // Handle reflection.
      const double *sgn = mirrored? up->sign_mirror : up->sign_plain;

      double w = up->weights[k+M];
      wsum += w;
      for (int c=0; c<num_basis; c++)
        ftar_c[c] += w*sgn[c]*fdo_c[c];
    }

    // Renormalize.
    for (int c=0; c<num_basis; c++)
      ftar_c[c] /= wsum;
  }
}

void
gkyl_dg_lowpass_filter_release(gkyl_dg_lowpass_filter *up)
{
  // Release memory associated with this updater.
  gkyl_free(up->weights);
  gkyl_free(up->sign_plain);
  gkyl_free(up->sign_mirror);
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_cu_free(up->weights_cu);
    gkyl_cu_free(up->sign_plain_cu);
    gkyl_cu_free(up->sign_mirror_cu);
  }
#endif
  gkyl_free(up);
}
