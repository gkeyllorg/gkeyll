#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type.
typedef struct gkyl_dg_lowpass_filter gkyl_dg_lowpass_filter;

/**
 * Create a new updater that low-pass filters a DG field along one direction
 * with a Blackman-windowed sinc kernel, normalized to preserve the zero mode.
 *
 * @param dir Direction along which to filter.
 * @param half_width Stencil half-width M in cells (stencil spans 2M+1 cells).
 * @param cutoff_wavelength Cutoff wavelength (physical units). Normalized
 *   cutoff f_c = dx/cutoff_wavelength must satisfy 0 < f_c <= 0.5.
 * @param basis DG basis of the filtered field.
 * @param grid Grid the filtered field is defined on.
 * @param range Range to filter in. The stencil is reflected or truncated at the edges.
 * @param use_gpu bool to determine if on GPU.
 * @return New filter updater.
 */
struct gkyl_dg_lowpass_filter*
gkyl_dg_lowpass_filter_new(int dir, int half_width, double cutoff_wavelength,
  const struct gkyl_basis *basis, const struct gkyl_rect_grid *grid,
  const struct gkyl_range *range, bool use_gpu);

/**
 * Run the filter updater. Cannot be used in-place; ftar cells outside the
 * range are left untouched.
 *
 * @param up Filter updater.
 * @param fdo Donor field.
 * @param ftar Target (filtered) field.
 */
void
gkyl_dg_lowpass_filter_advance(gkyl_dg_lowpass_filter *up,
  struct gkyl_array *fdo, struct gkyl_array *ftar);

/**
 * Release the memory associated with this filter updater.
 *
 * @param up Filter updater.
 */
void
gkyl_dg_lowpass_filter_release(gkyl_dg_lowpass_filter *up);
