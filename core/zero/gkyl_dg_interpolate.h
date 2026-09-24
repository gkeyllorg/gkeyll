#pragma once

#include <gkyl_basis.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_ops.h>
#include <gkyl_rect_decomp.h>
#include <string.h>
#include <gkyl_range.h>
#include <gkyl_array.h>

// Object type.
typedef struct gkyl_dg_interpolate gkyl_dg_interpolate;

/**
 * Create a new updater to interpolate a field from a donor grid to a target
 * grid with the same extents but different resolution.
 *
 * @param cdim Configuration-space dimensions.
 * @param basis DG basis.
 * @param grid_do Donor grid.
 * @param grid_tar Target grid.
 * @param range_do Range in donor grid.
 * @param range_tar Range in target grid.
 * @param nghost Number of ghost cells in each direction.
 * @param use_gpu bool to determine if on GPU.
 * @return New interpolation updater.
 */
struct gkyl_dg_interpolate*
gkyl_dg_interpolate_new(int cdim, const struct gkyl_basis *basis,
  const struct gkyl_rect_grid *grid_do, const struct gkyl_rect_grid *grid_tar,
  const struct gkyl_range *range_do, const struct gkyl_range *range_tar,
  const int *nghost, bool use_gpu);

/**
 * Run the interpolation updater in the indicated range.
 *
 * @param up Interpolation updater.
 * @param fdo Donor field.
 * @param ftar Target field.
 */
void
gkyl_dg_interpolate_advance(gkyl_dg_interpolate* up,
  struct gkyl_array *fdo, struct gkyl_array *ftar);

/**
 * Release the memory associated with this interpolating updater.
 *
 * @param up Interpolation updater.
 */
void
gkyl_dg_interpolate_release(gkyl_dg_interpolate* up);
