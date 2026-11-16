#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_deflate_geo gkyl_deflate_geo;

/**
 * Create an updater that restricts 3D geometry fields to a 1D or 2D grid
 * at the cell centers of the removed directions. Each removed direction
 * must contain one cell. Basis and grid objects must outlive the updater.
 *
 * @param cbasis Input 3D volume basis.
 * @param deflated_cbasis Output volume basis in the retained coordinates.
 * @param grid Input 3D configuration-space grid.
 * @param deflated_grid Output grid with the removed directions omitted.
 * @param rem_dirs Three-element array: 0 to retain a direction, 1 to remove it; copied.
 * @param use_gpu Reserved; currently unused.
 * @return New updater; release with gkyl_deflate_geo_release.
 */
gkyl_deflate_geo *gkyl_deflate_geo_new(
  const struct gkyl_basis *cbasis, const struct gkyl_basis *deflated_cbasis,
  const struct gkyl_rect_grid *grid, const struct gkyl_rect_grid *deflated_grid,
  const int *rem_dirs, bool use_gpu
);

/**
 * Create an updater that restricts geometry fields on 3D cell faces to
 * faces of a 1D or 2D grid, at the cell centers of the removed directions.
 * Each removed direction must contain one cell. Basis and grid objects
 * must outlive the updater.
 *
 * @param cbasis Input face basis.
 * @param deflated_num_basis Number of basis functions on an output face.
 * @param grid Input 3D configuration-space grid.
 * @param deflated_grid Output grid with the removed directions omitted.
 * @param rem_dirs Three-element array: 0 to retain a direction, 1 to remove it; copied.
 * @param dir Face-normal direction in the output grid.
 * @param use_gpu Reserved; currently unused.
 * @return New updater; release with gkyl_deflate_geo_surf_release.
 */
struct gkyl_deflate_geo_surf *gkyl_deflate_geo_surf_new(
  const struct gkyl_basis *cbasis, int deflated_num_basis, const struct gkyl_rect_grid *grid,
  const struct gkyl_rect_grid *deflated_grid, const int *rem_dirs, int dir, bool use_gpu
);

/**
 * Advance deflate_geo
 *
 * @param up deflate_geo updater object.
 * @param range Config-space range.
 * @param deflated_range range with one dimension removed
 * @param field 3d field
 * @param deflated_field 2d field on output
 * @param ncomp number of components
 */
void gkyl_deflate_geo_advance(
  const gkyl_deflate_geo *up, const struct gkyl_range *range,
  const struct gkyl_range *deflated_range, const struct gkyl_array *field,
  struct gkyl_array *deflated_field, int ncomp
);

/**
 * Restrict nodal quantities to the center of each removed coordinate.
 * Removed directions must contain one cell: two Gauss nodes are averaged,
 * while a three-node direction uses its center node.
 *
 * @param up deflate_geo updater object.
 * @param range Config-space range.
 * @param deflated_range range with one dimension removed
 * @param field 3d field
 * @param deflated_field 2d field on output
 * @param ncomp number of components
 */
void gkyl_deflate_geo_advance_nodal(
  const gkyl_deflate_geo *up, const struct gkyl_range *range,
  const struct gkyl_range *deflated_range, const struct gkyl_array *field,
  struct gkyl_array *deflated_field, int ncomp
);

/**
 * Advance deflate_geo_surf
 *
 * @param up deflate_geo updater object.
 * @param range Config-space range.
 * @param deflated_range range with one dimension removed
 * @param field 3d field
 * @param deflated_field 2d field on output
 * @param ncomp number of components
 */
void gkyl_deflate_geo_surf_advance(
  const struct gkyl_deflate_geo_surf *up, const struct gkyl_range *range,
  const struct gkyl_range *deflated_range, const struct gkyl_array *field,
  struct gkyl_array *deflated_field, int ncomp
);

/**
 * Restrict surface nodal quantities to the center of each removed coordinate.
 * Removed directions must contain one cell: two Gauss nodes are averaged,
 * while a three-node direction uses its center node.
 *
 * @param up deflate_geo updater object.
 * @param range Config-space nodal range.
 * @param deflated_range Nodal range with removed directions omitted.
 * @param field 3d field
 * @param deflated_field 2d field on output
 * @param ncomp number of components
 */
void gkyl_deflate_geo_surf_advance_nodal(
  const struct gkyl_deflate_geo_surf *up, const struct gkyl_range *range,
  const struct gkyl_range *deflated_range, const struct gkyl_array *field,
  struct gkyl_array *deflated_field, int ncomp
);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_deflate_geo_release(gkyl_deflate_geo *up);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_deflate_geo_surf_release(struct gkyl_deflate_geo_surf *up);
