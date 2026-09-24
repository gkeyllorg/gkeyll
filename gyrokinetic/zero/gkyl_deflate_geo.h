#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_deflate_geo gkyl_deflate_geo;

/**
 * Create new updater to compute deflated go 
 *
 * @param cbasis Basis object (configuration space).
 * @param deflated_cbasis Basis object (configuration space) with one less dimension
 * @param grid configuration space grid
 * @param defated_grid configuration space grid with onne less dimension
 * @param rem_dirs array of directions to remove. 0 to keep, 1 to remove
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */


gkyl_deflate_geo* gkyl_deflate_geo_new(const struct gkyl_basis *cbasis, const struct gkyl_basis *deflated_cbasis,
  const struct gkyl_rect_grid *grid, const struct gkyl_rect_grid *deflated_grid, const int *rem_dirs, bool use_gpu);
/**
 * Create new updater to compute the deflated surface geo
 *
 * @param cbasis Basis object (configuration space).
 * @param deflated_num_basis number of basis elements in deflated basis
 * @param grid configuration space grid
 * @param defated_grid configuration space grid with onne less dimension
 * @param rem_dirs array of directions to remove. 0 to keep, 1 to remove
 * @param dir direction of surface
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */

struct gkyl_deflate_geo_surf* gkyl_deflate_geo_surf_new(const struct gkyl_basis *cbasis, int deflated_num_basis, const struct gkyl_rect_grid *grid, const struct gkyl_rect_grid *deflated_grid, const int *rem_dirs, int dir, bool use_gpu);

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


void gkyl_deflate_geo_advance(const gkyl_deflate_geo *up, const struct gkyl_range *range, const struct gkyl_range* deflated_range, const struct gkyl_array *field, struct gkyl_array *deflated_field, int ncomp);


/**
 * Advance deflate_geo on nodal quantities
 *
 * @param up deflate_geo updater object.
 * @param range Config-space range.
 * @param deflated_range range with one dimension removed
 * @param field 3d field
 * @param deflated_field 2d field on output
 * @param ncomp number of components
 */


void gkyl_deflate_geo_advance_nodal(const gkyl_deflate_geo *up, const struct gkyl_range *range, const struct gkyl_range* deflated_range, const struct gkyl_array *field, struct gkyl_array *deflated_field, int ncomp);


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
void gkyl_deflate_geo_surf_advance(const struct gkyl_deflate_geo_surf *up, const struct gkyl_range *range, const struct gkyl_range* deflated_range, const struct gkyl_array *field, struct gkyl_array *deflated_field, int ncomp);

/**
 * Advance deflate_geo_surf on nodal quantities
 *
 * @param up deflate_geo updater object.
 * @param nrange Config-space range.
 * @param deflated_nrange range with one dimension removed
 * @param field 3d field
 * @param deflated_field 2d field on output
 * @param ncomp number of components
 */
void gkyl_deflate_geo_surf_advance_nodal(const struct gkyl_deflate_geo_surf *up, const struct gkyl_range *range, const struct gkyl_range* deflated_range, const struct gkyl_array *field, struct gkyl_array *deflated_field, int ncomp);


/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_deflate_geo_release(gkyl_deflate_geo* up);


/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_deflate_geo_surf_release(struct gkyl_deflate_geo_surf* up);
