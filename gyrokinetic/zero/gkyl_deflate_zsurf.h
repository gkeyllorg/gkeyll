#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_deflate_zsurf gkyl_deflate_zsurf;

/**
 * Create new updater to deflate a 2d (x,z)  or 3d (x,y,z) modal expansion to a 1d (x) modal expansion or 2d (x,y)
 *  at z=+/1 surfaces given a z index
 * 
 * @param cbasis Configuration-space basis.
 * @param deflated_cbasis Deflated configuration-space basis
 * @param edge Edge for surface (z=-1 or z=+1)
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */

struct gkyl_deflate_zsurf* 
gkyl_deflate_zsurf_new(const struct gkyl_basis *cbasis, const struct gkyl_basis *deflated_cbasis,
  int edge, bool use_gpu);

/**
 * Advance deflate_zsurf (compute the derived_zsurf coefficients).
 *
 * @param up deflate_zsurf updater object.
 * @param zidx Index for z surface
 * @param range Configuration-space range.
 * @param deflated_range Deflated configuration-space range.
 * @param field 2d (x,z) or 3d (x,y,z) field containing DG rep 
 * @param deflated_field 1d (x) or 2d (x,z) field containing DG rep
 * @param ncomp Number of components being deflated
 */

void gkyl_deflate_zsurf_advance(const struct gkyl_deflate_zsurf *up, int zidx, 
  const struct gkyl_range *range, const struct gkyl_range *deflated_range, 
  const struct gkyl_array *field, struct gkyl_array *deflated_field, int ncomp);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_deflate_zsurf_release(struct gkyl_deflate_zsurf* up);

/**
 * Host-side wrappers for deflation operations on device
 */

void gkyl_deflate_zsurf_advance_cu(const struct gkyl_deflate_zsurf *up, int zidx, 
  const struct gkyl_range *range, const struct gkyl_range *deflated_range, 
  const struct gkyl_array *field, struct gkyl_array *deflated_field, int ncomp);
