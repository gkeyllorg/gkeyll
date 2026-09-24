#pragma once

#include <gkyl_array.h>
#include <gkyl_evalf_def.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

// Geometry information for a single cell: recall a cell "owns" the
// faces on the lower side of cell
struct gkyl_wave_cell_geom {
  double kappa; // ratio of cell-volume in phy to comp space
  double lenr[GKYL_MAX_CDIM]; // ratio of face-area in phys to comp space for "lower" faces
  double norm[GKYL_MAX_CDIM][GKYL_MAX_CDIM]; // norm[d] is the normal to face perp to direction 'd'
  // tau1[d] X tau2[d] = norm[d] are tangents to face perp to direction 'd'
  double tau1[GKYL_MAX_CDIM][GKYL_MAX_CDIM];
  double tau2[GKYL_MAX_CDIM][GKYL_MAX_CDIM];
};

// geometry information over a range of cells
struct gkyl_wave_geom {
  struct gkyl_range range; // range over which geometry is defined
  struct gkyl_array *geom; // geometry in each cell

  uint32_t flags;
  struct gkyl_ref_count ref_count;  
  struct gkyl_wave_geom *on_dev; // pointer to itself or device object
};

/**
 * Create a new wave geometry object. 
 *
 * @param grid Grid on which geometry lives
 * @param range Range on which geometry should be constructed
 * @param mapc2p Mapping from computational to physical space
 * @param ctx Context for use in mapping
 */
struct gkyl_wave_geom*
gkyl_wave_geom_new(const struct gkyl_rect_grid *grid,
  struct gkyl_range *range, evalf_t mapc2p, void *ctx, bool use_gpu);

/**
 * Create a new wave geometry object that lives on NV-GPU: see new() method
 * above for documentation.
 */
struct gkyl_wave_geom*
gkyl_wave_geom_cu_dev_new(const struct gkyl_rect_grid *grid,
  struct gkyl_range *range, evalf_t mapc2p, void *ctx);

/**
 * Acquire pointer to geometry object. The pointer must be released
 * using gkyl_wave_geom_release method.
 *
 * @param wg Geometry to which a pointer is needed
 * @return Pointer to acquired geometry
 */
struct gkyl_wave_geom* gkyl_wave_geom_acquire(const struct gkyl_wave_geom* wg);

/**
 * Get pointer to geometry in a cell given by idx into the range over
 * which the geometry was constructed.
 *
 * @param wg Wave geometry object
 * @param idx Index into grid
 * @return cell geometry in cell @a idx
 */
GKYL_CU_DH
static inline const struct gkyl_wave_cell_geom*
gkyl_wave_geom_get(const struct gkyl_wave_geom *wg, const int *idx)
{
  return (const struct gkyl_wave_cell_geom*) gkyl_array_cfetch(wg->geom, gkyl_range_idx(&wg->range, idx));
}

/**
 * Release geometry object.
 *
 * @param wg Wave geometry object to release.
 */
void gkyl_wave_geom_release(const struct gkyl_wave_geom *wg);
