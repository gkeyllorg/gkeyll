#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type.
typedef struct gkyl_dg_fv_interpolate gkyl_dg_fv_interpolate;

/**
 * Create a new updater to interpolate a field from cell-average finite
 * volume representation to p=1 modal DG representation on a 2x coarsened grid
 * and vice-versa. 
 * 
 * @param basis 
 * @param grid_dg 
 * @param grid_fv 
 * @return struct gkyl_dg_fv_interpolate* 
 */
struct gkyl_dg_fv_interpolate*
gkyl_dg_fv_interpolate_new(struct gkyl_basis *basis, 
  struct gkyl_rect_grid *grid_dg, struct gkyl_rect_grid *grid_fv);

/**
  * Interpolate ncomp components from finite volume to p=1 modal DG.
  * Each DG cell requires data from the 2^ndim FV cells it covers.
  * The DG representation is constructed via recovery to match the respective
  * FV data when averaged over each FV cell.
  * 
  * @param up 
  * @param range_dg 
  * @param range_fv 
  * @param fv_in 
  * @param comp_in 
  * @param dg_out 
  * @param comp_out 
  * @param ncomp 
  */
void
gkyl_dg_fv_interpolate_fv_to_modal_advance(gkyl_dg_fv_interpolate *up, 
  struct gkyl_range *range_dg, struct gkyl_range *range_fv,
  const struct gkyl_array *fv_in, int comp_in,
  struct gkyl_array *dg_out, int comp_out, int ncomp);

/**
 * Interpolate ncomp components from p=1 modal DG to finite volume.
 * Each DG cell is broken into 2^ndim FV cells of equal size.
 * The FV data are computed by averaging the DG representation over
 * each FV cell.
 * 
 * @param up 
 * @param range_dg 
 * @param range_fv 
 * @param dg_in 
 * @param comp_in 
 * @param fv_out 
 * @param comp_out 
 * @param ncomp 
 */
void
gkyl_dg_fv_interpolate_modal_to_fv_advance(gkyl_dg_fv_interpolate *up, 
  struct gkyl_range *range_dg, struct gkyl_range *range_fv,
  const struct gkyl_array *dg_in, int comp_in,
  struct gkyl_array *fv_out, int comp_out, int ncomp);

/**
 * Release memory for interpolation object. 
 * 
 * @param up 
 */
void
gkyl_dg_fv_interpolate_release(gkyl_dg_fv_interpolate *up);