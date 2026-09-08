#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_rz_calc_derived_geo gkyl_rz_calc_derived_geo;

/**
 * Create new updater to compute the derived_geo coefficients
 *
 * @param cbasis Basis object (configuration space).
 * @param grid computational grid
 * @param node_type 0 for serendipity, 1 for gauss-legendre quadrature nodes
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_rz_calc_derived_geo *gkyl_rz_calc_derived_geo_new(
  const struct gkyl_basis *cbasis, const struct gkyl_rect_grid *grid, int node_type, bool use_gpu);

/**
 * Advance rz_calc_derived_geo (compute the derived_geo coefficients).
 *
 * @param up rz_calc_derived_geo updater object.
 * @param crange Config-space range.
 * @param g_ij input field containing DG rep of g_ij
 * @param bmag input field containing DG rep of B = |B|
 * @param jacobgeo input field with DG rep of jacobian (J)
 * @param jacobgeo_inv output field with DG rep of 1/J
 * @param gij output field with DG rep of  g^ij
 * @param jacobtot output field with DG rep of JB
 * @param gxxj output field with DG rep of Jg^xx
 * @param gxyj output field with DG rep of Jg^xy
 * @param gyyj output field with DG rep of Jg^yy
 * @param gxzj output field with DG rep of Jg^xz
 * @param eps2 output field with DG rep of eps2 = Jg^33 - J/g_33
 */
void gkyl_rz_calc_derived_geo_advance(const gkyl_rz_calc_derived_geo *up,
  const struct gkyl_range *crange, struct gkyl_array *g_ij, struct gkyl_array *bmag,
  struct gkyl_array *jacobgeo, struct gkyl_array *jacobgeo_inv, struct gkyl_array *gij,
  struct gkyl_array *b_i, struct gkyl_array *cmag, struct gkyl_array *jacobtot,
  struct gkyl_array *jacobtot_inv, struct gkyl_array *gxxj, struct gkyl_array *gxyj,
  struct gkyl_array *gyyj, struct gkyl_array *gxzj, struct gkyl_array *eps2);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_rz_calc_derived_geo_release(gkyl_rz_calc_derived_geo *up);
