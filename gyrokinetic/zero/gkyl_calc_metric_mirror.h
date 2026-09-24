#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_mirror_grid_gen.h>

// Object type
typedef struct gkyl_calc_metric_mirror gkyl_calc_metric_mirror;

/**
 * Create new updater to compute the metric coefficients
 *
 * @param cbasis Basis object (configuration space).
 * @param grid configuration space grid.
 * @param local, _ext computational local Config-space range and extended range.
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_calc_metric_mirror* gkyl_calc_metric_mirror_new(const struct gkyl_basis *cbasis,
  const struct gkyl_rect_grid *grid, const struct gkyl_range *local, 
  const struct gkyl_range *local_ext, bool use_gpu);

/**
 * Calculate metric coefficients and jacobian at nodes
 * Using the explicit in Eq. 66-73 of the GK coordinates document
 * Then convert to modal
 *
 * @param up calc_metric_mirror updater object.
 * @param gk_geom gk_geom objet to be populated
 * @param mirror_grid mirr_grid object with coordinates, tangents, duals, Jc, and B at nodes
 */
void gkyl_calc_metric_mirror_advance(gkyl_calc_metric_mirror *up, struct gk_geometry *gk_geom, struct gkyl_mirror_grid_gen *mirror_grid);

/**
 * calculate metric coefficients and jacobian at 
 * gauss-legendre quadrature nodes
 * Using the explicit in Eq. 66-73 of the GK coordinates document
 * Then convert to modal
 *
 * @param up calc_metric_mirror updater object.
 * @param gk_geom gk_geom objet to be populated
 * @param mirror_grid mirr_grid object with coordinates, tangents, duals, Jc, and B at nodes
 */
void gkyl_calc_metric_mirror_advance_interior(gkyl_calc_metric_mirror *up, struct gk_geometry *gk_geom, struct gkyl_mirror_grid_gen *mirror_grid);

/**
 * Use finite differences to calculate metric coefficients and jacobian at 
 * surface gauss-legendre quadrature nodes
 * Using the explicit in Eq. 66-73 of the GK coordinates document
 * Then convert to modal
 *
 * @param up calc_metric_mirror updater object.
 * @param dir direction of surface
 * @param gk_geom gk_geom objet to be populated
 * @param mirror_grid mirr_grid object with coordinates, tangents, duals, Jc, and B at nodes
 */
void gkyl_calc_metric_mirror_advance_surface( gkyl_calc_metric_mirror *up, int dir, struct gk_geometry* gk_geom, struct gkyl_mirror_grid_gen *mirror_grid);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_calc_metric_mirror_release(gkyl_calc_metric_mirror* up);
