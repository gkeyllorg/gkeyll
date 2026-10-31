#pragma once

#include <gkyl_gk_geometry.h>
/**
 * Create a new geometry object using mac2p. 
 *
 * @param geometry_inp geometry input struct containing grid, range, and other geo info
 */

struct gk_geometry *gkyl_gk_geometry_mapc2p_new(struct gkyl_gk_geometry_inp *geometry_inp);

/**
 * Check whether the lower and upper z boundaries of a mapc2p geometry are the same
 * poloidal location (same R,Z), as needed to compute deltats from mapc2p.
 *
 * @param geometry_inp geometry input struct (uses mapc2p, c2p_ctx, position_map and geo_grid)
 * @return true if (x,y,z_lo) and (x,y,z_up) map to the same R,Z.
 */
bool gkyl_gk_geometry_mapc2p_z_ends_match(const struct gkyl_gk_geometry_inp *geometry_inp);
