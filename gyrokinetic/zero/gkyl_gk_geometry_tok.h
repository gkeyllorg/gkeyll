#pragma once

#include <gkyl_gk_geometry.h>
/**
 * Create a new geometry object using tokamak input (efit)
 *
 * @param geometry_inp geometry input struct containing grid, range, and other geo info
 */
struct gk_geometry* gkyl_gk_geometry_tok_new(struct gkyl_gk_geometry_inp *geometry_inp);

/*
 * Set the lower and upper z-direction grid extents of a tokamak block
 *
 * @param efit info efit input info
 * @param grid_info computational grid input info
 * @param theta_lo on output lower extent in z direction
 * @param theta_up on output upper extent in z direction
 * */
void
gkyl_gk_geometry_tok_set_grid_extents(struct gkyl_efit_inp efit_info, 
  struct gkyl_tok_geo_grid_inp grid_info, double *theta_lo, double *theta_up);



