#pragma once

#include <gkyl_array.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_basis.h>
#include <gkyl_vlasov_velocity_map.h>

/**
 * Construct the velocity-space radiation drag force on the velocity-space grid
 * given the input radiation_id, either GKYL_VM_COMPTON_RADIATION (rad = -1/t_cool*(p - p0)) 
 * or GKYL_VM_CURVATURE_RADIATION (rad = -1/t_cool*sign(p)*p^4 = -1/t_cool*fabs(p)*p^3). 
 * 
 * @param vel_grid Velocity-space grid. 
 * @param vel_basis Velocity-space basis. 
 * @param vel_range Velocity-space range. 
 * @param radiation_id Radiation model indicator for which velocity-space radiation drag force to project.
 * @param vel_map Velocity-space mapping object for computing radiation drag force on mapped grids. NULL => uniform velocity grid.
 * @param t_cool Input cooling time parameter.
 * @param p0 Input (four-) velocity to relax to.
 * @param rad Output array storing the velocity-space radiation drag force.
 * @param use_gpu bool to determine if on GPU
 */
void gkyl_dg_vlasov_calc_radiation(
  const struct gkyl_rect_grid *vel_grid, const struct gkyl_basis *vel_basis,
  const struct gkyl_range *vel_range, enum gkyl_vlasov_radiation_id radiation_id,
  const struct gkyl_vlasov_velocity_map *vel_map, double t_cool, double p0, struct gkyl_array *rad,
  bool use_gpu
);

/**
 * Construct the velocity-space radiation drag force on the velocity-space grid on
 * NV-GPU. See new() method for documentation.
 */
void gkyl_dg_vlasov_calc_radiation_cu(
  const struct gkyl_rect_grid *vel_grid, const struct gkyl_basis *vel_basis,
  const struct gkyl_range *vel_range, enum gkyl_vlasov_radiation_id radiation_id,
  const struct gkyl_vlasov_velocity_map *vel_map, double t_cool, double p0, struct gkyl_array *rad
);