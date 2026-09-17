#pragma once

#include <gkyl_array.h>
#include <gkyl_fpo_vlasov_coeff_recovery.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_basis.h>

/**
 * Compute the drag coefficient a_i for the drag term of the FPO 
 * as the velocity space gradient of an arbitrary input potential H.
 * Calculates all 3 vector components of a_i. 
 *
 * @param grid Grid (for getting cell spacing)
 * @param cbasis Basis functions in configuration space
 * @param pbasis Basis functions in phase space
 * @param phase_range Range to calculate gradient
 * @param fpo_h Input potential
 * @param fpo_dhdv_surf Surface expansion of derivative of analytic potential
 * @param fpo_drag_coeff Output array of drag coefficient
 * @param fpo_drag_coeff_surf Output array for surface expansion of recovered drag coefficient
*/
void gkyl_calc_fpo_drag_coeff_recovery(const struct gkyl_fpo_vlasov_coeff_recovery* coeff_recovery, 
  const struct gkyl_rect_grid *grid, 
  struct gkyl_basis pbasis, const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *gamma, 
  const struct gkyl_array *fpo_h, const struct gkyl_array* fpo_dhdv_surf, 
  struct gkyl_array *fpo_drag_coeff, struct gkyl_array *fpo_drag_coeff_surf, bool use_gpu);

/**
 *  Host-side wrapper for drag coefficient computation.
 */
void gkyl_calc_fpo_drag_coeff_recovery_cu(const struct gkyl_fpo_vlasov_coeff_recovery* coeff_recovery,
  const struct gkyl_rect_grid *grid, 
  struct gkyl_basis pbasis, const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *gamma,
  const struct gkyl_array *fpo_h, const struct gkyl_array* fpo_dhdv_surf, 
  struct gkyl_array *fpo_drag_coeff, struct gkyl_array *fpo_drag_coeff_surf);

/**
 * Compute the sign of the surface expansion of the drag coefficient
 * at quadrature points for use in flux calculation for FPO drag update.
 *
 * @param pbasis Basis functions in phase space
 * @param phase_range Range over phase space
 * @param fpo_drag_coeff_surf Output array for surface expansion of recovered drag coefficient
 * @param sgn_drag_coeff_surf Sign of drag coefficient at quadrature points
 * @param const_sgn_drag_coeff_surf Flag that is true if the sign of the drag coefficient is constant at all quadrature points on a given surface
*/
void gkyl_calc_fpo_sgn_drag_coeff(const struct gkyl_fpo_vlasov_coeff_recovery *coeff_recovery,
  struct gkyl_basis pbasis, const struct gkyl_range *phase_range, 
  struct gkyl_array* fpo_drag_coeff_surf, 
  struct gkyl_array* sgn_drag_coeff_surf, struct gkyl_array* const_sgn_drag_coeff_surf, bool use_gpu);

void gkyl_calc_fpo_sgn_drag_coeff_cu(const struct gkyl_fpo_vlasov_coeff_recovery *coeff_recovery,
  struct gkyl_basis pbasis, const struct gkyl_range *phase_range,
  struct gkyl_array* fpo_drag_coeff_surf, 
  struct gkyl_array* sgn_drag_coeff_surf, struct gkyl_array* const_sgn_drag_coeff_surf);
