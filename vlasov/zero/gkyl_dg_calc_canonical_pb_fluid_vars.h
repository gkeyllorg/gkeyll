#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_dg_calc_canonical_pb_fluid_vars gkyl_dg_calc_canonical_pb_fluid_vars;

/**
 * Create new updater to compute canonical Poisson Bracket variables needed in
 * certain fluid systems, such as incompressible Euler, Hasegawa-Mima, and Hasegawa-Wakatani.  
 * Methods compute:
 * 1. alpha_surf : Surface expansion of configuration space flux alpha = {z, phi}
 * where phi is the potential given by a Poisson solve on (one of the) evolved quantit(ies)
 * such as in incompressible Euler, where phi is given by grad^2 phi = zeta, zeta is the vorticity. 
 * 2. source updates for different fluid systems, such as the adiabatic coupling and turbulence drive
 * in Hasegawa-Wakatani, source = alpha*(phi - n) + {phi, n0}, where n0 is the background density driving the turbulence, 
 * with the option to eliminate the zonal components of the adiabatic coupling with an average in y of phi and n 
 * and thus solve the modified Hasegawa-Wakatani system of equations. 
 * 
 * @param conf_grid Configuration-space grid (for getting cell spacing and cell center)
 * @param conf_basis Configuration-space basis functions
 * @param conf_range Configuration-space range. 
 * @param conf_ext_range Configuration-space extended range. 
 * @param wv_eqn Wave equation object which contains information about specific fluid system.
 *               For example, how many equations are we solving (1 for Hasegawa-Mima, 2 for Hasegawa-Wakatani), 
 *               and what the form of the source updates are for that fluid system. 
 * @param use_gpu bool to determine if on GPU
 * @return New updater pointer.
 */
struct gkyl_dg_calc_canonical_pb_fluid_vars* 
gkyl_dg_calc_canonical_pb_fluid_vars_new(const struct gkyl_rect_grid *conf_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_range *conf_range, const struct gkyl_range *conf_ext_range, 
  const struct gkyl_wv_eqn *wv_eqn,  bool use_gpu);

/**
 * Create new updater to compute canonical_pb general geometry variables on
 * NV-GPU. See new() method for documentation.
 */
struct gkyl_dg_calc_canonical_pb_fluid_vars* 
gkyl_dg_calc_canonical_pb_fluid_vars_cu_dev_new(const struct gkyl_rect_grid *conf_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_range *conf_range, const struct gkyl_range *conf_ext_range, 
  const struct gkyl_wv_eqn *wv_eqn);

/**
 * Compute surface expansion of configuration space flux alpha = {z, phi}
 * 
 * Note: Each cell stores the surface expansion on the *lower* edge of the cell
 * @param up Updater for computing canonical PB fluid variables 
 * @param conf_range Configuration space range 
 * @param conf_ext_range Extended Configuration space range (so we obtain alpha_surf at all the needed surfaces)
 * @param phi Potential expansion in a cell
 * @param alpha_surf Output surface expansion in a cell on the *lower* edge in each direction 
 * @param sgn_alpha_surf Output sign(alpha) at quadrature points along a surface 
 * @param const_sgn_alpha Output boolean array for if sign(alpha) is a constant on the surface
 *                        If sign(alpha) is a constant, kernels are simpler and we exploit this fact.
 */
void gkyl_dg_calc_canonical_pb_fluid_vars_alpha_surf(struct gkyl_dg_calc_canonical_pb_fluid_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *conf_ext_range, 
  const struct gkyl_array* phi,
  struct gkyl_array* alpha_surf, struct gkyl_array* sgn_alpha_surf, struct gkyl_array* const_sgn_alpha);

/**
 * Compute source update to canonical PB fluid system, such as the adiabatic coupling and turbulence drive
 * in Hasegawa-Wakatani, source = alpha*(phi - n) + {phi, background_n_gradient} where {.,.} is the Poisson bracket. 
 * 
 * @param up Updater for computing canonical PB fluid variables 
 * @param conf_range Configuration space range 
 * @param background_n_gradient Background density gradient for driving turbulence
 * @param phi Potential expansion in a cell
 * @param fluid Input array of fluid variables 
 * @param rhs Output increment to fluid variables from sources
 */
void gkyl_canonical_pb_fluid_vars_source(struct gkyl_dg_calc_canonical_pb_fluid_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_array *background_n_gradient, const struct gkyl_array *phi, 
  const struct gkyl_array *fluid, struct gkyl_array *rhs);

/**
 * Delete pointer to updater to compute canonical PB fluid variables.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_calc_canonical_pb_fluid_vars_release(struct gkyl_dg_calc_canonical_pb_fluid_vars *up);

/**
 * Host-side wrappers for canonical PB fluid variable operations on device
 */

void gkyl_dg_calc_canonical_pb_fluid_vars_alpha_surf_cu(struct gkyl_dg_calc_canonical_pb_fluid_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *conf_ext_range, 
  const struct gkyl_array* phi,
  struct gkyl_array* alpha_surf, struct gkyl_array* sgn_alpha_surf, struct gkyl_array* const_sgn_alpha);

void gkyl_canonical_pb_fluid_vars_source_cu(struct gkyl_dg_calc_canonical_pb_fluid_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_array *background_n_gradient, const struct gkyl_array *phi, 
  const struct gkyl_array *fluid, struct gkyl_array *rhs);
