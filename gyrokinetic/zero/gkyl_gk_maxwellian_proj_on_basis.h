#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_velocity_map.h>

// Object type
typedef struct gkyl_gk_maxwellian_proj_on_basis gkyl_gk_maxwellian_proj_on_basis;

// input packaged as a struct
struct gkyl_gk_maxwellian_proj_on_basis_inp {
  const struct gkyl_rect_grid *phase_grid; // Phase-space grid on which to compute moments
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions
  const struct gkyl_basis *phase_basis; // Phase-space basis functions
  const struct gkyl_range *conf_range; // Configuration-space range
  const struct gkyl_range *conf_range_ext; // Extended configuration-space range (for internal memory allocations)
  const struct gkyl_range *vel_range; // Velocity space range
  const struct gk_geometry *gk_geom; // Geometry object.
  const struct gkyl_velocity_map *vel_map; // Velocity space mapping object.
  enum gkyl_quad_type quad_type; // Quadrature type (default: Gauss-Legendre, see gkyl_eqn_type.h)
  double mass; // Mass factor
  bool bimaxwellian; // Bool for whether we are projecting a bi-Maxwellian instead of a Maxwellian.
  bool divide_jacobgeo; // Bool for whether to divide out the conf-space Jacobian from density.
  bool use_gpu; // bool for gpu useage
};

/**
 * Create new updater to project the Gyrokinetic Maxwellian distribution function 
 * onto basis functions. Also supports projecting a bi-Maxwellian distribution function. 
 * Free using gkyl_gk_maxwellian_proj_on_basis_release method.
 *
 * @param inp Input parameters defined in gkyl_gk_maxwellian_proj_on_basis_inp struct.
 * @return New updater pointer.
 */
struct gkyl_gk_maxwellian_proj_on_basis* 
gkyl_gk_maxwellian_proj_on_basis_inew(const struct gkyl_gk_maxwellian_proj_on_basis_inp *inp);

/**
 * Compute projection of Maxwellian (or bi-Maxwellian) distribution on basis. 
 * This method takes the moments as a single array moms_maxwellian = (n, upar, T/m)
 * or (n, upar, Tpar/m, Tperp/m) to compute the projection of the distribution function.
 * Further note: We correct the density moment of the projected distribution 
 * before we return f_maxwellian.
 *
 * @param up Project on basis updater to run
 * @param phase_rng Phase-space range
 * @param conf_rng Configuration-space range
 * @param moms_maxwellian Moments for computing distribution function, e.g., (n, upar, T/m)
 * @param use_jacobtot Bool to determine if we are using the total Jacobian or only
 *                     the velocity-space Jacobian in the exponential weighting.
 * @param f_maxwellian Output Maxwellian (or bi-Maxwellian) distribution function
 */
void gkyl_gk_maxwellian_proj_on_basis_advance(gkyl_gk_maxwellian_proj_on_basis *up,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *moms_maxwellian, bool use_jacobtot, 
  struct gkyl_array *f_maxwellian);

/**
 * Host-side wrapper for geometry variables (bmag and jacobtot) at quadrature points
 */
void gkyl_gk_maxwellian_proj_on_basis_geom_quad_vars_cu(gkyl_gk_maxwellian_proj_on_basis *up, 
  const struct gkyl_range *conf_range, 
  const struct gkyl_array *bmag, const struct gkyl_array *jacobtot);

/**
 * Host-side wrapper for projection of Maxwellian distribution function on device
 */
void gkyl_gk_maxwellian_proj_on_basis_advance_cu(gkyl_gk_maxwellian_proj_on_basis *up,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *moms_maxwellian, bool use_jacobtot, 
  struct gkyl_array *f_maxwellian);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_gk_maxwellian_proj_on_basis_release(gkyl_gk_maxwellian_proj_on_basis* up);
