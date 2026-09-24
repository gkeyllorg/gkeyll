#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h> 

// Object type
typedef struct gkyl_gk_maxwellian_moments gkyl_gk_maxwellian_moments;

// input packaged as a struct
struct gkyl_gk_maxwellian_moments_inp {
  const struct gkyl_rect_grid *phase_grid; // Phase-space grid on which to compute moments.
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions.
  const struct gkyl_basis *phase_basis; // Phase-space basis functions.
  const struct gkyl_range *conf_range; // Configuration-space range.
  const struct gkyl_range *conf_range_ext; // Extended configuration-space range (for internal memory allocations)
  const struct gk_geometry *gk_geom; // Geometry object.
  const struct gkyl_velocity_map *vel_map; // Velocity mapping object.
  double mass; // Mass factor. 
  bool divide_jacobgeo; // Whether to divide out the conf-space Jacobian from density.
  bool use_gpu; // bool for gpu useage
};

/**
 * Create new updater to compute the moments for the equivalent Maxwellian distribution
 * function in the gyrokinetic system of equations. 
 * Updater always returns (n, u_par, T/m) where T/m is the stationary frame 
 * temperature/mass (the frame moving at velocity u_par).
 * **Note**: if divide_jacobgeo = true, updater divides out the 
 * configuration-space Jacobian so the moments returned are the actual moments.
 * 
 * @param inp Input parameters defined in gkyl_gk_maxwellian_moments_inp struct.
 * @return New updater pointer.
 */
struct gkyl_gk_maxwellian_moments*
gkyl_gk_maxwellian_moments_inew(const struct gkyl_gk_maxwellian_moments_inp *inp);

/**
 * Compute the density moments of an arbitrary distribution function for the equivalent 
 * Maxwellian distribution function for the gyrokinetic system of equations.
 * Computes n, and if divide_jacobgeo = true *the configuration-space Jacobian factor is divided out*.
 *
 * @param up Maxwellian moments updater
 * @param phase_local Phase-space range on which to compute moments.
 * @param conf_local Configuration-space range on which to compute moments.
 * @param fin Input distribution function
 * @param density_out Output density
 */
void gkyl_gk_maxwellian_density_moment_advance(struct gkyl_gk_maxwellian_moments *up, 
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local, 
  const struct gkyl_array *fin, struct gkyl_array *density_out);

/**
 * Compute the moments of an arbitrary distribution function for the equivalent 
 * Maxwellian distribution function for the gyrokinetic system of equations.
 * Computes (n, u_par, T/m) where T/m is the stationary frame 
 * temperature/mass (the frame moving at velocity u_par).
 * **Note**: If divide_jacobgeo = true, updater divides out the 
 * configuration-space Jacobian so the moments returned are the actual moments.
 * 
 * @param up Maxwellian moments updater
 * @param phase_local Phase-space range on which to compute moments.
 * @param conf_local Configuration-space range on which to compute moments.
 * @param fin Input distribution function
 * @param moms_out Output Maxwellian moments (n, u_par, T/m)
 */
void gkyl_gk_maxwellian_moments_advance(struct gkyl_gk_maxwellian_moments *up, 
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local, 
  const struct gkyl_array *fin, struct gkyl_array *moms_out);

/**
 * Compute the moments of an arbitrary distribution function for the equivalent 
 * *Bi*-Maxwellian distribution function for the gyrokinetic system of equations.
 * Computes (n, u_par, T_par/m, T_perp/m) where T_par/m and T_perp/m are the stationary frame 
 * parallel and perpendicular temperature/mass (the frame moving at velocity u_par).
 * **Note**: If divide_jacobgeo = true, updater divides out the 
 * configuration-space Jacobian so the moments returned are the actual moments.
 * 
 * @param up Maxwellian moments updater
 * @param phase_local Phase-space range on which to compute moments.
 * @param conf_local Configuration-space range on which to compute moments.
 * @param fin Input distribution function
 * @param moms_out Output Maxwellian moments (n, u_par, T_par/m, T_perp/m)
 */
void gkyl_gk_bimaxwellian_moments_advance(struct gkyl_gk_maxwellian_moments *up, 
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local, 
  const struct gkyl_array *fin, struct gkyl_array *moms_out);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_gk_maxwellian_moments_release(gkyl_gk_maxwellian_moments* up);
