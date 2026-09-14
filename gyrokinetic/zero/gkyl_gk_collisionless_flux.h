#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_gk_bc_type.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_velocity_map.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_dg_geom.h>
#include <gkyl_gk_dg_geom.h>

// Object type
typedef struct gkyl_gk_collisionless_flux gkyl_gk_collisionless_flux;

/**
 * Create new updater to compute the surface expansion of phase space
 * flux due to the collisionless gyrokinetic terms.
 * 
 * @param phase_grid Phase space grid (for getting cell spacing and cell center)
 * @param conf_basis Configuration space basis functions
 * @param phase_basis Phase space basis functions
 * @param charge Species charge
 * @param mass Species mass
 * @param collless_type Type of collisionless terms.
 * @param no_by Whether to neglect the toroidal field (set b_y=0).
 * @param complete_em Whether it is meant to complete the RHS after partial EM update.
 * @param gk_geom Gyrokinetic geometry object.
 * @param gk_dg_geom DG geometry object.
 * @param vel_map Velocity space mapping object.
 * @param bctype_conf Type of conf-space BCs.
 * @param use_gpu bool to determine if on GPU
 * @return New updater pointer.
 */
struct gkyl_gk_collisionless_flux* 
gkyl_gk_collisionless_flux_new(const struct gkyl_rect_grid *phase_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, 
  const double charge, const double mass,
  enum gkyl_gk_collisionless_type collless_type,
  const bool no_by, const bool complete_em, const struct gk_geometry *gk_geom, 
  const struct gkyl_dg_geom *dg_geom, const struct gkyl_gk_dg_geom *gk_dg_geom, 
  const struct gkyl_velocity_map *vel_map, const enum gkyl_gyrokinetic_bc_type *bctype_conf, bool use_gpu);

/**
 * Compute surface expansion of phase space flux alpha
 * Computes the Poisson bracket of alpha = dz/dt = {z, H}, multiplies it by the upwindinded Jf,
 * and then evaluates the resulting polynomial expansion at a surface and projects the evaluated quantity 
 * onto the surface basis.
 * 
 * Note: Each cell stores the surface expansion on the *lower* edge of the cell
 * @param up Updater for computing gyrokinetic variables.
 * @param conf_range Configuration space range (should only be local range because geometry only defined on local range).
 * @param phase_range Phase space range.
 * @param conf_ext_range Extended configuration space range (so we obtain geo quantities at all the needed surfaces).
 * @param phase_ext_range Extended Phase space range (so we obtain alpha_surf at all the needed surfaces).
 * @param phi Electrostatic potential.
 * @param apar Parallel component of vector potential.
 * @param apardot Time derivative of parallel component of vector potential.
 * @param fin Distribution function.
 * @param flux_surf Output surface expansion in a cell on the *lower* edge in each direction.
 * @param clfrate Output CFL rate.
 */
void gkyl_gk_collisionless_flux_surf(struct gkyl_gk_collisionless_flux *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const struct gkyl_range *conf_ext_range, const struct gkyl_range *phase_ext_range, 
  const struct gkyl_array *phi, const struct gkyl_array *apar, const struct gkyl_array *apardot,
  const struct gkyl_array *fin, struct gkyl_array* flux_surf, struct gkyl_array *cflrate);

/**
 * Delete pointer to updater to compute gyrokinetic variables.
 *
 * @param up Updater to delete.
 */
void gkyl_gk_collisionless_flux_release(struct gkyl_gk_collisionless_flux *up);
