#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_velocity_map.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_dg_updater_gyrokinetic gkyl_dg_updater_gyrokinetic;

// return type for gyrokinetic timers
struct gkyl_dg_updater_gyrokinetic_tm {
  double gyrokinetic_tm; // time for gyrokinetic updates
};


/**
 * Create new updater to update gyrokinetic equations using hyper dg.
 *
 * @param grid Grid object
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis function
 * @param conf_range Configuration space range
 * @param phase_range Phase space range
 * @param is_zero_flux_bc[2*GKYL_MAX_DIM] True for boundaries with zero flux BCs.
 * @param charge Species charge
 * @param mass Species mass
 * @param collless_type Type of collisionless terms.
 * @param gk_geom Geometry struct 
 * @param vel_map Velocity space mapping object.
 * @param aux_inp Void pointer to auxiliary fields. Void to be flexible to different auxfields structs
 * @param use_gpu Boolean to determine if gyrokinetic equation object is on device
 * @return Pointer to updater object for Gyrokinetic equation
 */
gkyl_dg_updater_gyrokinetic* gkyl_dg_updater_gyrokinetic_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const bool *is_zero_flux_bc, double charge, double mass,
  enum gkyl_gk_collisionless_type collless_type, const struct gk_geometry *gk_geom,
  const struct gkyl_velocity_map *vel_map, void *aux_inp, bool use_gpu);

/**
 * Acquire gyrokinetic equation object.
 *
 * @param up gyrokinetic updater object.
 * 
 * @return gyrokinetic equation object.
 */
struct gkyl_dg_eqn* 
gkyl_dg_updater_gyrokinetic_acquire_eqn(const gkyl_dg_updater_gyrokinetic* up);

/**
 * Compute RHS of DG update. The update_rng MUST be a sub-range of the
 * range on which the array is defined. That is, it must be either the
 * same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param gyrokinetic gyrokinetic updater object
 * @param update_rng Range on which to compute.
 * @param fIn Input to updater
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T])
 * @param rhs RHS output
 */
void gkyl_dg_updater_gyrokinetic_advance(gkyl_dg_updater_gyrokinetic *gyrokinetic,
  const struct gkyl_range *update_rng, const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs);

/**
 * Return total time spent in gyrokinetic equation
 *
 * @param gyrokinetic Updater object
 * @return timers
 */
struct gkyl_dg_updater_gyrokinetic_tm gkyl_dg_updater_gyrokinetic_get_tm(const gkyl_dg_updater_gyrokinetic *gyrokinetic);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_updater_gyrokinetic_release(gkyl_dg_updater_gyrokinetic* up);
