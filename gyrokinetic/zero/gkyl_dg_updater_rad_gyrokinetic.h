#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_velocity_map.h>

// Object type
typedef struct gkyl_dg_updater_collisions gkyl_dg_updater_collisions;

// return type for drag and diffusion timers
struct gkyl_dg_updater_rad_gyrokinetic_tm {
  double drag_tm, diff_tm; // time for drag and diffusion updates
};

/**
 * Create new updater to update rad equations using hyper dg.
 *
 * @param grid Phase-space grid object
 * @param conf_basis Configuration-space basis functions
 * @param phase_basis Phase-space basis function
 * @param phase_range Phase-space range
 * @param conf_range Configuration-space range
 * @param vel_map Velocity space mapping object.
 * @param aux_inp Void pointer to auxiliary fields. Void to be flexible to different auxfields structs
 * @param use_gpu Boolean to determine if gyrokinetic equation object is on device
 * @return Pointer to updater object for Radiation operator in Gyrokinetic equation
 */
struct gkyl_dg_updater_collisions* 
gkyl_dg_updater_rad_gyrokinetic_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, 
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_velocity_map *vel_map, void *aux_inp, bool use_gpu);


/**
 * Compute RHS of DG update. The update_rng MUST be a sub-range of the
 * range on which the array is defined. That is, it must be either the
 * same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param rad Radiation operator updater object
 * @param update_rng Range on which to compute.
 * @param fIn Input to updater
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T])
 * @param rhs RHS output
 */
void gkyl_dg_updater_rad_gyrokinetic_advance(struct gkyl_dg_updater_collisions *rad,
  const struct gkyl_range *update_rng, const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs);

/**
 * Return total time spent in drag and diffusion terms
 *
 * @param rad Updater object
 * @return timers
 */
struct gkyl_dg_updater_rad_gyrokinetic_tm gkyl_dg_updater_rad_gyrokinetic_get_tm(const struct gkyl_dg_updater_collisions *coll);

/**
 * Delete updater.
 *
 * @param rad Updater to delete.
 */
void gkyl_dg_updater_rad_gyrokinetic_release(struct gkyl_dg_updater_collisions* coll);
