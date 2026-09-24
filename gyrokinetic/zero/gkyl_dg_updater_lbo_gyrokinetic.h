#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_lbo_gyrokinetic_diff.h>
#include <gkyl_dg_lbo_gyrokinetic_drag.h>
#include <gkyl_eqn_type.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_velocity_map.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_dg_updater_collisions gkyl_dg_updater_collisions;

// return type for drag and diffusion timers
struct gkyl_dg_updater_lbo_gyrokinetic_tm {
  double drag_tm, diff_tm; // time for drag and diffusion updates
};

/**
 * Create new updater to update lbo equations using hyper dg.
 *
 * @param phase_grid Phase space grid object
 * @param conf_basis Configuration space basis functions
 * @param phase_basis Phase space basis function
 * @param conf_range Configuration space range
 * @param drag_inp Input struct to gyrokinetic drag operator (see gkyl_dg_lbo_gyrokinetic_drag.h) 
 * @param diff_inp Input struct to gyrokinetic diffusion operator (see gkyl_dg_lbo_gyrokinetic_diff.h) 
 * @param mass Species mass.
 * @param gk_geom Gyrokinetic geometry object.
 * @param vel_map Velocity space mapping object.
 * @param use_gpu Bool for whether updater is on host or device
 * @return New gyrokinetic LBO updater object
 */
struct gkyl_dg_updater_collisions* 
gkyl_dg_updater_lbo_gyrokinetic_new(const struct gkyl_rect_grid *phase_grid,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, const struct gkyl_range *conf_range,
  struct gkyl_dg_lbo_gyrokinetic_drag_auxfields *drag_inp, struct gkyl_dg_lbo_gyrokinetic_diff_auxfields *diff_inp,
  double mass, const struct gk_geometry *gk_geom, const struct gkyl_velocity_map *vel_map,
  bool use_gpu);

/**
 * Compute RHS of DG update. The update_rng MUST be a sub-range of the
 * range on which the array is defined. That is, it must be either the
 * same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param lbo gyrokinetic LBO updater object
 * @param update_rng Range on which to compute.
 * @param fIn Input to updater
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T])
 * @param rhs RHS output
 */
void gkyl_dg_updater_lbo_gyrokinetic_advance(struct gkyl_dg_updater_collisions *lbo,
  const struct gkyl_range *update_rng, const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs);

/**
 * Return total time spent in drag and diffusion terms
 *
 * @param lbo Updater object
 * @return timers
 */
struct gkyl_dg_updater_lbo_gyrokinetic_tm gkyl_dg_updater_lbo_gyrokinetic_get_tm(const struct gkyl_dg_updater_collisions *coll);

/**
 * Delete updater.
 *
 * @param lbo Updater to delete.
 */
void gkyl_dg_updater_lbo_gyrokinetic_release(struct gkyl_dg_updater_collisions* coll);
