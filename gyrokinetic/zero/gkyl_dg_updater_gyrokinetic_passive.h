#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_velocity_map.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Fully separate type from gkyl_dg_updater_gyrokinetic.
typedef struct gkyl_dg_updater_gyrokinetic_passive gkyl_dg_updater_gyrokinetic_passive;

// Timer struct.
struct gkyl_dg_updater_gyrokinetic_passive_tm {
  double gyrokinetic_passive_tm;
};

/**
 * Create a new updater to advance the passive gyrokinetic equation using
 * hyper_dg.  Only conf-space directions are updated; no vpar surface kernel
 * is used.
 *
 * @param grid Phase space grid.
 * @param cbasis Configuration space basis functions.
 * @param pbasis Phase-space basis functions.
 * @param conf_range Configuration space range.
 * @param phase_range Phase space range.
 * @param is_zero_flux_bc[2*GKYL_MAX_DIM] True for zero-flux boundaries.
 * @param charge Species charge.
 * @param mass Species mass.
 * @param gk_geom Geometry struct.
 * @param vel_map Velocity space mapping object.
 * @param aux_inp Void pointer to gkyl_dg_gyrokinetic_passive_auxfields.
 * @param use_gpu Boolean to determine if on device.
 * @return Pointer to new updater object.
 */
gkyl_dg_updater_gyrokinetic_passive*
gkyl_dg_updater_gyrokinetic_passive_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis,
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const bool *is_zero_flux_bc,
  const double charge, const double mass,
  const struct gk_geometry *gk_geom,
  const struct gkyl_velocity_map *vel_map,
  void *aux_inp, bool use_gpu);

/**
 * Acquire passive gyrokinetic equation object.
 *
 * @param up Passyive gyrokinetic updater object.
 * 
 * @return Passive gyrokinetic equation object.
 */
struct gkyl_dg_eqn* 
gkyl_dg_updater_gyrokinetic_passive_acquire_eqn(const gkyl_dg_updater_gyrokinetic_passive* up);

/**
 * Compute the RHS DG update.
 *
 * @param up Passive gyrokinetic updater object.
 * @param update_rng Range on which to compute.
 * @param fIn Input distribution function.
 * @param cflrate CFL rate array (units of 1/[T]).
 * @param rhs RHS output.
 */
void gkyl_dg_updater_gyrokinetic_passive_advance(gkyl_dg_updater_gyrokinetic_passive *up,
  const struct gkyl_range *update_rng, const struct gkyl_array *fIn,
  struct gkyl_array *cflrate, struct gkyl_array *rhs);

/**
 * Return total time spent in the passive gyrokinetic equation.
 *
 * @param up Updater object.
 * @return Timers.
 */
struct gkyl_dg_updater_gyrokinetic_passive_tm
gkyl_dg_updater_gyrokinetic_passive_get_tm(const gkyl_dg_updater_gyrokinetic_passive *up);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_updater_gyrokinetic_passive_release(gkyl_dg_updater_gyrokinetic_passive *up);
