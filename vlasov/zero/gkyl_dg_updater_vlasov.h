#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_dg_updater_vlasov_timers.h>

// Object type
typedef struct gkyl_dg_updater_vlasov gkyl_dg_updater_vlasov;

/**
 * Create new updater to update vlasov equations using hyper dg.
 * Supports Vlasov-Maxwell, Vlasov-Poisson (with and without vector potential A)
 * and special relativistic Vlasov-Maxwell
 *
 * @param grid Grid object
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase space basis function
 * @param conf_range Configuration space range
 * @param vel_range Velocity space range
 * @param phase_range Phase space range
 * @param is_zero_flux_bc True for zero flux boundaries.
 * @param model_id Enum identifier for model type (e.g., SR, General Geometry, see gkyl_eqn_type.h)
 * @param field_id Enum identifier for field type (e.g., Maxwell's, Poisson, see gkyl_eqn_type.h)
 * @param aux_inp Void pointer to auxiliary fields. Void to be flexible to different auxfields structs
 * @param use_gpu Boolean to determine whether struct objects are on host or device
 * 
 * @return New vlasov updater object
 */
gkyl_dg_updater_vlasov* gkyl_dg_updater_vlasov_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis, 
  const struct gkyl_range *conf_range, const struct gkyl_range *vel_range, const struct gkyl_range *phase_range,
  const bool *is_zero_flux_bc, enum gkyl_model_id model_id, enum gkyl_field_id field_id, void *aux_inp, bool use_gpu);

/**
 * Acquire Vlasov equation object
 *
 * @param vlasov Vlasov updater object
 * 
 * @return Vlasov equation object
 */
struct gkyl_dg_eqn* 
gkyl_dg_updater_vlasov_acquire_eqn(const gkyl_dg_updater_vlasov* vlasov);

/**
 * Compute RHS of DG update. The update_rng MUST be a sub-range of the
 * range on which the array is defined. That is, it must be either the
 * same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param vlasov vlasov updater object
 * @param update_rng Range on which to compute.
 * @param fIn Input to updater
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T])
 * @param rhs RHS output
 */
void gkyl_dg_updater_vlasov_advance(gkyl_dg_updater_vlasov *vlasov,
  const struct gkyl_range *update_rng, const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs);

/**
 * Return total time spent in vlasov equation
 *
 * @param vlasov Updater object
 * @return timers
 */
struct gkyl_dg_updater_vlasov_tm gkyl_dg_updater_vlasov_get_tm(const gkyl_dg_updater_vlasov *vlasov);

/**
 * Delete updater.
 *
 * @param vlasov Updater to delete.
 */
void gkyl_dg_updater_vlasov_release(gkyl_dg_updater_vlasov* vlasov);
