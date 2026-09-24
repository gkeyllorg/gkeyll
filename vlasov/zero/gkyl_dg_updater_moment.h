#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_mom_type.h>

// Object type (defined in gkyl_dg_updater_moment_priv.h)
typedef struct gkyl_dg_updater_moment gkyl_dg_updater_moment;

// return type for moment timers (defined in gkyl_dg_updater_moment_priv.h)
typedef struct gkyl_dg_updater_moment_tm gkyl_dg_updater_moment_tm;

/**
 * Create new updater to compute moments of distribution function.
 * Supports Vlasov-Maxwell and special relativistic Vlasov-Maxwell,
 *
 * @param grid Grid object
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis function
 * @param conf_range Config space range
 * @param vel_range Velocity space range
 * @param phase_range phase space range
 * @param model_id Enum identifier for model type (e.g., SR, see gkyl_eqn_type.h)
 * @param aux_inp Void pointer to auxiliary fields. Void to be flexible to different auxfields structs
 * @param mom_type Name of moment.
 * @param is_integrated Boolean for if the moment is an integrated moment
 * @param mass Mass of species 
 * @param use_gpu Boolean to determine whether struct objects are on host or device
 * 
 * @return New moment updater object
 */
struct gkyl_dg_updater_moment*
gkyl_dg_updater_moment_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis, 
  const struct gkyl_range *conf_range, const struct gkyl_range *vel_range,
  const struct gkyl_range *phase_range,
  enum gkyl_model_id model_id, void *aux_inp, 
  enum gkyl_distribution_moments mom_type, bool is_integrated, bool use_gpu);

/**
 * Acquire moment object
 *
 * @param moment moment updater object
 * 
 * @return moment type object
 */
struct gkyl_mom_type* 
gkyl_dg_updater_moment_acquire_type(const struct gkyl_dg_updater_moment* moment);

/**
 * Acquire number of moments
 *
 * @param moment moment updater object
 * 
 * @return number of moments
 */
int 
gkyl_dg_updater_moment_num_mom(const struct gkyl_dg_updater_moment* moment);

/**
 * Compute moment. The update_phase_rng and update_conf_rng MUST be a sub-range of the
 * be a sub-range of the range on which the array is defined. 
 * That is, it must be either the same range as the array range, or one created using the
 * or one created using the gkyl_sub_range_init method.
 *
 * @param moment moment updater object
 * @param update_phase_rng Phase space range on which to compute.
 * @param update_conf_rng Configuration space range on which to compute.
 * @param fIn Input to updater
 * @param mout Output moment
 */
void
gkyl_dg_updater_moment_advance(struct gkyl_dg_updater_moment *moment,
  const struct gkyl_range *update_phase_rng, const struct gkyl_range *update_conf_rng,
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT mout);

/**
 * Return total time spent in computing moments
 *
 * @param moment Updater object
 * @return timers
 */
struct gkyl_dg_updater_moment_tm gkyl_dg_updater_moment_get_tm(const struct gkyl_dg_updater_moment *moment);

/**
 * Delete updater.
 *
 * @param moment Updater to delete.
 */
void gkyl_dg_updater_moment_release(struct gkyl_dg_updater_moment* moment);
