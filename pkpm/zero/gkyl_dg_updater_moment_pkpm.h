#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type (defined in gkyl_dg_updater_moment_priv.h)
typedef struct gkyl_dg_updater_moment gkyl_dg_updater_moment;

// return type for moment timers (defined in gkyl_dg_updater_moment_priv.h)
typedef struct gkyl_dg_updater_moment_tm gkyl_dg_updater_moment_tm;

// Input packaged as a struct
struct gkyl_dg_updater_moment_pkpm_inp {
  const struct gkyl_rect_grid *phase_grid; // Phase-space grid.
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions (p=1 serendipity).
  const struct gkyl_basis *phase_basis; // Phase-space basis functions (p=1 x p=2 hybrid).
  const struct gkyl_range *conf_range; // Configuration-space range (not read by the updater).
  const struct gkyl_range *vel_range; // Velocity-space range (not read by the updater).
  double mass; // Mass of species (the moments are mass weighted).
  bool is_diag; // bool to determine if computing the diagnostic moments.
  bool use_gpu; // bool to determine if on GPU.
};

/**
 * Create new updater to compute moments of distribution function
 * in parallel-kinetic-perpendicular-moment (PKPM) model.
 *
 * 
 * @param inp Input parameters defined in gkyl_dg_updater_moment_pkpm_inp struct.
 * @return New moment updater object
 */
struct gkyl_dg_updater_moment *gkyl_dg_updater_moment_pkpm_inew(
  const struct gkyl_dg_updater_moment_pkpm_inp *inp
);

/**
 * Acquire moment object
 *
 * @param moment moment updater object
 * 
 * @return moment type object
 */
struct gkyl_mom_type *gkyl_dg_updater_moment_pkpm_acquire_type(
  const struct gkyl_dg_updater_moment *moment
);

/**
 * Acquire number of moments
 *
 * @param moment moment updater object
 * 
 * @return number of moments
 */
int gkyl_dg_updater_moment_pkpm_num_mom(const struct gkyl_dg_updater_moment *moment);

/**
 * Compute moment. The update_phase_rng and update_conf_rng MUST each be a sub-range of the
 * range on which the corresponding array is defined. That is, it must be either the same range
 * as the array range, or one created using the gkyl_sub_range_init method.
 *
 * @param moment moment updater object
 * @param update_phase_rng Phase space range on which to compute.
 * @param update_conf_rng Configuration space range on which to compute.
 * @param fIn Input to updater
 * @param mout Output moment
 */
void gkyl_dg_updater_moment_pkpm_advance(
  struct gkyl_dg_updater_moment *moment, const struct gkyl_range *update_phase_rng,
  const struct gkyl_range *update_conf_rng, const struct gkyl_array *GKYL_RESTRICT fIn,
  struct gkyl_array *GKYL_RESTRICT mout
);

/**
 * Return total time spent in computing moments
 *
 * @param moment Updater object
 * @return timers
 */
struct gkyl_dg_updater_moment_tm gkyl_dg_updater_moment_pkpm_get_tm(
  const struct gkyl_dg_updater_moment *moment
);

/**
 * Delete updater.
 *
 * @param moment Updater to delete.
 */
void gkyl_dg_updater_moment_pkpm_release(struct gkyl_dg_updater_moment *moment);