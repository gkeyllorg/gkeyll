#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov_position_map.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wv_eqn.h>

// Object type
typedef struct gkyl_dg_updater_fluid gkyl_dg_updater_fluid;

// return type for drag and diffusion timers
struct gkyl_dg_updater_fluid_tm {
  double fluid_tm; // time for fluid updates
};

// Input packaged as a struct
struct gkyl_dg_updater_fluid_inp {
  const struct gkyl_rect_grid *grid; // Configuration-space grid.
  const struct gkyl_basis *cbasis; // Configuration-space basis functions.
  const struct gkyl_range *conf_range; // Configuration-space range.
  const struct gkyl_wv_eqn *wv_eqn; // Wave equation object: selects the fluid equation.
  const struct gkyl_wave_geom *geom; // Wave geometry of the faces (Euler and isothermal Euler).
  const struct gkyl_vlasov_position_map
    *pos_map; // Position map (REQUIRED for Euler and isothermal Euler).
  // Auxiliary arrays: struct gkyl_dg_advection_auxfields, gkyl_dg_euler_auxfields or
  // gkyl_dg_canonical_pb_fluid_auxfields, by the equation.
  void *aux_inp;
  // Edges whose face carries the absorbing-wall flux of an Euler or isothermal Euler fluid instead
  // of the interface flux: [d] the lower edge of direction d, [d + cdim] the upper one. Set only on
  // the ranks whose range holds that edge of the global domain.
  bool wall_flux[2 * GKYL_MAX_CDIM];
  bool use_gpu; // bool to determine if on GPU.
};

/**
 * Create new updater to update fluid equations using hyper dg.
 * Supports advection, isothermal Euler, Euler and canonical Poisson-bracket fluids.
 *
 * @param inp Input parameters defined in gkyl_dg_updater_fluid_inp struct.
 * @return New fluid updater object
 */
gkyl_dg_updater_fluid *gkyl_dg_updater_fluid_inew(const struct gkyl_dg_updater_fluid_inp *inp);

/**
 * Acquire fluid equation object
 *
 * @param fluid fluid updater object
 * 
 * @return fluid equation object
 */
struct gkyl_dg_eqn *gkyl_dg_updater_fluid_acquire_eqn(const gkyl_dg_updater_fluid *fluid);

/**
 * Compute RHS of DG update. The update_rng MUST be a sub-range of the
 * range on which the array is defined. That is, it must be either the
 * same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param fluid fluid updater object
 * @param update_rng Range on which to compute.
 * @param fIn Input to updater
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T])
 * @param rhs RHS output
 */
void gkyl_dg_updater_fluid_advance(
  gkyl_dg_updater_fluid *fluid, const struct gkyl_range *update_rng,
  const struct gkyl_array *GKYL_RESTRICT fIn, struct gkyl_array *GKYL_RESTRICT cflrate,
  struct gkyl_array *GKYL_RESTRICT rhs
);

/**
 * Return total time spent in drag and diffusion terms
 *
 * @param fluid Updater object
 * @return timers
 */
struct gkyl_dg_updater_fluid_tm gkyl_dg_updater_fluid_get_tm(const gkyl_dg_updater_fluid *fluid);

/**
 * Delete updater.
 *
 * @param fluid Updater to delete.
 */
void gkyl_dg_updater_fluid_release(gkyl_dg_updater_fluid *fluid);
