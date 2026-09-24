#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

// Identifiers for different Braginskii types
enum gkyl_braginskii_type
{
  NONE = 0,
  GKYL_BRAG_MAG = 1 << 0,
  GKYL_BRAG_VISC = 1 << 1,
  GKYL_BRAG_HEATFLUX = 1 << 2,

  GKYL_BRAG_UNMAG_FULL = GKYL_BRAG_VISC | GKYL_BRAG_HEATFLUX,
  GKYL_BRAG_MAG_FULL = GKYL_BRAG_MAG | GKYL_BRAG_UNMAG_FULL
};

struct gkyl_moment_braginskii_data {
  enum gkyl_eqn_type type_eqn; // equation type
  enum gkyl_braginskii_type type_brag; // which Braginskii equations
  double charge; // species charge
  double mass; // species mass
  double p_fac; // factor for obtaining pressure for Braginskii coefficients
};

struct gkyl_moment_braginskii_inp {
  const struct gkyl_rect_grid *grid; // grid on which to solve equations
  int nfluids; // number of fluids
  struct gkyl_moment_braginskii_data param[GKYL_MAX_SPECIES]; // species data
  double epsilon0; // permittivity of free space
  double coll_fac; // constant multiplicative factor for collision time to increase or decrease collisionality
};

// Object type
typedef struct gkyl_moment_braginskii gkyl_moment_braginskii;

/**
 * Create new updater to update fluid equations with Braginskii 
 * transport terms (viscosity, heat conduction, etc.).
 * Returns RHS for accumulation in a forward Euler method.
 *
 * @param inp Input parameters to updater
 */
gkyl_moment_braginskii* gkyl_moment_braginskii_new(struct gkyl_moment_braginskii_inp inp);

/**
 * Compute RHS contribution from Braginskii transport terms in the
 * multi-fluid system.  The update_rng MUST be a sub-range of the
 * range on which the array is defined.  That is, it must be either
 * the same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param bes Braginskii updater object.
 * @param brag_vars_rng Range on which to compute Braginskii variables (cell nodes)
 * @param update_rng Range on which to compute update.
 * @param fluid Input array of fluid variables (array size: nfluids)
 * @param em Total EM variables (plasma + external)
 * @param cflrate CFL scalar rate (frequency: units of 1/[T]) (array size: nfluids)
 * @param brag_vars Array for storing intermediate computation of Braginskii variables (cell nodes)
 * @param rhs RHS output (NOTE: Returns RHS output of all nfluids)
 */

void gkyl_moment_braginskii_advance(const gkyl_moment_braginskii *bes,   
  struct gkyl_range brag_vars_range, struct gkyl_range update_range,
  struct gkyl_array *fluid[GKYL_MAX_SPECIES], const struct gkyl_array *em_tot,
  struct gkyl_array *cflrate[GKYL_MAX_SPECIES], 
  struct gkyl_array *brag_vars[GKYL_MAX_SPECIES], struct gkyl_array *rhs[GKYL_MAX_SPECIES]);

/**
 * Delete updater.
 *
 * @param bes Updater to delete.
 */
void gkyl_moment_braginskii_release(gkyl_moment_braginskii *bes);
