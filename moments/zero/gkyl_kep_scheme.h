#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wv_eqn.h>

// Object type
typedef struct gkyl_kep_scheme gkyl_kep_scheme;

struct gkyl_kep_scheme_inp {
  const struct gkyl_rect_grid *grid; // grid on which to solve equations
  const struct gkyl_wv_eqn *equation; // equation to solve

  int num_up_dirs; // number of update directions  
  int update_dirs[GKYL_MAX_DIM]; // directions to update
  double cfl; // CFL number to use

  bool use_hybrid_flux; // should we use shock detector for hybrid flux

  const struct gkyl_wave_geom *geom; // geometry  
};

/**
 * Create new updater to update fluid equations using KEP scheme
 *
 * @param inp Input parameters to construct KEP scheme
 */
gkyl_kep_scheme* gkyl_kep_scheme_new(const struct gkyl_kep_scheme_inp *inp);

/**
 * Compute RHS of update. The update_rng MUST be a sub-range of the
 * range on which the array is defined. That is, it must be either the
 * same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param kep KEP scheme updater object.
 * @param update_rng Range on which to compute.
 * @param fIn Input to updater
 * @param alpha Array to store the shock indicator
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T])
 * @param rhs RHS output
 */
void gkyl_kep_scheme_advance(const gkyl_kep_scheme *kep,
  const struct gkyl_range *update_rng, const struct gkyl_array *fIn,
  struct gkyl_array *alpha,
  struct gkyl_array *cflrate, struct gkyl_array *rhs);

/**
 * Compute an estimate of maximum stable time-step for given input
 * state 'qin'
 *
 * @param kep Updater object
 * @param qin Input to compute dt for
 * @return maximum stable time-step
 */
double gkyl_kep_scheme_max_dt(const gkyl_kep_scheme *kep, const struct gkyl_range *update_range,
  const struct gkyl_array *qin);

/**
 * Delete updater.
 *
 * @param kep Updater to delete.
 */
void gkyl_kep_scheme_release(gkyl_kep_scheme* kep);
