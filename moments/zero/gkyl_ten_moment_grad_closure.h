#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

struct gkyl_ten_moment_grad_closure_inp {
  const struct gkyl_rect_grid *grid; // grid on which to solve equations
  double k0; // inverse damping coefficient
  double cfl; // CFL number to use
  const struct gkyl_range *update_range;
  const struct gkyl_range *heat_flux_range;

  struct gkyl_comm *comm;
};

struct gkyl_ten_moment_grad_closure_status {
  int success; // 1 if step worked, 0 otherwise
  double dt_suggested; // suggested time-step
};

// Object type
typedef struct gkyl_ten_moment_grad_closure gkyl_ten_moment_grad_closure;

/**
 * Create new updater to update pressure tensor in 10 moment equations
 * using a symmetrized gradient-based closure, q_ijk ~ \partial_i [_i T_{jk}]
 * Returns RHS for accumulation in a forward Euler method.
 *
 * @param inp Input parameters to updater
 */
gkyl_ten_moment_grad_closure *gkyl_ten_moment_grad_closure_new(
  const struct gkyl_ten_moment_grad_closure_inp *inp);

/**
 * Compute RHS contribution from symmetrized gradient-based closure
 * in 10 moment system.  The update_rng MUST be a sub-range of the
 * range on which the array is defined.  That is, it must be either
 * the same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param gces Gradient closure updater object.
 * @param heat_flux_rng Range on which to compute heat flux tensor (cell nodes)
 * @param update_rng Range on which to compute update.
 * @param fluid Input array of fluid variables 
 * @param em Total EM variables (plasma + external)
 * @param cflrate CFL scalar rate (frequency: units of 1/[T]) 
 * @param heat_flux Array for storing intermediate computation of heat flux tensor (cell nodes)
 * @param rhs RHS output (NOTE: Returns RHS output of all nfluids)
 */
struct gkyl_ten_moment_grad_closure_status gkyl_ten_moment_grad_closure_advance(
  const gkyl_ten_moment_grad_closure *gces, const struct gkyl_range *heat_flux_range,
  const struct gkyl_range *update_range, const struct gkyl_array *fluid,
  const struct gkyl_array *em_tot, struct gkyl_array *cflrate, double dt,
  struct gkyl_array *heat_flux, struct gkyl_array *rhs);

/**
 * Delete updater.
 *
 * @param gces Updater to delete.
 */
void gkyl_ten_moment_grad_closure_release(gkyl_ten_moment_grad_closure *gces);
