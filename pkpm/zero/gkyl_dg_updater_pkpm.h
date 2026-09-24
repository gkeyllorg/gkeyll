#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_euler_pkpm.h>
#include <gkyl_dg_vlasov_pkpm.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wv_eqn.h>

// Object type
typedef struct gkyl_dg_updater_pkpm gkyl_dg_updater_pkpm;

// return type for pkpm kinetic equations and fluids equations timers
struct gkyl_dg_updater_pkpm_tm {
  double vlasov_tm; // total time spent in computing kinetic equations
  double fluid_tm; // total time spent in computing fluid equations
};

/**
 * Create new updater to update PKPM equations using hyper dg.
 * Computes both the kinetic equation updates and fluid equation (momentum) updates.
 *
 * @param conf_grid Configuration space grid object
 * @param phase_grid Phase space grid object
 * @param conf_basis Configuration space basis functions
 * @param phase_basis Phase space basis function
 * @param conf_range Configuration space range
 * @param phase_range Phase space range
 * @param is_zero_flux_dir True in directions with (lower and upper) zero flux BCs.
 * @param vlasov_pkpm_inp Input struct to pkpm vlasov operator (see gkyl_dg_vlasov_pkpm.h)
 * @param euler_pkpm_inp Input struct to pkpm fluid operator (see gkyl_dg_euler_pkpm.h)
 * @param use_gpu Boolean to determine whether struct objects are on host or device
 * 
 * @return New PKPM updater object
 */
gkyl_dg_updater_pkpm* gkyl_dg_updater_pkpm_new(const struct gkyl_rect_grid *conf_grid, const struct gkyl_rect_grid *phase_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const bool *is_zero_flux_dir, 
  struct gkyl_dg_vlasov_pkpm_auxfields *vlasov_pkpm_inp, struct gkyl_dg_euler_pkpm_auxfields *euler_pkpm_inp, 
  bool use_gpu);

/**
 * Compute RHS of DG update. The update_rng MUST be a sub-range of the
 * range on which the array is defined. That is, it must be either the
 * same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param pkpm pkpm updater object
 * @param update_rng Range on which to compute.
 * @param fIn Input to updater
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T])
 * @param rhs RHS output
 */
void gkyl_dg_updater_pkpm_advance(gkyl_dg_updater_pkpm *pkpm,
  const struct gkyl_range *update_phase_rng, const struct gkyl_range *update_conf_rng, 
  const struct gkyl_array* GKYL_RESTRICT fIn, const struct gkyl_array* GKYL_RESTRICT fluidIn, 
  struct gkyl_array* GKYL_RESTRICT cflrate_f, struct gkyl_array* GKYL_RESTRICT cflrate_fluid, 
  struct gkyl_array* GKYL_RESTRICT rhs_f, struct gkyl_array* GKYL_RESTRICT rhs_fluid);

/**
 * Return total time spent in PKPM kinetic equations and fluid equations
 *
 * @param pkpm Updater object
 * @return timers
 */
struct gkyl_dg_updater_pkpm_tm gkyl_dg_updater_pkpm_get_tm(const gkyl_dg_updater_pkpm *pkpm);

/**
 * Delete updater.
 *
 * @param pkpm Updater to delete.
 */
void gkyl_dg_updater_pkpm_release(gkyl_dg_updater_pkpm* pkpm);
