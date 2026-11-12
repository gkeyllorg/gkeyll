#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>

// Object type
typedef struct gkyl_bgk_collisions gkyl_bgk_collisions;

/**
 * Create new updater to compute the increment due to a BGK
 * collision operator
 *   nu*f_M - nu*f)
 * where nu*f_M=sum_r nu_sr*f_Msr and nu=sum_r nu_sr, in
 * order to support multispecies collisions. The quantities
 * nu*f_M and nu must be computed elsewhere.
 *
 * @param cbasis Basis object (configuration space).
 * @param pbasis Basis object (phase space).
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_bgk_collisions *gkyl_bgk_collisions_new(
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis, bool use_gpu
);

/**
 * Set nu to zero for all points where density or temperature cell averages are negative/zero.
 * 
 * @param up BGK collision frequency updater object.
 * @param crange Config-space range.
 * @param marr Array of moments (density, momentum, energy).
 * @param nu_input Input collision frequency (nu).
 * @param actual_nu Output collision frequency (nu), corrected to be zero where density or temperature are negative/zero.
 */
void gkyl_bgk_collisions_correct_nu(
  const gkyl_bgk_collisions *up, const struct gkyl_range *crange, const struct gkyl_array *marr,
  const struct gkyl_array *nu_input, struct gkyl_array *actual_nu
);

/**
 * Advance BGK operator (compute the BGK contribution to df/dt).
 *
 * @param up BGK collision frequency updater object.
 * @param crange Config-space range.
 * @param prange Phase-space range.
 * @param nu Sum of collision frequencies.
 * @param nufM Sum of collision frequencies times their respective Maxwellian.
 * @param fin Input distribution function.
 * @param implicit_step  boolean of wheather or not to take an implicit step
 * @param dt timestep size used only in the implicit calc
 * @param out BGK contribution to df/dt.
 * @param cflfreq Output CFL frequency.
 */
void gkyl_bgk_collisions_advance(
  const gkyl_bgk_collisions *up, const struct gkyl_range *crange, const struct gkyl_range *prange,
  const struct gkyl_array *nu, const struct gkyl_array *nufM, const struct gkyl_array *fin,
  bool implicit_step, double dt, struct gkyl_array *out, struct gkyl_array *cflfreq
);

/**
 * Delete updater.
 *
 * @param pob Updater to delete.
 */
void gkyl_bgk_collisions_release(gkyl_bgk_collisions *up);
