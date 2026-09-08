#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>

// Type for storing preallocated memory
typedef struct gkyl_gyrokinetic_cross_prim_moms_bgk gkyl_gyrokinetic_cross_prim_moms_bgk;

/**
 * Allocate memory for use in cross moments calculation.
 *
 * @param phase_basis Phase space basis functions.
 * @param conf_basis Configuration space basis functions.
 * @param use_gpu Boolian to determine if on GPU.
 */
gkyl_gyrokinetic_cross_prim_moms_bgk *gkyl_gyrokinetic_cross_prim_moms_bgk_new(
  const struct gkyl_basis *phase_basis, const struct gkyl_basis *conf_basis, bool use_gpu);

/**
 * Compute the cross moments with moments of each species.
 *
 * @param up Cross moments updater.
 * @param conf_rng Configuration space range.
 * @param delta_sr Free parameter in relation between nu_sr and Morse's alpha_E.
 * @param betap1 Greene's beta parameter plus 1.
 * @param m_self Mass of the self species.
 * @param moms_self Moments of the self species.
 * @param m_other Mass of the other species.
 * @param moms_other Moments of the other species.
 * @param moms_cross Six output moments.
 */
void gkyl_gyrokinetic_cross_prim_moms_bgk_advance(gkyl_gyrokinetic_cross_prim_moms_bgk *up,
  const struct gkyl_range *conf_rng, double delta_sr, double betap1, double m_self,
  const struct gkyl_array *prim_moms_self, double m_other, const struct gkyl_array *prim_moms_other,
  struct gkyl_array *prim_moms_cross);

/**
 * Release memory needed in the cross moments calculation.
 *
 * @param up Memory to release.
 */
void gkyl_gyrokinetic_cross_prim_moms_bgk_release(gkyl_gyrokinetic_cross_prim_moms_bgk *up);
