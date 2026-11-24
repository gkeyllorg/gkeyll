#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>

// Updater for gyrokinetic BGK cross-species primitive moments.
typedef struct gkyl_gyrokinetic_cross_prim_moms_bgk gkyl_gyrokinetic_cross_prim_moms_bgk;

/**
 * Create an updater for gyrokinetic BGK cross-species primitive moments.
 *
 * @param phase_basis Phase-space basis functions.
 * @param conf_basis Configuration-space basis functions.
 * @param use_gpu Whether to perform calculations on the GPU.
 * @return New cross-species primitive moment updater.
 */
gkyl_gyrokinetic_cross_prim_moms_bgk *gkyl_gyrokinetic_cross_prim_moms_bgk_new(
  const struct gkyl_basis *phase_basis, const struct gkyl_basis *conf_basis, bool use_gpu
);

/**
 * Compute the BGK cross-species primitive moments from each species' primitive moments.
 *
 * Each moment array contains three configuration-space DG fields, ordered as
 * number density, parallel flow velocity, and squared thermal speed (T/m).
 * The output density equals the self-species density.
 *
 * @param up Cross moments updater.
 * @param conf_rng Configuration-space range over which to compute the moments.
 * @param delta_sr Free parameter in relation between nu_sr and Morse's alpha_E.
 * @param betap1 Greene's beta parameter plus 1.
 * @param m_self Mass of the self species.
 * @param prim_moms_self Input primitive moments of the self species.
 * @param m_other Mass of the other species.
 * @param prim_moms_other Input primitive moments of the other species.
 * @param prim_moms_cross Output cross-species primitive moments for the self species.
 */
void gkyl_gyrokinetic_cross_prim_moms_bgk_advance(
  gkyl_gyrokinetic_cross_prim_moms_bgk *up, const struct gkyl_range *conf_rng, double delta_sr,
  double betap1, double m_self, const struct gkyl_array *prim_moms_self, double m_other,
  const struct gkyl_array *prim_moms_other, struct gkyl_array *prim_moms_cross
);

/**
 * Release the cross-species primitive moment updater and its resources.
 *
 * @param up Updater to release.
 */
void gkyl_gyrokinetic_cross_prim_moms_bgk_release(gkyl_gyrokinetic_cross_prim_moms_bgk *up);

/**
 * Create a gyrokinetic BGK cross-species primitive moment updater without a
 * phase-space basis, allowing fluid collision partners.
 *
 * @param conf_basis Configuration-space basis functions.
 * @param vdim Number of gyrokinetic velocity coordinates (1 for v_parallel,
 *   2 for v_parallel and mu).
 * @param use_gpu Whether to perform calculations on the GPU.
 * @return New cross-species primitive moment updater.
 */
gkyl_gyrokinetic_cross_prim_moms_bgk *gkyl_gyrokinetic_cross_prim_moms_bgk_new_with_vdim(
  const struct gkyl_basis *conf_basis, int vdim, bool use_gpu
);
