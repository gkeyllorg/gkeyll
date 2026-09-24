#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>

// Object type
typedef struct gkyl_prim_cross_m0deltas gkyl_prim_cross_m0deltas;

/**
 * Create a new updater that computes the m0_s*delta_s*(beta+1) prefactor
 * in the calculation of cross-primitive moments for LBO and BGK
 * collisions. That is if the collision frequency is constant:
 *         n_s*delta_s*(beta+1) =         2*(beta+1) * n_s * m_r * n_r * nu_rs / (m_s * n_s * nu_sr + m_r * n_r * nu_rs)
 * and if the collision frequency varies in space and time:
 *   nu_sr*n_s*delta_s*(beta+1) = nu_sr * 2*(beta+1) * n_s * m_r * n_r * nu_rs / (m_s * n_s * nu_sr + m_r * n_r * nu_rs)
 *
 * @param normNu Whether we are using nu(x,t).
 * @param basis Basis object (configuration space).
 * @param range Range in which we'll compute m0_s*delta_s.
 * @param betap1 beta+1 parameter in Greene's formulism.
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_prim_cross_m0deltas* gkyl_prim_cross_m0deltas_new(bool normNu,
  const struct gkyl_basis *basis, const struct gkyl_range *range,
  double betap1, bool use_gpu);

/**
 * Compute 
 *         n_s*delta_s*(beta+1) =         2*(beta+1) * n_s * m_r * n_r * nu_rs / (m_s * n_s * nu_sr + m_r * n_r * nu_rs)
 * if collision frequency is constant, or
 *   nu_sr*n_s*delta_s*(beta+1) = nu_sr * 2*(beta+1) * n_s * m_r * n_r * nu_rs / (m_s * n_s * nu_sr + m_r * n_r * nu_rs)
 * if the collision frequency varies in space and time.
 *
 * @param up Struct defining this updater..
 * @param massself Mass of this species, m_s.
 * @param m0self Number density of this species, m0_s.
 * @param nuself Cross collision frequency of this species, nu_sr.
 * @param massother Mass of the other species, m_r.
 * @param m0other Number density of the other species, m0_r.
 * @param nuother Cross collision frequency of the other species, nu_rs.
 * @param out Output array.
 * @return New updater pointer.
 */
void gkyl_prim_cross_m0deltas_advance(gkyl_prim_cross_m0deltas *up,
  double massself, const struct gkyl_array* m0self, const struct gkyl_array* nuself,
  double massother, const struct gkyl_array* m0other, const struct gkyl_array* nuother,
  struct gkyl_array* out);

/**
 * Delete updater.
 *
 * @param pob Updater to delete.
 */
void gkyl_prim_cross_m0deltas_release(gkyl_prim_cross_m0deltas* up);
