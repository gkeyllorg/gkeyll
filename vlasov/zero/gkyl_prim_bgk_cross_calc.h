#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>

/**
 * Calculate the cross primitive moments, u_{sr,i} and v_{tsr}^2,
 * for cross-species collisions with the BGK operator.
 *
 *   u_sr = u_si-0.5*delta_s*(beta+1)*(u_si - u_ri)
 *
 *   v_tsr^2 = v_ts^2-(u_sri-u_ri).(u_sri-u_si)/vdim_phys
 *     -(delta_s*(beta+1)/(m_s+m_r))*
 *      (m_s*v_ts^2-m_r*v_tr^2+((m_s-m_r)/(2*vdim_phys))*(u_si-u_ri)^2)
 *
 * @param basis Basis object (configuration space).
 * @param vdim_phys Physical number of velocity dimensions represented.
 * @param m0sdeltas Prefactor m_0s*delta_s*(beta+1) (and m_0s should be 1 here).
 * @param massself Mass of this species.
 * @param primsself Primitive moments of this species, u_s and v_ts^2.
 * @param massother Mass of the other species.
 * @param primsother Primitive moments of the other species, u_r and v_tr^2.
 * @param range Range in which we'll compute m0_s*delta_s.
 * @return crossprims Cross primitive moments, u_sri and v_tsr^2.
 */
void gkyl_prim_bgk_cross_calc_advance(struct gkyl_basis basis,
  int vdim_phys, const struct gkyl_array* m0sdeltas,
  double massself, const struct gkyl_array* primsself,
  double massother, const struct gkyl_array* primsother,
  const struct gkyl_range *range, struct gkyl_array* crossprims);

void gkyl_prim_bgk_cross_calc_advance_cu(struct gkyl_basis basis,
  int vdim_phys, const struct gkyl_array* m0sdeltas,
  double massself, const struct gkyl_array* primsself,
  double massother, const struct gkyl_array* primsother,
  const struct gkyl_range *range, struct gkyl_array* crossprims);
