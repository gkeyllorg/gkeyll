#pragma once

#include <gkyl_basis.h>
#include <gkyl_prim_lbo_type.h>

/**
 * Create a gyrokinetic Lenard-Bernstein collision primitive moment object
 * for self- and cross-species collisions.
 *
 * @param cbasis Configuration-space basis functions.
 * @param pbasis Phase-space basis functions.
 * @param use_gpu Whether to perform calculations on the GPU.
 * @return New primitive moment object; release with gkyl_prim_lbo_type_release.
 */
struct gkyl_prim_lbo_type *gkyl_prim_lbo_gyrokinetic_new(
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis, bool use_gpu
);

/**
 * Create a gyrokinetic Lenard-Bernstein collision primitive moment object
 * without a phase-space basis, allowing fluid collision partners.
 *
 * @param cbasis Configuration-space basis functions.
 * @param vdim Number of gyrokinetic velocity coordinates (1 for v_parallel,
 *   2 for v_parallel and mu).
 * @param num_phase Number of phase-space basis functions; use 0 for fluid species.
 * @param use_gpu Whether to perform calculations on the GPU.
 * @return New primitive moment object; release with gkyl_prim_lbo_type_release.
 */
struct gkyl_prim_lbo_type *gkyl_prim_lbo_gyrokinetic_new_with_vdim(
  const struct gkyl_basis *cbasis, int vdim, int num_phase, bool use_gpu
);
