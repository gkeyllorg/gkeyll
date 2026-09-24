// Private header: not for direct use
#pragma once

#include <gkyl_array.h>
#include <gkyl_prim_lbo_type.h>
#include <gkyl_prim_lbo_pkpm_kernels.h>
#include <gkyl_mat.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

typedef void (*pkpm_self_prim_t)(struct gkyl_mat *A, struct gkyl_mat *rhs, 
  const double *moms, const double *boundary_corrections, const double *nu);

// for use in kernel tables
typedef struct { pkpm_self_prim_t kernels[3]; } gkyl_prim_lbo_pkpm_self_kern_list;

// PKPM self-primitive moment kernel list (Serendipity basis)
GKYL_CU_D
static const gkyl_prim_lbo_pkpm_self_kern_list ser_self_prim_kernels[] = {
  // 1x kernels
  { NULL, pkpm_self_prim_moments_1x1v_ser_p1, pkpm_self_prim_moments_1x1v_ser_p2 }, // 0
  // 2x kernels
  { NULL, pkpm_self_prim_moments_2x1v_ser_p1, NULL }, // 1
  // 3x kernels
  { NULL, pkpm_self_prim_moments_3x1v_ser_p1, NULL }, // 2
};

// PKPM self-primitive moment kernel list (Tensor basis)
GKYL_CU_D
static const gkyl_prim_lbo_pkpm_self_kern_list ten_self_prim_kernels[] = {
  // 1x kernels
  { NULL, pkpm_self_prim_moments_1x1v_ser_p1, pkpm_self_prim_moments_1x1v_tensor_p2 }, // 0
  // 2x kernels
  { NULL, pkpm_self_prim_moments_2x1v_ser_p1, pkpm_self_prim_moments_2x1v_tensor_p2 }, // 1
  // 3x kernels
  { NULL, pkpm_self_prim_moments_3x1v_ser_p1, NULL }, // 2
};

struct prim_lbo_type_pkpm {
  struct gkyl_prim_lbo_type prim; // Base object
  pkpm_self_prim_t self_prim; // Self-primitive moments kernel
  struct gkyl_range conf_range; // configuration space range
};

/**
 * Free primitive moment object.
 *
 * @param ref Reference counter for primitive moment to free
 */
void prim_lbo_pkpm_free(const struct gkyl_ref_count *ref);

GKYL_CU_D
static void
self_prim(const struct gkyl_prim_lbo_type *prim, struct gkyl_mat *A, struct gkyl_mat *rhs, 
  const int* idx, const double *moms, const double *boundary_corrections, const double *nu)
{
  struct prim_lbo_type_pkpm *prim_pkpm = container_of(prim, struct prim_lbo_type_pkpm, prim);

  long cidx = gkyl_range_idx(&prim_pkpm->conf_range, idx);
  return prim_pkpm->self_prim(A, rhs, moms, boundary_corrections, nu);
}
