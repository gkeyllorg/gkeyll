// Private header: not for direct use
#pragma once

#include <gkyl_prim_lbo_type.h>
#include <gkyl_prim_lbo_gyrokinetic_kernels.h>
#include <gkyl_mat.h>
#include <gkyl_util.h>

typedef void (*gyrokinetic_self_prim_t)(struct gkyl_mat *A, struct gkyl_mat *rhs, 
  const double *moms, const double *boundary_corrections, const double *nu);

typedef void (*gyrokinetic_cross_prim_t)(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *alpha_E,
  const double m_self, const double *moms_self, const double *prim_moms_self,
  const double m_other, const double *moms_other, const double *prim_moms_other,
  const double *boundary_corrections, const double *nu);

// for use in kernel tables
typedef struct { gyrokinetic_self_prim_t kernels[3]; } gkyl_prim_lbo_gyrokinetic_kern_list;
typedef struct { gyrokinetic_cross_prim_t kernels[3]; } gkyl_prim_lbo_gyrokinetic_cross_kern_list;

//
// Serendipity basis kernels
//

// self primitive moment kernel list
GKYL_CU_D
static const gkyl_prim_lbo_gyrokinetic_kern_list ser_self_prim_kernels[] = {
  // 1x kernels
  { NULL, gyrokinetic_self_prim_moments_1x1v_ser_p1, gyrokinetic_self_prim_moments_1x1v_ser_p2 }, // 0
  { NULL, gyrokinetic_self_prim_moments_1x2v_ser_p1, gyrokinetic_self_prim_moments_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, gyrokinetic_self_prim_moments_2x2v_ser_p1, gyrokinetic_self_prim_moments_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, gyrokinetic_self_prim_moments_3x2v_ser_p1, NULL }, // 3
};

// cross primitive moment kernel list
GKYL_CU_D
static const gkyl_prim_lbo_gyrokinetic_cross_kern_list ser_cross_prim_kernels[] = {
  // 1x kernels
  { NULL, gyrokinetic_cross_prim_moments_1x1v_ser_p1, gyrokinetic_cross_prim_moments_1x1v_ser_p2 }, // 0
  { NULL, gyrokinetic_cross_prim_moments_1x2v_ser_p1, gyrokinetic_cross_prim_moments_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, gyrokinetic_cross_prim_moments_2x2v_ser_p1, gyrokinetic_cross_prim_moments_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, gyrokinetic_cross_prim_moments_3x2v_ser_p1, NULL }, // 3
};

struct prim_lbo_type_gyrokinetic {
  struct gkyl_prim_lbo_type prim; // Base object
  gyrokinetic_self_prim_t self_prim; // Self-primitive moments kernel
  gyrokinetic_cross_prim_t cross_prim; // Cross-primitive moments kernels
};

/**
 * Free primitive moment object.
 *
 * @param ref Reference counter for primitive moment to free
 */
void prim_lbo_gyrokinetic_free(const struct gkyl_ref_count *ref);

GKYL_CU_D
static void
self_prim(const struct gkyl_prim_lbo_type *prim, struct gkyl_mat *A, struct gkyl_mat *rhs, 
  const int* idx, const double *moms, const double *boundary_corrections, const double *nu)
{
  struct prim_lbo_type_gyrokinetic *prim_gyrokinetic = container_of(prim, struct prim_lbo_type_gyrokinetic, prim);

  return prim_gyrokinetic->self_prim(A, rhs, moms, boundary_corrections, nu);
}

GKYL_CU_D
static void
cross_prim(const struct gkyl_prim_lbo_type *prim, struct gkyl_mat *A, struct gkyl_mat *rhs,
  const int *idx, const double *alpha_E,
  const double m_self, const double *moms_self, const double *prim_moms_self,
  const double m_other, const double *moms_other, const double *prim_moms_other,
  const double *boundary_corrections, const double *nu)
{
  struct prim_lbo_type_gyrokinetic *prim_gyrokinetic = container_of(prim, struct prim_lbo_type_gyrokinetic, prim);

  return prim_gyrokinetic->cross_prim(A, rhs, alpha_E, m_self, moms_self, prim_moms_self,
    m_other, moms_other, prim_moms_other, boundary_corrections, nu);
}

#ifdef GKYL_HAVE_CUDA

struct gkyl_prim_lbo_type*
gkyl_prim_lbo_gyrokinetic_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis);

#endif
