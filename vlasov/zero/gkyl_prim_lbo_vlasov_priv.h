// Private header: not for direct use
#pragma once

#include <gkyl_prim_lbo_type.h>
#include <gkyl_prim_lbo_vlasov_kernels.h>
#include <gkyl_mat.h>
#include <gkyl_util.h>

typedef void (*vlasov_self_prim_t)(struct gkyl_mat *A, struct gkyl_mat *rhs, 
  const double *moms, const double *boundary_corrections, const double *nu);

typedef void (*vlasov_cross_prim_t)(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *greene,
  const double m_self, const double *moms_self, const double *prim_moms_self,
  const double m_other, const double *moms_other, const double *prim_moms_other,
  const double *boundary_corrections, const double *nu);

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
static struct { int vdim[4]; } cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1,  0,  1,  2}, // 1x kernel indices
  {-1, -1,  3,  4}, // 2x kernel indices
  {-1, -1, -1,  5}, // 3x kernel indices  
};

// for use in kernel tables
typedef struct { vlasov_self_prim_t kernels[3]; } gkyl_prim_lbo_vlasov_self_kern_list;
typedef struct { vlasov_cross_prim_t kernels[3]; } gkyl_prim_lbo_vlasov_cross_kern_list;


//
// Serendipity basis kernels
//

// self-primitive moment kernel list
GKYL_CU_D
static const gkyl_prim_lbo_vlasov_self_kern_list ser_self_prim_kernels[] = {
  // 1x kernels
  { NULL, vlasov_self_prim_moments_1x1v_ser_p1, vlasov_self_prim_moments_1x1v_ser_p2 }, // 0
  { NULL, vlasov_self_prim_moments_1x2v_ser_p1, vlasov_self_prim_moments_1x2v_ser_p2 }, // 1
  { NULL, vlasov_self_prim_moments_1x3v_ser_p1, vlasov_self_prim_moments_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_self_prim_moments_2x2v_ser_p1, vlasov_self_prim_moments_2x2v_ser_p2 }, // 3
  { NULL, vlasov_self_prim_moments_2x3v_ser_p1, vlasov_self_prim_moments_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_self_prim_moments_3x3v_ser_p1, NULL }, // 5
};

// cross-primitive moment kernel list
GKYL_CU_D
static const gkyl_prim_lbo_vlasov_cross_kern_list ser_cross_prim_kernels[] = {
  // 1x kernels
  { NULL, vlasov_cross_prim_moments_1x1v_ser_p1, vlasov_cross_prim_moments_1x1v_ser_p2 }, // 0
  { NULL, vlasov_cross_prim_moments_1x2v_ser_p1, vlasov_cross_prim_moments_1x2v_ser_p2 }, // 1
  { NULL, vlasov_cross_prim_moments_1x3v_ser_p1, vlasov_cross_prim_moments_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_cross_prim_moments_2x2v_ser_p1, vlasov_cross_prim_moments_2x2v_ser_p2 }, // 3
  { NULL, vlasov_cross_prim_moments_2x3v_ser_p1, vlasov_cross_prim_moments_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_cross_prim_moments_3x3v_ser_p1, NULL }, // 5
};

struct prim_lbo_type_vlasov {
  struct gkyl_prim_lbo_type prim; // Base object
  vlasov_self_prim_t self_prim; // Self-primitive moments kernel  
  vlasov_cross_prim_t cross_prim; // Cross-primitive moments kernels
};

/**
 * Free primitive moment object.
 *
 * @param ref Reference counter for primitive moment to free
 */
void prim_lbo_vlasov_free(const struct gkyl_ref_count *ref);

GKYL_CU_D
static void
self_prim(const struct gkyl_prim_lbo_type *prim, struct gkyl_mat *A, struct gkyl_mat *rhs, 
  const int *idx, const double *moms, const double *boundary_corrections, const double *nu)
{
  struct prim_lbo_type_vlasov *prim_vlasov = container_of(prim, struct prim_lbo_type_vlasov, prim);

  return prim_vlasov->self_prim(A, rhs, moms, boundary_corrections, nu);
}

GKYL_CU_D
static void
cross_prim(const struct gkyl_prim_lbo_type *prim, struct gkyl_mat *A, struct gkyl_mat *rhs, 
  const int *idx, const double *greene,
  const double m_self, const double *moms_self, const double *prim_moms_self,
  const double m_other, const double *moms_other, const double *prim_moms_other,
  const double *boundary_corrections, const double *nu)
{
  struct prim_lbo_type_vlasov *prim_vlasov = container_of(prim, struct prim_lbo_type_vlasov, prim);

  return prim_vlasov->cross_prim(A, rhs, greene, m_self, moms_self, prim_moms_self,
    m_other, moms_other, prim_moms_other, boundary_corrections, nu);
}
