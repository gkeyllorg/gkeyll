#pragma once 
#include <math.h> 
#include <gkyl_mat.h> 
#include <gkyl_util.h> 
 
EXTERN_C_BEG 

GKYL_CU_DH void vlasov_self_prim_moments_1x1v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *moms, const double *boundary_corrections, const double *nu); 
GKYL_CU_DH void vlasov_cross_prim_moments_1x1v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *greene, const double m_self, const double *moms_self, const double *prim_mom_self, const double m_other, const double *moms_other, const double *prim_mom_other, const double *boundary_corrections, const double *nu); 

GKYL_CU_DH void vlasov_self_prim_moments_1x1v_ser_p3(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *moms, const double *boundary_corrections, const double *nu); 
GKYL_CU_DH void vlasov_cross_prim_moments_1x1v_ser_p3(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *greene, const double m_self, const double *moms_self, const double *prim_mom_self, const double m_other, const double *moms_other, const double *prim_mom_other, const double *boundary_corrections, const double *nu); 

GKYL_CU_DH void vlasov_self_prim_moments_1x2v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *moms, const double *boundary_corrections, const double *nu); 
GKYL_CU_DH void vlasov_cross_prim_moments_1x2v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *greene, const double m_self, const double *moms_self, const double *prim_mom_self, const double m_other, const double *moms_other, const double *prim_mom_other, const double *boundary_corrections, const double *nu); 

GKYL_CU_DH void vlasov_self_prim_moments_1x3v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *moms, const double *boundary_corrections, const double *nu); 
GKYL_CU_DH void vlasov_cross_prim_moments_1x3v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *greene, const double m_self, const double *moms_self, const double *prim_mom_self, const double m_other, const double *moms_other, const double *prim_mom_other, const double *boundary_corrections, const double *nu); 

GKYL_CU_DH void vlasov_self_prim_moments_2x1v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *moms, const double *boundary_corrections, const double *nu); 
GKYL_CU_DH void vlasov_cross_prim_moments_2x1v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *greene, const double m_self, const double *moms_self, const double *prim_mom_self, const double m_other, const double *moms_other, const double *prim_mom_other, const double *boundary_corrections, const double *nu); 

GKYL_CU_DH void vlasov_self_prim_moments_2x1v_ser_p3(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *moms, const double *boundary_corrections, const double *nu); 
GKYL_CU_DH void vlasov_cross_prim_moments_2x1v_ser_p3(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *greene, const double m_self, const double *moms_self, const double *prim_mom_self, const double m_other, const double *moms_other, const double *prim_mom_other, const double *boundary_corrections, const double *nu); 

GKYL_CU_DH void vlasov_self_prim_moments_2x2v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *moms, const double *boundary_corrections, const double *nu); 
GKYL_CU_DH void vlasov_cross_prim_moments_2x2v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *greene, const double m_self, const double *moms_self, const double *prim_mom_self, const double m_other, const double *moms_other, const double *prim_mom_other, const double *boundary_corrections, const double *nu); 

GKYL_CU_DH void vlasov_self_prim_moments_2x3v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *moms, const double *boundary_corrections, const double *nu); 
GKYL_CU_DH void vlasov_cross_prim_moments_2x3v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *greene, const double m_self, const double *moms_self, const double *prim_mom_self, const double m_other, const double *moms_other, const double *prim_mom_other, const double *boundary_corrections, const double *nu); 

GKYL_CU_DH void vlasov_self_prim_moments_1x1v_tensor_p1(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *moms, const double *boundary_corrections, const double *nu); 
GKYL_CU_DH void vlasov_cross_prim_moments_1x1v_tensor_p1(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *greene, const double m_self, const double *moms_self, const double *prim_mom_self, const double m_other, const double *moms_other, const double *prim_mom_other, const double *boundary_corrections, const double *nu); 

GKYL_CU_DH void vlasov_self_prim_moments_1x1v_tensor_p2(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *moms, const double *boundary_corrections, const double *nu); 
GKYL_CU_DH void vlasov_cross_prim_moments_1x1v_tensor_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *greene, const double m_self, const double *moms_self, const double *prim_mom_self, const double m_other, const double *moms_other, const double *prim_mom_other, const double *boundary_corrections, const double *nu); 

GKYL_CU_DH void vlasov_self_prim_moments_1x1v_tensor_p3(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *moms, const double *boundary_corrections, const double *nu); 
GKYL_CU_DH void vlasov_cross_prim_moments_1x1v_tensor_p3(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *greene, const double m_self, const double *moms_self, const double *prim_mom_self, const double m_other, const double *moms_other, const double *prim_mom_other, const double *boundary_corrections, const double *nu); 

GKYL_CU_DH void vlasov_self_prim_moments_1x2v_tensor_p1(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *moms, const double *boundary_corrections, const double *nu); 
GKYL_CU_DH void vlasov_cross_prim_moments_1x2v_tensor_p1(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *greene, const double m_self, const double *moms_self, const double *prim_mom_self, const double m_other, const double *moms_other, const double *prim_mom_other, const double *boundary_corrections, const double *nu); 

GKYL_CU_DH void vlasov_self_prim_moments_1x2v_tensor_p2(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *moms, const double *boundary_corrections, const double *nu); 
GKYL_CU_DH void vlasov_cross_prim_moments_1x2v_tensor_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *greene, const double m_self, const double *moms_self, const double *prim_mom_self, const double m_other, const double *moms_other, const double *prim_mom_other, const double *boundary_corrections, const double *nu); 

GKYL_CU_DH void vlasov_self_prim_moments_1x3v_tensor_p1(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *moms, const double *boundary_corrections, const double *nu); 
GKYL_CU_DH void vlasov_cross_prim_moments_1x3v_tensor_p1(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *greene, const double m_self, const double *moms_self, const double *prim_mom_self, const double m_other, const double *moms_other, const double *prim_mom_other, const double *boundary_corrections, const double *nu); 

GKYL_CU_DH void vlasov_self_prim_moments_1x3v_tensor_p2(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *moms, const double *boundary_corrections, const double *nu); 
GKYL_CU_DH void vlasov_cross_prim_moments_1x3v_tensor_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *greene, const double m_self, const double *moms_self, const double *prim_mom_self, const double m_other, const double *moms_other, const double *prim_mom_other, const double *boundary_corrections, const double *nu); 

GKYL_CU_DH void vlasov_self_prim_moments_2x1v_tensor_p1(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *moms, const double *boundary_corrections, const double *nu); 
GKYL_CU_DH void vlasov_cross_prim_moments_2x1v_tensor_p1(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *greene, const double m_self, const double *moms_self, const double *prim_mom_self, const double m_other, const double *moms_other, const double *prim_mom_other, const double *boundary_corrections, const double *nu); 

GKYL_CU_DH void vlasov_self_prim_moments_2x1v_tensor_p2(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *moms, const double *boundary_corrections, const double *nu); 
GKYL_CU_DH void vlasov_cross_prim_moments_2x1v_tensor_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *greene, const double m_self, const double *moms_self, const double *prim_mom_self, const double m_other, const double *moms_other, const double *prim_mom_other, const double *boundary_corrections, const double *nu); 

GKYL_CU_DH void vlasov_self_prim_moments_2x1v_tensor_p3(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *moms, const double *boundary_corrections, const double *nu); 
GKYL_CU_DH void vlasov_cross_prim_moments_2x1v_tensor_p3(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *greene, const double m_self, const double *moms_self, const double *prim_mom_self, const double m_other, const double *moms_other, const double *prim_mom_other, const double *boundary_corrections, const double *nu); 

GKYL_CU_DH void vlasov_self_prim_moments_2x2v_tensor_p1(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *moms, const double *boundary_corrections, const double *nu); 
GKYL_CU_DH void vlasov_cross_prim_moments_2x2v_tensor_p1(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *greene, const double m_self, const double *moms_self, const double *prim_mom_self, const double m_other, const double *moms_other, const double *prim_mom_other, const double *boundary_corrections, const double *nu); 

GKYL_CU_DH void vlasov_self_prim_moments_2x2v_tensor_p2(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *moms, const double *boundary_corrections, const double *nu); 
GKYL_CU_DH void vlasov_cross_prim_moments_2x2v_tensor_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *greene, const double m_self, const double *moms_self, const double *prim_mom_self, const double m_other, const double *moms_other, const double *prim_mom_other, const double *boundary_corrections, const double *nu); 

GKYL_CU_DH void vlasov_self_prim_moments_2x3v_tensor_p1(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *moms, const double *boundary_corrections, const double *nu); 
GKYL_CU_DH void vlasov_cross_prim_moments_2x3v_tensor_p1(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *greene, const double m_self, const double *moms_self, const double *prim_mom_self, const double m_other, const double *moms_other, const double *prim_mom_other, const double *boundary_corrections, const double *nu); 

GKYL_CU_DH void vlasov_self_prim_moments_2x3v_tensor_p2(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *moms, const double *boundary_corrections, const double *nu); 
GKYL_CU_DH void vlasov_cross_prim_moments_2x3v_tensor_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *greene, const double m_self, const double *moms_self, const double *prim_mom_self, const double m_other, const double *moms_other, const double *prim_mom_other, const double *boundary_corrections, const double *nu); 

GKYL_CU_DH void vlasov_self_prim_moments_3x3v_tensor_p1(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *moms, const double *boundary_corrections, const double *nu); 
GKYL_CU_DH void vlasov_cross_prim_moments_3x3v_tensor_p1(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *greene, const double m_self, const double *moms_self, const double *prim_mom_self, const double m_other, const double *moms_other, const double *prim_mom_other, const double *boundary_corrections, const double *nu); 

EXTERN_C_END 
