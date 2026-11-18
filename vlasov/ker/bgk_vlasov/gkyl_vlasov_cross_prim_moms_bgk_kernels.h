#pragma once 
#include <math.h> 
#include <gkyl_mat.h> 
#include <gkyl_util.h> 
 
EXTERN_C_BEG 

GKYL_CU_DH void vlasov_cross_prim_moms_bgk_1x1v_ser_p1(const double delta_sr, const double betaGreenep1, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, double *prim_moms_cross); 

GKYL_CU_DH void vlasov_cross_prim_moms_bgk_1x1v_ser_p2(const double delta_sr, const double betaGreenep1, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, double *prim_moms_cross); 

GKYL_CU_DH void vlasov_cross_prim_moms_bgk_1x1v_ser_p3(const double delta_sr, const double betaGreenep1, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, double *prim_moms_cross); 

GKYL_CU_DH void vlasov_cross_prim_moms_bgk_1x2v_ser_p1(const double delta_sr, const double betaGreenep1, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, double *prim_moms_cross); 

GKYL_CU_DH void vlasov_cross_prim_moms_bgk_1x2v_ser_p2(const double delta_sr, const double betaGreenep1, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, double *prim_moms_cross); 

GKYL_CU_DH void vlasov_cross_prim_moms_bgk_1x3v_ser_p1(const double delta_sr, const double betaGreenep1, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, double *prim_moms_cross); 

GKYL_CU_DH void vlasov_cross_prim_moms_bgk_1x3v_ser_p2(const double delta_sr, const double betaGreenep1, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, double *prim_moms_cross); 

GKYL_CU_DH void vlasov_cross_prim_moms_bgk_2x1v_ser_p1(const double delta_sr, const double betaGreenep1, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, double *prim_moms_cross); 

GKYL_CU_DH void vlasov_cross_prim_moms_bgk_2x1v_ser_p2(const double delta_sr, const double betaGreenep1, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, double *prim_moms_cross); 

GKYL_CU_DH void vlasov_cross_prim_moms_bgk_2x1v_ser_p3(const double delta_sr, const double betaGreenep1, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, double *prim_moms_cross); 

GKYL_CU_DH void vlasov_cross_prim_moms_bgk_2x2v_ser_p1(const double delta_sr, const double betaGreenep1, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, double *prim_moms_cross); 

GKYL_CU_DH void vlasov_cross_prim_moms_bgk_2x2v_ser_p2(const double delta_sr, const double betaGreenep1, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, double *prim_moms_cross); 

GKYL_CU_DH void vlasov_cross_prim_moms_bgk_2x3v_ser_p1(const double delta_sr, const double betaGreenep1, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, double *prim_moms_cross); 

GKYL_CU_DH void vlasov_cross_prim_moms_bgk_2x3v_ser_p2(const double delta_sr, const double betaGreenep1, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, double *prim_moms_cross); 

GKYL_CU_DH void vlasov_cross_prim_moms_bgk_3x3v_ser_p1(const double delta_sr, const double betaGreenep1, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, double *prim_moms_cross); 

GKYL_CU_DH void vlasov_cross_prim_moms_bgk_1x1v_tensor_p2(const double delta_sr, const double betaGreenep1, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, double *prim_moms_cross); 

GKYL_CU_DH void vlasov_cross_prim_moms_bgk_1x1v_tensor_p3(const double delta_sr, const double betaGreenep1, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, double *prim_moms_cross); 

GKYL_CU_DH void vlasov_cross_prim_moms_bgk_1x2v_tensor_p2(const double delta_sr, const double betaGreenep1, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, double *prim_moms_cross); 

GKYL_CU_DH void vlasov_cross_prim_moms_bgk_1x3v_tensor_p2(const double delta_sr, const double betaGreenep1, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, double *prim_moms_cross); 

GKYL_CU_DH void vlasov_cross_prim_moms_bgk_2x1v_tensor_p2(const double delta_sr, const double betaGreenep1, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, double *prim_moms_cross); 

GKYL_CU_DH void vlasov_cross_prim_moms_bgk_2x1v_tensor_p3(const double delta_sr, const double betaGreenep1, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, double *prim_moms_cross); 

GKYL_CU_DH void vlasov_cross_prim_moms_bgk_2x2v_tensor_p2(const double delta_sr, const double betaGreenep1, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, double *prim_moms_cross); 

GKYL_CU_DH void vlasov_cross_prim_moms_bgk_2x3v_tensor_p2(const double delta_sr, const double betaGreenep1, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, double *prim_moms_cross); 

EXTERN_C_END 
