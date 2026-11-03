#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH void vlasov_sr_M0_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M1i_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M2_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M3i_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_Ni_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_Tij_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_int_five_moments_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_sr_M0_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M1i_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M2_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M3i_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_Ni_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_Tij_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_int_five_moments_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_sr_M0_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M1i_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M2_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M3i_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_Ni_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_Tij_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_int_five_moments_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_sr_M0_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M1i_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M2_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M3i_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_Ni_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_Tij_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_int_five_moments_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_sr_M0_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M1i_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M2_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M3i_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_Ni_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_Tij_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_int_five_moments_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_sr_M0_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M1i_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M2_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M3i_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_Ni_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_Tij_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_int_five_moments_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_sr_M0_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M1i_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M2_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M3i_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_Ni_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_Tij_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_int_five_moments_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_sr_M0_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M1i_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M2_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M3i_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_Ni_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_Tij_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_int_five_moments_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_sr_M0_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M1i_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M2_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M3i_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_Ni_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_Tij_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_int_five_moments_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_sr_M0_2x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M1i_2x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M2_2x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M3i_2x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_Ni_2x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_Tij_2x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_int_five_moments_2x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out); 


EXTERN_C_END 
