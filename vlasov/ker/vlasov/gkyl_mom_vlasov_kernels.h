#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH void vlasov_M0_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M1i_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2ij_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3i_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3ijk_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_five_moments_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_int_five_moments_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 


GKYL_CU_DH void vlasov_M0_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M1i_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2ij_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3i_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3ijk_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_five_moments_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_int_five_moments_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 


GKYL_CU_DH void vlasov_M0_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M1i_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2ij_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3i_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3ijk_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_five_moments_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_int_five_moments_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 


GKYL_CU_DH void vlasov_M0_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M1i_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2ij_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3i_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3ijk_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_five_moments_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_int_five_moments_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 


GKYL_CU_DH void vlasov_M0_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M1i_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2ij_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3i_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3ijk_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_five_moments_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_int_five_moments_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 


GKYL_CU_DH void vlasov_M0_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M1i_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2ij_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3i_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3ijk_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_five_moments_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_int_five_moments_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 


GKYL_CU_DH void vlasov_M0_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M1i_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2ij_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3i_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3ijk_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_five_moments_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_int_five_moments_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 


GKYL_CU_DH void vlasov_M0_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M1i_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2ij_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3i_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3ijk_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_five_moments_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_int_five_moments_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 


GKYL_CU_DH void vlasov_M0_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M1i_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2ij_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3i_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3ijk_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_five_moments_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_int_five_moments_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 






GKYL_CU_DH void vlasov_M0_1x1v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M1i_1x1v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2_1x1v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2ij_1x1v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3i_1x1v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3ijk_1x1v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_five_moments_1x1v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_int_five_moments_1x1v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 


GKYL_CU_DH void vlasov_M0_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M1i_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2ij_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3i_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3ijk_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_five_moments_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_int_five_moments_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 


GKYL_CU_DH void vlasov_M0_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M1i_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2ij_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3i_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3ijk_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_five_moments_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_int_five_moments_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 


GKYL_CU_DH void vlasov_M0_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M1i_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2ij_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3i_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3ijk_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_five_moments_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_int_five_moments_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 


GKYL_CU_DH void vlasov_M0_1x3v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M1i_1x3v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2_1x3v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2ij_1x3v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3i_1x3v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3ijk_1x3v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_five_moments_1x3v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_int_five_moments_1x3v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 


GKYL_CU_DH void vlasov_M0_1x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M1i_1x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2_1x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2ij_1x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3i_1x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3ijk_1x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_five_moments_1x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_int_five_moments_1x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 


GKYL_CU_DH void vlasov_M0_2x2v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M1i_2x2v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2_2x2v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2ij_2x2v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3i_2x2v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3ijk_2x2v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_five_moments_2x2v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_int_five_moments_2x2v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 


GKYL_CU_DH void vlasov_M0_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M1i_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2ij_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3i_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3ijk_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_five_moments_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_int_five_moments_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 


GKYL_CU_DH void vlasov_M0_2x3v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M1i_2x3v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2_2x3v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2ij_2x3v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3i_2x3v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3ijk_2x3v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_five_moments_2x3v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_int_five_moments_2x3v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 






EXTERN_C_END 
