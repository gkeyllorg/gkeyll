#pragma once
#include <math.h>
#include <gkyl_eqn_type.h>
#include <gkyl_util.h>

EXTERN_C_BEG 

GKYL_CU_DH void mom_bcorr_lbo_vlasov_1x1v_ser_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_bcorr_lbo_vlasov_1x1v_ser_p3(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_bcorr_lbo_vlasov_1x2v_ser_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_bcorr_lbo_vlasov_1x3v_ser_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_bcorr_lbo_vlasov_2x1v_ser_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_bcorr_lbo_vlasov_2x1v_ser_p3(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_bcorr_lbo_vlasov_2x2v_ser_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_bcorr_lbo_vlasov_2x3v_ser_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_bcorr_lbo_vlasov_1x1v_tensor_p1(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_bcorr_lbo_vlasov_1x1v_tensor_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_bcorr_lbo_vlasov_1x1v_tensor_p3(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_bcorr_lbo_vlasov_1x2v_tensor_p1(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_bcorr_lbo_vlasov_1x2v_tensor_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_bcorr_lbo_vlasov_1x3v_tensor_p1(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_bcorr_lbo_vlasov_1x3v_tensor_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_bcorr_lbo_vlasov_2x1v_tensor_p1(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_bcorr_lbo_vlasov_2x1v_tensor_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_bcorr_lbo_vlasov_2x1v_tensor_p3(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_bcorr_lbo_vlasov_2x2v_tensor_p1(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_bcorr_lbo_vlasov_2x2v_tensor_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_bcorr_lbo_vlasov_2x3v_tensor_p1(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_bcorr_lbo_vlasov_3x3v_tensor_p1(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out); 

EXTERN_C_END 
