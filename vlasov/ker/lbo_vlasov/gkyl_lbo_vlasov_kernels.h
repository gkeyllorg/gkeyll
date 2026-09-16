#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH double lbo_vlasov_drag_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_diff_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_drag_vol_1x1v_ser_p3(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvx_1x1v_ser_p3(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvx_1x1v_ser_p3(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_diff_vol_1x1v_ser_p3(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvx_1x1v_ser_p3(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvx_1x1v_ser_p3(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_drag_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvy_1x2v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvy_1x2v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_diff_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvy_1x2v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvy_1x2v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_drag_vol_1x3v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvx_1x3v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvx_1x3v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvy_1x3v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvy_1x3v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvz_1x3v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvz_1x3v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_diff_vol_1x3v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvx_1x3v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvx_1x3v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvy_1x3v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvy_1x3v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvz_1x3v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvz_1x3v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_drag_vol_2x1v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvx_2x1v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvx_2x1v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_diff_vol_2x1v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvx_2x1v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvx_2x1v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_drag_vol_2x1v_ser_p3(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvx_2x1v_ser_p3(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvx_2x1v_ser_p3(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_diff_vol_2x1v_ser_p3(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvx_2x1v_ser_p3(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvx_2x1v_ser_p3(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_drag_vol_2x2v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_diff_vol_2x2v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_drag_vol_2x3v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvx_2x3v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvx_2x3v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvy_2x3v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvy_2x3v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvz_2x3v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvz_2x3v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_diff_vol_2x3v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvx_2x3v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvx_2x3v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvy_2x3v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvy_2x3v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvz_2x3v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvz_2x3v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_drag_vol_1x1v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvx_1x1v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvx_1x1v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_diff_vol_1x1v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvx_1x1v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvx_1x1v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_drag_vol_1x1v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvx_1x1v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvx_1x1v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_diff_vol_1x1v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvx_1x1v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvx_1x1v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_drag_vol_1x1v_tensor_p3(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvx_1x1v_tensor_p3(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvx_1x1v_tensor_p3(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_diff_vol_1x1v_tensor_p3(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvx_1x1v_tensor_p3(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvx_1x1v_tensor_p3(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_drag_vol_1x2v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvx_1x2v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvx_1x2v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvy_1x2v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvy_1x2v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_diff_vol_1x2v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvx_1x2v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvx_1x2v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvy_1x2v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvy_1x2v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_drag_vol_1x2v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvx_1x2v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvx_1x2v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvy_1x2v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvy_1x2v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_diff_vol_1x2v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvx_1x2v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvx_1x2v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvy_1x2v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvy_1x2v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_drag_vol_1x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvx_1x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvx_1x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvy_1x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvy_1x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvz_1x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvz_1x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_diff_vol_1x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvx_1x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvx_1x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvy_1x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvy_1x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvz_1x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvz_1x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_drag_vol_1x3v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvx_1x3v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvx_1x3v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvy_1x3v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvy_1x3v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvz_1x3v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvz_1x3v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_diff_vol_1x3v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvx_1x3v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvx_1x3v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvy_1x3v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvy_1x3v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvz_1x3v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvz_1x3v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_drag_vol_2x1v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvx_2x1v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvx_2x1v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_diff_vol_2x1v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvx_2x1v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvx_2x1v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_drag_vol_2x1v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvx_2x1v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvx_2x1v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_diff_vol_2x1v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvx_2x1v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvx_2x1v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_drag_vol_2x1v_tensor_p3(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvx_2x1v_tensor_p3(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvx_2x1v_tensor_p3(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_diff_vol_2x1v_tensor_p3(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvx_2x1v_tensor_p3(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvx_2x1v_tensor_p3(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_drag_vol_2x2v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvx_2x2v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvx_2x2v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvy_2x2v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvy_2x2v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_diff_vol_2x2v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvx_2x2v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvx_2x2v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvy_2x2v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvy_2x2v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_drag_vol_2x2v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvx_2x2v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvx_2x2v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvy_2x2v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvy_2x2v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_diff_vol_2x2v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvx_2x2v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvx_2x2v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvy_2x2v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvy_2x2v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_drag_vol_2x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvx_2x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvx_2x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvy_2x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvy_2x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvz_2x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvz_2x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_diff_vol_2x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvx_2x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvx_2x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvy_2x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvy_2x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvz_2x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvz_2x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_drag_vol_3x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvx_3x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvx_3x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvy_3x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvy_3x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvz_3x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_drag_surfvz_3x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_vlasov_diff_vol_3x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvx_3x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvx_3x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvy_3x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvy_3x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvz_3x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_vlasov_diff_surfvz_3x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

EXTERN_C_END 
