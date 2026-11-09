#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH double vlasov_sr_stream_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_surfx_1x1v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_boundary_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_sr_stream_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_surfx_1x1v_ser_p2(const double *w, const double *dxv, const double *gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_boundary_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_sr_stream_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_surfx_1x2v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_boundary_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_boundary_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_sr_stream_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_surfx_1x2v_ser_p2(const double *w, const double *dxv, const double *gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_boundary_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_surfvy_1x2v_ser_p2(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_boundary_surfvy_1x2v_ser_p2(const double *w, const double *dxv, const double *gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_sr_stream_vol_1x3v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_vol_1x3v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_surfx_1x3v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_surfvx_1x3v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_boundary_surfvx_1x3v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_surfvy_1x3v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_boundary_surfvy_1x3v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_surfvz_1x3v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_boundary_surfvz_1x3v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_sr_stream_vol_1x3v_ser_p2(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_vol_1x3v_ser_p2(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_surfx_1x3v_ser_p2(const double *w, const double *dxv, const double *gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_surfvx_1x3v_ser_p2(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_boundary_surfvx_1x3v_ser_p2(const double *w, const double *dxv, const double *gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_surfvy_1x3v_ser_p2(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_boundary_surfvy_1x3v_ser_p2(const double *w, const double *dxv, const double *gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_surfvz_1x3v_ser_p2(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_boundary_surfvz_1x3v_ser_p2(const double *w, const double *dxv, const double *gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_sr_stream_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_surfx_2x2v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_surfy_2x2v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_boundary_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_boundary_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_sr_stream_vol_2x2v_ser_p2(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_vol_2x2v_ser_p2(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_surfx_2x2v_ser_p2(const double *w, const double *dxv, const double *gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_surfy_2x2v_ser_p2(const double *w, const double *dxv, const double *gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_boundary_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_boundary_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_sr_stream_vol_2x3v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_vol_2x3v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_surfx_2x3v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_surfy_2x3v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_boundary_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_surfvy_2x3v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_boundary_surfvy_2x3v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_surfvz_2x3v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_boundary_surfvz_2x3v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 



EXTERN_C_END 
