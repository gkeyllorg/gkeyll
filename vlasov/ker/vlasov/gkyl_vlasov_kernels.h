#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH double vlasov_stream_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfx_1x1v_ser_p1(const double *w, const double *dxv, 
      const double *alpha_surf_l, const double *alpha_surf_r, 
      const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
      const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfx_1x1v_ser_p1(const double *w, const double *dxv, 
      const double *alpha_surf_edge, const double *alpha_surf_skin, 
      const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
      const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_stream_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfx_1x1v_ser_p2(const double *w, const double *dxv, 
      const double *alpha_surf_l, const double *alpha_surf_r, 
      const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
      const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfx_1x1v_ser_p2(const double *w, const double *dxv, 
      const double *alpha_surf_edge, const double *alpha_surf_skin, 
      const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
      const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_stream_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfx_1x2v_ser_p1(const double *w, const double *dxv, 
      const double *alpha_surf_l, const double *alpha_surf_r, 
      const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
      const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfx_1x2v_ser_p1(const double *w, const double *dxv, 
      const double *alpha_surf_edge, const double *alpha_surf_skin, 
      const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
      const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_stream_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfx_1x2v_ser_p2(const double *w, const double *dxv, 
      const double *alpha_surf_l, const double *alpha_surf_r, 
      const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
      const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfx_1x2v_ser_p2(const double *w, const double *dxv, 
      const double *alpha_surf_edge, const double *alpha_surf_skin, 
      const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
      const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvy_1x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvy_1x2v_ser_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_stream_vol_1x3v_ser_p1(const double *w, const double *dxv, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_vol_1x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfx_1x3v_ser_p1(const double *w, const double *dxv, 
      const double *alpha_surf_l, const double *alpha_surf_r, 
      const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
      const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfx_1x3v_ser_p1(const double *w, const double *dxv, 
      const double *alpha_surf_edge, const double *alpha_surf_skin, 
      const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
      const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_1x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_1x3v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvy_1x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvy_1x3v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvz_1x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvz_1x3v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_stream_vol_1x3v_ser_p2(const double *w, const double *dxv, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_vol_1x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfx_1x3v_ser_p2(const double *w, const double *dxv, 
      const double *alpha_surf_l, const double *alpha_surf_r, 
      const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
      const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfx_1x3v_ser_p2(const double *w, const double *dxv, 
      const double *alpha_surf_edge, const double *alpha_surf_skin, 
      const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
      const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_1x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_1x3v_ser_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvy_1x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvy_1x3v_ser_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvz_1x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvz_1x3v_ser_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_stream_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfx_2x2v_ser_p1(const double *w, const double *dxv, 
      const double *alpha_surf_l, const double *alpha_surf_r, 
      const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
      const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfx_2x2v_ser_p1(const double *w, const double *dxv, 
      const double *alpha_surf_edge, const double *alpha_surf_skin, 
      const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
      const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfy_2x2v_ser_p1(const double *w, const double *dxv, 
      const double *alpha_surf_l, const double *alpha_surf_r, 
      const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
      const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfy_2x2v_ser_p1(const double *w, const double *dxv, 
      const double *alpha_surf_edge, const double *alpha_surf_skin, 
      const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
      const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_stream_vol_2x2v_ser_p2(const double *w, const double *dxv, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_vol_2x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfx_2x2v_ser_p2(const double *w, const double *dxv, 
      const double *alpha_surf_l, const double *alpha_surf_r, 
      const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
      const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfx_2x2v_ser_p2(const double *w, const double *dxv, 
      const double *alpha_surf_edge, const double *alpha_surf_skin, 
      const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
      const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfy_2x2v_ser_p2(const double *w, const double *dxv, 
      const double *alpha_surf_l, const double *alpha_surf_r, 
      const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
      const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfy_2x2v_ser_p2(const double *w, const double *dxv, 
      const double *alpha_surf_edge, const double *alpha_surf_skin, 
      const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
      const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_stream_vol_2x3v_ser_p1(const double *w, const double *dxv, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_vol_2x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfx_2x3v_ser_p1(const double *w, const double *dxv, 
      const double *alpha_surf_l, const double *alpha_surf_r, 
      const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
      const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfx_2x3v_ser_p1(const double *w, const double *dxv, 
      const double *alpha_surf_edge, const double *alpha_surf_skin, 
      const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
      const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfy_2x3v_ser_p1(const double *w, const double *dxv, 
      const double *alpha_surf_l, const double *alpha_surf_r, 
      const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
      const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfy_2x3v_ser_p1(const double *w, const double *dxv, 
      const double *alpha_surf_edge, const double *alpha_surf_skin, 
      const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
      const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvy_2x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvy_2x3v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvz_2x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvz_2x3v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 



GKYL_CU_DH double vlasov_stream_vol_1x1v_tensor_p1(const double *w, const double *dxv, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_vol_1x1v_tensor_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfx_1x1v_tensor_p1(const double *w, const double *dxv, 
      const double *alpha_surf_l, const double *alpha_surf_r, 
      const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
      const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfx_1x1v_tensor_p1(const double *w, const double *dxv, 
      const double *alpha_surf_edge, const double *alpha_surf_skin, 
      const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
      const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_1x1v_tensor_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_1x1v_tensor_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_stream_vol_1x1v_tensor_p2(const double *w, const double *dxv, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_vol_1x1v_tensor_p2(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfx_1x1v_tensor_p2(const double *w, const double *dxv, 
      const double *alpha_surf_l, const double *alpha_surf_r, 
      const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
      const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfx_1x1v_tensor_p2(const double *w, const double *dxv, 
      const double *alpha_surf_edge, const double *alpha_surf_skin, 
      const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
      const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_1x1v_tensor_p2(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_1x1v_tensor_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_stream_vol_1x2v_tensor_p1(const double *w, const double *dxv, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_vol_1x2v_tensor_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfx_1x2v_tensor_p1(const double *w, const double *dxv, 
      const double *alpha_surf_l, const double *alpha_surf_r, 
      const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
      const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfx_1x2v_tensor_p1(const double *w, const double *dxv, 
      const double *alpha_surf_edge, const double *alpha_surf_skin, 
      const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
      const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_1x2v_tensor_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_1x2v_tensor_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvy_1x2v_tensor_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvy_1x2v_tensor_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_stream_vol_1x2v_tensor_p2(const double *w, const double *dxv, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_vol_1x2v_tensor_p2(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfx_1x2v_tensor_p2(const double *w, const double *dxv, 
      const double *alpha_surf_l, const double *alpha_surf_r, 
      const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
      const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfx_1x2v_tensor_p2(const double *w, const double *dxv, 
      const double *alpha_surf_edge, const double *alpha_surf_skin, 
      const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
      const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_1x2v_tensor_p2(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_1x2v_tensor_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvy_1x2v_tensor_p2(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvy_1x2v_tensor_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_stream_vol_1x3v_tensor_p1(const double *w, const double *dxv, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_vol_1x3v_tensor_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfx_1x3v_tensor_p1(const double *w, const double *dxv, 
      const double *alpha_surf_l, const double *alpha_surf_r, 
      const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
      const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfx_1x3v_tensor_p1(const double *w, const double *dxv, 
      const double *alpha_surf_edge, const double *alpha_surf_skin, 
      const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
      const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_1x3v_tensor_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_1x3v_tensor_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvy_1x3v_tensor_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvy_1x3v_tensor_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvz_1x3v_tensor_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvz_1x3v_tensor_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_stream_vol_1x3v_tensor_p2(const double *w, const double *dxv, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_vol_1x3v_tensor_p2(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfx_1x3v_tensor_p2(const double *w, const double *dxv, 
      const double *alpha_surf_l, const double *alpha_surf_r, 
      const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
      const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfx_1x3v_tensor_p2(const double *w, const double *dxv, 
      const double *alpha_surf_edge, const double *alpha_surf_skin, 
      const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
      const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_1x3v_tensor_p2(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_1x3v_tensor_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvy_1x3v_tensor_p2(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvy_1x3v_tensor_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvz_1x3v_tensor_p2(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvz_1x3v_tensor_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_stream_vol_2x2v_tensor_p1(const double *w, const double *dxv, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_vol_2x2v_tensor_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfx_2x2v_tensor_p1(const double *w, const double *dxv, 
      const double *alpha_surf_l, const double *alpha_surf_r, 
      const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
      const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfx_2x2v_tensor_p1(const double *w, const double *dxv, 
      const double *alpha_surf_edge, const double *alpha_surf_skin, 
      const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
      const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfy_2x2v_tensor_p1(const double *w, const double *dxv, 
      const double *alpha_surf_l, const double *alpha_surf_r, 
      const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
      const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfy_2x2v_tensor_p1(const double *w, const double *dxv, 
      const double *alpha_surf_edge, const double *alpha_surf_skin, 
      const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
      const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_2x2v_tensor_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_2x2v_tensor_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvy_2x2v_tensor_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvy_2x2v_tensor_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_stream_vol_2x2v_tensor_p2(const double *w, const double *dxv, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_vol_2x2v_tensor_p2(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfx_2x2v_tensor_p2(const double *w, const double *dxv, 
      const double *alpha_surf_l, const double *alpha_surf_r, 
      const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
      const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfx_2x2v_tensor_p2(const double *w, const double *dxv, 
      const double *alpha_surf_edge, const double *alpha_surf_skin, 
      const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
      const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfy_2x2v_tensor_p2(const double *w, const double *dxv, 
      const double *alpha_surf_l, const double *alpha_surf_r, 
      const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
      const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfy_2x2v_tensor_p2(const double *w, const double *dxv, 
      const double *alpha_surf_edge, const double *alpha_surf_skin, 
      const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
      const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_2x2v_tensor_p2(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_2x2v_tensor_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvy_2x2v_tensor_p2(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvy_2x2v_tensor_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_stream_vol_2x3v_tensor_p1(const double *w, const double *dxv, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_vol_2x3v_tensor_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfx_2x3v_tensor_p1(const double *w, const double *dxv, 
      const double *alpha_surf_l, const double *alpha_surf_r, 
      const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
      const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfx_2x3v_tensor_p1(const double *w, const double *dxv, 
      const double *alpha_surf_edge, const double *alpha_surf_skin, 
      const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
      const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfy_2x3v_tensor_p1(const double *w, const double *dxv, 
      const double *alpha_surf_l, const double *alpha_surf_r, 
      const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
      const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfy_2x3v_tensor_p1(const double *w, const double *dxv, 
      const double *alpha_surf_edge, const double *alpha_surf_skin, 
      const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
      const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_2x3v_tensor_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_2x3v_tensor_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvy_2x3v_tensor_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvy_2x3v_tensor_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvz_2x3v_tensor_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvz_2x3v_tensor_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 



EXTERN_C_END 
