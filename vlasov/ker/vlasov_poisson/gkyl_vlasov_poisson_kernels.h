#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH double vlasov_poisson_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfx_1x1v_ser_p1(const double *w, const double *dxv, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfx_1x1v_ser_p1(const double *w, const double *dxv, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_boundary_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_boundary_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_poisson_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfx_1x1v_ser_p2(const double *w, const double *dxv, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfx_1x1v_ser_p2(const double *w, const double *dxv, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_boundary_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_boundary_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_poisson_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfx_1x2v_ser_p1(const double *w, const double *dxv, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfx_1x2v_ser_p1(const double *w, const double *dxv, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_boundary_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_boundary_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_boundary_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_boundary_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_poisson_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfx_1x2v_ser_p2(const double *w, const double *dxv, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfx_1x2v_ser_p2(const double *w, const double *dxv, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_boundary_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_boundary_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvy_1x2v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_surfvy_1x2v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_surfvy_1x2v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_boundary_surfvy_1x2v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_boundary_surfvy_1x2v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_poisson_vol_1x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_vol_1x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_vol_1x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfx_1x3v_ser_p1(const double *w, const double *dxv, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfvx_1x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfx_1x3v_ser_p1(const double *w, const double *dxv, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvx_1x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_surfvx_1x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_surfvx_1x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_boundary_surfvx_1x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_boundary_surfvx_1x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvy_1x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_surfvy_1x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_surfvy_1x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_boundary_surfvy_1x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_boundary_surfvy_1x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvz_1x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_surfvz_1x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_surfvz_1x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_boundary_surfvz_1x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_boundary_surfvz_1x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_poisson_vol_1x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_vol_1x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_vol_1x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfx_1x3v_ser_p2(const double *w, const double *dxv, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfvx_1x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfx_1x3v_ser_p2(const double *w, const double *dxv, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvx_1x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_surfvx_1x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_surfvx_1x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_boundary_surfvx_1x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_boundary_surfvx_1x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvy_1x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_surfvy_1x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_surfvy_1x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_boundary_surfvy_1x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_boundary_surfvy_1x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvz_1x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_surfvz_1x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_surfvz_1x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_boundary_surfvz_1x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_boundary_surfvz_1x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_poisson_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfx_2x2v_ser_p1(const double *w, const double *dxv, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfx_2x2v_ser_p1(const double *w, const double *dxv, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfy_2x2v_ser_p1(const double *w, const double *dxv, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfy_2x2v_ser_p1(const double *w, const double *dxv, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_boundary_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_boundary_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_boundary_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_boundary_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_poisson_vol_2x2v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_vol_2x2v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_vol_2x2v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfx_2x2v_ser_p2(const double *w, const double *dxv, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfx_2x2v_ser_p2(const double *w, const double *dxv, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfy_2x2v_ser_p2(const double *w, const double *dxv, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfy_2x2v_ser_p2(const double *w, const double *dxv, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_boundary_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_boundary_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_boundary_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_boundary_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_poisson_vol_2x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_vol_2x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_vol_2x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfx_2x3v_ser_p1(const double *w, const double *dxv, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfx_2x3v_ser_p1(const double *w, const double *dxv, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfy_2x3v_ser_p1(const double *w, const double *dxv, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfvy_2x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfy_2x3v_ser_p1(const double *w, const double *dxv, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_boundary_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_boundary_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvy_2x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_surfvy_2x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_surfvy_2x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_boundary_surfvy_2x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_boundary_surfvy_2x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvz_2x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_surfvz_2x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_surfvz_2x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_boundary_surfvz_2x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_boundary_surfvz_2x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_poisson_vol_2x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_vol_2x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_vol_2x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfx_2x3v_ser_p2(const double *w, const double *dxv, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfvx_2x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfx_2x3v_ser_p2(const double *w, const double *dxv, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfy_2x3v_ser_p2(const double *w, const double *dxv, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfvy_2x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfy_2x3v_ser_p2(const double *w, const double *dxv, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvx_2x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_surfvx_2x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_surfvx_2x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_boundary_surfvx_2x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_boundary_surfvx_2x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvy_2x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_surfvy_2x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_surfvy_2x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_boundary_surfvy_2x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_boundary_surfvy_2x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvz_2x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_surfvz_2x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_surfvz_2x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_phiA_boundary_surfvz_2x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_ext_EB_boundary_surfvz_2x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 


EXTERN_C_END 
