#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH double advection_vol_1x_ser_p1(const double *w, const double *dxv, const double *u, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double advection_surfx_1x_ser_p1(const double *w, const double *dxv, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double advection_vol_1x_ser_p2(const double *w, const double *dxv, const double *u, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double advection_surfx_1x_ser_p2(const double *w, const double *dxv, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double advection_vol_1x_ser_p3(const double *w, const double *dxv, const double *u, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double advection_surfx_1x_ser_p3(const double *w, const double *dxv, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double advection_vol_2x_ser_p1(const double *w, const double *dxv, const double *u, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double advection_surfx_2x_ser_p1(const double *w, const double *dxv, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double advection_surfy_2x_ser_p1(const double *w, const double *dxv, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double advection_vol_2x_ser_p2(const double *w, const double *dxv, const double *u, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double advection_surfx_2x_ser_p2(const double *w, const double *dxv, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double advection_surfy_2x_ser_p2(const double *w, const double *dxv, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double advection_vol_2x_ser_p3(const double *w, const double *dxv, const double *u, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double advection_surfx_2x_ser_p3(const double *w, const double *dxv, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double advection_surfy_2x_ser_p3(const double *w, const double *dxv, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double advection_vol_3x_ser_p1(const double *w, const double *dxv, const double *u, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double advection_surfx_3x_ser_p1(const double *w, const double *dxv, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double advection_surfy_3x_ser_p1(const double *w, const double *dxv, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double advection_surfz_3x_ser_p1(const double *w, const double *dxv, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double advection_vol_3x_ser_p2(const double *w, const double *dxv, const double *u, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double advection_surfx_3x_ser_p2(const double *w, const double *dxv, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double advection_surfy_3x_ser_p2(const double *w, const double *dxv, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double advection_surfz_3x_ser_p2(const double *w, const double *dxv, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double advection_vol_2x_tensor_p2(const double *w, const double *dxv, const double *u, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double advection_surfx_2x_tensor_p2(const double *w, const double *dxv, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double advection_surfy_2x_tensor_p2(const double *w, const double *dxv, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double advection_vol_2x_tensor_p3(const double *w, const double *dxv, const double *u, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double advection_surfx_2x_tensor_p3(const double *w, const double *dxv, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double advection_surfy_2x_tensor_p3(const double *w, const double *dxv, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double advection_vol_3x_tensor_p2(const double *w, const double *dxv, const double *u, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double advection_surfx_3x_tensor_p2(const double *w, const double *dxv, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double advection_surfy_3x_tensor_p2(const double *w, const double *dxv, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double advection_surfz_3x_tensor_p2(const double *w, const double *dxv, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

EXTERN_C_END 
