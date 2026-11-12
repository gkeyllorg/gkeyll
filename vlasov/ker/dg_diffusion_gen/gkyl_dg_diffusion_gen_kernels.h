#pragma once 
#include <math.h> 
#include <gkyl_util.h> 

EXTERN_C_BEG 

GKYL_CU_DH double dg_diffusion_gen_vol_2x_ser_p1(const double* w, const double* dx, const double* D, const double* q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfxx_2x_ser_p1(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfxy_2x_ser_p1(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfyx_2x_ser_p1(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfyy_2x_ser_p1(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 

GKYL_CU_DH double dg_diffusion_gen_vol_2x_ser_p2(const double* w, const double* dx, const double* D, const double* q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfxx_2x_ser_p2(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfxy_2x_ser_p2(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfyx_2x_ser_p2(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfyy_2x_ser_p2(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 

GKYL_CU_DH double dg_diffusion_gen_vol_2x_ser_p3(const double* w, const double* dx, const double* D, const double* q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfxx_2x_ser_p3(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfxy_2x_ser_p3(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfyx_2x_ser_p3(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfyy_2x_ser_p3(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 

GKYL_CU_DH double dg_diffusion_gen_vol_3x_ser_p1(const double* w, const double* dx, const double* D, const double* q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfxx_3x_ser_p1(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfxy_3x_ser_p1(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfxz_3x_ser_p1(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfyx_3x_ser_p1(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfyy_3x_ser_p1(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfyz_3x_ser_p1(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfzx_3x_ser_p1(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfzy_3x_ser_p1(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfzz_3x_ser_p1(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 

GKYL_CU_DH double dg_diffusion_gen_vol_3x_ser_p2(const double* w, const double* dx, const double* D, const double* q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfxx_3x_ser_p2(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfxy_3x_ser_p2(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfxz_3x_ser_p2(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfyx_3x_ser_p2(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfyy_3x_ser_p2(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfyz_3x_ser_p2(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfzx_3x_ser_p2(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfzy_3x_ser_p2(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfzz_3x_ser_p2(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 

GKYL_CU_DH double dg_diffusion_gen_vol_2x_tensor_p2(const double* w, const double* dx, const double* D, const double* q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfxx_2x_tensor_p2(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfxy_2x_tensor_p2(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfyx_2x_tensor_p2(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfyy_2x_tensor_p2(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 

GKYL_CU_DH double dg_diffusion_gen_vol_2x_tensor_p3(const double* w, const double* dx, const double* D, const double* q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfxx_2x_tensor_p3(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfxy_2x_tensor_p3(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfyx_2x_tensor_p3(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion_gen_surfyy_2x_tensor_p3(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out); 

EXTERN_C_END 
