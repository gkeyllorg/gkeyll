#pragma once 
#include <math.h> 
#include <stdbool.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH double pkpm_fluid_vol_1x_ser_p1(const double *w, const double *dxv, const double *jacob_pos,
  const double *vlasov_pkpm_moms, const double *bvar, const double *u, const double *fluid,
  double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_fluid_surfx_1x_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH double pkpm_fluid_vol_2x_ser_p1(const double *w, const double *dxv, const double *jacob_pos,
  const double *vlasov_pkpm_moms, const double *bvar, const double *u, const double *fluid,
  double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_fluid_surfx_2x_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_fluid_surfy_2x_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH double pkpm_fluid_vol_3x_ser_p1(const double *w, const double *dxv, const double *jacob_pos,
  const double *vlasov_pkpm_moms, const double *bvar, const double *u, const double *fluid,
  double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_fluid_surfx_3x_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_fluid_surfy_3x_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_fluid_surfz_3x_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

EXTERN_C_END 
