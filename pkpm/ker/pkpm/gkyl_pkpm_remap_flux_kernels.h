#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH void pkpm_remap_flux_trace_x_1x1v_ser_p1(double jdinv_l, double jdinv_r, const double *f_l, const double *f_r, double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r); 
GKYL_CU_DH void pkpm_remap_flux_store_x_1x1v_ser_p1(const double *wpsi, const double *proj, double* GKYL_RESTRICT flux); 

GKYL_CU_DH void pkpm_remap_flux_trace_x_2x1v_ser_p1(double jdinv_l, double jdinv_r, const double *f_l, const double *f_r, double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r); 
GKYL_CU_DH void pkpm_remap_flux_store_x_2x1v_ser_p1(const double *wpsi, const double *proj, double* GKYL_RESTRICT flux); 
GKYL_CU_DH void pkpm_remap_flux_trace_y_2x1v_ser_p1(double jdinv_l, double jdinv_r, const double *f_l, const double *f_r, double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r); 
GKYL_CU_DH void pkpm_remap_flux_store_y_2x1v_ser_p1(const double *wpsi, const double *proj, double* GKYL_RESTRICT flux); 

GKYL_CU_DH void pkpm_remap_flux_trace_x_3x1v_ser_p1(double jdinv_l, double jdinv_r, const double *f_l, const double *f_r, double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r); 
GKYL_CU_DH void pkpm_remap_flux_store_x_3x1v_ser_p1(const double *wpsi, const double *proj, double* GKYL_RESTRICT flux); 
GKYL_CU_DH void pkpm_remap_flux_trace_y_3x1v_ser_p1(double jdinv_l, double jdinv_r, const double *f_l, const double *f_r, double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r); 
GKYL_CU_DH void pkpm_remap_flux_store_y_3x1v_ser_p1(const double *wpsi, const double *proj, double* GKYL_RESTRICT flux); 
GKYL_CU_DH void pkpm_remap_flux_trace_z_3x1v_ser_p1(double jdinv_l, double jdinv_r, const double *f_l, const double *f_r, double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r); 
GKYL_CU_DH void pkpm_remap_flux_store_z_3x1v_ser_p1(const double *wpsi, const double *proj, double* GKYL_RESTRICT flux); 

EXTERN_C_END 
