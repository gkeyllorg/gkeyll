#pragma once 
#include <math.h> 
#include <gkyl_mat.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH void sr_vars_n_set_1x1v_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *M0, const double *M1i); 
GKYL_CU_DH void sr_vars_n_copy_1x1v_ser_p1(int count, struct gkyl_nmat *x, const double *M0, double* GKYL_RESTRICT n); 
GKYL_CU_DH void sr_vars_GammaV_1x1v_ser_p1(const double *u_i, double* GKYL_RESTRICT u_i_sq, double* GKYL_RESTRICT GammaV, double* GKYL_RESTRICT GammaV_sq); 
GKYL_CU_DH void sr_vars_pressure_1x1v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure); 

GKYL_CU_DH void sr_vars_n_set_1x1v_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *M0, const double *M1i); 
GKYL_CU_DH void sr_vars_n_copy_1x1v_ser_p2(int count, struct gkyl_nmat *x, const double *M0, double* GKYL_RESTRICT n); 
GKYL_CU_DH void sr_vars_GammaV_1x1v_ser_p2(const double *u_i, double* GKYL_RESTRICT u_i_sq, double* GKYL_RESTRICT GammaV, double* GKYL_RESTRICT GammaV_sq); 
GKYL_CU_DH void sr_vars_pressure_1x1v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure); 

GKYL_CU_DH void sr_vars_n_set_1x1v_ser_p3(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *M0, const double *M1i); 
GKYL_CU_DH void sr_vars_n_copy_1x1v_ser_p3(int count, struct gkyl_nmat *x, const double *M0, double* GKYL_RESTRICT n); 
GKYL_CU_DH void sr_vars_GammaV_1x1v_ser_p3(const double *u_i, double* GKYL_RESTRICT u_i_sq, double* GKYL_RESTRICT GammaV, double* GKYL_RESTRICT GammaV_sq); 
GKYL_CU_DH void sr_vars_pressure_1x1v_ser_p3(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure); 

GKYL_CU_DH void sr_vars_n_set_1x2v_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *M0, const double *M1i); 
GKYL_CU_DH void sr_vars_n_copy_1x2v_ser_p1(int count, struct gkyl_nmat *x, const double *M0, double* GKYL_RESTRICT n); 
GKYL_CU_DH void sr_vars_GammaV_1x2v_ser_p1(const double *u_i, double* GKYL_RESTRICT u_i_sq, double* GKYL_RESTRICT GammaV, double* GKYL_RESTRICT GammaV_sq); 
GKYL_CU_DH void sr_vars_pressure_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure); 

GKYL_CU_DH void sr_vars_n_set_1x2v_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *M0, const double *M1i); 
GKYL_CU_DH void sr_vars_n_copy_1x2v_ser_p2(int count, struct gkyl_nmat *x, const double *M0, double* GKYL_RESTRICT n); 
GKYL_CU_DH void sr_vars_GammaV_1x2v_ser_p2(const double *u_i, double* GKYL_RESTRICT u_i_sq, double* GKYL_RESTRICT GammaV, double* GKYL_RESTRICT GammaV_sq); 
GKYL_CU_DH void sr_vars_pressure_1x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure); 

GKYL_CU_DH void sr_vars_n_set_1x3v_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *M0, const double *M1i); 
GKYL_CU_DH void sr_vars_n_copy_1x3v_ser_p1(int count, struct gkyl_nmat *x, const double *M0, double* GKYL_RESTRICT n); 
GKYL_CU_DH void sr_vars_GammaV_1x3v_ser_p1(const double *u_i, double* GKYL_RESTRICT u_i_sq, double* GKYL_RESTRICT GammaV, double* GKYL_RESTRICT GammaV_sq); 
GKYL_CU_DH void sr_vars_pressure_1x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure); 

GKYL_CU_DH void sr_vars_n_set_1x3v_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *M0, const double *M1i); 
GKYL_CU_DH void sr_vars_n_copy_1x3v_ser_p2(int count, struct gkyl_nmat *x, const double *M0, double* GKYL_RESTRICT n); 
GKYL_CU_DH void sr_vars_GammaV_1x3v_ser_p2(const double *u_i, double* GKYL_RESTRICT u_i_sq, double* GKYL_RESTRICT GammaV, double* GKYL_RESTRICT GammaV_sq); 
GKYL_CU_DH void sr_vars_pressure_1x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure); 

GKYL_CU_DH void sr_vars_n_set_2x1v_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *M0, const double *M1i); 
GKYL_CU_DH void sr_vars_n_copy_2x1v_ser_p1(int count, struct gkyl_nmat *x, const double *M0, double* GKYL_RESTRICT n); 
GKYL_CU_DH void sr_vars_GammaV_2x1v_ser_p1(const double *u_i, double* GKYL_RESTRICT u_i_sq, double* GKYL_RESTRICT GammaV, double* GKYL_RESTRICT GammaV_sq); 
GKYL_CU_DH void sr_vars_pressure_2x1v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure); 

GKYL_CU_DH void sr_vars_n_set_2x1v_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *M0, const double *M1i); 
GKYL_CU_DH void sr_vars_n_copy_2x1v_ser_p2(int count, struct gkyl_nmat *x, const double *M0, double* GKYL_RESTRICT n); 
GKYL_CU_DH void sr_vars_GammaV_2x1v_ser_p2(const double *u_i, double* GKYL_RESTRICT u_i_sq, double* GKYL_RESTRICT GammaV, double* GKYL_RESTRICT GammaV_sq); 
GKYL_CU_DH void sr_vars_pressure_2x1v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure); 

GKYL_CU_DH void sr_vars_n_set_2x1v_ser_p3(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *M0, const double *M1i); 
GKYL_CU_DH void sr_vars_n_copy_2x1v_ser_p3(int count, struct gkyl_nmat *x, const double *M0, double* GKYL_RESTRICT n); 
GKYL_CU_DH void sr_vars_GammaV_2x1v_ser_p3(const double *u_i, double* GKYL_RESTRICT u_i_sq, double* GKYL_RESTRICT GammaV, double* GKYL_RESTRICT GammaV_sq); 
GKYL_CU_DH void sr_vars_pressure_2x1v_ser_p3(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure); 

GKYL_CU_DH void sr_vars_n_set_2x2v_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *M0, const double *M1i); 
GKYL_CU_DH void sr_vars_n_copy_2x2v_ser_p1(int count, struct gkyl_nmat *x, const double *M0, double* GKYL_RESTRICT n); 
GKYL_CU_DH void sr_vars_GammaV_2x2v_ser_p1(const double *u_i, double* GKYL_RESTRICT u_i_sq, double* GKYL_RESTRICT GammaV, double* GKYL_RESTRICT GammaV_sq); 
GKYL_CU_DH void sr_vars_pressure_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure); 

GKYL_CU_DH void sr_vars_n_set_2x2v_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *M0, const double *M1i); 
GKYL_CU_DH void sr_vars_n_copy_2x2v_ser_p2(int count, struct gkyl_nmat *x, const double *M0, double* GKYL_RESTRICT n); 
GKYL_CU_DH void sr_vars_GammaV_2x2v_ser_p2(const double *u_i, double* GKYL_RESTRICT u_i_sq, double* GKYL_RESTRICT GammaV, double* GKYL_RESTRICT GammaV_sq); 
GKYL_CU_DH void sr_vars_pressure_2x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure); 

GKYL_CU_DH void sr_vars_n_set_2x3v_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *M0, const double *M1i); 
GKYL_CU_DH void sr_vars_n_copy_2x3v_ser_p1(int count, struct gkyl_nmat *x, const double *M0, double* GKYL_RESTRICT n); 
GKYL_CU_DH void sr_vars_GammaV_2x3v_ser_p1(const double *u_i, double* GKYL_RESTRICT u_i_sq, double* GKYL_RESTRICT GammaV, double* GKYL_RESTRICT GammaV_sq); 
GKYL_CU_DH void sr_vars_pressure_2x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure); 

GKYL_CU_DH void sr_vars_n_set_2x3v_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *M0, const double *M1i); 
GKYL_CU_DH void sr_vars_n_copy_2x3v_ser_p2(int count, struct gkyl_nmat *x, const double *M0, double* GKYL_RESTRICT n); 
GKYL_CU_DH void sr_vars_GammaV_2x3v_ser_p2(const double *u_i, double* GKYL_RESTRICT u_i_sq, double* GKYL_RESTRICT GammaV, double* GKYL_RESTRICT GammaV_sq); 
GKYL_CU_DH void sr_vars_pressure_2x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure); 

GKYL_CU_DH void sr_vars_n_set_3x3v_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *M0, const double *M1i); 
GKYL_CU_DH void sr_vars_n_copy_3x3v_ser_p1(int count, struct gkyl_nmat *x, const double *M0, double* GKYL_RESTRICT n); 
GKYL_CU_DH void sr_vars_GammaV_3x3v_ser_p1(const double *u_i, double* GKYL_RESTRICT u_i_sq, double* GKYL_RESTRICT GammaV, double* GKYL_RESTRICT GammaV_sq); 
GKYL_CU_DH void sr_vars_pressure_3x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure); 

GKYL_CU_DH void sr_vars_pressure_1x1v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure); 

GKYL_CU_DH void sr_vars_n_set_1x1v_tensor_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *M0, const double *M1i); 
GKYL_CU_DH void sr_vars_n_copy_1x1v_tensor_p2(int count, struct gkyl_nmat *x, const double *M0, double* GKYL_RESTRICT n); 
GKYL_CU_DH void sr_vars_GammaV_1x1v_tensor_p2(const double *u_i, double* GKYL_RESTRICT u_i_sq, double* GKYL_RESTRICT GammaV, double* GKYL_RESTRICT GammaV_sq); 
GKYL_CU_DH void sr_vars_pressure_1x1v_tensor_p2(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure); 

GKYL_CU_DH void sr_vars_n_set_1x1v_tensor_p3(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *M0, const double *M1i); 
GKYL_CU_DH void sr_vars_n_copy_1x1v_tensor_p3(int count, struct gkyl_nmat *x, const double *M0, double* GKYL_RESTRICT n); 
GKYL_CU_DH void sr_vars_GammaV_1x1v_tensor_p3(const double *u_i, double* GKYL_RESTRICT u_i_sq, double* GKYL_RESTRICT GammaV, double* GKYL_RESTRICT GammaV_sq); 
GKYL_CU_DH void sr_vars_pressure_1x1v_tensor_p3(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure); 

GKYL_CU_DH void sr_vars_pressure_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure); 

GKYL_CU_DH void sr_vars_n_set_1x2v_tensor_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *M0, const double *M1i); 
GKYL_CU_DH void sr_vars_n_copy_1x2v_tensor_p2(int count, struct gkyl_nmat *x, const double *M0, double* GKYL_RESTRICT n); 
GKYL_CU_DH void sr_vars_GammaV_1x2v_tensor_p2(const double *u_i, double* GKYL_RESTRICT u_i_sq, double* GKYL_RESTRICT GammaV, double* GKYL_RESTRICT GammaV_sq); 
GKYL_CU_DH void sr_vars_pressure_1x2v_tensor_p2(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure); 

GKYL_CU_DH void sr_vars_pressure_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure); 

GKYL_CU_DH void sr_vars_n_set_1x3v_tensor_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *M0, const double *M1i); 
GKYL_CU_DH void sr_vars_n_copy_1x3v_tensor_p2(int count, struct gkyl_nmat *x, const double *M0, double* GKYL_RESTRICT n); 
GKYL_CU_DH void sr_vars_GammaV_1x3v_tensor_p2(const double *u_i, double* GKYL_RESTRICT u_i_sq, double* GKYL_RESTRICT GammaV, double* GKYL_RESTRICT GammaV_sq); 
GKYL_CU_DH void sr_vars_pressure_1x3v_tensor_p2(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure); 

GKYL_CU_DH void sr_vars_pressure_2x1v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure); 

GKYL_CU_DH void sr_vars_n_set_2x1v_tensor_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *M0, const double *M1i); 
GKYL_CU_DH void sr_vars_n_copy_2x1v_tensor_p2(int count, struct gkyl_nmat *x, const double *M0, double* GKYL_RESTRICT n); 
GKYL_CU_DH void sr_vars_GammaV_2x1v_tensor_p2(const double *u_i, double* GKYL_RESTRICT u_i_sq, double* GKYL_RESTRICT GammaV, double* GKYL_RESTRICT GammaV_sq); 
GKYL_CU_DH void sr_vars_pressure_2x1v_tensor_p2(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure); 

GKYL_CU_DH void sr_vars_n_set_2x1v_tensor_p3(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *M0, const double *M1i); 
GKYL_CU_DH void sr_vars_n_copy_2x1v_tensor_p3(int count, struct gkyl_nmat *x, const double *M0, double* GKYL_RESTRICT n); 
GKYL_CU_DH void sr_vars_GammaV_2x1v_tensor_p3(const double *u_i, double* GKYL_RESTRICT u_i_sq, double* GKYL_RESTRICT GammaV, double* GKYL_RESTRICT GammaV_sq); 
GKYL_CU_DH void sr_vars_pressure_2x1v_tensor_p3(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure); 

GKYL_CU_DH void sr_vars_pressure_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure); 

GKYL_CU_DH void sr_vars_n_set_2x2v_tensor_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *M0, const double *M1i); 
GKYL_CU_DH void sr_vars_n_copy_2x2v_tensor_p2(int count, struct gkyl_nmat *x, const double *M0, double* GKYL_RESTRICT n); 
GKYL_CU_DH void sr_vars_GammaV_2x2v_tensor_p2(const double *u_i, double* GKYL_RESTRICT u_i_sq, double* GKYL_RESTRICT GammaV, double* GKYL_RESTRICT GammaV_sq); 
GKYL_CU_DH void sr_vars_pressure_2x2v_tensor_p2(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure); 

GKYL_CU_DH void sr_vars_pressure_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure); 

GKYL_CU_DH void sr_vars_n_set_2x3v_tensor_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *M0, const double *M1i); 
GKYL_CU_DH void sr_vars_n_copy_2x3v_tensor_p2(int count, struct gkyl_nmat *x, const double *M0, double* GKYL_RESTRICT n); 
GKYL_CU_DH void sr_vars_GammaV_2x3v_tensor_p2(const double *u_i, double* GKYL_RESTRICT u_i_sq, double* GKYL_RESTRICT GammaV, double* GKYL_RESTRICT GammaV_sq); 
GKYL_CU_DH void sr_vars_pressure_2x3v_tensor_p2(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure); 

GKYL_CU_DH void sr_vars_pressure_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure); 

EXTERN_C_END 
