#pragma once 
#include <math.h> 
#include <gkyl_mat.h> 
#include <gkyl_wave_geom.h> 
#include <gkyl_wv_eqn.h> 
#include <gkyl_wv_maxwell_priv.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 
typedef struct { double c, chi, gamma; } gkyl_maxwell_inp; 

GKYL_CU_DH double maxwell_vol_1x_ser_p1(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfx_1x_ser_p1(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double maxwell_vol_1x_ser_p2(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfx_1x_ser_p2(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double maxwell_vol_1x_ser_p3(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfx_1x_ser_p3(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void em_calc_BB_1x_ser_p1(const double *em, double* GKYL_RESTRICT out); 
GKYL_CU_DH void em_calc_num_ExB_1x_ser_p1(const double *em, double* GKYL_RESTRICT out); 
GKYL_CU_DH int em_set_bvar_1x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *BB); 
GKYL_CU_DH int em_set_ExB_1x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *num_ExB); 
GKYL_CU_DH void em_copy_bvar_1x_ser_p1(int count, struct gkyl_nmat *x, const double *em, int* cell_avg_magB2, 
    double* GKYL_RESTRICT bvar, double* GKYL_RESTRICT bvar_surf); 
GKYL_CU_DH void em_copy_ExB_1x_ser_p1(int count, struct gkyl_nmat *x, const double *em, int* cell_avg_magB2, 
    double* GKYL_RESTRICT ExB, double* GKYL_RESTRICT ExB_surf); 
GKYL_CU_DH void em_div_b_x_1x_ser_p1(const double *dxv, 
      const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
      const double *bvar_c, double* GKYL_RESTRICT max_b, double* GKYL_RESTRICT div_b); 
GKYL_CU_DH void em_vars_limiterx_1x_ser_p1(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, double *em_l, double *em_c, double *em_r); 

GKYL_CU_DH void em_calc_BB_1x_ser_p2(const double *em, double* GKYL_RESTRICT out); 
GKYL_CU_DH void em_calc_num_ExB_1x_ser_p2(const double *em, double* GKYL_RESTRICT out); 
GKYL_CU_DH int em_set_bvar_1x_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *BB); 
GKYL_CU_DH int em_set_ExB_1x_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *num_ExB); 
GKYL_CU_DH void em_copy_bvar_1x_ser_p2(int count, struct gkyl_nmat *x, const double *em, int* cell_avg_magB2, 
    double* GKYL_RESTRICT bvar, double* GKYL_RESTRICT bvar_surf); 
GKYL_CU_DH void em_copy_ExB_1x_ser_p2(int count, struct gkyl_nmat *x, const double *em, int* cell_avg_magB2, 
    double* GKYL_RESTRICT ExB, double* GKYL_RESTRICT ExB_surf); 
GKYL_CU_DH void em_div_b_x_1x_ser_p2(const double *dxv, 
      const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
      const double *bvar_c, double* GKYL_RESTRICT max_b, double* GKYL_RESTRICT div_b); 
GKYL_CU_DH void em_vars_limiterx_1x_ser_p2(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, double *em_l, double *em_c, double *em_r); 

GKYL_CU_DH void em_calc_BB_1x_ser_p3(const double *em, double* GKYL_RESTRICT out); 
GKYL_CU_DH void em_calc_num_ExB_1x_ser_p3(const double *em, double* GKYL_RESTRICT out); 
GKYL_CU_DH int em_set_bvar_1x_ser_p3(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *BB); 
GKYL_CU_DH int em_set_ExB_1x_ser_p3(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *num_ExB); 
GKYL_CU_DH void em_copy_bvar_1x_ser_p3(int count, struct gkyl_nmat *x, const double *em, int* cell_avg_magB2, 
    double* GKYL_RESTRICT bvar, double* GKYL_RESTRICT bvar_surf); 
GKYL_CU_DH void em_copy_ExB_1x_ser_p3(int count, struct gkyl_nmat *x, const double *em, int* cell_avg_magB2, 
    double* GKYL_RESTRICT ExB, double* GKYL_RESTRICT ExB_surf); 
GKYL_CU_DH void em_div_b_x_1x_ser_p3(const double *dxv, 
      const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
      const double *bvar_c, double* GKYL_RESTRICT max_b, double* GKYL_RESTRICT div_b); 
GKYL_CU_DH void em_vars_limiterx_1x_ser_p3(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, double *em_l, double *em_c, double *em_r); 

GKYL_CU_DH double maxwell_vol_2x_ser_p1(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfx_2x_ser_p1(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfy_2x_ser_p1(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double maxwell_vol_2x_ser_p2(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfx_2x_ser_p2(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfy_2x_ser_p2(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double maxwell_vol_2x_ser_p3(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfx_2x_ser_p3(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfy_2x_ser_p3(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void em_calc_BB_2x_ser_p1(const double *em, double* GKYL_RESTRICT out); 
GKYL_CU_DH void em_calc_num_ExB_2x_ser_p1(const double *em, double* GKYL_RESTRICT out); 
GKYL_CU_DH int em_set_bvar_2x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *BB); 
GKYL_CU_DH int em_set_ExB_2x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *num_ExB); 
GKYL_CU_DH void em_copy_bvar_2x_ser_p1(int count, struct gkyl_nmat *x, const double *em, int* cell_avg_magB2, 
    double* GKYL_RESTRICT bvar, double* GKYL_RESTRICT bvar_surf); 
GKYL_CU_DH void em_copy_ExB_2x_ser_p1(int count, struct gkyl_nmat *x, const double *em, int* cell_avg_magB2, 
    double* GKYL_RESTRICT ExB, double* GKYL_RESTRICT ExB_surf); 
GKYL_CU_DH void em_div_b_x_2x_ser_p1(const double *dxv, 
      const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
      const double *bvar_c, double* GKYL_RESTRICT max_b, double* GKYL_RESTRICT div_b); 
GKYL_CU_DH void em_vars_limiterx_2x_ser_p1(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, double *em_l, double *em_c, double *em_r); 
GKYL_CU_DH void em_div_b_y_2x_ser_p1(const double *dxv, 
      const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
      const double *bvar_c, double* GKYL_RESTRICT max_b, double* GKYL_RESTRICT div_b); 
GKYL_CU_DH void em_vars_limitery_2x_ser_p1(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, double *em_l, double *em_c, double *em_r); 

GKYL_CU_DH double maxwell_vol_3x_ser_p1(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfx_3x_ser_p1(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfy_3x_ser_p1(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfz_3x_ser_p1(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double maxwell_vol_3x_ser_p2(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfx_3x_ser_p2(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfy_3x_ser_p2(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfz_3x_ser_p2(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void em_calc_BB_3x_ser_p1(const double *em, double* GKYL_RESTRICT out); 
GKYL_CU_DH void em_calc_num_ExB_3x_ser_p1(const double *em, double* GKYL_RESTRICT out); 
GKYL_CU_DH int em_set_bvar_3x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *BB); 
GKYL_CU_DH int em_set_ExB_3x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *num_ExB); 
GKYL_CU_DH void em_copy_bvar_3x_ser_p1(int count, struct gkyl_nmat *x, const double *em, int* cell_avg_magB2, 
    double* GKYL_RESTRICT bvar, double* GKYL_RESTRICT bvar_surf); 
GKYL_CU_DH void em_copy_ExB_3x_ser_p1(int count, struct gkyl_nmat *x, const double *em, int* cell_avg_magB2, 
    double* GKYL_RESTRICT ExB, double* GKYL_RESTRICT ExB_surf); 
GKYL_CU_DH void em_div_b_x_3x_ser_p1(const double *dxv, 
      const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
      const double *bvar_c, double* GKYL_RESTRICT max_b, double* GKYL_RESTRICT div_b); 
GKYL_CU_DH void em_vars_limiterx_3x_ser_p1(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, double *em_l, double *em_c, double *em_r); 
GKYL_CU_DH void em_div_b_y_3x_ser_p1(const double *dxv, 
      const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
      const double *bvar_c, double* GKYL_RESTRICT max_b, double* GKYL_RESTRICT div_b); 
GKYL_CU_DH void em_vars_limitery_3x_ser_p1(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, double *em_l, double *em_c, double *em_r); 
GKYL_CU_DH void em_div_b_z_3x_ser_p1(const double *dxv, 
      const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
      const double *bvar_c, double* GKYL_RESTRICT max_b, double* GKYL_RESTRICT div_b); 
GKYL_CU_DH void em_vars_limiterz_3x_ser_p1(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, double *em_l, double *em_c, double *em_r); 

GKYL_CU_DH double maxwell_vol_2x_tensor_p2(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfx_2x_tensor_p2(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfy_2x_tensor_p2(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double maxwell_vol_2x_tensor_p3(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfx_2x_tensor_p3(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfy_2x_tensor_p3(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void em_calc_BB_2x_tensor_p2(const double *em, double* GKYL_RESTRICT out); 
GKYL_CU_DH void em_calc_num_ExB_2x_tensor_p2(const double *em, double* GKYL_RESTRICT out); 
GKYL_CU_DH int em_set_bvar_2x_tensor_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *BB); 
GKYL_CU_DH int em_set_ExB_2x_tensor_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *num_ExB); 
GKYL_CU_DH void em_copy_bvar_2x_tensor_p2(int count, struct gkyl_nmat *x, const double *em, int* cell_avg_magB2, 
    double* GKYL_RESTRICT bvar, double* GKYL_RESTRICT bvar_surf); 
GKYL_CU_DH void em_copy_ExB_2x_tensor_p2(int count, struct gkyl_nmat *x, const double *em, int* cell_avg_magB2, 
    double* GKYL_RESTRICT ExB, double* GKYL_RESTRICT ExB_surf); 
GKYL_CU_DH void em_div_b_x_2x_tensor_p2(const double *dxv, 
      const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
      const double *bvar_c, double* GKYL_RESTRICT max_b, double* GKYL_RESTRICT div_b); 
GKYL_CU_DH void em_vars_limiterx_2x_tensor_p2(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, double *em_l, double *em_c, double *em_r); 
GKYL_CU_DH void em_div_b_y_2x_tensor_p2(const double *dxv, 
      const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
      const double *bvar_c, double* GKYL_RESTRICT max_b, double* GKYL_RESTRICT div_b); 
GKYL_CU_DH void em_vars_limitery_2x_tensor_p2(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, double *em_l, double *em_c, double *em_r); 

GKYL_CU_DH void em_calc_BB_2x_tensor_p3(const double *em, double* GKYL_RESTRICT out); 
GKYL_CU_DH void em_calc_num_ExB_2x_tensor_p3(const double *em, double* GKYL_RESTRICT out); 
GKYL_CU_DH int em_set_bvar_2x_tensor_p3(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *BB); 
GKYL_CU_DH int em_set_ExB_2x_tensor_p3(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *num_ExB); 
GKYL_CU_DH void em_copy_bvar_2x_tensor_p3(int count, struct gkyl_nmat *x, const double *em, int* cell_avg_magB2, 
    double* GKYL_RESTRICT bvar, double* GKYL_RESTRICT bvar_surf); 
GKYL_CU_DH void em_copy_ExB_2x_tensor_p3(int count, struct gkyl_nmat *x, const double *em, int* cell_avg_magB2, 
    double* GKYL_RESTRICT ExB, double* GKYL_RESTRICT ExB_surf); 
GKYL_CU_DH void em_div_b_x_2x_tensor_p3(const double *dxv, 
      const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
      const double *bvar_c, double* GKYL_RESTRICT max_b, double* GKYL_RESTRICT div_b); 
GKYL_CU_DH void em_vars_limiterx_2x_tensor_p3(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, double *em_l, double *em_c, double *em_r); 
GKYL_CU_DH void em_div_b_y_2x_tensor_p3(const double *dxv, 
      const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
      const double *bvar_c, double* GKYL_RESTRICT max_b, double* GKYL_RESTRICT div_b); 
GKYL_CU_DH void em_vars_limitery_2x_tensor_p3(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, double *em_l, double *em_c, double *em_r); 

GKYL_CU_DH double maxwell_vol_3x_tensor_p2(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfx_3x_tensor_p2(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfy_3x_tensor_p2(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfz_3x_tensor_p2(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void em_calc_BB_3x_tensor_p2(const double *em, double* GKYL_RESTRICT out); 
GKYL_CU_DH void em_calc_num_ExB_3x_tensor_p2(const double *em, double* GKYL_RESTRICT out); 
GKYL_CU_DH int em_set_bvar_3x_tensor_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *BB); 
GKYL_CU_DH int em_set_ExB_3x_tensor_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *num_ExB); 
GKYL_CU_DH void em_copy_bvar_3x_tensor_p2(int count, struct gkyl_nmat *x, const double *em, int* cell_avg_magB2, 
    double* GKYL_RESTRICT bvar, double* GKYL_RESTRICT bvar_surf); 
GKYL_CU_DH void em_copy_ExB_3x_tensor_p2(int count, struct gkyl_nmat *x, const double *em, int* cell_avg_magB2, 
    double* GKYL_RESTRICT ExB, double* GKYL_RESTRICT ExB_surf); 
GKYL_CU_DH void em_div_b_x_3x_tensor_p2(const double *dxv, 
      const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
      const double *bvar_c, double* GKYL_RESTRICT max_b, double* GKYL_RESTRICT div_b); 
GKYL_CU_DH void em_vars_limiterx_3x_tensor_p2(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, double *em_l, double *em_c, double *em_r); 
GKYL_CU_DH void em_div_b_y_3x_tensor_p2(const double *dxv, 
      const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
      const double *bvar_c, double* GKYL_RESTRICT max_b, double* GKYL_RESTRICT div_b); 
GKYL_CU_DH void em_vars_limitery_3x_tensor_p2(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, double *em_l, double *em_c, double *em_r); 
GKYL_CU_DH void em_div_b_z_3x_tensor_p2(const double *dxv, 
      const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
      const double *bvar_c, double* GKYL_RESTRICT max_b, double* GKYL_RESTRICT div_b); 
GKYL_CU_DH void em_vars_limiterz_3x_tensor_p2(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, double *em_l, double *em_c, double *em_r); 

EXTERN_C_END 
