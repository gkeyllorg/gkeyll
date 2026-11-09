#pragma once 
#include <math.h> 
#include <stdbool.h> 
#include <gkyl_util.h> 
#include <gkyl_mat.h> 
EXTERN_C_BEG 

GKYL_CU_DH void pkpm_vars_grad_u_x_1x_ser_p1(const double *dxv, const double *jacob_pos, const double *u_l, const double *u_c, const double *u_r, double* GKYL_RESTRICT grad_u); 
GKYL_CU_DH void pkpm_vars_in_cell_grad_x_1x_ser_p1(const double *dxv, const double *jacob_pos, double blend_lo, double blend_up, const double *q_l, const double *q_c, const double *q_r, double* GKYL_RESTRICT grad); 
GKYL_CU_DH double pkpm_mom_diffusion_order2_x_1x_ser_p1(const double *dxv, double D, const double *jac_l, const double *jac_c,
  const double *jac_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_mom_diffusion_order4_x_1x_ser_p1(const double *dxv, double D, const double *jac_l, const double *jac_c,
  const double *jac_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH int pkpm_vars_closure_1x_ser_p1(double w_lo, double w_up, double dv, const double *norm_moms, const double *rhs_moms,
  const double *bvar, const double *wn, const double *theta, const double *alpha_0_prev, const double *c_h_prev,
  const double *f_lo, const double *f_up,
  double* GKYL_RESTRICT c_0, double* GKYL_RESTRICT alpha_0, double* GKYL_RESTRICT c_h, double* GKYL_RESTRICT a_perp,
  double* GKYL_RESTRICT flux_lo, double* GKYL_RESTRICT flux_up, double* GKYL_RESTRICT cfl, double* GKYL_RESTRICT cond); 
GKYL_CU_DH int pkpm_vars_grad_ln_n_1x_ser_p1(const double *rho, const double *grad_rho, double* GKYL_RESTRICT grad_ln_n); 
GKYL_CU_DH int pkpm_vars_prim_1x_ser_p1(const double *fluid, double* GKYL_RESTRICT u, double* GKYL_RESTRICT xi); 
GKYL_CU_DH int pkpm_dtB_over_magB_1x_ser_p1(const double *em, const double *dt_em, double* GKYL_RESTRICT out); 
GKYL_CU_DH int pkpm_fluid_positivity_1x_ser_p1(double eps, double* GKYL_RESTRICT fluid); 
GKYL_CU_DH void pkpm_vars_int_1x_ser_p1(const double *fluid, const double *vlasov_pkpm_moms, const double *u, double* GKYL_RESTRICT int_pkpm_vars); 
GKYL_CU_DH void pkpm_vars_source_1x_ser_p1(const double *qmem, const double *fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_diag_1x_ser_p1(const double *dxv, const double *vlasov_pkpm_moms, const double *euler_pkpm,
  const double *bvar, const double *grad_b, const double *u, const double *grad_u, const double *theta, const double *alpha_0,
  const double *xi, const double *c_h, double* GKYL_RESTRICT fluid_io, double* GKYL_RESTRICT pkpm_vars_io); 

GKYL_CU_DH void pkpm_vars_grad_u_x_2x_ser_p1(const double *dxv, const double *jacob_pos, const double *u_l, const double *u_c, const double *u_r, double* GKYL_RESTRICT grad_u); 
GKYL_CU_DH void pkpm_vars_grad_u_y_2x_ser_p1(const double *dxv, const double *jacob_pos, const double *u_l, const double *u_c, const double *u_r, double* GKYL_RESTRICT grad_u); 
GKYL_CU_DH void pkpm_vars_in_cell_grad_x_2x_ser_p1(const double *dxv, const double *jacob_pos, double blend_lo, double blend_up, const double *q_l, const double *q_c, const double *q_r, double* GKYL_RESTRICT grad); 
GKYL_CU_DH void pkpm_vars_in_cell_grad_y_2x_ser_p1(const double *dxv, const double *jacob_pos, double blend_lo, double blend_up, const double *q_l, const double *q_c, const double *q_r, double* GKYL_RESTRICT grad); 
GKYL_CU_DH double pkpm_mom_diffusion_order2_x_2x_ser_p1(const double *dxv, double D, const double *jac_l, const double *jac_c,
  const double *jac_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_mom_diffusion_order2_y_2x_ser_p1(const double *dxv, double D, const double *jac_l, const double *jac_c,
  const double *jac_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_mom_diffusion_order4_x_2x_ser_p1(const double *dxv, double D, const double *jac_l, const double *jac_c,
  const double *jac_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_mom_diffusion_order4_y_2x_ser_p1(const double *dxv, double D, const double *jac_l, const double *jac_c,
  const double *jac_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH int pkpm_vars_closure_2x_ser_p1(double w_lo, double w_up, double dv, const double *norm_moms, const double *rhs_moms,
  const double *bvar, const double *wn, const double *theta, const double *alpha_0_prev, const double *c_h_prev,
  const double *f_lo, const double *f_up,
  double* GKYL_RESTRICT c_0, double* GKYL_RESTRICT alpha_0, double* GKYL_RESTRICT c_h, double* GKYL_RESTRICT a_perp,
  double* GKYL_RESTRICT flux_lo, double* GKYL_RESTRICT flux_up, double* GKYL_RESTRICT cfl, double* GKYL_RESTRICT cond); 
GKYL_CU_DH int pkpm_vars_grad_ln_n_2x_ser_p1(const double *rho, const double *grad_rho, double* GKYL_RESTRICT grad_ln_n); 
GKYL_CU_DH int pkpm_vars_prim_2x_ser_p1(const double *fluid, double* GKYL_RESTRICT u, double* GKYL_RESTRICT xi); 
GKYL_CU_DH int pkpm_dtB_over_magB_2x_ser_p1(const double *em, const double *dt_em, double* GKYL_RESTRICT out); 
GKYL_CU_DH int pkpm_fluid_positivity_2x_ser_p1(double eps, double* GKYL_RESTRICT fluid); 
GKYL_CU_DH void pkpm_vars_int_2x_ser_p1(const double *fluid, const double *vlasov_pkpm_moms, const double *u, double* GKYL_RESTRICT int_pkpm_vars); 
GKYL_CU_DH void pkpm_vars_source_2x_ser_p1(const double *qmem, const double *fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_diag_2x_ser_p1(const double *dxv, const double *vlasov_pkpm_moms, const double *euler_pkpm,
  const double *bvar, const double *grad_b, const double *u, const double *grad_u, const double *theta, const double *alpha_0,
  const double *xi, const double *c_h, double* GKYL_RESTRICT fluid_io, double* GKYL_RESTRICT pkpm_vars_io); 

GKYL_CU_DH void pkpm_vars_grad_u_x_3x_ser_p1(const double *dxv, const double *jacob_pos, const double *u_l, const double *u_c, const double *u_r, double* GKYL_RESTRICT grad_u); 
GKYL_CU_DH void pkpm_vars_grad_u_y_3x_ser_p1(const double *dxv, const double *jacob_pos, const double *u_l, const double *u_c, const double *u_r, double* GKYL_RESTRICT grad_u); 
GKYL_CU_DH void pkpm_vars_grad_u_z_3x_ser_p1(const double *dxv, const double *jacob_pos, const double *u_l, const double *u_c, const double *u_r, double* GKYL_RESTRICT grad_u); 
GKYL_CU_DH void pkpm_vars_in_cell_grad_x_3x_ser_p1(const double *dxv, const double *jacob_pos, double blend_lo, double blend_up, const double *q_l, const double *q_c, const double *q_r, double* GKYL_RESTRICT grad); 
GKYL_CU_DH void pkpm_vars_in_cell_grad_y_3x_ser_p1(const double *dxv, const double *jacob_pos, double blend_lo, double blend_up, const double *q_l, const double *q_c, const double *q_r, double* GKYL_RESTRICT grad); 
GKYL_CU_DH void pkpm_vars_in_cell_grad_z_3x_ser_p1(const double *dxv, const double *jacob_pos, double blend_lo, double blend_up, const double *q_l, const double *q_c, const double *q_r, double* GKYL_RESTRICT grad); 
GKYL_CU_DH double pkpm_mom_diffusion_order2_x_3x_ser_p1(const double *dxv, double D, const double *jac_l, const double *jac_c,
  const double *jac_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_mom_diffusion_order2_y_3x_ser_p1(const double *dxv, double D, const double *jac_l, const double *jac_c,
  const double *jac_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_mom_diffusion_order2_z_3x_ser_p1(const double *dxv, double D, const double *jac_l, const double *jac_c,
  const double *jac_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_mom_diffusion_order4_x_3x_ser_p1(const double *dxv, double D, const double *jac_l, const double *jac_c,
  const double *jac_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_mom_diffusion_order4_y_3x_ser_p1(const double *dxv, double D, const double *jac_l, const double *jac_c,
  const double *jac_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_mom_diffusion_order4_z_3x_ser_p1(const double *dxv, double D, const double *jac_l, const double *jac_c,
  const double *jac_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH int pkpm_vars_closure_3x_ser_p1(double w_lo, double w_up, double dv, const double *norm_moms, const double *rhs_moms,
  const double *bvar, const double *wn, const double *theta, const double *alpha_0_prev, const double *c_h_prev,
  const double *f_lo, const double *f_up,
  double* GKYL_RESTRICT c_0, double* GKYL_RESTRICT alpha_0, double* GKYL_RESTRICT c_h, double* GKYL_RESTRICT a_perp,
  double* GKYL_RESTRICT flux_lo, double* GKYL_RESTRICT flux_up, double* GKYL_RESTRICT cfl, double* GKYL_RESTRICT cond); 
GKYL_CU_DH int pkpm_vars_grad_ln_n_3x_ser_p1(const double *rho, const double *grad_rho, double* GKYL_RESTRICT grad_ln_n); 
GKYL_CU_DH int pkpm_vars_prim_3x_ser_p1(const double *fluid, double* GKYL_RESTRICT u, double* GKYL_RESTRICT xi); 
GKYL_CU_DH int pkpm_dtB_over_magB_3x_ser_p1(const double *em, const double *dt_em, double* GKYL_RESTRICT out); 
GKYL_CU_DH int pkpm_fluid_positivity_3x_ser_p1(double eps, double* GKYL_RESTRICT fluid); 
GKYL_CU_DH void pkpm_vars_int_3x_ser_p1(const double *fluid, const double *vlasov_pkpm_moms, const double *u, double* GKYL_RESTRICT int_pkpm_vars); 
GKYL_CU_DH void pkpm_vars_source_3x_ser_p1(const double *qmem, const double *fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_diag_3x_ser_p1(const double *dxv, const double *vlasov_pkpm_moms, const double *euler_pkpm,
  const double *bvar, const double *grad_b, const double *u, const double *grad_u, const double *theta, const double *alpha_0,
  const double *xi, const double *c_h, double* GKYL_RESTRICT fluid_io, double* GKYL_RESTRICT pkpm_vars_io); 

EXTERN_C_END 
