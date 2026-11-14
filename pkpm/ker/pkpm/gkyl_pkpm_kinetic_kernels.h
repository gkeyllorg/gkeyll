#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH double pkpm_kinetic_conf_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *jacob_pos,
  const double *bvar, const double *u, const double *xi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH int pkpm_dist_normalize_1x1v_ser_p1(const double *nmoms, const double *f, double* GKYL_RESTRICT fbar); 
GKYL_CU_DH void pkpm_dist_transverse_1x1v_ser_p1(const double *bvar, double* GKYL_RESTRICT f); 
GKYL_CU_DH void pkpm_harm_exp_to_nodes_1x1v_ser_p1(const double *f, double* GKYL_RESTRICT fn); 
GKYL_CU_DH void pkpm_harm_exp_to_modal_1x1v_ser_p1(const double *fn, double* GKYL_RESTRICT f); 
GKYL_CU_DH void pkpm_harm_exp_fg_to_nodes_1x1v_ser_p1(const double *f, double* GKYL_RESTRICT fgn); 
GKYL_CU_DH void pkpm_harm_exp_wint_1x1v_ser_p1(const double *f, double* GKYL_RESTRICT wi); 
GKYL_CU_DH double pkpm_stream_surfx_1x1v_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_accel_surfvpar_1x1v_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vel_flux_nodes_1x_ser_p1(const double *wn, const double *theta, const double *bvar,
  const double *alpha_0, const double *c_h, double* GKYL_RESTRICT nv); 
GKYL_CU_DH void pkpm_vel_flux_vpar_1x1v_ser_p1(const double *w, const double *dxv, const double *nv,
  const double *f_l, const double *f_r, double* GKYL_RESTRICT flux); 
GKYL_CU_DH double pkpm_kinetic_vel_vol_1x1v_ser_p1(const double *w, const double *dxv,
  const double *c_0, const double *alpha_0, const double *c_h, const double *a_perp,
  const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_harm_build_1x_ser_p1(double qbym, const double *bvar, const double *grad_b, const double *u,
  const double *grad_u, const double *xi, const double *grad_xi, const double *grad_ln_n, const double *theta,
  const double *dtb, const double *em, double* GKYL_RESTRICT ho, double* GKYL_RESTRICT wn); 
GKYL_CU_DH void pkpm_harm_mat_1x_ser_p1(const double *dxv, const double *jacob_pos, const double *ho,
  double* GKYL_RESTRICT mat); 
GKYL_CU_DH void pkpm_harm_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *mat, const double *f,
  double* GKYL_RESTRICT out); 

GKYL_CU_DH double pkpm_kinetic_conf_vol_2x1v_ser_p1(const double *w, const double *dxv, const double *jacob_pos,
  const double *bvar, const double *u, const double *xi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH int pkpm_dist_normalize_2x1v_ser_p1(const double *nmoms, const double *f, double* GKYL_RESTRICT fbar); 
GKYL_CU_DH void pkpm_dist_transverse_2x1v_ser_p1(const double *bvar, double* GKYL_RESTRICT f); 
GKYL_CU_DH void pkpm_harm_exp_to_nodes_2x1v_ser_p1(const double *f, double* GKYL_RESTRICT fn); 
GKYL_CU_DH void pkpm_harm_exp_to_modal_2x1v_ser_p1(const double *fn, double* GKYL_RESTRICT f); 
GKYL_CU_DH void pkpm_harm_exp_fg_to_nodes_2x1v_ser_p1(const double *f, double* GKYL_RESTRICT fgn); 
GKYL_CU_DH void pkpm_harm_exp_wint_2x1v_ser_p1(const double *f, double* GKYL_RESTRICT wi); 
GKYL_CU_DH double pkpm_stream_surfx_2x1v_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_stream_surfy_2x1v_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_accel_surfvpar_2x1v_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vel_flux_nodes_2x_ser_p1(const double *wn, const double *theta, const double *bvar,
  const double *alpha_0, const double *c_h, double* GKYL_RESTRICT nv); 
GKYL_CU_DH void pkpm_vel_flux_vpar_2x1v_ser_p1(const double *w, const double *dxv, const double *nv,
  const double *f_l, const double *f_r, double* GKYL_RESTRICT flux); 
GKYL_CU_DH double pkpm_kinetic_vel_vol_2x1v_ser_p1(const double *w, const double *dxv,
  const double *c_0, const double *alpha_0, const double *c_h, const double *a_perp,
  const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_harm_build_2x_ser_p1(double qbym, const double *bvar, const double *grad_b, const double *u,
  const double *grad_u, const double *xi, const double *grad_xi, const double *grad_ln_n, const double *theta,
  const double *dtb, const double *em, double* GKYL_RESTRICT ho, double* GKYL_RESTRICT wn); 
GKYL_CU_DH void pkpm_harm_mat_2x_ser_p1(const double *dxv, const double *jacob_pos, const double *ho,
  double* GKYL_RESTRICT mat); 
GKYL_CU_DH void pkpm_harm_vol_2x1v_ser_p1(const double *w, const double *dxv, const double *mat, const double *f,
  double* GKYL_RESTRICT out); 

GKYL_CU_DH double pkpm_kinetic_conf_vol_3x1v_ser_p1(const double *w, const double *dxv, const double *jacob_pos,
  const double *bvar, const double *u, const double *xi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH int pkpm_dist_normalize_3x1v_ser_p1(const double *nmoms, const double *f, double* GKYL_RESTRICT fbar); 
GKYL_CU_DH void pkpm_dist_transverse_3x1v_ser_p1(const double *bvar, double* GKYL_RESTRICT f); 
GKYL_CU_DH void pkpm_harm_exp_to_nodes_3x1v_ser_p1(const double *f, double* GKYL_RESTRICT fn); 
GKYL_CU_DH void pkpm_harm_exp_to_modal_3x1v_ser_p1(const double *fn, double* GKYL_RESTRICT f); 
GKYL_CU_DH void pkpm_harm_exp_fg_to_nodes_3x1v_ser_p1(const double *f, double* GKYL_RESTRICT fgn); 
GKYL_CU_DH void pkpm_harm_exp_wint_3x1v_ser_p1(const double *f, double* GKYL_RESTRICT wi); 
GKYL_CU_DH double pkpm_stream_surfx_3x1v_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_stream_surfy_3x1v_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_stream_surfz_3x1v_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_accel_surfvpar_3x1v_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vel_flux_nodes_3x_ser_p1(const double *wn, const double *theta, const double *bvar,
  const double *alpha_0, const double *c_h, double* GKYL_RESTRICT nv); 
GKYL_CU_DH void pkpm_vel_flux_vpar_3x1v_ser_p1(const double *w, const double *dxv, const double *nv,
  const double *f_l, const double *f_r, double* GKYL_RESTRICT flux); 
GKYL_CU_DH double pkpm_kinetic_vel_vol_3x1v_ser_p1(const double *w, const double *dxv,
  const double *c_0, const double *alpha_0, const double *c_h, const double *a_perp,
  const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_harm_build_3x_ser_p1(double qbym, const double *bvar, const double *grad_b, const double *u,
  const double *grad_u, const double *xi, const double *grad_xi, const double *grad_ln_n, const double *theta,
  const double *dtb, const double *em, double* GKYL_RESTRICT ho, double* GKYL_RESTRICT wn); 
GKYL_CU_DH void pkpm_harm_mat_3x_ser_p1(const double *dxv, const double *jacob_pos, const double *ho,
  double* GKYL_RESTRICT mat); 
GKYL_CU_DH void pkpm_harm_vol_3x1v_ser_p1(const double *w, const double *dxv, const double *mat, const double *f,
  double* GKYL_RESTRICT out); 

EXTERN_C_END 
