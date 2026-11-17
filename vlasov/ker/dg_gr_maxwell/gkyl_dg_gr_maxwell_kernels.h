#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
#include <gkyl_maxwell_kernels.h> 
typedef struct { double chi, gamma, K_phi, K_psi; } gkyl_dg_gr_maxwell_inp; 
EXTERN_C_BEG 

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_1x_ser_p1(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_1x_ser_p1(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_1x_ser_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gr_maxwell_geom_source_1x_ser_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, 
      const double *geom_factor_con_nodal, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_x_1x_ser_p1(const double *dxv, const int theta_pole,
      const double *jacob_pos_l, const double *jacob_pos_r,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_1x_ser_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *jacob_pos_l, const double *jacob_pos_r,
      const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_1x_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_lorentz_conf_1x1v_ser_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, const double *vierb_cov_nodal,
      const double *vierb_con_nodal, double* GKYL_RESTRICT E_conf_for_force, double* GKYL_RESTRICT B_conf_for_force); 
GKYL_CU_DH void gr_maxwell_current_dep_1x_1v_ser_p1(const double q_over_eps0, const double *lapse_nodal, const double *shift_nodal,
      const double *vierb_con_nodal, const double *m0, const double *m1i, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_gr_maxwell_lorentz_conf_1x2v_ser_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, const double *vierb_cov_nodal,
      const double *vierb_con_nodal, double* GKYL_RESTRICT E_conf_for_force, double* GKYL_RESTRICT B_conf_for_force); 
GKYL_CU_DH void gr_maxwell_current_dep_1x_2v_ser_p1(const double q_over_eps0, const double *lapse_nodal, const double *shift_nodal,
      const double *vierb_con_nodal, const double *m0, const double *m1i, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_gr_maxwell_lorentz_conf_1x3v_ser_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, const double *vierb_cov_nodal,
      const double *vierb_con_nodal, double* GKYL_RESTRICT E_conf_for_force, double* GKYL_RESTRICT B_conf_for_force); 
GKYL_CU_DH void gr_maxwell_current_dep_1x_3v_ser_p1(const double q_over_eps0, const double *lapse_nodal, const double *shift_nodal,
      const double *vierb_con_nodal, const double *m0, const double *m1i, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_gr_maxwell_divide_Jc_1x_ser_p2(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_1x_ser_p2(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_1x_ser_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gr_maxwell_geom_source_1x_ser_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, 
      const double *geom_factor_con_nodal, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_x_1x_ser_p2(const double *dxv, const int theta_pole,
      const double *jacob_pos_l, const double *jacob_pos_r,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_1x_ser_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *jacob_pos_l, const double *jacob_pos_r,
      const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_1x_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_lorentz_conf_1x1v_ser_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, const double *vierb_cov_nodal,
      const double *vierb_con_nodal, double* GKYL_RESTRICT E_conf_for_force, double* GKYL_RESTRICT B_conf_for_force); 
GKYL_CU_DH void gr_maxwell_current_dep_1x_1v_ser_p2(const double q_over_eps0, const double *lapse_nodal, const double *shift_nodal,
      const double *vierb_con_nodal, const double *m0, const double *m1i, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_gr_maxwell_lorentz_conf_1x2v_ser_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, const double *vierb_cov_nodal,
      const double *vierb_con_nodal, double* GKYL_RESTRICT E_conf_for_force, double* GKYL_RESTRICT B_conf_for_force); 
GKYL_CU_DH void gr_maxwell_current_dep_1x_2v_ser_p2(const double q_over_eps0, const double *lapse_nodal, const double *shift_nodal,
      const double *vierb_con_nodal, const double *m0, const double *m1i, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_gr_maxwell_lorentz_conf_1x3v_ser_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, const double *vierb_cov_nodal,
      const double *vierb_con_nodal, double* GKYL_RESTRICT E_conf_for_force, double* GKYL_RESTRICT B_conf_for_force); 
GKYL_CU_DH void gr_maxwell_current_dep_1x_3v_ser_p2(const double q_over_eps0, const double *lapse_nodal, const double *shift_nodal,
      const double *vierb_con_nodal, const double *m0, const double *m1i, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_gr_maxwell_divide_Jc_2x_ser_p1(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_2x_ser_p1(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_2x_ser_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gr_maxwell_geom_source_2x_ser_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, 
      const double *geom_factor_con_nodal, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_x_2x_ser_p1(const double *dxv, const int theta_pole,
      const double *jacob_pos_l, const double *jacob_pos_r,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_2x_ser_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *jacob_pos_l, const double *jacob_pos_r,
      const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_2x_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_y_2x_ser_p1(const double *dxv, const int theta_pole,
      const double *jacob_pos_l, const double *jacob_pos_r,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_y_2x_ser_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *jacob_pos_l, const double *jacob_pos_r,
      const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfy_2x_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_lorentz_conf_2x2v_ser_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, const double *vierb_cov_nodal,
      const double *vierb_con_nodal, double* GKYL_RESTRICT E_conf_for_force, double* GKYL_RESTRICT B_conf_for_force); 
GKYL_CU_DH void gr_maxwell_current_dep_2x_2v_ser_p1(const double q_over_eps0, const double *lapse_nodal, const double *shift_nodal,
      const double *vierb_con_nodal, const double *m0, const double *m1i, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_gr_maxwell_lorentz_conf_2x3v_ser_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, const double *vierb_cov_nodal,
      const double *vierb_con_nodal, double* GKYL_RESTRICT E_conf_for_force, double* GKYL_RESTRICT B_conf_for_force); 
GKYL_CU_DH void gr_maxwell_current_dep_2x_3v_ser_p1(const double q_over_eps0, const double *lapse_nodal, const double *shift_nodal,
      const double *vierb_con_nodal, const double *m0, const double *m1i, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_gr_maxwell_divide_Jc_2x_ser_p2(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_2x_ser_p2(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_2x_ser_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gr_maxwell_geom_source_2x_ser_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, 
      const double *geom_factor_con_nodal, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_x_2x_ser_p2(const double *dxv, const int theta_pole,
      const double *jacob_pos_l, const double *jacob_pos_r,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_2x_ser_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *jacob_pos_l, const double *jacob_pos_r,
      const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_2x_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_y_2x_ser_p2(const double *dxv, const int theta_pole,
      const double *jacob_pos_l, const double *jacob_pos_r,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_y_2x_ser_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *jacob_pos_l, const double *jacob_pos_r,
      const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfy_2x_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_lorentz_conf_2x2v_ser_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, const double *vierb_cov_nodal,
      const double *vierb_con_nodal, double* GKYL_RESTRICT E_conf_for_force, double* GKYL_RESTRICT B_conf_for_force); 
GKYL_CU_DH void gr_maxwell_current_dep_2x_2v_ser_p2(const double q_over_eps0, const double *lapse_nodal, const double *shift_nodal,
      const double *vierb_con_nodal, const double *m0, const double *m1i, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_gr_maxwell_lorentz_conf_2x3v_ser_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, const double *vierb_cov_nodal,
      const double *vierb_con_nodal, double* GKYL_RESTRICT E_conf_for_force, double* GKYL_RESTRICT B_conf_for_force); 
GKYL_CU_DH void gr_maxwell_current_dep_2x_3v_ser_p2(const double q_over_eps0, const double *lapse_nodal, const double *shift_nodal,
      const double *vierb_con_nodal, const double *m0, const double *m1i, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_gr_maxwell_divide_Jc_3x_ser_p1(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_3x_ser_p1(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_3x_ser_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gr_maxwell_geom_source_3x_ser_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, 
      const double *geom_factor_con_nodal, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_x_3x_ser_p1(const double *dxv, const int theta_pole,
      const double *jacob_pos_l, const double *jacob_pos_r,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_3x_ser_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *jacob_pos_l, const double *jacob_pos_r,
      const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_3x_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_y_3x_ser_p1(const double *dxv, const int theta_pole,
      const double *jacob_pos_l, const double *jacob_pos_r,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_y_3x_ser_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *jacob_pos_l, const double *jacob_pos_r,
      const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfy_3x_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_z_3x_ser_p1(const double *dxv, const int theta_pole,
      const double *jacob_pos_l, const double *jacob_pos_r,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_z_3x_ser_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *jacob_pos_l, const double *jacob_pos_r,
      const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfz_3x_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_lorentz_conf_3x3v_ser_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, const double *vierb_cov_nodal,
      const double *vierb_con_nodal, double* GKYL_RESTRICT E_conf_for_force, double* GKYL_RESTRICT B_conf_for_force); 
GKYL_CU_DH void gr_maxwell_current_dep_3x_3v_ser_p1(const double q_over_eps0, const double *lapse_nodal, const double *shift_nodal,
      const double *vierb_con_nodal, const double *m0, const double *m1i, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_gr_maxwell_divide_Jc_3x_ser_p2(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_3x_ser_p2(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_3x_ser_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gr_maxwell_geom_source_3x_ser_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, 
      const double *geom_factor_con_nodal, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_x_3x_ser_p2(const double *dxv, const int theta_pole,
      const double *jacob_pos_l, const double *jacob_pos_r,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_3x_ser_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *jacob_pos_l, const double *jacob_pos_r,
      const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_3x_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_y_3x_ser_p2(const double *dxv, const int theta_pole,
      const double *jacob_pos_l, const double *jacob_pos_r,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_y_3x_ser_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *jacob_pos_l, const double *jacob_pos_r,
      const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfy_3x_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_z_3x_ser_p2(const double *dxv, const int theta_pole,
      const double *jacob_pos_l, const double *jacob_pos_r,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_z_3x_ser_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *jacob_pos_l, const double *jacob_pos_r,
      const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfz_3x_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_lorentz_conf_3x3v_ser_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, const double *vierb_cov_nodal,
      const double *vierb_con_nodal, double* GKYL_RESTRICT E_conf_for_force, double* GKYL_RESTRICT B_conf_for_force); 
GKYL_CU_DH void gr_maxwell_current_dep_3x_3v_ser_p2(const double q_over_eps0, const double *lapse_nodal, const double *shift_nodal,
      const double *vierb_con_nodal, const double *m0, const double *m1i, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_gr_maxwell_divide_Jc_1x_tensor_p1(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_1x_tensor_p1(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_1x_tensor_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gr_maxwell_geom_source_1x_tensor_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, 
      const double *geom_factor_con_nodal, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_x_1x_tensor_p1(const double *dxv, const int theta_pole,
      const double *jacob_pos_l, const double *jacob_pos_r,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_1x_tensor_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *jacob_pos_l, const double *jacob_pos_r,
      const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_1x_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_lorentz_conf_1x1v_tensor_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, const double *vierb_cov_nodal,
      const double *vierb_con_nodal, double* GKYL_RESTRICT E_conf_for_force, double* GKYL_RESTRICT B_conf_for_force); 
GKYL_CU_DH void gr_maxwell_current_dep_1x_1v_tensor_p1(const double q_over_eps0, const double *lapse_nodal, const double *shift_nodal,
      const double *vierb_con_nodal, const double *m0, const double *m1i, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_gr_maxwell_lorentz_conf_1x2v_tensor_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, const double *vierb_cov_nodal,
      const double *vierb_con_nodal, double* GKYL_RESTRICT E_conf_for_force, double* GKYL_RESTRICT B_conf_for_force); 
GKYL_CU_DH void gr_maxwell_current_dep_1x_2v_tensor_p1(const double q_over_eps0, const double *lapse_nodal, const double *shift_nodal,
      const double *vierb_con_nodal, const double *m0, const double *m1i, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_gr_maxwell_lorentz_conf_1x3v_tensor_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, const double *vierb_cov_nodal,
      const double *vierb_con_nodal, double* GKYL_RESTRICT E_conf_for_force, double* GKYL_RESTRICT B_conf_for_force); 
GKYL_CU_DH void gr_maxwell_current_dep_1x_3v_tensor_p1(const double q_over_eps0, const double *lapse_nodal, const double *shift_nodal,
      const double *vierb_con_nodal, const double *m0, const double *m1i, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_gr_maxwell_divide_Jc_1x_tensor_p2(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_1x_tensor_p2(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_1x_tensor_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gr_maxwell_geom_source_1x_tensor_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, 
      const double *geom_factor_con_nodal, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_x_1x_tensor_p2(const double *dxv, const int theta_pole,
      const double *jacob_pos_l, const double *jacob_pos_r,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_1x_tensor_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *jacob_pos_l, const double *jacob_pos_r,
      const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_1x_tensor_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_lorentz_conf_1x1v_tensor_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, const double *vierb_cov_nodal,
      const double *vierb_con_nodal, double* GKYL_RESTRICT E_conf_for_force, double* GKYL_RESTRICT B_conf_for_force); 
GKYL_CU_DH void gr_maxwell_current_dep_1x_1v_tensor_p2(const double q_over_eps0, const double *lapse_nodal, const double *shift_nodal,
      const double *vierb_con_nodal, const double *m0, const double *m1i, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_gr_maxwell_lorentz_conf_1x2v_tensor_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, const double *vierb_cov_nodal,
      const double *vierb_con_nodal, double* GKYL_RESTRICT E_conf_for_force, double* GKYL_RESTRICT B_conf_for_force); 
GKYL_CU_DH void gr_maxwell_current_dep_1x_2v_tensor_p2(const double q_over_eps0, const double *lapse_nodal, const double *shift_nodal,
      const double *vierb_con_nodal, const double *m0, const double *m1i, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_gr_maxwell_lorentz_conf_1x3v_tensor_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, const double *vierb_cov_nodal,
      const double *vierb_con_nodal, double* GKYL_RESTRICT E_conf_for_force, double* GKYL_RESTRICT B_conf_for_force); 
GKYL_CU_DH void gr_maxwell_current_dep_1x_3v_tensor_p2(const double q_over_eps0, const double *lapse_nodal, const double *shift_nodal,
      const double *vierb_con_nodal, const double *m0, const double *m1i, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_gr_maxwell_divide_Jc_2x_tensor_p1(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_2x_tensor_p1(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_2x_tensor_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gr_maxwell_geom_source_2x_tensor_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, 
      const double *geom_factor_con_nodal, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_x_2x_tensor_p1(const double *dxv, const int theta_pole,
      const double *jacob_pos_l, const double *jacob_pos_r,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_2x_tensor_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *jacob_pos_l, const double *jacob_pos_r,
      const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_2x_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_y_2x_tensor_p1(const double *dxv, const int theta_pole,
      const double *jacob_pos_l, const double *jacob_pos_r,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_y_2x_tensor_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *jacob_pos_l, const double *jacob_pos_r,
      const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfy_2x_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_lorentz_conf_2x2v_tensor_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, const double *vierb_cov_nodal,
      const double *vierb_con_nodal, double* GKYL_RESTRICT E_conf_for_force, double* GKYL_RESTRICT B_conf_for_force); 
GKYL_CU_DH void gr_maxwell_current_dep_2x_2v_tensor_p1(const double q_over_eps0, const double *lapse_nodal, const double *shift_nodal,
      const double *vierb_con_nodal, const double *m0, const double *m1i, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_gr_maxwell_lorentz_conf_2x3v_tensor_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, const double *vierb_cov_nodal,
      const double *vierb_con_nodal, double* GKYL_RESTRICT E_conf_for_force, double* GKYL_RESTRICT B_conf_for_force); 
GKYL_CU_DH void gr_maxwell_current_dep_2x_3v_tensor_p1(const double q_over_eps0, const double *lapse_nodal, const double *shift_nodal,
      const double *vierb_con_nodal, const double *m0, const double *m1i, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_gr_maxwell_divide_Jc_2x_tensor_p2(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_2x_tensor_p2(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_2x_tensor_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gr_maxwell_geom_source_2x_tensor_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, 
      const double *geom_factor_con_nodal, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_x_2x_tensor_p2(const double *dxv, const int theta_pole,
      const double *jacob_pos_l, const double *jacob_pos_r,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_2x_tensor_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *jacob_pos_l, const double *jacob_pos_r,
      const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_2x_tensor_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_y_2x_tensor_p2(const double *dxv, const int theta_pole,
      const double *jacob_pos_l, const double *jacob_pos_r,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_y_2x_tensor_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *jacob_pos_l, const double *jacob_pos_r,
      const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfy_2x_tensor_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_lorentz_conf_2x2v_tensor_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, const double *vierb_cov_nodal,
      const double *vierb_con_nodal, double* GKYL_RESTRICT E_conf_for_force, double* GKYL_RESTRICT B_conf_for_force); 
GKYL_CU_DH void gr_maxwell_current_dep_2x_2v_tensor_p2(const double q_over_eps0, const double *lapse_nodal, const double *shift_nodal,
      const double *vierb_con_nodal, const double *m0, const double *m1i, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_gr_maxwell_lorentz_conf_2x3v_tensor_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, const double *vierb_cov_nodal,
      const double *vierb_con_nodal, double* GKYL_RESTRICT E_conf_for_force, double* GKYL_RESTRICT B_conf_for_force); 
GKYL_CU_DH void gr_maxwell_current_dep_2x_3v_tensor_p2(const double q_over_eps0, const double *lapse_nodal, const double *shift_nodal,
      const double *vierb_con_nodal, const double *m0, const double *m1i, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_gr_maxwell_divide_Jc_3x_tensor_p1(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_3x_tensor_p1(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_3x_tensor_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gr_maxwell_geom_source_3x_tensor_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, 
      const double *geom_factor_con_nodal, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_x_3x_tensor_p1(const double *dxv, const int theta_pole,
      const double *jacob_pos_l, const double *jacob_pos_r,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_3x_tensor_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *jacob_pos_l, const double *jacob_pos_r,
      const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_3x_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_y_3x_tensor_p1(const double *dxv, const int theta_pole,
      const double *jacob_pos_l, const double *jacob_pos_r,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_y_3x_tensor_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *jacob_pos_l, const double *jacob_pos_r,
      const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfy_3x_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_z_3x_tensor_p1(const double *dxv, const int theta_pole,
      const double *jacob_pos_l, const double *jacob_pos_r,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_z_3x_tensor_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *jacob_pos_l, const double *jacob_pos_r,
      const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfz_3x_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_lorentz_conf_3x3v_tensor_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, const double *vierb_cov_nodal,
      const double *vierb_con_nodal, double* GKYL_RESTRICT E_conf_for_force, double* GKYL_RESTRICT B_conf_for_force); 
GKYL_CU_DH void gr_maxwell_current_dep_3x_3v_tensor_p1(const double q_over_eps0, const double *lapse_nodal, const double *shift_nodal,
      const double *vierb_con_nodal, const double *m0, const double *m1i, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_gr_maxwell_divide_Jc_3x_tensor_p2(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_3x_tensor_p2(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_3x_tensor_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gr_maxwell_geom_source_3x_tensor_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, 
      const double *geom_factor_con_nodal, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_x_3x_tensor_p2(const double *dxv, const int theta_pole,
      const double *jacob_pos_l, const double *jacob_pos_r,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_3x_tensor_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *jacob_pos_l, const double *jacob_pos_r,
      const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_3x_tensor_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_y_3x_tensor_p2(const double *dxv, const int theta_pole,
      const double *jacob_pos_l, const double *jacob_pos_r,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_y_3x_tensor_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *jacob_pos_l, const double *jacob_pos_r,
      const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfy_3x_tensor_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_z_3x_tensor_p2(const double *dxv, const int theta_pole,
      const double *jacob_pos_l, const double *jacob_pos_r,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_z_3x_tensor_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *jacob_pos_l, const double *jacob_pos_r,
      const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfz_3x_tensor_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_lorentz_conf_3x3v_tensor_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, const double *vierb_cov_nodal,
      const double *vierb_con_nodal, double* GKYL_RESTRICT E_conf_for_force, double* GKYL_RESTRICT B_conf_for_force); 
GKYL_CU_DH void gr_maxwell_current_dep_3x_3v_tensor_p2(const double q_over_eps0, const double *lapse_nodal, const double *shift_nodal,
      const double *vierb_con_nodal, const double *m0, const double *m1i, double* GKYL_RESTRICT out); 
EXTERN_C_END 
