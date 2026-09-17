#pragma once 
#include <math.h> 
#include <gkyl_eqn_type.h> 
#include <gkyl_mat.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

// 1x3v Kernels
GKYL_CU_DH void mom_fpo_vlasov_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *a_i, const double *D_ij,  const double *f, double* GKYL_RESTRICT out);
GKYL_CU_DH void mom_bcorr_fpo_vlasov_1x3v_ser_p1(const double *w, const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *D_ij, const double *fIn, double* GKYL_RESTRICT out);
GKYL_CU_DH void mom_fpo_vlasov_coeff_correct_mat_1x3v_ser_p1(struct gkyl_mat *lhs, struct gkyl_mat *rhs, const double *fpo_moms, const double *boundary_corrections, const double *moms);
GKYL_CU_DH void mom_fpo_vlasov_coeff_correct_accum_1x3v_ser_p1(const double *drag_diff_coeff_corrs, double *drag_coeff, double *drag_coeff_surf, double *diff_coeff, double *diff_coeff_surf);

// 1x3v Kernels
GKYL_CU_DH void mom_fpo_vlasov_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *a_i, const double *D_ij,  const double *f, double* GKYL_RESTRICT out);
GKYL_CU_DH void mom_bcorr_fpo_vlasov_1x3v_ser_p2(const double *w, const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *D_ij, const double *fIn, double* GKYL_RESTRICT out);
GKYL_CU_DH void mom_fpo_vlasov_coeff_correct_mat_1x3v_ser_p2(struct gkyl_mat *lhs, struct gkyl_mat *rhs, const double *fpo_moms, const double *boundary_corrections, const double *moms);
GKYL_CU_DH void mom_fpo_vlasov_coeff_correct_accum_1x3v_ser_p2(const double *drag_diff_coeff_corrs, double *drag_coeff, double *drag_coeff_surf, double *diff_coeff, double *diff_coeff_surf);

// 2x3v Kernels
GKYL_CU_DH void mom_fpo_vlasov_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *a_i, const double *D_ij,  const double *f, double* GKYL_RESTRICT out);
GKYL_CU_DH void mom_bcorr_fpo_vlasov_2x3v_ser_p1(const double *w, const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *D_ij, const double *fIn, double* GKYL_RESTRICT out);
GKYL_CU_DH void mom_fpo_vlasov_coeff_correct_mat_2x3v_ser_p1(struct gkyl_mat *lhs, struct gkyl_mat *rhs, const double *fpo_moms, const double *boundary_corrections, const double *moms);
GKYL_CU_DH void mom_fpo_vlasov_coeff_correct_accum_2x3v_ser_p1(const double *drag_diff_coeff_corrs, double *drag_coeff, double *drag_coeff_surf, double *diff_coeff, double *diff_coeff_surf);

EXTERN_C_END 
