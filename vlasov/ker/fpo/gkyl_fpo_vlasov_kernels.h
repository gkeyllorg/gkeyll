#pragma once
#include <math.h>
#include <gkyl_util.h>

EXTERN_C_BEG

// 1x3v p1 kernels
GKYL_CU_DH void fpo_drag_coeff_1x3v_vx_ser_p1_lovx(const double *dxv, const double *gamma, const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff, double *drag_coeff_surf);
GKYL_CU_DH void fpo_drag_coeff_1x3v_vx_ser_p1_invx(const double *dxv, const double *gamma, const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff, double *drag_coeff_surf);
GKYL_CU_DH void fpo_drag_coeff_1x3v_vx_ser_p1_upvx(const double *dxv, const double *gamma, const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff, double *drag_coeff_surf);
GKYL_CU_DH void fpo_drag_coeff_1x3v_vy_ser_p1_lovy(const double *dxv, const double *gamma, const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff, double *drag_coeff_surf);
GKYL_CU_DH void fpo_drag_coeff_1x3v_vy_ser_p1_invy(const double *dxv, const double *gamma, const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff, double *drag_coeff_surf);
GKYL_CU_DH void fpo_drag_coeff_1x3v_vy_ser_p1_upvy(const double *dxv, const double *gamma, const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff, double *drag_coeff_surf);
GKYL_CU_DH void fpo_drag_coeff_1x3v_vz_ser_p1_lovz(const double *dxv, const double *gamma, const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff, double *drag_coeff_surf);
GKYL_CU_DH void fpo_drag_coeff_1x3v_vz_ser_p1_invz(const double *dxv, const double *gamma, const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff, double *drag_coeff_surf);
GKYL_CU_DH void fpo_drag_coeff_1x3v_vz_ser_p1_upvz(const double *dxv, const double *gamma, const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff, double *drag_coeff_surf);

GKYL_CU_DH void fpo_sgn_drag_coeff_1x3v_vx_ser_p1_lovx(const double *drag_coeff_surf, double *sgn_drag_coeff_surf, int *const_sgn_drag_coeff_surf);
GKYL_CU_DH void fpo_sgn_drag_coeff_1x3v_vx_ser_p1_invx(const double *drag_coeff_surf, double *sgn_drag_coeff_surf, int *const_sgn_drag_coeff_surf);
GKYL_CU_DH void fpo_sgn_drag_coeff_1x3v_vx_ser_p1_upvx(const double *drag_coeff_surf, double *sgn_drag_coeff_surf, int *const_sgn_drag_coeff_surf);
GKYL_CU_DH void fpo_sgn_drag_coeff_1x3v_vy_ser_p1_lovy(const double *drag_coeff_surf, double *sgn_drag_coeff_surf, int *const_sgn_drag_coeff_surf);
GKYL_CU_DH void fpo_sgn_drag_coeff_1x3v_vy_ser_p1_invy(const double *drag_coeff_surf, double *sgn_drag_coeff_surf, int *const_sgn_drag_coeff_surf);
GKYL_CU_DH void fpo_sgn_drag_coeff_1x3v_vy_ser_p1_upvy(const double *drag_coeff_surf, double *sgn_drag_coeff_surf, int *const_sgn_drag_coeff_surf);
GKYL_CU_DH void fpo_sgn_drag_coeff_1x3v_vz_ser_p1_lovz(const double *drag_coeff_surf, double *sgn_drag_coeff_surf, int *const_sgn_drag_coeff_surf);
GKYL_CU_DH void fpo_sgn_drag_coeff_1x3v_vz_ser_p1_invz(const double *drag_coeff_surf, double *sgn_drag_coeff_surf, int *const_sgn_drag_coeff_surf);
GKYL_CU_DH void fpo_sgn_drag_coeff_1x3v_vz_ser_p1_upvz(const double *drag_coeff_surf, double *sgn_drag_coeff_surf, int *const_sgn_drag_coeff_surf);

GKYL_CU_DH void fpo_diff_coeff_diag_1x3v_vx_ser_p1_lovx(const double *dxv, const double *gamma, const double* fpo_g_stencil[3], const double* fpo_d2gdv2_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_diag_1x3v_vx_ser_p1_invx(const double *dxv, const double *gamma, const double* fpo_g_stencil[3], const double* fpo_d2gdv2_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_diag_1x3v_vx_ser_p1_upvx(const double *dxv, const double *gamma, const double* fpo_g_stencil[3], const double* fpo_d2gdv2_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_diag_1x3v_vy_ser_p1_lovy(const double *dxv, const double *gamma, const double* fpo_g_stencil[3], const double* fpo_d2gdv2_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_diag_1x3v_vy_ser_p1_invy(const double *dxv, const double *gamma, const double* fpo_g_stencil[3], const double* fpo_d2gdv2_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_diag_1x3v_vy_ser_p1_upvy(const double *dxv, const double *gamma, const double* fpo_g_stencil[3], const double* fpo_d2gdv2_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_diag_1x3v_vz_ser_p1_lovz(const double *dxv, const double *gamma, const double* fpo_g_stencil[3], const double* fpo_d2gdv2_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_diag_1x3v_vz_ser_p1_invz(const double *dxv, const double *gamma, const double* fpo_g_stencil[3], const double* fpo_d2gdv2_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_diag_1x3v_vz_ser_p1_upvz(const double *dxv, const double *gamma, const double* fpo_g_stencil[3], const double* fpo_d2gdv2_surf, double *diff_coeff);

GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvy_ser_p1_invx_invy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvy_ser_p1_invx_lovy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvy_ser_p1_invx_upvy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvy_ser_p1_lovx_invy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvy_ser_p1_lovx_lovy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvy_ser_p1_lovx_upvy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvy_ser_p1_upvx_invy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvy_ser_p1_upvx_lovy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvy_ser_p1_upvx_upvy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvz_ser_p1_invx_invz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvz_ser_p1_invx_lovz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvz_ser_p1_invx_upvz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvz_ser_p1_lovx_invz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvz_ser_p1_lovx_lovz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvz_ser_p1_lovx_upvz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvz_ser_p1_upvx_invz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvz_ser_p1_upvx_lovz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvz_ser_p1_upvx_upvz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvx_ser_p1_invy_invx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvx_ser_p1_invy_lovx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvx_ser_p1_invy_upvx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvx_ser_p1_lovy_invx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvx_ser_p1_lovy_lovx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvx_ser_p1_lovy_upvx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvx_ser_p1_upvy_invx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvx_ser_p1_upvy_lovx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvx_ser_p1_upvy_upvx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvz_ser_p1_invy_invz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvz_ser_p1_invy_lovz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvz_ser_p1_invy_upvz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvz_ser_p1_lovy_invz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvz_ser_p1_lovy_lovz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvz_ser_p1_lovy_upvz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvz_ser_p1_upvy_invz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvz_ser_p1_upvy_lovz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvz_ser_p1_upvy_upvz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvx_ser_p1_invz_invx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvx_ser_p1_invz_lovx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvx_ser_p1_invz_upvx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvx_ser_p1_lovz_invx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvx_ser_p1_lovz_lovx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvx_ser_p1_lovz_upvx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvx_ser_p1_upvz_invx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvx_ser_p1_upvz_lovx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvx_ser_p1_upvz_upvx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvy_ser_p1_invz_invy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvy_ser_p1_invz_lovy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvy_ser_p1_invz_upvy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvy_ser_p1_lovz_invy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvy_ser_p1_lovz_lovy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvy_ser_p1_lovz_upvy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvy_ser_p1_upvz_invy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvy_ser_p1_upvz_lovy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvy_ser_p1_upvz_upvy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);

GKYL_CU_DH void fpo_diff_coeff_surf_1x3v_vx_ser_p1(const double *diff_coeff_L, const double *diff_coeff_R, double *diff_coeff_surf_R);
GKYL_CU_DH void fpo_diff_coeff_surf_1x3v_vy_ser_p1(const double *diff_coeff_L, const double *diff_coeff_R, double *diff_coeff_surf_R);
GKYL_CU_DH void fpo_diff_coeff_surf_1x3v_vz_ser_p1(const double *diff_coeff_L, const double *diff_coeff_R, double *diff_coeff_surf_R);

GKYL_CU_DH double fpo_vlasov_drag_vol_1x3v_ser_p1(const double* dxv, const double* drag_coeff, const double* f, double* GKYL_RESTRICT out);

GKYL_CU_DH double fpo_vlasov_drag_surfvx_1x3v_ser_p1(const double* dxv, const double *alpha_surf_L, const double *alpha_surf_R, const double *sgn_alpha_surf_L, const double *sgn_alpha_surf_R, const int *const_sgn_alpha_L, const int *const_sgn_alpha_R, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_drag_surfvy_1x3v_ser_p1(const double* dxv, const double *alpha_surf_L, const double *alpha_surf_R, const double *sgn_alpha_surf_L, const double *sgn_alpha_surf_R, const int *const_sgn_alpha_L, const int *const_sgn_alpha_R, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_drag_surfvz_1x3v_ser_p1(const double* dxv, const double *alpha_surf_L, const double *alpha_surf_R, const double *sgn_alpha_surf_L, const double *sgn_alpha_surf_R, const int *const_sgn_alpha_L, const int *const_sgn_alpha_R, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out);

GKYL_CU_DH double fpo_vlasov_drag_boundary_surfvx_1x3v_ser_p1(const double* dxv, const double *alpha_surf_Edge, const double *alpha_surf_Skin, 
  const double *sgn_alpha_surf_Edge, const double *sgn_alpha_surf_Skin, const int *const_sgn_alpha_Edge, const int *const_sgn_alpha_Skin, const int Edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double fpo_vlasov_drag_boundary_surfvy_1x3v_ser_p1(const double* dxv, const double *alpha_surf_Edge, const double *alpha_surf_Skin, 
  const double *sgn_alpha_surf_Edge, const double *sgn_alpha_surf_Skin, const int *const_sgn_alpha_Edge, const int *const_sgn_alpha_Skin, const int Edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double fpo_vlasov_drag_boundary_surfvz_1x3v_ser_p1(const double* dxv, const double *alpha_surf_Edge, const double *alpha_surf_Skin, 
  const double *sgn_alpha_surf_Edge, const double *sgn_alpha_surf_Skin, const int *const_sgn_alpha_Edge, const int *const_sgn_alpha_Skin, const int Edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double fpo_vlasov_diff_vol_1x3v_ser_p1(const double* dxv, const double *diff_coeff, const double* f, double* GKYL_RESTRICT out);

GKYL_CU_DH double fpo_vlasov_diff_surfvxvx_1x3v_ser_p1(const double* dxv, const double* diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvxvy_1x3v_ser_p1(const double* dxv, const double* diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvxvz_1x3v_ser_p1(const double* dxv, const double* diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvyvx_1x3v_ser_p1(const double* dxv, const double* diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvyvy_1x3v_ser_p1(const double* dxv, const double* diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvyvz_1x3v_ser_p1(const double* dxv, const double* diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvzvx_1x3v_ser_p1(const double* dxv, const double* diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvzvy_1x3v_ser_p1(const double* dxv, const double* diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvzvz_1x3v_ser_p1(const double* dxv, const double* diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* GKYL_RESTRICT out);

GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvx_1x3v_ser_p1_lovx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvx_1x3v_ser_p1_upvx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvy_1x3v_ser_p1_lovy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvy_1x3v_ser_p1_upvy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvz_1x3v_ser_p1_lovz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvz_1x3v_ser_p1_upvz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);

GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvy_1x3v_ser_p1_invx_lovy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvy_1x3v_ser_p1_invx_upvy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvy_1x3v_ser_p1_lovx_invy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvy_1x3v_ser_p1_lovx_lovy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvy_1x3v_ser_p1_lovx_upvy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvy_1x3v_ser_p1_upvx_invy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvy_1x3v_ser_p1_upvx_lovy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvy_1x3v_ser_p1_upvx_upvy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvz_1x3v_ser_p1_invx_lovz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvz_1x3v_ser_p1_invx_upvz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvz_1x3v_ser_p1_lovx_invz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvz_1x3v_ser_p1_lovx_lovz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvz_1x3v_ser_p1_lovx_upvz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvz_1x3v_ser_p1_upvx_invz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvz_1x3v_ser_p1_upvx_lovz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvz_1x3v_ser_p1_upvx_upvz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvx_1x3v_ser_p1_invy_lovx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvx_1x3v_ser_p1_invy_upvx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvx_1x3v_ser_p1_lovy_invx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvx_1x3v_ser_p1_lovy_lovx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvx_1x3v_ser_p1_lovy_upvx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvx_1x3v_ser_p1_upvy_invx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvx_1x3v_ser_p1_upvy_lovx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvx_1x3v_ser_p1_upvy_upvx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvz_1x3v_ser_p1_invy_lovz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvz_1x3v_ser_p1_invy_upvz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvz_1x3v_ser_p1_lovy_invz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvz_1x3v_ser_p1_lovy_lovz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvz_1x3v_ser_p1_lovy_upvz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvz_1x3v_ser_p1_upvy_invz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvz_1x3v_ser_p1_upvy_lovz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvz_1x3v_ser_p1_upvy_upvz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvx_1x3v_ser_p1_invz_lovx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvx_1x3v_ser_p1_invz_upvx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvx_1x3v_ser_p1_lovz_invx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvx_1x3v_ser_p1_lovz_lovx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvx_1x3v_ser_p1_lovz_upvx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvx_1x3v_ser_p1_upvz_invx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvx_1x3v_ser_p1_upvz_lovx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvx_1x3v_ser_p1_upvz_upvx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvy_1x3v_ser_p1_invz_lovy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvy_1x3v_ser_p1_invz_upvy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvy_1x3v_ser_p1_lovz_invy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvy_1x3v_ser_p1_lovz_lovy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvy_1x3v_ser_p1_lovz_upvy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvy_1x3v_ser_p1_upvz_invy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvy_1x3v_ser_p1_upvz_lovy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvy_1x3v_ser_p1_upvz_upvy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);


// 1x3v p2 kernels
GKYL_CU_DH void fpo_drag_coeff_1x3v_vx_ser_p2_lovx(const double *dxv, const double *gamma, const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff, double *drag_coeff_surf);
GKYL_CU_DH void fpo_drag_coeff_1x3v_vx_ser_p2_invx(const double *dxv, const double *gamma, const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff, double *drag_coeff_surf);
GKYL_CU_DH void fpo_drag_coeff_1x3v_vx_ser_p2_upvx(const double *dxv, const double *gamma, const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff, double *drag_coeff_surf);
GKYL_CU_DH void fpo_drag_coeff_1x3v_vy_ser_p2_lovy(const double *dxv, const double *gamma, const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff, double *drag_coeff_surf);
GKYL_CU_DH void fpo_drag_coeff_1x3v_vy_ser_p2_invy(const double *dxv, const double *gamma, const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff, double *drag_coeff_surf);
GKYL_CU_DH void fpo_drag_coeff_1x3v_vy_ser_p2_upvy(const double *dxv, const double *gamma, const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff, double *drag_coeff_surf);
GKYL_CU_DH void fpo_drag_coeff_1x3v_vz_ser_p2_lovz(const double *dxv, const double *gamma, const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff, double *drag_coeff_surf);
GKYL_CU_DH void fpo_drag_coeff_1x3v_vz_ser_p2_invz(const double *dxv, const double *gamma, const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff, double *drag_coeff_surf);
GKYL_CU_DH void fpo_drag_coeff_1x3v_vz_ser_p2_upvz(const double *dxv, const double *gamma, const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff, double *drag_coeff_surf);

GKYL_CU_DH void fpo_sgn_drag_coeff_1x3v_vx_ser_p2_lovx(const double *drag_coeff_surf, double *sgn_drag_coeff_surf, int *const_sgn_drag_coeff_surf);
GKYL_CU_DH void fpo_sgn_drag_coeff_1x3v_vx_ser_p2_invx(const double *drag_coeff_surf, double *sgn_drag_coeff_surf, int *const_sgn_drag_coeff_surf);
GKYL_CU_DH void fpo_sgn_drag_coeff_1x3v_vx_ser_p2_upvx(const double *drag_coeff_surf, double *sgn_drag_coeff_surf, int *const_sgn_drag_coeff_surf);
GKYL_CU_DH void fpo_sgn_drag_coeff_1x3v_vy_ser_p2_lovy(const double *drag_coeff_surf, double *sgn_drag_coeff_surf, int *const_sgn_drag_coeff_surf);
GKYL_CU_DH void fpo_sgn_drag_coeff_1x3v_vy_ser_p2_invy(const double *drag_coeff_surf, double *sgn_drag_coeff_surf, int *const_sgn_drag_coeff_surf);
GKYL_CU_DH void fpo_sgn_drag_coeff_1x3v_vy_ser_p2_upvy(const double *drag_coeff_surf, double *sgn_drag_coeff_surf, int *const_sgn_drag_coeff_surf);
GKYL_CU_DH void fpo_sgn_drag_coeff_1x3v_vz_ser_p2_lovz(const double *drag_coeff_surf, double *sgn_drag_coeff_surf, int *const_sgn_drag_coeff_surf);
GKYL_CU_DH void fpo_sgn_drag_coeff_1x3v_vz_ser_p2_invz(const double *drag_coeff_surf, double *sgn_drag_coeff_surf, int *const_sgn_drag_coeff_surf);
GKYL_CU_DH void fpo_sgn_drag_coeff_1x3v_vz_ser_p2_upvz(const double *drag_coeff_surf, double *sgn_drag_coeff_surf, int *const_sgn_drag_coeff_surf);

GKYL_CU_DH void fpo_diff_coeff_diag_1x3v_vx_ser_p2_lovx(const double *dxv, const double *gamma, const double* fpo_g_stencil[3], const double* fpo_d2gdv2_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_diag_1x3v_vx_ser_p2_invx(const double *dxv, const double *gamma, const double* fpo_g_stencil[3], const double* fpo_d2gdv2_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_diag_1x3v_vx_ser_p2_upvx(const double *dxv, const double *gamma, const double* fpo_g_stencil[3], const double* fpo_d2gdv2_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_diag_1x3v_vy_ser_p2_lovy(const double *dxv, const double *gamma, const double* fpo_g_stencil[3], const double* fpo_d2gdv2_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_diag_1x3v_vy_ser_p2_invy(const double *dxv, const double *gamma, const double* fpo_g_stencil[3], const double* fpo_d2gdv2_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_diag_1x3v_vy_ser_p2_upvy(const double *dxv, const double *gamma, const double* fpo_g_stencil[3], const double* fpo_d2gdv2_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_diag_1x3v_vz_ser_p2_lovz(const double *dxv, const double *gamma, const double* fpo_g_stencil[3], const double* fpo_d2gdv2_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_diag_1x3v_vz_ser_p2_invz(const double *dxv, const double *gamma, const double* fpo_g_stencil[3], const double* fpo_d2gdv2_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_diag_1x3v_vz_ser_p2_upvz(const double *dxv, const double *gamma, const double* fpo_g_stencil[3], const double* fpo_d2gdv2_surf, double *diff_coeff);

GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvy_ser_p2_invx_invy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvy_ser_p2_invx_lovy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvy_ser_p2_invx_upvy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvy_ser_p2_lovx_invy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvy_ser_p2_lovx_lovy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvy_ser_p2_lovx_upvy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvy_ser_p2_upvx_invy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvy_ser_p2_upvx_lovy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvy_ser_p2_upvx_upvy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvz_ser_p2_invx_invz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvz_ser_p2_invx_lovz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvz_ser_p2_invx_upvz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvz_ser_p2_lovx_invz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvz_ser_p2_lovx_lovz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvz_ser_p2_lovx_upvz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvz_ser_p2_upvx_invz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvz_ser_p2_upvx_lovz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvz_ser_p2_upvx_upvz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvx_ser_p2_invy_invx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvx_ser_p2_invy_lovx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvx_ser_p2_invy_upvx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvx_ser_p2_lovy_invx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvx_ser_p2_lovy_lovx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvx_ser_p2_lovy_upvx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvx_ser_p2_upvy_invx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvx_ser_p2_upvy_lovx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvx_ser_p2_upvy_upvx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvz_ser_p2_invy_invz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvz_ser_p2_invy_lovz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvz_ser_p2_invy_upvz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvz_ser_p2_lovy_invz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvz_ser_p2_lovy_lovz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvz_ser_p2_lovy_upvz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvz_ser_p2_upvy_invz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvz_ser_p2_upvy_lovz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvz_ser_p2_upvy_upvz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvx_ser_p2_invz_invx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvx_ser_p2_invz_lovx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvx_ser_p2_invz_upvx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvx_ser_p2_lovz_invx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvx_ser_p2_lovz_lovx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvx_ser_p2_lovz_upvx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvx_ser_p2_upvz_invx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvx_ser_p2_upvz_lovx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvx_ser_p2_upvz_upvx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvy_ser_p2_invz_invy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvy_ser_p2_invz_lovy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvy_ser_p2_invz_upvy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvy_ser_p2_lovz_invy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvy_ser_p2_lovz_lovy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvy_ser_p2_lovz_upvy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvy_ser_p2_upvz_invy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvy_ser_p2_upvz_lovy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvy_ser_p2_upvz_upvy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);

GKYL_CU_DH void fpo_diff_coeff_surf_1x3v_vx_ser_p2(const double *diff_coeff_L, const double *diff_coeff_R, double *diff_coeff_surf_R);
GKYL_CU_DH void fpo_diff_coeff_surf_1x3v_vy_ser_p2(const double *diff_coeff_L, const double *diff_coeff_R, double *diff_coeff_surf_R);
GKYL_CU_DH void fpo_diff_coeff_surf_1x3v_vz_ser_p2(const double *diff_coeff_L, const double *diff_coeff_R, double *diff_coeff_surf_R);

GKYL_CU_DH double fpo_vlasov_drag_vol_1x3v_ser_p2(const double* dxv, const double* drag_coeff, const double* f, double* GKYL_RESTRICT out);

GKYL_CU_DH double fpo_vlasov_drag_surfvx_1x3v_ser_p2(const double* dxv, const double *alpha_surf_L, const double *alpha_surf_R, const double *sgn_alpha_surf_L, const double *sgn_alpha_surf_R, const int *const_sgn_alpha_L, const int *const_sgn_alpha_R, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_drag_surfvy_1x3v_ser_p2(const double* dxv, const double *alpha_surf_L, const double *alpha_surf_R, const double *sgn_alpha_surf_L, const double *sgn_alpha_surf_R, const int *const_sgn_alpha_L, const int *const_sgn_alpha_R, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_drag_surfvz_1x3v_ser_p2(const double* dxv, const double *alpha_surf_L, const double *alpha_surf_R, const double *sgn_alpha_surf_L, const double *sgn_alpha_surf_R, const int *const_sgn_alpha_L, const int *const_sgn_alpha_R, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out);

GKYL_CU_DH double fpo_vlasov_drag_boundary_surfvx_1x3v_ser_p2(const double* dxv, const double *alpha_surf_Edge, const double *alpha_surf_Skin, 
  const double *sgn_alpha_surf_Edge, const double *sgn_alpha_surf_Skin, const int *const_sgn_alpha_Edge, const int *const_sgn_alpha_Skin, const int Edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double fpo_vlasov_drag_boundary_surfvy_1x3v_ser_p2(const double* dxv, const double *alpha_surf_Edge, const double *alpha_surf_Skin, 
  const double *sgn_alpha_surf_Edge, const double *sgn_alpha_surf_Skin, const int *const_sgn_alpha_Edge, const int *const_sgn_alpha_Skin, const int Edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double fpo_vlasov_drag_boundary_surfvz_1x3v_ser_p2(const double* dxv, const double *alpha_surf_Edge, const double *alpha_surf_Skin, 
  const double *sgn_alpha_surf_Edge, const double *sgn_alpha_surf_Skin, const int *const_sgn_alpha_Edge, const int *const_sgn_alpha_Skin, const int Edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double fpo_vlasov_diff_vol_1x3v_ser_p2(const double* dxv, const double *diff_coeff, const double* f, double* GKYL_RESTRICT out);

GKYL_CU_DH double fpo_vlasov_diff_surfvxvx_1x3v_ser_p2(const double* dxv, const double* diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvxvy_1x3v_ser_p2(const double* dxv, const double* diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvxvz_1x3v_ser_p2(const double* dxv, const double* diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvyvx_1x3v_ser_p2(const double* dxv, const double* diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvyvy_1x3v_ser_p2(const double* dxv, const double* diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvyvz_1x3v_ser_p2(const double* dxv, const double* diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvzvx_1x3v_ser_p2(const double* dxv, const double* diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvzvy_1x3v_ser_p2(const double* dxv, const double* diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvzvz_1x3v_ser_p2(const double* dxv, const double* diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* GKYL_RESTRICT out);

GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvx_1x3v_ser_p2_lovx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvx_1x3v_ser_p2_upvx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvy_1x3v_ser_p2_lovy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvy_1x3v_ser_p2_upvy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvz_1x3v_ser_p2_lovz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvz_1x3v_ser_p2_upvz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);

GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvy_1x3v_ser_p2_invx_lovy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvy_1x3v_ser_p2_invx_upvy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvy_1x3v_ser_p2_lovx_invy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvy_1x3v_ser_p2_lovx_lovy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvy_1x3v_ser_p2_lovx_upvy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvy_1x3v_ser_p2_upvx_invy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvy_1x3v_ser_p2_upvx_lovy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvy_1x3v_ser_p2_upvx_upvy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvz_1x3v_ser_p2_invx_lovz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvz_1x3v_ser_p2_invx_upvz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvz_1x3v_ser_p2_lovx_invz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvz_1x3v_ser_p2_lovx_lovz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvz_1x3v_ser_p2_lovx_upvz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvz_1x3v_ser_p2_upvx_invz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvz_1x3v_ser_p2_upvx_lovz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvz_1x3v_ser_p2_upvx_upvz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvx_1x3v_ser_p2_invy_lovx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvx_1x3v_ser_p2_invy_upvx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvx_1x3v_ser_p2_lovy_invx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvx_1x3v_ser_p2_lovy_lovx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvx_1x3v_ser_p2_lovy_upvx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvx_1x3v_ser_p2_upvy_invx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvx_1x3v_ser_p2_upvy_lovx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvx_1x3v_ser_p2_upvy_upvx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvz_1x3v_ser_p2_invy_lovz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvz_1x3v_ser_p2_invy_upvz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvz_1x3v_ser_p2_lovy_invz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvz_1x3v_ser_p2_lovy_lovz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvz_1x3v_ser_p2_lovy_upvz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvz_1x3v_ser_p2_upvy_invz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvz_1x3v_ser_p2_upvy_lovz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvz_1x3v_ser_p2_upvy_upvz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvx_1x3v_ser_p2_invz_lovx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvx_1x3v_ser_p2_invz_upvx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvx_1x3v_ser_p2_lovz_invx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvx_1x3v_ser_p2_lovz_lovx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvx_1x3v_ser_p2_lovz_upvx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvx_1x3v_ser_p2_upvz_invx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvx_1x3v_ser_p2_upvz_lovx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvx_1x3v_ser_p2_upvz_upvx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvy_1x3v_ser_p2_invz_lovy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvy_1x3v_ser_p2_invz_upvy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvy_1x3v_ser_p2_lovz_invy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvy_1x3v_ser_p2_lovz_lovy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvy_1x3v_ser_p2_lovz_upvy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvy_1x3v_ser_p2_upvz_invy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvy_1x3v_ser_p2_upvz_lovy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvy_1x3v_ser_p2_upvz_upvy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);


// 2x3v p1 kernels
GKYL_CU_DH void fpo_drag_coeff_2x3v_vx_ser_p1_lovx(const double *dxv, const double *gamma, const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff, double *drag_coeff_surf);
GKYL_CU_DH void fpo_drag_coeff_2x3v_vx_ser_p1_invx(const double *dxv, const double *gamma, const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff, double *drag_coeff_surf);
GKYL_CU_DH void fpo_drag_coeff_2x3v_vx_ser_p1_upvx(const double *dxv, const double *gamma, const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff, double *drag_coeff_surf);
GKYL_CU_DH void fpo_drag_coeff_2x3v_vy_ser_p1_lovy(const double *dxv, const double *gamma, const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff, double *drag_coeff_surf);
GKYL_CU_DH void fpo_drag_coeff_2x3v_vy_ser_p1_invy(const double *dxv, const double *gamma, const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff, double *drag_coeff_surf);
GKYL_CU_DH void fpo_drag_coeff_2x3v_vy_ser_p1_upvy(const double *dxv, const double *gamma, const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff, double *drag_coeff_surf);
GKYL_CU_DH void fpo_drag_coeff_2x3v_vz_ser_p1_lovz(const double *dxv, const double *gamma, const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff, double *drag_coeff_surf);
GKYL_CU_DH void fpo_drag_coeff_2x3v_vz_ser_p1_invz(const double *dxv, const double *gamma, const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff, double *drag_coeff_surf);
GKYL_CU_DH void fpo_drag_coeff_2x3v_vz_ser_p1_upvz(const double *dxv, const double *gamma, const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff, double *drag_coeff_surf);

GKYL_CU_DH void fpo_sgn_drag_coeff_2x3v_vx_ser_p1_lovx(const double *drag_coeff_surf, double *sgn_drag_coeff_surf, int *const_sgn_drag_coeff_surf);
GKYL_CU_DH void fpo_sgn_drag_coeff_2x3v_vx_ser_p1_invx(const double *drag_coeff_surf, double *sgn_drag_coeff_surf, int *const_sgn_drag_coeff_surf);
GKYL_CU_DH void fpo_sgn_drag_coeff_2x3v_vx_ser_p1_upvx(const double *drag_coeff_surf, double *sgn_drag_coeff_surf, int *const_sgn_drag_coeff_surf);
GKYL_CU_DH void fpo_sgn_drag_coeff_2x3v_vy_ser_p1_lovy(const double *drag_coeff_surf, double *sgn_drag_coeff_surf, int *const_sgn_drag_coeff_surf);
GKYL_CU_DH void fpo_sgn_drag_coeff_2x3v_vy_ser_p1_invy(const double *drag_coeff_surf, double *sgn_drag_coeff_surf, int *const_sgn_drag_coeff_surf);
GKYL_CU_DH void fpo_sgn_drag_coeff_2x3v_vy_ser_p1_upvy(const double *drag_coeff_surf, double *sgn_drag_coeff_surf, int *const_sgn_drag_coeff_surf);
GKYL_CU_DH void fpo_sgn_drag_coeff_2x3v_vz_ser_p1_lovz(const double *drag_coeff_surf, double *sgn_drag_coeff_surf, int *const_sgn_drag_coeff_surf);
GKYL_CU_DH void fpo_sgn_drag_coeff_2x3v_vz_ser_p1_invz(const double *drag_coeff_surf, double *sgn_drag_coeff_surf, int *const_sgn_drag_coeff_surf);
GKYL_CU_DH void fpo_sgn_drag_coeff_2x3v_vz_ser_p1_upvz(const double *drag_coeff_surf, double *sgn_drag_coeff_surf, int *const_sgn_drag_coeff_surf);

GKYL_CU_DH void fpo_diff_coeff_diag_2x3v_vx_ser_p1_lovx(const double *dxv, const double *gamma, const double* fpo_g_stencil[3], const double* fpo_d2gdv2_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_diag_2x3v_vx_ser_p1_invx(const double *dxv, const double *gamma, const double* fpo_g_stencil[3], const double* fpo_d2gdv2_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_diag_2x3v_vx_ser_p1_upvx(const double *dxv, const double *gamma, const double* fpo_g_stencil[3], const double* fpo_d2gdv2_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_diag_2x3v_vy_ser_p1_lovy(const double *dxv, const double *gamma, const double* fpo_g_stencil[3], const double* fpo_d2gdv2_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_diag_2x3v_vy_ser_p1_invy(const double *dxv, const double *gamma, const double* fpo_g_stencil[3], const double* fpo_d2gdv2_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_diag_2x3v_vy_ser_p1_upvy(const double *dxv, const double *gamma, const double* fpo_g_stencil[3], const double* fpo_d2gdv2_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_diag_2x3v_vz_ser_p1_lovz(const double *dxv, const double *gamma, const double* fpo_g_stencil[3], const double* fpo_d2gdv2_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_diag_2x3v_vz_ser_p1_invz(const double *dxv, const double *gamma, const double* fpo_g_stencil[3], const double* fpo_d2gdv2_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_diag_2x3v_vz_ser_p1_upvz(const double *dxv, const double *gamma, const double* fpo_g_stencil[3], const double* fpo_d2gdv2_surf, double *diff_coeff);

GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vxvy_ser_p1_invx_invy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vxvy_ser_p1_invx_lovy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vxvy_ser_p1_invx_upvy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vxvy_ser_p1_lovx_invy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vxvy_ser_p1_lovx_lovy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vxvy_ser_p1_lovx_upvy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vxvy_ser_p1_upvx_invy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vxvy_ser_p1_upvx_lovy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vxvy_ser_p1_upvx_upvy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vxvz_ser_p1_invx_invz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vxvz_ser_p1_invx_lovz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vxvz_ser_p1_invx_upvz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vxvz_ser_p1_lovx_invz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vxvz_ser_p1_lovx_lovz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vxvz_ser_p1_lovx_upvz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vxvz_ser_p1_upvx_invz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vxvz_ser_p1_upvx_lovz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vxvz_ser_p1_upvx_upvz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vyvx_ser_p1_invy_invx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vyvx_ser_p1_invy_lovx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vyvx_ser_p1_invy_upvx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vyvx_ser_p1_lovy_invx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vyvx_ser_p1_lovy_lovx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vyvx_ser_p1_lovy_upvx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vyvx_ser_p1_upvy_invx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vyvx_ser_p1_upvy_lovx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vyvx_ser_p1_upvy_upvx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vyvz_ser_p1_invy_invz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vyvz_ser_p1_invy_lovz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vyvz_ser_p1_invy_upvz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vyvz_ser_p1_lovy_invz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vyvz_ser_p1_lovy_lovz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vyvz_ser_p1_lovy_upvz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vyvz_ser_p1_upvy_invz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vyvz_ser_p1_upvy_lovz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vyvz_ser_p1_upvy_upvz(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vzvx_ser_p1_invz_invx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vzvx_ser_p1_invz_lovx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vzvx_ser_p1_invz_upvx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vzvx_ser_p1_lovz_invx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vzvx_ser_p1_lovz_lovx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vzvx_ser_p1_lovz_upvx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vzvx_ser_p1_upvz_invx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vzvx_ser_p1_upvz_lovx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vzvx_ser_p1_upvz_upvx(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vzvy_ser_p1_invz_invy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vzvy_ser_p1_invz_lovy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vzvy_ser_p1_invz_upvy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vzvy_ser_p1_lovz_invy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vzvy_ser_p1_lovz_lovy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vzvy_ser_p1_lovz_upvy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vzvy_ser_p1_upvz_invy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vzvy_ser_p1_upvz_lovy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vzvy_ser_p1_upvz_upvy(const double *dxv, const double* gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff);

GKYL_CU_DH void fpo_diff_coeff_surf_2x3v_vx_ser_p1(const double *diff_coeff_L, const double *diff_coeff_R, double *diff_coeff_surf_R);
GKYL_CU_DH void fpo_diff_coeff_surf_2x3v_vy_ser_p1(const double *diff_coeff_L, const double *diff_coeff_R, double *diff_coeff_surf_R);
GKYL_CU_DH void fpo_diff_coeff_surf_2x3v_vz_ser_p1(const double *diff_coeff_L, const double *diff_coeff_R, double *diff_coeff_surf_R);

GKYL_CU_DH double fpo_vlasov_drag_vol_2x3v_ser_p1(const double* dxv, const double* drag_coeff, const double* f, double* GKYL_RESTRICT out);

GKYL_CU_DH double fpo_vlasov_drag_surfvx_2x3v_ser_p1(const double* dxv, const double *alpha_surf_L, const double *alpha_surf_R, const double *sgn_alpha_surf_L, const double *sgn_alpha_surf_R, const int *const_sgn_alpha_L, const int *const_sgn_alpha_R, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_drag_surfvy_2x3v_ser_p1(const double* dxv, const double *alpha_surf_L, const double *alpha_surf_R, const double *sgn_alpha_surf_L, const double *sgn_alpha_surf_R, const int *const_sgn_alpha_L, const int *const_sgn_alpha_R, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_drag_surfvz_2x3v_ser_p1(const double* dxv, const double *alpha_surf_L, const double *alpha_surf_R, const double *sgn_alpha_surf_L, const double *sgn_alpha_surf_R, const int *const_sgn_alpha_L, const int *const_sgn_alpha_R, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out);

GKYL_CU_DH double fpo_vlasov_drag_boundary_surfvx_2x3v_ser_p1(const double* dxv, const double *alpha_surf_Edge, const double *alpha_surf_Skin, 
  const double *sgn_alpha_surf_Edge, const double *sgn_alpha_surf_Skin, const int *const_sgn_alpha_Edge, const int *const_sgn_alpha_Skin, const int Edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double fpo_vlasov_drag_boundary_surfvy_2x3v_ser_p1(const double* dxv, const double *alpha_surf_Edge, const double *alpha_surf_Skin, 
  const double *sgn_alpha_surf_Edge, const double *sgn_alpha_surf_Skin, const int *const_sgn_alpha_Edge, const int *const_sgn_alpha_Skin, const int Edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double fpo_vlasov_drag_boundary_surfvz_2x3v_ser_p1(const double* dxv, const double *alpha_surf_Edge, const double *alpha_surf_Skin, 
  const double *sgn_alpha_surf_Edge, const double *sgn_alpha_surf_Skin, const int *const_sgn_alpha_Edge, const int *const_sgn_alpha_Skin, const int Edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double fpo_vlasov_diff_vol_2x3v_ser_p1(const double* dxv, const double *diff_coeff, const double* f, double* GKYL_RESTRICT out);

GKYL_CU_DH double fpo_vlasov_diff_surfvxvx_2x3v_ser_p1(const double* dxv, const double* diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvxvy_2x3v_ser_p1(const double* dxv, const double* diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvxvz_2x3v_ser_p1(const double* dxv, const double* diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvyvx_2x3v_ser_p1(const double* dxv, const double* diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvyvy_2x3v_ser_p1(const double* dxv, const double* diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvyvz_2x3v_ser_p1(const double* dxv, const double* diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvzvx_2x3v_ser_p1(const double* dxv, const double* diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvzvy_2x3v_ser_p1(const double* dxv, const double* diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvzvz_2x3v_ser_p1(const double* dxv, const double* diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* GKYL_RESTRICT out);

GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvx_2x3v_ser_p1_lovx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvx_2x3v_ser_p1_upvx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvy_2x3v_ser_p1_lovy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvy_2x3v_ser_p1_upvy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvz_2x3v_ser_p1_lovz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvz_2x3v_ser_p1_upvz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);

GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvy_2x3v_ser_p1_invx_lovy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvy_2x3v_ser_p1_invx_upvy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvy_2x3v_ser_p1_lovx_invy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvy_2x3v_ser_p1_lovx_lovy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvy_2x3v_ser_p1_lovx_upvy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvy_2x3v_ser_p1_upvx_invy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvy_2x3v_ser_p1_upvx_lovy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvy_2x3v_ser_p1_upvx_upvy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvz_2x3v_ser_p1_invx_lovz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvz_2x3v_ser_p1_invx_upvz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvz_2x3v_ser_p1_lovx_invz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvz_2x3v_ser_p1_lovx_lovz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvz_2x3v_ser_p1_lovx_upvz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvz_2x3v_ser_p1_upvx_invz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvz_2x3v_ser_p1_upvx_lovz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvz_2x3v_ser_p1_upvx_upvz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvx_2x3v_ser_p1_invy_lovx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvx_2x3v_ser_p1_invy_upvx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvx_2x3v_ser_p1_lovy_invx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvx_2x3v_ser_p1_lovy_lovx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvx_2x3v_ser_p1_lovy_upvx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvx_2x3v_ser_p1_upvy_invx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvx_2x3v_ser_p1_upvy_lovx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvx_2x3v_ser_p1_upvy_upvx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvz_2x3v_ser_p1_invy_lovz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvz_2x3v_ser_p1_invy_upvz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvz_2x3v_ser_p1_lovy_invz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvz_2x3v_ser_p1_lovy_lovz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvz_2x3v_ser_p1_lovy_upvz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvz_2x3v_ser_p1_upvy_invz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvz_2x3v_ser_p1_upvy_lovz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvz_2x3v_ser_p1_upvy_upvz(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvx_2x3v_ser_p1_invz_lovx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvx_2x3v_ser_p1_invz_upvx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvx_2x3v_ser_p1_lovz_invx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvx_2x3v_ser_p1_lovz_lovx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvx_2x3v_ser_p1_lovz_upvx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvx_2x3v_ser_p1_upvz_invx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvx_2x3v_ser_p1_upvz_lovx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvx_2x3v_ser_p1_upvz_upvx(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvy_2x3v_ser_p1_invz_lovy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvy_2x3v_ser_p1_invz_upvy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvy_2x3v_ser_p1_lovz_invy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvy_2x3v_ser_p1_lovz_lovy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvy_2x3v_ser_p1_lovz_upvy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvy_2x3v_ser_p1_upvz_invy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvy_2x3v_ser_p1_upvz_lovy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvzvy_2x3v_ser_p1_upvz_upvy(const double *dxv, const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out);



EXTERN_C_END