#pragma once 
#include <math.h> 
#include <gkyl_util.h> 

// Layouts of the inputs of the interval kernels pkpm_remap_flux_interval_* and of the moments. 
#define PKPM_RF_SIDE_R 0 // r: w_s = r wbar + c takes the face-frame label to the label of the side 
#define PKPM_RF_SIDE_C 1 // c, the part of the shift along the field direction of the face 
#define PKPM_RF_SIDE_KINV 2 // density factor of the map from the face frame back to the side 
#define PKPM_RF_SIDE_EG 3 // K_inv/r_inv^2, the factor of the block G mapped back to the side 
#define PKPM_RF_SIDE_CC 4 // |c_perp_inv|^2 
#define PKPM_RF_SIDE_CPERP 5 // first of the three components of c_perp_inv 
#define PKPM_RF_SIDE_NUM 8 // size of side_l/side_r 
#define PKPM_RF_CELL_ETA0 0 // eta_s = ETA0 + ETA1 wbar: cell coordinate of the side at the label wbar 
#define PKPM_RF_CELL_ETA1 1 // 2 r/dw of the velocity cell 
#define PKPM_RF_CELL_OUT 2 // 1 if the side is outside its velocity grid on the interval, else 0 
#define PKPM_RF_CELL_NUM 3 // size of cell_l/cell_r 
#define PKPM_RF_SHARE_MOM_OUT 0 // momentum share outside the grid of side s, at MOM_OUT + s 
#define PKPM_RF_SHARE_MOM_TOT 2 // total momentum share of side s, at MOM_TOT + s 
#define PKPM_RF_SHARE_EN_OUT 4 // energy share outside the grid of side s, at EN_OUT + s 
#define PKPM_RF_SHARE_EN_TOT 6 // total energy share of side s, at EN_TOT + s 
#define PKPM_RF_SHARE_NUM 8 // size of share 
#define PKPM_RF_NUM_MOMS 7 // face-frame moments of the flux: M0, M1, M2 of F, M0 of G, M0 of M0_perp (3) 

EXTERN_C_BEG 

GKYL_CU_DH void pkpm_remap_flux_interval_x_1v_p2(const double *fc, double w_lo, double w_up, const double *side_l, const double *side_r, const double *cell_l, const double *cell_r, const double *Y_l, const double *Y_r, double* GKYL_RESTRICT P_l, double* GKYL_RESTRICT P_r, double* GKYL_RESTRICT mom, double* GKYL_RESTRICT share); 
GKYL_CU_DH void pkpm_remap_flux_label_x_1v_p2(const double *fc, double w_center, double dw, const double *X_l, const double *X_r, double* GKYL_RESTRICT proj, double* GKYL_RESTRICT mom_l, double* GKYL_RESTRICT mom_r, double* GKYL_RESTRICT speed); 
GKYL_CU_DH void pkpm_remap_flux_interval_y_1v_p2(const double *fc, double w_lo, double w_up, const double *side_l, const double *side_r, const double *cell_l, const double *cell_r, const double *Y_l, const double *Y_r, double* GKYL_RESTRICT P_l, double* GKYL_RESTRICT P_r, double* GKYL_RESTRICT mom, double* GKYL_RESTRICT share); 
GKYL_CU_DH void pkpm_remap_flux_label_y_1v_p2(const double *fc, double w_center, double dw, const double *X_l, const double *X_r, double* GKYL_RESTRICT proj, double* GKYL_RESTRICT mom_l, double* GKYL_RESTRICT mom_r, double* GKYL_RESTRICT speed); 
GKYL_CU_DH void pkpm_remap_flux_interval_z_1v_p2(const double *fc, double w_lo, double w_up, const double *side_l, const double *side_r, const double *cell_l, const double *cell_r, const double *Y_l, const double *Y_r, double* GKYL_RESTRICT P_l, double* GKYL_RESTRICT P_r, double* GKYL_RESTRICT mom, double* GKYL_RESTRICT share); 
GKYL_CU_DH void pkpm_remap_flux_label_z_1v_p2(const double *fc, double w_center, double dw, const double *X_l, const double *X_r, double* GKYL_RESTRICT proj, double* GKYL_RESTRICT mom_l, double* GKYL_RESTRICT mom_r, double* GKYL_RESTRICT speed); 

GKYL_CU_DH void pkpm_remap_flux_to_face_1x1v_ser_p1(const double *psi, const double *fwd, const double *G, double* GKYL_RESTRICT X, double* GKYL_RESTRICT Y); 
GKYL_CU_DH void pkpm_remap_flux_from_face_x_1x1v_ser_p1(const double *wpsi, const double *bck, const double *P, double* GKYL_RESTRICT flux); 
GKYL_CU_DH void pkpm_remap_flux_trace_x_1x1v_ser_p1(double jdinv_l, double jdinv_r, const double *f_l, const double *f_r, double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r); 
GKYL_CU_DH void pkpm_remap_flux_store_x_1x1v_ser_p1(const double *wpsi, const double *proj, double* GKYL_RESTRICT flux); 

GKYL_CU_DH void pkpm_remap_flux_to_face_2x1v_ser_p1(const double *psi, const double *fwd, const double *G, double* GKYL_RESTRICT X, double* GKYL_RESTRICT Y); 
GKYL_CU_DH void pkpm_remap_flux_from_face_x_2x1v_ser_p1(const double *wpsi, const double *bck, const double *P, double* GKYL_RESTRICT flux); 
GKYL_CU_DH void pkpm_remap_flux_trace_x_2x1v_ser_p1(double jdinv_l, double jdinv_r, const double *f_l, const double *f_r, double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r); 
GKYL_CU_DH void pkpm_remap_flux_store_x_2x1v_ser_p1(const double *wpsi, const double *proj, double* GKYL_RESTRICT flux); 
GKYL_CU_DH void pkpm_remap_flux_from_face_y_2x1v_ser_p1(const double *wpsi, const double *bck, const double *P, double* GKYL_RESTRICT flux); 
GKYL_CU_DH void pkpm_remap_flux_trace_y_2x1v_ser_p1(double jdinv_l, double jdinv_r, const double *f_l, const double *f_r, double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r); 
GKYL_CU_DH void pkpm_remap_flux_store_y_2x1v_ser_p1(const double *wpsi, const double *proj, double* GKYL_RESTRICT flux); 

GKYL_CU_DH void pkpm_remap_flux_to_face_3x1v_ser_p1(const double *psi, const double *fwd, const double *G, double* GKYL_RESTRICT X, double* GKYL_RESTRICT Y); 
GKYL_CU_DH void pkpm_remap_flux_from_face_x_3x1v_ser_p1(const double *wpsi, const double *bck, const double *P, double* GKYL_RESTRICT flux); 
GKYL_CU_DH void pkpm_remap_flux_trace_x_3x1v_ser_p1(double jdinv_l, double jdinv_r, const double *f_l, const double *f_r, double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r); 
GKYL_CU_DH void pkpm_remap_flux_store_x_3x1v_ser_p1(const double *wpsi, const double *proj, double* GKYL_RESTRICT flux); 
GKYL_CU_DH void pkpm_remap_flux_from_face_y_3x1v_ser_p1(const double *wpsi, const double *bck, const double *P, double* GKYL_RESTRICT flux); 
GKYL_CU_DH void pkpm_remap_flux_trace_y_3x1v_ser_p1(double jdinv_l, double jdinv_r, const double *f_l, const double *f_r, double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r); 
GKYL_CU_DH void pkpm_remap_flux_store_y_3x1v_ser_p1(const double *wpsi, const double *proj, double* GKYL_RESTRICT flux); 
GKYL_CU_DH void pkpm_remap_flux_from_face_z_3x1v_ser_p1(const double *wpsi, const double *bck, const double *P, double* GKYL_RESTRICT flux); 
GKYL_CU_DH void pkpm_remap_flux_trace_z_3x1v_ser_p1(double jdinv_l, double jdinv_r, const double *f_l, const double *f_r, double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r); 
GKYL_CU_DH void pkpm_remap_flux_store_z_3x1v_ser_p1(const double *wpsi, const double *proj, double* GKYL_RESTRICT flux); 

EXTERN_C_END 
