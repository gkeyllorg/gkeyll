#pragma once

#include <math.h>
#include <gkyl_util.h>

EXTERN_C_BEG

GKYL_CU_DH double canonical_pb_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *hamil,  
              const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_1x1v_ser_p1(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_1x1v_ser_p1(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH int canonical_pb_alpha_surfx_1x1v_ser_p1(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfx_1x1v_ser_p1(const double *w, const double *dxv, const double *hamil,
                  double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfx_1x1v_ser_p1(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfx_1x1v_ser_p1(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double canonical_pb_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *hamil,  
              const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_1x1v_ser_p2(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_1x1v_ser_p2(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH int canonical_pb_alpha_surfx_1x1v_ser_p2(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfx_1x1v_ser_p2(const double *w, const double *dxv, const double *hamil,
                  double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfx_1x1v_ser_p2(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfx_1x1v_ser_p2(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double canonical_pb_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *hamil,  
              const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_1x2v_ser_p1(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_1x2v_ser_p1(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH int canonical_pb_alpha_surfx_1x2v_ser_p1(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfx_1x2v_ser_p1(const double *w, const double *dxv, const double *hamil,
                  double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfx_1x2v_ser_p1(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfx_1x2v_ser_p1(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double canonical_pb_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *hamil,  
              const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_1x2v_ser_p2(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_1x2v_ser_p2(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH int canonical_pb_alpha_surfx_1x2v_ser_p2(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfx_1x2v_ser_p2(const double *w, const double *dxv, const double *hamil,
                  double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfx_1x2v_ser_p2(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfx_1x2v_ser_p2(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvy_1x2v_ser_p2(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvy_1x2v_ser_p2(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvy_1x2v_ser_p2(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double canonical_pb_vol_1x3v_ser_p1(const double *w, const double *dxv, const double *hamil,  
              const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_1x3v_ser_p1(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_1x3v_ser_p1(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH int canonical_pb_alpha_surfx_1x3v_ser_p1(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfx_1x3v_ser_p1(const double *w, const double *dxv, const double *hamil,
                  double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfx_1x3v_ser_p1(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfx_1x3v_ser_p1(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvx_1x3v_ser_p1(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvx_1x3v_ser_p1(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvx_1x3v_ser_p1(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvy_1x3v_ser_p1(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvy_1x3v_ser_p1(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvy_1x3v_ser_p1(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvz_1x3v_ser_p1(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvz_1x3v_ser_p1(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvz_1x3v_ser_p1(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double canonical_pb_vol_1x3v_ser_p2(const double *w, const double *dxv, const double *hamil,  
              const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_1x3v_ser_p2(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_1x3v_ser_p2(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH int canonical_pb_alpha_surfx_1x3v_ser_p2(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfx_1x3v_ser_p2(const double *w, const double *dxv, const double *hamil,
                  double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfx_1x3v_ser_p2(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfx_1x3v_ser_p2(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvx_1x3v_ser_p2(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvx_1x3v_ser_p2(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvx_1x3v_ser_p2(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvy_1x3v_ser_p2(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvy_1x3v_ser_p2(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvy_1x3v_ser_p2(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvz_1x3v_ser_p2(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvz_1x3v_ser_p2(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvz_1x3v_ser_p2(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double canonical_pb_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *hamil,  
              const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_2x2v_ser_p1(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_2x2v_ser_p1(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH int canonical_pb_alpha_surfx_2x2v_ser_p1(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfx_2x2v_ser_p1(const double *w, const double *dxv, const double *hamil,
                  double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfx_2x2v_ser_p1(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfx_2x2v_ser_p1(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfy_2x2v_ser_p1(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfy_2x2v_ser_p1(const double *w, const double *dxv, const double *hamil,
                  double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfy_2x2v_ser_p1(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfy_2x2v_ser_p1(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double canonical_pb_vol_2x2v_ser_p2(const double *w, const double *dxv, const double *hamil,  
              const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_2x2v_ser_p2(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_2x2v_ser_p2(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH int canonical_pb_alpha_surfx_2x2v_ser_p2(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfx_2x2v_ser_p2(const double *w, const double *dxv, const double *hamil,
                  double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfx_2x2v_ser_p2(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfx_2x2v_ser_p2(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfy_2x2v_ser_p2(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfy_2x2v_ser_p2(const double *w, const double *dxv, const double *hamil,
                  double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfy_2x2v_ser_p2(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfy_2x2v_ser_p2(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double canonical_pb_vol_2x3v_ser_p1(const double *w, const double *dxv, const double *hamil,  
              const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_2x3v_ser_p1(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_2x3v_ser_p1(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH int canonical_pb_alpha_surfx_2x3v_ser_p1(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfx_2x3v_ser_p1(const double *w, const double *dxv, const double *hamil,
                  double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfx_2x3v_ser_p1(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfx_2x3v_ser_p1(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfy_2x3v_ser_p1(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfy_2x3v_ser_p1(const double *w, const double *dxv, const double *hamil,
                  double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfy_2x3v_ser_p1(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfy_2x3v_ser_p1(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvy_2x3v_ser_p1(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvy_2x3v_ser_p1(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvy_2x3v_ser_p1(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvz_2x3v_ser_p1(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvz_2x3v_ser_p1(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvz_2x3v_ser_p1(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double canonical_pb_vol_2x3v_ser_p2(const double *w, const double *dxv, const double *hamil,  
              const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_2x3v_ser_p2(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_2x3v_ser_p2(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH int canonical_pb_alpha_surfx_2x3v_ser_p2(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfx_2x3v_ser_p2(const double *w, const double *dxv, const double *hamil,
                  double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfx_2x3v_ser_p2(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfx_2x3v_ser_p2(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfy_2x3v_ser_p2(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfy_2x3v_ser_p2(const double *w, const double *dxv, const double *hamil,
                  double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfy_2x3v_ser_p2(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfy_2x3v_ser_p2(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvx_2x3v_ser_p2(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvx_2x3v_ser_p2(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvx_2x3v_ser_p2(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvy_2x3v_ser_p2(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvy_2x3v_ser_p2(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvy_2x3v_ser_p2(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvz_2x3v_ser_p2(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvz_2x3v_ser_p2(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvz_2x3v_ser_p2(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double canonical_pb_vol_1x1v_tensor_p1(const double *w, const double *dxv, const double *hamil,  
              const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_1x1v_tensor_p1(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_1x1v_tensor_p1(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH int canonical_pb_alpha_surfx_1x1v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfx_1x1v_tensor_p1(const double *w, const double *dxv, const double *hamil,
                  double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfx_1x1v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfx_1x1v_tensor_p1(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvx_1x1v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvx_1x1v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvx_1x1v_tensor_p1(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double canonical_pb_vol_1x1v_tensor_p2(const double *w, const double *dxv, const double *hamil,  
              const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_1x1v_tensor_p2(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_1x1v_tensor_p2(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH int canonical_pb_alpha_surfx_1x1v_tensor_p2(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfx_1x1v_tensor_p2(const double *w, const double *dxv, const double *hamil,
                  double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfx_1x1v_tensor_p2(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfx_1x1v_tensor_p2(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvx_1x1v_tensor_p2(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvx_1x1v_tensor_p2(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvx_1x1v_tensor_p2(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double canonical_pb_vol_1x2v_tensor_p1(const double *w, const double *dxv, const double *hamil,  
              const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_1x2v_tensor_p1(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_1x2v_tensor_p1(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH int canonical_pb_alpha_surfx_1x2v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfx_1x2v_tensor_p1(const double *w, const double *dxv, const double *hamil,
                  double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfx_1x2v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfx_1x2v_tensor_p1(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvx_1x2v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvx_1x2v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvx_1x2v_tensor_p1(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvy_1x2v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvy_1x2v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvy_1x2v_tensor_p1(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double canonical_pb_vol_1x2v_tensor_p2(const double *w, const double *dxv, const double *hamil,  
              const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_1x2v_tensor_p2(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_1x2v_tensor_p2(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH int canonical_pb_alpha_surfx_1x2v_tensor_p2(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfx_1x2v_tensor_p2(const double *w, const double *dxv, const double *hamil,
                  double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfx_1x2v_tensor_p2(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfx_1x2v_tensor_p2(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvx_1x2v_tensor_p2(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvx_1x2v_tensor_p2(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvx_1x2v_tensor_p2(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvy_1x2v_tensor_p2(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvy_1x2v_tensor_p2(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvy_1x2v_tensor_p2(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double canonical_pb_vol_1x3v_tensor_p1(const double *w, const double *dxv, const double *hamil,  
              const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_1x3v_tensor_p1(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_1x3v_tensor_p1(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH int canonical_pb_alpha_surfx_1x3v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfx_1x3v_tensor_p1(const double *w, const double *dxv, const double *hamil,
                  double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfx_1x3v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfx_1x3v_tensor_p1(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvx_1x3v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvx_1x3v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvx_1x3v_tensor_p1(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvy_1x3v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvy_1x3v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvy_1x3v_tensor_p1(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvz_1x3v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvz_1x3v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvz_1x3v_tensor_p1(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double canonical_pb_vol_1x3v_tensor_p2(const double *w, const double *dxv, const double *hamil,  
              const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_1x3v_tensor_p2(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_1x3v_tensor_p2(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH int canonical_pb_alpha_surfx_1x3v_tensor_p2(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfx_1x3v_tensor_p2(const double *w, const double *dxv, const double *hamil,
                  double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfx_1x3v_tensor_p2(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfx_1x3v_tensor_p2(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvx_1x3v_tensor_p2(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvx_1x3v_tensor_p2(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvx_1x3v_tensor_p2(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvy_1x3v_tensor_p2(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvy_1x3v_tensor_p2(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvy_1x3v_tensor_p2(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvz_1x3v_tensor_p2(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvz_1x3v_tensor_p2(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvz_1x3v_tensor_p2(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double canonical_pb_vol_2x2v_tensor_p1(const double *w, const double *dxv, const double *hamil,  
              const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_2x2v_tensor_p1(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_2x2v_tensor_p1(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH int canonical_pb_alpha_surfx_2x2v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfx_2x2v_tensor_p1(const double *w, const double *dxv, const double *hamil,
                  double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfx_2x2v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfx_2x2v_tensor_p1(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfy_2x2v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfy_2x2v_tensor_p1(const double *w, const double *dxv, const double *hamil,
                  double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfy_2x2v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfy_2x2v_tensor_p1(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvx_2x2v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvx_2x2v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvx_2x2v_tensor_p1(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvy_2x2v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvy_2x2v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvy_2x2v_tensor_p1(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double canonical_pb_vol_2x2v_tensor_p2(const double *w, const double *dxv, const double *hamil,  
              const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_2x2v_tensor_p2(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_2x2v_tensor_p2(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH int canonical_pb_alpha_surfx_2x2v_tensor_p2(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfx_2x2v_tensor_p2(const double *w, const double *dxv, const double *hamil,
                  double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfx_2x2v_tensor_p2(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfx_2x2v_tensor_p2(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfy_2x2v_tensor_p2(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfy_2x2v_tensor_p2(const double *w, const double *dxv, const double *hamil,
                  double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfy_2x2v_tensor_p2(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfy_2x2v_tensor_p2(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvx_2x2v_tensor_p2(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvx_2x2v_tensor_p2(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvx_2x2v_tensor_p2(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvy_2x2v_tensor_p2(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvy_2x2v_tensor_p2(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvy_2x2v_tensor_p2(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double canonical_pb_vol_2x3v_tensor_p1(const double *w, const double *dxv, const double *hamil,  
              const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_2x3v_tensor_p1(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_2x3v_tensor_p1(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH int canonical_pb_alpha_surfx_2x3v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfx_2x3v_tensor_p1(const double *w, const double *dxv, const double *hamil,
                  double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfx_2x3v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfx_2x3v_tensor_p1(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfy_2x3v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfy_2x3v_tensor_p1(const double *w, const double *dxv, const double *hamil,
                  double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfy_2x3v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfy_2x3v_tensor_p1(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvx_2x3v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvx_2x3v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvx_2x3v_tensor_p1(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvy_2x3v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvy_2x3v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvy_2x3v_tensor_p1(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvz_2x3v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
                double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvz_2x3v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvz_2x3v_tensor_p1(const double *w, const double *dxv, const double *hamil,
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 


GKYL_CU_DH double canonical_pb_vol_2x_ser_p1(const double *w, const double *dxv, const double *phi, const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_two_fluid_vol_2x_ser_p1(const double *w, const double *dxv, const double *phi, const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_fluid_hasegawa_mima_source_2x_ser_p1(const double *dxv, double alpha, const double *phi, const double *n0, const double *adiabatic_coupling_phi_n, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_fluid_hasegawa_wakatani_source_2x_ser_p1(const double *dxv, double alpha, const double *phi, const double *n0, const double *adiabatic_coupling_phi_n, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_fluid_subtract_zonal_2x_ser_p1(const double *phi_zonal, const double *n_zonal, double* GKYL_RESTRICT adiabatic_coupling_phi_n); 
GKYL_CU_DH int canonical_pb_alpha_surfx_2x_ser_p1(const double *w, const double *dxv, const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfx_2x_ser_p1(const double *w, const double *dxv, const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfx_2x_ser_p1(const double *w, const double *dxv, const double *phi, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_two_fluid_surfx_2x_ser_p1(const double *w, const double *dxv, const double *phi, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfy_2x_ser_p1(const double *w, const double *dxv, const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfy_2x_ser_p1(const double *w, const double *dxv, const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfy_2x_ser_p1(const double *w, const double *dxv, const double *phi, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_two_fluid_surfy_2x_ser_p1(const double *w, const double *dxv, const double *phi, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double canonical_pb_vol_2x_ser_p2(const double *w, const double *dxv, const double *phi, const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_two_fluid_vol_2x_ser_p2(const double *w, const double *dxv, const double *phi, const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_fluid_hasegawa_mima_source_2x_ser_p2(const double *dxv, double alpha, const double *phi, const double *n0, const double *adiabatic_coupling_phi_n, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_fluid_hasegawa_wakatani_source_2x_ser_p2(const double *dxv, double alpha, const double *phi, const double *n0, const double *adiabatic_coupling_phi_n, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_fluid_subtract_zonal_2x_ser_p2(const double *phi_zonal, const double *n_zonal, double* GKYL_RESTRICT adiabatic_coupling_phi_n); 
GKYL_CU_DH int canonical_pb_alpha_surfx_2x_ser_p2(const double *w, const double *dxv, const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfx_2x_ser_p2(const double *w, const double *dxv, const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfx_2x_ser_p2(const double *w, const double *dxv, const double *phi, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_two_fluid_surfx_2x_ser_p2(const double *w, const double *dxv, const double *phi, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfy_2x_ser_p2(const double *w, const double *dxv, const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfy_2x_ser_p2(const double *w, const double *dxv, const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfy_2x_ser_p2(const double *w, const double *dxv, const double *phi, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_two_fluid_surfy_2x_ser_p2(const double *w, const double *dxv, const double *phi, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double canonical_pb_vol_2x_tensor_p2(const double *w, const double *dxv, const double *phi, const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_two_fluid_vol_2x_tensor_p2(const double *w, const double *dxv, const double *phi, const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_fluid_hasegawa_mima_source_2x_tensor_p2(const double *dxv, double alpha, const double *phi, const double *n0, const double *adiabatic_coupling_phi_n, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_fluid_hasegawa_wakatani_source_2x_tensor_p2(const double *dxv, double alpha, const double *phi, const double *n0, const double *adiabatic_coupling_phi_n, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_fluid_subtract_zonal_2x_tensor_p2(const double *phi_zonal, const double *n_zonal, double* GKYL_RESTRICT adiabatic_coupling_phi_n); 
GKYL_CU_DH int canonical_pb_alpha_surfx_2x_tensor_p2(const double *w, const double *dxv, const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfx_2x_tensor_p2(const double *w, const double *dxv, const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfx_2x_tensor_p2(const double *w, const double *dxv, const double *phi, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_two_fluid_surfx_2x_tensor_p2(const double *w, const double *dxv, const double *phi, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfy_2x_tensor_p2(const double *w, const double *dxv, const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfy_2x_tensor_p2(const double *w, const double *dxv, const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfy_2x_tensor_p2(const double *w, const double *dxv, const double *phi, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_two_fluid_surfy_2x_tensor_p2(const double *w, const double *dxv, const double *phi, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void canonical_pb_fluid_default_source(const double *dxv, double alpha, const double *phi, const double *n0, const double *f, double* GKYL_RESTRICT out); 

EXTERN_C_END
