#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH void vlasov_triad_nodal_vel_dense_vol_1x1v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_alpha_quad_x_1x1v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_alpha_quad_x_1x1v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_alpha_quad_vx_1x1v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_alpha_quad_vx_1x1v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_ho_alpha_quad_x_1x1v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_ho_alpha_quad_x_1x1v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_ho_alpha_quad_vx_1x1v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_ho_alpha_quad_vx_1x1v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH void vlasov_triad_nodal_phase_vol_1x1v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 

GKYL_CU_DH int triad_nodal_hamil_phase_alpha_quad_x_1x1v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_alpha_quad_x_1x1v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_alpha_quad_vx_1x1v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_alpha_quad_vx_1x1v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_ho_alpha_quad_x_1x1v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_ho_alpha_quad_x_1x1v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_ho_alpha_quad_vx_1x1v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_ho_alpha_quad_vx_1x1v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH void vlasov_triad_nodal_vel_sparse_vol_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_alpha_quad_x_1x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_alpha_quad_x_1x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_alpha_quad_vx_1x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_alpha_quad_vx_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_alpha_quad_vy_1x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_alpha_quad_vy_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_ho_alpha_quad_x_1x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_ho_alpha_quad_x_1x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_ho_alpha_quad_vx_1x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_ho_alpha_quad_vx_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_ho_alpha_quad_vy_1x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_ho_alpha_quad_vy_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH void vlasov_triad_nodal_vel_dense_vol_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_alpha_quad_x_1x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_alpha_quad_x_1x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_alpha_quad_vx_1x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_alpha_quad_vx_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_alpha_quad_vy_1x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_alpha_quad_vy_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_ho_alpha_quad_x_1x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_ho_alpha_quad_x_1x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_ho_alpha_quad_vx_1x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_ho_alpha_quad_vx_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_ho_alpha_quad_vy_1x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_ho_alpha_quad_vy_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH void vlasov_triad_nodal_phase_vol_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 

GKYL_CU_DH int triad_nodal_hamil_phase_alpha_quad_x_1x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_alpha_quad_x_1x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_alpha_quad_vx_1x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_alpha_quad_vx_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_alpha_quad_vy_1x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_alpha_quad_vy_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_ho_alpha_quad_x_1x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_ho_alpha_quad_x_1x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_ho_alpha_quad_vx_1x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_ho_alpha_quad_vx_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_ho_alpha_quad_vy_1x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_ho_alpha_quad_vy_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH void vlasov_triad_nodal_vel_sparse_vol_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_alpha_quad_x_2x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_alpha_quad_x_2x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_alpha_quad_y_2x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_alpha_quad_y_2x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_alpha_quad_vx_2x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_alpha_quad_vx_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_alpha_quad_vy_2x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_alpha_quad_vy_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_ho_alpha_quad_x_2x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_ho_alpha_quad_x_2x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_ho_alpha_quad_y_2x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_ho_alpha_quad_y_2x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_ho_alpha_quad_vx_2x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_ho_alpha_quad_vx_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_ho_alpha_quad_vy_2x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_ho_alpha_quad_vy_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH void vlasov_triad_nodal_vel_dense_vol_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_alpha_quad_x_2x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_alpha_quad_x_2x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_alpha_quad_y_2x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_alpha_quad_y_2x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_alpha_quad_vx_2x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_alpha_quad_vx_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_alpha_quad_vy_2x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_alpha_quad_vy_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_ho_alpha_quad_x_2x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_ho_alpha_quad_x_2x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_ho_alpha_quad_y_2x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_ho_alpha_quad_y_2x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_ho_alpha_quad_vx_2x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_ho_alpha_quad_vx_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_ho_alpha_quad_vy_2x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_ho_alpha_quad_vy_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH void vlasov_triad_nodal_phase_vol_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 

GKYL_CU_DH int triad_nodal_hamil_phase_alpha_quad_x_2x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_alpha_quad_x_2x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_alpha_quad_y_2x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_alpha_quad_y_2x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_alpha_quad_vx_2x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_alpha_quad_vx_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_alpha_quad_vy_2x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_alpha_quad_vy_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_ho_alpha_quad_x_2x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_ho_alpha_quad_x_2x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_ho_alpha_quad_y_2x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_ho_alpha_quad_y_2x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_ho_alpha_quad_vx_2x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_ho_alpha_quad_vx_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_ho_alpha_quad_vy_2x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_ho_alpha_quad_vy_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH void vlasov_triad_nodal_vel_sparse_vol_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_alpha_quad_x_1x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_alpha_quad_x_1x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_alpha_quad_vx_1x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_alpha_quad_vx_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_alpha_quad_vy_1x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_alpha_quad_vy_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_alpha_quad_vz_1x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_alpha_quad_vz_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_ho_alpha_quad_x_1x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_ho_alpha_quad_x_1x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_ho_alpha_quad_vx_1x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_ho_alpha_quad_vx_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_ho_alpha_quad_vy_1x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_ho_alpha_quad_vy_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_ho_alpha_quad_vz_1x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_ho_alpha_quad_vz_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH void vlasov_triad_nodal_vel_dense_vol_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_alpha_quad_x_1x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_alpha_quad_x_1x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_alpha_quad_vx_1x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_alpha_quad_vx_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_alpha_quad_vy_1x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_alpha_quad_vy_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_alpha_quad_vz_1x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_alpha_quad_vz_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_ho_alpha_quad_x_1x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_ho_alpha_quad_x_1x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_ho_alpha_quad_vx_1x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_ho_alpha_quad_vx_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_ho_alpha_quad_vy_1x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_ho_alpha_quad_vy_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_ho_alpha_quad_vz_1x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_ho_alpha_quad_vz_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH void vlasov_triad_nodal_phase_vol_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 

GKYL_CU_DH int triad_nodal_hamil_phase_alpha_quad_x_1x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_alpha_quad_x_1x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_alpha_quad_vx_1x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_alpha_quad_vx_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_alpha_quad_vy_1x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_alpha_quad_vy_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_alpha_quad_vz_1x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_alpha_quad_vz_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_ho_alpha_quad_x_1x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_ho_alpha_quad_x_1x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_ho_alpha_quad_vx_1x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_ho_alpha_quad_vx_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_ho_alpha_quad_vy_1x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_ho_alpha_quad_vy_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_ho_alpha_quad_vz_1x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_ho_alpha_quad_vz_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH void vlasov_triad_nodal_vel_sparse_vol_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_alpha_quad_x_2x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_alpha_quad_x_2x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_alpha_quad_y_2x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_alpha_quad_y_2x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_alpha_quad_vx_2x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_alpha_quad_vx_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_alpha_quad_vy_2x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_alpha_quad_vy_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_alpha_quad_vz_2x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_alpha_quad_vz_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_ho_alpha_quad_x_2x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_ho_alpha_quad_x_2x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_ho_alpha_quad_y_2x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_ho_alpha_quad_y_2x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_ho_alpha_quad_vx_2x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_ho_alpha_quad_vx_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_ho_alpha_quad_vy_2x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_ho_alpha_quad_vy_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_ho_alpha_quad_vz_2x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_ho_alpha_quad_vz_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH void vlasov_triad_nodal_vel_dense_vol_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_alpha_quad_x_2x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_alpha_quad_x_2x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_alpha_quad_y_2x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_alpha_quad_y_2x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_alpha_quad_vx_2x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_alpha_quad_vx_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_alpha_quad_vy_2x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_alpha_quad_vy_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_alpha_quad_vz_2x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_alpha_quad_vz_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_ho_alpha_quad_x_2x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_ho_alpha_quad_x_2x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_ho_alpha_quad_y_2x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_ho_alpha_quad_y_2x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_ho_alpha_quad_vx_2x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_ho_alpha_quad_vx_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_ho_alpha_quad_vy_2x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_ho_alpha_quad_vy_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_ho_alpha_quad_vz_2x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_ho_alpha_quad_vz_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH void vlasov_triad_nodal_phase_vol_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 

GKYL_CU_DH int triad_nodal_hamil_phase_alpha_quad_x_2x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_alpha_quad_x_2x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_alpha_quad_y_2x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_alpha_quad_y_2x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_alpha_quad_vx_2x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_alpha_quad_vx_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_alpha_quad_vy_2x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_alpha_quad_vy_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_alpha_quad_vz_2x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_alpha_quad_vz_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_ho_alpha_quad_x_2x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_ho_alpha_quad_x_2x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_ho_alpha_quad_y_2x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_ho_alpha_quad_y_2x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_ho_alpha_quad_vx_2x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_ho_alpha_quad_vx_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_ho_alpha_quad_vy_2x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_ho_alpha_quad_vy_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_ho_alpha_quad_vz_2x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_ho_alpha_quad_vz_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH void vlasov_triad_nodal_vel_sparse_vol_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_alpha_quad_x_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_alpha_quad_x_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_alpha_quad_y_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_alpha_quad_y_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_alpha_quad_z_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_alpha_quad_z_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_alpha_quad_vx_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_alpha_quad_vx_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_alpha_quad_vy_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_alpha_quad_vy_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_alpha_quad_vz_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_alpha_quad_vz_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_ho_alpha_quad_x_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_ho_alpha_quad_x_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_ho_alpha_quad_y_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_ho_alpha_quad_y_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_ho_alpha_quad_z_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_ho_alpha_quad_z_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_ho_alpha_quad_vx_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_ho_alpha_quad_vx_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_ho_alpha_quad_vy_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_ho_alpha_quad_vy_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_sparse_ho_alpha_quad_vz_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_sparse_ho_alpha_quad_vz_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH void vlasov_triad_nodal_vel_dense_vol_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_alpha_quad_x_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_alpha_quad_x_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_alpha_quad_y_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_alpha_quad_y_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_alpha_quad_z_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_alpha_quad_z_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_alpha_quad_vx_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_alpha_quad_vx_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_alpha_quad_vy_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_alpha_quad_vy_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_alpha_quad_vz_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_alpha_quad_vz_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_ho_alpha_quad_x_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_ho_alpha_quad_x_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_ho_alpha_quad_y_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_ho_alpha_quad_y_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_ho_alpha_quad_z_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_ho_alpha_quad_z_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_ho_alpha_quad_vx_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_ho_alpha_quad_vx_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_ho_alpha_quad_vy_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_ho_alpha_quad_vy_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_vel_dense_ho_alpha_quad_vz_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_vel_dense_ho_alpha_quad_vz_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH void vlasov_triad_nodal_phase_vol_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 

GKYL_CU_DH int triad_nodal_hamil_phase_alpha_quad_x_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_alpha_quad_x_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_alpha_quad_y_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_alpha_quad_y_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_alpha_quad_z_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_alpha_quad_z_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_alpha_quad_vx_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_alpha_quad_vx_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_alpha_quad_vy_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_alpha_quad_vy_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_alpha_quad_vz_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_alpha_quad_vz_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_ho_alpha_quad_x_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_ho_alpha_quad_x_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_ho_alpha_quad_y_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_ho_alpha_quad_y_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_ho_alpha_quad_z_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_ho_alpha_quad_z_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_ho_alpha_quad_vx_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_ho_alpha_quad_vx_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_ho_alpha_quad_vy_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_ho_alpha_quad_vy_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH int triad_nodal_hamil_phase_ho_alpha_quad_vz_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT O, double* GKYL_RESTRICT I); 
GKYL_CU_DH void triad_nodal_hamil_phase_ho_alpha_quad_vz_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

EXTERN_C_END 
