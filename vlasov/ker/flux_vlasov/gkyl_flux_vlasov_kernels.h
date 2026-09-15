#pragma once
#include <math.h>
#include <gkyl_util.h>
// Declares the per-node surface kernels with C linkage so their definitions
// in the kernel .c files (compiled as C++ under nvcc) match the unmangled
// references from the dispatch code.
#include <gkyl_vlasov_surf_node_kernels.h>
EXTERN_C_BEG

GKYL_CU_DH double lax_flux_nodal_x_1x1v_ser_p1(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_x_1x1v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_x_1x1v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vx_1x1v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vx_1x1v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vx_1x1v_ser_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vx_1x1v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vx_1x1v_ser_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vx_1x1v_ser_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vx_1x1v_ser_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vx_1x1v_ser_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vx_1x1v_ser_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH double lax_flux_nodal_x_1x1v_ser_p2(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_x_1x1v_ser_p2(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_x_1x1v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_dense_ho_alpha_quad_x_1x1v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_x_1x1v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_x_1x1v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vx_1x1v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_ho_alpha_quad_vx_1x1v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vx_1x1v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_ho_alpha_quad_vx_1x1v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vx_1x1v_ser_p2(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_vx_1x1v_ser_p2(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vx_1x1v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vx_1x1v_ser_p2(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vx_1x1v_ser_p2(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vx_1x1v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vx_1x1v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vx_1x1v_ser_p2(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_vx_1x1v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_ho_alpha_quad_vx_1x1v_ser_p2(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_ho_alpha_quad_vx_1x1v_ser_p2(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_dense_alpha_quad_vx_1x1v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vx_1x1v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_ho_alpha_quad_vx_1x1v_ser_p2(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH double lax_flux_nodal_x_1x1v_ser_p3(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_x_1x1v_ser_p3(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_x_1x1v_ser_p3(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_dense_ho_alpha_quad_x_1x1v_ser_p3(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_x_1x1v_ser_p3(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_x_1x1v_ser_p3(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vx_1x1v_ser_p3(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_vx_1x1v_ser_p3(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vx_1x1v_ser_p3(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vx_1x1v_ser_p3(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vx_1x1v_ser_p3(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vx_1x1v_ser_p3(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vx_1x1v_ser_p3(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vx_1x1v_ser_p3(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_vx_1x1v_ser_p3(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_ho_alpha_quad_vx_1x1v_ser_p3(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_ho_alpha_quad_vx_1x1v_ser_p3(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_dense_alpha_quad_vx_1x1v_ser_p3(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vx_1x1v_ser_p3(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_ho_alpha_quad_vx_1x1v_ser_p3(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH double lax_flux_nodal_x_1x2v_ser_p1(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_x_1x2v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_alpha_quad_x_1x2v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_x_1x2v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vx_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vx_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vx_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vx_1x2v_ser_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vx_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vx_1x2v_ser_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vx_1x2v_ser_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vx_1x2v_ser_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vx_1x2v_ser_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vx_1x2v_ser_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vx_1x2v_ser_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vy_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vy_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vy_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vy_1x2v_ser_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vy_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vy_1x2v_ser_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vy_1x2v_ser_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vy_1x2v_ser_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vy_1x2v_ser_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vy_1x2v_ser_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vy_1x2v_ser_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH double lax_flux_nodal_x_1x2v_ser_p2(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_x_1x2v_ser_p2(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_x_1x2v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_alpha_quad_x_1x2v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_dense_ho_alpha_quad_x_1x2v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_ho_alpha_quad_x_1x2v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_x_1x2v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_x_1x2v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vx_1x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vx_1x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_ho_alpha_quad_vx_1x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_ho_alpha_quad_vx_1x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vx_1x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_ho_alpha_quad_vx_1x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vx_1x2v_ser_p2(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_vx_1x2v_ser_p2(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vx_1x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vx_1x2v_ser_p2(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vx_1x2v_ser_p2(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vx_1x2v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vx_1x2v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vx_1x2v_ser_p2(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vx_1x2v_ser_p2(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_vx_1x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_ho_alpha_quad_vx_1x2v_ser_p2(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_ho_alpha_quad_vx_1x2v_ser_p2(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_dense_alpha_quad_vx_1x2v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vx_1x2v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_sparse_alpha_quad_vx_1x2v_ser_p2(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_ho_alpha_quad_vx_1x2v_ser_p2(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vy_1x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vy_1x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_ho_alpha_quad_vy_1x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_ho_alpha_quad_vy_1x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vy_1x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_ho_alpha_quad_vy_1x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vy_1x2v_ser_p2(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_vy_1x2v_ser_p2(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vy_1x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vy_1x2v_ser_p2(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vy_1x2v_ser_p2(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vy_1x2v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vy_1x2v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vy_1x2v_ser_p2(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vy_1x2v_ser_p2(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_vy_1x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_ho_alpha_quad_vy_1x2v_ser_p2(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_ho_alpha_quad_vy_1x2v_ser_p2(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_dense_alpha_quad_vy_1x2v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vy_1x2v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_sparse_alpha_quad_vy_1x2v_ser_p2(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_ho_alpha_quad_vy_1x2v_ser_p2(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH double lax_flux_nodal_x_1x3v_ser_p1(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_x_1x3v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_alpha_quad_x_1x3v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_x_1x3v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vx_1x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vx_1x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vx_1x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vx_1x3v_ser_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vx_1x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vx_1x3v_ser_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vx_1x3v_ser_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vx_1x3v_ser_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vx_1x3v_ser_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vx_1x3v_ser_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vx_1x3v_ser_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vy_1x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vy_1x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vy_1x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vy_1x3v_ser_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vy_1x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vy_1x3v_ser_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vy_1x3v_ser_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vy_1x3v_ser_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vy_1x3v_ser_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vy_1x3v_ser_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vy_1x3v_ser_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vz_1x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vz_1x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vz_1x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vz_1x3v_ser_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vz_1x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vz_1x3v_ser_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vz_1x3v_ser_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vz_1x3v_ser_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vz_1x3v_ser_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vz_1x3v_ser_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vz_1x3v_ser_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH double lax_flux_nodal_x_1x3v_ser_p2(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_x_1x3v_ser_p2(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_x_1x3v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_alpha_quad_x_1x3v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_dense_ho_alpha_quad_x_1x3v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_ho_alpha_quad_x_1x3v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_x_1x3v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_x_1x3v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vx_1x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vx_1x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_ho_alpha_quad_vx_1x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_ho_alpha_quad_vx_1x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vx_1x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_ho_alpha_quad_vx_1x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vx_1x3v_ser_p2(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_vx_1x3v_ser_p2(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vx_1x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vx_1x3v_ser_p2(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vx_1x3v_ser_p2(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vx_1x3v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vx_1x3v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vx_1x3v_ser_p2(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vx_1x3v_ser_p2(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_vx_1x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_ho_alpha_quad_vx_1x3v_ser_p2(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_ho_alpha_quad_vx_1x3v_ser_p2(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_dense_alpha_quad_vx_1x3v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vx_1x3v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_sparse_alpha_quad_vx_1x3v_ser_p2(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_ho_alpha_quad_vx_1x3v_ser_p2(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vy_1x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vy_1x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_ho_alpha_quad_vy_1x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_ho_alpha_quad_vy_1x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vy_1x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_ho_alpha_quad_vy_1x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vy_1x3v_ser_p2(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_vy_1x3v_ser_p2(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vy_1x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vy_1x3v_ser_p2(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vy_1x3v_ser_p2(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vy_1x3v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vy_1x3v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vy_1x3v_ser_p2(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vy_1x3v_ser_p2(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_vy_1x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_ho_alpha_quad_vy_1x3v_ser_p2(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_ho_alpha_quad_vy_1x3v_ser_p2(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_dense_alpha_quad_vy_1x3v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vy_1x3v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_sparse_alpha_quad_vy_1x3v_ser_p2(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_ho_alpha_quad_vy_1x3v_ser_p2(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vz_1x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vz_1x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_ho_alpha_quad_vz_1x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_ho_alpha_quad_vz_1x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vz_1x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_ho_alpha_quad_vz_1x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vz_1x3v_ser_p2(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_vz_1x3v_ser_p2(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vz_1x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vz_1x3v_ser_p2(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vz_1x3v_ser_p2(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vz_1x3v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vz_1x3v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vz_1x3v_ser_p2(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vz_1x3v_ser_p2(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_vz_1x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_ho_alpha_quad_vz_1x3v_ser_p2(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_ho_alpha_quad_vz_1x3v_ser_p2(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_dense_alpha_quad_vz_1x3v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vz_1x3v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_sparse_alpha_quad_vz_1x3v_ser_p2(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_ho_alpha_quad_vz_1x3v_ser_p2(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH double lax_flux_nodal_x_2x1v_ser_p1(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_x_2x1v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_x_2x1v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_y_2x1v_ser_p1(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_y_2x1v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_y_2x1v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vx_2x1v_ser_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vx_2x1v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vx_2x1v_ser_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vx_2x1v_ser_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vx_2x1v_ser_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vx_2x1v_ser_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vx_2x1v_ser_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH double lax_flux_nodal_x_2x1v_ser_p2(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_x_2x1v_ser_p2(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_x_2x1v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_dense_ho_alpha_quad_x_2x1v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_x_2x1v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_x_2x1v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_y_2x1v_ser_p2(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_y_2x1v_ser_p2(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_y_2x1v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_dense_ho_alpha_quad_y_2x1v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_y_2x1v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_y_2x1v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vx_2x1v_ser_p2(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_vx_2x1v_ser_p2(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vx_2x1v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vx_2x1v_ser_p2(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vx_2x1v_ser_p2(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vx_2x1v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vx_2x1v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vx_2x1v_ser_p2(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_vx_2x1v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_ho_alpha_quad_vx_2x1v_ser_p2(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_ho_alpha_quad_vx_2x1v_ser_p2(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_dense_alpha_quad_vx_2x1v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vx_2x1v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_ho_alpha_quad_vx_2x1v_ser_p2(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH double lax_flux_nodal_x_2x1v_ser_p3(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_x_2x1v_ser_p3(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_x_2x1v_ser_p3(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_dense_ho_alpha_quad_x_2x1v_ser_p3(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_x_2x1v_ser_p3(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_x_2x1v_ser_p3(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_y_2x1v_ser_p3(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_y_2x1v_ser_p3(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_y_2x1v_ser_p3(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_dense_ho_alpha_quad_y_2x1v_ser_p3(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_y_2x1v_ser_p3(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_y_2x1v_ser_p3(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vx_2x1v_ser_p3(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_vx_2x1v_ser_p3(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vx_2x1v_ser_p3(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vx_2x1v_ser_p3(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vx_2x1v_ser_p3(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vx_2x1v_ser_p3(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vx_2x1v_ser_p3(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vx_2x1v_ser_p3(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_vx_2x1v_ser_p3(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_ho_alpha_quad_vx_2x1v_ser_p3(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_ho_alpha_quad_vx_2x1v_ser_p3(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_dense_alpha_quad_vx_2x1v_ser_p3(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vx_2x1v_ser_p3(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_ho_alpha_quad_vx_2x1v_ser_p3(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH double lax_flux_nodal_x_2x2v_ser_p1(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_x_2x2v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_alpha_quad_x_2x2v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_x_2x2v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_y_2x2v_ser_p1(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_y_2x2v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_alpha_quad_y_2x2v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_y_2x2v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vx_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vx_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vx_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vx_2x2v_ser_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vx_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vx_2x2v_ser_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vx_2x2v_ser_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vx_2x2v_ser_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vx_2x2v_ser_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vx_2x2v_ser_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vx_2x2v_ser_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vy_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vy_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vy_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vy_2x2v_ser_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vy_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vy_2x2v_ser_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vy_2x2v_ser_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vy_2x2v_ser_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vy_2x2v_ser_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vy_2x2v_ser_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vy_2x2v_ser_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH double lax_flux_nodal_x_2x2v_ser_p2(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_x_2x2v_ser_p2(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_x_2x2v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_alpha_quad_x_2x2v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_dense_ho_alpha_quad_x_2x2v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_ho_alpha_quad_x_2x2v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_x_2x2v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_x_2x2v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_y_2x2v_ser_p2(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_y_2x2v_ser_p2(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_y_2x2v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_alpha_quad_y_2x2v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_dense_ho_alpha_quad_y_2x2v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_ho_alpha_quad_y_2x2v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_y_2x2v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_y_2x2v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vx_2x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vx_2x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_ho_alpha_quad_vx_2x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_ho_alpha_quad_vx_2x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vx_2x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_ho_alpha_quad_vx_2x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vx_2x2v_ser_p2(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_vx_2x2v_ser_p2(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vx_2x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vx_2x2v_ser_p2(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vx_2x2v_ser_p2(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vx_2x2v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vx_2x2v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vx_2x2v_ser_p2(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vx_2x2v_ser_p2(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_vx_2x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_ho_alpha_quad_vx_2x2v_ser_p2(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_ho_alpha_quad_vx_2x2v_ser_p2(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_dense_alpha_quad_vx_2x2v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vx_2x2v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_sparse_alpha_quad_vx_2x2v_ser_p2(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_ho_alpha_quad_vx_2x2v_ser_p2(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vy_2x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vy_2x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_ho_alpha_quad_vy_2x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_ho_alpha_quad_vy_2x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vy_2x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_ho_alpha_quad_vy_2x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vy_2x2v_ser_p2(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_vy_2x2v_ser_p2(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vy_2x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vy_2x2v_ser_p2(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vy_2x2v_ser_p2(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vy_2x2v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vy_2x2v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vy_2x2v_ser_p2(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vy_2x2v_ser_p2(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_vy_2x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_ho_alpha_quad_vy_2x2v_ser_p2(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_ho_alpha_quad_vy_2x2v_ser_p2(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_dense_alpha_quad_vy_2x2v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vy_2x2v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_sparse_alpha_quad_vy_2x2v_ser_p2(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_ho_alpha_quad_vy_2x2v_ser_p2(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH double lax_flux_nodal_x_2x3v_ser_p1(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_x_2x3v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_alpha_quad_x_2x3v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_x_2x3v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_y_2x3v_ser_p1(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_y_2x3v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_alpha_quad_y_2x3v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_y_2x3v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vx_2x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vx_2x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vx_2x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vx_2x3v_ser_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vx_2x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vx_2x3v_ser_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vx_2x3v_ser_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vx_2x3v_ser_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vx_2x3v_ser_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vx_2x3v_ser_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vx_2x3v_ser_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vy_2x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vy_2x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vy_2x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vy_2x3v_ser_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vy_2x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vy_2x3v_ser_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vy_2x3v_ser_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vy_2x3v_ser_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vy_2x3v_ser_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vy_2x3v_ser_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vy_2x3v_ser_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vz_2x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vz_2x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vz_2x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vz_2x3v_ser_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vz_2x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vz_2x3v_ser_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vz_2x3v_ser_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vz_2x3v_ser_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vz_2x3v_ser_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vz_2x3v_ser_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vz_2x3v_ser_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH double lax_flux_nodal_x_2x3v_ser_p2(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_x_2x3v_ser_p2(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_x_2x3v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_alpha_quad_x_2x3v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_dense_ho_alpha_quad_x_2x3v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_ho_alpha_quad_x_2x3v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_x_2x3v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_x_2x3v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_y_2x3v_ser_p2(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_y_2x3v_ser_p2(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_y_2x3v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_alpha_quad_y_2x3v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_dense_ho_alpha_quad_y_2x3v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_ho_alpha_quad_y_2x3v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_y_2x3v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_y_2x3v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vx_2x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vx_2x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_ho_alpha_quad_vx_2x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_ho_alpha_quad_vx_2x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vx_2x3v_ser_p2(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_vx_2x3v_ser_p2(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void phi_alpha_quad_vx_2x3v_ser_p2(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vx_2x3v_ser_p2(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vx_2x3v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vx_2x3v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vx_2x3v_ser_p2(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vx_2x3v_ser_p2(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_ho_alpha_quad_vx_2x3v_ser_p2(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_ho_alpha_quad_vx_2x3v_ser_p2(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_dense_alpha_quad_vx_2x3v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vx_2x3v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_sparse_alpha_quad_vx_2x3v_ser_p2(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_ho_alpha_quad_vx_2x3v_ser_p2(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vy_2x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vy_2x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_ho_alpha_quad_vy_2x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_ho_alpha_quad_vy_2x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vy_2x3v_ser_p2(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_vy_2x3v_ser_p2(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void phi_alpha_quad_vy_2x3v_ser_p2(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vy_2x3v_ser_p2(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vy_2x3v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vy_2x3v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vy_2x3v_ser_p2(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vy_2x3v_ser_p2(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_ho_alpha_quad_vy_2x3v_ser_p2(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_ho_alpha_quad_vy_2x3v_ser_p2(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_dense_alpha_quad_vy_2x3v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vy_2x3v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_sparse_alpha_quad_vy_2x3v_ser_p2(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_ho_alpha_quad_vy_2x3v_ser_p2(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vz_2x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vz_2x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_ho_alpha_quad_vz_2x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_ho_alpha_quad_vz_2x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vz_2x3v_ser_p2(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_vz_2x3v_ser_p2(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void phi_alpha_quad_vz_2x3v_ser_p2(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vz_2x3v_ser_p2(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vz_2x3v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vz_2x3v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vz_2x3v_ser_p2(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vz_2x3v_ser_p2(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_ho_alpha_quad_vz_2x3v_ser_p2(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_ho_alpha_quad_vz_2x3v_ser_p2(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_dense_alpha_quad_vz_2x3v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vz_2x3v_ser_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_sparse_alpha_quad_vz_2x3v_ser_p2(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_ho_alpha_quad_vz_2x3v_ser_p2(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH double lax_flux_nodal_x_3x3v_ser_p1(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_x_3x3v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_alpha_quad_x_3x3v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_x_3x3v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_y_3x3v_ser_p1(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_y_3x3v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_alpha_quad_y_3x3v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_y_3x3v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_z_3x3v_ser_p1(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_z_3x3v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_alpha_quad_z_3x3v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_z_3x3v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vx_3x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vx_3x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vx_3x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vx_3x3v_ser_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vx_3x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vx_3x3v_ser_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vx_3x3v_ser_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vx_3x3v_ser_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vx_3x3v_ser_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vx_3x3v_ser_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vx_3x3v_ser_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vy_3x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vy_3x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vy_3x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vy_3x3v_ser_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vy_3x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vy_3x3v_ser_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vy_3x3v_ser_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vy_3x3v_ser_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vy_3x3v_ser_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vy_3x3v_ser_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vy_3x3v_ser_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vz_3x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vz_3x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vz_3x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
          const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vz_3x3v_ser_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vz_3x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vz_3x3v_ser_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vz_3x3v_ser_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vz_3x3v_ser_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vz_3x3v_ser_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vz_3x3v_ser_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vz_3x3v_ser_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH double lax_flux_nodal_x_1x1v_tensor_p1(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_x_1x1v_tensor_p1(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_x_1x1v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_dense_ho_alpha_quad_x_1x1v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_x_1x1v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_x_1x1v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vx_1x1v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_ho_alpha_quad_vx_1x1v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vx_1x1v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_ho_alpha_quad_vx_1x1v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vx_1x1v_tensor_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_vx_1x1v_tensor_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vx_1x1v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vx_1x1v_tensor_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vx_1x1v_tensor_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vx_1x1v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vx_1x1v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vx_1x1v_tensor_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_vx_1x1v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_ho_alpha_quad_vx_1x1v_tensor_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_ho_alpha_quad_vx_1x1v_tensor_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_dense_alpha_quad_vx_1x1v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vx_1x1v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_ho_alpha_quad_vx_1x1v_tensor_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH double lax_flux_nodal_x_1x1v_tensor_p2(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_x_1x1v_tensor_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_x_1x1v_tensor_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vx_1x1v_tensor_p2(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vx_1x1v_tensor_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vx_1x1v_tensor_p2(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vx_1x1v_tensor_p2(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vx_1x1v_tensor_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vx_1x1v_tensor_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vx_1x1v_tensor_p2(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH double lax_flux_nodal_x_1x1v_tensor_p3(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_x_1x1v_tensor_p3(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_x_1x1v_tensor_p3(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vx_1x1v_tensor_p3(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vx_1x1v_tensor_p3(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vx_1x1v_tensor_p3(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vx_1x1v_tensor_p3(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vx_1x1v_tensor_p3(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vx_1x1v_tensor_p3(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vx_1x1v_tensor_p3(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH double lax_flux_nodal_x_1x2v_tensor_p1(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_x_1x2v_tensor_p1(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_x_1x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_alpha_quad_x_1x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_dense_ho_alpha_quad_x_1x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_ho_alpha_quad_x_1x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_x_1x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_x_1x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vx_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vx_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_ho_alpha_quad_vx_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_ho_alpha_quad_vx_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vx_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_ho_alpha_quad_vx_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vx_1x2v_tensor_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_vx_1x2v_tensor_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vx_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vx_1x2v_tensor_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vx_1x2v_tensor_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vx_1x2v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vx_1x2v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vx_1x2v_tensor_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vx_1x2v_tensor_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_vx_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_ho_alpha_quad_vx_1x2v_tensor_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_ho_alpha_quad_vx_1x2v_tensor_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_dense_alpha_quad_vx_1x2v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vx_1x2v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_sparse_alpha_quad_vx_1x2v_tensor_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_ho_alpha_quad_vx_1x2v_tensor_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vy_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vy_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_ho_alpha_quad_vy_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_ho_alpha_quad_vy_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vy_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_ho_alpha_quad_vy_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vy_1x2v_tensor_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_vy_1x2v_tensor_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vy_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vy_1x2v_tensor_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vy_1x2v_tensor_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vy_1x2v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vy_1x2v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vy_1x2v_tensor_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vy_1x2v_tensor_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_vy_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_ho_alpha_quad_vy_1x2v_tensor_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_ho_alpha_quad_vy_1x2v_tensor_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_dense_alpha_quad_vy_1x2v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vy_1x2v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_sparse_alpha_quad_vy_1x2v_tensor_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_ho_alpha_quad_vy_1x2v_tensor_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH double lax_flux_nodal_x_1x2v_tensor_p2(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_x_1x2v_tensor_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_alpha_quad_x_1x2v_tensor_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_x_1x2v_tensor_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vx_1x2v_tensor_p2(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vx_1x2v_tensor_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vx_1x2v_tensor_p2(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vx_1x2v_tensor_p2(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vx_1x2v_tensor_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vx_1x2v_tensor_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vx_1x2v_tensor_p2(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vx_1x2v_tensor_p2(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vy_1x2v_tensor_p2(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vy_1x2v_tensor_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vy_1x2v_tensor_p2(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vy_1x2v_tensor_p2(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vy_1x2v_tensor_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vy_1x2v_tensor_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vy_1x2v_tensor_p2(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vy_1x2v_tensor_p2(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH double lax_flux_nodal_x_1x3v_tensor_p1(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_x_1x3v_tensor_p1(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_x_1x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_alpha_quad_x_1x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_dense_ho_alpha_quad_x_1x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_ho_alpha_quad_x_1x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_x_1x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_x_1x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vx_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vx_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_ho_alpha_quad_vx_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_ho_alpha_quad_vx_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vx_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_ho_alpha_quad_vx_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vx_1x3v_tensor_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_vx_1x3v_tensor_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vx_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vx_1x3v_tensor_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vx_1x3v_tensor_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vx_1x3v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vx_1x3v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vx_1x3v_tensor_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vx_1x3v_tensor_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_vx_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_ho_alpha_quad_vx_1x3v_tensor_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_ho_alpha_quad_vx_1x3v_tensor_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_dense_alpha_quad_vx_1x3v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vx_1x3v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_sparse_alpha_quad_vx_1x3v_tensor_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_ho_alpha_quad_vx_1x3v_tensor_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vy_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vy_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_ho_alpha_quad_vy_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_ho_alpha_quad_vy_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vy_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_ho_alpha_quad_vy_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vy_1x3v_tensor_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_vy_1x3v_tensor_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vy_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vy_1x3v_tensor_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vy_1x3v_tensor_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vy_1x3v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vy_1x3v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vy_1x3v_tensor_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vy_1x3v_tensor_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_vy_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_ho_alpha_quad_vy_1x3v_tensor_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_ho_alpha_quad_vy_1x3v_tensor_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_dense_alpha_quad_vy_1x3v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vy_1x3v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_sparse_alpha_quad_vy_1x3v_tensor_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_ho_alpha_quad_vy_1x3v_tensor_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vz_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vz_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_ho_alpha_quad_vz_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_ho_alpha_quad_vz_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vz_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_ho_alpha_quad_vz_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vz_1x3v_tensor_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_vz_1x3v_tensor_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vz_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vz_1x3v_tensor_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vz_1x3v_tensor_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vz_1x3v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vz_1x3v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vz_1x3v_tensor_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vz_1x3v_tensor_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_vz_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_ho_alpha_quad_vz_1x3v_tensor_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_ho_alpha_quad_vz_1x3v_tensor_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_dense_alpha_quad_vz_1x3v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vz_1x3v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_sparse_alpha_quad_vz_1x3v_tensor_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_ho_alpha_quad_vz_1x3v_tensor_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH double lax_flux_nodal_x_1x3v_tensor_p2(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_x_1x3v_tensor_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_alpha_quad_x_1x3v_tensor_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_x_1x3v_tensor_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vx_1x3v_tensor_p2(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vx_1x3v_tensor_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vx_1x3v_tensor_p2(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vx_1x3v_tensor_p2(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vx_1x3v_tensor_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vx_1x3v_tensor_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vx_1x3v_tensor_p2(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vx_1x3v_tensor_p2(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vy_1x3v_tensor_p2(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vy_1x3v_tensor_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vy_1x3v_tensor_p2(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vy_1x3v_tensor_p2(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vy_1x3v_tensor_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vy_1x3v_tensor_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vy_1x3v_tensor_p2(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vy_1x3v_tensor_p2(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vz_1x3v_tensor_p2(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vz_1x3v_tensor_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vz_1x3v_tensor_p2(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vz_1x3v_tensor_p2(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vz_1x3v_tensor_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vz_1x3v_tensor_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vz_1x3v_tensor_p2(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vz_1x3v_tensor_p2(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH double lax_flux_nodal_x_2x1v_tensor_p1(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_x_2x1v_tensor_p1(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_x_2x1v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_dense_ho_alpha_quad_x_2x1v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_x_2x1v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_x_2x1v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_y_2x1v_tensor_p1(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_y_2x1v_tensor_p1(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_y_2x1v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_dense_ho_alpha_quad_y_2x1v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_y_2x1v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_y_2x1v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vx_2x1v_tensor_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_vx_2x1v_tensor_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vx_2x1v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vx_2x1v_tensor_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vx_2x1v_tensor_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vx_2x1v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vx_2x1v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vx_2x1v_tensor_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_vx_2x1v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_ho_alpha_quad_vx_2x1v_tensor_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_ho_alpha_quad_vx_2x1v_tensor_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_dense_alpha_quad_vx_2x1v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vx_2x1v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_ho_alpha_quad_vx_2x1v_tensor_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH double lax_flux_nodal_x_2x1v_tensor_p2(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_x_2x1v_tensor_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_x_2x1v_tensor_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_y_2x1v_tensor_p2(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_y_2x1v_tensor_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_y_2x1v_tensor_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vx_2x1v_tensor_p2(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vx_2x1v_tensor_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vx_2x1v_tensor_p2(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vx_2x1v_tensor_p2(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vx_2x1v_tensor_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vx_2x1v_tensor_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vx_2x1v_tensor_p2(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH double lax_flux_nodal_x_2x1v_tensor_p3(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_x_2x1v_tensor_p3(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_x_2x1v_tensor_p3(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_y_2x1v_tensor_p3(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_y_2x1v_tensor_p3(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_y_2x1v_tensor_p3(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vx_2x1v_tensor_p3(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vx_2x1v_tensor_p3(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vx_2x1v_tensor_p3(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vx_2x1v_tensor_p3(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vx_2x1v_tensor_p3(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vx_2x1v_tensor_p3(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vx_2x1v_tensor_p3(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH double lax_flux_nodal_x_2x2v_tensor_p1(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_x_2x2v_tensor_p1(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_x_2x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_alpha_quad_x_2x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_dense_ho_alpha_quad_x_2x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_ho_alpha_quad_x_2x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_x_2x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_x_2x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_y_2x2v_tensor_p1(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_y_2x2v_tensor_p1(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_y_2x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_alpha_quad_y_2x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_dense_ho_alpha_quad_y_2x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_ho_alpha_quad_y_2x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_y_2x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_y_2x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vx_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vx_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_ho_alpha_quad_vx_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_ho_alpha_quad_vx_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vx_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_ho_alpha_quad_vx_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vx_2x2v_tensor_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_vx_2x2v_tensor_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vx_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vx_2x2v_tensor_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vx_2x2v_tensor_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vx_2x2v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vx_2x2v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vx_2x2v_tensor_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vx_2x2v_tensor_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_vx_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_ho_alpha_quad_vx_2x2v_tensor_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_ho_alpha_quad_vx_2x2v_tensor_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_dense_alpha_quad_vx_2x2v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vx_2x2v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_sparse_alpha_quad_vx_2x2v_tensor_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_ho_alpha_quad_vx_2x2v_tensor_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vy_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vy_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_ho_alpha_quad_vy_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_ho_alpha_quad_vy_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vy_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_ho_alpha_quad_vy_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vy_2x2v_tensor_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_vy_2x2v_tensor_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vy_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vy_2x2v_tensor_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vy_2x2v_tensor_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vy_2x2v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vy_2x2v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vy_2x2v_tensor_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vy_2x2v_tensor_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_vy_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_ho_alpha_quad_vy_2x2v_tensor_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_ho_alpha_quad_vy_2x2v_tensor_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_dense_alpha_quad_vy_2x2v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vy_2x2v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_sparse_alpha_quad_vy_2x2v_tensor_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_ho_alpha_quad_vy_2x2v_tensor_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH double lax_flux_nodal_x_2x2v_tensor_p2(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_x_2x2v_tensor_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_alpha_quad_x_2x2v_tensor_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_x_2x2v_tensor_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_y_2x2v_tensor_p2(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_y_2x2v_tensor_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_alpha_quad_y_2x2v_tensor_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_y_2x2v_tensor_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vx_2x2v_tensor_p2(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vx_2x2v_tensor_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vx_2x2v_tensor_p2(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vx_2x2v_tensor_p2(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vx_2x2v_tensor_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vx_2x2v_tensor_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vx_2x2v_tensor_p2(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vx_2x2v_tensor_p2(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vy_2x2v_tensor_p2(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vy_2x2v_tensor_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vy_2x2v_tensor_p2(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vy_2x2v_tensor_p2(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vy_2x2v_tensor_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vy_2x2v_tensor_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vy_2x2v_tensor_p2(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vy_2x2v_tensor_p2(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH double lax_flux_nodal_x_2x3v_tensor_p1(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_x_2x3v_tensor_p1(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_x_2x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_alpha_quad_x_2x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_dense_ho_alpha_quad_x_2x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_ho_alpha_quad_x_2x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_x_2x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_x_2x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_y_2x3v_tensor_p1(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_y_2x3v_tensor_p1(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_y_2x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_alpha_quad_y_2x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_dense_ho_alpha_quad_y_2x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_ho_alpha_quad_y_2x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_y_2x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_y_2x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vx_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vx_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_ho_alpha_quad_vx_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_ho_alpha_quad_vx_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vx_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_ho_alpha_quad_vx_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vx_2x3v_tensor_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_vx_2x3v_tensor_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vx_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vx_2x3v_tensor_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vx_2x3v_tensor_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vx_2x3v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vx_2x3v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vx_2x3v_tensor_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vx_2x3v_tensor_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_vx_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_ho_alpha_quad_vx_2x3v_tensor_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_ho_alpha_quad_vx_2x3v_tensor_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_dense_alpha_quad_vx_2x3v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vx_2x3v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_sparse_alpha_quad_vx_2x3v_tensor_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_ho_alpha_quad_vx_2x3v_tensor_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vy_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vy_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_ho_alpha_quad_vy_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_ho_alpha_quad_vy_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vy_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_ho_alpha_quad_vy_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vy_2x3v_tensor_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_vy_2x3v_tensor_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vy_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vy_2x3v_tensor_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vy_2x3v_tensor_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vy_2x3v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vy_2x3v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vy_2x3v_tensor_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vy_2x3v_tensor_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_vy_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_ho_alpha_quad_vy_2x3v_tensor_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_ho_alpha_quad_vy_2x3v_tensor_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_dense_alpha_quad_vy_2x3v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vy_2x3v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_sparse_alpha_quad_vy_2x3v_tensor_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_ho_alpha_quad_vy_2x3v_tensor_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vz_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vz_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_ho_alpha_quad_vz_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_ho_alpha_quad_vz_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vz_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_ho_alpha_quad_vz_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vz_2x3v_tensor_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_vz_2x3v_tensor_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vz_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vz_2x3v_tensor_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vz_2x3v_tensor_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vz_2x3v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vz_2x3v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vz_2x3v_tensor_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vz_2x3v_tensor_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_vz_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_ho_alpha_quad_vz_2x3v_tensor_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_ho_alpha_quad_vz_2x3v_tensor_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_dense_alpha_quad_vz_2x3v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vz_2x3v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_sparse_alpha_quad_vz_2x3v_tensor_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_ho_alpha_quad_vz_2x3v_tensor_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH double lax_flux_nodal_x_2x3v_tensor_p2(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_x_2x3v_tensor_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_alpha_quad_x_2x3v_tensor_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_x_2x3v_tensor_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_y_2x3v_tensor_p2(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_y_2x3v_tensor_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_alpha_quad_y_2x3v_tensor_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_y_2x3v_tensor_p2(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vx_2x3v_tensor_p2(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void phi_alpha_quad_vx_2x3v_tensor_p2(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vx_2x3v_tensor_p2(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vx_2x3v_tensor_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vx_2x3v_tensor_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vx_2x3v_tensor_p2(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vx_2x3v_tensor_p2(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vy_2x3v_tensor_p2(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void phi_alpha_quad_vy_2x3v_tensor_p2(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vy_2x3v_tensor_p2(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vy_2x3v_tensor_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vy_2x3v_tensor_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vy_2x3v_tensor_p2(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vy_2x3v_tensor_p2(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vz_2x3v_tensor_p2(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void phi_alpha_quad_vz_2x3v_tensor_p2(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vz_2x3v_tensor_p2(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vz_2x3v_tensor_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vz_2x3v_tensor_p2(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vz_2x3v_tensor_p2(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vz_2x3v_tensor_p2(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 

GKYL_CU_DH double lax_flux_nodal_x_3x3v_tensor_p1(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_x_3x3v_tensor_p1(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_x_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_alpha_quad_x_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_dense_ho_alpha_quad_x_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_ho_alpha_quad_x_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_x_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_x_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_y_3x3v_tensor_p1(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_y_3x3v_tensor_p1(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_y_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_alpha_quad_y_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_dense_ho_alpha_quad_y_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_ho_alpha_quad_y_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_y_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_y_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_z_3x3v_tensor_p1(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_z_3x3v_tensor_p1(const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_z_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_alpha_quad_z_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_dense_ho_alpha_quad_z_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_vel_sparse_ho_alpha_quad_z_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_alpha_quad_z_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_z_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
      const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
      const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vx_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vx_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_ho_alpha_quad_vx_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_ho_alpha_quad_vx_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vx_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_ho_alpha_quad_vx_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vx_3x3v_tensor_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_vx_3x3v_tensor_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vx_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vx_3x3v_tensor_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vx_3x3v_tensor_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vx_3x3v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vx_3x3v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vx_3x3v_tensor_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vx_3x3v_tensor_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_vx_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_ho_alpha_quad_vx_3x3v_tensor_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_ho_alpha_quad_vx_3x3v_tensor_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_dense_alpha_quad_vx_3x3v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vx_3x3v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_sparse_alpha_quad_vx_3x3v_tensor_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_ho_alpha_quad_vx_3x3v_tensor_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vy_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vy_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_ho_alpha_quad_vy_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_ho_alpha_quad_vy_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vy_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_ho_alpha_quad_vy_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vy_3x3v_tensor_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_vy_3x3v_tensor_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vy_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vy_3x3v_tensor_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vy_3x3v_tensor_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vy_3x3v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vy_3x3v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vy_3x3v_tensor_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vy_3x3v_tensor_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_vy_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_ho_alpha_quad_vy_3x3v_tensor_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_ho_alpha_quad_vy_3x3v_tensor_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_dense_alpha_quad_vy_3x3v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vy_3x3v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_sparse_alpha_quad_vy_3x3v_tensor_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_ho_alpha_quad_vy_3x3v_tensor_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vz_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vz_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_dense_ho_alpha_quad_vz_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_vel_sparse_ho_alpha_quad_vz_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vz_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void nc_hamil_phase_ho_alpha_quad_vz_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, 
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double lax_flux_nodal_vz_3x3v_tensor_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH double ho_lax_flux_nodal_vz_3x3v_tensor_p1(const double *dxv, 
      const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, const double *alpha_quad, 
      const double *f_l, const double *f_r, double* GKYL_RESTRICT vel_flux_surf) ; 
GKYL_CU_DH void hamil_phase_alpha_quad_vz_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_alpha_quad_vz_3x3v_tensor_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_alpha_quad_vz_3x3v_tensor_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vz_3x3v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vz_3x3v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vz_3x3v_tensor_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_alpha_quad_vz_3x3v_tensor_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_vz_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
        const double *hamil, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void phi_ho_alpha_quad_vz_3x3v_tensor_p1(const double *dxv, 
      const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void E_ho_alpha_quad_vz_3x3v_tensor_p1(const double *dxv, 
      const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_dense_alpha_quad_vz_3x3v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vz_3x3v_tensor_p1(const double *dxv, const double *jacob_vel,
      const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void B_ho_hamil_vel_sparse_alpha_quad_vz_3x3v_tensor_p1(const double *dxv, const double *jacob_vel,
        const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH void rad_ho_alpha_quad_vz_3x3v_tensor_p1(const double *dxv, 
      const double *rad, double* GKYL_RESTRICT alpha_quad); 

EXTERN_C_END 
