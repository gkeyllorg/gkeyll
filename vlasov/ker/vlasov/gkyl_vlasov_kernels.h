#pragma once
#include <math.h>
#include <gkyl_util.h>
EXTERN_C_BEG

GKYL_CU_DH void vlasov_divide_Jv_1x1v_ser_p1(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_rescale_Jv_1x1v_ser_p1(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *vmap, 
    const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_vel_dense_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *vmap,
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_phase_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *vmap,
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_phase_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_phi_vol_1x1v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_E_vol_1x1v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_dense_vol_1x1v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_dense_vol_1x1v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_dense_vol_1x1v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_phase_vol_1x1v_ser_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_phase_vol_1x1v_ser_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_phase_vol_1x1v_ser_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_rad_vol_1x1v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *rad, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfx_1x1v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfx_1x1v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfx_1x1v_ser_p1(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfx_1x1v_ser_p1(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_1x1v_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_1x1v_ser_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_divide_Jv_1x1v_ser_p2(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_rescale_Jv_1x1v_ser_p2(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *vmap, 
    const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_vel_dense_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *vmap,
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_phase_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *vmap,
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_phase_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_phi_vol_1x1v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_E_vol_1x1v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_dense_vol_1x1v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_dense_vol_1x1v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_dense_vol_1x1v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_phase_vol_1x1v_ser_p2(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_phase_vol_1x1v_ser_p2(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_phase_vol_1x1v_ser_p2(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_rad_vol_1x1v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *rad, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfx_1x1v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfx_1x1v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfx_1x1v_ser_p2(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfx_1x1v_ser_p2(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfx_1x1v_ser_p2(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfx_1x1v_ser_p2(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_1x1v_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_1x1v_ser_p2(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfvx_1x1v_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfvx_1x1v_ser_p2(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_divide_Jv_1x1v_ser_p3(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_rescale_Jv_1x1v_ser_p3(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_1x1v_ser_p3(const double *w, const double *dxv, const double *vmap, 
    const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_phase_vol_1x1v_ser_p3(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_phi_vol_1x1v_ser_p3(const double *w, const double *dxv, 
    const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_E_vol_1x1v_ser_p3(const double *w, const double *dxv, 
    const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_dense_vol_1x1v_ser_p3(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_dense_vol_1x1v_ser_p3(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_dense_vol_1x1v_ser_p3(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_phase_vol_1x1v_ser_p3(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_phase_vol_1x1v_ser_p3(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_phase_vol_1x1v_ser_p3(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_rad_vol_1x1v_ser_p3(const double *w, const double *dxv, 
    const double *jacob_vel, const double *rad, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfx_1x1v_ser_p3(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfx_1x1v_ser_p3(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_1x1v_ser_p3(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_1x1v_ser_p3(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfvx_1x1v_ser_p3(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfvx_1x1v_ser_p3(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_divide_Jv_1x2v_ser_p1(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_rescale_Jv_1x2v_ser_p1(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap, 
    const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_vel_sparse_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_vel_dense_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap,
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_vel_sparse_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap,
        const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_phase_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap,
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_phase_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_phi_vol_1x2v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_E_vol_1x2v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_dense_vol_1x2v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_dense_vol_1x2v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_dense_vol_1x2v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_sparse_vol_1x2v_ser_p1(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_sparse_vol_1x2v_ser_p1(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_sparse_vol_1x2v_ser_p1(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_phase_vol_1x2v_ser_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_phase_vol_1x2v_ser_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_phase_vol_1x2v_ser_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_rad_vol_1x2v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *rad, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfx_1x2v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfx_1x2v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_boundary_surfx_1x2v_ser_p1(const double *w, const double *dxv, 
      const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_surfx_1x2v_ser_p1(const double *w, const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfx_1x2v_ser_p1(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfx_1x2v_ser_p1(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_1x2v_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_1x2v_ser_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvy_1x2v_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvy_1x2v_ser_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_divide_Jv_1x2v_ser_p2(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_rescale_Jv_1x2v_ser_p2(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *vmap, 
    const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_vel_sparse_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_vel_dense_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *vmap,
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_vel_sparse_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *vmap,
        const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_phase_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *vmap,
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_phase_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_phi_vol_1x2v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_E_vol_1x2v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_dense_vol_1x2v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_dense_vol_1x2v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_dense_vol_1x2v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_sparse_vol_1x2v_ser_p2(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_sparse_vol_1x2v_ser_p2(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_sparse_vol_1x2v_ser_p2(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_phase_vol_1x2v_ser_p2(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_phase_vol_1x2v_ser_p2(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_phase_vol_1x2v_ser_p2(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_rad_vol_1x2v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *rad, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfx_1x2v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfx_1x2v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_boundary_surfx_1x2v_ser_p2(const double *w, const double *dxv, 
      const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_surfx_1x2v_ser_p2(const double *w, const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfx_1x2v_ser_p2(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfx_1x2v_ser_p2(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfx_1x2v_ser_p2(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfx_1x2v_ser_p2(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_1x2v_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_1x2v_ser_p2(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfvx_1x2v_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfvx_1x2v_ser_p2(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvy_1x2v_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvy_1x2v_ser_p2(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfvy_1x2v_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfvy_1x2v_ser_p2(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_divide_Jv_1x3v_ser_p1(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_rescale_Jv_1x3v_ser_p1(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_1x3v_ser_p1(const double *w, const double *dxv, const double *vmap, 
    const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_vel_sparse_vol_1x3v_ser_p1(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_vel_dense_vol_1x3v_ser_p1(const double *w, const double *dxv, const double *vmap,
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_vel_sparse_vol_1x3v_ser_p1(const double *w, const double *dxv, const double *vmap,
        const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_phase_vol_1x3v_ser_p1(const double *w, const double *dxv, const double *vmap,
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_phase_vol_1x3v_ser_p1(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_phi_vol_1x3v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_E_vol_1x3v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_dense_vol_1x3v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_dense_vol_1x3v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_dense_vol_1x3v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_sparse_vol_1x3v_ser_p1(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_sparse_vol_1x3v_ser_p1(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_sparse_vol_1x3v_ser_p1(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_phase_vol_1x3v_ser_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_phase_vol_1x3v_ser_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_phase_vol_1x3v_ser_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_rad_vol_1x3v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *rad, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfx_1x3v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfx_1x3v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_boundary_surfx_1x3v_ser_p1(const double *w, const double *dxv, 
      const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_surfx_1x3v_ser_p1(const double *w, const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfx_1x3v_ser_p1(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfx_1x3v_ser_p1(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_1x3v_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_1x3v_ser_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvy_1x3v_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvy_1x3v_ser_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvz_1x3v_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvz_1x3v_ser_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_divide_Jv_1x3v_ser_p2(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_rescale_Jv_1x3v_ser_p2(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_1x3v_ser_p2(const double *w, const double *dxv, const double *vmap, 
    const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_vel_sparse_vol_1x3v_ser_p2(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_vel_dense_vol_1x3v_ser_p2(const double *w, const double *dxv, const double *vmap,
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_vel_sparse_vol_1x3v_ser_p2(const double *w, const double *dxv, const double *vmap,
        const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_phase_vol_1x3v_ser_p2(const double *w, const double *dxv, const double *vmap,
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_phase_vol_1x3v_ser_p2(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_phi_vol_1x3v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_E_vol_1x3v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_dense_vol_1x3v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_dense_vol_1x3v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_dense_vol_1x3v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_sparse_vol_1x3v_ser_p2(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_sparse_vol_1x3v_ser_p2(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_sparse_vol_1x3v_ser_p2(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_phase_vol_1x3v_ser_p2(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_phase_vol_1x3v_ser_p2(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_phase_vol_1x3v_ser_p2(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_rad_vol_1x3v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *rad, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfx_1x3v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfx_1x3v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_boundary_surfx_1x3v_ser_p2(const double *w, const double *dxv, 
      const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_surfx_1x3v_ser_p2(const double *w, const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfx_1x3v_ser_p2(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfx_1x3v_ser_p2(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfx_1x3v_ser_p2(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfx_1x3v_ser_p2(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_1x3v_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_1x3v_ser_p2(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfvx_1x3v_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfvx_1x3v_ser_p2(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvy_1x3v_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvy_1x3v_ser_p2(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfvy_1x3v_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfvy_1x3v_ser_p2(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvz_1x3v_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvz_1x3v_ser_p2(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfvz_1x3v_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfvz_1x3v_ser_p2(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_divide_Jv_2x1v_ser_p1(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_rescale_Jv_2x1v_ser_p1(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_2x1v_ser_p1(const double *w, const double *dxv, const double *vmap, 
    const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_phase_vol_2x1v_ser_p1(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_phi_vol_2x1v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_E_vol_2x1v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_dense_vol_2x1v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_dense_vol_2x1v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_dense_vol_2x1v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_phase_vol_2x1v_ser_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_phase_vol_2x1v_ser_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_phase_vol_2x1v_ser_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_rad_vol_2x1v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *rad, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfx_2x1v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfx_2x1v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfy_2x1v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfy_2x1v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_2x1v_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_2x1v_ser_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_divide_Jv_2x1v_ser_p2(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_rescale_Jv_2x1v_ser_p2(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_2x1v_ser_p2(const double *w, const double *dxv, const double *vmap, 
    const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_phase_vol_2x1v_ser_p2(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_phi_vol_2x1v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_E_vol_2x1v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_dense_vol_2x1v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_dense_vol_2x1v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_dense_vol_2x1v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_phase_vol_2x1v_ser_p2(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_phase_vol_2x1v_ser_p2(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_phase_vol_2x1v_ser_p2(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_rad_vol_2x1v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *rad, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfx_2x1v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfx_2x1v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfy_2x1v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfy_2x1v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_2x1v_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_2x1v_ser_p2(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfvx_2x1v_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfvx_2x1v_ser_p2(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_divide_Jv_2x1v_ser_p3(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_rescale_Jv_2x1v_ser_p3(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_2x1v_ser_p3(const double *w, const double *dxv, const double *vmap, 
    const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_phase_vol_2x1v_ser_p3(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_phi_vol_2x1v_ser_p3(const double *w, const double *dxv, 
    const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_E_vol_2x1v_ser_p3(const double *w, const double *dxv, 
    const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_dense_vol_2x1v_ser_p3(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_dense_vol_2x1v_ser_p3(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_dense_vol_2x1v_ser_p3(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_phase_vol_2x1v_ser_p3(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_phase_vol_2x1v_ser_p3(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_phase_vol_2x1v_ser_p3(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_rad_vol_2x1v_ser_p3(const double *w, const double *dxv, 
    const double *jacob_vel, const double *rad, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfx_2x1v_ser_p3(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfx_2x1v_ser_p3(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfy_2x1v_ser_p3(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfy_2x1v_ser_p3(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_2x1v_ser_p3(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_2x1v_ser_p3(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfvx_2x1v_ser_p3(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfvx_2x1v_ser_p3(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_divide_Jv_2x2v_ser_p1(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_rescale_Jv_2x2v_ser_p1(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap, 
    const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_vel_sparse_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_vel_dense_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap,
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_vel_sparse_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap,
        const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_phase_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap,
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_phase_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_phi_vol_2x2v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_E_vol_2x2v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_dense_vol_2x2v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_dense_vol_2x2v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_dense_vol_2x2v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_sparse_vol_2x2v_ser_p1(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_sparse_vol_2x2v_ser_p1(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_sparse_vol_2x2v_ser_p1(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_phase_vol_2x2v_ser_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_phase_vol_2x2v_ser_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_phase_vol_2x2v_ser_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_rad_vol_2x2v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *rad, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfx_2x2v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfx_2x2v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_boundary_surfx_2x2v_ser_p1(const double *w, const double *dxv, 
      const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_surfx_2x2v_ser_p1(const double *w, const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfx_2x2v_ser_p1(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfx_2x2v_ser_p1(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfy_2x2v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfy_2x2v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_boundary_surfy_2x2v_ser_p1(const double *w, const double *dxv, 
      const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_surfy_2x2v_ser_p1(const double *w, const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfy_2x2v_ser_p1(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfy_2x2v_ser_p1(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_2x2v_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_2x2v_ser_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvy_2x2v_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvy_2x2v_ser_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_divide_Jv_2x2v_ser_p2(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_rescale_Jv_2x2v_ser_p2(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_2x2v_ser_p2(const double *w, const double *dxv, const double *vmap, 
    const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_vel_sparse_vol_2x2v_ser_p2(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_vel_dense_vol_2x2v_ser_p2(const double *w, const double *dxv, const double *vmap,
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_vel_sparse_vol_2x2v_ser_p2(const double *w, const double *dxv, const double *vmap,
        const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_phase_vol_2x2v_ser_p2(const double *w, const double *dxv, const double *vmap,
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_phase_vol_2x2v_ser_p2(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_phi_vol_2x2v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_E_vol_2x2v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_dense_vol_2x2v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_dense_vol_2x2v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_dense_vol_2x2v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_sparse_vol_2x2v_ser_p2(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_sparse_vol_2x2v_ser_p2(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_sparse_vol_2x2v_ser_p2(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_phase_vol_2x2v_ser_p2(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_phase_vol_2x2v_ser_p2(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_phase_vol_2x2v_ser_p2(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_rad_vol_2x2v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *rad, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfx_2x2v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfx_2x2v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_boundary_surfx_2x2v_ser_p2(const double *w, const double *dxv, 
      const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_surfx_2x2v_ser_p2(const double *w, const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfx_2x2v_ser_p2(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfx_2x2v_ser_p2(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfx_2x2v_ser_p2(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfx_2x2v_ser_p2(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfy_2x2v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfy_2x2v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_boundary_surfy_2x2v_ser_p2(const double *w, const double *dxv, 
      const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_surfy_2x2v_ser_p2(const double *w, const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfy_2x2v_ser_p2(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfy_2x2v_ser_p2(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfy_2x2v_ser_p2(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfy_2x2v_ser_p2(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_2x2v_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_2x2v_ser_p2(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfvx_2x2v_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfvx_2x2v_ser_p2(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvy_2x2v_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvy_2x2v_ser_p2(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfvy_2x2v_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfvy_2x2v_ser_p2(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_divide_Jv_2x3v_ser_p1(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_rescale_Jv_2x3v_ser_p1(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_2x3v_ser_p1(const double *w, const double *dxv, const double *vmap, 
    const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_vel_sparse_vol_2x3v_ser_p1(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_vel_dense_vol_2x3v_ser_p1(const double *w, const double *dxv, const double *vmap,
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_vel_sparse_vol_2x3v_ser_p1(const double *w, const double *dxv, const double *vmap,
        const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_phase_vol_2x3v_ser_p1(const double *w, const double *dxv, const double *vmap,
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_phase_vol_2x3v_ser_p1(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_phi_vol_2x3v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_E_vol_2x3v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_dense_vol_2x3v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_dense_vol_2x3v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_dense_vol_2x3v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_sparse_vol_2x3v_ser_p1(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_sparse_vol_2x3v_ser_p1(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_sparse_vol_2x3v_ser_p1(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_phase_vol_2x3v_ser_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_phase_vol_2x3v_ser_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_phase_vol_2x3v_ser_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_rad_vol_2x3v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *rad, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfx_2x3v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfx_2x3v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_boundary_surfx_2x3v_ser_p1(const double *w, const double *dxv, 
      const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_surfx_2x3v_ser_p1(const double *w, const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfx_2x3v_ser_p1(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfx_2x3v_ser_p1(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfy_2x3v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfy_2x3v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_boundary_surfy_2x3v_ser_p1(const double *w, const double *dxv, 
      const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_surfy_2x3v_ser_p1(const double *w, const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfy_2x3v_ser_p1(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfy_2x3v_ser_p1(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_2x3v_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_2x3v_ser_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvy_2x3v_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvy_2x3v_ser_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvz_2x3v_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvz_2x3v_ser_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_divide_Jv_2x3v_ser_p2(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_rescale_Jv_2x3v_ser_p2(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_2x3v_ser_p2(const double *w, const double *dxv, const double *vmap, 
    const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_vel_sparse_vol_2x3v_ser_p2(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_vel_dense_vol_2x3v_ser_p2(const double *w, const double *dxv, const double *vmap,
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_vel_sparse_vol_2x3v_ser_p2(const double *w, const double *dxv, const double *vmap,
        const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_phi_vol_2x3v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_E_vol_2x3v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_dense_vol_2x3v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_dense_vol_2x3v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_dense_vol_2x3v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_sparse_vol_2x3v_ser_p2(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_sparse_vol_2x3v_ser_p2(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_sparse_vol_2x3v_ser_p2(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_phase_vol_2x3v_ser_p2(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_phase_vol_2x3v_ser_p2(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_phase_vol_2x3v_ser_p2(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_rad_vol_2x3v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *rad, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfx_2x3v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfx_2x3v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_boundary_surfx_2x3v_ser_p2(const double *w, const double *dxv, 
      const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_surfx_2x3v_ser_p2(const double *w, const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfx_2x3v_ser_p2(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfx_2x3v_ser_p2(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfx_2x3v_ser_p2(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfx_2x3v_ser_p2(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfy_2x3v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfy_2x3v_ser_p2(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_boundary_surfy_2x3v_ser_p2(const double *w, const double *dxv, 
      const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_surfy_2x3v_ser_p2(const double *w, const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfy_2x3v_ser_p2(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfy_2x3v_ser_p2(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfy_2x3v_ser_p2(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfy_2x3v_ser_p2(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_2x3v_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_2x3v_ser_p2(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfvx_2x3v_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfvx_2x3v_ser_p2(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvy_2x3v_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvy_2x3v_ser_p2(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfvy_2x3v_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfvy_2x3v_ser_p2(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvz_2x3v_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvz_2x3v_ser_p2(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfvz_2x3v_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfvz_2x3v_ser_p2(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_divide_Jv_3x3v_ser_p1(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_rescale_Jv_3x3v_ser_p1(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_3x3v_ser_p1(const double *w, const double *dxv, const double *vmap, 
    const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_vel_sparse_vol_3x3v_ser_p1(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_vel_dense_vol_3x3v_ser_p1(const double *w, const double *dxv, const double *vmap,
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_vel_sparse_vol_3x3v_ser_p1(const double *w, const double *dxv, const double *vmap,
        const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_phase_vol_3x3v_ser_p1(const double *w, const double *dxv, const double *vmap,
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_phase_vol_3x3v_ser_p1(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_phi_vol_3x3v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_E_vol_3x3v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_dense_vol_3x3v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_dense_vol_3x3v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_dense_vol_3x3v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_sparse_vol_3x3v_ser_p1(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_sparse_vol_3x3v_ser_p1(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_sparse_vol_3x3v_ser_p1(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_phase_vol_3x3v_ser_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_phase_vol_3x3v_ser_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_phase_vol_3x3v_ser_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_rad_vol_3x3v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *rad, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfx_3x3v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfx_3x3v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_boundary_surfx_3x3v_ser_p1(const double *w, const double *dxv, 
      const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_surfx_3x3v_ser_p1(const double *w, const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfx_3x3v_ser_p1(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfx_3x3v_ser_p1(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfy_3x3v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfy_3x3v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_boundary_surfy_3x3v_ser_p1(const double *w, const double *dxv, 
      const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_surfy_3x3v_ser_p1(const double *w, const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfy_3x3v_ser_p1(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfy_3x3v_ser_p1(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfz_3x3v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfz_3x3v_ser_p1(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_boundary_surfz_3x3v_ser_p1(const double *w, const double *dxv, 
      const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_surfz_3x3v_ser_p1(const double *w, const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfz_3x3v_ser_p1(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfz_3x3v_ser_p1(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_3x3v_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_3x3v_ser_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvy_3x3v_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvy_3x3v_ser_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvz_3x3v_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvz_3x3v_ser_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 

GKYL_CU_DH void hamil_default_1v_ser_p1(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv); 
GKYL_CU_DH void hamil_sr_1v_ser_p1(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv); 

GKYL_CU_DH void rad_compton_1v_ser_p1(const double *w, const double *dxv, const double *vmap, double t_cool, double p0, double* GKYL_RESTRICT rad); 
GKYL_CU_DH void rad_curvature_1v_ser_p1(const double *w, const double *dxv, const double *vmap, double t_cool, double p0, double* GKYL_RESTRICT rad); 

GKYL_CU_DH void hamil_default_1v_ser_p2(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv); 
GKYL_CU_DH void hamil_sr_1v_ser_p2(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv); 

GKYL_CU_DH void rad_compton_1v_ser_p2(const double *w, const double *dxv, const double *vmap, double t_cool, double p0, double* GKYL_RESTRICT rad); 
GKYL_CU_DH void rad_curvature_1v_ser_p2(const double *w, const double *dxv, const double *vmap, double t_cool, double p0, double* GKYL_RESTRICT rad); 

GKYL_CU_DH void hamil_default_1v_ser_p3(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv); 
GKYL_CU_DH void hamil_sr_1v_ser_p3(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv); 

GKYL_CU_DH void rad_compton_1v_ser_p3(const double *w, const double *dxv, const double *vmap, double t_cool, double p0, double* GKYL_RESTRICT rad); 
GKYL_CU_DH void rad_curvature_1v_ser_p3(const double *w, const double *dxv, const double *vmap, double t_cool, double p0, double* GKYL_RESTRICT rad); 

GKYL_CU_DH void hamil_default_2v_ser_p1(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv); 
GKYL_CU_DH void hamil_sr_2v_ser_p1(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv); 

GKYL_CU_DH void rad_compton_2v_ser_p1(const double *w, const double *dxv, const double *vmap, double t_cool, double p0, double* GKYL_RESTRICT rad); 
GKYL_CU_DH void rad_curvature_2v_ser_p1(const double *w, const double *dxv, const double *vmap, double t_cool, double p0, double* GKYL_RESTRICT rad); 

GKYL_CU_DH void hamil_default_2v_ser_p2(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv); 
GKYL_CU_DH void hamil_sr_2v_ser_p2(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv); 

GKYL_CU_DH void rad_compton_2v_ser_p2(const double *w, const double *dxv, const double *vmap, double t_cool, double p0, double* GKYL_RESTRICT rad); 
GKYL_CU_DH void rad_curvature_2v_ser_p2(const double *w, const double *dxv, const double *vmap, double t_cool, double p0, double* GKYL_RESTRICT rad); 

GKYL_CU_DH void hamil_default_3v_ser_p1(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv); 
GKYL_CU_DH void hamil_sr_3v_ser_p1(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv); 

GKYL_CU_DH void rad_compton_3v_ser_p1(const double *w, const double *dxv, const double *vmap, double t_cool, double p0, double* GKYL_RESTRICT rad); 
GKYL_CU_DH void rad_curvature_3v_ser_p1(const double *w, const double *dxv, const double *vmap, double t_cool, double p0, double* GKYL_RESTRICT rad); 

GKYL_CU_DH void hamil_default_3v_ser_p2(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv); 
GKYL_CU_DH void hamil_sr_3v_ser_p2(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv); 

GKYL_CU_DH void rad_compton_3v_ser_p2(const double *w, const double *dxv, const double *vmap, double t_cool, double p0, double* GKYL_RESTRICT rad); 
GKYL_CU_DH void rad_curvature_3v_ser_p2(const double *w, const double *dxv, const double *vmap, double t_cool, double p0, double* GKYL_RESTRICT rad); 

GKYL_CU_DH void vlasov_divide_Jv_1x1v_tensor_p1(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_rescale_Jv_1x1v_tensor_p1(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_1x1v_tensor_p1(const double *w, const double *dxv, const double *vmap, 
    const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_vel_dense_vol_1x1v_tensor_p1(const double *w, const double *dxv, const double *vmap,
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_phase_vol_1x1v_tensor_p1(const double *w, const double *dxv, const double *vmap,
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_phase_vol_1x1v_tensor_p1(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_phi_vol_1x1v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_E_vol_1x1v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_dense_vol_1x1v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_dense_vol_1x1v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_dense_vol_1x1v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_phase_vol_1x1v_tensor_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_phase_vol_1x1v_tensor_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_phase_vol_1x1v_tensor_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_rad_vol_1x1v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *rad, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfx_1x1v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfx_1x1v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfx_1x1v_tensor_p1(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfx_1x1v_tensor_p1(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfx_1x1v_tensor_p1(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfx_1x1v_tensor_p1(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_1x1v_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_1x1v_tensor_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfvx_1x1v_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfvx_1x1v_tensor_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_divide_Jv_1x1v_tensor_p2(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_rescale_Jv_1x1v_tensor_p2(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_1x1v_tensor_p2(const double *w, const double *dxv, const double *vmap, 
    const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_phase_vol_1x1v_tensor_p2(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_phi_vol_1x1v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_E_vol_1x1v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_dense_vol_1x1v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_dense_vol_1x1v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_dense_vol_1x1v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_phase_vol_1x1v_tensor_p2(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_phase_vol_1x1v_tensor_p2(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_phase_vol_1x1v_tensor_p2(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_rad_vol_1x1v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *rad, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfx_1x1v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfx_1x1v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_1x1v_tensor_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_1x1v_tensor_p2(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_divide_Jv_1x1v_tensor_p3(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_rescale_Jv_1x1v_tensor_p3(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_1x1v_tensor_p3(const double *w, const double *dxv, const double *vmap, 
    const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_phase_vol_1x1v_tensor_p3(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_phi_vol_1x1v_tensor_p3(const double *w, const double *dxv, 
    const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_E_vol_1x1v_tensor_p3(const double *w, const double *dxv, 
    const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_dense_vol_1x1v_tensor_p3(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_dense_vol_1x1v_tensor_p3(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_dense_vol_1x1v_tensor_p3(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_phase_vol_1x1v_tensor_p3(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_phase_vol_1x1v_tensor_p3(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_phase_vol_1x1v_tensor_p3(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_rad_vol_1x1v_tensor_p3(const double *w, const double *dxv, 
    const double *jacob_vel, const double *rad, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfx_1x1v_tensor_p3(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfx_1x1v_tensor_p3(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_1x1v_tensor_p3(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_1x1v_tensor_p3(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_divide_Jv_1x2v_tensor_p1(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_rescale_Jv_1x2v_tensor_p1(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, 
    const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_vel_sparse_vol_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_vel_dense_vol_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap,
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_vel_sparse_vol_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap,
        const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_phase_vol_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap,
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_phase_vol_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_phi_vol_1x2v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_E_vol_1x2v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_dense_vol_1x2v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_dense_vol_1x2v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_dense_vol_1x2v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_sparse_vol_1x2v_tensor_p1(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_sparse_vol_1x2v_tensor_p1(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_sparse_vol_1x2v_tensor_p1(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_phase_vol_1x2v_tensor_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_phase_vol_1x2v_tensor_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_phase_vol_1x2v_tensor_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_rad_vol_1x2v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *rad, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfx_1x2v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfx_1x2v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_boundary_surfx_1x2v_tensor_p1(const double *w, const double *dxv, 
      const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_surfx_1x2v_tensor_p1(const double *w, const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfx_1x2v_tensor_p1(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfx_1x2v_tensor_p1(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfx_1x2v_tensor_p1(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfx_1x2v_tensor_p1(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_1x2v_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_1x2v_tensor_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfvx_1x2v_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfvx_1x2v_tensor_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvy_1x2v_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvy_1x2v_tensor_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfvy_1x2v_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfvy_1x2v_tensor_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_divide_Jv_1x2v_tensor_p2(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_rescale_Jv_1x2v_tensor_p2(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_1x2v_tensor_p2(const double *w, const double *dxv, const double *vmap, 
    const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_vel_sparse_vol_1x2v_tensor_p2(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_phase_vol_1x2v_tensor_p2(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_phi_vol_1x2v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_E_vol_1x2v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_dense_vol_1x2v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_dense_vol_1x2v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_dense_vol_1x2v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_sparse_vol_1x2v_tensor_p2(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_sparse_vol_1x2v_tensor_p2(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_sparse_vol_1x2v_tensor_p2(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_phase_vol_1x2v_tensor_p2(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_phase_vol_1x2v_tensor_p2(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_phase_vol_1x2v_tensor_p2(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_rad_vol_1x2v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *rad, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfx_1x2v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfx_1x2v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_boundary_surfx_1x2v_tensor_p2(const double *w, const double *dxv, 
      const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_surfx_1x2v_tensor_p2(const double *w, const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_1x2v_tensor_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_1x2v_tensor_p2(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvy_1x2v_tensor_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvy_1x2v_tensor_p2(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_divide_Jv_1x3v_tensor_p1(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_rescale_Jv_1x3v_tensor_p1(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, 
    const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_vel_sparse_vol_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_vel_dense_vol_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap,
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_vel_sparse_vol_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap,
        const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_phase_vol_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap,
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_phase_vol_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_phi_vol_1x3v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_E_vol_1x3v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_dense_vol_1x3v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_dense_vol_1x3v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_dense_vol_1x3v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_sparse_vol_1x3v_tensor_p1(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_sparse_vol_1x3v_tensor_p1(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_sparse_vol_1x3v_tensor_p1(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_phase_vol_1x3v_tensor_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_phase_vol_1x3v_tensor_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_phase_vol_1x3v_tensor_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_rad_vol_1x3v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *rad, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfx_1x3v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfx_1x3v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_boundary_surfx_1x3v_tensor_p1(const double *w, const double *dxv, 
      const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_surfx_1x3v_tensor_p1(const double *w, const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfx_1x3v_tensor_p1(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfx_1x3v_tensor_p1(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfx_1x3v_tensor_p1(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfx_1x3v_tensor_p1(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_1x3v_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_1x3v_tensor_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfvx_1x3v_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfvx_1x3v_tensor_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvy_1x3v_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvy_1x3v_tensor_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfvy_1x3v_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfvy_1x3v_tensor_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvz_1x3v_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvz_1x3v_tensor_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfvz_1x3v_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfvz_1x3v_tensor_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_divide_Jv_1x3v_tensor_p2(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_rescale_Jv_1x3v_tensor_p2(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_1x3v_tensor_p2(const double *w, const double *dxv, const double *vmap, 
    const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_vel_sparse_vol_1x3v_tensor_p2(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_phase_vol_1x3v_tensor_p2(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_phi_vol_1x3v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_E_vol_1x3v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_dense_vol_1x3v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_dense_vol_1x3v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_dense_vol_1x3v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_sparse_vol_1x3v_tensor_p2(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_sparse_vol_1x3v_tensor_p2(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_sparse_vol_1x3v_tensor_p2(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_phase_vol_1x3v_tensor_p2(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_phase_vol_1x3v_tensor_p2(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_phase_vol_1x3v_tensor_p2(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_rad_vol_1x3v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *rad, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfx_1x3v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfx_1x3v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_boundary_surfx_1x3v_tensor_p2(const double *w, const double *dxv, 
      const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_surfx_1x3v_tensor_p2(const double *w, const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_1x3v_tensor_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_1x3v_tensor_p2(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvy_1x3v_tensor_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvy_1x3v_tensor_p2(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvz_1x3v_tensor_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvz_1x3v_tensor_p2(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_divide_Jv_2x1v_tensor_p1(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_rescale_Jv_2x1v_tensor_p1(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_2x1v_tensor_p1(const double *w, const double *dxv, const double *vmap, 
    const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_phase_vol_2x1v_tensor_p1(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_phi_vol_2x1v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_E_vol_2x1v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_dense_vol_2x1v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_dense_vol_2x1v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_dense_vol_2x1v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_phase_vol_2x1v_tensor_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_phase_vol_2x1v_tensor_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_phase_vol_2x1v_tensor_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_rad_vol_2x1v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *rad, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfx_2x1v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfx_2x1v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfy_2x1v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfy_2x1v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_2x1v_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_2x1v_tensor_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfvx_2x1v_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfvx_2x1v_tensor_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_divide_Jv_2x1v_tensor_p2(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_rescale_Jv_2x1v_tensor_p2(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_2x1v_tensor_p2(const double *w, const double *dxv, const double *vmap, 
    const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_phase_vol_2x1v_tensor_p2(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_phi_vol_2x1v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_E_vol_2x1v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_dense_vol_2x1v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_dense_vol_2x1v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_dense_vol_2x1v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_phase_vol_2x1v_tensor_p2(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_phase_vol_2x1v_tensor_p2(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_phase_vol_2x1v_tensor_p2(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_rad_vol_2x1v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *rad, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfx_2x1v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfx_2x1v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfy_2x1v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfy_2x1v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_2x1v_tensor_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_2x1v_tensor_p2(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_divide_Jv_2x1v_tensor_p3(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_rescale_Jv_2x1v_tensor_p3(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_2x1v_tensor_p3(const double *w, const double *dxv, const double *vmap, 
    const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_phase_vol_2x1v_tensor_p3(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_phi_vol_2x1v_tensor_p3(const double *w, const double *dxv, 
    const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_E_vol_2x1v_tensor_p3(const double *w, const double *dxv, 
    const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_dense_vol_2x1v_tensor_p3(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_dense_vol_2x1v_tensor_p3(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_dense_vol_2x1v_tensor_p3(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_phase_vol_2x1v_tensor_p3(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_phase_vol_2x1v_tensor_p3(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_phase_vol_2x1v_tensor_p3(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_rad_vol_2x1v_tensor_p3(const double *w, const double *dxv, 
    const double *jacob_vel, const double *rad, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfx_2x1v_tensor_p3(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfx_2x1v_tensor_p3(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfy_2x1v_tensor_p3(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfy_2x1v_tensor_p3(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_2x1v_tensor_p3(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_2x1v_tensor_p3(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_divide_Jv_2x2v_tensor_p1(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_rescale_Jv_2x2v_tensor_p1(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, 
    const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_vel_sparse_vol_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_vel_dense_vol_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap,
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_vel_sparse_vol_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap,
        const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_phase_vol_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap,
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_phase_vol_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_phi_vol_2x2v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_E_vol_2x2v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_dense_vol_2x2v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_dense_vol_2x2v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_dense_vol_2x2v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_sparse_vol_2x2v_tensor_p1(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_sparse_vol_2x2v_tensor_p1(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_sparse_vol_2x2v_tensor_p1(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_phase_vol_2x2v_tensor_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_phase_vol_2x2v_tensor_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_phase_vol_2x2v_tensor_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_rad_vol_2x2v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *rad, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfx_2x2v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfx_2x2v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_boundary_surfx_2x2v_tensor_p1(const double *w, const double *dxv, 
      const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_surfx_2x2v_tensor_p1(const double *w, const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfx_2x2v_tensor_p1(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfx_2x2v_tensor_p1(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfx_2x2v_tensor_p1(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfx_2x2v_tensor_p1(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfy_2x2v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfy_2x2v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_boundary_surfy_2x2v_tensor_p1(const double *w, const double *dxv, 
      const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_surfy_2x2v_tensor_p1(const double *w, const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfy_2x2v_tensor_p1(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfy_2x2v_tensor_p1(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfy_2x2v_tensor_p1(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfy_2x2v_tensor_p1(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_2x2v_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_2x2v_tensor_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfvx_2x2v_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfvx_2x2v_tensor_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvy_2x2v_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvy_2x2v_tensor_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfvy_2x2v_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfvy_2x2v_tensor_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_divide_Jv_2x2v_tensor_p2(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_rescale_Jv_2x2v_tensor_p2(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_2x2v_tensor_p2(const double *w, const double *dxv, const double *vmap, 
    const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_vel_sparse_vol_2x2v_tensor_p2(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_phase_vol_2x2v_tensor_p2(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_phi_vol_2x2v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_E_vol_2x2v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_dense_vol_2x2v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_dense_vol_2x2v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_dense_vol_2x2v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_sparse_vol_2x2v_tensor_p2(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_sparse_vol_2x2v_tensor_p2(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_sparse_vol_2x2v_tensor_p2(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_phase_vol_2x2v_tensor_p2(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_phase_vol_2x2v_tensor_p2(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_phase_vol_2x2v_tensor_p2(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_rad_vol_2x2v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *rad, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfx_2x2v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfx_2x2v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_boundary_surfx_2x2v_tensor_p2(const double *w, const double *dxv, 
      const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_surfx_2x2v_tensor_p2(const double *w, const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfy_2x2v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfy_2x2v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_boundary_surfy_2x2v_tensor_p2(const double *w, const double *dxv, 
      const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_surfy_2x2v_tensor_p2(const double *w, const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_2x2v_tensor_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_2x2v_tensor_p2(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvy_2x2v_tensor_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvy_2x2v_tensor_p2(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_divide_Jv_2x3v_tensor_p1(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_rescale_Jv_2x3v_tensor_p1(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, 
    const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_vel_sparse_vol_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_vel_dense_vol_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap,
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_vel_sparse_vol_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap,
        const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_phase_vol_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap,
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_phase_vol_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_phi_vol_2x3v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_E_vol_2x3v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_dense_vol_2x3v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_dense_vol_2x3v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_dense_vol_2x3v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_sparse_vol_2x3v_tensor_p1(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_sparse_vol_2x3v_tensor_p1(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_sparse_vol_2x3v_tensor_p1(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_phase_vol_2x3v_tensor_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_phase_vol_2x3v_tensor_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_phase_vol_2x3v_tensor_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_rad_vol_2x3v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *rad, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfx_2x3v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfx_2x3v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_boundary_surfx_2x3v_tensor_p1(const double *w, const double *dxv, 
      const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_surfx_2x3v_tensor_p1(const double *w, const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfx_2x3v_tensor_p1(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfx_2x3v_tensor_p1(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfx_2x3v_tensor_p1(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfx_2x3v_tensor_p1(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfy_2x3v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfy_2x3v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_boundary_surfy_2x3v_tensor_p1(const double *w, const double *dxv, 
      const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_surfy_2x3v_tensor_p1(const double *w, const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfy_2x3v_tensor_p1(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfy_2x3v_tensor_p1(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfy_2x3v_tensor_p1(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfy_2x3v_tensor_p1(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_2x3v_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_2x3v_tensor_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfvx_2x3v_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfvx_2x3v_tensor_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvy_2x3v_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvy_2x3v_tensor_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfvy_2x3v_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfvy_2x3v_tensor_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvz_2x3v_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvz_2x3v_tensor_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfvz_2x3v_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfvz_2x3v_tensor_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_divide_Jv_2x3v_tensor_p2(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_rescale_Jv_2x3v_tensor_p2(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_2x3v_tensor_p2(const double *w, const double *dxv, const double *vmap, 
    const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_vel_sparse_vol_2x3v_tensor_p2(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_phi_vol_2x3v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_E_vol_2x3v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_dense_vol_2x3v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_dense_vol_2x3v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_dense_vol_2x3v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_sparse_vol_2x3v_tensor_p2(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_sparse_vol_2x3v_tensor_p2(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_sparse_vol_2x3v_tensor_p2(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_phase_vol_2x3v_tensor_p2(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_phase_vol_2x3v_tensor_p2(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_phase_vol_2x3v_tensor_p2(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_rad_vol_2x3v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_vel, const double *rad, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfx_2x3v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfx_2x3v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_boundary_surfx_2x3v_tensor_p2(const double *w, const double *dxv, 
      const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_surfx_2x3v_tensor_p2(const double *w, const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfy_2x3v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfy_2x3v_tensor_p2(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_boundary_surfy_2x3v_tensor_p2(const double *w, const double *dxv, 
      const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_surfy_2x3v_tensor_p2(const double *w, const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_2x3v_tensor_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_2x3v_tensor_p2(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvy_2x3v_tensor_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvy_2x3v_tensor_p2(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvz_2x3v_tensor_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvz_2x3v_tensor_p2(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_divide_Jv_3x3v_tensor_p1(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_rescale_Jv_3x3v_tensor_p1(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J); 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, 
    const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_vel_sparse_vol_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_vel_dense_vol_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap,
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_vel_sparse_vol_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap,
        const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_phase_vol_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap,
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_hamil_phase_vol_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, 
      const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, const double *f_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_phi_vol_3x3v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_E_vol_3x3v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_dense_vol_3x3v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_dense_vol_3x3v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_dense_vol_3x3v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_vel_sparse_vol_3x3v_tensor_p1(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_vel_sparse_vol_3x3v_tensor_p1(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_vel_sparse_vol_3x3v_tensor_p1(const double *w, const double *dxv, 
      const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bx_hamil_phase_vol_3x3v_tensor_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_By_hamil_phase_vol_3x3v_tensor_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_Bz_hamil_phase_vol_3x3v_tensor_p1(const double *w, const double *dxv,
    const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_rad_vol_3x3v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_vel, const double *rad, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfx_3x3v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfx_3x3v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_boundary_surfx_3x3v_tensor_p1(const double *w, const double *dxv, 
      const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_surfx_3x3v_tensor_p1(const double *w, const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfx_3x3v_tensor_p1(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfx_3x3v_tensor_p1(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfx_3x3v_tensor_p1(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfx_3x3v_tensor_p1(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfy_3x3v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfy_3x3v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_boundary_surfy_3x3v_tensor_p1(const double *w, const double *dxv, 
      const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_surfy_3x3v_tensor_p1(const double *w, const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfy_3x3v_tensor_p1(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfy_3x3v_tensor_p1(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfy_3x3v_tensor_p1(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfy_3x3v_tensor_p1(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfz_3x3v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfz_3x3v_tensor_p1(const double *w, const double *dxv, 
    const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_boundary_surfz_3x3v_tensor_p1(const double *w, const double *dxv, 
      const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_hamil_vel_sparse_surfz_3x3v_tensor_p1(const double *w, const double *dxv, 
      const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfz_3x3v_tensor_p1(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfz_3x3v_tensor_p1(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfz_3x3v_tensor_p1(const double *w, const double *dxv,
        const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfz_3x3v_tensor_p1(const double *w, const double *dxv,
        const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvx_3x3v_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvx_3x3v_tensor_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfvx_3x3v_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfvx_3x3v_tensor_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvy_3x3v_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvy_3x3v_tensor_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfvy_3x3v_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfvy_3x3v_tensor_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_surfvz_3x3v_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_surfvz_3x3v_tensor_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_ho_surfvz_3x3v_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_boundary_ho_surfvz_3x3v_tensor_p1(const double *w, const double *dxv,
      const int edge, const double *flux, double* GKYL_RESTRICT out); 

GKYL_CU_DH void hamil_default_1v_tensor_p1(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv); 
GKYL_CU_DH void hamil_sr_1v_tensor_p1(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv); 

GKYL_CU_DH void rad_compton_1v_tensor_p1(const double *w, const double *dxv, const double *vmap, double t_cool, double p0, double* GKYL_RESTRICT rad); 
GKYL_CU_DH void rad_curvature_1v_tensor_p1(const double *w, const double *dxv, const double *vmap, double t_cool, double p0, double* GKYL_RESTRICT rad); 

GKYL_CU_DH void hamil_default_1v_tensor_p2(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv); 
GKYL_CU_DH void hamil_sr_1v_tensor_p2(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv); 

GKYL_CU_DH void rad_compton_1v_tensor_p2(const double *w, const double *dxv, const double *vmap, double t_cool, double p0, double* GKYL_RESTRICT rad); 
GKYL_CU_DH void rad_curvature_1v_tensor_p2(const double *w, const double *dxv, const double *vmap, double t_cool, double p0, double* GKYL_RESTRICT rad); 

GKYL_CU_DH void hamil_default_1v_tensor_p3(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv); 
GKYL_CU_DH void hamil_sr_1v_tensor_p3(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv); 

GKYL_CU_DH void rad_compton_1v_tensor_p3(const double *w, const double *dxv, const double *vmap, double t_cool, double p0, double* GKYL_RESTRICT rad); 
GKYL_CU_DH void rad_curvature_1v_tensor_p3(const double *w, const double *dxv, const double *vmap, double t_cool, double p0, double* GKYL_RESTRICT rad); 

GKYL_CU_DH void hamil_default_2v_tensor_p1(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv); 
GKYL_CU_DH void hamil_sr_2v_tensor_p1(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv); 

GKYL_CU_DH void rad_compton_2v_tensor_p1(const double *w, const double *dxv, const double *vmap, double t_cool, double p0, double* GKYL_RESTRICT rad); 
GKYL_CU_DH void rad_curvature_2v_tensor_p1(const double *w, const double *dxv, const double *vmap, double t_cool, double p0, double* GKYL_RESTRICT rad); 

GKYL_CU_DH void hamil_default_2v_tensor_p2(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv); 
GKYL_CU_DH void hamil_sr_2v_tensor_p2(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv); 

GKYL_CU_DH void rad_compton_2v_tensor_p2(const double *w, const double *dxv, const double *vmap, double t_cool, double p0, double* GKYL_RESTRICT rad); 
GKYL_CU_DH void rad_curvature_2v_tensor_p2(const double *w, const double *dxv, const double *vmap, double t_cool, double p0, double* GKYL_RESTRICT rad); 

GKYL_CU_DH void hamil_default_3v_tensor_p1(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv); 
GKYL_CU_DH void hamil_sr_3v_tensor_p1(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv); 

GKYL_CU_DH void rad_compton_3v_tensor_p1(const double *w, const double *dxv, const double *vmap, double t_cool, double p0, double* GKYL_RESTRICT rad); 
GKYL_CU_DH void rad_curvature_3v_tensor_p1(const double *w, const double *dxv, const double *vmap, double t_cool, double p0, double* GKYL_RESTRICT rad); 

GKYL_CU_DH void hamil_default_3v_tensor_p2(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv); 
GKYL_CU_DH void hamil_sr_3v_tensor_p2(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv); 

GKYL_CU_DH void rad_compton_3v_tensor_p2(const double *w, const double *dxv, const double *vmap, double t_cool, double p0, double* GKYL_RESTRICT rad); 
GKYL_CU_DH void rad_curvature_3v_tensor_p2(const double *w, const double *dxv, const double *vmap, double t_cool, double p0, double* GKYL_RESTRICT rad); 

EXTERN_C_END 
