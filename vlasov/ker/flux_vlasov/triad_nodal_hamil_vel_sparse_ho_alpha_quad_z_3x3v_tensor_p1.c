#include <gkyl_vlasov_triad_nodal_kernels.h> 
#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_3x3v_tensor_p1.h> 
GKYL_CU_DH int triad_nodal_hamil_vel_sparse_ho_alpha_quad_z_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 3; 
  O += off*4; 
  I += off*64; 
  const double dv10 = 2.0/dxv[3]; 
  const double dv11 = 2.0/dxv[4]; 
  const double dv12 = 2.0/dxv[5]; 
  const double *pt_face = &poisson_tensor_conf[(4 + (hamil_pt_edge == -1 ? 0 : 1))*12]; 
  if (tid < 4) { 
    const int i = tid; 
    O[0*4 + i] = pt_face[i*3 + 0]; 
    O[1*4 + i] = pt_face[i*3 + 1]; 
    O[2*4 + i] = pt_face[i*3 + 2]; 
  } 
  if (tid < 64) { 
    const int j = tid; 
    double dH_dv0 = 0.0; 
    for (int s = 0; s < 7; ++s) { 
      const int b = vst_3x3v_tensor_p1_ho_vel_sparse_idx[s]; 
      dH_dv0 += vst_3x3v_tensor_p1_ho_vel_vol_dv0[j*27 + b]*hamil[b]; 
    } 
    I[0*64 + j] = dH_dv0*dv10*(1.0/jacob_vel_surf[0 + j/16]); 
    double dH_dv1 = 0.0; 
    for (int s = 0; s < 7; ++s) { 
      const int b = vst_3x3v_tensor_p1_ho_vel_sparse_idx[s]; 
      dH_dv1 += vst_3x3v_tensor_p1_ho_vel_vol_dv1[j*27 + b]*hamil[b]; 
    } 
    I[1*64 + j] = dH_dv1*dv11*(1.0/jacob_vel_surf[4 + j/4%4]); 
    double dH_dv2 = 0.0; 
    for (int s = 0; s < 7; ++s) { 
      const int b = vst_3x3v_tensor_p1_ho_vel_sparse_idx[s]; 
      dH_dv2 += vst_3x3v_tensor_p1_ho_vel_vol_dv2[j*27 + b]*hamil[b]; 
    } 
    I[2*64 + j] = dH_dv2*dv12*(1.0/jacob_vel_surf[8 + j%4]); 
  } 
  return 3; 
} 

GKYL_CU_DH void triad_nodal_hamil_vel_sparse_ho_alpha_quad_z_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[12]; 
  double I[192]; 
  for (int tid = 0; tid < 64; ++tid) triad_nodal_hamil_vel_sparse_ho_alpha_quad_z_3x3v_tensor_p1_shared(tid, 0, w, dxv, hamil_pt_edge, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 64; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 3; ++t) alpha += O[t*4 + i]*I[t*64 + j]; 
      alpha_quad[i*64 + j] += alpha; 
    } 
  } 
} 
