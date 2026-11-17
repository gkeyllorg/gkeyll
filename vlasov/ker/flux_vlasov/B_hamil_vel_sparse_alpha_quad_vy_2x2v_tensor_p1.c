#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x2v_tensor_p1.h> 
GKYL_CU_DH int B_hamil_vel_sparse_alpha_quad_vy_2x2v_tensor_p1_shared(int tid, int off, const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 1; 
  O += off*4; 
  I += off*3; 
  const double *Bz = &qmem[20]; 
  if (tid < 4) { 
    const int i = tid; 
    double Bz_quad = 0.0; 
    for (int a = 0; a < 4; ++a) Bz_quad += vst_2x2v_tensor_p1_conf_ev[i*4 + a]*Bz[a]; 
    O[0*4 + i] = -Bz_quad; 
  } 
  if (tid < 3) { 
    const int j = tid; 
    double dH_dvx = 0.0; 
    for (int s = 0; s < 5; ++s) { 
      const int b = vst_2x2v_tensor_p1_vel_sparse_idx[s]; 
      dH_dvx += vst_2x2v_tensor_p1_vel_dv0_v1[j*9 + b]*hamil[b]; 
    } 
    I[0*3 + j] = 2.0/(dxv[2]*jacob_vel_surf[0 + j])*dH_dvx; 
  } 
  return 1; 
} 

GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vy_2x2v_tensor_p1(const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[4]; 
  double I[3]; 
  for (int tid = 0; tid < 4; ++tid) B_hamil_vel_sparse_alpha_quad_vy_2x2v_tensor_p1_shared(tid, 0, dxv, jacob_vel_surf, hamil, qmem, O, I); 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 3; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 1; ++t) alpha += O[t*4 + i]*I[t*3 + j]; 
      alpha_quad[i*3 + j] += alpha; 
    } 
  } 
} 
