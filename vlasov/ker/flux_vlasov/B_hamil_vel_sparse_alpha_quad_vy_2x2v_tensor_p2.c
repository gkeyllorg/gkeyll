#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x2v_tensor_p2.h> 
GKYL_CU_DH int B_hamil_vel_sparse_alpha_quad_vy_2x2v_tensor_p2_shared(int tid, int nthreads, const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 1; 
  const double *Bz = &qmem[45]; 
  for (int i = tid; i < 16; i += nthreads) { 
    double Bz_quad = 0.0; 
    for (int a = 0; a < 9; ++a) Bz_quad += vst_2x2v_tensor_p2_conf_ev[i*9 + a]*Bz[a]; 
    O[0*16 + i] = -Bz_quad; 
  } 
  for (int j = tid; j < 4; j += nthreads) { 
    double dH_dvx = 0.0; 
    for (int s = 0; s < 5; ++s) { 
      const int b = vst_2x2v_tensor_p2_vel_sparse_idx[s]; 
      dH_dvx += vst_2x2v_tensor_p2_vel_dv0_v1[j*9 + b]*hamil[b]; 
    } 
    I[0*4 + j] = 2.0/(dxv[2]*jacob_vel_surf[0 + j])*dH_dvx; 
  } 
  return 1; 
} 

GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vy_2x2v_tensor_p2(const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[16]; 
  double I[4]; 
  B_hamil_vel_sparse_alpha_quad_vy_2x2v_tensor_p2_shared(0, 1, dxv, jacob_vel_surf, hamil, qmem, O, I); 
  for (int i = 0; i < 16; ++i) { 
    for (int j = 0; j < 4; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 1; ++t) alpha += O[t*16 + i]*I[t*4 + j]; 
      alpha_quad[i*4 + j] += alpha; 
    } 
  } 
} 
