#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x2v_ser_p2.h> 
GKYL_CU_DH int phi_ho_alpha_quad_vy_2x2v_ser_p2_shared(int tid, int nthreads, const double *dxv, const double *jacob_pos, const double *phi,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 1; 
  const double dx11 = 2.0/dxv[1]; 
  const double jacob_cy_inv = 1.0/jacob_pos[3]; 
  for (int i = tid; i < 16; i += nthreads) { 
    double force_quad = 0.0; 
    for (int a = 0; a < 8; ++a) force_quad += vst_2x2v_ser_p2_ho_conf_dx1[i*8 + a]*phi[a]; 
    O[i] = -dx11*(force_quad*jacob_cy_inv); 
  } 
  for (int j = tid; j < 4; j += nthreads) { 
    I[j] = 1.0; 
  } 
  return 1; 
} 

GKYL_CU_DH void phi_ho_alpha_quad_vy_2x2v_ser_p2(const double *dxv, const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[16]; 
  double I[4]; 
  phi_ho_alpha_quad_vy_2x2v_ser_p2_shared(0, 1, dxv, jacob_pos, phi, O, I); 
  for (int i = 0; i < 16; ++i) { 
    for (int j = 0; j < 4; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 1; ++t) alpha += O[t*16 + i]*I[t*4 + j]; 
      alpha_quad[i*4 + j] += alpha; 
    } 
  } 
} 
