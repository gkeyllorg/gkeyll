#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x2v_ser_p2.h> 
GKYL_CU_DH int phi_alpha_quad_vx_1x2v_ser_p2_shared(int tid, int nthreads, const double *dxv, const double *jacob_pos, const double *phi,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 1; 
  const double dx10 = 2.0/dxv[0]; 
  const double jacob_cx_inv = 1.0/jacob_pos[0]; 
  for (int i = tid; i < 3; i += nthreads) { 
    double force_quad = 0.0; 
    for (int a = 0; a < 3; ++a) force_quad += vst_1x2v_ser_p2_conf_dx0[i*3 + a]*phi[a]; 
    O[i] = -dx10*(force_quad*jacob_cx_inv); 
  } 
  for (int j = tid; j < 3; j += nthreads) { 
    I[j] = 1.0; 
  } 
  return 1; 
} 

GKYL_CU_DH void phi_alpha_quad_vx_1x2v_ser_p2(const double *dxv, const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[3]; 
  double I[3]; 
  phi_alpha_quad_vx_1x2v_ser_p2_shared(0, 1, dxv, jacob_pos, phi, O, I); 
  for (int i = 0; i < 3; ++i) { 
    for (int j = 0; j < 3; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 1; ++t) alpha += O[t*3 + i]*I[t*3 + j]; 
      alpha_quad[i*3 + j] += alpha; 
    } 
  } 
} 
