#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_3x3v_tensor_p1.h> 
GKYL_CU_DH int phi_alpha_quad_vx_3x3v_tensor_p1_shared(int tid, int off, const double *dxv, const double *jacob_pos, const double *phi,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 1; 
  O += off*8; 
  I += off*9; 
  const double dx10 = 2.0/dxv[0]; 
  const double jacob_cx_inv = 1.0/jacob_pos[0]; 
  if (tid < 8) { 
    const int i = tid; 
    double force_quad = 0.0; 
    for (int a = 0; a < 8; ++a) force_quad += vst_3x3v_tensor_p1_conf_dx0[i*8 + a]*phi[a]; 
    O[i] = -dx10*(force_quad*jacob_cx_inv); 
  } 
  if (tid < 9) { 
    const int j = tid; 
    I[j] = 1.0; 
  } 
  return 1; 
} 

GKYL_CU_DH void phi_alpha_quad_vx_3x3v_tensor_p1(const double *dxv, const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[8]; 
  double I[9]; 
  for (int tid = 0; tid < 9; ++tid) phi_alpha_quad_vx_3x3v_tensor_p1_shared(tid, 0, dxv, jacob_pos, phi, O, I); 
  for (int i = 0; i < 8; ++i) { 
    for (int j = 0; j < 9; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 1; ++t) alpha += O[t*8 + i]*I[t*9 + j]; 
      alpha_quad[i*9 + j] += alpha; 
    } 
  } 
} 
