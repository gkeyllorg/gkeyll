#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_3x3v_tensor_p1.h> 
GKYL_CU_DH int phi_ho_alpha_quad_vz_3x3v_tensor_p1_shared(int tid, int off, const double *dxv, const double *jacob_pos, const double *phi,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 1; 
  O += off*8; 
  I += off*16; 
  const double dx12 = 2.0/dxv[2]; 
  const double jacob_cz_inv = 1.0/jacob_pos[4]; 
  if (tid < 8) { 
    const int i = tid; 
    double force_quad = 0.0; 
    for (int a = 0; a < 8; ++a) force_quad += vst_3x3v_tensor_p1_ho_conf_dx2[i*8 + a]*phi[a]; 
    O[i] = -dx12*(force_quad*jacob_cz_inv); 
  } 
  if (tid < 16) { 
    const int j = tid; 
    I[j] = 1.0; 
  } 
  return 1; 
} 

GKYL_CU_DH void phi_ho_alpha_quad_vz_3x3v_tensor_p1(const double *dxv, const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[8]; 
  double I[16]; 
  for (int tid = 0; tid < 16; ++tid) phi_ho_alpha_quad_vz_3x3v_tensor_p1_shared(tid, 0, dxv, jacob_pos, phi, O, I); 
  for (int i = 0; i < 8; ++i) { 
    for (int j = 0; j < 16; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 1; ++t) alpha += O[t*8 + i]*I[t*16 + j]; 
      alpha_quad[i*16 + j] += alpha; 
    } 
  } 
} 
