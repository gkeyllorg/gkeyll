#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x2v_tensor_p1.h> 
GKYL_CU_DH int nc_hamil_phase_alpha_quad_vx_2x2v_tensor_p1_shared(int tid, int nthreads, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 12; 
  const double dx10 = 2.0/dxv[0]; 
  const double dx11 = 2.0/dxv[1]; 
  const double jacob_cx_inv = 1.0/jacob_pos[0]; 
  const double jacob_cy_inv = 1.0/jacob_pos[2]; 
  const double *poisson_tensor_conf_x0 = &poisson_tensor_conf[0]; 
  const double *poisson_tensor_conf_x2 = &poisson_tensor_conf[8]; 
  const double *vmap_v0 = &vmap[0]; 
  const double *vmap_v1 = &vmap[4]; 
  const double dv11 = 2.0/dxv[3]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[16]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[20]; 
  for (int i = tid; i < 4; i += nthreads) { 
    double px0_q = 0.0; 
    double px1_q = 0.0; 
    for (int a = 0; a < 4; ++a) { 
      px0_q += vst_2x2v_tensor_p1_conf_ev[i*4 + a]*poisson_tensor_conf_x0[a]; 
      px1_q += vst_2x2v_tensor_p1_conf_ev[i*4 + a]*poisson_tensor_conf_x2[a]; 
    } 
    for (int a = 0; a < 4; ++a) O[a*4 + i] = -(px0_q*vst_2x2v_tensor_p1_ph_v0_CmDx0[i*4 + a]*dx10*jacob_cx_inv + px1_q*vst_2x2v_tensor_p1_ph_v0_CmDx1[i*4 + a]*dx11*jacob_cy_inv); 
    double p0_q = 0.0; 
    double p1_q = 0.0; 
    for (int a = 0; a < 4; ++a) { 
      p0_q += vst_2x2v_tensor_p1_conf_ev[i*4 + a]*poisson_tensor_conf_0[a]; 
      p1_q += vst_2x2v_tensor_p1_conf_ev[i*4 + a]*poisson_tensor_conf_1[a]; 
    } 
    for (int a = 0; a < 4; ++a) O[(4 + a)*4 + i] = p0_q*vst_2x2v_tensor_p1_ph_v0_Cm[i*4 + a]; 
    for (int a = 0; a < 4; ++a) O[(8 + a)*4 + i] = p1_q*vst_2x2v_tensor_p1_ph_v0_Cm[i*4 + a]; 
  } 
  for (int j = tid; j < 3; j += nthreads) { 
    double G[4]; 
    double Gd1[4]; 
    for (int a = 0; a < 4; ++a) { G[a] = 0.0; Gd1[a] = 0.0; } 
    for (int k = 0; k < 36; ++k) { 
      const int a = vst_2x2v_tensor_p1_ph_v0_cmap[k]; 
      G[a] += vst_2x2v_tensor_p1_ph_v0_V[j*9 + vst_2x2v_tensor_p1_ph_v0_vrmap[k]]*(vst_2x2v_tensor_p1_ph_v0_coefr[k]*hamil[k]); 
      Gd1[a] += vst_2x2v_tensor_p1_ph_v0_Vd1[j*7 + vst_2x2v_tensor_p1_ph_v0_vrd1map[k]]*(vst_2x2v_tensor_p1_ph_v0_dcoefr1[k]*hamil[k]); 
    } 
    for (int a = 0; a < 4; ++a) I[a*3 + j] = G[a]; 
    const double jacob_vy_inv = 1.0/jacob_vel_surf[3 + j]; 
    const double vt1 = -(1.8708286933869707*vmap_v0[3])+1.5811388300841895*vmap_v0[2]-1.224744871391589*vmap_v0[1]+0.7071067811865475*vmap_v0[0]; 
    const double xn1 = vst_2x2v_tensor_p1_vel_nodes_v0[j*1 + 0]; 
    const double xn1_sq = xn1*xn1; 
    const double vt2 = 4.677071733467426*vmap_v1[3]*xn1*xn1_sq+2.371708245126284*vmap_v1[2]*xn1_sq-2.806243040080455*vmap_v1[3]*xn1+1.224744871391589*vmap_v1[1]*xn1-0.7905694150420947*vmap_v1[2]+0.7071067811865475*vmap_v1[0]; 
    for (int a = 0; a < 4; ++a) I[(4 + a)*3 + j] = vt1*Gd1[a]*dv11*jacob_vy_inv; 
    for (int a = 0; a < 4; ++a) I[(8 + a)*3 + j] = vt2*Gd1[a]*dv11*jacob_vy_inv; 
  } 
  return 12; 
} 

GKYL_CU_DH void nc_hamil_phase_alpha_quad_vx_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[48]; 
  double I[36]; 
  nc_hamil_phase_alpha_quad_vx_2x2v_tensor_p1_shared(0, 1, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 3; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 12; ++t) alpha += O[t*4 + i]*I[t*3 + j]; 
      alpha_quad[i*3 + j] += alpha; 
    } 
  } 
} 
