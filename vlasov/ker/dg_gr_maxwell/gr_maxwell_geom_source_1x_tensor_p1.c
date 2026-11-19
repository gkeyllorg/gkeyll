#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH void gr_maxwell_geom_source_1x_tensor_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, 
      const double *geom_factor_con, const double *fields, double* GKYL_RESTRICT out) 
{ 
  const double chi = meq->chi, gamma = meq->gamma; 
  const double K_phi = meq->K_phi, K_psi = meq->K_psi; 
  double dx0 = 2.0/dx[0]; 
  
  double Jphi_n[2] = {0.0};
  double Jpsi_n[2] = {0.0};
  double geom_sourceD1_n[2] = {0.0};
  double geom_sourceD2_n[2] = {0.0};
  double geom_sourceD3_n[2] = {0.0};
  double geom_sourceB1_n[2] = {0.0};
  double geom_sourceB2_n[2] = {0.0};
  double geom_sourceB3_n[2] = {0.0};
  double geom_sourcephi_n[2] = {0.0};
  double geom_sourcepsi_n[2] = {0.0};
  
  const double *Jphi = &fields[12]; 
  const double *Jpsi = &fields[14]; 
  
  double *outJDx = &out[0]; 
  double *outJDy = &out[2]; 
  double *outJDz = &out[4]; 
  double *outJBx = &out[6]; 
  double *outJBy = &out[8]; 
  double *outJBz = &out[10]; 
  double *outJphi = &out[12]; 
  double *outJpsi = &out[14]; 
  
  Jphi_n[0] = 0.7071067811865475*Jphi[0]-0.7071067811865475*Jphi[1];
  Jpsi_n[0] = 0.7071067811865475*Jpsi[0]-0.7071067811865475*Jpsi[1];
  
  Jphi_n[1] = 0.7071067811865475*Jphi[1]+0.7071067811865475*Jphi[0];
  Jpsi_n[1] = 0.7071067811865475*Jpsi[1]+0.7071067811865475*Jpsi[0];
  
  const double *geom_factor_con_x = &geom_factor_con[0]; 
  const double *geom_factor_con_y = &geom_factor_con[2]; 
  const double *geom_factor_con_z = &geom_factor_con[4]; 
  
  for (int i=0; i<2; ++i) {
    geom_sourceD1_n[i] = chi * Jphi_n[i] * geom_factor_con_x[i]; 
    geom_sourceD2_n[i] = chi * Jphi_n[i] * geom_factor_con_y[i]; 
    geom_sourceD3_n[i] = chi * Jphi_n[i] * geom_factor_con_z[i]; 
    geom_sourceB1_n[i] = gamma * Jpsi_n[i] * geom_factor_con_x[i]; 
    geom_sourceB2_n[i] = gamma * Jpsi_n[i] * geom_factor_con_y[i]; 
    geom_sourceB3_n[i] = gamma * Jpsi_n[i] * geom_factor_con_z[i]; 
    geom_sourcephi_n[i] = - K_phi * Jphi_n[i]; 
    geom_sourcepsi_n[i] = - K_psi * Jpsi_n[i]; 
  }
  
  outJDx[0] += 0.7071067811865475*geom_sourceD1_n[1]+0.7071067811865475*geom_sourceD1_n[0]; 
  outJDx[1] += 0.7071067811865475*geom_sourceD1_n[1]-0.7071067811865475*geom_sourceD1_n[0]; 

  outJDy[0] += 0.7071067811865475*geom_sourceD2_n[1]+0.7071067811865475*geom_sourceD2_n[0]; 
  outJDy[1] += 0.7071067811865475*geom_sourceD2_n[1]-0.7071067811865475*geom_sourceD2_n[0]; 

  outJDz[0] += 0.7071067811865475*geom_sourceD3_n[1]+0.7071067811865475*geom_sourceD3_n[0]; 
  outJDz[1] += 0.7071067811865475*geom_sourceD3_n[1]-0.7071067811865475*geom_sourceD3_n[0]; 

  outJBx[0] += 0.7071067811865475*geom_sourceB1_n[1]+0.7071067811865475*geom_sourceB1_n[0]; 
  outJBx[1] += 0.7071067811865475*geom_sourceB1_n[1]-0.7071067811865475*geom_sourceB1_n[0]; 

  outJBy[0] += 0.7071067811865475*geom_sourceB2_n[1]+0.7071067811865475*geom_sourceB2_n[0]; 
  outJBy[1] += 0.7071067811865475*geom_sourceB2_n[1]-0.7071067811865475*geom_sourceB2_n[0]; 

  outJBz[0] += 0.7071067811865475*geom_sourceB3_n[1]+0.7071067811865475*geom_sourceB3_n[0]; 
  outJBz[1] += 0.7071067811865475*geom_sourceB3_n[1]-0.7071067811865475*geom_sourceB3_n[0]; 

  outJphi[0] += 0.7071067811865475*geom_sourcephi_n[1]+0.7071067811865475*geom_sourcephi_n[0]; 
  outJphi[1] += 0.7071067811865475*geom_sourcephi_n[1]-0.7071067811865475*geom_sourcephi_n[0]; 

  outJpsi[0] += 0.7071067811865475*geom_sourcepsi_n[1]+0.7071067811865475*geom_sourcepsi_n[0]; 
  outJpsi[1] += 0.7071067811865475*geom_sourcepsi_n[1]-0.7071067811865475*geom_sourcepsi_n[0]; 

} 
