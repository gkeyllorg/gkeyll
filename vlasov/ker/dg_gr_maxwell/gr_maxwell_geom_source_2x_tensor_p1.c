#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH void gr_maxwell_geom_source_2x_tensor_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, 
      const double *geom_factor_con, const double *fields, double* GKYL_RESTRICT out) 
{ 
  const double chi = meq->chi, gamma = meq->gamma; 
  const double K_phi = meq->K_phi, K_psi = meq->K_psi; 
  double dx0 = 2.0/dx[0]; 
  double dx1 = 2.0/dx[1]; 
  
  double Jphi_n[4] = {0.0};
  double Jpsi_n[4] = {0.0};
  double geom_sourceD1_n[4] = {0.0};
  double geom_sourceD2_n[4] = {0.0};
  double geom_sourceD3_n[4] = {0.0};
  double geom_sourceB1_n[4] = {0.0};
  double geom_sourceB2_n[4] = {0.0};
  double geom_sourceB3_n[4] = {0.0};
  double geom_sourcephi_n[4] = {0.0};
  double geom_sourcepsi_n[4] = {0.0};
  
  const double *Jphi = &fields[24]; 
  const double *Jpsi = &fields[28]; 
  
  double *outJDx = &out[0]; 
  double *outJDy = &out[4]; 
  double *outJDz = &out[8]; 
  double *outJBx = &out[12]; 
  double *outJBy = &out[16]; 
  double *outJBz = &out[20]; 
  double *outJphi = &out[24]; 
  double *outJpsi = &out[28]; 
  
  Jphi_n[0] = 0.5*Jphi[3]-0.5*Jphi[2]-0.5*Jphi[1]+0.5*Jphi[0];
  Jpsi_n[0] = 0.5*Jpsi[3]-0.5*Jpsi[2]-0.5*Jpsi[1]+0.5*Jpsi[0];
  
  Jphi_n[1] = -(0.5*Jphi[3])+0.5*Jphi[2]-0.5*Jphi[1]+0.5*Jphi[0];
  Jpsi_n[1] = -(0.5*Jpsi[3])+0.5*Jpsi[2]-0.5*Jpsi[1]+0.5*Jpsi[0];
  
  Jphi_n[2] = -(0.5*Jphi[3])-0.5*Jphi[2]+0.5*Jphi[1]+0.5*Jphi[0];
  Jpsi_n[2] = -(0.5*Jpsi[3])-0.5*Jpsi[2]+0.5*Jpsi[1]+0.5*Jpsi[0];
  
  Jphi_n[3] = 0.5*Jphi[3]+0.5*Jphi[2]+0.5*Jphi[1]+0.5*Jphi[0];
  Jpsi_n[3] = 0.5*Jpsi[3]+0.5*Jpsi[2]+0.5*Jpsi[1]+0.5*Jpsi[0];
  
  const double *geom_factor_con_x = &geom_factor_con[0]; 
  const double *geom_factor_con_y = &geom_factor_con[4]; 
  const double *geom_factor_con_z = &geom_factor_con[8]; 
  
  for (int i=0; i<4; ++i) {
    geom_sourceD1_n[i] = chi * Jphi_n[i] * geom_factor_con_x[i]; 
    geom_sourceD2_n[i] = chi * Jphi_n[i] * geom_factor_con_y[i]; 
    geom_sourceD3_n[i] = chi * Jphi_n[i] * geom_factor_con_z[i]; 
    geom_sourceB1_n[i] = gamma * Jpsi_n[i] * geom_factor_con_x[i]; 
    geom_sourceB2_n[i] = gamma * Jpsi_n[i] * geom_factor_con_y[i]; 
    geom_sourceB3_n[i] = gamma * Jpsi_n[i] * geom_factor_con_z[i]; 
    geom_sourcephi_n[i] = - K_phi * Jphi_n[i]; 
    geom_sourcepsi_n[i] = - K_psi * Jpsi_n[i]; 
  }
  
  outJDx[0] += 0.5*geom_sourceD1_n[3]+0.5*geom_sourceD1_n[2]+0.5*geom_sourceD1_n[1]+0.5*geom_sourceD1_n[0]; 
  outJDx[1] += 0.5*geom_sourceD1_n[3]+0.5*geom_sourceD1_n[2]-0.5*geom_sourceD1_n[1]-0.5*geom_sourceD1_n[0]; 
  outJDx[2] += 0.5*geom_sourceD1_n[3]-0.5*geom_sourceD1_n[2]+0.5*geom_sourceD1_n[1]-0.5*geom_sourceD1_n[0]; 
  outJDx[3] += 0.5*geom_sourceD1_n[3]-0.5*geom_sourceD1_n[2]-0.5*geom_sourceD1_n[1]+0.5*geom_sourceD1_n[0]; 

  outJDy[0] += 0.5*geom_sourceD2_n[3]+0.5*geom_sourceD2_n[2]+0.5*geom_sourceD2_n[1]+0.5*geom_sourceD2_n[0]; 
  outJDy[1] += 0.5*geom_sourceD2_n[3]+0.5*geom_sourceD2_n[2]-0.5*geom_sourceD2_n[1]-0.5*geom_sourceD2_n[0]; 
  outJDy[2] += 0.5*geom_sourceD2_n[3]-0.5*geom_sourceD2_n[2]+0.5*geom_sourceD2_n[1]-0.5*geom_sourceD2_n[0]; 
  outJDy[3] += 0.5*geom_sourceD2_n[3]-0.5*geom_sourceD2_n[2]-0.5*geom_sourceD2_n[1]+0.5*geom_sourceD2_n[0]; 

  outJDz[0] += 0.5*geom_sourceD3_n[3]+0.5*geom_sourceD3_n[2]+0.5*geom_sourceD3_n[1]+0.5*geom_sourceD3_n[0]; 
  outJDz[1] += 0.5*geom_sourceD3_n[3]+0.5*geom_sourceD3_n[2]-0.5*geom_sourceD3_n[1]-0.5*geom_sourceD3_n[0]; 
  outJDz[2] += 0.5*geom_sourceD3_n[3]-0.5*geom_sourceD3_n[2]+0.5*geom_sourceD3_n[1]-0.5*geom_sourceD3_n[0]; 
  outJDz[3] += 0.5*geom_sourceD3_n[3]-0.5*geom_sourceD3_n[2]-0.5*geom_sourceD3_n[1]+0.5*geom_sourceD3_n[0]; 

  outJBx[0] += 0.5*geom_sourceB1_n[3]+0.5*geom_sourceB1_n[2]+0.5*geom_sourceB1_n[1]+0.5*geom_sourceB1_n[0]; 
  outJBx[1] += 0.5*geom_sourceB1_n[3]+0.5*geom_sourceB1_n[2]-0.5*geom_sourceB1_n[1]-0.5*geom_sourceB1_n[0]; 
  outJBx[2] += 0.5*geom_sourceB1_n[3]-0.5*geom_sourceB1_n[2]+0.5*geom_sourceB1_n[1]-0.5*geom_sourceB1_n[0]; 
  outJBx[3] += 0.5*geom_sourceB1_n[3]-0.5*geom_sourceB1_n[2]-0.5*geom_sourceB1_n[1]+0.5*geom_sourceB1_n[0]; 

  outJBy[0] += 0.5*geom_sourceB2_n[3]+0.5*geom_sourceB2_n[2]+0.5*geom_sourceB2_n[1]+0.5*geom_sourceB2_n[0]; 
  outJBy[1] += 0.5*geom_sourceB2_n[3]+0.5*geom_sourceB2_n[2]-0.5*geom_sourceB2_n[1]-0.5*geom_sourceB2_n[0]; 
  outJBy[2] += 0.5*geom_sourceB2_n[3]-0.5*geom_sourceB2_n[2]+0.5*geom_sourceB2_n[1]-0.5*geom_sourceB2_n[0]; 
  outJBy[3] += 0.5*geom_sourceB2_n[3]-0.5*geom_sourceB2_n[2]-0.5*geom_sourceB2_n[1]+0.5*geom_sourceB2_n[0]; 

  outJBz[0] += 0.5*geom_sourceB3_n[3]+0.5*geom_sourceB3_n[2]+0.5*geom_sourceB3_n[1]+0.5*geom_sourceB3_n[0]; 
  outJBz[1] += 0.5*geom_sourceB3_n[3]+0.5*geom_sourceB3_n[2]-0.5*geom_sourceB3_n[1]-0.5*geom_sourceB3_n[0]; 
  outJBz[2] += 0.5*geom_sourceB3_n[3]-0.5*geom_sourceB3_n[2]+0.5*geom_sourceB3_n[1]-0.5*geom_sourceB3_n[0]; 
  outJBz[3] += 0.5*geom_sourceB3_n[3]-0.5*geom_sourceB3_n[2]-0.5*geom_sourceB3_n[1]+0.5*geom_sourceB3_n[0]; 

  outJphi[0] += 0.5*geom_sourcephi_n[3]+0.5*geom_sourcephi_n[2]+0.5*geom_sourcephi_n[1]+0.5*geom_sourcephi_n[0]; 
  outJphi[1] += 0.5*geom_sourcephi_n[3]+0.5*geom_sourcephi_n[2]-0.5*geom_sourcephi_n[1]-0.5*geom_sourcephi_n[0]; 
  outJphi[2] += 0.5*geom_sourcephi_n[3]-0.5*geom_sourcephi_n[2]+0.5*geom_sourcephi_n[1]-0.5*geom_sourcephi_n[0]; 
  outJphi[3] += 0.5*geom_sourcephi_n[3]-0.5*geom_sourcephi_n[2]-0.5*geom_sourcephi_n[1]+0.5*geom_sourcephi_n[0]; 

  outJpsi[0] += 0.5*geom_sourcepsi_n[3]+0.5*geom_sourcepsi_n[2]+0.5*geom_sourcepsi_n[1]+0.5*geom_sourcepsi_n[0]; 
  outJpsi[1] += 0.5*geom_sourcepsi_n[3]+0.5*geom_sourcepsi_n[2]-0.5*geom_sourcepsi_n[1]-0.5*geom_sourcepsi_n[0]; 
  outJpsi[2] += 0.5*geom_sourcepsi_n[3]-0.5*geom_sourcepsi_n[2]+0.5*geom_sourcepsi_n[1]-0.5*geom_sourcepsi_n[0]; 
  outJpsi[3] += 0.5*geom_sourcepsi_n[3]-0.5*geom_sourcepsi_n[2]-0.5*geom_sourcepsi_n[1]+0.5*geom_sourcepsi_n[0]; 

} 
