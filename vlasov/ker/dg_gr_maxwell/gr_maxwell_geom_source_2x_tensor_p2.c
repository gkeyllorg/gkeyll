#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH void gr_maxwell_geom_source_2x_tensor_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, 
      const double *geom_factor_con, const double *fields, double* GKYL_RESTRICT out) 
{ 
  const double chi = meq->chi, gamma = meq->gamma; 
  const double K_phi = meq->K_phi, K_psi = meq->K_psi; 
  double dx0 = 2.0/dx[0]; 
  double dx1 = 2.0/dx[1]; 
  
  double Jphi_n[9] = {0.0};
  double Jpsi_n[9] = {0.0};
  double geom_sourceD1_n[9] = {0.0};
  double geom_sourceD2_n[9] = {0.0};
  double geom_sourceD3_n[9] = {0.0};
  double geom_sourceB1_n[9] = {0.0};
  double geom_sourceB2_n[9] = {0.0};
  double geom_sourceB3_n[9] = {0.0};
  double geom_sourcephi_n[9] = {0.0};
  double geom_sourcepsi_n[9] = {0.0};
  
  const double *Jphi = &fields[54]; 
  const double *Jpsi = &fields[63]; 
  
  double *outJDx = &out[0]; 
  double *outJDy = &out[9]; 
  double *outJDz = &out[18]; 
  double *outJBx = &out[27]; 
  double *outJBy = &out[36]; 
  double *outJBz = &out[45]; 
  double *outJphi = &out[54]; 
  double *outJpsi = &out[63]; 
  
  Jphi_n[0] = 0.4*Jphi[8]-0.5999999999999995*Jphi[7]-0.5999999999999999*Jphi[6]+0.4472135954999579*Jphi[5]+0.4472135954999579*Jphi[4]+0.9*Jphi[3]-0.6708203932499369*Jphi[2]-0.6708203932499369*Jphi[1]+0.5*Jphi[0];
  Jpsi_n[0] = 0.4*Jpsi[8]-0.5999999999999995*Jpsi[7]-0.5999999999999999*Jpsi[6]+0.4472135954999579*Jpsi[5]+0.4472135954999579*Jpsi[4]+0.9*Jpsi[3]-0.6708203932499369*Jpsi[2]-0.6708203932499369*Jpsi[1]+0.5*Jpsi[0];
  
  Jphi_n[1] = -(0.5*Jphi[8])+0.75*Jphi[7]-0.5590169943749475*Jphi[5]+0.4472135954999579*Jphi[4]-0.6708203932499369*Jphi[1]+0.5*Jphi[0];
  Jpsi_n[1] = -(0.5*Jpsi[8])+0.75*Jpsi[7]-0.5590169943749475*Jpsi[5]+0.4472135954999579*Jpsi[4]-0.6708203932499369*Jpsi[1]+0.5*Jpsi[0];
  
  Jphi_n[2] = 0.4*Jphi[8]-0.5999999999999995*Jphi[7]+0.5999999999999999*Jphi[6]+0.4472135954999579*Jphi[5]+0.4472135954999579*Jphi[4]-0.9*Jphi[3]+0.6708203932499369*Jphi[2]-0.6708203932499369*Jphi[1]+0.5*Jphi[0];
  Jpsi_n[2] = 0.4*Jpsi[8]-0.5999999999999995*Jpsi[7]+0.5999999999999999*Jpsi[6]+0.4472135954999579*Jpsi[5]+0.4472135954999579*Jpsi[4]-0.9*Jpsi[3]+0.6708203932499369*Jpsi[2]-0.6708203932499369*Jpsi[1]+0.5*Jpsi[0];
  
  Jphi_n[3] = -(0.5*Jphi[8])+0.75*Jphi[6]+0.4472135954999579*Jphi[5]-0.5590169943749475*Jphi[4]-0.6708203932499369*Jphi[2]+0.5*Jphi[0];
  Jpsi_n[3] = -(0.5*Jpsi[8])+0.75*Jpsi[6]+0.4472135954999579*Jpsi[5]-0.5590169943749475*Jpsi[4]-0.6708203932499369*Jpsi[2]+0.5*Jpsi[0];
  
  Jphi_n[4] = 0.625*Jphi[8]-0.5590169943749475*Jphi[5]-0.5590169943749475*Jphi[4]+0.5*Jphi[0];
  Jpsi_n[4] = 0.625*Jpsi[8]-0.5590169943749475*Jpsi[5]-0.5590169943749475*Jpsi[4]+0.5*Jpsi[0];
  
  Jphi_n[5] = -(0.5*Jphi[8])-0.75*Jphi[6]+0.4472135954999579*Jphi[5]-0.5590169943749475*Jphi[4]+0.6708203932499369*Jphi[2]+0.5*Jphi[0];
  Jpsi_n[5] = -(0.5*Jpsi[8])-0.75*Jpsi[6]+0.4472135954999579*Jpsi[5]-0.5590169943749475*Jpsi[4]+0.6708203932499369*Jpsi[2]+0.5*Jpsi[0];
  
  Jphi_n[6] = 0.4*Jphi[8]+0.5999999999999995*Jphi[7]-0.5999999999999999*Jphi[6]+0.4472135954999579*Jphi[5]+0.4472135954999579*Jphi[4]-0.9*Jphi[3]-0.6708203932499369*Jphi[2]+0.6708203932499369*Jphi[1]+0.5*Jphi[0];
  Jpsi_n[6] = 0.4*Jpsi[8]+0.5999999999999995*Jpsi[7]-0.5999999999999999*Jpsi[6]+0.4472135954999579*Jpsi[5]+0.4472135954999579*Jpsi[4]-0.9*Jpsi[3]-0.6708203932499369*Jpsi[2]+0.6708203932499369*Jpsi[1]+0.5*Jpsi[0];
  
  Jphi_n[7] = -(0.5*Jphi[8])-0.75*Jphi[7]-0.5590169943749475*Jphi[5]+0.4472135954999579*Jphi[4]+0.6708203932499369*Jphi[1]+0.5*Jphi[0];
  Jpsi_n[7] = -(0.5*Jpsi[8])-0.75*Jpsi[7]-0.5590169943749475*Jpsi[5]+0.4472135954999579*Jpsi[4]+0.6708203932499369*Jpsi[1]+0.5*Jpsi[0];
  
  Jphi_n[8] = 0.4*Jphi[8]+0.5999999999999995*Jphi[7]+0.5999999999999999*Jphi[6]+0.4472135954999579*Jphi[5]+0.4472135954999579*Jphi[4]+0.9*Jphi[3]+0.6708203932499369*Jphi[2]+0.6708203932499369*Jphi[1]+0.5*Jphi[0];
  Jpsi_n[8] = 0.4*Jpsi[8]+0.5999999999999995*Jpsi[7]+0.5999999999999999*Jpsi[6]+0.4472135954999579*Jpsi[5]+0.4472135954999579*Jpsi[4]+0.9*Jpsi[3]+0.6708203932499369*Jpsi[2]+0.6708203932499369*Jpsi[1]+0.5*Jpsi[0];
  
  const double *geom_factor_con_x = &geom_factor_con[0]; 
  const double *geom_factor_con_y = &geom_factor_con[9]; 
  const double *geom_factor_con_z = &geom_factor_con[18]; 
  
  for (int i=0; i<9; ++i) {
    geom_sourceD1_n[i] = chi * Jphi_n[i] * geom_factor_con_x[i]; 
    geom_sourceD2_n[i] = chi * Jphi_n[i] * geom_factor_con_y[i]; 
    geom_sourceD3_n[i] = chi * Jphi_n[i] * geom_factor_con_z[i]; 
    geom_sourceB1_n[i] = gamma * Jpsi_n[i] * geom_factor_con_x[i]; 
    geom_sourceB2_n[i] = gamma * Jpsi_n[i] * geom_factor_con_y[i]; 
    geom_sourceB3_n[i] = gamma * Jpsi_n[i] * geom_factor_con_z[i]; 
    geom_sourcephi_n[i] = - K_phi * Jphi_n[i]; 
    geom_sourcepsi_n[i] = - K_psi * Jpsi_n[i]; 
  }
  
  outJDx[0] += 0.15432098765432098*geom_sourceD1_n[8]+0.24691358024691357*geom_sourceD1_n[7]+0.15432098765432098*geom_sourceD1_n[6]+0.24691358024691357*geom_sourceD1_n[5]+0.3950617283950617*geom_sourceD1_n[4]+0.24691358024691357*geom_sourceD1_n[3]+0.15432098765432098*geom_sourceD1_n[2]+0.24691358024691357*geom_sourceD1_n[1]+0.15432098765432098*geom_sourceD1_n[0]; 
  outJDx[1] += 0.20704333124998056*geom_sourceD1_n[8]+0.33126932999996883*geom_sourceD1_n[7]+0.20704333124998056*geom_sourceD1_n[6]-0.20704333124998056*geom_sourceD1_n[2]-0.33126932999996883*geom_sourceD1_n[1]-0.20704333124998056*geom_sourceD1_n[0]; 
  outJDx[2] += 0.20704333124998056*geom_sourceD1_n[8]-0.20704333124998056*geom_sourceD1_n[6]+0.33126932999996883*geom_sourceD1_n[5]-0.33126932999996883*geom_sourceD1_n[3]+0.20704333124998056*geom_sourceD1_n[2]-0.20704333124998056*geom_sourceD1_n[0]; 
  outJDx[3] += 0.2777777777777778*geom_sourceD1_n[8]-0.2777777777777778*geom_sourceD1_n[6]-0.2777777777777778*geom_sourceD1_n[2]+0.2777777777777778*geom_sourceD1_n[0]; 
  outJDx[4] += 0.13802888749998704*geom_sourceD1_n[8]+0.22084621999997922*geom_sourceD1_n[7]+0.13802888749998704*geom_sourceD1_n[6]-0.2760577749999741*geom_sourceD1_n[5]-0.44169243999995844*geom_sourceD1_n[4]-0.2760577749999741*geom_sourceD1_n[3]+0.13802888749998704*geom_sourceD1_n[2]+0.22084621999997922*geom_sourceD1_n[1]+0.13802888749998704*geom_sourceD1_n[0]; 
  outJDx[5] += 0.13802888749998704*geom_sourceD1_n[8]-0.2760577749999741*geom_sourceD1_n[7]+0.13802888749998704*geom_sourceD1_n[6]+0.22084621999997922*geom_sourceD1_n[5]-0.44169243999995844*geom_sourceD1_n[4]+0.22084621999997922*geom_sourceD1_n[3]+0.13802888749998704*geom_sourceD1_n[2]-0.2760577749999741*geom_sourceD1_n[1]+0.13802888749998704*geom_sourceD1_n[0]; 
  outJDx[6] += 0.18518518518518526*geom_sourceD1_n[8]-0.18518518518518526*geom_sourceD1_n[6]-0.3703703703703705*geom_sourceD1_n[5]+0.3703703703703705*geom_sourceD1_n[3]+0.18518518518518526*geom_sourceD1_n[2]-0.18518518518518526*geom_sourceD1_n[0]; 
  outJDx[7] += 0.18518518518518526*geom_sourceD1_n[8]-0.3703703703703705*geom_sourceD1_n[7]+0.18518518518518526*geom_sourceD1_n[6]-0.18518518518518526*geom_sourceD1_n[2]+0.3703703703703705*geom_sourceD1_n[1]-0.18518518518518526*geom_sourceD1_n[0]; 
  outJDx[8] += 0.12345679012345678*geom_sourceD1_n[8]-0.24691358024691357*geom_sourceD1_n[7]+0.12345679012345678*geom_sourceD1_n[6]-0.24691358024691357*geom_sourceD1_n[5]+0.49382716049382713*geom_sourceD1_n[4]-0.24691358024691357*geom_sourceD1_n[3]+0.12345679012345678*geom_sourceD1_n[2]-0.24691358024691357*geom_sourceD1_n[1]+0.12345679012345678*geom_sourceD1_n[0]; 

  outJDy[0] += 0.15432098765432098*geom_sourceD2_n[8]+0.24691358024691357*geom_sourceD2_n[7]+0.15432098765432098*geom_sourceD2_n[6]+0.24691358024691357*geom_sourceD2_n[5]+0.3950617283950617*geom_sourceD2_n[4]+0.24691358024691357*geom_sourceD2_n[3]+0.15432098765432098*geom_sourceD2_n[2]+0.24691358024691357*geom_sourceD2_n[1]+0.15432098765432098*geom_sourceD2_n[0]; 
  outJDy[1] += 0.20704333124998056*geom_sourceD2_n[8]+0.33126932999996883*geom_sourceD2_n[7]+0.20704333124998056*geom_sourceD2_n[6]-0.20704333124998056*geom_sourceD2_n[2]-0.33126932999996883*geom_sourceD2_n[1]-0.20704333124998056*geom_sourceD2_n[0]; 
  outJDy[2] += 0.20704333124998056*geom_sourceD2_n[8]-0.20704333124998056*geom_sourceD2_n[6]+0.33126932999996883*geom_sourceD2_n[5]-0.33126932999996883*geom_sourceD2_n[3]+0.20704333124998056*geom_sourceD2_n[2]-0.20704333124998056*geom_sourceD2_n[0]; 
  outJDy[3] += 0.2777777777777778*geom_sourceD2_n[8]-0.2777777777777778*geom_sourceD2_n[6]-0.2777777777777778*geom_sourceD2_n[2]+0.2777777777777778*geom_sourceD2_n[0]; 
  outJDy[4] += 0.13802888749998704*geom_sourceD2_n[8]+0.22084621999997922*geom_sourceD2_n[7]+0.13802888749998704*geom_sourceD2_n[6]-0.2760577749999741*geom_sourceD2_n[5]-0.44169243999995844*geom_sourceD2_n[4]-0.2760577749999741*geom_sourceD2_n[3]+0.13802888749998704*geom_sourceD2_n[2]+0.22084621999997922*geom_sourceD2_n[1]+0.13802888749998704*geom_sourceD2_n[0]; 
  outJDy[5] += 0.13802888749998704*geom_sourceD2_n[8]-0.2760577749999741*geom_sourceD2_n[7]+0.13802888749998704*geom_sourceD2_n[6]+0.22084621999997922*geom_sourceD2_n[5]-0.44169243999995844*geom_sourceD2_n[4]+0.22084621999997922*geom_sourceD2_n[3]+0.13802888749998704*geom_sourceD2_n[2]-0.2760577749999741*geom_sourceD2_n[1]+0.13802888749998704*geom_sourceD2_n[0]; 
  outJDy[6] += 0.18518518518518526*geom_sourceD2_n[8]-0.18518518518518526*geom_sourceD2_n[6]-0.3703703703703705*geom_sourceD2_n[5]+0.3703703703703705*geom_sourceD2_n[3]+0.18518518518518526*geom_sourceD2_n[2]-0.18518518518518526*geom_sourceD2_n[0]; 
  outJDy[7] += 0.18518518518518526*geom_sourceD2_n[8]-0.3703703703703705*geom_sourceD2_n[7]+0.18518518518518526*geom_sourceD2_n[6]-0.18518518518518526*geom_sourceD2_n[2]+0.3703703703703705*geom_sourceD2_n[1]-0.18518518518518526*geom_sourceD2_n[0]; 
  outJDy[8] += 0.12345679012345678*geom_sourceD2_n[8]-0.24691358024691357*geom_sourceD2_n[7]+0.12345679012345678*geom_sourceD2_n[6]-0.24691358024691357*geom_sourceD2_n[5]+0.49382716049382713*geom_sourceD2_n[4]-0.24691358024691357*geom_sourceD2_n[3]+0.12345679012345678*geom_sourceD2_n[2]-0.24691358024691357*geom_sourceD2_n[1]+0.12345679012345678*geom_sourceD2_n[0]; 

  outJDz[0] += 0.15432098765432098*geom_sourceD3_n[8]+0.24691358024691357*geom_sourceD3_n[7]+0.15432098765432098*geom_sourceD3_n[6]+0.24691358024691357*geom_sourceD3_n[5]+0.3950617283950617*geom_sourceD3_n[4]+0.24691358024691357*geom_sourceD3_n[3]+0.15432098765432098*geom_sourceD3_n[2]+0.24691358024691357*geom_sourceD3_n[1]+0.15432098765432098*geom_sourceD3_n[0]; 
  outJDz[1] += 0.20704333124998056*geom_sourceD3_n[8]+0.33126932999996883*geom_sourceD3_n[7]+0.20704333124998056*geom_sourceD3_n[6]-0.20704333124998056*geom_sourceD3_n[2]-0.33126932999996883*geom_sourceD3_n[1]-0.20704333124998056*geom_sourceD3_n[0]; 
  outJDz[2] += 0.20704333124998056*geom_sourceD3_n[8]-0.20704333124998056*geom_sourceD3_n[6]+0.33126932999996883*geom_sourceD3_n[5]-0.33126932999996883*geom_sourceD3_n[3]+0.20704333124998056*geom_sourceD3_n[2]-0.20704333124998056*geom_sourceD3_n[0]; 
  outJDz[3] += 0.2777777777777778*geom_sourceD3_n[8]-0.2777777777777778*geom_sourceD3_n[6]-0.2777777777777778*geom_sourceD3_n[2]+0.2777777777777778*geom_sourceD3_n[0]; 
  outJDz[4] += 0.13802888749998704*geom_sourceD3_n[8]+0.22084621999997922*geom_sourceD3_n[7]+0.13802888749998704*geom_sourceD3_n[6]-0.2760577749999741*geom_sourceD3_n[5]-0.44169243999995844*geom_sourceD3_n[4]-0.2760577749999741*geom_sourceD3_n[3]+0.13802888749998704*geom_sourceD3_n[2]+0.22084621999997922*geom_sourceD3_n[1]+0.13802888749998704*geom_sourceD3_n[0]; 
  outJDz[5] += 0.13802888749998704*geom_sourceD3_n[8]-0.2760577749999741*geom_sourceD3_n[7]+0.13802888749998704*geom_sourceD3_n[6]+0.22084621999997922*geom_sourceD3_n[5]-0.44169243999995844*geom_sourceD3_n[4]+0.22084621999997922*geom_sourceD3_n[3]+0.13802888749998704*geom_sourceD3_n[2]-0.2760577749999741*geom_sourceD3_n[1]+0.13802888749998704*geom_sourceD3_n[0]; 
  outJDz[6] += 0.18518518518518526*geom_sourceD3_n[8]-0.18518518518518526*geom_sourceD3_n[6]-0.3703703703703705*geom_sourceD3_n[5]+0.3703703703703705*geom_sourceD3_n[3]+0.18518518518518526*geom_sourceD3_n[2]-0.18518518518518526*geom_sourceD3_n[0]; 
  outJDz[7] += 0.18518518518518526*geom_sourceD3_n[8]-0.3703703703703705*geom_sourceD3_n[7]+0.18518518518518526*geom_sourceD3_n[6]-0.18518518518518526*geom_sourceD3_n[2]+0.3703703703703705*geom_sourceD3_n[1]-0.18518518518518526*geom_sourceD3_n[0]; 
  outJDz[8] += 0.12345679012345678*geom_sourceD3_n[8]-0.24691358024691357*geom_sourceD3_n[7]+0.12345679012345678*geom_sourceD3_n[6]-0.24691358024691357*geom_sourceD3_n[5]+0.49382716049382713*geom_sourceD3_n[4]-0.24691358024691357*geom_sourceD3_n[3]+0.12345679012345678*geom_sourceD3_n[2]-0.24691358024691357*geom_sourceD3_n[1]+0.12345679012345678*geom_sourceD3_n[0]; 

  outJBx[0] += 0.15432098765432098*geom_sourceB1_n[8]+0.24691358024691357*geom_sourceB1_n[7]+0.15432098765432098*geom_sourceB1_n[6]+0.24691358024691357*geom_sourceB1_n[5]+0.3950617283950617*geom_sourceB1_n[4]+0.24691358024691357*geom_sourceB1_n[3]+0.15432098765432098*geom_sourceB1_n[2]+0.24691358024691357*geom_sourceB1_n[1]+0.15432098765432098*geom_sourceB1_n[0]; 
  outJBx[1] += 0.20704333124998056*geom_sourceB1_n[8]+0.33126932999996883*geom_sourceB1_n[7]+0.20704333124998056*geom_sourceB1_n[6]-0.20704333124998056*geom_sourceB1_n[2]-0.33126932999996883*geom_sourceB1_n[1]-0.20704333124998056*geom_sourceB1_n[0]; 
  outJBx[2] += 0.20704333124998056*geom_sourceB1_n[8]-0.20704333124998056*geom_sourceB1_n[6]+0.33126932999996883*geom_sourceB1_n[5]-0.33126932999996883*geom_sourceB1_n[3]+0.20704333124998056*geom_sourceB1_n[2]-0.20704333124998056*geom_sourceB1_n[0]; 
  outJBx[3] += 0.2777777777777778*geom_sourceB1_n[8]-0.2777777777777778*geom_sourceB1_n[6]-0.2777777777777778*geom_sourceB1_n[2]+0.2777777777777778*geom_sourceB1_n[0]; 
  outJBx[4] += 0.13802888749998704*geom_sourceB1_n[8]+0.22084621999997922*geom_sourceB1_n[7]+0.13802888749998704*geom_sourceB1_n[6]-0.2760577749999741*geom_sourceB1_n[5]-0.44169243999995844*geom_sourceB1_n[4]-0.2760577749999741*geom_sourceB1_n[3]+0.13802888749998704*geom_sourceB1_n[2]+0.22084621999997922*geom_sourceB1_n[1]+0.13802888749998704*geom_sourceB1_n[0]; 
  outJBx[5] += 0.13802888749998704*geom_sourceB1_n[8]-0.2760577749999741*geom_sourceB1_n[7]+0.13802888749998704*geom_sourceB1_n[6]+0.22084621999997922*geom_sourceB1_n[5]-0.44169243999995844*geom_sourceB1_n[4]+0.22084621999997922*geom_sourceB1_n[3]+0.13802888749998704*geom_sourceB1_n[2]-0.2760577749999741*geom_sourceB1_n[1]+0.13802888749998704*geom_sourceB1_n[0]; 
  outJBx[6] += 0.18518518518518526*geom_sourceB1_n[8]-0.18518518518518526*geom_sourceB1_n[6]-0.3703703703703705*geom_sourceB1_n[5]+0.3703703703703705*geom_sourceB1_n[3]+0.18518518518518526*geom_sourceB1_n[2]-0.18518518518518526*geom_sourceB1_n[0]; 
  outJBx[7] += 0.18518518518518526*geom_sourceB1_n[8]-0.3703703703703705*geom_sourceB1_n[7]+0.18518518518518526*geom_sourceB1_n[6]-0.18518518518518526*geom_sourceB1_n[2]+0.3703703703703705*geom_sourceB1_n[1]-0.18518518518518526*geom_sourceB1_n[0]; 
  outJBx[8] += 0.12345679012345678*geom_sourceB1_n[8]-0.24691358024691357*geom_sourceB1_n[7]+0.12345679012345678*geom_sourceB1_n[6]-0.24691358024691357*geom_sourceB1_n[5]+0.49382716049382713*geom_sourceB1_n[4]-0.24691358024691357*geom_sourceB1_n[3]+0.12345679012345678*geom_sourceB1_n[2]-0.24691358024691357*geom_sourceB1_n[1]+0.12345679012345678*geom_sourceB1_n[0]; 

  outJBy[0] += 0.15432098765432098*geom_sourceB2_n[8]+0.24691358024691357*geom_sourceB2_n[7]+0.15432098765432098*geom_sourceB2_n[6]+0.24691358024691357*geom_sourceB2_n[5]+0.3950617283950617*geom_sourceB2_n[4]+0.24691358024691357*geom_sourceB2_n[3]+0.15432098765432098*geom_sourceB2_n[2]+0.24691358024691357*geom_sourceB2_n[1]+0.15432098765432098*geom_sourceB2_n[0]; 
  outJBy[1] += 0.20704333124998056*geom_sourceB2_n[8]+0.33126932999996883*geom_sourceB2_n[7]+0.20704333124998056*geom_sourceB2_n[6]-0.20704333124998056*geom_sourceB2_n[2]-0.33126932999996883*geom_sourceB2_n[1]-0.20704333124998056*geom_sourceB2_n[0]; 
  outJBy[2] += 0.20704333124998056*geom_sourceB2_n[8]-0.20704333124998056*geom_sourceB2_n[6]+0.33126932999996883*geom_sourceB2_n[5]-0.33126932999996883*geom_sourceB2_n[3]+0.20704333124998056*geom_sourceB2_n[2]-0.20704333124998056*geom_sourceB2_n[0]; 
  outJBy[3] += 0.2777777777777778*geom_sourceB2_n[8]-0.2777777777777778*geom_sourceB2_n[6]-0.2777777777777778*geom_sourceB2_n[2]+0.2777777777777778*geom_sourceB2_n[0]; 
  outJBy[4] += 0.13802888749998704*geom_sourceB2_n[8]+0.22084621999997922*geom_sourceB2_n[7]+0.13802888749998704*geom_sourceB2_n[6]-0.2760577749999741*geom_sourceB2_n[5]-0.44169243999995844*geom_sourceB2_n[4]-0.2760577749999741*geom_sourceB2_n[3]+0.13802888749998704*geom_sourceB2_n[2]+0.22084621999997922*geom_sourceB2_n[1]+0.13802888749998704*geom_sourceB2_n[0]; 
  outJBy[5] += 0.13802888749998704*geom_sourceB2_n[8]-0.2760577749999741*geom_sourceB2_n[7]+0.13802888749998704*geom_sourceB2_n[6]+0.22084621999997922*geom_sourceB2_n[5]-0.44169243999995844*geom_sourceB2_n[4]+0.22084621999997922*geom_sourceB2_n[3]+0.13802888749998704*geom_sourceB2_n[2]-0.2760577749999741*geom_sourceB2_n[1]+0.13802888749998704*geom_sourceB2_n[0]; 
  outJBy[6] += 0.18518518518518526*geom_sourceB2_n[8]-0.18518518518518526*geom_sourceB2_n[6]-0.3703703703703705*geom_sourceB2_n[5]+0.3703703703703705*geom_sourceB2_n[3]+0.18518518518518526*geom_sourceB2_n[2]-0.18518518518518526*geom_sourceB2_n[0]; 
  outJBy[7] += 0.18518518518518526*geom_sourceB2_n[8]-0.3703703703703705*geom_sourceB2_n[7]+0.18518518518518526*geom_sourceB2_n[6]-0.18518518518518526*geom_sourceB2_n[2]+0.3703703703703705*geom_sourceB2_n[1]-0.18518518518518526*geom_sourceB2_n[0]; 
  outJBy[8] += 0.12345679012345678*geom_sourceB2_n[8]-0.24691358024691357*geom_sourceB2_n[7]+0.12345679012345678*geom_sourceB2_n[6]-0.24691358024691357*geom_sourceB2_n[5]+0.49382716049382713*geom_sourceB2_n[4]-0.24691358024691357*geom_sourceB2_n[3]+0.12345679012345678*geom_sourceB2_n[2]-0.24691358024691357*geom_sourceB2_n[1]+0.12345679012345678*geom_sourceB2_n[0]; 

  outJBz[0] += 0.15432098765432098*geom_sourceB3_n[8]+0.24691358024691357*geom_sourceB3_n[7]+0.15432098765432098*geom_sourceB3_n[6]+0.24691358024691357*geom_sourceB3_n[5]+0.3950617283950617*geom_sourceB3_n[4]+0.24691358024691357*geom_sourceB3_n[3]+0.15432098765432098*geom_sourceB3_n[2]+0.24691358024691357*geom_sourceB3_n[1]+0.15432098765432098*geom_sourceB3_n[0]; 
  outJBz[1] += 0.20704333124998056*geom_sourceB3_n[8]+0.33126932999996883*geom_sourceB3_n[7]+0.20704333124998056*geom_sourceB3_n[6]-0.20704333124998056*geom_sourceB3_n[2]-0.33126932999996883*geom_sourceB3_n[1]-0.20704333124998056*geom_sourceB3_n[0]; 
  outJBz[2] += 0.20704333124998056*geom_sourceB3_n[8]-0.20704333124998056*geom_sourceB3_n[6]+0.33126932999996883*geom_sourceB3_n[5]-0.33126932999996883*geom_sourceB3_n[3]+0.20704333124998056*geom_sourceB3_n[2]-0.20704333124998056*geom_sourceB3_n[0]; 
  outJBz[3] += 0.2777777777777778*geom_sourceB3_n[8]-0.2777777777777778*geom_sourceB3_n[6]-0.2777777777777778*geom_sourceB3_n[2]+0.2777777777777778*geom_sourceB3_n[0]; 
  outJBz[4] += 0.13802888749998704*geom_sourceB3_n[8]+0.22084621999997922*geom_sourceB3_n[7]+0.13802888749998704*geom_sourceB3_n[6]-0.2760577749999741*geom_sourceB3_n[5]-0.44169243999995844*geom_sourceB3_n[4]-0.2760577749999741*geom_sourceB3_n[3]+0.13802888749998704*geom_sourceB3_n[2]+0.22084621999997922*geom_sourceB3_n[1]+0.13802888749998704*geom_sourceB3_n[0]; 
  outJBz[5] += 0.13802888749998704*geom_sourceB3_n[8]-0.2760577749999741*geom_sourceB3_n[7]+0.13802888749998704*geom_sourceB3_n[6]+0.22084621999997922*geom_sourceB3_n[5]-0.44169243999995844*geom_sourceB3_n[4]+0.22084621999997922*geom_sourceB3_n[3]+0.13802888749998704*geom_sourceB3_n[2]-0.2760577749999741*geom_sourceB3_n[1]+0.13802888749998704*geom_sourceB3_n[0]; 
  outJBz[6] += 0.18518518518518526*geom_sourceB3_n[8]-0.18518518518518526*geom_sourceB3_n[6]-0.3703703703703705*geom_sourceB3_n[5]+0.3703703703703705*geom_sourceB3_n[3]+0.18518518518518526*geom_sourceB3_n[2]-0.18518518518518526*geom_sourceB3_n[0]; 
  outJBz[7] += 0.18518518518518526*geom_sourceB3_n[8]-0.3703703703703705*geom_sourceB3_n[7]+0.18518518518518526*geom_sourceB3_n[6]-0.18518518518518526*geom_sourceB3_n[2]+0.3703703703703705*geom_sourceB3_n[1]-0.18518518518518526*geom_sourceB3_n[0]; 
  outJBz[8] += 0.12345679012345678*geom_sourceB3_n[8]-0.24691358024691357*geom_sourceB3_n[7]+0.12345679012345678*geom_sourceB3_n[6]-0.24691358024691357*geom_sourceB3_n[5]+0.49382716049382713*geom_sourceB3_n[4]-0.24691358024691357*geom_sourceB3_n[3]+0.12345679012345678*geom_sourceB3_n[2]-0.24691358024691357*geom_sourceB3_n[1]+0.12345679012345678*geom_sourceB3_n[0]; 

  outJphi[0] += 0.15432098765432098*geom_sourcephi_n[8]+0.24691358024691357*geom_sourcephi_n[7]+0.15432098765432098*geom_sourcephi_n[6]+0.24691358024691357*geom_sourcephi_n[5]+0.3950617283950617*geom_sourcephi_n[4]+0.24691358024691357*geom_sourcephi_n[3]+0.15432098765432098*geom_sourcephi_n[2]+0.24691358024691357*geom_sourcephi_n[1]+0.15432098765432098*geom_sourcephi_n[0]; 
  outJphi[1] += 0.20704333124998056*geom_sourcephi_n[8]+0.33126932999996883*geom_sourcephi_n[7]+0.20704333124998056*geom_sourcephi_n[6]-0.20704333124998056*geom_sourcephi_n[2]-0.33126932999996883*geom_sourcephi_n[1]-0.20704333124998056*geom_sourcephi_n[0]; 
  outJphi[2] += 0.20704333124998056*geom_sourcephi_n[8]-0.20704333124998056*geom_sourcephi_n[6]+0.33126932999996883*geom_sourcephi_n[5]-0.33126932999996883*geom_sourcephi_n[3]+0.20704333124998056*geom_sourcephi_n[2]-0.20704333124998056*geom_sourcephi_n[0]; 
  outJphi[3] += 0.2777777777777778*geom_sourcephi_n[8]-0.2777777777777778*geom_sourcephi_n[6]-0.2777777777777778*geom_sourcephi_n[2]+0.2777777777777778*geom_sourcephi_n[0]; 
  outJphi[4] += 0.13802888749998704*geom_sourcephi_n[8]+0.22084621999997922*geom_sourcephi_n[7]+0.13802888749998704*geom_sourcephi_n[6]-0.2760577749999741*geom_sourcephi_n[5]-0.44169243999995844*geom_sourcephi_n[4]-0.2760577749999741*geom_sourcephi_n[3]+0.13802888749998704*geom_sourcephi_n[2]+0.22084621999997922*geom_sourcephi_n[1]+0.13802888749998704*geom_sourcephi_n[0]; 
  outJphi[5] += 0.13802888749998704*geom_sourcephi_n[8]-0.2760577749999741*geom_sourcephi_n[7]+0.13802888749998704*geom_sourcephi_n[6]+0.22084621999997922*geom_sourcephi_n[5]-0.44169243999995844*geom_sourcephi_n[4]+0.22084621999997922*geom_sourcephi_n[3]+0.13802888749998704*geom_sourcephi_n[2]-0.2760577749999741*geom_sourcephi_n[1]+0.13802888749998704*geom_sourcephi_n[0]; 
  outJphi[6] += 0.18518518518518526*geom_sourcephi_n[8]-0.18518518518518526*geom_sourcephi_n[6]-0.3703703703703705*geom_sourcephi_n[5]+0.3703703703703705*geom_sourcephi_n[3]+0.18518518518518526*geom_sourcephi_n[2]-0.18518518518518526*geom_sourcephi_n[0]; 
  outJphi[7] += 0.18518518518518526*geom_sourcephi_n[8]-0.3703703703703705*geom_sourcephi_n[7]+0.18518518518518526*geom_sourcephi_n[6]-0.18518518518518526*geom_sourcephi_n[2]+0.3703703703703705*geom_sourcephi_n[1]-0.18518518518518526*geom_sourcephi_n[0]; 
  outJphi[8] += 0.12345679012345678*geom_sourcephi_n[8]-0.24691358024691357*geom_sourcephi_n[7]+0.12345679012345678*geom_sourcephi_n[6]-0.24691358024691357*geom_sourcephi_n[5]+0.49382716049382713*geom_sourcephi_n[4]-0.24691358024691357*geom_sourcephi_n[3]+0.12345679012345678*geom_sourcephi_n[2]-0.24691358024691357*geom_sourcephi_n[1]+0.12345679012345678*geom_sourcephi_n[0]; 

  outJpsi[0] += 0.15432098765432098*geom_sourcepsi_n[8]+0.24691358024691357*geom_sourcepsi_n[7]+0.15432098765432098*geom_sourcepsi_n[6]+0.24691358024691357*geom_sourcepsi_n[5]+0.3950617283950617*geom_sourcepsi_n[4]+0.24691358024691357*geom_sourcepsi_n[3]+0.15432098765432098*geom_sourcepsi_n[2]+0.24691358024691357*geom_sourcepsi_n[1]+0.15432098765432098*geom_sourcepsi_n[0]; 
  outJpsi[1] += 0.20704333124998056*geom_sourcepsi_n[8]+0.33126932999996883*geom_sourcepsi_n[7]+0.20704333124998056*geom_sourcepsi_n[6]-0.20704333124998056*geom_sourcepsi_n[2]-0.33126932999996883*geom_sourcepsi_n[1]-0.20704333124998056*geom_sourcepsi_n[0]; 
  outJpsi[2] += 0.20704333124998056*geom_sourcepsi_n[8]-0.20704333124998056*geom_sourcepsi_n[6]+0.33126932999996883*geom_sourcepsi_n[5]-0.33126932999996883*geom_sourcepsi_n[3]+0.20704333124998056*geom_sourcepsi_n[2]-0.20704333124998056*geom_sourcepsi_n[0]; 
  outJpsi[3] += 0.2777777777777778*geom_sourcepsi_n[8]-0.2777777777777778*geom_sourcepsi_n[6]-0.2777777777777778*geom_sourcepsi_n[2]+0.2777777777777778*geom_sourcepsi_n[0]; 
  outJpsi[4] += 0.13802888749998704*geom_sourcepsi_n[8]+0.22084621999997922*geom_sourcepsi_n[7]+0.13802888749998704*geom_sourcepsi_n[6]-0.2760577749999741*geom_sourcepsi_n[5]-0.44169243999995844*geom_sourcepsi_n[4]-0.2760577749999741*geom_sourcepsi_n[3]+0.13802888749998704*geom_sourcepsi_n[2]+0.22084621999997922*geom_sourcepsi_n[1]+0.13802888749998704*geom_sourcepsi_n[0]; 
  outJpsi[5] += 0.13802888749998704*geom_sourcepsi_n[8]-0.2760577749999741*geom_sourcepsi_n[7]+0.13802888749998704*geom_sourcepsi_n[6]+0.22084621999997922*geom_sourcepsi_n[5]-0.44169243999995844*geom_sourcepsi_n[4]+0.22084621999997922*geom_sourcepsi_n[3]+0.13802888749998704*geom_sourcepsi_n[2]-0.2760577749999741*geom_sourcepsi_n[1]+0.13802888749998704*geom_sourcepsi_n[0]; 
  outJpsi[6] += 0.18518518518518526*geom_sourcepsi_n[8]-0.18518518518518526*geom_sourcepsi_n[6]-0.3703703703703705*geom_sourcepsi_n[5]+0.3703703703703705*geom_sourcepsi_n[3]+0.18518518518518526*geom_sourcepsi_n[2]-0.18518518518518526*geom_sourcepsi_n[0]; 
  outJpsi[7] += 0.18518518518518526*geom_sourcepsi_n[8]-0.3703703703703705*geom_sourcepsi_n[7]+0.18518518518518526*geom_sourcepsi_n[6]-0.18518518518518526*geom_sourcepsi_n[2]+0.3703703703703705*geom_sourcepsi_n[1]-0.18518518518518526*geom_sourcepsi_n[0]; 
  outJpsi[8] += 0.12345679012345678*geom_sourcepsi_n[8]-0.24691358024691357*geom_sourcepsi_n[7]+0.12345679012345678*geom_sourcepsi_n[6]-0.24691358024691357*geom_sourcepsi_n[5]+0.49382716049382713*geom_sourcepsi_n[4]-0.24691358024691357*geom_sourcepsi_n[3]+0.12345679012345678*geom_sourcepsi_n[2]-0.24691358024691357*geom_sourcepsi_n[1]+0.12345679012345678*geom_sourcepsi_n[0]; 

} 
