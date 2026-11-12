#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfx_2x1v_ser_p3(const double *w, const double *dxv,
  const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel,
  const double *poisson_tensor_conf, const double *hamil,
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[2]; 
  double wv = w[2]; 

  double fUpwind[12]; 
  double Ghat[12]; 

  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 

  const double *jacob_cx_edge = &jacob_pos_edge[0]; 
  const double *jacob_cx_skin = &jacob_pos_skin[0]; 
  const double jacob_cx_edge_inv = 1.0/jacob_cx_edge[0]; 
  const double jacob_cx_skin_inv = 1.0/jacob_cx_skin[0]; 

  if (edge == -1) { 

  if (wv>0) { 

  fUpwind[0] = (1.8708286933869707*fskin[17]+1.5811388300841895*fskin[7]+1.224744871391589*fskin[1]+0.7071067811865475*fskin[0])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[1] = (1.8708286933869707*fskin[23]+1.5811388300841898*fskin[11]+1.224744871391589*fskin[4]+0.7071067811865475*fskin[2])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[2] = (1.8708286933869707*fskin[25]+1.5811388300841898*fskin[13]+1.224744871391589*fskin[5]+0.7071067811865475*fskin[3])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[3] = (1.8708286933869707*fskin[29]+1.5811388300841895*fskin[20]+1.224744871391589*fskin[10]+0.7071067811865475*fskin[6])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[4] = (1.224744871391589*fskin[12]+0.7071067811865475*fskin[8])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[5] = (1.224744871391589*fskin[15]+0.7071067811865475*fskin[9])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[6] = (1.224744871391589*fskin[21]+0.7071067811865475*fskin[14])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[7] = (1.224744871391589*fskin[22]+0.7071067811865475*fskin[16])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[8] = (1.2247448713915887*fskin[24]+0.7071067811865475*fskin[18])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[9] = (1.2247448713915887*fskin[27]+0.7071067811865475*fskin[19])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[10] = (1.2247448713915887*fskin[30]+0.7071067811865475*fskin[26])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[11] = (1.2247448713915887*fskin[31]+0.7071067811865475*fskin[28])*jacob_cx_skin_inv*jacob_vx_inv; 

  } else { 

  fUpwind[0] = (-(1.8708286933869707*fedge[17])+1.5811388300841895*fedge[7]-1.224744871391589*fedge[1]+0.7071067811865475*fedge[0])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[1] = (-(1.8708286933869707*fedge[23])+1.5811388300841898*fedge[11]-1.224744871391589*fedge[4]+0.7071067811865475*fedge[2])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[2] = (-(1.8708286933869707*fedge[25])+1.5811388300841898*fedge[13]-1.224744871391589*fedge[5]+0.7071067811865475*fedge[3])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[3] = (-(1.8708286933869707*fedge[29])+1.5811388300841895*fedge[20]-1.224744871391589*fedge[10]+0.7071067811865475*fedge[6])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[4] = (0.7071067811865475*fedge[8]-1.224744871391589*fedge[12])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[5] = (0.7071067811865475*fedge[9]-1.224744871391589*fedge[15])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[6] = (0.7071067811865475*fedge[14]-1.224744871391589*fedge[21])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[7] = (0.7071067811865475*fedge[16]-1.224744871391589*fedge[22])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[8] = (0.7071067811865475*fedge[18]-1.2247448713915887*fedge[24])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[9] = (0.7071067811865475*fedge[19]-1.2247448713915887*fedge[27])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[10] = (0.7071067811865475*fedge[26]-1.2247448713915887*fedge[30])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[11] = (0.7071067811865475*fedge[28]-1.2247448713915887*fedge[31])*jacob_cx_edge_inv*jacob_vx_inv; 

  } 
  Ghat[0] = 4.183300132670378*hamil[3]*fUpwind[5]+1.8708286933869707*fUpwind[0]*hamil[3]+2.7386127875258306*fUpwind[2]*hamil[2]+1.224744871391589*fUpwind[0]*hamil[1]; 
  Ghat[1] = 4.183300132670379*hamil[3]*fUpwind[7]+1.8708286933869707*fUpwind[1]*hamil[3]+2.7386127875258306*hamil[2]*fUpwind[3]+1.224744871391589*fUpwind[1]*hamil[1]; 
  Ghat[2] = 3.6742346141747664*hamil[3]*fUpwind[9]+2.4494897427831783*hamil[2]*fUpwind[5]+5.612486080160912*fUpwind[2]*hamil[3]+2.7386127875258306*fUpwind[0]*hamil[2]+1.224744871391589*hamil[1]*fUpwind[2]; 
  Ghat[3] = 3.674234614174766*hamil[3]*fUpwind[11]+2.4494897427831783*hamil[2]*fUpwind[7]+5.612486080160912*fUpwind[3]*hamil[3]+1.224744871391589*hamil[1]*fUpwind[3]+2.7386127875258306*fUpwind[1]*hamil[2]; 
  Ghat[4] = 2.7386127875258306*hamil[2]*fUpwind[6]+1.8708286933869707*hamil[3]*fUpwind[4]+1.224744871391589*hamil[1]*fUpwind[4]; 
  Ghat[5] = 2.405351177211819*hamil[2]*fUpwind[9]+4.543441112511214*hamil[3]*fUpwind[5]+1.224744871391589*hamil[1]*fUpwind[5]+4.183300132670378*fUpwind[0]*hamil[3]+2.4494897427831783*fUpwind[2]*hamil[2]; 
  Ghat[6] = 5.612486080160912*hamil[3]*fUpwind[6]+1.224744871391589*hamil[1]*fUpwind[6]+2.7386127875258306*hamil[2]*fUpwind[4]; 
  Ghat[7] = 2.405351177211819*hamil[2]*fUpwind[11]+4.543441112511214*hamil[3]*fUpwind[7]+1.224744871391589*hamil[1]*fUpwind[7]+4.183300132670379*fUpwind[1]*hamil[3]+2.4494897427831783*hamil[2]*fUpwind[3]; 
  Ghat[8] = 2.7386127875258297*hamil[2]*fUpwind[10]+1.8708286933869707*hamil[3]*fUpwind[8]+1.224744871391589*hamil[1]*fUpwind[8]; 
  Ghat[9] = 4.365266951236265*hamil[3]*fUpwind[9]+1.224744871391589*hamil[1]*fUpwind[9]+2.405351177211819*hamil[2]*fUpwind[5]+3.6742346141747664*fUpwind[2]*hamil[3]; 
  Ghat[10] = 5.612486080160912*hamil[3]*fUpwind[10]+1.224744871391589*hamil[1]*fUpwind[10]+2.7386127875258297*hamil[2]*fUpwind[8]; 
  Ghat[11] = 4.365266951236265*hamil[3]*fUpwind[11]+1.224744871391589*hamil[1]*fUpwind[11]+2.405351177211819*hamil[2]*fUpwind[7]+3.674234614174766*fUpwind[3]*hamil[3]; 

  out[0] += -(0.7071067811865475*Ghat[0]*dv10*dx10); 
  out[1] += -(1.224744871391589*Ghat[0]*dv10*dx10); 
  out[2] += -(0.7071067811865475*Ghat[1]*dv10*dx10); 
  out[3] += -(0.7071067811865475*Ghat[2]*dv10*dx10); 
  out[4] += -(1.224744871391589*Ghat[1]*dv10*dx10); 
  out[5] += -(1.224744871391589*Ghat[2]*dv10*dx10); 
  out[6] += -(0.7071067811865475*Ghat[3]*dv10*dx10); 
  out[7] += -(1.5811388300841895*Ghat[0]*dv10*dx10); 
  out[8] += -(0.7071067811865475*Ghat[4]*dv10*dx10); 
  out[9] += -(0.7071067811865475*Ghat[5]*dv10*dx10); 
  out[10] += -(1.224744871391589*Ghat[3]*dv10*dx10); 
  out[11] += -(1.5811388300841898*Ghat[1]*dv10*dx10); 
  out[12] += -(1.224744871391589*Ghat[4]*dv10*dx10); 
  out[13] += -(1.5811388300841898*Ghat[2]*dv10*dx10); 
  out[14] += -(0.7071067811865475*Ghat[6]*dv10*dx10); 
  out[15] += -(1.224744871391589*Ghat[5]*dv10*dx10); 
  out[16] += -(0.7071067811865475*Ghat[7]*dv10*dx10); 
  out[17] += -(1.8708286933869707*Ghat[0]*dv10*dx10); 
  out[18] += -(0.7071067811865475*Ghat[8]*dv10*dx10); 
  out[19] += -(0.7071067811865475*Ghat[9]*dv10*dx10); 
  out[20] += -(1.5811388300841895*Ghat[3]*dv10*dx10); 
  out[21] += -(1.224744871391589*Ghat[6]*dv10*dx10); 
  out[22] += -(1.224744871391589*Ghat[7]*dv10*dx10); 
  out[23] += -(1.8708286933869707*Ghat[1]*dv10*dx10); 
  out[24] += -(1.2247448713915887*Ghat[8]*dv10*dx10); 
  out[25] += -(1.8708286933869707*Ghat[2]*dv10*dx10); 
  out[26] += -(0.7071067811865475*Ghat[10]*dv10*dx10); 
  out[27] += -(1.2247448713915887*Ghat[9]*dv10*dx10); 
  out[28] += -(0.7071067811865475*Ghat[11]*dv10*dx10); 
  out[29] += -(1.8708286933869707*Ghat[3]*dv10*dx10); 
  out[30] += -(1.2247448713915887*Ghat[10]*dv10*dx10); 
  out[31] += -(1.2247448713915887*Ghat[11]*dv10*dx10); 

  } else { 

  if (wv>0) { 

  fUpwind[0] = (1.8708286933869707*fedge[17]+1.5811388300841895*fedge[7]+1.224744871391589*fedge[1]+0.7071067811865475*fedge[0])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[1] = (1.8708286933869707*fedge[23]+1.5811388300841898*fedge[11]+1.224744871391589*fedge[4]+0.7071067811865475*fedge[2])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[2] = (1.8708286933869707*fedge[25]+1.5811388300841898*fedge[13]+1.224744871391589*fedge[5]+0.7071067811865475*fedge[3])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[3] = (1.8708286933869707*fedge[29]+1.5811388300841895*fedge[20]+1.224744871391589*fedge[10]+0.7071067811865475*fedge[6])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[4] = (1.224744871391589*fedge[12]+0.7071067811865475*fedge[8])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[5] = (1.224744871391589*fedge[15]+0.7071067811865475*fedge[9])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[6] = (1.224744871391589*fedge[21]+0.7071067811865475*fedge[14])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[7] = (1.224744871391589*fedge[22]+0.7071067811865475*fedge[16])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[8] = (1.2247448713915887*fedge[24]+0.7071067811865475*fedge[18])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[9] = (1.2247448713915887*fedge[27]+0.7071067811865475*fedge[19])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[10] = (1.2247448713915887*fedge[30]+0.7071067811865475*fedge[26])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[11] = (1.2247448713915887*fedge[31]+0.7071067811865475*fedge[28])*jacob_cx_edge_inv*jacob_vx_inv; 

  } else { 

  fUpwind[0] = (-(1.8708286933869707*fskin[17])+1.5811388300841895*fskin[7]-1.224744871391589*fskin[1]+0.7071067811865475*fskin[0])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[1] = (-(1.8708286933869707*fskin[23])+1.5811388300841898*fskin[11]-1.224744871391589*fskin[4]+0.7071067811865475*fskin[2])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[2] = (-(1.8708286933869707*fskin[25])+1.5811388300841898*fskin[13]-1.224744871391589*fskin[5]+0.7071067811865475*fskin[3])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[3] = (-(1.8708286933869707*fskin[29])+1.5811388300841895*fskin[20]-1.224744871391589*fskin[10]+0.7071067811865475*fskin[6])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[4] = (0.7071067811865475*fskin[8]-1.224744871391589*fskin[12])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[5] = (0.7071067811865475*fskin[9]-1.224744871391589*fskin[15])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[6] = (0.7071067811865475*fskin[14]-1.224744871391589*fskin[21])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[7] = (0.7071067811865475*fskin[16]-1.224744871391589*fskin[22])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[8] = (0.7071067811865475*fskin[18]-1.2247448713915887*fskin[24])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[9] = (0.7071067811865475*fskin[19]-1.2247448713915887*fskin[27])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[10] = (0.7071067811865475*fskin[26]-1.2247448713915887*fskin[30])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[11] = (0.7071067811865475*fskin[28]-1.2247448713915887*fskin[31])*jacob_cx_skin_inv*jacob_vx_inv; 

  } 
  Ghat[0] = 4.183300132670378*hamil[3]*fUpwind[5]+1.8708286933869707*fUpwind[0]*hamil[3]+2.7386127875258306*fUpwind[2]*hamil[2]+1.224744871391589*fUpwind[0]*hamil[1]; 
  Ghat[1] = 4.183300132670379*hamil[3]*fUpwind[7]+1.8708286933869707*fUpwind[1]*hamil[3]+2.7386127875258306*hamil[2]*fUpwind[3]+1.224744871391589*fUpwind[1]*hamil[1]; 
  Ghat[2] = 3.6742346141747664*hamil[3]*fUpwind[9]+2.4494897427831783*hamil[2]*fUpwind[5]+5.612486080160912*fUpwind[2]*hamil[3]+2.7386127875258306*fUpwind[0]*hamil[2]+1.224744871391589*hamil[1]*fUpwind[2]; 
  Ghat[3] = 3.674234614174766*hamil[3]*fUpwind[11]+2.4494897427831783*hamil[2]*fUpwind[7]+5.612486080160912*fUpwind[3]*hamil[3]+1.224744871391589*hamil[1]*fUpwind[3]+2.7386127875258306*fUpwind[1]*hamil[2]; 
  Ghat[4] = 2.7386127875258306*hamil[2]*fUpwind[6]+1.8708286933869707*hamil[3]*fUpwind[4]+1.224744871391589*hamil[1]*fUpwind[4]; 
  Ghat[5] = 2.405351177211819*hamil[2]*fUpwind[9]+4.543441112511214*hamil[3]*fUpwind[5]+1.224744871391589*hamil[1]*fUpwind[5]+4.183300132670378*fUpwind[0]*hamil[3]+2.4494897427831783*fUpwind[2]*hamil[2]; 
  Ghat[6] = 5.612486080160912*hamil[3]*fUpwind[6]+1.224744871391589*hamil[1]*fUpwind[6]+2.7386127875258306*hamil[2]*fUpwind[4]; 
  Ghat[7] = 2.405351177211819*hamil[2]*fUpwind[11]+4.543441112511214*hamil[3]*fUpwind[7]+1.224744871391589*hamil[1]*fUpwind[7]+4.183300132670379*fUpwind[1]*hamil[3]+2.4494897427831783*hamil[2]*fUpwind[3]; 
  Ghat[8] = 2.7386127875258297*hamil[2]*fUpwind[10]+1.8708286933869707*hamil[3]*fUpwind[8]+1.224744871391589*hamil[1]*fUpwind[8]; 
  Ghat[9] = 4.365266951236265*hamil[3]*fUpwind[9]+1.224744871391589*hamil[1]*fUpwind[9]+2.405351177211819*hamil[2]*fUpwind[5]+3.6742346141747664*fUpwind[2]*hamil[3]; 
  Ghat[10] = 5.612486080160912*hamil[3]*fUpwind[10]+1.224744871391589*hamil[1]*fUpwind[10]+2.7386127875258297*hamil[2]*fUpwind[8]; 
  Ghat[11] = 4.365266951236265*hamil[3]*fUpwind[11]+1.224744871391589*hamil[1]*fUpwind[11]+2.405351177211819*hamil[2]*fUpwind[7]+3.674234614174766*fUpwind[3]*hamil[3]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv10*dx10; 
  out[1] += -(1.224744871391589*Ghat[0]*dv10*dx10); 
  out[2] += 0.7071067811865475*Ghat[1]*dv10*dx10; 
  out[3] += 0.7071067811865475*Ghat[2]*dv10*dx10; 
  out[4] += -(1.224744871391589*Ghat[1]*dv10*dx10); 
  out[5] += -(1.224744871391589*Ghat[2]*dv10*dx10); 
  out[6] += 0.7071067811865475*Ghat[3]*dv10*dx10; 
  out[7] += 1.5811388300841895*Ghat[0]*dv10*dx10; 
  out[8] += 0.7071067811865475*Ghat[4]*dv10*dx10; 
  out[9] += 0.7071067811865475*Ghat[5]*dv10*dx10; 
  out[10] += -(1.224744871391589*Ghat[3]*dv10*dx10); 
  out[11] += 1.5811388300841898*Ghat[1]*dv10*dx10; 
  out[12] += -(1.224744871391589*Ghat[4]*dv10*dx10); 
  out[13] += 1.5811388300841898*Ghat[2]*dv10*dx10; 
  out[14] += 0.7071067811865475*Ghat[6]*dv10*dx10; 
  out[15] += -(1.224744871391589*Ghat[5]*dv10*dx10); 
  out[16] += 0.7071067811865475*Ghat[7]*dv10*dx10; 
  out[17] += -(1.8708286933869707*Ghat[0]*dv10*dx10); 
  out[18] += 0.7071067811865475*Ghat[8]*dv10*dx10; 
  out[19] += 0.7071067811865475*Ghat[9]*dv10*dx10; 
  out[20] += 1.5811388300841895*Ghat[3]*dv10*dx10; 
  out[21] += -(1.224744871391589*Ghat[6]*dv10*dx10); 
  out[22] += -(1.224744871391589*Ghat[7]*dv10*dx10); 
  out[23] += -(1.8708286933869707*Ghat[1]*dv10*dx10); 
  out[24] += -(1.2247448713915887*Ghat[8]*dv10*dx10); 
  out[25] += -(1.8708286933869707*Ghat[2]*dv10*dx10); 
  out[26] += 0.7071067811865475*Ghat[10]*dv10*dx10; 
  out[27] += -(1.2247448713915887*Ghat[9]*dv10*dx10); 
  out[28] += 0.7071067811865475*Ghat[11]*dv10*dx10; 
  out[29] += -(1.8708286933869707*Ghat[3]*dv10*dx10); 
  out[30] += -(1.2247448713915887*Ghat[10]*dv10*dx10); 
  out[31] += -(1.2247448713915887*Ghat[11]*dv10*dx10); 

  } 
  return 0.0;

} 
