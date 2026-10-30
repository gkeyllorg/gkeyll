#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfx_1x1v_ser_p3(const double *w, const double *dxv,
  const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel,
  const double *poisson_tensor_conf, const double *hamil,
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[1]; 
  double wv = w[1]; 

  double fUpwind[4]; 
  double Ghat[4]; 

  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 

  const double *jacob_cx_edge = &jacob_pos_edge[0]; 
  const double *jacob_cx_skin = &jacob_pos_skin[0]; 
  const double jacob_cx_edge_inv = 1.0/jacob_cx_edge[0]; 
  const double jacob_cx_skin_inv = 1.0/jacob_cx_skin[0]; 

  if (edge == -1) { 

  if (wv>0) { 

  fUpwind[0] = (1.8708286933869707*fskin[8]+1.5811388300841895*fskin[4]+1.224744871391589*fskin[1]+0.7071067811865475*fskin[0])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[1] = (1.8708286933869707*fskin[10]+1.5811388300841898*fskin[6]+1.224744871391589*fskin[3]+0.7071067811865475*fskin[2])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[2] = (1.224744871391589*fskin[7]+0.7071067811865475*fskin[5])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[3] = (1.2247448713915887*fskin[11]+0.7071067811865475*fskin[9])*jacob_cx_skin_inv*jacob_vx_inv; 

  } else { 

  fUpwind[0] = (-(1.8708286933869707*fedge[8])+1.5811388300841895*fedge[4]-1.224744871391589*fedge[1]+0.7071067811865475*fedge[0])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[1] = (-(1.8708286933869707*fedge[10])+1.5811388300841898*fedge[6]-1.224744871391589*fedge[3]+0.7071067811865475*fedge[2])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[2] = (0.7071067811865475*fedge[5]-1.224744871391589*fedge[7])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[3] = (0.7071067811865475*fedge[9]-1.2247448713915887*fedge[11])*jacob_cx_edge_inv*jacob_vx_inv; 

  } 
  Ghat[0] = 4.183300132670378*fUpwind[2]*hamil[3]+1.8708286933869707*fUpwind[0]*hamil[3]+2.7386127875258306*fUpwind[1]*hamil[2]+1.224744871391589*fUpwind[0]*hamil[1]; 
  Ghat[1] = 3.6742346141747664*fUpwind[3]*hamil[3]+5.612486080160912*fUpwind[1]*hamil[3]+2.4494897427831783*fUpwind[2]*hamil[2]+2.7386127875258306*fUpwind[0]*hamil[2]+1.224744871391589*fUpwind[1]*hamil[1]; 
  Ghat[2] = 4.543441112511214*fUpwind[2]*hamil[3]+4.183300132670378*fUpwind[0]*hamil[3]+2.405351177211819*hamil[2]*fUpwind[3]+2.4494897427831783*fUpwind[1]*hamil[2]+1.224744871391589*hamil[1]*fUpwind[2]; 
  Ghat[3] = 4.365266951236265*fUpwind[3]*hamil[3]+3.6742346141747664*fUpwind[1]*hamil[3]+1.224744871391589*hamil[1]*fUpwind[3]+2.405351177211819*fUpwind[2]*hamil[2]; 

  out[0] += -(0.7071067811865475*Ghat[0]*dv10*dx10); 
  out[1] += -(1.224744871391589*Ghat[0]*dv10*dx10); 
  out[2] += -(0.7071067811865475*Ghat[1]*dv10*dx10); 
  out[3] += -(1.224744871391589*Ghat[1]*dv10*dx10); 
  out[4] += -(1.5811388300841895*Ghat[0]*dv10*dx10); 
  out[5] += -(0.7071067811865475*Ghat[2]*dv10*dx10); 
  out[6] += -(1.5811388300841898*Ghat[1]*dv10*dx10); 
  out[7] += -(1.224744871391589*Ghat[2]*dv10*dx10); 
  out[8] += -(1.8708286933869707*Ghat[0]*dv10*dx10); 
  out[9] += -(0.7071067811865475*Ghat[3]*dv10*dx10); 
  out[10] += -(1.8708286933869707*Ghat[1]*dv10*dx10); 
  out[11] += -(1.2247448713915887*Ghat[3]*dv10*dx10); 

  } else { 

  if (wv>0) { 

  fUpwind[0] = (1.8708286933869707*fedge[8]+1.5811388300841895*fedge[4]+1.224744871391589*fedge[1]+0.7071067811865475*fedge[0])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[1] = (1.8708286933869707*fedge[10]+1.5811388300841898*fedge[6]+1.224744871391589*fedge[3]+0.7071067811865475*fedge[2])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[2] = (1.224744871391589*fedge[7]+0.7071067811865475*fedge[5])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[3] = (1.2247448713915887*fedge[11]+0.7071067811865475*fedge[9])*jacob_cx_edge_inv*jacob_vx_inv; 

  } else { 

  fUpwind[0] = (-(1.8708286933869707*fskin[8])+1.5811388300841895*fskin[4]-1.224744871391589*fskin[1]+0.7071067811865475*fskin[0])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[1] = (-(1.8708286933869707*fskin[10])+1.5811388300841898*fskin[6]-1.224744871391589*fskin[3]+0.7071067811865475*fskin[2])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[2] = (0.7071067811865475*fskin[5]-1.224744871391589*fskin[7])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[3] = (0.7071067811865475*fskin[9]-1.2247448713915887*fskin[11])*jacob_cx_skin_inv*jacob_vx_inv; 

  } 
  Ghat[0] = 4.183300132670378*fUpwind[2]*hamil[3]+1.8708286933869707*fUpwind[0]*hamil[3]+2.7386127875258306*fUpwind[1]*hamil[2]+1.224744871391589*fUpwind[0]*hamil[1]; 
  Ghat[1] = 3.6742346141747664*fUpwind[3]*hamil[3]+5.612486080160912*fUpwind[1]*hamil[3]+2.4494897427831783*fUpwind[2]*hamil[2]+2.7386127875258306*fUpwind[0]*hamil[2]+1.224744871391589*fUpwind[1]*hamil[1]; 
  Ghat[2] = 4.543441112511214*fUpwind[2]*hamil[3]+4.183300132670378*fUpwind[0]*hamil[3]+2.405351177211819*hamil[2]*fUpwind[3]+2.4494897427831783*fUpwind[1]*hamil[2]+1.224744871391589*hamil[1]*fUpwind[2]; 
  Ghat[3] = 4.365266951236265*fUpwind[3]*hamil[3]+3.6742346141747664*fUpwind[1]*hamil[3]+1.224744871391589*hamil[1]*fUpwind[3]+2.405351177211819*fUpwind[2]*hamil[2]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv10*dx10; 
  out[1] += -(1.224744871391589*Ghat[0]*dv10*dx10); 
  out[2] += 0.7071067811865475*Ghat[1]*dv10*dx10; 
  out[3] += -(1.224744871391589*Ghat[1]*dv10*dx10); 
  out[4] += 1.5811388300841895*Ghat[0]*dv10*dx10; 
  out[5] += 0.7071067811865475*Ghat[2]*dv10*dx10; 
  out[6] += 1.5811388300841898*Ghat[1]*dv10*dx10; 
  out[7] += -(1.224744871391589*Ghat[2]*dv10*dx10); 
  out[8] += -(1.8708286933869707*Ghat[0]*dv10*dx10); 
  out[9] += 0.7071067811865475*Ghat[3]*dv10*dx10; 
  out[10] += -(1.8708286933869707*Ghat[1]*dv10*dx10); 
  out[11] += -(1.2247448713915887*Ghat[3]*dv10*dx10); 

  } 
  return 0.0;

} 
