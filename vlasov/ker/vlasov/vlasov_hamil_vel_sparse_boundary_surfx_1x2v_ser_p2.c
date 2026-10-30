#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_vel_sparse_boundary_surfx_1x2v_ser_p2(const double *w, const double *dxv,
  const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel,
  const double *poisson_tensor_conf, const double *hamil,
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[1]; 
  double wv = w[1]; 

  double fUpwind[8]; 
  double Ghat[8]; 

  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 

  const double *jacob_cx_edge = &jacob_pos_edge[0]; 
  const double *jacob_cx_skin = &jacob_pos_skin[0]; 
  const double jacob_cx_edge_inv = 1.0/jacob_cx_edge[0]; 
  const double jacob_cx_skin_inv = 1.0/jacob_cx_skin[0]; 

  if (edge == -1) { 

  if (wv>0) { 

  fUpwind[0] = (1.5811388300841895*fskin[7]+1.224744871391589*fskin[1]+0.7071067811865475*fskin[0])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[1] = (1.5811388300841898*fskin[11]+1.224744871391589*fskin[4]+0.7071067811865475*fskin[2])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[2] = (1.5811388300841898*fskin[13]+1.224744871391589*fskin[5]+0.7071067811865475*fskin[3])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[3] = (1.5811388300841895*fskin[17]+1.224744871391589*fskin[10]+0.7071067811865475*fskin[6])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[4] = (1.224744871391589*fskin[12]+0.7071067811865475*fskin[8])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[5] = (1.224744871391589*fskin[15]+0.7071067811865475*fskin[9])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[6] = (1.224744871391589*fskin[18]+0.7071067811865475*fskin[14])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[7] = (1.224744871391589*fskin[19]+0.7071067811865475*fskin[16])*jacob_cx_skin_inv*jacob_vx_inv; 

  } else { 

  fUpwind[0] = (1.5811388300841895*fedge[7]-1.224744871391589*fedge[1]+0.7071067811865475*fedge[0])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[1] = (1.5811388300841898*fedge[11]-1.224744871391589*fedge[4]+0.7071067811865475*fedge[2])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[2] = (1.5811388300841898*fedge[13]-1.224744871391589*fedge[5]+0.7071067811865475*fedge[3])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[3] = (1.5811388300841895*fedge[17]-1.224744871391589*fedge[10]+0.7071067811865475*fedge[6])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[4] = (0.7071067811865475*fedge[8]-1.224744871391589*fedge[12])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[5] = (0.7071067811865475*fedge[9]-1.224744871391589*fedge[15])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[6] = (0.7071067811865475*fedge[14]-1.224744871391589*fedge[18])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[7] = (0.7071067811865475*fedge[16]-1.224744871391589*fedge[19])*jacob_cx_edge_inv*jacob_vx_inv; 

  } 
  Ghat[0] = 1.9364916731037085*fUpwind[1]*hamil[4]+0.8660254037844386*fUpwind[0]*hamil[1]; 
  Ghat[1] = 1.7320508075688772*fUpwind[4]*hamil[4]+1.9364916731037085*fUpwind[0]*hamil[4]+0.8660254037844386*fUpwind[1]*hamil[1]; 
  Ghat[2] = 1.9364916731037085*fUpwind[3]*hamil[4]+0.8660254037844386*hamil[1]*fUpwind[2]; 
  Ghat[3] = 1.7320508075688774*hamil[4]*fUpwind[6]+1.9364916731037085*fUpwind[2]*hamil[4]+0.8660254037844386*hamil[1]*fUpwind[3]; 
  Ghat[4] = 1.7320508075688772*fUpwind[1]*hamil[4]+0.8660254037844386*hamil[1]*fUpwind[4]; 
  Ghat[5] = 1.9364916731037085*hamil[4]*fUpwind[7]+0.8660254037844386*hamil[1]*fUpwind[5]; 
  Ghat[6] = 0.8660254037844386*hamil[1]*fUpwind[6]+1.7320508075688774*fUpwind[3]*hamil[4]; 
  Ghat[7] = 0.8660254037844386*hamil[1]*fUpwind[7]+1.9364916731037085*hamil[4]*fUpwind[5]; 

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
  out[17] += -(1.5811388300841895*Ghat[3]*dv10*dx10); 
  out[18] += -(1.224744871391589*Ghat[6]*dv10*dx10); 
  out[19] += -(1.224744871391589*Ghat[7]*dv10*dx10); 

  } else { 

  if (wv>0) { 

  fUpwind[0] = (1.5811388300841895*fedge[7]+1.224744871391589*fedge[1]+0.7071067811865475*fedge[0])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[1] = (1.5811388300841898*fedge[11]+1.224744871391589*fedge[4]+0.7071067811865475*fedge[2])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[2] = (1.5811388300841898*fedge[13]+1.224744871391589*fedge[5]+0.7071067811865475*fedge[3])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[3] = (1.5811388300841895*fedge[17]+1.224744871391589*fedge[10]+0.7071067811865475*fedge[6])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[4] = (1.224744871391589*fedge[12]+0.7071067811865475*fedge[8])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[5] = (1.224744871391589*fedge[15]+0.7071067811865475*fedge[9])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[6] = (1.224744871391589*fedge[18]+0.7071067811865475*fedge[14])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[7] = (1.224744871391589*fedge[19]+0.7071067811865475*fedge[16])*jacob_cx_edge_inv*jacob_vx_inv; 

  } else { 

  fUpwind[0] = (1.5811388300841895*fskin[7]-1.224744871391589*fskin[1]+0.7071067811865475*fskin[0])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[1] = (1.5811388300841898*fskin[11]-1.224744871391589*fskin[4]+0.7071067811865475*fskin[2])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[2] = (1.5811388300841898*fskin[13]-1.224744871391589*fskin[5]+0.7071067811865475*fskin[3])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[3] = (1.5811388300841895*fskin[17]-1.224744871391589*fskin[10]+0.7071067811865475*fskin[6])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[4] = (0.7071067811865475*fskin[8]-1.224744871391589*fskin[12])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[5] = (0.7071067811865475*fskin[9]-1.224744871391589*fskin[15])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[6] = (0.7071067811865475*fskin[14]-1.224744871391589*fskin[18])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[7] = (0.7071067811865475*fskin[16]-1.224744871391589*fskin[19])*jacob_cx_skin_inv*jacob_vx_inv; 

  } 
  Ghat[0] = 1.9364916731037085*fUpwind[1]*hamil[4]+0.8660254037844386*fUpwind[0]*hamil[1]; 
  Ghat[1] = 1.7320508075688772*fUpwind[4]*hamil[4]+1.9364916731037085*fUpwind[0]*hamil[4]+0.8660254037844386*fUpwind[1]*hamil[1]; 
  Ghat[2] = 1.9364916731037085*fUpwind[3]*hamil[4]+0.8660254037844386*hamil[1]*fUpwind[2]; 
  Ghat[3] = 1.7320508075688774*hamil[4]*fUpwind[6]+1.9364916731037085*fUpwind[2]*hamil[4]+0.8660254037844386*hamil[1]*fUpwind[3]; 
  Ghat[4] = 1.7320508075688772*fUpwind[1]*hamil[4]+0.8660254037844386*hamil[1]*fUpwind[4]; 
  Ghat[5] = 1.9364916731037085*hamil[4]*fUpwind[7]+0.8660254037844386*hamil[1]*fUpwind[5]; 
  Ghat[6] = 0.8660254037844386*hamil[1]*fUpwind[6]+1.7320508075688774*fUpwind[3]*hamil[4]; 
  Ghat[7] = 0.8660254037844386*hamil[1]*fUpwind[7]+1.9364916731037085*hamil[4]*fUpwind[5]; 

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
  out[17] += 1.5811388300841895*Ghat[3]*dv10*dx10; 
  out[18] += -(1.224744871391589*Ghat[6]*dv10*dx10); 
  out[19] += -(1.224744871391589*Ghat[7]*dv10*dx10); 

  } 
  return 0.0;

} 
