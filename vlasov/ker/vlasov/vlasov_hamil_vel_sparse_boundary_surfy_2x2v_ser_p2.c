#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_vel_sparse_boundary_surfy_2x2v_ser_p2(const double *w, const double *dxv,
  const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel,
  const double *poisson_tensor_conf, const double *hamil,
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  double dx11 = 2.0/dxv[1]; 
  double dv11 = 2.0/dxv[3]; 
  double wv = w[3]; 

  double fUpwind[20]; 
  double Ghat[20]; 

  const double *jacob_vy = &jacob_vel[3]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 

  const double *jacob_cy_edge = &jacob_pos_edge[3]; 
  const double *jacob_cy_skin = &jacob_pos_skin[3]; 
  const double jacob_cy_edge_inv = 1.0/jacob_cy_edge[0]; 
  const double jacob_cy_skin_inv = 1.0/jacob_cy_skin[0]; 

  if (edge == -1) { 

  if (wv>0) { 

  fUpwind[0] = (1.5811388300841895*fskin[12]+1.224744871391589*fskin[2]+0.7071067811865475*fskin[0])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[1] = (1.5811388300841898*fskin[20]+1.224744871391589*fskin[5]+0.7071067811865475*fskin[1])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[2] = (1.5811388300841898*fskin[22]+1.224744871391589*fskin[7]+0.7071067811865475*fskin[3])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[3] = (1.5811388300841898*fskin[26]+1.224744871391589*fskin[9]+0.7071067811865475*fskin[4])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[4] = (1.5811388300841895*fskin[33]+1.224744871391589*fskin[15]+0.7071067811865475*fskin[6])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[5] = (1.5811388300841895*fskin[36]+1.224744871391589*fskin[16]+0.7071067811865475*fskin[8])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[6] = (1.5811388300841895*fskin[38]+1.224744871391589*fskin[18]+0.7071067811865475*fskin[10])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[7] = (1.224744871391589*fskin[19]+0.7071067811865475*fskin[11])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[8] = (1.224744871391589*fskin[24]+0.7071067811865475*fskin[13])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[9] = (1.224744871391589*fskin[29]+0.7071067811865475*fskin[14])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[10] = (1.5811388300841898*fskin[45]+1.224744871391589*fskin[31]+0.7071067811865475*fskin[17])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[11] = (1.224744871391589*fskin[32]+0.7071067811865475*fskin[21])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[12] = (1.224744871391589*fskin[34]+0.7071067811865475*fskin[23])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[13] = (1.224744871391589*fskin[35]+0.7071067811865475*fskin[25])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[14] = (1.224744871391589*fskin[40]+0.7071067811865475*fskin[27])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[15] = (1.224744871391589*fskin[41]+0.7071067811865475*fskin[28])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[16] = (1.224744871391589*fskin[43]+0.7071067811865475*fskin[30])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[17] = (1.224744871391589*fskin[44]+0.7071067811865475*fskin[37])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[18] = (1.224744871391589*fskin[46]+0.7071067811865475*fskin[39])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[19] = (1.224744871391589*fskin[47]+0.7071067811865475*fskin[42])*jacob_cy_skin_inv*jacob_vy_inv; 

  } else { 

  fUpwind[0] = (1.5811388300841895*fedge[12]-1.224744871391589*fedge[2]+0.7071067811865475*fedge[0])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[1] = (1.5811388300841898*fedge[20]-1.224744871391589*fedge[5]+0.7071067811865475*fedge[1])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[2] = (1.5811388300841898*fedge[22]-1.224744871391589*fedge[7]+0.7071067811865475*fedge[3])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[3] = (1.5811388300841898*fedge[26]-1.224744871391589*fedge[9]+0.7071067811865475*fedge[4])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[4] = (1.5811388300841895*fedge[33]-1.224744871391589*fedge[15]+0.7071067811865475*fedge[6])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[5] = (1.5811388300841895*fedge[36]-1.224744871391589*fedge[16]+0.7071067811865475*fedge[8])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[6] = (1.5811388300841895*fedge[38]-1.224744871391589*fedge[18]+0.7071067811865475*fedge[10])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[7] = (0.7071067811865475*fedge[11]-1.224744871391589*fedge[19])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[8] = (0.7071067811865475*fedge[13]-1.224744871391589*fedge[24])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[9] = (0.7071067811865475*fedge[14]-1.224744871391589*fedge[29])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[10] = (1.5811388300841898*fedge[45]-1.224744871391589*fedge[31]+0.7071067811865475*fedge[17])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[11] = (0.7071067811865475*fedge[21]-1.224744871391589*fedge[32])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[12] = (0.7071067811865475*fedge[23]-1.224744871391589*fedge[34])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[13] = (0.7071067811865475*fedge[25]-1.224744871391589*fedge[35])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[14] = (0.7071067811865475*fedge[27]-1.224744871391589*fedge[40])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[15] = (0.7071067811865475*fedge[28]-1.224744871391589*fedge[41])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[16] = (0.7071067811865475*fedge[30]-1.224744871391589*fedge[43])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[17] = (0.7071067811865475*fedge[37]-1.224744871391589*fedge[44])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[18] = (0.7071067811865475*fedge[39]-1.224744871391589*fedge[46])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[19] = (0.7071067811865475*fedge[42]-1.224744871391589*fedge[47])*jacob_cy_edge_inv*jacob_vy_inv; 

  } 
  Ghat[0] = 1.9364916731037085*fUpwind[3]*hamil[5]+0.8660254037844386*fUpwind[0]*hamil[2]; 
  Ghat[1] = 1.9364916731037085*fUpwind[5]*hamil[5]+0.8660254037844386*fUpwind[1]*hamil[2]; 
  Ghat[2] = 1.9364916731037085*hamil[5]*fUpwind[6]+0.8660254037844386*fUpwind[2]*hamil[2]; 
  Ghat[3] = 1.7320508075688772*hamil[5]*fUpwind[9]+1.9364916731037085*fUpwind[0]*hamil[5]+0.8660254037844386*hamil[2]*fUpwind[3]; 
  Ghat[4] = 1.9364916731037085*hamil[5]*fUpwind[10]+0.8660254037844386*hamil[2]*fUpwind[4]; 
  Ghat[5] = 1.7320508075688774*hamil[5]*fUpwind[15]+1.9364916731037085*fUpwind[1]*hamil[5]+0.8660254037844386*hamil[2]*fUpwind[5]; 
  Ghat[6] = 1.7320508075688774*hamil[5]*fUpwind[16]+0.8660254037844386*hamil[2]*fUpwind[6]+1.9364916731037085*fUpwind[2]*hamil[5]; 
  Ghat[7] = 1.9364916731037085*hamil[5]*fUpwind[13]+0.8660254037844386*hamil[2]*fUpwind[7]; 
  Ghat[8] = 1.9364916731037085*hamil[5]*fUpwind[14]+0.8660254037844386*hamil[2]*fUpwind[8]; 
  Ghat[9] = 0.8660254037844386*hamil[2]*fUpwind[9]+1.7320508075688772*fUpwind[3]*hamil[5]; 
  Ghat[10] = 1.7320508075688772*hamil[5]*fUpwind[19]+0.8660254037844386*hamil[2]*fUpwind[10]+1.9364916731037085*fUpwind[4]*hamil[5]; 
  Ghat[11] = 1.9364916731037085*hamil[5]*fUpwind[17]+0.8660254037844386*hamil[2]*fUpwind[11]; 
  Ghat[12] = 1.9364916731037085*hamil[5]*fUpwind[18]+0.8660254037844386*hamil[2]*fUpwind[12]; 
  Ghat[13] = 0.8660254037844386*hamil[2]*fUpwind[13]+1.9364916731037085*hamil[5]*fUpwind[7]; 
  Ghat[14] = 0.8660254037844386*hamil[2]*fUpwind[14]+1.9364916731037085*hamil[5]*fUpwind[8]; 
  Ghat[15] = 0.8660254037844386*hamil[2]*fUpwind[15]+1.7320508075688774*fUpwind[5]*hamil[5]; 
  Ghat[16] = 0.8660254037844386*hamil[2]*fUpwind[16]+1.7320508075688774*hamil[5]*fUpwind[6]; 
  Ghat[17] = 0.8660254037844386*hamil[2]*fUpwind[17]+1.9364916731037085*hamil[5]*fUpwind[11]; 
  Ghat[18] = 0.8660254037844386*hamil[2]*fUpwind[18]+1.9364916731037085*hamil[5]*fUpwind[12]; 
  Ghat[19] = 0.8660254037844386*hamil[2]*fUpwind[19]+1.7320508075688772*hamil[5]*fUpwind[10]; 

  out[0] += -(0.7071067811865475*Ghat[0]*dv11*dx11); 
  out[1] += -(0.7071067811865475*Ghat[1]*dv11*dx11); 
  out[2] += -(1.224744871391589*Ghat[0]*dv11*dx11); 
  out[3] += -(0.7071067811865475*Ghat[2]*dv11*dx11); 
  out[4] += -(0.7071067811865475*Ghat[3]*dv11*dx11); 
  out[5] += -(1.224744871391589*Ghat[1]*dv11*dx11); 
  out[6] += -(0.7071067811865475*Ghat[4]*dv11*dx11); 
  out[7] += -(1.224744871391589*Ghat[2]*dv11*dx11); 
  out[8] += -(0.7071067811865475*Ghat[5]*dv11*dx11); 
  out[9] += -(1.224744871391589*Ghat[3]*dv11*dx11); 
  out[10] += -(0.7071067811865475*Ghat[6]*dv11*dx11); 
  out[11] += -(0.7071067811865475*Ghat[7]*dv11*dx11); 
  out[12] += -(1.5811388300841895*Ghat[0]*dv11*dx11); 
  out[13] += -(0.7071067811865475*Ghat[8]*dv11*dx11); 
  out[14] += -(0.7071067811865475*Ghat[9]*dv11*dx11); 
  out[15] += -(1.224744871391589*Ghat[4]*dv11*dx11); 
  out[16] += -(1.224744871391589*Ghat[5]*dv11*dx11); 
  out[17] += -(0.7071067811865475*Ghat[10]*dv11*dx11); 
  out[18] += -(1.224744871391589*Ghat[6]*dv11*dx11); 
  out[19] += -(1.224744871391589*Ghat[7]*dv11*dx11); 
  out[20] += -(1.5811388300841898*Ghat[1]*dv11*dx11); 
  out[21] += -(0.7071067811865475*Ghat[11]*dv11*dx11); 
  out[22] += -(1.5811388300841898*Ghat[2]*dv11*dx11); 
  out[23] += -(0.7071067811865475*Ghat[12]*dv11*dx11); 
  out[24] += -(1.224744871391589*Ghat[8]*dv11*dx11); 
  out[25] += -(0.7071067811865475*Ghat[13]*dv11*dx11); 
  out[26] += -(1.5811388300841898*Ghat[3]*dv11*dx11); 
  out[27] += -(0.7071067811865475*Ghat[14]*dv11*dx11); 
  out[28] += -(0.7071067811865475*Ghat[15]*dv11*dx11); 
  out[29] += -(1.224744871391589*Ghat[9]*dv11*dx11); 
  out[30] += -(0.7071067811865475*Ghat[16]*dv11*dx11); 
  out[31] += -(1.224744871391589*Ghat[10]*dv11*dx11); 
  out[32] += -(1.224744871391589*Ghat[11]*dv11*dx11); 
  out[33] += -(1.5811388300841895*Ghat[4]*dv11*dx11); 
  out[34] += -(1.224744871391589*Ghat[12]*dv11*dx11); 
  out[35] += -(1.224744871391589*Ghat[13]*dv11*dx11); 
  out[36] += -(1.5811388300841895*Ghat[5]*dv11*dx11); 
  out[37] += -(0.7071067811865475*Ghat[17]*dv11*dx11); 
  out[38] += -(1.5811388300841895*Ghat[6]*dv11*dx11); 
  out[39] += -(0.7071067811865475*Ghat[18]*dv11*dx11); 
  out[40] += -(1.224744871391589*Ghat[14]*dv11*dx11); 
  out[41] += -(1.224744871391589*Ghat[15]*dv11*dx11); 
  out[42] += -(0.7071067811865475*Ghat[19]*dv11*dx11); 
  out[43] += -(1.224744871391589*Ghat[16]*dv11*dx11); 
  out[44] += -(1.224744871391589*Ghat[17]*dv11*dx11); 
  out[45] += -(1.5811388300841898*Ghat[10]*dv11*dx11); 
  out[46] += -(1.224744871391589*Ghat[18]*dv11*dx11); 
  out[47] += -(1.224744871391589*Ghat[19]*dv11*dx11); 

  } else { 

  if (wv>0) { 

  fUpwind[0] = (1.5811388300841895*fedge[12]+1.224744871391589*fedge[2]+0.7071067811865475*fedge[0])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[1] = (1.5811388300841898*fedge[20]+1.224744871391589*fedge[5]+0.7071067811865475*fedge[1])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[2] = (1.5811388300841898*fedge[22]+1.224744871391589*fedge[7]+0.7071067811865475*fedge[3])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[3] = (1.5811388300841898*fedge[26]+1.224744871391589*fedge[9]+0.7071067811865475*fedge[4])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[4] = (1.5811388300841895*fedge[33]+1.224744871391589*fedge[15]+0.7071067811865475*fedge[6])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[5] = (1.5811388300841895*fedge[36]+1.224744871391589*fedge[16]+0.7071067811865475*fedge[8])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[6] = (1.5811388300841895*fedge[38]+1.224744871391589*fedge[18]+0.7071067811865475*fedge[10])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[7] = (1.224744871391589*fedge[19]+0.7071067811865475*fedge[11])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[8] = (1.224744871391589*fedge[24]+0.7071067811865475*fedge[13])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[9] = (1.224744871391589*fedge[29]+0.7071067811865475*fedge[14])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[10] = (1.5811388300841898*fedge[45]+1.224744871391589*fedge[31]+0.7071067811865475*fedge[17])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[11] = (1.224744871391589*fedge[32]+0.7071067811865475*fedge[21])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[12] = (1.224744871391589*fedge[34]+0.7071067811865475*fedge[23])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[13] = (1.224744871391589*fedge[35]+0.7071067811865475*fedge[25])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[14] = (1.224744871391589*fedge[40]+0.7071067811865475*fedge[27])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[15] = (1.224744871391589*fedge[41]+0.7071067811865475*fedge[28])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[16] = (1.224744871391589*fedge[43]+0.7071067811865475*fedge[30])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[17] = (1.224744871391589*fedge[44]+0.7071067811865475*fedge[37])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[18] = (1.224744871391589*fedge[46]+0.7071067811865475*fedge[39])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[19] = (1.224744871391589*fedge[47]+0.7071067811865475*fedge[42])*jacob_cy_edge_inv*jacob_vy_inv; 

  } else { 

  fUpwind[0] = (1.5811388300841895*fskin[12]-1.224744871391589*fskin[2]+0.7071067811865475*fskin[0])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[1] = (1.5811388300841898*fskin[20]-1.224744871391589*fskin[5]+0.7071067811865475*fskin[1])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[2] = (1.5811388300841898*fskin[22]-1.224744871391589*fskin[7]+0.7071067811865475*fskin[3])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[3] = (1.5811388300841898*fskin[26]-1.224744871391589*fskin[9]+0.7071067811865475*fskin[4])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[4] = (1.5811388300841895*fskin[33]-1.224744871391589*fskin[15]+0.7071067811865475*fskin[6])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[5] = (1.5811388300841895*fskin[36]-1.224744871391589*fskin[16]+0.7071067811865475*fskin[8])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[6] = (1.5811388300841895*fskin[38]-1.224744871391589*fskin[18]+0.7071067811865475*fskin[10])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[7] = (0.7071067811865475*fskin[11]-1.224744871391589*fskin[19])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[8] = (0.7071067811865475*fskin[13]-1.224744871391589*fskin[24])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[9] = (0.7071067811865475*fskin[14]-1.224744871391589*fskin[29])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[10] = (1.5811388300841898*fskin[45]-1.224744871391589*fskin[31]+0.7071067811865475*fskin[17])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[11] = (0.7071067811865475*fskin[21]-1.224744871391589*fskin[32])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[12] = (0.7071067811865475*fskin[23]-1.224744871391589*fskin[34])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[13] = (0.7071067811865475*fskin[25]-1.224744871391589*fskin[35])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[14] = (0.7071067811865475*fskin[27]-1.224744871391589*fskin[40])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[15] = (0.7071067811865475*fskin[28]-1.224744871391589*fskin[41])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[16] = (0.7071067811865475*fskin[30]-1.224744871391589*fskin[43])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[17] = (0.7071067811865475*fskin[37]-1.224744871391589*fskin[44])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[18] = (0.7071067811865475*fskin[39]-1.224744871391589*fskin[46])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[19] = (0.7071067811865475*fskin[42]-1.224744871391589*fskin[47])*jacob_cy_skin_inv*jacob_vy_inv; 

  } 
  Ghat[0] = 1.9364916731037085*fUpwind[3]*hamil[5]+0.8660254037844386*fUpwind[0]*hamil[2]; 
  Ghat[1] = 1.9364916731037085*fUpwind[5]*hamil[5]+0.8660254037844386*fUpwind[1]*hamil[2]; 
  Ghat[2] = 1.9364916731037085*hamil[5]*fUpwind[6]+0.8660254037844386*fUpwind[2]*hamil[2]; 
  Ghat[3] = 1.7320508075688772*hamil[5]*fUpwind[9]+1.9364916731037085*fUpwind[0]*hamil[5]+0.8660254037844386*hamil[2]*fUpwind[3]; 
  Ghat[4] = 1.9364916731037085*hamil[5]*fUpwind[10]+0.8660254037844386*hamil[2]*fUpwind[4]; 
  Ghat[5] = 1.7320508075688774*hamil[5]*fUpwind[15]+1.9364916731037085*fUpwind[1]*hamil[5]+0.8660254037844386*hamil[2]*fUpwind[5]; 
  Ghat[6] = 1.7320508075688774*hamil[5]*fUpwind[16]+0.8660254037844386*hamil[2]*fUpwind[6]+1.9364916731037085*fUpwind[2]*hamil[5]; 
  Ghat[7] = 1.9364916731037085*hamil[5]*fUpwind[13]+0.8660254037844386*hamil[2]*fUpwind[7]; 
  Ghat[8] = 1.9364916731037085*hamil[5]*fUpwind[14]+0.8660254037844386*hamil[2]*fUpwind[8]; 
  Ghat[9] = 0.8660254037844386*hamil[2]*fUpwind[9]+1.7320508075688772*fUpwind[3]*hamil[5]; 
  Ghat[10] = 1.7320508075688772*hamil[5]*fUpwind[19]+0.8660254037844386*hamil[2]*fUpwind[10]+1.9364916731037085*fUpwind[4]*hamil[5]; 
  Ghat[11] = 1.9364916731037085*hamil[5]*fUpwind[17]+0.8660254037844386*hamil[2]*fUpwind[11]; 
  Ghat[12] = 1.9364916731037085*hamil[5]*fUpwind[18]+0.8660254037844386*hamil[2]*fUpwind[12]; 
  Ghat[13] = 0.8660254037844386*hamil[2]*fUpwind[13]+1.9364916731037085*hamil[5]*fUpwind[7]; 
  Ghat[14] = 0.8660254037844386*hamil[2]*fUpwind[14]+1.9364916731037085*hamil[5]*fUpwind[8]; 
  Ghat[15] = 0.8660254037844386*hamil[2]*fUpwind[15]+1.7320508075688774*fUpwind[5]*hamil[5]; 
  Ghat[16] = 0.8660254037844386*hamil[2]*fUpwind[16]+1.7320508075688774*hamil[5]*fUpwind[6]; 
  Ghat[17] = 0.8660254037844386*hamil[2]*fUpwind[17]+1.9364916731037085*hamil[5]*fUpwind[11]; 
  Ghat[18] = 0.8660254037844386*hamil[2]*fUpwind[18]+1.9364916731037085*hamil[5]*fUpwind[12]; 
  Ghat[19] = 0.8660254037844386*hamil[2]*fUpwind[19]+1.7320508075688772*hamil[5]*fUpwind[10]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv11*dx11; 
  out[1] += 0.7071067811865475*Ghat[1]*dv11*dx11; 
  out[2] += -(1.224744871391589*Ghat[0]*dv11*dx11); 
  out[3] += 0.7071067811865475*Ghat[2]*dv11*dx11; 
  out[4] += 0.7071067811865475*Ghat[3]*dv11*dx11; 
  out[5] += -(1.224744871391589*Ghat[1]*dv11*dx11); 
  out[6] += 0.7071067811865475*Ghat[4]*dv11*dx11; 
  out[7] += -(1.224744871391589*Ghat[2]*dv11*dx11); 
  out[8] += 0.7071067811865475*Ghat[5]*dv11*dx11; 
  out[9] += -(1.224744871391589*Ghat[3]*dv11*dx11); 
  out[10] += 0.7071067811865475*Ghat[6]*dv11*dx11; 
  out[11] += 0.7071067811865475*Ghat[7]*dv11*dx11; 
  out[12] += 1.5811388300841895*Ghat[0]*dv11*dx11; 
  out[13] += 0.7071067811865475*Ghat[8]*dv11*dx11; 
  out[14] += 0.7071067811865475*Ghat[9]*dv11*dx11; 
  out[15] += -(1.224744871391589*Ghat[4]*dv11*dx11); 
  out[16] += -(1.224744871391589*Ghat[5]*dv11*dx11); 
  out[17] += 0.7071067811865475*Ghat[10]*dv11*dx11; 
  out[18] += -(1.224744871391589*Ghat[6]*dv11*dx11); 
  out[19] += -(1.224744871391589*Ghat[7]*dv11*dx11); 
  out[20] += 1.5811388300841898*Ghat[1]*dv11*dx11; 
  out[21] += 0.7071067811865475*Ghat[11]*dv11*dx11; 
  out[22] += 1.5811388300841898*Ghat[2]*dv11*dx11; 
  out[23] += 0.7071067811865475*Ghat[12]*dv11*dx11; 
  out[24] += -(1.224744871391589*Ghat[8]*dv11*dx11); 
  out[25] += 0.7071067811865475*Ghat[13]*dv11*dx11; 
  out[26] += 1.5811388300841898*Ghat[3]*dv11*dx11; 
  out[27] += 0.7071067811865475*Ghat[14]*dv11*dx11; 
  out[28] += 0.7071067811865475*Ghat[15]*dv11*dx11; 
  out[29] += -(1.224744871391589*Ghat[9]*dv11*dx11); 
  out[30] += 0.7071067811865475*Ghat[16]*dv11*dx11; 
  out[31] += -(1.224744871391589*Ghat[10]*dv11*dx11); 
  out[32] += -(1.224744871391589*Ghat[11]*dv11*dx11); 
  out[33] += 1.5811388300841895*Ghat[4]*dv11*dx11; 
  out[34] += -(1.224744871391589*Ghat[12]*dv11*dx11); 
  out[35] += -(1.224744871391589*Ghat[13]*dv11*dx11); 
  out[36] += 1.5811388300841895*Ghat[5]*dv11*dx11; 
  out[37] += 0.7071067811865475*Ghat[17]*dv11*dx11; 
  out[38] += 1.5811388300841895*Ghat[6]*dv11*dx11; 
  out[39] += 0.7071067811865475*Ghat[18]*dv11*dx11; 
  out[40] += -(1.224744871391589*Ghat[14]*dv11*dx11); 
  out[41] += -(1.224744871391589*Ghat[15]*dv11*dx11); 
  out[42] += 0.7071067811865475*Ghat[19]*dv11*dx11; 
  out[43] += -(1.224744871391589*Ghat[16]*dv11*dx11); 
  out[44] += -(1.224744871391589*Ghat[17]*dv11*dx11); 
  out[45] += 1.5811388300841898*Ghat[10]*dv11*dx11; 
  out[46] += -(1.224744871391589*Ghat[18]*dv11*dx11); 
  out[47] += -(1.224744871391589*Ghat[19]*dv11*dx11); 

  } 
  return 0.0;

} 
