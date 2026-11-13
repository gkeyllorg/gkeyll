#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_vel_sparse_boundary_surfx_2x2v_ser_p2(const double *w, const double *dxv,
  const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel,
  const double *poisson_tensor_conf, const double *hamil,
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[2]; 
  double wv = w[2]; 

  double fUpwind[20]; 
  double Ghat[20]; 

  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 

  const double *jacob_cx_edge = &jacob_pos_edge[0]; 
  const double *jacob_cx_skin = &jacob_pos_skin[0]; 
  const double jacob_cx_edge_inv = 1.0/jacob_cx_edge[0]; 
  const double jacob_cx_skin_inv = 1.0/jacob_cx_skin[0]; 

  if (edge == -1) { 

  if (wv>0) { 

  fUpwind[0] = (1.5811388300841895*fskin[11]+1.224744871391589*fskin[1]+0.7071067811865475*fskin[0])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[1] = (1.5811388300841898*fskin[19]+1.224744871391589*fskin[5]+0.7071067811865475*fskin[2])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[2] = (1.5811388300841898*fskin[21]+1.224744871391589*fskin[6]+0.7071067811865475*fskin[3])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[3] = (1.5811388300841898*fskin[25]+1.224744871391589*fskin[8]+0.7071067811865475*fskin[4])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[4] = (1.5811388300841895*fskin[32]+1.224744871391589*fskin[15]+0.7071067811865475*fskin[7])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[5] = (1.5811388300841895*fskin[35]+1.224744871391589*fskin[16]+0.7071067811865475*fskin[9])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[6] = (1.5811388300841895*fskin[37]+1.224744871391589*fskin[17]+0.7071067811865475*fskin[10])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[7] = (1.224744871391589*fskin[20]+0.7071067811865475*fskin[12])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[8] = (1.224744871391589*fskin[23]+0.7071067811865475*fskin[13])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[9] = (1.224744871391589*fskin[28]+0.7071067811865475*fskin[14])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[10] = (1.5811388300841898*fskin[44]+1.224744871391589*fskin[31]+0.7071067811865475*fskin[18])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[11] = (1.224744871391589*fskin[33]+0.7071067811865475*fskin[22])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[12] = (1.224744871391589*fskin[34]+0.7071067811865475*fskin[24])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[13] = (1.224744871391589*fskin[36]+0.7071067811865475*fskin[26])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[14] = (1.224744871391589*fskin[39]+0.7071067811865475*fskin[27])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[15] = (1.224744871391589*fskin[41]+0.7071067811865475*fskin[29])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[16] = (1.224744871391589*fskin[42]+0.7071067811865475*fskin[30])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[17] = (1.224744871391589*fskin[45]+0.7071067811865475*fskin[38])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[18] = (1.224744871391589*fskin[46]+0.7071067811865475*fskin[40])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[19] = (1.224744871391589*fskin[47]+0.7071067811865475*fskin[43])*jacob_cx_skin_inv*jacob_vx_inv; 

  } else { 

  fUpwind[0] = (1.5811388300841895*fedge[11]-1.224744871391589*fedge[1]+0.7071067811865475*fedge[0])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[1] = (1.5811388300841898*fedge[19]-1.224744871391589*fedge[5]+0.7071067811865475*fedge[2])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[2] = (1.5811388300841898*fedge[21]-1.224744871391589*fedge[6]+0.7071067811865475*fedge[3])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[3] = (1.5811388300841898*fedge[25]-1.224744871391589*fedge[8]+0.7071067811865475*fedge[4])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[4] = (1.5811388300841895*fedge[32]-1.224744871391589*fedge[15]+0.7071067811865475*fedge[7])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[5] = (1.5811388300841895*fedge[35]-1.224744871391589*fedge[16]+0.7071067811865475*fedge[9])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[6] = (1.5811388300841895*fedge[37]-1.224744871391589*fedge[17]+0.7071067811865475*fedge[10])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[7] = (0.7071067811865475*fedge[12]-1.224744871391589*fedge[20])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[8] = (0.7071067811865475*fedge[13]-1.224744871391589*fedge[23])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[9] = (0.7071067811865475*fedge[14]-1.224744871391589*fedge[28])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[10] = (1.5811388300841898*fedge[44]-1.224744871391589*fedge[31]+0.7071067811865475*fedge[18])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[11] = (0.7071067811865475*fedge[22]-1.224744871391589*fedge[33])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[12] = (0.7071067811865475*fedge[24]-1.224744871391589*fedge[34])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[13] = (0.7071067811865475*fedge[26]-1.224744871391589*fedge[36])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[14] = (0.7071067811865475*fedge[27]-1.224744871391589*fedge[39])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[15] = (0.7071067811865475*fedge[29]-1.224744871391589*fedge[41])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[16] = (0.7071067811865475*fedge[30]-1.224744871391589*fedge[42])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[17] = (0.7071067811865475*fedge[38]-1.224744871391589*fedge[45])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[18] = (0.7071067811865475*fedge[40]-1.224744871391589*fedge[46])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[19] = (0.7071067811865475*fedge[43]-1.224744871391589*fedge[47])*jacob_cx_edge_inv*jacob_vx_inv; 

  } 
  Ghat[0] = 1.9364916731037085*fUpwind[2]*hamil[4]+0.8660254037844386*fUpwind[0]*hamil[1]; 
  Ghat[1] = 1.9364916731037085*fUpwind[4]*hamil[4]+0.8660254037844386*fUpwind[1]*hamil[1]; 
  Ghat[2] = 1.7320508075688772*hamil[4]*fUpwind[8]+1.9364916731037085*fUpwind[0]*hamil[4]+0.8660254037844386*hamil[1]*fUpwind[2]; 
  Ghat[3] = 1.9364916731037085*hamil[4]*fUpwind[6]+0.8660254037844386*hamil[1]*fUpwind[3]; 
  Ghat[4] = 1.7320508075688774*hamil[4]*fUpwind[12]+1.9364916731037085*fUpwind[1]*hamil[4]+0.8660254037844386*hamil[1]*fUpwind[4]; 
  Ghat[5] = 1.9364916731037085*hamil[4]*fUpwind[10]+0.8660254037844386*hamil[1]*fUpwind[5]; 
  Ghat[6] = 1.7320508075688774*hamil[4]*fUpwind[14]+0.8660254037844386*hamil[1]*fUpwind[6]+1.9364916731037085*fUpwind[3]*hamil[4]; 
  Ghat[7] = 1.9364916731037085*hamil[4]*fUpwind[11]+0.8660254037844386*hamil[1]*fUpwind[7]; 
  Ghat[8] = 0.8660254037844386*hamil[1]*fUpwind[8]+1.7320508075688772*fUpwind[2]*hamil[4]; 
  Ghat[9] = 1.9364916731037085*hamil[4]*fUpwind[16]+0.8660254037844386*hamil[1]*fUpwind[9]; 
  Ghat[10] = 1.7320508075688772*hamil[4]*fUpwind[18]+0.8660254037844386*hamil[1]*fUpwind[10]+1.9364916731037085*hamil[4]*fUpwind[5]; 
  Ghat[11] = 0.8660254037844386*hamil[1]*fUpwind[11]+1.9364916731037085*hamil[4]*fUpwind[7]; 
  Ghat[12] = 0.8660254037844386*hamil[1]*fUpwind[12]+1.7320508075688774*fUpwind[4]*hamil[4]; 
  Ghat[13] = 1.9364916731037085*hamil[4]*fUpwind[17]+0.8660254037844386*hamil[1]*fUpwind[13]; 
  Ghat[14] = 0.8660254037844386*hamil[1]*fUpwind[14]+1.7320508075688774*hamil[4]*fUpwind[6]; 
  Ghat[15] = 1.9364916731037085*hamil[4]*fUpwind[19]+0.8660254037844386*hamil[1]*fUpwind[15]; 
  Ghat[16] = 0.8660254037844386*hamil[1]*fUpwind[16]+1.9364916731037085*hamil[4]*fUpwind[9]; 
  Ghat[17] = 0.8660254037844386*hamil[1]*fUpwind[17]+1.9364916731037085*hamil[4]*fUpwind[13]; 
  Ghat[18] = 0.8660254037844386*hamil[1]*fUpwind[18]+1.7320508075688772*hamil[4]*fUpwind[10]; 
  Ghat[19] = 0.8660254037844386*hamil[1]*fUpwind[19]+1.9364916731037085*hamil[4]*fUpwind[15]; 

  out[0] += -(0.7071067811865475*Ghat[0]*dv10*dx10); 
  out[1] += -(1.224744871391589*Ghat[0]*dv10*dx10); 
  out[2] += -(0.7071067811865475*Ghat[1]*dv10*dx10); 
  out[3] += -(0.7071067811865475*Ghat[2]*dv10*dx10); 
  out[4] += -(0.7071067811865475*Ghat[3]*dv10*dx10); 
  out[5] += -(1.224744871391589*Ghat[1]*dv10*dx10); 
  out[6] += -(1.224744871391589*Ghat[2]*dv10*dx10); 
  out[7] += -(0.7071067811865475*Ghat[4]*dv10*dx10); 
  out[8] += -(1.224744871391589*Ghat[3]*dv10*dx10); 
  out[9] += -(0.7071067811865475*Ghat[5]*dv10*dx10); 
  out[10] += -(0.7071067811865475*Ghat[6]*dv10*dx10); 
  out[11] += -(1.5811388300841895*Ghat[0]*dv10*dx10); 
  out[12] += -(0.7071067811865475*Ghat[7]*dv10*dx10); 
  out[13] += -(0.7071067811865475*Ghat[8]*dv10*dx10); 
  out[14] += -(0.7071067811865475*Ghat[9]*dv10*dx10); 
  out[15] += -(1.224744871391589*Ghat[4]*dv10*dx10); 
  out[16] += -(1.224744871391589*Ghat[5]*dv10*dx10); 
  out[17] += -(1.224744871391589*Ghat[6]*dv10*dx10); 
  out[18] += -(0.7071067811865475*Ghat[10]*dv10*dx10); 
  out[19] += -(1.5811388300841898*Ghat[1]*dv10*dx10); 
  out[20] += -(1.224744871391589*Ghat[7]*dv10*dx10); 
  out[21] += -(1.5811388300841898*Ghat[2]*dv10*dx10); 
  out[22] += -(0.7071067811865475*Ghat[11]*dv10*dx10); 
  out[23] += -(1.224744871391589*Ghat[8]*dv10*dx10); 
  out[24] += -(0.7071067811865475*Ghat[12]*dv10*dx10); 
  out[25] += -(1.5811388300841898*Ghat[3]*dv10*dx10); 
  out[26] += -(0.7071067811865475*Ghat[13]*dv10*dx10); 
  out[27] += -(0.7071067811865475*Ghat[14]*dv10*dx10); 
  out[28] += -(1.224744871391589*Ghat[9]*dv10*dx10); 
  out[29] += -(0.7071067811865475*Ghat[15]*dv10*dx10); 
  out[30] += -(0.7071067811865475*Ghat[16]*dv10*dx10); 
  out[31] += -(1.224744871391589*Ghat[10]*dv10*dx10); 
  out[32] += -(1.5811388300841895*Ghat[4]*dv10*dx10); 
  out[33] += -(1.224744871391589*Ghat[11]*dv10*dx10); 
  out[34] += -(1.224744871391589*Ghat[12]*dv10*dx10); 
  out[35] += -(1.5811388300841895*Ghat[5]*dv10*dx10); 
  out[36] += -(1.224744871391589*Ghat[13]*dv10*dx10); 
  out[37] += -(1.5811388300841895*Ghat[6]*dv10*dx10); 
  out[38] += -(0.7071067811865475*Ghat[17]*dv10*dx10); 
  out[39] += -(1.224744871391589*Ghat[14]*dv10*dx10); 
  out[40] += -(0.7071067811865475*Ghat[18]*dv10*dx10); 
  out[41] += -(1.224744871391589*Ghat[15]*dv10*dx10); 
  out[42] += -(1.224744871391589*Ghat[16]*dv10*dx10); 
  out[43] += -(0.7071067811865475*Ghat[19]*dv10*dx10); 
  out[44] += -(1.5811388300841898*Ghat[10]*dv10*dx10); 
  out[45] += -(1.224744871391589*Ghat[17]*dv10*dx10); 
  out[46] += -(1.224744871391589*Ghat[18]*dv10*dx10); 
  out[47] += -(1.224744871391589*Ghat[19]*dv10*dx10); 

  } else { 

  if (wv>0) { 

  fUpwind[0] = (1.5811388300841895*fedge[11]+1.224744871391589*fedge[1]+0.7071067811865475*fedge[0])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[1] = (1.5811388300841898*fedge[19]+1.224744871391589*fedge[5]+0.7071067811865475*fedge[2])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[2] = (1.5811388300841898*fedge[21]+1.224744871391589*fedge[6]+0.7071067811865475*fedge[3])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[3] = (1.5811388300841898*fedge[25]+1.224744871391589*fedge[8]+0.7071067811865475*fedge[4])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[4] = (1.5811388300841895*fedge[32]+1.224744871391589*fedge[15]+0.7071067811865475*fedge[7])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[5] = (1.5811388300841895*fedge[35]+1.224744871391589*fedge[16]+0.7071067811865475*fedge[9])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[6] = (1.5811388300841895*fedge[37]+1.224744871391589*fedge[17]+0.7071067811865475*fedge[10])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[7] = (1.224744871391589*fedge[20]+0.7071067811865475*fedge[12])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[8] = (1.224744871391589*fedge[23]+0.7071067811865475*fedge[13])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[9] = (1.224744871391589*fedge[28]+0.7071067811865475*fedge[14])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[10] = (1.5811388300841898*fedge[44]+1.224744871391589*fedge[31]+0.7071067811865475*fedge[18])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[11] = (1.224744871391589*fedge[33]+0.7071067811865475*fedge[22])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[12] = (1.224744871391589*fedge[34]+0.7071067811865475*fedge[24])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[13] = (1.224744871391589*fedge[36]+0.7071067811865475*fedge[26])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[14] = (1.224744871391589*fedge[39]+0.7071067811865475*fedge[27])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[15] = (1.224744871391589*fedge[41]+0.7071067811865475*fedge[29])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[16] = (1.224744871391589*fedge[42]+0.7071067811865475*fedge[30])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[17] = (1.224744871391589*fedge[45]+0.7071067811865475*fedge[38])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[18] = (1.224744871391589*fedge[46]+0.7071067811865475*fedge[40])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[19] = (1.224744871391589*fedge[47]+0.7071067811865475*fedge[43])*jacob_cx_edge_inv*jacob_vx_inv; 

  } else { 

  fUpwind[0] = (1.5811388300841895*fskin[11]-1.224744871391589*fskin[1]+0.7071067811865475*fskin[0])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[1] = (1.5811388300841898*fskin[19]-1.224744871391589*fskin[5]+0.7071067811865475*fskin[2])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[2] = (1.5811388300841898*fskin[21]-1.224744871391589*fskin[6]+0.7071067811865475*fskin[3])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[3] = (1.5811388300841898*fskin[25]-1.224744871391589*fskin[8]+0.7071067811865475*fskin[4])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[4] = (1.5811388300841895*fskin[32]-1.224744871391589*fskin[15]+0.7071067811865475*fskin[7])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[5] = (1.5811388300841895*fskin[35]-1.224744871391589*fskin[16]+0.7071067811865475*fskin[9])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[6] = (1.5811388300841895*fskin[37]-1.224744871391589*fskin[17]+0.7071067811865475*fskin[10])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[7] = (0.7071067811865475*fskin[12]-1.224744871391589*fskin[20])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[8] = (0.7071067811865475*fskin[13]-1.224744871391589*fskin[23])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[9] = (0.7071067811865475*fskin[14]-1.224744871391589*fskin[28])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[10] = (1.5811388300841898*fskin[44]-1.224744871391589*fskin[31]+0.7071067811865475*fskin[18])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[11] = (0.7071067811865475*fskin[22]-1.224744871391589*fskin[33])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[12] = (0.7071067811865475*fskin[24]-1.224744871391589*fskin[34])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[13] = (0.7071067811865475*fskin[26]-1.224744871391589*fskin[36])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[14] = (0.7071067811865475*fskin[27]-1.224744871391589*fskin[39])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[15] = (0.7071067811865475*fskin[29]-1.224744871391589*fskin[41])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[16] = (0.7071067811865475*fskin[30]-1.224744871391589*fskin[42])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[17] = (0.7071067811865475*fskin[38]-1.224744871391589*fskin[45])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[18] = (0.7071067811865475*fskin[40]-1.224744871391589*fskin[46])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[19] = (0.7071067811865475*fskin[43]-1.224744871391589*fskin[47])*jacob_cx_skin_inv*jacob_vx_inv; 

  } 
  Ghat[0] = 1.9364916731037085*fUpwind[2]*hamil[4]+0.8660254037844386*fUpwind[0]*hamil[1]; 
  Ghat[1] = 1.9364916731037085*fUpwind[4]*hamil[4]+0.8660254037844386*fUpwind[1]*hamil[1]; 
  Ghat[2] = 1.7320508075688772*hamil[4]*fUpwind[8]+1.9364916731037085*fUpwind[0]*hamil[4]+0.8660254037844386*hamil[1]*fUpwind[2]; 
  Ghat[3] = 1.9364916731037085*hamil[4]*fUpwind[6]+0.8660254037844386*hamil[1]*fUpwind[3]; 
  Ghat[4] = 1.7320508075688774*hamil[4]*fUpwind[12]+1.9364916731037085*fUpwind[1]*hamil[4]+0.8660254037844386*hamil[1]*fUpwind[4]; 
  Ghat[5] = 1.9364916731037085*hamil[4]*fUpwind[10]+0.8660254037844386*hamil[1]*fUpwind[5]; 
  Ghat[6] = 1.7320508075688774*hamil[4]*fUpwind[14]+0.8660254037844386*hamil[1]*fUpwind[6]+1.9364916731037085*fUpwind[3]*hamil[4]; 
  Ghat[7] = 1.9364916731037085*hamil[4]*fUpwind[11]+0.8660254037844386*hamil[1]*fUpwind[7]; 
  Ghat[8] = 0.8660254037844386*hamil[1]*fUpwind[8]+1.7320508075688772*fUpwind[2]*hamil[4]; 
  Ghat[9] = 1.9364916731037085*hamil[4]*fUpwind[16]+0.8660254037844386*hamil[1]*fUpwind[9]; 
  Ghat[10] = 1.7320508075688772*hamil[4]*fUpwind[18]+0.8660254037844386*hamil[1]*fUpwind[10]+1.9364916731037085*hamil[4]*fUpwind[5]; 
  Ghat[11] = 0.8660254037844386*hamil[1]*fUpwind[11]+1.9364916731037085*hamil[4]*fUpwind[7]; 
  Ghat[12] = 0.8660254037844386*hamil[1]*fUpwind[12]+1.7320508075688774*fUpwind[4]*hamil[4]; 
  Ghat[13] = 1.9364916731037085*hamil[4]*fUpwind[17]+0.8660254037844386*hamil[1]*fUpwind[13]; 
  Ghat[14] = 0.8660254037844386*hamil[1]*fUpwind[14]+1.7320508075688774*hamil[4]*fUpwind[6]; 
  Ghat[15] = 1.9364916731037085*hamil[4]*fUpwind[19]+0.8660254037844386*hamil[1]*fUpwind[15]; 
  Ghat[16] = 0.8660254037844386*hamil[1]*fUpwind[16]+1.9364916731037085*hamil[4]*fUpwind[9]; 
  Ghat[17] = 0.8660254037844386*hamil[1]*fUpwind[17]+1.9364916731037085*hamil[4]*fUpwind[13]; 
  Ghat[18] = 0.8660254037844386*hamil[1]*fUpwind[18]+1.7320508075688772*hamil[4]*fUpwind[10]; 
  Ghat[19] = 0.8660254037844386*hamil[1]*fUpwind[19]+1.9364916731037085*hamil[4]*fUpwind[15]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv10*dx10; 
  out[1] += -(1.224744871391589*Ghat[0]*dv10*dx10); 
  out[2] += 0.7071067811865475*Ghat[1]*dv10*dx10; 
  out[3] += 0.7071067811865475*Ghat[2]*dv10*dx10; 
  out[4] += 0.7071067811865475*Ghat[3]*dv10*dx10; 
  out[5] += -(1.224744871391589*Ghat[1]*dv10*dx10); 
  out[6] += -(1.224744871391589*Ghat[2]*dv10*dx10); 
  out[7] += 0.7071067811865475*Ghat[4]*dv10*dx10; 
  out[8] += -(1.224744871391589*Ghat[3]*dv10*dx10); 
  out[9] += 0.7071067811865475*Ghat[5]*dv10*dx10; 
  out[10] += 0.7071067811865475*Ghat[6]*dv10*dx10; 
  out[11] += 1.5811388300841895*Ghat[0]*dv10*dx10; 
  out[12] += 0.7071067811865475*Ghat[7]*dv10*dx10; 
  out[13] += 0.7071067811865475*Ghat[8]*dv10*dx10; 
  out[14] += 0.7071067811865475*Ghat[9]*dv10*dx10; 
  out[15] += -(1.224744871391589*Ghat[4]*dv10*dx10); 
  out[16] += -(1.224744871391589*Ghat[5]*dv10*dx10); 
  out[17] += -(1.224744871391589*Ghat[6]*dv10*dx10); 
  out[18] += 0.7071067811865475*Ghat[10]*dv10*dx10; 
  out[19] += 1.5811388300841898*Ghat[1]*dv10*dx10; 
  out[20] += -(1.224744871391589*Ghat[7]*dv10*dx10); 
  out[21] += 1.5811388300841898*Ghat[2]*dv10*dx10; 
  out[22] += 0.7071067811865475*Ghat[11]*dv10*dx10; 
  out[23] += -(1.224744871391589*Ghat[8]*dv10*dx10); 
  out[24] += 0.7071067811865475*Ghat[12]*dv10*dx10; 
  out[25] += 1.5811388300841898*Ghat[3]*dv10*dx10; 
  out[26] += 0.7071067811865475*Ghat[13]*dv10*dx10; 
  out[27] += 0.7071067811865475*Ghat[14]*dv10*dx10; 
  out[28] += -(1.224744871391589*Ghat[9]*dv10*dx10); 
  out[29] += 0.7071067811865475*Ghat[15]*dv10*dx10; 
  out[30] += 0.7071067811865475*Ghat[16]*dv10*dx10; 
  out[31] += -(1.224744871391589*Ghat[10]*dv10*dx10); 
  out[32] += 1.5811388300841895*Ghat[4]*dv10*dx10; 
  out[33] += -(1.224744871391589*Ghat[11]*dv10*dx10); 
  out[34] += -(1.224744871391589*Ghat[12]*dv10*dx10); 
  out[35] += 1.5811388300841895*Ghat[5]*dv10*dx10; 
  out[36] += -(1.224744871391589*Ghat[13]*dv10*dx10); 
  out[37] += 1.5811388300841895*Ghat[6]*dv10*dx10; 
  out[38] += 0.7071067811865475*Ghat[17]*dv10*dx10; 
  out[39] += -(1.224744871391589*Ghat[14]*dv10*dx10); 
  out[40] += 0.7071067811865475*Ghat[18]*dv10*dx10; 
  out[41] += -(1.224744871391589*Ghat[15]*dv10*dx10); 
  out[42] += -(1.224744871391589*Ghat[16]*dv10*dx10); 
  out[43] += 0.7071067811865475*Ghat[19]*dv10*dx10; 
  out[44] += 1.5811388300841898*Ghat[10]*dv10*dx10; 
  out[45] += -(1.224744871391589*Ghat[17]*dv10*dx10); 
  out[46] += -(1.224744871391589*Ghat[18]*dv10*dx10); 
  out[47] += -(1.224744871391589*Ghat[19]*dv10*dx10); 

  } 
  return 0.0;

} 
