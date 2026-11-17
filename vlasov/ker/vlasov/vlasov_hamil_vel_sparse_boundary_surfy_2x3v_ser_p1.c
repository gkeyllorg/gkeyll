#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_vel_sparse_boundary_surfy_2x3v_ser_p1(const double *w, const double *dxv,
  const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel,
  const double *poisson_tensor_conf, const double *hamil,
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  double dx11 = 2.0/dxv[1]; 
  double dv11 = 2.0/dxv[3]; 
  double wv = w[3]; 

  double fUpwind[16]; 
  double Ghat[16]; 

  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 

  const double *jacob_cy_edge = &jacob_pos_edge[2]; 
  const double *jacob_cy_skin = &jacob_pos_skin[2]; 
  const double jacob_cy_edge_inv = 1.0/jacob_cy_edge[0]; 
  const double jacob_cy_skin_inv = 1.0/jacob_cy_skin[0]; 

  if (edge == -1) { 

  if (wv>0) { 

  fUpwind[0] = (1.224744871391589*fskin[2]+0.7071067811865475*fskin[0])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[1] = (1.224744871391589*fskin[6]+0.7071067811865475*fskin[1])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[2] = (1.224744871391589*fskin[8]+0.7071067811865475*fskin[3])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[3] = (1.224744871391589*fskin[10]+0.7071067811865475*fskin[4])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[4] = (1.224744871391589*fskin[13]+0.7071067811865475*fskin[5])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[5] = (1.224744871391589*fskin[16]+0.7071067811865475*fskin[7])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[6] = (1.224744871391589*fskin[17]+0.7071067811865475*fskin[9])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[7] = (1.224744871391589*fskin[19]+0.7071067811865475*fskin[11])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[8] = (1.224744871391589*fskin[20]+0.7071067811865475*fskin[12])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[9] = (1.224744871391589*fskin[22]+0.7071067811865475*fskin[14])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[10] = (1.224744871391589*fskin[24]+0.7071067811865475*fskin[15])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[11] = (1.224744871391589*fskin[26]+0.7071067811865475*fskin[18])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[12] = (1.224744871391589*fskin[27]+0.7071067811865475*fskin[21])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[13] = (1.224744871391589*fskin[28]+0.7071067811865475*fskin[23])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[14] = (1.224744871391589*fskin[30]+0.7071067811865475*fskin[25])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[15] = (1.224744871391589*fskin[31]+0.7071067811865475*fskin[29])*jacob_cy_skin_inv*jacob_vy_inv; 

  } else { 

  fUpwind[0] = (0.7071067811865475*fedge[0]-1.224744871391589*fedge[2])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[1] = (0.7071067811865475*fedge[1]-1.224744871391589*fedge[6])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[2] = (0.7071067811865475*fedge[3]-1.224744871391589*fedge[8])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[3] = (0.7071067811865475*fedge[4]-1.224744871391589*fedge[10])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[4] = (0.7071067811865475*fedge[5]-1.224744871391589*fedge[13])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[5] = (0.7071067811865475*fedge[7]-1.224744871391589*fedge[16])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[6] = (0.7071067811865475*fedge[9]-1.224744871391589*fedge[17])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[7] = (0.7071067811865475*fedge[11]-1.224744871391589*fedge[19])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[8] = (0.7071067811865475*fedge[12]-1.224744871391589*fedge[20])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[9] = (0.7071067811865475*fedge[14]-1.224744871391589*fedge[22])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[10] = (0.7071067811865475*fedge[15]-1.224744871391589*fedge[24])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[11] = (0.7071067811865475*fedge[18]-1.224744871391589*fedge[26])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[12] = (0.7071067811865475*fedge[21]-1.224744871391589*fedge[27])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[13] = (0.7071067811865475*fedge[23]-1.224744871391589*fedge[28])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[14] = (0.7071067811865475*fedge[25]-1.224744871391589*fedge[30])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[15] = (0.7071067811865475*fedge[29]-1.224744871391589*fedge[31])*jacob_cy_edge_inv*jacob_vy_inv; 

  } 
  Ghat[0] = 0.6123724356957944*fUpwind[0]*hamil[2]; 
  Ghat[1] = 0.6123724356957944*fUpwind[1]*hamil[2]; 
  Ghat[2] = 0.6123724356957944*fUpwind[2]*hamil[2]; 
  Ghat[3] = 0.6123724356957944*hamil[2]*fUpwind[3]; 
  Ghat[4] = 0.6123724356957944*hamil[2]*fUpwind[4]; 
  Ghat[5] = 0.6123724356957944*hamil[2]*fUpwind[5]; 
  Ghat[6] = 0.6123724356957944*hamil[2]*fUpwind[6]; 
  Ghat[7] = 0.6123724356957944*hamil[2]*fUpwind[7]; 
  Ghat[8] = 0.6123724356957944*hamil[2]*fUpwind[8]; 
  Ghat[9] = 0.6123724356957944*hamil[2]*fUpwind[9]; 
  Ghat[10] = 0.6123724356957944*hamil[2]*fUpwind[10]; 
  Ghat[11] = 0.6123724356957944*hamil[2]*fUpwind[11]; 
  Ghat[12] = 0.6123724356957944*hamil[2]*fUpwind[12]; 
  Ghat[13] = 0.6123724356957944*hamil[2]*fUpwind[13]; 
  Ghat[14] = 0.6123724356957944*hamil[2]*fUpwind[14]; 
  Ghat[15] = 0.6123724356957944*hamil[2]*fUpwind[15]; 

  out[0] += -(0.7071067811865475*Ghat[0]*dv11*dx11); 
  out[1] += -(0.7071067811865475*Ghat[1]*dv11*dx11); 
  out[2] += -(1.224744871391589*Ghat[0]*dv11*dx11); 
  out[3] += -(0.7071067811865475*Ghat[2]*dv11*dx11); 
  out[4] += -(0.7071067811865475*Ghat[3]*dv11*dx11); 
  out[5] += -(0.7071067811865475*Ghat[4]*dv11*dx11); 
  out[6] += -(1.224744871391589*Ghat[1]*dv11*dx11); 
  out[7] += -(0.7071067811865475*Ghat[5]*dv11*dx11); 
  out[8] += -(1.224744871391589*Ghat[2]*dv11*dx11); 
  out[9] += -(0.7071067811865475*Ghat[6]*dv11*dx11); 
  out[10] += -(1.224744871391589*Ghat[3]*dv11*dx11); 
  out[11] += -(0.7071067811865475*Ghat[7]*dv11*dx11); 
  out[12] += -(0.7071067811865475*Ghat[8]*dv11*dx11); 
  out[13] += -(1.224744871391589*Ghat[4]*dv11*dx11); 
  out[14] += -(0.7071067811865475*Ghat[9]*dv11*dx11); 
  out[15] += -(0.7071067811865475*Ghat[10]*dv11*dx11); 
  out[16] += -(1.224744871391589*Ghat[5]*dv11*dx11); 
  out[17] += -(1.224744871391589*Ghat[6]*dv11*dx11); 
  out[18] += -(0.7071067811865475*Ghat[11]*dv11*dx11); 
  out[19] += -(1.224744871391589*Ghat[7]*dv11*dx11); 
  out[20] += -(1.224744871391589*Ghat[8]*dv11*dx11); 
  out[21] += -(0.7071067811865475*Ghat[12]*dv11*dx11); 
  out[22] += -(1.224744871391589*Ghat[9]*dv11*dx11); 
  out[23] += -(0.7071067811865475*Ghat[13]*dv11*dx11); 
  out[24] += -(1.224744871391589*Ghat[10]*dv11*dx11); 
  out[25] += -(0.7071067811865475*Ghat[14]*dv11*dx11); 
  out[26] += -(1.224744871391589*Ghat[11]*dv11*dx11); 
  out[27] += -(1.224744871391589*Ghat[12]*dv11*dx11); 
  out[28] += -(1.224744871391589*Ghat[13]*dv11*dx11); 
  out[29] += -(0.7071067811865475*Ghat[15]*dv11*dx11); 
  out[30] += -(1.224744871391589*Ghat[14]*dv11*dx11); 
  out[31] += -(1.224744871391589*Ghat[15]*dv11*dx11); 

  } else { 

  if (wv>0) { 

  fUpwind[0] = (1.224744871391589*fedge[2]+0.7071067811865475*fedge[0])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[1] = (1.224744871391589*fedge[6]+0.7071067811865475*fedge[1])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[2] = (1.224744871391589*fedge[8]+0.7071067811865475*fedge[3])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[3] = (1.224744871391589*fedge[10]+0.7071067811865475*fedge[4])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[4] = (1.224744871391589*fedge[13]+0.7071067811865475*fedge[5])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[5] = (1.224744871391589*fedge[16]+0.7071067811865475*fedge[7])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[6] = (1.224744871391589*fedge[17]+0.7071067811865475*fedge[9])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[7] = (1.224744871391589*fedge[19]+0.7071067811865475*fedge[11])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[8] = (1.224744871391589*fedge[20]+0.7071067811865475*fedge[12])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[9] = (1.224744871391589*fedge[22]+0.7071067811865475*fedge[14])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[10] = (1.224744871391589*fedge[24]+0.7071067811865475*fedge[15])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[11] = (1.224744871391589*fedge[26]+0.7071067811865475*fedge[18])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[12] = (1.224744871391589*fedge[27]+0.7071067811865475*fedge[21])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[13] = (1.224744871391589*fedge[28]+0.7071067811865475*fedge[23])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[14] = (1.224744871391589*fedge[30]+0.7071067811865475*fedge[25])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[15] = (1.224744871391589*fedge[31]+0.7071067811865475*fedge[29])*jacob_cy_edge_inv*jacob_vy_inv; 

  } else { 

  fUpwind[0] = (0.7071067811865475*fskin[0]-1.224744871391589*fskin[2])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[1] = (0.7071067811865475*fskin[1]-1.224744871391589*fskin[6])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[2] = (0.7071067811865475*fskin[3]-1.224744871391589*fskin[8])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[3] = (0.7071067811865475*fskin[4]-1.224744871391589*fskin[10])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[4] = (0.7071067811865475*fskin[5]-1.224744871391589*fskin[13])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[5] = (0.7071067811865475*fskin[7]-1.224744871391589*fskin[16])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[6] = (0.7071067811865475*fskin[9]-1.224744871391589*fskin[17])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[7] = (0.7071067811865475*fskin[11]-1.224744871391589*fskin[19])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[8] = (0.7071067811865475*fskin[12]-1.224744871391589*fskin[20])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[9] = (0.7071067811865475*fskin[14]-1.224744871391589*fskin[22])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[10] = (0.7071067811865475*fskin[15]-1.224744871391589*fskin[24])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[11] = (0.7071067811865475*fskin[18]-1.224744871391589*fskin[26])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[12] = (0.7071067811865475*fskin[21]-1.224744871391589*fskin[27])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[13] = (0.7071067811865475*fskin[23]-1.224744871391589*fskin[28])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[14] = (0.7071067811865475*fskin[25]-1.224744871391589*fskin[30])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[15] = (0.7071067811865475*fskin[29]-1.224744871391589*fskin[31])*jacob_cy_skin_inv*jacob_vy_inv; 

  } 
  Ghat[0] = 0.6123724356957944*fUpwind[0]*hamil[2]; 
  Ghat[1] = 0.6123724356957944*fUpwind[1]*hamil[2]; 
  Ghat[2] = 0.6123724356957944*fUpwind[2]*hamil[2]; 
  Ghat[3] = 0.6123724356957944*hamil[2]*fUpwind[3]; 
  Ghat[4] = 0.6123724356957944*hamil[2]*fUpwind[4]; 
  Ghat[5] = 0.6123724356957944*hamil[2]*fUpwind[5]; 
  Ghat[6] = 0.6123724356957944*hamil[2]*fUpwind[6]; 
  Ghat[7] = 0.6123724356957944*hamil[2]*fUpwind[7]; 
  Ghat[8] = 0.6123724356957944*hamil[2]*fUpwind[8]; 
  Ghat[9] = 0.6123724356957944*hamil[2]*fUpwind[9]; 
  Ghat[10] = 0.6123724356957944*hamil[2]*fUpwind[10]; 
  Ghat[11] = 0.6123724356957944*hamil[2]*fUpwind[11]; 
  Ghat[12] = 0.6123724356957944*hamil[2]*fUpwind[12]; 
  Ghat[13] = 0.6123724356957944*hamil[2]*fUpwind[13]; 
  Ghat[14] = 0.6123724356957944*hamil[2]*fUpwind[14]; 
  Ghat[15] = 0.6123724356957944*hamil[2]*fUpwind[15]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv11*dx11; 
  out[1] += 0.7071067811865475*Ghat[1]*dv11*dx11; 
  out[2] += -(1.224744871391589*Ghat[0]*dv11*dx11); 
  out[3] += 0.7071067811865475*Ghat[2]*dv11*dx11; 
  out[4] += 0.7071067811865475*Ghat[3]*dv11*dx11; 
  out[5] += 0.7071067811865475*Ghat[4]*dv11*dx11; 
  out[6] += -(1.224744871391589*Ghat[1]*dv11*dx11); 
  out[7] += 0.7071067811865475*Ghat[5]*dv11*dx11; 
  out[8] += -(1.224744871391589*Ghat[2]*dv11*dx11); 
  out[9] += 0.7071067811865475*Ghat[6]*dv11*dx11; 
  out[10] += -(1.224744871391589*Ghat[3]*dv11*dx11); 
  out[11] += 0.7071067811865475*Ghat[7]*dv11*dx11; 
  out[12] += 0.7071067811865475*Ghat[8]*dv11*dx11; 
  out[13] += -(1.224744871391589*Ghat[4]*dv11*dx11); 
  out[14] += 0.7071067811865475*Ghat[9]*dv11*dx11; 
  out[15] += 0.7071067811865475*Ghat[10]*dv11*dx11; 
  out[16] += -(1.224744871391589*Ghat[5]*dv11*dx11); 
  out[17] += -(1.224744871391589*Ghat[6]*dv11*dx11); 
  out[18] += 0.7071067811865475*Ghat[11]*dv11*dx11; 
  out[19] += -(1.224744871391589*Ghat[7]*dv11*dx11); 
  out[20] += -(1.224744871391589*Ghat[8]*dv11*dx11); 
  out[21] += 0.7071067811865475*Ghat[12]*dv11*dx11; 
  out[22] += -(1.224744871391589*Ghat[9]*dv11*dx11); 
  out[23] += 0.7071067811865475*Ghat[13]*dv11*dx11; 
  out[24] += -(1.224744871391589*Ghat[10]*dv11*dx11); 
  out[25] += 0.7071067811865475*Ghat[14]*dv11*dx11; 
  out[26] += -(1.224744871391589*Ghat[11]*dv11*dx11); 
  out[27] += -(1.224744871391589*Ghat[12]*dv11*dx11); 
  out[28] += -(1.224744871391589*Ghat[13]*dv11*dx11); 
  out[29] += 0.7071067811865475*Ghat[15]*dv11*dx11; 
  out[30] += -(1.224744871391589*Ghat[14]*dv11*dx11); 
  out[31] += -(1.224744871391589*Ghat[15]*dv11*dx11); 

  } 
  return 0.0;

} 
