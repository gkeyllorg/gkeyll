#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfy_2x2v_ser_p1(const double *w, const double *dxv,
  const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel,
  const double *poisson_tensor_conf, const double *hamil,
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  double dx11 = 2.0/dxv[1]; 
  double dv11 = 2.0/dxv[3]; 
  double wv = w[3]; 

  double fUpwind[8]; 
  double Ghat[8]; 

  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 

  const double *jacob_cy_edge = &jacob_pos_edge[2]; 
  const double *jacob_cy_skin = &jacob_pos_skin[2]; 
  const double jacob_cy_edge_inv = 1.0/jacob_cy_edge[0]; 
  const double jacob_cy_skin_inv = 1.0/jacob_cy_skin[0]; 

  if (edge == -1) { 

  if (wv>0) { 

  fUpwind[0] = (1.224744871391589*fskin[2]+0.7071067811865475*fskin[0])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[1] = (1.224744871391589*fskin[5]+0.7071067811865475*fskin[1])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[2] = (1.224744871391589*fskin[7]+0.7071067811865475*fskin[3])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[3] = (1.224744871391589*fskin[9]+0.7071067811865475*fskin[4])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[4] = (1.224744871391589*fskin[11]+0.7071067811865475*fskin[6])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[5] = (1.224744871391589*fskin[12]+0.7071067811865475*fskin[8])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[6] = (1.224744871391589*fskin[14]+0.7071067811865475*fskin[10])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[7] = (1.224744871391589*fskin[15]+0.7071067811865475*fskin[13])*jacob_cy_skin_inv*jacob_vy_inv; 

  } else { 

  fUpwind[0] = (0.7071067811865475*fedge[0]-1.224744871391589*fedge[2])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[1] = (0.7071067811865475*fedge[1]-1.224744871391589*fedge[5])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[2] = (0.7071067811865475*fedge[3]-1.224744871391589*fedge[7])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[3] = (0.7071067811865475*fedge[4]-1.224744871391589*fedge[9])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[4] = (0.7071067811865475*fedge[6]-1.224744871391589*fedge[11])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[5] = (0.7071067811865475*fedge[8]-1.224744871391589*fedge[12])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[6] = (0.7071067811865475*fedge[10]-1.224744871391589*fedge[14])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[7] = (0.7071067811865475*fedge[13]-1.224744871391589*fedge[15])*jacob_cy_edge_inv*jacob_vy_inv; 

  } 
  Ghat[0] = 0.8660254037844386*fUpwind[2]*hamil[3]+0.8660254037844386*fUpwind[0]*hamil[2]; 
  Ghat[1] = 0.8660254037844386*hamil[3]*fUpwind[4]+0.8660254037844386*fUpwind[1]*hamil[2]; 
  Ghat[2] = 0.8660254037844386*fUpwind[0]*hamil[3]+0.8660254037844386*fUpwind[2]*hamil[2]; 
  Ghat[3] = 0.8660254037844386*hamil[3]*fUpwind[6]+0.8660254037844386*hamil[2]*fUpwind[3]; 
  Ghat[4] = 0.8660254037844386*hamil[2]*fUpwind[4]+0.8660254037844386*fUpwind[1]*hamil[3]; 
  Ghat[5] = 0.8660254037844386*hamil[3]*fUpwind[7]+0.8660254037844386*hamil[2]*fUpwind[5]; 
  Ghat[6] = 0.8660254037844386*hamil[2]*fUpwind[6]+0.8660254037844386*fUpwind[3]*hamil[3]; 
  Ghat[7] = 0.8660254037844386*hamil[2]*fUpwind[7]+0.8660254037844386*hamil[3]*fUpwind[5]; 

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
  out[11] += -(1.224744871391589*Ghat[4]*dv11*dx11); 
  out[12] += -(1.224744871391589*Ghat[5]*dv11*dx11); 
  out[13] += -(0.7071067811865475*Ghat[7]*dv11*dx11); 
  out[14] += -(1.224744871391589*Ghat[6]*dv11*dx11); 
  out[15] += -(1.224744871391589*Ghat[7]*dv11*dx11); 

  } else { 

  if (wv>0) { 

  fUpwind[0] = (1.224744871391589*fedge[2]+0.7071067811865475*fedge[0])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[1] = (1.224744871391589*fedge[5]+0.7071067811865475*fedge[1])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[2] = (1.224744871391589*fedge[7]+0.7071067811865475*fedge[3])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[3] = (1.224744871391589*fedge[9]+0.7071067811865475*fedge[4])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[4] = (1.224744871391589*fedge[11]+0.7071067811865475*fedge[6])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[5] = (1.224744871391589*fedge[12]+0.7071067811865475*fedge[8])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[6] = (1.224744871391589*fedge[14]+0.7071067811865475*fedge[10])*jacob_cy_edge_inv*jacob_vy_inv; 
  fUpwind[7] = (1.224744871391589*fedge[15]+0.7071067811865475*fedge[13])*jacob_cy_edge_inv*jacob_vy_inv; 

  } else { 

  fUpwind[0] = (0.7071067811865475*fskin[0]-1.224744871391589*fskin[2])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[1] = (0.7071067811865475*fskin[1]-1.224744871391589*fskin[5])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[2] = (0.7071067811865475*fskin[3]-1.224744871391589*fskin[7])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[3] = (0.7071067811865475*fskin[4]-1.224744871391589*fskin[9])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[4] = (0.7071067811865475*fskin[6]-1.224744871391589*fskin[11])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[5] = (0.7071067811865475*fskin[8]-1.224744871391589*fskin[12])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[6] = (0.7071067811865475*fskin[10]-1.224744871391589*fskin[14])*jacob_cy_skin_inv*jacob_vy_inv; 
  fUpwind[7] = (0.7071067811865475*fskin[13]-1.224744871391589*fskin[15])*jacob_cy_skin_inv*jacob_vy_inv; 

  } 
  Ghat[0] = 0.8660254037844386*fUpwind[2]*hamil[3]+0.8660254037844386*fUpwind[0]*hamil[2]; 
  Ghat[1] = 0.8660254037844386*hamil[3]*fUpwind[4]+0.8660254037844386*fUpwind[1]*hamil[2]; 
  Ghat[2] = 0.8660254037844386*fUpwind[0]*hamil[3]+0.8660254037844386*fUpwind[2]*hamil[2]; 
  Ghat[3] = 0.8660254037844386*hamil[3]*fUpwind[6]+0.8660254037844386*hamil[2]*fUpwind[3]; 
  Ghat[4] = 0.8660254037844386*hamil[2]*fUpwind[4]+0.8660254037844386*fUpwind[1]*hamil[3]; 
  Ghat[5] = 0.8660254037844386*hamil[3]*fUpwind[7]+0.8660254037844386*hamil[2]*fUpwind[5]; 
  Ghat[6] = 0.8660254037844386*hamil[2]*fUpwind[6]+0.8660254037844386*fUpwind[3]*hamil[3]; 
  Ghat[7] = 0.8660254037844386*hamil[2]*fUpwind[7]+0.8660254037844386*hamil[3]*fUpwind[5]; 

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
  out[11] += -(1.224744871391589*Ghat[4]*dv11*dx11); 
  out[12] += -(1.224744871391589*Ghat[5]*dv11*dx11); 
  out[13] += 0.7071067811865475*Ghat[7]*dv11*dx11; 
  out[14] += -(1.224744871391589*Ghat[6]*dv11*dx11); 
  out[15] += -(1.224744871391589*Ghat[7]*dv11*dx11); 

  } 
  return 0.0;

} 
