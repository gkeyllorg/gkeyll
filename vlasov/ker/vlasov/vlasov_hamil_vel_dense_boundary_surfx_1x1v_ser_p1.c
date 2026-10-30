#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfx_1x1v_ser_p1(const double *w, const double *dxv,
  const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel,
  const double *poisson_tensor_conf, const double *hamil,
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[1]; 
  double wv = w[1]; 

  double fUpwind[2]; 
  double Ghat[2]; 

  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 

  const double *jacob_cx_edge = &jacob_pos_edge[0]; 
  const double *jacob_cx_skin = &jacob_pos_skin[0]; 
  const double jacob_cx_edge_inv = 1.0/jacob_cx_edge[0]; 
  const double jacob_cx_skin_inv = 1.0/jacob_cx_skin[0]; 

  if (edge == -1) { 

  if (wv>0) { 

  fUpwind[0] = (1.224744871391589*fskin[1]+0.7071067811865475*fskin[0])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[1] = (1.224744871391589*fskin[3]+0.7071067811865475*fskin[2])*jacob_cx_skin_inv*jacob_vx_inv; 

  } else { 

  fUpwind[0] = (0.7071067811865475*fedge[0]-1.224744871391589*fedge[1])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[1] = (0.7071067811865475*fedge[2]-1.224744871391589*fedge[3])*jacob_cx_edge_inv*jacob_vx_inv; 

  } 
  Ghat[0] = 1.224744871391589*fUpwind[0]*hamil[1]; 
  Ghat[1] = 1.224744871391589*fUpwind[1]*hamil[1]; 

  out[0] += -(0.7071067811865475*Ghat[0]*dv10*dx10); 
  out[1] += -(1.224744871391589*Ghat[0]*dv10*dx10); 
  out[2] += -(0.7071067811865475*Ghat[1]*dv10*dx10); 
  out[3] += -(1.224744871391589*Ghat[1]*dv10*dx10); 

  } else { 

  if (wv>0) { 

  fUpwind[0] = (1.224744871391589*fedge[1]+0.7071067811865475*fedge[0])*jacob_cx_edge_inv*jacob_vx_inv; 
  fUpwind[1] = (1.224744871391589*fedge[3]+0.7071067811865475*fedge[2])*jacob_cx_edge_inv*jacob_vx_inv; 

  } else { 

  fUpwind[0] = (0.7071067811865475*fskin[0]-1.224744871391589*fskin[1])*jacob_cx_skin_inv*jacob_vx_inv; 
  fUpwind[1] = (0.7071067811865475*fskin[2]-1.224744871391589*fskin[3])*jacob_cx_skin_inv*jacob_vx_inv; 

  } 
  Ghat[0] = 1.224744871391589*fUpwind[0]*hamil[1]; 
  Ghat[1] = 1.224744871391589*fUpwind[1]*hamil[1]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv10*dx10; 
  out[1] += -(1.224744871391589*Ghat[0]*dv10*dx10); 
  out[2] += 0.7071067811865475*Ghat[1]*dv10*dx10; 
  out[3] += -(1.224744871391589*Ghat[1]*dv10*dx10); 

  } 
  return 0.0;

} 
