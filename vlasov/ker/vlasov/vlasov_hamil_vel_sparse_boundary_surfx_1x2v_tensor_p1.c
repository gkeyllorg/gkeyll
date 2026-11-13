#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_vel_sparse_boundary_surfx_1x2v_tensor_p1(const double *w, const double *dxv,
  const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel,
  const double *poisson_tensor_conf, const double *hamil,
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  const double dx10 = 2.0/dxv[0]; 
  const double dv10 = 2.0/dxv[1]; 
  const double wv = w[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  double bdiv_vx[9]; 
  bdiv_vx[0] = 0.2777777777777778/jacob_vx[2]+0.4444444444444444/jacob_vx[1]+0.2777777777777778/jacob_vx[0]; 
  bdiv_vx[1] = 0.37267799624996495/jacob_vx[2]-0.37267799624996495/jacob_vx[0]; 
  bdiv_vx[2] = 0.24845199749997662/jacob_vx[2]-0.49690399499995325/jacob_vx[1]+0.24845199749997662/jacob_vx[0]; 
  bdiv_vx[3] = 0.37267799624996495/jacob_vx[2]-0.37267799624996495/jacob_vx[0]; 
  bdiv_vx[4] = 0.5/jacob_vx[2]+0.5/jacob_vx[0]; 
  bdiv_vx[5] = 0.3333333333333333/jacob_vx[2]-0.3333333333333333/jacob_vx[0]; 
  bdiv_vx[6] = 0.24845199749997662/jacob_vx[2]-0.49690399499995325/jacob_vx[1]+0.24845199749997662/jacob_vx[0]; 
  bdiv_vx[7] = 0.3333333333333333/jacob_vx[2]-0.3333333333333333/jacob_vx[0]; 
  bdiv_vx[8] = 0.2222222222222222/jacob_vx[2]+0.5555555555555556/jacob_vx[1]+0.2222222222222222/jacob_vx[0]; 
  double mb_vx[9]; 
  mb_vx[0] = 1.9364916731037085*bdiv_vx[3]*hamil[4]+0.8660254037844386*bdiv_vx[0]*hamil[1]; 
  mb_vx[1] = 1.9364916731037085*bdiv_vx[4]*hamil[4]+0.8660254037844386*bdiv_vx[1]*hamil[1]; 
  mb_vx[2] = 1.9364916731037085*hamil[4]*bdiv_vx[5]+0.8660254037844386*hamil[1]*bdiv_vx[2]; 
  mb_vx[3] = 1.7320508075688772*hamil[4]*bdiv_vx[6]+1.9364916731037085*bdiv_vx[0]*hamil[4]+0.8660254037844386*hamil[1]*bdiv_vx[3]; 
  mb_vx[4] = 1.7320508075688772*hamil[4]*bdiv_vx[7]+1.9364916731037085*bdiv_vx[1]*hamil[4]+0.8660254037844386*hamil[1]*bdiv_vx[4]; 
  mb_vx[5] = 1.7320508075688772*hamil[4]*bdiv_vx[8]+0.8660254037844386*hamil[1]*bdiv_vx[5]+1.9364916731037085*bdiv_vx[2]*hamil[4]; 
  mb_vx[6] = 0.8660254037844386*hamil[1]*bdiv_vx[6]+1.7320508075688772*bdiv_vx[3]*hamil[4]; 
  mb_vx[7] = 0.8660254037844386*hamil[1]*bdiv_vx[7]+1.7320508075688772*bdiv_vx[4]*hamil[4]; 
  mb_vx[8] = 0.8660254037844386*hamil[1]*bdiv_vx[8]+1.7320508075688772*hamil[4]*bdiv_vx[5]; 
  const double jacob_cx_edge_inv = 1.0/jacob_pos_edge[0]; 
  const double jacob_cx_skin_inv = 1.0/jacob_pos_skin[0]; 
  double u[9]; 
  if (edge == -1) { 
  if (wv>0) { 
  u[0] = 1.224744871391589*fskin[1]*jacob_cx_skin_inv+0.7071067811865475*fskin[0]*jacob_cx_skin_inv; 
  u[1] = 1.224744871391589*fskin[5]*jacob_cx_skin_inv+0.7071067811865475*fskin[3]*jacob_cx_skin_inv; 
  u[2] = 1.224744871391589*fskin[12]*jacob_cx_skin_inv+0.7071067811865475*fskin[8]*jacob_cx_skin_inv; 
  u[3] = 1.224744871391589*fskin[4]*jacob_cx_skin_inv+0.7071067811865475*fskin[2]*jacob_cx_skin_inv; 
  u[4] = 1.224744871391589*fskin[9]*jacob_cx_skin_inv+0.7071067811865475*fskin[6]*jacob_cx_skin_inv; 
  u[5] = 1.224744871391589*fskin[15]*jacob_cx_skin_inv+0.7071067811865475*fskin[13]*jacob_cx_skin_inv; 
  u[6] = 1.224744871391589*fskin[10]*jacob_cx_skin_inv+0.7071067811865475*fskin[7]*jacob_cx_skin_inv; 
  u[7] = 1.224744871391589*fskin[14]*jacob_cx_skin_inv+0.7071067811865475*fskin[11]*jacob_cx_skin_inv; 
  u[8] = 1.224744871391589*fskin[17]*jacob_cx_skin_inv+0.7071067811865475*fskin[16]*jacob_cx_skin_inv; 
  } else { 
  u[0] = 0.7071067811865475*fedge[0]*jacob_cx_edge_inv-1.224744871391589*fedge[1]*jacob_cx_edge_inv; 
  u[1] = 0.7071067811865475*fedge[3]*jacob_cx_edge_inv-1.224744871391589*fedge[5]*jacob_cx_edge_inv; 
  u[2] = 0.7071067811865475*fedge[8]*jacob_cx_edge_inv-1.224744871391589*fedge[12]*jacob_cx_edge_inv; 
  u[3] = 0.7071067811865475*fedge[2]*jacob_cx_edge_inv-1.224744871391589*fedge[4]*jacob_cx_edge_inv; 
  u[4] = 0.7071067811865475*fedge[6]*jacob_cx_edge_inv-1.224744871391589*fedge[9]*jacob_cx_edge_inv; 
  u[5] = 0.7071067811865475*fedge[13]*jacob_cx_edge_inv-1.224744871391589*fedge[15]*jacob_cx_edge_inv; 
  u[6] = 0.7071067811865475*fedge[7]*jacob_cx_edge_inv-1.224744871391589*fedge[10]*jacob_cx_edge_inv; 
  u[7] = 0.7071067811865475*fedge[11]*jacob_cx_edge_inv-1.224744871391589*fedge[14]*jacob_cx_edge_inv; 
  u[8] = 0.7071067811865475*fedge[16]*jacob_cx_edge_inv-1.224744871391589*fedge[17]*jacob_cx_edge_inv; 
  } 
  double Ghat[9]; 
  Ghat[0] = mb_vx[2]*u[6]+mb_vx[1]*u[3]+mb_vx[0]*u[0]; 
  Ghat[1] = mb_vx[2]*u[7]+mb_vx[1]*u[4]+mb_vx[0]*u[1]; 
  Ghat[2] = mb_vx[2]*u[8]+mb_vx[1]*u[5]+mb_vx[0]*u[2]; 
  Ghat[3] = mb_vx[5]*u[6]+u[3]*mb_vx[4]+u[0]*mb_vx[3]; 
  Ghat[4] = mb_vx[5]*u[7]+mb_vx[4]*u[4]+u[1]*mb_vx[3]; 
  Ghat[5] = mb_vx[5]*u[8]+mb_vx[4]*u[5]+u[2]*mb_vx[3]; 
  Ghat[6] = u[6]*mb_vx[8]+u[3]*mb_vx[7]+u[0]*mb_vx[6]; 
  Ghat[7] = u[7]*mb_vx[8]+u[4]*mb_vx[7]+u[1]*mb_vx[6]; 
  Ghat[8] = mb_vx[8]*u[8]+u[5]*mb_vx[7]+u[2]*mb_vx[6]; 
  out[0] += -(0.7071067811865475*Ghat[0]*dv10*dx10); 
  out[1] += -(1.224744871391589*Ghat[0]*dv10*dx10); 
  out[2] += -(0.7071067811865475*Ghat[3]*dv10*dx10); 
  out[3] += -(0.7071067811865475*Ghat[1]*dv10*dx10); 
  out[4] += -(1.224744871391589*Ghat[3]*dv10*dx10); 
  out[5] += -(1.224744871391589*Ghat[1]*dv10*dx10); 
  out[6] += -(0.7071067811865475*Ghat[4]*dv10*dx10); 
  out[7] += -(0.7071067811865475*Ghat[6]*dv10*dx10); 
  out[8] += -(0.7071067811865475*Ghat[2]*dv10*dx10); 
  out[9] += -(1.224744871391589*Ghat[4]*dv10*dx10); 
  out[10] += -(1.224744871391589*Ghat[6]*dv10*dx10); 
  out[11] += -(0.7071067811865475*Ghat[7]*dv10*dx10); 
  out[12] += -(1.224744871391589*Ghat[2]*dv10*dx10); 
  out[13] += -(0.7071067811865475*Ghat[5]*dv10*dx10); 
  out[14] += -(1.224744871391589*Ghat[7]*dv10*dx10); 
  out[15] += -(1.224744871391589*Ghat[5]*dv10*dx10); 
  out[16] += -(0.7071067811865475*Ghat[8]*dv10*dx10); 
  out[17] += -(1.224744871391589*Ghat[8]*dv10*dx10); 
  } else { 
  if (wv>0) { 
  u[0] = 1.224744871391589*fedge[1]*jacob_cx_edge_inv+0.7071067811865475*fedge[0]*jacob_cx_edge_inv; 
  u[1] = 1.224744871391589*fedge[5]*jacob_cx_edge_inv+0.7071067811865475*fedge[3]*jacob_cx_edge_inv; 
  u[2] = 1.224744871391589*fedge[12]*jacob_cx_edge_inv+0.7071067811865475*fedge[8]*jacob_cx_edge_inv; 
  u[3] = 1.224744871391589*fedge[4]*jacob_cx_edge_inv+0.7071067811865475*fedge[2]*jacob_cx_edge_inv; 
  u[4] = 1.224744871391589*fedge[9]*jacob_cx_edge_inv+0.7071067811865475*fedge[6]*jacob_cx_edge_inv; 
  u[5] = 1.224744871391589*fedge[15]*jacob_cx_edge_inv+0.7071067811865475*fedge[13]*jacob_cx_edge_inv; 
  u[6] = 1.224744871391589*fedge[10]*jacob_cx_edge_inv+0.7071067811865475*fedge[7]*jacob_cx_edge_inv; 
  u[7] = 1.224744871391589*fedge[14]*jacob_cx_edge_inv+0.7071067811865475*fedge[11]*jacob_cx_edge_inv; 
  u[8] = 1.224744871391589*fedge[17]*jacob_cx_edge_inv+0.7071067811865475*fedge[16]*jacob_cx_edge_inv; 
  } else { 
  u[0] = 0.7071067811865475*fskin[0]*jacob_cx_skin_inv-1.224744871391589*fskin[1]*jacob_cx_skin_inv; 
  u[1] = 0.7071067811865475*fskin[3]*jacob_cx_skin_inv-1.224744871391589*fskin[5]*jacob_cx_skin_inv; 
  u[2] = 0.7071067811865475*fskin[8]*jacob_cx_skin_inv-1.224744871391589*fskin[12]*jacob_cx_skin_inv; 
  u[3] = 0.7071067811865475*fskin[2]*jacob_cx_skin_inv-1.224744871391589*fskin[4]*jacob_cx_skin_inv; 
  u[4] = 0.7071067811865475*fskin[6]*jacob_cx_skin_inv-1.224744871391589*fskin[9]*jacob_cx_skin_inv; 
  u[5] = 0.7071067811865475*fskin[13]*jacob_cx_skin_inv-1.224744871391589*fskin[15]*jacob_cx_skin_inv; 
  u[6] = 0.7071067811865475*fskin[7]*jacob_cx_skin_inv-1.224744871391589*fskin[10]*jacob_cx_skin_inv; 
  u[7] = 0.7071067811865475*fskin[11]*jacob_cx_skin_inv-1.224744871391589*fskin[14]*jacob_cx_skin_inv; 
  u[8] = 0.7071067811865475*fskin[16]*jacob_cx_skin_inv-1.224744871391589*fskin[17]*jacob_cx_skin_inv; 
  } 
  double Ghat[9]; 
  Ghat[0] = mb_vx[2]*u[6]+mb_vx[1]*u[3]+mb_vx[0]*u[0]; 
  Ghat[1] = mb_vx[2]*u[7]+mb_vx[1]*u[4]+mb_vx[0]*u[1]; 
  Ghat[2] = mb_vx[2]*u[8]+mb_vx[1]*u[5]+mb_vx[0]*u[2]; 
  Ghat[3] = mb_vx[5]*u[6]+u[3]*mb_vx[4]+u[0]*mb_vx[3]; 
  Ghat[4] = mb_vx[5]*u[7]+mb_vx[4]*u[4]+u[1]*mb_vx[3]; 
  Ghat[5] = mb_vx[5]*u[8]+mb_vx[4]*u[5]+u[2]*mb_vx[3]; 
  Ghat[6] = u[6]*mb_vx[8]+u[3]*mb_vx[7]+u[0]*mb_vx[6]; 
  Ghat[7] = u[7]*mb_vx[8]+u[4]*mb_vx[7]+u[1]*mb_vx[6]; 
  Ghat[8] = mb_vx[8]*u[8]+u[5]*mb_vx[7]+u[2]*mb_vx[6]; 
  out[0] += 0.7071067811865475*Ghat[0]*dv10*dx10; 
  out[1] += -(1.224744871391589*Ghat[0]*dv10*dx10); 
  out[2] += 0.7071067811865475*Ghat[3]*dv10*dx10; 
  out[3] += 0.7071067811865475*Ghat[1]*dv10*dx10; 
  out[4] += -(1.224744871391589*Ghat[3]*dv10*dx10); 
  out[5] += -(1.224744871391589*Ghat[1]*dv10*dx10); 
  out[6] += 0.7071067811865475*Ghat[4]*dv10*dx10; 
  out[7] += 0.7071067811865475*Ghat[6]*dv10*dx10; 
  out[8] += 0.7071067811865475*Ghat[2]*dv10*dx10; 
  out[9] += -(1.224744871391589*Ghat[4]*dv10*dx10); 
  out[10] += -(1.224744871391589*Ghat[6]*dv10*dx10); 
  out[11] += 0.7071067811865475*Ghat[7]*dv10*dx10; 
  out[12] += -(1.224744871391589*Ghat[2]*dv10*dx10); 
  out[13] += 0.7071067811865475*Ghat[5]*dv10*dx10; 
  out[14] += -(1.224744871391589*Ghat[7]*dv10*dx10); 
  out[15] += -(1.224744871391589*Ghat[5]*dv10*dx10); 
  out[16] += 0.7071067811865475*Ghat[8]*dv10*dx10; 
  out[17] += -(1.224744871391589*Ghat[8]*dv10*dx10); 
  } 
  return 0.0;

} 
