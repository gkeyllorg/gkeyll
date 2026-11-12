#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_vel_sparse_surfx_1x2v_tensor_p1(const double *w, const double *dxv,
  const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel,
  const double *poisson_tensor_conf, const double *hamil,
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // Upwinded f on the two faces, divided by the velocity-map Jacobian along v_x and the 
  // configuration-space one along x, times dH/dv_x (velocity-space operations factored by direction). 
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
  const double jacob_cx_l_inv = 1.0/jacob_pos_l[0]; 
  const double jacob_cx_c_inv = 1.0/jacob_pos_c[0]; 
  const double jacob_cx_r_inv = 1.0/jacob_pos_r[0]; 
  double u_l[9], u_r[9]; 
  if (wv>0) { 
  u_r[0] = 1.224744871391589*fc[1]*jacob_cx_c_inv+0.7071067811865475*fc[0]*jacob_cx_c_inv; 
  u_r[1] = 1.224744871391589*fc[5]*jacob_cx_c_inv+0.7071067811865475*fc[3]*jacob_cx_c_inv; 
  u_r[2] = 1.224744871391589*fc[12]*jacob_cx_c_inv+0.7071067811865475*fc[8]*jacob_cx_c_inv; 
  u_r[3] = 1.224744871391589*fc[4]*jacob_cx_c_inv+0.7071067811865475*fc[2]*jacob_cx_c_inv; 
  u_r[4] = 1.224744871391589*fc[9]*jacob_cx_c_inv+0.7071067811865475*fc[6]*jacob_cx_c_inv; 
  u_r[5] = 1.224744871391589*fc[15]*jacob_cx_c_inv+0.7071067811865475*fc[13]*jacob_cx_c_inv; 
  u_r[6] = 1.224744871391589*fc[10]*jacob_cx_c_inv+0.7071067811865475*fc[7]*jacob_cx_c_inv; 
  u_r[7] = 1.224744871391589*fc[14]*jacob_cx_c_inv+0.7071067811865475*fc[11]*jacob_cx_c_inv; 
  u_r[8] = 1.224744871391589*fc[17]*jacob_cx_c_inv+0.7071067811865475*fc[16]*jacob_cx_c_inv; 
  u_l[0] = 1.224744871391589*fl[1]*jacob_cx_l_inv+0.7071067811865475*fl[0]*jacob_cx_l_inv; 
  u_l[1] = 1.224744871391589*fl[5]*jacob_cx_l_inv+0.7071067811865475*fl[3]*jacob_cx_l_inv; 
  u_l[2] = 1.224744871391589*fl[12]*jacob_cx_l_inv+0.7071067811865475*fl[8]*jacob_cx_l_inv; 
  u_l[3] = 1.224744871391589*fl[4]*jacob_cx_l_inv+0.7071067811865475*fl[2]*jacob_cx_l_inv; 
  u_l[4] = 1.224744871391589*fl[9]*jacob_cx_l_inv+0.7071067811865475*fl[6]*jacob_cx_l_inv; 
  u_l[5] = 1.224744871391589*fl[15]*jacob_cx_l_inv+0.7071067811865475*fl[13]*jacob_cx_l_inv; 
  u_l[6] = 1.224744871391589*fl[10]*jacob_cx_l_inv+0.7071067811865475*fl[7]*jacob_cx_l_inv; 
  u_l[7] = 1.224744871391589*fl[14]*jacob_cx_l_inv+0.7071067811865475*fl[11]*jacob_cx_l_inv; 
  u_l[8] = 1.224744871391589*fl[17]*jacob_cx_l_inv+0.7071067811865475*fl[16]*jacob_cx_l_inv; 
  } else { 
  u_r[0] = 0.7071067811865475*fr[0]*jacob_cx_r_inv-1.224744871391589*fr[1]*jacob_cx_r_inv; 
  u_r[1] = 0.7071067811865475*fr[3]*jacob_cx_r_inv-1.224744871391589*fr[5]*jacob_cx_r_inv; 
  u_r[2] = 0.7071067811865475*fr[8]*jacob_cx_r_inv-1.224744871391589*fr[12]*jacob_cx_r_inv; 
  u_r[3] = 0.7071067811865475*fr[2]*jacob_cx_r_inv-1.224744871391589*fr[4]*jacob_cx_r_inv; 
  u_r[4] = 0.7071067811865475*fr[6]*jacob_cx_r_inv-1.224744871391589*fr[9]*jacob_cx_r_inv; 
  u_r[5] = 0.7071067811865475*fr[13]*jacob_cx_r_inv-1.224744871391589*fr[15]*jacob_cx_r_inv; 
  u_r[6] = 0.7071067811865475*fr[7]*jacob_cx_r_inv-1.224744871391589*fr[10]*jacob_cx_r_inv; 
  u_r[7] = 0.7071067811865475*fr[11]*jacob_cx_r_inv-1.224744871391589*fr[14]*jacob_cx_r_inv; 
  u_r[8] = 0.7071067811865475*fr[16]*jacob_cx_r_inv-1.224744871391589*fr[17]*jacob_cx_r_inv; 
  u_l[0] = 0.7071067811865475*fc[0]*jacob_cx_c_inv-1.224744871391589*fc[1]*jacob_cx_c_inv; 
  u_l[1] = 0.7071067811865475*fc[3]*jacob_cx_c_inv-1.224744871391589*fc[5]*jacob_cx_c_inv; 
  u_l[2] = 0.7071067811865475*fc[8]*jacob_cx_c_inv-1.224744871391589*fc[12]*jacob_cx_c_inv; 
  u_l[3] = 0.7071067811865475*fc[2]*jacob_cx_c_inv-1.224744871391589*fc[4]*jacob_cx_c_inv; 
  u_l[4] = 0.7071067811865475*fc[6]*jacob_cx_c_inv-1.224744871391589*fc[9]*jacob_cx_c_inv; 
  u_l[5] = 0.7071067811865475*fc[13]*jacob_cx_c_inv-1.224744871391589*fc[15]*jacob_cx_c_inv; 
  u_l[6] = 0.7071067811865475*fc[7]*jacob_cx_c_inv-1.224744871391589*fc[10]*jacob_cx_c_inv; 
  u_l[7] = 0.7071067811865475*fc[11]*jacob_cx_c_inv-1.224744871391589*fc[14]*jacob_cx_c_inv; 
  u_l[8] = 0.7071067811865475*fc[16]*jacob_cx_c_inv-1.224744871391589*fc[17]*jacob_cx_c_inv; 
  } 
  double Ghat_l[9]; 
  Ghat_l[0] = mb_vx[2]*u_l[6]+mb_vx[1]*u_l[3]+mb_vx[0]*u_l[0]; 
  Ghat_l[1] = mb_vx[2]*u_l[7]+mb_vx[1]*u_l[4]+mb_vx[0]*u_l[1]; 
  Ghat_l[2] = mb_vx[2]*u_l[8]+mb_vx[1]*u_l[5]+mb_vx[0]*u_l[2]; 
  Ghat_l[3] = mb_vx[5]*u_l[6]+u_l[3]*mb_vx[4]+u_l[0]*mb_vx[3]; 
  Ghat_l[4] = mb_vx[5]*u_l[7]+mb_vx[4]*u_l[4]+u_l[1]*mb_vx[3]; 
  Ghat_l[5] = mb_vx[5]*u_l[8]+mb_vx[4]*u_l[5]+u_l[2]*mb_vx[3]; 
  Ghat_l[6] = u_l[6]*mb_vx[8]+u_l[3]*mb_vx[7]+u_l[0]*mb_vx[6]; 
  Ghat_l[7] = u_l[7]*mb_vx[8]+u_l[4]*mb_vx[7]+u_l[1]*mb_vx[6]; 
  Ghat_l[8] = mb_vx[8]*u_l[8]+u_l[5]*mb_vx[7]+u_l[2]*mb_vx[6]; 
  double Ghat_r[9]; 
  Ghat_r[0] = mb_vx[2]*u_r[6]+mb_vx[1]*u_r[3]+mb_vx[0]*u_r[0]; 
  Ghat_r[1] = mb_vx[2]*u_r[7]+mb_vx[1]*u_r[4]+mb_vx[0]*u_r[1]; 
  Ghat_r[2] = mb_vx[2]*u_r[8]+mb_vx[1]*u_r[5]+mb_vx[0]*u_r[2]; 
  Ghat_r[3] = mb_vx[5]*u_r[6]+u_r[3]*mb_vx[4]+u_r[0]*mb_vx[3]; 
  Ghat_r[4] = mb_vx[5]*u_r[7]+mb_vx[4]*u_r[4]+u_r[1]*mb_vx[3]; 
  Ghat_r[5] = mb_vx[5]*u_r[8]+mb_vx[4]*u_r[5]+u_r[2]*mb_vx[3]; 
  Ghat_r[6] = u_r[6]*mb_vx[8]+u_r[3]*mb_vx[7]+u_r[0]*mb_vx[6]; 
  Ghat_r[7] = u_r[7]*mb_vx[8]+u_r[4]*mb_vx[7]+u_r[1]*mb_vx[6]; 
  Ghat_r[8] = mb_vx[8]*u_r[8]+u_r[5]*mb_vx[7]+u_r[2]*mb_vx[6]; 
  out[0] += 0.7071067811865475*Ghat_l[0]*dv10*dx10-0.7071067811865475*Ghat_r[0]*dv10*dx10; 
  out[1] += -(1.224744871391589*Ghat_r[0]*dv10*dx10)-1.224744871391589*Ghat_l[0]*dv10*dx10; 
  out[2] += 0.7071067811865475*Ghat_l[3]*dv10*dx10-0.7071067811865475*Ghat_r[3]*dv10*dx10; 
  out[3] += 0.7071067811865475*Ghat_l[1]*dv10*dx10-0.7071067811865475*Ghat_r[1]*dv10*dx10; 
  out[4] += -(1.224744871391589*Ghat_r[3]*dv10*dx10)-1.224744871391589*Ghat_l[3]*dv10*dx10; 
  out[5] += -(1.224744871391589*Ghat_r[1]*dv10*dx10)-1.224744871391589*Ghat_l[1]*dv10*dx10; 
  out[6] += 0.7071067811865475*Ghat_l[4]*dv10*dx10-0.7071067811865475*Ghat_r[4]*dv10*dx10; 
  out[7] += 0.7071067811865475*Ghat_l[6]*dv10*dx10-0.7071067811865475*Ghat_r[6]*dv10*dx10; 
  out[8] += 0.7071067811865475*Ghat_l[2]*dv10*dx10-0.7071067811865475*Ghat_r[2]*dv10*dx10; 
  out[9] += -(1.224744871391589*Ghat_r[4]*dv10*dx10)-1.224744871391589*Ghat_l[4]*dv10*dx10; 
  out[10] += -(1.224744871391589*Ghat_r[6]*dv10*dx10)-1.224744871391589*Ghat_l[6]*dv10*dx10; 
  out[11] += 0.7071067811865475*Ghat_l[7]*dv10*dx10-0.7071067811865475*Ghat_r[7]*dv10*dx10; 
  out[12] += -(1.224744871391589*Ghat_r[2]*dv10*dx10)-1.224744871391589*Ghat_l[2]*dv10*dx10; 
  out[13] += 0.7071067811865475*Ghat_l[5]*dv10*dx10-0.7071067811865475*Ghat_r[5]*dv10*dx10; 
  out[14] += -(1.224744871391589*Ghat_r[7]*dv10*dx10)-1.224744871391589*Ghat_l[7]*dv10*dx10; 
  out[15] += -(1.224744871391589*Ghat_r[5]*dv10*dx10)-1.224744871391589*Ghat_l[5]*dv10*dx10; 
  out[16] += 0.7071067811865475*Ghat_l[8]*dv10*dx10-0.7071067811865475*Ghat_r[8]*dv10*dx10; 
  out[17] += -(1.224744871391589*Ghat_r[8]*dv10*dx10)-1.224744871391589*Ghat_l[8]*dv10*dx10; 
  return fabs(0.75*dv10*dx10*jacob_cx_c_inv*((1.4433756729740645*hamil[4])/jacob_vx[2]-(1.4433756729740645*hamil[4])/jacob_vx[0]+(0.48112522432468824*hamil[1])/jacob_vx[2]+(0.7698003589195012*hamil[1])/jacob_vx[1]+(0.48112522432468824*hamil[1])/jacob_vx[0]));

} 
