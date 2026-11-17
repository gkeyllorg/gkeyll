#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfx_1x1v_tensor_p1(const double *w, const double *dxv,
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
  double eb_vx[12]; 
  eb_vx[0] = 0.9681844646844027*bdiv_vx[6]-1.0546722811938851*bdiv_vx[3]+0.7071067811865475*bdiv_vx[0]; 
  eb_vx[1] = 0.9681844646844027*bdiv_vx[7]-1.0546722811938851*bdiv_vx[4]+0.7071067811865475*bdiv_vx[1]; 
  eb_vx[2] = 0.9681844646844027*bdiv_vx[8]-1.0546722811938851*bdiv_vx[5]+0.7071067811865475*bdiv_vx[2]; 
  eb_vx[3] = -(0.516430513231777*bdiv_vx[6])-0.4163900395009129*bdiv_vx[3]+0.7071067811865475*bdiv_vx[0]; 
  eb_vx[4] = -(0.516430513231777*bdiv_vx[7])-0.4163900395009129*bdiv_vx[4]+0.7071067811865475*bdiv_vx[1]; 
  eb_vx[5] = -(0.516430513231777*bdiv_vx[8])-0.4163900395009129*bdiv_vx[5]+0.7071067811865475*bdiv_vx[2]; 
  eb_vx[6] = -(0.516430513231777*bdiv_vx[6])+0.4163900395009129*bdiv_vx[3]+0.7071067811865475*bdiv_vx[0]; 
  eb_vx[7] = -(0.516430513231777*bdiv_vx[7])+0.4163900395009129*bdiv_vx[4]+0.7071067811865475*bdiv_vx[1]; 
  eb_vx[8] = -(0.516430513231777*bdiv_vx[8])+0.4163900395009129*bdiv_vx[5]+0.7071067811865475*bdiv_vx[2]; 
  eb_vx[9] = 0.9681844646844027*bdiv_vx[6]+1.0546722811938851*bdiv_vx[3]+0.7071067811865475*bdiv_vx[0]; 
  eb_vx[10] = 0.9681844646844027*bdiv_vx[7]+1.0546722811938851*bdiv_vx[4]+0.7071067811865475*bdiv_vx[1]; 
  eb_vx[11] = 0.9681844646844027*bdiv_vx[8]+1.0546722811938851*bdiv_vx[5]+0.7071067811865475*bdiv_vx[2]; 
  double gn_vx[4]; 
  gn_vx[0] = 1.224744871391589*hamil[1]-4.084728180770504*hamil[2]; 
  gn_vx[1] = 1.224744871391589*hamil[1]-1.6126716885136843*hamil[2]; 
  gn_vx[2] = 1.6126716885136843*hamil[2]+1.224744871391589*hamil[1]; 
  gn_vx[3] = 4.084728180770504*hamil[2]+1.224744871391589*hamil[1]; 
  const double jacob_cx_l_inv = 1.0/jacob_pos_l[0]; 
  const double jacob_cx_c_inv = 1.0/jacob_pos_c[0]; 
  const double jacob_cx_r_inv = 1.0/jacob_pos_r[0]; 
  double u_l[3], u_r[3]; 
  if (wv>0) { 
  u_r[0] = 1.224744871391589*fc[1]*jacob_cx_c_inv+0.7071067811865475*fc[0]*jacob_cx_c_inv; 
  u_r[1] = 1.224744871391589*fc[3]*jacob_cx_c_inv+0.7071067811865475*fc[2]*jacob_cx_c_inv; 
  u_r[2] = 1.224744871391589*fc[5]*jacob_cx_c_inv+0.7071067811865475*fc[4]*jacob_cx_c_inv; 
  u_l[0] = 1.224744871391589*fl[1]*jacob_cx_l_inv+0.7071067811865475*fl[0]*jacob_cx_l_inv; 
  u_l[1] = 1.224744871391589*fl[3]*jacob_cx_l_inv+0.7071067811865475*fl[2]*jacob_cx_l_inv; 
  u_l[2] = 1.224744871391589*fl[5]*jacob_cx_l_inv+0.7071067811865475*fl[4]*jacob_cx_l_inv; 
  } else { 
  u_r[0] = 0.7071067811865475*fr[0]*jacob_cx_r_inv-1.224744871391589*fr[1]*jacob_cx_r_inv; 
  u_r[1] = 0.7071067811865475*fr[2]*jacob_cx_r_inv-1.224744871391589*fr[3]*jacob_cx_r_inv; 
  u_r[2] = 0.7071067811865475*fr[4]*jacob_cx_r_inv-1.224744871391589*fr[5]*jacob_cx_r_inv; 
  u_l[0] = 0.7071067811865475*fc[0]*jacob_cx_c_inv-1.224744871391589*fc[1]*jacob_cx_c_inv; 
  u_l[1] = 0.7071067811865475*fc[2]*jacob_cx_c_inv-1.224744871391589*fc[3]*jacob_cx_c_inv; 
  u_l[2] = 0.7071067811865475*fc[4]*jacob_cx_c_inv-1.224744871391589*fc[5]*jacob_cx_c_inv; 
  } 
  double Ghat_l_s0_a[4]; 
  Ghat_l_s0_a[0] = eb_vx[2]*u_l[2]+eb_vx[1]*u_l[1]+eb_vx[0]*u_l[0]; 
  Ghat_l_s0_a[1] = u_l[2]*eb_vx[5]+u_l[1]*eb_vx[4]+u_l[0]*eb_vx[3]; 
  Ghat_l_s0_a[2] = u_l[2]*eb_vx[8]+u_l[1]*eb_vx[7]+u_l[0]*eb_vx[6]; 
  Ghat_l_s0_a[3] = u_l[2]*eb_vx[11]+u_l[1]*eb_vx[10]+u_l[0]*eb_vx[9]; 
  double Ghat_l_s0_n[4]; 
  Ghat_l_s0_n[0] = Ghat_l_s0_a[0]*gn_vx[0]; 
  Ghat_l_s0_n[1] = Ghat_l_s0_a[1]*gn_vx[1]; 
  Ghat_l_s0_n[2] = Ghat_l_s0_a[2]*gn_vx[2]; 
  Ghat_l_s0_n[3] = Ghat_l_s0_a[3]*gn_vx[3]; 
  double Ghat_l[3]; 
  Ghat_l[0] = 0.24597051986528992*Ghat_l_s0_n[3]+0.4611362613212574*Ghat_l_s0_n[2]+0.4611362613212574*Ghat_l_s0_n[1]+0.24597051986528992*Ghat_l_s0_n[0]; 
  Ghat_l[1] = 0.36687286304546396*Ghat_l_s0_n[3]+0.2715467467935445*Ghat_l_s0_n[2]-0.2715467467935445*Ghat_l_s0_n[1]-0.36687286304546396*Ghat_l_s0_n[0]; 
  Ghat_l[2] = 0.33678765702728136*Ghat_l_s0_n[3]-0.33678765702728136*Ghat_l_s0_n[2]-0.33678765702728136*Ghat_l_s0_n[1]+0.33678765702728136*Ghat_l_s0_n[0]; 
  double Ghat_r_s0_a[4]; 
  Ghat_r_s0_a[0] = eb_vx[2]*u_r[2]+eb_vx[1]*u_r[1]+eb_vx[0]*u_r[0]; 
  Ghat_r_s0_a[1] = u_r[2]*eb_vx[5]+u_r[1]*eb_vx[4]+u_r[0]*eb_vx[3]; 
  Ghat_r_s0_a[2] = u_r[2]*eb_vx[8]+u_r[1]*eb_vx[7]+u_r[0]*eb_vx[6]; 
  Ghat_r_s0_a[3] = u_r[2]*eb_vx[11]+u_r[1]*eb_vx[10]+u_r[0]*eb_vx[9]; 
  double Ghat_r_s0_n[4]; 
  Ghat_r_s0_n[0] = Ghat_r_s0_a[0]*gn_vx[0]; 
  Ghat_r_s0_n[1] = Ghat_r_s0_a[1]*gn_vx[1]; 
  Ghat_r_s0_n[2] = Ghat_r_s0_a[2]*gn_vx[2]; 
  Ghat_r_s0_n[3] = Ghat_r_s0_a[3]*gn_vx[3]; 
  double Ghat_r[3]; 
  Ghat_r[0] = 0.24597051986528992*Ghat_r_s0_n[3]+0.4611362613212574*Ghat_r_s0_n[2]+0.4611362613212574*Ghat_r_s0_n[1]+0.24597051986528992*Ghat_r_s0_n[0]; 
  Ghat_r[1] = 0.36687286304546396*Ghat_r_s0_n[3]+0.2715467467935445*Ghat_r_s0_n[2]-0.2715467467935445*Ghat_r_s0_n[1]-0.36687286304546396*Ghat_r_s0_n[0]; 
  Ghat_r[2] = 0.33678765702728136*Ghat_r_s0_n[3]-0.33678765702728136*Ghat_r_s0_n[2]-0.33678765702728136*Ghat_r_s0_n[1]+0.33678765702728136*Ghat_r_s0_n[0]; 
  out[0] += 0.7071067811865475*Ghat_l[0]*dv10*dx10-0.7071067811865475*Ghat_r[0]*dv10*dx10; 
  out[1] += -(1.224744871391589*Ghat_r[0]*dv10*dx10)-1.224744871391589*Ghat_l[0]*dv10*dx10; 
  out[2] += 0.7071067811865475*Ghat_l[1]*dv10*dx10-0.7071067811865475*Ghat_r[1]*dv10*dx10; 
  out[3] += -(1.224744871391589*Ghat_r[1]*dv10*dx10)-1.224744871391589*Ghat_l[1]*dv10*dx10; 
  out[4] += 0.7071067811865475*Ghat_l[2]*dv10*dx10-0.7071067811865475*Ghat_r[2]*dv10*dx10; 
  out[5] += -(1.224744871391589*Ghat_r[2]*dv10*dx10)-1.224744871391589*Ghat_l[2]*dv10*dx10; 
  return fabs(1.0606601717798212*dv10*dx10*jacob_cx_c_inv*((1.4433756729740645*hamil[2])/jacob_vx[2]+(0.48112522432468824*hamil[1])/jacob_vx[2]-(1.4433756729740645*hamil[2])/jacob_vx[0]+(0.7698003589195012*hamil[1])/jacob_vx[1]+(0.48112522432468824*hamil[1])/jacob_vx[0]));

} 
