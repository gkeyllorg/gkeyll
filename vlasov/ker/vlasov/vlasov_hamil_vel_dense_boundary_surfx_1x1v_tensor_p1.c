#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_vel_dense_boundary_surfx_1x1v_tensor_p1(const double *w, const double *dxv,
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
  const double jacob_cx_edge_inv = 1.0/jacob_pos_edge[0]; 
  const double jacob_cx_skin_inv = 1.0/jacob_pos_skin[0]; 
  double u[3]; 
  if (edge == -1) { 
  if (wv>0) { 
  u[0] = 1.224744871391589*fskin[1]*jacob_cx_skin_inv+0.7071067811865475*fskin[0]*jacob_cx_skin_inv; 
  u[1] = 1.224744871391589*fskin[3]*jacob_cx_skin_inv+0.7071067811865475*fskin[2]*jacob_cx_skin_inv; 
  u[2] = 1.224744871391589*fskin[5]*jacob_cx_skin_inv+0.7071067811865475*fskin[4]*jacob_cx_skin_inv; 
  } else { 
  u[0] = 0.7071067811865475*fedge[0]*jacob_cx_edge_inv-1.224744871391589*fedge[1]*jacob_cx_edge_inv; 
  u[1] = 0.7071067811865475*fedge[2]*jacob_cx_edge_inv-1.224744871391589*fedge[3]*jacob_cx_edge_inv; 
  u[2] = 0.7071067811865475*fedge[4]*jacob_cx_edge_inv-1.224744871391589*fedge[5]*jacob_cx_edge_inv; 
  } 
  double Ghat_s0_a[4]; 
  Ghat_s0_a[0] = eb_vx[2]*u[2]+eb_vx[1]*u[1]+eb_vx[0]*u[0]; 
  Ghat_s0_a[1] = u[2]*eb_vx[5]+u[1]*eb_vx[4]+u[0]*eb_vx[3]; 
  Ghat_s0_a[2] = u[2]*eb_vx[8]+u[1]*eb_vx[7]+u[0]*eb_vx[6]; 
  Ghat_s0_a[3] = u[2]*eb_vx[11]+u[1]*eb_vx[10]+u[0]*eb_vx[9]; 
  double Ghat_s0_n[4]; 
  Ghat_s0_n[0] = Ghat_s0_a[0]*gn_vx[0]; 
  Ghat_s0_n[1] = Ghat_s0_a[1]*gn_vx[1]; 
  Ghat_s0_n[2] = Ghat_s0_a[2]*gn_vx[2]; 
  Ghat_s0_n[3] = Ghat_s0_a[3]*gn_vx[3]; 
  double Ghat[3]; 
  Ghat[0] = 0.24597051986528992*Ghat_s0_n[3]+0.4611362613212574*Ghat_s0_n[2]+0.4611362613212574*Ghat_s0_n[1]+0.24597051986528992*Ghat_s0_n[0]; 
  Ghat[1] = 0.36687286304546396*Ghat_s0_n[3]+0.2715467467935445*Ghat_s0_n[2]-0.2715467467935445*Ghat_s0_n[1]-0.36687286304546396*Ghat_s0_n[0]; 
  Ghat[2] = 0.33678765702728136*Ghat_s0_n[3]-0.33678765702728136*Ghat_s0_n[2]-0.33678765702728136*Ghat_s0_n[1]+0.33678765702728136*Ghat_s0_n[0]; 
  out[0] += -(0.7071067811865475*Ghat[0]*dv10*dx10); 
  out[1] += -(1.224744871391589*Ghat[0]*dv10*dx10); 
  out[2] += -(0.7071067811865475*Ghat[1]*dv10*dx10); 
  out[3] += -(1.224744871391589*Ghat[1]*dv10*dx10); 
  out[4] += -(0.7071067811865475*Ghat[2]*dv10*dx10); 
  out[5] += -(1.224744871391589*Ghat[2]*dv10*dx10); 
  } else { 
  if (wv>0) { 
  u[0] = 1.224744871391589*fedge[1]*jacob_cx_edge_inv+0.7071067811865475*fedge[0]*jacob_cx_edge_inv; 
  u[1] = 1.224744871391589*fedge[3]*jacob_cx_edge_inv+0.7071067811865475*fedge[2]*jacob_cx_edge_inv; 
  u[2] = 1.224744871391589*fedge[5]*jacob_cx_edge_inv+0.7071067811865475*fedge[4]*jacob_cx_edge_inv; 
  } else { 
  u[0] = 0.7071067811865475*fskin[0]*jacob_cx_skin_inv-1.224744871391589*fskin[1]*jacob_cx_skin_inv; 
  u[1] = 0.7071067811865475*fskin[2]*jacob_cx_skin_inv-1.224744871391589*fskin[3]*jacob_cx_skin_inv; 
  u[2] = 0.7071067811865475*fskin[4]*jacob_cx_skin_inv-1.224744871391589*fskin[5]*jacob_cx_skin_inv; 
  } 
  double Ghat_s0_a[4]; 
  Ghat_s0_a[0] = eb_vx[2]*u[2]+eb_vx[1]*u[1]+eb_vx[0]*u[0]; 
  Ghat_s0_a[1] = u[2]*eb_vx[5]+u[1]*eb_vx[4]+u[0]*eb_vx[3]; 
  Ghat_s0_a[2] = u[2]*eb_vx[8]+u[1]*eb_vx[7]+u[0]*eb_vx[6]; 
  Ghat_s0_a[3] = u[2]*eb_vx[11]+u[1]*eb_vx[10]+u[0]*eb_vx[9]; 
  double Ghat_s0_n[4]; 
  Ghat_s0_n[0] = Ghat_s0_a[0]*gn_vx[0]; 
  Ghat_s0_n[1] = Ghat_s0_a[1]*gn_vx[1]; 
  Ghat_s0_n[2] = Ghat_s0_a[2]*gn_vx[2]; 
  Ghat_s0_n[3] = Ghat_s0_a[3]*gn_vx[3]; 
  double Ghat[3]; 
  Ghat[0] = 0.24597051986528992*Ghat_s0_n[3]+0.4611362613212574*Ghat_s0_n[2]+0.4611362613212574*Ghat_s0_n[1]+0.24597051986528992*Ghat_s0_n[0]; 
  Ghat[1] = 0.36687286304546396*Ghat_s0_n[3]+0.2715467467935445*Ghat_s0_n[2]-0.2715467467935445*Ghat_s0_n[1]-0.36687286304546396*Ghat_s0_n[0]; 
  Ghat[2] = 0.33678765702728136*Ghat_s0_n[3]-0.33678765702728136*Ghat_s0_n[2]-0.33678765702728136*Ghat_s0_n[1]+0.33678765702728136*Ghat_s0_n[0]; 
  out[0] += 0.7071067811865475*Ghat[0]*dv10*dx10; 
  out[1] += -(1.224744871391589*Ghat[0]*dv10*dx10); 
  out[2] += 0.7071067811865475*Ghat[1]*dv10*dx10; 
  out[3] += -(1.224744871391589*Ghat[1]*dv10*dx10); 
  out[4] += 0.7071067811865475*Ghat[2]*dv10*dx10; 
  out[5] += -(1.224744871391589*Ghat[2]*dv10*dx10); 
  } 
  return 0.0;

} 
