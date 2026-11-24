#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_2x1v_tensor_p1(const double *w, const double *dxv, const double *vmap,
  const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil,
  const double *f,
  const double *f_no_J, double* GKYL_RESTRICT out) 
{ 
  const double dx10 = 2.0/dxv[0]; 
  const double dv10 = 2.0/dxv[2]; 
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
  const double jacob_cx_inv = 1.0/jacob_pos[0]; 
  double G_vx_s0_a[4]; 
  G_vx_s0_a[0] = eb_vx[2]*f[7]*jacob_cx_inv+eb_vx[1]*f[3]*jacob_cx_inv+eb_vx[0]*f[0]*jacob_cx_inv; 
  G_vx_s0_a[1] = eb_vx[5]*f[7]*jacob_cx_inv+f[3]*eb_vx[4]*jacob_cx_inv+f[0]*eb_vx[3]*jacob_cx_inv; 
  G_vx_s0_a[2] = f[7]*eb_vx[8]*jacob_cx_inv+f[3]*eb_vx[7]*jacob_cx_inv+f[0]*eb_vx[6]*jacob_cx_inv; 
  G_vx_s0_a[3] = f[7]*eb_vx[11]*jacob_cx_inv+f[3]*eb_vx[10]*jacob_cx_inv+f[0]*eb_vx[9]*jacob_cx_inv; 
  double G_vx_s0_n[4]; 
  G_vx_s0_n[0] = G_vx_s0_a[0]*gn_vx[0]; 
  G_vx_s0_n[1] = G_vx_s0_a[1]*gn_vx[1]; 
  G_vx_s0_n[2] = G_vx_s0_a[2]*gn_vx[2]; 
  G_vx_s0_n[3] = G_vx_s0_a[3]*gn_vx[3]; 
  double G_vx_s1_a[4]; 
  G_vx_s1_a[0] = eb_vx[2]*f[10]*jacob_cx_inv+eb_vx[1]*f[6]*jacob_cx_inv+eb_vx[0]*f[2]*jacob_cx_inv; 
  G_vx_s1_a[1] = eb_vx[5]*f[10]*jacob_cx_inv+eb_vx[4]*f[6]*jacob_cx_inv+f[2]*eb_vx[3]*jacob_cx_inv; 
  G_vx_s1_a[2] = eb_vx[8]*f[10]*jacob_cx_inv+f[6]*eb_vx[7]*jacob_cx_inv+f[2]*eb_vx[6]*jacob_cx_inv; 
  G_vx_s1_a[3] = f[10]*eb_vx[11]*jacob_cx_inv+f[6]*eb_vx[10]*jacob_cx_inv+f[2]*eb_vx[9]*jacob_cx_inv; 
  double G_vx_s1_n[4]; 
  G_vx_s1_n[0] = G_vx_s1_a[0]*gn_vx[0]; 
  G_vx_s1_n[1] = G_vx_s1_a[1]*gn_vx[1]; 
  G_vx_s1_n[2] = G_vx_s1_a[2]*gn_vx[2]; 
  G_vx_s1_n[3] = G_vx_s1_a[3]*gn_vx[3]; 
  double G_vx[6]; 
  G_vx[0] = 0.24597051986528992*G_vx_s0_n[3]+0.4611362613212574*G_vx_s0_n[2]+0.4611362613212574*G_vx_s0_n[1]+0.24597051986528992*G_vx_s0_n[0]; 
  G_vx[1] = 0.36687286304546396*G_vx_s0_n[3]+0.2715467467935445*G_vx_s0_n[2]-0.2715467467935445*G_vx_s0_n[1]-0.36687286304546396*G_vx_s0_n[0]; 
  G_vx[2] = 0.33678765702728136*G_vx_s0_n[3]-0.33678765702728136*G_vx_s0_n[2]-0.33678765702728136*G_vx_s0_n[1]+0.33678765702728136*G_vx_s0_n[0]; 
  G_vx[3] = 0.24597051986528992*G_vx_s1_n[3]+0.4611362613212574*G_vx_s1_n[2]+0.4611362613212574*G_vx_s1_n[1]+0.24597051986528992*G_vx_s1_n[0]; 
  G_vx[4] = 0.36687286304546396*G_vx_s1_n[3]+0.2715467467935445*G_vx_s1_n[2]-0.2715467467935445*G_vx_s1_n[1]-0.36687286304546396*G_vx_s1_n[0]; 
  G_vx[5] = 0.33678765702728136*G_vx_s1_n[3]-0.33678765702728136*G_vx_s1_n[2]-0.33678765702728136*G_vx_s1_n[1]+0.33678765702728136*G_vx_s1_n[0]; 
  out[1] += 1.7320508075688772*G_vx[0]*dv10*dx10; 
  out[4] += 1.7320508075688772*G_vx[3]*dv10*dx10; 
  out[5] += 1.7320508075688772*G_vx[1]*dv10*dx10; 
  out[8] += 1.7320508075688772*G_vx[4]*dv10*dx10; 
  out[9] += 1.7320508075688772*G_vx[2]*dv10*dx10; 
  out[11] += 1.7320508075688772*G_vx[5]*dv10*dx10; 
} 
