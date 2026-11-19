#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_hamil_vel_sparse_vol_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap,
  const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil,
  const double *f,
  const double *f_no_J, double* GKYL_RESTRICT out) 
{ 
  const double dx10 = 2.0/dxv[0]; 
  const double dv10 = 2.0/dxv[1]; 
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
  const double jacob_cx_inv = 1.0/jacob_pos[0]; 
  double G_vx[9]; 
  G_vx[0] = mb_vx[2]*f[7]*jacob_cx_inv+mb_vx[1]*f[2]*jacob_cx_inv+f[0]*mb_vx[0]*jacob_cx_inv; 
  G_vx[1] = mb_vx[2]*f[11]*jacob_cx_inv+mb_vx[1]*f[6]*jacob_cx_inv+mb_vx[0]*f[3]*jacob_cx_inv; 
  G_vx[2] = mb_vx[2]*f[16]*jacob_cx_inv+mb_vx[1]*f[13]*jacob_cx_inv+mb_vx[0]*f[8]*jacob_cx_inv; 
  G_vx[3] = mb_vx[5]*f[7]*jacob_cx_inv+f[2]*mb_vx[4]*jacob_cx_inv+f[0]*mb_vx[3]*jacob_cx_inv; 
  G_vx[4] = mb_vx[5]*f[11]*jacob_cx_inv+mb_vx[4]*f[6]*jacob_cx_inv+f[3]*mb_vx[3]*jacob_cx_inv; 
  G_vx[5] = mb_vx[5]*f[16]*jacob_cx_inv+mb_vx[4]*f[13]*jacob_cx_inv+mb_vx[3]*f[8]*jacob_cx_inv; 
  G_vx[6] = f[7]*mb_vx[8]*jacob_cx_inv+f[2]*mb_vx[7]*jacob_cx_inv+f[0]*mb_vx[6]*jacob_cx_inv; 
  G_vx[7] = mb_vx[8]*f[11]*jacob_cx_inv+f[6]*mb_vx[7]*jacob_cx_inv+f[3]*mb_vx[6]*jacob_cx_inv; 
  G_vx[8] = mb_vx[8]*f[16]*jacob_cx_inv+mb_vx[7]*f[13]*jacob_cx_inv+mb_vx[6]*f[8]*jacob_cx_inv; 
  out[1] += 1.7320508075688772*G_vx[0]*dv10*dx10; 
  out[4] += 1.7320508075688772*G_vx[3]*dv10*dx10; 
  out[5] += 1.7320508075688772*G_vx[1]*dv10*dx10; 
  out[9] += 1.7320508075688772*G_vx[4]*dv10*dx10; 
  out[10] += 1.7320508075688772*G_vx[6]*dv10*dx10; 
  out[12] += 1.7320508075688772*G_vx[2]*dv10*dx10; 
  out[14] += 1.7320508075688772*G_vx[7]*dv10*dx10; 
  out[15] += 1.7320508075688772*G_vx[5]*dv10*dx10; 
  out[17] += 1.7320508075688772*G_vx[8]*dv10*dx10; 
} 
