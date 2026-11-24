#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_hamil_vel_sparse_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *vmap,
  const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil,
  const double *f,
  const double *f_no_J, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  const double *jacob_cx = &jacob_pos[0]; 
  const double jacob_cx_inv = 1.0/jacob_cx[0]; 
  out[1] += (3.3541019662496847*f[2]*hamil[4]+1.5*f[0]*hamil[1])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[4] += (3.0*hamil[4]*f[8]+3.3541019662496847*f[0]*hamil[4]+1.5*hamil[1]*f[2])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[5] += (3.3541019662496847*hamil[4]*f[6]+1.5*hamil[1]*f[3])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[7] += (7.5*f[4]*hamil[4]+3.3541019662496847*f[1]*hamil[1])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[10] += (3.0*hamil[4]*f[14]+1.5*hamil[1]*f[6]+3.3541019662496847*f[3]*hamil[4])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[11] += (6.708203932499369*hamil[4]*f[12]+7.500000000000001*f[1]*hamil[4]+3.3541019662496843*hamil[1]*f[4])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[12] += (1.5*hamil[1]*f[8]+3.0*f[2]*hamil[4])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[13] += (7.500000000000001*hamil[4]*f[10]+3.3541019662496843*hamil[1]*f[5])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[15] += (3.3541019662496847*hamil[4]*f[16]+1.5*hamil[1]*f[9])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[17] += (6.708203932499369*hamil[4]*f[18]+3.3541019662496847*hamil[1]*f[10]+7.5*hamil[4]*f[5])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[18] += (1.5*hamil[1]*f[14]+3.0*hamil[4]*f[6])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[19] += (1.5*hamil[1]*f[16]+3.3541019662496847*hamil[4]*f[9])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 

} 
