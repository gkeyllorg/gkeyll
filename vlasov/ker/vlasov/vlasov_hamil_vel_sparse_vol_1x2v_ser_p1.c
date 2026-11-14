#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_hamil_vel_sparse_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap,
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
  out[1] += 1.5*f[0]*hamil[1]*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[4] += 1.5*hamil[1]*f[2]*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[5] += 1.5*hamil[1]*f[3]*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[7] += 1.5*hamil[1]*f[6]*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 

} 
