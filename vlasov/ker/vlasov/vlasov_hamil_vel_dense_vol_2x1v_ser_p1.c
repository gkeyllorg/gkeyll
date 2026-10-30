#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_2x1v_ser_p1(const double *w, const double *dxv, const double *vmap,
  const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil,
  const double *f,
  const double *f_no_J, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  const double *jacob_cx = &jacob_pos[0]; 
  const double jacob_cx_inv = 1.0/jacob_cx[0]; 
  out[1] += 2.1213203435596424*f[0]*hamil[1]*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[4] += 2.1213203435596424*hamil[1]*f[2]*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[5] += 2.1213203435596424*hamil[1]*f[3]*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[7] += 2.1213203435596424*hamil[1]*f[6]*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 

} 
