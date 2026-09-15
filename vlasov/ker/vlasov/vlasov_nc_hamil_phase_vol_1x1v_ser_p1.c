#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_phase_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *vmap,
  const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f,
  const double *f_no_J, double* GKYL_RESTRICT out) 
{ 
  const double dx10 = 2.0/dxv[0]; 
  const double dv0 = dxv[1]; 
  const double dv10 = 2.0/dxv[1]; 
  const double jacob_cx_inv = 1.0/jacob_pos[0]; 
  const double jacob_vx_inv = 1.0/jacob_vel[0]; 
  const double wx1 = w[1]; 
  const double dv10_sq = dv10*dv10; 
  
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[0]; 
  out[1] += (1.060660171779821*f[0]*poisson_tensor_conf_0[1]*hamil[3]+1.060660171779821*poisson_tensor_conf_0[0]*f[1]*hamil[3]+1.060660171779821*f[1]*poisson_tensor_conf_0[1]*hamil[2]+1.060660171779821*f[0]*poisson_tensor_conf_0[0]*hamil[2])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[2] += (-(1.060660171779821*poisson_tensor_conf_0[1]*f[3]*hamil[3])-1.060660171779821*poisson_tensor_conf_0[0]*f[2]*hamil[3]-1.060660171779821*f[1]*hamil[1]*poisson_tensor_conf_0[1]-1.060660171779821*f[0]*poisson_tensor_conf_0[0]*hamil[1])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[3] += (1.060660171779821*poisson_tensor_conf_0[1]*hamil[2]*f[3]+1.060660171779821*poisson_tensor_conf_0[0]*f[2]*hamil[2]-1.060660171779821*f[0]*hamil[1]*poisson_tensor_conf_0[1]-1.060660171779821*poisson_tensor_conf_0[0]*f[1]*hamil[1])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 

} 
