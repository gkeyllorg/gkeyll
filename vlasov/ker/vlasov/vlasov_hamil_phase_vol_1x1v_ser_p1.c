#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_hamil_phase_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *vmap,
  const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil,
  const double *f,
  const double *f_no_J, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[1]; 
  out[1] += (1.5*f[1]*hamil[3]+1.5*f[0]*hamil[2])*dv10*dx10; 
  out[2] += (-(1.5*f[2]*hamil[3])-1.5*f[0]*hamil[1])*dv10*dx10; 
  out[3] += (1.5*f[2]*hamil[2]-1.5*f[1]*hamil[1])*dv10*dx10; 

} 
