#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_hamil_phase_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *vmap,
  const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil,
  const double *f,
  const double *f_no_J, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[1]; 
  out[1] += (3.3541019662496843*f[3]*hamil[7]+1.5*f[4]*hamil[6]+3.3541019662496847*f[2]*hamil[5]+1.5*f[1]*hamil[3]+1.5*f[0]*hamil[2])*dv10*dx10; 
  out[2] += (-(1.5*f[5]*hamil[7])-3.3541019662496843*f[3]*hamil[6]-3.3541019662496847*f[1]*hamil[4]-1.5*f[2]*hamil[3]-1.5*f[0]*hamil[1])*dv10*dx10; 
  out[3] += (1.5*f[7]*hamil[7]+3.3541019662496843*f[1]*hamil[7]-1.5*f[6]*hamil[6]-3.3541019662496843*f[2]*hamil[6]+3.0*f[5]*hamil[5]+3.3541019662496847*f[0]*hamil[5]-3.0*f[4]*hamil[4]-3.3541019662496847*f[0]*hamil[4]+1.5*f[2]*hamil[2]-1.5*f[1]*hamil[1])*dv10*dx10; 
  out[4] += (6.708203932499369*f[6]*hamil[7]+7.500000000000001*f[2]*hamil[7]+3.0*f[1]*hamil[6]+7.5*f[3]*hamil[5]+3.0*hamil[3]*f[4]+3.3541019662496847*f[0]*hamil[3]+3.3541019662496847*f[1]*hamil[2])*dv10*dx10; 
  out[5] += (-(3.0*f[2]*hamil[7])-6.708203932499369*hamil[6]*f[7]-7.500000000000001*f[1]*hamil[6]-3.0*hamil[3]*f[5]-7.5*f[3]*hamil[4]-3.3541019662496847*f[0]*hamil[3]-3.3541019662496847*hamil[1]*f[2])*dv10*dx10; 
  out[6] += (6.708203932499369*f[5]*hamil[7]+6.708203932499369*f[4]*hamil[7]+7.5*f[0]*hamil[7]+6.708203932499369*hamil[5]*f[7]+1.5*hamil[3]*f[6]+7.500000000000001*f[1]*hamil[5]-3.0*f[1]*hamil[4]-1.5*hamil[1]*f[4]+3.3541019662496843*f[2]*hamil[3]+3.3541019662496843*hamil[2]*f[3])*dv10*dx10; 
  out[7] += (-(1.5*hamil[3]*f[7])-6.708203932499369*f[5]*hamil[6]-6.708203932499369*f[4]*hamil[6]-7.5*f[0]*hamil[6]-6.708203932499369*hamil[4]*f[6]+3.0*f[2]*hamil[5]+1.5*hamil[2]*f[5]-7.500000000000001*f[2]*hamil[4]-3.3541019662496843*f[1]*hamil[3]-3.3541019662496843*hamil[1]*f[3])*dv10*dx10; 

} 
