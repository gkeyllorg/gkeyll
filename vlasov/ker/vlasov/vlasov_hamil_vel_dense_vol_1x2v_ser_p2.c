#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *vmap,
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
  out[1] += (1.5*hamil[7]*f[9]+3.3541019662496843*f[6]*hamil[6]+3.3541019662496847*f[2]*hamil[4]+1.5*f[3]*hamil[3]+1.5*f[0]*hamil[1])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[4] += (1.5*hamil[7]*f[16]+3.0*hamil[6]*f[14]+3.0*hamil[4]*f[8]+3.3541019662496843*f[3]*hamil[6]+1.5*hamil[3]*f[6]+3.3541019662496847*f[0]*hamil[4]+1.5*hamil[1]*f[2])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[5] += (3.0*hamil[6]*f[16]+1.3416407864998738*hamil[3]*f[9]+1.3416407864998738*f[3]*hamil[7]+3.3541019662496843*f[2]*hamil[6]+3.3541019662496847*hamil[4]*f[6]+1.5*f[0]*hamil[3]+1.5*hamil[1]*f[3])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[7] += (3.3541019662496847*hamil[7]*f[15]+7.500000000000001*hamil[6]*f[10]+3.3541019662496847*hamil[3]*f[5]+7.5*f[4]*hamil[4]+3.3541019662496847*f[1]*hamil[1])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[10] += (1.3416407864998738*hamil[3]*f[16]+3.0*hamil[4]*f[14]+3.0*hamil[6]*f[9]+3.0*hamil[6]*f[8]+1.3416407864998738*f[6]*hamil[7]+3.3541019662496843*f[0]*hamil[6]+1.5*hamil[1]*f[6]+3.3541019662496847*f[3]*hamil[4]+1.5*f[2]*hamil[3])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[11] += (3.3541019662496847*hamil[7]*f[19]+6.708203932499369*hamil[6]*f[18]+6.708203932499369*hamil[4]*f[12]+3.3541019662496843*hamil[3]*f[10]+7.5*f[5]*hamil[6]+7.500000000000001*f[1]*hamil[4]+3.3541019662496843*hamil[1]*f[4])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[12] += (1.5*hamil[3]*f[14]+1.5*hamil[1]*f[8]+3.0*f[6]*hamil[6]+3.0*f[2]*hamil[4])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[13] += (6.708203932499369*hamil[6]*f[19]+3.0*hamil[3]*f[15]+7.500000000000001*hamil[4]*f[10]+3.0*f[5]*hamil[7]+7.5*f[4]*hamil[6]+3.3541019662496843*hamil[1]*f[5]+3.3541019662496843*f[1]*hamil[3])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[15] += (3.3541019662496847*hamil[4]*f[16]+0.9583148474999099*hamil[7]*f[9]+1.5*hamil[1]*f[9]+1.5*f[0]*hamil[7]+3.0*f[6]*hamil[6]+1.3416407864998738*f[3]*hamil[3])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[17] += (3.0*hamil[3]*f[19]+6.708203932499369*hamil[4]*f[18]+6.708203932499369*hamil[6]*f[15]+6.708203932499369*hamil[6]*f[12]+3.0*hamil[7]*f[10]+3.3541019662496847*hamil[1]*f[10]+7.500000000000001*f[1]*hamil[6]+7.5*hamil[4]*f[5]+3.3541019662496847*hamil[3]*f[4])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[18] += (2.6832815729997477*hamil[6]*f[16]+1.3416407864998738*hamil[7]*f[14]+1.5*hamil[1]*f[14]+1.5*hamil[3]*f[8]+3.0*f[2]*hamil[6]+3.0*hamil[4]*f[6])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[19] += (0.9583148474999099*hamil[7]*f[16]+1.5*hamil[1]*f[16]+2.6832815729997477*hamil[6]*f[14]+3.3541019662496847*hamil[4]*f[9]+1.5*f[2]*hamil[7]+3.0*f[3]*hamil[6]+1.3416407864998738*hamil[3]*f[6])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 

} 
