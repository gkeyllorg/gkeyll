#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_1x1v_ser_p3(const double *w, const double *dxv, const double *vmap,
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
  out[1] += (7.24568837309472*hamil[3]*f[5]+3.24037034920393*f[0]*hamil[3]+4.743416490252569*f[2]*hamil[2]+2.1213203435596424*f[0]*hamil[1])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[3] += (6.363961030678928*hamil[3]*f[9]+4.242640687119286*hamil[2]*f[5]+9.721111047611789*f[2]*hamil[3]+4.743416490252569*f[0]*hamil[2]+2.1213203435596424*hamil[1]*f[2])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[4] += (16.201851746019653*hamil[3]*f[7]+7.24568837309472*f[1]*hamil[3]+10.606601717798211*hamil[2]*f[3]+4.743416490252569*f[1]*hamil[1])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[6] += (14.230249470757704*hamil[3]*f[11]+9.48683298050514*hamil[2]*f[7]+21.73706511928416*f[3]*hamil[3]+4.743416490252569*hamil[1]*f[3]+10.606601717798213*f[1]*hamil[2])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[7] += (4.16619044897648*hamil[2]*f[9]+7.869470848066687*hamil[3]*f[5]+2.1213203435596424*hamil[1]*f[5]+7.245688373094719*f[0]*hamil[3]+4.242640687119286*f[2]*hamil[2])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[8] += (16.201851746019653*hamil[2]*f[6]+11.067971810589325*hamil[3]*f[5]+11.067971810589325*hamil[3]*f[4]+7.24568837309472*hamil[1]*f[4]+4.949747468305832*f[0]*hamil[3]+7.24568837309472*f[2]*hamil[2]+3.24037034920393*f[0]*hamil[1])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[10] += (9.721111047611789*hamil[3]*f[9]+33.20391543176798*hamil[3]*f[6]+7.245688373094719*hamil[1]*f[6]+6.48074069840786*hamil[2]*f[5]+16.20185174601965*hamil[2]*f[4]+14.849242404917497*f[2]*hamil[3]+7.245688373094718*f[0]*hamil[2]+3.2403703492039297*hamil[1]*f[2])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[11] += (7.560864148142503*hamil[3]*f[9]+2.121320343559642*hamil[1]*f[9]+4.166190448976481*hamil[2]*f[5]+6.363961030678925*f[2]*hamil[3])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 

} 
