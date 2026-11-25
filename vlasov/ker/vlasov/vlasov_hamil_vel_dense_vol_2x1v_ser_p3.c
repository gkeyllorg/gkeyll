#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_2x1v_ser_p3(const double *w, const double *dxv, const double *vmap,
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
  out[1] += (7.24568837309472*hamil[3]*f[9]+3.24037034920393*f[0]*hamil[3]+4.743416490252569*hamil[2]*f[3]+2.1213203435596424*f[0]*hamil[1])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[4] += (7.245688373094719*hamil[3]*f[16]+4.743416490252569*hamil[2]*f[6]+3.24037034920393*f[2]*hamil[3]+2.1213203435596424*hamil[1]*f[2])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[5] += (6.363961030678928*hamil[3]*f[19]+4.242640687119286*hamil[2]*f[9]+9.721111047611789*f[3]*hamil[3]+2.1213203435596424*hamil[1]*f[3]+4.743416490252569*f[0]*hamil[2])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[7] += (16.201851746019653*hamil[3]*f[15]+10.606601717798211*hamil[2]*f[5]+7.24568837309472*f[1]*hamil[3]+4.743416490252569*f[1]*hamil[1])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[10] += (6.363961030678925*hamil[3]*f[28]+4.242640687119286*hamil[2]*f[16]+9.721111047611789*hamil[3]*f[6]+2.1213203435596424*hamil[1]*f[6]+4.743416490252569*f[2]*hamil[2])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[11] += (16.201851746019653*hamil[3]*f[22]+10.606601717798213*hamil[2]*f[10]+7.245688373094719*hamil[3]*f[4]+4.743416490252569*hamil[1]*f[4])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[12] += (4.743416490252569*hamil[2]*f[14]+3.2403703492039297*hamil[3]*f[8]+2.1213203435596424*hamil[1]*f[8])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[13] += (14.230249470757704*hamil[3]*f[27]+9.48683298050514*hamil[2]*f[15]+21.73706511928416*hamil[3]*f[5]+4.743416490252569*hamil[1]*f[5]+10.606601717798213*f[1]*hamil[2])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[15] += (4.16619044897648*hamil[2]*f[19]+7.869470848066687*hamil[3]*f[9]+2.1213203435596424*hamil[1]*f[9]+7.245688373094719*f[0]*hamil[3]+4.242640687119286*hamil[2]*f[3])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[17] += (16.201851746019653*hamil[2]*f[13]+11.067971810589325*hamil[3]*f[9]+11.067971810589325*hamil[3]*f[7]+7.24568837309472*hamil[1]*f[7]+4.949747468305832*f[0]*hamil[3]+7.24568837309472*hamil[2]*f[3]+3.24037034920393*f[0]*hamil[1])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[20] += (14.230249470757707*hamil[3]*f[31]+9.48683298050514*hamil[2]*f[22]+21.737065119284154*hamil[3]*f[10]+4.743416490252569*hamil[1]*f[10]+10.606601717798211*hamil[2]*f[4])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[21] += (9.721111047611789*hamil[3]*f[14]+2.1213203435596424*hamil[1]*f[14]+4.743416490252569*hamil[2]*f[8])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[22] += (4.166190448976481*hamil[2]*f[28]+7.869470848066687*hamil[3]*f[16]+2.1213203435596424*hamil[1]*f[16]+4.242640687119286*hamil[2]*f[6]+7.24568837309472*f[2]*hamil[3])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[23] += (16.20185174601965*hamil[2]*f[20]+11.067971810589327*hamil[3]*f[16]+11.067971810589327*hamil[3]*f[11]+7.245688373094719*hamil[1]*f[11]+7.245688373094718*hamil[2]*f[6]+4.949747468305833*f[2]*hamil[3]+3.2403703492039297*hamil[1]*f[2])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[24] += (4.743416490252569*hamil[2]*f[26]+3.2403703492039297*hamil[3]*f[18]+2.121320343559642*hamil[1]*f[18])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[25] += (9.721111047611789*hamil[3]*f[19]+33.20391543176798*hamil[3]*f[13]+7.245688373094719*hamil[1]*f[13]+6.48074069840786*hamil[2]*f[9]+16.20185174601965*hamil[2]*f[7]+14.849242404917497*f[3]*hamil[3]+3.2403703492039297*hamil[1]*f[3]+7.245688373094718*f[0]*hamil[2])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[27] += (7.560864148142503*hamil[3]*f[19]+2.121320343559642*hamil[1]*f[19]+4.166190448976481*hamil[2]*f[9]+6.363961030678925*f[3]*hamil[3])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[29] += (9.721111047611789*hamil[3]*f[28]+33.20391543176798*hamil[3]*f[20]+7.24568837309472*hamil[1]*f[20]+6.48074069840786*hamil[2]*f[16]+16.201851746019653*hamil[2]*f[11]+14.849242404917497*hamil[3]*f[6]+3.24037034920393*hamil[1]*f[6]+7.24568837309472*f[2]*hamil[2])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[30] += (9.721111047611787*hamil[3]*f[26]+2.121320343559642*hamil[1]*f[26]+4.743416490252569*hamil[2]*f[18])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[31] += (7.560864148142503*hamil[3]*f[28]+2.121320343559642*hamil[1]*f[28]+4.16619044897648*hamil[2]*f[16]+6.363961030678928*hamil[3]*f[6])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 

} 
