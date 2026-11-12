#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_2x1v_ser_p2(const double *w, const double *dxv, const double *vmap,
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
  out[1] += (4.743416490252569*hamil[2]*f[3]+2.1213203435596424*f[0]*hamil[1])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[4] += (4.743416490252569*hamil[2]*f[6]+2.1213203435596424*hamil[1]*f[2])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[5] += (4.242640687119286*hamil[2]*f[9]+2.1213203435596424*hamil[1]*f[3]+4.743416490252569*f[0]*hamil[2])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[7] += (10.606601717798211*hamil[2]*f[5]+4.743416490252569*f[1]*hamil[1])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[10] += (4.242640687119286*hamil[2]*f[16]+2.1213203435596424*hamil[1]*f[6]+4.743416490252569*f[2]*hamil[2])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[11] += (10.606601717798213*hamil[2]*f[10]+4.743416490252569*hamil[1]*f[4])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[12] += (4.743416490252569*hamil[2]*f[14]+2.1213203435596424*hamil[1]*f[8])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[13] += (9.48683298050514*hamil[2]*f[15]+4.743416490252569*hamil[1]*f[5]+10.606601717798213*f[1]*hamil[2])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[15] += (2.1213203435596424*hamil[1]*f[9]+4.242640687119286*hamil[2]*f[3])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[17] += (9.48683298050514*hamil[2]*f[19]+4.743416490252569*hamil[1]*f[10]+10.606601717798211*hamil[2]*f[4])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[18] += (2.1213203435596424*hamil[1]*f[14]+4.743416490252569*hamil[2]*f[8])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[19] += (2.1213203435596424*hamil[1]*f[16]+4.242640687119286*hamil[2]*f[6])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 

} 
