#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_1x3v_ser_p1(const double *w, const double *dxv, const double *vmap,
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
  out[1] += (1.060660171779821*hamil[7]*f[10]+1.060660171779821*f[4]*hamil[5]+1.060660171779821*f[3]*hamil[4]+1.060660171779821*f[0]*hamil[1])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[5] += (1.060660171779821*hamil[7]*f[14]+1.060660171779821*hamil[5]*f[9]+1.060660171779821*hamil[4]*f[7]+1.060660171779821*hamil[1]*f[2])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[6] += (1.060660171779821*hamil[5]*f[10]+1.060660171779821*f[4]*hamil[7]+1.060660171779821*f[0]*hamil[4]+1.060660171779821*hamil[1]*f[3])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[8] += (1.060660171779821*hamil[4]*f[10]+1.060660171779821*f[3]*hamil[7]+1.060660171779821*f[0]*hamil[5]+1.060660171779821*hamil[1]*f[4])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[11] += (1.060660171779821*hamil[5]*f[14]+1.060660171779821*hamil[7]*f[9]+1.060660171779821*hamil[1]*f[7]+1.060660171779821*f[2]*hamil[4])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[12] += (1.060660171779821*hamil[4]*f[14]+1.060660171779821*hamil[1]*f[9]+1.060660171779821*f[7]*hamil[7]+1.060660171779821*f[2]*hamil[5])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[13] += (1.060660171779821*hamil[1]*f[10]+1.060660171779821*f[0]*hamil[7]+1.060660171779821*f[3]*hamil[5]+1.060660171779821*f[4]*hamil[4])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[15] += (1.060660171779821*hamil[1]*f[14]+1.060660171779821*hamil[4]*f[9]+1.060660171779821*f[2]*hamil[7]+1.060660171779821*hamil[5]*f[7])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 

} 
