#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_hamil_vel_dense_vol_2x3v_ser_p1(const double *w, const double *dxv, const double *vmap,
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
  out[1] += (1.060660171779821*hamil[7]*f[15]+1.060660171779821*f[5]*hamil[5]+1.060660171779821*f[4]*hamil[4]+1.060660171779821*f[0]*hamil[1])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[6] += (1.060660171779821*hamil[7]*f[24]+1.060660171779821*hamil[5]*f[13]+1.060660171779821*hamil[4]*f[10]+1.060660171779821*hamil[1]*f[2])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[7] += (1.060660171779821*hamil[7]*f[25]+1.060660171779821*hamil[5]*f[14]+1.060660171779821*hamil[4]*f[11]+1.060660171779821*hamil[1]*f[3])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[9] += (1.060660171779821*hamil[5]*f[15]+1.060660171779821*f[5]*hamil[7]+1.060660171779821*f[0]*hamil[4]+1.060660171779821*hamil[1]*f[4])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[12] += (1.060660171779821*hamil[4]*f[15]+1.060660171779821*f[4]*hamil[7]+1.060660171779821*f[0]*hamil[5]+1.060660171779821*hamil[1]*f[5])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[16] += (1.060660171779821*hamil[7]*f[30]+1.060660171779821*hamil[5]*f[22]+1.060660171779821*hamil[4]*f[19]+1.060660171779821*hamil[1]*f[8])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[17] += (1.060660171779821*hamil[5]*f[24]+1.060660171779821*hamil[7]*f[13]+1.060660171779821*hamil[1]*f[10]+1.060660171779821*f[2]*hamil[4])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[18] += (1.060660171779821*hamil[5]*f[25]+1.060660171779821*hamil[7]*f[14]+1.060660171779821*hamil[1]*f[11]+1.060660171779821*f[3]*hamil[4])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[20] += (1.060660171779821*hamil[4]*f[24]+1.060660171779821*hamil[1]*f[13]+1.060660171779821*hamil[7]*f[10]+1.060660171779821*f[2]*hamil[5])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[21] += (1.060660171779821*hamil[4]*f[25]+1.060660171779821*hamil[1]*f[14]+1.060660171779821*hamil[7]*f[11]+1.060660171779821*f[3]*hamil[5])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[23] += (1.060660171779821*hamil[1]*f[15]+1.060660171779821*f[0]*hamil[7]+1.060660171779821*f[4]*hamil[5]+1.060660171779821*hamil[4]*f[5])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[26] += (1.060660171779821*hamil[5]*f[30]+1.060660171779821*hamil[7]*f[22]+1.060660171779821*hamil[1]*f[19]+1.060660171779821*hamil[4]*f[8])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[27] += (1.060660171779821*hamil[4]*f[30]+1.060660171779821*hamil[1]*f[22]+1.060660171779821*hamil[7]*f[19]+1.060660171779821*hamil[5]*f[8])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[28] += (1.060660171779821*hamil[1]*f[24]+1.060660171779821*hamil[4]*f[13]+1.060660171779821*hamil[5]*f[10]+1.060660171779821*f[2]*hamil[7])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[29] += (1.060660171779821*hamil[1]*f[25]+1.060660171779821*hamil[4]*f[14]+1.060660171779821*hamil[5]*f[11]+1.060660171779821*f[3]*hamil[7])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[31] += (1.060660171779821*hamil[1]*f[30]+1.060660171779821*hamil[4]*f[22]+1.060660171779821*hamil[5]*f[19]+1.060660171779821*hamil[7]*f[8])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 

  double dx11 = 2.0/dxv[1]; 
  double dv11 = 2.0/dxv[3]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  const double *jacob_cy = &jacob_pos[2]; 
  const double jacob_cy_inv = 1.0/jacob_cy[0]; 
  out[2] += (1.060660171779821*hamil[7]*f[14]+1.060660171779821*f[5]*hamil[6]+1.060660171779821*f[3]*hamil[4]+1.060660171779821*f[0]*hamil[2])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[6] += (1.060660171779821*hamil[7]*f[21]+1.060660171779821*hamil[6]*f[12]+1.060660171779821*hamil[4]*f[7]+1.060660171779821*f[1]*hamil[2])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[8] += (1.060660171779821*hamil[6]*f[14]+1.060660171779821*f[5]*hamil[7]+1.060660171779821*f[0]*hamil[4]+1.060660171779821*hamil[2]*f[3])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[10] += (1.060660171779821*hamil[7]*f[25]+1.060660171779821*hamil[6]*f[15]+1.060660171779821*hamil[4]*f[11]+1.060660171779821*hamil[2]*f[4])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[13] += (1.060660171779821*hamil[4]*f[14]+1.060660171779821*f[3]*hamil[7]+1.060660171779821*f[0]*hamil[6]+1.060660171779821*hamil[2]*f[5])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[16] += (1.060660171779821*hamil[6]*f[21]+1.060660171779821*hamil[7]*f[12]+1.060660171779821*hamil[2]*f[7]+1.060660171779821*f[1]*hamil[4])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[17] += (1.060660171779821*hamil[7]*f[29]+1.060660171779821*hamil[6]*f[23]+1.060660171779821*hamil[4]*f[18]+1.060660171779821*hamil[2]*f[9])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[19] += (1.060660171779821*hamil[6]*f[25]+1.060660171779821*hamil[7]*f[15]+1.060660171779821*hamil[2]*f[11]+1.060660171779821*f[4]*hamil[4])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[20] += (1.060660171779821*hamil[4]*f[21]+1.060660171779821*hamil[2]*f[12]+1.060660171779821*f[7]*hamil[7]+1.060660171779821*f[1]*hamil[6])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[22] += (1.060660171779821*hamil[2]*f[14]+1.060660171779821*f[0]*hamil[7]+1.060660171779821*f[3]*hamil[6]+1.060660171779821*hamil[4]*f[5])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[24] += (1.060660171779821*hamil[4]*f[25]+1.060660171779821*hamil[2]*f[15]+1.060660171779821*hamil[7]*f[11]+1.060660171779821*f[4]*hamil[6])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[26] += (1.060660171779821*hamil[6]*f[29]+1.060660171779821*hamil[7]*f[23]+1.060660171779821*hamil[2]*f[18]+1.060660171779821*hamil[4]*f[9])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[27] += (1.060660171779821*hamil[2]*f[21]+1.060660171779821*hamil[4]*f[12]+1.060660171779821*f[1]*hamil[7]+1.060660171779821*hamil[6]*f[7])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[28] += (1.060660171779821*hamil[4]*f[29]+1.060660171779821*hamil[2]*f[23]+1.060660171779821*hamil[7]*f[18]+1.060660171779821*hamil[6]*f[9])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[30] += (1.060660171779821*hamil[2]*f[25]+1.060660171779821*hamil[4]*f[15]+1.060660171779821*hamil[6]*f[11]+1.060660171779821*f[4]*hamil[7])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[31] += (1.060660171779821*hamil[2]*f[29]+1.060660171779821*hamil[4]*f[23]+1.060660171779821*hamil[6]*f[18]+1.060660171779821*hamil[7]*f[9])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 

} 
