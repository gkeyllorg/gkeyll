#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_hamil_vel_sparse_vol_2x2v_ser_p2(const double *w, const double *dxv, const double *vmap,
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
  out[1] += (3.3541019662496847*f[3]*hamil[4]+1.5*f[0]*hamil[1])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[5] += (3.3541019662496847*hamil[4]*f[7]+1.5*hamil[1]*f[2])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[6] += (3.0*hamil[4]*f[13]+3.3541019662496847*f[0]*hamil[4]+1.5*hamil[1]*f[3])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[8] += (3.3541019662496847*hamil[4]*f[10]+1.5*hamil[1]*f[4])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[11] += (7.5*hamil[4]*f[6]+3.3541019662496847*f[1]*hamil[1])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[15] += (3.0*hamil[4]*f[24]+1.5*hamil[1]*f[7]+3.3541019662496847*f[2]*hamil[4])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[16] += (3.3541019662496847*hamil[4]*f[18]+1.5*hamil[1]*f[9])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[17] += (3.0*hamil[4]*f[27]+1.5*hamil[1]*f[10]+3.3541019662496847*f[4]*hamil[4])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[19] += (7.500000000000001*hamil[4]*f[15]+3.3541019662496843*hamil[1]*f[5])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[20] += (3.3541019662496847*hamil[4]*f[22]+1.5*hamil[1]*f[12])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[21] += (6.708203932499369*hamil[4]*f[23]+3.3541019662496843*hamil[1]*f[6]+7.500000000000001*f[1]*hamil[4])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[23] += (1.5*hamil[1]*f[13]+3.0*f[3]*hamil[4])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[25] += (7.500000000000001*hamil[4]*f[17]+3.3541019662496843*hamil[1]*f[8])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[28] += (3.3541019662496847*hamil[4]*f[30]+1.5*hamil[1]*f[14])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[31] += (3.0*hamil[4]*f[40]+1.5*hamil[1]*f[18]+3.3541019662496847*hamil[4]*f[9])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[32] += (6.708203932499369*hamil[4]*f[34]+3.3541019662496847*hamil[1]*f[15]+7.5*hamil[4]*f[5])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[33] += (1.5*hamil[1]*f[22]+3.3541019662496847*hamil[4]*f[12])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[34] += (1.5*hamil[1]*f[24]+3.0*hamil[4]*f[7])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[35] += (7.5*hamil[4]*f[31]+3.3541019662496847*hamil[1]*f[16])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[36] += (3.3541019662496847*hamil[4]*f[38]+1.5*hamil[1]*f[26])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[37] += (6.708203932499369*hamil[4]*f[39]+3.3541019662496847*hamil[1]*f[17]+7.5*hamil[4]*f[8])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[39] += (1.5*hamil[1]*f[27]+3.0*hamil[4]*f[10])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[41] += (3.3541019662496847*hamil[4]*f[43]+1.5*hamil[1]*f[29])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[42] += (1.5*hamil[1]*f[30]+3.3541019662496847*hamil[4]*f[14])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[44] += (6.708203932499369*hamil[4]*f[46]+3.3541019662496843*hamil[1]*f[31]+7.500000000000001*hamil[4]*f[16])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[45] += (1.5*hamil[1]*f[38]+3.3541019662496847*hamil[4]*f[26])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[46] += (1.5*hamil[1]*f[40]+3.0*hamil[4]*f[18])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[47] += (1.5*hamil[1]*f[43]+3.3541019662496847*hamil[4]*f[29])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 

  double dx11 = 2.0/dxv[1]; 
  double dv11 = 2.0/dxv[3]; 
  const double *jacob_vy = &jacob_vel[3]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  const double *jacob_cy = &jacob_pos[3]; 
  const double jacob_cy_inv = 1.0/jacob_cy[0]; 
  out[2] += (3.3541019662496847*f[4]*hamil[5]+1.5*f[0]*hamil[2])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[5] += (3.3541019662496847*hamil[5]*f[8]+1.5*f[1]*hamil[2])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[7] += (3.3541019662496847*hamil[5]*f[10]+1.5*hamil[2]*f[3])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[9] += (3.0*hamil[5]*f[14]+3.3541019662496847*f[0]*hamil[5]+1.5*hamil[2]*f[4])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[12] += (7.5*hamil[5]*f[9]+3.3541019662496847*f[2]*hamil[2])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[15] += (3.3541019662496847*hamil[5]*f[17]+1.5*hamil[2]*f[6])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[16] += (3.0*hamil[5]*f[28]+1.5*hamil[2]*f[8]+3.3541019662496847*f[1]*hamil[5])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[18] += (3.0*hamil[5]*f[30]+1.5*hamil[2]*f[10]+3.3541019662496847*f[3]*hamil[5])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[19] += (3.3541019662496847*hamil[5]*f[25]+1.5*hamil[2]*f[11])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[20] += (7.500000000000001*hamil[5]*f[16]+3.3541019662496843*hamil[2]*f[5])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[22] += (7.500000000000001*hamil[5]*f[18]+3.3541019662496843*hamil[2]*f[7])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[24] += (3.3541019662496847*hamil[5]*f[27]+1.5*hamil[2]*f[13])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[26] += (6.708203932499369*hamil[5]*f[29]+3.3541019662496843*hamil[2]*f[9]+7.500000000000001*f[2]*hamil[5])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[29] += (1.5*hamil[2]*f[14]+3.0*f[4]*hamil[5])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[31] += (3.0*hamil[5]*f[42]+1.5*hamil[2]*f[17]+3.3541019662496847*hamil[5]*f[6])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[32] += (3.3541019662496847*hamil[5]*f[37]+1.5*hamil[2]*f[21])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[33] += (7.5*hamil[5]*f[31]+3.3541019662496847*hamil[2]*f[15])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[34] += (3.3541019662496847*hamil[5]*f[39]+1.5*hamil[2]*f[23])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[35] += (1.5*hamil[2]*f[25]+3.3541019662496847*hamil[5]*f[11])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[36] += (6.708203932499369*hamil[5]*f[41]+3.3541019662496847*hamil[2]*f[16]+7.5*f[5]*hamil[5])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[38] += (6.708203932499369*hamil[5]*f[43]+3.3541019662496847*hamil[2]*f[18]+7.5*hamil[5]*f[7])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[40] += (1.5*hamil[2]*f[27]+3.3541019662496847*hamil[5]*f[13])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[41] += (1.5*hamil[2]*f[28]+3.0*hamil[5]*f[8])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[43] += (1.5*hamil[2]*f[30]+3.0*hamil[5]*f[10])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[44] += (1.5*hamil[2]*f[37]+3.3541019662496847*hamil[5]*f[21])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[45] += (6.708203932499369*hamil[5]*f[47]+3.3541019662496843*hamil[2]*f[31]+7.500000000000001*hamil[5]*f[15])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[46] += (1.5*hamil[2]*f[39]+3.3541019662496847*hamil[5]*f[23])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[47] += (1.5*hamil[2]*f[42]+3.0*hamil[5]*f[17])*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 

} 
