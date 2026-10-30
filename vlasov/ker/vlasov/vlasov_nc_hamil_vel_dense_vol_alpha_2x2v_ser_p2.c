#include <gkyl_nc_hamil_gen_vol_comps_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_vel_dense_vol_alpha_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
   const double *vmap, const double *jacob_pos, const double *jacob_vel,
   const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT out) 
{ 
  const double dx10 = 2.0/dxv[0]; 
  const double dx11 = 2.0/dxv[1]; 
  const double dv0 = dxv[2]; 
  const double dv1 = dxv[3]; 
  const double dv10 = 2.0/dxv[2]; 
  const double dv11 = 2.0/dxv[3]; 
  const double wx1 = w[2]; 
  const double wx2 = w[3]; 
  const double jacob_cx_inv = 1.0/jacob_pos[0]; 
  const double jacob_cy_inv = 1.0/jacob_pos[3]; 
  const double jacob_vx_inv = 1.0/jacob_vel[0]; 
  const double jacob_vy_inv = 1.0/jacob_vel[3]; 
  const double *vmap_vx = &vmap[0]; 
  const double *vmap_vy = &vmap[4]; 
  
  const double *pt_conf_0 = &poisson_tensor_conf[0]; 
  if (dir == 0) { 
  out[0] += 1.7320508075688772*pt_conf_0[0]*hamil[1]*dv10*jacob_vx_inv; 
  out[1] += 1.7320508075688772*hamil[1]*pt_conf_0[1]*dv10*jacob_vx_inv; 
  out[2] += 1.7320508075688772*hamil[1]*pt_conf_0[2]*dv10*jacob_vx_inv; 
  out[3] += 3.872983346207417*pt_conf_0[0]*hamil[4]*dv10*jacob_vx_inv; 
  out[4] += 1.7320508075688772*pt_conf_0[0]*hamil[3]*dv10*jacob_vx_inv; 
  out[5] += 1.7320508075688772*hamil[1]*pt_conf_0[3]*dv10*jacob_vx_inv; 
  out[6] += 3.872983346207417*pt_conf_0[1]*hamil[4]*dv10*jacob_vx_inv; 
  out[7] += 3.872983346207417*pt_conf_0[2]*hamil[4]*dv10*jacob_vx_inv; 
  out[8] += 1.7320508075688772*pt_conf_0[1]*hamil[3]*dv10*jacob_vx_inv; 
  out[9] += 1.7320508075688772*pt_conf_0[2]*hamil[3]*dv10*jacob_vx_inv; 
  out[10] += 3.872983346207417*pt_conf_0[0]*hamil[6]*dv10*jacob_vx_inv; 
  out[11] += 1.7320508075688772*hamil[1]*pt_conf_0[4]*dv10*jacob_vx_inv; 
  out[12] += 1.7320508075688772*hamil[1]*pt_conf_0[5]*dv10*jacob_vx_inv; 
  out[14] += 1.7320508075688774*pt_conf_0[0]*hamil[7]*dv10*jacob_vx_inv; 
  out[15] += 3.872983346207417*pt_conf_0[3]*hamil[4]*dv10*jacob_vx_inv; 
  out[16] += 1.7320508075688772*hamil[3]*pt_conf_0[3]*dv10*jacob_vx_inv; 
  out[17] += 3.872983346207417*pt_conf_0[1]*hamil[6]*dv10*jacob_vx_inv; 
  out[18] += 3.872983346207417*pt_conf_0[2]*hamil[6]*dv10*jacob_vx_inv; 
  out[19] += 1.7320508075688772*hamil[1]*pt_conf_0[6]*dv10*jacob_vx_inv; 
  out[20] += 1.7320508075688772*hamil[1]*pt_conf_0[7]*dv10*jacob_vx_inv; 
  out[21] += 3.872983346207417*hamil[4]*pt_conf_0[4]*dv10*jacob_vx_inv; 
  out[22] += 3.872983346207417*hamil[4]*pt_conf_0[5]*dv10*jacob_vx_inv; 
  out[25] += 1.7320508075688774*hamil[3]*pt_conf_0[4]*dv10*jacob_vx_inv; 
  out[26] += 1.7320508075688774*hamil[3]*pt_conf_0[5]*dv10*jacob_vx_inv; 
  out[28] += 1.7320508075688772*pt_conf_0[1]*hamil[7]*dv10*jacob_vx_inv; 
  out[29] += 1.7320508075688772*pt_conf_0[2]*hamil[7]*dv10*jacob_vx_inv; 
  out[31] += 3.872983346207417*pt_conf_0[3]*hamil[6]*dv10*jacob_vx_inv; 
  out[32] += 3.872983346207417*hamil[4]*pt_conf_0[6]*dv10*jacob_vx_inv; 
  out[33] += 3.872983346207417*hamil[4]*pt_conf_0[7]*dv10*jacob_vx_inv; 
  out[35] += 1.7320508075688774*hamil[3]*pt_conf_0[6]*dv10*jacob_vx_inv; 
  out[36] += 1.7320508075688774*hamil[3]*pt_conf_0[7]*dv10*jacob_vx_inv; 
  out[37] += 3.872983346207417*pt_conf_0[4]*hamil[6]*dv10*jacob_vx_inv; 
  out[38] += 3.872983346207417*pt_conf_0[5]*hamil[6]*dv10*jacob_vx_inv; 
  out[41] += 1.7320508075688774*pt_conf_0[3]*hamil[7]*dv10*jacob_vx_inv; 
  out[44] += 3.872983346207417*hamil[6]*pt_conf_0[6]*dv10*jacob_vx_inv; 
  out[45] += 3.872983346207417*hamil[6]*pt_conf_0[7]*dv10*jacob_vx_inv; 
  } 
  if (dir == 2) { 
  } 

  const double *pt_conf_1 = &poisson_tensor_conf[8]; 
  if (dir == 0) { 
  out[0] += 1.7320508075688772*pt_conf_1[0]*hamil[2]*dv11*jacob_vy_inv; 
  out[1] += 1.7320508075688772*pt_conf_1[1]*hamil[2]*dv11*jacob_vy_inv; 
  out[2] += 1.7320508075688772*hamil[2]*pt_conf_1[2]*dv11*jacob_vy_inv; 
  out[3] += 1.7320508075688772*pt_conf_1[0]*hamil[3]*dv11*jacob_vy_inv; 
  out[4] += 3.872983346207417*pt_conf_1[0]*hamil[5]*dv11*jacob_vy_inv; 
  out[5] += 1.7320508075688772*hamil[2]*pt_conf_1[3]*dv11*jacob_vy_inv; 
  out[6] += 1.7320508075688772*pt_conf_1[1]*hamil[3]*dv11*jacob_vy_inv; 
  out[7] += 1.7320508075688772*pt_conf_1[2]*hamil[3]*dv11*jacob_vy_inv; 
  out[8] += 3.872983346207417*pt_conf_1[1]*hamil[5]*dv11*jacob_vy_inv; 
  out[9] += 3.872983346207417*pt_conf_1[2]*hamil[5]*dv11*jacob_vy_inv; 
  out[10] += 3.872983346207417*pt_conf_1[0]*hamil[7]*dv11*jacob_vy_inv; 
  out[11] += 1.7320508075688772*hamil[2]*pt_conf_1[4]*dv11*jacob_vy_inv; 
  out[12] += 1.7320508075688772*hamil[2]*pt_conf_1[5]*dv11*jacob_vy_inv; 
  out[13] += 1.7320508075688772*pt_conf_1[0]*hamil[6]*dv11*jacob_vy_inv; 
  out[15] += 1.7320508075688772*hamil[3]*pt_conf_1[3]*dv11*jacob_vy_inv; 
  out[16] += 3.872983346207417*pt_conf_1[3]*hamil[5]*dv11*jacob_vy_inv; 
  out[17] += 3.872983346207417*pt_conf_1[1]*hamil[7]*dv11*jacob_vy_inv; 
  out[18] += 3.872983346207417*pt_conf_1[2]*hamil[7]*dv11*jacob_vy_inv; 
  out[19] += 1.7320508075688772*hamil[2]*pt_conf_1[6]*dv11*jacob_vy_inv; 
  out[20] += 1.7320508075688772*hamil[2]*pt_conf_1[7]*dv11*jacob_vy_inv; 
  out[21] += 1.7320508075688774*hamil[3]*pt_conf_1[4]*dv11*jacob_vy_inv; 
  out[22] += 1.7320508075688774*hamil[3]*pt_conf_1[5]*dv11*jacob_vy_inv; 
  out[23] += 1.7320508075688772*pt_conf_1[1]*hamil[6]*dv11*jacob_vy_inv; 
  out[24] += 1.7320508075688772*pt_conf_1[2]*hamil[6]*dv11*jacob_vy_inv; 
  out[25] += 3.872983346207417*pt_conf_1[4]*hamil[5]*dv11*jacob_vy_inv; 
  out[26] += 3.872983346207417*hamil[5]*pt_conf_1[5]*dv11*jacob_vy_inv; 
  out[31] += 3.872983346207417*pt_conf_1[3]*hamil[7]*dv11*jacob_vy_inv; 
  out[32] += 1.7320508075688774*hamil[3]*pt_conf_1[6]*dv11*jacob_vy_inv; 
  out[33] += 1.7320508075688774*hamil[3]*pt_conf_1[7]*dv11*jacob_vy_inv; 
  out[34] += 1.7320508075688772*pt_conf_1[3]*hamil[6]*dv11*jacob_vy_inv; 
  out[35] += 3.872983346207417*hamil[5]*pt_conf_1[6]*dv11*jacob_vy_inv; 
  out[36] += 3.872983346207417*hamil[5]*pt_conf_1[7]*dv11*jacob_vy_inv; 
  out[37] += 3.872983346207417*pt_conf_1[4]*hamil[7]*dv11*jacob_vy_inv; 
  out[38] += 3.872983346207417*pt_conf_1[5]*hamil[7]*dv11*jacob_vy_inv; 
  out[44] += 3.872983346207417*pt_conf_1[6]*hamil[7]*dv11*jacob_vy_inv; 
  out[45] += 3.872983346207417*hamil[7]*pt_conf_1[7]*dv11*jacob_vy_inv; 
  } 
  if (dir == 3) { 
  } 

  const double *pt_conf_2 = &poisson_tensor_conf[16]; 
  if (dir == 1) { 
  out[0] += 1.7320508075688772*pt_conf_2[0]*hamil[1]*dv10*jacob_vx_inv; 
  out[1] += 1.7320508075688772*hamil[1]*pt_conf_2[1]*dv10*jacob_vx_inv; 
  out[2] += 1.7320508075688772*hamil[1]*pt_conf_2[2]*dv10*jacob_vx_inv; 
  out[3] += 3.872983346207417*pt_conf_2[0]*hamil[4]*dv10*jacob_vx_inv; 
  out[4] += 1.7320508075688772*pt_conf_2[0]*hamil[3]*dv10*jacob_vx_inv; 
  out[5] += 1.7320508075688772*hamil[1]*pt_conf_2[3]*dv10*jacob_vx_inv; 
  out[6] += 3.872983346207417*pt_conf_2[1]*hamil[4]*dv10*jacob_vx_inv; 
  out[7] += 3.872983346207417*pt_conf_2[2]*hamil[4]*dv10*jacob_vx_inv; 
  out[8] += 1.7320508075688772*pt_conf_2[1]*hamil[3]*dv10*jacob_vx_inv; 
  out[9] += 1.7320508075688772*pt_conf_2[2]*hamil[3]*dv10*jacob_vx_inv; 
  out[10] += 3.872983346207417*pt_conf_2[0]*hamil[6]*dv10*jacob_vx_inv; 
  out[11] += 1.7320508075688772*hamil[1]*pt_conf_2[4]*dv10*jacob_vx_inv; 
  out[12] += 1.7320508075688772*hamil[1]*pt_conf_2[5]*dv10*jacob_vx_inv; 
  out[14] += 1.7320508075688774*pt_conf_2[0]*hamil[7]*dv10*jacob_vx_inv; 
  out[15] += 3.872983346207417*pt_conf_2[3]*hamil[4]*dv10*jacob_vx_inv; 
  out[16] += 1.7320508075688772*hamil[3]*pt_conf_2[3]*dv10*jacob_vx_inv; 
  out[17] += 3.872983346207417*pt_conf_2[1]*hamil[6]*dv10*jacob_vx_inv; 
  out[18] += 3.872983346207417*pt_conf_2[2]*hamil[6]*dv10*jacob_vx_inv; 
  out[19] += 1.7320508075688772*hamil[1]*pt_conf_2[6]*dv10*jacob_vx_inv; 
  out[20] += 1.7320508075688772*hamil[1]*pt_conf_2[7]*dv10*jacob_vx_inv; 
  out[21] += 3.872983346207417*hamil[4]*pt_conf_2[4]*dv10*jacob_vx_inv; 
  out[22] += 3.872983346207417*hamil[4]*pt_conf_2[5]*dv10*jacob_vx_inv; 
  out[25] += 1.7320508075688774*hamil[3]*pt_conf_2[4]*dv10*jacob_vx_inv; 
  out[26] += 1.7320508075688774*hamil[3]*pt_conf_2[5]*dv10*jacob_vx_inv; 
  out[28] += 1.7320508075688772*pt_conf_2[1]*hamil[7]*dv10*jacob_vx_inv; 
  out[29] += 1.7320508075688772*pt_conf_2[2]*hamil[7]*dv10*jacob_vx_inv; 
  out[31] += 3.872983346207417*pt_conf_2[3]*hamil[6]*dv10*jacob_vx_inv; 
  out[32] += 3.872983346207417*hamil[4]*pt_conf_2[6]*dv10*jacob_vx_inv; 
  out[33] += 3.872983346207417*hamil[4]*pt_conf_2[7]*dv10*jacob_vx_inv; 
  out[35] += 1.7320508075688774*hamil[3]*pt_conf_2[6]*dv10*jacob_vx_inv; 
  out[36] += 1.7320508075688774*hamil[3]*pt_conf_2[7]*dv10*jacob_vx_inv; 
  out[37] += 3.872983346207417*pt_conf_2[4]*hamil[6]*dv10*jacob_vx_inv; 
  out[38] += 3.872983346207417*pt_conf_2[5]*hamil[6]*dv10*jacob_vx_inv; 
  out[41] += 1.7320508075688774*pt_conf_2[3]*hamil[7]*dv10*jacob_vx_inv; 
  out[44] += 3.872983346207417*hamil[6]*pt_conf_2[6]*dv10*jacob_vx_inv; 
  out[45] += 3.872983346207417*hamil[6]*pt_conf_2[7]*dv10*jacob_vx_inv; 
  } 
  if (dir == 2) { 
  } 

  const double *pt_conf_3 = &poisson_tensor_conf[24]; 
  if (dir == 1) { 
  out[0] += 1.7320508075688772*pt_conf_3[0]*hamil[2]*dv11*jacob_vy_inv; 
  out[1] += 1.7320508075688772*pt_conf_3[1]*hamil[2]*dv11*jacob_vy_inv; 
  out[2] += 1.7320508075688772*hamil[2]*pt_conf_3[2]*dv11*jacob_vy_inv; 
  out[3] += 1.7320508075688772*pt_conf_3[0]*hamil[3]*dv11*jacob_vy_inv; 
  out[4] += 3.872983346207417*pt_conf_3[0]*hamil[5]*dv11*jacob_vy_inv; 
  out[5] += 1.7320508075688772*hamil[2]*pt_conf_3[3]*dv11*jacob_vy_inv; 
  out[6] += 1.7320508075688772*pt_conf_3[1]*hamil[3]*dv11*jacob_vy_inv; 
  out[7] += 1.7320508075688772*pt_conf_3[2]*hamil[3]*dv11*jacob_vy_inv; 
  out[8] += 3.872983346207417*pt_conf_3[1]*hamil[5]*dv11*jacob_vy_inv; 
  out[9] += 3.872983346207417*pt_conf_3[2]*hamil[5]*dv11*jacob_vy_inv; 
  out[10] += 3.872983346207417*pt_conf_3[0]*hamil[7]*dv11*jacob_vy_inv; 
  out[11] += 1.7320508075688772*hamil[2]*pt_conf_3[4]*dv11*jacob_vy_inv; 
  out[12] += 1.7320508075688772*hamil[2]*pt_conf_3[5]*dv11*jacob_vy_inv; 
  out[13] += 1.7320508075688772*pt_conf_3[0]*hamil[6]*dv11*jacob_vy_inv; 
  out[15] += 1.7320508075688772*hamil[3]*pt_conf_3[3]*dv11*jacob_vy_inv; 
  out[16] += 3.872983346207417*pt_conf_3[3]*hamil[5]*dv11*jacob_vy_inv; 
  out[17] += 3.872983346207417*pt_conf_3[1]*hamil[7]*dv11*jacob_vy_inv; 
  out[18] += 3.872983346207417*pt_conf_3[2]*hamil[7]*dv11*jacob_vy_inv; 
  out[19] += 1.7320508075688772*hamil[2]*pt_conf_3[6]*dv11*jacob_vy_inv; 
  out[20] += 1.7320508075688772*hamil[2]*pt_conf_3[7]*dv11*jacob_vy_inv; 
  out[21] += 1.7320508075688774*hamil[3]*pt_conf_3[4]*dv11*jacob_vy_inv; 
  out[22] += 1.7320508075688774*hamil[3]*pt_conf_3[5]*dv11*jacob_vy_inv; 
  out[23] += 1.7320508075688772*pt_conf_3[1]*hamil[6]*dv11*jacob_vy_inv; 
  out[24] += 1.7320508075688772*pt_conf_3[2]*hamil[6]*dv11*jacob_vy_inv; 
  out[25] += 3.872983346207417*pt_conf_3[4]*hamil[5]*dv11*jacob_vy_inv; 
  out[26] += 3.872983346207417*hamil[5]*pt_conf_3[5]*dv11*jacob_vy_inv; 
  out[31] += 3.872983346207417*pt_conf_3[3]*hamil[7]*dv11*jacob_vy_inv; 
  out[32] += 1.7320508075688774*hamil[3]*pt_conf_3[6]*dv11*jacob_vy_inv; 
  out[33] += 1.7320508075688774*hamil[3]*pt_conf_3[7]*dv11*jacob_vy_inv; 
  out[34] += 1.7320508075688772*pt_conf_3[3]*hamil[6]*dv11*jacob_vy_inv; 
  out[35] += 3.872983346207417*hamil[5]*pt_conf_3[6]*dv11*jacob_vy_inv; 
  out[36] += 3.872983346207417*hamil[5]*pt_conf_3[7]*dv11*jacob_vy_inv; 
  out[37] += 3.872983346207417*pt_conf_3[4]*hamil[7]*dv11*jacob_vy_inv; 
  out[38] += 3.872983346207417*pt_conf_3[5]*hamil[7]*dv11*jacob_vy_inv; 
  out[44] += 3.872983346207417*pt_conf_3[6]*hamil[7]*dv11*jacob_vy_inv; 
  out[45] += 3.872983346207417*hamil[7]*pt_conf_3[7]*dv11*jacob_vy_inv; 
  } 
  if (dir == 3) { 
  } 

  const double *pt_conf_4 = &poisson_tensor_conf[32]; 
  if (dir == 2) { 
  out[0] += (1.224744871391589*pt_conf_4[0]*vmap_vx[1]*hamil[3]+1.224744871391589*pt_conf_4[0]*vmap_vx[0]*hamil[2])*dv11*jacob_vy_inv; 
  out[1] += (1.224744871391589*pt_conf_4[1]*vmap_vx[1]*hamil[3]+1.224744871391589*vmap_vx[0]*pt_conf_4[1]*hamil[2])*dv11*jacob_vy_inv; 
  out[2] += (1.224744871391589*vmap_vx[1]*pt_conf_4[2]*hamil[3]+1.224744871391589*vmap_vx[0]*hamil[2]*pt_conf_4[2])*dv11*jacob_vy_inv; 
  out[3] += (1.0954451150103324*pt_conf_4[0]*vmap_vx[1]*hamil[6]+1.224744871391589*pt_conf_4[0]*vmap_vx[0]*hamil[3]+1.224744871391589*pt_conf_4[0]*vmap_vx[1]*hamil[2])*dv11*jacob_vy_inv; 
  out[4] += (2.7386127875258306*pt_conf_4[0]*vmap_vx[1]*hamil[7]+2.7386127875258306*pt_conf_4[0]*vmap_vx[0]*hamil[5])*dv11*jacob_vy_inv; 
  out[5] += (1.224744871391589*vmap_vx[1]*hamil[3]*pt_conf_4[3]+1.224744871391589*vmap_vx[0]*hamil[2]*pt_conf_4[3])*dv11*jacob_vy_inv; 
  out[6] += (1.0954451150103324*pt_conf_4[1]*vmap_vx[1]*hamil[6]+1.224744871391589*vmap_vx[0]*pt_conf_4[1]*hamil[3]+1.224744871391589*pt_conf_4[1]*vmap_vx[1]*hamil[2])*dv11*jacob_vy_inv; 
  out[7] += (1.0954451150103324*vmap_vx[1]*pt_conf_4[2]*hamil[6]+1.224744871391589*vmap_vx[0]*pt_conf_4[2]*hamil[3]+1.224744871391589*vmap_vx[1]*hamil[2]*pt_conf_4[2])*dv11*jacob_vy_inv; 
  out[8] += (2.7386127875258306*pt_conf_4[1]*vmap_vx[1]*hamil[7]+2.7386127875258306*vmap_vx[0]*pt_conf_4[1]*hamil[5])*dv11*jacob_vy_inv; 
  out[9] += (2.7386127875258306*vmap_vx[1]*pt_conf_4[2]*hamil[7]+2.7386127875258306*vmap_vx[0]*pt_conf_4[2]*hamil[5])*dv11*jacob_vy_inv; 
  out[10] += (2.7386127875258306*pt_conf_4[0]*vmap_vx[0]*hamil[7]+2.7386127875258306*pt_conf_4[0]*vmap_vx[1]*hamil[5])*dv11*jacob_vy_inv; 
  out[11] += (1.224744871391589*vmap_vx[1]*hamil[3]*pt_conf_4[4]+1.224744871391589*vmap_vx[0]*hamil[2]*pt_conf_4[4])*dv11*jacob_vy_inv; 
  out[12] += (1.224744871391589*vmap_vx[1]*hamil[3]*pt_conf_4[5]+1.224744871391589*vmap_vx[0]*hamil[2]*pt_conf_4[5])*dv11*jacob_vy_inv; 
  out[13] += (1.224744871391589*pt_conf_4[0]*vmap_vx[0]*hamil[6]+1.0954451150103324*pt_conf_4[0]*vmap_vx[1]*hamil[3])*dv11*jacob_vy_inv; 
  out[15] += (1.0954451150103324*vmap_vx[1]*pt_conf_4[3]*hamil[6]+1.224744871391589*vmap_vx[0]*hamil[3]*pt_conf_4[3]+1.224744871391589*vmap_vx[1]*hamil[2]*pt_conf_4[3])*dv11*jacob_vy_inv; 
  out[16] += (2.7386127875258306*vmap_vx[1]*pt_conf_4[3]*hamil[7]+2.7386127875258306*vmap_vx[0]*pt_conf_4[3]*hamil[5])*dv11*jacob_vy_inv; 
  out[17] += (2.7386127875258306*vmap_vx[0]*pt_conf_4[1]*hamil[7]+2.7386127875258306*pt_conf_4[1]*vmap_vx[1]*hamil[5])*dv11*jacob_vy_inv; 
  out[18] += (2.7386127875258306*vmap_vx[0]*pt_conf_4[2]*hamil[7]+2.7386127875258306*vmap_vx[1]*pt_conf_4[2]*hamil[5])*dv11*jacob_vy_inv; 
  out[19] += (1.224744871391589*vmap_vx[1]*hamil[3]*pt_conf_4[6]+1.224744871391589*vmap_vx[0]*hamil[2]*pt_conf_4[6])*dv11*jacob_vy_inv; 
  out[20] += (1.224744871391589*vmap_vx[1]*hamil[3]*pt_conf_4[7]+1.224744871391589*vmap_vx[0]*hamil[2]*pt_conf_4[7])*dv11*jacob_vy_inv; 
  out[21] += (1.0954451150103324*vmap_vx[1]*pt_conf_4[4]*hamil[6]+1.224744871391589*vmap_vx[0]*hamil[3]*pt_conf_4[4]+1.224744871391589*vmap_vx[1]*hamil[2]*pt_conf_4[4])*dv11*jacob_vy_inv; 
  out[22] += (1.0954451150103324*vmap_vx[1]*pt_conf_4[5]*hamil[6]+1.224744871391589*vmap_vx[0]*hamil[3]*pt_conf_4[5]+1.224744871391589*vmap_vx[1]*hamil[2]*pt_conf_4[5])*dv11*jacob_vy_inv; 
  out[23] += (1.224744871391589*vmap_vx[0]*pt_conf_4[1]*hamil[6]+1.0954451150103324*pt_conf_4[1]*vmap_vx[1]*hamil[3])*dv11*jacob_vy_inv; 
  out[24] += (1.224744871391589*vmap_vx[0]*pt_conf_4[2]*hamil[6]+1.0954451150103324*vmap_vx[1]*pt_conf_4[2]*hamil[3])*dv11*jacob_vy_inv; 
  out[25] += (2.7386127875258306*vmap_vx[1]*pt_conf_4[4]*hamil[7]+2.7386127875258306*vmap_vx[0]*pt_conf_4[4]*hamil[5])*dv11*jacob_vy_inv; 
  out[26] += (2.7386127875258306*vmap_vx[1]*pt_conf_4[5]*hamil[7]+2.7386127875258306*vmap_vx[0]*hamil[5]*pt_conf_4[5])*dv11*jacob_vy_inv; 
  out[27] += 2.4494897427831783*pt_conf_4[0]*vmap_vx[1]*hamil[7]*dv11*jacob_vy_inv; 
  out[31] += (2.7386127875258306*vmap_vx[0]*pt_conf_4[3]*hamil[7]+2.7386127875258306*vmap_vx[1]*pt_conf_4[3]*hamil[5])*dv11*jacob_vy_inv; 
  out[32] += (1.0954451150103324*vmap_vx[1]*hamil[6]*pt_conf_4[6]+1.224744871391589*vmap_vx[0]*hamil[3]*pt_conf_4[6]+1.224744871391589*vmap_vx[1]*hamil[2]*pt_conf_4[6])*dv11*jacob_vy_inv; 
  out[33] += (1.0954451150103324*vmap_vx[1]*hamil[6]*pt_conf_4[7]+1.224744871391589*vmap_vx[0]*hamil[3]*pt_conf_4[7]+1.224744871391589*vmap_vx[1]*hamil[2]*pt_conf_4[7])*dv11*jacob_vy_inv; 
  out[34] += (1.224744871391589*vmap_vx[0]*pt_conf_4[3]*hamil[6]+1.0954451150103324*vmap_vx[1]*hamil[3]*pt_conf_4[3])*dv11*jacob_vy_inv; 
  out[35] += (2.7386127875258306*vmap_vx[1]*pt_conf_4[6]*hamil[7]+2.7386127875258306*vmap_vx[0]*hamil[5]*pt_conf_4[6])*dv11*jacob_vy_inv; 
  out[36] += (2.7386127875258306*vmap_vx[1]*hamil[7]*pt_conf_4[7]+2.7386127875258306*vmap_vx[0]*hamil[5]*pt_conf_4[7])*dv11*jacob_vy_inv; 
  out[37] += (2.7386127875258306*vmap_vx[0]*pt_conf_4[4]*hamil[7]+2.7386127875258306*vmap_vx[1]*pt_conf_4[4]*hamil[5])*dv11*jacob_vy_inv; 
  out[38] += (2.7386127875258306*vmap_vx[0]*pt_conf_4[5]*hamil[7]+2.7386127875258306*vmap_vx[1]*hamil[5]*pt_conf_4[5])*dv11*jacob_vy_inv; 
  out[39] += 2.4494897427831783*pt_conf_4[1]*vmap_vx[1]*hamil[7]*dv11*jacob_vy_inv; 
  out[40] += 2.4494897427831783*vmap_vx[1]*pt_conf_4[2]*hamil[7]*dv11*jacob_vy_inv; 
  out[44] += (2.7386127875258306*vmap_vx[0]*pt_conf_4[6]*hamil[7]+2.7386127875258306*vmap_vx[1]*hamil[5]*pt_conf_4[6])*dv11*jacob_vy_inv; 
  out[45] += (2.7386127875258306*vmap_vx[0]*hamil[7]*pt_conf_4[7]+2.7386127875258306*vmap_vx[1]*hamil[5]*pt_conf_4[7])*dv11*jacob_vy_inv; 
  out[46] += 2.4494897427831783*vmap_vx[1]*pt_conf_4[3]*hamil[7]*dv11*jacob_vy_inv; 
  } 
  if (dir == 3) { 
  out[0] += (-(2.7386127875258306*pt_conf_4[0]*vmap_vx[1]*hamil[4])-1.224744871391589*pt_conf_4[0]*vmap_vx[0]*hamil[1])*dv10*jacob_vx_inv; 
  out[1] += (-(2.7386127875258306*pt_conf_4[1]*vmap_vx[1]*hamil[4])-1.224744871391589*vmap_vx[0]*hamil[1]*pt_conf_4[1])*dv10*jacob_vx_inv; 
  out[2] += (-(2.7386127875258306*vmap_vx[1]*pt_conf_4[2]*hamil[4])-1.224744871391589*vmap_vx[0]*hamil[1]*pt_conf_4[2])*dv10*jacob_vx_inv; 
  out[3] += (-(2.7386127875258306*pt_conf_4[0]*vmap_vx[0]*hamil[4])-1.224744871391589*pt_conf_4[0]*hamil[1]*vmap_vx[1])*dv10*jacob_vx_inv; 
  out[4] += (-(2.7386127875258306*pt_conf_4[0]*vmap_vx[1]*hamil[6])-1.224744871391589*pt_conf_4[0]*vmap_vx[0]*hamil[3])*dv10*jacob_vx_inv; 
  out[5] += (-(2.7386127875258306*vmap_vx[1]*pt_conf_4[3]*hamil[4])-1.224744871391589*vmap_vx[0]*hamil[1]*pt_conf_4[3])*dv10*jacob_vx_inv; 
  out[6] += (-(2.7386127875258306*vmap_vx[0]*pt_conf_4[1]*hamil[4])-1.224744871391589*hamil[1]*pt_conf_4[1]*vmap_vx[1])*dv10*jacob_vx_inv; 
  out[7] += (-(2.7386127875258306*vmap_vx[0]*pt_conf_4[2]*hamil[4])-1.224744871391589*hamil[1]*vmap_vx[1]*pt_conf_4[2])*dv10*jacob_vx_inv; 
  out[8] += (-(2.7386127875258306*pt_conf_4[1]*vmap_vx[1]*hamil[6])-1.224744871391589*vmap_vx[0]*pt_conf_4[1]*hamil[3])*dv10*jacob_vx_inv; 
  out[9] += (-(2.7386127875258306*vmap_vx[1]*pt_conf_4[2]*hamil[6])-1.224744871391589*vmap_vx[0]*pt_conf_4[2]*hamil[3])*dv10*jacob_vx_inv; 
  out[10] += (-(2.7386127875258306*pt_conf_4[0]*vmap_vx[0]*hamil[6])-1.224744871391589*pt_conf_4[0]*vmap_vx[1]*hamil[3])*dv10*jacob_vx_inv; 
  out[11] += (-(2.7386127875258306*vmap_vx[1]*hamil[4]*pt_conf_4[4])-1.224744871391589*vmap_vx[0]*hamil[1]*pt_conf_4[4])*dv10*jacob_vx_inv; 
  out[12] += (-(2.7386127875258306*vmap_vx[1]*hamil[4]*pt_conf_4[5])-1.224744871391589*vmap_vx[0]*hamil[1]*pt_conf_4[5])*dv10*jacob_vx_inv; 
  out[13] += -(2.4494897427831783*pt_conf_4[0]*vmap_vx[1]*hamil[4]*dv10*jacob_vx_inv); 
  out[14] += -(1.224744871391589*pt_conf_4[0]*vmap_vx[0]*hamil[7]*dv10*jacob_vx_inv); 
  out[15] += (-(2.7386127875258306*vmap_vx[0]*pt_conf_4[3]*hamil[4])-1.224744871391589*hamil[1]*vmap_vx[1]*pt_conf_4[3])*dv10*jacob_vx_inv; 
  out[16] += (-(2.7386127875258306*vmap_vx[1]*pt_conf_4[3]*hamil[6])-1.224744871391589*vmap_vx[0]*hamil[3]*pt_conf_4[3])*dv10*jacob_vx_inv; 
  out[17] += (-(2.7386127875258306*vmap_vx[0]*pt_conf_4[1]*hamil[6])-1.224744871391589*pt_conf_4[1]*vmap_vx[1]*hamil[3])*dv10*jacob_vx_inv; 
  out[18] += (-(2.7386127875258306*vmap_vx[0]*pt_conf_4[2]*hamil[6])-1.224744871391589*vmap_vx[1]*pt_conf_4[2]*hamil[3])*dv10*jacob_vx_inv; 
  out[19] += (-(2.7386127875258306*vmap_vx[1]*hamil[4]*pt_conf_4[6])-1.224744871391589*vmap_vx[0]*hamil[1]*pt_conf_4[6])*dv10*jacob_vx_inv; 
  out[20] += (-(2.7386127875258306*vmap_vx[1]*hamil[4]*pt_conf_4[7])-1.224744871391589*vmap_vx[0]*hamil[1]*pt_conf_4[7])*dv10*jacob_vx_inv; 
  out[21] += (-(2.7386127875258306*vmap_vx[0]*hamil[4]*pt_conf_4[4])-1.224744871391589*hamil[1]*vmap_vx[1]*pt_conf_4[4])*dv10*jacob_vx_inv; 
  out[22] += (-(2.7386127875258306*vmap_vx[0]*hamil[4]*pt_conf_4[5])-1.224744871391589*hamil[1]*vmap_vx[1]*pt_conf_4[5])*dv10*jacob_vx_inv; 
  out[23] += -(2.4494897427831783*pt_conf_4[1]*vmap_vx[1]*hamil[4]*dv10*jacob_vx_inv); 
  out[24] += -(2.4494897427831783*vmap_vx[1]*pt_conf_4[2]*hamil[4]*dv10*jacob_vx_inv); 
  out[25] += (-(2.7386127875258306*vmap_vx[1]*pt_conf_4[4]*hamil[6])-1.224744871391589*vmap_vx[0]*hamil[3]*pt_conf_4[4])*dv10*jacob_vx_inv; 
  out[26] += (-(2.7386127875258306*vmap_vx[1]*pt_conf_4[5]*hamil[6])-1.224744871391589*vmap_vx[0]*hamil[3]*pt_conf_4[5])*dv10*jacob_vx_inv; 
  out[27] += -(2.4494897427831783*pt_conf_4[0]*vmap_vx[1]*hamil[6]*dv10*jacob_vx_inv); 
  out[28] += -(1.224744871391589*vmap_vx[0]*pt_conf_4[1]*hamil[7]*dv10*jacob_vx_inv); 
  out[29] += -(1.224744871391589*vmap_vx[0]*pt_conf_4[2]*hamil[7]*dv10*jacob_vx_inv); 
  out[30] += -(1.224744871391589*pt_conf_4[0]*vmap_vx[1]*hamil[7]*dv10*jacob_vx_inv); 
  out[31] += (-(2.7386127875258306*vmap_vx[0]*pt_conf_4[3]*hamil[6])-1.224744871391589*vmap_vx[1]*hamil[3]*pt_conf_4[3])*dv10*jacob_vx_inv; 
  out[32] += (-(2.7386127875258306*vmap_vx[0]*hamil[4]*pt_conf_4[6])-1.224744871391589*hamil[1]*vmap_vx[1]*pt_conf_4[6])*dv10*jacob_vx_inv; 
  out[33] += (-(2.7386127875258306*vmap_vx[0]*hamil[4]*pt_conf_4[7])-1.224744871391589*hamil[1]*vmap_vx[1]*pt_conf_4[7])*dv10*jacob_vx_inv; 
  out[34] += -(2.4494897427831783*vmap_vx[1]*pt_conf_4[3]*hamil[4]*dv10*jacob_vx_inv); 
  out[35] += (-(2.7386127875258306*vmap_vx[1]*hamil[6]*pt_conf_4[6])-1.224744871391589*vmap_vx[0]*hamil[3]*pt_conf_4[6])*dv10*jacob_vx_inv; 
  out[36] += (-(2.7386127875258306*vmap_vx[1]*hamil[6]*pt_conf_4[7])-1.224744871391589*vmap_vx[0]*hamil[3]*pt_conf_4[7])*dv10*jacob_vx_inv; 
  out[37] += (-(2.7386127875258306*vmap_vx[0]*pt_conf_4[4]*hamil[6])-1.224744871391589*vmap_vx[1]*hamil[3]*pt_conf_4[4])*dv10*jacob_vx_inv; 
  out[38] += (-(2.7386127875258306*vmap_vx[0]*pt_conf_4[5]*hamil[6])-1.224744871391589*vmap_vx[1]*hamil[3]*pt_conf_4[5])*dv10*jacob_vx_inv; 
  out[39] += -(2.4494897427831783*pt_conf_4[1]*vmap_vx[1]*hamil[6]*dv10*jacob_vx_inv); 
  out[40] += -(2.4494897427831783*vmap_vx[1]*pt_conf_4[2]*hamil[6]*dv10*jacob_vx_inv); 
  out[41] += -(1.224744871391589*vmap_vx[0]*pt_conf_4[3]*hamil[7]*dv10*jacob_vx_inv); 
  out[42] += -(1.224744871391589*pt_conf_4[1]*vmap_vx[1]*hamil[7]*dv10*jacob_vx_inv); 
  out[43] += -(1.224744871391589*vmap_vx[1]*pt_conf_4[2]*hamil[7]*dv10*jacob_vx_inv); 
  out[44] += (-(2.7386127875258306*vmap_vx[0]*hamil[6]*pt_conf_4[6])-1.224744871391589*vmap_vx[1]*hamil[3]*pt_conf_4[6])*dv10*jacob_vx_inv; 
  out[45] += (-(2.7386127875258306*vmap_vx[0]*hamil[6]*pt_conf_4[7])-1.224744871391589*vmap_vx[1]*hamil[3]*pt_conf_4[7])*dv10*jacob_vx_inv; 
  out[46] += -(2.4494897427831783*vmap_vx[1]*pt_conf_4[3]*hamil[6]*dv10*jacob_vx_inv); 
  out[47] += -(1.224744871391589*vmap_vx[1]*pt_conf_4[3]*hamil[7]*dv10*jacob_vx_inv); 
  } 

  const double *pt_conf_5 = &poisson_tensor_conf[40]; 
  if (dir == 2) { 
  out[0] += (2.7386127875258306*pt_conf_5[0]*vmap_vy[1]*hamil[5]+1.224744871391589*pt_conf_5[0]*vmap_vy[0]*hamil[2])*dv11*jacob_vy_inv; 
  out[1] += (2.7386127875258306*pt_conf_5[1]*vmap_vy[1]*hamil[5]+1.224744871391589*vmap_vy[0]*pt_conf_5[1]*hamil[2])*dv11*jacob_vy_inv; 
  out[2] += (2.7386127875258306*vmap_vy[1]*pt_conf_5[2]*hamil[5]+1.224744871391589*vmap_vy[0]*hamil[2]*pt_conf_5[2])*dv11*jacob_vy_inv; 
  out[3] += (2.7386127875258306*pt_conf_5[0]*vmap_vy[1]*hamil[7]+1.224744871391589*pt_conf_5[0]*vmap_vy[0]*hamil[3])*dv11*jacob_vy_inv; 
  out[4] += (2.7386127875258306*pt_conf_5[0]*vmap_vy[0]*hamil[5]+1.224744871391589*pt_conf_5[0]*vmap_vy[1]*hamil[2])*dv11*jacob_vy_inv; 
  out[5] += (2.7386127875258306*vmap_vy[1]*pt_conf_5[3]*hamil[5]+1.224744871391589*vmap_vy[0]*hamil[2]*pt_conf_5[3])*dv11*jacob_vy_inv; 
  out[6] += (2.7386127875258306*pt_conf_5[1]*vmap_vy[1]*hamil[7]+1.224744871391589*vmap_vy[0]*pt_conf_5[1]*hamil[3])*dv11*jacob_vy_inv; 
  out[7] += (2.7386127875258306*vmap_vy[1]*pt_conf_5[2]*hamil[7]+1.224744871391589*vmap_vy[0]*pt_conf_5[2]*hamil[3])*dv11*jacob_vy_inv; 
  out[8] += (2.7386127875258306*vmap_vy[0]*pt_conf_5[1]*hamil[5]+1.224744871391589*pt_conf_5[1]*vmap_vy[1]*hamil[2])*dv11*jacob_vy_inv; 
  out[9] += (2.7386127875258306*vmap_vy[0]*pt_conf_5[2]*hamil[5]+1.224744871391589*vmap_vy[1]*hamil[2]*pt_conf_5[2])*dv11*jacob_vy_inv; 
  out[10] += (2.7386127875258306*pt_conf_5[0]*vmap_vy[0]*hamil[7]+1.224744871391589*pt_conf_5[0]*vmap_vy[1]*hamil[3])*dv11*jacob_vy_inv; 
  out[11] += (2.7386127875258306*vmap_vy[1]*pt_conf_5[4]*hamil[5]+1.224744871391589*vmap_vy[0]*hamil[2]*pt_conf_5[4])*dv11*jacob_vy_inv; 
  out[12] += (2.7386127875258306*vmap_vy[1]*hamil[5]*pt_conf_5[5]+1.224744871391589*vmap_vy[0]*hamil[2]*pt_conf_5[5])*dv11*jacob_vy_inv; 
  out[13] += 1.224744871391589*pt_conf_5[0]*vmap_vy[0]*hamil[6]*dv11*jacob_vy_inv; 
  out[14] += 2.4494897427831783*pt_conf_5[0]*vmap_vy[1]*hamil[5]*dv11*jacob_vy_inv; 
  out[15] += (2.7386127875258306*vmap_vy[1]*pt_conf_5[3]*hamil[7]+1.224744871391589*vmap_vy[0]*hamil[3]*pt_conf_5[3])*dv11*jacob_vy_inv; 
  out[16] += (2.7386127875258306*vmap_vy[0]*pt_conf_5[3]*hamil[5]+1.224744871391589*vmap_vy[1]*hamil[2]*pt_conf_5[3])*dv11*jacob_vy_inv; 
  out[17] += (2.7386127875258306*vmap_vy[0]*pt_conf_5[1]*hamil[7]+1.224744871391589*pt_conf_5[1]*vmap_vy[1]*hamil[3])*dv11*jacob_vy_inv; 
  out[18] += (2.7386127875258306*vmap_vy[0]*pt_conf_5[2]*hamil[7]+1.224744871391589*vmap_vy[1]*pt_conf_5[2]*hamil[3])*dv11*jacob_vy_inv; 
  out[19] += (2.7386127875258306*vmap_vy[1]*hamil[5]*pt_conf_5[6]+1.224744871391589*vmap_vy[0]*hamil[2]*pt_conf_5[6])*dv11*jacob_vy_inv; 
  out[20] += (2.7386127875258306*vmap_vy[1]*hamil[5]*pt_conf_5[7]+1.224744871391589*vmap_vy[0]*hamil[2]*pt_conf_5[7])*dv11*jacob_vy_inv; 
  out[21] += (2.7386127875258306*vmap_vy[1]*pt_conf_5[4]*hamil[7]+1.224744871391589*vmap_vy[0]*hamil[3]*pt_conf_5[4])*dv11*jacob_vy_inv; 
  out[22] += (2.7386127875258306*vmap_vy[1]*pt_conf_5[5]*hamil[7]+1.224744871391589*vmap_vy[0]*hamil[3]*pt_conf_5[5])*dv11*jacob_vy_inv; 
  out[23] += 1.224744871391589*vmap_vy[0]*pt_conf_5[1]*hamil[6]*dv11*jacob_vy_inv; 
  out[24] += 1.224744871391589*vmap_vy[0]*pt_conf_5[2]*hamil[6]*dv11*jacob_vy_inv; 
  out[25] += (2.7386127875258306*vmap_vy[0]*pt_conf_5[4]*hamil[5]+1.224744871391589*vmap_vy[1]*hamil[2]*pt_conf_5[4])*dv11*jacob_vy_inv; 
  out[26] += (2.7386127875258306*vmap_vy[0]*hamil[5]*pt_conf_5[5]+1.224744871391589*vmap_vy[1]*hamil[2]*pt_conf_5[5])*dv11*jacob_vy_inv; 
  out[27] += 1.224744871391589*pt_conf_5[0]*vmap_vy[1]*hamil[6]*dv11*jacob_vy_inv; 
  out[28] += 2.4494897427831783*pt_conf_5[1]*vmap_vy[1]*hamil[5]*dv11*jacob_vy_inv; 
  out[29] += 2.4494897427831783*vmap_vy[1]*pt_conf_5[2]*hamil[5]*dv11*jacob_vy_inv; 
  out[30] += 2.4494897427831783*pt_conf_5[0]*vmap_vy[1]*hamil[7]*dv11*jacob_vy_inv; 
  out[31] += (2.7386127875258306*vmap_vy[0]*pt_conf_5[3]*hamil[7]+1.224744871391589*vmap_vy[1]*hamil[3]*pt_conf_5[3])*dv11*jacob_vy_inv; 
  out[32] += (2.7386127875258306*vmap_vy[1]*pt_conf_5[6]*hamil[7]+1.224744871391589*vmap_vy[0]*hamil[3]*pt_conf_5[6])*dv11*jacob_vy_inv; 
  out[33] += (2.7386127875258306*vmap_vy[1]*hamil[7]*pt_conf_5[7]+1.224744871391589*vmap_vy[0]*hamil[3]*pt_conf_5[7])*dv11*jacob_vy_inv; 
  out[34] += 1.224744871391589*vmap_vy[0]*pt_conf_5[3]*hamil[6]*dv11*jacob_vy_inv; 
  out[35] += (2.7386127875258306*vmap_vy[0]*hamil[5]*pt_conf_5[6]+1.224744871391589*vmap_vy[1]*hamil[2]*pt_conf_5[6])*dv11*jacob_vy_inv; 
  out[36] += (2.7386127875258306*vmap_vy[0]*hamil[5]*pt_conf_5[7]+1.224744871391589*vmap_vy[1]*hamil[2]*pt_conf_5[7])*dv11*jacob_vy_inv; 
  out[37] += (2.7386127875258306*vmap_vy[0]*pt_conf_5[4]*hamil[7]+1.224744871391589*vmap_vy[1]*hamil[3]*pt_conf_5[4])*dv11*jacob_vy_inv; 
  out[38] += (2.7386127875258306*vmap_vy[0]*pt_conf_5[5]*hamil[7]+1.224744871391589*vmap_vy[1]*hamil[3]*pt_conf_5[5])*dv11*jacob_vy_inv; 
  out[39] += 1.224744871391589*pt_conf_5[1]*vmap_vy[1]*hamil[6]*dv11*jacob_vy_inv; 
  out[40] += 1.224744871391589*vmap_vy[1]*pt_conf_5[2]*hamil[6]*dv11*jacob_vy_inv; 
  out[41] += 2.4494897427831783*vmap_vy[1]*pt_conf_5[3]*hamil[5]*dv11*jacob_vy_inv; 
  out[42] += 2.4494897427831783*pt_conf_5[1]*vmap_vy[1]*hamil[7]*dv11*jacob_vy_inv; 
  out[43] += 2.4494897427831783*vmap_vy[1]*pt_conf_5[2]*hamil[7]*dv11*jacob_vy_inv; 
  out[44] += (2.7386127875258306*vmap_vy[0]*pt_conf_5[6]*hamil[7]+1.224744871391589*vmap_vy[1]*hamil[3]*pt_conf_5[6])*dv11*jacob_vy_inv; 
  out[45] += (2.7386127875258306*vmap_vy[0]*hamil[7]*pt_conf_5[7]+1.224744871391589*vmap_vy[1]*hamil[3]*pt_conf_5[7])*dv11*jacob_vy_inv; 
  out[46] += 1.224744871391589*vmap_vy[1]*pt_conf_5[3]*hamil[6]*dv11*jacob_vy_inv; 
  out[47] += 2.4494897427831783*vmap_vy[1]*pt_conf_5[3]*hamil[7]*dv11*jacob_vy_inv; 
  } 
  if (dir == 3) { 
  out[0] += (-(1.224744871391589*pt_conf_5[0]*vmap_vy[1]*hamil[3])-1.224744871391589*pt_conf_5[0]*vmap_vy[0]*hamil[1])*dv10*jacob_vx_inv; 
  out[1] += (-(1.224744871391589*pt_conf_5[1]*vmap_vy[1]*hamil[3])-1.224744871391589*vmap_vy[0]*hamil[1]*pt_conf_5[1])*dv10*jacob_vx_inv; 
  out[2] += (-(1.224744871391589*vmap_vy[1]*pt_conf_5[2]*hamil[3])-1.224744871391589*vmap_vy[0]*hamil[1]*pt_conf_5[2])*dv10*jacob_vx_inv; 
  out[3] += (-(2.7386127875258306*pt_conf_5[0]*vmap_vy[1]*hamil[6])-2.7386127875258306*pt_conf_5[0]*vmap_vy[0]*hamil[4])*dv10*jacob_vx_inv; 
  out[4] += (-(1.0954451150103324*pt_conf_5[0]*vmap_vy[1]*hamil[7])-1.224744871391589*pt_conf_5[0]*vmap_vy[0]*hamil[3]-1.224744871391589*pt_conf_5[0]*hamil[1]*vmap_vy[1])*dv10*jacob_vx_inv; 
  out[5] += (-(1.224744871391589*vmap_vy[1]*hamil[3]*pt_conf_5[3])-1.224744871391589*vmap_vy[0]*hamil[1]*pt_conf_5[3])*dv10*jacob_vx_inv; 
  out[6] += (-(2.7386127875258306*pt_conf_5[1]*vmap_vy[1]*hamil[6])-2.7386127875258306*vmap_vy[0]*pt_conf_5[1]*hamil[4])*dv10*jacob_vx_inv; 
  out[7] += (-(2.7386127875258306*vmap_vy[1]*pt_conf_5[2]*hamil[6])-2.7386127875258306*vmap_vy[0]*pt_conf_5[2]*hamil[4])*dv10*jacob_vx_inv; 
  out[8] += (-(1.0954451150103324*pt_conf_5[1]*vmap_vy[1]*hamil[7])-1.224744871391589*vmap_vy[0]*pt_conf_5[1]*hamil[3]-1.224744871391589*hamil[1]*pt_conf_5[1]*vmap_vy[1])*dv10*jacob_vx_inv; 
  out[9] += (-(1.0954451150103324*vmap_vy[1]*pt_conf_5[2]*hamil[7])-1.224744871391589*vmap_vy[0]*pt_conf_5[2]*hamil[3]-1.224744871391589*hamil[1]*vmap_vy[1]*pt_conf_5[2])*dv10*jacob_vx_inv; 
  out[10] += (-(2.7386127875258306*pt_conf_5[0]*vmap_vy[0]*hamil[6])-2.7386127875258306*pt_conf_5[0]*vmap_vy[1]*hamil[4])*dv10*jacob_vx_inv; 
  out[11] += (-(1.224744871391589*vmap_vy[1]*hamil[3]*pt_conf_5[4])-1.224744871391589*vmap_vy[0]*hamil[1]*pt_conf_5[4])*dv10*jacob_vx_inv; 
  out[12] += (-(1.224744871391589*vmap_vy[1]*hamil[3]*pt_conf_5[5])-1.224744871391589*vmap_vy[0]*hamil[1]*pt_conf_5[5])*dv10*jacob_vx_inv; 
  out[14] += (-(1.224744871391589*pt_conf_5[0]*vmap_vy[0]*hamil[7])-1.0954451150103324*pt_conf_5[0]*vmap_vy[1]*hamil[3])*dv10*jacob_vx_inv; 
  out[15] += (-(2.7386127875258306*vmap_vy[1]*pt_conf_5[3]*hamil[6])-2.7386127875258306*vmap_vy[0]*pt_conf_5[3]*hamil[4])*dv10*jacob_vx_inv; 
  out[16] += (-(1.0954451150103324*vmap_vy[1]*pt_conf_5[3]*hamil[7])-1.224744871391589*vmap_vy[0]*hamil[3]*pt_conf_5[3]-1.224744871391589*hamil[1]*vmap_vy[1]*pt_conf_5[3])*dv10*jacob_vx_inv; 
  out[17] += (-(2.7386127875258306*vmap_vy[0]*pt_conf_5[1]*hamil[6])-2.7386127875258306*pt_conf_5[1]*vmap_vy[1]*hamil[4])*dv10*jacob_vx_inv; 
  out[18] += (-(2.7386127875258306*vmap_vy[0]*pt_conf_5[2]*hamil[6])-2.7386127875258306*vmap_vy[1]*pt_conf_5[2]*hamil[4])*dv10*jacob_vx_inv; 
  out[19] += (-(1.224744871391589*vmap_vy[1]*hamil[3]*pt_conf_5[6])-1.224744871391589*vmap_vy[0]*hamil[1]*pt_conf_5[6])*dv10*jacob_vx_inv; 
  out[20] += (-(1.224744871391589*vmap_vy[1]*hamil[3]*pt_conf_5[7])-1.224744871391589*vmap_vy[0]*hamil[1]*pt_conf_5[7])*dv10*jacob_vx_inv; 
  out[21] += (-(2.7386127875258306*vmap_vy[1]*pt_conf_5[4]*hamil[6])-2.7386127875258306*vmap_vy[0]*hamil[4]*pt_conf_5[4])*dv10*jacob_vx_inv; 
  out[22] += (-(2.7386127875258306*vmap_vy[1]*pt_conf_5[5]*hamil[6])-2.7386127875258306*vmap_vy[0]*hamil[4]*pt_conf_5[5])*dv10*jacob_vx_inv; 
  out[25] += (-(1.0954451150103324*vmap_vy[1]*pt_conf_5[4]*hamil[7])-1.224744871391589*vmap_vy[0]*hamil[3]*pt_conf_5[4]-1.224744871391589*hamil[1]*vmap_vy[1]*pt_conf_5[4])*dv10*jacob_vx_inv; 
  out[26] += (-(1.0954451150103324*vmap_vy[1]*pt_conf_5[5]*hamil[7])-1.224744871391589*vmap_vy[0]*hamil[3]*pt_conf_5[5]-1.224744871391589*hamil[1]*vmap_vy[1]*pt_conf_5[5])*dv10*jacob_vx_inv; 
  out[28] += (-(1.224744871391589*vmap_vy[0]*pt_conf_5[1]*hamil[7])-1.0954451150103324*pt_conf_5[1]*vmap_vy[1]*hamil[3])*dv10*jacob_vx_inv; 
  out[29] += (-(1.224744871391589*vmap_vy[0]*pt_conf_5[2]*hamil[7])-1.0954451150103324*vmap_vy[1]*pt_conf_5[2]*hamil[3])*dv10*jacob_vx_inv; 
  out[30] += -(2.4494897427831783*pt_conf_5[0]*vmap_vy[1]*hamil[6]*dv10*jacob_vx_inv); 
  out[31] += (-(2.7386127875258306*vmap_vy[0]*pt_conf_5[3]*hamil[6])-2.7386127875258306*vmap_vy[1]*pt_conf_5[3]*hamil[4])*dv10*jacob_vx_inv; 
  out[32] += (-(2.7386127875258306*vmap_vy[1]*hamil[6]*pt_conf_5[6])-2.7386127875258306*vmap_vy[0]*hamil[4]*pt_conf_5[6])*dv10*jacob_vx_inv; 
  out[33] += (-(2.7386127875258306*vmap_vy[1]*hamil[6]*pt_conf_5[7])-2.7386127875258306*vmap_vy[0]*hamil[4]*pt_conf_5[7])*dv10*jacob_vx_inv; 
  out[35] += (-(1.0954451150103324*vmap_vy[1]*pt_conf_5[6]*hamil[7])-1.224744871391589*vmap_vy[0]*hamil[3]*pt_conf_5[6]-1.224744871391589*hamil[1]*vmap_vy[1]*pt_conf_5[6])*dv10*jacob_vx_inv; 
  out[36] += (-(1.0954451150103324*vmap_vy[1]*hamil[7]*pt_conf_5[7])-1.224744871391589*vmap_vy[0]*hamil[3]*pt_conf_5[7]-1.224744871391589*hamil[1]*vmap_vy[1]*pt_conf_5[7])*dv10*jacob_vx_inv; 
  out[37] += (-(2.7386127875258306*vmap_vy[0]*pt_conf_5[4]*hamil[6])-2.7386127875258306*vmap_vy[1]*hamil[4]*pt_conf_5[4])*dv10*jacob_vx_inv; 
  out[38] += (-(2.7386127875258306*vmap_vy[0]*pt_conf_5[5]*hamil[6])-2.7386127875258306*vmap_vy[1]*hamil[4]*pt_conf_5[5])*dv10*jacob_vx_inv; 
  out[41] += (-(1.224744871391589*vmap_vy[0]*pt_conf_5[3]*hamil[7])-1.0954451150103324*vmap_vy[1]*hamil[3]*pt_conf_5[3])*dv10*jacob_vx_inv; 
  out[42] += -(2.4494897427831783*pt_conf_5[1]*vmap_vy[1]*hamil[6]*dv10*jacob_vx_inv); 
  out[43] += -(2.4494897427831783*vmap_vy[1]*pt_conf_5[2]*hamil[6]*dv10*jacob_vx_inv); 
  out[44] += (-(2.7386127875258306*vmap_vy[0]*hamil[6]*pt_conf_5[6])-2.7386127875258306*vmap_vy[1]*hamil[4]*pt_conf_5[6])*dv10*jacob_vx_inv; 
  out[45] += (-(2.7386127875258306*vmap_vy[0]*hamil[6]*pt_conf_5[7])-2.7386127875258306*vmap_vy[1]*hamil[4]*pt_conf_5[7])*dv10*jacob_vx_inv; 
  out[47] += -(2.4494897427831783*vmap_vy[1]*pt_conf_5[3]*hamil[6]*dv10*jacob_vx_inv); 
  } 

} 
