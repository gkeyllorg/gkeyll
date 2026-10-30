#include <gkyl_nc_hamil_gen_vol_comps_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_vel_sparse_vol_alpha_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
   const double *vmap, const double *jacob_pos, const double *jacob_vel,
   const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT out) 
{ 
  const double dx10 = 2.0/dxv[0]; 
  const double dx11 = 2.0/dxv[1]; 
  const double dv0 = dxv[2]; 
  const double dv1 = dxv[3]; 
  const double dv2 = dxv[4]; 
  const double dv10 = 2.0/dxv[2]; 
  const double dv11 = 2.0/dxv[3]; 
  const double dv12 = 2.0/dxv[4]; 
  const double wx1 = w[2]; 
  const double wx2 = w[3]; 
  const double wx3 = w[4]; 
  const double jacob_cx_inv = 1.0/jacob_pos[0]; 
  const double jacob_cy_inv = 1.0/jacob_pos[3]; 
  const double jacob_vx_inv = 1.0/jacob_vel[0]; 
  const double jacob_vy_inv = 1.0/jacob_vel[3]; 
  const double jacob_vz_inv = 1.0/jacob_vel[6]; 
  const double *vmap_vx = &vmap[0]; 
  const double *vmap_vy = &vmap[4]; 
  const double *vmap_vz = &vmap[8]; 
  
  const double *pt_conf_0 = &poisson_tensor_conf[0]; 
  if (dir == 0) { 
  out[0] += 1.7320508075688772*pt_conf_0[0]*hamil[1]*dv10*jacob_vx_inv; 
  out[1] += 1.7320508075688772*hamil[1]*pt_conf_0[1]*dv10*jacob_vx_inv; 
  out[2] += 1.7320508075688772*hamil[1]*pt_conf_0[2]*dv10*jacob_vx_inv; 
  out[3] += 3.872983346207417*pt_conf_0[0]*hamil[7]*dv10*jacob_vx_inv; 
  out[6] += 1.7320508075688772*hamil[1]*pt_conf_0[3]*dv10*jacob_vx_inv; 
  out[7] += 3.872983346207417*pt_conf_0[1]*hamil[7]*dv10*jacob_vx_inv; 
  out[8] += 3.872983346207417*pt_conf_0[2]*hamil[7]*dv10*jacob_vx_inv; 
  out[16] += 1.7320508075688772*hamil[1]*pt_conf_0[4]*dv10*jacob_vx_inv; 
  out[17] += 1.7320508075688772*hamil[1]*pt_conf_0[5]*dv10*jacob_vx_inv; 
  out[21] += 3.872983346207417*pt_conf_0[3]*hamil[7]*dv10*jacob_vx_inv; 
  out[31] += 1.7320508075688772*hamil[1]*pt_conf_0[6]*dv10*jacob_vx_inv; 
  out[32] += 1.7320508075688772*hamil[1]*pt_conf_0[7]*dv10*jacob_vx_inv; 
  out[33] += 3.872983346207417*pt_conf_0[4]*hamil[7]*dv10*jacob_vx_inv; 
  out[34] += 3.872983346207417*pt_conf_0[5]*hamil[7]*dv10*jacob_vx_inv; 
  out[56] += 3.872983346207417*pt_conf_0[6]*hamil[7]*dv10*jacob_vx_inv; 
  out[57] += 3.872983346207417*hamil[7]*pt_conf_0[7]*dv10*jacob_vx_inv; 
  } 
  if (dir == 2) { 
  } 

  const double *pt_conf_1 = &poisson_tensor_conf[8]; 
  if (dir == 0) { 
  out[0] += 1.7320508075688772*pt_conf_1[0]*hamil[2]*dv11*jacob_vy_inv; 
  out[1] += 1.7320508075688772*pt_conf_1[1]*hamil[2]*dv11*jacob_vy_inv; 
  out[2] += 1.7320508075688772*hamil[2]*pt_conf_1[2]*dv11*jacob_vy_inv; 
  out[4] += 3.872983346207417*pt_conf_1[0]*hamil[8]*dv11*jacob_vy_inv; 
  out[6] += 1.7320508075688772*hamil[2]*pt_conf_1[3]*dv11*jacob_vy_inv; 
  out[9] += 3.872983346207417*pt_conf_1[1]*hamil[8]*dv11*jacob_vy_inv; 
  out[10] += 3.872983346207417*pt_conf_1[2]*hamil[8]*dv11*jacob_vy_inv; 
  out[16] += 1.7320508075688772*hamil[2]*pt_conf_1[4]*dv11*jacob_vy_inv; 
  out[17] += 1.7320508075688772*hamil[2]*pt_conf_1[5]*dv11*jacob_vy_inv; 
  out[22] += 3.872983346207417*pt_conf_1[3]*hamil[8]*dv11*jacob_vy_inv; 
  out[31] += 1.7320508075688772*hamil[2]*pt_conf_1[6]*dv11*jacob_vy_inv; 
  out[32] += 1.7320508075688772*hamil[2]*pt_conf_1[7]*dv11*jacob_vy_inv; 
  out[37] += 3.872983346207417*pt_conf_1[4]*hamil[8]*dv11*jacob_vy_inv; 
  out[38] += 3.872983346207417*pt_conf_1[5]*hamil[8]*dv11*jacob_vy_inv; 
  out[59] += 3.872983346207417*pt_conf_1[6]*hamil[8]*dv11*jacob_vy_inv; 
  out[60] += 3.872983346207417*pt_conf_1[7]*hamil[8]*dv11*jacob_vy_inv; 
  } 
  if (dir == 3) { 
  } 

  const double *pt_conf_2 = &poisson_tensor_conf[16]; 
  if (dir == 0) { 
  out[0] += 1.7320508075688772*pt_conf_2[0]*hamil[3]*dv12*jacob_vz_inv; 
  out[1] += 1.7320508075688772*pt_conf_2[1]*hamil[3]*dv12*jacob_vz_inv; 
  out[2] += 1.7320508075688772*pt_conf_2[2]*hamil[3]*dv12*jacob_vz_inv; 
  out[5] += 3.872983346207417*pt_conf_2[0]*hamil[9]*dv12*jacob_vz_inv; 
  out[6] += 1.7320508075688772*hamil[3]*pt_conf_2[3]*dv12*jacob_vz_inv; 
  out[12] += 3.872983346207417*pt_conf_2[1]*hamil[9]*dv12*jacob_vz_inv; 
  out[13] += 3.872983346207417*pt_conf_2[2]*hamil[9]*dv12*jacob_vz_inv; 
  out[16] += 1.7320508075688772*hamil[3]*pt_conf_2[4]*dv12*jacob_vz_inv; 
  out[17] += 1.7320508075688772*hamil[3]*pt_conf_2[5]*dv12*jacob_vz_inv; 
  out[25] += 3.872983346207417*pt_conf_2[3]*hamil[9]*dv12*jacob_vz_inv; 
  out[31] += 1.7320508075688772*hamil[3]*pt_conf_2[6]*dv12*jacob_vz_inv; 
  out[32] += 1.7320508075688772*hamil[3]*pt_conf_2[7]*dv12*jacob_vz_inv; 
  out[43] += 3.872983346207417*pt_conf_2[4]*hamil[9]*dv12*jacob_vz_inv; 
  out[44] += 3.872983346207417*pt_conf_2[5]*hamil[9]*dv12*jacob_vz_inv; 
  out[68] += 3.872983346207417*pt_conf_2[6]*hamil[9]*dv12*jacob_vz_inv; 
  out[69] += 3.872983346207417*pt_conf_2[7]*hamil[9]*dv12*jacob_vz_inv; 
  } 
  if (dir == 4) { 
  } 

  const double *pt_conf_3 = &poisson_tensor_conf[24]; 
  if (dir == 1) { 
  out[0] += 1.7320508075688772*pt_conf_3[0]*hamil[1]*dv10*jacob_vx_inv; 
  out[1] += 1.7320508075688772*hamil[1]*pt_conf_3[1]*dv10*jacob_vx_inv; 
  out[2] += 1.7320508075688772*hamil[1]*pt_conf_3[2]*dv10*jacob_vx_inv; 
  out[3] += 3.872983346207417*pt_conf_3[0]*hamil[7]*dv10*jacob_vx_inv; 
  out[6] += 1.7320508075688772*hamil[1]*pt_conf_3[3]*dv10*jacob_vx_inv; 
  out[7] += 3.872983346207417*pt_conf_3[1]*hamil[7]*dv10*jacob_vx_inv; 
  out[8] += 3.872983346207417*pt_conf_3[2]*hamil[7]*dv10*jacob_vx_inv; 
  out[16] += 1.7320508075688772*hamil[1]*pt_conf_3[4]*dv10*jacob_vx_inv; 
  out[17] += 1.7320508075688772*hamil[1]*pt_conf_3[5]*dv10*jacob_vx_inv; 
  out[21] += 3.872983346207417*pt_conf_3[3]*hamil[7]*dv10*jacob_vx_inv; 
  out[31] += 1.7320508075688772*hamil[1]*pt_conf_3[6]*dv10*jacob_vx_inv; 
  out[32] += 1.7320508075688772*hamil[1]*pt_conf_3[7]*dv10*jacob_vx_inv; 
  out[33] += 3.872983346207417*pt_conf_3[4]*hamil[7]*dv10*jacob_vx_inv; 
  out[34] += 3.872983346207417*pt_conf_3[5]*hamil[7]*dv10*jacob_vx_inv; 
  out[56] += 3.872983346207417*pt_conf_3[6]*hamil[7]*dv10*jacob_vx_inv; 
  out[57] += 3.872983346207417*hamil[7]*pt_conf_3[7]*dv10*jacob_vx_inv; 
  } 
  if (dir == 2) { 
  } 

  const double *pt_conf_4 = &poisson_tensor_conf[32]; 
  if (dir == 1) { 
  out[0] += 1.7320508075688772*pt_conf_4[0]*hamil[2]*dv11*jacob_vy_inv; 
  out[1] += 1.7320508075688772*pt_conf_4[1]*hamil[2]*dv11*jacob_vy_inv; 
  out[2] += 1.7320508075688772*hamil[2]*pt_conf_4[2]*dv11*jacob_vy_inv; 
  out[4] += 3.872983346207417*pt_conf_4[0]*hamil[8]*dv11*jacob_vy_inv; 
  out[6] += 1.7320508075688772*hamil[2]*pt_conf_4[3]*dv11*jacob_vy_inv; 
  out[9] += 3.872983346207417*pt_conf_4[1]*hamil[8]*dv11*jacob_vy_inv; 
  out[10] += 3.872983346207417*pt_conf_4[2]*hamil[8]*dv11*jacob_vy_inv; 
  out[16] += 1.7320508075688772*hamil[2]*pt_conf_4[4]*dv11*jacob_vy_inv; 
  out[17] += 1.7320508075688772*hamil[2]*pt_conf_4[5]*dv11*jacob_vy_inv; 
  out[22] += 3.872983346207417*pt_conf_4[3]*hamil[8]*dv11*jacob_vy_inv; 
  out[31] += 1.7320508075688772*hamil[2]*pt_conf_4[6]*dv11*jacob_vy_inv; 
  out[32] += 1.7320508075688772*hamil[2]*pt_conf_4[7]*dv11*jacob_vy_inv; 
  out[37] += 3.872983346207417*pt_conf_4[4]*hamil[8]*dv11*jacob_vy_inv; 
  out[38] += 3.872983346207417*pt_conf_4[5]*hamil[8]*dv11*jacob_vy_inv; 
  out[59] += 3.872983346207417*pt_conf_4[6]*hamil[8]*dv11*jacob_vy_inv; 
  out[60] += 3.872983346207417*pt_conf_4[7]*hamil[8]*dv11*jacob_vy_inv; 
  } 
  if (dir == 3) { 
  } 

  const double *pt_conf_5 = &poisson_tensor_conf[40]; 
  if (dir == 1) { 
  out[0] += 1.7320508075688772*pt_conf_5[0]*hamil[3]*dv12*jacob_vz_inv; 
  out[1] += 1.7320508075688772*pt_conf_5[1]*hamil[3]*dv12*jacob_vz_inv; 
  out[2] += 1.7320508075688772*pt_conf_5[2]*hamil[3]*dv12*jacob_vz_inv; 
  out[5] += 3.872983346207417*pt_conf_5[0]*hamil[9]*dv12*jacob_vz_inv; 
  out[6] += 1.7320508075688772*hamil[3]*pt_conf_5[3]*dv12*jacob_vz_inv; 
  out[12] += 3.872983346207417*pt_conf_5[1]*hamil[9]*dv12*jacob_vz_inv; 
  out[13] += 3.872983346207417*pt_conf_5[2]*hamil[9]*dv12*jacob_vz_inv; 
  out[16] += 1.7320508075688772*hamil[3]*pt_conf_5[4]*dv12*jacob_vz_inv; 
  out[17] += 1.7320508075688772*hamil[3]*pt_conf_5[5]*dv12*jacob_vz_inv; 
  out[25] += 3.872983346207417*pt_conf_5[3]*hamil[9]*dv12*jacob_vz_inv; 
  out[31] += 1.7320508075688772*hamil[3]*pt_conf_5[6]*dv12*jacob_vz_inv; 
  out[32] += 1.7320508075688772*hamil[3]*pt_conf_5[7]*dv12*jacob_vz_inv; 
  out[43] += 3.872983346207417*pt_conf_5[4]*hamil[9]*dv12*jacob_vz_inv; 
  out[44] += 3.872983346207417*pt_conf_5[5]*hamil[9]*dv12*jacob_vz_inv; 
  out[68] += 3.872983346207417*pt_conf_5[6]*hamil[9]*dv12*jacob_vz_inv; 
  out[69] += 3.872983346207417*pt_conf_5[7]*hamil[9]*dv12*jacob_vz_inv; 
  } 
  if (dir == 4) { 
  } 

  const double *pt_conf_9 = &poisson_tensor_conf[72]; 
  if (dir == 2) { 
  out[0] += 1.224744871391589*pt_conf_9[0]*vmap_vx[0]*hamil[2]*dv11*jacob_vy_inv; 
  out[1] += 1.224744871391589*vmap_vx[0]*pt_conf_9[1]*hamil[2]*dv11*jacob_vy_inv; 
  out[2] += 1.224744871391589*vmap_vx[0]*hamil[2]*pt_conf_9[2]*dv11*jacob_vy_inv; 
  out[3] += 1.224744871391589*pt_conf_9[0]*vmap_vx[1]*hamil[2]*dv11*jacob_vy_inv; 
  out[4] += 2.7386127875258306*pt_conf_9[0]*vmap_vx[0]*hamil[8]*dv11*jacob_vy_inv; 
  out[6] += 1.224744871391589*vmap_vx[0]*hamil[2]*pt_conf_9[3]*dv11*jacob_vy_inv; 
  out[7] += 1.224744871391589*pt_conf_9[1]*vmap_vx[1]*hamil[2]*dv11*jacob_vy_inv; 
  out[8] += 1.224744871391589*vmap_vx[1]*hamil[2]*pt_conf_9[2]*dv11*jacob_vy_inv; 
  out[9] += 2.7386127875258306*vmap_vx[0]*pt_conf_9[1]*hamil[8]*dv11*jacob_vy_inv; 
  out[10] += 2.7386127875258306*vmap_vx[0]*pt_conf_9[2]*hamil[8]*dv11*jacob_vy_inv; 
  out[11] += 2.7386127875258306*pt_conf_9[0]*vmap_vx[1]*hamil[8]*dv11*jacob_vy_inv; 
  out[16] += 1.224744871391589*vmap_vx[0]*hamil[2]*pt_conf_9[4]*dv11*jacob_vy_inv; 
  out[17] += 1.224744871391589*vmap_vx[0]*hamil[2]*pt_conf_9[5]*dv11*jacob_vy_inv; 
  out[21] += 1.224744871391589*vmap_vx[1]*hamil[2]*pt_conf_9[3]*dv11*jacob_vy_inv; 
  out[22] += 2.7386127875258306*vmap_vx[0]*pt_conf_9[3]*hamil[8]*dv11*jacob_vy_inv; 
  out[23] += 2.7386127875258306*pt_conf_9[1]*vmap_vx[1]*hamil[8]*dv11*jacob_vy_inv; 
  out[24] += 2.7386127875258306*vmap_vx[1]*pt_conf_9[2]*hamil[8]*dv11*jacob_vy_inv; 
  out[31] += 1.224744871391589*vmap_vx[0]*hamil[2]*pt_conf_9[6]*dv11*jacob_vy_inv; 
  out[32] += 1.224744871391589*vmap_vx[0]*hamil[2]*pt_conf_9[7]*dv11*jacob_vy_inv; 
  out[33] += 1.224744871391589*vmap_vx[1]*hamil[2]*pt_conf_9[4]*dv11*jacob_vy_inv; 
  out[34] += 1.224744871391589*vmap_vx[1]*hamil[2]*pt_conf_9[5]*dv11*jacob_vy_inv; 
  out[37] += 2.7386127875258306*vmap_vx[0]*pt_conf_9[4]*hamil[8]*dv11*jacob_vy_inv; 
  out[38] += 2.7386127875258306*vmap_vx[0]*pt_conf_9[5]*hamil[8]*dv11*jacob_vy_inv; 
  out[51] += 2.7386127875258306*vmap_vx[1]*pt_conf_9[3]*hamil[8]*dv11*jacob_vy_inv; 
  out[56] += 1.224744871391589*vmap_vx[1]*hamil[2]*pt_conf_9[6]*dv11*jacob_vy_inv; 
  out[57] += 1.224744871391589*vmap_vx[1]*hamil[2]*pt_conf_9[7]*dv11*jacob_vy_inv; 
  out[59] += 2.7386127875258306*vmap_vx[0]*pt_conf_9[6]*hamil[8]*dv11*jacob_vy_inv; 
  out[60] += 2.7386127875258306*vmap_vx[0]*pt_conf_9[7]*hamil[8]*dv11*jacob_vy_inv; 
  out[61] += 2.7386127875258306*vmap_vx[1]*pt_conf_9[4]*hamil[8]*dv11*jacob_vy_inv; 
  out[62] += 2.7386127875258306*vmap_vx[1]*pt_conf_9[5]*hamil[8]*dv11*jacob_vy_inv; 
  out[87] += 2.7386127875258306*vmap_vx[1]*pt_conf_9[6]*hamil[8]*dv11*jacob_vy_inv; 
  out[88] += 2.7386127875258306*vmap_vx[1]*pt_conf_9[7]*hamil[8]*dv11*jacob_vy_inv; 
  } 
  if (dir == 3) { 
  out[0] += (-(2.7386127875258306*pt_conf_9[0]*vmap_vx[1]*hamil[7])-1.224744871391589*pt_conf_9[0]*vmap_vx[0]*hamil[1])*dv10*jacob_vx_inv; 
  out[1] += (-(2.7386127875258306*pt_conf_9[1]*vmap_vx[1]*hamil[7])-1.224744871391589*vmap_vx[0]*hamil[1]*pt_conf_9[1])*dv10*jacob_vx_inv; 
  out[2] += (-(2.7386127875258306*vmap_vx[1]*pt_conf_9[2]*hamil[7])-1.224744871391589*vmap_vx[0]*hamil[1]*pt_conf_9[2])*dv10*jacob_vx_inv; 
  out[3] += (-(2.7386127875258306*pt_conf_9[0]*vmap_vx[0]*hamil[7])-1.224744871391589*pt_conf_9[0]*hamil[1]*vmap_vx[1])*dv10*jacob_vx_inv; 
  out[6] += (-(2.7386127875258306*vmap_vx[1]*pt_conf_9[3]*hamil[7])-1.224744871391589*vmap_vx[0]*hamil[1]*pt_conf_9[3])*dv10*jacob_vx_inv; 
  out[7] += (-(2.7386127875258306*vmap_vx[0]*pt_conf_9[1]*hamil[7])-1.224744871391589*hamil[1]*pt_conf_9[1]*vmap_vx[1])*dv10*jacob_vx_inv; 
  out[8] += (-(2.7386127875258306*vmap_vx[0]*pt_conf_9[2]*hamil[7])-1.224744871391589*hamil[1]*vmap_vx[1]*pt_conf_9[2])*dv10*jacob_vx_inv; 
  out[16] += (-(2.7386127875258306*vmap_vx[1]*pt_conf_9[4]*hamil[7])-1.224744871391589*vmap_vx[0]*hamil[1]*pt_conf_9[4])*dv10*jacob_vx_inv; 
  out[17] += (-(2.7386127875258306*vmap_vx[1]*pt_conf_9[5]*hamil[7])-1.224744871391589*vmap_vx[0]*hamil[1]*pt_conf_9[5])*dv10*jacob_vx_inv; 
  out[18] += -(2.4494897427831783*pt_conf_9[0]*vmap_vx[1]*hamil[7]*dv10*jacob_vx_inv); 
  out[21] += (-(2.7386127875258306*vmap_vx[0]*pt_conf_9[3]*hamil[7])-1.224744871391589*hamil[1]*vmap_vx[1]*pt_conf_9[3])*dv10*jacob_vx_inv; 
  out[31] += (-(2.7386127875258306*vmap_vx[1]*pt_conf_9[6]*hamil[7])-1.224744871391589*vmap_vx[0]*hamil[1]*pt_conf_9[6])*dv10*jacob_vx_inv; 
  out[32] += (-(2.7386127875258306*vmap_vx[1]*hamil[7]*pt_conf_9[7])-1.224744871391589*vmap_vx[0]*hamil[1]*pt_conf_9[7])*dv10*jacob_vx_inv; 
  out[33] += (-(2.7386127875258306*vmap_vx[0]*pt_conf_9[4]*hamil[7])-1.224744871391589*hamil[1]*vmap_vx[1]*pt_conf_9[4])*dv10*jacob_vx_inv; 
  out[34] += (-(2.7386127875258306*vmap_vx[0]*pt_conf_9[5]*hamil[7])-1.224744871391589*hamil[1]*vmap_vx[1]*pt_conf_9[5])*dv10*jacob_vx_inv; 
  out[35] += -(2.4494897427831783*pt_conf_9[1]*vmap_vx[1]*hamil[7]*dv10*jacob_vx_inv); 
  out[36] += -(2.4494897427831783*vmap_vx[1]*pt_conf_9[2]*hamil[7]*dv10*jacob_vx_inv); 
  out[56] += (-(2.7386127875258306*vmap_vx[0]*pt_conf_9[6]*hamil[7])-1.224744871391589*hamil[1]*vmap_vx[1]*pt_conf_9[6])*dv10*jacob_vx_inv; 
  out[57] += (-(2.7386127875258306*vmap_vx[0]*hamil[7]*pt_conf_9[7])-1.224744871391589*hamil[1]*vmap_vx[1]*pt_conf_9[7])*dv10*jacob_vx_inv; 
  out[58] += -(2.4494897427831783*vmap_vx[1]*pt_conf_9[3]*hamil[7]*dv10*jacob_vx_inv); 
  } 

  const double *pt_conf_10 = &poisson_tensor_conf[80]; 
  if (dir == 2) { 
  out[0] += 1.224744871391589*pt_conf_10[0]*vmap_vx[0]*hamil[3]*dv12*jacob_vz_inv; 
  out[1] += 1.224744871391589*vmap_vx[0]*pt_conf_10[1]*hamil[3]*dv12*jacob_vz_inv; 
  out[2] += 1.224744871391589*vmap_vx[0]*pt_conf_10[2]*hamil[3]*dv12*jacob_vz_inv; 
  out[3] += 1.224744871391589*pt_conf_10[0]*vmap_vx[1]*hamil[3]*dv12*jacob_vz_inv; 
  out[5] += 2.7386127875258306*pt_conf_10[0]*vmap_vx[0]*hamil[9]*dv12*jacob_vz_inv; 
  out[6] += 1.224744871391589*vmap_vx[0]*hamil[3]*pt_conf_10[3]*dv12*jacob_vz_inv; 
  out[7] += 1.224744871391589*pt_conf_10[1]*vmap_vx[1]*hamil[3]*dv12*jacob_vz_inv; 
  out[8] += 1.224744871391589*vmap_vx[1]*pt_conf_10[2]*hamil[3]*dv12*jacob_vz_inv; 
  out[12] += 2.7386127875258306*vmap_vx[0]*pt_conf_10[1]*hamil[9]*dv12*jacob_vz_inv; 
  out[13] += 2.7386127875258306*vmap_vx[0]*pt_conf_10[2]*hamil[9]*dv12*jacob_vz_inv; 
  out[14] += 2.7386127875258306*pt_conf_10[0]*vmap_vx[1]*hamil[9]*dv12*jacob_vz_inv; 
  out[16] += 1.224744871391589*vmap_vx[0]*hamil[3]*pt_conf_10[4]*dv12*jacob_vz_inv; 
  out[17] += 1.224744871391589*vmap_vx[0]*hamil[3]*pt_conf_10[5]*dv12*jacob_vz_inv; 
  out[21] += 1.224744871391589*vmap_vx[1]*hamil[3]*pt_conf_10[3]*dv12*jacob_vz_inv; 
  out[25] += 2.7386127875258306*vmap_vx[0]*pt_conf_10[3]*hamil[9]*dv12*jacob_vz_inv; 
  out[26] += 2.7386127875258306*pt_conf_10[1]*vmap_vx[1]*hamil[9]*dv12*jacob_vz_inv; 
  out[27] += 2.7386127875258306*vmap_vx[1]*pt_conf_10[2]*hamil[9]*dv12*jacob_vz_inv; 
  out[31] += 1.224744871391589*vmap_vx[0]*hamil[3]*pt_conf_10[6]*dv12*jacob_vz_inv; 
  out[32] += 1.224744871391589*vmap_vx[0]*hamil[3]*pt_conf_10[7]*dv12*jacob_vz_inv; 
  out[33] += 1.224744871391589*vmap_vx[1]*hamil[3]*pt_conf_10[4]*dv12*jacob_vz_inv; 
  out[34] += 1.224744871391589*vmap_vx[1]*hamil[3]*pt_conf_10[5]*dv12*jacob_vz_inv; 
  out[43] += 2.7386127875258306*vmap_vx[0]*pt_conf_10[4]*hamil[9]*dv12*jacob_vz_inv; 
  out[44] += 2.7386127875258306*vmap_vx[0]*pt_conf_10[5]*hamil[9]*dv12*jacob_vz_inv; 
  out[52] += 2.7386127875258306*vmap_vx[1]*pt_conf_10[3]*hamil[9]*dv12*jacob_vz_inv; 
  out[56] += 1.224744871391589*vmap_vx[1]*hamil[3]*pt_conf_10[6]*dv12*jacob_vz_inv; 
  out[57] += 1.224744871391589*vmap_vx[1]*hamil[3]*pt_conf_10[7]*dv12*jacob_vz_inv; 
  out[68] += 2.7386127875258306*vmap_vx[0]*pt_conf_10[6]*hamil[9]*dv12*jacob_vz_inv; 
  out[69] += 2.7386127875258306*vmap_vx[0]*pt_conf_10[7]*hamil[9]*dv12*jacob_vz_inv; 
  out[70] += 2.7386127875258306*vmap_vx[1]*pt_conf_10[4]*hamil[9]*dv12*jacob_vz_inv; 
  out[71] += 2.7386127875258306*vmap_vx[1]*pt_conf_10[5]*hamil[9]*dv12*jacob_vz_inv; 
  out[91] += 2.7386127875258306*vmap_vx[1]*pt_conf_10[6]*hamil[9]*dv12*jacob_vz_inv; 
  out[92] += 2.7386127875258306*vmap_vx[1]*pt_conf_10[7]*hamil[9]*dv12*jacob_vz_inv; 
  } 
  if (dir == 4) { 
  out[0] += (-(2.7386127875258306*pt_conf_10[0]*vmap_vx[1]*hamil[7])-1.224744871391589*pt_conf_10[0]*vmap_vx[0]*hamil[1])*dv10*jacob_vx_inv; 
  out[1] += (-(2.7386127875258306*pt_conf_10[1]*vmap_vx[1]*hamil[7])-1.224744871391589*vmap_vx[0]*hamil[1]*pt_conf_10[1])*dv10*jacob_vx_inv; 
  out[2] += (-(2.7386127875258306*vmap_vx[1]*pt_conf_10[2]*hamil[7])-1.224744871391589*vmap_vx[0]*hamil[1]*pt_conf_10[2])*dv10*jacob_vx_inv; 
  out[3] += (-(2.7386127875258306*pt_conf_10[0]*vmap_vx[0]*hamil[7])-1.224744871391589*pt_conf_10[0]*hamil[1]*vmap_vx[1])*dv10*jacob_vx_inv; 
  out[6] += (-(2.7386127875258306*vmap_vx[1]*pt_conf_10[3]*hamil[7])-1.224744871391589*vmap_vx[0]*hamil[1]*pt_conf_10[3])*dv10*jacob_vx_inv; 
  out[7] += (-(2.7386127875258306*vmap_vx[0]*pt_conf_10[1]*hamil[7])-1.224744871391589*hamil[1]*pt_conf_10[1]*vmap_vx[1])*dv10*jacob_vx_inv; 
  out[8] += (-(2.7386127875258306*vmap_vx[0]*pt_conf_10[2]*hamil[7])-1.224744871391589*hamil[1]*vmap_vx[1]*pt_conf_10[2])*dv10*jacob_vx_inv; 
  out[16] += (-(2.7386127875258306*vmap_vx[1]*pt_conf_10[4]*hamil[7])-1.224744871391589*vmap_vx[0]*hamil[1]*pt_conf_10[4])*dv10*jacob_vx_inv; 
  out[17] += (-(2.7386127875258306*vmap_vx[1]*pt_conf_10[5]*hamil[7])-1.224744871391589*vmap_vx[0]*hamil[1]*pt_conf_10[5])*dv10*jacob_vx_inv; 
  out[18] += -(2.4494897427831783*pt_conf_10[0]*vmap_vx[1]*hamil[7]*dv10*jacob_vx_inv); 
  out[21] += (-(2.7386127875258306*vmap_vx[0]*pt_conf_10[3]*hamil[7])-1.224744871391589*hamil[1]*vmap_vx[1]*pt_conf_10[3])*dv10*jacob_vx_inv; 
  out[31] += (-(2.7386127875258306*vmap_vx[1]*pt_conf_10[6]*hamil[7])-1.224744871391589*vmap_vx[0]*hamil[1]*pt_conf_10[6])*dv10*jacob_vx_inv; 
  out[32] += (-(2.7386127875258306*vmap_vx[1]*hamil[7]*pt_conf_10[7])-1.224744871391589*vmap_vx[0]*hamil[1]*pt_conf_10[7])*dv10*jacob_vx_inv; 
  out[33] += (-(2.7386127875258306*vmap_vx[0]*pt_conf_10[4]*hamil[7])-1.224744871391589*hamil[1]*vmap_vx[1]*pt_conf_10[4])*dv10*jacob_vx_inv; 
  out[34] += (-(2.7386127875258306*vmap_vx[0]*pt_conf_10[5]*hamil[7])-1.224744871391589*hamil[1]*vmap_vx[1]*pt_conf_10[5])*dv10*jacob_vx_inv; 
  out[35] += -(2.4494897427831783*pt_conf_10[1]*vmap_vx[1]*hamil[7]*dv10*jacob_vx_inv); 
  out[36] += -(2.4494897427831783*vmap_vx[1]*pt_conf_10[2]*hamil[7]*dv10*jacob_vx_inv); 
  out[56] += (-(2.7386127875258306*vmap_vx[0]*pt_conf_10[6]*hamil[7])-1.224744871391589*hamil[1]*vmap_vx[1]*pt_conf_10[6])*dv10*jacob_vx_inv; 
  out[57] += (-(2.7386127875258306*vmap_vx[0]*hamil[7]*pt_conf_10[7])-1.224744871391589*hamil[1]*vmap_vx[1]*pt_conf_10[7])*dv10*jacob_vx_inv; 
  out[58] += -(2.4494897427831783*vmap_vx[1]*pt_conf_10[3]*hamil[7]*dv10*jacob_vx_inv); 
  } 

  const double *pt_conf_11 = &poisson_tensor_conf[88]; 
  if (dir == 3) { 
  out[0] += 1.224744871391589*pt_conf_11[0]*vmap_vx[0]*hamil[3]*dv12*jacob_vz_inv; 
  out[1] += 1.224744871391589*vmap_vx[0]*pt_conf_11[1]*hamil[3]*dv12*jacob_vz_inv; 
  out[2] += 1.224744871391589*vmap_vx[0]*pt_conf_11[2]*hamil[3]*dv12*jacob_vz_inv; 
  out[3] += 1.224744871391589*pt_conf_11[0]*vmap_vx[1]*hamil[3]*dv12*jacob_vz_inv; 
  out[5] += 2.7386127875258306*pt_conf_11[0]*vmap_vx[0]*hamil[9]*dv12*jacob_vz_inv; 
  out[6] += 1.224744871391589*vmap_vx[0]*hamil[3]*pt_conf_11[3]*dv12*jacob_vz_inv; 
  out[7] += 1.224744871391589*pt_conf_11[1]*vmap_vx[1]*hamil[3]*dv12*jacob_vz_inv; 
  out[8] += 1.224744871391589*vmap_vx[1]*pt_conf_11[2]*hamil[3]*dv12*jacob_vz_inv; 
  out[12] += 2.7386127875258306*vmap_vx[0]*pt_conf_11[1]*hamil[9]*dv12*jacob_vz_inv; 
  out[13] += 2.7386127875258306*vmap_vx[0]*pt_conf_11[2]*hamil[9]*dv12*jacob_vz_inv; 
  out[14] += 2.7386127875258306*pt_conf_11[0]*vmap_vx[1]*hamil[9]*dv12*jacob_vz_inv; 
  out[16] += 1.224744871391589*vmap_vx[0]*hamil[3]*pt_conf_11[4]*dv12*jacob_vz_inv; 
  out[17] += 1.224744871391589*vmap_vx[0]*hamil[3]*pt_conf_11[5]*dv12*jacob_vz_inv; 
  out[21] += 1.224744871391589*vmap_vx[1]*hamil[3]*pt_conf_11[3]*dv12*jacob_vz_inv; 
  out[25] += 2.7386127875258306*vmap_vx[0]*pt_conf_11[3]*hamil[9]*dv12*jacob_vz_inv; 
  out[26] += 2.7386127875258306*pt_conf_11[1]*vmap_vx[1]*hamil[9]*dv12*jacob_vz_inv; 
  out[27] += 2.7386127875258306*vmap_vx[1]*pt_conf_11[2]*hamil[9]*dv12*jacob_vz_inv; 
  out[31] += 1.224744871391589*vmap_vx[0]*hamil[3]*pt_conf_11[6]*dv12*jacob_vz_inv; 
  out[32] += 1.224744871391589*vmap_vx[0]*hamil[3]*pt_conf_11[7]*dv12*jacob_vz_inv; 
  out[33] += 1.224744871391589*vmap_vx[1]*hamil[3]*pt_conf_11[4]*dv12*jacob_vz_inv; 
  out[34] += 1.224744871391589*vmap_vx[1]*hamil[3]*pt_conf_11[5]*dv12*jacob_vz_inv; 
  out[43] += 2.7386127875258306*vmap_vx[0]*pt_conf_11[4]*hamil[9]*dv12*jacob_vz_inv; 
  out[44] += 2.7386127875258306*vmap_vx[0]*pt_conf_11[5]*hamil[9]*dv12*jacob_vz_inv; 
  out[52] += 2.7386127875258306*vmap_vx[1]*pt_conf_11[3]*hamil[9]*dv12*jacob_vz_inv; 
  out[56] += 1.224744871391589*vmap_vx[1]*hamil[3]*pt_conf_11[6]*dv12*jacob_vz_inv; 
  out[57] += 1.224744871391589*vmap_vx[1]*hamil[3]*pt_conf_11[7]*dv12*jacob_vz_inv; 
  out[68] += 2.7386127875258306*vmap_vx[0]*pt_conf_11[6]*hamil[9]*dv12*jacob_vz_inv; 
  out[69] += 2.7386127875258306*vmap_vx[0]*pt_conf_11[7]*hamil[9]*dv12*jacob_vz_inv; 
  out[70] += 2.7386127875258306*vmap_vx[1]*pt_conf_11[4]*hamil[9]*dv12*jacob_vz_inv; 
  out[71] += 2.7386127875258306*vmap_vx[1]*pt_conf_11[5]*hamil[9]*dv12*jacob_vz_inv; 
  out[91] += 2.7386127875258306*vmap_vx[1]*pt_conf_11[6]*hamil[9]*dv12*jacob_vz_inv; 
  out[92] += 2.7386127875258306*vmap_vx[1]*pt_conf_11[7]*hamil[9]*dv12*jacob_vz_inv; 
  } 
  if (dir == 4) { 
  out[0] += -(1.224744871391589*pt_conf_11[0]*vmap_vx[0]*hamil[2]*dv11*jacob_vy_inv); 
  out[1] += -(1.224744871391589*vmap_vx[0]*pt_conf_11[1]*hamil[2]*dv11*jacob_vy_inv); 
  out[2] += -(1.224744871391589*vmap_vx[0]*hamil[2]*pt_conf_11[2]*dv11*jacob_vy_inv); 
  out[3] += -(1.224744871391589*pt_conf_11[0]*vmap_vx[1]*hamil[2]*dv11*jacob_vy_inv); 
  out[4] += -(2.7386127875258306*pt_conf_11[0]*vmap_vx[0]*hamil[8]*dv11*jacob_vy_inv); 
  out[6] += -(1.224744871391589*vmap_vx[0]*hamil[2]*pt_conf_11[3]*dv11*jacob_vy_inv); 
  out[7] += -(1.224744871391589*pt_conf_11[1]*vmap_vx[1]*hamil[2]*dv11*jacob_vy_inv); 
  out[8] += -(1.224744871391589*vmap_vx[1]*hamil[2]*pt_conf_11[2]*dv11*jacob_vy_inv); 
  out[9] += -(2.7386127875258306*vmap_vx[0]*pt_conf_11[1]*hamil[8]*dv11*jacob_vy_inv); 
  out[10] += -(2.7386127875258306*vmap_vx[0]*pt_conf_11[2]*hamil[8]*dv11*jacob_vy_inv); 
  out[11] += -(2.7386127875258306*pt_conf_11[0]*vmap_vx[1]*hamil[8]*dv11*jacob_vy_inv); 
  out[16] += -(1.224744871391589*vmap_vx[0]*hamil[2]*pt_conf_11[4]*dv11*jacob_vy_inv); 
  out[17] += -(1.224744871391589*vmap_vx[0]*hamil[2]*pt_conf_11[5]*dv11*jacob_vy_inv); 
  out[21] += -(1.224744871391589*vmap_vx[1]*hamil[2]*pt_conf_11[3]*dv11*jacob_vy_inv); 
  out[22] += -(2.7386127875258306*vmap_vx[0]*pt_conf_11[3]*hamil[8]*dv11*jacob_vy_inv); 
  out[23] += -(2.7386127875258306*pt_conf_11[1]*vmap_vx[1]*hamil[8]*dv11*jacob_vy_inv); 
  out[24] += -(2.7386127875258306*vmap_vx[1]*pt_conf_11[2]*hamil[8]*dv11*jacob_vy_inv); 
  out[31] += -(1.224744871391589*vmap_vx[0]*hamil[2]*pt_conf_11[6]*dv11*jacob_vy_inv); 
  out[32] += -(1.224744871391589*vmap_vx[0]*hamil[2]*pt_conf_11[7]*dv11*jacob_vy_inv); 
  out[33] += -(1.224744871391589*vmap_vx[1]*hamil[2]*pt_conf_11[4]*dv11*jacob_vy_inv); 
  out[34] += -(1.224744871391589*vmap_vx[1]*hamil[2]*pt_conf_11[5]*dv11*jacob_vy_inv); 
  out[37] += -(2.7386127875258306*vmap_vx[0]*pt_conf_11[4]*hamil[8]*dv11*jacob_vy_inv); 
  out[38] += -(2.7386127875258306*vmap_vx[0]*pt_conf_11[5]*hamil[8]*dv11*jacob_vy_inv); 
  out[51] += -(2.7386127875258306*vmap_vx[1]*pt_conf_11[3]*hamil[8]*dv11*jacob_vy_inv); 
  out[56] += -(1.224744871391589*vmap_vx[1]*hamil[2]*pt_conf_11[6]*dv11*jacob_vy_inv); 
  out[57] += -(1.224744871391589*vmap_vx[1]*hamil[2]*pt_conf_11[7]*dv11*jacob_vy_inv); 
  out[59] += -(2.7386127875258306*vmap_vx[0]*pt_conf_11[6]*hamil[8]*dv11*jacob_vy_inv); 
  out[60] += -(2.7386127875258306*vmap_vx[0]*pt_conf_11[7]*hamil[8]*dv11*jacob_vy_inv); 
  out[61] += -(2.7386127875258306*vmap_vx[1]*pt_conf_11[4]*hamil[8]*dv11*jacob_vy_inv); 
  out[62] += -(2.7386127875258306*vmap_vx[1]*pt_conf_11[5]*hamil[8]*dv11*jacob_vy_inv); 
  out[87] += -(2.7386127875258306*vmap_vx[1]*pt_conf_11[6]*hamil[8]*dv11*jacob_vy_inv); 
  out[88] += -(2.7386127875258306*vmap_vx[1]*pt_conf_11[7]*hamil[8]*dv11*jacob_vy_inv); 
  } 

  const double *pt_conf_12 = &poisson_tensor_conf[96]; 
  if (dir == 2) { 
  out[0] += (2.7386127875258306*pt_conf_12[0]*vmap_vy[1]*hamil[8]+1.224744871391589*pt_conf_12[0]*vmap_vy[0]*hamil[2])*dv11*jacob_vy_inv; 
  out[1] += (2.7386127875258306*pt_conf_12[1]*vmap_vy[1]*hamil[8]+1.224744871391589*vmap_vy[0]*pt_conf_12[1]*hamil[2])*dv11*jacob_vy_inv; 
  out[2] += (2.7386127875258306*vmap_vy[1]*pt_conf_12[2]*hamil[8]+1.224744871391589*vmap_vy[0]*hamil[2]*pt_conf_12[2])*dv11*jacob_vy_inv; 
  out[4] += (2.7386127875258306*pt_conf_12[0]*vmap_vy[0]*hamil[8]+1.224744871391589*pt_conf_12[0]*vmap_vy[1]*hamil[2])*dv11*jacob_vy_inv; 
  out[6] += (2.7386127875258306*vmap_vy[1]*pt_conf_12[3]*hamil[8]+1.224744871391589*vmap_vy[0]*hamil[2]*pt_conf_12[3])*dv11*jacob_vy_inv; 
  out[9] += (2.7386127875258306*vmap_vy[0]*pt_conf_12[1]*hamil[8]+1.224744871391589*pt_conf_12[1]*vmap_vy[1]*hamil[2])*dv11*jacob_vy_inv; 
  out[10] += (2.7386127875258306*vmap_vy[0]*pt_conf_12[2]*hamil[8]+1.224744871391589*vmap_vy[1]*hamil[2]*pt_conf_12[2])*dv11*jacob_vy_inv; 
  out[16] += (2.7386127875258306*vmap_vy[1]*pt_conf_12[4]*hamil[8]+1.224744871391589*vmap_vy[0]*hamil[2]*pt_conf_12[4])*dv11*jacob_vy_inv; 
  out[17] += (2.7386127875258306*vmap_vy[1]*pt_conf_12[5]*hamil[8]+1.224744871391589*vmap_vy[0]*hamil[2]*pt_conf_12[5])*dv11*jacob_vy_inv; 
  out[19] += 2.4494897427831783*pt_conf_12[0]*vmap_vy[1]*hamil[8]*dv11*jacob_vy_inv; 
  out[22] += (2.7386127875258306*vmap_vy[0]*pt_conf_12[3]*hamil[8]+1.224744871391589*vmap_vy[1]*hamil[2]*pt_conf_12[3])*dv11*jacob_vy_inv; 
  out[31] += (2.7386127875258306*vmap_vy[1]*pt_conf_12[6]*hamil[8]+1.224744871391589*vmap_vy[0]*hamil[2]*pt_conf_12[6])*dv11*jacob_vy_inv; 
  out[32] += (2.7386127875258306*vmap_vy[1]*pt_conf_12[7]*hamil[8]+1.224744871391589*vmap_vy[0]*hamil[2]*pt_conf_12[7])*dv11*jacob_vy_inv; 
  out[37] += (2.7386127875258306*vmap_vy[0]*pt_conf_12[4]*hamil[8]+1.224744871391589*vmap_vy[1]*hamil[2]*pt_conf_12[4])*dv11*jacob_vy_inv; 
  out[38] += (2.7386127875258306*vmap_vy[0]*pt_conf_12[5]*hamil[8]+1.224744871391589*vmap_vy[1]*hamil[2]*pt_conf_12[5])*dv11*jacob_vy_inv; 
  out[40] += 2.4494897427831783*pt_conf_12[1]*vmap_vy[1]*hamil[8]*dv11*jacob_vy_inv; 
  out[41] += 2.4494897427831783*vmap_vy[1]*pt_conf_12[2]*hamil[8]*dv11*jacob_vy_inv; 
  out[59] += (2.7386127875258306*vmap_vy[0]*pt_conf_12[6]*hamil[8]+1.224744871391589*vmap_vy[1]*hamil[2]*pt_conf_12[6])*dv11*jacob_vy_inv; 
  out[60] += (2.7386127875258306*vmap_vy[0]*pt_conf_12[7]*hamil[8]+1.224744871391589*vmap_vy[1]*hamil[2]*pt_conf_12[7])*dv11*jacob_vy_inv; 
  out[65] += 2.4494897427831783*vmap_vy[1]*pt_conf_12[3]*hamil[8]*dv11*jacob_vy_inv; 
  } 
  if (dir == 3) { 
  out[0] += -(1.224744871391589*pt_conf_12[0]*vmap_vy[0]*hamil[1]*dv10*jacob_vx_inv); 
  out[1] += -(1.224744871391589*vmap_vy[0]*hamil[1]*pt_conf_12[1]*dv10*jacob_vx_inv); 
  out[2] += -(1.224744871391589*vmap_vy[0]*hamil[1]*pt_conf_12[2]*dv10*jacob_vx_inv); 
  out[3] += -(2.7386127875258306*pt_conf_12[0]*vmap_vy[0]*hamil[7]*dv10*jacob_vx_inv); 
  out[4] += -(1.224744871391589*pt_conf_12[0]*hamil[1]*vmap_vy[1]*dv10*jacob_vx_inv); 
  out[6] += -(1.224744871391589*vmap_vy[0]*hamil[1]*pt_conf_12[3]*dv10*jacob_vx_inv); 
  out[7] += -(2.7386127875258306*vmap_vy[0]*pt_conf_12[1]*hamil[7]*dv10*jacob_vx_inv); 
  out[8] += -(2.7386127875258306*vmap_vy[0]*pt_conf_12[2]*hamil[7]*dv10*jacob_vx_inv); 
  out[9] += -(1.224744871391589*hamil[1]*pt_conf_12[1]*vmap_vy[1]*dv10*jacob_vx_inv); 
  out[10] += -(1.224744871391589*hamil[1]*vmap_vy[1]*pt_conf_12[2]*dv10*jacob_vx_inv); 
  out[11] += -(2.7386127875258306*pt_conf_12[0]*vmap_vy[1]*hamil[7]*dv10*jacob_vx_inv); 
  out[16] += -(1.224744871391589*vmap_vy[0]*hamil[1]*pt_conf_12[4]*dv10*jacob_vx_inv); 
  out[17] += -(1.224744871391589*vmap_vy[0]*hamil[1]*pt_conf_12[5]*dv10*jacob_vx_inv); 
  out[21] += -(2.7386127875258306*vmap_vy[0]*pt_conf_12[3]*hamil[7]*dv10*jacob_vx_inv); 
  out[22] += -(1.224744871391589*hamil[1]*vmap_vy[1]*pt_conf_12[3]*dv10*jacob_vx_inv); 
  out[23] += -(2.7386127875258306*pt_conf_12[1]*vmap_vy[1]*hamil[7]*dv10*jacob_vx_inv); 
  out[24] += -(2.7386127875258306*vmap_vy[1]*pt_conf_12[2]*hamil[7]*dv10*jacob_vx_inv); 
  out[31] += -(1.224744871391589*vmap_vy[0]*hamil[1]*pt_conf_12[6]*dv10*jacob_vx_inv); 
  out[32] += -(1.224744871391589*vmap_vy[0]*hamil[1]*pt_conf_12[7]*dv10*jacob_vx_inv); 
  out[33] += -(2.7386127875258306*vmap_vy[0]*pt_conf_12[4]*hamil[7]*dv10*jacob_vx_inv); 
  out[34] += -(2.7386127875258306*vmap_vy[0]*pt_conf_12[5]*hamil[7]*dv10*jacob_vx_inv); 
  out[37] += -(1.224744871391589*hamil[1]*vmap_vy[1]*pt_conf_12[4]*dv10*jacob_vx_inv); 
  out[38] += -(1.224744871391589*hamil[1]*vmap_vy[1]*pt_conf_12[5]*dv10*jacob_vx_inv); 
  out[51] += -(2.7386127875258306*vmap_vy[1]*pt_conf_12[3]*hamil[7]*dv10*jacob_vx_inv); 
  out[56] += -(2.7386127875258306*vmap_vy[0]*pt_conf_12[6]*hamil[7]*dv10*jacob_vx_inv); 
  out[57] += -(2.7386127875258306*vmap_vy[0]*hamil[7]*pt_conf_12[7]*dv10*jacob_vx_inv); 
  out[59] += -(1.224744871391589*hamil[1]*vmap_vy[1]*pt_conf_12[6]*dv10*jacob_vx_inv); 
  out[60] += -(1.224744871391589*hamil[1]*vmap_vy[1]*pt_conf_12[7]*dv10*jacob_vx_inv); 
  out[61] += -(2.7386127875258306*vmap_vy[1]*pt_conf_12[4]*hamil[7]*dv10*jacob_vx_inv); 
  out[62] += -(2.7386127875258306*vmap_vy[1]*pt_conf_12[5]*hamil[7]*dv10*jacob_vx_inv); 
  out[87] += -(2.7386127875258306*vmap_vy[1]*pt_conf_12[6]*hamil[7]*dv10*jacob_vx_inv); 
  out[88] += -(2.7386127875258306*vmap_vy[1]*hamil[7]*pt_conf_12[7]*dv10*jacob_vx_inv); 
  } 

  const double *pt_conf_13 = &poisson_tensor_conf[104]; 
  if (dir == 2) { 
  out[0] += 1.224744871391589*pt_conf_13[0]*vmap_vy[0]*hamil[3]*dv12*jacob_vz_inv; 
  out[1] += 1.224744871391589*vmap_vy[0]*pt_conf_13[1]*hamil[3]*dv12*jacob_vz_inv; 
  out[2] += 1.224744871391589*vmap_vy[0]*pt_conf_13[2]*hamil[3]*dv12*jacob_vz_inv; 
  out[4] += 1.224744871391589*pt_conf_13[0]*vmap_vy[1]*hamil[3]*dv12*jacob_vz_inv; 
  out[5] += 2.7386127875258306*pt_conf_13[0]*vmap_vy[0]*hamil[9]*dv12*jacob_vz_inv; 
  out[6] += 1.224744871391589*vmap_vy[0]*hamil[3]*pt_conf_13[3]*dv12*jacob_vz_inv; 
  out[9] += 1.224744871391589*pt_conf_13[1]*vmap_vy[1]*hamil[3]*dv12*jacob_vz_inv; 
  out[10] += 1.224744871391589*vmap_vy[1]*pt_conf_13[2]*hamil[3]*dv12*jacob_vz_inv; 
  out[12] += 2.7386127875258306*vmap_vy[0]*pt_conf_13[1]*hamil[9]*dv12*jacob_vz_inv; 
  out[13] += 2.7386127875258306*vmap_vy[0]*pt_conf_13[2]*hamil[9]*dv12*jacob_vz_inv; 
  out[15] += 2.7386127875258306*pt_conf_13[0]*vmap_vy[1]*hamil[9]*dv12*jacob_vz_inv; 
  out[16] += 1.224744871391589*vmap_vy[0]*hamil[3]*pt_conf_13[4]*dv12*jacob_vz_inv; 
  out[17] += 1.224744871391589*vmap_vy[0]*hamil[3]*pt_conf_13[5]*dv12*jacob_vz_inv; 
  out[22] += 1.224744871391589*vmap_vy[1]*hamil[3]*pt_conf_13[3]*dv12*jacob_vz_inv; 
  out[25] += 2.7386127875258306*vmap_vy[0]*pt_conf_13[3]*hamil[9]*dv12*jacob_vz_inv; 
  out[28] += 2.7386127875258306*pt_conf_13[1]*vmap_vy[1]*hamil[9]*dv12*jacob_vz_inv; 
  out[29] += 2.7386127875258306*vmap_vy[1]*pt_conf_13[2]*hamil[9]*dv12*jacob_vz_inv; 
  out[31] += 1.224744871391589*vmap_vy[0]*hamil[3]*pt_conf_13[6]*dv12*jacob_vz_inv; 
  out[32] += 1.224744871391589*vmap_vy[0]*hamil[3]*pt_conf_13[7]*dv12*jacob_vz_inv; 
  out[37] += 1.224744871391589*vmap_vy[1]*hamil[3]*pt_conf_13[4]*dv12*jacob_vz_inv; 
  out[38] += 1.224744871391589*vmap_vy[1]*hamil[3]*pt_conf_13[5]*dv12*jacob_vz_inv; 
  out[43] += 2.7386127875258306*vmap_vy[0]*pt_conf_13[4]*hamil[9]*dv12*jacob_vz_inv; 
  out[44] += 2.7386127875258306*vmap_vy[0]*pt_conf_13[5]*hamil[9]*dv12*jacob_vz_inv; 
  out[53] += 2.7386127875258306*vmap_vy[1]*pt_conf_13[3]*hamil[9]*dv12*jacob_vz_inv; 
  out[59] += 1.224744871391589*vmap_vy[1]*hamil[3]*pt_conf_13[6]*dv12*jacob_vz_inv; 
  out[60] += 1.224744871391589*vmap_vy[1]*hamil[3]*pt_conf_13[7]*dv12*jacob_vz_inv; 
  out[68] += 2.7386127875258306*vmap_vy[0]*pt_conf_13[6]*hamil[9]*dv12*jacob_vz_inv; 
  out[69] += 2.7386127875258306*vmap_vy[0]*pt_conf_13[7]*hamil[9]*dv12*jacob_vz_inv; 
  out[74] += 2.7386127875258306*vmap_vy[1]*pt_conf_13[4]*hamil[9]*dv12*jacob_vz_inv; 
  out[75] += 2.7386127875258306*vmap_vy[1]*pt_conf_13[5]*hamil[9]*dv12*jacob_vz_inv; 
  out[94] += 2.7386127875258306*vmap_vy[1]*pt_conf_13[6]*hamil[9]*dv12*jacob_vz_inv; 
  out[95] += 2.7386127875258306*vmap_vy[1]*pt_conf_13[7]*hamil[9]*dv12*jacob_vz_inv; 
  } 
  if (dir == 4) { 
  out[0] += -(1.224744871391589*pt_conf_13[0]*vmap_vy[0]*hamil[1]*dv10*jacob_vx_inv); 
  out[1] += -(1.224744871391589*vmap_vy[0]*hamil[1]*pt_conf_13[1]*dv10*jacob_vx_inv); 
  out[2] += -(1.224744871391589*vmap_vy[0]*hamil[1]*pt_conf_13[2]*dv10*jacob_vx_inv); 
  out[3] += -(2.7386127875258306*pt_conf_13[0]*vmap_vy[0]*hamil[7]*dv10*jacob_vx_inv); 
  out[4] += -(1.224744871391589*pt_conf_13[0]*hamil[1]*vmap_vy[1]*dv10*jacob_vx_inv); 
  out[6] += -(1.224744871391589*vmap_vy[0]*hamil[1]*pt_conf_13[3]*dv10*jacob_vx_inv); 
  out[7] += -(2.7386127875258306*vmap_vy[0]*pt_conf_13[1]*hamil[7]*dv10*jacob_vx_inv); 
  out[8] += -(2.7386127875258306*vmap_vy[0]*pt_conf_13[2]*hamil[7]*dv10*jacob_vx_inv); 
  out[9] += -(1.224744871391589*hamil[1]*pt_conf_13[1]*vmap_vy[1]*dv10*jacob_vx_inv); 
  out[10] += -(1.224744871391589*hamil[1]*vmap_vy[1]*pt_conf_13[2]*dv10*jacob_vx_inv); 
  out[11] += -(2.7386127875258306*pt_conf_13[0]*vmap_vy[1]*hamil[7]*dv10*jacob_vx_inv); 
  out[16] += -(1.224744871391589*vmap_vy[0]*hamil[1]*pt_conf_13[4]*dv10*jacob_vx_inv); 
  out[17] += -(1.224744871391589*vmap_vy[0]*hamil[1]*pt_conf_13[5]*dv10*jacob_vx_inv); 
  out[21] += -(2.7386127875258306*vmap_vy[0]*pt_conf_13[3]*hamil[7]*dv10*jacob_vx_inv); 
  out[22] += -(1.224744871391589*hamil[1]*vmap_vy[1]*pt_conf_13[3]*dv10*jacob_vx_inv); 
  out[23] += -(2.7386127875258306*pt_conf_13[1]*vmap_vy[1]*hamil[7]*dv10*jacob_vx_inv); 
  out[24] += -(2.7386127875258306*vmap_vy[1]*pt_conf_13[2]*hamil[7]*dv10*jacob_vx_inv); 
  out[31] += -(1.224744871391589*vmap_vy[0]*hamil[1]*pt_conf_13[6]*dv10*jacob_vx_inv); 
  out[32] += -(1.224744871391589*vmap_vy[0]*hamil[1]*pt_conf_13[7]*dv10*jacob_vx_inv); 
  out[33] += -(2.7386127875258306*vmap_vy[0]*pt_conf_13[4]*hamil[7]*dv10*jacob_vx_inv); 
  out[34] += -(2.7386127875258306*vmap_vy[0]*pt_conf_13[5]*hamil[7]*dv10*jacob_vx_inv); 
  out[37] += -(1.224744871391589*hamil[1]*vmap_vy[1]*pt_conf_13[4]*dv10*jacob_vx_inv); 
  out[38] += -(1.224744871391589*hamil[1]*vmap_vy[1]*pt_conf_13[5]*dv10*jacob_vx_inv); 
  out[51] += -(2.7386127875258306*vmap_vy[1]*pt_conf_13[3]*hamil[7]*dv10*jacob_vx_inv); 
  out[56] += -(2.7386127875258306*vmap_vy[0]*pt_conf_13[6]*hamil[7]*dv10*jacob_vx_inv); 
  out[57] += -(2.7386127875258306*vmap_vy[0]*hamil[7]*pt_conf_13[7]*dv10*jacob_vx_inv); 
  out[59] += -(1.224744871391589*hamil[1]*vmap_vy[1]*pt_conf_13[6]*dv10*jacob_vx_inv); 
  out[60] += -(1.224744871391589*hamil[1]*vmap_vy[1]*pt_conf_13[7]*dv10*jacob_vx_inv); 
  out[61] += -(2.7386127875258306*vmap_vy[1]*pt_conf_13[4]*hamil[7]*dv10*jacob_vx_inv); 
  out[62] += -(2.7386127875258306*vmap_vy[1]*pt_conf_13[5]*hamil[7]*dv10*jacob_vx_inv); 
  out[87] += -(2.7386127875258306*vmap_vy[1]*pt_conf_13[6]*hamil[7]*dv10*jacob_vx_inv); 
  out[88] += -(2.7386127875258306*vmap_vy[1]*hamil[7]*pt_conf_13[7]*dv10*jacob_vx_inv); 
  } 

  const double *pt_conf_14 = &poisson_tensor_conf[112]; 
  if (dir == 3) { 
  out[0] += 1.224744871391589*pt_conf_14[0]*vmap_vy[0]*hamil[3]*dv12*jacob_vz_inv; 
  out[1] += 1.224744871391589*vmap_vy[0]*pt_conf_14[1]*hamil[3]*dv12*jacob_vz_inv; 
  out[2] += 1.224744871391589*vmap_vy[0]*pt_conf_14[2]*hamil[3]*dv12*jacob_vz_inv; 
  out[4] += 1.224744871391589*pt_conf_14[0]*vmap_vy[1]*hamil[3]*dv12*jacob_vz_inv; 
  out[5] += 2.7386127875258306*pt_conf_14[0]*vmap_vy[0]*hamil[9]*dv12*jacob_vz_inv; 
  out[6] += 1.224744871391589*vmap_vy[0]*hamil[3]*pt_conf_14[3]*dv12*jacob_vz_inv; 
  out[9] += 1.224744871391589*pt_conf_14[1]*vmap_vy[1]*hamil[3]*dv12*jacob_vz_inv; 
  out[10] += 1.224744871391589*vmap_vy[1]*pt_conf_14[2]*hamil[3]*dv12*jacob_vz_inv; 
  out[12] += 2.7386127875258306*vmap_vy[0]*pt_conf_14[1]*hamil[9]*dv12*jacob_vz_inv; 
  out[13] += 2.7386127875258306*vmap_vy[0]*pt_conf_14[2]*hamil[9]*dv12*jacob_vz_inv; 
  out[15] += 2.7386127875258306*pt_conf_14[0]*vmap_vy[1]*hamil[9]*dv12*jacob_vz_inv; 
  out[16] += 1.224744871391589*vmap_vy[0]*hamil[3]*pt_conf_14[4]*dv12*jacob_vz_inv; 
  out[17] += 1.224744871391589*vmap_vy[0]*hamil[3]*pt_conf_14[5]*dv12*jacob_vz_inv; 
  out[22] += 1.224744871391589*vmap_vy[1]*hamil[3]*pt_conf_14[3]*dv12*jacob_vz_inv; 
  out[25] += 2.7386127875258306*vmap_vy[0]*pt_conf_14[3]*hamil[9]*dv12*jacob_vz_inv; 
  out[28] += 2.7386127875258306*pt_conf_14[1]*vmap_vy[1]*hamil[9]*dv12*jacob_vz_inv; 
  out[29] += 2.7386127875258306*vmap_vy[1]*pt_conf_14[2]*hamil[9]*dv12*jacob_vz_inv; 
  out[31] += 1.224744871391589*vmap_vy[0]*hamil[3]*pt_conf_14[6]*dv12*jacob_vz_inv; 
  out[32] += 1.224744871391589*vmap_vy[0]*hamil[3]*pt_conf_14[7]*dv12*jacob_vz_inv; 
  out[37] += 1.224744871391589*vmap_vy[1]*hamil[3]*pt_conf_14[4]*dv12*jacob_vz_inv; 
  out[38] += 1.224744871391589*vmap_vy[1]*hamil[3]*pt_conf_14[5]*dv12*jacob_vz_inv; 
  out[43] += 2.7386127875258306*vmap_vy[0]*pt_conf_14[4]*hamil[9]*dv12*jacob_vz_inv; 
  out[44] += 2.7386127875258306*vmap_vy[0]*pt_conf_14[5]*hamil[9]*dv12*jacob_vz_inv; 
  out[53] += 2.7386127875258306*vmap_vy[1]*pt_conf_14[3]*hamil[9]*dv12*jacob_vz_inv; 
  out[59] += 1.224744871391589*vmap_vy[1]*hamil[3]*pt_conf_14[6]*dv12*jacob_vz_inv; 
  out[60] += 1.224744871391589*vmap_vy[1]*hamil[3]*pt_conf_14[7]*dv12*jacob_vz_inv; 
  out[68] += 2.7386127875258306*vmap_vy[0]*pt_conf_14[6]*hamil[9]*dv12*jacob_vz_inv; 
  out[69] += 2.7386127875258306*vmap_vy[0]*pt_conf_14[7]*hamil[9]*dv12*jacob_vz_inv; 
  out[74] += 2.7386127875258306*vmap_vy[1]*pt_conf_14[4]*hamil[9]*dv12*jacob_vz_inv; 
  out[75] += 2.7386127875258306*vmap_vy[1]*pt_conf_14[5]*hamil[9]*dv12*jacob_vz_inv; 
  out[94] += 2.7386127875258306*vmap_vy[1]*pt_conf_14[6]*hamil[9]*dv12*jacob_vz_inv; 
  out[95] += 2.7386127875258306*vmap_vy[1]*pt_conf_14[7]*hamil[9]*dv12*jacob_vz_inv; 
  } 
  if (dir == 4) { 
  out[0] += (-(2.7386127875258306*pt_conf_14[0]*vmap_vy[1]*hamil[8])-1.224744871391589*pt_conf_14[0]*vmap_vy[0]*hamil[2])*dv11*jacob_vy_inv; 
  out[1] += (-(2.7386127875258306*pt_conf_14[1]*vmap_vy[1]*hamil[8])-1.224744871391589*vmap_vy[0]*pt_conf_14[1]*hamil[2])*dv11*jacob_vy_inv; 
  out[2] += (-(2.7386127875258306*vmap_vy[1]*pt_conf_14[2]*hamil[8])-1.224744871391589*vmap_vy[0]*hamil[2]*pt_conf_14[2])*dv11*jacob_vy_inv; 
  out[4] += (-(2.7386127875258306*pt_conf_14[0]*vmap_vy[0]*hamil[8])-1.224744871391589*pt_conf_14[0]*vmap_vy[1]*hamil[2])*dv11*jacob_vy_inv; 
  out[6] += (-(2.7386127875258306*vmap_vy[1]*pt_conf_14[3]*hamil[8])-1.224744871391589*vmap_vy[0]*hamil[2]*pt_conf_14[3])*dv11*jacob_vy_inv; 
  out[9] += (-(2.7386127875258306*vmap_vy[0]*pt_conf_14[1]*hamil[8])-1.224744871391589*pt_conf_14[1]*vmap_vy[1]*hamil[2])*dv11*jacob_vy_inv; 
  out[10] += (-(2.7386127875258306*vmap_vy[0]*pt_conf_14[2]*hamil[8])-1.224744871391589*vmap_vy[1]*hamil[2]*pt_conf_14[2])*dv11*jacob_vy_inv; 
  out[16] += (-(2.7386127875258306*vmap_vy[1]*pt_conf_14[4]*hamil[8])-1.224744871391589*vmap_vy[0]*hamil[2]*pt_conf_14[4])*dv11*jacob_vy_inv; 
  out[17] += (-(2.7386127875258306*vmap_vy[1]*pt_conf_14[5]*hamil[8])-1.224744871391589*vmap_vy[0]*hamil[2]*pt_conf_14[5])*dv11*jacob_vy_inv; 
  out[19] += -(2.4494897427831783*pt_conf_14[0]*vmap_vy[1]*hamil[8]*dv11*jacob_vy_inv); 
  out[22] += (-(2.7386127875258306*vmap_vy[0]*pt_conf_14[3]*hamil[8])-1.224744871391589*vmap_vy[1]*hamil[2]*pt_conf_14[3])*dv11*jacob_vy_inv; 
  out[31] += (-(2.7386127875258306*vmap_vy[1]*pt_conf_14[6]*hamil[8])-1.224744871391589*vmap_vy[0]*hamil[2]*pt_conf_14[6])*dv11*jacob_vy_inv; 
  out[32] += (-(2.7386127875258306*vmap_vy[1]*pt_conf_14[7]*hamil[8])-1.224744871391589*vmap_vy[0]*hamil[2]*pt_conf_14[7])*dv11*jacob_vy_inv; 
  out[37] += (-(2.7386127875258306*vmap_vy[0]*pt_conf_14[4]*hamil[8])-1.224744871391589*vmap_vy[1]*hamil[2]*pt_conf_14[4])*dv11*jacob_vy_inv; 
  out[38] += (-(2.7386127875258306*vmap_vy[0]*pt_conf_14[5]*hamil[8])-1.224744871391589*vmap_vy[1]*hamil[2]*pt_conf_14[5])*dv11*jacob_vy_inv; 
  out[40] += -(2.4494897427831783*pt_conf_14[1]*vmap_vy[1]*hamil[8]*dv11*jacob_vy_inv); 
  out[41] += -(2.4494897427831783*vmap_vy[1]*pt_conf_14[2]*hamil[8]*dv11*jacob_vy_inv); 
  out[59] += (-(2.7386127875258306*vmap_vy[0]*pt_conf_14[6]*hamil[8])-1.224744871391589*vmap_vy[1]*hamil[2]*pt_conf_14[6])*dv11*jacob_vy_inv; 
  out[60] += (-(2.7386127875258306*vmap_vy[0]*pt_conf_14[7]*hamil[8])-1.224744871391589*vmap_vy[1]*hamil[2]*pt_conf_14[7])*dv11*jacob_vy_inv; 
  out[65] += -(2.4494897427831783*vmap_vy[1]*pt_conf_14[3]*hamil[8]*dv11*jacob_vy_inv); 
  } 

  const double *pt_conf_15 = &poisson_tensor_conf[120]; 
  if (dir == 2) { 
  out[0] += 1.224744871391589*pt_conf_15[0]*vmap_vz[0]*hamil[2]*dv11*jacob_vy_inv; 
  out[1] += 1.224744871391589*vmap_vz[0]*pt_conf_15[1]*hamil[2]*dv11*jacob_vy_inv; 
  out[2] += 1.224744871391589*vmap_vz[0]*hamil[2]*pt_conf_15[2]*dv11*jacob_vy_inv; 
  out[4] += 2.7386127875258306*pt_conf_15[0]*vmap_vz[0]*hamil[8]*dv11*jacob_vy_inv; 
  out[5] += 1.224744871391589*pt_conf_15[0]*vmap_vz[1]*hamil[2]*dv11*jacob_vy_inv; 
  out[6] += 1.224744871391589*vmap_vz[0]*hamil[2]*pt_conf_15[3]*dv11*jacob_vy_inv; 
  out[9] += 2.7386127875258306*vmap_vz[0]*pt_conf_15[1]*hamil[8]*dv11*jacob_vy_inv; 
  out[10] += 2.7386127875258306*vmap_vz[0]*pt_conf_15[2]*hamil[8]*dv11*jacob_vy_inv; 
  out[12] += 1.224744871391589*pt_conf_15[1]*vmap_vz[1]*hamil[2]*dv11*jacob_vy_inv; 
  out[13] += 1.224744871391589*vmap_vz[1]*hamil[2]*pt_conf_15[2]*dv11*jacob_vy_inv; 
  out[15] += 2.7386127875258306*pt_conf_15[0]*vmap_vz[1]*hamil[8]*dv11*jacob_vy_inv; 
  out[16] += 1.224744871391589*vmap_vz[0]*hamil[2]*pt_conf_15[4]*dv11*jacob_vy_inv; 
  out[17] += 1.224744871391589*vmap_vz[0]*hamil[2]*pt_conf_15[5]*dv11*jacob_vy_inv; 
  out[22] += 2.7386127875258306*vmap_vz[0]*pt_conf_15[3]*hamil[8]*dv11*jacob_vy_inv; 
  out[25] += 1.224744871391589*vmap_vz[1]*hamil[2]*pt_conf_15[3]*dv11*jacob_vy_inv; 
  out[28] += 2.7386127875258306*pt_conf_15[1]*vmap_vz[1]*hamil[8]*dv11*jacob_vy_inv; 
  out[29] += 2.7386127875258306*vmap_vz[1]*pt_conf_15[2]*hamil[8]*dv11*jacob_vy_inv; 
  out[31] += 1.224744871391589*vmap_vz[0]*hamil[2]*pt_conf_15[6]*dv11*jacob_vy_inv; 
  out[32] += 1.224744871391589*vmap_vz[0]*hamil[2]*pt_conf_15[7]*dv11*jacob_vy_inv; 
  out[37] += 2.7386127875258306*vmap_vz[0]*pt_conf_15[4]*hamil[8]*dv11*jacob_vy_inv; 
  out[38] += 2.7386127875258306*vmap_vz[0]*pt_conf_15[5]*hamil[8]*dv11*jacob_vy_inv; 
  out[43] += 1.224744871391589*vmap_vz[1]*hamil[2]*pt_conf_15[4]*dv11*jacob_vy_inv; 
  out[44] += 1.224744871391589*vmap_vz[1]*hamil[2]*pt_conf_15[5]*dv11*jacob_vy_inv; 
  out[53] += 2.7386127875258306*vmap_vz[1]*pt_conf_15[3]*hamil[8]*dv11*jacob_vy_inv; 
  out[59] += 2.7386127875258306*vmap_vz[0]*pt_conf_15[6]*hamil[8]*dv11*jacob_vy_inv; 
  out[60] += 2.7386127875258306*vmap_vz[0]*pt_conf_15[7]*hamil[8]*dv11*jacob_vy_inv; 
  out[68] += 1.224744871391589*vmap_vz[1]*hamil[2]*pt_conf_15[6]*dv11*jacob_vy_inv; 
  out[69] += 1.224744871391589*vmap_vz[1]*hamil[2]*pt_conf_15[7]*dv11*jacob_vy_inv; 
  out[74] += 2.7386127875258306*vmap_vz[1]*pt_conf_15[4]*hamil[8]*dv11*jacob_vy_inv; 
  out[75] += 2.7386127875258306*vmap_vz[1]*pt_conf_15[5]*hamil[8]*dv11*jacob_vy_inv; 
  out[94] += 2.7386127875258306*vmap_vz[1]*pt_conf_15[6]*hamil[8]*dv11*jacob_vy_inv; 
  out[95] += 2.7386127875258306*vmap_vz[1]*pt_conf_15[7]*hamil[8]*dv11*jacob_vy_inv; 
  } 
  if (dir == 3) { 
  out[0] += -(1.224744871391589*pt_conf_15[0]*vmap_vz[0]*hamil[1]*dv10*jacob_vx_inv); 
  out[1] += -(1.224744871391589*vmap_vz[0]*hamil[1]*pt_conf_15[1]*dv10*jacob_vx_inv); 
  out[2] += -(1.224744871391589*vmap_vz[0]*hamil[1]*pt_conf_15[2]*dv10*jacob_vx_inv); 
  out[3] += -(2.7386127875258306*pt_conf_15[0]*vmap_vz[0]*hamil[7]*dv10*jacob_vx_inv); 
  out[5] += -(1.224744871391589*pt_conf_15[0]*hamil[1]*vmap_vz[1]*dv10*jacob_vx_inv); 
  out[6] += -(1.224744871391589*vmap_vz[0]*hamil[1]*pt_conf_15[3]*dv10*jacob_vx_inv); 
  out[7] += -(2.7386127875258306*vmap_vz[0]*pt_conf_15[1]*hamil[7]*dv10*jacob_vx_inv); 
  out[8] += -(2.7386127875258306*vmap_vz[0]*pt_conf_15[2]*hamil[7]*dv10*jacob_vx_inv); 
  out[12] += -(1.224744871391589*hamil[1]*pt_conf_15[1]*vmap_vz[1]*dv10*jacob_vx_inv); 
  out[13] += -(1.224744871391589*hamil[1]*vmap_vz[1]*pt_conf_15[2]*dv10*jacob_vx_inv); 
  out[14] += -(2.7386127875258306*pt_conf_15[0]*vmap_vz[1]*hamil[7]*dv10*jacob_vx_inv); 
  out[16] += -(1.224744871391589*vmap_vz[0]*hamil[1]*pt_conf_15[4]*dv10*jacob_vx_inv); 
  out[17] += -(1.224744871391589*vmap_vz[0]*hamil[1]*pt_conf_15[5]*dv10*jacob_vx_inv); 
  out[21] += -(2.7386127875258306*vmap_vz[0]*pt_conf_15[3]*hamil[7]*dv10*jacob_vx_inv); 
  out[25] += -(1.224744871391589*hamil[1]*vmap_vz[1]*pt_conf_15[3]*dv10*jacob_vx_inv); 
  out[26] += -(2.7386127875258306*pt_conf_15[1]*vmap_vz[1]*hamil[7]*dv10*jacob_vx_inv); 
  out[27] += -(2.7386127875258306*vmap_vz[1]*pt_conf_15[2]*hamil[7]*dv10*jacob_vx_inv); 
  out[31] += -(1.224744871391589*vmap_vz[0]*hamil[1]*pt_conf_15[6]*dv10*jacob_vx_inv); 
  out[32] += -(1.224744871391589*vmap_vz[0]*hamil[1]*pt_conf_15[7]*dv10*jacob_vx_inv); 
  out[33] += -(2.7386127875258306*vmap_vz[0]*pt_conf_15[4]*hamil[7]*dv10*jacob_vx_inv); 
  out[34] += -(2.7386127875258306*vmap_vz[0]*pt_conf_15[5]*hamil[7]*dv10*jacob_vx_inv); 
  out[43] += -(1.224744871391589*hamil[1]*vmap_vz[1]*pt_conf_15[4]*dv10*jacob_vx_inv); 
  out[44] += -(1.224744871391589*hamil[1]*vmap_vz[1]*pt_conf_15[5]*dv10*jacob_vx_inv); 
  out[52] += -(2.7386127875258306*vmap_vz[1]*pt_conf_15[3]*hamil[7]*dv10*jacob_vx_inv); 
  out[56] += -(2.7386127875258306*vmap_vz[0]*pt_conf_15[6]*hamil[7]*dv10*jacob_vx_inv); 
  out[57] += -(2.7386127875258306*vmap_vz[0]*hamil[7]*pt_conf_15[7]*dv10*jacob_vx_inv); 
  out[68] += -(1.224744871391589*hamil[1]*vmap_vz[1]*pt_conf_15[6]*dv10*jacob_vx_inv); 
  out[69] += -(1.224744871391589*hamil[1]*vmap_vz[1]*pt_conf_15[7]*dv10*jacob_vx_inv); 
  out[70] += -(2.7386127875258306*vmap_vz[1]*pt_conf_15[4]*hamil[7]*dv10*jacob_vx_inv); 
  out[71] += -(2.7386127875258306*vmap_vz[1]*pt_conf_15[5]*hamil[7]*dv10*jacob_vx_inv); 
  out[91] += -(2.7386127875258306*vmap_vz[1]*pt_conf_15[6]*hamil[7]*dv10*jacob_vx_inv); 
  out[92] += -(2.7386127875258306*vmap_vz[1]*hamil[7]*pt_conf_15[7]*dv10*jacob_vx_inv); 
  } 

  const double *pt_conf_16 = &poisson_tensor_conf[128]; 
  if (dir == 2) { 
  out[0] += (2.7386127875258306*pt_conf_16[0]*vmap_vz[1]*hamil[9]+1.224744871391589*pt_conf_16[0]*vmap_vz[0]*hamil[3])*dv12*jacob_vz_inv; 
  out[1] += (2.7386127875258306*pt_conf_16[1]*vmap_vz[1]*hamil[9]+1.224744871391589*vmap_vz[0]*pt_conf_16[1]*hamil[3])*dv12*jacob_vz_inv; 
  out[2] += (2.7386127875258306*vmap_vz[1]*pt_conf_16[2]*hamil[9]+1.224744871391589*vmap_vz[0]*pt_conf_16[2]*hamil[3])*dv12*jacob_vz_inv; 
  out[5] += (2.7386127875258306*pt_conf_16[0]*vmap_vz[0]*hamil[9]+1.224744871391589*pt_conf_16[0]*vmap_vz[1]*hamil[3])*dv12*jacob_vz_inv; 
  out[6] += (2.7386127875258306*vmap_vz[1]*pt_conf_16[3]*hamil[9]+1.224744871391589*vmap_vz[0]*hamil[3]*pt_conf_16[3])*dv12*jacob_vz_inv; 
  out[12] += (2.7386127875258306*vmap_vz[0]*pt_conf_16[1]*hamil[9]+1.224744871391589*pt_conf_16[1]*vmap_vz[1]*hamil[3])*dv12*jacob_vz_inv; 
  out[13] += (2.7386127875258306*vmap_vz[0]*pt_conf_16[2]*hamil[9]+1.224744871391589*vmap_vz[1]*pt_conf_16[2]*hamil[3])*dv12*jacob_vz_inv; 
  out[16] += (2.7386127875258306*vmap_vz[1]*pt_conf_16[4]*hamil[9]+1.224744871391589*vmap_vz[0]*hamil[3]*pt_conf_16[4])*dv12*jacob_vz_inv; 
  out[17] += (2.7386127875258306*vmap_vz[1]*pt_conf_16[5]*hamil[9]+1.224744871391589*vmap_vz[0]*hamil[3]*pt_conf_16[5])*dv12*jacob_vz_inv; 
  out[20] += 2.4494897427831783*pt_conf_16[0]*vmap_vz[1]*hamil[9]*dv12*jacob_vz_inv; 
  out[25] += (2.7386127875258306*vmap_vz[0]*pt_conf_16[3]*hamil[9]+1.224744871391589*vmap_vz[1]*hamil[3]*pt_conf_16[3])*dv12*jacob_vz_inv; 
  out[31] += (2.7386127875258306*vmap_vz[1]*pt_conf_16[6]*hamil[9]+1.224744871391589*vmap_vz[0]*hamil[3]*pt_conf_16[6])*dv12*jacob_vz_inv; 
  out[32] += (2.7386127875258306*vmap_vz[1]*pt_conf_16[7]*hamil[9]+1.224744871391589*vmap_vz[0]*hamil[3]*pt_conf_16[7])*dv12*jacob_vz_inv; 
  out[43] += (2.7386127875258306*vmap_vz[0]*pt_conf_16[4]*hamil[9]+1.224744871391589*vmap_vz[1]*hamil[3]*pt_conf_16[4])*dv12*jacob_vz_inv; 
  out[44] += (2.7386127875258306*vmap_vz[0]*pt_conf_16[5]*hamil[9]+1.224744871391589*vmap_vz[1]*hamil[3]*pt_conf_16[5])*dv12*jacob_vz_inv; 
  out[47] += 2.4494897427831783*pt_conf_16[1]*vmap_vz[1]*hamil[9]*dv12*jacob_vz_inv; 
  out[48] += 2.4494897427831783*vmap_vz[1]*pt_conf_16[2]*hamil[9]*dv12*jacob_vz_inv; 
  out[68] += (2.7386127875258306*vmap_vz[0]*pt_conf_16[6]*hamil[9]+1.224744871391589*vmap_vz[1]*hamil[3]*pt_conf_16[6])*dv12*jacob_vz_inv; 
  out[69] += (2.7386127875258306*vmap_vz[0]*pt_conf_16[7]*hamil[9]+1.224744871391589*vmap_vz[1]*hamil[3]*pt_conf_16[7])*dv12*jacob_vz_inv; 
  out[80] += 2.4494897427831783*vmap_vz[1]*pt_conf_16[3]*hamil[9]*dv12*jacob_vz_inv; 
  } 
  if (dir == 4) { 
  out[0] += -(1.224744871391589*pt_conf_16[0]*vmap_vz[0]*hamil[1]*dv10*jacob_vx_inv); 
  out[1] += -(1.224744871391589*vmap_vz[0]*hamil[1]*pt_conf_16[1]*dv10*jacob_vx_inv); 
  out[2] += -(1.224744871391589*vmap_vz[0]*hamil[1]*pt_conf_16[2]*dv10*jacob_vx_inv); 
  out[3] += -(2.7386127875258306*pt_conf_16[0]*vmap_vz[0]*hamil[7]*dv10*jacob_vx_inv); 
  out[5] += -(1.224744871391589*pt_conf_16[0]*hamil[1]*vmap_vz[1]*dv10*jacob_vx_inv); 
  out[6] += -(1.224744871391589*vmap_vz[0]*hamil[1]*pt_conf_16[3]*dv10*jacob_vx_inv); 
  out[7] += -(2.7386127875258306*vmap_vz[0]*pt_conf_16[1]*hamil[7]*dv10*jacob_vx_inv); 
  out[8] += -(2.7386127875258306*vmap_vz[0]*pt_conf_16[2]*hamil[7]*dv10*jacob_vx_inv); 
  out[12] += -(1.224744871391589*hamil[1]*pt_conf_16[1]*vmap_vz[1]*dv10*jacob_vx_inv); 
  out[13] += -(1.224744871391589*hamil[1]*vmap_vz[1]*pt_conf_16[2]*dv10*jacob_vx_inv); 
  out[14] += -(2.7386127875258306*pt_conf_16[0]*vmap_vz[1]*hamil[7]*dv10*jacob_vx_inv); 
  out[16] += -(1.224744871391589*vmap_vz[0]*hamil[1]*pt_conf_16[4]*dv10*jacob_vx_inv); 
  out[17] += -(1.224744871391589*vmap_vz[0]*hamil[1]*pt_conf_16[5]*dv10*jacob_vx_inv); 
  out[21] += -(2.7386127875258306*vmap_vz[0]*pt_conf_16[3]*hamil[7]*dv10*jacob_vx_inv); 
  out[25] += -(1.224744871391589*hamil[1]*vmap_vz[1]*pt_conf_16[3]*dv10*jacob_vx_inv); 
  out[26] += -(2.7386127875258306*pt_conf_16[1]*vmap_vz[1]*hamil[7]*dv10*jacob_vx_inv); 
  out[27] += -(2.7386127875258306*vmap_vz[1]*pt_conf_16[2]*hamil[7]*dv10*jacob_vx_inv); 
  out[31] += -(1.224744871391589*vmap_vz[0]*hamil[1]*pt_conf_16[6]*dv10*jacob_vx_inv); 
  out[32] += -(1.224744871391589*vmap_vz[0]*hamil[1]*pt_conf_16[7]*dv10*jacob_vx_inv); 
  out[33] += -(2.7386127875258306*vmap_vz[0]*pt_conf_16[4]*hamil[7]*dv10*jacob_vx_inv); 
  out[34] += -(2.7386127875258306*vmap_vz[0]*pt_conf_16[5]*hamil[7]*dv10*jacob_vx_inv); 
  out[43] += -(1.224744871391589*hamil[1]*vmap_vz[1]*pt_conf_16[4]*dv10*jacob_vx_inv); 
  out[44] += -(1.224744871391589*hamil[1]*vmap_vz[1]*pt_conf_16[5]*dv10*jacob_vx_inv); 
  out[52] += -(2.7386127875258306*vmap_vz[1]*pt_conf_16[3]*hamil[7]*dv10*jacob_vx_inv); 
  out[56] += -(2.7386127875258306*vmap_vz[0]*pt_conf_16[6]*hamil[7]*dv10*jacob_vx_inv); 
  out[57] += -(2.7386127875258306*vmap_vz[0]*hamil[7]*pt_conf_16[7]*dv10*jacob_vx_inv); 
  out[68] += -(1.224744871391589*hamil[1]*vmap_vz[1]*pt_conf_16[6]*dv10*jacob_vx_inv); 
  out[69] += -(1.224744871391589*hamil[1]*vmap_vz[1]*pt_conf_16[7]*dv10*jacob_vx_inv); 
  out[70] += -(2.7386127875258306*vmap_vz[1]*pt_conf_16[4]*hamil[7]*dv10*jacob_vx_inv); 
  out[71] += -(2.7386127875258306*vmap_vz[1]*pt_conf_16[5]*hamil[7]*dv10*jacob_vx_inv); 
  out[91] += -(2.7386127875258306*vmap_vz[1]*pt_conf_16[6]*hamil[7]*dv10*jacob_vx_inv); 
  out[92] += -(2.7386127875258306*vmap_vz[1]*hamil[7]*pt_conf_16[7]*dv10*jacob_vx_inv); 
  } 

  const double *pt_conf_17 = &poisson_tensor_conf[136]; 
  if (dir == 3) { 
  out[0] += (2.7386127875258306*pt_conf_17[0]*vmap_vz[1]*hamil[9]+1.224744871391589*pt_conf_17[0]*vmap_vz[0]*hamil[3])*dv12*jacob_vz_inv; 
  out[1] += (2.7386127875258306*pt_conf_17[1]*vmap_vz[1]*hamil[9]+1.224744871391589*vmap_vz[0]*pt_conf_17[1]*hamil[3])*dv12*jacob_vz_inv; 
  out[2] += (2.7386127875258306*vmap_vz[1]*pt_conf_17[2]*hamil[9]+1.224744871391589*vmap_vz[0]*pt_conf_17[2]*hamil[3])*dv12*jacob_vz_inv; 
  out[5] += (2.7386127875258306*pt_conf_17[0]*vmap_vz[0]*hamil[9]+1.224744871391589*pt_conf_17[0]*vmap_vz[1]*hamil[3])*dv12*jacob_vz_inv; 
  out[6] += (2.7386127875258306*vmap_vz[1]*pt_conf_17[3]*hamil[9]+1.224744871391589*vmap_vz[0]*hamil[3]*pt_conf_17[3])*dv12*jacob_vz_inv; 
  out[12] += (2.7386127875258306*vmap_vz[0]*pt_conf_17[1]*hamil[9]+1.224744871391589*pt_conf_17[1]*vmap_vz[1]*hamil[3])*dv12*jacob_vz_inv; 
  out[13] += (2.7386127875258306*vmap_vz[0]*pt_conf_17[2]*hamil[9]+1.224744871391589*vmap_vz[1]*pt_conf_17[2]*hamil[3])*dv12*jacob_vz_inv; 
  out[16] += (2.7386127875258306*vmap_vz[1]*pt_conf_17[4]*hamil[9]+1.224744871391589*vmap_vz[0]*hamil[3]*pt_conf_17[4])*dv12*jacob_vz_inv; 
  out[17] += (2.7386127875258306*vmap_vz[1]*pt_conf_17[5]*hamil[9]+1.224744871391589*vmap_vz[0]*hamil[3]*pt_conf_17[5])*dv12*jacob_vz_inv; 
  out[20] += 2.4494897427831783*pt_conf_17[0]*vmap_vz[1]*hamil[9]*dv12*jacob_vz_inv; 
  out[25] += (2.7386127875258306*vmap_vz[0]*pt_conf_17[3]*hamil[9]+1.224744871391589*vmap_vz[1]*hamil[3]*pt_conf_17[3])*dv12*jacob_vz_inv; 
  out[31] += (2.7386127875258306*vmap_vz[1]*pt_conf_17[6]*hamil[9]+1.224744871391589*vmap_vz[0]*hamil[3]*pt_conf_17[6])*dv12*jacob_vz_inv; 
  out[32] += (2.7386127875258306*vmap_vz[1]*pt_conf_17[7]*hamil[9]+1.224744871391589*vmap_vz[0]*hamil[3]*pt_conf_17[7])*dv12*jacob_vz_inv; 
  out[43] += (2.7386127875258306*vmap_vz[0]*pt_conf_17[4]*hamil[9]+1.224744871391589*vmap_vz[1]*hamil[3]*pt_conf_17[4])*dv12*jacob_vz_inv; 
  out[44] += (2.7386127875258306*vmap_vz[0]*pt_conf_17[5]*hamil[9]+1.224744871391589*vmap_vz[1]*hamil[3]*pt_conf_17[5])*dv12*jacob_vz_inv; 
  out[47] += 2.4494897427831783*pt_conf_17[1]*vmap_vz[1]*hamil[9]*dv12*jacob_vz_inv; 
  out[48] += 2.4494897427831783*vmap_vz[1]*pt_conf_17[2]*hamil[9]*dv12*jacob_vz_inv; 
  out[68] += (2.7386127875258306*vmap_vz[0]*pt_conf_17[6]*hamil[9]+1.224744871391589*vmap_vz[1]*hamil[3]*pt_conf_17[6])*dv12*jacob_vz_inv; 
  out[69] += (2.7386127875258306*vmap_vz[0]*pt_conf_17[7]*hamil[9]+1.224744871391589*vmap_vz[1]*hamil[3]*pt_conf_17[7])*dv12*jacob_vz_inv; 
  out[80] += 2.4494897427831783*vmap_vz[1]*pt_conf_17[3]*hamil[9]*dv12*jacob_vz_inv; 
  } 
  if (dir == 4) { 
  out[0] += -(1.224744871391589*pt_conf_17[0]*vmap_vz[0]*hamil[2]*dv11*jacob_vy_inv); 
  out[1] += -(1.224744871391589*vmap_vz[0]*pt_conf_17[1]*hamil[2]*dv11*jacob_vy_inv); 
  out[2] += -(1.224744871391589*vmap_vz[0]*hamil[2]*pt_conf_17[2]*dv11*jacob_vy_inv); 
  out[4] += -(2.7386127875258306*pt_conf_17[0]*vmap_vz[0]*hamil[8]*dv11*jacob_vy_inv); 
  out[5] += -(1.224744871391589*pt_conf_17[0]*vmap_vz[1]*hamil[2]*dv11*jacob_vy_inv); 
  out[6] += -(1.224744871391589*vmap_vz[0]*hamil[2]*pt_conf_17[3]*dv11*jacob_vy_inv); 
  out[9] += -(2.7386127875258306*vmap_vz[0]*pt_conf_17[1]*hamil[8]*dv11*jacob_vy_inv); 
  out[10] += -(2.7386127875258306*vmap_vz[0]*pt_conf_17[2]*hamil[8]*dv11*jacob_vy_inv); 
  out[12] += -(1.224744871391589*pt_conf_17[1]*vmap_vz[1]*hamil[2]*dv11*jacob_vy_inv); 
  out[13] += -(1.224744871391589*vmap_vz[1]*hamil[2]*pt_conf_17[2]*dv11*jacob_vy_inv); 
  out[15] += -(2.7386127875258306*pt_conf_17[0]*vmap_vz[1]*hamil[8]*dv11*jacob_vy_inv); 
  out[16] += -(1.224744871391589*vmap_vz[0]*hamil[2]*pt_conf_17[4]*dv11*jacob_vy_inv); 
  out[17] += -(1.224744871391589*vmap_vz[0]*hamil[2]*pt_conf_17[5]*dv11*jacob_vy_inv); 
  out[22] += -(2.7386127875258306*vmap_vz[0]*pt_conf_17[3]*hamil[8]*dv11*jacob_vy_inv); 
  out[25] += -(1.224744871391589*vmap_vz[1]*hamil[2]*pt_conf_17[3]*dv11*jacob_vy_inv); 
  out[28] += -(2.7386127875258306*pt_conf_17[1]*vmap_vz[1]*hamil[8]*dv11*jacob_vy_inv); 
  out[29] += -(2.7386127875258306*vmap_vz[1]*pt_conf_17[2]*hamil[8]*dv11*jacob_vy_inv); 
  out[31] += -(1.224744871391589*vmap_vz[0]*hamil[2]*pt_conf_17[6]*dv11*jacob_vy_inv); 
  out[32] += -(1.224744871391589*vmap_vz[0]*hamil[2]*pt_conf_17[7]*dv11*jacob_vy_inv); 
  out[37] += -(2.7386127875258306*vmap_vz[0]*pt_conf_17[4]*hamil[8]*dv11*jacob_vy_inv); 
  out[38] += -(2.7386127875258306*vmap_vz[0]*pt_conf_17[5]*hamil[8]*dv11*jacob_vy_inv); 
  out[43] += -(1.224744871391589*vmap_vz[1]*hamil[2]*pt_conf_17[4]*dv11*jacob_vy_inv); 
  out[44] += -(1.224744871391589*vmap_vz[1]*hamil[2]*pt_conf_17[5]*dv11*jacob_vy_inv); 
  out[53] += -(2.7386127875258306*vmap_vz[1]*pt_conf_17[3]*hamil[8]*dv11*jacob_vy_inv); 
  out[59] += -(2.7386127875258306*vmap_vz[0]*pt_conf_17[6]*hamil[8]*dv11*jacob_vy_inv); 
  out[60] += -(2.7386127875258306*vmap_vz[0]*pt_conf_17[7]*hamil[8]*dv11*jacob_vy_inv); 
  out[68] += -(1.224744871391589*vmap_vz[1]*hamil[2]*pt_conf_17[6]*dv11*jacob_vy_inv); 
  out[69] += -(1.224744871391589*vmap_vz[1]*hamil[2]*pt_conf_17[7]*dv11*jacob_vy_inv); 
  out[74] += -(2.7386127875258306*vmap_vz[1]*pt_conf_17[4]*hamil[8]*dv11*jacob_vy_inv); 
  out[75] += -(2.7386127875258306*vmap_vz[1]*pt_conf_17[5]*hamil[8]*dv11*jacob_vy_inv); 
  out[94] += -(2.7386127875258306*vmap_vz[1]*pt_conf_17[6]*hamil[8]*dv11*jacob_vy_inv); 
  out[95] += -(2.7386127875258306*vmap_vz[1]*pt_conf_17[7]*hamil[8]*dv11*jacob_vy_inv); 
  } 

} 
