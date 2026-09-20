#include <gkyl_nc_hamil_gen_vol_comps_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_vel_sparse_vol_alpha_2x3v_tensor_p1(const double *w, const double *dxv, const int dir,
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
  const double jacob_cy_inv = 1.0/jacob_pos[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double *jacob_vy = &jacob_vel[3]; 
  const double *jacob_vz = &jacob_vel[6]; 
  const double *vmap_vx = &vmap[0]; 
  const double *vmap_vy = &vmap[4]; 
  const double *vmap_vz = &vmap[8]; 
  
  const double *pt_conf_0 = &poisson_tensor_conf[0]; 
  if (dir == 0) { 
  double *outb2 = &out[216]; 
  outb2[0] += 1.7320508075688772*pt_conf_0[0]*hamil[1]*dv10; 
  outb2[1] += 1.7320508075688772*hamil[1]*pt_conf_0[1]*dv10; 
  outb2[2] += 1.7320508075688772*hamil[1]*pt_conf_0[2]*dv10; 
  outb2[3] += 3.872983346207417*pt_conf_0[0]*hamil[7]*dv10; 
  outb2[6] += 1.7320508075688772*hamil[1]*pt_conf_0[3]*dv10; 
  outb2[7] += 3.872983346207417*pt_conf_0[1]*hamil[7]*dv10; 
  outb2[8] += 3.872983346207417*pt_conf_0[2]*hamil[7]*dv10; 
  outb2[19] += 3.872983346207417*pt_conf_0[3]*hamil[7]*dv10; 
  } 
  if (dir == 2) { 
  double *outb0 = &out[0]; 
  } 

  const double *pt_conf_1 = &poisson_tensor_conf[4]; 
  if (dir == 0) { 
  double *outb3 = &out[324]; 
  outb3[0] += 1.7320508075688772*pt_conf_1[0]*hamil[2]*dv11; 
  outb3[1] += 1.7320508075688772*pt_conf_1[1]*hamil[2]*dv11; 
  outb3[2] += 1.7320508075688772*hamil[2]*pt_conf_1[2]*dv11; 
  outb3[4] += 3.872983346207417*pt_conf_1[0]*hamil[8]*dv11; 
  outb3[6] += 1.7320508075688772*hamil[2]*pt_conf_1[3]*dv11; 
  outb3[9] += 3.872983346207417*pt_conf_1[1]*hamil[8]*dv11; 
  outb3[10] += 3.872983346207417*pt_conf_1[2]*hamil[8]*dv11; 
  outb3[20] += 3.872983346207417*pt_conf_1[3]*hamil[8]*dv11; 
  } 
  if (dir == 3) { 
  double *outb0 = &out[0]; 
  } 

  const double *pt_conf_2 = &poisson_tensor_conf[8]; 
  if (dir == 0) { 
  double *outb4 = &out[432]; 
  outb4[0] += 1.7320508075688772*pt_conf_2[0]*hamil[3]*dv12; 
  outb4[1] += 1.7320508075688772*pt_conf_2[1]*hamil[3]*dv12; 
  outb4[2] += 1.7320508075688772*pt_conf_2[2]*hamil[3]*dv12; 
  outb4[5] += 3.872983346207417*pt_conf_2[0]*hamil[9]*dv12; 
  outb4[6] += 1.7320508075688772*hamil[3]*pt_conf_2[3]*dv12; 
  outb4[12] += 3.872983346207417*pt_conf_2[1]*hamil[9]*dv12; 
  outb4[13] += 3.872983346207417*pt_conf_2[2]*hamil[9]*dv12; 
  outb4[23] += 3.872983346207417*pt_conf_2[3]*hamil[9]*dv12; 
  } 
  if (dir == 4) { 
  double *outb0 = &out[0]; 
  } 

  const double *pt_conf_3 = &poisson_tensor_conf[12]; 
  if (dir == 1) { 
  double *outb2 = &out[216]; 
  outb2[0] += 1.7320508075688772*pt_conf_3[0]*hamil[1]*dv10; 
  outb2[1] += 1.7320508075688772*hamil[1]*pt_conf_3[1]*dv10; 
  outb2[2] += 1.7320508075688772*hamil[1]*pt_conf_3[2]*dv10; 
  outb2[3] += 3.872983346207417*pt_conf_3[0]*hamil[7]*dv10; 
  outb2[6] += 1.7320508075688772*hamil[1]*pt_conf_3[3]*dv10; 
  outb2[7] += 3.872983346207417*pt_conf_3[1]*hamil[7]*dv10; 
  outb2[8] += 3.872983346207417*pt_conf_3[2]*hamil[7]*dv10; 
  outb2[19] += 3.872983346207417*pt_conf_3[3]*hamil[7]*dv10; 
  } 
  if (dir == 2) { 
  double *outb1 = &out[108]; 
  } 

  const double *pt_conf_4 = &poisson_tensor_conf[16]; 
  if (dir == 1) { 
  double *outb3 = &out[324]; 
  outb3[0] += 1.7320508075688772*pt_conf_4[0]*hamil[2]*dv11; 
  outb3[1] += 1.7320508075688772*pt_conf_4[1]*hamil[2]*dv11; 
  outb3[2] += 1.7320508075688772*hamil[2]*pt_conf_4[2]*dv11; 
  outb3[4] += 3.872983346207417*pt_conf_4[0]*hamil[8]*dv11; 
  outb3[6] += 1.7320508075688772*hamil[2]*pt_conf_4[3]*dv11; 
  outb3[9] += 3.872983346207417*pt_conf_4[1]*hamil[8]*dv11; 
  outb3[10] += 3.872983346207417*pt_conf_4[2]*hamil[8]*dv11; 
  outb3[20] += 3.872983346207417*pt_conf_4[3]*hamil[8]*dv11; 
  } 
  if (dir == 3) { 
  double *outb1 = &out[108]; 
  } 

  const double *pt_conf_5 = &poisson_tensor_conf[20]; 
  if (dir == 1) { 
  double *outb4 = &out[432]; 
  outb4[0] += 1.7320508075688772*pt_conf_5[0]*hamil[3]*dv12; 
  outb4[1] += 1.7320508075688772*pt_conf_5[1]*hamil[3]*dv12; 
  outb4[2] += 1.7320508075688772*pt_conf_5[2]*hamil[3]*dv12; 
  outb4[5] += 3.872983346207417*pt_conf_5[0]*hamil[9]*dv12; 
  outb4[6] += 1.7320508075688772*hamil[3]*pt_conf_5[3]*dv12; 
  outb4[12] += 3.872983346207417*pt_conf_5[1]*hamil[9]*dv12; 
  outb4[13] += 3.872983346207417*pt_conf_5[2]*hamil[9]*dv12; 
  outb4[23] += 3.872983346207417*pt_conf_5[3]*hamil[9]*dv12; 
  } 
  if (dir == 4) { 
  double *outb1 = &out[108]; 
  } 

  const double *pt_conf_9 = &poisson_tensor_conf[36]; 
  if (dir == 2) { 
  double *outb3 = &out[324]; 
  outb3[0] += 1.224744871391589*pt_conf_9[0]*vmap_vx[0]*hamil[2]*dv11; 
  outb3[1] += 1.224744871391589*vmap_vx[0]*pt_conf_9[1]*hamil[2]*dv11; 
  outb3[2] += 1.224744871391589*vmap_vx[0]*hamil[2]*pt_conf_9[2]*dv11; 
  outb3[3] += 1.224744871391589*pt_conf_9[0]*vmap_vx[1]*hamil[2]*dv11; 
  outb3[4] += 2.7386127875258306*pt_conf_9[0]*vmap_vx[0]*hamil[8]*dv11; 
  outb3[6] += 1.224744871391589*vmap_vx[0]*hamil[2]*pt_conf_9[3]*dv11; 
  outb3[7] += 1.224744871391589*pt_conf_9[1]*vmap_vx[1]*hamil[2]*dv11; 
  outb3[8] += 1.224744871391589*vmap_vx[1]*hamil[2]*pt_conf_9[2]*dv11; 
  outb3[9] += 2.7386127875258306*vmap_vx[0]*pt_conf_9[1]*hamil[8]*dv11; 
  outb3[10] += 2.7386127875258306*vmap_vx[0]*pt_conf_9[2]*hamil[8]*dv11; 
  outb3[11] += 2.7386127875258306*pt_conf_9[0]*vmap_vx[1]*hamil[8]*dv11; 
  outb3[16] += 1.224744871391589*pt_conf_9[0]*hamil[2]*vmap_vx[2]*dv11; 
  outb3[19] += 1.224744871391589*vmap_vx[1]*hamil[2]*pt_conf_9[3]*dv11; 
  outb3[20] += 2.7386127875258306*vmap_vx[0]*pt_conf_9[3]*hamil[8]*dv11; 
  outb3[21] += 2.7386127875258306*pt_conf_9[1]*vmap_vx[1]*hamil[8]*dv11; 
  outb3[22] += 2.7386127875258306*vmap_vx[1]*pt_conf_9[2]*hamil[8]*dv11; 
  outb3[29] += 1.224744871391589*pt_conf_9[1]*hamil[2]*vmap_vx[2]*dv11; 
  outb3[30] += 1.224744871391589*hamil[2]*pt_conf_9[2]*vmap_vx[2]*dv11; 
  outb3[31] += 2.7386127875258306*pt_conf_9[0]*vmap_vx[2]*hamil[8]*dv11; 
  outb3[41] += 2.7386127875258306*vmap_vx[1]*pt_conf_9[3]*hamil[8]*dv11; 
  outb3[46] += 1.224744871391589*hamil[2]*vmap_vx[2]*pt_conf_9[3]*dv11; 
  outb3[47] += 2.7386127875258306*pt_conf_9[1]*vmap_vx[2]*hamil[8]*dv11; 
  outb3[48] += 2.7386127875258306*pt_conf_9[2]*vmap_vx[2]*hamil[8]*dv11; 
  outb3[68] += 2.7386127875258306*vmap_vx[2]*pt_conf_9[3]*hamil[8]*dv11; 
  } 
  if (dir == 3) { 
  double *outb2 = &out[216]; 
  outb2[0] += (-(2.7386127875258306*pt_conf_9[0]*vmap_vx[1]*hamil[7])-1.224744871391589*pt_conf_9[0]*vmap_vx[0]*hamil[1])*dv10; 
  outb2[1] += (-(2.7386127875258306*pt_conf_9[1]*vmap_vx[1]*hamil[7])-1.224744871391589*vmap_vx[0]*hamil[1]*pt_conf_9[1])*dv10; 
  outb2[2] += (-(2.7386127875258306*vmap_vx[1]*pt_conf_9[2]*hamil[7])-1.224744871391589*vmap_vx[0]*hamil[1]*pt_conf_9[2])*dv10; 
  outb2[3] += (-(2.4494897427831783*pt_conf_9[0]*vmap_vx[2]*hamil[7])-2.7386127875258306*pt_conf_9[0]*vmap_vx[0]*hamil[7]-1.224744871391589*pt_conf_9[0]*hamil[1]*vmap_vx[1])*dv10; 
  outb2[6] += (-(2.7386127875258306*vmap_vx[1]*pt_conf_9[3]*hamil[7])-1.224744871391589*vmap_vx[0]*hamil[1]*pt_conf_9[3])*dv10; 
  outb2[7] += (-(2.4494897427831783*pt_conf_9[1]*vmap_vx[2]*hamil[7])-2.7386127875258306*vmap_vx[0]*pt_conf_9[1]*hamil[7]-1.224744871391589*hamil[1]*pt_conf_9[1]*vmap_vx[1])*dv10; 
  outb2[8] += (-(2.4494897427831783*pt_conf_9[2]*vmap_vx[2]*hamil[7])-2.7386127875258306*vmap_vx[0]*pt_conf_9[2]*hamil[7]-1.224744871391589*hamil[1]*vmap_vx[1]*pt_conf_9[2])*dv10; 
  outb2[16] += (-(2.4494897427831783*pt_conf_9[0]*vmap_vx[1]*hamil[7])-1.224744871391589*pt_conf_9[0]*hamil[1]*vmap_vx[2])*dv10; 
  outb2[19] += (-(2.4494897427831783*vmap_vx[2]*pt_conf_9[3]*hamil[7])-2.7386127875258306*vmap_vx[0]*pt_conf_9[3]*hamil[7]-1.224744871391589*hamil[1]*vmap_vx[1]*pt_conf_9[3])*dv10; 
  outb2[29] += (-(2.4494897427831783*pt_conf_9[1]*vmap_vx[1]*hamil[7])-1.224744871391589*hamil[1]*pt_conf_9[1]*vmap_vx[2])*dv10; 
  outb2[30] += (-(2.4494897427831783*vmap_vx[1]*pt_conf_9[2]*hamil[7])-1.224744871391589*hamil[1]*pt_conf_9[2]*vmap_vx[2])*dv10; 
  outb2[46] += (-(2.4494897427831783*vmap_vx[1]*pt_conf_9[3]*hamil[7])-1.224744871391589*hamil[1]*vmap_vx[2]*pt_conf_9[3])*dv10; 
  } 

  const double *pt_conf_10 = &poisson_tensor_conf[40]; 
  if (dir == 2) { 
  double *outb4 = &out[432]; 
  outb4[0] += 1.224744871391589*pt_conf_10[0]*vmap_vx[0]*hamil[3]*dv12; 
  outb4[1] += 1.224744871391589*vmap_vx[0]*pt_conf_10[1]*hamil[3]*dv12; 
  outb4[2] += 1.224744871391589*vmap_vx[0]*pt_conf_10[2]*hamil[3]*dv12; 
  outb4[3] += 1.224744871391589*pt_conf_10[0]*vmap_vx[1]*hamil[3]*dv12; 
  outb4[5] += 2.7386127875258306*pt_conf_10[0]*vmap_vx[0]*hamil[9]*dv12; 
  outb4[6] += 1.224744871391589*vmap_vx[0]*hamil[3]*pt_conf_10[3]*dv12; 
  outb4[7] += 1.224744871391589*pt_conf_10[1]*vmap_vx[1]*hamil[3]*dv12; 
  outb4[8] += 1.224744871391589*vmap_vx[1]*pt_conf_10[2]*hamil[3]*dv12; 
  outb4[12] += 2.7386127875258306*vmap_vx[0]*pt_conf_10[1]*hamil[9]*dv12; 
  outb4[13] += 2.7386127875258306*vmap_vx[0]*pt_conf_10[2]*hamil[9]*dv12; 
  outb4[14] += 2.7386127875258306*pt_conf_10[0]*vmap_vx[1]*hamil[9]*dv12; 
  outb4[16] += 1.224744871391589*pt_conf_10[0]*vmap_vx[2]*hamil[3]*dv12; 
  outb4[19] += 1.224744871391589*vmap_vx[1]*hamil[3]*pt_conf_10[3]*dv12; 
  outb4[23] += 2.7386127875258306*vmap_vx[0]*pt_conf_10[3]*hamil[9]*dv12; 
  outb4[24] += 2.7386127875258306*pt_conf_10[1]*vmap_vx[1]*hamil[9]*dv12; 
  outb4[25] += 2.7386127875258306*vmap_vx[1]*pt_conf_10[2]*hamil[9]*dv12; 
  outb4[29] += 1.224744871391589*pt_conf_10[1]*vmap_vx[2]*hamil[3]*dv12; 
  outb4[30] += 1.224744871391589*pt_conf_10[2]*vmap_vx[2]*hamil[3]*dv12; 
  outb4[35] += 2.7386127875258306*pt_conf_10[0]*vmap_vx[2]*hamil[9]*dv12; 
  outb4[42] += 2.7386127875258306*vmap_vx[1]*pt_conf_10[3]*hamil[9]*dv12; 
  outb4[46] += 1.224744871391589*vmap_vx[2]*hamil[3]*pt_conf_10[3]*dv12; 
  outb4[52] += 2.7386127875258306*pt_conf_10[1]*vmap_vx[2]*hamil[9]*dv12; 
  outb4[53] += 2.7386127875258306*pt_conf_10[2]*vmap_vx[2]*hamil[9]*dv12; 
  outb4[70] += 2.7386127875258306*vmap_vx[2]*pt_conf_10[3]*hamil[9]*dv12; 
  } 
  if (dir == 4) { 
  double *outb2 = &out[216]; 
  outb2[0] += (-(2.7386127875258306*pt_conf_10[0]*vmap_vx[1]*hamil[7])-1.224744871391589*pt_conf_10[0]*vmap_vx[0]*hamil[1])*dv10; 
  outb2[1] += (-(2.7386127875258306*pt_conf_10[1]*vmap_vx[1]*hamil[7])-1.224744871391589*vmap_vx[0]*hamil[1]*pt_conf_10[1])*dv10; 
  outb2[2] += (-(2.7386127875258306*vmap_vx[1]*pt_conf_10[2]*hamil[7])-1.224744871391589*vmap_vx[0]*hamil[1]*pt_conf_10[2])*dv10; 
  outb2[3] += (-(2.4494897427831783*pt_conf_10[0]*vmap_vx[2]*hamil[7])-2.7386127875258306*pt_conf_10[0]*vmap_vx[0]*hamil[7]-1.224744871391589*pt_conf_10[0]*hamil[1]*vmap_vx[1])*dv10; 
  outb2[6] += (-(2.7386127875258306*vmap_vx[1]*pt_conf_10[3]*hamil[7])-1.224744871391589*vmap_vx[0]*hamil[1]*pt_conf_10[3])*dv10; 
  outb2[7] += (-(2.4494897427831783*pt_conf_10[1]*vmap_vx[2]*hamil[7])-2.7386127875258306*vmap_vx[0]*pt_conf_10[1]*hamil[7]-1.224744871391589*hamil[1]*pt_conf_10[1]*vmap_vx[1])*dv10; 
  outb2[8] += (-(2.4494897427831783*pt_conf_10[2]*vmap_vx[2]*hamil[7])-2.7386127875258306*vmap_vx[0]*pt_conf_10[2]*hamil[7]-1.224744871391589*hamil[1]*vmap_vx[1]*pt_conf_10[2])*dv10; 
  outb2[16] += (-(2.4494897427831783*pt_conf_10[0]*vmap_vx[1]*hamil[7])-1.224744871391589*pt_conf_10[0]*hamil[1]*vmap_vx[2])*dv10; 
  outb2[19] += (-(2.4494897427831783*vmap_vx[2]*pt_conf_10[3]*hamil[7])-2.7386127875258306*vmap_vx[0]*pt_conf_10[3]*hamil[7]-1.224744871391589*hamil[1]*vmap_vx[1]*pt_conf_10[3])*dv10; 
  outb2[29] += (-(2.4494897427831783*pt_conf_10[1]*vmap_vx[1]*hamil[7])-1.224744871391589*hamil[1]*pt_conf_10[1]*vmap_vx[2])*dv10; 
  outb2[30] += (-(2.4494897427831783*vmap_vx[1]*pt_conf_10[2]*hamil[7])-1.224744871391589*hamil[1]*pt_conf_10[2]*vmap_vx[2])*dv10; 
  outb2[46] += (-(2.4494897427831783*vmap_vx[1]*pt_conf_10[3]*hamil[7])-1.224744871391589*hamil[1]*vmap_vx[2]*pt_conf_10[3])*dv10; 
  } 

  const double *pt_conf_11 = &poisson_tensor_conf[44]; 
  if (dir == 3) { 
  double *outb4 = &out[432]; 
  outb4[0] += 1.224744871391589*pt_conf_11[0]*vmap_vx[0]*hamil[3]*dv12; 
  outb4[1] += 1.224744871391589*vmap_vx[0]*pt_conf_11[1]*hamil[3]*dv12; 
  outb4[2] += 1.224744871391589*vmap_vx[0]*pt_conf_11[2]*hamil[3]*dv12; 
  outb4[3] += 1.224744871391589*pt_conf_11[0]*vmap_vx[1]*hamil[3]*dv12; 
  outb4[5] += 2.7386127875258306*pt_conf_11[0]*vmap_vx[0]*hamil[9]*dv12; 
  outb4[6] += 1.224744871391589*vmap_vx[0]*hamil[3]*pt_conf_11[3]*dv12; 
  outb4[7] += 1.224744871391589*pt_conf_11[1]*vmap_vx[1]*hamil[3]*dv12; 
  outb4[8] += 1.224744871391589*vmap_vx[1]*pt_conf_11[2]*hamil[3]*dv12; 
  outb4[12] += 2.7386127875258306*vmap_vx[0]*pt_conf_11[1]*hamil[9]*dv12; 
  outb4[13] += 2.7386127875258306*vmap_vx[0]*pt_conf_11[2]*hamil[9]*dv12; 
  outb4[14] += 2.7386127875258306*pt_conf_11[0]*vmap_vx[1]*hamil[9]*dv12; 
  outb4[16] += 1.224744871391589*pt_conf_11[0]*vmap_vx[2]*hamil[3]*dv12; 
  outb4[19] += 1.224744871391589*vmap_vx[1]*hamil[3]*pt_conf_11[3]*dv12; 
  outb4[23] += 2.7386127875258306*vmap_vx[0]*pt_conf_11[3]*hamil[9]*dv12; 
  outb4[24] += 2.7386127875258306*pt_conf_11[1]*vmap_vx[1]*hamil[9]*dv12; 
  outb4[25] += 2.7386127875258306*vmap_vx[1]*pt_conf_11[2]*hamil[9]*dv12; 
  outb4[29] += 1.224744871391589*pt_conf_11[1]*vmap_vx[2]*hamil[3]*dv12; 
  outb4[30] += 1.224744871391589*pt_conf_11[2]*vmap_vx[2]*hamil[3]*dv12; 
  outb4[35] += 2.7386127875258306*pt_conf_11[0]*vmap_vx[2]*hamil[9]*dv12; 
  outb4[42] += 2.7386127875258306*vmap_vx[1]*pt_conf_11[3]*hamil[9]*dv12; 
  outb4[46] += 1.224744871391589*vmap_vx[2]*hamil[3]*pt_conf_11[3]*dv12; 
  outb4[52] += 2.7386127875258306*pt_conf_11[1]*vmap_vx[2]*hamil[9]*dv12; 
  outb4[53] += 2.7386127875258306*pt_conf_11[2]*vmap_vx[2]*hamil[9]*dv12; 
  outb4[70] += 2.7386127875258306*vmap_vx[2]*pt_conf_11[3]*hamil[9]*dv12; 
  } 
  if (dir == 4) { 
  double *outb3 = &out[324]; 
  outb3[0] += -(1.224744871391589*pt_conf_11[0]*vmap_vx[0]*hamil[2]*dv11); 
  outb3[1] += -(1.224744871391589*vmap_vx[0]*pt_conf_11[1]*hamil[2]*dv11); 
  outb3[2] += -(1.224744871391589*vmap_vx[0]*hamil[2]*pt_conf_11[2]*dv11); 
  outb3[3] += -(1.224744871391589*pt_conf_11[0]*vmap_vx[1]*hamil[2]*dv11); 
  outb3[4] += -(2.7386127875258306*pt_conf_11[0]*vmap_vx[0]*hamil[8]*dv11); 
  outb3[6] += -(1.224744871391589*vmap_vx[0]*hamil[2]*pt_conf_11[3]*dv11); 
  outb3[7] += -(1.224744871391589*pt_conf_11[1]*vmap_vx[1]*hamil[2]*dv11); 
  outb3[8] += -(1.224744871391589*vmap_vx[1]*hamil[2]*pt_conf_11[2]*dv11); 
  outb3[9] += -(2.7386127875258306*vmap_vx[0]*pt_conf_11[1]*hamil[8]*dv11); 
  outb3[10] += -(2.7386127875258306*vmap_vx[0]*pt_conf_11[2]*hamil[8]*dv11); 
  outb3[11] += -(2.7386127875258306*pt_conf_11[0]*vmap_vx[1]*hamil[8]*dv11); 
  outb3[16] += -(1.224744871391589*pt_conf_11[0]*hamil[2]*vmap_vx[2]*dv11); 
  outb3[19] += -(1.224744871391589*vmap_vx[1]*hamil[2]*pt_conf_11[3]*dv11); 
  outb3[20] += -(2.7386127875258306*vmap_vx[0]*pt_conf_11[3]*hamil[8]*dv11); 
  outb3[21] += -(2.7386127875258306*pt_conf_11[1]*vmap_vx[1]*hamil[8]*dv11); 
  outb3[22] += -(2.7386127875258306*vmap_vx[1]*pt_conf_11[2]*hamil[8]*dv11); 
  outb3[29] += -(1.224744871391589*pt_conf_11[1]*hamil[2]*vmap_vx[2]*dv11); 
  outb3[30] += -(1.224744871391589*hamil[2]*pt_conf_11[2]*vmap_vx[2]*dv11); 
  outb3[31] += -(2.7386127875258306*pt_conf_11[0]*vmap_vx[2]*hamil[8]*dv11); 
  outb3[41] += -(2.7386127875258306*vmap_vx[1]*pt_conf_11[3]*hamil[8]*dv11); 
  outb3[46] += -(1.224744871391589*hamil[2]*vmap_vx[2]*pt_conf_11[3]*dv11); 
  outb3[47] += -(2.7386127875258306*pt_conf_11[1]*vmap_vx[2]*hamil[8]*dv11); 
  outb3[48] += -(2.7386127875258306*pt_conf_11[2]*vmap_vx[2]*hamil[8]*dv11); 
  outb3[68] += -(2.7386127875258306*vmap_vx[2]*pt_conf_11[3]*hamil[8]*dv11); 
  } 

  const double *pt_conf_12 = &poisson_tensor_conf[48]; 
  if (dir == 2) { 
  double *outb3 = &out[324]; 
  outb3[0] += (2.7386127875258306*pt_conf_12[0]*vmap_vy[1]*hamil[8]+1.224744871391589*pt_conf_12[0]*vmap_vy[0]*hamil[2])*dv11; 
  outb3[1] += (2.7386127875258306*pt_conf_12[1]*vmap_vy[1]*hamil[8]+1.224744871391589*vmap_vy[0]*pt_conf_12[1]*hamil[2])*dv11; 
  outb3[2] += (2.7386127875258306*vmap_vy[1]*pt_conf_12[2]*hamil[8]+1.224744871391589*vmap_vy[0]*hamil[2]*pt_conf_12[2])*dv11; 
  outb3[4] += (2.4494897427831783*pt_conf_12[0]*vmap_vy[2]*hamil[8]+2.7386127875258306*pt_conf_12[0]*vmap_vy[0]*hamil[8]+1.224744871391589*pt_conf_12[0]*vmap_vy[1]*hamil[2])*dv11; 
  outb3[6] += (2.7386127875258306*vmap_vy[1]*pt_conf_12[3]*hamil[8]+1.224744871391589*vmap_vy[0]*hamil[2]*pt_conf_12[3])*dv11; 
  outb3[9] += (2.4494897427831783*pt_conf_12[1]*vmap_vy[2]*hamil[8]+2.7386127875258306*vmap_vy[0]*pt_conf_12[1]*hamil[8]+1.224744871391589*pt_conf_12[1]*vmap_vy[1]*hamil[2])*dv11; 
  outb3[10] += (2.4494897427831783*pt_conf_12[2]*vmap_vy[2]*hamil[8]+2.7386127875258306*vmap_vy[0]*pt_conf_12[2]*hamil[8]+1.224744871391589*vmap_vy[1]*hamil[2]*pt_conf_12[2])*dv11; 
  outb3[17] += (2.4494897427831783*pt_conf_12[0]*vmap_vy[1]*hamil[8]+1.224744871391589*pt_conf_12[0]*hamil[2]*vmap_vy[2])*dv11; 
  outb3[20] += (2.4494897427831783*vmap_vy[2]*pt_conf_12[3]*hamil[8]+2.7386127875258306*vmap_vy[0]*pt_conf_12[3]*hamil[8]+1.224744871391589*vmap_vy[1]*hamil[2]*pt_conf_12[3])*dv11; 
  outb3[32] += (2.4494897427831783*pt_conf_12[1]*vmap_vy[1]*hamil[8]+1.224744871391589*pt_conf_12[1]*hamil[2]*vmap_vy[2])*dv11; 
  outb3[33] += (2.4494897427831783*vmap_vy[1]*pt_conf_12[2]*hamil[8]+1.224744871391589*hamil[2]*pt_conf_12[2]*vmap_vy[2])*dv11; 
  outb3[49] += (2.4494897427831783*vmap_vy[1]*pt_conf_12[3]*hamil[8]+1.224744871391589*hamil[2]*vmap_vy[2]*pt_conf_12[3])*dv11; 
  } 
  if (dir == 3) { 
  double *outb2 = &out[216]; 
  outb2[0] += -(1.224744871391589*pt_conf_12[0]*vmap_vy[0]*hamil[1]*dv10); 
  outb2[1] += -(1.224744871391589*vmap_vy[0]*hamil[1]*pt_conf_12[1]*dv10); 
  outb2[2] += -(1.224744871391589*vmap_vy[0]*hamil[1]*pt_conf_12[2]*dv10); 
  outb2[3] += -(2.7386127875258306*pt_conf_12[0]*vmap_vy[0]*hamil[7]*dv10); 
  outb2[4] += -(1.224744871391589*pt_conf_12[0]*hamil[1]*vmap_vy[1]*dv10); 
  outb2[6] += -(1.224744871391589*vmap_vy[0]*hamil[1]*pt_conf_12[3]*dv10); 
  outb2[7] += -(2.7386127875258306*vmap_vy[0]*pt_conf_12[1]*hamil[7]*dv10); 
  outb2[8] += -(2.7386127875258306*vmap_vy[0]*pt_conf_12[2]*hamil[7]*dv10); 
  outb2[9] += -(1.224744871391589*hamil[1]*pt_conf_12[1]*vmap_vy[1]*dv10); 
  outb2[10] += -(1.224744871391589*hamil[1]*vmap_vy[1]*pt_conf_12[2]*dv10); 
  outb2[11] += -(2.7386127875258306*pt_conf_12[0]*vmap_vy[1]*hamil[7]*dv10); 
  outb2[17] += -(1.224744871391589*pt_conf_12[0]*hamil[1]*vmap_vy[2]*dv10); 
  outb2[19] += -(2.7386127875258306*vmap_vy[0]*pt_conf_12[3]*hamil[7]*dv10); 
  outb2[20] += -(1.224744871391589*hamil[1]*vmap_vy[1]*pt_conf_12[3]*dv10); 
  outb2[21] += -(2.7386127875258306*pt_conf_12[1]*vmap_vy[1]*hamil[7]*dv10); 
  outb2[22] += -(2.7386127875258306*vmap_vy[1]*pt_conf_12[2]*hamil[7]*dv10); 
  outb2[32] += -(1.224744871391589*hamil[1]*pt_conf_12[1]*vmap_vy[2]*dv10); 
  outb2[33] += -(1.224744871391589*hamil[1]*pt_conf_12[2]*vmap_vy[2]*dv10); 
  outb2[34] += -(2.7386127875258306*pt_conf_12[0]*vmap_vy[2]*hamil[7]*dv10); 
  outb2[41] += -(2.7386127875258306*vmap_vy[1]*pt_conf_12[3]*hamil[7]*dv10); 
  outb2[49] += -(1.224744871391589*hamil[1]*vmap_vy[2]*pt_conf_12[3]*dv10); 
  outb2[50] += -(2.7386127875258306*pt_conf_12[1]*vmap_vy[2]*hamil[7]*dv10); 
  outb2[51] += -(2.7386127875258306*pt_conf_12[2]*vmap_vy[2]*hamil[7]*dv10); 
  outb2[69] += -(2.7386127875258306*vmap_vy[2]*pt_conf_12[3]*hamil[7]*dv10); 
  } 

  const double *pt_conf_13 = &poisson_tensor_conf[52]; 
  if (dir == 2) { 
  double *outb4 = &out[432]; 
  outb4[0] += 1.224744871391589*pt_conf_13[0]*vmap_vy[0]*hamil[3]*dv12; 
  outb4[1] += 1.224744871391589*vmap_vy[0]*pt_conf_13[1]*hamil[3]*dv12; 
  outb4[2] += 1.224744871391589*vmap_vy[0]*pt_conf_13[2]*hamil[3]*dv12; 
  outb4[4] += 1.224744871391589*pt_conf_13[0]*vmap_vy[1]*hamil[3]*dv12; 
  outb4[5] += 2.7386127875258306*pt_conf_13[0]*vmap_vy[0]*hamil[9]*dv12; 
  outb4[6] += 1.224744871391589*vmap_vy[0]*hamil[3]*pt_conf_13[3]*dv12; 
  outb4[9] += 1.224744871391589*pt_conf_13[1]*vmap_vy[1]*hamil[3]*dv12; 
  outb4[10] += 1.224744871391589*vmap_vy[1]*pt_conf_13[2]*hamil[3]*dv12; 
  outb4[12] += 2.7386127875258306*vmap_vy[0]*pt_conf_13[1]*hamil[9]*dv12; 
  outb4[13] += 2.7386127875258306*vmap_vy[0]*pt_conf_13[2]*hamil[9]*dv12; 
  outb4[15] += 2.7386127875258306*pt_conf_13[0]*vmap_vy[1]*hamil[9]*dv12; 
  outb4[17] += 1.224744871391589*pt_conf_13[0]*vmap_vy[2]*hamil[3]*dv12; 
  outb4[20] += 1.224744871391589*vmap_vy[1]*hamil[3]*pt_conf_13[3]*dv12; 
  outb4[23] += 2.7386127875258306*vmap_vy[0]*pt_conf_13[3]*hamil[9]*dv12; 
  outb4[26] += 2.7386127875258306*pt_conf_13[1]*vmap_vy[1]*hamil[9]*dv12; 
  outb4[27] += 2.7386127875258306*vmap_vy[1]*pt_conf_13[2]*hamil[9]*dv12; 
  outb4[32] += 1.224744871391589*pt_conf_13[1]*vmap_vy[2]*hamil[3]*dv12; 
  outb4[33] += 1.224744871391589*pt_conf_13[2]*vmap_vy[2]*hamil[3]*dv12; 
  outb4[36] += 2.7386127875258306*pt_conf_13[0]*vmap_vy[2]*hamil[9]*dv12; 
  outb4[43] += 2.7386127875258306*vmap_vy[1]*pt_conf_13[3]*hamil[9]*dv12; 
  outb4[49] += 1.224744871391589*vmap_vy[2]*hamil[3]*pt_conf_13[3]*dv12; 
  outb4[55] += 2.7386127875258306*pt_conf_13[1]*vmap_vy[2]*hamil[9]*dv12; 
  outb4[56] += 2.7386127875258306*pt_conf_13[2]*vmap_vy[2]*hamil[9]*dv12; 
  outb4[73] += 2.7386127875258306*vmap_vy[2]*pt_conf_13[3]*hamil[9]*dv12; 
  } 
  if (dir == 4) { 
  double *outb2 = &out[216]; 
  outb2[0] += -(1.224744871391589*pt_conf_13[0]*vmap_vy[0]*hamil[1]*dv10); 
  outb2[1] += -(1.224744871391589*vmap_vy[0]*hamil[1]*pt_conf_13[1]*dv10); 
  outb2[2] += -(1.224744871391589*vmap_vy[0]*hamil[1]*pt_conf_13[2]*dv10); 
  outb2[3] += -(2.7386127875258306*pt_conf_13[0]*vmap_vy[0]*hamil[7]*dv10); 
  outb2[4] += -(1.224744871391589*pt_conf_13[0]*hamil[1]*vmap_vy[1]*dv10); 
  outb2[6] += -(1.224744871391589*vmap_vy[0]*hamil[1]*pt_conf_13[3]*dv10); 
  outb2[7] += -(2.7386127875258306*vmap_vy[0]*pt_conf_13[1]*hamil[7]*dv10); 
  outb2[8] += -(2.7386127875258306*vmap_vy[0]*pt_conf_13[2]*hamil[7]*dv10); 
  outb2[9] += -(1.224744871391589*hamil[1]*pt_conf_13[1]*vmap_vy[1]*dv10); 
  outb2[10] += -(1.224744871391589*hamil[1]*vmap_vy[1]*pt_conf_13[2]*dv10); 
  outb2[11] += -(2.7386127875258306*pt_conf_13[0]*vmap_vy[1]*hamil[7]*dv10); 
  outb2[17] += -(1.224744871391589*pt_conf_13[0]*hamil[1]*vmap_vy[2]*dv10); 
  outb2[19] += -(2.7386127875258306*vmap_vy[0]*pt_conf_13[3]*hamil[7]*dv10); 
  outb2[20] += -(1.224744871391589*hamil[1]*vmap_vy[1]*pt_conf_13[3]*dv10); 
  outb2[21] += -(2.7386127875258306*pt_conf_13[1]*vmap_vy[1]*hamil[7]*dv10); 
  outb2[22] += -(2.7386127875258306*vmap_vy[1]*pt_conf_13[2]*hamil[7]*dv10); 
  outb2[32] += -(1.224744871391589*hamil[1]*pt_conf_13[1]*vmap_vy[2]*dv10); 
  outb2[33] += -(1.224744871391589*hamil[1]*pt_conf_13[2]*vmap_vy[2]*dv10); 
  outb2[34] += -(2.7386127875258306*pt_conf_13[0]*vmap_vy[2]*hamil[7]*dv10); 
  outb2[41] += -(2.7386127875258306*vmap_vy[1]*pt_conf_13[3]*hamil[7]*dv10); 
  outb2[49] += -(1.224744871391589*hamil[1]*vmap_vy[2]*pt_conf_13[3]*dv10); 
  outb2[50] += -(2.7386127875258306*pt_conf_13[1]*vmap_vy[2]*hamil[7]*dv10); 
  outb2[51] += -(2.7386127875258306*pt_conf_13[2]*vmap_vy[2]*hamil[7]*dv10); 
  outb2[69] += -(2.7386127875258306*vmap_vy[2]*pt_conf_13[3]*hamil[7]*dv10); 
  } 

  const double *pt_conf_14 = &poisson_tensor_conf[56]; 
  if (dir == 3) { 
  double *outb4 = &out[432]; 
  outb4[0] += 1.224744871391589*pt_conf_14[0]*vmap_vy[0]*hamil[3]*dv12; 
  outb4[1] += 1.224744871391589*vmap_vy[0]*pt_conf_14[1]*hamil[3]*dv12; 
  outb4[2] += 1.224744871391589*vmap_vy[0]*pt_conf_14[2]*hamil[3]*dv12; 
  outb4[4] += 1.224744871391589*pt_conf_14[0]*vmap_vy[1]*hamil[3]*dv12; 
  outb4[5] += 2.7386127875258306*pt_conf_14[0]*vmap_vy[0]*hamil[9]*dv12; 
  outb4[6] += 1.224744871391589*vmap_vy[0]*hamil[3]*pt_conf_14[3]*dv12; 
  outb4[9] += 1.224744871391589*pt_conf_14[1]*vmap_vy[1]*hamil[3]*dv12; 
  outb4[10] += 1.224744871391589*vmap_vy[1]*pt_conf_14[2]*hamil[3]*dv12; 
  outb4[12] += 2.7386127875258306*vmap_vy[0]*pt_conf_14[1]*hamil[9]*dv12; 
  outb4[13] += 2.7386127875258306*vmap_vy[0]*pt_conf_14[2]*hamil[9]*dv12; 
  outb4[15] += 2.7386127875258306*pt_conf_14[0]*vmap_vy[1]*hamil[9]*dv12; 
  outb4[17] += 1.224744871391589*pt_conf_14[0]*vmap_vy[2]*hamil[3]*dv12; 
  outb4[20] += 1.224744871391589*vmap_vy[1]*hamil[3]*pt_conf_14[3]*dv12; 
  outb4[23] += 2.7386127875258306*vmap_vy[0]*pt_conf_14[3]*hamil[9]*dv12; 
  outb4[26] += 2.7386127875258306*pt_conf_14[1]*vmap_vy[1]*hamil[9]*dv12; 
  outb4[27] += 2.7386127875258306*vmap_vy[1]*pt_conf_14[2]*hamil[9]*dv12; 
  outb4[32] += 1.224744871391589*pt_conf_14[1]*vmap_vy[2]*hamil[3]*dv12; 
  outb4[33] += 1.224744871391589*pt_conf_14[2]*vmap_vy[2]*hamil[3]*dv12; 
  outb4[36] += 2.7386127875258306*pt_conf_14[0]*vmap_vy[2]*hamil[9]*dv12; 
  outb4[43] += 2.7386127875258306*vmap_vy[1]*pt_conf_14[3]*hamil[9]*dv12; 
  outb4[49] += 1.224744871391589*vmap_vy[2]*hamil[3]*pt_conf_14[3]*dv12; 
  outb4[55] += 2.7386127875258306*pt_conf_14[1]*vmap_vy[2]*hamil[9]*dv12; 
  outb4[56] += 2.7386127875258306*pt_conf_14[2]*vmap_vy[2]*hamil[9]*dv12; 
  outb4[73] += 2.7386127875258306*vmap_vy[2]*pt_conf_14[3]*hamil[9]*dv12; 
  } 
  if (dir == 4) { 
  double *outb3 = &out[324]; 
  outb3[0] += (-(2.7386127875258306*pt_conf_14[0]*vmap_vy[1]*hamil[8])-1.224744871391589*pt_conf_14[0]*vmap_vy[0]*hamil[2])*dv11; 
  outb3[1] += (-(2.7386127875258306*pt_conf_14[1]*vmap_vy[1]*hamil[8])-1.224744871391589*vmap_vy[0]*pt_conf_14[1]*hamil[2])*dv11; 
  outb3[2] += (-(2.7386127875258306*vmap_vy[1]*pt_conf_14[2]*hamil[8])-1.224744871391589*vmap_vy[0]*hamil[2]*pt_conf_14[2])*dv11; 
  outb3[4] += (-(2.4494897427831783*pt_conf_14[0]*vmap_vy[2]*hamil[8])-2.7386127875258306*pt_conf_14[0]*vmap_vy[0]*hamil[8]-1.224744871391589*pt_conf_14[0]*vmap_vy[1]*hamil[2])*dv11; 
  outb3[6] += (-(2.7386127875258306*vmap_vy[1]*pt_conf_14[3]*hamil[8])-1.224744871391589*vmap_vy[0]*hamil[2]*pt_conf_14[3])*dv11; 
  outb3[9] += (-(2.4494897427831783*pt_conf_14[1]*vmap_vy[2]*hamil[8])-2.7386127875258306*vmap_vy[0]*pt_conf_14[1]*hamil[8]-1.224744871391589*pt_conf_14[1]*vmap_vy[1]*hamil[2])*dv11; 
  outb3[10] += (-(2.4494897427831783*pt_conf_14[2]*vmap_vy[2]*hamil[8])-2.7386127875258306*vmap_vy[0]*pt_conf_14[2]*hamil[8]-1.224744871391589*vmap_vy[1]*hamil[2]*pt_conf_14[2])*dv11; 
  outb3[17] += (-(2.4494897427831783*pt_conf_14[0]*vmap_vy[1]*hamil[8])-1.224744871391589*pt_conf_14[0]*hamil[2]*vmap_vy[2])*dv11; 
  outb3[20] += (-(2.4494897427831783*vmap_vy[2]*pt_conf_14[3]*hamil[8])-2.7386127875258306*vmap_vy[0]*pt_conf_14[3]*hamil[8]-1.224744871391589*vmap_vy[1]*hamil[2]*pt_conf_14[3])*dv11; 
  outb3[32] += (-(2.4494897427831783*pt_conf_14[1]*vmap_vy[1]*hamil[8])-1.224744871391589*pt_conf_14[1]*hamil[2]*vmap_vy[2])*dv11; 
  outb3[33] += (-(2.4494897427831783*vmap_vy[1]*pt_conf_14[2]*hamil[8])-1.224744871391589*hamil[2]*pt_conf_14[2]*vmap_vy[2])*dv11; 
  outb3[49] += (-(2.4494897427831783*vmap_vy[1]*pt_conf_14[3]*hamil[8])-1.224744871391589*hamil[2]*vmap_vy[2]*pt_conf_14[3])*dv11; 
  } 

  const double *pt_conf_15 = &poisson_tensor_conf[60]; 
  if (dir == 2) { 
  double *outb3 = &out[324]; 
  outb3[0] += 1.224744871391589*pt_conf_15[0]*vmap_vz[0]*hamil[2]*dv11; 
  outb3[1] += 1.224744871391589*vmap_vz[0]*pt_conf_15[1]*hamil[2]*dv11; 
  outb3[2] += 1.224744871391589*vmap_vz[0]*hamil[2]*pt_conf_15[2]*dv11; 
  outb3[4] += 2.7386127875258306*pt_conf_15[0]*vmap_vz[0]*hamil[8]*dv11; 
  outb3[5] += 1.224744871391589*pt_conf_15[0]*vmap_vz[1]*hamil[2]*dv11; 
  outb3[6] += 1.224744871391589*vmap_vz[0]*hamil[2]*pt_conf_15[3]*dv11; 
  outb3[9] += 2.7386127875258306*vmap_vz[0]*pt_conf_15[1]*hamil[8]*dv11; 
  outb3[10] += 2.7386127875258306*vmap_vz[0]*pt_conf_15[2]*hamil[8]*dv11; 
  outb3[12] += 1.224744871391589*pt_conf_15[1]*vmap_vz[1]*hamil[2]*dv11; 
  outb3[13] += 1.224744871391589*vmap_vz[1]*hamil[2]*pt_conf_15[2]*dv11; 
  outb3[15] += 2.7386127875258306*pt_conf_15[0]*vmap_vz[1]*hamil[8]*dv11; 
  outb3[18] += 1.224744871391589*pt_conf_15[0]*hamil[2]*vmap_vz[2]*dv11; 
  outb3[20] += 2.7386127875258306*vmap_vz[0]*pt_conf_15[3]*hamil[8]*dv11; 
  outb3[23] += 1.224744871391589*vmap_vz[1]*hamil[2]*pt_conf_15[3]*dv11; 
  outb3[26] += 2.7386127875258306*pt_conf_15[1]*vmap_vz[1]*hamil[8]*dv11; 
  outb3[27] += 2.7386127875258306*vmap_vz[1]*pt_conf_15[2]*hamil[8]*dv11; 
  outb3[37] += 1.224744871391589*pt_conf_15[1]*hamil[2]*vmap_vz[2]*dv11; 
  outb3[38] += 1.224744871391589*hamil[2]*pt_conf_15[2]*vmap_vz[2]*dv11; 
  outb3[40] += 2.7386127875258306*pt_conf_15[0]*vmap_vz[2]*hamil[8]*dv11; 
  outb3[43] += 2.7386127875258306*vmap_vz[1]*pt_conf_15[3]*hamil[8]*dv11; 
  outb3[58] += 1.224744871391589*hamil[2]*vmap_vz[2]*pt_conf_15[3]*dv11; 
  outb3[61] += 2.7386127875258306*pt_conf_15[1]*vmap_vz[2]*hamil[8]*dv11; 
  outb3[62] += 2.7386127875258306*pt_conf_15[2]*vmap_vz[2]*hamil[8]*dv11; 
  outb3[77] += 2.7386127875258306*vmap_vz[2]*pt_conf_15[3]*hamil[8]*dv11; 
  } 
  if (dir == 3) { 
  double *outb2 = &out[216]; 
  outb2[0] += -(1.224744871391589*pt_conf_15[0]*vmap_vz[0]*hamil[1]*dv10); 
  outb2[1] += -(1.224744871391589*vmap_vz[0]*hamil[1]*pt_conf_15[1]*dv10); 
  outb2[2] += -(1.224744871391589*vmap_vz[0]*hamil[1]*pt_conf_15[2]*dv10); 
  outb2[3] += -(2.7386127875258306*pt_conf_15[0]*vmap_vz[0]*hamil[7]*dv10); 
  outb2[5] += -(1.224744871391589*pt_conf_15[0]*hamil[1]*vmap_vz[1]*dv10); 
  outb2[6] += -(1.224744871391589*vmap_vz[0]*hamil[1]*pt_conf_15[3]*dv10); 
  outb2[7] += -(2.7386127875258306*vmap_vz[0]*pt_conf_15[1]*hamil[7]*dv10); 
  outb2[8] += -(2.7386127875258306*vmap_vz[0]*pt_conf_15[2]*hamil[7]*dv10); 
  outb2[12] += -(1.224744871391589*hamil[1]*pt_conf_15[1]*vmap_vz[1]*dv10); 
  outb2[13] += -(1.224744871391589*hamil[1]*vmap_vz[1]*pt_conf_15[2]*dv10); 
  outb2[14] += -(2.7386127875258306*pt_conf_15[0]*vmap_vz[1]*hamil[7]*dv10); 
  outb2[18] += -(1.224744871391589*pt_conf_15[0]*hamil[1]*vmap_vz[2]*dv10); 
  outb2[19] += -(2.7386127875258306*vmap_vz[0]*pt_conf_15[3]*hamil[7]*dv10); 
  outb2[23] += -(1.224744871391589*hamil[1]*vmap_vz[1]*pt_conf_15[3]*dv10); 
  outb2[24] += -(2.7386127875258306*pt_conf_15[1]*vmap_vz[1]*hamil[7]*dv10); 
  outb2[25] += -(2.7386127875258306*vmap_vz[1]*pt_conf_15[2]*hamil[7]*dv10); 
  outb2[37] += -(1.224744871391589*hamil[1]*pt_conf_15[1]*vmap_vz[2]*dv10); 
  outb2[38] += -(1.224744871391589*hamil[1]*pt_conf_15[2]*vmap_vz[2]*dv10); 
  outb2[39] += -(2.7386127875258306*pt_conf_15[0]*vmap_vz[2]*hamil[7]*dv10); 
  outb2[42] += -(2.7386127875258306*vmap_vz[1]*pt_conf_15[3]*hamil[7]*dv10); 
  outb2[58] += -(1.224744871391589*hamil[1]*vmap_vz[2]*pt_conf_15[3]*dv10); 
  outb2[59] += -(2.7386127875258306*pt_conf_15[1]*vmap_vz[2]*hamil[7]*dv10); 
  outb2[60] += -(2.7386127875258306*pt_conf_15[2]*vmap_vz[2]*hamil[7]*dv10); 
  outb2[76] += -(2.7386127875258306*vmap_vz[2]*pt_conf_15[3]*hamil[7]*dv10); 
  } 

  const double *pt_conf_16 = &poisson_tensor_conf[64]; 
  if (dir == 2) { 
  double *outb4 = &out[432]; 
  outb4[0] += (2.7386127875258306*pt_conf_16[0]*vmap_vz[1]*hamil[9]+1.224744871391589*pt_conf_16[0]*vmap_vz[0]*hamil[3])*dv12; 
  outb4[1] += (2.7386127875258306*pt_conf_16[1]*vmap_vz[1]*hamil[9]+1.224744871391589*vmap_vz[0]*pt_conf_16[1]*hamil[3])*dv12; 
  outb4[2] += (2.7386127875258306*vmap_vz[1]*pt_conf_16[2]*hamil[9]+1.224744871391589*vmap_vz[0]*pt_conf_16[2]*hamil[3])*dv12; 
  outb4[5] += (2.4494897427831783*pt_conf_16[0]*vmap_vz[2]*hamil[9]+2.7386127875258306*pt_conf_16[0]*vmap_vz[0]*hamil[9]+1.224744871391589*pt_conf_16[0]*vmap_vz[1]*hamil[3])*dv12; 
  outb4[6] += (2.7386127875258306*vmap_vz[1]*pt_conf_16[3]*hamil[9]+1.224744871391589*vmap_vz[0]*hamil[3]*pt_conf_16[3])*dv12; 
  outb4[12] += (2.4494897427831783*pt_conf_16[1]*vmap_vz[2]*hamil[9]+2.7386127875258306*vmap_vz[0]*pt_conf_16[1]*hamil[9]+1.224744871391589*pt_conf_16[1]*vmap_vz[1]*hamil[3])*dv12; 
  outb4[13] += (2.4494897427831783*pt_conf_16[2]*vmap_vz[2]*hamil[9]+2.7386127875258306*vmap_vz[0]*pt_conf_16[2]*hamil[9]+1.224744871391589*vmap_vz[1]*pt_conf_16[2]*hamil[3])*dv12; 
  outb4[18] += (2.4494897427831783*pt_conf_16[0]*vmap_vz[1]*hamil[9]+1.224744871391589*pt_conf_16[0]*vmap_vz[2]*hamil[3])*dv12; 
  outb4[23] += (2.4494897427831783*vmap_vz[2]*pt_conf_16[3]*hamil[9]+2.7386127875258306*vmap_vz[0]*pt_conf_16[3]*hamil[9]+1.224744871391589*vmap_vz[1]*hamil[3]*pt_conf_16[3])*dv12; 
  outb4[37] += (2.4494897427831783*pt_conf_16[1]*vmap_vz[1]*hamil[9]+1.224744871391589*pt_conf_16[1]*vmap_vz[2]*hamil[3])*dv12; 
  outb4[38] += (2.4494897427831783*vmap_vz[1]*pt_conf_16[2]*hamil[9]+1.224744871391589*pt_conf_16[2]*vmap_vz[2]*hamil[3])*dv12; 
  outb4[58] += (2.4494897427831783*vmap_vz[1]*pt_conf_16[3]*hamil[9]+1.224744871391589*vmap_vz[2]*hamil[3]*pt_conf_16[3])*dv12; 
  } 
  if (dir == 4) { 
  double *outb2 = &out[216]; 
  outb2[0] += -(1.224744871391589*pt_conf_16[0]*vmap_vz[0]*hamil[1]*dv10); 
  outb2[1] += -(1.224744871391589*vmap_vz[0]*hamil[1]*pt_conf_16[1]*dv10); 
  outb2[2] += -(1.224744871391589*vmap_vz[0]*hamil[1]*pt_conf_16[2]*dv10); 
  outb2[3] += -(2.7386127875258306*pt_conf_16[0]*vmap_vz[0]*hamil[7]*dv10); 
  outb2[5] += -(1.224744871391589*pt_conf_16[0]*hamil[1]*vmap_vz[1]*dv10); 
  outb2[6] += -(1.224744871391589*vmap_vz[0]*hamil[1]*pt_conf_16[3]*dv10); 
  outb2[7] += -(2.7386127875258306*vmap_vz[0]*pt_conf_16[1]*hamil[7]*dv10); 
  outb2[8] += -(2.7386127875258306*vmap_vz[0]*pt_conf_16[2]*hamil[7]*dv10); 
  outb2[12] += -(1.224744871391589*hamil[1]*pt_conf_16[1]*vmap_vz[1]*dv10); 
  outb2[13] += -(1.224744871391589*hamil[1]*vmap_vz[1]*pt_conf_16[2]*dv10); 
  outb2[14] += -(2.7386127875258306*pt_conf_16[0]*vmap_vz[1]*hamil[7]*dv10); 
  outb2[18] += -(1.224744871391589*pt_conf_16[0]*hamil[1]*vmap_vz[2]*dv10); 
  outb2[19] += -(2.7386127875258306*vmap_vz[0]*pt_conf_16[3]*hamil[7]*dv10); 
  outb2[23] += -(1.224744871391589*hamil[1]*vmap_vz[1]*pt_conf_16[3]*dv10); 
  outb2[24] += -(2.7386127875258306*pt_conf_16[1]*vmap_vz[1]*hamil[7]*dv10); 
  outb2[25] += -(2.7386127875258306*vmap_vz[1]*pt_conf_16[2]*hamil[7]*dv10); 
  outb2[37] += -(1.224744871391589*hamil[1]*pt_conf_16[1]*vmap_vz[2]*dv10); 
  outb2[38] += -(1.224744871391589*hamil[1]*pt_conf_16[2]*vmap_vz[2]*dv10); 
  outb2[39] += -(2.7386127875258306*pt_conf_16[0]*vmap_vz[2]*hamil[7]*dv10); 
  outb2[42] += -(2.7386127875258306*vmap_vz[1]*pt_conf_16[3]*hamil[7]*dv10); 
  outb2[58] += -(1.224744871391589*hamil[1]*vmap_vz[2]*pt_conf_16[3]*dv10); 
  outb2[59] += -(2.7386127875258306*pt_conf_16[1]*vmap_vz[2]*hamil[7]*dv10); 
  outb2[60] += -(2.7386127875258306*pt_conf_16[2]*vmap_vz[2]*hamil[7]*dv10); 
  outb2[76] += -(2.7386127875258306*vmap_vz[2]*pt_conf_16[3]*hamil[7]*dv10); 
  } 

  const double *pt_conf_17 = &poisson_tensor_conf[68]; 
  if (dir == 3) { 
  double *outb4 = &out[432]; 
  outb4[0] += (2.7386127875258306*pt_conf_17[0]*vmap_vz[1]*hamil[9]+1.224744871391589*pt_conf_17[0]*vmap_vz[0]*hamil[3])*dv12; 
  outb4[1] += (2.7386127875258306*pt_conf_17[1]*vmap_vz[1]*hamil[9]+1.224744871391589*vmap_vz[0]*pt_conf_17[1]*hamil[3])*dv12; 
  outb4[2] += (2.7386127875258306*vmap_vz[1]*pt_conf_17[2]*hamil[9]+1.224744871391589*vmap_vz[0]*pt_conf_17[2]*hamil[3])*dv12; 
  outb4[5] += (2.4494897427831783*pt_conf_17[0]*vmap_vz[2]*hamil[9]+2.7386127875258306*pt_conf_17[0]*vmap_vz[0]*hamil[9]+1.224744871391589*pt_conf_17[0]*vmap_vz[1]*hamil[3])*dv12; 
  outb4[6] += (2.7386127875258306*vmap_vz[1]*pt_conf_17[3]*hamil[9]+1.224744871391589*vmap_vz[0]*hamil[3]*pt_conf_17[3])*dv12; 
  outb4[12] += (2.4494897427831783*pt_conf_17[1]*vmap_vz[2]*hamil[9]+2.7386127875258306*vmap_vz[0]*pt_conf_17[1]*hamil[9]+1.224744871391589*pt_conf_17[1]*vmap_vz[1]*hamil[3])*dv12; 
  outb4[13] += (2.4494897427831783*pt_conf_17[2]*vmap_vz[2]*hamil[9]+2.7386127875258306*vmap_vz[0]*pt_conf_17[2]*hamil[9]+1.224744871391589*vmap_vz[1]*pt_conf_17[2]*hamil[3])*dv12; 
  outb4[18] += (2.4494897427831783*pt_conf_17[0]*vmap_vz[1]*hamil[9]+1.224744871391589*pt_conf_17[0]*vmap_vz[2]*hamil[3])*dv12; 
  outb4[23] += (2.4494897427831783*vmap_vz[2]*pt_conf_17[3]*hamil[9]+2.7386127875258306*vmap_vz[0]*pt_conf_17[3]*hamil[9]+1.224744871391589*vmap_vz[1]*hamil[3]*pt_conf_17[3])*dv12; 
  outb4[37] += (2.4494897427831783*pt_conf_17[1]*vmap_vz[1]*hamil[9]+1.224744871391589*pt_conf_17[1]*vmap_vz[2]*hamil[3])*dv12; 
  outb4[38] += (2.4494897427831783*vmap_vz[1]*pt_conf_17[2]*hamil[9]+1.224744871391589*pt_conf_17[2]*vmap_vz[2]*hamil[3])*dv12; 
  outb4[58] += (2.4494897427831783*vmap_vz[1]*pt_conf_17[3]*hamil[9]+1.224744871391589*vmap_vz[2]*hamil[3]*pt_conf_17[3])*dv12; 
  } 
  if (dir == 4) { 
  double *outb3 = &out[324]; 
  outb3[0] += -(1.224744871391589*pt_conf_17[0]*vmap_vz[0]*hamil[2]*dv11); 
  outb3[1] += -(1.224744871391589*vmap_vz[0]*pt_conf_17[1]*hamil[2]*dv11); 
  outb3[2] += -(1.224744871391589*vmap_vz[0]*hamil[2]*pt_conf_17[2]*dv11); 
  outb3[4] += -(2.7386127875258306*pt_conf_17[0]*vmap_vz[0]*hamil[8]*dv11); 
  outb3[5] += -(1.224744871391589*pt_conf_17[0]*vmap_vz[1]*hamil[2]*dv11); 
  outb3[6] += -(1.224744871391589*vmap_vz[0]*hamil[2]*pt_conf_17[3]*dv11); 
  outb3[9] += -(2.7386127875258306*vmap_vz[0]*pt_conf_17[1]*hamil[8]*dv11); 
  outb3[10] += -(2.7386127875258306*vmap_vz[0]*pt_conf_17[2]*hamil[8]*dv11); 
  outb3[12] += -(1.224744871391589*pt_conf_17[1]*vmap_vz[1]*hamil[2]*dv11); 
  outb3[13] += -(1.224744871391589*vmap_vz[1]*hamil[2]*pt_conf_17[2]*dv11); 
  outb3[15] += -(2.7386127875258306*pt_conf_17[0]*vmap_vz[1]*hamil[8]*dv11); 
  outb3[18] += -(1.224744871391589*pt_conf_17[0]*hamil[2]*vmap_vz[2]*dv11); 
  outb3[20] += -(2.7386127875258306*vmap_vz[0]*pt_conf_17[3]*hamil[8]*dv11); 
  outb3[23] += -(1.224744871391589*vmap_vz[1]*hamil[2]*pt_conf_17[3]*dv11); 
  outb3[26] += -(2.7386127875258306*pt_conf_17[1]*vmap_vz[1]*hamil[8]*dv11); 
  outb3[27] += -(2.7386127875258306*vmap_vz[1]*pt_conf_17[2]*hamil[8]*dv11); 
  outb3[37] += -(1.224744871391589*pt_conf_17[1]*hamil[2]*vmap_vz[2]*dv11); 
  outb3[38] += -(1.224744871391589*hamil[2]*pt_conf_17[2]*vmap_vz[2]*dv11); 
  outb3[40] += -(2.7386127875258306*pt_conf_17[0]*vmap_vz[2]*hamil[8]*dv11); 
  outb3[43] += -(2.7386127875258306*vmap_vz[1]*pt_conf_17[3]*hamil[8]*dv11); 
  outb3[58] += -(1.224744871391589*hamil[2]*vmap_vz[2]*pt_conf_17[3]*dv11); 
  outb3[61] += -(2.7386127875258306*pt_conf_17[1]*vmap_vz[2]*hamil[8]*dv11); 
  outb3[62] += -(2.7386127875258306*pt_conf_17[2]*vmap_vz[2]*hamil[8]*dv11); 
  outb3[77] += -(2.7386127875258306*vmap_vz[2]*pt_conf_17[3]*hamil[8]*dv11); 
  } 

} 
