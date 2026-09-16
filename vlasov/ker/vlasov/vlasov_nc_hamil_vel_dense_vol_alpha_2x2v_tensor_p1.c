#include <gkyl_nc_hamil_gen_vol_comps_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_vel_dense_vol_alpha_2x2v_tensor_p1(const double *w, const double *dxv, const int dir,
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
  const double jacob_cy_inv = 1.0/jacob_pos[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double *jacob_vy = &jacob_vel[3]; 
  const double *vmap_vx = &vmap[0]; 
  const double *vmap_vy = &vmap[4]; 
  
  const double *pt_conf_0 = &poisson_tensor_conf[0]; 
  if (dir == 0) { 
  double *outb2 = &out[72]; 
  outb2[0] += 1.7320508075688772*pt_conf_0[0]*hamil[1]*dv10; 
  outb2[1] += 1.7320508075688772*hamil[1]*pt_conf_0[1]*dv10; 
  outb2[2] += 1.7320508075688772*hamil[1]*pt_conf_0[2]*dv10; 
  outb2[3] += 3.872983346207417*pt_conf_0[0]*hamil[4]*dv10; 
  outb2[4] += 1.7320508075688772*pt_conf_0[0]*hamil[3]*dv10; 
  outb2[5] += 1.7320508075688772*hamil[1]*pt_conf_0[3]*dv10; 
  outb2[6] += 3.872983346207417*pt_conf_0[1]*hamil[4]*dv10; 
  outb2[7] += 3.872983346207417*pt_conf_0[2]*hamil[4]*dv10; 
  outb2[8] += 1.7320508075688772*pt_conf_0[1]*hamil[3]*dv10; 
  outb2[9] += 1.7320508075688772*pt_conf_0[2]*hamil[3]*dv10; 
  outb2[10] += 3.872983346207417*pt_conf_0[0]*hamil[6]*dv10; 
  outb2[12] += 1.7320508075688774*pt_conf_0[0]*hamil[7]*dv10; 
  outb2[13] += 3.872983346207417*pt_conf_0[3]*hamil[4]*dv10; 
  outb2[14] += 1.7320508075688772*hamil[3]*pt_conf_0[3]*dv10; 
  outb2[15] += 3.872983346207417*pt_conf_0[1]*hamil[6]*dv10; 
  outb2[16] += 3.872983346207417*pt_conf_0[2]*hamil[6]*dv10; 
  outb2[20] += 1.7320508075688772*pt_conf_0[1]*hamil[7]*dv10; 
  outb2[21] += 1.7320508075688772*pt_conf_0[2]*hamil[7]*dv10; 
  outb2[22] += 3.872983346207417*pt_conf_0[0]*hamil[8]*dv10; 
  outb2[23] += 3.872983346207417*pt_conf_0[3]*hamil[6]*dv10; 
  outb2[27] += 1.7320508075688774*pt_conf_0[3]*hamil[7]*dv10; 
  outb2[28] += 3.872983346207417*pt_conf_0[1]*hamil[8]*dv10; 
  outb2[29] += 3.872983346207417*pt_conf_0[2]*hamil[8]*dv10; 
  outb2[32] += 3.872983346207417*pt_conf_0[3]*hamil[8]*dv10; 
  } 
  if (dir == 2) { 
  double *outb0 = &out[0]; 
  } 

  const double *pt_conf_1 = &poisson_tensor_conf[4]; 
  if (dir == 0) { 
  double *outb3 = &out[108]; 
  outb3[0] += 1.7320508075688772*pt_conf_1[0]*hamil[2]*dv11; 
  outb3[1] += 1.7320508075688772*pt_conf_1[1]*hamil[2]*dv11; 
  outb3[2] += 1.7320508075688772*hamil[2]*pt_conf_1[2]*dv11; 
  outb3[3] += 1.7320508075688772*pt_conf_1[0]*hamil[3]*dv11; 
  outb3[4] += 3.872983346207417*pt_conf_1[0]*hamil[5]*dv11; 
  outb3[5] += 1.7320508075688772*hamil[2]*pt_conf_1[3]*dv11; 
  outb3[6] += 1.7320508075688772*pt_conf_1[1]*hamil[3]*dv11; 
  outb3[7] += 1.7320508075688772*pt_conf_1[2]*hamil[3]*dv11; 
  outb3[8] += 3.872983346207417*pt_conf_1[1]*hamil[5]*dv11; 
  outb3[9] += 3.872983346207417*pt_conf_1[2]*hamil[5]*dv11; 
  outb3[10] += 3.872983346207417*pt_conf_1[0]*hamil[7]*dv11; 
  outb3[11] += 1.7320508075688774*pt_conf_1[0]*hamil[6]*dv11; 
  outb3[13] += 1.7320508075688772*hamil[3]*pt_conf_1[3]*dv11; 
  outb3[14] += 3.872983346207417*pt_conf_1[3]*hamil[5]*dv11; 
  outb3[15] += 3.872983346207417*pt_conf_1[1]*hamil[7]*dv11; 
  outb3[16] += 3.872983346207417*pt_conf_1[2]*hamil[7]*dv11; 
  outb3[17] += 1.7320508075688772*pt_conf_1[1]*hamil[6]*dv11; 
  outb3[18] += 1.7320508075688772*pt_conf_1[2]*hamil[6]*dv11; 
  outb3[19] += 3.872983346207417*pt_conf_1[0]*hamil[8]*dv11; 
  outb3[23] += 3.872983346207417*pt_conf_1[3]*hamil[7]*dv11; 
  outb3[24] += 1.7320508075688774*pt_conf_1[3]*hamil[6]*dv11; 
  outb3[25] += 3.872983346207417*pt_conf_1[1]*hamil[8]*dv11; 
  outb3[26] += 3.872983346207417*pt_conf_1[2]*hamil[8]*dv11; 
  outb3[31] += 3.872983346207417*pt_conf_1[3]*hamil[8]*dv11; 
  } 
  if (dir == 3) { 
  double *outb0 = &out[0]; 
  } 

  const double *pt_conf_2 = &poisson_tensor_conf[8]; 
  if (dir == 1) { 
  double *outb2 = &out[72]; 
  outb2[0] += 1.7320508075688772*pt_conf_2[0]*hamil[1]*dv10; 
  outb2[1] += 1.7320508075688772*hamil[1]*pt_conf_2[1]*dv10; 
  outb2[2] += 1.7320508075688772*hamil[1]*pt_conf_2[2]*dv10; 
  outb2[3] += 3.872983346207417*pt_conf_2[0]*hamil[4]*dv10; 
  outb2[4] += 1.7320508075688772*pt_conf_2[0]*hamil[3]*dv10; 
  outb2[5] += 1.7320508075688772*hamil[1]*pt_conf_2[3]*dv10; 
  outb2[6] += 3.872983346207417*pt_conf_2[1]*hamil[4]*dv10; 
  outb2[7] += 3.872983346207417*pt_conf_2[2]*hamil[4]*dv10; 
  outb2[8] += 1.7320508075688772*pt_conf_2[1]*hamil[3]*dv10; 
  outb2[9] += 1.7320508075688772*pt_conf_2[2]*hamil[3]*dv10; 
  outb2[10] += 3.872983346207417*pt_conf_2[0]*hamil[6]*dv10; 
  outb2[12] += 1.7320508075688774*pt_conf_2[0]*hamil[7]*dv10; 
  outb2[13] += 3.872983346207417*pt_conf_2[3]*hamil[4]*dv10; 
  outb2[14] += 1.7320508075688772*hamil[3]*pt_conf_2[3]*dv10; 
  outb2[15] += 3.872983346207417*pt_conf_2[1]*hamil[6]*dv10; 
  outb2[16] += 3.872983346207417*pt_conf_2[2]*hamil[6]*dv10; 
  outb2[20] += 1.7320508075688772*pt_conf_2[1]*hamil[7]*dv10; 
  outb2[21] += 1.7320508075688772*pt_conf_2[2]*hamil[7]*dv10; 
  outb2[22] += 3.872983346207417*pt_conf_2[0]*hamil[8]*dv10; 
  outb2[23] += 3.872983346207417*pt_conf_2[3]*hamil[6]*dv10; 
  outb2[27] += 1.7320508075688774*pt_conf_2[3]*hamil[7]*dv10; 
  outb2[28] += 3.872983346207417*pt_conf_2[1]*hamil[8]*dv10; 
  outb2[29] += 3.872983346207417*pt_conf_2[2]*hamil[8]*dv10; 
  outb2[32] += 3.872983346207417*pt_conf_2[3]*hamil[8]*dv10; 
  } 
  if (dir == 2) { 
  double *outb1 = &out[36]; 
  } 

  const double *pt_conf_3 = &poisson_tensor_conf[12]; 
  if (dir == 1) { 
  double *outb3 = &out[108]; 
  outb3[0] += 1.7320508075688772*pt_conf_3[0]*hamil[2]*dv11; 
  outb3[1] += 1.7320508075688772*pt_conf_3[1]*hamil[2]*dv11; 
  outb3[2] += 1.7320508075688772*hamil[2]*pt_conf_3[2]*dv11; 
  outb3[3] += 1.7320508075688772*pt_conf_3[0]*hamil[3]*dv11; 
  outb3[4] += 3.872983346207417*pt_conf_3[0]*hamil[5]*dv11; 
  outb3[5] += 1.7320508075688772*hamil[2]*pt_conf_3[3]*dv11; 
  outb3[6] += 1.7320508075688772*pt_conf_3[1]*hamil[3]*dv11; 
  outb3[7] += 1.7320508075688772*pt_conf_3[2]*hamil[3]*dv11; 
  outb3[8] += 3.872983346207417*pt_conf_3[1]*hamil[5]*dv11; 
  outb3[9] += 3.872983346207417*pt_conf_3[2]*hamil[5]*dv11; 
  outb3[10] += 3.872983346207417*pt_conf_3[0]*hamil[7]*dv11; 
  outb3[11] += 1.7320508075688774*pt_conf_3[0]*hamil[6]*dv11; 
  outb3[13] += 1.7320508075688772*hamil[3]*pt_conf_3[3]*dv11; 
  outb3[14] += 3.872983346207417*pt_conf_3[3]*hamil[5]*dv11; 
  outb3[15] += 3.872983346207417*pt_conf_3[1]*hamil[7]*dv11; 
  outb3[16] += 3.872983346207417*pt_conf_3[2]*hamil[7]*dv11; 
  outb3[17] += 1.7320508075688772*pt_conf_3[1]*hamil[6]*dv11; 
  outb3[18] += 1.7320508075688772*pt_conf_3[2]*hamil[6]*dv11; 
  outb3[19] += 3.872983346207417*pt_conf_3[0]*hamil[8]*dv11; 
  outb3[23] += 3.872983346207417*pt_conf_3[3]*hamil[7]*dv11; 
  outb3[24] += 1.7320508075688774*pt_conf_3[3]*hamil[6]*dv11; 
  outb3[25] += 3.872983346207417*pt_conf_3[1]*hamil[8]*dv11; 
  outb3[26] += 3.872983346207417*pt_conf_3[2]*hamil[8]*dv11; 
  outb3[31] += 3.872983346207417*pt_conf_3[3]*hamil[8]*dv11; 
  } 
  if (dir == 3) { 
  double *outb1 = &out[36]; 
  } 

  const double *pt_conf_4 = &poisson_tensor_conf[16]; 
  if (dir == 2) { 
  double *outb3 = &out[108]; 
  outb3[0] += (1.224744871391589*pt_conf_4[0]*vmap_vx[2]*hamil[6]+1.224744871391589*pt_conf_4[0]*vmap_vx[1]*hamil[3]+1.224744871391589*pt_conf_4[0]*vmap_vx[0]*hamil[2])*dv11; 
  outb3[1] += (1.224744871391589*pt_conf_4[1]*vmap_vx[2]*hamil[6]+1.224744871391589*pt_conf_4[1]*vmap_vx[1]*hamil[3]+1.224744871391589*vmap_vx[0]*pt_conf_4[1]*hamil[2])*dv11; 
  outb3[2] += (1.224744871391589*pt_conf_4[2]*vmap_vx[2]*hamil[6]+1.224744871391589*vmap_vx[1]*pt_conf_4[2]*hamil[3]+1.224744871391589*vmap_vx[0]*hamil[2]*pt_conf_4[2])*dv11; 
  outb3[3] += (1.0954451150103324*pt_conf_4[0]*vmap_vx[1]*hamil[6]+1.0954451150103324*pt_conf_4[0]*vmap_vx[2]*hamil[3]+1.224744871391589*pt_conf_4[0]*vmap_vx[0]*hamil[3]+1.224744871391589*pt_conf_4[0]*vmap_vx[1]*hamil[2])*dv11; 
  outb3[4] += (2.7386127875258306*pt_conf_4[0]*vmap_vx[2]*hamil[8]+2.7386127875258306*pt_conf_4[0]*vmap_vx[1]*hamil[7]+2.7386127875258306*pt_conf_4[0]*vmap_vx[0]*hamil[5])*dv11; 
  outb3[5] += (1.224744871391589*vmap_vx[2]*pt_conf_4[3]*hamil[6]+1.224744871391589*vmap_vx[1]*hamil[3]*pt_conf_4[3]+1.224744871391589*vmap_vx[0]*hamil[2]*pt_conf_4[3])*dv11; 
  outb3[6] += (1.0954451150103324*pt_conf_4[1]*vmap_vx[1]*hamil[6]+1.0954451150103324*pt_conf_4[1]*vmap_vx[2]*hamil[3]+1.224744871391589*vmap_vx[0]*pt_conf_4[1]*hamil[3]+1.224744871391589*pt_conf_4[1]*vmap_vx[1]*hamil[2])*dv11; 
  outb3[7] += (1.0954451150103324*vmap_vx[1]*pt_conf_4[2]*hamil[6]+1.0954451150103324*pt_conf_4[2]*vmap_vx[2]*hamil[3]+1.224744871391589*vmap_vx[0]*pt_conf_4[2]*hamil[3]+1.224744871391589*vmap_vx[1]*hamil[2]*pt_conf_4[2])*dv11; 
  outb3[8] += (2.7386127875258306*pt_conf_4[1]*vmap_vx[2]*hamil[8]+2.7386127875258306*pt_conf_4[1]*vmap_vx[1]*hamil[7]+2.7386127875258306*vmap_vx[0]*pt_conf_4[1]*hamil[5])*dv11; 
  outb3[9] += (2.7386127875258306*pt_conf_4[2]*vmap_vx[2]*hamil[8]+2.7386127875258306*vmap_vx[1]*pt_conf_4[2]*hamil[7]+2.7386127875258306*vmap_vx[0]*pt_conf_4[2]*hamil[5])*dv11; 
  outb3[10] += (2.4494897427831783*pt_conf_4[0]*vmap_vx[1]*hamil[8]+2.4494897427831783*pt_conf_4[0]*vmap_vx[2]*hamil[7]+2.7386127875258306*pt_conf_4[0]*vmap_vx[0]*hamil[7]+2.7386127875258306*pt_conf_4[0]*vmap_vx[1]*hamil[5])*dv11; 
  outb3[11] += (0.7824607964359517*pt_conf_4[0]*vmap_vx[2]*hamil[6]+1.224744871391589*pt_conf_4[0]*vmap_vx[0]*hamil[6]+1.0954451150103324*pt_conf_4[0]*vmap_vx[1]*hamil[3]+1.224744871391589*pt_conf_4[0]*hamil[2]*vmap_vx[2])*dv11; 
  outb3[13] += (1.0954451150103324*vmap_vx[1]*pt_conf_4[3]*hamil[6]+1.0954451150103324*vmap_vx[2]*hamil[3]*pt_conf_4[3]+1.224744871391589*vmap_vx[0]*hamil[3]*pt_conf_4[3]+1.224744871391589*vmap_vx[1]*hamil[2]*pt_conf_4[3])*dv11; 
  outb3[14] += (2.7386127875258306*vmap_vx[2]*pt_conf_4[3]*hamil[8]+2.7386127875258306*vmap_vx[1]*pt_conf_4[3]*hamil[7]+2.7386127875258306*vmap_vx[0]*pt_conf_4[3]*hamil[5])*dv11; 
  outb3[15] += (2.4494897427831783*pt_conf_4[1]*vmap_vx[1]*hamil[8]+2.4494897427831783*pt_conf_4[1]*vmap_vx[2]*hamil[7]+2.7386127875258306*vmap_vx[0]*pt_conf_4[1]*hamil[7]+2.7386127875258306*pt_conf_4[1]*vmap_vx[1]*hamil[5])*dv11; 
  outb3[16] += (2.4494897427831783*vmap_vx[1]*pt_conf_4[2]*hamil[8]+2.4494897427831783*pt_conf_4[2]*vmap_vx[2]*hamil[7]+2.7386127875258306*vmap_vx[0]*pt_conf_4[2]*hamil[7]+2.7386127875258306*vmap_vx[1]*pt_conf_4[2]*hamil[5])*dv11; 
  outb3[17] += (0.7824607964359517*pt_conf_4[1]*vmap_vx[2]*hamil[6]+1.224744871391589*vmap_vx[0]*pt_conf_4[1]*hamil[6]+1.0954451150103324*pt_conf_4[1]*vmap_vx[1]*hamil[3]+1.224744871391589*pt_conf_4[1]*hamil[2]*vmap_vx[2])*dv11; 
  outb3[18] += (0.7824607964359517*pt_conf_4[2]*vmap_vx[2]*hamil[6]+1.224744871391589*vmap_vx[0]*pt_conf_4[2]*hamil[6]+1.0954451150103324*vmap_vx[1]*pt_conf_4[2]*hamil[3]+1.224744871391589*hamil[2]*pt_conf_4[2]*vmap_vx[2])*dv11; 
  outb3[19] += (1.7496355305594131*pt_conf_4[0]*vmap_vx[2]*hamil[8]+2.7386127875258306*pt_conf_4[0]*vmap_vx[0]*hamil[8]+2.4494897427831783*pt_conf_4[0]*vmap_vx[1]*hamil[7]+2.7386127875258306*pt_conf_4[0]*vmap_vx[2]*hamil[5])*dv11; 
  outb3[23] += (2.4494897427831783*vmap_vx[1]*pt_conf_4[3]*hamil[8]+2.4494897427831783*vmap_vx[2]*pt_conf_4[3]*hamil[7]+2.7386127875258306*vmap_vx[0]*pt_conf_4[3]*hamil[7]+2.7386127875258306*vmap_vx[1]*pt_conf_4[3]*hamil[5])*dv11; 
  outb3[24] += (0.7824607964359517*vmap_vx[2]*pt_conf_4[3]*hamil[6]+1.224744871391589*vmap_vx[0]*pt_conf_4[3]*hamil[6]+1.0954451150103324*vmap_vx[1]*hamil[3]*pt_conf_4[3]+1.224744871391589*hamil[2]*vmap_vx[2]*pt_conf_4[3])*dv11; 
  outb3[25] += (1.749635530559413*pt_conf_4[1]*vmap_vx[2]*hamil[8]+2.7386127875258306*vmap_vx[0]*pt_conf_4[1]*hamil[8]+2.4494897427831783*pt_conf_4[1]*vmap_vx[1]*hamil[7]+2.7386127875258306*pt_conf_4[1]*vmap_vx[2]*hamil[5])*dv11; 
  outb3[26] += (1.749635530559413*pt_conf_4[2]*vmap_vx[2]*hamil[8]+2.7386127875258306*vmap_vx[0]*pt_conf_4[2]*hamil[8]+2.4494897427831783*vmap_vx[1]*pt_conf_4[2]*hamil[7]+2.7386127875258306*pt_conf_4[2]*vmap_vx[2]*hamil[5])*dv11; 
  outb3[31] += (1.7496355305594131*vmap_vx[2]*pt_conf_4[3]*hamil[8]+2.7386127875258306*vmap_vx[0]*pt_conf_4[3]*hamil[8]+2.4494897427831783*vmap_vx[1]*pt_conf_4[3]*hamil[7]+2.7386127875258306*vmap_vx[2]*pt_conf_4[3]*hamil[5])*dv11; 
  } 
  if (dir == 3) { 
  double *outb2 = &out[72]; 
  outb2[0] += (-(2.7386127875258306*pt_conf_4[0]*vmap_vx[1]*hamil[4])-1.224744871391589*pt_conf_4[0]*vmap_vx[0]*hamil[1])*dv10; 
  outb2[1] += (-(2.7386127875258306*pt_conf_4[1]*vmap_vx[1]*hamil[4])-1.224744871391589*vmap_vx[0]*hamil[1]*pt_conf_4[1])*dv10; 
  outb2[2] += (-(2.7386127875258306*vmap_vx[1]*pt_conf_4[2]*hamil[4])-1.224744871391589*vmap_vx[0]*hamil[1]*pt_conf_4[2])*dv10; 
  outb2[3] += (-(2.4494897427831783*pt_conf_4[0]*vmap_vx[2]*hamil[4])-2.7386127875258306*pt_conf_4[0]*vmap_vx[0]*hamil[4]-1.224744871391589*pt_conf_4[0]*hamil[1]*vmap_vx[1])*dv10; 
  outb2[4] += (-(2.7386127875258306*pt_conf_4[0]*vmap_vx[1]*hamil[6])-1.224744871391589*pt_conf_4[0]*vmap_vx[0]*hamil[3])*dv10; 
  outb2[5] += (-(2.7386127875258306*vmap_vx[1]*pt_conf_4[3]*hamil[4])-1.224744871391589*vmap_vx[0]*hamil[1]*pt_conf_4[3])*dv10; 
  outb2[6] += (-(2.4494897427831783*pt_conf_4[1]*vmap_vx[2]*hamil[4])-2.7386127875258306*vmap_vx[0]*pt_conf_4[1]*hamil[4]-1.224744871391589*hamil[1]*pt_conf_4[1]*vmap_vx[1])*dv10; 
  outb2[7] += (-(2.4494897427831783*pt_conf_4[2]*vmap_vx[2]*hamil[4])-2.7386127875258306*vmap_vx[0]*pt_conf_4[2]*hamil[4]-1.224744871391589*hamil[1]*vmap_vx[1]*pt_conf_4[2])*dv10; 
  outb2[8] += (-(2.7386127875258306*pt_conf_4[1]*vmap_vx[1]*hamil[6])-1.224744871391589*vmap_vx[0]*pt_conf_4[1]*hamil[3])*dv10; 
  outb2[9] += (-(2.7386127875258306*vmap_vx[1]*pt_conf_4[2]*hamil[6])-1.224744871391589*vmap_vx[0]*pt_conf_4[2]*hamil[3])*dv10; 
  outb2[10] += (-(2.4494897427831783*pt_conf_4[0]*vmap_vx[2]*hamil[6])-2.7386127875258306*pt_conf_4[0]*vmap_vx[0]*hamil[6]-1.224744871391589*pt_conf_4[0]*vmap_vx[1]*hamil[3])*dv10; 
  outb2[11] += (-(2.4494897427831783*pt_conf_4[0]*vmap_vx[1]*hamil[4])-1.224744871391589*pt_conf_4[0]*hamil[1]*vmap_vx[2])*dv10; 
  outb2[12] += (-(2.7386127875258306*pt_conf_4[0]*vmap_vx[1]*hamil[8])-1.224744871391589*pt_conf_4[0]*vmap_vx[0]*hamil[7])*dv10; 
  outb2[13] += (-(2.4494897427831783*vmap_vx[2]*pt_conf_4[3]*hamil[4])-2.7386127875258306*vmap_vx[0]*pt_conf_4[3]*hamil[4]-1.224744871391589*hamil[1]*vmap_vx[1]*pt_conf_4[3])*dv10; 
  outb2[14] += (-(2.7386127875258306*vmap_vx[1]*pt_conf_4[3]*hamil[6])-1.224744871391589*vmap_vx[0]*hamil[3]*pt_conf_4[3])*dv10; 
  outb2[15] += (-(2.4494897427831783*pt_conf_4[1]*vmap_vx[2]*hamil[6])-2.7386127875258306*vmap_vx[0]*pt_conf_4[1]*hamil[6]-1.224744871391589*pt_conf_4[1]*vmap_vx[1]*hamil[3])*dv10; 
  outb2[16] += (-(2.4494897427831783*pt_conf_4[2]*vmap_vx[2]*hamil[6])-2.7386127875258306*vmap_vx[0]*pt_conf_4[2]*hamil[6]-1.224744871391589*vmap_vx[1]*pt_conf_4[2]*hamil[3])*dv10; 
  outb2[17] += (-(2.4494897427831783*pt_conf_4[1]*vmap_vx[1]*hamil[4])-1.224744871391589*hamil[1]*pt_conf_4[1]*vmap_vx[2])*dv10; 
  outb2[18] += (-(2.4494897427831783*vmap_vx[1]*pt_conf_4[2]*hamil[4])-1.224744871391589*hamil[1]*pt_conf_4[2]*vmap_vx[2])*dv10; 
  outb2[19] += (-(2.4494897427831783*pt_conf_4[0]*vmap_vx[1]*hamil[6])-1.224744871391589*pt_conf_4[0]*vmap_vx[2]*hamil[3])*dv10; 
  outb2[20] += (-(2.7386127875258306*pt_conf_4[1]*vmap_vx[1]*hamil[8])-1.224744871391589*vmap_vx[0]*pt_conf_4[1]*hamil[7])*dv10; 
  outb2[21] += (-(2.7386127875258306*vmap_vx[1]*pt_conf_4[2]*hamil[8])-1.224744871391589*vmap_vx[0]*pt_conf_4[2]*hamil[7])*dv10; 
  outb2[22] += (-(2.4494897427831783*pt_conf_4[0]*vmap_vx[2]*hamil[8])-2.7386127875258306*pt_conf_4[0]*vmap_vx[0]*hamil[8]-1.224744871391589*pt_conf_4[0]*vmap_vx[1]*hamil[7])*dv10; 
  outb2[23] += (-(2.4494897427831783*vmap_vx[2]*pt_conf_4[3]*hamil[6])-2.7386127875258306*vmap_vx[0]*pt_conf_4[3]*hamil[6]-1.224744871391589*vmap_vx[1]*hamil[3]*pt_conf_4[3])*dv10; 
  outb2[24] += (-(2.4494897427831783*vmap_vx[1]*pt_conf_4[3]*hamil[4])-1.224744871391589*hamil[1]*vmap_vx[2]*pt_conf_4[3])*dv10; 
  outb2[25] += (-(2.4494897427831783*pt_conf_4[1]*vmap_vx[1]*hamil[6])-1.224744871391589*pt_conf_4[1]*vmap_vx[2]*hamil[3])*dv10; 
  outb2[26] += (-(2.4494897427831783*vmap_vx[1]*pt_conf_4[2]*hamil[6])-1.224744871391589*pt_conf_4[2]*vmap_vx[2]*hamil[3])*dv10; 
  outb2[27] += (-(2.7386127875258306*vmap_vx[1]*pt_conf_4[3]*hamil[8])-1.224744871391589*vmap_vx[0]*pt_conf_4[3]*hamil[7])*dv10; 
  outb2[28] += (-(2.4494897427831783*pt_conf_4[1]*vmap_vx[2]*hamil[8])-2.7386127875258306*vmap_vx[0]*pt_conf_4[1]*hamil[8]-1.224744871391589*pt_conf_4[1]*vmap_vx[1]*hamil[7])*dv10; 
  outb2[29] += (-(2.4494897427831783*pt_conf_4[2]*vmap_vx[2]*hamil[8])-2.7386127875258306*vmap_vx[0]*pt_conf_4[2]*hamil[8]-1.224744871391589*vmap_vx[1]*pt_conf_4[2]*hamil[7])*dv10; 
  outb2[30] += (-(2.4494897427831783*pt_conf_4[0]*vmap_vx[1]*hamil[8])-1.224744871391589*pt_conf_4[0]*vmap_vx[2]*hamil[7])*dv10; 
  outb2[31] += (-(2.4494897427831783*vmap_vx[1]*pt_conf_4[3]*hamil[6])-1.224744871391589*vmap_vx[2]*hamil[3]*pt_conf_4[3])*dv10; 
  outb2[32] += (-(2.4494897427831783*vmap_vx[2]*pt_conf_4[3]*hamil[8])-2.7386127875258306*vmap_vx[0]*pt_conf_4[3]*hamil[8]-1.224744871391589*vmap_vx[1]*pt_conf_4[3]*hamil[7])*dv10; 
  outb2[33] += (-(2.4494897427831783*pt_conf_4[1]*vmap_vx[1]*hamil[8])-1.224744871391589*pt_conf_4[1]*vmap_vx[2]*hamil[7])*dv10; 
  outb2[34] += (-(2.4494897427831783*vmap_vx[1]*pt_conf_4[2]*hamil[8])-1.224744871391589*pt_conf_4[2]*vmap_vx[2]*hamil[7])*dv10; 
  outb2[35] += (-(2.4494897427831783*vmap_vx[1]*pt_conf_4[3]*hamil[8])-1.224744871391589*vmap_vx[2]*pt_conf_4[3]*hamil[7])*dv10; 
  } 

  const double *pt_conf_5 = &poisson_tensor_conf[20]; 
  if (dir == 2) { 
  double *outb3 = &out[108]; 
  outb3[0] += (2.7386127875258306*pt_conf_5[0]*vmap_vy[1]*hamil[5]+1.224744871391589*pt_conf_5[0]*vmap_vy[0]*hamil[2])*dv11; 
  outb3[1] += (2.7386127875258306*pt_conf_5[1]*vmap_vy[1]*hamil[5]+1.224744871391589*vmap_vy[0]*pt_conf_5[1]*hamil[2])*dv11; 
  outb3[2] += (2.7386127875258306*vmap_vy[1]*pt_conf_5[2]*hamil[5]+1.224744871391589*vmap_vy[0]*hamil[2]*pt_conf_5[2])*dv11; 
  outb3[3] += (2.7386127875258306*pt_conf_5[0]*vmap_vy[1]*hamil[7]+1.224744871391589*pt_conf_5[0]*vmap_vy[0]*hamil[3])*dv11; 
  outb3[4] += (2.4494897427831783*pt_conf_5[0]*vmap_vy[2]*hamil[5]+2.7386127875258306*pt_conf_5[0]*vmap_vy[0]*hamil[5]+1.224744871391589*pt_conf_5[0]*vmap_vy[1]*hamil[2])*dv11; 
  outb3[5] += (2.7386127875258306*vmap_vy[1]*pt_conf_5[3]*hamil[5]+1.224744871391589*vmap_vy[0]*hamil[2]*pt_conf_5[3])*dv11; 
  outb3[6] += (2.7386127875258306*pt_conf_5[1]*vmap_vy[1]*hamil[7]+1.224744871391589*vmap_vy[0]*pt_conf_5[1]*hamil[3])*dv11; 
  outb3[7] += (2.7386127875258306*vmap_vy[1]*pt_conf_5[2]*hamil[7]+1.224744871391589*vmap_vy[0]*pt_conf_5[2]*hamil[3])*dv11; 
  outb3[8] += (2.4494897427831783*pt_conf_5[1]*vmap_vy[2]*hamil[5]+2.7386127875258306*vmap_vy[0]*pt_conf_5[1]*hamil[5]+1.224744871391589*pt_conf_5[1]*vmap_vy[1]*hamil[2])*dv11; 
  outb3[9] += (2.4494897427831783*pt_conf_5[2]*vmap_vy[2]*hamil[5]+2.7386127875258306*vmap_vy[0]*pt_conf_5[2]*hamil[5]+1.224744871391589*vmap_vy[1]*hamil[2]*pt_conf_5[2])*dv11; 
  outb3[10] += (2.4494897427831783*pt_conf_5[0]*vmap_vy[2]*hamil[7]+2.7386127875258306*pt_conf_5[0]*vmap_vy[0]*hamil[7]+1.224744871391589*pt_conf_5[0]*vmap_vy[1]*hamil[3])*dv11; 
  outb3[11] += (2.7386127875258306*pt_conf_5[0]*vmap_vy[1]*hamil[8]+1.224744871391589*pt_conf_5[0]*vmap_vy[0]*hamil[6])*dv11; 
  outb3[12] += (2.4494897427831783*pt_conf_5[0]*vmap_vy[1]*hamil[5]+1.224744871391589*pt_conf_5[0]*hamil[2]*vmap_vy[2])*dv11; 
  outb3[13] += (2.7386127875258306*vmap_vy[1]*pt_conf_5[3]*hamil[7]+1.224744871391589*vmap_vy[0]*hamil[3]*pt_conf_5[3])*dv11; 
  outb3[14] += (2.4494897427831783*vmap_vy[2]*pt_conf_5[3]*hamil[5]+2.7386127875258306*vmap_vy[0]*pt_conf_5[3]*hamil[5]+1.224744871391589*vmap_vy[1]*hamil[2]*pt_conf_5[3])*dv11; 
  outb3[15] += (2.4494897427831783*pt_conf_5[1]*vmap_vy[2]*hamil[7]+2.7386127875258306*vmap_vy[0]*pt_conf_5[1]*hamil[7]+1.224744871391589*pt_conf_5[1]*vmap_vy[1]*hamil[3])*dv11; 
  outb3[16] += (2.4494897427831783*pt_conf_5[2]*vmap_vy[2]*hamil[7]+2.7386127875258306*vmap_vy[0]*pt_conf_5[2]*hamil[7]+1.224744871391589*vmap_vy[1]*pt_conf_5[2]*hamil[3])*dv11; 
  outb3[17] += (2.7386127875258306*pt_conf_5[1]*vmap_vy[1]*hamil[8]+1.224744871391589*vmap_vy[0]*pt_conf_5[1]*hamil[6])*dv11; 
  outb3[18] += (2.7386127875258306*vmap_vy[1]*pt_conf_5[2]*hamil[8]+1.224744871391589*vmap_vy[0]*pt_conf_5[2]*hamil[6])*dv11; 
  outb3[19] += (2.4494897427831783*pt_conf_5[0]*vmap_vy[2]*hamil[8]+2.7386127875258306*pt_conf_5[0]*vmap_vy[0]*hamil[8]+1.224744871391589*pt_conf_5[0]*vmap_vy[1]*hamil[6])*dv11; 
  outb3[20] += (2.4494897427831783*pt_conf_5[1]*vmap_vy[1]*hamil[5]+1.224744871391589*pt_conf_5[1]*hamil[2]*vmap_vy[2])*dv11; 
  outb3[21] += (2.4494897427831783*vmap_vy[1]*pt_conf_5[2]*hamil[5]+1.224744871391589*hamil[2]*pt_conf_5[2]*vmap_vy[2])*dv11; 
  outb3[22] += (2.4494897427831783*pt_conf_5[0]*vmap_vy[1]*hamil[7]+1.224744871391589*pt_conf_5[0]*vmap_vy[2]*hamil[3])*dv11; 
  outb3[23] += (2.4494897427831783*vmap_vy[2]*pt_conf_5[3]*hamil[7]+2.7386127875258306*vmap_vy[0]*pt_conf_5[3]*hamil[7]+1.224744871391589*vmap_vy[1]*hamil[3]*pt_conf_5[3])*dv11; 
  outb3[24] += (2.7386127875258306*vmap_vy[1]*pt_conf_5[3]*hamil[8]+1.224744871391589*vmap_vy[0]*pt_conf_5[3]*hamil[6])*dv11; 
  outb3[25] += (2.4494897427831783*pt_conf_5[1]*vmap_vy[2]*hamil[8]+2.7386127875258306*vmap_vy[0]*pt_conf_5[1]*hamil[8]+1.224744871391589*pt_conf_5[1]*vmap_vy[1]*hamil[6])*dv11; 
  outb3[26] += (2.4494897427831783*pt_conf_5[2]*vmap_vy[2]*hamil[8]+2.7386127875258306*vmap_vy[0]*pt_conf_5[2]*hamil[8]+1.224744871391589*vmap_vy[1]*pt_conf_5[2]*hamil[6])*dv11; 
  outb3[27] += (2.4494897427831783*vmap_vy[1]*pt_conf_5[3]*hamil[5]+1.224744871391589*hamil[2]*vmap_vy[2]*pt_conf_5[3])*dv11; 
  outb3[28] += (2.4494897427831783*pt_conf_5[1]*vmap_vy[1]*hamil[7]+1.224744871391589*pt_conf_5[1]*vmap_vy[2]*hamil[3])*dv11; 
  outb3[29] += (2.4494897427831783*vmap_vy[1]*pt_conf_5[2]*hamil[7]+1.224744871391589*pt_conf_5[2]*vmap_vy[2]*hamil[3])*dv11; 
  outb3[30] += (2.4494897427831783*pt_conf_5[0]*vmap_vy[1]*hamil[8]+1.224744871391589*pt_conf_5[0]*vmap_vy[2]*hamil[6])*dv11; 
  outb3[31] += (2.4494897427831783*vmap_vy[2]*pt_conf_5[3]*hamil[8]+2.7386127875258306*vmap_vy[0]*pt_conf_5[3]*hamil[8]+1.224744871391589*vmap_vy[1]*pt_conf_5[3]*hamil[6])*dv11; 
  outb3[32] += (2.4494897427831783*vmap_vy[1]*pt_conf_5[3]*hamil[7]+1.224744871391589*vmap_vy[2]*hamil[3]*pt_conf_5[3])*dv11; 
  outb3[33] += (2.4494897427831783*pt_conf_5[1]*vmap_vy[1]*hamil[8]+1.224744871391589*pt_conf_5[1]*vmap_vy[2]*hamil[6])*dv11; 
  outb3[34] += (2.4494897427831783*vmap_vy[1]*pt_conf_5[2]*hamil[8]+1.224744871391589*pt_conf_5[2]*vmap_vy[2]*hamil[6])*dv11; 
  outb3[35] += (2.4494897427831783*vmap_vy[1]*pt_conf_5[3]*hamil[8]+1.224744871391589*vmap_vy[2]*pt_conf_5[3]*hamil[6])*dv11; 
  } 
  if (dir == 3) { 
  double *outb2 = &out[72]; 
  outb2[0] += (-(1.224744871391589*pt_conf_5[0]*vmap_vy[2]*hamil[7])-1.224744871391589*pt_conf_5[0]*vmap_vy[1]*hamil[3]-1.224744871391589*pt_conf_5[0]*vmap_vy[0]*hamil[1])*dv10; 
  outb2[1] += (-(1.224744871391589*pt_conf_5[1]*vmap_vy[2]*hamil[7])-1.224744871391589*pt_conf_5[1]*vmap_vy[1]*hamil[3]-1.224744871391589*vmap_vy[0]*hamil[1]*pt_conf_5[1])*dv10; 
  outb2[2] += (-(1.224744871391589*pt_conf_5[2]*vmap_vy[2]*hamil[7])-1.224744871391589*vmap_vy[1]*pt_conf_5[2]*hamil[3]-1.224744871391589*vmap_vy[0]*hamil[1]*pt_conf_5[2])*dv10; 
  outb2[3] += (-(2.7386127875258306*pt_conf_5[0]*vmap_vy[2]*hamil[8])-2.7386127875258306*pt_conf_5[0]*vmap_vy[1]*hamil[6]-2.7386127875258306*pt_conf_5[0]*vmap_vy[0]*hamil[4])*dv10; 
  outb2[4] += (-(1.0954451150103324*pt_conf_5[0]*vmap_vy[1]*hamil[7])-1.0954451150103324*pt_conf_5[0]*vmap_vy[2]*hamil[3]-1.224744871391589*pt_conf_5[0]*vmap_vy[0]*hamil[3]-1.224744871391589*pt_conf_5[0]*hamil[1]*vmap_vy[1])*dv10; 
  outb2[5] += (-(1.224744871391589*vmap_vy[2]*pt_conf_5[3]*hamil[7])-1.224744871391589*vmap_vy[1]*hamil[3]*pt_conf_5[3]-1.224744871391589*vmap_vy[0]*hamil[1]*pt_conf_5[3])*dv10; 
  outb2[6] += (-(2.7386127875258306*pt_conf_5[1]*vmap_vy[2]*hamil[8])-2.7386127875258306*pt_conf_5[1]*vmap_vy[1]*hamil[6]-2.7386127875258306*vmap_vy[0]*pt_conf_5[1]*hamil[4])*dv10; 
  outb2[7] += (-(2.7386127875258306*pt_conf_5[2]*vmap_vy[2]*hamil[8])-2.7386127875258306*vmap_vy[1]*pt_conf_5[2]*hamil[6]-2.7386127875258306*vmap_vy[0]*pt_conf_5[2]*hamil[4])*dv10; 
  outb2[8] += (-(1.0954451150103324*pt_conf_5[1]*vmap_vy[1]*hamil[7])-1.0954451150103324*pt_conf_5[1]*vmap_vy[2]*hamil[3]-1.224744871391589*vmap_vy[0]*pt_conf_5[1]*hamil[3]-1.224744871391589*hamil[1]*pt_conf_5[1]*vmap_vy[1])*dv10; 
  outb2[9] += (-(1.0954451150103324*vmap_vy[1]*pt_conf_5[2]*hamil[7])-1.0954451150103324*pt_conf_5[2]*vmap_vy[2]*hamil[3]-1.224744871391589*vmap_vy[0]*pt_conf_5[2]*hamil[3]-1.224744871391589*hamil[1]*vmap_vy[1]*pt_conf_5[2])*dv10; 
  outb2[10] += (-(2.4494897427831783*pt_conf_5[0]*vmap_vy[1]*hamil[8])-2.4494897427831783*pt_conf_5[0]*vmap_vy[2]*hamil[6]-2.7386127875258306*pt_conf_5[0]*vmap_vy[0]*hamil[6]-2.7386127875258306*pt_conf_5[0]*vmap_vy[1]*hamil[4])*dv10; 
  outb2[12] += (-(0.7824607964359517*pt_conf_5[0]*vmap_vy[2]*hamil[7])-1.224744871391589*pt_conf_5[0]*vmap_vy[0]*hamil[7]-1.0954451150103324*pt_conf_5[0]*vmap_vy[1]*hamil[3]-1.224744871391589*pt_conf_5[0]*hamil[1]*vmap_vy[2])*dv10; 
  outb2[13] += (-(2.7386127875258306*vmap_vy[2]*pt_conf_5[3]*hamil[8])-2.7386127875258306*vmap_vy[1]*pt_conf_5[3]*hamil[6]-2.7386127875258306*vmap_vy[0]*pt_conf_5[3]*hamil[4])*dv10; 
  outb2[14] += (-(1.0954451150103324*vmap_vy[1]*pt_conf_5[3]*hamil[7])-1.0954451150103324*vmap_vy[2]*hamil[3]*pt_conf_5[3]-1.224744871391589*vmap_vy[0]*hamil[3]*pt_conf_5[3]-1.224744871391589*hamil[1]*vmap_vy[1]*pt_conf_5[3])*dv10; 
  outb2[15] += (-(2.4494897427831783*pt_conf_5[1]*vmap_vy[1]*hamil[8])-2.4494897427831783*pt_conf_5[1]*vmap_vy[2]*hamil[6]-2.7386127875258306*vmap_vy[0]*pt_conf_5[1]*hamil[6]-2.7386127875258306*pt_conf_5[1]*vmap_vy[1]*hamil[4])*dv10; 
  outb2[16] += (-(2.4494897427831783*vmap_vy[1]*pt_conf_5[2]*hamil[8])-2.4494897427831783*pt_conf_5[2]*vmap_vy[2]*hamil[6]-2.7386127875258306*vmap_vy[0]*pt_conf_5[2]*hamil[6]-2.7386127875258306*vmap_vy[1]*pt_conf_5[2]*hamil[4])*dv10; 
  outb2[20] += (-(0.7824607964359517*pt_conf_5[1]*vmap_vy[2]*hamil[7])-1.224744871391589*vmap_vy[0]*pt_conf_5[1]*hamil[7]-1.0954451150103324*pt_conf_5[1]*vmap_vy[1]*hamil[3]-1.224744871391589*hamil[1]*pt_conf_5[1]*vmap_vy[2])*dv10; 
  outb2[21] += (-(0.7824607964359517*pt_conf_5[2]*vmap_vy[2]*hamil[7])-1.224744871391589*vmap_vy[0]*pt_conf_5[2]*hamil[7]-1.0954451150103324*vmap_vy[1]*pt_conf_5[2]*hamil[3]-1.224744871391589*hamil[1]*pt_conf_5[2]*vmap_vy[2])*dv10; 
  outb2[22] += (-(1.7496355305594131*pt_conf_5[0]*vmap_vy[2]*hamil[8])-2.7386127875258306*pt_conf_5[0]*vmap_vy[0]*hamil[8]-2.4494897427831783*pt_conf_5[0]*vmap_vy[1]*hamil[6]-2.7386127875258306*pt_conf_5[0]*vmap_vy[2]*hamil[4])*dv10; 
  outb2[23] += (-(2.4494897427831783*vmap_vy[1]*pt_conf_5[3]*hamil[8])-2.4494897427831783*vmap_vy[2]*pt_conf_5[3]*hamil[6]-2.7386127875258306*vmap_vy[0]*pt_conf_5[3]*hamil[6]-2.7386127875258306*vmap_vy[1]*pt_conf_5[3]*hamil[4])*dv10; 
  outb2[27] += (-(0.7824607964359517*vmap_vy[2]*pt_conf_5[3]*hamil[7])-1.224744871391589*vmap_vy[0]*pt_conf_5[3]*hamil[7]-1.0954451150103324*vmap_vy[1]*hamil[3]*pt_conf_5[3]-1.224744871391589*hamil[1]*vmap_vy[2]*pt_conf_5[3])*dv10; 
  outb2[28] += (-(1.749635530559413*pt_conf_5[1]*vmap_vy[2]*hamil[8])-2.7386127875258306*vmap_vy[0]*pt_conf_5[1]*hamil[8]-2.4494897427831783*pt_conf_5[1]*vmap_vy[1]*hamil[6]-2.7386127875258306*pt_conf_5[1]*vmap_vy[2]*hamil[4])*dv10; 
  outb2[29] += (-(1.749635530559413*pt_conf_5[2]*vmap_vy[2]*hamil[8])-2.7386127875258306*vmap_vy[0]*pt_conf_5[2]*hamil[8]-2.4494897427831783*vmap_vy[1]*pt_conf_5[2]*hamil[6]-2.7386127875258306*pt_conf_5[2]*vmap_vy[2]*hamil[4])*dv10; 
  outb2[32] += (-(1.7496355305594131*vmap_vy[2]*pt_conf_5[3]*hamil[8])-2.7386127875258306*vmap_vy[0]*pt_conf_5[3]*hamil[8]-2.4494897427831783*vmap_vy[1]*pt_conf_5[3]*hamil[6]-2.7386127875258306*vmap_vy[2]*pt_conf_5[3]*hamil[4])*dv10; 
  } 

} 
