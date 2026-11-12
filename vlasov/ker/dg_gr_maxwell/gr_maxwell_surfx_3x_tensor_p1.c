#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH double gr_maxwell_surfx_3x_tensor_p1(const double *w, const double *dxv,
  const double *Fhat_l_nodal, const double *Fhat_r_nodal, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 

  const double *Fhat_l_nodal_1 = &Fhat_l_nodal[0]; 
  const double *Fhat_r_nodal_1 = &Fhat_r_nodal[0]; 
  double *out_1 = &out[0]; 
  out_1[0] += (-(0.3535533905932737*Fhat_r_nodal_1[3])+0.3535533905932737*Fhat_l_nodal_1[3]-0.3535533905932737*Fhat_r_nodal_1[2]+0.3535533905932737*Fhat_l_nodal_1[2]-0.3535533905932737*Fhat_r_nodal_1[1]+0.3535533905932737*Fhat_l_nodal_1[1]-0.3535533905932737*Fhat_r_nodal_1[0]+0.3535533905932737*Fhat_l_nodal_1[0])*dx10; 
  out_1[1] += -(0.6123724356957944*(Fhat_r_nodal_1[3]+Fhat_l_nodal_1[3]+Fhat_r_nodal_1[2]+Fhat_l_nodal_1[2]+Fhat_r_nodal_1[1]+Fhat_l_nodal_1[1]+Fhat_r_nodal_1[0]+Fhat_l_nodal_1[0])*dx10); 
  out_1[2] += (-(0.3535533905932737*Fhat_r_nodal_1[3])+0.3535533905932737*Fhat_l_nodal_1[3]-0.3535533905932737*Fhat_r_nodal_1[2]+0.3535533905932737*(Fhat_l_nodal_1[2]+Fhat_r_nodal_1[1])-0.3535533905932737*Fhat_l_nodal_1[1]+0.3535533905932737*Fhat_r_nodal_1[0]-0.3535533905932737*Fhat_l_nodal_1[0])*dx10; 
  out_1[3] += (-(0.3535533905932737*Fhat_r_nodal_1[3])+0.3535533905932737*(Fhat_l_nodal_1[3]+Fhat_r_nodal_1[2])-0.3535533905932737*(Fhat_l_nodal_1[2]+Fhat_r_nodal_1[1])+0.3535533905932737*(Fhat_l_nodal_1[1]+Fhat_r_nodal_1[0])-0.3535533905932737*Fhat_l_nodal_1[0])*dx10; 
  out_1[4] += (0.6123724356957944*(Fhat_r_nodal_1[1]+Fhat_l_nodal_1[1]+Fhat_r_nodal_1[0]+Fhat_l_nodal_1[0])-0.6123724356957944*(Fhat_r_nodal_1[3]+Fhat_l_nodal_1[3]+Fhat_r_nodal_1[2]+Fhat_l_nodal_1[2]))*dx10; 
  out_1[5] += (-(0.6123724356957944*(Fhat_r_nodal_1[3]+Fhat_l_nodal_1[3]))+0.6123724356957944*(Fhat_r_nodal_1[2]+Fhat_l_nodal_1[2])-0.6123724356957944*(Fhat_r_nodal_1[1]+Fhat_l_nodal_1[1])+0.6123724356957944*(Fhat_r_nodal_1[0]+Fhat_l_nodal_1[0]))*dx10; 
  out_1[6] += (-(0.3535533905932737*Fhat_r_nodal_1[3])+0.3535533905932737*(Fhat_l_nodal_1[3]+Fhat_r_nodal_1[2])-0.3535533905932737*Fhat_l_nodal_1[2]+0.3535533905932737*Fhat_r_nodal_1[1]-0.3535533905932737*(Fhat_l_nodal_1[1]+Fhat_r_nodal_1[0])+0.3535533905932737*Fhat_l_nodal_1[0])*dx10; 
  out_1[7] += (-(0.6123724356957944*(Fhat_r_nodal_1[3]+Fhat_l_nodal_1[3]))+0.6123724356957944*(Fhat_r_nodal_1[2]+Fhat_l_nodal_1[2]+Fhat_r_nodal_1[1]+Fhat_l_nodal_1[1])-0.6123724356957944*(Fhat_r_nodal_1[0]+Fhat_l_nodal_1[0]))*dx10; 

  const double *Fhat_l_nodal_2 = &Fhat_l_nodal[4]; 
  const double *Fhat_r_nodal_2 = &Fhat_r_nodal[4]; 
  double *out_2 = &out[8]; 
  out_2[0] += (-(0.3535533905932737*Fhat_r_nodal_2[3])+0.3535533905932737*Fhat_l_nodal_2[3]-0.3535533905932737*Fhat_r_nodal_2[2]+0.3535533905932737*Fhat_l_nodal_2[2]-0.3535533905932737*Fhat_r_nodal_2[1]+0.3535533905932737*Fhat_l_nodal_2[1]-0.3535533905932737*Fhat_r_nodal_2[0]+0.3535533905932737*Fhat_l_nodal_2[0])*dx10; 
  out_2[1] += -(0.6123724356957944*(Fhat_r_nodal_2[3]+Fhat_l_nodal_2[3]+Fhat_r_nodal_2[2]+Fhat_l_nodal_2[2]+Fhat_r_nodal_2[1]+Fhat_l_nodal_2[1]+Fhat_r_nodal_2[0]+Fhat_l_nodal_2[0])*dx10); 
  out_2[2] += (-(0.3535533905932737*Fhat_r_nodal_2[3])+0.3535533905932737*Fhat_l_nodal_2[3]-0.3535533905932737*Fhat_r_nodal_2[2]+0.3535533905932737*(Fhat_l_nodal_2[2]+Fhat_r_nodal_2[1])-0.3535533905932737*Fhat_l_nodal_2[1]+0.3535533905932737*Fhat_r_nodal_2[0]-0.3535533905932737*Fhat_l_nodal_2[0])*dx10; 
  out_2[3] += (-(0.3535533905932737*Fhat_r_nodal_2[3])+0.3535533905932737*(Fhat_l_nodal_2[3]+Fhat_r_nodal_2[2])-0.3535533905932737*(Fhat_l_nodal_2[2]+Fhat_r_nodal_2[1])+0.3535533905932737*(Fhat_l_nodal_2[1]+Fhat_r_nodal_2[0])-0.3535533905932737*Fhat_l_nodal_2[0])*dx10; 
  out_2[4] += (0.6123724356957944*(Fhat_r_nodal_2[1]+Fhat_l_nodal_2[1]+Fhat_r_nodal_2[0]+Fhat_l_nodal_2[0])-0.6123724356957944*(Fhat_r_nodal_2[3]+Fhat_l_nodal_2[3]+Fhat_r_nodal_2[2]+Fhat_l_nodal_2[2]))*dx10; 
  out_2[5] += (-(0.6123724356957944*(Fhat_r_nodal_2[3]+Fhat_l_nodal_2[3]))+0.6123724356957944*(Fhat_r_nodal_2[2]+Fhat_l_nodal_2[2])-0.6123724356957944*(Fhat_r_nodal_2[1]+Fhat_l_nodal_2[1])+0.6123724356957944*(Fhat_r_nodal_2[0]+Fhat_l_nodal_2[0]))*dx10; 
  out_2[6] += (-(0.3535533905932737*Fhat_r_nodal_2[3])+0.3535533905932737*(Fhat_l_nodal_2[3]+Fhat_r_nodal_2[2])-0.3535533905932737*Fhat_l_nodal_2[2]+0.3535533905932737*Fhat_r_nodal_2[1]-0.3535533905932737*(Fhat_l_nodal_2[1]+Fhat_r_nodal_2[0])+0.3535533905932737*Fhat_l_nodal_2[0])*dx10; 
  out_2[7] += (-(0.6123724356957944*(Fhat_r_nodal_2[3]+Fhat_l_nodal_2[3]))+0.6123724356957944*(Fhat_r_nodal_2[2]+Fhat_l_nodal_2[2]+Fhat_r_nodal_2[1]+Fhat_l_nodal_2[1])-0.6123724356957944*(Fhat_r_nodal_2[0]+Fhat_l_nodal_2[0]))*dx10; 

  const double *Fhat_l_nodal_3 = &Fhat_l_nodal[8]; 
  const double *Fhat_r_nodal_3 = &Fhat_r_nodal[8]; 
  double *out_3 = &out[16]; 
  out_3[0] += (-(0.3535533905932737*Fhat_r_nodal_3[3])+0.3535533905932737*Fhat_l_nodal_3[3]-0.3535533905932737*Fhat_r_nodal_3[2]+0.3535533905932737*Fhat_l_nodal_3[2]-0.3535533905932737*Fhat_r_nodal_3[1]+0.3535533905932737*Fhat_l_nodal_3[1]-0.3535533905932737*Fhat_r_nodal_3[0]+0.3535533905932737*Fhat_l_nodal_3[0])*dx10; 
  out_3[1] += -(0.6123724356957944*(Fhat_r_nodal_3[3]+Fhat_l_nodal_3[3]+Fhat_r_nodal_3[2]+Fhat_l_nodal_3[2]+Fhat_r_nodal_3[1]+Fhat_l_nodal_3[1]+Fhat_r_nodal_3[0]+Fhat_l_nodal_3[0])*dx10); 
  out_3[2] += (-(0.3535533905932737*Fhat_r_nodal_3[3])+0.3535533905932737*Fhat_l_nodal_3[3]-0.3535533905932737*Fhat_r_nodal_3[2]+0.3535533905932737*(Fhat_l_nodal_3[2]+Fhat_r_nodal_3[1])-0.3535533905932737*Fhat_l_nodal_3[1]+0.3535533905932737*Fhat_r_nodal_3[0]-0.3535533905932737*Fhat_l_nodal_3[0])*dx10; 
  out_3[3] += (-(0.3535533905932737*Fhat_r_nodal_3[3])+0.3535533905932737*(Fhat_l_nodal_3[3]+Fhat_r_nodal_3[2])-0.3535533905932737*(Fhat_l_nodal_3[2]+Fhat_r_nodal_3[1])+0.3535533905932737*(Fhat_l_nodal_3[1]+Fhat_r_nodal_3[0])-0.3535533905932737*Fhat_l_nodal_3[0])*dx10; 
  out_3[4] += (0.6123724356957944*(Fhat_r_nodal_3[1]+Fhat_l_nodal_3[1]+Fhat_r_nodal_3[0]+Fhat_l_nodal_3[0])-0.6123724356957944*(Fhat_r_nodal_3[3]+Fhat_l_nodal_3[3]+Fhat_r_nodal_3[2]+Fhat_l_nodal_3[2]))*dx10; 
  out_3[5] += (-(0.6123724356957944*(Fhat_r_nodal_3[3]+Fhat_l_nodal_3[3]))+0.6123724356957944*(Fhat_r_nodal_3[2]+Fhat_l_nodal_3[2])-0.6123724356957944*(Fhat_r_nodal_3[1]+Fhat_l_nodal_3[1])+0.6123724356957944*(Fhat_r_nodal_3[0]+Fhat_l_nodal_3[0]))*dx10; 
  out_3[6] += (-(0.3535533905932737*Fhat_r_nodal_3[3])+0.3535533905932737*(Fhat_l_nodal_3[3]+Fhat_r_nodal_3[2])-0.3535533905932737*Fhat_l_nodal_3[2]+0.3535533905932737*Fhat_r_nodal_3[1]-0.3535533905932737*(Fhat_l_nodal_3[1]+Fhat_r_nodal_3[0])+0.3535533905932737*Fhat_l_nodal_3[0])*dx10; 
  out_3[7] += (-(0.6123724356957944*(Fhat_r_nodal_3[3]+Fhat_l_nodal_3[3]))+0.6123724356957944*(Fhat_r_nodal_3[2]+Fhat_l_nodal_3[2]+Fhat_r_nodal_3[1]+Fhat_l_nodal_3[1])-0.6123724356957944*(Fhat_r_nodal_3[0]+Fhat_l_nodal_3[0]))*dx10; 

  const double *Fhat_l_nodal_4 = &Fhat_l_nodal[12]; 
  const double *Fhat_r_nodal_4 = &Fhat_r_nodal[12]; 
  double *out_4 = &out[24]; 
  out_4[0] += (-(0.3535533905932737*Fhat_r_nodal_4[3])+0.3535533905932737*Fhat_l_nodal_4[3]-0.3535533905932737*Fhat_r_nodal_4[2]+0.3535533905932737*Fhat_l_nodal_4[2]-0.3535533905932737*Fhat_r_nodal_4[1]+0.3535533905932737*Fhat_l_nodal_4[1]-0.3535533905932737*Fhat_r_nodal_4[0]+0.3535533905932737*Fhat_l_nodal_4[0])*dx10; 
  out_4[1] += -(0.6123724356957944*(Fhat_r_nodal_4[3]+Fhat_l_nodal_4[3]+Fhat_r_nodal_4[2]+Fhat_l_nodal_4[2]+Fhat_r_nodal_4[1]+Fhat_l_nodal_4[1]+Fhat_r_nodal_4[0]+Fhat_l_nodal_4[0])*dx10); 
  out_4[2] += (-(0.3535533905932737*Fhat_r_nodal_4[3])+0.3535533905932737*Fhat_l_nodal_4[3]-0.3535533905932737*Fhat_r_nodal_4[2]+0.3535533905932737*(Fhat_l_nodal_4[2]+Fhat_r_nodal_4[1])-0.3535533905932737*Fhat_l_nodal_4[1]+0.3535533905932737*Fhat_r_nodal_4[0]-0.3535533905932737*Fhat_l_nodal_4[0])*dx10; 
  out_4[3] += (-(0.3535533905932737*Fhat_r_nodal_4[3])+0.3535533905932737*(Fhat_l_nodal_4[3]+Fhat_r_nodal_4[2])-0.3535533905932737*(Fhat_l_nodal_4[2]+Fhat_r_nodal_4[1])+0.3535533905932737*(Fhat_l_nodal_4[1]+Fhat_r_nodal_4[0])-0.3535533905932737*Fhat_l_nodal_4[0])*dx10; 
  out_4[4] += (0.6123724356957944*(Fhat_r_nodal_4[1]+Fhat_l_nodal_4[1]+Fhat_r_nodal_4[0]+Fhat_l_nodal_4[0])-0.6123724356957944*(Fhat_r_nodal_4[3]+Fhat_l_nodal_4[3]+Fhat_r_nodal_4[2]+Fhat_l_nodal_4[2]))*dx10; 
  out_4[5] += (-(0.6123724356957944*(Fhat_r_nodal_4[3]+Fhat_l_nodal_4[3]))+0.6123724356957944*(Fhat_r_nodal_4[2]+Fhat_l_nodal_4[2])-0.6123724356957944*(Fhat_r_nodal_4[1]+Fhat_l_nodal_4[1])+0.6123724356957944*(Fhat_r_nodal_4[0]+Fhat_l_nodal_4[0]))*dx10; 
  out_4[6] += (-(0.3535533905932737*Fhat_r_nodal_4[3])+0.3535533905932737*(Fhat_l_nodal_4[3]+Fhat_r_nodal_4[2])-0.3535533905932737*Fhat_l_nodal_4[2]+0.3535533905932737*Fhat_r_nodal_4[1]-0.3535533905932737*(Fhat_l_nodal_4[1]+Fhat_r_nodal_4[0])+0.3535533905932737*Fhat_l_nodal_4[0])*dx10; 
  out_4[7] += (-(0.6123724356957944*(Fhat_r_nodal_4[3]+Fhat_l_nodal_4[3]))+0.6123724356957944*(Fhat_r_nodal_4[2]+Fhat_l_nodal_4[2]+Fhat_r_nodal_4[1]+Fhat_l_nodal_4[1])-0.6123724356957944*(Fhat_r_nodal_4[0]+Fhat_l_nodal_4[0]))*dx10; 

  const double *Fhat_l_nodal_5 = &Fhat_l_nodal[16]; 
  const double *Fhat_r_nodal_5 = &Fhat_r_nodal[16]; 
  double *out_5 = &out[32]; 
  out_5[0] += (-(0.3535533905932737*Fhat_r_nodal_5[3])+0.3535533905932737*Fhat_l_nodal_5[3]-0.3535533905932737*Fhat_r_nodal_5[2]+0.3535533905932737*Fhat_l_nodal_5[2]-0.3535533905932737*Fhat_r_nodal_5[1]+0.3535533905932737*Fhat_l_nodal_5[1]-0.3535533905932737*Fhat_r_nodal_5[0]+0.3535533905932737*Fhat_l_nodal_5[0])*dx10; 
  out_5[1] += -(0.6123724356957944*(Fhat_r_nodal_5[3]+Fhat_l_nodal_5[3]+Fhat_r_nodal_5[2]+Fhat_l_nodal_5[2]+Fhat_r_nodal_5[1]+Fhat_l_nodal_5[1]+Fhat_r_nodal_5[0]+Fhat_l_nodal_5[0])*dx10); 
  out_5[2] += (-(0.3535533905932737*Fhat_r_nodal_5[3])+0.3535533905932737*Fhat_l_nodal_5[3]-0.3535533905932737*Fhat_r_nodal_5[2]+0.3535533905932737*(Fhat_l_nodal_5[2]+Fhat_r_nodal_5[1])-0.3535533905932737*Fhat_l_nodal_5[1]+0.3535533905932737*Fhat_r_nodal_5[0]-0.3535533905932737*Fhat_l_nodal_5[0])*dx10; 
  out_5[3] += (-(0.3535533905932737*Fhat_r_nodal_5[3])+0.3535533905932737*(Fhat_l_nodal_5[3]+Fhat_r_nodal_5[2])-0.3535533905932737*(Fhat_l_nodal_5[2]+Fhat_r_nodal_5[1])+0.3535533905932737*(Fhat_l_nodal_5[1]+Fhat_r_nodal_5[0])-0.3535533905932737*Fhat_l_nodal_5[0])*dx10; 
  out_5[4] += (0.6123724356957944*(Fhat_r_nodal_5[1]+Fhat_l_nodal_5[1]+Fhat_r_nodal_5[0]+Fhat_l_nodal_5[0])-0.6123724356957944*(Fhat_r_nodal_5[3]+Fhat_l_nodal_5[3]+Fhat_r_nodal_5[2]+Fhat_l_nodal_5[2]))*dx10; 
  out_5[5] += (-(0.6123724356957944*(Fhat_r_nodal_5[3]+Fhat_l_nodal_5[3]))+0.6123724356957944*(Fhat_r_nodal_5[2]+Fhat_l_nodal_5[2])-0.6123724356957944*(Fhat_r_nodal_5[1]+Fhat_l_nodal_5[1])+0.6123724356957944*(Fhat_r_nodal_5[0]+Fhat_l_nodal_5[0]))*dx10; 
  out_5[6] += (-(0.3535533905932737*Fhat_r_nodal_5[3])+0.3535533905932737*(Fhat_l_nodal_5[3]+Fhat_r_nodal_5[2])-0.3535533905932737*Fhat_l_nodal_5[2]+0.3535533905932737*Fhat_r_nodal_5[1]-0.3535533905932737*(Fhat_l_nodal_5[1]+Fhat_r_nodal_5[0])+0.3535533905932737*Fhat_l_nodal_5[0])*dx10; 
  out_5[7] += (-(0.6123724356957944*(Fhat_r_nodal_5[3]+Fhat_l_nodal_5[3]))+0.6123724356957944*(Fhat_r_nodal_5[2]+Fhat_l_nodal_5[2]+Fhat_r_nodal_5[1]+Fhat_l_nodal_5[1])-0.6123724356957944*(Fhat_r_nodal_5[0]+Fhat_l_nodal_5[0]))*dx10; 

  const double *Fhat_l_nodal_6 = &Fhat_l_nodal[20]; 
  const double *Fhat_r_nodal_6 = &Fhat_r_nodal[20]; 
  double *out_6 = &out[40]; 
  out_6[0] += (-(0.3535533905932737*Fhat_r_nodal_6[3])+0.3535533905932737*Fhat_l_nodal_6[3]-0.3535533905932737*Fhat_r_nodal_6[2]+0.3535533905932737*Fhat_l_nodal_6[2]-0.3535533905932737*Fhat_r_nodal_6[1]+0.3535533905932737*Fhat_l_nodal_6[1]-0.3535533905932737*Fhat_r_nodal_6[0]+0.3535533905932737*Fhat_l_nodal_6[0])*dx10; 
  out_6[1] += -(0.6123724356957944*(Fhat_r_nodal_6[3]+Fhat_l_nodal_6[3]+Fhat_r_nodal_6[2]+Fhat_l_nodal_6[2]+Fhat_r_nodal_6[1]+Fhat_l_nodal_6[1]+Fhat_r_nodal_6[0]+Fhat_l_nodal_6[0])*dx10); 
  out_6[2] += (-(0.3535533905932737*Fhat_r_nodal_6[3])+0.3535533905932737*Fhat_l_nodal_6[3]-0.3535533905932737*Fhat_r_nodal_6[2]+0.3535533905932737*(Fhat_l_nodal_6[2]+Fhat_r_nodal_6[1])-0.3535533905932737*Fhat_l_nodal_6[1]+0.3535533905932737*Fhat_r_nodal_6[0]-0.3535533905932737*Fhat_l_nodal_6[0])*dx10; 
  out_6[3] += (-(0.3535533905932737*Fhat_r_nodal_6[3])+0.3535533905932737*(Fhat_l_nodal_6[3]+Fhat_r_nodal_6[2])-0.3535533905932737*(Fhat_l_nodal_6[2]+Fhat_r_nodal_6[1])+0.3535533905932737*(Fhat_l_nodal_6[1]+Fhat_r_nodal_6[0])-0.3535533905932737*Fhat_l_nodal_6[0])*dx10; 
  out_6[4] += (0.6123724356957944*(Fhat_r_nodal_6[1]+Fhat_l_nodal_6[1]+Fhat_r_nodal_6[0]+Fhat_l_nodal_6[0])-0.6123724356957944*(Fhat_r_nodal_6[3]+Fhat_l_nodal_6[3]+Fhat_r_nodal_6[2]+Fhat_l_nodal_6[2]))*dx10; 
  out_6[5] += (-(0.6123724356957944*(Fhat_r_nodal_6[3]+Fhat_l_nodal_6[3]))+0.6123724356957944*(Fhat_r_nodal_6[2]+Fhat_l_nodal_6[2])-0.6123724356957944*(Fhat_r_nodal_6[1]+Fhat_l_nodal_6[1])+0.6123724356957944*(Fhat_r_nodal_6[0]+Fhat_l_nodal_6[0]))*dx10; 
  out_6[6] += (-(0.3535533905932737*Fhat_r_nodal_6[3])+0.3535533905932737*(Fhat_l_nodal_6[3]+Fhat_r_nodal_6[2])-0.3535533905932737*Fhat_l_nodal_6[2]+0.3535533905932737*Fhat_r_nodal_6[1]-0.3535533905932737*(Fhat_l_nodal_6[1]+Fhat_r_nodal_6[0])+0.3535533905932737*Fhat_l_nodal_6[0])*dx10; 
  out_6[7] += (-(0.6123724356957944*(Fhat_r_nodal_6[3]+Fhat_l_nodal_6[3]))+0.6123724356957944*(Fhat_r_nodal_6[2]+Fhat_l_nodal_6[2]+Fhat_r_nodal_6[1]+Fhat_l_nodal_6[1])-0.6123724356957944*(Fhat_r_nodal_6[0]+Fhat_l_nodal_6[0]))*dx10; 

  const double *Fhat_l_nodal_7 = &Fhat_l_nodal[24]; 
  const double *Fhat_r_nodal_7 = &Fhat_r_nodal[24]; 
  double *out_7 = &out[48]; 
  out_7[0] += (-(0.3535533905932737*Fhat_r_nodal_7[3])+0.3535533905932737*Fhat_l_nodal_7[3]-0.3535533905932737*Fhat_r_nodal_7[2]+0.3535533905932737*Fhat_l_nodal_7[2]-0.3535533905932737*Fhat_r_nodal_7[1]+0.3535533905932737*Fhat_l_nodal_7[1]-0.3535533905932737*Fhat_r_nodal_7[0]+0.3535533905932737*Fhat_l_nodal_7[0])*dx10; 
  out_7[1] += -(0.6123724356957944*(Fhat_r_nodal_7[3]+Fhat_l_nodal_7[3]+Fhat_r_nodal_7[2]+Fhat_l_nodal_7[2]+Fhat_r_nodal_7[1]+Fhat_l_nodal_7[1]+Fhat_r_nodal_7[0]+Fhat_l_nodal_7[0])*dx10); 
  out_7[2] += (-(0.3535533905932737*Fhat_r_nodal_7[3])+0.3535533905932737*Fhat_l_nodal_7[3]-0.3535533905932737*Fhat_r_nodal_7[2]+0.3535533905932737*(Fhat_l_nodal_7[2]+Fhat_r_nodal_7[1])-0.3535533905932737*Fhat_l_nodal_7[1]+0.3535533905932737*Fhat_r_nodal_7[0]-0.3535533905932737*Fhat_l_nodal_7[0])*dx10; 
  out_7[3] += (-(0.3535533905932737*Fhat_r_nodal_7[3])+0.3535533905932737*(Fhat_l_nodal_7[3]+Fhat_r_nodal_7[2])-0.3535533905932737*(Fhat_l_nodal_7[2]+Fhat_r_nodal_7[1])+0.3535533905932737*(Fhat_l_nodal_7[1]+Fhat_r_nodal_7[0])-0.3535533905932737*Fhat_l_nodal_7[0])*dx10; 
  out_7[4] += (0.6123724356957944*(Fhat_r_nodal_7[1]+Fhat_l_nodal_7[1]+Fhat_r_nodal_7[0]+Fhat_l_nodal_7[0])-0.6123724356957944*(Fhat_r_nodal_7[3]+Fhat_l_nodal_7[3]+Fhat_r_nodal_7[2]+Fhat_l_nodal_7[2]))*dx10; 
  out_7[5] += (-(0.6123724356957944*(Fhat_r_nodal_7[3]+Fhat_l_nodal_7[3]))+0.6123724356957944*(Fhat_r_nodal_7[2]+Fhat_l_nodal_7[2])-0.6123724356957944*(Fhat_r_nodal_7[1]+Fhat_l_nodal_7[1])+0.6123724356957944*(Fhat_r_nodal_7[0]+Fhat_l_nodal_7[0]))*dx10; 
  out_7[6] += (-(0.3535533905932737*Fhat_r_nodal_7[3])+0.3535533905932737*(Fhat_l_nodal_7[3]+Fhat_r_nodal_7[2])-0.3535533905932737*Fhat_l_nodal_7[2]+0.3535533905932737*Fhat_r_nodal_7[1]-0.3535533905932737*(Fhat_l_nodal_7[1]+Fhat_r_nodal_7[0])+0.3535533905932737*Fhat_l_nodal_7[0])*dx10; 
  out_7[7] += (-(0.6123724356957944*(Fhat_r_nodal_7[3]+Fhat_l_nodal_7[3]))+0.6123724356957944*(Fhat_r_nodal_7[2]+Fhat_l_nodal_7[2]+Fhat_r_nodal_7[1]+Fhat_l_nodal_7[1])-0.6123724356957944*(Fhat_r_nodal_7[0]+Fhat_l_nodal_7[0]))*dx10; 

  const double *Fhat_l_nodal_8 = &Fhat_l_nodal[28]; 
  const double *Fhat_r_nodal_8 = &Fhat_r_nodal[28]; 
  double *out_8 = &out[56]; 
  out_8[0] += (-(0.3535533905932737*Fhat_r_nodal_8[3])+0.3535533905932737*Fhat_l_nodal_8[3]-0.3535533905932737*Fhat_r_nodal_8[2]+0.3535533905932737*Fhat_l_nodal_8[2]-0.3535533905932737*Fhat_r_nodal_8[1]+0.3535533905932737*Fhat_l_nodal_8[1]-0.3535533905932737*Fhat_r_nodal_8[0]+0.3535533905932737*Fhat_l_nodal_8[0])*dx10; 
  out_8[1] += -(0.6123724356957944*(Fhat_r_nodal_8[3]+Fhat_l_nodal_8[3]+Fhat_r_nodal_8[2]+Fhat_l_nodal_8[2]+Fhat_r_nodal_8[1]+Fhat_l_nodal_8[1]+Fhat_r_nodal_8[0]+Fhat_l_nodal_8[0])*dx10); 
  out_8[2] += (-(0.3535533905932737*Fhat_r_nodal_8[3])+0.3535533905932737*Fhat_l_nodal_8[3]-0.3535533905932737*Fhat_r_nodal_8[2]+0.3535533905932737*(Fhat_l_nodal_8[2]+Fhat_r_nodal_8[1])-0.3535533905932737*Fhat_l_nodal_8[1]+0.3535533905932737*Fhat_r_nodal_8[0]-0.3535533905932737*Fhat_l_nodal_8[0])*dx10; 
  out_8[3] += (-(0.3535533905932737*Fhat_r_nodal_8[3])+0.3535533905932737*(Fhat_l_nodal_8[3]+Fhat_r_nodal_8[2])-0.3535533905932737*(Fhat_l_nodal_8[2]+Fhat_r_nodal_8[1])+0.3535533905932737*(Fhat_l_nodal_8[1]+Fhat_r_nodal_8[0])-0.3535533905932737*Fhat_l_nodal_8[0])*dx10; 
  out_8[4] += (0.6123724356957944*(Fhat_r_nodal_8[1]+Fhat_l_nodal_8[1]+Fhat_r_nodal_8[0]+Fhat_l_nodal_8[0])-0.6123724356957944*(Fhat_r_nodal_8[3]+Fhat_l_nodal_8[3]+Fhat_r_nodal_8[2]+Fhat_l_nodal_8[2]))*dx10; 
  out_8[5] += (-(0.6123724356957944*(Fhat_r_nodal_8[3]+Fhat_l_nodal_8[3]))+0.6123724356957944*(Fhat_r_nodal_8[2]+Fhat_l_nodal_8[2])-0.6123724356957944*(Fhat_r_nodal_8[1]+Fhat_l_nodal_8[1])+0.6123724356957944*(Fhat_r_nodal_8[0]+Fhat_l_nodal_8[0]))*dx10; 
  out_8[6] += (-(0.3535533905932737*Fhat_r_nodal_8[3])+0.3535533905932737*(Fhat_l_nodal_8[3]+Fhat_r_nodal_8[2])-0.3535533905932737*Fhat_l_nodal_8[2]+0.3535533905932737*Fhat_r_nodal_8[1]-0.3535533905932737*(Fhat_l_nodal_8[1]+Fhat_r_nodal_8[0])+0.3535533905932737*Fhat_l_nodal_8[0])*dx10; 
  out_8[7] += (-(0.6123724356957944*(Fhat_r_nodal_8[3]+Fhat_l_nodal_8[3]))+0.6123724356957944*(Fhat_r_nodal_8[2]+Fhat_l_nodal_8[2]+Fhat_r_nodal_8[1]+Fhat_l_nodal_8[1])-0.6123724356957944*(Fhat_r_nodal_8[0]+Fhat_l_nodal_8[0]))*dx10; 


  return 0.0;

} 
