#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH double gr_maxwell_surfx_1x_ser_p1(const double *w, const double *dxv,
  const double *Fhat_l_nodal, const double *Fhat_r_nodal, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 

  const double *Fhat_l_nodal_1 = &Fhat_l_nodal[0]; 
  const double *Fhat_r_nodal_1 = &Fhat_r_nodal[0]; 
  double *out_1 = &out[0]; 
  out_1[0] += (0.7071067811865475*Fhat_l_nodal_1[0]-0.7071067811865475*Fhat_r_nodal_1[0])*dx10; 
  out_1[1] += -(1.224744871391589*(Fhat_r_nodal_1[0]+Fhat_l_nodal_1[0])*dx10); 

  const double *Fhat_l_nodal_2 = &Fhat_l_nodal[1]; 
  const double *Fhat_r_nodal_2 = &Fhat_r_nodal[1]; 
  double *out_2 = &out[2]; 
  out_2[0] += (0.7071067811865475*Fhat_l_nodal_2[0]-0.7071067811865475*Fhat_r_nodal_2[0])*dx10; 
  out_2[1] += -(1.224744871391589*(Fhat_r_nodal_2[0]+Fhat_l_nodal_2[0])*dx10); 

  const double *Fhat_l_nodal_3 = &Fhat_l_nodal[2]; 
  const double *Fhat_r_nodal_3 = &Fhat_r_nodal[2]; 
  double *out_3 = &out[4]; 
  out_3[0] += (0.7071067811865475*Fhat_l_nodal_3[0]-0.7071067811865475*Fhat_r_nodal_3[0])*dx10; 
  out_3[1] += -(1.224744871391589*(Fhat_r_nodal_3[0]+Fhat_l_nodal_3[0])*dx10); 

  const double *Fhat_l_nodal_4 = &Fhat_l_nodal[3]; 
  const double *Fhat_r_nodal_4 = &Fhat_r_nodal[3]; 
  double *out_4 = &out[6]; 
  out_4[0] += (0.7071067811865475*Fhat_l_nodal_4[0]-0.7071067811865475*Fhat_r_nodal_4[0])*dx10; 
  out_4[1] += -(1.224744871391589*(Fhat_r_nodal_4[0]+Fhat_l_nodal_4[0])*dx10); 

  const double *Fhat_l_nodal_5 = &Fhat_l_nodal[4]; 
  const double *Fhat_r_nodal_5 = &Fhat_r_nodal[4]; 
  double *out_5 = &out[8]; 
  out_5[0] += (0.7071067811865475*Fhat_l_nodal_5[0]-0.7071067811865475*Fhat_r_nodal_5[0])*dx10; 
  out_5[1] += -(1.224744871391589*(Fhat_r_nodal_5[0]+Fhat_l_nodal_5[0])*dx10); 

  const double *Fhat_l_nodal_6 = &Fhat_l_nodal[5]; 
  const double *Fhat_r_nodal_6 = &Fhat_r_nodal[5]; 
  double *out_6 = &out[10]; 
  out_6[0] += (0.7071067811865475*Fhat_l_nodal_6[0]-0.7071067811865475*Fhat_r_nodal_6[0])*dx10; 
  out_6[1] += -(1.224744871391589*(Fhat_r_nodal_6[0]+Fhat_l_nodal_6[0])*dx10); 

  const double *Fhat_l_nodal_7 = &Fhat_l_nodal[6]; 
  const double *Fhat_r_nodal_7 = &Fhat_r_nodal[6]; 
  double *out_7 = &out[12]; 
  out_7[0] += (0.7071067811865475*Fhat_l_nodal_7[0]-0.7071067811865475*Fhat_r_nodal_7[0])*dx10; 
  out_7[1] += -(1.224744871391589*(Fhat_r_nodal_7[0]+Fhat_l_nodal_7[0])*dx10); 

  const double *Fhat_l_nodal_8 = &Fhat_l_nodal[7]; 
  const double *Fhat_r_nodal_8 = &Fhat_r_nodal[7]; 
  double *out_8 = &out[14]; 
  out_8[0] += (0.7071067811865475*Fhat_l_nodal_8[0]-0.7071067811865475*Fhat_r_nodal_8[0])*dx10; 
  out_8[1] += -(1.224744871391589*(Fhat_r_nodal_8[0]+Fhat_l_nodal_8[0])*dx10); 


  return 0.0;

} 
