#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH double gr_maxwell_surfx_2x_ser_p1(const double *w, const double *dxv,
  const double *Fhat_l_nodal, const double *Fhat_r_nodal, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 

  const double *Fhat_l_nodal_1 = &Fhat_l_nodal[0]; 
  const double *Fhat_r_nodal_1 = &Fhat_r_nodal[0]; 
  double *out_1 = &out[0]; 
  out_1[0] += (-(0.5*Fhat_r_nodal_1[1])+0.5*Fhat_l_nodal_1[1]-0.5*Fhat_r_nodal_1[0]+0.5*Fhat_l_nodal_1[0])*dx10; 
  out_1[1] += -(0.8660254037844386*(Fhat_r_nodal_1[1]+Fhat_l_nodal_1[1]+Fhat_r_nodal_1[0]+Fhat_l_nodal_1[0])*dx10); 
  out_1[2] += (-(0.5*Fhat_r_nodal_1[1])+0.5*(Fhat_l_nodal_1[1]+Fhat_r_nodal_1[0])-0.5*Fhat_l_nodal_1[0])*dx10; 
  out_1[3] += (0.8660254037844386*(Fhat_r_nodal_1[0]+Fhat_l_nodal_1[0])-0.8660254037844386*(Fhat_r_nodal_1[1]+Fhat_l_nodal_1[1]))*dx10; 

  const double *Fhat_l_nodal_2 = &Fhat_l_nodal[2]; 
  const double *Fhat_r_nodal_2 = &Fhat_r_nodal[2]; 
  double *out_2 = &out[4]; 
  out_2[0] += (-(0.5*Fhat_r_nodal_2[1])+0.5*Fhat_l_nodal_2[1]-0.5*Fhat_r_nodal_2[0]+0.5*Fhat_l_nodal_2[0])*dx10; 
  out_2[1] += -(0.8660254037844386*(Fhat_r_nodal_2[1]+Fhat_l_nodal_2[1]+Fhat_r_nodal_2[0]+Fhat_l_nodal_2[0])*dx10); 
  out_2[2] += (-(0.5*Fhat_r_nodal_2[1])+0.5*(Fhat_l_nodal_2[1]+Fhat_r_nodal_2[0])-0.5*Fhat_l_nodal_2[0])*dx10; 
  out_2[3] += (0.8660254037844386*(Fhat_r_nodal_2[0]+Fhat_l_nodal_2[0])-0.8660254037844386*(Fhat_r_nodal_2[1]+Fhat_l_nodal_2[1]))*dx10; 

  const double *Fhat_l_nodal_3 = &Fhat_l_nodal[4]; 
  const double *Fhat_r_nodal_3 = &Fhat_r_nodal[4]; 
  double *out_3 = &out[8]; 
  out_3[0] += (-(0.5*Fhat_r_nodal_3[1])+0.5*Fhat_l_nodal_3[1]-0.5*Fhat_r_nodal_3[0]+0.5*Fhat_l_nodal_3[0])*dx10; 
  out_3[1] += -(0.8660254037844386*(Fhat_r_nodal_3[1]+Fhat_l_nodal_3[1]+Fhat_r_nodal_3[0]+Fhat_l_nodal_3[0])*dx10); 
  out_3[2] += (-(0.5*Fhat_r_nodal_3[1])+0.5*(Fhat_l_nodal_3[1]+Fhat_r_nodal_3[0])-0.5*Fhat_l_nodal_3[0])*dx10; 
  out_3[3] += (0.8660254037844386*(Fhat_r_nodal_3[0]+Fhat_l_nodal_3[0])-0.8660254037844386*(Fhat_r_nodal_3[1]+Fhat_l_nodal_3[1]))*dx10; 

  const double *Fhat_l_nodal_4 = &Fhat_l_nodal[6]; 
  const double *Fhat_r_nodal_4 = &Fhat_r_nodal[6]; 
  double *out_4 = &out[12]; 
  out_4[0] += (-(0.5*Fhat_r_nodal_4[1])+0.5*Fhat_l_nodal_4[1]-0.5*Fhat_r_nodal_4[0]+0.5*Fhat_l_nodal_4[0])*dx10; 
  out_4[1] += -(0.8660254037844386*(Fhat_r_nodal_4[1]+Fhat_l_nodal_4[1]+Fhat_r_nodal_4[0]+Fhat_l_nodal_4[0])*dx10); 
  out_4[2] += (-(0.5*Fhat_r_nodal_4[1])+0.5*(Fhat_l_nodal_4[1]+Fhat_r_nodal_4[0])-0.5*Fhat_l_nodal_4[0])*dx10; 
  out_4[3] += (0.8660254037844386*(Fhat_r_nodal_4[0]+Fhat_l_nodal_4[0])-0.8660254037844386*(Fhat_r_nodal_4[1]+Fhat_l_nodal_4[1]))*dx10; 

  const double *Fhat_l_nodal_5 = &Fhat_l_nodal[8]; 
  const double *Fhat_r_nodal_5 = &Fhat_r_nodal[8]; 
  double *out_5 = &out[16]; 
  out_5[0] += (-(0.5*Fhat_r_nodal_5[1])+0.5*Fhat_l_nodal_5[1]-0.5*Fhat_r_nodal_5[0]+0.5*Fhat_l_nodal_5[0])*dx10; 
  out_5[1] += -(0.8660254037844386*(Fhat_r_nodal_5[1]+Fhat_l_nodal_5[1]+Fhat_r_nodal_5[0]+Fhat_l_nodal_5[0])*dx10); 
  out_5[2] += (-(0.5*Fhat_r_nodal_5[1])+0.5*(Fhat_l_nodal_5[1]+Fhat_r_nodal_5[0])-0.5*Fhat_l_nodal_5[0])*dx10; 
  out_5[3] += (0.8660254037844386*(Fhat_r_nodal_5[0]+Fhat_l_nodal_5[0])-0.8660254037844386*(Fhat_r_nodal_5[1]+Fhat_l_nodal_5[1]))*dx10; 

  const double *Fhat_l_nodal_6 = &Fhat_l_nodal[10]; 
  const double *Fhat_r_nodal_6 = &Fhat_r_nodal[10]; 
  double *out_6 = &out[20]; 
  out_6[0] += (-(0.5*Fhat_r_nodal_6[1])+0.5*Fhat_l_nodal_6[1]-0.5*Fhat_r_nodal_6[0]+0.5*Fhat_l_nodal_6[0])*dx10; 
  out_6[1] += -(0.8660254037844386*(Fhat_r_nodal_6[1]+Fhat_l_nodal_6[1]+Fhat_r_nodal_6[0]+Fhat_l_nodal_6[0])*dx10); 
  out_6[2] += (-(0.5*Fhat_r_nodal_6[1])+0.5*(Fhat_l_nodal_6[1]+Fhat_r_nodal_6[0])-0.5*Fhat_l_nodal_6[0])*dx10; 
  out_6[3] += (0.8660254037844386*(Fhat_r_nodal_6[0]+Fhat_l_nodal_6[0])-0.8660254037844386*(Fhat_r_nodal_6[1]+Fhat_l_nodal_6[1]))*dx10; 

  const double *Fhat_l_nodal_7 = &Fhat_l_nodal[12]; 
  const double *Fhat_r_nodal_7 = &Fhat_r_nodal[12]; 
  double *out_7 = &out[24]; 
  out_7[0] += (-(0.5*Fhat_r_nodal_7[1])+0.5*Fhat_l_nodal_7[1]-0.5*Fhat_r_nodal_7[0]+0.5*Fhat_l_nodal_7[0])*dx10; 
  out_7[1] += -(0.8660254037844386*(Fhat_r_nodal_7[1]+Fhat_l_nodal_7[1]+Fhat_r_nodal_7[0]+Fhat_l_nodal_7[0])*dx10); 
  out_7[2] += (-(0.5*Fhat_r_nodal_7[1])+0.5*(Fhat_l_nodal_7[1]+Fhat_r_nodal_7[0])-0.5*Fhat_l_nodal_7[0])*dx10; 
  out_7[3] += (0.8660254037844386*(Fhat_r_nodal_7[0]+Fhat_l_nodal_7[0])-0.8660254037844386*(Fhat_r_nodal_7[1]+Fhat_l_nodal_7[1]))*dx10; 

  const double *Fhat_l_nodal_8 = &Fhat_l_nodal[14]; 
  const double *Fhat_r_nodal_8 = &Fhat_r_nodal[14]; 
  double *out_8 = &out[28]; 
  out_8[0] += (-(0.5*Fhat_r_nodal_8[1])+0.5*Fhat_l_nodal_8[1]-0.5*Fhat_r_nodal_8[0]+0.5*Fhat_l_nodal_8[0])*dx10; 
  out_8[1] += -(0.8660254037844386*(Fhat_r_nodal_8[1]+Fhat_l_nodal_8[1]+Fhat_r_nodal_8[0]+Fhat_l_nodal_8[0])*dx10); 
  out_8[2] += (-(0.5*Fhat_r_nodal_8[1])+0.5*(Fhat_l_nodal_8[1]+Fhat_r_nodal_8[0])-0.5*Fhat_l_nodal_8[0])*dx10; 
  out_8[3] += (0.8660254037844386*(Fhat_r_nodal_8[0]+Fhat_l_nodal_8[0])-0.8660254037844386*(Fhat_r_nodal_8[1]+Fhat_l_nodal_8[1]))*dx10; 


  return 0.0;

} 
