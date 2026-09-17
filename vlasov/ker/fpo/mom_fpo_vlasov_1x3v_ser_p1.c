#include <gkyl_mom_fpo_vlasov_kernels.h> 
GKYL_CU_DH void mom_fpo_vlasov_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *a_i, const double *D_ij,  const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[4]:   cell center in each direction. 
  // dxv[4]: cell length in each direction. 
  // idx[4]: cell index in each direction. 
  // a_i:     Drag coefficient. 
  // D_ij:    Diffusion tensor. 
  // f[16]:   Input distribution function. 
  // out:     Output moments a + div(D) (3-components), a . v + div(D . v). 
 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  const double *ax = &a_i[0]; 
  const double *ay = &a_i[16]; 
  const double *az = &a_i[32]; 
  const double *Dxx = &D_ij[0]; 
  const double *Dxy = &D_ij[16]; 
  const double *Dxz = &D_ij[32]; 
  const double *Dyx = &D_ij[48]; 
  const double *Dyy = &D_ij[64]; 
  const double *Dyz = &D_ij[80]; 
  const double *Dzx = &D_ij[96]; 
  const double *Dzy = &D_ij[112]; 
  const double *Dzz = &D_ij[128]; 

  double *ax_corr = &out[0]; 
  double *ay_corr = &out[2]; 
  double *az_corr = &out[4]; 
  double *energy_corr = &out[6]; 

  double temp_1[2] = {0.0}; 
  double temp_2[2] = {0.0}; 
  double temp_3[2] = {0.0}; 
  double temp_4[2] = {0.0}; 
  double temp_5[2] = {0.0}; 
  double temp_6[2] = {0.0}; 
  double temp_7[2] = {0.0}; 
  double temp_8[2] = {0.0}; 
  double temp_9[2] = {0.0}; 
  double temp_10[2] = {0.0}; 
  double temp_11[2] = {0.0}; 
  double temp_12[2] = {0.0}; 
  temp_1[0] = 0.7071067811865475*(ax[15]*f[15]+ax[14]*f[14]+ax[13]*f[13]+ax[12]*f[12]+ax[11]*f[11]+ax[10]*f[10]+ax[9]*f[9]+ax[8]*f[8]+ax[7]*f[7]+ax[6]*f[6]+ax[5]*f[5]+ax[4]*f[4]+ax[3]*f[3]+ax[2]*f[2]+ax[1]*f[1]+ax[0]*f[0]); 
  temp_1[1] = 0.7071067811865475*(ax[14]*f[15]+f[14]*ax[15]+ax[10]*f[13]+f[10]*ax[13]+ax[9]*f[12]+f[9]*ax[12]+ax[7]*f[11]+f[7]*ax[11]+ax[4]*f[8]+f[4]*ax[8]+ax[3]*f[6]+f[3]*ax[6]+ax[2]*f[5]+f[2]*ax[5]+ax[0]*f[1]+f[0]*ax[1]); 
  temp_2[0] = (2.4494897427831783*(f[13]*Dxx[15]+f[10]*Dxx[14]+f[8]*Dxx[12]+f[6]*Dxx[11]+f[4]*Dxx[9]+f[3]*Dxx[7]+f[1]*Dxx[5]+f[0]*Dxx[2]))/dv1; 
  temp_2[1] = (2.4494897427831783*(f[10]*Dxx[15]+f[13]*Dxx[14]+f[4]*Dxx[12]+f[3]*Dxx[11]+f[8]*Dxx[9]+f[6]*Dxx[7]+f[0]*Dxx[5]+f[1]*Dxx[2]))/dv1; 
  temp_3[0] = (2.4494897427831783*(f[12]*Dyx[15]+f[9]*Dyx[14]+f[8]*Dyx[13]+f[5]*Dyx[11]+f[4]*Dyx[10]+f[2]*Dyx[7]+f[1]*Dyx[6]+f[0]*Dyx[3]))/dv2; 
  temp_3[1] = (2.4494897427831783*(f[9]*Dyx[15]+f[12]*Dyx[14]+f[4]*Dyx[13]+f[2]*Dyx[11]+f[8]*Dyx[10]+f[5]*Dyx[7]+f[0]*Dyx[6]+f[1]*Dyx[3]))/dv2; 
  temp_4[0] = (2.4494897427831783*(f[11]*Dzx[15]+f[7]*Dzx[14]+f[6]*Dzx[13]+f[5]*Dzx[12]+f[3]*Dzx[10]+f[2]*Dzx[9]+f[1]*Dzx[8]+f[0]*Dzx[4]))/dv3; 
  temp_4[1] = (2.4494897427831783*(f[7]*Dzx[15]+f[11]*Dzx[14]+f[3]*Dzx[13]+f[2]*Dzx[12]+f[6]*Dzx[10]+f[5]*Dzx[9]+f[0]*Dzx[8]+f[1]*Dzx[4]))/dv3; 
  ax_corr[0] += (temp_4[0]+temp_3[0]+temp_2[0]+temp_1[0])*volFact; 
  ax_corr[1] += (temp_4[1]+temp_3[1]+temp_2[1]+temp_1[1])*volFact; 

  temp_1[0] = 0.7071067811865475*(ay[15]*f[15]+ay[14]*f[14]+ay[13]*f[13]+ay[12]*f[12]+ay[11]*f[11]+ay[10]*f[10]+ay[9]*f[9]+ay[8]*f[8]+ay[7]*f[7]+ay[6]*f[6]+ay[5]*f[5]+ay[4]*f[4]+ay[3]*f[3]+ay[2]*f[2]+ay[1]*f[1]+ay[0]*f[0]); 
  temp_1[1] = 0.7071067811865475*(ay[14]*f[15]+f[14]*ay[15]+ay[10]*f[13]+f[10]*ay[13]+ay[9]*f[12]+f[9]*ay[12]+ay[7]*f[11]+f[7]*ay[11]+ay[4]*f[8]+f[4]*ay[8]+ay[3]*f[6]+f[3]*ay[6]+ay[2]*f[5]+f[2]*ay[5]+ay[0]*f[1]+f[0]*ay[1]); 
  temp_2[0] = (2.4494897427831783*(f[13]*Dxy[15]+f[10]*Dxy[14]+f[8]*Dxy[12]+f[6]*Dxy[11]+f[4]*Dxy[9]+f[3]*Dxy[7]+f[1]*Dxy[5]+f[0]*Dxy[2]))/dv1; 
  temp_2[1] = (2.4494897427831783*(f[10]*Dxy[15]+f[13]*Dxy[14]+f[4]*Dxy[12]+f[3]*Dxy[11]+f[8]*Dxy[9]+f[6]*Dxy[7]+f[0]*Dxy[5]+f[1]*Dxy[2]))/dv1; 
  temp_3[0] = (2.4494897427831783*(f[12]*Dyy[15]+f[9]*Dyy[14]+f[8]*Dyy[13]+f[5]*Dyy[11]+f[4]*Dyy[10]+f[2]*Dyy[7]+f[1]*Dyy[6]+f[0]*Dyy[3]))/dv2; 
  temp_3[1] = (2.4494897427831783*(f[9]*Dyy[15]+f[12]*Dyy[14]+f[4]*Dyy[13]+f[2]*Dyy[11]+f[8]*Dyy[10]+f[5]*Dyy[7]+f[0]*Dyy[6]+f[1]*Dyy[3]))/dv2; 
  temp_4[0] = (2.4494897427831783*(f[11]*Dzy[15]+f[7]*Dzy[14]+f[6]*Dzy[13]+f[5]*Dzy[12]+f[3]*Dzy[10]+f[2]*Dzy[9]+f[1]*Dzy[8]+f[0]*Dzy[4]))/dv3; 
  temp_4[1] = (2.4494897427831783*(f[7]*Dzy[15]+f[11]*Dzy[14]+f[3]*Dzy[13]+f[2]*Dzy[12]+f[6]*Dzy[10]+f[5]*Dzy[9]+f[0]*Dzy[8]+f[1]*Dzy[4]))/dv3; 
  ay_corr[0] += (temp_4[0]+temp_3[0]+temp_2[0]+temp_1[0])*volFact; 
  ay_corr[1] += (temp_4[1]+temp_3[1]+temp_2[1]+temp_1[1])*volFact; 

  temp_1[0] = 0.7071067811865475*(az[15]*f[15]+az[14]*f[14]+az[13]*f[13]+az[12]*f[12]+az[11]*f[11]+az[10]*f[10]+az[9]*f[9]+az[8]*f[8]+az[7]*f[7]+az[6]*f[6]+az[5]*f[5]+az[4]*f[4]+az[3]*f[3]+az[2]*f[2]+az[1]*f[1]+az[0]*f[0]); 
  temp_1[1] = 0.7071067811865475*(az[14]*f[15]+f[14]*az[15]+az[10]*f[13]+f[10]*az[13]+az[9]*f[12]+f[9]*az[12]+az[7]*f[11]+f[7]*az[11]+az[4]*f[8]+f[4]*az[8]+az[3]*f[6]+f[3]*az[6]+az[2]*f[5]+f[2]*az[5]+az[0]*f[1]+f[0]*az[1]); 
  temp_2[0] = (2.4494897427831783*(f[13]*Dxz[15]+f[10]*Dxz[14]+f[8]*Dxz[12]+f[6]*Dxz[11]+f[4]*Dxz[9]+f[3]*Dxz[7]+f[1]*Dxz[5]+f[0]*Dxz[2]))/dv1; 
  temp_2[1] = (2.4494897427831783*(f[10]*Dxz[15]+f[13]*Dxz[14]+f[4]*Dxz[12]+f[3]*Dxz[11]+f[8]*Dxz[9]+f[6]*Dxz[7]+f[0]*Dxz[5]+f[1]*Dxz[2]))/dv1; 
  temp_3[0] = (2.4494897427831783*(f[12]*Dyz[15]+f[9]*Dyz[14]+f[8]*Dyz[13]+f[5]*Dyz[11]+f[4]*Dyz[10]+f[2]*Dyz[7]+f[1]*Dyz[6]+f[0]*Dyz[3]))/dv2; 
  temp_3[1] = (2.4494897427831783*(f[9]*Dyz[15]+f[12]*Dyz[14]+f[4]*Dyz[13]+f[2]*Dyz[11]+f[8]*Dyz[10]+f[5]*Dyz[7]+f[0]*Dyz[6]+f[1]*Dyz[3]))/dv2; 
  temp_4[0] = (2.4494897427831783*(f[11]*Dzz[15]+f[7]*Dzz[14]+f[6]*Dzz[13]+f[5]*Dzz[12]+f[3]*Dzz[10]+f[2]*Dzz[9]+f[1]*Dzz[8]+f[0]*Dzz[4]))/dv3; 
  temp_4[1] = (2.4494897427831783*(f[7]*Dzz[15]+f[11]*Dzz[14]+f[3]*Dzz[13]+f[2]*Dzz[12]+f[6]*Dzz[10]+f[5]*Dzz[9]+f[0]*Dzz[8]+f[1]*Dzz[4]))/dv3; 
  az_corr[0] += (temp_4[0]+temp_3[0]+temp_2[0]+temp_1[0])*volFact; 
  az_corr[1] += (temp_4[1]+temp_3[1]+temp_2[1]+temp_1[1])*volFact; 

  temp_1[0] = 0.7071067811865475*(ax[15]*f[15]+ax[14]*f[14]+ax[13]*f[13]+ax[12]*f[12]+ax[11]*f[11]+ax[10]*f[10]+ax[9]*f[9]+ax[8]*f[8]+ax[7]*f[7]+ax[6]*f[6]+ax[5]*f[5]+ax[4]*f[4]+ax[3]*f[3]+ax[2]*f[2]+ax[1]*f[1]+ax[0]*f[0])*wx1+0.20412414523193148*(ax[13]*f[15]+f[13]*ax[15]+ax[10]*f[14]+f[10]*ax[14]+ax[8]*f[12]+f[8]*ax[12]+ax[6]*f[11]+f[6]*ax[11]+ax[4]*f[9]+f[4]*ax[9]+ax[3]*f[7]+f[3]*ax[7]+ax[1]*f[5]+f[1]*ax[5]+ax[0]*f[2]+f[0]*ax[2])*dv1; 
  temp_1[1] = 0.7071067811865475*(ax[14]*f[15]+f[14]*ax[15]+ax[10]*f[13]+f[10]*ax[13]+ax[9]*f[12]+f[9]*ax[12]+ax[7]*f[11]+f[7]*ax[11]+ax[4]*f[8]+f[4]*ax[8]+ax[3]*f[6]+f[3]*ax[6]+ax[2]*f[5]+f[2]*ax[5]+ax[0]*f[1]+f[0]*ax[1])*wx1+0.20412414523193148*(ax[10]*f[15]+f[10]*ax[15]+ax[13]*f[14]+f[13]*ax[14]+ax[4]*f[12]+f[4]*ax[12]+ax[3]*f[11]+f[3]*ax[11]+ax[8]*f[9]+f[8]*ax[9]+ax[6]*f[7]+f[6]*ax[7]+ax[0]*f[5]+f[0]*ax[5]+ax[1]*f[2]+f[1]*ax[2])*dv1; 
  temp_2[0] = 0.7071067811865475*(ay[15]*f[15]+ay[14]*f[14]+ay[13]*f[13]+ay[12]*f[12]+ay[11]*f[11]+ay[10]*f[10]+ay[9]*f[9]+ay[8]*f[8]+ay[7]*f[7]+ay[6]*f[6]+ay[5]*f[5]+ay[4]*f[4]+ay[3]*f[3]+ay[2]*f[2]+ay[1]*f[1]+ay[0]*f[0])*wx2+0.20412414523193148*(ay[12]*f[15]+f[12]*ay[15]+ay[9]*f[14]+f[9]*ay[14]+ay[8]*f[13]+f[8]*ay[13]+ay[5]*f[11]+f[5]*ay[11]+ay[4]*f[10]+f[4]*ay[10]+ay[2]*f[7]+f[2]*ay[7]+ay[1]*f[6]+f[1]*ay[6]+ay[0]*f[3]+f[0]*ay[3])*dv2; 
  temp_2[1] = 0.7071067811865475*(ay[14]*f[15]+f[14]*ay[15]+ay[10]*f[13]+f[10]*ay[13]+ay[9]*f[12]+f[9]*ay[12]+ay[7]*f[11]+f[7]*ay[11]+ay[4]*f[8]+f[4]*ay[8]+ay[3]*f[6]+f[3]*ay[6]+ay[2]*f[5]+f[2]*ay[5]+ay[0]*f[1]+f[0]*ay[1])*wx2+0.20412414523193148*(ay[9]*f[15]+f[9]*ay[15]+ay[12]*f[14]+f[12]*ay[14]+ay[4]*f[13]+f[4]*ay[13]+ay[2]*f[11]+f[2]*ay[11]+ay[8]*f[10]+f[8]*ay[10]+ay[5]*f[7]+f[5]*ay[7]+ay[0]*f[6]+f[0]*ay[6]+ay[1]*f[3]+f[1]*ay[3])*dv2; 
  temp_3[0] = 0.7071067811865475*(az[15]*f[15]+az[14]*f[14]+az[13]*f[13]+az[12]*f[12]+az[11]*f[11]+az[10]*f[10]+az[9]*f[9]+az[8]*f[8]+az[7]*f[7]+az[6]*f[6]+az[5]*f[5]+az[4]*f[4]+az[3]*f[3]+az[2]*f[2]+az[1]*f[1]+az[0]*f[0])*wx3+0.20412414523193148*(az[11]*f[15]+f[11]*az[15]+az[7]*f[14]+f[7]*az[14]+az[6]*f[13]+f[6]*az[13]+az[5]*f[12]+f[5]*az[12]+az[3]*f[10]+f[3]*az[10]+az[2]*f[9]+f[2]*az[9]+az[1]*f[8]+f[1]*az[8]+az[0]*f[4]+f[0]*az[4])*dv3; 
  temp_3[1] = 0.7071067811865475*(az[14]*f[15]+f[14]*az[15]+az[10]*f[13]+f[10]*az[13]+az[9]*f[12]+f[9]*az[12]+az[7]*f[11]+f[7]*az[11]+az[4]*f[8]+f[4]*az[8]+az[3]*f[6]+f[3]*az[6]+az[2]*f[5]+f[2]*az[5]+az[0]*f[1]+f[0]*az[1])*wx3+0.20412414523193148*(az[7]*f[15]+f[7]*az[15]+az[11]*f[14]+f[11]*az[14]+az[3]*f[13]+f[3]*az[13]+az[2]*f[12]+f[2]*az[12]+az[6]*f[10]+f[6]*az[10]+az[5]*f[9]+f[5]*az[9]+az[0]*f[8]+f[0]*az[8]+az[1]*f[4]+f[1]*az[4])*dv3; 
  temp_4[0] = (2.4494897427831783*(f[13]*Dxx[15]+f[10]*Dxx[14]+f[8]*Dxx[12]+f[6]*Dxx[11]+f[4]*Dxx[9]+f[3]*Dxx[7]+f[1]*Dxx[5]+f[0]*Dxx[2])*wx1)/dv1+1.4142135623730951*(Dxx[15]*f[15]+Dxx[14]*f[14]+Dxx[12]*f[12]+Dxx[11]*f[11]+Dxx[9]*f[9]+Dxx[7]*f[7]+Dxx[5]*f[5]+Dxx[2]*f[2])+0.7071067811865475*(Dxx[13]*f[13]+Dxx[10]*f[10]+Dxx[8]*f[8]+Dxx[6]*f[6]+Dxx[4]*f[4]+Dxx[3]*f[3]+Dxx[1]*f[1]+Dxx[0]*f[0]); 
  temp_4[1] = (2.4494897427831783*(f[10]*Dxx[15]+f[13]*Dxx[14]+f[4]*Dxx[12]+f[3]*Dxx[11]+f[8]*Dxx[9]+f[6]*Dxx[7]+f[0]*Dxx[5]+f[1]*Dxx[2])*wx1)/dv1+1.4142135623730951*(Dxx[14]*f[15]+f[14]*Dxx[15]+Dxx[9]*f[12]+f[9]*Dxx[12]+Dxx[7]*f[11]+f[7]*Dxx[11]+Dxx[2]*f[5]+f[2]*Dxx[5])+0.7071067811865475*(Dxx[10]*f[13]+f[10]*Dxx[13]+Dxx[4]*f[8]+f[4]*Dxx[8]+Dxx[3]*f[6]+f[3]*Dxx[6]+Dxx[0]*f[1]+f[0]*Dxx[1]); 
  temp_5[0] = (2.4494897427831783*(f[13]*Dxy[15]+f[10]*Dxy[14]+f[8]*Dxy[12]+f[6]*Dxy[11]+f[4]*Dxy[9]+f[3]*Dxy[7]+f[1]*Dxy[5]+f[0]*Dxy[2])*wx2+0.7071067811865475*(f[8]*Dxy[15]+f[4]*Dxy[14]+Dxy[12]*f[13]+f[1]*Dxy[11]+Dxy[9]*f[10]+f[0]*Dxy[7]+Dxy[5]*f[6]+Dxy[2]*f[3])*dv2)/dv1; 
  temp_5[1] = (2.4494897427831783*(f[10]*Dxy[15]+f[13]*Dxy[14]+f[4]*Dxy[12]+f[3]*Dxy[11]+f[8]*Dxy[9]+f[6]*Dxy[7]+f[0]*Dxy[5]+f[1]*Dxy[2])*wx2+0.7071067811865475*(f[4]*Dxy[15]+f[8]*Dxy[14]+Dxy[9]*f[13]+f[10]*Dxy[12]+f[0]*Dxy[11]+f[1]*Dxy[7]+Dxy[2]*f[6]+f[3]*Dxy[5])*dv2)/dv1; 
  temp_6[0] = (2.4494897427831783*(f[13]*Dxz[15]+f[10]*Dxz[14]+f[8]*Dxz[12]+f[6]*Dxz[11]+f[4]*Dxz[9]+f[3]*Dxz[7]+f[1]*Dxz[5]+f[0]*Dxz[2])*wx3+0.7071067811865475*(f[6]*Dxz[15]+f[3]*Dxz[14]+Dxz[11]*f[13]+f[1]*Dxz[12]+Dxz[7]*f[10]+f[0]*Dxz[9]+Dxz[5]*f[8]+Dxz[2]*f[4])*dv3)/dv1; 
  temp_6[1] = (2.4494897427831783*(f[10]*Dxz[15]+f[13]*Dxz[14]+f[4]*Dxz[12]+f[3]*Dxz[11]+f[8]*Dxz[9]+f[6]*Dxz[7]+f[0]*Dxz[5]+f[1]*Dxz[2])*wx3+0.7071067811865475*(f[3]*Dxz[15]+f[6]*Dxz[14]+Dxz[7]*f[13]+f[0]*Dxz[12]+f[10]*Dxz[11]+f[1]*Dxz[9]+Dxz[2]*f[8]+f[4]*Dxz[5])*dv3)/dv1; 
  temp_7[0] = (2.4494897427831783*(f[12]*Dyx[15]+f[9]*Dyx[14]+f[8]*Dyx[13]+f[5]*Dyx[11]+f[4]*Dyx[10]+f[2]*Dyx[7]+f[1]*Dyx[6]+f[0]*Dyx[3])*wx1+0.7071067811865475*(f[8]*Dyx[15]+f[4]*Dyx[14]+f[12]*Dyx[13]+f[1]*Dyx[11]+f[9]*Dyx[10]+f[0]*Dyx[7]+f[5]*Dyx[6]+f[2]*Dyx[3])*dv1)/dv2; 
  temp_7[1] = (2.4494897427831783*(f[9]*Dyx[15]+f[12]*Dyx[14]+f[4]*Dyx[13]+f[2]*Dyx[11]+f[8]*Dyx[10]+f[5]*Dyx[7]+f[0]*Dyx[6]+f[1]*Dyx[3])*wx1+0.7071067811865475*(f[4]*Dyx[15]+f[8]*Dyx[14]+f[9]*Dyx[13]+Dyx[10]*f[12]+f[0]*Dyx[11]+f[1]*Dyx[7]+f[2]*Dyx[6]+Dyx[3]*f[5])*dv1)/dv2; 
  temp_8[0] = (2.4494897427831783*(f[12]*Dyy[15]+f[9]*Dyy[14]+f[8]*Dyy[13]+f[5]*Dyy[11]+f[4]*Dyy[10]+f[2]*Dyy[7]+f[1]*Dyy[6]+f[0]*Dyy[3])*wx2)/dv2+1.4142135623730951*(Dyy[15]*f[15]+Dyy[14]*f[14]+Dyy[13]*f[13]+Dyy[11]*f[11]+Dyy[10]*f[10]+Dyy[7]*f[7]+Dyy[6]*f[6]+Dyy[3]*f[3])+0.7071067811865475*(Dyy[12]*f[12]+Dyy[9]*f[9]+Dyy[8]*f[8]+Dyy[5]*f[5]+Dyy[4]*f[4]+Dyy[2]*f[2]+Dyy[1]*f[1]+Dyy[0]*f[0]); 
  temp_8[1] = (2.4494897427831783*(f[9]*Dyy[15]+f[12]*Dyy[14]+f[4]*Dyy[13]+f[2]*Dyy[11]+f[8]*Dyy[10]+f[5]*Dyy[7]+f[0]*Dyy[6]+f[1]*Dyy[3])*wx2)/dv2+1.4142135623730951*(Dyy[14]*f[15]+f[14]*Dyy[15]+Dyy[10]*f[13]+f[10]*Dyy[13]+Dyy[7]*f[11]+f[7]*Dyy[11]+Dyy[3]*f[6]+f[3]*Dyy[6])+0.7071067811865475*(Dyy[9]*f[12]+f[9]*Dyy[12]+Dyy[4]*f[8]+f[4]*Dyy[8]+Dyy[2]*f[5]+f[2]*Dyy[5]+Dyy[0]*f[1]+f[0]*Dyy[1]); 
  temp_9[0] = (2.4494897427831783*(f[12]*Dyz[15]+f[9]*Dyz[14]+f[8]*Dyz[13]+f[5]*Dyz[11]+f[4]*Dyz[10]+f[2]*Dyz[7]+f[1]*Dyz[6]+f[0]*Dyz[3])*wx3+0.7071067811865475*(f[5]*Dyz[15]+f[2]*Dyz[14]+f[1]*Dyz[13]+Dyz[11]*f[12]+f[0]*Dyz[10]+Dyz[7]*f[9]+Dyz[6]*f[8]+Dyz[3]*f[4])*dv3)/dv2; 
  temp_9[1] = (2.4494897427831783*(f[9]*Dyz[15]+f[12]*Dyz[14]+f[4]*Dyz[13]+f[2]*Dyz[11]+f[8]*Dyz[10]+f[5]*Dyz[7]+f[0]*Dyz[6]+f[1]*Dyz[3])*wx3+0.7071067811865475*(f[2]*Dyz[15]+f[5]*Dyz[14]+f[0]*Dyz[13]+Dyz[7]*f[12]+f[9]*Dyz[11]+f[1]*Dyz[10]+Dyz[3]*f[8]+f[4]*Dyz[6])*dv3)/dv2; 
  temp_10[0] = (2.4494897427831783*(f[11]*Dzx[15]+f[7]*Dzx[14]+f[6]*Dzx[13]+f[5]*Dzx[12]+f[3]*Dzx[10]+f[2]*Dzx[9]+f[1]*Dzx[8]+f[0]*Dzx[4])*wx1+0.7071067811865475*(f[6]*Dzx[15]+f[3]*Dzx[14]+f[11]*Dzx[13]+f[1]*Dzx[12]+f[7]*Dzx[10]+f[0]*Dzx[9]+f[5]*Dzx[8]+f[2]*Dzx[4])*dv1)/dv3; 
  temp_10[1] = (2.4494897427831783*(f[7]*Dzx[15]+f[11]*Dzx[14]+f[3]*Dzx[13]+f[2]*Dzx[12]+f[6]*Dzx[10]+f[5]*Dzx[9]+f[0]*Dzx[8]+f[1]*Dzx[4])*wx1+0.7071067811865475*(f[3]*Dzx[15]+f[6]*Dzx[14]+f[7]*Dzx[13]+f[0]*Dzx[12]+Dzx[10]*f[11]+f[1]*Dzx[9]+f[2]*Dzx[8]+Dzx[4]*f[5])*dv1)/dv3; 
  temp_11[0] = (2.4494897427831783*(f[11]*Dzy[15]+f[7]*Dzy[14]+f[6]*Dzy[13]+f[5]*Dzy[12]+f[3]*Dzy[10]+f[2]*Dzy[9]+f[1]*Dzy[8]+f[0]*Dzy[4])*wx2+0.7071067811865475*(f[5]*Dzy[15]+f[2]*Dzy[14]+f[1]*Dzy[13]+f[11]*Dzy[12]+f[0]*Dzy[10]+f[7]*Dzy[9]+f[6]*Dzy[8]+f[3]*Dzy[4])*dv2)/dv3; 
  temp_11[1] = (2.4494897427831783*(f[7]*Dzy[15]+f[11]*Dzy[14]+f[3]*Dzy[13]+f[2]*Dzy[12]+f[6]*Dzy[10]+f[5]*Dzy[9]+f[0]*Dzy[8]+f[1]*Dzy[4])*wx2+0.7071067811865475*(f[2]*Dzy[15]+f[5]*Dzy[14]+f[0]*Dzy[13]+f[7]*Dzy[12]+Dzy[9]*f[11]+f[1]*Dzy[10]+f[3]*Dzy[8]+Dzy[4]*f[6])*dv2)/dv3; 
  temp_12[0] = (2.4494897427831783*(f[11]*Dzz[15]+f[7]*Dzz[14]+f[6]*Dzz[13]+f[5]*Dzz[12]+f[3]*Dzz[10]+f[2]*Dzz[9]+f[1]*Dzz[8]+f[0]*Dzz[4])*wx3)/dv3+1.4142135623730951*(Dzz[15]*f[15]+Dzz[14]*f[14]+Dzz[13]*f[13]+Dzz[12]*f[12]+Dzz[10]*f[10]+Dzz[9]*f[9]+Dzz[8]*f[8]+Dzz[4]*f[4])+0.7071067811865475*(Dzz[11]*f[11]+Dzz[7]*f[7]+Dzz[6]*f[6]+Dzz[5]*f[5]+Dzz[3]*f[3]+Dzz[2]*f[2]+Dzz[1]*f[1]+Dzz[0]*f[0]); 
  temp_12[1] = (2.4494897427831783*(f[7]*Dzz[15]+f[11]*Dzz[14]+f[3]*Dzz[13]+f[2]*Dzz[12]+f[6]*Dzz[10]+f[5]*Dzz[9]+f[0]*Dzz[8]+f[1]*Dzz[4])*wx3)/dv3+1.4142135623730951*(Dzz[14]*f[15]+f[14]*Dzz[15]+Dzz[10]*f[13]+f[10]*Dzz[13]+Dzz[9]*f[12]+f[9]*Dzz[12]+Dzz[4]*f[8]+f[4]*Dzz[8])+0.7071067811865475*(Dzz[7]*f[11]+f[7]*Dzz[11]+Dzz[3]*f[6]+f[3]*Dzz[6]+Dzz[2]*f[5]+f[2]*Dzz[5]+Dzz[0]*f[1]+f[0]*Dzz[1]); 
  energy_corr[0] += (temp_9[0]+temp_8[0]+temp_7[0]+temp_6[0]+temp_5[0]+temp_4[0]+temp_3[0]+temp_2[0]+temp_12[0]+temp_11[0]+temp_10[0]+temp_1[0])*volFact; 
  energy_corr[1] += (temp_9[1]+temp_8[1]+temp_7[1]+temp_6[1]+temp_5[1]+temp_4[1]+temp_3[1]+temp_2[1]+temp_12[1]+temp_11[1]+temp_10[1]+temp_1[1])*volFact; 

} 
