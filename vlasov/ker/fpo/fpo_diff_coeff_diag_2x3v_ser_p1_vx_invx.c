#include <gkyl_fpo_vlasov_kernels.h> 

GKYL_CU_DH void fpo_diff_coeff_diag_2x3v_vx_ser_p1_invx(const double *dxv, const double *gamma, const double* fpo_g_stencil[3], const double* fpo_d2gdv2_surf, double *diff_coeff) {
  // dxv[NDIM]: Cell spacing in each direction. 
  // gamma: Scalar factor gamma. 
  // fpo_g_stencil[3]: 3 cell stencil of Rosenbluth potential G. 
  // fpo_d2gdv2_surf: Surface projection of d2G/dv2 in center cell. 
  // diff_coeff: Output array for diffusion tensor. 

  // Use cell-average value for gamma. 
  double gamma_avg = gamma[0]/sqrt(pow(2, 2)); 
  double dv1_sq = 4.0/dxv[2]/dxv[2]; 

  const double* G_L = fpo_g_stencil[0]; 
  const double* G_C = fpo_g_stencil[1]; 
  const double* G_R = fpo_g_stencil[2]; 
  
  const double* d2G_surf_C = &fpo_d2gdv2_surf[0]; 
  
  double *out = &diff_coeff[0]; 
  
  out[0] = 0.5625*(G_R[0]*dv1_sq*gamma_avg+G_L[0]*dv1_sq*gamma_avg)-0.5412658773652741*G_R[3]*dv1_sq*gamma_avg+0.5412658773652741*G_L[3]*dv1_sq*gamma_avg-1.125*G_C[0]*dv1_sq*gamma_avg; 
  out[1] = 0.5625*(G_R[1]*dv1_sq*gamma_avg+G_L[1]*dv1_sq*gamma_avg)-0.5412658773652741*G_R[7]*dv1_sq*gamma_avg+0.5412658773652741*G_L[7]*dv1_sq*gamma_avg-1.125*G_C[1]*dv1_sq*gamma_avg; 
  out[2] = 0.5625*(G_R[2]*dv1_sq*gamma_avg+G_L[2]*dv1_sq*gamma_avg)-0.5412658773652741*G_R[8]*dv1_sq*gamma_avg+0.5412658773652741*G_L[8]*dv1_sq*gamma_avg-1.125*G_C[2]*dv1_sq*gamma_avg; 
  out[3] = -(0.4375*(G_R[3]*dv1_sq*gamma_avg+G_L[3]*dv1_sq*gamma_avg))-2.875*G_C[3]*dv1_sq*gamma_avg+0.5412658773652741*G_R[0]*dv1_sq*gamma_avg-0.5412658773652741*G_L[0]*dv1_sq*gamma_avg; 
  out[4] = 0.5625*(G_R[4]*dv1_sq*gamma_avg+G_L[4]*dv1_sq*gamma_avg)-0.5412658773652741*G_R[11]*dv1_sq*gamma_avg+0.5412658773652741*G_L[11]*dv1_sq*gamma_avg-1.125*G_C[4]*dv1_sq*gamma_avg; 
  out[5] = 0.5625*(G_R[5]*dv1_sq*gamma_avg+G_L[5]*dv1_sq*gamma_avg)-0.5412658773652741*G_R[14]*dv1_sq*gamma_avg+0.5412658773652741*G_L[14]*dv1_sq*gamma_avg-1.125*G_C[5]*dv1_sq*gamma_avg; 
  out[6] = 0.5625*(G_R[6]*dv1_sq*gamma_avg+G_L[6]*dv1_sq*gamma_avg)-0.5412658773652741*G_R[16]*dv1_sq*gamma_avg+0.5412658773652741*G_L[16]*dv1_sq*gamma_avg-1.125*G_C[6]*dv1_sq*gamma_avg; 
  out[7] = -(0.4375*(G_R[7]*dv1_sq*gamma_avg+G_L[7]*dv1_sq*gamma_avg))-2.875*G_C[7]*dv1_sq*gamma_avg+0.5412658773652741*G_R[1]*dv1_sq*gamma_avg-0.5412658773652741*G_L[1]*dv1_sq*gamma_avg; 
  out[8] = -(0.4375*(G_R[8]*dv1_sq*gamma_avg+G_L[8]*dv1_sq*gamma_avg))-2.875*G_C[8]*dv1_sq*gamma_avg+0.5412658773652741*G_R[2]*dv1_sq*gamma_avg-0.5412658773652741*G_L[2]*dv1_sq*gamma_avg; 
  out[9] = 0.5625*(G_R[9]*dv1_sq*gamma_avg+G_L[9]*dv1_sq*gamma_avg)-0.5412658773652741*G_R[18]*dv1_sq*gamma_avg+0.5412658773652741*G_L[18]*dv1_sq*gamma_avg-1.125*G_C[9]*dv1_sq*gamma_avg; 
  out[10] = 0.5625*(G_R[10]*dv1_sq*gamma_avg+G_L[10]*dv1_sq*gamma_avg)-0.5412658773652741*G_R[19]*dv1_sq*gamma_avg+0.5412658773652741*G_L[19]*dv1_sq*gamma_avg-1.125*G_C[10]*dv1_sq*gamma_avg; 
  out[11] = -(0.4375*(G_R[11]*dv1_sq*gamma_avg+G_L[11]*dv1_sq*gamma_avg))-2.875*G_C[11]*dv1_sq*gamma_avg+0.5412658773652741*G_R[4]*dv1_sq*gamma_avg-0.5412658773652741*G_L[4]*dv1_sq*gamma_avg; 
  out[12] = 0.5625*(G_R[12]*dv1_sq*gamma_avg+G_L[12]*dv1_sq*gamma_avg)-0.5412658773652741*G_R[21]*dv1_sq*gamma_avg+0.5412658773652741*G_L[21]*dv1_sq*gamma_avg-1.125*G_C[12]*dv1_sq*gamma_avg; 
  out[13] = 0.5625*(G_R[13]*dv1_sq*gamma_avg+G_L[13]*dv1_sq*gamma_avg)-0.5412658773652741*G_R[22]*dv1_sq*gamma_avg+0.5412658773652741*G_L[22]*dv1_sq*gamma_avg-1.125*G_C[13]*dv1_sq*gamma_avg; 
  out[14] = -(0.4375*(G_R[14]*dv1_sq*gamma_avg+G_L[14]*dv1_sq*gamma_avg))-2.875*G_C[14]*dv1_sq*gamma_avg+0.5412658773652741*G_R[5]*dv1_sq*gamma_avg-0.5412658773652741*G_L[5]*dv1_sq*gamma_avg; 
  out[15] = 0.5625*(G_R[15]*dv1_sq*gamma_avg+G_L[15]*dv1_sq*gamma_avg)-0.5412658773652741*G_R[25]*dv1_sq*gamma_avg+0.5412658773652741*G_L[25]*dv1_sq*gamma_avg-1.125*G_C[15]*dv1_sq*gamma_avg; 
  out[16] = -(0.4375*(G_R[16]*dv1_sq*gamma_avg+G_L[16]*dv1_sq*gamma_avg))-2.875*G_C[16]*dv1_sq*gamma_avg+0.5412658773652741*G_R[6]*dv1_sq*gamma_avg-0.5412658773652741*G_L[6]*dv1_sq*gamma_avg; 
  out[17] = 0.5625*(G_R[17]*dv1_sq*gamma_avg+G_L[17]*dv1_sq*gamma_avg)-0.5412658773652741*G_R[26]*dv1_sq*gamma_avg+0.5412658773652741*G_L[26]*dv1_sq*gamma_avg-1.125*G_C[17]*dv1_sq*gamma_avg; 
  out[18] = -(0.4375*(G_R[18]*dv1_sq*gamma_avg+G_L[18]*dv1_sq*gamma_avg))-2.875*G_C[18]*dv1_sq*gamma_avg+0.5412658773652741*G_R[9]*dv1_sq*gamma_avg-0.5412658773652741*G_L[9]*dv1_sq*gamma_avg; 
  out[19] = -(0.4375*(G_R[19]*dv1_sq*gamma_avg+G_L[19]*dv1_sq*gamma_avg))-2.875*G_C[19]*dv1_sq*gamma_avg+0.5412658773652741*G_R[10]*dv1_sq*gamma_avg-0.5412658773652741*G_L[10]*dv1_sq*gamma_avg; 
  out[20] = 0.5625*(G_R[20]*dv1_sq*gamma_avg+G_L[20]*dv1_sq*gamma_avg)-0.5412658773652741*G_R[27]*dv1_sq*gamma_avg+0.5412658773652741*G_L[27]*dv1_sq*gamma_avg-1.125*G_C[20]*dv1_sq*gamma_avg; 
  out[21] = -(0.4375*(G_R[21]*dv1_sq*gamma_avg+G_L[21]*dv1_sq*gamma_avg))-2.875*G_C[21]*dv1_sq*gamma_avg+0.5412658773652741*G_R[12]*dv1_sq*gamma_avg-0.5412658773652741*G_L[12]*dv1_sq*gamma_avg; 
  out[22] = -(0.4375*(G_R[22]*dv1_sq*gamma_avg+G_L[22]*dv1_sq*gamma_avg))-2.875*G_C[22]*dv1_sq*gamma_avg+0.5412658773652741*G_R[13]*dv1_sq*gamma_avg-0.5412658773652741*G_L[13]*dv1_sq*gamma_avg; 
  out[23] = 0.5625*(G_R[23]*dv1_sq*gamma_avg+G_L[23]*dv1_sq*gamma_avg)-0.5412658773652741*G_R[29]*dv1_sq*gamma_avg+0.5412658773652741*G_L[29]*dv1_sq*gamma_avg-1.125*G_C[23]*dv1_sq*gamma_avg; 
  out[24] = 0.5625*(G_R[24]*dv1_sq*gamma_avg+G_L[24]*dv1_sq*gamma_avg)-0.5412658773652741*G_R[30]*dv1_sq*gamma_avg+0.5412658773652741*G_L[30]*dv1_sq*gamma_avg-1.125*G_C[24]*dv1_sq*gamma_avg; 
  out[25] = -(0.4375*(G_R[25]*dv1_sq*gamma_avg+G_L[25]*dv1_sq*gamma_avg))-2.875*G_C[25]*dv1_sq*gamma_avg+0.5412658773652741*G_R[15]*dv1_sq*gamma_avg-0.5412658773652741*G_L[15]*dv1_sq*gamma_avg; 
  out[26] = -(0.4375*(G_R[26]*dv1_sq*gamma_avg+G_L[26]*dv1_sq*gamma_avg))-2.875*G_C[26]*dv1_sq*gamma_avg+0.5412658773652741*G_R[17]*dv1_sq*gamma_avg-0.5412658773652741*G_L[17]*dv1_sq*gamma_avg; 
  out[27] = -(0.4375*(G_R[27]*dv1_sq*gamma_avg+G_L[27]*dv1_sq*gamma_avg))-2.875*G_C[27]*dv1_sq*gamma_avg+0.5412658773652741*G_R[20]*dv1_sq*gamma_avg-0.5412658773652741*G_L[20]*dv1_sq*gamma_avg; 
  out[28] = 0.5625*(G_R[28]*dv1_sq*gamma_avg+G_L[28]*dv1_sq*gamma_avg)-0.5412658773652741*G_R[31]*dv1_sq*gamma_avg+0.5412658773652741*G_L[31]*dv1_sq*gamma_avg-1.125*G_C[28]*dv1_sq*gamma_avg; 
  out[29] = -(0.4375*(G_R[29]*dv1_sq*gamma_avg+G_L[29]*dv1_sq*gamma_avg))-2.875*G_C[29]*dv1_sq*gamma_avg+0.5412658773652741*G_R[23]*dv1_sq*gamma_avg-0.5412658773652741*G_L[23]*dv1_sq*gamma_avg; 
  out[30] = -(0.4375*(G_R[30]*dv1_sq*gamma_avg+G_L[30]*dv1_sq*gamma_avg))-2.875*G_C[30]*dv1_sq*gamma_avg+0.5412658773652741*G_R[24]*dv1_sq*gamma_avg-0.5412658773652741*G_L[24]*dv1_sq*gamma_avg; 
  out[31] = -(0.4375*(G_R[31]*dv1_sq*gamma_avg+G_L[31]*dv1_sq*gamma_avg))-2.875*G_C[31]*dv1_sq*gamma_avg+0.5412658773652741*G_R[28]*dv1_sq*gamma_avg-0.5412658773652741*G_L[28]*dv1_sq*gamma_avg; 
} 

