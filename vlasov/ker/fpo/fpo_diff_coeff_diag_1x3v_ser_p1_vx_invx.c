#include <gkyl_fpo_vlasov_kernels.h> 

GKYL_CU_DH void fpo_diff_coeff_diag_1x3v_vx_ser_p1_invx(const double *dxv, const double *gamma, const double* fpo_g_stencil[3], const double* fpo_d2gdv2_surf, double *diff_coeff) {
  // dxv[NDIM]: Cell spacing in each direction. 
  // gamma: Scalar factor gamma. 
  // fpo_g_stencil[3]: 3 cell stencil of Rosenbluth potential G. 
  // fpo_d2gdv2_surf: Surface projection of d2G/dv2 in center cell. 
  // diff_coeff: Output array for diffusion tensor. 

  // Use cell-average value for gamma. 
  double gamma_avg = gamma[0]/sqrt(pow(2, 1)); 
  double dv1_sq = 4.0/dxv[1]/dxv[1]; 

  const double* G_L = fpo_g_stencil[0]; 
  const double* G_C = fpo_g_stencil[1]; 
  const double* G_R = fpo_g_stencil[2]; 
  
  const double* d2G_surf_C = &fpo_d2gdv2_surf[0]; 
  
  double *out = &diff_coeff[0]; 
  
  out[0] = 0.5625*(G_R[0]*dv1_sq*gamma_avg+G_L[0]*dv1_sq*gamma_avg)-0.5412658773652741*G_R[2]*dv1_sq*gamma_avg+0.5412658773652741*G_L[2]*dv1_sq*gamma_avg-1.125*G_C[0]*dv1_sq*gamma_avg; 
  out[1] = 0.5625*(G_R[1]*dv1_sq*gamma_avg+G_L[1]*dv1_sq*gamma_avg)-0.5412658773652741*G_R[5]*dv1_sq*gamma_avg+0.5412658773652741*G_L[5]*dv1_sq*gamma_avg-1.125*G_C[1]*dv1_sq*gamma_avg; 
  out[2] = -(0.4375*(G_R[2]*dv1_sq*gamma_avg+G_L[2]*dv1_sq*gamma_avg))-2.875*G_C[2]*dv1_sq*gamma_avg+0.5412658773652741*G_R[0]*dv1_sq*gamma_avg-0.5412658773652741*G_L[0]*dv1_sq*gamma_avg; 
  out[3] = 0.5625*(G_R[3]*dv1_sq*gamma_avg+G_L[3]*dv1_sq*gamma_avg)-0.5412658773652741*G_R[7]*dv1_sq*gamma_avg+0.5412658773652741*G_L[7]*dv1_sq*gamma_avg-1.125*G_C[3]*dv1_sq*gamma_avg; 
  out[4] = 0.5625*(G_R[4]*dv1_sq*gamma_avg+G_L[4]*dv1_sq*gamma_avg)-0.5412658773652741*G_R[9]*dv1_sq*gamma_avg+0.5412658773652741*G_L[9]*dv1_sq*gamma_avg-1.125*G_C[4]*dv1_sq*gamma_avg; 
  out[5] = -(0.4375*(G_R[5]*dv1_sq*gamma_avg+G_L[5]*dv1_sq*gamma_avg))-2.875*G_C[5]*dv1_sq*gamma_avg+0.5412658773652741*G_R[1]*dv1_sq*gamma_avg-0.5412658773652741*G_L[1]*dv1_sq*gamma_avg; 
  out[6] = 0.5625*(G_R[6]*dv1_sq*gamma_avg+G_L[6]*dv1_sq*gamma_avg)-0.5412658773652741*G_R[11]*dv1_sq*gamma_avg+0.5412658773652741*G_L[11]*dv1_sq*gamma_avg-1.125*G_C[6]*dv1_sq*gamma_avg; 
  out[7] = -(0.4375*(G_R[7]*dv1_sq*gamma_avg+G_L[7]*dv1_sq*gamma_avg))-2.875*G_C[7]*dv1_sq*gamma_avg+0.5412658773652741*G_R[3]*dv1_sq*gamma_avg-0.5412658773652741*G_L[3]*dv1_sq*gamma_avg; 
  out[8] = 0.5625*(G_R[8]*dv1_sq*gamma_avg+G_L[8]*dv1_sq*gamma_avg)-0.5412658773652741*G_R[12]*dv1_sq*gamma_avg+0.5412658773652741*G_L[12]*dv1_sq*gamma_avg-1.125*G_C[8]*dv1_sq*gamma_avg; 
  out[9] = -(0.4375*(G_R[9]*dv1_sq*gamma_avg+G_L[9]*dv1_sq*gamma_avg))-2.875*G_C[9]*dv1_sq*gamma_avg+0.5412658773652741*G_R[4]*dv1_sq*gamma_avg-0.5412658773652741*G_L[4]*dv1_sq*gamma_avg; 
  out[10] = 0.5625*(G_R[10]*dv1_sq*gamma_avg+G_L[10]*dv1_sq*gamma_avg)-0.5412658773652741*G_R[14]*dv1_sq*gamma_avg+0.5412658773652741*G_L[14]*dv1_sq*gamma_avg-1.125*G_C[10]*dv1_sq*gamma_avg; 
  out[11] = -(0.4375*(G_R[11]*dv1_sq*gamma_avg+G_L[11]*dv1_sq*gamma_avg))-2.875*G_C[11]*dv1_sq*gamma_avg+0.5412658773652741*G_R[6]*dv1_sq*gamma_avg-0.5412658773652741*G_L[6]*dv1_sq*gamma_avg; 
  out[12] = -(0.4375*(G_R[12]*dv1_sq*gamma_avg+G_L[12]*dv1_sq*gamma_avg))-2.875*G_C[12]*dv1_sq*gamma_avg+0.5412658773652741*G_R[8]*dv1_sq*gamma_avg-0.5412658773652741*G_L[8]*dv1_sq*gamma_avg; 
  out[13] = 0.5625*(G_R[13]*dv1_sq*gamma_avg+G_L[13]*dv1_sq*gamma_avg)-0.5412658773652741*G_R[15]*dv1_sq*gamma_avg+0.5412658773652741*G_L[15]*dv1_sq*gamma_avg-1.125*G_C[13]*dv1_sq*gamma_avg; 
  out[14] = -(0.4375*(G_R[14]*dv1_sq*gamma_avg+G_L[14]*dv1_sq*gamma_avg))-2.875*G_C[14]*dv1_sq*gamma_avg+0.5412658773652741*G_R[10]*dv1_sq*gamma_avg-0.5412658773652741*G_L[10]*dv1_sq*gamma_avg; 
  out[15] = -(0.4375*(G_R[15]*dv1_sq*gamma_avg+G_L[15]*dv1_sq*gamma_avg))-2.875*G_C[15]*dv1_sq*gamma_avg+0.5412658773652741*G_R[13]*dv1_sq*gamma_avg-0.5412658773652741*G_L[13]*dv1_sq*gamma_avg; 
} 

