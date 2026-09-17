#include <gkyl_fpo_vlasov_kernels.h> 

GKYL_CU_DH void fpo_diff_coeff_diag_2x3v_vy_ser_p1_upvy(const double *dxv, const double *gamma, const double* fpo_g_stencil[3], const double* fpo_d2gdv2_surf, double *diff_coeff) {
  // dxv[NDIM]: Cell spacing in each direction. 
  // gamma: Scalar factor gamma. 
  // fpo_g_stencil[3]: 3 cell stencil of Rosenbluth potential G. 
  // fpo_d2gdv2_surf: Surface projection of d2G/dv2 in center cell. 
  // diff_coeff: Output array for diffusion tensor. 

  // Use cell-average value for gamma. 
  double gamma_avg = gamma[0]/sqrt(pow(2, 2)); 
  double dv1_sq = 4.0/dxv[3]/dxv[3]; 

  const double* G_L = fpo_g_stencil[0]; 
  const double* G_C = fpo_g_stencil[1]; 
  const double* G_R = fpo_g_stencil[2]; 
  
  const double* d2G_surf_C = &fpo_d2gdv2_surf[16]; 
  
  double *out = &diff_coeff[128]; 
  
  out[0] = 1.4433756729740645*G_L[4]*dv1_sq*gamma_avg+2.886751345948129*G_C[4]*dv1_sq*gamma_avg+1.25*G_L[0]*dv1_sq*gamma_avg-1.25*G_C[0]*dv1_sq*gamma_avg+0.2357022603955158*d2G_surf_C[0]*gamma_avg; 
  out[1] = 1.4433756729740645*G_L[9]*dv1_sq*gamma_avg+2.886751345948129*G_C[9]*dv1_sq*gamma_avg+1.25*G_L[1]*dv1_sq*gamma_avg-1.25*G_C[1]*dv1_sq*gamma_avg+0.2357022603955158*d2G_surf_C[1]*gamma_avg; 
  out[2] = 1.4433756729740645*G_L[10]*dv1_sq*gamma_avg+2.886751345948129*G_C[10]*dv1_sq*gamma_avg+1.25*G_L[2]*dv1_sq*gamma_avg-1.25*G_C[2]*dv1_sq*gamma_avg+0.2357022603955158*d2G_surf_C[2]*gamma_avg; 
  out[3] = 1.4433756729740645*G_L[11]*dv1_sq*gamma_avg+2.886751345948129*G_C[11]*dv1_sq*gamma_avg+1.25*G_L[3]*dv1_sq*gamma_avg-1.25*G_C[3]*dv1_sq*gamma_avg+0.2357022603955158*d2G_surf_C[3]*gamma_avg; 
  out[4] = 0.5208333333333334*G_L[4]*dv1_sq*gamma_avg+0.10416666666666667*G_C[4]*dv1_sq*gamma_avg+0.18042195912175807*G_L[0]*dv1_sq*gamma_avg-0.18042195912175807*G_C[0]*dv1_sq*gamma_avg+0.34020690871988585*d2G_surf_C[0]*gamma_avg; 
  out[5] = 1.4433756729740645*G_L[15]*dv1_sq*gamma_avg+2.886751345948129*G_C[15]*dv1_sq*gamma_avg+1.25*G_L[5]*dv1_sq*gamma_avg-1.25*G_C[5]*dv1_sq*gamma_avg+0.2357022603955158*d2G_surf_C[4]*gamma_avg; 
  out[6] = 1.4433756729740645*G_L[17]*dv1_sq*gamma_avg+2.886751345948129*G_C[17]*dv1_sq*gamma_avg+1.25*G_L[6]*dv1_sq*gamma_avg-1.25*G_C[6]*dv1_sq*gamma_avg+0.2357022603955158*d2G_surf_C[5]*gamma_avg; 
  out[7] = 1.4433756729740645*G_L[18]*dv1_sq*gamma_avg+2.886751345948129*G_C[18]*dv1_sq*gamma_avg+1.25*G_L[7]*dv1_sq*gamma_avg-1.25*G_C[7]*dv1_sq*gamma_avg+0.2357022603955158*d2G_surf_C[6]*gamma_avg; 
  out[8] = 1.4433756729740645*G_L[19]*dv1_sq*gamma_avg+2.886751345948129*G_C[19]*dv1_sq*gamma_avg+1.25*G_L[8]*dv1_sq*gamma_avg-1.25*G_C[8]*dv1_sq*gamma_avg+0.2357022603955158*d2G_surf_C[7]*gamma_avg; 
  out[9] = 0.5208333333333334*G_L[9]*dv1_sq*gamma_avg+0.10416666666666667*G_C[9]*dv1_sq*gamma_avg+0.18042195912175807*G_L[1]*dv1_sq*gamma_avg-0.18042195912175807*G_C[1]*dv1_sq*gamma_avg+0.34020690871988585*d2G_surf_C[1]*gamma_avg; 
  out[10] = 0.5208333333333334*G_L[10]*dv1_sq*gamma_avg+0.10416666666666667*G_C[10]*dv1_sq*gamma_avg+0.18042195912175807*G_L[2]*dv1_sq*gamma_avg-0.18042195912175807*G_C[2]*dv1_sq*gamma_avg+0.34020690871988585*d2G_surf_C[2]*gamma_avg; 
  out[11] = 0.5208333333333334*G_L[11]*dv1_sq*gamma_avg+0.10416666666666667*G_C[11]*dv1_sq*gamma_avg+0.18042195912175807*G_L[3]*dv1_sq*gamma_avg-0.18042195912175807*G_C[3]*dv1_sq*gamma_avg+0.34020690871988585*d2G_surf_C[3]*gamma_avg; 
  out[12] = 1.4433756729740645*G_L[23]*dv1_sq*gamma_avg+2.886751345948129*G_C[23]*dv1_sq*gamma_avg+1.25*G_L[12]*dv1_sq*gamma_avg-1.25*G_C[12]*dv1_sq*gamma_avg+0.2357022603955158*d2G_surf_C[8]*gamma_avg; 
  out[13] = 1.4433756729740645*G_L[24]*dv1_sq*gamma_avg+2.886751345948129*G_C[24]*dv1_sq*gamma_avg+1.25*G_L[13]*dv1_sq*gamma_avg-1.25*G_C[13]*dv1_sq*gamma_avg+0.2357022603955158*d2G_surf_C[9]*gamma_avg; 
  out[14] = 1.4433756729740645*G_L[25]*dv1_sq*gamma_avg+2.886751345948129*G_C[25]*dv1_sq*gamma_avg+1.25*G_L[14]*dv1_sq*gamma_avg-1.25*G_C[14]*dv1_sq*gamma_avg+0.2357022603955158*d2G_surf_C[10]*gamma_avg; 
  out[15] = 0.5208333333333334*G_L[15]*dv1_sq*gamma_avg+0.10416666666666667*G_C[15]*dv1_sq*gamma_avg+0.18042195912175807*G_L[5]*dv1_sq*gamma_avg-0.18042195912175807*G_C[5]*dv1_sq*gamma_avg+0.34020690871988585*d2G_surf_C[4]*gamma_avg; 
  out[16] = 1.4433756729740645*G_L[26]*dv1_sq*gamma_avg+2.886751345948129*G_C[26]*dv1_sq*gamma_avg+1.25*G_L[16]*dv1_sq*gamma_avg-1.25*G_C[16]*dv1_sq*gamma_avg+0.2357022603955158*d2G_surf_C[11]*gamma_avg; 
  out[17] = 0.5208333333333334*G_L[17]*dv1_sq*gamma_avg+0.10416666666666667*G_C[17]*dv1_sq*gamma_avg+0.18042195912175807*G_L[6]*dv1_sq*gamma_avg-0.18042195912175807*G_C[6]*dv1_sq*gamma_avg+0.34020690871988585*d2G_surf_C[5]*gamma_avg; 
  out[18] = 0.5208333333333334*G_L[18]*dv1_sq*gamma_avg+0.10416666666666667*G_C[18]*dv1_sq*gamma_avg+0.18042195912175807*G_L[7]*dv1_sq*gamma_avg-0.18042195912175807*G_C[7]*dv1_sq*gamma_avg+0.34020690871988585*d2G_surf_C[6]*gamma_avg; 
  out[19] = 0.5208333333333334*G_L[19]*dv1_sq*gamma_avg+0.10416666666666667*G_C[19]*dv1_sq*gamma_avg+0.18042195912175807*G_L[8]*dv1_sq*gamma_avg-0.18042195912175807*G_C[8]*dv1_sq*gamma_avg+0.34020690871988585*d2G_surf_C[7]*gamma_avg; 
  out[20] = 1.4433756729740645*G_L[28]*dv1_sq*gamma_avg+2.886751345948129*G_C[28]*dv1_sq*gamma_avg+1.25*G_L[20]*dv1_sq*gamma_avg-1.25*G_C[20]*dv1_sq*gamma_avg+0.2357022603955158*d2G_surf_C[12]*gamma_avg; 
  out[21] = 1.4433756729740645*G_L[29]*dv1_sq*gamma_avg+2.886751345948129*G_C[29]*dv1_sq*gamma_avg+1.25*G_L[21]*dv1_sq*gamma_avg-1.25*G_C[21]*dv1_sq*gamma_avg+0.2357022603955158*d2G_surf_C[13]*gamma_avg; 
  out[22] = 1.4433756729740645*G_L[30]*dv1_sq*gamma_avg+2.886751345948129*G_C[30]*dv1_sq*gamma_avg+1.25*G_L[22]*dv1_sq*gamma_avg-1.25*G_C[22]*dv1_sq*gamma_avg+0.2357022603955158*d2G_surf_C[14]*gamma_avg; 
  out[23] = 0.5208333333333334*G_L[23]*dv1_sq*gamma_avg+0.10416666666666667*G_C[23]*dv1_sq*gamma_avg+0.18042195912175807*G_L[12]*dv1_sq*gamma_avg-0.18042195912175807*G_C[12]*dv1_sq*gamma_avg+0.34020690871988585*d2G_surf_C[8]*gamma_avg; 
  out[24] = 0.5208333333333334*G_L[24]*dv1_sq*gamma_avg+0.10416666666666667*G_C[24]*dv1_sq*gamma_avg+0.18042195912175807*G_L[13]*dv1_sq*gamma_avg-0.18042195912175807*G_C[13]*dv1_sq*gamma_avg+0.34020690871988585*d2G_surf_C[9]*gamma_avg; 
  out[25] = 0.5208333333333334*G_L[25]*dv1_sq*gamma_avg+0.10416666666666667*G_C[25]*dv1_sq*gamma_avg+0.18042195912175807*G_L[14]*dv1_sq*gamma_avg-0.18042195912175807*G_C[14]*dv1_sq*gamma_avg+0.34020690871988585*d2G_surf_C[10]*gamma_avg; 
  out[26] = 0.5208333333333334*G_L[26]*dv1_sq*gamma_avg+0.10416666666666667*G_C[26]*dv1_sq*gamma_avg+0.18042195912175807*G_L[16]*dv1_sq*gamma_avg-0.18042195912175807*G_C[16]*dv1_sq*gamma_avg+0.34020690871988585*d2G_surf_C[11]*gamma_avg; 
  out[27] = 1.4433756729740645*G_L[31]*dv1_sq*gamma_avg+2.886751345948129*G_C[31]*dv1_sq*gamma_avg+1.25*G_L[27]*dv1_sq*gamma_avg-1.25*G_C[27]*dv1_sq*gamma_avg+0.2357022603955158*d2G_surf_C[15]*gamma_avg; 
  out[28] = 0.5208333333333334*G_L[28]*dv1_sq*gamma_avg+0.10416666666666667*G_C[28]*dv1_sq*gamma_avg+0.18042195912175807*G_L[20]*dv1_sq*gamma_avg-0.18042195912175807*G_C[20]*dv1_sq*gamma_avg+0.34020690871988585*d2G_surf_C[12]*gamma_avg; 
  out[29] = 0.5208333333333334*G_L[29]*dv1_sq*gamma_avg+0.10416666666666667*G_C[29]*dv1_sq*gamma_avg+0.18042195912175807*G_L[21]*dv1_sq*gamma_avg-0.18042195912175807*G_C[21]*dv1_sq*gamma_avg+0.34020690871988585*d2G_surf_C[13]*gamma_avg; 
  out[30] = 0.5208333333333334*G_L[30]*dv1_sq*gamma_avg+0.10416666666666667*G_C[30]*dv1_sq*gamma_avg+0.18042195912175807*G_L[22]*dv1_sq*gamma_avg-0.18042195912175807*G_C[22]*dv1_sq*gamma_avg+0.34020690871988585*d2G_surf_C[14]*gamma_avg; 
  out[31] = 0.5208333333333334*G_L[31]*dv1_sq*gamma_avg+0.10416666666666667*G_C[31]*dv1_sq*gamma_avg+0.18042195912175807*G_L[27]*dv1_sq*gamma_avg-0.18042195912175807*G_C[27]*dv1_sq*gamma_avg+0.34020690871988585*d2G_surf_C[15]*gamma_avg; 
} 

