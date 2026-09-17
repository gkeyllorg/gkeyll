#include <gkyl_fpo_vlasov_kernels.h> 

GKYL_CU_DH void fpo_diff_coeff_diag_1x3v_vz_ser_p2_invz(const double *dxv, const double *gamma, const double* fpo_g_stencil[3], const double* fpo_d2gdv2_surf, double *diff_coeff) {
  // dxv[NDIM]: Cell spacing in each direction. 
  // gamma: Scalar factor gamma. 
  // fpo_g_stencil[3]: 3 cell stencil of Rosenbluth potential G. 
  // fpo_d2gdv2_surf: Surface projection of d2G/dv2 in center cell. 
  // diff_coeff: Output array for diffusion tensor. 

  // Use cell-average value for gamma. 
  double gamma_avg = gamma[0]/sqrt(pow(2, 1)); 
  double dv1_sq = 4.0/dxv[3]/dxv[3]; 

  const double* G_L = fpo_g_stencil[0]; 
  const double* G_C = fpo_g_stencil[1]; 
  const double* G_R = fpo_g_stencil[2]; 
  
  const double* d2G_surf_C = &fpo_d2gdv2_surf[40]; 
  
  double *out = &diff_coeff[384]; 
  
  out[0] = 0.6708203932499369*(G_R[14]*dv1_sq*gamma_avg+G_L[14]*dv1_sq*gamma_avg)+0.9375*(G_R[0]*dv1_sq*gamma_avg+G_L[0]*dv1_sq*gamma_avg)-1.3416407864998738*G_C[14]*dv1_sq*gamma_avg-1.190784930203603*G_R[4]*dv1_sq*gamma_avg+1.190784930203603*G_L[4]*dv1_sq*gamma_avg-1.875*G_C[0]*dv1_sq*gamma_avg; 
  out[1] = 0.6708203932499369*(G_R[28]*dv1_sq*gamma_avg+G_L[28]*dv1_sq*gamma_avg)+0.9375*(G_R[1]*dv1_sq*gamma_avg+G_L[1]*dv1_sq*gamma_avg)-1.3416407864998738*G_C[28]*dv1_sq*gamma_avg-1.190784930203603*G_R[8]*dv1_sq*gamma_avg+1.190784930203603*G_L[8]*dv1_sq*gamma_avg-1.875*G_C[1]*dv1_sq*gamma_avg; 
  out[2] = 0.6708203932499369*(G_R[29]*dv1_sq*gamma_avg+G_L[29]*dv1_sq*gamma_avg)+0.9375*(G_R[2]*dv1_sq*gamma_avg+G_L[2]*dv1_sq*gamma_avg)-1.3416407864998738*G_C[29]*dv1_sq*gamma_avg-1.190784930203603*G_R[9]*dv1_sq*gamma_avg+1.190784930203603*G_L[9]*dv1_sq*gamma_avg-1.875*G_C[2]*dv1_sq*gamma_avg; 
  out[3] = 0.6708203932499369*(G_R[30]*dv1_sq*gamma_avg+G_L[30]*dv1_sq*gamma_avg)+0.9375*(G_R[3]*dv1_sq*gamma_avg+G_L[3]*dv1_sq*gamma_avg)-1.3416407864998738*G_C[30]*dv1_sq*gamma_avg-1.190784930203603*G_R[10]*dv1_sq*gamma_avg+1.190784930203603*G_L[10]*dv1_sq*gamma_avg-1.875*G_C[3]*dv1_sq*gamma_avg; 
  out[4] = -(1.453125*(G_R[4]*dv1_sq*gamma_avg+G_L[4]*dv1_sq*gamma_avg))+0.7382874503707888*G_R[14]*dv1_sq*gamma_avg-0.7382874503707888*G_L[14]*dv1_sq*gamma_avg-5.34375*G_C[4]*dv1_sq*gamma_avg+1.190784930203603*G_R[0]*dv1_sq*gamma_avg-1.190784930203603*G_L[0]*dv1_sq*gamma_avg; 
  out[5] = 0.6708203932499369*(G_R[41]*dv1_sq*gamma_avg+G_L[41]*dv1_sq*gamma_avg)+0.9375*(G_R[5]*dv1_sq*gamma_avg+G_L[5]*dv1_sq*gamma_avg)-1.3416407864998738*G_C[41]*dv1_sq*gamma_avg-1.190784930203603*G_R[16]*dv1_sq*gamma_avg+1.190784930203603*G_L[16]*dv1_sq*gamma_avg-1.875*G_C[5]*dv1_sq*gamma_avg; 
  out[6] = 0.6708203932499369*(G_R[42]*dv1_sq*gamma_avg+G_L[42]*dv1_sq*gamma_avg)+0.9375*(G_R[6]*dv1_sq*gamma_avg+G_L[6]*dv1_sq*gamma_avg)-1.3416407864998738*G_C[42]*dv1_sq*gamma_avg-1.190784930203603*G_R[17]*dv1_sq*gamma_avg+1.190784930203603*G_L[17]*dv1_sq*gamma_avg-1.875*G_C[6]*dv1_sq*gamma_avg; 
  out[7] = 0.6708203932499369*(G_R[43]*dv1_sq*gamma_avg+G_L[43]*dv1_sq*gamma_avg)+0.9375*(G_R[7]*dv1_sq*gamma_avg+G_L[7]*dv1_sq*gamma_avg)-1.3416407864998738*G_C[43]*dv1_sq*gamma_avg-1.190784930203603*G_R[18]*dv1_sq*gamma_avg+1.190784930203603*G_L[18]*dv1_sq*gamma_avg-1.875*G_C[7]*dv1_sq*gamma_avg; 
  out[8] = -(1.453125*(G_R[8]*dv1_sq*gamma_avg+G_L[8]*dv1_sq*gamma_avg))+0.7382874503707888*G_R[28]*dv1_sq*gamma_avg-0.7382874503707888*G_L[28]*dv1_sq*gamma_avg-5.34375*G_C[8]*dv1_sq*gamma_avg+1.190784930203603*G_R[1]*dv1_sq*gamma_avg-1.190784930203603*G_L[1]*dv1_sq*gamma_avg; 
  out[9] = -(1.453125*(G_R[9]*dv1_sq*gamma_avg+G_L[9]*dv1_sq*gamma_avg))+0.7382874503707888*G_R[29]*dv1_sq*gamma_avg-0.7382874503707888*G_L[29]*dv1_sq*gamma_avg-5.34375*G_C[9]*dv1_sq*gamma_avg+1.190784930203603*G_R[2]*dv1_sq*gamma_avg-1.190784930203603*G_L[2]*dv1_sq*gamma_avg; 
  out[10] = -(1.453125*(G_R[10]*dv1_sq*gamma_avg+G_L[10]*dv1_sq*gamma_avg))+0.7382874503707888*G_R[30]*dv1_sq*gamma_avg-0.7382874503707888*G_L[30]*dv1_sq*gamma_avg-5.34375*G_C[10]*dv1_sq*gamma_avg+1.190784930203603*G_R[3]*dv1_sq*gamma_avg-1.190784930203603*G_L[3]*dv1_sq*gamma_avg; 
  out[11] = 0.9375*(G_R[11]*dv1_sq*gamma_avg+G_L[11]*dv1_sq*gamma_avg)-1.190784930203603*G_R[25]*dv1_sq*gamma_avg+1.190784930203603*G_L[25]*dv1_sq*gamma_avg-1.875*G_C[11]*dv1_sq*gamma_avg; 
  out[12] = 0.9375*(G_R[12]*dv1_sq*gamma_avg+G_L[12]*dv1_sq*gamma_avg)-1.190784930203603*G_R[26]*dv1_sq*gamma_avg+1.190784930203603*G_L[26]*dv1_sq*gamma_avg-1.875*G_C[12]*dv1_sq*gamma_avg; 
  out[13] = 0.9375*(G_R[13]*dv1_sq*gamma_avg+G_L[13]*dv1_sq*gamma_avg)-1.190784930203603*G_R[27]*dv1_sq*gamma_avg+1.190784930203603*G_L[27]*dv1_sq*gamma_avg-1.875*G_C[13]*dv1_sq*gamma_avg; 
  out[14] = -(0.140625*(G_R[14]*dv1_sq*gamma_avg+G_L[14]*dv1_sq*gamma_avg))+0.4192627457812106*(G_R[0]*dv1_sq*gamma_avg+G_L[0]*dv1_sq*gamma_avg)-6.28125*G_C[14]*dv1_sq*gamma_avg-0.3025768239224545*G_R[4]*dv1_sq*gamma_avg+0.3025768239224545*G_L[4]*dv1_sq*gamma_avg-0.8385254915624212*G_C[0]*dv1_sq*gamma_avg; 
  out[15] = 0.6708203932499369*(G_R[47]*dv1_sq*gamma_avg+G_L[47]*dv1_sq*gamma_avg)+0.9375*(G_R[15]*dv1_sq*gamma_avg+G_L[15]*dv1_sq*gamma_avg)-1.3416407864998738*G_C[47]*dv1_sq*gamma_avg-1.190784930203603*G_R[31]*dv1_sq*gamma_avg+1.190784930203603*G_L[31]*dv1_sq*gamma_avg-1.875*G_C[15]*dv1_sq*gamma_avg; 
  out[16] = -(1.453125*(G_R[16]*dv1_sq*gamma_avg+G_L[16]*dv1_sq*gamma_avg))+0.7382874503707888*G_R[41]*dv1_sq*gamma_avg-0.7382874503707888*G_L[41]*dv1_sq*gamma_avg-5.34375*G_C[16]*dv1_sq*gamma_avg+1.190784930203603*G_R[5]*dv1_sq*gamma_avg-1.190784930203603*G_L[5]*dv1_sq*gamma_avg; 
  out[17] = -(1.453125*(G_R[17]*dv1_sq*gamma_avg+G_L[17]*dv1_sq*gamma_avg))+0.7382874503707888*G_R[42]*dv1_sq*gamma_avg-0.7382874503707888*G_L[42]*dv1_sq*gamma_avg-5.34375*G_C[17]*dv1_sq*gamma_avg+1.190784930203603*G_R[6]*dv1_sq*gamma_avg-1.190784930203603*G_L[6]*dv1_sq*gamma_avg; 
  out[18] = -(1.453125*(G_R[18]*dv1_sq*gamma_avg+G_L[18]*dv1_sq*gamma_avg))+0.7382874503707888*G_R[43]*dv1_sq*gamma_avg-0.7382874503707888*G_L[43]*dv1_sq*gamma_avg-5.34375*G_C[18]*dv1_sq*gamma_avg+1.190784930203603*G_R[7]*dv1_sq*gamma_avg-1.190784930203603*G_L[7]*dv1_sq*gamma_avg; 
  out[19] = 0.9375*(G_R[19]*dv1_sq*gamma_avg+G_L[19]*dv1_sq*gamma_avg)-1.190784930203603*G_R[35]*dv1_sq*gamma_avg+1.190784930203603*G_L[35]*dv1_sq*gamma_avg-1.875*G_C[19]*dv1_sq*gamma_avg; 
  out[20] = 0.9375*(G_R[20]*dv1_sq*gamma_avg+G_L[20]*dv1_sq*gamma_avg)-1.190784930203603*G_R[36]*dv1_sq*gamma_avg+1.190784930203603*G_L[36]*dv1_sq*gamma_avg-1.875*G_C[20]*dv1_sq*gamma_avg; 
  out[21] = 0.9375*(G_R[21]*dv1_sq*gamma_avg+G_L[21]*dv1_sq*gamma_avg)-1.190784930203603*G_R[37]*dv1_sq*gamma_avg+1.190784930203603*G_L[37]*dv1_sq*gamma_avg-1.875*G_C[21]*dv1_sq*gamma_avg; 
  out[22] = 0.9375*(G_R[22]*dv1_sq*gamma_avg+G_L[22]*dv1_sq*gamma_avg)-1.190784930203603*G_R[38]*dv1_sq*gamma_avg+1.190784930203603*G_L[38]*dv1_sq*gamma_avg-1.875*G_C[22]*dv1_sq*gamma_avg; 
  out[23] = 0.9375*(G_R[23]*dv1_sq*gamma_avg+G_L[23]*dv1_sq*gamma_avg)-1.190784930203603*G_R[39]*dv1_sq*gamma_avg+1.190784930203603*G_L[39]*dv1_sq*gamma_avg-1.875*G_C[23]*dv1_sq*gamma_avg; 
  out[24] = 0.9375*(G_R[24]*dv1_sq*gamma_avg+G_L[24]*dv1_sq*gamma_avg)-1.190784930203603*G_R[40]*dv1_sq*gamma_avg+1.190784930203603*G_L[40]*dv1_sq*gamma_avg-1.875*G_C[24]*dv1_sq*gamma_avg; 
  out[25] = -(1.453125*(G_R[25]*dv1_sq*gamma_avg+G_L[25]*dv1_sq*gamma_avg))-5.34375*G_C[25]*dv1_sq*gamma_avg+1.190784930203603*G_R[11]*dv1_sq*gamma_avg-1.190784930203603*G_L[11]*dv1_sq*gamma_avg; 
  out[26] = -(1.453125*(G_R[26]*dv1_sq*gamma_avg+G_L[26]*dv1_sq*gamma_avg))-5.34375*G_C[26]*dv1_sq*gamma_avg+1.190784930203603*G_R[12]*dv1_sq*gamma_avg-1.190784930203603*G_L[12]*dv1_sq*gamma_avg; 
  out[27] = -(1.453125*(G_R[27]*dv1_sq*gamma_avg+G_L[27]*dv1_sq*gamma_avg))-5.34375*G_C[27]*dv1_sq*gamma_avg+1.190784930203603*G_R[13]*dv1_sq*gamma_avg-1.190784930203603*G_L[13]*dv1_sq*gamma_avg; 
  out[28] = -(0.140625*(G_R[28]*dv1_sq*gamma_avg+G_L[28]*dv1_sq*gamma_avg))+0.41926274578121053*(G_R[1]*dv1_sq*gamma_avg+G_L[1]*dv1_sq*gamma_avg)-6.28125*G_C[28]*dv1_sq*gamma_avg-0.30257682392245444*G_R[8]*dv1_sq*gamma_avg+0.30257682392245444*G_L[8]*dv1_sq*gamma_avg-0.8385254915624211*G_C[1]*dv1_sq*gamma_avg; 
  out[29] = -(0.140625*(G_R[29]*dv1_sq*gamma_avg+G_L[29]*dv1_sq*gamma_avg))+0.41926274578121053*(G_R[2]*dv1_sq*gamma_avg+G_L[2]*dv1_sq*gamma_avg)-6.28125*G_C[29]*dv1_sq*gamma_avg-0.30257682392245444*G_R[9]*dv1_sq*gamma_avg+0.30257682392245444*G_L[9]*dv1_sq*gamma_avg-0.8385254915624211*G_C[2]*dv1_sq*gamma_avg; 
  out[30] = -(0.140625*(G_R[30]*dv1_sq*gamma_avg+G_L[30]*dv1_sq*gamma_avg))+0.41926274578121053*(G_R[3]*dv1_sq*gamma_avg+G_L[3]*dv1_sq*gamma_avg)-6.28125*G_C[30]*dv1_sq*gamma_avg-0.30257682392245444*G_R[10]*dv1_sq*gamma_avg+0.30257682392245444*G_L[10]*dv1_sq*gamma_avg-0.8385254915624211*G_C[3]*dv1_sq*gamma_avg; 
  out[31] = -(1.453125*(G_R[31]*dv1_sq*gamma_avg+G_L[31]*dv1_sq*gamma_avg))+0.7382874503707888*G_R[47]*dv1_sq*gamma_avg-0.7382874503707888*G_L[47]*dv1_sq*gamma_avg-5.34375*G_C[31]*dv1_sq*gamma_avg+1.190784930203603*G_R[15]*dv1_sq*gamma_avg-1.190784930203603*G_L[15]*dv1_sq*gamma_avg; 
  out[32] = 0.9375*(G_R[32]*dv1_sq*gamma_avg+G_L[32]*dv1_sq*gamma_avg)-1.190784930203603*G_R[44]*dv1_sq*gamma_avg+1.190784930203603*G_L[44]*dv1_sq*gamma_avg-1.875*G_C[32]*dv1_sq*gamma_avg; 
  out[33] = 0.9375*(G_R[33]*dv1_sq*gamma_avg+G_L[33]*dv1_sq*gamma_avg)-1.190784930203603*G_R[45]*dv1_sq*gamma_avg+1.190784930203603*G_L[45]*dv1_sq*gamma_avg-1.875*G_C[33]*dv1_sq*gamma_avg; 
  out[34] = 0.9375*(G_R[34]*dv1_sq*gamma_avg+G_L[34]*dv1_sq*gamma_avg)-1.190784930203603*G_R[46]*dv1_sq*gamma_avg+1.190784930203603*G_L[46]*dv1_sq*gamma_avg-1.875*G_C[34]*dv1_sq*gamma_avg; 
  out[35] = -(1.453125*(G_R[35]*dv1_sq*gamma_avg+G_L[35]*dv1_sq*gamma_avg))-5.34375*G_C[35]*dv1_sq*gamma_avg+1.190784930203603*G_R[19]*dv1_sq*gamma_avg-1.190784930203603*G_L[19]*dv1_sq*gamma_avg; 
  out[36] = -(1.453125*(G_R[36]*dv1_sq*gamma_avg+G_L[36]*dv1_sq*gamma_avg))-5.34375*G_C[36]*dv1_sq*gamma_avg+1.190784930203603*G_R[20]*dv1_sq*gamma_avg-1.190784930203603*G_L[20]*dv1_sq*gamma_avg; 
  out[37] = -(1.453125*(G_R[37]*dv1_sq*gamma_avg+G_L[37]*dv1_sq*gamma_avg))-5.34375*G_C[37]*dv1_sq*gamma_avg+1.190784930203603*G_R[21]*dv1_sq*gamma_avg-1.190784930203603*G_L[21]*dv1_sq*gamma_avg; 
  out[38] = -(1.453125*(G_R[38]*dv1_sq*gamma_avg+G_L[38]*dv1_sq*gamma_avg))-5.34375*G_C[38]*dv1_sq*gamma_avg+1.190784930203603*G_R[22]*dv1_sq*gamma_avg-1.190784930203603*G_L[22]*dv1_sq*gamma_avg; 
  out[39] = -(1.453125*(G_R[39]*dv1_sq*gamma_avg+G_L[39]*dv1_sq*gamma_avg))-5.34375*G_C[39]*dv1_sq*gamma_avg+1.190784930203603*G_R[23]*dv1_sq*gamma_avg-1.190784930203603*G_L[23]*dv1_sq*gamma_avg; 
  out[40] = -(1.453125*(G_R[40]*dv1_sq*gamma_avg+G_L[40]*dv1_sq*gamma_avg))-5.34375*G_C[40]*dv1_sq*gamma_avg+1.190784930203603*G_R[24]*dv1_sq*gamma_avg-1.190784930203603*G_L[24]*dv1_sq*gamma_avg; 
  out[41] = -(0.140625*(G_R[41]*dv1_sq*gamma_avg+G_L[41]*dv1_sq*gamma_avg))+0.4192627457812106*(G_R[5]*dv1_sq*gamma_avg+G_L[5]*dv1_sq*gamma_avg)-6.28125*G_C[41]*dv1_sq*gamma_avg-0.3025768239224545*G_R[16]*dv1_sq*gamma_avg+0.3025768239224545*G_L[16]*dv1_sq*gamma_avg-0.8385254915624212*G_C[5]*dv1_sq*gamma_avg; 
  out[42] = -(0.140625*(G_R[42]*dv1_sq*gamma_avg+G_L[42]*dv1_sq*gamma_avg))+0.4192627457812106*(G_R[6]*dv1_sq*gamma_avg+G_L[6]*dv1_sq*gamma_avg)-6.28125*G_C[42]*dv1_sq*gamma_avg-0.3025768239224545*G_R[17]*dv1_sq*gamma_avg+0.3025768239224545*G_L[17]*dv1_sq*gamma_avg-0.8385254915624212*G_C[6]*dv1_sq*gamma_avg; 
  out[43] = -(0.140625*(G_R[43]*dv1_sq*gamma_avg+G_L[43]*dv1_sq*gamma_avg))+0.4192627457812106*(G_R[7]*dv1_sq*gamma_avg+G_L[7]*dv1_sq*gamma_avg)-6.28125*G_C[43]*dv1_sq*gamma_avg-0.3025768239224545*G_R[18]*dv1_sq*gamma_avg+0.3025768239224545*G_L[18]*dv1_sq*gamma_avg-0.8385254915624212*G_C[7]*dv1_sq*gamma_avg; 
  out[44] = -(1.453125*(G_R[44]*dv1_sq*gamma_avg+G_L[44]*dv1_sq*gamma_avg))-5.34375*G_C[44]*dv1_sq*gamma_avg+1.190784930203603*G_R[32]*dv1_sq*gamma_avg-1.190784930203603*G_L[32]*dv1_sq*gamma_avg; 
  out[45] = -(1.453125*(G_R[45]*dv1_sq*gamma_avg+G_L[45]*dv1_sq*gamma_avg))-5.34375*G_C[45]*dv1_sq*gamma_avg+1.190784930203603*G_R[33]*dv1_sq*gamma_avg-1.190784930203603*G_L[33]*dv1_sq*gamma_avg; 
  out[46] = -(1.453125*(G_R[46]*dv1_sq*gamma_avg+G_L[46]*dv1_sq*gamma_avg))-5.34375*G_C[46]*dv1_sq*gamma_avg+1.190784930203603*G_R[34]*dv1_sq*gamma_avg-1.190784930203603*G_L[34]*dv1_sq*gamma_avg; 
  out[47] = -(0.140625*(G_R[47]*dv1_sq*gamma_avg+G_L[47]*dv1_sq*gamma_avg))+0.41926274578121053*(G_R[15]*dv1_sq*gamma_avg+G_L[15]*dv1_sq*gamma_avg)-6.28125*G_C[47]*dv1_sq*gamma_avg-0.30257682392245444*G_R[31]*dv1_sq*gamma_avg+0.30257682392245444*G_L[31]*dv1_sq*gamma_avg-0.8385254915624211*G_C[15]*dv1_sq*gamma_avg; 
} 

