#include <gkyl_fpo_vlasov_kernels.h> 
#include <math.h> 
 
void fpo_diff_coeff_cross_2x3v_vyvz_ser_p1_lovy_invz(const double *dxv, const double *gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff) { 
  // dxv[NDIM]: Cell spacing in each direction. 
  // gamma: Scalar factor gamma. 
  // fpo_g_stencil[9]: 9 cell stencil of Rosenbluth potential G. 
  // fpo_g_surf_stencil[9]: 9 cell stencil of surface projection of G. 
  // fpo_dgdv_surf: Surface expansion of dG/dv in center cell. 
  // diff_coeff: Output array for diffusion tensor. 

  // Use cell-average value for gamma. 
  double gamma_avg = gamma[0]/sqrt(pow(2, 2)); 
  double dv1 = 2.0/dxv[3]; 
  double dv2 = 2.0/dxv[4]; 
  double dv1_sq = 4.0/dxv[3]/dxv[4]; 
 
  double surft1_lo[16] = {0.0}; 
  double surft1_up[16] = {0.0}; 
  double surft2_lo[16] = {0.0}; 
  double surft2_up[16] = {0.0}; 
  double vol[32] = {0.0}; 
  double *out = &diff_coeff[160]; 

  const double* GBL = fpo_g_stencil[0]; 
  const double* GCL = fpo_g_stencil[1]; 
  const double* GTL = fpo_g_stencil[2]; 
  const double* GBC = fpo_g_stencil[3]; 
  const double* GCC = fpo_g_stencil[4]; 
  const double* G_surf_CC_vy = &fpo_g_surf_stencil[4][16]; 
  const double* GTC = fpo_g_stencil[5]; 
  const double* GBR = fpo_g_stencil[6]; 
  const double* GCR = fpo_g_stencil[7]; 
  const double* GTR = fpo_g_stencil[8]; 
  const double* dGdvy_surf_CC_vz = &fpo_dgdv_surf[80]; 

  surft1_lo[0] = -(0.7071067811865475*GCC[15])+0.7071067811865475*GBC[15]+0.6123724356957944*(GCC[4]+GBC[4]); 
  surft1_lo[1] = -(0.7071067811865475*GCC[23])+0.7071067811865475*GBC[23]+0.6123724356957944*(GCC[9]+GBC[9]); 
  surft1_lo[2] = -(0.7071067811865475*GCC[24])+0.7071067811865475*GBC[24]+0.6123724356957944*(GCC[10]+GBC[10]); 
  surft1_lo[3] = -(0.7071067811865475*GCC[25])+0.7071067811865475*GBC[25]+0.6123724356957944*(GCC[11]+GBC[11]); 
  surft1_lo[5] = -(0.7071067811865475*GCC[28])+0.7071067811865475*GBC[28]+0.6123724356957944*(GCC[17]+GBC[17]); 
  surft1_lo[6] = -(0.7071067811865475*GCC[29])+0.7071067811865475*GBC[29]+0.6123724356957944*(GCC[18]+GBC[18]); 
  surft1_lo[7] = -(0.7071067811865475*GCC[30])+0.7071067811865475*GBC[30]+0.6123724356957944*(GCC[19]+GBC[19]); 
  surft1_lo[11] = -(0.7071067811865475*GCC[31])+0.7071067811865475*GBC[31]+0.6123724356957944*(GCC[26]+GBC[26]); 
  surft1_up[0] = -(0.7071067811865475*GTC[15])+0.7071067811865475*GCC[15]+0.6123724356957944*(GTC[4]+GCC[4]); 
  surft1_up[1] = -(0.7071067811865475*GTC[23])+0.7071067811865475*GCC[23]+0.6123724356957944*(GTC[9]+GCC[9]); 
  surft1_up[2] = -(0.7071067811865475*GTC[24])+0.7071067811865475*GCC[24]+0.6123724356957944*(GTC[10]+GCC[10]); 
  surft1_up[3] = -(0.7071067811865475*GTC[25])+0.7071067811865475*GCC[25]+0.6123724356957944*(GTC[11]+GCC[11]); 
  surft1_up[5] = -(0.7071067811865475*GTC[28])+0.7071067811865475*GCC[28]+0.6123724356957944*(GTC[17]+GCC[17]); 
  surft1_up[6] = -(0.7071067811865475*GTC[29])+0.7071067811865475*GCC[29]+0.6123724356957944*(GTC[18]+GCC[18]); 
  surft1_up[7] = -(0.7071067811865475*GTC[30])+0.7071067811865475*GCC[30]+0.6123724356957944*(GTC[19]+GCC[19]); 
  surft1_up[11] = -(0.7071067811865475*GTC[31])+0.7071067811865475*GCC[31]+0.6123724356957944*(GTC[26]+GCC[26]); 

  surft2_lo[0] = G_surf_CC_vy[0]; 
  surft2_lo[1] = G_surf_CC_vy[1]; 
  surft2_lo[2] = G_surf_CC_vy[2]; 
  surft2_lo[3] = G_surf_CC_vy[3]; 
  surft2_lo[4] = G_surf_CC_vy[4]; 
  surft2_lo[5] = G_surf_CC_vy[5]; 
  surft2_lo[6] = G_surf_CC_vy[6]; 
  surft2_lo[7] = G_surf_CC_vy[7]; 
  surft2_lo[8] = G_surf_CC_vy[8]; 
  surft2_lo[9] = G_surf_CC_vy[9]; 
  surft2_lo[10] = G_surf_CC_vy[10]; 
  surft2_lo[11] = G_surf_CC_vy[11]; 
  surft2_lo[12] = G_surf_CC_vy[12]; 
  surft2_lo[13] = G_surf_CC_vy[13]; 
  surft2_lo[14] = G_surf_CC_vy[14]; 
  surft2_lo[15] = G_surf_CC_vy[15]; 
  surft2_up[0] = -(0.408248290463863*GCR[4])+0.408248290463863*GCC[4]+0.3535533905932737*(GCR[0]+GCC[0]); 
  surft2_up[1] = -(0.408248290463863*GCR[9])+0.408248290463863*GCC[9]+0.3535533905932737*(GCR[1]+GCC[1]); 
  surft2_up[2] = -(0.408248290463863*GCR[10])+0.408248290463863*GCC[10]+0.3535533905932737*(GCR[2]+GCC[2]); 
  surft2_up[3] = -(0.408248290463863*GCR[11])+0.408248290463863*GCC[11]+0.3535533905932737*(GCR[3]+GCC[3]); 
  surft2_up[4] = -(0.408248290463863*GCR[15])+0.408248290463863*GCC[15]+0.3535533905932737*(GCR[5]+GCC[5]); 
  surft2_up[5] = -(0.408248290463863*GCR[17])+0.408248290463863*GCC[17]+0.3535533905932737*(GCR[6]+GCC[6]); 
  surft2_up[6] = -(0.408248290463863*GCR[18])+0.408248290463863*GCC[18]+0.3535533905932737*(GCR[7]+GCC[7]); 
  surft2_up[7] = -(0.408248290463863*GCR[19])+0.408248290463863*GCC[19]+0.3535533905932737*(GCR[8]+GCC[8]); 
  surft2_up[8] = -(0.408248290463863*GCR[23])+0.408248290463863*GCC[23]+0.3535533905932737*(GCR[12]+GCC[12]); 
  surft2_up[9] = -(0.408248290463863*GCR[24])+0.408248290463863*GCC[24]+0.3535533905932737*(GCR[13]+GCC[13]); 
  surft2_up[10] = -(0.408248290463863*GCR[25])+0.408248290463863*GCC[25]+0.3535533905932737*(GCR[14]+GCC[14]); 
  surft2_up[11] = -(0.408248290463863*GCR[26])+0.408248290463863*GCC[26]+0.3535533905932737*(GCR[16]+GCC[16]); 
  surft2_up[12] = -(0.408248290463863*GCR[28])+0.408248290463863*GCC[28]+0.3535533905932737*(GCR[20]+GCC[20]); 
  surft2_up[13] = -(0.408248290463863*GCR[29])+0.408248290463863*GCC[29]+0.3535533905932737*(GCR[21]+GCC[21]); 
  surft2_up[14] = -(0.408248290463863*GCR[30])+0.408248290463863*GCC[30]+0.3535533905932737*(GCR[22]+GCC[22]); 
  surft2_up[15] = -(0.408248290463863*GCR[31])+0.408248290463863*GCC[31]+0.3535533905932737*(GCR[27]+GCC[27]); 

  vol[15] = 3.0*GCC[0]; 
  vol[23] = 3.0*GCC[1]; 
  vol[24] = 3.0*GCC[2]; 
  vol[25] = 3.0*GCC[3]; 
  vol[28] = 3.0*GCC[6]; 
  vol[29] = 3.0*GCC[7]; 
  vol[30] = 3.0*GCC[8]; 
  vol[31] = 3.0*GCC[16]; 

  out[0] = (vol[0]+0.7071067811865475*surft1_up[0]-0.7071067811865475*surft1_lo[0])*dv1_sq*gamma_avg; 
  out[1] = (vol[1]+0.7071067811865475*surft1_up[1]-0.7071067811865475*surft1_lo[1])*dv1_sq*gamma_avg; 
  out[2] = (vol[2]+0.7071067811865475*surft1_up[2]-0.7071067811865475*surft1_lo[2])*dv1_sq*gamma_avg; 
  out[3] = (vol[3]+0.7071067811865475*surft1_up[3]-0.7071067811865475*surft1_lo[3])*dv1_sq*gamma_avg; 
  out[4] = (vol[4]+0.7071067811865475*surft1_up[4]-0.7071067811865475*surft1_lo[4])*dv1_sq*gamma_avg; 
  out[5] = (vol[5]-1.224744871391589*surft2_up[0]+1.224744871391589*(surft2_lo[0]+surft1_up[0]+surft1_lo[0]))*dv1_sq*gamma_avg; 
  out[6] = (vol[6]+0.7071067811865475*surft1_up[5]-0.7071067811865475*surft1_lo[5])*dv1_sq*gamma_avg; 
  out[7] = (vol[7]+0.7071067811865475*surft1_up[6]-0.7071067811865475*surft1_lo[6])*dv1_sq*gamma_avg; 
  out[8] = (vol[8]+0.7071067811865475*surft1_up[7]-0.7071067811865475*surft1_lo[7])*dv1_sq*gamma_avg; 
  out[9] = (vol[9]+0.7071067811865475*surft1_up[8]-0.7071067811865475*surft1_lo[8])*dv1_sq*gamma_avg; 
  out[10] = (vol[10]+0.7071067811865475*surft1_up[9]-0.7071067811865475*surft1_lo[9])*dv1_sq*gamma_avg; 
  out[11] = (vol[11]+0.7071067811865475*surft1_up[10]-0.7071067811865475*surft1_lo[10])*dv1_sq*gamma_avg; 
  out[12] = (vol[12]-1.224744871391589*surft2_up[1]+1.224744871391589*(surft2_lo[1]+surft1_up[1]+surft1_lo[1]))*dv1_sq*gamma_avg; 
  out[13] = (vol[13]-1.224744871391589*surft2_up[2]+1.224744871391589*(surft2_lo[2]+surft1_up[2]+surft1_lo[2]))*dv1_sq*gamma_avg; 
  out[14] = (vol[14]-1.224744871391589*surft2_up[3]+1.224744871391589*(surft2_lo[3]+surft1_up[3]+surft1_lo[3]))*dv1_sq*gamma_avg; 
  out[15] = (vol[15]+1.224744871391589*(surft1_up[4]+surft1_lo[4])-2.1213203435596424*(surft2_up[0]+surft2_lo[0]))*dv1_sq*gamma_avg; 
  out[16] = (vol[16]+0.7071067811865475*surft1_up[11]-0.7071067811865475*surft1_lo[11])*dv1_sq*gamma_avg; 
  out[17] = (vol[17]+0.7071067811865475*surft1_up[12]-0.7071067811865475*surft1_lo[12])*dv1_sq*gamma_avg; 
  out[18] = (vol[18]+0.7071067811865475*surft1_up[13]-0.7071067811865475*surft1_lo[13])*dv1_sq*gamma_avg; 
  out[19] = (vol[19]+0.7071067811865475*surft1_up[14]-0.7071067811865475*surft1_lo[14])*dv1_sq*gamma_avg; 
  out[20] = (vol[20]-1.224744871391589*surft2_up[5]+1.224744871391589*(surft2_lo[5]+surft1_up[5]+surft1_lo[5]))*dv1_sq*gamma_avg; 
  out[21] = (vol[21]-1.224744871391589*surft2_up[6]+1.224744871391589*(surft2_lo[6]+surft1_up[6]+surft1_lo[6]))*dv1_sq*gamma_avg; 
  out[22] = (vol[22]-1.224744871391589*surft2_up[7]+1.224744871391589*(surft2_lo[7]+surft1_up[7]+surft1_lo[7]))*dv1_sq*gamma_avg; 
  out[23] = (vol[23]+1.224744871391589*(surft1_up[8]+surft1_lo[8])-2.1213203435596424*(surft2_up[1]+surft2_lo[1]))*dv1_sq*gamma_avg; 
  out[24] = (vol[24]+1.224744871391589*(surft1_up[9]+surft1_lo[9])-2.1213203435596424*(surft2_up[2]+surft2_lo[2]))*dv1_sq*gamma_avg; 
  out[25] = (vol[25]+1.224744871391589*(surft1_up[10]+surft1_lo[10])-2.1213203435596424*(surft2_up[3]+surft2_lo[3]))*dv1_sq*gamma_avg; 
  out[26] = (vol[26]+0.7071067811865475*surft1_up[15]-0.7071067811865475*surft1_lo[15])*dv1_sq*gamma_avg; 
  out[27] = (vol[27]-1.224744871391589*surft2_up[11]+1.224744871391589*(surft2_lo[11]+surft1_up[11]+surft1_lo[11]))*dv1_sq*gamma_avg; 
  out[28] = (vol[28]+1.224744871391589*(surft1_up[12]+surft1_lo[12])-2.1213203435596424*(surft2_up[5]+surft2_lo[5]))*dv1_sq*gamma_avg; 
  out[29] = (vol[29]+1.224744871391589*(surft1_up[13]+surft1_lo[13])-2.1213203435596424*(surft2_up[6]+surft2_lo[6]))*dv1_sq*gamma_avg; 
  out[30] = (vol[30]+1.224744871391589*(surft1_up[14]+surft1_lo[14])-2.1213203435596424*(surft2_up[7]+surft2_lo[7]))*dv1_sq*gamma_avg; 
  out[31] = (vol[31]+1.224744871391589*(surft1_up[15]+surft1_lo[15])-2.1213203435596424*(surft2_up[11]+surft2_lo[11]))*dv1_sq*gamma_avg; 
} 
