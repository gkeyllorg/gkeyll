#include <gkyl_fpo_vlasov_kernels.h> 
 
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vyvx_ser_p1_upvy_upvx(const double *dxv, const double *gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff) { 
  // dxv[NDIM]: Cell spacing in each direction. 
  // gamma: Scalar factor gamma. 
  // fpo_g_stencil[9]: 9 cell stencil of Rosenbluth potential G. 
  // fpo_g_surf_stencil[9]: 9 cell stencil of surface projection of G. 
  // fpo_dgdv_surf: Surface expansion of dG/dv in center cell. 
  // diff_coeff: Output array for diffusion tensor. 

  // Use cell-average value for gamma. 
  double gamma_avg = gamma[0]/sqrt(pow(2, 2)); 
  double dv1 = 2.0/dxv[3]; 
  double dv2 = 2.0/dxv[2]; 
  double dv1_sq = 4.0/dxv[3]/dxv[2]; 
 
  double surft1_lo[16] = {0.0}; 
  double surft1_up[16] = {0.0}; 
  double surft2_lo[16] = {0.0}; 
  double surft2_up[16] = {0.0}; 
  double vol[32] = {0.0}; 
  double *out = &diff_coeff[96]; 

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
  const double* dGdvy_surf_CC_vx = &fpo_dgdv_surf[0]; 

  surft1_lo[0] = 0.7071067811865475*GCL[11]-0.7071067811865475*GCC[11]+0.6123724356957944*(GCL[4]+GCC[4]); 
  surft1_lo[1] = 0.7071067811865475*GCL[18]-0.7071067811865475*GCC[18]+0.6123724356957944*(GCL[9]+GCC[9]); 
  surft1_lo[2] = 0.7071067811865475*GCL[19]-0.7071067811865475*GCC[19]+0.6123724356957944*(GCL[10]+GCC[10]); 
  surft1_lo[4] = 0.7071067811865475*GCL[25]-0.7071067811865475*GCC[25]+0.6123724356957944*(GCL[15]+GCC[15]); 
  surft1_lo[5] = 0.7071067811865475*GCL[26]-0.7071067811865475*GCC[26]+0.6123724356957944*(GCL[17]+GCC[17]); 
  surft1_lo[8] = 0.7071067811865475*GCL[29]-0.7071067811865475*GCC[29]+0.6123724356957944*(GCL[23]+GCC[23]); 
  surft1_lo[9] = 0.7071067811865475*GCL[30]-0.7071067811865475*GCC[30]+0.6123724356957944*(GCL[24]+GCC[24]); 
  surft1_lo[12] = 0.7071067811865475*GCL[31]-0.7071067811865475*GCC[31]+0.6123724356957944*(GCL[28]+GCC[28]); 
  surft1_up[0] = dGdvy_surf_CC_vx[0]/dv1; 
  surft1_up[1] = dGdvy_surf_CC_vx[1]/dv1; 
  surft1_up[2] = dGdvy_surf_CC_vx[2]/dv1; 
  surft1_up[3] = dGdvy_surf_CC_vx[3]/dv1; 
  surft1_up[4] = dGdvy_surf_CC_vx[4]/dv1; 
  surft1_up[5] = dGdvy_surf_CC_vx[5]/dv1; 
  surft1_up[6] = dGdvy_surf_CC_vx[6]/dv1; 
  surft1_up[7] = dGdvy_surf_CC_vx[7]/dv1; 
  surft1_up[8] = dGdvy_surf_CC_vx[8]/dv1; 
  surft1_up[9] = dGdvy_surf_CC_vx[9]/dv1; 
  surft1_up[10] = dGdvy_surf_CC_vx[10]/dv1; 
  surft1_up[11] = dGdvy_surf_CC_vx[11]/dv1; 
  surft1_up[12] = dGdvy_surf_CC_vx[12]/dv1; 
  surft1_up[13] = dGdvy_surf_CC_vx[13]/dv1; 
  surft1_up[14] = dGdvy_surf_CC_vx[14]/dv1; 
  surft1_up[15] = dGdvy_surf_CC_vx[15]/dv1; 

  surft2_lo[0] = -(0.408248290463863*GCC[4])+0.408248290463863*GBC[4]+0.3535533905932737*(GCC[0]+GBC[0]); 
  surft2_lo[1] = -(0.408248290463863*GCC[9])+0.408248290463863*GBC[9]+0.3535533905932737*(GCC[1]+GBC[1]); 
  surft2_lo[2] = -(0.408248290463863*GCC[10])+0.408248290463863*GBC[10]+0.3535533905932737*(GCC[2]+GBC[2]); 
  surft2_lo[3] = -(0.408248290463863*GCC[11])+0.408248290463863*GBC[11]+0.3535533905932737*(GCC[3]+GBC[3]); 
  surft2_lo[4] = -(0.408248290463863*GCC[15])+0.408248290463863*GBC[15]+0.3535533905932737*(GCC[5]+GBC[5]); 
  surft2_lo[5] = -(0.408248290463863*GCC[17])+0.408248290463863*GBC[17]+0.3535533905932737*(GCC[6]+GBC[6]); 
  surft2_lo[6] = -(0.408248290463863*GCC[18])+0.408248290463863*GBC[18]+0.3535533905932737*(GCC[7]+GBC[7]); 
  surft2_lo[7] = -(0.408248290463863*GCC[19])+0.408248290463863*GBC[19]+0.3535533905932737*(GCC[8]+GBC[8]); 
  surft2_lo[8] = -(0.408248290463863*GCC[23])+0.408248290463863*GBC[23]+0.3535533905932737*(GCC[12]+GBC[12]); 
  surft2_lo[9] = -(0.408248290463863*GCC[24])+0.408248290463863*GBC[24]+0.3535533905932737*(GCC[13]+GBC[13]); 
  surft2_lo[10] = -(0.408248290463863*GCC[25])+0.408248290463863*GBC[25]+0.3535533905932737*(GCC[14]+GBC[14]); 
  surft2_lo[11] = -(0.408248290463863*GCC[26])+0.408248290463863*GBC[26]+0.3535533905932737*(GCC[16]+GBC[16]); 
  surft2_lo[12] = -(0.408248290463863*GCC[28])+0.408248290463863*GBC[28]+0.3535533905932737*(GCC[20]+GBC[20]); 
  surft2_lo[13] = -(0.408248290463863*GCC[29])+0.408248290463863*GBC[29]+0.3535533905932737*(GCC[21]+GBC[21]); 
  surft2_lo[14] = -(0.408248290463863*GCC[30])+0.408248290463863*GBC[30]+0.3535533905932737*(GCC[22]+GBC[22]); 
  surft2_lo[15] = -(0.408248290463863*GCC[31])+0.408248290463863*GBC[31]+0.3535533905932737*(GCC[27]+GBC[27]); 
  surft2_up[0] = G_surf_CC_vy[0]; 
  surft2_up[1] = G_surf_CC_vy[1]; 
  surft2_up[2] = G_surf_CC_vy[2]; 
  surft2_up[3] = G_surf_CC_vy[3]; 
  surft2_up[4] = G_surf_CC_vy[4]; 
  surft2_up[5] = G_surf_CC_vy[5]; 
  surft2_up[6] = G_surf_CC_vy[6]; 
  surft2_up[7] = G_surf_CC_vy[7]; 
  surft2_up[8] = G_surf_CC_vy[8]; 
  surft2_up[9] = G_surf_CC_vy[9]; 
  surft2_up[10] = G_surf_CC_vy[10]; 
  surft2_up[11] = G_surf_CC_vy[11]; 
  surft2_up[12] = G_surf_CC_vy[12]; 
  surft2_up[13] = G_surf_CC_vy[13]; 
  surft2_up[14] = G_surf_CC_vy[14]; 
  surft2_up[15] = G_surf_CC_vy[15]; 

  vol[11] = 3.0*GCC[0]; 
  vol[18] = 3.0*GCC[1]; 
  vol[19] = 3.0*GCC[2]; 
  vol[25] = 3.0*GCC[5]; 
  vol[26] = 3.0*GCC[6]; 
  vol[29] = 3.0*GCC[12]; 
  vol[30] = 3.0*GCC[13]; 
  vol[31] = 3.0*GCC[20]; 

  out[0] = (vol[0]+0.7071067811865475*surft1_up[0]-0.7071067811865475*surft1_lo[0])*dv1_sq*gamma_avg; 
  out[1] = (vol[1]+0.7071067811865475*surft1_up[1]-0.7071067811865475*surft1_lo[1])*dv1_sq*gamma_avg; 
  out[2] = (vol[2]+0.7071067811865475*surft1_up[2]-0.7071067811865475*surft1_lo[2])*dv1_sq*gamma_avg; 
  out[3] = (vol[3]-1.224744871391589*surft2_up[0]+1.224744871391589*(surft2_lo[0]+surft1_up[0]+surft1_lo[0]))*dv1_sq*gamma_avg; 
  out[4] = (vol[4]+0.7071067811865475*surft1_up[3]-0.7071067811865475*surft1_lo[3])*dv1_sq*gamma_avg; 
  out[5] = (vol[5]+0.7071067811865475*surft1_up[4]-0.7071067811865475*surft1_lo[4])*dv1_sq*gamma_avg; 
  out[6] = (vol[6]+0.7071067811865475*surft1_up[5]-0.7071067811865475*surft1_lo[5])*dv1_sq*gamma_avg; 
  out[7] = (vol[7]-1.224744871391589*surft2_up[1]+1.224744871391589*(surft2_lo[1]+surft1_up[1]+surft1_lo[1]))*dv1_sq*gamma_avg; 
  out[8] = (vol[8]-1.224744871391589*surft2_up[2]+1.224744871391589*(surft2_lo[2]+surft1_up[2]+surft1_lo[2]))*dv1_sq*gamma_avg; 
  out[9] = (vol[9]+0.7071067811865475*surft1_up[6]-0.7071067811865475*surft1_lo[6])*dv1_sq*gamma_avg; 
  out[10] = (vol[10]+0.7071067811865475*surft1_up[7]-0.7071067811865475*surft1_lo[7])*dv1_sq*gamma_avg; 
  out[11] = (vol[11]+1.224744871391589*(surft1_up[3]+surft1_lo[3])-2.1213203435596424*(surft2_up[0]+surft2_lo[0]))*dv1_sq*gamma_avg; 
  out[12] = (vol[12]+0.7071067811865475*surft1_up[8]-0.7071067811865475*surft1_lo[8])*dv1_sq*gamma_avg; 
  out[13] = (vol[13]+0.7071067811865475*surft1_up[9]-0.7071067811865475*surft1_lo[9])*dv1_sq*gamma_avg; 
  out[14] = (vol[14]-1.224744871391589*surft2_up[4]+1.224744871391589*(surft2_lo[4]+surft1_up[4]+surft1_lo[4]))*dv1_sq*gamma_avg; 
  out[15] = (vol[15]+0.7071067811865475*surft1_up[10]-0.7071067811865475*surft1_lo[10])*dv1_sq*gamma_avg; 
  out[16] = (vol[16]-1.224744871391589*surft2_up[5]+1.224744871391589*(surft2_lo[5]+surft1_up[5]+surft1_lo[5]))*dv1_sq*gamma_avg; 
  out[17] = (vol[17]+0.7071067811865475*surft1_up[11]-0.7071067811865475*surft1_lo[11])*dv1_sq*gamma_avg; 
  out[18] = (vol[18]+1.224744871391589*(surft1_up[6]+surft1_lo[6])-2.1213203435596424*(surft2_up[1]+surft2_lo[1]))*dv1_sq*gamma_avg; 
  out[19] = (vol[19]+1.224744871391589*(surft1_up[7]+surft1_lo[7])-2.1213203435596424*(surft2_up[2]+surft2_lo[2]))*dv1_sq*gamma_avg; 
  out[20] = (vol[20]+0.7071067811865475*surft1_up[12]-0.7071067811865475*surft1_lo[12])*dv1_sq*gamma_avg; 
  out[21] = (vol[21]-1.224744871391589*surft2_up[8]+1.224744871391589*(surft2_lo[8]+surft1_up[8]+surft1_lo[8]))*dv1_sq*gamma_avg; 
  out[22] = (vol[22]-1.224744871391589*surft2_up[9]+1.224744871391589*(surft2_lo[9]+surft1_up[9]+surft1_lo[9]))*dv1_sq*gamma_avg; 
  out[23] = (vol[23]+0.7071067811865475*surft1_up[13]-0.7071067811865475*surft1_lo[13])*dv1_sq*gamma_avg; 
  out[24] = (vol[24]+0.7071067811865475*surft1_up[14]-0.7071067811865475*surft1_lo[14])*dv1_sq*gamma_avg; 
  out[25] = (vol[25]+1.224744871391589*(surft1_up[10]+surft1_lo[10])-2.1213203435596424*(surft2_up[4]+surft2_lo[4]))*dv1_sq*gamma_avg; 
  out[26] = (vol[26]+1.224744871391589*(surft1_up[11]+surft1_lo[11])-2.1213203435596424*(surft2_up[5]+surft2_lo[5]))*dv1_sq*gamma_avg; 
  out[27] = (vol[27]-1.224744871391589*surft2_up[12]+1.224744871391589*(surft2_lo[12]+surft1_up[12]+surft1_lo[12]))*dv1_sq*gamma_avg; 
  out[28] = (vol[28]+0.7071067811865475*surft1_up[15]-0.7071067811865475*surft1_lo[15])*dv1_sq*gamma_avg; 
  out[29] = (vol[29]+1.224744871391589*(surft1_up[13]+surft1_lo[13])-2.1213203435596424*(surft2_up[8]+surft2_lo[8]))*dv1_sq*gamma_avg; 
  out[30] = (vol[30]+1.224744871391589*(surft1_up[14]+surft1_lo[14])-2.1213203435596424*(surft2_up[9]+surft2_lo[9]))*dv1_sq*gamma_avg; 
  out[31] = (vol[31]+1.224744871391589*(surft1_up[15]+surft1_lo[15])-2.1213203435596424*(surft2_up[12]+surft2_lo[12]))*dv1_sq*gamma_avg; 
} 
