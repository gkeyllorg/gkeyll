#include <gkyl_fpo_vlasov_kernels.h> 
 
GKYL_CU_DH void fpo_diff_coeff_cross_2x3v_vzvx_ser_p1_lovz_lovx(const double *dxv, const double *gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff) { 
  // dxv[NDIM]: Cell spacing in each direction. 
  // gamma: Scalar factor gamma. 
  // fpo_g_stencil[9]: 9 cell stencil of Rosenbluth potential G. 
  // fpo_g_surf_stencil[9]: 9 cell stencil of surface projection of G. 
  // fpo_dgdv_surf: Surface expansion of dG/dv in center cell. 
  // diff_coeff: Output array for diffusion tensor. 

  // Use cell-average value for gamma. 
  double gamma_avg = gamma[0]/sqrt(pow(2, 2)); 
  double dv1 = 2.0/dxv[4]; 
  double dv2 = 2.0/dxv[2]; 
  double dv1_sq = 4.0/dxv[4]/dxv[2]; 
 
  double surft1_lo[16] = {0.0}; 
  double surft1_up[16] = {0.0}; 
  double surft2_lo[16] = {0.0}; 
  double surft2_up[16] = {0.0}; 
  double vol[32] = {0.0}; 
  double *out = &diff_coeff[192]; 

  const double* GBL = fpo_g_stencil[0]; 
  const double* GCL = fpo_g_stencil[1]; 
  const double* GTL = fpo_g_stencil[2]; 
  const double* GBC = fpo_g_stencil[3]; 
  const double* GCC = fpo_g_stencil[4]; 
  const double* G_surf_CC_vz = &fpo_g_surf_stencil[4][32]; 
  const double* GTC = fpo_g_stencil[5]; 
  const double* GBR = fpo_g_stencil[6]; 
  const double* GCR = fpo_g_stencil[7]; 
  const double* GTR = fpo_g_stencil[8]; 
  const double* dGdvz_surf_CC_vx = &fpo_dgdv_surf[16]; 

  surft1_lo[0] = dGdvz_surf_CC_vx[0]/dv1; 
  surft1_lo[1] = dGdvz_surf_CC_vx[1]/dv1; 
  surft1_lo[2] = dGdvz_surf_CC_vx[2]/dv1; 
  surft1_lo[3] = dGdvz_surf_CC_vx[3]/dv1; 
  surft1_lo[4] = dGdvz_surf_CC_vx[4]/dv1; 
  surft1_lo[5] = dGdvz_surf_CC_vx[5]/dv1; 
  surft1_lo[6] = dGdvz_surf_CC_vx[6]/dv1; 
  surft1_lo[7] = dGdvz_surf_CC_vx[7]/dv1; 
  surft1_lo[8] = dGdvz_surf_CC_vx[8]/dv1; 
  surft1_lo[9] = dGdvz_surf_CC_vx[9]/dv1; 
  surft1_lo[10] = dGdvz_surf_CC_vx[10]/dv1; 
  surft1_lo[11] = dGdvz_surf_CC_vx[11]/dv1; 
  surft1_lo[12] = dGdvz_surf_CC_vx[12]/dv1; 
  surft1_lo[13] = dGdvz_surf_CC_vx[13]/dv1; 
  surft1_lo[14] = dGdvz_surf_CC_vx[14]/dv1; 
  surft1_lo[15] = dGdvz_surf_CC_vx[15]/dv1; 
  surft1_up[0] = -(0.7071067811865475*GCR[14])+0.7071067811865475*GCC[14]+0.6123724356957944*(GCR[5]+GCC[5]); 
  surft1_up[1] = -(0.7071067811865475*GCR[21])+0.7071067811865475*GCC[21]+0.6123724356957944*(GCR[12]+GCC[12]); 
  surft1_up[2] = -(0.7071067811865475*GCR[22])+0.7071067811865475*GCC[22]+0.6123724356957944*(GCR[13]+GCC[13]); 
  surft1_up[3] = -(0.7071067811865475*GCR[25])+0.7071067811865475*GCC[25]+0.6123724356957944*(GCR[15]+GCC[15]); 
  surft1_up[5] = -(0.7071067811865475*GCR[27])+0.7071067811865475*GCC[27]+0.6123724356957944*(GCR[20]+GCC[20]); 
  surft1_up[6] = -(0.7071067811865475*GCR[29])+0.7071067811865475*GCC[29]+0.6123724356957944*(GCR[23]+GCC[23]); 
  surft1_up[7] = -(0.7071067811865475*GCR[30])+0.7071067811865475*GCC[30]+0.6123724356957944*(GCR[24]+GCC[24]); 
  surft1_up[11] = -(0.7071067811865475*GCR[31])+0.7071067811865475*GCC[31]+0.6123724356957944*(GCR[28]+GCC[28]); 

  surft2_lo[0] = G_surf_CC_vz[0]; 
  surft2_lo[1] = G_surf_CC_vz[1]; 
  surft2_lo[2] = G_surf_CC_vz[2]; 
  surft2_lo[3] = G_surf_CC_vz[3]; 
  surft2_lo[4] = G_surf_CC_vz[4]; 
  surft2_lo[5] = G_surf_CC_vz[5]; 
  surft2_lo[6] = G_surf_CC_vz[6]; 
  surft2_lo[7] = G_surf_CC_vz[7]; 
  surft2_lo[8] = G_surf_CC_vz[8]; 
  surft2_lo[9] = G_surf_CC_vz[9]; 
  surft2_lo[10] = G_surf_CC_vz[10]; 
  surft2_lo[11] = G_surf_CC_vz[11]; 
  surft2_lo[12] = G_surf_CC_vz[12]; 
  surft2_lo[13] = G_surf_CC_vz[13]; 
  surft2_lo[14] = G_surf_CC_vz[14]; 
  surft2_lo[15] = G_surf_CC_vz[15]; 
  surft2_up[0] = -(0.408248290463863*GTC[5])+0.408248290463863*GCC[5]+0.3535533905932737*(GTC[0]+GCC[0]); 
  surft2_up[1] = -(0.408248290463863*GTC[12])+0.408248290463863*GCC[12]+0.3535533905932737*(GTC[1]+GCC[1]); 
  surft2_up[2] = -(0.408248290463863*GTC[13])+0.408248290463863*GCC[13]+0.3535533905932737*(GTC[2]+GCC[2]); 
  surft2_up[3] = -(0.408248290463863*GTC[14])+0.408248290463863*GCC[14]+0.3535533905932737*(GTC[3]+GCC[3]); 
  surft2_up[4] = -(0.408248290463863*GTC[15])+0.408248290463863*GCC[15]+0.3535533905932737*(GTC[4]+GCC[4]); 
  surft2_up[5] = -(0.408248290463863*GTC[20])+0.408248290463863*GCC[20]+0.3535533905932737*(GTC[6]+GCC[6]); 
  surft2_up[6] = -(0.408248290463863*GTC[21])+0.408248290463863*GCC[21]+0.3535533905932737*(GTC[7]+GCC[7]); 
  surft2_up[7] = -(0.408248290463863*GTC[22])+0.408248290463863*GCC[22]+0.3535533905932737*(GTC[8]+GCC[8]); 
  surft2_up[8] = -(0.408248290463863*GTC[23])+0.408248290463863*GCC[23]+0.3535533905932737*(GTC[9]+GCC[9]); 
  surft2_up[9] = -(0.408248290463863*GTC[24])+0.408248290463863*GCC[24]+0.3535533905932737*(GTC[10]+GCC[10]); 
  surft2_up[10] = -(0.408248290463863*GTC[25])+0.408248290463863*GCC[25]+0.3535533905932737*(GTC[11]+GCC[11]); 
  surft2_up[11] = -(0.408248290463863*GTC[27])+0.408248290463863*GCC[27]+0.3535533905932737*(GTC[16]+GCC[16]); 
  surft2_up[12] = -(0.408248290463863*GTC[28])+0.408248290463863*GCC[28]+0.3535533905932737*(GTC[17]+GCC[17]); 
  surft2_up[13] = -(0.408248290463863*GTC[29])+0.408248290463863*GCC[29]+0.3535533905932737*(GTC[18]+GCC[18]); 
  surft2_up[14] = -(0.408248290463863*GTC[30])+0.408248290463863*GCC[30]+0.3535533905932737*(GTC[19]+GCC[19]); 
  surft2_up[15] = -(0.408248290463863*GTC[31])+0.408248290463863*GCC[31]+0.3535533905932737*(GTC[26]+GCC[26]); 

  vol[14] = 3.0*GCC[0]; 
  vol[21] = 3.0*GCC[1]; 
  vol[22] = 3.0*GCC[2]; 
  vol[25] = 3.0*GCC[4]; 
  vol[27] = 3.0*GCC[6]; 
  vol[29] = 3.0*GCC[9]; 
  vol[30] = 3.0*GCC[10]; 
  vol[31] = 3.0*GCC[17]; 

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
  out[11] = (vol[11]-1.224744871391589*surft2_up[4]+1.224744871391589*(surft2_lo[4]+surft1_up[3]+surft1_lo[3]))*dv1_sq*gamma_avg; 
  out[12] = (vol[12]+0.7071067811865475*surft1_up[8]-0.7071067811865475*surft1_lo[8])*dv1_sq*gamma_avg; 
  out[13] = (vol[13]+0.7071067811865475*surft1_up[9]-0.7071067811865475*surft1_lo[9])*dv1_sq*gamma_avg; 
  out[14] = (vol[14]+1.224744871391589*(surft1_up[4]+surft1_lo[4])-2.1213203435596424*(surft2_up[0]+surft2_lo[0]))*dv1_sq*gamma_avg; 
  out[15] = (vol[15]+0.7071067811865475*surft1_up[10]-0.7071067811865475*surft1_lo[10])*dv1_sq*gamma_avg; 
  out[16] = (vol[16]-1.224744871391589*surft2_up[5]+1.224744871391589*(surft2_lo[5]+surft1_up[5]+surft1_lo[5]))*dv1_sq*gamma_avg; 
  out[17] = (vol[17]+0.7071067811865475*surft1_up[11]-0.7071067811865475*surft1_lo[11])*dv1_sq*gamma_avg; 
  out[18] = (vol[18]-1.224744871391589*surft2_up[8]+1.224744871391589*(surft2_lo[8]+surft1_up[6]+surft1_lo[6]))*dv1_sq*gamma_avg; 
  out[19] = (vol[19]-1.224744871391589*surft2_up[9]+1.224744871391589*(surft2_lo[9]+surft1_up[7]+surft1_lo[7]))*dv1_sq*gamma_avg; 
  out[20] = (vol[20]+0.7071067811865475*surft1_up[12]-0.7071067811865475*surft1_lo[12])*dv1_sq*gamma_avg; 
  out[21] = (vol[21]+1.224744871391589*(surft1_up[8]+surft1_lo[8])-2.1213203435596424*(surft2_up[1]+surft2_lo[1]))*dv1_sq*gamma_avg; 
  out[22] = (vol[22]+1.224744871391589*(surft1_up[9]+surft1_lo[9])-2.1213203435596424*(surft2_up[2]+surft2_lo[2]))*dv1_sq*gamma_avg; 
  out[23] = (vol[23]+0.7071067811865475*surft1_up[13]-0.7071067811865475*surft1_lo[13])*dv1_sq*gamma_avg; 
  out[24] = (vol[24]+0.7071067811865475*surft1_up[14]-0.7071067811865475*surft1_lo[14])*dv1_sq*gamma_avg; 
  out[25] = (vol[25]+1.224744871391589*(surft1_up[10]+surft1_lo[10])-2.1213203435596424*(surft2_up[4]+surft2_lo[4]))*dv1_sq*gamma_avg; 
  out[26] = (vol[26]-1.224744871391589*surft2_up[12]+1.224744871391589*(surft2_lo[12]+surft1_up[11]+surft1_lo[11]))*dv1_sq*gamma_avg; 
  out[27] = (vol[27]+1.224744871391589*(surft1_up[12]+surft1_lo[12])-2.1213203435596424*(surft2_up[5]+surft2_lo[5]))*dv1_sq*gamma_avg; 
  out[28] = (vol[28]+0.7071067811865475*surft1_up[15]-0.7071067811865475*surft1_lo[15])*dv1_sq*gamma_avg; 
  out[29] = (vol[29]+1.224744871391589*(surft1_up[13]+surft1_lo[13])-2.1213203435596424*(surft2_up[8]+surft2_lo[8]))*dv1_sq*gamma_avg; 
  out[30] = (vol[30]+1.224744871391589*(surft1_up[14]+surft1_lo[14])-2.1213203435596424*(surft2_up[9]+surft2_lo[9]))*dv1_sq*gamma_avg; 
  out[31] = (vol[31]+1.224744871391589*(surft1_up[15]+surft1_lo[15])-2.1213203435596424*(surft2_up[12]+surft2_lo[12]))*dv1_sq*gamma_avg; 
} 
