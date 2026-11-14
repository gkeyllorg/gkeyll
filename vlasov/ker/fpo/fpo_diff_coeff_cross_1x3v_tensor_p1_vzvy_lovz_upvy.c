#include <gkyl_fpo_vlasov_kernels.h> 
 
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvy_tensor_p1_lovz_upvy(const double *dxv, const double *gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff) { 
  // dxv[NDIM]: Cell spacing in each direction. 
  // gamma: Scalar factor gamma. 
  // fpo_g_stencil[9]: 9 cell stencil of Rosenbluth potential G. 
  // fpo_g_surf_stencil[9]: 9 cell stencil of surface projection of G. 
  // fpo_dgdv_surf: Surface expansion of dG/dv in center cell. 
  // diff_coeff: Output array for diffusion tensor. 

  // Use cell-average value for gamma. 
  double gamma_avg = gamma[0]/sqrt(pow(2, 1)); 
  double dv1 = 2.0/dxv[3]; 
  double dv2 = 2.0/dxv[2]; 
  double dv1_sq = 4.0/dxv[3]/dxv[2]; 
 
  double surft1_lo[18] = {0.0}; 
  double surft1_up[18] = {0.0}; 
  double surft2_lo[18] = {0.0}; 
  double surft2_up[18] = {0.0}; 
  double vol[54] = {0.0}; 
  double *out = &diff_coeff[378]; 

  const double* GBL = fpo_g_stencil[0]; 
  const double* GCL = fpo_g_stencil[1]; 
  const double* GTL = fpo_g_stencil[2]; 
  const double* GBC = fpo_g_stencil[3]; 
  const double* GCC = fpo_g_stencil[4]; 
  const double* G_surf_CC_vz = &fpo_g_surf_stencil[4][36]; 
  const double* GTC = fpo_g_stencil[5]; 
  const double* GBR = fpo_g_stencil[6]; 
  const double* GCR = fpo_g_stencil[7]; 
  const double* GTR = fpo_g_stencil[8]; 
  const double* dGdvz_surf_CC_vy = &fpo_dgdv_surf[54]; 

  surft1_lo[0] = 0.599071547271275*(GCL[23]+GCC[23])+0.8617863895711042*GCL[10]-0.8617863895711042*GCC[10]+0.6123724356957944*(GCL[4]+GCC[4]); 
  surft1_lo[1] = 0.599071547271275*(GCL[32]+GCC[32])+0.8617863895711042*GCL[16]-0.8617863895711042*GCC[16]+0.6123724356957944*(GCL[8]+GCC[8]); 
  surft1_lo[2] = 0.599071547271275*(GCL[33]+GCC[33])+0.8617863895711042*GCL[17]-0.8617863895711042*GCC[17]+0.6123724356957944*(GCL[9]+GCC[9]); 
  surft1_lo[3] = 1.3395647030845494*(GCL[39]+GCC[39])+1.9270129491651047*GCL[26]-1.9270129491651047*GCC[26]+1.369306393762915*(GCL[13]+GCC[13]); 
  surft1_lo[4] = 0.599071547271275*(GCL[41]+GCC[41])+0.8617863895711042*GCL[27]-0.8617863895711042*GCC[27]+0.6123724356957944*(GCL[15]+GCC[15]); 
  surft1_lo[5] = 1.3395647030845494*(GCL[47]+GCC[47])+1.927012949165105*GCL[35]-1.927012949165105*GCC[35]+1.369306393762915*(GCL[24]+GCC[24]); 
  surft1_lo[6] = 1.3395647030845494*(GCL[48]+GCC[48])+1.927012949165105*GCL[36]-1.927012949165105*GCC[36]+1.369306393762915*(GCL[25]+GCC[25]); 
  surft1_lo[7] = 0.599071547271275*(GCL[44]+GCC[44])+0.8617863895711042*GCL[31]-0.8617863895711042*GCC[31]+0.6123724356957944*(GCL[22]+GCC[22]); 
  surft1_lo[9] = 1.3395647030845494*(GCL[51]+GCC[51])+1.9270129491651047*GCL[42]-1.9270129491651047*GCC[42]+1.369306393762915*(GCL[34]+GCC[34]); 
  surft1_lo[10] = 0.599071547271275*(GCL[49]+GCC[49])+0.8617863895711042*GCL[40]-0.8617863895711042*GCC[40]+0.6123724356957944*(GCL[30]+GCC[30]); 
  surft1_lo[11] = 1.3395647030845497*(GCL[52]+GCC[52])+1.9270129491651047*GCL[46]-1.9270129491651047*GCC[46]+1.369306393762915*(GCL[38]+GCC[38]); 
  surft1_lo[14] = 1.3395647030845497*(GCL[53]+GCC[53])+1.927012949165105*GCL[50]-1.927012949165105*GCC[50]+1.369306393762915*(GCL[45]+GCC[45]); 
  surft1_up[0] = dGdvz_surf_CC_vy[0]/dv1; 
  surft1_up[1] = dGdvz_surf_CC_vy[1]/dv1; 
  surft1_up[2] = dGdvz_surf_CC_vy[2]/dv1; 
  surft1_up[3] = dGdvz_surf_CC_vy[3]/dv1; 
  surft1_up[4] = dGdvz_surf_CC_vy[4]/dv1; 
  surft1_up[5] = dGdvz_surf_CC_vy[5]/dv1; 
  surft1_up[6] = dGdvz_surf_CC_vy[6]/dv1; 
  surft1_up[7] = dGdvz_surf_CC_vy[7]/dv1; 
  surft1_up[8] = dGdvz_surf_CC_vy[8]/dv1; 
  surft1_up[9] = dGdvz_surf_CC_vy[9]/dv1; 
  surft1_up[10] = dGdvz_surf_CC_vy[10]/dv1; 
  surft1_up[11] = dGdvz_surf_CC_vy[11]/dv1; 
  surft1_up[12] = dGdvz_surf_CC_vy[12]/dv1; 
  surft1_up[13] = dGdvz_surf_CC_vy[13]/dv1; 
  surft1_up[14] = dGdvz_surf_CC_vy[14]/dv1; 
  surft1_up[15] = dGdvz_surf_CC_vy[15]/dv1; 
  surft1_up[16] = dGdvz_surf_CC_vy[16]/dv1; 
  surft1_up[17] = dGdvz_surf_CC_vy[17]/dv1; 

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
  surft2_lo[16] = G_surf_CC_vz[16]; 
  surft2_lo[17] = G_surf_CC_vz[17]; 
  surft2_up[0] = 0.34587411908091625*(GTC[13]+GCC[13])-0.49755260400283263*GTC[4]+0.49755260400283263*GCC[4]+0.3535533905932737*(GTC[0]+GCC[0]); 
  surft2_up[1] = 0.34587411908091625*(GTC[24]+GCC[24])-0.49755260400283263*GTC[8]+0.49755260400283263*GCC[8]+0.3535533905932737*(GTC[1]+GCC[1]); 
  surft2_up[2] = 0.34587411908091625*(GTC[25]+GCC[25])-0.49755260400283263*GTC[9]+0.49755260400283263*GCC[9]+0.3535533905932737*(GTC[2]+GCC[2]); 
  surft2_up[3] = 0.34587411908091625*(GTC[26]+GCC[26])-0.49755260400283263*GTC[10]+0.49755260400283263*GCC[10]+0.3535533905932737*(GTC[3]+GCC[3]); 
  surft2_up[4] = 0.34587411908091625*(GTC[34]+GCC[34])-0.49755260400283263*GTC[15]+0.49755260400283263*GCC[15]+0.3535533905932737*(GTC[5]+GCC[5]); 
  surft2_up[5] = 0.34587411908091625*(GTC[35]+GCC[35])-0.49755260400283263*GTC[16]+0.49755260400283263*GCC[16]+0.3535533905932737*(GTC[6]+GCC[6]); 
  surft2_up[6] = 0.34587411908091625*(GTC[36]+GCC[36])-0.49755260400283263*GTC[17]+0.49755260400283263*GCC[17]+0.3535533905932737*(GTC[7]+GCC[7]); 
  surft2_up[7] = 0.34587411908091625*(GTC[38]+GCC[38])-0.49755260400283263*GTC[22]+0.49755260400283263*GCC[22]+0.3535533905932737*(GTC[11]+GCC[11]); 
  surft2_up[8] = 0.34587411908091625*(GTC[39]+GCC[39])-0.49755260400283263*GTC[23]+0.49755260400283263*GCC[23]+0.3535533905932737*(GTC[12]+GCC[12]); 
  surft2_up[9] = 0.34587411908091625*(GTC[42]+GCC[42])-0.49755260400283263*GTC[27]+0.49755260400283263*GCC[27]+0.3535533905932737*(GTC[14]+GCC[14]); 
  surft2_up[10] = 0.34587411908091625*(GTC[45]+GCC[45])-0.49755260400283263*GTC[30]+0.49755260400283263*GCC[30]+0.3535533905932737*(GTC[18]+GCC[18]); 
  surft2_up[11] = 0.34587411908091625*(GTC[46]+GCC[46])-0.49755260400283263*GTC[31]+0.49755260400283263*GCC[31]+0.3535533905932737*(GTC[19]+GCC[19]); 
  surft2_up[12] = 0.34587411908091625*(GTC[47]+GCC[47])-0.49755260400283263*GTC[32]+0.49755260400283263*GCC[32]+0.3535533905932737*(GTC[20]+GCC[20]); 
  surft2_up[13] = 0.34587411908091625*(GTC[48]+GCC[48])-0.49755260400283263*GTC[33]+0.49755260400283263*GCC[33]+0.3535533905932737*(GTC[21]+GCC[21]); 
  surft2_up[14] = 0.34587411908091625*(GTC[50]+GCC[50])-0.49755260400283263*GTC[40]+0.49755260400283263*GCC[40]+0.3535533905932737*(GTC[28]+GCC[28]); 
  surft2_up[15] = 0.34587411908091625*(GTC[51]+GCC[51])-0.49755260400283263*GTC[41]+0.49755260400283263*GCC[41]+0.3535533905932737*(GTC[29]+GCC[29]); 
  surft2_up[16] = 0.34587411908091625*(GTC[52]+GCC[52])-0.49755260400283263*GTC[44]+0.49755260400283263*GCC[44]+0.3535533905932737*(GTC[37]+GCC[37]); 
  surft2_up[17] = 0.34587411908091625*(GTC[53]+GCC[53])-0.49755260400283263*GTC[49]+0.49755260400283263*GCC[49]+0.3535533905932737*(GTC[43]+GCC[43]); 

  vol[10] = 3.0*GCC[0]; 
  vol[16] = 3.0*GCC[1]; 
  vol[17] = 3.0*GCC[2]; 
  vol[23] = 6.7082039324993685*GCC[3]; 
  vol[26] = 6.7082039324993685*GCC[4]; 
  vol[27] = 3.0*GCC[5]; 
  vol[31] = 3.0*GCC[11]; 
  vol[32] = 6.708203932499369*GCC[6]; 
  vol[33] = 6.708203932499369*GCC[7]; 
  vol[35] = 6.708203932499369*GCC[8]; 
  vol[36] = 6.708203932499369*GCC[9]; 
  vol[39] = 15.0*GCC[10]; 
  vol[40] = 3.0*GCC[18]; 
  vol[41] = 6.7082039324993685*GCC[14]; 
  vol[42] = 6.7082039324993685*GCC[15]; 
  vol[44] = 6.7082039324993685*GCC[19]; 
  vol[46] = 6.7082039324993685*GCC[22]; 
  vol[47] = 15.0*GCC[16]; 
  vol[48] = 15.0*GCC[17]; 
  vol[49] = 6.708203932499369*GCC[28]; 
  vol[50] = 6.708203932499369*GCC[30]; 
  vol[51] = 15.0*GCC[27]; 
  vol[52] = 15.0*GCC[31]; 
  vol[53] = 15.0*GCC[40]; 

  out[0] = (vol[0]+0.7071067811865475*surft1_up[0]-0.7071067811865475*surft1_lo[0])*dv1_sq*gamma_avg; 
  out[1] = (vol[1]+0.7071067811865475*surft1_up[1]-0.7071067811865475*surft1_lo[1])*dv1_sq*gamma_avg; 
  out[2] = (vol[2]+0.7071067811865475*surft1_up[2]-0.7071067811865475*surft1_lo[2])*dv1_sq*gamma_avg; 
  out[3] = (vol[3]-1.224744871391589*surft2_up[0]+1.224744871391589*(surft2_lo[0]+surft1_up[0]+surft1_lo[0]))*dv1_sq*gamma_avg; 
  out[4] = (vol[4]+0.7071067811865475*surft1_up[3]-0.7071067811865475*surft1_lo[3])*dv1_sq*gamma_avg; 
  out[5] = (vol[5]+0.7071067811865475*surft1_up[4]-0.7071067811865475*surft1_lo[4])*dv1_sq*gamma_avg; 
  out[6] = (vol[6]-1.224744871391589*surft2_up[1]+1.224744871391589*(surft2_lo[1]+surft1_up[1]+surft1_lo[1]))*dv1_sq*gamma_avg; 
  out[7] = (vol[7]-1.224744871391589*surft2_up[2]+1.224744871391589*(surft2_lo[2]+surft1_up[2]+surft1_lo[2]))*dv1_sq*gamma_avg; 
  out[8] = (vol[8]+0.7071067811865475*surft1_up[5]-0.7071067811865475*surft1_lo[5])*dv1_sq*gamma_avg; 
  out[9] = (vol[9]+0.7071067811865475*surft1_up[6]-0.7071067811865475*surft1_lo[6])*dv1_sq*gamma_avg; 
  out[10] = (vol[10]+1.224744871391589*(surft1_up[3]+surft1_lo[3])-2.1213203435596424*(surft2_up[0]+surft2_lo[0]))*dv1_sq*gamma_avg; 
  out[11] = (vol[11]+0.7071067811865475*surft1_up[7]-0.7071067811865475*surft1_lo[7])*dv1_sq*gamma_avg; 
  out[12] = (vol[12]-2.7386127875258306*surft2_up[3]+2.7386127875258306*surft2_lo[3]+1.5811388300841895*surft1_up[0]-1.5811388300841895*surft1_lo[0])*dv1_sq*gamma_avg; 
  out[13] = (vol[13]+0.7071067811865475*surft1_up[8]-0.7071067811865475*surft1_lo[8])*dv1_sq*gamma_avg; 
  out[14] = (vol[14]-1.224744871391589*surft2_up[4]+1.224744871391589*(surft2_lo[4]+surft1_up[4]+surft1_lo[4]))*dv1_sq*gamma_avg; 
  out[15] = (vol[15]+0.7071067811865475*surft1_up[9]-0.7071067811865475*surft1_lo[9])*dv1_sq*gamma_avg; 
  out[16] = (vol[16]+1.224744871391589*(surft1_up[5]+surft1_lo[5])-2.1213203435596424*(surft2_up[1]+surft2_lo[1]))*dv1_sq*gamma_avg; 
  out[17] = (vol[17]+1.224744871391589*(surft1_up[6]+surft1_lo[6])-2.1213203435596424*(surft2_up[2]+surft2_lo[2]))*dv1_sq*gamma_avg; 
  out[18] = (vol[18]+0.7071067811865475*surft1_up[10]-0.7071067811865475*surft1_lo[10])*dv1_sq*gamma_avg; 
  out[19] = (vol[19]-1.224744871391589*surft2_up[7]+1.224744871391589*(surft2_lo[7]+surft1_up[7]+surft1_lo[7]))*dv1_sq*gamma_avg; 
  out[20] = (vol[20]-2.7386127875258306*surft2_up[5]+2.7386127875258306*surft2_lo[5]+1.5811388300841898*surft1_up[1]-1.5811388300841898*surft1_lo[1])*dv1_sq*gamma_avg; 
  out[21] = (vol[21]-2.7386127875258306*surft2_up[6]+2.7386127875258306*surft2_lo[6]+1.5811388300841898*surft1_up[2]-1.5811388300841898*surft1_lo[2])*dv1_sq*gamma_avg; 
  out[22] = (vol[22]+0.7071067811865475*surft1_up[11]-0.7071067811865475*surft1_lo[11])*dv1_sq*gamma_avg; 
  out[23] = (vol[23]-4.743416490252569*(surft2_up[3]+surft2_lo[3])+1.5811388300841898*surft1_up[3]-1.5811388300841898*surft1_lo[3])*dv1_sq*gamma_avg; 
  out[24] = (vol[24]+0.7071067811865475*surft1_up[12]-0.7071067811865475*surft1_lo[12])*dv1_sq*gamma_avg; 
  out[25] = (vol[25]+0.7071067811865475*surft1_up[13]-0.7071067811865475*surft1_lo[13])*dv1_sq*gamma_avg; 
  out[26] = (vol[26]+1.224744871391589*(surft1_up[8]+surft1_lo[8])-2.7386127875258306*surft2_up[0]+2.7386127875258306*surft2_lo[0])*dv1_sq*gamma_avg; 
  out[27] = (vol[27]+1.224744871391589*(surft1_up[9]+surft1_lo[9])-2.1213203435596424*(surft2_up[4]+surft2_lo[4]))*dv1_sq*gamma_avg; 
  out[28] = (vol[28]-1.224744871391589*surft2_up[10]+1.224744871391589*(surft2_lo[10]+surft1_up[10]+surft1_lo[10]))*dv1_sq*gamma_avg; 
  out[29] = (vol[29]-2.7386127875258306*surft2_up[9]+2.7386127875258306*surft2_lo[9]+1.5811388300841895*surft1_up[4]-1.5811388300841895*surft1_lo[4])*dv1_sq*gamma_avg; 
  out[30] = (vol[30]+0.7071067811865475*surft1_up[14]-0.7071067811865475*surft1_lo[14])*dv1_sq*gamma_avg; 
  out[31] = (vol[31]+1.224744871391589*(surft1_up[11]+surft1_lo[11])-2.1213203435596424*(surft2_up[7]+surft2_lo[7]))*dv1_sq*gamma_avg; 
  out[32] = (vol[32]-4.743416490252569*(surft2_up[5]+surft2_lo[5])+1.5811388300841895*surft1_up[5]-1.5811388300841895*surft1_lo[5])*dv1_sq*gamma_avg; 
  out[33] = (vol[33]-4.743416490252569*(surft2_up[6]+surft2_lo[6])+1.5811388300841895*surft1_up[6]-1.5811388300841895*surft1_lo[6])*dv1_sq*gamma_avg; 
  out[34] = (vol[34]+0.7071067811865475*surft1_up[15]-0.7071067811865475*surft1_lo[15])*dv1_sq*gamma_avg; 
  out[35] = (vol[35]+1.224744871391589*(surft1_up[12]+surft1_lo[12])-2.7386127875258306*surft2_up[1]+2.7386127875258306*surft2_lo[1])*dv1_sq*gamma_avg; 
  out[36] = (vol[36]+1.224744871391589*(surft1_up[13]+surft1_lo[13])-2.7386127875258306*surft2_up[2]+2.7386127875258306*surft2_lo[2])*dv1_sq*gamma_avg; 
  out[37] = (vol[37]-2.7386127875258306*surft2_up[11]+2.7386127875258306*surft2_lo[11]+1.5811388300841895*surft1_up[7]-1.5811388300841895*surft1_lo[7])*dv1_sq*gamma_avg; 
  out[38] = (vol[38]+0.7071067811865475*surft1_up[16]-0.7071067811865475*surft1_lo[16])*dv1_sq*gamma_avg; 
  out[39] = (vol[39]+1.5811388300841895*surft1_up[8]-1.5811388300841895*surft1_lo[8]-6.123724356957944*surft2_up[3]+6.123724356957944*surft2_lo[3])*dv1_sq*gamma_avg; 
  out[40] = (vol[40]+1.224744871391589*(surft1_up[14]+surft1_lo[14])-2.1213203435596424*(surft2_up[10]+surft2_lo[10]))*dv1_sq*gamma_avg; 
  out[41] = (vol[41]-4.743416490252569*(surft2_up[9]+surft2_lo[9])+1.5811388300841898*surft1_up[9]-1.5811388300841898*surft1_lo[9])*dv1_sq*gamma_avg; 
  out[42] = (vol[42]+1.224744871391589*(surft1_up[15]+surft1_lo[15])-2.7386127875258306*surft2_up[4]+2.7386127875258306*surft2_lo[4])*dv1_sq*gamma_avg; 
  out[43] = (vol[43]-2.7386127875258306*surft2_up[14]+2.7386127875258306*surft2_lo[14]+1.5811388300841898*surft1_up[10]-1.5811388300841898*surft1_lo[10])*dv1_sq*gamma_avg; 
  out[44] = (vol[44]-4.743416490252569*(surft2_up[11]+surft2_lo[11])+1.5811388300841898*surft1_up[11]-1.5811388300841898*surft1_lo[11])*dv1_sq*gamma_avg; 
  out[45] = (vol[45]+0.7071067811865475*surft1_up[17]-0.7071067811865475*surft1_lo[17])*dv1_sq*gamma_avg; 
  out[46] = (vol[46]+1.224744871391589*(surft1_up[16]+surft1_lo[16])-2.7386127875258306*surft2_up[7]+2.7386127875258306*surft2_lo[7])*dv1_sq*gamma_avg; 
  out[47] = (vol[47]+1.5811388300841898*surft1_up[12]-1.5811388300841898*surft1_lo[12]-6.123724356957944*surft2_up[5]+6.123724356957944*surft2_lo[5])*dv1_sq*gamma_avg; 
  out[48] = (vol[48]+1.5811388300841898*surft1_up[13]-1.5811388300841898*surft1_lo[13]-6.123724356957944*surft2_up[6]+6.123724356957944*surft2_lo[6])*dv1_sq*gamma_avg; 
  out[49] = (vol[49]-4.743416490252569*(surft2_up[14]+surft2_lo[14])+1.5811388300841895*surft1_up[14]-1.5811388300841895*surft1_lo[14])*dv1_sq*gamma_avg; 
  out[50] = (vol[50]+1.224744871391589*(surft1_up[17]+surft1_lo[17])-2.7386127875258306*surft2_up[10]+2.7386127875258306*surft2_lo[10])*dv1_sq*gamma_avg; 
  out[51] = (vol[51]+1.5811388300841895*surft1_up[15]-1.5811388300841895*surft1_lo[15]-6.123724356957944*surft2_up[9]+6.123724356957944*surft2_lo[9])*dv1_sq*gamma_avg; 
  out[52] = (vol[52]+1.5811388300841895*surft1_up[16]-1.5811388300841895*surft1_lo[16]-6.123724356957945*surft2_up[11]+6.123724356957945*surft2_lo[11])*dv1_sq*gamma_avg; 
  out[53] = (vol[53]+1.5811388300841898*surft1_up[17]-1.5811388300841898*surft1_lo[17]-6.123724356957945*surft2_up[14]+6.123724356957945*surft2_lo[14])*dv1_sq*gamma_avg; 
} 
