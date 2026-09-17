#include <gkyl_fpo_vlasov_kernels.h> 
#include <math.h> 
 
void fpo_diff_coeff_cross_1x3v_vxvz_ser_p1_upvx_lovz(const double *dxv, const double *gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff) { 
  // dxv[NDIM]: Cell spacing in each direction. 
  // gamma: Scalar factor gamma. 
  // fpo_g_stencil[9]: 9 cell stencil of Rosenbluth potential G. 
  // fpo_g_surf_stencil[9]: 9 cell stencil of surface projection of G. 
  // fpo_dgdv_surf: Surface expansion of dG/dv in center cell. 
  // diff_coeff: Output array for diffusion tensor. 

  // Use cell-average value for gamma. 
  double gamma_avg = gamma[0]/sqrt(pow(2, 1)); 
  double dv1 = 2.0/dxv[1]; 
  double dv2 = 2.0/dxv[3]; 
  double dv1_sq = 4.0/dxv[1]/dxv[3]; 
 
  double surft1_lo[8] = {0.0}; 
  double surft1_up[8] = {0.0}; 
  double surft2_lo[8] = {0.0}; 
  double surft2_up[8] = {0.0}; 
  double vol[16] = {0.0}; 
  double *out = &diff_coeff[32]; 

  const double* GBL = fpo_g_stencil[0]; 
  const double* GCL = fpo_g_stencil[1]; 
  const double* GTL = fpo_g_stencil[2]; 
  const double* GBC = fpo_g_stencil[3]; 
  const double* GCC = fpo_g_stencil[4]; 
  const double* G_surf_CC_vx = &fpo_g_surf_stencil[4][0]; 
  const double* GTC = fpo_g_stencil[5]; 
  const double* GBR = fpo_g_stencil[6]; 
  const double* GCR = fpo_g_stencil[7]; 
  const double* GTR = fpo_g_stencil[8]; 
  const double* dGdvx_surf_CC_vz = &fpo_dgdv_surf[32]; 

  surft1_lo[0] = dGdvx_surf_CC_vz[0]/dv1; 
  surft1_lo[1] = dGdvx_surf_CC_vz[1]/dv1; 
  surft1_lo[2] = dGdvx_surf_CC_vz[2]/dv1; 
  surft1_lo[3] = dGdvx_surf_CC_vz[3]/dv1; 
  surft1_lo[4] = dGdvx_surf_CC_vz[4]/dv1; 
  surft1_lo[5] = dGdvx_surf_CC_vz[5]/dv1; 
  surft1_lo[6] = dGdvx_surf_CC_vz[6]/dv1; 
  surft1_lo[7] = dGdvx_surf_CC_vz[7]/dv1; 
  surft1_up[0] = -(0.7071067811865475*GTC[9])+0.7071067811865475*GCC[9]+0.6123724356957944*(GTC[2]+GCC[2]); 
  surft1_up[1] = -(0.7071067811865475*GTC[12])+0.7071067811865475*GCC[12]+0.6123724356957944*(GTC[5]+GCC[5]); 
  surft1_up[3] = -(0.7071067811865475*GTC[14])+0.7071067811865475*GCC[14]+0.6123724356957944*(GTC[7]+GCC[7]); 
  surft1_up[5] = -(0.7071067811865475*GTC[15])+0.7071067811865475*GCC[15]+0.6123724356957944*(GTC[11]+GCC[11]); 

  surft2_lo[0] = 0.408248290463863*GCL[2]-0.408248290463863*GCC[2]+0.3535533905932737*(GCL[0]+GCC[0]); 
  surft2_lo[1] = 0.408248290463863*GCL[5]-0.408248290463863*GCC[5]+0.3535533905932737*(GCL[1]+GCC[1]); 
  surft2_lo[2] = 0.408248290463863*GCL[7]-0.408248290463863*GCC[7]+0.3535533905932737*(GCL[3]+GCC[3]); 
  surft2_lo[3] = 0.408248290463863*GCL[9]-0.408248290463863*GCC[9]+0.3535533905932737*(GCL[4]+GCC[4]); 
  surft2_lo[4] = 0.408248290463863*GCL[11]-0.408248290463863*GCC[11]+0.3535533905932737*(GCL[6]+GCC[6]); 
  surft2_lo[5] = 0.408248290463863*GCL[12]-0.408248290463863*GCC[12]+0.3535533905932737*(GCL[8]+GCC[8]); 
  surft2_lo[6] = 0.408248290463863*GCL[14]-0.408248290463863*GCC[14]+0.3535533905932737*(GCL[10]+GCC[10]); 
  surft2_lo[7] = 0.408248290463863*GCL[15]-0.408248290463863*GCC[15]+0.3535533905932737*(GCL[13]+GCC[13]); 
  surft2_up[0] = G_surf_CC_vx[0]; 
  surft2_up[1] = G_surf_CC_vx[1]; 
  surft2_up[2] = G_surf_CC_vx[2]; 
  surft2_up[3] = G_surf_CC_vx[3]; 
  surft2_up[4] = G_surf_CC_vx[4]; 
  surft2_up[5] = G_surf_CC_vx[5]; 
  surft2_up[6] = G_surf_CC_vx[6]; 
  surft2_up[7] = G_surf_CC_vx[7]; 

  vol[9] = 3.0*GCC[0]; 
  vol[12] = 3.0*GCC[1]; 
  vol[14] = 3.0*GCC[3]; 
  vol[15] = 3.0*GCC[6]; 

  out[0] = (vol[0]+0.7071067811865475*surft1_up[0]-0.7071067811865475*surft1_lo[0])*dv1_sq*gamma_avg; 
  out[1] = (vol[1]+0.7071067811865475*surft1_up[1]-0.7071067811865475*surft1_lo[1])*dv1_sq*gamma_avg; 
  out[2] = (vol[2]+0.7071067811865475*surft1_up[2]-0.7071067811865475*surft1_lo[2])*dv1_sq*gamma_avg; 
  out[3] = (vol[3]+0.7071067811865475*surft1_up[3]-0.7071067811865475*surft1_lo[3])*dv1_sq*gamma_avg; 
  out[4] = (vol[4]-1.224744871391589*surft2_up[0]+1.224744871391589*(surft2_lo[0]+surft1_up[0]+surft1_lo[0]))*dv1_sq*gamma_avg; 
  out[5] = (vol[5]+0.7071067811865475*surft1_up[4]-0.7071067811865475*surft1_lo[4])*dv1_sq*gamma_avg; 
  out[6] = (vol[6]+0.7071067811865475*surft1_up[5]-0.7071067811865475*surft1_lo[5])*dv1_sq*gamma_avg; 
  out[7] = (vol[7]+0.7071067811865475*surft1_up[6]-0.7071067811865475*surft1_lo[6])*dv1_sq*gamma_avg; 
  out[8] = (vol[8]-1.224744871391589*surft2_up[1]+1.224744871391589*(surft2_lo[1]+surft1_up[1]+surft1_lo[1]))*dv1_sq*gamma_avg; 
  out[9] = (vol[9]+1.224744871391589*(surft1_up[2]+surft1_lo[2])-2.1213203435596424*(surft2_up[0]+surft2_lo[0]))*dv1_sq*gamma_avg; 
  out[10] = (vol[10]+1.224744871391589*(surft1_up[3]+surft1_lo[3]+surft2_lo[2])-1.224744871391589*surft2_up[2])*dv1_sq*gamma_avg; 
  out[11] = (vol[11]+0.7071067811865475*surft1_up[7]-0.7071067811865475*surft1_lo[7])*dv1_sq*gamma_avg; 
  out[12] = (vol[12]+1.224744871391589*(surft1_up[4]+surft1_lo[4])-2.1213203435596424*(surft2_up[1]+surft2_lo[1]))*dv1_sq*gamma_avg; 
  out[13] = (vol[13]+1.224744871391589*(surft1_up[5]+surft1_lo[5]+surft2_lo[4])-1.224744871391589*surft2_up[4])*dv1_sq*gamma_avg; 
  out[14] = (vol[14]+1.224744871391589*(surft1_up[6]+surft1_lo[6])-2.1213203435596424*(surft2_up[2]+surft2_lo[2]))*dv1_sq*gamma_avg; 
  out[15] = (vol[15]+1.224744871391589*(surft1_up[7]+surft1_lo[7])-2.1213203435596424*(surft2_up[4]+surft2_lo[4]))*dv1_sq*gamma_avg; 
} 
