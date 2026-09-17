#include <gkyl_fpo_vlasov_kernels.h> 
 
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvx_ser_p1_invz_upvx(const double *dxv, const double *gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff) { 
  // dxv[NDIM]: Cell spacing in each direction. 
  // gamma: Scalar factor gamma. 
  // fpo_g_stencil[9]: 9 cell stencil of Rosenbluth potential G. 
  // fpo_g_surf_stencil[9]: 9 cell stencil of surface projection of G. 
  // fpo_dgdv_surf: Surface expansion of dG/dv in center cell. 
  // diff_coeff: Output array for diffusion tensor. 

  // Use cell-average value for gamma. 
  double gamma_avg = gamma[0]/sqrt(pow(2, 1)); 
  double dv1 = 2.0/dxv[3]; 
  double dv2 = 2.0/dxv[1]; 
  double dv1_sq = 4.0/dxv[3]/dxv[1]; 
 
  double surft1_lo[8] = {0.0}; 
  double surft1_up[8] = {0.0}; 
  double surft2_lo[8] = {0.0}; 
  double surft2_up[8] = {0.0}; 
  double vol[16] = {0.0}; 
  double *out = &diff_coeff[96]; 

  const double* GBL = fpo_g_stencil[0]; 
  const double* GCL = fpo_g_stencil[1]; 
  const double* GTL = fpo_g_stencil[2]; 
  const double* GBC = fpo_g_stencil[3]; 
  const double* GCC = fpo_g_stencil[4]; 
  const double* G_surf_CC_vz = &fpo_g_surf_stencil[4][16]; 
  const double* GTC = fpo_g_stencil[5]; 
  const double* GBR = fpo_g_stencil[6]; 
  const double* GCR = fpo_g_stencil[7]; 
  const double* GTR = fpo_g_stencil[8]; 
  const double* dGdvz_surf_CC_vx = &fpo_dgdv_surf[8]; 

  surft1_lo[0] = -(0.11785113019775789*(GTL[9]+GBL[9]))+0.11785113019775789*(GTC[9]+GBC[9])+0.2357022603955158*GCL[9]-0.2357022603955158*GCC[9]-0.10206207261596573*(GTL[4]+GTC[4]+GBL[4]+GBC[4]+GTC[2]+GBL[2])+0.20412414523193148*(GCL[4]+GCC[4])+0.10206207261596573*(GTL[2]+GBC[2])+0.0883883476483184*(GTL[0]+GTC[0])-0.0883883476483184*(GBL[0]+GBC[0]); 
  surft1_lo[1] = -(0.11785113019775789*(GTL[12]+GBL[12]))+0.11785113019775789*(GTC[12]+GBC[12])+0.2357022603955158*GCL[12]-0.2357022603955158*GCC[12]-0.10206207261596573*(GTL[8]+GTC[8]+GBL[8]+GBC[8]+GTC[5]+GBL[5])+0.20412414523193148*(GCL[8]+GCC[8])+0.10206207261596573*(GTL[5]+GBC[5])+0.0883883476483184*(GTL[1]+GTC[1])-0.0883883476483184*(GBL[1]+GBC[1]); 
  surft1_lo[2] = -(0.11785113019775789*(GTL[14]+GBL[14]))+0.11785113019775789*(GTC[14]+GBC[14])+0.2357022603955158*GCL[14]-0.2357022603955158*GCC[14]-0.10206207261596573*(GTL[10]+GTC[10]+GBL[10]+GBC[10]+GTC[7]+GBL[7])+0.20412414523193148*(GCL[10]+GCC[10])+0.10206207261596573*(GTL[7]+GBC[7])+0.0883883476483184*(GTL[3]+GTC[3])-0.0883883476483184*(GBL[3]+GBC[3]); 
  surft1_lo[3] = -(0.20412414523193148*(GTL[9]+GBC[9]))+0.20412414523193148*(GTC[9]+GBL[9])-0.1767766952966368*(GTL[4]+GTC[4]+GTC[2]+GBC[2])+0.1767766952966368*(GBL[4]+GBC[4]+GTL[2]+GBL[2])-0.3535533905932737*GCL[2]+0.3535533905932737*GCC[2]+0.15309310892394856*(GTL[0]+GTC[0]+GBL[0]+GBC[0])-0.3061862178478971*(GCL[0]+GCC[0]); 
  surft1_lo[4] = -(0.11785113019775789*(GTL[15]+GBL[15]))+0.11785113019775789*(GTC[15]+GBC[15])+0.2357022603955158*GCL[15]-0.2357022603955158*GCC[15]-0.10206207261596573*(GTL[13]+GTC[13]+GBL[13]+GBC[13]+GTC[11]+GBL[11])+0.20412414523193148*(GCL[13]+GCC[13])+0.10206207261596573*(GTL[11]+GBC[11])+0.0883883476483184*(GTL[6]+GTC[6])-0.0883883476483184*(GBL[6]+GBC[6]); 
  surft1_lo[5] = -(0.20412414523193148*(GTL[12]+GBC[12]))+0.20412414523193148*(GTC[12]+GBL[12])-0.1767766952966368*(GTL[8]+GTC[8]+GTC[5]+GBC[5])+0.1767766952966368*(GBL[8]+GBC[8]+GTL[5]+GBL[5])-0.3535533905932737*GCL[5]+0.3535533905932737*GCC[5]+0.15309310892394856*(GTL[1]+GTC[1]+GBL[1]+GBC[1])-0.3061862178478971*(GCL[1]+GCC[1]); 
  surft1_lo[6] = -(0.20412414523193148*(GTL[14]+GBC[14]))+0.20412414523193148*(GTC[14]+GBL[14])-0.1767766952966368*(GTL[10]+GTC[10]+GTC[7]+GBC[7])+0.1767766952966368*(GBL[10]+GBC[10]+GTL[7]+GBL[7])-0.3535533905932737*GCL[7]+0.3535533905932737*GCC[7]+0.15309310892394856*(GTL[3]+GTC[3]+GBL[3]+GBC[3])-0.3061862178478971*(GCL[3]+GCC[3]); 
  surft1_lo[7] = -(0.20412414523193148*(GTL[15]+GBC[15]))+0.20412414523193148*(GTC[15]+GBL[15])-0.1767766952966368*(GTL[13]+GTC[13]+GTC[11]+GBC[11])+0.1767766952966368*(GBL[13]+GBC[13]+GTL[11]+GBL[11])-0.3535533905932737*GCL[11]+0.3535533905932737*GCC[11]+0.15309310892394856*(GTL[6]+GTC[6]+GBL[6]+GBC[6])-0.3061862178478971*(GCL[6]+GCC[6]); 
  surft1_up[0] = dGdvz_surf_CC_vx[0]/dv1; 
  surft1_up[1] = dGdvz_surf_CC_vx[1]/dv1; 
  surft1_up[2] = dGdvz_surf_CC_vx[2]/dv1; 
  surft1_up[3] = dGdvz_surf_CC_vx[3]/dv1; 
  surft1_up[4] = dGdvz_surf_CC_vx[4]/dv1; 
  surft1_up[5] = dGdvz_surf_CC_vx[5]/dv1; 
  surft1_up[6] = dGdvz_surf_CC_vx[6]/dv1; 
  surft1_up[7] = dGdvz_surf_CC_vx[7]/dv1; 

  surft2_lo[0] = -(0.408248290463863*GCC[4])+0.408248290463863*GBC[4]+0.3535533905932737*(GCC[0]+GBC[0]); 
  surft2_lo[1] = -(0.408248290463863*GCC[8])+0.408248290463863*GBC[8]+0.3535533905932737*(GCC[1]+GBC[1]); 
  surft2_lo[2] = -(0.408248290463863*GCC[9])+0.408248290463863*GBC[9]+0.3535533905932737*(GCC[2]+GBC[2]); 
  surft2_lo[3] = -(0.408248290463863*GCC[10])+0.408248290463863*GBC[10]+0.3535533905932737*(GCC[3]+GBC[3]); 
  surft2_lo[4] = -(0.408248290463863*GCC[12])+0.408248290463863*GBC[12]+0.3535533905932737*(GCC[5]+GBC[5]); 
  surft2_lo[5] = -(0.408248290463863*GCC[13])+0.408248290463863*GBC[13]+0.3535533905932737*(GCC[6]+GBC[6]); 
  surft2_lo[6] = -(0.408248290463863*GCC[14])+0.408248290463863*GBC[14]+0.3535533905932737*(GCC[7]+GBC[7]); 
  surft2_lo[7] = -(0.408248290463863*GCC[15])+0.408248290463863*GBC[15]+0.3535533905932737*(GCC[11]+GBC[11]); 
  surft2_up[0] = -(0.408248290463863*GTC[4])+0.408248290463863*GCC[4]+0.3535533905932737*(GTC[0]+GCC[0]); 
  surft2_up[1] = -(0.408248290463863*GTC[8])+0.408248290463863*GCC[8]+0.3535533905932737*(GTC[1]+GCC[1]); 
  surft2_up[2] = -(0.408248290463863*GTC[9])+0.408248290463863*GCC[9]+0.3535533905932737*(GTC[2]+GCC[2]); 
  surft2_up[3] = -(0.408248290463863*GTC[10])+0.408248290463863*GCC[10]+0.3535533905932737*(GTC[3]+GCC[3]); 
  surft2_up[4] = -(0.408248290463863*GTC[12])+0.408248290463863*GCC[12]+0.3535533905932737*(GTC[5]+GCC[5]); 
  surft2_up[5] = -(0.408248290463863*GTC[13])+0.408248290463863*GCC[13]+0.3535533905932737*(GTC[6]+GCC[6]); 
  surft2_up[6] = -(0.408248290463863*GTC[14])+0.408248290463863*GCC[14]+0.3535533905932737*(GTC[7]+GCC[7]); 
  surft2_up[7] = -(0.408248290463863*GTC[15])+0.408248290463863*GCC[15]+0.3535533905932737*(GTC[11]+GCC[11]); 

  vol[9] = 3.0*GCC[0]; 
  vol[12] = 3.0*GCC[1]; 
  vol[14] = 3.0*GCC[3]; 
  vol[15] = 3.0*GCC[6]; 

  out[0] = (vol[0]+0.7071067811865475*surft1_up[0]-0.7071067811865475*surft1_lo[0])*dv1_sq*gamma_avg; 
  out[1] = (vol[1]+0.7071067811865475*surft1_up[1]-0.7071067811865475*surft1_lo[1])*dv1_sq*gamma_avg; 
  out[2] = (vol[2]-1.224744871391589*surft2_up[0]+1.224744871391589*(surft2_lo[0]+surft1_up[0]+surft1_lo[0]))*dv1_sq*gamma_avg; 
  out[3] = (vol[3]+0.7071067811865475*surft1_up[2]-0.7071067811865475*surft1_lo[2])*dv1_sq*gamma_avg; 
  out[4] = (vol[4]+0.7071067811865475*surft1_up[3]-0.7071067811865475*surft1_lo[3])*dv1_sq*gamma_avg; 
  out[5] = (vol[5]-1.224744871391589*surft2_up[1]+1.224744871391589*(surft2_lo[1]+surft1_up[1]+surft1_lo[1]))*dv1_sq*gamma_avg; 
  out[6] = (vol[6]+0.7071067811865475*surft1_up[4]-0.7071067811865475*surft1_lo[4])*dv1_sq*gamma_avg; 
  out[7] = (vol[7]-1.224744871391589*surft2_up[3]+1.224744871391589*(surft2_lo[3]+surft1_up[2]+surft1_lo[2]))*dv1_sq*gamma_avg; 
  out[8] = (vol[8]+0.7071067811865475*surft1_up[5]-0.7071067811865475*surft1_lo[5])*dv1_sq*gamma_avg; 
  out[9] = (vol[9]+1.224744871391589*(surft1_up[3]+surft1_lo[3])-2.1213203435596424*(surft2_up[0]+surft2_lo[0]))*dv1_sq*gamma_avg; 
  out[10] = (vol[10]+0.7071067811865475*surft1_up[6]-0.7071067811865475*surft1_lo[6])*dv1_sq*gamma_avg; 
  out[11] = (vol[11]-1.224744871391589*surft2_up[5]+1.224744871391589*(surft2_lo[5]+surft1_up[4]+surft1_lo[4]))*dv1_sq*gamma_avg; 
  out[12] = (vol[12]+1.224744871391589*(surft1_up[5]+surft1_lo[5])-2.1213203435596424*(surft2_up[1]+surft2_lo[1]))*dv1_sq*gamma_avg; 
  out[13] = (vol[13]+0.7071067811865475*surft1_up[7]-0.7071067811865475*surft1_lo[7])*dv1_sq*gamma_avg; 
  out[14] = (vol[14]+1.224744871391589*(surft1_up[6]+surft1_lo[6])-2.1213203435596424*(surft2_up[3]+surft2_lo[3]))*dv1_sq*gamma_avg; 
  out[15] = (vol[15]+1.224744871391589*(surft1_up[7]+surft1_lo[7])-2.1213203435596424*(surft2_up[5]+surft2_lo[5]))*dv1_sq*gamma_avg; 
} 
