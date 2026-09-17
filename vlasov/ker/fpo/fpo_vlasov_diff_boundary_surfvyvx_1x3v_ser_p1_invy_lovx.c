#include <gkyl_fpo_vlasov_kernels.h> 
 
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvyvx_1x3v_ser_p1_invy_lovx(const double *dxv, const double *diff_coeff_C,
      const double *diff_coeff_surf_stencil[9], const double *f_stencil[9], double* GKYL_RESTRICT out) { 
  // Stencil indices for this kernel: [-1, 0] 

  // dxv[NDIM]: Cell spacing in each direction. 
  // diff_coeff_stencil[3]: 3-cell stencil of diffusion tensor. 
  // f_stencil[9]: 9-cell stencil of distribution function. 
  // out: Incremented output. 


  double dv_inv_sq = 4.0/dxv[2]/dxv[1]; 
 
  double f_rec_lo[8] = {0.0}; 
  double f_rec_up[8] = {0.0}; 
  double df_rec_lo[8] = {0.0}; 
  double df_rec_up[8] = {0.0}; 
  double surft1_lo[8] = {0.0}; 
  double surft1_up[8] = {0.0}; 
  double surft2_lo[8] = {0.0}; 
  double surft2_up[8] = {0.0}; 
  double vol[16] = {0.0}; 

  // Index into D and f stencils. 
  const double *DCC = &diff_coeff_C[48]; 
  const double *fBL = f_stencil[0]; 
  const double *fCL = f_stencil[1]; 
  const double *fTL = f_stencil[2]; 
  const double *fBC = f_stencil[3]; 
  const double *fCC = f_stencil[4]; 
  const double *Dsurf_CC_vy = &diff_coeff_surf_stencil[4][48]; 
  const double *Dsurf_CC_vx = &diff_coeff_surf_stencil[4][56]; 
  const double *fTC = f_stencil[5]; 
  const double* Dsurf_TC_vy = &diff_coeff_surf_stencil[5][48]; 
  const double *fBR = f_stencil[6]; 
  const double *fCR = f_stencil[7]; 
  const double* Dsurf_CR_vx = &diff_coeff_surf_stencil[7][56]; 
  const double *fTR = f_stencil[8]; 

  f_rec_lo[0] = 0.7071067811865475*fCC[0]-1.224744871391589*fCC[2]; 
  f_rec_lo[1] = 0.7071067811865475*fCC[1]-1.224744871391589*fCC[5]; 
  f_rec_lo[2] = 0.7071067811865475*fCC[3]-1.224744871391589*fCC[7]; 
  f_rec_lo[3] = 0.7071067811865475*fCC[4]-1.224744871391589*fCC[9]; 
  f_rec_lo[4] = 0.7071067811865475*fCC[6]-1.224744871391589*fCC[11]; 
  f_rec_lo[5] = 0.7071067811865475*fCC[8]-1.224744871391589*fCC[12]; 
  f_rec_lo[6] = 0.7071067811865475*fCC[10]-1.224744871391589*fCC[14]; 
  f_rec_lo[7] = 0.7071067811865475*fCC[13]-1.224744871391589*fCC[15]; 
  f_rec_up[0] = -(0.408248290463863*fCR[2])+0.408248290463863*fCC[2]+0.3535533905932737*(fCR[0]+fCC[0]); 
  f_rec_up[1] = -(0.408248290463863*fCR[5])+0.408248290463863*fCC[5]+0.3535533905932737*(fCR[1]+fCC[1]); 
  f_rec_up[2] = -(0.408248290463863*fCR[7])+0.408248290463863*fCC[7]+0.3535533905932737*(fCR[3]+fCC[3]); 
  f_rec_up[3] = -(0.408248290463863*fCR[9])+0.408248290463863*fCC[9]+0.3535533905932737*(fCR[4]+fCC[4]); 
  f_rec_up[4] = -(0.408248290463863*fCR[11])+0.408248290463863*fCC[11]+0.3535533905932737*(fCR[6]+fCC[6]); 
  f_rec_up[5] = -(0.408248290463863*fCR[12])+0.408248290463863*fCC[12]+0.3535533905932737*(fCR[8]+fCC[8]); 
  f_rec_up[6] = -(0.408248290463863*fCR[14])+0.408248290463863*fCC[14]+0.3535533905932737*(fCR[10]+fCC[10]); 
  f_rec_up[7] = -(0.408248290463863*fCR[15])+0.408248290463863*fCC[15]+0.3535533905932737*(fCR[13]+fCC[13]); 

  df_rec_lo[0] = -(0.7071067811865475*fCC[7])+0.7071067811865475*fBC[7]+0.6123724356957944*(fCC[2]+fBC[2]); 
  df_rec_lo[1] = -(0.7071067811865475*fCC[11])+0.7071067811865475*fBC[11]+0.6123724356957944*(fCC[5]+fBC[5]); 
  df_rec_lo[3] = -(0.7071067811865475*fCC[14])+0.7071067811865475*fBC[14]+0.6123724356957944*(fCC[9]+fBC[9]); 
  df_rec_lo[5] = -(0.7071067811865475*fCC[15])+0.7071067811865475*fBC[15]+0.6123724356957944*(fCC[12]+fBC[12]); 
  df_rec_up[0] = -(0.7071067811865475*fTC[7])+0.7071067811865475*fCC[7]+0.6123724356957944*(fTC[2]+fCC[2]); 
  df_rec_up[1] = -(0.7071067811865475*fTC[11])+0.7071067811865475*fCC[11]+0.6123724356957944*(fTC[5]+fCC[5]); 
  df_rec_up[3] = -(0.7071067811865475*fTC[14])+0.7071067811865475*fCC[14]+0.6123724356957944*(fTC[9]+fCC[9]); 
  df_rec_up[5] = -(0.7071067811865475*fTC[15])+0.7071067811865475*fCC[15]+0.6123724356957944*(fTC[12]+fCC[12]); 

  surft1_lo[0] = 0.3535533905932737*(Dsurf_CC_vy[7]*df_rec_lo[7]+Dsurf_CC_vy[6]*df_rec_lo[6]+Dsurf_CC_vy[5]*df_rec_lo[5]+Dsurf_CC_vy[4]*df_rec_lo[4]+Dsurf_CC_vy[3]*df_rec_lo[3]+Dsurf_CC_vy[2]*df_rec_lo[2]+Dsurf_CC_vy[1]*df_rec_lo[1]+Dsurf_CC_vy[0]*df_rec_lo[0]); 
  surft1_lo[1] = 0.3535533905932737*(Dsurf_CC_vy[6]*df_rec_lo[7]+df_rec_lo[6]*Dsurf_CC_vy[7]+Dsurf_CC_vy[3]*df_rec_lo[5]+df_rec_lo[3]*Dsurf_CC_vy[5]+Dsurf_CC_vy[2]*df_rec_lo[4]+df_rec_lo[2]*Dsurf_CC_vy[4]+Dsurf_CC_vy[0]*df_rec_lo[1]+df_rec_lo[0]*Dsurf_CC_vy[1]); 
  surft1_lo[2] = 0.3535533905932737*(Dsurf_CC_vy[5]*df_rec_lo[7]+df_rec_lo[5]*Dsurf_CC_vy[7]+Dsurf_CC_vy[3]*df_rec_lo[6]+df_rec_lo[3]*Dsurf_CC_vy[6]+Dsurf_CC_vy[1]*df_rec_lo[4]+df_rec_lo[1]*Dsurf_CC_vy[4]+Dsurf_CC_vy[0]*df_rec_lo[2]+df_rec_lo[0]*Dsurf_CC_vy[2]); 
  surft1_lo[3] = 0.3535533905932737*(Dsurf_CC_vy[4]*df_rec_lo[7]+df_rec_lo[4]*Dsurf_CC_vy[7]+Dsurf_CC_vy[2]*df_rec_lo[6]+df_rec_lo[2]*Dsurf_CC_vy[6]+Dsurf_CC_vy[1]*df_rec_lo[5]+df_rec_lo[1]*Dsurf_CC_vy[5]+Dsurf_CC_vy[0]*df_rec_lo[3]+df_rec_lo[0]*Dsurf_CC_vy[3]); 
  surft1_lo[4] = 0.3535533905932737*(Dsurf_CC_vy[3]*df_rec_lo[7]+df_rec_lo[3]*Dsurf_CC_vy[7]+Dsurf_CC_vy[5]*df_rec_lo[6]+df_rec_lo[5]*Dsurf_CC_vy[6]+Dsurf_CC_vy[0]*df_rec_lo[4]+df_rec_lo[0]*Dsurf_CC_vy[4]+Dsurf_CC_vy[1]*df_rec_lo[2]+df_rec_lo[1]*Dsurf_CC_vy[2]); 
  surft1_lo[5] = 0.3535533905932737*(Dsurf_CC_vy[2]*df_rec_lo[7]+df_rec_lo[2]*Dsurf_CC_vy[7]+Dsurf_CC_vy[4]*df_rec_lo[6]+df_rec_lo[4]*Dsurf_CC_vy[6]+Dsurf_CC_vy[0]*df_rec_lo[5]+df_rec_lo[0]*Dsurf_CC_vy[5]+Dsurf_CC_vy[1]*df_rec_lo[3]+df_rec_lo[1]*Dsurf_CC_vy[3]); 
  surft1_lo[6] = 0.3535533905932737*(Dsurf_CC_vy[1]*df_rec_lo[7]+df_rec_lo[1]*Dsurf_CC_vy[7]+Dsurf_CC_vy[0]*df_rec_lo[6]+df_rec_lo[0]*Dsurf_CC_vy[6]+Dsurf_CC_vy[4]*df_rec_lo[5]+df_rec_lo[4]*Dsurf_CC_vy[5]+Dsurf_CC_vy[2]*df_rec_lo[3]+df_rec_lo[2]*Dsurf_CC_vy[3]); 
  surft1_lo[7] = 0.3535533905932737*(Dsurf_CC_vy[0]*df_rec_lo[7]+df_rec_lo[0]*Dsurf_CC_vy[7]+Dsurf_CC_vy[1]*df_rec_lo[6]+df_rec_lo[1]*Dsurf_CC_vy[6]+Dsurf_CC_vy[2]*df_rec_lo[5]+df_rec_lo[2]*Dsurf_CC_vy[5]+Dsurf_CC_vy[3]*df_rec_lo[4]+df_rec_lo[3]*Dsurf_CC_vy[4]); 
  surft1_up[0] = 0.3535533905932737*(Dsurf_TC_vy[7]*df_rec_up[7]+Dsurf_TC_vy[6]*df_rec_up[6]+Dsurf_TC_vy[5]*df_rec_up[5]+Dsurf_TC_vy[4]*df_rec_up[4]+Dsurf_TC_vy[3]*df_rec_up[3]+Dsurf_TC_vy[2]*df_rec_up[2]+Dsurf_TC_vy[1]*df_rec_up[1]+Dsurf_TC_vy[0]*df_rec_up[0]); 
  surft1_up[1] = 0.3535533905932737*(Dsurf_TC_vy[6]*df_rec_up[7]+df_rec_up[6]*Dsurf_TC_vy[7]+Dsurf_TC_vy[3]*df_rec_up[5]+df_rec_up[3]*Dsurf_TC_vy[5]+Dsurf_TC_vy[2]*df_rec_up[4]+df_rec_up[2]*Dsurf_TC_vy[4]+Dsurf_TC_vy[0]*df_rec_up[1]+df_rec_up[0]*Dsurf_TC_vy[1]); 
  surft1_up[2] = 0.3535533905932737*(Dsurf_TC_vy[5]*df_rec_up[7]+df_rec_up[5]*Dsurf_TC_vy[7]+Dsurf_TC_vy[3]*df_rec_up[6]+df_rec_up[3]*Dsurf_TC_vy[6]+Dsurf_TC_vy[1]*df_rec_up[4]+df_rec_up[1]*Dsurf_TC_vy[4]+Dsurf_TC_vy[0]*df_rec_up[2]+df_rec_up[0]*Dsurf_TC_vy[2]); 
  surft1_up[3] = 0.3535533905932737*(Dsurf_TC_vy[4]*df_rec_up[7]+df_rec_up[4]*Dsurf_TC_vy[7]+Dsurf_TC_vy[2]*df_rec_up[6]+df_rec_up[2]*Dsurf_TC_vy[6]+Dsurf_TC_vy[1]*df_rec_up[5]+df_rec_up[1]*Dsurf_TC_vy[5]+Dsurf_TC_vy[0]*df_rec_up[3]+df_rec_up[0]*Dsurf_TC_vy[3]); 
  surft1_up[4] = 0.3535533905932737*(Dsurf_TC_vy[3]*df_rec_up[7]+df_rec_up[3]*Dsurf_TC_vy[7]+Dsurf_TC_vy[5]*df_rec_up[6]+df_rec_up[5]*Dsurf_TC_vy[6]+Dsurf_TC_vy[0]*df_rec_up[4]+df_rec_up[0]*Dsurf_TC_vy[4]+Dsurf_TC_vy[1]*df_rec_up[2]+df_rec_up[1]*Dsurf_TC_vy[2]); 
  surft1_up[5] = 0.3535533905932737*(Dsurf_TC_vy[2]*df_rec_up[7]+df_rec_up[2]*Dsurf_TC_vy[7]+Dsurf_TC_vy[4]*df_rec_up[6]+df_rec_up[4]*Dsurf_TC_vy[6]+Dsurf_TC_vy[0]*df_rec_up[5]+df_rec_up[0]*Dsurf_TC_vy[5]+Dsurf_TC_vy[1]*df_rec_up[3]+df_rec_up[1]*Dsurf_TC_vy[3]); 
  surft1_up[6] = 0.3535533905932737*(Dsurf_TC_vy[1]*df_rec_up[7]+df_rec_up[1]*Dsurf_TC_vy[7]+Dsurf_TC_vy[0]*df_rec_up[6]+df_rec_up[0]*Dsurf_TC_vy[6]+Dsurf_TC_vy[4]*df_rec_up[5]+df_rec_up[4]*Dsurf_TC_vy[5]+Dsurf_TC_vy[2]*df_rec_up[3]+df_rec_up[2]*Dsurf_TC_vy[3]); 
  surft1_up[7] = 0.3535533905932737*(Dsurf_TC_vy[0]*df_rec_up[7]+df_rec_up[0]*Dsurf_TC_vy[7]+Dsurf_TC_vy[1]*df_rec_up[6]+df_rec_up[1]*Dsurf_TC_vy[6]+Dsurf_TC_vy[2]*df_rec_up[5]+df_rec_up[2]*Dsurf_TC_vy[5]+Dsurf_TC_vy[3]*df_rec_up[4]+df_rec_up[3]*Dsurf_TC_vy[4]); 

  surft2_lo[0] = 0.25*(f_rec_lo[7]*DCC[13]+f_rec_lo[6]*DCC[10]+f_rec_lo[5]*DCC[8]+f_rec_lo[4]*DCC[6]+f_rec_lo[3]*DCC[4]+f_rec_lo[2]*DCC[3]+DCC[1]*f_rec_lo[1]+DCC[0]*f_rec_lo[0])-0.4330127018922193*(f_rec_lo[7]*DCC[15]+f_rec_lo[6]*DCC[14]+f_rec_lo[5]*DCC[12]+f_rec_lo[4]*DCC[11]+f_rec_lo[3]*DCC[9]+f_rec_lo[2]*DCC[7]+f_rec_lo[1]*DCC[5]+f_rec_lo[0]*DCC[2]); 
  surft2_lo[1] = 0.25*(f_rec_lo[6]*DCC[13]+f_rec_lo[7]*DCC[10]+f_rec_lo[3]*DCC[8]+f_rec_lo[2]*DCC[6]+DCC[4]*f_rec_lo[5]+DCC[3]*f_rec_lo[4]+DCC[0]*f_rec_lo[1]+f_rec_lo[0]*DCC[1])-0.4330127018922193*(f_rec_lo[6]*DCC[15]+f_rec_lo[7]*DCC[14]+f_rec_lo[3]*DCC[12]+f_rec_lo[2]*DCC[11]+f_rec_lo[5]*DCC[9]+f_rec_lo[4]*DCC[7]+f_rec_lo[0]*DCC[5]+f_rec_lo[1]*DCC[2]); 
  surft2_lo[2] = 0.25*(f_rec_lo[5]*DCC[13]+f_rec_lo[3]*DCC[10]+f_rec_lo[7]*DCC[8]+DCC[4]*f_rec_lo[6]+f_rec_lo[1]*DCC[6]+DCC[1]*f_rec_lo[4]+f_rec_lo[0]*DCC[3]+DCC[0]*f_rec_lo[2])-0.4330127018922193*(f_rec_lo[5]*DCC[15]+f_rec_lo[3]*DCC[14]+f_rec_lo[7]*DCC[12]+f_rec_lo[1]*DCC[11]+f_rec_lo[6]*DCC[9]+f_rec_lo[0]*DCC[7]+f_rec_lo[4]*DCC[5]+DCC[2]*f_rec_lo[2]); 
  surft2_lo[3] = 0.25*(f_rec_lo[4]*DCC[13]+f_rec_lo[2]*DCC[10]+f_rec_lo[1]*DCC[8]+DCC[6]*f_rec_lo[7]+DCC[3]*f_rec_lo[6]+DCC[1]*f_rec_lo[5]+f_rec_lo[0]*DCC[4]+DCC[0]*f_rec_lo[3])-0.4330127018922193*(f_rec_lo[4]*DCC[15]+f_rec_lo[2]*DCC[14]+f_rec_lo[1]*DCC[12]+f_rec_lo[7]*DCC[11]+f_rec_lo[0]*DCC[9]+f_rec_lo[6]*DCC[7]+DCC[5]*f_rec_lo[5]+DCC[2]*f_rec_lo[3]); 
  surft2_lo[4] = 0.25*(f_rec_lo[3]*DCC[13]+f_rec_lo[5]*DCC[10]+f_rec_lo[6]*DCC[8]+DCC[4]*f_rec_lo[7]+f_rec_lo[0]*DCC[6]+DCC[0]*f_rec_lo[4]+f_rec_lo[1]*DCC[3]+DCC[1]*f_rec_lo[2])-0.4330127018922193*(f_rec_lo[3]*DCC[15]+f_rec_lo[5]*DCC[14]+f_rec_lo[6]*DCC[12]+f_rec_lo[0]*DCC[11]+f_rec_lo[7]*DCC[9]+f_rec_lo[1]*DCC[7]+f_rec_lo[2]*DCC[5]+DCC[2]*f_rec_lo[4]); 
  surft2_lo[5] = 0.25*(f_rec_lo[2]*DCC[13]+f_rec_lo[4]*DCC[10]+f_rec_lo[0]*DCC[8]+DCC[3]*f_rec_lo[7]+DCC[6]*f_rec_lo[6]+DCC[0]*f_rec_lo[5]+f_rec_lo[1]*DCC[4]+DCC[1]*f_rec_lo[3])-0.4330127018922193*(f_rec_lo[2]*DCC[15]+f_rec_lo[4]*DCC[14]+f_rec_lo[0]*DCC[12]+f_rec_lo[6]*DCC[11]+f_rec_lo[1]*DCC[9]+DCC[7]*f_rec_lo[7]+DCC[2]*f_rec_lo[5]+f_rec_lo[3]*DCC[5]); 
  surft2_lo[6] = 0.25*(f_rec_lo[1]*DCC[13]+f_rec_lo[0]*DCC[10]+f_rec_lo[4]*DCC[8]+DCC[1]*f_rec_lo[7]+DCC[0]*f_rec_lo[6]+f_rec_lo[5]*DCC[6]+f_rec_lo[2]*DCC[4]+DCC[3]*f_rec_lo[3])-0.4330127018922193*(f_rec_lo[1]*DCC[15]+f_rec_lo[0]*DCC[14]+f_rec_lo[4]*DCC[12]+f_rec_lo[5]*DCC[11]+f_rec_lo[2]*DCC[9]+DCC[5]*f_rec_lo[7]+f_rec_lo[3]*DCC[7]+DCC[2]*f_rec_lo[6]); 
  surft2_lo[7] = 0.25*(f_rec_lo[0]*DCC[13]+f_rec_lo[1]*DCC[10]+f_rec_lo[2]*DCC[8]+DCC[0]*f_rec_lo[7]+DCC[1]*f_rec_lo[6]+f_rec_lo[3]*DCC[6]+DCC[3]*f_rec_lo[5]+DCC[4]*f_rec_lo[4])-0.4330127018922193*(f_rec_lo[0]*DCC[15]+f_rec_lo[1]*DCC[14]+f_rec_lo[2]*DCC[12]+f_rec_lo[3]*DCC[11]+f_rec_lo[4]*DCC[9]+DCC[2]*f_rec_lo[7]+f_rec_lo[5]*DCC[7]+DCC[5]*f_rec_lo[6]); 
  surft2_up[0] = 0.3535533905932737*(Dsurf_CR_vx[7]*f_rec_up[7]+Dsurf_CR_vx[6]*f_rec_up[6]+Dsurf_CR_vx[5]*f_rec_up[5]+Dsurf_CR_vx[4]*f_rec_up[4]+Dsurf_CR_vx[3]*f_rec_up[3]+Dsurf_CR_vx[2]*f_rec_up[2]+Dsurf_CR_vx[1]*f_rec_up[1]+Dsurf_CR_vx[0]*f_rec_up[0]); 
  surft2_up[1] = 0.3535533905932737*(Dsurf_CR_vx[6]*f_rec_up[7]+f_rec_up[6]*Dsurf_CR_vx[7]+Dsurf_CR_vx[3]*f_rec_up[5]+f_rec_up[3]*Dsurf_CR_vx[5]+Dsurf_CR_vx[2]*f_rec_up[4]+f_rec_up[2]*Dsurf_CR_vx[4]+Dsurf_CR_vx[0]*f_rec_up[1]+f_rec_up[0]*Dsurf_CR_vx[1]); 
  surft2_up[2] = 0.3535533905932737*(Dsurf_CR_vx[5]*f_rec_up[7]+f_rec_up[5]*Dsurf_CR_vx[7]+Dsurf_CR_vx[3]*f_rec_up[6]+f_rec_up[3]*Dsurf_CR_vx[6]+Dsurf_CR_vx[1]*f_rec_up[4]+f_rec_up[1]*Dsurf_CR_vx[4]+Dsurf_CR_vx[0]*f_rec_up[2]+f_rec_up[0]*Dsurf_CR_vx[2]); 
  surft2_up[3] = 0.3535533905932737*(Dsurf_CR_vx[4]*f_rec_up[7]+f_rec_up[4]*Dsurf_CR_vx[7]+Dsurf_CR_vx[2]*f_rec_up[6]+f_rec_up[2]*Dsurf_CR_vx[6]+Dsurf_CR_vx[1]*f_rec_up[5]+f_rec_up[1]*Dsurf_CR_vx[5]+Dsurf_CR_vx[0]*f_rec_up[3]+f_rec_up[0]*Dsurf_CR_vx[3]); 
  surft2_up[4] = 0.3535533905932737*(Dsurf_CR_vx[3]*f_rec_up[7]+f_rec_up[3]*Dsurf_CR_vx[7]+Dsurf_CR_vx[5]*f_rec_up[6]+f_rec_up[5]*Dsurf_CR_vx[6]+Dsurf_CR_vx[0]*f_rec_up[4]+f_rec_up[0]*Dsurf_CR_vx[4]+Dsurf_CR_vx[1]*f_rec_up[2]+f_rec_up[1]*Dsurf_CR_vx[2]); 
  surft2_up[5] = 0.3535533905932737*(Dsurf_CR_vx[2]*f_rec_up[7]+f_rec_up[2]*Dsurf_CR_vx[7]+Dsurf_CR_vx[4]*f_rec_up[6]+f_rec_up[4]*Dsurf_CR_vx[6]+Dsurf_CR_vx[0]*f_rec_up[5]+f_rec_up[0]*Dsurf_CR_vx[5]+Dsurf_CR_vx[1]*f_rec_up[3]+f_rec_up[1]*Dsurf_CR_vx[3]); 
  surft2_up[6] = 0.3535533905932737*(Dsurf_CR_vx[1]*f_rec_up[7]+f_rec_up[1]*Dsurf_CR_vx[7]+Dsurf_CR_vx[0]*f_rec_up[6]+f_rec_up[0]*Dsurf_CR_vx[6]+Dsurf_CR_vx[4]*f_rec_up[5]+f_rec_up[4]*Dsurf_CR_vx[5]+Dsurf_CR_vx[2]*f_rec_up[3]+f_rec_up[2]*Dsurf_CR_vx[3]); 
  surft2_up[7] = 0.3535533905932737*(Dsurf_CR_vx[0]*f_rec_up[7]+f_rec_up[0]*Dsurf_CR_vx[7]+Dsurf_CR_vx[1]*f_rec_up[6]+f_rec_up[1]*Dsurf_CR_vx[6]+Dsurf_CR_vx[2]*f_rec_up[5]+f_rec_up[2]*Dsurf_CR_vx[5]+Dsurf_CR_vx[3]*f_rec_up[4]+f_rec_up[3]*Dsurf_CR_vx[4]); 

  vol[2] = 0.75*(fCC[12]*DCC[15]+fCC[9]*DCC[14]+fCC[8]*DCC[13]+fCC[5]*DCC[11]+fCC[4]*DCC[10]+fCC[2]*DCC[7]+fCC[1]*DCC[6]+fCC[0]*DCC[3]); 
  vol[5] = 0.75*(fCC[9]*DCC[15]+fCC[12]*DCC[14]+fCC[4]*DCC[13]+fCC[2]*DCC[11]+fCC[8]*DCC[10]+fCC[5]*DCC[7]+fCC[0]*DCC[6]+fCC[1]*DCC[3]); 
  vol[7] = 1.5*(DCC[15]*fCC[15]+DCC[14]*fCC[14]+DCC[13]*fCC[13]+DCC[11]*fCC[11]+DCC[10]*fCC[10]+DCC[7]*fCC[7]+DCC[6]*fCC[6]+DCC[3]*fCC[3])+0.75*(DCC[12]*fCC[12]+DCC[9]*fCC[9]+DCC[8]*fCC[8]+DCC[5]*fCC[5]+DCC[4]*fCC[4]+DCC[2]*fCC[2]+DCC[1]*fCC[1]+DCC[0]*fCC[0]); 
  vol[9] = 0.75*(fCC[5]*DCC[15]+fCC[2]*DCC[14]+fCC[1]*DCC[13]+DCC[11]*fCC[12]+fCC[0]*DCC[10]+DCC[7]*fCC[9]+DCC[6]*fCC[8]+DCC[3]*fCC[4]); 
  vol[11] = 1.5*(DCC[14]*fCC[15]+fCC[14]*DCC[15]+DCC[10]*fCC[13]+fCC[10]*DCC[13]+DCC[7]*fCC[11]+fCC[7]*DCC[11]+DCC[3]*fCC[6]+fCC[3]*DCC[6])+0.75*(DCC[9]*fCC[12]+fCC[9]*DCC[12]+DCC[4]*fCC[8]+fCC[4]*DCC[8]+DCC[2]*fCC[5]+fCC[2]*DCC[5]+DCC[0]*fCC[1]+fCC[0]*DCC[1]); 
  vol[12] = 0.75*(fCC[2]*DCC[15]+fCC[5]*DCC[14]+fCC[0]*DCC[13]+DCC[7]*fCC[12]+fCC[9]*DCC[11]+fCC[1]*DCC[10]+DCC[3]*fCC[8]+fCC[4]*DCC[6]); 
  vol[14] = 1.5*(DCC[11]*fCC[15]+fCC[11]*DCC[15]+DCC[7]*fCC[14]+fCC[7]*DCC[14]+DCC[6]*fCC[13]+fCC[6]*DCC[13]+DCC[3]*fCC[10]+fCC[3]*DCC[10])+0.75*(DCC[5]*fCC[12]+fCC[5]*DCC[12]+DCC[2]*fCC[9]+fCC[2]*DCC[9]+DCC[1]*fCC[8]+fCC[1]*DCC[8]+DCC[0]*fCC[4]+fCC[0]*DCC[4]); 
  vol[15] = 1.5*(DCC[7]*fCC[15]+fCC[7]*DCC[15]+DCC[11]*fCC[14]+fCC[11]*DCC[14]+DCC[3]*fCC[13]+fCC[3]*DCC[13]+DCC[6]*fCC[10]+fCC[6]*DCC[10])+0.75*(DCC[2]*fCC[12]+fCC[2]*DCC[12]+DCC[5]*fCC[9]+fCC[5]*DCC[9]+DCC[0]*fCC[8]+fCC[0]*DCC[8]+DCC[1]*fCC[4]+fCC[1]*DCC[4]); 

  out[0] += (0.5*vol[0]+0.35355339059327373*surft1_up[0]-0.35355339059327373*surft1_lo[0])*dv_inv_sq; 
  out[1] += (0.5*vol[1]+0.35355339059327373*surft1_up[1]-0.35355339059327373*surft1_lo[1])*dv_inv_sq; 
  out[2] += (0.5*vol[2]+0.35355339059327373*surft1_up[2]-0.35355339059327373*surft1_lo[2])*dv_inv_sq; 
  out[3] += (0.5*vol[3]-0.6123724356957945*surft2_up[0]+0.6123724356957945*(surft2_lo[0]+surft1_up[0]+surft1_lo[0]))*dv_inv_sq; 
  out[4] += (0.5*vol[4]+0.35355339059327373*surft1_up[3]-0.35355339059327373*surft1_lo[3])*dv_inv_sq; 
  out[5] += (0.5*vol[5]+0.35355339059327373*surft1_up[4]-0.35355339059327373*surft1_lo[4])*dv_inv_sq; 
  out[6] += (0.5*vol[6]-0.6123724356957945*surft2_up[1]+0.6123724356957945*(surft2_lo[1]+surft1_up[1]+surft1_lo[1]))*dv_inv_sq; 
  out[7] += (0.5*vol[7]+0.6123724356957945*(surft1_up[2]+surft1_lo[2])-1.0606601717798212*(surft2_up[0]+surft2_lo[0]))*dv_inv_sq; 
  out[8] += (0.5*vol[8]+0.35355339059327373*surft1_up[5]-0.35355339059327373*surft1_lo[5])*dv_inv_sq; 
  out[9] += (0.5*vol[9]+0.35355339059327373*surft1_up[6]-0.35355339059327373*surft1_lo[6])*dv_inv_sq; 
  out[10] += (0.5*vol[10]-0.6123724356957945*surft2_up[3]+0.6123724356957945*(surft2_lo[3]+surft1_up[3]+surft1_lo[3]))*dv_inv_sq; 
  out[11] += (0.5*vol[11]+0.6123724356957945*(surft1_up[4]+surft1_lo[4])-1.0606601717798212*(surft2_up[1]+surft2_lo[1]))*dv_inv_sq; 
  out[12] += (0.5*vol[12]+0.35355339059327373*surft1_up[7]-0.35355339059327373*surft1_lo[7])*dv_inv_sq; 
  out[13] += (0.5*vol[13]-0.6123724356957945*surft2_up[5]+0.6123724356957945*(surft2_lo[5]+surft1_up[5]+surft1_lo[5]))*dv_inv_sq; 
  out[14] += (0.5*vol[14]+0.6123724356957945*(surft1_up[6]+surft1_lo[6])-1.0606601717798212*(surft2_up[3]+surft2_lo[3]))*dv_inv_sq; 
  out[15] += (0.5*vol[15]+0.6123724356957945*(surft1_up[7]+surft1_lo[7])-1.0606601717798212*(surft2_up[5]+surft2_lo[5]))*dv_inv_sq; 
  double cflFreq = fmax(fabs(Dsurf_CC_vy[0]), fabs(Dsurf_TC_vy[0])); 

  return 1.4142135623730947*dv_inv_sq*cflFreq; 
} 
