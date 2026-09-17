#include <gkyl_fpo_vlasov_kernels.h> 
 
  // Stencil indices for this kernel: [-1, 0] 

GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvz_1x3v_ser_p1_lovx_invz(const double *dxv, const double *diff_coeff_C,
      const double *diff_coeff_surf_stencil[9], const double *f_stencil[9], double* GKYL_RESTRICT out) { 
  // dxv[NDIM]: Cell spacing in each direction. 
  // diff_coeff_stencil[3]: 3-cell stencil of diffusion tensor. 
  // f_stencil[9]: 9-cell stencil of distribution function. 
  // out: Incremented output. 


  double dv_inv_sq = 4.0/dxv[1]/dxv[3]; 
 
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
  const double *DCC = &diff_coeff_C[32]; 
  const double *fBL = f_stencil[0]; 
  const double *fCL = f_stencil[1]; 
  const double *fTL = f_stencil[2]; 
  const double *fBC = f_stencil[3]; 
  const double *fCC = f_stencil[4]; 
  const double *Dsurf_CC_vx = &diff_coeff_surf_stencil[4][32]; 
  const double *Dsurf_CC_vz = &diff_coeff_surf_stencil[4][40]; 
  const double *fTC = f_stencil[5]; 
  const double* Dsurf_TC_vz = &diff_coeff_surf_stencil[5][40]; 
  const double *fBR = f_stencil[6]; 
  const double *fCR = f_stencil[7]; 
  const double* Dsurf_CR_vx = &diff_coeff_surf_stencil[7][32]; 
  const double *fTR = f_stencil[8]; 

  f_rec_lo[0] = -(0.408248290463863*fCC[4])+0.408248290463863*fBC[4]+0.3535533905932737*(fCC[0]+fBC[0]); 
  f_rec_lo[1] = -(0.408248290463863*fCC[8])+0.408248290463863*fBC[8]+0.3535533905932737*(fCC[1]+fBC[1]); 
  f_rec_lo[2] = -(0.408248290463863*fCC[9])+0.408248290463863*fBC[9]+0.3535533905932737*(fCC[2]+fBC[2]); 
  f_rec_lo[3] = -(0.408248290463863*fCC[10])+0.408248290463863*fBC[10]+0.3535533905932737*(fCC[3]+fBC[3]); 
  f_rec_lo[4] = -(0.408248290463863*fCC[12])+0.408248290463863*fBC[12]+0.3535533905932737*(fCC[5]+fBC[5]); 
  f_rec_lo[5] = -(0.408248290463863*fCC[13])+0.408248290463863*fBC[13]+0.3535533905932737*(fCC[6]+fBC[6]); 
  f_rec_lo[6] = -(0.408248290463863*fCC[14])+0.408248290463863*fBC[14]+0.3535533905932737*(fCC[7]+fBC[7]); 
  f_rec_lo[7] = -(0.408248290463863*fCC[15])+0.408248290463863*fBC[15]+0.3535533905932737*(fCC[11]+fBC[11]); 
  f_rec_up[0] = -(0.408248290463863*fTC[4])+0.408248290463863*fCC[4]+0.3535533905932737*(fTC[0]+fCC[0]); 
  f_rec_up[1] = -(0.408248290463863*fTC[8])+0.408248290463863*fCC[8]+0.3535533905932737*(fTC[1]+fCC[1]); 
  f_rec_up[2] = -(0.408248290463863*fTC[9])+0.408248290463863*fCC[9]+0.3535533905932737*(fTC[2]+fCC[2]); 
  f_rec_up[3] = -(0.408248290463863*fTC[10])+0.408248290463863*fCC[10]+0.3535533905932737*(fTC[3]+fCC[3]); 
  f_rec_up[4] = -(0.408248290463863*fTC[12])+0.408248290463863*fCC[12]+0.3535533905932737*(fTC[5]+fCC[5]); 
  f_rec_up[5] = -(0.408248290463863*fTC[13])+0.408248290463863*fCC[13]+0.3535533905932737*(fTC[6]+fCC[6]); 
  f_rec_up[6] = -(0.408248290463863*fTC[14])+0.408248290463863*fCC[14]+0.3535533905932737*(fTC[7]+fCC[7]); 
  f_rec_up[7] = -(0.408248290463863*fTC[15])+0.408248290463863*fCC[15]+0.3535533905932737*(fTC[11]+fCC[11]); 

  df_rec_up[0] = 0.11785113019775789*(fTR[9]+fBR[9])-0.11785113019775789*(fTC[9]+fBC[9])-0.2357022603955158*fCR[9]+0.2357022603955158*fCC[9]-0.10206207261596573*(fTR[4]+fTC[4]+fBR[4]+fBC[4]+fTR[2]+fBC[2])+0.20412414523193148*(fCR[4]+fCC[4])+0.10206207261596573*(fTC[2]+fBR[2])+0.0883883476483184*(fTR[0]+fTC[0])-0.0883883476483184*(fBR[0]+fBC[0]); 
  df_rec_up[1] = 0.11785113019775789*(fTR[12]+fBR[12])-0.11785113019775789*(fTC[12]+fBC[12])-0.2357022603955158*fCR[12]+0.2357022603955158*fCC[12]-0.10206207261596573*(fTR[8]+fTC[8]+fBR[8]+fBC[8]+fTR[5]+fBC[5])+0.20412414523193148*(fCR[8]+fCC[8])+0.10206207261596573*(fTC[5]+fBR[5])+0.0883883476483184*(fTR[1]+fTC[1])-0.0883883476483184*(fBR[1]+fBC[1]); 
  df_rec_up[2] = 0.11785113019775789*(fTR[14]+fBR[14])-0.11785113019775789*(fTC[14]+fBC[14])-0.2357022603955158*fCR[14]+0.2357022603955158*fCC[14]-0.10206207261596573*(fTR[10]+fTC[10]+fBR[10]+fBC[10]+fTR[7]+fBC[7])+0.20412414523193148*(fCR[10]+fCC[10])+0.10206207261596573*(fTC[7]+fBR[7])+0.0883883476483184*(fTR[3]+fTC[3])-0.0883883476483184*(fBR[3]+fBC[3]); 
  df_rec_up[3] = 0.20412414523193148*(fTR[9]+fBC[9])-0.20412414523193148*(fTC[9]+fBR[9])-0.1767766952966368*(fTR[4]+fTC[4]+fTR[2]+fBR[2])+0.1767766952966368*(fBR[4]+fBC[4]+fTC[2]+fBC[2])+0.3535533905932737*fCR[2]-0.3535533905932737*fCC[2]+0.15309310892394856*(fTR[0]+fTC[0]+fBR[0]+fBC[0])-0.3061862178478971*(fCR[0]+fCC[0]); 
  df_rec_up[4] = 0.11785113019775789*(fTR[15]+fBR[15])-0.11785113019775789*(fTC[15]+fBC[15])-0.2357022603955158*fCR[15]+0.2357022603955158*fCC[15]-0.10206207261596573*(fTR[13]+fTC[13]+fBR[13]+fBC[13]+fTR[11]+fBC[11])+0.20412414523193148*(fCR[13]+fCC[13])+0.10206207261596573*(fTC[11]+fBR[11])+0.0883883476483184*(fTR[6]+fTC[6])-0.0883883476483184*(fBR[6]+fBC[6]); 
  df_rec_up[5] = 0.20412414523193148*(fTR[12]+fBC[12])-0.20412414523193148*(fTC[12]+fBR[12])-0.1767766952966368*(fTR[8]+fTC[8]+fTR[5]+fBR[5])+0.1767766952966368*(fBR[8]+fBC[8]+fTC[5]+fBC[5])+0.3535533905932737*fCR[5]-0.3535533905932737*fCC[5]+0.15309310892394856*(fTR[1]+fTC[1]+fBR[1]+fBC[1])-0.3061862178478971*(fCR[1]+fCC[1]); 
  df_rec_up[6] = 0.20412414523193148*(fTR[14]+fBC[14])-0.20412414523193148*(fTC[14]+fBR[14])-0.1767766952966368*(fTR[10]+fTC[10]+fTR[7]+fBR[7])+0.1767766952966368*(fBR[10]+fBC[10]+fTC[7]+fBC[7])+0.3535533905932737*fCR[7]-0.3535533905932737*fCC[7]+0.15309310892394856*(fTR[3]+fTC[3]+fBR[3]+fBC[3])-0.3061862178478971*(fCR[3]+fCC[3]); 
  df_rec_up[7] = 0.20412414523193148*(fTR[15]+fBC[15])-0.20412414523193148*(fTC[15]+fBR[15])-0.1767766952966368*(fTR[13]+fTC[13]+fTR[11]+fBR[11])+0.1767766952966368*(fBR[13]+fBC[13]+fTC[11]+fBC[11])+0.3535533905932737*fCR[11]-0.3535533905932737*fCC[11]+0.15309310892394856*(fTR[6]+fTC[6]+fBR[6]+fBC[6])-0.3061862178478971*(fCR[6]+fCC[6]); 

  surft1_lo[0] = 0.3535533905932737*(Dsurf_CC_vx[7]*df_rec_lo[7]+Dsurf_CC_vx[6]*df_rec_lo[6]+Dsurf_CC_vx[5]*df_rec_lo[5]+Dsurf_CC_vx[4]*df_rec_lo[4]+Dsurf_CC_vx[3]*df_rec_lo[3]+Dsurf_CC_vx[2]*df_rec_lo[2]+Dsurf_CC_vx[1]*df_rec_lo[1]+Dsurf_CC_vx[0]*df_rec_lo[0]); 
  surft1_lo[1] = 0.3535533905932737*(Dsurf_CC_vx[6]*df_rec_lo[7]+df_rec_lo[6]*Dsurf_CC_vx[7]+Dsurf_CC_vx[3]*df_rec_lo[5]+df_rec_lo[3]*Dsurf_CC_vx[5]+Dsurf_CC_vx[2]*df_rec_lo[4]+df_rec_lo[2]*Dsurf_CC_vx[4]+Dsurf_CC_vx[0]*df_rec_lo[1]+df_rec_lo[0]*Dsurf_CC_vx[1]); 
  surft1_lo[2] = 0.3535533905932737*(Dsurf_CC_vx[5]*df_rec_lo[7]+df_rec_lo[5]*Dsurf_CC_vx[7]+Dsurf_CC_vx[3]*df_rec_lo[6]+df_rec_lo[3]*Dsurf_CC_vx[6]+Dsurf_CC_vx[1]*df_rec_lo[4]+df_rec_lo[1]*Dsurf_CC_vx[4]+Dsurf_CC_vx[0]*df_rec_lo[2]+df_rec_lo[0]*Dsurf_CC_vx[2]); 
  surft1_lo[3] = 0.3535533905932737*(Dsurf_CC_vx[4]*df_rec_lo[7]+df_rec_lo[4]*Dsurf_CC_vx[7]+Dsurf_CC_vx[2]*df_rec_lo[6]+df_rec_lo[2]*Dsurf_CC_vx[6]+Dsurf_CC_vx[1]*df_rec_lo[5]+df_rec_lo[1]*Dsurf_CC_vx[5]+Dsurf_CC_vx[0]*df_rec_lo[3]+df_rec_lo[0]*Dsurf_CC_vx[3]); 
  surft1_lo[4] = 0.3535533905932737*(Dsurf_CC_vx[3]*df_rec_lo[7]+df_rec_lo[3]*Dsurf_CC_vx[7]+Dsurf_CC_vx[5]*df_rec_lo[6]+df_rec_lo[5]*Dsurf_CC_vx[6]+Dsurf_CC_vx[0]*df_rec_lo[4]+df_rec_lo[0]*Dsurf_CC_vx[4]+Dsurf_CC_vx[1]*df_rec_lo[2]+df_rec_lo[1]*Dsurf_CC_vx[2]); 
  surft1_lo[5] = 0.3535533905932737*(Dsurf_CC_vx[2]*df_rec_lo[7]+df_rec_lo[2]*Dsurf_CC_vx[7]+Dsurf_CC_vx[4]*df_rec_lo[6]+df_rec_lo[4]*Dsurf_CC_vx[6]+Dsurf_CC_vx[0]*df_rec_lo[5]+df_rec_lo[0]*Dsurf_CC_vx[5]+Dsurf_CC_vx[1]*df_rec_lo[3]+df_rec_lo[1]*Dsurf_CC_vx[3]); 
  surft1_lo[6] = 0.3535533905932737*(Dsurf_CC_vx[1]*df_rec_lo[7]+df_rec_lo[1]*Dsurf_CC_vx[7]+Dsurf_CC_vx[0]*df_rec_lo[6]+df_rec_lo[0]*Dsurf_CC_vx[6]+Dsurf_CC_vx[4]*df_rec_lo[5]+df_rec_lo[4]*Dsurf_CC_vx[5]+Dsurf_CC_vx[2]*df_rec_lo[3]+df_rec_lo[2]*Dsurf_CC_vx[3]); 
  surft1_lo[7] = 0.3535533905932737*(Dsurf_CC_vx[0]*df_rec_lo[7]+df_rec_lo[0]*Dsurf_CC_vx[7]+Dsurf_CC_vx[1]*df_rec_lo[6]+df_rec_lo[1]*Dsurf_CC_vx[6]+Dsurf_CC_vx[2]*df_rec_lo[5]+df_rec_lo[2]*Dsurf_CC_vx[5]+Dsurf_CC_vx[3]*df_rec_lo[4]+df_rec_lo[3]*Dsurf_CC_vx[4]); 
  surft1_up[0] = 0.3535533905932737*(Dsurf_CR_vx[7]*df_rec_up[7]+Dsurf_CR_vx[6]*df_rec_up[6]+Dsurf_CR_vx[5]*df_rec_up[5]+Dsurf_CR_vx[4]*df_rec_up[4]+Dsurf_CR_vx[3]*df_rec_up[3]+Dsurf_CR_vx[2]*df_rec_up[2]+Dsurf_CR_vx[1]*df_rec_up[1]+Dsurf_CR_vx[0]*df_rec_up[0]); 
  surft1_up[1] = 0.3535533905932737*(Dsurf_CR_vx[6]*df_rec_up[7]+df_rec_up[6]*Dsurf_CR_vx[7]+Dsurf_CR_vx[3]*df_rec_up[5]+df_rec_up[3]*Dsurf_CR_vx[5]+Dsurf_CR_vx[2]*df_rec_up[4]+df_rec_up[2]*Dsurf_CR_vx[4]+Dsurf_CR_vx[0]*df_rec_up[1]+df_rec_up[0]*Dsurf_CR_vx[1]); 
  surft1_up[2] = 0.3535533905932737*(Dsurf_CR_vx[5]*df_rec_up[7]+df_rec_up[5]*Dsurf_CR_vx[7]+Dsurf_CR_vx[3]*df_rec_up[6]+df_rec_up[3]*Dsurf_CR_vx[6]+Dsurf_CR_vx[1]*df_rec_up[4]+df_rec_up[1]*Dsurf_CR_vx[4]+Dsurf_CR_vx[0]*df_rec_up[2]+df_rec_up[0]*Dsurf_CR_vx[2]); 
  surft1_up[3] = 0.3535533905932737*(Dsurf_CR_vx[4]*df_rec_up[7]+df_rec_up[4]*Dsurf_CR_vx[7]+Dsurf_CR_vx[2]*df_rec_up[6]+df_rec_up[2]*Dsurf_CR_vx[6]+Dsurf_CR_vx[1]*df_rec_up[5]+df_rec_up[1]*Dsurf_CR_vx[5]+Dsurf_CR_vx[0]*df_rec_up[3]+df_rec_up[0]*Dsurf_CR_vx[3]); 
  surft1_up[4] = 0.3535533905932737*(Dsurf_CR_vx[3]*df_rec_up[7]+df_rec_up[3]*Dsurf_CR_vx[7]+Dsurf_CR_vx[5]*df_rec_up[6]+df_rec_up[5]*Dsurf_CR_vx[6]+Dsurf_CR_vx[0]*df_rec_up[4]+df_rec_up[0]*Dsurf_CR_vx[4]+Dsurf_CR_vx[1]*df_rec_up[2]+df_rec_up[1]*Dsurf_CR_vx[2]); 
  surft1_up[5] = 0.3535533905932737*(Dsurf_CR_vx[2]*df_rec_up[7]+df_rec_up[2]*Dsurf_CR_vx[7]+Dsurf_CR_vx[4]*df_rec_up[6]+df_rec_up[4]*Dsurf_CR_vx[6]+Dsurf_CR_vx[0]*df_rec_up[5]+df_rec_up[0]*Dsurf_CR_vx[5]+Dsurf_CR_vx[1]*df_rec_up[3]+df_rec_up[1]*Dsurf_CR_vx[3]); 
  surft1_up[6] = 0.3535533905932737*(Dsurf_CR_vx[1]*df_rec_up[7]+df_rec_up[1]*Dsurf_CR_vx[7]+Dsurf_CR_vx[0]*df_rec_up[6]+df_rec_up[0]*Dsurf_CR_vx[6]+Dsurf_CR_vx[4]*df_rec_up[5]+df_rec_up[4]*Dsurf_CR_vx[5]+Dsurf_CR_vx[2]*df_rec_up[3]+df_rec_up[2]*Dsurf_CR_vx[3]); 
  surft1_up[7] = 0.3535533905932737*(Dsurf_CR_vx[0]*df_rec_up[7]+df_rec_up[0]*Dsurf_CR_vx[7]+Dsurf_CR_vx[1]*df_rec_up[6]+df_rec_up[1]*Dsurf_CR_vx[6]+Dsurf_CR_vx[2]*df_rec_up[5]+df_rec_up[2]*Dsurf_CR_vx[5]+Dsurf_CR_vx[3]*df_rec_up[4]+df_rec_up[3]*Dsurf_CR_vx[4]); 

  surft2_lo[0] = 0.3535533905932737*(Dsurf_CC_vz[7]*f_rec_lo[7]+Dsurf_CC_vz[6]*f_rec_lo[6]+Dsurf_CC_vz[5]*f_rec_lo[5]+Dsurf_CC_vz[4]*f_rec_lo[4]+Dsurf_CC_vz[3]*f_rec_lo[3]+Dsurf_CC_vz[2]*f_rec_lo[2]+Dsurf_CC_vz[1]*f_rec_lo[1]+Dsurf_CC_vz[0]*f_rec_lo[0]); 
  surft2_lo[1] = 0.3535533905932737*(Dsurf_CC_vz[6]*f_rec_lo[7]+f_rec_lo[6]*Dsurf_CC_vz[7]+Dsurf_CC_vz[3]*f_rec_lo[5]+f_rec_lo[3]*Dsurf_CC_vz[5]+Dsurf_CC_vz[2]*f_rec_lo[4]+f_rec_lo[2]*Dsurf_CC_vz[4]+Dsurf_CC_vz[0]*f_rec_lo[1]+f_rec_lo[0]*Dsurf_CC_vz[1]); 
  surft2_lo[2] = 0.3535533905932737*(Dsurf_CC_vz[5]*f_rec_lo[7]+f_rec_lo[5]*Dsurf_CC_vz[7]+Dsurf_CC_vz[3]*f_rec_lo[6]+f_rec_lo[3]*Dsurf_CC_vz[6]+Dsurf_CC_vz[1]*f_rec_lo[4]+f_rec_lo[1]*Dsurf_CC_vz[4]+Dsurf_CC_vz[0]*f_rec_lo[2]+f_rec_lo[0]*Dsurf_CC_vz[2]); 
  surft2_lo[3] = 0.3535533905932737*(Dsurf_CC_vz[4]*f_rec_lo[7]+f_rec_lo[4]*Dsurf_CC_vz[7]+Dsurf_CC_vz[2]*f_rec_lo[6]+f_rec_lo[2]*Dsurf_CC_vz[6]+Dsurf_CC_vz[1]*f_rec_lo[5]+f_rec_lo[1]*Dsurf_CC_vz[5]+Dsurf_CC_vz[0]*f_rec_lo[3]+f_rec_lo[0]*Dsurf_CC_vz[3]); 
  surft2_lo[4] = 0.3535533905932737*(Dsurf_CC_vz[3]*f_rec_lo[7]+f_rec_lo[3]*Dsurf_CC_vz[7]+Dsurf_CC_vz[5]*f_rec_lo[6]+f_rec_lo[5]*Dsurf_CC_vz[6]+Dsurf_CC_vz[0]*f_rec_lo[4]+f_rec_lo[0]*Dsurf_CC_vz[4]+Dsurf_CC_vz[1]*f_rec_lo[2]+f_rec_lo[1]*Dsurf_CC_vz[2]); 
  surft2_lo[5] = 0.3535533905932737*(Dsurf_CC_vz[2]*f_rec_lo[7]+f_rec_lo[2]*Dsurf_CC_vz[7]+Dsurf_CC_vz[4]*f_rec_lo[6]+f_rec_lo[4]*Dsurf_CC_vz[6]+Dsurf_CC_vz[0]*f_rec_lo[5]+f_rec_lo[0]*Dsurf_CC_vz[5]+Dsurf_CC_vz[1]*f_rec_lo[3]+f_rec_lo[1]*Dsurf_CC_vz[3]); 
  surft2_lo[6] = 0.3535533905932737*(Dsurf_CC_vz[1]*f_rec_lo[7]+f_rec_lo[1]*Dsurf_CC_vz[7]+Dsurf_CC_vz[0]*f_rec_lo[6]+f_rec_lo[0]*Dsurf_CC_vz[6]+Dsurf_CC_vz[4]*f_rec_lo[5]+f_rec_lo[4]*Dsurf_CC_vz[5]+Dsurf_CC_vz[2]*f_rec_lo[3]+f_rec_lo[2]*Dsurf_CC_vz[3]); 
  surft2_lo[7] = 0.3535533905932737*(Dsurf_CC_vz[0]*f_rec_lo[7]+f_rec_lo[0]*Dsurf_CC_vz[7]+Dsurf_CC_vz[1]*f_rec_lo[6]+f_rec_lo[1]*Dsurf_CC_vz[6]+Dsurf_CC_vz[2]*f_rec_lo[5]+f_rec_lo[2]*Dsurf_CC_vz[5]+Dsurf_CC_vz[3]*f_rec_lo[4]+f_rec_lo[3]*Dsurf_CC_vz[4]); 
  surft2_up[0] = 0.3535533905932737*(Dsurf_TC_vz[7]*f_rec_up[7]+Dsurf_TC_vz[6]*f_rec_up[6]+Dsurf_TC_vz[5]*f_rec_up[5]+Dsurf_TC_vz[4]*f_rec_up[4]+Dsurf_TC_vz[3]*f_rec_up[3]+Dsurf_TC_vz[2]*f_rec_up[2]+Dsurf_TC_vz[1]*f_rec_up[1]+Dsurf_TC_vz[0]*f_rec_up[0]); 
  surft2_up[1] = 0.3535533905932737*(Dsurf_TC_vz[6]*f_rec_up[7]+f_rec_up[6]*Dsurf_TC_vz[7]+Dsurf_TC_vz[3]*f_rec_up[5]+f_rec_up[3]*Dsurf_TC_vz[5]+Dsurf_TC_vz[2]*f_rec_up[4]+f_rec_up[2]*Dsurf_TC_vz[4]+Dsurf_TC_vz[0]*f_rec_up[1]+f_rec_up[0]*Dsurf_TC_vz[1]); 
  surft2_up[2] = 0.3535533905932737*(Dsurf_TC_vz[5]*f_rec_up[7]+f_rec_up[5]*Dsurf_TC_vz[7]+Dsurf_TC_vz[3]*f_rec_up[6]+f_rec_up[3]*Dsurf_TC_vz[6]+Dsurf_TC_vz[1]*f_rec_up[4]+f_rec_up[1]*Dsurf_TC_vz[4]+Dsurf_TC_vz[0]*f_rec_up[2]+f_rec_up[0]*Dsurf_TC_vz[2]); 
  surft2_up[3] = 0.3535533905932737*(Dsurf_TC_vz[4]*f_rec_up[7]+f_rec_up[4]*Dsurf_TC_vz[7]+Dsurf_TC_vz[2]*f_rec_up[6]+f_rec_up[2]*Dsurf_TC_vz[6]+Dsurf_TC_vz[1]*f_rec_up[5]+f_rec_up[1]*Dsurf_TC_vz[5]+Dsurf_TC_vz[0]*f_rec_up[3]+f_rec_up[0]*Dsurf_TC_vz[3]); 
  surft2_up[4] = 0.3535533905932737*(Dsurf_TC_vz[3]*f_rec_up[7]+f_rec_up[3]*Dsurf_TC_vz[7]+Dsurf_TC_vz[5]*f_rec_up[6]+f_rec_up[5]*Dsurf_TC_vz[6]+Dsurf_TC_vz[0]*f_rec_up[4]+f_rec_up[0]*Dsurf_TC_vz[4]+Dsurf_TC_vz[1]*f_rec_up[2]+f_rec_up[1]*Dsurf_TC_vz[2]); 
  surft2_up[5] = 0.3535533905932737*(Dsurf_TC_vz[2]*f_rec_up[7]+f_rec_up[2]*Dsurf_TC_vz[7]+Dsurf_TC_vz[4]*f_rec_up[6]+f_rec_up[4]*Dsurf_TC_vz[6]+Dsurf_TC_vz[0]*f_rec_up[5]+f_rec_up[0]*Dsurf_TC_vz[5]+Dsurf_TC_vz[1]*f_rec_up[3]+f_rec_up[1]*Dsurf_TC_vz[3]); 
  surft2_up[6] = 0.3535533905932737*(Dsurf_TC_vz[1]*f_rec_up[7]+f_rec_up[1]*Dsurf_TC_vz[7]+Dsurf_TC_vz[0]*f_rec_up[6]+f_rec_up[0]*Dsurf_TC_vz[6]+Dsurf_TC_vz[4]*f_rec_up[5]+f_rec_up[4]*Dsurf_TC_vz[5]+Dsurf_TC_vz[2]*f_rec_up[3]+f_rec_up[2]*Dsurf_TC_vz[3]); 
  surft2_up[7] = 0.3535533905932737*(Dsurf_TC_vz[0]*f_rec_up[7]+f_rec_up[0]*Dsurf_TC_vz[7]+Dsurf_TC_vz[1]*f_rec_up[6]+f_rec_up[1]*Dsurf_TC_vz[6]+Dsurf_TC_vz[2]*f_rec_up[5]+f_rec_up[2]*Dsurf_TC_vz[5]+Dsurf_TC_vz[3]*f_rec_up[4]+f_rec_up[3]*Dsurf_TC_vz[4]); 

  vol[4] = 0.75*(fCC[13]*DCC[15]+fCC[10]*DCC[14]+fCC[8]*DCC[12]+fCC[6]*DCC[11]+fCC[4]*DCC[9]+fCC[3]*DCC[7]+fCC[1]*DCC[5]+fCC[0]*DCC[2]); 
  vol[8] = 0.75*(fCC[10]*DCC[15]+fCC[13]*DCC[14]+fCC[4]*DCC[12]+fCC[3]*DCC[11]+fCC[8]*DCC[9]+fCC[6]*DCC[7]+fCC[0]*DCC[5]+fCC[1]*DCC[2]); 
  vol[9] = 1.5*(DCC[15]*fCC[15]+DCC[14]*fCC[14]+DCC[12]*fCC[12]+DCC[11]*fCC[11]+DCC[9]*fCC[9]+DCC[7]*fCC[7]+DCC[5]*fCC[5]+DCC[2]*fCC[2])+0.75*(DCC[13]*fCC[13]+DCC[10]*fCC[10]+DCC[8]*fCC[8]+DCC[6]*fCC[6]+DCC[4]*fCC[4]+DCC[3]*fCC[3]+DCC[1]*fCC[1]+DCC[0]*fCC[0]); 
  vol[10] = 0.75*(fCC[8]*DCC[15]+fCC[4]*DCC[14]+DCC[12]*fCC[13]+fCC[1]*DCC[11]+DCC[9]*fCC[10]+fCC[0]*DCC[7]+DCC[5]*fCC[6]+DCC[2]*fCC[3]); 
  vol[12] = 1.5*(DCC[14]*fCC[15]+fCC[14]*DCC[15]+DCC[9]*fCC[12]+fCC[9]*DCC[12]+DCC[7]*fCC[11]+fCC[7]*DCC[11]+DCC[2]*fCC[5]+fCC[2]*DCC[5])+0.75*(DCC[10]*fCC[13]+fCC[10]*DCC[13]+DCC[4]*fCC[8]+fCC[4]*DCC[8]+DCC[3]*fCC[6]+fCC[3]*DCC[6]+DCC[0]*fCC[1]+fCC[0]*DCC[1]); 
  vol[13] = 0.75*(fCC[4]*DCC[15]+fCC[8]*DCC[14]+DCC[9]*fCC[13]+fCC[10]*DCC[12]+fCC[0]*DCC[11]+fCC[1]*DCC[7]+DCC[2]*fCC[6]+fCC[3]*DCC[5]); 
  vol[14] = 1.5*(DCC[12]*fCC[15]+fCC[12]*DCC[15]+DCC[9]*fCC[14]+fCC[9]*DCC[14]+DCC[5]*fCC[11]+fCC[5]*DCC[11]+DCC[2]*fCC[7]+fCC[2]*DCC[7])+0.75*(DCC[8]*fCC[13]+fCC[8]*DCC[13]+DCC[4]*fCC[10]+fCC[4]*DCC[10]+DCC[1]*fCC[6]+fCC[1]*DCC[6]+DCC[0]*fCC[3]+fCC[0]*DCC[3]); 
  vol[15] = 1.5*(DCC[9]*fCC[15]+fCC[9]*DCC[15]+DCC[12]*fCC[14]+fCC[12]*DCC[14]+DCC[2]*fCC[11]+fCC[2]*DCC[11]+DCC[5]*fCC[7]+fCC[5]*DCC[7])+0.75*(DCC[4]*fCC[13]+fCC[4]*DCC[13]+DCC[8]*fCC[10]+fCC[8]*DCC[10]+DCC[0]*fCC[6]+fCC[0]*DCC[6]+DCC[1]*fCC[3]+fCC[1]*DCC[3]); 

  out[0] += (0.5*vol[0]+0.35355339059327373*surft1_up[0]-0.35355339059327373*surft1_lo[0])*dv_inv_sq; 
  out[1] += (0.5*vol[1]+0.35355339059327373*surft1_up[1]-0.35355339059327373*surft1_lo[1])*dv_inv_sq; 
  out[2] += (0.5*vol[2]-0.6123724356957945*surft2_up[0]+0.6123724356957945*(surft2_lo[0]+surft1_up[0]+surft1_lo[0]))*dv_inv_sq; 
  out[3] += (0.5*vol[3]+0.35355339059327373*surft1_up[2]-0.35355339059327373*surft1_lo[2])*dv_inv_sq; 
  out[4] += (0.5*vol[4]+0.35355339059327373*surft1_up[3]-0.35355339059327373*surft1_lo[3])*dv_inv_sq; 
  out[5] += (0.5*vol[5]-0.6123724356957945*surft2_up[1]+0.6123724356957945*(surft2_lo[1]+surft1_up[1]+surft1_lo[1]))*dv_inv_sq; 
  out[6] += (0.5*vol[6]+0.35355339059327373*surft1_up[4]-0.35355339059327373*surft1_lo[4])*dv_inv_sq; 
  out[7] += (0.5*vol[7]-0.6123724356957945*surft2_up[3]+0.6123724356957945*(surft2_lo[3]+surft1_up[2]+surft1_lo[2]))*dv_inv_sq; 
  out[8] += (0.5*vol[8]+0.35355339059327373*surft1_up[5]-0.35355339059327373*surft1_lo[5])*dv_inv_sq; 
  out[9] += (0.5*vol[9]+0.6123724356957945*(surft1_up[3]+surft1_lo[3])-1.0606601717798212*(surft2_up[0]+surft2_lo[0]))*dv_inv_sq; 
  out[10] += (0.5*vol[10]+0.35355339059327373*surft1_up[6]-0.35355339059327373*surft1_lo[6])*dv_inv_sq; 
  out[11] += (0.5*vol[11]-0.6123724356957945*surft2_up[5]+0.6123724356957945*(surft2_lo[5]+surft1_up[4]+surft1_lo[4]))*dv_inv_sq; 
  out[12] += (0.5*vol[12]+0.6123724356957945*(surft1_up[5]+surft1_lo[5])-1.0606601717798212*(surft2_up[1]+surft2_lo[1]))*dv_inv_sq; 
  out[13] += (0.5*vol[13]+0.35355339059327373*surft1_up[7]-0.35355339059327373*surft1_lo[7])*dv_inv_sq; 
  out[14] += (0.5*vol[14]+0.6123724356957945*(surft1_up[6]+surft1_lo[6])-1.0606601717798212*(surft2_up[3]+surft2_lo[3]))*dv_inv_sq; 
  out[15] += (0.5*vol[15]+0.6123724356957945*(surft1_up[7]+surft1_lo[7])-1.0606601717798212*(surft2_up[5]+surft2_lo[5]))*dv_inv_sq; 
  double cflFreq = fabs(Dsurf_CR_vx[0]); 

  return 1.4142135623730947*dv_inv_sq*cflFreq; 
} 
