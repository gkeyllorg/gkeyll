#include <gkyl_fpo_vlasov_kernels.h> 

GKYL_CU_DH double fpo_vlasov_diff_surfvyvz_1x3v_ser_p1(const double* dxv, const double* diff_coeff_C, 
  const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* GKYL_RESTRICT out) { 
  // dxv[NDIM]: Cell spacing in each direction. 
  // diff_coeff_C: Diffusion tensor in center cell. 
  // diff_coeff_surf_stencil[9]: 9-cell stencil of surface expansion of recovered diffusion tensor. 
  // f_stencil[9]: 9-cell stencil of distribution function. 
  // out: Incremented output. 

  double dv_inv_sq = 4.0/dxv[2]/dxv[3]; 

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
  const double *fBL = f_stencil[0]; 
  const double *fCL = f_stencil[1]; 
  const double *fTL = f_stencil[2]; 
  const double *fBC = f_stencil[3]; 
  const double *fCC = f_stencil[4]; 
  const double *fTC = f_stencil[5]; 
  const double *fBR = f_stencil[6]; 
  const double *fCR = f_stencil[7]; 
  const double *fTR = f_stencil[8]; 

  const double *DCC = &diff_coeff_C[80]; 

  const double *Dsurf_CC_vy = &diff_coeff_surf_stencil[4][80]; 
  const double *Dsurf_CC_vz = &diff_coeff_surf_stencil[4][88]; 
  const double* Dsurf_TC_vz = &diff_coeff_surf_stencil[5][88]; 
  const double* Dsurf_CR_vy = &diff_coeff_surf_stencil[7][80]; 
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

  df_rec_lo[0] = -(0.11785113019775789*(fTL[10]+fBL[10]))+0.11785113019775789*(fTC[10]+fBC[10])+0.2357022603955158*fCL[10]-0.2357022603955158*fCC[10]-0.10206207261596573*(fTL[4]+fTC[4]+fBL[4]+fBC[4]+fTC[3]+fBL[3])+0.20412414523193148*(fCL[4]+fCC[4])+0.10206207261596573*(fTL[3]+fBC[3])+0.0883883476483184*(fTL[0]+fTC[0])-0.0883883476483184*(fBL[0]+fBC[0]); 
  df_rec_lo[1] = -(0.11785113019775789*(fTL[13]+fBL[13]))+0.11785113019775789*(fTC[13]+fBC[13])+0.2357022603955158*fCL[13]-0.2357022603955158*fCC[13]-0.10206207261596573*(fTL[8]+fTC[8]+fBL[8]+fBC[8]+fTC[6]+fBL[6])+0.20412414523193148*(fCL[8]+fCC[8])+0.10206207261596573*(fTL[6]+fBC[6])+0.0883883476483184*(fTL[1]+fTC[1])-0.0883883476483184*(fBL[1]+fBC[1]); 
  df_rec_lo[2] = -(0.11785113019775789*(fTL[14]+fBL[14]))+0.11785113019775789*(fTC[14]+fBC[14])+0.2357022603955158*fCL[14]-0.2357022603955158*fCC[14]-0.10206207261596573*(fTL[9]+fTC[9]+fBL[9]+fBC[9]+fTC[7]+fBL[7])+0.20412414523193148*(fCL[9]+fCC[9])+0.10206207261596573*(fTL[7]+fBC[7])+0.0883883476483184*(fTL[2]+fTC[2])-0.0883883476483184*(fBL[2]+fBC[2]); 
  df_rec_lo[3] = -(0.20412414523193148*(fTL[10]+fBC[10]))+0.20412414523193148*(fTC[10]+fBL[10])-0.1767766952966368*(fTL[4]+fTC[4]+fTC[3]+fBC[3])+0.1767766952966368*(fBL[4]+fBC[4]+fTL[3]+fBL[3])-0.3535533905932737*fCL[3]+0.3535533905932737*fCC[3]+0.15309310892394856*(fTL[0]+fTC[0]+fBL[0]+fBC[0])-0.3061862178478971*(fCL[0]+fCC[0]); 
  df_rec_lo[4] = -(0.11785113019775789*(fTL[15]+fBL[15]))+0.11785113019775789*(fTC[15]+fBC[15])+0.2357022603955158*fCL[15]-0.2357022603955158*fCC[15]-0.10206207261596573*(fTL[12]+fTC[12]+fBL[12]+fBC[12]+fTC[11]+fBL[11])+0.20412414523193148*(fCL[12]+fCC[12])+0.10206207261596573*(fTL[11]+fBC[11])+0.0883883476483184*(fTL[5]+fTC[5])-0.0883883476483184*(fBL[5]+fBC[5]); 
  df_rec_lo[5] = -(0.20412414523193148*(fTL[13]+fBC[13]))+0.20412414523193148*(fTC[13]+fBL[13])-0.1767766952966368*(fTL[8]+fTC[8]+fTC[6]+fBC[6])+0.1767766952966368*(fBL[8]+fBC[8]+fTL[6]+fBL[6])-0.3535533905932737*fCL[6]+0.3535533905932737*fCC[6]+0.15309310892394856*(fTL[1]+fTC[1]+fBL[1]+fBC[1])-0.3061862178478971*(fCL[1]+fCC[1]); 
  df_rec_lo[6] = -(0.20412414523193148*(fTL[14]+fBC[14]))+0.20412414523193148*(fTC[14]+fBL[14])-0.1767766952966368*(fTL[9]+fTC[9]+fTC[7]+fBC[7])+0.1767766952966368*(fBL[9]+fBC[9]+fTL[7]+fBL[7])-0.3535533905932737*fCL[7]+0.3535533905932737*fCC[7]+0.15309310892394856*(fTL[2]+fTC[2]+fBL[2]+fBC[2])-0.3061862178478971*(fCL[2]+fCC[2]); 
  df_rec_lo[7] = -(0.20412414523193148*(fTL[15]+fBC[15]))+0.20412414523193148*(fTC[15]+fBL[15])-0.1767766952966368*(fTL[12]+fTC[12]+fTC[11]+fBC[11])+0.1767766952966368*(fBL[12]+fBC[12]+fTL[11]+fBL[11])-0.3535533905932737*fCL[11]+0.3535533905932737*fCC[11]+0.15309310892394856*(fTL[5]+fTC[5]+fBL[5]+fBC[5])-0.3061862178478971*(fCL[5]+fCC[5]); 
  df_rec_up[0] = 0.11785113019775789*(fTR[10]+fBR[10])-0.11785113019775789*(fTC[10]+fBC[10])-0.2357022603955158*fCR[10]+0.2357022603955158*fCC[10]-0.10206207261596573*(fTR[4]+fTC[4]+fBR[4]+fBC[4]+fTR[3]+fBC[3])+0.20412414523193148*(fCR[4]+fCC[4])+0.10206207261596573*(fTC[3]+fBR[3])+0.0883883476483184*(fTR[0]+fTC[0])-0.0883883476483184*(fBR[0]+fBC[0]); 
  df_rec_up[1] = 0.11785113019775789*(fTR[13]+fBR[13])-0.11785113019775789*(fTC[13]+fBC[13])-0.2357022603955158*fCR[13]+0.2357022603955158*fCC[13]-0.10206207261596573*(fTR[8]+fTC[8]+fBR[8]+fBC[8]+fTR[6]+fBC[6])+0.20412414523193148*(fCR[8]+fCC[8])+0.10206207261596573*(fTC[6]+fBR[6])+0.0883883476483184*(fTR[1]+fTC[1])-0.0883883476483184*(fBR[1]+fBC[1]); 
  df_rec_up[2] = 0.11785113019775789*(fTR[14]+fBR[14])-0.11785113019775789*(fTC[14]+fBC[14])-0.2357022603955158*fCR[14]+0.2357022603955158*fCC[14]-0.10206207261596573*(fTR[9]+fTC[9]+fBR[9]+fBC[9]+fTR[7]+fBC[7])+0.20412414523193148*(fCR[9]+fCC[9])+0.10206207261596573*(fTC[7]+fBR[7])+0.0883883476483184*(fTR[2]+fTC[2])-0.0883883476483184*(fBR[2]+fBC[2]); 
  df_rec_up[3] = 0.20412414523193148*(fTR[10]+fBC[10])-0.20412414523193148*(fTC[10]+fBR[10])-0.1767766952966368*(fTR[4]+fTC[4]+fTR[3]+fBR[3])+0.1767766952966368*(fBR[4]+fBC[4]+fTC[3]+fBC[3])+0.3535533905932737*fCR[3]-0.3535533905932737*fCC[3]+0.15309310892394856*(fTR[0]+fTC[0]+fBR[0]+fBC[0])-0.3061862178478971*(fCR[0]+fCC[0]); 
  df_rec_up[4] = 0.11785113019775789*(fTR[15]+fBR[15])-0.11785113019775789*(fTC[15]+fBC[15])-0.2357022603955158*fCR[15]+0.2357022603955158*fCC[15]-0.10206207261596573*(fTR[12]+fTC[12]+fBR[12]+fBC[12]+fTR[11]+fBC[11])+0.20412414523193148*(fCR[12]+fCC[12])+0.10206207261596573*(fTC[11]+fBR[11])+0.0883883476483184*(fTR[5]+fTC[5])-0.0883883476483184*(fBR[5]+fBC[5]); 
  df_rec_up[5] = 0.20412414523193148*(fTR[13]+fBC[13])-0.20412414523193148*(fTC[13]+fBR[13])-0.1767766952966368*(fTR[8]+fTC[8]+fTR[6]+fBR[6])+0.1767766952966368*(fBR[8]+fBC[8]+fTC[6]+fBC[6])+0.3535533905932737*fCR[6]-0.3535533905932737*fCC[6]+0.15309310892394856*(fTR[1]+fTC[1]+fBR[1]+fBC[1])-0.3061862178478971*(fCR[1]+fCC[1]); 
  df_rec_up[6] = 0.20412414523193148*(fTR[14]+fBC[14])-0.20412414523193148*(fTC[14]+fBR[14])-0.1767766952966368*(fTR[9]+fTC[9]+fTR[7]+fBR[7])+0.1767766952966368*(fBR[9]+fBC[9]+fTC[7]+fBC[7])+0.3535533905932737*fCR[7]-0.3535533905932737*fCC[7]+0.15309310892394856*(fTR[2]+fTC[2]+fBR[2]+fBC[2])-0.3061862178478971*(fCR[2]+fCC[2]); 
  df_rec_up[7] = 0.20412414523193148*(fTR[15]+fBC[15])-0.20412414523193148*(fTC[15]+fBR[15])-0.1767766952966368*(fTR[12]+fTC[12]+fTR[11]+fBR[11])+0.1767766952966368*(fBR[12]+fBC[12]+fTC[11]+fBC[11])+0.3535533905932737*fCR[11]-0.3535533905932737*fCC[11]+0.15309310892394856*(fTR[5]+fTC[5]+fBR[5]+fBC[5])-0.3061862178478971*(fCR[5]+fCC[5]); 

  surft1_lo[0] = 0.3535533905932737*(Dsurf_CC_vy[7]*df_rec_lo[7]+Dsurf_CC_vy[6]*df_rec_lo[6]+Dsurf_CC_vy[5]*df_rec_lo[5]+Dsurf_CC_vy[4]*df_rec_lo[4]+Dsurf_CC_vy[3]*df_rec_lo[3]+Dsurf_CC_vy[2]*df_rec_lo[2]+Dsurf_CC_vy[1]*df_rec_lo[1]+Dsurf_CC_vy[0]*df_rec_lo[0]); 
  surft1_lo[1] = 0.3535533905932737*(Dsurf_CC_vy[6]*df_rec_lo[7]+df_rec_lo[6]*Dsurf_CC_vy[7]+Dsurf_CC_vy[3]*df_rec_lo[5]+df_rec_lo[3]*Dsurf_CC_vy[5]+Dsurf_CC_vy[2]*df_rec_lo[4]+df_rec_lo[2]*Dsurf_CC_vy[4]+Dsurf_CC_vy[0]*df_rec_lo[1]+df_rec_lo[0]*Dsurf_CC_vy[1]); 
  surft1_lo[2] = 0.3535533905932737*(Dsurf_CC_vy[5]*df_rec_lo[7]+df_rec_lo[5]*Dsurf_CC_vy[7]+Dsurf_CC_vy[3]*df_rec_lo[6]+df_rec_lo[3]*Dsurf_CC_vy[6]+Dsurf_CC_vy[1]*df_rec_lo[4]+df_rec_lo[1]*Dsurf_CC_vy[4]+Dsurf_CC_vy[0]*df_rec_lo[2]+df_rec_lo[0]*Dsurf_CC_vy[2]); 
  surft1_lo[3] = 0.3535533905932737*(Dsurf_CC_vy[4]*df_rec_lo[7]+df_rec_lo[4]*Dsurf_CC_vy[7]+Dsurf_CC_vy[2]*df_rec_lo[6]+df_rec_lo[2]*Dsurf_CC_vy[6]+Dsurf_CC_vy[1]*df_rec_lo[5]+df_rec_lo[1]*Dsurf_CC_vy[5]+Dsurf_CC_vy[0]*df_rec_lo[3]+df_rec_lo[0]*Dsurf_CC_vy[3]); 
  surft1_lo[4] = 0.3535533905932737*(Dsurf_CC_vy[3]*df_rec_lo[7]+df_rec_lo[3]*Dsurf_CC_vy[7]+Dsurf_CC_vy[5]*df_rec_lo[6]+df_rec_lo[5]*Dsurf_CC_vy[6]+Dsurf_CC_vy[0]*df_rec_lo[4]+df_rec_lo[0]*Dsurf_CC_vy[4]+Dsurf_CC_vy[1]*df_rec_lo[2]+df_rec_lo[1]*Dsurf_CC_vy[2]); 
  surft1_lo[5] = 0.3535533905932737*(Dsurf_CC_vy[2]*df_rec_lo[7]+df_rec_lo[2]*Dsurf_CC_vy[7]+Dsurf_CC_vy[4]*df_rec_lo[6]+df_rec_lo[4]*Dsurf_CC_vy[6]+Dsurf_CC_vy[0]*df_rec_lo[5]+df_rec_lo[0]*Dsurf_CC_vy[5]+Dsurf_CC_vy[1]*df_rec_lo[3]+df_rec_lo[1]*Dsurf_CC_vy[3]); 
  surft1_lo[6] = 0.3535533905932737*(Dsurf_CC_vy[1]*df_rec_lo[7]+df_rec_lo[1]*Dsurf_CC_vy[7]+Dsurf_CC_vy[0]*df_rec_lo[6]+df_rec_lo[0]*Dsurf_CC_vy[6]+Dsurf_CC_vy[4]*df_rec_lo[5]+df_rec_lo[4]*Dsurf_CC_vy[5]+Dsurf_CC_vy[2]*df_rec_lo[3]+df_rec_lo[2]*Dsurf_CC_vy[3]); 
  surft1_lo[7] = 0.3535533905932737*(Dsurf_CC_vy[0]*df_rec_lo[7]+df_rec_lo[0]*Dsurf_CC_vy[7]+Dsurf_CC_vy[1]*df_rec_lo[6]+df_rec_lo[1]*Dsurf_CC_vy[6]+Dsurf_CC_vy[2]*df_rec_lo[5]+df_rec_lo[2]*Dsurf_CC_vy[5]+Dsurf_CC_vy[3]*df_rec_lo[4]+df_rec_lo[3]*Dsurf_CC_vy[4]); 
  surft1_up[0] = 0.3535533905932737*(Dsurf_CR_vy[7]*df_rec_up[7]+Dsurf_CR_vy[6]*df_rec_up[6]+Dsurf_CR_vy[5]*df_rec_up[5]+Dsurf_CR_vy[4]*df_rec_up[4]+Dsurf_CR_vy[3]*df_rec_up[3]+Dsurf_CR_vy[2]*df_rec_up[2]+Dsurf_CR_vy[1]*df_rec_up[1]+Dsurf_CR_vy[0]*df_rec_up[0]); 
  surft1_up[1] = 0.3535533905932737*(Dsurf_CR_vy[6]*df_rec_up[7]+df_rec_up[6]*Dsurf_CR_vy[7]+Dsurf_CR_vy[3]*df_rec_up[5]+df_rec_up[3]*Dsurf_CR_vy[5]+Dsurf_CR_vy[2]*df_rec_up[4]+df_rec_up[2]*Dsurf_CR_vy[4]+Dsurf_CR_vy[0]*df_rec_up[1]+df_rec_up[0]*Dsurf_CR_vy[1]); 
  surft1_up[2] = 0.3535533905932737*(Dsurf_CR_vy[5]*df_rec_up[7]+df_rec_up[5]*Dsurf_CR_vy[7]+Dsurf_CR_vy[3]*df_rec_up[6]+df_rec_up[3]*Dsurf_CR_vy[6]+Dsurf_CR_vy[1]*df_rec_up[4]+df_rec_up[1]*Dsurf_CR_vy[4]+Dsurf_CR_vy[0]*df_rec_up[2]+df_rec_up[0]*Dsurf_CR_vy[2]); 
  surft1_up[3] = 0.3535533905932737*(Dsurf_CR_vy[4]*df_rec_up[7]+df_rec_up[4]*Dsurf_CR_vy[7]+Dsurf_CR_vy[2]*df_rec_up[6]+df_rec_up[2]*Dsurf_CR_vy[6]+Dsurf_CR_vy[1]*df_rec_up[5]+df_rec_up[1]*Dsurf_CR_vy[5]+Dsurf_CR_vy[0]*df_rec_up[3]+df_rec_up[0]*Dsurf_CR_vy[3]); 
  surft1_up[4] = 0.3535533905932737*(Dsurf_CR_vy[3]*df_rec_up[7]+df_rec_up[3]*Dsurf_CR_vy[7]+Dsurf_CR_vy[5]*df_rec_up[6]+df_rec_up[5]*Dsurf_CR_vy[6]+Dsurf_CR_vy[0]*df_rec_up[4]+df_rec_up[0]*Dsurf_CR_vy[4]+Dsurf_CR_vy[1]*df_rec_up[2]+df_rec_up[1]*Dsurf_CR_vy[2]); 
  surft1_up[5] = 0.3535533905932737*(Dsurf_CR_vy[2]*df_rec_up[7]+df_rec_up[2]*Dsurf_CR_vy[7]+Dsurf_CR_vy[4]*df_rec_up[6]+df_rec_up[4]*Dsurf_CR_vy[6]+Dsurf_CR_vy[0]*df_rec_up[5]+df_rec_up[0]*Dsurf_CR_vy[5]+Dsurf_CR_vy[1]*df_rec_up[3]+df_rec_up[1]*Dsurf_CR_vy[3]); 
  surft1_up[6] = 0.3535533905932737*(Dsurf_CR_vy[1]*df_rec_up[7]+df_rec_up[1]*Dsurf_CR_vy[7]+Dsurf_CR_vy[0]*df_rec_up[6]+df_rec_up[0]*Dsurf_CR_vy[6]+Dsurf_CR_vy[4]*df_rec_up[5]+df_rec_up[4]*Dsurf_CR_vy[5]+Dsurf_CR_vy[2]*df_rec_up[3]+df_rec_up[2]*Dsurf_CR_vy[3]); 
  surft1_up[7] = 0.3535533905932737*(Dsurf_CR_vy[0]*df_rec_up[7]+df_rec_up[0]*Dsurf_CR_vy[7]+Dsurf_CR_vy[1]*df_rec_up[6]+df_rec_up[1]*Dsurf_CR_vy[6]+Dsurf_CR_vy[2]*df_rec_up[5]+df_rec_up[2]*Dsurf_CR_vy[5]+Dsurf_CR_vy[3]*df_rec_up[4]+df_rec_up[3]*Dsurf_CR_vy[4]); 

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

  vol[4] = 0.75*(fCC[12]*DCC[15]+fCC[9]*DCC[14]+fCC[8]*DCC[13]+fCC[5]*DCC[11]+fCC[4]*DCC[10]+fCC[2]*DCC[7]+fCC[1]*DCC[6]+fCC[0]*DCC[3]); 
  vol[8] = 0.75*(fCC[9]*DCC[15]+fCC[12]*DCC[14]+fCC[4]*DCC[13]+fCC[2]*DCC[11]+fCC[8]*DCC[10]+fCC[5]*DCC[7]+fCC[0]*DCC[6]+fCC[1]*DCC[3]); 
  vol[9] = 0.75*(fCC[8]*DCC[15]+fCC[4]*DCC[14]+fCC[12]*DCC[13]+fCC[1]*DCC[11]+fCC[9]*DCC[10]+fCC[0]*DCC[7]+fCC[5]*DCC[6]+fCC[2]*DCC[3]); 
  vol[10] = 1.5*(DCC[15]*fCC[15]+DCC[14]*fCC[14]+DCC[13]*fCC[13]+DCC[11]*fCC[11]+DCC[10]*fCC[10]+DCC[7]*fCC[7]+DCC[6]*fCC[6]+DCC[3]*fCC[3])+0.75*(DCC[12]*fCC[12]+DCC[9]*fCC[9]+DCC[8]*fCC[8]+DCC[5]*fCC[5]+DCC[4]*fCC[4]+DCC[2]*fCC[2]+DCC[1]*fCC[1]+DCC[0]*fCC[0]); 
  vol[12] = 0.75*(fCC[4]*DCC[15]+fCC[8]*DCC[14]+fCC[9]*DCC[13]+DCC[10]*fCC[12]+fCC[0]*DCC[11]+fCC[1]*DCC[7]+fCC[2]*DCC[6]+DCC[3]*fCC[5]); 
  vol[13] = 1.5*(DCC[14]*fCC[15]+fCC[14]*DCC[15]+DCC[10]*fCC[13]+fCC[10]*DCC[13]+DCC[7]*fCC[11]+fCC[7]*DCC[11]+DCC[3]*fCC[6]+fCC[3]*DCC[6])+0.75*(DCC[9]*fCC[12]+fCC[9]*DCC[12]+DCC[4]*fCC[8]+fCC[4]*DCC[8]+DCC[2]*fCC[5]+fCC[2]*DCC[5]+DCC[0]*fCC[1]+fCC[0]*DCC[1]); 
  vol[14] = 1.5*(DCC[13]*fCC[15]+fCC[13]*DCC[15]+DCC[10]*fCC[14]+fCC[10]*DCC[14]+DCC[6]*fCC[11]+fCC[6]*DCC[11]+DCC[3]*fCC[7]+fCC[3]*DCC[7])+0.75*(DCC[8]*fCC[12]+fCC[8]*DCC[12]+DCC[4]*fCC[9]+fCC[4]*DCC[9]+DCC[1]*fCC[5]+fCC[1]*DCC[5]+DCC[0]*fCC[2]+fCC[0]*DCC[2]); 
  vol[15] = 1.5*(DCC[10]*fCC[15]+fCC[10]*DCC[15]+DCC[13]*fCC[14]+fCC[13]*DCC[14]+DCC[3]*fCC[11]+fCC[3]*DCC[11]+DCC[6]*fCC[7]+fCC[6]*DCC[7])+0.75*(DCC[4]*fCC[12]+fCC[4]*DCC[12]+DCC[8]*fCC[9]+fCC[8]*DCC[9]+DCC[0]*fCC[5]+fCC[0]*DCC[5]+DCC[1]*fCC[2]+fCC[1]*DCC[2]); 

  out[0] += (0.5*vol[0]+0.35355339059327373*surft1_up[0]-0.35355339059327373*surft1_lo[0])*dv_inv_sq; 
  out[1] += (0.5*vol[1]+0.35355339059327373*surft1_up[1]-0.35355339059327373*surft1_lo[1])*dv_inv_sq; 
  out[2] += (0.5*vol[2]+0.35355339059327373*surft1_up[2]-0.35355339059327373*surft1_lo[2])*dv_inv_sq; 
  out[3] += (0.5*vol[3]-0.6123724356957945*surft2_up[0]+0.6123724356957945*(surft2_lo[0]+surft1_up[0]+surft1_lo[0]))*dv_inv_sq; 
  out[4] += (0.5*vol[4]+0.35355339059327373*surft1_up[3]-0.35355339059327373*surft1_lo[3])*dv_inv_sq; 
  out[5] += (0.5*vol[5]+0.35355339059327373*surft1_up[4]-0.35355339059327373*surft1_lo[4])*dv_inv_sq; 
  out[6] += (0.5*vol[6]-0.6123724356957945*surft2_up[1]+0.6123724356957945*(surft2_lo[1]+surft1_up[1]+surft1_lo[1]))*dv_inv_sq; 
  out[7] += (0.5*vol[7]-0.6123724356957945*surft2_up[2]+0.6123724356957945*(surft2_lo[2]+surft1_up[2]+surft1_lo[2]))*dv_inv_sq; 
  out[8] += (0.5*vol[8]+0.35355339059327373*surft1_up[5]-0.35355339059327373*surft1_lo[5])*dv_inv_sq; 
  out[9] += (0.5*vol[9]+0.35355339059327373*surft1_up[6]-0.35355339059327373*surft1_lo[6])*dv_inv_sq; 
  out[10] += (0.5*vol[10]+0.6123724356957945*(surft1_up[3]+surft1_lo[3])-1.0606601717798212*(surft2_up[0]+surft2_lo[0]))*dv_inv_sq; 
  out[11] += (0.5*vol[11]-0.6123724356957945*surft2_up[4]+0.6123724356957945*(surft2_lo[4]+surft1_up[4]+surft1_lo[4]))*dv_inv_sq; 
  out[12] += (0.5*vol[12]+0.35355339059327373*surft1_up[7]-0.35355339059327373*surft1_lo[7])*dv_inv_sq; 
  out[13] += (0.5*vol[13]+0.6123724356957945*(surft1_up[5]+surft1_lo[5])-1.0606601717798212*(surft2_up[1]+surft2_lo[1]))*dv_inv_sq; 
  out[14] += (0.5*vol[14]+0.6123724356957945*(surft1_up[6]+surft1_lo[6])-1.0606601717798212*(surft2_up[2]+surft2_lo[2]))*dv_inv_sq; 
  out[15] += (0.5*vol[15]+0.6123724356957945*(surft1_up[7]+surft1_lo[7])-1.0606601717798212*(surft2_up[4]+surft2_lo[4]))*dv_inv_sq; 

  double cflFreq = fmax(fabs(Dsurf_CC_vy[0]), fabs(Dsurf_CR_vy[0])); 

  return 1.4142135623730947*dv_inv_sq*cflFreq; 
} 
