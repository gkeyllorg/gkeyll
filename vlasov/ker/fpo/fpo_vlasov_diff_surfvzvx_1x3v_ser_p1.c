#include <gkyl_fpo_vlasov_kernels.h> 

GKYL_CU_DH double fpo_vlasov_diff_surfvzvx_1x3v_ser_p1(const double* dxv, const double* diff_coeff_C, 
  const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* GKYL_RESTRICT out) { 
  // dxv[NDIM]: Cell spacing in each direction. 
  // diff_coeff_C: Diffusion tensor in center cell. 
  // diff_coeff_surf_stencil[9]: 9-cell stencil of surface expansion of recovered diffusion tensor. 
  // f_stencil[9]: 9-cell stencil of distribution function. 
  // out: Incremented output. 

  double dv_inv_sq = 4.0/dxv[3]/dxv[1]; 

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

  const double *DCC = &diff_coeff_C[96]; 

  const double *Dsurf_CC_vz = &diff_coeff_surf_stencil[4][96]; 
  const double *Dsurf_CC_vx = &diff_coeff_surf_stencil[4][104]; 
  const double* Dsurf_TC_vz = &diff_coeff_surf_stencil[5][96]; 
  const double* Dsurf_CR_vx = &diff_coeff_surf_stencil[7][104]; 
  f_rec_lo[0] = 0.408248290463863*fCL[2]-0.408248290463863*fCC[2]+0.3535533905932737*(fCL[0]+fCC[0]); 
  f_rec_lo[1] = 0.408248290463863*fCL[5]-0.408248290463863*fCC[5]+0.3535533905932737*(fCL[1]+fCC[1]); 
  f_rec_lo[2] = 0.408248290463863*fCL[7]-0.408248290463863*fCC[7]+0.3535533905932737*(fCL[3]+fCC[3]); 
  f_rec_lo[3] = 0.408248290463863*fCL[9]-0.408248290463863*fCC[9]+0.3535533905932737*(fCL[4]+fCC[4]); 
  f_rec_lo[4] = 0.408248290463863*fCL[11]-0.408248290463863*fCC[11]+0.3535533905932737*(fCL[6]+fCC[6]); 
  f_rec_lo[5] = 0.408248290463863*fCL[12]-0.408248290463863*fCC[12]+0.3535533905932737*(fCL[8]+fCC[8]); 
  f_rec_lo[6] = 0.408248290463863*fCL[14]-0.408248290463863*fCC[14]+0.3535533905932737*(fCL[10]+fCC[10]); 
  f_rec_lo[7] = 0.408248290463863*fCL[15]-0.408248290463863*fCC[15]+0.3535533905932737*(fCL[13]+fCC[13]); 
  f_rec_up[0] = -(0.408248290463863*fCR[2])+0.408248290463863*fCC[2]+0.3535533905932737*(fCR[0]+fCC[0]); 
  f_rec_up[1] = -(0.408248290463863*fCR[5])+0.408248290463863*fCC[5]+0.3535533905932737*(fCR[1]+fCC[1]); 
  f_rec_up[2] = -(0.408248290463863*fCR[7])+0.408248290463863*fCC[7]+0.3535533905932737*(fCR[3]+fCC[3]); 
  f_rec_up[3] = -(0.408248290463863*fCR[9])+0.408248290463863*fCC[9]+0.3535533905932737*(fCR[4]+fCC[4]); 
  f_rec_up[4] = -(0.408248290463863*fCR[11])+0.408248290463863*fCC[11]+0.3535533905932737*(fCR[6]+fCC[6]); 
  f_rec_up[5] = -(0.408248290463863*fCR[12])+0.408248290463863*fCC[12]+0.3535533905932737*(fCR[8]+fCC[8]); 
  f_rec_up[6] = -(0.408248290463863*fCR[14])+0.408248290463863*fCC[14]+0.3535533905932737*(fCR[10]+fCC[10]); 
  f_rec_up[7] = -(0.408248290463863*fCR[15])+0.408248290463863*fCC[15]+0.3535533905932737*(fCR[13]+fCC[13]); 

  df_rec_lo[0] = 0.11785113019775789*(fCR[9]+fCL[9])-0.2357022603955158*fCC[9]-0.11785113019775789*(fBR[9]+fBL[9])+0.2357022603955158*fBC[9]-0.10206207261596573*(fCR[4]+fBL[4]+fCR[2]+fCL[2]+fBR[2]+fBL[2])+0.10206207261596573*(fCL[4]+fBR[4])+0.20412414523193148*(fCC[2]+fBC[2])+0.0883883476483184*(fCR[0]+fBR[0])-0.0883883476483184*(fCL[0]+fBL[0]); 
  df_rec_lo[1] = 0.11785113019775789*(fCR[12]+fCL[12])-0.2357022603955158*fCC[12]-0.11785113019775789*(fBR[12]+fBL[12])+0.2357022603955158*fBC[12]-0.10206207261596573*(fCR[8]+fBL[8]+fCR[5]+fCL[5]+fBR[5]+fBL[5])+0.10206207261596573*(fCL[8]+fBR[8])+0.20412414523193148*(fCC[5]+fBC[5])+0.0883883476483184*(fCR[1]+fBR[1])-0.0883883476483184*(fCL[1]+fBL[1]); 
  df_rec_lo[2] = 0.20412414523193148*(fCR[9]+fBL[9])-0.20412414523193148*(fCL[9]+fBR[9])-0.1767766952966368*(fCR[4]+fCL[4]+fCR[2]+fBR[2])+0.3535533905932737*fCC[4]+0.1767766952966368*(fBR[4]+fBL[4]+fCL[2]+fBL[2])-0.3535533905932737*fBC[4]+0.15309310892394856*(fCR[0]+fCL[0]+fBR[0]+fBL[0])-0.3061862178478971*(fCC[0]+fBC[0]); 
  df_rec_lo[3] = 0.11785113019775789*(fCR[14]+fCL[14])-0.2357022603955158*fCC[14]-0.11785113019775789*(fBR[14]+fBL[14])+0.2357022603955158*fBC[14]-0.10206207261596573*(fCR[10]+fBL[10]+fCR[7]+fCL[7]+fBR[7]+fBL[7])+0.10206207261596573*(fCL[10]+fBR[10])+0.20412414523193148*(fCC[7]+fBC[7])+0.0883883476483184*(fCR[3]+fBR[3])-0.0883883476483184*(fCL[3]+fBL[3]); 
  df_rec_lo[4] = 0.20412414523193148*(fCR[12]+fBL[12])-0.20412414523193148*(fCL[12]+fBR[12])-0.1767766952966368*(fCR[8]+fCL[8]+fCR[5]+fBR[5])+0.3535533905932737*fCC[8]+0.1767766952966368*(fBR[8]+fBL[8]+fCL[5]+fBL[5])-0.3535533905932737*fBC[8]+0.15309310892394856*(fCR[1]+fCL[1]+fBR[1]+fBL[1])-0.3061862178478971*(fCC[1]+fBC[1]); 
  df_rec_lo[5] = 0.11785113019775789*(fCR[15]+fCL[15])-0.2357022603955158*fCC[15]-0.11785113019775789*(fBR[15]+fBL[15])+0.2357022603955158*fBC[15]-0.10206207261596573*(fCR[13]+fBL[13]+fCR[11]+fCL[11]+fBR[11]+fBL[11])+0.10206207261596573*(fCL[13]+fBR[13])+0.20412414523193148*(fCC[11]+fBC[11])+0.0883883476483184*(fCR[6]+fBR[6])-0.0883883476483184*(fCL[6]+fBL[6]); 
  df_rec_lo[6] = 0.20412414523193148*(fCR[14]+fBL[14])-0.20412414523193148*(fCL[14]+fBR[14])-0.1767766952966368*(fCR[10]+fCL[10]+fCR[7]+fBR[7])+0.3535533905932737*fCC[10]+0.1767766952966368*(fBR[10]+fBL[10]+fCL[7]+fBL[7])-0.3535533905932737*fBC[10]+0.15309310892394856*(fCR[3]+fCL[3]+fBR[3]+fBL[3])-0.3061862178478971*(fCC[3]+fBC[3]); 
  df_rec_lo[7] = 0.20412414523193148*(fCR[15]+fBL[15])-0.20412414523193148*(fCL[15]+fBR[15])-0.1767766952966368*(fCR[13]+fCL[13]+fCR[11]+fBR[11])+0.3535533905932737*fCC[13]+0.1767766952966368*(fBR[13]+fBL[13]+fCL[11]+fBL[11])-0.3535533905932737*fBC[13]+0.15309310892394856*(fCR[6]+fCL[6]+fBR[6]+fBL[6])-0.3061862178478971*(fCC[6]+fBC[6]); 
  df_rec_up[0] = 0.11785113019775789*(fTR[9]+fTL[9])-0.2357022603955158*fTC[9]-0.11785113019775789*(fCR[9]+fCL[9])+0.2357022603955158*fCC[9]-0.10206207261596573*(fTR[4]+fCL[4]+fTR[2]+fTL[2]+fCR[2]+fCL[2])+0.10206207261596573*(fTL[4]+fCR[4])+0.20412414523193148*(fTC[2]+fCC[2])+0.0883883476483184*(fTR[0]+fCR[0])-0.0883883476483184*(fTL[0]+fCL[0]); 
  df_rec_up[1] = 0.11785113019775789*(fTR[12]+fTL[12])-0.2357022603955158*fTC[12]-0.11785113019775789*(fCR[12]+fCL[12])+0.2357022603955158*fCC[12]-0.10206207261596573*(fTR[8]+fCL[8]+fTR[5]+fTL[5]+fCR[5]+fCL[5])+0.10206207261596573*(fTL[8]+fCR[8])+0.20412414523193148*(fTC[5]+fCC[5])+0.0883883476483184*(fTR[1]+fCR[1])-0.0883883476483184*(fTL[1]+fCL[1]); 
  df_rec_up[2] = 0.20412414523193148*(fTR[9]+fCL[9])-0.20412414523193148*(fTL[9]+fCR[9])-0.1767766952966368*(fTR[4]+fTL[4]+fTR[2]+fCR[2])+0.3535533905932737*fTC[4]+0.1767766952966368*(fCR[4]+fCL[4]+fTL[2]+fCL[2])-0.3535533905932737*fCC[4]+0.15309310892394856*(fTR[0]+fTL[0]+fCR[0]+fCL[0])-0.3061862178478971*(fTC[0]+fCC[0]); 
  df_rec_up[3] = 0.11785113019775789*(fTR[14]+fTL[14])-0.2357022603955158*fTC[14]-0.11785113019775789*(fCR[14]+fCL[14])+0.2357022603955158*fCC[14]-0.10206207261596573*(fTR[10]+fCL[10]+fTR[7]+fTL[7]+fCR[7]+fCL[7])+0.10206207261596573*(fTL[10]+fCR[10])+0.20412414523193148*(fTC[7]+fCC[7])+0.0883883476483184*(fTR[3]+fCR[3])-0.0883883476483184*(fTL[3]+fCL[3]); 
  df_rec_up[4] = 0.20412414523193148*(fTR[12]+fCL[12])-0.20412414523193148*(fTL[12]+fCR[12])-0.1767766952966368*(fTR[8]+fTL[8]+fTR[5]+fCR[5])+0.3535533905932737*fTC[8]+0.1767766952966368*(fCR[8]+fCL[8]+fTL[5]+fCL[5])-0.3535533905932737*fCC[8]+0.15309310892394856*(fTR[1]+fTL[1]+fCR[1]+fCL[1])-0.3061862178478971*(fTC[1]+fCC[1]); 
  df_rec_up[5] = 0.11785113019775789*(fTR[15]+fTL[15])-0.2357022603955158*fTC[15]-0.11785113019775789*(fCR[15]+fCL[15])+0.2357022603955158*fCC[15]-0.10206207261596573*(fTR[13]+fCL[13]+fTR[11]+fTL[11]+fCR[11]+fCL[11])+0.10206207261596573*(fTL[13]+fCR[13])+0.20412414523193148*(fTC[11]+fCC[11])+0.0883883476483184*(fTR[6]+fCR[6])-0.0883883476483184*(fTL[6]+fCL[6]); 
  df_rec_up[6] = 0.20412414523193148*(fTR[14]+fCL[14])-0.20412414523193148*(fTL[14]+fCR[14])-0.1767766952966368*(fTR[10]+fTL[10]+fTR[7]+fCR[7])+0.3535533905932737*fTC[10]+0.1767766952966368*(fCR[10]+fCL[10]+fTL[7]+fCL[7])-0.3535533905932737*fCC[10]+0.15309310892394856*(fTR[3]+fTL[3]+fCR[3]+fCL[3])-0.3061862178478971*(fTC[3]+fCC[3]); 
  df_rec_up[7] = 0.20412414523193148*(fTR[15]+fCL[15])-0.20412414523193148*(fTL[15]+fCR[15])-0.1767766952966368*(fTR[13]+fTL[13]+fTR[11]+fCR[11])+0.3535533905932737*fTC[13]+0.1767766952966368*(fCR[13]+fCL[13]+fTL[11]+fCL[11])-0.3535533905932737*fCC[13]+0.15309310892394856*(fTR[6]+fTL[6]+fCR[6]+fCL[6])-0.3061862178478971*(fTC[6]+fCC[6]); 

  surft1_lo[0] = 0.3535533905932737*(Dsurf_CC_vz[7]*df_rec_lo[7]+Dsurf_CC_vz[6]*df_rec_lo[6]+Dsurf_CC_vz[5]*df_rec_lo[5]+Dsurf_CC_vz[4]*df_rec_lo[4]+Dsurf_CC_vz[3]*df_rec_lo[3]+Dsurf_CC_vz[2]*df_rec_lo[2]+Dsurf_CC_vz[1]*df_rec_lo[1]+Dsurf_CC_vz[0]*df_rec_lo[0]); 
  surft1_lo[1] = 0.3535533905932737*(Dsurf_CC_vz[6]*df_rec_lo[7]+df_rec_lo[6]*Dsurf_CC_vz[7]+Dsurf_CC_vz[3]*df_rec_lo[5]+df_rec_lo[3]*Dsurf_CC_vz[5]+Dsurf_CC_vz[2]*df_rec_lo[4]+df_rec_lo[2]*Dsurf_CC_vz[4]+Dsurf_CC_vz[0]*df_rec_lo[1]+df_rec_lo[0]*Dsurf_CC_vz[1]); 
  surft1_lo[2] = 0.3535533905932737*(Dsurf_CC_vz[5]*df_rec_lo[7]+df_rec_lo[5]*Dsurf_CC_vz[7]+Dsurf_CC_vz[3]*df_rec_lo[6]+df_rec_lo[3]*Dsurf_CC_vz[6]+Dsurf_CC_vz[1]*df_rec_lo[4]+df_rec_lo[1]*Dsurf_CC_vz[4]+Dsurf_CC_vz[0]*df_rec_lo[2]+df_rec_lo[0]*Dsurf_CC_vz[2]); 
  surft1_lo[3] = 0.3535533905932737*(Dsurf_CC_vz[4]*df_rec_lo[7]+df_rec_lo[4]*Dsurf_CC_vz[7]+Dsurf_CC_vz[2]*df_rec_lo[6]+df_rec_lo[2]*Dsurf_CC_vz[6]+Dsurf_CC_vz[1]*df_rec_lo[5]+df_rec_lo[1]*Dsurf_CC_vz[5]+Dsurf_CC_vz[0]*df_rec_lo[3]+df_rec_lo[0]*Dsurf_CC_vz[3]); 
  surft1_lo[4] = 0.3535533905932737*(Dsurf_CC_vz[3]*df_rec_lo[7]+df_rec_lo[3]*Dsurf_CC_vz[7]+Dsurf_CC_vz[5]*df_rec_lo[6]+df_rec_lo[5]*Dsurf_CC_vz[6]+Dsurf_CC_vz[0]*df_rec_lo[4]+df_rec_lo[0]*Dsurf_CC_vz[4]+Dsurf_CC_vz[1]*df_rec_lo[2]+df_rec_lo[1]*Dsurf_CC_vz[2]); 
  surft1_lo[5] = 0.3535533905932737*(Dsurf_CC_vz[2]*df_rec_lo[7]+df_rec_lo[2]*Dsurf_CC_vz[7]+Dsurf_CC_vz[4]*df_rec_lo[6]+df_rec_lo[4]*Dsurf_CC_vz[6]+Dsurf_CC_vz[0]*df_rec_lo[5]+df_rec_lo[0]*Dsurf_CC_vz[5]+Dsurf_CC_vz[1]*df_rec_lo[3]+df_rec_lo[1]*Dsurf_CC_vz[3]); 
  surft1_lo[6] = 0.3535533905932737*(Dsurf_CC_vz[1]*df_rec_lo[7]+df_rec_lo[1]*Dsurf_CC_vz[7]+Dsurf_CC_vz[0]*df_rec_lo[6]+df_rec_lo[0]*Dsurf_CC_vz[6]+Dsurf_CC_vz[4]*df_rec_lo[5]+df_rec_lo[4]*Dsurf_CC_vz[5]+Dsurf_CC_vz[2]*df_rec_lo[3]+df_rec_lo[2]*Dsurf_CC_vz[3]); 
  surft1_lo[7] = 0.3535533905932737*(Dsurf_CC_vz[0]*df_rec_lo[7]+df_rec_lo[0]*Dsurf_CC_vz[7]+Dsurf_CC_vz[1]*df_rec_lo[6]+df_rec_lo[1]*Dsurf_CC_vz[6]+Dsurf_CC_vz[2]*df_rec_lo[5]+df_rec_lo[2]*Dsurf_CC_vz[5]+Dsurf_CC_vz[3]*df_rec_lo[4]+df_rec_lo[3]*Dsurf_CC_vz[4]); 
  surft1_up[0] = 0.3535533905932737*(Dsurf_TC_vz[7]*df_rec_up[7]+Dsurf_TC_vz[6]*df_rec_up[6]+Dsurf_TC_vz[5]*df_rec_up[5]+Dsurf_TC_vz[4]*df_rec_up[4]+Dsurf_TC_vz[3]*df_rec_up[3]+Dsurf_TC_vz[2]*df_rec_up[2]+Dsurf_TC_vz[1]*df_rec_up[1]+Dsurf_TC_vz[0]*df_rec_up[0]); 
  surft1_up[1] = 0.3535533905932737*(Dsurf_TC_vz[6]*df_rec_up[7]+df_rec_up[6]*Dsurf_TC_vz[7]+Dsurf_TC_vz[3]*df_rec_up[5]+df_rec_up[3]*Dsurf_TC_vz[5]+Dsurf_TC_vz[2]*df_rec_up[4]+df_rec_up[2]*Dsurf_TC_vz[4]+Dsurf_TC_vz[0]*df_rec_up[1]+df_rec_up[0]*Dsurf_TC_vz[1]); 
  surft1_up[2] = 0.3535533905932737*(Dsurf_TC_vz[5]*df_rec_up[7]+df_rec_up[5]*Dsurf_TC_vz[7]+Dsurf_TC_vz[3]*df_rec_up[6]+df_rec_up[3]*Dsurf_TC_vz[6]+Dsurf_TC_vz[1]*df_rec_up[4]+df_rec_up[1]*Dsurf_TC_vz[4]+Dsurf_TC_vz[0]*df_rec_up[2]+df_rec_up[0]*Dsurf_TC_vz[2]); 
  surft1_up[3] = 0.3535533905932737*(Dsurf_TC_vz[4]*df_rec_up[7]+df_rec_up[4]*Dsurf_TC_vz[7]+Dsurf_TC_vz[2]*df_rec_up[6]+df_rec_up[2]*Dsurf_TC_vz[6]+Dsurf_TC_vz[1]*df_rec_up[5]+df_rec_up[1]*Dsurf_TC_vz[5]+Dsurf_TC_vz[0]*df_rec_up[3]+df_rec_up[0]*Dsurf_TC_vz[3]); 
  surft1_up[4] = 0.3535533905932737*(Dsurf_TC_vz[3]*df_rec_up[7]+df_rec_up[3]*Dsurf_TC_vz[7]+Dsurf_TC_vz[5]*df_rec_up[6]+df_rec_up[5]*Dsurf_TC_vz[6]+Dsurf_TC_vz[0]*df_rec_up[4]+df_rec_up[0]*Dsurf_TC_vz[4]+Dsurf_TC_vz[1]*df_rec_up[2]+df_rec_up[1]*Dsurf_TC_vz[2]); 
  surft1_up[5] = 0.3535533905932737*(Dsurf_TC_vz[2]*df_rec_up[7]+df_rec_up[2]*Dsurf_TC_vz[7]+Dsurf_TC_vz[4]*df_rec_up[6]+df_rec_up[4]*Dsurf_TC_vz[6]+Dsurf_TC_vz[0]*df_rec_up[5]+df_rec_up[0]*Dsurf_TC_vz[5]+Dsurf_TC_vz[1]*df_rec_up[3]+df_rec_up[1]*Dsurf_TC_vz[3]); 
  surft1_up[6] = 0.3535533905932737*(Dsurf_TC_vz[1]*df_rec_up[7]+df_rec_up[1]*Dsurf_TC_vz[7]+Dsurf_TC_vz[0]*df_rec_up[6]+df_rec_up[0]*Dsurf_TC_vz[6]+Dsurf_TC_vz[4]*df_rec_up[5]+df_rec_up[4]*Dsurf_TC_vz[5]+Dsurf_TC_vz[2]*df_rec_up[3]+df_rec_up[2]*Dsurf_TC_vz[3]); 
  surft1_up[7] = 0.3535533905932737*(Dsurf_TC_vz[0]*df_rec_up[7]+df_rec_up[0]*Dsurf_TC_vz[7]+Dsurf_TC_vz[1]*df_rec_up[6]+df_rec_up[1]*Dsurf_TC_vz[6]+Dsurf_TC_vz[2]*df_rec_up[5]+df_rec_up[2]*Dsurf_TC_vz[5]+Dsurf_TC_vz[3]*df_rec_up[4]+df_rec_up[3]*Dsurf_TC_vz[4]); 

  surft2_lo[0] = 0.3535533905932737*(Dsurf_CC_vx[7]*f_rec_lo[7]+Dsurf_CC_vx[6]*f_rec_lo[6]+Dsurf_CC_vx[5]*f_rec_lo[5]+Dsurf_CC_vx[4]*f_rec_lo[4]+Dsurf_CC_vx[3]*f_rec_lo[3]+Dsurf_CC_vx[2]*f_rec_lo[2]+Dsurf_CC_vx[1]*f_rec_lo[1]+Dsurf_CC_vx[0]*f_rec_lo[0]); 
  surft2_lo[1] = 0.3535533905932737*(Dsurf_CC_vx[6]*f_rec_lo[7]+f_rec_lo[6]*Dsurf_CC_vx[7]+Dsurf_CC_vx[3]*f_rec_lo[5]+f_rec_lo[3]*Dsurf_CC_vx[5]+Dsurf_CC_vx[2]*f_rec_lo[4]+f_rec_lo[2]*Dsurf_CC_vx[4]+Dsurf_CC_vx[0]*f_rec_lo[1]+f_rec_lo[0]*Dsurf_CC_vx[1]); 
  surft2_lo[2] = 0.3535533905932737*(Dsurf_CC_vx[5]*f_rec_lo[7]+f_rec_lo[5]*Dsurf_CC_vx[7]+Dsurf_CC_vx[3]*f_rec_lo[6]+f_rec_lo[3]*Dsurf_CC_vx[6]+Dsurf_CC_vx[1]*f_rec_lo[4]+f_rec_lo[1]*Dsurf_CC_vx[4]+Dsurf_CC_vx[0]*f_rec_lo[2]+f_rec_lo[0]*Dsurf_CC_vx[2]); 
  surft2_lo[3] = 0.3535533905932737*(Dsurf_CC_vx[4]*f_rec_lo[7]+f_rec_lo[4]*Dsurf_CC_vx[7]+Dsurf_CC_vx[2]*f_rec_lo[6]+f_rec_lo[2]*Dsurf_CC_vx[6]+Dsurf_CC_vx[1]*f_rec_lo[5]+f_rec_lo[1]*Dsurf_CC_vx[5]+Dsurf_CC_vx[0]*f_rec_lo[3]+f_rec_lo[0]*Dsurf_CC_vx[3]); 
  surft2_lo[4] = 0.3535533905932737*(Dsurf_CC_vx[3]*f_rec_lo[7]+f_rec_lo[3]*Dsurf_CC_vx[7]+Dsurf_CC_vx[5]*f_rec_lo[6]+f_rec_lo[5]*Dsurf_CC_vx[6]+Dsurf_CC_vx[0]*f_rec_lo[4]+f_rec_lo[0]*Dsurf_CC_vx[4]+Dsurf_CC_vx[1]*f_rec_lo[2]+f_rec_lo[1]*Dsurf_CC_vx[2]); 
  surft2_lo[5] = 0.3535533905932737*(Dsurf_CC_vx[2]*f_rec_lo[7]+f_rec_lo[2]*Dsurf_CC_vx[7]+Dsurf_CC_vx[4]*f_rec_lo[6]+f_rec_lo[4]*Dsurf_CC_vx[6]+Dsurf_CC_vx[0]*f_rec_lo[5]+f_rec_lo[0]*Dsurf_CC_vx[5]+Dsurf_CC_vx[1]*f_rec_lo[3]+f_rec_lo[1]*Dsurf_CC_vx[3]); 
  surft2_lo[6] = 0.3535533905932737*(Dsurf_CC_vx[1]*f_rec_lo[7]+f_rec_lo[1]*Dsurf_CC_vx[7]+Dsurf_CC_vx[0]*f_rec_lo[6]+f_rec_lo[0]*Dsurf_CC_vx[6]+Dsurf_CC_vx[4]*f_rec_lo[5]+f_rec_lo[4]*Dsurf_CC_vx[5]+Dsurf_CC_vx[2]*f_rec_lo[3]+f_rec_lo[2]*Dsurf_CC_vx[3]); 
  surft2_lo[7] = 0.3535533905932737*(Dsurf_CC_vx[0]*f_rec_lo[7]+f_rec_lo[0]*Dsurf_CC_vx[7]+Dsurf_CC_vx[1]*f_rec_lo[6]+f_rec_lo[1]*Dsurf_CC_vx[6]+Dsurf_CC_vx[2]*f_rec_lo[5]+f_rec_lo[2]*Dsurf_CC_vx[5]+Dsurf_CC_vx[3]*f_rec_lo[4]+f_rec_lo[3]*Dsurf_CC_vx[4]); 
  surft2_up[0] = 0.3535533905932737*(Dsurf_CR_vx[7]*f_rec_up[7]+Dsurf_CR_vx[6]*f_rec_up[6]+Dsurf_CR_vx[5]*f_rec_up[5]+Dsurf_CR_vx[4]*f_rec_up[4]+Dsurf_CR_vx[3]*f_rec_up[3]+Dsurf_CR_vx[2]*f_rec_up[2]+Dsurf_CR_vx[1]*f_rec_up[1]+Dsurf_CR_vx[0]*f_rec_up[0]); 
  surft2_up[1] = 0.3535533905932737*(Dsurf_CR_vx[6]*f_rec_up[7]+f_rec_up[6]*Dsurf_CR_vx[7]+Dsurf_CR_vx[3]*f_rec_up[5]+f_rec_up[3]*Dsurf_CR_vx[5]+Dsurf_CR_vx[2]*f_rec_up[4]+f_rec_up[2]*Dsurf_CR_vx[4]+Dsurf_CR_vx[0]*f_rec_up[1]+f_rec_up[0]*Dsurf_CR_vx[1]); 
  surft2_up[2] = 0.3535533905932737*(Dsurf_CR_vx[5]*f_rec_up[7]+f_rec_up[5]*Dsurf_CR_vx[7]+Dsurf_CR_vx[3]*f_rec_up[6]+f_rec_up[3]*Dsurf_CR_vx[6]+Dsurf_CR_vx[1]*f_rec_up[4]+f_rec_up[1]*Dsurf_CR_vx[4]+Dsurf_CR_vx[0]*f_rec_up[2]+f_rec_up[0]*Dsurf_CR_vx[2]); 
  surft2_up[3] = 0.3535533905932737*(Dsurf_CR_vx[4]*f_rec_up[7]+f_rec_up[4]*Dsurf_CR_vx[7]+Dsurf_CR_vx[2]*f_rec_up[6]+f_rec_up[2]*Dsurf_CR_vx[6]+Dsurf_CR_vx[1]*f_rec_up[5]+f_rec_up[1]*Dsurf_CR_vx[5]+Dsurf_CR_vx[0]*f_rec_up[3]+f_rec_up[0]*Dsurf_CR_vx[3]); 
  surft2_up[4] = 0.3535533905932737*(Dsurf_CR_vx[3]*f_rec_up[7]+f_rec_up[3]*Dsurf_CR_vx[7]+Dsurf_CR_vx[5]*f_rec_up[6]+f_rec_up[5]*Dsurf_CR_vx[6]+Dsurf_CR_vx[0]*f_rec_up[4]+f_rec_up[0]*Dsurf_CR_vx[4]+Dsurf_CR_vx[1]*f_rec_up[2]+f_rec_up[1]*Dsurf_CR_vx[2]); 
  surft2_up[5] = 0.3535533905932737*(Dsurf_CR_vx[2]*f_rec_up[7]+f_rec_up[2]*Dsurf_CR_vx[7]+Dsurf_CR_vx[4]*f_rec_up[6]+f_rec_up[4]*Dsurf_CR_vx[6]+Dsurf_CR_vx[0]*f_rec_up[5]+f_rec_up[0]*Dsurf_CR_vx[5]+Dsurf_CR_vx[1]*f_rec_up[3]+f_rec_up[1]*Dsurf_CR_vx[3]); 
  surft2_up[6] = 0.3535533905932737*(Dsurf_CR_vx[1]*f_rec_up[7]+f_rec_up[1]*Dsurf_CR_vx[7]+Dsurf_CR_vx[0]*f_rec_up[6]+f_rec_up[0]*Dsurf_CR_vx[6]+Dsurf_CR_vx[4]*f_rec_up[5]+f_rec_up[4]*Dsurf_CR_vx[5]+Dsurf_CR_vx[2]*f_rec_up[3]+f_rec_up[2]*Dsurf_CR_vx[3]); 
  surft2_up[7] = 0.3535533905932737*(Dsurf_CR_vx[0]*f_rec_up[7]+f_rec_up[0]*Dsurf_CR_vx[7]+Dsurf_CR_vx[1]*f_rec_up[6]+f_rec_up[1]*Dsurf_CR_vx[6]+Dsurf_CR_vx[2]*f_rec_up[5]+f_rec_up[2]*Dsurf_CR_vx[5]+Dsurf_CR_vx[3]*f_rec_up[4]+f_rec_up[3]*Dsurf_CR_vx[4]); 

  vol[2] = 0.75*(fCC[11]*DCC[15]+fCC[7]*DCC[14]+fCC[6]*DCC[13]+fCC[5]*DCC[12]+fCC[3]*DCC[10]+fCC[2]*DCC[9]+fCC[1]*DCC[8]+fCC[0]*DCC[4]); 
  vol[5] = 0.75*(fCC[7]*DCC[15]+fCC[11]*DCC[14]+fCC[3]*DCC[13]+fCC[2]*DCC[12]+fCC[6]*DCC[10]+fCC[5]*DCC[9]+fCC[0]*DCC[8]+fCC[1]*DCC[4]); 
  vol[7] = 0.75*(fCC[5]*DCC[15]+fCC[2]*DCC[14]+fCC[1]*DCC[13]+fCC[11]*DCC[12]+fCC[0]*DCC[10]+fCC[7]*DCC[9]+fCC[6]*DCC[8]+fCC[3]*DCC[4]); 
  vol[9] = 1.5*(DCC[15]*fCC[15]+DCC[14]*fCC[14]+DCC[13]*fCC[13]+DCC[12]*fCC[12]+DCC[10]*fCC[10]+DCC[9]*fCC[9]+DCC[8]*fCC[8]+DCC[4]*fCC[4])+0.75*(DCC[11]*fCC[11]+DCC[7]*fCC[7]+DCC[6]*fCC[6]+DCC[5]*fCC[5]+DCC[3]*fCC[3]+DCC[2]*fCC[2]+DCC[1]*fCC[1]+DCC[0]*fCC[0]); 
  vol[11] = 0.75*(fCC[2]*DCC[15]+fCC[5]*DCC[14]+fCC[0]*DCC[13]+fCC[7]*DCC[12]+DCC[9]*fCC[11]+fCC[1]*DCC[10]+fCC[3]*DCC[8]+DCC[4]*fCC[6]); 
  vol[12] = 1.5*(DCC[14]*fCC[15]+fCC[14]*DCC[15]+DCC[10]*fCC[13]+fCC[10]*DCC[13]+DCC[9]*fCC[12]+fCC[9]*DCC[12]+DCC[4]*fCC[8]+fCC[4]*DCC[8])+0.75*(DCC[7]*fCC[11]+fCC[7]*DCC[11]+DCC[3]*fCC[6]+fCC[3]*DCC[6]+DCC[2]*fCC[5]+fCC[2]*DCC[5]+DCC[0]*fCC[1]+fCC[0]*DCC[1]); 
  vol[14] = 1.5*(DCC[12]*fCC[15]+fCC[12]*DCC[15]+DCC[9]*fCC[14]+fCC[9]*DCC[14]+DCC[8]*fCC[13]+fCC[8]*DCC[13]+DCC[4]*fCC[10]+fCC[4]*DCC[10])+0.75*(DCC[5]*fCC[11]+fCC[5]*DCC[11]+DCC[2]*fCC[7]+fCC[2]*DCC[7]+DCC[1]*fCC[6]+fCC[1]*DCC[6]+DCC[0]*fCC[3]+fCC[0]*DCC[3]); 
  vol[15] = 1.5*(DCC[9]*fCC[15]+fCC[9]*DCC[15]+DCC[12]*fCC[14]+fCC[12]*DCC[14]+DCC[4]*fCC[13]+fCC[4]*DCC[13]+DCC[8]*fCC[10]+fCC[8]*DCC[10])+0.75*(DCC[2]*fCC[11]+fCC[2]*DCC[11]+DCC[5]*fCC[7]+fCC[5]*DCC[7]+DCC[0]*fCC[6]+fCC[0]*DCC[6]+DCC[1]*fCC[3]+fCC[1]*DCC[3]); 

  out[0] += (0.5*vol[0]+0.35355339059327373*surft1_up[0]-0.35355339059327373*surft1_lo[0])*dv_inv_sq; 
  out[1] += (0.5*vol[1]+0.35355339059327373*surft1_up[1]-0.35355339059327373*surft1_lo[1])*dv_inv_sq; 
  out[2] += (0.5*vol[2]+0.35355339059327373*surft1_up[2]-0.35355339059327373*surft1_lo[2])*dv_inv_sq; 
  out[3] += (0.5*vol[3]+0.35355339059327373*surft1_up[3]-0.35355339059327373*surft1_lo[3])*dv_inv_sq; 
  out[4] += (0.5*vol[4]-0.6123724356957945*surft2_up[0]+0.6123724356957945*(surft2_lo[0]+surft1_up[0]+surft1_lo[0]))*dv_inv_sq; 
  out[5] += (0.5*vol[5]+0.35355339059327373*surft1_up[4]-0.35355339059327373*surft1_lo[4])*dv_inv_sq; 
  out[6] += (0.5*vol[6]+0.35355339059327373*surft1_up[5]-0.35355339059327373*surft1_lo[5])*dv_inv_sq; 
  out[7] += (0.5*vol[7]+0.35355339059327373*surft1_up[6]-0.35355339059327373*surft1_lo[6])*dv_inv_sq; 
  out[8] += (0.5*vol[8]-0.6123724356957945*surft2_up[1]+0.6123724356957945*(surft2_lo[1]+surft1_up[1]+surft1_lo[1]))*dv_inv_sq; 
  out[9] += (0.5*vol[9]+0.6123724356957945*(surft1_up[2]+surft1_lo[2])-1.0606601717798212*(surft2_up[0]+surft2_lo[0]))*dv_inv_sq; 
  out[10] += (0.5*vol[10]+0.6123724356957945*(surft1_up[3]+surft1_lo[3]+surft2_lo[2])-0.6123724356957945*surft2_up[2])*dv_inv_sq; 
  out[11] += (0.5*vol[11]+0.35355339059327373*surft1_up[7]-0.35355339059327373*surft1_lo[7])*dv_inv_sq; 
  out[12] += (0.5*vol[12]+0.6123724356957945*(surft1_up[4]+surft1_lo[4])-1.0606601717798212*(surft2_up[1]+surft2_lo[1]))*dv_inv_sq; 
  out[13] += (0.5*vol[13]+0.6123724356957945*(surft1_up[5]+surft1_lo[5]+surft2_lo[4])-0.6123724356957945*surft2_up[4])*dv_inv_sq; 
  out[14] += (0.5*vol[14]+0.6123724356957945*(surft1_up[6]+surft1_lo[6])-1.0606601717798212*(surft2_up[2]+surft2_lo[2]))*dv_inv_sq; 
  out[15] += (0.5*vol[15]+0.6123724356957945*(surft1_up[7]+surft1_lo[7])-1.0606601717798212*(surft2_up[4]+surft2_lo[4]))*dv_inv_sq; 

  double cflFreq = fmax(fabs(Dsurf_CC_vz[0]), fabs(Dsurf_TC_vz[0])); 

  return 1.4142135623730947*dv_inv_sq*cflFreq; 
} 
