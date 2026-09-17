#include <gkyl_fpo_vlasov_kernels.h> 

GKYL_CU_DH double fpo_vlasov_diff_surfvzvz_1x3v_ser_p1(const double* dxv, const double* diff_coeff_C, 
  const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* GKYL_RESTRICT out) { 
  // dxv[NDIM]: Cell spacing in each direction. 
  // diff_coeff_C: Diffusion tensor in center cell. 
  // diff_coeff_surf_stencil[9]: 9-cell stencil of surface expansion of recovered diffusion tensor. 
  // f_stencil[9]: 9-cell stencil of distribution function. 
  // out: Incremented output. 

  double dv_inv_sq = 4.0/dxv[3]/dxv[3]; 

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
  const double* fL = f_stencil[0]; 
  const double* fC = f_stencil[1]; 
  const double* fR = f_stencil[2]; 

  const double* DC = &diff_coeff_C[128]; 
  const double* D_surf_C = &diff_coeff_surf_stencil[1][128]; 
  const double* D_surf_R = &diff_coeff_surf_stencil[2][128]; 

  f_rec_lo[0] = 0.408248290463863*fL[4]-0.408248290463863*fC[4]+0.3535533905932737*(fL[0]+fC[0]); 
  f_rec_lo[1] = 0.408248290463863*fL[8]-0.408248290463863*fC[8]+0.3535533905932737*(fL[1]+fC[1]); 
  f_rec_lo[2] = 0.408248290463863*fL[9]-0.408248290463863*fC[9]+0.3535533905932737*(fL[2]+fC[2]); 
  f_rec_lo[3] = 0.408248290463863*fL[10]-0.408248290463863*fC[10]+0.3535533905932737*(fL[3]+fC[3]); 
  f_rec_lo[4] = 0.408248290463863*fL[12]-0.408248290463863*fC[12]+0.3535533905932737*(fL[5]+fC[5]); 
  f_rec_lo[5] = 0.408248290463863*fL[13]-0.408248290463863*fC[13]+0.3535533905932737*(fL[6]+fC[6]); 
  f_rec_lo[6] = 0.408248290463863*fL[14]-0.408248290463863*fC[14]+0.3535533905932737*(fL[7]+fC[7]); 
  f_rec_lo[7] = 0.408248290463863*fL[15]-0.408248290463863*fC[15]+0.3535533905932737*(fL[11]+fC[11]); 
  f_rec_up[0] = -(0.408248290463863*fR[4])+0.408248290463863*fC[4]+0.3535533905932737*(fR[0]+fC[0]); 
  f_rec_up[1] = -(0.408248290463863*fR[8])+0.408248290463863*fC[8]+0.3535533905932737*(fR[1]+fC[1]); 
  f_rec_up[2] = -(0.408248290463863*fR[9])+0.408248290463863*fC[9]+0.3535533905932737*(fR[2]+fC[2]); 
  f_rec_up[3] = -(0.408248290463863*fR[10])+0.408248290463863*fC[10]+0.3535533905932737*(fR[3]+fC[3]); 
  f_rec_up[4] = -(0.408248290463863*fR[12])+0.408248290463863*fC[12]+0.3535533905932737*(fR[5]+fC[5]); 
  f_rec_up[5] = -(0.408248290463863*fR[13])+0.408248290463863*fC[13]+0.3535533905932737*(fR[6]+fC[6]); 
  f_rec_up[6] = -(0.408248290463863*fR[14])+0.408248290463863*fC[14]+0.3535533905932737*(fR[7]+fC[7]); 
  f_rec_up[7] = -(0.408248290463863*fR[15])+0.408248290463863*fC[15]+0.3535533905932737*(fR[11]+fC[11]); 

  df_rec_lo[0] = -(0.7654655446197428*(fL[4]+fC[4]))-0.7954951288348656*fL[0]+0.7954951288348656*fC[0]; 
  df_rec_lo[1] = -(0.7654655446197428*(fL[8]+fC[8]))-0.7954951288348656*fL[1]+0.7954951288348656*fC[1]; 
  df_rec_lo[2] = -(0.7654655446197428*(fL[9]+fC[9]))-0.7954951288348656*fL[2]+0.7954951288348656*fC[2]; 
  df_rec_lo[3] = -(0.7654655446197428*(fL[10]+fC[10]))-0.7954951288348656*fL[3]+0.7954951288348656*fC[3]; 
  df_rec_lo[4] = -(0.7654655446197428*(fL[12]+fC[12]))-0.7954951288348656*fL[5]+0.7954951288348656*fC[5]; 
  df_rec_lo[5] = -(0.7654655446197428*(fL[13]+fC[13]))-0.7954951288348656*fL[6]+0.7954951288348656*fC[6]; 
  df_rec_lo[6] = -(0.7654655446197428*(fL[14]+fC[14]))-0.7954951288348656*fL[7]+0.7954951288348656*fC[7]; 
  df_rec_lo[7] = -(0.7654655446197428*(fL[15]+fC[15]))-0.7954951288348656*fL[11]+0.7954951288348656*fC[11]; 
  df_rec_up[0] = -(0.7654655446197428*(fR[4]+fC[4]))+0.7954951288348656*fR[0]-0.7954951288348656*fC[0]; 
  df_rec_up[1] = -(0.7654655446197428*(fR[8]+fC[8]))+0.7954951288348656*fR[1]-0.7954951288348656*fC[1]; 
  df_rec_up[2] = -(0.7654655446197428*(fR[9]+fC[9]))+0.7954951288348656*fR[2]-0.7954951288348656*fC[2]; 
  df_rec_up[3] = -(0.7654655446197428*(fR[10]+fC[10]))+0.7954951288348656*fR[3]-0.7954951288348656*fC[3]; 
  df_rec_up[4] = -(0.7654655446197428*(fR[12]+fC[12]))+0.7954951288348656*fR[5]-0.7954951288348656*fC[5]; 
  df_rec_up[5] = -(0.7654655446197428*(fR[13]+fC[13]))+0.7954951288348656*fR[6]-0.7954951288348656*fC[6]; 
  df_rec_up[6] = -(0.7654655446197428*(fR[14]+fC[14]))+0.7954951288348656*fR[7]-0.7954951288348656*fC[7]; 
  df_rec_up[7] = -(0.7654655446197428*(fR[15]+fC[15]))+0.7954951288348656*fR[11]-0.7954951288348656*fC[11]; 

  surft1_lo[0] = 0.3535533905932737*(D_surf_C[7]*df_rec_lo[7]+D_surf_C[6]*df_rec_lo[6]+D_surf_C[5]*df_rec_lo[5]+D_surf_C[4]*df_rec_lo[4]+D_surf_C[3]*df_rec_lo[3]+D_surf_C[2]*df_rec_lo[2]+D_surf_C[1]*df_rec_lo[1]+D_surf_C[0]*df_rec_lo[0]); 
  surft1_lo[1] = 0.3535533905932737*(D_surf_C[6]*df_rec_lo[7]+df_rec_lo[6]*D_surf_C[7]+D_surf_C[3]*df_rec_lo[5]+df_rec_lo[3]*D_surf_C[5]+D_surf_C[2]*df_rec_lo[4]+df_rec_lo[2]*D_surf_C[4]+D_surf_C[0]*df_rec_lo[1]+df_rec_lo[0]*D_surf_C[1]); 
  surft1_lo[2] = 0.3535533905932737*(D_surf_C[5]*df_rec_lo[7]+df_rec_lo[5]*D_surf_C[7]+D_surf_C[3]*df_rec_lo[6]+df_rec_lo[3]*D_surf_C[6]+D_surf_C[1]*df_rec_lo[4]+df_rec_lo[1]*D_surf_C[4]+D_surf_C[0]*df_rec_lo[2]+df_rec_lo[0]*D_surf_C[2]); 
  surft1_lo[3] = 0.3535533905932737*(D_surf_C[4]*df_rec_lo[7]+df_rec_lo[4]*D_surf_C[7]+D_surf_C[2]*df_rec_lo[6]+df_rec_lo[2]*D_surf_C[6]+D_surf_C[1]*df_rec_lo[5]+df_rec_lo[1]*D_surf_C[5]+D_surf_C[0]*df_rec_lo[3]+df_rec_lo[0]*D_surf_C[3]); 
  surft1_lo[4] = 0.3535533905932737*(D_surf_C[3]*df_rec_lo[7]+df_rec_lo[3]*D_surf_C[7]+D_surf_C[5]*df_rec_lo[6]+df_rec_lo[5]*D_surf_C[6]+D_surf_C[0]*df_rec_lo[4]+df_rec_lo[0]*D_surf_C[4]+D_surf_C[1]*df_rec_lo[2]+df_rec_lo[1]*D_surf_C[2]); 
  surft1_lo[5] = 0.3535533905932737*(D_surf_C[2]*df_rec_lo[7]+df_rec_lo[2]*D_surf_C[7]+D_surf_C[4]*df_rec_lo[6]+df_rec_lo[4]*D_surf_C[6]+D_surf_C[0]*df_rec_lo[5]+df_rec_lo[0]*D_surf_C[5]+D_surf_C[1]*df_rec_lo[3]+df_rec_lo[1]*D_surf_C[3]); 
  surft1_lo[6] = 0.3535533905932737*(D_surf_C[1]*df_rec_lo[7]+df_rec_lo[1]*D_surf_C[7]+D_surf_C[0]*df_rec_lo[6]+df_rec_lo[0]*D_surf_C[6]+D_surf_C[4]*df_rec_lo[5]+df_rec_lo[4]*D_surf_C[5]+D_surf_C[2]*df_rec_lo[3]+df_rec_lo[2]*D_surf_C[3]); 
  surft1_lo[7] = 0.3535533905932737*(D_surf_C[0]*df_rec_lo[7]+df_rec_lo[0]*D_surf_C[7]+D_surf_C[1]*df_rec_lo[6]+df_rec_lo[1]*D_surf_C[6]+D_surf_C[2]*df_rec_lo[5]+df_rec_lo[2]*D_surf_C[5]+D_surf_C[3]*df_rec_lo[4]+df_rec_lo[3]*D_surf_C[4]); 
  surft1_up[0] = 0.3535533905932737*(D_surf_R[7]*df_rec_up[7]+D_surf_R[6]*df_rec_up[6]+D_surf_R[5]*df_rec_up[5]+D_surf_R[4]*df_rec_up[4]+D_surf_R[3]*df_rec_up[3]+D_surf_R[2]*df_rec_up[2]+D_surf_R[1]*df_rec_up[1]+D_surf_R[0]*df_rec_up[0]); 
  surft1_up[1] = 0.3535533905932737*(D_surf_R[6]*df_rec_up[7]+df_rec_up[6]*D_surf_R[7]+D_surf_R[3]*df_rec_up[5]+df_rec_up[3]*D_surf_R[5]+D_surf_R[2]*df_rec_up[4]+df_rec_up[2]*D_surf_R[4]+D_surf_R[0]*df_rec_up[1]+df_rec_up[0]*D_surf_R[1]); 
  surft1_up[2] = 0.3535533905932737*(D_surf_R[5]*df_rec_up[7]+df_rec_up[5]*D_surf_R[7]+D_surf_R[3]*df_rec_up[6]+df_rec_up[3]*D_surf_R[6]+D_surf_R[1]*df_rec_up[4]+df_rec_up[1]*D_surf_R[4]+D_surf_R[0]*df_rec_up[2]+df_rec_up[0]*D_surf_R[2]); 
  surft1_up[3] = 0.3535533905932737*(D_surf_R[4]*df_rec_up[7]+df_rec_up[4]*D_surf_R[7]+D_surf_R[2]*df_rec_up[6]+df_rec_up[2]*D_surf_R[6]+D_surf_R[1]*df_rec_up[5]+df_rec_up[1]*D_surf_R[5]+D_surf_R[0]*df_rec_up[3]+df_rec_up[0]*D_surf_R[3]); 
  surft1_up[4] = 0.3535533905932737*(D_surf_R[3]*df_rec_up[7]+df_rec_up[3]*D_surf_R[7]+D_surf_R[5]*df_rec_up[6]+df_rec_up[5]*D_surf_R[6]+D_surf_R[0]*df_rec_up[4]+df_rec_up[0]*D_surf_R[4]+D_surf_R[1]*df_rec_up[2]+df_rec_up[1]*D_surf_R[2]); 
  surft1_up[5] = 0.3535533905932737*(D_surf_R[2]*df_rec_up[7]+df_rec_up[2]*D_surf_R[7]+D_surf_R[4]*df_rec_up[6]+df_rec_up[4]*D_surf_R[6]+D_surf_R[0]*df_rec_up[5]+df_rec_up[0]*D_surf_R[5]+D_surf_R[1]*df_rec_up[3]+df_rec_up[1]*D_surf_R[3]); 
  surft1_up[6] = 0.3535533905932737*(D_surf_R[1]*df_rec_up[7]+df_rec_up[1]*D_surf_R[7]+D_surf_R[0]*df_rec_up[6]+df_rec_up[0]*D_surf_R[6]+D_surf_R[4]*df_rec_up[5]+df_rec_up[4]*D_surf_R[5]+D_surf_R[2]*df_rec_up[3]+df_rec_up[2]*D_surf_R[3]); 
  surft1_up[7] = 0.3535533905932737*(D_surf_R[0]*df_rec_up[7]+df_rec_up[0]*D_surf_R[7]+D_surf_R[1]*df_rec_up[6]+df_rec_up[1]*D_surf_R[6]+D_surf_R[2]*df_rec_up[5]+df_rec_up[2]*D_surf_R[5]+D_surf_R[3]*df_rec_up[4]+df_rec_up[3]*D_surf_R[4]); 

  surft2_lo[0] = 0.3535533905932737*(D_surf_C[7]*f_rec_lo[7]+D_surf_C[6]*f_rec_lo[6]+D_surf_C[5]*f_rec_lo[5]+D_surf_C[4]*f_rec_lo[4]+D_surf_C[3]*f_rec_lo[3]+D_surf_C[2]*f_rec_lo[2]+D_surf_C[1]*f_rec_lo[1]+D_surf_C[0]*f_rec_lo[0]); 
  surft2_lo[1] = 0.3535533905932737*(D_surf_C[6]*f_rec_lo[7]+f_rec_lo[6]*D_surf_C[7]+D_surf_C[3]*f_rec_lo[5]+f_rec_lo[3]*D_surf_C[5]+D_surf_C[2]*f_rec_lo[4]+f_rec_lo[2]*D_surf_C[4]+D_surf_C[0]*f_rec_lo[1]+f_rec_lo[0]*D_surf_C[1]); 
  surft2_lo[2] = 0.3535533905932737*(D_surf_C[5]*f_rec_lo[7]+f_rec_lo[5]*D_surf_C[7]+D_surf_C[3]*f_rec_lo[6]+f_rec_lo[3]*D_surf_C[6]+D_surf_C[1]*f_rec_lo[4]+f_rec_lo[1]*D_surf_C[4]+D_surf_C[0]*f_rec_lo[2]+f_rec_lo[0]*D_surf_C[2]); 
  surft2_lo[3] = 0.3535533905932737*(D_surf_C[4]*f_rec_lo[7]+f_rec_lo[4]*D_surf_C[7]+D_surf_C[2]*f_rec_lo[6]+f_rec_lo[2]*D_surf_C[6]+D_surf_C[1]*f_rec_lo[5]+f_rec_lo[1]*D_surf_C[5]+D_surf_C[0]*f_rec_lo[3]+f_rec_lo[0]*D_surf_C[3]); 
  surft2_lo[4] = 0.3535533905932737*(D_surf_C[3]*f_rec_lo[7]+f_rec_lo[3]*D_surf_C[7]+D_surf_C[5]*f_rec_lo[6]+f_rec_lo[5]*D_surf_C[6]+D_surf_C[0]*f_rec_lo[4]+f_rec_lo[0]*D_surf_C[4]+D_surf_C[1]*f_rec_lo[2]+f_rec_lo[1]*D_surf_C[2]); 
  surft2_lo[5] = 0.3535533905932737*(D_surf_C[2]*f_rec_lo[7]+f_rec_lo[2]*D_surf_C[7]+D_surf_C[4]*f_rec_lo[6]+f_rec_lo[4]*D_surf_C[6]+D_surf_C[0]*f_rec_lo[5]+f_rec_lo[0]*D_surf_C[5]+D_surf_C[1]*f_rec_lo[3]+f_rec_lo[1]*D_surf_C[3]); 
  surft2_lo[6] = 0.3535533905932737*(D_surf_C[1]*f_rec_lo[7]+f_rec_lo[1]*D_surf_C[7]+D_surf_C[0]*f_rec_lo[6]+f_rec_lo[0]*D_surf_C[6]+D_surf_C[4]*f_rec_lo[5]+f_rec_lo[4]*D_surf_C[5]+D_surf_C[2]*f_rec_lo[3]+f_rec_lo[2]*D_surf_C[3]); 
  surft2_lo[7] = 0.3535533905932737*(D_surf_C[0]*f_rec_lo[7]+f_rec_lo[0]*D_surf_C[7]+D_surf_C[1]*f_rec_lo[6]+f_rec_lo[1]*D_surf_C[6]+D_surf_C[2]*f_rec_lo[5]+f_rec_lo[2]*D_surf_C[5]+D_surf_C[3]*f_rec_lo[4]+f_rec_lo[3]*D_surf_C[4]); 
  surft2_up[0] = 0.3535533905932737*(D_surf_R[7]*f_rec_up[7]+D_surf_R[6]*f_rec_up[6]+D_surf_R[5]*f_rec_up[5]+D_surf_R[4]*f_rec_up[4]+D_surf_R[3]*f_rec_up[3]+D_surf_R[2]*f_rec_up[2]+D_surf_R[1]*f_rec_up[1]+D_surf_R[0]*f_rec_up[0]); 
  surft2_up[1] = 0.3535533905932737*(D_surf_R[6]*f_rec_up[7]+f_rec_up[6]*D_surf_R[7]+D_surf_R[3]*f_rec_up[5]+f_rec_up[3]*D_surf_R[5]+D_surf_R[2]*f_rec_up[4]+f_rec_up[2]*D_surf_R[4]+D_surf_R[0]*f_rec_up[1]+f_rec_up[0]*D_surf_R[1]); 
  surft2_up[2] = 0.3535533905932737*(D_surf_R[5]*f_rec_up[7]+f_rec_up[5]*D_surf_R[7]+D_surf_R[3]*f_rec_up[6]+f_rec_up[3]*D_surf_R[6]+D_surf_R[1]*f_rec_up[4]+f_rec_up[1]*D_surf_R[4]+D_surf_R[0]*f_rec_up[2]+f_rec_up[0]*D_surf_R[2]); 
  surft2_up[3] = 0.3535533905932737*(D_surf_R[4]*f_rec_up[7]+f_rec_up[4]*D_surf_R[7]+D_surf_R[2]*f_rec_up[6]+f_rec_up[2]*D_surf_R[6]+D_surf_R[1]*f_rec_up[5]+f_rec_up[1]*D_surf_R[5]+D_surf_R[0]*f_rec_up[3]+f_rec_up[0]*D_surf_R[3]); 
  surft2_up[4] = 0.3535533905932737*(D_surf_R[3]*f_rec_up[7]+f_rec_up[3]*D_surf_R[7]+D_surf_R[5]*f_rec_up[6]+f_rec_up[5]*D_surf_R[6]+D_surf_R[0]*f_rec_up[4]+f_rec_up[0]*D_surf_R[4]+D_surf_R[1]*f_rec_up[2]+f_rec_up[1]*D_surf_R[2]); 
  surft2_up[5] = 0.3535533905932737*(D_surf_R[2]*f_rec_up[7]+f_rec_up[2]*D_surf_R[7]+D_surf_R[4]*f_rec_up[6]+f_rec_up[4]*D_surf_R[6]+D_surf_R[0]*f_rec_up[5]+f_rec_up[0]*D_surf_R[5]+D_surf_R[1]*f_rec_up[3]+f_rec_up[1]*D_surf_R[3]); 
  surft2_up[6] = 0.3535533905932737*(D_surf_R[1]*f_rec_up[7]+f_rec_up[1]*D_surf_R[7]+D_surf_R[0]*f_rec_up[6]+f_rec_up[0]*D_surf_R[6]+D_surf_R[4]*f_rec_up[5]+f_rec_up[4]*D_surf_R[5]+D_surf_R[2]*f_rec_up[3]+f_rec_up[2]*D_surf_R[3]); 
  surft2_up[7] = 0.3535533905932737*(D_surf_R[0]*f_rec_up[7]+f_rec_up[0]*D_surf_R[7]+D_surf_R[1]*f_rec_up[6]+f_rec_up[1]*D_surf_R[6]+D_surf_R[2]*f_rec_up[5]+f_rec_up[2]*D_surf_R[5]+D_surf_R[3]*f_rec_up[4]+f_rec_up[3]*D_surf_R[4]); 

  vol[4] = 0.75*(fC[11]*DC[15]+fC[7]*DC[14]+fC[6]*DC[13]+fC[5]*DC[12]+fC[3]*DC[10]+fC[2]*DC[9]+fC[1]*DC[8]+fC[0]*DC[4]); 
  vol[8] = 0.75*(fC[7]*DC[15]+fC[11]*DC[14]+fC[3]*DC[13]+fC[2]*DC[12]+fC[6]*DC[10]+fC[5]*DC[9]+fC[0]*DC[8]+fC[1]*DC[4]); 
  vol[9] = 0.75*(fC[6]*DC[15]+fC[3]*DC[14]+fC[11]*DC[13]+fC[1]*DC[12]+fC[7]*DC[10]+fC[0]*DC[9]+fC[5]*DC[8]+fC[2]*DC[4]); 
  vol[10] = 0.75*(fC[5]*DC[15]+fC[2]*DC[14]+fC[1]*DC[13]+fC[11]*DC[12]+fC[0]*DC[10]+fC[7]*DC[9]+fC[6]*DC[8]+fC[3]*DC[4]); 
  vol[12] = 0.75*(fC[3]*DC[15]+fC[6]*DC[14]+fC[7]*DC[13]+fC[0]*DC[12]+DC[10]*fC[11]+fC[1]*DC[9]+fC[2]*DC[8]+DC[4]*fC[5]); 
  vol[13] = 0.75*(fC[2]*DC[15]+fC[5]*DC[14]+fC[0]*DC[13]+fC[7]*DC[12]+DC[9]*fC[11]+fC[1]*DC[10]+fC[3]*DC[8]+DC[4]*fC[6]); 
  vol[14] = 0.75*(fC[1]*DC[15]+fC[0]*DC[14]+fC[5]*DC[13]+fC[6]*DC[12]+DC[8]*fC[11]+fC[2]*DC[10]+fC[3]*DC[9]+DC[4]*fC[7]); 
  vol[15] = 0.75*(fC[0]*DC[15]+fC[1]*DC[14]+fC[2]*DC[13]+fC[3]*DC[12]+DC[4]*fC[11]+fC[5]*DC[10]+fC[6]*DC[9]+fC[7]*DC[8]); 

  out[0] += (0.5*vol[0]+0.35355339059327373*surft1_up[0]-0.35355339059327373*surft1_lo[0])*dv_inv_sq; 
  out[1] += (0.5*vol[1]+0.35355339059327373*surft1_up[1]-0.35355339059327373*surft1_lo[1])*dv_inv_sq; 
  out[2] += (0.5*vol[2]+0.35355339059327373*surft1_up[2]-0.35355339059327373*surft1_lo[2])*dv_inv_sq; 
  out[3] += (0.5*vol[3]+0.35355339059327373*surft1_up[3]-0.35355339059327373*surft1_lo[3])*dv_inv_sq; 
  out[4] += (0.5*vol[4]-0.6123724356957945*surft2_up[0]+0.6123724356957945*(surft2_lo[0]+surft1_up[0]+surft1_lo[0]))*dv_inv_sq; 
  out[5] += (0.5*vol[5]+0.35355339059327373*surft1_up[4]-0.35355339059327373*surft1_lo[4])*dv_inv_sq; 
  out[6] += (0.5*vol[6]+0.35355339059327373*surft1_up[5]-0.35355339059327373*surft1_lo[5])*dv_inv_sq; 
  out[7] += (0.5*vol[7]+0.35355339059327373*surft1_up[6]-0.35355339059327373*surft1_lo[6])*dv_inv_sq; 
  out[8] += (0.5*vol[8]-0.6123724356957945*surft2_up[1]+0.6123724356957945*(surft2_lo[1]+surft1_up[1]+surft1_lo[1]))*dv_inv_sq; 
  out[9] += (0.5*vol[9]-0.6123724356957945*surft2_up[2]+0.6123724356957945*(surft2_lo[2]+surft1_up[2]+surft1_lo[2]))*dv_inv_sq; 
  out[10] += (0.5*vol[10]-0.6123724356957945*surft2_up[3]+0.6123724356957945*(surft2_lo[3]+surft1_up[3]+surft1_lo[3]))*dv_inv_sq; 
  out[11] += (0.5*vol[11]+0.35355339059327373*surft1_up[7]-0.35355339059327373*surft1_lo[7])*dv_inv_sq; 
  out[12] += (0.5*vol[12]-0.6123724356957945*surft2_up[4]+0.6123724356957945*(surft2_lo[4]+surft1_up[4]+surft1_lo[4]))*dv_inv_sq; 
  out[13] += (0.5*vol[13]-0.6123724356957945*surft2_up[5]+0.6123724356957945*(surft2_lo[5]+surft1_up[5]+surft1_lo[5]))*dv_inv_sq; 
  out[14] += (0.5*vol[14]-0.6123724356957945*surft2_up[6]+0.6123724356957945*(surft2_lo[6]+surft1_up[6]+surft1_lo[6]))*dv_inv_sq; 
  out[15] += (0.5*vol[15]-0.6123724356957945*surft2_up[7]+0.6123724356957945*(surft2_lo[7]+surft1_up[7]+surft1_lo[7]))*dv_inv_sq; 

  double cflFreq = fmax(fabs(D_surf_C[0]), fabs(D_surf_R[0])); 

  return 1.4142135623730947*dv_inv_sq*cflFreq; 
} 
