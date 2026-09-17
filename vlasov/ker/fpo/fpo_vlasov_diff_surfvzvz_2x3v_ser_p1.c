#include <gkyl_fpo_vlasov_kernels.h> 

GKYL_CU_DH double fpo_vlasov_diff_surfvzvz_2x3v_ser_p1(const double* dxv, const double* diff_coeff_C, 
  const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* GKYL_RESTRICT out) { 
  // dxv[NDIM]: Cell spacing in each direction. 
  // diff_coeff_C: Diffusion tensor in center cell. 
  // diff_coeff_surf_stencil[9]: 9-cell stencil of surface expansion of recovered diffusion tensor. 
  // f_stencil[9]: 9-cell stencil of distribution function. 
  // out: Incremented output. 

  double dv_inv_sq = 4.0/dxv[4]/dxv[4]; 

  double f_rec_lo[16] = {0.0}; 
  double f_rec_up[16] = {0.0}; 
  double df_rec_lo[16] = {0.0}; 
  double df_rec_up[16] = {0.0}; 
  double surft1_lo[16] = {0.0}; 
  double surft1_up[16] = {0.0}; 
  double surft2_lo[16] = {0.0}; 
  double surft2_up[16] = {0.0}; 
  double vol[32] = {0.0}; 

  // Index into D and f stencils. 
  const double* fL = f_stencil[0]; 
  const double* fC = f_stencil[1]; 
  const double* fR = f_stencil[2]; 

  const double* DC = &diff_coeff_C[256]; 
  const double* D_surf_C = &diff_coeff_surf_stencil[1][256]; 
  const double* D_surf_R = &diff_coeff_surf_stencil[2][256]; 

  f_rec_lo[0] = 0.408248290463863*fL[5]-0.408248290463863*fC[5]+0.3535533905932737*(fL[0]+fC[0]); 
  f_rec_lo[1] = 0.408248290463863*fL[12]-0.408248290463863*fC[12]+0.3535533905932737*(fL[1]+fC[1]); 
  f_rec_lo[2] = 0.408248290463863*fL[13]-0.408248290463863*fC[13]+0.3535533905932737*(fL[2]+fC[2]); 
  f_rec_lo[3] = 0.408248290463863*fL[14]-0.408248290463863*fC[14]+0.3535533905932737*(fL[3]+fC[3]); 
  f_rec_lo[4] = 0.408248290463863*fL[15]-0.408248290463863*fC[15]+0.3535533905932737*(fL[4]+fC[4]); 
  f_rec_lo[5] = 0.408248290463863*fL[20]-0.408248290463863*fC[20]+0.3535533905932737*(fL[6]+fC[6]); 
  f_rec_lo[6] = 0.408248290463863*fL[21]-0.408248290463863*fC[21]+0.3535533905932737*(fL[7]+fC[7]); 
  f_rec_lo[7] = 0.408248290463863*fL[22]-0.408248290463863*fC[22]+0.3535533905932737*(fL[8]+fC[8]); 
  f_rec_lo[8] = 0.408248290463863*fL[23]-0.408248290463863*fC[23]+0.3535533905932737*(fL[9]+fC[9]); 
  f_rec_lo[9] = 0.408248290463863*fL[24]-0.408248290463863*fC[24]+0.3535533905932737*(fL[10]+fC[10]); 
  f_rec_lo[10] = 0.408248290463863*fL[25]-0.408248290463863*fC[25]+0.3535533905932737*(fL[11]+fC[11]); 
  f_rec_lo[11] = 0.408248290463863*fL[27]-0.408248290463863*fC[27]+0.3535533905932737*(fL[16]+fC[16]); 
  f_rec_lo[12] = 0.408248290463863*fL[28]-0.408248290463863*fC[28]+0.3535533905932737*(fL[17]+fC[17]); 
  f_rec_lo[13] = 0.408248290463863*fL[29]-0.408248290463863*fC[29]+0.3535533905932737*(fL[18]+fC[18]); 
  f_rec_lo[14] = 0.408248290463863*fL[30]-0.408248290463863*fC[30]+0.3535533905932737*(fL[19]+fC[19]); 
  f_rec_lo[15] = 0.408248290463863*fL[31]-0.408248290463863*fC[31]+0.3535533905932737*(fL[26]+fC[26]); 
  f_rec_up[0] = -(0.408248290463863*fR[5])+0.408248290463863*fC[5]+0.3535533905932737*(fR[0]+fC[0]); 
  f_rec_up[1] = -(0.408248290463863*fR[12])+0.408248290463863*fC[12]+0.3535533905932737*(fR[1]+fC[1]); 
  f_rec_up[2] = -(0.408248290463863*fR[13])+0.408248290463863*fC[13]+0.3535533905932737*(fR[2]+fC[2]); 
  f_rec_up[3] = -(0.408248290463863*fR[14])+0.408248290463863*fC[14]+0.3535533905932737*(fR[3]+fC[3]); 
  f_rec_up[4] = -(0.408248290463863*fR[15])+0.408248290463863*fC[15]+0.3535533905932737*(fR[4]+fC[4]); 
  f_rec_up[5] = -(0.408248290463863*fR[20])+0.408248290463863*fC[20]+0.3535533905932737*(fR[6]+fC[6]); 
  f_rec_up[6] = -(0.408248290463863*fR[21])+0.408248290463863*fC[21]+0.3535533905932737*(fR[7]+fC[7]); 
  f_rec_up[7] = -(0.408248290463863*fR[22])+0.408248290463863*fC[22]+0.3535533905932737*(fR[8]+fC[8]); 
  f_rec_up[8] = -(0.408248290463863*fR[23])+0.408248290463863*fC[23]+0.3535533905932737*(fR[9]+fC[9]); 
  f_rec_up[9] = -(0.408248290463863*fR[24])+0.408248290463863*fC[24]+0.3535533905932737*(fR[10]+fC[10]); 
  f_rec_up[10] = -(0.408248290463863*fR[25])+0.408248290463863*fC[25]+0.3535533905932737*(fR[11]+fC[11]); 
  f_rec_up[11] = -(0.408248290463863*fR[27])+0.408248290463863*fC[27]+0.3535533905932737*(fR[16]+fC[16]); 
  f_rec_up[12] = -(0.408248290463863*fR[28])+0.408248290463863*fC[28]+0.3535533905932737*(fR[17]+fC[17]); 
  f_rec_up[13] = -(0.408248290463863*fR[29])+0.408248290463863*fC[29]+0.3535533905932737*(fR[18]+fC[18]); 
  f_rec_up[14] = -(0.408248290463863*fR[30])+0.408248290463863*fC[30]+0.3535533905932737*(fR[19]+fC[19]); 
  f_rec_up[15] = -(0.408248290463863*fR[31])+0.408248290463863*fC[31]+0.3535533905932737*(fR[26]+fC[26]); 

  df_rec_lo[0] = -(0.7654655446197428*(fL[5]+fC[5]))-0.7954951288348656*fL[0]+0.7954951288348656*fC[0]; 
  df_rec_lo[1] = -(0.7654655446197428*(fL[12]+fC[12]))-0.7954951288348656*fL[1]+0.7954951288348656*fC[1]; 
  df_rec_lo[2] = -(0.7654655446197428*(fL[13]+fC[13]))-0.7954951288348656*fL[2]+0.7954951288348656*fC[2]; 
  df_rec_lo[3] = -(0.7654655446197428*(fL[14]+fC[14]))-0.7954951288348656*fL[3]+0.7954951288348656*fC[3]; 
  df_rec_lo[4] = -(0.7654655446197428*(fL[15]+fC[15]))-0.7954951288348656*fL[4]+0.7954951288348656*fC[4]; 
  df_rec_lo[5] = -(0.7654655446197428*(fL[20]+fC[20]))-0.7954951288348656*fL[6]+0.7954951288348656*fC[6]; 
  df_rec_lo[6] = -(0.7654655446197428*(fL[21]+fC[21]))-0.7954951288348656*fL[7]+0.7954951288348656*fC[7]; 
  df_rec_lo[7] = -(0.7654655446197428*(fL[22]+fC[22]))-0.7954951288348656*fL[8]+0.7954951288348656*fC[8]; 
  df_rec_lo[8] = -(0.7654655446197428*(fL[23]+fC[23]))-0.7954951288348656*fL[9]+0.7954951288348656*fC[9]; 
  df_rec_lo[9] = -(0.7654655446197428*(fL[24]+fC[24]))-0.7954951288348656*fL[10]+0.7954951288348656*fC[10]; 
  df_rec_lo[10] = -(0.7654655446197428*(fL[25]+fC[25]))-0.7954951288348656*fL[11]+0.7954951288348656*fC[11]; 
  df_rec_lo[11] = -(0.7654655446197428*(fL[27]+fC[27]))-0.7954951288348656*fL[16]+0.7954951288348656*fC[16]; 
  df_rec_lo[12] = -(0.7654655446197428*(fL[28]+fC[28]))-0.7954951288348656*fL[17]+0.7954951288348656*fC[17]; 
  df_rec_lo[13] = -(0.7654655446197428*(fL[29]+fC[29]))-0.7954951288348656*fL[18]+0.7954951288348656*fC[18]; 
  df_rec_lo[14] = -(0.7654655446197428*(fL[30]+fC[30]))-0.7954951288348656*fL[19]+0.7954951288348656*fC[19]; 
  df_rec_lo[15] = -(0.7654655446197428*(fL[31]+fC[31]))-0.7954951288348656*fL[26]+0.7954951288348656*fC[26]; 
  df_rec_up[0] = -(0.7654655446197428*(fR[5]+fC[5]))+0.7954951288348656*fR[0]-0.7954951288348656*fC[0]; 
  df_rec_up[1] = -(0.7654655446197428*(fR[12]+fC[12]))+0.7954951288348656*fR[1]-0.7954951288348656*fC[1]; 
  df_rec_up[2] = -(0.7654655446197428*(fR[13]+fC[13]))+0.7954951288348656*fR[2]-0.7954951288348656*fC[2]; 
  df_rec_up[3] = -(0.7654655446197428*(fR[14]+fC[14]))+0.7954951288348656*fR[3]-0.7954951288348656*fC[3]; 
  df_rec_up[4] = -(0.7654655446197428*(fR[15]+fC[15]))+0.7954951288348656*fR[4]-0.7954951288348656*fC[4]; 
  df_rec_up[5] = -(0.7654655446197428*(fR[20]+fC[20]))+0.7954951288348656*fR[6]-0.7954951288348656*fC[6]; 
  df_rec_up[6] = -(0.7654655446197428*(fR[21]+fC[21]))+0.7954951288348656*fR[7]-0.7954951288348656*fC[7]; 
  df_rec_up[7] = -(0.7654655446197428*(fR[22]+fC[22]))+0.7954951288348656*fR[8]-0.7954951288348656*fC[8]; 
  df_rec_up[8] = -(0.7654655446197428*(fR[23]+fC[23]))+0.7954951288348656*fR[9]-0.7954951288348656*fC[9]; 
  df_rec_up[9] = -(0.7654655446197428*(fR[24]+fC[24]))+0.7954951288348656*fR[10]-0.7954951288348656*fC[10]; 
  df_rec_up[10] = -(0.7654655446197428*(fR[25]+fC[25]))+0.7954951288348656*fR[11]-0.7954951288348656*fC[11]; 
  df_rec_up[11] = -(0.7654655446197428*(fR[27]+fC[27]))+0.7954951288348656*fR[16]-0.7954951288348656*fC[16]; 
  df_rec_up[12] = -(0.7654655446197428*(fR[28]+fC[28]))+0.7954951288348656*fR[17]-0.7954951288348656*fC[17]; 
  df_rec_up[13] = -(0.7654655446197428*(fR[29]+fC[29]))+0.7954951288348656*fR[18]-0.7954951288348656*fC[18]; 
  df_rec_up[14] = -(0.7654655446197428*(fR[30]+fC[30]))+0.7954951288348656*fR[19]-0.7954951288348656*fC[19]; 
  df_rec_up[15] = -(0.7654655446197428*(fR[31]+fC[31]))+0.7954951288348656*fR[26]-0.7954951288348656*fC[26]; 

  surft1_lo[0] = 0.25*(D_surf_C[15]*df_rec_lo[15]+D_surf_C[14]*df_rec_lo[14]+D_surf_C[13]*df_rec_lo[13]+D_surf_C[12]*df_rec_lo[12]+D_surf_C[11]*df_rec_lo[11]+D_surf_C[10]*df_rec_lo[10]+D_surf_C[9]*df_rec_lo[9]+D_surf_C[8]*df_rec_lo[8]+D_surf_C[7]*df_rec_lo[7]+D_surf_C[6]*df_rec_lo[6]+D_surf_C[5]*df_rec_lo[5]+D_surf_C[4]*df_rec_lo[4]+D_surf_C[3]*df_rec_lo[3]+D_surf_C[2]*df_rec_lo[2]+D_surf_C[1]*df_rec_lo[1]+D_surf_C[0]*df_rec_lo[0]); 
  surft1_lo[1] = 0.25*(D_surf_C[14]*df_rec_lo[15]+df_rec_lo[14]*D_surf_C[15]+D_surf_C[10]*df_rec_lo[13]+df_rec_lo[10]*D_surf_C[13]+D_surf_C[9]*df_rec_lo[12]+df_rec_lo[9]*D_surf_C[12]+D_surf_C[7]*df_rec_lo[11]+df_rec_lo[7]*D_surf_C[11]+D_surf_C[4]*df_rec_lo[8]+df_rec_lo[4]*D_surf_C[8]+D_surf_C[3]*df_rec_lo[6]+df_rec_lo[3]*D_surf_C[6]+D_surf_C[2]*df_rec_lo[5]+df_rec_lo[2]*D_surf_C[5]+D_surf_C[0]*df_rec_lo[1]+df_rec_lo[0]*D_surf_C[1]); 
  surft1_lo[2] = 0.25*(D_surf_C[13]*df_rec_lo[15]+df_rec_lo[13]*D_surf_C[15]+D_surf_C[10]*df_rec_lo[14]+df_rec_lo[10]*D_surf_C[14]+D_surf_C[8]*df_rec_lo[12]+df_rec_lo[8]*D_surf_C[12]+D_surf_C[6]*df_rec_lo[11]+df_rec_lo[6]*D_surf_C[11]+D_surf_C[4]*df_rec_lo[9]+df_rec_lo[4]*D_surf_C[9]+D_surf_C[3]*df_rec_lo[7]+df_rec_lo[3]*D_surf_C[7]+D_surf_C[1]*df_rec_lo[5]+df_rec_lo[1]*D_surf_C[5]+D_surf_C[0]*df_rec_lo[2]+df_rec_lo[0]*D_surf_C[2]); 
  surft1_lo[3] = 0.25*(D_surf_C[12]*df_rec_lo[15]+df_rec_lo[12]*D_surf_C[15]+D_surf_C[9]*df_rec_lo[14]+df_rec_lo[9]*D_surf_C[14]+D_surf_C[8]*df_rec_lo[13]+df_rec_lo[8]*D_surf_C[13]+D_surf_C[5]*df_rec_lo[11]+df_rec_lo[5]*D_surf_C[11]+D_surf_C[4]*df_rec_lo[10]+df_rec_lo[4]*D_surf_C[10]+D_surf_C[2]*df_rec_lo[7]+df_rec_lo[2]*D_surf_C[7]+D_surf_C[1]*df_rec_lo[6]+df_rec_lo[1]*D_surf_C[6]+D_surf_C[0]*df_rec_lo[3]+df_rec_lo[0]*D_surf_C[3]); 
  surft1_lo[4] = 0.25*(D_surf_C[11]*df_rec_lo[15]+df_rec_lo[11]*D_surf_C[15]+D_surf_C[7]*df_rec_lo[14]+df_rec_lo[7]*D_surf_C[14]+D_surf_C[6]*df_rec_lo[13]+df_rec_lo[6]*D_surf_C[13]+D_surf_C[5]*df_rec_lo[12]+df_rec_lo[5]*D_surf_C[12]+D_surf_C[3]*df_rec_lo[10]+df_rec_lo[3]*D_surf_C[10]+D_surf_C[2]*df_rec_lo[9]+df_rec_lo[2]*D_surf_C[9]+D_surf_C[1]*df_rec_lo[8]+df_rec_lo[1]*D_surf_C[8]+D_surf_C[0]*df_rec_lo[4]+df_rec_lo[0]*D_surf_C[4]); 
  surft1_lo[5] = 0.25*(D_surf_C[10]*df_rec_lo[15]+df_rec_lo[10]*D_surf_C[15]+D_surf_C[13]*df_rec_lo[14]+df_rec_lo[13]*D_surf_C[14]+D_surf_C[4]*df_rec_lo[12]+df_rec_lo[4]*D_surf_C[12]+D_surf_C[3]*df_rec_lo[11]+df_rec_lo[3]*D_surf_C[11]+D_surf_C[8]*df_rec_lo[9]+df_rec_lo[8]*D_surf_C[9]+D_surf_C[6]*df_rec_lo[7]+df_rec_lo[6]*D_surf_C[7]+D_surf_C[0]*df_rec_lo[5]+df_rec_lo[0]*D_surf_C[5]+D_surf_C[1]*df_rec_lo[2]+df_rec_lo[1]*D_surf_C[2]); 
  surft1_lo[6] = 0.25*(D_surf_C[9]*df_rec_lo[15]+df_rec_lo[9]*D_surf_C[15]+D_surf_C[12]*df_rec_lo[14]+df_rec_lo[12]*D_surf_C[14]+D_surf_C[4]*df_rec_lo[13]+df_rec_lo[4]*D_surf_C[13]+D_surf_C[2]*df_rec_lo[11]+df_rec_lo[2]*D_surf_C[11]+D_surf_C[8]*df_rec_lo[10]+df_rec_lo[8]*D_surf_C[10]+D_surf_C[5]*df_rec_lo[7]+df_rec_lo[5]*D_surf_C[7]+D_surf_C[0]*df_rec_lo[6]+df_rec_lo[0]*D_surf_C[6]+D_surf_C[1]*df_rec_lo[3]+df_rec_lo[1]*D_surf_C[3]); 
  surft1_lo[7] = 0.25*(D_surf_C[8]*df_rec_lo[15]+df_rec_lo[8]*D_surf_C[15]+D_surf_C[4]*df_rec_lo[14]+df_rec_lo[4]*D_surf_C[14]+D_surf_C[12]*df_rec_lo[13]+df_rec_lo[12]*D_surf_C[13]+D_surf_C[1]*df_rec_lo[11]+df_rec_lo[1]*D_surf_C[11]+D_surf_C[9]*df_rec_lo[10]+df_rec_lo[9]*D_surf_C[10]+D_surf_C[0]*df_rec_lo[7]+df_rec_lo[0]*D_surf_C[7]+D_surf_C[5]*df_rec_lo[6]+df_rec_lo[5]*D_surf_C[6]+D_surf_C[2]*df_rec_lo[3]+df_rec_lo[2]*D_surf_C[3]); 
  surft1_lo[8] = 0.25*(D_surf_C[7]*df_rec_lo[15]+df_rec_lo[7]*D_surf_C[15]+D_surf_C[11]*df_rec_lo[14]+df_rec_lo[11]*D_surf_C[14]+D_surf_C[3]*df_rec_lo[13]+df_rec_lo[3]*D_surf_C[13]+D_surf_C[2]*df_rec_lo[12]+df_rec_lo[2]*D_surf_C[12]+D_surf_C[6]*df_rec_lo[10]+df_rec_lo[6]*D_surf_C[10]+D_surf_C[5]*df_rec_lo[9]+df_rec_lo[5]*D_surf_C[9]+D_surf_C[0]*df_rec_lo[8]+df_rec_lo[0]*D_surf_C[8]+D_surf_C[1]*df_rec_lo[4]+df_rec_lo[1]*D_surf_C[4]); 
  surft1_lo[9] = 0.25*(D_surf_C[6]*df_rec_lo[15]+df_rec_lo[6]*D_surf_C[15]+D_surf_C[3]*df_rec_lo[14]+df_rec_lo[3]*D_surf_C[14]+D_surf_C[11]*df_rec_lo[13]+df_rec_lo[11]*D_surf_C[13]+D_surf_C[1]*df_rec_lo[12]+df_rec_lo[1]*D_surf_C[12]+D_surf_C[7]*df_rec_lo[10]+df_rec_lo[7]*D_surf_C[10]+D_surf_C[0]*df_rec_lo[9]+df_rec_lo[0]*D_surf_C[9]+D_surf_C[5]*df_rec_lo[8]+df_rec_lo[5]*D_surf_C[8]+D_surf_C[2]*df_rec_lo[4]+df_rec_lo[2]*D_surf_C[4]); 
  surft1_lo[10] = 0.25*(D_surf_C[5]*df_rec_lo[15]+df_rec_lo[5]*D_surf_C[15]+D_surf_C[2]*df_rec_lo[14]+df_rec_lo[2]*D_surf_C[14]+D_surf_C[1]*df_rec_lo[13]+df_rec_lo[1]*D_surf_C[13]+D_surf_C[11]*df_rec_lo[12]+df_rec_lo[11]*D_surf_C[12]+D_surf_C[0]*df_rec_lo[10]+df_rec_lo[0]*D_surf_C[10]+D_surf_C[7]*df_rec_lo[9]+df_rec_lo[7]*D_surf_C[9]+D_surf_C[6]*df_rec_lo[8]+df_rec_lo[6]*D_surf_C[8]+D_surf_C[3]*df_rec_lo[4]+df_rec_lo[3]*D_surf_C[4]); 
  surft1_lo[11] = 0.25*(D_surf_C[4]*df_rec_lo[15]+df_rec_lo[4]*D_surf_C[15]+D_surf_C[8]*df_rec_lo[14]+df_rec_lo[8]*D_surf_C[14]+D_surf_C[9]*df_rec_lo[13]+df_rec_lo[9]*D_surf_C[13]+D_surf_C[10]*df_rec_lo[12]+df_rec_lo[10]*D_surf_C[12]+D_surf_C[0]*df_rec_lo[11]+df_rec_lo[0]*D_surf_C[11]+D_surf_C[1]*df_rec_lo[7]+df_rec_lo[1]*D_surf_C[7]+D_surf_C[2]*df_rec_lo[6]+df_rec_lo[2]*D_surf_C[6]+D_surf_C[3]*df_rec_lo[5]+df_rec_lo[3]*D_surf_C[5]); 
  surft1_lo[12] = 0.25*(D_surf_C[3]*df_rec_lo[15]+df_rec_lo[3]*D_surf_C[15]+D_surf_C[6]*df_rec_lo[14]+df_rec_lo[6]*D_surf_C[14]+D_surf_C[7]*df_rec_lo[13]+df_rec_lo[7]*D_surf_C[13]+D_surf_C[0]*df_rec_lo[12]+df_rec_lo[0]*D_surf_C[12]+D_surf_C[10]*df_rec_lo[11]+df_rec_lo[10]*D_surf_C[11]+D_surf_C[1]*df_rec_lo[9]+df_rec_lo[1]*D_surf_C[9]+D_surf_C[2]*df_rec_lo[8]+df_rec_lo[2]*D_surf_C[8]+D_surf_C[4]*df_rec_lo[5]+df_rec_lo[4]*D_surf_C[5]); 
  surft1_lo[13] = 0.25*(D_surf_C[2]*df_rec_lo[15]+df_rec_lo[2]*D_surf_C[15]+D_surf_C[5]*df_rec_lo[14]+df_rec_lo[5]*D_surf_C[14]+D_surf_C[0]*df_rec_lo[13]+df_rec_lo[0]*D_surf_C[13]+D_surf_C[7]*df_rec_lo[12]+df_rec_lo[7]*D_surf_C[12]+D_surf_C[9]*df_rec_lo[11]+df_rec_lo[9]*D_surf_C[11]+D_surf_C[1]*df_rec_lo[10]+df_rec_lo[1]*D_surf_C[10]+D_surf_C[3]*df_rec_lo[8]+df_rec_lo[3]*D_surf_C[8]+D_surf_C[4]*df_rec_lo[6]+df_rec_lo[4]*D_surf_C[6]); 
  surft1_lo[14] = 0.25*(D_surf_C[1]*df_rec_lo[15]+df_rec_lo[1]*D_surf_C[15]+D_surf_C[0]*df_rec_lo[14]+df_rec_lo[0]*D_surf_C[14]+D_surf_C[5]*df_rec_lo[13]+df_rec_lo[5]*D_surf_C[13]+D_surf_C[6]*df_rec_lo[12]+df_rec_lo[6]*D_surf_C[12]+D_surf_C[8]*df_rec_lo[11]+df_rec_lo[8]*D_surf_C[11]+D_surf_C[2]*df_rec_lo[10]+df_rec_lo[2]*D_surf_C[10]+D_surf_C[3]*df_rec_lo[9]+df_rec_lo[3]*D_surf_C[9]+D_surf_C[4]*df_rec_lo[7]+df_rec_lo[4]*D_surf_C[7]); 
  surft1_lo[15] = 0.25*(D_surf_C[0]*df_rec_lo[15]+df_rec_lo[0]*D_surf_C[15]+D_surf_C[1]*df_rec_lo[14]+df_rec_lo[1]*D_surf_C[14]+D_surf_C[2]*df_rec_lo[13]+df_rec_lo[2]*D_surf_C[13]+D_surf_C[3]*df_rec_lo[12]+df_rec_lo[3]*D_surf_C[12]+D_surf_C[4]*df_rec_lo[11]+df_rec_lo[4]*D_surf_C[11]+D_surf_C[5]*df_rec_lo[10]+df_rec_lo[5]*D_surf_C[10]+D_surf_C[6]*df_rec_lo[9]+df_rec_lo[6]*D_surf_C[9]+D_surf_C[7]*df_rec_lo[8]+df_rec_lo[7]*D_surf_C[8]); 
  surft1_up[0] = 0.25*(D_surf_R[15]*df_rec_up[15]+D_surf_R[14]*df_rec_up[14]+D_surf_R[13]*df_rec_up[13]+D_surf_R[12]*df_rec_up[12]+D_surf_R[11]*df_rec_up[11]+D_surf_R[10]*df_rec_up[10]+D_surf_R[9]*df_rec_up[9]+D_surf_R[8]*df_rec_up[8]+D_surf_R[7]*df_rec_up[7]+D_surf_R[6]*df_rec_up[6]+D_surf_R[5]*df_rec_up[5]+D_surf_R[4]*df_rec_up[4]+D_surf_R[3]*df_rec_up[3]+D_surf_R[2]*df_rec_up[2]+D_surf_R[1]*df_rec_up[1]+D_surf_R[0]*df_rec_up[0]); 
  surft1_up[1] = 0.25*(D_surf_R[14]*df_rec_up[15]+df_rec_up[14]*D_surf_R[15]+D_surf_R[10]*df_rec_up[13]+df_rec_up[10]*D_surf_R[13]+D_surf_R[9]*df_rec_up[12]+df_rec_up[9]*D_surf_R[12]+D_surf_R[7]*df_rec_up[11]+df_rec_up[7]*D_surf_R[11]+D_surf_R[4]*df_rec_up[8]+df_rec_up[4]*D_surf_R[8]+D_surf_R[3]*df_rec_up[6]+df_rec_up[3]*D_surf_R[6]+D_surf_R[2]*df_rec_up[5]+df_rec_up[2]*D_surf_R[5]+D_surf_R[0]*df_rec_up[1]+df_rec_up[0]*D_surf_R[1]); 
  surft1_up[2] = 0.25*(D_surf_R[13]*df_rec_up[15]+df_rec_up[13]*D_surf_R[15]+D_surf_R[10]*df_rec_up[14]+df_rec_up[10]*D_surf_R[14]+D_surf_R[8]*df_rec_up[12]+df_rec_up[8]*D_surf_R[12]+D_surf_R[6]*df_rec_up[11]+df_rec_up[6]*D_surf_R[11]+D_surf_R[4]*df_rec_up[9]+df_rec_up[4]*D_surf_R[9]+D_surf_R[3]*df_rec_up[7]+df_rec_up[3]*D_surf_R[7]+D_surf_R[1]*df_rec_up[5]+df_rec_up[1]*D_surf_R[5]+D_surf_R[0]*df_rec_up[2]+df_rec_up[0]*D_surf_R[2]); 
  surft1_up[3] = 0.25*(D_surf_R[12]*df_rec_up[15]+df_rec_up[12]*D_surf_R[15]+D_surf_R[9]*df_rec_up[14]+df_rec_up[9]*D_surf_R[14]+D_surf_R[8]*df_rec_up[13]+df_rec_up[8]*D_surf_R[13]+D_surf_R[5]*df_rec_up[11]+df_rec_up[5]*D_surf_R[11]+D_surf_R[4]*df_rec_up[10]+df_rec_up[4]*D_surf_R[10]+D_surf_R[2]*df_rec_up[7]+df_rec_up[2]*D_surf_R[7]+D_surf_R[1]*df_rec_up[6]+df_rec_up[1]*D_surf_R[6]+D_surf_R[0]*df_rec_up[3]+df_rec_up[0]*D_surf_R[3]); 
  surft1_up[4] = 0.25*(D_surf_R[11]*df_rec_up[15]+df_rec_up[11]*D_surf_R[15]+D_surf_R[7]*df_rec_up[14]+df_rec_up[7]*D_surf_R[14]+D_surf_R[6]*df_rec_up[13]+df_rec_up[6]*D_surf_R[13]+D_surf_R[5]*df_rec_up[12]+df_rec_up[5]*D_surf_R[12]+D_surf_R[3]*df_rec_up[10]+df_rec_up[3]*D_surf_R[10]+D_surf_R[2]*df_rec_up[9]+df_rec_up[2]*D_surf_R[9]+D_surf_R[1]*df_rec_up[8]+df_rec_up[1]*D_surf_R[8]+D_surf_R[0]*df_rec_up[4]+df_rec_up[0]*D_surf_R[4]); 
  surft1_up[5] = 0.25*(D_surf_R[10]*df_rec_up[15]+df_rec_up[10]*D_surf_R[15]+D_surf_R[13]*df_rec_up[14]+df_rec_up[13]*D_surf_R[14]+D_surf_R[4]*df_rec_up[12]+df_rec_up[4]*D_surf_R[12]+D_surf_R[3]*df_rec_up[11]+df_rec_up[3]*D_surf_R[11]+D_surf_R[8]*df_rec_up[9]+df_rec_up[8]*D_surf_R[9]+D_surf_R[6]*df_rec_up[7]+df_rec_up[6]*D_surf_R[7]+D_surf_R[0]*df_rec_up[5]+df_rec_up[0]*D_surf_R[5]+D_surf_R[1]*df_rec_up[2]+df_rec_up[1]*D_surf_R[2]); 
  surft1_up[6] = 0.25*(D_surf_R[9]*df_rec_up[15]+df_rec_up[9]*D_surf_R[15]+D_surf_R[12]*df_rec_up[14]+df_rec_up[12]*D_surf_R[14]+D_surf_R[4]*df_rec_up[13]+df_rec_up[4]*D_surf_R[13]+D_surf_R[2]*df_rec_up[11]+df_rec_up[2]*D_surf_R[11]+D_surf_R[8]*df_rec_up[10]+df_rec_up[8]*D_surf_R[10]+D_surf_R[5]*df_rec_up[7]+df_rec_up[5]*D_surf_R[7]+D_surf_R[0]*df_rec_up[6]+df_rec_up[0]*D_surf_R[6]+D_surf_R[1]*df_rec_up[3]+df_rec_up[1]*D_surf_R[3]); 
  surft1_up[7] = 0.25*(D_surf_R[8]*df_rec_up[15]+df_rec_up[8]*D_surf_R[15]+D_surf_R[4]*df_rec_up[14]+df_rec_up[4]*D_surf_R[14]+D_surf_R[12]*df_rec_up[13]+df_rec_up[12]*D_surf_R[13]+D_surf_R[1]*df_rec_up[11]+df_rec_up[1]*D_surf_R[11]+D_surf_R[9]*df_rec_up[10]+df_rec_up[9]*D_surf_R[10]+D_surf_R[0]*df_rec_up[7]+df_rec_up[0]*D_surf_R[7]+D_surf_R[5]*df_rec_up[6]+df_rec_up[5]*D_surf_R[6]+D_surf_R[2]*df_rec_up[3]+df_rec_up[2]*D_surf_R[3]); 
  surft1_up[8] = 0.25*(D_surf_R[7]*df_rec_up[15]+df_rec_up[7]*D_surf_R[15]+D_surf_R[11]*df_rec_up[14]+df_rec_up[11]*D_surf_R[14]+D_surf_R[3]*df_rec_up[13]+df_rec_up[3]*D_surf_R[13]+D_surf_R[2]*df_rec_up[12]+df_rec_up[2]*D_surf_R[12]+D_surf_R[6]*df_rec_up[10]+df_rec_up[6]*D_surf_R[10]+D_surf_R[5]*df_rec_up[9]+df_rec_up[5]*D_surf_R[9]+D_surf_R[0]*df_rec_up[8]+df_rec_up[0]*D_surf_R[8]+D_surf_R[1]*df_rec_up[4]+df_rec_up[1]*D_surf_R[4]); 
  surft1_up[9] = 0.25*(D_surf_R[6]*df_rec_up[15]+df_rec_up[6]*D_surf_R[15]+D_surf_R[3]*df_rec_up[14]+df_rec_up[3]*D_surf_R[14]+D_surf_R[11]*df_rec_up[13]+df_rec_up[11]*D_surf_R[13]+D_surf_R[1]*df_rec_up[12]+df_rec_up[1]*D_surf_R[12]+D_surf_R[7]*df_rec_up[10]+df_rec_up[7]*D_surf_R[10]+D_surf_R[0]*df_rec_up[9]+df_rec_up[0]*D_surf_R[9]+D_surf_R[5]*df_rec_up[8]+df_rec_up[5]*D_surf_R[8]+D_surf_R[2]*df_rec_up[4]+df_rec_up[2]*D_surf_R[4]); 
  surft1_up[10] = 0.25*(D_surf_R[5]*df_rec_up[15]+df_rec_up[5]*D_surf_R[15]+D_surf_R[2]*df_rec_up[14]+df_rec_up[2]*D_surf_R[14]+D_surf_R[1]*df_rec_up[13]+df_rec_up[1]*D_surf_R[13]+D_surf_R[11]*df_rec_up[12]+df_rec_up[11]*D_surf_R[12]+D_surf_R[0]*df_rec_up[10]+df_rec_up[0]*D_surf_R[10]+D_surf_R[7]*df_rec_up[9]+df_rec_up[7]*D_surf_R[9]+D_surf_R[6]*df_rec_up[8]+df_rec_up[6]*D_surf_R[8]+D_surf_R[3]*df_rec_up[4]+df_rec_up[3]*D_surf_R[4]); 
  surft1_up[11] = 0.25*(D_surf_R[4]*df_rec_up[15]+df_rec_up[4]*D_surf_R[15]+D_surf_R[8]*df_rec_up[14]+df_rec_up[8]*D_surf_R[14]+D_surf_R[9]*df_rec_up[13]+df_rec_up[9]*D_surf_R[13]+D_surf_R[10]*df_rec_up[12]+df_rec_up[10]*D_surf_R[12]+D_surf_R[0]*df_rec_up[11]+df_rec_up[0]*D_surf_R[11]+D_surf_R[1]*df_rec_up[7]+df_rec_up[1]*D_surf_R[7]+D_surf_R[2]*df_rec_up[6]+df_rec_up[2]*D_surf_R[6]+D_surf_R[3]*df_rec_up[5]+df_rec_up[3]*D_surf_R[5]); 
  surft1_up[12] = 0.25*(D_surf_R[3]*df_rec_up[15]+df_rec_up[3]*D_surf_R[15]+D_surf_R[6]*df_rec_up[14]+df_rec_up[6]*D_surf_R[14]+D_surf_R[7]*df_rec_up[13]+df_rec_up[7]*D_surf_R[13]+D_surf_R[0]*df_rec_up[12]+df_rec_up[0]*D_surf_R[12]+D_surf_R[10]*df_rec_up[11]+df_rec_up[10]*D_surf_R[11]+D_surf_R[1]*df_rec_up[9]+df_rec_up[1]*D_surf_R[9]+D_surf_R[2]*df_rec_up[8]+df_rec_up[2]*D_surf_R[8]+D_surf_R[4]*df_rec_up[5]+df_rec_up[4]*D_surf_R[5]); 
  surft1_up[13] = 0.25*(D_surf_R[2]*df_rec_up[15]+df_rec_up[2]*D_surf_R[15]+D_surf_R[5]*df_rec_up[14]+df_rec_up[5]*D_surf_R[14]+D_surf_R[0]*df_rec_up[13]+df_rec_up[0]*D_surf_R[13]+D_surf_R[7]*df_rec_up[12]+df_rec_up[7]*D_surf_R[12]+D_surf_R[9]*df_rec_up[11]+df_rec_up[9]*D_surf_R[11]+D_surf_R[1]*df_rec_up[10]+df_rec_up[1]*D_surf_R[10]+D_surf_R[3]*df_rec_up[8]+df_rec_up[3]*D_surf_R[8]+D_surf_R[4]*df_rec_up[6]+df_rec_up[4]*D_surf_R[6]); 
  surft1_up[14] = 0.25*(D_surf_R[1]*df_rec_up[15]+df_rec_up[1]*D_surf_R[15]+D_surf_R[0]*df_rec_up[14]+df_rec_up[0]*D_surf_R[14]+D_surf_R[5]*df_rec_up[13]+df_rec_up[5]*D_surf_R[13]+D_surf_R[6]*df_rec_up[12]+df_rec_up[6]*D_surf_R[12]+D_surf_R[8]*df_rec_up[11]+df_rec_up[8]*D_surf_R[11]+D_surf_R[2]*df_rec_up[10]+df_rec_up[2]*D_surf_R[10]+D_surf_R[3]*df_rec_up[9]+df_rec_up[3]*D_surf_R[9]+D_surf_R[4]*df_rec_up[7]+df_rec_up[4]*D_surf_R[7]); 
  surft1_up[15] = 0.25*(D_surf_R[0]*df_rec_up[15]+df_rec_up[0]*D_surf_R[15]+D_surf_R[1]*df_rec_up[14]+df_rec_up[1]*D_surf_R[14]+D_surf_R[2]*df_rec_up[13]+df_rec_up[2]*D_surf_R[13]+D_surf_R[3]*df_rec_up[12]+df_rec_up[3]*D_surf_R[12]+D_surf_R[4]*df_rec_up[11]+df_rec_up[4]*D_surf_R[11]+D_surf_R[5]*df_rec_up[10]+df_rec_up[5]*D_surf_R[10]+D_surf_R[6]*df_rec_up[9]+df_rec_up[6]*D_surf_R[9]+D_surf_R[7]*df_rec_up[8]+df_rec_up[7]*D_surf_R[8]); 

  surft2_lo[0] = 0.25*(D_surf_C[15]*f_rec_lo[15]+D_surf_C[14]*f_rec_lo[14]+D_surf_C[13]*f_rec_lo[13]+D_surf_C[12]*f_rec_lo[12]+D_surf_C[11]*f_rec_lo[11]+D_surf_C[10]*f_rec_lo[10]+D_surf_C[9]*f_rec_lo[9]+D_surf_C[8]*f_rec_lo[8]+D_surf_C[7]*f_rec_lo[7]+D_surf_C[6]*f_rec_lo[6]+D_surf_C[5]*f_rec_lo[5]+D_surf_C[4]*f_rec_lo[4]+D_surf_C[3]*f_rec_lo[3]+D_surf_C[2]*f_rec_lo[2]+D_surf_C[1]*f_rec_lo[1]+D_surf_C[0]*f_rec_lo[0]); 
  surft2_lo[1] = 0.25*(D_surf_C[14]*f_rec_lo[15]+f_rec_lo[14]*D_surf_C[15]+D_surf_C[10]*f_rec_lo[13]+f_rec_lo[10]*D_surf_C[13]+D_surf_C[9]*f_rec_lo[12]+f_rec_lo[9]*D_surf_C[12]+D_surf_C[7]*f_rec_lo[11]+f_rec_lo[7]*D_surf_C[11]+D_surf_C[4]*f_rec_lo[8]+f_rec_lo[4]*D_surf_C[8]+D_surf_C[3]*f_rec_lo[6]+f_rec_lo[3]*D_surf_C[6]+D_surf_C[2]*f_rec_lo[5]+f_rec_lo[2]*D_surf_C[5]+D_surf_C[0]*f_rec_lo[1]+f_rec_lo[0]*D_surf_C[1]); 
  surft2_lo[2] = 0.25*(D_surf_C[13]*f_rec_lo[15]+f_rec_lo[13]*D_surf_C[15]+D_surf_C[10]*f_rec_lo[14]+f_rec_lo[10]*D_surf_C[14]+D_surf_C[8]*f_rec_lo[12]+f_rec_lo[8]*D_surf_C[12]+D_surf_C[6]*f_rec_lo[11]+f_rec_lo[6]*D_surf_C[11]+D_surf_C[4]*f_rec_lo[9]+f_rec_lo[4]*D_surf_C[9]+D_surf_C[3]*f_rec_lo[7]+f_rec_lo[3]*D_surf_C[7]+D_surf_C[1]*f_rec_lo[5]+f_rec_lo[1]*D_surf_C[5]+D_surf_C[0]*f_rec_lo[2]+f_rec_lo[0]*D_surf_C[2]); 
  surft2_lo[3] = 0.25*(D_surf_C[12]*f_rec_lo[15]+f_rec_lo[12]*D_surf_C[15]+D_surf_C[9]*f_rec_lo[14]+f_rec_lo[9]*D_surf_C[14]+D_surf_C[8]*f_rec_lo[13]+f_rec_lo[8]*D_surf_C[13]+D_surf_C[5]*f_rec_lo[11]+f_rec_lo[5]*D_surf_C[11]+D_surf_C[4]*f_rec_lo[10]+f_rec_lo[4]*D_surf_C[10]+D_surf_C[2]*f_rec_lo[7]+f_rec_lo[2]*D_surf_C[7]+D_surf_C[1]*f_rec_lo[6]+f_rec_lo[1]*D_surf_C[6]+D_surf_C[0]*f_rec_lo[3]+f_rec_lo[0]*D_surf_C[3]); 
  surft2_lo[4] = 0.25*(D_surf_C[11]*f_rec_lo[15]+f_rec_lo[11]*D_surf_C[15]+D_surf_C[7]*f_rec_lo[14]+f_rec_lo[7]*D_surf_C[14]+D_surf_C[6]*f_rec_lo[13]+f_rec_lo[6]*D_surf_C[13]+D_surf_C[5]*f_rec_lo[12]+f_rec_lo[5]*D_surf_C[12]+D_surf_C[3]*f_rec_lo[10]+f_rec_lo[3]*D_surf_C[10]+D_surf_C[2]*f_rec_lo[9]+f_rec_lo[2]*D_surf_C[9]+D_surf_C[1]*f_rec_lo[8]+f_rec_lo[1]*D_surf_C[8]+D_surf_C[0]*f_rec_lo[4]+f_rec_lo[0]*D_surf_C[4]); 
  surft2_lo[5] = 0.25*(D_surf_C[10]*f_rec_lo[15]+f_rec_lo[10]*D_surf_C[15]+D_surf_C[13]*f_rec_lo[14]+f_rec_lo[13]*D_surf_C[14]+D_surf_C[4]*f_rec_lo[12]+f_rec_lo[4]*D_surf_C[12]+D_surf_C[3]*f_rec_lo[11]+f_rec_lo[3]*D_surf_C[11]+D_surf_C[8]*f_rec_lo[9]+f_rec_lo[8]*D_surf_C[9]+D_surf_C[6]*f_rec_lo[7]+f_rec_lo[6]*D_surf_C[7]+D_surf_C[0]*f_rec_lo[5]+f_rec_lo[0]*D_surf_C[5]+D_surf_C[1]*f_rec_lo[2]+f_rec_lo[1]*D_surf_C[2]); 
  surft2_lo[6] = 0.25*(D_surf_C[9]*f_rec_lo[15]+f_rec_lo[9]*D_surf_C[15]+D_surf_C[12]*f_rec_lo[14]+f_rec_lo[12]*D_surf_C[14]+D_surf_C[4]*f_rec_lo[13]+f_rec_lo[4]*D_surf_C[13]+D_surf_C[2]*f_rec_lo[11]+f_rec_lo[2]*D_surf_C[11]+D_surf_C[8]*f_rec_lo[10]+f_rec_lo[8]*D_surf_C[10]+D_surf_C[5]*f_rec_lo[7]+f_rec_lo[5]*D_surf_C[7]+D_surf_C[0]*f_rec_lo[6]+f_rec_lo[0]*D_surf_C[6]+D_surf_C[1]*f_rec_lo[3]+f_rec_lo[1]*D_surf_C[3]); 
  surft2_lo[7] = 0.25*(D_surf_C[8]*f_rec_lo[15]+f_rec_lo[8]*D_surf_C[15]+D_surf_C[4]*f_rec_lo[14]+f_rec_lo[4]*D_surf_C[14]+D_surf_C[12]*f_rec_lo[13]+f_rec_lo[12]*D_surf_C[13]+D_surf_C[1]*f_rec_lo[11]+f_rec_lo[1]*D_surf_C[11]+D_surf_C[9]*f_rec_lo[10]+f_rec_lo[9]*D_surf_C[10]+D_surf_C[0]*f_rec_lo[7]+f_rec_lo[0]*D_surf_C[7]+D_surf_C[5]*f_rec_lo[6]+f_rec_lo[5]*D_surf_C[6]+D_surf_C[2]*f_rec_lo[3]+f_rec_lo[2]*D_surf_C[3]); 
  surft2_lo[8] = 0.25*(D_surf_C[7]*f_rec_lo[15]+f_rec_lo[7]*D_surf_C[15]+D_surf_C[11]*f_rec_lo[14]+f_rec_lo[11]*D_surf_C[14]+D_surf_C[3]*f_rec_lo[13]+f_rec_lo[3]*D_surf_C[13]+D_surf_C[2]*f_rec_lo[12]+f_rec_lo[2]*D_surf_C[12]+D_surf_C[6]*f_rec_lo[10]+f_rec_lo[6]*D_surf_C[10]+D_surf_C[5]*f_rec_lo[9]+f_rec_lo[5]*D_surf_C[9]+D_surf_C[0]*f_rec_lo[8]+f_rec_lo[0]*D_surf_C[8]+D_surf_C[1]*f_rec_lo[4]+f_rec_lo[1]*D_surf_C[4]); 
  surft2_lo[9] = 0.25*(D_surf_C[6]*f_rec_lo[15]+f_rec_lo[6]*D_surf_C[15]+D_surf_C[3]*f_rec_lo[14]+f_rec_lo[3]*D_surf_C[14]+D_surf_C[11]*f_rec_lo[13]+f_rec_lo[11]*D_surf_C[13]+D_surf_C[1]*f_rec_lo[12]+f_rec_lo[1]*D_surf_C[12]+D_surf_C[7]*f_rec_lo[10]+f_rec_lo[7]*D_surf_C[10]+D_surf_C[0]*f_rec_lo[9]+f_rec_lo[0]*D_surf_C[9]+D_surf_C[5]*f_rec_lo[8]+f_rec_lo[5]*D_surf_C[8]+D_surf_C[2]*f_rec_lo[4]+f_rec_lo[2]*D_surf_C[4]); 
  surft2_lo[10] = 0.25*(D_surf_C[5]*f_rec_lo[15]+f_rec_lo[5]*D_surf_C[15]+D_surf_C[2]*f_rec_lo[14]+f_rec_lo[2]*D_surf_C[14]+D_surf_C[1]*f_rec_lo[13]+f_rec_lo[1]*D_surf_C[13]+D_surf_C[11]*f_rec_lo[12]+f_rec_lo[11]*D_surf_C[12]+D_surf_C[0]*f_rec_lo[10]+f_rec_lo[0]*D_surf_C[10]+D_surf_C[7]*f_rec_lo[9]+f_rec_lo[7]*D_surf_C[9]+D_surf_C[6]*f_rec_lo[8]+f_rec_lo[6]*D_surf_C[8]+D_surf_C[3]*f_rec_lo[4]+f_rec_lo[3]*D_surf_C[4]); 
  surft2_lo[11] = 0.25*(D_surf_C[4]*f_rec_lo[15]+f_rec_lo[4]*D_surf_C[15]+D_surf_C[8]*f_rec_lo[14]+f_rec_lo[8]*D_surf_C[14]+D_surf_C[9]*f_rec_lo[13]+f_rec_lo[9]*D_surf_C[13]+D_surf_C[10]*f_rec_lo[12]+f_rec_lo[10]*D_surf_C[12]+D_surf_C[0]*f_rec_lo[11]+f_rec_lo[0]*D_surf_C[11]+D_surf_C[1]*f_rec_lo[7]+f_rec_lo[1]*D_surf_C[7]+D_surf_C[2]*f_rec_lo[6]+f_rec_lo[2]*D_surf_C[6]+D_surf_C[3]*f_rec_lo[5]+f_rec_lo[3]*D_surf_C[5]); 
  surft2_lo[12] = 0.25*(D_surf_C[3]*f_rec_lo[15]+f_rec_lo[3]*D_surf_C[15]+D_surf_C[6]*f_rec_lo[14]+f_rec_lo[6]*D_surf_C[14]+D_surf_C[7]*f_rec_lo[13]+f_rec_lo[7]*D_surf_C[13]+D_surf_C[0]*f_rec_lo[12]+f_rec_lo[0]*D_surf_C[12]+D_surf_C[10]*f_rec_lo[11]+f_rec_lo[10]*D_surf_C[11]+D_surf_C[1]*f_rec_lo[9]+f_rec_lo[1]*D_surf_C[9]+D_surf_C[2]*f_rec_lo[8]+f_rec_lo[2]*D_surf_C[8]+D_surf_C[4]*f_rec_lo[5]+f_rec_lo[4]*D_surf_C[5]); 
  surft2_lo[13] = 0.25*(D_surf_C[2]*f_rec_lo[15]+f_rec_lo[2]*D_surf_C[15]+D_surf_C[5]*f_rec_lo[14]+f_rec_lo[5]*D_surf_C[14]+D_surf_C[0]*f_rec_lo[13]+f_rec_lo[0]*D_surf_C[13]+D_surf_C[7]*f_rec_lo[12]+f_rec_lo[7]*D_surf_C[12]+D_surf_C[9]*f_rec_lo[11]+f_rec_lo[9]*D_surf_C[11]+D_surf_C[1]*f_rec_lo[10]+f_rec_lo[1]*D_surf_C[10]+D_surf_C[3]*f_rec_lo[8]+f_rec_lo[3]*D_surf_C[8]+D_surf_C[4]*f_rec_lo[6]+f_rec_lo[4]*D_surf_C[6]); 
  surft2_lo[14] = 0.25*(D_surf_C[1]*f_rec_lo[15]+f_rec_lo[1]*D_surf_C[15]+D_surf_C[0]*f_rec_lo[14]+f_rec_lo[0]*D_surf_C[14]+D_surf_C[5]*f_rec_lo[13]+f_rec_lo[5]*D_surf_C[13]+D_surf_C[6]*f_rec_lo[12]+f_rec_lo[6]*D_surf_C[12]+D_surf_C[8]*f_rec_lo[11]+f_rec_lo[8]*D_surf_C[11]+D_surf_C[2]*f_rec_lo[10]+f_rec_lo[2]*D_surf_C[10]+D_surf_C[3]*f_rec_lo[9]+f_rec_lo[3]*D_surf_C[9]+D_surf_C[4]*f_rec_lo[7]+f_rec_lo[4]*D_surf_C[7]); 
  surft2_lo[15] = 0.25*(D_surf_C[0]*f_rec_lo[15]+f_rec_lo[0]*D_surf_C[15]+D_surf_C[1]*f_rec_lo[14]+f_rec_lo[1]*D_surf_C[14]+D_surf_C[2]*f_rec_lo[13]+f_rec_lo[2]*D_surf_C[13]+D_surf_C[3]*f_rec_lo[12]+f_rec_lo[3]*D_surf_C[12]+D_surf_C[4]*f_rec_lo[11]+f_rec_lo[4]*D_surf_C[11]+D_surf_C[5]*f_rec_lo[10]+f_rec_lo[5]*D_surf_C[10]+D_surf_C[6]*f_rec_lo[9]+f_rec_lo[6]*D_surf_C[9]+D_surf_C[7]*f_rec_lo[8]+f_rec_lo[7]*D_surf_C[8]); 
  surft2_up[0] = 0.25*(D_surf_R[15]*f_rec_up[15]+D_surf_R[14]*f_rec_up[14]+D_surf_R[13]*f_rec_up[13]+D_surf_R[12]*f_rec_up[12]+D_surf_R[11]*f_rec_up[11]+D_surf_R[10]*f_rec_up[10]+D_surf_R[9]*f_rec_up[9]+D_surf_R[8]*f_rec_up[8]+D_surf_R[7]*f_rec_up[7]+D_surf_R[6]*f_rec_up[6]+D_surf_R[5]*f_rec_up[5]+D_surf_R[4]*f_rec_up[4]+D_surf_R[3]*f_rec_up[3]+D_surf_R[2]*f_rec_up[2]+D_surf_R[1]*f_rec_up[1]+D_surf_R[0]*f_rec_up[0]); 
  surft2_up[1] = 0.25*(D_surf_R[14]*f_rec_up[15]+f_rec_up[14]*D_surf_R[15]+D_surf_R[10]*f_rec_up[13]+f_rec_up[10]*D_surf_R[13]+D_surf_R[9]*f_rec_up[12]+f_rec_up[9]*D_surf_R[12]+D_surf_R[7]*f_rec_up[11]+f_rec_up[7]*D_surf_R[11]+D_surf_R[4]*f_rec_up[8]+f_rec_up[4]*D_surf_R[8]+D_surf_R[3]*f_rec_up[6]+f_rec_up[3]*D_surf_R[6]+D_surf_R[2]*f_rec_up[5]+f_rec_up[2]*D_surf_R[5]+D_surf_R[0]*f_rec_up[1]+f_rec_up[0]*D_surf_R[1]); 
  surft2_up[2] = 0.25*(D_surf_R[13]*f_rec_up[15]+f_rec_up[13]*D_surf_R[15]+D_surf_R[10]*f_rec_up[14]+f_rec_up[10]*D_surf_R[14]+D_surf_R[8]*f_rec_up[12]+f_rec_up[8]*D_surf_R[12]+D_surf_R[6]*f_rec_up[11]+f_rec_up[6]*D_surf_R[11]+D_surf_R[4]*f_rec_up[9]+f_rec_up[4]*D_surf_R[9]+D_surf_R[3]*f_rec_up[7]+f_rec_up[3]*D_surf_R[7]+D_surf_R[1]*f_rec_up[5]+f_rec_up[1]*D_surf_R[5]+D_surf_R[0]*f_rec_up[2]+f_rec_up[0]*D_surf_R[2]); 
  surft2_up[3] = 0.25*(D_surf_R[12]*f_rec_up[15]+f_rec_up[12]*D_surf_R[15]+D_surf_R[9]*f_rec_up[14]+f_rec_up[9]*D_surf_R[14]+D_surf_R[8]*f_rec_up[13]+f_rec_up[8]*D_surf_R[13]+D_surf_R[5]*f_rec_up[11]+f_rec_up[5]*D_surf_R[11]+D_surf_R[4]*f_rec_up[10]+f_rec_up[4]*D_surf_R[10]+D_surf_R[2]*f_rec_up[7]+f_rec_up[2]*D_surf_R[7]+D_surf_R[1]*f_rec_up[6]+f_rec_up[1]*D_surf_R[6]+D_surf_R[0]*f_rec_up[3]+f_rec_up[0]*D_surf_R[3]); 
  surft2_up[4] = 0.25*(D_surf_R[11]*f_rec_up[15]+f_rec_up[11]*D_surf_R[15]+D_surf_R[7]*f_rec_up[14]+f_rec_up[7]*D_surf_R[14]+D_surf_R[6]*f_rec_up[13]+f_rec_up[6]*D_surf_R[13]+D_surf_R[5]*f_rec_up[12]+f_rec_up[5]*D_surf_R[12]+D_surf_R[3]*f_rec_up[10]+f_rec_up[3]*D_surf_R[10]+D_surf_R[2]*f_rec_up[9]+f_rec_up[2]*D_surf_R[9]+D_surf_R[1]*f_rec_up[8]+f_rec_up[1]*D_surf_R[8]+D_surf_R[0]*f_rec_up[4]+f_rec_up[0]*D_surf_R[4]); 
  surft2_up[5] = 0.25*(D_surf_R[10]*f_rec_up[15]+f_rec_up[10]*D_surf_R[15]+D_surf_R[13]*f_rec_up[14]+f_rec_up[13]*D_surf_R[14]+D_surf_R[4]*f_rec_up[12]+f_rec_up[4]*D_surf_R[12]+D_surf_R[3]*f_rec_up[11]+f_rec_up[3]*D_surf_R[11]+D_surf_R[8]*f_rec_up[9]+f_rec_up[8]*D_surf_R[9]+D_surf_R[6]*f_rec_up[7]+f_rec_up[6]*D_surf_R[7]+D_surf_R[0]*f_rec_up[5]+f_rec_up[0]*D_surf_R[5]+D_surf_R[1]*f_rec_up[2]+f_rec_up[1]*D_surf_R[2]); 
  surft2_up[6] = 0.25*(D_surf_R[9]*f_rec_up[15]+f_rec_up[9]*D_surf_R[15]+D_surf_R[12]*f_rec_up[14]+f_rec_up[12]*D_surf_R[14]+D_surf_R[4]*f_rec_up[13]+f_rec_up[4]*D_surf_R[13]+D_surf_R[2]*f_rec_up[11]+f_rec_up[2]*D_surf_R[11]+D_surf_R[8]*f_rec_up[10]+f_rec_up[8]*D_surf_R[10]+D_surf_R[5]*f_rec_up[7]+f_rec_up[5]*D_surf_R[7]+D_surf_R[0]*f_rec_up[6]+f_rec_up[0]*D_surf_R[6]+D_surf_R[1]*f_rec_up[3]+f_rec_up[1]*D_surf_R[3]); 
  surft2_up[7] = 0.25*(D_surf_R[8]*f_rec_up[15]+f_rec_up[8]*D_surf_R[15]+D_surf_R[4]*f_rec_up[14]+f_rec_up[4]*D_surf_R[14]+D_surf_R[12]*f_rec_up[13]+f_rec_up[12]*D_surf_R[13]+D_surf_R[1]*f_rec_up[11]+f_rec_up[1]*D_surf_R[11]+D_surf_R[9]*f_rec_up[10]+f_rec_up[9]*D_surf_R[10]+D_surf_R[0]*f_rec_up[7]+f_rec_up[0]*D_surf_R[7]+D_surf_R[5]*f_rec_up[6]+f_rec_up[5]*D_surf_R[6]+D_surf_R[2]*f_rec_up[3]+f_rec_up[2]*D_surf_R[3]); 
  surft2_up[8] = 0.25*(D_surf_R[7]*f_rec_up[15]+f_rec_up[7]*D_surf_R[15]+D_surf_R[11]*f_rec_up[14]+f_rec_up[11]*D_surf_R[14]+D_surf_R[3]*f_rec_up[13]+f_rec_up[3]*D_surf_R[13]+D_surf_R[2]*f_rec_up[12]+f_rec_up[2]*D_surf_R[12]+D_surf_R[6]*f_rec_up[10]+f_rec_up[6]*D_surf_R[10]+D_surf_R[5]*f_rec_up[9]+f_rec_up[5]*D_surf_R[9]+D_surf_R[0]*f_rec_up[8]+f_rec_up[0]*D_surf_R[8]+D_surf_R[1]*f_rec_up[4]+f_rec_up[1]*D_surf_R[4]); 
  surft2_up[9] = 0.25*(D_surf_R[6]*f_rec_up[15]+f_rec_up[6]*D_surf_R[15]+D_surf_R[3]*f_rec_up[14]+f_rec_up[3]*D_surf_R[14]+D_surf_R[11]*f_rec_up[13]+f_rec_up[11]*D_surf_R[13]+D_surf_R[1]*f_rec_up[12]+f_rec_up[1]*D_surf_R[12]+D_surf_R[7]*f_rec_up[10]+f_rec_up[7]*D_surf_R[10]+D_surf_R[0]*f_rec_up[9]+f_rec_up[0]*D_surf_R[9]+D_surf_R[5]*f_rec_up[8]+f_rec_up[5]*D_surf_R[8]+D_surf_R[2]*f_rec_up[4]+f_rec_up[2]*D_surf_R[4]); 
  surft2_up[10] = 0.25*(D_surf_R[5]*f_rec_up[15]+f_rec_up[5]*D_surf_R[15]+D_surf_R[2]*f_rec_up[14]+f_rec_up[2]*D_surf_R[14]+D_surf_R[1]*f_rec_up[13]+f_rec_up[1]*D_surf_R[13]+D_surf_R[11]*f_rec_up[12]+f_rec_up[11]*D_surf_R[12]+D_surf_R[0]*f_rec_up[10]+f_rec_up[0]*D_surf_R[10]+D_surf_R[7]*f_rec_up[9]+f_rec_up[7]*D_surf_R[9]+D_surf_R[6]*f_rec_up[8]+f_rec_up[6]*D_surf_R[8]+D_surf_R[3]*f_rec_up[4]+f_rec_up[3]*D_surf_R[4]); 
  surft2_up[11] = 0.25*(D_surf_R[4]*f_rec_up[15]+f_rec_up[4]*D_surf_R[15]+D_surf_R[8]*f_rec_up[14]+f_rec_up[8]*D_surf_R[14]+D_surf_R[9]*f_rec_up[13]+f_rec_up[9]*D_surf_R[13]+D_surf_R[10]*f_rec_up[12]+f_rec_up[10]*D_surf_R[12]+D_surf_R[0]*f_rec_up[11]+f_rec_up[0]*D_surf_R[11]+D_surf_R[1]*f_rec_up[7]+f_rec_up[1]*D_surf_R[7]+D_surf_R[2]*f_rec_up[6]+f_rec_up[2]*D_surf_R[6]+D_surf_R[3]*f_rec_up[5]+f_rec_up[3]*D_surf_R[5]); 
  surft2_up[12] = 0.25*(D_surf_R[3]*f_rec_up[15]+f_rec_up[3]*D_surf_R[15]+D_surf_R[6]*f_rec_up[14]+f_rec_up[6]*D_surf_R[14]+D_surf_R[7]*f_rec_up[13]+f_rec_up[7]*D_surf_R[13]+D_surf_R[0]*f_rec_up[12]+f_rec_up[0]*D_surf_R[12]+D_surf_R[10]*f_rec_up[11]+f_rec_up[10]*D_surf_R[11]+D_surf_R[1]*f_rec_up[9]+f_rec_up[1]*D_surf_R[9]+D_surf_R[2]*f_rec_up[8]+f_rec_up[2]*D_surf_R[8]+D_surf_R[4]*f_rec_up[5]+f_rec_up[4]*D_surf_R[5]); 
  surft2_up[13] = 0.25*(D_surf_R[2]*f_rec_up[15]+f_rec_up[2]*D_surf_R[15]+D_surf_R[5]*f_rec_up[14]+f_rec_up[5]*D_surf_R[14]+D_surf_R[0]*f_rec_up[13]+f_rec_up[0]*D_surf_R[13]+D_surf_R[7]*f_rec_up[12]+f_rec_up[7]*D_surf_R[12]+D_surf_R[9]*f_rec_up[11]+f_rec_up[9]*D_surf_R[11]+D_surf_R[1]*f_rec_up[10]+f_rec_up[1]*D_surf_R[10]+D_surf_R[3]*f_rec_up[8]+f_rec_up[3]*D_surf_R[8]+D_surf_R[4]*f_rec_up[6]+f_rec_up[4]*D_surf_R[6]); 
  surft2_up[14] = 0.25*(D_surf_R[1]*f_rec_up[15]+f_rec_up[1]*D_surf_R[15]+D_surf_R[0]*f_rec_up[14]+f_rec_up[0]*D_surf_R[14]+D_surf_R[5]*f_rec_up[13]+f_rec_up[5]*D_surf_R[13]+D_surf_R[6]*f_rec_up[12]+f_rec_up[6]*D_surf_R[12]+D_surf_R[8]*f_rec_up[11]+f_rec_up[8]*D_surf_R[11]+D_surf_R[2]*f_rec_up[10]+f_rec_up[2]*D_surf_R[10]+D_surf_R[3]*f_rec_up[9]+f_rec_up[3]*D_surf_R[9]+D_surf_R[4]*f_rec_up[7]+f_rec_up[4]*D_surf_R[7]); 
  surft2_up[15] = 0.25*(D_surf_R[0]*f_rec_up[15]+f_rec_up[0]*D_surf_R[15]+D_surf_R[1]*f_rec_up[14]+f_rec_up[1]*D_surf_R[14]+D_surf_R[2]*f_rec_up[13]+f_rec_up[2]*D_surf_R[13]+D_surf_R[3]*f_rec_up[12]+f_rec_up[3]*D_surf_R[12]+D_surf_R[4]*f_rec_up[11]+f_rec_up[4]*D_surf_R[11]+D_surf_R[5]*f_rec_up[10]+f_rec_up[5]*D_surf_R[10]+D_surf_R[6]*f_rec_up[9]+f_rec_up[6]*D_surf_R[9]+D_surf_R[7]*f_rec_up[8]+f_rec_up[7]*D_surf_R[8]); 

  vol[5] = 0.5303300858899105*(fC[26]*DC[31]+fC[19]*DC[30]+fC[18]*DC[29]+fC[17]*DC[28]+fC[16]*DC[27]+fC[11]*DC[25]+fC[10]*DC[24]+fC[9]*DC[23]+fC[8]*DC[22]+fC[7]*DC[21]+fC[6]*DC[20]+fC[4]*DC[15]+fC[3]*DC[14]+fC[2]*DC[13]+fC[1]*DC[12]+fC[0]*DC[5]); 
  vol[12] = 0.5303300858899105*(fC[19]*DC[31]+fC[26]*DC[30]+fC[11]*DC[29]+fC[10]*DC[28]+fC[8]*DC[27]+fC[18]*DC[25]+fC[17]*DC[24]+fC[4]*DC[23]+fC[16]*DC[22]+fC[3]*DC[21]+fC[2]*DC[20]+fC[9]*DC[15]+fC[7]*DC[14]+fC[6]*DC[13]+fC[0]*DC[12]+fC[1]*DC[5]); 
  vol[13] = 0.5303300858899105*(fC[18]*DC[31]+fC[11]*DC[30]+fC[26]*DC[29]+fC[9]*DC[28]+fC[7]*DC[27]+fC[19]*DC[25]+fC[4]*DC[24]+fC[17]*DC[23]+fC[3]*DC[22]+fC[16]*DC[21]+fC[1]*DC[20]+fC[10]*DC[15]+fC[8]*DC[14]+fC[0]*DC[13]+fC[6]*DC[12]+fC[2]*DC[5]); 
  vol[14] = 0.5303300858899105*(fC[17]*DC[31]+fC[10]*DC[30]+fC[9]*DC[29]+fC[26]*DC[28]+fC[6]*DC[27]+fC[4]*DC[25]+fC[19]*DC[24]+fC[18]*DC[23]+fC[2]*DC[22]+fC[1]*DC[21]+fC[16]*DC[20]+fC[11]*DC[15]+fC[0]*DC[14]+fC[8]*DC[13]+fC[7]*DC[12]+fC[3]*DC[5]); 
  vol[15] = 0.5303300858899105*(fC[16]*DC[31]+fC[8]*DC[30]+fC[7]*DC[29]+fC[6]*DC[28]+fC[26]*DC[27]+fC[3]*DC[25]+fC[2]*DC[24]+fC[1]*DC[23]+fC[19]*DC[22]+fC[18]*DC[21]+fC[17]*DC[20]+fC[0]*DC[15]+fC[11]*DC[14]+fC[10]*DC[13]+fC[9]*DC[12]+fC[4]*DC[5]); 
  vol[20] = 0.5303300858899105*(fC[11]*DC[31]+fC[18]*DC[30]+fC[19]*DC[29]+fC[4]*DC[28]+fC[3]*DC[27]+DC[25]*fC[26]+fC[9]*DC[24]+fC[10]*DC[23]+fC[7]*DC[22]+fC[8]*DC[21]+fC[0]*DC[20]+DC[15]*fC[17]+DC[14]*fC[16]+fC[1]*DC[13]+fC[2]*DC[12]+DC[5]*fC[6]); 
  vol[21] = 0.5303300858899105*(fC[10]*DC[31]+fC[17]*DC[30]+fC[4]*DC[29]+fC[19]*DC[28]+fC[2]*DC[27]+DC[24]*fC[26]+fC[9]*DC[25]+fC[11]*DC[23]+fC[6]*DC[22]+fC[0]*DC[21]+fC[8]*DC[20]+DC[15]*fC[18]+DC[13]*fC[16]+fC[1]*DC[14]+fC[3]*DC[12]+DC[5]*fC[7]); 
  vol[22] = 0.5303300858899105*(fC[9]*DC[31]+fC[4]*DC[30]+fC[17]*DC[29]+fC[18]*DC[28]+fC[1]*DC[27]+DC[23]*fC[26]+fC[10]*DC[25]+fC[11]*DC[24]+fC[0]*DC[22]+fC[6]*DC[21]+fC[7]*DC[20]+DC[15]*fC[19]+DC[12]*fC[16]+fC[2]*DC[14]+fC[3]*DC[13]+DC[5]*fC[8]); 
  vol[23] = 0.5303300858899105*(fC[8]*DC[31]+fC[16]*DC[30]+fC[3]*DC[29]+fC[2]*DC[28]+fC[19]*DC[27]+DC[22]*fC[26]+fC[7]*DC[25]+fC[6]*DC[24]+fC[0]*DC[23]+fC[11]*DC[21]+fC[10]*DC[20]+DC[14]*fC[18]+DC[13]*fC[17]+fC[1]*DC[15]+fC[4]*DC[12]+DC[5]*fC[9]); 
  vol[24] = 0.5303300858899105*(fC[7]*DC[31]+fC[3]*DC[30]+fC[16]*DC[29]+fC[1]*DC[28]+fC[18]*DC[27]+DC[21]*fC[26]+fC[8]*DC[25]+fC[0]*DC[24]+fC[6]*DC[23]+fC[11]*DC[22]+fC[9]*DC[20]+DC[14]*fC[19]+DC[12]*fC[17]+fC[2]*DC[15]+fC[4]*DC[13]+DC[5]*fC[10]); 
  vol[25] = 0.5303300858899105*(fC[6]*DC[31]+fC[2]*DC[30]+fC[1]*DC[29]+fC[16]*DC[28]+fC[17]*DC[27]+DC[20]*fC[26]+fC[0]*DC[25]+fC[8]*DC[24]+fC[7]*DC[23]+fC[10]*DC[22]+fC[9]*DC[21]+DC[13]*fC[19]+DC[12]*fC[18]+fC[3]*DC[15]+fC[4]*DC[14]+DC[5]*fC[11]); 
  vol[27] = 0.5303300858899105*(fC[4]*DC[31]+fC[9]*DC[30]+fC[10]*DC[29]+fC[11]*DC[28]+fC[0]*DC[27]+DC[15]*fC[26]+fC[17]*DC[25]+fC[18]*DC[24]+fC[19]*DC[23]+fC[1]*DC[22]+fC[2]*DC[21]+fC[3]*DC[20]+DC[5]*fC[16]+fC[6]*DC[14]+fC[7]*DC[13]+fC[8]*DC[12]); 
  vol[28] = 0.5303300858899105*(fC[3]*DC[31]+fC[7]*DC[30]+fC[8]*DC[29]+fC[0]*DC[28]+fC[11]*DC[27]+DC[14]*fC[26]+fC[16]*DC[25]+fC[1]*DC[24]+fC[2]*DC[23]+fC[18]*DC[22]+fC[19]*DC[21]+fC[4]*DC[20]+DC[5]*fC[17]+fC[6]*DC[15]+fC[9]*DC[13]+fC[10]*DC[12]); 
  vol[29] = 0.5303300858899105*(fC[2]*DC[31]+fC[6]*DC[30]+fC[0]*DC[29]+fC[8]*DC[28]+fC[10]*DC[27]+DC[13]*fC[26]+fC[1]*DC[25]+fC[16]*DC[24]+fC[3]*DC[23]+fC[17]*DC[22]+fC[4]*DC[21]+fC[19]*DC[20]+DC[5]*fC[18]+fC[7]*DC[15]+fC[9]*DC[14]+fC[11]*DC[12]); 
  vol[30] = 0.5303300858899105*(fC[1]*DC[31]+fC[0]*DC[30]+fC[6]*DC[29]+fC[7]*DC[28]+fC[9]*DC[27]+DC[12]*fC[26]+fC[2]*DC[25]+fC[3]*DC[24]+fC[16]*DC[23]+fC[4]*DC[22]+fC[17]*DC[21]+fC[18]*DC[20]+DC[5]*fC[19]+fC[8]*DC[15]+fC[10]*DC[14]+fC[11]*DC[13]); 
  vol[31] = 0.5303300858899105*(fC[0]*DC[31]+fC[1]*DC[30]+fC[2]*DC[29]+fC[3]*DC[28]+fC[4]*DC[27]+DC[5]*fC[26]+fC[6]*DC[25]+fC[7]*DC[24]+fC[8]*DC[23]+fC[9]*DC[22]+fC[10]*DC[21]+fC[11]*DC[20]+DC[12]*fC[19]+DC[13]*fC[18]+DC[14]*fC[17]+DC[15]*fC[16]); 

  out[0] += (0.5*vol[0]+0.35355339059327373*surft1_up[0]-0.35355339059327373*surft1_lo[0])*dv_inv_sq; 
  out[1] += (0.5*vol[1]+0.35355339059327373*surft1_up[1]-0.35355339059327373*surft1_lo[1])*dv_inv_sq; 
  out[2] += (0.5*vol[2]+0.35355339059327373*surft1_up[2]-0.35355339059327373*surft1_lo[2])*dv_inv_sq; 
  out[3] += (0.5*vol[3]+0.35355339059327373*surft1_up[3]-0.35355339059327373*surft1_lo[3])*dv_inv_sq; 
  out[4] += (0.5*vol[4]+0.35355339059327373*surft1_up[4]-0.35355339059327373*surft1_lo[4])*dv_inv_sq; 
  out[5] += (0.5*vol[5]-0.6123724356957945*surft2_up[0]+0.6123724356957945*(surft2_lo[0]+surft1_up[0]+surft1_lo[0]))*dv_inv_sq; 
  out[6] += (0.5*vol[6]+0.35355339059327373*surft1_up[5]-0.35355339059327373*surft1_lo[5])*dv_inv_sq; 
  out[7] += (0.5*vol[7]+0.35355339059327373*surft1_up[6]-0.35355339059327373*surft1_lo[6])*dv_inv_sq; 
  out[8] += (0.5*vol[8]+0.35355339059327373*surft1_up[7]-0.35355339059327373*surft1_lo[7])*dv_inv_sq; 
  out[9] += (0.5*vol[9]+0.35355339059327373*surft1_up[8]-0.35355339059327373*surft1_lo[8])*dv_inv_sq; 
  out[10] += (0.5*vol[10]+0.35355339059327373*surft1_up[9]-0.35355339059327373*surft1_lo[9])*dv_inv_sq; 
  out[11] += (0.5*vol[11]+0.35355339059327373*surft1_up[10]-0.35355339059327373*surft1_lo[10])*dv_inv_sq; 
  out[12] += (0.5*vol[12]-0.6123724356957945*surft2_up[1]+0.6123724356957945*(surft2_lo[1]+surft1_up[1]+surft1_lo[1]))*dv_inv_sq; 
  out[13] += (0.5*vol[13]-0.6123724356957945*surft2_up[2]+0.6123724356957945*(surft2_lo[2]+surft1_up[2]+surft1_lo[2]))*dv_inv_sq; 
  out[14] += (0.5*vol[14]-0.6123724356957945*surft2_up[3]+0.6123724356957945*(surft2_lo[3]+surft1_up[3]+surft1_lo[3]))*dv_inv_sq; 
  out[15] += (0.5*vol[15]-0.6123724356957945*surft2_up[4]+0.6123724356957945*(surft2_lo[4]+surft1_up[4]+surft1_lo[4]))*dv_inv_sq; 
  out[16] += (0.5*vol[16]+0.35355339059327373*surft1_up[11]-0.35355339059327373*surft1_lo[11])*dv_inv_sq; 
  out[17] += (0.5*vol[17]+0.35355339059327373*surft1_up[12]-0.35355339059327373*surft1_lo[12])*dv_inv_sq; 
  out[18] += (0.5*vol[18]+0.35355339059327373*surft1_up[13]-0.35355339059327373*surft1_lo[13])*dv_inv_sq; 
  out[19] += (0.5*vol[19]+0.35355339059327373*surft1_up[14]-0.35355339059327373*surft1_lo[14])*dv_inv_sq; 
  out[20] += (0.5*vol[20]-0.6123724356957945*surft2_up[5]+0.6123724356957945*(surft2_lo[5]+surft1_up[5]+surft1_lo[5]))*dv_inv_sq; 
  out[21] += (0.5*vol[21]-0.6123724356957945*surft2_up[6]+0.6123724356957945*(surft2_lo[6]+surft1_up[6]+surft1_lo[6]))*dv_inv_sq; 
  out[22] += (0.5*vol[22]-0.6123724356957945*surft2_up[7]+0.6123724356957945*(surft2_lo[7]+surft1_up[7]+surft1_lo[7]))*dv_inv_sq; 
  out[23] += (0.5*vol[23]-0.6123724356957945*surft2_up[8]+0.6123724356957945*(surft2_lo[8]+surft1_up[8]+surft1_lo[8]))*dv_inv_sq; 
  out[24] += (0.5*vol[24]-0.6123724356957945*surft2_up[9]+0.6123724356957945*(surft2_lo[9]+surft1_up[9]+surft1_lo[9]))*dv_inv_sq; 
  out[25] += (0.5*vol[25]-0.6123724356957945*surft2_up[10]+0.6123724356957945*(surft2_lo[10]+surft1_up[10]+surft1_lo[10]))*dv_inv_sq; 
  out[26] += (0.5*vol[26]+0.35355339059327373*surft1_up[15]-0.35355339059327373*surft1_lo[15])*dv_inv_sq; 
  out[27] += (0.5*vol[27]-0.6123724356957945*surft2_up[11]+0.6123724356957945*(surft2_lo[11]+surft1_up[11]+surft1_lo[11]))*dv_inv_sq; 
  out[28] += (0.5*vol[28]-0.6123724356957945*surft2_up[12]+0.6123724356957945*(surft2_lo[12]+surft1_up[12]+surft1_lo[12]))*dv_inv_sq; 
  out[29] += (0.5*vol[29]-0.6123724356957945*surft2_up[13]+0.6123724356957945*(surft2_lo[13]+surft1_up[13]+surft1_lo[13]))*dv_inv_sq; 
  out[30] += (0.5*vol[30]-0.6123724356957945*surft2_up[14]+0.6123724356957945*(surft2_lo[14]+surft1_up[14]+surft1_lo[14]))*dv_inv_sq; 
  out[31] += (0.5*vol[31]-0.6123724356957945*surft2_up[15]+0.6123724356957945*(surft2_lo[15]+surft1_up[15]+surft1_lo[15]))*dv_inv_sq; 

  double cflFreq = fmax(fabs(D_surf_C[0]), fabs(D_surf_R[0])); 

  return 1.0*dv_inv_sq*cflFreq; 
} 
