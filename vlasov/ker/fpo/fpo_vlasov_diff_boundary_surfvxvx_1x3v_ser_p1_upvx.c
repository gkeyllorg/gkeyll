#include <gkyl_fpo_vlasov_kernels.h> 
 
GKYL_CU_DH double fpo_vlasov_diff_boundary_surfvxvx_1x3v_ser_p1_upvx(const double *dxv, const double *diff_coeff_C,
      const double *diff_coeff_surf_stencil[9], const double* f_stencil[9], double* out) { 
  // dxv[NDIM]: Cell spacing in each direction. 
  // diff_coeff_stencil[9]: 9-cell stencil of diffusion tensor. 
  // f_stencil[9]: 9-cell stencil of distribution function. 
  // out: Incremented output. 


  double dv_inv_sq = 4.0/dxv[1]/dxv[1]; 
 
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
  const double *DC = &diff_coeff_C[0]; 
  const double* fL = f_stencil[0]; 
  const double *D_surf_C = &diff_coeff_surf_stencil[1][0]; 
  const double* fC = f_stencil[1]; 
  const double *D_surf_R = &diff_coeff_surf_stencil[2][0]; 
  const double* fR = f_stencil[2]; 

  f_rec_lo[0] = 0.408248290463863*fL[2]-0.408248290463863*fC[2]+0.3535533905932737*(fL[0]+fC[0]); 
  f_rec_lo[1] = 0.408248290463863*fL[5]-0.408248290463863*fC[5]+0.3535533905932737*(fL[1]+fC[1]); 
  f_rec_lo[2] = 0.408248290463863*fL[7]-0.408248290463863*fC[7]+0.3535533905932737*(fL[3]+fC[3]); 
  f_rec_lo[3] = 0.408248290463863*fL[9]-0.408248290463863*fC[9]+0.3535533905932737*(fL[4]+fC[4]); 
  f_rec_lo[4] = 0.408248290463863*fL[11]-0.408248290463863*fC[11]+0.3535533905932737*(fL[6]+fC[6]); 
  f_rec_lo[5] = 0.408248290463863*fL[12]-0.408248290463863*fC[12]+0.3535533905932737*(fL[8]+fC[8]); 
  f_rec_lo[6] = 0.408248290463863*fL[14]-0.408248290463863*fC[14]+0.3535533905932737*(fL[10]+fC[10]); 
  f_rec_lo[7] = 0.408248290463863*fL[15]-0.408248290463863*fC[15]+0.3535533905932737*(fL[13]+fC[13]); 
  f_rec_up[0] = 1.224744871391589*fC[2]+0.7071067811865475*fC[0]; 
  f_rec_up[1] = 1.224744871391589*fC[5]+0.7071067811865475*fC[1]; 
  f_rec_up[2] = 1.224744871391589*fC[7]+0.7071067811865475*fC[3]; 
  f_rec_up[3] = 1.224744871391589*fC[9]+0.7071067811865475*fC[4]; 
  f_rec_up[4] = 1.224744871391589*fC[11]+0.7071067811865475*fC[6]; 
  f_rec_up[5] = 1.224744871391589*fC[12]+0.7071067811865475*fC[8]; 
  f_rec_up[6] = 1.224744871391589*fC[14]+0.7071067811865475*fC[10]; 
  f_rec_up[7] = 1.224744871391589*fC[15]+0.7071067811865475*fC[13]; 

  df_rec_lo[0] = -(0.7654655446197428*(fL[2]+fC[2]))-0.7954951288348656*fL[0]+0.7954951288348656*fC[0]; 
  df_rec_lo[1] = -(0.7654655446197428*(fL[5]+fC[5]))-0.7954951288348656*fL[1]+0.7954951288348656*fC[1]; 
  df_rec_lo[2] = -(0.7654655446197428*(fL[7]+fC[7]))-0.7954951288348656*fL[3]+0.7954951288348656*fC[3]; 
  df_rec_lo[3] = -(0.7654655446197428*(fL[9]+fC[9]))-0.7954951288348656*fL[4]+0.7954951288348656*fC[4]; 
  df_rec_lo[4] = -(0.7654655446197428*(fL[11]+fC[11]))-0.7954951288348656*fL[6]+0.7954951288348656*fC[6]; 
  df_rec_lo[5] = -(0.7654655446197428*(fL[12]+fC[12]))-0.7954951288348656*fL[8]+0.7954951288348656*fC[8]; 
  df_rec_lo[6] = -(0.7654655446197428*(fL[14]+fC[14]))-0.7954951288348656*fL[10]+0.7954951288348656*fC[10]; 
  df_rec_lo[7] = -(0.7654655446197428*(fL[15]+fC[15]))-0.7954951288348656*fL[13]+0.7954951288348656*fC[13]; 

  surft1_lo[0] = 0.3535533905932737*(D_surf_C[7]*df_rec_lo[7]+D_surf_C[6]*df_rec_lo[6]+D_surf_C[5]*df_rec_lo[5]+D_surf_C[4]*df_rec_lo[4]+D_surf_C[3]*df_rec_lo[3]+D_surf_C[2]*df_rec_lo[2]+D_surf_C[1]*df_rec_lo[1]+D_surf_C[0]*df_rec_lo[0]); 
  surft1_lo[1] = 0.3535533905932737*(D_surf_C[6]*df_rec_lo[7]+df_rec_lo[6]*D_surf_C[7]+D_surf_C[3]*df_rec_lo[5]+df_rec_lo[3]*D_surf_C[5]+D_surf_C[2]*df_rec_lo[4]+df_rec_lo[2]*D_surf_C[4]+D_surf_C[0]*df_rec_lo[1]+df_rec_lo[0]*D_surf_C[1]); 
  surft1_lo[2] = 0.3535533905932737*(D_surf_C[5]*df_rec_lo[7]+df_rec_lo[5]*D_surf_C[7]+D_surf_C[3]*df_rec_lo[6]+df_rec_lo[3]*D_surf_C[6]+D_surf_C[1]*df_rec_lo[4]+df_rec_lo[1]*D_surf_C[4]+D_surf_C[0]*df_rec_lo[2]+df_rec_lo[0]*D_surf_C[2]); 
  surft1_lo[3] = 0.3535533905932737*(D_surf_C[4]*df_rec_lo[7]+df_rec_lo[4]*D_surf_C[7]+D_surf_C[2]*df_rec_lo[6]+df_rec_lo[2]*D_surf_C[6]+D_surf_C[1]*df_rec_lo[5]+df_rec_lo[1]*D_surf_C[5]+D_surf_C[0]*df_rec_lo[3]+df_rec_lo[0]*D_surf_C[3]); 
  surft1_lo[4] = 0.3535533905932737*(D_surf_C[3]*df_rec_lo[7]+df_rec_lo[3]*D_surf_C[7]+D_surf_C[5]*df_rec_lo[6]+df_rec_lo[5]*D_surf_C[6]+D_surf_C[0]*df_rec_lo[4]+df_rec_lo[0]*D_surf_C[4]+D_surf_C[1]*df_rec_lo[2]+df_rec_lo[1]*D_surf_C[2]); 
  surft1_lo[5] = 0.3535533905932737*(D_surf_C[2]*df_rec_lo[7]+df_rec_lo[2]*D_surf_C[7]+D_surf_C[4]*df_rec_lo[6]+df_rec_lo[4]*D_surf_C[6]+D_surf_C[0]*df_rec_lo[5]+df_rec_lo[0]*D_surf_C[5]+D_surf_C[1]*df_rec_lo[3]+df_rec_lo[1]*D_surf_C[3]); 
  surft1_lo[6] = 0.3535533905932737*(D_surf_C[1]*df_rec_lo[7]+df_rec_lo[1]*D_surf_C[7]+D_surf_C[0]*df_rec_lo[6]+df_rec_lo[0]*D_surf_C[6]+D_surf_C[4]*df_rec_lo[5]+df_rec_lo[4]*D_surf_C[5]+D_surf_C[2]*df_rec_lo[3]+df_rec_lo[2]*D_surf_C[3]); 
  surft1_lo[7] = 0.3535533905932737*(D_surf_C[0]*df_rec_lo[7]+df_rec_lo[0]*D_surf_C[7]+D_surf_C[1]*df_rec_lo[6]+df_rec_lo[1]*D_surf_C[6]+D_surf_C[2]*df_rec_lo[5]+df_rec_lo[2]*D_surf_C[5]+D_surf_C[3]*df_rec_lo[4]+df_rec_lo[3]*D_surf_C[4]); 

  surft2_lo[0] = 0.3535533905932737*(D_surf_C[7]*f_rec_lo[7]+D_surf_C[6]*f_rec_lo[6]+D_surf_C[5]*f_rec_lo[5]+D_surf_C[4]*f_rec_lo[4]+D_surf_C[3]*f_rec_lo[3]+D_surf_C[2]*f_rec_lo[2]+D_surf_C[1]*f_rec_lo[1]+D_surf_C[0]*f_rec_lo[0]); 
  surft2_lo[1] = 0.3535533905932737*(D_surf_C[6]*f_rec_lo[7]+f_rec_lo[6]*D_surf_C[7]+D_surf_C[3]*f_rec_lo[5]+f_rec_lo[3]*D_surf_C[5]+D_surf_C[2]*f_rec_lo[4]+f_rec_lo[2]*D_surf_C[4]+D_surf_C[0]*f_rec_lo[1]+f_rec_lo[0]*D_surf_C[1]); 
  surft2_lo[2] = 0.3535533905932737*(D_surf_C[5]*f_rec_lo[7]+f_rec_lo[5]*D_surf_C[7]+D_surf_C[3]*f_rec_lo[6]+f_rec_lo[3]*D_surf_C[6]+D_surf_C[1]*f_rec_lo[4]+f_rec_lo[1]*D_surf_C[4]+D_surf_C[0]*f_rec_lo[2]+f_rec_lo[0]*D_surf_C[2]); 
  surft2_lo[3] = 0.3535533905932737*(D_surf_C[4]*f_rec_lo[7]+f_rec_lo[4]*D_surf_C[7]+D_surf_C[2]*f_rec_lo[6]+f_rec_lo[2]*D_surf_C[6]+D_surf_C[1]*f_rec_lo[5]+f_rec_lo[1]*D_surf_C[5]+D_surf_C[0]*f_rec_lo[3]+f_rec_lo[0]*D_surf_C[3]); 
  surft2_lo[4] = 0.3535533905932737*(D_surf_C[3]*f_rec_lo[7]+f_rec_lo[3]*D_surf_C[7]+D_surf_C[5]*f_rec_lo[6]+f_rec_lo[5]*D_surf_C[6]+D_surf_C[0]*f_rec_lo[4]+f_rec_lo[0]*D_surf_C[4]+D_surf_C[1]*f_rec_lo[2]+f_rec_lo[1]*D_surf_C[2]); 
  surft2_lo[5] = 0.3535533905932737*(D_surf_C[2]*f_rec_lo[7]+f_rec_lo[2]*D_surf_C[7]+D_surf_C[4]*f_rec_lo[6]+f_rec_lo[4]*D_surf_C[6]+D_surf_C[0]*f_rec_lo[5]+f_rec_lo[0]*D_surf_C[5]+D_surf_C[1]*f_rec_lo[3]+f_rec_lo[1]*D_surf_C[3]); 
  surft2_lo[6] = 0.3535533905932737*(D_surf_C[1]*f_rec_lo[7]+f_rec_lo[1]*D_surf_C[7]+D_surf_C[0]*f_rec_lo[6]+f_rec_lo[0]*D_surf_C[6]+D_surf_C[4]*f_rec_lo[5]+f_rec_lo[4]*D_surf_C[5]+D_surf_C[2]*f_rec_lo[3]+f_rec_lo[2]*D_surf_C[3]); 
  surft2_lo[7] = 0.3535533905932737*(D_surf_C[0]*f_rec_lo[7]+f_rec_lo[0]*D_surf_C[7]+D_surf_C[1]*f_rec_lo[6]+f_rec_lo[1]*D_surf_C[6]+D_surf_C[2]*f_rec_lo[5]+f_rec_lo[2]*D_surf_C[5]+D_surf_C[3]*f_rec_lo[4]+f_rec_lo[3]*D_surf_C[4]); 
  surft2_up[0] = 0.4330127018922193*(f_rec_up[7]*DC[15]+f_rec_up[6]*DC[14]+f_rec_up[5]*DC[12]+f_rec_up[4]*DC[11]+f_rec_up[3]*DC[9]+f_rec_up[2]*DC[7]+f_rec_up[1]*DC[5]+f_rec_up[0]*DC[2])+0.25*(f_rec_up[7]*DC[13]+f_rec_up[6]*DC[10]+f_rec_up[5]*DC[8]+f_rec_up[4]*DC[6]+f_rec_up[3]*DC[4]+f_rec_up[2]*DC[3]+DC[1]*f_rec_up[1]+DC[0]*f_rec_up[0]); 
  surft2_up[1] = 0.4330127018922193*(f_rec_up[6]*DC[15]+f_rec_up[7]*DC[14]+f_rec_up[3]*DC[12]+f_rec_up[2]*DC[11]+f_rec_up[5]*DC[9]+f_rec_up[4]*DC[7]+f_rec_up[0]*DC[5]+f_rec_up[1]*DC[2])+0.25*(f_rec_up[6]*DC[13]+f_rec_up[7]*DC[10]+f_rec_up[3]*DC[8]+f_rec_up[2]*DC[6]+DC[4]*f_rec_up[5]+DC[3]*f_rec_up[4]+DC[0]*f_rec_up[1]+f_rec_up[0]*DC[1]); 
  surft2_up[2] = 0.4330127018922193*(f_rec_up[5]*DC[15]+f_rec_up[3]*DC[14]+f_rec_up[7]*DC[12]+f_rec_up[1]*DC[11]+f_rec_up[6]*DC[9]+f_rec_up[0]*DC[7]+f_rec_up[4]*DC[5]+DC[2]*f_rec_up[2])+0.25*(f_rec_up[5]*DC[13]+f_rec_up[3]*DC[10]+f_rec_up[7]*DC[8]+DC[4]*f_rec_up[6]+f_rec_up[1]*DC[6]+DC[1]*f_rec_up[4]+f_rec_up[0]*DC[3]+DC[0]*f_rec_up[2]); 
  surft2_up[3] = 0.4330127018922193*(f_rec_up[4]*DC[15]+f_rec_up[2]*DC[14]+f_rec_up[1]*DC[12]+f_rec_up[7]*DC[11]+f_rec_up[0]*DC[9]+f_rec_up[6]*DC[7]+DC[5]*f_rec_up[5]+DC[2]*f_rec_up[3])+0.25*(f_rec_up[4]*DC[13]+f_rec_up[2]*DC[10]+f_rec_up[1]*DC[8]+DC[6]*f_rec_up[7]+DC[3]*f_rec_up[6]+DC[1]*f_rec_up[5]+f_rec_up[0]*DC[4]+DC[0]*f_rec_up[3]); 
  surft2_up[4] = 0.4330127018922193*(f_rec_up[3]*DC[15]+f_rec_up[5]*DC[14]+f_rec_up[6]*DC[12]+f_rec_up[0]*DC[11]+f_rec_up[7]*DC[9]+f_rec_up[1]*DC[7]+f_rec_up[2]*DC[5]+DC[2]*f_rec_up[4])+0.25*(f_rec_up[3]*DC[13]+f_rec_up[5]*DC[10]+f_rec_up[6]*DC[8]+DC[4]*f_rec_up[7]+f_rec_up[0]*DC[6]+DC[0]*f_rec_up[4]+f_rec_up[1]*DC[3]+DC[1]*f_rec_up[2]); 
  surft2_up[5] = 0.4330127018922193*(f_rec_up[2]*DC[15]+f_rec_up[4]*DC[14]+f_rec_up[0]*DC[12]+f_rec_up[6]*DC[11]+f_rec_up[1]*DC[9]+DC[7]*f_rec_up[7]+DC[2]*f_rec_up[5]+f_rec_up[3]*DC[5])+0.25*(f_rec_up[2]*DC[13]+f_rec_up[4]*DC[10]+f_rec_up[0]*DC[8]+DC[3]*f_rec_up[7]+DC[6]*f_rec_up[6]+DC[0]*f_rec_up[5]+f_rec_up[1]*DC[4]+DC[1]*f_rec_up[3]); 
  surft2_up[6] = 0.4330127018922193*(f_rec_up[1]*DC[15]+f_rec_up[0]*DC[14]+f_rec_up[4]*DC[12]+f_rec_up[5]*DC[11]+f_rec_up[2]*DC[9]+DC[5]*f_rec_up[7]+f_rec_up[3]*DC[7]+DC[2]*f_rec_up[6])+0.25*(f_rec_up[1]*DC[13]+f_rec_up[0]*DC[10]+f_rec_up[4]*DC[8]+DC[1]*f_rec_up[7]+DC[0]*f_rec_up[6]+f_rec_up[5]*DC[6]+f_rec_up[2]*DC[4]+DC[3]*f_rec_up[3]); 
  surft2_up[7] = 0.4330127018922193*(f_rec_up[0]*DC[15]+f_rec_up[1]*DC[14]+f_rec_up[2]*DC[12]+f_rec_up[3]*DC[11]+f_rec_up[4]*DC[9]+DC[2]*f_rec_up[7]+f_rec_up[5]*DC[7]+DC[5]*f_rec_up[6])+0.25*(f_rec_up[0]*DC[13]+f_rec_up[1]*DC[10]+f_rec_up[2]*DC[8]+DC[0]*f_rec_up[7]+DC[1]*f_rec_up[6]+f_rec_up[3]*DC[6]+DC[3]*f_rec_up[5]+DC[4]*f_rec_up[4]); 

  vol[2] = 0.75*(fC[13]*DC[15]+fC[10]*DC[14]+fC[8]*DC[12]+fC[6]*DC[11]+fC[4]*DC[9]+fC[3]*DC[7]+fC[1]*DC[5]+fC[0]*DC[2]); 
  vol[5] = 0.75*(fC[10]*DC[15]+fC[13]*DC[14]+fC[4]*DC[12]+fC[3]*DC[11]+fC[8]*DC[9]+fC[6]*DC[7]+fC[0]*DC[5]+fC[1]*DC[2]); 
  vol[7] = 0.75*(fC[8]*DC[15]+fC[4]*DC[14]+DC[12]*fC[13]+fC[1]*DC[11]+DC[9]*fC[10]+fC[0]*DC[7]+DC[5]*fC[6]+DC[2]*fC[3]); 
  vol[9] = 0.75*(fC[6]*DC[15]+fC[3]*DC[14]+DC[11]*fC[13]+fC[1]*DC[12]+DC[7]*fC[10]+fC[0]*DC[9]+DC[5]*fC[8]+DC[2]*fC[4]); 
  vol[11] = 0.75*(fC[4]*DC[15]+fC[8]*DC[14]+DC[9]*fC[13]+fC[10]*DC[12]+fC[0]*DC[11]+fC[1]*DC[7]+DC[2]*fC[6]+fC[3]*DC[5]); 
  vol[12] = 0.75*(fC[3]*DC[15]+fC[6]*DC[14]+DC[7]*fC[13]+fC[0]*DC[12]+fC[10]*DC[11]+fC[1]*DC[9]+DC[2]*fC[8]+fC[4]*DC[5]); 
  vol[14] = 0.75*(fC[1]*DC[15]+fC[0]*DC[14]+DC[5]*fC[13]+fC[6]*DC[12]+fC[8]*DC[11]+DC[2]*fC[10]+fC[3]*DC[9]+fC[4]*DC[7]); 
  vol[15] = 0.75*(fC[0]*DC[15]+fC[1]*DC[14]+DC[2]*fC[13]+fC[3]*DC[12]+fC[4]*DC[11]+DC[5]*fC[10]+fC[6]*DC[9]+DC[7]*fC[8]); 

  out[0] += (0.5*vol[0]-0.35355339059327373*surft1_lo[0])*dv_inv_sq; 
  out[1] += (0.5*vol[1]-0.35355339059327373*surft1_lo[1])*dv_inv_sq; 
  out[2] += (0.5*vol[2]-0.6123724356957945*surft2_up[0]+0.6123724356957945*(surft2_lo[0]+surft1_lo[0]))*dv_inv_sq; 
  out[3] += (0.5*vol[3]-0.35355339059327373*surft1_lo[2])*dv_inv_sq; 
  out[4] += (0.5*vol[4]-0.35355339059327373*surft1_lo[3])*dv_inv_sq; 
  out[5] += (0.5*vol[5]-0.6123724356957945*surft2_up[1]+0.6123724356957945*(surft2_lo[1]+surft1_lo[1]))*dv_inv_sq; 
  out[6] += (0.5*vol[6]-0.35355339059327373*surft1_lo[4])*dv_inv_sq; 
  out[7] += (0.5*vol[7]-0.6123724356957945*surft2_up[2]+0.6123724356957945*(surft2_lo[2]+surft1_lo[2]))*dv_inv_sq; 
  out[8] += (0.5*vol[8]-0.35355339059327373*surft1_lo[5])*dv_inv_sq; 
  out[9] += (0.5*vol[9]-0.6123724356957945*surft2_up[3]+0.6123724356957945*(surft2_lo[3]+surft1_lo[3]))*dv_inv_sq; 
  out[10] += (0.5*vol[10]-0.35355339059327373*surft1_lo[6])*dv_inv_sq; 
  out[11] += (0.5*vol[11]-0.6123724356957945*surft2_up[4]+0.6123724356957945*(surft2_lo[4]+surft1_lo[4]))*dv_inv_sq; 
  out[12] += (0.5*vol[12]-0.6123724356957945*surft2_up[5]+0.6123724356957945*(surft2_lo[5]+surft1_lo[5]))*dv_inv_sq; 
  out[13] += (0.5*vol[13]-0.35355339059327373*surft1_lo[7])*dv_inv_sq; 
  out[14] += (0.5*vol[14]-0.6123724356957945*surft2_up[6]+0.6123724356957945*(surft2_lo[6]+surft1_lo[6]))*dv_inv_sq; 
  out[15] += (0.5*vol[15]-0.6123724356957945*surft2_up[7]+0.6123724356957945*(surft2_lo[7]+surft1_lo[7]))*dv_inv_sq; 

  double cflFreq = fabs(D_surf_C[0]); 

  return 1.4142135623730947*dv_inv_sq*cflFreq; 
} 
