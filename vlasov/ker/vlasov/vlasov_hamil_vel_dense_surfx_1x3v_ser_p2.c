#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfx_1x3v_ser_p2(const double *w, const double *dxv,
  const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel,
  const double *poisson_tensor_conf, const double *hamil,
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[1]; 
  double wv = w[1]; 

  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  const double *jacob_cx_l = &jacob_pos_l[0]; 
  const double *jacob_cx_c = &jacob_pos_c[0]; 
  const double *jacob_cx_r = &jacob_pos_r[0]; 
  const double jacob_cx_l_inv = 1.0/jacob_cx_l[0]; 
  const double jacob_cx_c_inv = 1.0/jacob_cx_c[0]; 
  const double jacob_cx_r_inv = 1.0/jacob_cx_r[0]; 
  double fUpwind_r[20] = {0.0}; 
  double fUpwind_l[20] = {0.0}; 
  double Ghat_r[20] = {0.0}; 
  double Ghat_l[20] = {0.0}; 
  if (wv>0) { 

  fUpwind_r[0] = (1.5811388300841895*fc[11]+1.224744871391589*fc[1]+0.7071067811865475*fc[0])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[1] = (1.5811388300841898*fc[19]+1.224744871391589*fc[5]+0.7071067811865475*fc[2])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[2] = (1.5811388300841898*fc[21]+1.224744871391589*fc[6]+0.7071067811865475*fc[3])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[3] = (1.5811388300841898*fc[25]+1.224744871391589*fc[8]+0.7071067811865475*fc[4])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[4] = (1.5811388300841895*fc[32]+1.224744871391589*fc[15]+0.7071067811865475*fc[7])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[5] = (1.5811388300841895*fc[35]+1.224744871391589*fc[16]+0.7071067811865475*fc[9])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[6] = (1.5811388300841895*fc[37]+1.224744871391589*fc[17]+0.7071067811865475*fc[10])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[7] = (1.224744871391589*fc[20]+0.7071067811865475*fc[12])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[8] = (1.224744871391589*fc[23]+0.7071067811865475*fc[13])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[9] = (1.224744871391589*fc[28]+0.7071067811865475*fc[14])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[10] = (1.5811388300841898*fc[44]+1.224744871391589*fc[31]+0.7071067811865475*fc[18])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[11] = (1.224744871391589*fc[33]+0.7071067811865475*fc[22])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[12] = (1.224744871391589*fc[34]+0.7071067811865475*fc[24])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[13] = (1.224744871391589*fc[36]+0.7071067811865475*fc[26])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[14] = (1.224744871391589*fc[39]+0.7071067811865475*fc[27])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[15] = (1.224744871391589*fc[41]+0.7071067811865475*fc[29])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[16] = (1.224744871391589*fc[42]+0.7071067811865475*fc[30])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[17] = (1.224744871391589*fc[45]+0.7071067811865475*fc[38])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[18] = (1.224744871391589*fc[46]+0.7071067811865475*fc[40])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[19] = (1.224744871391589*fc[47]+0.7071067811865475*fc[43])*jacob_cx_c_inv*jacob_vx_inv; 

  fUpwind_l[0] = (1.5811388300841895*fl[11]+1.224744871391589*fl[1]+0.7071067811865475*fl[0])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[1] = (1.5811388300841898*fl[19]+1.224744871391589*fl[5]+0.7071067811865475*fl[2])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[2] = (1.5811388300841898*fl[21]+1.224744871391589*fl[6]+0.7071067811865475*fl[3])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[3] = (1.5811388300841898*fl[25]+1.224744871391589*fl[8]+0.7071067811865475*fl[4])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[4] = (1.5811388300841895*fl[32]+1.224744871391589*fl[15]+0.7071067811865475*fl[7])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[5] = (1.5811388300841895*fl[35]+1.224744871391589*fl[16]+0.7071067811865475*fl[9])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[6] = (1.5811388300841895*fl[37]+1.224744871391589*fl[17]+0.7071067811865475*fl[10])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[7] = (1.224744871391589*fl[20]+0.7071067811865475*fl[12])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[8] = (1.224744871391589*fl[23]+0.7071067811865475*fl[13])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[9] = (1.224744871391589*fl[28]+0.7071067811865475*fl[14])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[10] = (1.5811388300841898*fl[44]+1.224744871391589*fl[31]+0.7071067811865475*fl[18])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[11] = (1.224744871391589*fl[33]+0.7071067811865475*fl[22])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[12] = (1.224744871391589*fl[34]+0.7071067811865475*fl[24])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[13] = (1.224744871391589*fl[36]+0.7071067811865475*fl[26])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[14] = (1.224744871391589*fl[39]+0.7071067811865475*fl[27])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[15] = (1.224744871391589*fl[41]+0.7071067811865475*fl[29])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[16] = (1.224744871391589*fl[42]+0.7071067811865475*fl[30])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[17] = (1.224744871391589*fl[45]+0.7071067811865475*fl[38])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[18] = (1.224744871391589*fl[46]+0.7071067811865475*fl[40])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[19] = (1.224744871391589*fl[47]+0.7071067811865475*fl[43])*jacob_cx_l_inv*jacob_vx_inv; 

  } else { 

  fUpwind_r[0] = (1.5811388300841895*fr[11]-1.224744871391589*fr[1]+0.7071067811865475*fr[0])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[1] = (1.5811388300841898*fr[19]-1.224744871391589*fr[5]+0.7071067811865475*fr[2])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[2] = (1.5811388300841898*fr[21]-1.224744871391589*fr[6]+0.7071067811865475*fr[3])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[3] = (1.5811388300841898*fr[25]-1.224744871391589*fr[8]+0.7071067811865475*fr[4])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[4] = (1.5811388300841895*fr[32]-1.224744871391589*fr[15]+0.7071067811865475*fr[7])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[5] = (1.5811388300841895*fr[35]-1.224744871391589*fr[16]+0.7071067811865475*fr[9])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[6] = (1.5811388300841895*fr[37]-1.224744871391589*fr[17]+0.7071067811865475*fr[10])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[7] = (0.7071067811865475*fr[12]-1.224744871391589*fr[20])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[8] = (0.7071067811865475*fr[13]-1.224744871391589*fr[23])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[9] = (0.7071067811865475*fr[14]-1.224744871391589*fr[28])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[10] = (1.5811388300841898*fr[44]-1.224744871391589*fr[31]+0.7071067811865475*fr[18])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[11] = (0.7071067811865475*fr[22]-1.224744871391589*fr[33])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[12] = (0.7071067811865475*fr[24]-1.224744871391589*fr[34])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[13] = (0.7071067811865475*fr[26]-1.224744871391589*fr[36])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[14] = (0.7071067811865475*fr[27]-1.224744871391589*fr[39])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[15] = (0.7071067811865475*fr[29]-1.224744871391589*fr[41])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[16] = (0.7071067811865475*fr[30]-1.224744871391589*fr[42])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[17] = (0.7071067811865475*fr[38]-1.224744871391589*fr[45])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[18] = (0.7071067811865475*fr[40]-1.224744871391589*fr[46])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[19] = (0.7071067811865475*fr[43]-1.224744871391589*fr[47])*jacob_cx_r_inv*jacob_vx_inv; 

  fUpwind_l[0] = (1.5811388300841895*fc[11]-1.224744871391589*fc[1]+0.7071067811865475*fc[0])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[1] = (1.5811388300841898*fc[19]-1.224744871391589*fc[5]+0.7071067811865475*fc[2])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[2] = (1.5811388300841898*fc[21]-1.224744871391589*fc[6]+0.7071067811865475*fc[3])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[3] = (1.5811388300841898*fc[25]-1.224744871391589*fc[8]+0.7071067811865475*fc[4])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[4] = (1.5811388300841895*fc[32]-1.224744871391589*fc[15]+0.7071067811865475*fc[7])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[5] = (1.5811388300841895*fc[35]-1.224744871391589*fc[16]+0.7071067811865475*fc[9])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[6] = (1.5811388300841895*fc[37]-1.224744871391589*fc[17]+0.7071067811865475*fc[10])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[7] = (0.7071067811865475*fc[12]-1.224744871391589*fc[20])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[8] = (0.7071067811865475*fc[13]-1.224744871391589*fc[23])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[9] = (0.7071067811865475*fc[14]-1.224744871391589*fc[28])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[10] = (1.5811388300841898*fc[44]-1.224744871391589*fc[31]+0.7071067811865475*fc[18])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[11] = (0.7071067811865475*fc[22]-1.224744871391589*fc[33])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[12] = (0.7071067811865475*fc[24]-1.224744871391589*fc[34])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[13] = (0.7071067811865475*fc[26]-1.224744871391589*fc[36])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[14] = (0.7071067811865475*fc[27]-1.224744871391589*fc[39])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[15] = (0.7071067811865475*fc[29]-1.224744871391589*fc[41])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[16] = (0.7071067811865475*fc[30]-1.224744871391589*fc[42])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[17] = (0.7071067811865475*fc[38]-1.224744871391589*fc[45])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[18] = (0.7071067811865475*fc[40]-1.224744871391589*fc[46])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[19] = (0.7071067811865475*fc[43]-1.224744871391589*fc[47])*jacob_cx_c_inv*jacob_vx_inv; 

  } 
  Ghat_l[0] = 0.6123724356957944*(fUpwind_l[16]*hamil[19]+fUpwind_l[14]*hamil[18])+1.369306393762915*fUpwind_l[10]*hamil[17]+0.6123724356957944*fUpwind_l[9]*hamil[15]+1.369306393762915*fUpwind_l[5]*hamil[13]+0.6123724356957944*fUpwind_l[8]*hamil[12]+1.369306393762915*fUpwind_l[4]*hamil[11]+0.6123724356957944*fUpwind_l[6]*hamil[10]+1.369306393762915*fUpwind_l[1]*hamil[7]+0.6123724356957944*(fUpwind_l[3]*hamil[5]+fUpwind_l[2]*hamil[4]+fUpwind_l[0]*hamil[1]); 
  Ghat_l[1] = 0.6123724356957944*(fUpwind_l[19]*hamil[19]+fUpwind_l[18]*hamil[18])+(1.224744871391589*fUpwind_l[17]+1.369306393762915*fUpwind_l[6])*hamil[17]+0.6123724356957944*fUpwind_l[15]*hamil[15]+(1.224744871391589*fUpwind_l[13]+1.369306393762915*fUpwind_l[3])*hamil[13]+0.6123724356957944*fUpwind_l[12]*hamil[12]+(1.224744871391589*fUpwind_l[11]+1.369306393762915*fUpwind_l[2])*hamil[11]+0.6123724356957944*fUpwind_l[10]*hamil[10]+(1.224744871391589*fUpwind_l[7]+1.369306393762915*fUpwind_l[0])*hamil[7]+0.6123724356957944*(fUpwind_l[5]*hamil[5]+fUpwind_l[4]*hamil[4]+fUpwind_l[1]*hamil[1]); 
  Ghat_l[2] = 0.6123724356957944*fUpwind_l[9]*hamil[19]+0.5477225575051661*fUpwind_l[6]*hamil[18]+hamil[17]*(1.224744871391589*fUpwind_l[18]+1.369306393762915*fUpwind_l[5])+0.6123724356957944*hamil[15]*fUpwind_l[16]+0.5477225575051661*hamil[10]*fUpwind_l[14]+1.369306393762915*fUpwind_l[10]*hamil[13]+0.5477225575051661*fUpwind_l[2]*hamil[12]+hamil[11]*(1.224744871391589*fUpwind_l[12]+1.369306393762915*fUpwind_l[1])+0.6123724356957944*fUpwind_l[3]*hamil[10]+0.5477225575051661*hamil[4]*fUpwind_l[8]+1.369306393762915*fUpwind_l[4]*hamil[7]+0.6123724356957944*(hamil[5]*fUpwind_l[6]+fUpwind_l[0]*hamil[4]+hamil[1]*fUpwind_l[2]); 
  Ghat_l[3] = 0.5477225575051661*fUpwind_l[6]*hamil[19]+1.224744871391589*hamil[17]*fUpwind_l[19]+0.6123724356957944*fUpwind_l[8]*hamil[18]+1.369306393762915*fUpwind_l[4]*hamil[17]+0.5477225575051661*(hamil[10]*fUpwind_l[16]+fUpwind_l[3]*hamil[15])+1.224744871391589*hamil[13]*fUpwind_l[15]+0.6123724356957944*hamil[12]*fUpwind_l[14]+1.369306393762915*(fUpwind_l[1]*hamil[13]+fUpwind_l[10]*hamil[11])+0.6123724356957944*fUpwind_l[2]*hamil[10]+0.5477225575051661*hamil[5]*fUpwind_l[9]+1.369306393762915*fUpwind_l[5]*hamil[7]+0.6123724356957944*(hamil[4]*fUpwind_l[6]+fUpwind_l[0]*hamil[5]+hamil[1]*fUpwind_l[3]); 
  Ghat_l[4] = 0.6123724356957944*(fUpwind_l[15]*hamil[19]+hamil[15]*fUpwind_l[19])+0.5477225575051661*(fUpwind_l[10]*hamil[18]+hamil[10]*fUpwind_l[18])+(1.224744871391589*(fUpwind_l[14]+fUpwind_l[13])+1.369306393762915*fUpwind_l[3])*hamil[17]+hamil[13]*(1.224744871391589*fUpwind_l[17]+1.369306393762915*fUpwind_l[6])+0.5477225575051661*(fUpwind_l[4]*hamil[12]+hamil[4]*fUpwind_l[12])+(1.224744871391589*(fUpwind_l[8]+fUpwind_l[7])+1.369306393762915*fUpwind_l[0])*hamil[11]+1.224744871391589*hamil[7]*fUpwind_l[11]+0.6123724356957944*(fUpwind_l[5]*hamil[10]+hamil[5]*fUpwind_l[10])+1.369306393762915*fUpwind_l[2]*hamil[7]+0.6123724356957944*(fUpwind_l[1]*hamil[4]+hamil[1]*fUpwind_l[4]); 
  Ghat_l[5] = 0.5477225575051661*(fUpwind_l[10]*hamil[19]+hamil[10]*fUpwind_l[19])+0.6123724356957944*(fUpwind_l[12]*hamil[18]+hamil[12]*fUpwind_l[18])+(1.224744871391589*(fUpwind_l[16]+fUpwind_l[11])+1.369306393762915*fUpwind_l[2])*hamil[17]+1.224744871391589*hamil[11]*fUpwind_l[17]+0.5477225575051661*(fUpwind_l[5]*hamil[15]+hamil[5]*fUpwind_l[15])+(1.224744871391589*(fUpwind_l[9]+fUpwind_l[7])+1.369306393762915*fUpwind_l[0])*hamil[13]+1.224744871391589*hamil[7]*fUpwind_l[13]+1.369306393762915*fUpwind_l[6]*hamil[11]+0.6123724356957944*(fUpwind_l[4]*hamil[10]+hamil[4]*fUpwind_l[10])+1.369306393762915*fUpwind_l[3]*hamil[7]+0.6123724356957944*(fUpwind_l[1]*hamil[5]+hamil[1]*fUpwind_l[5]); 
  Ghat_l[6] = (0.4898979485566356*fUpwind_l[14]+0.5477225575051661*fUpwind_l[3])*hamil[19]+1.224744871391589*hamil[13]*fUpwind_l[19]+(0.4898979485566356*fUpwind_l[16]+0.5477225575051661*fUpwind_l[2])*hamil[18]+1.224744871391589*hamil[11]*fUpwind_l[18]+(1.224744871391589*(fUpwind_l[15]+fUpwind_l[12])+1.369306393762915*fUpwind_l[1])*hamil[17]+0.5477225575051661*(hamil[5]*fUpwind_l[16]+fUpwind_l[6]*hamil[15]+hamil[4]*fUpwind_l[14])+1.369306393762915*fUpwind_l[4]*hamil[13]+0.5477225575051661*fUpwind_l[6]*hamil[12]+1.369306393762915*fUpwind_l[5]*hamil[11]+(0.5477225575051661*(fUpwind_l[9]+fUpwind_l[8])+0.6123724356957944*fUpwind_l[0])*hamil[10]+1.369306393762915*hamil[7]*fUpwind_l[10]+0.6123724356957944*(hamil[1]*fUpwind_l[6]+fUpwind_l[2]*hamil[5]+fUpwind_l[3]*hamil[4]); 
  Ghat_l[7] = 1.224744871391589*fUpwind_l[10]*hamil[17]+0.6123724356957944*hamil[10]*fUpwind_l[17]+1.224744871391589*fUpwind_l[5]*hamil[13]+0.6123724356957944*hamil[5]*fUpwind_l[13]+1.224744871391589*fUpwind_l[4]*hamil[11]+0.6123724356957944*hamil[4]*fUpwind_l[11]+1.224744871391589*fUpwind_l[1]*hamil[7]+0.6123724356957944*hamil[1]*fUpwind_l[7]; 
  Ghat_l[8] = 0.5477225575051661*fUpwind_l[16]*hamil[19]+(0.39123039821797573*fUpwind_l[14]+0.6123724356957944*fUpwind_l[3])*hamil[18]+1.369306393762915*hamil[13]*fUpwind_l[18]+1.224744871391589*fUpwind_l[10]*hamil[17]+0.6123724356957944*hamil[5]*fUpwind_l[14]+(0.39123039821797573*fUpwind_l[8]+0.6123724356957944*fUpwind_l[0])*hamil[12]+1.369306393762915*hamil[7]*fUpwind_l[12]+1.224744871391589*fUpwind_l[4]*hamil[11]+0.5477225575051661*fUpwind_l[6]*hamil[10]+0.6123724356957944*hamil[1]*fUpwind_l[8]+0.5477225575051661*fUpwind_l[2]*hamil[4]; 
  Ghat_l[9] = (0.39123039821797573*fUpwind_l[16]+0.6123724356957944*fUpwind_l[2])*hamil[19]+1.369306393762915*hamil[11]*fUpwind_l[19]+0.5477225575051661*fUpwind_l[14]*hamil[18]+1.224744871391589*fUpwind_l[10]*hamil[17]+0.6123724356957944*hamil[4]*fUpwind_l[16]+(0.39123039821797573*fUpwind_l[9]+0.6123724356957944*fUpwind_l[0])*hamil[15]+1.369306393762915*hamil[7]*fUpwind_l[15]+1.224744871391589*fUpwind_l[5]*hamil[13]+0.5477225575051661*fUpwind_l[6]*hamil[10]+0.6123724356957944*hamil[1]*fUpwind_l[9]+0.5477225575051661*fUpwind_l[3]*hamil[5]; 
  Ghat_l[10] = (0.4898979485566357*fUpwind_l[18]+0.5477225575051661*fUpwind_l[5])*hamil[19]+0.4898979485566357*hamil[18]*fUpwind_l[19]+0.5477225575051661*(hamil[5]*fUpwind_l[19]+fUpwind_l[4]*hamil[18]+hamil[4]*fUpwind_l[18])+(1.224744871391589*(fUpwind_l[9]+fUpwind_l[8]+fUpwind_l[7])+1.369306393762915*fUpwind_l[0])*hamil[17]+1.224744871391589*(hamil[7]*fUpwind_l[17]+hamil[13]*fUpwind_l[16])+0.5477225575051661*(fUpwind_l[10]*hamil[15]+hamil[10]*fUpwind_l[15])+1.224744871391589*hamil[11]*fUpwind_l[14]+(1.224744871391589*fUpwind_l[11]+1.369306393762915*fUpwind_l[2])*hamil[13]+1.224744871391589*hamil[11]*fUpwind_l[13]+0.5477225575051661*(fUpwind_l[10]*hamil[12]+hamil[10]*fUpwind_l[12])+1.369306393762915*fUpwind_l[3]*hamil[11]+0.6123724356957944*(fUpwind_l[1]*hamil[10]+hamil[1]*fUpwind_l[10])+1.369306393762915*fUpwind_l[6]*hamil[7]+0.6123724356957944*(fUpwind_l[4]*hamil[5]+hamil[4]*fUpwind_l[5]); 
  Ghat_l[11] = 0.5477225575051661*fUpwind_l[17]*hamil[18]+hamil[17]*(1.0954451150103324*fUpwind_l[18]+1.224744871391589*fUpwind_l[5])+0.6123724356957944*hamil[5]*fUpwind_l[17]+1.224744871391589*fUpwind_l[10]*hamil[13]+0.6123724356957944*hamil[10]*fUpwind_l[13]+0.5477225575051661*fUpwind_l[11]*hamil[12]+hamil[11]*(1.0954451150103324*fUpwind_l[12]+1.224744871391589*fUpwind_l[1])+0.6123724356957944*hamil[1]*fUpwind_l[11]+1.224744871391589*fUpwind_l[4]*hamil[7]+0.6123724356957944*hamil[4]*fUpwind_l[7]; 
  Ghat_l[12] = 0.5477225575051661*fUpwind_l[19]*hamil[19]+0.39123039821797573*fUpwind_l[18]*hamil[18]+0.6123724356957944*(fUpwind_l[5]*hamil[18]+hamil[5]*fUpwind_l[18])+(1.0954451150103324*fUpwind_l[17]+1.224744871391589*fUpwind_l[6])*hamil[17]+1.369306393762915*hamil[13]*fUpwind_l[14]+0.39123039821797573*fUpwind_l[12]*hamil[12]+0.6123724356957944*(fUpwind_l[1]*hamil[12]+hamil[1]*fUpwind_l[12])+(1.0954451150103324*fUpwind_l[11]+1.224744871391589*fUpwind_l[2])*hamil[11]+0.5477225575051661*fUpwind_l[10]*hamil[10]+1.369306393762915*hamil[7]*fUpwind_l[8]+0.5477225575051661*fUpwind_l[4]*hamil[4]; 
  Ghat_l[13] = 0.5477225575051661*fUpwind_l[17]*hamil[19]+hamil[17]*(1.0954451150103324*fUpwind_l[19]+1.224744871391589*fUpwind_l[4])+0.6123724356957944*hamil[4]*fUpwind_l[17]+0.5477225575051661*fUpwind_l[13]*hamil[15]+hamil[13]*(1.0954451150103324*fUpwind_l[15]+1.224744871391589*fUpwind_l[1])+0.6123724356957944*hamil[1]*fUpwind_l[13]+1.224744871391589*fUpwind_l[10]*hamil[11]+0.6123724356957944*hamil[10]*fUpwind_l[11]+1.224744871391589*fUpwind_l[5]*hamil[7]+0.6123724356957944*hamil[5]*fUpwind_l[7]; 
  Ghat_l[14] = 0.4898979485566356*fUpwind_l[6]*hamil[19]+1.0954451150103324*hamil[17]*fUpwind_l[19]+(0.5477225575051661*fUpwind_l[9]+0.39123039821797573*fUpwind_l[8]+0.6123724356957944*fUpwind_l[0])*hamil[18]+1.369306393762915*hamil[7]*fUpwind_l[18]+1.224744871391589*fUpwind_l[4]*hamil[17]+0.4898979485566357*hamil[10]*fUpwind_l[16]+fUpwind_l[14]*(0.5477225575051661*hamil[15]+0.39123039821797573*hamil[12]+0.6123724356957944*hamil[1])+1.369306393762915*fUpwind_l[12]*hamil[13]+0.6123724356957944*fUpwind_l[3]*hamil[12]+1.224744871391589*fUpwind_l[10]*hamil[11]+0.5477225575051661*fUpwind_l[2]*hamil[10]+0.6123724356957944*hamil[5]*fUpwind_l[8]+0.5477225575051661*hamil[4]*fUpwind_l[6]; 
  Ghat_l[15] = 0.39123039821797573*fUpwind_l[19]*hamil[19]+0.6123724356957944*(fUpwind_l[4]*hamil[19]+hamil[4]*fUpwind_l[19])+0.5477225575051661*fUpwind_l[18]*hamil[18]+(1.0954451150103324*fUpwind_l[17]+1.224744871391589*fUpwind_l[6])*hamil[17]+1.369306393762915*hamil[11]*fUpwind_l[16]+0.39123039821797573*fUpwind_l[15]*hamil[15]+0.6123724356957944*(fUpwind_l[1]*hamil[15]+hamil[1]*fUpwind_l[15])+(1.0954451150103324*fUpwind_l[13]+1.224744871391589*fUpwind_l[3])*hamil[13]+0.5477225575051661*fUpwind_l[10]*hamil[10]+1.369306393762915*hamil[7]*fUpwind_l[9]+0.5477225575051661*fUpwind_l[5]*hamil[5]; 
  Ghat_l[16] = (0.39123039821797573*fUpwind_l[9]+0.5477225575051661*fUpwind_l[8]+0.6123724356957944*fUpwind_l[0])*hamil[19]+1.369306393762915*hamil[7]*fUpwind_l[19]+0.4898979485566356*fUpwind_l[6]*hamil[18]+hamil[17]*(1.0954451150103324*fUpwind_l[18]+1.224744871391589*fUpwind_l[5])+(0.39123039821797573*hamil[15]+0.5477225575051661*hamil[12])*fUpwind_l[16]+0.6123724356957944*(hamil[1]*fUpwind_l[16]+fUpwind_l[2]*hamil[15])+1.369306393762915*hamil[11]*fUpwind_l[15]+0.4898979485566357*hamil[10]*fUpwind_l[14]+1.224744871391589*fUpwind_l[10]*hamil[13]+0.5477225575051661*fUpwind_l[3]*hamil[10]+0.6123724356957944*hamil[4]*fUpwind_l[9]+0.5477225575051661*hamil[5]*fUpwind_l[6]; 
  Ghat_l[17] = 0.5477225575051661*fUpwind_l[13]*hamil[19]+1.0954451150103324*hamil[13]*fUpwind_l[19]+0.5477225575051661*fUpwind_l[11]*hamil[18]+1.0954451150103324*hamil[11]*fUpwind_l[18]+(1.0954451150103324*(fUpwind_l[15]+fUpwind_l[12])+1.224744871391589*fUpwind_l[1])*hamil[17]+(0.5477225575051661*(hamil[15]+hamil[12])+0.6123724356957944*hamil[1])*fUpwind_l[17]+1.224744871391589*fUpwind_l[4]*hamil[13]+0.6123724356957944*hamil[4]*fUpwind_l[13]+1.224744871391589*fUpwind_l[5]*hamil[11]+0.6123724356957944*(hamil[5]*fUpwind_l[11]+fUpwind_l[7]*hamil[10])+1.224744871391589*hamil[7]*fUpwind_l[10]; 
  Ghat_l[18] = 0.4898979485566357*(fUpwind_l[10]*hamil[19]+hamil[10]*fUpwind_l[19])+(0.5477225575051661*fUpwind_l[15]+0.39123039821797573*fUpwind_l[12]+0.6123724356957944*fUpwind_l[1])*hamil[18]+(0.5477225575051661*hamil[15]+0.39123039821797573*hamil[12]+0.6123724356957944*hamil[1])*fUpwind_l[18]+(1.0954451150103324*(fUpwind_l[16]+fUpwind_l[11])+1.224744871391589*fUpwind_l[2])*hamil[17]+1.0954451150103324*hamil[11]*fUpwind_l[17]+1.369306393762915*(hamil[7]*fUpwind_l[14]+fUpwind_l[8]*hamil[13])+0.6123724356957944*(fUpwind_l[5]*hamil[12]+hamil[5]*fUpwind_l[12])+1.224744871391589*fUpwind_l[6]*hamil[11]+0.5477225575051661*(fUpwind_l[4]*hamil[10]+hamil[4]*fUpwind_l[10]); 
  Ghat_l[19] = (0.39123039821797573*fUpwind_l[15]+0.5477225575051661*fUpwind_l[12]+0.6123724356957944*fUpwind_l[1])*hamil[19]+(0.39123039821797573*hamil[15]+0.5477225575051661*hamil[12]+0.6123724356957944*hamil[1])*fUpwind_l[19]+0.4898979485566357*(fUpwind_l[10]*hamil[18]+hamil[10]*fUpwind_l[18])+(1.0954451150103324*(fUpwind_l[14]+fUpwind_l[13])+1.224744871391589*fUpwind_l[3])*hamil[17]+1.0954451150103324*hamil[13]*fUpwind_l[17]+1.369306393762915*hamil[7]*fUpwind_l[16]+0.6123724356957944*(fUpwind_l[4]*hamil[15]+hamil[4]*fUpwind_l[15])+1.224744871391589*fUpwind_l[6]*hamil[13]+1.369306393762915*fUpwind_l[9]*hamil[11]+0.5477225575051661*(fUpwind_l[5]*hamil[10]+hamil[5]*fUpwind_l[10]); 

  Ghat_r[0] = 0.6123724356957944*(fUpwind_r[16]*hamil[19]+fUpwind_r[14]*hamil[18])+1.369306393762915*fUpwind_r[10]*hamil[17]+0.6123724356957944*fUpwind_r[9]*hamil[15]+1.369306393762915*fUpwind_r[5]*hamil[13]+0.6123724356957944*fUpwind_r[8]*hamil[12]+1.369306393762915*fUpwind_r[4]*hamil[11]+0.6123724356957944*fUpwind_r[6]*hamil[10]+1.369306393762915*fUpwind_r[1]*hamil[7]+0.6123724356957944*(fUpwind_r[3]*hamil[5]+fUpwind_r[2]*hamil[4]+fUpwind_r[0]*hamil[1]); 
  Ghat_r[1] = 0.6123724356957944*(fUpwind_r[19]*hamil[19]+fUpwind_r[18]*hamil[18])+(1.224744871391589*fUpwind_r[17]+1.369306393762915*fUpwind_r[6])*hamil[17]+0.6123724356957944*fUpwind_r[15]*hamil[15]+(1.224744871391589*fUpwind_r[13]+1.369306393762915*fUpwind_r[3])*hamil[13]+0.6123724356957944*fUpwind_r[12]*hamil[12]+(1.224744871391589*fUpwind_r[11]+1.369306393762915*fUpwind_r[2])*hamil[11]+0.6123724356957944*fUpwind_r[10]*hamil[10]+(1.224744871391589*fUpwind_r[7]+1.369306393762915*fUpwind_r[0])*hamil[7]+0.6123724356957944*(fUpwind_r[5]*hamil[5]+fUpwind_r[4]*hamil[4]+fUpwind_r[1]*hamil[1]); 
  Ghat_r[2] = 0.6123724356957944*fUpwind_r[9]*hamil[19]+0.5477225575051661*fUpwind_r[6]*hamil[18]+hamil[17]*(1.224744871391589*fUpwind_r[18]+1.369306393762915*fUpwind_r[5])+0.6123724356957944*hamil[15]*fUpwind_r[16]+0.5477225575051661*hamil[10]*fUpwind_r[14]+1.369306393762915*fUpwind_r[10]*hamil[13]+0.5477225575051661*fUpwind_r[2]*hamil[12]+hamil[11]*(1.224744871391589*fUpwind_r[12]+1.369306393762915*fUpwind_r[1])+0.6123724356957944*fUpwind_r[3]*hamil[10]+0.5477225575051661*hamil[4]*fUpwind_r[8]+1.369306393762915*fUpwind_r[4]*hamil[7]+0.6123724356957944*(hamil[5]*fUpwind_r[6]+fUpwind_r[0]*hamil[4]+hamil[1]*fUpwind_r[2]); 
  Ghat_r[3] = 0.5477225575051661*fUpwind_r[6]*hamil[19]+1.224744871391589*hamil[17]*fUpwind_r[19]+0.6123724356957944*fUpwind_r[8]*hamil[18]+1.369306393762915*fUpwind_r[4]*hamil[17]+0.5477225575051661*(hamil[10]*fUpwind_r[16]+fUpwind_r[3]*hamil[15])+1.224744871391589*hamil[13]*fUpwind_r[15]+0.6123724356957944*hamil[12]*fUpwind_r[14]+1.369306393762915*(fUpwind_r[1]*hamil[13]+fUpwind_r[10]*hamil[11])+0.6123724356957944*fUpwind_r[2]*hamil[10]+0.5477225575051661*hamil[5]*fUpwind_r[9]+1.369306393762915*fUpwind_r[5]*hamil[7]+0.6123724356957944*(hamil[4]*fUpwind_r[6]+fUpwind_r[0]*hamil[5]+hamil[1]*fUpwind_r[3]); 
  Ghat_r[4] = 0.6123724356957944*(fUpwind_r[15]*hamil[19]+hamil[15]*fUpwind_r[19])+0.5477225575051661*(fUpwind_r[10]*hamil[18]+hamil[10]*fUpwind_r[18])+(1.224744871391589*(fUpwind_r[14]+fUpwind_r[13])+1.369306393762915*fUpwind_r[3])*hamil[17]+hamil[13]*(1.224744871391589*fUpwind_r[17]+1.369306393762915*fUpwind_r[6])+0.5477225575051661*(fUpwind_r[4]*hamil[12]+hamil[4]*fUpwind_r[12])+(1.224744871391589*(fUpwind_r[8]+fUpwind_r[7])+1.369306393762915*fUpwind_r[0])*hamil[11]+1.224744871391589*hamil[7]*fUpwind_r[11]+0.6123724356957944*(fUpwind_r[5]*hamil[10]+hamil[5]*fUpwind_r[10])+1.369306393762915*fUpwind_r[2]*hamil[7]+0.6123724356957944*(fUpwind_r[1]*hamil[4]+hamil[1]*fUpwind_r[4]); 
  Ghat_r[5] = 0.5477225575051661*(fUpwind_r[10]*hamil[19]+hamil[10]*fUpwind_r[19])+0.6123724356957944*(fUpwind_r[12]*hamil[18]+hamil[12]*fUpwind_r[18])+(1.224744871391589*(fUpwind_r[16]+fUpwind_r[11])+1.369306393762915*fUpwind_r[2])*hamil[17]+1.224744871391589*hamil[11]*fUpwind_r[17]+0.5477225575051661*(fUpwind_r[5]*hamil[15]+hamil[5]*fUpwind_r[15])+(1.224744871391589*(fUpwind_r[9]+fUpwind_r[7])+1.369306393762915*fUpwind_r[0])*hamil[13]+1.224744871391589*hamil[7]*fUpwind_r[13]+1.369306393762915*fUpwind_r[6]*hamil[11]+0.6123724356957944*(fUpwind_r[4]*hamil[10]+hamil[4]*fUpwind_r[10])+1.369306393762915*fUpwind_r[3]*hamil[7]+0.6123724356957944*(fUpwind_r[1]*hamil[5]+hamil[1]*fUpwind_r[5]); 
  Ghat_r[6] = (0.4898979485566356*fUpwind_r[14]+0.5477225575051661*fUpwind_r[3])*hamil[19]+1.224744871391589*hamil[13]*fUpwind_r[19]+(0.4898979485566356*fUpwind_r[16]+0.5477225575051661*fUpwind_r[2])*hamil[18]+1.224744871391589*hamil[11]*fUpwind_r[18]+(1.224744871391589*(fUpwind_r[15]+fUpwind_r[12])+1.369306393762915*fUpwind_r[1])*hamil[17]+0.5477225575051661*(hamil[5]*fUpwind_r[16]+fUpwind_r[6]*hamil[15]+hamil[4]*fUpwind_r[14])+1.369306393762915*fUpwind_r[4]*hamil[13]+0.5477225575051661*fUpwind_r[6]*hamil[12]+1.369306393762915*fUpwind_r[5]*hamil[11]+(0.5477225575051661*(fUpwind_r[9]+fUpwind_r[8])+0.6123724356957944*fUpwind_r[0])*hamil[10]+1.369306393762915*hamil[7]*fUpwind_r[10]+0.6123724356957944*(hamil[1]*fUpwind_r[6]+fUpwind_r[2]*hamil[5]+fUpwind_r[3]*hamil[4]); 
  Ghat_r[7] = 1.224744871391589*fUpwind_r[10]*hamil[17]+0.6123724356957944*hamil[10]*fUpwind_r[17]+1.224744871391589*fUpwind_r[5]*hamil[13]+0.6123724356957944*hamil[5]*fUpwind_r[13]+1.224744871391589*fUpwind_r[4]*hamil[11]+0.6123724356957944*hamil[4]*fUpwind_r[11]+1.224744871391589*fUpwind_r[1]*hamil[7]+0.6123724356957944*hamil[1]*fUpwind_r[7]; 
  Ghat_r[8] = 0.5477225575051661*fUpwind_r[16]*hamil[19]+(0.39123039821797573*fUpwind_r[14]+0.6123724356957944*fUpwind_r[3])*hamil[18]+1.369306393762915*hamil[13]*fUpwind_r[18]+1.224744871391589*fUpwind_r[10]*hamil[17]+0.6123724356957944*hamil[5]*fUpwind_r[14]+(0.39123039821797573*fUpwind_r[8]+0.6123724356957944*fUpwind_r[0])*hamil[12]+1.369306393762915*hamil[7]*fUpwind_r[12]+1.224744871391589*fUpwind_r[4]*hamil[11]+0.5477225575051661*fUpwind_r[6]*hamil[10]+0.6123724356957944*hamil[1]*fUpwind_r[8]+0.5477225575051661*fUpwind_r[2]*hamil[4]; 
  Ghat_r[9] = (0.39123039821797573*fUpwind_r[16]+0.6123724356957944*fUpwind_r[2])*hamil[19]+1.369306393762915*hamil[11]*fUpwind_r[19]+0.5477225575051661*fUpwind_r[14]*hamil[18]+1.224744871391589*fUpwind_r[10]*hamil[17]+0.6123724356957944*hamil[4]*fUpwind_r[16]+(0.39123039821797573*fUpwind_r[9]+0.6123724356957944*fUpwind_r[0])*hamil[15]+1.369306393762915*hamil[7]*fUpwind_r[15]+1.224744871391589*fUpwind_r[5]*hamil[13]+0.5477225575051661*fUpwind_r[6]*hamil[10]+0.6123724356957944*hamil[1]*fUpwind_r[9]+0.5477225575051661*fUpwind_r[3]*hamil[5]; 
  Ghat_r[10] = (0.4898979485566357*fUpwind_r[18]+0.5477225575051661*fUpwind_r[5])*hamil[19]+0.4898979485566357*hamil[18]*fUpwind_r[19]+0.5477225575051661*(hamil[5]*fUpwind_r[19]+fUpwind_r[4]*hamil[18]+hamil[4]*fUpwind_r[18])+(1.224744871391589*(fUpwind_r[9]+fUpwind_r[8]+fUpwind_r[7])+1.369306393762915*fUpwind_r[0])*hamil[17]+1.224744871391589*(hamil[7]*fUpwind_r[17]+hamil[13]*fUpwind_r[16])+0.5477225575051661*(fUpwind_r[10]*hamil[15]+hamil[10]*fUpwind_r[15])+1.224744871391589*hamil[11]*fUpwind_r[14]+(1.224744871391589*fUpwind_r[11]+1.369306393762915*fUpwind_r[2])*hamil[13]+1.224744871391589*hamil[11]*fUpwind_r[13]+0.5477225575051661*(fUpwind_r[10]*hamil[12]+hamil[10]*fUpwind_r[12])+1.369306393762915*fUpwind_r[3]*hamil[11]+0.6123724356957944*(fUpwind_r[1]*hamil[10]+hamil[1]*fUpwind_r[10])+1.369306393762915*fUpwind_r[6]*hamil[7]+0.6123724356957944*(fUpwind_r[4]*hamil[5]+hamil[4]*fUpwind_r[5]); 
  Ghat_r[11] = 0.5477225575051661*fUpwind_r[17]*hamil[18]+hamil[17]*(1.0954451150103324*fUpwind_r[18]+1.224744871391589*fUpwind_r[5])+0.6123724356957944*hamil[5]*fUpwind_r[17]+1.224744871391589*fUpwind_r[10]*hamil[13]+0.6123724356957944*hamil[10]*fUpwind_r[13]+0.5477225575051661*fUpwind_r[11]*hamil[12]+hamil[11]*(1.0954451150103324*fUpwind_r[12]+1.224744871391589*fUpwind_r[1])+0.6123724356957944*hamil[1]*fUpwind_r[11]+1.224744871391589*fUpwind_r[4]*hamil[7]+0.6123724356957944*hamil[4]*fUpwind_r[7]; 
  Ghat_r[12] = 0.5477225575051661*fUpwind_r[19]*hamil[19]+0.39123039821797573*fUpwind_r[18]*hamil[18]+0.6123724356957944*(fUpwind_r[5]*hamil[18]+hamil[5]*fUpwind_r[18])+(1.0954451150103324*fUpwind_r[17]+1.224744871391589*fUpwind_r[6])*hamil[17]+1.369306393762915*hamil[13]*fUpwind_r[14]+0.39123039821797573*fUpwind_r[12]*hamil[12]+0.6123724356957944*(fUpwind_r[1]*hamil[12]+hamil[1]*fUpwind_r[12])+(1.0954451150103324*fUpwind_r[11]+1.224744871391589*fUpwind_r[2])*hamil[11]+0.5477225575051661*fUpwind_r[10]*hamil[10]+1.369306393762915*hamil[7]*fUpwind_r[8]+0.5477225575051661*fUpwind_r[4]*hamil[4]; 
  Ghat_r[13] = 0.5477225575051661*fUpwind_r[17]*hamil[19]+hamil[17]*(1.0954451150103324*fUpwind_r[19]+1.224744871391589*fUpwind_r[4])+0.6123724356957944*hamil[4]*fUpwind_r[17]+0.5477225575051661*fUpwind_r[13]*hamil[15]+hamil[13]*(1.0954451150103324*fUpwind_r[15]+1.224744871391589*fUpwind_r[1])+0.6123724356957944*hamil[1]*fUpwind_r[13]+1.224744871391589*fUpwind_r[10]*hamil[11]+0.6123724356957944*hamil[10]*fUpwind_r[11]+1.224744871391589*fUpwind_r[5]*hamil[7]+0.6123724356957944*hamil[5]*fUpwind_r[7]; 
  Ghat_r[14] = 0.4898979485566356*fUpwind_r[6]*hamil[19]+1.0954451150103324*hamil[17]*fUpwind_r[19]+(0.5477225575051661*fUpwind_r[9]+0.39123039821797573*fUpwind_r[8]+0.6123724356957944*fUpwind_r[0])*hamil[18]+1.369306393762915*hamil[7]*fUpwind_r[18]+1.224744871391589*fUpwind_r[4]*hamil[17]+0.4898979485566357*hamil[10]*fUpwind_r[16]+fUpwind_r[14]*(0.5477225575051661*hamil[15]+0.39123039821797573*hamil[12]+0.6123724356957944*hamil[1])+1.369306393762915*fUpwind_r[12]*hamil[13]+0.6123724356957944*fUpwind_r[3]*hamil[12]+1.224744871391589*fUpwind_r[10]*hamil[11]+0.5477225575051661*fUpwind_r[2]*hamil[10]+0.6123724356957944*hamil[5]*fUpwind_r[8]+0.5477225575051661*hamil[4]*fUpwind_r[6]; 
  Ghat_r[15] = 0.39123039821797573*fUpwind_r[19]*hamil[19]+0.6123724356957944*(fUpwind_r[4]*hamil[19]+hamil[4]*fUpwind_r[19])+0.5477225575051661*fUpwind_r[18]*hamil[18]+(1.0954451150103324*fUpwind_r[17]+1.224744871391589*fUpwind_r[6])*hamil[17]+1.369306393762915*hamil[11]*fUpwind_r[16]+0.39123039821797573*fUpwind_r[15]*hamil[15]+0.6123724356957944*(fUpwind_r[1]*hamil[15]+hamil[1]*fUpwind_r[15])+(1.0954451150103324*fUpwind_r[13]+1.224744871391589*fUpwind_r[3])*hamil[13]+0.5477225575051661*fUpwind_r[10]*hamil[10]+1.369306393762915*hamil[7]*fUpwind_r[9]+0.5477225575051661*fUpwind_r[5]*hamil[5]; 
  Ghat_r[16] = (0.39123039821797573*fUpwind_r[9]+0.5477225575051661*fUpwind_r[8]+0.6123724356957944*fUpwind_r[0])*hamil[19]+1.369306393762915*hamil[7]*fUpwind_r[19]+0.4898979485566356*fUpwind_r[6]*hamil[18]+hamil[17]*(1.0954451150103324*fUpwind_r[18]+1.224744871391589*fUpwind_r[5])+(0.39123039821797573*hamil[15]+0.5477225575051661*hamil[12])*fUpwind_r[16]+0.6123724356957944*(hamil[1]*fUpwind_r[16]+fUpwind_r[2]*hamil[15])+1.369306393762915*hamil[11]*fUpwind_r[15]+0.4898979485566357*hamil[10]*fUpwind_r[14]+1.224744871391589*fUpwind_r[10]*hamil[13]+0.5477225575051661*fUpwind_r[3]*hamil[10]+0.6123724356957944*hamil[4]*fUpwind_r[9]+0.5477225575051661*hamil[5]*fUpwind_r[6]; 
  Ghat_r[17] = 0.5477225575051661*fUpwind_r[13]*hamil[19]+1.0954451150103324*hamil[13]*fUpwind_r[19]+0.5477225575051661*fUpwind_r[11]*hamil[18]+1.0954451150103324*hamil[11]*fUpwind_r[18]+(1.0954451150103324*(fUpwind_r[15]+fUpwind_r[12])+1.224744871391589*fUpwind_r[1])*hamil[17]+(0.5477225575051661*(hamil[15]+hamil[12])+0.6123724356957944*hamil[1])*fUpwind_r[17]+1.224744871391589*fUpwind_r[4]*hamil[13]+0.6123724356957944*hamil[4]*fUpwind_r[13]+1.224744871391589*fUpwind_r[5]*hamil[11]+0.6123724356957944*(hamil[5]*fUpwind_r[11]+fUpwind_r[7]*hamil[10])+1.224744871391589*hamil[7]*fUpwind_r[10]; 
  Ghat_r[18] = 0.4898979485566357*(fUpwind_r[10]*hamil[19]+hamil[10]*fUpwind_r[19])+(0.5477225575051661*fUpwind_r[15]+0.39123039821797573*fUpwind_r[12]+0.6123724356957944*fUpwind_r[1])*hamil[18]+(0.5477225575051661*hamil[15]+0.39123039821797573*hamil[12]+0.6123724356957944*hamil[1])*fUpwind_r[18]+(1.0954451150103324*(fUpwind_r[16]+fUpwind_r[11])+1.224744871391589*fUpwind_r[2])*hamil[17]+1.0954451150103324*hamil[11]*fUpwind_r[17]+1.369306393762915*(hamil[7]*fUpwind_r[14]+fUpwind_r[8]*hamil[13])+0.6123724356957944*(fUpwind_r[5]*hamil[12]+hamil[5]*fUpwind_r[12])+1.224744871391589*fUpwind_r[6]*hamil[11]+0.5477225575051661*(fUpwind_r[4]*hamil[10]+hamil[4]*fUpwind_r[10]); 
  Ghat_r[19] = (0.39123039821797573*fUpwind_r[15]+0.5477225575051661*fUpwind_r[12]+0.6123724356957944*fUpwind_r[1])*hamil[19]+(0.39123039821797573*hamil[15]+0.5477225575051661*hamil[12]+0.6123724356957944*hamil[1])*fUpwind_r[19]+0.4898979485566357*(fUpwind_r[10]*hamil[18]+hamil[10]*fUpwind_r[18])+(1.0954451150103324*(fUpwind_r[14]+fUpwind_r[13])+1.224744871391589*fUpwind_r[3])*hamil[17]+1.0954451150103324*hamil[13]*fUpwind_r[17]+1.369306393762915*hamil[7]*fUpwind_r[16]+0.6123724356957944*(fUpwind_r[4]*hamil[15]+hamil[4]*fUpwind_r[15])+1.224744871391589*fUpwind_r[6]*hamil[13]+1.369306393762915*fUpwind_r[9]*hamil[11]+0.5477225575051661*(fUpwind_r[5]*hamil[10]+hamil[5]*fUpwind_r[10]); 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv10*dx10; 
  out[1] += -(1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv10*dx10); 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv10*dx10; 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv10*dx10; 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv10*dx10; 
  out[5] += -(1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv10*dx10); 
  out[6] += -(1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv10*dx10); 
  out[7] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dv10*dx10; 
  out[8] += -(1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv10*dx10); 
  out[9] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dv10*dx10; 
  out[10] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dv10*dx10; 
  out[11] += (1.5811388300841895*Ghat_l[0]-1.5811388300841895*Ghat_r[0])*dv10*dx10; 
  out[12] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dv10*dx10; 
  out[13] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dv10*dx10; 
  out[14] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*dv10*dx10; 
  out[15] += -(1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv10*dx10); 
  out[16] += -(1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv10*dx10); 
  out[17] += -(1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dv10*dx10); 
  out[18] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*dv10*dx10; 
  out[19] += (1.5811388300841898*Ghat_l[1]-1.5811388300841898*Ghat_r[1])*dv10*dx10; 
  out[20] += -(1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dv10*dx10); 
  out[21] += (1.5811388300841898*Ghat_l[2]-1.5811388300841898*Ghat_r[2])*dv10*dx10; 
  out[22] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*dv10*dx10; 
  out[23] += -(1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dv10*dx10); 
  out[24] += (0.7071067811865475*Ghat_l[12]-0.7071067811865475*Ghat_r[12])*dv10*dx10; 
  out[25] += (1.5811388300841898*Ghat_l[3]-1.5811388300841898*Ghat_r[3])*dv10*dx10; 
  out[26] += (0.7071067811865475*Ghat_l[13]-0.7071067811865475*Ghat_r[13])*dv10*dx10; 
  out[27] += (0.7071067811865475*Ghat_l[14]-0.7071067811865475*Ghat_r[14])*dv10*dx10; 
  out[28] += -(1.224744871391589*(Ghat_r[9]+Ghat_l[9])*dv10*dx10); 
  out[29] += (0.7071067811865475*Ghat_l[15]-0.7071067811865475*Ghat_r[15])*dv10*dx10; 
  out[30] += (0.7071067811865475*Ghat_l[16]-0.7071067811865475*Ghat_r[16])*dv10*dx10; 
  out[31] += -(1.224744871391589*(Ghat_r[10]+Ghat_l[10])*dv10*dx10); 
  out[32] += (1.5811388300841895*Ghat_l[4]-1.5811388300841895*Ghat_r[4])*dv10*dx10; 
  out[33] += -(1.224744871391589*(Ghat_r[11]+Ghat_l[11])*dv10*dx10); 
  out[34] += -(1.224744871391589*(Ghat_r[12]+Ghat_l[12])*dv10*dx10); 
  out[35] += (1.5811388300841895*Ghat_l[5]-1.5811388300841895*Ghat_r[5])*dv10*dx10; 
  out[36] += -(1.224744871391589*(Ghat_r[13]+Ghat_l[13])*dv10*dx10); 
  out[37] += (1.5811388300841895*Ghat_l[6]-1.5811388300841895*Ghat_r[6])*dv10*dx10; 
  out[38] += (0.7071067811865475*Ghat_l[17]-0.7071067811865475*Ghat_r[17])*dv10*dx10; 
  out[39] += -(1.224744871391589*(Ghat_r[14]+Ghat_l[14])*dv10*dx10); 
  out[40] += (0.7071067811865475*Ghat_l[18]-0.7071067811865475*Ghat_r[18])*dv10*dx10; 
  out[41] += -(1.224744871391589*(Ghat_r[15]+Ghat_l[15])*dv10*dx10); 
  out[42] += -(1.224744871391589*(Ghat_r[16]+Ghat_l[16])*dv10*dx10); 
  out[43] += (0.7071067811865475*Ghat_l[19]-0.7071067811865475*Ghat_r[19])*dv10*dx10; 
  out[44] += (1.5811388300841898*Ghat_l[10]-1.5811388300841898*Ghat_r[10])*dv10*dx10; 
  out[45] += -(1.224744871391589*(Ghat_r[17]+Ghat_l[17])*dv10*dx10); 
  out[46] += -(1.224744871391589*(Ghat_r[18]+Ghat_l[18])*dv10*dx10); 
  out[47] += -(1.224744871391589*(Ghat_r[19]+Ghat_l[19])*dv10*dx10); 

  return fabs(0.8838834764831842*dv10*dx10*jacob_cx_c_inv*(1.7320508075688772*hamil[1]*jacob_vx_inv));

} 
