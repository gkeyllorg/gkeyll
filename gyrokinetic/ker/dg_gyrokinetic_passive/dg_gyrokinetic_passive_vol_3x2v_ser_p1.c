#include <gkyl_dg_gyrokinetic_passive_kernels.h> 
GKYL_CU_DH double dg_gyrokinetic_passive_vol_3x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *vmapSq,
    const double q_, const double m_, const double *bmag, const double *speeds,
    const double *dualcurlbhatoverB, const double *rtg33inv, const double *bioverJB,
    const double *fin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // vmap: velocity space mapping.
  // vmapSq: velocity space mapping squared.
  // q_,m_: species charge and mass.
  // bmag: magnetic field amplitude.
  // speeds: passive advection speeds.
  // fin: Distribution function.
  // out: output increment.

  double rdx2 = 2.0/dxv[0];
  double rdy2 = 2.0/dxv[1];
  double rdz2 = 2.0/dxv[2];
  double rdvpar2 = 2.0/dxv[3];
  double rdmu2 = 2.0/dxv[4];

  double rdvpar2Sq = rdvpar2*rdvpar2;
  double dvparSq = dxv[3]*dxv[3];

  const double *bioverJB_x = &bioverJB[0]; 
  const double *bioverJB_y = &bioverJB[8]; 
  const double *bioverJB_z = &bioverJB[16]; 

  const double *dualcurlbhatoverB_x = &dualcurlbhatoverB[0]; 
  const double *dualcurlbhatoverB_y = &dualcurlbhatoverB[8]; 
  const double *dualcurlbhatoverB_z = &dualcurlbhatoverB[16]; 

  double alphax[48] = {0.}; 
  alphax[0] = 2.0*speeds[0]*rdx2; 
  alphax[1] = 2.0*speeds[1]*rdx2; 
  alphax[2] = 2.0*speeds[2]*rdx2; 
  alphax[3] = 2.0*speeds[3]*rdx2; 
  alphax[6] = 2.0*speeds[4]*rdx2; 
  alphax[7] = 2.0*speeds[5]*rdx2; 
  alphax[8] = 2.0*speeds[6]*rdx2; 
  alphax[16] = 2.0*speeds[7]*rdx2; 


  out[1] += 0.3061862178478971*(alphax[16]*fin[16]+alphax[8]*fin[8]+alphax[7]*fin[7]+alphax[6]*fin[6]+alphax[3]*fin[3]+alphax[2]*fin[2]+alphax[1]*fin[1]+alphax[0]*fin[0]); 
  out[6] += 0.3061862178478971*(alphax[7]*fin[16]+fin[7]*alphax[16]+alphax[3]*fin[8]+fin[3]*alphax[8]+alphax[1]*fin[6]+fin[1]*alphax[6]+alphax[0]*fin[2]+fin[0]*alphax[2]); 
  out[7] += 0.3061862178478971*(alphax[6]*fin[16]+fin[6]*alphax[16]+alphax[2]*fin[8]+fin[2]*alphax[8]+alphax[1]*fin[7]+fin[1]*alphax[7]+alphax[0]*fin[3]+fin[0]*alphax[3]); 
  out[9] += 0.3061862178478971*(alphax[16]*fin[26]+alphax[8]*fin[19]+alphax[7]*fin[18]+alphax[6]*fin[17]+alphax[3]*fin[11]+alphax[2]*fin[10]+alphax[1]*fin[9]+alphax[0]*fin[4]); 
  out[12] += 0.3061862178478971*(alphax[16]*fin[27]+alphax[8]*fin[22]+alphax[7]*fin[21]+alphax[6]*fin[20]+alphax[3]*fin[14]+alphax[2]*fin[13]+alphax[1]*fin[12]+alphax[0]*fin[5]); 
  out[16] += 0.3061862178478971*(alphax[1]*fin[16]+fin[1]*alphax[16]+alphax[0]*fin[8]+fin[0]*alphax[8]+alphax[6]*fin[7]+fin[6]*alphax[7]+alphax[2]*fin[3]+fin[2]*alphax[3]); 
  out[17] += 0.3061862178478971*(alphax[7]*fin[26]+alphax[3]*fin[19]+alphax[16]*fin[18]+alphax[1]*fin[17]+alphax[8]*fin[11]+alphax[0]*fin[10]+alphax[6]*fin[9]+alphax[2]*fin[4]); 
  out[18] += 0.3061862178478971*(alphax[6]*fin[26]+alphax[2]*fin[19]+alphax[1]*fin[18]+alphax[16]*fin[17]+alphax[0]*fin[11]+alphax[8]*fin[10]+alphax[7]*fin[9]+alphax[3]*fin[4]); 
  out[20] += 0.3061862178478971*(alphax[7]*fin[27]+alphax[3]*fin[22]+alphax[16]*fin[21]+alphax[1]*fin[20]+alphax[8]*fin[14]+alphax[0]*fin[13]+alphax[6]*fin[12]+alphax[2]*fin[5]); 
  out[21] += 0.3061862178478971*(alphax[6]*fin[27]+alphax[2]*fin[22]+alphax[1]*fin[21]+alphax[16]*fin[20]+alphax[0]*fin[14]+alphax[8]*fin[13]+alphax[7]*fin[12]+alphax[3]*fin[5]); 
  out[23] += 0.3061862178478971*(alphax[16]*fin[31]+alphax[8]*fin[30]+alphax[7]*fin[29]+alphax[6]*fin[28]+alphax[3]*fin[25]+alphax[2]*fin[24]+alphax[1]*fin[23]+alphax[0]*fin[15]); 
  out[26] += 0.3061862178478971*(alphax[1]*fin[26]+alphax[0]*fin[19]+alphax[6]*fin[18]+alphax[7]*fin[17]+fin[9]*alphax[16]+alphax[2]*fin[11]+alphax[3]*fin[10]+fin[4]*alphax[8]); 
  out[27] += 0.3061862178478971*(alphax[1]*fin[27]+alphax[0]*fin[22]+alphax[6]*fin[21]+alphax[7]*fin[20]+fin[12]*alphax[16]+alphax[2]*fin[14]+alphax[3]*fin[13]+fin[5]*alphax[8]); 
  out[28] += 0.3061862178478971*(alphax[7]*fin[31]+alphax[3]*fin[30]+alphax[16]*fin[29]+alphax[1]*fin[28]+alphax[8]*fin[25]+alphax[0]*fin[24]+alphax[6]*fin[23]+alphax[2]*fin[15]); 
  out[29] += 0.3061862178478971*(alphax[6]*fin[31]+alphax[2]*fin[30]+alphax[1]*fin[29]+alphax[16]*fin[28]+alphax[0]*fin[25]+alphax[8]*fin[24]+alphax[7]*fin[23]+alphax[3]*fin[15]); 
  out[31] += 0.3061862178478971*(alphax[1]*fin[31]+alphax[0]*fin[30]+alphax[6]*fin[29]+alphax[7]*fin[28]+alphax[2]*fin[25]+alphax[3]*fin[24]+alphax[16]*fin[23]+alphax[8]*fin[15]); 
  out[33] += 0.3061862178478971*(alphax[16]*fin[43]+alphax[8]*fin[39]+alphax[7]*fin[38]+alphax[6]*fin[37]+alphax[3]*fin[35]+alphax[2]*fin[34]+alphax[1]*fin[33]+alphax[0]*fin[32]); 
  out[37] += 0.3061862178478971*(alphax[7]*fin[43]+alphax[3]*fin[39]+alphax[16]*fin[38]+alphax[1]*fin[37]+alphax[8]*fin[35]+alphax[0]*fin[34]+alphax[6]*fin[33]+alphax[2]*fin[32]); 
  out[38] += 0.3061862178478971*(alphax[6]*fin[43]+alphax[2]*fin[39]+alphax[1]*fin[38]+alphax[16]*fin[37]+alphax[0]*fin[35]+alphax[8]*fin[34]+alphax[7]*fin[33]+alphax[3]*fin[32]); 
  out[40] += 0.3061862178478971*(alphax[16]*fin[47]+alphax[8]*fin[46]+alphax[7]*fin[45]+alphax[6]*fin[44]+alphax[3]*fin[42]+alphax[2]*fin[41]+alphax[1]*fin[40]+alphax[0]*fin[36]); 
  out[43] += 0.3061862178478971*(alphax[1]*fin[43]+alphax[0]*fin[39]+alphax[6]*fin[38]+alphax[7]*fin[37]+alphax[2]*fin[35]+alphax[3]*fin[34]+alphax[16]*fin[33]+alphax[8]*fin[32]); 
  out[44] += 0.3061862178478971*(alphax[7]*fin[47]+alphax[3]*fin[46]+alphax[16]*fin[45]+alphax[1]*fin[44]+alphax[8]*fin[42]+alphax[0]*fin[41]+alphax[6]*fin[40]+alphax[2]*fin[36]); 
  out[45] += 0.3061862178478971*(alphax[6]*fin[47]+alphax[2]*fin[46]+alphax[1]*fin[45]+alphax[16]*fin[44]+alphax[0]*fin[42]+alphax[8]*fin[41]+alphax[7]*fin[40]+alphax[3]*fin[36]); 
  out[47] += 0.3061862178478971*(alphax[1]*fin[47]+alphax[0]*fin[46]+alphax[6]*fin[45]+alphax[7]*fin[44]+alphax[2]*fin[42]+alphax[3]*fin[41]+alphax[16]*fin[40]+alphax[8]*fin[36]); 

  double alphay[48] = {0.}; 
  alphay[0] = 2.0*speeds[8]*rdy2; 
  alphay[1] = 2.0*speeds[9]*rdy2; 
  alphay[2] = 2.0*speeds[10]*rdy2; 
  alphay[3] = 2.0*speeds[11]*rdy2; 
  alphay[6] = 2.0*speeds[12]*rdy2; 
  alphay[7] = 2.0*speeds[13]*rdy2; 
  alphay[8] = 2.0*speeds[14]*rdy2; 
  alphay[16] = 2.0*speeds[15]*rdy2; 


  out[2] += 0.3061862178478971*(alphay[16]*fin[16]+alphay[8]*fin[8]+alphay[7]*fin[7]+alphay[6]*fin[6]+alphay[3]*fin[3]+alphay[2]*fin[2]+alphay[1]*fin[1]+alphay[0]*fin[0]); 
  out[6] += 0.3061862178478971*(alphay[8]*fin[16]+fin[8]*alphay[16]+alphay[3]*fin[7]+fin[3]*alphay[7]+alphay[2]*fin[6]+fin[2]*alphay[6]+alphay[0]*fin[1]+fin[0]*alphay[1]); 
  out[8] += 0.3061862178478971*(alphay[6]*fin[16]+fin[6]*alphay[16]+alphay[2]*fin[8]+fin[2]*alphay[8]+alphay[1]*fin[7]+fin[1]*alphay[7]+alphay[0]*fin[3]+fin[0]*alphay[3]); 
  out[10] += 0.3061862178478971*(alphay[16]*fin[26]+alphay[8]*fin[19]+alphay[7]*fin[18]+alphay[6]*fin[17]+alphay[3]*fin[11]+alphay[2]*fin[10]+alphay[1]*fin[9]+alphay[0]*fin[4]); 
  out[13] += 0.3061862178478971*(alphay[16]*fin[27]+alphay[8]*fin[22]+alphay[7]*fin[21]+alphay[6]*fin[20]+alphay[3]*fin[14]+alphay[2]*fin[13]+alphay[1]*fin[12]+alphay[0]*fin[5]); 
  out[16] += 0.3061862178478971*(alphay[2]*fin[16]+fin[2]*alphay[16]+alphay[6]*fin[8]+fin[6]*alphay[8]+alphay[0]*fin[7]+fin[0]*alphay[7]+alphay[1]*fin[3]+fin[1]*alphay[3]); 
  out[17] += 0.3061862178478971*(alphay[8]*fin[26]+alphay[16]*fin[19]+alphay[3]*fin[18]+alphay[2]*fin[17]+alphay[7]*fin[11]+alphay[6]*fin[10]+alphay[0]*fin[9]+alphay[1]*fin[4]); 
  out[19] += 0.3061862178478971*(alphay[6]*fin[26]+alphay[2]*fin[19]+alphay[1]*fin[18]+alphay[16]*fin[17]+alphay[0]*fin[11]+alphay[8]*fin[10]+alphay[7]*fin[9]+alphay[3]*fin[4]); 
  out[20] += 0.3061862178478971*(alphay[8]*fin[27]+alphay[16]*fin[22]+alphay[3]*fin[21]+alphay[2]*fin[20]+alphay[7]*fin[14]+alphay[6]*fin[13]+alphay[0]*fin[12]+alphay[1]*fin[5]); 
  out[22] += 0.3061862178478971*(alphay[6]*fin[27]+alphay[2]*fin[22]+alphay[1]*fin[21]+alphay[16]*fin[20]+alphay[0]*fin[14]+alphay[8]*fin[13]+alphay[7]*fin[12]+alphay[3]*fin[5]); 
  out[24] += 0.3061862178478971*(alphay[16]*fin[31]+alphay[8]*fin[30]+alphay[7]*fin[29]+alphay[6]*fin[28]+alphay[3]*fin[25]+alphay[2]*fin[24]+alphay[1]*fin[23]+alphay[0]*fin[15]); 
  out[26] += 0.3061862178478971*(alphay[2]*fin[26]+alphay[6]*fin[19]+alphay[0]*fin[18]+alphay[8]*fin[17]+fin[10]*alphay[16]+alphay[1]*fin[11]+alphay[3]*fin[9]+fin[4]*alphay[7]); 
  out[27] += 0.3061862178478971*(alphay[2]*fin[27]+alphay[6]*fin[22]+alphay[0]*fin[21]+alphay[8]*fin[20]+fin[13]*alphay[16]+alphay[1]*fin[14]+alphay[3]*fin[12]+fin[5]*alphay[7]); 
  out[28] += 0.3061862178478971*(alphay[8]*fin[31]+alphay[16]*fin[30]+alphay[3]*fin[29]+alphay[2]*fin[28]+alphay[7]*fin[25]+alphay[6]*fin[24]+alphay[0]*fin[23]+alphay[1]*fin[15]); 
  out[30] += 0.3061862178478971*(alphay[6]*fin[31]+alphay[2]*fin[30]+alphay[1]*fin[29]+alphay[16]*fin[28]+alphay[0]*fin[25]+alphay[8]*fin[24]+alphay[7]*fin[23]+alphay[3]*fin[15]); 
  out[31] += 0.3061862178478971*(alphay[2]*fin[31]+alphay[6]*fin[30]+alphay[0]*fin[29]+alphay[8]*fin[28]+alphay[1]*fin[25]+alphay[16]*fin[24]+alphay[3]*fin[23]+alphay[7]*fin[15]); 
  out[34] += 0.3061862178478971*(alphay[16]*fin[43]+alphay[8]*fin[39]+alphay[7]*fin[38]+alphay[6]*fin[37]+alphay[3]*fin[35]+alphay[2]*fin[34]+alphay[1]*fin[33]+alphay[0]*fin[32]); 
  out[37] += 0.3061862178478971*(alphay[8]*fin[43]+alphay[16]*fin[39]+alphay[3]*fin[38]+alphay[2]*fin[37]+alphay[7]*fin[35]+alphay[6]*fin[34]+alphay[0]*fin[33]+alphay[1]*fin[32]); 
  out[39] += 0.3061862178478971*(alphay[6]*fin[43]+alphay[2]*fin[39]+alphay[1]*fin[38]+alphay[16]*fin[37]+alphay[0]*fin[35]+alphay[8]*fin[34]+alphay[7]*fin[33]+alphay[3]*fin[32]); 
  out[41] += 0.3061862178478971*(alphay[16]*fin[47]+alphay[8]*fin[46]+alphay[7]*fin[45]+alphay[6]*fin[44]+alphay[3]*fin[42]+alphay[2]*fin[41]+alphay[1]*fin[40]+alphay[0]*fin[36]); 
  out[43] += 0.3061862178478971*(alphay[2]*fin[43]+alphay[6]*fin[39]+alphay[0]*fin[38]+alphay[8]*fin[37]+alphay[1]*fin[35]+alphay[16]*fin[34]+alphay[3]*fin[33]+alphay[7]*fin[32]); 
  out[44] += 0.3061862178478971*(alphay[8]*fin[47]+alphay[16]*fin[46]+alphay[3]*fin[45]+alphay[2]*fin[44]+alphay[7]*fin[42]+alphay[6]*fin[41]+alphay[0]*fin[40]+alphay[1]*fin[36]); 
  out[46] += 0.3061862178478971*(alphay[6]*fin[47]+alphay[2]*fin[46]+alphay[1]*fin[45]+alphay[16]*fin[44]+alphay[0]*fin[42]+alphay[8]*fin[41]+alphay[7]*fin[40]+alphay[3]*fin[36]); 
  out[47] += 0.3061862178478971*(alphay[2]*fin[47]+alphay[6]*fin[46]+alphay[0]*fin[45]+alphay[8]*fin[44]+alphay[1]*fin[42]+alphay[16]*fin[41]+alphay[3]*fin[40]+alphay[7]*fin[36]); 

  double alphaz[48] = {0.}; 
  alphaz[0] = 2.0*speeds[16]*rdz2; 
  alphaz[1] = 2.0*speeds[17]*rdz2; 
  alphaz[2] = 2.0*speeds[18]*rdz2; 
  alphaz[3] = 2.0*speeds[19]*rdz2; 
  alphaz[6] = 2.0*speeds[20]*rdz2; 
  alphaz[7] = 2.0*speeds[21]*rdz2; 
  alphaz[8] = 2.0*speeds[22]*rdz2; 
  alphaz[16] = 2.0*speeds[23]*rdz2; 


  out[3] += 0.3061862178478971*(alphaz[16]*fin[16]+alphaz[8]*fin[8]+alphaz[7]*fin[7]+alphaz[6]*fin[6]+alphaz[3]*fin[3]+alphaz[2]*fin[2]+alphaz[1]*fin[1]+alphaz[0]*fin[0]); 
  out[7] += 0.3061862178478971*(alphaz[8]*fin[16]+fin[8]*alphaz[16]+alphaz[3]*fin[7]+fin[3]*alphaz[7]+alphaz[2]*fin[6]+fin[2]*alphaz[6]+alphaz[0]*fin[1]+fin[0]*alphaz[1]); 
  out[8] += 0.3061862178478971*(alphaz[7]*fin[16]+fin[7]*alphaz[16]+alphaz[3]*fin[8]+fin[3]*alphaz[8]+alphaz[1]*fin[6]+fin[1]*alphaz[6]+alphaz[0]*fin[2]+fin[0]*alphaz[2]); 
  out[11] += 0.3061862178478971*(alphaz[16]*fin[26]+alphaz[8]*fin[19]+alphaz[7]*fin[18]+alphaz[6]*fin[17]+alphaz[3]*fin[11]+alphaz[2]*fin[10]+alphaz[1]*fin[9]+alphaz[0]*fin[4]); 
  out[14] += 0.3061862178478971*(alphaz[16]*fin[27]+alphaz[8]*fin[22]+alphaz[7]*fin[21]+alphaz[6]*fin[20]+alphaz[3]*fin[14]+alphaz[2]*fin[13]+alphaz[1]*fin[12]+alphaz[0]*fin[5]); 
  out[16] += 0.3061862178478971*(alphaz[3]*fin[16]+fin[3]*alphaz[16]+alphaz[7]*fin[8]+fin[7]*alphaz[8]+alphaz[0]*fin[6]+fin[0]*alphaz[6]+alphaz[1]*fin[2]+fin[1]*alphaz[2]); 
  out[18] += 0.3061862178478971*(alphaz[8]*fin[26]+alphaz[16]*fin[19]+alphaz[3]*fin[18]+alphaz[2]*fin[17]+alphaz[7]*fin[11]+alphaz[6]*fin[10]+alphaz[0]*fin[9]+alphaz[1]*fin[4]); 
  out[19] += 0.3061862178478971*(alphaz[7]*fin[26]+alphaz[3]*fin[19]+alphaz[16]*fin[18]+alphaz[1]*fin[17]+alphaz[8]*fin[11]+alphaz[0]*fin[10]+alphaz[6]*fin[9]+alphaz[2]*fin[4]); 
  out[21] += 0.3061862178478971*(alphaz[8]*fin[27]+alphaz[16]*fin[22]+alphaz[3]*fin[21]+alphaz[2]*fin[20]+alphaz[7]*fin[14]+alphaz[6]*fin[13]+alphaz[0]*fin[12]+alphaz[1]*fin[5]); 
  out[22] += 0.3061862178478971*(alphaz[7]*fin[27]+alphaz[3]*fin[22]+alphaz[16]*fin[21]+alphaz[1]*fin[20]+alphaz[8]*fin[14]+alphaz[0]*fin[13]+alphaz[6]*fin[12]+alphaz[2]*fin[5]); 
  out[25] += 0.3061862178478971*(alphaz[16]*fin[31]+alphaz[8]*fin[30]+alphaz[7]*fin[29]+alphaz[6]*fin[28]+alphaz[3]*fin[25]+alphaz[2]*fin[24]+alphaz[1]*fin[23]+alphaz[0]*fin[15]); 
  out[26] += 0.3061862178478971*(alphaz[3]*fin[26]+alphaz[7]*fin[19]+alphaz[8]*fin[18]+alphaz[0]*fin[17]+fin[11]*alphaz[16]+alphaz[1]*fin[10]+alphaz[2]*fin[9]+fin[4]*alphaz[6]); 
  out[27] += 0.3061862178478971*(alphaz[3]*fin[27]+alphaz[7]*fin[22]+alphaz[8]*fin[21]+alphaz[0]*fin[20]+fin[14]*alphaz[16]+alphaz[1]*fin[13]+alphaz[2]*fin[12]+fin[5]*alphaz[6]); 
  out[29] += 0.3061862178478971*(alphaz[8]*fin[31]+alphaz[16]*fin[30]+alphaz[3]*fin[29]+alphaz[2]*fin[28]+alphaz[7]*fin[25]+alphaz[6]*fin[24]+alphaz[0]*fin[23]+alphaz[1]*fin[15]); 
  out[30] += 0.3061862178478971*(alphaz[7]*fin[31]+alphaz[3]*fin[30]+alphaz[16]*fin[29]+alphaz[1]*fin[28]+alphaz[8]*fin[25]+alphaz[0]*fin[24]+alphaz[6]*fin[23]+alphaz[2]*fin[15]); 
  out[31] += 0.3061862178478971*(alphaz[3]*fin[31]+alphaz[7]*fin[30]+alphaz[8]*fin[29]+alphaz[0]*fin[28]+alphaz[16]*fin[25]+alphaz[1]*fin[24]+alphaz[2]*fin[23]+alphaz[6]*fin[15]); 
  out[35] += 0.3061862178478971*(alphaz[16]*fin[43]+alphaz[8]*fin[39]+alphaz[7]*fin[38]+alphaz[6]*fin[37]+alphaz[3]*fin[35]+alphaz[2]*fin[34]+alphaz[1]*fin[33]+alphaz[0]*fin[32]); 
  out[38] += 0.3061862178478971*(alphaz[8]*fin[43]+alphaz[16]*fin[39]+alphaz[3]*fin[38]+alphaz[2]*fin[37]+alphaz[7]*fin[35]+alphaz[6]*fin[34]+alphaz[0]*fin[33]+alphaz[1]*fin[32]); 
  out[39] += 0.3061862178478971*(alphaz[7]*fin[43]+alphaz[3]*fin[39]+alphaz[16]*fin[38]+alphaz[1]*fin[37]+alphaz[8]*fin[35]+alphaz[0]*fin[34]+alphaz[6]*fin[33]+alphaz[2]*fin[32]); 
  out[42] += 0.3061862178478971*(alphaz[16]*fin[47]+alphaz[8]*fin[46]+alphaz[7]*fin[45]+alphaz[6]*fin[44]+alphaz[3]*fin[42]+alphaz[2]*fin[41]+alphaz[1]*fin[40]+alphaz[0]*fin[36]); 
  out[43] += 0.3061862178478971*(alphaz[3]*fin[43]+alphaz[7]*fin[39]+alphaz[8]*fin[38]+alphaz[0]*fin[37]+alphaz[16]*fin[35]+alphaz[1]*fin[34]+alphaz[2]*fin[33]+alphaz[6]*fin[32]); 
  out[45] += 0.3061862178478971*(alphaz[8]*fin[47]+alphaz[16]*fin[46]+alphaz[3]*fin[45]+alphaz[2]*fin[44]+alphaz[7]*fin[42]+alphaz[6]*fin[41]+alphaz[0]*fin[40]+alphaz[1]*fin[36]); 
  out[46] += 0.3061862178478971*(alphaz[7]*fin[47]+alphaz[3]*fin[46]+alphaz[16]*fin[45]+alphaz[1]*fin[44]+alphaz[8]*fin[42]+alphaz[0]*fin[41]+alphaz[6]*fin[40]+alphaz[2]*fin[36]); 
  out[47] += 0.3061862178478971*(alphaz[3]*fin[47]+alphaz[7]*fin[46]+alphaz[8]*fin[45]+alphaz[0]*fin[44]+alphaz[16]*fin[42]+alphaz[1]*fin[41]+alphaz[2]*fin[40]+alphaz[6]*fin[36]); 

  return 0.; 
} 
