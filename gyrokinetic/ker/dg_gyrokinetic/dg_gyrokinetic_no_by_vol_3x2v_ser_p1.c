#include <gkyl_dg_gyrokinetic_kernels.h> 
GKYL_CU_DH double dg_gyrokinetic_no_by_vol_3x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *vmapSq,
    const double q_, const double m_, const double *bmag, const double *phi,
    const double *dualcurlbhatoverB, const double *rtg33inv, const double *bioverJB,
    const double *fin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // vmap: velocity space mapping.
  // vmapSq: velocity space mapping squared.
  // q_,m_: species charge and mass.
  // bmag: magnetic field amplitude.
  // phi: electrostatic potential .
  // fin: Distribution function.
  // out: output increment.

  double rdx2 = 2.0/dxv[0];
  double rdy2 = 2.0/dxv[1];
  double rdz2 = 2.0/dxv[2];
  double rdvpar2 = 2.0/dxv[3];
  double rdmu2 = 2.0/dxv[4];

  double rdvpar2Sq = rdvpar2*rdvpar2;
  double dvparSq = dxv[3]*dxv[3];

  const double *bioverJB_1 = &bioverJB[0]; 
  const double *bioverJB_2 = &bioverJB[8]; 
  const double *bioverJB_3 = &bioverJB[16]; 

  const double *dualcurlbhatoverB_1 = &dualcurlbhatoverB[0]; 
  const double *dualcurlbhatoverB_2 = &dualcurlbhatoverB[8]; 
  const double *dualcurlbhatoverB_3 = &dualcurlbhatoverB[16]; 

  double hamil[48] = {0.}; 
  hamil[0] = 2.0*(phi[0]*q_+vmapSq[0]*m_)+1.4142135623730951*bmag[0]*vmap[2]; 
  hamil[1] = 2.0*phi[1]*q_+1.4142135623730951*bmag[1]*vmap[2]; 
  hamil[2] = 2.0*phi[2]*q_; 
  hamil[3] = 2.0*phi[3]*q_+1.4142135623730951*vmap[2]*bmag[3]; 
  hamil[4] = 2.0*vmapSq[1]*m_; 
  hamil[5] = 1.4142135623730951*bmag[0]*vmap[3]; 
  hamil[6] = 2.0*phi[4]*q_; 
  hamil[7] = 2.0*phi[5]*q_+1.4142135623730951*vmap[2]*bmag[5]; 
  hamil[8] = 2.0*phi[6]*q_; 
  hamil[12] = 1.4142135623730951*bmag[1]*vmap[3]; 
  hamil[14] = 1.4142135623730951*bmag[3]*vmap[3]; 
  hamil[16] = 2.0*phi[7]*q_; 
  hamil[21] = 1.4142135623730951*vmap[3]*bmag[5]; 
  hamil[32] = 2.0*vmapSq[2]*m_; 

  double vmap2 = vmap[1]*vmap[1]; 

  double hamil2[2] = {0.}; 
  hamil2[0] = hamil[4]*hamil[4]; 
  hamil2[1] = hamil[32]*hamil[32]; 

  double alphax[48] = {0.}; 



  double alphay[48] = {0.}; 



  double alphaz[48] = {0.}; 
  alphaz[0] = (0.5*rtg33inv[0]*hamil[4]*rdz2)/(vmap[1]*m_); 
  alphaz[1] = (0.5*rtg33inv[1]*hamil[4]*rdz2)/(vmap[1]*m_); 
  alphaz[3] = (0.5*rtg33inv[3]*hamil[4]*rdz2)/(vmap[1]*m_); 
  alphaz[4] = (1.118033988749895*rtg33inv[0]*hamil[32]*rdz2)/(vmap[1]*m_); 
  alphaz[7] = (0.5*hamil[4]*rtg33inv[5]*rdz2)/(vmap[1]*m_); 
  alphaz[9] = (1.118033988749895*rtg33inv[1]*hamil[32]*rdz2)/(vmap[1]*m_); 
  alphaz[11] = (1.118033988749895*rtg33inv[3]*hamil[32]*rdz2)/(vmap[1]*m_); 
  alphaz[18] = (1.118033988749895*rtg33inv[5]*hamil[32]*rdz2)/(vmap[1]*m_); 


  out[3] += 0.3061862178478971*(alphaz[18]*fin[18]+alphaz[11]*fin[11]+alphaz[9]*fin[9]+alphaz[7]*fin[7]+alphaz[4]*fin[4]+alphaz[3]*fin[3]+alphaz[1]*fin[1]+alphaz[0]*fin[0]); 
  out[7] += 0.3061862178478971*(alphaz[11]*fin[18]+fin[11]*alphaz[18]+alphaz[4]*fin[9]+fin[4]*alphaz[9]+alphaz[3]*fin[7]+fin[3]*alphaz[7]+alphaz[0]*fin[1]+fin[0]*alphaz[1]); 
  out[8] += 0.3061862178478971*(alphaz[18]*fin[26]+alphaz[11]*fin[19]+alphaz[9]*fin[17]+alphaz[7]*fin[16]+alphaz[4]*fin[10]+alphaz[3]*fin[8]+alphaz[1]*fin[6]+alphaz[0]*fin[2]); 
  out[11] += 0.273861278752583*alphaz[18]*fin[38]+0.27386127875258304*(alphaz[11]*fin[35]+alphaz[9]*fin[33])+0.273861278752583*alphaz[4]*fin[32]+0.3061862178478971*(alphaz[7]*fin[18]+fin[7]*alphaz[18]+alphaz[3]*fin[11]+fin[3]*alphaz[11]+alphaz[1]*fin[9]+fin[1]*alphaz[9]+alphaz[0]*fin[4]+fin[0]*alphaz[4]); 
  out[14] += 0.3061862178478971*(alphaz[18]*fin[29]+alphaz[11]*fin[25]+alphaz[9]*fin[23]+alphaz[7]*fin[21]+alphaz[4]*fin[15]+alphaz[3]*fin[14]+alphaz[1]*fin[12]+alphaz[0]*fin[5]); 
  out[16] += 0.3061862178478971*(alphaz[11]*fin[26]+alphaz[18]*fin[19]+alphaz[4]*fin[17]+alphaz[3]*fin[16]+alphaz[9]*fin[10]+alphaz[7]*fin[8]+alphaz[0]*fin[6]+alphaz[1]*fin[2]); 
  out[18] += 0.273861278752583*alphaz[11]*fin[38]+0.27386127875258304*(alphaz[18]*fin[35]+alphaz[4]*fin[33])+0.273861278752583*alphaz[9]*fin[32]+0.3061862178478971*(alphaz[3]*fin[18]+fin[3]*alphaz[18]+alphaz[7]*fin[11]+fin[7]*alphaz[11]+alphaz[0]*fin[9]+fin[0]*alphaz[9]+alphaz[1]*fin[4]+fin[1]*alphaz[4]); 
  out[19] += 0.27386127875258304*alphaz[18]*fin[43]+0.273861278752583*(alphaz[11]*fin[39]+alphaz[9]*fin[37])+0.27386127875258304*alphaz[4]*fin[34]+0.3061862178478971*(alphaz[7]*fin[26]+alphaz[3]*fin[19]+fin[16]*alphaz[18]+alphaz[1]*fin[17]+fin[8]*alphaz[11]+alphaz[0]*fin[10]+fin[6]*alphaz[9]+fin[2]*alphaz[4]); 
  out[21] += 0.3061862178478971*(alphaz[11]*fin[29]+alphaz[18]*fin[25]+alphaz[4]*fin[23]+alphaz[3]*fin[21]+alphaz[9]*fin[15]+alphaz[7]*fin[14]+alphaz[0]*fin[12]+alphaz[1]*fin[5]); 
  out[22] += 0.3061862178478971*(alphaz[18]*fin[31]+alphaz[11]*fin[30]+alphaz[9]*fin[28]+alphaz[7]*fin[27]+alphaz[4]*fin[24]+alphaz[3]*fin[22]+alphaz[1]*fin[20]+alphaz[0]*fin[13]); 
  out[25] += 0.27386127875258304*alphaz[18]*fin[45]+0.273861278752583*(alphaz[11]*fin[42]+alphaz[9]*fin[40])+0.27386127875258304*alphaz[4]*fin[36]+0.3061862178478971*(alphaz[7]*fin[29]+alphaz[3]*fin[25]+alphaz[1]*fin[23]+alphaz[18]*fin[21]+alphaz[0]*fin[15]+alphaz[11]*fin[14]+alphaz[9]*fin[12]+alphaz[4]*fin[5]); 
  out[26] += 0.27386127875258304*alphaz[11]*fin[43]+0.273861278752583*(alphaz[18]*fin[39]+alphaz[4]*fin[37])+0.27386127875258304*alphaz[9]*fin[34]+0.3061862178478971*(alphaz[3]*fin[26]+alphaz[7]*fin[19]+fin[8]*alphaz[18]+alphaz[0]*fin[17]+alphaz[11]*fin[16]+alphaz[1]*fin[10]+fin[2]*alphaz[9]+alphaz[4]*fin[6]); 
  out[27] += 0.3061862178478971*(alphaz[11]*fin[31]+alphaz[18]*fin[30]+alphaz[4]*fin[28]+alphaz[3]*fin[27]+alphaz[9]*fin[24]+alphaz[7]*fin[22]+alphaz[0]*fin[20]+alphaz[1]*fin[13]); 
  out[29] += 0.27386127875258304*alphaz[11]*fin[45]+0.273861278752583*(alphaz[18]*fin[42]+alphaz[4]*fin[40])+0.27386127875258304*alphaz[9]*fin[36]+0.3061862178478971*(alphaz[3]*fin[29]+alphaz[7]*fin[25]+alphaz[0]*fin[23]+alphaz[11]*fin[21]+fin[14]*alphaz[18]+alphaz[1]*fin[15]+alphaz[4]*fin[12]+fin[5]*alphaz[9]); 
  out[30] += 0.273861278752583*alphaz[18]*fin[47]+0.27386127875258304*(alphaz[11]*fin[46]+alphaz[9]*fin[44])+0.273861278752583*alphaz[4]*fin[41]+0.3061862178478971*(alphaz[7]*fin[31]+alphaz[3]*fin[30]+alphaz[1]*fin[28]+alphaz[18]*fin[27]+alphaz[0]*fin[24]+alphaz[11]*fin[22]+alphaz[9]*fin[20]+alphaz[4]*fin[13]); 
  out[31] += 0.273861278752583*alphaz[11]*fin[47]+0.27386127875258304*(alphaz[18]*fin[46]+alphaz[4]*fin[44])+0.273861278752583*alphaz[9]*fin[41]+0.3061862178478971*(alphaz[3]*fin[31]+alphaz[7]*fin[30]+alphaz[0]*fin[28]+alphaz[11]*fin[27]+alphaz[1]*fin[24]+alphaz[18]*fin[22]+alphaz[4]*fin[20]+alphaz[9]*fin[13]); 
  out[35] += 0.3061862178478971*(alphaz[7]*fin[38]+alphaz[3]*fin[35]+alphaz[1]*fin[33]+alphaz[0]*fin[32])+0.27386127875258304*(alphaz[18]*fin[18]+alphaz[11]*fin[11]+alphaz[9]*fin[9]+alphaz[4]*fin[4]); 
  out[38] += 0.3061862178478971*(alphaz[3]*fin[38]+alphaz[7]*fin[35]+alphaz[0]*fin[33]+alphaz[1]*fin[32])+0.273861278752583*(alphaz[11]*fin[18]+fin[11]*alphaz[18]+alphaz[4]*fin[9]+fin[4]*alphaz[9]); 
  out[39] += 0.3061862178478971*(alphaz[7]*fin[43]+alphaz[3]*fin[39]+alphaz[1]*fin[37]+alphaz[0]*fin[34])+0.273861278752583*(alphaz[18]*fin[26]+alphaz[11]*fin[19]+alphaz[9]*fin[17]+alphaz[4]*fin[10]); 
  out[42] += 0.3061862178478971*(alphaz[7]*fin[45]+alphaz[3]*fin[42]+alphaz[1]*fin[40]+alphaz[0]*fin[36])+0.273861278752583*(alphaz[18]*fin[29]+alphaz[11]*fin[25]+alphaz[9]*fin[23]+alphaz[4]*fin[15]); 
  out[43] += 0.3061862178478971*(alphaz[3]*fin[43]+alphaz[7]*fin[39]+alphaz[0]*fin[37]+alphaz[1]*fin[34])+0.27386127875258304*(alphaz[11]*fin[26]+alphaz[18]*fin[19]+alphaz[4]*fin[17]+alphaz[9]*fin[10]); 
  out[45] += 0.3061862178478971*(alphaz[3]*fin[45]+alphaz[7]*fin[42]+alphaz[0]*fin[40]+alphaz[1]*fin[36])+0.27386127875258304*(alphaz[11]*fin[29]+alphaz[18]*fin[25]+alphaz[4]*fin[23]+alphaz[9]*fin[15]); 
  out[46] += 0.3061862178478971*(alphaz[7]*fin[47]+alphaz[3]*fin[46]+alphaz[1]*fin[44]+alphaz[0]*fin[41])+0.27386127875258304*(alphaz[18]*fin[31]+alphaz[11]*fin[30]+alphaz[9]*fin[28]+alphaz[4]*fin[24]); 
  out[47] += 0.3061862178478971*(alphaz[3]*fin[47]+alphaz[7]*fin[46]+alphaz[0]*fin[44]+alphaz[1]*fin[41])+0.273861278752583*(alphaz[11]*fin[31]+alphaz[18]*fin[30]+alphaz[4]*fin[28]+alphaz[9]*fin[24]); 

  double alphavpar[48] = {0.}; 
  alphavpar[0] = ((-(0.5*rtg33inv[1]*hamil[7])-0.5*rtg33inv[0]*hamil[3])*rdz2)/(vmap[1]*m_); 
  alphavpar[1] = ((-(0.5*rtg33inv[0]*hamil[7])-0.5*rtg33inv[1]*hamil[3])*rdz2)/(vmap[1]*m_); 
  alphavpar[2] = ((-(0.5*rtg33inv[1]*hamil[16])-0.5*rtg33inv[0]*hamil[8])*rdz2)/(vmap[1]*m_); 
  alphavpar[3] = ((-(0.5*rtg33inv[5]*hamil[7])-0.5*hamil[3]*rtg33inv[3])*rdz2)/(vmap[1]*m_); 
  alphavpar[5] = ((-(0.5*rtg33inv[1]*hamil[21])-0.5*rtg33inv[0]*hamil[14])*rdz2)/(vmap[1]*m_); 
  alphavpar[6] = ((-(0.5*rtg33inv[0]*hamil[16])-0.5*rtg33inv[1]*hamil[8])*rdz2)/(vmap[1]*m_); 
  alphavpar[7] = ((-(0.5*rtg33inv[3]*hamil[7])-0.5*hamil[3]*rtg33inv[5])*rdz2)/(vmap[1]*m_); 
  alphavpar[8] = ((-(0.5*rtg33inv[5]*hamil[16])-0.5*rtg33inv[3]*hamil[8])*rdz2)/(vmap[1]*m_); 
  alphavpar[12] = ((-(0.5*rtg33inv[0]*hamil[21])-0.5*rtg33inv[1]*hamil[14])*rdz2)/(vmap[1]*m_); 
  alphavpar[14] = ((-(0.5*rtg33inv[5]*hamil[21])-0.5*rtg33inv[3]*hamil[14])*rdz2)/(vmap[1]*m_); 
  alphavpar[16] = ((-(0.5*rtg33inv[3]*hamil[16])-0.5*rtg33inv[5]*hamil[8])*rdz2)/(vmap[1]*m_); 
  alphavpar[21] = ((-(0.5*rtg33inv[3]*hamil[21])-0.5*rtg33inv[5]*hamil[14])*rdz2)/(vmap[1]*m_); 


  out[4] += 0.3061862178478971*(alphavpar[21]*fin[21]+alphavpar[16]*fin[16]+alphavpar[14]*fin[14]+alphavpar[12]*fin[12]+alphavpar[8]*fin[8]+alphavpar[7]*fin[7]+alphavpar[6]*fin[6]+alphavpar[5]*fin[5]+alphavpar[3]*fin[3]+alphavpar[2]*fin[2]+alphavpar[1]*fin[1]+alphavpar[0]*fin[0]); 
  out[9] += 0.3061862178478971*(alphavpar[14]*fin[21]+fin[14]*alphavpar[21]+alphavpar[8]*fin[16]+fin[8]*alphavpar[16]+alphavpar[5]*fin[12]+fin[5]*alphavpar[12]+alphavpar[3]*fin[7]+fin[3]*alphavpar[7]+alphavpar[2]*fin[6]+fin[2]*alphavpar[6]+alphavpar[0]*fin[1]+fin[0]*alphavpar[1]); 
  out[10] += 0.3061862178478971*(alphavpar[21]*fin[27]+alphavpar[14]*fin[22]+alphavpar[12]*fin[20]+alphavpar[7]*fin[16]+fin[7]*alphavpar[16]+alphavpar[5]*fin[13]+alphavpar[3]*fin[8]+fin[3]*alphavpar[8]+alphavpar[1]*fin[6]+fin[1]*alphavpar[6]+alphavpar[0]*fin[2]+fin[0]*alphavpar[2]); 
  out[11] += 0.3061862178478971*(alphavpar[12]*fin[21]+fin[12]*alphavpar[21]+alphavpar[6]*fin[16]+fin[6]*alphavpar[16]+alphavpar[5]*fin[14]+fin[5]*alphavpar[14]+alphavpar[2]*fin[8]+fin[2]*alphavpar[8]+alphavpar[1]*fin[7]+fin[1]*alphavpar[7]+alphavpar[0]*fin[3]+fin[0]*alphavpar[3]); 
  out[15] += 0.3061862178478971*(alphavpar[16]*fin[27]+alphavpar[8]*fin[22]+alphavpar[7]*fin[21]+fin[7]*alphavpar[21]+alphavpar[6]*fin[20]+alphavpar[3]*fin[14]+fin[3]*alphavpar[14]+alphavpar[2]*fin[13]+alphavpar[1]*fin[12]+fin[1]*alphavpar[12]+alphavpar[0]*fin[5]+fin[0]*alphavpar[5]); 
  out[17] += 0.3061862178478971*(alphavpar[14]*fin[27]+alphavpar[21]*fin[22]+alphavpar[5]*fin[20]+alphavpar[3]*fin[16]+fin[3]*alphavpar[16]+alphavpar[12]*fin[13]+alphavpar[7]*fin[8]+fin[7]*alphavpar[8]+alphavpar[0]*fin[6]+fin[0]*alphavpar[6]+alphavpar[1]*fin[2]+fin[1]*alphavpar[2]); 
  out[18] += 0.3061862178478971*(alphavpar[5]*fin[21]+fin[5]*alphavpar[21]+alphavpar[2]*fin[16]+fin[2]*alphavpar[16]+alphavpar[12]*fin[14]+fin[12]*alphavpar[14]+alphavpar[6]*fin[8]+fin[6]*alphavpar[8]+alphavpar[0]*fin[7]+fin[0]*alphavpar[7]+alphavpar[1]*fin[3]+fin[1]*alphavpar[3]); 
  out[19] += 0.3061862178478971*(alphavpar[12]*fin[27]+alphavpar[5]*fin[22]+fin[20]*alphavpar[21]+alphavpar[1]*fin[16]+fin[1]*alphavpar[16]+fin[13]*alphavpar[14]+alphavpar[0]*fin[8]+fin[0]*alphavpar[8]+alphavpar[6]*fin[7]+fin[6]*alphavpar[7]+alphavpar[2]*fin[3]+fin[2]*alphavpar[3]); 
  out[23] += 0.3061862178478971*(alphavpar[8]*fin[27]+alphavpar[16]*fin[22]+alphavpar[3]*fin[21]+fin[3]*alphavpar[21]+alphavpar[2]*fin[20]+alphavpar[7]*fin[14]+fin[7]*alphavpar[14]+alphavpar[6]*fin[13]+alphavpar[0]*fin[12]+fin[0]*alphavpar[12]+alphavpar[1]*fin[5]+fin[1]*alphavpar[5]); 
  out[24] += 0.3061862178478971*(alphavpar[7]*fin[27]+alphavpar[3]*fin[22]+alphavpar[16]*fin[21]+fin[16]*alphavpar[21]+alphavpar[1]*fin[20]+alphavpar[8]*fin[14]+fin[8]*alphavpar[14]+alphavpar[0]*fin[13]+alphavpar[6]*fin[12]+fin[6]*alphavpar[12]+alphavpar[2]*fin[5]+fin[2]*alphavpar[5]); 
  out[25] += 0.3061862178478971*(alphavpar[6]*fin[27]+alphavpar[2]*fin[22]+alphavpar[1]*fin[21]+fin[1]*alphavpar[21]+alphavpar[16]*fin[20]+alphavpar[0]*fin[14]+fin[0]*alphavpar[14]+alphavpar[8]*fin[13]+alphavpar[7]*fin[12]+fin[7]*alphavpar[12]+alphavpar[3]*fin[5]+fin[3]*alphavpar[5]); 
  out[26] += 0.3061862178478971*(alphavpar[5]*fin[27]+alphavpar[12]*fin[22]+fin[13]*alphavpar[21]+alphavpar[14]*fin[20]+alphavpar[0]*fin[16]+fin[0]*alphavpar[16]+alphavpar[1]*fin[8]+fin[1]*alphavpar[8]+alphavpar[2]*fin[7]+fin[2]*alphavpar[7]+alphavpar[3]*fin[6]+fin[3]*alphavpar[6]); 
  out[28] += 0.3061862178478971*(alphavpar[3]*fin[27]+alphavpar[7]*fin[22]+alphavpar[8]*fin[21]+fin[8]*alphavpar[21]+alphavpar[0]*fin[20]+alphavpar[14]*fin[16]+fin[14]*alphavpar[16]+alphavpar[1]*fin[13]+alphavpar[2]*fin[12]+fin[2]*alphavpar[12]+alphavpar[5]*fin[6]+fin[5]*alphavpar[6]); 
  out[29] += 0.3061862178478971*(alphavpar[2]*fin[27]+alphavpar[6]*fin[22]+alphavpar[0]*fin[21]+fin[0]*alphavpar[21]+alphavpar[8]*fin[20]+fin[13]*alphavpar[16]+alphavpar[1]*fin[14]+fin[1]*alphavpar[14]+alphavpar[3]*fin[12]+fin[3]*alphavpar[12]+alphavpar[5]*fin[7]+fin[5]*alphavpar[7]); 
  out[30] += 0.3061862178478971*(alphavpar[1]*fin[27]+alphavpar[0]*fin[22]+alphavpar[6]*fin[21]+fin[6]*alphavpar[21]+alphavpar[7]*fin[20]+alphavpar[12]*fin[16]+fin[12]*alphavpar[16]+alphavpar[2]*fin[14]+fin[2]*alphavpar[14]+alphavpar[3]*fin[13]+alphavpar[5]*fin[8]+fin[5]*alphavpar[8]); 
  out[31] += 0.3061862178478971*(alphavpar[0]*fin[27]+alphavpar[1]*fin[22]+alphavpar[2]*fin[21]+fin[2]*alphavpar[21]+alphavpar[3]*fin[20]+alphavpar[5]*fin[16]+fin[5]*alphavpar[16]+alphavpar[6]*fin[14]+fin[6]*alphavpar[14]+alphavpar[7]*fin[13]+alphavpar[8]*fin[12]+fin[8]*alphavpar[12]); 
  out[32] += 0.6846531968814573*(alphavpar[21]*fin[29]+alphavpar[16]*fin[26]+alphavpar[14]*fin[25]+alphavpar[12]*fin[23]+alphavpar[8]*fin[19]+alphavpar[7]*fin[18]+alphavpar[6]*fin[17]+alphavpar[5]*fin[15]+alphavpar[3]*fin[11]+alphavpar[2]*fin[10]+alphavpar[1]*fin[9]+alphavpar[0]*fin[4]); 
  out[33] += 0.6846531968814574*(alphavpar[14]*fin[29]+alphavpar[8]*fin[26]+alphavpar[21]*fin[25]+alphavpar[5]*fin[23]+alphavpar[16]*fin[19]+alphavpar[3]*fin[18]+alphavpar[2]*fin[17]+alphavpar[12]*fin[15]+alphavpar[7]*fin[11]+alphavpar[6]*fin[10]+alphavpar[0]*fin[9]+alphavpar[1]*fin[4]); 
  out[34] += 0.6846531968814574*(alphavpar[21]*fin[31]+alphavpar[14]*fin[30]+alphavpar[12]*fin[28]+alphavpar[7]*fin[26]+alphavpar[5]*fin[24]+alphavpar[3]*fin[19]+alphavpar[16]*fin[18]+alphavpar[1]*fin[17]+alphavpar[8]*fin[11]+alphavpar[0]*fin[10]+alphavpar[6]*fin[9]+alphavpar[2]*fin[4]); 
  out[35] += 0.6846531968814574*(alphavpar[12]*fin[29]+alphavpar[6]*fin[26]+alphavpar[5]*fin[25]+alphavpar[21]*fin[23]+alphavpar[2]*fin[19]+alphavpar[1]*fin[18]+alphavpar[16]*fin[17]+alphavpar[14]*fin[15]+alphavpar[0]*fin[11]+alphavpar[8]*fin[10]+alphavpar[7]*fin[9]+alphavpar[3]*fin[4]); 
  out[36] += 0.6846531968814574*(alphavpar[16]*fin[31]+alphavpar[8]*fin[30]+alphavpar[7]*fin[29]+alphavpar[6]*fin[28]+alphavpar[3]*fin[25]+alphavpar[2]*fin[24]+alphavpar[1]*fin[23]+fin[18]*alphavpar[21]+alphavpar[0]*fin[15]+fin[11]*alphavpar[14]+fin[9]*alphavpar[12]+fin[4]*alphavpar[5]); 
  out[37] += 0.6846531968814573*(alphavpar[14]*fin[31]+alphavpar[21]*fin[30]+alphavpar[5]*fin[28]+alphavpar[3]*fin[26]+alphavpar[12]*fin[24]+alphavpar[7]*fin[19]+alphavpar[8]*fin[18]+alphavpar[0]*fin[17]+fin[11]*alphavpar[16]+alphavpar[1]*fin[10]+alphavpar[2]*fin[9]+fin[4]*alphavpar[6]); 
  out[38] += 0.6846531968814573*(alphavpar[5]*fin[29]+alphavpar[2]*fin[26]+alphavpar[12]*fin[25]+alphavpar[14]*fin[23]+fin[15]*alphavpar[21]+alphavpar[6]*fin[19]+alphavpar[0]*fin[18]+alphavpar[8]*fin[17]+fin[10]*alphavpar[16]+alphavpar[1]*fin[11]+alphavpar[3]*fin[9]+fin[4]*alphavpar[7]); 
  out[39] += 0.6846531968814573*(alphavpar[12]*fin[31]+alphavpar[5]*fin[30]+alphavpar[21]*fin[28]+alphavpar[1]*fin[26]+alphavpar[14]*fin[24]+alphavpar[0]*fin[19]+alphavpar[6]*fin[18]+alphavpar[7]*fin[17]+fin[9]*alphavpar[16]+alphavpar[2]*fin[11]+alphavpar[3]*fin[10]+fin[4]*alphavpar[8]); 
  out[40] += 0.6846531968814573*(alphavpar[8]*fin[31]+alphavpar[16]*fin[30]+alphavpar[3]*fin[29]+alphavpar[2]*fin[28]+alphavpar[7]*fin[25]+alphavpar[6]*fin[24]+alphavpar[0]*fin[23]+fin[11]*alphavpar[21]+alphavpar[14]*fin[18]+alphavpar[1]*fin[15]+fin[4]*alphavpar[12]+alphavpar[5]*fin[9]); 
  out[41] += 0.6846531968814573*(alphavpar[7]*fin[31]+alphavpar[3]*fin[30]+alphavpar[16]*fin[29]+alphavpar[1]*fin[28]+alphavpar[21]*fin[26]+alphavpar[8]*fin[25]+alphavpar[0]*fin[24]+alphavpar[6]*fin[23]+alphavpar[14]*fin[19]+alphavpar[12]*fin[17]+alphavpar[2]*fin[15]+alphavpar[5]*fin[10]); 
  out[42] += 0.6846531968814573*(alphavpar[6]*fin[31]+alphavpar[2]*fin[30]+alphavpar[1]*fin[29]+alphavpar[16]*fin[28]+alphavpar[0]*fin[25]+alphavpar[8]*fin[24]+alphavpar[7]*fin[23]+fin[9]*alphavpar[21]+alphavpar[12]*fin[18]+alphavpar[3]*fin[15]+fin[4]*alphavpar[14]+alphavpar[5]*fin[11]); 
  out[43] += 0.6846531968814574*(alphavpar[5]*fin[31]+alphavpar[12]*fin[30]+alphavpar[14]*fin[28]+alphavpar[0]*fin[26]+alphavpar[21]*fin[24]+alphavpar[1]*fin[19]+alphavpar[2]*fin[18]+alphavpar[3]*fin[17]+fin[4]*alphavpar[16]+alphavpar[6]*fin[11]+alphavpar[7]*fin[10]+alphavpar[8]*fin[9]); 
  out[44] += 0.6846531968814574*(alphavpar[3]*fin[31]+alphavpar[7]*fin[30]+alphavpar[8]*fin[29]+alphavpar[0]*fin[28]+alphavpar[14]*fin[26]+alphavpar[16]*fin[25]+alphavpar[1]*fin[24]+alphavpar[2]*fin[23]+fin[19]*alphavpar[21]+alphavpar[5]*fin[17]+alphavpar[6]*fin[15]+fin[10]*alphavpar[12]); 
  out[45] += 0.6846531968814574*(alphavpar[2]*fin[31]+alphavpar[6]*fin[30]+alphavpar[0]*fin[29]+alphavpar[8]*fin[28]+alphavpar[1]*fin[25]+alphavpar[16]*fin[24]+alphavpar[3]*fin[23]+fin[4]*alphavpar[21]+alphavpar[5]*fin[18]+alphavpar[7]*fin[15]+fin[9]*alphavpar[14]+fin[11]*alphavpar[12]); 
  out[46] += 0.6846531968814574*(alphavpar[1]*fin[31]+alphavpar[0]*fin[30]+alphavpar[6]*fin[29]+alphavpar[7]*fin[28]+alphavpar[12]*fin[26]+alphavpar[2]*fin[25]+alphavpar[3]*fin[24]+alphavpar[16]*fin[23]+fin[17]*alphavpar[21]+alphavpar[5]*fin[19]+alphavpar[8]*fin[15]+fin[10]*alphavpar[14]); 
  out[47] += 0.6846531968814573*(alphavpar[0]*fin[31]+alphavpar[1]*fin[30]+alphavpar[2]*fin[29]+alphavpar[3]*fin[28]+alphavpar[5]*fin[26]+alphavpar[6]*fin[25]+alphavpar[7]*fin[24]+alphavpar[8]*fin[23]+fin[10]*alphavpar[21]+alphavpar[12]*fin[19]+alphavpar[14]*fin[17]+fin[15]*alphavpar[16]); 

  return 0.; 
} 
