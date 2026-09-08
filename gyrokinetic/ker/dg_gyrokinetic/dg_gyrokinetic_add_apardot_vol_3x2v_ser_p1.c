#include <gkyl_dg_gyrokinetic_kernels.h> 
GKYL_CU_DH double dg_gyrokinetic_add_apardot_vol_3x2v_ser_p1(const double *vmap, const double q_, const double m_, const double *apardot,
    const double *fin, double* GKYL_RESTRICT out) 
{ 
  // q_,m_: species charge and mass.
  // apardot: time derivative of parallel component of magnetic vector potential.
  // fin: Distribution function.
  // out: output increment.

  double alphavpar[48] = {0.}; 
  alphavpar[0] = -((1.6329931618554527*apardot[0]*q_)/(vmap[1]*m_)); 
  alphavpar[1] = -((1.6329931618554527*apardot[1]*q_)/(vmap[1]*m_)); 
  alphavpar[2] = -((1.6329931618554527*apardot[2]*q_)/(vmap[1]*m_)); 
  alphavpar[3] = -((1.6329931618554527*apardot[3]*q_)/(vmap[1]*m_)); 
  alphavpar[6] = -((1.6329931618554527*apardot[4]*q_)/(vmap[1]*m_)); 
  alphavpar[7] = -((1.6329931618554527*apardot[5]*q_)/(vmap[1]*m_)); 
  alphavpar[8] = -((1.6329931618554527*apardot[6]*q_)/(vmap[1]*m_)); 
  alphavpar[16] = -((1.6329931618554527*apardot[7]*q_)/(vmap[1]*m_)); 


  out[4] += 0.3061862178478971*(alphavpar[16]*fin[16]+alphavpar[8]*fin[8]+alphavpar[7]*fin[7]+alphavpar[6]*fin[6]+alphavpar[3]*fin[3]+alphavpar[2]*fin[2]+alphavpar[1]*fin[1]+alphavpar[0]*fin[0]); 
  out[9] += 0.3061862178478971*(alphavpar[8]*fin[16]+fin[8]*alphavpar[16]+alphavpar[3]*fin[7]+fin[3]*alphavpar[7]+alphavpar[2]*fin[6]+fin[2]*alphavpar[6]+alphavpar[0]*fin[1]+fin[0]*alphavpar[1]); 
  out[10] += 0.3061862178478971*(alphavpar[7]*fin[16]+fin[7]*alphavpar[16]+alphavpar[3]*fin[8]+fin[3]*alphavpar[8]+alphavpar[1]*fin[6]+fin[1]*alphavpar[6]+alphavpar[0]*fin[2]+fin[0]*alphavpar[2]); 
  out[11] += 0.3061862178478971*(alphavpar[6]*fin[16]+fin[6]*alphavpar[16]+alphavpar[2]*fin[8]+fin[2]*alphavpar[8]+alphavpar[1]*fin[7]+fin[1]*alphavpar[7]+alphavpar[0]*fin[3]+fin[0]*alphavpar[3]); 
  out[15] += 0.3061862178478971*(alphavpar[16]*fin[27]+alphavpar[8]*fin[22]+alphavpar[7]*fin[21]+alphavpar[6]*fin[20]+alphavpar[3]*fin[14]+alphavpar[2]*fin[13]+alphavpar[1]*fin[12]+alphavpar[0]*fin[5]); 
  out[17] += 0.3061862178478971*(alphavpar[3]*fin[16]+fin[3]*alphavpar[16]+alphavpar[7]*fin[8]+fin[7]*alphavpar[8]+alphavpar[0]*fin[6]+fin[0]*alphavpar[6]+alphavpar[1]*fin[2]+fin[1]*alphavpar[2]); 
  out[18] += 0.3061862178478971*(alphavpar[2]*fin[16]+fin[2]*alphavpar[16]+alphavpar[6]*fin[8]+fin[6]*alphavpar[8]+alphavpar[0]*fin[7]+fin[0]*alphavpar[7]+alphavpar[1]*fin[3]+fin[1]*alphavpar[3]); 
  out[19] += 0.3061862178478971*(alphavpar[1]*fin[16]+fin[1]*alphavpar[16]+alphavpar[0]*fin[8]+fin[0]*alphavpar[8]+alphavpar[6]*fin[7]+fin[6]*alphavpar[7]+alphavpar[2]*fin[3]+fin[2]*alphavpar[3]); 
  out[23] += 0.3061862178478971*(alphavpar[8]*fin[27]+alphavpar[16]*fin[22]+alphavpar[3]*fin[21]+alphavpar[2]*fin[20]+alphavpar[7]*fin[14]+alphavpar[6]*fin[13]+alphavpar[0]*fin[12]+alphavpar[1]*fin[5]); 
  out[24] += 0.3061862178478971*(alphavpar[7]*fin[27]+alphavpar[3]*fin[22]+alphavpar[16]*fin[21]+alphavpar[1]*fin[20]+alphavpar[8]*fin[14]+alphavpar[0]*fin[13]+alphavpar[6]*fin[12]+alphavpar[2]*fin[5]); 
  out[25] += 0.3061862178478971*(alphavpar[6]*fin[27]+alphavpar[2]*fin[22]+alphavpar[1]*fin[21]+alphavpar[16]*fin[20]+alphavpar[0]*fin[14]+alphavpar[8]*fin[13]+alphavpar[7]*fin[12]+alphavpar[3]*fin[5]); 
  out[26] += 0.3061862178478971*(alphavpar[0]*fin[16]+fin[0]*alphavpar[16]+alphavpar[1]*fin[8]+fin[1]*alphavpar[8]+alphavpar[2]*fin[7]+fin[2]*alphavpar[7]+alphavpar[3]*fin[6]+fin[3]*alphavpar[6]); 
  out[28] += 0.3061862178478971*(alphavpar[3]*fin[27]+alphavpar[7]*fin[22]+alphavpar[8]*fin[21]+alphavpar[0]*fin[20]+fin[14]*alphavpar[16]+alphavpar[1]*fin[13]+alphavpar[2]*fin[12]+fin[5]*alphavpar[6]); 
  out[29] += 0.3061862178478971*(alphavpar[2]*fin[27]+alphavpar[6]*fin[22]+alphavpar[0]*fin[21]+alphavpar[8]*fin[20]+fin[13]*alphavpar[16]+alphavpar[1]*fin[14]+alphavpar[3]*fin[12]+fin[5]*alphavpar[7]); 
  out[30] += 0.3061862178478971*(alphavpar[1]*fin[27]+alphavpar[0]*fin[22]+alphavpar[6]*fin[21]+alphavpar[7]*fin[20]+fin[12]*alphavpar[16]+alphavpar[2]*fin[14]+alphavpar[3]*fin[13]+fin[5]*alphavpar[8]); 
  out[31] += 0.3061862178478971*(alphavpar[0]*fin[27]+alphavpar[1]*fin[22]+alphavpar[2]*fin[21]+alphavpar[3]*fin[20]+fin[5]*alphavpar[16]+alphavpar[6]*fin[14]+alphavpar[7]*fin[13]+alphavpar[8]*fin[12]); 
  out[32] += 0.6846531968814573*(alphavpar[16]*fin[26]+alphavpar[8]*fin[19]+alphavpar[7]*fin[18]+alphavpar[6]*fin[17]+alphavpar[3]*fin[11]+alphavpar[2]*fin[10]+alphavpar[1]*fin[9]+alphavpar[0]*fin[4]); 
  out[33] += 0.6846531968814574*(alphavpar[8]*fin[26]+alphavpar[16]*fin[19]+alphavpar[3]*fin[18]+alphavpar[2]*fin[17]+alphavpar[7]*fin[11]+alphavpar[6]*fin[10]+alphavpar[0]*fin[9]+alphavpar[1]*fin[4]); 
  out[34] += 0.6846531968814574*(alphavpar[7]*fin[26]+alphavpar[3]*fin[19]+alphavpar[16]*fin[18]+alphavpar[1]*fin[17]+alphavpar[8]*fin[11]+alphavpar[0]*fin[10]+alphavpar[6]*fin[9]+alphavpar[2]*fin[4]); 
  out[35] += 0.6846531968814574*(alphavpar[6]*fin[26]+alphavpar[2]*fin[19]+alphavpar[1]*fin[18]+alphavpar[16]*fin[17]+alphavpar[0]*fin[11]+alphavpar[8]*fin[10]+alphavpar[7]*fin[9]+alphavpar[3]*fin[4]); 
  out[36] += 0.6846531968814574*(alphavpar[16]*fin[31]+alphavpar[8]*fin[30]+alphavpar[7]*fin[29]+alphavpar[6]*fin[28]+alphavpar[3]*fin[25]+alphavpar[2]*fin[24]+alphavpar[1]*fin[23]+alphavpar[0]*fin[15]); 
  out[37] += 0.6846531968814573*(alphavpar[3]*fin[26]+alphavpar[7]*fin[19]+alphavpar[8]*fin[18]+alphavpar[0]*fin[17]+fin[11]*alphavpar[16]+alphavpar[1]*fin[10]+alphavpar[2]*fin[9]+fin[4]*alphavpar[6]); 
  out[38] += 0.6846531968814573*(alphavpar[2]*fin[26]+alphavpar[6]*fin[19]+alphavpar[0]*fin[18]+alphavpar[8]*fin[17]+fin[10]*alphavpar[16]+alphavpar[1]*fin[11]+alphavpar[3]*fin[9]+fin[4]*alphavpar[7]); 
  out[39] += 0.6846531968814573*(alphavpar[1]*fin[26]+alphavpar[0]*fin[19]+alphavpar[6]*fin[18]+alphavpar[7]*fin[17]+fin[9]*alphavpar[16]+alphavpar[2]*fin[11]+alphavpar[3]*fin[10]+fin[4]*alphavpar[8]); 
  out[40] += 0.6846531968814573*(alphavpar[8]*fin[31]+alphavpar[16]*fin[30]+alphavpar[3]*fin[29]+alphavpar[2]*fin[28]+alphavpar[7]*fin[25]+alphavpar[6]*fin[24]+alphavpar[0]*fin[23]+alphavpar[1]*fin[15]); 
  out[41] += 0.6846531968814573*(alphavpar[7]*fin[31]+alphavpar[3]*fin[30]+alphavpar[16]*fin[29]+alphavpar[1]*fin[28]+alphavpar[8]*fin[25]+alphavpar[0]*fin[24]+alphavpar[6]*fin[23]+alphavpar[2]*fin[15]); 
  out[42] += 0.6846531968814573*(alphavpar[6]*fin[31]+alphavpar[2]*fin[30]+alphavpar[1]*fin[29]+alphavpar[16]*fin[28]+alphavpar[0]*fin[25]+alphavpar[8]*fin[24]+alphavpar[7]*fin[23]+alphavpar[3]*fin[15]); 
  out[43] += 0.6846531968814574*(alphavpar[0]*fin[26]+alphavpar[1]*fin[19]+alphavpar[2]*fin[18]+alphavpar[3]*fin[17]+fin[4]*alphavpar[16]+alphavpar[6]*fin[11]+alphavpar[7]*fin[10]+alphavpar[8]*fin[9]); 
  out[44] += 0.6846531968814574*(alphavpar[3]*fin[31]+alphavpar[7]*fin[30]+alphavpar[8]*fin[29]+alphavpar[0]*fin[28]+alphavpar[16]*fin[25]+alphavpar[1]*fin[24]+alphavpar[2]*fin[23]+alphavpar[6]*fin[15]); 
  out[45] += 0.6846531968814574*(alphavpar[2]*fin[31]+alphavpar[6]*fin[30]+alphavpar[0]*fin[29]+alphavpar[8]*fin[28]+alphavpar[1]*fin[25]+alphavpar[16]*fin[24]+alphavpar[3]*fin[23]+alphavpar[7]*fin[15]); 
  out[46] += 0.6846531968814574*(alphavpar[1]*fin[31]+alphavpar[0]*fin[30]+alphavpar[6]*fin[29]+alphavpar[7]*fin[28]+alphavpar[2]*fin[25]+alphavpar[3]*fin[24]+alphavpar[16]*fin[23]+alphavpar[8]*fin[15]); 
  out[47] += 0.6846531968814573*(alphavpar[0]*fin[31]+alphavpar[1]*fin[30]+alphavpar[2]*fin[29]+alphavpar[3]*fin[28]+alphavpar[6]*fin[25]+alphavpar[7]*fin[24]+alphavpar[8]*fin[23]+fin[15]*alphavpar[16]); 

  return 0.; 
} 
