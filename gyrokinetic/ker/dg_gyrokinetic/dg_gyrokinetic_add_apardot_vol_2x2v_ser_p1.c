#include <gkyl_dg_gyrokinetic_kernels.h> 
GKYL_CU_DH double dg_gyrokinetic_add_apardot_vol_2x2v_ser_p1(const double *vmap, const double q_, const double m_, const double *apardot,
    const double *fin, double* GKYL_RESTRICT out) 
{ 
  // q_,m_: species charge and mass.
  // apardot: time derivative of parallel component of magnetic vector potential.
  // fin: Distribution function.
  // out: output increment.

  double alphavpar[24] = {0.}; 
  alphavpar[0] = -((1.6329931618554527*apardot[0]*q_)/(vmap[1]*m_)); 
  alphavpar[1] = -((1.6329931618554527*apardot[1]*q_)/(vmap[1]*m_)); 
  alphavpar[2] = -((1.6329931618554527*apardot[2]*q_)/(vmap[1]*m_)); 
  alphavpar[5] = -((1.6329931618554527*apardot[3]*q_)/(vmap[1]*m_)); 


  out[3] += 0.4330127018922193*(alphavpar[5]*fin[5]+alphavpar[2]*fin[2]+alphavpar[1]*fin[1]+alphavpar[0]*fin[0]); 
  out[6] += 0.4330127018922193*(alphavpar[2]*fin[5]+fin[2]*alphavpar[5]+alphavpar[0]*fin[1]+fin[0]*alphavpar[1]); 
  out[7] += 0.4330127018922193*(alphavpar[1]*fin[5]+fin[1]*alphavpar[5]+alphavpar[0]*fin[2]+fin[0]*alphavpar[2]); 
  out[10] += 0.4330127018922193*(alphavpar[5]*fin[12]+alphavpar[2]*fin[9]+alphavpar[1]*fin[8]+alphavpar[0]*fin[4]); 
  out[11] += 0.4330127018922193*(alphavpar[0]*fin[5]+fin[0]*alphavpar[5]+alphavpar[1]*fin[2]+fin[1]*alphavpar[2]); 
  out[13] += 0.4330127018922193*(alphavpar[2]*fin[12]+alphavpar[5]*fin[9]+alphavpar[0]*fin[8]+alphavpar[1]*fin[4]); 
  out[14] += 0.4330127018922193*(alphavpar[1]*fin[12]+alphavpar[0]*fin[9]+alphavpar[5]*fin[8]+alphavpar[2]*fin[4]); 
  out[15] += 0.4330127018922193*(alphavpar[0]*fin[12]+alphavpar[1]*fin[9]+alphavpar[2]*fin[8]+fin[4]*alphavpar[5]); 
  out[16] += 0.9682458365518543*(alphavpar[5]*fin[11]+alphavpar[2]*fin[7]+alphavpar[1]*fin[6]+alphavpar[0]*fin[3]); 
  out[17] += 0.9682458365518543*(alphavpar[2]*fin[11]+alphavpar[5]*fin[7]+alphavpar[0]*fin[6]+alphavpar[1]*fin[3]); 
  out[18] += 0.9682458365518543*(alphavpar[1]*fin[11]+alphavpar[0]*fin[7]+alphavpar[5]*fin[6]+alphavpar[2]*fin[3]); 
  out[19] += 0.9682458365518543*(alphavpar[5]*fin[15]+alphavpar[2]*fin[14]+alphavpar[1]*fin[13]+alphavpar[0]*fin[10]); 
  out[20] += 0.9682458365518543*(alphavpar[0]*fin[11]+alphavpar[1]*fin[7]+alphavpar[2]*fin[6]+fin[3]*alphavpar[5]); 
  out[21] += 0.9682458365518543*(alphavpar[2]*fin[15]+alphavpar[5]*fin[14]+alphavpar[0]*fin[13]+alphavpar[1]*fin[10]); 
  out[22] += 0.9682458365518543*(alphavpar[1]*fin[15]+alphavpar[0]*fin[14]+alphavpar[5]*fin[13]+alphavpar[2]*fin[10]); 
  out[23] += 0.9682458365518543*(alphavpar[0]*fin[15]+alphavpar[1]*fin[14]+alphavpar[2]*fin[13]+alphavpar[5]*fin[10]); 

  return 0.; 
} 
