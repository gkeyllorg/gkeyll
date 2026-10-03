#include <gkyl_dg_gyrokinetic_kernels.h> 
GKYL_CU_DH double dg_gyrokinetic_add_apardot_vol_1x2v_ser_p1(const double *vmap, const double q_, const double m_, const double *apardot,
    const double *fin, double* GKYL_RESTRICT out) 
{ 
  // q_,m_: species charge and mass.
  // apardot: time derivative of parallel component of magnetic vector potential.
  // fin: Distribution function.
  // out: output increment.

  double alphavpar[12] = {0.}; 
  alphavpar[0] = -((1.6329931618554527*apardot[0]*q_)/(vmap[1]*m_)); 
  alphavpar[1] = -((1.6329931618554527*apardot[1]*q_)/(vmap[1]*m_)); 


  out[2] += 0.6123724356957944*(alphavpar[1]*fin[1]+alphavpar[0]*fin[0]); 
  out[4] += 0.6123724356957944*(alphavpar[0]*fin[1]+fin[0]*alphavpar[1]); 
  out[6] += 0.6123724356957944*(alphavpar[1]*fin[5]+alphavpar[0]*fin[3]); 
  out[7] += 0.6123724356957944*(alphavpar[0]*fin[5]+alphavpar[1]*fin[3]); 
  out[8] += 1.369306393762915*(alphavpar[1]*fin[4]+alphavpar[0]*fin[2]); 
  out[9] += 1.369306393762915*(alphavpar[0]*fin[4]+alphavpar[1]*fin[2]); 
  out[10] += 1.369306393762915*(alphavpar[1]*fin[7]+alphavpar[0]*fin[6]); 
  out[11] += 1.369306393762915*(alphavpar[0]*fin[7]+alphavpar[1]*fin[6]); 

  return 0.; 
} 
