#include <gkyl_dg_gyrokinetic_kernels.h> 
GKYL_CU_DH double dg_gyrokinetic_add_apardot_vol_1x1v_ser_p1(const double *vmap, const double q_, const double m_, const double *apardot,
    const double *fin, double* GKYL_RESTRICT out) 
{ 
  // q_,m_: species charge and mass.
  // apardot: time derivative of parallel component of magnetic vector potential.
  // fin: Distribution function.
  // out: output increment.

  double alphavpar[6] = {0.}; 
  alphavpar[0] = -((1.1547005383792517*apardot[0]*q_)/(vmap[1]*m_)); 
  alphavpar[1] = -((1.1547005383792517*apardot[1]*q_)/(vmap[1]*m_)); 


  out[2] += 0.8660254037844386*(alphavpar[1]*fin[1]+alphavpar[0]*fin[0]); 
  out[3] += 0.8660254037844386*(alphavpar[0]*fin[1]+fin[0]*alphavpar[1]); 
  out[4] += 1.9364916731037085*(alphavpar[1]*fin[3]+alphavpar[0]*fin[2]); 
  out[5] += 1.9364916731037085*(alphavpar[0]*fin[3]+alphavpar[1]*fin[2]); 

  return 0.; 
} 
