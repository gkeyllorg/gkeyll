#include <gkyl_dg_gyrokinetic_kernels.h> 
GKYL_CU_DH double dg_gyrokinetic_add_apar_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *vmapSq,
    const double q_, const double m_, const double *bmag, const double *jacobtot_inv, const double *dualcurlbhatoverB, const double *bioverJB,
    const double *b_i, const double *phi, const double *apar, const double *fin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // vmap: velocity space mapping.
  // vmapSq: velocity space mapping squared.
  // q_,m_: species charge and mass.
  // bmag: magnetic field amplitude.
  // jacobtot_inv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // dualcurlbhatoverB: dual curl of bhat over B.
  // bioverJB: b_i over J times B.
  // b_i: covariant components of the field aligned unit vector.
  // apar: parallel component of magnetic vector potential.
  // phi: electrostatic potential .
  // fin: Distribution function.
  // out: output increment.

  double rdx2 = 2.0/dxv[0];
  double rdvpar2 = 2.0/dxv[1];


  const double *b_1 = &b_i[0];
  const double *b_2 = &b_i[2];
  const double *b_3 = &b_i[4];

  const double *bioverJB_1 = &bioverJB[0]; 
  const double *bioverJB_2 = &bioverJB[2]; 
  const double *bioverJB_3 = &bioverJB[4]; 

  const double *dualcurlbhatoverB_1 = &dualcurlbhatoverB[0]; 
  const double *dualcurlbhatoverB_2 = &dualcurlbhatoverB[2]; 
  const double *dualcurlbhatoverB_3 = &dualcurlbhatoverB[4]; 

  double hamil[6] = {0.}; 
  hamil[0] = 1.4142135623730951*phi[0]*q_+0.7071067811865475*vmapSq[0]*m_; 
  hamil[1] = 1.4142135623730951*phi[1]*q_; 
  hamil[2] = 0.7071067811865475*vmapSq[1]*m_; 
  hamil[4] = 0.7071067811865475*vmapSq[2]*m_; 

  double vmap2 = vmap[1]*vmap[1]; 

  double hamil2[2] = {0.}; 
  hamil2[0] = hamil[2]*hamil[2]; 
  hamil2[1] = hamil[4]*hamil[4]; 

  double alphax[6] = {0.}; 
  alphax[0] = ((0.7071067811865475*apar[1]*dualcurlbhatoverB_3[1]*hamil[2]+0.7071067811865475*apar[0]*dualcurlbhatoverB_3[0]*hamil[2])*rdx2)/(vmap[1]*m_); 
  alphax[1] = ((0.7071067811865475*apar[0]*dualcurlbhatoverB_3[1]*hamil[2]+0.7071067811865475*dualcurlbhatoverB_3[0]*apar[1]*hamil[2])*rdx2)/(vmap[1]*m_); 
  alphax[2] = ((1.5811388300841895*apar[1]*dualcurlbhatoverB_3[1]*hamil[4]+1.5811388300841895*apar[0]*dualcurlbhatoverB_3[0]*hamil[4])*rdx2)/(vmap[1]*m_); 
  alphax[3] = ((1.5811388300841895*apar[0]*dualcurlbhatoverB_3[1]*hamil[4]+1.5811388300841895*dualcurlbhatoverB_3[0]*apar[1]*hamil[4])*rdx2)/(vmap[1]*m_); 


  out[1] += 0.8660254037844386*(alphax[3]*fin[3]+alphax[2]*fin[2]+alphax[1]*fin[1]+alphax[0]*fin[0]); 
  out[3] += 0.7745966692414834*alphax[3]*fin[5]+0.7745966692414833*alphax[2]*fin[4]+0.8660254037844386*(alphax[1]*fin[3]+fin[1]*alphax[3]+alphax[0]*fin[2]+fin[0]*alphax[2]); 
  out[5] += 0.8660254037844386*alphax[1]*fin[5]+0.8660254037844387*alphax[0]*fin[4]+0.7745966692414834*(alphax[3]*fin[3]+alphax[2]*fin[2]); 

  double alphavpar[6] = {0.}; 
  alphavpar[0] = ((-(0.7071067811865475*apar[1]*dualcurlbhatoverB_3[1]*hamil[1])-0.7071067811865475*apar[0]*dualcurlbhatoverB_3[0]*hamil[1])*rdx2)/(vmap[1]*m_); 
  alphavpar[1] = ((-(0.7071067811865475*apar[0]*dualcurlbhatoverB_3[1]*hamil[1])-0.7071067811865475*dualcurlbhatoverB_3[0]*apar[1]*hamil[1])*rdx2)/(vmap[1]*m_); 


  out[2] += 0.8660254037844386*(alphavpar[1]*fin[1]+alphavpar[0]*fin[0]); 
  out[3] += 0.8660254037844386*(alphavpar[0]*fin[1]+fin[0]*alphavpar[1]); 
  out[4] += 1.9364916731037085*(alphavpar[1]*fin[3]+alphavpar[0]*fin[2]); 
  out[5] += 1.9364916731037085*(alphavpar[0]*fin[3]+alphavpar[1]*fin[2]); 

  return 0.; 
} 
