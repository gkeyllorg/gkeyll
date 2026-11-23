#include <gkyl_dg_gyrokinetic_kernels.h> 
GKYL_CU_DH double dg_gyrokinetic_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *vmapSq,
    const double q_, const double m_, const double *bmag, const double *phi,
    const double *dualcurlbhatoverB, const double *B3overB, const double *bioverJB,
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
  double rdvpar2 = 2.0/dxv[1];
  double rdmu2 = 2.0/dxv[2];

  double rdvpar2Sq = rdvpar2*rdvpar2;
  double dvparSq = dxv[1]*dxv[1];

  const double *bioverJB_x = &bioverJB[0]; 
  const double *bioverJB_y = &bioverJB[2]; 
  const double *bioverJB_z = &bioverJB[4]; 

  const double *dualcurlbhatoverB_x = &dualcurlbhatoverB[0]; 
  const double *dualcurlbhatoverB_y = &dualcurlbhatoverB[2]; 
  const double *dualcurlbhatoverB_z = &dualcurlbhatoverB[4]; 

  double hamil[12] = {0.}; 
  hamil[0] = 2.0*phi[0]*q_+vmapSq[0]*m_+1.4142135623730951*bmag[0]*vmap[2]; 
  hamil[1] = 2.0*phi[1]*q_+1.4142135623730951*bmag[1]*vmap[2]; 
  hamil[2] = vmapSq[1]*m_; 
  hamil[3] = 1.4142135623730951*bmag[0]*vmap[3]; 
  hamil[5] = 1.4142135623730951*bmag[1]*vmap[3]; 
  hamil[8] = vmapSq[2]*m_; 

  double vmap2 = vmap[1]*vmap[1]; 

  double hamil2[2] = {0.}; 
  hamil2[0] = hamil[2]*hamil[2]; 
  hamil2[1] = hamil[8]*hamil[8]; 

  double alphax[12] = {0.}; 
  alphax[0] = (B3overB[0]*hamil[2]*rdx2)/(vmap[1]*m_); 
  alphax[1] = (B3overB[1]*hamil[2]*rdx2)/(vmap[1]*m_); 
  alphax[2] = (2.23606797749979*B3overB[0]*hamil[8]*rdx2)/(vmap[1]*m_); 
  alphax[4] = (2.23606797749979*B3overB[1]*hamil[8]*rdx2)/(vmap[1]*m_); 


  out[1] += 0.6123724356957944*(alphax[4]*fin[4]+alphax[2]*fin[2]+alphax[1]*fin[1]+alphax[0]*fin[0]); 
  out[4] += 0.5477225575051661*(alphax[4]*fin[9]+alphax[2]*fin[8])+0.6123724356957944*(alphax[1]*fin[4]+fin[1]*alphax[4]+alphax[0]*fin[2]+fin[0]*alphax[2]); 
  out[5] += 0.6123724356957944*(alphax[4]*fin[7]+alphax[2]*fin[6]+alphax[1]*fin[5]+alphax[0]*fin[3]); 
  out[7] += 0.5477225575051661*(alphax[4]*fin[11]+alphax[2]*fin[10])+0.6123724356957944*(alphax[1]*fin[7]+alphax[0]*fin[6]+alphax[4]*fin[5]+alphax[2]*fin[3]); 
  out[9] += 0.6123724356957944*(alphax[1]*fin[9]+alphax[0]*fin[8])+0.5477225575051661*(alphax[4]*fin[4]+alphax[2]*fin[2]); 
  out[11] += 0.6123724356957944*(alphax[1]*fin[11]+alphax[0]*fin[10])+0.5477225575051661*(alphax[4]*fin[7]+alphax[2]*fin[6]); 

  double alphavpar[12] = {0.}; 
  alphavpar[0] = -((1.0*B3overB[0]*hamil[1]*rdx2)/(vmap[1]*m_)); 
  alphavpar[1] = -((1.0*B3overB[1]*hamil[1]*rdx2)/(vmap[1]*m_)); 
  alphavpar[3] = -((1.0*B3overB[0]*hamil[5]*rdx2)/(vmap[1]*m_)); 
  alphavpar[5] = -((1.0*B3overB[1]*hamil[5]*rdx2)/(vmap[1]*m_)); 


  out[2] += 0.6123724356957944*(alphavpar[5]*fin[5]+alphavpar[3]*fin[3]+alphavpar[1]*fin[1]+alphavpar[0]*fin[0]); 
  out[4] += 0.6123724356957944*(alphavpar[3]*fin[5]+fin[3]*alphavpar[5]+alphavpar[0]*fin[1]+fin[0]*alphavpar[1]); 
  out[6] += 0.6123724356957944*(alphavpar[1]*fin[5]+fin[1]*alphavpar[5]+alphavpar[0]*fin[3]+fin[0]*alphavpar[3]); 
  out[7] += 0.6123724356957944*(alphavpar[0]*fin[5]+fin[0]*alphavpar[5]+alphavpar[1]*fin[3]+fin[1]*alphavpar[3]); 
  out[8] += 1.369306393762915*(alphavpar[5]*fin[7]+alphavpar[3]*fin[6]+alphavpar[1]*fin[4]+alphavpar[0]*fin[2]); 
  out[9] += 1.369306393762915*(alphavpar[3]*fin[7]+alphavpar[5]*fin[6]+alphavpar[0]*fin[4]+alphavpar[1]*fin[2]); 
  out[10] += 1.369306393762915*(alphavpar[1]*fin[7]+alphavpar[0]*fin[6]+fin[4]*alphavpar[5]+fin[2]*alphavpar[3]); 
  out[11] += 1.369306393762915*(alphavpar[0]*fin[7]+alphavpar[1]*fin[6]+fin[2]*alphavpar[5]+alphavpar[3]*fin[4]); 

  return 0.; 
} 
