#include <gkyl_dg_gyrokinetic_passive_kernels.h> 
GKYL_CU_DH double dg_gyrokinetic_passive_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *vmapSq,
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

  double alphax[12] = {0.}; 
  alphax[0] = 2.0*speeds[0]*rdx2; 
  alphax[1] = 2.0*speeds[1]*rdx2; 


  out[1] += 0.6123724356957944*(alphax[1]*fin[1]+alphax[0]*fin[0]); 
  out[4] += 0.6123724356957944*(alphax[1]*fin[4]+alphax[0]*fin[2]); 
  out[5] += 0.6123724356957944*(alphax[1]*fin[5]+alphax[0]*fin[3]); 
  out[7] += 0.6123724356957944*(alphax[1]*fin[7]+alphax[0]*fin[6]); 
  out[9] += 0.6123724356957944*(alphax[1]*fin[9]+alphax[0]*fin[8]); 
  out[11] += 0.6123724356957944*(alphax[1]*fin[11]+alphax[0]*fin[10]); 

  return 0.; 
} 
