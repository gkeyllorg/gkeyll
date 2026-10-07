#include <gkyl_dg_gyrokinetic_kernels.h> 
GKYL_CU_DH double dg_gyrokinetic_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *vmapSq,
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
  double rdvpar2 = 2.0/dxv[1];
  double rdmu2 = 2.0/dxv[2];

  double rdvpar2Sq = rdvpar2*rdvpar2;
  double dvparSq = dxv[1]*dxv[1];

  const double *bioverJB_1 = &bioverJB[0]; 
  const double *bioverJB_2 = &bioverJB[2]; 
  const double *bioverJB_3 = &bioverJB[4]; 

  const double *dualcurlbhatoverB_1 = &dualcurlbhatoverB[0]; 
  const double *dualcurlbhatoverB_2 = &dualcurlbhatoverB[2]; 
  const double *dualcurlbhatoverB_3 = &dualcurlbhatoverB[4]; 

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
  alphax[0] = (rdx2*((2.5*dualcurlbhatoverB_3[0]*hamil2[1])/vmap2+(0.5*dualcurlbhatoverB_3[0]*hamil2[0])/vmap2))/(m_*q_)+(rtg33inv[0]*vmap[1]*hamil[2]*rdx2)/(m_*vmap2); 
  alphax[1] = (rdx2*((2.5*dualcurlbhatoverB_3[1]*hamil2[1])/vmap2+(0.5*hamil2[0]*dualcurlbhatoverB_3[1])/vmap2))/(m_*q_)+(rtg33inv[1]*vmap[1]*hamil[2]*rdx2)/(m_*vmap2); 
  alphax[2] = (2.23606797749979*dualcurlbhatoverB_3[0]*hamil[2]*hamil[8]*rdx2)/(m_*q_*vmap2)+(2.23606797749979*rtg33inv[0]*vmap[1]*hamil[8]*rdx2)/(m_*vmap2); 
  alphax[4] = (2.23606797749979*dualcurlbhatoverB_3[1]*hamil[2]*hamil[8]*rdx2)/(m_*q_*vmap2)+(2.23606797749979*rtg33inv[1]*vmap[1]*hamil[8]*rdx2)/(m_*vmap2); 
  alphax[8] = (2.23606797749979*dualcurlbhatoverB_3[0]*hamil2[1]*rdx2)/(m_*q_*vmap2); 
  alphax[9] = (2.2360679774997902*dualcurlbhatoverB_3[1]*hamil2[1]*rdx2)/(m_*q_*vmap2); 


  out[1] += 0.6123724356957944*(alphax[9]*fin[9]+alphax[8]*fin[8]+alphax[4]*fin[4]+alphax[2]*fin[2]+alphax[1]*fin[1]+alphax[0]*fin[0]); 
  out[4] += 0.5477225575051661*(alphax[4]*fin[9]+fin[4]*alphax[9]+alphax[2]*fin[8]+fin[2]*alphax[8])+0.6123724356957944*(alphax[1]*fin[4]+fin[1]*alphax[4]+alphax[0]*fin[2]+fin[0]*alphax[2]); 
  out[5] += 0.6123724356957944*(alphax[9]*fin[11]+alphax[8]*fin[10]+alphax[4]*fin[7]+alphax[2]*fin[6]+alphax[1]*fin[5]+alphax[0]*fin[3]); 
  out[7] += 0.5477225575051661*(alphax[4]*fin[11]+alphax[2]*fin[10]+fin[7]*alphax[9]+fin[6]*alphax[8])+0.6123724356957944*(alphax[1]*fin[7]+alphax[0]*fin[6]+alphax[4]*fin[5]+alphax[2]*fin[3]); 
  out[9] += 0.39123039821797573*alphax[9]*fin[9]+0.6123724356957944*(alphax[1]*fin[9]+fin[1]*alphax[9])+0.39123039821797573*alphax[8]*fin[8]+0.6123724356957944*(alphax[0]*fin[8]+fin[0]*alphax[8])+0.5477225575051661*(alphax[4]*fin[4]+alphax[2]*fin[2]); 
  out[11] += (0.39123039821797573*alphax[9]+0.6123724356957944*alphax[1])*fin[11]+0.39123039821797573*alphax[8]*fin[10]+0.6123724356957944*(alphax[0]*fin[10]+fin[5]*alphax[9]+fin[3]*alphax[8])+0.5477225575051661*(alphax[4]*fin[7]+alphax[2]*fin[6]); 

  double alphavpar[12] = {0.}; 
  alphavpar[0] = -((0.5*dualcurlbhatoverB_3[0]*hamil[1]*hamil[2]*rdx2)/(m_*q_*vmap2))-(1.0*rtg33inv[0]*hamil[1]*vmap[1]*rdx2)/(m_*vmap2); 
  alphavpar[1] = -((0.5*dualcurlbhatoverB_3[1]*hamil[1]*hamil[2]*rdx2)/(m_*q_*vmap2))-(1.0*hamil[1]*rtg33inv[1]*vmap[1]*rdx2)/(m_*vmap2); 
  alphavpar[2] = -((1.118033988749895*dualcurlbhatoverB_3[0]*hamil[1]*hamil[8]*rdx2)/(m_*q_*vmap2)); 
  alphavpar[3] = -((0.5*dualcurlbhatoverB_3[0]*hamil[2]*hamil[5]*rdx2)/(m_*q_*vmap2))-(1.0*rtg33inv[0]*vmap[1]*hamil[5]*rdx2)/(m_*vmap2); 
  alphavpar[4] = -((1.118033988749895*dualcurlbhatoverB_3[1]*hamil[1]*hamil[8]*rdx2)/(m_*q_*vmap2)); 
  alphavpar[5] = -((0.5*dualcurlbhatoverB_3[1]*hamil[2]*hamil[5]*rdx2)/(m_*q_*vmap2))-(1.0*rtg33inv[1]*vmap[1]*hamil[5]*rdx2)/(m_*vmap2); 
  alphavpar[6] = -((1.118033988749895*dualcurlbhatoverB_3[0]*hamil[5]*hamil[8]*rdx2)/(m_*q_*vmap2)); 
  alphavpar[7] = -((1.118033988749895*dualcurlbhatoverB_3[1]*hamil[5]*hamil[8]*rdx2)/(m_*q_*vmap2)); 


  out[2] += 0.6123724356957944*(alphavpar[7]*fin[7]+alphavpar[6]*fin[6]+alphavpar[5]*fin[5]+alphavpar[4]*fin[4]+alphavpar[3]*fin[3]+alphavpar[2]*fin[2]+alphavpar[1]*fin[1]+alphavpar[0]*fin[0]); 
  out[4] += 0.6123724356957944*(alphavpar[6]*fin[7]+fin[6]*alphavpar[7]+alphavpar[3]*fin[5]+fin[3]*alphavpar[5]+alphavpar[2]*fin[4]+fin[2]*alphavpar[4]+alphavpar[0]*fin[1]+fin[0]*alphavpar[1]); 
  out[6] += 0.6123724356957944*(alphavpar[4]*fin[7]+fin[4]*alphavpar[7]+alphavpar[2]*fin[6]+fin[2]*alphavpar[6]+alphavpar[1]*fin[5]+fin[1]*alphavpar[5]+alphavpar[0]*fin[3]+fin[0]*alphavpar[3]); 
  out[7] += 0.6123724356957944*(alphavpar[2]*fin[7]+fin[2]*alphavpar[7]+alphavpar[4]*fin[6]+fin[4]*alphavpar[6]+alphavpar[0]*fin[5]+fin[0]*alphavpar[5]+alphavpar[1]*fin[3]+fin[1]*alphavpar[3]); 
  out[8] += 1.224744871391589*(alphavpar[7]*fin[11]+alphavpar[6]*fin[10]+alphavpar[4]*fin[9]+alphavpar[2]*fin[8])+1.369306393762915*(alphavpar[5]*fin[7]+fin[5]*alphavpar[7]+alphavpar[3]*fin[6]+fin[3]*alphavpar[6]+alphavpar[1]*fin[4]+fin[1]*alphavpar[4]+alphavpar[0]*fin[2]+fin[0]*alphavpar[2]); 
  out[9] += 1.224744871391589*(alphavpar[6]*fin[11]+alphavpar[7]*fin[10]+alphavpar[2]*fin[9]+alphavpar[4]*fin[8])+1.369306393762915*(alphavpar[3]*fin[7]+fin[3]*alphavpar[7]+alphavpar[5]*fin[6]+fin[5]*alphavpar[6]+alphavpar[0]*fin[4]+fin[0]*alphavpar[4]+alphavpar[1]*fin[2]+fin[1]*alphavpar[2]); 
  out[10] += 1.224744871391589*(alphavpar[4]*fin[11]+alphavpar[2]*fin[10]+alphavpar[7]*fin[9]+alphavpar[6]*fin[8])+1.369306393762915*(alphavpar[1]*fin[7]+fin[1]*alphavpar[7]+alphavpar[0]*fin[6]+fin[0]*alphavpar[6]+alphavpar[4]*fin[5]+fin[4]*alphavpar[5]+alphavpar[2]*fin[3]+fin[2]*alphavpar[3]); 
  out[11] += 1.224744871391589*(alphavpar[2]*fin[11]+alphavpar[4]*fin[10]+alphavpar[6]*fin[9]+alphavpar[7]*fin[8])+1.369306393762915*(alphavpar[0]*fin[7]+fin[0]*alphavpar[7]+alphavpar[1]*fin[6]+fin[1]*alphavpar[6]+alphavpar[2]*fin[5]+fin[2]*alphavpar[5]+alphavpar[3]*fin[4]+fin[3]*alphavpar[4]); 

  return 0.; 
} 
