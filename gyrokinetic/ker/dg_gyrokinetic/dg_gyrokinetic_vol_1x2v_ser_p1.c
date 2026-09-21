#include <gkyl_dg_gyrokinetic_kernels.h> 
GKYL_CU_DH double dg_gyrokinetic_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *vmapSq,
    const double q_, const double m_, const double *bmag, const double *yfield,
    const double *dualcurlbhatoverB, const double *rtg33inv, const double *bioverJB,
    const double *fin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // vmap: velocity space mapping.
  // vmapSq: velocity space mapping squared.
  // q_,m_: species charge and mass.
  // bmag: magnetic field amplitude.
  // yfield: Yushmanov field (gradient of Yushmanov potential).
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
  hamil[0] = vmapSq[0]*m_; 
  hamil[2] = vmapSq[1]*m_; 
  hamil[8] = vmapSq[2]*m_; 

  double vmap2 = vmap[1]*vmap[1]; 

  double hamil2[2] = {0.}; 
  hamil2[0] = hamil[2]*hamil[2]; 
  hamil2[1] = hamil[8]*hamil[8]; 

  double alphax[12] = {0.}; 
  alphax[0] = (rtg33inv[0]*hamil[2]*rdx2)/(vmap[1]*m_); 
  alphax[1] = (rtg33inv[1]*hamil[2]*rdx2)/(vmap[1]*m_); 
  alphax[2] = (2.23606797749979*rtg33inv[0]*hamil[8]*rdx2)/(vmap[1]*m_); 
  alphax[4] = (2.23606797749979*rtg33inv[1]*hamil[8]*rdx2)/(vmap[1]*m_); 


  out[1] += 0.6123724356957944*(alphax[4]*fin[4]+alphax[2]*fin[2]+alphax[1]*fin[1]+alphax[0]*fin[0]); 
  out[4] += 0.5477225575051661*(alphax[4]*fin[9]+alphax[2]*fin[8])+0.6123724356957944*(alphax[1]*fin[4]+fin[1]*alphax[4]+alphax[0]*fin[2]+fin[0]*alphax[2]); 
  out[5] += 0.6123724356957944*(alphax[4]*fin[7]+alphax[2]*fin[6]+alphax[1]*fin[5]+alphax[0]*fin[3]); 
  out[7] += 0.5477225575051661*(alphax[4]*fin[11]+alphax[2]*fin[10])+0.6123724356957944*(alphax[1]*fin[7]+alphax[0]*fin[6]+alphax[4]*fin[5]+alphax[2]*fin[3]); 
  out[9] += 0.6123724356957944*(alphax[1]*fin[9]+alphax[0]*fin[8])+0.5477225575051661*(alphax[4]*fin[4]+alphax[2]*fin[2]); 
  out[11] += 0.6123724356957944*(alphax[1]*fin[11]+alphax[0]*fin[10])+0.5477225575051661*(alphax[4]*fin[7]+alphax[2]*fin[6]); 

  double alphavpar[12] = {0.}; 
  alphavpar[0] = (-(0.5773502691896258*rtg33inv[1]*yfield[1])-0.5773502691896258*rtg33inv[0]*yfield[0])/(vmap[1]*m_); 
  alphavpar[1] = (-(0.5773502691896258*rtg33inv[0]*yfield[1])-0.5773502691896258*yfield[0]*rtg33inv[1])/(vmap[1]*m_); 
  alphavpar[2] = (-(0.5773502691896258*rtg33inv[1]*yfield[4])-0.5773502691896258*rtg33inv[0]*yfield[2])/(vmap[1]*m_); 
  alphavpar[3] = (-(0.5773502691896258*rtg33inv[1]*yfield[5])-0.5773502691896258*rtg33inv[0]*yfield[3])/(vmap[1]*m_); 
  alphavpar[4] = (-(0.5773502691896258*rtg33inv[0]*yfield[4])-0.5773502691896258*rtg33inv[1]*yfield[2])/(vmap[1]*m_); 
  alphavpar[5] = (-(0.5773502691896258*rtg33inv[0]*yfield[5])-0.5773502691896258*rtg33inv[1]*yfield[3])/(vmap[1]*m_); 
  alphavpar[6] = (-(0.5773502691896258*rtg33inv[1]*yfield[7])-0.5773502691896258*rtg33inv[0]*yfield[6])/(vmap[1]*m_); 
  alphavpar[7] = (-(0.5773502691896258*rtg33inv[0]*yfield[7])-0.5773502691896258*rtg33inv[1]*yfield[6])/(vmap[1]*m_); 
  alphavpar[8] = (-(0.5773502691896257*rtg33inv[1]*yfield[9])-0.5773502691896258*rtg33inv[0]*yfield[8])/(vmap[1]*m_); 
  alphavpar[9] = (-(0.5773502691896258*rtg33inv[0]*yfield[9])-0.5773502691896257*rtg33inv[1]*yfield[8])/(vmap[1]*m_); 
  alphavpar[10] = (-(0.5773502691896257*rtg33inv[1]*yfield[11])-0.5773502691896258*rtg33inv[0]*yfield[10])/(vmap[1]*m_); 
  alphavpar[11] = (-(0.5773502691896258*rtg33inv[0]*yfield[11])-0.5773502691896257*rtg33inv[1]*yfield[10])/(vmap[1]*m_); 


  out[2] += 0.6123724356957944*(alphavpar[11]*fin[11]+alphavpar[10]*fin[10]+alphavpar[9]*fin[9]+alphavpar[8]*fin[8]+alphavpar[7]*fin[7]+alphavpar[6]*fin[6]+alphavpar[5]*fin[5]+alphavpar[4]*fin[4]+alphavpar[3]*fin[3]+alphavpar[2]*fin[2]+alphavpar[1]*fin[1]+alphavpar[0]*fin[0]); 
  out[4] += 0.6123724356957944*(alphavpar[10]*fin[11]+fin[10]*alphavpar[11]+alphavpar[8]*fin[9]+fin[8]*alphavpar[9]+alphavpar[6]*fin[7]+fin[6]*alphavpar[7]+alphavpar[3]*fin[5]+fin[3]*alphavpar[5]+alphavpar[2]*fin[4]+fin[2]*alphavpar[4]+alphavpar[0]*fin[1]+fin[0]*alphavpar[1]); 
  out[6] += 0.6123724356957944*(alphavpar[9]*fin[11]+fin[9]*alphavpar[11]+alphavpar[8]*fin[10]+fin[8]*alphavpar[10]+alphavpar[4]*fin[7]+fin[4]*alphavpar[7]+alphavpar[2]*fin[6]+fin[2]*alphavpar[6]+alphavpar[1]*fin[5]+fin[1]*alphavpar[5]+alphavpar[0]*fin[3]+fin[0]*alphavpar[3]); 
  out[7] += 0.6123724356957944*(alphavpar[8]*fin[11]+fin[8]*alphavpar[11]+alphavpar[9]*fin[10]+fin[9]*alphavpar[10]+alphavpar[2]*fin[7]+fin[2]*alphavpar[7]+alphavpar[4]*fin[6]+fin[4]*alphavpar[6]+alphavpar[0]*fin[5]+fin[0]*alphavpar[5]+alphavpar[1]*fin[3]+fin[1]*alphavpar[3]); 
  out[8] += 1.224744871391589*(alphavpar[7]*fin[11]+fin[7]*alphavpar[11]+alphavpar[6]*fin[10]+fin[6]*alphavpar[10]+alphavpar[4]*fin[9]+fin[4]*alphavpar[9]+alphavpar[2]*fin[8]+fin[2]*alphavpar[8])+1.369306393762915*(alphavpar[5]*fin[7]+fin[5]*alphavpar[7]+alphavpar[3]*fin[6]+fin[3]*alphavpar[6]+alphavpar[1]*fin[4]+fin[1]*alphavpar[4]+alphavpar[0]*fin[2]+fin[0]*alphavpar[2]); 
  out[9] += 1.224744871391589*(alphavpar[6]*fin[11]+fin[6]*alphavpar[11]+alphavpar[7]*fin[10]+fin[7]*alphavpar[10]+alphavpar[2]*fin[9]+fin[2]*alphavpar[9]+alphavpar[4]*fin[8]+fin[4]*alphavpar[8])+1.369306393762915*(alphavpar[3]*fin[7]+fin[3]*alphavpar[7]+alphavpar[5]*fin[6]+fin[5]*alphavpar[6]+alphavpar[0]*fin[4]+fin[0]*alphavpar[4]+alphavpar[1]*fin[2]+fin[1]*alphavpar[2]); 
  out[10] += 1.224744871391589*(alphavpar[4]*fin[11]+fin[4]*alphavpar[11]+alphavpar[2]*fin[10]+fin[2]*alphavpar[10]+alphavpar[7]*fin[9]+fin[7]*alphavpar[9]+alphavpar[6]*fin[8]+fin[6]*alphavpar[8])+1.369306393762915*(alphavpar[1]*fin[7]+fin[1]*alphavpar[7]+alphavpar[0]*fin[6]+fin[0]*alphavpar[6]+alphavpar[4]*fin[5]+fin[4]*alphavpar[5]+alphavpar[2]*fin[3]+fin[2]*alphavpar[3]); 
  out[11] += 1.224744871391589*(alphavpar[2]*fin[11]+fin[2]*alphavpar[11]+alphavpar[4]*fin[10]+fin[4]*alphavpar[10]+alphavpar[6]*fin[9]+fin[6]*alphavpar[9]+alphavpar[7]*fin[8]+fin[7]*alphavpar[8])+1.369306393762915*(alphavpar[0]*fin[7]+fin[0]*alphavpar[7]+alphavpar[1]*fin[6]+fin[1]*alphavpar[6]+alphavpar[2]*fin[5]+fin[2]*alphavpar[5]+alphavpar[3]*fin[4]+fin[3]*alphavpar[4]); 

  return 0.; 
} 
