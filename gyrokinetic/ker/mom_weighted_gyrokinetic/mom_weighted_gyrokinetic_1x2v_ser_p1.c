#include <gkyl_mom_weighted_gyrokinetic_kernels.h> 
GKYL_CU_DH void mom_weighted_gyrokinetic_M0_WEIGHTX_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[1]*dxv[2]/m_; 
  out[0] += (1.4142135623730951*f[1]*wgt[1]+1.4142135623730951*f[0]*wgt[0])*volFact; 
  out[1] += (1.4142135623730951*f[0]*wgt[1]+1.4142135623730951*wgt[0]*f[1])*volFact; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M0_WEIGHTV_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[1]*dxv[2]/m_; 
  out[0] += (wgt[5]*f[10]+wgt[4]*f[8]+wgt[3]*f[6]+wgt[2]*f[3]+wgt[1]*f[2]+f[0]*wgt[0])*volFact; 
  out[1] += (1.0000000000000002*wgt[5]*f[11]+1.0000000000000002*wgt[4]*f[9]+wgt[3]*f[7]+wgt[2]*f[5]+wgt[1]*f[4]+wgt[0]*f[1])*volFact; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M1_WEIGHTX_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[1]*dxv[2]/m_; 
  double tmp[2]; 
  tmp[0] = 1.4142135623730951*vmap[1]*f[2]+1.4142135623730951*f[0]*vmap[0]; 
  tmp[1] = 1.4142135623730951*vmap[1]*f[4]+1.4142135623730951*vmap[0]*f[1]; 
  out[0] += (0.7071067811865475*tmp[1]*wgt[1]+0.7071067811865475*tmp[0]*wgt[0])*volFact; 
  out[1] += (0.7071067811865475*tmp[0]*wgt[1]+0.7071067811865475*wgt[0]*tmp[1])*volFact; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M1_WEIGHTV_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[1]*dxv[2]/m_; 
  double tmp[2]; 
  tmp[0] = 0.7071067811865475*vmap[0]*wgt[5]*f[10]+0.632455532033676*vmap[1]*wgt[3]*f[10]+0.7071067811865475*vmap[0]*wgt[4]*f[8]+0.6324555320336759*vmap[1]*wgt[1]*f[8]+0.632455532033676*vmap[1]*wgt[5]*f[6]+0.7071067811865475*vmap[0]*wgt[3]*f[6]+0.7071067811865475*vmap[1]*wgt[2]*f[6]+0.6324555320336759*vmap[1]*f[2]*wgt[4]+0.7071067811865475*vmap[1]*f[3]*wgt[3]+0.7071067811865475*vmap[0]*wgt[2]*f[3]+0.7071067811865475*vmap[0]*wgt[1]*f[2]+0.7071067811865475*wgt[0]*vmap[1]*f[2]+0.7071067811865475*f[0]*vmap[1]*wgt[1]+0.7071067811865475*f[0]*vmap[0]*wgt[0]; 
  tmp[1] = 0.7071067811865475*vmap[0]*wgt[5]*f[11]+0.6324555320336759*vmap[1]*wgt[3]*f[11]+0.7071067811865475*vmap[0]*wgt[4]*f[9]+0.632455532033676*vmap[1]*wgt[1]*f[9]+0.632455532033676*vmap[1]*wgt[5]*f[7]+0.7071067811865475*vmap[0]*wgt[3]*f[7]+0.7071067811865475*vmap[1]*wgt[2]*f[7]+0.7071067811865475*vmap[1]*wgt[3]*f[5]+0.7071067811865475*vmap[0]*wgt[2]*f[5]+0.6324555320336759*vmap[1]*f[4]*wgt[4]+0.7071067811865475*vmap[0]*wgt[1]*f[4]+0.7071067811865475*wgt[0]*vmap[1]*f[4]+0.7071067811865475*f[1]*vmap[1]*wgt[1]+0.7071067811865475*vmap[0]*wgt[0]*f[1]; 
  out[0] += tmp[0]*volFact; 
  out[1] += tmp[1]*volFact; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M2_WEIGHTX_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[1]*dxv[2]/m_; 
  double tmpPar[2]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  tmpPar[0] = 0.8944271909999159*vmap1R2*f[8]+2.0*vmap[0]*vmap[1]*f[2]+f[0]*vmap1R2+f[0]*vmap0R2; 
  tmpPar[1] = 0.8944271909999161*vmap1R2*f[9]+2.0*vmap[0]*vmap[1]*f[4]+f[1]*vmap1R2+vmap0R2*f[1]; 

  out[0] += (0.7071067811865475*tmpPar[1]*wgt[1]+0.7071067811865475*tmpPar[0]*wgt[0])*volFact; 
  out[1] += (0.7071067811865475*tmpPar[0]*wgt[1]+0.7071067811865475*wgt[0]*tmpPar[1])*volFact; 
  double tmpPerp[2]; 
  tmpPerp[0] = 1.4142135623730951*f[3]*vmap[3]+1.4142135623730951*f[0]*vmap[2]; 
  tmpPerp[1] = 1.4142135623730951*vmap[3]*f[5]+1.4142135623730951*f[1]*vmap[2]; 
  out[0] += (2.0*(0.5*bmag[0]*tmpPerp[1]*wgt[1]+0.5*tmpPerp[0]*bmag[1]*wgt[1]+0.5*wgt[0]*bmag[1]*tmpPerp[1]+0.5*bmag[0]*tmpPerp[0]*wgt[0])*volFact)/m_; 
  out[1] += (2.0*(0.9*bmag[1]*tmpPerp[1]*wgt[1]+0.5*bmag[0]*tmpPerp[0]*wgt[1]+0.5*bmag[0]*wgt[0]*tmpPerp[1]+0.5*tmpPerp[0]*wgt[0]*bmag[1])*volFact)/m_; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M2_WEIGHTV_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[1]*dxv[2]/m_; 
  double tmpPar[2]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  tmpPar[0] = 0.7857142857142857*vmap1R2*wgt[5]*f[10]+0.5*vmap0R2*wgt[5]*f[10]+0.8944271909999161*vmap[0]*vmap[1]*wgt[3]*f[10]+0.44721359549995804*vmap1R2*wgt[2]*f[10]+0.7857142857142857*vmap1R2*wgt[4]*f[8]+0.5*vmap0R2*wgt[4]*f[8]+0.8944271909999159*vmap[0]*vmap[1]*wgt[1]*f[8]+0.4472135954999579*wgt[0]*vmap1R2*f[8]+0.8944271909999161*vmap[0]*vmap[1]*wgt[5]*f[6]+0.9*vmap1R2*wgt[3]*f[6]+0.5*vmap0R2*wgt[3]*f[6]+vmap[0]*vmap[1]*wgt[2]*f[6]+0.44721359549995804*vmap1R2*f[3]*wgt[5]+0.8944271909999159*vmap[0]*vmap[1]*f[2]*wgt[4]+0.4472135954999579*f[0]*vmap1R2*wgt[4]+vmap[0]*vmap[1]*f[3]*wgt[3]+0.5*vmap1R2*wgt[2]*f[3]+0.5*vmap0R2*wgt[2]*f[3]+0.9*vmap1R2*wgt[1]*f[2]+0.5*vmap0R2*wgt[1]*f[2]+vmap[0]*wgt[0]*vmap[1]*f[2]+f[0]*vmap[0]*vmap[1]*wgt[1]+0.5*f[0]*wgt[0]*vmap1R2+0.5*f[0]*vmap0R2*wgt[0]; 
  tmpPar[1] = 0.7857142857142858*vmap1R2*wgt[5]*f[11]+0.5000000000000001*vmap0R2*wgt[5]*f[11]+0.8944271909999159*vmap[0]*vmap[1]*wgt[3]*f[11]+0.4472135954999579*vmap1R2*wgt[2]*f[11]+0.7857142857142858*vmap1R2*wgt[4]*f[9]+0.5000000000000001*vmap0R2*wgt[4]*f[9]+0.8944271909999161*vmap[0]*vmap[1]*wgt[1]*f[9]+0.44721359549995804*wgt[0]*vmap1R2*f[9]+0.8944271909999161*vmap[0]*vmap[1]*wgt[5]*f[7]+0.9*vmap1R2*wgt[3]*f[7]+0.5*vmap0R2*wgt[3]*f[7]+vmap[0]*vmap[1]*wgt[2]*f[7]+0.44721359549995804*vmap1R2*f[5]*wgt[5]+vmap[0]*vmap[1]*wgt[3]*f[5]+0.5*vmap1R2*wgt[2]*f[5]+0.5*vmap0R2*wgt[2]*f[5]+0.8944271909999159*vmap[0]*vmap[1]*f[4]*wgt[4]+0.4472135954999579*f[1]*vmap1R2*wgt[4]+0.9*vmap1R2*wgt[1]*f[4]+0.5*vmap0R2*wgt[1]*f[4]+vmap[0]*wgt[0]*vmap[1]*f[4]+vmap[0]*f[1]*vmap[1]*wgt[1]+0.5*wgt[0]*f[1]*vmap1R2+0.5*vmap0R2*wgt[0]*f[1]; 

  out[0] += tmpPar[0]*volFact; 
  out[1] += tmpPar[1]*volFact; 
  double tmpPerp[2]; 
  tmpPerp[0] = 0.7071067811865475*vmap[2]*wgt[5]*f[10]+0.7071067811865475*vmap[3]*wgt[4]*f[10]+0.7071067811865475*vmap[3]*wgt[5]*f[8]+0.7071067811865475*vmap[2]*wgt[4]*f[8]+0.7071067811865475*vmap[2]*wgt[3]*f[6]+0.7071067811865475*wgt[1]*vmap[3]*f[6]+0.7071067811865475*f[2]*vmap[3]*wgt[3]+0.7071067811865475*wgt[0]*f[3]*vmap[3]+0.7071067811865475*f[0]*wgt[2]*vmap[3]+0.7071067811865475*vmap[2]*wgt[2]*f[3]+0.7071067811865475*wgt[1]*f[2]*vmap[2]+0.7071067811865475*f[0]*wgt[0]*vmap[2]; 
  tmpPerp[1] = 0.7071067811865475*vmap[2]*wgt[5]*f[11]+0.7071067811865475*vmap[3]*wgt[4]*f[11]+0.7071067811865475*vmap[3]*wgt[5]*f[9]+0.7071067811865475*vmap[2]*wgt[4]*f[9]+0.7071067811865475*vmap[2]*wgt[3]*f[7]+0.7071067811865475*wgt[1]*vmap[3]*f[7]+0.7071067811865475*wgt[0]*vmap[3]*f[5]+0.7071067811865475*vmap[2]*wgt[2]*f[5]+0.7071067811865475*vmap[3]*wgt[3]*f[4]+0.7071067811865475*wgt[1]*vmap[2]*f[4]+0.7071067811865475*f[1]*wgt[2]*vmap[3]+0.7071067811865475*wgt[0]*f[1]*vmap[2]; 
  out[0] += (2.0*(0.7071067811865475*bmag[1]*tmpPerp[1]+0.7071067811865475*bmag[0]*tmpPerp[0])*volFact)/m_; 
  out[1] += (2.0*(0.7071067811865475*bmag[0]*tmpPerp[1]+0.7071067811865475*tmpPerp[0]*bmag[1])*volFact)/m_; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M2PAR_WEIGHTX_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[1]*dxv[2]/m_; 
  double tmp[2]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  tmp[0] = 0.8944271909999159*vmap1R2*f[8]+2.0*vmap[0]*vmap[1]*f[2]+f[0]*vmap1R2+f[0]*vmap0R2; 
  tmp[1] = 0.8944271909999161*vmap1R2*f[9]+2.0*vmap[0]*vmap[1]*f[4]+f[1]*vmap1R2+vmap0R2*f[1]; 

  out[0] += (0.7071067811865475*tmp[1]*wgt[1]+0.7071067811865475*tmp[0]*wgt[0])*volFact; 
  out[1] += (0.7071067811865475*tmp[0]*wgt[1]+0.7071067811865475*wgt[0]*tmp[1])*volFact; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M2PAR_WEIGHTV_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[1]*dxv[2]/m_; 
  double tmp[2]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  tmp[0] = 0.7857142857142857*vmap1R2*wgt[5]*f[10]+0.5*vmap0R2*wgt[5]*f[10]+0.8944271909999161*vmap[0]*vmap[1]*wgt[3]*f[10]+0.44721359549995804*vmap1R2*wgt[2]*f[10]+0.7857142857142857*vmap1R2*wgt[4]*f[8]+0.5*vmap0R2*wgt[4]*f[8]+0.8944271909999159*vmap[0]*vmap[1]*wgt[1]*f[8]+0.4472135954999579*wgt[0]*vmap1R2*f[8]+0.8944271909999161*vmap[0]*vmap[1]*wgt[5]*f[6]+0.9*vmap1R2*wgt[3]*f[6]+0.5*vmap0R2*wgt[3]*f[6]+vmap[0]*vmap[1]*wgt[2]*f[6]+0.44721359549995804*vmap1R2*f[3]*wgt[5]+0.8944271909999159*vmap[0]*vmap[1]*f[2]*wgt[4]+0.4472135954999579*f[0]*vmap1R2*wgt[4]+vmap[0]*vmap[1]*f[3]*wgt[3]+0.5*vmap1R2*wgt[2]*f[3]+0.5*vmap0R2*wgt[2]*f[3]+0.9*vmap1R2*wgt[1]*f[2]+0.5*vmap0R2*wgt[1]*f[2]+vmap[0]*wgt[0]*vmap[1]*f[2]+f[0]*vmap[0]*vmap[1]*wgt[1]+0.5*f[0]*wgt[0]*vmap1R2+0.5*f[0]*vmap0R2*wgt[0]; 
  tmp[1] = 0.7857142857142858*vmap1R2*wgt[5]*f[11]+0.5000000000000001*vmap0R2*wgt[5]*f[11]+0.8944271909999159*vmap[0]*vmap[1]*wgt[3]*f[11]+0.4472135954999579*vmap1R2*wgt[2]*f[11]+0.7857142857142858*vmap1R2*wgt[4]*f[9]+0.5000000000000001*vmap0R2*wgt[4]*f[9]+0.8944271909999161*vmap[0]*vmap[1]*wgt[1]*f[9]+0.44721359549995804*wgt[0]*vmap1R2*f[9]+0.8944271909999161*vmap[0]*vmap[1]*wgt[5]*f[7]+0.9*vmap1R2*wgt[3]*f[7]+0.5*vmap0R2*wgt[3]*f[7]+vmap[0]*vmap[1]*wgt[2]*f[7]+0.44721359549995804*vmap1R2*f[5]*wgt[5]+vmap[0]*vmap[1]*wgt[3]*f[5]+0.5*vmap1R2*wgt[2]*f[5]+0.5*vmap0R2*wgt[2]*f[5]+0.8944271909999159*vmap[0]*vmap[1]*f[4]*wgt[4]+0.4472135954999579*f[1]*vmap1R2*wgt[4]+0.9*vmap1R2*wgt[1]*f[4]+0.5*vmap0R2*wgt[1]*f[4]+vmap[0]*wgt[0]*vmap[1]*f[4]+vmap[0]*f[1]*vmap[1]*wgt[1]+0.5*wgt[0]*f[1]*vmap1R2+0.5*vmap0R2*wgt[0]*f[1]; 

  out[0] += tmp[0]*volFact; 
  out[1] += tmp[1]*volFact; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M0M1M2_WEIGHTX_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[1]*dxv[2]/m_; 
  double tmp1[2]; 
  tmp1[0] = 1.4142135623730951*vmap[1]*f[2]+1.4142135623730951*f[0]*vmap[0]; 
  tmp1[1] = 1.4142135623730951*vmap[1]*f[4]+1.4142135623730951*vmap[0]*f[1]; 
  double tmp2par[2]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  tmp2par[0] = 0.8944271909999159*vmap1R2*f[8]+2.0*vmap[0]*vmap[1]*f[2]+f[0]*vmap1R2+f[0]*vmap0R2; 
  tmp2par[1] = 0.8944271909999161*vmap1R2*f[9]+2.0*vmap[0]*vmap[1]*f[4]+f[1]*vmap1R2+vmap0R2*f[1]; 
  double tmp2perp[2]; 
  tmp2perp[0] = 1.4142135623730951*f[3]*vmap[3]+1.4142135623730951*f[0]*vmap[2]; 
  tmp2perp[1] = 1.4142135623730951*vmap[3]*f[5]+1.4142135623730951*f[1]*vmap[2]; 

  out[0] += (1.4142135623730951*f[1]*wgt[1]+1.4142135623730951*f[0]*wgt[0])*volFact; 
  out[1] += (1.4142135623730951*f[0]*wgt[1]+1.4142135623730951*wgt[0]*f[1])*volFact; 
  out[2] += (0.7071067811865475*tmp1[1]*wgt[1]+0.7071067811865475*tmp1[0]*wgt[0])*volFact; 
  out[3] += (0.7071067811865475*tmp1[0]*wgt[1]+0.7071067811865475*wgt[0]*tmp1[1])*volFact; 
  out[4] += ((bmag[0]*tmp2perp[1]*wgt[1]+tmp2perp[0]*bmag[1]*wgt[1]+wgt[0]*bmag[1]*tmp2perp[1]+bmag[0]*tmp2perp[0]*wgt[0])*volFact)/m_+(0.7071067811865475*tmp2par[1]*wgt[1]+0.7071067811865475*tmp2par[0]*wgt[0])*volFact; 
  out[5] += ((1.8*bmag[1]*tmp2perp[1]*wgt[1]+bmag[0]*tmp2perp[0]*wgt[1]+bmag[0]*wgt[0]*tmp2perp[1]+tmp2perp[0]*wgt[0]*bmag[1])*volFact)/m_+(0.7071067811865475*tmp2par[0]*wgt[1]+0.7071067811865475*wgt[0]*tmp2par[1])*volFact; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M0M1M2_WEIGHTV_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[1]*dxv[2]/m_; 
  double tmp1[2]; 
  tmp1[0] = 0.7071067811865475*vmap[0]*wgt[5]*f[10]+0.632455532033676*vmap[1]*wgt[3]*f[10]+0.7071067811865475*vmap[0]*wgt[4]*f[8]+0.6324555320336759*vmap[1]*wgt[1]*f[8]+0.632455532033676*vmap[1]*wgt[5]*f[6]+0.7071067811865475*vmap[0]*wgt[3]*f[6]+0.7071067811865475*vmap[1]*wgt[2]*f[6]+0.6324555320336759*vmap[1]*f[2]*wgt[4]+0.7071067811865475*vmap[1]*f[3]*wgt[3]+0.7071067811865475*vmap[0]*wgt[2]*f[3]+0.7071067811865475*vmap[0]*wgt[1]*f[2]+0.7071067811865475*wgt[0]*vmap[1]*f[2]+0.7071067811865475*f[0]*vmap[1]*wgt[1]+0.7071067811865475*f[0]*vmap[0]*wgt[0]; 
  tmp1[1] = 0.7071067811865475*vmap[0]*wgt[5]*f[11]+0.6324555320336759*vmap[1]*wgt[3]*f[11]+0.7071067811865475*vmap[0]*wgt[4]*f[9]+0.632455532033676*vmap[1]*wgt[1]*f[9]+0.632455532033676*vmap[1]*wgt[5]*f[7]+0.7071067811865475*vmap[0]*wgt[3]*f[7]+0.7071067811865475*vmap[1]*wgt[2]*f[7]+0.7071067811865475*vmap[1]*wgt[3]*f[5]+0.7071067811865475*vmap[0]*wgt[2]*f[5]+0.6324555320336759*vmap[1]*f[4]*wgt[4]+0.7071067811865475*vmap[0]*wgt[1]*f[4]+0.7071067811865475*wgt[0]*vmap[1]*f[4]+0.7071067811865475*f[1]*vmap[1]*wgt[1]+0.7071067811865475*vmap[0]*wgt[0]*f[1]; 
  double tmp2par[2]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  tmp2par[0] = 0.7857142857142857*vmap1R2*wgt[5]*f[10]+0.5*vmap0R2*wgt[5]*f[10]+0.8944271909999161*vmap[0]*vmap[1]*wgt[3]*f[10]+0.44721359549995804*vmap1R2*wgt[2]*f[10]+0.7857142857142857*vmap1R2*wgt[4]*f[8]+0.5*vmap0R2*wgt[4]*f[8]+0.8944271909999159*vmap[0]*vmap[1]*wgt[1]*f[8]+0.4472135954999579*wgt[0]*vmap1R2*f[8]+0.8944271909999161*vmap[0]*vmap[1]*wgt[5]*f[6]+0.9*vmap1R2*wgt[3]*f[6]+0.5*vmap0R2*wgt[3]*f[6]+vmap[0]*vmap[1]*wgt[2]*f[6]+0.44721359549995804*vmap1R2*f[3]*wgt[5]+0.8944271909999159*vmap[0]*vmap[1]*f[2]*wgt[4]+0.4472135954999579*f[0]*vmap1R2*wgt[4]+vmap[0]*vmap[1]*f[3]*wgt[3]+0.5*vmap1R2*wgt[2]*f[3]+0.5*vmap0R2*wgt[2]*f[3]+0.9*vmap1R2*wgt[1]*f[2]+0.5*vmap0R2*wgt[1]*f[2]+vmap[0]*wgt[0]*vmap[1]*f[2]+f[0]*vmap[0]*vmap[1]*wgt[1]+0.5*f[0]*wgt[0]*vmap1R2+0.5*f[0]*vmap0R2*wgt[0]; 
  tmp2par[1] = 0.7857142857142858*vmap1R2*wgt[5]*f[11]+0.5000000000000001*vmap0R2*wgt[5]*f[11]+0.8944271909999159*vmap[0]*vmap[1]*wgt[3]*f[11]+0.4472135954999579*vmap1R2*wgt[2]*f[11]+0.7857142857142858*vmap1R2*wgt[4]*f[9]+0.5000000000000001*vmap0R2*wgt[4]*f[9]+0.8944271909999161*vmap[0]*vmap[1]*wgt[1]*f[9]+0.44721359549995804*wgt[0]*vmap1R2*f[9]+0.8944271909999161*vmap[0]*vmap[1]*wgt[5]*f[7]+0.9*vmap1R2*wgt[3]*f[7]+0.5*vmap0R2*wgt[3]*f[7]+vmap[0]*vmap[1]*wgt[2]*f[7]+0.44721359549995804*vmap1R2*f[5]*wgt[5]+vmap[0]*vmap[1]*wgt[3]*f[5]+0.5*vmap1R2*wgt[2]*f[5]+0.5*vmap0R2*wgt[2]*f[5]+0.8944271909999159*vmap[0]*vmap[1]*f[4]*wgt[4]+0.4472135954999579*f[1]*vmap1R2*wgt[4]+0.9*vmap1R2*wgt[1]*f[4]+0.5*vmap0R2*wgt[1]*f[4]+vmap[0]*wgt[0]*vmap[1]*f[4]+vmap[0]*f[1]*vmap[1]*wgt[1]+0.5*wgt[0]*f[1]*vmap1R2+0.5*vmap0R2*wgt[0]*f[1]; 
  double tmp2perp[2]; 
  tmp2perp[0] = 0.7071067811865475*vmap[2]*wgt[5]*f[10]+0.7071067811865475*vmap[3]*wgt[4]*f[10]+0.7071067811865475*vmap[3]*wgt[5]*f[8]+0.7071067811865475*vmap[2]*wgt[4]*f[8]+0.7071067811865475*vmap[2]*wgt[3]*f[6]+0.7071067811865475*wgt[1]*vmap[3]*f[6]+0.7071067811865475*f[2]*vmap[3]*wgt[3]+0.7071067811865475*wgt[0]*f[3]*vmap[3]+0.7071067811865475*f[0]*wgt[2]*vmap[3]+0.7071067811865475*vmap[2]*wgt[2]*f[3]+0.7071067811865475*wgt[1]*f[2]*vmap[2]+0.7071067811865475*f[0]*wgt[0]*vmap[2]; 
  tmp2perp[1] = 0.7071067811865475*vmap[2]*wgt[5]*f[11]+0.7071067811865475*vmap[3]*wgt[4]*f[11]+0.7071067811865475*vmap[3]*wgt[5]*f[9]+0.7071067811865475*vmap[2]*wgt[4]*f[9]+0.7071067811865475*vmap[2]*wgt[3]*f[7]+0.7071067811865475*wgt[1]*vmap[3]*f[7]+0.7071067811865475*wgt[0]*vmap[3]*f[5]+0.7071067811865475*vmap[2]*wgt[2]*f[5]+0.7071067811865475*vmap[3]*wgt[3]*f[4]+0.7071067811865475*wgt[1]*vmap[2]*f[4]+0.7071067811865475*f[1]*wgt[2]*vmap[3]+0.7071067811865475*wgt[0]*f[1]*vmap[2]; 

  out[0] += (wgt[5]*f[10]+wgt[4]*f[8]+wgt[3]*f[6]+wgt[2]*f[3]+wgt[1]*f[2]+f[0]*wgt[0])*volFact; 
  out[1] += (1.0000000000000002*wgt[5]*f[11]+1.0000000000000002*wgt[4]*f[9]+wgt[3]*f[7]+wgt[2]*f[5]+wgt[1]*f[4]+wgt[0]*f[1])*volFact; 
  out[2] += tmp1[0]*volFact; 
  out[3] += tmp1[1]*volFact; 
  out[4] += ((1.4142135623730951*bmag[1]*tmp2perp[1]+1.4142135623730951*bmag[0]*tmp2perp[0])*volFact)/m_+tmp2par[0]*volFact; 
  out[5] += ((1.4142135623730951*bmag[0]*tmp2perp[1]+1.4142135623730951*tmp2perp[0]*bmag[1])*volFact)/m_+tmp2par[1]*volFact; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M2PERP_WEIGHTX_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[1]*dxv[2]/m_; 
  double tmp[2]; 
  tmp[0] = 1.4142135623730951*f[3]*vmap[3]+1.4142135623730951*f[0]*vmap[2]; 
  tmp[1] = 1.4142135623730951*vmap[3]*f[5]+1.4142135623730951*f[1]*vmap[2]; 
  out[0] += (2.0*(0.5*bmag[0]*tmp[1]*wgt[1]+0.5*tmp[0]*bmag[1]*wgt[1]+0.5*wgt[0]*bmag[1]*tmp[1]+0.5*bmag[0]*tmp[0]*wgt[0])*volFact)/m_; 
  out[1] += (2.0*(0.9*bmag[1]*tmp[1]*wgt[1]+0.5*bmag[0]*tmp[0]*wgt[1]+0.5*bmag[0]*wgt[0]*tmp[1]+0.5*tmp[0]*wgt[0]*bmag[1])*volFact)/m_; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M2PERP_WEIGHTV_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[1]*dxv[2]/m_; 
  double tmp[2]; 
  tmp[0] = 0.7071067811865475*vmap[2]*wgt[5]*f[10]+0.7071067811865475*vmap[3]*wgt[4]*f[10]+0.7071067811865475*vmap[3]*wgt[5]*f[8]+0.7071067811865475*vmap[2]*wgt[4]*f[8]+0.7071067811865475*vmap[2]*wgt[3]*f[6]+0.7071067811865475*wgt[1]*vmap[3]*f[6]+0.7071067811865475*f[2]*vmap[3]*wgt[3]+0.7071067811865475*wgt[0]*f[3]*vmap[3]+0.7071067811865475*f[0]*wgt[2]*vmap[3]+0.7071067811865475*vmap[2]*wgt[2]*f[3]+0.7071067811865475*wgt[1]*f[2]*vmap[2]+0.7071067811865475*f[0]*wgt[0]*vmap[2]; 
  tmp[1] = 0.7071067811865475*vmap[2]*wgt[5]*f[11]+0.7071067811865475*vmap[3]*wgt[4]*f[11]+0.7071067811865475*vmap[3]*wgt[5]*f[9]+0.7071067811865475*vmap[2]*wgt[4]*f[9]+0.7071067811865475*vmap[2]*wgt[3]*f[7]+0.7071067811865475*wgt[1]*vmap[3]*f[7]+0.7071067811865475*wgt[0]*vmap[3]*f[5]+0.7071067811865475*vmap[2]*wgt[2]*f[5]+0.7071067811865475*vmap[3]*wgt[3]*f[4]+0.7071067811865475*wgt[1]*vmap[2]*f[4]+0.7071067811865475*f[1]*wgt[2]*vmap[3]+0.7071067811865475*wgt[0]*f[1]*vmap[2]; 
  out[0] += (2.0*(0.7071067811865475*bmag[1]*tmp[1]+0.7071067811865475*bmag[0]*tmp[0])*volFact)/m_; 
  out[1] += (2.0*(0.7071067811865475*bmag[0]*tmp[1]+0.7071067811865475*tmp[0]*bmag[1])*volFact)/m_; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M0M1M2PARM2PERP_WEIGHTX_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[1]*dxv[2]/m_; 
  double tmp1[2]; 
  tmp1[0] = 1.4142135623730951*vmap[1]*f[2]+1.4142135623730951*f[0]*vmap[0]; 
  tmp1[1] = 1.4142135623730951*vmap[1]*f[4]+1.4142135623730951*vmap[0]*f[1]; 
  double tmp2par[2]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  tmp2par[0] = 0.8944271909999159*vmap1R2*f[8]+2.0*vmap[0]*vmap[1]*f[2]+f[0]*vmap1R2+f[0]*vmap0R2; 
  tmp2par[1] = 0.8944271909999161*vmap1R2*f[9]+2.0*vmap[0]*vmap[1]*f[4]+f[1]*vmap1R2+vmap0R2*f[1]; 
  double tmp2perp[2]; 
  tmp2perp[0] = 1.4142135623730951*f[3]*vmap[3]+1.4142135623730951*f[0]*vmap[2]; 
  tmp2perp[1] = 1.4142135623730951*vmap[3]*f[5]+1.4142135623730951*f[1]*vmap[2]; 

  out[0] += (1.4142135623730951*f[1]*wgt[1]+1.4142135623730951*f[0]*wgt[0])*volFact; 
  out[1] += (1.4142135623730951*f[0]*wgt[1]+1.4142135623730951*wgt[0]*f[1])*volFact; 
  out[2] += (0.7071067811865475*tmp1[1]*wgt[1]+0.7071067811865475*tmp1[0]*wgt[0])*volFact; 
  out[3] += (0.7071067811865475*tmp1[0]*wgt[1]+0.7071067811865475*wgt[0]*tmp1[1])*volFact; 
  out[4] += (0.7071067811865475*tmp2par[1]*wgt[1]+0.7071067811865475*tmp2par[0]*wgt[0])*volFact; 
  out[5] += (0.7071067811865475*tmp2par[0]*wgt[1]+0.7071067811865475*wgt[0]*tmp2par[1])*volFact; 
  out[6] += ((bmag[0]*tmp2perp[1]*wgt[1]+tmp2perp[0]*bmag[1]*wgt[1]+wgt[0]*bmag[1]*tmp2perp[1]+bmag[0]*tmp2perp[0]*wgt[0])*volFact)/m_; 
  out[7] += ((1.8*bmag[1]*tmp2perp[1]*wgt[1]+bmag[0]*tmp2perp[0]*wgt[1]+bmag[0]*wgt[0]*tmp2perp[1]+tmp2perp[0]*wgt[0]*bmag[1])*volFact)/m_; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M0M1M2PARM2PERP_WEIGHTV_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[1]*dxv[2]/m_; 
  double tmp1[2]; 
  tmp1[0] = 0.7071067811865475*vmap[0]*wgt[5]*f[10]+0.632455532033676*vmap[1]*wgt[3]*f[10]+0.7071067811865475*vmap[0]*wgt[4]*f[8]+0.6324555320336759*vmap[1]*wgt[1]*f[8]+0.632455532033676*vmap[1]*wgt[5]*f[6]+0.7071067811865475*vmap[0]*wgt[3]*f[6]+0.7071067811865475*vmap[1]*wgt[2]*f[6]+0.6324555320336759*vmap[1]*f[2]*wgt[4]+0.7071067811865475*vmap[1]*f[3]*wgt[3]+0.7071067811865475*vmap[0]*wgt[2]*f[3]+0.7071067811865475*vmap[0]*wgt[1]*f[2]+0.7071067811865475*wgt[0]*vmap[1]*f[2]+0.7071067811865475*f[0]*vmap[1]*wgt[1]+0.7071067811865475*f[0]*vmap[0]*wgt[0]; 
  tmp1[1] = 0.7071067811865475*vmap[0]*wgt[5]*f[11]+0.6324555320336759*vmap[1]*wgt[3]*f[11]+0.7071067811865475*vmap[0]*wgt[4]*f[9]+0.632455532033676*vmap[1]*wgt[1]*f[9]+0.632455532033676*vmap[1]*wgt[5]*f[7]+0.7071067811865475*vmap[0]*wgt[3]*f[7]+0.7071067811865475*vmap[1]*wgt[2]*f[7]+0.7071067811865475*vmap[1]*wgt[3]*f[5]+0.7071067811865475*vmap[0]*wgt[2]*f[5]+0.6324555320336759*vmap[1]*f[4]*wgt[4]+0.7071067811865475*vmap[0]*wgt[1]*f[4]+0.7071067811865475*wgt[0]*vmap[1]*f[4]+0.7071067811865475*f[1]*vmap[1]*wgt[1]+0.7071067811865475*vmap[0]*wgt[0]*f[1]; 
  double tmp2par[2]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  tmp2par[0] = 0.7857142857142857*vmap1R2*wgt[5]*f[10]+0.5*vmap0R2*wgt[5]*f[10]+0.8944271909999161*vmap[0]*vmap[1]*wgt[3]*f[10]+0.44721359549995804*vmap1R2*wgt[2]*f[10]+0.7857142857142857*vmap1R2*wgt[4]*f[8]+0.5*vmap0R2*wgt[4]*f[8]+0.8944271909999159*vmap[0]*vmap[1]*wgt[1]*f[8]+0.4472135954999579*wgt[0]*vmap1R2*f[8]+0.8944271909999161*vmap[0]*vmap[1]*wgt[5]*f[6]+0.9*vmap1R2*wgt[3]*f[6]+0.5*vmap0R2*wgt[3]*f[6]+vmap[0]*vmap[1]*wgt[2]*f[6]+0.44721359549995804*vmap1R2*f[3]*wgt[5]+0.8944271909999159*vmap[0]*vmap[1]*f[2]*wgt[4]+0.4472135954999579*f[0]*vmap1R2*wgt[4]+vmap[0]*vmap[1]*f[3]*wgt[3]+0.5*vmap1R2*wgt[2]*f[3]+0.5*vmap0R2*wgt[2]*f[3]+0.9*vmap1R2*wgt[1]*f[2]+0.5*vmap0R2*wgt[1]*f[2]+vmap[0]*wgt[0]*vmap[1]*f[2]+f[0]*vmap[0]*vmap[1]*wgt[1]+0.5*f[0]*wgt[0]*vmap1R2+0.5*f[0]*vmap0R2*wgt[0]; 
  tmp2par[1] = 0.7857142857142858*vmap1R2*wgt[5]*f[11]+0.5000000000000001*vmap0R2*wgt[5]*f[11]+0.8944271909999159*vmap[0]*vmap[1]*wgt[3]*f[11]+0.4472135954999579*vmap1R2*wgt[2]*f[11]+0.7857142857142858*vmap1R2*wgt[4]*f[9]+0.5000000000000001*vmap0R2*wgt[4]*f[9]+0.8944271909999161*vmap[0]*vmap[1]*wgt[1]*f[9]+0.44721359549995804*wgt[0]*vmap1R2*f[9]+0.8944271909999161*vmap[0]*vmap[1]*wgt[5]*f[7]+0.9*vmap1R2*wgt[3]*f[7]+0.5*vmap0R2*wgt[3]*f[7]+vmap[0]*vmap[1]*wgt[2]*f[7]+0.44721359549995804*vmap1R2*f[5]*wgt[5]+vmap[0]*vmap[1]*wgt[3]*f[5]+0.5*vmap1R2*wgt[2]*f[5]+0.5*vmap0R2*wgt[2]*f[5]+0.8944271909999159*vmap[0]*vmap[1]*f[4]*wgt[4]+0.4472135954999579*f[1]*vmap1R2*wgt[4]+0.9*vmap1R2*wgt[1]*f[4]+0.5*vmap0R2*wgt[1]*f[4]+vmap[0]*wgt[0]*vmap[1]*f[4]+vmap[0]*f[1]*vmap[1]*wgt[1]+0.5*wgt[0]*f[1]*vmap1R2+0.5*vmap0R2*wgt[0]*f[1]; 
  double tmp2perp[2]; 
  tmp2perp[0] = 0.7071067811865475*vmap[2]*wgt[5]*f[10]+0.7071067811865475*vmap[3]*wgt[4]*f[10]+0.7071067811865475*vmap[3]*wgt[5]*f[8]+0.7071067811865475*vmap[2]*wgt[4]*f[8]+0.7071067811865475*vmap[2]*wgt[3]*f[6]+0.7071067811865475*wgt[1]*vmap[3]*f[6]+0.7071067811865475*f[2]*vmap[3]*wgt[3]+0.7071067811865475*wgt[0]*f[3]*vmap[3]+0.7071067811865475*f[0]*wgt[2]*vmap[3]+0.7071067811865475*vmap[2]*wgt[2]*f[3]+0.7071067811865475*wgt[1]*f[2]*vmap[2]+0.7071067811865475*f[0]*wgt[0]*vmap[2]; 
  tmp2perp[1] = 0.7071067811865475*vmap[2]*wgt[5]*f[11]+0.7071067811865475*vmap[3]*wgt[4]*f[11]+0.7071067811865475*vmap[3]*wgt[5]*f[9]+0.7071067811865475*vmap[2]*wgt[4]*f[9]+0.7071067811865475*vmap[2]*wgt[3]*f[7]+0.7071067811865475*wgt[1]*vmap[3]*f[7]+0.7071067811865475*wgt[0]*vmap[3]*f[5]+0.7071067811865475*vmap[2]*wgt[2]*f[5]+0.7071067811865475*vmap[3]*wgt[3]*f[4]+0.7071067811865475*wgt[1]*vmap[2]*f[4]+0.7071067811865475*f[1]*wgt[2]*vmap[3]+0.7071067811865475*wgt[0]*f[1]*vmap[2]; 

  out[0] += (wgt[5]*f[10]+wgt[4]*f[8]+wgt[3]*f[6]+wgt[2]*f[3]+wgt[1]*f[2]+f[0]*wgt[0])*volFact; 
  out[1] += (1.0000000000000002*wgt[5]*f[11]+1.0000000000000002*wgt[4]*f[9]+wgt[3]*f[7]+wgt[2]*f[5]+wgt[1]*f[4]+wgt[0]*f[1])*volFact; 
  out[2] += tmp1[0]*volFact; 
  out[3] += tmp1[1]*volFact; 
  out[4] += tmp2par[0]*volFact; 
  out[5] += tmp2par[1]*volFact; 
  out[6] += ((1.4142135623730951*bmag[1]*tmp2perp[1]+1.4142135623730951*bmag[0]*tmp2perp[0])*volFact)/m_; 
  out[7] += ((1.4142135623730951*bmag[0]*tmp2perp[1]+1.4142135623730951*tmp2perp[0]*bmag[1])*volFact)/m_; 
} 

