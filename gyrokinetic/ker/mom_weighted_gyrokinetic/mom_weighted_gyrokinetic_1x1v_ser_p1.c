#include <gkyl_mom_weighted_gyrokinetic_kernels.h> 
GKYL_CU_DH void mom_weighted_gyrokinetic_M0_WEIGHTX_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.5*dxv[1]; 
  out[0] += (f[1]*wgt[1]+f[0]*wgt[0])*volFact; 
  out[1] += (f[0]*wgt[1]+wgt[0]*f[1])*volFact; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M0_WEIGHTV_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.5*dxv[1]; 
  out[0] += (wgt[2]*f[4]+wgt[1]*f[2]+f[0]*wgt[0])*volFact; 
  out[1] += (1.0000000000000002*wgt[2]*f[5]+wgt[1]*f[3]+wgt[0]*f[1])*volFact; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M1_WEIGHTX_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.5*dxv[1]; 
  double tmp[2]; 
  tmp[0] = vmap[1]*f[2]+f[0]*vmap[0]; 
  tmp[1] = vmap[1]*f[3]+vmap[0]*f[1]; 
  out[0] += (0.7071067811865475*tmp[1]*wgt[1]+0.7071067811865475*tmp[0]*wgt[0])*volFact; 
  out[1] += (0.7071067811865475*tmp[0]*wgt[1]+0.7071067811865475*wgt[0]*tmp[1])*volFact; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M1_WEIGHTV_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.5*dxv[1]; 
  double tmp[2]; 
  tmp[0] = 0.7071067811865475*vmap[0]*wgt[2]*f[4]+0.6324555320336759*vmap[1]*wgt[1]*f[4]+0.6324555320336759*vmap[1]*f[2]*wgt[2]+0.7071067811865475*vmap[0]*wgt[1]*f[2]+0.7071067811865475*wgt[0]*vmap[1]*f[2]+0.7071067811865475*f[0]*vmap[1]*wgt[1]+0.7071067811865475*f[0]*vmap[0]*wgt[0]; 
  tmp[1] = 0.7071067811865475*vmap[0]*wgt[2]*f[5]+0.632455532033676*vmap[1]*wgt[1]*f[5]+0.6324555320336759*vmap[1]*wgt[2]*f[3]+0.7071067811865475*vmap[0]*wgt[1]*f[3]+0.7071067811865475*wgt[0]*vmap[1]*f[3]+0.7071067811865475*f[1]*vmap[1]*wgt[1]+0.7071067811865475*vmap[0]*wgt[0]*f[1]; 
  out[0] += tmp[0]*volFact; 
  out[1] += tmp[1]*volFact; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M2_WEIGHTX_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.5*dxv[1]; 
  double tmpPar[2]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  tmpPar[0] = 0.6324555320336759*vmap1R2*f[4]+1.4142135623730951*vmap[0]*vmap[1]*f[2]+0.7071067811865475*f[0]*vmap1R2+0.7071067811865475*f[0]*vmap0R2; 
  tmpPar[1] = 0.632455532033676*vmap1R2*f[5]+1.4142135623730951*vmap[0]*vmap[1]*f[3]+0.7071067811865475*f[1]*vmap1R2+0.7071067811865475*vmap0R2*f[1]; 

  out[0] += (0.7071067811865475*tmpPar[1]*wgt[1]+0.7071067811865475*tmpPar[0]*wgt[0])*volFact; 
  out[1] += (0.7071067811865475*tmpPar[0]*wgt[1]+0.7071067811865475*wgt[0]*tmpPar[1])*volFact; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M2_WEIGHTV_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.5*dxv[1]; 
  double tmpPar[2]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  tmpPar[0] = 0.7857142857142857*vmap1R2*wgt[2]*f[4]+0.5*vmap0R2*wgt[2]*f[4]+0.8944271909999159*vmap[0]*vmap[1]*wgt[1]*f[4]+0.4472135954999579*wgt[0]*vmap1R2*f[4]+0.8944271909999159*vmap[0]*vmap[1]*f[2]*wgt[2]+0.4472135954999579*f[0]*vmap1R2*wgt[2]+0.9*vmap1R2*wgt[1]*f[2]+0.5*vmap0R2*wgt[1]*f[2]+vmap[0]*wgt[0]*vmap[1]*f[2]+f[0]*vmap[0]*vmap[1]*wgt[1]+0.5*f[0]*wgt[0]*vmap1R2+0.5*f[0]*vmap0R2*wgt[0]; 
  tmpPar[1] = 0.7857142857142858*vmap1R2*wgt[2]*f[5]+0.5000000000000001*vmap0R2*wgt[2]*f[5]+0.8944271909999161*vmap[0]*vmap[1]*wgt[1]*f[5]+0.44721359549995804*wgt[0]*vmap1R2*f[5]+0.8944271909999159*vmap[0]*vmap[1]*wgt[2]*f[3]+0.9*vmap1R2*wgt[1]*f[3]+0.5*vmap0R2*wgt[1]*f[3]+vmap[0]*wgt[0]*vmap[1]*f[3]+0.4472135954999579*f[1]*vmap1R2*wgt[2]+vmap[0]*f[1]*vmap[1]*wgt[1]+0.5*wgt[0]*f[1]*vmap1R2+0.5*vmap0R2*wgt[0]*f[1]; 

  out[0] += tmpPar[0]*volFact; 
  out[1] += tmpPar[1]*volFact; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M2PAR_WEIGHTX_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.5*dxv[1]; 
  double tmp[2]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  tmp[0] = 0.6324555320336759*vmap1R2*f[4]+1.4142135623730951*vmap[0]*vmap[1]*f[2]+0.7071067811865475*f[0]*vmap1R2+0.7071067811865475*f[0]*vmap0R2; 
  tmp[1] = 0.632455532033676*vmap1R2*f[5]+1.4142135623730951*vmap[0]*vmap[1]*f[3]+0.7071067811865475*f[1]*vmap1R2+0.7071067811865475*vmap0R2*f[1]; 

  out[0] += (0.7071067811865475*tmp[1]*wgt[1]+0.7071067811865475*tmp[0]*wgt[0])*volFact; 
  out[1] += (0.7071067811865475*tmp[0]*wgt[1]+0.7071067811865475*wgt[0]*tmp[1])*volFact; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M2PAR_WEIGHTV_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.5*dxv[1]; 
  double tmp[2]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  tmp[0] = 0.7857142857142857*vmap1R2*wgt[2]*f[4]+0.5*vmap0R2*wgt[2]*f[4]+0.8944271909999159*vmap[0]*vmap[1]*wgt[1]*f[4]+0.4472135954999579*wgt[0]*vmap1R2*f[4]+0.8944271909999159*vmap[0]*vmap[1]*f[2]*wgt[2]+0.4472135954999579*f[0]*vmap1R2*wgt[2]+0.9*vmap1R2*wgt[1]*f[2]+0.5*vmap0R2*wgt[1]*f[2]+vmap[0]*wgt[0]*vmap[1]*f[2]+f[0]*vmap[0]*vmap[1]*wgt[1]+0.5*f[0]*wgt[0]*vmap1R2+0.5*f[0]*vmap0R2*wgt[0]; 
  tmp[1] = 0.7857142857142858*vmap1R2*wgt[2]*f[5]+0.5000000000000001*vmap0R2*wgt[2]*f[5]+0.8944271909999161*vmap[0]*vmap[1]*wgt[1]*f[5]+0.44721359549995804*wgt[0]*vmap1R2*f[5]+0.8944271909999159*vmap[0]*vmap[1]*wgt[2]*f[3]+0.9*vmap1R2*wgt[1]*f[3]+0.5*vmap0R2*wgt[1]*f[3]+vmap[0]*wgt[0]*vmap[1]*f[3]+0.4472135954999579*f[1]*vmap1R2*wgt[2]+vmap[0]*f[1]*vmap[1]*wgt[1]+0.5*wgt[0]*f[1]*vmap1R2+0.5*vmap0R2*wgt[0]*f[1]; 

  out[0] += tmp[0]*volFact; 
  out[1] += tmp[1]*volFact; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M0M1M2_WEIGHTX_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.5*dxv[1]; 
  double tmp1[2]; 
  tmp1[0] = vmap[1]*f[2]+f[0]*vmap[0]; 
  tmp1[1] = vmap[1]*f[3]+vmap[0]*f[1]; 
  double tmp2par[2]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  tmp2par[0] = 0.6324555320336759*vmap1R2*f[4]+1.4142135623730951*vmap[0]*vmap[1]*f[2]+0.7071067811865475*f[0]*vmap1R2+0.7071067811865475*f[0]*vmap0R2; 
  tmp2par[1] = 0.632455532033676*vmap1R2*f[5]+1.4142135623730951*vmap[0]*vmap[1]*f[3]+0.7071067811865475*f[1]*vmap1R2+0.7071067811865475*vmap0R2*f[1]; 

  out[0] += (f[1]*wgt[1]+f[0]*wgt[0])*volFact; 
  out[1] += (f[0]*wgt[1]+wgt[0]*f[1])*volFact; 
  out[2] += (0.7071067811865475*tmp1[1]*wgt[1]+0.7071067811865475*tmp1[0]*wgt[0])*volFact; 
  out[3] += (0.7071067811865475*tmp1[0]*wgt[1]+0.7071067811865475*wgt[0]*tmp1[1])*volFact; 
  out[4] += (0.7071067811865475*tmp2par[1]*wgt[1]+0.7071067811865475*tmp2par[0]*wgt[0])*volFact; 
  out[5] += (0.7071067811865475*tmp2par[0]*wgt[1]+0.7071067811865475*wgt[0]*tmp2par[1])*volFact; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M0M1M2_WEIGHTV_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.5*dxv[1]; 
  double tmp1[2]; 
  tmp1[0] = 0.7071067811865475*vmap[0]*wgt[2]*f[4]+0.6324555320336759*vmap[1]*wgt[1]*f[4]+0.6324555320336759*vmap[1]*f[2]*wgt[2]+0.7071067811865475*vmap[0]*wgt[1]*f[2]+0.7071067811865475*wgt[0]*vmap[1]*f[2]+0.7071067811865475*f[0]*vmap[1]*wgt[1]+0.7071067811865475*f[0]*vmap[0]*wgt[0]; 
  tmp1[1] = 0.7071067811865475*vmap[0]*wgt[2]*f[5]+0.632455532033676*vmap[1]*wgt[1]*f[5]+0.6324555320336759*vmap[1]*wgt[2]*f[3]+0.7071067811865475*vmap[0]*wgt[1]*f[3]+0.7071067811865475*wgt[0]*vmap[1]*f[3]+0.7071067811865475*f[1]*vmap[1]*wgt[1]+0.7071067811865475*vmap[0]*wgt[0]*f[1]; 
  double tmp2par[2]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  tmp2par[0] = 0.7857142857142857*vmap1R2*wgt[2]*f[4]+0.5*vmap0R2*wgt[2]*f[4]+0.8944271909999159*vmap[0]*vmap[1]*wgt[1]*f[4]+0.4472135954999579*wgt[0]*vmap1R2*f[4]+0.8944271909999159*vmap[0]*vmap[1]*f[2]*wgt[2]+0.4472135954999579*f[0]*vmap1R2*wgt[2]+0.9*vmap1R2*wgt[1]*f[2]+0.5*vmap0R2*wgt[1]*f[2]+vmap[0]*wgt[0]*vmap[1]*f[2]+f[0]*vmap[0]*vmap[1]*wgt[1]+0.5*f[0]*wgt[0]*vmap1R2+0.5*f[0]*vmap0R2*wgt[0]; 
  tmp2par[1] = 0.7857142857142858*vmap1R2*wgt[2]*f[5]+0.5000000000000001*vmap0R2*wgt[2]*f[5]+0.8944271909999161*vmap[0]*vmap[1]*wgt[1]*f[5]+0.44721359549995804*wgt[0]*vmap1R2*f[5]+0.8944271909999159*vmap[0]*vmap[1]*wgt[2]*f[3]+0.9*vmap1R2*wgt[1]*f[3]+0.5*vmap0R2*wgt[1]*f[3]+vmap[0]*wgt[0]*vmap[1]*f[3]+0.4472135954999579*f[1]*vmap1R2*wgt[2]+vmap[0]*f[1]*vmap[1]*wgt[1]+0.5*wgt[0]*f[1]*vmap1R2+0.5*vmap0R2*wgt[0]*f[1]; 

  out[0] += (wgt[2]*f[4]+wgt[1]*f[2]+f[0]*wgt[0])*volFact; 
  out[1] += (1.0000000000000002*wgt[2]*f[5]+wgt[1]*f[3]+wgt[0]*f[1])*volFact; 
  out[2] += tmp1[0]*volFact; 
  out[3] += tmp1[1]*volFact; 
  out[4] += tmp2par[0]*volFact; 
  out[5] += tmp2par[1]*volFact; 
} 

