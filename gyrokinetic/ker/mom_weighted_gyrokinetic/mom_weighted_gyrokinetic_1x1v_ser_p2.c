#include <gkyl_mom_weighted_gyrokinetic_kernels.h> 
GKYL_CU_DH void mom_weighted_gyrokinetic_M0_WEIGHTX_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.5*dxv[1]; 
  out[0] += (wgt[2]*f[4]+f[1]*wgt[1]+f[0]*wgt[0])*volFact; 
  out[1] += (0.8944271909999159*wgt[1]*f[4]+0.8944271909999159*f[1]*wgt[2]+f[0]*wgt[1]+wgt[0]*f[1])*volFact; 
  out[2] += (0.6388765649999399*wgt[2]*f[4]+wgt[0]*f[4]+f[0]*wgt[2]+0.8944271909999159*f[1]*wgt[1])*volFact; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M0_WEIGHTV_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.5*dxv[1]; 
  out[0] += (wgt[2]*f[5]+wgt[1]*f[2]+f[0]*wgt[0])*volFact; 
  out[1] += (1.0000000000000002*wgt[2]*f[7]+wgt[1]*f[3]+wgt[0]*f[1])*volFact; 
  out[2] += (1.0000000000000002*wgt[1]*f[6]+wgt[0]*f[4])*volFact; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M1_WEIGHTX_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.5*dxv[1]; 
  double tmp[3]; 
  tmp[0] = vmap[1]*f[2]+f[0]*vmap[0]; 
  tmp[1] = vmap[1]*f[3]+vmap[0]*f[1]; 
  tmp[2] = 1.0000000000000002*vmap[1]*f[6]+vmap[0]*f[4]; 
  out[0] += (0.7071067811865475*tmp[2]*wgt[2]+0.7071067811865475*tmp[1]*wgt[1]+0.7071067811865475*tmp[0]*wgt[0])*volFact; 
  out[1] += (0.6324555320336759*tmp[1]*wgt[2]+0.6324555320336759*wgt[1]*tmp[2]+0.7071067811865475*tmp[0]*wgt[1]+0.7071067811865475*wgt[0]*tmp[1])*volFact; 
  out[2] += (0.45175395145262565*tmp[2]*wgt[2]+0.7071067811865475*tmp[0]*wgt[2]+0.7071067811865475*wgt[0]*tmp[2]+0.6324555320336759*tmp[1]*wgt[1])*volFact; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M1_WEIGHTV_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.5*dxv[1]; 
  double tmp[3]; 
  tmp[0] = 0.7071067811865475*vmap[0]*wgt[2]*f[5]+0.6324555320336759*vmap[1]*wgt[1]*f[5]+0.6324555320336759*vmap[1]*f[2]*wgt[2]+0.7071067811865475*vmap[0]*wgt[1]*f[2]+0.7071067811865475*wgt[0]*vmap[1]*f[2]+0.7071067811865475*f[0]*vmap[1]*wgt[1]+0.7071067811865475*f[0]*vmap[0]*wgt[0]; 
  tmp[1] = 0.7071067811865475*vmap[0]*wgt[2]*f[7]+0.632455532033676*vmap[1]*wgt[1]*f[7]+0.6324555320336759*vmap[1]*wgt[2]*f[3]+0.7071067811865475*vmap[0]*wgt[1]*f[3]+0.7071067811865475*wgt[0]*vmap[1]*f[3]+0.7071067811865475*f[1]*vmap[1]*wgt[1]+0.7071067811865475*vmap[0]*wgt[0]*f[1]; 
  tmp[2] = 0.632455532033676*vmap[1]*wgt[2]*f[6]+0.7071067811865475*vmap[0]*wgt[1]*f[6]+0.7071067811865475*wgt[0]*vmap[1]*f[6]+0.7071067811865475*vmap[1]*wgt[1]*f[4]+0.7071067811865475*vmap[0]*wgt[0]*f[4]; 
  out[0] += tmp[0]*volFact; 
  out[1] += tmp[1]*volFact; 
  out[2] += tmp[2]*volFact; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M2_WEIGHTX_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.5*dxv[1]; 
  double tmpPar[3]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  tmpPar[0] = 0.6324555320336759*vmap1R2*f[5]+1.4142135623730951*vmap[0]*vmap[1]*f[2]+0.7071067811865475*f[0]*vmap1R2+0.7071067811865475*f[0]*vmap0R2; 
  tmpPar[1] = 0.632455532033676*vmap1R2*f[7]+1.4142135623730951*vmap[0]*vmap[1]*f[3]+0.7071067811865475*f[1]*vmap1R2+0.7071067811865475*vmap0R2*f[1]; 
  tmpPar[2] = 1.4142135623730951*vmap[0]*vmap[1]*f[6]+0.7071067811865475*vmap1R2*f[4]+0.7071067811865475*vmap0R2*f[4]; 

  out[0] += (0.7071067811865475*tmpPar[2]*wgt[2]+0.7071067811865475*tmpPar[1]*wgt[1]+0.7071067811865475*tmpPar[0]*wgt[0])*volFact; 
  out[1] += (0.6324555320336759*tmpPar[1]*wgt[2]+0.6324555320336759*wgt[1]*tmpPar[2]+0.7071067811865475*tmpPar[0]*wgt[1]+0.7071067811865475*wgt[0]*tmpPar[1])*volFact; 
  out[2] += (0.45175395145262565*tmpPar[2]*wgt[2]+0.7071067811865475*tmpPar[0]*wgt[2]+0.7071067811865475*wgt[0]*tmpPar[2]+0.6324555320336759*tmpPar[1]*wgt[1])*volFact; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M2_WEIGHTV_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.5*dxv[1]; 
  double tmpPar[3]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  tmpPar[0] = 0.7857142857142857*vmap1R2*wgt[2]*f[5]+0.5*vmap0R2*wgt[2]*f[5]+0.8944271909999159*vmap[0]*vmap[1]*wgt[1]*f[5]+0.4472135954999579*wgt[0]*vmap1R2*f[5]+0.8944271909999159*vmap[0]*vmap[1]*f[2]*wgt[2]+0.4472135954999579*f[0]*vmap1R2*wgt[2]+0.9*vmap1R2*wgt[1]*f[2]+0.5*vmap0R2*wgt[1]*f[2]+vmap[0]*wgt[0]*vmap[1]*f[2]+f[0]*vmap[0]*vmap[1]*wgt[1]+0.5*f[0]*wgt[0]*vmap1R2+0.5*f[0]*vmap0R2*wgt[0]; 
  tmpPar[1] = 0.7857142857142858*vmap1R2*wgt[2]*f[7]+0.5000000000000001*vmap0R2*wgt[2]*f[7]+0.8944271909999161*vmap[0]*vmap[1]*wgt[1]*f[7]+0.44721359549995804*wgt[0]*vmap1R2*f[7]+0.8944271909999159*vmap[0]*vmap[1]*wgt[2]*f[3]+0.9*vmap1R2*wgt[1]*f[3]+0.5*vmap0R2*wgt[1]*f[3]+vmap[0]*wgt[0]*vmap[1]*f[3]+0.4472135954999579*f[1]*vmap1R2*wgt[2]+vmap[0]*f[1]*vmap[1]*wgt[1]+0.5*wgt[0]*f[1]*vmap1R2+0.5*vmap0R2*wgt[0]*f[1]; 
  tmpPar[2] = 0.8944271909999161*vmap[0]*vmap[1]*wgt[2]*f[6]+0.8999999999999998*vmap1R2*wgt[1]*f[6]+0.5000000000000001*vmap0R2*wgt[1]*f[6]+1.0000000000000002*vmap[0]*wgt[0]*vmap[1]*f[6]+0.4472135954999579*vmap1R2*wgt[2]*f[4]+vmap[0]*vmap[1]*wgt[1]*f[4]+0.5*wgt[0]*vmap1R2*f[4]+0.5*vmap0R2*wgt[0]*f[4]; 

  out[0] += tmpPar[0]*volFact; 
  out[1] += tmpPar[1]*volFact; 
  out[2] += tmpPar[2]*volFact; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M2PAR_WEIGHTX_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.5*dxv[1]; 
  double tmp[3]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  tmp[0] = 0.6324555320336759*vmap1R2*f[5]+1.4142135623730951*vmap[0]*vmap[1]*f[2]+0.7071067811865475*f[0]*vmap1R2+0.7071067811865475*f[0]*vmap0R2; 
  tmp[1] = 0.632455532033676*vmap1R2*f[7]+1.4142135623730951*vmap[0]*vmap[1]*f[3]+0.7071067811865475*f[1]*vmap1R2+0.7071067811865475*vmap0R2*f[1]; 
  tmp[2] = 1.4142135623730951*vmap[0]*vmap[1]*f[6]+0.7071067811865475*vmap1R2*f[4]+0.7071067811865475*vmap0R2*f[4]; 

  out[0] += (0.7071067811865475*tmp[2]*wgt[2]+0.7071067811865475*tmp[1]*wgt[1]+0.7071067811865475*tmp[0]*wgt[0])*volFact; 
  out[1] += (0.6324555320336759*tmp[1]*wgt[2]+0.6324555320336759*wgt[1]*tmp[2]+0.7071067811865475*tmp[0]*wgt[1]+0.7071067811865475*wgt[0]*tmp[1])*volFact; 
  out[2] += (0.45175395145262565*tmp[2]*wgt[2]+0.7071067811865475*tmp[0]*wgt[2]+0.7071067811865475*wgt[0]*tmp[2]+0.6324555320336759*tmp[1]*wgt[1])*volFact; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M2PAR_WEIGHTV_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.5*dxv[1]; 
  double tmp[3]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  tmp[0] = 0.7857142857142857*vmap1R2*wgt[2]*f[5]+0.5*vmap0R2*wgt[2]*f[5]+0.8944271909999159*vmap[0]*vmap[1]*wgt[1]*f[5]+0.4472135954999579*wgt[0]*vmap1R2*f[5]+0.8944271909999159*vmap[0]*vmap[1]*f[2]*wgt[2]+0.4472135954999579*f[0]*vmap1R2*wgt[2]+0.9*vmap1R2*wgt[1]*f[2]+0.5*vmap0R2*wgt[1]*f[2]+vmap[0]*wgt[0]*vmap[1]*f[2]+f[0]*vmap[0]*vmap[1]*wgt[1]+0.5*f[0]*wgt[0]*vmap1R2+0.5*f[0]*vmap0R2*wgt[0]; 
  tmp[1] = 0.7857142857142858*vmap1R2*wgt[2]*f[7]+0.5000000000000001*vmap0R2*wgt[2]*f[7]+0.8944271909999161*vmap[0]*vmap[1]*wgt[1]*f[7]+0.44721359549995804*wgt[0]*vmap1R2*f[7]+0.8944271909999159*vmap[0]*vmap[1]*wgt[2]*f[3]+0.9*vmap1R2*wgt[1]*f[3]+0.5*vmap0R2*wgt[1]*f[3]+vmap[0]*wgt[0]*vmap[1]*f[3]+0.4472135954999579*f[1]*vmap1R2*wgt[2]+vmap[0]*f[1]*vmap[1]*wgt[1]+0.5*wgt[0]*f[1]*vmap1R2+0.5*vmap0R2*wgt[0]*f[1]; 
  tmp[2] = 0.8944271909999161*vmap[0]*vmap[1]*wgt[2]*f[6]+0.8999999999999998*vmap1R2*wgt[1]*f[6]+0.5000000000000001*vmap0R2*wgt[1]*f[6]+1.0000000000000002*vmap[0]*wgt[0]*vmap[1]*f[6]+0.4472135954999579*vmap1R2*wgt[2]*f[4]+vmap[0]*vmap[1]*wgt[1]*f[4]+0.5*wgt[0]*vmap1R2*f[4]+0.5*vmap0R2*wgt[0]*f[4]; 

  out[0] += tmp[0]*volFact; 
  out[1] += tmp[1]*volFact; 
  out[2] += tmp[2]*volFact; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M0M1M2_WEIGHTX_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.5*dxv[1]; 
  double tmp1[3]; 
  tmp1[0] = vmap[1]*f[2]+f[0]*vmap[0]; 
  tmp1[1] = vmap[1]*f[3]+vmap[0]*f[1]; 
  tmp1[2] = 1.0000000000000002*vmap[1]*f[6]+vmap[0]*f[4]; 
  double tmp2par[3]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  tmp2par[0] = 0.6324555320336759*vmap1R2*f[5]+1.4142135623730951*vmap[0]*vmap[1]*f[2]+0.7071067811865475*f[0]*vmap1R2+0.7071067811865475*f[0]*vmap0R2; 
  tmp2par[1] = 0.632455532033676*vmap1R2*f[7]+1.4142135623730951*vmap[0]*vmap[1]*f[3]+0.7071067811865475*f[1]*vmap1R2+0.7071067811865475*vmap0R2*f[1]; 
  tmp2par[2] = 1.4142135623730951*vmap[0]*vmap[1]*f[6]+0.7071067811865475*vmap1R2*f[4]+0.7071067811865475*vmap0R2*f[4]; 

  out[0] += (wgt[2]*f[4]+f[1]*wgt[1]+f[0]*wgt[0])*volFact; 
  out[1] += (0.8944271909999159*wgt[1]*f[4]+0.8944271909999159*f[1]*wgt[2]+f[0]*wgt[1]+wgt[0]*f[1])*volFact; 
  out[2] += (0.6388765649999399*wgt[2]*f[4]+wgt[0]*f[4]+f[0]*wgt[2]+0.8944271909999159*f[1]*wgt[1])*volFact; 
  out[3] += (0.7071067811865475*tmp1[2]*wgt[2]+0.7071067811865475*tmp1[1]*wgt[1]+0.7071067811865475*tmp1[0]*wgt[0])*volFact; 
  out[4] += (0.6324555320336759*tmp1[1]*wgt[2]+0.6324555320336759*wgt[1]*tmp1[2]+0.7071067811865475*tmp1[0]*wgt[1]+0.7071067811865475*wgt[0]*tmp1[1])*volFact; 
  out[5] += (0.45175395145262565*tmp1[2]*wgt[2]+0.7071067811865475*tmp1[0]*wgt[2]+0.7071067811865475*wgt[0]*tmp1[2]+0.6324555320336759*tmp1[1]*wgt[1])*volFact; 
  out[6] += (0.7071067811865475*tmp2par[2]*wgt[2]+0.7071067811865475*tmp2par[1]*wgt[1]+0.7071067811865475*tmp2par[0]*wgt[0])*volFact; 
  out[7] += (0.6324555320336759*tmp2par[1]*wgt[2]+0.6324555320336759*wgt[1]*tmp2par[2]+0.7071067811865475*tmp2par[0]*wgt[1]+0.7071067811865475*wgt[0]*tmp2par[1])*volFact; 
  out[8] += (0.45175395145262565*tmp2par[2]*wgt[2]+0.7071067811865475*tmp2par[0]*wgt[2]+0.7071067811865475*wgt[0]*tmp2par[2]+0.6324555320336759*tmp2par[1]*wgt[1])*volFact; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M0M1M2_WEIGHTV_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.5*dxv[1]; 
  double tmp1[3]; 
  tmp1[0] = 0.7071067811865475*vmap[0]*wgt[2]*f[5]+0.6324555320336759*vmap[1]*wgt[1]*f[5]+0.6324555320336759*vmap[1]*f[2]*wgt[2]+0.7071067811865475*vmap[0]*wgt[1]*f[2]+0.7071067811865475*wgt[0]*vmap[1]*f[2]+0.7071067811865475*f[0]*vmap[1]*wgt[1]+0.7071067811865475*f[0]*vmap[0]*wgt[0]; 
  tmp1[1] = 0.7071067811865475*vmap[0]*wgt[2]*f[7]+0.632455532033676*vmap[1]*wgt[1]*f[7]+0.6324555320336759*vmap[1]*wgt[2]*f[3]+0.7071067811865475*vmap[0]*wgt[1]*f[3]+0.7071067811865475*wgt[0]*vmap[1]*f[3]+0.7071067811865475*f[1]*vmap[1]*wgt[1]+0.7071067811865475*vmap[0]*wgt[0]*f[1]; 
  tmp1[2] = 0.632455532033676*vmap[1]*wgt[2]*f[6]+0.7071067811865475*vmap[0]*wgt[1]*f[6]+0.7071067811865475*wgt[0]*vmap[1]*f[6]+0.7071067811865475*vmap[1]*wgt[1]*f[4]+0.7071067811865475*vmap[0]*wgt[0]*f[4]; 
  double tmp2par[3]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  tmp2par[0] = 0.7857142857142857*vmap1R2*wgt[2]*f[5]+0.5*vmap0R2*wgt[2]*f[5]+0.8944271909999159*vmap[0]*vmap[1]*wgt[1]*f[5]+0.4472135954999579*wgt[0]*vmap1R2*f[5]+0.8944271909999159*vmap[0]*vmap[1]*f[2]*wgt[2]+0.4472135954999579*f[0]*vmap1R2*wgt[2]+0.9*vmap1R2*wgt[1]*f[2]+0.5*vmap0R2*wgt[1]*f[2]+vmap[0]*wgt[0]*vmap[1]*f[2]+f[0]*vmap[0]*vmap[1]*wgt[1]+0.5*f[0]*wgt[0]*vmap1R2+0.5*f[0]*vmap0R2*wgt[0]; 
  tmp2par[1] = 0.7857142857142858*vmap1R2*wgt[2]*f[7]+0.5000000000000001*vmap0R2*wgt[2]*f[7]+0.8944271909999161*vmap[0]*vmap[1]*wgt[1]*f[7]+0.44721359549995804*wgt[0]*vmap1R2*f[7]+0.8944271909999159*vmap[0]*vmap[1]*wgt[2]*f[3]+0.9*vmap1R2*wgt[1]*f[3]+0.5*vmap0R2*wgt[1]*f[3]+vmap[0]*wgt[0]*vmap[1]*f[3]+0.4472135954999579*f[1]*vmap1R2*wgt[2]+vmap[0]*f[1]*vmap[1]*wgt[1]+0.5*wgt[0]*f[1]*vmap1R2+0.5*vmap0R2*wgt[0]*f[1]; 
  tmp2par[2] = 0.8944271909999161*vmap[0]*vmap[1]*wgt[2]*f[6]+0.8999999999999998*vmap1R2*wgt[1]*f[6]+0.5000000000000001*vmap0R2*wgt[1]*f[6]+1.0000000000000002*vmap[0]*wgt[0]*vmap[1]*f[6]+0.4472135954999579*vmap1R2*wgt[2]*f[4]+vmap[0]*vmap[1]*wgt[1]*f[4]+0.5*wgt[0]*vmap1R2*f[4]+0.5*vmap0R2*wgt[0]*f[4]; 

  out[0] += (wgt[2]*f[5]+wgt[1]*f[2]+f[0]*wgt[0])*volFact; 
  out[1] += (1.0000000000000002*wgt[2]*f[7]+wgt[1]*f[3]+wgt[0]*f[1])*volFact; 
  out[2] += (1.0000000000000002*wgt[1]*f[6]+wgt[0]*f[4])*volFact; 
  out[3] += tmp1[0]*volFact; 
  out[4] += tmp1[1]*volFact; 
  out[5] += tmp1[2]*volFact; 
  out[6] += tmp2par[0]*volFact; 
  out[7] += tmp2par[1]*volFact; 
  out[8] += tmp2par[2]*volFact; 
} 

