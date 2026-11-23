#include <gkyl_mom_weighted_gyrokinetic_kernels.h> 
GKYL_CU_DH void int_mom_weighted_gyrokinetic_M0_WEIGHTX_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.39269908169872414*dxv[0]*dxv[1]*dxv[2]*dxv[3]/m_; 
 

  out[0] += (2.0*wgt[3]*f[5]+2.0*f[2]*wgt[2]+2.0*f[1]*wgt[1]+2.0*f[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_weighted_gyrokinetic_M0_WEIGHTV_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.39269908169872414*dxv[0]*dxv[1]*dxv[2]*dxv[3]/m_; 
 

  out[0] += (2.0*wgt[5]*f[19]+2.0*wgt[4]*f[16]+2.0*wgt[3]*f[10]+2.0*wgt[2]*f[4]+2.0*wgt[1]*f[3]+2.0*f[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_weighted_gyrokinetic_M1_WEIGHTX_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.39269908169872414*dxv[0]*dxv[1]*dxv[2]*dxv[3]/m_; 
 

  out[0] += (1.4142135623730951*vmap[1]*wgt[3]*f[11]+1.4142135623730951*vmap[1]*wgt[2]*f[7]+1.4142135623730951*vmap[1]*wgt[1]*f[6]+1.4142135623730951*vmap[0]*wgt[3]*f[5]+1.4142135623730951*wgt[0]*vmap[1]*f[3]+1.4142135623730951*vmap[0]*f[2]*wgt[2]+1.4142135623730951*vmap[0]*f[1]*wgt[1]+1.4142135623730951*f[0]*vmap[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_weighted_gyrokinetic_M1_WEIGHTV_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.39269908169872414*dxv[0]*dxv[1]*dxv[2]*dxv[3]/m_; 
 

  out[0] += (1.4142135623730951*vmap[0]*wgt[5]*f[19]+1.264911064067352*vmap[1]*wgt[3]*f[19]+1.4142135623730951*vmap[0]*wgt[4]*f[16]+1.264911064067352*vmap[1]*wgt[1]*f[16]+1.264911064067352*vmap[1]*wgt[5]*f[10]+1.4142135623730951*vmap[0]*wgt[3]*f[10]+1.4142135623730951*vmap[1]*wgt[2]*f[10]+1.264911064067352*vmap[1]*f[3]*wgt[4]+1.4142135623730951*vmap[1]*wgt[3]*f[4]+1.4142135623730951*vmap[0]*wgt[2]*f[4]+1.4142135623730951*vmap[0]*wgt[1]*f[3]+1.4142135623730951*wgt[0]*vmap[1]*f[3]+1.4142135623730951*f[0]*vmap[1]*wgt[1]+1.4142135623730951*f[0]*vmap[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_weighted_gyrokinetic_M2PAR_WEIGHTX_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.39269908169872414*dxv[0]*dxv[1]*dxv[2]*dxv[3]/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (0.8944271909999159*vmap1R2*wgt[3]*f[20]+0.8944271909999161*vmap1R2*wgt[2]*f[18]+0.8944271909999161*vmap1R2*wgt[1]*f[17]+0.8944271909999159*wgt[0]*vmap1R2*f[16]+2.0*vmap[0]*vmap[1]*wgt[3]*f[11]+2.0*vmap[0]*vmap[1]*wgt[2]*f[7]+2.0*vmap[0]*vmap[1]*wgt[1]*f[6]+vmap1R2*wgt[3]*f[5]+vmap0R2*wgt[3]*f[5]+2.0*vmap[0]*wgt[0]*vmap[1]*f[3]+vmap1R2*f[2]*wgt[2]+vmap0R2*f[2]*wgt[2]+f[1]*vmap1R2*wgt[1]+vmap0R2*f[1]*wgt[1]+f[0]*wgt[0]*vmap1R2+f[0]*vmap0R2*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_weighted_gyrokinetic_M2PAR_WEIGHTV_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.39269908169872414*dxv[0]*dxv[1]*dxv[2]*dxv[3]/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (1.5714285714285714*vmap1R2*wgt[5]*f[19]+vmap0R2*wgt[5]*f[19]+1.7888543819998322*vmap[0]*vmap[1]*wgt[3]*f[19]+0.8944271909999161*vmap1R2*wgt[2]*f[19]+1.5714285714285714*vmap1R2*wgt[4]*f[16]+vmap0R2*wgt[4]*f[16]+1.7888543819998317*vmap[0]*vmap[1]*wgt[1]*f[16]+0.8944271909999159*wgt[0]*vmap1R2*f[16]+1.7888543819998322*vmap[0]*vmap[1]*wgt[5]*f[10]+1.8*vmap1R2*wgt[3]*f[10]+vmap0R2*wgt[3]*f[10]+2.0*vmap[0]*vmap[1]*wgt[2]*f[10]+0.8944271909999161*vmap1R2*f[4]*wgt[5]+1.7888543819998317*vmap[0]*vmap[1]*f[3]*wgt[4]+0.8944271909999159*f[0]*vmap1R2*wgt[4]+2.0*vmap[0]*vmap[1]*wgt[3]*f[4]+vmap1R2*wgt[2]*f[4]+vmap0R2*wgt[2]*f[4]+1.8*vmap1R2*wgt[1]*f[3]+vmap0R2*wgt[1]*f[3]+2.0*vmap[0]*wgt[0]*vmap[1]*f[3]+2.0*f[0]*vmap[0]*vmap[1]*wgt[1]+f[0]*wgt[0]*vmap1R2+f[0]*vmap0R2*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_weighted_gyrokinetic_M2_WEIGHTX_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.39269908169872414*dxv[0]*dxv[1]*dxv[2]*dxv[3]/m_; 
 
  double tmp[4]; 
  tmp[0] = (2.8284271247461907*vmap[3]*f[4])/m_+(2.8284271247461907*f[0]*vmap[2])/m_; 
  tmp[1] = (2.8284271247461907*vmap[3]*f[8])/m_+(2.8284271247461907*f[1]*vmap[2])/m_; 
  tmp[2] = (2.8284271247461907*vmap[3]*f[9])/m_+(2.8284271247461907*f[2]*vmap[2])/m_; 
  tmp[3] = (2.8284271247461907*vmap[3]*f[12])/m_+(2.8284271247461907*vmap[2]*f[5])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (0.8944271909999159*vmap1R2*wgt[3]*f[20]+0.8944271909999161*vmap1R2*wgt[2]*f[18]+0.8944271909999161*vmap1R2*wgt[1]*f[17]+0.8944271909999159*wgt[0]*vmap1R2*f[16]+2.0*vmap[0]*vmap[1]*wgt[3]*f[11]+2.0*vmap[0]*vmap[1]*wgt[2]*f[7]+2.0*vmap[0]*vmap[1]*wgt[1]*f[6]+vmap1R2*wgt[3]*f[5]+vmap0R2*wgt[3]*f[5]+0.5*bmag[0]*tmp[3]*wgt[3]+0.5*tmp[0]*bmag[3]*wgt[3]+0.5*bmag[1]*tmp[2]*wgt[3]+0.5*tmp[1]*bmag[2]*wgt[3]+0.5*wgt[0]*bmag[3]*tmp[3]+0.5*bmag[1]*wgt[2]*tmp[3]+0.5*wgt[1]*bmag[2]*tmp[3]+2.0*vmap[0]*wgt[0]*vmap[1]*f[3]+0.5*tmp[1]*wgt[2]*bmag[3]+0.5*wgt[1]*tmp[2]*bmag[3]+0.5*bmag[0]*tmp[2]*wgt[2]+vmap1R2*f[2]*wgt[2]+vmap0R2*f[2]*wgt[2]+0.5*tmp[0]*bmag[2]*wgt[2]+0.5*wgt[0]*bmag[2]*tmp[2]+f[1]*vmap1R2*wgt[1]+0.5*bmag[0]*tmp[1]*wgt[1]+vmap0R2*f[1]*wgt[1]+0.5*tmp[0]*bmag[1]*wgt[1]+f[0]*wgt[0]*vmap1R2+0.5*wgt[0]*bmag[1]*tmp[1]+f[0]*vmap0R2*wgt[0]+0.5*bmag[0]*tmp[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_weighted_gyrokinetic_M2_WEIGHTV_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.39269908169872414*dxv[0]*dxv[1]*dxv[2]*dxv[3]/m_; 
 
  double tmp[4]; 
  tmp[0] = (1.4142135623730951*vmap[2]*wgt[5]*f[19])/m_+(1.4142135623730951*vmap[3]*wgt[4]*f[19])/m_+(1.4142135623730951*vmap[3]*wgt[5]*f[16])/m_+(1.4142135623730951*vmap[2]*wgt[4]*f[16])/m_+(1.4142135623730951*vmap[2]*wgt[3]*f[10])/m_+(1.4142135623730951*wgt[1]*vmap[3]*f[10])/m_+(1.4142135623730951*wgt[0]*vmap[3]*f[4])/m_+(1.4142135623730951*vmap[2]*wgt[2]*f[4])/m_+(1.4142135623730951*f[3]*vmap[3]*wgt[3])/m_+(1.4142135623730951*f[0]*wgt[2]*vmap[3])/m_+(1.4142135623730951*wgt[1]*vmap[2]*f[3])/m_+(1.4142135623730951*f[0]*wgt[0]*vmap[2])/m_; 
  tmp[1] = (1.4142135623730951*vmap[2]*wgt[5]*f[21])/m_+(1.4142135623730951*vmap[3]*wgt[4]*f[21])/m_+(1.4142135623730951*vmap[3]*wgt[5]*f[17])/m_+(1.4142135623730951*vmap[2]*wgt[4]*f[17])/m_+(1.4142135623730951*vmap[2]*wgt[3]*f[13])/m_+(1.4142135623730951*wgt[1]*vmap[3]*f[13])/m_+(1.4142135623730951*wgt[0]*vmap[3]*f[8])/m_+(1.4142135623730951*vmap[2]*wgt[2]*f[8])/m_+(1.4142135623730951*vmap[3]*wgt[3]*f[6])/m_+(1.4142135623730951*wgt[1]*vmap[2]*f[6])/m_+(1.4142135623730951*f[1]*wgt[2]*vmap[3])/m_+(1.4142135623730951*wgt[0]*f[1]*vmap[2])/m_; 
  tmp[2] = (1.4142135623730951*vmap[2]*wgt[5]*f[22])/m_+(1.4142135623730951*vmap[3]*wgt[4]*f[22])/m_+(1.4142135623730951*vmap[3]*wgt[5]*f[18])/m_+(1.4142135623730951*vmap[2]*wgt[4]*f[18])/m_+(1.4142135623730951*vmap[2]*wgt[3]*f[14])/m_+(1.4142135623730951*wgt[1]*vmap[3]*f[14])/m_+(1.4142135623730951*wgt[0]*vmap[3]*f[9])/m_+(1.4142135623730951*vmap[2]*wgt[2]*f[9])/m_+(1.4142135623730951*vmap[3]*wgt[3]*f[7])/m_+(1.4142135623730951*wgt[1]*vmap[2]*f[7])/m_+(1.4142135623730951*f[2]*wgt[2]*vmap[3])/m_+(1.4142135623730951*wgt[0]*f[2]*vmap[2])/m_; 
  tmp[3] = (1.4142135623730951*vmap[2]*wgt[5]*f[23])/m_+(1.4142135623730951*vmap[3]*wgt[4]*f[23])/m_+(1.4142135623730951*vmap[3]*wgt[5]*f[20])/m_+(1.4142135623730951*vmap[2]*wgt[4]*f[20])/m_+(1.4142135623730951*vmap[2]*wgt[3]*f[15])/m_+(1.4142135623730951*wgt[1]*vmap[3]*f[15])/m_+(1.4142135623730951*wgt[0]*vmap[3]*f[12])/m_+(1.4142135623730951*vmap[2]*wgt[2]*f[12])/m_+(1.4142135623730951*vmap[3]*wgt[3]*f[11])/m_+(1.4142135623730951*wgt[1]*vmap[2]*f[11])/m_+(1.4142135623730951*wgt[2]*vmap[3]*f[5])/m_+(1.4142135623730951*wgt[0]*vmap[2]*f[5])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (1.5714285714285714*vmap1R2*wgt[5]*f[19]+vmap0R2*wgt[5]*f[19]+1.7888543819998322*vmap[0]*vmap[1]*wgt[3]*f[19]+0.8944271909999161*vmap1R2*wgt[2]*f[19]+1.5714285714285714*vmap1R2*wgt[4]*f[16]+vmap0R2*wgt[4]*f[16]+1.7888543819998317*vmap[0]*vmap[1]*wgt[1]*f[16]+0.8944271909999159*wgt[0]*vmap1R2*f[16]+1.7888543819998322*vmap[0]*vmap[1]*wgt[5]*f[10]+1.8*vmap1R2*wgt[3]*f[10]+vmap0R2*wgt[3]*f[10]+2.0*vmap[0]*vmap[1]*wgt[2]*f[10]+0.8944271909999161*vmap1R2*f[4]*wgt[5]+1.7888543819998317*vmap[0]*vmap[1]*f[3]*wgt[4]+0.8944271909999159*f[0]*vmap1R2*wgt[4]+2.0*vmap[0]*vmap[1]*wgt[3]*f[4]+vmap1R2*wgt[2]*f[4]+vmap0R2*wgt[2]*f[4]+bmag[3]*tmp[3]+1.8*vmap1R2*wgt[1]*f[3]+vmap0R2*wgt[1]*f[3]+2.0*vmap[0]*wgt[0]*vmap[1]*f[3]+bmag[2]*tmp[2]+2.0*f[0]*vmap[0]*vmap[1]*wgt[1]+f[0]*wgt[0]*vmap1R2+bmag[1]*tmp[1]+f[0]*vmap0R2*wgt[0]+bmag[0]*tmp[0])*volFact; 
} 

GKYL_CU_DH void int_mom_weighted_gyrokinetic_M0M1M2_WEIGHTX_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.39269908169872414*dxv[0]*dxv[1]*dxv[2]*dxv[3]/m_; 
 
  double tmp[4]; 
  tmp[0] = (2.8284271247461907*vmap[3]*f[4])/m_+(2.8284271247461907*f[0]*vmap[2])/m_; 
  tmp[1] = (2.8284271247461907*vmap[3]*f[8])/m_+(2.8284271247461907*f[1]*vmap[2])/m_; 
  tmp[2] = (2.8284271247461907*vmap[3]*f[9])/m_+(2.8284271247461907*f[2]*vmap[2])/m_; 
  tmp[3] = (2.8284271247461907*vmap[3]*f[12])/m_+(2.8284271247461907*vmap[2]*f[5])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (2.0*wgt[3]*f[5]+2.0*f[2]*wgt[2]+2.0*f[1]*wgt[1]+2.0*f[0]*wgt[0])*volFact; 
  out[1] += (1.4142135623730951*vmap[1]*wgt[3]*f[11]+1.4142135623730951*vmap[1]*wgt[2]*f[7]+1.4142135623730951*vmap[1]*wgt[1]*f[6]+1.4142135623730951*vmap[0]*wgt[3]*f[5]+1.4142135623730951*wgt[0]*vmap[1]*f[3]+1.4142135623730951*vmap[0]*f[2]*wgt[2]+1.4142135623730951*vmap[0]*f[1]*wgt[1]+1.4142135623730951*f[0]*vmap[0]*wgt[0])*volFact; 
  out[2] += (0.8944271909999159*vmap1R2*wgt[3]*f[20]+0.8944271909999161*vmap1R2*wgt[2]*f[18]+0.8944271909999161*vmap1R2*wgt[1]*f[17]+0.8944271909999159*wgt[0]*vmap1R2*f[16]+2.0*vmap[0]*vmap[1]*wgt[3]*f[11]+2.0*vmap[0]*vmap[1]*wgt[2]*f[7]+2.0*vmap[0]*vmap[1]*wgt[1]*f[6]+vmap1R2*wgt[3]*f[5]+vmap0R2*wgt[3]*f[5]+0.5*bmag[0]*tmp[3]*wgt[3]+0.5*tmp[0]*bmag[3]*wgt[3]+0.5*bmag[1]*tmp[2]*wgt[3]+0.5*tmp[1]*bmag[2]*wgt[3]+0.5*wgt[0]*bmag[3]*tmp[3]+0.5*bmag[1]*wgt[2]*tmp[3]+0.5*wgt[1]*bmag[2]*tmp[3]+2.0*vmap[0]*wgt[0]*vmap[1]*f[3]+0.5*tmp[1]*wgt[2]*bmag[3]+0.5*wgt[1]*tmp[2]*bmag[3]+0.5*bmag[0]*tmp[2]*wgt[2]+vmap1R2*f[2]*wgt[2]+vmap0R2*f[2]*wgt[2]+0.5*tmp[0]*bmag[2]*wgt[2]+0.5*wgt[0]*bmag[2]*tmp[2]+f[1]*vmap1R2*wgt[1]+0.5*bmag[0]*tmp[1]*wgt[1]+vmap0R2*f[1]*wgt[1]+0.5*tmp[0]*bmag[1]*wgt[1]+f[0]*wgt[0]*vmap1R2+0.5*wgt[0]*bmag[1]*tmp[1]+f[0]*vmap0R2*wgt[0]+0.5*bmag[0]*tmp[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_weighted_gyrokinetic_M0M1M2_WEIGHTV_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.39269908169872414*dxv[0]*dxv[1]*dxv[2]*dxv[3]/m_; 
 
  double tmp[4]; 
  tmp[0] = (1.4142135623730951*vmap[2]*wgt[5]*f[19])/m_+(1.4142135623730951*vmap[3]*wgt[4]*f[19])/m_+(1.4142135623730951*vmap[3]*wgt[5]*f[16])/m_+(1.4142135623730951*vmap[2]*wgt[4]*f[16])/m_+(1.4142135623730951*vmap[2]*wgt[3]*f[10])/m_+(1.4142135623730951*wgt[1]*vmap[3]*f[10])/m_+(1.4142135623730951*wgt[0]*vmap[3]*f[4])/m_+(1.4142135623730951*vmap[2]*wgt[2]*f[4])/m_+(1.4142135623730951*f[3]*vmap[3]*wgt[3])/m_+(1.4142135623730951*f[0]*wgt[2]*vmap[3])/m_+(1.4142135623730951*wgt[1]*vmap[2]*f[3])/m_+(1.4142135623730951*f[0]*wgt[0]*vmap[2])/m_; 
  tmp[1] = (1.4142135623730951*vmap[2]*wgt[5]*f[21])/m_+(1.4142135623730951*vmap[3]*wgt[4]*f[21])/m_+(1.4142135623730951*vmap[3]*wgt[5]*f[17])/m_+(1.4142135623730951*vmap[2]*wgt[4]*f[17])/m_+(1.4142135623730951*vmap[2]*wgt[3]*f[13])/m_+(1.4142135623730951*wgt[1]*vmap[3]*f[13])/m_+(1.4142135623730951*wgt[0]*vmap[3]*f[8])/m_+(1.4142135623730951*vmap[2]*wgt[2]*f[8])/m_+(1.4142135623730951*vmap[3]*wgt[3]*f[6])/m_+(1.4142135623730951*wgt[1]*vmap[2]*f[6])/m_+(1.4142135623730951*f[1]*wgt[2]*vmap[3])/m_+(1.4142135623730951*wgt[0]*f[1]*vmap[2])/m_; 
  tmp[2] = (1.4142135623730951*vmap[2]*wgt[5]*f[22])/m_+(1.4142135623730951*vmap[3]*wgt[4]*f[22])/m_+(1.4142135623730951*vmap[3]*wgt[5]*f[18])/m_+(1.4142135623730951*vmap[2]*wgt[4]*f[18])/m_+(1.4142135623730951*vmap[2]*wgt[3]*f[14])/m_+(1.4142135623730951*wgt[1]*vmap[3]*f[14])/m_+(1.4142135623730951*wgt[0]*vmap[3]*f[9])/m_+(1.4142135623730951*vmap[2]*wgt[2]*f[9])/m_+(1.4142135623730951*vmap[3]*wgt[3]*f[7])/m_+(1.4142135623730951*wgt[1]*vmap[2]*f[7])/m_+(1.4142135623730951*f[2]*wgt[2]*vmap[3])/m_+(1.4142135623730951*wgt[0]*f[2]*vmap[2])/m_; 
  tmp[3] = (1.4142135623730951*vmap[2]*wgt[5]*f[23])/m_+(1.4142135623730951*vmap[3]*wgt[4]*f[23])/m_+(1.4142135623730951*vmap[3]*wgt[5]*f[20])/m_+(1.4142135623730951*vmap[2]*wgt[4]*f[20])/m_+(1.4142135623730951*vmap[2]*wgt[3]*f[15])/m_+(1.4142135623730951*wgt[1]*vmap[3]*f[15])/m_+(1.4142135623730951*wgt[0]*vmap[3]*f[12])/m_+(1.4142135623730951*vmap[2]*wgt[2]*f[12])/m_+(1.4142135623730951*vmap[3]*wgt[3]*f[11])/m_+(1.4142135623730951*wgt[1]*vmap[2]*f[11])/m_+(1.4142135623730951*wgt[2]*vmap[3]*f[5])/m_+(1.4142135623730951*wgt[0]*vmap[2]*f[5])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (2.0*wgt[5]*f[19]+2.0*wgt[4]*f[16]+2.0*wgt[3]*f[10]+2.0*wgt[2]*f[4]+2.0*wgt[1]*f[3]+2.0*f[0]*wgt[0])*volFact; 
  out[1] += (1.4142135623730951*vmap[0]*wgt[5]*f[19]+1.264911064067352*vmap[1]*wgt[3]*f[19]+1.4142135623730951*vmap[0]*wgt[4]*f[16]+1.264911064067352*vmap[1]*wgt[1]*f[16]+1.264911064067352*vmap[1]*wgt[5]*f[10]+1.4142135623730951*vmap[0]*wgt[3]*f[10]+1.4142135623730951*vmap[1]*wgt[2]*f[10]+1.264911064067352*vmap[1]*f[3]*wgt[4]+1.4142135623730951*vmap[1]*wgt[3]*f[4]+1.4142135623730951*vmap[0]*wgt[2]*f[4]+1.4142135623730951*vmap[0]*wgt[1]*f[3]+1.4142135623730951*wgt[0]*vmap[1]*f[3]+1.4142135623730951*f[0]*vmap[1]*wgt[1]+1.4142135623730951*f[0]*vmap[0]*wgt[0])*volFact; 
  out[2] += (1.5714285714285714*vmap1R2*wgt[5]*f[19]+vmap0R2*wgt[5]*f[19]+1.7888543819998322*vmap[0]*vmap[1]*wgt[3]*f[19]+0.8944271909999161*vmap1R2*wgt[2]*f[19]+1.5714285714285714*vmap1R2*wgt[4]*f[16]+vmap0R2*wgt[4]*f[16]+1.7888543819998317*vmap[0]*vmap[1]*wgt[1]*f[16]+0.8944271909999159*wgt[0]*vmap1R2*f[16]+1.7888543819998322*vmap[0]*vmap[1]*wgt[5]*f[10]+1.8*vmap1R2*wgt[3]*f[10]+vmap0R2*wgt[3]*f[10]+2.0*vmap[0]*vmap[1]*wgt[2]*f[10]+0.8944271909999161*vmap1R2*f[4]*wgt[5]+1.7888543819998317*vmap[0]*vmap[1]*f[3]*wgt[4]+0.8944271909999159*f[0]*vmap1R2*wgt[4]+2.0*vmap[0]*vmap[1]*wgt[3]*f[4]+vmap1R2*wgt[2]*f[4]+vmap0R2*wgt[2]*f[4]+bmag[3]*tmp[3]+1.8*vmap1R2*wgt[1]*f[3]+vmap0R2*wgt[1]*f[3]+2.0*vmap[0]*wgt[0]*vmap[1]*f[3]+bmag[2]*tmp[2]+2.0*f[0]*vmap[0]*vmap[1]*wgt[1]+f[0]*wgt[0]*vmap1R2+bmag[1]*tmp[1]+f[0]*vmap0R2*wgt[0]+bmag[0]*tmp[0])*volFact; 
} 

GKYL_CU_DH void int_mom_weighted_gyrokinetic_M2PERP_WEIGHTX_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.39269908169872414*dxv[0]*dxv[1]*dxv[2]*dxv[3]/m_; 
 
  double tmp[4]; 
  tmp[0] = (2.8284271247461907*vmap[3]*f[4])/m_+(2.8284271247461907*f[0]*vmap[2])/m_; 
  tmp[1] = (2.8284271247461907*vmap[3]*f[8])/m_+(2.8284271247461907*f[1]*vmap[2])/m_; 
  tmp[2] = (2.8284271247461907*vmap[3]*f[9])/m_+(2.8284271247461907*f[2]*vmap[2])/m_; 
  tmp[3] = (2.8284271247461907*vmap[3]*f[12])/m_+(2.8284271247461907*vmap[2]*f[5])/m_; 
 

  out[0] += (0.5*bmag[0]*tmp[3]*wgt[3]+0.5*tmp[0]*bmag[3]*wgt[3]+0.5*bmag[1]*tmp[2]*wgt[3]+0.5*tmp[1]*bmag[2]*wgt[3]+0.5*wgt[0]*bmag[3]*tmp[3]+0.5*bmag[1]*wgt[2]*tmp[3]+0.5*wgt[1]*bmag[2]*tmp[3]+0.5*tmp[1]*wgt[2]*bmag[3]+0.5*wgt[1]*tmp[2]*bmag[3]+0.5*bmag[0]*tmp[2]*wgt[2]+0.5*tmp[0]*bmag[2]*wgt[2]+0.5*wgt[0]*bmag[2]*tmp[2]+0.5*bmag[0]*tmp[1]*wgt[1]+0.5*tmp[0]*bmag[1]*wgt[1]+0.5*wgt[0]*bmag[1]*tmp[1]+0.5*bmag[0]*tmp[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_weighted_gyrokinetic_M2PERP_WEIGHTV_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.39269908169872414*dxv[0]*dxv[1]*dxv[2]*dxv[3]/m_; 
 
  double tmp[4]; 
  tmp[0] = (1.4142135623730951*vmap[2]*wgt[5]*f[19])/m_+(1.4142135623730951*vmap[3]*wgt[4]*f[19])/m_+(1.4142135623730951*vmap[3]*wgt[5]*f[16])/m_+(1.4142135623730951*vmap[2]*wgt[4]*f[16])/m_+(1.4142135623730951*vmap[2]*wgt[3]*f[10])/m_+(1.4142135623730951*wgt[1]*vmap[3]*f[10])/m_+(1.4142135623730951*wgt[0]*vmap[3]*f[4])/m_+(1.4142135623730951*vmap[2]*wgt[2]*f[4])/m_+(1.4142135623730951*f[3]*vmap[3]*wgt[3])/m_+(1.4142135623730951*f[0]*wgt[2]*vmap[3])/m_+(1.4142135623730951*wgt[1]*vmap[2]*f[3])/m_+(1.4142135623730951*f[0]*wgt[0]*vmap[2])/m_; 
  tmp[1] = (1.4142135623730951*vmap[2]*wgt[5]*f[21])/m_+(1.4142135623730951*vmap[3]*wgt[4]*f[21])/m_+(1.4142135623730951*vmap[3]*wgt[5]*f[17])/m_+(1.4142135623730951*vmap[2]*wgt[4]*f[17])/m_+(1.4142135623730951*vmap[2]*wgt[3]*f[13])/m_+(1.4142135623730951*wgt[1]*vmap[3]*f[13])/m_+(1.4142135623730951*wgt[0]*vmap[3]*f[8])/m_+(1.4142135623730951*vmap[2]*wgt[2]*f[8])/m_+(1.4142135623730951*vmap[3]*wgt[3]*f[6])/m_+(1.4142135623730951*wgt[1]*vmap[2]*f[6])/m_+(1.4142135623730951*f[1]*wgt[2]*vmap[3])/m_+(1.4142135623730951*wgt[0]*f[1]*vmap[2])/m_; 
  tmp[2] = (1.4142135623730951*vmap[2]*wgt[5]*f[22])/m_+(1.4142135623730951*vmap[3]*wgt[4]*f[22])/m_+(1.4142135623730951*vmap[3]*wgt[5]*f[18])/m_+(1.4142135623730951*vmap[2]*wgt[4]*f[18])/m_+(1.4142135623730951*vmap[2]*wgt[3]*f[14])/m_+(1.4142135623730951*wgt[1]*vmap[3]*f[14])/m_+(1.4142135623730951*wgt[0]*vmap[3]*f[9])/m_+(1.4142135623730951*vmap[2]*wgt[2]*f[9])/m_+(1.4142135623730951*vmap[3]*wgt[3]*f[7])/m_+(1.4142135623730951*wgt[1]*vmap[2]*f[7])/m_+(1.4142135623730951*f[2]*wgt[2]*vmap[3])/m_+(1.4142135623730951*wgt[0]*f[2]*vmap[2])/m_; 
  tmp[3] = (1.4142135623730951*vmap[2]*wgt[5]*f[23])/m_+(1.4142135623730951*vmap[3]*wgt[4]*f[23])/m_+(1.4142135623730951*vmap[3]*wgt[5]*f[20])/m_+(1.4142135623730951*vmap[2]*wgt[4]*f[20])/m_+(1.4142135623730951*vmap[2]*wgt[3]*f[15])/m_+(1.4142135623730951*wgt[1]*vmap[3]*f[15])/m_+(1.4142135623730951*wgt[0]*vmap[3]*f[12])/m_+(1.4142135623730951*vmap[2]*wgt[2]*f[12])/m_+(1.4142135623730951*vmap[3]*wgt[3]*f[11])/m_+(1.4142135623730951*wgt[1]*vmap[2]*f[11])/m_+(1.4142135623730951*wgt[2]*vmap[3]*f[5])/m_+(1.4142135623730951*wgt[0]*vmap[2]*f[5])/m_; 
 

  out[0] += (bmag[3]*tmp[3]+bmag[2]*tmp[2]+bmag[1]*tmp[1]+bmag[0]*tmp[0])*volFact; 
} 

GKYL_CU_DH void int_mom_weighted_gyrokinetic_M0M1M2PARM2PERP_WEIGHTX_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.39269908169872414*dxv[0]*dxv[1]*dxv[2]*dxv[3]/m_; 
 
  double tmp[4]; 
  tmp[0] = (2.8284271247461907*vmap[3]*f[4])/m_+(2.8284271247461907*f[0]*vmap[2])/m_; 
  tmp[1] = (2.8284271247461907*vmap[3]*f[8])/m_+(2.8284271247461907*f[1]*vmap[2])/m_; 
  tmp[2] = (2.8284271247461907*vmap[3]*f[9])/m_+(2.8284271247461907*f[2]*vmap[2])/m_; 
  tmp[3] = (2.8284271247461907*vmap[3]*f[12])/m_+(2.8284271247461907*vmap[2]*f[5])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (2.0*wgt[3]*f[5]+2.0*f[2]*wgt[2]+2.0*f[1]*wgt[1]+2.0*f[0]*wgt[0])*volFact; 
  out[1] += (1.4142135623730951*vmap[1]*wgt[3]*f[11]+1.4142135623730951*vmap[1]*wgt[2]*f[7]+1.4142135623730951*vmap[1]*wgt[1]*f[6]+1.4142135623730951*vmap[0]*wgt[3]*f[5]+1.4142135623730951*wgt[0]*vmap[1]*f[3]+1.4142135623730951*vmap[0]*f[2]*wgt[2]+1.4142135623730951*vmap[0]*f[1]*wgt[1]+1.4142135623730951*f[0]*vmap[0]*wgt[0])*volFact; 
  out[2] += (0.8944271909999159*vmap1R2*wgt[3]*f[20]+0.8944271909999161*vmap1R2*wgt[2]*f[18]+0.8944271909999161*vmap1R2*wgt[1]*f[17]+0.8944271909999159*wgt[0]*vmap1R2*f[16]+2.0*vmap[0]*vmap[1]*wgt[3]*f[11]+2.0*vmap[0]*vmap[1]*wgt[2]*f[7]+2.0*vmap[0]*vmap[1]*wgt[1]*f[6]+vmap1R2*wgt[3]*f[5]+vmap0R2*wgt[3]*f[5]+2.0*vmap[0]*wgt[0]*vmap[1]*f[3]+vmap1R2*f[2]*wgt[2]+vmap0R2*f[2]*wgt[2]+f[1]*vmap1R2*wgt[1]+vmap0R2*f[1]*wgt[1]+f[0]*wgt[0]*vmap1R2+f[0]*vmap0R2*wgt[0])*volFact; 
  out[3] += (0.5*bmag[0]*tmp[3]*wgt[3]+0.5*tmp[0]*bmag[3]*wgt[3]+0.5*bmag[1]*tmp[2]*wgt[3]+0.5*tmp[1]*bmag[2]*wgt[3]+0.5*wgt[0]*bmag[3]*tmp[3]+0.5*bmag[1]*wgt[2]*tmp[3]+0.5*wgt[1]*bmag[2]*tmp[3]+0.5*tmp[1]*wgt[2]*bmag[3]+0.5*wgt[1]*tmp[2]*bmag[3]+0.5*bmag[0]*tmp[2]*wgt[2]+0.5*tmp[0]*bmag[2]*wgt[2]+0.5*wgt[0]*bmag[2]*tmp[2]+0.5*bmag[0]*tmp[1]*wgt[1]+0.5*tmp[0]*bmag[1]*wgt[1]+0.5*wgt[0]*bmag[1]*tmp[1]+0.5*bmag[0]*tmp[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_weighted_gyrokinetic_M0M1M2PARM2PERP_WEIGHTV_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.39269908169872414*dxv[0]*dxv[1]*dxv[2]*dxv[3]/m_; 
 
  double tmp[4]; 
  tmp[0] = (1.4142135623730951*vmap[2]*wgt[5]*f[19])/m_+(1.4142135623730951*vmap[3]*wgt[4]*f[19])/m_+(1.4142135623730951*vmap[3]*wgt[5]*f[16])/m_+(1.4142135623730951*vmap[2]*wgt[4]*f[16])/m_+(1.4142135623730951*vmap[2]*wgt[3]*f[10])/m_+(1.4142135623730951*wgt[1]*vmap[3]*f[10])/m_+(1.4142135623730951*wgt[0]*vmap[3]*f[4])/m_+(1.4142135623730951*vmap[2]*wgt[2]*f[4])/m_+(1.4142135623730951*f[3]*vmap[3]*wgt[3])/m_+(1.4142135623730951*f[0]*wgt[2]*vmap[3])/m_+(1.4142135623730951*wgt[1]*vmap[2]*f[3])/m_+(1.4142135623730951*f[0]*wgt[0]*vmap[2])/m_; 
  tmp[1] = (1.4142135623730951*vmap[2]*wgt[5]*f[21])/m_+(1.4142135623730951*vmap[3]*wgt[4]*f[21])/m_+(1.4142135623730951*vmap[3]*wgt[5]*f[17])/m_+(1.4142135623730951*vmap[2]*wgt[4]*f[17])/m_+(1.4142135623730951*vmap[2]*wgt[3]*f[13])/m_+(1.4142135623730951*wgt[1]*vmap[3]*f[13])/m_+(1.4142135623730951*wgt[0]*vmap[3]*f[8])/m_+(1.4142135623730951*vmap[2]*wgt[2]*f[8])/m_+(1.4142135623730951*vmap[3]*wgt[3]*f[6])/m_+(1.4142135623730951*wgt[1]*vmap[2]*f[6])/m_+(1.4142135623730951*f[1]*wgt[2]*vmap[3])/m_+(1.4142135623730951*wgt[0]*f[1]*vmap[2])/m_; 
  tmp[2] = (1.4142135623730951*vmap[2]*wgt[5]*f[22])/m_+(1.4142135623730951*vmap[3]*wgt[4]*f[22])/m_+(1.4142135623730951*vmap[3]*wgt[5]*f[18])/m_+(1.4142135623730951*vmap[2]*wgt[4]*f[18])/m_+(1.4142135623730951*vmap[2]*wgt[3]*f[14])/m_+(1.4142135623730951*wgt[1]*vmap[3]*f[14])/m_+(1.4142135623730951*wgt[0]*vmap[3]*f[9])/m_+(1.4142135623730951*vmap[2]*wgt[2]*f[9])/m_+(1.4142135623730951*vmap[3]*wgt[3]*f[7])/m_+(1.4142135623730951*wgt[1]*vmap[2]*f[7])/m_+(1.4142135623730951*f[2]*wgt[2]*vmap[3])/m_+(1.4142135623730951*wgt[0]*f[2]*vmap[2])/m_; 
  tmp[3] = (1.4142135623730951*vmap[2]*wgt[5]*f[23])/m_+(1.4142135623730951*vmap[3]*wgt[4]*f[23])/m_+(1.4142135623730951*vmap[3]*wgt[5]*f[20])/m_+(1.4142135623730951*vmap[2]*wgt[4]*f[20])/m_+(1.4142135623730951*vmap[2]*wgt[3]*f[15])/m_+(1.4142135623730951*wgt[1]*vmap[3]*f[15])/m_+(1.4142135623730951*wgt[0]*vmap[3]*f[12])/m_+(1.4142135623730951*vmap[2]*wgt[2]*f[12])/m_+(1.4142135623730951*vmap[3]*wgt[3]*f[11])/m_+(1.4142135623730951*wgt[1]*vmap[2]*f[11])/m_+(1.4142135623730951*wgt[2]*vmap[3]*f[5])/m_+(1.4142135623730951*wgt[0]*vmap[2]*f[5])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (2.0*wgt[5]*f[19]+2.0*wgt[4]*f[16]+2.0*wgt[3]*f[10]+2.0*wgt[2]*f[4]+2.0*wgt[1]*f[3]+2.0*f[0]*wgt[0])*volFact; 
  out[1] += (1.4142135623730951*vmap[0]*wgt[5]*f[19]+1.264911064067352*vmap[1]*wgt[3]*f[19]+1.4142135623730951*vmap[0]*wgt[4]*f[16]+1.264911064067352*vmap[1]*wgt[1]*f[16]+1.264911064067352*vmap[1]*wgt[5]*f[10]+1.4142135623730951*vmap[0]*wgt[3]*f[10]+1.4142135623730951*vmap[1]*wgt[2]*f[10]+1.264911064067352*vmap[1]*f[3]*wgt[4]+1.4142135623730951*vmap[1]*wgt[3]*f[4]+1.4142135623730951*vmap[0]*wgt[2]*f[4]+1.4142135623730951*vmap[0]*wgt[1]*f[3]+1.4142135623730951*wgt[0]*vmap[1]*f[3]+1.4142135623730951*f[0]*vmap[1]*wgt[1]+1.4142135623730951*f[0]*vmap[0]*wgt[0])*volFact; 
  out[2] += (1.5714285714285714*vmap1R2*wgt[5]*f[19]+vmap0R2*wgt[5]*f[19]+1.7888543819998322*vmap[0]*vmap[1]*wgt[3]*f[19]+0.8944271909999161*vmap1R2*wgt[2]*f[19]+1.5714285714285714*vmap1R2*wgt[4]*f[16]+vmap0R2*wgt[4]*f[16]+1.7888543819998317*vmap[0]*vmap[1]*wgt[1]*f[16]+0.8944271909999159*wgt[0]*vmap1R2*f[16]+1.7888543819998322*vmap[0]*vmap[1]*wgt[5]*f[10]+1.8*vmap1R2*wgt[3]*f[10]+vmap0R2*wgt[3]*f[10]+2.0*vmap[0]*vmap[1]*wgt[2]*f[10]+0.8944271909999161*vmap1R2*f[4]*wgt[5]+1.7888543819998317*vmap[0]*vmap[1]*f[3]*wgt[4]+0.8944271909999159*f[0]*vmap1R2*wgt[4]+2.0*vmap[0]*vmap[1]*wgt[3]*f[4]+vmap1R2*wgt[2]*f[4]+vmap0R2*wgt[2]*f[4]+1.8*vmap1R2*wgt[1]*f[3]+vmap0R2*wgt[1]*f[3]+2.0*vmap[0]*wgt[0]*vmap[1]*f[3]+2.0*f[0]*vmap[0]*vmap[1]*wgt[1]+f[0]*wgt[0]*vmap1R2+f[0]*vmap0R2*wgt[0])*volFact; 
  out[3] += (bmag[3]*tmp[3]+bmag[2]*tmp[2]+bmag[1]*tmp[1]+bmag[0]*tmp[0])*volFact; 
} 

