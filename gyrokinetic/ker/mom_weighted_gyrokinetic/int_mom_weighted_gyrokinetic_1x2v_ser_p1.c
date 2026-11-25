#include <gkyl_mom_weighted_gyrokinetic_kernels.h> 
GKYL_CU_DH void int_mom_weighted_gyrokinetic_M0_WEIGHTX_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.7853981633974483*dxv[0]*dxv[1]*dxv[2]/m_; 
 

  out[0] += (2.0*f[1]*wgt[1]+2.0*f[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_weighted_gyrokinetic_M0_WEIGHTV_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.7853981633974483*dxv[0]*dxv[1]*dxv[2]/m_; 
 

  out[0] += (1.4142135623730951*wgt[5]*f[10]+1.4142135623730951*wgt[4]*f[8]+1.4142135623730951*wgt[3]*f[6]+1.4142135623730951*wgt[2]*f[3]+1.4142135623730951*wgt[1]*f[2]+1.4142135623730951*f[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_weighted_gyrokinetic_M1_WEIGHTX_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.7853981633974483*dxv[0]*dxv[1]*dxv[2]/m_; 
 

  out[0] += (1.4142135623730951*vmap[1]*wgt[1]*f[4]+1.4142135623730951*wgt[0]*vmap[1]*f[2]+1.4142135623730951*vmap[0]*f[1]*wgt[1]+1.4142135623730951*f[0]*vmap[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_weighted_gyrokinetic_M1_WEIGHTV_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.7853981633974483*dxv[0]*dxv[1]*dxv[2]/m_; 
 

  out[0] += (vmap[0]*wgt[5]*f[10]+0.8944271909999161*vmap[1]*wgt[3]*f[10]+vmap[0]*wgt[4]*f[8]+0.8944271909999159*vmap[1]*wgt[1]*f[8]+0.8944271909999161*vmap[1]*wgt[5]*f[6]+vmap[0]*wgt[3]*f[6]+vmap[1]*wgt[2]*f[6]+0.8944271909999159*vmap[1]*f[2]*wgt[4]+vmap[1]*f[3]*wgt[3]+vmap[0]*wgt[2]*f[3]+vmap[0]*wgt[1]*f[2]+wgt[0]*vmap[1]*f[2]+f[0]*vmap[1]*wgt[1]+f[0]*vmap[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_weighted_gyrokinetic_M2PAR_WEIGHTX_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.7853981633974483*dxv[0]*dxv[1]*dxv[2]/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (0.8944271909999161*vmap1R2*wgt[1]*f[9]+0.8944271909999159*wgt[0]*vmap1R2*f[8]+2.0*vmap[0]*vmap[1]*wgt[1]*f[4]+2.0*vmap[0]*wgt[0]*vmap[1]*f[2]+f[1]*vmap1R2*wgt[1]+vmap0R2*f[1]*wgt[1]+f[0]*wgt[0]*vmap1R2+f[0]*vmap0R2*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_weighted_gyrokinetic_M2PAR_WEIGHTV_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.7853981633974483*dxv[0]*dxv[1]*dxv[2]/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (1.1111677990074318*vmap1R2*wgt[5]*f[10]+0.7071067811865475*vmap0R2*wgt[5]*f[10]+1.264911064067352*vmap[0]*vmap[1]*wgt[3]*f[10]+0.632455532033676*vmap1R2*wgt[2]*f[10]+1.1111677990074318*vmap1R2*wgt[4]*f[8]+0.7071067811865475*vmap0R2*wgt[4]*f[8]+1.264911064067352*vmap[0]*vmap[1]*wgt[1]*f[8]+0.6324555320336759*wgt[0]*vmap1R2*f[8]+1.264911064067352*vmap[0]*vmap[1]*wgt[5]*f[6]+1.2727922061357855*vmap1R2*wgt[3]*f[6]+0.7071067811865475*vmap0R2*wgt[3]*f[6]+1.4142135623730951*vmap[0]*vmap[1]*wgt[2]*f[6]+0.632455532033676*vmap1R2*f[3]*wgt[5]+1.264911064067352*vmap[0]*vmap[1]*f[2]*wgt[4]+0.6324555320336759*f[0]*vmap1R2*wgt[4]+1.4142135623730951*vmap[0]*vmap[1]*f[3]*wgt[3]+0.7071067811865475*vmap1R2*wgt[2]*f[3]+0.7071067811865475*vmap0R2*wgt[2]*f[3]+1.2727922061357855*vmap1R2*wgt[1]*f[2]+0.7071067811865475*vmap0R2*wgt[1]*f[2]+1.4142135623730951*vmap[0]*wgt[0]*vmap[1]*f[2]+1.4142135623730951*f[0]*vmap[0]*vmap[1]*wgt[1]+0.7071067811865475*f[0]*wgt[0]*vmap1R2+0.7071067811865475*f[0]*vmap0R2*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_weighted_gyrokinetic_M2_WEIGHTX_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.7853981633974483*dxv[0]*dxv[1]*dxv[2]/m_; 
 
  double tmp[2]; 
  tmp[0] = (2.8284271247461907*f[3]*vmap[3])/m_+(2.8284271247461907*f[0]*vmap[2])/m_; 
  tmp[1] = (2.8284271247461907*vmap[3]*f[5])/m_+(2.8284271247461907*f[1]*vmap[2])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (0.8944271909999161*vmap1R2*wgt[1]*f[9]+0.8944271909999159*wgt[0]*vmap1R2*f[8]+2.0*vmap[0]*vmap[1]*wgt[1]*f[4]+2.0*vmap[0]*wgt[0]*vmap[1]*f[2]+f[1]*vmap1R2*wgt[1]+0.7071067811865475*bmag[0]*tmp[1]*wgt[1]+vmap0R2*f[1]*wgt[1]+0.7071067811865475*tmp[0]*bmag[1]*wgt[1]+f[0]*wgt[0]*vmap1R2+0.7071067811865475*wgt[0]*bmag[1]*tmp[1]+f[0]*vmap0R2*wgt[0]+0.7071067811865475*bmag[0]*tmp[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_weighted_gyrokinetic_M2_WEIGHTV_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.7853981633974483*dxv[0]*dxv[1]*dxv[2]/m_; 
 
  double tmp[2]; 
  tmp[0] = (1.4142135623730951*vmap[2]*wgt[5]*f[10])/m_+(1.4142135623730951*vmap[3]*wgt[4]*f[10])/m_+(1.4142135623730951*vmap[3]*wgt[5]*f[8])/m_+(1.4142135623730951*vmap[2]*wgt[4]*f[8])/m_+(1.4142135623730951*vmap[2]*wgt[3]*f[6])/m_+(1.4142135623730951*wgt[1]*vmap[3]*f[6])/m_+(1.4142135623730951*f[2]*vmap[3]*wgt[3])/m_+(1.4142135623730951*wgt[0]*f[3]*vmap[3])/m_+(1.4142135623730951*f[0]*wgt[2]*vmap[3])/m_+(1.4142135623730951*vmap[2]*wgt[2]*f[3])/m_+(1.4142135623730951*wgt[1]*f[2]*vmap[2])/m_+(1.4142135623730951*f[0]*wgt[0]*vmap[2])/m_; 
  tmp[1] = (1.4142135623730951*vmap[2]*wgt[5]*f[11])/m_+(1.4142135623730951*vmap[3]*wgt[4]*f[11])/m_+(1.4142135623730951*vmap[3]*wgt[5]*f[9])/m_+(1.4142135623730951*vmap[2]*wgt[4]*f[9])/m_+(1.4142135623730951*vmap[2]*wgt[3]*f[7])/m_+(1.4142135623730951*wgt[1]*vmap[3]*f[7])/m_+(1.4142135623730951*wgt[0]*vmap[3]*f[5])/m_+(1.4142135623730951*vmap[2]*wgt[2]*f[5])/m_+(1.4142135623730951*vmap[3]*wgt[3]*f[4])/m_+(1.4142135623730951*wgt[1]*vmap[2]*f[4])/m_+(1.4142135623730951*f[1]*wgt[2]*vmap[3])/m_+(1.4142135623730951*wgt[0]*f[1]*vmap[2])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (1.1111677990074318*vmap1R2*wgt[5]*f[10]+0.7071067811865475*vmap0R2*wgt[5]*f[10]+1.264911064067352*vmap[0]*vmap[1]*wgt[3]*f[10]+0.632455532033676*vmap1R2*wgt[2]*f[10]+1.1111677990074318*vmap1R2*wgt[4]*f[8]+0.7071067811865475*vmap0R2*wgt[4]*f[8]+1.264911064067352*vmap[0]*vmap[1]*wgt[1]*f[8]+0.6324555320336759*wgt[0]*vmap1R2*f[8]+1.264911064067352*vmap[0]*vmap[1]*wgt[5]*f[6]+1.2727922061357855*vmap1R2*wgt[3]*f[6]+0.7071067811865475*vmap0R2*wgt[3]*f[6]+1.4142135623730951*vmap[0]*vmap[1]*wgt[2]*f[6]+0.632455532033676*vmap1R2*f[3]*wgt[5]+1.264911064067352*vmap[0]*vmap[1]*f[2]*wgt[4]+0.6324555320336759*f[0]*vmap1R2*wgt[4]+1.4142135623730951*vmap[0]*vmap[1]*f[3]*wgt[3]+0.7071067811865475*vmap1R2*wgt[2]*f[3]+0.7071067811865475*vmap0R2*wgt[2]*f[3]+1.2727922061357855*vmap1R2*wgt[1]*f[2]+0.7071067811865475*vmap0R2*wgt[1]*f[2]+1.4142135623730951*vmap[0]*wgt[0]*vmap[1]*f[2]+1.4142135623730951*f[0]*vmap[0]*vmap[1]*wgt[1]+0.7071067811865475*f[0]*wgt[0]*vmap1R2+bmag[1]*tmp[1]+0.7071067811865475*f[0]*vmap0R2*wgt[0]+bmag[0]*tmp[0])*volFact; 
} 

GKYL_CU_DH void int_mom_weighted_gyrokinetic_M0M1M2_WEIGHTX_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.7853981633974483*dxv[0]*dxv[1]*dxv[2]/m_; 
 
  double tmp[2]; 
  tmp[0] = (2.8284271247461907*f[3]*vmap[3])/m_+(2.8284271247461907*f[0]*vmap[2])/m_; 
  tmp[1] = (2.8284271247461907*vmap[3]*f[5])/m_+(2.8284271247461907*f[1]*vmap[2])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (2.0*f[1]*wgt[1]+2.0*f[0]*wgt[0])*volFact; 
  out[1] += (1.4142135623730951*vmap[1]*wgt[1]*f[4]+1.4142135623730951*wgt[0]*vmap[1]*f[2]+1.4142135623730951*vmap[0]*f[1]*wgt[1]+1.4142135623730951*f[0]*vmap[0]*wgt[0])*volFact; 
  out[2] += (0.8944271909999161*vmap1R2*wgt[1]*f[9]+0.8944271909999159*wgt[0]*vmap1R2*f[8]+2.0*vmap[0]*vmap[1]*wgt[1]*f[4]+2.0*vmap[0]*wgt[0]*vmap[1]*f[2]+f[1]*vmap1R2*wgt[1]+0.7071067811865475*bmag[0]*tmp[1]*wgt[1]+vmap0R2*f[1]*wgt[1]+0.7071067811865475*tmp[0]*bmag[1]*wgt[1]+f[0]*wgt[0]*vmap1R2+0.7071067811865475*wgt[0]*bmag[1]*tmp[1]+f[0]*vmap0R2*wgt[0]+0.7071067811865475*bmag[0]*tmp[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_weighted_gyrokinetic_M0M1M2_WEIGHTV_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.7853981633974483*dxv[0]*dxv[1]*dxv[2]/m_; 
 
  double tmp[2]; 
  tmp[0] = (1.4142135623730951*vmap[2]*wgt[5]*f[10])/m_+(1.4142135623730951*vmap[3]*wgt[4]*f[10])/m_+(1.4142135623730951*vmap[3]*wgt[5]*f[8])/m_+(1.4142135623730951*vmap[2]*wgt[4]*f[8])/m_+(1.4142135623730951*vmap[2]*wgt[3]*f[6])/m_+(1.4142135623730951*wgt[1]*vmap[3]*f[6])/m_+(1.4142135623730951*f[2]*vmap[3]*wgt[3])/m_+(1.4142135623730951*wgt[0]*f[3]*vmap[3])/m_+(1.4142135623730951*f[0]*wgt[2]*vmap[3])/m_+(1.4142135623730951*vmap[2]*wgt[2]*f[3])/m_+(1.4142135623730951*wgt[1]*f[2]*vmap[2])/m_+(1.4142135623730951*f[0]*wgt[0]*vmap[2])/m_; 
  tmp[1] = (1.4142135623730951*vmap[2]*wgt[5]*f[11])/m_+(1.4142135623730951*vmap[3]*wgt[4]*f[11])/m_+(1.4142135623730951*vmap[3]*wgt[5]*f[9])/m_+(1.4142135623730951*vmap[2]*wgt[4]*f[9])/m_+(1.4142135623730951*vmap[2]*wgt[3]*f[7])/m_+(1.4142135623730951*wgt[1]*vmap[3]*f[7])/m_+(1.4142135623730951*wgt[0]*vmap[3]*f[5])/m_+(1.4142135623730951*vmap[2]*wgt[2]*f[5])/m_+(1.4142135623730951*vmap[3]*wgt[3]*f[4])/m_+(1.4142135623730951*wgt[1]*vmap[2]*f[4])/m_+(1.4142135623730951*f[1]*wgt[2]*vmap[3])/m_+(1.4142135623730951*wgt[0]*f[1]*vmap[2])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (1.4142135623730951*wgt[5]*f[10]+1.4142135623730951*wgt[4]*f[8]+1.4142135623730951*wgt[3]*f[6]+1.4142135623730951*wgt[2]*f[3]+1.4142135623730951*wgt[1]*f[2]+1.4142135623730951*f[0]*wgt[0])*volFact; 
  out[1] += (vmap[0]*wgt[5]*f[10]+0.8944271909999161*vmap[1]*wgt[3]*f[10]+vmap[0]*wgt[4]*f[8]+0.8944271909999159*vmap[1]*wgt[1]*f[8]+0.8944271909999161*vmap[1]*wgt[5]*f[6]+vmap[0]*wgt[3]*f[6]+vmap[1]*wgt[2]*f[6]+0.8944271909999159*vmap[1]*f[2]*wgt[4]+vmap[1]*f[3]*wgt[3]+vmap[0]*wgt[2]*f[3]+vmap[0]*wgt[1]*f[2]+wgt[0]*vmap[1]*f[2]+f[0]*vmap[1]*wgt[1]+f[0]*vmap[0]*wgt[0])*volFact; 
  out[2] += (1.1111677990074318*vmap1R2*wgt[5]*f[10]+0.7071067811865475*vmap0R2*wgt[5]*f[10]+1.264911064067352*vmap[0]*vmap[1]*wgt[3]*f[10]+0.632455532033676*vmap1R2*wgt[2]*f[10]+1.1111677990074318*vmap1R2*wgt[4]*f[8]+0.7071067811865475*vmap0R2*wgt[4]*f[8]+1.264911064067352*vmap[0]*vmap[1]*wgt[1]*f[8]+0.6324555320336759*wgt[0]*vmap1R2*f[8]+1.264911064067352*vmap[0]*vmap[1]*wgt[5]*f[6]+1.2727922061357855*vmap1R2*wgt[3]*f[6]+0.7071067811865475*vmap0R2*wgt[3]*f[6]+1.4142135623730951*vmap[0]*vmap[1]*wgt[2]*f[6]+0.632455532033676*vmap1R2*f[3]*wgt[5]+1.264911064067352*vmap[0]*vmap[1]*f[2]*wgt[4]+0.6324555320336759*f[0]*vmap1R2*wgt[4]+1.4142135623730951*vmap[0]*vmap[1]*f[3]*wgt[3]+0.7071067811865475*vmap1R2*wgt[2]*f[3]+0.7071067811865475*vmap0R2*wgt[2]*f[3]+1.2727922061357855*vmap1R2*wgt[1]*f[2]+0.7071067811865475*vmap0R2*wgt[1]*f[2]+1.4142135623730951*vmap[0]*wgt[0]*vmap[1]*f[2]+1.4142135623730951*f[0]*vmap[0]*vmap[1]*wgt[1]+0.7071067811865475*f[0]*wgt[0]*vmap1R2+bmag[1]*tmp[1]+0.7071067811865475*f[0]*vmap0R2*wgt[0]+bmag[0]*tmp[0])*volFact; 
} 

GKYL_CU_DH void int_mom_weighted_gyrokinetic_M2PERP_WEIGHTX_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.7853981633974483*dxv[0]*dxv[1]*dxv[2]/m_; 
 
  double tmp[2]; 
  tmp[0] = (2.8284271247461907*f[3]*vmap[3])/m_+(2.8284271247461907*f[0]*vmap[2])/m_; 
  tmp[1] = (2.8284271247461907*vmap[3]*f[5])/m_+(2.8284271247461907*f[1]*vmap[2])/m_; 
 

  out[0] += (0.7071067811865475*bmag[0]*tmp[1]*wgt[1]+0.7071067811865475*tmp[0]*bmag[1]*wgt[1]+0.7071067811865475*wgt[0]*bmag[1]*tmp[1]+0.7071067811865475*bmag[0]*tmp[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_weighted_gyrokinetic_M2PERP_WEIGHTV_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.7853981633974483*dxv[0]*dxv[1]*dxv[2]/m_; 
 
  double tmp[2]; 
  tmp[0] = (1.4142135623730951*vmap[2]*wgt[5]*f[10])/m_+(1.4142135623730951*vmap[3]*wgt[4]*f[10])/m_+(1.4142135623730951*vmap[3]*wgt[5]*f[8])/m_+(1.4142135623730951*vmap[2]*wgt[4]*f[8])/m_+(1.4142135623730951*vmap[2]*wgt[3]*f[6])/m_+(1.4142135623730951*wgt[1]*vmap[3]*f[6])/m_+(1.4142135623730951*f[2]*vmap[3]*wgt[3])/m_+(1.4142135623730951*wgt[0]*f[3]*vmap[3])/m_+(1.4142135623730951*f[0]*wgt[2]*vmap[3])/m_+(1.4142135623730951*vmap[2]*wgt[2]*f[3])/m_+(1.4142135623730951*wgt[1]*f[2]*vmap[2])/m_+(1.4142135623730951*f[0]*wgt[0]*vmap[2])/m_; 
  tmp[1] = (1.4142135623730951*vmap[2]*wgt[5]*f[11])/m_+(1.4142135623730951*vmap[3]*wgt[4]*f[11])/m_+(1.4142135623730951*vmap[3]*wgt[5]*f[9])/m_+(1.4142135623730951*vmap[2]*wgt[4]*f[9])/m_+(1.4142135623730951*vmap[2]*wgt[3]*f[7])/m_+(1.4142135623730951*wgt[1]*vmap[3]*f[7])/m_+(1.4142135623730951*wgt[0]*vmap[3]*f[5])/m_+(1.4142135623730951*vmap[2]*wgt[2]*f[5])/m_+(1.4142135623730951*vmap[3]*wgt[3]*f[4])/m_+(1.4142135623730951*wgt[1]*vmap[2]*f[4])/m_+(1.4142135623730951*f[1]*wgt[2]*vmap[3])/m_+(1.4142135623730951*wgt[0]*f[1]*vmap[2])/m_; 
 

  out[0] += (bmag[1]*tmp[1]+bmag[0]*tmp[0])*volFact; 
} 

GKYL_CU_DH void int_mom_weighted_gyrokinetic_M0M1M2PARM2PERP_WEIGHTX_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.7853981633974483*dxv[0]*dxv[1]*dxv[2]/m_; 
 
  double tmp[2]; 
  tmp[0] = (2.8284271247461907*f[3]*vmap[3])/m_+(2.8284271247461907*f[0]*vmap[2])/m_; 
  tmp[1] = (2.8284271247461907*vmap[3]*f[5])/m_+(2.8284271247461907*f[1]*vmap[2])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (2.0*f[1]*wgt[1]+2.0*f[0]*wgt[0])*volFact; 
  out[1] += (1.4142135623730951*vmap[1]*wgt[1]*f[4]+1.4142135623730951*wgt[0]*vmap[1]*f[2]+1.4142135623730951*vmap[0]*f[1]*wgt[1]+1.4142135623730951*f[0]*vmap[0]*wgt[0])*volFact; 
  out[2] += (0.8944271909999161*vmap1R2*wgt[1]*f[9]+0.8944271909999159*wgt[0]*vmap1R2*f[8]+2.0*vmap[0]*vmap[1]*wgt[1]*f[4]+2.0*vmap[0]*wgt[0]*vmap[1]*f[2]+f[1]*vmap1R2*wgt[1]+vmap0R2*f[1]*wgt[1]+f[0]*wgt[0]*vmap1R2+f[0]*vmap0R2*wgt[0])*volFact; 
  out[3] += (0.7071067811865475*bmag[0]*tmp[1]*wgt[1]+0.7071067811865475*tmp[0]*bmag[1]*wgt[1]+0.7071067811865475*wgt[0]*bmag[1]*tmp[1]+0.7071067811865475*bmag[0]*tmp[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_weighted_gyrokinetic_M0M1M2PARM2PERP_WEIGHTV_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.7853981633974483*dxv[0]*dxv[1]*dxv[2]/m_; 
 
  double tmp[2]; 
  tmp[0] = (1.4142135623730951*vmap[2]*wgt[5]*f[10])/m_+(1.4142135623730951*vmap[3]*wgt[4]*f[10])/m_+(1.4142135623730951*vmap[3]*wgt[5]*f[8])/m_+(1.4142135623730951*vmap[2]*wgt[4]*f[8])/m_+(1.4142135623730951*vmap[2]*wgt[3]*f[6])/m_+(1.4142135623730951*wgt[1]*vmap[3]*f[6])/m_+(1.4142135623730951*f[2]*vmap[3]*wgt[3])/m_+(1.4142135623730951*wgt[0]*f[3]*vmap[3])/m_+(1.4142135623730951*f[0]*wgt[2]*vmap[3])/m_+(1.4142135623730951*vmap[2]*wgt[2]*f[3])/m_+(1.4142135623730951*wgt[1]*f[2]*vmap[2])/m_+(1.4142135623730951*f[0]*wgt[0]*vmap[2])/m_; 
  tmp[1] = (1.4142135623730951*vmap[2]*wgt[5]*f[11])/m_+(1.4142135623730951*vmap[3]*wgt[4]*f[11])/m_+(1.4142135623730951*vmap[3]*wgt[5]*f[9])/m_+(1.4142135623730951*vmap[2]*wgt[4]*f[9])/m_+(1.4142135623730951*vmap[2]*wgt[3]*f[7])/m_+(1.4142135623730951*wgt[1]*vmap[3]*f[7])/m_+(1.4142135623730951*wgt[0]*vmap[3]*f[5])/m_+(1.4142135623730951*vmap[2]*wgt[2]*f[5])/m_+(1.4142135623730951*vmap[3]*wgt[3]*f[4])/m_+(1.4142135623730951*wgt[1]*vmap[2]*f[4])/m_+(1.4142135623730951*f[1]*wgt[2]*vmap[3])/m_+(1.4142135623730951*wgt[0]*f[1]*vmap[2])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (1.4142135623730951*wgt[5]*f[10]+1.4142135623730951*wgt[4]*f[8]+1.4142135623730951*wgt[3]*f[6]+1.4142135623730951*wgt[2]*f[3]+1.4142135623730951*wgt[1]*f[2]+1.4142135623730951*f[0]*wgt[0])*volFact; 
  out[1] += (vmap[0]*wgt[5]*f[10]+0.8944271909999161*vmap[1]*wgt[3]*f[10]+vmap[0]*wgt[4]*f[8]+0.8944271909999159*vmap[1]*wgt[1]*f[8]+0.8944271909999161*vmap[1]*wgt[5]*f[6]+vmap[0]*wgt[3]*f[6]+vmap[1]*wgt[2]*f[6]+0.8944271909999159*vmap[1]*f[2]*wgt[4]+vmap[1]*f[3]*wgt[3]+vmap[0]*wgt[2]*f[3]+vmap[0]*wgt[1]*f[2]+wgt[0]*vmap[1]*f[2]+f[0]*vmap[1]*wgt[1]+f[0]*vmap[0]*wgt[0])*volFact; 
  out[2] += (1.1111677990074318*vmap1R2*wgt[5]*f[10]+0.7071067811865475*vmap0R2*wgt[5]*f[10]+1.264911064067352*vmap[0]*vmap[1]*wgt[3]*f[10]+0.632455532033676*vmap1R2*wgt[2]*f[10]+1.1111677990074318*vmap1R2*wgt[4]*f[8]+0.7071067811865475*vmap0R2*wgt[4]*f[8]+1.264911064067352*vmap[0]*vmap[1]*wgt[1]*f[8]+0.6324555320336759*wgt[0]*vmap1R2*f[8]+1.264911064067352*vmap[0]*vmap[1]*wgt[5]*f[6]+1.2727922061357855*vmap1R2*wgt[3]*f[6]+0.7071067811865475*vmap0R2*wgt[3]*f[6]+1.4142135623730951*vmap[0]*vmap[1]*wgt[2]*f[6]+0.632455532033676*vmap1R2*f[3]*wgt[5]+1.264911064067352*vmap[0]*vmap[1]*f[2]*wgt[4]+0.6324555320336759*f[0]*vmap1R2*wgt[4]+1.4142135623730951*vmap[0]*vmap[1]*f[3]*wgt[3]+0.7071067811865475*vmap1R2*wgt[2]*f[3]+0.7071067811865475*vmap0R2*wgt[2]*f[3]+1.2727922061357855*vmap1R2*wgt[1]*f[2]+0.7071067811865475*vmap0R2*wgt[1]*f[2]+1.4142135623730951*vmap[0]*wgt[0]*vmap[1]*f[2]+1.4142135623730951*f[0]*vmap[0]*vmap[1]*wgt[1]+0.7071067811865475*f[0]*wgt[0]*vmap1R2+0.7071067811865475*f[0]*vmap0R2*wgt[0])*volFact; 
  out[3] += (bmag[1]*tmp[1]+bmag[0]*tmp[0])*volFact; 
} 

