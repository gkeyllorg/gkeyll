#include <gkyl_mom_gyrokinetic_kernels.h> 
GKYL_CU_DH void mom_gyrokinetic_M0_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[2]*dxv[3]/m_; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[5]*volFact; 
  out[4] += 2.0*f[11]*volFact; 
  out[5] += 2.0*f[12]*volFact; 
  out[6] += 2.0*f[19]*volFact; 
  out[7] += 2.0*f[20]*volFact; 
} 

GKYL_CU_DH void mom_gyrokinetic_M1_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[2]*dxv[3]/m_; 
  out[0] += (1.4142135623730951*vmap[1]*f[3]+1.4142135623730951*f[0]*vmap[0])*volFact; 
  out[1] += (1.4142135623730951*vmap[1]*f[6]+1.4142135623730951*vmap[0]*f[1])*volFact; 
  out[2] += (1.4142135623730951*vmap[1]*f[7]+1.4142135623730951*vmap[0]*f[2])*volFact; 
  out[3] += (1.4142135623730951*vmap[1]*f[15]+1.4142135623730951*vmap[0]*f[5])*volFact; 
  out[4] += (1.4142135623730951*vmap[1]*f[21]+1.4142135623730951*vmap[0]*f[11])*volFact; 
  out[5] += (1.4142135623730951*vmap[1]*f[22]+1.4142135623730951*vmap[0]*f[12])*volFact; 
  out[6] += (1.4142135623730951*vmap[1]*f[32]+1.4142135623730951*vmap[0]*f[19])*volFact; 
  out[7] += (1.4142135623730951*vmap[1]*f[33]+1.4142135623730951*vmap[0]*f[20])*volFact; 
} 

GKYL_CU_DH void mom_gyrokinetic_M2_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[2]*dxv[3]/m_; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (0.8944271909999159*vmap1R2*f[13]+2.0*vmap[0]*vmap[1]*f[3]+f[0]*vmap1R2+f[0]*vmap0R2)*volFact; 
  out[1] += (0.8944271909999161*vmap1R2*f[23]+2.0*vmap[0]*vmap[1]*f[6]+f[1]*vmap1R2+vmap0R2*f[1])*volFact; 
  out[2] += (0.8944271909999161*vmap1R2*f[24]+2.0*vmap[0]*vmap[1]*f[7]+vmap1R2*f[2]+vmap0R2*f[2])*volFact; 
  out[3] += (0.8944271909999159*vmap1R2*f[34]+2.0*vmap[0]*vmap[1]*f[15]+vmap1R2*f[5]+vmap0R2*f[5])*volFact; 
  out[4] += (2.0000000000000004*vmap[0]*vmap[1]*f[21]+vmap1R2*f[11]+vmap0R2*f[11])*volFact; 
  out[5] += (2.0000000000000004*vmap[0]*vmap[1]*f[22]+vmap1R2*f[12]+vmap0R2*f[12])*volFact; 
  out[6] += (2.0000000000000004*vmap[0]*vmap[1]*f[32]+vmap1R2*f[19]+vmap0R2*f[19])*volFact; 
  out[7] += (2.0000000000000004*vmap[0]*vmap[1]*f[33]+vmap1R2*f[20]+vmap0R2*f[20])*volFact; 
  double tmp[8]; 
  tmp[0] = 1.4142135623730951*vmap[3]*f[4]+1.4142135623730951*f[0]*vmap[2]; 
  tmp[1] = 1.4142135623730951*vmap[3]*f[8]+1.4142135623730951*f[1]*vmap[2]; 
  tmp[2] = 1.4142135623730951*vmap[3]*f[9]+1.4142135623730951*f[2]*vmap[2]; 
  tmp[3] = 1.4142135623730951*vmap[3]*f[16]+1.4142135623730951*vmap[2]*f[5]; 
  tmp[4] = 1.4142135623730951*vmap[3]*f[25]+1.4142135623730951*vmap[2]*f[11]; 
  tmp[5] = 1.4142135623730951*vmap[3]*f[26]+1.4142135623730951*vmap[2]*f[12]; 
  tmp[6] = 1.4142135623730951*vmap[3]*f[35]+1.4142135623730951*vmap[2]*f[19]; 
  tmp[7] = 1.4142135623730951*vmap[3]*f[36]+1.4142135623730951*vmap[2]*f[20]; 
  out[0] += (2.0*(0.5*bmag[7]*tmp[7]+0.5*bmag[6]*tmp[6]+0.5*bmag[5]*tmp[5]+0.5*bmag[4]*tmp[4]+0.5*bmag[3]*tmp[3]+0.5*bmag[2]*tmp[2]+0.5*bmag[1]*tmp[1]+0.5*bmag[0]*tmp[0])*volFact)/m_; 
  out[1] += (2.0*(0.5000000000000001*bmag[5]*tmp[7]+0.5000000000000001*tmp[5]*bmag[7]+0.44721359549995804*bmag[3]*tmp[6]+0.44721359549995804*tmp[3]*bmag[6]+0.4472135954999579*bmag[1]*tmp[4]+0.4472135954999579*tmp[1]*bmag[4]+0.5*bmag[2]*tmp[3]+0.5*tmp[2]*bmag[3]+0.5*bmag[0]*tmp[1]+0.5*tmp[0]*bmag[1])*volFact)/m_; 
  out[2] += (2.0*(0.44721359549995804*bmag[3]*tmp[7]+0.44721359549995804*tmp[3]*bmag[7]+0.5000000000000001*bmag[4]*tmp[6]+0.5000000000000001*tmp[4]*bmag[6]+0.4472135954999579*bmag[2]*tmp[5]+0.4472135954999579*tmp[2]*bmag[5]+0.5*bmag[1]*tmp[3]+0.5*tmp[1]*bmag[3]+0.5*bmag[0]*tmp[2]+0.5*tmp[0]*bmag[2])*volFact)/m_; 
  out[3] += (2.0*(0.4*bmag[6]*tmp[7]+0.44721359549995804*bmag[2]*tmp[7]+0.4*tmp[6]*bmag[7]+0.44721359549995804*tmp[2]*bmag[7]+0.44721359549995804*bmag[1]*tmp[6]+0.44721359549995804*tmp[1]*bmag[6]+0.4472135954999579*bmag[3]*tmp[5]+0.4472135954999579*tmp[3]*bmag[5]+0.4472135954999579*bmag[3]*tmp[4]+0.4472135954999579*tmp[3]*bmag[4]+0.5*bmag[0]*tmp[3]+0.5*tmp[0]*bmag[3]+0.5*bmag[1]*tmp[2]+0.5*tmp[1]*bmag[2])*volFact)/m_; 
  out[4] += (2.0*(0.4472135954999579*bmag[7]*tmp[7]+0.31943828249996997*bmag[6]*tmp[6]+0.5000000000000001*bmag[2]*tmp[6]+0.5000000000000001*tmp[2]*bmag[6]+0.31943828249996997*bmag[4]*tmp[4]+0.5*bmag[0]*tmp[4]+0.5*tmp[0]*bmag[4]+0.4472135954999579*bmag[3]*tmp[3]+0.4472135954999579*bmag[1]*tmp[1])*volFact)/m_; 
  out[5] += (2.0*(0.31943828249996997*bmag[7]*tmp[7]+0.5000000000000001*bmag[1]*tmp[7]+0.5000000000000001*tmp[1]*bmag[7]+0.4472135954999579*bmag[6]*tmp[6]+0.31943828249996997*bmag[5]*tmp[5]+0.5*bmag[0]*tmp[5]+0.5*tmp[0]*bmag[5]+0.4472135954999579*bmag[3]*tmp[3]+0.4472135954999579*bmag[2]*tmp[2])*volFact)/m_; 
  out[6] += (2.0*(0.4*bmag[3]*tmp[7]+0.4*tmp[3]*bmag[7]+0.4472135954999579*bmag[5]*tmp[6]+0.31943828249996997*bmag[4]*tmp[6]+0.5*bmag[0]*tmp[6]+0.4472135954999579*tmp[5]*bmag[6]+0.31943828249996997*tmp[4]*bmag[6]+0.5*tmp[0]*bmag[6]+0.5000000000000001*bmag[2]*tmp[4]+0.5000000000000001*tmp[2]*bmag[4]+0.44721359549995804*bmag[1]*tmp[3]+0.44721359549995804*tmp[1]*bmag[3])*volFact)/m_; 
  out[7] += (2.0*(0.31943828249996997*bmag[5]*tmp[7]+0.4472135954999579*bmag[4]*tmp[7]+0.5*bmag[0]*tmp[7]+0.31943828249996997*tmp[5]*bmag[7]+0.4472135954999579*tmp[4]*bmag[7]+0.5*tmp[0]*bmag[7]+0.4*bmag[3]*tmp[6]+0.4*tmp[3]*bmag[6]+0.5000000000000001*bmag[1]*tmp[5]+0.5000000000000001*tmp[1]*bmag[5]+0.44721359549995804*bmag[2]*tmp[3]+0.44721359549995804*tmp[2]*bmag[3])*volFact)/m_; 
} 

GKYL_CU_DH void mom_gyrokinetic_M3_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[2]*dxv[3]/m_; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap0R3 = pow(vmap[0],3);
  const double vmap1R2 = pow(vmap[1],2);
  const double vmap1R3 = pow(vmap[1],3);

  out[0] += ((vmap[1]*vmap[3]*bmag[7]*f[45]+vmap[1]*vmap[3]*bmag[6]*f[44]+vmap[1]*vmap[3]*bmag[5]*f[38]+vmap[1]*vmap[3]*bmag[4]*f[37]+1.0000000000000002*vmap[0]*vmap[3]*bmag[7]*f[36]+1.0000000000000002*vmap[0]*vmap[3]*bmag[6]*f[35]+1.0000000000000002*vmap[1]*vmap[2]*bmag[7]*f[33]+1.0000000000000002*vmap[1]*vmap[2]*bmag[6]*f[32]+vmap[1]*bmag[3]*vmap[3]*f[31]+1.0000000000000002*vmap[0]*vmap[3]*bmag[5]*f[26]+1.0000000000000002*vmap[0]*vmap[3]*bmag[4]*f[25]+1.0000000000000002*vmap[1]*vmap[2]*bmag[5]*f[22]+1.0000000000000002*vmap[1]*vmap[2]*bmag[4]*f[21]+vmap[0]*vmap[2]*bmag[7]*f[20]+vmap[0]*vmap[2]*bmag[6]*f[19]+vmap[1]*bmag[2]*vmap[3]*f[18]+bmag[1]*vmap[1]*vmap[3]*f[17]+vmap[0]*bmag[3]*vmap[3]*f[16]+vmap[1]*vmap[2]*bmag[3]*f[15]+vmap[0]*vmap[2]*bmag[5]*f[12]+vmap[0]*vmap[2]*bmag[4]*f[11]+bmag[0]*vmap[1]*vmap[3]*f[10]+vmap[0]*bmag[2]*vmap[3]*f[9]+vmap[0]*bmag[1]*vmap[3]*f[8]+vmap[1]*bmag[2]*vmap[2]*f[7]+bmag[1]*vmap[1]*vmap[2]*f[6]+vmap[0]*vmap[2]*bmag[3]*f[5]+bmag[0]*vmap[0]*vmap[3]*f[4]+bmag[0]*vmap[1]*vmap[2]*f[3]+vmap[0]*bmag[2]*f[2]*vmap[2]+vmap[0]*bmag[1]*f[1]*vmap[2]+bmag[0]*f[0]*vmap[0]*vmap[2])*volFact)/m_+(1.8973665961010278*vmap[0]*vmap1R2*f[13]+1.2727922061357855*vmap1R3*f[3]+2.1213203435596424*vmap0R2*vmap[1]*f[3]+2.1213203435596424*f[0]*vmap[0]*vmap1R2+0.7071067811865475*f[0]*vmap0R3)*volFact; 
  out[1] += ((1.0000000000000002*vmap[1]*vmap[3]*bmag[5]*f[45]+0.8944271909999161*vmap[1]*bmag[3]*vmap[3]*f[44]+1.0000000000000002*vmap[1]*vmap[3]*bmag[7]*f[38]+0.8944271909999159*bmag[1]*vmap[1]*vmap[3]*f[37]+vmap[0]*vmap[3]*bmag[5]*f[36]+0.8944271909999159*vmap[0]*bmag[3]*vmap[3]*f[35]+vmap[1]*vmap[2]*bmag[5]*f[33]+0.8944271909999159*vmap[1]*vmap[2]*bmag[3]*f[32]+0.8944271909999161*vmap[1]*vmap[3]*bmag[6]*f[31]+vmap[1]*bmag[2]*vmap[3]*f[31]+vmap[0]*vmap[3]*bmag[7]*f[26]+0.8944271909999161*vmap[0]*bmag[1]*vmap[3]*f[25]+vmap[1]*vmap[2]*bmag[7]*f[22]+0.8944271909999161*bmag[1]*vmap[1]*vmap[2]*f[21]+1.0000000000000002*vmap[0]*vmap[2]*bmag[5]*f[20]+0.8944271909999161*vmap[0]*vmap[2]*bmag[3]*f[19]+vmap[1]*bmag[3]*vmap[3]*f[18]+0.8944271909999159*vmap[1]*vmap[3]*bmag[4]*f[17]+bmag[0]*vmap[1]*vmap[3]*f[17]+0.8944271909999161*vmap[0]*vmap[3]*bmag[6]*f[16]+vmap[0]*bmag[2]*vmap[3]*f[16]+0.8944271909999161*vmap[1]*vmap[2]*bmag[6]*f[15]+vmap[1]*bmag[2]*vmap[2]*f[15]+1.0000000000000002*vmap[0]*vmap[2]*bmag[7]*f[12]+0.8944271909999159*vmap[0]*bmag[1]*vmap[2]*f[11]+bmag[1]*vmap[1]*vmap[3]*f[10]+vmap[0]*bmag[3]*vmap[3]*f[9]+0.8944271909999159*vmap[0]*vmap[3]*bmag[4]*f[8]+bmag[0]*vmap[0]*vmap[3]*f[8]+vmap[1]*vmap[2]*bmag[3]*f[7]+0.8944271909999159*vmap[1]*vmap[2]*bmag[4]*f[6]+bmag[0]*vmap[1]*vmap[2]*f[6]+0.8944271909999161*vmap[0]*vmap[2]*f[5]*bmag[6]+vmap[0]*bmag[2]*vmap[2]*f[5]+vmap[0]*bmag[1]*vmap[3]*f[4]+0.8944271909999159*vmap[0]*f[1]*vmap[2]*bmag[4]+bmag[1]*vmap[1]*vmap[2]*f[3]+vmap[0]*f[2]*vmap[2]*bmag[3]+bmag[0]*vmap[0]*f[1]*vmap[2]+f[0]*vmap[0]*bmag[1]*vmap[2])*volFact)/m_+(1.897366596101028*vmap[0]*vmap1R2*f[23]+1.2727922061357855*vmap1R3*f[6]+2.1213203435596424*vmap0R2*vmap[1]*f[6]+2.1213203435596424*vmap[0]*f[1]*vmap1R2+0.7071067811865475*vmap0R3*f[1])*volFact; 
  out[2] += ((0.8944271909999161*vmap[1]*bmag[3]*vmap[3]*f[45]+1.0000000000000002*vmap[1]*vmap[3]*bmag[4]*f[44]+0.8944271909999159*vmap[1]*bmag[2]*vmap[3]*f[38]+1.0000000000000002*vmap[1]*vmap[3]*bmag[6]*f[37]+0.8944271909999159*vmap[0]*bmag[3]*vmap[3]*f[36]+vmap[0]*vmap[3]*bmag[4]*f[35]+0.8944271909999159*vmap[1]*vmap[2]*bmag[3]*f[33]+vmap[1]*vmap[2]*bmag[4]*f[32]+0.8944271909999161*vmap[1]*vmap[3]*bmag[7]*f[31]+bmag[1]*vmap[1]*vmap[3]*f[31]+0.8944271909999161*vmap[0]*bmag[2]*vmap[3]*f[26]+vmap[0]*vmap[3]*bmag[6]*f[25]+0.8944271909999161*vmap[1]*bmag[2]*vmap[2]*f[22]+vmap[1]*vmap[2]*bmag[6]*f[21]+0.8944271909999161*vmap[0]*vmap[2]*bmag[3]*f[20]+1.0000000000000002*vmap[0]*vmap[2]*bmag[4]*f[19]+0.8944271909999159*vmap[1]*vmap[3]*bmag[5]*f[18]+bmag[0]*vmap[1]*vmap[3]*f[18]+vmap[1]*bmag[3]*vmap[3]*f[17]+0.8944271909999161*vmap[0]*vmap[3]*bmag[7]*f[16]+vmap[0]*bmag[1]*vmap[3]*f[16]+0.8944271909999161*vmap[1]*vmap[2]*bmag[7]*f[15]+bmag[1]*vmap[1]*vmap[2]*f[15]+0.8944271909999159*vmap[0]*bmag[2]*vmap[2]*f[12]+1.0000000000000002*vmap[0]*vmap[2]*bmag[6]*f[11]+vmap[1]*bmag[2]*vmap[3]*f[10]+0.8944271909999159*vmap[0]*vmap[3]*bmag[5]*f[9]+bmag[0]*vmap[0]*vmap[3]*f[9]+vmap[0]*bmag[3]*vmap[3]*f[8]+0.8944271909999159*vmap[1]*vmap[2]*bmag[5]*f[7]+bmag[0]*vmap[1]*vmap[2]*f[7]+0.8944271909999161*vmap[0]*vmap[2]*f[5]*bmag[7]+vmap[1]*vmap[2]*bmag[3]*f[6]+vmap[0]*bmag[1]*vmap[2]*f[5]+0.8944271909999159*vmap[0]*f[2]*vmap[2]*bmag[5]+vmap[0]*bmag[2]*vmap[3]*f[4]+vmap[1]*bmag[2]*vmap[2]*f[3]+vmap[0]*f[1]*vmap[2]*bmag[3]+bmag[0]*vmap[0]*f[2]*vmap[2]+f[0]*vmap[0]*bmag[2]*vmap[2])*volFact)/m_+(1.897366596101028*vmap[0]*vmap1R2*f[24]+1.2727922061357855*vmap1R3*f[7]+2.1213203435596424*vmap0R2*vmap[1]*f[7]+2.1213203435596424*vmap[0]*vmap1R2*f[2]+0.7071067811865475*vmap0R3*f[2])*volFact; 
  out[3] += ((0.8*vmap[1]*vmap[3]*bmag[6]*f[45]+0.8944271909999161*vmap[1]*bmag[2]*vmap[3]*f[45]+0.8*vmap[1]*vmap[3]*bmag[7]*f[44]+0.8944271909999161*bmag[1]*vmap[1]*vmap[3]*f[44]+0.8944271909999159*vmap[1]*bmag[3]*vmap[3]*f[38]+0.8944271909999159*vmap[1]*bmag[3]*vmap[3]*f[37]+0.8*vmap[0]*vmap[3]*bmag[6]*f[36]+0.8944271909999159*vmap[0]*bmag[2]*vmap[3]*f[36]+0.8*vmap[0]*vmap[3]*bmag[7]*f[35]+0.8944271909999159*vmap[0]*bmag[1]*vmap[3]*f[35]+0.8*vmap[1]*vmap[2]*bmag[6]*f[33]+0.8944271909999159*vmap[1]*bmag[2]*vmap[2]*f[33]+0.8*vmap[1]*vmap[2]*bmag[7]*f[32]+0.8944271909999159*bmag[1]*vmap[1]*vmap[2]*f[32]+0.8944271909999159*vmap[1]*vmap[3]*bmag[5]*f[31]+0.8944271909999159*vmap[1]*vmap[3]*bmag[4]*f[31]+bmag[0]*vmap[1]*vmap[3]*f[31]+0.8944271909999161*vmap[0]*bmag[3]*vmap[3]*f[26]+0.8944271909999161*vmap[0]*bmag[3]*vmap[3]*f[25]+0.8944271909999161*vmap[1]*vmap[2]*bmag[3]*f[22]+0.8944271909999161*vmap[1]*vmap[2]*bmag[3]*f[21]+0.8*vmap[0]*vmap[2]*bmag[6]*f[20]+0.8944271909999161*vmap[0]*bmag[2]*vmap[2]*f[20]+0.8*vmap[0]*vmap[2]*bmag[7]*f[19]+0.8944271909999161*vmap[0]*bmag[1]*vmap[2]*f[19]+0.8944271909999161*vmap[1]*vmap[3]*bmag[7]*f[18]+bmag[1]*vmap[1]*vmap[3]*f[18]+0.8944271909999161*vmap[1]*vmap[3]*bmag[6]*f[17]+vmap[1]*bmag[2]*vmap[3]*f[17]+0.8944271909999159*vmap[0]*vmap[3]*bmag[5]*f[16]+0.8944271909999159*vmap[0]*vmap[3]*bmag[4]*f[16]+bmag[0]*vmap[0]*vmap[3]*f[16]+0.8944271909999159*vmap[1]*vmap[2]*bmag[5]*f[15]+0.8944271909999159*vmap[1]*vmap[2]*bmag[4]*f[15]+bmag[0]*vmap[1]*vmap[2]*f[15]+0.8944271909999159*vmap[0]*vmap[2]*bmag[3]*f[12]+0.8944271909999159*vmap[0]*vmap[2]*bmag[3]*f[11]+vmap[1]*bmag[3]*vmap[3]*f[10]+0.8944271909999161*vmap[0]*vmap[3]*bmag[7]*f[9]+vmap[0]*bmag[1]*vmap[3]*f[9]+0.8944271909999161*vmap[0]*vmap[3]*bmag[6]*f[8]+vmap[0]*bmag[2]*vmap[3]*f[8]+0.8944271909999161*vmap[1]*vmap[2]*bmag[7]*f[7]+bmag[1]*vmap[1]*vmap[2]*f[7]+0.8944271909999161*vmap[0]*f[2]*vmap[2]*bmag[7]+0.8944271909999161*vmap[1]*vmap[2]*bmag[6]*f[6]+vmap[1]*bmag[2]*vmap[2]*f[6]+0.8944271909999161*vmap[0]*f[1]*vmap[2]*bmag[6]+0.8944271909999159*vmap[0]*vmap[2]*bmag[5]*f[5]+0.8944271909999159*vmap[0]*vmap[2]*bmag[4]*f[5]+bmag[0]*vmap[0]*vmap[2]*f[5]+vmap[0]*bmag[3]*vmap[3]*f[4]+vmap[1]*vmap[2]*bmag[3]*f[3]+f[0]*vmap[0]*vmap[2]*bmag[3]+vmap[0]*bmag[1]*f[2]*vmap[2]+vmap[0]*f[1]*bmag[2]*vmap[2])*volFact)/m_+(1.8973665961010278*vmap[0]*vmap1R2*f[34]+1.2727922061357855*vmap1R3*f[15]+2.1213203435596424*vmap0R2*vmap[1]*f[15]+2.1213203435596424*vmap[0]*vmap1R2*f[5]+0.7071067811865475*vmap0R3*f[5])*volFact; 
  out[4] += ((0.8944271909999159*vmap[1]*vmap[3]*bmag[7]*f[45]+0.6388765649999399*vmap[1]*vmap[3]*bmag[6]*f[44]+1.0000000000000002*vmap[1]*bmag[2]*vmap[3]*f[44]+0.6388765649999399*vmap[1]*vmap[3]*bmag[4]*f[37]+bmag[0]*vmap[1]*vmap[3]*f[37]+0.8944271909999161*vmap[0]*vmap[3]*bmag[7]*f[36]+0.63887656499994*vmap[0]*vmap[3]*bmag[6]*f[35]+vmap[0]*bmag[2]*vmap[3]*f[35]+0.8944271909999161*vmap[1]*vmap[2]*bmag[7]*f[33]+0.63887656499994*vmap[1]*vmap[2]*bmag[6]*f[32]+vmap[1]*bmag[2]*vmap[2]*f[32]+0.8944271909999159*vmap[1]*bmag[3]*vmap[3]*f[31]+0.63887656499994*vmap[0]*vmap[3]*bmag[4]*f[25]+1.0000000000000002*bmag[0]*vmap[0]*vmap[3]*f[25]+0.63887656499994*vmap[1]*vmap[2]*bmag[4]*f[21]+1.0000000000000002*bmag[0]*vmap[1]*vmap[2]*f[21]+0.8944271909999159*vmap[0]*vmap[2]*bmag[7]*f[20]+0.6388765649999399*vmap[0]*vmap[2]*bmag[6]*f[19]+1.0000000000000002*vmap[0]*bmag[2]*vmap[2]*f[19]+1.0000000000000002*vmap[1]*vmap[3]*bmag[6]*f[18]+0.8944271909999159*bmag[1]*vmap[1]*vmap[3]*f[17]+0.8944271909999159*vmap[0]*bmag[3]*vmap[3]*f[16]+0.8944271909999159*vmap[1]*vmap[2]*bmag[3]*f[15]+0.6388765649999399*vmap[0]*vmap[2]*bmag[4]*f[11]+bmag[0]*vmap[0]*vmap[2]*f[11]+vmap[1]*vmap[3]*bmag[4]*f[10]+1.0000000000000002*vmap[0]*vmap[3]*bmag[6]*f[9]+0.8944271909999159*vmap[0]*bmag[1]*vmap[3]*f[8]+1.0000000000000002*vmap[1]*vmap[2]*bmag[6]*f[7]+0.8944271909999159*bmag[1]*vmap[1]*vmap[2]*f[6]+1.0000000000000002*vmap[0]*f[2]*vmap[2]*bmag[6]+0.8944271909999159*vmap[0]*vmap[2]*bmag[3]*f[5]+vmap[0]*vmap[3]*bmag[4]*f[4]+vmap[1]*vmap[2]*f[3]*bmag[4]+f[0]*vmap[0]*vmap[2]*bmag[4]+0.8944271909999159*vmap[0]*bmag[1]*f[1]*vmap[2])*volFact)/m_+(1.2727922061357853*vmap1R3*f[21]+2.1213203435596424*vmap0R2*vmap[1]*f[21]+2.1213203435596424*vmap[0]*vmap1R2*f[11]+0.7071067811865475*vmap0R3*f[11])*volFact; 
  out[5] += ((0.6388765649999399*vmap[1]*vmap[3]*bmag[7]*f[45]+1.0000000000000002*bmag[1]*vmap[1]*vmap[3]*f[45]+0.8944271909999159*vmap[1]*vmap[3]*bmag[6]*f[44]+0.6388765649999399*vmap[1]*vmap[3]*bmag[5]*f[38]+bmag[0]*vmap[1]*vmap[3]*f[38]+0.63887656499994*vmap[0]*vmap[3]*bmag[7]*f[36]+vmap[0]*bmag[1]*vmap[3]*f[36]+0.8944271909999161*vmap[0]*vmap[3]*bmag[6]*f[35]+0.63887656499994*vmap[1]*vmap[2]*bmag[7]*f[33]+bmag[1]*vmap[1]*vmap[2]*f[33]+0.8944271909999161*vmap[1]*vmap[2]*bmag[6]*f[32]+0.8944271909999159*vmap[1]*bmag[3]*vmap[3]*f[31]+0.63887656499994*vmap[0]*vmap[3]*bmag[5]*f[26]+1.0000000000000002*bmag[0]*vmap[0]*vmap[3]*f[26]+0.63887656499994*vmap[1]*vmap[2]*bmag[5]*f[22]+1.0000000000000002*bmag[0]*vmap[1]*vmap[2]*f[22]+0.6388765649999399*vmap[0]*vmap[2]*bmag[7]*f[20]+1.0000000000000002*vmap[0]*bmag[1]*vmap[2]*f[20]+0.8944271909999159*vmap[0]*vmap[2]*bmag[6]*f[19]+0.8944271909999159*vmap[1]*bmag[2]*vmap[3]*f[18]+1.0000000000000002*vmap[1]*vmap[3]*bmag[7]*f[17]+0.8944271909999159*vmap[0]*bmag[3]*vmap[3]*f[16]+0.8944271909999159*vmap[1]*vmap[2]*bmag[3]*f[15]+0.6388765649999399*vmap[0]*vmap[2]*bmag[5]*f[12]+bmag[0]*vmap[0]*vmap[2]*f[12]+vmap[1]*vmap[3]*bmag[5]*f[10]+0.8944271909999159*vmap[0]*bmag[2]*vmap[3]*f[9]+1.0000000000000002*vmap[0]*vmap[3]*bmag[7]*f[8]+0.8944271909999159*vmap[1]*bmag[2]*vmap[2]*f[7]+1.0000000000000002*vmap[1]*vmap[2]*f[6]*bmag[7]+1.0000000000000002*vmap[0]*f[1]*vmap[2]*bmag[7]+0.8944271909999159*vmap[0]*vmap[2]*bmag[3]*f[5]+vmap[0]*vmap[3]*f[4]*bmag[5]+vmap[1]*vmap[2]*f[3]*bmag[5]+f[0]*vmap[0]*vmap[2]*bmag[5]+0.8944271909999159*vmap[0]*bmag[2]*f[2]*vmap[2])*volFact)/m_+(1.2727922061357853*vmap1R3*f[22]+2.1213203435596424*vmap0R2*vmap[1]*f[22]+2.1213203435596424*vmap[0]*vmap1R2*f[12]+0.7071067811865475*vmap0R3*f[12])*volFact; 
  out[6] += ((0.8*vmap[1]*bmag[3]*vmap[3]*f[45]+0.8944271909999159*vmap[1]*vmap[3]*bmag[5]*f[44]+0.6388765649999399*vmap[1]*vmap[3]*bmag[4]*f[44]+bmag[0]*vmap[1]*vmap[3]*f[44]+0.8944271909999159*vmap[1]*vmap[3]*bmag[6]*f[38]+0.6388765649999399*vmap[1]*vmap[3]*bmag[6]*f[37]+1.0000000000000002*vmap[1]*bmag[2]*vmap[3]*f[37]+0.8*vmap[0]*bmag[3]*vmap[3]*f[36]+0.8944271909999161*vmap[0]*vmap[3]*bmag[5]*f[35]+0.63887656499994*vmap[0]*vmap[3]*bmag[4]*f[35]+1.0000000000000002*bmag[0]*vmap[0]*vmap[3]*f[35]+0.8*vmap[1]*vmap[2]*bmag[3]*f[33]+0.8944271909999161*vmap[1]*vmap[2]*bmag[5]*f[32]+0.63887656499994*vmap[1]*vmap[2]*bmag[4]*f[32]+1.0000000000000002*bmag[0]*vmap[1]*vmap[2]*f[32]+0.8*vmap[1]*vmap[3]*bmag[7]*f[31]+0.8944271909999161*bmag[1]*vmap[1]*vmap[3]*f[31]+0.8944271909999161*vmap[0]*vmap[3]*bmag[6]*f[26]+0.63887656499994*vmap[0]*vmap[3]*bmag[6]*f[25]+vmap[0]*bmag[2]*vmap[3]*f[25]+0.8944271909999161*vmap[1]*vmap[2]*bmag[6]*f[22]+0.63887656499994*vmap[1]*vmap[2]*bmag[6]*f[21]+vmap[1]*bmag[2]*vmap[2]*f[21]+0.8*vmap[0]*vmap[2]*bmag[3]*f[20]+0.8944271909999159*vmap[0]*vmap[2]*bmag[5]*f[19]+0.6388765649999399*vmap[0]*vmap[2]*bmag[4]*f[19]+bmag[0]*vmap[0]*vmap[2]*f[19]+1.0000000000000002*vmap[1]*vmap[3]*bmag[4]*f[18]+0.8944271909999161*vmap[1]*bmag[3]*vmap[3]*f[17]+0.8*vmap[0]*vmap[3]*bmag[7]*f[16]+0.8944271909999161*vmap[0]*bmag[1]*vmap[3]*f[16]+0.8*vmap[1]*vmap[2]*bmag[7]*f[15]+0.8944271909999161*bmag[1]*vmap[1]*vmap[2]*f[15]+0.8944271909999159*vmap[0]*vmap[2]*bmag[6]*f[12]+0.6388765649999399*vmap[0]*vmap[2]*bmag[6]*f[11]+1.0000000000000002*vmap[0]*bmag[2]*vmap[2]*f[11]+vmap[1]*vmap[3]*bmag[6]*f[10]+1.0000000000000002*vmap[0]*vmap[3]*bmag[4]*f[9]+0.8944271909999161*vmap[0]*bmag[3]*vmap[3]*f[8]+1.0000000000000002*vmap[1]*vmap[2]*bmag[4]*f[7]+0.8*vmap[0]*vmap[2]*f[5]*bmag[7]+0.8944271909999161*vmap[1]*vmap[2]*bmag[3]*f[6]+vmap[0]*vmap[3]*f[4]*bmag[6]+vmap[1]*vmap[2]*f[3]*bmag[6]+f[0]*vmap[0]*vmap[2]*bmag[6]+0.8944271909999161*vmap[0]*bmag[1]*vmap[2]*f[5]+1.0000000000000002*vmap[0]*f[2]*vmap[2]*bmag[4]+0.8944271909999161*vmap[0]*f[1]*vmap[2]*bmag[3])*volFact)/m_+(1.2727922061357853*vmap1R3*f[32]+2.1213203435596424*vmap0R2*vmap[1]*f[32]+2.1213203435596424*vmap[0]*vmap1R2*f[19]+0.7071067811865475*vmap0R3*f[19])*volFact; 
  out[7] += ((0.6388765649999399*vmap[1]*vmap[3]*bmag[5]*f[45]+0.8944271909999159*vmap[1]*vmap[3]*bmag[4]*f[45]+bmag[0]*vmap[1]*vmap[3]*f[45]+0.8*vmap[1]*bmag[3]*vmap[3]*f[44]+0.6388765649999399*vmap[1]*vmap[3]*bmag[7]*f[38]+1.0000000000000002*bmag[1]*vmap[1]*vmap[3]*f[38]+0.8944271909999159*vmap[1]*vmap[3]*bmag[7]*f[37]+0.63887656499994*vmap[0]*vmap[3]*bmag[5]*f[36]+0.8944271909999161*vmap[0]*vmap[3]*bmag[4]*f[36]+1.0000000000000002*bmag[0]*vmap[0]*vmap[3]*f[36]+0.8*vmap[0]*bmag[3]*vmap[3]*f[35]+0.63887656499994*vmap[1]*vmap[2]*bmag[5]*f[33]+0.8944271909999161*vmap[1]*vmap[2]*bmag[4]*f[33]+1.0000000000000002*bmag[0]*vmap[1]*vmap[2]*f[33]+0.8*vmap[1]*vmap[2]*bmag[3]*f[32]+0.8*vmap[1]*vmap[3]*bmag[6]*f[31]+0.8944271909999161*vmap[1]*bmag[2]*vmap[3]*f[31]+0.63887656499994*vmap[0]*vmap[3]*bmag[7]*f[26]+vmap[0]*bmag[1]*vmap[3]*f[26]+0.8944271909999161*vmap[0]*vmap[3]*bmag[7]*f[25]+0.63887656499994*vmap[1]*vmap[2]*bmag[7]*f[22]+bmag[1]*vmap[1]*vmap[2]*f[22]+0.8944271909999161*vmap[1]*vmap[2]*bmag[7]*f[21]+0.6388765649999399*vmap[0]*vmap[2]*bmag[5]*f[20]+0.8944271909999159*vmap[0]*vmap[2]*bmag[4]*f[20]+bmag[0]*vmap[0]*vmap[2]*f[20]+0.8*vmap[0]*vmap[2]*bmag[3]*f[19]+0.8944271909999161*vmap[1]*bmag[3]*vmap[3]*f[18]+1.0000000000000002*vmap[1]*vmap[3]*bmag[5]*f[17]+0.8*vmap[0]*vmap[3]*bmag[6]*f[16]+0.8944271909999161*vmap[0]*bmag[2]*vmap[3]*f[16]+0.8*vmap[1]*vmap[2]*bmag[6]*f[15]+0.8944271909999161*vmap[1]*bmag[2]*vmap[2]*f[15]+0.6388765649999399*vmap[0]*vmap[2]*bmag[7]*f[12]+1.0000000000000002*vmap[0]*bmag[1]*vmap[2]*f[12]+0.8944271909999159*vmap[0]*vmap[2]*bmag[7]*f[11]+vmap[1]*vmap[3]*bmag[7]*f[10]+0.8944271909999161*vmap[0]*bmag[3]*vmap[3]*f[9]+1.0000000000000002*vmap[0]*vmap[3]*bmag[5]*f[8]+0.8944271909999161*vmap[1]*vmap[2]*bmag[3]*f[7]+vmap[0]*vmap[3]*f[4]*bmag[7]+vmap[1]*vmap[2]*f[3]*bmag[7]+f[0]*vmap[0]*vmap[2]*bmag[7]+1.0000000000000002*vmap[1]*vmap[2]*bmag[5]*f[6]+0.8*vmap[0]*vmap[2]*f[5]*bmag[6]+0.8944271909999161*vmap[0]*bmag[2]*vmap[2]*f[5]+1.0000000000000002*vmap[0]*f[1]*vmap[2]*bmag[5]+0.8944271909999161*vmap[0]*f[2]*vmap[2]*bmag[3])*volFact)/m_+(1.2727922061357853*vmap1R3*f[33]+2.1213203435596424*vmap0R2*vmap[1]*f[33]+2.1213203435596424*vmap[0]*vmap1R2*f[20]+0.7071067811865475*vmap0R3*f[20])*volFact; 
} 

GKYL_CU_DH void mom_gyrokinetic_M2PAR_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[2]*dxv[3]/m_; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (0.8944271909999159*vmap1R2*f[13]+2.0*vmap[0]*vmap[1]*f[3]+f[0]*vmap1R2+f[0]*vmap0R2)*volFact; 
  out[1] += (0.8944271909999161*vmap1R2*f[23]+2.0*vmap[0]*vmap[1]*f[6]+f[1]*vmap1R2+vmap0R2*f[1])*volFact; 
  out[2] += (0.8944271909999161*vmap1R2*f[24]+2.0*vmap[0]*vmap[1]*f[7]+vmap1R2*f[2]+vmap0R2*f[2])*volFact; 
  out[3] += (0.8944271909999159*vmap1R2*f[34]+2.0*vmap[0]*vmap[1]*f[15]+vmap1R2*f[5]+vmap0R2*f[5])*volFact; 
  out[4] += (2.0000000000000004*vmap[0]*vmap[1]*f[21]+vmap1R2*f[11]+vmap0R2*f[11])*volFact; 
  out[5] += (2.0000000000000004*vmap[0]*vmap[1]*f[22]+vmap1R2*f[12]+vmap0R2*f[12])*volFact; 
  out[6] += (2.0000000000000004*vmap[0]*vmap[1]*f[32]+vmap1R2*f[19]+vmap0R2*f[19])*volFact; 
  out[7] += (2.0000000000000004*vmap[0]*vmap[1]*f[33]+vmap1R2*f[20]+vmap0R2*f[20])*volFact; 
} 

GKYL_CU_DH void mom_gyrokinetic_M3PAR_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[2]*dxv[3]/m_; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap0R3 = pow(vmap[0],3);
  const double vmap1R2 = pow(vmap[1],2);
  const double vmap1R3 = pow(vmap[1],3);

  out[0] += (1.8973665961010278*vmap[0]*vmap1R2*f[13]+1.2727922061357855*vmap1R3*f[3]+2.1213203435596424*vmap0R2*vmap[1]*f[3]+2.1213203435596424*f[0]*vmap[0]*vmap1R2+0.7071067811865475*f[0]*vmap0R3)*volFact; 
  out[1] += (1.897366596101028*vmap[0]*vmap1R2*f[23]+1.2727922061357855*vmap1R3*f[6]+2.1213203435596424*vmap0R2*vmap[1]*f[6]+2.1213203435596424*vmap[0]*f[1]*vmap1R2+0.7071067811865475*vmap0R3*f[1])*volFact; 
  out[2] += (1.897366596101028*vmap[0]*vmap1R2*f[24]+1.2727922061357855*vmap1R3*f[7]+2.1213203435596424*vmap0R2*vmap[1]*f[7]+2.1213203435596424*vmap[0]*vmap1R2*f[2]+0.7071067811865475*vmap0R3*f[2])*volFact; 
  out[3] += (1.8973665961010278*vmap[0]*vmap1R2*f[34]+1.2727922061357855*vmap1R3*f[15]+2.1213203435596424*vmap0R2*vmap[1]*f[15]+2.1213203435596424*vmap[0]*vmap1R2*f[5]+0.7071067811865475*vmap0R3*f[5])*volFact; 
  out[4] += (1.2727922061357853*vmap1R3*f[21]+2.1213203435596424*vmap0R2*vmap[1]*f[21]+2.1213203435596424*vmap[0]*vmap1R2*f[11]+0.7071067811865475*vmap0R3*f[11])*volFact; 
  out[5] += (1.2727922061357853*vmap1R3*f[22]+2.1213203435596424*vmap0R2*vmap[1]*f[22]+2.1213203435596424*vmap[0]*vmap1R2*f[12]+0.7071067811865475*vmap0R3*f[12])*volFact; 
  out[6] += (1.2727922061357853*vmap1R3*f[32]+2.1213203435596424*vmap0R2*vmap[1]*f[32]+2.1213203435596424*vmap[0]*vmap1R2*f[19]+0.7071067811865475*vmap0R3*f[19])*volFact; 
  out[7] += (1.2727922061357853*vmap1R3*f[33]+2.1213203435596424*vmap0R2*vmap[1]*f[33]+2.1213203435596424*vmap[0]*vmap1R2*f[20]+0.7071067811865475*vmap0R3*f[20])*volFact; 
} 

GKYL_CU_DH void mom_gyrokinetic_M0M1M2_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[2]*dxv[3]/m_; 
  double tmp[8]; 
  tmp[0] = 1.4142135623730951*vmap[3]*f[4]+1.4142135623730951*f[0]*vmap[2]; 
  tmp[1] = 1.4142135623730951*vmap[3]*f[8]+1.4142135623730951*f[1]*vmap[2]; 
  tmp[2] = 1.4142135623730951*vmap[3]*f[9]+1.4142135623730951*f[2]*vmap[2]; 
  tmp[3] = 1.4142135623730951*vmap[3]*f[16]+1.4142135623730951*vmap[2]*f[5]; 
  tmp[4] = 1.4142135623730951*vmap[3]*f[25]+1.4142135623730951*vmap[2]*f[11]; 
  tmp[5] = 1.4142135623730951*vmap[3]*f[26]+1.4142135623730951*vmap[2]*f[12]; 
  tmp[6] = 1.4142135623730951*vmap[3]*f[35]+1.4142135623730951*vmap[2]*f[19]; 
  tmp[7] = 1.4142135623730951*vmap[3]*f[36]+1.4142135623730951*vmap[2]*f[20]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[5]*volFact; 
  out[4] += 2.0*f[11]*volFact; 
  out[5] += 2.0*f[12]*volFact; 
  out[6] += 2.0*f[19]*volFact; 
  out[7] += 2.0*f[20]*volFact; 
  out[8] += (1.4142135623730951*vmap[1]*f[3]+1.4142135623730951*f[0]*vmap[0])*volFact; 
  out[9] += (1.4142135623730951*vmap[1]*f[6]+1.4142135623730951*vmap[0]*f[1])*volFact; 
  out[10] += (1.4142135623730951*vmap[1]*f[7]+1.4142135623730951*vmap[0]*f[2])*volFact; 
  out[11] += (1.4142135623730951*vmap[1]*f[15]+1.4142135623730951*vmap[0]*f[5])*volFact; 
  out[12] += (1.4142135623730951*vmap[1]*f[21]+1.4142135623730951*vmap[0]*f[11])*volFact; 
  out[13] += (1.4142135623730951*vmap[1]*f[22]+1.4142135623730951*vmap[0]*f[12])*volFact; 
  out[14] += (1.4142135623730951*vmap[1]*f[32]+1.4142135623730951*vmap[0]*f[19])*volFact; 
  out[15] += (1.4142135623730951*vmap[1]*f[33]+1.4142135623730951*vmap[0]*f[20])*volFact; 
  out[16] += ((bmag[7]*tmp[7]+bmag[6]*tmp[6]+bmag[5]*tmp[5]+bmag[4]*tmp[4]+bmag[3]*tmp[3]+bmag[2]*tmp[2]+bmag[1]*tmp[1]+bmag[0]*tmp[0])*volFact)/m_+(0.8944271909999159*vmap1R2*f[13]+2.0*vmap[0]*vmap[1]*f[3]+f[0]*vmap1R2+f[0]*vmap0R2)*volFact; 
  out[17] += ((1.0000000000000002*bmag[5]*tmp[7]+1.0000000000000002*tmp[5]*bmag[7]+0.8944271909999161*bmag[3]*tmp[6]+0.8944271909999161*tmp[3]*bmag[6]+0.8944271909999159*bmag[1]*tmp[4]+0.8944271909999159*tmp[1]*bmag[4]+bmag[2]*tmp[3]+tmp[2]*bmag[3]+bmag[0]*tmp[1]+tmp[0]*bmag[1])*volFact)/m_+(0.8944271909999161*vmap1R2*f[23]+2.0*vmap[0]*vmap[1]*f[6]+f[1]*vmap1R2+vmap0R2*f[1])*volFact; 
  out[18] += ((0.8944271909999161*bmag[3]*tmp[7]+0.8944271909999161*tmp[3]*bmag[7]+1.0000000000000002*bmag[4]*tmp[6]+1.0000000000000002*tmp[4]*bmag[6]+0.8944271909999159*bmag[2]*tmp[5]+0.8944271909999159*tmp[2]*bmag[5]+bmag[1]*tmp[3]+tmp[1]*bmag[3]+bmag[0]*tmp[2]+tmp[0]*bmag[2])*volFact)/m_+(0.8944271909999161*vmap1R2*f[24]+2.0*vmap[0]*vmap[1]*f[7]+vmap1R2*f[2]+vmap0R2*f[2])*volFact; 
  out[19] += ((0.8*bmag[6]*tmp[7]+0.8944271909999161*bmag[2]*tmp[7]+0.8*tmp[6]*bmag[7]+0.8944271909999161*tmp[2]*bmag[7]+0.8944271909999161*bmag[1]*tmp[6]+0.8944271909999161*tmp[1]*bmag[6]+0.8944271909999159*bmag[3]*tmp[5]+0.8944271909999159*tmp[3]*bmag[5]+0.8944271909999159*bmag[3]*tmp[4]+0.8944271909999159*tmp[3]*bmag[4]+bmag[0]*tmp[3]+tmp[0]*bmag[3]+bmag[1]*tmp[2]+tmp[1]*bmag[2])*volFact)/m_+(0.8944271909999159*vmap1R2*f[34]+2.0*vmap[0]*vmap[1]*f[15]+vmap1R2*f[5]+vmap0R2*f[5])*volFact; 
  out[20] += ((0.8944271909999159*bmag[7]*tmp[7]+0.6388765649999399*bmag[6]*tmp[6]+1.0000000000000002*bmag[2]*tmp[6]+1.0000000000000002*tmp[2]*bmag[6]+0.6388765649999399*bmag[4]*tmp[4]+bmag[0]*tmp[4]+tmp[0]*bmag[4]+0.8944271909999159*bmag[3]*tmp[3]+0.8944271909999159*bmag[1]*tmp[1])*volFact)/m_+(2.0000000000000004*vmap[0]*vmap[1]*f[21]+vmap1R2*f[11]+vmap0R2*f[11])*volFact; 
  out[21] += ((0.6388765649999399*bmag[7]*tmp[7]+1.0000000000000002*bmag[1]*tmp[7]+1.0000000000000002*tmp[1]*bmag[7]+0.8944271909999159*bmag[6]*tmp[6]+0.6388765649999399*bmag[5]*tmp[5]+bmag[0]*tmp[5]+tmp[0]*bmag[5]+0.8944271909999159*bmag[3]*tmp[3]+0.8944271909999159*bmag[2]*tmp[2])*volFact)/m_+(2.0000000000000004*vmap[0]*vmap[1]*f[22]+vmap1R2*f[12]+vmap0R2*f[12])*volFact; 
  out[22] += ((0.8*bmag[3]*tmp[7]+0.8*tmp[3]*bmag[7]+0.8944271909999159*bmag[5]*tmp[6]+0.6388765649999399*bmag[4]*tmp[6]+bmag[0]*tmp[6]+0.8944271909999159*tmp[5]*bmag[6]+0.6388765649999399*tmp[4]*bmag[6]+tmp[0]*bmag[6]+1.0000000000000002*bmag[2]*tmp[4]+1.0000000000000002*tmp[2]*bmag[4]+0.8944271909999161*bmag[1]*tmp[3]+0.8944271909999161*tmp[1]*bmag[3])*volFact)/m_+(2.0000000000000004*vmap[0]*vmap[1]*f[32]+vmap1R2*f[19]+vmap0R2*f[19])*volFact; 
  out[23] += ((0.6388765649999399*bmag[5]*tmp[7]+0.8944271909999159*bmag[4]*tmp[7]+bmag[0]*tmp[7]+0.6388765649999399*tmp[5]*bmag[7]+0.8944271909999159*tmp[4]*bmag[7]+tmp[0]*bmag[7]+0.8*bmag[3]*tmp[6]+0.8*tmp[3]*bmag[6]+1.0000000000000002*bmag[1]*tmp[5]+1.0000000000000002*tmp[1]*bmag[5]+0.8944271909999161*bmag[2]*tmp[3]+0.8944271909999161*tmp[2]*bmag[3])*volFact)/m_+(2.0000000000000004*vmap[0]*vmap[1]*f[33]+vmap1R2*f[20]+vmap0R2*f[20])*volFact; 
} 

GKYL_CU_DH void mom_gyrokinetic_M0M1M2M3_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[2]*dxv[3]/m_; 
  double tmp2[8]; 
  tmp2[0] = 1.4142135623730951*vmap[3]*f[4]+1.4142135623730951*f[0]*vmap[2]; 
  tmp2[1] = 1.4142135623730951*vmap[3]*f[8]+1.4142135623730951*f[1]*vmap[2]; 
  tmp2[2] = 1.4142135623730951*vmap[3]*f[9]+1.4142135623730951*f[2]*vmap[2]; 
  tmp2[3] = 1.4142135623730951*vmap[3]*f[16]+1.4142135623730951*vmap[2]*f[5]; 
  tmp2[4] = 1.4142135623730951*vmap[3]*f[25]+1.4142135623730951*vmap[2]*f[11]; 
  tmp2[5] = 1.4142135623730951*vmap[3]*f[26]+1.4142135623730951*vmap[2]*f[12]; 
  tmp2[6] = 1.4142135623730951*vmap[3]*f[35]+1.4142135623730951*vmap[2]*f[19]; 
  tmp2[7] = 1.4142135623730951*vmap[3]*f[36]+1.4142135623730951*vmap[2]*f[20]; 
  double tmp3[8]; 
  tmp3[0] = vmap[1]*vmap[3]*f[10]+vmap[0]*vmap[3]*f[4]+vmap[1]*vmap[2]*f[3]+f[0]*vmap[0]*vmap[2]; 
  tmp3[1] = vmap[1]*vmap[3]*f[17]+vmap[0]*vmap[3]*f[8]+vmap[1]*vmap[2]*f[6]+vmap[0]*f[1]*vmap[2]; 
  tmp3[2] = vmap[1]*vmap[3]*f[18]+vmap[0]*vmap[3]*f[9]+vmap[1]*vmap[2]*f[7]+vmap[0]*f[2]*vmap[2]; 
  tmp3[3] = vmap[1]*vmap[3]*f[31]+vmap[0]*vmap[3]*f[16]+vmap[1]*vmap[2]*f[15]+vmap[0]*vmap[2]*f[5]; 
  tmp3[4] = vmap[1]*vmap[3]*f[37]+1.0000000000000002*vmap[0]*vmap[3]*f[25]+1.0000000000000002*vmap[1]*vmap[2]*f[21]+vmap[0]*vmap[2]*f[11]; 
  tmp3[5] = vmap[1]*vmap[3]*f[38]+1.0000000000000002*vmap[0]*vmap[3]*f[26]+1.0000000000000002*vmap[1]*vmap[2]*f[22]+vmap[0]*vmap[2]*f[12]; 
  tmp3[6] = vmap[1]*vmap[3]*f[44]+1.0000000000000002*vmap[0]*vmap[3]*f[35]+1.0000000000000002*vmap[1]*vmap[2]*f[32]+vmap[0]*vmap[2]*f[19]; 
  tmp3[7] = vmap[1]*vmap[3]*f[45]+1.0000000000000002*vmap[0]*vmap[3]*f[36]+1.0000000000000002*vmap[1]*vmap[2]*f[33]+vmap[0]*vmap[2]*f[20]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap0R3 = pow(vmap[0],3);
  const double vmap1R2 = pow(vmap[1],2);
  const double vmap1R3 = pow(vmap[1],3);

  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[5]*volFact; 
  out[4] += 2.0*f[11]*volFact; 
  out[5] += 2.0*f[12]*volFact; 
  out[6] += 2.0*f[19]*volFact; 
  out[7] += 2.0*f[20]*volFact; 
  out[8] += (1.4142135623730951*vmap[1]*f[3]+1.4142135623730951*f[0]*vmap[0])*volFact; 
  out[9] += (1.4142135623730951*vmap[1]*f[6]+1.4142135623730951*vmap[0]*f[1])*volFact; 
  out[10] += (1.4142135623730951*vmap[1]*f[7]+1.4142135623730951*vmap[0]*f[2])*volFact; 
  out[11] += (1.4142135623730951*vmap[1]*f[15]+1.4142135623730951*vmap[0]*f[5])*volFact; 
  out[12] += (1.4142135623730951*vmap[1]*f[21]+1.4142135623730951*vmap[0]*f[11])*volFact; 
  out[13] += (1.4142135623730951*vmap[1]*f[22]+1.4142135623730951*vmap[0]*f[12])*volFact; 
  out[14] += (1.4142135623730951*vmap[1]*f[32]+1.4142135623730951*vmap[0]*f[19])*volFact; 
  out[15] += (1.4142135623730951*vmap[1]*f[33]+1.4142135623730951*vmap[0]*f[20])*volFact; 
  out[16] += ((bmag[7]*tmp2[7]+bmag[6]*tmp2[6]+bmag[5]*tmp2[5]+bmag[4]*tmp2[4]+bmag[3]*tmp2[3]+bmag[2]*tmp2[2]+bmag[1]*tmp2[1]+bmag[0]*tmp2[0])*volFact)/m_+(0.8944271909999159*vmap1R2*f[13]+2.0*vmap[0]*vmap[1]*f[3]+f[0]*vmap1R2+f[0]*vmap0R2)*volFact; 
  out[17] += ((1.0000000000000002*bmag[5]*tmp2[7]+1.0000000000000002*tmp2[5]*bmag[7]+0.8944271909999161*bmag[3]*tmp2[6]+0.8944271909999161*tmp2[3]*bmag[6]+0.8944271909999159*bmag[1]*tmp2[4]+0.8944271909999159*tmp2[1]*bmag[4]+bmag[2]*tmp2[3]+tmp2[2]*bmag[3]+bmag[0]*tmp2[1]+tmp2[0]*bmag[1])*volFact)/m_+(0.8944271909999161*vmap1R2*f[23]+2.0*vmap[0]*vmap[1]*f[6]+f[1]*vmap1R2+vmap0R2*f[1])*volFact; 
  out[18] += ((0.8944271909999161*bmag[3]*tmp2[7]+0.8944271909999161*tmp2[3]*bmag[7]+1.0000000000000002*bmag[4]*tmp2[6]+1.0000000000000002*tmp2[4]*bmag[6]+0.8944271909999159*bmag[2]*tmp2[5]+0.8944271909999159*tmp2[2]*bmag[5]+bmag[1]*tmp2[3]+tmp2[1]*bmag[3]+bmag[0]*tmp2[2]+tmp2[0]*bmag[2])*volFact)/m_+(0.8944271909999161*vmap1R2*f[24]+2.0*vmap[0]*vmap[1]*f[7]+vmap1R2*f[2]+vmap0R2*f[2])*volFact; 
  out[19] += ((0.8*bmag[6]*tmp2[7]+0.8944271909999161*bmag[2]*tmp2[7]+0.8*tmp2[6]*bmag[7]+0.8944271909999161*tmp2[2]*bmag[7]+0.8944271909999161*bmag[1]*tmp2[6]+0.8944271909999161*tmp2[1]*bmag[6]+0.8944271909999159*bmag[3]*tmp2[5]+0.8944271909999159*tmp2[3]*bmag[5]+0.8944271909999159*bmag[3]*tmp2[4]+0.8944271909999159*tmp2[3]*bmag[4]+bmag[0]*tmp2[3]+tmp2[0]*bmag[3]+bmag[1]*tmp2[2]+tmp2[1]*bmag[2])*volFact)/m_+(0.8944271909999159*vmap1R2*f[34]+2.0*vmap[0]*vmap[1]*f[15]+vmap1R2*f[5]+vmap0R2*f[5])*volFact; 
  out[20] += ((0.8944271909999159*bmag[7]*tmp2[7]+0.6388765649999399*bmag[6]*tmp2[6]+1.0000000000000002*bmag[2]*tmp2[6]+1.0000000000000002*tmp2[2]*bmag[6]+0.6388765649999399*bmag[4]*tmp2[4]+bmag[0]*tmp2[4]+tmp2[0]*bmag[4]+0.8944271909999159*bmag[3]*tmp2[3]+0.8944271909999159*bmag[1]*tmp2[1])*volFact)/m_+(2.0000000000000004*vmap[0]*vmap[1]*f[21]+vmap1R2*f[11]+vmap0R2*f[11])*volFact; 
  out[21] += ((0.6388765649999399*bmag[7]*tmp2[7]+1.0000000000000002*bmag[1]*tmp2[7]+1.0000000000000002*tmp2[1]*bmag[7]+0.8944271909999159*bmag[6]*tmp2[6]+0.6388765649999399*bmag[5]*tmp2[5]+bmag[0]*tmp2[5]+tmp2[0]*bmag[5]+0.8944271909999159*bmag[3]*tmp2[3]+0.8944271909999159*bmag[2]*tmp2[2])*volFact)/m_+(2.0000000000000004*vmap[0]*vmap[1]*f[22]+vmap1R2*f[12]+vmap0R2*f[12])*volFact; 
  out[22] += ((0.8*bmag[3]*tmp2[7]+0.8*tmp2[3]*bmag[7]+0.8944271909999159*bmag[5]*tmp2[6]+0.6388765649999399*bmag[4]*tmp2[6]+bmag[0]*tmp2[6]+0.8944271909999159*tmp2[5]*bmag[6]+0.6388765649999399*tmp2[4]*bmag[6]+tmp2[0]*bmag[6]+1.0000000000000002*bmag[2]*tmp2[4]+1.0000000000000002*tmp2[2]*bmag[4]+0.8944271909999161*bmag[1]*tmp2[3]+0.8944271909999161*tmp2[1]*bmag[3])*volFact)/m_+(2.0000000000000004*vmap[0]*vmap[1]*f[32]+vmap1R2*f[19]+vmap0R2*f[19])*volFact; 
  out[23] += ((0.6388765649999399*bmag[5]*tmp2[7]+0.8944271909999159*bmag[4]*tmp2[7]+bmag[0]*tmp2[7]+0.6388765649999399*tmp2[5]*bmag[7]+0.8944271909999159*tmp2[4]*bmag[7]+tmp2[0]*bmag[7]+0.8*bmag[3]*tmp2[6]+0.8*tmp2[3]*bmag[6]+1.0000000000000002*bmag[1]*tmp2[5]+1.0000000000000002*tmp2[1]*bmag[5]+0.8944271909999161*bmag[2]*tmp2[3]+0.8944271909999161*tmp2[2]*bmag[3])*volFact)/m_+(2.0000000000000004*vmap[0]*vmap[1]*f[33]+vmap1R2*f[20]+vmap0R2*f[20])*volFact; 
  out[24] += ((bmag[7]*tmp3[7]+bmag[6]*tmp3[6]+bmag[5]*tmp3[5]+bmag[4]*tmp3[4]+bmag[3]*tmp3[3]+bmag[2]*tmp3[2]+bmag[1]*tmp3[1]+bmag[0]*tmp3[0])*volFact)/m_+(1.8973665961010278*vmap[0]*vmap1R2*f[13]+1.2727922061357855*vmap1R3*f[3]+2.1213203435596424*vmap0R2*vmap[1]*f[3]+2.1213203435596424*f[0]*vmap[0]*vmap1R2+0.7071067811865475*f[0]*vmap0R3)*volFact; 
  out[25] += ((1.0000000000000002*bmag[5]*tmp3[7]+1.0000000000000002*tmp3[5]*bmag[7]+0.8944271909999161*bmag[3]*tmp3[6]+0.8944271909999161*tmp3[3]*bmag[6]+0.8944271909999159*bmag[1]*tmp3[4]+0.8944271909999159*tmp3[1]*bmag[4]+bmag[2]*tmp3[3]+tmp3[2]*bmag[3]+bmag[0]*tmp3[1]+tmp3[0]*bmag[1])*volFact)/m_+(1.897366596101028*vmap[0]*vmap1R2*f[23]+1.2727922061357855*vmap1R3*f[6]+2.1213203435596424*vmap0R2*vmap[1]*f[6]+2.1213203435596424*vmap[0]*f[1]*vmap1R2+0.7071067811865475*vmap0R3*f[1])*volFact; 
  out[26] += ((0.8944271909999161*bmag[3]*tmp3[7]+0.8944271909999161*tmp3[3]*bmag[7]+1.0000000000000002*bmag[4]*tmp3[6]+1.0000000000000002*tmp3[4]*bmag[6]+0.8944271909999159*bmag[2]*tmp3[5]+0.8944271909999159*tmp3[2]*bmag[5]+bmag[1]*tmp3[3]+tmp3[1]*bmag[3]+bmag[0]*tmp3[2]+tmp3[0]*bmag[2])*volFact)/m_+(1.897366596101028*vmap[0]*vmap1R2*f[24]+1.2727922061357855*vmap1R3*f[7]+2.1213203435596424*vmap0R2*vmap[1]*f[7]+2.1213203435596424*vmap[0]*vmap1R2*f[2]+0.7071067811865475*vmap0R3*f[2])*volFact; 
  out[27] += ((0.8*bmag[6]*tmp3[7]+0.8944271909999161*bmag[2]*tmp3[7]+0.8*tmp3[6]*bmag[7]+0.8944271909999161*tmp3[2]*bmag[7]+0.8944271909999161*bmag[1]*tmp3[6]+0.8944271909999161*tmp3[1]*bmag[6]+0.8944271909999159*bmag[3]*tmp3[5]+0.8944271909999159*tmp3[3]*bmag[5]+0.8944271909999159*bmag[3]*tmp3[4]+0.8944271909999159*tmp3[3]*bmag[4]+bmag[0]*tmp3[3]+tmp3[0]*bmag[3]+bmag[1]*tmp3[2]+tmp3[1]*bmag[2])*volFact)/m_+(1.8973665961010278*vmap[0]*vmap1R2*f[34]+1.2727922061357855*vmap1R3*f[15]+2.1213203435596424*vmap0R2*vmap[1]*f[15]+2.1213203435596424*vmap[0]*vmap1R2*f[5]+0.7071067811865475*vmap0R3*f[5])*volFact; 
  out[28] += ((0.8944271909999159*bmag[7]*tmp3[7]+0.6388765649999399*bmag[6]*tmp3[6]+1.0000000000000002*bmag[2]*tmp3[6]+1.0000000000000002*tmp3[2]*bmag[6]+0.6388765649999399*bmag[4]*tmp3[4]+bmag[0]*tmp3[4]+tmp3[0]*bmag[4]+0.8944271909999159*bmag[3]*tmp3[3]+0.8944271909999159*bmag[1]*tmp3[1])*volFact)/m_+(1.2727922061357853*vmap1R3*f[21]+2.1213203435596424*vmap0R2*vmap[1]*f[21]+2.1213203435596424*vmap[0]*vmap1R2*f[11]+0.7071067811865475*vmap0R3*f[11])*volFact; 
  out[29] += ((0.6388765649999399*bmag[7]*tmp3[7]+1.0000000000000002*bmag[1]*tmp3[7]+1.0000000000000002*tmp3[1]*bmag[7]+0.8944271909999159*bmag[6]*tmp3[6]+0.6388765649999399*bmag[5]*tmp3[5]+bmag[0]*tmp3[5]+tmp3[0]*bmag[5]+0.8944271909999159*bmag[3]*tmp3[3]+0.8944271909999159*bmag[2]*tmp3[2])*volFact)/m_+(1.2727922061357853*vmap1R3*f[22]+2.1213203435596424*vmap0R2*vmap[1]*f[22]+2.1213203435596424*vmap[0]*vmap1R2*f[12]+0.7071067811865475*vmap0R3*f[12])*volFact; 
  out[30] += ((0.8*bmag[3]*tmp3[7]+0.8*tmp3[3]*bmag[7]+0.8944271909999159*bmag[5]*tmp3[6]+0.6388765649999399*bmag[4]*tmp3[6]+bmag[0]*tmp3[6]+0.8944271909999159*tmp3[5]*bmag[6]+0.6388765649999399*tmp3[4]*bmag[6]+tmp3[0]*bmag[6]+1.0000000000000002*bmag[2]*tmp3[4]+1.0000000000000002*tmp3[2]*bmag[4]+0.8944271909999161*bmag[1]*tmp3[3]+0.8944271909999161*tmp3[1]*bmag[3])*volFact)/m_+(1.2727922061357853*vmap1R3*f[32]+2.1213203435596424*vmap0R2*vmap[1]*f[32]+2.1213203435596424*vmap[0]*vmap1R2*f[19]+0.7071067811865475*vmap0R3*f[19])*volFact; 
  out[31] += ((0.6388765649999399*bmag[5]*tmp3[7]+0.8944271909999159*bmag[4]*tmp3[7]+bmag[0]*tmp3[7]+0.6388765649999399*tmp3[5]*bmag[7]+0.8944271909999159*tmp3[4]*bmag[7]+tmp3[0]*bmag[7]+0.8*bmag[3]*tmp3[6]+0.8*tmp3[3]*bmag[6]+1.0000000000000002*bmag[1]*tmp3[5]+1.0000000000000002*tmp3[1]*bmag[5]+0.8944271909999161*bmag[2]*tmp3[3]+0.8944271909999161*tmp3[2]*bmag[3])*volFact)/m_+(1.2727922061357853*vmap1R3*f[33]+2.1213203435596424*vmap0R2*vmap[1]*f[33]+2.1213203435596424*vmap[0]*vmap1R2*f[20]+0.7071067811865475*vmap0R3*f[20])*volFact; 
} 

GKYL_CU_DH void mom_gyrokinetic_HAMILTONIAN_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[2]*dxv[3]/m_; 
  double tmp[8]; 
  tmp[0] = 1.4142135623730951*vmap[3]*f[4]+1.4142135623730951*f[0]*vmap[2]; 
  tmp[1] = 1.4142135623730951*vmap[3]*f[8]+1.4142135623730951*f[1]*vmap[2]; 
  tmp[2] = 1.4142135623730951*vmap[3]*f[9]+1.4142135623730951*f[2]*vmap[2]; 
  tmp[3] = 1.4142135623730951*vmap[3]*f[16]+1.4142135623730951*vmap[2]*f[5]; 
  tmp[4] = 1.4142135623730951*vmap[3]*f[25]+1.4142135623730951*vmap[2]*f[11]; 
  tmp[5] = 1.4142135623730951*vmap[3]*f[26]+1.4142135623730951*vmap[2]*f[12]; 
  tmp[6] = 1.4142135623730951*vmap[3]*f[35]+1.4142135623730951*vmap[2]*f[19]; 
  tmp[7] = 1.4142135623730951*vmap[3]*f[36]+1.4142135623730951*vmap[2]*f[20]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[5]*volFact; 
  out[4] += 2.0*f[11]*volFact; 
  out[5] += 2.0*f[12]*volFact; 
  out[6] += 2.0*f[19]*volFact; 
  out[7] += 2.0*f[20]*volFact; 
  out[8] += (1.4142135623730951*vmap[1]*f[3]+1.4142135623730951*f[0]*vmap[0])*m_*volFact; 
  out[9] += (1.4142135623730951*vmap[1]*f[6]+1.4142135623730951*vmap[0]*f[1])*m_*volFact; 
  out[10] += (1.4142135623730951*vmap[1]*f[7]+1.4142135623730951*vmap[0]*f[2])*m_*volFact; 
  out[11] += (1.4142135623730951*vmap[1]*f[15]+1.4142135623730951*vmap[0]*f[5])*m_*volFact; 
  out[12] += (1.4142135623730951*vmap[1]*f[21]+1.4142135623730951*vmap[0]*f[11])*m_*volFact; 
  out[13] += (1.4142135623730951*vmap[1]*f[22]+1.4142135623730951*vmap[0]*f[12])*m_*volFact; 
  out[14] += (1.4142135623730951*vmap[1]*f[32]+1.4142135623730951*vmap[0]*f[19])*m_*volFact; 
  out[15] += (1.4142135623730951*vmap[1]*f[33]+1.4142135623730951*vmap[0]*f[20])*m_*volFact; 
  out[16] += (phi[7]*f[20]+phi[6]*f[19]+phi[5]*f[12]+phi[4]*f[11]+phi[3]*f[5]+f[2]*phi[2]+f[1]*phi[1]+f[0]*phi[0])*q_*volFact+(0.4472135954999579*vmap1R2*f[13]+vmap[0]*vmap[1]*f[3]+0.5*f[0]*vmap1R2+0.5*f[0]*vmap0R2)*m_*volFact+(0.5*bmag[7]*tmp[7]+0.5*bmag[6]*tmp[6]+0.5*bmag[5]*tmp[5]+0.5*bmag[4]*tmp[4]+0.5*bmag[3]*tmp[3]+0.5*bmag[2]*tmp[2]+0.5*bmag[1]*tmp[1]+0.5*bmag[0]*tmp[0])*volFact; 
  out[17] += (1.0000000000000002*phi[5]*f[20]+0.8944271909999161*phi[3]*f[19]+1.0000000000000002*phi[7]*f[12]+0.8944271909999159*phi[1]*f[11]+0.8944271909999161*f[5]*phi[6]+phi[2]*f[5]+0.8944271909999159*f[1]*phi[4]+f[2]*phi[3]+f[0]*phi[1]+phi[0]*f[1])*q_*volFact+(0.44721359549995804*vmap1R2*f[23]+vmap[0]*vmap[1]*f[6]+0.5*f[1]*vmap1R2+0.5*vmap0R2*f[1])*m_*volFact+(0.5000000000000001*bmag[5]*tmp[7]+0.5000000000000001*tmp[5]*bmag[7]+0.44721359549995804*bmag[3]*tmp[6]+0.44721359549995804*tmp[3]*bmag[6]+0.4472135954999579*bmag[1]*tmp[4]+0.4472135954999579*tmp[1]*bmag[4]+0.5*bmag[2]*tmp[3]+0.5*tmp[2]*bmag[3]+0.5*bmag[0]*tmp[1]+0.5*tmp[0]*bmag[1])*volFact; 
  out[18] += (0.8944271909999161*phi[3]*f[20]+1.0000000000000002*phi[4]*f[19]+0.8944271909999159*phi[2]*f[12]+1.0000000000000002*phi[6]*f[11]+0.8944271909999161*f[5]*phi[7]+0.8944271909999159*f[2]*phi[5]+phi[1]*f[5]+f[1]*phi[3]+f[0]*phi[2]+phi[0]*f[2])*q_*volFact+(0.44721359549995804*vmap1R2*f[24]+vmap[0]*vmap[1]*f[7]+0.5*vmap1R2*f[2]+0.5*vmap0R2*f[2])*m_*volFact+(0.44721359549995804*bmag[3]*tmp[7]+0.44721359549995804*tmp[3]*bmag[7]+0.5000000000000001*bmag[4]*tmp[6]+0.5000000000000001*tmp[4]*bmag[6]+0.4472135954999579*bmag[2]*tmp[5]+0.4472135954999579*tmp[2]*bmag[5]+0.5*bmag[1]*tmp[3]+0.5*tmp[1]*bmag[3]+0.5*bmag[0]*tmp[2]+0.5*tmp[0]*bmag[2])*volFact; 
  out[19] += (0.8*phi[6]*f[20]+0.8944271909999161*phi[2]*f[20]+0.8*phi[7]*f[19]+0.8944271909999161*phi[1]*f[19]+0.8944271909999159*phi[3]*f[12]+0.8944271909999159*phi[3]*f[11]+0.8944271909999161*f[2]*phi[7]+0.8944271909999161*f[1]*phi[6]+0.8944271909999159*f[5]*phi[5]+0.8944271909999159*phi[4]*f[5]+phi[0]*f[5]+f[0]*phi[3]+f[1]*phi[2]+phi[1]*f[2])*q_*volFact+(0.4472135954999579*vmap1R2*f[34]+vmap[0]*vmap[1]*f[15]+0.5*vmap1R2*f[5]+0.5*vmap0R2*f[5])*m_*volFact+(0.4*bmag[6]*tmp[7]+0.44721359549995804*bmag[2]*tmp[7]+0.4*tmp[6]*bmag[7]+0.44721359549995804*tmp[2]*bmag[7]+0.44721359549995804*bmag[1]*tmp[6]+0.44721359549995804*tmp[1]*bmag[6]+0.4472135954999579*bmag[3]*tmp[5]+0.4472135954999579*tmp[3]*bmag[5]+0.4472135954999579*bmag[3]*tmp[4]+0.4472135954999579*tmp[3]*bmag[4]+0.5*bmag[0]*tmp[3]+0.5*tmp[0]*bmag[3]+0.5*bmag[1]*tmp[2]+0.5*tmp[1]*bmag[2])*volFact; 
  out[20] += (0.8944271909999159*phi[7]*f[20]+0.6388765649999399*phi[6]*f[19]+1.0000000000000002*phi[2]*f[19]+0.6388765649999399*phi[4]*f[11]+phi[0]*f[11]+1.0000000000000002*f[2]*phi[6]+0.8944271909999159*phi[3]*f[5]+f[0]*phi[4]+0.8944271909999159*f[1]*phi[1])*q_*volFact+(1.0000000000000002*vmap[0]*vmap[1]*f[21]+0.5*vmap1R2*f[11]+0.5*vmap0R2*f[11])*m_*volFact+(0.4472135954999579*bmag[7]*tmp[7]+0.31943828249996997*bmag[6]*tmp[6]+0.5000000000000001*bmag[2]*tmp[6]+0.5000000000000001*tmp[2]*bmag[6]+0.31943828249996997*bmag[4]*tmp[4]+0.5*bmag[0]*tmp[4]+0.5*tmp[0]*bmag[4]+0.4472135954999579*bmag[3]*tmp[3]+0.4472135954999579*bmag[1]*tmp[1])*volFact; 
  out[21] += (0.6388765649999399*phi[7]*f[20]+1.0000000000000002*phi[1]*f[20]+0.8944271909999159*phi[6]*f[19]+0.6388765649999399*phi[5]*f[12]+phi[0]*f[12]+1.0000000000000002*f[1]*phi[7]+f[0]*phi[5]+0.8944271909999159*phi[3]*f[5]+0.8944271909999159*f[2]*phi[2])*q_*volFact+(1.0000000000000002*vmap[0]*vmap[1]*f[22]+0.5*vmap1R2*f[12]+0.5*vmap0R2*f[12])*m_*volFact+(0.31943828249996997*bmag[7]*tmp[7]+0.5000000000000001*bmag[1]*tmp[7]+0.5000000000000001*tmp[1]*bmag[7]+0.4472135954999579*bmag[6]*tmp[6]+0.31943828249996997*bmag[5]*tmp[5]+0.5*bmag[0]*tmp[5]+0.5*tmp[0]*bmag[5]+0.4472135954999579*bmag[3]*tmp[3]+0.4472135954999579*bmag[2]*tmp[2])*volFact; 
  out[22] += (0.8*phi[3]*f[20]+0.8944271909999159*phi[5]*f[19]+0.6388765649999399*phi[4]*f[19]+phi[0]*f[19]+0.8944271909999159*phi[6]*f[12]+0.6388765649999399*phi[6]*f[11]+1.0000000000000002*phi[2]*f[11]+0.8*f[5]*phi[7]+f[0]*phi[6]+0.8944271909999161*phi[1]*f[5]+1.0000000000000002*f[2]*phi[4]+0.8944271909999161*f[1]*phi[3])*q_*volFact+(1.0000000000000002*vmap[0]*vmap[1]*f[32]+0.5*vmap1R2*f[19]+0.5*vmap0R2*f[19])*m_*volFact+(0.4*bmag[3]*tmp[7]+0.4*tmp[3]*bmag[7]+0.4472135954999579*bmag[5]*tmp[6]+0.31943828249996997*bmag[4]*tmp[6]+0.5*bmag[0]*tmp[6]+0.4472135954999579*tmp[5]*bmag[6]+0.31943828249996997*tmp[4]*bmag[6]+0.5*tmp[0]*bmag[6]+0.5000000000000001*bmag[2]*tmp[4]+0.5000000000000001*tmp[2]*bmag[4]+0.44721359549995804*bmag[1]*tmp[3]+0.44721359549995804*tmp[1]*bmag[3])*volFact; 
  out[23] += (0.6388765649999399*phi[5]*f[20]+0.8944271909999159*phi[4]*f[20]+phi[0]*f[20]+0.8*phi[3]*f[19]+0.6388765649999399*phi[7]*f[12]+1.0000000000000002*phi[1]*f[12]+0.8944271909999159*phi[7]*f[11]+f[0]*phi[7]+0.8*f[5]*phi[6]+1.0000000000000002*f[1]*phi[5]+0.8944271909999161*phi[2]*f[5]+0.8944271909999161*f[2]*phi[3])*q_*volFact+(1.0000000000000002*vmap[0]*vmap[1]*f[33]+0.5*vmap1R2*f[20]+0.5*vmap0R2*f[20])*m_*volFact+(0.31943828249996997*bmag[5]*tmp[7]+0.4472135954999579*bmag[4]*tmp[7]+0.5*bmag[0]*tmp[7]+0.31943828249996997*tmp[5]*bmag[7]+0.4472135954999579*tmp[4]*bmag[7]+0.5*tmp[0]*bmag[7]+0.4*bmag[3]*tmp[6]+0.4*tmp[3]*bmag[6]+0.5000000000000001*bmag[1]*tmp[5]+0.5000000000000001*tmp[1]*bmag[5]+0.44721359549995804*bmag[2]*tmp[3]+0.44721359549995804*tmp[2]*bmag[3])*volFact; 
} 

GKYL_CU_DH void mom_gyrokinetic_M2PERP_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[2]*dxv[3]/m_; 
  double tmp[8]; 
  tmp[0] = 1.4142135623730951*vmap[3]*f[4]+1.4142135623730951*f[0]*vmap[2]; 
  tmp[1] = 1.4142135623730951*vmap[3]*f[8]+1.4142135623730951*f[1]*vmap[2]; 
  tmp[2] = 1.4142135623730951*vmap[3]*f[9]+1.4142135623730951*f[2]*vmap[2]; 
  tmp[3] = 1.4142135623730951*vmap[3]*f[16]+1.4142135623730951*vmap[2]*f[5]; 
  tmp[4] = 1.4142135623730951*vmap[3]*f[25]+1.4142135623730951*vmap[2]*f[11]; 
  tmp[5] = 1.4142135623730951*vmap[3]*f[26]+1.4142135623730951*vmap[2]*f[12]; 
  tmp[6] = 1.4142135623730951*vmap[3]*f[35]+1.4142135623730951*vmap[2]*f[19]; 
  tmp[7] = 1.4142135623730951*vmap[3]*f[36]+1.4142135623730951*vmap[2]*f[20]; 
  out[0] += (2.0*(0.5*bmag[7]*tmp[7]+0.5*bmag[6]*tmp[6]+0.5*bmag[5]*tmp[5]+0.5*bmag[4]*tmp[4]+0.5*bmag[3]*tmp[3]+0.5*bmag[2]*tmp[2]+0.5*bmag[1]*tmp[1]+0.5*bmag[0]*tmp[0])*volFact)/m_; 
  out[1] += (2.0*(0.5000000000000001*bmag[5]*tmp[7]+0.5000000000000001*tmp[5]*bmag[7]+0.44721359549995804*bmag[3]*tmp[6]+0.44721359549995804*tmp[3]*bmag[6]+0.4472135954999579*bmag[1]*tmp[4]+0.4472135954999579*tmp[1]*bmag[4]+0.5*bmag[2]*tmp[3]+0.5*tmp[2]*bmag[3]+0.5*bmag[0]*tmp[1]+0.5*tmp[0]*bmag[1])*volFact)/m_; 
  out[2] += (2.0*(0.44721359549995804*bmag[3]*tmp[7]+0.44721359549995804*tmp[3]*bmag[7]+0.5000000000000001*bmag[4]*tmp[6]+0.5000000000000001*tmp[4]*bmag[6]+0.4472135954999579*bmag[2]*tmp[5]+0.4472135954999579*tmp[2]*bmag[5]+0.5*bmag[1]*tmp[3]+0.5*tmp[1]*bmag[3]+0.5*bmag[0]*tmp[2]+0.5*tmp[0]*bmag[2])*volFact)/m_; 
  out[3] += (2.0*(0.4*bmag[6]*tmp[7]+0.44721359549995804*bmag[2]*tmp[7]+0.4*tmp[6]*bmag[7]+0.44721359549995804*tmp[2]*bmag[7]+0.44721359549995804*bmag[1]*tmp[6]+0.44721359549995804*tmp[1]*bmag[6]+0.4472135954999579*bmag[3]*tmp[5]+0.4472135954999579*tmp[3]*bmag[5]+0.4472135954999579*bmag[3]*tmp[4]+0.4472135954999579*tmp[3]*bmag[4]+0.5*bmag[0]*tmp[3]+0.5*tmp[0]*bmag[3]+0.5*bmag[1]*tmp[2]+0.5*tmp[1]*bmag[2])*volFact)/m_; 
  out[4] += (2.0*(0.4472135954999579*bmag[7]*tmp[7]+0.31943828249996997*bmag[6]*tmp[6]+0.5000000000000001*bmag[2]*tmp[6]+0.5000000000000001*tmp[2]*bmag[6]+0.31943828249996997*bmag[4]*tmp[4]+0.5*bmag[0]*tmp[4]+0.5*tmp[0]*bmag[4]+0.4472135954999579*bmag[3]*tmp[3]+0.4472135954999579*bmag[1]*tmp[1])*volFact)/m_; 
  out[5] += (2.0*(0.31943828249996997*bmag[7]*tmp[7]+0.5000000000000001*bmag[1]*tmp[7]+0.5000000000000001*tmp[1]*bmag[7]+0.4472135954999579*bmag[6]*tmp[6]+0.31943828249996997*bmag[5]*tmp[5]+0.5*bmag[0]*tmp[5]+0.5*tmp[0]*bmag[5]+0.4472135954999579*bmag[3]*tmp[3]+0.4472135954999579*bmag[2]*tmp[2])*volFact)/m_; 
  out[6] += (2.0*(0.4*bmag[3]*tmp[7]+0.4*tmp[3]*bmag[7]+0.4472135954999579*bmag[5]*tmp[6]+0.31943828249996997*bmag[4]*tmp[6]+0.5*bmag[0]*tmp[6]+0.4472135954999579*tmp[5]*bmag[6]+0.31943828249996997*tmp[4]*bmag[6]+0.5*tmp[0]*bmag[6]+0.5000000000000001*bmag[2]*tmp[4]+0.5000000000000001*tmp[2]*bmag[4]+0.44721359549995804*bmag[1]*tmp[3]+0.44721359549995804*tmp[1]*bmag[3])*volFact)/m_; 
  out[7] += (2.0*(0.31943828249996997*bmag[5]*tmp[7]+0.4472135954999579*bmag[4]*tmp[7]+0.5*bmag[0]*tmp[7]+0.31943828249996997*tmp[5]*bmag[7]+0.4472135954999579*tmp[4]*bmag[7]+0.5*tmp[0]*bmag[7]+0.4*bmag[3]*tmp[6]+0.4*tmp[3]*bmag[6]+0.5000000000000001*bmag[1]*tmp[5]+0.5000000000000001*tmp[1]*bmag[5]+0.44721359549995804*bmag[2]*tmp[3]+0.44721359549995804*tmp[2]*bmag[3])*volFact)/m_; 
} 

GKYL_CU_DH void mom_gyrokinetic_M3PERP_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[2]*dxv[3]/m_; 

  out[0] += ((vmap[1]*vmap[3]*bmag[7]*f[45]+vmap[1]*vmap[3]*bmag[6]*f[44]+vmap[1]*vmap[3]*bmag[5]*f[38]+vmap[1]*vmap[3]*bmag[4]*f[37]+1.0000000000000002*vmap[0]*vmap[3]*bmag[7]*f[36]+1.0000000000000002*vmap[0]*vmap[3]*bmag[6]*f[35]+1.0000000000000002*vmap[1]*vmap[2]*bmag[7]*f[33]+1.0000000000000002*vmap[1]*vmap[2]*bmag[6]*f[32]+vmap[1]*bmag[3]*vmap[3]*f[31]+1.0000000000000002*vmap[0]*vmap[3]*bmag[5]*f[26]+1.0000000000000002*vmap[0]*vmap[3]*bmag[4]*f[25]+1.0000000000000002*vmap[1]*vmap[2]*bmag[5]*f[22]+1.0000000000000002*vmap[1]*vmap[2]*bmag[4]*f[21]+vmap[0]*vmap[2]*bmag[7]*f[20]+vmap[0]*vmap[2]*bmag[6]*f[19]+vmap[1]*bmag[2]*vmap[3]*f[18]+bmag[1]*vmap[1]*vmap[3]*f[17]+vmap[0]*bmag[3]*vmap[3]*f[16]+vmap[1]*vmap[2]*bmag[3]*f[15]+vmap[0]*vmap[2]*bmag[5]*f[12]+vmap[0]*vmap[2]*bmag[4]*f[11]+bmag[0]*vmap[1]*vmap[3]*f[10]+vmap[0]*bmag[2]*vmap[3]*f[9]+vmap[0]*bmag[1]*vmap[3]*f[8]+vmap[1]*bmag[2]*vmap[2]*f[7]+bmag[1]*vmap[1]*vmap[2]*f[6]+vmap[0]*vmap[2]*bmag[3]*f[5]+bmag[0]*vmap[0]*vmap[3]*f[4]+bmag[0]*vmap[1]*vmap[2]*f[3]+vmap[0]*bmag[2]*f[2]*vmap[2]+vmap[0]*bmag[1]*f[1]*vmap[2]+bmag[0]*f[0]*vmap[0]*vmap[2])*volFact)/m_; 
  out[1] += ((1.0000000000000002*vmap[1]*vmap[3]*bmag[5]*f[45]+0.8944271909999161*vmap[1]*bmag[3]*vmap[3]*f[44]+1.0000000000000002*vmap[1]*vmap[3]*bmag[7]*f[38]+0.8944271909999159*bmag[1]*vmap[1]*vmap[3]*f[37]+vmap[0]*vmap[3]*bmag[5]*f[36]+0.8944271909999159*vmap[0]*bmag[3]*vmap[3]*f[35]+vmap[1]*vmap[2]*bmag[5]*f[33]+0.8944271909999159*vmap[1]*vmap[2]*bmag[3]*f[32]+0.8944271909999161*vmap[1]*vmap[3]*bmag[6]*f[31]+vmap[1]*bmag[2]*vmap[3]*f[31]+vmap[0]*vmap[3]*bmag[7]*f[26]+0.8944271909999161*vmap[0]*bmag[1]*vmap[3]*f[25]+vmap[1]*vmap[2]*bmag[7]*f[22]+0.8944271909999161*bmag[1]*vmap[1]*vmap[2]*f[21]+1.0000000000000002*vmap[0]*vmap[2]*bmag[5]*f[20]+0.8944271909999161*vmap[0]*vmap[2]*bmag[3]*f[19]+vmap[1]*bmag[3]*vmap[3]*f[18]+0.8944271909999159*vmap[1]*vmap[3]*bmag[4]*f[17]+bmag[0]*vmap[1]*vmap[3]*f[17]+0.8944271909999161*vmap[0]*vmap[3]*bmag[6]*f[16]+vmap[0]*bmag[2]*vmap[3]*f[16]+0.8944271909999161*vmap[1]*vmap[2]*bmag[6]*f[15]+vmap[1]*bmag[2]*vmap[2]*f[15]+1.0000000000000002*vmap[0]*vmap[2]*bmag[7]*f[12]+0.8944271909999159*vmap[0]*bmag[1]*vmap[2]*f[11]+bmag[1]*vmap[1]*vmap[3]*f[10]+vmap[0]*bmag[3]*vmap[3]*f[9]+0.8944271909999159*vmap[0]*vmap[3]*bmag[4]*f[8]+bmag[0]*vmap[0]*vmap[3]*f[8]+vmap[1]*vmap[2]*bmag[3]*f[7]+0.8944271909999159*vmap[1]*vmap[2]*bmag[4]*f[6]+bmag[0]*vmap[1]*vmap[2]*f[6]+0.8944271909999161*vmap[0]*vmap[2]*f[5]*bmag[6]+vmap[0]*bmag[2]*vmap[2]*f[5]+vmap[0]*bmag[1]*vmap[3]*f[4]+0.8944271909999159*vmap[0]*f[1]*vmap[2]*bmag[4]+bmag[1]*vmap[1]*vmap[2]*f[3]+vmap[0]*f[2]*vmap[2]*bmag[3]+bmag[0]*vmap[0]*f[1]*vmap[2]+f[0]*vmap[0]*bmag[1]*vmap[2])*volFact)/m_; 
  out[2] += ((0.8944271909999161*vmap[1]*bmag[3]*vmap[3]*f[45]+1.0000000000000002*vmap[1]*vmap[3]*bmag[4]*f[44]+0.8944271909999159*vmap[1]*bmag[2]*vmap[3]*f[38]+1.0000000000000002*vmap[1]*vmap[3]*bmag[6]*f[37]+0.8944271909999159*vmap[0]*bmag[3]*vmap[3]*f[36]+vmap[0]*vmap[3]*bmag[4]*f[35]+0.8944271909999159*vmap[1]*vmap[2]*bmag[3]*f[33]+vmap[1]*vmap[2]*bmag[4]*f[32]+0.8944271909999161*vmap[1]*vmap[3]*bmag[7]*f[31]+bmag[1]*vmap[1]*vmap[3]*f[31]+0.8944271909999161*vmap[0]*bmag[2]*vmap[3]*f[26]+vmap[0]*vmap[3]*bmag[6]*f[25]+0.8944271909999161*vmap[1]*bmag[2]*vmap[2]*f[22]+vmap[1]*vmap[2]*bmag[6]*f[21]+0.8944271909999161*vmap[0]*vmap[2]*bmag[3]*f[20]+1.0000000000000002*vmap[0]*vmap[2]*bmag[4]*f[19]+0.8944271909999159*vmap[1]*vmap[3]*bmag[5]*f[18]+bmag[0]*vmap[1]*vmap[3]*f[18]+vmap[1]*bmag[3]*vmap[3]*f[17]+0.8944271909999161*vmap[0]*vmap[3]*bmag[7]*f[16]+vmap[0]*bmag[1]*vmap[3]*f[16]+0.8944271909999161*vmap[1]*vmap[2]*bmag[7]*f[15]+bmag[1]*vmap[1]*vmap[2]*f[15]+0.8944271909999159*vmap[0]*bmag[2]*vmap[2]*f[12]+1.0000000000000002*vmap[0]*vmap[2]*bmag[6]*f[11]+vmap[1]*bmag[2]*vmap[3]*f[10]+0.8944271909999159*vmap[0]*vmap[3]*bmag[5]*f[9]+bmag[0]*vmap[0]*vmap[3]*f[9]+vmap[0]*bmag[3]*vmap[3]*f[8]+0.8944271909999159*vmap[1]*vmap[2]*bmag[5]*f[7]+bmag[0]*vmap[1]*vmap[2]*f[7]+0.8944271909999161*vmap[0]*vmap[2]*f[5]*bmag[7]+vmap[1]*vmap[2]*bmag[3]*f[6]+vmap[0]*bmag[1]*vmap[2]*f[5]+0.8944271909999159*vmap[0]*f[2]*vmap[2]*bmag[5]+vmap[0]*bmag[2]*vmap[3]*f[4]+vmap[1]*bmag[2]*vmap[2]*f[3]+vmap[0]*f[1]*vmap[2]*bmag[3]+bmag[0]*vmap[0]*f[2]*vmap[2]+f[0]*vmap[0]*bmag[2]*vmap[2])*volFact)/m_; 
  out[3] += ((0.8*vmap[1]*vmap[3]*bmag[6]*f[45]+0.8944271909999161*vmap[1]*bmag[2]*vmap[3]*f[45]+0.8*vmap[1]*vmap[3]*bmag[7]*f[44]+0.8944271909999161*bmag[1]*vmap[1]*vmap[3]*f[44]+0.8944271909999159*vmap[1]*bmag[3]*vmap[3]*f[38]+0.8944271909999159*vmap[1]*bmag[3]*vmap[3]*f[37]+0.8*vmap[0]*vmap[3]*bmag[6]*f[36]+0.8944271909999159*vmap[0]*bmag[2]*vmap[3]*f[36]+0.8*vmap[0]*vmap[3]*bmag[7]*f[35]+0.8944271909999159*vmap[0]*bmag[1]*vmap[3]*f[35]+0.8*vmap[1]*vmap[2]*bmag[6]*f[33]+0.8944271909999159*vmap[1]*bmag[2]*vmap[2]*f[33]+0.8*vmap[1]*vmap[2]*bmag[7]*f[32]+0.8944271909999159*bmag[1]*vmap[1]*vmap[2]*f[32]+0.8944271909999159*vmap[1]*vmap[3]*bmag[5]*f[31]+0.8944271909999159*vmap[1]*vmap[3]*bmag[4]*f[31]+bmag[0]*vmap[1]*vmap[3]*f[31]+0.8944271909999161*vmap[0]*bmag[3]*vmap[3]*f[26]+0.8944271909999161*vmap[0]*bmag[3]*vmap[3]*f[25]+0.8944271909999161*vmap[1]*vmap[2]*bmag[3]*f[22]+0.8944271909999161*vmap[1]*vmap[2]*bmag[3]*f[21]+0.8*vmap[0]*vmap[2]*bmag[6]*f[20]+0.8944271909999161*vmap[0]*bmag[2]*vmap[2]*f[20]+0.8*vmap[0]*vmap[2]*bmag[7]*f[19]+0.8944271909999161*vmap[0]*bmag[1]*vmap[2]*f[19]+0.8944271909999161*vmap[1]*vmap[3]*bmag[7]*f[18]+bmag[1]*vmap[1]*vmap[3]*f[18]+0.8944271909999161*vmap[1]*vmap[3]*bmag[6]*f[17]+vmap[1]*bmag[2]*vmap[3]*f[17]+0.8944271909999159*vmap[0]*vmap[3]*bmag[5]*f[16]+0.8944271909999159*vmap[0]*vmap[3]*bmag[4]*f[16]+bmag[0]*vmap[0]*vmap[3]*f[16]+0.8944271909999159*vmap[1]*vmap[2]*bmag[5]*f[15]+0.8944271909999159*vmap[1]*vmap[2]*bmag[4]*f[15]+bmag[0]*vmap[1]*vmap[2]*f[15]+0.8944271909999159*vmap[0]*vmap[2]*bmag[3]*f[12]+0.8944271909999159*vmap[0]*vmap[2]*bmag[3]*f[11]+vmap[1]*bmag[3]*vmap[3]*f[10]+0.8944271909999161*vmap[0]*vmap[3]*bmag[7]*f[9]+vmap[0]*bmag[1]*vmap[3]*f[9]+0.8944271909999161*vmap[0]*vmap[3]*bmag[6]*f[8]+vmap[0]*bmag[2]*vmap[3]*f[8]+0.8944271909999161*vmap[1]*vmap[2]*bmag[7]*f[7]+bmag[1]*vmap[1]*vmap[2]*f[7]+0.8944271909999161*vmap[0]*f[2]*vmap[2]*bmag[7]+0.8944271909999161*vmap[1]*vmap[2]*bmag[6]*f[6]+vmap[1]*bmag[2]*vmap[2]*f[6]+0.8944271909999161*vmap[0]*f[1]*vmap[2]*bmag[6]+0.8944271909999159*vmap[0]*vmap[2]*bmag[5]*f[5]+0.8944271909999159*vmap[0]*vmap[2]*bmag[4]*f[5]+bmag[0]*vmap[0]*vmap[2]*f[5]+vmap[0]*bmag[3]*vmap[3]*f[4]+vmap[1]*vmap[2]*bmag[3]*f[3]+f[0]*vmap[0]*vmap[2]*bmag[3]+vmap[0]*bmag[1]*f[2]*vmap[2]+vmap[0]*f[1]*bmag[2]*vmap[2])*volFact)/m_; 
  out[4] += ((0.8944271909999159*vmap[1]*vmap[3]*bmag[7]*f[45]+0.6388765649999399*vmap[1]*vmap[3]*bmag[6]*f[44]+1.0000000000000002*vmap[1]*bmag[2]*vmap[3]*f[44]+0.6388765649999399*vmap[1]*vmap[3]*bmag[4]*f[37]+bmag[0]*vmap[1]*vmap[3]*f[37]+0.8944271909999161*vmap[0]*vmap[3]*bmag[7]*f[36]+0.63887656499994*vmap[0]*vmap[3]*bmag[6]*f[35]+vmap[0]*bmag[2]*vmap[3]*f[35]+0.8944271909999161*vmap[1]*vmap[2]*bmag[7]*f[33]+0.63887656499994*vmap[1]*vmap[2]*bmag[6]*f[32]+vmap[1]*bmag[2]*vmap[2]*f[32]+0.8944271909999159*vmap[1]*bmag[3]*vmap[3]*f[31]+0.63887656499994*vmap[0]*vmap[3]*bmag[4]*f[25]+1.0000000000000002*bmag[0]*vmap[0]*vmap[3]*f[25]+0.63887656499994*vmap[1]*vmap[2]*bmag[4]*f[21]+1.0000000000000002*bmag[0]*vmap[1]*vmap[2]*f[21]+0.8944271909999159*vmap[0]*vmap[2]*bmag[7]*f[20]+0.6388765649999399*vmap[0]*vmap[2]*bmag[6]*f[19]+1.0000000000000002*vmap[0]*bmag[2]*vmap[2]*f[19]+1.0000000000000002*vmap[1]*vmap[3]*bmag[6]*f[18]+0.8944271909999159*bmag[1]*vmap[1]*vmap[3]*f[17]+0.8944271909999159*vmap[0]*bmag[3]*vmap[3]*f[16]+0.8944271909999159*vmap[1]*vmap[2]*bmag[3]*f[15]+0.6388765649999399*vmap[0]*vmap[2]*bmag[4]*f[11]+bmag[0]*vmap[0]*vmap[2]*f[11]+vmap[1]*vmap[3]*bmag[4]*f[10]+1.0000000000000002*vmap[0]*vmap[3]*bmag[6]*f[9]+0.8944271909999159*vmap[0]*bmag[1]*vmap[3]*f[8]+1.0000000000000002*vmap[1]*vmap[2]*bmag[6]*f[7]+0.8944271909999159*bmag[1]*vmap[1]*vmap[2]*f[6]+1.0000000000000002*vmap[0]*f[2]*vmap[2]*bmag[6]+0.8944271909999159*vmap[0]*vmap[2]*bmag[3]*f[5]+vmap[0]*vmap[3]*bmag[4]*f[4]+vmap[1]*vmap[2]*f[3]*bmag[4]+f[0]*vmap[0]*vmap[2]*bmag[4]+0.8944271909999159*vmap[0]*bmag[1]*f[1]*vmap[2])*volFact)/m_; 
  out[5] += ((0.6388765649999399*vmap[1]*vmap[3]*bmag[7]*f[45]+1.0000000000000002*bmag[1]*vmap[1]*vmap[3]*f[45]+0.8944271909999159*vmap[1]*vmap[3]*bmag[6]*f[44]+0.6388765649999399*vmap[1]*vmap[3]*bmag[5]*f[38]+bmag[0]*vmap[1]*vmap[3]*f[38]+0.63887656499994*vmap[0]*vmap[3]*bmag[7]*f[36]+vmap[0]*bmag[1]*vmap[3]*f[36]+0.8944271909999161*vmap[0]*vmap[3]*bmag[6]*f[35]+0.63887656499994*vmap[1]*vmap[2]*bmag[7]*f[33]+bmag[1]*vmap[1]*vmap[2]*f[33]+0.8944271909999161*vmap[1]*vmap[2]*bmag[6]*f[32]+0.8944271909999159*vmap[1]*bmag[3]*vmap[3]*f[31]+0.63887656499994*vmap[0]*vmap[3]*bmag[5]*f[26]+1.0000000000000002*bmag[0]*vmap[0]*vmap[3]*f[26]+0.63887656499994*vmap[1]*vmap[2]*bmag[5]*f[22]+1.0000000000000002*bmag[0]*vmap[1]*vmap[2]*f[22]+0.6388765649999399*vmap[0]*vmap[2]*bmag[7]*f[20]+1.0000000000000002*vmap[0]*bmag[1]*vmap[2]*f[20]+0.8944271909999159*vmap[0]*vmap[2]*bmag[6]*f[19]+0.8944271909999159*vmap[1]*bmag[2]*vmap[3]*f[18]+1.0000000000000002*vmap[1]*vmap[3]*bmag[7]*f[17]+0.8944271909999159*vmap[0]*bmag[3]*vmap[3]*f[16]+0.8944271909999159*vmap[1]*vmap[2]*bmag[3]*f[15]+0.6388765649999399*vmap[0]*vmap[2]*bmag[5]*f[12]+bmag[0]*vmap[0]*vmap[2]*f[12]+vmap[1]*vmap[3]*bmag[5]*f[10]+0.8944271909999159*vmap[0]*bmag[2]*vmap[3]*f[9]+1.0000000000000002*vmap[0]*vmap[3]*bmag[7]*f[8]+0.8944271909999159*vmap[1]*bmag[2]*vmap[2]*f[7]+1.0000000000000002*vmap[1]*vmap[2]*f[6]*bmag[7]+1.0000000000000002*vmap[0]*f[1]*vmap[2]*bmag[7]+0.8944271909999159*vmap[0]*vmap[2]*bmag[3]*f[5]+vmap[0]*vmap[3]*f[4]*bmag[5]+vmap[1]*vmap[2]*f[3]*bmag[5]+f[0]*vmap[0]*vmap[2]*bmag[5]+0.8944271909999159*vmap[0]*bmag[2]*f[2]*vmap[2])*volFact)/m_; 
  out[6] += ((0.8*vmap[1]*bmag[3]*vmap[3]*f[45]+0.8944271909999159*vmap[1]*vmap[3]*bmag[5]*f[44]+0.6388765649999399*vmap[1]*vmap[3]*bmag[4]*f[44]+bmag[0]*vmap[1]*vmap[3]*f[44]+0.8944271909999159*vmap[1]*vmap[3]*bmag[6]*f[38]+0.6388765649999399*vmap[1]*vmap[3]*bmag[6]*f[37]+1.0000000000000002*vmap[1]*bmag[2]*vmap[3]*f[37]+0.8*vmap[0]*bmag[3]*vmap[3]*f[36]+0.8944271909999161*vmap[0]*vmap[3]*bmag[5]*f[35]+0.63887656499994*vmap[0]*vmap[3]*bmag[4]*f[35]+1.0000000000000002*bmag[0]*vmap[0]*vmap[3]*f[35]+0.8*vmap[1]*vmap[2]*bmag[3]*f[33]+0.8944271909999161*vmap[1]*vmap[2]*bmag[5]*f[32]+0.63887656499994*vmap[1]*vmap[2]*bmag[4]*f[32]+1.0000000000000002*bmag[0]*vmap[1]*vmap[2]*f[32]+0.8*vmap[1]*vmap[3]*bmag[7]*f[31]+0.8944271909999161*bmag[1]*vmap[1]*vmap[3]*f[31]+0.8944271909999161*vmap[0]*vmap[3]*bmag[6]*f[26]+0.63887656499994*vmap[0]*vmap[3]*bmag[6]*f[25]+vmap[0]*bmag[2]*vmap[3]*f[25]+0.8944271909999161*vmap[1]*vmap[2]*bmag[6]*f[22]+0.63887656499994*vmap[1]*vmap[2]*bmag[6]*f[21]+vmap[1]*bmag[2]*vmap[2]*f[21]+0.8*vmap[0]*vmap[2]*bmag[3]*f[20]+0.8944271909999159*vmap[0]*vmap[2]*bmag[5]*f[19]+0.6388765649999399*vmap[0]*vmap[2]*bmag[4]*f[19]+bmag[0]*vmap[0]*vmap[2]*f[19]+1.0000000000000002*vmap[1]*vmap[3]*bmag[4]*f[18]+0.8944271909999161*vmap[1]*bmag[3]*vmap[3]*f[17]+0.8*vmap[0]*vmap[3]*bmag[7]*f[16]+0.8944271909999161*vmap[0]*bmag[1]*vmap[3]*f[16]+0.8*vmap[1]*vmap[2]*bmag[7]*f[15]+0.8944271909999161*bmag[1]*vmap[1]*vmap[2]*f[15]+0.8944271909999159*vmap[0]*vmap[2]*bmag[6]*f[12]+0.6388765649999399*vmap[0]*vmap[2]*bmag[6]*f[11]+1.0000000000000002*vmap[0]*bmag[2]*vmap[2]*f[11]+vmap[1]*vmap[3]*bmag[6]*f[10]+1.0000000000000002*vmap[0]*vmap[3]*bmag[4]*f[9]+0.8944271909999161*vmap[0]*bmag[3]*vmap[3]*f[8]+1.0000000000000002*vmap[1]*vmap[2]*bmag[4]*f[7]+0.8*vmap[0]*vmap[2]*f[5]*bmag[7]+0.8944271909999161*vmap[1]*vmap[2]*bmag[3]*f[6]+vmap[0]*vmap[3]*f[4]*bmag[6]+vmap[1]*vmap[2]*f[3]*bmag[6]+f[0]*vmap[0]*vmap[2]*bmag[6]+0.8944271909999161*vmap[0]*bmag[1]*vmap[2]*f[5]+1.0000000000000002*vmap[0]*f[2]*vmap[2]*bmag[4]+0.8944271909999161*vmap[0]*f[1]*vmap[2]*bmag[3])*volFact)/m_; 
  out[7] += ((0.6388765649999399*vmap[1]*vmap[3]*bmag[5]*f[45]+0.8944271909999159*vmap[1]*vmap[3]*bmag[4]*f[45]+bmag[0]*vmap[1]*vmap[3]*f[45]+0.8*vmap[1]*bmag[3]*vmap[3]*f[44]+0.6388765649999399*vmap[1]*vmap[3]*bmag[7]*f[38]+1.0000000000000002*bmag[1]*vmap[1]*vmap[3]*f[38]+0.8944271909999159*vmap[1]*vmap[3]*bmag[7]*f[37]+0.63887656499994*vmap[0]*vmap[3]*bmag[5]*f[36]+0.8944271909999161*vmap[0]*vmap[3]*bmag[4]*f[36]+1.0000000000000002*bmag[0]*vmap[0]*vmap[3]*f[36]+0.8*vmap[0]*bmag[3]*vmap[3]*f[35]+0.63887656499994*vmap[1]*vmap[2]*bmag[5]*f[33]+0.8944271909999161*vmap[1]*vmap[2]*bmag[4]*f[33]+1.0000000000000002*bmag[0]*vmap[1]*vmap[2]*f[33]+0.8*vmap[1]*vmap[2]*bmag[3]*f[32]+0.8*vmap[1]*vmap[3]*bmag[6]*f[31]+0.8944271909999161*vmap[1]*bmag[2]*vmap[3]*f[31]+0.63887656499994*vmap[0]*vmap[3]*bmag[7]*f[26]+vmap[0]*bmag[1]*vmap[3]*f[26]+0.8944271909999161*vmap[0]*vmap[3]*bmag[7]*f[25]+0.63887656499994*vmap[1]*vmap[2]*bmag[7]*f[22]+bmag[1]*vmap[1]*vmap[2]*f[22]+0.8944271909999161*vmap[1]*vmap[2]*bmag[7]*f[21]+0.6388765649999399*vmap[0]*vmap[2]*bmag[5]*f[20]+0.8944271909999159*vmap[0]*vmap[2]*bmag[4]*f[20]+bmag[0]*vmap[0]*vmap[2]*f[20]+0.8*vmap[0]*vmap[2]*bmag[3]*f[19]+0.8944271909999161*vmap[1]*bmag[3]*vmap[3]*f[18]+1.0000000000000002*vmap[1]*vmap[3]*bmag[5]*f[17]+0.8*vmap[0]*vmap[3]*bmag[6]*f[16]+0.8944271909999161*vmap[0]*bmag[2]*vmap[3]*f[16]+0.8*vmap[1]*vmap[2]*bmag[6]*f[15]+0.8944271909999161*vmap[1]*bmag[2]*vmap[2]*f[15]+0.6388765649999399*vmap[0]*vmap[2]*bmag[7]*f[12]+1.0000000000000002*vmap[0]*bmag[1]*vmap[2]*f[12]+0.8944271909999159*vmap[0]*vmap[2]*bmag[7]*f[11]+vmap[1]*vmap[3]*bmag[7]*f[10]+0.8944271909999161*vmap[0]*bmag[3]*vmap[3]*f[9]+1.0000000000000002*vmap[0]*vmap[3]*bmag[5]*f[8]+0.8944271909999161*vmap[1]*vmap[2]*bmag[3]*f[7]+vmap[0]*vmap[3]*f[4]*bmag[7]+vmap[1]*vmap[2]*f[3]*bmag[7]+f[0]*vmap[0]*vmap[2]*bmag[7]+1.0000000000000002*vmap[1]*vmap[2]*bmag[5]*f[6]+0.8*vmap[0]*vmap[2]*f[5]*bmag[6]+0.8944271909999161*vmap[0]*bmag[2]*vmap[2]*f[5]+1.0000000000000002*vmap[0]*f[1]*vmap[2]*bmag[5]+0.8944271909999161*vmap[0]*f[2]*vmap[2]*bmag[3])*volFact)/m_; 
} 

GKYL_CU_DH void mom_gyrokinetic_M0M1M2PARM2PERP_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[2]*dxv[3]/m_; 
  double tmp[8]; 
  tmp[0] = 1.4142135623730951*vmap[3]*f[4]+1.4142135623730951*f[0]*vmap[2]; 
  tmp[1] = 1.4142135623730951*vmap[3]*f[8]+1.4142135623730951*f[1]*vmap[2]; 
  tmp[2] = 1.4142135623730951*vmap[3]*f[9]+1.4142135623730951*f[2]*vmap[2]; 
  tmp[3] = 1.4142135623730951*vmap[3]*f[16]+1.4142135623730951*vmap[2]*f[5]; 
  tmp[4] = 1.4142135623730951*vmap[3]*f[25]+1.4142135623730951*vmap[2]*f[11]; 
  tmp[5] = 1.4142135623730951*vmap[3]*f[26]+1.4142135623730951*vmap[2]*f[12]; 
  tmp[6] = 1.4142135623730951*vmap[3]*f[35]+1.4142135623730951*vmap[2]*f[19]; 
  tmp[7] = 1.4142135623730951*vmap[3]*f[36]+1.4142135623730951*vmap[2]*f[20]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[5]*volFact; 
  out[4] += 2.0*f[11]*volFact; 
  out[5] += 2.0*f[12]*volFact; 
  out[6] += 2.0*f[19]*volFact; 
  out[7] += 2.0*f[20]*volFact; 
  out[8] += (1.4142135623730951*vmap[1]*f[3]+1.4142135623730951*f[0]*vmap[0])*volFact; 
  out[9] += (1.4142135623730951*vmap[1]*f[6]+1.4142135623730951*vmap[0]*f[1])*volFact; 
  out[10] += (1.4142135623730951*vmap[1]*f[7]+1.4142135623730951*vmap[0]*f[2])*volFact; 
  out[11] += (1.4142135623730951*vmap[1]*f[15]+1.4142135623730951*vmap[0]*f[5])*volFact; 
  out[12] += (1.4142135623730951*vmap[1]*f[21]+1.4142135623730951*vmap[0]*f[11])*volFact; 
  out[13] += (1.4142135623730951*vmap[1]*f[22]+1.4142135623730951*vmap[0]*f[12])*volFact; 
  out[14] += (1.4142135623730951*vmap[1]*f[32]+1.4142135623730951*vmap[0]*f[19])*volFact; 
  out[15] += (1.4142135623730951*vmap[1]*f[33]+1.4142135623730951*vmap[0]*f[20])*volFact; 
  out[16] += (0.8944271909999159*vmap1R2*f[13]+2.0*vmap[0]*vmap[1]*f[3]+f[0]*vmap1R2+f[0]*vmap0R2)*volFact; 
  out[17] += (0.8944271909999161*vmap1R2*f[23]+2.0*vmap[0]*vmap[1]*f[6]+f[1]*vmap1R2+vmap0R2*f[1])*volFact; 
  out[18] += (0.8944271909999161*vmap1R2*f[24]+2.0*vmap[0]*vmap[1]*f[7]+vmap1R2*f[2]+vmap0R2*f[2])*volFact; 
  out[19] += (0.8944271909999159*vmap1R2*f[34]+2.0*vmap[0]*vmap[1]*f[15]+vmap1R2*f[5]+vmap0R2*f[5])*volFact; 
  out[20] += (2.0000000000000004*vmap[0]*vmap[1]*f[21]+vmap1R2*f[11]+vmap0R2*f[11])*volFact; 
  out[21] += (2.0000000000000004*vmap[0]*vmap[1]*f[22]+vmap1R2*f[12]+vmap0R2*f[12])*volFact; 
  out[22] += (2.0000000000000004*vmap[0]*vmap[1]*f[32]+vmap1R2*f[19]+vmap0R2*f[19])*volFact; 
  out[23] += (2.0000000000000004*vmap[0]*vmap[1]*f[33]+vmap1R2*f[20]+vmap0R2*f[20])*volFact; 
  out[24] += ((bmag[7]*tmp[7]+bmag[6]*tmp[6]+bmag[5]*tmp[5]+bmag[4]*tmp[4]+bmag[3]*tmp[3]+bmag[2]*tmp[2]+bmag[1]*tmp[1]+bmag[0]*tmp[0])*volFact)/m_; 
  out[25] += ((1.0000000000000002*bmag[5]*tmp[7]+1.0000000000000002*tmp[5]*bmag[7]+0.8944271909999161*bmag[3]*tmp[6]+0.8944271909999161*tmp[3]*bmag[6]+0.8944271909999159*bmag[1]*tmp[4]+0.8944271909999159*tmp[1]*bmag[4]+bmag[2]*tmp[3]+tmp[2]*bmag[3]+bmag[0]*tmp[1]+tmp[0]*bmag[1])*volFact)/m_; 
  out[26] += ((0.8944271909999161*bmag[3]*tmp[7]+0.8944271909999161*tmp[3]*bmag[7]+1.0000000000000002*bmag[4]*tmp[6]+1.0000000000000002*tmp[4]*bmag[6]+0.8944271909999159*bmag[2]*tmp[5]+0.8944271909999159*tmp[2]*bmag[5]+bmag[1]*tmp[3]+tmp[1]*bmag[3]+bmag[0]*tmp[2]+tmp[0]*bmag[2])*volFact)/m_; 
  out[27] += ((0.8*bmag[6]*tmp[7]+0.8944271909999161*bmag[2]*tmp[7]+0.8*tmp[6]*bmag[7]+0.8944271909999161*tmp[2]*bmag[7]+0.8944271909999161*bmag[1]*tmp[6]+0.8944271909999161*tmp[1]*bmag[6]+0.8944271909999159*bmag[3]*tmp[5]+0.8944271909999159*tmp[3]*bmag[5]+0.8944271909999159*bmag[3]*tmp[4]+0.8944271909999159*tmp[3]*bmag[4]+bmag[0]*tmp[3]+tmp[0]*bmag[3]+bmag[1]*tmp[2]+tmp[1]*bmag[2])*volFact)/m_; 
  out[28] += ((0.8944271909999159*bmag[7]*tmp[7]+0.6388765649999399*bmag[6]*tmp[6]+1.0000000000000002*bmag[2]*tmp[6]+1.0000000000000002*tmp[2]*bmag[6]+0.6388765649999399*bmag[4]*tmp[4]+bmag[0]*tmp[4]+tmp[0]*bmag[4]+0.8944271909999159*bmag[3]*tmp[3]+0.8944271909999159*bmag[1]*tmp[1])*volFact)/m_; 
  out[29] += ((0.6388765649999399*bmag[7]*tmp[7]+1.0000000000000002*bmag[1]*tmp[7]+1.0000000000000002*tmp[1]*bmag[7]+0.8944271909999159*bmag[6]*tmp[6]+0.6388765649999399*bmag[5]*tmp[5]+bmag[0]*tmp[5]+tmp[0]*bmag[5]+0.8944271909999159*bmag[3]*tmp[3]+0.8944271909999159*bmag[2]*tmp[2])*volFact)/m_; 
  out[30] += ((0.8*bmag[3]*tmp[7]+0.8*tmp[3]*bmag[7]+0.8944271909999159*bmag[5]*tmp[6]+0.6388765649999399*bmag[4]*tmp[6]+bmag[0]*tmp[6]+0.8944271909999159*tmp[5]*bmag[6]+0.6388765649999399*tmp[4]*bmag[6]+tmp[0]*bmag[6]+1.0000000000000002*bmag[2]*tmp[4]+1.0000000000000002*tmp[2]*bmag[4]+0.8944271909999161*bmag[1]*tmp[3]+0.8944271909999161*tmp[1]*bmag[3])*volFact)/m_; 
  out[31] += ((0.6388765649999399*bmag[5]*tmp[7]+0.8944271909999159*bmag[4]*tmp[7]+bmag[0]*tmp[7]+0.6388765649999399*tmp[5]*bmag[7]+0.8944271909999159*tmp[4]*bmag[7]+tmp[0]*bmag[7]+0.8*bmag[3]*tmp[6]+0.8*tmp[3]*bmag[6]+1.0000000000000002*bmag[1]*tmp[5]+1.0000000000000002*tmp[1]*bmag[5]+0.8944271909999161*bmag[2]*tmp[3]+0.8944271909999161*tmp[2]*bmag[3])*volFact)/m_; 
} 

GKYL_CU_DH void mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[2]*dxv[3]/m_; 
  double tmp2[8]; 
  tmp2[0] = 1.4142135623730951*vmap[3]*f[4]+1.4142135623730951*f[0]*vmap[2]; 
  tmp2[1] = 1.4142135623730951*vmap[3]*f[8]+1.4142135623730951*f[1]*vmap[2]; 
  tmp2[2] = 1.4142135623730951*vmap[3]*f[9]+1.4142135623730951*f[2]*vmap[2]; 
  tmp2[3] = 1.4142135623730951*vmap[3]*f[16]+1.4142135623730951*vmap[2]*f[5]; 
  tmp2[4] = 1.4142135623730951*vmap[3]*f[25]+1.4142135623730951*vmap[2]*f[11]; 
  tmp2[5] = 1.4142135623730951*vmap[3]*f[26]+1.4142135623730951*vmap[2]*f[12]; 
  tmp2[6] = 1.4142135623730951*vmap[3]*f[35]+1.4142135623730951*vmap[2]*f[19]; 
  tmp2[7] = 1.4142135623730951*vmap[3]*f[36]+1.4142135623730951*vmap[2]*f[20]; 
  double tmp3[8]; 
  tmp3[0] = vmap[1]*vmap[3]*f[10]+vmap[0]*vmap[3]*f[4]+vmap[1]*vmap[2]*f[3]+f[0]*vmap[0]*vmap[2]; 
  tmp3[1] = vmap[1]*vmap[3]*f[17]+vmap[0]*vmap[3]*f[8]+vmap[1]*vmap[2]*f[6]+vmap[0]*f[1]*vmap[2]; 
  tmp3[2] = vmap[1]*vmap[3]*f[18]+vmap[0]*vmap[3]*f[9]+vmap[1]*vmap[2]*f[7]+vmap[0]*f[2]*vmap[2]; 
  tmp3[3] = vmap[1]*vmap[3]*f[31]+vmap[0]*vmap[3]*f[16]+vmap[1]*vmap[2]*f[15]+vmap[0]*vmap[2]*f[5]; 
  tmp3[4] = vmap[1]*vmap[3]*f[37]+1.0000000000000002*vmap[0]*vmap[3]*f[25]+1.0000000000000002*vmap[1]*vmap[2]*f[21]+vmap[0]*vmap[2]*f[11]; 
  tmp3[5] = vmap[1]*vmap[3]*f[38]+1.0000000000000002*vmap[0]*vmap[3]*f[26]+1.0000000000000002*vmap[1]*vmap[2]*f[22]+vmap[0]*vmap[2]*f[12]; 
  tmp3[6] = vmap[1]*vmap[3]*f[44]+1.0000000000000002*vmap[0]*vmap[3]*f[35]+1.0000000000000002*vmap[1]*vmap[2]*f[32]+vmap[0]*vmap[2]*f[19]; 
  tmp3[7] = vmap[1]*vmap[3]*f[45]+1.0000000000000002*vmap[0]*vmap[3]*f[36]+1.0000000000000002*vmap[1]*vmap[2]*f[33]+vmap[0]*vmap[2]*f[20]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap0R3 = pow(vmap[0],3);
  const double vmap1R2 = pow(vmap[1],2);
  const double vmap1R3 = pow(vmap[1],3);

  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[5]*volFact; 
  out[4] += 2.0*f[11]*volFact; 
  out[5] += 2.0*f[12]*volFact; 
  out[6] += 2.0*f[19]*volFact; 
  out[7] += 2.0*f[20]*volFact; 
  out[8] += (1.4142135623730951*vmap[1]*f[3]+1.4142135623730951*f[0]*vmap[0])*volFact; 
  out[9] += (1.4142135623730951*vmap[1]*f[6]+1.4142135623730951*vmap[0]*f[1])*volFact; 
  out[10] += (1.4142135623730951*vmap[1]*f[7]+1.4142135623730951*vmap[0]*f[2])*volFact; 
  out[11] += (1.4142135623730951*vmap[1]*f[15]+1.4142135623730951*vmap[0]*f[5])*volFact; 
  out[12] += (1.4142135623730951*vmap[1]*f[21]+1.4142135623730951*vmap[0]*f[11])*volFact; 
  out[13] += (1.4142135623730951*vmap[1]*f[22]+1.4142135623730951*vmap[0]*f[12])*volFact; 
  out[14] += (1.4142135623730951*vmap[1]*f[32]+1.4142135623730951*vmap[0]*f[19])*volFact; 
  out[15] += (1.4142135623730951*vmap[1]*f[33]+1.4142135623730951*vmap[0]*f[20])*volFact; 
  out[16] += (0.8944271909999159*vmap1R2*f[13]+2.0*vmap[0]*vmap[1]*f[3]+f[0]*vmap1R2+f[0]*vmap0R2)*volFact; 
  out[17] += (0.8944271909999161*vmap1R2*f[23]+2.0*vmap[0]*vmap[1]*f[6]+f[1]*vmap1R2+vmap0R2*f[1])*volFact; 
  out[18] += (0.8944271909999161*vmap1R2*f[24]+2.0*vmap[0]*vmap[1]*f[7]+vmap1R2*f[2]+vmap0R2*f[2])*volFact; 
  out[19] += (0.8944271909999159*vmap1R2*f[34]+2.0*vmap[0]*vmap[1]*f[15]+vmap1R2*f[5]+vmap0R2*f[5])*volFact; 
  out[20] += (2.0000000000000004*vmap[0]*vmap[1]*f[21]+vmap1R2*f[11]+vmap0R2*f[11])*volFact; 
  out[21] += (2.0000000000000004*vmap[0]*vmap[1]*f[22]+vmap1R2*f[12]+vmap0R2*f[12])*volFact; 
  out[22] += (2.0000000000000004*vmap[0]*vmap[1]*f[32]+vmap1R2*f[19]+vmap0R2*f[19])*volFact; 
  out[23] += (2.0000000000000004*vmap[0]*vmap[1]*f[33]+vmap1R2*f[20]+vmap0R2*f[20])*volFact; 
  out[24] += ((bmag[7]*tmp2[7]+bmag[6]*tmp2[6]+bmag[5]*tmp2[5]+bmag[4]*tmp2[4]+bmag[3]*tmp2[3]+bmag[2]*tmp2[2]+bmag[1]*tmp2[1]+bmag[0]*tmp2[0])*volFact)/m_; 
  out[25] += ((1.0000000000000002*bmag[5]*tmp2[7]+1.0000000000000002*tmp2[5]*bmag[7]+0.8944271909999161*bmag[3]*tmp2[6]+0.8944271909999161*tmp2[3]*bmag[6]+0.8944271909999159*bmag[1]*tmp2[4]+0.8944271909999159*tmp2[1]*bmag[4]+bmag[2]*tmp2[3]+tmp2[2]*bmag[3]+bmag[0]*tmp2[1]+tmp2[0]*bmag[1])*volFact)/m_; 
  out[26] += ((0.8944271909999161*bmag[3]*tmp2[7]+0.8944271909999161*tmp2[3]*bmag[7]+1.0000000000000002*bmag[4]*tmp2[6]+1.0000000000000002*tmp2[4]*bmag[6]+0.8944271909999159*bmag[2]*tmp2[5]+0.8944271909999159*tmp2[2]*bmag[5]+bmag[1]*tmp2[3]+tmp2[1]*bmag[3]+bmag[0]*tmp2[2]+tmp2[0]*bmag[2])*volFact)/m_; 
  out[27] += ((0.8*bmag[6]*tmp2[7]+0.8944271909999161*bmag[2]*tmp2[7]+0.8*tmp2[6]*bmag[7]+0.8944271909999161*tmp2[2]*bmag[7]+0.8944271909999161*bmag[1]*tmp2[6]+0.8944271909999161*tmp2[1]*bmag[6]+0.8944271909999159*bmag[3]*tmp2[5]+0.8944271909999159*tmp2[3]*bmag[5]+0.8944271909999159*bmag[3]*tmp2[4]+0.8944271909999159*tmp2[3]*bmag[4]+bmag[0]*tmp2[3]+tmp2[0]*bmag[3]+bmag[1]*tmp2[2]+tmp2[1]*bmag[2])*volFact)/m_; 
  out[28] += ((0.8944271909999159*bmag[7]*tmp2[7]+0.6388765649999399*bmag[6]*tmp2[6]+1.0000000000000002*bmag[2]*tmp2[6]+1.0000000000000002*tmp2[2]*bmag[6]+0.6388765649999399*bmag[4]*tmp2[4]+bmag[0]*tmp2[4]+tmp2[0]*bmag[4]+0.8944271909999159*bmag[3]*tmp2[3]+0.8944271909999159*bmag[1]*tmp2[1])*volFact)/m_; 
  out[29] += ((0.6388765649999399*bmag[7]*tmp2[7]+1.0000000000000002*bmag[1]*tmp2[7]+1.0000000000000002*tmp2[1]*bmag[7]+0.8944271909999159*bmag[6]*tmp2[6]+0.6388765649999399*bmag[5]*tmp2[5]+bmag[0]*tmp2[5]+tmp2[0]*bmag[5]+0.8944271909999159*bmag[3]*tmp2[3]+0.8944271909999159*bmag[2]*tmp2[2])*volFact)/m_; 
  out[30] += ((0.8*bmag[3]*tmp2[7]+0.8*tmp2[3]*bmag[7]+0.8944271909999159*bmag[5]*tmp2[6]+0.6388765649999399*bmag[4]*tmp2[6]+bmag[0]*tmp2[6]+0.8944271909999159*tmp2[5]*bmag[6]+0.6388765649999399*tmp2[4]*bmag[6]+tmp2[0]*bmag[6]+1.0000000000000002*bmag[2]*tmp2[4]+1.0000000000000002*tmp2[2]*bmag[4]+0.8944271909999161*bmag[1]*tmp2[3]+0.8944271909999161*tmp2[1]*bmag[3])*volFact)/m_; 
  out[31] += ((0.6388765649999399*bmag[5]*tmp2[7]+0.8944271909999159*bmag[4]*tmp2[7]+bmag[0]*tmp2[7]+0.6388765649999399*tmp2[5]*bmag[7]+0.8944271909999159*tmp2[4]*bmag[7]+tmp2[0]*bmag[7]+0.8*bmag[3]*tmp2[6]+0.8*tmp2[3]*bmag[6]+1.0000000000000002*bmag[1]*tmp2[5]+1.0000000000000002*tmp2[1]*bmag[5]+0.8944271909999161*bmag[2]*tmp2[3]+0.8944271909999161*tmp2[2]*bmag[3])*volFact)/m_; 
  out[32] += (1.8973665961010278*vmap[0]*vmap1R2*f[13]+1.2727922061357855*vmap1R3*f[3]+2.1213203435596424*vmap0R2*vmap[1]*f[3]+2.1213203435596424*f[0]*vmap[0]*vmap1R2+0.7071067811865475*f[0]*vmap0R3)*volFact; 
  out[33] += (1.897366596101028*vmap[0]*vmap1R2*f[23]+1.2727922061357855*vmap1R3*f[6]+2.1213203435596424*vmap0R2*vmap[1]*f[6]+2.1213203435596424*vmap[0]*f[1]*vmap1R2+0.7071067811865475*vmap0R3*f[1])*volFact; 
  out[34] += (1.897366596101028*vmap[0]*vmap1R2*f[24]+1.2727922061357855*vmap1R3*f[7]+2.1213203435596424*vmap0R2*vmap[1]*f[7]+2.1213203435596424*vmap[0]*vmap1R2*f[2]+0.7071067811865475*vmap0R3*f[2])*volFact; 
  out[35] += (1.8973665961010278*vmap[0]*vmap1R2*f[34]+1.2727922061357855*vmap1R3*f[15]+2.1213203435596424*vmap0R2*vmap[1]*f[15]+2.1213203435596424*vmap[0]*vmap1R2*f[5]+0.7071067811865475*vmap0R3*f[5])*volFact; 
  out[36] += (1.2727922061357853*vmap1R3*f[21]+2.1213203435596424*vmap0R2*vmap[1]*f[21]+2.1213203435596424*vmap[0]*vmap1R2*f[11]+0.7071067811865475*vmap0R3*f[11])*volFact; 
  out[37] += (1.2727922061357853*vmap1R3*f[22]+2.1213203435596424*vmap0R2*vmap[1]*f[22]+2.1213203435596424*vmap[0]*vmap1R2*f[12]+0.7071067811865475*vmap0R3*f[12])*volFact; 
  out[38] += (1.2727922061357853*vmap1R3*f[32]+2.1213203435596424*vmap0R2*vmap[1]*f[32]+2.1213203435596424*vmap[0]*vmap1R2*f[19]+0.7071067811865475*vmap0R3*f[19])*volFact; 
  out[39] += (1.2727922061357853*vmap1R3*f[33]+2.1213203435596424*vmap0R2*vmap[1]*f[33]+2.1213203435596424*vmap[0]*vmap1R2*f[20]+0.7071067811865475*vmap0R3*f[20])*volFact; 
  out[40] += ((bmag[7]*tmp3[7]+bmag[6]*tmp3[6]+bmag[5]*tmp3[5]+bmag[4]*tmp3[4]+bmag[3]*tmp3[3]+bmag[2]*tmp3[2]+bmag[1]*tmp3[1]+bmag[0]*tmp3[0])*volFact)/m_; 
  out[41] += ((1.0000000000000002*bmag[5]*tmp3[7]+1.0000000000000002*tmp3[5]*bmag[7]+0.8944271909999161*bmag[3]*tmp3[6]+0.8944271909999161*tmp3[3]*bmag[6]+0.8944271909999159*bmag[1]*tmp3[4]+0.8944271909999159*tmp3[1]*bmag[4]+bmag[2]*tmp3[3]+tmp3[2]*bmag[3]+bmag[0]*tmp3[1]+tmp3[0]*bmag[1])*volFact)/m_; 
  out[42] += ((0.8944271909999161*bmag[3]*tmp3[7]+0.8944271909999161*tmp3[3]*bmag[7]+1.0000000000000002*bmag[4]*tmp3[6]+1.0000000000000002*tmp3[4]*bmag[6]+0.8944271909999159*bmag[2]*tmp3[5]+0.8944271909999159*tmp3[2]*bmag[5]+bmag[1]*tmp3[3]+tmp3[1]*bmag[3]+bmag[0]*tmp3[2]+tmp3[0]*bmag[2])*volFact)/m_; 
  out[43] += ((0.8*bmag[6]*tmp3[7]+0.8944271909999161*bmag[2]*tmp3[7]+0.8*tmp3[6]*bmag[7]+0.8944271909999161*tmp3[2]*bmag[7]+0.8944271909999161*bmag[1]*tmp3[6]+0.8944271909999161*tmp3[1]*bmag[6]+0.8944271909999159*bmag[3]*tmp3[5]+0.8944271909999159*tmp3[3]*bmag[5]+0.8944271909999159*bmag[3]*tmp3[4]+0.8944271909999159*tmp3[3]*bmag[4]+bmag[0]*tmp3[3]+tmp3[0]*bmag[3]+bmag[1]*tmp3[2]+tmp3[1]*bmag[2])*volFact)/m_; 
  out[44] += ((0.8944271909999159*bmag[7]*tmp3[7]+0.6388765649999399*bmag[6]*tmp3[6]+1.0000000000000002*bmag[2]*tmp3[6]+1.0000000000000002*tmp3[2]*bmag[6]+0.6388765649999399*bmag[4]*tmp3[4]+bmag[0]*tmp3[4]+tmp3[0]*bmag[4]+0.8944271909999159*bmag[3]*tmp3[3]+0.8944271909999159*bmag[1]*tmp3[1])*volFact)/m_; 
  out[45] += ((0.6388765649999399*bmag[7]*tmp3[7]+1.0000000000000002*bmag[1]*tmp3[7]+1.0000000000000002*tmp3[1]*bmag[7]+0.8944271909999159*bmag[6]*tmp3[6]+0.6388765649999399*bmag[5]*tmp3[5]+bmag[0]*tmp3[5]+tmp3[0]*bmag[5]+0.8944271909999159*bmag[3]*tmp3[3]+0.8944271909999159*bmag[2]*tmp3[2])*volFact)/m_; 
  out[46] += ((0.8*bmag[3]*tmp3[7]+0.8*tmp3[3]*bmag[7]+0.8944271909999159*bmag[5]*tmp3[6]+0.6388765649999399*bmag[4]*tmp3[6]+bmag[0]*tmp3[6]+0.8944271909999159*tmp3[5]*bmag[6]+0.6388765649999399*tmp3[4]*bmag[6]+tmp3[0]*bmag[6]+1.0000000000000002*bmag[2]*tmp3[4]+1.0000000000000002*tmp3[2]*bmag[4]+0.8944271909999161*bmag[1]*tmp3[3]+0.8944271909999161*tmp3[1]*bmag[3])*volFact)/m_; 
  out[47] += ((0.6388765649999399*bmag[5]*tmp3[7]+0.8944271909999159*bmag[4]*tmp3[7]+bmag[0]*tmp3[7]+0.6388765649999399*tmp3[5]*bmag[7]+0.8944271909999159*tmp3[4]*bmag[7]+tmp3[0]*bmag[7]+0.8*bmag[3]*tmp3[6]+0.8*tmp3[3]*bmag[6]+1.0000000000000002*bmag[1]*tmp3[5]+1.0000000000000002*tmp3[1]*bmag[5]+0.8944271909999161*bmag[2]*tmp3[3]+0.8944271909999161*tmp3[2]*bmag[3])*volFact)/m_; 
} 

