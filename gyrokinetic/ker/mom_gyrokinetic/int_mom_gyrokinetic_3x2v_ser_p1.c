#include <gkyl_mom_gyrokinetic_kernels.h> 
GKYL_CU_DH void int_mom_gyrokinetic_M0_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.19634954084936207*dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]/m_; 
 

  out[0] += 5.656854249492382*f[0]*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M1_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.19634954084936207*dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]/m_; 
 

  out[0] += (4.0*vmap[1]*f[4]+4.0*f[0]*vmap[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M2PAR_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.19634954084936207*dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (2.5298221281347044*vmap1R2*f[32]+5.656854249492382*vmap[0]*vmap[1]*f[4]+2.8284271247461907*f[0]*vmap1R2+2.8284271247461907*f[0]*vmap0R2)*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M2_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.19634954084936207*dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]/m_; 
 
  double tmp[8]; 
  tmp[0] = (2.8284271247461907*vmap[3]*f[5])/m_+(2.8284271247461907*f[0]*vmap[2])/m_; 
  tmp[1] = (2.8284271247461907*vmap[3]*f[12])/m_+(2.8284271247461907*f[1]*vmap[2])/m_; 
  tmp[2] = (2.8284271247461907*vmap[3]*f[13])/m_+(2.8284271247461907*f[2]*vmap[2])/m_; 
  tmp[3] = (2.8284271247461907*vmap[3]*f[14])/m_+(2.8284271247461907*vmap[2]*f[3])/m_; 
  tmp[4] = (2.8284271247461907*vmap[3]*f[20])/m_+(2.8284271247461907*vmap[2]*f[6])/m_; 
  tmp[5] = (2.8284271247461907*vmap[3]*f[21])/m_+(2.8284271247461907*vmap[2]*f[7])/m_; 
  tmp[6] = (2.8284271247461907*vmap[3]*f[22])/m_+(2.8284271247461907*vmap[2]*f[8])/m_; 
  tmp[7] = (2.8284271247461907*vmap[3]*f[27])/m_+(2.8284271247461907*vmap[2]*f[16])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (2.5298221281347044*vmap1R2*f[32]+bmag[7]*tmp[7]+bmag[6]*tmp[6]+bmag[5]*tmp[5]+bmag[4]*tmp[4]+5.656854249492382*vmap[0]*vmap[1]*f[4]+bmag[3]*tmp[3]+bmag[2]*tmp[2]+2.8284271247461907*f[0]*vmap1R2+bmag[1]*tmp[1]+2.8284271247461907*f[0]*vmap0R2+bmag[0]*tmp[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M3PAR_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.19634954084936207*dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap0R3 = pow(vmap[0],3);
  const double vmap1R2 = pow(vmap[1],2);
  const double vmap1R3 = pow(vmap[1],3);

  out[0] += (5.366563145999495*vmap[0]*vmap1R2*f[32]+3.6*vmap1R3*f[4]+6.0*vmap0R2*vmap[1]*f[4]+6.0*f[0]*vmap[0]*vmap1R2+2.0*f[0]*vmap0R3)*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M3_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.19634954084936207*dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]/m_; 
 
  double tmp[8]; 
  tmp[0] = (2.0*vmap[1]*vmap[3]*f[15])/m_+(2.0*vmap[0]*vmap[3]*f[5])/m_+(2.0*vmap[1]*vmap[2]*f[4])/m_+(2.0*f[0]*vmap[0]*vmap[2])/m_; 
  tmp[1] = (2.0*vmap[1]*vmap[3]*f[23])/m_+(2.0*vmap[0]*vmap[3]*f[12])/m_+(2.0*vmap[1]*vmap[2]*f[9])/m_+(2.0*vmap[0]*f[1]*vmap[2])/m_; 
  tmp[2] = (2.0*vmap[1]*vmap[3]*f[24])/m_+(2.0*vmap[0]*vmap[3]*f[13])/m_+(2.0*vmap[1]*vmap[2]*f[10])/m_+(2.0*vmap[0]*f[2]*vmap[2])/m_; 
  tmp[3] = (2.0*vmap[1]*vmap[3]*f[25])/m_+(2.0*vmap[0]*vmap[3]*f[14])/m_+(2.0*vmap[1]*vmap[2]*f[11])/m_+(2.0*vmap[0]*vmap[2]*f[3])/m_; 
  tmp[4] = (2.0*vmap[1]*vmap[3]*f[28])/m_+(2.0*vmap[0]*vmap[3]*f[20])/m_+(2.0*vmap[1]*vmap[2]*f[17])/m_+(2.0*vmap[0]*vmap[2]*f[6])/m_; 
  tmp[5] = (2.0*vmap[1]*vmap[3]*f[29])/m_+(2.0*vmap[0]*vmap[3]*f[21])/m_+(2.0*vmap[1]*vmap[2]*f[18])/m_+(2.0*vmap[0]*vmap[2]*f[7])/m_; 
  tmp[6] = (2.0*vmap[1]*vmap[3]*f[30])/m_+(2.0*vmap[0]*vmap[3]*f[22])/m_+(2.0*vmap[1]*vmap[2]*f[19])/m_+(2.0*vmap[0]*vmap[2]*f[8])/m_; 
  tmp[7] = (2.0*vmap[1]*vmap[3]*f[31])/m_+(2.0*vmap[0]*vmap[3]*f[27])/m_+(2.0*vmap[1]*vmap[2]*f[26])/m_+(2.0*vmap[0]*vmap[2]*f[16])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap0R3 = pow(vmap[0],3);
  const double vmap1R2 = pow(vmap[1],2);
  const double vmap1R3 = pow(vmap[1],3);

  out[0] += (5.366563145999495*vmap[0]*vmap1R2*f[32]+bmag[7]*tmp[7]+bmag[6]*tmp[6]+bmag[5]*tmp[5]+bmag[4]*tmp[4]+3.6*vmap1R3*f[4]+6.0*vmap0R2*vmap[1]*f[4]+bmag[3]*tmp[3]+bmag[2]*tmp[2]+6.0*f[0]*vmap[0]*vmap1R2+bmag[1]*tmp[1]+2.0*f[0]*vmap0R3+bmag[0]*tmp[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.19634954084936207*dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]/m_; 
 
  double tmp[8]; 
  tmp[0] = (2.8284271247461907*vmap[3]*f[5])/m_+(2.8284271247461907*f[0]*vmap[2])/m_; 
  tmp[1] = (2.8284271247461907*vmap[3]*f[12])/m_+(2.8284271247461907*f[1]*vmap[2])/m_; 
  tmp[2] = (2.8284271247461907*vmap[3]*f[13])/m_+(2.8284271247461907*f[2]*vmap[2])/m_; 
  tmp[3] = (2.8284271247461907*vmap[3]*f[14])/m_+(2.8284271247461907*vmap[2]*f[3])/m_; 
  tmp[4] = (2.8284271247461907*vmap[3]*f[20])/m_+(2.8284271247461907*vmap[2]*f[6])/m_; 
  tmp[5] = (2.8284271247461907*vmap[3]*f[21])/m_+(2.8284271247461907*vmap[2]*f[7])/m_; 
  tmp[6] = (2.8284271247461907*vmap[3]*f[22])/m_+(2.8284271247461907*vmap[2]*f[8])/m_; 
  tmp[7] = (2.8284271247461907*vmap[3]*f[27])/m_+(2.8284271247461907*vmap[2]*f[16])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += 5.656854249492382*f[0]*volFact; 
  out[1] += (4.0*vmap[1]*f[4]+4.0*f[0]*vmap[0])*volFact; 
  out[2] += (2.5298221281347044*vmap1R2*f[32]+bmag[7]*tmp[7]+bmag[6]*tmp[6]+bmag[5]*tmp[5]+bmag[4]*tmp[4]+5.656854249492382*vmap[0]*vmap[1]*f[4]+bmag[3]*tmp[3]+bmag[2]*tmp[2]+2.8284271247461907*f[0]*vmap1R2+bmag[1]*tmp[1]+2.8284271247461907*f[0]*vmap0R2+bmag[0]*tmp[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2M3_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.19634954084936207*dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]/m_; 
 
  double tmp2[8]; 
  tmp2[0] = (2.8284271247461907*vmap[3]*f[5])/m_+(2.8284271247461907*f[0]*vmap[2])/m_; 
  tmp2[1] = (2.8284271247461907*vmap[3]*f[12])/m_+(2.8284271247461907*f[1]*vmap[2])/m_; 
  tmp2[2] = (2.8284271247461907*vmap[3]*f[13])/m_+(2.8284271247461907*f[2]*vmap[2])/m_; 
  tmp2[3] = (2.8284271247461907*vmap[3]*f[14])/m_+(2.8284271247461907*vmap[2]*f[3])/m_; 
  tmp2[4] = (2.8284271247461907*vmap[3]*f[20])/m_+(2.8284271247461907*vmap[2]*f[6])/m_; 
  tmp2[5] = (2.8284271247461907*vmap[3]*f[21])/m_+(2.8284271247461907*vmap[2]*f[7])/m_; 
  tmp2[6] = (2.8284271247461907*vmap[3]*f[22])/m_+(2.8284271247461907*vmap[2]*f[8])/m_; 
  tmp2[7] = (2.8284271247461907*vmap[3]*f[27])/m_+(2.8284271247461907*vmap[2]*f[16])/m_; 
 
  double tmp3[8]; 
  tmp3[0] = (2.0*vmap[1]*vmap[3]*f[15])/m_+(2.0*vmap[0]*vmap[3]*f[5])/m_+(2.0*vmap[1]*vmap[2]*f[4])/m_+(2.0*f[0]*vmap[0]*vmap[2])/m_; 
  tmp3[1] = (2.0*vmap[1]*vmap[3]*f[23])/m_+(2.0*vmap[0]*vmap[3]*f[12])/m_+(2.0*vmap[1]*vmap[2]*f[9])/m_+(2.0*vmap[0]*f[1]*vmap[2])/m_; 
  tmp3[2] = (2.0*vmap[1]*vmap[3]*f[24])/m_+(2.0*vmap[0]*vmap[3]*f[13])/m_+(2.0*vmap[1]*vmap[2]*f[10])/m_+(2.0*vmap[0]*f[2]*vmap[2])/m_; 
  tmp3[3] = (2.0*vmap[1]*vmap[3]*f[25])/m_+(2.0*vmap[0]*vmap[3]*f[14])/m_+(2.0*vmap[1]*vmap[2]*f[11])/m_+(2.0*vmap[0]*vmap[2]*f[3])/m_; 
  tmp3[4] = (2.0*vmap[1]*vmap[3]*f[28])/m_+(2.0*vmap[0]*vmap[3]*f[20])/m_+(2.0*vmap[1]*vmap[2]*f[17])/m_+(2.0*vmap[0]*vmap[2]*f[6])/m_; 
  tmp3[5] = (2.0*vmap[1]*vmap[3]*f[29])/m_+(2.0*vmap[0]*vmap[3]*f[21])/m_+(2.0*vmap[1]*vmap[2]*f[18])/m_+(2.0*vmap[0]*vmap[2]*f[7])/m_; 
  tmp3[6] = (2.0*vmap[1]*vmap[3]*f[30])/m_+(2.0*vmap[0]*vmap[3]*f[22])/m_+(2.0*vmap[1]*vmap[2]*f[19])/m_+(2.0*vmap[0]*vmap[2]*f[8])/m_; 
  tmp3[7] = (2.0*vmap[1]*vmap[3]*f[31])/m_+(2.0*vmap[0]*vmap[3]*f[27])/m_+(2.0*vmap[1]*vmap[2]*f[26])/m_+(2.0*vmap[0]*vmap[2]*f[16])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap0R3 = pow(vmap[0],3);
  const double vmap1R2 = pow(vmap[1],2);
  const double vmap1R3 = pow(vmap[1],3);

  out[0] += 5.656854249492382*f[0]*volFact; 
  out[1] += (4.0*vmap[1]*f[4]+4.0*f[0]*vmap[0])*volFact; 
  out[2] += (2.5298221281347044*vmap1R2*f[32]+bmag[7]*tmp2[7]+bmag[6]*tmp2[6]+bmag[5]*tmp2[5]+bmag[4]*tmp2[4]+5.656854249492382*vmap[0]*vmap[1]*f[4]+bmag[3]*tmp2[3]+bmag[2]*tmp2[2]+2.8284271247461907*f[0]*vmap1R2+bmag[1]*tmp2[1]+2.8284271247461907*f[0]*vmap0R2+bmag[0]*tmp2[0])*volFact; 
  out[3] += (5.366563145999495*vmap[0]*vmap1R2*f[32]+bmag[7]*tmp3[7]+bmag[6]*tmp3[6]+bmag[5]*tmp3[5]+bmag[4]*tmp3[4]+3.6*vmap1R3*f[4]+6.0*vmap0R2*vmap[1]*f[4]+bmag[3]*tmp3[3]+bmag[2]*tmp3[2]+6.0*f[0]*vmap[0]*vmap1R2+bmag[1]*tmp3[1]+2.0*f[0]*vmap0R3+bmag[0]*tmp3[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_HAMILTONIAN_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.19634954084936207*dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]/m_; 
 
  double tmp[8]; 
  tmp[0] = 1.4142135623730951*vmap[3]*f[5]+1.4142135623730951*f[0]*vmap[2]; 
  tmp[1] = 1.4142135623730951*vmap[3]*f[12]+1.4142135623730951*f[1]*vmap[2]; 
  tmp[2] = 1.4142135623730951*vmap[3]*f[13]+1.4142135623730951*f[2]*vmap[2]; 
  tmp[3] = 1.4142135623730951*vmap[3]*f[14]+1.4142135623730951*vmap[2]*f[3]; 
  tmp[4] = 1.4142135623730951*vmap[3]*f[20]+1.4142135623730951*vmap[2]*f[6]; 
  tmp[5] = 1.4142135623730951*vmap[3]*f[21]+1.4142135623730951*vmap[2]*f[7]; 
  tmp[6] = 1.4142135623730951*vmap[3]*f[22]+1.4142135623730951*vmap[2]*f[8]; 
  tmp[7] = 1.4142135623730951*vmap[3]*f[27]+1.4142135623730951*vmap[2]*f[16]; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += 5.656854249492382*f[0]*volFact; 
  out[1] += (4.0*vmap[1]*f[4]+4.0*f[0]*vmap[0])*m_*volFact; 
  out[2] += (2.0*phi[7]*f[16]*q_+2.0*phi[6]*f[8]*q_+2.0*phi[5]*f[7]*q_+2.0*phi[4]*f[6]*q_+2.0*f[3]*phi[3]*q_+2.0*f[2]*phi[2]*q_+2.0*f[1]*phi[1]*q_+2.0*f[0]*phi[0]*q_+bmag[7]*tmp[7]+bmag[6]*tmp[6]+bmag[5]*tmp[5]+bmag[4]*tmp[4]+bmag[3]*tmp[3]+bmag[2]*tmp[2]+bmag[1]*tmp[1]+bmag[0]*tmp[0])*volFact+(1.264911064067352*vmap1R2*f[32]+2.8284271247461907*vmap[0]*vmap[1]*f[4]+1.4142135623730951*f[0]*vmap1R2+1.4142135623730951*f[0]*vmap0R2)*m_*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M2PERP_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.19634954084936207*dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]/m_; 
 
  double tmp[8]; 
  tmp[0] = (2.8284271247461907*vmap[3]*f[5])/m_+(2.8284271247461907*f[0]*vmap[2])/m_; 
  tmp[1] = (2.8284271247461907*vmap[3]*f[12])/m_+(2.8284271247461907*f[1]*vmap[2])/m_; 
  tmp[2] = (2.8284271247461907*vmap[3]*f[13])/m_+(2.8284271247461907*f[2]*vmap[2])/m_; 
  tmp[3] = (2.8284271247461907*vmap[3]*f[14])/m_+(2.8284271247461907*vmap[2]*f[3])/m_; 
  tmp[4] = (2.8284271247461907*vmap[3]*f[20])/m_+(2.8284271247461907*vmap[2]*f[6])/m_; 
  tmp[5] = (2.8284271247461907*vmap[3]*f[21])/m_+(2.8284271247461907*vmap[2]*f[7])/m_; 
  tmp[6] = (2.8284271247461907*vmap[3]*f[22])/m_+(2.8284271247461907*vmap[2]*f[8])/m_; 
  tmp[7] = (2.8284271247461907*vmap[3]*f[27])/m_+(2.8284271247461907*vmap[2]*f[16])/m_; 
 

  out[0] += (bmag[7]*tmp[7]+bmag[6]*tmp[6]+bmag[5]*tmp[5]+bmag[4]*tmp[4]+bmag[3]*tmp[3]+bmag[2]*tmp[2]+bmag[1]*tmp[1]+bmag[0]*tmp[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M3PERP_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.19634954084936207*dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]/m_; 
 
  double tmp[8]; 
  tmp[0] = (2.0*vmap[1]*vmap[3]*f[15])/m_+(2.0*vmap[0]*vmap[3]*f[5])/m_+(2.0*vmap[1]*vmap[2]*f[4])/m_+(2.0*f[0]*vmap[0]*vmap[2])/m_; 
  tmp[1] = (2.0*vmap[1]*vmap[3]*f[23])/m_+(2.0*vmap[0]*vmap[3]*f[12])/m_+(2.0*vmap[1]*vmap[2]*f[9])/m_+(2.0*vmap[0]*f[1]*vmap[2])/m_; 
  tmp[2] = (2.0*vmap[1]*vmap[3]*f[24])/m_+(2.0*vmap[0]*vmap[3]*f[13])/m_+(2.0*vmap[1]*vmap[2]*f[10])/m_+(2.0*vmap[0]*f[2]*vmap[2])/m_; 
  tmp[3] = (2.0*vmap[1]*vmap[3]*f[25])/m_+(2.0*vmap[0]*vmap[3]*f[14])/m_+(2.0*vmap[1]*vmap[2]*f[11])/m_+(2.0*vmap[0]*vmap[2]*f[3])/m_; 
  tmp[4] = (2.0*vmap[1]*vmap[3]*f[28])/m_+(2.0*vmap[0]*vmap[3]*f[20])/m_+(2.0*vmap[1]*vmap[2]*f[17])/m_+(2.0*vmap[0]*vmap[2]*f[6])/m_; 
  tmp[5] = (2.0*vmap[1]*vmap[3]*f[29])/m_+(2.0*vmap[0]*vmap[3]*f[21])/m_+(2.0*vmap[1]*vmap[2]*f[18])/m_+(2.0*vmap[0]*vmap[2]*f[7])/m_; 
  tmp[6] = (2.0*vmap[1]*vmap[3]*f[30])/m_+(2.0*vmap[0]*vmap[3]*f[22])/m_+(2.0*vmap[1]*vmap[2]*f[19])/m_+(2.0*vmap[0]*vmap[2]*f[8])/m_; 
  tmp[7] = (2.0*vmap[1]*vmap[3]*f[31])/m_+(2.0*vmap[0]*vmap[3]*f[27])/m_+(2.0*vmap[1]*vmap[2]*f[26])/m_+(2.0*vmap[0]*vmap[2]*f[16])/m_; 
 

  out[0] += (bmag[7]*tmp[7]+bmag[6]*tmp[6]+bmag[5]*tmp[5]+bmag[4]*tmp[4]+bmag[3]*tmp[3]+bmag[2]*tmp[2]+bmag[1]*tmp[1]+bmag[0]*tmp[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2PARM2PERP_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.19634954084936207*dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]/m_; 
 
  double tmp[8]; 
  tmp[0] = (2.8284271247461907*vmap[3]*f[5])/m_+(2.8284271247461907*f[0]*vmap[2])/m_; 
  tmp[1] = (2.8284271247461907*vmap[3]*f[12])/m_+(2.8284271247461907*f[1]*vmap[2])/m_; 
  tmp[2] = (2.8284271247461907*vmap[3]*f[13])/m_+(2.8284271247461907*f[2]*vmap[2])/m_; 
  tmp[3] = (2.8284271247461907*vmap[3]*f[14])/m_+(2.8284271247461907*vmap[2]*f[3])/m_; 
  tmp[4] = (2.8284271247461907*vmap[3]*f[20])/m_+(2.8284271247461907*vmap[2]*f[6])/m_; 
  tmp[5] = (2.8284271247461907*vmap[3]*f[21])/m_+(2.8284271247461907*vmap[2]*f[7])/m_; 
  tmp[6] = (2.8284271247461907*vmap[3]*f[22])/m_+(2.8284271247461907*vmap[2]*f[8])/m_; 
  tmp[7] = (2.8284271247461907*vmap[3]*f[27])/m_+(2.8284271247461907*vmap[2]*f[16])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += 5.656854249492382*f[0]*volFact; 
  out[1] += (4.0*vmap[1]*f[4]+4.0*f[0]*vmap[0])*volFact; 
  out[2] += (2.5298221281347044*vmap1R2*f[32]+5.656854249492382*vmap[0]*vmap[1]*f[4]+2.8284271247461907*f[0]*vmap1R2+2.8284271247461907*f[0]*vmap0R2)*volFact; 
  out[3] += (bmag[7]*tmp[7]+bmag[6]*tmp[6]+bmag[5]*tmp[5]+bmag[4]*tmp[4]+bmag[3]*tmp[3]+bmag[2]*tmp[2]+bmag[1]*tmp[1]+bmag[0]*tmp[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.19634954084936207*dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]/m_; 
 
  double tmp2[8]; 
  tmp2[0] = (2.8284271247461907*vmap[3]*f[5])/m_+(2.8284271247461907*f[0]*vmap[2])/m_; 
  tmp2[1] = (2.8284271247461907*vmap[3]*f[12])/m_+(2.8284271247461907*f[1]*vmap[2])/m_; 
  tmp2[2] = (2.8284271247461907*vmap[3]*f[13])/m_+(2.8284271247461907*f[2]*vmap[2])/m_; 
  tmp2[3] = (2.8284271247461907*vmap[3]*f[14])/m_+(2.8284271247461907*vmap[2]*f[3])/m_; 
  tmp2[4] = (2.8284271247461907*vmap[3]*f[20])/m_+(2.8284271247461907*vmap[2]*f[6])/m_; 
  tmp2[5] = (2.8284271247461907*vmap[3]*f[21])/m_+(2.8284271247461907*vmap[2]*f[7])/m_; 
  tmp2[6] = (2.8284271247461907*vmap[3]*f[22])/m_+(2.8284271247461907*vmap[2]*f[8])/m_; 
  tmp2[7] = (2.8284271247461907*vmap[3]*f[27])/m_+(2.8284271247461907*vmap[2]*f[16])/m_; 
 
  double tmp3[8]; 
  tmp3[0] = (2.0*vmap[1]*vmap[3]*f[15])/m_+(2.0*vmap[0]*vmap[3]*f[5])/m_+(2.0*vmap[1]*vmap[2]*f[4])/m_+(2.0*f[0]*vmap[0]*vmap[2])/m_; 
  tmp3[1] = (2.0*vmap[1]*vmap[3]*f[23])/m_+(2.0*vmap[0]*vmap[3]*f[12])/m_+(2.0*vmap[1]*vmap[2]*f[9])/m_+(2.0*vmap[0]*f[1]*vmap[2])/m_; 
  tmp3[2] = (2.0*vmap[1]*vmap[3]*f[24])/m_+(2.0*vmap[0]*vmap[3]*f[13])/m_+(2.0*vmap[1]*vmap[2]*f[10])/m_+(2.0*vmap[0]*f[2]*vmap[2])/m_; 
  tmp3[3] = (2.0*vmap[1]*vmap[3]*f[25])/m_+(2.0*vmap[0]*vmap[3]*f[14])/m_+(2.0*vmap[1]*vmap[2]*f[11])/m_+(2.0*vmap[0]*vmap[2]*f[3])/m_; 
  tmp3[4] = (2.0*vmap[1]*vmap[3]*f[28])/m_+(2.0*vmap[0]*vmap[3]*f[20])/m_+(2.0*vmap[1]*vmap[2]*f[17])/m_+(2.0*vmap[0]*vmap[2]*f[6])/m_; 
  tmp3[5] = (2.0*vmap[1]*vmap[3]*f[29])/m_+(2.0*vmap[0]*vmap[3]*f[21])/m_+(2.0*vmap[1]*vmap[2]*f[18])/m_+(2.0*vmap[0]*vmap[2]*f[7])/m_; 
  tmp3[6] = (2.0*vmap[1]*vmap[3]*f[30])/m_+(2.0*vmap[0]*vmap[3]*f[22])/m_+(2.0*vmap[1]*vmap[2]*f[19])/m_+(2.0*vmap[0]*vmap[2]*f[8])/m_; 
  tmp3[7] = (2.0*vmap[1]*vmap[3]*f[31])/m_+(2.0*vmap[0]*vmap[3]*f[27])/m_+(2.0*vmap[1]*vmap[2]*f[26])/m_+(2.0*vmap[0]*vmap[2]*f[16])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap0R3 = pow(vmap[0],3);
  const double vmap1R2 = pow(vmap[1],2);
  const double vmap1R3 = pow(vmap[1],3);

  out[0] += 5.656854249492382*f[0]*volFact; 
  out[1] += (4.0*vmap[1]*f[4]+4.0*f[0]*vmap[0])*volFact; 
  out[2] += (2.5298221281347044*vmap1R2*f[32]+5.656854249492382*vmap[0]*vmap[1]*f[4]+2.8284271247461907*f[0]*vmap1R2+2.8284271247461907*f[0]*vmap0R2)*volFact; 
  out[3] += (bmag[7]*tmp2[7]+bmag[6]*tmp2[6]+bmag[5]*tmp2[5]+bmag[4]*tmp2[4]+bmag[3]*tmp2[3]+bmag[2]*tmp2[2]+bmag[1]*tmp2[1]+bmag[0]*tmp2[0])*volFact; 
  out[4] += (5.366563145999495*vmap[0]*vmap1R2*f[32]+3.6*vmap1R3*f[4]+6.0*vmap0R2*vmap[1]*f[4]+6.0*f[0]*vmap[0]*vmap1R2+2.0*f[0]*vmap0R3)*volFact; 
  out[5] += (bmag[7]*tmp3[7]+bmag[6]*tmp3[6]+bmag[5]*tmp3[5]+bmag[4]*tmp3[4]+bmag[3]*tmp3[3]+bmag[2]*tmp3[2]+bmag[1]*tmp3[1]+bmag[0]*tmp3[0])*volFact; 
} 

