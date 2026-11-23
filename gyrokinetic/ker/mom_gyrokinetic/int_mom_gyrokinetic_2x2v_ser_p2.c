#include <gkyl_mom_gyrokinetic_kernels.h> 
GKYL_CU_DH void int_mom_gyrokinetic_M0_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.39269908169872414*dxv[0]*dxv[1]*dxv[2]*dxv[3]/m_; 
 

  out[0] += 4.0*f[0]*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M1_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.39269908169872414*dxv[0]*dxv[1]*dxv[2]*dxv[3]/m_; 
 

  out[0] += (2.8284271247461907*vmap[1]*f[3]+2.8284271247461907*f[0]*vmap[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M2PAR_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.39269908169872414*dxv[0]*dxv[1]*dxv[2]*dxv[3]/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (1.7888543819998317*vmap1R2*f[13]+4.0*vmap[0]*vmap[1]*f[3]+2.0*f[0]*vmap1R2+2.0*f[0]*vmap0R2)*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M2_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.39269908169872414*dxv[0]*dxv[1]*dxv[2]*dxv[3]/m_; 
 
  double tmp[8]; 
  tmp[0] = (2.8284271247461907*vmap[3]*f[4])/m_+(2.8284271247461907*f[0]*vmap[2])/m_; 
  tmp[1] = (2.8284271247461907*vmap[3]*f[8])/m_+(2.8284271247461907*f[1]*vmap[2])/m_; 
  tmp[2] = (2.8284271247461907*vmap[3]*f[9])/m_+(2.8284271247461907*f[2]*vmap[2])/m_; 
  tmp[3] = (2.8284271247461907*vmap[3]*f[16])/m_+(2.8284271247461907*vmap[2]*f[5])/m_; 
  tmp[4] = (2.828427124746191*vmap[3]*f[25])/m_+(2.8284271247461907*vmap[2]*f[11])/m_; 
  tmp[5] = (2.828427124746191*vmap[3]*f[26])/m_+(2.8284271247461907*vmap[2]*f[12])/m_; 
  tmp[6] = (2.828427124746191*vmap[3]*f[35])/m_+(2.8284271247461907*vmap[2]*f[19])/m_; 
  tmp[7] = (2.828427124746191*vmap[3]*f[36])/m_+(2.8284271247461907*vmap[2]*f[20])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (1.7888543819998317*vmap1R2*f[13]+bmag[7]*tmp[7]+bmag[6]*tmp[6]+bmag[5]*tmp[5]+bmag[4]*tmp[4]+bmag[3]*tmp[3]+4.0*vmap[0]*vmap[1]*f[3]+bmag[2]*tmp[2]+2.0*f[0]*vmap1R2+bmag[1]*tmp[1]+2.0*f[0]*vmap0R2+bmag[0]*tmp[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M3PAR_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.39269908169872414*dxv[0]*dxv[1]*dxv[2]*dxv[3]/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap0R3 = pow(vmap[0],3);
  const double vmap1R2 = pow(vmap[1],2);
  const double vmap1R3 = pow(vmap[1],3);

  out[0] += (3.7947331922020555*vmap[0]*vmap1R2*f[13]+2.5455844122715714*vmap1R3*f[3]+4.242640687119286*vmap0R2*vmap[1]*f[3]+4.242640687119286*f[0]*vmap[0]*vmap1R2+1.4142135623730951*f[0]*vmap0R3)*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M3_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.39269908169872414*dxv[0]*dxv[1]*dxv[2]*dxv[3]/m_; 
 
  double tmp[8]; 
  tmp[0] = (2.0*vmap[1]*vmap[3]*f[10])/m_+(2.0*vmap[0]*vmap[3]*f[4])/m_+(2.0*vmap[1]*vmap[2]*f[3])/m_+(2.0*f[0]*vmap[0]*vmap[2])/m_; 
  tmp[1] = (2.0*vmap[1]*vmap[3]*f[17])/m_+(2.0*vmap[0]*vmap[3]*f[8])/m_+(2.0*vmap[1]*vmap[2]*f[6])/m_+(2.0*vmap[0]*f[1]*vmap[2])/m_; 
  tmp[2] = (2.0*vmap[1]*vmap[3]*f[18])/m_+(2.0*vmap[0]*vmap[3]*f[9])/m_+(2.0*vmap[1]*vmap[2]*f[7])/m_+(2.0*vmap[0]*f[2]*vmap[2])/m_; 
  tmp[3] = (2.0*vmap[1]*vmap[3]*f[31])/m_+(2.0*vmap[0]*vmap[3]*f[16])/m_+(2.0*vmap[1]*vmap[2]*f[15])/m_+(2.0*vmap[0]*vmap[2]*f[5])/m_; 
  tmp[4] = (2.0*vmap[1]*vmap[3]*f[37])/m_+(2.0000000000000004*vmap[0]*vmap[3]*f[25])/m_+(2.0000000000000004*vmap[1]*vmap[2]*f[21])/m_+(2.0*vmap[0]*vmap[2]*f[11])/m_; 
  tmp[5] = (2.0*vmap[1]*vmap[3]*f[38])/m_+(2.0000000000000004*vmap[0]*vmap[3]*f[26])/m_+(2.0000000000000004*vmap[1]*vmap[2]*f[22])/m_+(2.0*vmap[0]*vmap[2]*f[12])/m_; 
  tmp[6] = (2.0*vmap[1]*vmap[3]*f[44])/m_+(2.0000000000000004*vmap[0]*vmap[3]*f[35])/m_+(2.0000000000000004*vmap[1]*vmap[2]*f[32])/m_+(2.0*vmap[0]*vmap[2]*f[19])/m_; 
  tmp[7] = (2.0*vmap[1]*vmap[3]*f[45])/m_+(2.0000000000000004*vmap[0]*vmap[3]*f[36])/m_+(2.0000000000000004*vmap[1]*vmap[2]*f[33])/m_+(2.0*vmap[0]*vmap[2]*f[20])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap0R3 = pow(vmap[0],3);
  const double vmap1R2 = pow(vmap[1],2);
  const double vmap1R3 = pow(vmap[1],3);

  out[0] += (3.7947331922020555*vmap[0]*vmap1R2*f[13]+bmag[7]*tmp[7]+bmag[6]*tmp[6]+bmag[5]*tmp[5]+bmag[4]*tmp[4]+bmag[3]*tmp[3]+2.5455844122715714*vmap1R3*f[3]+4.242640687119286*vmap0R2*vmap[1]*f[3]+bmag[2]*tmp[2]+4.242640687119286*f[0]*vmap[0]*vmap1R2+bmag[1]*tmp[1]+1.4142135623730951*f[0]*vmap0R3+bmag[0]*tmp[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.39269908169872414*dxv[0]*dxv[1]*dxv[2]*dxv[3]/m_; 
 
  double tmp[8]; 
  tmp[0] = (2.8284271247461907*vmap[3]*f[4])/m_+(2.8284271247461907*f[0]*vmap[2])/m_; 
  tmp[1] = (2.8284271247461907*vmap[3]*f[8])/m_+(2.8284271247461907*f[1]*vmap[2])/m_; 
  tmp[2] = (2.8284271247461907*vmap[3]*f[9])/m_+(2.8284271247461907*f[2]*vmap[2])/m_; 
  tmp[3] = (2.8284271247461907*vmap[3]*f[16])/m_+(2.8284271247461907*vmap[2]*f[5])/m_; 
  tmp[4] = (2.828427124746191*vmap[3]*f[25])/m_+(2.8284271247461907*vmap[2]*f[11])/m_; 
  tmp[5] = (2.828427124746191*vmap[3]*f[26])/m_+(2.8284271247461907*vmap[2]*f[12])/m_; 
  tmp[6] = (2.828427124746191*vmap[3]*f[35])/m_+(2.8284271247461907*vmap[2]*f[19])/m_; 
  tmp[7] = (2.828427124746191*vmap[3]*f[36])/m_+(2.8284271247461907*vmap[2]*f[20])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += 4.0*f[0]*volFact; 
  out[1] += (2.8284271247461907*vmap[1]*f[3]+2.8284271247461907*f[0]*vmap[0])*volFact; 
  out[2] += (1.7888543819998317*vmap1R2*f[13]+bmag[7]*tmp[7]+bmag[6]*tmp[6]+bmag[5]*tmp[5]+bmag[4]*tmp[4]+bmag[3]*tmp[3]+4.0*vmap[0]*vmap[1]*f[3]+bmag[2]*tmp[2]+2.0*f[0]*vmap1R2+bmag[1]*tmp[1]+2.0*f[0]*vmap0R2+bmag[0]*tmp[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2M3_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.39269908169872414*dxv[0]*dxv[1]*dxv[2]*dxv[3]/m_; 
 
  double tmp2[8]; 
  tmp2[0] = (2.8284271247461907*vmap[3]*f[4])/m_+(2.8284271247461907*f[0]*vmap[2])/m_; 
  tmp2[1] = (2.8284271247461907*vmap[3]*f[8])/m_+(2.8284271247461907*f[1]*vmap[2])/m_; 
  tmp2[2] = (2.8284271247461907*vmap[3]*f[9])/m_+(2.8284271247461907*f[2]*vmap[2])/m_; 
  tmp2[3] = (2.8284271247461907*vmap[3]*f[16])/m_+(2.8284271247461907*vmap[2]*f[5])/m_; 
  tmp2[4] = (2.828427124746191*vmap[3]*f[25])/m_+(2.8284271247461907*vmap[2]*f[11])/m_; 
  tmp2[5] = (2.828427124746191*vmap[3]*f[26])/m_+(2.8284271247461907*vmap[2]*f[12])/m_; 
  tmp2[6] = (2.828427124746191*vmap[3]*f[35])/m_+(2.8284271247461907*vmap[2]*f[19])/m_; 
  tmp2[7] = (2.828427124746191*vmap[3]*f[36])/m_+(2.8284271247461907*vmap[2]*f[20])/m_; 
 
  double tmp3[8]; 
  tmp3[0] = (2.0*vmap[1]*vmap[3]*f[10])/m_+(2.0*vmap[0]*vmap[3]*f[4])/m_+(2.0*vmap[1]*vmap[2]*f[3])/m_+(2.0*f[0]*vmap[0]*vmap[2])/m_; 
  tmp3[1] = (2.0*vmap[1]*vmap[3]*f[17])/m_+(2.0*vmap[0]*vmap[3]*f[8])/m_+(2.0*vmap[1]*vmap[2]*f[6])/m_+(2.0*vmap[0]*f[1]*vmap[2])/m_; 
  tmp3[2] = (2.0*vmap[1]*vmap[3]*f[18])/m_+(2.0*vmap[0]*vmap[3]*f[9])/m_+(2.0*vmap[1]*vmap[2]*f[7])/m_+(2.0*vmap[0]*f[2]*vmap[2])/m_; 
  tmp3[3] = (2.0*vmap[1]*vmap[3]*f[31])/m_+(2.0*vmap[0]*vmap[3]*f[16])/m_+(2.0*vmap[1]*vmap[2]*f[15])/m_+(2.0*vmap[0]*vmap[2]*f[5])/m_; 
  tmp3[4] = (2.0*vmap[1]*vmap[3]*f[37])/m_+(2.0000000000000004*vmap[0]*vmap[3]*f[25])/m_+(2.0000000000000004*vmap[1]*vmap[2]*f[21])/m_+(2.0*vmap[0]*vmap[2]*f[11])/m_; 
  tmp3[5] = (2.0*vmap[1]*vmap[3]*f[38])/m_+(2.0000000000000004*vmap[0]*vmap[3]*f[26])/m_+(2.0000000000000004*vmap[1]*vmap[2]*f[22])/m_+(2.0*vmap[0]*vmap[2]*f[12])/m_; 
  tmp3[6] = (2.0*vmap[1]*vmap[3]*f[44])/m_+(2.0000000000000004*vmap[0]*vmap[3]*f[35])/m_+(2.0000000000000004*vmap[1]*vmap[2]*f[32])/m_+(2.0*vmap[0]*vmap[2]*f[19])/m_; 
  tmp3[7] = (2.0*vmap[1]*vmap[3]*f[45])/m_+(2.0000000000000004*vmap[0]*vmap[3]*f[36])/m_+(2.0000000000000004*vmap[1]*vmap[2]*f[33])/m_+(2.0*vmap[0]*vmap[2]*f[20])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap0R3 = pow(vmap[0],3);
  const double vmap1R2 = pow(vmap[1],2);
  const double vmap1R3 = pow(vmap[1],3);

  out[0] += 4.0*f[0]*volFact; 
  out[1] += (2.8284271247461907*vmap[1]*f[3]+2.8284271247461907*f[0]*vmap[0])*volFact; 
  out[2] += (1.7888543819998317*vmap1R2*f[13]+bmag[7]*tmp2[7]+bmag[6]*tmp2[6]+bmag[5]*tmp2[5]+bmag[4]*tmp2[4]+bmag[3]*tmp2[3]+4.0*vmap[0]*vmap[1]*f[3]+bmag[2]*tmp2[2]+2.0*f[0]*vmap1R2+bmag[1]*tmp2[1]+2.0*f[0]*vmap0R2+bmag[0]*tmp2[0])*volFact; 
  out[3] += (3.7947331922020555*vmap[0]*vmap1R2*f[13]+bmag[7]*tmp3[7]+bmag[6]*tmp3[6]+bmag[5]*tmp3[5]+bmag[4]*tmp3[4]+bmag[3]*tmp3[3]+2.5455844122715714*vmap1R3*f[3]+4.242640687119286*vmap0R2*vmap[1]*f[3]+bmag[2]*tmp3[2]+4.242640687119286*f[0]*vmap[0]*vmap1R2+bmag[1]*tmp3[1]+1.4142135623730951*f[0]*vmap0R3+bmag[0]*tmp3[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_HAMILTONIAN_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.39269908169872414*dxv[0]*dxv[1]*dxv[2]*dxv[3]/m_; 
 
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

  out[0] += 4.0*f[0]*volFact; 
  out[1] += (2.8284271247461907*vmap[1]*f[3]+2.8284271247461907*f[0]*vmap[0])*m_*volFact; 
  out[2] += (2.0*phi[7]*f[20]*q_+2.0*phi[6]*f[19]*q_+2.0*phi[5]*f[12]*q_+2.0*phi[4]*f[11]*q_+2.0*phi[3]*f[5]*q_+2.0*f[2]*phi[2]*q_+2.0*f[1]*phi[1]*q_+2.0*f[0]*phi[0]*q_+bmag[7]*tmp[7]+bmag[6]*tmp[6]+bmag[5]*tmp[5]+bmag[4]*tmp[4]+bmag[3]*tmp[3]+bmag[2]*tmp[2]+bmag[1]*tmp[1]+bmag[0]*tmp[0])*volFact+(0.8944271909999159*vmap1R2*f[13]+2.0*vmap[0]*vmap[1]*f[3]+f[0]*vmap1R2+f[0]*vmap0R2)*m_*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M2PERP_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.39269908169872414*dxv[0]*dxv[1]*dxv[2]*dxv[3]/m_; 
 
  double tmp[8]; 
  tmp[0] = (2.8284271247461907*vmap[3]*f[4])/m_+(2.8284271247461907*f[0]*vmap[2])/m_; 
  tmp[1] = (2.8284271247461907*vmap[3]*f[8])/m_+(2.8284271247461907*f[1]*vmap[2])/m_; 
  tmp[2] = (2.8284271247461907*vmap[3]*f[9])/m_+(2.8284271247461907*f[2]*vmap[2])/m_; 
  tmp[3] = (2.8284271247461907*vmap[3]*f[16])/m_+(2.8284271247461907*vmap[2]*f[5])/m_; 
  tmp[4] = (2.828427124746191*vmap[3]*f[25])/m_+(2.8284271247461907*vmap[2]*f[11])/m_; 
  tmp[5] = (2.828427124746191*vmap[3]*f[26])/m_+(2.8284271247461907*vmap[2]*f[12])/m_; 
  tmp[6] = (2.828427124746191*vmap[3]*f[35])/m_+(2.8284271247461907*vmap[2]*f[19])/m_; 
  tmp[7] = (2.828427124746191*vmap[3]*f[36])/m_+(2.8284271247461907*vmap[2]*f[20])/m_; 
 

  out[0] += (bmag[7]*tmp[7]+bmag[6]*tmp[6]+bmag[5]*tmp[5]+bmag[4]*tmp[4]+bmag[3]*tmp[3]+bmag[2]*tmp[2]+bmag[1]*tmp[1]+bmag[0]*tmp[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M3PERP_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.39269908169872414*dxv[0]*dxv[1]*dxv[2]*dxv[3]/m_; 
 
  double tmp[8]; 
  tmp[0] = (2.0*vmap[1]*vmap[3]*f[10])/m_+(2.0*vmap[0]*vmap[3]*f[4])/m_+(2.0*vmap[1]*vmap[2]*f[3])/m_+(2.0*f[0]*vmap[0]*vmap[2])/m_; 
  tmp[1] = (2.0*vmap[1]*vmap[3]*f[17])/m_+(2.0*vmap[0]*vmap[3]*f[8])/m_+(2.0*vmap[1]*vmap[2]*f[6])/m_+(2.0*vmap[0]*f[1]*vmap[2])/m_; 
  tmp[2] = (2.0*vmap[1]*vmap[3]*f[18])/m_+(2.0*vmap[0]*vmap[3]*f[9])/m_+(2.0*vmap[1]*vmap[2]*f[7])/m_+(2.0*vmap[0]*f[2]*vmap[2])/m_; 
  tmp[3] = (2.0*vmap[1]*vmap[3]*f[31])/m_+(2.0*vmap[0]*vmap[3]*f[16])/m_+(2.0*vmap[1]*vmap[2]*f[15])/m_+(2.0*vmap[0]*vmap[2]*f[5])/m_; 
  tmp[4] = (2.0*vmap[1]*vmap[3]*f[37])/m_+(2.0000000000000004*vmap[0]*vmap[3]*f[25])/m_+(2.0000000000000004*vmap[1]*vmap[2]*f[21])/m_+(2.0*vmap[0]*vmap[2]*f[11])/m_; 
  tmp[5] = (2.0*vmap[1]*vmap[3]*f[38])/m_+(2.0000000000000004*vmap[0]*vmap[3]*f[26])/m_+(2.0000000000000004*vmap[1]*vmap[2]*f[22])/m_+(2.0*vmap[0]*vmap[2]*f[12])/m_; 
  tmp[6] = (2.0*vmap[1]*vmap[3]*f[44])/m_+(2.0000000000000004*vmap[0]*vmap[3]*f[35])/m_+(2.0000000000000004*vmap[1]*vmap[2]*f[32])/m_+(2.0*vmap[0]*vmap[2]*f[19])/m_; 
  tmp[7] = (2.0*vmap[1]*vmap[3]*f[45])/m_+(2.0000000000000004*vmap[0]*vmap[3]*f[36])/m_+(2.0000000000000004*vmap[1]*vmap[2]*f[33])/m_+(2.0*vmap[0]*vmap[2]*f[20])/m_; 
 

  out[0] += (bmag[7]*tmp[7]+bmag[6]*tmp[6]+bmag[5]*tmp[5]+bmag[4]*tmp[4]+bmag[3]*tmp[3]+bmag[2]*tmp[2]+bmag[1]*tmp[1]+bmag[0]*tmp[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2PARM2PERP_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.39269908169872414*dxv[0]*dxv[1]*dxv[2]*dxv[3]/m_; 
 
  double tmp[8]; 
  tmp[0] = (2.8284271247461907*vmap[3]*f[4])/m_+(2.8284271247461907*f[0]*vmap[2])/m_; 
  tmp[1] = (2.8284271247461907*vmap[3]*f[8])/m_+(2.8284271247461907*f[1]*vmap[2])/m_; 
  tmp[2] = (2.8284271247461907*vmap[3]*f[9])/m_+(2.8284271247461907*f[2]*vmap[2])/m_; 
  tmp[3] = (2.8284271247461907*vmap[3]*f[16])/m_+(2.8284271247461907*vmap[2]*f[5])/m_; 
  tmp[4] = (2.828427124746191*vmap[3]*f[25])/m_+(2.8284271247461907*vmap[2]*f[11])/m_; 
  tmp[5] = (2.828427124746191*vmap[3]*f[26])/m_+(2.8284271247461907*vmap[2]*f[12])/m_; 
  tmp[6] = (2.828427124746191*vmap[3]*f[35])/m_+(2.8284271247461907*vmap[2]*f[19])/m_; 
  tmp[7] = (2.828427124746191*vmap[3]*f[36])/m_+(2.8284271247461907*vmap[2]*f[20])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += 4.0*f[0]*volFact; 
  out[1] += (2.8284271247461907*vmap[1]*f[3]+2.8284271247461907*f[0]*vmap[0])*volFact; 
  out[2] += (1.7888543819998317*vmap1R2*f[13]+4.0*vmap[0]*vmap[1]*f[3]+2.0*f[0]*vmap1R2+2.0*f[0]*vmap0R2)*volFact; 
  out[3] += (bmag[7]*tmp[7]+bmag[6]*tmp[6]+bmag[5]*tmp[5]+bmag[4]*tmp[4]+bmag[3]*tmp[3]+bmag[2]*tmp[2]+bmag[1]*tmp[1]+bmag[0]*tmp[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.39269908169872414*dxv[0]*dxv[1]*dxv[2]*dxv[3]/m_; 
 
  double tmp2[8]; 
  tmp2[0] = (2.8284271247461907*vmap[3]*f[4])/m_+(2.8284271247461907*f[0]*vmap[2])/m_; 
  tmp2[1] = (2.8284271247461907*vmap[3]*f[8])/m_+(2.8284271247461907*f[1]*vmap[2])/m_; 
  tmp2[2] = (2.8284271247461907*vmap[3]*f[9])/m_+(2.8284271247461907*f[2]*vmap[2])/m_; 
  tmp2[3] = (2.8284271247461907*vmap[3]*f[16])/m_+(2.8284271247461907*vmap[2]*f[5])/m_; 
  tmp2[4] = (2.828427124746191*vmap[3]*f[25])/m_+(2.8284271247461907*vmap[2]*f[11])/m_; 
  tmp2[5] = (2.828427124746191*vmap[3]*f[26])/m_+(2.8284271247461907*vmap[2]*f[12])/m_; 
  tmp2[6] = (2.828427124746191*vmap[3]*f[35])/m_+(2.8284271247461907*vmap[2]*f[19])/m_; 
  tmp2[7] = (2.828427124746191*vmap[3]*f[36])/m_+(2.8284271247461907*vmap[2]*f[20])/m_; 
 
  double tmp3[8]; 
  tmp3[0] = (2.0*vmap[1]*vmap[3]*f[10])/m_+(2.0*vmap[0]*vmap[3]*f[4])/m_+(2.0*vmap[1]*vmap[2]*f[3])/m_+(2.0*f[0]*vmap[0]*vmap[2])/m_; 
  tmp3[1] = (2.0*vmap[1]*vmap[3]*f[17])/m_+(2.0*vmap[0]*vmap[3]*f[8])/m_+(2.0*vmap[1]*vmap[2]*f[6])/m_+(2.0*vmap[0]*f[1]*vmap[2])/m_; 
  tmp3[2] = (2.0*vmap[1]*vmap[3]*f[18])/m_+(2.0*vmap[0]*vmap[3]*f[9])/m_+(2.0*vmap[1]*vmap[2]*f[7])/m_+(2.0*vmap[0]*f[2]*vmap[2])/m_; 
  tmp3[3] = (2.0*vmap[1]*vmap[3]*f[31])/m_+(2.0*vmap[0]*vmap[3]*f[16])/m_+(2.0*vmap[1]*vmap[2]*f[15])/m_+(2.0*vmap[0]*vmap[2]*f[5])/m_; 
  tmp3[4] = (2.0*vmap[1]*vmap[3]*f[37])/m_+(2.0000000000000004*vmap[0]*vmap[3]*f[25])/m_+(2.0000000000000004*vmap[1]*vmap[2]*f[21])/m_+(2.0*vmap[0]*vmap[2]*f[11])/m_; 
  tmp3[5] = (2.0*vmap[1]*vmap[3]*f[38])/m_+(2.0000000000000004*vmap[0]*vmap[3]*f[26])/m_+(2.0000000000000004*vmap[1]*vmap[2]*f[22])/m_+(2.0*vmap[0]*vmap[2]*f[12])/m_; 
  tmp3[6] = (2.0*vmap[1]*vmap[3]*f[44])/m_+(2.0000000000000004*vmap[0]*vmap[3]*f[35])/m_+(2.0000000000000004*vmap[1]*vmap[2]*f[32])/m_+(2.0*vmap[0]*vmap[2]*f[19])/m_; 
  tmp3[7] = (2.0*vmap[1]*vmap[3]*f[45])/m_+(2.0000000000000004*vmap[0]*vmap[3]*f[36])/m_+(2.0000000000000004*vmap[1]*vmap[2]*f[33])/m_+(2.0*vmap[0]*vmap[2]*f[20])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap0R3 = pow(vmap[0],3);
  const double vmap1R2 = pow(vmap[1],2);
  const double vmap1R3 = pow(vmap[1],3);

  out[0] += 4.0*f[0]*volFact; 
  out[1] += (2.8284271247461907*vmap[1]*f[3]+2.8284271247461907*f[0]*vmap[0])*volFact; 
  out[2] += (1.7888543819998317*vmap1R2*f[13]+4.0*vmap[0]*vmap[1]*f[3]+2.0*f[0]*vmap1R2+2.0*f[0]*vmap0R2)*volFact; 
  out[3] += (bmag[7]*tmp2[7]+bmag[6]*tmp2[6]+bmag[5]*tmp2[5]+bmag[4]*tmp2[4]+bmag[3]*tmp2[3]+bmag[2]*tmp2[2]+bmag[1]*tmp2[1]+bmag[0]*tmp2[0])*volFact; 
  out[4] += (3.7947331922020555*vmap[0]*vmap1R2*f[13]+2.5455844122715714*vmap1R3*f[3]+4.242640687119286*vmap0R2*vmap[1]*f[3]+4.242640687119286*f[0]*vmap[0]*vmap1R2+1.4142135623730951*f[0]*vmap0R3)*volFact; 
  out[5] += (bmag[7]*tmp3[7]+bmag[6]*tmp3[6]+bmag[5]*tmp3[5]+bmag[4]*tmp3[4]+bmag[3]*tmp3[3]+bmag[2]*tmp3[2]+bmag[1]*tmp3[1]+bmag[0]*tmp3[0])*volFact; 
} 

