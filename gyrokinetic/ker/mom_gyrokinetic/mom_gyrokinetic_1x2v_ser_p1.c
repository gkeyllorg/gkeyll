#include <gkyl_mom_gyrokinetic_kernels.h> 
GKYL_CU_DH void mom_gyrokinetic_M0_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[1]*dxv[2]/m_; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
} 

GKYL_CU_DH void mom_gyrokinetic_M1_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[1]*dxv[2]/m_; 
  out[0] += (1.4142135623730951*vmap[1]*f[2]+1.4142135623730951*f[0]*vmap[0])*volFact; 
  out[1] += (1.4142135623730951*vmap[1]*f[4]+1.4142135623730951*vmap[0]*f[1])*volFact; 
} 

GKYL_CU_DH void mom_gyrokinetic_M2_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[1]*dxv[2]/m_; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (0.8944271909999159*vmap1R2*f[8]+2.0*vmap[0]*vmap[1]*f[2]+f[0]*vmap1R2+f[0]*vmap0R2)*volFact; 
  out[1] += (0.8944271909999161*vmap1R2*f[9]+2.0*vmap[0]*vmap[1]*f[4]+f[1]*vmap1R2+vmap0R2*f[1])*volFact; 
  double tmp[2]; 
  tmp[0] = 1.4142135623730951*f[3]*vmap[3]+1.4142135623730951*f[0]*vmap[2]; 
  tmp[1] = 1.4142135623730951*vmap[3]*f[5]+1.4142135623730951*f[1]*vmap[2]; 
  out[0] += (2.0*(0.7071067811865475*bmag[1]*tmp[1]+0.7071067811865475*bmag[0]*tmp[0])*volFact)/m_; 
  out[1] += (2.0*(0.7071067811865475*bmag[0]*tmp[1]+0.7071067811865475*tmp[0]*bmag[1])*volFact)/m_; 
} 

GKYL_CU_DH void mom_gyrokinetic_M3_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[1]*dxv[2]/m_; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap0R3 = pow(vmap[0],3);
  const double vmap1R2 = pow(vmap[1],2);
  const double vmap1R3 = pow(vmap[1],3);

  out[0] += ((1.4142135623730951*bmag[1]*vmap[1]*vmap[3]*f[7]+1.4142135623730951*bmag[0]*vmap[1]*vmap[3]*f[6]+1.4142135623730951*vmap[0]*bmag[1]*vmap[3]*f[5]+1.4142135623730951*bmag[1]*vmap[1]*vmap[2]*f[4]+1.4142135623730951*bmag[0]*vmap[0]*f[3]*vmap[3]+1.4142135623730951*bmag[0]*vmap[1]*f[2]*vmap[2]+1.4142135623730951*vmap[0]*bmag[1]*f[1]*vmap[2]+1.4142135623730951*bmag[0]*f[0]*vmap[0]*vmap[2])*volFact)/m_+(1.8973665961010278*vmap[0]*vmap1R2*f[8]+1.2727922061357855*vmap1R3*f[2]+2.1213203435596424*vmap0R2*vmap[1]*f[2]+2.1213203435596424*f[0]*vmap[0]*vmap1R2+0.7071067811865475*f[0]*vmap0R3)*volFact; 
  out[1] += ((1.4142135623730951*bmag[0]*vmap[1]*vmap[3]*f[7]+1.4142135623730951*bmag[1]*vmap[1]*vmap[3]*f[6]+1.4142135623730951*bmag[0]*vmap[0]*vmap[3]*f[5]+1.4142135623730951*bmag[0]*vmap[1]*vmap[2]*f[4]+1.4142135623730951*vmap[0]*bmag[1]*f[3]*vmap[3]+1.4142135623730951*bmag[1]*vmap[1]*f[2]*vmap[2]+1.4142135623730951*bmag[0]*vmap[0]*f[1]*vmap[2]+1.4142135623730951*f[0]*vmap[0]*bmag[1]*vmap[2])*volFact)/m_+(1.897366596101028*vmap[0]*vmap1R2*f[9]+1.2727922061357855*vmap1R3*f[4]+2.1213203435596424*vmap0R2*vmap[1]*f[4]+2.1213203435596424*vmap[0]*f[1]*vmap1R2+0.7071067811865475*vmap0R3*f[1])*volFact; 
} 

GKYL_CU_DH void mom_gyrokinetic_M2PAR_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[1]*dxv[2]/m_; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (0.8944271909999159*vmap1R2*f[8]+2.0*vmap[0]*vmap[1]*f[2]+f[0]*vmap1R2+f[0]*vmap0R2)*volFact; 
  out[1] += (0.8944271909999161*vmap1R2*f[9]+2.0*vmap[0]*vmap[1]*f[4]+f[1]*vmap1R2+vmap0R2*f[1])*volFact; 
} 

GKYL_CU_DH void mom_gyrokinetic_M3PAR_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[1]*dxv[2]/m_; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap0R3 = pow(vmap[0],3);
  const double vmap1R2 = pow(vmap[1],2);
  const double vmap1R3 = pow(vmap[1],3);

  out[0] += (1.8973665961010278*vmap[0]*vmap1R2*f[8]+1.2727922061357855*vmap1R3*f[2]+2.1213203435596424*vmap0R2*vmap[1]*f[2]+2.1213203435596424*f[0]*vmap[0]*vmap1R2+0.7071067811865475*f[0]*vmap0R3)*volFact; 
  out[1] += (1.897366596101028*vmap[0]*vmap1R2*f[9]+1.2727922061357855*vmap1R3*f[4]+2.1213203435596424*vmap0R2*vmap[1]*f[4]+2.1213203435596424*vmap[0]*f[1]*vmap1R2+0.7071067811865475*vmap0R3*f[1])*volFact; 
} 

GKYL_CU_DH void mom_gyrokinetic_M0M1M2_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[1]*dxv[2]/m_; 
  double tmp[2]; 
  tmp[0] = 1.4142135623730951*f[3]*vmap[3]+1.4142135623730951*f[0]*vmap[2]; 
  tmp[1] = 1.4142135623730951*vmap[3]*f[5]+1.4142135623730951*f[1]*vmap[2]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += (1.4142135623730951*vmap[1]*f[2]+1.4142135623730951*f[0]*vmap[0])*volFact; 
  out[3] += (1.4142135623730951*vmap[1]*f[4]+1.4142135623730951*vmap[0]*f[1])*volFact; 
  out[4] += ((1.4142135623730951*bmag[1]*tmp[1]+1.4142135623730951*bmag[0]*tmp[0])*volFact)/m_+(0.8944271909999159*vmap1R2*f[8]+2.0*vmap[0]*vmap[1]*f[2]+f[0]*vmap1R2+f[0]*vmap0R2)*volFact; 
  out[5] += ((1.4142135623730951*bmag[0]*tmp[1]+1.4142135623730951*tmp[0]*bmag[1])*volFact)/m_+(0.8944271909999161*vmap1R2*f[9]+2.0*vmap[0]*vmap[1]*f[4]+f[1]*vmap1R2+vmap0R2*f[1])*volFact; 
} 

GKYL_CU_DH void mom_gyrokinetic_M0M1M2M3_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[1]*dxv[2]/m_; 
  double tmp2[2]; 
  tmp2[0] = 1.4142135623730951*f[3]*vmap[3]+1.4142135623730951*f[0]*vmap[2]; 
  tmp2[1] = 1.4142135623730951*vmap[3]*f[5]+1.4142135623730951*f[1]*vmap[2]; 
  double tmp3[2]; 
  tmp3[0] = vmap[1]*vmap[3]*f[6]+vmap[0]*f[3]*vmap[3]+vmap[1]*f[2]*vmap[2]+f[0]*vmap[0]*vmap[2]; 
  tmp3[1] = vmap[1]*vmap[3]*f[7]+vmap[0]*vmap[3]*f[5]+vmap[1]*vmap[2]*f[4]+vmap[0]*f[1]*vmap[2]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap0R3 = pow(vmap[0],3);
  const double vmap1R2 = pow(vmap[1],2);
  const double vmap1R3 = pow(vmap[1],3);

  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += (1.4142135623730951*vmap[1]*f[2]+1.4142135623730951*f[0]*vmap[0])*volFact; 
  out[3] += (1.4142135623730951*vmap[1]*f[4]+1.4142135623730951*vmap[0]*f[1])*volFact; 
  out[4] += ((1.4142135623730951*bmag[1]*tmp2[1]+1.4142135623730951*bmag[0]*tmp2[0])*volFact)/m_+(0.8944271909999159*vmap1R2*f[8]+2.0*vmap[0]*vmap[1]*f[2]+f[0]*vmap1R2+f[0]*vmap0R2)*volFact; 
  out[5] += ((1.4142135623730951*bmag[0]*tmp2[1]+1.4142135623730951*tmp2[0]*bmag[1])*volFact)/m_+(0.8944271909999161*vmap1R2*f[9]+2.0*vmap[0]*vmap[1]*f[4]+f[1]*vmap1R2+vmap0R2*f[1])*volFact; 
  out[6] += ((1.4142135623730951*bmag[1]*tmp3[1]+1.4142135623730951*bmag[0]*tmp3[0])*volFact)/m_+(1.8973665961010278*vmap[0]*vmap1R2*f[8]+1.2727922061357855*vmap1R3*f[2]+2.1213203435596424*vmap0R2*vmap[1]*f[2]+2.1213203435596424*f[0]*vmap[0]*vmap1R2+0.7071067811865475*f[0]*vmap0R3)*volFact; 
  out[7] += ((1.4142135623730951*bmag[0]*tmp3[1]+1.4142135623730951*tmp3[0]*bmag[1])*volFact)/m_+(1.897366596101028*vmap[0]*vmap1R2*f[9]+1.2727922061357855*vmap1R3*f[4]+2.1213203435596424*vmap0R2*vmap[1]*f[4]+2.1213203435596424*vmap[0]*f[1]*vmap1R2+0.7071067811865475*vmap0R3*f[1])*volFact; 
} 

GKYL_CU_DH void mom_gyrokinetic_HAMILTONIAN_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[1]*dxv[2]/m_; 
  double tmp[2]; 
  tmp[0] = 1.4142135623730951*f[3]*vmap[3]+1.4142135623730951*f[0]*vmap[2]; 
  tmp[1] = 1.4142135623730951*vmap[3]*f[5]+1.4142135623730951*f[1]*vmap[2]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += (1.4142135623730951*vmap[1]*f[2]+1.4142135623730951*f[0]*vmap[0])*m_*volFact; 
  out[3] += (1.4142135623730951*vmap[1]*f[4]+1.4142135623730951*vmap[0]*f[1])*m_*volFact; 
  out[4] += (1.4142135623730951*f[1]*phi[1]+1.4142135623730951*f[0]*phi[0])*q_*volFact+(0.4472135954999579*vmap1R2*f[8]+vmap[0]*vmap[1]*f[2]+0.5*f[0]*vmap1R2+0.5*f[0]*vmap0R2)*m_*volFact+(0.7071067811865475*bmag[1]*tmp[1]+0.7071067811865475*bmag[0]*tmp[0])*volFact; 
  out[5] += (1.4142135623730951*f[0]*phi[1]+1.4142135623730951*phi[0]*f[1])*q_*volFact+(0.44721359549995804*vmap1R2*f[9]+vmap[0]*vmap[1]*f[4]+0.5*f[1]*vmap1R2+0.5*vmap0R2*f[1])*m_*volFact+(0.7071067811865475*bmag[0]*tmp[1]+0.7071067811865475*tmp[0]*bmag[1])*volFact; 
} 

GKYL_CU_DH void mom_gyrokinetic_M2PERP_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[1]*dxv[2]/m_; 
  double tmp[2]; 
  tmp[0] = 1.4142135623730951*f[3]*vmap[3]+1.4142135623730951*f[0]*vmap[2]; 
  tmp[1] = 1.4142135623730951*vmap[3]*f[5]+1.4142135623730951*f[1]*vmap[2]; 
  out[0] += (2.0*(0.7071067811865475*bmag[1]*tmp[1]+0.7071067811865475*bmag[0]*tmp[0])*volFact)/m_; 
  out[1] += (2.0*(0.7071067811865475*bmag[0]*tmp[1]+0.7071067811865475*tmp[0]*bmag[1])*volFact)/m_; 
} 

GKYL_CU_DH void mom_gyrokinetic_M3PERP_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[1]*dxv[2]/m_; 

  out[0] += ((1.4142135623730951*bmag[1]*vmap[1]*vmap[3]*f[7]+1.4142135623730951*bmag[0]*vmap[1]*vmap[3]*f[6]+1.4142135623730951*vmap[0]*bmag[1]*vmap[3]*f[5]+1.4142135623730951*bmag[1]*vmap[1]*vmap[2]*f[4]+1.4142135623730951*bmag[0]*vmap[0]*f[3]*vmap[3]+1.4142135623730951*bmag[0]*vmap[1]*f[2]*vmap[2]+1.4142135623730951*vmap[0]*bmag[1]*f[1]*vmap[2]+1.4142135623730951*bmag[0]*f[0]*vmap[0]*vmap[2])*volFact)/m_; 
  out[1] += ((1.4142135623730951*bmag[0]*vmap[1]*vmap[3]*f[7]+1.4142135623730951*bmag[1]*vmap[1]*vmap[3]*f[6]+1.4142135623730951*bmag[0]*vmap[0]*vmap[3]*f[5]+1.4142135623730951*bmag[0]*vmap[1]*vmap[2]*f[4]+1.4142135623730951*vmap[0]*bmag[1]*f[3]*vmap[3]+1.4142135623730951*bmag[1]*vmap[1]*f[2]*vmap[2]+1.4142135623730951*bmag[0]*vmap[0]*f[1]*vmap[2]+1.4142135623730951*f[0]*vmap[0]*bmag[1]*vmap[2])*volFact)/m_; 
} 

GKYL_CU_DH void mom_gyrokinetic_M0M1M2PARM2PERP_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[1]*dxv[2]/m_; 
  double tmp[2]; 
  tmp[0] = 1.4142135623730951*f[3]*vmap[3]+1.4142135623730951*f[0]*vmap[2]; 
  tmp[1] = 1.4142135623730951*vmap[3]*f[5]+1.4142135623730951*f[1]*vmap[2]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += (1.4142135623730951*vmap[1]*f[2]+1.4142135623730951*f[0]*vmap[0])*volFact; 
  out[3] += (1.4142135623730951*vmap[1]*f[4]+1.4142135623730951*vmap[0]*f[1])*volFact; 
  out[4] += (0.8944271909999159*vmap1R2*f[8]+2.0*vmap[0]*vmap[1]*f[2]+f[0]*vmap1R2+f[0]*vmap0R2)*volFact; 
  out[5] += (0.8944271909999161*vmap1R2*f[9]+2.0*vmap[0]*vmap[1]*f[4]+f[1]*vmap1R2+vmap0R2*f[1])*volFact; 
  out[6] += ((1.4142135623730951*bmag[1]*tmp[1]+1.4142135623730951*bmag[0]*tmp[0])*volFact)/m_; 
  out[7] += ((1.4142135623730951*bmag[0]*tmp[1]+1.4142135623730951*tmp[0]*bmag[1])*volFact)/m_; 
} 

GKYL_CU_DH void mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[1]*dxv[2]/m_; 
  double tmp2[2]; 
  tmp2[0] = 1.4142135623730951*f[3]*vmap[3]+1.4142135623730951*f[0]*vmap[2]; 
  tmp2[1] = 1.4142135623730951*vmap[3]*f[5]+1.4142135623730951*f[1]*vmap[2]; 
  double tmp3[2]; 
  tmp3[0] = vmap[1]*vmap[3]*f[6]+vmap[0]*f[3]*vmap[3]+vmap[1]*f[2]*vmap[2]+f[0]*vmap[0]*vmap[2]; 
  tmp3[1] = vmap[1]*vmap[3]*f[7]+vmap[0]*vmap[3]*f[5]+vmap[1]*vmap[2]*f[4]+vmap[0]*f[1]*vmap[2]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap0R3 = pow(vmap[0],3);
  const double vmap1R2 = pow(vmap[1],2);
  const double vmap1R3 = pow(vmap[1],3);

  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += (1.4142135623730951*vmap[1]*f[2]+1.4142135623730951*f[0]*vmap[0])*volFact; 
  out[3] += (1.4142135623730951*vmap[1]*f[4]+1.4142135623730951*vmap[0]*f[1])*volFact; 
  out[4] += (0.8944271909999159*vmap1R2*f[8]+2.0*vmap[0]*vmap[1]*f[2]+f[0]*vmap1R2+f[0]*vmap0R2)*volFact; 
  out[5] += (0.8944271909999161*vmap1R2*f[9]+2.0*vmap[0]*vmap[1]*f[4]+f[1]*vmap1R2+vmap0R2*f[1])*volFact; 
  out[6] += ((1.4142135623730951*bmag[1]*tmp2[1]+1.4142135623730951*bmag[0]*tmp2[0])*volFact)/m_; 
  out[7] += ((1.4142135623730951*bmag[0]*tmp2[1]+1.4142135623730951*tmp2[0]*bmag[1])*volFact)/m_; 
  out[8] += (1.8973665961010278*vmap[0]*vmap1R2*f[8]+1.2727922061357855*vmap1R3*f[2]+2.1213203435596424*vmap0R2*vmap[1]*f[2]+2.1213203435596424*f[0]*vmap[0]*vmap1R2+0.7071067811865475*f[0]*vmap0R3)*volFact; 
  out[9] += (1.897366596101028*vmap[0]*vmap1R2*f[9]+1.2727922061357855*vmap1R3*f[4]+2.1213203435596424*vmap0R2*vmap[1]*f[4]+2.1213203435596424*vmap[0]*f[1]*vmap1R2+0.7071067811865475*vmap0R3*f[1])*volFact; 
  out[10] += ((1.4142135623730951*bmag[1]*tmp3[1]+1.4142135623730951*bmag[0]*tmp3[0])*volFact)/m_; 
  out[11] += ((1.4142135623730951*bmag[0]*tmp3[1]+1.4142135623730951*tmp3[0]*bmag[1])*volFact)/m_; 
} 

