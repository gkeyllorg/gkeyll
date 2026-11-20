#include <gkyl_mom_gyrokinetic_kernels.h> 
GKYL_CU_DH void int_mom_gyrokinetic_M0_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.25*dxv[0]*dxv[1]; 
 

  out[0] += 2.0*f[0]*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M1_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.25*dxv[0]*dxv[1]; 
 

  out[0] += (1.4142135623730951*vmap[1]*f[2]+1.4142135623730951*f[0]*vmap[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M2PAR_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.25*dxv[0]*dxv[1]; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (0.8944271909999159*vmap1R2*f[5]+2.0*vmap[0]*vmap[1]*f[2]+f[0]*vmap1R2+f[0]*vmap0R2)*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M2_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.25*dxv[0]*dxv[1]; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (0.8944271909999159*vmap1R2*f[5]+2.0*vmap[0]*vmap[1]*f[2]+f[0]*vmap1R2+f[0]*vmap0R2)*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M3PAR_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.25*dxv[0]*dxv[1]; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap0R3 = pow(vmap[0],3);
  const double vmap1R2 = pow(vmap[1],2);
  const double vmap1R3 = pow(vmap[1],3);

  out[0] += (1.8973665961010278*vmap[0]*vmap1R2*f[5]+1.2727922061357855*vmap1R3*f[2]+2.1213203435596424*vmap0R2*vmap[1]*f[2]+2.1213203435596424*f[0]*vmap[0]*vmap1R2+0.7071067811865475*f[0]*vmap0R3)*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M3_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.25*dxv[0]*dxv[1]; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap0R3 = pow(vmap[0],3);
  const double vmap1R2 = pow(vmap[1],2);
  const double vmap1R3 = pow(vmap[1],3);

  out[0] += (1.8973665961010278*vmap[0]*vmap1R2*f[5]+1.2727922061357855*vmap1R3*f[2]+2.1213203435596424*vmap0R2*vmap[1]*f[2]+2.1213203435596424*f[0]*vmap[0]*vmap1R2+0.7071067811865475*f[0]*vmap0R3)*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.25*dxv[0]*dxv[1]; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += 2.0*f[0]*volFact; 
  out[1] += (1.4142135623730951*vmap[1]*f[2]+1.4142135623730951*f[0]*vmap[0])*volFact; 
  out[2] += (0.8944271909999159*vmap1R2*f[5]+2.0*vmap[0]*vmap[1]*f[2]+f[0]*vmap1R2+f[0]*vmap0R2)*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2M3_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.25*dxv[0]*dxv[1]; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap0R3 = pow(vmap[0],3);
  const double vmap1R2 = pow(vmap[1],2);
  const double vmap1R3 = pow(vmap[1],3);

  out[0] += 2.0*f[0]*volFact; 
  out[1] += (1.4142135623730951*vmap[1]*f[2]+1.4142135623730951*f[0]*vmap[0])*volFact; 
  out[2] += (0.8944271909999159*vmap1R2*f[5]+2.0*vmap[0]*vmap[1]*f[2]+f[0]*vmap1R2+f[0]*vmap0R2)*volFact; 
  out[3] += (1.8973665961010278*vmap[0]*vmap1R2*f[5]+1.2727922061357855*vmap1R3*f[2]+2.1213203435596424*vmap0R2*vmap[1]*f[2]+2.1213203435596424*f[0]*vmap[0]*vmap1R2+0.7071067811865475*f[0]*vmap0R3)*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_HAMILTONIAN_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.25*dxv[0]*dxv[1]; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += 2.0*f[0]*volFact; 
  out[1] += (1.4142135623730951*vmap[1]*f[2]+1.4142135623730951*f[0]*vmap[0])*m_*volFact; 
  out[2] += (1.4142135623730951*phi[2]*f[4]*q_+1.4142135623730951*f[1]*phi[1]*q_+1.4142135623730951*f[0]*phi[0]*q_)*volFact+(0.4472135954999579*vmap1R2*f[5]+vmap[0]*vmap[1]*f[2]+0.5*f[0]*vmap1R2+0.5*f[0]*vmap0R2)*m_*volFact; 
} 

