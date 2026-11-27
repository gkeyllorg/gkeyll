#include <gkyl_mom_weighted_gyrokinetic_kernels.h> 
GKYL_CU_DH void int_mom_weighted_gyrokinetic_M0_WEIGHTX_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.25*dxv[0]*dxv[1]; 
 

  out[0] += (1.4142135623730951*wgt[2]*f[4]+1.4142135623730951*f[1]*wgt[1]+1.4142135623730951*f[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_weighted_gyrokinetic_M0_WEIGHTV_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.25*dxv[0]*dxv[1]; 
 

  out[0] += (1.4142135623730951*wgt[2]*f[5]+1.4142135623730951*wgt[1]*f[2]+1.4142135623730951*f[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_weighted_gyrokinetic_M1_WEIGHTX_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.25*dxv[0]*dxv[1]; 
 

  out[0] += (1.0000000000000002*vmap[1]*wgt[2]*f[6]+vmap[0]*wgt[2]*f[4]+vmap[1]*wgt[1]*f[3]+wgt[0]*vmap[1]*f[2]+vmap[0]*f[1]*wgt[1]+f[0]*vmap[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_weighted_gyrokinetic_M1_WEIGHTV_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.25*dxv[0]*dxv[1]; 
 

  out[0] += (vmap[0]*wgt[2]*f[5]+0.8944271909999159*vmap[1]*wgt[1]*f[5]+0.8944271909999159*vmap[1]*f[2]*wgt[2]+vmap[0]*wgt[1]*f[2]+wgt[0]*vmap[1]*f[2]+f[0]*vmap[1]*wgt[1]+f[0]*vmap[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_weighted_gyrokinetic_M2PAR_WEIGHTX_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.25*dxv[0]*dxv[1]; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (0.632455532033676*vmap1R2*wgt[1]*f[7]+1.4142135623730951*vmap[0]*vmap[1]*wgt[2]*f[6]+0.6324555320336759*wgt[0]*vmap1R2*f[5]+0.7071067811865475*vmap1R2*wgt[2]*f[4]+0.7071067811865475*vmap0R2*wgt[2]*f[4]+1.4142135623730951*vmap[0]*vmap[1]*wgt[1]*f[3]+1.4142135623730951*vmap[0]*wgt[0]*vmap[1]*f[2]+0.7071067811865475*f[1]*vmap1R2*wgt[1]+0.7071067811865475*vmap0R2*f[1]*wgt[1]+0.7071067811865475*f[0]*wgt[0]*vmap1R2+0.7071067811865475*f[0]*vmap0R2*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_weighted_gyrokinetic_M2PAR_WEIGHTV_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.25*dxv[0]*dxv[1]; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (1.1111677990074318*vmap1R2*wgt[2]*f[5]+0.7071067811865475*vmap0R2*wgt[2]*f[5]+1.264911064067352*vmap[0]*vmap[1]*wgt[1]*f[5]+0.6324555320336759*wgt[0]*vmap1R2*f[5]+1.264911064067352*vmap[0]*vmap[1]*f[2]*wgt[2]+0.6324555320336759*f[0]*vmap1R2*wgt[2]+1.2727922061357855*vmap1R2*wgt[1]*f[2]+0.7071067811865475*vmap0R2*wgt[1]*f[2]+1.4142135623730951*vmap[0]*wgt[0]*vmap[1]*f[2]+1.4142135623730951*f[0]*vmap[0]*vmap[1]*wgt[1]+0.7071067811865475*f[0]*wgt[0]*vmap1R2+0.7071067811865475*f[0]*vmap0R2*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_weighted_gyrokinetic_M2_WEIGHTX_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.25*dxv[0]*dxv[1]; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (0.632455532033676*vmap1R2*wgt[1]*f[7]+1.4142135623730951*vmap[0]*vmap[1]*wgt[2]*f[6]+0.6324555320336759*wgt[0]*vmap1R2*f[5]+0.7071067811865475*vmap1R2*wgt[2]*f[4]+0.7071067811865475*vmap0R2*wgt[2]*f[4]+1.4142135623730951*vmap[0]*vmap[1]*wgt[1]*f[3]+1.4142135623730951*vmap[0]*wgt[0]*vmap[1]*f[2]+0.7071067811865475*f[1]*vmap1R2*wgt[1]+0.7071067811865475*vmap0R2*f[1]*wgt[1]+0.7071067811865475*f[0]*wgt[0]*vmap1R2+0.7071067811865475*f[0]*vmap0R2*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_weighted_gyrokinetic_M2_WEIGHTV_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.25*dxv[0]*dxv[1]; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (1.1111677990074318*vmap1R2*wgt[2]*f[5]+0.7071067811865475*vmap0R2*wgt[2]*f[5]+1.264911064067352*vmap[0]*vmap[1]*wgt[1]*f[5]+0.6324555320336759*wgt[0]*vmap1R2*f[5]+1.264911064067352*vmap[0]*vmap[1]*f[2]*wgt[2]+0.6324555320336759*f[0]*vmap1R2*wgt[2]+1.2727922061357855*vmap1R2*wgt[1]*f[2]+0.7071067811865475*vmap0R2*wgt[1]*f[2]+1.4142135623730951*vmap[0]*wgt[0]*vmap[1]*f[2]+1.4142135623730951*f[0]*vmap[0]*vmap[1]*wgt[1]+0.7071067811865475*f[0]*wgt[0]*vmap1R2+0.7071067811865475*f[0]*vmap0R2*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_weighted_gyrokinetic_M0M1M2_WEIGHTX_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.25*dxv[0]*dxv[1]; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (1.4142135623730951*wgt[2]*f[4]+1.4142135623730951*f[1]*wgt[1]+1.4142135623730951*f[0]*wgt[0])*volFact; 
  out[1] += (1.0000000000000002*vmap[1]*wgt[2]*f[6]+vmap[0]*wgt[2]*f[4]+vmap[1]*wgt[1]*f[3]+wgt[0]*vmap[1]*f[2]+vmap[0]*f[1]*wgt[1]+f[0]*vmap[0]*wgt[0])*volFact; 
  out[2] += (0.632455532033676*vmap1R2*wgt[1]*f[7]+1.4142135623730951*vmap[0]*vmap[1]*wgt[2]*f[6]+0.6324555320336759*wgt[0]*vmap1R2*f[5]+0.7071067811865475*vmap1R2*wgt[2]*f[4]+0.7071067811865475*vmap0R2*wgt[2]*f[4]+1.4142135623730951*vmap[0]*vmap[1]*wgt[1]*f[3]+1.4142135623730951*vmap[0]*wgt[0]*vmap[1]*f[2]+0.7071067811865475*f[1]*vmap1R2*wgt[1]+0.7071067811865475*vmap0R2*f[1]*wgt[1]+0.7071067811865475*f[0]*wgt[0]*vmap1R2+0.7071067811865475*f[0]*vmap0R2*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_weighted_gyrokinetic_M0M1M2_WEIGHTV_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.25*dxv[0]*dxv[1]; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (1.4142135623730951*wgt[2]*f[5]+1.4142135623730951*wgt[1]*f[2]+1.4142135623730951*f[0]*wgt[0])*volFact; 
  out[1] += (vmap[0]*wgt[2]*f[5]+0.8944271909999159*vmap[1]*wgt[1]*f[5]+0.8944271909999159*vmap[1]*f[2]*wgt[2]+vmap[0]*wgt[1]*f[2]+wgt[0]*vmap[1]*f[2]+f[0]*vmap[1]*wgt[1]+f[0]*vmap[0]*wgt[0])*volFact; 
  out[2] += (1.1111677990074318*vmap1R2*wgt[2]*f[5]+0.7071067811865475*vmap0R2*wgt[2]*f[5]+1.264911064067352*vmap[0]*vmap[1]*wgt[1]*f[5]+0.6324555320336759*wgt[0]*vmap1R2*f[5]+1.264911064067352*vmap[0]*vmap[1]*f[2]*wgt[2]+0.6324555320336759*f[0]*vmap1R2*wgt[2]+1.2727922061357855*vmap1R2*wgt[1]*f[2]+0.7071067811865475*vmap0R2*wgt[1]*f[2]+1.4142135623730951*vmap[0]*wgt[0]*vmap[1]*f[2]+1.4142135623730951*f[0]*vmap[0]*vmap[1]*wgt[1]+0.7071067811865475*f[0]*wgt[0]*vmap1R2+0.7071067811865475*f[0]*vmap0R2*wgt[0])*volFact; 
} 

