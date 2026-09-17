#include <gkyl_fpo_vlasov_kernels.h> 

GKYL_CU_DH double fpo_vlasov_drag_vol_1x3v_ser_p1(const double* dxv, const double* drag_coeff, const double* f, double* GKYL_RESTRICT out) { 
  // dxv[NDIM]: Cell spacing in each direction. 
  // drag_coeff: Drag coefficient. 
  // f: Distribution function. 
  // out: Incremented output. 

  // dv1 for each direction. 
  double dv1_vx = 2.0/dxv[1]; 
  double dv1_vy = 2.0/dxv[2]; 
  double dv1_vz = 2.0/dxv[3]; 
  
  // Index into drag coefficient. 
  const double* a_vx = &drag_coeff[0]; 
  const double* a_vy = &drag_coeff[16]; 
  const double* a_vz = &drag_coeff[32]; 
  
  // vx contribution 
  out[2] += 0.21650635094610965*(a_vx[15]*f[15]+a_vx[14]*f[14]+a_vx[13]*f[13]+a_vx[12]*f[12]+a_vx[11]*f[11]+a_vx[10]*f[10]+a_vx[9]*f[9]+a_vx[8]*f[8]+a_vx[7]*f[7]+a_vx[6]*f[6]+a_vx[5]*f[5]+a_vx[4]*f[4]+a_vx[3]*f[3]+a_vx[2]*f[2]+a_vx[1]*f[1]+a_vx[0]*f[0])*dv1_vx; 
  out[5] += 0.21650635094610965*(a_vx[14]*f[15]+f[14]*a_vx[15]+a_vx[10]*f[13]+f[10]*a_vx[13]+a_vx[9]*f[12]+f[9]*a_vx[12]+a_vx[7]*f[11]+f[7]*a_vx[11]+a_vx[4]*f[8]+f[4]*a_vx[8]+a_vx[3]*f[6]+f[3]*a_vx[6]+a_vx[2]*f[5]+f[2]*a_vx[5]+a_vx[0]*f[1]+f[0]*a_vx[1])*dv1_vx; 
  out[7] += 0.21650635094610965*(a_vx[12]*f[15]+f[12]*a_vx[15]+a_vx[9]*f[14]+f[9]*a_vx[14]+a_vx[8]*f[13]+f[8]*a_vx[13]+a_vx[5]*f[11]+f[5]*a_vx[11]+a_vx[4]*f[10]+f[4]*a_vx[10]+a_vx[2]*f[7]+f[2]*a_vx[7]+a_vx[1]*f[6]+f[1]*a_vx[6]+a_vx[0]*f[3]+f[0]*a_vx[3])*dv1_vx; 
  out[9] += 0.21650635094610965*(a_vx[11]*f[15]+f[11]*a_vx[15]+a_vx[7]*f[14]+f[7]*a_vx[14]+a_vx[6]*f[13]+f[6]*a_vx[13]+a_vx[5]*f[12]+f[5]*a_vx[12]+a_vx[3]*f[10]+f[3]*a_vx[10]+a_vx[2]*f[9]+f[2]*a_vx[9]+a_vx[1]*f[8]+f[1]*a_vx[8]+a_vx[0]*f[4]+f[0]*a_vx[4])*dv1_vx; 
  out[11] += 0.21650635094610965*(a_vx[9]*f[15]+f[9]*a_vx[15]+a_vx[12]*f[14]+f[12]*a_vx[14]+a_vx[4]*f[13]+f[4]*a_vx[13]+a_vx[2]*f[11]+f[2]*a_vx[11]+a_vx[8]*f[10]+f[8]*a_vx[10]+a_vx[5]*f[7]+f[5]*a_vx[7]+a_vx[0]*f[6]+f[0]*a_vx[6]+a_vx[1]*f[3]+f[1]*a_vx[3])*dv1_vx; 
  out[12] += 0.21650635094610965*(a_vx[7]*f[15]+f[7]*a_vx[15]+a_vx[11]*f[14]+f[11]*a_vx[14]+a_vx[3]*f[13]+f[3]*a_vx[13]+a_vx[2]*f[12]+f[2]*a_vx[12]+a_vx[6]*f[10]+f[6]*a_vx[10]+a_vx[5]*f[9]+f[5]*a_vx[9]+a_vx[0]*f[8]+f[0]*a_vx[8]+a_vx[1]*f[4]+f[1]*a_vx[4])*dv1_vx; 
  out[14] += 0.21650635094610965*(a_vx[5]*f[15]+f[5]*a_vx[15]+a_vx[2]*f[14]+f[2]*a_vx[14]+a_vx[1]*f[13]+f[1]*a_vx[13]+a_vx[11]*f[12]+f[11]*a_vx[12]+a_vx[0]*f[10]+f[0]*a_vx[10]+a_vx[7]*f[9]+f[7]*a_vx[9]+a_vx[6]*f[8]+f[6]*a_vx[8]+a_vx[3]*f[4]+f[3]*a_vx[4])*dv1_vx; 
  out[15] += 0.21650635094610965*(a_vx[2]*f[15]+f[2]*a_vx[15]+a_vx[5]*f[14]+f[5]*a_vx[14]+a_vx[0]*f[13]+f[0]*a_vx[13]+a_vx[7]*f[12]+f[7]*a_vx[12]+a_vx[9]*f[11]+f[9]*a_vx[11]+a_vx[1]*f[10]+f[1]*a_vx[10]+a_vx[3]*f[8]+f[3]*a_vx[8]+a_vx[4]*f[6]+f[4]*a_vx[6])*dv1_vx; 

  // vy contribution 
  out[3] += 0.21650635094610965*(a_vy[15]*f[15]+a_vy[14]*f[14]+a_vy[13]*f[13]+a_vy[12]*f[12]+a_vy[11]*f[11]+a_vy[10]*f[10]+a_vy[9]*f[9]+a_vy[8]*f[8]+a_vy[7]*f[7]+a_vy[6]*f[6]+a_vy[5]*f[5]+a_vy[4]*f[4]+a_vy[3]*f[3]+a_vy[2]*f[2]+a_vy[1]*f[1]+a_vy[0]*f[0])*dv1_vy; 
  out[6] += 0.21650635094610965*(a_vy[14]*f[15]+f[14]*a_vy[15]+a_vy[10]*f[13]+f[10]*a_vy[13]+a_vy[9]*f[12]+f[9]*a_vy[12]+a_vy[7]*f[11]+f[7]*a_vy[11]+a_vy[4]*f[8]+f[4]*a_vy[8]+a_vy[3]*f[6]+f[3]*a_vy[6]+a_vy[2]*f[5]+f[2]*a_vy[5]+a_vy[0]*f[1]+f[0]*a_vy[1])*dv1_vy; 
  out[7] += 0.21650635094610965*(a_vy[13]*f[15]+f[13]*a_vy[15]+a_vy[10]*f[14]+f[10]*a_vy[14]+a_vy[8]*f[12]+f[8]*a_vy[12]+a_vy[6]*f[11]+f[6]*a_vy[11]+a_vy[4]*f[9]+f[4]*a_vy[9]+a_vy[3]*f[7]+f[3]*a_vy[7]+a_vy[1]*f[5]+f[1]*a_vy[5]+a_vy[0]*f[2]+f[0]*a_vy[2])*dv1_vy; 
  out[10] += 0.21650635094610965*(a_vy[11]*f[15]+f[11]*a_vy[15]+a_vy[7]*f[14]+f[7]*a_vy[14]+a_vy[6]*f[13]+f[6]*a_vy[13]+a_vy[5]*f[12]+f[5]*a_vy[12]+a_vy[3]*f[10]+f[3]*a_vy[10]+a_vy[2]*f[9]+f[2]*a_vy[9]+a_vy[1]*f[8]+f[1]*a_vy[8]+a_vy[0]*f[4]+f[0]*a_vy[4])*dv1_vy; 
  out[11] += 0.21650635094610965*(a_vy[10]*f[15]+f[10]*a_vy[15]+a_vy[13]*f[14]+f[13]*a_vy[14]+a_vy[4]*f[12]+f[4]*a_vy[12]+a_vy[3]*f[11]+f[3]*a_vy[11]+a_vy[8]*f[9]+f[8]*a_vy[9]+a_vy[6]*f[7]+f[6]*a_vy[7]+a_vy[0]*f[5]+f[0]*a_vy[5]+a_vy[1]*f[2]+f[1]*a_vy[2])*dv1_vy; 
  out[13] += 0.21650635094610965*(a_vy[7]*f[15]+f[7]*a_vy[15]+a_vy[11]*f[14]+f[11]*a_vy[14]+a_vy[3]*f[13]+f[3]*a_vy[13]+a_vy[2]*f[12]+f[2]*a_vy[12]+a_vy[6]*f[10]+f[6]*a_vy[10]+a_vy[5]*f[9]+f[5]*a_vy[9]+a_vy[0]*f[8]+f[0]*a_vy[8]+a_vy[1]*f[4]+f[1]*a_vy[4])*dv1_vy; 
  out[14] += 0.21650635094610965*(a_vy[6]*f[15]+f[6]*a_vy[15]+a_vy[3]*f[14]+f[3]*a_vy[14]+a_vy[11]*f[13]+f[11]*a_vy[13]+a_vy[1]*f[12]+f[1]*a_vy[12]+a_vy[7]*f[10]+f[7]*a_vy[10]+a_vy[0]*f[9]+f[0]*a_vy[9]+a_vy[5]*f[8]+f[5]*a_vy[8]+a_vy[2]*f[4]+f[2]*a_vy[4])*dv1_vy; 
  out[15] += 0.21650635094610965*(a_vy[3]*f[15]+f[3]*a_vy[15]+a_vy[6]*f[14]+f[6]*a_vy[14]+a_vy[7]*f[13]+f[7]*a_vy[13]+a_vy[0]*f[12]+f[0]*a_vy[12]+a_vy[10]*f[11]+f[10]*a_vy[11]+a_vy[1]*f[9]+f[1]*a_vy[9]+a_vy[2]*f[8]+f[2]*a_vy[8]+a_vy[4]*f[5]+f[4]*a_vy[5])*dv1_vy; 

  // vz contribution 
  out[4] += 0.21650635094610965*(a_vz[15]*f[15]+a_vz[14]*f[14]+a_vz[13]*f[13]+a_vz[12]*f[12]+a_vz[11]*f[11]+a_vz[10]*f[10]+a_vz[9]*f[9]+a_vz[8]*f[8]+a_vz[7]*f[7]+a_vz[6]*f[6]+a_vz[5]*f[5]+a_vz[4]*f[4]+a_vz[3]*f[3]+a_vz[2]*f[2]+a_vz[1]*f[1]+a_vz[0]*f[0])*dv1_vz; 
  out[8] += 0.21650635094610965*(a_vz[14]*f[15]+f[14]*a_vz[15]+a_vz[10]*f[13]+f[10]*a_vz[13]+a_vz[9]*f[12]+f[9]*a_vz[12]+a_vz[7]*f[11]+f[7]*a_vz[11]+a_vz[4]*f[8]+f[4]*a_vz[8]+a_vz[3]*f[6]+f[3]*a_vz[6]+a_vz[2]*f[5]+f[2]*a_vz[5]+a_vz[0]*f[1]+f[0]*a_vz[1])*dv1_vz; 
  out[9] += 0.21650635094610965*(a_vz[13]*f[15]+f[13]*a_vz[15]+a_vz[10]*f[14]+f[10]*a_vz[14]+a_vz[8]*f[12]+f[8]*a_vz[12]+a_vz[6]*f[11]+f[6]*a_vz[11]+a_vz[4]*f[9]+f[4]*a_vz[9]+a_vz[3]*f[7]+f[3]*a_vz[7]+a_vz[1]*f[5]+f[1]*a_vz[5]+a_vz[0]*f[2]+f[0]*a_vz[2])*dv1_vz; 
  out[10] += 0.21650635094610965*(a_vz[12]*f[15]+f[12]*a_vz[15]+a_vz[9]*f[14]+f[9]*a_vz[14]+a_vz[8]*f[13]+f[8]*a_vz[13]+a_vz[5]*f[11]+f[5]*a_vz[11]+a_vz[4]*f[10]+f[4]*a_vz[10]+a_vz[2]*f[7]+f[2]*a_vz[7]+a_vz[1]*f[6]+f[1]*a_vz[6]+a_vz[0]*f[3]+f[0]*a_vz[3])*dv1_vz; 
  out[12] += 0.21650635094610965*(a_vz[10]*f[15]+f[10]*a_vz[15]+a_vz[13]*f[14]+f[13]*a_vz[14]+a_vz[4]*f[12]+f[4]*a_vz[12]+a_vz[3]*f[11]+f[3]*a_vz[11]+a_vz[8]*f[9]+f[8]*a_vz[9]+a_vz[6]*f[7]+f[6]*a_vz[7]+a_vz[0]*f[5]+f[0]*a_vz[5]+a_vz[1]*f[2]+f[1]*a_vz[2])*dv1_vz; 
  out[13] += 0.21650635094610965*(a_vz[9]*f[15]+f[9]*a_vz[15]+a_vz[12]*f[14]+f[12]*a_vz[14]+a_vz[4]*f[13]+f[4]*a_vz[13]+a_vz[2]*f[11]+f[2]*a_vz[11]+a_vz[8]*f[10]+f[8]*a_vz[10]+a_vz[5]*f[7]+f[5]*a_vz[7]+a_vz[0]*f[6]+f[0]*a_vz[6]+a_vz[1]*f[3]+f[1]*a_vz[3])*dv1_vz; 
  out[14] += 0.21650635094610965*(a_vz[8]*f[15]+f[8]*a_vz[15]+a_vz[4]*f[14]+f[4]*a_vz[14]+a_vz[12]*f[13]+f[12]*a_vz[13]+a_vz[1]*f[11]+f[1]*a_vz[11]+a_vz[9]*f[10]+f[9]*a_vz[10]+a_vz[0]*f[7]+f[0]*a_vz[7]+a_vz[5]*f[6]+f[5]*a_vz[6]+a_vz[2]*f[3]+f[2]*a_vz[3])*dv1_vz; 
  out[15] += 0.21650635094610965*(a_vz[4]*f[15]+f[4]*a_vz[15]+a_vz[8]*f[14]+f[8]*a_vz[14]+a_vz[9]*f[13]+f[9]*a_vz[13]+a_vz[10]*f[12]+f[10]*a_vz[12]+a_vz[0]*f[11]+f[0]*a_vz[11]+a_vz[1]*f[7]+f[1]*a_vz[7]+a_vz[2]*f[6]+f[2]*a_vz[6]+a_vz[3]*f[5]+f[3]*a_vz[5])*dv1_vz; 



  return 0.0; 

}
