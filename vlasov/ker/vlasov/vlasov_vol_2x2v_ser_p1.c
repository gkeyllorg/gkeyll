#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // field:      q/m*EM fields.
  // cot_vec:   Only used in gen geo.
  // f:         Input distribution function.
  // out:       Incremented output.
  double dv0dx0 = dxv[2]/dxv[0]; 
  double w0dx0 = w[2]/dxv[0]; 
  double dv1dx1 = dxv[3]/dxv[1]; 
  double w1dx1 = w[3]/dxv[1]; 
  const double dv10 = 2/dxv[2]; 
  const double *E0 = &field[0]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv11 = 2/dxv[3]; 
  const double *E1 = &field[4]; 
  const double dv2 = dxv[3], wv2 = w[3]; 

  const double *B2 = &field[20]; 
  double cflFreq_mid = 0.0; 
  double alpha_vdim[32] = {0.0}; 

  cflFreq_mid += 3.0*(fabs(w0dx0)+0.5*dv0dx0); 

  cflFreq_mid += 3.0*(fabs(w1dx1)+0.5*dv1dx1); 

  out[1] += 3.464101615137754*f[0]*w0dx0+f[3]*dv0dx0; 
  out[2] += 3.464101615137754*f[0]*w1dx1+f[4]*dv1dx1; 
  out[5] += 3.464101615137754*(f[1]*w1dx1+f[2]*w0dx0)+f[8]*dv1dx1+f[7]*dv0dx0; 
  out[6] += 3.464101615137754*f[3]*w0dx0+(0.8944271909999159*f[16]+f[0])*dv0dx0; 
  out[7] += 3.464101615137754*f[3]*w1dx1+f[10]*dv1dx1; 
  out[8] += 3.464101615137754*f[4]*w0dx0+f[10]*dv0dx0; 
  out[9] += 3.464101615137754*f[4]*w1dx1+(0.8944271909999159*f[24]+f[0])*dv1dx1; 
  out[11] += 3.464101615137754*(f[6]*w1dx1+f[7]*w0dx0)+f[13]*dv1dx1+(0.8944271909999161*f[18]+f[2])*dv0dx0; 
  out[12] += 3.464101615137754*(f[8]*w1dx1+f[9]*w0dx0)+(0.8944271909999161*f[25]+f[1])*dv1dx1+f[14]*dv0dx0; 
  out[13] += 3.464101615137754*f[10]*w0dx0+(0.8944271909999161*f[19]+f[4])*dv0dx0; 
  out[14] += 3.464101615137754*f[10]*w1dx1+(0.8944271909999161*f[27]+f[3])*dv1dx1; 
  out[15] += 3.464101615137754*(f[13]*w1dx1+f[14]*w0dx0)+(0.8944271909999159*f[29]+f[6])*dv1dx1+(0.8944271909999159*f[22]+f[9])*dv0dx0; 
  out[17] += 3.464101615137755*f[16]*w0dx0+0.8944271909999161*f[3]*dv0dx0; 
  out[18] += 3.464101615137755*f[16]*w1dx1+f[19]*dv1dx1; 
  out[20] += 3.464101615137755*(f[17]*w1dx1+f[18]*w0dx0)+f[21]*dv1dx1+0.8944271909999159*f[7]*dv0dx0; 
  out[21] += 3.464101615137755*f[19]*w0dx0+0.8944271909999159*f[10]*dv0dx0; 
  out[22] += 3.464101615137755*f[19]*w1dx1+f[16]*dv1dx1; 
  out[23] += 3.464101615137755*(f[21]*w1dx1+f[22]*w0dx0)+f[17]*dv1dx1+0.8944271909999161*f[14]*dv0dx0; 
  out[25] += 3.464101615137755*f[24]*w0dx0+f[27]*dv0dx0; 
  out[26] += 3.464101615137755*f[24]*w1dx1+0.8944271909999161*f[4]*dv1dx1; 
  out[28] += 3.464101615137755*(f[25]*w1dx1+f[26]*w0dx0)+0.8944271909999159*f[8]*dv1dx1+f[30]*dv0dx0; 
  out[29] += 3.464101615137755*f[27]*w0dx0+f[24]*dv0dx0; 
  out[30] += 3.464101615137755*f[27]*w1dx1+0.8944271909999159*f[10]*dv1dx1; 
  out[31] += 3.464101615137755*(f[29]*w1dx1+f[30]*w0dx0)+0.8944271909999161*f[13]*dv1dx1+f[26]*dv0dx0; 

  alpha_vdim[0] = 2.0*dv10*(B2[0]*wv2+E0[0]); 
  alpha_vdim[1] = 2.0*dv10*(B2[1]*wv2+E0[1]); 
  alpha_vdim[2] = 2.0*dv10*(B2[2]*wv2+E0[2]); 
  alpha_vdim[3] = 0.0; 
  alpha_vdim[4] = 0.5773502691896258*B2[0]*dv10*dv2; 
  alpha_vdim[5] = 2.0*dv10*(B2[3]*wv2+E0[3]); 
  alpha_vdim[6] = 0.0; 
  alpha_vdim[7] = 0.0; 
  alpha_vdim[8] = 0.5773502691896258*B2[1]*dv10*dv2; 
  alpha_vdim[9] = 0.5773502691896258*B2[2]*dv10*dv2; 
  alpha_vdim[10] = 0.0; 
  alpha_vdim[11] = 0.0; 
  alpha_vdim[12] = 0.5773502691896258*B2[3]*dv10*dv2; 
  alpha_vdim[13] = 0.0; 
  alpha_vdim[14] = 0.0; 
  alpha_vdim[15] = 0.0; 
  alpha_vdim[16] = 0.0; 
  alpha_vdim[17] = 0.0; 
  alpha_vdim[18] = 0.0; 
  alpha_vdim[19] = 0.0; 
  alpha_vdim[20] = 0.0; 
  alpha_vdim[21] = 0.0; 
  alpha_vdim[22] = 0.0; 
  alpha_vdim[23] = 0.0; 
  alpha_vdim[24] = 0.0; 
  alpha_vdim[25] = 0.0; 
  alpha_vdim[26] = 0.0; 
  alpha_vdim[27] = 0.0; 
  alpha_vdim[28] = 0.0; 
  alpha_vdim[29] = 0.0; 
  alpha_vdim[30] = 0.0; 
  alpha_vdim[31] = 0.0; 
  cflFreq_mid += 5.0*fabs(0.125*alpha_vdim[0]); 

  out[3] += 0.4330127018922193*(alpha_vdim[12]*f[12]+alpha_vdim[9]*f[9]+alpha_vdim[8]*f[8]+alpha_vdim[5]*f[5]+alpha_vdim[4]*f[4]+alpha_vdim[2]*f[2]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[6] += 0.4330127018922193*(alpha_vdim[9]*f[12]+f[9]*alpha_vdim[12]+alpha_vdim[4]*f[8]+f[4]*alpha_vdim[8]+alpha_vdim[2]*f[5]+f[2]*alpha_vdim[5]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[7] += 0.4330127018922193*(alpha_vdim[8]*f[12]+f[8]*alpha_vdim[12]+alpha_vdim[4]*f[9]+f[4]*alpha_vdim[9]+alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5]+alpha_vdim[0]*f[2]+f[0]*alpha_vdim[2]); 
  out[10] += 0.3872983346207416*alpha_vdim[12]*f[28]+0.3872983346207417*(alpha_vdim[9]*f[26]+alpha_vdim[8]*f[25])+0.3872983346207416*alpha_vdim[4]*f[24]+0.4330127018922193*(alpha_vdim[5]*f[12]+f[5]*alpha_vdim[12]+alpha_vdim[2]*f[9]+f[2]*alpha_vdim[9]+alpha_vdim[1]*f[8]+f[1]*alpha_vdim[8]+alpha_vdim[0]*f[4]+f[0]*alpha_vdim[4]); 
  out[11] += 0.4330127018922193*(alpha_vdim[4]*f[12]+f[4]*alpha_vdim[12]+alpha_vdim[8]*f[9]+f[8]*alpha_vdim[9]+alpha_vdim[0]*f[5]+f[0]*alpha_vdim[5]+alpha_vdim[1]*f[2]+f[1]*alpha_vdim[2]); 
  out[13] += 0.3872983346207416*alpha_vdim[9]*f[28]+0.3872983346207417*(alpha_vdim[12]*f[26]+alpha_vdim[4]*f[25])+0.3872983346207416*alpha_vdim[8]*f[24]+0.4330127018922193*(alpha_vdim[2]*f[12]+f[2]*alpha_vdim[12]+alpha_vdim[5]*f[9]+f[5]*alpha_vdim[9]+alpha_vdim[0]*f[8]+f[0]*alpha_vdim[8]+alpha_vdim[1]*f[4]+f[1]*alpha_vdim[4]); 
  out[14] += 0.3872983346207416*alpha_vdim[8]*f[28]+0.3872983346207417*(alpha_vdim[4]*f[26]+alpha_vdim[12]*f[25])+0.3872983346207416*alpha_vdim[9]*f[24]+0.4330127018922193*(alpha_vdim[1]*f[12]+f[1]*alpha_vdim[12]+alpha_vdim[0]*f[9]+f[0]*alpha_vdim[9]+alpha_vdim[5]*f[8]+f[5]*alpha_vdim[8]+alpha_vdim[2]*f[4]+f[2]*alpha_vdim[4]); 
  out[15] += 0.3872983346207416*alpha_vdim[4]*f[28]+0.3872983346207417*(alpha_vdim[8]*f[26]+alpha_vdim[9]*f[25])+0.3872983346207416*alpha_vdim[12]*f[24]+0.4330127018922193*(alpha_vdim[0]*f[12]+f[0]*alpha_vdim[12]+alpha_vdim[1]*f[9]+f[1]*alpha_vdim[9]+alpha_vdim[2]*f[8]+f[2]*alpha_vdim[8]+alpha_vdim[4]*f[5]+f[4]*alpha_vdim[5]); 
  out[16] += 0.9682458365518543*(alpha_vdim[12]*f[15]+alpha_vdim[9]*f[14]+alpha_vdim[8]*f[13]+alpha_vdim[5]*f[11]+alpha_vdim[4]*f[10]+alpha_vdim[2]*f[7]+alpha_vdim[1]*f[6]+alpha_vdim[0]*f[3]); 
  out[17] += 0.9682458365518543*(alpha_vdim[9]*f[15]+alpha_vdim[12]*f[14]+alpha_vdim[4]*f[13]+alpha_vdim[2]*f[11]+alpha_vdim[8]*f[10]+alpha_vdim[5]*f[7]+alpha_vdim[0]*f[6]+alpha_vdim[1]*f[3]); 
  out[18] += 0.9682458365518543*(alpha_vdim[8]*f[15]+alpha_vdim[4]*f[14]+alpha_vdim[12]*f[13]+alpha_vdim[1]*f[11]+alpha_vdim[9]*f[10]+alpha_vdim[0]*f[7]+alpha_vdim[5]*f[6]+alpha_vdim[2]*f[3]); 
  out[19] += 0.8660254037844386*alpha_vdim[12]*f[31]+0.8660254037844387*(alpha_vdim[9]*f[30]+alpha_vdim[8]*f[29])+0.8660254037844386*alpha_vdim[4]*f[27]+0.9682458365518543*(alpha_vdim[5]*f[15]+alpha_vdim[2]*f[14]+alpha_vdim[1]*f[13]+f[11]*alpha_vdim[12]+alpha_vdim[0]*f[10]+f[7]*alpha_vdim[9]+f[6]*alpha_vdim[8]+f[3]*alpha_vdim[4]); 
  out[20] += 0.9682458365518543*(alpha_vdim[4]*f[15]+alpha_vdim[8]*f[14]+alpha_vdim[9]*f[13]+f[10]*alpha_vdim[12]+alpha_vdim[0]*f[11]+alpha_vdim[1]*f[7]+alpha_vdim[2]*f[6]+f[3]*alpha_vdim[5]); 
  out[21] += 0.8660254037844387*alpha_vdim[9]*f[31]+0.8660254037844386*(alpha_vdim[12]*f[30]+alpha_vdim[4]*f[29])+0.8660254037844387*alpha_vdim[8]*f[27]+0.9682458365518543*(alpha_vdim[2]*f[15]+alpha_vdim[5]*f[14]+alpha_vdim[0]*f[13]+f[7]*alpha_vdim[12]+alpha_vdim[9]*f[11]+alpha_vdim[1]*f[10]+f[3]*alpha_vdim[8]+alpha_vdim[4]*f[6]); 
  out[22] += 0.8660254037844387*alpha_vdim[8]*f[31]+0.8660254037844386*(alpha_vdim[4]*f[30]+alpha_vdim[12]*f[29])+0.8660254037844387*alpha_vdim[9]*f[27]+0.9682458365518543*(alpha_vdim[1]*f[15]+alpha_vdim[0]*f[14]+alpha_vdim[5]*f[13]+f[6]*alpha_vdim[12]+alpha_vdim[8]*f[11]+alpha_vdim[2]*f[10]+f[3]*alpha_vdim[9]+alpha_vdim[4]*f[7]); 
  out[23] += 0.8660254037844386*alpha_vdim[4]*f[31]+0.8660254037844387*(alpha_vdim[8]*f[30]+alpha_vdim[9]*f[29])+0.8660254037844386*alpha_vdim[12]*f[27]+0.9682458365518543*(alpha_vdim[0]*f[15]+alpha_vdim[1]*f[14]+alpha_vdim[2]*f[13]+f[3]*alpha_vdim[12]+alpha_vdim[4]*f[11]+alpha_vdim[5]*f[10]+f[6]*alpha_vdim[9]+f[7]*alpha_vdim[8]); 
  out[27] += 0.4330127018922194*alpha_vdim[5]*f[28]+0.4330127018922193*(alpha_vdim[2]*f[26]+alpha_vdim[1]*f[25])+0.4330127018922194*alpha_vdim[0]*f[24]+0.3872983346207417*(alpha_vdim[12]*f[12]+alpha_vdim[9]*f[9]+alpha_vdim[8]*f[8]+alpha_vdim[4]*f[4]); 
  out[29] += 0.4330127018922193*alpha_vdim[2]*f[28]+0.4330127018922194*(alpha_vdim[5]*f[26]+alpha_vdim[0]*f[25])+0.4330127018922193*alpha_vdim[1]*f[24]+0.3872983346207416*(alpha_vdim[9]*f[12]+f[9]*alpha_vdim[12]+alpha_vdim[4]*f[8]+f[4]*alpha_vdim[8]); 
  out[30] += 0.4330127018922193*alpha_vdim[1]*f[28]+0.4330127018922194*(alpha_vdim[0]*f[26]+alpha_vdim[5]*f[25])+0.4330127018922193*alpha_vdim[2]*f[24]+0.3872983346207416*(alpha_vdim[8]*f[12]+f[8]*alpha_vdim[12]+alpha_vdim[4]*f[9]+f[4]*alpha_vdim[9]); 
  out[31] += 0.4330127018922194*alpha_vdim[0]*f[28]+0.4330127018922193*(alpha_vdim[1]*f[26]+alpha_vdim[2]*f[25])+0.4330127018922194*alpha_vdim[5]*f[24]+0.3872983346207417*(alpha_vdim[4]*f[12]+f[4]*alpha_vdim[12]+alpha_vdim[8]*f[9]+f[8]*alpha_vdim[9]); 

  alpha_vdim[0] = dv11*(2.0*E1[0]-2.0*B2[0]*wv1); 
  alpha_vdim[1] = dv11*(2.0*E1[1]-2.0*B2[1]*wv1); 
  alpha_vdim[2] = dv11*(2.0*E1[2]-2.0*B2[2]*wv1); 
  alpha_vdim[3] = -0.5773502691896258*B2[0]*dv1*dv11; 
  alpha_vdim[4] = 0.0; 
  alpha_vdim[5] = dv11*(2.0*E1[3]-2.0*B2[3]*wv1); 
  alpha_vdim[6] = -0.5773502691896258*B2[1]*dv1*dv11; 
  alpha_vdim[7] = -0.5773502691896258*B2[2]*dv1*dv11; 
  alpha_vdim[8] = 0.0; 
  alpha_vdim[9] = 0.0; 
  alpha_vdim[10] = 0.0; 
  alpha_vdim[11] = -0.5773502691896258*B2[3]*dv1*dv11; 
  alpha_vdim[12] = 0.0; 
  alpha_vdim[13] = 0.0; 
  alpha_vdim[14] = 0.0; 
  alpha_vdim[15] = 0.0; 
  alpha_vdim[16] = 0.0; 
  alpha_vdim[17] = 0.0; 
  alpha_vdim[18] = 0.0; 
  alpha_vdim[19] = 0.0; 
  alpha_vdim[20] = 0.0; 
  alpha_vdim[21] = 0.0; 
  alpha_vdim[22] = 0.0; 
  alpha_vdim[23] = 0.0; 
  alpha_vdim[24] = 0.0; 
  alpha_vdim[25] = 0.0; 
  alpha_vdim[26] = 0.0; 
  alpha_vdim[27] = 0.0; 
  alpha_vdim[28] = 0.0; 
  alpha_vdim[29] = 0.0; 
  alpha_vdim[30] = 0.0; 
  alpha_vdim[31] = 0.0; 
  cflFreq_mid += 5.0*fabs(0.125*alpha_vdim[0]); 

  out[4] += 0.4330127018922193*(alpha_vdim[11]*f[11]+alpha_vdim[7]*f[7]+alpha_vdim[6]*f[6]+alpha_vdim[5]*f[5]+alpha_vdim[3]*f[3]+alpha_vdim[2]*f[2]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[8] += 0.4330127018922193*(alpha_vdim[7]*f[11]+f[7]*alpha_vdim[11]+alpha_vdim[3]*f[6]+f[3]*alpha_vdim[6]+alpha_vdim[2]*f[5]+f[2]*alpha_vdim[5]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[9] += 0.4330127018922193*(alpha_vdim[6]*f[11]+f[6]*alpha_vdim[11]+alpha_vdim[3]*f[7]+f[3]*alpha_vdim[7]+alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5]+alpha_vdim[0]*f[2]+f[0]*alpha_vdim[2]); 
  out[10] += 0.3872983346207416*alpha_vdim[11]*f[20]+0.3872983346207417*(alpha_vdim[7]*f[18]+alpha_vdim[6]*f[17])+0.3872983346207416*alpha_vdim[3]*f[16]+0.4330127018922193*(alpha_vdim[5]*f[11]+f[5]*alpha_vdim[11]+alpha_vdim[2]*f[7]+f[2]*alpha_vdim[7]+alpha_vdim[1]*f[6]+f[1]*alpha_vdim[6]+alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3]); 
  out[12] += 0.4330127018922193*(alpha_vdim[3]*f[11]+f[3]*alpha_vdim[11]+alpha_vdim[6]*f[7]+f[6]*alpha_vdim[7]+alpha_vdim[0]*f[5]+f[0]*alpha_vdim[5]+alpha_vdim[1]*f[2]+f[1]*alpha_vdim[2]); 
  out[13] += 0.3872983346207416*alpha_vdim[7]*f[20]+0.3872983346207417*(alpha_vdim[11]*f[18]+alpha_vdim[3]*f[17])+0.3872983346207416*alpha_vdim[6]*f[16]+0.4330127018922193*(alpha_vdim[2]*f[11]+f[2]*alpha_vdim[11]+alpha_vdim[5]*f[7]+f[5]*alpha_vdim[7]+alpha_vdim[0]*f[6]+f[0]*alpha_vdim[6]+alpha_vdim[1]*f[3]+f[1]*alpha_vdim[3]); 
  out[14] += 0.3872983346207416*alpha_vdim[6]*f[20]+0.3872983346207417*(alpha_vdim[3]*f[18]+alpha_vdim[11]*f[17])+0.3872983346207416*alpha_vdim[7]*f[16]+0.4330127018922193*(alpha_vdim[1]*f[11]+f[1]*alpha_vdim[11]+alpha_vdim[0]*f[7]+f[0]*alpha_vdim[7]+alpha_vdim[5]*f[6]+f[5]*alpha_vdim[6]+alpha_vdim[2]*f[3]+f[2]*alpha_vdim[3]); 
  out[15] += 0.3872983346207416*alpha_vdim[3]*f[20]+0.3872983346207417*(alpha_vdim[6]*f[18]+alpha_vdim[7]*f[17])+0.3872983346207416*alpha_vdim[11]*f[16]+0.4330127018922193*(alpha_vdim[0]*f[11]+f[0]*alpha_vdim[11]+alpha_vdim[1]*f[7]+f[1]*alpha_vdim[7]+alpha_vdim[2]*f[6]+f[2]*alpha_vdim[6]+alpha_vdim[3]*f[5]+f[3]*alpha_vdim[5]); 
  out[19] += 0.4330127018922194*alpha_vdim[5]*f[20]+0.4330127018922193*(alpha_vdim[2]*f[18]+alpha_vdim[1]*f[17])+0.4330127018922194*alpha_vdim[0]*f[16]+0.3872983346207417*(alpha_vdim[11]*f[11]+alpha_vdim[7]*f[7]+alpha_vdim[6]*f[6]+alpha_vdim[3]*f[3]); 
  out[21] += 0.4330127018922193*alpha_vdim[2]*f[20]+0.4330127018922194*(alpha_vdim[5]*f[18]+alpha_vdim[0]*f[17])+0.4330127018922193*alpha_vdim[1]*f[16]+0.3872983346207416*(alpha_vdim[7]*f[11]+f[7]*alpha_vdim[11]+alpha_vdim[3]*f[6]+f[3]*alpha_vdim[6]); 
  out[22] += 0.4330127018922193*alpha_vdim[1]*f[20]+0.4330127018922194*(alpha_vdim[0]*f[18]+alpha_vdim[5]*f[17])+0.4330127018922193*alpha_vdim[2]*f[16]+0.3872983346207416*(alpha_vdim[6]*f[11]+f[6]*alpha_vdim[11]+alpha_vdim[3]*f[7]+f[3]*alpha_vdim[7]); 
  out[23] += 0.4330127018922194*alpha_vdim[0]*f[20]+0.4330127018922193*(alpha_vdim[1]*f[18]+alpha_vdim[2]*f[17])+0.4330127018922194*alpha_vdim[5]*f[16]+0.3872983346207417*(alpha_vdim[3]*f[11]+f[3]*alpha_vdim[11]+alpha_vdim[6]*f[7]+f[6]*alpha_vdim[7]); 
  out[24] += 0.9682458365518543*(alpha_vdim[11]*f[15]+alpha_vdim[7]*f[14]+alpha_vdim[6]*f[13]+alpha_vdim[5]*f[12]+alpha_vdim[3]*f[10]+alpha_vdim[2]*f[9]+alpha_vdim[1]*f[8]+alpha_vdim[0]*f[4]); 
  out[25] += 0.9682458365518543*(alpha_vdim[7]*f[15]+alpha_vdim[11]*f[14]+alpha_vdim[3]*f[13]+alpha_vdim[2]*f[12]+alpha_vdim[6]*f[10]+alpha_vdim[5]*f[9]+alpha_vdim[0]*f[8]+alpha_vdim[1]*f[4]); 
  out[26] += 0.9682458365518543*(alpha_vdim[6]*f[15]+alpha_vdim[3]*f[14]+alpha_vdim[11]*f[13]+alpha_vdim[1]*f[12]+alpha_vdim[7]*f[10]+alpha_vdim[0]*f[9]+alpha_vdim[5]*f[8]+alpha_vdim[2]*f[4]); 
  out[27] += 0.8660254037844386*alpha_vdim[11]*f[23]+0.8660254037844387*(alpha_vdim[7]*f[22]+alpha_vdim[6]*f[21])+0.8660254037844386*alpha_vdim[3]*f[19]+0.9682458365518543*(alpha_vdim[5]*f[15]+alpha_vdim[2]*f[14]+alpha_vdim[1]*f[13]+alpha_vdim[11]*f[12]+alpha_vdim[0]*f[10]+alpha_vdim[7]*f[9]+alpha_vdim[6]*f[8]+alpha_vdim[3]*f[4]); 
  out[28] += 0.9682458365518543*(alpha_vdim[3]*f[15]+alpha_vdim[6]*f[14]+alpha_vdim[7]*f[13]+alpha_vdim[0]*f[12]+f[10]*alpha_vdim[11]+alpha_vdim[1]*f[9]+alpha_vdim[2]*f[8]+f[4]*alpha_vdim[5]); 
  out[29] += 0.8660254037844387*alpha_vdim[7]*f[23]+0.8660254037844386*(alpha_vdim[11]*f[22]+alpha_vdim[3]*f[21])+0.8660254037844387*alpha_vdim[6]*f[19]+0.9682458365518543*(alpha_vdim[2]*f[15]+alpha_vdim[5]*f[14]+alpha_vdim[0]*f[13]+alpha_vdim[7]*f[12]+f[9]*alpha_vdim[11]+alpha_vdim[1]*f[10]+alpha_vdim[3]*f[8]+f[4]*alpha_vdim[6]); 
  out[30] += 0.8660254037844387*alpha_vdim[6]*f[23]+0.8660254037844386*(alpha_vdim[3]*f[22]+alpha_vdim[11]*f[21])+0.8660254037844387*alpha_vdim[7]*f[19]+0.9682458365518543*(alpha_vdim[1]*f[15]+alpha_vdim[0]*f[14]+alpha_vdim[5]*f[13]+alpha_vdim[6]*f[12]+f[8]*alpha_vdim[11]+alpha_vdim[2]*f[10]+alpha_vdim[3]*f[9]+f[4]*alpha_vdim[7]); 
  out[31] += 0.8660254037844386*alpha_vdim[3]*f[23]+0.8660254037844387*(alpha_vdim[6]*f[22]+alpha_vdim[7]*f[21])+0.8660254037844386*alpha_vdim[11]*f[19]+0.9682458365518543*(alpha_vdim[0]*f[15]+alpha_vdim[1]*f[14]+alpha_vdim[2]*f[13]+alpha_vdim[3]*f[12]+f[4]*alpha_vdim[11]+alpha_vdim[5]*f[10]+alpha_vdim[6]*f[9]+alpha_vdim[7]*f[8]); 

  return cflFreq_mid; 
} 
