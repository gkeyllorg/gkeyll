#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // field:      q/m*EM fields.
  // cot_vec:   Only used in gen geo.
  // f:         Input distribution function.
  // out:       Incremented output.
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  const double dv10 = 2/dxv[1]; 
  const double *E0 = &field[0]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv11 = 2/dxv[2]; 
  const double *E1 = &field[2]; 
  const double dv2 = dxv[2], wv2 = w[2]; 

  const double *B2 = &field[10]; 
  double cflFreq_mid = 0.0; 
  double alpha_vdim[16] = {0.0}; 

  cflFreq_mid += 3.0*(fabs(w0dx0)+0.5*dv0dx0); 

  out[1] += 3.464101615137754*f[0]*w0dx0+f[2]*dv0dx0; 
  out[4] += 3.464101615137754*f[2]*w0dx0+(0.8944271909999159*f[8]+f[0])*dv0dx0; 
  out[5] += 3.464101615137754*f[3]*w0dx0+f[6]*dv0dx0; 
  out[7] += 3.464101615137754*f[6]*w0dx0+(0.8944271909999161*f[10]+f[3])*dv0dx0; 
  out[9] += 3.464101615137755*f[8]*w0dx0+0.8944271909999161*f[2]*dv0dx0; 
  out[11] += 3.464101615137755*f[10]*w0dx0+0.8944271909999159*f[6]*dv0dx0; 
  out[13] += 3.464101615137755*f[12]*w0dx0+f[14]*dv0dx0; 
  out[15] += 3.464101615137755*f[14]*w0dx0+f[12]*dv0dx0; 

  alpha_vdim[0] = 2.0*dv10*(B2[0]*wv2+E0[0]); 
  alpha_vdim[1] = 2.0*dv10*(B2[1]*wv2+E0[1]); 
  alpha_vdim[2] = 0.0; 
  alpha_vdim[3] = 0.5773502691896258*B2[0]*dv10*dv2; 
  alpha_vdim[4] = 0.0; 
  alpha_vdim[5] = 0.5773502691896258*B2[1]*dv10*dv2; 
  alpha_vdim[6] = 0.0; 
  alpha_vdim[7] = 0.0; 
  alpha_vdim[8] = 0.0; 
  alpha_vdim[9] = 0.0; 
  alpha_vdim[10] = 0.0; 
  alpha_vdim[11] = 0.0; 
  alpha_vdim[12] = 0.0; 
  alpha_vdim[13] = 0.0; 
  alpha_vdim[14] = 0.0; 
  alpha_vdim[15] = 0.0; 
  cflFreq_mid += 5.0*fabs(0.1767766952966368*alpha_vdim[0]); 

  out[2] += 0.6123724356957944*(alpha_vdim[5]*f[5]+alpha_vdim[3]*f[3]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[4] += 0.6123724356957944*(alpha_vdim[3]*f[5]+f[3]*alpha_vdim[5]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[6] += 0.5477225575051661*(alpha_vdim[5]*f[13]+alpha_vdim[3]*f[12])+0.6123724356957944*(alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5]+alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3]); 
  out[7] += 0.5477225575051661*(alpha_vdim[3]*f[13]+alpha_vdim[5]*f[12])+0.6123724356957944*(alpha_vdim[0]*f[5]+f[0]*alpha_vdim[5]+alpha_vdim[1]*f[3]+f[1]*alpha_vdim[3]); 
  out[8] += 1.369306393762915*(alpha_vdim[5]*f[7]+alpha_vdim[3]*f[6]+alpha_vdim[1]*f[4]+alpha_vdim[0]*f[2]); 
  out[9] += 1.369306393762915*(alpha_vdim[3]*f[7]+alpha_vdim[5]*f[6]+alpha_vdim[0]*f[4]+alpha_vdim[1]*f[2]); 
  out[10] += 1.224744871391589*(alpha_vdim[5]*f[15]+alpha_vdim[3]*f[14])+1.369306393762915*(alpha_vdim[1]*f[7]+alpha_vdim[0]*f[6]+f[4]*alpha_vdim[5]+f[2]*alpha_vdim[3]); 
  out[11] += 1.224744871391589*(alpha_vdim[3]*f[15]+alpha_vdim[5]*f[14])+1.369306393762915*(alpha_vdim[0]*f[7]+alpha_vdim[1]*f[6]+f[2]*alpha_vdim[5]+alpha_vdim[3]*f[4]); 
  out[14] += 0.6123724356957944*(alpha_vdim[1]*f[13]+alpha_vdim[0]*f[12])+0.5477225575051661*(alpha_vdim[5]*f[5]+alpha_vdim[3]*f[3]); 
  out[15] += 0.6123724356957944*(alpha_vdim[0]*f[13]+alpha_vdim[1]*f[12])+0.5477225575051661*(alpha_vdim[3]*f[5]+f[3]*alpha_vdim[5]); 

  alpha_vdim[0] = dv11*(2.0*E1[0]-2.0*B2[0]*wv1); 
  alpha_vdim[1] = dv11*(2.0*E1[1]-2.0*B2[1]*wv1); 
  alpha_vdim[2] = -0.5773502691896258*B2[0]*dv1*dv11; 
  alpha_vdim[3] = 0.0; 
  alpha_vdim[4] = -0.5773502691896258*B2[1]*dv1*dv11; 
  alpha_vdim[5] = 0.0; 
  alpha_vdim[6] = 0.0; 
  alpha_vdim[7] = 0.0; 
  alpha_vdim[8] = 0.0; 
  alpha_vdim[9] = 0.0; 
  alpha_vdim[10] = 0.0; 
  alpha_vdim[11] = 0.0; 
  alpha_vdim[12] = 0.0; 
  alpha_vdim[13] = 0.0; 
  alpha_vdim[14] = 0.0; 
  alpha_vdim[15] = 0.0; 
  cflFreq_mid += 5.0*fabs(0.1767766952966368*alpha_vdim[0]); 

  out[3] += 0.6123724356957944*(alpha_vdim[4]*f[4]+alpha_vdim[2]*f[2]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[5] += 0.6123724356957944*(alpha_vdim[2]*f[4]+f[2]*alpha_vdim[4]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[6] += 0.5477225575051661*(alpha_vdim[4]*f[9]+alpha_vdim[2]*f[8])+0.6123724356957944*(alpha_vdim[1]*f[4]+f[1]*alpha_vdim[4]+alpha_vdim[0]*f[2]+f[0]*alpha_vdim[2]); 
  out[7] += 0.5477225575051661*(alpha_vdim[2]*f[9]+alpha_vdim[4]*f[8])+0.6123724356957944*(alpha_vdim[0]*f[4]+f[0]*alpha_vdim[4]+alpha_vdim[1]*f[2]+f[1]*alpha_vdim[2]); 
  out[10] += 0.6123724356957944*(alpha_vdim[1]*f[9]+alpha_vdim[0]*f[8])+0.5477225575051661*(alpha_vdim[4]*f[4]+alpha_vdim[2]*f[2]); 
  out[11] += 0.6123724356957944*(alpha_vdim[0]*f[9]+alpha_vdim[1]*f[8])+0.5477225575051661*(alpha_vdim[2]*f[4]+f[2]*alpha_vdim[4]); 
  out[12] += 1.369306393762915*(alpha_vdim[4]*f[7]+alpha_vdim[2]*f[6]+alpha_vdim[1]*f[5]+alpha_vdim[0]*f[3]); 
  out[13] += 1.369306393762915*(alpha_vdim[2]*f[7]+alpha_vdim[4]*f[6]+alpha_vdim[0]*f[5]+alpha_vdim[1]*f[3]); 
  out[14] += 1.224744871391589*(alpha_vdim[4]*f[11]+alpha_vdim[2]*f[10])+1.369306393762915*(alpha_vdim[1]*f[7]+alpha_vdim[0]*f[6]+alpha_vdim[4]*f[5]+alpha_vdim[2]*f[3]); 
  out[15] += 1.224744871391589*(alpha_vdim[2]*f[11]+alpha_vdim[4]*f[10])+1.369306393762915*(alpha_vdim[0]*f[7]+alpha_vdim[1]*f[6]+alpha_vdim[2]*f[5]+f[3]*alpha_vdim[4]); 

  return cflFreq_mid; 
} 
