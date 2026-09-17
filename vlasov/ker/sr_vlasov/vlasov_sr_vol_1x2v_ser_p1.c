#include <gkyl_vlasov_sr_kernels.h> 
GKYL_CU_DH double vlasov_sr_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // gamma:     Particle Lorentz boost factor sqrt(1 + p^2).
  // qmem:      q/m*EM fields.
  // f:         Input distribution function.
  // out:       Incremented output.
  const double dx10 = 2.0/dxv[0]; 
  const double dv10 = 2.0/dxv[1]; 
  const double *E0 = &qmem[0]; 
  const double dv11 = 2.0/dxv[2]; 
  const double *E1 = &qmem[2]; 
  const double *B2 = &qmem[10]; 
  double p0_over_gamma[8] = {0.0}; 
  p0_over_gamma[0] = 1.732050807568877*gamma[1]*dv10; 
  p0_over_gamma[1] = 3.872983346207417*gamma[4]*dv10; 
  p0_over_gamma[2] = 1.732050807568877*gamma[3]*dv10; 
  p0_over_gamma[3] = 3.872983346207417*gamma[6]*dv10; 
  p0_over_gamma[5] = 1.732050807568877*gamma[7]*dv10; 

  double p1_over_gamma[8] = {0.0}; 
  p1_over_gamma[0] = 1.732050807568877*gamma[2]*dv11; 
  p1_over_gamma[1] = 1.732050807568877*gamma[3]*dv11; 
  p1_over_gamma[2] = 3.872983346207417*gamma[5]*dv11; 
  p1_over_gamma[3] = 3.872983346207417*gamma[7]*dv11; 
  p1_over_gamma[4] = 1.732050807568877*gamma[6]*dv11; 

  double cflFreq_mid = 0.0; 
  double alpha_vdim[16] = {0.0}; 

  cflFreq_mid += 3.0*fabs((0.25*p0_over_gamma[0]-0.2795084971874737*p0_over_gamma[5])*dx10); 

  out[1] += 0.8660254037844386*(p0_over_gamma[5]*f[12]+p0_over_gamma[3]*f[6]+p0_over_gamma[2]*f[3]+p0_over_gamma[1]*f[2]+f[0]*p0_over_gamma[0])*dx10; 
  out[4] += (0.8660254037844387*p0_over_gamma[5]*f[14]+0.7745966692414834*p0_over_gamma[3]*f[10]+0.7745966692414833*p0_over_gamma[1]*f[8]+0.8660254037844386*(p0_over_gamma[2]*f[6]+f[3]*p0_over_gamma[3]+p0_over_gamma[0]*f[2]+f[0]*p0_over_gamma[1]))*dx10; 
  out[5] += (0.7745966692414834*p0_over_gamma[3]*f[14]+0.7745966692414833*p0_over_gamma[2]*f[12]+0.8660254037844386*p0_over_gamma[1]*f[6]+0.7745966692414833*f[3]*p0_over_gamma[5]+0.8660254037844386*(f[2]*p0_over_gamma[3]+p0_over_gamma[0]*f[3]+f[0]*p0_over_gamma[2]))*dx10; 
  out[7] += (0.7745966692414834*p0_over_gamma[2]*f[14]+0.7745966692414833*p0_over_gamma[3]*f[12]+0.7745966692414834*p0_over_gamma[1]*f[10]+0.7745966692414833*(p0_over_gamma[3]*f[8]+p0_over_gamma[5]*f[6])+0.8660254037844386*(p0_over_gamma[0]*f[6]+f[0]*p0_over_gamma[3]+p0_over_gamma[1]*f[3]+f[2]*p0_over_gamma[2]))*dx10; 
  out[9] += (0.8660254037844386*p0_over_gamma[2]*f[10]+0.8660254037844387*p0_over_gamma[0]*f[8]+0.7745966692414834*(p0_over_gamma[3]*f[6]+p0_over_gamma[1]*f[2]))*dx10; 
  out[11] += (0.6928203230275508*p0_over_gamma[3]*f[14]+(0.7745966692414834*p0_over_gamma[5]+0.8660254037844387*p0_over_gamma[0])*f[10]+0.8660254037844386*p0_over_gamma[2]*f[8]+0.7745966692414833*(p0_over_gamma[1]*f[6]+f[2]*p0_over_gamma[3]))*dx10; 
  out[13] += (0.8660254037844386*p0_over_gamma[1]*f[14]+(0.5532833351724881*p0_over_gamma[5]+0.8660254037844387*p0_over_gamma[0])*f[12]+0.7745966692414834*p0_over_gamma[3]*f[6]+0.8660254037844387*f[0]*p0_over_gamma[5]+0.7745966692414834*p0_over_gamma[2]*f[3])*dx10; 
  out[15] += ((0.5532833351724881*p0_over_gamma[5]+0.8660254037844387*p0_over_gamma[0])*f[14]+0.8660254037844386*p0_over_gamma[1]*f[12]+0.6928203230275508*p0_over_gamma[3]*f[10]+0.7745966692414833*p0_over_gamma[2]*f[6]+0.8660254037844386*f[2]*p0_over_gamma[5]+0.7745966692414833*f[3]*p0_over_gamma[3])*dx10; 

  alpha_vdim[0] = (B2[0]*p1_over_gamma[0]+2.0*E0[0])*dv10; 
  alpha_vdim[1] = (2.0*E0[1]+p1_over_gamma[0]*B2[1])*dv10; 
  alpha_vdim[2] = B2[0]*p1_over_gamma[1]*dv10; 
  alpha_vdim[3] = B2[0]*p1_over_gamma[2]*dv10; 
  alpha_vdim[4] = B2[1]*p1_over_gamma[1]*dv10; 
  alpha_vdim[5] = B2[1]*p1_over_gamma[2]*dv10; 
  alpha_vdim[6] = B2[0]*p1_over_gamma[3]*dv10; 
  alpha_vdim[7] = B2[1]*p1_over_gamma[3]*dv10; 
  alpha_vdim[8] = B2[0]*p1_over_gamma[4]*dv10; 
  alpha_vdim[9] = B2[1]*p1_over_gamma[4]*dv10; 
  alpha_vdim[10] = 0.0; 
  alpha_vdim[11] = 0.0; 
  alpha_vdim[12] = 0.0; 
  alpha_vdim[13] = 0.0; 
  alpha_vdim[14] = 0.0; 
  alpha_vdim[15] = 0.0; 
  cflFreq_mid += 5.0*fabs(0.1767766952966368*alpha_vdim[0]-0.1976423537605236*alpha_vdim[8]); 

  out[2] += 0.6123724356957944*(alpha_vdim[9]*f[9]+alpha_vdim[8]*f[8]+alpha_vdim[7]*f[7]+alpha_vdim[6]*f[6]+alpha_vdim[5]*f[5]+alpha_vdim[4]*f[4]+alpha_vdim[3]*f[3]+alpha_vdim[2]*f[2]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[4] += 0.6123724356957944*(alpha_vdim[8]*f[9]+f[8]*alpha_vdim[9]+alpha_vdim[6]*f[7]+f[6]*alpha_vdim[7]+alpha_vdim[3]*f[5]+f[3]*alpha_vdim[5]+alpha_vdim[2]*f[4]+f[2]*alpha_vdim[4]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[6] += 0.5477225575051661*(alpha_vdim[7]*f[15]+alpha_vdim[6]*f[14]+alpha_vdim[5]*f[13]+alpha_vdim[3]*f[12])+0.6123724356957944*(alpha_vdim[9]*f[11]+alpha_vdim[8]*f[10]+alpha_vdim[4]*f[7]+f[4]*alpha_vdim[7]+alpha_vdim[2]*f[6]+f[2]*alpha_vdim[6]+alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5]+alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3]); 
  out[7] += 0.5477225575051661*(alpha_vdim[6]*f[15]+alpha_vdim[7]*f[14]+alpha_vdim[3]*f[13]+alpha_vdim[5]*f[12])+0.6123724356957944*(alpha_vdim[8]*f[11]+alpha_vdim[9]*f[10]+alpha_vdim[2]*f[7]+f[2]*alpha_vdim[7]+alpha_vdim[4]*f[6]+f[4]*alpha_vdim[6]+alpha_vdim[0]*f[5]+f[0]*alpha_vdim[5]+alpha_vdim[1]*f[3]+f[1]*alpha_vdim[3]); 
  out[8] += 1.224744871391589*(alpha_vdim[7]*f[11]+alpha_vdim[6]*f[10]+alpha_vdim[4]*f[9]+f[4]*alpha_vdim[9]+alpha_vdim[2]*f[8]+f[2]*alpha_vdim[8])+1.369306393762915*(alpha_vdim[5]*f[7]+f[5]*alpha_vdim[7]+alpha_vdim[3]*f[6]+f[3]*alpha_vdim[6]+alpha_vdim[1]*f[4]+f[1]*alpha_vdim[4]+alpha_vdim[0]*f[2]+f[0]*alpha_vdim[2]); 
  out[9] += 1.224744871391589*(alpha_vdim[6]*f[11]+alpha_vdim[7]*f[10]+alpha_vdim[2]*f[9]+f[2]*alpha_vdim[9]+alpha_vdim[4]*f[8]+f[4]*alpha_vdim[8])+1.369306393762915*(alpha_vdim[3]*f[7]+f[3]*alpha_vdim[7]+alpha_vdim[5]*f[6]+f[5]*alpha_vdim[6]+alpha_vdim[0]*f[4]+f[0]*alpha_vdim[4]+alpha_vdim[1]*f[2]+f[1]*alpha_vdim[2]); 
  out[10] += 1.224744871391589*(alpha_vdim[5]*f[15]+alpha_vdim[3]*f[14]+alpha_vdim[7]*f[13]+alpha_vdim[6]*f[12]+alpha_vdim[4]*f[11]+alpha_vdim[2]*f[10]+alpha_vdim[7]*f[9]+f[7]*alpha_vdim[9]+alpha_vdim[6]*f[8]+f[6]*alpha_vdim[8])+1.369306393762915*(alpha_vdim[1]*f[7]+f[1]*alpha_vdim[7]+alpha_vdim[0]*f[6]+f[0]*alpha_vdim[6]+alpha_vdim[4]*f[5]+f[4]*alpha_vdim[5]+alpha_vdim[2]*f[3]+f[2]*alpha_vdim[3]); 
  out[11] += 1.224744871391589*(alpha_vdim[3]*f[15]+alpha_vdim[5]*f[14]+alpha_vdim[6]*f[13]+alpha_vdim[7]*f[12]+alpha_vdim[2]*f[11]+alpha_vdim[4]*f[10]+alpha_vdim[6]*f[9]+f[6]*alpha_vdim[9]+alpha_vdim[7]*f[8]+f[7]*alpha_vdim[8])+1.369306393762915*(alpha_vdim[0]*f[7]+f[0]*alpha_vdim[7]+alpha_vdim[1]*f[6]+f[1]*alpha_vdim[6]+alpha_vdim[2]*f[5]+f[2]*alpha_vdim[5]+alpha_vdim[3]*f[4]+f[3]*alpha_vdim[4]); 
  out[14] += 0.6123724356957944*(alpha_vdim[4]*f[15]+alpha_vdim[2]*f[14]+alpha_vdim[1]*f[13]+alpha_vdim[0]*f[12])+0.5477225575051661*(alpha_vdim[7]*f[7]+alpha_vdim[6]*f[6]+alpha_vdim[5]*f[5]+alpha_vdim[3]*f[3]); 
  out[15] += 0.6123724356957944*(alpha_vdim[2]*f[15]+alpha_vdim[4]*f[14]+alpha_vdim[0]*f[13]+alpha_vdim[1]*f[12])+0.5477225575051661*(alpha_vdim[6]*f[7]+f[6]*alpha_vdim[7]+alpha_vdim[3]*f[5]+f[3]*alpha_vdim[5]); 

  alpha_vdim[0] = (2.0*E1[0]-1.0*B2[0]*p0_over_gamma[0])*dv11; 
  alpha_vdim[1] = (2.0*E1[1]-1.0*p0_over_gamma[0]*B2[1])*dv11; 
  alpha_vdim[2] = -1.0*B2[0]*p0_over_gamma[1]*dv11; 
  alpha_vdim[3] = -1.0*B2[0]*p0_over_gamma[2]*dv11; 
  alpha_vdim[4] = -1.0*B2[1]*p0_over_gamma[1]*dv11; 
  alpha_vdim[5] = -1.0*B2[1]*p0_over_gamma[2]*dv11; 
  alpha_vdim[6] = -1.0*B2[0]*p0_over_gamma[3]*dv11; 
  alpha_vdim[7] = -1.0*B2[1]*p0_over_gamma[3]*dv11; 
  alpha_vdim[8] = 0.0; 
  alpha_vdim[9] = 0.0; 
  alpha_vdim[10] = 0.0; 
  alpha_vdim[11] = 0.0; 
  alpha_vdim[12] = -1.0*B2[0]*p0_over_gamma[5]*dv11; 
  alpha_vdim[13] = -1.0*B2[1]*p0_over_gamma[5]*dv11; 
  alpha_vdim[14] = 0.0; 
  alpha_vdim[15] = 0.0; 
  cflFreq_mid += 5.0*fabs(0.1767766952966368*alpha_vdim[0]-0.1976423537605236*alpha_vdim[12]); 

  out[3] += 0.6123724356957944*(alpha_vdim[13]*f[13]+alpha_vdim[12]*f[12]+alpha_vdim[7]*f[7]+alpha_vdim[6]*f[6]+alpha_vdim[5]*f[5]+alpha_vdim[4]*f[4]+alpha_vdim[3]*f[3]+alpha_vdim[2]*f[2]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[5] += 0.6123724356957944*(alpha_vdim[12]*f[13]+f[12]*alpha_vdim[13]+alpha_vdim[6]*f[7]+f[6]*alpha_vdim[7]+alpha_vdim[3]*f[5]+f[3]*alpha_vdim[5]+alpha_vdim[2]*f[4]+f[2]*alpha_vdim[4]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[6] += 0.6123724356957944*(alpha_vdim[13]*f[15]+alpha_vdim[12]*f[14])+0.5477225575051661*(alpha_vdim[7]*f[11]+alpha_vdim[6]*f[10]+alpha_vdim[4]*f[9]+alpha_vdim[2]*f[8])+0.6123724356957944*(alpha_vdim[5]*f[7]+f[5]*alpha_vdim[7]+alpha_vdim[3]*f[6]+f[3]*alpha_vdim[6]+alpha_vdim[1]*f[4]+f[1]*alpha_vdim[4]+alpha_vdim[0]*f[2]+f[0]*alpha_vdim[2]); 
  out[7] += 0.6123724356957944*(alpha_vdim[12]*f[15]+alpha_vdim[13]*f[14])+0.5477225575051661*(alpha_vdim[6]*f[11]+alpha_vdim[7]*f[10]+alpha_vdim[2]*f[9]+alpha_vdim[4]*f[8])+0.6123724356957944*(alpha_vdim[3]*f[7]+f[3]*alpha_vdim[7]+alpha_vdim[5]*f[6]+f[5]*alpha_vdim[6]+alpha_vdim[0]*f[4]+f[0]*alpha_vdim[4]+alpha_vdim[1]*f[2]+f[1]*alpha_vdim[2]); 
  out[10] += 0.6123724356957944*(alpha_vdim[5]*f[11]+alpha_vdim[3]*f[10]+alpha_vdim[1]*f[9]+alpha_vdim[0]*f[8])+0.5477225575051661*(alpha_vdim[7]*f[7]+alpha_vdim[6]*f[6]+alpha_vdim[4]*f[4]+alpha_vdim[2]*f[2]); 
  out[11] += 0.6123724356957944*(alpha_vdim[3]*f[11]+alpha_vdim[5]*f[10]+alpha_vdim[0]*f[9]+alpha_vdim[1]*f[8])+0.5477225575051661*(alpha_vdim[6]*f[7]+f[6]*alpha_vdim[7]+alpha_vdim[2]*f[4]+f[2]*alpha_vdim[4]); 
  out[12] += 1.224744871391589*(alpha_vdim[7]*f[15]+alpha_vdim[6]*f[14]+alpha_vdim[5]*f[13]+f[5]*alpha_vdim[13]+alpha_vdim[3]*f[12]+f[3]*alpha_vdim[12])+1.369306393762915*(alpha_vdim[4]*f[7]+f[4]*alpha_vdim[7]+alpha_vdim[2]*f[6]+f[2]*alpha_vdim[6]+alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5]+alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3]); 
  out[13] += 1.224744871391589*(alpha_vdim[6]*f[15]+alpha_vdim[7]*f[14]+alpha_vdim[3]*f[13]+f[3]*alpha_vdim[13]+alpha_vdim[5]*f[12]+f[5]*alpha_vdim[12])+1.369306393762915*(alpha_vdim[2]*f[7]+f[2]*alpha_vdim[7]+alpha_vdim[4]*f[6]+f[4]*alpha_vdim[6]+alpha_vdim[0]*f[5]+f[0]*alpha_vdim[5]+alpha_vdim[1]*f[3]+f[1]*alpha_vdim[3]); 
  out[14] += 1.224744871391589*(alpha_vdim[5]*f[15]+alpha_vdim[3]*f[14]+alpha_vdim[7]*f[13]+f[7]*alpha_vdim[13]+alpha_vdim[6]*f[12]+f[6]*alpha_vdim[12]+alpha_vdim[4]*f[11]+alpha_vdim[2]*f[10]+alpha_vdim[7]*f[9]+alpha_vdim[6]*f[8])+1.369306393762915*(alpha_vdim[1]*f[7]+f[1]*alpha_vdim[7]+alpha_vdim[0]*f[6]+f[0]*alpha_vdim[6]+alpha_vdim[4]*f[5]+f[4]*alpha_vdim[5]+alpha_vdim[2]*f[3]+f[2]*alpha_vdim[3]); 
  out[15] += 1.224744871391589*(alpha_vdim[3]*f[15]+alpha_vdim[5]*f[14]+alpha_vdim[6]*f[13]+f[6]*alpha_vdim[13]+alpha_vdim[7]*f[12]+f[7]*alpha_vdim[12]+alpha_vdim[2]*f[11]+alpha_vdim[4]*f[10]+alpha_vdim[6]*f[9]+alpha_vdim[7]*f[8])+1.369306393762915*(alpha_vdim[0]*f[7]+f[0]*alpha_vdim[7]+alpha_vdim[1]*f[6]+f[1]*alpha_vdim[6]+alpha_vdim[2]*f[5]+f[2]*alpha_vdim[5]+alpha_vdim[3]*f[4]+f[3]*alpha_vdim[4]); 

  return cflFreq_mid; 
} 
