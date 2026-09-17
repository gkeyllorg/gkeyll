#include <gkyl_vlasov_sr_kernels.h> 
GKYL_CU_DH double vlasov_sr_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out) 
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
  double p0_over_gamma[3] = {0.0}; 
  p0_over_gamma[0] = 1.732050807568877*gamma[1]*dv10; 
  p0_over_gamma[1] = 3.872983346207417*gamma[2]*dv10; 

  double cflFreq_mid = 0.0; 
  double alpha_vdim[6] = {0.0}; 

  cflFreq_mid += 3.0*fabs(0.3535533905932737*p0_over_gamma[0]*dx10); 

  out[1] += 1.224744871391589*(p0_over_gamma[1]*f[2]+f[0]*p0_over_gamma[0])*dx10; 
  out[3] += (1.095445115010332*p0_over_gamma[1]*f[4]+1.224744871391589*(p0_over_gamma[0]*f[2]+f[0]*p0_over_gamma[1]))*dx10; 
  out[5] += (1.224744871391589*p0_over_gamma[0]*f[4]+1.095445115010332*p0_over_gamma[1]*f[2])*dx10; 

  alpha_vdim[0] = 1.414213562373095*E0[0]*dv10; 
  alpha_vdim[1] = 1.414213562373095*E0[1]*dv10; 
  alpha_vdim[2] = 0.0; 
  alpha_vdim[3] = 0.0; 
  alpha_vdim[4] = 0.0; 
  alpha_vdim[5] = 0.0; 
  cflFreq_mid += 5.0*fabs(0.25*alpha_vdim[0]); 

  out[2] += 0.8660254037844386*(alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.8660254037844386*(alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[4] += 1.936491673103709*(alpha_vdim[1]*f[3]+alpha_vdim[0]*f[2]); 
  out[5] += 1.936491673103709*(alpha_vdim[0]*f[3]+alpha_vdim[1]*f[2]); 

  return cflFreq_mid; 
} 
