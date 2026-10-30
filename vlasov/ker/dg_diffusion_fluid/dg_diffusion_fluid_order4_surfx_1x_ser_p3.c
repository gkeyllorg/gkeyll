#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order4_surfx_1x_ser_p3_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],4.);

  out[0] += -(0.00390625*(3928.940696930917*coeff[0]*qr[3]-3928.940696930917*coeff[0]*ql[3]-6808.826991486861*coeff[0]*qr[2]-6808.826991486861*coeff[0]*ql[2]+13617.653982973721*coeff[0]*qc[2]+6365.286717815623*coeff[0]*qr[1]-6365.286717815623*coeff[0]*ql[1]-3675.0*coeff[0]*qr[0]-3675.0*coeff[0]*ql[0]+7350.0*coeff[0]*qc[0])*Jfac); 
  out[1] += -(5.580357142857143e-4*(38837.329014750736*coeff[0]*qr[3]+38837.329014750736*coeff[0]*ql[3]+112868.83936676232*coeff[0]*qc[3]-73118.05259304981*coeff[0]*qr[2]+73118.05259304981*coeff[0]*ql[2]+73395.0*coeff[0]*qr[1]+73395.0*coeff[0]*ql[1]+161910.0*coeff[0]*qc[1]-44557.00702470936*coeff[0]*qr[0]+44557.00702470936*coeff[0]*ql[0])*Jfac); 
  out[2] += -(5.580357142857143e-4*(19221.343215290653*coeff[0]*qr[3]-19221.343215290653*coeff[0]*ql[3]-51345.0*coeff[0]*qr[2]-51345.0*coeff[0]*ql[2]+248850.0*coeff[0]*qc[2]+61731.4815552*coeff[0]*qr[1]-61731.4815552*coeff[0]*ql[1]-41087.749086558644*coeff[0]*qr[0]-41087.749086558644*coeff[0]*ql[0]+82175.49817311729*coeff[0]*qc[0])*Jfac); 
  out[3] += 0.00390625*(6135.0*coeff[0]*qr[3]+6135.0*coeff[0]*ql[3]-26130.0*coeff[0]*qc[3]-6584.596798589874*coeff[0]*qr[2]+6584.596798589874*coeff[0]*ql[2]+3918.102219187242*coeff[0]*qr[1]+3918.102219187242*coeff[0]*ql[1]+3436.9317712168795*coeff[0]*qc[1]-1627.1370563047233*coeff[0]*qr[0]+1627.1370563047233*coeff[0]*ql[0])*Jfac; 

  return 0.;

}

