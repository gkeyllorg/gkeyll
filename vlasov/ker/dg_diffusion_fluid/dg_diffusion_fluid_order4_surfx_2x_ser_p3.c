#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order4_surfx_2x_ser_p3_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],4.);

  out[0] += -(0.00390625*(3928.940696930917*coeff[0]*qr[8]-3928.940696930917*coeff[0]*ql[8]-6808.826991486861*coeff[0]*qr[4]-6808.826991486861*coeff[0]*ql[4]+13617.653982973721*coeff[0]*qc[4]+6365.286717815623*coeff[0]*qr[1]-6365.286717815623*coeff[0]*ql[1]-3675.0*coeff[0]*qr[0]-3675.0*coeff[0]*ql[0]+7350.0*coeff[0]*qc[0])*Jfac); 
  out[1] += -(5.580357142857143e-4*(38837.329014750736*coeff[0]*qr[8]+38837.329014750736*coeff[0]*ql[8]+112868.83936676232*coeff[0]*qc[8]-73118.05259304981*coeff[0]*qr[4]+73118.05259304981*coeff[0]*ql[4]+73395.0*coeff[0]*qr[1]+73395.0*coeff[0]*ql[1]+161910.0*coeff[0]*qc[1]-44557.00702470936*coeff[0]*qr[0]+44557.00702470936*coeff[0]*ql[0])*Jfac); 
  out[2] += -(0.00390625*(3928.9406969309152*coeff[0]*qr[10]-3928.9406969309152*coeff[0]*ql[10]-6808.82699148686*coeff[0]*qr[6]-6808.82699148686*coeff[0]*ql[6]+13617.65398297372*coeff[0]*qc[6]+6365.286717815623*coeff[0]*qr[3]-6365.286717815623*coeff[0]*ql[3]-3675.0*coeff[0]*qr[2]-3675.0*coeff[0]*ql[2]+7350.0*coeff[0]*qc[2])*Jfac); 
  out[3] += -(5.580357142857143e-4*(38837.32901475074*coeff[0]*qr[10]+38837.32901475074*coeff[0]*ql[10]+112868.83936676233*coeff[0]*qc[10]-73118.05259304983*coeff[0]*qr[6]+73118.05259304983*coeff[0]*ql[6]+73395.0*coeff[0]*qr[3]+73395.0*coeff[0]*ql[3]+161910.0*coeff[0]*qc[3]-44557.00702470936*coeff[0]*qr[2]+44557.00702470936*coeff[0]*ql[2])*Jfac); 
  out[4] += -(5.580357142857143e-4*(19221.343215290653*coeff[0]*qr[8]-19221.343215290653*coeff[0]*ql[8]-51345.0*coeff[0]*qr[4]-51345.0*coeff[0]*ql[4]+248850.0*coeff[0]*qc[4]+61731.4815552*coeff[0]*qr[1]-61731.4815552*coeff[0]*ql[1]-41087.749086558644*coeff[0]*qr[0]-41087.749086558644*coeff[0]*ql[0]+82175.49817311729*coeff[0]*qc[0])*Jfac); 
  out[5] += -(0.00390625*(6365.286717815626*coeff[0]*qr[7]-6365.286717815626*coeff[0]*ql[7]-3675.0*coeff[0]*qr[5]-3675.0*coeff[0]*ql[5]+7350.0*coeff[0]*qc[5])*Jfac); 
  out[6] += -(5.580357142857143e-4*(19221.343215290653*coeff[0]*qr[10]-19221.343215290653*coeff[0]*ql[10]-51345.0*coeff[0]*qr[6]-51345.0*coeff[0]*ql[6]+248850.0*coeff[0]*qc[6]+61731.48155520002*coeff[0]*qr[3]-61731.48155520002*coeff[0]*ql[3]-41087.74908655864*coeff[0]*qr[2]-41087.74908655864*coeff[0]*ql[2]+82175.49817311727*coeff[0]*qc[2])*Jfac); 
  out[7] += -(0.00390625*(10485.0*coeff[0]*qr[7]+10485.0*coeff[0]*ql[7]+23130.0*coeff[0]*qc[7]-6365.286717815626*coeff[0]*qr[5]+6365.286717815626*coeff[0]*ql[5])*Jfac); 
  out[8] += 0.00390625*(6135.0*coeff[0]*qr[8]+6135.0*coeff[0]*ql[8]-26130.0*coeff[0]*qc[8]-6584.596798589874*coeff[0]*qr[4]+6584.596798589874*coeff[0]*ql[4]+3918.102219187242*coeff[0]*qr[1]+3918.102219187242*coeff[0]*ql[1]+3436.9317712168795*coeff[0]*qc[1]-1627.1370563047233*coeff[0]*qr[0]+1627.1370563047233*coeff[0]*ql[0])*Jfac; 
  out[9] += -(0.00390625*(6365.286717815624*coeff[0]*qr[11]-6365.286717815624*coeff[0]*ql[11]-3675.0*coeff[0]*qr[9]-3675.0*coeff[0]*ql[9]+7350.0*coeff[0]*qc[9])*Jfac); 
  out[10] += 0.00390625*(6135.0*coeff[0]*qr[10]+6135.0*coeff[0]*ql[10]-26130.0*coeff[0]*qc[10]-6584.596798589872*coeff[0]*qr[6]+6584.596798589872*coeff[0]*ql[6]+3918.102219187243*coeff[0]*qr[3]+3918.102219187243*coeff[0]*ql[3]+3436.93177121688*coeff[0]*qc[3]-1627.137056304723*coeff[0]*qr[2]+1627.137056304723*coeff[0]*ql[2])*Jfac; 
  out[11] += -(0.00390625*(10485.0*coeff[0]*qr[11]+10485.0*coeff[0]*ql[11]+23130.0*coeff[0]*qc[11]-6365.286717815624*coeff[0]*qr[9]+6365.286717815624*coeff[0]*ql[9])*Jfac); 

  return 0.;

}

