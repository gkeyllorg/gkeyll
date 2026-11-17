#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order4_surfy_2x_ser_p3_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[1],4.);

  out[0] += -(0.00390625*(3928.940696930917*coeff[1]*qr[9]-3928.940696930917*coeff[1]*ql[9]-6808.826991486861*coeff[1]*qr[5]-6808.826991486861*coeff[1]*ql[5]+13617.653982973721*coeff[1]*qc[5]+6365.286717815623*coeff[1]*qr[2]-6365.286717815623*coeff[1]*ql[2]+(-(3675.0*qr[0])-3675.0*ql[0]+7350.0*qc[0])*coeff[1])*Jfac); 
  out[1] += -(0.00390625*(3928.9406969309152*coeff[1]*qr[11]-3928.9406969309152*coeff[1]*ql[11]-6808.82699148686*coeff[1]*qr[7]-6808.82699148686*coeff[1]*ql[7]+13617.65398297372*coeff[1]*qc[7]+6365.286717815623*coeff[1]*qr[3]-6365.286717815623*coeff[1]*ql[3]-3675.0*coeff[1]*qr[1]-3675.0*coeff[1]*ql[1]+7350.0*coeff[1]*qc[1])*Jfac); 
  out[2] += -(5.580357142857143e-4*(38837.329014750736*coeff[1]*qr[9]+38837.329014750736*coeff[1]*ql[9]+112868.83936676232*coeff[1]*qc[9]-73118.05259304981*coeff[1]*qr[5]+73118.05259304981*coeff[1]*ql[5]+73395.0*coeff[1]*qr[2]+73395.0*coeff[1]*ql[2]+161910.0*coeff[1]*qc[2]+(44557.00702470936*ql[0]-44557.00702470936*qr[0])*coeff[1])*Jfac); 
  out[3] += -(5.580357142857143e-4*(38837.32901475074*coeff[1]*qr[11]+38837.32901475074*coeff[1]*ql[11]+112868.83936676233*coeff[1]*qc[11]-73118.05259304983*coeff[1]*qr[7]+73118.05259304983*coeff[1]*ql[7]+73395.0*coeff[1]*qr[3]+73395.0*coeff[1]*ql[3]+161910.0*coeff[1]*qc[3]-44557.00702470936*coeff[1]*qr[1]+44557.00702470936*coeff[1]*ql[1])*Jfac); 
  out[4] += -(0.00390625*(6365.286717815626*coeff[1]*qr[6]-6365.286717815626*coeff[1]*ql[6]-3675.0*coeff[1]*qr[4]-3675.0*coeff[1]*ql[4]+7350.0*coeff[1]*qc[4])*Jfac); 
  out[5] += -(5.580357142857143e-4*(19221.343215290653*coeff[1]*qr[9]-19221.343215290653*coeff[1]*ql[9]-51345.0*coeff[1]*qr[5]-51345.0*coeff[1]*ql[5]+248850.0*coeff[1]*qc[5]+61731.4815552*coeff[1]*qr[2]-61731.4815552*coeff[1]*ql[2]+(-(41087.749086558644*qr[0])-41087.749086558644*ql[0]+82175.49817311729*qc[0])*coeff[1])*Jfac); 
  out[6] += -(0.00390625*(10485.0*coeff[1]*qr[6]+10485.0*coeff[1]*ql[6]+23130.0*coeff[1]*qc[6]-6365.286717815626*coeff[1]*qr[4]+6365.286717815626*coeff[1]*ql[4])*Jfac); 
  out[7] += -(5.580357142857143e-4*(19221.343215290653*coeff[1]*qr[11]-19221.343215290653*coeff[1]*ql[11]-51345.0*coeff[1]*qr[7]-51345.0*coeff[1]*ql[7]+248850.0*coeff[1]*qc[7]+61731.48155520002*coeff[1]*qr[3]-61731.48155520002*coeff[1]*ql[3]-41087.74908655864*coeff[1]*qr[1]-41087.74908655864*coeff[1]*ql[1]+82175.49817311727*coeff[1]*qc[1])*Jfac); 
  out[8] += -(0.00390625*(6365.286717815624*coeff[1]*qr[10]-6365.286717815624*coeff[1]*ql[10]-3675.0*coeff[1]*qr[8]-3675.0*coeff[1]*ql[8]+7350.0*coeff[1]*qc[8])*Jfac); 
  out[9] += 0.00390625*(6135.0*coeff[1]*qr[9]+6135.0*coeff[1]*ql[9]-26130.0*coeff[1]*qc[9]-6584.596798589874*coeff[1]*qr[5]+6584.596798589874*coeff[1]*ql[5]+3918.102219187242*coeff[1]*qr[2]+3918.102219187242*coeff[1]*ql[2]+3436.9317712168795*coeff[1]*qc[2]+(1627.1370563047233*ql[0]-1627.1370563047233*qr[0])*coeff[1])*Jfac; 
  out[10] += -(0.00390625*(10485.0*coeff[1]*qr[10]+10485.0*coeff[1]*ql[10]+23130.0*coeff[1]*qc[10]-6365.286717815624*coeff[1]*qr[8]+6365.286717815624*coeff[1]*ql[8])*Jfac); 
  out[11] += 0.00390625*(6135.0*coeff[1]*qr[11]+6135.0*coeff[1]*ql[11]-26130.0*coeff[1]*qc[11]-6584.596798589872*coeff[1]*qr[7]+6584.596798589872*coeff[1]*ql[7]+3918.102219187243*coeff[1]*qr[3]+3918.102219187243*coeff[1]*ql[3]+3436.93177121688*coeff[1]*qc[3]-1627.137056304723*coeff[1]*qr[1]+1627.137056304723*coeff[1]*ql[1])*Jfac; 

  return 0.;

}

