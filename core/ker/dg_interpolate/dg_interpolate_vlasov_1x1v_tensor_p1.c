#include <gkyl_dg_interpolate_kernels.h> 
 
GKYL_CU_DH void dg_interpolate_vlasov_1x1v_tensor_p1_x(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar) 
{ 
  // wDo: cell center of donor cell.
  // wTar: cell center of target cell.
  // dxDo: cell length of donor cell.
  // dxTar: cell length of target cell.
  // fldDo: donor field.
  // fldTar: target field in cells pointed to by the stencil.

  double eLo = fmax(-1.0,1.0-(2.0*(wTar[0]-1.0*wDo[0]+0.5*dxTar[0]+0.5*dxDo[0]))/dxTar[0]);
  double eUp = fmin( 1.0,(2.0*(-(1.0*wTar[0])+wDo[0]+0.5*dxTar[0]+0.5*dxDo[0]))/dxTar[0]-1.0);

#ifdef __CUDA_ARCH__
  const double eLoR2 = pow(eLo,2);
  const double eLoR3 = pow(eLo,3);
  const double eUpR2 = pow(eUp,2);
  const double eUpR3 = pow(eUp,3);

  atomicAdd(&fldTar[0], (dxTar[0]*(0.4330127018922193*fldDo[1]*eUpR2-0.4330127018922193*fldDo[1]*eLoR2))/dxDo[0]+(1.7320508075688772*wTar[0]*fldDo[1]*eUp-1.7320508075688772*wDo[0]*fldDo[1]*eUp-1.7320508075688772*wTar[0]*fldDo[1]*eLo+1.7320508075688772*wDo[0]*fldDo[1]*eLo)/dxDo[0]+0.5*fldDo[0]*eUp-0.5*fldDo[0]*eLo); 
  atomicAdd(&fldTar[1], (dxTar[0]*(0.5*fldDo[1]*eUpR3-0.5*fldDo[1]*eLoR3))/dxDo[0]+(1.5*wTar[0]*fldDo[1]*eUpR2-1.5*wDo[0]*fldDo[1]*eUpR2-1.5*wTar[0]*fldDo[1]*eLoR2+1.5*wDo[0]*fldDo[1]*eLoR2)/dxDo[0]+0.4330127018922193*fldDo[0]*eUpR2-0.4330127018922193*fldDo[0]*eLoR2); 
  atomicAdd(&fldTar[2], (dxTar[0]*(0.4330127018922193*fldDo[3]*eUpR2-0.4330127018922193*fldDo[3]*eLoR2))/dxDo[0]+(1.7320508075688772*wTar[0]*fldDo[3]*eUp-1.7320508075688772*wDo[0]*fldDo[3]*eUp-1.7320508075688772*wTar[0]*fldDo[3]*eLo+1.7320508075688772*wDo[0]*fldDo[3]*eLo)/dxDo[0]+0.5*fldDo[2]*eUp-0.5*fldDo[2]*eLo); 
  atomicAdd(&fldTar[3], (dxTar[0]*(0.5*fldDo[3]*eUpR3-0.5*fldDo[3]*eLoR3))/dxDo[0]+(1.5*wTar[0]*fldDo[3]*eUpR2-1.5*wDo[0]*fldDo[3]*eUpR2-1.5*wTar[0]*fldDo[3]*eLoR2+1.5*wDo[0]*fldDo[3]*eLoR2)/dxDo[0]+0.4330127018922193*fldDo[2]*eUpR2-0.4330127018922193*fldDo[2]*eLoR2); 
  atomicAdd(&fldTar[4], (dxTar[0]*(0.43301270189221935*fldDo[5]*eUpR2-0.43301270189221935*fldDo[5]*eLoR2))/dxDo[0]+(1.7320508075688774*wTar[0]*fldDo[5]*eUp-1.7320508075688774*wDo[0]*fldDo[5]*eUp-1.7320508075688774*wTar[0]*fldDo[5]*eLo+1.7320508075688774*wDo[0]*fldDo[5]*eLo)/dxDo[0]+0.5*fldDo[4]*eUp-0.5*fldDo[4]*eLo); 
  atomicAdd(&fldTar[5], (dxTar[0]*(0.5*fldDo[5]*eUpR3-0.5*fldDo[5]*eLoR3))/dxDo[0]+(1.5*wTar[0]*fldDo[5]*eUpR2-1.5*wDo[0]*fldDo[5]*eUpR2-1.5*wTar[0]*fldDo[5]*eLoR2+1.5*wDo[0]*fldDo[5]*eLoR2)/dxDo[0]+0.43301270189221935*fldDo[4]*eUpR2-0.43301270189221935*fldDo[4]*eLoR2); 
#else
  const double eLoR2 = pow(eLo,2);
  const double eLoR3 = pow(eLo,3);
  const double eUpR2 = pow(eUp,2);
  const double eUpR3 = pow(eUp,3);

  fldTar[0] += (dxTar[0]*(0.4330127018922193*fldDo[1]*eUpR2-0.4330127018922193*fldDo[1]*eLoR2))/dxDo[0]+(1.7320508075688772*wTar[0]*fldDo[1]*eUp-1.7320508075688772*wDo[0]*fldDo[1]*eUp-1.7320508075688772*wTar[0]*fldDo[1]*eLo+1.7320508075688772*wDo[0]*fldDo[1]*eLo)/dxDo[0]+0.5*fldDo[0]*eUp-0.5*fldDo[0]*eLo; 
  fldTar[1] += (dxTar[0]*(0.5*fldDo[1]*eUpR3-0.5*fldDo[1]*eLoR3))/dxDo[0]+(1.5*wTar[0]*fldDo[1]*eUpR2-1.5*wDo[0]*fldDo[1]*eUpR2-1.5*wTar[0]*fldDo[1]*eLoR2+1.5*wDo[0]*fldDo[1]*eLoR2)/dxDo[0]+0.4330127018922193*fldDo[0]*eUpR2-0.4330127018922193*fldDo[0]*eLoR2; 
  fldTar[2] += (dxTar[0]*(0.4330127018922193*fldDo[3]*eUpR2-0.4330127018922193*fldDo[3]*eLoR2))/dxDo[0]+(1.7320508075688772*wTar[0]*fldDo[3]*eUp-1.7320508075688772*wDo[0]*fldDo[3]*eUp-1.7320508075688772*wTar[0]*fldDo[3]*eLo+1.7320508075688772*wDo[0]*fldDo[3]*eLo)/dxDo[0]+0.5*fldDo[2]*eUp-0.5*fldDo[2]*eLo; 
  fldTar[3] += (dxTar[0]*(0.5*fldDo[3]*eUpR3-0.5*fldDo[3]*eLoR3))/dxDo[0]+(1.5*wTar[0]*fldDo[3]*eUpR2-1.5*wDo[0]*fldDo[3]*eUpR2-1.5*wTar[0]*fldDo[3]*eLoR2+1.5*wDo[0]*fldDo[3]*eLoR2)/dxDo[0]+0.4330127018922193*fldDo[2]*eUpR2-0.4330127018922193*fldDo[2]*eLoR2; 
  fldTar[4] += (dxTar[0]*(0.43301270189221935*fldDo[5]*eUpR2-0.43301270189221935*fldDo[5]*eLoR2))/dxDo[0]+(1.7320508075688774*wTar[0]*fldDo[5]*eUp-1.7320508075688774*wDo[0]*fldDo[5]*eUp-1.7320508075688774*wTar[0]*fldDo[5]*eLo+1.7320508075688774*wDo[0]*fldDo[5]*eLo)/dxDo[0]+0.5*fldDo[4]*eUp-0.5*fldDo[4]*eLo; 
  fldTar[5] += (dxTar[0]*(0.5*fldDo[5]*eUpR3-0.5*fldDo[5]*eLoR3))/dxDo[0]+(1.5*wTar[0]*fldDo[5]*eUpR2-1.5*wDo[0]*fldDo[5]*eUpR2-1.5*wTar[0]*fldDo[5]*eLoR2+1.5*wDo[0]*fldDo[5]*eLoR2)/dxDo[0]+0.43301270189221935*fldDo[4]*eUpR2-0.43301270189221935*fldDo[4]*eLoR2; 
#endif

}

GKYL_CU_DH void dg_interpolate_vlasov_1x1v_tensor_p1_vx(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar) 
{ 
  // wDo: cell center of donor cell.
  // wTar: cell center of target cell.
  // dxDo: cell length of donor cell.
  // dxTar: cell length of target cell.
  // fldDo: donor field.
  // fldTar: target field in cells pointed to by the stencil.

  double eLo = fmax(-1.0,1.0-(2.0*(wTar[1]-1.0*wDo[1]+0.5*dxTar[1]+0.5*dxDo[1]))/dxTar[1]);
  double eUp = fmin( 1.0,(2.0*(-(1.0*wTar[1])+wDo[1]+0.5*dxTar[1]+0.5*dxDo[1]))/dxTar[1]-1.0);

#ifdef __CUDA_ARCH__
  const double wDo1R2 = pow(wDo[1],2);
  const double wTar1R2 = pow(wTar[1],2);
  const double dxDo1R2 = pow(dxDo[1],2);
  const double dxTar1R2 = pow(dxTar[1],2);
  const double eLoR2 = pow(eLo,2);
  const double eLoR3 = pow(eLo,3);
  const double eLoR4 = pow(eLo,4);
  const double eLoR5 = pow(eLo,5);
  const double eUpR2 = pow(eUp,2);
  const double eUpR3 = pow(eUp,3);
  const double eUpR4 = pow(eUp,4);
  const double eUpR5 = pow(eUp,5);

  atomicAdd(&fldTar[0], (dxTar1R2*(0.5590169943749475*fldDo[4]*eUpR3-0.5590169943749475*fldDo[4]*eLoR3))/dxDo1R2+(dxTar[1]*(3.3541019662496847*wTar[1]*fldDo[4]*eUpR2-3.3541019662496847*wDo[1]*fldDo[4]*eUpR2-3.3541019662496847*wTar[1]*fldDo[4]*eLoR2+3.3541019662496847*wDo[1]*fldDo[4]*eLoR2))/dxDo1R2+(dxTar[1]*(0.4330127018922193*fldDo[2]*eUpR2-0.4330127018922193*fldDo[2]*eLoR2))/dxDo[1]+(6.708203932499369*wTar1R2*fldDo[4]*eUp-13.416407864998739*wDo[1]*wTar[1]*fldDo[4]*eUp+6.708203932499369*wDo1R2*fldDo[4]*eUp-6.708203932499369*wTar1R2*fldDo[4]*eLo+13.416407864998739*wDo[1]*wTar[1]*fldDo[4]*eLo-6.708203932499369*wDo1R2*fldDo[4]*eLo)/dxDo1R2+(1.7320508075688772*wTar[1]*fldDo[2]*eUp-1.7320508075688772*wDo[1]*fldDo[2]*eUp-1.7320508075688772*wTar[1]*fldDo[2]*eLo+1.7320508075688772*wDo[1]*fldDo[2]*eLo)/dxDo[1]-0.5590169943749475*fldDo[4]*eUp+0.5*fldDo[0]*eUp+0.5590169943749475*fldDo[4]*eLo-0.5*fldDo[0]*eLo); 
  atomicAdd(&fldTar[1], (dxTar1R2*(0.5590169943749476*fldDo[5]*eUpR3-0.5590169943749476*fldDo[5]*eLoR3))/dxDo1R2+(dxTar[1]*(3.3541019662496843*wTar[1]*fldDo[5]*eUpR2-3.3541019662496843*wDo[1]*fldDo[5]*eUpR2-3.3541019662496843*wTar[1]*fldDo[5]*eLoR2+3.3541019662496843*wDo[1]*fldDo[5]*eLoR2))/dxDo1R2+(dxTar[1]*(0.4330127018922193*fldDo[3]*eUpR2-0.4330127018922193*fldDo[3]*eLoR2))/dxDo[1]+(6.7082039324993685*wTar1R2*fldDo[5]*eUp-13.416407864998737*wDo[1]*wTar[1]*fldDo[5]*eUp+6.7082039324993685*wDo1R2*fldDo[5]*eUp-6.7082039324993685*wTar1R2*fldDo[5]*eLo+13.416407864998737*wDo[1]*wTar[1]*fldDo[5]*eLo-6.7082039324993685*wDo1R2*fldDo[5]*eLo)/dxDo1R2+(1.7320508075688772*wTar[1]*fldDo[3]*eUp-1.7320508075688772*wDo[1]*fldDo[3]*eUp-1.7320508075688772*wTar[1]*fldDo[3]*eLo+1.7320508075688772*wDo[1]*fldDo[3]*eLo)/dxDo[1]-0.5590169943749476*fldDo[5]*eUp+0.5*fldDo[1]*eUp+0.5590169943749476*fldDo[5]*eLo-0.5*fldDo[1]*eLo); 
  atomicAdd(&fldTar[2], (dxTar1R2*(0.7261843774138906*fldDo[4]*eUpR4-0.7261843774138906*fldDo[4]*eLoR4))/dxDo1R2+(dxTar[1]*(3.872983346207417*wTar[1]*fldDo[4]*eUpR3-3.872983346207417*wDo[1]*fldDo[4]*eUpR3-3.872983346207417*wTar[1]*fldDo[4]*eLoR3+3.872983346207417*wDo[1]*fldDo[4]*eLoR3))/dxDo1R2+(dxTar[1]*(0.5*fldDo[2]*eUpR3-0.5*fldDo[2]*eLoR3))/dxDo[1]+(5.809475019311125*wTar1R2*fldDo[4]*eUpR2-11.61895003862225*wDo[1]*wTar[1]*fldDo[4]*eUpR2+5.809475019311125*wDo1R2*fldDo[4]*eUpR2-5.809475019311125*wTar1R2*fldDo[4]*eLoR2+11.61895003862225*wDo[1]*wTar[1]*fldDo[4]*eLoR2-5.809475019311125*wDo1R2*fldDo[4]*eLoR2)/dxDo1R2+(1.5*wTar[1]*fldDo[2]*eUpR2-1.5*wDo[1]*fldDo[2]*eUpR2-1.5*wTar[1]*fldDo[2]*eLoR2+1.5*wDo[1]*fldDo[2]*eLoR2)/dxDo[1]-0.4841229182759271*fldDo[4]*eUpR2+0.4330127018922193*fldDo[0]*eUpR2+0.4841229182759271*fldDo[4]*eLoR2-0.4330127018922193*fldDo[0]*eLoR2); 
  atomicAdd(&fldTar[3], (dxTar1R2*(0.7261843774138907*fldDo[5]*eUpR4-0.7261843774138907*fldDo[5]*eLoR4))/dxDo1R2+(dxTar[1]*(3.872983346207417*wTar[1]*fldDo[5]*eUpR3-3.872983346207417*wDo[1]*fldDo[5]*eUpR3-3.872983346207417*wTar[1]*fldDo[5]*eLoR3+3.872983346207417*wDo[1]*fldDo[5]*eLoR3))/dxDo1R2+(dxTar[1]*(0.5*fldDo[3]*eUpR3-0.5*fldDo[3]*eLoR3))/dxDo[1]+(5.809475019311126*wTar1R2*fldDo[5]*eUpR2-11.618950038622252*wDo[1]*wTar[1]*fldDo[5]*eUpR2+5.809475019311126*wDo1R2*fldDo[5]*eUpR2-5.809475019311126*wTar1R2*fldDo[5]*eLoR2+11.618950038622252*wDo[1]*wTar[1]*fldDo[5]*eLoR2-5.809475019311126*wDo1R2*fldDo[5]*eLoR2)/dxDo1R2+(1.5*wTar[1]*fldDo[3]*eUpR2-1.5*wDo[1]*fldDo[3]*eUpR2-1.5*wTar[1]*fldDo[3]*eLoR2+1.5*wDo[1]*fldDo[3]*eLoR2)/dxDo[1]-0.4841229182759271*fldDo[5]*eUpR2+0.4330127018922193*fldDo[1]*eUpR2+0.4841229182759271*fldDo[5]*eLoR2-0.4330127018922193*fldDo[1]*eLoR2); 
  atomicAdd(&fldTar[4], (dxTar1R2*(1.125*fldDo[4]*eUpR5-0.625*fldDo[4]*eUpR3-1.125*fldDo[4]*eLoR5+0.625*fldDo[4]*eLoR3))/dxDo1R2+(dxTar[1]*(5.625*wTar[1]*fldDo[4]*eUpR4-5.625*wDo[1]*fldDo[4]*eUpR4-3.75*wTar[1]*fldDo[4]*eUpR2+3.75*wDo[1]*fldDo[4]*eUpR2-5.625*wTar[1]*fldDo[4]*eLoR4+5.625*wDo[1]*fldDo[4]*eLoR4+3.75*wTar[1]*fldDo[4]*eLoR2-3.75*wDo[1]*fldDo[4]*eLoR2))/dxDo1R2+(dxTar[1]*(0.7261843774138906*fldDo[2]*eUpR4-0.4841229182759271*fldDo[2]*eUpR2-0.7261843774138906*fldDo[2]*eLoR4+0.4841229182759271*fldDo[2]*eLoR2))/dxDo[1]+(7.5*wTar1R2*fldDo[4]*eUpR3-15.0*wDo[1]*wTar[1]*fldDo[4]*eUpR3+7.5*wDo1R2*fldDo[4]*eUpR3-7.5*wTar1R2*fldDo[4]*eUp+15.0*wDo[1]*wTar[1]*fldDo[4]*eUp-7.5*wDo1R2*fldDo[4]*eUp-7.5*wTar1R2*fldDo[4]*eLoR3+15.0*wDo[1]*wTar[1]*fldDo[4]*eLoR3-7.5*wDo1R2*fldDo[4]*eLoR3+7.5*wTar1R2*fldDo[4]*eLo-15.0*wDo[1]*wTar[1]*fldDo[4]*eLo+7.5*wDo1R2*fldDo[4]*eLo)/dxDo1R2+(1.9364916731037085*wTar[1]*fldDo[2]*eUpR3-1.9364916731037085*wDo[1]*fldDo[2]*eUpR3-1.9364916731037085*wTar[1]*fldDo[2]*eUp+1.9364916731037085*wDo[1]*fldDo[2]*eUp-1.9364916731037085*wTar[1]*fldDo[2]*eLoR3+1.9364916731037085*wDo[1]*fldDo[2]*eLoR3+1.9364916731037085*wTar[1]*fldDo[2]*eLo-1.9364916731037085*wDo[1]*fldDo[2]*eLo)/dxDo[1]-0.625*fldDo[4]*eUpR3+0.5590169943749475*fldDo[0]*eUpR3+0.625*fldDo[4]*eUp-0.5590169943749475*fldDo[0]*eUp+0.625*fldDo[4]*eLoR3-0.5590169943749475*fldDo[0]*eLoR3-0.625*fldDo[4]*eLo+0.5590169943749475*fldDo[0]*eLo); 
  atomicAdd(&fldTar[5], (dxTar1R2*(1.125*fldDo[5]*eUpR5-0.625*fldDo[5]*eUpR3-1.125*fldDo[5]*eLoR5+0.625*fldDo[5]*eLoR3))/dxDo1R2+(dxTar[1]*(5.625*wTar[1]*fldDo[5]*eUpR4-5.625*wDo[1]*fldDo[5]*eUpR4-3.75*wTar[1]*fldDo[5]*eUpR2+3.75*wDo[1]*fldDo[5]*eUpR2-5.625*wTar[1]*fldDo[5]*eLoR4+5.625*wDo[1]*fldDo[5]*eLoR4+3.75*wTar[1]*fldDo[5]*eLoR2-3.75*wDo[1]*fldDo[5]*eLoR2))/dxDo1R2+(dxTar[1]*(0.7261843774138907*fldDo[3]*eUpR4-0.4841229182759271*fldDo[3]*eUpR2-0.7261843774138907*fldDo[3]*eLoR4+0.4841229182759271*fldDo[3]*eLoR2))/dxDo[1]+(7.5*wTar1R2*fldDo[5]*eUpR3-15.0*wDo[1]*wTar[1]*fldDo[5]*eUpR3+7.5*wDo1R2*fldDo[5]*eUpR3-7.5*wTar1R2*fldDo[5]*eUp+15.0*wDo[1]*wTar[1]*fldDo[5]*eUp-7.5*wDo1R2*fldDo[5]*eUp-7.5*wTar1R2*fldDo[5]*eLoR3+15.0*wDo[1]*wTar[1]*fldDo[5]*eLoR3-7.5*wDo1R2*fldDo[5]*eLoR3+7.5*wTar1R2*fldDo[5]*eLo-15.0*wDo[1]*wTar[1]*fldDo[5]*eLo+7.5*wDo1R2*fldDo[5]*eLo)/dxDo1R2+(1.9364916731037085*wTar[1]*fldDo[3]*eUpR3-1.9364916731037085*wDo[1]*fldDo[3]*eUpR3-1.9364916731037085*wTar[1]*fldDo[3]*eUp+1.9364916731037085*wDo[1]*fldDo[3]*eUp-1.9364916731037085*wTar[1]*fldDo[3]*eLoR3+1.9364916731037085*wDo[1]*fldDo[3]*eLoR3+1.9364916731037085*wTar[1]*fldDo[3]*eLo-1.9364916731037085*wDo[1]*fldDo[3]*eLo)/dxDo[1]-0.625*fldDo[5]*eUpR3+0.5590169943749476*fldDo[1]*eUpR3+0.625*fldDo[5]*eUp-0.5590169943749476*fldDo[1]*eUp+0.625*fldDo[5]*eLoR3-0.5590169943749476*fldDo[1]*eLoR3-0.625*fldDo[5]*eLo+0.5590169943749476*fldDo[1]*eLo); 
#else
  const double wDo1R2 = pow(wDo[1],2);
  const double wTar1R2 = pow(wTar[1],2);
  const double dxDo1R2 = pow(dxDo[1],2);
  const double dxTar1R2 = pow(dxTar[1],2);
  const double eLoR2 = pow(eLo,2);
  const double eLoR3 = pow(eLo,3);
  const double eLoR4 = pow(eLo,4);
  const double eLoR5 = pow(eLo,5);
  const double eUpR2 = pow(eUp,2);
  const double eUpR3 = pow(eUp,3);
  const double eUpR4 = pow(eUp,4);
  const double eUpR5 = pow(eUp,5);

  fldTar[0] += (dxTar1R2*(0.5590169943749475*fldDo[4]*eUpR3-0.5590169943749475*fldDo[4]*eLoR3))/dxDo1R2+(dxTar[1]*(3.3541019662496847*wTar[1]*fldDo[4]*eUpR2-3.3541019662496847*wDo[1]*fldDo[4]*eUpR2-3.3541019662496847*wTar[1]*fldDo[4]*eLoR2+3.3541019662496847*wDo[1]*fldDo[4]*eLoR2))/dxDo1R2+(dxTar[1]*(0.4330127018922193*fldDo[2]*eUpR2-0.4330127018922193*fldDo[2]*eLoR2))/dxDo[1]+(6.708203932499369*wTar1R2*fldDo[4]*eUp-13.416407864998739*wDo[1]*wTar[1]*fldDo[4]*eUp+6.708203932499369*wDo1R2*fldDo[4]*eUp-6.708203932499369*wTar1R2*fldDo[4]*eLo+13.416407864998739*wDo[1]*wTar[1]*fldDo[4]*eLo-6.708203932499369*wDo1R2*fldDo[4]*eLo)/dxDo1R2+(1.7320508075688772*wTar[1]*fldDo[2]*eUp-1.7320508075688772*wDo[1]*fldDo[2]*eUp-1.7320508075688772*wTar[1]*fldDo[2]*eLo+1.7320508075688772*wDo[1]*fldDo[2]*eLo)/dxDo[1]-0.5590169943749475*fldDo[4]*eUp+0.5*fldDo[0]*eUp+0.5590169943749475*fldDo[4]*eLo-0.5*fldDo[0]*eLo; 
  fldTar[1] += (dxTar1R2*(0.5590169943749476*fldDo[5]*eUpR3-0.5590169943749476*fldDo[5]*eLoR3))/dxDo1R2+(dxTar[1]*(3.3541019662496843*wTar[1]*fldDo[5]*eUpR2-3.3541019662496843*wDo[1]*fldDo[5]*eUpR2-3.3541019662496843*wTar[1]*fldDo[5]*eLoR2+3.3541019662496843*wDo[1]*fldDo[5]*eLoR2))/dxDo1R2+(dxTar[1]*(0.4330127018922193*fldDo[3]*eUpR2-0.4330127018922193*fldDo[3]*eLoR2))/dxDo[1]+(6.7082039324993685*wTar1R2*fldDo[5]*eUp-13.416407864998737*wDo[1]*wTar[1]*fldDo[5]*eUp+6.7082039324993685*wDo1R2*fldDo[5]*eUp-6.7082039324993685*wTar1R2*fldDo[5]*eLo+13.416407864998737*wDo[1]*wTar[1]*fldDo[5]*eLo-6.7082039324993685*wDo1R2*fldDo[5]*eLo)/dxDo1R2+(1.7320508075688772*wTar[1]*fldDo[3]*eUp-1.7320508075688772*wDo[1]*fldDo[3]*eUp-1.7320508075688772*wTar[1]*fldDo[3]*eLo+1.7320508075688772*wDo[1]*fldDo[3]*eLo)/dxDo[1]-0.5590169943749476*fldDo[5]*eUp+0.5*fldDo[1]*eUp+0.5590169943749476*fldDo[5]*eLo-0.5*fldDo[1]*eLo; 
  fldTar[2] += (dxTar1R2*(0.7261843774138906*fldDo[4]*eUpR4-0.7261843774138906*fldDo[4]*eLoR4))/dxDo1R2+(dxTar[1]*(3.872983346207417*wTar[1]*fldDo[4]*eUpR3-3.872983346207417*wDo[1]*fldDo[4]*eUpR3-3.872983346207417*wTar[1]*fldDo[4]*eLoR3+3.872983346207417*wDo[1]*fldDo[4]*eLoR3))/dxDo1R2+(dxTar[1]*(0.5*fldDo[2]*eUpR3-0.5*fldDo[2]*eLoR3))/dxDo[1]+(5.809475019311125*wTar1R2*fldDo[4]*eUpR2-11.61895003862225*wDo[1]*wTar[1]*fldDo[4]*eUpR2+5.809475019311125*wDo1R2*fldDo[4]*eUpR2-5.809475019311125*wTar1R2*fldDo[4]*eLoR2+11.61895003862225*wDo[1]*wTar[1]*fldDo[4]*eLoR2-5.809475019311125*wDo1R2*fldDo[4]*eLoR2)/dxDo1R2+(1.5*wTar[1]*fldDo[2]*eUpR2-1.5*wDo[1]*fldDo[2]*eUpR2-1.5*wTar[1]*fldDo[2]*eLoR2+1.5*wDo[1]*fldDo[2]*eLoR2)/dxDo[1]-0.4841229182759271*fldDo[4]*eUpR2+0.4330127018922193*fldDo[0]*eUpR2+0.4841229182759271*fldDo[4]*eLoR2-0.4330127018922193*fldDo[0]*eLoR2; 
  fldTar[3] += (dxTar1R2*(0.7261843774138907*fldDo[5]*eUpR4-0.7261843774138907*fldDo[5]*eLoR4))/dxDo1R2+(dxTar[1]*(3.872983346207417*wTar[1]*fldDo[5]*eUpR3-3.872983346207417*wDo[1]*fldDo[5]*eUpR3-3.872983346207417*wTar[1]*fldDo[5]*eLoR3+3.872983346207417*wDo[1]*fldDo[5]*eLoR3))/dxDo1R2+(dxTar[1]*(0.5*fldDo[3]*eUpR3-0.5*fldDo[3]*eLoR3))/dxDo[1]+(5.809475019311126*wTar1R2*fldDo[5]*eUpR2-11.618950038622252*wDo[1]*wTar[1]*fldDo[5]*eUpR2+5.809475019311126*wDo1R2*fldDo[5]*eUpR2-5.809475019311126*wTar1R2*fldDo[5]*eLoR2+11.618950038622252*wDo[1]*wTar[1]*fldDo[5]*eLoR2-5.809475019311126*wDo1R2*fldDo[5]*eLoR2)/dxDo1R2+(1.5*wTar[1]*fldDo[3]*eUpR2-1.5*wDo[1]*fldDo[3]*eUpR2-1.5*wTar[1]*fldDo[3]*eLoR2+1.5*wDo[1]*fldDo[3]*eLoR2)/dxDo[1]-0.4841229182759271*fldDo[5]*eUpR2+0.4330127018922193*fldDo[1]*eUpR2+0.4841229182759271*fldDo[5]*eLoR2-0.4330127018922193*fldDo[1]*eLoR2; 
  fldTar[4] += (dxTar1R2*(1.125*fldDo[4]*eUpR5-0.625*fldDo[4]*eUpR3-1.125*fldDo[4]*eLoR5+0.625*fldDo[4]*eLoR3))/dxDo1R2+(dxTar[1]*(5.625*wTar[1]*fldDo[4]*eUpR4-5.625*wDo[1]*fldDo[4]*eUpR4-3.75*wTar[1]*fldDo[4]*eUpR2+3.75*wDo[1]*fldDo[4]*eUpR2-5.625*wTar[1]*fldDo[4]*eLoR4+5.625*wDo[1]*fldDo[4]*eLoR4+3.75*wTar[1]*fldDo[4]*eLoR2-3.75*wDo[1]*fldDo[4]*eLoR2))/dxDo1R2+(dxTar[1]*(0.7261843774138906*fldDo[2]*eUpR4-0.4841229182759271*fldDo[2]*eUpR2-0.7261843774138906*fldDo[2]*eLoR4+0.4841229182759271*fldDo[2]*eLoR2))/dxDo[1]+(7.5*wTar1R2*fldDo[4]*eUpR3-15.0*wDo[1]*wTar[1]*fldDo[4]*eUpR3+7.5*wDo1R2*fldDo[4]*eUpR3-7.5*wTar1R2*fldDo[4]*eUp+15.0*wDo[1]*wTar[1]*fldDo[4]*eUp-7.5*wDo1R2*fldDo[4]*eUp-7.5*wTar1R2*fldDo[4]*eLoR3+15.0*wDo[1]*wTar[1]*fldDo[4]*eLoR3-7.5*wDo1R2*fldDo[4]*eLoR3+7.5*wTar1R2*fldDo[4]*eLo-15.0*wDo[1]*wTar[1]*fldDo[4]*eLo+7.5*wDo1R2*fldDo[4]*eLo)/dxDo1R2+(1.9364916731037085*wTar[1]*fldDo[2]*eUpR3-1.9364916731037085*wDo[1]*fldDo[2]*eUpR3-1.9364916731037085*wTar[1]*fldDo[2]*eUp+1.9364916731037085*wDo[1]*fldDo[2]*eUp-1.9364916731037085*wTar[1]*fldDo[2]*eLoR3+1.9364916731037085*wDo[1]*fldDo[2]*eLoR3+1.9364916731037085*wTar[1]*fldDo[2]*eLo-1.9364916731037085*wDo[1]*fldDo[2]*eLo)/dxDo[1]-0.625*fldDo[4]*eUpR3+0.5590169943749475*fldDo[0]*eUpR3+0.625*fldDo[4]*eUp-0.5590169943749475*fldDo[0]*eUp+0.625*fldDo[4]*eLoR3-0.5590169943749475*fldDo[0]*eLoR3-0.625*fldDo[4]*eLo+0.5590169943749475*fldDo[0]*eLo; 
  fldTar[5] += (dxTar1R2*(1.125*fldDo[5]*eUpR5-0.625*fldDo[5]*eUpR3-1.125*fldDo[5]*eLoR5+0.625*fldDo[5]*eLoR3))/dxDo1R2+(dxTar[1]*(5.625*wTar[1]*fldDo[5]*eUpR4-5.625*wDo[1]*fldDo[5]*eUpR4-3.75*wTar[1]*fldDo[5]*eUpR2+3.75*wDo[1]*fldDo[5]*eUpR2-5.625*wTar[1]*fldDo[5]*eLoR4+5.625*wDo[1]*fldDo[5]*eLoR4+3.75*wTar[1]*fldDo[5]*eLoR2-3.75*wDo[1]*fldDo[5]*eLoR2))/dxDo1R2+(dxTar[1]*(0.7261843774138907*fldDo[3]*eUpR4-0.4841229182759271*fldDo[3]*eUpR2-0.7261843774138907*fldDo[3]*eLoR4+0.4841229182759271*fldDo[3]*eLoR2))/dxDo[1]+(7.5*wTar1R2*fldDo[5]*eUpR3-15.0*wDo[1]*wTar[1]*fldDo[5]*eUpR3+7.5*wDo1R2*fldDo[5]*eUpR3-7.5*wTar1R2*fldDo[5]*eUp+15.0*wDo[1]*wTar[1]*fldDo[5]*eUp-7.5*wDo1R2*fldDo[5]*eUp-7.5*wTar1R2*fldDo[5]*eLoR3+15.0*wDo[1]*wTar[1]*fldDo[5]*eLoR3-7.5*wDo1R2*fldDo[5]*eLoR3+7.5*wTar1R2*fldDo[5]*eLo-15.0*wDo[1]*wTar[1]*fldDo[5]*eLo+7.5*wDo1R2*fldDo[5]*eLo)/dxDo1R2+(1.9364916731037085*wTar[1]*fldDo[3]*eUpR3-1.9364916731037085*wDo[1]*fldDo[3]*eUpR3-1.9364916731037085*wTar[1]*fldDo[3]*eUp+1.9364916731037085*wDo[1]*fldDo[3]*eUp-1.9364916731037085*wTar[1]*fldDo[3]*eLoR3+1.9364916731037085*wDo[1]*fldDo[3]*eLoR3+1.9364916731037085*wTar[1]*fldDo[3]*eLo-1.9364916731037085*wDo[1]*fldDo[3]*eLo)/dxDo[1]-0.625*fldDo[5]*eUpR3+0.5590169943749476*fldDo[1]*eUpR3+0.625*fldDo[5]*eUp-0.5590169943749476*fldDo[1]*eUp+0.625*fldDo[5]*eLoR3-0.5590169943749476*fldDo[1]*eLoR3-0.625*fldDo[5]*eLo+0.5590169943749476*fldDo[1]*eLo; 
#endif

}

