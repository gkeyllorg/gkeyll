#include <gkyl_dg_interpolate_kernels.h> 
 
GKYL_CU_DH void dg_interpolate_vlasov_2x2v_ser_p1_x(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar) 
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
  atomicAdd(&fldTar[2], (dxTar[0]*(0.4330127018922193*fldDo[5]*eUpR2-0.4330127018922193*fldDo[5]*eLoR2))/dxDo[0]+(1.7320508075688772*wTar[0]*fldDo[5]*eUp-1.7320508075688772*wDo[0]*fldDo[5]*eUp-1.7320508075688772*wTar[0]*fldDo[5]*eLo+1.7320508075688772*wDo[0]*fldDo[5]*eLo)/dxDo[0]+0.5*fldDo[2]*eUp-0.5*fldDo[2]*eLo); 
  atomicAdd(&fldTar[3], (dxTar[0]*(0.4330127018922193*fldDo[6]*eUpR2-0.4330127018922193*fldDo[6]*eLoR2))/dxDo[0]+(1.7320508075688772*wTar[0]*fldDo[6]*eUp-1.7320508075688772*wDo[0]*fldDo[6]*eUp-1.7320508075688772*wTar[0]*fldDo[6]*eLo+1.7320508075688772*wDo[0]*fldDo[6]*eLo)/dxDo[0]+0.5*fldDo[3]*eUp-0.5*fldDo[3]*eLo); 
  atomicAdd(&fldTar[4], (dxTar[0]*(0.4330127018922193*fldDo[8]*eUpR2-0.4330127018922193*fldDo[8]*eLoR2))/dxDo[0]+(1.7320508075688772*wTar[0]*fldDo[8]*eUp-1.7320508075688772*wDo[0]*fldDo[8]*eUp-1.7320508075688772*wTar[0]*fldDo[8]*eLo+1.7320508075688772*wDo[0]*fldDo[8]*eLo)/dxDo[0]+0.5*fldDo[4]*eUp-0.5*fldDo[4]*eLo); 
  atomicAdd(&fldTar[5], (dxTar[0]*(0.5*fldDo[5]*eUpR3-0.5*fldDo[5]*eLoR3))/dxDo[0]+(1.5*wTar[0]*fldDo[5]*eUpR2-1.5*wDo[0]*fldDo[5]*eUpR2-1.5*wTar[0]*fldDo[5]*eLoR2+1.5*wDo[0]*fldDo[5]*eLoR2)/dxDo[0]+0.4330127018922193*fldDo[2]*eUpR2-0.4330127018922193*fldDo[2]*eLoR2); 
  atomicAdd(&fldTar[6], (dxTar[0]*(0.5*fldDo[6]*eUpR3-0.5*fldDo[6]*eLoR3))/dxDo[0]+(1.5*wTar[0]*fldDo[6]*eUpR2-1.5*wDo[0]*fldDo[6]*eUpR2-1.5*wTar[0]*fldDo[6]*eLoR2+1.5*wDo[0]*fldDo[6]*eLoR2)/dxDo[0]+0.4330127018922193*fldDo[3]*eUpR2-0.4330127018922193*fldDo[3]*eLoR2); 
  atomicAdd(&fldTar[7], (dxTar[0]*(0.4330127018922193*fldDo[11]*eUpR2-0.4330127018922193*fldDo[11]*eLoR2))/dxDo[0]+(1.7320508075688772*wTar[0]*fldDo[11]*eUp-1.7320508075688772*wDo[0]*fldDo[11]*eUp-1.7320508075688772*wTar[0]*fldDo[11]*eLo+1.7320508075688772*wDo[0]*fldDo[11]*eLo)/dxDo[0]+0.5*fldDo[7]*eUp-0.5*fldDo[7]*eLo); 
  atomicAdd(&fldTar[8], (dxTar[0]*(0.5*fldDo[8]*eUpR3-0.5*fldDo[8]*eLoR3))/dxDo[0]+(1.5*wTar[0]*fldDo[8]*eUpR2-1.5*wDo[0]*fldDo[8]*eUpR2-1.5*wTar[0]*fldDo[8]*eLoR2+1.5*wDo[0]*fldDo[8]*eLoR2)/dxDo[0]+0.4330127018922193*fldDo[4]*eUpR2-0.4330127018922193*fldDo[4]*eLoR2); 
  atomicAdd(&fldTar[9], (dxTar[0]*(0.4330127018922193*fldDo[12]*eUpR2-0.4330127018922193*fldDo[12]*eLoR2))/dxDo[0]+(1.7320508075688772*wTar[0]*fldDo[12]*eUp-1.7320508075688772*wDo[0]*fldDo[12]*eUp-1.7320508075688772*wTar[0]*fldDo[12]*eLo+1.7320508075688772*wDo[0]*fldDo[12]*eLo)/dxDo[0]+0.5*fldDo[9]*eUp-0.5*fldDo[9]*eLo); 
  atomicAdd(&fldTar[10], (dxTar[0]*(0.4330127018922193*fldDo[13]*eUpR2-0.4330127018922193*fldDo[13]*eLoR2))/dxDo[0]+(1.7320508075688772*wTar[0]*fldDo[13]*eUp-1.7320508075688772*wDo[0]*fldDo[13]*eUp-1.7320508075688772*wTar[0]*fldDo[13]*eLo+1.7320508075688772*wDo[0]*fldDo[13]*eLo)/dxDo[0]+0.5*fldDo[10]*eUp-0.5*fldDo[10]*eLo); 
  atomicAdd(&fldTar[11], (dxTar[0]*(0.5*fldDo[11]*eUpR3-0.5*fldDo[11]*eLoR3))/dxDo[0]+(1.5*wTar[0]*fldDo[11]*eUpR2-1.5*wDo[0]*fldDo[11]*eUpR2-1.5*wTar[0]*fldDo[11]*eLoR2+1.5*wDo[0]*fldDo[11]*eLoR2)/dxDo[0]+0.4330127018922193*fldDo[7]*eUpR2-0.4330127018922193*fldDo[7]*eLoR2); 
  atomicAdd(&fldTar[12], (dxTar[0]*(0.5*fldDo[12]*eUpR3-0.5*fldDo[12]*eLoR3))/dxDo[0]+(1.5*wTar[0]*fldDo[12]*eUpR2-1.5*wDo[0]*fldDo[12]*eUpR2-1.5*wTar[0]*fldDo[12]*eLoR2+1.5*wDo[0]*fldDo[12]*eLoR2)/dxDo[0]+0.4330127018922193*fldDo[9]*eUpR2-0.4330127018922193*fldDo[9]*eLoR2); 
  atomicAdd(&fldTar[13], (dxTar[0]*(0.5*fldDo[13]*eUpR3-0.5*fldDo[13]*eLoR3))/dxDo[0]+(1.5*wTar[0]*fldDo[13]*eUpR2-1.5*wDo[0]*fldDo[13]*eUpR2-1.5*wTar[0]*fldDo[13]*eLoR2+1.5*wDo[0]*fldDo[13]*eLoR2)/dxDo[0]+0.4330127018922193*fldDo[10]*eUpR2-0.4330127018922193*fldDo[10]*eLoR2); 
  atomicAdd(&fldTar[14], (dxTar[0]*(0.4330127018922193*fldDo[15]*eUpR2-0.4330127018922193*fldDo[15]*eLoR2))/dxDo[0]+(1.7320508075688772*wTar[0]*fldDo[15]*eUp-1.7320508075688772*wDo[0]*fldDo[15]*eUp-1.7320508075688772*wTar[0]*fldDo[15]*eLo+1.7320508075688772*wDo[0]*fldDo[15]*eLo)/dxDo[0]+0.5*fldDo[14]*eUp-0.5*fldDo[14]*eLo); 
  atomicAdd(&fldTar[15], (dxTar[0]*(0.5*fldDo[15]*eUpR3-0.5*fldDo[15]*eLoR3))/dxDo[0]+(1.5*wTar[0]*fldDo[15]*eUpR2-1.5*wDo[0]*fldDo[15]*eUpR2-1.5*wTar[0]*fldDo[15]*eLoR2+1.5*wDo[0]*fldDo[15]*eLoR2)/dxDo[0]+0.4330127018922193*fldDo[14]*eUpR2-0.4330127018922193*fldDo[14]*eLoR2); 
#else
  const double eLoR2 = pow(eLo,2);
  const double eLoR3 = pow(eLo,3);
  const double eUpR2 = pow(eUp,2);
  const double eUpR3 = pow(eUp,3);

  fldTar[0] += (dxTar[0]*(0.4330127018922193*fldDo[1]*eUpR2-0.4330127018922193*fldDo[1]*eLoR2))/dxDo[0]+(1.7320508075688772*wTar[0]*fldDo[1]*eUp-1.7320508075688772*wDo[0]*fldDo[1]*eUp-1.7320508075688772*wTar[0]*fldDo[1]*eLo+1.7320508075688772*wDo[0]*fldDo[1]*eLo)/dxDo[0]+0.5*fldDo[0]*eUp-0.5*fldDo[0]*eLo; 
  fldTar[1] += (dxTar[0]*(0.5*fldDo[1]*eUpR3-0.5*fldDo[1]*eLoR3))/dxDo[0]+(1.5*wTar[0]*fldDo[1]*eUpR2-1.5*wDo[0]*fldDo[1]*eUpR2-1.5*wTar[0]*fldDo[1]*eLoR2+1.5*wDo[0]*fldDo[1]*eLoR2)/dxDo[0]+0.4330127018922193*fldDo[0]*eUpR2-0.4330127018922193*fldDo[0]*eLoR2; 
  fldTar[2] += (dxTar[0]*(0.4330127018922193*fldDo[5]*eUpR2-0.4330127018922193*fldDo[5]*eLoR2))/dxDo[0]+(1.7320508075688772*wTar[0]*fldDo[5]*eUp-1.7320508075688772*wDo[0]*fldDo[5]*eUp-1.7320508075688772*wTar[0]*fldDo[5]*eLo+1.7320508075688772*wDo[0]*fldDo[5]*eLo)/dxDo[0]+0.5*fldDo[2]*eUp-0.5*fldDo[2]*eLo; 
  fldTar[3] += (dxTar[0]*(0.4330127018922193*fldDo[6]*eUpR2-0.4330127018922193*fldDo[6]*eLoR2))/dxDo[0]+(1.7320508075688772*wTar[0]*fldDo[6]*eUp-1.7320508075688772*wDo[0]*fldDo[6]*eUp-1.7320508075688772*wTar[0]*fldDo[6]*eLo+1.7320508075688772*wDo[0]*fldDo[6]*eLo)/dxDo[0]+0.5*fldDo[3]*eUp-0.5*fldDo[3]*eLo; 
  fldTar[4] += (dxTar[0]*(0.4330127018922193*fldDo[8]*eUpR2-0.4330127018922193*fldDo[8]*eLoR2))/dxDo[0]+(1.7320508075688772*wTar[0]*fldDo[8]*eUp-1.7320508075688772*wDo[0]*fldDo[8]*eUp-1.7320508075688772*wTar[0]*fldDo[8]*eLo+1.7320508075688772*wDo[0]*fldDo[8]*eLo)/dxDo[0]+0.5*fldDo[4]*eUp-0.5*fldDo[4]*eLo; 
  fldTar[5] += (dxTar[0]*(0.5*fldDo[5]*eUpR3-0.5*fldDo[5]*eLoR3))/dxDo[0]+(1.5*wTar[0]*fldDo[5]*eUpR2-1.5*wDo[0]*fldDo[5]*eUpR2-1.5*wTar[0]*fldDo[5]*eLoR2+1.5*wDo[0]*fldDo[5]*eLoR2)/dxDo[0]+0.4330127018922193*fldDo[2]*eUpR2-0.4330127018922193*fldDo[2]*eLoR2; 
  fldTar[6] += (dxTar[0]*(0.5*fldDo[6]*eUpR3-0.5*fldDo[6]*eLoR3))/dxDo[0]+(1.5*wTar[0]*fldDo[6]*eUpR2-1.5*wDo[0]*fldDo[6]*eUpR2-1.5*wTar[0]*fldDo[6]*eLoR2+1.5*wDo[0]*fldDo[6]*eLoR2)/dxDo[0]+0.4330127018922193*fldDo[3]*eUpR2-0.4330127018922193*fldDo[3]*eLoR2; 
  fldTar[7] += (dxTar[0]*(0.4330127018922193*fldDo[11]*eUpR2-0.4330127018922193*fldDo[11]*eLoR2))/dxDo[0]+(1.7320508075688772*wTar[0]*fldDo[11]*eUp-1.7320508075688772*wDo[0]*fldDo[11]*eUp-1.7320508075688772*wTar[0]*fldDo[11]*eLo+1.7320508075688772*wDo[0]*fldDo[11]*eLo)/dxDo[0]+0.5*fldDo[7]*eUp-0.5*fldDo[7]*eLo; 
  fldTar[8] += (dxTar[0]*(0.5*fldDo[8]*eUpR3-0.5*fldDo[8]*eLoR3))/dxDo[0]+(1.5*wTar[0]*fldDo[8]*eUpR2-1.5*wDo[0]*fldDo[8]*eUpR2-1.5*wTar[0]*fldDo[8]*eLoR2+1.5*wDo[0]*fldDo[8]*eLoR2)/dxDo[0]+0.4330127018922193*fldDo[4]*eUpR2-0.4330127018922193*fldDo[4]*eLoR2; 
  fldTar[9] += (dxTar[0]*(0.4330127018922193*fldDo[12]*eUpR2-0.4330127018922193*fldDo[12]*eLoR2))/dxDo[0]+(1.7320508075688772*wTar[0]*fldDo[12]*eUp-1.7320508075688772*wDo[0]*fldDo[12]*eUp-1.7320508075688772*wTar[0]*fldDo[12]*eLo+1.7320508075688772*wDo[0]*fldDo[12]*eLo)/dxDo[0]+0.5*fldDo[9]*eUp-0.5*fldDo[9]*eLo; 
  fldTar[10] += (dxTar[0]*(0.4330127018922193*fldDo[13]*eUpR2-0.4330127018922193*fldDo[13]*eLoR2))/dxDo[0]+(1.7320508075688772*wTar[0]*fldDo[13]*eUp-1.7320508075688772*wDo[0]*fldDo[13]*eUp-1.7320508075688772*wTar[0]*fldDo[13]*eLo+1.7320508075688772*wDo[0]*fldDo[13]*eLo)/dxDo[0]+0.5*fldDo[10]*eUp-0.5*fldDo[10]*eLo; 
  fldTar[11] += (dxTar[0]*(0.5*fldDo[11]*eUpR3-0.5*fldDo[11]*eLoR3))/dxDo[0]+(1.5*wTar[0]*fldDo[11]*eUpR2-1.5*wDo[0]*fldDo[11]*eUpR2-1.5*wTar[0]*fldDo[11]*eLoR2+1.5*wDo[0]*fldDo[11]*eLoR2)/dxDo[0]+0.4330127018922193*fldDo[7]*eUpR2-0.4330127018922193*fldDo[7]*eLoR2; 
  fldTar[12] += (dxTar[0]*(0.5*fldDo[12]*eUpR3-0.5*fldDo[12]*eLoR3))/dxDo[0]+(1.5*wTar[0]*fldDo[12]*eUpR2-1.5*wDo[0]*fldDo[12]*eUpR2-1.5*wTar[0]*fldDo[12]*eLoR2+1.5*wDo[0]*fldDo[12]*eLoR2)/dxDo[0]+0.4330127018922193*fldDo[9]*eUpR2-0.4330127018922193*fldDo[9]*eLoR2; 
  fldTar[13] += (dxTar[0]*(0.5*fldDo[13]*eUpR3-0.5*fldDo[13]*eLoR3))/dxDo[0]+(1.5*wTar[0]*fldDo[13]*eUpR2-1.5*wDo[0]*fldDo[13]*eUpR2-1.5*wTar[0]*fldDo[13]*eLoR2+1.5*wDo[0]*fldDo[13]*eLoR2)/dxDo[0]+0.4330127018922193*fldDo[10]*eUpR2-0.4330127018922193*fldDo[10]*eLoR2; 
  fldTar[14] += (dxTar[0]*(0.4330127018922193*fldDo[15]*eUpR2-0.4330127018922193*fldDo[15]*eLoR2))/dxDo[0]+(1.7320508075688772*wTar[0]*fldDo[15]*eUp-1.7320508075688772*wDo[0]*fldDo[15]*eUp-1.7320508075688772*wTar[0]*fldDo[15]*eLo+1.7320508075688772*wDo[0]*fldDo[15]*eLo)/dxDo[0]+0.5*fldDo[14]*eUp-0.5*fldDo[14]*eLo; 
  fldTar[15] += (dxTar[0]*(0.5*fldDo[15]*eUpR3-0.5*fldDo[15]*eLoR3))/dxDo[0]+(1.5*wTar[0]*fldDo[15]*eUpR2-1.5*wDo[0]*fldDo[15]*eUpR2-1.5*wTar[0]*fldDo[15]*eLoR2+1.5*wDo[0]*fldDo[15]*eLoR2)/dxDo[0]+0.4330127018922193*fldDo[14]*eUpR2-0.4330127018922193*fldDo[14]*eLoR2; 
#endif

}

GKYL_CU_DH void dg_interpolate_vlasov_2x2v_ser_p1_y(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar) 
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
  const double eLoR2 = pow(eLo,2);
  const double eLoR3 = pow(eLo,3);
  const double eUpR2 = pow(eUp,2);
  const double eUpR3 = pow(eUp,3);

  atomicAdd(&fldTar[0], (dxTar[1]*(0.4330127018922193*fldDo[2]*eUpR2-0.4330127018922193*fldDo[2]*eLoR2))/dxDo[1]+(1.7320508075688772*wTar[1]*fldDo[2]*eUp-1.7320508075688772*wDo[1]*fldDo[2]*eUp-1.7320508075688772*wTar[1]*fldDo[2]*eLo+1.7320508075688772*wDo[1]*fldDo[2]*eLo)/dxDo[1]+0.5*fldDo[0]*eUp-0.5*fldDo[0]*eLo); 
  atomicAdd(&fldTar[1], (dxTar[1]*(0.4330127018922193*fldDo[5]*eUpR2-0.4330127018922193*fldDo[5]*eLoR2))/dxDo[1]+(1.7320508075688772*wTar[1]*fldDo[5]*eUp-1.7320508075688772*wDo[1]*fldDo[5]*eUp-1.7320508075688772*wTar[1]*fldDo[5]*eLo+1.7320508075688772*wDo[1]*fldDo[5]*eLo)/dxDo[1]+0.5*fldDo[1]*eUp-0.5*fldDo[1]*eLo); 
  atomicAdd(&fldTar[2], (dxTar[1]*(0.5*fldDo[2]*eUpR3-0.5*fldDo[2]*eLoR3))/dxDo[1]+(1.5*wTar[1]*fldDo[2]*eUpR2-1.5*wDo[1]*fldDo[2]*eUpR2-1.5*wTar[1]*fldDo[2]*eLoR2+1.5*wDo[1]*fldDo[2]*eLoR2)/dxDo[1]+0.4330127018922193*fldDo[0]*eUpR2-0.4330127018922193*fldDo[0]*eLoR2); 
  atomicAdd(&fldTar[3], (dxTar[1]*(0.4330127018922193*fldDo[7]*eUpR2-0.4330127018922193*fldDo[7]*eLoR2))/dxDo[1]+(1.7320508075688772*wTar[1]*fldDo[7]*eUp-1.7320508075688772*wDo[1]*fldDo[7]*eUp-1.7320508075688772*wTar[1]*fldDo[7]*eLo+1.7320508075688772*wDo[1]*fldDo[7]*eLo)/dxDo[1]+0.5*fldDo[3]*eUp-0.5*fldDo[3]*eLo); 
  atomicAdd(&fldTar[4], (dxTar[1]*(0.4330127018922193*fldDo[9]*eUpR2-0.4330127018922193*fldDo[9]*eLoR2))/dxDo[1]+(1.7320508075688772*wTar[1]*fldDo[9]*eUp-1.7320508075688772*wDo[1]*fldDo[9]*eUp-1.7320508075688772*wTar[1]*fldDo[9]*eLo+1.7320508075688772*wDo[1]*fldDo[9]*eLo)/dxDo[1]+0.5*fldDo[4]*eUp-0.5*fldDo[4]*eLo); 
  atomicAdd(&fldTar[5], (dxTar[1]*(0.5*fldDo[5]*eUpR3-0.5*fldDo[5]*eLoR3))/dxDo[1]+(1.5*wTar[1]*fldDo[5]*eUpR2-1.5*wDo[1]*fldDo[5]*eUpR2-1.5*wTar[1]*fldDo[5]*eLoR2+1.5*wDo[1]*fldDo[5]*eLoR2)/dxDo[1]+0.4330127018922193*fldDo[1]*eUpR2-0.4330127018922193*fldDo[1]*eLoR2); 
  atomicAdd(&fldTar[6], (dxTar[1]*(0.4330127018922193*fldDo[11]*eUpR2-0.4330127018922193*fldDo[11]*eLoR2))/dxDo[1]+(1.7320508075688772*wTar[1]*fldDo[11]*eUp-1.7320508075688772*wDo[1]*fldDo[11]*eUp-1.7320508075688772*wTar[1]*fldDo[11]*eLo+1.7320508075688772*wDo[1]*fldDo[11]*eLo)/dxDo[1]+0.5*fldDo[6]*eUp-0.5*fldDo[6]*eLo); 
  atomicAdd(&fldTar[7], (dxTar[1]*(0.5*fldDo[7]*eUpR3-0.5*fldDo[7]*eLoR3))/dxDo[1]+(1.5*wTar[1]*fldDo[7]*eUpR2-1.5*wDo[1]*fldDo[7]*eUpR2-1.5*wTar[1]*fldDo[7]*eLoR2+1.5*wDo[1]*fldDo[7]*eLoR2)/dxDo[1]+0.4330127018922193*fldDo[3]*eUpR2-0.4330127018922193*fldDo[3]*eLoR2); 
  atomicAdd(&fldTar[8], (dxTar[1]*(0.4330127018922193*fldDo[12]*eUpR2-0.4330127018922193*fldDo[12]*eLoR2))/dxDo[1]+(1.7320508075688772*wTar[1]*fldDo[12]*eUp-1.7320508075688772*wDo[1]*fldDo[12]*eUp-1.7320508075688772*wTar[1]*fldDo[12]*eLo+1.7320508075688772*wDo[1]*fldDo[12]*eLo)/dxDo[1]+0.5*fldDo[8]*eUp-0.5*fldDo[8]*eLo); 
  atomicAdd(&fldTar[9], (dxTar[1]*(0.5*fldDo[9]*eUpR3-0.5*fldDo[9]*eLoR3))/dxDo[1]+(1.5*wTar[1]*fldDo[9]*eUpR2-1.5*wDo[1]*fldDo[9]*eUpR2-1.5*wTar[1]*fldDo[9]*eLoR2+1.5*wDo[1]*fldDo[9]*eLoR2)/dxDo[1]+0.4330127018922193*fldDo[4]*eUpR2-0.4330127018922193*fldDo[4]*eLoR2); 
  atomicAdd(&fldTar[10], (dxTar[1]*(0.4330127018922193*fldDo[14]*eUpR2-0.4330127018922193*fldDo[14]*eLoR2))/dxDo[1]+(1.7320508075688772*wTar[1]*fldDo[14]*eUp-1.7320508075688772*wDo[1]*fldDo[14]*eUp-1.7320508075688772*wTar[1]*fldDo[14]*eLo+1.7320508075688772*wDo[1]*fldDo[14]*eLo)/dxDo[1]+0.5*fldDo[10]*eUp-0.5*fldDo[10]*eLo); 
  atomicAdd(&fldTar[11], (dxTar[1]*(0.5*fldDo[11]*eUpR3-0.5*fldDo[11]*eLoR3))/dxDo[1]+(1.5*wTar[1]*fldDo[11]*eUpR2-1.5*wDo[1]*fldDo[11]*eUpR2-1.5*wTar[1]*fldDo[11]*eLoR2+1.5*wDo[1]*fldDo[11]*eLoR2)/dxDo[1]+0.4330127018922193*fldDo[6]*eUpR2-0.4330127018922193*fldDo[6]*eLoR2); 
  atomicAdd(&fldTar[12], (dxTar[1]*(0.5*fldDo[12]*eUpR3-0.5*fldDo[12]*eLoR3))/dxDo[1]+(1.5*wTar[1]*fldDo[12]*eUpR2-1.5*wDo[1]*fldDo[12]*eUpR2-1.5*wTar[1]*fldDo[12]*eLoR2+1.5*wDo[1]*fldDo[12]*eLoR2)/dxDo[1]+0.4330127018922193*fldDo[8]*eUpR2-0.4330127018922193*fldDo[8]*eLoR2); 
  atomicAdd(&fldTar[13], (dxTar[1]*(0.4330127018922193*fldDo[15]*eUpR2-0.4330127018922193*fldDo[15]*eLoR2))/dxDo[1]+(1.7320508075688772*wTar[1]*fldDo[15]*eUp-1.7320508075688772*wDo[1]*fldDo[15]*eUp-1.7320508075688772*wTar[1]*fldDo[15]*eLo+1.7320508075688772*wDo[1]*fldDo[15]*eLo)/dxDo[1]+0.5*fldDo[13]*eUp-0.5*fldDo[13]*eLo); 
  atomicAdd(&fldTar[14], (dxTar[1]*(0.5*fldDo[14]*eUpR3-0.5*fldDo[14]*eLoR3))/dxDo[1]+(1.5*wTar[1]*fldDo[14]*eUpR2-1.5*wDo[1]*fldDo[14]*eUpR2-1.5*wTar[1]*fldDo[14]*eLoR2+1.5*wDo[1]*fldDo[14]*eLoR2)/dxDo[1]+0.4330127018922193*fldDo[10]*eUpR2-0.4330127018922193*fldDo[10]*eLoR2); 
  atomicAdd(&fldTar[15], (dxTar[1]*(0.5*fldDo[15]*eUpR3-0.5*fldDo[15]*eLoR3))/dxDo[1]+(1.5*wTar[1]*fldDo[15]*eUpR2-1.5*wDo[1]*fldDo[15]*eUpR2-1.5*wTar[1]*fldDo[15]*eLoR2+1.5*wDo[1]*fldDo[15]*eLoR2)/dxDo[1]+0.4330127018922193*fldDo[13]*eUpR2-0.4330127018922193*fldDo[13]*eLoR2); 
#else
  const double eLoR2 = pow(eLo,2);
  const double eLoR3 = pow(eLo,3);
  const double eUpR2 = pow(eUp,2);
  const double eUpR3 = pow(eUp,3);

  fldTar[0] += (dxTar[1]*(0.4330127018922193*fldDo[2]*eUpR2-0.4330127018922193*fldDo[2]*eLoR2))/dxDo[1]+(1.7320508075688772*wTar[1]*fldDo[2]*eUp-1.7320508075688772*wDo[1]*fldDo[2]*eUp-1.7320508075688772*wTar[1]*fldDo[2]*eLo+1.7320508075688772*wDo[1]*fldDo[2]*eLo)/dxDo[1]+0.5*fldDo[0]*eUp-0.5*fldDo[0]*eLo; 
  fldTar[1] += (dxTar[1]*(0.4330127018922193*fldDo[5]*eUpR2-0.4330127018922193*fldDo[5]*eLoR2))/dxDo[1]+(1.7320508075688772*wTar[1]*fldDo[5]*eUp-1.7320508075688772*wDo[1]*fldDo[5]*eUp-1.7320508075688772*wTar[1]*fldDo[5]*eLo+1.7320508075688772*wDo[1]*fldDo[5]*eLo)/dxDo[1]+0.5*fldDo[1]*eUp-0.5*fldDo[1]*eLo; 
  fldTar[2] += (dxTar[1]*(0.5*fldDo[2]*eUpR3-0.5*fldDo[2]*eLoR3))/dxDo[1]+(1.5*wTar[1]*fldDo[2]*eUpR2-1.5*wDo[1]*fldDo[2]*eUpR2-1.5*wTar[1]*fldDo[2]*eLoR2+1.5*wDo[1]*fldDo[2]*eLoR2)/dxDo[1]+0.4330127018922193*fldDo[0]*eUpR2-0.4330127018922193*fldDo[0]*eLoR2; 
  fldTar[3] += (dxTar[1]*(0.4330127018922193*fldDo[7]*eUpR2-0.4330127018922193*fldDo[7]*eLoR2))/dxDo[1]+(1.7320508075688772*wTar[1]*fldDo[7]*eUp-1.7320508075688772*wDo[1]*fldDo[7]*eUp-1.7320508075688772*wTar[1]*fldDo[7]*eLo+1.7320508075688772*wDo[1]*fldDo[7]*eLo)/dxDo[1]+0.5*fldDo[3]*eUp-0.5*fldDo[3]*eLo; 
  fldTar[4] += (dxTar[1]*(0.4330127018922193*fldDo[9]*eUpR2-0.4330127018922193*fldDo[9]*eLoR2))/dxDo[1]+(1.7320508075688772*wTar[1]*fldDo[9]*eUp-1.7320508075688772*wDo[1]*fldDo[9]*eUp-1.7320508075688772*wTar[1]*fldDo[9]*eLo+1.7320508075688772*wDo[1]*fldDo[9]*eLo)/dxDo[1]+0.5*fldDo[4]*eUp-0.5*fldDo[4]*eLo; 
  fldTar[5] += (dxTar[1]*(0.5*fldDo[5]*eUpR3-0.5*fldDo[5]*eLoR3))/dxDo[1]+(1.5*wTar[1]*fldDo[5]*eUpR2-1.5*wDo[1]*fldDo[5]*eUpR2-1.5*wTar[1]*fldDo[5]*eLoR2+1.5*wDo[1]*fldDo[5]*eLoR2)/dxDo[1]+0.4330127018922193*fldDo[1]*eUpR2-0.4330127018922193*fldDo[1]*eLoR2; 
  fldTar[6] += (dxTar[1]*(0.4330127018922193*fldDo[11]*eUpR2-0.4330127018922193*fldDo[11]*eLoR2))/dxDo[1]+(1.7320508075688772*wTar[1]*fldDo[11]*eUp-1.7320508075688772*wDo[1]*fldDo[11]*eUp-1.7320508075688772*wTar[1]*fldDo[11]*eLo+1.7320508075688772*wDo[1]*fldDo[11]*eLo)/dxDo[1]+0.5*fldDo[6]*eUp-0.5*fldDo[6]*eLo; 
  fldTar[7] += (dxTar[1]*(0.5*fldDo[7]*eUpR3-0.5*fldDo[7]*eLoR3))/dxDo[1]+(1.5*wTar[1]*fldDo[7]*eUpR2-1.5*wDo[1]*fldDo[7]*eUpR2-1.5*wTar[1]*fldDo[7]*eLoR2+1.5*wDo[1]*fldDo[7]*eLoR2)/dxDo[1]+0.4330127018922193*fldDo[3]*eUpR2-0.4330127018922193*fldDo[3]*eLoR2; 
  fldTar[8] += (dxTar[1]*(0.4330127018922193*fldDo[12]*eUpR2-0.4330127018922193*fldDo[12]*eLoR2))/dxDo[1]+(1.7320508075688772*wTar[1]*fldDo[12]*eUp-1.7320508075688772*wDo[1]*fldDo[12]*eUp-1.7320508075688772*wTar[1]*fldDo[12]*eLo+1.7320508075688772*wDo[1]*fldDo[12]*eLo)/dxDo[1]+0.5*fldDo[8]*eUp-0.5*fldDo[8]*eLo; 
  fldTar[9] += (dxTar[1]*(0.5*fldDo[9]*eUpR3-0.5*fldDo[9]*eLoR3))/dxDo[1]+(1.5*wTar[1]*fldDo[9]*eUpR2-1.5*wDo[1]*fldDo[9]*eUpR2-1.5*wTar[1]*fldDo[9]*eLoR2+1.5*wDo[1]*fldDo[9]*eLoR2)/dxDo[1]+0.4330127018922193*fldDo[4]*eUpR2-0.4330127018922193*fldDo[4]*eLoR2; 
  fldTar[10] += (dxTar[1]*(0.4330127018922193*fldDo[14]*eUpR2-0.4330127018922193*fldDo[14]*eLoR2))/dxDo[1]+(1.7320508075688772*wTar[1]*fldDo[14]*eUp-1.7320508075688772*wDo[1]*fldDo[14]*eUp-1.7320508075688772*wTar[1]*fldDo[14]*eLo+1.7320508075688772*wDo[1]*fldDo[14]*eLo)/dxDo[1]+0.5*fldDo[10]*eUp-0.5*fldDo[10]*eLo; 
  fldTar[11] += (dxTar[1]*(0.5*fldDo[11]*eUpR3-0.5*fldDo[11]*eLoR3))/dxDo[1]+(1.5*wTar[1]*fldDo[11]*eUpR2-1.5*wDo[1]*fldDo[11]*eUpR2-1.5*wTar[1]*fldDo[11]*eLoR2+1.5*wDo[1]*fldDo[11]*eLoR2)/dxDo[1]+0.4330127018922193*fldDo[6]*eUpR2-0.4330127018922193*fldDo[6]*eLoR2; 
  fldTar[12] += (dxTar[1]*(0.5*fldDo[12]*eUpR3-0.5*fldDo[12]*eLoR3))/dxDo[1]+(1.5*wTar[1]*fldDo[12]*eUpR2-1.5*wDo[1]*fldDo[12]*eUpR2-1.5*wTar[1]*fldDo[12]*eLoR2+1.5*wDo[1]*fldDo[12]*eLoR2)/dxDo[1]+0.4330127018922193*fldDo[8]*eUpR2-0.4330127018922193*fldDo[8]*eLoR2; 
  fldTar[13] += (dxTar[1]*(0.4330127018922193*fldDo[15]*eUpR2-0.4330127018922193*fldDo[15]*eLoR2))/dxDo[1]+(1.7320508075688772*wTar[1]*fldDo[15]*eUp-1.7320508075688772*wDo[1]*fldDo[15]*eUp-1.7320508075688772*wTar[1]*fldDo[15]*eLo+1.7320508075688772*wDo[1]*fldDo[15]*eLo)/dxDo[1]+0.5*fldDo[13]*eUp-0.5*fldDo[13]*eLo; 
  fldTar[14] += (dxTar[1]*(0.5*fldDo[14]*eUpR3-0.5*fldDo[14]*eLoR3))/dxDo[1]+(1.5*wTar[1]*fldDo[14]*eUpR2-1.5*wDo[1]*fldDo[14]*eUpR2-1.5*wTar[1]*fldDo[14]*eLoR2+1.5*wDo[1]*fldDo[14]*eLoR2)/dxDo[1]+0.4330127018922193*fldDo[10]*eUpR2-0.4330127018922193*fldDo[10]*eLoR2; 
  fldTar[15] += (dxTar[1]*(0.5*fldDo[15]*eUpR3-0.5*fldDo[15]*eLoR3))/dxDo[1]+(1.5*wTar[1]*fldDo[15]*eUpR2-1.5*wDo[1]*fldDo[15]*eUpR2-1.5*wTar[1]*fldDo[15]*eLoR2+1.5*wDo[1]*fldDo[15]*eLoR2)/dxDo[1]+0.4330127018922193*fldDo[13]*eUpR2-0.4330127018922193*fldDo[13]*eLoR2; 
#endif

}

GKYL_CU_DH void dg_interpolate_vlasov_2x2v_ser_p1_vx(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar) 
{ 
  // wDo: cell center of donor cell.
  // wTar: cell center of target cell.
  // dxDo: cell length of donor cell.
  // dxTar: cell length of target cell.
  // fldDo: donor field.
  // fldTar: target field in cells pointed to by the stencil.

  double eLo = fmax(-1.0,1.0-(2.0*(wTar[2]-1.0*wDo[2]+0.5*dxTar[2]+0.5*dxDo[2]))/dxTar[2]);
  double eUp = fmin( 1.0,(2.0*(-(1.0*wTar[2])+wDo[2]+0.5*dxTar[2]+0.5*dxDo[2]))/dxTar[2]-1.0);

#ifdef __CUDA_ARCH__
  const double eLoR2 = pow(eLo,2);
  const double eLoR3 = pow(eLo,3);
  const double eUpR2 = pow(eUp,2);
  const double eUpR3 = pow(eUp,3);

  atomicAdd(&fldTar[0], (dxTar[2]*(0.4330127018922193*fldDo[3]*eUpR2-0.4330127018922193*fldDo[3]*eLoR2))/dxDo[2]+(1.7320508075688772*wTar[2]*fldDo[3]*eUp-1.7320508075688772*wDo[2]*fldDo[3]*eUp-1.7320508075688772*wTar[2]*fldDo[3]*eLo+1.7320508075688772*wDo[2]*fldDo[3]*eLo)/dxDo[2]+0.5*fldDo[0]*eUp-0.5*fldDo[0]*eLo); 
  atomicAdd(&fldTar[1], (dxTar[2]*(0.4330127018922193*fldDo[6]*eUpR2-0.4330127018922193*fldDo[6]*eLoR2))/dxDo[2]+(1.7320508075688772*wTar[2]*fldDo[6]*eUp-1.7320508075688772*wDo[2]*fldDo[6]*eUp-1.7320508075688772*wTar[2]*fldDo[6]*eLo+1.7320508075688772*wDo[2]*fldDo[6]*eLo)/dxDo[2]+0.5*fldDo[1]*eUp-0.5*fldDo[1]*eLo); 
  atomicAdd(&fldTar[2], (dxTar[2]*(0.4330127018922193*fldDo[7]*eUpR2-0.4330127018922193*fldDo[7]*eLoR2))/dxDo[2]+(1.7320508075688772*wTar[2]*fldDo[7]*eUp-1.7320508075688772*wDo[2]*fldDo[7]*eUp-1.7320508075688772*wTar[2]*fldDo[7]*eLo+1.7320508075688772*wDo[2]*fldDo[7]*eLo)/dxDo[2]+0.5*fldDo[2]*eUp-0.5*fldDo[2]*eLo); 
  atomicAdd(&fldTar[3], (dxTar[2]*(0.5*fldDo[3]*eUpR3-0.5*fldDo[3]*eLoR3))/dxDo[2]+(1.5*wTar[2]*fldDo[3]*eUpR2-1.5*wDo[2]*fldDo[3]*eUpR2-1.5*wTar[2]*fldDo[3]*eLoR2+1.5*wDo[2]*fldDo[3]*eLoR2)/dxDo[2]+0.4330127018922193*fldDo[0]*eUpR2-0.4330127018922193*fldDo[0]*eLoR2); 
  atomicAdd(&fldTar[4], (dxTar[2]*(0.4330127018922193*fldDo[10]*eUpR2-0.4330127018922193*fldDo[10]*eLoR2))/dxDo[2]+(1.7320508075688772*wTar[2]*fldDo[10]*eUp-1.7320508075688772*wDo[2]*fldDo[10]*eUp-1.7320508075688772*wTar[2]*fldDo[10]*eLo+1.7320508075688772*wDo[2]*fldDo[10]*eLo)/dxDo[2]+0.5*fldDo[4]*eUp-0.5*fldDo[4]*eLo); 
  atomicAdd(&fldTar[5], (dxTar[2]*(0.4330127018922193*fldDo[11]*eUpR2-0.4330127018922193*fldDo[11]*eLoR2))/dxDo[2]+(1.7320508075688772*wTar[2]*fldDo[11]*eUp-1.7320508075688772*wDo[2]*fldDo[11]*eUp-1.7320508075688772*wTar[2]*fldDo[11]*eLo+1.7320508075688772*wDo[2]*fldDo[11]*eLo)/dxDo[2]+0.5*fldDo[5]*eUp-0.5*fldDo[5]*eLo); 
  atomicAdd(&fldTar[6], (dxTar[2]*(0.5*fldDo[6]*eUpR3-0.5*fldDo[6]*eLoR3))/dxDo[2]+(1.5*wTar[2]*fldDo[6]*eUpR2-1.5*wDo[2]*fldDo[6]*eUpR2-1.5*wTar[2]*fldDo[6]*eLoR2+1.5*wDo[2]*fldDo[6]*eLoR2)/dxDo[2]+0.4330127018922193*fldDo[1]*eUpR2-0.4330127018922193*fldDo[1]*eLoR2); 
  atomicAdd(&fldTar[7], (dxTar[2]*(0.5*fldDo[7]*eUpR3-0.5*fldDo[7]*eLoR3))/dxDo[2]+(1.5*wTar[2]*fldDo[7]*eUpR2-1.5*wDo[2]*fldDo[7]*eUpR2-1.5*wTar[2]*fldDo[7]*eLoR2+1.5*wDo[2]*fldDo[7]*eLoR2)/dxDo[2]+0.4330127018922193*fldDo[2]*eUpR2-0.4330127018922193*fldDo[2]*eLoR2); 
  atomicAdd(&fldTar[8], (dxTar[2]*(0.4330127018922193*fldDo[13]*eUpR2-0.4330127018922193*fldDo[13]*eLoR2))/dxDo[2]+(1.7320508075688772*wTar[2]*fldDo[13]*eUp-1.7320508075688772*wDo[2]*fldDo[13]*eUp-1.7320508075688772*wTar[2]*fldDo[13]*eLo+1.7320508075688772*wDo[2]*fldDo[13]*eLo)/dxDo[2]+0.5*fldDo[8]*eUp-0.5*fldDo[8]*eLo); 
  atomicAdd(&fldTar[9], (dxTar[2]*(0.4330127018922193*fldDo[14]*eUpR2-0.4330127018922193*fldDo[14]*eLoR2))/dxDo[2]+(1.7320508075688772*wTar[2]*fldDo[14]*eUp-1.7320508075688772*wDo[2]*fldDo[14]*eUp-1.7320508075688772*wTar[2]*fldDo[14]*eLo+1.7320508075688772*wDo[2]*fldDo[14]*eLo)/dxDo[2]+0.5*fldDo[9]*eUp-0.5*fldDo[9]*eLo); 
  atomicAdd(&fldTar[10], (dxTar[2]*(0.5*fldDo[10]*eUpR3-0.5*fldDo[10]*eLoR3))/dxDo[2]+(1.5*wTar[2]*fldDo[10]*eUpR2-1.5*wDo[2]*fldDo[10]*eUpR2-1.5*wTar[2]*fldDo[10]*eLoR2+1.5*wDo[2]*fldDo[10]*eLoR2)/dxDo[2]+0.4330127018922193*fldDo[4]*eUpR2-0.4330127018922193*fldDo[4]*eLoR2); 
  atomicAdd(&fldTar[11], (dxTar[2]*(0.5*fldDo[11]*eUpR3-0.5*fldDo[11]*eLoR3))/dxDo[2]+(1.5*wTar[2]*fldDo[11]*eUpR2-1.5*wDo[2]*fldDo[11]*eUpR2-1.5*wTar[2]*fldDo[11]*eLoR2+1.5*wDo[2]*fldDo[11]*eLoR2)/dxDo[2]+0.4330127018922193*fldDo[5]*eUpR2-0.4330127018922193*fldDo[5]*eLoR2); 
  atomicAdd(&fldTar[12], (dxTar[2]*(0.4330127018922193*fldDo[15]*eUpR2-0.4330127018922193*fldDo[15]*eLoR2))/dxDo[2]+(1.7320508075688772*wTar[2]*fldDo[15]*eUp-1.7320508075688772*wDo[2]*fldDo[15]*eUp-1.7320508075688772*wTar[2]*fldDo[15]*eLo+1.7320508075688772*wDo[2]*fldDo[15]*eLo)/dxDo[2]+0.5*fldDo[12]*eUp-0.5*fldDo[12]*eLo); 
  atomicAdd(&fldTar[13], (dxTar[2]*(0.5*fldDo[13]*eUpR3-0.5*fldDo[13]*eLoR3))/dxDo[2]+(1.5*wTar[2]*fldDo[13]*eUpR2-1.5*wDo[2]*fldDo[13]*eUpR2-1.5*wTar[2]*fldDo[13]*eLoR2+1.5*wDo[2]*fldDo[13]*eLoR2)/dxDo[2]+0.4330127018922193*fldDo[8]*eUpR2-0.4330127018922193*fldDo[8]*eLoR2); 
  atomicAdd(&fldTar[14], (dxTar[2]*(0.5*fldDo[14]*eUpR3-0.5*fldDo[14]*eLoR3))/dxDo[2]+(1.5*wTar[2]*fldDo[14]*eUpR2-1.5*wDo[2]*fldDo[14]*eUpR2-1.5*wTar[2]*fldDo[14]*eLoR2+1.5*wDo[2]*fldDo[14]*eLoR2)/dxDo[2]+0.4330127018922193*fldDo[9]*eUpR2-0.4330127018922193*fldDo[9]*eLoR2); 
  atomicAdd(&fldTar[15], (dxTar[2]*(0.5*fldDo[15]*eUpR3-0.5*fldDo[15]*eLoR3))/dxDo[2]+(1.5*wTar[2]*fldDo[15]*eUpR2-1.5*wDo[2]*fldDo[15]*eUpR2-1.5*wTar[2]*fldDo[15]*eLoR2+1.5*wDo[2]*fldDo[15]*eLoR2)/dxDo[2]+0.4330127018922193*fldDo[12]*eUpR2-0.4330127018922193*fldDo[12]*eLoR2); 
#else
  const double eLoR2 = pow(eLo,2);
  const double eLoR3 = pow(eLo,3);
  const double eUpR2 = pow(eUp,2);
  const double eUpR3 = pow(eUp,3);

  fldTar[0] += (dxTar[2]*(0.4330127018922193*fldDo[3]*eUpR2-0.4330127018922193*fldDo[3]*eLoR2))/dxDo[2]+(1.7320508075688772*wTar[2]*fldDo[3]*eUp-1.7320508075688772*wDo[2]*fldDo[3]*eUp-1.7320508075688772*wTar[2]*fldDo[3]*eLo+1.7320508075688772*wDo[2]*fldDo[3]*eLo)/dxDo[2]+0.5*fldDo[0]*eUp-0.5*fldDo[0]*eLo; 
  fldTar[1] += (dxTar[2]*(0.4330127018922193*fldDo[6]*eUpR2-0.4330127018922193*fldDo[6]*eLoR2))/dxDo[2]+(1.7320508075688772*wTar[2]*fldDo[6]*eUp-1.7320508075688772*wDo[2]*fldDo[6]*eUp-1.7320508075688772*wTar[2]*fldDo[6]*eLo+1.7320508075688772*wDo[2]*fldDo[6]*eLo)/dxDo[2]+0.5*fldDo[1]*eUp-0.5*fldDo[1]*eLo; 
  fldTar[2] += (dxTar[2]*(0.4330127018922193*fldDo[7]*eUpR2-0.4330127018922193*fldDo[7]*eLoR2))/dxDo[2]+(1.7320508075688772*wTar[2]*fldDo[7]*eUp-1.7320508075688772*wDo[2]*fldDo[7]*eUp-1.7320508075688772*wTar[2]*fldDo[7]*eLo+1.7320508075688772*wDo[2]*fldDo[7]*eLo)/dxDo[2]+0.5*fldDo[2]*eUp-0.5*fldDo[2]*eLo; 
  fldTar[3] += (dxTar[2]*(0.5*fldDo[3]*eUpR3-0.5*fldDo[3]*eLoR3))/dxDo[2]+(1.5*wTar[2]*fldDo[3]*eUpR2-1.5*wDo[2]*fldDo[3]*eUpR2-1.5*wTar[2]*fldDo[3]*eLoR2+1.5*wDo[2]*fldDo[3]*eLoR2)/dxDo[2]+0.4330127018922193*fldDo[0]*eUpR2-0.4330127018922193*fldDo[0]*eLoR2; 
  fldTar[4] += (dxTar[2]*(0.4330127018922193*fldDo[10]*eUpR2-0.4330127018922193*fldDo[10]*eLoR2))/dxDo[2]+(1.7320508075688772*wTar[2]*fldDo[10]*eUp-1.7320508075688772*wDo[2]*fldDo[10]*eUp-1.7320508075688772*wTar[2]*fldDo[10]*eLo+1.7320508075688772*wDo[2]*fldDo[10]*eLo)/dxDo[2]+0.5*fldDo[4]*eUp-0.5*fldDo[4]*eLo; 
  fldTar[5] += (dxTar[2]*(0.4330127018922193*fldDo[11]*eUpR2-0.4330127018922193*fldDo[11]*eLoR2))/dxDo[2]+(1.7320508075688772*wTar[2]*fldDo[11]*eUp-1.7320508075688772*wDo[2]*fldDo[11]*eUp-1.7320508075688772*wTar[2]*fldDo[11]*eLo+1.7320508075688772*wDo[2]*fldDo[11]*eLo)/dxDo[2]+0.5*fldDo[5]*eUp-0.5*fldDo[5]*eLo; 
  fldTar[6] += (dxTar[2]*(0.5*fldDo[6]*eUpR3-0.5*fldDo[6]*eLoR3))/dxDo[2]+(1.5*wTar[2]*fldDo[6]*eUpR2-1.5*wDo[2]*fldDo[6]*eUpR2-1.5*wTar[2]*fldDo[6]*eLoR2+1.5*wDo[2]*fldDo[6]*eLoR2)/dxDo[2]+0.4330127018922193*fldDo[1]*eUpR2-0.4330127018922193*fldDo[1]*eLoR2; 
  fldTar[7] += (dxTar[2]*(0.5*fldDo[7]*eUpR3-0.5*fldDo[7]*eLoR3))/dxDo[2]+(1.5*wTar[2]*fldDo[7]*eUpR2-1.5*wDo[2]*fldDo[7]*eUpR2-1.5*wTar[2]*fldDo[7]*eLoR2+1.5*wDo[2]*fldDo[7]*eLoR2)/dxDo[2]+0.4330127018922193*fldDo[2]*eUpR2-0.4330127018922193*fldDo[2]*eLoR2; 
  fldTar[8] += (dxTar[2]*(0.4330127018922193*fldDo[13]*eUpR2-0.4330127018922193*fldDo[13]*eLoR2))/dxDo[2]+(1.7320508075688772*wTar[2]*fldDo[13]*eUp-1.7320508075688772*wDo[2]*fldDo[13]*eUp-1.7320508075688772*wTar[2]*fldDo[13]*eLo+1.7320508075688772*wDo[2]*fldDo[13]*eLo)/dxDo[2]+0.5*fldDo[8]*eUp-0.5*fldDo[8]*eLo; 
  fldTar[9] += (dxTar[2]*(0.4330127018922193*fldDo[14]*eUpR2-0.4330127018922193*fldDo[14]*eLoR2))/dxDo[2]+(1.7320508075688772*wTar[2]*fldDo[14]*eUp-1.7320508075688772*wDo[2]*fldDo[14]*eUp-1.7320508075688772*wTar[2]*fldDo[14]*eLo+1.7320508075688772*wDo[2]*fldDo[14]*eLo)/dxDo[2]+0.5*fldDo[9]*eUp-0.5*fldDo[9]*eLo; 
  fldTar[10] += (dxTar[2]*(0.5*fldDo[10]*eUpR3-0.5*fldDo[10]*eLoR3))/dxDo[2]+(1.5*wTar[2]*fldDo[10]*eUpR2-1.5*wDo[2]*fldDo[10]*eUpR2-1.5*wTar[2]*fldDo[10]*eLoR2+1.5*wDo[2]*fldDo[10]*eLoR2)/dxDo[2]+0.4330127018922193*fldDo[4]*eUpR2-0.4330127018922193*fldDo[4]*eLoR2; 
  fldTar[11] += (dxTar[2]*(0.5*fldDo[11]*eUpR3-0.5*fldDo[11]*eLoR3))/dxDo[2]+(1.5*wTar[2]*fldDo[11]*eUpR2-1.5*wDo[2]*fldDo[11]*eUpR2-1.5*wTar[2]*fldDo[11]*eLoR2+1.5*wDo[2]*fldDo[11]*eLoR2)/dxDo[2]+0.4330127018922193*fldDo[5]*eUpR2-0.4330127018922193*fldDo[5]*eLoR2; 
  fldTar[12] += (dxTar[2]*(0.4330127018922193*fldDo[15]*eUpR2-0.4330127018922193*fldDo[15]*eLoR2))/dxDo[2]+(1.7320508075688772*wTar[2]*fldDo[15]*eUp-1.7320508075688772*wDo[2]*fldDo[15]*eUp-1.7320508075688772*wTar[2]*fldDo[15]*eLo+1.7320508075688772*wDo[2]*fldDo[15]*eLo)/dxDo[2]+0.5*fldDo[12]*eUp-0.5*fldDo[12]*eLo; 
  fldTar[13] += (dxTar[2]*(0.5*fldDo[13]*eUpR3-0.5*fldDo[13]*eLoR3))/dxDo[2]+(1.5*wTar[2]*fldDo[13]*eUpR2-1.5*wDo[2]*fldDo[13]*eUpR2-1.5*wTar[2]*fldDo[13]*eLoR2+1.5*wDo[2]*fldDo[13]*eLoR2)/dxDo[2]+0.4330127018922193*fldDo[8]*eUpR2-0.4330127018922193*fldDo[8]*eLoR2; 
  fldTar[14] += (dxTar[2]*(0.5*fldDo[14]*eUpR3-0.5*fldDo[14]*eLoR3))/dxDo[2]+(1.5*wTar[2]*fldDo[14]*eUpR2-1.5*wDo[2]*fldDo[14]*eUpR2-1.5*wTar[2]*fldDo[14]*eLoR2+1.5*wDo[2]*fldDo[14]*eLoR2)/dxDo[2]+0.4330127018922193*fldDo[9]*eUpR2-0.4330127018922193*fldDo[9]*eLoR2; 
  fldTar[15] += (dxTar[2]*(0.5*fldDo[15]*eUpR3-0.5*fldDo[15]*eLoR3))/dxDo[2]+(1.5*wTar[2]*fldDo[15]*eUpR2-1.5*wDo[2]*fldDo[15]*eUpR2-1.5*wTar[2]*fldDo[15]*eLoR2+1.5*wDo[2]*fldDo[15]*eLoR2)/dxDo[2]+0.4330127018922193*fldDo[12]*eUpR2-0.4330127018922193*fldDo[12]*eLoR2; 
#endif

}

GKYL_CU_DH void dg_interpolate_vlasov_2x2v_ser_p1_vy(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar) 
{ 
  // wDo: cell center of donor cell.
  // wTar: cell center of target cell.
  // dxDo: cell length of donor cell.
  // dxTar: cell length of target cell.
  // fldDo: donor field.
  // fldTar: target field in cells pointed to by the stencil.

  double eLo = fmax(-1.0,1.0-(2.0*(wTar[3]-1.0*wDo[3]+0.5*dxTar[3]+0.5*dxDo[3]))/dxTar[3]);
  double eUp = fmin( 1.0,(2.0*(-(1.0*wTar[3])+wDo[3]+0.5*dxTar[3]+0.5*dxDo[3]))/dxTar[3]-1.0);

#ifdef __CUDA_ARCH__
  const double eLoR2 = pow(eLo,2);
  const double eLoR3 = pow(eLo,3);
  const double eUpR2 = pow(eUp,2);
  const double eUpR3 = pow(eUp,3);

  atomicAdd(&fldTar[0], (dxTar[3]*(0.4330127018922193*fldDo[4]*eUpR2-0.4330127018922193*fldDo[4]*eLoR2))/dxDo[3]+(1.7320508075688772*wTar[3]*fldDo[4]*eUp-1.7320508075688772*wDo[3]*fldDo[4]*eUp-1.7320508075688772*wTar[3]*fldDo[4]*eLo+1.7320508075688772*wDo[3]*fldDo[4]*eLo)/dxDo[3]+0.5*fldDo[0]*eUp-0.5*fldDo[0]*eLo); 
  atomicAdd(&fldTar[1], (dxTar[3]*(0.4330127018922193*fldDo[8]*eUpR2-0.4330127018922193*fldDo[8]*eLoR2))/dxDo[3]+(1.7320508075688772*wTar[3]*fldDo[8]*eUp-1.7320508075688772*wDo[3]*fldDo[8]*eUp-1.7320508075688772*wTar[3]*fldDo[8]*eLo+1.7320508075688772*wDo[3]*fldDo[8]*eLo)/dxDo[3]+0.5*fldDo[1]*eUp-0.5*fldDo[1]*eLo); 
  atomicAdd(&fldTar[2], (dxTar[3]*(0.4330127018922193*fldDo[9]*eUpR2-0.4330127018922193*fldDo[9]*eLoR2))/dxDo[3]+(1.7320508075688772*wTar[3]*fldDo[9]*eUp-1.7320508075688772*wDo[3]*fldDo[9]*eUp-1.7320508075688772*wTar[3]*fldDo[9]*eLo+1.7320508075688772*wDo[3]*fldDo[9]*eLo)/dxDo[3]+0.5*fldDo[2]*eUp-0.5*fldDo[2]*eLo); 
  atomicAdd(&fldTar[3], (dxTar[3]*(0.4330127018922193*fldDo[10]*eUpR2-0.4330127018922193*fldDo[10]*eLoR2))/dxDo[3]+(1.7320508075688772*wTar[3]*fldDo[10]*eUp-1.7320508075688772*wDo[3]*fldDo[10]*eUp-1.7320508075688772*wTar[3]*fldDo[10]*eLo+1.7320508075688772*wDo[3]*fldDo[10]*eLo)/dxDo[3]+0.5*fldDo[3]*eUp-0.5*fldDo[3]*eLo); 
  atomicAdd(&fldTar[4], (dxTar[3]*(0.5*fldDo[4]*eUpR3-0.5*fldDo[4]*eLoR3))/dxDo[3]+(1.5*wTar[3]*fldDo[4]*eUpR2-1.5*wDo[3]*fldDo[4]*eUpR2-1.5*wTar[3]*fldDo[4]*eLoR2+1.5*wDo[3]*fldDo[4]*eLoR2)/dxDo[3]+0.4330127018922193*fldDo[0]*eUpR2-0.4330127018922193*fldDo[0]*eLoR2); 
  atomicAdd(&fldTar[5], (dxTar[3]*(0.4330127018922193*fldDo[12]*eUpR2-0.4330127018922193*fldDo[12]*eLoR2))/dxDo[3]+(1.7320508075688772*wTar[3]*fldDo[12]*eUp-1.7320508075688772*wDo[3]*fldDo[12]*eUp-1.7320508075688772*wTar[3]*fldDo[12]*eLo+1.7320508075688772*wDo[3]*fldDo[12]*eLo)/dxDo[3]+0.5*fldDo[5]*eUp-0.5*fldDo[5]*eLo); 
  atomicAdd(&fldTar[6], (dxTar[3]*(0.4330127018922193*fldDo[13]*eUpR2-0.4330127018922193*fldDo[13]*eLoR2))/dxDo[3]+(1.7320508075688772*wTar[3]*fldDo[13]*eUp-1.7320508075688772*wDo[3]*fldDo[13]*eUp-1.7320508075688772*wTar[3]*fldDo[13]*eLo+1.7320508075688772*wDo[3]*fldDo[13]*eLo)/dxDo[3]+0.5*fldDo[6]*eUp-0.5*fldDo[6]*eLo); 
  atomicAdd(&fldTar[7], (dxTar[3]*(0.4330127018922193*fldDo[14]*eUpR2-0.4330127018922193*fldDo[14]*eLoR2))/dxDo[3]+(1.7320508075688772*wTar[3]*fldDo[14]*eUp-1.7320508075688772*wDo[3]*fldDo[14]*eUp-1.7320508075688772*wTar[3]*fldDo[14]*eLo+1.7320508075688772*wDo[3]*fldDo[14]*eLo)/dxDo[3]+0.5*fldDo[7]*eUp-0.5*fldDo[7]*eLo); 
  atomicAdd(&fldTar[8], (dxTar[3]*(0.5*fldDo[8]*eUpR3-0.5*fldDo[8]*eLoR3))/dxDo[3]+(1.5*wTar[3]*fldDo[8]*eUpR2-1.5*wDo[3]*fldDo[8]*eUpR2-1.5*wTar[3]*fldDo[8]*eLoR2+1.5*wDo[3]*fldDo[8]*eLoR2)/dxDo[3]+0.4330127018922193*fldDo[1]*eUpR2-0.4330127018922193*fldDo[1]*eLoR2); 
  atomicAdd(&fldTar[9], (dxTar[3]*(0.5*fldDo[9]*eUpR3-0.5*fldDo[9]*eLoR3))/dxDo[3]+(1.5*wTar[3]*fldDo[9]*eUpR2-1.5*wDo[3]*fldDo[9]*eUpR2-1.5*wTar[3]*fldDo[9]*eLoR2+1.5*wDo[3]*fldDo[9]*eLoR2)/dxDo[3]+0.4330127018922193*fldDo[2]*eUpR2-0.4330127018922193*fldDo[2]*eLoR2); 
  atomicAdd(&fldTar[10], (dxTar[3]*(0.5*fldDo[10]*eUpR3-0.5*fldDo[10]*eLoR3))/dxDo[3]+(1.5*wTar[3]*fldDo[10]*eUpR2-1.5*wDo[3]*fldDo[10]*eUpR2-1.5*wTar[3]*fldDo[10]*eLoR2+1.5*wDo[3]*fldDo[10]*eLoR2)/dxDo[3]+0.4330127018922193*fldDo[3]*eUpR2-0.4330127018922193*fldDo[3]*eLoR2); 
  atomicAdd(&fldTar[11], (dxTar[3]*(0.4330127018922193*fldDo[15]*eUpR2-0.4330127018922193*fldDo[15]*eLoR2))/dxDo[3]+(1.7320508075688772*wTar[3]*fldDo[15]*eUp-1.7320508075688772*wDo[3]*fldDo[15]*eUp-1.7320508075688772*wTar[3]*fldDo[15]*eLo+1.7320508075688772*wDo[3]*fldDo[15]*eLo)/dxDo[3]+0.5*fldDo[11]*eUp-0.5*fldDo[11]*eLo); 
  atomicAdd(&fldTar[12], (dxTar[3]*(0.5*fldDo[12]*eUpR3-0.5*fldDo[12]*eLoR3))/dxDo[3]+(1.5*wTar[3]*fldDo[12]*eUpR2-1.5*wDo[3]*fldDo[12]*eUpR2-1.5*wTar[3]*fldDo[12]*eLoR2+1.5*wDo[3]*fldDo[12]*eLoR2)/dxDo[3]+0.4330127018922193*fldDo[5]*eUpR2-0.4330127018922193*fldDo[5]*eLoR2); 
  atomicAdd(&fldTar[13], (dxTar[3]*(0.5*fldDo[13]*eUpR3-0.5*fldDo[13]*eLoR3))/dxDo[3]+(1.5*wTar[3]*fldDo[13]*eUpR2-1.5*wDo[3]*fldDo[13]*eUpR2-1.5*wTar[3]*fldDo[13]*eLoR2+1.5*wDo[3]*fldDo[13]*eLoR2)/dxDo[3]+0.4330127018922193*fldDo[6]*eUpR2-0.4330127018922193*fldDo[6]*eLoR2); 
  atomicAdd(&fldTar[14], (dxTar[3]*(0.5*fldDo[14]*eUpR3-0.5*fldDo[14]*eLoR3))/dxDo[3]+(1.5*wTar[3]*fldDo[14]*eUpR2-1.5*wDo[3]*fldDo[14]*eUpR2-1.5*wTar[3]*fldDo[14]*eLoR2+1.5*wDo[3]*fldDo[14]*eLoR2)/dxDo[3]+0.4330127018922193*fldDo[7]*eUpR2-0.4330127018922193*fldDo[7]*eLoR2); 
  atomicAdd(&fldTar[15], (dxTar[3]*(0.5*fldDo[15]*eUpR3-0.5*fldDo[15]*eLoR3))/dxDo[3]+(1.5*wTar[3]*fldDo[15]*eUpR2-1.5*wDo[3]*fldDo[15]*eUpR2-1.5*wTar[3]*fldDo[15]*eLoR2+1.5*wDo[3]*fldDo[15]*eLoR2)/dxDo[3]+0.4330127018922193*fldDo[11]*eUpR2-0.4330127018922193*fldDo[11]*eLoR2); 
#else
  const double eLoR2 = pow(eLo,2);
  const double eLoR3 = pow(eLo,3);
  const double eUpR2 = pow(eUp,2);
  const double eUpR3 = pow(eUp,3);

  fldTar[0] += (dxTar[3]*(0.4330127018922193*fldDo[4]*eUpR2-0.4330127018922193*fldDo[4]*eLoR2))/dxDo[3]+(1.7320508075688772*wTar[3]*fldDo[4]*eUp-1.7320508075688772*wDo[3]*fldDo[4]*eUp-1.7320508075688772*wTar[3]*fldDo[4]*eLo+1.7320508075688772*wDo[3]*fldDo[4]*eLo)/dxDo[3]+0.5*fldDo[0]*eUp-0.5*fldDo[0]*eLo; 
  fldTar[1] += (dxTar[3]*(0.4330127018922193*fldDo[8]*eUpR2-0.4330127018922193*fldDo[8]*eLoR2))/dxDo[3]+(1.7320508075688772*wTar[3]*fldDo[8]*eUp-1.7320508075688772*wDo[3]*fldDo[8]*eUp-1.7320508075688772*wTar[3]*fldDo[8]*eLo+1.7320508075688772*wDo[3]*fldDo[8]*eLo)/dxDo[3]+0.5*fldDo[1]*eUp-0.5*fldDo[1]*eLo; 
  fldTar[2] += (dxTar[3]*(0.4330127018922193*fldDo[9]*eUpR2-0.4330127018922193*fldDo[9]*eLoR2))/dxDo[3]+(1.7320508075688772*wTar[3]*fldDo[9]*eUp-1.7320508075688772*wDo[3]*fldDo[9]*eUp-1.7320508075688772*wTar[3]*fldDo[9]*eLo+1.7320508075688772*wDo[3]*fldDo[9]*eLo)/dxDo[3]+0.5*fldDo[2]*eUp-0.5*fldDo[2]*eLo; 
  fldTar[3] += (dxTar[3]*(0.4330127018922193*fldDo[10]*eUpR2-0.4330127018922193*fldDo[10]*eLoR2))/dxDo[3]+(1.7320508075688772*wTar[3]*fldDo[10]*eUp-1.7320508075688772*wDo[3]*fldDo[10]*eUp-1.7320508075688772*wTar[3]*fldDo[10]*eLo+1.7320508075688772*wDo[3]*fldDo[10]*eLo)/dxDo[3]+0.5*fldDo[3]*eUp-0.5*fldDo[3]*eLo; 
  fldTar[4] += (dxTar[3]*(0.5*fldDo[4]*eUpR3-0.5*fldDo[4]*eLoR3))/dxDo[3]+(1.5*wTar[3]*fldDo[4]*eUpR2-1.5*wDo[3]*fldDo[4]*eUpR2-1.5*wTar[3]*fldDo[4]*eLoR2+1.5*wDo[3]*fldDo[4]*eLoR2)/dxDo[3]+0.4330127018922193*fldDo[0]*eUpR2-0.4330127018922193*fldDo[0]*eLoR2; 
  fldTar[5] += (dxTar[3]*(0.4330127018922193*fldDo[12]*eUpR2-0.4330127018922193*fldDo[12]*eLoR2))/dxDo[3]+(1.7320508075688772*wTar[3]*fldDo[12]*eUp-1.7320508075688772*wDo[3]*fldDo[12]*eUp-1.7320508075688772*wTar[3]*fldDo[12]*eLo+1.7320508075688772*wDo[3]*fldDo[12]*eLo)/dxDo[3]+0.5*fldDo[5]*eUp-0.5*fldDo[5]*eLo; 
  fldTar[6] += (dxTar[3]*(0.4330127018922193*fldDo[13]*eUpR2-0.4330127018922193*fldDo[13]*eLoR2))/dxDo[3]+(1.7320508075688772*wTar[3]*fldDo[13]*eUp-1.7320508075688772*wDo[3]*fldDo[13]*eUp-1.7320508075688772*wTar[3]*fldDo[13]*eLo+1.7320508075688772*wDo[3]*fldDo[13]*eLo)/dxDo[3]+0.5*fldDo[6]*eUp-0.5*fldDo[6]*eLo; 
  fldTar[7] += (dxTar[3]*(0.4330127018922193*fldDo[14]*eUpR2-0.4330127018922193*fldDo[14]*eLoR2))/dxDo[3]+(1.7320508075688772*wTar[3]*fldDo[14]*eUp-1.7320508075688772*wDo[3]*fldDo[14]*eUp-1.7320508075688772*wTar[3]*fldDo[14]*eLo+1.7320508075688772*wDo[3]*fldDo[14]*eLo)/dxDo[3]+0.5*fldDo[7]*eUp-0.5*fldDo[7]*eLo; 
  fldTar[8] += (dxTar[3]*(0.5*fldDo[8]*eUpR3-0.5*fldDo[8]*eLoR3))/dxDo[3]+(1.5*wTar[3]*fldDo[8]*eUpR2-1.5*wDo[3]*fldDo[8]*eUpR2-1.5*wTar[3]*fldDo[8]*eLoR2+1.5*wDo[3]*fldDo[8]*eLoR2)/dxDo[3]+0.4330127018922193*fldDo[1]*eUpR2-0.4330127018922193*fldDo[1]*eLoR2; 
  fldTar[9] += (dxTar[3]*(0.5*fldDo[9]*eUpR3-0.5*fldDo[9]*eLoR3))/dxDo[3]+(1.5*wTar[3]*fldDo[9]*eUpR2-1.5*wDo[3]*fldDo[9]*eUpR2-1.5*wTar[3]*fldDo[9]*eLoR2+1.5*wDo[3]*fldDo[9]*eLoR2)/dxDo[3]+0.4330127018922193*fldDo[2]*eUpR2-0.4330127018922193*fldDo[2]*eLoR2; 
  fldTar[10] += (dxTar[3]*(0.5*fldDo[10]*eUpR3-0.5*fldDo[10]*eLoR3))/dxDo[3]+(1.5*wTar[3]*fldDo[10]*eUpR2-1.5*wDo[3]*fldDo[10]*eUpR2-1.5*wTar[3]*fldDo[10]*eLoR2+1.5*wDo[3]*fldDo[10]*eLoR2)/dxDo[3]+0.4330127018922193*fldDo[3]*eUpR2-0.4330127018922193*fldDo[3]*eLoR2; 
  fldTar[11] += (dxTar[3]*(0.4330127018922193*fldDo[15]*eUpR2-0.4330127018922193*fldDo[15]*eLoR2))/dxDo[3]+(1.7320508075688772*wTar[3]*fldDo[15]*eUp-1.7320508075688772*wDo[3]*fldDo[15]*eUp-1.7320508075688772*wTar[3]*fldDo[15]*eLo+1.7320508075688772*wDo[3]*fldDo[15]*eLo)/dxDo[3]+0.5*fldDo[11]*eUp-0.5*fldDo[11]*eLo; 
  fldTar[12] += (dxTar[3]*(0.5*fldDo[12]*eUpR3-0.5*fldDo[12]*eLoR3))/dxDo[3]+(1.5*wTar[3]*fldDo[12]*eUpR2-1.5*wDo[3]*fldDo[12]*eUpR2-1.5*wTar[3]*fldDo[12]*eLoR2+1.5*wDo[3]*fldDo[12]*eLoR2)/dxDo[3]+0.4330127018922193*fldDo[5]*eUpR2-0.4330127018922193*fldDo[5]*eLoR2; 
  fldTar[13] += (dxTar[3]*(0.5*fldDo[13]*eUpR3-0.5*fldDo[13]*eLoR3))/dxDo[3]+(1.5*wTar[3]*fldDo[13]*eUpR2-1.5*wDo[3]*fldDo[13]*eUpR2-1.5*wTar[3]*fldDo[13]*eLoR2+1.5*wDo[3]*fldDo[13]*eLoR2)/dxDo[3]+0.4330127018922193*fldDo[6]*eUpR2-0.4330127018922193*fldDo[6]*eLoR2; 
  fldTar[14] += (dxTar[3]*(0.5*fldDo[14]*eUpR3-0.5*fldDo[14]*eLoR3))/dxDo[3]+(1.5*wTar[3]*fldDo[14]*eUpR2-1.5*wDo[3]*fldDo[14]*eUpR2-1.5*wTar[3]*fldDo[14]*eLoR2+1.5*wDo[3]*fldDo[14]*eLoR2)/dxDo[3]+0.4330127018922193*fldDo[7]*eUpR2-0.4330127018922193*fldDo[7]*eLoR2; 
  fldTar[15] += (dxTar[3]*(0.5*fldDo[15]*eUpR3-0.5*fldDo[15]*eLoR3))/dxDo[3]+(1.5*wTar[3]*fldDo[15]*eUpR2-1.5*wDo[3]*fldDo[15]*eUpR2-1.5*wTar[3]*fldDo[15]*eLoR2+1.5*wDo[3]*fldDo[15]*eLoR2)/dxDo[3]+0.4330127018922193*fldDo[11]*eUpR2-0.4330127018922193*fldDo[11]*eLoR2; 
#endif

}

