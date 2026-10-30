#include <gkyl_dg_interpolate_kernels.h> 
 
GKYL_CU_DH void dg_interpolate_vlasov_1x2v_ser_p1_x(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar) 
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
  atomicAdd(&fldTar[2], (dxTar[0]*(0.4330127018922193*fldDo[4]*eUpR2-0.4330127018922193*fldDo[4]*eLoR2))/dxDo[0]+(1.7320508075688772*wTar[0]*fldDo[4]*eUp-1.7320508075688772*wDo[0]*fldDo[4]*eUp-1.7320508075688772*wTar[0]*fldDo[4]*eLo+1.7320508075688772*wDo[0]*fldDo[4]*eLo)/dxDo[0]+0.5*fldDo[2]*eUp-0.5*fldDo[2]*eLo); 
  atomicAdd(&fldTar[3], (dxTar[0]*(0.4330127018922193*fldDo[5]*eUpR2-0.4330127018922193*fldDo[5]*eLoR2))/dxDo[0]+(1.7320508075688772*wTar[0]*fldDo[5]*eUp-1.7320508075688772*wDo[0]*fldDo[5]*eUp-1.7320508075688772*wTar[0]*fldDo[5]*eLo+1.7320508075688772*wDo[0]*fldDo[5]*eLo)/dxDo[0]+0.5*fldDo[3]*eUp-0.5*fldDo[3]*eLo); 
  atomicAdd(&fldTar[4], (dxTar[0]*(0.5*fldDo[4]*eUpR3-0.5*fldDo[4]*eLoR3))/dxDo[0]+(1.5*wTar[0]*fldDo[4]*eUpR2-1.5*wDo[0]*fldDo[4]*eUpR2-1.5*wTar[0]*fldDo[4]*eLoR2+1.5*wDo[0]*fldDo[4]*eLoR2)/dxDo[0]+0.4330127018922193*fldDo[2]*eUpR2-0.4330127018922193*fldDo[2]*eLoR2); 
  atomicAdd(&fldTar[5], (dxTar[0]*(0.5*fldDo[5]*eUpR3-0.5*fldDo[5]*eLoR3))/dxDo[0]+(1.5*wTar[0]*fldDo[5]*eUpR2-1.5*wDo[0]*fldDo[5]*eUpR2-1.5*wTar[0]*fldDo[5]*eLoR2+1.5*wDo[0]*fldDo[5]*eLoR2)/dxDo[0]+0.4330127018922193*fldDo[3]*eUpR2-0.4330127018922193*fldDo[3]*eLoR2); 
  atomicAdd(&fldTar[6], (dxTar[0]*(0.4330127018922193*fldDo[7]*eUpR2-0.4330127018922193*fldDo[7]*eLoR2))/dxDo[0]+(1.7320508075688772*wTar[0]*fldDo[7]*eUp-1.7320508075688772*wDo[0]*fldDo[7]*eUp-1.7320508075688772*wTar[0]*fldDo[7]*eLo+1.7320508075688772*wDo[0]*fldDo[7]*eLo)/dxDo[0]+0.5*fldDo[6]*eUp-0.5*fldDo[6]*eLo); 
  atomicAdd(&fldTar[7], (dxTar[0]*(0.5*fldDo[7]*eUpR3-0.5*fldDo[7]*eLoR3))/dxDo[0]+(1.5*wTar[0]*fldDo[7]*eUpR2-1.5*wDo[0]*fldDo[7]*eUpR2-1.5*wTar[0]*fldDo[7]*eLoR2+1.5*wDo[0]*fldDo[7]*eLoR2)/dxDo[0]+0.4330127018922193*fldDo[6]*eUpR2-0.4330127018922193*fldDo[6]*eLoR2); 
#else
  const double eLoR2 = pow(eLo,2);
  const double eLoR3 = pow(eLo,3);
  const double eUpR2 = pow(eUp,2);
  const double eUpR3 = pow(eUp,3);

  fldTar[0] += (dxTar[0]*(0.4330127018922193*fldDo[1]*eUpR2-0.4330127018922193*fldDo[1]*eLoR2))/dxDo[0]+(1.7320508075688772*wTar[0]*fldDo[1]*eUp-1.7320508075688772*wDo[0]*fldDo[1]*eUp-1.7320508075688772*wTar[0]*fldDo[1]*eLo+1.7320508075688772*wDo[0]*fldDo[1]*eLo)/dxDo[0]+0.5*fldDo[0]*eUp-0.5*fldDo[0]*eLo; 
  fldTar[1] += (dxTar[0]*(0.5*fldDo[1]*eUpR3-0.5*fldDo[1]*eLoR3))/dxDo[0]+(1.5*wTar[0]*fldDo[1]*eUpR2-1.5*wDo[0]*fldDo[1]*eUpR2-1.5*wTar[0]*fldDo[1]*eLoR2+1.5*wDo[0]*fldDo[1]*eLoR2)/dxDo[0]+0.4330127018922193*fldDo[0]*eUpR2-0.4330127018922193*fldDo[0]*eLoR2; 
  fldTar[2] += (dxTar[0]*(0.4330127018922193*fldDo[4]*eUpR2-0.4330127018922193*fldDo[4]*eLoR2))/dxDo[0]+(1.7320508075688772*wTar[0]*fldDo[4]*eUp-1.7320508075688772*wDo[0]*fldDo[4]*eUp-1.7320508075688772*wTar[0]*fldDo[4]*eLo+1.7320508075688772*wDo[0]*fldDo[4]*eLo)/dxDo[0]+0.5*fldDo[2]*eUp-0.5*fldDo[2]*eLo; 
  fldTar[3] += (dxTar[0]*(0.4330127018922193*fldDo[5]*eUpR2-0.4330127018922193*fldDo[5]*eLoR2))/dxDo[0]+(1.7320508075688772*wTar[0]*fldDo[5]*eUp-1.7320508075688772*wDo[0]*fldDo[5]*eUp-1.7320508075688772*wTar[0]*fldDo[5]*eLo+1.7320508075688772*wDo[0]*fldDo[5]*eLo)/dxDo[0]+0.5*fldDo[3]*eUp-0.5*fldDo[3]*eLo; 
  fldTar[4] += (dxTar[0]*(0.5*fldDo[4]*eUpR3-0.5*fldDo[4]*eLoR3))/dxDo[0]+(1.5*wTar[0]*fldDo[4]*eUpR2-1.5*wDo[0]*fldDo[4]*eUpR2-1.5*wTar[0]*fldDo[4]*eLoR2+1.5*wDo[0]*fldDo[4]*eLoR2)/dxDo[0]+0.4330127018922193*fldDo[2]*eUpR2-0.4330127018922193*fldDo[2]*eLoR2; 
  fldTar[5] += (dxTar[0]*(0.5*fldDo[5]*eUpR3-0.5*fldDo[5]*eLoR3))/dxDo[0]+(1.5*wTar[0]*fldDo[5]*eUpR2-1.5*wDo[0]*fldDo[5]*eUpR2-1.5*wTar[0]*fldDo[5]*eLoR2+1.5*wDo[0]*fldDo[5]*eLoR2)/dxDo[0]+0.4330127018922193*fldDo[3]*eUpR2-0.4330127018922193*fldDo[3]*eLoR2; 
  fldTar[6] += (dxTar[0]*(0.4330127018922193*fldDo[7]*eUpR2-0.4330127018922193*fldDo[7]*eLoR2))/dxDo[0]+(1.7320508075688772*wTar[0]*fldDo[7]*eUp-1.7320508075688772*wDo[0]*fldDo[7]*eUp-1.7320508075688772*wTar[0]*fldDo[7]*eLo+1.7320508075688772*wDo[0]*fldDo[7]*eLo)/dxDo[0]+0.5*fldDo[6]*eUp-0.5*fldDo[6]*eLo; 
  fldTar[7] += (dxTar[0]*(0.5*fldDo[7]*eUpR3-0.5*fldDo[7]*eLoR3))/dxDo[0]+(1.5*wTar[0]*fldDo[7]*eUpR2-1.5*wDo[0]*fldDo[7]*eUpR2-1.5*wTar[0]*fldDo[7]*eLoR2+1.5*wDo[0]*fldDo[7]*eLoR2)/dxDo[0]+0.4330127018922193*fldDo[6]*eUpR2-0.4330127018922193*fldDo[6]*eLoR2; 
#endif

}

GKYL_CU_DH void dg_interpolate_vlasov_1x2v_ser_p1_vx(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar) 
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
  atomicAdd(&fldTar[1], (dxTar[1]*(0.4330127018922193*fldDo[4]*eUpR2-0.4330127018922193*fldDo[4]*eLoR2))/dxDo[1]+(1.7320508075688772*wTar[1]*fldDo[4]*eUp-1.7320508075688772*wDo[1]*fldDo[4]*eUp-1.7320508075688772*wTar[1]*fldDo[4]*eLo+1.7320508075688772*wDo[1]*fldDo[4]*eLo)/dxDo[1]+0.5*fldDo[1]*eUp-0.5*fldDo[1]*eLo); 
  atomicAdd(&fldTar[2], (dxTar[1]*(0.5*fldDo[2]*eUpR3-0.5*fldDo[2]*eLoR3))/dxDo[1]+(1.5*wTar[1]*fldDo[2]*eUpR2-1.5*wDo[1]*fldDo[2]*eUpR2-1.5*wTar[1]*fldDo[2]*eLoR2+1.5*wDo[1]*fldDo[2]*eLoR2)/dxDo[1]+0.4330127018922193*fldDo[0]*eUpR2-0.4330127018922193*fldDo[0]*eLoR2); 
  atomicAdd(&fldTar[3], (dxTar[1]*(0.4330127018922193*fldDo[6]*eUpR2-0.4330127018922193*fldDo[6]*eLoR2))/dxDo[1]+(1.7320508075688772*wTar[1]*fldDo[6]*eUp-1.7320508075688772*wDo[1]*fldDo[6]*eUp-1.7320508075688772*wTar[1]*fldDo[6]*eLo+1.7320508075688772*wDo[1]*fldDo[6]*eLo)/dxDo[1]+0.5*fldDo[3]*eUp-0.5*fldDo[3]*eLo); 
  atomicAdd(&fldTar[4], (dxTar[1]*(0.5*fldDo[4]*eUpR3-0.5*fldDo[4]*eLoR3))/dxDo[1]+(1.5*wTar[1]*fldDo[4]*eUpR2-1.5*wDo[1]*fldDo[4]*eUpR2-1.5*wTar[1]*fldDo[4]*eLoR2+1.5*wDo[1]*fldDo[4]*eLoR2)/dxDo[1]+0.4330127018922193*fldDo[1]*eUpR2-0.4330127018922193*fldDo[1]*eLoR2); 
  atomicAdd(&fldTar[5], (dxTar[1]*(0.4330127018922193*fldDo[7]*eUpR2-0.4330127018922193*fldDo[7]*eLoR2))/dxDo[1]+(1.7320508075688772*wTar[1]*fldDo[7]*eUp-1.7320508075688772*wDo[1]*fldDo[7]*eUp-1.7320508075688772*wTar[1]*fldDo[7]*eLo+1.7320508075688772*wDo[1]*fldDo[7]*eLo)/dxDo[1]+0.5*fldDo[5]*eUp-0.5*fldDo[5]*eLo); 
  atomicAdd(&fldTar[6], (dxTar[1]*(0.5*fldDo[6]*eUpR3-0.5*fldDo[6]*eLoR3))/dxDo[1]+(1.5*wTar[1]*fldDo[6]*eUpR2-1.5*wDo[1]*fldDo[6]*eUpR2-1.5*wTar[1]*fldDo[6]*eLoR2+1.5*wDo[1]*fldDo[6]*eLoR2)/dxDo[1]+0.4330127018922193*fldDo[3]*eUpR2-0.4330127018922193*fldDo[3]*eLoR2); 
  atomicAdd(&fldTar[7], (dxTar[1]*(0.5*fldDo[7]*eUpR3-0.5*fldDo[7]*eLoR3))/dxDo[1]+(1.5*wTar[1]*fldDo[7]*eUpR2-1.5*wDo[1]*fldDo[7]*eUpR2-1.5*wTar[1]*fldDo[7]*eLoR2+1.5*wDo[1]*fldDo[7]*eLoR2)/dxDo[1]+0.4330127018922193*fldDo[5]*eUpR2-0.4330127018922193*fldDo[5]*eLoR2); 
#else
  const double eLoR2 = pow(eLo,2);
  const double eLoR3 = pow(eLo,3);
  const double eUpR2 = pow(eUp,2);
  const double eUpR3 = pow(eUp,3);

  fldTar[0] += (dxTar[1]*(0.4330127018922193*fldDo[2]*eUpR2-0.4330127018922193*fldDo[2]*eLoR2))/dxDo[1]+(1.7320508075688772*wTar[1]*fldDo[2]*eUp-1.7320508075688772*wDo[1]*fldDo[2]*eUp-1.7320508075688772*wTar[1]*fldDo[2]*eLo+1.7320508075688772*wDo[1]*fldDo[2]*eLo)/dxDo[1]+0.5*fldDo[0]*eUp-0.5*fldDo[0]*eLo; 
  fldTar[1] += (dxTar[1]*(0.4330127018922193*fldDo[4]*eUpR2-0.4330127018922193*fldDo[4]*eLoR2))/dxDo[1]+(1.7320508075688772*wTar[1]*fldDo[4]*eUp-1.7320508075688772*wDo[1]*fldDo[4]*eUp-1.7320508075688772*wTar[1]*fldDo[4]*eLo+1.7320508075688772*wDo[1]*fldDo[4]*eLo)/dxDo[1]+0.5*fldDo[1]*eUp-0.5*fldDo[1]*eLo; 
  fldTar[2] += (dxTar[1]*(0.5*fldDo[2]*eUpR3-0.5*fldDo[2]*eLoR3))/dxDo[1]+(1.5*wTar[1]*fldDo[2]*eUpR2-1.5*wDo[1]*fldDo[2]*eUpR2-1.5*wTar[1]*fldDo[2]*eLoR2+1.5*wDo[1]*fldDo[2]*eLoR2)/dxDo[1]+0.4330127018922193*fldDo[0]*eUpR2-0.4330127018922193*fldDo[0]*eLoR2; 
  fldTar[3] += (dxTar[1]*(0.4330127018922193*fldDo[6]*eUpR2-0.4330127018922193*fldDo[6]*eLoR2))/dxDo[1]+(1.7320508075688772*wTar[1]*fldDo[6]*eUp-1.7320508075688772*wDo[1]*fldDo[6]*eUp-1.7320508075688772*wTar[1]*fldDo[6]*eLo+1.7320508075688772*wDo[1]*fldDo[6]*eLo)/dxDo[1]+0.5*fldDo[3]*eUp-0.5*fldDo[3]*eLo; 
  fldTar[4] += (dxTar[1]*(0.5*fldDo[4]*eUpR3-0.5*fldDo[4]*eLoR3))/dxDo[1]+(1.5*wTar[1]*fldDo[4]*eUpR2-1.5*wDo[1]*fldDo[4]*eUpR2-1.5*wTar[1]*fldDo[4]*eLoR2+1.5*wDo[1]*fldDo[4]*eLoR2)/dxDo[1]+0.4330127018922193*fldDo[1]*eUpR2-0.4330127018922193*fldDo[1]*eLoR2; 
  fldTar[5] += (dxTar[1]*(0.4330127018922193*fldDo[7]*eUpR2-0.4330127018922193*fldDo[7]*eLoR2))/dxDo[1]+(1.7320508075688772*wTar[1]*fldDo[7]*eUp-1.7320508075688772*wDo[1]*fldDo[7]*eUp-1.7320508075688772*wTar[1]*fldDo[7]*eLo+1.7320508075688772*wDo[1]*fldDo[7]*eLo)/dxDo[1]+0.5*fldDo[5]*eUp-0.5*fldDo[5]*eLo; 
  fldTar[6] += (dxTar[1]*(0.5*fldDo[6]*eUpR3-0.5*fldDo[6]*eLoR3))/dxDo[1]+(1.5*wTar[1]*fldDo[6]*eUpR2-1.5*wDo[1]*fldDo[6]*eUpR2-1.5*wTar[1]*fldDo[6]*eLoR2+1.5*wDo[1]*fldDo[6]*eLoR2)/dxDo[1]+0.4330127018922193*fldDo[3]*eUpR2-0.4330127018922193*fldDo[3]*eLoR2; 
  fldTar[7] += (dxTar[1]*(0.5*fldDo[7]*eUpR3-0.5*fldDo[7]*eLoR3))/dxDo[1]+(1.5*wTar[1]*fldDo[7]*eUpR2-1.5*wDo[1]*fldDo[7]*eUpR2-1.5*wTar[1]*fldDo[7]*eLoR2+1.5*wDo[1]*fldDo[7]*eLoR2)/dxDo[1]+0.4330127018922193*fldDo[5]*eUpR2-0.4330127018922193*fldDo[5]*eLoR2; 
#endif

}

GKYL_CU_DH void dg_interpolate_vlasov_1x2v_ser_p1_vy(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar) 
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
  atomicAdd(&fldTar[1], (dxTar[2]*(0.4330127018922193*fldDo[5]*eUpR2-0.4330127018922193*fldDo[5]*eLoR2))/dxDo[2]+(1.7320508075688772*wTar[2]*fldDo[5]*eUp-1.7320508075688772*wDo[2]*fldDo[5]*eUp-1.7320508075688772*wTar[2]*fldDo[5]*eLo+1.7320508075688772*wDo[2]*fldDo[5]*eLo)/dxDo[2]+0.5*fldDo[1]*eUp-0.5*fldDo[1]*eLo); 
  atomicAdd(&fldTar[2], (dxTar[2]*(0.4330127018922193*fldDo[6]*eUpR2-0.4330127018922193*fldDo[6]*eLoR2))/dxDo[2]+(1.7320508075688772*wTar[2]*fldDo[6]*eUp-1.7320508075688772*wDo[2]*fldDo[6]*eUp-1.7320508075688772*wTar[2]*fldDo[6]*eLo+1.7320508075688772*wDo[2]*fldDo[6]*eLo)/dxDo[2]+0.5*fldDo[2]*eUp-0.5*fldDo[2]*eLo); 
  atomicAdd(&fldTar[3], (dxTar[2]*(0.5*fldDo[3]*eUpR3-0.5*fldDo[3]*eLoR3))/dxDo[2]+(1.5*wTar[2]*fldDo[3]*eUpR2-1.5*wDo[2]*fldDo[3]*eUpR2-1.5*wTar[2]*fldDo[3]*eLoR2+1.5*wDo[2]*fldDo[3]*eLoR2)/dxDo[2]+0.4330127018922193*fldDo[0]*eUpR2-0.4330127018922193*fldDo[0]*eLoR2); 
  atomicAdd(&fldTar[4], (dxTar[2]*(0.4330127018922193*fldDo[7]*eUpR2-0.4330127018922193*fldDo[7]*eLoR2))/dxDo[2]+(1.7320508075688772*wTar[2]*fldDo[7]*eUp-1.7320508075688772*wDo[2]*fldDo[7]*eUp-1.7320508075688772*wTar[2]*fldDo[7]*eLo+1.7320508075688772*wDo[2]*fldDo[7]*eLo)/dxDo[2]+0.5*fldDo[4]*eUp-0.5*fldDo[4]*eLo); 
  atomicAdd(&fldTar[5], (dxTar[2]*(0.5*fldDo[5]*eUpR3-0.5*fldDo[5]*eLoR3))/dxDo[2]+(1.5*wTar[2]*fldDo[5]*eUpR2-1.5*wDo[2]*fldDo[5]*eUpR2-1.5*wTar[2]*fldDo[5]*eLoR2+1.5*wDo[2]*fldDo[5]*eLoR2)/dxDo[2]+0.4330127018922193*fldDo[1]*eUpR2-0.4330127018922193*fldDo[1]*eLoR2); 
  atomicAdd(&fldTar[6], (dxTar[2]*(0.5*fldDo[6]*eUpR3-0.5*fldDo[6]*eLoR3))/dxDo[2]+(1.5*wTar[2]*fldDo[6]*eUpR2-1.5*wDo[2]*fldDo[6]*eUpR2-1.5*wTar[2]*fldDo[6]*eLoR2+1.5*wDo[2]*fldDo[6]*eLoR2)/dxDo[2]+0.4330127018922193*fldDo[2]*eUpR2-0.4330127018922193*fldDo[2]*eLoR2); 
  atomicAdd(&fldTar[7], (dxTar[2]*(0.5*fldDo[7]*eUpR3-0.5*fldDo[7]*eLoR3))/dxDo[2]+(1.5*wTar[2]*fldDo[7]*eUpR2-1.5*wDo[2]*fldDo[7]*eUpR2-1.5*wTar[2]*fldDo[7]*eLoR2+1.5*wDo[2]*fldDo[7]*eLoR2)/dxDo[2]+0.4330127018922193*fldDo[4]*eUpR2-0.4330127018922193*fldDo[4]*eLoR2); 
#else
  const double eLoR2 = pow(eLo,2);
  const double eLoR3 = pow(eLo,3);
  const double eUpR2 = pow(eUp,2);
  const double eUpR3 = pow(eUp,3);

  fldTar[0] += (dxTar[2]*(0.4330127018922193*fldDo[3]*eUpR2-0.4330127018922193*fldDo[3]*eLoR2))/dxDo[2]+(1.7320508075688772*wTar[2]*fldDo[3]*eUp-1.7320508075688772*wDo[2]*fldDo[3]*eUp-1.7320508075688772*wTar[2]*fldDo[3]*eLo+1.7320508075688772*wDo[2]*fldDo[3]*eLo)/dxDo[2]+0.5*fldDo[0]*eUp-0.5*fldDo[0]*eLo; 
  fldTar[1] += (dxTar[2]*(0.4330127018922193*fldDo[5]*eUpR2-0.4330127018922193*fldDo[5]*eLoR2))/dxDo[2]+(1.7320508075688772*wTar[2]*fldDo[5]*eUp-1.7320508075688772*wDo[2]*fldDo[5]*eUp-1.7320508075688772*wTar[2]*fldDo[5]*eLo+1.7320508075688772*wDo[2]*fldDo[5]*eLo)/dxDo[2]+0.5*fldDo[1]*eUp-0.5*fldDo[1]*eLo; 
  fldTar[2] += (dxTar[2]*(0.4330127018922193*fldDo[6]*eUpR2-0.4330127018922193*fldDo[6]*eLoR2))/dxDo[2]+(1.7320508075688772*wTar[2]*fldDo[6]*eUp-1.7320508075688772*wDo[2]*fldDo[6]*eUp-1.7320508075688772*wTar[2]*fldDo[6]*eLo+1.7320508075688772*wDo[2]*fldDo[6]*eLo)/dxDo[2]+0.5*fldDo[2]*eUp-0.5*fldDo[2]*eLo; 
  fldTar[3] += (dxTar[2]*(0.5*fldDo[3]*eUpR3-0.5*fldDo[3]*eLoR3))/dxDo[2]+(1.5*wTar[2]*fldDo[3]*eUpR2-1.5*wDo[2]*fldDo[3]*eUpR2-1.5*wTar[2]*fldDo[3]*eLoR2+1.5*wDo[2]*fldDo[3]*eLoR2)/dxDo[2]+0.4330127018922193*fldDo[0]*eUpR2-0.4330127018922193*fldDo[0]*eLoR2; 
  fldTar[4] += (dxTar[2]*(0.4330127018922193*fldDo[7]*eUpR2-0.4330127018922193*fldDo[7]*eLoR2))/dxDo[2]+(1.7320508075688772*wTar[2]*fldDo[7]*eUp-1.7320508075688772*wDo[2]*fldDo[7]*eUp-1.7320508075688772*wTar[2]*fldDo[7]*eLo+1.7320508075688772*wDo[2]*fldDo[7]*eLo)/dxDo[2]+0.5*fldDo[4]*eUp-0.5*fldDo[4]*eLo; 
  fldTar[5] += (dxTar[2]*(0.5*fldDo[5]*eUpR3-0.5*fldDo[5]*eLoR3))/dxDo[2]+(1.5*wTar[2]*fldDo[5]*eUpR2-1.5*wDo[2]*fldDo[5]*eUpR2-1.5*wTar[2]*fldDo[5]*eLoR2+1.5*wDo[2]*fldDo[5]*eLoR2)/dxDo[2]+0.4330127018922193*fldDo[1]*eUpR2-0.4330127018922193*fldDo[1]*eLoR2; 
  fldTar[6] += (dxTar[2]*(0.5*fldDo[6]*eUpR3-0.5*fldDo[6]*eLoR3))/dxDo[2]+(1.5*wTar[2]*fldDo[6]*eUpR2-1.5*wDo[2]*fldDo[6]*eUpR2-1.5*wTar[2]*fldDo[6]*eLoR2+1.5*wDo[2]*fldDo[6]*eLoR2)/dxDo[2]+0.4330127018922193*fldDo[2]*eUpR2-0.4330127018922193*fldDo[2]*eLoR2; 
  fldTar[7] += (dxTar[2]*(0.5*fldDo[7]*eUpR3-0.5*fldDo[7]*eLoR3))/dxDo[2]+(1.5*wTar[2]*fldDo[7]*eUpR2-1.5*wDo[2]*fldDo[7]*eUpR2-1.5*wTar[2]*fldDo[7]*eLoR2+1.5*wDo[2]*fldDo[7]*eLoR2)/dxDo[2]+0.4330127018922193*fldDo[4]*eUpR2-0.4330127018922193*fldDo[4]*eLoR2; 
#endif

}

