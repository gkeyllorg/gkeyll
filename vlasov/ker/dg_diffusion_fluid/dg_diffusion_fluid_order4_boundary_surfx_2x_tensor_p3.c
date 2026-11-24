#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order4_boundary_surfx_2x_tensor_p3_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],4.);

  double vol_incr[16] = {0.0}; 

  double edgeSurf_incr[16] = {0.0}; 
  double boundSurf_incr[16] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 15.347424597386395*coeff[0]*fSkin[8]+15.347424597386395*coeff[0]*fEdge[8]+26.59698043549555*coeff[0]*fSkin[4]-26.59698043549555*coeff[0]*fEdge[4]+24.864401241467277*coeff[0]*fSkin[1]+24.864401241467277*coeff[0]*fEdge[1]+14.35546875*coeff[0]*fSkin[0]-14.35546875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 31.492421698315376*coeff[0]*fSkin[8]+21.672616637695718*coeff[0]*fEdge[8]+51.3321581784444*coeff[0]*fSkin[4]-40.80248470594299*coeff[0]*fEdge[4]+45.17578125*coeff[0]*fSkin[1]+40.95703125*coeff[0]*fEdge[1]+24.864401241467277*coeff[0]*fSkin[0]-24.864401241467277*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 15.347424597386388*coeff[0]*fSkin[11]+15.347424597386388*coeff[0]*fEdge[11]+26.596980435495546*coeff[0]*fSkin[6]-26.596980435495546*coeff[0]*fEdge[6]+24.864401241467277*coeff[0]*fSkin[3]+24.864401241467277*coeff[0]*fEdge[3]+14.35546875*coeff[0]*fSkin[2]-14.35546875*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 31.492421698315383*coeff[0]*fSkin[11]+21.67261663769573*coeff[0]*fEdge[11]+51.3321581784444*coeff[0]*fSkin[6]-40.80248470594299*coeff[0]*fEdge[6]+45.17578125*coeff[0]*fSkin[3]+40.95703125*coeff[0]*fEdge[3]+24.864401241467277*coeff[0]*fSkin[2]-24.864401241467277*coeff[0]*fEdge[2]; 
  edgeSurf_incr[4] = 48.75813745345884*coeff[0]*fSkin[8]+10.726195990675592*coeff[0]*fEdge[8]+69.43359375*coeff[0]*fSkin[4]-28.65234375*coeff[0]*fEdge[4]+50.78751989538398*coeff[0]*fSkin[1]+34.44837140357144*coeff[0]*fEdge[1]+22.928431409909958*coeff[0]*fSkin[0]-22.928431409909958*coeff[0]*fEdge[0]; 
  edgeSurf_incr[5] = 15.347424597386397*coeff[0]*fSkin[13]+15.347424597386397*coeff[0]*fEdge[13]+26.59698043549555*coeff[0]*fSkin[10]-26.59698043549555*coeff[0]*fEdge[10]+24.864401241467288*coeff[0]*fSkin[7]+24.864401241467288*coeff[0]*fEdge[7]+14.35546875*coeff[0]*fSkin[5]-14.35546875*coeff[0]*fEdge[5]; 
  edgeSurf_incr[6] = 48.758137453458836*coeff[0]*fSkin[11]+10.726195990675588*coeff[0]*fEdge[11]+69.43359375*coeff[0]*fSkin[6]-28.65234375*coeff[0]*fEdge[6]+50.78751989538398*coeff[0]*fSkin[3]+34.44837140357144*coeff[0]*fEdge[3]+22.92843140990995*coeff[0]*fSkin[2]-22.92843140990995*coeff[0]*fEdge[2]; 
  edgeSurf_incr[7] = 31.492421698315383*coeff[0]*fSkin[13]+21.67261663769573*coeff[0]*fEdge[13]+51.3321581784444*coeff[0]*fSkin[10]-40.80248470594299*coeff[0]*fEdge[10]+45.17578125*coeff[0]*fSkin[7]+40.95703125*coeff[0]*fEdge[7]+24.864401241467288*coeff[0]*fSkin[5]-24.864401241467288*coeff[0]*fEdge[5]; 
  edgeSurf_incr[8] = 51.03515625*coeff[0]*fSkin[8]-23.96484375*coeff[0]*fEdge[8]+43.053346234041356*coeff[0]*fSkin[4]+25.72108124449167*coeff[0]*fEdge[4]-6.712757365657978*coeff[0]*fSkin[1]-15.305086793700156*coeff[0]*fEdge[1]-26.199138959174753*coeff[0]*fSkin[0]+6.356004126190321*coeff[0]*fEdge[0]; 
  edgeSurf_incr[9] = 15.347424597386395*coeff[0]*fSkin[15]+15.347424597386395*coeff[0]*fEdge[15]+26.59698043549555*coeff[0]*fSkin[14]-26.59698043549555*coeff[0]*fEdge[14]+24.86440124146728*coeff[0]*fSkin[12]+24.86440124146728*coeff[0]*fEdge[12]+14.35546875*coeff[0]*fSkin[9]-14.35546875*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = 48.75813745345883*coeff[0]*fSkin[13]+10.726195990675588*coeff[0]*fEdge[13]+69.43359375*coeff[0]*fSkin[10]-28.65234375*coeff[0]*fEdge[10]+50.78751989538398*coeff[0]*fSkin[7]+34.44837140357144*coeff[0]*fEdge[7]+22.928431409909958*coeff[0]*fSkin[5]-22.928431409909958*coeff[0]*fEdge[5]; 
  edgeSurf_incr[11] = 51.03515625*coeff[0]*fSkin[11]-23.96484375*coeff[0]*fEdge[11]+43.05334623404135*coeff[0]*fSkin[6]+25.721081244491685*coeff[0]*fEdge[6]-6.712757365657971*coeff[0]*fSkin[3]-15.30508679370017*coeff[0]*fEdge[3]-26.19913895917474*coeff[0]*fSkin[2]+6.356004126190324*coeff[0]*fEdge[2]; 
  edgeSurf_incr[12] = 31.492421698315383*coeff[0]*fSkin[15]+21.67261663769573*coeff[0]*fEdge[15]+51.332158178444395*coeff[0]*fSkin[14]-40.80248470594298*coeff[0]*fEdge[14]+45.17578125*coeff[0]*fSkin[12]+40.95703125*coeff[0]*fEdge[12]+24.86440124146728*coeff[0]*fSkin[9]-24.86440124146728*coeff[0]*fEdge[9]; 
  edgeSurf_incr[13] = 51.03515625*coeff[0]*fSkin[13]-23.96484375*coeff[0]*fEdge[13]+43.05334623404134*coeff[0]*fSkin[10]+25.721081244491685*coeff[0]*fEdge[10]-6.7127573656579695*coeff[0]*fSkin[7]-15.30508679370017*coeff[0]*fEdge[7]-26.19913895917476*coeff[0]*fSkin[5]+6.356004126190328*coeff[0]*fEdge[5]; 
  edgeSurf_incr[14] = 48.75813745345883*coeff[0]*fSkin[15]+10.726195990675588*coeff[0]*fEdge[15]+69.43359375*coeff[0]*fSkin[14]-28.65234375*coeff[0]*fEdge[14]+50.78751989538397*coeff[0]*fSkin[12]+34.44837140357143*coeff[0]*fEdge[12]+22.928431409909958*coeff[0]*fSkin[9]-22.928431409909958*coeff[0]*fEdge[9]; 
  edgeSurf_incr[15] = 51.03515625*coeff[0]*fSkin[15]-23.96484375*coeff[0]*fEdge[15]+43.05334623404134*coeff[0]*fSkin[14]+25.721081244491685*coeff[0]*fEdge[14]-6.712757365657971*coeff[0]*fSkin[12]-15.30508679370017*coeff[0]*fEdge[12]-26.199138959174753*coeff[0]*fSkin[9]+6.356004126190321*coeff[0]*fEdge[9]; 

  boundSurf_incr[1] = 11.783766072743584*coeff[0]*fSkin[8]-11.61895003862225*coeff[0]*fSkin[4]+4.5*coeff[0]*fSkin[1]; 
  boundSurf_incr[3] = 11.783766072743589*coeff[0]*fSkin[11]-11.618950038622252*coeff[0]*fSkin[6]+4.5*coeff[0]*fSkin[3]; 
  boundSurf_incr[4] = -(45.638329755339896*coeff[0]*fSkin[8])+45.0*coeff[0]*fSkin[4]-17.42842505793337*coeff[0]*fSkin[1]; 
  boundSurf_incr[6] = -(45.6383297553399*coeff[0]*fSkin[11])+45.0*coeff[0]*fSkin[6]-17.428425057933378*coeff[0]*fSkin[3]; 
  boundSurf_incr[7] = 11.783766072743589*coeff[0]*fSkin[13]-11.618950038622252*coeff[0]*fSkin[10]+4.5*coeff[0]*fSkin[7]; 
  boundSurf_incr[8] = 91.8*coeff[0]*fSkin[8]-78.09225313691493*coeff[0]*fSkin[4]+10.998181667894002*coeff[0]*fSkin[1]+19.84313483298443*coeff[0]*fSkin[0]; 
  boundSurf_incr[10] = -(45.638329755339896*coeff[0]*fSkin[13])+45.0*coeff[0]*fSkin[10]-17.428425057933378*coeff[0]*fSkin[7]; 
  boundSurf_incr[11] = 91.8*coeff[0]*fSkin[11]-78.09225313691493*coeff[0]*fSkin[6]+10.998181667894015*coeff[0]*fSkin[3]+19.843134832984425*coeff[0]*fSkin[2]; 
  boundSurf_incr[12] = 11.783766072743589*coeff[0]*fSkin[15]-11.61895003862225*coeff[0]*fSkin[14]+4.5*coeff[0]*fSkin[12]; 
  boundSurf_incr[13] = 91.8*coeff[0]*fSkin[13]-78.09225313691493*coeff[0]*fSkin[10]+10.998181667894016*coeff[0]*fSkin[7]+19.84313483298443*coeff[0]*fSkin[5]; 
  boundSurf_incr[14] = -(45.638329755339896*coeff[0]*fSkin[15])+45.0*coeff[0]*fSkin[14]-17.428425057933378*coeff[0]*fSkin[12]; 
  boundSurf_incr[15] = 91.8*coeff[0]*fSkin[15]-78.09225313691493*coeff[0]*fSkin[14]+10.998181667894015*coeff[0]*fSkin[12]+19.84313483298443*coeff[0]*fSkin[9]; 

  } else { 

  edgeSurf_incr[0] = -(15.347424597386395*coeff[0]*fSkin[8])-15.347424597386395*coeff[0]*fEdge[8]+26.59698043549555*coeff[0]*fSkin[4]-26.59698043549555*coeff[0]*fEdge[4]-24.864401241467277*coeff[0]*fSkin[1]-24.864401241467277*coeff[0]*fEdge[1]+14.35546875*coeff[0]*fSkin[0]-14.35546875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 31.492421698315376*coeff[0]*fSkin[8]+21.672616637695718*coeff[0]*fEdge[8]-51.3321581784444*coeff[0]*fSkin[4]+40.80248470594299*coeff[0]*fEdge[4]+45.17578125*coeff[0]*fSkin[1]+40.95703125*coeff[0]*fEdge[1]-24.864401241467277*coeff[0]*fSkin[0]+24.864401241467277*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = -(15.347424597386388*coeff[0]*fSkin[11])-15.347424597386388*coeff[0]*fEdge[11]+26.596980435495546*coeff[0]*fSkin[6]-26.596980435495546*coeff[0]*fEdge[6]-24.864401241467277*coeff[0]*fSkin[3]-24.864401241467277*coeff[0]*fEdge[3]+14.35546875*coeff[0]*fSkin[2]-14.35546875*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 31.492421698315383*coeff[0]*fSkin[11]+21.67261663769573*coeff[0]*fEdge[11]-51.3321581784444*coeff[0]*fSkin[6]+40.80248470594299*coeff[0]*fEdge[6]+45.17578125*coeff[0]*fSkin[3]+40.95703125*coeff[0]*fEdge[3]-24.864401241467277*coeff[0]*fSkin[2]+24.864401241467277*coeff[0]*fEdge[2]; 
  edgeSurf_incr[4] = -(48.75813745345884*coeff[0]*fSkin[8])-10.726195990675592*coeff[0]*fEdge[8]+69.43359375*coeff[0]*fSkin[4]-28.65234375*coeff[0]*fEdge[4]-50.78751989538398*coeff[0]*fSkin[1]-34.44837140357144*coeff[0]*fEdge[1]+22.928431409909958*coeff[0]*fSkin[0]-22.928431409909958*coeff[0]*fEdge[0]; 
  edgeSurf_incr[5] = -(15.347424597386397*coeff[0]*fSkin[13])-15.347424597386397*coeff[0]*fEdge[13]+26.59698043549555*coeff[0]*fSkin[10]-26.59698043549555*coeff[0]*fEdge[10]-24.864401241467288*coeff[0]*fSkin[7]-24.864401241467288*coeff[0]*fEdge[7]+14.35546875*coeff[0]*fSkin[5]-14.35546875*coeff[0]*fEdge[5]; 
  edgeSurf_incr[6] = -(48.758137453458836*coeff[0]*fSkin[11])-10.726195990675588*coeff[0]*fEdge[11]+69.43359375*coeff[0]*fSkin[6]-28.65234375*coeff[0]*fEdge[6]-50.78751989538398*coeff[0]*fSkin[3]-34.44837140357144*coeff[0]*fEdge[3]+22.92843140990995*coeff[0]*fSkin[2]-22.92843140990995*coeff[0]*fEdge[2]; 
  edgeSurf_incr[7] = 31.492421698315383*coeff[0]*fSkin[13]+21.67261663769573*coeff[0]*fEdge[13]-51.3321581784444*coeff[0]*fSkin[10]+40.80248470594299*coeff[0]*fEdge[10]+45.17578125*coeff[0]*fSkin[7]+40.95703125*coeff[0]*fEdge[7]-24.864401241467288*coeff[0]*fSkin[5]+24.864401241467288*coeff[0]*fEdge[5]; 
  edgeSurf_incr[8] = 51.03515625*coeff[0]*fSkin[8]-23.96484375*coeff[0]*fEdge[8]-43.053346234041356*coeff[0]*fSkin[4]-25.72108124449167*coeff[0]*fEdge[4]-6.712757365657978*coeff[0]*fSkin[1]-15.305086793700156*coeff[0]*fEdge[1]+26.199138959174753*coeff[0]*fSkin[0]-6.356004126190321*coeff[0]*fEdge[0]; 
  edgeSurf_incr[9] = -(15.347424597386395*coeff[0]*fSkin[15])-15.347424597386395*coeff[0]*fEdge[15]+26.59698043549555*coeff[0]*fSkin[14]-26.59698043549555*coeff[0]*fEdge[14]-24.86440124146728*coeff[0]*fSkin[12]-24.86440124146728*coeff[0]*fEdge[12]+14.35546875*coeff[0]*fSkin[9]-14.35546875*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = -(48.75813745345883*coeff[0]*fSkin[13])-10.726195990675588*coeff[0]*fEdge[13]+69.43359375*coeff[0]*fSkin[10]-28.65234375*coeff[0]*fEdge[10]-50.78751989538398*coeff[0]*fSkin[7]-34.44837140357144*coeff[0]*fEdge[7]+22.928431409909958*coeff[0]*fSkin[5]-22.928431409909958*coeff[0]*fEdge[5]; 
  edgeSurf_incr[11] = 51.03515625*coeff[0]*fSkin[11]-23.96484375*coeff[0]*fEdge[11]-43.05334623404135*coeff[0]*fSkin[6]-25.721081244491685*coeff[0]*fEdge[6]-6.712757365657971*coeff[0]*fSkin[3]-15.30508679370017*coeff[0]*fEdge[3]+26.19913895917474*coeff[0]*fSkin[2]-6.356004126190324*coeff[0]*fEdge[2]; 
  edgeSurf_incr[12] = 31.492421698315383*coeff[0]*fSkin[15]+21.67261663769573*coeff[0]*fEdge[15]-51.332158178444395*coeff[0]*fSkin[14]+40.80248470594298*coeff[0]*fEdge[14]+45.17578125*coeff[0]*fSkin[12]+40.95703125*coeff[0]*fEdge[12]-24.86440124146728*coeff[0]*fSkin[9]+24.86440124146728*coeff[0]*fEdge[9]; 
  edgeSurf_incr[13] = 51.03515625*coeff[0]*fSkin[13]-23.96484375*coeff[0]*fEdge[13]-43.05334623404134*coeff[0]*fSkin[10]-25.721081244491685*coeff[0]*fEdge[10]-6.7127573656579695*coeff[0]*fSkin[7]-15.30508679370017*coeff[0]*fEdge[7]+26.19913895917476*coeff[0]*fSkin[5]-6.356004126190328*coeff[0]*fEdge[5]; 
  edgeSurf_incr[14] = -(48.75813745345883*coeff[0]*fSkin[15])-10.726195990675588*coeff[0]*fEdge[15]+69.43359375*coeff[0]*fSkin[14]-28.65234375*coeff[0]*fEdge[14]-50.78751989538397*coeff[0]*fSkin[12]-34.44837140357143*coeff[0]*fEdge[12]+22.928431409909958*coeff[0]*fSkin[9]-22.928431409909958*coeff[0]*fEdge[9]; 
  edgeSurf_incr[15] = 51.03515625*coeff[0]*fSkin[15]-23.96484375*coeff[0]*fEdge[15]-43.05334623404134*coeff[0]*fSkin[14]-25.721081244491685*coeff[0]*fEdge[14]-6.712757365657971*coeff[0]*fSkin[12]-15.30508679370017*coeff[0]*fEdge[12]+26.199138959174753*coeff[0]*fSkin[9]-6.356004126190321*coeff[0]*fEdge[9]; 

  boundSurf_incr[1] = 11.783766072743584*coeff[0]*fSkin[8]+11.61895003862225*coeff[0]*fSkin[4]+4.5*coeff[0]*fSkin[1]; 
  boundSurf_incr[3] = 11.783766072743589*coeff[0]*fSkin[11]+11.618950038622252*coeff[0]*fSkin[6]+4.5*coeff[0]*fSkin[3]; 
  boundSurf_incr[4] = 45.638329755339896*coeff[0]*fSkin[8]+45.0*coeff[0]*fSkin[4]+17.42842505793337*coeff[0]*fSkin[1]; 
  boundSurf_incr[6] = 45.6383297553399*coeff[0]*fSkin[11]+45.0*coeff[0]*fSkin[6]+17.428425057933378*coeff[0]*fSkin[3]; 
  boundSurf_incr[7] = 11.783766072743589*coeff[0]*fSkin[13]+11.618950038622252*coeff[0]*fSkin[10]+4.5*coeff[0]*fSkin[7]; 
  boundSurf_incr[8] = 91.8*coeff[0]*fSkin[8]+78.09225313691493*coeff[0]*fSkin[4]+10.998181667894002*coeff[0]*fSkin[1]-19.84313483298443*coeff[0]*fSkin[0]; 
  boundSurf_incr[10] = 45.638329755339896*coeff[0]*fSkin[13]+45.0*coeff[0]*fSkin[10]+17.428425057933378*coeff[0]*fSkin[7]; 
  boundSurf_incr[11] = 91.8*coeff[0]*fSkin[11]+78.09225313691493*coeff[0]*fSkin[6]+10.998181667894015*coeff[0]*fSkin[3]-19.843134832984425*coeff[0]*fSkin[2]; 
  boundSurf_incr[12] = 11.783766072743589*coeff[0]*fSkin[15]+11.61895003862225*coeff[0]*fSkin[14]+4.5*coeff[0]*fSkin[12]; 
  boundSurf_incr[13] = 91.8*coeff[0]*fSkin[13]+78.09225313691493*coeff[0]*fSkin[10]+10.998181667894016*coeff[0]*fSkin[7]-19.84313483298443*coeff[0]*fSkin[5]; 
  boundSurf_incr[14] = 45.638329755339896*coeff[0]*fSkin[15]+45.0*coeff[0]*fSkin[14]+17.428425057933378*coeff[0]*fSkin[12]; 
  boundSurf_incr[15] = 91.8*coeff[0]*fSkin[15]+78.09225313691493*coeff[0]*fSkin[14]+10.998181667894015*coeff[0]*fSkin[12]-19.84313483298443*coeff[0]*fSkin[9]; 

  }

  out[0] += -(1.0*(vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac); 
  out[1] += -(1.0*(vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac); 
  out[2] += -(1.0*(vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac); 
  out[3] += -(1.0*(vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*Jfac); 
  out[4] += -(1.0*(vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*Jfac); 
  out[5] += -(1.0*(vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*Jfac); 
  out[6] += -(1.0*(vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*Jfac); 
  out[7] += -(1.0*(vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*Jfac); 
  out[8] += -(1.0*(vol_incr[8]+edgeSurf_incr[8]+boundSurf_incr[8])*Jfac); 
  out[9] += -(1.0*(vol_incr[9]+edgeSurf_incr[9]+boundSurf_incr[9])*Jfac); 
  out[10] += -(1.0*(vol_incr[10]+edgeSurf_incr[10]+boundSurf_incr[10])*Jfac); 
  out[11] += -(1.0*(vol_incr[11]+edgeSurf_incr[11]+boundSurf_incr[11])*Jfac); 
  out[12] += -(1.0*(vol_incr[12]+edgeSurf_incr[12]+boundSurf_incr[12])*Jfac); 
  out[13] += -(1.0*(vol_incr[13]+edgeSurf_incr[13]+boundSurf_incr[13])*Jfac); 
  out[14] += -(1.0*(vol_incr[14]+edgeSurf_incr[14]+boundSurf_incr[14])*Jfac); 
  out[15] += -(1.0*(vol_incr[15]+edgeSurf_incr[15]+boundSurf_incr[15])*Jfac); 

  return 0.;
}

