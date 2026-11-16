#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order4_boundary_surfy_2x_ser_p3_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[1],4.);

  double vol_incr[12] = {0.0}; 

  double edgeSurf_incr[12] = {0.0}; 
  double boundSurf_incr[12] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 15.347424597386395*coeff[1]*fSkin[9]+15.347424597386395*coeff[1]*fEdge[9]+26.59698043549555*coeff[1]*fSkin[5]-26.59698043549555*coeff[1]*fEdge[5]+24.864401241467277*coeff[1]*fSkin[2]+24.864401241467277*coeff[1]*fEdge[2]+14.35546875*fSkin[0]*coeff[1]-14.35546875*fEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = 15.347424597386388*coeff[1]*fSkin[11]+15.347424597386388*coeff[1]*fEdge[11]+26.596980435495546*coeff[1]*fSkin[7]-26.596980435495546*coeff[1]*fEdge[7]+24.864401241467277*coeff[1]*fSkin[3]+24.864401241467277*coeff[1]*fEdge[3]+14.35546875*coeff[1]*fSkin[1]-14.35546875*coeff[1]*fEdge[1]; 
  edgeSurf_incr[2] = 31.492421698315376*coeff[1]*fSkin[9]+21.672616637695718*coeff[1]*fEdge[9]+51.3321581784444*coeff[1]*fSkin[5]-40.80248470594299*coeff[1]*fEdge[5]+45.17578125*coeff[1]*fSkin[2]+40.95703125*coeff[1]*fEdge[2]+24.864401241467277*fSkin[0]*coeff[1]-24.864401241467277*fEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = 31.492421698315383*coeff[1]*fSkin[11]+21.67261663769573*coeff[1]*fEdge[11]+51.3321581784444*coeff[1]*fSkin[7]-40.80248470594299*coeff[1]*fEdge[7]+45.17578125*coeff[1]*fSkin[3]+40.95703125*coeff[1]*fEdge[3]+24.864401241467277*coeff[1]*fSkin[1]-24.864401241467277*coeff[1]*fEdge[1]; 
  edgeSurf_incr[4] = 24.864401241467288*coeff[1]*fSkin[6]+24.864401241467288*coeff[1]*fEdge[6]+14.35546875*coeff[1]*fSkin[4]-14.35546875*coeff[1]*fEdge[4]; 
  edgeSurf_incr[5] = 48.75813745345884*coeff[1]*fSkin[9]+10.726195990675592*coeff[1]*fEdge[9]+69.43359375*coeff[1]*fSkin[5]-28.65234375*coeff[1]*fEdge[5]+50.78751989538398*coeff[1]*fSkin[2]+34.44837140357144*coeff[1]*fEdge[2]+22.928431409909958*fSkin[0]*coeff[1]-22.928431409909958*fEdge[0]*coeff[1]; 
  edgeSurf_incr[6] = 45.17578125*coeff[1]*fSkin[6]+40.95703125*coeff[1]*fEdge[6]+24.864401241467288*coeff[1]*fSkin[4]-24.864401241467288*coeff[1]*fEdge[4]; 
  edgeSurf_incr[7] = 48.758137453458836*coeff[1]*fSkin[11]+10.726195990675588*coeff[1]*fEdge[11]+69.43359375*coeff[1]*fSkin[7]-28.65234375*coeff[1]*fEdge[7]+50.78751989538398*coeff[1]*fSkin[3]+34.44837140357144*coeff[1]*fEdge[3]+22.92843140990995*coeff[1]*fSkin[1]-22.92843140990995*coeff[1]*fEdge[1]; 
  edgeSurf_incr[8] = 24.86440124146728*coeff[1]*fSkin[10]+24.86440124146728*coeff[1]*fEdge[10]+14.35546875*coeff[1]*fSkin[8]-14.35546875*coeff[1]*fEdge[8]; 
  edgeSurf_incr[9] = 51.03515625*coeff[1]*fSkin[9]-23.96484375*coeff[1]*fEdge[9]+43.053346234041356*coeff[1]*fSkin[5]+25.72108124449167*coeff[1]*fEdge[5]-6.712757365657978*coeff[1]*fSkin[2]-15.305086793700156*coeff[1]*fEdge[2]-26.199138959174753*fSkin[0]*coeff[1]+6.356004126190321*fEdge[0]*coeff[1]; 
  edgeSurf_incr[10] = 45.17578125*coeff[1]*fSkin[10]+40.95703125*coeff[1]*fEdge[10]+24.86440124146728*coeff[1]*fSkin[8]-24.86440124146728*coeff[1]*fEdge[8]; 
  edgeSurf_incr[11] = 51.03515625*coeff[1]*fSkin[11]-23.96484375*coeff[1]*fEdge[11]+43.05334623404135*coeff[1]*fSkin[7]+25.721081244491685*coeff[1]*fEdge[7]-6.712757365657971*coeff[1]*fSkin[3]-15.30508679370017*coeff[1]*fEdge[3]-26.19913895917474*coeff[1]*fSkin[1]+6.356004126190324*coeff[1]*fEdge[1]; 

  boundSurf_incr[2] = 11.783766072743584*coeff[1]*fSkin[9]-11.61895003862225*coeff[1]*fSkin[5]+4.5*coeff[1]*fSkin[2]; 
  boundSurf_incr[3] = 11.783766072743589*coeff[1]*fSkin[11]-11.618950038622252*coeff[1]*fSkin[7]+4.5*coeff[1]*fSkin[3]; 
  boundSurf_incr[5] = -(45.638329755339896*coeff[1]*fSkin[9])+45.0*coeff[1]*fSkin[5]-17.42842505793337*coeff[1]*fSkin[2]; 
  boundSurf_incr[6] = 4.5*coeff[1]*fSkin[6]; 
  boundSurf_incr[7] = -(45.6383297553399*coeff[1]*fSkin[11])+45.0*coeff[1]*fSkin[7]-17.428425057933378*coeff[1]*fSkin[3]; 
  boundSurf_incr[9] = 91.8*coeff[1]*fSkin[9]-78.09225313691493*coeff[1]*fSkin[5]+10.998181667894002*coeff[1]*fSkin[2]+19.84313483298443*fSkin[0]*coeff[1]; 
  boundSurf_incr[10] = 4.5*coeff[1]*fSkin[10]; 
  boundSurf_incr[11] = 91.8*coeff[1]*fSkin[11]-78.09225313691493*coeff[1]*fSkin[7]+10.998181667894015*coeff[1]*fSkin[3]+19.843134832984425*coeff[1]*fSkin[1]; 

  } else { 

  edgeSurf_incr[0] = -(15.347424597386395*coeff[1]*fSkin[9])-15.347424597386395*coeff[1]*fEdge[9]+26.59698043549555*coeff[1]*fSkin[5]-26.59698043549555*coeff[1]*fEdge[5]-24.864401241467277*coeff[1]*fSkin[2]-24.864401241467277*coeff[1]*fEdge[2]+14.35546875*fSkin[0]*coeff[1]-14.35546875*fEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = -(15.347424597386388*coeff[1]*fSkin[11])-15.347424597386388*coeff[1]*fEdge[11]+26.596980435495546*coeff[1]*fSkin[7]-26.596980435495546*coeff[1]*fEdge[7]-24.864401241467277*coeff[1]*fSkin[3]-24.864401241467277*coeff[1]*fEdge[3]+14.35546875*coeff[1]*fSkin[1]-14.35546875*coeff[1]*fEdge[1]; 
  edgeSurf_incr[2] = 31.492421698315376*coeff[1]*fSkin[9]+21.672616637695718*coeff[1]*fEdge[9]-51.3321581784444*coeff[1]*fSkin[5]+40.80248470594299*coeff[1]*fEdge[5]+45.17578125*coeff[1]*fSkin[2]+40.95703125*coeff[1]*fEdge[2]-24.864401241467277*fSkin[0]*coeff[1]+24.864401241467277*fEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = 31.492421698315383*coeff[1]*fSkin[11]+21.67261663769573*coeff[1]*fEdge[11]-51.3321581784444*coeff[1]*fSkin[7]+40.80248470594299*coeff[1]*fEdge[7]+45.17578125*coeff[1]*fSkin[3]+40.95703125*coeff[1]*fEdge[3]-24.864401241467277*coeff[1]*fSkin[1]+24.864401241467277*coeff[1]*fEdge[1]; 
  edgeSurf_incr[4] = -(24.864401241467288*coeff[1]*fSkin[6])-24.864401241467288*coeff[1]*fEdge[6]+14.35546875*coeff[1]*fSkin[4]-14.35546875*coeff[1]*fEdge[4]; 
  edgeSurf_incr[5] = -(48.75813745345884*coeff[1]*fSkin[9])-10.726195990675592*coeff[1]*fEdge[9]+69.43359375*coeff[1]*fSkin[5]-28.65234375*coeff[1]*fEdge[5]-50.78751989538398*coeff[1]*fSkin[2]-34.44837140357144*coeff[1]*fEdge[2]+22.928431409909958*fSkin[0]*coeff[1]-22.928431409909958*fEdge[0]*coeff[1]; 
  edgeSurf_incr[6] = 45.17578125*coeff[1]*fSkin[6]+40.95703125*coeff[1]*fEdge[6]-24.864401241467288*coeff[1]*fSkin[4]+24.864401241467288*coeff[1]*fEdge[4]; 
  edgeSurf_incr[7] = -(48.758137453458836*coeff[1]*fSkin[11])-10.726195990675588*coeff[1]*fEdge[11]+69.43359375*coeff[1]*fSkin[7]-28.65234375*coeff[1]*fEdge[7]-50.78751989538398*coeff[1]*fSkin[3]-34.44837140357144*coeff[1]*fEdge[3]+22.92843140990995*coeff[1]*fSkin[1]-22.92843140990995*coeff[1]*fEdge[1]; 
  edgeSurf_incr[8] = -(24.86440124146728*coeff[1]*fSkin[10])-24.86440124146728*coeff[1]*fEdge[10]+14.35546875*coeff[1]*fSkin[8]-14.35546875*coeff[1]*fEdge[8]; 
  edgeSurf_incr[9] = 51.03515625*coeff[1]*fSkin[9]-23.96484375*coeff[1]*fEdge[9]-43.053346234041356*coeff[1]*fSkin[5]-25.72108124449167*coeff[1]*fEdge[5]-6.712757365657978*coeff[1]*fSkin[2]-15.305086793700156*coeff[1]*fEdge[2]+26.199138959174753*fSkin[0]*coeff[1]-6.356004126190321*fEdge[0]*coeff[1]; 
  edgeSurf_incr[10] = 45.17578125*coeff[1]*fSkin[10]+40.95703125*coeff[1]*fEdge[10]-24.86440124146728*coeff[1]*fSkin[8]+24.86440124146728*coeff[1]*fEdge[8]; 
  edgeSurf_incr[11] = 51.03515625*coeff[1]*fSkin[11]-23.96484375*coeff[1]*fEdge[11]-43.05334623404135*coeff[1]*fSkin[7]-25.721081244491685*coeff[1]*fEdge[7]-6.712757365657971*coeff[1]*fSkin[3]-15.30508679370017*coeff[1]*fEdge[3]+26.19913895917474*coeff[1]*fSkin[1]-6.356004126190324*coeff[1]*fEdge[1]; 

  boundSurf_incr[2] = 11.783766072743584*coeff[1]*fSkin[9]+11.61895003862225*coeff[1]*fSkin[5]+4.5*coeff[1]*fSkin[2]; 
  boundSurf_incr[3] = 11.783766072743589*coeff[1]*fSkin[11]+11.618950038622252*coeff[1]*fSkin[7]+4.5*coeff[1]*fSkin[3]; 
  boundSurf_incr[5] = 45.638329755339896*coeff[1]*fSkin[9]+45.0*coeff[1]*fSkin[5]+17.42842505793337*coeff[1]*fSkin[2]; 
  boundSurf_incr[6] = 4.5*coeff[1]*fSkin[6]; 
  boundSurf_incr[7] = 45.6383297553399*coeff[1]*fSkin[11]+45.0*coeff[1]*fSkin[7]+17.428425057933378*coeff[1]*fSkin[3]; 
  boundSurf_incr[9] = 91.8*coeff[1]*fSkin[9]+78.09225313691493*coeff[1]*fSkin[5]+10.998181667894002*coeff[1]*fSkin[2]-19.84313483298443*fSkin[0]*coeff[1]; 
  boundSurf_incr[10] = 4.5*coeff[1]*fSkin[10]; 
  boundSurf_incr[11] = 91.8*coeff[1]*fSkin[11]+78.09225313691493*coeff[1]*fSkin[7]+10.998181667894015*coeff[1]*fSkin[3]-19.843134832984425*coeff[1]*fSkin[1]; 

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

  return 0.;
}

