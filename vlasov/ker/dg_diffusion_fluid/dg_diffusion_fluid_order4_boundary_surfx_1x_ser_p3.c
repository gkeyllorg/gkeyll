#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order4_boundary_surfx_1x_ser_p3_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],4.);

  double vol_incr[4] = {0.0}; 

  double edgeSurf_incr[4] = {0.0}; 
  double boundSurf_incr[4] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 15.347424597386395*coeff[0]*fSkin[3]+15.347424597386395*coeff[0]*fEdge[3]+26.59698043549555*coeff[0]*fSkin[2]-26.59698043549555*coeff[0]*fEdge[2]+24.864401241467277*coeff[0]*fSkin[1]+24.864401241467277*coeff[0]*fEdge[1]+14.35546875*coeff[0]*fSkin[0]-14.35546875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 31.492421698315376*coeff[0]*fSkin[3]+21.672616637695718*coeff[0]*fEdge[3]+51.3321581784444*coeff[0]*fSkin[2]-40.80248470594299*coeff[0]*fEdge[2]+45.17578125*coeff[0]*fSkin[1]+40.95703125*coeff[0]*fEdge[1]+24.864401241467277*coeff[0]*fSkin[0]-24.864401241467277*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 48.75813745345884*coeff[0]*fSkin[3]+10.726195990675592*coeff[0]*fEdge[3]+69.43359375*coeff[0]*fSkin[2]-28.65234375*coeff[0]*fEdge[2]+50.78751989538398*coeff[0]*fSkin[1]+34.44837140357144*coeff[0]*fEdge[1]+22.928431409909958*coeff[0]*fSkin[0]-22.928431409909958*coeff[0]*fEdge[0]; 
  edgeSurf_incr[3] = 51.03515625*coeff[0]*fSkin[3]-23.96484375*coeff[0]*fEdge[3]+43.053346234041356*coeff[0]*fSkin[2]+25.72108124449167*coeff[0]*fEdge[2]-6.712757365657978*coeff[0]*fSkin[1]-15.305086793700156*coeff[0]*fEdge[1]-26.199138959174753*coeff[0]*fSkin[0]+6.356004126190321*coeff[0]*fEdge[0]; 

  boundSurf_incr[1] = 11.783766072743584*coeff[0]*fSkin[3]-11.61895003862225*coeff[0]*fSkin[2]+4.5*coeff[0]*fSkin[1]; 
  boundSurf_incr[2] = -(45.638329755339896*coeff[0]*fSkin[3])+45.0*coeff[0]*fSkin[2]-17.42842505793337*coeff[0]*fSkin[1]; 
  boundSurf_incr[3] = 91.8*coeff[0]*fSkin[3]-78.09225313691493*coeff[0]*fSkin[2]+10.998181667894002*coeff[0]*fSkin[1]+19.84313483298443*coeff[0]*fSkin[0]; 

  } else { 

  edgeSurf_incr[0] = -(15.347424597386395*coeff[0]*fSkin[3])-15.347424597386395*coeff[0]*fEdge[3]+26.59698043549555*coeff[0]*fSkin[2]-26.59698043549555*coeff[0]*fEdge[2]-24.864401241467277*coeff[0]*fSkin[1]-24.864401241467277*coeff[0]*fEdge[1]+14.35546875*coeff[0]*fSkin[0]-14.35546875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 31.492421698315376*coeff[0]*fSkin[3]+21.672616637695718*coeff[0]*fEdge[3]-51.3321581784444*coeff[0]*fSkin[2]+40.80248470594299*coeff[0]*fEdge[2]+45.17578125*coeff[0]*fSkin[1]+40.95703125*coeff[0]*fEdge[1]-24.864401241467277*coeff[0]*fSkin[0]+24.864401241467277*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = -(48.75813745345884*coeff[0]*fSkin[3])-10.726195990675592*coeff[0]*fEdge[3]+69.43359375*coeff[0]*fSkin[2]-28.65234375*coeff[0]*fEdge[2]-50.78751989538398*coeff[0]*fSkin[1]-34.44837140357144*coeff[0]*fEdge[1]+22.928431409909958*coeff[0]*fSkin[0]-22.928431409909958*coeff[0]*fEdge[0]; 
  edgeSurf_incr[3] = 51.03515625*coeff[0]*fSkin[3]-23.96484375*coeff[0]*fEdge[3]-43.053346234041356*coeff[0]*fSkin[2]-25.72108124449167*coeff[0]*fEdge[2]-6.712757365657978*coeff[0]*fSkin[1]-15.305086793700156*coeff[0]*fEdge[1]+26.199138959174753*coeff[0]*fSkin[0]-6.356004126190321*coeff[0]*fEdge[0]; 

  boundSurf_incr[1] = 11.783766072743584*coeff[0]*fSkin[3]+11.61895003862225*coeff[0]*fSkin[2]+4.5*coeff[0]*fSkin[1]; 
  boundSurf_incr[2] = 45.638329755339896*coeff[0]*fSkin[3]+45.0*coeff[0]*fSkin[2]+17.42842505793337*coeff[0]*fSkin[1]; 
  boundSurf_incr[3] = 91.8*coeff[0]*fSkin[3]+78.09225313691493*coeff[0]*fSkin[2]+10.998181667894002*coeff[0]*fSkin[1]-19.84313483298443*coeff[0]*fSkin[0]; 

  }

  out[0] += -(1.0*(vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac); 
  out[1] += -(1.0*(vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac); 
  out[2] += -(1.0*(vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac); 
  out[3] += -(1.0*(vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*Jfac); 

  return 0.;
}

