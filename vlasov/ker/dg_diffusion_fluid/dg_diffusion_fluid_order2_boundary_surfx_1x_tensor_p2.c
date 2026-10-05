#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfx_1x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],2.);

  double vol_incr[3] = {0.0}; 
  vol_incr[2] = 6.708203932499369*coeff[0]*fSkin[0]; 

  double edgeSurf_incr[3] = {0.0}; 
  double boundSurf_incr[3] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = -(0.6708203932499369*coeff[0]*fSkin[2])+0.6708203932499369*coeff[0]*fEdge[2]-1.190784930203603*coeff[0]*fSkin[1]-1.190784930203603*coeff[0]*fEdge[1]-0.9375*coeff[0]*fSkin[0]+0.9375*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = -(1.5855025573536612*coeff[0]*fSkin[2])+0.7382874503707886*coeff[0]*fEdge[2]-2.671875*coeff[0]*fSkin[1]-1.453125*coeff[0]*fEdge[1]-2.0568103339880417*coeff[0]*fSkin[0]+1.1907849302036029*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = -(3.140625*coeff[0]*fSkin[2])-0.140625*coeff[0]*fEdge[2]-5.022775277112744*coeff[0]*fSkin[1]-0.3025768239224549*coeff[0]*fEdge[1]-3.7733647120308955*coeff[0]*fSkin[0]+0.4192627457812108*coeff[0]*fEdge[0]; 

  boundSurf_incr[1] = 0.9682458365518543*coeff[0]*fSkin[2]-1.25*coeff[0]*fSkin[1]+0.8660254037844386*coeff[0]*fSkin[0]; 
  boundSurf_incr[2] = -(3.75*coeff[0]*fSkin[2])+4.841229182759272*coeff[0]*fSkin[1]-3.3541019662496847*coeff[0]*fSkin[0]; 

  } else { 

  edgeSurf_incr[0] = -(0.6708203932499369*coeff[0]*fSkin[2])+0.6708203932499369*coeff[0]*fEdge[2]+1.190784930203603*coeff[0]*fSkin[1]+1.190784930203603*coeff[0]*fEdge[1]-0.9375*coeff[0]*fSkin[0]+0.9375*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 1.5855025573536612*coeff[0]*fSkin[2]-0.7382874503707886*coeff[0]*fEdge[2]-2.671875*coeff[0]*fSkin[1]-1.453125*coeff[0]*fEdge[1]+2.0568103339880417*coeff[0]*fSkin[0]-1.1907849302036029*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = -(3.140625*coeff[0]*fSkin[2])-0.140625*coeff[0]*fEdge[2]+5.022775277112744*coeff[0]*fSkin[1]+0.3025768239224549*coeff[0]*fEdge[1]-3.7733647120308955*coeff[0]*fSkin[0]+0.4192627457812108*coeff[0]*fEdge[0]; 

  boundSurf_incr[1] = -(0.9682458365518543*coeff[0]*fSkin[2])-1.25*coeff[0]*fSkin[1]-0.8660254037844386*coeff[0]*fSkin[0]; 
  boundSurf_incr[2] = -(3.75*coeff[0]*fSkin[2])-4.841229182759272*coeff[0]*fSkin[1]-3.3541019662496847*coeff[0]*fSkin[0]; 

  }

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 

  return 0.;
}

