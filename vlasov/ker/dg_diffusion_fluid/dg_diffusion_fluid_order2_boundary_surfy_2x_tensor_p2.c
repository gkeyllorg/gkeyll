#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfy_2x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[1],2.);

  double vol_incr[9] = {0.0}; 
  vol_incr[5] = 6.708203932499369*fSkin[0]*coeff[1]; 
  vol_incr[7] = 6.7082039324993685*coeff[1]*fSkin[1]; 
  vol_incr[8] = 6.708203932499369*coeff[1]*fSkin[4]; 

  double edgeSurf_incr[9] = {0.0}; 
  double boundSurf_incr[9] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = -(0.6708203932499369*coeff[1]*fSkin[5])+0.6708203932499369*coeff[1]*fEdge[5]-1.190784930203603*coeff[1]*fSkin[2]-1.190784930203603*coeff[1]*fEdge[2]-0.9375*fSkin[0]*coeff[1]+0.9375*fEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = -(0.6708203932499369*coeff[1]*fSkin[7])+0.6708203932499369*coeff[1]*fEdge[7]-1.190784930203603*coeff[1]*fSkin[3]-1.190784930203603*coeff[1]*fEdge[3]-0.9375*coeff[1]*fSkin[1]+0.9375*coeff[1]*fEdge[1]; 
  edgeSurf_incr[2] = -(1.5855025573536612*coeff[1]*fSkin[5])+0.7382874503707886*coeff[1]*fEdge[5]-2.671875*coeff[1]*fSkin[2]-1.453125*coeff[1]*fEdge[2]-2.0568103339880417*fSkin[0]*coeff[1]+1.1907849302036029*fEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = -(1.5855025573536612*coeff[1]*fSkin[7])+0.7382874503707888*coeff[1]*fEdge[7]-2.671875*coeff[1]*fSkin[3]-1.453125*coeff[1]*fEdge[3]-2.0568103339880417*coeff[1]*fSkin[1]+1.1907849302036029*coeff[1]*fEdge[1]; 
  edgeSurf_incr[4] = -(0.6708203932499369*coeff[1]*fSkin[8])+0.6708203932499369*coeff[1]*fEdge[8]-1.190784930203603*coeff[1]*fSkin[6]-1.190784930203603*coeff[1]*fEdge[6]-0.9375*coeff[1]*fSkin[4]+0.9375*coeff[1]*fEdge[4]; 
  edgeSurf_incr[5] = -(3.140625*coeff[1]*fSkin[5])-0.140625*coeff[1]*fEdge[5]-5.022775277112744*coeff[1]*fSkin[2]-0.3025768239224549*coeff[1]*fEdge[2]-3.7733647120308955*fSkin[0]*coeff[1]+0.4192627457812108*fEdge[0]*coeff[1]; 
  edgeSurf_incr[6] = -(1.5855025573536612*coeff[1]*fSkin[8])+0.7382874503707888*coeff[1]*fEdge[8]-2.671875*coeff[1]*fSkin[6]-1.453125*coeff[1]*fEdge[6]-2.0568103339880417*coeff[1]*fSkin[4]+1.190784930203603*coeff[1]*fEdge[4]; 
  edgeSurf_incr[7] = -(3.140625*coeff[1]*fSkin[7])-0.140625*coeff[1]*fEdge[7]-5.022775277112744*coeff[1]*fSkin[3]-0.30257682392245444*coeff[1]*fEdge[3]-3.773364712030894*coeff[1]*fSkin[1]+0.41926274578121053*coeff[1]*fEdge[1]; 
  edgeSurf_incr[8] = -(3.140625*coeff[1]*fSkin[8])-0.140625*coeff[1]*fEdge[8]-5.022775277112744*coeff[1]*fSkin[6]-0.30257682392245444*coeff[1]*fEdge[6]-3.7733647120308955*coeff[1]*fSkin[4]+0.4192627457812108*coeff[1]*fEdge[4]; 

  boundSurf_incr[2] = 0.9682458365518543*coeff[1]*fSkin[5]-1.25*coeff[1]*fSkin[2]+0.8660254037844386*fSkin[0]*coeff[1]; 
  boundSurf_incr[3] = 0.9682458365518543*coeff[1]*fSkin[7]-1.25*coeff[1]*fSkin[3]+0.8660254037844386*coeff[1]*fSkin[1]; 
  boundSurf_incr[5] = -(3.75*coeff[1]*fSkin[5])+4.841229182759272*coeff[1]*fSkin[2]-3.3541019662496847*fSkin[0]*coeff[1]; 
  boundSurf_incr[6] = 0.9682458365518543*coeff[1]*fSkin[8]-1.25*coeff[1]*fSkin[6]+0.8660254037844387*coeff[1]*fSkin[4]; 
  boundSurf_incr[7] = -(3.75*coeff[1]*fSkin[7])+4.841229182759271*coeff[1]*fSkin[3]-3.3541019662496843*coeff[1]*fSkin[1]; 
  boundSurf_incr[8] = -(3.75*coeff[1]*fSkin[8])+4.841229182759271*coeff[1]*fSkin[6]-3.3541019662496847*coeff[1]*fSkin[4]; 

  } else { 

  edgeSurf_incr[0] = -(0.6708203932499369*coeff[1]*fSkin[5])+0.6708203932499369*coeff[1]*fEdge[5]+1.190784930203603*coeff[1]*fSkin[2]+1.190784930203603*coeff[1]*fEdge[2]-0.9375*fSkin[0]*coeff[1]+0.9375*fEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = -(0.6708203932499369*coeff[1]*fSkin[7])+0.6708203932499369*coeff[1]*fEdge[7]+1.190784930203603*coeff[1]*fSkin[3]+1.190784930203603*coeff[1]*fEdge[3]-0.9375*coeff[1]*fSkin[1]+0.9375*coeff[1]*fEdge[1]; 
  edgeSurf_incr[2] = 1.5855025573536612*coeff[1]*fSkin[5]-0.7382874503707886*coeff[1]*fEdge[5]-2.671875*coeff[1]*fSkin[2]-1.453125*coeff[1]*fEdge[2]+2.0568103339880417*fSkin[0]*coeff[1]-1.1907849302036029*fEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = 1.5855025573536612*coeff[1]*fSkin[7]-0.7382874503707888*coeff[1]*fEdge[7]-2.671875*coeff[1]*fSkin[3]-1.453125*coeff[1]*fEdge[3]+2.0568103339880417*coeff[1]*fSkin[1]-1.1907849302036029*coeff[1]*fEdge[1]; 
  edgeSurf_incr[4] = -(0.6708203932499369*coeff[1]*fSkin[8])+0.6708203932499369*coeff[1]*fEdge[8]+1.190784930203603*coeff[1]*fSkin[6]+1.190784930203603*coeff[1]*fEdge[6]-0.9375*coeff[1]*fSkin[4]+0.9375*coeff[1]*fEdge[4]; 
  edgeSurf_incr[5] = -(3.140625*coeff[1]*fSkin[5])-0.140625*coeff[1]*fEdge[5]+5.022775277112744*coeff[1]*fSkin[2]+0.3025768239224549*coeff[1]*fEdge[2]-3.7733647120308955*fSkin[0]*coeff[1]+0.4192627457812108*fEdge[0]*coeff[1]; 
  edgeSurf_incr[6] = 1.5855025573536612*coeff[1]*fSkin[8]-0.7382874503707888*coeff[1]*fEdge[8]-2.671875*coeff[1]*fSkin[6]-1.453125*coeff[1]*fEdge[6]+2.0568103339880417*coeff[1]*fSkin[4]-1.190784930203603*coeff[1]*fEdge[4]; 
  edgeSurf_incr[7] = -(3.140625*coeff[1]*fSkin[7])-0.140625*coeff[1]*fEdge[7]+5.022775277112744*coeff[1]*fSkin[3]+0.30257682392245444*coeff[1]*fEdge[3]-3.773364712030894*coeff[1]*fSkin[1]+0.41926274578121053*coeff[1]*fEdge[1]; 
  edgeSurf_incr[8] = -(3.140625*coeff[1]*fSkin[8])-0.140625*coeff[1]*fEdge[8]+5.022775277112744*coeff[1]*fSkin[6]+0.30257682392245444*coeff[1]*fEdge[6]-3.7733647120308955*coeff[1]*fSkin[4]+0.4192627457812108*coeff[1]*fEdge[4]; 

  boundSurf_incr[2] = -(0.9682458365518543*coeff[1]*fSkin[5])-1.25*coeff[1]*fSkin[2]-0.8660254037844386*fSkin[0]*coeff[1]; 
  boundSurf_incr[3] = -(0.9682458365518543*coeff[1]*fSkin[7])-1.25*coeff[1]*fSkin[3]-0.8660254037844386*coeff[1]*fSkin[1]; 
  boundSurf_incr[5] = -(3.75*coeff[1]*fSkin[5])-4.841229182759272*coeff[1]*fSkin[2]-3.3541019662496847*fSkin[0]*coeff[1]; 
  boundSurf_incr[6] = -(0.9682458365518543*coeff[1]*fSkin[8])-1.25*coeff[1]*fSkin[6]-0.8660254037844387*coeff[1]*fSkin[4]; 
  boundSurf_incr[7] = -(3.75*coeff[1]*fSkin[7])-4.841229182759271*coeff[1]*fSkin[3]-3.3541019662496843*coeff[1]*fSkin[1]; 
  boundSurf_incr[8] = -(3.75*coeff[1]*fSkin[8])-4.841229182759271*coeff[1]*fSkin[6]-3.3541019662496847*coeff[1]*fSkin[4]; 

  }

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*Jfac; 
  out[4] += (vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*Jfac; 
  out[5] += (vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*Jfac; 
  out[6] += (vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*Jfac; 
  out[7] += (vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*Jfac; 
  out[8] += (vol_incr[8]+edgeSurf_incr[8]+boundSurf_incr[8])*Jfac; 

  return 0.;
}

