#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfz_3x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[2],2.);

  double vol_incr[20] = {0.0}; 
  vol_incr[9] = 6.708203932499369*fSkin[0]*coeff[2]; 
  vol_incr[15] = 6.7082039324993685*fSkin[1]*coeff[2]; 
  vol_incr[16] = 6.7082039324993685*coeff[2]*fSkin[2]; 
  vol_incr[19] = 6.708203932499369*coeff[2]*fSkin[4]; 

  double edgeSurf_incr[20] = {0.0}; 
  double boundSurf_incr[20] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = -(0.6708203932499369*coeff[2]*fSkin[9])+0.6708203932499369*coeff[2]*fEdge[9]-1.190784930203603*coeff[2]*fSkin[3]-1.190784930203603*coeff[2]*fEdge[3]-0.9375*fSkin[0]*coeff[2]+0.9375*fEdge[0]*coeff[2]; 
  edgeSurf_incr[1] = -(0.6708203932499369*coeff[2]*fSkin[15])+0.6708203932499369*coeff[2]*fEdge[15]-1.190784930203603*coeff[2]*fSkin[5]-1.190784930203603*coeff[2]*fEdge[5]-0.9375*fSkin[1]*coeff[2]+0.9375*fEdge[1]*coeff[2]; 
  edgeSurf_incr[2] = -(0.6708203932499369*coeff[2]*fSkin[16])+0.6708203932499369*coeff[2]*fEdge[16]-1.190784930203603*coeff[2]*fSkin[6]-1.190784930203603*coeff[2]*fEdge[6]-0.9375*coeff[2]*fSkin[2]+0.9375*coeff[2]*fEdge[2]; 
  edgeSurf_incr[3] = -(1.5855025573536612*coeff[2]*fSkin[9])+0.7382874503707886*coeff[2]*fEdge[9]-2.671875*coeff[2]*fSkin[3]-1.453125*coeff[2]*fEdge[3]-2.0568103339880417*fSkin[0]*coeff[2]+1.1907849302036029*fEdge[0]*coeff[2]; 
  edgeSurf_incr[4] = -(0.6708203932499369*coeff[2]*fSkin[19])+0.6708203932499369*coeff[2]*fEdge[19]-1.190784930203603*coeff[2]*fSkin[10]-1.190784930203603*coeff[2]*fEdge[10]-0.9375*coeff[2]*fSkin[4]+0.9375*coeff[2]*fEdge[4]; 
  edgeSurf_incr[5] = -(1.5855025573536612*coeff[2]*fSkin[15])+0.7382874503707888*coeff[2]*fEdge[15]-2.671875*coeff[2]*fSkin[5]-1.453125*coeff[2]*fEdge[5]-2.0568103339880417*fSkin[1]*coeff[2]+1.1907849302036029*fEdge[1]*coeff[2]; 
  edgeSurf_incr[6] = -(1.5855025573536612*coeff[2]*fSkin[16])+0.7382874503707888*coeff[2]*fEdge[16]-2.671875*coeff[2]*fSkin[6]-1.453125*coeff[2]*fEdge[6]-2.0568103339880417*coeff[2]*fSkin[2]+1.1907849302036029*coeff[2]*fEdge[2]; 
  edgeSurf_incr[7] = -(1.190784930203603*coeff[2]*fSkin[13])-1.190784930203603*coeff[2]*fEdge[13]-0.9375*coeff[2]*fSkin[7]+0.9375*coeff[2]*fEdge[7]; 
  edgeSurf_incr[8] = -(1.190784930203603*coeff[2]*fSkin[14])-1.190784930203603*coeff[2]*fEdge[14]-0.9375*coeff[2]*fSkin[8]+0.9375*coeff[2]*fEdge[8]; 
  edgeSurf_incr[9] = -(3.140625*coeff[2]*fSkin[9])-0.140625*coeff[2]*fEdge[9]-5.022775277112744*coeff[2]*fSkin[3]-0.3025768239224549*coeff[2]*fEdge[3]-3.7733647120308955*fSkin[0]*coeff[2]+0.4192627457812108*fEdge[0]*coeff[2]; 
  edgeSurf_incr[10] = -(1.5855025573536612*coeff[2]*fSkin[19])+0.7382874503707886*coeff[2]*fEdge[19]-2.671875*coeff[2]*fSkin[10]-1.453125*coeff[2]*fEdge[10]-2.0568103339880417*coeff[2]*fSkin[4]+1.1907849302036029*coeff[2]*fEdge[4]; 
  edgeSurf_incr[11] = -(1.190784930203603*coeff[2]*fSkin[17])-1.190784930203603*coeff[2]*fEdge[17]-0.9375*coeff[2]*fSkin[11]+0.9375*coeff[2]*fEdge[11]; 
  edgeSurf_incr[12] = -(1.190784930203603*coeff[2]*fSkin[18])-1.190784930203603*coeff[2]*fEdge[18]-0.9375*coeff[2]*fSkin[12]+0.9375*coeff[2]*fEdge[12]; 
  edgeSurf_incr[13] = -(2.671875*coeff[2]*fSkin[13])-1.453125*coeff[2]*fEdge[13]-2.0568103339880417*coeff[2]*fSkin[7]+1.190784930203603*coeff[2]*fEdge[7]; 
  edgeSurf_incr[14] = -(2.671875*coeff[2]*fSkin[14])-1.453125*coeff[2]*fEdge[14]-2.0568103339880417*coeff[2]*fSkin[8]+1.190784930203603*coeff[2]*fEdge[8]; 
  edgeSurf_incr[15] = -(3.140625*coeff[2]*fSkin[15])-0.140625*coeff[2]*fEdge[15]-5.022775277112744*coeff[2]*fSkin[5]-0.30257682392245444*coeff[2]*fEdge[5]-3.773364712030894*fSkin[1]*coeff[2]+0.41926274578121053*fEdge[1]*coeff[2]; 
  edgeSurf_incr[16] = -(3.140625*coeff[2]*fSkin[16])-0.140625*coeff[2]*fEdge[16]-5.022775277112744*coeff[2]*fSkin[6]-0.30257682392245444*coeff[2]*fEdge[6]-3.773364712030894*coeff[2]*fSkin[2]+0.41926274578121053*coeff[2]*fEdge[2]; 
  edgeSurf_incr[17] = -(2.671875*coeff[2]*fSkin[17])-1.453125*coeff[2]*fEdge[17]-2.0568103339880417*coeff[2]*fSkin[11]+1.190784930203603*coeff[2]*fEdge[11]; 
  edgeSurf_incr[18] = -(2.671875*coeff[2]*fSkin[18])-1.453125*coeff[2]*fEdge[18]-2.0568103339880417*coeff[2]*fSkin[12]+1.190784930203603*coeff[2]*fEdge[12]; 
  edgeSurf_incr[19] = -(3.140625*coeff[2]*fSkin[19])-0.140625*coeff[2]*fEdge[19]-5.022775277112744*coeff[2]*fSkin[10]-0.3025768239224549*coeff[2]*fEdge[10]-3.7733647120308955*coeff[2]*fSkin[4]+0.4192627457812108*coeff[2]*fEdge[4]; 

  boundSurf_incr[3] = 0.9682458365518543*coeff[2]*fSkin[9]-1.25*coeff[2]*fSkin[3]+0.8660254037844386*fSkin[0]*coeff[2]; 
  boundSurf_incr[5] = 0.9682458365518543*coeff[2]*fSkin[15]-1.25*coeff[2]*fSkin[5]+0.8660254037844386*fSkin[1]*coeff[2]; 
  boundSurf_incr[6] = 0.9682458365518543*coeff[2]*fSkin[16]-1.25*coeff[2]*fSkin[6]+0.8660254037844386*coeff[2]*fSkin[2]; 
  boundSurf_incr[9] = -(3.75*coeff[2]*fSkin[9])+4.841229182759272*coeff[2]*fSkin[3]-3.3541019662496847*fSkin[0]*coeff[2]; 
  boundSurf_incr[10] = 0.9682458365518543*coeff[2]*fSkin[19]-1.25*coeff[2]*fSkin[10]+0.8660254037844386*coeff[2]*fSkin[4]; 
  boundSurf_incr[13] = 0.8660254037844387*coeff[2]*fSkin[7]-1.25*coeff[2]*fSkin[13]; 
  boundSurf_incr[14] = 0.8660254037844387*coeff[2]*fSkin[8]-1.25*coeff[2]*fSkin[14]; 
  boundSurf_incr[15] = -(3.75*coeff[2]*fSkin[15])+4.841229182759271*coeff[2]*fSkin[5]-3.3541019662496843*fSkin[1]*coeff[2]; 
  boundSurf_incr[16] = -(3.75*coeff[2]*fSkin[16])+4.841229182759271*coeff[2]*fSkin[6]-3.3541019662496843*coeff[2]*fSkin[2]; 
  boundSurf_incr[17] = 0.8660254037844387*coeff[2]*fSkin[11]-1.25*coeff[2]*fSkin[17]; 
  boundSurf_incr[18] = 0.8660254037844387*coeff[2]*fSkin[12]-1.25*coeff[2]*fSkin[18]; 
  boundSurf_incr[19] = -(3.75*coeff[2]*fSkin[19])+4.841229182759272*coeff[2]*fSkin[10]-3.3541019662496847*coeff[2]*fSkin[4]; 

  } else { 

  edgeSurf_incr[0] = -(0.6708203932499369*coeff[2]*fSkin[9])+0.6708203932499369*coeff[2]*fEdge[9]+1.190784930203603*coeff[2]*fSkin[3]+1.190784930203603*coeff[2]*fEdge[3]-0.9375*fSkin[0]*coeff[2]+0.9375*fEdge[0]*coeff[2]; 
  edgeSurf_incr[1] = -(0.6708203932499369*coeff[2]*fSkin[15])+0.6708203932499369*coeff[2]*fEdge[15]+1.190784930203603*coeff[2]*fSkin[5]+1.190784930203603*coeff[2]*fEdge[5]-0.9375*fSkin[1]*coeff[2]+0.9375*fEdge[1]*coeff[2]; 
  edgeSurf_incr[2] = -(0.6708203932499369*coeff[2]*fSkin[16])+0.6708203932499369*coeff[2]*fEdge[16]+1.190784930203603*coeff[2]*fSkin[6]+1.190784930203603*coeff[2]*fEdge[6]-0.9375*coeff[2]*fSkin[2]+0.9375*coeff[2]*fEdge[2]; 
  edgeSurf_incr[3] = 1.5855025573536612*coeff[2]*fSkin[9]-0.7382874503707886*coeff[2]*fEdge[9]-2.671875*coeff[2]*fSkin[3]-1.453125*coeff[2]*fEdge[3]+2.0568103339880417*fSkin[0]*coeff[2]-1.1907849302036029*fEdge[0]*coeff[2]; 
  edgeSurf_incr[4] = -(0.6708203932499369*coeff[2]*fSkin[19])+0.6708203932499369*coeff[2]*fEdge[19]+1.190784930203603*coeff[2]*fSkin[10]+1.190784930203603*coeff[2]*fEdge[10]-0.9375*coeff[2]*fSkin[4]+0.9375*coeff[2]*fEdge[4]; 
  edgeSurf_incr[5] = 1.5855025573536612*coeff[2]*fSkin[15]-0.7382874503707888*coeff[2]*fEdge[15]-2.671875*coeff[2]*fSkin[5]-1.453125*coeff[2]*fEdge[5]+2.0568103339880417*fSkin[1]*coeff[2]-1.1907849302036029*fEdge[1]*coeff[2]; 
  edgeSurf_incr[6] = 1.5855025573536612*coeff[2]*fSkin[16]-0.7382874503707888*coeff[2]*fEdge[16]-2.671875*coeff[2]*fSkin[6]-1.453125*coeff[2]*fEdge[6]+2.0568103339880417*coeff[2]*fSkin[2]-1.1907849302036029*coeff[2]*fEdge[2]; 
  edgeSurf_incr[7] = 1.190784930203603*coeff[2]*fSkin[13]+1.190784930203603*coeff[2]*fEdge[13]-0.9375*coeff[2]*fSkin[7]+0.9375*coeff[2]*fEdge[7]; 
  edgeSurf_incr[8] = 1.190784930203603*coeff[2]*fSkin[14]+1.190784930203603*coeff[2]*fEdge[14]-0.9375*coeff[2]*fSkin[8]+0.9375*coeff[2]*fEdge[8]; 
  edgeSurf_incr[9] = -(3.140625*coeff[2]*fSkin[9])-0.140625*coeff[2]*fEdge[9]+5.022775277112744*coeff[2]*fSkin[3]+0.3025768239224549*coeff[2]*fEdge[3]-3.7733647120308955*fSkin[0]*coeff[2]+0.4192627457812108*fEdge[0]*coeff[2]; 
  edgeSurf_incr[10] = 1.5855025573536612*coeff[2]*fSkin[19]-0.7382874503707886*coeff[2]*fEdge[19]-2.671875*coeff[2]*fSkin[10]-1.453125*coeff[2]*fEdge[10]+2.0568103339880417*coeff[2]*fSkin[4]-1.1907849302036029*coeff[2]*fEdge[4]; 
  edgeSurf_incr[11] = 1.190784930203603*coeff[2]*fSkin[17]+1.190784930203603*coeff[2]*fEdge[17]-0.9375*coeff[2]*fSkin[11]+0.9375*coeff[2]*fEdge[11]; 
  edgeSurf_incr[12] = 1.190784930203603*coeff[2]*fSkin[18]+1.190784930203603*coeff[2]*fEdge[18]-0.9375*coeff[2]*fSkin[12]+0.9375*coeff[2]*fEdge[12]; 
  edgeSurf_incr[13] = -(2.671875*coeff[2]*fSkin[13])-1.453125*coeff[2]*fEdge[13]+2.0568103339880417*coeff[2]*fSkin[7]-1.190784930203603*coeff[2]*fEdge[7]; 
  edgeSurf_incr[14] = -(2.671875*coeff[2]*fSkin[14])-1.453125*coeff[2]*fEdge[14]+2.0568103339880417*coeff[2]*fSkin[8]-1.190784930203603*coeff[2]*fEdge[8]; 
  edgeSurf_incr[15] = -(3.140625*coeff[2]*fSkin[15])-0.140625*coeff[2]*fEdge[15]+5.022775277112744*coeff[2]*fSkin[5]+0.30257682392245444*coeff[2]*fEdge[5]-3.773364712030894*fSkin[1]*coeff[2]+0.41926274578121053*fEdge[1]*coeff[2]; 
  edgeSurf_incr[16] = -(3.140625*coeff[2]*fSkin[16])-0.140625*coeff[2]*fEdge[16]+5.022775277112744*coeff[2]*fSkin[6]+0.30257682392245444*coeff[2]*fEdge[6]-3.773364712030894*coeff[2]*fSkin[2]+0.41926274578121053*coeff[2]*fEdge[2]; 
  edgeSurf_incr[17] = -(2.671875*coeff[2]*fSkin[17])-1.453125*coeff[2]*fEdge[17]+2.0568103339880417*coeff[2]*fSkin[11]-1.190784930203603*coeff[2]*fEdge[11]; 
  edgeSurf_incr[18] = -(2.671875*coeff[2]*fSkin[18])-1.453125*coeff[2]*fEdge[18]+2.0568103339880417*coeff[2]*fSkin[12]-1.190784930203603*coeff[2]*fEdge[12]; 
  edgeSurf_incr[19] = -(3.140625*coeff[2]*fSkin[19])-0.140625*coeff[2]*fEdge[19]+5.022775277112744*coeff[2]*fSkin[10]+0.3025768239224549*coeff[2]*fEdge[10]-3.7733647120308955*coeff[2]*fSkin[4]+0.4192627457812108*coeff[2]*fEdge[4]; 

  boundSurf_incr[3] = -(0.9682458365518543*coeff[2]*fSkin[9])-1.25*coeff[2]*fSkin[3]-0.8660254037844386*fSkin[0]*coeff[2]; 
  boundSurf_incr[5] = -(0.9682458365518543*coeff[2]*fSkin[15])-1.25*coeff[2]*fSkin[5]-0.8660254037844386*fSkin[1]*coeff[2]; 
  boundSurf_incr[6] = -(0.9682458365518543*coeff[2]*fSkin[16])-1.25*coeff[2]*fSkin[6]-0.8660254037844386*coeff[2]*fSkin[2]; 
  boundSurf_incr[9] = -(3.75*coeff[2]*fSkin[9])-4.841229182759272*coeff[2]*fSkin[3]-3.3541019662496847*fSkin[0]*coeff[2]; 
  boundSurf_incr[10] = -(0.9682458365518543*coeff[2]*fSkin[19])-1.25*coeff[2]*fSkin[10]-0.8660254037844386*coeff[2]*fSkin[4]; 
  boundSurf_incr[13] = -(1.25*coeff[2]*fSkin[13])-0.8660254037844387*coeff[2]*fSkin[7]; 
  boundSurf_incr[14] = -(1.25*coeff[2]*fSkin[14])-0.8660254037844387*coeff[2]*fSkin[8]; 
  boundSurf_incr[15] = -(3.75*coeff[2]*fSkin[15])-4.841229182759271*coeff[2]*fSkin[5]-3.3541019662496843*fSkin[1]*coeff[2]; 
  boundSurf_incr[16] = -(3.75*coeff[2]*fSkin[16])-4.841229182759271*coeff[2]*fSkin[6]-3.3541019662496843*coeff[2]*fSkin[2]; 
  boundSurf_incr[17] = -(1.25*coeff[2]*fSkin[17])-0.8660254037844387*coeff[2]*fSkin[11]; 
  boundSurf_incr[18] = -(1.25*coeff[2]*fSkin[18])-0.8660254037844387*coeff[2]*fSkin[12]; 
  boundSurf_incr[19] = -(3.75*coeff[2]*fSkin[19])-4.841229182759272*coeff[2]*fSkin[10]-3.3541019662496847*coeff[2]*fSkin[4]; 

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
  out[9] += (vol_incr[9]+edgeSurf_incr[9]+boundSurf_incr[9])*Jfac; 
  out[10] += (vol_incr[10]+edgeSurf_incr[10]+boundSurf_incr[10])*Jfac; 
  out[11] += (vol_incr[11]+edgeSurf_incr[11]+boundSurf_incr[11])*Jfac; 
  out[12] += (vol_incr[12]+edgeSurf_incr[12]+boundSurf_incr[12])*Jfac; 
  out[13] += (vol_incr[13]+edgeSurf_incr[13]+boundSurf_incr[13])*Jfac; 
  out[14] += (vol_incr[14]+edgeSurf_incr[14]+boundSurf_incr[14])*Jfac; 
  out[15] += (vol_incr[15]+edgeSurf_incr[15]+boundSurf_incr[15])*Jfac; 
  out[16] += (vol_incr[16]+edgeSurf_incr[16]+boundSurf_incr[16])*Jfac; 
  out[17] += (vol_incr[17]+edgeSurf_incr[17]+boundSurf_incr[17])*Jfac; 
  out[18] += (vol_incr[18]+edgeSurf_incr[18]+boundSurf_incr[18])*Jfac; 
  out[19] += (vol_incr[19]+edgeSurf_incr[19]+boundSurf_incr[19])*Jfac; 

  return 0.;
}

