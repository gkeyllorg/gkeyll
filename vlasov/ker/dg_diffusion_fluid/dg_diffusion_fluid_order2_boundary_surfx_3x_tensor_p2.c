#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfx_3x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],2.);

  double vol_incr[27] = {0.0}; 
  vol_incr[7] = 6.708203932499369*coeff[0]*fSkin[0]; 
  vol_incr[11] = 6.7082039324993685*coeff[0]*fSkin[2]; 
  vol_incr[13] = 6.7082039324993685*coeff[0]*fSkin[3]; 
  vol_incr[17] = 6.708203932499369*coeff[0]*fSkin[6]; 
  vol_incr[20] = 6.708203932499369*coeff[0]*fSkin[8]; 
  vol_incr[21] = 6.708203932499369*coeff[0]*fSkin[9]; 
  vol_incr[23] = 6.7082039324993685*coeff[0]*fSkin[14]; 
  vol_incr[24] = 6.7082039324993685*coeff[0]*fSkin[16]; 
  vol_incr[26] = 6.708203932499369*coeff[0]*fSkin[22]; 

  double edgeSurf_incr[27] = {0.0}; 
  double boundSurf_incr[27] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = -(0.6708203932499369*coeff[0]*fSkin[7])+0.6708203932499369*coeff[0]*fEdge[7]-1.190784930203603*coeff[0]*fSkin[1]-1.190784930203603*coeff[0]*fEdge[1]-0.9375*coeff[0]*fSkin[0]+0.9375*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = -(1.5855025573536612*coeff[0]*fSkin[7])+0.7382874503707886*coeff[0]*fEdge[7]-2.671875*coeff[0]*fSkin[1]-1.453125*coeff[0]*fEdge[1]-2.0568103339880417*coeff[0]*fSkin[0]+1.1907849302036029*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = -(0.6708203932499369*coeff[0]*fSkin[11])+0.6708203932499369*coeff[0]*fEdge[11]-1.190784930203603*coeff[0]*fSkin[4]-1.190784930203603*coeff[0]*fEdge[4]-0.9375*coeff[0]*fSkin[2]+0.9375*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = -(0.6708203932499369*coeff[0]*fSkin[13])+0.6708203932499369*coeff[0]*fEdge[13]-1.190784930203603*coeff[0]*fSkin[5]-1.190784930203603*coeff[0]*fEdge[5]-0.9375*coeff[0]*fSkin[3]+0.9375*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = -(1.5855025573536612*coeff[0]*fSkin[11])+0.7382874503707888*coeff[0]*fEdge[11]-2.671875*coeff[0]*fSkin[4]-1.453125*coeff[0]*fEdge[4]-2.0568103339880417*coeff[0]*fSkin[2]+1.1907849302036029*coeff[0]*fEdge[2]; 
  edgeSurf_incr[5] = -(1.5855025573536612*coeff[0]*fSkin[13])+0.7382874503707888*coeff[0]*fEdge[13]-2.671875*coeff[0]*fSkin[5]-1.453125*coeff[0]*fEdge[5]-2.0568103339880417*coeff[0]*fSkin[3]+1.1907849302036029*coeff[0]*fEdge[3]; 
  edgeSurf_incr[6] = -(0.6708203932499369*coeff[0]*fSkin[17])+0.6708203932499369*coeff[0]*fEdge[17]-1.190784930203603*coeff[0]*fSkin[10]-1.190784930203603*coeff[0]*fEdge[10]-0.9375*coeff[0]*fSkin[6]+0.9375*coeff[0]*fEdge[6]; 
  edgeSurf_incr[7] = -(3.140625*coeff[0]*fSkin[7])-0.140625*coeff[0]*fEdge[7]-5.022775277112744*coeff[0]*fSkin[1]-0.3025768239224549*coeff[0]*fEdge[1]-3.7733647120308955*coeff[0]*fSkin[0]+0.4192627457812108*coeff[0]*fEdge[0]; 
  edgeSurf_incr[8] = -(0.6708203932499369*coeff[0]*fSkin[20])+0.6708203932499369*coeff[0]*fEdge[20]-1.190784930203603*coeff[0]*fSkin[12]-1.190784930203603*coeff[0]*fEdge[12]-0.9375*coeff[0]*fSkin[8]+0.9375*coeff[0]*fEdge[8]; 
  edgeSurf_incr[9] = -(0.6708203932499369*coeff[0]*fSkin[21])+0.6708203932499369*coeff[0]*fEdge[21]-1.190784930203603*coeff[0]*fSkin[15]-1.190784930203603*coeff[0]*fEdge[15]-0.9375*coeff[0]*fSkin[9]+0.9375*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = -(1.5855025573536612*coeff[0]*fSkin[17])+0.7382874503707886*coeff[0]*fEdge[17]-2.671875*coeff[0]*fSkin[10]-1.453125*coeff[0]*fEdge[10]-2.0568103339880417*coeff[0]*fSkin[6]+1.1907849302036029*coeff[0]*fEdge[6]; 
  edgeSurf_incr[11] = -(3.140625*coeff[0]*fSkin[11])-0.140625*coeff[0]*fEdge[11]-5.022775277112744*coeff[0]*fSkin[4]-0.30257682392245444*coeff[0]*fEdge[4]-3.773364712030894*coeff[0]*fSkin[2]+0.41926274578121053*coeff[0]*fEdge[2]; 
  edgeSurf_incr[12] = -(1.5855025573536612*coeff[0]*fSkin[20])+0.7382874503707888*coeff[0]*fEdge[20]-2.671875*coeff[0]*fSkin[12]-1.453125*coeff[0]*fEdge[12]-2.0568103339880417*coeff[0]*fSkin[8]+1.190784930203603*coeff[0]*fEdge[8]; 
  edgeSurf_incr[13] = -(3.140625*coeff[0]*fSkin[13])-0.140625*coeff[0]*fEdge[13]-5.022775277112744*coeff[0]*fSkin[5]-0.30257682392245444*coeff[0]*fEdge[5]-3.773364712030894*coeff[0]*fSkin[3]+0.41926274578121053*coeff[0]*fEdge[3]; 
  edgeSurf_incr[14] = -(0.6708203932499369*coeff[0]*fSkin[23])+0.6708203932499369*coeff[0]*fEdge[23]-1.190784930203603*coeff[0]*fSkin[18]-1.190784930203603*coeff[0]*fEdge[18]-0.9375*coeff[0]*fSkin[14]+0.9375*coeff[0]*fEdge[14]; 
  edgeSurf_incr[15] = -(1.5855025573536612*coeff[0]*fSkin[21])+0.7382874503707888*coeff[0]*fEdge[21]-2.671875*coeff[0]*fSkin[15]-1.453125*coeff[0]*fEdge[15]-2.0568103339880417*coeff[0]*fSkin[9]+1.190784930203603*coeff[0]*fEdge[9]; 
  edgeSurf_incr[16] = -(0.6708203932499369*coeff[0]*fSkin[24])+0.6708203932499369*coeff[0]*fEdge[24]-1.190784930203603*coeff[0]*fSkin[19]-1.190784930203603*coeff[0]*fEdge[19]-0.9375*coeff[0]*fSkin[16]+0.9375*coeff[0]*fEdge[16]; 
  edgeSurf_incr[17] = -(3.140625*coeff[0]*fSkin[17])-0.140625*coeff[0]*fEdge[17]-5.022775277112744*coeff[0]*fSkin[10]-0.3025768239224549*coeff[0]*fEdge[10]-3.7733647120308955*coeff[0]*fSkin[6]+0.4192627457812108*coeff[0]*fEdge[6]; 
  edgeSurf_incr[18] = -(1.5855025573536612*coeff[0]*fSkin[23])+0.7382874503707886*coeff[0]*fEdge[23]-2.671875*coeff[0]*fSkin[18]-1.453125*coeff[0]*fEdge[18]-2.0568103339880417*coeff[0]*fSkin[14]+1.190784930203603*coeff[0]*fEdge[14]; 
  edgeSurf_incr[19] = -(1.5855025573536612*coeff[0]*fSkin[24])+0.7382874503707886*coeff[0]*fEdge[24]-2.671875*coeff[0]*fSkin[19]-1.453125*coeff[0]*fEdge[19]-2.0568103339880417*coeff[0]*fSkin[16]+1.190784930203603*coeff[0]*fEdge[16]; 
  edgeSurf_incr[20] = -(3.140625*coeff[0]*fSkin[20])-0.140625*coeff[0]*fEdge[20]-5.022775277112744*coeff[0]*fSkin[12]-0.30257682392245444*coeff[0]*fEdge[12]-3.7733647120308955*coeff[0]*fSkin[8]+0.4192627457812108*coeff[0]*fEdge[8]; 
  edgeSurf_incr[21] = -(3.140625*coeff[0]*fSkin[21])-0.140625*coeff[0]*fEdge[21]-5.022775277112744*coeff[0]*fSkin[15]-0.30257682392245444*coeff[0]*fEdge[15]-3.7733647120308955*coeff[0]*fSkin[9]+0.4192627457812108*coeff[0]*fEdge[9]; 
  edgeSurf_incr[22] = -(0.6708203932499369*coeff[0]*fSkin[26])+0.6708203932499369*coeff[0]*fEdge[26]-1.190784930203603*coeff[0]*fSkin[25]-1.190784930203603*coeff[0]*fEdge[25]-0.9375*coeff[0]*fSkin[22]+0.9375*coeff[0]*fEdge[22]; 
  edgeSurf_incr[23] = -(3.140625*coeff[0]*fSkin[23])-0.140625*coeff[0]*fEdge[23]-5.022775277112744*coeff[0]*fSkin[18]-0.3025768239224549*coeff[0]*fEdge[18]-3.773364712030894*coeff[0]*fSkin[14]+0.41926274578121053*coeff[0]*fEdge[14]; 
  edgeSurf_incr[24] = -(3.140625*coeff[0]*fSkin[24])-0.140625*coeff[0]*fEdge[24]-5.022775277112744*coeff[0]*fSkin[19]-0.3025768239224549*coeff[0]*fEdge[19]-3.773364712030894*coeff[0]*fSkin[16]+0.41926274578121053*coeff[0]*fEdge[16]; 
  edgeSurf_incr[25] = -(1.5855025573536612*coeff[0]*fSkin[26])+0.7382874503707886*coeff[0]*fEdge[26]-2.671875*coeff[0]*fSkin[25]-1.453125*coeff[0]*fEdge[25]-2.0568103339880417*coeff[0]*fSkin[22]+1.1907849302036029*coeff[0]*fEdge[22]; 
  edgeSurf_incr[26] = -(3.140625*coeff[0]*fSkin[26])-0.140625*coeff[0]*fEdge[26]-5.022775277112744*coeff[0]*fSkin[25]-0.3025768239224549*coeff[0]*fEdge[25]-3.7733647120308955*coeff[0]*fSkin[22]+0.4192627457812108*coeff[0]*fEdge[22]; 

  boundSurf_incr[1] = 0.9682458365518543*coeff[0]*fSkin[7]-1.25*coeff[0]*fSkin[1]+0.8660254037844386*coeff[0]*fSkin[0]; 
  boundSurf_incr[4] = 0.9682458365518543*coeff[0]*fSkin[11]-1.25*coeff[0]*fSkin[4]+0.8660254037844386*coeff[0]*fSkin[2]; 
  boundSurf_incr[5] = 0.9682458365518543*coeff[0]*fSkin[13]-1.25*coeff[0]*fSkin[5]+0.8660254037844386*coeff[0]*fSkin[3]; 
  boundSurf_incr[7] = -(3.75*coeff[0]*fSkin[7])+4.841229182759272*coeff[0]*fSkin[1]-3.3541019662496847*coeff[0]*fSkin[0]; 
  boundSurf_incr[10] = 0.9682458365518543*coeff[0]*fSkin[17]-1.25*coeff[0]*fSkin[10]+0.8660254037844386*coeff[0]*fSkin[6]; 
  boundSurf_incr[11] = -(3.75*coeff[0]*fSkin[11])+4.841229182759271*coeff[0]*fSkin[4]-3.3541019662496843*coeff[0]*fSkin[2]; 
  boundSurf_incr[12] = 0.9682458365518543*coeff[0]*fSkin[20]-1.25*coeff[0]*fSkin[12]+0.8660254037844387*coeff[0]*fSkin[8]; 
  boundSurf_incr[13] = -(3.75*coeff[0]*fSkin[13])+4.841229182759271*coeff[0]*fSkin[5]-3.3541019662496843*coeff[0]*fSkin[3]; 
  boundSurf_incr[15] = 0.9682458365518543*coeff[0]*fSkin[21]-1.25*coeff[0]*fSkin[15]+0.8660254037844387*coeff[0]*fSkin[9]; 
  boundSurf_incr[17] = -(3.75*coeff[0]*fSkin[17])+4.841229182759272*coeff[0]*fSkin[10]-3.3541019662496847*coeff[0]*fSkin[6]; 
  boundSurf_incr[18] = 0.9682458365518543*coeff[0]*fSkin[23]-1.25*coeff[0]*fSkin[18]+0.8660254037844387*coeff[0]*fSkin[14]; 
  boundSurf_incr[19] = 0.9682458365518543*coeff[0]*fSkin[24]-1.25*coeff[0]*fSkin[19]+0.8660254037844387*coeff[0]*fSkin[16]; 
  boundSurf_incr[20] = -(3.75*coeff[0]*fSkin[20])+4.841229182759271*coeff[0]*fSkin[12]-3.3541019662496847*coeff[0]*fSkin[8]; 
  boundSurf_incr[21] = -(3.75*coeff[0]*fSkin[21])+4.841229182759271*coeff[0]*fSkin[15]-3.3541019662496847*coeff[0]*fSkin[9]; 
  boundSurf_incr[23] = -(3.75*coeff[0]*fSkin[23])+4.841229182759272*coeff[0]*fSkin[18]-3.3541019662496843*coeff[0]*fSkin[14]; 
  boundSurf_incr[24] = -(3.75*coeff[0]*fSkin[24])+4.841229182759272*coeff[0]*fSkin[19]-3.3541019662496843*coeff[0]*fSkin[16]; 
  boundSurf_incr[25] = 0.9682458365518543*coeff[0]*fSkin[26]-1.25*coeff[0]*fSkin[25]+0.8660254037844386*coeff[0]*fSkin[22]; 
  boundSurf_incr[26] = -(3.75*coeff[0]*fSkin[26])+4.841229182759272*coeff[0]*fSkin[25]-3.3541019662496847*coeff[0]*fSkin[22]; 

  } else { 

  edgeSurf_incr[0] = -(0.6708203932499369*coeff[0]*fSkin[7])+0.6708203932499369*coeff[0]*fEdge[7]+1.190784930203603*coeff[0]*fSkin[1]+1.190784930203603*coeff[0]*fEdge[1]-0.9375*coeff[0]*fSkin[0]+0.9375*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 1.5855025573536612*coeff[0]*fSkin[7]-0.7382874503707886*coeff[0]*fEdge[7]-2.671875*coeff[0]*fSkin[1]-1.453125*coeff[0]*fEdge[1]+2.0568103339880417*coeff[0]*fSkin[0]-1.1907849302036029*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = -(0.6708203932499369*coeff[0]*fSkin[11])+0.6708203932499369*coeff[0]*fEdge[11]+1.190784930203603*coeff[0]*fSkin[4]+1.190784930203603*coeff[0]*fEdge[4]-0.9375*coeff[0]*fSkin[2]+0.9375*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = -(0.6708203932499369*coeff[0]*fSkin[13])+0.6708203932499369*coeff[0]*fEdge[13]+1.190784930203603*coeff[0]*fSkin[5]+1.190784930203603*coeff[0]*fEdge[5]-0.9375*coeff[0]*fSkin[3]+0.9375*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = 1.5855025573536612*coeff[0]*fSkin[11]-0.7382874503707888*coeff[0]*fEdge[11]-2.671875*coeff[0]*fSkin[4]-1.453125*coeff[0]*fEdge[4]+2.0568103339880417*coeff[0]*fSkin[2]-1.1907849302036029*coeff[0]*fEdge[2]; 
  edgeSurf_incr[5] = 1.5855025573536612*coeff[0]*fSkin[13]-0.7382874503707888*coeff[0]*fEdge[13]-2.671875*coeff[0]*fSkin[5]-1.453125*coeff[0]*fEdge[5]+2.0568103339880417*coeff[0]*fSkin[3]-1.1907849302036029*coeff[0]*fEdge[3]; 
  edgeSurf_incr[6] = -(0.6708203932499369*coeff[0]*fSkin[17])+0.6708203932499369*coeff[0]*fEdge[17]+1.190784930203603*coeff[0]*fSkin[10]+1.190784930203603*coeff[0]*fEdge[10]-0.9375*coeff[0]*fSkin[6]+0.9375*coeff[0]*fEdge[6]; 
  edgeSurf_incr[7] = -(3.140625*coeff[0]*fSkin[7])-0.140625*coeff[0]*fEdge[7]+5.022775277112744*coeff[0]*fSkin[1]+0.3025768239224549*coeff[0]*fEdge[1]-3.7733647120308955*coeff[0]*fSkin[0]+0.4192627457812108*coeff[0]*fEdge[0]; 
  edgeSurf_incr[8] = -(0.6708203932499369*coeff[0]*fSkin[20])+0.6708203932499369*coeff[0]*fEdge[20]+1.190784930203603*coeff[0]*fSkin[12]+1.190784930203603*coeff[0]*fEdge[12]-0.9375*coeff[0]*fSkin[8]+0.9375*coeff[0]*fEdge[8]; 
  edgeSurf_incr[9] = -(0.6708203932499369*coeff[0]*fSkin[21])+0.6708203932499369*coeff[0]*fEdge[21]+1.190784930203603*coeff[0]*fSkin[15]+1.190784930203603*coeff[0]*fEdge[15]-0.9375*coeff[0]*fSkin[9]+0.9375*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = 1.5855025573536612*coeff[0]*fSkin[17]-0.7382874503707886*coeff[0]*fEdge[17]-2.671875*coeff[0]*fSkin[10]-1.453125*coeff[0]*fEdge[10]+2.0568103339880417*coeff[0]*fSkin[6]-1.1907849302036029*coeff[0]*fEdge[6]; 
  edgeSurf_incr[11] = -(3.140625*coeff[0]*fSkin[11])-0.140625*coeff[0]*fEdge[11]+5.022775277112744*coeff[0]*fSkin[4]+0.30257682392245444*coeff[0]*fEdge[4]-3.773364712030894*coeff[0]*fSkin[2]+0.41926274578121053*coeff[0]*fEdge[2]; 
  edgeSurf_incr[12] = 1.5855025573536612*coeff[0]*fSkin[20]-0.7382874503707888*coeff[0]*fEdge[20]-2.671875*coeff[0]*fSkin[12]-1.453125*coeff[0]*fEdge[12]+2.0568103339880417*coeff[0]*fSkin[8]-1.190784930203603*coeff[0]*fEdge[8]; 
  edgeSurf_incr[13] = -(3.140625*coeff[0]*fSkin[13])-0.140625*coeff[0]*fEdge[13]+5.022775277112744*coeff[0]*fSkin[5]+0.30257682392245444*coeff[0]*fEdge[5]-3.773364712030894*coeff[0]*fSkin[3]+0.41926274578121053*coeff[0]*fEdge[3]; 
  edgeSurf_incr[14] = -(0.6708203932499369*coeff[0]*fSkin[23])+0.6708203932499369*coeff[0]*fEdge[23]+1.190784930203603*coeff[0]*fSkin[18]+1.190784930203603*coeff[0]*fEdge[18]-0.9375*coeff[0]*fSkin[14]+0.9375*coeff[0]*fEdge[14]; 
  edgeSurf_incr[15] = 1.5855025573536612*coeff[0]*fSkin[21]-0.7382874503707888*coeff[0]*fEdge[21]-2.671875*coeff[0]*fSkin[15]-1.453125*coeff[0]*fEdge[15]+2.0568103339880417*coeff[0]*fSkin[9]-1.190784930203603*coeff[0]*fEdge[9]; 
  edgeSurf_incr[16] = -(0.6708203932499369*coeff[0]*fSkin[24])+0.6708203932499369*coeff[0]*fEdge[24]+1.190784930203603*coeff[0]*fSkin[19]+1.190784930203603*coeff[0]*fEdge[19]-0.9375*coeff[0]*fSkin[16]+0.9375*coeff[0]*fEdge[16]; 
  edgeSurf_incr[17] = -(3.140625*coeff[0]*fSkin[17])-0.140625*coeff[0]*fEdge[17]+5.022775277112744*coeff[0]*fSkin[10]+0.3025768239224549*coeff[0]*fEdge[10]-3.7733647120308955*coeff[0]*fSkin[6]+0.4192627457812108*coeff[0]*fEdge[6]; 
  edgeSurf_incr[18] = 1.5855025573536612*coeff[0]*fSkin[23]-0.7382874503707886*coeff[0]*fEdge[23]-2.671875*coeff[0]*fSkin[18]-1.453125*coeff[0]*fEdge[18]+2.0568103339880417*coeff[0]*fSkin[14]-1.190784930203603*coeff[0]*fEdge[14]; 
  edgeSurf_incr[19] = 1.5855025573536612*coeff[0]*fSkin[24]-0.7382874503707886*coeff[0]*fEdge[24]-2.671875*coeff[0]*fSkin[19]-1.453125*coeff[0]*fEdge[19]+2.0568103339880417*coeff[0]*fSkin[16]-1.190784930203603*coeff[0]*fEdge[16]; 
  edgeSurf_incr[20] = -(3.140625*coeff[0]*fSkin[20])-0.140625*coeff[0]*fEdge[20]+5.022775277112744*coeff[0]*fSkin[12]+0.30257682392245444*coeff[0]*fEdge[12]-3.7733647120308955*coeff[0]*fSkin[8]+0.4192627457812108*coeff[0]*fEdge[8]; 
  edgeSurf_incr[21] = -(3.140625*coeff[0]*fSkin[21])-0.140625*coeff[0]*fEdge[21]+5.022775277112744*coeff[0]*fSkin[15]+0.30257682392245444*coeff[0]*fEdge[15]-3.7733647120308955*coeff[0]*fSkin[9]+0.4192627457812108*coeff[0]*fEdge[9]; 
  edgeSurf_incr[22] = -(0.6708203932499369*coeff[0]*fSkin[26])+0.6708203932499369*coeff[0]*fEdge[26]+1.190784930203603*coeff[0]*fSkin[25]+1.190784930203603*coeff[0]*fEdge[25]-0.9375*coeff[0]*fSkin[22]+0.9375*coeff[0]*fEdge[22]; 
  edgeSurf_incr[23] = -(3.140625*coeff[0]*fSkin[23])-0.140625*coeff[0]*fEdge[23]+5.022775277112744*coeff[0]*fSkin[18]+0.3025768239224549*coeff[0]*fEdge[18]-3.773364712030894*coeff[0]*fSkin[14]+0.41926274578121053*coeff[0]*fEdge[14]; 
  edgeSurf_incr[24] = -(3.140625*coeff[0]*fSkin[24])-0.140625*coeff[0]*fEdge[24]+5.022775277112744*coeff[0]*fSkin[19]+0.3025768239224549*coeff[0]*fEdge[19]-3.773364712030894*coeff[0]*fSkin[16]+0.41926274578121053*coeff[0]*fEdge[16]; 
  edgeSurf_incr[25] = 1.5855025573536612*coeff[0]*fSkin[26]-0.7382874503707886*coeff[0]*fEdge[26]-2.671875*coeff[0]*fSkin[25]-1.453125*coeff[0]*fEdge[25]+2.0568103339880417*coeff[0]*fSkin[22]-1.1907849302036029*coeff[0]*fEdge[22]; 
  edgeSurf_incr[26] = -(3.140625*coeff[0]*fSkin[26])-0.140625*coeff[0]*fEdge[26]+5.022775277112744*coeff[0]*fSkin[25]+0.3025768239224549*coeff[0]*fEdge[25]-3.7733647120308955*coeff[0]*fSkin[22]+0.4192627457812108*coeff[0]*fEdge[22]; 

  boundSurf_incr[1] = -(0.9682458365518543*coeff[0]*fSkin[7])-1.25*coeff[0]*fSkin[1]-0.8660254037844386*coeff[0]*fSkin[0]; 
  boundSurf_incr[4] = -(0.9682458365518543*coeff[0]*fSkin[11])-1.25*coeff[0]*fSkin[4]-0.8660254037844386*coeff[0]*fSkin[2]; 
  boundSurf_incr[5] = -(0.9682458365518543*coeff[0]*fSkin[13])-1.25*coeff[0]*fSkin[5]-0.8660254037844386*coeff[0]*fSkin[3]; 
  boundSurf_incr[7] = -(3.75*coeff[0]*fSkin[7])-4.841229182759272*coeff[0]*fSkin[1]-3.3541019662496847*coeff[0]*fSkin[0]; 
  boundSurf_incr[10] = -(0.9682458365518543*coeff[0]*fSkin[17])-1.25*coeff[0]*fSkin[10]-0.8660254037844386*coeff[0]*fSkin[6]; 
  boundSurf_incr[11] = -(3.75*coeff[0]*fSkin[11])-4.841229182759271*coeff[0]*fSkin[4]-3.3541019662496843*coeff[0]*fSkin[2]; 
  boundSurf_incr[12] = -(0.9682458365518543*coeff[0]*fSkin[20])-1.25*coeff[0]*fSkin[12]-0.8660254037844387*coeff[0]*fSkin[8]; 
  boundSurf_incr[13] = -(3.75*coeff[0]*fSkin[13])-4.841229182759271*coeff[0]*fSkin[5]-3.3541019662496843*coeff[0]*fSkin[3]; 
  boundSurf_incr[15] = -(0.9682458365518543*coeff[0]*fSkin[21])-1.25*coeff[0]*fSkin[15]-0.8660254037844387*coeff[0]*fSkin[9]; 
  boundSurf_incr[17] = -(3.75*coeff[0]*fSkin[17])-4.841229182759272*coeff[0]*fSkin[10]-3.3541019662496847*coeff[0]*fSkin[6]; 
  boundSurf_incr[18] = -(0.9682458365518543*coeff[0]*fSkin[23])-1.25*coeff[0]*fSkin[18]-0.8660254037844387*coeff[0]*fSkin[14]; 
  boundSurf_incr[19] = -(0.9682458365518543*coeff[0]*fSkin[24])-1.25*coeff[0]*fSkin[19]-0.8660254037844387*coeff[0]*fSkin[16]; 
  boundSurf_incr[20] = -(3.75*coeff[0]*fSkin[20])-4.841229182759271*coeff[0]*fSkin[12]-3.3541019662496847*coeff[0]*fSkin[8]; 
  boundSurf_incr[21] = -(3.75*coeff[0]*fSkin[21])-4.841229182759271*coeff[0]*fSkin[15]-3.3541019662496847*coeff[0]*fSkin[9]; 
  boundSurf_incr[23] = -(3.75*coeff[0]*fSkin[23])-4.841229182759272*coeff[0]*fSkin[18]-3.3541019662496843*coeff[0]*fSkin[14]; 
  boundSurf_incr[24] = -(3.75*coeff[0]*fSkin[24])-4.841229182759272*coeff[0]*fSkin[19]-3.3541019662496843*coeff[0]*fSkin[16]; 
  boundSurf_incr[25] = -(0.9682458365518543*coeff[0]*fSkin[26])-1.25*coeff[0]*fSkin[25]-0.8660254037844386*coeff[0]*fSkin[22]; 
  boundSurf_incr[26] = -(3.75*coeff[0]*fSkin[26])-4.841229182759272*coeff[0]*fSkin[25]-3.3541019662496847*coeff[0]*fSkin[22]; 

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
  out[20] += (vol_incr[20]+edgeSurf_incr[20]+boundSurf_incr[20])*Jfac; 
  out[21] += (vol_incr[21]+edgeSurf_incr[21]+boundSurf_incr[21])*Jfac; 
  out[22] += (vol_incr[22]+edgeSurf_incr[22]+boundSurf_incr[22])*Jfac; 
  out[23] += (vol_incr[23]+edgeSurf_incr[23]+boundSurf_incr[23])*Jfac; 
  out[24] += (vol_incr[24]+edgeSurf_incr[24]+boundSurf_incr[24])*Jfac; 
  out[25] += (vol_incr[25]+edgeSurf_incr[25]+boundSurf_incr[25])*Jfac; 
  out[26] += (vol_incr[26]+edgeSurf_incr[26]+boundSurf_incr[26])*Jfac; 

  return 0.;
}

