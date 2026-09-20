#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvx_1x1v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[2]: Cell-center coordinates. 
  // dxv[2]: Cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[4]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fSkin/Edge: Distribution function in cells 
  // out: Incremented distribution function in cell 
  const double *nuVtSqSum = &nuPrimMomsSum[2];

  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 

  double vol_incr[6] = {0.0}; 
  vol_incr[4] = 4.743416490252569*fSkin[1]*nuVtSqSum[1]+4.743416490252569*fSkin[0]*nuVtSqSum[0]; 
  vol_incr[5] = 4.743416490252569*fSkin[0]*nuVtSqSum[1]+4.743416490252569*nuVtSqSum[0]*fSkin[1]; 

  double edgeSurf_incr[6] = {0.0}; 
  double boundSurf_incr[6] = {0.0}; 

  if (edge == -1) { 

  double edgeSurf[6] = {0.0}; 
  edgeSurf[0] = -(0.6708203932499369*fSkin[4])+0.6708203932499369*fEdge[4]-1.190784930203603*fSkin[2]-1.190784930203603*fEdge[2]-0.9375*fSkin[0]+0.9375*fEdge[0]; 
  edgeSurf[1] = -(0.6708203932499369*fSkin[5])+0.6708203932499369*fEdge[5]-1.190784930203603*fSkin[3]-1.190784930203603*fEdge[3]-0.9375*fSkin[1]+0.9375*fEdge[1]; 
  edgeSurf[2] = -(1.5855025573536612*fSkin[4])+0.7382874503707888*fEdge[4]-2.671875*fSkin[2]-1.453125*fEdge[2]-2.0568103339880417*fSkin[0]+1.190784930203603*fEdge[0]; 
  edgeSurf[3] = -(1.5855025573536612*fSkin[5])+0.7382874503707888*fEdge[5]-2.671875*fSkin[3]-1.453125*fEdge[3]-2.0568103339880417*fSkin[1]+1.190784930203603*fEdge[1]; 
  edgeSurf[4] = -(3.140625*fSkin[4])-0.140625*fEdge[4]-5.022775277112744*fSkin[2]-0.3025768239224545*fEdge[2]-3.7733647120308955*fSkin[0]+0.4192627457812106*fEdge[0]; 
  edgeSurf[5] = -(3.140625*fSkin[5])-0.140625*fEdge[5]-5.022775277112744*fSkin[3]-0.30257682392245444*fEdge[3]-3.773364712030894*fSkin[1]+0.41926274578121053*fEdge[1]; 

  double boundSurf[6] = {0.0}; 
  boundSurf[2] = 1.9364916731037085*fSkin[4]-1.5*fSkin[2]+0.8660254037844386*fSkin[0]; 
  boundSurf[3] = 1.9364916731037085*fSkin[5]-1.5*fSkin[3]+0.8660254037844386*fSkin[1]; 
  boundSurf[4] = -(7.5*fSkin[4])+5.809475019311125*fSkin[2]-3.3541019662496847*fSkin[0]; 
  boundSurf[5] = -(7.5*fSkin[5])+5.809475019311126*fSkin[3]-3.3541019662496843*fSkin[1]; 

  edgeSurf_incr[0] = 0.7071067811865475*edgeSurf[1]*nuVtSqSum[1]+0.7071067811865475*edgeSurf[0]*nuVtSqSum[0]; 
  edgeSurf_incr[1] = 0.7071067811865475*edgeSurf[0]*nuVtSqSum[1]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[1]; 
  edgeSurf_incr[2] = 0.7071067811865475*nuVtSqSum[1]*edgeSurf[3]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[2]; 
  edgeSurf_incr[3] = 0.7071067811865475*nuVtSqSum[0]*edgeSurf[3]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[2]; 
  edgeSurf_incr[4] = 0.7071067811865475*nuVtSqSum[1]*edgeSurf[5]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[4]; 
  edgeSurf_incr[5] = 0.7071067811865475*nuVtSqSum[0]*edgeSurf[5]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[4]; 

  boundSurf_incr[0] = 0.7071067811865475*boundSurf[1]*nuVtSqSum[1]+0.7071067811865475*boundSurf[0]*nuVtSqSum[0]; 
  boundSurf_incr[1] = 0.7071067811865475*boundSurf[0]*nuVtSqSum[1]+0.7071067811865475*nuVtSqSum[0]*boundSurf[1]; 
  boundSurf_incr[2] = 0.7071067811865475*nuVtSqSum[1]*boundSurf[3]+0.7071067811865475*nuVtSqSum[0]*boundSurf[2]; 
  boundSurf_incr[3] = 0.7071067811865475*nuVtSqSum[0]*boundSurf[3]+0.7071067811865475*nuVtSqSum[1]*boundSurf[2]; 
  boundSurf_incr[4] = 0.7071067811865475*nuVtSqSum[1]*boundSurf[5]+0.7071067811865475*nuVtSqSum[0]*boundSurf[4]; 
  boundSurf_incr[5] = 0.7071067811865475*nuVtSqSum[0]*boundSurf[5]+0.7071067811865475*nuVtSqSum[1]*boundSurf[4]; 


  } else { 

  double edgeSurf[6] = {0.0}; 
  edgeSurf[0] = -(0.6708203932499369*fSkin[4])+0.6708203932499369*fEdge[4]+1.190784930203603*fSkin[2]+1.190784930203603*fEdge[2]-0.9375*fSkin[0]+0.9375*fEdge[0]; 
  edgeSurf[1] = -(0.6708203932499369*fSkin[5])+0.6708203932499369*fEdge[5]+1.190784930203603*fSkin[3]+1.190784930203603*fEdge[3]-0.9375*fSkin[1]+0.9375*fEdge[1]; 
  edgeSurf[2] = 1.5855025573536612*fSkin[4]-0.7382874503707888*fEdge[4]-2.671875*fSkin[2]-1.453125*fEdge[2]+2.0568103339880417*fSkin[0]-1.190784930203603*fEdge[0]; 
  edgeSurf[3] = 1.5855025573536612*fSkin[5]-0.7382874503707888*fEdge[5]-2.671875*fSkin[3]-1.453125*fEdge[3]+2.0568103339880417*fSkin[1]-1.190784930203603*fEdge[1]; 
  edgeSurf[4] = -(3.140625*fSkin[4])-0.140625*fEdge[4]+5.022775277112744*fSkin[2]+0.3025768239224545*fEdge[2]-3.7733647120308955*fSkin[0]+0.4192627457812106*fEdge[0]; 
  edgeSurf[5] = -(3.140625*fSkin[5])-0.140625*fEdge[5]+5.022775277112744*fSkin[3]+0.30257682392245444*fEdge[3]-3.773364712030894*fSkin[1]+0.41926274578121053*fEdge[1]; 

  double boundSurf[6] = {0.0}; 
  boundSurf[2] = -(1.9364916731037085*fSkin[4])-1.5*fSkin[2]-0.8660254037844386*fSkin[0]; 
  boundSurf[3] = -(1.9364916731037085*fSkin[5])-1.5*fSkin[3]-0.8660254037844386*fSkin[1]; 
  boundSurf[4] = -(7.5*fSkin[4])-5.809475019311125*fSkin[2]-3.3541019662496847*fSkin[0]; 
  boundSurf[5] = -(7.5*fSkin[5])-5.809475019311126*fSkin[3]-3.3541019662496843*fSkin[1]; 

  edgeSurf_incr[0] = 0.7071067811865475*edgeSurf[1]*nuVtSqSum[1]+0.7071067811865475*edgeSurf[0]*nuVtSqSum[0]; 
  edgeSurf_incr[1] = 0.7071067811865475*edgeSurf[0]*nuVtSqSum[1]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[1]; 
  edgeSurf_incr[2] = 0.7071067811865475*nuVtSqSum[1]*edgeSurf[3]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[2]; 
  edgeSurf_incr[3] = 0.7071067811865475*nuVtSqSum[0]*edgeSurf[3]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[2]; 
  edgeSurf_incr[4] = 0.7071067811865475*nuVtSqSum[1]*edgeSurf[5]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[4]; 
  edgeSurf_incr[5] = 0.7071067811865475*nuVtSqSum[0]*edgeSurf[5]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[4]; 

  boundSurf_incr[0] = 0.7071067811865475*boundSurf[1]*nuVtSqSum[1]+0.7071067811865475*boundSurf[0]*nuVtSqSum[0]; 
  boundSurf_incr[1] = 0.7071067811865475*boundSurf[0]*nuVtSqSum[1]+0.7071067811865475*nuVtSqSum[0]*boundSurf[1]; 
  boundSurf_incr[2] = 0.7071067811865475*nuVtSqSum[1]*boundSurf[3]+0.7071067811865475*nuVtSqSum[0]*boundSurf[2]; 
  boundSurf_incr[3] = 0.7071067811865475*nuVtSqSum[0]*boundSurf[3]+0.7071067811865475*nuVtSqSum[1]*boundSurf[2]; 
  boundSurf_incr[4] = 0.7071067811865475*nuVtSqSum[1]*boundSurf[5]+0.7071067811865475*nuVtSqSum[0]*boundSurf[4]; 
  boundSurf_incr[5] = 0.7071067811865475*nuVtSqSum[0]*boundSurf[5]+0.7071067811865475*nuVtSqSum[1]*boundSurf[4]; 

  } 

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*rdvSq4; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*rdvSq4; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*rdvSq4; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*rdvSq4; 
  out[4] += (vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*rdvSq4; 
  out[5] += (vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*rdvSq4; 

  return 0.;

} 
