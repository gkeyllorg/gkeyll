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

  double g[6] = {0.0}; 
  g[4] = 6.708203932499369*fSkin[0]; 
  g[5] = 6.7082039324993685*fSkin[1]; 

  if (edge == -1) { 

  g[0] += -(0.6708203932499369*fSkin[4])+0.6708203932499369*fEdge[4]-1.190784930203603*fSkin[2]-1.190784930203603*fEdge[2]-0.9375*fSkin[0]+0.9375*fEdge[0]; 
  g[1] += -(0.6708203932499369*fSkin[5])+0.6708203932499369*fEdge[5]-1.190784930203603*fSkin[3]-1.190784930203603*fEdge[3]-0.9375*fSkin[1]+0.9375*fEdge[1]; 
  g[2] += 0.3509891157500473*fSkin[4]+0.7382874503707888*fEdge[4]-4.171875*fSkin[2]-1.453125*fEdge[2]-1.190784930203603*fSkin[0]+1.190784930203603*fEdge[0]; 
  g[3] += 0.35098911575004715*fSkin[5]+0.7382874503707888*fEdge[5]-4.171875*fSkin[3]-1.453125*fEdge[3]-1.190784930203603*fSkin[1]+1.190784930203603*fEdge[1]; 
  g[4] += -(10.640625*fSkin[4])-0.140625*fEdge[4]+0.7866997421983815*fSkin[2]-0.3025768239224545*fEdge[2]-7.12746667828058*fSkin[0]+0.4192627457812106*fEdge[0]; 
  g[5] += -(10.640625*fSkin[5])-0.140625*fEdge[5]+0.7866997421983816*fSkin[3]-0.30257682392245444*fEdge[3]-7.1274666782805784*fSkin[1]+0.41926274578121053*fEdge[1]; 

  } else { 

  g[0] += -(0.6708203932499369*fSkin[4])+0.6708203932499369*fEdge[4]+1.190784930203603*fSkin[2]+1.190784930203603*fEdge[2]-0.9375*fSkin[0]+0.9375*fEdge[0]; 
  g[1] += -(0.6708203932499369*fSkin[5])+0.6708203932499369*fEdge[5]+1.190784930203603*fSkin[3]+1.190784930203603*fEdge[3]-0.9375*fSkin[1]+0.9375*fEdge[1]; 
  g[2] += -(0.3509891157500473*fSkin[4])-0.7382874503707888*fEdge[4]-4.171875*fSkin[2]-1.453125*fEdge[2]+1.190784930203603*fSkin[0]-1.190784930203603*fEdge[0]; 
  g[3] += -(0.35098911575004715*fSkin[5])-0.7382874503707888*fEdge[5]-4.171875*fSkin[3]-1.453125*fEdge[3]+1.190784930203603*fSkin[1]-1.190784930203603*fEdge[1]; 
  g[4] += -(10.640625*fSkin[4])-0.140625*fEdge[4]-0.7866997421983815*fSkin[2]+0.3025768239224545*fEdge[2]-7.12746667828058*fSkin[0]+0.4192627457812106*fEdge[0]; 
  g[5] += -(10.640625*fSkin[5])-0.140625*fEdge[5]-0.7866997421983816*fSkin[3]+0.30257682392245444*fEdge[3]-7.1274666782805784*fSkin[1]+0.41926274578121053*fEdge[1]; 

  } 

  out[0] += 0.7071067811865475*g[1]*nuVtSqSum[1]*rdvSq4+0.7071067811865475*g[0]*nuVtSqSum[0]*rdvSq4; 
  out[1] += 0.7071067811865475*g[0]*nuVtSqSum[1]*rdvSq4+0.7071067811865475*nuVtSqSum[0]*g[1]*rdvSq4; 
  out[2] += 0.7071067811865475*nuVtSqSum[1]*g[3]*rdvSq4+0.7071067811865475*nuVtSqSum[0]*g[2]*rdvSq4; 
  out[3] += 0.7071067811865475*nuVtSqSum[0]*g[3]*rdvSq4+0.7071067811865475*nuVtSqSum[1]*g[2]*rdvSq4; 
  out[4] += 0.7071067811865475*nuVtSqSum[1]*g[5]*rdvSq4+0.7071067811865475*nuVtSqSum[0]*g[4]*rdvSq4; 
  out[5] += 0.7071067811865475*nuVtSqSum[0]*g[5]*rdvSq4+0.7071067811865475*nuVtSqSum[1]*g[4]*rdvSq4; 

  return 0.;

} 
