#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvx_2x1v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[3]: Cell-center coordinates. 
  // dxv[3]: Cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[8]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fSkin/Edge: Distribution function in cells 
  // out: Incremented distribution function in cell 
  const double *nuVtSqSum = &nuPrimMomsSum[4];

  double rdvSq4 = 4.0/(dxv[2]*dxv[2]); 

  double g[12] = {0.0}; 
  g[7] = 6.708203932499369*fSkin[0]; 
  g[9] = 6.7082039324993685*fSkin[1]; 
  g[10] = 6.7082039324993685*fSkin[2]; 
  g[11] = 6.708203932499369*fSkin[4]; 

  if (edge == -1) { 

  g[0] += -(0.6708203932499369*fSkin[7])+0.6708203932499369*fEdge[7]-1.190784930203603*fSkin[3]-1.190784930203603*fEdge[3]-0.9375*fSkin[0]+0.9375*fEdge[0]; 
  g[1] += -(0.6708203932499369*fSkin[9])+0.6708203932499369*fEdge[9]-1.190784930203603*fSkin[5]-1.190784930203603*fEdge[5]-0.9375*fSkin[1]+0.9375*fEdge[1]; 
  g[2] += -(0.6708203932499369*fSkin[10])+0.6708203932499369*fEdge[10]-1.190784930203603*fSkin[6]-1.190784930203603*fEdge[6]-0.9375*fSkin[2]+0.9375*fEdge[2]; 
  g[3] += 0.3509891157500473*fSkin[7]+0.7382874503707888*fEdge[7]-4.171875*fSkin[3]-1.453125*fEdge[3]-1.190784930203603*fSkin[0]+1.190784930203603*fEdge[0]; 
  g[4] += -(0.6708203932499369*fSkin[11])+0.6708203932499369*fEdge[11]-1.190784930203603*fSkin[8]-1.190784930203603*fEdge[8]-0.9375*fSkin[4]+0.9375*fEdge[4]; 
  g[5] += 0.35098911575004715*fSkin[9]+0.7382874503707888*fEdge[9]-4.171875*fSkin[5]-1.453125*fEdge[5]-1.190784930203603*fSkin[1]+1.190784930203603*fEdge[1]; 
  g[6] += 0.35098911575004715*fSkin[10]+0.7382874503707888*fEdge[10]-4.171875*fSkin[6]-1.453125*fEdge[6]-1.190784930203603*fSkin[2]+1.190784930203603*fEdge[2]; 
  g[7] += -(10.640625*fSkin[7])-0.140625*fEdge[7]+0.7866997421983815*fSkin[3]-0.3025768239224545*fEdge[3]-7.12746667828058*fSkin[0]+0.4192627457812106*fEdge[0]; 
  g[8] += 0.3509891157500473*fSkin[11]+0.7382874503707888*fEdge[11]-4.171875*fSkin[8]-1.453125*fEdge[8]-1.190784930203603*fSkin[4]+1.190784930203603*fEdge[4]; 
  g[9] += -(10.640625*fSkin[9])-0.140625*fEdge[9]+0.7866997421983816*fSkin[5]-0.30257682392245444*fEdge[5]-7.1274666782805784*fSkin[1]+0.41926274578121053*fEdge[1]; 
  g[10] += -(10.640625*fSkin[10])-0.140625*fEdge[10]+0.7866997421983816*fSkin[6]-0.30257682392245444*fEdge[6]-7.1274666782805784*fSkin[2]+0.41926274578121053*fEdge[2]; 
  g[11] += -(10.640625*fSkin[11])-0.140625*fEdge[11]+0.7866997421983815*fSkin[8]-0.3025768239224545*fEdge[8]-7.12746667828058*fSkin[4]+0.4192627457812106*fEdge[4]; 

  } else { 

  g[0] += -(0.6708203932499369*fSkin[7])+0.6708203932499369*fEdge[7]+1.190784930203603*fSkin[3]+1.190784930203603*fEdge[3]-0.9375*fSkin[0]+0.9375*fEdge[0]; 
  g[1] += -(0.6708203932499369*fSkin[9])+0.6708203932499369*fEdge[9]+1.190784930203603*fSkin[5]+1.190784930203603*fEdge[5]-0.9375*fSkin[1]+0.9375*fEdge[1]; 
  g[2] += -(0.6708203932499369*fSkin[10])+0.6708203932499369*fEdge[10]+1.190784930203603*fSkin[6]+1.190784930203603*fEdge[6]-0.9375*fSkin[2]+0.9375*fEdge[2]; 
  g[3] += -(0.3509891157500473*fSkin[7])-0.7382874503707888*fEdge[7]-4.171875*fSkin[3]-1.453125*fEdge[3]+1.190784930203603*fSkin[0]-1.190784930203603*fEdge[0]; 
  g[4] += -(0.6708203932499369*fSkin[11])+0.6708203932499369*fEdge[11]+1.190784930203603*fSkin[8]+1.190784930203603*fEdge[8]-0.9375*fSkin[4]+0.9375*fEdge[4]; 
  g[5] += -(0.35098911575004715*fSkin[9])-0.7382874503707888*fEdge[9]-4.171875*fSkin[5]-1.453125*fEdge[5]+1.190784930203603*fSkin[1]-1.190784930203603*fEdge[1]; 
  g[6] += -(0.35098911575004715*fSkin[10])-0.7382874503707888*fEdge[10]-4.171875*fSkin[6]-1.453125*fEdge[6]+1.190784930203603*fSkin[2]-1.190784930203603*fEdge[2]; 
  g[7] += -(10.640625*fSkin[7])-0.140625*fEdge[7]-0.7866997421983815*fSkin[3]+0.3025768239224545*fEdge[3]-7.12746667828058*fSkin[0]+0.4192627457812106*fEdge[0]; 
  g[8] += -(0.3509891157500473*fSkin[11])-0.7382874503707888*fEdge[11]-4.171875*fSkin[8]-1.453125*fEdge[8]+1.190784930203603*fSkin[4]-1.190784930203603*fEdge[4]; 
  g[9] += -(10.640625*fSkin[9])-0.140625*fEdge[9]-0.7866997421983816*fSkin[5]+0.30257682392245444*fEdge[5]-7.1274666782805784*fSkin[1]+0.41926274578121053*fEdge[1]; 
  g[10] += -(10.640625*fSkin[10])-0.140625*fEdge[10]-0.7866997421983816*fSkin[6]+0.30257682392245444*fEdge[6]-7.1274666782805784*fSkin[2]+0.41926274578121053*fEdge[2]; 
  g[11] += -(10.640625*fSkin[11])-0.140625*fEdge[11]-0.7866997421983815*fSkin[8]+0.3025768239224545*fEdge[8]-7.12746667828058*fSkin[4]+0.4192627457812106*fEdge[4]; 

  } 

  out[0] += 0.5*nuVtSqSum[3]*g[4]*rdvSq4+0.5*g[2]*nuVtSqSum[2]*rdvSq4+0.5*g[1]*nuVtSqSum[1]*rdvSq4+0.5*g[0]*nuVtSqSum[0]*rdvSq4; 
  out[1] += 0.5*nuVtSqSum[2]*g[4]*rdvSq4+0.5*g[2]*nuVtSqSum[3]*rdvSq4+0.5*g[0]*nuVtSqSum[1]*rdvSq4+0.5*nuVtSqSum[0]*g[1]*rdvSq4; 
  out[2] += 0.5*nuVtSqSum[1]*g[4]*rdvSq4+0.5*g[1]*nuVtSqSum[3]*rdvSq4+0.5*g[0]*nuVtSqSum[2]*rdvSq4+0.5*nuVtSqSum[0]*g[2]*rdvSq4; 
  out[3] += 0.5*nuVtSqSum[3]*g[8]*rdvSq4+0.5*nuVtSqSum[2]*g[6]*rdvSq4+0.5*nuVtSqSum[1]*g[5]*rdvSq4+0.5*nuVtSqSum[0]*g[3]*rdvSq4; 
  out[4] += 0.5*nuVtSqSum[0]*g[4]*rdvSq4+0.5*g[0]*nuVtSqSum[3]*rdvSq4+0.5*g[1]*nuVtSqSum[2]*rdvSq4+0.5*nuVtSqSum[1]*g[2]*rdvSq4; 
  out[5] += 0.5*nuVtSqSum[2]*g[8]*rdvSq4+0.5*nuVtSqSum[3]*g[6]*rdvSq4+0.5*nuVtSqSum[0]*g[5]*rdvSq4+0.5*nuVtSqSum[1]*g[3]*rdvSq4; 
  out[6] += 0.5*nuVtSqSum[1]*g[8]*rdvSq4+0.5*nuVtSqSum[0]*g[6]*rdvSq4+0.5*nuVtSqSum[3]*g[5]*rdvSq4+0.5*nuVtSqSum[2]*g[3]*rdvSq4; 
  out[7] += 0.5*nuVtSqSum[3]*g[11]*rdvSq4+0.5000000000000001*nuVtSqSum[2]*g[10]*rdvSq4+0.5000000000000001*nuVtSqSum[1]*g[9]*rdvSq4+0.5*nuVtSqSum[0]*g[7]*rdvSq4; 
  out[8] += 0.5*nuVtSqSum[0]*g[8]*rdvSq4+0.5*nuVtSqSum[1]*g[6]*rdvSq4+0.5*nuVtSqSum[2]*g[5]*rdvSq4+0.5*g[3]*nuVtSqSum[3]*rdvSq4; 
  out[9] += 0.5000000000000001*nuVtSqSum[2]*g[11]*rdvSq4+0.5*nuVtSqSum[3]*g[10]*rdvSq4+0.5*nuVtSqSum[0]*g[9]*rdvSq4+0.5000000000000001*nuVtSqSum[1]*g[7]*rdvSq4; 
  out[10] += 0.5000000000000001*nuVtSqSum[1]*g[11]*rdvSq4+0.5*nuVtSqSum[0]*g[10]*rdvSq4+0.5*nuVtSqSum[3]*g[9]*rdvSq4+0.5000000000000001*nuVtSqSum[2]*g[7]*rdvSq4; 
  out[11] += 0.5*nuVtSqSum[0]*g[11]*rdvSq4+0.5000000000000001*nuVtSqSum[1]*g[10]*rdvSq4+0.5000000000000001*nuVtSqSum[2]*g[9]*rdvSq4+0.5*nuVtSqSum[3]*g[7]*rdvSq4; 

  return 0.;

} 
