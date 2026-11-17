#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_diff_surfvx_2x1v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[3]: cell-center coordinates. 
  // dxv[3]: cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[8]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fl/fc/fr: distribution function in cells 
  // out: incremented distribution function in cell 

  const double *nuVtSqSum = &nuPrimMomsSum[4];

  double rdvSq4 = 4.0/(dxv[2]*dxv[2]); 
  double incr[12]; 

  double f_xx[12] = {0.0}; 
  f_xx[0] = 0.6708203932499369*fr[7]+0.6708203932499369*fl[7]-1.3416407864998738*fc[7]-1.190784930203603*fr[3]+1.190784930203603*fl[3]+0.9375*fr[0]+0.9375*fl[0]-1.875*fc[0]; 
  f_xx[1] = 0.6708203932499369*fr[9]+0.6708203932499369*fl[9]-1.3416407864998738*fc[9]-1.190784930203603*fr[5]+1.190784930203603*fl[5]+0.9375*fr[1]+0.9375*fl[1]-1.875*fc[1]; 
  f_xx[2] = 0.6708203932499369*fr[10]+0.6708203932499369*fl[10]-1.3416407864998738*fc[10]-1.190784930203603*fr[6]+1.190784930203603*fl[6]+0.9375*fr[2]+0.9375*fl[2]-1.875*fc[2]; 
  f_xx[3] = 0.7382874503707888*fr[7]-0.7382874503707888*fl[7]-1.453125*fr[3]-1.453125*fl[3]-5.34375*fc[3]+1.190784930203603*fr[0]-1.190784930203603*fl[0]; 
  f_xx[4] = 0.6708203932499369*fr[11]+0.6708203932499369*fl[11]-1.3416407864998738*fc[11]-1.190784930203603*fr[8]+1.190784930203603*fl[8]+0.9375*fr[4]+0.9375*fl[4]-1.875*fc[4]; 
  f_xx[5] = 0.7382874503707888*fr[9]-0.7382874503707888*fl[9]-1.453125*fr[5]-1.453125*fl[5]-5.34375*fc[5]+1.190784930203603*fr[1]-1.190784930203603*fl[1]; 
  f_xx[6] = 0.7382874503707888*fr[10]-0.7382874503707888*fl[10]-1.453125*fr[6]-1.453125*fl[6]-5.34375*fc[6]+1.190784930203603*fr[2]-1.190784930203603*fl[2]; 
  f_xx[7] = -(0.140625*fr[7])-0.140625*fl[7]-6.28125*fc[7]-0.3025768239224545*fr[3]+0.3025768239224545*fl[3]+0.4192627457812106*fr[0]+0.4192627457812106*fl[0]-0.8385254915624212*fc[0]; 
  f_xx[8] = 0.7382874503707888*fr[11]-0.7382874503707888*fl[11]-1.453125*fr[8]-1.453125*fl[8]-5.34375*fc[8]+1.190784930203603*fr[4]-1.190784930203603*fl[4]; 
  f_xx[9] = -(0.140625*fr[9])-0.140625*fl[9]-6.28125*fc[9]-0.30257682392245444*fr[5]+0.30257682392245444*fl[5]+0.41926274578121053*fr[1]+0.41926274578121053*fl[1]-0.8385254915624211*fc[1]; 
  f_xx[10] = -(0.140625*fr[10])-0.140625*fl[10]-6.28125*fc[10]-0.30257682392245444*fr[6]+0.30257682392245444*fl[6]+0.41926274578121053*fr[2]+0.41926274578121053*fl[2]-0.8385254915624211*fc[2]; 
  f_xx[11] = -(0.140625*fr[11])-0.140625*fl[11]-6.28125*fc[11]-0.3025768239224545*fr[8]+0.3025768239224545*fl[8]+0.4192627457812106*fr[4]+0.4192627457812106*fl[4]-0.8385254915624212*fc[4]; 

  incr[0] = 0.5*nuVtSqSum[3]*f_xx[4]+0.5*f_xx[2]*nuVtSqSum[2]+0.5*f_xx[1]*nuVtSqSum[1]+0.5*f_xx[0]*nuVtSqSum[0]; 
  incr[1] = 0.5*nuVtSqSum[2]*f_xx[4]+0.5*f_xx[2]*nuVtSqSum[3]+0.5*f_xx[0]*nuVtSqSum[1]+0.5*nuVtSqSum[0]*f_xx[1]; 
  incr[2] = 0.5*nuVtSqSum[1]*f_xx[4]+0.5*f_xx[1]*nuVtSqSum[3]+0.5*f_xx[0]*nuVtSqSum[2]+0.5*nuVtSqSum[0]*f_xx[2]; 
  incr[3] = 0.5*nuVtSqSum[3]*f_xx[8]+0.5*nuVtSqSum[2]*f_xx[6]+0.5*nuVtSqSum[1]*f_xx[5]+0.5*nuVtSqSum[0]*f_xx[3]; 
  incr[4] = 0.5*nuVtSqSum[0]*f_xx[4]+0.5*f_xx[0]*nuVtSqSum[3]+0.5*f_xx[1]*nuVtSqSum[2]+0.5*nuVtSqSum[1]*f_xx[2]; 
  incr[5] = 0.5*nuVtSqSum[2]*f_xx[8]+0.5*nuVtSqSum[3]*f_xx[6]+0.5*nuVtSqSum[0]*f_xx[5]+0.5*nuVtSqSum[1]*f_xx[3]; 
  incr[6] = 0.5*nuVtSqSum[1]*f_xx[8]+0.5*nuVtSqSum[0]*f_xx[6]+0.5*nuVtSqSum[3]*f_xx[5]+0.5*nuVtSqSum[2]*f_xx[3]; 
  incr[7] = 0.5*nuVtSqSum[3]*f_xx[11]+0.5000000000000001*nuVtSqSum[2]*f_xx[10]+0.5000000000000001*nuVtSqSum[1]*f_xx[9]+0.5*nuVtSqSum[0]*f_xx[7]; 
  incr[8] = 0.5*nuVtSqSum[0]*f_xx[8]+0.5*nuVtSqSum[1]*f_xx[6]+0.5*nuVtSqSum[2]*f_xx[5]+0.5*f_xx[3]*nuVtSqSum[3]; 
  incr[9] = 0.5000000000000001*nuVtSqSum[2]*f_xx[11]+0.5*nuVtSqSum[3]*f_xx[10]+0.5*nuVtSqSum[0]*f_xx[9]+0.5000000000000001*nuVtSqSum[1]*f_xx[7]; 
  incr[10] = 0.5000000000000001*nuVtSqSum[1]*f_xx[11]+0.5*nuVtSqSum[0]*f_xx[10]+0.5*nuVtSqSum[3]*f_xx[9]+0.5000000000000001*nuVtSqSum[2]*f_xx[7]; 
  incr[11] = 0.5*nuVtSqSum[0]*f_xx[11]+0.5000000000000001*nuVtSqSum[1]*f_xx[10]+0.5000000000000001*nuVtSqSum[2]*f_xx[9]+0.5*nuVtSqSum[3]*f_xx[7]; 

  out[0] += incr[0]*rdvSq4; 
  out[1] += incr[1]*rdvSq4; 
  out[2] += incr[2]*rdvSq4; 
  out[3] += incr[3]*rdvSq4; 
  out[4] += incr[4]*rdvSq4; 
  out[5] += incr[5]*rdvSq4; 
  out[6] += incr[6]*rdvSq4; 
  out[7] += incr[7]*rdvSq4; 
  out[8] += incr[8]*rdvSq4; 
  out[9] += incr[9]*rdvSq4; 
  out[10] += incr[10]*rdvSq4; 
  out[11] += incr[11]*rdvSq4; 

  return 0.;

} 
