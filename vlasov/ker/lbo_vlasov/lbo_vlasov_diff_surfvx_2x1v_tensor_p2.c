#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_diff_surfvx_2x1v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[3]: cell-center coordinates. 
  // dxv[3]: cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[18]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fl/fc/fr: distribution function in cells 
  // out: incremented distribution function in cell 

  const double *nuVtSqSum = &nuPrimMomsSum[9];

  double rdvSq4 = 4.0/(dxv[2]*dxv[2]); 
  double incr[27]; 

  double f_xx[27] = {0.0}; 
  f_xx[0] = 0.6708203932499369*fr[9]+0.6708203932499369*fl[9]-1.3416407864998738*fc[9]-1.190784930203603*fr[3]+1.190784930203603*fl[3]+0.9375*fr[0]+0.9375*fl[0]-1.875*fc[0]; 
  f_xx[1] = 0.6708203932499369*fr[15]+0.6708203932499369*fl[15]-1.3416407864998738*fc[15]-1.190784930203603*fr[5]+1.190784930203603*fl[5]+0.9375*fr[1]+0.9375*fl[1]-1.875*fc[1]; 
  f_xx[2] = 0.6708203932499369*fr[16]+0.6708203932499369*fl[16]-1.3416407864998738*fc[16]-1.190784930203603*fr[6]+1.190784930203603*fl[6]+0.9375*fr[2]+0.9375*fl[2]-1.875*fc[2]; 
  f_xx[3] = 0.7382874503707888*fr[9]-0.7382874503707888*fl[9]-1.453125*fr[3]-1.453125*fl[3]-5.34375*fc[3]+1.190784930203603*fr[0]-1.190784930203603*fl[0]; 
  f_xx[4] = 0.6708203932499369*fr[19]+0.6708203932499369*fl[19]-1.3416407864998738*fc[19]-1.190784930203603*fr[10]+1.190784930203603*fl[10]+0.9375*fr[4]+0.9375*fl[4]-1.875*fc[4]; 
  f_xx[5] = 0.7382874503707888*fr[15]-0.7382874503707888*fl[15]-1.453125*fr[5]-1.453125*fl[5]-5.34375*fc[5]+1.190784930203603*fr[1]-1.190784930203603*fl[1]; 
  f_xx[6] = 0.7382874503707888*fr[16]-0.7382874503707888*fl[16]-1.453125*fr[6]-1.453125*fl[6]-5.34375*fc[6]+1.190784930203603*fr[2]-1.190784930203603*fl[2]; 
  f_xx[7] = 0.6708203932499369*fr[21]+0.6708203932499369*fl[21]-1.3416407864998738*fc[21]-1.190784930203603*fr[13]+1.190784930203603*fl[13]+0.9375*fr[7]+0.9375*fl[7]-1.875*fc[7]; 
  f_xx[8] = 0.6708203932499369*fr[22]+0.6708203932499369*fl[22]-1.3416407864998738*fc[22]-1.190784930203603*fr[14]+1.190784930203603*fl[14]+0.9375*fr[8]+0.9375*fl[8]-1.875*fc[8]; 
  f_xx[9] = -(0.140625*fr[9])-0.140625*fl[9]-6.28125*fc[9]-0.3025768239224545*fr[3]+0.3025768239224545*fl[3]+0.4192627457812106*fr[0]+0.4192627457812106*fl[0]-0.8385254915624212*fc[0]; 
  f_xx[10] = 0.7382874503707888*fr[19]-0.7382874503707888*fl[19]-1.453125*fr[10]-1.453125*fl[10]-5.34375*fc[10]+1.190784930203603*fr[4]-1.190784930203603*fl[4]; 
  f_xx[11] = 0.6708203932499369*fr[24]+0.6708203932499369*fl[24]-1.3416407864998738*fc[24]-1.190784930203603*fr[17]+1.190784930203603*fl[17]+0.9375*fr[11]+0.9375*fl[11]-1.875*fc[11]; 
  f_xx[12] = 0.6708203932499369*fr[25]+0.6708203932499369*fl[25]-1.3416407864998738*fc[25]-1.190784930203603*fr[18]+1.190784930203603*fl[18]+0.9375*fr[12]+0.9375*fl[12]-1.875*fc[12]; 
  f_xx[13] = 0.7382874503707888*fr[21]-0.7382874503707888*fl[21]-1.453125*fr[13]-1.453125*fl[13]-5.34375*fc[13]+1.190784930203603*fr[7]-1.190784930203603*fl[7]; 
  f_xx[14] = 0.7382874503707888*fr[22]-0.7382874503707888*fl[22]-1.453125*fr[14]-1.453125*fl[14]-5.34375*fc[14]+1.190784930203603*fr[8]-1.190784930203603*fl[8]; 
  f_xx[15] = -(0.140625*fr[15])-0.140625*fl[15]-6.28125*fc[15]-0.30257682392245444*fr[5]+0.30257682392245444*fl[5]+0.41926274578121053*fr[1]+0.41926274578121053*fl[1]-0.8385254915624211*fc[1]; 
  f_xx[16] = -(0.140625*fr[16])-0.140625*fl[16]-6.28125*fc[16]-0.30257682392245444*fr[6]+0.30257682392245444*fl[6]+0.41926274578121053*fr[2]+0.41926274578121053*fl[2]-0.8385254915624211*fc[2]; 
  f_xx[17] = 0.7382874503707888*fr[24]-0.7382874503707888*fl[24]-1.453125*fr[17]-1.453125*fl[17]-5.34375*fc[17]+1.190784930203603*fr[11]-1.190784930203603*fl[11]; 
  f_xx[18] = 0.7382874503707888*fr[25]-0.7382874503707888*fl[25]-1.453125*fr[18]-1.453125*fl[18]-5.34375*fc[18]+1.190784930203603*fr[12]-1.190784930203603*fl[12]; 
  f_xx[19] = -(0.140625*fr[19])-0.140625*fl[19]-6.28125*fc[19]-0.3025768239224545*fr[10]+0.3025768239224545*fl[10]+0.4192627457812106*fr[4]+0.4192627457812106*fl[4]-0.8385254915624212*fc[4]; 
  f_xx[20] = 0.6708203932499369*fr[26]+0.6708203932499369*fl[26]-1.3416407864998738*fc[26]-1.190784930203603*fr[23]+1.190784930203603*fl[23]+0.9375*fr[20]+0.9375*fl[20]-1.875*fc[20]; 
  f_xx[21] = -(0.140625*fr[21])-0.140625*fl[21]-6.28125*fc[21]-0.30257682392245444*fr[13]+0.30257682392245444*fl[13]+0.4192627457812106*fr[7]+0.4192627457812106*fl[7]-0.8385254915624212*fc[7]; 
  f_xx[22] = -(0.140625*fr[22])-0.140625*fl[22]-6.28125*fc[22]-0.30257682392245444*fr[14]+0.30257682392245444*fl[14]+0.4192627457812106*fr[8]+0.4192627457812106*fl[8]-0.8385254915624212*fc[8]; 
  f_xx[23] = 0.7382874503707888*fr[26]-0.7382874503707888*fl[26]-1.453125*fr[23]-1.453125*fl[23]-5.34375*fc[23]+1.190784930203603*fr[20]-1.190784930203603*fl[20]; 
  f_xx[24] = -(0.140625*fr[24])-0.140625*fl[24]-6.28125*fc[24]-0.3025768239224545*fr[17]+0.3025768239224545*fl[17]+0.41926274578121053*fr[11]+0.41926274578121053*fl[11]-0.8385254915624211*fc[11]; 
  f_xx[25] = -(0.140625*fr[25])-0.140625*fl[25]-6.28125*fc[25]-0.3025768239224545*fr[18]+0.3025768239224545*fl[18]+0.41926274578121053*fr[12]+0.41926274578121053*fl[12]-0.8385254915624211*fc[12]; 
  f_xx[26] = -(0.140625*fr[26])-0.140625*fl[26]-6.28125*fc[26]-0.3025768239224545*fr[23]+0.3025768239224545*fl[23]+0.4192627457812106*fr[20]+0.4192627457812106*fl[20]-0.8385254915624212*fc[20]; 

  incr[0] = 0.5*nuVtSqSum[8]*f_xx[20]+0.5*nuVtSqSum[7]*f_xx[12]+0.5*nuVtSqSum[6]*f_xx[11]+0.5*nuVtSqSum[5]*f_xx[8]+0.5*nuVtSqSum[4]*f_xx[7]+0.5*nuVtSqSum[3]*f_xx[4]+0.5*f_xx[2]*nuVtSqSum[2]+0.5*f_xx[1]*nuVtSqSum[1]+0.5*f_xx[0]*nuVtSqSum[0]; 
  incr[1] = 0.44721359549995804*nuVtSqSum[7]*f_xx[20]+0.44721359549995804*nuVtSqSum[8]*f_xx[12]+0.5000000000000001*nuVtSqSum[5]*f_xx[12]+0.44721359549995804*nuVtSqSum[3]*f_xx[11]+0.5000000000000001*nuVtSqSum[7]*f_xx[8]+0.4472135954999579*nuVtSqSum[1]*f_xx[7]+0.44721359549995804*f_xx[4]*nuVtSqSum[6]+0.4472135954999579*f_xx[1]*nuVtSqSum[4]+0.5*nuVtSqSum[2]*f_xx[4]+0.5*f_xx[2]*nuVtSqSum[3]+0.5*f_xx[0]*nuVtSqSum[1]+0.5*nuVtSqSum[0]*f_xx[1]; 
  incr[2] = 0.44721359549995804*nuVtSqSum[6]*f_xx[20]+0.44721359549995804*nuVtSqSum[3]*f_xx[12]+0.44721359549995804*nuVtSqSum[8]*f_xx[11]+0.5000000000000001*nuVtSqSum[4]*f_xx[11]+0.4472135954999579*nuVtSqSum[2]*f_xx[8]+0.44721359549995804*f_xx[4]*nuVtSqSum[7]+0.5000000000000001*nuVtSqSum[6]*f_xx[7]+0.4472135954999579*f_xx[2]*nuVtSqSum[5]+0.5*nuVtSqSum[1]*f_xx[4]+0.5*f_xx[1]*nuVtSqSum[3]+0.5*f_xx[0]*nuVtSqSum[2]+0.5*nuVtSqSum[0]*f_xx[2]; 
  incr[3] = 0.5*nuVtSqSum[8]*f_xx[23]+0.5000000000000001*nuVtSqSum[7]*f_xx[18]+0.5000000000000001*nuVtSqSum[6]*f_xx[17]+0.5000000000000001*nuVtSqSum[5]*f_xx[14]+0.5000000000000001*nuVtSqSum[4]*f_xx[13]+0.5*nuVtSqSum[3]*f_xx[10]+0.5*nuVtSqSum[2]*f_xx[6]+0.5*nuVtSqSum[1]*f_xx[5]+0.5*nuVtSqSum[0]*f_xx[3]; 
  incr[4] = 0.4*nuVtSqSum[3]*f_xx[20]+0.4*nuVtSqSum[6]*f_xx[12]+0.44721359549995804*nuVtSqSum[2]*f_xx[12]+0.4*nuVtSqSum[7]*f_xx[11]+0.44721359549995804*nuVtSqSum[1]*f_xx[11]+0.4*f_xx[4]*nuVtSqSum[8]+0.4472135954999579*nuVtSqSum[3]*f_xx[8]+0.44721359549995804*f_xx[2]*nuVtSqSum[7]+0.4472135954999579*nuVtSqSum[3]*f_xx[7]+0.44721359549995804*f_xx[1]*nuVtSqSum[6]+0.4472135954999579*f_xx[4]*nuVtSqSum[5]+0.4472135954999579*f_xx[4]*nuVtSqSum[4]+0.5*nuVtSqSum[0]*f_xx[4]+0.5*f_xx[0]*nuVtSqSum[3]+0.5*f_xx[1]*nuVtSqSum[2]+0.5*nuVtSqSum[1]*f_xx[2]; 
  incr[5] = 0.44721359549995804*nuVtSqSum[7]*f_xx[23]+0.4472135954999579*nuVtSqSum[8]*f_xx[18]+0.5*nuVtSqSum[5]*f_xx[18]+0.4472135954999579*nuVtSqSum[3]*f_xx[17]+0.5*nuVtSqSum[7]*f_xx[14]+0.44721359549995804*nuVtSqSum[1]*f_xx[13]+0.44721359549995804*nuVtSqSum[6]*f_xx[10]+0.5*nuVtSqSum[2]*f_xx[10]+0.5*nuVtSqSum[3]*f_xx[6]+0.4472135954999579*nuVtSqSum[4]*f_xx[5]+0.5*nuVtSqSum[0]*f_xx[5]+0.5*nuVtSqSum[1]*f_xx[3]; 
  incr[6] = 0.44721359549995804*nuVtSqSum[6]*f_xx[23]+0.4472135954999579*nuVtSqSum[3]*f_xx[18]+0.4472135954999579*nuVtSqSum[8]*f_xx[17]+0.5*nuVtSqSum[4]*f_xx[17]+0.44721359549995804*nuVtSqSum[2]*f_xx[14]+0.5*nuVtSqSum[6]*f_xx[13]+0.44721359549995804*nuVtSqSum[7]*f_xx[10]+0.5*nuVtSqSum[1]*f_xx[10]+0.4472135954999579*nuVtSqSum[5]*f_xx[6]+0.5*nuVtSqSum[0]*f_xx[6]+0.5*nuVtSqSum[3]*f_xx[5]+0.5*nuVtSqSum[2]*f_xx[3]; 
  incr[7] = 0.31943828249996997*nuVtSqSum[8]*f_xx[20]+0.5*nuVtSqSum[5]*f_xx[20]+0.4472135954999579*nuVtSqSum[7]*f_xx[12]+0.31943828249996997*nuVtSqSum[6]*f_xx[11]+0.5000000000000001*nuVtSqSum[2]*f_xx[11]+0.5*f_xx[8]*nuVtSqSum[8]+0.31943828249996997*nuVtSqSum[4]*f_xx[7]+0.5*nuVtSqSum[0]*f_xx[7]+0.5000000000000001*f_xx[2]*nuVtSqSum[6]+0.5*f_xx[0]*nuVtSqSum[4]+0.4472135954999579*nuVtSqSum[3]*f_xx[4]+0.4472135954999579*f_xx[1]*nuVtSqSum[1]; 
  incr[8] = 0.31943828249996997*nuVtSqSum[8]*f_xx[20]+0.5*nuVtSqSum[4]*f_xx[20]+0.31943828249996997*nuVtSqSum[7]*f_xx[12]+0.5000000000000001*nuVtSqSum[1]*f_xx[12]+0.4472135954999579*nuVtSqSum[6]*f_xx[11]+0.5*f_xx[7]*nuVtSqSum[8]+0.31943828249996997*nuVtSqSum[5]*f_xx[8]+0.5*nuVtSqSum[0]*f_xx[8]+0.5000000000000001*f_xx[1]*nuVtSqSum[7]+0.5*f_xx[0]*nuVtSqSum[5]+0.4472135954999579*nuVtSqSum[3]*f_xx[4]+0.4472135954999579*f_xx[2]*nuVtSqSum[2]; 
  incr[9] = 0.5*nuVtSqSum[8]*f_xx[26]+0.5000000000000001*nuVtSqSum[7]*f_xx[25]+0.5000000000000001*nuVtSqSum[6]*f_xx[24]+0.5*nuVtSqSum[5]*f_xx[22]+0.5*nuVtSqSum[4]*f_xx[21]+0.5*nuVtSqSum[3]*f_xx[19]+0.5000000000000001*nuVtSqSum[2]*f_xx[16]+0.5000000000000001*nuVtSqSum[1]*f_xx[15]+0.5*nuVtSqSum[0]*f_xx[9]; 
  incr[10] = 0.4*nuVtSqSum[3]*f_xx[23]+0.4*nuVtSqSum[6]*f_xx[18]+0.4472135954999579*nuVtSqSum[2]*f_xx[18]+0.4*nuVtSqSum[7]*f_xx[17]+0.4472135954999579*nuVtSqSum[1]*f_xx[17]+0.44721359549995804*nuVtSqSum[3]*f_xx[14]+0.44721359549995804*nuVtSqSum[3]*f_xx[13]+0.4*nuVtSqSum[8]*f_xx[10]+0.4472135954999579*nuVtSqSum[5]*f_xx[10]+0.4472135954999579*nuVtSqSum[4]*f_xx[10]+0.5*nuVtSqSum[0]*f_xx[10]+0.44721359549995804*f_xx[6]*nuVtSqSum[7]+0.44721359549995804*f_xx[5]*nuVtSqSum[6]+0.5*nuVtSqSum[1]*f_xx[6]+0.5*nuVtSqSum[2]*f_xx[5]+0.5*f_xx[3]*nuVtSqSum[3]; 
  incr[11] = 0.2857142857142857*nuVtSqSum[6]*f_xx[20]+0.44721359549995804*nuVtSqSum[2]*f_xx[20]+0.4*nuVtSqSum[3]*f_xx[12]+0.2857142857142857*nuVtSqSum[8]*f_xx[11]+0.4472135954999579*nuVtSqSum[5]*f_xx[11]+0.31943828249996997*nuVtSqSum[4]*f_xx[11]+0.5*nuVtSqSum[0]*f_xx[11]+0.44721359549995804*f_xx[2]*nuVtSqSum[8]+0.4472135954999579*nuVtSqSum[6]*f_xx[8]+0.4*f_xx[4]*nuVtSqSum[7]+0.31943828249996997*nuVtSqSum[6]*f_xx[7]+0.5000000000000001*nuVtSqSum[2]*f_xx[7]+0.5*f_xx[0]*nuVtSqSum[6]+0.5000000000000001*f_xx[2]*nuVtSqSum[4]+0.44721359549995804*nuVtSqSum[1]*f_xx[4]+0.44721359549995804*f_xx[1]*nuVtSqSum[3]; 
  incr[12] = 0.2857142857142857*nuVtSqSum[7]*f_xx[20]+0.44721359549995804*nuVtSqSum[1]*f_xx[20]+0.2857142857142857*nuVtSqSum[8]*f_xx[12]+0.31943828249996997*nuVtSqSum[5]*f_xx[12]+0.4472135954999579*nuVtSqSum[4]*f_xx[12]+0.5*nuVtSqSum[0]*f_xx[12]+0.4*nuVtSqSum[3]*f_xx[11]+0.44721359549995804*f_xx[1]*nuVtSqSum[8]+0.31943828249996997*nuVtSqSum[7]*f_xx[8]+0.5000000000000001*nuVtSqSum[1]*f_xx[8]+0.4472135954999579*f_xx[7]*nuVtSqSum[7]+0.5*f_xx[0]*nuVtSqSum[7]+0.4*f_xx[4]*nuVtSqSum[6]+0.5000000000000001*f_xx[1]*nuVtSqSum[5]+0.44721359549995804*nuVtSqSum[2]*f_xx[4]+0.44721359549995804*f_xx[2]*nuVtSqSum[3]; 
  incr[13] = 0.31943828249997*nuVtSqSum[8]*f_xx[23]+0.5000000000000001*nuVtSqSum[5]*f_xx[23]+0.4472135954999579*nuVtSqSum[7]*f_xx[18]+0.31943828249996997*nuVtSqSum[6]*f_xx[17]+0.5000000000000001*nuVtSqSum[2]*f_xx[17]+0.5*nuVtSqSum[8]*f_xx[14]+0.31943828249996997*nuVtSqSum[4]*f_xx[13]+0.5*nuVtSqSum[0]*f_xx[13]+0.44721359549995804*nuVtSqSum[3]*f_xx[10]+0.5*f_xx[6]*nuVtSqSum[6]+0.44721359549995804*nuVtSqSum[1]*f_xx[5]+0.5000000000000001*f_xx[3]*nuVtSqSum[4]; 
  incr[14] = 0.31943828249997*nuVtSqSum[8]*f_xx[23]+0.5000000000000001*nuVtSqSum[4]*f_xx[23]+0.31943828249996997*nuVtSqSum[7]*f_xx[18]+0.5000000000000001*nuVtSqSum[1]*f_xx[18]+0.4472135954999579*nuVtSqSum[6]*f_xx[17]+0.31943828249996997*nuVtSqSum[5]*f_xx[14]+0.5*nuVtSqSum[0]*f_xx[14]+0.5*nuVtSqSum[8]*f_xx[13]+0.44721359549995804*nuVtSqSum[3]*f_xx[10]+0.5*f_xx[5]*nuVtSqSum[7]+0.44721359549995804*nuVtSqSum[2]*f_xx[6]+0.5000000000000001*f_xx[3]*nuVtSqSum[5]; 
  incr[15] = 0.4472135954999579*nuVtSqSum[7]*f_xx[26]+0.44721359549995804*nuVtSqSum[8]*f_xx[25]+0.5000000000000001*nuVtSqSum[5]*f_xx[25]+0.44721359549995804*nuVtSqSum[3]*f_xx[24]+0.5*nuVtSqSum[7]*f_xx[22]+0.44721359549995804*nuVtSqSum[1]*f_xx[21]+0.4472135954999579*nuVtSqSum[6]*f_xx[19]+0.5000000000000001*nuVtSqSum[2]*f_xx[19]+0.5*nuVtSqSum[3]*f_xx[16]+0.4472135954999579*nuVtSqSum[4]*f_xx[15]+0.5*nuVtSqSum[0]*f_xx[15]+0.5000000000000001*nuVtSqSum[1]*f_xx[9]; 
  incr[16] = 0.4472135954999579*nuVtSqSum[6]*f_xx[26]+0.44721359549995804*nuVtSqSum[3]*f_xx[25]+0.44721359549995804*nuVtSqSum[8]*f_xx[24]+0.5000000000000001*nuVtSqSum[4]*f_xx[24]+0.44721359549995804*nuVtSqSum[2]*f_xx[22]+0.5*nuVtSqSum[6]*f_xx[21]+0.4472135954999579*nuVtSqSum[7]*f_xx[19]+0.5000000000000001*nuVtSqSum[1]*f_xx[19]+0.4472135954999579*nuVtSqSum[5]*f_xx[16]+0.5*nuVtSqSum[0]*f_xx[16]+0.5*nuVtSqSum[3]*f_xx[15]+0.5000000000000001*nuVtSqSum[2]*f_xx[9]; 
  incr[17] = 0.28571428571428575*nuVtSqSum[6]*f_xx[23]+0.4472135954999579*nuVtSqSum[2]*f_xx[23]+0.4*nuVtSqSum[3]*f_xx[18]+0.2857142857142857*nuVtSqSum[8]*f_xx[17]+0.4472135954999579*nuVtSqSum[5]*f_xx[17]+0.31943828249996997*nuVtSqSum[4]*f_xx[17]+0.5*nuVtSqSum[0]*f_xx[17]+0.4472135954999579*nuVtSqSum[6]*f_xx[14]+0.31943828249996997*nuVtSqSum[6]*f_xx[13]+0.5000000000000001*nuVtSqSum[2]*f_xx[13]+0.4*nuVtSqSum[7]*f_xx[10]+0.4472135954999579*nuVtSqSum[1]*f_xx[10]+0.4472135954999579*f_xx[6]*nuVtSqSum[8]+0.5000000000000001*f_xx[3]*nuVtSqSum[6]+0.5*nuVtSqSum[4]*f_xx[6]+0.4472135954999579*nuVtSqSum[3]*f_xx[5]; 
  incr[18] = 0.28571428571428575*nuVtSqSum[7]*f_xx[23]+0.4472135954999579*nuVtSqSum[1]*f_xx[23]+0.2857142857142857*nuVtSqSum[8]*f_xx[18]+0.31943828249996997*nuVtSqSum[5]*f_xx[18]+0.4472135954999579*nuVtSqSum[4]*f_xx[18]+0.5*nuVtSqSum[0]*f_xx[18]+0.4*nuVtSqSum[3]*f_xx[17]+0.31943828249996997*nuVtSqSum[7]*f_xx[14]+0.5000000000000001*nuVtSqSum[1]*f_xx[14]+0.4472135954999579*nuVtSqSum[7]*f_xx[13]+0.4*nuVtSqSum[6]*f_xx[10]+0.4472135954999579*nuVtSqSum[2]*f_xx[10]+0.4472135954999579*f_xx[5]*nuVtSqSum[8]+0.5000000000000001*f_xx[3]*nuVtSqSum[7]+0.4472135954999579*nuVtSqSum[3]*f_xx[6]+0.5*f_xx[5]*nuVtSqSum[5]; 
  incr[19] = 0.4*nuVtSqSum[3]*f_xx[26]+0.4*nuVtSqSum[6]*f_xx[25]+0.4472135954999579*nuVtSqSum[2]*f_xx[25]+0.4*nuVtSqSum[7]*f_xx[24]+0.4472135954999579*nuVtSqSum[1]*f_xx[24]+0.4472135954999579*nuVtSqSum[3]*f_xx[22]+0.4472135954999579*nuVtSqSum[3]*f_xx[21]+0.4*nuVtSqSum[8]*f_xx[19]+0.4472135954999579*nuVtSqSum[5]*f_xx[19]+0.4472135954999579*nuVtSqSum[4]*f_xx[19]+0.5*nuVtSqSum[0]*f_xx[19]+0.4472135954999579*nuVtSqSum[7]*f_xx[16]+0.5000000000000001*nuVtSqSum[1]*f_xx[16]+0.4472135954999579*nuVtSqSum[6]*f_xx[15]+0.5000000000000001*nuVtSqSum[2]*f_xx[15]+0.5*nuVtSqSum[3]*f_xx[9]; 
  incr[20] = 0.20408163265306123*nuVtSqSum[8]*f_xx[20]+0.31943828249996997*nuVtSqSum[5]*f_xx[20]+0.31943828249996997*nuVtSqSum[4]*f_xx[20]+0.5*nuVtSqSum[0]*f_xx[20]+0.2857142857142857*nuVtSqSum[7]*f_xx[12]+0.44721359549995804*nuVtSqSum[1]*f_xx[12]+0.2857142857142857*nuVtSqSum[6]*f_xx[11]+0.44721359549995804*nuVtSqSum[2]*f_xx[11]+0.31943828249996997*f_xx[8]*nuVtSqSum[8]+0.31943828249996997*f_xx[7]*nuVtSqSum[8]+0.5*f_xx[0]*nuVtSqSum[8]+0.5*nuVtSqSum[4]*f_xx[8]+0.44721359549995804*f_xx[1]*nuVtSqSum[7]+0.5*nuVtSqSum[5]*f_xx[7]+0.44721359549995804*f_xx[2]*nuVtSqSum[6]+0.4*nuVtSqSum[3]*f_xx[4]; 
  incr[21] = 0.31943828249996997*nuVtSqSum[8]*f_xx[26]+0.5*nuVtSqSum[5]*f_xx[26]+0.44721359549995804*nuVtSqSum[7]*f_xx[25]+0.31943828249997*nuVtSqSum[6]*f_xx[24]+0.5*nuVtSqSum[2]*f_xx[24]+0.5*nuVtSqSum[8]*f_xx[22]+0.31943828249996997*nuVtSqSum[4]*f_xx[21]+0.5*nuVtSqSum[0]*f_xx[21]+0.4472135954999579*nuVtSqSum[3]*f_xx[19]+0.5*nuVtSqSum[6]*f_xx[16]+0.44721359549995804*nuVtSqSum[1]*f_xx[15]+0.5*nuVtSqSum[4]*f_xx[9]; 
  incr[22] = 0.31943828249996997*nuVtSqSum[8]*f_xx[26]+0.5*nuVtSqSum[4]*f_xx[26]+0.31943828249997*nuVtSqSum[7]*f_xx[25]+0.5*nuVtSqSum[1]*f_xx[25]+0.44721359549995804*nuVtSqSum[6]*f_xx[24]+0.31943828249996997*nuVtSqSum[5]*f_xx[22]+0.5*nuVtSqSum[0]*f_xx[22]+0.5*nuVtSqSum[8]*f_xx[21]+0.4472135954999579*nuVtSqSum[3]*f_xx[19]+0.44721359549995804*nuVtSqSum[2]*f_xx[16]+0.5*nuVtSqSum[7]*f_xx[15]+0.5*nuVtSqSum[5]*f_xx[9]; 
  incr[23] = 0.20408163265306123*nuVtSqSum[8]*f_xx[23]+0.31943828249996997*nuVtSqSum[5]*f_xx[23]+0.31943828249996997*nuVtSqSum[4]*f_xx[23]+0.5*nuVtSqSum[0]*f_xx[23]+0.28571428571428575*nuVtSqSum[7]*f_xx[18]+0.4472135954999579*nuVtSqSum[1]*f_xx[18]+0.28571428571428575*nuVtSqSum[6]*f_xx[17]+0.4472135954999579*nuVtSqSum[2]*f_xx[17]+0.31943828249997*nuVtSqSum[8]*f_xx[14]+0.5000000000000001*nuVtSqSum[4]*f_xx[14]+0.31943828249997*nuVtSqSum[8]*f_xx[13]+0.5000000000000001*nuVtSqSum[5]*f_xx[13]+0.4*nuVtSqSum[3]*f_xx[10]+0.5*f_xx[3]*nuVtSqSum[8]+0.44721359549995804*f_xx[5]*nuVtSqSum[7]+0.44721359549995804*f_xx[6]*nuVtSqSum[6]; 
  incr[24] = 0.28571428571428575*nuVtSqSum[6]*f_xx[26]+0.4472135954999579*nuVtSqSum[2]*f_xx[26]+0.4*nuVtSqSum[3]*f_xx[25]+0.2857142857142857*nuVtSqSum[8]*f_xx[24]+0.4472135954999579*nuVtSqSum[5]*f_xx[24]+0.31943828249996997*nuVtSqSum[4]*f_xx[24]+0.5*nuVtSqSum[0]*f_xx[24]+0.44721359549995804*nuVtSqSum[6]*f_xx[22]+0.31943828249997*nuVtSqSum[6]*f_xx[21]+0.5*nuVtSqSum[2]*f_xx[21]+0.4*nuVtSqSum[7]*f_xx[19]+0.4472135954999579*nuVtSqSum[1]*f_xx[19]+0.44721359549995804*nuVtSqSum[8]*f_xx[16]+0.5000000000000001*nuVtSqSum[4]*f_xx[16]+0.44721359549995804*nuVtSqSum[3]*f_xx[15]+0.5000000000000001*nuVtSqSum[6]*f_xx[9]; 
  incr[25] = 0.28571428571428575*nuVtSqSum[7]*f_xx[26]+0.4472135954999579*nuVtSqSum[1]*f_xx[26]+0.2857142857142857*nuVtSqSum[8]*f_xx[25]+0.31943828249996997*nuVtSqSum[5]*f_xx[25]+0.4472135954999579*nuVtSqSum[4]*f_xx[25]+0.5*nuVtSqSum[0]*f_xx[25]+0.4*nuVtSqSum[3]*f_xx[24]+0.31943828249997*nuVtSqSum[7]*f_xx[22]+0.5*nuVtSqSum[1]*f_xx[22]+0.44721359549995804*nuVtSqSum[7]*f_xx[21]+0.4*nuVtSqSum[6]*f_xx[19]+0.4472135954999579*nuVtSqSum[2]*f_xx[19]+0.44721359549995804*nuVtSqSum[3]*f_xx[16]+0.44721359549995804*nuVtSqSum[8]*f_xx[15]+0.5000000000000001*nuVtSqSum[5]*f_xx[15]+0.5000000000000001*nuVtSqSum[7]*f_xx[9]; 
  incr[26] = 0.20408163265306123*nuVtSqSum[8]*f_xx[26]+0.31943828249996997*nuVtSqSum[5]*f_xx[26]+0.31943828249996997*nuVtSqSum[4]*f_xx[26]+0.5*nuVtSqSum[0]*f_xx[26]+0.28571428571428575*nuVtSqSum[7]*f_xx[25]+0.4472135954999579*nuVtSqSum[1]*f_xx[25]+0.28571428571428575*nuVtSqSum[6]*f_xx[24]+0.4472135954999579*nuVtSqSum[2]*f_xx[24]+0.31943828249996997*nuVtSqSum[8]*f_xx[22]+0.5*nuVtSqSum[4]*f_xx[22]+0.31943828249996997*nuVtSqSum[8]*f_xx[21]+0.5*nuVtSqSum[5]*f_xx[21]+0.4*nuVtSqSum[3]*f_xx[19]+0.4472135954999579*nuVtSqSum[6]*f_xx[16]+0.4472135954999579*nuVtSqSum[7]*f_xx[15]+0.5*nuVtSqSum[8]*f_xx[9]; 

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
  out[12] += incr[12]*rdvSq4; 
  out[13] += incr[13]*rdvSq4; 
  out[14] += incr[14]*rdvSq4; 
  out[15] += incr[15]*rdvSq4; 
  out[16] += incr[16]*rdvSq4; 
  out[17] += incr[17]*rdvSq4; 
  out[18] += incr[18]*rdvSq4; 
  out[19] += incr[19]*rdvSq4; 
  out[20] += incr[20]*rdvSq4; 
  out[21] += incr[21]*rdvSq4; 
  out[22] += incr[22]*rdvSq4; 
  out[23] += incr[23]*rdvSq4; 
  out[24] += incr[24]*rdvSq4; 
  out[25] += incr[25]*rdvSq4; 
  out[26] += incr[26]*rdvSq4; 

  return 0.;

} 
