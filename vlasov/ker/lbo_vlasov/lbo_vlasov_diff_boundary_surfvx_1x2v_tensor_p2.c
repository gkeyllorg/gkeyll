#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvx_1x2v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[3]: Cell-center coordinates. 
  // dxv[3]: Cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[9]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fSkin/Edge: Distribution function in cells 
  // out: Incremented distribution function in cell 
  const double *nuVtSqSum = &nuPrimMomsSum[6];

  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 

  double facDiff[3]; 
  // Expand nuVtSqSum in phase basis.
  facDiff[0] = nuVtSqSum[0]; 
  facDiff[1] = nuVtSqSum[1]; 
  facDiff[2] = nuVtSqSum[2]; 

  double vol_incr[27] = {0.0}; 
  vol_incr[8] = 4.743416490252569*facDiff[2]*fSkin[7]+4.743416490252569*fSkin[1]*facDiff[1]+4.743416490252569*fSkin[0]*facDiff[0]; 
  vol_incr[12] = 4.242640687119286*facDiff[1]*fSkin[7]+4.242640687119286*fSkin[1]*facDiff[2]+4.743416490252569*fSkin[0]*facDiff[1]+4.743416490252569*facDiff[0]*fSkin[1]; 
  vol_incr[14] = 4.743416490252569*facDiff[2]*fSkin[13]+4.743416490252569*facDiff[1]*fSkin[5]+4.743416490252569*facDiff[0]*fSkin[3]; 
  vol_incr[18] = 4.242640687119286*facDiff[1]*fSkin[13]+4.242640687119286*facDiff[2]*fSkin[5]+4.743416490252569*facDiff[0]*fSkin[5]+4.743416490252569*facDiff[1]*fSkin[3]; 
  vol_incr[20] = 3.0304576336566322*facDiff[2]*fSkin[7]+4.743416490252569*facDiff[0]*fSkin[7]+4.743416490252569*fSkin[0]*facDiff[2]+4.242640687119286*fSkin[1]*facDiff[1]; 
  vol_incr[22] = 4.743416490252569*facDiff[2]*fSkin[21]+4.743416490252569*facDiff[1]*fSkin[15]+4.743416490252569*facDiff[0]*fSkin[9]; 
  vol_incr[23] = 3.0304576336566327*facDiff[2]*fSkin[13]+4.743416490252569*facDiff[0]*fSkin[13]+4.242640687119286*facDiff[1]*fSkin[5]+4.743416490252569*facDiff[2]*fSkin[3]; 
  vol_incr[25] = 4.242640687119286*facDiff[1]*fSkin[21]+4.242640687119286*facDiff[2]*fSkin[15]+4.743416490252569*facDiff[0]*fSkin[15]+4.743416490252569*facDiff[1]*fSkin[9]; 
  vol_incr[26] = 3.0304576336566322*facDiff[2]*fSkin[21]+4.743416490252569*facDiff[0]*fSkin[21]+4.242640687119286*facDiff[1]*fSkin[15]+4.743416490252569*facDiff[2]*fSkin[9]; 

  double edgeSurf_incr[27] = {0.0}; 
  double boundSurf_incr[27] = {0.0}; 

  if (edge == -1) { 

  double edgeSurf[27] = {0.0}; 
  edgeSurf[0] = -(0.6708203932499369*fSkin[8])+0.6708203932499369*fEdge[8]-1.190784930203603*fSkin[2]-1.190784930203603*fEdge[2]-0.9375*fSkin[0]+0.9375*fEdge[0]; 
  edgeSurf[1] = -(0.6708203932499369*fSkin[12])+0.6708203932499369*fEdge[12]-1.190784930203603*fSkin[4]-1.190784930203603*fEdge[4]-0.9375*fSkin[1]+0.9375*fEdge[1]; 
  edgeSurf[2] = -(1.5855025573536612*fSkin[8])+0.7382874503707888*fEdge[8]-2.671875*fSkin[2]-1.453125*fEdge[2]-2.0568103339880417*fSkin[0]+1.190784930203603*fEdge[0]; 
  edgeSurf[3] = -(0.6708203932499369*fSkin[14])+0.6708203932499369*fEdge[14]-1.190784930203603*fSkin[6]-1.190784930203603*fEdge[6]-0.9375*fSkin[3]+0.9375*fEdge[3]; 
  edgeSurf[4] = -(1.5855025573536612*fSkin[12])+0.7382874503707888*fEdge[12]-2.671875*fSkin[4]-1.453125*fEdge[4]-2.0568103339880417*fSkin[1]+1.190784930203603*fEdge[1]; 
  edgeSurf[5] = -(0.6708203932499369*fSkin[18])+0.6708203932499369*fEdge[18]-1.190784930203603*fSkin[10]-1.190784930203603*fEdge[10]-0.9375*fSkin[5]+0.9375*fEdge[5]; 
  edgeSurf[6] = -(1.5855025573536612*fSkin[14])+0.7382874503707888*fEdge[14]-2.671875*fSkin[6]-1.453125*fEdge[6]-2.0568103339880417*fSkin[3]+1.190784930203603*fEdge[3]; 
  edgeSurf[7] = -(0.6708203932499369*fSkin[20])+0.6708203932499369*fEdge[20]-1.190784930203603*fSkin[11]-1.190784930203603*fEdge[11]-0.9375*fSkin[7]+0.9375*fEdge[7]; 
  edgeSurf[8] = -(3.140625*fSkin[8])-0.140625*fEdge[8]-5.022775277112744*fSkin[2]-0.3025768239224545*fEdge[2]-3.7733647120308955*fSkin[0]+0.4192627457812106*fEdge[0]; 
  edgeSurf[9] = -(0.6708203932499369*fSkin[22])+0.6708203932499369*fEdge[22]-1.190784930203603*fSkin[16]-1.190784930203603*fEdge[16]-0.9375*fSkin[9]+0.9375*fEdge[9]; 
  edgeSurf[10] = -(1.5855025573536612*fSkin[18])+0.7382874503707888*fEdge[18]-2.671875*fSkin[10]-1.453125*fEdge[10]-2.0568103339880417*fSkin[5]+1.190784930203603*fEdge[5]; 
  edgeSurf[11] = -(1.5855025573536612*fSkin[20])+0.7382874503707888*fEdge[20]-2.671875*fSkin[11]-1.453125*fEdge[11]-2.0568103339880417*fSkin[7]+1.190784930203603*fEdge[7]; 
  edgeSurf[12] = -(3.140625*fSkin[12])-0.140625*fEdge[12]-5.022775277112744*fSkin[4]-0.30257682392245444*fEdge[4]-3.773364712030894*fSkin[1]+0.41926274578121053*fEdge[1]; 
  edgeSurf[13] = -(0.6708203932499369*fSkin[23])+0.6708203932499369*fEdge[23]-1.190784930203603*fSkin[17]-1.190784930203603*fEdge[17]-0.9375*fSkin[13]+0.9375*fEdge[13]; 
  edgeSurf[14] = -(3.140625*fSkin[14])-0.140625*fEdge[14]-5.022775277112744*fSkin[6]-0.30257682392245444*fEdge[6]-3.773364712030894*fSkin[3]+0.41926274578121053*fEdge[3]; 
  edgeSurf[15] = -(0.6708203932499369*fSkin[25])+0.6708203932499369*fEdge[25]-1.190784930203603*fSkin[19]-1.190784930203603*fEdge[19]-0.9375*fSkin[15]+0.9375*fEdge[15]; 
  edgeSurf[16] = -(1.5855025573536612*fSkin[22])+0.7382874503707888*fEdge[22]-2.671875*fSkin[16]-1.453125*fEdge[16]-2.0568103339880417*fSkin[9]+1.190784930203603*fEdge[9]; 
  edgeSurf[17] = -(1.5855025573536612*fSkin[23])+0.7382874503707888*fEdge[23]-2.671875*fSkin[17]-1.453125*fEdge[17]-2.0568103339880417*fSkin[13]+1.190784930203603*fEdge[13]; 
  edgeSurf[18] = -(3.140625*fSkin[18])-0.140625*fEdge[18]-5.022775277112744*fSkin[10]-0.3025768239224545*fEdge[10]-3.7733647120308955*fSkin[5]+0.4192627457812106*fEdge[5]; 
  edgeSurf[19] = -(1.5855025573536612*fSkin[25])+0.7382874503707888*fEdge[25]-2.671875*fSkin[19]-1.453125*fEdge[19]-2.0568103339880417*fSkin[15]+1.190784930203603*fEdge[15]; 
  edgeSurf[20] = -(3.140625*fSkin[20])-0.140625*fEdge[20]-5.022775277112744*fSkin[11]-0.30257682392245444*fEdge[11]-3.7733647120308955*fSkin[7]+0.4192627457812106*fEdge[7]; 
  edgeSurf[21] = -(0.6708203932499369*fSkin[26])+0.6708203932499369*fEdge[26]-1.190784930203603*fSkin[24]-1.190784930203603*fEdge[24]-0.9375*fSkin[21]+0.9375*fEdge[21]; 
  edgeSurf[22] = -(3.140625*fSkin[22])-0.140625*fEdge[22]-5.022775277112744*fSkin[16]-0.30257682392245444*fEdge[16]-3.7733647120308955*fSkin[9]+0.4192627457812106*fEdge[9]; 
  edgeSurf[23] = -(3.140625*fSkin[23])-0.140625*fEdge[23]-5.022775277112744*fSkin[17]-0.3025768239224545*fEdge[17]-3.773364712030894*fSkin[13]+0.41926274578121053*fEdge[13]; 
  edgeSurf[24] = -(1.5855025573536612*fSkin[26])+0.7382874503707888*fEdge[26]-2.671875*fSkin[24]-1.453125*fEdge[24]-2.0568103339880417*fSkin[21]+1.190784930203603*fEdge[21]; 
  edgeSurf[25] = -(3.140625*fSkin[25])-0.140625*fEdge[25]-5.022775277112744*fSkin[19]-0.3025768239224545*fEdge[19]-3.773364712030894*fSkin[15]+0.41926274578121053*fEdge[15]; 
  edgeSurf[26] = -(3.140625*fSkin[26])-0.140625*fEdge[26]-5.022775277112744*fSkin[24]-0.3025768239224545*fEdge[24]-3.7733647120308955*fSkin[21]+0.4192627457812106*fEdge[21]; 

  double boundSurf[27] = {0.0}; 
  boundSurf[2] = 1.9364916731037085*fSkin[8]-1.5*fSkin[2]+0.8660254037844386*fSkin[0]; 
  boundSurf[4] = 1.9364916731037085*fSkin[12]-1.5*fSkin[4]+0.8660254037844386*fSkin[1]; 
  boundSurf[6] = 1.9364916731037085*fSkin[14]-1.5*fSkin[6]+0.8660254037844386*fSkin[3]; 
  boundSurf[8] = -(7.5*fSkin[8])+5.809475019311125*fSkin[2]-3.3541019662496847*fSkin[0]; 
  boundSurf[10] = 1.9364916731037085*fSkin[18]-1.5*fSkin[10]+0.8660254037844386*fSkin[5]; 
  boundSurf[11] = 1.9364916731037085*fSkin[20]-1.5*fSkin[11]+0.8660254037844387*fSkin[7]; 
  boundSurf[12] = -(7.5*fSkin[12])+5.809475019311126*fSkin[4]-3.3541019662496843*fSkin[1]; 
  boundSurf[14] = -(7.5*fSkin[14])+5.809475019311126*fSkin[6]-3.3541019662496843*fSkin[3]; 
  boundSurf[16] = 1.9364916731037085*fSkin[22]-1.5*fSkin[16]+0.8660254037844387*fSkin[9]; 
  boundSurf[17] = 1.9364916731037085*fSkin[23]-1.5*fSkin[17]+0.8660254037844387*fSkin[13]; 
  boundSurf[18] = -(7.5*fSkin[18])+5.809475019311125*fSkin[10]-3.3541019662496847*fSkin[5]; 
  boundSurf[19] = 1.9364916731037085*fSkin[25]-1.5*fSkin[19]+0.8660254037844387*fSkin[15]; 
  boundSurf[20] = -(7.5*fSkin[20])+5.809475019311126*fSkin[11]-3.3541019662496847*fSkin[7]; 
  boundSurf[22] = -(7.5*fSkin[22])+5.809475019311126*fSkin[16]-3.3541019662496847*fSkin[9]; 
  boundSurf[23] = -(7.5*fSkin[23])+5.809475019311125*fSkin[17]-3.3541019662496843*fSkin[13]; 
  boundSurf[24] = 1.9364916731037085*fSkin[26]-1.5*fSkin[24]+0.8660254037844386*fSkin[21]; 
  boundSurf[25] = -(7.5*fSkin[25])+5.809475019311125*fSkin[19]-3.3541019662496843*fSkin[15]; 
  boundSurf[26] = -(7.5*fSkin[26])+5.809475019311125*fSkin[24]-3.3541019662496847*fSkin[21]; 

  edgeSurf_incr[0] = 0.7071067811865475*nuVtSqSum[2]*edgeSurf[7]+0.7071067811865475*edgeSurf[1]*nuVtSqSum[1]+0.7071067811865475*edgeSurf[0]*nuVtSqSum[0]; 
  edgeSurf_incr[1] = 0.6324555320336759*nuVtSqSum[1]*edgeSurf[7]+0.6324555320336759*edgeSurf[1]*nuVtSqSum[2]+0.7071067811865475*edgeSurf[0]*nuVtSqSum[1]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[1]; 
  edgeSurf_incr[2] = 0.7071067811865475*nuVtSqSum[2]*edgeSurf[11]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[4]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[2]; 
  edgeSurf_incr[3] = 0.7071067811865475*nuVtSqSum[2]*edgeSurf[13]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[5]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[3]; 
  edgeSurf_incr[4] = 0.632455532033676*nuVtSqSum[1]*edgeSurf[11]+0.6324555320336759*nuVtSqSum[2]*edgeSurf[4]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[4]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[2]; 
  edgeSurf_incr[5] = 0.632455532033676*nuVtSqSum[1]*edgeSurf[13]+0.6324555320336759*nuVtSqSum[2]*edgeSurf[5]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[5]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[3]; 
  edgeSurf_incr[6] = 0.7071067811865475*nuVtSqSum[2]*edgeSurf[17]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[10]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[6]; 
  edgeSurf_incr[7] = 0.45175395145262565*nuVtSqSum[2]*edgeSurf[7]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[7]+0.7071067811865475*edgeSurf[0]*nuVtSqSum[2]+0.6324555320336759*edgeSurf[1]*nuVtSqSum[1]; 
  edgeSurf_incr[8] = 0.7071067811865475*nuVtSqSum[2]*edgeSurf[20]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[12]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[8]; 
  edgeSurf_incr[9] = 0.7071067811865475*nuVtSqSum[2]*edgeSurf[21]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[15]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[9]; 
  edgeSurf_incr[10] = 0.6324555320336759*nuVtSqSum[1]*edgeSurf[17]+0.6324555320336759*nuVtSqSum[2]*edgeSurf[10]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[10]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[6]; 
  edgeSurf_incr[11] = 0.45175395145262565*nuVtSqSum[2]*edgeSurf[11]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[11]+0.632455532033676*nuVtSqSum[1]*edgeSurf[4]+0.7071067811865475*edgeSurf[2]*nuVtSqSum[2]; 
  edgeSurf_incr[12] = 0.632455532033676*nuVtSqSum[1]*edgeSurf[20]+0.6324555320336759*nuVtSqSum[2]*edgeSurf[12]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[12]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[8]; 
  edgeSurf_incr[13] = 0.45175395145262565*nuVtSqSum[2]*edgeSurf[13]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[13]+0.632455532033676*nuVtSqSum[1]*edgeSurf[5]+0.7071067811865475*nuVtSqSum[2]*edgeSurf[3]; 
  edgeSurf_incr[14] = 0.7071067811865475*nuVtSqSum[2]*edgeSurf[23]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[18]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[14]; 
  edgeSurf_incr[15] = 0.632455532033676*nuVtSqSum[1]*edgeSurf[21]+0.6324555320336759*nuVtSqSum[2]*edgeSurf[15]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[15]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[9]; 
  edgeSurf_incr[16] = 0.7071067811865475*nuVtSqSum[2]*edgeSurf[24]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[19]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[16]; 
  edgeSurf_incr[17] = 0.45175395145262565*nuVtSqSum[2]*edgeSurf[17]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[17]+0.6324555320336759*nuVtSqSum[1]*edgeSurf[10]+0.7071067811865475*nuVtSqSum[2]*edgeSurf[6]; 
  edgeSurf_incr[18] = 0.6324555320336759*nuVtSqSum[1]*edgeSurf[23]+0.6324555320336759*nuVtSqSum[2]*edgeSurf[18]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[18]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[14]; 
  edgeSurf_incr[19] = 0.6324555320336759*nuVtSqSum[1]*edgeSurf[24]+0.6324555320336759*nuVtSqSum[2]*edgeSurf[19]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[19]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[16]; 
  edgeSurf_incr[20] = 0.45175395145262565*nuVtSqSum[2]*edgeSurf[20]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[20]+0.632455532033676*nuVtSqSum[1]*edgeSurf[12]+0.7071067811865475*nuVtSqSum[2]*edgeSurf[8]; 
  edgeSurf_incr[21] = 0.45175395145262565*nuVtSqSum[2]*edgeSurf[21]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[21]+0.632455532033676*nuVtSqSum[1]*edgeSurf[15]+0.7071067811865475*nuVtSqSum[2]*edgeSurf[9]; 
  edgeSurf_incr[22] = 0.7071067811865475*nuVtSqSum[2]*edgeSurf[26]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[25]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[22]; 
  edgeSurf_incr[23] = 0.45175395145262565*nuVtSqSum[2]*edgeSurf[23]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[23]+0.6324555320336759*nuVtSqSum[1]*edgeSurf[18]+0.7071067811865475*nuVtSqSum[2]*edgeSurf[14]; 
  edgeSurf_incr[24] = 0.45175395145262565*nuVtSqSum[2]*edgeSurf[24]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[24]+0.6324555320336759*nuVtSqSum[1]*edgeSurf[19]+0.7071067811865475*nuVtSqSum[2]*edgeSurf[16]; 
  edgeSurf_incr[25] = 0.6324555320336759*nuVtSqSum[1]*edgeSurf[26]+0.6324555320336759*nuVtSqSum[2]*edgeSurf[25]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[25]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[22]; 
  edgeSurf_incr[26] = 0.45175395145262565*nuVtSqSum[2]*edgeSurf[26]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[26]+0.6324555320336759*nuVtSqSum[1]*edgeSurf[25]+0.7071067811865475*nuVtSqSum[2]*edgeSurf[22]; 

  boundSurf_incr[0] = 0.7071067811865475*nuVtSqSum[2]*boundSurf[7]+0.7071067811865475*boundSurf[1]*nuVtSqSum[1]+0.7071067811865475*boundSurf[0]*nuVtSqSum[0]; 
  boundSurf_incr[1] = 0.6324555320336759*nuVtSqSum[1]*boundSurf[7]+0.6324555320336759*boundSurf[1]*nuVtSqSum[2]+0.7071067811865475*boundSurf[0]*nuVtSqSum[1]+0.7071067811865475*nuVtSqSum[0]*boundSurf[1]; 
  boundSurf_incr[2] = 0.7071067811865475*nuVtSqSum[2]*boundSurf[11]+0.7071067811865475*nuVtSqSum[1]*boundSurf[4]+0.7071067811865475*nuVtSqSum[0]*boundSurf[2]; 
  boundSurf_incr[3] = 0.7071067811865475*nuVtSqSum[2]*boundSurf[13]+0.7071067811865475*nuVtSqSum[1]*boundSurf[5]+0.7071067811865475*nuVtSqSum[0]*boundSurf[3]; 
  boundSurf_incr[4] = 0.632455532033676*nuVtSqSum[1]*boundSurf[11]+0.6324555320336759*nuVtSqSum[2]*boundSurf[4]+0.7071067811865475*nuVtSqSum[0]*boundSurf[4]+0.7071067811865475*nuVtSqSum[1]*boundSurf[2]; 
  boundSurf_incr[5] = 0.632455532033676*nuVtSqSum[1]*boundSurf[13]+0.6324555320336759*nuVtSqSum[2]*boundSurf[5]+0.7071067811865475*nuVtSqSum[0]*boundSurf[5]+0.7071067811865475*nuVtSqSum[1]*boundSurf[3]; 
  boundSurf_incr[6] = 0.7071067811865475*nuVtSqSum[2]*boundSurf[17]+0.7071067811865475*nuVtSqSum[1]*boundSurf[10]+0.7071067811865475*nuVtSqSum[0]*boundSurf[6]; 
  boundSurf_incr[7] = 0.45175395145262565*nuVtSqSum[2]*boundSurf[7]+0.7071067811865475*nuVtSqSum[0]*boundSurf[7]+0.7071067811865475*boundSurf[0]*nuVtSqSum[2]+0.6324555320336759*boundSurf[1]*nuVtSqSum[1]; 
  boundSurf_incr[8] = 0.7071067811865475*nuVtSqSum[2]*boundSurf[20]+0.7071067811865475*nuVtSqSum[1]*boundSurf[12]+0.7071067811865475*nuVtSqSum[0]*boundSurf[8]; 
  boundSurf_incr[9] = 0.7071067811865475*nuVtSqSum[2]*boundSurf[21]+0.7071067811865475*nuVtSqSum[1]*boundSurf[15]+0.7071067811865475*nuVtSqSum[0]*boundSurf[9]; 
  boundSurf_incr[10] = 0.6324555320336759*nuVtSqSum[1]*boundSurf[17]+0.6324555320336759*nuVtSqSum[2]*boundSurf[10]+0.7071067811865475*nuVtSqSum[0]*boundSurf[10]+0.7071067811865475*nuVtSqSum[1]*boundSurf[6]; 
  boundSurf_incr[11] = 0.45175395145262565*nuVtSqSum[2]*boundSurf[11]+0.7071067811865475*nuVtSqSum[0]*boundSurf[11]+0.632455532033676*nuVtSqSum[1]*boundSurf[4]+0.7071067811865475*boundSurf[2]*nuVtSqSum[2]; 
  boundSurf_incr[12] = 0.632455532033676*nuVtSqSum[1]*boundSurf[20]+0.6324555320336759*nuVtSqSum[2]*boundSurf[12]+0.7071067811865475*nuVtSqSum[0]*boundSurf[12]+0.7071067811865475*nuVtSqSum[1]*boundSurf[8]; 
  boundSurf_incr[13] = 0.45175395145262565*nuVtSqSum[2]*boundSurf[13]+0.7071067811865475*nuVtSqSum[0]*boundSurf[13]+0.632455532033676*nuVtSqSum[1]*boundSurf[5]+0.7071067811865475*nuVtSqSum[2]*boundSurf[3]; 
  boundSurf_incr[14] = 0.7071067811865475*nuVtSqSum[2]*boundSurf[23]+0.7071067811865475*nuVtSqSum[1]*boundSurf[18]+0.7071067811865475*nuVtSqSum[0]*boundSurf[14]; 
  boundSurf_incr[15] = 0.632455532033676*nuVtSqSum[1]*boundSurf[21]+0.6324555320336759*nuVtSqSum[2]*boundSurf[15]+0.7071067811865475*nuVtSqSum[0]*boundSurf[15]+0.7071067811865475*nuVtSqSum[1]*boundSurf[9]; 
  boundSurf_incr[16] = 0.7071067811865475*nuVtSqSum[2]*boundSurf[24]+0.7071067811865475*nuVtSqSum[1]*boundSurf[19]+0.7071067811865475*nuVtSqSum[0]*boundSurf[16]; 
  boundSurf_incr[17] = 0.45175395145262565*nuVtSqSum[2]*boundSurf[17]+0.7071067811865475*nuVtSqSum[0]*boundSurf[17]+0.6324555320336759*nuVtSqSum[1]*boundSurf[10]+0.7071067811865475*nuVtSqSum[2]*boundSurf[6]; 
  boundSurf_incr[18] = 0.6324555320336759*nuVtSqSum[1]*boundSurf[23]+0.6324555320336759*nuVtSqSum[2]*boundSurf[18]+0.7071067811865475*nuVtSqSum[0]*boundSurf[18]+0.7071067811865475*nuVtSqSum[1]*boundSurf[14]; 
  boundSurf_incr[19] = 0.6324555320336759*nuVtSqSum[1]*boundSurf[24]+0.6324555320336759*nuVtSqSum[2]*boundSurf[19]+0.7071067811865475*nuVtSqSum[0]*boundSurf[19]+0.7071067811865475*nuVtSqSum[1]*boundSurf[16]; 
  boundSurf_incr[20] = 0.45175395145262565*nuVtSqSum[2]*boundSurf[20]+0.7071067811865475*nuVtSqSum[0]*boundSurf[20]+0.632455532033676*nuVtSqSum[1]*boundSurf[12]+0.7071067811865475*nuVtSqSum[2]*boundSurf[8]; 
  boundSurf_incr[21] = 0.45175395145262565*nuVtSqSum[2]*boundSurf[21]+0.7071067811865475*nuVtSqSum[0]*boundSurf[21]+0.632455532033676*nuVtSqSum[1]*boundSurf[15]+0.7071067811865475*nuVtSqSum[2]*boundSurf[9]; 
  boundSurf_incr[22] = 0.7071067811865475*nuVtSqSum[2]*boundSurf[26]+0.7071067811865475*nuVtSqSum[1]*boundSurf[25]+0.7071067811865475*nuVtSqSum[0]*boundSurf[22]; 
  boundSurf_incr[23] = 0.45175395145262565*nuVtSqSum[2]*boundSurf[23]+0.7071067811865475*nuVtSqSum[0]*boundSurf[23]+0.6324555320336759*nuVtSqSum[1]*boundSurf[18]+0.7071067811865475*nuVtSqSum[2]*boundSurf[14]; 
  boundSurf_incr[24] = 0.45175395145262565*nuVtSqSum[2]*boundSurf[24]+0.7071067811865475*nuVtSqSum[0]*boundSurf[24]+0.6324555320336759*nuVtSqSum[1]*boundSurf[19]+0.7071067811865475*nuVtSqSum[2]*boundSurf[16]; 
  boundSurf_incr[25] = 0.6324555320336759*nuVtSqSum[1]*boundSurf[26]+0.6324555320336759*nuVtSqSum[2]*boundSurf[25]+0.7071067811865475*nuVtSqSum[0]*boundSurf[25]+0.7071067811865475*nuVtSqSum[1]*boundSurf[22]; 
  boundSurf_incr[26] = 0.45175395145262565*nuVtSqSum[2]*boundSurf[26]+0.7071067811865475*nuVtSqSum[0]*boundSurf[26]+0.6324555320336759*nuVtSqSum[1]*boundSurf[25]+0.7071067811865475*nuVtSqSum[2]*boundSurf[22]; 


  } else { 

  double edgeSurf[27] = {0.0}; 
  edgeSurf[0] = -(0.6708203932499369*fSkin[8])+0.6708203932499369*fEdge[8]+1.190784930203603*fSkin[2]+1.190784930203603*fEdge[2]-0.9375*fSkin[0]+0.9375*fEdge[0]; 
  edgeSurf[1] = -(0.6708203932499369*fSkin[12])+0.6708203932499369*fEdge[12]+1.190784930203603*fSkin[4]+1.190784930203603*fEdge[4]-0.9375*fSkin[1]+0.9375*fEdge[1]; 
  edgeSurf[2] = 1.5855025573536612*fSkin[8]-0.7382874503707888*fEdge[8]-2.671875*fSkin[2]-1.453125*fEdge[2]+2.0568103339880417*fSkin[0]-1.190784930203603*fEdge[0]; 
  edgeSurf[3] = -(0.6708203932499369*fSkin[14])+0.6708203932499369*fEdge[14]+1.190784930203603*fSkin[6]+1.190784930203603*fEdge[6]-0.9375*fSkin[3]+0.9375*fEdge[3]; 
  edgeSurf[4] = 1.5855025573536612*fSkin[12]-0.7382874503707888*fEdge[12]-2.671875*fSkin[4]-1.453125*fEdge[4]+2.0568103339880417*fSkin[1]-1.190784930203603*fEdge[1]; 
  edgeSurf[5] = -(0.6708203932499369*fSkin[18])+0.6708203932499369*fEdge[18]+1.190784930203603*fSkin[10]+1.190784930203603*fEdge[10]-0.9375*fSkin[5]+0.9375*fEdge[5]; 
  edgeSurf[6] = 1.5855025573536612*fSkin[14]-0.7382874503707888*fEdge[14]-2.671875*fSkin[6]-1.453125*fEdge[6]+2.0568103339880417*fSkin[3]-1.190784930203603*fEdge[3]; 
  edgeSurf[7] = -(0.6708203932499369*fSkin[20])+0.6708203932499369*fEdge[20]+1.190784930203603*fSkin[11]+1.190784930203603*fEdge[11]-0.9375*fSkin[7]+0.9375*fEdge[7]; 
  edgeSurf[8] = -(3.140625*fSkin[8])-0.140625*fEdge[8]+5.022775277112744*fSkin[2]+0.3025768239224545*fEdge[2]-3.7733647120308955*fSkin[0]+0.4192627457812106*fEdge[0]; 
  edgeSurf[9] = -(0.6708203932499369*fSkin[22])+0.6708203932499369*fEdge[22]+1.190784930203603*fSkin[16]+1.190784930203603*fEdge[16]-0.9375*fSkin[9]+0.9375*fEdge[9]; 
  edgeSurf[10] = 1.5855025573536612*fSkin[18]-0.7382874503707888*fEdge[18]-2.671875*fSkin[10]-1.453125*fEdge[10]+2.0568103339880417*fSkin[5]-1.190784930203603*fEdge[5]; 
  edgeSurf[11] = 1.5855025573536612*fSkin[20]-0.7382874503707888*fEdge[20]-2.671875*fSkin[11]-1.453125*fEdge[11]+2.0568103339880417*fSkin[7]-1.190784930203603*fEdge[7]; 
  edgeSurf[12] = -(3.140625*fSkin[12])-0.140625*fEdge[12]+5.022775277112744*fSkin[4]+0.30257682392245444*fEdge[4]-3.773364712030894*fSkin[1]+0.41926274578121053*fEdge[1]; 
  edgeSurf[13] = -(0.6708203932499369*fSkin[23])+0.6708203932499369*fEdge[23]+1.190784930203603*fSkin[17]+1.190784930203603*fEdge[17]-0.9375*fSkin[13]+0.9375*fEdge[13]; 
  edgeSurf[14] = -(3.140625*fSkin[14])-0.140625*fEdge[14]+5.022775277112744*fSkin[6]+0.30257682392245444*fEdge[6]-3.773364712030894*fSkin[3]+0.41926274578121053*fEdge[3]; 
  edgeSurf[15] = -(0.6708203932499369*fSkin[25])+0.6708203932499369*fEdge[25]+1.190784930203603*fSkin[19]+1.190784930203603*fEdge[19]-0.9375*fSkin[15]+0.9375*fEdge[15]; 
  edgeSurf[16] = 1.5855025573536612*fSkin[22]-0.7382874503707888*fEdge[22]-2.671875*fSkin[16]-1.453125*fEdge[16]+2.0568103339880417*fSkin[9]-1.190784930203603*fEdge[9]; 
  edgeSurf[17] = 1.5855025573536612*fSkin[23]-0.7382874503707888*fEdge[23]-2.671875*fSkin[17]-1.453125*fEdge[17]+2.0568103339880417*fSkin[13]-1.190784930203603*fEdge[13]; 
  edgeSurf[18] = -(3.140625*fSkin[18])-0.140625*fEdge[18]+5.022775277112744*fSkin[10]+0.3025768239224545*fEdge[10]-3.7733647120308955*fSkin[5]+0.4192627457812106*fEdge[5]; 
  edgeSurf[19] = 1.5855025573536612*fSkin[25]-0.7382874503707888*fEdge[25]-2.671875*fSkin[19]-1.453125*fEdge[19]+2.0568103339880417*fSkin[15]-1.190784930203603*fEdge[15]; 
  edgeSurf[20] = -(3.140625*fSkin[20])-0.140625*fEdge[20]+5.022775277112744*fSkin[11]+0.30257682392245444*fEdge[11]-3.7733647120308955*fSkin[7]+0.4192627457812106*fEdge[7]; 
  edgeSurf[21] = -(0.6708203932499369*fSkin[26])+0.6708203932499369*fEdge[26]+1.190784930203603*fSkin[24]+1.190784930203603*fEdge[24]-0.9375*fSkin[21]+0.9375*fEdge[21]; 
  edgeSurf[22] = -(3.140625*fSkin[22])-0.140625*fEdge[22]+5.022775277112744*fSkin[16]+0.30257682392245444*fEdge[16]-3.7733647120308955*fSkin[9]+0.4192627457812106*fEdge[9]; 
  edgeSurf[23] = -(3.140625*fSkin[23])-0.140625*fEdge[23]+5.022775277112744*fSkin[17]+0.3025768239224545*fEdge[17]-3.773364712030894*fSkin[13]+0.41926274578121053*fEdge[13]; 
  edgeSurf[24] = 1.5855025573536612*fSkin[26]-0.7382874503707888*fEdge[26]-2.671875*fSkin[24]-1.453125*fEdge[24]+2.0568103339880417*fSkin[21]-1.190784930203603*fEdge[21]; 
  edgeSurf[25] = -(3.140625*fSkin[25])-0.140625*fEdge[25]+5.022775277112744*fSkin[19]+0.3025768239224545*fEdge[19]-3.773364712030894*fSkin[15]+0.41926274578121053*fEdge[15]; 
  edgeSurf[26] = -(3.140625*fSkin[26])-0.140625*fEdge[26]+5.022775277112744*fSkin[24]+0.3025768239224545*fEdge[24]-3.7733647120308955*fSkin[21]+0.4192627457812106*fEdge[21]; 

  double boundSurf[27] = {0.0}; 
  boundSurf[2] = -(1.9364916731037085*fSkin[8])-1.5*fSkin[2]-0.8660254037844386*fSkin[0]; 
  boundSurf[4] = -(1.9364916731037085*fSkin[12])-1.5*fSkin[4]-0.8660254037844386*fSkin[1]; 
  boundSurf[6] = -(1.9364916731037085*fSkin[14])-1.5*fSkin[6]-0.8660254037844386*fSkin[3]; 
  boundSurf[8] = -(7.5*fSkin[8])-5.809475019311125*fSkin[2]-3.3541019662496847*fSkin[0]; 
  boundSurf[10] = -(1.9364916731037085*fSkin[18])-1.5*fSkin[10]-0.8660254037844386*fSkin[5]; 
  boundSurf[11] = -(1.9364916731037085*fSkin[20])-1.5*fSkin[11]-0.8660254037844387*fSkin[7]; 
  boundSurf[12] = -(7.5*fSkin[12])-5.809475019311126*fSkin[4]-3.3541019662496843*fSkin[1]; 
  boundSurf[14] = -(7.5*fSkin[14])-5.809475019311126*fSkin[6]-3.3541019662496843*fSkin[3]; 
  boundSurf[16] = -(1.9364916731037085*fSkin[22])-1.5*fSkin[16]-0.8660254037844387*fSkin[9]; 
  boundSurf[17] = -(1.9364916731037085*fSkin[23])-1.5*fSkin[17]-0.8660254037844387*fSkin[13]; 
  boundSurf[18] = -(7.5*fSkin[18])-5.809475019311125*fSkin[10]-3.3541019662496847*fSkin[5]; 
  boundSurf[19] = -(1.9364916731037085*fSkin[25])-1.5*fSkin[19]-0.8660254037844387*fSkin[15]; 
  boundSurf[20] = -(7.5*fSkin[20])-5.809475019311126*fSkin[11]-3.3541019662496847*fSkin[7]; 
  boundSurf[22] = -(7.5*fSkin[22])-5.809475019311126*fSkin[16]-3.3541019662496847*fSkin[9]; 
  boundSurf[23] = -(7.5*fSkin[23])-5.809475019311125*fSkin[17]-3.3541019662496843*fSkin[13]; 
  boundSurf[24] = -(1.9364916731037085*fSkin[26])-1.5*fSkin[24]-0.8660254037844386*fSkin[21]; 
  boundSurf[25] = -(7.5*fSkin[25])-5.809475019311125*fSkin[19]-3.3541019662496843*fSkin[15]; 
  boundSurf[26] = -(7.5*fSkin[26])-5.809475019311125*fSkin[24]-3.3541019662496847*fSkin[21]; 

  edgeSurf_incr[0] = 0.7071067811865475*nuVtSqSum[2]*edgeSurf[7]+0.7071067811865475*edgeSurf[1]*nuVtSqSum[1]+0.7071067811865475*edgeSurf[0]*nuVtSqSum[0]; 
  edgeSurf_incr[1] = 0.6324555320336759*nuVtSqSum[1]*edgeSurf[7]+0.6324555320336759*edgeSurf[1]*nuVtSqSum[2]+0.7071067811865475*edgeSurf[0]*nuVtSqSum[1]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[1]; 
  edgeSurf_incr[2] = 0.7071067811865475*nuVtSqSum[2]*edgeSurf[11]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[4]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[2]; 
  edgeSurf_incr[3] = 0.7071067811865475*nuVtSqSum[2]*edgeSurf[13]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[5]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[3]; 
  edgeSurf_incr[4] = 0.632455532033676*nuVtSqSum[1]*edgeSurf[11]+0.6324555320336759*nuVtSqSum[2]*edgeSurf[4]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[4]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[2]; 
  edgeSurf_incr[5] = 0.632455532033676*nuVtSqSum[1]*edgeSurf[13]+0.6324555320336759*nuVtSqSum[2]*edgeSurf[5]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[5]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[3]; 
  edgeSurf_incr[6] = 0.7071067811865475*nuVtSqSum[2]*edgeSurf[17]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[10]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[6]; 
  edgeSurf_incr[7] = 0.45175395145262565*nuVtSqSum[2]*edgeSurf[7]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[7]+0.7071067811865475*edgeSurf[0]*nuVtSqSum[2]+0.6324555320336759*edgeSurf[1]*nuVtSqSum[1]; 
  edgeSurf_incr[8] = 0.7071067811865475*nuVtSqSum[2]*edgeSurf[20]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[12]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[8]; 
  edgeSurf_incr[9] = 0.7071067811865475*nuVtSqSum[2]*edgeSurf[21]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[15]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[9]; 
  edgeSurf_incr[10] = 0.6324555320336759*nuVtSqSum[1]*edgeSurf[17]+0.6324555320336759*nuVtSqSum[2]*edgeSurf[10]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[10]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[6]; 
  edgeSurf_incr[11] = 0.45175395145262565*nuVtSqSum[2]*edgeSurf[11]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[11]+0.632455532033676*nuVtSqSum[1]*edgeSurf[4]+0.7071067811865475*edgeSurf[2]*nuVtSqSum[2]; 
  edgeSurf_incr[12] = 0.632455532033676*nuVtSqSum[1]*edgeSurf[20]+0.6324555320336759*nuVtSqSum[2]*edgeSurf[12]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[12]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[8]; 
  edgeSurf_incr[13] = 0.45175395145262565*nuVtSqSum[2]*edgeSurf[13]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[13]+0.632455532033676*nuVtSqSum[1]*edgeSurf[5]+0.7071067811865475*nuVtSqSum[2]*edgeSurf[3]; 
  edgeSurf_incr[14] = 0.7071067811865475*nuVtSqSum[2]*edgeSurf[23]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[18]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[14]; 
  edgeSurf_incr[15] = 0.632455532033676*nuVtSqSum[1]*edgeSurf[21]+0.6324555320336759*nuVtSqSum[2]*edgeSurf[15]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[15]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[9]; 
  edgeSurf_incr[16] = 0.7071067811865475*nuVtSqSum[2]*edgeSurf[24]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[19]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[16]; 
  edgeSurf_incr[17] = 0.45175395145262565*nuVtSqSum[2]*edgeSurf[17]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[17]+0.6324555320336759*nuVtSqSum[1]*edgeSurf[10]+0.7071067811865475*nuVtSqSum[2]*edgeSurf[6]; 
  edgeSurf_incr[18] = 0.6324555320336759*nuVtSqSum[1]*edgeSurf[23]+0.6324555320336759*nuVtSqSum[2]*edgeSurf[18]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[18]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[14]; 
  edgeSurf_incr[19] = 0.6324555320336759*nuVtSqSum[1]*edgeSurf[24]+0.6324555320336759*nuVtSqSum[2]*edgeSurf[19]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[19]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[16]; 
  edgeSurf_incr[20] = 0.45175395145262565*nuVtSqSum[2]*edgeSurf[20]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[20]+0.632455532033676*nuVtSqSum[1]*edgeSurf[12]+0.7071067811865475*nuVtSqSum[2]*edgeSurf[8]; 
  edgeSurf_incr[21] = 0.45175395145262565*nuVtSqSum[2]*edgeSurf[21]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[21]+0.632455532033676*nuVtSqSum[1]*edgeSurf[15]+0.7071067811865475*nuVtSqSum[2]*edgeSurf[9]; 
  edgeSurf_incr[22] = 0.7071067811865475*nuVtSqSum[2]*edgeSurf[26]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[25]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[22]; 
  edgeSurf_incr[23] = 0.45175395145262565*nuVtSqSum[2]*edgeSurf[23]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[23]+0.6324555320336759*nuVtSqSum[1]*edgeSurf[18]+0.7071067811865475*nuVtSqSum[2]*edgeSurf[14]; 
  edgeSurf_incr[24] = 0.45175395145262565*nuVtSqSum[2]*edgeSurf[24]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[24]+0.6324555320336759*nuVtSqSum[1]*edgeSurf[19]+0.7071067811865475*nuVtSqSum[2]*edgeSurf[16]; 
  edgeSurf_incr[25] = 0.6324555320336759*nuVtSqSum[1]*edgeSurf[26]+0.6324555320336759*nuVtSqSum[2]*edgeSurf[25]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[25]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[22]; 
  edgeSurf_incr[26] = 0.45175395145262565*nuVtSqSum[2]*edgeSurf[26]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[26]+0.6324555320336759*nuVtSqSum[1]*edgeSurf[25]+0.7071067811865475*nuVtSqSum[2]*edgeSurf[22]; 

  boundSurf_incr[0] = 0.7071067811865475*nuVtSqSum[2]*boundSurf[7]+0.7071067811865475*boundSurf[1]*nuVtSqSum[1]+0.7071067811865475*boundSurf[0]*nuVtSqSum[0]; 
  boundSurf_incr[1] = 0.6324555320336759*nuVtSqSum[1]*boundSurf[7]+0.6324555320336759*boundSurf[1]*nuVtSqSum[2]+0.7071067811865475*boundSurf[0]*nuVtSqSum[1]+0.7071067811865475*nuVtSqSum[0]*boundSurf[1]; 
  boundSurf_incr[2] = 0.7071067811865475*nuVtSqSum[2]*boundSurf[11]+0.7071067811865475*nuVtSqSum[1]*boundSurf[4]+0.7071067811865475*nuVtSqSum[0]*boundSurf[2]; 
  boundSurf_incr[3] = 0.7071067811865475*nuVtSqSum[2]*boundSurf[13]+0.7071067811865475*nuVtSqSum[1]*boundSurf[5]+0.7071067811865475*nuVtSqSum[0]*boundSurf[3]; 
  boundSurf_incr[4] = 0.632455532033676*nuVtSqSum[1]*boundSurf[11]+0.6324555320336759*nuVtSqSum[2]*boundSurf[4]+0.7071067811865475*nuVtSqSum[0]*boundSurf[4]+0.7071067811865475*nuVtSqSum[1]*boundSurf[2]; 
  boundSurf_incr[5] = 0.632455532033676*nuVtSqSum[1]*boundSurf[13]+0.6324555320336759*nuVtSqSum[2]*boundSurf[5]+0.7071067811865475*nuVtSqSum[0]*boundSurf[5]+0.7071067811865475*nuVtSqSum[1]*boundSurf[3]; 
  boundSurf_incr[6] = 0.7071067811865475*nuVtSqSum[2]*boundSurf[17]+0.7071067811865475*nuVtSqSum[1]*boundSurf[10]+0.7071067811865475*nuVtSqSum[0]*boundSurf[6]; 
  boundSurf_incr[7] = 0.45175395145262565*nuVtSqSum[2]*boundSurf[7]+0.7071067811865475*nuVtSqSum[0]*boundSurf[7]+0.7071067811865475*boundSurf[0]*nuVtSqSum[2]+0.6324555320336759*boundSurf[1]*nuVtSqSum[1]; 
  boundSurf_incr[8] = 0.7071067811865475*nuVtSqSum[2]*boundSurf[20]+0.7071067811865475*nuVtSqSum[1]*boundSurf[12]+0.7071067811865475*nuVtSqSum[0]*boundSurf[8]; 
  boundSurf_incr[9] = 0.7071067811865475*nuVtSqSum[2]*boundSurf[21]+0.7071067811865475*nuVtSqSum[1]*boundSurf[15]+0.7071067811865475*nuVtSqSum[0]*boundSurf[9]; 
  boundSurf_incr[10] = 0.6324555320336759*nuVtSqSum[1]*boundSurf[17]+0.6324555320336759*nuVtSqSum[2]*boundSurf[10]+0.7071067811865475*nuVtSqSum[0]*boundSurf[10]+0.7071067811865475*nuVtSqSum[1]*boundSurf[6]; 
  boundSurf_incr[11] = 0.45175395145262565*nuVtSqSum[2]*boundSurf[11]+0.7071067811865475*nuVtSqSum[0]*boundSurf[11]+0.632455532033676*nuVtSqSum[1]*boundSurf[4]+0.7071067811865475*boundSurf[2]*nuVtSqSum[2]; 
  boundSurf_incr[12] = 0.632455532033676*nuVtSqSum[1]*boundSurf[20]+0.6324555320336759*nuVtSqSum[2]*boundSurf[12]+0.7071067811865475*nuVtSqSum[0]*boundSurf[12]+0.7071067811865475*nuVtSqSum[1]*boundSurf[8]; 
  boundSurf_incr[13] = 0.45175395145262565*nuVtSqSum[2]*boundSurf[13]+0.7071067811865475*nuVtSqSum[0]*boundSurf[13]+0.632455532033676*nuVtSqSum[1]*boundSurf[5]+0.7071067811865475*nuVtSqSum[2]*boundSurf[3]; 
  boundSurf_incr[14] = 0.7071067811865475*nuVtSqSum[2]*boundSurf[23]+0.7071067811865475*nuVtSqSum[1]*boundSurf[18]+0.7071067811865475*nuVtSqSum[0]*boundSurf[14]; 
  boundSurf_incr[15] = 0.632455532033676*nuVtSqSum[1]*boundSurf[21]+0.6324555320336759*nuVtSqSum[2]*boundSurf[15]+0.7071067811865475*nuVtSqSum[0]*boundSurf[15]+0.7071067811865475*nuVtSqSum[1]*boundSurf[9]; 
  boundSurf_incr[16] = 0.7071067811865475*nuVtSqSum[2]*boundSurf[24]+0.7071067811865475*nuVtSqSum[1]*boundSurf[19]+0.7071067811865475*nuVtSqSum[0]*boundSurf[16]; 
  boundSurf_incr[17] = 0.45175395145262565*nuVtSqSum[2]*boundSurf[17]+0.7071067811865475*nuVtSqSum[0]*boundSurf[17]+0.6324555320336759*nuVtSqSum[1]*boundSurf[10]+0.7071067811865475*nuVtSqSum[2]*boundSurf[6]; 
  boundSurf_incr[18] = 0.6324555320336759*nuVtSqSum[1]*boundSurf[23]+0.6324555320336759*nuVtSqSum[2]*boundSurf[18]+0.7071067811865475*nuVtSqSum[0]*boundSurf[18]+0.7071067811865475*nuVtSqSum[1]*boundSurf[14]; 
  boundSurf_incr[19] = 0.6324555320336759*nuVtSqSum[1]*boundSurf[24]+0.6324555320336759*nuVtSqSum[2]*boundSurf[19]+0.7071067811865475*nuVtSqSum[0]*boundSurf[19]+0.7071067811865475*nuVtSqSum[1]*boundSurf[16]; 
  boundSurf_incr[20] = 0.45175395145262565*nuVtSqSum[2]*boundSurf[20]+0.7071067811865475*nuVtSqSum[0]*boundSurf[20]+0.632455532033676*nuVtSqSum[1]*boundSurf[12]+0.7071067811865475*nuVtSqSum[2]*boundSurf[8]; 
  boundSurf_incr[21] = 0.45175395145262565*nuVtSqSum[2]*boundSurf[21]+0.7071067811865475*nuVtSqSum[0]*boundSurf[21]+0.632455532033676*nuVtSqSum[1]*boundSurf[15]+0.7071067811865475*nuVtSqSum[2]*boundSurf[9]; 
  boundSurf_incr[22] = 0.7071067811865475*nuVtSqSum[2]*boundSurf[26]+0.7071067811865475*nuVtSqSum[1]*boundSurf[25]+0.7071067811865475*nuVtSqSum[0]*boundSurf[22]; 
  boundSurf_incr[23] = 0.45175395145262565*nuVtSqSum[2]*boundSurf[23]+0.7071067811865475*nuVtSqSum[0]*boundSurf[23]+0.6324555320336759*nuVtSqSum[1]*boundSurf[18]+0.7071067811865475*nuVtSqSum[2]*boundSurf[14]; 
  boundSurf_incr[24] = 0.45175395145262565*nuVtSqSum[2]*boundSurf[24]+0.7071067811865475*nuVtSqSum[0]*boundSurf[24]+0.6324555320336759*nuVtSqSum[1]*boundSurf[19]+0.7071067811865475*nuVtSqSum[2]*boundSurf[16]; 
  boundSurf_incr[25] = 0.6324555320336759*nuVtSqSum[1]*boundSurf[26]+0.6324555320336759*nuVtSqSum[2]*boundSurf[25]+0.7071067811865475*nuVtSqSum[0]*boundSurf[25]+0.7071067811865475*nuVtSqSum[1]*boundSurf[22]; 
  boundSurf_incr[26] = 0.45175395145262565*nuVtSqSum[2]*boundSurf[26]+0.7071067811865475*nuVtSqSum[0]*boundSurf[26]+0.6324555320336759*nuVtSqSum[1]*boundSurf[25]+0.7071067811865475*nuVtSqSum[2]*boundSurf[22]; 

  } 

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*rdvSq4; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*rdvSq4; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*rdvSq4; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*rdvSq4; 
  out[4] += (vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*rdvSq4; 
  out[5] += (vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*rdvSq4; 
  out[6] += (vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*rdvSq4; 
  out[7] += (vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*rdvSq4; 
  out[8] += (vol_incr[8]+edgeSurf_incr[8]+boundSurf_incr[8])*rdvSq4; 
  out[9] += (vol_incr[9]+edgeSurf_incr[9]+boundSurf_incr[9])*rdvSq4; 
  out[10] += (vol_incr[10]+edgeSurf_incr[10]+boundSurf_incr[10])*rdvSq4; 
  out[11] += (vol_incr[11]+edgeSurf_incr[11]+boundSurf_incr[11])*rdvSq4; 
  out[12] += (vol_incr[12]+edgeSurf_incr[12]+boundSurf_incr[12])*rdvSq4; 
  out[13] += (vol_incr[13]+edgeSurf_incr[13]+boundSurf_incr[13])*rdvSq4; 
  out[14] += (vol_incr[14]+edgeSurf_incr[14]+boundSurf_incr[14])*rdvSq4; 
  out[15] += (vol_incr[15]+edgeSurf_incr[15]+boundSurf_incr[15])*rdvSq4; 
  out[16] += (vol_incr[16]+edgeSurf_incr[16]+boundSurf_incr[16])*rdvSq4; 
  out[17] += (vol_incr[17]+edgeSurf_incr[17]+boundSurf_incr[17])*rdvSq4; 
  out[18] += (vol_incr[18]+edgeSurf_incr[18]+boundSurf_incr[18])*rdvSq4; 
  out[19] += (vol_incr[19]+edgeSurf_incr[19]+boundSurf_incr[19])*rdvSq4; 
  out[20] += (vol_incr[20]+edgeSurf_incr[20]+boundSurf_incr[20])*rdvSq4; 
  out[21] += (vol_incr[21]+edgeSurf_incr[21]+boundSurf_incr[21])*rdvSq4; 
  out[22] += (vol_incr[22]+edgeSurf_incr[22]+boundSurf_incr[22])*rdvSq4; 
  out[23] += (vol_incr[23]+edgeSurf_incr[23]+boundSurf_incr[23])*rdvSq4; 
  out[24] += (vol_incr[24]+edgeSurf_incr[24]+boundSurf_incr[24])*rdvSq4; 
  out[25] += (vol_incr[25]+edgeSurf_incr[25]+boundSurf_incr[25])*rdvSq4; 
  out[26] += (vol_incr[26]+edgeSurf_incr[26]+boundSurf_incr[26])*rdvSq4; 

  return 0.;

} 
