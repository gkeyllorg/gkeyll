#include <gkyl_lbo_pkpm_kernels.h> 
GKYL_CU_DH double lbo_pkpm_diff_boundary_surfvpar_2x1v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:       Cell-center coordinates. 
  // dxv[NDIM]:     Cell spacing. 
  // nuSum:         Collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: Sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // fSkin/fEdge:   Input distribution functions [F_0, T_perp G = T_perp (F_1 - F_0)] in skin cell/last edge cell. 
  // out:           Incremented distribution functions in skin cell. 
  double rdvSq4 = 4.0/(dxv[2]*dxv[2]); 
  const double *nuVtSqSum = &nuPrimMomsSum[4];

  const double *F_0Skin = &fSkin[0]; 
  const double *G_1Skin = &fSkin[12]; 
  const double *F_0Edge = &fEdge[0]; 
  const double *G_1Edge = &fEdge[12]; 
  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[12]; 

  double vol_incr_F_0[12] = {0.0}; 
  double vol_incr_G_1[12] = {0.0}; 
  vol_incr_F_0[7] = 3.3541019662496847*nuVtSqSum[3]*F_0Skin[4]*rdvSq4+3.3541019662496847*F_0Skin[2]*nuVtSqSum[2]*rdvSq4+3.3541019662496847*F_0Skin[1]*nuVtSqSum[1]*rdvSq4+3.3541019662496847*F_0Skin[0]*nuVtSqSum[0]*rdvSq4; 
  vol_incr_F_0[9] = 3.3541019662496843*nuVtSqSum[2]*F_0Skin[4]*rdvSq4+3.3541019662496843*F_0Skin[2]*nuVtSqSum[3]*rdvSq4+3.3541019662496843*F_0Skin[0]*nuVtSqSum[1]*rdvSq4+3.3541019662496843*nuVtSqSum[0]*F_0Skin[1]*rdvSq4; 
  vol_incr_F_0[10] = 3.3541019662496843*nuVtSqSum[1]*F_0Skin[4]*rdvSq4+3.3541019662496843*F_0Skin[1]*nuVtSqSum[3]*rdvSq4+3.3541019662496843*F_0Skin[0]*nuVtSqSum[2]*rdvSq4+3.3541019662496843*nuVtSqSum[0]*F_0Skin[2]*rdvSq4; 
  vol_incr_F_0[11] = 3.3541019662496847*nuVtSqSum[0]*F_0Skin[4]*rdvSq4+3.3541019662496847*F_0Skin[0]*nuVtSqSum[3]*rdvSq4+3.3541019662496847*F_0Skin[1]*nuVtSqSum[2]*rdvSq4+3.3541019662496847*nuVtSqSum[1]*F_0Skin[2]*rdvSq4; 
  vol_incr_G_1[7] = 3.3541019662496847*nuVtSqSum[3]*G_1Skin[4]*rdvSq4+3.3541019662496847*G_1Skin[2]*nuVtSqSum[2]*rdvSq4+3.3541019662496847*G_1Skin[1]*nuVtSqSum[1]*rdvSq4+3.3541019662496847*G_1Skin[0]*nuVtSqSum[0]*rdvSq4; 
  vol_incr_G_1[9] = 3.3541019662496843*nuVtSqSum[2]*G_1Skin[4]*rdvSq4+3.3541019662496843*G_1Skin[2]*nuVtSqSum[3]*rdvSq4+3.3541019662496843*G_1Skin[0]*nuVtSqSum[1]*rdvSq4+3.3541019662496843*nuVtSqSum[0]*G_1Skin[1]*rdvSq4; 
  vol_incr_G_1[10] = 3.3541019662496843*nuVtSqSum[1]*G_1Skin[4]*rdvSq4+3.3541019662496843*G_1Skin[1]*nuVtSqSum[3]*rdvSq4+3.3541019662496843*G_1Skin[0]*nuVtSqSum[2]*rdvSq4+3.3541019662496843*nuVtSqSum[0]*G_1Skin[2]*rdvSq4; 
  vol_incr_G_1[11] = 3.3541019662496847*nuVtSqSum[0]*G_1Skin[4]*rdvSq4+3.3541019662496847*G_1Skin[0]*nuVtSqSum[3]*rdvSq4+3.3541019662496847*G_1Skin[1]*nuVtSqSum[2]*rdvSq4+3.3541019662496847*nuVtSqSum[1]*G_1Skin[2]*rdvSq4; 

  double temp_F_0_diff[12] = {0.0}; 
  double temp_F_0_edge[12] = {0.0}; 
  double diff_F_0_incr[12] = {0.0}; 
  double edge_F_0_incr[12] = {0.0}; 

  double temp_G_1_diff[12] = {0.0}; 
  double temp_G_1_edge[12] = {0.0}; 
  double diff_G_1_incr[12] = {0.0}; 
  double edge_G_1_incr[12] = {0.0}; 

  if (edge == -1) { 

  temp_F_0_diff[0] = -(0.6708203932499369*F_0Skin[7])+0.6708203932499369*F_0Edge[7]-1.190784930203603*F_0Skin[3]-1.190784930203603*F_0Edge[3]-0.9375*F_0Skin[0]+0.9375*F_0Edge[0]; 
  temp_F_0_diff[1] = -(0.6708203932499369*F_0Skin[9])+0.6708203932499369*F_0Edge[9]-1.190784930203603*F_0Skin[5]-1.190784930203603*F_0Edge[5]-0.9375*F_0Skin[1]+0.9375*F_0Edge[1]; 
  temp_F_0_diff[2] = -(0.6708203932499369*F_0Skin[10])+0.6708203932499369*F_0Edge[10]-1.190784930203603*F_0Skin[6]-1.190784930203603*F_0Edge[6]-0.9375*F_0Skin[2]+0.9375*F_0Edge[2]; 
  temp_F_0_diff[3] = -(1.5855025573536612*F_0Skin[7])+0.7382874503707888*F_0Edge[7]-2.671875*F_0Skin[3]-1.453125*F_0Edge[3]-2.0568103339880417*F_0Skin[0]+1.190784930203603*F_0Edge[0]; 
  temp_F_0_diff[4] = -(0.6708203932499369*F_0Skin[11])+0.6708203932499369*F_0Edge[11]-1.190784930203603*F_0Skin[8]-1.190784930203603*F_0Edge[8]-0.9375*F_0Skin[4]+0.9375*F_0Edge[4]; 
  temp_F_0_diff[5] = -(1.5855025573536612*F_0Skin[9])+0.7382874503707888*F_0Edge[9]-2.671875*F_0Skin[5]-1.453125*F_0Edge[5]-2.0568103339880417*F_0Skin[1]+1.190784930203603*F_0Edge[1]; 
  temp_F_0_diff[6] = -(1.5855025573536612*F_0Skin[10])+0.7382874503707888*F_0Edge[10]-2.671875*F_0Skin[6]-1.453125*F_0Edge[6]-2.0568103339880417*F_0Skin[2]+1.190784930203603*F_0Edge[2]; 
  temp_F_0_diff[7] = -(3.140625*F_0Skin[7])-0.140625*F_0Edge[7]-5.022775277112744*F_0Skin[3]-0.3025768239224545*F_0Edge[3]-3.7733647120308955*F_0Skin[0]+0.4192627457812106*F_0Edge[0]; 
  temp_F_0_diff[8] = -(1.5855025573536612*F_0Skin[11])+0.7382874503707888*F_0Edge[11]-2.671875*F_0Skin[8]-1.453125*F_0Edge[8]-2.0568103339880417*F_0Skin[4]+1.190784930203603*F_0Edge[4]; 
  temp_F_0_diff[9] = -(3.140625*F_0Skin[9])-0.140625*F_0Edge[9]-5.022775277112744*F_0Skin[5]-0.30257682392245444*F_0Edge[5]-3.773364712030894*F_0Skin[1]+0.41926274578121053*F_0Edge[1]; 
  temp_F_0_diff[10] = -(3.140625*F_0Skin[10])-0.140625*F_0Edge[10]-5.022775277112744*F_0Skin[6]-0.30257682392245444*F_0Edge[6]-3.773364712030894*F_0Skin[2]+0.41926274578121053*F_0Edge[2]; 
  temp_F_0_diff[11] = -(3.140625*F_0Skin[11])-0.140625*F_0Edge[11]-5.022775277112744*F_0Skin[8]-0.3025768239224545*F_0Edge[8]-3.7733647120308955*F_0Skin[4]+0.4192627457812106*F_0Edge[4]; 
  temp_G_1_diff[0] = -(0.6708203932499369*G_1Skin[7])+0.6708203932499369*G_1Edge[7]-1.190784930203603*G_1Skin[3]-1.190784930203603*G_1Edge[3]-0.9375*G_1Skin[0]+0.9375*G_1Edge[0]; 
  temp_G_1_diff[1] = -(0.6708203932499369*G_1Skin[9])+0.6708203932499369*G_1Edge[9]-1.190784930203603*G_1Skin[5]-1.190784930203603*G_1Edge[5]-0.9375*G_1Skin[1]+0.9375*G_1Edge[1]; 
  temp_G_1_diff[2] = -(0.6708203932499369*G_1Skin[10])+0.6708203932499369*G_1Edge[10]-1.190784930203603*G_1Skin[6]-1.190784930203603*G_1Edge[6]-0.9375*G_1Skin[2]+0.9375*G_1Edge[2]; 
  temp_G_1_diff[3] = -(1.5855025573536612*G_1Skin[7])+0.7382874503707888*G_1Edge[7]-2.671875*G_1Skin[3]-1.453125*G_1Edge[3]-2.0568103339880417*G_1Skin[0]+1.190784930203603*G_1Edge[0]; 
  temp_G_1_diff[4] = -(0.6708203932499369*G_1Skin[11])+0.6708203932499369*G_1Edge[11]-1.190784930203603*G_1Skin[8]-1.190784930203603*G_1Edge[8]-0.9375*G_1Skin[4]+0.9375*G_1Edge[4]; 
  temp_G_1_diff[5] = -(1.5855025573536612*G_1Skin[9])+0.7382874503707888*G_1Edge[9]-2.671875*G_1Skin[5]-1.453125*G_1Edge[5]-2.0568103339880417*G_1Skin[1]+1.190784930203603*G_1Edge[1]; 
  temp_G_1_diff[6] = -(1.5855025573536612*G_1Skin[10])+0.7382874503707888*G_1Edge[10]-2.671875*G_1Skin[6]-1.453125*G_1Edge[6]-2.0568103339880417*G_1Skin[2]+1.190784930203603*G_1Edge[2]; 
  temp_G_1_diff[7] = -(3.140625*G_1Skin[7])-0.140625*G_1Edge[7]-5.022775277112744*G_1Skin[3]-0.3025768239224545*G_1Edge[3]-3.7733647120308955*G_1Skin[0]+0.4192627457812106*G_1Edge[0]; 
  temp_G_1_diff[8] = -(1.5855025573536612*G_1Skin[11])+0.7382874503707888*G_1Edge[11]-2.671875*G_1Skin[8]-1.453125*G_1Edge[8]-2.0568103339880417*G_1Skin[4]+1.190784930203603*G_1Edge[4]; 
  temp_G_1_diff[9] = -(3.140625*G_1Skin[9])-0.140625*G_1Edge[9]-5.022775277112744*G_1Skin[5]-0.30257682392245444*G_1Edge[5]-3.773364712030894*G_1Skin[1]+0.41926274578121053*G_1Edge[1]; 
  temp_G_1_diff[10] = -(3.140625*G_1Skin[10])-0.140625*G_1Edge[10]-5.022775277112744*G_1Skin[6]-0.30257682392245444*G_1Edge[6]-3.773364712030894*G_1Skin[2]+0.41926274578121053*G_1Edge[2]; 
  temp_G_1_diff[11] = -(3.140625*G_1Skin[11])-0.140625*G_1Edge[11]-5.022775277112744*G_1Skin[8]-0.3025768239224545*G_1Edge[8]-3.7733647120308955*G_1Skin[4]+0.4192627457812106*G_1Edge[4]; 

  temp_F_0_edge[3] = 1.9364916731037085*F_0Skin[7]-1.5*F_0Skin[3]+0.8660254037844386*F_0Skin[0]; 
  temp_F_0_edge[5] = 1.9364916731037085*F_0Skin[9]-1.5*F_0Skin[5]+0.8660254037844386*F_0Skin[1]; 
  temp_F_0_edge[6] = 1.9364916731037085*F_0Skin[10]-1.5*F_0Skin[6]+0.8660254037844386*F_0Skin[2]; 
  temp_F_0_edge[7] = -(7.5*F_0Skin[7])+5.809475019311125*F_0Skin[3]-3.3541019662496847*F_0Skin[0]; 
  temp_F_0_edge[8] = 1.9364916731037085*F_0Skin[11]-1.5*F_0Skin[8]+0.8660254037844386*F_0Skin[4]; 
  temp_F_0_edge[9] = -(7.5*F_0Skin[9])+5.809475019311126*F_0Skin[5]-3.3541019662496843*F_0Skin[1]; 
  temp_F_0_edge[10] = -(7.5*F_0Skin[10])+5.809475019311126*F_0Skin[6]-3.3541019662496843*F_0Skin[2]; 
  temp_F_0_edge[11] = -(7.5*F_0Skin[11])+5.809475019311125*F_0Skin[8]-3.3541019662496847*F_0Skin[4]; 
  temp_G_1_edge[3] = 1.9364916731037085*G_1Skin[7]-1.5*G_1Skin[3]+0.8660254037844386*G_1Skin[0]; 
  temp_G_1_edge[5] = 1.9364916731037085*G_1Skin[9]-1.5*G_1Skin[5]+0.8660254037844386*G_1Skin[1]; 
  temp_G_1_edge[6] = 1.9364916731037085*G_1Skin[10]-1.5*G_1Skin[6]+0.8660254037844386*G_1Skin[2]; 
  temp_G_1_edge[7] = -(7.5*G_1Skin[7])+5.809475019311125*G_1Skin[3]-3.3541019662496847*G_1Skin[0]; 
  temp_G_1_edge[8] = 1.9364916731037085*G_1Skin[11]-1.5*G_1Skin[8]+0.8660254037844386*G_1Skin[4]; 
  temp_G_1_edge[9] = -(7.5*G_1Skin[9])+5.809475019311126*G_1Skin[5]-3.3541019662496843*G_1Skin[1]; 
  temp_G_1_edge[10] = -(7.5*G_1Skin[10])+5.809475019311126*G_1Skin[6]-3.3541019662496843*G_1Skin[2]; 
  temp_G_1_edge[11] = -(7.5*G_1Skin[11])+5.809475019311125*G_1Skin[8]-3.3541019662496847*G_1Skin[4]; 

  diff_F_0_incr[0] = 0.5*nuVtSqSum[3]*temp_F_0_diff[4]+0.5*nuVtSqSum[2]*temp_F_0_diff[2]+0.5*nuVtSqSum[1]*temp_F_0_diff[1]+0.5*nuVtSqSum[0]*temp_F_0_diff[0]; 
  diff_F_0_incr[1] = 0.5*nuVtSqSum[2]*temp_F_0_diff[4]+0.5*temp_F_0_diff[2]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_F_0_diff[1]+0.5*temp_F_0_diff[0]*nuVtSqSum[1]; 
  diff_F_0_incr[2] = 0.5*nuVtSqSum[1]*temp_F_0_diff[4]+0.5*temp_F_0_diff[1]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_F_0_diff[2]+0.5*temp_F_0_diff[0]*nuVtSqSum[2]; 
  diff_F_0_incr[3] = 0.5*nuVtSqSum[3]*temp_F_0_diff[8]+0.5*nuVtSqSum[2]*temp_F_0_diff[6]+0.5*nuVtSqSum[1]*temp_F_0_diff[5]+0.5*nuVtSqSum[0]*temp_F_0_diff[3]; 
  diff_F_0_incr[4] = 0.5*nuVtSqSum[0]*temp_F_0_diff[4]+0.5*temp_F_0_diff[0]*nuVtSqSum[3]+0.5*nuVtSqSum[1]*temp_F_0_diff[2]+0.5*temp_F_0_diff[1]*nuVtSqSum[2]; 
  diff_F_0_incr[5] = 0.5*nuVtSqSum[2]*temp_F_0_diff[8]+0.5*nuVtSqSum[3]*temp_F_0_diff[6]+0.5*nuVtSqSum[0]*temp_F_0_diff[5]+0.5*nuVtSqSum[1]*temp_F_0_diff[3]; 
  diff_F_0_incr[6] = 0.5*nuVtSqSum[1]*temp_F_0_diff[8]+0.5*nuVtSqSum[0]*temp_F_0_diff[6]+0.5*nuVtSqSum[3]*temp_F_0_diff[5]+0.5*nuVtSqSum[2]*temp_F_0_diff[3]; 
  diff_F_0_incr[7] = 0.5*nuVtSqSum[3]*temp_F_0_diff[11]+0.5000000000000001*nuVtSqSum[2]*temp_F_0_diff[10]+0.5000000000000001*nuVtSqSum[1]*temp_F_0_diff[9]+0.5*nuVtSqSum[0]*temp_F_0_diff[7]; 
  diff_F_0_incr[8] = 0.5*nuVtSqSum[0]*temp_F_0_diff[8]+0.5*nuVtSqSum[1]*temp_F_0_diff[6]+0.5*nuVtSqSum[2]*temp_F_0_diff[5]+0.5*nuVtSqSum[3]*temp_F_0_diff[3]; 
  diff_F_0_incr[9] = 0.5000000000000001*nuVtSqSum[2]*temp_F_0_diff[11]+0.5*nuVtSqSum[3]*temp_F_0_diff[10]+0.5*nuVtSqSum[0]*temp_F_0_diff[9]+0.5000000000000001*nuVtSqSum[1]*temp_F_0_diff[7]; 
  diff_F_0_incr[10] = 0.5000000000000001*nuVtSqSum[1]*temp_F_0_diff[11]+0.5*nuVtSqSum[0]*temp_F_0_diff[10]+0.5*nuVtSqSum[3]*temp_F_0_diff[9]+0.5000000000000001*nuVtSqSum[2]*temp_F_0_diff[7]; 
  diff_F_0_incr[11] = 0.5*nuVtSqSum[0]*temp_F_0_diff[11]+0.5000000000000001*nuVtSqSum[1]*temp_F_0_diff[10]+0.5000000000000001*nuVtSqSum[2]*temp_F_0_diff[9]+0.5*nuVtSqSum[3]*temp_F_0_diff[7]; 
  diff_G_1_incr[0] = 0.5*nuVtSqSum[3]*temp_G_1_diff[4]+0.5*nuVtSqSum[2]*temp_G_1_diff[2]+0.5*nuVtSqSum[1]*temp_G_1_diff[1]+0.5*nuVtSqSum[0]*temp_G_1_diff[0]; 
  diff_G_1_incr[1] = 0.5*nuVtSqSum[2]*temp_G_1_diff[4]+0.5*temp_G_1_diff[2]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_G_1_diff[1]+0.5*temp_G_1_diff[0]*nuVtSqSum[1]; 
  diff_G_1_incr[2] = 0.5*nuVtSqSum[1]*temp_G_1_diff[4]+0.5*temp_G_1_diff[1]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_G_1_diff[2]+0.5*temp_G_1_diff[0]*nuVtSqSum[2]; 
  diff_G_1_incr[3] = 0.5*nuVtSqSum[3]*temp_G_1_diff[8]+0.5*nuVtSqSum[2]*temp_G_1_diff[6]+0.5*nuVtSqSum[1]*temp_G_1_diff[5]+0.5*nuVtSqSum[0]*temp_G_1_diff[3]; 
  diff_G_1_incr[4] = 0.5*nuVtSqSum[0]*temp_G_1_diff[4]+0.5*temp_G_1_diff[0]*nuVtSqSum[3]+0.5*nuVtSqSum[1]*temp_G_1_diff[2]+0.5*temp_G_1_diff[1]*nuVtSqSum[2]; 
  diff_G_1_incr[5] = 0.5*nuVtSqSum[2]*temp_G_1_diff[8]+0.5*nuVtSqSum[3]*temp_G_1_diff[6]+0.5*nuVtSqSum[0]*temp_G_1_diff[5]+0.5*nuVtSqSum[1]*temp_G_1_diff[3]; 
  diff_G_1_incr[6] = 0.5*nuVtSqSum[1]*temp_G_1_diff[8]+0.5*nuVtSqSum[0]*temp_G_1_diff[6]+0.5*nuVtSqSum[3]*temp_G_1_diff[5]+0.5*nuVtSqSum[2]*temp_G_1_diff[3]; 
  diff_G_1_incr[7] = 0.5*nuVtSqSum[3]*temp_G_1_diff[11]+0.5000000000000001*nuVtSqSum[2]*temp_G_1_diff[10]+0.5000000000000001*nuVtSqSum[1]*temp_G_1_diff[9]+0.5*nuVtSqSum[0]*temp_G_1_diff[7]; 
  diff_G_1_incr[8] = 0.5*nuVtSqSum[0]*temp_G_1_diff[8]+0.5*nuVtSqSum[1]*temp_G_1_diff[6]+0.5*nuVtSqSum[2]*temp_G_1_diff[5]+0.5*nuVtSqSum[3]*temp_G_1_diff[3]; 
  diff_G_1_incr[9] = 0.5000000000000001*nuVtSqSum[2]*temp_G_1_diff[11]+0.5*nuVtSqSum[3]*temp_G_1_diff[10]+0.5*nuVtSqSum[0]*temp_G_1_diff[9]+0.5000000000000001*nuVtSqSum[1]*temp_G_1_diff[7]; 
  diff_G_1_incr[10] = 0.5000000000000001*nuVtSqSum[1]*temp_G_1_diff[11]+0.5*nuVtSqSum[0]*temp_G_1_diff[10]+0.5*nuVtSqSum[3]*temp_G_1_diff[9]+0.5000000000000001*nuVtSqSum[2]*temp_G_1_diff[7]; 
  diff_G_1_incr[11] = 0.5*nuVtSqSum[0]*temp_G_1_diff[11]+0.5000000000000001*nuVtSqSum[1]*temp_G_1_diff[10]+0.5000000000000001*nuVtSqSum[2]*temp_G_1_diff[9]+0.5*nuVtSqSum[3]*temp_G_1_diff[7]; 

  edge_F_0_incr[0] = 0.5*nuVtSqSum[3]*temp_F_0_edge[4]+0.5*nuVtSqSum[2]*temp_F_0_edge[2]+0.5*nuVtSqSum[1]*temp_F_0_edge[1]+0.5*nuVtSqSum[0]*temp_F_0_edge[0]; 
  edge_F_0_incr[1] = 0.5*nuVtSqSum[2]*temp_F_0_edge[4]+0.5*temp_F_0_edge[2]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_F_0_edge[1]+0.5*temp_F_0_edge[0]*nuVtSqSum[1]; 
  edge_F_0_incr[2] = 0.5*nuVtSqSum[1]*temp_F_0_edge[4]+0.5*temp_F_0_edge[1]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_F_0_edge[2]+0.5*temp_F_0_edge[0]*nuVtSqSum[2]; 
  edge_F_0_incr[3] = 0.5*nuVtSqSum[3]*temp_F_0_edge[8]+0.5*nuVtSqSum[2]*temp_F_0_edge[6]+0.5*nuVtSqSum[1]*temp_F_0_edge[5]+0.5*nuVtSqSum[0]*temp_F_0_edge[3]; 
  edge_F_0_incr[4] = 0.5*nuVtSqSum[0]*temp_F_0_edge[4]+0.5*temp_F_0_edge[0]*nuVtSqSum[3]+0.5*nuVtSqSum[1]*temp_F_0_edge[2]+0.5*temp_F_0_edge[1]*nuVtSqSum[2]; 
  edge_F_0_incr[5] = 0.5*nuVtSqSum[2]*temp_F_0_edge[8]+0.5*nuVtSqSum[3]*temp_F_0_edge[6]+0.5*nuVtSqSum[0]*temp_F_0_edge[5]+0.5*nuVtSqSum[1]*temp_F_0_edge[3]; 
  edge_F_0_incr[6] = 0.5*nuVtSqSum[1]*temp_F_0_edge[8]+0.5*nuVtSqSum[0]*temp_F_0_edge[6]+0.5*nuVtSqSum[3]*temp_F_0_edge[5]+0.5*nuVtSqSum[2]*temp_F_0_edge[3]; 
  edge_F_0_incr[7] = 0.5*nuVtSqSum[3]*temp_F_0_edge[11]+0.5000000000000001*nuVtSqSum[2]*temp_F_0_edge[10]+0.5000000000000001*nuVtSqSum[1]*temp_F_0_edge[9]+0.5*nuVtSqSum[0]*temp_F_0_edge[7]; 
  edge_F_0_incr[8] = 0.5*nuVtSqSum[0]*temp_F_0_edge[8]+0.5*nuVtSqSum[1]*temp_F_0_edge[6]+0.5*nuVtSqSum[2]*temp_F_0_edge[5]+0.5*nuVtSqSum[3]*temp_F_0_edge[3]; 
  edge_F_0_incr[9] = 0.5000000000000001*nuVtSqSum[2]*temp_F_0_edge[11]+0.5*nuVtSqSum[3]*temp_F_0_edge[10]+0.5*nuVtSqSum[0]*temp_F_0_edge[9]+0.5000000000000001*nuVtSqSum[1]*temp_F_0_edge[7]; 
  edge_F_0_incr[10] = 0.5000000000000001*nuVtSqSum[1]*temp_F_0_edge[11]+0.5*nuVtSqSum[0]*temp_F_0_edge[10]+0.5*nuVtSqSum[3]*temp_F_0_edge[9]+0.5000000000000001*nuVtSqSum[2]*temp_F_0_edge[7]; 
  edge_F_0_incr[11] = 0.5*nuVtSqSum[0]*temp_F_0_edge[11]+0.5000000000000001*nuVtSqSum[1]*temp_F_0_edge[10]+0.5000000000000001*nuVtSqSum[2]*temp_F_0_edge[9]+0.5*nuVtSqSum[3]*temp_F_0_edge[7]; 
  edge_G_1_incr[0] = 0.5*nuVtSqSum[3]*temp_G_1_edge[4]+0.5*nuVtSqSum[2]*temp_G_1_edge[2]+0.5*nuVtSqSum[1]*temp_G_1_edge[1]+0.5*nuVtSqSum[0]*temp_G_1_edge[0]; 
  edge_G_1_incr[1] = 0.5*nuVtSqSum[2]*temp_G_1_edge[4]+0.5*temp_G_1_edge[2]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_G_1_edge[1]+0.5*temp_G_1_edge[0]*nuVtSqSum[1]; 
  edge_G_1_incr[2] = 0.5*nuVtSqSum[1]*temp_G_1_edge[4]+0.5*temp_G_1_edge[1]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_G_1_edge[2]+0.5*temp_G_1_edge[0]*nuVtSqSum[2]; 
  edge_G_1_incr[3] = 0.5*nuVtSqSum[3]*temp_G_1_edge[8]+0.5*nuVtSqSum[2]*temp_G_1_edge[6]+0.5*nuVtSqSum[1]*temp_G_1_edge[5]+0.5*nuVtSqSum[0]*temp_G_1_edge[3]; 
  edge_G_1_incr[4] = 0.5*nuVtSqSum[0]*temp_G_1_edge[4]+0.5*temp_G_1_edge[0]*nuVtSqSum[3]+0.5*nuVtSqSum[1]*temp_G_1_edge[2]+0.5*temp_G_1_edge[1]*nuVtSqSum[2]; 
  edge_G_1_incr[5] = 0.5*nuVtSqSum[2]*temp_G_1_edge[8]+0.5*nuVtSqSum[3]*temp_G_1_edge[6]+0.5*nuVtSqSum[0]*temp_G_1_edge[5]+0.5*nuVtSqSum[1]*temp_G_1_edge[3]; 
  edge_G_1_incr[6] = 0.5*nuVtSqSum[1]*temp_G_1_edge[8]+0.5*nuVtSqSum[0]*temp_G_1_edge[6]+0.5*nuVtSqSum[3]*temp_G_1_edge[5]+0.5*nuVtSqSum[2]*temp_G_1_edge[3]; 
  edge_G_1_incr[7] = 0.5*nuVtSqSum[3]*temp_G_1_edge[11]+0.5000000000000001*nuVtSqSum[2]*temp_G_1_edge[10]+0.5000000000000001*nuVtSqSum[1]*temp_G_1_edge[9]+0.5*nuVtSqSum[0]*temp_G_1_edge[7]; 
  edge_G_1_incr[8] = 0.5*nuVtSqSum[0]*temp_G_1_edge[8]+0.5*nuVtSqSum[1]*temp_G_1_edge[6]+0.5*nuVtSqSum[2]*temp_G_1_edge[5]+0.5*nuVtSqSum[3]*temp_G_1_edge[3]; 
  edge_G_1_incr[9] = 0.5000000000000001*nuVtSqSum[2]*temp_G_1_edge[11]+0.5*nuVtSqSum[3]*temp_G_1_edge[10]+0.5*nuVtSqSum[0]*temp_G_1_edge[9]+0.5000000000000001*nuVtSqSum[1]*temp_G_1_edge[7]; 
  edge_G_1_incr[10] = 0.5000000000000001*nuVtSqSum[1]*temp_G_1_edge[11]+0.5*nuVtSqSum[0]*temp_G_1_edge[10]+0.5*nuVtSqSum[3]*temp_G_1_edge[9]+0.5000000000000001*nuVtSqSum[2]*temp_G_1_edge[7]; 
  edge_G_1_incr[11] = 0.5*nuVtSqSum[0]*temp_G_1_edge[11]+0.5000000000000001*nuVtSqSum[1]*temp_G_1_edge[10]+0.5000000000000001*nuVtSqSum[2]*temp_G_1_edge[9]+0.5*nuVtSqSum[3]*temp_G_1_edge[7]; 


  } else { 

  temp_F_0_diff[0] = -(0.6708203932499369*F_0Skin[7])+0.6708203932499369*F_0Edge[7]+1.190784930203603*F_0Skin[3]+1.190784930203603*F_0Edge[3]-0.9375*F_0Skin[0]+0.9375*F_0Edge[0]; 
  temp_F_0_diff[1] = -(0.6708203932499369*F_0Skin[9])+0.6708203932499369*F_0Edge[9]+1.190784930203603*F_0Skin[5]+1.190784930203603*F_0Edge[5]-0.9375*F_0Skin[1]+0.9375*F_0Edge[1]; 
  temp_F_0_diff[2] = -(0.6708203932499369*F_0Skin[10])+0.6708203932499369*F_0Edge[10]+1.190784930203603*F_0Skin[6]+1.190784930203603*F_0Edge[6]-0.9375*F_0Skin[2]+0.9375*F_0Edge[2]; 
  temp_F_0_diff[3] = 1.5855025573536612*F_0Skin[7]-0.7382874503707888*F_0Edge[7]-2.671875*F_0Skin[3]-1.453125*F_0Edge[3]+2.0568103339880417*F_0Skin[0]-1.190784930203603*F_0Edge[0]; 
  temp_F_0_diff[4] = -(0.6708203932499369*F_0Skin[11])+0.6708203932499369*F_0Edge[11]+1.190784930203603*F_0Skin[8]+1.190784930203603*F_0Edge[8]-0.9375*F_0Skin[4]+0.9375*F_0Edge[4]; 
  temp_F_0_diff[5] = 1.5855025573536612*F_0Skin[9]-0.7382874503707888*F_0Edge[9]-2.671875*F_0Skin[5]-1.453125*F_0Edge[5]+2.0568103339880417*F_0Skin[1]-1.190784930203603*F_0Edge[1]; 
  temp_F_0_diff[6] = 1.5855025573536612*F_0Skin[10]-0.7382874503707888*F_0Edge[10]-2.671875*F_0Skin[6]-1.453125*F_0Edge[6]+2.0568103339880417*F_0Skin[2]-1.190784930203603*F_0Edge[2]; 
  temp_F_0_diff[7] = -(3.140625*F_0Skin[7])-0.140625*F_0Edge[7]+5.022775277112744*F_0Skin[3]+0.3025768239224545*F_0Edge[3]-3.7733647120308955*F_0Skin[0]+0.4192627457812106*F_0Edge[0]; 
  temp_F_0_diff[8] = 1.5855025573536612*F_0Skin[11]-0.7382874503707888*F_0Edge[11]-2.671875*F_0Skin[8]-1.453125*F_0Edge[8]+2.0568103339880417*F_0Skin[4]-1.190784930203603*F_0Edge[4]; 
  temp_F_0_diff[9] = -(3.140625*F_0Skin[9])-0.140625*F_0Edge[9]+5.022775277112744*F_0Skin[5]+0.30257682392245444*F_0Edge[5]-3.773364712030894*F_0Skin[1]+0.41926274578121053*F_0Edge[1]; 
  temp_F_0_diff[10] = -(3.140625*F_0Skin[10])-0.140625*F_0Edge[10]+5.022775277112744*F_0Skin[6]+0.30257682392245444*F_0Edge[6]-3.773364712030894*F_0Skin[2]+0.41926274578121053*F_0Edge[2]; 
  temp_F_0_diff[11] = -(3.140625*F_0Skin[11])-0.140625*F_0Edge[11]+5.022775277112744*F_0Skin[8]+0.3025768239224545*F_0Edge[8]-3.7733647120308955*F_0Skin[4]+0.4192627457812106*F_0Edge[4]; 
  temp_G_1_diff[0] = -(0.6708203932499369*G_1Skin[7])+0.6708203932499369*G_1Edge[7]+1.190784930203603*G_1Skin[3]+1.190784930203603*G_1Edge[3]-0.9375*G_1Skin[0]+0.9375*G_1Edge[0]; 
  temp_G_1_diff[1] = -(0.6708203932499369*G_1Skin[9])+0.6708203932499369*G_1Edge[9]+1.190784930203603*G_1Skin[5]+1.190784930203603*G_1Edge[5]-0.9375*G_1Skin[1]+0.9375*G_1Edge[1]; 
  temp_G_1_diff[2] = -(0.6708203932499369*G_1Skin[10])+0.6708203932499369*G_1Edge[10]+1.190784930203603*G_1Skin[6]+1.190784930203603*G_1Edge[6]-0.9375*G_1Skin[2]+0.9375*G_1Edge[2]; 
  temp_G_1_diff[3] = 1.5855025573536612*G_1Skin[7]-0.7382874503707888*G_1Edge[7]-2.671875*G_1Skin[3]-1.453125*G_1Edge[3]+2.0568103339880417*G_1Skin[0]-1.190784930203603*G_1Edge[0]; 
  temp_G_1_diff[4] = -(0.6708203932499369*G_1Skin[11])+0.6708203932499369*G_1Edge[11]+1.190784930203603*G_1Skin[8]+1.190784930203603*G_1Edge[8]-0.9375*G_1Skin[4]+0.9375*G_1Edge[4]; 
  temp_G_1_diff[5] = 1.5855025573536612*G_1Skin[9]-0.7382874503707888*G_1Edge[9]-2.671875*G_1Skin[5]-1.453125*G_1Edge[5]+2.0568103339880417*G_1Skin[1]-1.190784930203603*G_1Edge[1]; 
  temp_G_1_diff[6] = 1.5855025573536612*G_1Skin[10]-0.7382874503707888*G_1Edge[10]-2.671875*G_1Skin[6]-1.453125*G_1Edge[6]+2.0568103339880417*G_1Skin[2]-1.190784930203603*G_1Edge[2]; 
  temp_G_1_diff[7] = -(3.140625*G_1Skin[7])-0.140625*G_1Edge[7]+5.022775277112744*G_1Skin[3]+0.3025768239224545*G_1Edge[3]-3.7733647120308955*G_1Skin[0]+0.4192627457812106*G_1Edge[0]; 
  temp_G_1_diff[8] = 1.5855025573536612*G_1Skin[11]-0.7382874503707888*G_1Edge[11]-2.671875*G_1Skin[8]-1.453125*G_1Edge[8]+2.0568103339880417*G_1Skin[4]-1.190784930203603*G_1Edge[4]; 
  temp_G_1_diff[9] = -(3.140625*G_1Skin[9])-0.140625*G_1Edge[9]+5.022775277112744*G_1Skin[5]+0.30257682392245444*G_1Edge[5]-3.773364712030894*G_1Skin[1]+0.41926274578121053*G_1Edge[1]; 
  temp_G_1_diff[10] = -(3.140625*G_1Skin[10])-0.140625*G_1Edge[10]+5.022775277112744*G_1Skin[6]+0.30257682392245444*G_1Edge[6]-3.773364712030894*G_1Skin[2]+0.41926274578121053*G_1Edge[2]; 
  temp_G_1_diff[11] = -(3.140625*G_1Skin[11])-0.140625*G_1Edge[11]+5.022775277112744*G_1Skin[8]+0.3025768239224545*G_1Edge[8]-3.7733647120308955*G_1Skin[4]+0.4192627457812106*G_1Edge[4]; 

  temp_F_0_edge[3] = -(1.9364916731037085*F_0Skin[7])-1.5*F_0Skin[3]-0.8660254037844386*F_0Skin[0]; 
  temp_F_0_edge[5] = -(1.9364916731037085*F_0Skin[9])-1.5*F_0Skin[5]-0.8660254037844386*F_0Skin[1]; 
  temp_F_0_edge[6] = -(1.9364916731037085*F_0Skin[10])-1.5*F_0Skin[6]-0.8660254037844386*F_0Skin[2]; 
  temp_F_0_edge[7] = -(7.5*F_0Skin[7])-5.809475019311125*F_0Skin[3]-3.3541019662496847*F_0Skin[0]; 
  temp_F_0_edge[8] = -(1.9364916731037085*F_0Skin[11])-1.5*F_0Skin[8]-0.8660254037844386*F_0Skin[4]; 
  temp_F_0_edge[9] = -(7.5*F_0Skin[9])-5.809475019311126*F_0Skin[5]-3.3541019662496843*F_0Skin[1]; 
  temp_F_0_edge[10] = -(7.5*F_0Skin[10])-5.809475019311126*F_0Skin[6]-3.3541019662496843*F_0Skin[2]; 
  temp_F_0_edge[11] = -(7.5*F_0Skin[11])-5.809475019311125*F_0Skin[8]-3.3541019662496847*F_0Skin[4]; 
  temp_G_1_edge[3] = -(1.9364916731037085*G_1Skin[7])-1.5*G_1Skin[3]-0.8660254037844386*G_1Skin[0]; 
  temp_G_1_edge[5] = -(1.9364916731037085*G_1Skin[9])-1.5*G_1Skin[5]-0.8660254037844386*G_1Skin[1]; 
  temp_G_1_edge[6] = -(1.9364916731037085*G_1Skin[10])-1.5*G_1Skin[6]-0.8660254037844386*G_1Skin[2]; 
  temp_G_1_edge[7] = -(7.5*G_1Skin[7])-5.809475019311125*G_1Skin[3]-3.3541019662496847*G_1Skin[0]; 
  temp_G_1_edge[8] = -(1.9364916731037085*G_1Skin[11])-1.5*G_1Skin[8]-0.8660254037844386*G_1Skin[4]; 
  temp_G_1_edge[9] = -(7.5*G_1Skin[9])-5.809475019311126*G_1Skin[5]-3.3541019662496843*G_1Skin[1]; 
  temp_G_1_edge[10] = -(7.5*G_1Skin[10])-5.809475019311126*G_1Skin[6]-3.3541019662496843*G_1Skin[2]; 
  temp_G_1_edge[11] = -(7.5*G_1Skin[11])-5.809475019311125*G_1Skin[8]-3.3541019662496847*G_1Skin[4]; 

  diff_F_0_incr[0] = 0.5*nuVtSqSum[3]*temp_F_0_diff[4]+0.5*nuVtSqSum[2]*temp_F_0_diff[2]+0.5*nuVtSqSum[1]*temp_F_0_diff[1]+0.5*nuVtSqSum[0]*temp_F_0_diff[0]; 
  diff_F_0_incr[1] = 0.5*nuVtSqSum[2]*temp_F_0_diff[4]+0.5*temp_F_0_diff[2]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_F_0_diff[1]+0.5*temp_F_0_diff[0]*nuVtSqSum[1]; 
  diff_F_0_incr[2] = 0.5*nuVtSqSum[1]*temp_F_0_diff[4]+0.5*temp_F_0_diff[1]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_F_0_diff[2]+0.5*temp_F_0_diff[0]*nuVtSqSum[2]; 
  diff_F_0_incr[3] = 0.5*nuVtSqSum[3]*temp_F_0_diff[8]+0.5*nuVtSqSum[2]*temp_F_0_diff[6]+0.5*nuVtSqSum[1]*temp_F_0_diff[5]+0.5*nuVtSqSum[0]*temp_F_0_diff[3]; 
  diff_F_0_incr[4] = 0.5*nuVtSqSum[0]*temp_F_0_diff[4]+0.5*temp_F_0_diff[0]*nuVtSqSum[3]+0.5*nuVtSqSum[1]*temp_F_0_diff[2]+0.5*temp_F_0_diff[1]*nuVtSqSum[2]; 
  diff_F_0_incr[5] = 0.5*nuVtSqSum[2]*temp_F_0_diff[8]+0.5*nuVtSqSum[3]*temp_F_0_diff[6]+0.5*nuVtSqSum[0]*temp_F_0_diff[5]+0.5*nuVtSqSum[1]*temp_F_0_diff[3]; 
  diff_F_0_incr[6] = 0.5*nuVtSqSum[1]*temp_F_0_diff[8]+0.5*nuVtSqSum[0]*temp_F_0_diff[6]+0.5*nuVtSqSum[3]*temp_F_0_diff[5]+0.5*nuVtSqSum[2]*temp_F_0_diff[3]; 
  diff_F_0_incr[7] = 0.5*nuVtSqSum[3]*temp_F_0_diff[11]+0.5000000000000001*nuVtSqSum[2]*temp_F_0_diff[10]+0.5000000000000001*nuVtSqSum[1]*temp_F_0_diff[9]+0.5*nuVtSqSum[0]*temp_F_0_diff[7]; 
  diff_F_0_incr[8] = 0.5*nuVtSqSum[0]*temp_F_0_diff[8]+0.5*nuVtSqSum[1]*temp_F_0_diff[6]+0.5*nuVtSqSum[2]*temp_F_0_diff[5]+0.5*nuVtSqSum[3]*temp_F_0_diff[3]; 
  diff_F_0_incr[9] = 0.5000000000000001*nuVtSqSum[2]*temp_F_0_diff[11]+0.5*nuVtSqSum[3]*temp_F_0_diff[10]+0.5*nuVtSqSum[0]*temp_F_0_diff[9]+0.5000000000000001*nuVtSqSum[1]*temp_F_0_diff[7]; 
  diff_F_0_incr[10] = 0.5000000000000001*nuVtSqSum[1]*temp_F_0_diff[11]+0.5*nuVtSqSum[0]*temp_F_0_diff[10]+0.5*nuVtSqSum[3]*temp_F_0_diff[9]+0.5000000000000001*nuVtSqSum[2]*temp_F_0_diff[7]; 
  diff_F_0_incr[11] = 0.5*nuVtSqSum[0]*temp_F_0_diff[11]+0.5000000000000001*nuVtSqSum[1]*temp_F_0_diff[10]+0.5000000000000001*nuVtSqSum[2]*temp_F_0_diff[9]+0.5*nuVtSqSum[3]*temp_F_0_diff[7]; 
  diff_G_1_incr[0] = 0.5*nuVtSqSum[3]*temp_G_1_diff[4]+0.5*nuVtSqSum[2]*temp_G_1_diff[2]+0.5*nuVtSqSum[1]*temp_G_1_diff[1]+0.5*nuVtSqSum[0]*temp_G_1_diff[0]; 
  diff_G_1_incr[1] = 0.5*nuVtSqSum[2]*temp_G_1_diff[4]+0.5*temp_G_1_diff[2]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_G_1_diff[1]+0.5*temp_G_1_diff[0]*nuVtSqSum[1]; 
  diff_G_1_incr[2] = 0.5*nuVtSqSum[1]*temp_G_1_diff[4]+0.5*temp_G_1_diff[1]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_G_1_diff[2]+0.5*temp_G_1_diff[0]*nuVtSqSum[2]; 
  diff_G_1_incr[3] = 0.5*nuVtSqSum[3]*temp_G_1_diff[8]+0.5*nuVtSqSum[2]*temp_G_1_diff[6]+0.5*nuVtSqSum[1]*temp_G_1_diff[5]+0.5*nuVtSqSum[0]*temp_G_1_diff[3]; 
  diff_G_1_incr[4] = 0.5*nuVtSqSum[0]*temp_G_1_diff[4]+0.5*temp_G_1_diff[0]*nuVtSqSum[3]+0.5*nuVtSqSum[1]*temp_G_1_diff[2]+0.5*temp_G_1_diff[1]*nuVtSqSum[2]; 
  diff_G_1_incr[5] = 0.5*nuVtSqSum[2]*temp_G_1_diff[8]+0.5*nuVtSqSum[3]*temp_G_1_diff[6]+0.5*nuVtSqSum[0]*temp_G_1_diff[5]+0.5*nuVtSqSum[1]*temp_G_1_diff[3]; 
  diff_G_1_incr[6] = 0.5*nuVtSqSum[1]*temp_G_1_diff[8]+0.5*nuVtSqSum[0]*temp_G_1_diff[6]+0.5*nuVtSqSum[3]*temp_G_1_diff[5]+0.5*nuVtSqSum[2]*temp_G_1_diff[3]; 
  diff_G_1_incr[7] = 0.5*nuVtSqSum[3]*temp_G_1_diff[11]+0.5000000000000001*nuVtSqSum[2]*temp_G_1_diff[10]+0.5000000000000001*nuVtSqSum[1]*temp_G_1_diff[9]+0.5*nuVtSqSum[0]*temp_G_1_diff[7]; 
  diff_G_1_incr[8] = 0.5*nuVtSqSum[0]*temp_G_1_diff[8]+0.5*nuVtSqSum[1]*temp_G_1_diff[6]+0.5*nuVtSqSum[2]*temp_G_1_diff[5]+0.5*nuVtSqSum[3]*temp_G_1_diff[3]; 
  diff_G_1_incr[9] = 0.5000000000000001*nuVtSqSum[2]*temp_G_1_diff[11]+0.5*nuVtSqSum[3]*temp_G_1_diff[10]+0.5*nuVtSqSum[0]*temp_G_1_diff[9]+0.5000000000000001*nuVtSqSum[1]*temp_G_1_diff[7]; 
  diff_G_1_incr[10] = 0.5000000000000001*nuVtSqSum[1]*temp_G_1_diff[11]+0.5*nuVtSqSum[0]*temp_G_1_diff[10]+0.5*nuVtSqSum[3]*temp_G_1_diff[9]+0.5000000000000001*nuVtSqSum[2]*temp_G_1_diff[7]; 
  diff_G_1_incr[11] = 0.5*nuVtSqSum[0]*temp_G_1_diff[11]+0.5000000000000001*nuVtSqSum[1]*temp_G_1_diff[10]+0.5000000000000001*nuVtSqSum[2]*temp_G_1_diff[9]+0.5*nuVtSqSum[3]*temp_G_1_diff[7]; 

  edge_F_0_incr[0] = 0.5*nuVtSqSum[3]*temp_F_0_edge[4]+0.5*nuVtSqSum[2]*temp_F_0_edge[2]+0.5*nuVtSqSum[1]*temp_F_0_edge[1]+0.5*nuVtSqSum[0]*temp_F_0_edge[0]; 
  edge_F_0_incr[1] = 0.5*nuVtSqSum[2]*temp_F_0_edge[4]+0.5*temp_F_0_edge[2]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_F_0_edge[1]+0.5*temp_F_0_edge[0]*nuVtSqSum[1]; 
  edge_F_0_incr[2] = 0.5*nuVtSqSum[1]*temp_F_0_edge[4]+0.5*temp_F_0_edge[1]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_F_0_edge[2]+0.5*temp_F_0_edge[0]*nuVtSqSum[2]; 
  edge_F_0_incr[3] = 0.5*nuVtSqSum[3]*temp_F_0_edge[8]+0.5*nuVtSqSum[2]*temp_F_0_edge[6]+0.5*nuVtSqSum[1]*temp_F_0_edge[5]+0.5*nuVtSqSum[0]*temp_F_0_edge[3]; 
  edge_F_0_incr[4] = 0.5*nuVtSqSum[0]*temp_F_0_edge[4]+0.5*temp_F_0_edge[0]*nuVtSqSum[3]+0.5*nuVtSqSum[1]*temp_F_0_edge[2]+0.5*temp_F_0_edge[1]*nuVtSqSum[2]; 
  edge_F_0_incr[5] = 0.5*nuVtSqSum[2]*temp_F_0_edge[8]+0.5*nuVtSqSum[3]*temp_F_0_edge[6]+0.5*nuVtSqSum[0]*temp_F_0_edge[5]+0.5*nuVtSqSum[1]*temp_F_0_edge[3]; 
  edge_F_0_incr[6] = 0.5*nuVtSqSum[1]*temp_F_0_edge[8]+0.5*nuVtSqSum[0]*temp_F_0_edge[6]+0.5*nuVtSqSum[3]*temp_F_0_edge[5]+0.5*nuVtSqSum[2]*temp_F_0_edge[3]; 
  edge_F_0_incr[7] = 0.5*nuVtSqSum[3]*temp_F_0_edge[11]+0.5000000000000001*nuVtSqSum[2]*temp_F_0_edge[10]+0.5000000000000001*nuVtSqSum[1]*temp_F_0_edge[9]+0.5*nuVtSqSum[0]*temp_F_0_edge[7]; 
  edge_F_0_incr[8] = 0.5*nuVtSqSum[0]*temp_F_0_edge[8]+0.5*nuVtSqSum[1]*temp_F_0_edge[6]+0.5*nuVtSqSum[2]*temp_F_0_edge[5]+0.5*nuVtSqSum[3]*temp_F_0_edge[3]; 
  edge_F_0_incr[9] = 0.5000000000000001*nuVtSqSum[2]*temp_F_0_edge[11]+0.5*nuVtSqSum[3]*temp_F_0_edge[10]+0.5*nuVtSqSum[0]*temp_F_0_edge[9]+0.5000000000000001*nuVtSqSum[1]*temp_F_0_edge[7]; 
  edge_F_0_incr[10] = 0.5000000000000001*nuVtSqSum[1]*temp_F_0_edge[11]+0.5*nuVtSqSum[0]*temp_F_0_edge[10]+0.5*nuVtSqSum[3]*temp_F_0_edge[9]+0.5000000000000001*nuVtSqSum[2]*temp_F_0_edge[7]; 
  edge_F_0_incr[11] = 0.5*nuVtSqSum[0]*temp_F_0_edge[11]+0.5000000000000001*nuVtSqSum[1]*temp_F_0_edge[10]+0.5000000000000001*nuVtSqSum[2]*temp_F_0_edge[9]+0.5*nuVtSqSum[3]*temp_F_0_edge[7]; 
  edge_G_1_incr[0] = 0.5*nuVtSqSum[3]*temp_G_1_edge[4]+0.5*nuVtSqSum[2]*temp_G_1_edge[2]+0.5*nuVtSqSum[1]*temp_G_1_edge[1]+0.5*nuVtSqSum[0]*temp_G_1_edge[0]; 
  edge_G_1_incr[1] = 0.5*nuVtSqSum[2]*temp_G_1_edge[4]+0.5*temp_G_1_edge[2]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_G_1_edge[1]+0.5*temp_G_1_edge[0]*nuVtSqSum[1]; 
  edge_G_1_incr[2] = 0.5*nuVtSqSum[1]*temp_G_1_edge[4]+0.5*temp_G_1_edge[1]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_G_1_edge[2]+0.5*temp_G_1_edge[0]*nuVtSqSum[2]; 
  edge_G_1_incr[3] = 0.5*nuVtSqSum[3]*temp_G_1_edge[8]+0.5*nuVtSqSum[2]*temp_G_1_edge[6]+0.5*nuVtSqSum[1]*temp_G_1_edge[5]+0.5*nuVtSqSum[0]*temp_G_1_edge[3]; 
  edge_G_1_incr[4] = 0.5*nuVtSqSum[0]*temp_G_1_edge[4]+0.5*temp_G_1_edge[0]*nuVtSqSum[3]+0.5*nuVtSqSum[1]*temp_G_1_edge[2]+0.5*temp_G_1_edge[1]*nuVtSqSum[2]; 
  edge_G_1_incr[5] = 0.5*nuVtSqSum[2]*temp_G_1_edge[8]+0.5*nuVtSqSum[3]*temp_G_1_edge[6]+0.5*nuVtSqSum[0]*temp_G_1_edge[5]+0.5*nuVtSqSum[1]*temp_G_1_edge[3]; 
  edge_G_1_incr[6] = 0.5*nuVtSqSum[1]*temp_G_1_edge[8]+0.5*nuVtSqSum[0]*temp_G_1_edge[6]+0.5*nuVtSqSum[3]*temp_G_1_edge[5]+0.5*nuVtSqSum[2]*temp_G_1_edge[3]; 
  edge_G_1_incr[7] = 0.5*nuVtSqSum[3]*temp_G_1_edge[11]+0.5000000000000001*nuVtSqSum[2]*temp_G_1_edge[10]+0.5000000000000001*nuVtSqSum[1]*temp_G_1_edge[9]+0.5*nuVtSqSum[0]*temp_G_1_edge[7]; 
  edge_G_1_incr[8] = 0.5*nuVtSqSum[0]*temp_G_1_edge[8]+0.5*nuVtSqSum[1]*temp_G_1_edge[6]+0.5*nuVtSqSum[2]*temp_G_1_edge[5]+0.5*nuVtSqSum[3]*temp_G_1_edge[3]; 
  edge_G_1_incr[9] = 0.5000000000000001*nuVtSqSum[2]*temp_G_1_edge[11]+0.5*nuVtSqSum[3]*temp_G_1_edge[10]+0.5*nuVtSqSum[0]*temp_G_1_edge[9]+0.5000000000000001*nuVtSqSum[1]*temp_G_1_edge[7]; 
  edge_G_1_incr[10] = 0.5000000000000001*nuVtSqSum[1]*temp_G_1_edge[11]+0.5*nuVtSqSum[0]*temp_G_1_edge[10]+0.5*nuVtSqSum[3]*temp_G_1_edge[9]+0.5000000000000001*nuVtSqSum[2]*temp_G_1_edge[7]; 
  edge_G_1_incr[11] = 0.5*nuVtSqSum[0]*temp_G_1_edge[11]+0.5000000000000001*nuVtSqSum[1]*temp_G_1_edge[10]+0.5000000000000001*nuVtSqSum[2]*temp_G_1_edge[9]+0.5*nuVtSqSum[3]*temp_G_1_edge[7]; 

  } 

  out_F_0[0] += edge_F_0_incr[0]*rdvSq4+diff_F_0_incr[0]*rdvSq4+vol_incr_F_0[0]; 
  out_F_0[1] += edge_F_0_incr[1]*rdvSq4+diff_F_0_incr[1]*rdvSq4+vol_incr_F_0[1]; 
  out_F_0[2] += edge_F_0_incr[2]*rdvSq4+diff_F_0_incr[2]*rdvSq4+vol_incr_F_0[2]; 
  out_F_0[3] += edge_F_0_incr[3]*rdvSq4+diff_F_0_incr[3]*rdvSq4+vol_incr_F_0[3]; 
  out_F_0[4] += edge_F_0_incr[4]*rdvSq4+diff_F_0_incr[4]*rdvSq4+vol_incr_F_0[4]; 
  out_F_0[5] += edge_F_0_incr[5]*rdvSq4+diff_F_0_incr[5]*rdvSq4+vol_incr_F_0[5]; 
  out_F_0[6] += edge_F_0_incr[6]*rdvSq4+diff_F_0_incr[6]*rdvSq4+vol_incr_F_0[6]; 
  out_F_0[7] += edge_F_0_incr[7]*rdvSq4+diff_F_0_incr[7]*rdvSq4+vol_incr_F_0[7]; 
  out_F_0[8] += edge_F_0_incr[8]*rdvSq4+diff_F_0_incr[8]*rdvSq4+vol_incr_F_0[8]; 
  out_F_0[9] += edge_F_0_incr[9]*rdvSq4+diff_F_0_incr[9]*rdvSq4+vol_incr_F_0[9]; 
  out_F_0[10] += edge_F_0_incr[10]*rdvSq4+diff_F_0_incr[10]*rdvSq4+vol_incr_F_0[10]; 
  out_F_0[11] += edge_F_0_incr[11]*rdvSq4+diff_F_0_incr[11]*rdvSq4+vol_incr_F_0[11]; 
  out_G_1[0] += edge_G_1_incr[0]*rdvSq4+diff_G_1_incr[0]*rdvSq4+vol_incr_G_1[0]; 
  out_G_1[1] += edge_G_1_incr[1]*rdvSq4+diff_G_1_incr[1]*rdvSq4+vol_incr_G_1[1]; 
  out_G_1[2] += edge_G_1_incr[2]*rdvSq4+diff_G_1_incr[2]*rdvSq4+vol_incr_G_1[2]; 
  out_G_1[3] += edge_G_1_incr[3]*rdvSq4+diff_G_1_incr[3]*rdvSq4+vol_incr_G_1[3]; 
  out_G_1[4] += edge_G_1_incr[4]*rdvSq4+diff_G_1_incr[4]*rdvSq4+vol_incr_G_1[4]; 
  out_G_1[5] += edge_G_1_incr[5]*rdvSq4+diff_G_1_incr[5]*rdvSq4+vol_incr_G_1[5]; 
  out_G_1[6] += edge_G_1_incr[6]*rdvSq4+diff_G_1_incr[6]*rdvSq4+vol_incr_G_1[6]; 
  out_G_1[7] += edge_G_1_incr[7]*rdvSq4+diff_G_1_incr[7]*rdvSq4+vol_incr_G_1[7]; 
  out_G_1[8] += edge_G_1_incr[8]*rdvSq4+diff_G_1_incr[8]*rdvSq4+vol_incr_G_1[8]; 
  out_G_1[9] += edge_G_1_incr[9]*rdvSq4+diff_G_1_incr[9]*rdvSq4+vol_incr_G_1[9]; 
  out_G_1[10] += edge_G_1_incr[10]*rdvSq4+diff_G_1_incr[10]*rdvSq4+vol_incr_G_1[10]; 
  out_G_1[11] += edge_G_1_incr[11]*rdvSq4+diff_G_1_incr[11]*rdvSq4+vol_incr_G_1[11]; 

  return 0.;

} 
