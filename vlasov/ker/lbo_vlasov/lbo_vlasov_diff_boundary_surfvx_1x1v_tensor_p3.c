#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvx_1x1v_tensor_p3(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[2]: Cell-center coordinates. 
  // dxv[2]: Cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[8]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fSkin/Edge: Distribution function in cells 
  // out: Incremented distribution function in cell 
  const double *nuVtSqSum = &nuPrimMomsSum[4];

  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 

  double facDiff[4]; 
  // Expand nuVtSqSum in phase basis.
  facDiff[0] = nuVtSqSum[0]; 
  facDiff[1] = nuVtSqSum[1]; 
  facDiff[2] = nuVtSqSum[2]; 
  facDiff[3] = nuVtSqSum[3]; 

  double vol_incr[16] = {0.0}; 
  vol_incr[5] = 4.743416490252569*facDiff[3]*fSkin[8]+4.743416490252569*facDiff[2]*fSkin[4]+4.743416490252569*fSkin[1]*facDiff[1]+4.743416490252569*fSkin[0]*facDiff[0]; 
  vol_incr[7] = 4.16619044897648*facDiff[2]*fSkin[8]+4.16619044897648*facDiff[3]*fSkin[4]+4.242640687119286*facDiff[1]*fSkin[4]+4.242640687119286*fSkin[1]*facDiff[2]+4.743416490252569*fSkin[0]*facDiff[1]+4.743416490252569*facDiff[0]*fSkin[1]; 
  vol_incr[9] = 16.20185174601965*facDiff[3]*fSkin[11]+16.201851746019653*facDiff[2]*fSkin[6]+16.20185174601965*facDiff[1]*fSkin[3]+16.20185174601965*facDiff[0]*fSkin[2]; 
  vol_incr[10] = 2.8284271247461907*facDiff[3]*fSkin[8]+4.166190448976479*facDiff[1]*fSkin[8]+3.0304576336566322*facDiff[2]*fSkin[4]+4.743416490252569*facDiff[0]*fSkin[4]+4.166190448976479*fSkin[1]*facDiff[3]+4.743416490252569*fSkin[0]*facDiff[2]+4.242640687119286*fSkin[1]*facDiff[1]; 
  vol_incr[12] = 14.230249470757707*facDiff[2]*fSkin[11]+14.230249470757704*facDiff[3]*fSkin[6]+14.491376746189442*facDiff[1]*fSkin[6]+14.491376746189438*facDiff[2]*fSkin[3]+16.20185174601965*facDiff[0]*fSkin[3]+16.20185174601965*facDiff[1]*fSkin[2]; 
  vol_incr[13] = 2.8284271247461903*facDiff[2]*fSkin[8]+4.743416490252568*facDiff[0]*fSkin[8]+2.8284271247461903*facDiff[3]*fSkin[4]+4.166190448976479*facDiff[1]*fSkin[4]+4.743416490252568*fSkin[0]*facDiff[3]+4.166190448976479*fSkin[1]*facDiff[2]; 
  vol_incr[14] = 9.660917830792961*facDiff[3]*fSkin[11]+14.230249470757702*facDiff[1]*fSkin[11]+10.350983390135315*facDiff[2]*fSkin[6]+16.201851746019653*facDiff[0]*fSkin[6]+14.230249470757705*fSkin[3]*facDiff[3]+14.49137674618944*facDiff[1]*fSkin[3]+16.201851746019653*fSkin[2]*facDiff[2]; 
  vol_incr[15] = 9.660917830792961*facDiff[2]*fSkin[11]+16.20185174601965*facDiff[0]*fSkin[11]+9.660917830792963*facDiff[3]*fSkin[6]+14.230249470757704*facDiff[1]*fSkin[6]+16.20185174601965*fSkin[2]*facDiff[3]+14.230249470757707*facDiff[2]*fSkin[3]; 

  double edgeSurf_incr[16] = {0.0}; 
  double boundSurf_incr[16] = {0.0}; 

  if (edge == -1) { 

  double edgeSurf[16] = {0.0}; 
  edgeSurf[0] = -(0.6821077598838398*fSkin[9])-0.6821077598838398*fEdge[9]-1.5547660156053227*fSkin[5]+1.5547660156053227*fEdge[5]-1.935025511580855*fSkin[2]-1.935025511580855*fEdge[2]-1.3671875*fSkin[0]+1.3671875*fEdge[0]; 
  edgeSurf[1] = -(0.6821077598838398*fSkin[12])-0.6821077598838398*fEdge[12]-1.5547660156053227*fSkin[7]+1.5547660156053227*fEdge[7]-1.935025511580855*fSkin[3]-1.935025511580855*fEdge[3]-1.3671875*fSkin[1]+1.3671875*fEdge[1]; 
  edgeSurf[2] = -(1.5087721317097909*fSkin[9])-0.8541184610018139*fEdge[9]-3.2980873807547537*fSkin[5]+2.0877800850649355*fEdge[5]-4.0078125*fSkin[2]-2.6953125*fEdge[2]-2.801050915365293*fSkin[0]+1.935025511580855*fEdge[0]; 
  edgeSurf[3] = -(1.508772131709791*fSkin[12])-0.8541184610018139*fEdge[12]-3.2980873807547537*fSkin[7]+2.087780085064936*fEdge[7]-4.0078125*fSkin[3]-2.6953125*fEdge[3]-2.801050915365293*fSkin[1]+1.935025511580855*fEdge[1]; 
  edgeSurf[4] = -(0.6821077598838398*fSkin[14])-0.6821077598838398*fEdge[14]-1.5547660156053227*fSkin[10]+1.5547660156053227*fEdge[10]-1.935025511580855*fSkin[6]-1.935025511580855*fEdge[6]-1.3671875*fSkin[4]+1.3671875*fEdge[4]; 
  edgeSurf[5] = -(2.792970701173145*fSkin[9])-0.2575079369875949*fEdge[9]-5.8203125*fSkin[5]+1.1328125*fEdge[5]-6.868493903039716*fSkin[2]-1.7852032611424813*fEdge[2]-4.734175171112836*fSkin[0]+1.3800732048631514*fEdge[0]; 
  edgeSurf[6] = -(1.508772131709791*fSkin[14])-0.8541184610018139*fEdge[14]-3.2980873807547537*fSkin[10]+2.087780085064936*fEdge[10]-4.0078125*fSkin[6]-2.6953125*fEdge[6]-2.8010509153652934*fSkin[4]+1.935025511580855*fEdge[4]; 
  edgeSurf[7] = -(2.792970701173145*fSkin[12])-0.25750793698759494*fEdge[12]-5.8203125*fSkin[7]+1.1328125*fEdge[7]-6.868493903039716*fSkin[3]-1.7852032611424813*fEdge[3]-4.734175171112837*fSkin[1]+1.3800732048631517*fEdge[1]; 
  edgeSurf[8] = -(0.6821077598838398*fSkin[15])-0.6821077598838398*fEdge[15]-1.5547660156053225*fSkin[13]+1.5547660156053225*fEdge[13]-1.9350255115808548*fSkin[11]-1.9350255115808548*fEdge[11]-1.3671875*fSkin[8]+1.3671875*fEdge[8]; 
  edgeSurf[9] = -(4.8046875*fSkin[9])+1.1953125*fEdge[9]-9.659849020842342*fSkin[5]-1.4328005724694384*fEdge[5]-11.134226883838016*fSkin[2]+0.8950343154210625*fEdge[2]-7.585865087193006*fSkin[0]-0.351388846000766*fEdge[0]; 
  edgeSurf[10] = -(2.792970701173145*fSkin[14])-0.2575079369875949*fEdge[14]-5.8203125*fSkin[10]+1.1328125*fEdge[10]-6.868493903039716*fSkin[6]-1.7852032611424813*fEdge[6]-4.734175171112836*fSkin[4]+1.3800732048631514*fEdge[4]; 
  edgeSurf[11] = -(1.508772131709791*fSkin[15])-0.8541184610018139*fEdge[15]-3.2980873807547533*fSkin[13]+2.087780085064936*fEdge[13]-4.0078125*fSkin[11]-2.6953125*fEdge[11]-2.8010509153652934*fSkin[8]+1.9350255115808548*fEdge[8]; 
  edgeSurf[12] = -(4.8046875*fSkin[12])+1.1953125*fEdge[12]-9.659849020842342*fSkin[7]-1.4328005724694384*fEdge[7]-11.134226883838018*fSkin[3]+0.8950343154210625*fEdge[3]-7.585865087193007*fSkin[1]-0.351388846000766*fEdge[1]; 
  edgeSurf[13] = -(2.792970701173145*fSkin[15])-0.2575079369875949*fEdge[15]-5.8203125*fSkin[13]+1.1328125*fEdge[13]-6.868493903039717*fSkin[11]-1.7852032611424813*fEdge[11]-4.7341751711128355*fSkin[8]+1.3800732048631514*fEdge[8]; 
  edgeSurf[14] = -(4.8046875*fSkin[14])+1.1953125*fEdge[14]-9.659849020842342*fSkin[10]-1.4328005724694384*fEdge[10]-11.13422688383802*fSkin[6]+0.8950343154210625*fEdge[6]-7.585865087193006*fSkin[4]-0.351388846000766*fEdge[4]; 
  edgeSurf[15] = -(4.8046875*fSkin[15])+1.1953125*fEdge[15]-9.659849020842342*fSkin[13]-1.4328005724694384*fEdge[13]-11.134226883838018*fSkin[11]+0.8950343154210625*fEdge[11]-7.585865087193006*fSkin[8]-0.351388846000766*fEdge[8]; 

  double boundSurf[16] = {0.0}; 
  boundSurf[2] = -(2.29128784747792*fSkin[9])+1.9364916731037085*fSkin[5]-1.5*fSkin[2]+0.8660254037844386*fSkin[0]; 
  boundSurf[3] = -(2.29128784747792*fSkin[12])+1.9364916731037085*fSkin[7]-1.5*fSkin[3]+0.8660254037844386*fSkin[1]; 
  boundSurf[5] = 8.874119674649425*fSkin[9]-7.5*fSkin[5]+5.809475019311125*fSkin[2]-3.3541019662496847*fSkin[0]; 
  boundSurf[6] = -(2.29128784747792*fSkin[14])+1.9364916731037085*fSkin[10]-1.5*fSkin[6]+0.8660254037844387*fSkin[4]; 
  boundSurf[7] = 8.874119674649425*fSkin[12]-7.5*fSkin[7]+5.809475019311126*fSkin[3]-3.3541019662496843*fSkin[1]; 
  boundSurf[9] = -(21.0*fSkin[9])+17.74823934929885*fSkin[5]-13.747727084867519*fSkin[2]+7.937253933193772*fSkin[0]; 
  boundSurf[10] = 8.874119674649425*fSkin[14]-7.5*fSkin[10]+5.809475019311126*fSkin[6]-3.3541019662496847*fSkin[4]; 
  boundSurf[11] = -(2.29128784747792*fSkin[15])+1.9364916731037083*fSkin[13]-1.5*fSkin[11]+0.8660254037844386*fSkin[8]; 
  boundSurf[12] = -(21.0*fSkin[12])+17.74823934929885*fSkin[7]-13.74772708486752*fSkin[3]+7.937253933193771*fSkin[1]; 
  boundSurf[13] = 8.874119674649425*fSkin[15]-7.5*fSkin[13]+5.809475019311125*fSkin[11]-3.354101966249684*fSkin[8]; 
  boundSurf[14] = -(21.0*fSkin[14])+17.74823934929885*fSkin[10]-13.74772708486752*fSkin[6]+7.937253933193772*fSkin[4]; 
  boundSurf[15] = -(21.0*fSkin[15])+17.74823934929885*fSkin[13]-13.74772708486752*fSkin[11]+7.937253933193772*fSkin[8]; 

  edgeSurf_incr[0] = 0.7071067811865475*nuVtSqSum[3]*edgeSurf[8]+0.7071067811865475*nuVtSqSum[2]*edgeSurf[4]+0.7071067811865475*edgeSurf[1]*nuVtSqSum[1]+0.7071067811865475*edgeSurf[0]*nuVtSqSum[0]; 
  edgeSurf_incr[1] = 0.6210590034081186*nuVtSqSum[2]*edgeSurf[8]+0.6210590034081186*nuVtSqSum[3]*edgeSurf[4]+0.6324555320336759*nuVtSqSum[1]*edgeSurf[4]+0.6324555320336759*edgeSurf[1]*nuVtSqSum[2]+0.7071067811865475*edgeSurf[0]*nuVtSqSum[1]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[1]; 
  edgeSurf_incr[2] = 0.7071067811865474*nuVtSqSum[3]*edgeSurf[11]+0.7071067811865475*nuVtSqSum[2]*edgeSurf[6]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[3]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[2]; 
  edgeSurf_incr[3] = 0.6210590034081187*nuVtSqSum[2]*edgeSurf[11]+0.6210590034081187*nuVtSqSum[3]*edgeSurf[6]+0.632455532033676*nuVtSqSum[1]*edgeSurf[6]+0.6324555320336759*nuVtSqSum[2]*edgeSurf[3]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[3]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[2]; 
  edgeSurf_incr[4] = 0.421637021355784*nuVtSqSum[3]*edgeSurf[8]+0.6210590034081186*nuVtSqSum[1]*edgeSurf[8]+0.45175395145262565*nuVtSqSum[2]*edgeSurf[4]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[4]+0.6210590034081186*edgeSurf[1]*nuVtSqSum[3]+0.7071067811865475*edgeSurf[0]*nuVtSqSum[2]+0.6324555320336759*edgeSurf[1]*nuVtSqSum[1]; 
  edgeSurf_incr[5] = 0.7071067811865474*nuVtSqSum[3]*edgeSurf[13]+0.7071067811865475*nuVtSqSum[2]*edgeSurf[10]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[7]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[5]; 
  edgeSurf_incr[6] = 0.4216370213557839*nuVtSqSum[3]*edgeSurf[11]+0.6210590034081187*nuVtSqSum[1]*edgeSurf[11]+0.45175395145262565*nuVtSqSum[2]*edgeSurf[6]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[6]+0.6210590034081187*edgeSurf[3]*nuVtSqSum[3]+0.632455532033676*nuVtSqSum[1]*edgeSurf[3]+0.7071067811865475*edgeSurf[2]*nuVtSqSum[2]; 
  edgeSurf_incr[7] = 0.6210590034081186*nuVtSqSum[2]*edgeSurf[13]+0.6210590034081188*nuVtSqSum[3]*edgeSurf[10]+0.632455532033676*nuVtSqSum[1]*edgeSurf[10]+0.6324555320336759*nuVtSqSum[2]*edgeSurf[7]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[7]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[5]; 
  edgeSurf_incr[8] = 0.421637021355784*nuVtSqSum[2]*edgeSurf[8]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[8]+0.421637021355784*nuVtSqSum[3]*edgeSurf[4]+0.6210590034081186*nuVtSqSum[1]*edgeSurf[4]+0.7071067811865475*edgeSurf[0]*nuVtSqSum[3]+0.6210590034081186*edgeSurf[1]*nuVtSqSum[2]; 
  edgeSurf_incr[9] = 0.7071067811865475*nuVtSqSum[3]*edgeSurf[15]+0.7071067811865474*nuVtSqSum[2]*edgeSurf[14]+0.7071067811865474*nuVtSqSum[1]*edgeSurf[12]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[9]; 
  edgeSurf_incr[10] = 0.42163702135578396*nuVtSqSum[3]*edgeSurf[13]+0.6210590034081186*nuVtSqSum[1]*edgeSurf[13]+0.45175395145262565*nuVtSqSum[2]*edgeSurf[10]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[10]+0.6210590034081187*nuVtSqSum[3]*edgeSurf[7]+0.632455532033676*nuVtSqSum[1]*edgeSurf[7]+0.7071067811865475*nuVtSqSum[2]*edgeSurf[5]; 
  edgeSurf_incr[11] = 0.421637021355784*nuVtSqSum[2]*edgeSurf[11]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[11]+0.4216370213557839*nuVtSqSum[3]*edgeSurf[6]+0.6210590034081187*nuVtSqSum[1]*edgeSurf[6]+0.7071067811865474*edgeSurf[2]*nuVtSqSum[3]+0.6210590034081187*nuVtSqSum[2]*edgeSurf[3]; 
  edgeSurf_incr[12] = 0.6210590034081186*nuVtSqSum[2]*edgeSurf[15]+0.6210590034081186*nuVtSqSum[3]*edgeSurf[14]+0.6324555320336759*nuVtSqSum[1]*edgeSurf[14]+0.6324555320336759*nuVtSqSum[2]*edgeSurf[12]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[12]+0.7071067811865474*nuVtSqSum[1]*edgeSurf[9]; 
  edgeSurf_incr[13] = 0.421637021355784*nuVtSqSum[2]*edgeSurf[13]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[13]+0.42163702135578396*nuVtSqSum[3]*edgeSurf[10]+0.6210590034081186*nuVtSqSum[1]*edgeSurf[10]+0.6210590034081186*nuVtSqSum[2]*edgeSurf[7]+0.7071067811865474*nuVtSqSum[3]*edgeSurf[5]; 
  edgeSurf_incr[14] = 0.42163702135578396*nuVtSqSum[3]*edgeSurf[15]+0.6210590034081186*nuVtSqSum[1]*edgeSurf[15]+0.45175395145262565*nuVtSqSum[2]*edgeSurf[14]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[14]+0.6210590034081186*nuVtSqSum[3]*edgeSurf[12]+0.6324555320336759*nuVtSqSum[1]*edgeSurf[12]+0.7071067811865474*nuVtSqSum[2]*edgeSurf[9]; 
  edgeSurf_incr[15] = 0.421637021355784*nuVtSqSum[2]*edgeSurf[15]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[15]+0.42163702135578396*nuVtSqSum[3]*edgeSurf[14]+0.6210590034081186*nuVtSqSum[1]*edgeSurf[14]+0.6210590034081187*nuVtSqSum[2]*edgeSurf[12]+0.7071067811865475*nuVtSqSum[3]*edgeSurf[9]; 

  boundSurf_incr[0] = 0.7071067811865475*nuVtSqSum[3]*boundSurf[8]+0.7071067811865475*nuVtSqSum[2]*boundSurf[4]+0.7071067811865475*boundSurf[1]*nuVtSqSum[1]+0.7071067811865475*boundSurf[0]*nuVtSqSum[0]; 
  boundSurf_incr[1] = 0.6210590034081186*nuVtSqSum[2]*boundSurf[8]+0.6210590034081186*nuVtSqSum[3]*boundSurf[4]+0.6324555320336759*nuVtSqSum[1]*boundSurf[4]+0.6324555320336759*boundSurf[1]*nuVtSqSum[2]+0.7071067811865475*boundSurf[0]*nuVtSqSum[1]+0.7071067811865475*nuVtSqSum[0]*boundSurf[1]; 
  boundSurf_incr[2] = 0.7071067811865474*nuVtSqSum[3]*boundSurf[11]+0.7071067811865475*nuVtSqSum[2]*boundSurf[6]+0.7071067811865475*nuVtSqSum[1]*boundSurf[3]+0.7071067811865475*nuVtSqSum[0]*boundSurf[2]; 
  boundSurf_incr[3] = 0.6210590034081187*nuVtSqSum[2]*boundSurf[11]+0.6210590034081187*nuVtSqSum[3]*boundSurf[6]+0.632455532033676*nuVtSqSum[1]*boundSurf[6]+0.6324555320336759*nuVtSqSum[2]*boundSurf[3]+0.7071067811865475*nuVtSqSum[0]*boundSurf[3]+0.7071067811865475*nuVtSqSum[1]*boundSurf[2]; 
  boundSurf_incr[4] = 0.421637021355784*nuVtSqSum[3]*boundSurf[8]+0.6210590034081186*nuVtSqSum[1]*boundSurf[8]+0.45175395145262565*nuVtSqSum[2]*boundSurf[4]+0.7071067811865475*nuVtSqSum[0]*boundSurf[4]+0.6210590034081186*boundSurf[1]*nuVtSqSum[3]+0.7071067811865475*boundSurf[0]*nuVtSqSum[2]+0.6324555320336759*boundSurf[1]*nuVtSqSum[1]; 
  boundSurf_incr[5] = 0.7071067811865474*nuVtSqSum[3]*boundSurf[13]+0.7071067811865475*nuVtSqSum[2]*boundSurf[10]+0.7071067811865475*nuVtSqSum[1]*boundSurf[7]+0.7071067811865475*nuVtSqSum[0]*boundSurf[5]; 
  boundSurf_incr[6] = 0.4216370213557839*nuVtSqSum[3]*boundSurf[11]+0.6210590034081187*nuVtSqSum[1]*boundSurf[11]+0.45175395145262565*nuVtSqSum[2]*boundSurf[6]+0.7071067811865475*nuVtSqSum[0]*boundSurf[6]+0.6210590034081187*boundSurf[3]*nuVtSqSum[3]+0.632455532033676*nuVtSqSum[1]*boundSurf[3]+0.7071067811865475*boundSurf[2]*nuVtSqSum[2]; 
  boundSurf_incr[7] = 0.6210590034081186*nuVtSqSum[2]*boundSurf[13]+0.6210590034081188*nuVtSqSum[3]*boundSurf[10]+0.632455532033676*nuVtSqSum[1]*boundSurf[10]+0.6324555320336759*nuVtSqSum[2]*boundSurf[7]+0.7071067811865475*nuVtSqSum[0]*boundSurf[7]+0.7071067811865475*nuVtSqSum[1]*boundSurf[5]; 
  boundSurf_incr[8] = 0.421637021355784*nuVtSqSum[2]*boundSurf[8]+0.7071067811865475*nuVtSqSum[0]*boundSurf[8]+0.421637021355784*nuVtSqSum[3]*boundSurf[4]+0.6210590034081186*nuVtSqSum[1]*boundSurf[4]+0.7071067811865475*boundSurf[0]*nuVtSqSum[3]+0.6210590034081186*boundSurf[1]*nuVtSqSum[2]; 
  boundSurf_incr[9] = 0.7071067811865475*nuVtSqSum[3]*boundSurf[15]+0.7071067811865474*nuVtSqSum[2]*boundSurf[14]+0.7071067811865474*nuVtSqSum[1]*boundSurf[12]+0.7071067811865475*nuVtSqSum[0]*boundSurf[9]; 
  boundSurf_incr[10] = 0.42163702135578396*nuVtSqSum[3]*boundSurf[13]+0.6210590034081186*nuVtSqSum[1]*boundSurf[13]+0.45175395145262565*nuVtSqSum[2]*boundSurf[10]+0.7071067811865475*nuVtSqSum[0]*boundSurf[10]+0.6210590034081187*nuVtSqSum[3]*boundSurf[7]+0.632455532033676*nuVtSqSum[1]*boundSurf[7]+0.7071067811865475*nuVtSqSum[2]*boundSurf[5]; 
  boundSurf_incr[11] = 0.421637021355784*nuVtSqSum[2]*boundSurf[11]+0.7071067811865475*nuVtSqSum[0]*boundSurf[11]+0.4216370213557839*nuVtSqSum[3]*boundSurf[6]+0.6210590034081187*nuVtSqSum[1]*boundSurf[6]+0.7071067811865474*boundSurf[2]*nuVtSqSum[3]+0.6210590034081187*nuVtSqSum[2]*boundSurf[3]; 
  boundSurf_incr[12] = 0.6210590034081186*nuVtSqSum[2]*boundSurf[15]+0.6210590034081186*nuVtSqSum[3]*boundSurf[14]+0.6324555320336759*nuVtSqSum[1]*boundSurf[14]+0.6324555320336759*nuVtSqSum[2]*boundSurf[12]+0.7071067811865475*nuVtSqSum[0]*boundSurf[12]+0.7071067811865474*nuVtSqSum[1]*boundSurf[9]; 
  boundSurf_incr[13] = 0.421637021355784*nuVtSqSum[2]*boundSurf[13]+0.7071067811865475*nuVtSqSum[0]*boundSurf[13]+0.42163702135578396*nuVtSqSum[3]*boundSurf[10]+0.6210590034081186*nuVtSqSum[1]*boundSurf[10]+0.6210590034081186*nuVtSqSum[2]*boundSurf[7]+0.7071067811865474*nuVtSqSum[3]*boundSurf[5]; 
  boundSurf_incr[14] = 0.42163702135578396*nuVtSqSum[3]*boundSurf[15]+0.6210590034081186*nuVtSqSum[1]*boundSurf[15]+0.45175395145262565*nuVtSqSum[2]*boundSurf[14]+0.7071067811865475*nuVtSqSum[0]*boundSurf[14]+0.6210590034081186*nuVtSqSum[3]*boundSurf[12]+0.6324555320336759*nuVtSqSum[1]*boundSurf[12]+0.7071067811865474*nuVtSqSum[2]*boundSurf[9]; 
  boundSurf_incr[15] = 0.421637021355784*nuVtSqSum[2]*boundSurf[15]+0.7071067811865475*nuVtSqSum[0]*boundSurf[15]+0.42163702135578396*nuVtSqSum[3]*boundSurf[14]+0.6210590034081186*nuVtSqSum[1]*boundSurf[14]+0.6210590034081187*nuVtSqSum[2]*boundSurf[12]+0.7071067811865475*nuVtSqSum[3]*boundSurf[9]; 


  } else { 

  double edgeSurf[16] = {0.0}; 
  edgeSurf[0] = 0.6821077598838398*fSkin[9]+0.6821077598838398*fEdge[9]-1.5547660156053227*fSkin[5]+1.5547660156053227*fEdge[5]+1.935025511580855*fSkin[2]+1.935025511580855*fEdge[2]-1.3671875*fSkin[0]+1.3671875*fEdge[0]; 
  edgeSurf[1] = 0.6821077598838398*fSkin[12]+0.6821077598838398*fEdge[12]-1.5547660156053227*fSkin[7]+1.5547660156053227*fEdge[7]+1.935025511580855*fSkin[3]+1.935025511580855*fEdge[3]-1.3671875*fSkin[1]+1.3671875*fEdge[1]; 
  edgeSurf[2] = -(1.5087721317097909*fSkin[9])-0.8541184610018139*fEdge[9]+3.2980873807547537*fSkin[5]-2.0877800850649355*fEdge[5]-4.0078125*fSkin[2]-2.6953125*fEdge[2]+2.801050915365293*fSkin[0]-1.935025511580855*fEdge[0]; 
  edgeSurf[3] = -(1.508772131709791*fSkin[12])-0.8541184610018139*fEdge[12]+3.2980873807547537*fSkin[7]-2.087780085064936*fEdge[7]-4.0078125*fSkin[3]-2.6953125*fEdge[3]+2.801050915365293*fSkin[1]-1.935025511580855*fEdge[1]; 
  edgeSurf[4] = 0.6821077598838398*fSkin[14]+0.6821077598838398*fEdge[14]-1.5547660156053227*fSkin[10]+1.5547660156053227*fEdge[10]+1.935025511580855*fSkin[6]+1.935025511580855*fEdge[6]-1.3671875*fSkin[4]+1.3671875*fEdge[4]; 
  edgeSurf[5] = 2.792970701173145*fSkin[9]+0.2575079369875949*fEdge[9]-5.8203125*fSkin[5]+1.1328125*fEdge[5]+6.868493903039716*fSkin[2]+1.7852032611424813*fEdge[2]-4.734175171112836*fSkin[0]+1.3800732048631514*fEdge[0]; 
  edgeSurf[6] = -(1.508772131709791*fSkin[14])-0.8541184610018139*fEdge[14]+3.2980873807547537*fSkin[10]-2.087780085064936*fEdge[10]-4.0078125*fSkin[6]-2.6953125*fEdge[6]+2.8010509153652934*fSkin[4]-1.935025511580855*fEdge[4]; 
  edgeSurf[7] = 2.792970701173145*fSkin[12]+0.25750793698759494*fEdge[12]-5.8203125*fSkin[7]+1.1328125*fEdge[7]+6.868493903039716*fSkin[3]+1.7852032611424813*fEdge[3]-4.734175171112837*fSkin[1]+1.3800732048631517*fEdge[1]; 
  edgeSurf[8] = 0.6821077598838398*fSkin[15]+0.6821077598838398*fEdge[15]-1.5547660156053225*fSkin[13]+1.5547660156053225*fEdge[13]+1.9350255115808548*fSkin[11]+1.9350255115808548*fEdge[11]-1.3671875*fSkin[8]+1.3671875*fEdge[8]; 
  edgeSurf[9] = -(4.8046875*fSkin[9])+1.1953125*fEdge[9]+9.659849020842342*fSkin[5]+1.4328005724694384*fEdge[5]-11.134226883838016*fSkin[2]+0.8950343154210625*fEdge[2]+7.585865087193006*fSkin[0]+0.351388846000766*fEdge[0]; 
  edgeSurf[10] = 2.792970701173145*fSkin[14]+0.2575079369875949*fEdge[14]-5.8203125*fSkin[10]+1.1328125*fEdge[10]+6.868493903039716*fSkin[6]+1.7852032611424813*fEdge[6]-4.734175171112836*fSkin[4]+1.3800732048631514*fEdge[4]; 
  edgeSurf[11] = -(1.508772131709791*fSkin[15])-0.8541184610018139*fEdge[15]+3.2980873807547533*fSkin[13]-2.087780085064936*fEdge[13]-4.0078125*fSkin[11]-2.6953125*fEdge[11]+2.8010509153652934*fSkin[8]-1.9350255115808548*fEdge[8]; 
  edgeSurf[12] = -(4.8046875*fSkin[12])+1.1953125*fEdge[12]+9.659849020842342*fSkin[7]+1.4328005724694384*fEdge[7]-11.134226883838018*fSkin[3]+0.8950343154210625*fEdge[3]+7.585865087193007*fSkin[1]+0.351388846000766*fEdge[1]; 
  edgeSurf[13] = 2.792970701173145*fSkin[15]+0.2575079369875949*fEdge[15]-5.8203125*fSkin[13]+1.1328125*fEdge[13]+6.868493903039717*fSkin[11]+1.7852032611424813*fEdge[11]-4.7341751711128355*fSkin[8]+1.3800732048631514*fEdge[8]; 
  edgeSurf[14] = -(4.8046875*fSkin[14])+1.1953125*fEdge[14]+9.659849020842342*fSkin[10]+1.4328005724694384*fEdge[10]-11.13422688383802*fSkin[6]+0.8950343154210625*fEdge[6]+7.585865087193006*fSkin[4]+0.351388846000766*fEdge[4]; 
  edgeSurf[15] = -(4.8046875*fSkin[15])+1.1953125*fEdge[15]+9.659849020842342*fSkin[13]+1.4328005724694384*fEdge[13]-11.134226883838018*fSkin[11]+0.8950343154210625*fEdge[11]+7.585865087193006*fSkin[8]+0.351388846000766*fEdge[8]; 

  double boundSurf[16] = {0.0}; 
  boundSurf[2] = -(2.29128784747792*fSkin[9])-1.9364916731037085*fSkin[5]-1.5*fSkin[2]-0.8660254037844386*fSkin[0]; 
  boundSurf[3] = -(2.29128784747792*fSkin[12])-1.9364916731037085*fSkin[7]-1.5*fSkin[3]-0.8660254037844386*fSkin[1]; 
  boundSurf[5] = -(8.874119674649425*fSkin[9])-7.5*fSkin[5]-5.809475019311125*fSkin[2]-3.3541019662496847*fSkin[0]; 
  boundSurf[6] = -(2.29128784747792*fSkin[14])-1.9364916731037085*fSkin[10]-1.5*fSkin[6]-0.8660254037844387*fSkin[4]; 
  boundSurf[7] = -(8.874119674649425*fSkin[12])-7.5*fSkin[7]-5.809475019311126*fSkin[3]-3.3541019662496843*fSkin[1]; 
  boundSurf[9] = -(21.0*fSkin[9])-17.74823934929885*fSkin[5]-13.747727084867519*fSkin[2]-7.937253933193772*fSkin[0]; 
  boundSurf[10] = -(8.874119674649425*fSkin[14])-7.5*fSkin[10]-5.809475019311126*fSkin[6]-3.3541019662496847*fSkin[4]; 
  boundSurf[11] = -(2.29128784747792*fSkin[15])-1.9364916731037083*fSkin[13]-1.5*fSkin[11]-0.8660254037844386*fSkin[8]; 
  boundSurf[12] = -(21.0*fSkin[12])-17.74823934929885*fSkin[7]-13.74772708486752*fSkin[3]-7.937253933193771*fSkin[1]; 
  boundSurf[13] = -(8.874119674649425*fSkin[15])-7.5*fSkin[13]-5.809475019311125*fSkin[11]-3.354101966249684*fSkin[8]; 
  boundSurf[14] = -(21.0*fSkin[14])-17.74823934929885*fSkin[10]-13.74772708486752*fSkin[6]-7.937253933193772*fSkin[4]; 
  boundSurf[15] = -(21.0*fSkin[15])-17.74823934929885*fSkin[13]-13.74772708486752*fSkin[11]-7.937253933193772*fSkin[8]; 

  edgeSurf_incr[0] = 0.7071067811865475*nuVtSqSum[3]*edgeSurf[8]+0.7071067811865475*nuVtSqSum[2]*edgeSurf[4]+0.7071067811865475*edgeSurf[1]*nuVtSqSum[1]+0.7071067811865475*edgeSurf[0]*nuVtSqSum[0]; 
  edgeSurf_incr[1] = 0.6210590034081186*nuVtSqSum[2]*edgeSurf[8]+0.6210590034081186*nuVtSqSum[3]*edgeSurf[4]+0.6324555320336759*nuVtSqSum[1]*edgeSurf[4]+0.6324555320336759*edgeSurf[1]*nuVtSqSum[2]+0.7071067811865475*edgeSurf[0]*nuVtSqSum[1]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[1]; 
  edgeSurf_incr[2] = 0.7071067811865474*nuVtSqSum[3]*edgeSurf[11]+0.7071067811865475*nuVtSqSum[2]*edgeSurf[6]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[3]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[2]; 
  edgeSurf_incr[3] = 0.6210590034081187*nuVtSqSum[2]*edgeSurf[11]+0.6210590034081187*nuVtSqSum[3]*edgeSurf[6]+0.632455532033676*nuVtSqSum[1]*edgeSurf[6]+0.6324555320336759*nuVtSqSum[2]*edgeSurf[3]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[3]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[2]; 
  edgeSurf_incr[4] = 0.421637021355784*nuVtSqSum[3]*edgeSurf[8]+0.6210590034081186*nuVtSqSum[1]*edgeSurf[8]+0.45175395145262565*nuVtSqSum[2]*edgeSurf[4]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[4]+0.6210590034081186*edgeSurf[1]*nuVtSqSum[3]+0.7071067811865475*edgeSurf[0]*nuVtSqSum[2]+0.6324555320336759*edgeSurf[1]*nuVtSqSum[1]; 
  edgeSurf_incr[5] = 0.7071067811865474*nuVtSqSum[3]*edgeSurf[13]+0.7071067811865475*nuVtSqSum[2]*edgeSurf[10]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[7]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[5]; 
  edgeSurf_incr[6] = 0.4216370213557839*nuVtSqSum[3]*edgeSurf[11]+0.6210590034081187*nuVtSqSum[1]*edgeSurf[11]+0.45175395145262565*nuVtSqSum[2]*edgeSurf[6]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[6]+0.6210590034081187*edgeSurf[3]*nuVtSqSum[3]+0.632455532033676*nuVtSqSum[1]*edgeSurf[3]+0.7071067811865475*edgeSurf[2]*nuVtSqSum[2]; 
  edgeSurf_incr[7] = 0.6210590034081186*nuVtSqSum[2]*edgeSurf[13]+0.6210590034081188*nuVtSqSum[3]*edgeSurf[10]+0.632455532033676*nuVtSqSum[1]*edgeSurf[10]+0.6324555320336759*nuVtSqSum[2]*edgeSurf[7]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[7]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[5]; 
  edgeSurf_incr[8] = 0.421637021355784*nuVtSqSum[2]*edgeSurf[8]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[8]+0.421637021355784*nuVtSqSum[3]*edgeSurf[4]+0.6210590034081186*nuVtSqSum[1]*edgeSurf[4]+0.7071067811865475*edgeSurf[0]*nuVtSqSum[3]+0.6210590034081186*edgeSurf[1]*nuVtSqSum[2]; 
  edgeSurf_incr[9] = 0.7071067811865475*nuVtSqSum[3]*edgeSurf[15]+0.7071067811865474*nuVtSqSum[2]*edgeSurf[14]+0.7071067811865474*nuVtSqSum[1]*edgeSurf[12]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[9]; 
  edgeSurf_incr[10] = 0.42163702135578396*nuVtSqSum[3]*edgeSurf[13]+0.6210590034081186*nuVtSqSum[1]*edgeSurf[13]+0.45175395145262565*nuVtSqSum[2]*edgeSurf[10]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[10]+0.6210590034081187*nuVtSqSum[3]*edgeSurf[7]+0.632455532033676*nuVtSqSum[1]*edgeSurf[7]+0.7071067811865475*nuVtSqSum[2]*edgeSurf[5]; 
  edgeSurf_incr[11] = 0.421637021355784*nuVtSqSum[2]*edgeSurf[11]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[11]+0.4216370213557839*nuVtSqSum[3]*edgeSurf[6]+0.6210590034081187*nuVtSqSum[1]*edgeSurf[6]+0.7071067811865474*edgeSurf[2]*nuVtSqSum[3]+0.6210590034081187*nuVtSqSum[2]*edgeSurf[3]; 
  edgeSurf_incr[12] = 0.6210590034081186*nuVtSqSum[2]*edgeSurf[15]+0.6210590034081186*nuVtSqSum[3]*edgeSurf[14]+0.6324555320336759*nuVtSqSum[1]*edgeSurf[14]+0.6324555320336759*nuVtSqSum[2]*edgeSurf[12]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[12]+0.7071067811865474*nuVtSqSum[1]*edgeSurf[9]; 
  edgeSurf_incr[13] = 0.421637021355784*nuVtSqSum[2]*edgeSurf[13]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[13]+0.42163702135578396*nuVtSqSum[3]*edgeSurf[10]+0.6210590034081186*nuVtSqSum[1]*edgeSurf[10]+0.6210590034081186*nuVtSqSum[2]*edgeSurf[7]+0.7071067811865474*nuVtSqSum[3]*edgeSurf[5]; 
  edgeSurf_incr[14] = 0.42163702135578396*nuVtSqSum[3]*edgeSurf[15]+0.6210590034081186*nuVtSqSum[1]*edgeSurf[15]+0.45175395145262565*nuVtSqSum[2]*edgeSurf[14]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[14]+0.6210590034081186*nuVtSqSum[3]*edgeSurf[12]+0.6324555320336759*nuVtSqSum[1]*edgeSurf[12]+0.7071067811865474*nuVtSqSum[2]*edgeSurf[9]; 
  edgeSurf_incr[15] = 0.421637021355784*nuVtSqSum[2]*edgeSurf[15]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[15]+0.42163702135578396*nuVtSqSum[3]*edgeSurf[14]+0.6210590034081186*nuVtSqSum[1]*edgeSurf[14]+0.6210590034081187*nuVtSqSum[2]*edgeSurf[12]+0.7071067811865475*nuVtSqSum[3]*edgeSurf[9]; 

  boundSurf_incr[0] = 0.7071067811865475*nuVtSqSum[3]*boundSurf[8]+0.7071067811865475*nuVtSqSum[2]*boundSurf[4]+0.7071067811865475*boundSurf[1]*nuVtSqSum[1]+0.7071067811865475*boundSurf[0]*nuVtSqSum[0]; 
  boundSurf_incr[1] = 0.6210590034081186*nuVtSqSum[2]*boundSurf[8]+0.6210590034081186*nuVtSqSum[3]*boundSurf[4]+0.6324555320336759*nuVtSqSum[1]*boundSurf[4]+0.6324555320336759*boundSurf[1]*nuVtSqSum[2]+0.7071067811865475*boundSurf[0]*nuVtSqSum[1]+0.7071067811865475*nuVtSqSum[0]*boundSurf[1]; 
  boundSurf_incr[2] = 0.7071067811865474*nuVtSqSum[3]*boundSurf[11]+0.7071067811865475*nuVtSqSum[2]*boundSurf[6]+0.7071067811865475*nuVtSqSum[1]*boundSurf[3]+0.7071067811865475*nuVtSqSum[0]*boundSurf[2]; 
  boundSurf_incr[3] = 0.6210590034081187*nuVtSqSum[2]*boundSurf[11]+0.6210590034081187*nuVtSqSum[3]*boundSurf[6]+0.632455532033676*nuVtSqSum[1]*boundSurf[6]+0.6324555320336759*nuVtSqSum[2]*boundSurf[3]+0.7071067811865475*nuVtSqSum[0]*boundSurf[3]+0.7071067811865475*nuVtSqSum[1]*boundSurf[2]; 
  boundSurf_incr[4] = 0.421637021355784*nuVtSqSum[3]*boundSurf[8]+0.6210590034081186*nuVtSqSum[1]*boundSurf[8]+0.45175395145262565*nuVtSqSum[2]*boundSurf[4]+0.7071067811865475*nuVtSqSum[0]*boundSurf[4]+0.6210590034081186*boundSurf[1]*nuVtSqSum[3]+0.7071067811865475*boundSurf[0]*nuVtSqSum[2]+0.6324555320336759*boundSurf[1]*nuVtSqSum[1]; 
  boundSurf_incr[5] = 0.7071067811865474*nuVtSqSum[3]*boundSurf[13]+0.7071067811865475*nuVtSqSum[2]*boundSurf[10]+0.7071067811865475*nuVtSqSum[1]*boundSurf[7]+0.7071067811865475*nuVtSqSum[0]*boundSurf[5]; 
  boundSurf_incr[6] = 0.4216370213557839*nuVtSqSum[3]*boundSurf[11]+0.6210590034081187*nuVtSqSum[1]*boundSurf[11]+0.45175395145262565*nuVtSqSum[2]*boundSurf[6]+0.7071067811865475*nuVtSqSum[0]*boundSurf[6]+0.6210590034081187*boundSurf[3]*nuVtSqSum[3]+0.632455532033676*nuVtSqSum[1]*boundSurf[3]+0.7071067811865475*boundSurf[2]*nuVtSqSum[2]; 
  boundSurf_incr[7] = 0.6210590034081186*nuVtSqSum[2]*boundSurf[13]+0.6210590034081188*nuVtSqSum[3]*boundSurf[10]+0.632455532033676*nuVtSqSum[1]*boundSurf[10]+0.6324555320336759*nuVtSqSum[2]*boundSurf[7]+0.7071067811865475*nuVtSqSum[0]*boundSurf[7]+0.7071067811865475*nuVtSqSum[1]*boundSurf[5]; 
  boundSurf_incr[8] = 0.421637021355784*nuVtSqSum[2]*boundSurf[8]+0.7071067811865475*nuVtSqSum[0]*boundSurf[8]+0.421637021355784*nuVtSqSum[3]*boundSurf[4]+0.6210590034081186*nuVtSqSum[1]*boundSurf[4]+0.7071067811865475*boundSurf[0]*nuVtSqSum[3]+0.6210590034081186*boundSurf[1]*nuVtSqSum[2]; 
  boundSurf_incr[9] = 0.7071067811865475*nuVtSqSum[3]*boundSurf[15]+0.7071067811865474*nuVtSqSum[2]*boundSurf[14]+0.7071067811865474*nuVtSqSum[1]*boundSurf[12]+0.7071067811865475*nuVtSqSum[0]*boundSurf[9]; 
  boundSurf_incr[10] = 0.42163702135578396*nuVtSqSum[3]*boundSurf[13]+0.6210590034081186*nuVtSqSum[1]*boundSurf[13]+0.45175395145262565*nuVtSqSum[2]*boundSurf[10]+0.7071067811865475*nuVtSqSum[0]*boundSurf[10]+0.6210590034081187*nuVtSqSum[3]*boundSurf[7]+0.632455532033676*nuVtSqSum[1]*boundSurf[7]+0.7071067811865475*nuVtSqSum[2]*boundSurf[5]; 
  boundSurf_incr[11] = 0.421637021355784*nuVtSqSum[2]*boundSurf[11]+0.7071067811865475*nuVtSqSum[0]*boundSurf[11]+0.4216370213557839*nuVtSqSum[3]*boundSurf[6]+0.6210590034081187*nuVtSqSum[1]*boundSurf[6]+0.7071067811865474*boundSurf[2]*nuVtSqSum[3]+0.6210590034081187*nuVtSqSum[2]*boundSurf[3]; 
  boundSurf_incr[12] = 0.6210590034081186*nuVtSqSum[2]*boundSurf[15]+0.6210590034081186*nuVtSqSum[3]*boundSurf[14]+0.6324555320336759*nuVtSqSum[1]*boundSurf[14]+0.6324555320336759*nuVtSqSum[2]*boundSurf[12]+0.7071067811865475*nuVtSqSum[0]*boundSurf[12]+0.7071067811865474*nuVtSqSum[1]*boundSurf[9]; 
  boundSurf_incr[13] = 0.421637021355784*nuVtSqSum[2]*boundSurf[13]+0.7071067811865475*nuVtSqSum[0]*boundSurf[13]+0.42163702135578396*nuVtSqSum[3]*boundSurf[10]+0.6210590034081186*nuVtSqSum[1]*boundSurf[10]+0.6210590034081186*nuVtSqSum[2]*boundSurf[7]+0.7071067811865474*nuVtSqSum[3]*boundSurf[5]; 
  boundSurf_incr[14] = 0.42163702135578396*nuVtSqSum[3]*boundSurf[15]+0.6210590034081186*nuVtSqSum[1]*boundSurf[15]+0.45175395145262565*nuVtSqSum[2]*boundSurf[14]+0.7071067811865475*nuVtSqSum[0]*boundSurf[14]+0.6210590034081186*nuVtSqSum[3]*boundSurf[12]+0.6324555320336759*nuVtSqSum[1]*boundSurf[12]+0.7071067811865474*nuVtSqSum[2]*boundSurf[9]; 
  boundSurf_incr[15] = 0.421637021355784*nuVtSqSum[2]*boundSurf[15]+0.7071067811865475*nuVtSqSum[0]*boundSurf[15]+0.42163702135578396*nuVtSqSum[3]*boundSurf[14]+0.6210590034081186*nuVtSqSum[1]*boundSurf[14]+0.6210590034081187*nuVtSqSum[2]*boundSurf[12]+0.7071067811865475*nuVtSqSum[3]*boundSurf[9]; 

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

  return 0.;

} 
