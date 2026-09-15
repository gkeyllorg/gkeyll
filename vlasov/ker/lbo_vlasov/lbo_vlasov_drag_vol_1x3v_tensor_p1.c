#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_drag_vol_1x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[4]: Cell-center coordinates. 
  // dxv[4]: Cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // f: Input distribution function.
  // out: Incremented output 
  const double *nuUSum = nuPrimMomsSum;

  const double rdvx2 = 2.0/dxv[1]; 
  const double rdvy2 = 2.0/dxv[2]; 
  const double rdvz2 = 2.0/dxv[3]; 

  double alphaDrag[162]; 
  // Expand rdv2*(nu*vx-nuUSumx) in phase basis.
  alphaDrag[0] = (2.8284271247461907*nuUSum[0]-2.8284271247461907*nuSum[0]*w[1])*rdvx2; 
  alphaDrag[1] = (2.8284271247461907*nuUSum[1]-2.8284271247461907*nuSum[1]*w[1])*rdvx2; 
  alphaDrag[2] = -(0.8164965809277261*nuSum[0]*dxv[1]*rdvx2); 
  alphaDrag[5] = -(0.8164965809277261*dxv[1]*nuSum[1]*rdvx2); 

  // Expand rdv2*(nu*vy-nuUSumy) in phase basis.
  alphaDrag[54] = (2.8284271247461907*nuUSum[2]-2.8284271247461907*nuSum[0]*w[2])*rdvy2; 
  alphaDrag[55] = (2.8284271247461907*nuUSum[3]-2.8284271247461907*nuSum[1]*w[2])*rdvy2; 
  alphaDrag[57] = -(0.8164965809277261*nuSum[0]*dxv[2]*rdvy2); 
  alphaDrag[60] = -(0.8164965809277261*nuSum[1]*dxv[2]*rdvy2); 

  // Expand rdv2*(nu*vz-nuUSumz) in phase basis.
  alphaDrag[108] = (2.8284271247461907*nuUSum[4]-2.8284271247461907*nuSum[0]*w[3])*rdvz2; 
  alphaDrag[109] = (2.8284271247461907*nuUSum[5]-2.8284271247461907*nuSum[1]*w[3])*rdvz2; 
  alphaDrag[112] = -(0.8164965809277261*nuSum[0]*dxv[3]*rdvz2); 
  alphaDrag[116] = -(0.8164965809277261*nuSum[1]*dxv[3]*rdvz2); 

  out[2] += 0.4330127018922193*(alphaDrag[5]*f[5]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.4330127018922193*(f[6]*alphaDrag[60]+f[3]*alphaDrag[57]+f[1]*alphaDrag[55]+f[0]*alphaDrag[54]); 
  out[4] += 0.4330127018922193*(f[8]*alphaDrag[116]+f[4]*alphaDrag[112]+f[1]*alphaDrag[109]+f[0]*alphaDrag[108]); 
  out[5] += 0.4330127018922193*(alphaDrag[2]*f[5]+f[2]*alphaDrag[5]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[6] += 0.4330127018922193*(f[3]*alphaDrag[60]+f[6]*alphaDrag[57]+f[0]*alphaDrag[55]+f[1]*alphaDrag[54]); 
  out[7] += 0.4330127018922193*(f[14]*alphaDrag[60]+f[7]*alphaDrag[57]+f[5]*alphaDrag[55]+f[2]*alphaDrag[54]+alphaDrag[5]*f[14]+alphaDrag[2]*f[7]+alphaDrag[1]*f[6]+alphaDrag[0]*f[3]); 
  out[8] += 0.4330127018922193*(f[4]*alphaDrag[116]+f[8]*alphaDrag[112]+f[0]*alphaDrag[109]+f[1]*alphaDrag[108]); 
  out[9] += 0.4330127018922193*(f[15]*alphaDrag[116]+f[9]*alphaDrag[112]+f[5]*alphaDrag[109]+f[2]*alphaDrag[108]+alphaDrag[5]*f[15]+alphaDrag[2]*f[9]+alphaDrag[1]*f[8]+alphaDrag[0]*f[4]); 
  out[10] += 0.4330127018922193*(f[16]*alphaDrag[116]+f[10]*alphaDrag[112]+f[6]*alphaDrag[109]+f[3]*alphaDrag[108]+f[16]*alphaDrag[60]+f[10]*alphaDrag[57]+f[8]*alphaDrag[55]+f[4]*alphaDrag[54]); 
  out[11] += 0.8660254037844386*(alphaDrag[5]*f[18]+alphaDrag[2]*f[11])+0.9682458365518543*(alphaDrag[1]*f[5]+f[1]*alphaDrag[5]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[12] += (0.8660254037844386*f[20]+0.9682458365518543*f[1])*alphaDrag[60]+0.8660254037844386*f[12]*alphaDrag[57]+0.9682458365518543*(f[0]*alphaDrag[57]+f[6]*alphaDrag[55]+f[3]*alphaDrag[54]); 
  out[13] += (0.8660254037844386*f[24]+0.9682458365518543*f[1])*alphaDrag[116]+0.8660254037844386*f[13]*alphaDrag[112]+0.9682458365518543*(f[0]*alphaDrag[112]+f[8]*alphaDrag[109]+f[4]*alphaDrag[108]); 
  out[14] += 0.4330127018922193*(f[7]*alphaDrag[60]+f[14]*alphaDrag[57]+f[2]*alphaDrag[55]+f[5]*alphaDrag[54]+alphaDrag[2]*f[14]+alphaDrag[5]*f[7]+alphaDrag[0]*f[6]+alphaDrag[1]*f[3]); 
  out[15] += 0.4330127018922193*(f[9]*alphaDrag[116]+f[15]*alphaDrag[112]+f[2]*alphaDrag[109]+f[5]*alphaDrag[108]+alphaDrag[2]*f[15]+alphaDrag[5]*f[9]+alphaDrag[0]*f[8]+alphaDrag[1]*f[4]); 
  out[16] += 0.4330127018922193*(f[10]*alphaDrag[116]+f[16]*alphaDrag[112]+f[3]*alphaDrag[109]+f[6]*alphaDrag[108]+f[10]*alphaDrag[60]+f[16]*alphaDrag[57]+f[4]*alphaDrag[55]+f[8]*alphaDrag[54]); 
  out[17] += 0.4330127018922193*(f[27]*alphaDrag[116]+f[17]*alphaDrag[112]+f[14]*alphaDrag[109]+f[7]*alphaDrag[108]+f[27]*alphaDrag[60]+f[17]*alphaDrag[57]+f[15]*alphaDrag[55]+f[9]*alphaDrag[54]+alphaDrag[5]*f[27]+alphaDrag[2]*f[17]+alphaDrag[1]*f[16]+alphaDrag[0]*f[10]); 
  out[18] += 0.8660254037844386*(alphaDrag[2]*f[18]+alphaDrag[5]*f[11])+0.9682458365518543*(alphaDrag[0]*f[5]+f[0]*alphaDrag[5]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[19] += 0.4330127018922193*(f[28]*alphaDrag[60]+f[19]*alphaDrag[57]+f[18]*alphaDrag[55]+f[11]*alphaDrag[54])+0.8660254037844386*(alphaDrag[5]*f[28]+alphaDrag[2]*f[19])+0.9682458365518543*(alphaDrag[1]*f[14]+alphaDrag[0]*f[7]+alphaDrag[5]*f[6]+alphaDrag[2]*f[3]); 
  out[20] += (0.8660254037844386*f[12]+0.9682458365518543*f[0])*alphaDrag[60]+0.8660254037844386*f[20]*alphaDrag[57]+0.9682458365518543*(f[1]*alphaDrag[57]+f[3]*alphaDrag[55]+f[6]*alphaDrag[54]); 
  out[21] += (0.8660254037844386*f[29]+0.9682458365518543*f[5])*alphaDrag[60]+0.8660254037844386*f[21]*alphaDrag[57]+0.9682458365518543*(f[2]*alphaDrag[57]+f[14]*alphaDrag[55]+f[7]*alphaDrag[54])+0.4330127018922193*(alphaDrag[5]*f[29]+alphaDrag[2]*f[21]+alphaDrag[1]*f[20]+alphaDrag[0]*f[12]); 
  out[22] += 0.4330127018922193*(f[30]*alphaDrag[116]+f[22]*alphaDrag[112]+f[18]*alphaDrag[109]+f[11]*alphaDrag[108])+0.8660254037844386*(alphaDrag[5]*f[30]+alphaDrag[2]*f[22])+0.9682458365518543*(alphaDrag[1]*f[15]+alphaDrag[0]*f[9]+alphaDrag[5]*f[8]+alphaDrag[2]*f[4]); 
  out[23] += 0.4330127018922193*(f[32]*alphaDrag[116]+f[23]*alphaDrag[112]+f[20]*alphaDrag[109]+f[12]*alphaDrag[108])+(0.8660254037844386*f[32]+0.9682458365518543*f[8])*alphaDrag[60]+0.8660254037844386*f[23]*alphaDrag[57]+0.9682458365518543*(f[4]*alphaDrag[57]+f[16]*alphaDrag[55]+f[10]*alphaDrag[54]); 
  out[24] += (0.8660254037844386*f[13]+0.9682458365518543*f[0])*alphaDrag[116]+0.8660254037844386*f[24]*alphaDrag[112]+0.9682458365518543*(f[1]*alphaDrag[112]+f[4]*alphaDrag[109]+f[8]*alphaDrag[108]); 
  out[25] += (0.8660254037844386*f[34]+0.9682458365518543*f[5])*alphaDrag[116]+0.8660254037844386*f[25]*alphaDrag[112]+0.9682458365518543*(f[2]*alphaDrag[112]+f[15]*alphaDrag[109]+f[9]*alphaDrag[108])+0.4330127018922193*(alphaDrag[5]*f[34]+alphaDrag[2]*f[25]+alphaDrag[1]*f[24]+alphaDrag[0]*f[13]); 
  out[26] += (0.8660254037844386*f[35]+0.9682458365518543*f[6])*alphaDrag[116]+0.8660254037844386*f[26]*alphaDrag[112]+0.9682458365518543*(f[3]*alphaDrag[112]+f[16]*alphaDrag[109]+f[10]*alphaDrag[108])+0.4330127018922193*(f[35]*alphaDrag[60]+f[26]*alphaDrag[57]+f[24]*alphaDrag[55]+f[13]*alphaDrag[54]); 
  out[27] += 0.4330127018922193*(f[17]*alphaDrag[116]+f[27]*alphaDrag[112]+f[7]*alphaDrag[109]+f[14]*alphaDrag[108]+f[17]*alphaDrag[60]+f[27]*alphaDrag[57]+f[9]*alphaDrag[55]+f[15]*alphaDrag[54]+alphaDrag[2]*f[27]+alphaDrag[5]*f[17]+alphaDrag[0]*f[16]+alphaDrag[1]*f[10]); 
  out[28] += 0.4330127018922193*(f[19]*alphaDrag[60]+f[28]*alphaDrag[57]+f[11]*alphaDrag[55]+f[18]*alphaDrag[54])+0.8660254037844386*(alphaDrag[2]*f[28]+alphaDrag[5]*f[19])+0.9682458365518543*(alphaDrag[0]*f[14]+alphaDrag[1]*f[7]+alphaDrag[2]*f[6]+f[3]*alphaDrag[5]); 
  out[29] += (0.8660254037844386*f[21]+0.9682458365518543*f[2])*alphaDrag[60]+0.8660254037844386*f[29]*alphaDrag[57]+0.9682458365518543*(f[5]*alphaDrag[57]+f[7]*alphaDrag[55]+f[14]*alphaDrag[54])+0.4330127018922193*(alphaDrag[2]*f[29]+alphaDrag[5]*f[21]+alphaDrag[0]*f[20]+alphaDrag[1]*f[12]); 
  out[30] += 0.4330127018922193*(f[22]*alphaDrag[116]+f[30]*alphaDrag[112]+f[11]*alphaDrag[109]+f[18]*alphaDrag[108])+0.8660254037844386*(alphaDrag[2]*f[30]+alphaDrag[5]*f[22])+0.9682458365518543*(alphaDrag[0]*f[15]+alphaDrag[1]*f[9]+alphaDrag[2]*f[8]+f[4]*alphaDrag[5]); 
  out[31] += 0.4330127018922193*(f[40]*alphaDrag[116]+f[31]*alphaDrag[112]+f[28]*alphaDrag[109]+f[19]*alphaDrag[108]+f[40]*alphaDrag[60]+f[31]*alphaDrag[57]+f[30]*alphaDrag[55]+f[22]*alphaDrag[54])+0.8660254037844386*(alphaDrag[5]*f[40]+alphaDrag[2]*f[31])+0.9682458365518543*(alphaDrag[1]*f[27]+alphaDrag[0]*f[17]+alphaDrag[5]*f[16]+alphaDrag[2]*f[10]); 
  out[32] += 0.4330127018922193*(f[23]*alphaDrag[116]+f[32]*alphaDrag[112]+f[12]*alphaDrag[109]+f[20]*alphaDrag[108])+(0.8660254037844386*f[23]+0.9682458365518543*f[4])*alphaDrag[60]+0.8660254037844386*f[32]*alphaDrag[57]+0.9682458365518543*(f[8]*alphaDrag[57]+f[10]*alphaDrag[55]+f[16]*alphaDrag[54]); 
  out[33] += 0.4330127018922193*(f[41]*alphaDrag[116]+f[33]*alphaDrag[112]+f[29]*alphaDrag[109]+f[21]*alphaDrag[108])+(0.8660254037844386*f[41]+0.9682458365518543*f[15])*alphaDrag[60]+0.8660254037844386*f[33]*alphaDrag[57]+0.9682458365518543*(f[9]*alphaDrag[57]+f[27]*alphaDrag[55]+f[17]*alphaDrag[54])+0.4330127018922193*(alphaDrag[5]*f[41]+alphaDrag[2]*f[33]+alphaDrag[1]*f[32]+alphaDrag[0]*f[23]); 
  out[34] += (0.8660254037844386*f[25]+0.9682458365518543*f[2])*alphaDrag[116]+0.8660254037844386*f[34]*alphaDrag[112]+0.9682458365518543*(f[5]*alphaDrag[112]+f[9]*alphaDrag[109]+f[15]*alphaDrag[108])+0.4330127018922193*(alphaDrag[2]*f[34]+alphaDrag[5]*f[25]+alphaDrag[0]*f[24]+alphaDrag[1]*f[13]); 
  out[35] += (0.8660254037844386*f[26]+0.9682458365518543*f[3])*alphaDrag[116]+0.8660254037844386*f[35]*alphaDrag[112]+0.9682458365518543*(f[6]*alphaDrag[112]+f[10]*alphaDrag[109]+f[16]*alphaDrag[108])+0.4330127018922193*(f[26]*alphaDrag[60]+f[35]*alphaDrag[57]+f[13]*alphaDrag[55]+f[24]*alphaDrag[54]); 
  out[36] += (0.8660254037844386*f[42]+0.9682458365518543*f[14])*alphaDrag[116]+0.8660254037844386*f[36]*alphaDrag[112]+0.9682458365518543*(f[7]*alphaDrag[112]+f[27]*alphaDrag[109]+f[17]*alphaDrag[108])+0.4330127018922193*(f[42]*alphaDrag[60]+f[36]*alphaDrag[57]+f[34]*alphaDrag[55]+f[25]*alphaDrag[54]+alphaDrag[5]*f[42]+alphaDrag[2]*f[36]+alphaDrag[1]*f[35]+alphaDrag[0]*f[26]); 
  out[37] += (0.8660254037844386*f[43]+0.9682458365518543*f[18])*alphaDrag[60]+0.8660254037844386*f[37]*alphaDrag[57]+0.9682458365518543*(f[11]*alphaDrag[57]+f[28]*alphaDrag[55]+f[19]*alphaDrag[54])+0.8660254037844386*(alphaDrag[5]*f[43]+alphaDrag[2]*f[37])+0.9682458365518543*(alphaDrag[1]*f[29]+alphaDrag[0]*f[21]+alphaDrag[5]*f[20]+alphaDrag[2]*f[12]); 
  out[38] += (0.8660254037844386*f[45]+0.9682458365518543*f[18])*alphaDrag[116]+0.8660254037844386*f[38]*alphaDrag[112]+0.9682458365518543*(f[11]*alphaDrag[112]+f[30]*alphaDrag[109]+f[22]*alphaDrag[108])+0.8660254037844386*(alphaDrag[5]*f[45]+alphaDrag[2]*f[38])+0.9682458365518543*(alphaDrag[1]*f[34]+alphaDrag[0]*f[25]+alphaDrag[5]*f[24]+alphaDrag[2]*f[13]); 
  out[39] += (0.8660254037844386*f[47]+0.9682458365518543*f[20])*alphaDrag[116]+0.8660254037844386*f[39]*alphaDrag[112]+0.9682458365518543*(f[12]*alphaDrag[112]+f[32]*alphaDrag[109]+f[23]*alphaDrag[108])+(0.8660254037844386*f[47]+0.9682458365518543*f[24])*alphaDrag[60]+0.8660254037844386*f[39]*alphaDrag[57]+0.9682458365518543*(f[13]*alphaDrag[57]+f[35]*alphaDrag[55]+f[26]*alphaDrag[54]); 
  out[40] += 0.4330127018922193*(f[31]*alphaDrag[116]+f[40]*alphaDrag[112]+f[19]*alphaDrag[109]+f[28]*alphaDrag[108]+f[31]*alphaDrag[60]+f[40]*alphaDrag[57]+f[22]*alphaDrag[55]+f[30]*alphaDrag[54])+0.8660254037844386*(alphaDrag[2]*f[40]+alphaDrag[5]*f[31])+0.9682458365518543*(alphaDrag[0]*f[27]+alphaDrag[1]*f[17]+alphaDrag[2]*f[16]+alphaDrag[5]*f[10]); 
  out[41] += 0.4330127018922193*(f[33]*alphaDrag[116]+f[41]*alphaDrag[112]+f[21]*alphaDrag[109]+f[29]*alphaDrag[108])+(0.8660254037844386*f[33]+0.9682458365518543*f[9])*alphaDrag[60]+0.8660254037844386*f[41]*alphaDrag[57]+0.9682458365518543*(f[15]*alphaDrag[57]+f[17]*alphaDrag[55]+f[27]*alphaDrag[54])+0.4330127018922193*(alphaDrag[2]*f[41]+alphaDrag[5]*f[33]+alphaDrag[0]*f[32]+alphaDrag[1]*f[23]); 
  out[42] += (0.8660254037844386*f[36]+0.9682458365518543*f[7])*alphaDrag[116]+0.8660254037844386*f[42]*alphaDrag[112]+0.9682458365518543*(f[14]*alphaDrag[112]+f[17]*alphaDrag[109]+f[27]*alphaDrag[108])+0.4330127018922193*(f[36]*alphaDrag[60]+f[42]*alphaDrag[57]+f[25]*alphaDrag[55]+f[34]*alphaDrag[54]+alphaDrag[2]*f[42]+alphaDrag[5]*f[36]+alphaDrag[0]*f[35]+alphaDrag[1]*f[26]); 
  out[43] += (0.8660254037844386*f[37]+0.9682458365518543*f[11])*alphaDrag[60]+0.8660254037844386*f[43]*alphaDrag[57]+0.9682458365518543*(f[18]*alphaDrag[57]+f[19]*alphaDrag[55]+f[28]*alphaDrag[54])+0.8660254037844386*(alphaDrag[2]*f[43]+alphaDrag[5]*f[37])+0.9682458365518543*(alphaDrag[0]*f[29]+alphaDrag[1]*f[21]+alphaDrag[2]*f[20]+alphaDrag[5]*f[12]); 
  out[44] += 0.4330127018922193*(f[49]*alphaDrag[116]+f[44]*alphaDrag[112]+f[43]*alphaDrag[109]+f[37]*alphaDrag[108])+(0.8660254037844386*f[49]+0.9682458365518543*f[30])*alphaDrag[60]+0.8660254037844386*f[44]*alphaDrag[57]+0.9682458365518543*(f[22]*alphaDrag[57]+f[40]*alphaDrag[55]+f[31]*alphaDrag[54])+0.8660254037844386*(alphaDrag[5]*f[49]+alphaDrag[2]*f[44])+0.9682458365518543*(alphaDrag[1]*f[41]+alphaDrag[0]*f[33]+alphaDrag[5]*f[32]+alphaDrag[2]*f[23]); 
  out[45] += (0.8660254037844386*f[38]+0.9682458365518543*f[11])*alphaDrag[116]+0.8660254037844386*f[45]*alphaDrag[112]+0.9682458365518543*(f[18]*alphaDrag[112]+f[22]*alphaDrag[109]+f[30]*alphaDrag[108])+0.8660254037844386*(alphaDrag[2]*f[45]+alphaDrag[5]*f[38])+0.9682458365518543*(alphaDrag[0]*f[34]+alphaDrag[1]*f[25]+alphaDrag[2]*f[24]+alphaDrag[5]*f[13]); 
  out[46] += (0.8660254037844386*f[50]+0.9682458365518543*f[28])*alphaDrag[116]+0.8660254037844386*f[46]*alphaDrag[112]+0.9682458365518543*(f[19]*alphaDrag[112]+f[40]*alphaDrag[109]+f[31]*alphaDrag[108])+0.4330127018922193*(f[50]*alphaDrag[60]+f[46]*alphaDrag[57]+f[45]*alphaDrag[55]+f[38]*alphaDrag[54])+0.8660254037844386*(alphaDrag[5]*f[50]+alphaDrag[2]*f[46])+0.9682458365518543*(alphaDrag[1]*f[42]+alphaDrag[0]*f[36]+alphaDrag[5]*f[35]+alphaDrag[2]*f[26]); 
  out[47] += (0.8660254037844386*f[39]+0.9682458365518543*f[12])*alphaDrag[116]+0.8660254037844386*f[47]*alphaDrag[112]+0.9682458365518543*(f[20]*alphaDrag[112]+f[23]*alphaDrag[109]+f[32]*alphaDrag[108])+(0.8660254037844386*f[39]+0.9682458365518543*f[13])*alphaDrag[60]+0.8660254037844386*f[47]*alphaDrag[57]+0.9682458365518543*(f[24]*alphaDrag[57]+f[26]*alphaDrag[55]+f[35]*alphaDrag[54]); 
  out[48] += (0.8660254037844386*f[51]+0.9682458365518543*f[29])*alphaDrag[116]+0.8660254037844386*f[48]*alphaDrag[112]+0.9682458365518543*(f[21]*alphaDrag[112]+f[41]*alphaDrag[109]+f[33]*alphaDrag[108])+(0.8660254037844386*f[51]+0.9682458365518543*f[34])*alphaDrag[60]+0.8660254037844386*f[48]*alphaDrag[57]+0.9682458365518543*(f[25]*alphaDrag[57]+f[42]*alphaDrag[55]+f[36]*alphaDrag[54])+0.4330127018922193*(alphaDrag[5]*f[51]+alphaDrag[2]*f[48]+alphaDrag[1]*f[47]+alphaDrag[0]*f[39]); 
  out[49] += 0.4330127018922193*(f[44]*alphaDrag[116]+f[49]*alphaDrag[112]+f[37]*alphaDrag[109]+f[43]*alphaDrag[108])+(0.8660254037844386*f[44]+0.9682458365518543*f[22])*alphaDrag[60]+0.8660254037844386*f[49]*alphaDrag[57]+0.9682458365518543*(f[30]*alphaDrag[57]+f[31]*alphaDrag[55]+f[40]*alphaDrag[54])+0.8660254037844386*(alphaDrag[2]*f[49]+alphaDrag[5]*f[44])+0.9682458365518543*(alphaDrag[0]*f[41]+alphaDrag[1]*f[33]+alphaDrag[2]*f[32]+alphaDrag[5]*f[23]); 
  out[50] += (0.8660254037844386*f[46]+0.9682458365518543*f[19])*alphaDrag[116]+0.8660254037844386*f[50]*alphaDrag[112]+0.9682458365518543*(f[28]*alphaDrag[112]+f[31]*alphaDrag[109]+f[40]*alphaDrag[108])+0.4330127018922193*(f[46]*alphaDrag[60]+f[50]*alphaDrag[57]+f[38]*alphaDrag[55]+f[45]*alphaDrag[54])+0.8660254037844386*(alphaDrag[2]*f[50]+alphaDrag[5]*f[46])+0.9682458365518543*(alphaDrag[0]*f[42]+alphaDrag[1]*f[36]+alphaDrag[2]*f[35]+alphaDrag[5]*f[26]); 
  out[51] += (0.8660254037844386*f[48]+0.9682458365518543*f[21])*alphaDrag[116]+0.8660254037844386*f[51]*alphaDrag[112]+0.9682458365518543*(f[29]*alphaDrag[112]+f[33]*alphaDrag[109]+f[41]*alphaDrag[108])+(0.8660254037844386*f[48]+0.9682458365518543*f[25])*alphaDrag[60]+0.8660254037844386*f[51]*alphaDrag[57]+0.9682458365518543*(f[34]*alphaDrag[57]+f[36]*alphaDrag[55]+f[42]*alphaDrag[54])+0.4330127018922193*(alphaDrag[2]*f[51]+alphaDrag[5]*f[48]+alphaDrag[0]*f[47]+alphaDrag[1]*f[39]); 
  out[52] += (0.8660254037844386*f[53]+0.9682458365518543*f[43])*alphaDrag[116]+0.8660254037844386*f[52]*alphaDrag[112]+0.9682458365518543*(f[37]*alphaDrag[112]+f[49]*alphaDrag[109]+f[44]*alphaDrag[108])+(0.8660254037844386*f[53]+0.9682458365518543*f[45])*alphaDrag[60]+0.8660254037844386*f[52]*alphaDrag[57]+0.9682458365518543*(f[38]*alphaDrag[57]+f[50]*alphaDrag[55]+f[46]*alphaDrag[54])+0.8660254037844386*(alphaDrag[5]*f[53]+alphaDrag[2]*f[52])+0.9682458365518543*(alphaDrag[1]*f[51]+alphaDrag[0]*f[48]+alphaDrag[5]*f[47]+alphaDrag[2]*f[39]); 
  out[53] += (0.8660254037844386*f[52]+0.9682458365518543*f[37])*alphaDrag[116]+0.8660254037844386*f[53]*alphaDrag[112]+0.9682458365518543*(f[43]*alphaDrag[112]+f[44]*alphaDrag[109]+f[49]*alphaDrag[108])+(0.8660254037844386*f[52]+0.9682458365518543*f[38])*alphaDrag[60]+0.8660254037844386*f[53]*alphaDrag[57]+0.9682458365518543*(f[45]*alphaDrag[57]+f[46]*alphaDrag[55]+f[50]*alphaDrag[54])+0.8660254037844386*(alphaDrag[2]*f[53]+alphaDrag[5]*f[52])+0.9682458365518543*(alphaDrag[0]*f[51]+alphaDrag[1]*f[48]+alphaDrag[2]*f[47]+alphaDrag[5]*f[39]); 

  return fabs(0.625*alphaDrag[0])+fabs(0.625*alphaDrag[54])+fabs(0.625*alphaDrag[108]); 

} 
