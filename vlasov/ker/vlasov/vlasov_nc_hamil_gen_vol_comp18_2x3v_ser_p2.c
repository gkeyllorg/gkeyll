#include <gkyl_nc_hamil_gen_vol_comps_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp18_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
   const double *jacob_pos, const double *jacob_vel,
   const double *alpha, double* GKYL_RESTRICT out) 
{ 
  const double dx10 = 2.0/dxv[0]; 
  const double dx11 = 2.0/dxv[1]; 
  const double dv10 = 2.0/dxv[2]; 
  const double dv11 = 2.0/dxv[3]; 
  const double dv12 = 2.0/dxv[4]; 
  const double jacob_cx_inv = 1.0/jacob_pos[0]; 
  const double jacob_cy_inv = 1.0/jacob_pos[3]; 
  const double jacob_vx_inv = 1.0/jacob_vel[0]; 
  const double jacob_vy_inv = 1.0/jacob_vel[3]; 
  const double jacob_vz_inv = 1.0/jacob_vel[6]; 
  
  if (dir == 0) { 
  } 

  if (dir == 1) { 
  } 

  if (dir == 2) { 
  out[0] += 0.6846531968814573*alpha[3]*dv10*jacob_vx_inv; 
  out[1] += 0.6846531968814573*alpha[7]*dv10*jacob_vx_inv; 
  out[2] += 0.6846531968814573*alpha[8]*dv10*jacob_vx_inv; 
  out[3] += (0.6123724356957944*alpha[18]+0.6846531968814573*alpha[0])*dv10*jacob_vx_inv; 
  out[4] += 0.6846531968814573*alpha[11]*dv10*jacob_vx_inv; 
  out[5] += 0.6846531968814573*alpha[14]*dv10*jacob_vx_inv; 
  out[6] += 0.6846531968814573*alpha[21]*dv10*jacob_vx_inv; 
  out[7] += (0.6123724356957944*alpha[35]+0.6846531968814573*alpha[1])*dv10*jacob_vx_inv; 
  out[8] += (0.6123724356957944*alpha[36]+0.6846531968814573*alpha[2])*dv10*jacob_vx_inv; 
  out[9] += 0.6846531968814573*alpha[23]*dv10*jacob_vx_inv; 
  out[10] += 0.6846531968814573*alpha[24]*dv10*jacob_vx_inv; 
  out[11] += (0.6123724356957944*alpha[39]+0.6846531968814573*alpha[4])*dv10*jacob_vx_inv; 
  out[12] += 0.6846531968814573*alpha[26]*dv10*jacob_vx_inv; 
  out[13] += 0.6846531968814573*alpha[27]*dv10*jacob_vx_inv; 
  out[14] += (0.6123724356957944*alpha[45]+0.6846531968814573*alpha[5])*dv10*jacob_vx_inv; 
  out[15] += 0.6846531968814573*alpha[30]*dv10*jacob_vx_inv; 
  out[16] += 0.6846531968814574*alpha[33]*dv10*jacob_vx_inv; 
  out[17] += 0.6846531968814574*alpha[34]*dv10*jacob_vx_inv; 
  out[18] += 0.6123724356957944*alpha[3]*dv10*jacob_vx_inv; 
  out[19] += 0.6846531968814574*alpha[42]*dv10*jacob_vx_inv; 
  out[20] += 0.6846531968814574*alpha[49]*dv10*jacob_vx_inv; 
  out[21] += (0.6123724356957944*alpha[58]+0.6846531968814573*alpha[6])*dv10*jacob_vx_inv; 
  out[22] += 0.6846531968814573*alpha[51]*dv10*jacob_vx_inv; 
  out[23] += (0.6123724356957944*alpha[63]+0.6846531968814573*alpha[9])*dv10*jacob_vx_inv; 
  out[24] += (0.6123724356957944*alpha[64]+0.6846531968814573*alpha[10])*dv10*jacob_vx_inv; 
  out[25] += 0.6846531968814573*alpha[52]*dv10*jacob_vx_inv; 
  out[26] += (0.6123724356957944*alpha[72]+0.6846531968814573*alpha[12])*dv10*jacob_vx_inv; 
  out[27] += (0.6123724356957944*alpha[73]+0.6846531968814573*alpha[13])*dv10*jacob_vx_inv; 
  out[28] += 0.6846531968814573*alpha[54]*dv10*jacob_vx_inv; 
  out[29] += 0.6846531968814573*alpha[55]*dv10*jacob_vx_inv; 
  out[30] += (0.6123724356957944*alpha[76]+0.6846531968814573*alpha[15])*dv10*jacob_vx_inv; 
  out[31] += 0.6846531968814574*alpha[56]*dv10*jacob_vx_inv; 
  out[32] += 0.6846531968814574*alpha[57]*dv10*jacob_vx_inv; 
  out[33] += 0.6846531968814574*alpha[16]*dv10*jacob_vx_inv; 
  out[34] += 0.6846531968814574*alpha[17]*dv10*jacob_vx_inv; 
  out[35] += 0.6123724356957944*alpha[7]*dv10*jacob_vx_inv; 
  out[36] += 0.6123724356957944*alpha[8]*dv10*jacob_vx_inv; 
  out[37] += 0.6846531968814574*alpha[61]*dv10*jacob_vx_inv; 
  out[38] += 0.6846531968814574*alpha[62]*dv10*jacob_vx_inv; 
  out[39] += 0.6123724356957944*alpha[11]*dv10*jacob_vx_inv; 
  out[40] += 0.6846531968814574*alpha[66]*dv10*jacob_vx_inv; 
  out[41] += 0.6846531968814574*alpha[67]*dv10*jacob_vx_inv; 
  out[42] += 0.6846531968814574*alpha[19]*dv10*jacob_vx_inv; 
  out[43] += 0.6846531968814574*alpha[70]*dv10*jacob_vx_inv; 
  out[44] += 0.6846531968814574*alpha[71]*dv10*jacob_vx_inv; 
  out[45] += 0.6123724356957944*alpha[14]*dv10*jacob_vx_inv; 
  out[46] += 0.6846531968814574*alpha[79]*dv10*jacob_vx_inv; 
  out[47] += 0.6846531968814574*alpha[81]*dv10*jacob_vx_inv; 
  out[48] += 0.6846531968814574*alpha[82]*dv10*jacob_vx_inv; 
  out[49] += 0.6846531968814574*alpha[20]*dv10*jacob_vx_inv; 
  out[50] += 0.6846531968814574*alpha[85]*dv10*jacob_vx_inv; 
  out[51] += (0.6123724356957944*alpha[89]+0.6846531968814573*alpha[22])*dv10*jacob_vx_inv; 
  out[52] += (0.6123724356957944*alpha[93]+0.6846531968814573*alpha[25])*dv10*jacob_vx_inv; 
  out[53] += 0.6846531968814573*alpha[86]*dv10*jacob_vx_inv; 
  out[54] += (0.6123724356957944*alpha[98]+0.6846531968814573*alpha[28])*dv10*jacob_vx_inv; 
  out[55] += (0.6123724356957944*alpha[99]+0.6846531968814573*alpha[29])*dv10*jacob_vx_inv; 
  out[56] += 0.6846531968814574*alpha[31]*dv10*jacob_vx_inv; 
  out[57] += 0.6846531968814574*alpha[32]*dv10*jacob_vx_inv; 
  out[58] += 0.6123724356957944*alpha[21]*dv10*jacob_vx_inv; 
  out[59] += 0.6846531968814574*alpha[87]*dv10*jacob_vx_inv; 
  out[60] += 0.6846531968814574*alpha[88]*dv10*jacob_vx_inv; 
  out[61] += 0.6846531968814574*alpha[37]*dv10*jacob_vx_inv; 
  out[62] += 0.6846531968814574*alpha[38]*dv10*jacob_vx_inv; 
  out[63] += 0.6123724356957944*alpha[23]*dv10*jacob_vx_inv; 
  out[64] += 0.6123724356957944*alpha[24]*dv10*jacob_vx_inv; 
  out[65] += 0.6846531968814574*alpha[90]*dv10*jacob_vx_inv; 
  out[66] += 0.6846531968814574*alpha[40]*dv10*jacob_vx_inv; 
  out[67] += 0.6846531968814574*alpha[41]*dv10*jacob_vx_inv; 
  out[68] += 0.6846531968814574*alpha[91]*dv10*jacob_vx_inv; 
  out[69] += 0.6846531968814574*alpha[92]*dv10*jacob_vx_inv; 
  out[70] += 0.6846531968814574*alpha[43]*dv10*jacob_vx_inv; 
  out[71] += 0.6846531968814574*alpha[44]*dv10*jacob_vx_inv; 
  out[72] += 0.6123724356957944*alpha[26]*dv10*jacob_vx_inv; 
  out[73] += 0.6123724356957944*alpha[27]*dv10*jacob_vx_inv; 
  out[74] += 0.6846531968814574*alpha[96]*dv10*jacob_vx_inv; 
  out[75] += 0.6846531968814574*alpha[97]*dv10*jacob_vx_inv; 
  out[76] += 0.6123724356957944*alpha[30]*dv10*jacob_vx_inv; 
  out[77] += 0.6846531968814574*alpha[101]*dv10*jacob_vx_inv; 
  out[78] += 0.6846531968814574*alpha[102]*dv10*jacob_vx_inv; 
  out[79] += 0.6846531968814574*alpha[46]*dv10*jacob_vx_inv; 
  out[80] += 0.6846531968814574*alpha[103]*dv10*jacob_vx_inv; 
  out[81] += 0.6846531968814574*alpha[47]*dv10*jacob_vx_inv; 
  out[82] += 0.6846531968814574*alpha[48]*dv10*jacob_vx_inv; 
  out[83] += 0.6846531968814574*alpha[105]*dv10*jacob_vx_inv; 
  out[84] += 0.6846531968814574*alpha[106]*dv10*jacob_vx_inv; 
  out[85] += 0.6846531968814574*alpha[50]*dv10*jacob_vx_inv; 
  out[86] += (0.6123724356957944*alpha[109]+0.6846531968814573*alpha[53])*dv10*jacob_vx_inv; 
  out[87] += 0.6846531968814574*alpha[59]*dv10*jacob_vx_inv; 
  out[88] += 0.6846531968814574*alpha[60]*dv10*jacob_vx_inv; 
  out[89] += 0.6123724356957944*alpha[51]*dv10*jacob_vx_inv; 
  out[90] += 0.6846531968814574*alpha[65]*dv10*jacob_vx_inv; 
  out[91] += 0.6846531968814574*alpha[68]*dv10*jacob_vx_inv; 
  out[92] += 0.6846531968814574*alpha[69]*dv10*jacob_vx_inv; 
  out[93] += 0.6123724356957944*alpha[52]*dv10*jacob_vx_inv; 
  out[94] += 0.6846531968814574*alpha[107]*dv10*jacob_vx_inv; 
  out[95] += 0.6846531968814574*alpha[108]*dv10*jacob_vx_inv; 
  out[96] += 0.6846531968814574*alpha[74]*dv10*jacob_vx_inv; 
  out[97] += 0.6846531968814574*alpha[75]*dv10*jacob_vx_inv; 
  out[98] += 0.6123724356957944*alpha[54]*dv10*jacob_vx_inv; 
  out[99] += 0.6123724356957944*alpha[55]*dv10*jacob_vx_inv; 
  out[100] += 0.6846531968814574*alpha[110]*dv10*jacob_vx_inv; 
  out[101] += 0.6846531968814574*alpha[77]*dv10*jacob_vx_inv; 
  out[102] += 0.6846531968814574*alpha[78]*dv10*jacob_vx_inv; 
  out[103] += 0.6846531968814574*alpha[80]*dv10*jacob_vx_inv; 
  out[104] += 0.6846531968814574*alpha[111]*dv10*jacob_vx_inv; 
  out[105] += 0.6846531968814574*alpha[83]*dv10*jacob_vx_inv; 
  out[106] += 0.6846531968814574*alpha[84]*dv10*jacob_vx_inv; 
  out[107] += 0.6846531968814574*alpha[94]*dv10*jacob_vx_inv; 
  out[108] += 0.6846531968814574*alpha[95]*dv10*jacob_vx_inv; 
  out[109] += 0.6123724356957944*alpha[86]*dv10*jacob_vx_inv; 
  out[110] += 0.6846531968814574*alpha[100]*dv10*jacob_vx_inv; 
  out[111] += 0.6846531968814574*alpha[104]*dv10*jacob_vx_inv; 
  } 

  if (dir == 3) { 
  } 

  if (dir == 4) { 
  } 

} 
