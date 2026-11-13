#include <gkyl_nc_hamil_gen_vol_comps_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp45_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
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
  out[0] += 0.6846531968814574*alpha[14]*dv10*jacob_vx_inv; 
  out[1] += 0.6846531968814574*alpha[26]*dv10*jacob_vx_inv; 
  out[2] += 0.6846531968814574*alpha[27]*dv10*jacob_vx_inv; 
  out[3] += (0.6123724356957944*alpha[45]+0.6846531968814574*alpha[5])*dv10*jacob_vx_inv; 
  out[4] += 0.6846531968814574*alpha[30]*dv10*jacob_vx_inv; 
  out[5] += (0.6123724356957944*alpha[49]+0.6846531968814574*alpha[3])*dv10*jacob_vx_inv; 
  out[6] += 0.6846531968814574*alpha[52]*dv10*jacob_vx_inv; 
  out[7] += (0.6123724356957944*alpha[72]+0.6846531968814574*alpha[12])*dv10*jacob_vx_inv; 
  out[8] += (0.6123724356957944*alpha[73]+0.6846531968814574*alpha[13])*dv10*jacob_vx_inv; 
  out[9] += 0.6846531968814574*alpha[54]*dv10*jacob_vx_inv; 
  out[10] += 0.6846531968814574*alpha[55]*dv10*jacob_vx_inv; 
  out[11] += (0.6123724356957944*alpha[76]+0.6846531968814574*alpha[15])*dv10*jacob_vx_inv; 
  out[12] += (0.6123724356957944*alpha[81]+0.6846531968814574*alpha[7])*dv10*jacob_vx_inv; 
  out[13] += (0.6123724356957944*alpha[82]+0.6846531968814574*alpha[8])*dv10*jacob_vx_inv; 
  out[14] += (0.6123724356957944*alpha[20]+0.6123724356957944*alpha[18]+0.6846531968814574*alpha[0])*dv10*jacob_vx_inv; 
  out[15] += (0.6123724356957944*alpha[85]+0.6846531968814574*alpha[11])*dv10*jacob_vx_inv; 
  out[16] += 0.6846531968814574*alpha[70]*dv10*jacob_vx_inv; 
  out[17] += 0.6846531968814574*alpha[71]*dv10*jacob_vx_inv; 
  out[18] += 0.6123724356957944*alpha[14]*dv10*jacob_vx_inv; 
  out[19] += 0.6846531968814574*alpha[79]*dv10*jacob_vx_inv; 
  out[20] += 0.6123724356957944*alpha[14]*dv10*jacob_vx_inv; 
  out[21] += (0.6123724356957944*alpha[93]+0.6846531968814574*alpha[25])*dv10*jacob_vx_inv; 
  out[22] += 0.6846531968814574*alpha[86]*dv10*jacob_vx_inv; 
  out[23] += (0.6123724356957944*alpha[98]+0.6846531968814574*alpha[28])*dv10*jacob_vx_inv; 
  out[24] += (0.6123724356957944*alpha[99]+0.6846531968814574*alpha[29])*dv10*jacob_vx_inv; 
  out[25] += (0.6123724356957944*alpha[103]+0.6846531968814574*alpha[21])*dv10*jacob_vx_inv; 
  out[26] += (0.6123724356957944*alpha[47]+0.6123724356957944*alpha[35]+0.6846531968814574*alpha[1])*dv10*jacob_vx_inv; 
  out[27] += (0.6123724356957944*alpha[48]+0.6123724356957944*alpha[36]+0.6846531968814574*alpha[2])*dv10*jacob_vx_inv; 
  out[28] += (0.6123724356957944*alpha[105]+0.6846531968814574*alpha[23])*dv10*jacob_vx_inv; 
  out[29] += (0.6123724356957944*alpha[106]+0.6846531968814574*alpha[24])*dv10*jacob_vx_inv; 
  out[30] += (0.6123724356957944*alpha[50]+0.6123724356957944*alpha[39]+0.6846531968814574*alpha[4])*dv10*jacob_vx_inv; 
  out[31] += 0.6846531968814574*alpha[91]*dv10*jacob_vx_inv; 
  out[32] += 0.6846531968814574*alpha[92]*dv10*jacob_vx_inv; 
  out[33] += 0.6846531968814574*alpha[43]*dv10*jacob_vx_inv; 
  out[34] += 0.6846531968814574*alpha[44]*dv10*jacob_vx_inv; 
  out[35] += 0.6123724356957944*alpha[26]*dv10*jacob_vx_inv; 
  out[36] += 0.6123724356957944*alpha[27]*dv10*jacob_vx_inv; 
  out[37] += 0.6846531968814574*alpha[96]*dv10*jacob_vx_inv; 
  out[38] += 0.6846531968814574*alpha[97]*dv10*jacob_vx_inv; 
  out[39] += 0.6123724356957944*alpha[30]*dv10*jacob_vx_inv; 
  out[40] += 0.6846531968814574*alpha[101]*dv10*jacob_vx_inv; 
  out[41] += 0.6846531968814574*alpha[102]*dv10*jacob_vx_inv; 
  out[42] += 0.6846531968814574*alpha[46]*dv10*jacob_vx_inv; 
  out[43] += 0.6846531968814574*alpha[33]*dv10*jacob_vx_inv; 
  out[44] += 0.6846531968814574*alpha[34]*dv10*jacob_vx_inv; 
  out[45] += (0.5477225575051661*alpha[49]+0.6123724356957944*alpha[3])*dv10*jacob_vx_inv; 
  out[46] += 0.6846531968814574*alpha[42]*dv10*jacob_vx_inv; 
  out[47] += 0.6123724356957944*alpha[26]*dv10*jacob_vx_inv; 
  out[48] += 0.6123724356957944*alpha[27]*dv10*jacob_vx_inv; 
  out[49] += (0.5477225575051661*alpha[45]+0.6123724356957944*alpha[5])*dv10*jacob_vx_inv; 
  out[50] += 0.6123724356957944*alpha[30]*dv10*jacob_vx_inv; 
  out[51] += (0.6123724356957944*alpha[109]+0.6846531968814574*alpha[53])*dv10*jacob_vx_inv; 
  out[52] += (0.6123724356957944*alpha[80]+0.6123724356957944*alpha[58]+0.6846531968814574*alpha[6])*dv10*jacob_vx_inv; 
  out[53] += (0.6123724356957944*alpha[111]+0.6846531968814574*alpha[51])*dv10*jacob_vx_inv; 
  out[54] += (0.6123724356957944*alpha[83]+0.6123724356957944*alpha[63]+0.6846531968814574*alpha[9])*dv10*jacob_vx_inv; 
  out[55] += (0.6123724356957944*alpha[84]+0.6123724356957944*alpha[64]+0.6846531968814574*alpha[10])*dv10*jacob_vx_inv; 
  out[56] += 0.6846531968814574*alpha[68]*dv10*jacob_vx_inv; 
  out[57] += 0.6846531968814574*alpha[69]*dv10*jacob_vx_inv; 
  out[58] += 0.6123724356957944*alpha[52]*dv10*jacob_vx_inv; 
  out[59] += 0.6846531968814574*alpha[107]*dv10*jacob_vx_inv; 
  out[60] += 0.6846531968814574*alpha[108]*dv10*jacob_vx_inv; 
  out[61] += 0.6846531968814574*alpha[74]*dv10*jacob_vx_inv; 
  out[62] += 0.6846531968814574*alpha[75]*dv10*jacob_vx_inv; 
  out[63] += 0.6123724356957944*alpha[54]*dv10*jacob_vx_inv; 
  out[64] += 0.6123724356957944*alpha[55]*dv10*jacob_vx_inv; 
  out[65] += 0.6846531968814574*alpha[110]*dv10*jacob_vx_inv; 
  out[66] += 0.6846531968814574*alpha[77]*dv10*jacob_vx_inv; 
  out[67] += 0.6846531968814574*alpha[78]*dv10*jacob_vx_inv; 
  out[68] += 0.6846531968814574*alpha[56]*dv10*jacob_vx_inv; 
  out[69] += 0.6846531968814574*alpha[57]*dv10*jacob_vx_inv; 
  out[70] += 0.6846531968814574*alpha[16]*dv10*jacob_vx_inv; 
  out[71] += 0.6846531968814574*alpha[17]*dv10*jacob_vx_inv; 
  out[72] += (0.5477225575051661*alpha[81]+0.6123724356957944*alpha[7])*dv10*jacob_vx_inv; 
  out[73] += (0.5477225575051661*alpha[82]+0.6123724356957944*alpha[8])*dv10*jacob_vx_inv; 
  out[74] += 0.6846531968814574*alpha[61]*dv10*jacob_vx_inv; 
  out[75] += 0.6846531968814574*alpha[62]*dv10*jacob_vx_inv; 
  out[76] += (0.5477225575051661*alpha[85]+0.6123724356957944*alpha[11])*dv10*jacob_vx_inv; 
  out[77] += 0.6846531968814574*alpha[66]*dv10*jacob_vx_inv; 
  out[78] += 0.6846531968814574*alpha[67]*dv10*jacob_vx_inv; 
  out[79] += 0.6846531968814574*alpha[19]*dv10*jacob_vx_inv; 
  out[80] += 0.6123724356957944*alpha[52]*dv10*jacob_vx_inv; 
  out[81] += (0.5477225575051661*alpha[72]+0.6123724356957944*alpha[12])*dv10*jacob_vx_inv; 
  out[82] += (0.5477225575051661*alpha[73]+0.6123724356957944*alpha[13])*dv10*jacob_vx_inv; 
  out[83] += 0.6123724356957944*alpha[54]*dv10*jacob_vx_inv; 
  out[84] += 0.6123724356957944*alpha[55]*dv10*jacob_vx_inv; 
  out[85] += (0.5477225575051661*alpha[76]+0.6123724356957944*alpha[15])*dv10*jacob_vx_inv; 
  out[86] += (0.6123724356957944*alpha[104]+0.6123724356957944*alpha[89]+0.6846531968814574*alpha[22])*dv10*jacob_vx_inv; 
  out[87] += 0.6846531968814574*alpha[94]*dv10*jacob_vx_inv; 
  out[88] += 0.6846531968814574*alpha[95]*dv10*jacob_vx_inv; 
  out[89] += 0.6123724356957944*alpha[86]*dv10*jacob_vx_inv; 
  out[90] += 0.6846531968814574*alpha[100]*dv10*jacob_vx_inv; 
  out[91] += 0.6846531968814574*alpha[31]*dv10*jacob_vx_inv; 
  out[92] += 0.6846531968814574*alpha[32]*dv10*jacob_vx_inv; 
  out[93] += (0.5477225575051661*alpha[103]+0.6123724356957944*alpha[21])*dv10*jacob_vx_inv; 
  out[94] += 0.6846531968814574*alpha[87]*dv10*jacob_vx_inv; 
  out[95] += 0.6846531968814574*alpha[88]*dv10*jacob_vx_inv; 
  out[96] += 0.6846531968814574*alpha[37]*dv10*jacob_vx_inv; 
  out[97] += 0.6846531968814574*alpha[38]*dv10*jacob_vx_inv; 
  out[98] += (0.5477225575051661*alpha[105]+0.6123724356957944*alpha[23])*dv10*jacob_vx_inv; 
  out[99] += (0.5477225575051661*alpha[106]+0.6123724356957944*alpha[24])*dv10*jacob_vx_inv; 
  out[100] += 0.6846531968814574*alpha[90]*dv10*jacob_vx_inv; 
  out[101] += 0.6846531968814574*alpha[40]*dv10*jacob_vx_inv; 
  out[102] += 0.6846531968814574*alpha[41]*dv10*jacob_vx_inv; 
  out[103] += (0.5477225575051661*alpha[93]+0.6123724356957944*alpha[25])*dv10*jacob_vx_inv; 
  out[104] += 0.6123724356957944*alpha[86]*dv10*jacob_vx_inv; 
  out[105] += (0.5477225575051661*alpha[98]+0.6123724356957944*alpha[28])*dv10*jacob_vx_inv; 
  out[106] += (0.5477225575051661*alpha[99]+0.6123724356957944*alpha[29])*dv10*jacob_vx_inv; 
  out[107] += 0.6846531968814574*alpha[59]*dv10*jacob_vx_inv; 
  out[108] += 0.6846531968814574*alpha[60]*dv10*jacob_vx_inv; 
  out[109] += (0.5477225575051661*alpha[111]+0.6123724356957944*alpha[51])*dv10*jacob_vx_inv; 
  out[110] += 0.6846531968814574*alpha[65]*dv10*jacob_vx_inv; 
  out[111] += (0.5477225575051661*alpha[109]+0.6123724356957944*alpha[53])*dv10*jacob_vx_inv; 
  } 

  if (dir == 3) { 
  } 

  if (dir == 4) { 
  out[0] += 0.3061862178478971*alpha[18]*dv12*jacob_vz_inv; 
  out[1] += 0.3061862178478971*alpha[35]*dv12*jacob_vz_inv; 
  out[2] += 0.3061862178478971*alpha[36]*dv12*jacob_vz_inv; 
  out[3] += 0.27386127875258304*alpha[3]*dv12*jacob_vz_inv; 
  out[4] += 0.3061862178478971*alpha[39]*dv12*jacob_vz_inv; 
  out[5] += 0.3061862178478971*alpha[45]*dv12*jacob_vz_inv; 
  out[6] += 0.3061862178478971*alpha[58]*dv12*jacob_vz_inv; 
  out[7] += 0.27386127875258304*alpha[7]*dv12*jacob_vz_inv; 
  out[8] += 0.27386127875258304*alpha[8]*dv12*jacob_vz_inv; 
  out[9] += 0.3061862178478971*alpha[63]*dv12*jacob_vz_inv; 
  out[10] += 0.3061862178478971*alpha[64]*dv12*jacob_vz_inv; 
  out[11] += 0.27386127875258304*alpha[11]*dv12*jacob_vz_inv; 
  out[12] += 0.3061862178478971*alpha[72]*dv12*jacob_vz_inv; 
  out[13] += 0.3061862178478971*alpha[73]*dv12*jacob_vz_inv; 
  out[14] += 0.27386127875258304*alpha[14]*dv12*jacob_vz_inv; 
  out[15] += 0.3061862178478971*alpha[76]*dv12*jacob_vz_inv; 
  out[18] += (0.19561519910898784*alpha[18]+0.3061862178478971*alpha[0])*dv12*jacob_vz_inv; 
  out[21] += 0.27386127875258304*alpha[21]*dv12*jacob_vz_inv; 
  out[22] += 0.3061862178478971*alpha[89]*dv12*jacob_vz_inv; 
  out[23] += 0.27386127875258304*alpha[23]*dv12*jacob_vz_inv; 
  out[24] += 0.27386127875258304*alpha[24]*dv12*jacob_vz_inv; 
  out[25] += 0.3061862178478971*alpha[93]*dv12*jacob_vz_inv; 
  out[26] += 0.27386127875258304*alpha[26]*dv12*jacob_vz_inv; 
  out[27] += 0.27386127875258304*alpha[27]*dv12*jacob_vz_inv; 
  out[28] += 0.3061862178478971*alpha[98]*dv12*jacob_vz_inv; 
  out[29] += 0.3061862178478971*alpha[99]*dv12*jacob_vz_inv; 
  out[30] += 0.27386127875258304*alpha[30]*dv12*jacob_vz_inv; 
  out[33] += 0.27386127875258304*alpha[33]*dv12*jacob_vz_inv; 
  out[34] += 0.27386127875258304*alpha[34]*dv12*jacob_vz_inv; 
  out[35] += (0.19561519910898784*alpha[35]+0.3061862178478971*alpha[1])*dv12*jacob_vz_inv; 
  out[36] += (0.19561519910898784*alpha[36]+0.3061862178478971*alpha[2])*dv12*jacob_vz_inv; 
  out[39] += (0.19561519910898784*alpha[39]+0.3061862178478971*alpha[4])*dv12*jacob_vz_inv; 
  out[42] += 0.27386127875258304*alpha[42]*dv12*jacob_vz_inv; 
  out[45] += (0.19561519910898784*alpha[45]+0.3061862178478971*alpha[5])*dv12*jacob_vz_inv; 
  out[49] += 0.27386127875258304*alpha[49]*dv12*jacob_vz_inv; 
  out[51] += 0.27386127875258304*alpha[51]*dv12*jacob_vz_inv; 
  out[52] += 0.27386127875258304*alpha[52]*dv12*jacob_vz_inv; 
  out[53] += 0.3061862178478971*alpha[109]*dv12*jacob_vz_inv; 
  out[54] += 0.27386127875258304*alpha[54]*dv12*jacob_vz_inv; 
  out[55] += 0.27386127875258304*alpha[55]*dv12*jacob_vz_inv; 
  out[56] += 0.27386127875258304*alpha[56]*dv12*jacob_vz_inv; 
  out[57] += 0.27386127875258304*alpha[57]*dv12*jacob_vz_inv; 
  out[58] += (0.19561519910898784*alpha[58]+0.3061862178478971*alpha[6])*dv12*jacob_vz_inv; 
  out[61] += 0.27386127875258304*alpha[61]*dv12*jacob_vz_inv; 
  out[62] += 0.27386127875258304*alpha[62]*dv12*jacob_vz_inv; 
  out[63] += (0.19561519910898784*alpha[63]+0.3061862178478971*alpha[9])*dv12*jacob_vz_inv; 
  out[64] += (0.19561519910898784*alpha[64]+0.3061862178478971*alpha[10])*dv12*jacob_vz_inv; 
  out[66] += 0.27386127875258304*alpha[66]*dv12*jacob_vz_inv; 
  out[67] += 0.27386127875258304*alpha[67]*dv12*jacob_vz_inv; 
  out[70] += 0.27386127875258304*alpha[70]*dv12*jacob_vz_inv; 
  out[71] += 0.27386127875258304*alpha[71]*dv12*jacob_vz_inv; 
  out[72] += (0.19561519910898784*alpha[72]+0.3061862178478971*alpha[12])*dv12*jacob_vz_inv; 
  out[73] += (0.19561519910898784*alpha[73]+0.3061862178478971*alpha[13])*dv12*jacob_vz_inv; 
  out[76] += (0.19561519910898784*alpha[76]+0.3061862178478971*alpha[15])*dv12*jacob_vz_inv; 
  out[79] += 0.27386127875258304*alpha[79]*dv12*jacob_vz_inv; 
  out[81] += 0.27386127875258304*alpha[81]*dv12*jacob_vz_inv; 
  out[82] += 0.27386127875258304*alpha[82]*dv12*jacob_vz_inv; 
  out[85] += 0.27386127875258304*alpha[85]*dv12*jacob_vz_inv; 
  out[86] += 0.27386127875258304*alpha[86]*dv12*jacob_vz_inv; 
  out[87] += 0.27386127875258304*alpha[87]*dv12*jacob_vz_inv; 
  out[88] += 0.27386127875258304*alpha[88]*dv12*jacob_vz_inv; 
  out[89] += (0.19561519910898784*alpha[89]+0.3061862178478971*alpha[22])*dv12*jacob_vz_inv; 
  out[90] += 0.27386127875258304*alpha[90]*dv12*jacob_vz_inv; 
  out[91] += 0.27386127875258304*alpha[91]*dv12*jacob_vz_inv; 
  out[92] += 0.27386127875258304*alpha[92]*dv12*jacob_vz_inv; 
  out[93] += (0.19561519910898784*alpha[93]+0.3061862178478971*alpha[25])*dv12*jacob_vz_inv; 
  out[96] += 0.27386127875258304*alpha[96]*dv12*jacob_vz_inv; 
  out[97] += 0.27386127875258304*alpha[97]*dv12*jacob_vz_inv; 
  out[98] += (0.19561519910898784*alpha[98]+0.3061862178478971*alpha[28])*dv12*jacob_vz_inv; 
  out[99] += (0.19561519910898784*alpha[99]+0.3061862178478971*alpha[29])*dv12*jacob_vz_inv; 
  out[101] += 0.27386127875258304*alpha[101]*dv12*jacob_vz_inv; 
  out[102] += 0.27386127875258304*alpha[102]*dv12*jacob_vz_inv; 
  out[103] += 0.27386127875258304*alpha[103]*dv12*jacob_vz_inv; 
  out[105] += 0.27386127875258304*alpha[105]*dv12*jacob_vz_inv; 
  out[106] += 0.27386127875258304*alpha[106]*dv12*jacob_vz_inv; 
  out[107] += 0.27386127875258304*alpha[107]*dv12*jacob_vz_inv; 
  out[108] += 0.27386127875258304*alpha[108]*dv12*jacob_vz_inv; 
  out[109] += (0.19561519910898784*alpha[109]+0.3061862178478971*alpha[53])*dv12*jacob_vz_inv; 
  out[110] += 0.27386127875258304*alpha[110]*dv12*jacob_vz_inv; 
  out[111] += 0.27386127875258304*alpha[111]*dv12*jacob_vz_inv; 
  } 

} 
