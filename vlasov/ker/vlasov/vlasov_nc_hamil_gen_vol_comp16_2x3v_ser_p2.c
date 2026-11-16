#include <gkyl_nc_hamil_gen_vol_comps_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp16_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
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
  out[0] += 0.6846531968814573*alpha[1]*dx10*jacob_cx_inv; 
  out[1] += (0.6123724356957944*alpha[16]+0.6846531968814573*alpha[0])*dx10*jacob_cx_inv; 
  out[2] += 0.6846531968814573*alpha[6]*dx10*jacob_cx_inv; 
  out[3] += 0.6846531968814573*alpha[7]*dx10*jacob_cx_inv; 
  out[4] += 0.6846531968814573*alpha[9]*dx10*jacob_cx_inv; 
  out[5] += 0.6846531968814573*alpha[12]*dx10*jacob_cx_inv; 
  out[6] += (0.6123724356957944*alpha[31]+0.6846531968814573*alpha[2])*dx10*jacob_cx_inv; 
  out[7] += (0.6123724356957944*alpha[33]+0.6846531968814573*alpha[3])*dx10*jacob_cx_inv; 
  out[8] += 0.6846531968814573*alpha[21]*dx10*jacob_cx_inv; 
  out[9] += (0.6123724356957944*alpha[37]+0.6846531968814573*alpha[4])*dx10*jacob_cx_inv; 
  out[10] += 0.6846531968814573*alpha[22]*dx10*jacob_cx_inv; 
  out[11] += 0.6846531968814573*alpha[23]*dx10*jacob_cx_inv; 
  out[12] += (0.6123724356957944*alpha[43]+0.6846531968814573*alpha[5])*dx10*jacob_cx_inv; 
  out[13] += 0.6846531968814573*alpha[25]*dx10*jacob_cx_inv; 
  out[14] += 0.6846531968814573*alpha[26]*dx10*jacob_cx_inv; 
  out[15] += 0.6846531968814573*alpha[28]*dx10*jacob_cx_inv; 
  out[16] += 0.6123724356957944*alpha[1]*dx10*jacob_cx_inv; 
  out[17] += 0.6846531968814574*alpha[32]*dx10*jacob_cx_inv; 
  out[18] += 0.6846531968814574*alpha[35]*dx10*jacob_cx_inv; 
  out[19] += 0.6846531968814574*alpha[40]*dx10*jacob_cx_inv; 
  out[20] += 0.6846531968814574*alpha[47]*dx10*jacob_cx_inv; 
  out[21] += (0.6123724356957944*alpha[56]+0.6846531968814573*alpha[8])*dx10*jacob_cx_inv; 
  out[22] += (0.6123724356957944*alpha[59]+0.6846531968814573*alpha[10])*dx10*jacob_cx_inv; 
  out[23] += (0.6123724356957944*alpha[61]+0.6846531968814573*alpha[11])*dx10*jacob_cx_inv; 
  out[24] += 0.6846531968814573*alpha[51]*dx10*jacob_cx_inv; 
  out[25] += (0.6123724356957944*alpha[68]+0.6846531968814573*alpha[13])*dx10*jacob_cx_inv; 
  out[26] += (0.6123724356957944*alpha[70]+0.6846531968814573*alpha[14])*dx10*jacob_cx_inv; 
  out[27] += 0.6846531968814573*alpha[52]*dx10*jacob_cx_inv; 
  out[28] += (0.6123724356957944*alpha[74]+0.6846531968814573*alpha[15])*dx10*jacob_cx_inv; 
  out[29] += 0.6846531968814573*alpha[53]*dx10*jacob_cx_inv; 
  out[30] += 0.6846531968814573*alpha[54]*dx10*jacob_cx_inv; 
  out[31] += 0.6123724356957944*alpha[6]*dx10*jacob_cx_inv; 
  out[32] += 0.6846531968814574*alpha[17]*dx10*jacob_cx_inv; 
  out[33] += 0.6123724356957944*alpha[7]*dx10*jacob_cx_inv; 
  out[34] += 0.6846531968814574*alpha[57]*dx10*jacob_cx_inv; 
  out[35] += 0.6846531968814574*alpha[18]*dx10*jacob_cx_inv; 
  out[36] += 0.6846531968814574*alpha[58]*dx10*jacob_cx_inv; 
  out[37] += 0.6123724356957944*alpha[9]*dx10*jacob_cx_inv; 
  out[38] += 0.6846531968814574*alpha[60]*dx10*jacob_cx_inv; 
  out[39] += 0.6846531968814574*alpha[63]*dx10*jacob_cx_inv; 
  out[40] += 0.6846531968814574*alpha[19]*dx10*jacob_cx_inv; 
  out[41] += 0.6846531968814574*alpha[65]*dx10*jacob_cx_inv; 
  out[42] += 0.6846531968814574*alpha[66]*dx10*jacob_cx_inv; 
  out[43] += 0.6123724356957944*alpha[12]*dx10*jacob_cx_inv; 
  out[44] += 0.6846531968814574*alpha[69]*dx10*jacob_cx_inv; 
  out[45] += 0.6846531968814574*alpha[72]*dx10*jacob_cx_inv; 
  out[46] += 0.6846531968814574*alpha[77]*dx10*jacob_cx_inv; 
  out[47] += 0.6846531968814574*alpha[20]*dx10*jacob_cx_inv; 
  out[48] += 0.6846531968814574*alpha[80]*dx10*jacob_cx_inv; 
  out[49] += 0.6846531968814574*alpha[81]*dx10*jacob_cx_inv; 
  out[50] += 0.6846531968814574*alpha[83]*dx10*jacob_cx_inv; 
  out[51] += (0.6123724356957944*alpha[87]+0.6846531968814573*alpha[24])*dx10*jacob_cx_inv; 
  out[52] += (0.6123724356957944*alpha[91]+0.6846531968814573*alpha[27])*dx10*jacob_cx_inv; 
  out[53] += (0.6123724356957944*alpha[94]+0.6846531968814573*alpha[29])*dx10*jacob_cx_inv; 
  out[54] += (0.6123724356957944*alpha[96]+0.6846531968814573*alpha[30])*dx10*jacob_cx_inv; 
  out[55] += 0.6846531968814573*alpha[86]*dx10*jacob_cx_inv; 
  out[56] += 0.6123724356957944*alpha[21]*dx10*jacob_cx_inv; 
  out[57] += 0.6846531968814574*alpha[34]*dx10*jacob_cx_inv; 
  out[58] += 0.6846531968814574*alpha[36]*dx10*jacob_cx_inv; 
  out[59] += 0.6123724356957944*alpha[22]*dx10*jacob_cx_inv; 
  out[60] += 0.6846531968814574*alpha[38]*dx10*jacob_cx_inv; 
  out[61] += 0.6123724356957944*alpha[23]*dx10*jacob_cx_inv; 
  out[62] += 0.6846531968814574*alpha[88]*dx10*jacob_cx_inv; 
  out[63] += 0.6846531968814574*alpha[39]*dx10*jacob_cx_inv; 
  out[64] += 0.6846531968814574*alpha[89]*dx10*jacob_cx_inv; 
  out[65] += 0.6846531968814574*alpha[41]*dx10*jacob_cx_inv; 
  out[66] += 0.6846531968814574*alpha[42]*dx10*jacob_cx_inv; 
  out[67] += 0.6846531968814574*alpha[90]*dx10*jacob_cx_inv; 
  out[68] += 0.6123724356957944*alpha[25]*dx10*jacob_cx_inv; 
  out[69] += 0.6846531968814574*alpha[44]*dx10*jacob_cx_inv; 
  out[70] += 0.6123724356957944*alpha[26]*dx10*jacob_cx_inv; 
  out[71] += 0.6846531968814574*alpha[92]*dx10*jacob_cx_inv; 
  out[72] += 0.6846531968814574*alpha[45]*dx10*jacob_cx_inv; 
  out[73] += 0.6846531968814574*alpha[93]*dx10*jacob_cx_inv; 
  out[74] += 0.6123724356957944*alpha[28]*dx10*jacob_cx_inv; 
  out[75] += 0.6846531968814574*alpha[95]*dx10*jacob_cx_inv; 
  out[76] += 0.6846531968814574*alpha[98]*dx10*jacob_cx_inv; 
  out[77] += 0.6846531968814574*alpha[46]*dx10*jacob_cx_inv; 
  out[78] += 0.6846531968814574*alpha[100]*dx10*jacob_cx_inv; 
  out[79] += 0.6846531968814574*alpha[101]*dx10*jacob_cx_inv; 
  out[80] += 0.6846531968814574*alpha[48]*dx10*jacob_cx_inv; 
  out[81] += 0.6846531968814574*alpha[49]*dx10*jacob_cx_inv; 
  out[82] += 0.6846531968814574*alpha[103]*dx10*jacob_cx_inv; 
  out[83] += 0.6846531968814574*alpha[50]*dx10*jacob_cx_inv; 
  out[84] += 0.6846531968814574*alpha[104]*dx10*jacob_cx_inv; 
  out[85] += 0.6846531968814574*alpha[105]*dx10*jacob_cx_inv; 
  out[86] += (0.6123724356957944*alpha[107]+0.6846531968814573*alpha[55])*dx10*jacob_cx_inv; 
  out[87] += 0.6123724356957944*alpha[51]*dx10*jacob_cx_inv; 
  out[88] += 0.6846531968814574*alpha[62]*dx10*jacob_cx_inv; 
  out[89] += 0.6846531968814574*alpha[64]*dx10*jacob_cx_inv; 
  out[90] += 0.6846531968814574*alpha[67]*dx10*jacob_cx_inv; 
  out[91] += 0.6123724356957944*alpha[52]*dx10*jacob_cx_inv; 
  out[92] += 0.6846531968814574*alpha[71]*dx10*jacob_cx_inv; 
  out[93] += 0.6846531968814574*alpha[73]*dx10*jacob_cx_inv; 
  out[94] += 0.6123724356957944*alpha[53]*dx10*jacob_cx_inv; 
  out[95] += 0.6846531968814574*alpha[75]*dx10*jacob_cx_inv; 
  out[96] += 0.6123724356957944*alpha[54]*dx10*jacob_cx_inv; 
  out[97] += 0.6846531968814574*alpha[108]*dx10*jacob_cx_inv; 
  out[98] += 0.6846531968814574*alpha[76]*dx10*jacob_cx_inv; 
  out[99] += 0.6846531968814574*alpha[109]*dx10*jacob_cx_inv; 
  out[100] += 0.6846531968814574*alpha[78]*dx10*jacob_cx_inv; 
  out[101] += 0.6846531968814574*alpha[79]*dx10*jacob_cx_inv; 
  out[102] += 0.6846531968814574*alpha[110]*dx10*jacob_cx_inv; 
  out[103] += 0.6846531968814574*alpha[82]*dx10*jacob_cx_inv; 
  out[104] += 0.6846531968814574*alpha[84]*dx10*jacob_cx_inv; 
  out[105] += 0.6846531968814574*alpha[85]*dx10*jacob_cx_inv; 
  out[106] += 0.6846531968814574*alpha[111]*dx10*jacob_cx_inv; 
  out[107] += 0.6123724356957944*alpha[86]*dx10*jacob_cx_inv; 
  out[108] += 0.6846531968814574*alpha[97]*dx10*jacob_cx_inv; 
  out[109] += 0.6846531968814574*alpha[99]*dx10*jacob_cx_inv; 
  out[110] += 0.6846531968814574*alpha[102]*dx10*jacob_cx_inv; 
  out[111] += 0.6846531968814574*alpha[106]*dx10*jacob_cx_inv; 
  } 

  if (dir == 1) { 
  } 

  if (dir == 2) { 
  } 

  if (dir == 3) { 
  } 

  if (dir == 4) { 
  } 

} 
