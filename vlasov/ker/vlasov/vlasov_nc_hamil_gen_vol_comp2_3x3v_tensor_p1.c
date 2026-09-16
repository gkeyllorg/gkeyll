#include <gkyl_nc_hamil_gen_vol_comps_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp2_3x3v_tensor_p1(const double *w, const double *dxv, const int dir,
   const double *jacob_pos, const double *jacob_vel,
   const double *alpha, double* GKYL_RESTRICT out) 
{ 
  const double dx10 = 2.0/dxv[0]; 
  const double dx11 = 2.0/dxv[1]; 
  const double dx12 = 2.0/dxv[2]; 
  const double dv10 = 2.0/dxv[3]; 
  const double dv11 = 2.0/dxv[4]; 
  const double dv12 = 2.0/dxv[5]; 
  const double jacob_cx_inv = 1.0/jacob_pos[0]; 
  const double jacob_cy_inv = 1.0/jacob_pos[2]; 
  const double jacob_cz_inv = 1.0/jacob_pos[4]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double *jacob_vy = &jacob_vel[3]; 
  const double *jacob_vz = &jacob_vel[6]; 
  
  if (dir == 0) { 
  } 

  if (dir == 1) { 
  out[0] += 0.21650635094610965*alpha[0]*dx11*jacob_cy_inv; 
  out[1] += 0.21650635094610965*alpha[1]*dx11*jacob_cy_inv; 
  out[2] += 0.21650635094610965*alpha[2]*dx11*jacob_cy_inv; 
  out[3] += 0.21650635094610965*alpha[3]*dx11*jacob_cy_inv; 
  out[4] += 0.21650635094610965*alpha[4]*dx11*jacob_cy_inv; 
  out[5] += 0.21650635094610965*alpha[5]*dx11*jacob_cy_inv; 
  out[6] += 0.21650635094610965*alpha[6]*dx11*jacob_cy_inv; 
  out[7] += 0.21650635094610965*alpha[7]*dx11*jacob_cy_inv; 
  out[8] += 0.21650635094610965*alpha[8]*dx11*jacob_cy_inv; 
  out[9] += 0.21650635094610965*alpha[9]*dx11*jacob_cy_inv; 
  out[10] += 0.21650635094610965*alpha[10]*dx11*jacob_cy_inv; 
  out[11] += 0.21650635094610965*alpha[11]*dx11*jacob_cy_inv; 
  out[12] += 0.21650635094610965*alpha[12]*dx11*jacob_cy_inv; 
  out[13] += 0.21650635094610965*alpha[13]*dx11*jacob_cy_inv; 
  out[14] += 0.21650635094610965*alpha[14]*dx11*jacob_cy_inv; 
  out[15] += 0.21650635094610965*alpha[15]*dx11*jacob_cy_inv; 
  out[16] += 0.21650635094610965*alpha[16]*dx11*jacob_cy_inv; 
  out[17] += 0.21650635094610965*alpha[17]*dx11*jacob_cy_inv; 
  out[18] += 0.21650635094610965*alpha[18]*dx11*jacob_cy_inv; 
  out[19] += 0.21650635094610965*alpha[19]*dx11*jacob_cy_inv; 
  out[20] += 0.21650635094610965*alpha[20]*dx11*jacob_cy_inv; 
  out[21] += 0.21650635094610965*alpha[21]*dx11*jacob_cy_inv; 
  out[22] += 0.21650635094610965*alpha[22]*dx11*jacob_cy_inv; 
  out[23] += 0.21650635094610965*alpha[23]*dx11*jacob_cy_inv; 
  out[24] += 0.21650635094610965*alpha[24]*dx11*jacob_cy_inv; 
  out[25] += 0.21650635094610965*alpha[25]*dx11*jacob_cy_inv; 
  out[26] += 0.21650635094610965*alpha[26]*dx11*jacob_cy_inv; 
  out[27] += 0.21650635094610965*alpha[27]*dx11*jacob_cy_inv; 
  out[28] += 0.21650635094610965*alpha[28]*dx11*jacob_cy_inv; 
  out[29] += 0.21650635094610965*alpha[29]*dx11*jacob_cy_inv; 
  out[30] += 0.21650635094610965*alpha[30]*dx11*jacob_cy_inv; 
  out[31] += 0.21650635094610965*alpha[31]*dx11*jacob_cy_inv; 
  out[32] += 0.21650635094610965*alpha[32]*dx11*jacob_cy_inv; 
  out[33] += 0.21650635094610965*alpha[33]*dx11*jacob_cy_inv; 
  out[34] += 0.21650635094610965*alpha[34]*dx11*jacob_cy_inv; 
  out[35] += 0.21650635094610965*alpha[35]*dx11*jacob_cy_inv; 
  out[36] += 0.21650635094610965*alpha[36]*dx11*jacob_cy_inv; 
  out[37] += 0.21650635094610965*alpha[37]*dx11*jacob_cy_inv; 
  out[38] += 0.21650635094610965*alpha[38]*dx11*jacob_cy_inv; 
  out[39] += 0.21650635094610965*alpha[39]*dx11*jacob_cy_inv; 
  out[40] += 0.21650635094610965*alpha[40]*dx11*jacob_cy_inv; 
  out[41] += 0.21650635094610965*alpha[41]*dx11*jacob_cy_inv; 
  out[42] += 0.21650635094610965*alpha[42]*dx11*jacob_cy_inv; 
  out[43] += 0.21650635094610965*alpha[43]*dx11*jacob_cy_inv; 
  out[44] += 0.21650635094610965*alpha[44]*dx11*jacob_cy_inv; 
  out[45] += 0.21650635094610965*alpha[45]*dx11*jacob_cy_inv; 
  out[46] += 0.21650635094610965*alpha[46]*dx11*jacob_cy_inv; 
  out[47] += 0.21650635094610965*alpha[47]*dx11*jacob_cy_inv; 
  out[48] += 0.21650635094610965*alpha[48]*dx11*jacob_cy_inv; 
  out[49] += 0.21650635094610965*alpha[49]*dx11*jacob_cy_inv; 
  out[50] += 0.21650635094610965*alpha[50]*dx11*jacob_cy_inv; 
  out[51] += 0.21650635094610965*alpha[51]*dx11*jacob_cy_inv; 
  out[52] += 0.21650635094610965*alpha[52]*dx11*jacob_cy_inv; 
  out[53] += 0.21650635094610965*alpha[53]*dx11*jacob_cy_inv; 
  out[54] += 0.21650635094610965*alpha[54]*dx11*jacob_cy_inv; 
  out[55] += 0.21650635094610965*alpha[55]*dx11*jacob_cy_inv; 
  out[56] += 0.21650635094610965*alpha[56]*dx11*jacob_cy_inv; 
  out[57] += 0.21650635094610965*alpha[57]*dx11*jacob_cy_inv; 
  out[58] += 0.21650635094610965*alpha[58]*dx11*jacob_cy_inv; 
  out[59] += 0.21650635094610965*alpha[59]*dx11*jacob_cy_inv; 
  out[60] += 0.21650635094610965*alpha[60]*dx11*jacob_cy_inv; 
  out[61] += 0.21650635094610965*alpha[61]*dx11*jacob_cy_inv; 
  out[62] += 0.21650635094610965*alpha[62]*dx11*jacob_cy_inv; 
  out[63] += 0.21650635094610965*alpha[63]*dx11*jacob_cy_inv; 
  out[64] += 0.21650635094610965*alpha[64]*dx11*jacob_cy_inv; 
  out[65] += 0.21650635094610965*alpha[65]*dx11*jacob_cy_inv; 
  out[66] += 0.21650635094610965*alpha[66]*dx11*jacob_cy_inv; 
  out[67] += 0.21650635094610965*alpha[67]*dx11*jacob_cy_inv; 
  out[68] += 0.21650635094610965*alpha[68]*dx11*jacob_cy_inv; 
  out[69] += 0.21650635094610965*alpha[69]*dx11*jacob_cy_inv; 
  out[70] += 0.21650635094610965*alpha[70]*dx11*jacob_cy_inv; 
  out[71] += 0.21650635094610965*alpha[71]*dx11*jacob_cy_inv; 
  out[72] += 0.21650635094610965*alpha[72]*dx11*jacob_cy_inv; 
  out[73] += 0.21650635094610965*alpha[73]*dx11*jacob_cy_inv; 
  out[74] += 0.21650635094610965*alpha[74]*dx11*jacob_cy_inv; 
  out[75] += 0.21650635094610965*alpha[75]*dx11*jacob_cy_inv; 
  out[76] += 0.21650635094610965*alpha[76]*dx11*jacob_cy_inv; 
  out[77] += 0.21650635094610965*alpha[77]*dx11*jacob_cy_inv; 
  out[78] += 0.21650635094610965*alpha[78]*dx11*jacob_cy_inv; 
  out[79] += 0.21650635094610965*alpha[79]*dx11*jacob_cy_inv; 
  out[80] += 0.21650635094610965*alpha[80]*dx11*jacob_cy_inv; 
  out[81] += 0.21650635094610965*alpha[81]*dx11*jacob_cy_inv; 
  out[82] += 0.21650635094610965*alpha[82]*dx11*jacob_cy_inv; 
  out[83] += 0.21650635094610965*alpha[83]*dx11*jacob_cy_inv; 
  out[84] += 0.21650635094610965*alpha[84]*dx11*jacob_cy_inv; 
  out[85] += 0.21650635094610965*alpha[85]*dx11*jacob_cy_inv; 
  out[86] += 0.21650635094610965*alpha[86]*dx11*jacob_cy_inv; 
  out[87] += 0.21650635094610965*alpha[87]*dx11*jacob_cy_inv; 
  out[88] += 0.21650635094610965*alpha[88]*dx11*jacob_cy_inv; 
  out[89] += 0.21650635094610965*alpha[89]*dx11*jacob_cy_inv; 
  out[90] += 0.21650635094610965*alpha[90]*dx11*jacob_cy_inv; 
  out[91] += 0.21650635094610965*alpha[91]*dx11*jacob_cy_inv; 
  out[92] += 0.21650635094610965*alpha[92]*dx11*jacob_cy_inv; 
  out[93] += 0.21650635094610965*alpha[93]*dx11*jacob_cy_inv; 
  out[94] += 0.21650635094610965*alpha[94]*dx11*jacob_cy_inv; 
  out[95] += 0.21650635094610965*alpha[95]*dx11*jacob_cy_inv; 
  out[96] += 0.21650635094610965*alpha[96]*dx11*jacob_cy_inv; 
  out[97] += 0.21650635094610965*alpha[97]*dx11*jacob_cy_inv; 
  out[98] += 0.21650635094610965*alpha[98]*dx11*jacob_cy_inv; 
  out[99] += 0.21650635094610965*alpha[99]*dx11*jacob_cy_inv; 
  out[100] += 0.21650635094610965*alpha[100]*dx11*jacob_cy_inv; 
  out[101] += 0.21650635094610965*alpha[101]*dx11*jacob_cy_inv; 
  out[102] += 0.21650635094610965*alpha[102]*dx11*jacob_cy_inv; 
  out[103] += 0.21650635094610965*alpha[103]*dx11*jacob_cy_inv; 
  out[104] += 0.21650635094610965*alpha[104]*dx11*jacob_cy_inv; 
  out[105] += 0.21650635094610965*alpha[105]*dx11*jacob_cy_inv; 
  out[106] += 0.21650635094610965*alpha[106]*dx11*jacob_cy_inv; 
  out[107] += 0.21650635094610965*alpha[107]*dx11*jacob_cy_inv; 
  out[108] += 0.21650635094610965*alpha[108]*dx11*jacob_cy_inv; 
  out[109] += 0.21650635094610965*alpha[109]*dx11*jacob_cy_inv; 
  out[110] += 0.21650635094610965*alpha[110]*dx11*jacob_cy_inv; 
  out[111] += 0.21650635094610965*alpha[111]*dx11*jacob_cy_inv; 
  out[112] += 0.21650635094610965*alpha[112]*dx11*jacob_cy_inv; 
  out[113] += 0.21650635094610965*alpha[113]*dx11*jacob_cy_inv; 
  out[114] += 0.21650635094610965*alpha[114]*dx11*jacob_cy_inv; 
  out[115] += 0.21650635094610965*alpha[115]*dx11*jacob_cy_inv; 
  out[116] += 0.21650635094610965*alpha[116]*dx11*jacob_cy_inv; 
  out[117] += 0.21650635094610965*alpha[117]*dx11*jacob_cy_inv; 
  out[118] += 0.21650635094610965*alpha[118]*dx11*jacob_cy_inv; 
  out[119] += 0.21650635094610965*alpha[119]*dx11*jacob_cy_inv; 
  out[120] += 0.21650635094610965*alpha[120]*dx11*jacob_cy_inv; 
  out[121] += 0.21650635094610965*alpha[121]*dx11*jacob_cy_inv; 
  out[122] += 0.21650635094610965*alpha[122]*dx11*jacob_cy_inv; 
  out[123] += 0.21650635094610965*alpha[123]*dx11*jacob_cy_inv; 
  out[124] += 0.21650635094610965*alpha[124]*dx11*jacob_cy_inv; 
  out[125] += 0.21650635094610965*alpha[125]*dx11*jacob_cy_inv; 
  out[126] += 0.21650635094610965*alpha[126]*dx11*jacob_cy_inv; 
  out[127] += 0.21650635094610965*alpha[127]*dx11*jacob_cy_inv; 
  out[128] += 0.21650635094610965*alpha[128]*dx11*jacob_cy_inv; 
  out[129] += 0.21650635094610965*alpha[129]*dx11*jacob_cy_inv; 
  out[130] += 0.21650635094610965*alpha[130]*dx11*jacob_cy_inv; 
  out[131] += 0.21650635094610965*alpha[131]*dx11*jacob_cy_inv; 
  out[132] += 0.21650635094610965*alpha[132]*dx11*jacob_cy_inv; 
  out[133] += 0.21650635094610965*alpha[133]*dx11*jacob_cy_inv; 
  out[134] += 0.21650635094610965*alpha[134]*dx11*jacob_cy_inv; 
  out[135] += 0.21650635094610965*alpha[135]*dx11*jacob_cy_inv; 
  out[136] += 0.21650635094610965*alpha[136]*dx11*jacob_cy_inv; 
  out[137] += 0.21650635094610965*alpha[137]*dx11*jacob_cy_inv; 
  out[138] += 0.21650635094610965*alpha[138]*dx11*jacob_cy_inv; 
  out[139] += 0.21650635094610965*alpha[139]*dx11*jacob_cy_inv; 
  out[140] += 0.21650635094610965*alpha[140]*dx11*jacob_cy_inv; 
  out[141] += 0.21650635094610965*alpha[141]*dx11*jacob_cy_inv; 
  out[142] += 0.21650635094610965*alpha[142]*dx11*jacob_cy_inv; 
  out[143] += 0.21650635094610965*alpha[143]*dx11*jacob_cy_inv; 
  out[144] += 0.21650635094610965*alpha[144]*dx11*jacob_cy_inv; 
  out[145] += 0.21650635094610965*alpha[145]*dx11*jacob_cy_inv; 
  out[146] += 0.21650635094610965*alpha[146]*dx11*jacob_cy_inv; 
  out[147] += 0.21650635094610965*alpha[147]*dx11*jacob_cy_inv; 
  out[148] += 0.21650635094610965*alpha[148]*dx11*jacob_cy_inv; 
  out[149] += 0.21650635094610965*alpha[149]*dx11*jacob_cy_inv; 
  out[150] += 0.21650635094610965*alpha[150]*dx11*jacob_cy_inv; 
  out[151] += 0.21650635094610965*alpha[151]*dx11*jacob_cy_inv; 
  out[152] += 0.21650635094610965*alpha[152]*dx11*jacob_cy_inv; 
  out[153] += 0.21650635094610965*alpha[153]*dx11*jacob_cy_inv; 
  out[154] += 0.21650635094610965*alpha[154]*dx11*jacob_cy_inv; 
  out[155] += 0.21650635094610965*alpha[155]*dx11*jacob_cy_inv; 
  out[156] += 0.21650635094610965*alpha[156]*dx11*jacob_cy_inv; 
  out[157] += 0.21650635094610965*alpha[157]*dx11*jacob_cy_inv; 
  out[158] += 0.21650635094610965*alpha[158]*dx11*jacob_cy_inv; 
  out[159] += 0.21650635094610965*alpha[159]*dx11*jacob_cy_inv; 
  out[160] += 0.21650635094610965*alpha[160]*dx11*jacob_cy_inv; 
  out[161] += 0.21650635094610965*alpha[161]*dx11*jacob_cy_inv; 
  out[162] += 0.21650635094610965*alpha[162]*dx11*jacob_cy_inv; 
  out[163] += 0.21650635094610965*alpha[163]*dx11*jacob_cy_inv; 
  out[164] += 0.21650635094610965*alpha[164]*dx11*jacob_cy_inv; 
  out[165] += 0.21650635094610965*alpha[165]*dx11*jacob_cy_inv; 
  out[166] += 0.21650635094610965*alpha[166]*dx11*jacob_cy_inv; 
  out[167] += 0.21650635094610965*alpha[167]*dx11*jacob_cy_inv; 
  out[168] += 0.21650635094610965*alpha[168]*dx11*jacob_cy_inv; 
  out[169] += 0.21650635094610965*alpha[169]*dx11*jacob_cy_inv; 
  out[170] += 0.21650635094610965*alpha[170]*dx11*jacob_cy_inv; 
  out[171] += 0.21650635094610965*alpha[171]*dx11*jacob_cy_inv; 
  out[172] += 0.21650635094610965*alpha[172]*dx11*jacob_cy_inv; 
  out[173] += 0.21650635094610965*alpha[173]*dx11*jacob_cy_inv; 
  out[174] += 0.21650635094610965*alpha[174]*dx11*jacob_cy_inv; 
  out[175] += 0.21650635094610965*alpha[175]*dx11*jacob_cy_inv; 
  out[176] += 0.21650635094610965*alpha[176]*dx11*jacob_cy_inv; 
  out[177] += 0.21650635094610965*alpha[177]*dx11*jacob_cy_inv; 
  out[178] += 0.21650635094610965*alpha[178]*dx11*jacob_cy_inv; 
  out[179] += 0.21650635094610965*alpha[179]*dx11*jacob_cy_inv; 
  out[180] += 0.21650635094610965*alpha[180]*dx11*jacob_cy_inv; 
  out[181] += 0.21650635094610965*alpha[181]*dx11*jacob_cy_inv; 
  out[182] += 0.21650635094610965*alpha[182]*dx11*jacob_cy_inv; 
  out[183] += 0.21650635094610965*alpha[183]*dx11*jacob_cy_inv; 
  out[184] += 0.21650635094610965*alpha[184]*dx11*jacob_cy_inv; 
  out[185] += 0.21650635094610965*alpha[185]*dx11*jacob_cy_inv; 
  out[186] += 0.21650635094610965*alpha[186]*dx11*jacob_cy_inv; 
  out[187] += 0.21650635094610965*alpha[187]*dx11*jacob_cy_inv; 
  out[188] += 0.21650635094610965*alpha[188]*dx11*jacob_cy_inv; 
  out[189] += 0.21650635094610965*alpha[189]*dx11*jacob_cy_inv; 
  out[190] += 0.21650635094610965*alpha[190]*dx11*jacob_cy_inv; 
  out[191] += 0.21650635094610965*alpha[191]*dx11*jacob_cy_inv; 
  out[192] += 0.21650635094610965*alpha[192]*dx11*jacob_cy_inv; 
  out[193] += 0.21650635094610965*alpha[193]*dx11*jacob_cy_inv; 
  out[194] += 0.21650635094610965*alpha[194]*dx11*jacob_cy_inv; 
  out[195] += 0.21650635094610965*alpha[195]*dx11*jacob_cy_inv; 
  out[196] += 0.21650635094610965*alpha[196]*dx11*jacob_cy_inv; 
  out[197] += 0.21650635094610965*alpha[197]*dx11*jacob_cy_inv; 
  out[198] += 0.21650635094610965*alpha[198]*dx11*jacob_cy_inv; 
  out[199] += 0.21650635094610965*alpha[199]*dx11*jacob_cy_inv; 
  out[200] += 0.21650635094610965*alpha[200]*dx11*jacob_cy_inv; 
  out[201] += 0.21650635094610965*alpha[201]*dx11*jacob_cy_inv; 
  out[202] += 0.21650635094610965*alpha[202]*dx11*jacob_cy_inv; 
  out[203] += 0.21650635094610965*alpha[203]*dx11*jacob_cy_inv; 
  out[204] += 0.21650635094610965*alpha[204]*dx11*jacob_cy_inv; 
  out[205] += 0.21650635094610965*alpha[205]*dx11*jacob_cy_inv; 
  out[206] += 0.21650635094610965*alpha[206]*dx11*jacob_cy_inv; 
  out[207] += 0.21650635094610965*alpha[207]*dx11*jacob_cy_inv; 
  out[208] += 0.21650635094610965*alpha[208]*dx11*jacob_cy_inv; 
  out[209] += 0.21650635094610965*alpha[209]*dx11*jacob_cy_inv; 
  out[210] += 0.21650635094610965*alpha[210]*dx11*jacob_cy_inv; 
  out[211] += 0.21650635094610965*alpha[211]*dx11*jacob_cy_inv; 
  out[212] += 0.21650635094610965*alpha[212]*dx11*jacob_cy_inv; 
  out[213] += 0.21650635094610965*alpha[213]*dx11*jacob_cy_inv; 
  out[214] += 0.21650635094610965*alpha[214]*dx11*jacob_cy_inv; 
  out[215] += 0.21650635094610965*alpha[215]*dx11*jacob_cy_inv; 
  } 

  if (dir == 2) { 
  } 

  if (dir == 3) { 
  } 

  if (dir == 4) { 
  } 

  if (dir == 5) { 
  } 

} 
