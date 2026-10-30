#include <gkyl_array_average_kernels.h>

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.1767766952966368*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.1767766952966368*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.1767766952966368*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.1767766952966368*win[0]*fin[4]*subvol); 
  atomicAdd(&out[4], 0.1767766952966368*win[0]*fin[5]*subvol); 
  atomicAdd(&out[5], 0.1767766952966368*win[0]*fin[6]*subvol); 
  atomicAdd(&out[6], 0.1767766952966368*win[0]*fin[9]*subvol); 
  atomicAdd(&out[7], 0.1767766952966368*win[0]*fin[11]*subvol); 
  atomicAdd(&out[8], 0.1767766952966368*win[0]*fin[12]*subvol); 
  atomicAdd(&out[9], 0.1767766952966368*win[0]*fin[14]*subvol); 
  atomicAdd(&out[10], 0.1767766952966368*win[0]*fin[15]*subvol); 
  atomicAdd(&out[11], 0.1767766952966368*win[0]*fin[16]*subvol); 
  atomicAdd(&out[12], 0.1767766952966368*win[0]*fin[18]*subvol); 
  atomicAdd(&out[13], 0.1767766952966368*win[0]*fin[19]*subvol); 
  atomicAdd(&out[14], 0.1767766952966368*win[0]*fin[20]*subvol); 
  atomicAdd(&out[15], 0.1767766952966368*win[0]*fin[21]*subvol); 
  atomicAdd(&out[16], 0.1767766952966368*win[0]*fin[22]*subvol); 
  atomicAdd(&out[17], 0.1767766952966368*win[0]*fin[23]*subvol); 
  atomicAdd(&out[18], 0.1767766952966368*win[0]*fin[24]*subvol); 
  atomicAdd(&out[19], 0.1767766952966368*win[0]*fin[28]*subvol); 
  atomicAdd(&out[20], 0.1767766952966368*win[0]*fin[31]*subvol); 
  atomicAdd(&out[21], 0.1767766952966368*win[0]*fin[33]*subvol); 
  atomicAdd(&out[22], 0.1767766952966368*win[0]*fin[34]*subvol); 
  atomicAdd(&out[23], 0.1767766952966368*win[0]*fin[37]*subvol); 
  atomicAdd(&out[24], 0.1767766952966368*win[0]*fin[39]*subvol); 
  atomicAdd(&out[25], 0.1767766952966368*win[0]*fin[40]*subvol); 
  atomicAdd(&out[26], 0.1767766952966368*win[0]*fin[42]*subvol); 
  atomicAdd(&out[27], 0.1767766952966368*win[0]*fin[43]*subvol); 
  atomicAdd(&out[28], 0.1767766952966368*win[0]*fin[44]*subvol); 
  atomicAdd(&out[29], 0.1767766952966368*win[0]*fin[46]*subvol); 
  atomicAdd(&out[30], 0.1767766952966368*win[0]*fin[47]*subvol); 
  atomicAdd(&out[31], 0.1767766952966368*win[0]*fin[48]*subvol); 
  atomicAdd(&out[32], 0.1767766952966368*win[0]*fin[50]*subvol); 
  atomicAdd(&out[33], 0.1767766952966368*win[0]*fin[51]*subvol); 
  atomicAdd(&out[34], 0.1767766952966368*win[0]*fin[52]*subvol); 
  atomicAdd(&out[35], 0.1767766952966368*win[0]*fin[53]*subvol); 
  atomicAdd(&out[36], 0.1767766952966368*win[0]*fin[54]*subvol); 
  atomicAdd(&out[37], 0.1767766952966368*win[0]*fin[56]*subvol); 
  atomicAdd(&out[38], 0.1767766952966368*win[0]*fin[57]*subvol); 
  atomicAdd(&out[39], 0.1767766952966368*win[0]*fin[58]*subvol); 
  atomicAdd(&out[40], 0.1767766952966368*win[0]*fin[59]*subvol); 
  atomicAdd(&out[41], 0.1767766952966368*win[0]*fin[64]*subvol); 
  atomicAdd(&out[42], 0.1767766952966368*win[0]*fin[68]*subvol); 
  atomicAdd(&out[43], 0.1767766952966368*win[0]*fin[71]*subvol); 
  atomicAdd(&out[44], 0.1767766952966368*win[0]*fin[73]*subvol); 
  atomicAdd(&out[45], 0.1767766952966368*win[0]*fin[74]*subvol); 
  atomicAdd(&out[46], 0.1767766952966368*win[0]*fin[77]*subvol); 
  atomicAdd(&out[47], 0.1767766952966368*win[0]*fin[79]*subvol); 
  atomicAdd(&out[48], 0.1767766952966368*win[0]*fin[80]*subvol); 
  atomicAdd(&out[49], 0.1767766952966368*win[0]*fin[83]*subvol); 
  atomicAdd(&out[50], 0.1767766952966368*win[0]*fin[85]*subvol); 
  atomicAdd(&out[51], 0.1767766952966368*win[0]*fin[86]*subvol); 
  atomicAdd(&out[52], 0.1767766952966368*win[0]*fin[88]*subvol); 
  atomicAdd(&out[53], 0.1767766952966368*win[0]*fin[89]*subvol); 
  atomicAdd(&out[54], 0.1767766952966368*win[0]*fin[90]*subvol); 
  atomicAdd(&out[55], 0.1767766952966368*win[0]*fin[92]*subvol); 
  atomicAdd(&out[56], 0.1767766952966368*win[0]*fin[93]*subvol); 
  atomicAdd(&out[57], 0.1767766952966368*win[0]*fin[94]*subvol); 
  atomicAdd(&out[58], 0.1767766952966368*win[0]*fin[97]*subvol); 
  atomicAdd(&out[59], 0.1767766952966368*win[0]*fin[99]*subvol); 
  atomicAdd(&out[60], 0.1767766952966368*win[0]*fin[100]*subvol); 
  atomicAdd(&out[61], 0.1767766952966368*win[0]*fin[102]*subvol); 
  atomicAdd(&out[62], 0.1767766952966368*win[0]*fin[103]*subvol); 
  atomicAdd(&out[63], 0.1767766952966368*win[0]*fin[104]*subvol); 
  atomicAdd(&out[64], 0.1767766952966368*win[0]*fin[105]*subvol); 
  atomicAdd(&out[65], 0.1767766952966368*win[0]*fin[106]*subvol); 
  atomicAdd(&out[66], 0.1767766952966368*win[0]*fin[107]*subvol); 
  atomicAdd(&out[67], 0.1767766952966368*win[0]*fin[113]*subvol); 
  atomicAdd(&out[68], 0.1767766952966368*win[0]*fin[117]*subvol); 
  atomicAdd(&out[69], 0.1767766952966368*win[0]*fin[121]*subvol); 
  atomicAdd(&out[70], 0.1767766952966368*win[0]*fin[124]*subvol); 
  atomicAdd(&out[71], 0.1767766952966368*win[0]*fin[126]*subvol); 
  atomicAdd(&out[72], 0.1767766952966368*win[0]*fin[127]*subvol); 
  atomicAdd(&out[73], 0.1767766952966368*win[0]*fin[130]*subvol); 
  atomicAdd(&out[74], 0.1767766952966368*win[0]*fin[132]*subvol); 
  atomicAdd(&out[75], 0.1767766952966368*win[0]*fin[133]*subvol); 
  atomicAdd(&out[76], 0.1767766952966368*win[0]*fin[137]*subvol); 
  atomicAdd(&out[77], 0.1767766952966368*win[0]*fin[140]*subvol); 
  atomicAdd(&out[78], 0.1767766952966368*win[0]*fin[142]*subvol); 
  atomicAdd(&out[79], 0.1767766952966368*win[0]*fin[143]*subvol); 
  atomicAdd(&out[80], 0.1767766952966368*win[0]*fin[145]*subvol); 
  atomicAdd(&out[81], 0.1767766952966368*win[0]*fin[146]*subvol); 
  atomicAdd(&out[82], 0.1767766952966368*win[0]*fin[147]*subvol); 
  atomicAdd(&out[83], 0.1767766952966368*win[0]*fin[149]*subvol); 
  atomicAdd(&out[84], 0.1767766952966368*win[0]*fin[150]*subvol); 
  atomicAdd(&out[85], 0.1767766952966368*win[0]*fin[151]*subvol); 
  atomicAdd(&out[86], 0.1767766952966368*win[0]*fin[153]*subvol); 
  atomicAdd(&out[87], 0.1767766952966368*win[0]*fin[154]*subvol); 
  atomicAdd(&out[88], 0.1767766952966368*win[0]*fin[155]*subvol); 
  atomicAdd(&out[89], 0.1767766952966368*win[0]*fin[162]*subvol); 
  atomicAdd(&out[90], 0.1767766952966368*win[0]*fin[166]*subvol); 
  atomicAdd(&out[91], 0.1767766952966368*win[0]*fin[171]*subvol); 
  atomicAdd(&out[92], 0.1767766952966368*win[0]*fin[174]*subvol); 
  atomicAdd(&out[93], 0.1767766952966368*win[0]*fin[176]*subvol); 
  atomicAdd(&out[94], 0.1767766952966368*win[0]*fin[177]*subvol); 
  atomicAdd(&out[95], 0.1767766952966368*win[0]*fin[180]*subvol); 
  atomicAdd(&out[96], 0.1767766952966368*win[0]*fin[182]*subvol); 
  atomicAdd(&out[97], 0.1767766952966368*win[0]*fin[183]*subvol); 
  atomicAdd(&out[98], 0.1767766952966368*win[0]*fin[186]*subvol); 
  atomicAdd(&out[99], 0.1767766952966368*win[0]*fin[188]*subvol); 
  atomicAdd(&out[100], 0.1767766952966368*win[0]*fin[189]*subvol); 
  atomicAdd(&out[101], 0.1767766952966368*win[0]*fin[190]*subvol); 
  atomicAdd(&out[102], 0.1767766952966368*win[0]*fin[197]*subvol); 
  atomicAdd(&out[103], 0.1767766952966368*win[0]*fin[201]*subvol); 
  atomicAdd(&out[104], 0.1767766952966368*win[0]*fin[205]*subvol); 
  atomicAdd(&out[105], 0.1767766952966368*win[0]*fin[207]*subvol); 
  atomicAdd(&out[106], 0.1767766952966368*win[0]*fin[208]*subvol); 
  atomicAdd(&out[107], 0.1767766952966368*win[0]*fin[214]*subvol); 
#else
  out[0] += 0.1767766952966368*fin[0]*win[0]*subvol; 
  out[1] += 0.1767766952966368*win[0]*fin[2]*subvol; 
  out[2] += 0.1767766952966368*win[0]*fin[3]*subvol; 
  out[3] += 0.1767766952966368*win[0]*fin[4]*subvol; 
  out[4] += 0.1767766952966368*win[0]*fin[5]*subvol; 
  out[5] += 0.1767766952966368*win[0]*fin[6]*subvol; 
  out[6] += 0.1767766952966368*win[0]*fin[9]*subvol; 
  out[7] += 0.1767766952966368*win[0]*fin[11]*subvol; 
  out[8] += 0.1767766952966368*win[0]*fin[12]*subvol; 
  out[9] += 0.1767766952966368*win[0]*fin[14]*subvol; 
  out[10] += 0.1767766952966368*win[0]*fin[15]*subvol; 
  out[11] += 0.1767766952966368*win[0]*fin[16]*subvol; 
  out[12] += 0.1767766952966368*win[0]*fin[18]*subvol; 
  out[13] += 0.1767766952966368*win[0]*fin[19]*subvol; 
  out[14] += 0.1767766952966368*win[0]*fin[20]*subvol; 
  out[15] += 0.1767766952966368*win[0]*fin[21]*subvol; 
  out[16] += 0.1767766952966368*win[0]*fin[22]*subvol; 
  out[17] += 0.1767766952966368*win[0]*fin[23]*subvol; 
  out[18] += 0.1767766952966368*win[0]*fin[24]*subvol; 
  out[19] += 0.1767766952966368*win[0]*fin[28]*subvol; 
  out[20] += 0.1767766952966368*win[0]*fin[31]*subvol; 
  out[21] += 0.1767766952966368*win[0]*fin[33]*subvol; 
  out[22] += 0.1767766952966368*win[0]*fin[34]*subvol; 
  out[23] += 0.1767766952966368*win[0]*fin[37]*subvol; 
  out[24] += 0.1767766952966368*win[0]*fin[39]*subvol; 
  out[25] += 0.1767766952966368*win[0]*fin[40]*subvol; 
  out[26] += 0.1767766952966368*win[0]*fin[42]*subvol; 
  out[27] += 0.1767766952966368*win[0]*fin[43]*subvol; 
  out[28] += 0.1767766952966368*win[0]*fin[44]*subvol; 
  out[29] += 0.1767766952966368*win[0]*fin[46]*subvol; 
  out[30] += 0.1767766952966368*win[0]*fin[47]*subvol; 
  out[31] += 0.1767766952966368*win[0]*fin[48]*subvol; 
  out[32] += 0.1767766952966368*win[0]*fin[50]*subvol; 
  out[33] += 0.1767766952966368*win[0]*fin[51]*subvol; 
  out[34] += 0.1767766952966368*win[0]*fin[52]*subvol; 
  out[35] += 0.1767766952966368*win[0]*fin[53]*subvol; 
  out[36] += 0.1767766952966368*win[0]*fin[54]*subvol; 
  out[37] += 0.1767766952966368*win[0]*fin[56]*subvol; 
  out[38] += 0.1767766952966368*win[0]*fin[57]*subvol; 
  out[39] += 0.1767766952966368*win[0]*fin[58]*subvol; 
  out[40] += 0.1767766952966368*win[0]*fin[59]*subvol; 
  out[41] += 0.1767766952966368*win[0]*fin[64]*subvol; 
  out[42] += 0.1767766952966368*win[0]*fin[68]*subvol; 
  out[43] += 0.1767766952966368*win[0]*fin[71]*subvol; 
  out[44] += 0.1767766952966368*win[0]*fin[73]*subvol; 
  out[45] += 0.1767766952966368*win[0]*fin[74]*subvol; 
  out[46] += 0.1767766952966368*win[0]*fin[77]*subvol; 
  out[47] += 0.1767766952966368*win[0]*fin[79]*subvol; 
  out[48] += 0.1767766952966368*win[0]*fin[80]*subvol; 
  out[49] += 0.1767766952966368*win[0]*fin[83]*subvol; 
  out[50] += 0.1767766952966368*win[0]*fin[85]*subvol; 
  out[51] += 0.1767766952966368*win[0]*fin[86]*subvol; 
  out[52] += 0.1767766952966368*win[0]*fin[88]*subvol; 
  out[53] += 0.1767766952966368*win[0]*fin[89]*subvol; 
  out[54] += 0.1767766952966368*win[0]*fin[90]*subvol; 
  out[55] += 0.1767766952966368*win[0]*fin[92]*subvol; 
  out[56] += 0.1767766952966368*win[0]*fin[93]*subvol; 
  out[57] += 0.1767766952966368*win[0]*fin[94]*subvol; 
  out[58] += 0.1767766952966368*win[0]*fin[97]*subvol; 
  out[59] += 0.1767766952966368*win[0]*fin[99]*subvol; 
  out[60] += 0.1767766952966368*win[0]*fin[100]*subvol; 
  out[61] += 0.1767766952966368*win[0]*fin[102]*subvol; 
  out[62] += 0.1767766952966368*win[0]*fin[103]*subvol; 
  out[63] += 0.1767766952966368*win[0]*fin[104]*subvol; 
  out[64] += 0.1767766952966368*win[0]*fin[105]*subvol; 
  out[65] += 0.1767766952966368*win[0]*fin[106]*subvol; 
  out[66] += 0.1767766952966368*win[0]*fin[107]*subvol; 
  out[67] += 0.1767766952966368*win[0]*fin[113]*subvol; 
  out[68] += 0.1767766952966368*win[0]*fin[117]*subvol; 
  out[69] += 0.1767766952966368*win[0]*fin[121]*subvol; 
  out[70] += 0.1767766952966368*win[0]*fin[124]*subvol; 
  out[71] += 0.1767766952966368*win[0]*fin[126]*subvol; 
  out[72] += 0.1767766952966368*win[0]*fin[127]*subvol; 
  out[73] += 0.1767766952966368*win[0]*fin[130]*subvol; 
  out[74] += 0.1767766952966368*win[0]*fin[132]*subvol; 
  out[75] += 0.1767766952966368*win[0]*fin[133]*subvol; 
  out[76] += 0.1767766952966368*win[0]*fin[137]*subvol; 
  out[77] += 0.1767766952966368*win[0]*fin[140]*subvol; 
  out[78] += 0.1767766952966368*win[0]*fin[142]*subvol; 
  out[79] += 0.1767766952966368*win[0]*fin[143]*subvol; 
  out[80] += 0.1767766952966368*win[0]*fin[145]*subvol; 
  out[81] += 0.1767766952966368*win[0]*fin[146]*subvol; 
  out[82] += 0.1767766952966368*win[0]*fin[147]*subvol; 
  out[83] += 0.1767766952966368*win[0]*fin[149]*subvol; 
  out[84] += 0.1767766952966368*win[0]*fin[150]*subvol; 
  out[85] += 0.1767766952966368*win[0]*fin[151]*subvol; 
  out[86] += 0.1767766952966368*win[0]*fin[153]*subvol; 
  out[87] += 0.1767766952966368*win[0]*fin[154]*subvol; 
  out[88] += 0.1767766952966368*win[0]*fin[155]*subvol; 
  out[89] += 0.1767766952966368*win[0]*fin[162]*subvol; 
  out[90] += 0.1767766952966368*win[0]*fin[166]*subvol; 
  out[91] += 0.1767766952966368*win[0]*fin[171]*subvol; 
  out[92] += 0.1767766952966368*win[0]*fin[174]*subvol; 
  out[93] += 0.1767766952966368*win[0]*fin[176]*subvol; 
  out[94] += 0.1767766952966368*win[0]*fin[177]*subvol; 
  out[95] += 0.1767766952966368*win[0]*fin[180]*subvol; 
  out[96] += 0.1767766952966368*win[0]*fin[182]*subvol; 
  out[97] += 0.1767766952966368*win[0]*fin[183]*subvol; 
  out[98] += 0.1767766952966368*win[0]*fin[186]*subvol; 
  out[99] += 0.1767766952966368*win[0]*fin[188]*subvol; 
  out[100] += 0.1767766952966368*win[0]*fin[189]*subvol; 
  out[101] += 0.1767766952966368*win[0]*fin[190]*subvol; 
  out[102] += 0.1767766952966368*win[0]*fin[197]*subvol; 
  out[103] += 0.1767766952966368*win[0]*fin[201]*subvol; 
  out[104] += 0.1767766952966368*win[0]*fin[205]*subvol; 
  out[105] += 0.1767766952966368*win[0]*fin[207]*subvol; 
  out[106] += 0.1767766952966368*win[0]*fin[208]*subvol; 
  out[107] += 0.1767766952966368*win[0]*fin[214]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.25*fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], 0.25*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.25*win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], 0.25*win[0]*fin[4]*subvol); 
  atomicAdd(&out[5], 0.25*win[0]*fin[5]*subvol); 
  atomicAdd(&out[6], 0.25*win[0]*fin[6]*subvol); 
  atomicAdd(&out[9], 0.25*win[0]*fin[9]*subvol); 
  atomicAdd(&out[11], 0.25*win[0]*fin[11]*subvol); 
  atomicAdd(&out[12], 0.25*win[0]*fin[12]*subvol); 
  atomicAdd(&out[14], 0.25*win[0]*fin[14]*subvol); 
  atomicAdd(&out[15], 0.25*win[0]*fin[15]*subvol); 
  atomicAdd(&out[16], 0.25*win[0]*fin[16]*subvol); 
  atomicAdd(&out[18], 0.25*win[0]*fin[18]*subvol); 
  atomicAdd(&out[19], 0.25*win[0]*fin[19]*subvol); 
  atomicAdd(&out[20], 0.25*win[0]*fin[20]*subvol); 
  atomicAdd(&out[21], 0.25*win[0]*fin[21]*subvol); 
  atomicAdd(&out[22], 0.25*win[0]*fin[22]*subvol); 
  atomicAdd(&out[23], 0.25*win[0]*fin[23]*subvol); 
  atomicAdd(&out[24], 0.25*win[0]*fin[24]*subvol); 
  atomicAdd(&out[28], 0.25*win[0]*fin[28]*subvol); 
  atomicAdd(&out[31], 0.25*win[0]*fin[31]*subvol); 
  atomicAdd(&out[33], 0.25*win[0]*fin[33]*subvol); 
  atomicAdd(&out[34], 0.25*win[0]*fin[34]*subvol); 
  atomicAdd(&out[37], 0.25*win[0]*fin[37]*subvol); 
  atomicAdd(&out[39], 0.25*win[0]*fin[39]*subvol); 
  atomicAdd(&out[40], 0.25*win[0]*fin[40]*subvol); 
  atomicAdd(&out[42], 0.25*win[0]*fin[42]*subvol); 
  atomicAdd(&out[43], 0.25*win[0]*fin[43]*subvol); 
  atomicAdd(&out[44], 0.25*win[0]*fin[44]*subvol); 
  atomicAdd(&out[46], 0.25*win[0]*fin[46]*subvol); 
  atomicAdd(&out[47], 0.25*win[0]*fin[47]*subvol); 
  atomicAdd(&out[48], 0.25*win[0]*fin[48]*subvol); 
  atomicAdd(&out[50], 0.25*win[0]*fin[50]*subvol); 
  atomicAdd(&out[51], 0.25*win[0]*fin[51]*subvol); 
  atomicAdd(&out[52], 0.25*win[0]*fin[52]*subvol); 
  atomicAdd(&out[53], 0.25*win[0]*fin[53]*subvol); 
  atomicAdd(&out[54], 0.25*win[0]*fin[54]*subvol); 
  atomicAdd(&out[56], 0.25*win[0]*fin[56]*subvol); 
  atomicAdd(&out[57], 0.25*win[0]*fin[57]*subvol); 
  atomicAdd(&out[58], 0.25*win[0]*fin[58]*subvol); 
  atomicAdd(&out[59], 0.25*win[0]*fin[59]*subvol); 
  atomicAdd(&out[64], 0.25*win[0]*fin[64]*subvol); 
  atomicAdd(&out[68], 0.25*win[0]*fin[68]*subvol); 
  atomicAdd(&out[71], 0.25*win[0]*fin[71]*subvol); 
  atomicAdd(&out[73], 0.25*win[0]*fin[73]*subvol); 
  atomicAdd(&out[74], 0.25*win[0]*fin[74]*subvol); 
  atomicAdd(&out[77], 0.25*win[0]*fin[77]*subvol); 
  atomicAdd(&out[79], 0.25*win[0]*fin[79]*subvol); 
  atomicAdd(&out[80], 0.25*win[0]*fin[80]*subvol); 
  atomicAdd(&out[83], 0.25*win[0]*fin[83]*subvol); 
  atomicAdd(&out[85], 0.25*win[0]*fin[85]*subvol); 
  atomicAdd(&out[86], 0.25*win[0]*fin[86]*subvol); 
  atomicAdd(&out[88], 0.25*win[0]*fin[88]*subvol); 
  atomicAdd(&out[89], 0.25*win[0]*fin[89]*subvol); 
  atomicAdd(&out[90], 0.25*win[0]*fin[90]*subvol); 
  atomicAdd(&out[92], 0.25*win[0]*fin[92]*subvol); 
  atomicAdd(&out[93], 0.25*win[0]*fin[93]*subvol); 
  atomicAdd(&out[94], 0.25*win[0]*fin[94]*subvol); 
  atomicAdd(&out[97], 0.25*win[0]*fin[97]*subvol); 
  atomicAdd(&out[99], 0.25*win[0]*fin[99]*subvol); 
  atomicAdd(&out[100], 0.25*win[0]*fin[100]*subvol); 
  atomicAdd(&out[102], 0.25*win[0]*fin[102]*subvol); 
  atomicAdd(&out[103], 0.25*win[0]*fin[103]*subvol); 
  atomicAdd(&out[104], 0.25*win[0]*fin[104]*subvol); 
  atomicAdd(&out[105], 0.25*win[0]*fin[105]*subvol); 
  atomicAdd(&out[106], 0.25*win[0]*fin[106]*subvol); 
  atomicAdd(&out[107], 0.25*win[0]*fin[107]*subvol); 
  atomicAdd(&out[113], 0.25*win[0]*fin[113]*subvol); 
  atomicAdd(&out[117], 0.25*win[0]*fin[117]*subvol); 
  atomicAdd(&out[121], 0.25*win[0]*fin[121]*subvol); 
  atomicAdd(&out[124], 0.25*win[0]*fin[124]*subvol); 
  atomicAdd(&out[126], 0.25*win[0]*fin[126]*subvol); 
  atomicAdd(&out[127], 0.25*win[0]*fin[127]*subvol); 
  atomicAdd(&out[130], 0.25*win[0]*fin[130]*subvol); 
  atomicAdd(&out[132], 0.25*win[0]*fin[132]*subvol); 
  atomicAdd(&out[133], 0.25*win[0]*fin[133]*subvol); 
  atomicAdd(&out[137], 0.25*win[0]*fin[137]*subvol); 
  atomicAdd(&out[140], 0.25*win[0]*fin[140]*subvol); 
  atomicAdd(&out[142], 0.25*win[0]*fin[142]*subvol); 
  atomicAdd(&out[143], 0.25*win[0]*fin[143]*subvol); 
  atomicAdd(&out[145], 0.25*win[0]*fin[145]*subvol); 
  atomicAdd(&out[146], 0.25*win[0]*fin[146]*subvol); 
  atomicAdd(&out[147], 0.25*win[0]*fin[147]*subvol); 
  atomicAdd(&out[149], 0.25*win[0]*fin[149]*subvol); 
  atomicAdd(&out[150], 0.25*win[0]*fin[150]*subvol); 
  atomicAdd(&out[151], 0.25*win[0]*fin[151]*subvol); 
  atomicAdd(&out[153], 0.25*win[0]*fin[153]*subvol); 
  atomicAdd(&out[154], 0.25*win[0]*fin[154]*subvol); 
  atomicAdd(&out[155], 0.25*win[0]*fin[155]*subvol); 
  atomicAdd(&out[162], 0.25*win[0]*fin[162]*subvol); 
  atomicAdd(&out[166], 0.25*win[0]*fin[166]*subvol); 
  atomicAdd(&out[171], 0.25*win[0]*fin[171]*subvol); 
  atomicAdd(&out[174], 0.25*win[0]*fin[174]*subvol); 
  atomicAdd(&out[176], 0.25*win[0]*fin[176]*subvol); 
  atomicAdd(&out[177], 0.25*win[0]*fin[177]*subvol); 
  atomicAdd(&out[180], 0.25*win[0]*fin[180]*subvol); 
  atomicAdd(&out[182], 0.25*win[0]*fin[182]*subvol); 
  atomicAdd(&out[183], 0.25*win[0]*fin[183]*subvol); 
  atomicAdd(&out[186], 0.25*win[0]*fin[186]*subvol); 
  atomicAdd(&out[188], 0.25*win[0]*fin[188]*subvol); 
  atomicAdd(&out[189], 0.25*win[0]*fin[189]*subvol); 
  atomicAdd(&out[190], 0.25*win[0]*fin[190]*subvol); 
  atomicAdd(&out[197], 0.25*win[0]*fin[197]*subvol); 
  atomicAdd(&out[201], 0.25*win[0]*fin[201]*subvol); 
  atomicAdd(&out[205], 0.25*win[0]*fin[205]*subvol); 
  atomicAdd(&out[207], 0.25*win[0]*fin[207]*subvol); 
  atomicAdd(&out[208], 0.25*win[0]*fin[208]*subvol); 
  atomicAdd(&out[214], 0.25*win[0]*fin[214]*subvol); 
#else
  out[0] += 0.25*fin[0]*win[0]*subvol; 
  out[2] += 0.25*win[0]*fin[2]*subvol; 
  out[3] += 0.25*win[0]*fin[3]*subvol; 
  out[4] += 0.25*win[0]*fin[4]*subvol; 
  out[5] += 0.25*win[0]*fin[5]*subvol; 
  out[6] += 0.25*win[0]*fin[6]*subvol; 
  out[9] += 0.25*win[0]*fin[9]*subvol; 
  out[11] += 0.25*win[0]*fin[11]*subvol; 
  out[12] += 0.25*win[0]*fin[12]*subvol; 
  out[14] += 0.25*win[0]*fin[14]*subvol; 
  out[15] += 0.25*win[0]*fin[15]*subvol; 
  out[16] += 0.25*win[0]*fin[16]*subvol; 
  out[18] += 0.25*win[0]*fin[18]*subvol; 
  out[19] += 0.25*win[0]*fin[19]*subvol; 
  out[20] += 0.25*win[0]*fin[20]*subvol; 
  out[21] += 0.25*win[0]*fin[21]*subvol; 
  out[22] += 0.25*win[0]*fin[22]*subvol; 
  out[23] += 0.25*win[0]*fin[23]*subvol; 
  out[24] += 0.25*win[0]*fin[24]*subvol; 
  out[28] += 0.25*win[0]*fin[28]*subvol; 
  out[31] += 0.25*win[0]*fin[31]*subvol; 
  out[33] += 0.25*win[0]*fin[33]*subvol; 
  out[34] += 0.25*win[0]*fin[34]*subvol; 
  out[37] += 0.25*win[0]*fin[37]*subvol; 
  out[39] += 0.25*win[0]*fin[39]*subvol; 
  out[40] += 0.25*win[0]*fin[40]*subvol; 
  out[42] += 0.25*win[0]*fin[42]*subvol; 
  out[43] += 0.25*win[0]*fin[43]*subvol; 
  out[44] += 0.25*win[0]*fin[44]*subvol; 
  out[46] += 0.25*win[0]*fin[46]*subvol; 
  out[47] += 0.25*win[0]*fin[47]*subvol; 
  out[48] += 0.25*win[0]*fin[48]*subvol; 
  out[50] += 0.25*win[0]*fin[50]*subvol; 
  out[51] += 0.25*win[0]*fin[51]*subvol; 
  out[52] += 0.25*win[0]*fin[52]*subvol; 
  out[53] += 0.25*win[0]*fin[53]*subvol; 
  out[54] += 0.25*win[0]*fin[54]*subvol; 
  out[56] += 0.25*win[0]*fin[56]*subvol; 
  out[57] += 0.25*win[0]*fin[57]*subvol; 
  out[58] += 0.25*win[0]*fin[58]*subvol; 
  out[59] += 0.25*win[0]*fin[59]*subvol; 
  out[64] += 0.25*win[0]*fin[64]*subvol; 
  out[68] += 0.25*win[0]*fin[68]*subvol; 
  out[71] += 0.25*win[0]*fin[71]*subvol; 
  out[73] += 0.25*win[0]*fin[73]*subvol; 
  out[74] += 0.25*win[0]*fin[74]*subvol; 
  out[77] += 0.25*win[0]*fin[77]*subvol; 
  out[79] += 0.25*win[0]*fin[79]*subvol; 
  out[80] += 0.25*win[0]*fin[80]*subvol; 
  out[83] += 0.25*win[0]*fin[83]*subvol; 
  out[85] += 0.25*win[0]*fin[85]*subvol; 
  out[86] += 0.25*win[0]*fin[86]*subvol; 
  out[88] += 0.25*win[0]*fin[88]*subvol; 
  out[89] += 0.25*win[0]*fin[89]*subvol; 
  out[90] += 0.25*win[0]*fin[90]*subvol; 
  out[92] += 0.25*win[0]*fin[92]*subvol; 
  out[93] += 0.25*win[0]*fin[93]*subvol; 
  out[94] += 0.25*win[0]*fin[94]*subvol; 
  out[97] += 0.25*win[0]*fin[97]*subvol; 
  out[99] += 0.25*win[0]*fin[99]*subvol; 
  out[100] += 0.25*win[0]*fin[100]*subvol; 
  out[102] += 0.25*win[0]*fin[102]*subvol; 
  out[103] += 0.25*win[0]*fin[103]*subvol; 
  out[104] += 0.25*win[0]*fin[104]*subvol; 
  out[105] += 0.25*win[0]*fin[105]*subvol; 
  out[106] += 0.25*win[0]*fin[106]*subvol; 
  out[107] += 0.25*win[0]*fin[107]*subvol; 
  out[113] += 0.25*win[0]*fin[113]*subvol; 
  out[117] += 0.25*win[0]*fin[117]*subvol; 
  out[121] += 0.25*win[0]*fin[121]*subvol; 
  out[124] += 0.25*win[0]*fin[124]*subvol; 
  out[126] += 0.25*win[0]*fin[126]*subvol; 
  out[127] += 0.25*win[0]*fin[127]*subvol; 
  out[130] += 0.25*win[0]*fin[130]*subvol; 
  out[132] += 0.25*win[0]*fin[132]*subvol; 
  out[133] += 0.25*win[0]*fin[133]*subvol; 
  out[137] += 0.25*win[0]*fin[137]*subvol; 
  out[140] += 0.25*win[0]*fin[140]*subvol; 
  out[142] += 0.25*win[0]*fin[142]*subvol; 
  out[143] += 0.25*win[0]*fin[143]*subvol; 
  out[145] += 0.25*win[0]*fin[145]*subvol; 
  out[146] += 0.25*win[0]*fin[146]*subvol; 
  out[147] += 0.25*win[0]*fin[147]*subvol; 
  out[149] += 0.25*win[0]*fin[149]*subvol; 
  out[150] += 0.25*win[0]*fin[150]*subvol; 
  out[151] += 0.25*win[0]*fin[151]*subvol; 
  out[153] += 0.25*win[0]*fin[153]*subvol; 
  out[154] += 0.25*win[0]*fin[154]*subvol; 
  out[155] += 0.25*win[0]*fin[155]*subvol; 
  out[162] += 0.25*win[0]*fin[162]*subvol; 
  out[166] += 0.25*win[0]*fin[166]*subvol; 
  out[171] += 0.25*win[0]*fin[171]*subvol; 
  out[174] += 0.25*win[0]*fin[174]*subvol; 
  out[176] += 0.25*win[0]*fin[176]*subvol; 
  out[177] += 0.25*win[0]*fin[177]*subvol; 
  out[180] += 0.25*win[0]*fin[180]*subvol; 
  out[182] += 0.25*win[0]*fin[182]*subvol; 
  out[183] += 0.25*win[0]*fin[183]*subvol; 
  out[186] += 0.25*win[0]*fin[186]*subvol; 
  out[188] += 0.25*win[0]*fin[188]*subvol; 
  out[189] += 0.25*win[0]*fin[189]*subvol; 
  out[190] += 0.25*win[0]*fin[190]*subvol; 
  out[197] += 0.25*win[0]*fin[197]*subvol; 
  out[201] += 0.25*win[0]*fin[201]*subvol; 
  out[205] += 0.25*win[0]*fin[205]*subvol; 
  out[207] += 0.25*win[0]*fin[207]*subvol; 
  out[208] += 0.25*win[0]*fin[208]*subvol; 
  out[214] += 0.25*win[0]*fin[214]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.1767766952966368*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.1767766952966368*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.1767766952966368*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.1767766952966368*win[0]*fin[4]*subvol); 
  atomicAdd(&out[4], 0.1767766952966368*win[0]*fin[5]*subvol); 
  atomicAdd(&out[5], 0.1767766952966368*win[0]*fin[6]*subvol); 
  atomicAdd(&out[6], 0.1767766952966368*win[0]*fin[8]*subvol); 
  atomicAdd(&out[7], 0.1767766952966368*win[0]*fin[10]*subvol); 
  atomicAdd(&out[8], 0.1767766952966368*win[0]*fin[12]*subvol); 
  atomicAdd(&out[9], 0.1767766952966368*win[0]*fin[13]*subvol); 
  atomicAdd(&out[10], 0.1767766952966368*win[0]*fin[15]*subvol); 
  atomicAdd(&out[11], 0.1767766952966368*win[0]*fin[16]*subvol); 
  atomicAdd(&out[12], 0.1767766952966368*win[0]*fin[17]*subvol); 
  atomicAdd(&out[13], 0.1767766952966368*win[0]*fin[19]*subvol); 
  atomicAdd(&out[14], 0.1767766952966368*win[0]*fin[20]*subvol); 
  atomicAdd(&out[15], 0.1767766952966368*win[0]*fin[21]*subvol); 
  atomicAdd(&out[16], 0.1767766952966368*win[0]*fin[22]*subvol); 
  atomicAdd(&out[17], 0.1767766952966368*win[0]*fin[23]*subvol); 
  atomicAdd(&out[18], 0.1767766952966368*win[0]*fin[24]*subvol); 
  atomicAdd(&out[19], 0.1767766952966368*win[0]*fin[27]*subvol); 
  atomicAdd(&out[20], 0.1767766952966368*win[0]*fin[30]*subvol); 
  atomicAdd(&out[21], 0.1767766952966368*win[0]*fin[32]*subvol); 
  atomicAdd(&out[22], 0.1767766952966368*win[0]*fin[34]*subvol); 
  atomicAdd(&out[23], 0.1767766952966368*win[0]*fin[36]*subvol); 
  atomicAdd(&out[24], 0.1767766952966368*win[0]*fin[38]*subvol); 
  atomicAdd(&out[25], 0.1767766952966368*win[0]*fin[40]*subvol); 
  atomicAdd(&out[26], 0.1767766952966368*win[0]*fin[41]*subvol); 
  atomicAdd(&out[27], 0.1767766952966368*win[0]*fin[43]*subvol); 
  atomicAdd(&out[28], 0.1767766952966368*win[0]*fin[44]*subvol); 
  atomicAdd(&out[29], 0.1767766952966368*win[0]*fin[45]*subvol); 
  atomicAdd(&out[30], 0.1767766952966368*win[0]*fin[47]*subvol); 
  atomicAdd(&out[31], 0.1767766952966368*win[0]*fin[48]*subvol); 
  atomicAdd(&out[32], 0.1767766952966368*win[0]*fin[49]*subvol); 
  atomicAdd(&out[33], 0.1767766952966368*win[0]*fin[51]*subvol); 
  atomicAdd(&out[34], 0.1767766952966368*win[0]*fin[52]*subvol); 
  atomicAdd(&out[35], 0.1767766952966368*win[0]*fin[53]*subvol); 
  atomicAdd(&out[36], 0.1767766952966368*win[0]*fin[54]*subvol); 
  atomicAdd(&out[37], 0.1767766952966368*win[0]*fin[55]*subvol); 
  atomicAdd(&out[38], 0.1767766952966368*win[0]*fin[57]*subvol); 
  atomicAdd(&out[39], 0.1767766952966368*win[0]*fin[58]*subvol); 
  atomicAdd(&out[40], 0.1767766952966368*win[0]*fin[59]*subvol); 
  atomicAdd(&out[41], 0.1767766952966368*win[0]*fin[63]*subvol); 
  atomicAdd(&out[42], 0.1767766952966368*win[0]*fin[67]*subvol); 
  atomicAdd(&out[43], 0.1767766952966368*win[0]*fin[70]*subvol); 
  atomicAdd(&out[44], 0.1767766952966368*win[0]*fin[72]*subvol); 
  atomicAdd(&out[45], 0.1767766952966368*win[0]*fin[74]*subvol); 
  atomicAdd(&out[46], 0.1767766952966368*win[0]*fin[76]*subvol); 
  atomicAdd(&out[47], 0.1767766952966368*win[0]*fin[78]*subvol); 
  atomicAdd(&out[48], 0.1767766952966368*win[0]*fin[80]*subvol); 
  atomicAdd(&out[49], 0.1767766952966368*win[0]*fin[82]*subvol); 
  atomicAdd(&out[50], 0.1767766952966368*win[0]*fin[84]*subvol); 
  atomicAdd(&out[51], 0.1767766952966368*win[0]*fin[86]*subvol); 
  atomicAdd(&out[52], 0.1767766952966368*win[0]*fin[87]*subvol); 
  atomicAdd(&out[53], 0.1767766952966368*win[0]*fin[89]*subvol); 
  atomicAdd(&out[54], 0.1767766952966368*win[0]*fin[90]*subvol); 
  atomicAdd(&out[55], 0.1767766952966368*win[0]*fin[91]*subvol); 
  atomicAdd(&out[56], 0.1767766952966368*win[0]*fin[93]*subvol); 
  atomicAdd(&out[57], 0.1767766952966368*win[0]*fin[94]*subvol); 
  atomicAdd(&out[58], 0.1767766952966368*win[0]*fin[96]*subvol); 
  atomicAdd(&out[59], 0.1767766952966368*win[0]*fin[98]*subvol); 
  atomicAdd(&out[60], 0.1767766952966368*win[0]*fin[100]*subvol); 
  atomicAdd(&out[61], 0.1767766952966368*win[0]*fin[101]*subvol); 
  atomicAdd(&out[62], 0.1767766952966368*win[0]*fin[103]*subvol); 
  atomicAdd(&out[63], 0.1767766952966368*win[0]*fin[104]*subvol); 
  atomicAdd(&out[64], 0.1767766952966368*win[0]*fin[105]*subvol); 
  atomicAdd(&out[65], 0.1767766952966368*win[0]*fin[106]*subvol); 
  atomicAdd(&out[66], 0.1767766952966368*win[0]*fin[107]*subvol); 
  atomicAdd(&out[67], 0.1767766952966368*win[0]*fin[112]*subvol); 
  atomicAdd(&out[68], 0.1767766952966368*win[0]*fin[116]*subvol); 
  atomicAdd(&out[69], 0.1767766952966368*win[0]*fin[120]*subvol); 
  atomicAdd(&out[70], 0.1767766952966368*win[0]*fin[123]*subvol); 
  atomicAdd(&out[71], 0.1767766952966368*win[0]*fin[125]*subvol); 
  atomicAdd(&out[72], 0.1767766952966368*win[0]*fin[127]*subvol); 
  atomicAdd(&out[73], 0.1767766952966368*win[0]*fin[129]*subvol); 
  atomicAdd(&out[74], 0.1767766952966368*win[0]*fin[131]*subvol); 
  atomicAdd(&out[75], 0.1767766952966368*win[0]*fin[133]*subvol); 
  atomicAdd(&out[76], 0.1767766952966368*win[0]*fin[136]*subvol); 
  atomicAdd(&out[77], 0.1767766952966368*win[0]*fin[139]*subvol); 
  atomicAdd(&out[78], 0.1767766952966368*win[0]*fin[141]*subvol); 
  atomicAdd(&out[79], 0.1767766952966368*win[0]*fin[143]*subvol); 
  atomicAdd(&out[80], 0.1767766952966368*win[0]*fin[144]*subvol); 
  atomicAdd(&out[81], 0.1767766952966368*win[0]*fin[146]*subvol); 
  atomicAdd(&out[82], 0.1767766952966368*win[0]*fin[147]*subvol); 
  atomicAdd(&out[83], 0.1767766952966368*win[0]*fin[148]*subvol); 
  atomicAdd(&out[84], 0.1767766952966368*win[0]*fin[150]*subvol); 
  atomicAdd(&out[85], 0.1767766952966368*win[0]*fin[151]*subvol); 
  atomicAdd(&out[86], 0.1767766952966368*win[0]*fin[152]*subvol); 
  atomicAdd(&out[87], 0.1767766952966368*win[0]*fin[154]*subvol); 
  atomicAdd(&out[88], 0.1767766952966368*win[0]*fin[155]*subvol); 
  atomicAdd(&out[89], 0.1767766952966368*win[0]*fin[161]*subvol); 
  atomicAdd(&out[90], 0.1767766952966368*win[0]*fin[165]*subvol); 
  atomicAdd(&out[91], 0.1767766952966368*win[0]*fin[170]*subvol); 
  atomicAdd(&out[92], 0.1767766952966368*win[0]*fin[173]*subvol); 
  atomicAdd(&out[93], 0.1767766952966368*win[0]*fin[175]*subvol); 
  atomicAdd(&out[94], 0.1767766952966368*win[0]*fin[177]*subvol); 
  atomicAdd(&out[95], 0.1767766952966368*win[0]*fin[179]*subvol); 
  atomicAdd(&out[96], 0.1767766952966368*win[0]*fin[181]*subvol); 
  atomicAdd(&out[97], 0.1767766952966368*win[0]*fin[183]*subvol); 
  atomicAdd(&out[98], 0.1767766952966368*win[0]*fin[185]*subvol); 
  atomicAdd(&out[99], 0.1767766952966368*win[0]*fin[187]*subvol); 
  atomicAdd(&out[100], 0.1767766952966368*win[0]*fin[189]*subvol); 
  atomicAdd(&out[101], 0.1767766952966368*win[0]*fin[190]*subvol); 
  atomicAdd(&out[102], 0.1767766952966368*win[0]*fin[196]*subvol); 
  atomicAdd(&out[103], 0.1767766952966368*win[0]*fin[200]*subvol); 
  atomicAdd(&out[104], 0.1767766952966368*win[0]*fin[204]*subvol); 
  atomicAdd(&out[105], 0.1767766952966368*win[0]*fin[206]*subvol); 
  atomicAdd(&out[106], 0.1767766952966368*win[0]*fin[208]*subvol); 
  atomicAdd(&out[107], 0.1767766952966368*win[0]*fin[213]*subvol); 
#else
  out[0] += 0.1767766952966368*fin[0]*win[0]*subvol; 
  out[1] += 0.1767766952966368*win[0]*fin[1]*subvol; 
  out[2] += 0.1767766952966368*win[0]*fin[3]*subvol; 
  out[3] += 0.1767766952966368*win[0]*fin[4]*subvol; 
  out[4] += 0.1767766952966368*win[0]*fin[5]*subvol; 
  out[5] += 0.1767766952966368*win[0]*fin[6]*subvol; 
  out[6] += 0.1767766952966368*win[0]*fin[8]*subvol; 
  out[7] += 0.1767766952966368*win[0]*fin[10]*subvol; 
  out[8] += 0.1767766952966368*win[0]*fin[12]*subvol; 
  out[9] += 0.1767766952966368*win[0]*fin[13]*subvol; 
  out[10] += 0.1767766952966368*win[0]*fin[15]*subvol; 
  out[11] += 0.1767766952966368*win[0]*fin[16]*subvol; 
  out[12] += 0.1767766952966368*win[0]*fin[17]*subvol; 
  out[13] += 0.1767766952966368*win[0]*fin[19]*subvol; 
  out[14] += 0.1767766952966368*win[0]*fin[20]*subvol; 
  out[15] += 0.1767766952966368*win[0]*fin[21]*subvol; 
  out[16] += 0.1767766952966368*win[0]*fin[22]*subvol; 
  out[17] += 0.1767766952966368*win[0]*fin[23]*subvol; 
  out[18] += 0.1767766952966368*win[0]*fin[24]*subvol; 
  out[19] += 0.1767766952966368*win[0]*fin[27]*subvol; 
  out[20] += 0.1767766952966368*win[0]*fin[30]*subvol; 
  out[21] += 0.1767766952966368*win[0]*fin[32]*subvol; 
  out[22] += 0.1767766952966368*win[0]*fin[34]*subvol; 
  out[23] += 0.1767766952966368*win[0]*fin[36]*subvol; 
  out[24] += 0.1767766952966368*win[0]*fin[38]*subvol; 
  out[25] += 0.1767766952966368*win[0]*fin[40]*subvol; 
  out[26] += 0.1767766952966368*win[0]*fin[41]*subvol; 
  out[27] += 0.1767766952966368*win[0]*fin[43]*subvol; 
  out[28] += 0.1767766952966368*win[0]*fin[44]*subvol; 
  out[29] += 0.1767766952966368*win[0]*fin[45]*subvol; 
  out[30] += 0.1767766952966368*win[0]*fin[47]*subvol; 
  out[31] += 0.1767766952966368*win[0]*fin[48]*subvol; 
  out[32] += 0.1767766952966368*win[0]*fin[49]*subvol; 
  out[33] += 0.1767766952966368*win[0]*fin[51]*subvol; 
  out[34] += 0.1767766952966368*win[0]*fin[52]*subvol; 
  out[35] += 0.1767766952966368*win[0]*fin[53]*subvol; 
  out[36] += 0.1767766952966368*win[0]*fin[54]*subvol; 
  out[37] += 0.1767766952966368*win[0]*fin[55]*subvol; 
  out[38] += 0.1767766952966368*win[0]*fin[57]*subvol; 
  out[39] += 0.1767766952966368*win[0]*fin[58]*subvol; 
  out[40] += 0.1767766952966368*win[0]*fin[59]*subvol; 
  out[41] += 0.1767766952966368*win[0]*fin[63]*subvol; 
  out[42] += 0.1767766952966368*win[0]*fin[67]*subvol; 
  out[43] += 0.1767766952966368*win[0]*fin[70]*subvol; 
  out[44] += 0.1767766952966368*win[0]*fin[72]*subvol; 
  out[45] += 0.1767766952966368*win[0]*fin[74]*subvol; 
  out[46] += 0.1767766952966368*win[0]*fin[76]*subvol; 
  out[47] += 0.1767766952966368*win[0]*fin[78]*subvol; 
  out[48] += 0.1767766952966368*win[0]*fin[80]*subvol; 
  out[49] += 0.1767766952966368*win[0]*fin[82]*subvol; 
  out[50] += 0.1767766952966368*win[0]*fin[84]*subvol; 
  out[51] += 0.1767766952966368*win[0]*fin[86]*subvol; 
  out[52] += 0.1767766952966368*win[0]*fin[87]*subvol; 
  out[53] += 0.1767766952966368*win[0]*fin[89]*subvol; 
  out[54] += 0.1767766952966368*win[0]*fin[90]*subvol; 
  out[55] += 0.1767766952966368*win[0]*fin[91]*subvol; 
  out[56] += 0.1767766952966368*win[0]*fin[93]*subvol; 
  out[57] += 0.1767766952966368*win[0]*fin[94]*subvol; 
  out[58] += 0.1767766952966368*win[0]*fin[96]*subvol; 
  out[59] += 0.1767766952966368*win[0]*fin[98]*subvol; 
  out[60] += 0.1767766952966368*win[0]*fin[100]*subvol; 
  out[61] += 0.1767766952966368*win[0]*fin[101]*subvol; 
  out[62] += 0.1767766952966368*win[0]*fin[103]*subvol; 
  out[63] += 0.1767766952966368*win[0]*fin[104]*subvol; 
  out[64] += 0.1767766952966368*win[0]*fin[105]*subvol; 
  out[65] += 0.1767766952966368*win[0]*fin[106]*subvol; 
  out[66] += 0.1767766952966368*win[0]*fin[107]*subvol; 
  out[67] += 0.1767766952966368*win[0]*fin[112]*subvol; 
  out[68] += 0.1767766952966368*win[0]*fin[116]*subvol; 
  out[69] += 0.1767766952966368*win[0]*fin[120]*subvol; 
  out[70] += 0.1767766952966368*win[0]*fin[123]*subvol; 
  out[71] += 0.1767766952966368*win[0]*fin[125]*subvol; 
  out[72] += 0.1767766952966368*win[0]*fin[127]*subvol; 
  out[73] += 0.1767766952966368*win[0]*fin[129]*subvol; 
  out[74] += 0.1767766952966368*win[0]*fin[131]*subvol; 
  out[75] += 0.1767766952966368*win[0]*fin[133]*subvol; 
  out[76] += 0.1767766952966368*win[0]*fin[136]*subvol; 
  out[77] += 0.1767766952966368*win[0]*fin[139]*subvol; 
  out[78] += 0.1767766952966368*win[0]*fin[141]*subvol; 
  out[79] += 0.1767766952966368*win[0]*fin[143]*subvol; 
  out[80] += 0.1767766952966368*win[0]*fin[144]*subvol; 
  out[81] += 0.1767766952966368*win[0]*fin[146]*subvol; 
  out[82] += 0.1767766952966368*win[0]*fin[147]*subvol; 
  out[83] += 0.1767766952966368*win[0]*fin[148]*subvol; 
  out[84] += 0.1767766952966368*win[0]*fin[150]*subvol; 
  out[85] += 0.1767766952966368*win[0]*fin[151]*subvol; 
  out[86] += 0.1767766952966368*win[0]*fin[152]*subvol; 
  out[87] += 0.1767766952966368*win[0]*fin[154]*subvol; 
  out[88] += 0.1767766952966368*win[0]*fin[155]*subvol; 
  out[89] += 0.1767766952966368*win[0]*fin[161]*subvol; 
  out[90] += 0.1767766952966368*win[0]*fin[165]*subvol; 
  out[91] += 0.1767766952966368*win[0]*fin[170]*subvol; 
  out[92] += 0.1767766952966368*win[0]*fin[173]*subvol; 
  out[93] += 0.1767766952966368*win[0]*fin[175]*subvol; 
  out[94] += 0.1767766952966368*win[0]*fin[177]*subvol; 
  out[95] += 0.1767766952966368*win[0]*fin[179]*subvol; 
  out[96] += 0.1767766952966368*win[0]*fin[181]*subvol; 
  out[97] += 0.1767766952966368*win[0]*fin[183]*subvol; 
  out[98] += 0.1767766952966368*win[0]*fin[185]*subvol; 
  out[99] += 0.1767766952966368*win[0]*fin[187]*subvol; 
  out[100] += 0.1767766952966368*win[0]*fin[189]*subvol; 
  out[101] += 0.1767766952966368*win[0]*fin[190]*subvol; 
  out[102] += 0.1767766952966368*win[0]*fin[196]*subvol; 
  out[103] += 0.1767766952966368*win[0]*fin[200]*subvol; 
  out[104] += 0.1767766952966368*win[0]*fin[204]*subvol; 
  out[105] += 0.1767766952966368*win[0]*fin[206]*subvol; 
  out[106] += 0.1767766952966368*win[0]*fin[208]*subvol; 
  out[107] += 0.1767766952966368*win[0]*fin[213]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.25*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.25*win[0]*fin[1]*subvol); 
  atomicAdd(&out[3], 0.25*win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], 0.25*win[0]*fin[4]*subvol); 
  atomicAdd(&out[5], 0.25*win[0]*fin[5]*subvol); 
  atomicAdd(&out[6], 0.25*win[0]*fin[6]*subvol); 
  atomicAdd(&out[8], 0.25*win[0]*fin[8]*subvol); 
  atomicAdd(&out[10], 0.25*win[0]*fin[10]*subvol); 
  atomicAdd(&out[12], 0.25*win[0]*fin[12]*subvol); 
  atomicAdd(&out[13], 0.25*win[0]*fin[13]*subvol); 
  atomicAdd(&out[15], 0.25*win[0]*fin[15]*subvol); 
  atomicAdd(&out[16], 0.25*win[0]*fin[16]*subvol); 
  atomicAdd(&out[17], 0.25*win[0]*fin[17]*subvol); 
  atomicAdd(&out[19], 0.25*win[0]*fin[19]*subvol); 
  atomicAdd(&out[20], 0.25*win[0]*fin[20]*subvol); 
  atomicAdd(&out[21], 0.25*win[0]*fin[21]*subvol); 
  atomicAdd(&out[22], 0.25*win[0]*fin[22]*subvol); 
  atomicAdd(&out[23], 0.25*win[0]*fin[23]*subvol); 
  atomicAdd(&out[24], 0.25*win[0]*fin[24]*subvol); 
  atomicAdd(&out[27], 0.25*win[0]*fin[27]*subvol); 
  atomicAdd(&out[30], 0.25*win[0]*fin[30]*subvol); 
  atomicAdd(&out[32], 0.25*win[0]*fin[32]*subvol); 
  atomicAdd(&out[34], 0.25*win[0]*fin[34]*subvol); 
  atomicAdd(&out[36], 0.25*win[0]*fin[36]*subvol); 
  atomicAdd(&out[38], 0.25*win[0]*fin[38]*subvol); 
  atomicAdd(&out[40], 0.25*win[0]*fin[40]*subvol); 
  atomicAdd(&out[41], 0.25*win[0]*fin[41]*subvol); 
  atomicAdd(&out[43], 0.25*win[0]*fin[43]*subvol); 
  atomicAdd(&out[44], 0.25*win[0]*fin[44]*subvol); 
  atomicAdd(&out[45], 0.25*win[0]*fin[45]*subvol); 
  atomicAdd(&out[47], 0.25*win[0]*fin[47]*subvol); 
  atomicAdd(&out[48], 0.25*win[0]*fin[48]*subvol); 
  atomicAdd(&out[49], 0.25*win[0]*fin[49]*subvol); 
  atomicAdd(&out[51], 0.25*win[0]*fin[51]*subvol); 
  atomicAdd(&out[52], 0.25*win[0]*fin[52]*subvol); 
  atomicAdd(&out[53], 0.25*win[0]*fin[53]*subvol); 
  atomicAdd(&out[54], 0.25*win[0]*fin[54]*subvol); 
  atomicAdd(&out[55], 0.25*win[0]*fin[55]*subvol); 
  atomicAdd(&out[57], 0.25*win[0]*fin[57]*subvol); 
  atomicAdd(&out[58], 0.25*win[0]*fin[58]*subvol); 
  atomicAdd(&out[59], 0.25*win[0]*fin[59]*subvol); 
  atomicAdd(&out[63], 0.25*win[0]*fin[63]*subvol); 
  atomicAdd(&out[67], 0.25*win[0]*fin[67]*subvol); 
  atomicAdd(&out[70], 0.25*win[0]*fin[70]*subvol); 
  atomicAdd(&out[72], 0.25*win[0]*fin[72]*subvol); 
  atomicAdd(&out[74], 0.25*win[0]*fin[74]*subvol); 
  atomicAdd(&out[76], 0.25*win[0]*fin[76]*subvol); 
  atomicAdd(&out[78], 0.25*win[0]*fin[78]*subvol); 
  atomicAdd(&out[80], 0.25*win[0]*fin[80]*subvol); 
  atomicAdd(&out[82], 0.25*win[0]*fin[82]*subvol); 
  atomicAdd(&out[84], 0.25*win[0]*fin[84]*subvol); 
  atomicAdd(&out[86], 0.25*win[0]*fin[86]*subvol); 
  atomicAdd(&out[87], 0.25*win[0]*fin[87]*subvol); 
  atomicAdd(&out[89], 0.25*win[0]*fin[89]*subvol); 
  atomicAdd(&out[90], 0.25*win[0]*fin[90]*subvol); 
  atomicAdd(&out[91], 0.25*win[0]*fin[91]*subvol); 
  atomicAdd(&out[93], 0.25*win[0]*fin[93]*subvol); 
  atomicAdd(&out[94], 0.25*win[0]*fin[94]*subvol); 
  atomicAdd(&out[96], 0.25*win[0]*fin[96]*subvol); 
  atomicAdd(&out[98], 0.25*win[0]*fin[98]*subvol); 
  atomicAdd(&out[100], 0.25*win[0]*fin[100]*subvol); 
  atomicAdd(&out[101], 0.25*win[0]*fin[101]*subvol); 
  atomicAdd(&out[103], 0.25*win[0]*fin[103]*subvol); 
  atomicAdd(&out[104], 0.25*win[0]*fin[104]*subvol); 
  atomicAdd(&out[105], 0.25*win[0]*fin[105]*subvol); 
  atomicAdd(&out[106], 0.25*win[0]*fin[106]*subvol); 
  atomicAdd(&out[107], 0.25*win[0]*fin[107]*subvol); 
  atomicAdd(&out[112], 0.25*win[0]*fin[112]*subvol); 
  atomicAdd(&out[116], 0.25*win[0]*fin[116]*subvol); 
  atomicAdd(&out[120], 0.25*win[0]*fin[120]*subvol); 
  atomicAdd(&out[123], 0.25*win[0]*fin[123]*subvol); 
  atomicAdd(&out[125], 0.25*win[0]*fin[125]*subvol); 
  atomicAdd(&out[127], 0.25*win[0]*fin[127]*subvol); 
  atomicAdd(&out[129], 0.25*win[0]*fin[129]*subvol); 
  atomicAdd(&out[131], 0.25*win[0]*fin[131]*subvol); 
  atomicAdd(&out[133], 0.25*win[0]*fin[133]*subvol); 
  atomicAdd(&out[136], 0.25*win[0]*fin[136]*subvol); 
  atomicAdd(&out[139], 0.25*win[0]*fin[139]*subvol); 
  atomicAdd(&out[141], 0.25*win[0]*fin[141]*subvol); 
  atomicAdd(&out[143], 0.25*win[0]*fin[143]*subvol); 
  atomicAdd(&out[144], 0.25*win[0]*fin[144]*subvol); 
  atomicAdd(&out[146], 0.25*win[0]*fin[146]*subvol); 
  atomicAdd(&out[147], 0.25*win[0]*fin[147]*subvol); 
  atomicAdd(&out[148], 0.25*win[0]*fin[148]*subvol); 
  atomicAdd(&out[150], 0.25*win[0]*fin[150]*subvol); 
  atomicAdd(&out[151], 0.25*win[0]*fin[151]*subvol); 
  atomicAdd(&out[152], 0.25*win[0]*fin[152]*subvol); 
  atomicAdd(&out[154], 0.25*win[0]*fin[154]*subvol); 
  atomicAdd(&out[155], 0.25*win[0]*fin[155]*subvol); 
  atomicAdd(&out[161], 0.25*win[0]*fin[161]*subvol); 
  atomicAdd(&out[165], 0.25*win[0]*fin[165]*subvol); 
  atomicAdd(&out[170], 0.25*win[0]*fin[170]*subvol); 
  atomicAdd(&out[173], 0.25*win[0]*fin[173]*subvol); 
  atomicAdd(&out[175], 0.25*win[0]*fin[175]*subvol); 
  atomicAdd(&out[177], 0.25*win[0]*fin[177]*subvol); 
  atomicAdd(&out[179], 0.25*win[0]*fin[179]*subvol); 
  atomicAdd(&out[181], 0.25*win[0]*fin[181]*subvol); 
  atomicAdd(&out[183], 0.25*win[0]*fin[183]*subvol); 
  atomicAdd(&out[185], 0.25*win[0]*fin[185]*subvol); 
  atomicAdd(&out[187], 0.25*win[0]*fin[187]*subvol); 
  atomicAdd(&out[189], 0.25*win[0]*fin[189]*subvol); 
  atomicAdd(&out[190], 0.25*win[0]*fin[190]*subvol); 
  atomicAdd(&out[196], 0.25*win[0]*fin[196]*subvol); 
  atomicAdd(&out[200], 0.25*win[0]*fin[200]*subvol); 
  atomicAdd(&out[204], 0.25*win[0]*fin[204]*subvol); 
  atomicAdd(&out[206], 0.25*win[0]*fin[206]*subvol); 
  atomicAdd(&out[208], 0.25*win[0]*fin[208]*subvol); 
  atomicAdd(&out[213], 0.25*win[0]*fin[213]*subvol); 
#else
  out[0] += 0.25*fin[0]*win[0]*subvol; 
  out[1] += 0.25*win[0]*fin[1]*subvol; 
  out[3] += 0.25*win[0]*fin[3]*subvol; 
  out[4] += 0.25*win[0]*fin[4]*subvol; 
  out[5] += 0.25*win[0]*fin[5]*subvol; 
  out[6] += 0.25*win[0]*fin[6]*subvol; 
  out[8] += 0.25*win[0]*fin[8]*subvol; 
  out[10] += 0.25*win[0]*fin[10]*subvol; 
  out[12] += 0.25*win[0]*fin[12]*subvol; 
  out[13] += 0.25*win[0]*fin[13]*subvol; 
  out[15] += 0.25*win[0]*fin[15]*subvol; 
  out[16] += 0.25*win[0]*fin[16]*subvol; 
  out[17] += 0.25*win[0]*fin[17]*subvol; 
  out[19] += 0.25*win[0]*fin[19]*subvol; 
  out[20] += 0.25*win[0]*fin[20]*subvol; 
  out[21] += 0.25*win[0]*fin[21]*subvol; 
  out[22] += 0.25*win[0]*fin[22]*subvol; 
  out[23] += 0.25*win[0]*fin[23]*subvol; 
  out[24] += 0.25*win[0]*fin[24]*subvol; 
  out[27] += 0.25*win[0]*fin[27]*subvol; 
  out[30] += 0.25*win[0]*fin[30]*subvol; 
  out[32] += 0.25*win[0]*fin[32]*subvol; 
  out[34] += 0.25*win[0]*fin[34]*subvol; 
  out[36] += 0.25*win[0]*fin[36]*subvol; 
  out[38] += 0.25*win[0]*fin[38]*subvol; 
  out[40] += 0.25*win[0]*fin[40]*subvol; 
  out[41] += 0.25*win[0]*fin[41]*subvol; 
  out[43] += 0.25*win[0]*fin[43]*subvol; 
  out[44] += 0.25*win[0]*fin[44]*subvol; 
  out[45] += 0.25*win[0]*fin[45]*subvol; 
  out[47] += 0.25*win[0]*fin[47]*subvol; 
  out[48] += 0.25*win[0]*fin[48]*subvol; 
  out[49] += 0.25*win[0]*fin[49]*subvol; 
  out[51] += 0.25*win[0]*fin[51]*subvol; 
  out[52] += 0.25*win[0]*fin[52]*subvol; 
  out[53] += 0.25*win[0]*fin[53]*subvol; 
  out[54] += 0.25*win[0]*fin[54]*subvol; 
  out[55] += 0.25*win[0]*fin[55]*subvol; 
  out[57] += 0.25*win[0]*fin[57]*subvol; 
  out[58] += 0.25*win[0]*fin[58]*subvol; 
  out[59] += 0.25*win[0]*fin[59]*subvol; 
  out[63] += 0.25*win[0]*fin[63]*subvol; 
  out[67] += 0.25*win[0]*fin[67]*subvol; 
  out[70] += 0.25*win[0]*fin[70]*subvol; 
  out[72] += 0.25*win[0]*fin[72]*subvol; 
  out[74] += 0.25*win[0]*fin[74]*subvol; 
  out[76] += 0.25*win[0]*fin[76]*subvol; 
  out[78] += 0.25*win[0]*fin[78]*subvol; 
  out[80] += 0.25*win[0]*fin[80]*subvol; 
  out[82] += 0.25*win[0]*fin[82]*subvol; 
  out[84] += 0.25*win[0]*fin[84]*subvol; 
  out[86] += 0.25*win[0]*fin[86]*subvol; 
  out[87] += 0.25*win[0]*fin[87]*subvol; 
  out[89] += 0.25*win[0]*fin[89]*subvol; 
  out[90] += 0.25*win[0]*fin[90]*subvol; 
  out[91] += 0.25*win[0]*fin[91]*subvol; 
  out[93] += 0.25*win[0]*fin[93]*subvol; 
  out[94] += 0.25*win[0]*fin[94]*subvol; 
  out[96] += 0.25*win[0]*fin[96]*subvol; 
  out[98] += 0.25*win[0]*fin[98]*subvol; 
  out[100] += 0.25*win[0]*fin[100]*subvol; 
  out[101] += 0.25*win[0]*fin[101]*subvol; 
  out[103] += 0.25*win[0]*fin[103]*subvol; 
  out[104] += 0.25*win[0]*fin[104]*subvol; 
  out[105] += 0.25*win[0]*fin[105]*subvol; 
  out[106] += 0.25*win[0]*fin[106]*subvol; 
  out[107] += 0.25*win[0]*fin[107]*subvol; 
  out[112] += 0.25*win[0]*fin[112]*subvol; 
  out[116] += 0.25*win[0]*fin[116]*subvol; 
  out[120] += 0.25*win[0]*fin[120]*subvol; 
  out[123] += 0.25*win[0]*fin[123]*subvol; 
  out[125] += 0.25*win[0]*fin[125]*subvol; 
  out[127] += 0.25*win[0]*fin[127]*subvol; 
  out[129] += 0.25*win[0]*fin[129]*subvol; 
  out[131] += 0.25*win[0]*fin[131]*subvol; 
  out[133] += 0.25*win[0]*fin[133]*subvol; 
  out[136] += 0.25*win[0]*fin[136]*subvol; 
  out[139] += 0.25*win[0]*fin[139]*subvol; 
  out[141] += 0.25*win[0]*fin[141]*subvol; 
  out[143] += 0.25*win[0]*fin[143]*subvol; 
  out[144] += 0.25*win[0]*fin[144]*subvol; 
  out[146] += 0.25*win[0]*fin[146]*subvol; 
  out[147] += 0.25*win[0]*fin[147]*subvol; 
  out[148] += 0.25*win[0]*fin[148]*subvol; 
  out[150] += 0.25*win[0]*fin[150]*subvol; 
  out[151] += 0.25*win[0]*fin[151]*subvol; 
  out[152] += 0.25*win[0]*fin[152]*subvol; 
  out[154] += 0.25*win[0]*fin[154]*subvol; 
  out[155] += 0.25*win[0]*fin[155]*subvol; 
  out[161] += 0.25*win[0]*fin[161]*subvol; 
  out[165] += 0.25*win[0]*fin[165]*subvol; 
  out[170] += 0.25*win[0]*fin[170]*subvol; 
  out[173] += 0.25*win[0]*fin[173]*subvol; 
  out[175] += 0.25*win[0]*fin[175]*subvol; 
  out[177] += 0.25*win[0]*fin[177]*subvol; 
  out[179] += 0.25*win[0]*fin[179]*subvol; 
  out[181] += 0.25*win[0]*fin[181]*subvol; 
  out[183] += 0.25*win[0]*fin[183]*subvol; 
  out[185] += 0.25*win[0]*fin[185]*subvol; 
  out[187] += 0.25*win[0]*fin[187]*subvol; 
  out[189] += 0.25*win[0]*fin[189]*subvol; 
  out[190] += 0.25*win[0]*fin[190]*subvol; 
  out[196] += 0.25*win[0]*fin[196]*subvol; 
  out[200] += 0.25*win[0]*fin[200]*subvol; 
  out[204] += 0.25*win[0]*fin[204]*subvol; 
  out[206] += 0.25*win[0]*fin[206]*subvol; 
  out[208] += 0.25*win[0]*fin[208]*subvol; 
  out[213] += 0.25*win[0]*fin[213]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.25*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.25*win[0]*fin[3]*subvol); 
  atomicAdd(&out[2], 0.25*win[0]*fin[4]*subvol); 
  atomicAdd(&out[3], 0.25*win[0]*fin[5]*subvol); 
  atomicAdd(&out[4], 0.25*win[0]*fin[6]*subvol); 
  atomicAdd(&out[5], 0.25*win[0]*fin[12]*subvol); 
  atomicAdd(&out[6], 0.25*win[0]*fin[15]*subvol); 
  atomicAdd(&out[7], 0.25*win[0]*fin[16]*subvol); 
  atomicAdd(&out[8], 0.25*win[0]*fin[19]*subvol); 
  atomicAdd(&out[9], 0.25*win[0]*fin[20]*subvol); 
  atomicAdd(&out[10], 0.25*win[0]*fin[21]*subvol); 
  atomicAdd(&out[11], 0.25*win[0]*fin[22]*subvol); 
  atomicAdd(&out[12], 0.25*win[0]*fin[23]*subvol); 
  atomicAdd(&out[13], 0.25*win[0]*fin[24]*subvol); 
  atomicAdd(&out[14], 0.25*win[0]*fin[34]*subvol); 
  atomicAdd(&out[15], 0.25*win[0]*fin[40]*subvol); 
  atomicAdd(&out[16], 0.25*win[0]*fin[43]*subvol); 
  atomicAdd(&out[17], 0.25*win[0]*fin[44]*subvol); 
  atomicAdd(&out[18], 0.25*win[0]*fin[47]*subvol); 
  atomicAdd(&out[19], 0.25*win[0]*fin[48]*subvol); 
  atomicAdd(&out[20], 0.25*win[0]*fin[51]*subvol); 
  atomicAdd(&out[21], 0.25*win[0]*fin[52]*subvol); 
  atomicAdd(&out[22], 0.25*win[0]*fin[53]*subvol); 
  atomicAdd(&out[23], 0.25*win[0]*fin[54]*subvol); 
  atomicAdd(&out[24], 0.25*win[0]*fin[57]*subvol); 
  atomicAdd(&out[25], 0.25*win[0]*fin[58]*subvol); 
  atomicAdd(&out[26], 0.25*win[0]*fin[59]*subvol); 
  atomicAdd(&out[27], 0.25*win[0]*fin[74]*subvol); 
  atomicAdd(&out[28], 0.25*win[0]*fin[80]*subvol); 
  atomicAdd(&out[29], 0.25*win[0]*fin[86]*subvol); 
  atomicAdd(&out[30], 0.25*win[0]*fin[89]*subvol); 
  atomicAdd(&out[31], 0.25*win[0]*fin[90]*subvol); 
  atomicAdd(&out[32], 0.25*win[0]*fin[93]*subvol); 
  atomicAdd(&out[33], 0.25*win[0]*fin[94]*subvol); 
  atomicAdd(&out[34], 0.25*win[0]*fin[100]*subvol); 
  atomicAdd(&out[35], 0.25*win[0]*fin[103]*subvol); 
  atomicAdd(&out[36], 0.25*win[0]*fin[104]*subvol); 
  atomicAdd(&out[37], 0.25*win[0]*fin[105]*subvol); 
  atomicAdd(&out[38], 0.25*win[0]*fin[106]*subvol); 
  atomicAdd(&out[39], 0.25*win[0]*fin[107]*subvol); 
  atomicAdd(&out[40], 0.25*win[0]*fin[127]*subvol); 
  atomicAdd(&out[41], 0.25*win[0]*fin[133]*subvol); 
  atomicAdd(&out[42], 0.25*win[0]*fin[143]*subvol); 
  atomicAdd(&out[43], 0.25*win[0]*fin[146]*subvol); 
  atomicAdd(&out[44], 0.25*win[0]*fin[147]*subvol); 
  atomicAdd(&out[45], 0.25*win[0]*fin[150]*subvol); 
  atomicAdd(&out[46], 0.25*win[0]*fin[151]*subvol); 
  atomicAdd(&out[47], 0.25*win[0]*fin[154]*subvol); 
  atomicAdd(&out[48], 0.25*win[0]*fin[155]*subvol); 
  atomicAdd(&out[49], 0.25*win[0]*fin[177]*subvol); 
  atomicAdd(&out[50], 0.25*win[0]*fin[183]*subvol); 
  atomicAdd(&out[51], 0.25*win[0]*fin[189]*subvol); 
  atomicAdd(&out[52], 0.25*win[0]*fin[190]*subvol); 
  atomicAdd(&out[53], 0.25*win[0]*fin[208]*subvol); 
#else
  out[0] += 0.25*fin[0]*win[0]*subvol; 
  out[1] += 0.25*win[0]*fin[3]*subvol; 
  out[2] += 0.25*win[0]*fin[4]*subvol; 
  out[3] += 0.25*win[0]*fin[5]*subvol; 
  out[4] += 0.25*win[0]*fin[6]*subvol; 
  out[5] += 0.25*win[0]*fin[12]*subvol; 
  out[6] += 0.25*win[0]*fin[15]*subvol; 
  out[7] += 0.25*win[0]*fin[16]*subvol; 
  out[8] += 0.25*win[0]*fin[19]*subvol; 
  out[9] += 0.25*win[0]*fin[20]*subvol; 
  out[10] += 0.25*win[0]*fin[21]*subvol; 
  out[11] += 0.25*win[0]*fin[22]*subvol; 
  out[12] += 0.25*win[0]*fin[23]*subvol; 
  out[13] += 0.25*win[0]*fin[24]*subvol; 
  out[14] += 0.25*win[0]*fin[34]*subvol; 
  out[15] += 0.25*win[0]*fin[40]*subvol; 
  out[16] += 0.25*win[0]*fin[43]*subvol; 
  out[17] += 0.25*win[0]*fin[44]*subvol; 
  out[18] += 0.25*win[0]*fin[47]*subvol; 
  out[19] += 0.25*win[0]*fin[48]*subvol; 
  out[20] += 0.25*win[0]*fin[51]*subvol; 
  out[21] += 0.25*win[0]*fin[52]*subvol; 
  out[22] += 0.25*win[0]*fin[53]*subvol; 
  out[23] += 0.25*win[0]*fin[54]*subvol; 
  out[24] += 0.25*win[0]*fin[57]*subvol; 
  out[25] += 0.25*win[0]*fin[58]*subvol; 
  out[26] += 0.25*win[0]*fin[59]*subvol; 
  out[27] += 0.25*win[0]*fin[74]*subvol; 
  out[28] += 0.25*win[0]*fin[80]*subvol; 
  out[29] += 0.25*win[0]*fin[86]*subvol; 
  out[30] += 0.25*win[0]*fin[89]*subvol; 
  out[31] += 0.25*win[0]*fin[90]*subvol; 
  out[32] += 0.25*win[0]*fin[93]*subvol; 
  out[33] += 0.25*win[0]*fin[94]*subvol; 
  out[34] += 0.25*win[0]*fin[100]*subvol; 
  out[35] += 0.25*win[0]*fin[103]*subvol; 
  out[36] += 0.25*win[0]*fin[104]*subvol; 
  out[37] += 0.25*win[0]*fin[105]*subvol; 
  out[38] += 0.25*win[0]*fin[106]*subvol; 
  out[39] += 0.25*win[0]*fin[107]*subvol; 
  out[40] += 0.25*win[0]*fin[127]*subvol; 
  out[41] += 0.25*win[0]*fin[133]*subvol; 
  out[42] += 0.25*win[0]*fin[143]*subvol; 
  out[43] += 0.25*win[0]*fin[146]*subvol; 
  out[44] += 0.25*win[0]*fin[147]*subvol; 
  out[45] += 0.25*win[0]*fin[150]*subvol; 
  out[46] += 0.25*win[0]*fin[151]*subvol; 
  out[47] += 0.25*win[0]*fin[154]*subvol; 
  out[48] += 0.25*win[0]*fin[155]*subvol; 
  out[49] += 0.25*win[0]*fin[177]*subvol; 
  out[50] += 0.25*win[0]*fin[183]*subvol; 
  out[51] += 0.25*win[0]*fin[189]*subvol; 
  out[52] += 0.25*win[0]*fin[190]*subvol; 
  out[53] += 0.25*win[0]*fin[208]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[5]*subvol); 
  atomicAdd(&out[6], 0.5*win[0]*fin[6]*subvol); 
  atomicAdd(&out[12], 0.5*win[0]*fin[12]*subvol); 
  atomicAdd(&out[15], 0.5*win[0]*fin[15]*subvol); 
  atomicAdd(&out[16], 0.5*win[0]*fin[16]*subvol); 
  atomicAdd(&out[19], 0.5*win[0]*fin[19]*subvol); 
  atomicAdd(&out[20], 0.5*win[0]*fin[20]*subvol); 
  atomicAdd(&out[21], 0.5*win[0]*fin[21]*subvol); 
  atomicAdd(&out[22], 0.5*win[0]*fin[22]*subvol); 
  atomicAdd(&out[23], 0.5*win[0]*fin[23]*subvol); 
  atomicAdd(&out[24], 0.5*win[0]*fin[24]*subvol); 
  atomicAdd(&out[34], 0.5*win[0]*fin[34]*subvol); 
  atomicAdd(&out[40], 0.5*win[0]*fin[40]*subvol); 
  atomicAdd(&out[43], 0.5*win[0]*fin[43]*subvol); 
  atomicAdd(&out[44], 0.5*win[0]*fin[44]*subvol); 
  atomicAdd(&out[47], 0.5*win[0]*fin[47]*subvol); 
  atomicAdd(&out[48], 0.5*win[0]*fin[48]*subvol); 
  atomicAdd(&out[51], 0.5*win[0]*fin[51]*subvol); 
  atomicAdd(&out[52], 0.5*win[0]*fin[52]*subvol); 
  atomicAdd(&out[53], 0.5*win[0]*fin[53]*subvol); 
  atomicAdd(&out[54], 0.5*win[0]*fin[54]*subvol); 
  atomicAdd(&out[57], 0.5*win[0]*fin[57]*subvol); 
  atomicAdd(&out[58], 0.5*win[0]*fin[58]*subvol); 
  atomicAdd(&out[59], 0.5*win[0]*fin[59]*subvol); 
  atomicAdd(&out[74], 0.5*win[0]*fin[74]*subvol); 
  atomicAdd(&out[80], 0.5*win[0]*fin[80]*subvol); 
  atomicAdd(&out[86], 0.5*win[0]*fin[86]*subvol); 
  atomicAdd(&out[89], 0.5*win[0]*fin[89]*subvol); 
  atomicAdd(&out[90], 0.5*win[0]*fin[90]*subvol); 
  atomicAdd(&out[93], 0.5*win[0]*fin[93]*subvol); 
  atomicAdd(&out[94], 0.5*win[0]*fin[94]*subvol); 
  atomicAdd(&out[100], 0.5*win[0]*fin[100]*subvol); 
  atomicAdd(&out[103], 0.5*win[0]*fin[103]*subvol); 
  atomicAdd(&out[104], 0.5*win[0]*fin[104]*subvol); 
  atomicAdd(&out[105], 0.5*win[0]*fin[105]*subvol); 
  atomicAdd(&out[106], 0.5*win[0]*fin[106]*subvol); 
  atomicAdd(&out[107], 0.5*win[0]*fin[107]*subvol); 
  atomicAdd(&out[127], 0.5*win[0]*fin[127]*subvol); 
  atomicAdd(&out[133], 0.5*win[0]*fin[133]*subvol); 
  atomicAdd(&out[143], 0.5*win[0]*fin[143]*subvol); 
  atomicAdd(&out[146], 0.5*win[0]*fin[146]*subvol); 
  atomicAdd(&out[147], 0.5*win[0]*fin[147]*subvol); 
  atomicAdd(&out[150], 0.5*win[0]*fin[150]*subvol); 
  atomicAdd(&out[151], 0.5*win[0]*fin[151]*subvol); 
  atomicAdd(&out[154], 0.5*win[0]*fin[154]*subvol); 
  atomicAdd(&out[155], 0.5*win[0]*fin[155]*subvol); 
  atomicAdd(&out[177], 0.5*win[0]*fin[177]*subvol); 
  atomicAdd(&out[183], 0.5*win[0]*fin[183]*subvol); 
  atomicAdd(&out[189], 0.5*win[0]*fin[189]*subvol); 
  atomicAdd(&out[190], 0.5*win[0]*fin[190]*subvol); 
  atomicAdd(&out[208], 0.5*win[0]*fin[208]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[3] += 0.5*win[0]*fin[3]*subvol; 
  out[4] += 0.5*win[0]*fin[4]*subvol; 
  out[5] += 0.5*win[0]*fin[5]*subvol; 
  out[6] += 0.5*win[0]*fin[6]*subvol; 
  out[12] += 0.5*win[0]*fin[12]*subvol; 
  out[15] += 0.5*win[0]*fin[15]*subvol; 
  out[16] += 0.5*win[0]*fin[16]*subvol; 
  out[19] += 0.5*win[0]*fin[19]*subvol; 
  out[20] += 0.5*win[0]*fin[20]*subvol; 
  out[21] += 0.5*win[0]*fin[21]*subvol; 
  out[22] += 0.5*win[0]*fin[22]*subvol; 
  out[23] += 0.5*win[0]*fin[23]*subvol; 
  out[24] += 0.5*win[0]*fin[24]*subvol; 
  out[34] += 0.5*win[0]*fin[34]*subvol; 
  out[40] += 0.5*win[0]*fin[40]*subvol; 
  out[43] += 0.5*win[0]*fin[43]*subvol; 
  out[44] += 0.5*win[0]*fin[44]*subvol; 
  out[47] += 0.5*win[0]*fin[47]*subvol; 
  out[48] += 0.5*win[0]*fin[48]*subvol; 
  out[51] += 0.5*win[0]*fin[51]*subvol; 
  out[52] += 0.5*win[0]*fin[52]*subvol; 
  out[53] += 0.5*win[0]*fin[53]*subvol; 
  out[54] += 0.5*win[0]*fin[54]*subvol; 
  out[57] += 0.5*win[0]*fin[57]*subvol; 
  out[58] += 0.5*win[0]*fin[58]*subvol; 
  out[59] += 0.5*win[0]*fin[59]*subvol; 
  out[74] += 0.5*win[0]*fin[74]*subvol; 
  out[80] += 0.5*win[0]*fin[80]*subvol; 
  out[86] += 0.5*win[0]*fin[86]*subvol; 
  out[89] += 0.5*win[0]*fin[89]*subvol; 
  out[90] += 0.5*win[0]*fin[90]*subvol; 
  out[93] += 0.5*win[0]*fin[93]*subvol; 
  out[94] += 0.5*win[0]*fin[94]*subvol; 
  out[100] += 0.5*win[0]*fin[100]*subvol; 
  out[103] += 0.5*win[0]*fin[103]*subvol; 
  out[104] += 0.5*win[0]*fin[104]*subvol; 
  out[105] += 0.5*win[0]*fin[105]*subvol; 
  out[106] += 0.5*win[0]*fin[106]*subvol; 
  out[107] += 0.5*win[0]*fin[107]*subvol; 
  out[127] += 0.5*win[0]*fin[127]*subvol; 
  out[133] += 0.5*win[0]*fin[133]*subvol; 
  out[143] += 0.5*win[0]*fin[143]*subvol; 
  out[146] += 0.5*win[0]*fin[146]*subvol; 
  out[147] += 0.5*win[0]*fin[147]*subvol; 
  out[150] += 0.5*win[0]*fin[150]*subvol; 
  out[151] += 0.5*win[0]*fin[151]*subvol; 
  out[154] += 0.5*win[0]*fin[154]*subvol; 
  out[155] += 0.5*win[0]*fin[155]*subvol; 
  out[177] += 0.5*win[0]*fin[177]*subvol; 
  out[183] += 0.5*win[0]*fin[183]*subvol; 
  out[189] += 0.5*win[0]*fin[189]*subvol; 
  out[190] += 0.5*win[0]*fin[190]*subvol; 
  out[208] += 0.5*win[0]*fin[208]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.1767766952966368*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.1767766952966368*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.1767766952966368*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.1767766952966368*win[0]*fin[4]*subvol); 
  atomicAdd(&out[4], 0.1767766952966368*win[0]*fin[5]*subvol); 
  atomicAdd(&out[5], 0.1767766952966368*win[0]*fin[6]*subvol); 
  atomicAdd(&out[6], 0.1767766952966368*win[0]*fin[7]*subvol); 
  atomicAdd(&out[7], 0.1767766952966368*win[0]*fin[10]*subvol); 
  atomicAdd(&out[8], 0.1767766952966368*win[0]*fin[11]*subvol); 
  atomicAdd(&out[9], 0.1767766952966368*win[0]*fin[13]*subvol); 
  atomicAdd(&out[10], 0.1767766952966368*win[0]*fin[14]*subvol); 
  atomicAdd(&out[11], 0.1767766952966368*win[0]*fin[16]*subvol); 
  atomicAdd(&out[12], 0.1767766952966368*win[0]*fin[17]*subvol); 
  atomicAdd(&out[13], 0.1767766952966368*win[0]*fin[18]*subvol); 
  atomicAdd(&out[14], 0.1767766952966368*win[0]*fin[20]*subvol); 
  atomicAdd(&out[15], 0.1767766952966368*win[0]*fin[21]*subvol); 
  atomicAdd(&out[16], 0.1767766952966368*win[0]*fin[22]*subvol); 
  atomicAdd(&out[17], 0.1767766952966368*win[0]*fin[23]*subvol); 
  atomicAdd(&out[18], 0.1767766952966368*win[0]*fin[24]*subvol); 
  atomicAdd(&out[19], 0.1767766952966368*win[0]*fin[26]*subvol); 
  atomicAdd(&out[20], 0.1767766952966368*win[0]*fin[29]*subvol); 
  atomicAdd(&out[21], 0.1767766952966368*win[0]*fin[32]*subvol); 
  atomicAdd(&out[22], 0.1767766952966368*win[0]*fin[33]*subvol); 
  atomicAdd(&out[23], 0.1767766952966368*win[0]*fin[35]*subvol); 
  atomicAdd(&out[24], 0.1767766952966368*win[0]*fin[38]*subvol); 
  atomicAdd(&out[25], 0.1767766952966368*win[0]*fin[39]*subvol); 
  atomicAdd(&out[26], 0.1767766952966368*win[0]*fin[41]*subvol); 
  atomicAdd(&out[27], 0.1767766952966368*win[0]*fin[42]*subvol); 
  atomicAdd(&out[28], 0.1767766952966368*win[0]*fin[44]*subvol); 
  atomicAdd(&out[29], 0.1767766952966368*win[0]*fin[45]*subvol); 
  atomicAdd(&out[30], 0.1767766952966368*win[0]*fin[46]*subvol); 
  atomicAdd(&out[31], 0.1767766952966368*win[0]*fin[48]*subvol); 
  atomicAdd(&out[32], 0.1767766952966368*win[0]*fin[49]*subvol); 
  atomicAdd(&out[33], 0.1767766952966368*win[0]*fin[50]*subvol); 
  atomicAdd(&out[34], 0.1767766952966368*win[0]*fin[52]*subvol); 
  atomicAdd(&out[35], 0.1767766952966368*win[0]*fin[53]*subvol); 
  atomicAdd(&out[36], 0.1767766952966368*win[0]*fin[54]*subvol); 
  atomicAdd(&out[37], 0.1767766952966368*win[0]*fin[55]*subvol); 
  atomicAdd(&out[38], 0.1767766952966368*win[0]*fin[56]*subvol); 
  atomicAdd(&out[39], 0.1767766952966368*win[0]*fin[58]*subvol); 
  atomicAdd(&out[40], 0.1767766952966368*win[0]*fin[59]*subvol); 
  atomicAdd(&out[41], 0.1767766952966368*win[0]*fin[62]*subvol); 
  atomicAdd(&out[42], 0.1767766952966368*win[0]*fin[66]*subvol); 
  atomicAdd(&out[43], 0.1767766952966368*win[0]*fin[69]*subvol); 
  atomicAdd(&out[44], 0.1767766952966368*win[0]*fin[72]*subvol); 
  atomicAdd(&out[45], 0.1767766952966368*win[0]*fin[73]*subvol); 
  atomicAdd(&out[46], 0.1767766952966368*win[0]*fin[75]*subvol); 
  atomicAdd(&out[47], 0.1767766952966368*win[0]*fin[78]*subvol); 
  atomicAdd(&out[48], 0.1767766952966368*win[0]*fin[79]*subvol); 
  atomicAdd(&out[49], 0.1767766952966368*win[0]*fin[81]*subvol); 
  atomicAdd(&out[50], 0.1767766952966368*win[0]*fin[84]*subvol); 
  atomicAdd(&out[51], 0.1767766952966368*win[0]*fin[85]*subvol); 
  atomicAdd(&out[52], 0.1767766952966368*win[0]*fin[87]*subvol); 
  atomicAdd(&out[53], 0.1767766952966368*win[0]*fin[88]*subvol); 
  atomicAdd(&out[54], 0.1767766952966368*win[0]*fin[90]*subvol); 
  atomicAdd(&out[55], 0.1767766952966368*win[0]*fin[91]*subvol); 
  atomicAdd(&out[56], 0.1767766952966368*win[0]*fin[92]*subvol); 
  atomicAdd(&out[57], 0.1767766952966368*win[0]*fin[94]*subvol); 
  atomicAdd(&out[58], 0.1767766952966368*win[0]*fin[95]*subvol); 
  atomicAdd(&out[59], 0.1767766952966368*win[0]*fin[98]*subvol); 
  atomicAdd(&out[60], 0.1767766952966368*win[0]*fin[99]*subvol); 
  atomicAdd(&out[61], 0.1767766952966368*win[0]*fin[101]*subvol); 
  atomicAdd(&out[62], 0.1767766952966368*win[0]*fin[102]*subvol); 
  atomicAdd(&out[63], 0.1767766952966368*win[0]*fin[104]*subvol); 
  atomicAdd(&out[64], 0.1767766952966368*win[0]*fin[105]*subvol); 
  atomicAdd(&out[65], 0.1767766952966368*win[0]*fin[106]*subvol); 
  atomicAdd(&out[66], 0.1767766952966368*win[0]*fin[107]*subvol); 
  atomicAdd(&out[67], 0.1767766952966368*win[0]*fin[111]*subvol); 
  atomicAdd(&out[68], 0.1767766952966368*win[0]*fin[115]*subvol); 
  atomicAdd(&out[69], 0.1767766952966368*win[0]*fin[119]*subvol); 
  atomicAdd(&out[70], 0.1767766952966368*win[0]*fin[122]*subvol); 
  atomicAdd(&out[71], 0.1767766952966368*win[0]*fin[125]*subvol); 
  atomicAdd(&out[72], 0.1767766952966368*win[0]*fin[126]*subvol); 
  atomicAdd(&out[73], 0.1767766952966368*win[0]*fin[128]*subvol); 
  atomicAdd(&out[74], 0.1767766952966368*win[0]*fin[131]*subvol); 
  atomicAdd(&out[75], 0.1767766952966368*win[0]*fin[132]*subvol); 
  atomicAdd(&out[76], 0.1767766952966368*win[0]*fin[135]*subvol); 
  atomicAdd(&out[77], 0.1767766952966368*win[0]*fin[138]*subvol); 
  atomicAdd(&out[78], 0.1767766952966368*win[0]*fin[141]*subvol); 
  atomicAdd(&out[79], 0.1767766952966368*win[0]*fin[142]*subvol); 
  atomicAdd(&out[80], 0.1767766952966368*win[0]*fin[144]*subvol); 
  atomicAdd(&out[81], 0.1767766952966368*win[0]*fin[145]*subvol); 
  atomicAdd(&out[82], 0.1767766952966368*win[0]*fin[147]*subvol); 
  atomicAdd(&out[83], 0.1767766952966368*win[0]*fin[148]*subvol); 
  atomicAdd(&out[84], 0.1767766952966368*win[0]*fin[149]*subvol); 
  atomicAdd(&out[85], 0.1767766952966368*win[0]*fin[151]*subvol); 
  atomicAdd(&out[86], 0.1767766952966368*win[0]*fin[152]*subvol); 
  atomicAdd(&out[87], 0.1767766952966368*win[0]*fin[153]*subvol); 
  atomicAdd(&out[88], 0.1767766952966368*win[0]*fin[155]*subvol); 
  atomicAdd(&out[89], 0.1767766952966368*win[0]*fin[160]*subvol); 
  atomicAdd(&out[90], 0.1767766952966368*win[0]*fin[164]*subvol); 
  atomicAdd(&out[91], 0.1767766952966368*win[0]*fin[169]*subvol); 
  atomicAdd(&out[92], 0.1767766952966368*win[0]*fin[172]*subvol); 
  atomicAdd(&out[93], 0.1767766952966368*win[0]*fin[175]*subvol); 
  atomicAdd(&out[94], 0.1767766952966368*win[0]*fin[176]*subvol); 
  atomicAdd(&out[95], 0.1767766952966368*win[0]*fin[178]*subvol); 
  atomicAdd(&out[96], 0.1767766952966368*win[0]*fin[181]*subvol); 
  atomicAdd(&out[97], 0.1767766952966368*win[0]*fin[182]*subvol); 
  atomicAdd(&out[98], 0.1767766952966368*win[0]*fin[184]*subvol); 
  atomicAdd(&out[99], 0.1767766952966368*win[0]*fin[187]*subvol); 
  atomicAdd(&out[100], 0.1767766952966368*win[0]*fin[188]*subvol); 
  atomicAdd(&out[101], 0.1767766952966368*win[0]*fin[190]*subvol); 
  atomicAdd(&out[102], 0.1767766952966368*win[0]*fin[195]*subvol); 
  atomicAdd(&out[103], 0.1767766952966368*win[0]*fin[199]*subvol); 
  atomicAdd(&out[104], 0.1767766952966368*win[0]*fin[203]*subvol); 
  atomicAdd(&out[105], 0.1767766952966368*win[0]*fin[206]*subvol); 
  atomicAdd(&out[106], 0.1767766952966368*win[0]*fin[207]*subvol); 
  atomicAdd(&out[107], 0.1767766952966368*win[0]*fin[212]*subvol); 
#else
  out[0] += 0.1767766952966368*fin[0]*win[0]*subvol; 
  out[1] += 0.1767766952966368*win[0]*fin[1]*subvol; 
  out[2] += 0.1767766952966368*win[0]*fin[2]*subvol; 
  out[3] += 0.1767766952966368*win[0]*fin[4]*subvol; 
  out[4] += 0.1767766952966368*win[0]*fin[5]*subvol; 
  out[5] += 0.1767766952966368*win[0]*fin[6]*subvol; 
  out[6] += 0.1767766952966368*win[0]*fin[7]*subvol; 
  out[7] += 0.1767766952966368*win[0]*fin[10]*subvol; 
  out[8] += 0.1767766952966368*win[0]*fin[11]*subvol; 
  out[9] += 0.1767766952966368*win[0]*fin[13]*subvol; 
  out[10] += 0.1767766952966368*win[0]*fin[14]*subvol; 
  out[11] += 0.1767766952966368*win[0]*fin[16]*subvol; 
  out[12] += 0.1767766952966368*win[0]*fin[17]*subvol; 
  out[13] += 0.1767766952966368*win[0]*fin[18]*subvol; 
  out[14] += 0.1767766952966368*win[0]*fin[20]*subvol; 
  out[15] += 0.1767766952966368*win[0]*fin[21]*subvol; 
  out[16] += 0.1767766952966368*win[0]*fin[22]*subvol; 
  out[17] += 0.1767766952966368*win[0]*fin[23]*subvol; 
  out[18] += 0.1767766952966368*win[0]*fin[24]*subvol; 
  out[19] += 0.1767766952966368*win[0]*fin[26]*subvol; 
  out[20] += 0.1767766952966368*win[0]*fin[29]*subvol; 
  out[21] += 0.1767766952966368*win[0]*fin[32]*subvol; 
  out[22] += 0.1767766952966368*win[0]*fin[33]*subvol; 
  out[23] += 0.1767766952966368*win[0]*fin[35]*subvol; 
  out[24] += 0.1767766952966368*win[0]*fin[38]*subvol; 
  out[25] += 0.1767766952966368*win[0]*fin[39]*subvol; 
  out[26] += 0.1767766952966368*win[0]*fin[41]*subvol; 
  out[27] += 0.1767766952966368*win[0]*fin[42]*subvol; 
  out[28] += 0.1767766952966368*win[0]*fin[44]*subvol; 
  out[29] += 0.1767766952966368*win[0]*fin[45]*subvol; 
  out[30] += 0.1767766952966368*win[0]*fin[46]*subvol; 
  out[31] += 0.1767766952966368*win[0]*fin[48]*subvol; 
  out[32] += 0.1767766952966368*win[0]*fin[49]*subvol; 
  out[33] += 0.1767766952966368*win[0]*fin[50]*subvol; 
  out[34] += 0.1767766952966368*win[0]*fin[52]*subvol; 
  out[35] += 0.1767766952966368*win[0]*fin[53]*subvol; 
  out[36] += 0.1767766952966368*win[0]*fin[54]*subvol; 
  out[37] += 0.1767766952966368*win[0]*fin[55]*subvol; 
  out[38] += 0.1767766952966368*win[0]*fin[56]*subvol; 
  out[39] += 0.1767766952966368*win[0]*fin[58]*subvol; 
  out[40] += 0.1767766952966368*win[0]*fin[59]*subvol; 
  out[41] += 0.1767766952966368*win[0]*fin[62]*subvol; 
  out[42] += 0.1767766952966368*win[0]*fin[66]*subvol; 
  out[43] += 0.1767766952966368*win[0]*fin[69]*subvol; 
  out[44] += 0.1767766952966368*win[0]*fin[72]*subvol; 
  out[45] += 0.1767766952966368*win[0]*fin[73]*subvol; 
  out[46] += 0.1767766952966368*win[0]*fin[75]*subvol; 
  out[47] += 0.1767766952966368*win[0]*fin[78]*subvol; 
  out[48] += 0.1767766952966368*win[0]*fin[79]*subvol; 
  out[49] += 0.1767766952966368*win[0]*fin[81]*subvol; 
  out[50] += 0.1767766952966368*win[0]*fin[84]*subvol; 
  out[51] += 0.1767766952966368*win[0]*fin[85]*subvol; 
  out[52] += 0.1767766952966368*win[0]*fin[87]*subvol; 
  out[53] += 0.1767766952966368*win[0]*fin[88]*subvol; 
  out[54] += 0.1767766952966368*win[0]*fin[90]*subvol; 
  out[55] += 0.1767766952966368*win[0]*fin[91]*subvol; 
  out[56] += 0.1767766952966368*win[0]*fin[92]*subvol; 
  out[57] += 0.1767766952966368*win[0]*fin[94]*subvol; 
  out[58] += 0.1767766952966368*win[0]*fin[95]*subvol; 
  out[59] += 0.1767766952966368*win[0]*fin[98]*subvol; 
  out[60] += 0.1767766952966368*win[0]*fin[99]*subvol; 
  out[61] += 0.1767766952966368*win[0]*fin[101]*subvol; 
  out[62] += 0.1767766952966368*win[0]*fin[102]*subvol; 
  out[63] += 0.1767766952966368*win[0]*fin[104]*subvol; 
  out[64] += 0.1767766952966368*win[0]*fin[105]*subvol; 
  out[65] += 0.1767766952966368*win[0]*fin[106]*subvol; 
  out[66] += 0.1767766952966368*win[0]*fin[107]*subvol; 
  out[67] += 0.1767766952966368*win[0]*fin[111]*subvol; 
  out[68] += 0.1767766952966368*win[0]*fin[115]*subvol; 
  out[69] += 0.1767766952966368*win[0]*fin[119]*subvol; 
  out[70] += 0.1767766952966368*win[0]*fin[122]*subvol; 
  out[71] += 0.1767766952966368*win[0]*fin[125]*subvol; 
  out[72] += 0.1767766952966368*win[0]*fin[126]*subvol; 
  out[73] += 0.1767766952966368*win[0]*fin[128]*subvol; 
  out[74] += 0.1767766952966368*win[0]*fin[131]*subvol; 
  out[75] += 0.1767766952966368*win[0]*fin[132]*subvol; 
  out[76] += 0.1767766952966368*win[0]*fin[135]*subvol; 
  out[77] += 0.1767766952966368*win[0]*fin[138]*subvol; 
  out[78] += 0.1767766952966368*win[0]*fin[141]*subvol; 
  out[79] += 0.1767766952966368*win[0]*fin[142]*subvol; 
  out[80] += 0.1767766952966368*win[0]*fin[144]*subvol; 
  out[81] += 0.1767766952966368*win[0]*fin[145]*subvol; 
  out[82] += 0.1767766952966368*win[0]*fin[147]*subvol; 
  out[83] += 0.1767766952966368*win[0]*fin[148]*subvol; 
  out[84] += 0.1767766952966368*win[0]*fin[149]*subvol; 
  out[85] += 0.1767766952966368*win[0]*fin[151]*subvol; 
  out[86] += 0.1767766952966368*win[0]*fin[152]*subvol; 
  out[87] += 0.1767766952966368*win[0]*fin[153]*subvol; 
  out[88] += 0.1767766952966368*win[0]*fin[155]*subvol; 
  out[89] += 0.1767766952966368*win[0]*fin[160]*subvol; 
  out[90] += 0.1767766952966368*win[0]*fin[164]*subvol; 
  out[91] += 0.1767766952966368*win[0]*fin[169]*subvol; 
  out[92] += 0.1767766952966368*win[0]*fin[172]*subvol; 
  out[93] += 0.1767766952966368*win[0]*fin[175]*subvol; 
  out[94] += 0.1767766952966368*win[0]*fin[176]*subvol; 
  out[95] += 0.1767766952966368*win[0]*fin[178]*subvol; 
  out[96] += 0.1767766952966368*win[0]*fin[181]*subvol; 
  out[97] += 0.1767766952966368*win[0]*fin[182]*subvol; 
  out[98] += 0.1767766952966368*win[0]*fin[184]*subvol; 
  out[99] += 0.1767766952966368*win[0]*fin[187]*subvol; 
  out[100] += 0.1767766952966368*win[0]*fin[188]*subvol; 
  out[101] += 0.1767766952966368*win[0]*fin[190]*subvol; 
  out[102] += 0.1767766952966368*win[0]*fin[195]*subvol; 
  out[103] += 0.1767766952966368*win[0]*fin[199]*subvol; 
  out[104] += 0.1767766952966368*win[0]*fin[203]*subvol; 
  out[105] += 0.1767766952966368*win[0]*fin[206]*subvol; 
  out[106] += 0.1767766952966368*win[0]*fin[207]*subvol; 
  out[107] += 0.1767766952966368*win[0]*fin[212]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.25*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.25*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.25*win[0]*fin[2]*subvol); 
  atomicAdd(&out[4], 0.25*win[0]*fin[4]*subvol); 
  atomicAdd(&out[5], 0.25*win[0]*fin[5]*subvol); 
  atomicAdd(&out[6], 0.25*win[0]*fin[6]*subvol); 
  atomicAdd(&out[7], 0.25*win[0]*fin[7]*subvol); 
  atomicAdd(&out[10], 0.25*win[0]*fin[10]*subvol); 
  atomicAdd(&out[11], 0.25*win[0]*fin[11]*subvol); 
  atomicAdd(&out[13], 0.25*win[0]*fin[13]*subvol); 
  atomicAdd(&out[14], 0.25*win[0]*fin[14]*subvol); 
  atomicAdd(&out[16], 0.25*win[0]*fin[16]*subvol); 
  atomicAdd(&out[17], 0.25*win[0]*fin[17]*subvol); 
  atomicAdd(&out[18], 0.25*win[0]*fin[18]*subvol); 
  atomicAdd(&out[20], 0.25*win[0]*fin[20]*subvol); 
  atomicAdd(&out[21], 0.25*win[0]*fin[21]*subvol); 
  atomicAdd(&out[22], 0.25*win[0]*fin[22]*subvol); 
  atomicAdd(&out[23], 0.25*win[0]*fin[23]*subvol); 
  atomicAdd(&out[24], 0.25*win[0]*fin[24]*subvol); 
  atomicAdd(&out[26], 0.25*win[0]*fin[26]*subvol); 
  atomicAdd(&out[29], 0.25*win[0]*fin[29]*subvol); 
  atomicAdd(&out[32], 0.25*win[0]*fin[32]*subvol); 
  atomicAdd(&out[33], 0.25*win[0]*fin[33]*subvol); 
  atomicAdd(&out[35], 0.25*win[0]*fin[35]*subvol); 
  atomicAdd(&out[38], 0.25*win[0]*fin[38]*subvol); 
  atomicAdd(&out[39], 0.25*win[0]*fin[39]*subvol); 
  atomicAdd(&out[41], 0.25*win[0]*fin[41]*subvol); 
  atomicAdd(&out[42], 0.25*win[0]*fin[42]*subvol); 
  atomicAdd(&out[44], 0.25*win[0]*fin[44]*subvol); 
  atomicAdd(&out[45], 0.25*win[0]*fin[45]*subvol); 
  atomicAdd(&out[46], 0.25*win[0]*fin[46]*subvol); 
  atomicAdd(&out[48], 0.25*win[0]*fin[48]*subvol); 
  atomicAdd(&out[49], 0.25*win[0]*fin[49]*subvol); 
  atomicAdd(&out[50], 0.25*win[0]*fin[50]*subvol); 
  atomicAdd(&out[52], 0.25*win[0]*fin[52]*subvol); 
  atomicAdd(&out[53], 0.25*win[0]*fin[53]*subvol); 
  atomicAdd(&out[54], 0.25*win[0]*fin[54]*subvol); 
  atomicAdd(&out[55], 0.25*win[0]*fin[55]*subvol); 
  atomicAdd(&out[56], 0.25*win[0]*fin[56]*subvol); 
  atomicAdd(&out[58], 0.25*win[0]*fin[58]*subvol); 
  atomicAdd(&out[59], 0.25*win[0]*fin[59]*subvol); 
  atomicAdd(&out[62], 0.25*win[0]*fin[62]*subvol); 
  atomicAdd(&out[66], 0.25*win[0]*fin[66]*subvol); 
  atomicAdd(&out[69], 0.25*win[0]*fin[69]*subvol); 
  atomicAdd(&out[72], 0.25*win[0]*fin[72]*subvol); 
  atomicAdd(&out[73], 0.25*win[0]*fin[73]*subvol); 
  atomicAdd(&out[75], 0.25*win[0]*fin[75]*subvol); 
  atomicAdd(&out[78], 0.25*win[0]*fin[78]*subvol); 
  atomicAdd(&out[79], 0.25*win[0]*fin[79]*subvol); 
  atomicAdd(&out[81], 0.25*win[0]*fin[81]*subvol); 
  atomicAdd(&out[84], 0.25*win[0]*fin[84]*subvol); 
  atomicAdd(&out[85], 0.25*win[0]*fin[85]*subvol); 
  atomicAdd(&out[87], 0.25*win[0]*fin[87]*subvol); 
  atomicAdd(&out[88], 0.25*win[0]*fin[88]*subvol); 
  atomicAdd(&out[90], 0.25*win[0]*fin[90]*subvol); 
  atomicAdd(&out[91], 0.25*win[0]*fin[91]*subvol); 
  atomicAdd(&out[92], 0.25*win[0]*fin[92]*subvol); 
  atomicAdd(&out[94], 0.25*win[0]*fin[94]*subvol); 
  atomicAdd(&out[95], 0.25*win[0]*fin[95]*subvol); 
  atomicAdd(&out[98], 0.25*win[0]*fin[98]*subvol); 
  atomicAdd(&out[99], 0.25*win[0]*fin[99]*subvol); 
  atomicAdd(&out[101], 0.25*win[0]*fin[101]*subvol); 
  atomicAdd(&out[102], 0.25*win[0]*fin[102]*subvol); 
  atomicAdd(&out[104], 0.25*win[0]*fin[104]*subvol); 
  atomicAdd(&out[105], 0.25*win[0]*fin[105]*subvol); 
  atomicAdd(&out[106], 0.25*win[0]*fin[106]*subvol); 
  atomicAdd(&out[107], 0.25*win[0]*fin[107]*subvol); 
  atomicAdd(&out[111], 0.25*win[0]*fin[111]*subvol); 
  atomicAdd(&out[115], 0.25*win[0]*fin[115]*subvol); 
  atomicAdd(&out[119], 0.25*win[0]*fin[119]*subvol); 
  atomicAdd(&out[122], 0.25*win[0]*fin[122]*subvol); 
  atomicAdd(&out[125], 0.25*win[0]*fin[125]*subvol); 
  atomicAdd(&out[126], 0.25*win[0]*fin[126]*subvol); 
  atomicAdd(&out[128], 0.25*win[0]*fin[128]*subvol); 
  atomicAdd(&out[131], 0.25*win[0]*fin[131]*subvol); 
  atomicAdd(&out[132], 0.25*win[0]*fin[132]*subvol); 
  atomicAdd(&out[135], 0.25*win[0]*fin[135]*subvol); 
  atomicAdd(&out[138], 0.25*win[0]*fin[138]*subvol); 
  atomicAdd(&out[141], 0.25*win[0]*fin[141]*subvol); 
  atomicAdd(&out[142], 0.25*win[0]*fin[142]*subvol); 
  atomicAdd(&out[144], 0.25*win[0]*fin[144]*subvol); 
  atomicAdd(&out[145], 0.25*win[0]*fin[145]*subvol); 
  atomicAdd(&out[147], 0.25*win[0]*fin[147]*subvol); 
  atomicAdd(&out[148], 0.25*win[0]*fin[148]*subvol); 
  atomicAdd(&out[149], 0.25*win[0]*fin[149]*subvol); 
  atomicAdd(&out[151], 0.25*win[0]*fin[151]*subvol); 
  atomicAdd(&out[152], 0.25*win[0]*fin[152]*subvol); 
  atomicAdd(&out[153], 0.25*win[0]*fin[153]*subvol); 
  atomicAdd(&out[155], 0.25*win[0]*fin[155]*subvol); 
  atomicAdd(&out[160], 0.25*win[0]*fin[160]*subvol); 
  atomicAdd(&out[164], 0.25*win[0]*fin[164]*subvol); 
  atomicAdd(&out[169], 0.25*win[0]*fin[169]*subvol); 
  atomicAdd(&out[172], 0.25*win[0]*fin[172]*subvol); 
  atomicAdd(&out[175], 0.25*win[0]*fin[175]*subvol); 
  atomicAdd(&out[176], 0.25*win[0]*fin[176]*subvol); 
  atomicAdd(&out[178], 0.25*win[0]*fin[178]*subvol); 
  atomicAdd(&out[181], 0.25*win[0]*fin[181]*subvol); 
  atomicAdd(&out[182], 0.25*win[0]*fin[182]*subvol); 
  atomicAdd(&out[184], 0.25*win[0]*fin[184]*subvol); 
  atomicAdd(&out[187], 0.25*win[0]*fin[187]*subvol); 
  atomicAdd(&out[188], 0.25*win[0]*fin[188]*subvol); 
  atomicAdd(&out[190], 0.25*win[0]*fin[190]*subvol); 
  atomicAdd(&out[195], 0.25*win[0]*fin[195]*subvol); 
  atomicAdd(&out[199], 0.25*win[0]*fin[199]*subvol); 
  atomicAdd(&out[203], 0.25*win[0]*fin[203]*subvol); 
  atomicAdd(&out[206], 0.25*win[0]*fin[206]*subvol); 
  atomicAdd(&out[207], 0.25*win[0]*fin[207]*subvol); 
  atomicAdd(&out[212], 0.25*win[0]*fin[212]*subvol); 
#else
  out[0] += 0.25*fin[0]*win[0]*subvol; 
  out[1] += 0.25*win[0]*fin[1]*subvol; 
  out[2] += 0.25*win[0]*fin[2]*subvol; 
  out[4] += 0.25*win[0]*fin[4]*subvol; 
  out[5] += 0.25*win[0]*fin[5]*subvol; 
  out[6] += 0.25*win[0]*fin[6]*subvol; 
  out[7] += 0.25*win[0]*fin[7]*subvol; 
  out[10] += 0.25*win[0]*fin[10]*subvol; 
  out[11] += 0.25*win[0]*fin[11]*subvol; 
  out[13] += 0.25*win[0]*fin[13]*subvol; 
  out[14] += 0.25*win[0]*fin[14]*subvol; 
  out[16] += 0.25*win[0]*fin[16]*subvol; 
  out[17] += 0.25*win[0]*fin[17]*subvol; 
  out[18] += 0.25*win[0]*fin[18]*subvol; 
  out[20] += 0.25*win[0]*fin[20]*subvol; 
  out[21] += 0.25*win[0]*fin[21]*subvol; 
  out[22] += 0.25*win[0]*fin[22]*subvol; 
  out[23] += 0.25*win[0]*fin[23]*subvol; 
  out[24] += 0.25*win[0]*fin[24]*subvol; 
  out[26] += 0.25*win[0]*fin[26]*subvol; 
  out[29] += 0.25*win[0]*fin[29]*subvol; 
  out[32] += 0.25*win[0]*fin[32]*subvol; 
  out[33] += 0.25*win[0]*fin[33]*subvol; 
  out[35] += 0.25*win[0]*fin[35]*subvol; 
  out[38] += 0.25*win[0]*fin[38]*subvol; 
  out[39] += 0.25*win[0]*fin[39]*subvol; 
  out[41] += 0.25*win[0]*fin[41]*subvol; 
  out[42] += 0.25*win[0]*fin[42]*subvol; 
  out[44] += 0.25*win[0]*fin[44]*subvol; 
  out[45] += 0.25*win[0]*fin[45]*subvol; 
  out[46] += 0.25*win[0]*fin[46]*subvol; 
  out[48] += 0.25*win[0]*fin[48]*subvol; 
  out[49] += 0.25*win[0]*fin[49]*subvol; 
  out[50] += 0.25*win[0]*fin[50]*subvol; 
  out[52] += 0.25*win[0]*fin[52]*subvol; 
  out[53] += 0.25*win[0]*fin[53]*subvol; 
  out[54] += 0.25*win[0]*fin[54]*subvol; 
  out[55] += 0.25*win[0]*fin[55]*subvol; 
  out[56] += 0.25*win[0]*fin[56]*subvol; 
  out[58] += 0.25*win[0]*fin[58]*subvol; 
  out[59] += 0.25*win[0]*fin[59]*subvol; 
  out[62] += 0.25*win[0]*fin[62]*subvol; 
  out[66] += 0.25*win[0]*fin[66]*subvol; 
  out[69] += 0.25*win[0]*fin[69]*subvol; 
  out[72] += 0.25*win[0]*fin[72]*subvol; 
  out[73] += 0.25*win[0]*fin[73]*subvol; 
  out[75] += 0.25*win[0]*fin[75]*subvol; 
  out[78] += 0.25*win[0]*fin[78]*subvol; 
  out[79] += 0.25*win[0]*fin[79]*subvol; 
  out[81] += 0.25*win[0]*fin[81]*subvol; 
  out[84] += 0.25*win[0]*fin[84]*subvol; 
  out[85] += 0.25*win[0]*fin[85]*subvol; 
  out[87] += 0.25*win[0]*fin[87]*subvol; 
  out[88] += 0.25*win[0]*fin[88]*subvol; 
  out[90] += 0.25*win[0]*fin[90]*subvol; 
  out[91] += 0.25*win[0]*fin[91]*subvol; 
  out[92] += 0.25*win[0]*fin[92]*subvol; 
  out[94] += 0.25*win[0]*fin[94]*subvol; 
  out[95] += 0.25*win[0]*fin[95]*subvol; 
  out[98] += 0.25*win[0]*fin[98]*subvol; 
  out[99] += 0.25*win[0]*fin[99]*subvol; 
  out[101] += 0.25*win[0]*fin[101]*subvol; 
  out[102] += 0.25*win[0]*fin[102]*subvol; 
  out[104] += 0.25*win[0]*fin[104]*subvol; 
  out[105] += 0.25*win[0]*fin[105]*subvol; 
  out[106] += 0.25*win[0]*fin[106]*subvol; 
  out[107] += 0.25*win[0]*fin[107]*subvol; 
  out[111] += 0.25*win[0]*fin[111]*subvol; 
  out[115] += 0.25*win[0]*fin[115]*subvol; 
  out[119] += 0.25*win[0]*fin[119]*subvol; 
  out[122] += 0.25*win[0]*fin[122]*subvol; 
  out[125] += 0.25*win[0]*fin[125]*subvol; 
  out[126] += 0.25*win[0]*fin[126]*subvol; 
  out[128] += 0.25*win[0]*fin[128]*subvol; 
  out[131] += 0.25*win[0]*fin[131]*subvol; 
  out[132] += 0.25*win[0]*fin[132]*subvol; 
  out[135] += 0.25*win[0]*fin[135]*subvol; 
  out[138] += 0.25*win[0]*fin[138]*subvol; 
  out[141] += 0.25*win[0]*fin[141]*subvol; 
  out[142] += 0.25*win[0]*fin[142]*subvol; 
  out[144] += 0.25*win[0]*fin[144]*subvol; 
  out[145] += 0.25*win[0]*fin[145]*subvol; 
  out[147] += 0.25*win[0]*fin[147]*subvol; 
  out[148] += 0.25*win[0]*fin[148]*subvol; 
  out[149] += 0.25*win[0]*fin[149]*subvol; 
  out[151] += 0.25*win[0]*fin[151]*subvol; 
  out[152] += 0.25*win[0]*fin[152]*subvol; 
  out[153] += 0.25*win[0]*fin[153]*subvol; 
  out[155] += 0.25*win[0]*fin[155]*subvol; 
  out[160] += 0.25*win[0]*fin[160]*subvol; 
  out[164] += 0.25*win[0]*fin[164]*subvol; 
  out[169] += 0.25*win[0]*fin[169]*subvol; 
  out[172] += 0.25*win[0]*fin[172]*subvol; 
  out[175] += 0.25*win[0]*fin[175]*subvol; 
  out[176] += 0.25*win[0]*fin[176]*subvol; 
  out[178] += 0.25*win[0]*fin[178]*subvol; 
  out[181] += 0.25*win[0]*fin[181]*subvol; 
  out[182] += 0.25*win[0]*fin[182]*subvol; 
  out[184] += 0.25*win[0]*fin[184]*subvol; 
  out[187] += 0.25*win[0]*fin[187]*subvol; 
  out[188] += 0.25*win[0]*fin[188]*subvol; 
  out[190] += 0.25*win[0]*fin[190]*subvol; 
  out[195] += 0.25*win[0]*fin[195]*subvol; 
  out[199] += 0.25*win[0]*fin[199]*subvol; 
  out[203] += 0.25*win[0]*fin[203]*subvol; 
  out[206] += 0.25*win[0]*fin[206]*subvol; 
  out[207] += 0.25*win[0]*fin[207]*subvol; 
  out[212] += 0.25*win[0]*fin[212]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.25*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.25*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.25*win[0]*fin[4]*subvol); 
  atomicAdd(&out[3], 0.25*win[0]*fin[5]*subvol); 
  atomicAdd(&out[4], 0.25*win[0]*fin[6]*subvol); 
  atomicAdd(&out[5], 0.25*win[0]*fin[11]*subvol); 
  atomicAdd(&out[6], 0.25*win[0]*fin[14]*subvol); 
  atomicAdd(&out[7], 0.25*win[0]*fin[16]*subvol); 
  atomicAdd(&out[8], 0.25*win[0]*fin[18]*subvol); 
  atomicAdd(&out[9], 0.25*win[0]*fin[20]*subvol); 
  atomicAdd(&out[10], 0.25*win[0]*fin[21]*subvol); 
  atomicAdd(&out[11], 0.25*win[0]*fin[22]*subvol); 
  atomicAdd(&out[12], 0.25*win[0]*fin[23]*subvol); 
  atomicAdd(&out[13], 0.25*win[0]*fin[24]*subvol); 
  atomicAdd(&out[14], 0.25*win[0]*fin[33]*subvol); 
  atomicAdd(&out[15], 0.25*win[0]*fin[39]*subvol); 
  atomicAdd(&out[16], 0.25*win[0]*fin[42]*subvol); 
  atomicAdd(&out[17], 0.25*win[0]*fin[44]*subvol); 
  atomicAdd(&out[18], 0.25*win[0]*fin[46]*subvol); 
  atomicAdd(&out[19], 0.25*win[0]*fin[48]*subvol); 
  atomicAdd(&out[20], 0.25*win[0]*fin[50]*subvol); 
  atomicAdd(&out[21], 0.25*win[0]*fin[52]*subvol); 
  atomicAdd(&out[22], 0.25*win[0]*fin[53]*subvol); 
  atomicAdd(&out[23], 0.25*win[0]*fin[54]*subvol); 
  atomicAdd(&out[24], 0.25*win[0]*fin[56]*subvol); 
  atomicAdd(&out[25], 0.25*win[0]*fin[58]*subvol); 
  atomicAdd(&out[26], 0.25*win[0]*fin[59]*subvol); 
  atomicAdd(&out[27], 0.25*win[0]*fin[73]*subvol); 
  atomicAdd(&out[28], 0.25*win[0]*fin[79]*subvol); 
  atomicAdd(&out[29], 0.25*win[0]*fin[85]*subvol); 
  atomicAdd(&out[30], 0.25*win[0]*fin[88]*subvol); 
  atomicAdd(&out[31], 0.25*win[0]*fin[90]*subvol); 
  atomicAdd(&out[32], 0.25*win[0]*fin[92]*subvol); 
  atomicAdd(&out[33], 0.25*win[0]*fin[94]*subvol); 
  atomicAdd(&out[34], 0.25*win[0]*fin[99]*subvol); 
  atomicAdd(&out[35], 0.25*win[0]*fin[102]*subvol); 
  atomicAdd(&out[36], 0.25*win[0]*fin[104]*subvol); 
  atomicAdd(&out[37], 0.25*win[0]*fin[105]*subvol); 
  atomicAdd(&out[38], 0.25*win[0]*fin[106]*subvol); 
  atomicAdd(&out[39], 0.25*win[0]*fin[107]*subvol); 
  atomicAdd(&out[40], 0.25*win[0]*fin[126]*subvol); 
  atomicAdd(&out[41], 0.25*win[0]*fin[132]*subvol); 
  atomicAdd(&out[42], 0.25*win[0]*fin[142]*subvol); 
  atomicAdd(&out[43], 0.25*win[0]*fin[145]*subvol); 
  atomicAdd(&out[44], 0.25*win[0]*fin[147]*subvol); 
  atomicAdd(&out[45], 0.25*win[0]*fin[149]*subvol); 
  atomicAdd(&out[46], 0.25*win[0]*fin[151]*subvol); 
  atomicAdd(&out[47], 0.25*win[0]*fin[153]*subvol); 
  atomicAdd(&out[48], 0.25*win[0]*fin[155]*subvol); 
  atomicAdd(&out[49], 0.25*win[0]*fin[176]*subvol); 
  atomicAdd(&out[50], 0.25*win[0]*fin[182]*subvol); 
  atomicAdd(&out[51], 0.25*win[0]*fin[188]*subvol); 
  atomicAdd(&out[52], 0.25*win[0]*fin[190]*subvol); 
  atomicAdd(&out[53], 0.25*win[0]*fin[207]*subvol); 
#else
  out[0] += 0.25*fin[0]*win[0]*subvol; 
  out[1] += 0.25*win[0]*fin[2]*subvol; 
  out[2] += 0.25*win[0]*fin[4]*subvol; 
  out[3] += 0.25*win[0]*fin[5]*subvol; 
  out[4] += 0.25*win[0]*fin[6]*subvol; 
  out[5] += 0.25*win[0]*fin[11]*subvol; 
  out[6] += 0.25*win[0]*fin[14]*subvol; 
  out[7] += 0.25*win[0]*fin[16]*subvol; 
  out[8] += 0.25*win[0]*fin[18]*subvol; 
  out[9] += 0.25*win[0]*fin[20]*subvol; 
  out[10] += 0.25*win[0]*fin[21]*subvol; 
  out[11] += 0.25*win[0]*fin[22]*subvol; 
  out[12] += 0.25*win[0]*fin[23]*subvol; 
  out[13] += 0.25*win[0]*fin[24]*subvol; 
  out[14] += 0.25*win[0]*fin[33]*subvol; 
  out[15] += 0.25*win[0]*fin[39]*subvol; 
  out[16] += 0.25*win[0]*fin[42]*subvol; 
  out[17] += 0.25*win[0]*fin[44]*subvol; 
  out[18] += 0.25*win[0]*fin[46]*subvol; 
  out[19] += 0.25*win[0]*fin[48]*subvol; 
  out[20] += 0.25*win[0]*fin[50]*subvol; 
  out[21] += 0.25*win[0]*fin[52]*subvol; 
  out[22] += 0.25*win[0]*fin[53]*subvol; 
  out[23] += 0.25*win[0]*fin[54]*subvol; 
  out[24] += 0.25*win[0]*fin[56]*subvol; 
  out[25] += 0.25*win[0]*fin[58]*subvol; 
  out[26] += 0.25*win[0]*fin[59]*subvol; 
  out[27] += 0.25*win[0]*fin[73]*subvol; 
  out[28] += 0.25*win[0]*fin[79]*subvol; 
  out[29] += 0.25*win[0]*fin[85]*subvol; 
  out[30] += 0.25*win[0]*fin[88]*subvol; 
  out[31] += 0.25*win[0]*fin[90]*subvol; 
  out[32] += 0.25*win[0]*fin[92]*subvol; 
  out[33] += 0.25*win[0]*fin[94]*subvol; 
  out[34] += 0.25*win[0]*fin[99]*subvol; 
  out[35] += 0.25*win[0]*fin[102]*subvol; 
  out[36] += 0.25*win[0]*fin[104]*subvol; 
  out[37] += 0.25*win[0]*fin[105]*subvol; 
  out[38] += 0.25*win[0]*fin[106]*subvol; 
  out[39] += 0.25*win[0]*fin[107]*subvol; 
  out[40] += 0.25*win[0]*fin[126]*subvol; 
  out[41] += 0.25*win[0]*fin[132]*subvol; 
  out[42] += 0.25*win[0]*fin[142]*subvol; 
  out[43] += 0.25*win[0]*fin[145]*subvol; 
  out[44] += 0.25*win[0]*fin[147]*subvol; 
  out[45] += 0.25*win[0]*fin[149]*subvol; 
  out[46] += 0.25*win[0]*fin[151]*subvol; 
  out[47] += 0.25*win[0]*fin[153]*subvol; 
  out[48] += 0.25*win[0]*fin[155]*subvol; 
  out[49] += 0.25*win[0]*fin[176]*subvol; 
  out[50] += 0.25*win[0]*fin[182]*subvol; 
  out[51] += 0.25*win[0]*fin[188]*subvol; 
  out[52] += 0.25*win[0]*fin[190]*subvol; 
  out[53] += 0.25*win[0]*fin[207]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[5]*subvol); 
  atomicAdd(&out[6], 0.5*win[0]*fin[6]*subvol); 
  atomicAdd(&out[11], 0.5*win[0]*fin[11]*subvol); 
  atomicAdd(&out[14], 0.5*win[0]*fin[14]*subvol); 
  atomicAdd(&out[16], 0.5*win[0]*fin[16]*subvol); 
  atomicAdd(&out[18], 0.5*win[0]*fin[18]*subvol); 
  atomicAdd(&out[20], 0.5*win[0]*fin[20]*subvol); 
  atomicAdd(&out[21], 0.5*win[0]*fin[21]*subvol); 
  atomicAdd(&out[22], 0.5*win[0]*fin[22]*subvol); 
  atomicAdd(&out[23], 0.5*win[0]*fin[23]*subvol); 
  atomicAdd(&out[24], 0.5*win[0]*fin[24]*subvol); 
  atomicAdd(&out[33], 0.5*win[0]*fin[33]*subvol); 
  atomicAdd(&out[39], 0.5*win[0]*fin[39]*subvol); 
  atomicAdd(&out[42], 0.5*win[0]*fin[42]*subvol); 
  atomicAdd(&out[44], 0.5*win[0]*fin[44]*subvol); 
  atomicAdd(&out[46], 0.5*win[0]*fin[46]*subvol); 
  atomicAdd(&out[48], 0.5*win[0]*fin[48]*subvol); 
  atomicAdd(&out[50], 0.5*win[0]*fin[50]*subvol); 
  atomicAdd(&out[52], 0.5*win[0]*fin[52]*subvol); 
  atomicAdd(&out[53], 0.5*win[0]*fin[53]*subvol); 
  atomicAdd(&out[54], 0.5*win[0]*fin[54]*subvol); 
  atomicAdd(&out[56], 0.5*win[0]*fin[56]*subvol); 
  atomicAdd(&out[58], 0.5*win[0]*fin[58]*subvol); 
  atomicAdd(&out[59], 0.5*win[0]*fin[59]*subvol); 
  atomicAdd(&out[73], 0.5*win[0]*fin[73]*subvol); 
  atomicAdd(&out[79], 0.5*win[0]*fin[79]*subvol); 
  atomicAdd(&out[85], 0.5*win[0]*fin[85]*subvol); 
  atomicAdd(&out[88], 0.5*win[0]*fin[88]*subvol); 
  atomicAdd(&out[90], 0.5*win[0]*fin[90]*subvol); 
  atomicAdd(&out[92], 0.5*win[0]*fin[92]*subvol); 
  atomicAdd(&out[94], 0.5*win[0]*fin[94]*subvol); 
  atomicAdd(&out[99], 0.5*win[0]*fin[99]*subvol); 
  atomicAdd(&out[102], 0.5*win[0]*fin[102]*subvol); 
  atomicAdd(&out[104], 0.5*win[0]*fin[104]*subvol); 
  atomicAdd(&out[105], 0.5*win[0]*fin[105]*subvol); 
  atomicAdd(&out[106], 0.5*win[0]*fin[106]*subvol); 
  atomicAdd(&out[107], 0.5*win[0]*fin[107]*subvol); 
  atomicAdd(&out[126], 0.5*win[0]*fin[126]*subvol); 
  atomicAdd(&out[132], 0.5*win[0]*fin[132]*subvol); 
  atomicAdd(&out[142], 0.5*win[0]*fin[142]*subvol); 
  atomicAdd(&out[145], 0.5*win[0]*fin[145]*subvol); 
  atomicAdd(&out[147], 0.5*win[0]*fin[147]*subvol); 
  atomicAdd(&out[149], 0.5*win[0]*fin[149]*subvol); 
  atomicAdd(&out[151], 0.5*win[0]*fin[151]*subvol); 
  atomicAdd(&out[153], 0.5*win[0]*fin[153]*subvol); 
  atomicAdd(&out[155], 0.5*win[0]*fin[155]*subvol); 
  atomicAdd(&out[176], 0.5*win[0]*fin[176]*subvol); 
  atomicAdd(&out[182], 0.5*win[0]*fin[182]*subvol); 
  atomicAdd(&out[188], 0.5*win[0]*fin[188]*subvol); 
  atomicAdd(&out[190], 0.5*win[0]*fin[190]*subvol); 
  atomicAdd(&out[207], 0.5*win[0]*fin[207]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[2] += 0.5*win[0]*fin[2]*subvol; 
  out[4] += 0.5*win[0]*fin[4]*subvol; 
  out[5] += 0.5*win[0]*fin[5]*subvol; 
  out[6] += 0.5*win[0]*fin[6]*subvol; 
  out[11] += 0.5*win[0]*fin[11]*subvol; 
  out[14] += 0.5*win[0]*fin[14]*subvol; 
  out[16] += 0.5*win[0]*fin[16]*subvol; 
  out[18] += 0.5*win[0]*fin[18]*subvol; 
  out[20] += 0.5*win[0]*fin[20]*subvol; 
  out[21] += 0.5*win[0]*fin[21]*subvol; 
  out[22] += 0.5*win[0]*fin[22]*subvol; 
  out[23] += 0.5*win[0]*fin[23]*subvol; 
  out[24] += 0.5*win[0]*fin[24]*subvol; 
  out[33] += 0.5*win[0]*fin[33]*subvol; 
  out[39] += 0.5*win[0]*fin[39]*subvol; 
  out[42] += 0.5*win[0]*fin[42]*subvol; 
  out[44] += 0.5*win[0]*fin[44]*subvol; 
  out[46] += 0.5*win[0]*fin[46]*subvol; 
  out[48] += 0.5*win[0]*fin[48]*subvol; 
  out[50] += 0.5*win[0]*fin[50]*subvol; 
  out[52] += 0.5*win[0]*fin[52]*subvol; 
  out[53] += 0.5*win[0]*fin[53]*subvol; 
  out[54] += 0.5*win[0]*fin[54]*subvol; 
  out[56] += 0.5*win[0]*fin[56]*subvol; 
  out[58] += 0.5*win[0]*fin[58]*subvol; 
  out[59] += 0.5*win[0]*fin[59]*subvol; 
  out[73] += 0.5*win[0]*fin[73]*subvol; 
  out[79] += 0.5*win[0]*fin[79]*subvol; 
  out[85] += 0.5*win[0]*fin[85]*subvol; 
  out[88] += 0.5*win[0]*fin[88]*subvol; 
  out[90] += 0.5*win[0]*fin[90]*subvol; 
  out[92] += 0.5*win[0]*fin[92]*subvol; 
  out[94] += 0.5*win[0]*fin[94]*subvol; 
  out[99] += 0.5*win[0]*fin[99]*subvol; 
  out[102] += 0.5*win[0]*fin[102]*subvol; 
  out[104] += 0.5*win[0]*fin[104]*subvol; 
  out[105] += 0.5*win[0]*fin[105]*subvol; 
  out[106] += 0.5*win[0]*fin[106]*subvol; 
  out[107] += 0.5*win[0]*fin[107]*subvol; 
  out[126] += 0.5*win[0]*fin[126]*subvol; 
  out[132] += 0.5*win[0]*fin[132]*subvol; 
  out[142] += 0.5*win[0]*fin[142]*subvol; 
  out[145] += 0.5*win[0]*fin[145]*subvol; 
  out[147] += 0.5*win[0]*fin[147]*subvol; 
  out[149] += 0.5*win[0]*fin[149]*subvol; 
  out[151] += 0.5*win[0]*fin[151]*subvol; 
  out[153] += 0.5*win[0]*fin[153]*subvol; 
  out[155] += 0.5*win[0]*fin[155]*subvol; 
  out[176] += 0.5*win[0]*fin[176]*subvol; 
  out[182] += 0.5*win[0]*fin[182]*subvol; 
  out[188] += 0.5*win[0]*fin[188]*subvol; 
  out[190] += 0.5*win[0]*fin[190]*subvol; 
  out[207] += 0.5*win[0]*fin[207]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgyz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.25*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.25*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.25*win[0]*fin[4]*subvol); 
  atomicAdd(&out[3], 0.25*win[0]*fin[5]*subvol); 
  atomicAdd(&out[4], 0.25*win[0]*fin[6]*subvol); 
  atomicAdd(&out[5], 0.25*win[0]*fin[10]*subvol); 
  atomicAdd(&out[6], 0.25*win[0]*fin[13]*subvol); 
  atomicAdd(&out[7], 0.25*win[0]*fin[16]*subvol); 
  atomicAdd(&out[8], 0.25*win[0]*fin[17]*subvol); 
  atomicAdd(&out[9], 0.25*win[0]*fin[20]*subvol); 
  atomicAdd(&out[10], 0.25*win[0]*fin[21]*subvol); 
  atomicAdd(&out[11], 0.25*win[0]*fin[22]*subvol); 
  atomicAdd(&out[12], 0.25*win[0]*fin[23]*subvol); 
  atomicAdd(&out[13], 0.25*win[0]*fin[24]*subvol); 
  atomicAdd(&out[14], 0.25*win[0]*fin[32]*subvol); 
  atomicAdd(&out[15], 0.25*win[0]*fin[38]*subvol); 
  atomicAdd(&out[16], 0.25*win[0]*fin[41]*subvol); 
  atomicAdd(&out[17], 0.25*win[0]*fin[44]*subvol); 
  atomicAdd(&out[18], 0.25*win[0]*fin[45]*subvol); 
  atomicAdd(&out[19], 0.25*win[0]*fin[48]*subvol); 
  atomicAdd(&out[20], 0.25*win[0]*fin[49]*subvol); 
  atomicAdd(&out[21], 0.25*win[0]*fin[52]*subvol); 
  atomicAdd(&out[22], 0.25*win[0]*fin[53]*subvol); 
  atomicAdd(&out[23], 0.25*win[0]*fin[54]*subvol); 
  atomicAdd(&out[24], 0.25*win[0]*fin[55]*subvol); 
  atomicAdd(&out[25], 0.25*win[0]*fin[58]*subvol); 
  atomicAdd(&out[26], 0.25*win[0]*fin[59]*subvol); 
  atomicAdd(&out[27], 0.25*win[0]*fin[72]*subvol); 
  atomicAdd(&out[28], 0.25*win[0]*fin[78]*subvol); 
  atomicAdd(&out[29], 0.25*win[0]*fin[84]*subvol); 
  atomicAdd(&out[30], 0.25*win[0]*fin[87]*subvol); 
  atomicAdd(&out[31], 0.25*win[0]*fin[90]*subvol); 
  atomicAdd(&out[32], 0.25*win[0]*fin[91]*subvol); 
  atomicAdd(&out[33], 0.25*win[0]*fin[94]*subvol); 
  atomicAdd(&out[34], 0.25*win[0]*fin[98]*subvol); 
  atomicAdd(&out[35], 0.25*win[0]*fin[101]*subvol); 
  atomicAdd(&out[36], 0.25*win[0]*fin[104]*subvol); 
  atomicAdd(&out[37], 0.25*win[0]*fin[105]*subvol); 
  atomicAdd(&out[38], 0.25*win[0]*fin[106]*subvol); 
  atomicAdd(&out[39], 0.25*win[0]*fin[107]*subvol); 
  atomicAdd(&out[40], 0.25*win[0]*fin[125]*subvol); 
  atomicAdd(&out[41], 0.25*win[0]*fin[131]*subvol); 
  atomicAdd(&out[42], 0.25*win[0]*fin[141]*subvol); 
  atomicAdd(&out[43], 0.25*win[0]*fin[144]*subvol); 
  atomicAdd(&out[44], 0.25*win[0]*fin[147]*subvol); 
  atomicAdd(&out[45], 0.25*win[0]*fin[148]*subvol); 
  atomicAdd(&out[46], 0.25*win[0]*fin[151]*subvol); 
  atomicAdd(&out[47], 0.25*win[0]*fin[152]*subvol); 
  atomicAdd(&out[48], 0.25*win[0]*fin[155]*subvol); 
  atomicAdd(&out[49], 0.25*win[0]*fin[175]*subvol); 
  atomicAdd(&out[50], 0.25*win[0]*fin[181]*subvol); 
  atomicAdd(&out[51], 0.25*win[0]*fin[187]*subvol); 
  atomicAdd(&out[52], 0.25*win[0]*fin[190]*subvol); 
  atomicAdd(&out[53], 0.25*win[0]*fin[206]*subvol); 
#else
  out[0] += 0.25*fin[0]*win[0]*subvol; 
  out[1] += 0.25*win[0]*fin[1]*subvol; 
  out[2] += 0.25*win[0]*fin[4]*subvol; 
  out[3] += 0.25*win[0]*fin[5]*subvol; 
  out[4] += 0.25*win[0]*fin[6]*subvol; 
  out[5] += 0.25*win[0]*fin[10]*subvol; 
  out[6] += 0.25*win[0]*fin[13]*subvol; 
  out[7] += 0.25*win[0]*fin[16]*subvol; 
  out[8] += 0.25*win[0]*fin[17]*subvol; 
  out[9] += 0.25*win[0]*fin[20]*subvol; 
  out[10] += 0.25*win[0]*fin[21]*subvol; 
  out[11] += 0.25*win[0]*fin[22]*subvol; 
  out[12] += 0.25*win[0]*fin[23]*subvol; 
  out[13] += 0.25*win[0]*fin[24]*subvol; 
  out[14] += 0.25*win[0]*fin[32]*subvol; 
  out[15] += 0.25*win[0]*fin[38]*subvol; 
  out[16] += 0.25*win[0]*fin[41]*subvol; 
  out[17] += 0.25*win[0]*fin[44]*subvol; 
  out[18] += 0.25*win[0]*fin[45]*subvol; 
  out[19] += 0.25*win[0]*fin[48]*subvol; 
  out[20] += 0.25*win[0]*fin[49]*subvol; 
  out[21] += 0.25*win[0]*fin[52]*subvol; 
  out[22] += 0.25*win[0]*fin[53]*subvol; 
  out[23] += 0.25*win[0]*fin[54]*subvol; 
  out[24] += 0.25*win[0]*fin[55]*subvol; 
  out[25] += 0.25*win[0]*fin[58]*subvol; 
  out[26] += 0.25*win[0]*fin[59]*subvol; 
  out[27] += 0.25*win[0]*fin[72]*subvol; 
  out[28] += 0.25*win[0]*fin[78]*subvol; 
  out[29] += 0.25*win[0]*fin[84]*subvol; 
  out[30] += 0.25*win[0]*fin[87]*subvol; 
  out[31] += 0.25*win[0]*fin[90]*subvol; 
  out[32] += 0.25*win[0]*fin[91]*subvol; 
  out[33] += 0.25*win[0]*fin[94]*subvol; 
  out[34] += 0.25*win[0]*fin[98]*subvol; 
  out[35] += 0.25*win[0]*fin[101]*subvol; 
  out[36] += 0.25*win[0]*fin[104]*subvol; 
  out[37] += 0.25*win[0]*fin[105]*subvol; 
  out[38] += 0.25*win[0]*fin[106]*subvol; 
  out[39] += 0.25*win[0]*fin[107]*subvol; 
  out[40] += 0.25*win[0]*fin[125]*subvol; 
  out[41] += 0.25*win[0]*fin[131]*subvol; 
  out[42] += 0.25*win[0]*fin[141]*subvol; 
  out[43] += 0.25*win[0]*fin[144]*subvol; 
  out[44] += 0.25*win[0]*fin[147]*subvol; 
  out[45] += 0.25*win[0]*fin[148]*subvol; 
  out[46] += 0.25*win[0]*fin[151]*subvol; 
  out[47] += 0.25*win[0]*fin[152]*subvol; 
  out[48] += 0.25*win[0]*fin[155]*subvol; 
  out[49] += 0.25*win[0]*fin[175]*subvol; 
  out[50] += 0.25*win[0]*fin[181]*subvol; 
  out[51] += 0.25*win[0]*fin[187]*subvol; 
  out[52] += 0.25*win[0]*fin[190]*subvol; 
  out[53] += 0.25*win[0]*fin[206]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgyz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[5]*subvol); 
  atomicAdd(&out[6], 0.5*win[0]*fin[6]*subvol); 
  atomicAdd(&out[10], 0.5*win[0]*fin[10]*subvol); 
  atomicAdd(&out[13], 0.5*win[0]*fin[13]*subvol); 
  atomicAdd(&out[16], 0.5*win[0]*fin[16]*subvol); 
  atomicAdd(&out[17], 0.5*win[0]*fin[17]*subvol); 
  atomicAdd(&out[20], 0.5*win[0]*fin[20]*subvol); 
  atomicAdd(&out[21], 0.5*win[0]*fin[21]*subvol); 
  atomicAdd(&out[22], 0.5*win[0]*fin[22]*subvol); 
  atomicAdd(&out[23], 0.5*win[0]*fin[23]*subvol); 
  atomicAdd(&out[24], 0.5*win[0]*fin[24]*subvol); 
  atomicAdd(&out[32], 0.5*win[0]*fin[32]*subvol); 
  atomicAdd(&out[38], 0.5*win[0]*fin[38]*subvol); 
  atomicAdd(&out[41], 0.5*win[0]*fin[41]*subvol); 
  atomicAdd(&out[44], 0.5*win[0]*fin[44]*subvol); 
  atomicAdd(&out[45], 0.5*win[0]*fin[45]*subvol); 
  atomicAdd(&out[48], 0.5*win[0]*fin[48]*subvol); 
  atomicAdd(&out[49], 0.5*win[0]*fin[49]*subvol); 
  atomicAdd(&out[52], 0.5*win[0]*fin[52]*subvol); 
  atomicAdd(&out[53], 0.5*win[0]*fin[53]*subvol); 
  atomicAdd(&out[54], 0.5*win[0]*fin[54]*subvol); 
  atomicAdd(&out[55], 0.5*win[0]*fin[55]*subvol); 
  atomicAdd(&out[58], 0.5*win[0]*fin[58]*subvol); 
  atomicAdd(&out[59], 0.5*win[0]*fin[59]*subvol); 
  atomicAdd(&out[72], 0.5*win[0]*fin[72]*subvol); 
  atomicAdd(&out[78], 0.5*win[0]*fin[78]*subvol); 
  atomicAdd(&out[84], 0.5*win[0]*fin[84]*subvol); 
  atomicAdd(&out[87], 0.5*win[0]*fin[87]*subvol); 
  atomicAdd(&out[90], 0.5*win[0]*fin[90]*subvol); 
  atomicAdd(&out[91], 0.5*win[0]*fin[91]*subvol); 
  atomicAdd(&out[94], 0.5*win[0]*fin[94]*subvol); 
  atomicAdd(&out[98], 0.5*win[0]*fin[98]*subvol); 
  atomicAdd(&out[101], 0.5*win[0]*fin[101]*subvol); 
  atomicAdd(&out[104], 0.5*win[0]*fin[104]*subvol); 
  atomicAdd(&out[105], 0.5*win[0]*fin[105]*subvol); 
  atomicAdd(&out[106], 0.5*win[0]*fin[106]*subvol); 
  atomicAdd(&out[107], 0.5*win[0]*fin[107]*subvol); 
  atomicAdd(&out[125], 0.5*win[0]*fin[125]*subvol); 
  atomicAdd(&out[131], 0.5*win[0]*fin[131]*subvol); 
  atomicAdd(&out[141], 0.5*win[0]*fin[141]*subvol); 
  atomicAdd(&out[144], 0.5*win[0]*fin[144]*subvol); 
  atomicAdd(&out[147], 0.5*win[0]*fin[147]*subvol); 
  atomicAdd(&out[148], 0.5*win[0]*fin[148]*subvol); 
  atomicAdd(&out[151], 0.5*win[0]*fin[151]*subvol); 
  atomicAdd(&out[152], 0.5*win[0]*fin[152]*subvol); 
  atomicAdd(&out[155], 0.5*win[0]*fin[155]*subvol); 
  atomicAdd(&out[175], 0.5*win[0]*fin[175]*subvol); 
  atomicAdd(&out[181], 0.5*win[0]*fin[181]*subvol); 
  atomicAdd(&out[187], 0.5*win[0]*fin[187]*subvol); 
  atomicAdd(&out[190], 0.5*win[0]*fin[190]*subvol); 
  atomicAdd(&out[206], 0.5*win[0]*fin[206]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[4] += 0.5*win[0]*fin[4]*subvol; 
  out[5] += 0.5*win[0]*fin[5]*subvol; 
  out[6] += 0.5*win[0]*fin[6]*subvol; 
  out[10] += 0.5*win[0]*fin[10]*subvol; 
  out[13] += 0.5*win[0]*fin[13]*subvol; 
  out[16] += 0.5*win[0]*fin[16]*subvol; 
  out[17] += 0.5*win[0]*fin[17]*subvol; 
  out[20] += 0.5*win[0]*fin[20]*subvol; 
  out[21] += 0.5*win[0]*fin[21]*subvol; 
  out[22] += 0.5*win[0]*fin[22]*subvol; 
  out[23] += 0.5*win[0]*fin[23]*subvol; 
  out[24] += 0.5*win[0]*fin[24]*subvol; 
  out[32] += 0.5*win[0]*fin[32]*subvol; 
  out[38] += 0.5*win[0]*fin[38]*subvol; 
  out[41] += 0.5*win[0]*fin[41]*subvol; 
  out[44] += 0.5*win[0]*fin[44]*subvol; 
  out[45] += 0.5*win[0]*fin[45]*subvol; 
  out[48] += 0.5*win[0]*fin[48]*subvol; 
  out[49] += 0.5*win[0]*fin[49]*subvol; 
  out[52] += 0.5*win[0]*fin[52]*subvol; 
  out[53] += 0.5*win[0]*fin[53]*subvol; 
  out[54] += 0.5*win[0]*fin[54]*subvol; 
  out[55] += 0.5*win[0]*fin[55]*subvol; 
  out[58] += 0.5*win[0]*fin[58]*subvol; 
  out[59] += 0.5*win[0]*fin[59]*subvol; 
  out[72] += 0.5*win[0]*fin[72]*subvol; 
  out[78] += 0.5*win[0]*fin[78]*subvol; 
  out[84] += 0.5*win[0]*fin[84]*subvol; 
  out[87] += 0.5*win[0]*fin[87]*subvol; 
  out[90] += 0.5*win[0]*fin[90]*subvol; 
  out[91] += 0.5*win[0]*fin[91]*subvol; 
  out[94] += 0.5*win[0]*fin[94]*subvol; 
  out[98] += 0.5*win[0]*fin[98]*subvol; 
  out[101] += 0.5*win[0]*fin[101]*subvol; 
  out[104] += 0.5*win[0]*fin[104]*subvol; 
  out[105] += 0.5*win[0]*fin[105]*subvol; 
  out[106] += 0.5*win[0]*fin[106]*subvol; 
  out[107] += 0.5*win[0]*fin[107]*subvol; 
  out[125] += 0.5*win[0]*fin[125]*subvol; 
  out[131] += 0.5*win[0]*fin[131]*subvol; 
  out[141] += 0.5*win[0]*fin[141]*subvol; 
  out[144] += 0.5*win[0]*fin[144]*subvol; 
  out[147] += 0.5*win[0]*fin[147]*subvol; 
  out[148] += 0.5*win[0]*fin[148]*subvol; 
  out[151] += 0.5*win[0]*fin[151]*subvol; 
  out[152] += 0.5*win[0]*fin[152]*subvol; 
  out[155] += 0.5*win[0]*fin[155]*subvol; 
  out[175] += 0.5*win[0]*fin[175]*subvol; 
  out[181] += 0.5*win[0]*fin[181]*subvol; 
  out[187] += 0.5*win[0]*fin[187]*subvol; 
  out[190] += 0.5*win[0]*fin[190]*subvol; 
  out[206] += 0.5*win[0]*fin[206]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxyz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[4]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[5]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[6]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[16]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[20]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[21]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[22]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[23]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[24]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[44]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[48]*subvol); 
  atomicAdd(&out[12], 0.3535533905932737*win[0]*fin[52]*subvol); 
  atomicAdd(&out[13], 0.3535533905932737*win[0]*fin[53]*subvol); 
  atomicAdd(&out[14], 0.3535533905932737*win[0]*fin[54]*subvol); 
  atomicAdd(&out[15], 0.3535533905932737*win[0]*fin[58]*subvol); 
  atomicAdd(&out[16], 0.3535533905932737*win[0]*fin[59]*subvol); 
  atomicAdd(&out[17], 0.3535533905932737*win[0]*fin[90]*subvol); 
  atomicAdd(&out[18], 0.3535533905932737*win[0]*fin[94]*subvol); 
  atomicAdd(&out[19], 0.3535533905932737*win[0]*fin[104]*subvol); 
  atomicAdd(&out[20], 0.3535533905932737*win[0]*fin[105]*subvol); 
  atomicAdd(&out[21], 0.3535533905932737*win[0]*fin[106]*subvol); 
  atomicAdd(&out[22], 0.3535533905932737*win[0]*fin[107]*subvol); 
  atomicAdd(&out[23], 0.3535533905932737*win[0]*fin[147]*subvol); 
  atomicAdd(&out[24], 0.3535533905932737*win[0]*fin[151]*subvol); 
  atomicAdd(&out[25], 0.3535533905932737*win[0]*fin[155]*subvol); 
  atomicAdd(&out[26], 0.3535533905932737*win[0]*fin[190]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[4]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[5]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[6]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[16]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[20]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[21]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[22]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[23]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[24]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[44]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[48]*subvol; 
  out[12] += 0.3535533905932737*win[0]*fin[52]*subvol; 
  out[13] += 0.3535533905932737*win[0]*fin[53]*subvol; 
  out[14] += 0.3535533905932737*win[0]*fin[54]*subvol; 
  out[15] += 0.3535533905932737*win[0]*fin[58]*subvol; 
  out[16] += 0.3535533905932737*win[0]*fin[59]*subvol; 
  out[17] += 0.3535533905932737*win[0]*fin[90]*subvol; 
  out[18] += 0.3535533905932737*win[0]*fin[94]*subvol; 
  out[19] += 0.3535533905932737*win[0]*fin[104]*subvol; 
  out[20] += 0.3535533905932737*win[0]*fin[105]*subvol; 
  out[21] += 0.3535533905932737*win[0]*fin[106]*subvol; 
  out[22] += 0.3535533905932737*win[0]*fin[107]*subvol; 
  out[23] += 0.3535533905932737*win[0]*fin[147]*subvol; 
  out[24] += 0.3535533905932737*win[0]*fin[151]*subvol; 
  out[25] += 0.3535533905932737*win[0]*fin[155]*subvol; 
  out[26] += 0.3535533905932737*win[0]*fin[190]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxyz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
  atomicAdd(&out[4], win[0]*fin[4]*subvol); 
  atomicAdd(&out[5], win[0]*fin[5]*subvol); 
  atomicAdd(&out[6], win[0]*fin[6]*subvol); 
  atomicAdd(&out[16], win[0]*fin[16]*subvol); 
  atomicAdd(&out[20], win[0]*fin[20]*subvol); 
  atomicAdd(&out[21], win[0]*fin[21]*subvol); 
  atomicAdd(&out[22], win[0]*fin[22]*subvol); 
  atomicAdd(&out[23], win[0]*fin[23]*subvol); 
  atomicAdd(&out[24], win[0]*fin[24]*subvol); 
  atomicAdd(&out[44], win[0]*fin[44]*subvol); 
  atomicAdd(&out[48], win[0]*fin[48]*subvol); 
  atomicAdd(&out[52], win[0]*fin[52]*subvol); 
  atomicAdd(&out[53], win[0]*fin[53]*subvol); 
  atomicAdd(&out[54], win[0]*fin[54]*subvol); 
  atomicAdd(&out[58], win[0]*fin[58]*subvol); 
  atomicAdd(&out[59], win[0]*fin[59]*subvol); 
  atomicAdd(&out[90], win[0]*fin[90]*subvol); 
  atomicAdd(&out[94], win[0]*fin[94]*subvol); 
  atomicAdd(&out[104], win[0]*fin[104]*subvol); 
  atomicAdd(&out[105], win[0]*fin[105]*subvol); 
  atomicAdd(&out[106], win[0]*fin[106]*subvol); 
  atomicAdd(&out[107], win[0]*fin[107]*subvol); 
  atomicAdd(&out[147], win[0]*fin[147]*subvol); 
  atomicAdd(&out[151], win[0]*fin[151]*subvol); 
  atomicAdd(&out[155], win[0]*fin[155]*subvol); 
  atomicAdd(&out[190], win[0]*fin[190]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
  out[4] += win[0]*fin[4]*subvol; 
  out[5] += win[0]*fin[5]*subvol; 
  out[6] += win[0]*fin[6]*subvol; 
  out[16] += win[0]*fin[16]*subvol; 
  out[20] += win[0]*fin[20]*subvol; 
  out[21] += win[0]*fin[21]*subvol; 
  out[22] += win[0]*fin[22]*subvol; 
  out[23] += win[0]*fin[23]*subvol; 
  out[24] += win[0]*fin[24]*subvol; 
  out[44] += win[0]*fin[44]*subvol; 
  out[48] += win[0]*fin[48]*subvol; 
  out[52] += win[0]*fin[52]*subvol; 
  out[53] += win[0]*fin[53]*subvol; 
  out[54] += win[0]*fin[54]*subvol; 
  out[58] += win[0]*fin[58]*subvol; 
  out[59] += win[0]*fin[59]*subvol; 
  out[90] += win[0]*fin[90]*subvol; 
  out[94] += win[0]*fin[94]*subvol; 
  out[104] += win[0]*fin[104]*subvol; 
  out[105] += win[0]*fin[105]*subvol; 
  out[106] += win[0]*fin[106]*subvol; 
  out[107] += win[0]*fin[107]*subvol; 
  out[147] += win[0]*fin[147]*subvol; 
  out[151] += win[0]*fin[151]*subvol; 
  out[155] += win[0]*fin[155]*subvol; 
  out[190] += win[0]*fin[190]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgvx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.1767766952966368*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.1767766952966368*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.1767766952966368*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.1767766952966368*win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], 0.1767766952966368*win[0]*fin[5]*subvol); 
  atomicAdd(&out[5], 0.1767766952966368*win[0]*fin[6]*subvol); 
  atomicAdd(&out[6], 0.1767766952966368*win[0]*fin[7]*subvol); 
  atomicAdd(&out[7], 0.1767766952966368*win[0]*fin[8]*subvol); 
  atomicAdd(&out[8], 0.1767766952966368*win[0]*fin[9]*subvol); 
  atomicAdd(&out[9], 0.1767766952966368*win[0]*fin[13]*subvol); 
  atomicAdd(&out[10], 0.1767766952966368*win[0]*fin[14]*subvol); 
  atomicAdd(&out[11], 0.1767766952966368*win[0]*fin[15]*subvol); 
  atomicAdd(&out[12], 0.1767766952966368*win[0]*fin[17]*subvol); 
  atomicAdd(&out[13], 0.1767766952966368*win[0]*fin[18]*subvol); 
  atomicAdd(&out[14], 0.1767766952966368*win[0]*fin[19]*subvol); 
  atomicAdd(&out[15], 0.1767766952966368*win[0]*fin[21]*subvol); 
  atomicAdd(&out[16], 0.1767766952966368*win[0]*fin[23]*subvol); 
  atomicAdd(&out[17], 0.1767766952966368*win[0]*fin[24]*subvol); 
  atomicAdd(&out[18], 0.1767766952966368*win[0]*fin[25]*subvol); 
  atomicAdd(&out[19], 0.1767766952966368*win[0]*fin[29]*subvol); 
  atomicAdd(&out[20], 0.1767766952966368*win[0]*fin[30]*subvol); 
  atomicAdd(&out[21], 0.1767766952966368*win[0]*fin[31]*subvol); 
  atomicAdd(&out[22], 0.1767766952966368*win[0]*fin[35]*subvol); 
  atomicAdd(&out[23], 0.1767766952966368*win[0]*fin[36]*subvol); 
  atomicAdd(&out[24], 0.1767766952966368*win[0]*fin[37]*subvol); 
  atomicAdd(&out[25], 0.1767766952966368*win[0]*fin[41]*subvol); 
  atomicAdd(&out[26], 0.1767766952966368*win[0]*fin[42]*subvol); 
  atomicAdd(&out[27], 0.1767766952966368*win[0]*fin[43]*subvol); 
  atomicAdd(&out[28], 0.1767766952966368*win[0]*fin[49]*subvol); 
  atomicAdd(&out[29], 0.1767766952966368*win[0]*fin[50]*subvol); 
  atomicAdd(&out[30], 0.1767766952966368*win[0]*fin[51]*subvol); 
  atomicAdd(&out[31], 0.1767766952966368*win[0]*fin[54]*subvol); 
  atomicAdd(&out[32], 0.1767766952966368*win[0]*fin[55]*subvol); 
  atomicAdd(&out[33], 0.1767766952966368*win[0]*fin[56]*subvol); 
  atomicAdd(&out[34], 0.1767766952966368*win[0]*fin[57]*subvol); 
  atomicAdd(&out[35], 0.1767766952966368*win[0]*fin[59]*subvol); 
  atomicAdd(&out[36], 0.1767766952966368*win[0]*fin[61]*subvol); 
  atomicAdd(&out[37], 0.1767766952966368*win[0]*fin[65]*subvol); 
  atomicAdd(&out[38], 0.1767766952966368*win[0]*fin[69]*subvol); 
  atomicAdd(&out[39], 0.1767766952966368*win[0]*fin[70]*subvol); 
  atomicAdd(&out[40], 0.1767766952966368*win[0]*fin[71]*subvol); 
  atomicAdd(&out[41], 0.1767766952966368*win[0]*fin[81]*subvol); 
  atomicAdd(&out[42], 0.1767766952966368*win[0]*fin[82]*subvol); 
  atomicAdd(&out[43], 0.1767766952966368*win[0]*fin[83]*subvol); 
  atomicAdd(&out[44], 0.1767766952966368*win[0]*fin[91]*subvol); 
  atomicAdd(&out[45], 0.1767766952966368*win[0]*fin[92]*subvol); 
  atomicAdd(&out[46], 0.1767766952966368*win[0]*fin[93]*subvol); 
  atomicAdd(&out[47], 0.1767766952966368*win[0]*fin[95]*subvol); 
  atomicAdd(&out[48], 0.1767766952966368*win[0]*fin[96]*subvol); 
  atomicAdd(&out[49], 0.1767766952966368*win[0]*fin[97]*subvol); 
  atomicAdd(&out[50], 0.1767766952966368*win[0]*fin[101]*subvol); 
  atomicAdd(&out[51], 0.1767766952966368*win[0]*fin[102]*subvol); 
  atomicAdd(&out[52], 0.1767766952966368*win[0]*fin[103]*subvol); 
  atomicAdd(&out[53], 0.1767766952966368*win[0]*fin[107]*subvol); 
  atomicAdd(&out[54], 0.1767766952966368*win[0]*fin[110]*subvol); 
  atomicAdd(&out[55], 0.1767766952966368*win[0]*fin[118]*subvol); 
  atomicAdd(&out[56], 0.1767766952966368*win[0]*fin[128]*subvol); 
  atomicAdd(&out[57], 0.1767766952966368*win[0]*fin[129]*subvol); 
  atomicAdd(&out[58], 0.1767766952966368*win[0]*fin[130]*subvol); 
  atomicAdd(&out[59], 0.1767766952966368*win[0]*fin[134]*subvol); 
  atomicAdd(&out[60], 0.1767766952966368*win[0]*fin[138]*subvol); 
  atomicAdd(&out[61], 0.1767766952966368*win[0]*fin[139]*subvol); 
  atomicAdd(&out[62], 0.1767766952966368*win[0]*fin[140]*subvol); 
  atomicAdd(&out[63], 0.1767766952966368*win[0]*fin[152]*subvol); 
  atomicAdd(&out[64], 0.1767766952966368*win[0]*fin[153]*subvol); 
  atomicAdd(&out[65], 0.1767766952966368*win[0]*fin[154]*subvol); 
  atomicAdd(&out[66], 0.1767766952966368*win[0]*fin[163]*subvol); 
  atomicAdd(&out[67], 0.1767766952966368*win[0]*fin[168]*subvol); 
  atomicAdd(&out[68], 0.1767766952966368*win[0]*fin[184]*subvol); 
  atomicAdd(&out[69], 0.1767766952966368*win[0]*fin[185]*subvol); 
  atomicAdd(&out[70], 0.1767766952966368*win[0]*fin[186]*subvol); 
  atomicAdd(&out[71], 0.1767766952966368*win[0]*fin[202]*subvol); 
#else
  out[0] += 0.1767766952966368*fin[0]*win[0]*subvol; 
  out[1] += 0.1767766952966368*win[0]*fin[1]*subvol; 
  out[2] += 0.1767766952966368*win[0]*fin[2]*subvol; 
  out[3] += 0.1767766952966368*win[0]*fin[3]*subvol; 
  out[4] += 0.1767766952966368*win[0]*fin[5]*subvol; 
  out[5] += 0.1767766952966368*win[0]*fin[6]*subvol; 
  out[6] += 0.1767766952966368*win[0]*fin[7]*subvol; 
  out[7] += 0.1767766952966368*win[0]*fin[8]*subvol; 
  out[8] += 0.1767766952966368*win[0]*fin[9]*subvol; 
  out[9] += 0.1767766952966368*win[0]*fin[13]*subvol; 
  out[10] += 0.1767766952966368*win[0]*fin[14]*subvol; 
  out[11] += 0.1767766952966368*win[0]*fin[15]*subvol; 
  out[12] += 0.1767766952966368*win[0]*fin[17]*subvol; 
  out[13] += 0.1767766952966368*win[0]*fin[18]*subvol; 
  out[14] += 0.1767766952966368*win[0]*fin[19]*subvol; 
  out[15] += 0.1767766952966368*win[0]*fin[21]*subvol; 
  out[16] += 0.1767766952966368*win[0]*fin[23]*subvol; 
  out[17] += 0.1767766952966368*win[0]*fin[24]*subvol; 
  out[18] += 0.1767766952966368*win[0]*fin[25]*subvol; 
  out[19] += 0.1767766952966368*win[0]*fin[29]*subvol; 
  out[20] += 0.1767766952966368*win[0]*fin[30]*subvol; 
  out[21] += 0.1767766952966368*win[0]*fin[31]*subvol; 
  out[22] += 0.1767766952966368*win[0]*fin[35]*subvol; 
  out[23] += 0.1767766952966368*win[0]*fin[36]*subvol; 
  out[24] += 0.1767766952966368*win[0]*fin[37]*subvol; 
  out[25] += 0.1767766952966368*win[0]*fin[41]*subvol; 
  out[26] += 0.1767766952966368*win[0]*fin[42]*subvol; 
  out[27] += 0.1767766952966368*win[0]*fin[43]*subvol; 
  out[28] += 0.1767766952966368*win[0]*fin[49]*subvol; 
  out[29] += 0.1767766952966368*win[0]*fin[50]*subvol; 
  out[30] += 0.1767766952966368*win[0]*fin[51]*subvol; 
  out[31] += 0.1767766952966368*win[0]*fin[54]*subvol; 
  out[32] += 0.1767766952966368*win[0]*fin[55]*subvol; 
  out[33] += 0.1767766952966368*win[0]*fin[56]*subvol; 
  out[34] += 0.1767766952966368*win[0]*fin[57]*subvol; 
  out[35] += 0.1767766952966368*win[0]*fin[59]*subvol; 
  out[36] += 0.1767766952966368*win[0]*fin[61]*subvol; 
  out[37] += 0.1767766952966368*win[0]*fin[65]*subvol; 
  out[38] += 0.1767766952966368*win[0]*fin[69]*subvol; 
  out[39] += 0.1767766952966368*win[0]*fin[70]*subvol; 
  out[40] += 0.1767766952966368*win[0]*fin[71]*subvol; 
  out[41] += 0.1767766952966368*win[0]*fin[81]*subvol; 
  out[42] += 0.1767766952966368*win[0]*fin[82]*subvol; 
  out[43] += 0.1767766952966368*win[0]*fin[83]*subvol; 
  out[44] += 0.1767766952966368*win[0]*fin[91]*subvol; 
  out[45] += 0.1767766952966368*win[0]*fin[92]*subvol; 
  out[46] += 0.1767766952966368*win[0]*fin[93]*subvol; 
  out[47] += 0.1767766952966368*win[0]*fin[95]*subvol; 
  out[48] += 0.1767766952966368*win[0]*fin[96]*subvol; 
  out[49] += 0.1767766952966368*win[0]*fin[97]*subvol; 
  out[50] += 0.1767766952966368*win[0]*fin[101]*subvol; 
  out[51] += 0.1767766952966368*win[0]*fin[102]*subvol; 
  out[52] += 0.1767766952966368*win[0]*fin[103]*subvol; 
  out[53] += 0.1767766952966368*win[0]*fin[107]*subvol; 
  out[54] += 0.1767766952966368*win[0]*fin[110]*subvol; 
  out[55] += 0.1767766952966368*win[0]*fin[118]*subvol; 
  out[56] += 0.1767766952966368*win[0]*fin[128]*subvol; 
  out[57] += 0.1767766952966368*win[0]*fin[129]*subvol; 
  out[58] += 0.1767766952966368*win[0]*fin[130]*subvol; 
  out[59] += 0.1767766952966368*win[0]*fin[134]*subvol; 
  out[60] += 0.1767766952966368*win[0]*fin[138]*subvol; 
  out[61] += 0.1767766952966368*win[0]*fin[139]*subvol; 
  out[62] += 0.1767766952966368*win[0]*fin[140]*subvol; 
  out[63] += 0.1767766952966368*win[0]*fin[152]*subvol; 
  out[64] += 0.1767766952966368*win[0]*fin[153]*subvol; 
  out[65] += 0.1767766952966368*win[0]*fin[154]*subvol; 
  out[66] += 0.1767766952966368*win[0]*fin[163]*subvol; 
  out[67] += 0.1767766952966368*win[0]*fin[168]*subvol; 
  out[68] += 0.1767766952966368*win[0]*fin[184]*subvol; 
  out[69] += 0.1767766952966368*win[0]*fin[185]*subvol; 
  out[70] += 0.1767766952966368*win[0]*fin[186]*subvol; 
  out[71] += 0.1767766952966368*win[0]*fin[202]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgvx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.25*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.25*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.25*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.25*win[0]*fin[3]*subvol); 
  atomicAdd(&out[5], 0.25*win[0]*fin[5]*subvol); 
  atomicAdd(&out[6], 0.25*win[0]*fin[6]*subvol); 
  atomicAdd(&out[7], 0.25*win[0]*fin[7]*subvol); 
  atomicAdd(&out[8], 0.25*win[0]*fin[8]*subvol); 
  atomicAdd(&out[9], 0.25*win[0]*fin[9]*subvol); 
  atomicAdd(&out[13], 0.25*win[0]*fin[13]*subvol); 
  atomicAdd(&out[14], 0.25*win[0]*fin[14]*subvol); 
  atomicAdd(&out[15], 0.25*win[0]*fin[15]*subvol); 
  atomicAdd(&out[17], 0.25*win[0]*fin[17]*subvol); 
  atomicAdd(&out[18], 0.25*win[0]*fin[18]*subvol); 
  atomicAdd(&out[19], 0.25*win[0]*fin[19]*subvol); 
  atomicAdd(&out[21], 0.25*win[0]*fin[21]*subvol); 
  atomicAdd(&out[23], 0.25*win[0]*fin[23]*subvol); 
  atomicAdd(&out[24], 0.25*win[0]*fin[24]*subvol); 
  atomicAdd(&out[25], 0.25*win[0]*fin[25]*subvol); 
  atomicAdd(&out[29], 0.25*win[0]*fin[29]*subvol); 
  atomicAdd(&out[30], 0.25*win[0]*fin[30]*subvol); 
  atomicAdd(&out[31], 0.25*win[0]*fin[31]*subvol); 
  atomicAdd(&out[35], 0.25*win[0]*fin[35]*subvol); 
  atomicAdd(&out[36], 0.25*win[0]*fin[36]*subvol); 
  atomicAdd(&out[37], 0.25*win[0]*fin[37]*subvol); 
  atomicAdd(&out[41], 0.25*win[0]*fin[41]*subvol); 
  atomicAdd(&out[42], 0.25*win[0]*fin[42]*subvol); 
  atomicAdd(&out[43], 0.25*win[0]*fin[43]*subvol); 
  atomicAdd(&out[49], 0.25*win[0]*fin[49]*subvol); 
  atomicAdd(&out[50], 0.25*win[0]*fin[50]*subvol); 
  atomicAdd(&out[51], 0.25*win[0]*fin[51]*subvol); 
  atomicAdd(&out[54], 0.25*win[0]*fin[54]*subvol); 
  atomicAdd(&out[55], 0.25*win[0]*fin[55]*subvol); 
  atomicAdd(&out[56], 0.25*win[0]*fin[56]*subvol); 
  atomicAdd(&out[57], 0.25*win[0]*fin[57]*subvol); 
  atomicAdd(&out[59], 0.25*win[0]*fin[59]*subvol); 
  atomicAdd(&out[61], 0.25*win[0]*fin[61]*subvol); 
  atomicAdd(&out[65], 0.25*win[0]*fin[65]*subvol); 
  atomicAdd(&out[69], 0.25*win[0]*fin[69]*subvol); 
  atomicAdd(&out[70], 0.25*win[0]*fin[70]*subvol); 
  atomicAdd(&out[71], 0.25*win[0]*fin[71]*subvol); 
  atomicAdd(&out[81], 0.25*win[0]*fin[81]*subvol); 
  atomicAdd(&out[82], 0.25*win[0]*fin[82]*subvol); 
  atomicAdd(&out[83], 0.25*win[0]*fin[83]*subvol); 
  atomicAdd(&out[91], 0.25*win[0]*fin[91]*subvol); 
  atomicAdd(&out[92], 0.25*win[0]*fin[92]*subvol); 
  atomicAdd(&out[93], 0.25*win[0]*fin[93]*subvol); 
  atomicAdd(&out[95], 0.25*win[0]*fin[95]*subvol); 
  atomicAdd(&out[96], 0.25*win[0]*fin[96]*subvol); 
  atomicAdd(&out[97], 0.25*win[0]*fin[97]*subvol); 
  atomicAdd(&out[101], 0.25*win[0]*fin[101]*subvol); 
  atomicAdd(&out[102], 0.25*win[0]*fin[102]*subvol); 
  atomicAdd(&out[103], 0.25*win[0]*fin[103]*subvol); 
  atomicAdd(&out[107], 0.25*win[0]*fin[107]*subvol); 
  atomicAdd(&out[110], 0.25*win[0]*fin[110]*subvol); 
  atomicAdd(&out[118], 0.25*win[0]*fin[118]*subvol); 
  atomicAdd(&out[128], 0.25*win[0]*fin[128]*subvol); 
  atomicAdd(&out[129], 0.25*win[0]*fin[129]*subvol); 
  atomicAdd(&out[130], 0.25*win[0]*fin[130]*subvol); 
  atomicAdd(&out[134], 0.25*win[0]*fin[134]*subvol); 
  atomicAdd(&out[138], 0.25*win[0]*fin[138]*subvol); 
  atomicAdd(&out[139], 0.25*win[0]*fin[139]*subvol); 
  atomicAdd(&out[140], 0.25*win[0]*fin[140]*subvol); 
  atomicAdd(&out[152], 0.25*win[0]*fin[152]*subvol); 
  atomicAdd(&out[153], 0.25*win[0]*fin[153]*subvol); 
  atomicAdd(&out[154], 0.25*win[0]*fin[154]*subvol); 
  atomicAdd(&out[163], 0.25*win[0]*fin[163]*subvol); 
  atomicAdd(&out[168], 0.25*win[0]*fin[168]*subvol); 
  atomicAdd(&out[184], 0.25*win[0]*fin[184]*subvol); 
  atomicAdd(&out[185], 0.25*win[0]*fin[185]*subvol); 
  atomicAdd(&out[186], 0.25*win[0]*fin[186]*subvol); 
  atomicAdd(&out[202], 0.25*win[0]*fin[202]*subvol); 
#else
  out[0] += 0.25*fin[0]*win[0]*subvol; 
  out[1] += 0.25*win[0]*fin[1]*subvol; 
  out[2] += 0.25*win[0]*fin[2]*subvol; 
  out[3] += 0.25*win[0]*fin[3]*subvol; 
  out[5] += 0.25*win[0]*fin[5]*subvol; 
  out[6] += 0.25*win[0]*fin[6]*subvol; 
  out[7] += 0.25*win[0]*fin[7]*subvol; 
  out[8] += 0.25*win[0]*fin[8]*subvol; 
  out[9] += 0.25*win[0]*fin[9]*subvol; 
  out[13] += 0.25*win[0]*fin[13]*subvol; 
  out[14] += 0.25*win[0]*fin[14]*subvol; 
  out[15] += 0.25*win[0]*fin[15]*subvol; 
  out[17] += 0.25*win[0]*fin[17]*subvol; 
  out[18] += 0.25*win[0]*fin[18]*subvol; 
  out[19] += 0.25*win[0]*fin[19]*subvol; 
  out[21] += 0.25*win[0]*fin[21]*subvol; 
  out[23] += 0.25*win[0]*fin[23]*subvol; 
  out[24] += 0.25*win[0]*fin[24]*subvol; 
  out[25] += 0.25*win[0]*fin[25]*subvol; 
  out[29] += 0.25*win[0]*fin[29]*subvol; 
  out[30] += 0.25*win[0]*fin[30]*subvol; 
  out[31] += 0.25*win[0]*fin[31]*subvol; 
  out[35] += 0.25*win[0]*fin[35]*subvol; 
  out[36] += 0.25*win[0]*fin[36]*subvol; 
  out[37] += 0.25*win[0]*fin[37]*subvol; 
  out[41] += 0.25*win[0]*fin[41]*subvol; 
  out[42] += 0.25*win[0]*fin[42]*subvol; 
  out[43] += 0.25*win[0]*fin[43]*subvol; 
  out[49] += 0.25*win[0]*fin[49]*subvol; 
  out[50] += 0.25*win[0]*fin[50]*subvol; 
  out[51] += 0.25*win[0]*fin[51]*subvol; 
  out[54] += 0.25*win[0]*fin[54]*subvol; 
  out[55] += 0.25*win[0]*fin[55]*subvol; 
  out[56] += 0.25*win[0]*fin[56]*subvol; 
  out[57] += 0.25*win[0]*fin[57]*subvol; 
  out[59] += 0.25*win[0]*fin[59]*subvol; 
  out[61] += 0.25*win[0]*fin[61]*subvol; 
  out[65] += 0.25*win[0]*fin[65]*subvol; 
  out[69] += 0.25*win[0]*fin[69]*subvol; 
  out[70] += 0.25*win[0]*fin[70]*subvol; 
  out[71] += 0.25*win[0]*fin[71]*subvol; 
  out[81] += 0.25*win[0]*fin[81]*subvol; 
  out[82] += 0.25*win[0]*fin[82]*subvol; 
  out[83] += 0.25*win[0]*fin[83]*subvol; 
  out[91] += 0.25*win[0]*fin[91]*subvol; 
  out[92] += 0.25*win[0]*fin[92]*subvol; 
  out[93] += 0.25*win[0]*fin[93]*subvol; 
  out[95] += 0.25*win[0]*fin[95]*subvol; 
  out[96] += 0.25*win[0]*fin[96]*subvol; 
  out[97] += 0.25*win[0]*fin[97]*subvol; 
  out[101] += 0.25*win[0]*fin[101]*subvol; 
  out[102] += 0.25*win[0]*fin[102]*subvol; 
  out[103] += 0.25*win[0]*fin[103]*subvol; 
  out[107] += 0.25*win[0]*fin[107]*subvol; 
  out[110] += 0.25*win[0]*fin[110]*subvol; 
  out[118] += 0.25*win[0]*fin[118]*subvol; 
  out[128] += 0.25*win[0]*fin[128]*subvol; 
  out[129] += 0.25*win[0]*fin[129]*subvol; 
  out[130] += 0.25*win[0]*fin[130]*subvol; 
  out[134] += 0.25*win[0]*fin[134]*subvol; 
  out[138] += 0.25*win[0]*fin[138]*subvol; 
  out[139] += 0.25*win[0]*fin[139]*subvol; 
  out[140] += 0.25*win[0]*fin[140]*subvol; 
  out[152] += 0.25*win[0]*fin[152]*subvol; 
  out[153] += 0.25*win[0]*fin[153]*subvol; 
  out[154] += 0.25*win[0]*fin[154]*subvol; 
  out[163] += 0.25*win[0]*fin[163]*subvol; 
  out[168] += 0.25*win[0]*fin[168]*subvol; 
  out[184] += 0.25*win[0]*fin[184]*subvol; 
  out[185] += 0.25*win[0]*fin[185]*subvol; 
  out[186] += 0.25*win[0]*fin[186]*subvol; 
  out[202] += 0.25*win[0]*fin[202]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxvx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.25*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.25*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.25*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.25*win[0]*fin[5]*subvol); 
  atomicAdd(&out[4], 0.25*win[0]*fin[6]*subvol); 
  atomicAdd(&out[5], 0.25*win[0]*fin[9]*subvol); 
  atomicAdd(&out[6], 0.25*win[0]*fin[14]*subvol); 
  atomicAdd(&out[7], 0.25*win[0]*fin[15]*subvol); 
  atomicAdd(&out[8], 0.25*win[0]*fin[18]*subvol); 
  atomicAdd(&out[9], 0.25*win[0]*fin[19]*subvol); 
  atomicAdd(&out[10], 0.25*win[0]*fin[21]*subvol); 
  atomicAdd(&out[11], 0.25*win[0]*fin[23]*subvol); 
  atomicAdd(&out[12], 0.25*win[0]*fin[24]*subvol); 
  atomicAdd(&out[13], 0.25*win[0]*fin[31]*subvol); 
  atomicAdd(&out[14], 0.25*win[0]*fin[37]*subvol); 
  atomicAdd(&out[15], 0.25*win[0]*fin[42]*subvol); 
  atomicAdd(&out[16], 0.25*win[0]*fin[43]*subvol); 
  atomicAdd(&out[17], 0.25*win[0]*fin[50]*subvol); 
  atomicAdd(&out[18], 0.25*win[0]*fin[51]*subvol); 
  atomicAdd(&out[19], 0.25*win[0]*fin[54]*subvol); 
  atomicAdd(&out[20], 0.25*win[0]*fin[56]*subvol); 
  atomicAdd(&out[21], 0.25*win[0]*fin[57]*subvol); 
  atomicAdd(&out[22], 0.25*win[0]*fin[59]*subvol); 
  atomicAdd(&out[23], 0.25*win[0]*fin[71]*subvol); 
  atomicAdd(&out[24], 0.25*win[0]*fin[83]*subvol); 
  atomicAdd(&out[25], 0.25*win[0]*fin[92]*subvol); 
  atomicAdd(&out[26], 0.25*win[0]*fin[93]*subvol); 
  atomicAdd(&out[27], 0.25*win[0]*fin[97]*subvol); 
  atomicAdd(&out[28], 0.25*win[0]*fin[102]*subvol); 
  atomicAdd(&out[29], 0.25*win[0]*fin[103]*subvol); 
  atomicAdd(&out[30], 0.25*win[0]*fin[107]*subvol); 
  atomicAdd(&out[31], 0.25*win[0]*fin[130]*subvol); 
  atomicAdd(&out[32], 0.25*win[0]*fin[140]*subvol); 
  atomicAdd(&out[33], 0.25*win[0]*fin[153]*subvol); 
  atomicAdd(&out[34], 0.25*win[0]*fin[154]*subvol); 
  atomicAdd(&out[35], 0.25*win[0]*fin[186]*subvol); 
#else
  out[0] += 0.25*fin[0]*win[0]*subvol; 
  out[1] += 0.25*win[0]*fin[2]*subvol; 
  out[2] += 0.25*win[0]*fin[3]*subvol; 
  out[3] += 0.25*win[0]*fin[5]*subvol; 
  out[4] += 0.25*win[0]*fin[6]*subvol; 
  out[5] += 0.25*win[0]*fin[9]*subvol; 
  out[6] += 0.25*win[0]*fin[14]*subvol; 
  out[7] += 0.25*win[0]*fin[15]*subvol; 
  out[8] += 0.25*win[0]*fin[18]*subvol; 
  out[9] += 0.25*win[0]*fin[19]*subvol; 
  out[10] += 0.25*win[0]*fin[21]*subvol; 
  out[11] += 0.25*win[0]*fin[23]*subvol; 
  out[12] += 0.25*win[0]*fin[24]*subvol; 
  out[13] += 0.25*win[0]*fin[31]*subvol; 
  out[14] += 0.25*win[0]*fin[37]*subvol; 
  out[15] += 0.25*win[0]*fin[42]*subvol; 
  out[16] += 0.25*win[0]*fin[43]*subvol; 
  out[17] += 0.25*win[0]*fin[50]*subvol; 
  out[18] += 0.25*win[0]*fin[51]*subvol; 
  out[19] += 0.25*win[0]*fin[54]*subvol; 
  out[20] += 0.25*win[0]*fin[56]*subvol; 
  out[21] += 0.25*win[0]*fin[57]*subvol; 
  out[22] += 0.25*win[0]*fin[59]*subvol; 
  out[23] += 0.25*win[0]*fin[71]*subvol; 
  out[24] += 0.25*win[0]*fin[83]*subvol; 
  out[25] += 0.25*win[0]*fin[92]*subvol; 
  out[26] += 0.25*win[0]*fin[93]*subvol; 
  out[27] += 0.25*win[0]*fin[97]*subvol; 
  out[28] += 0.25*win[0]*fin[102]*subvol; 
  out[29] += 0.25*win[0]*fin[103]*subvol; 
  out[30] += 0.25*win[0]*fin[107]*subvol; 
  out[31] += 0.25*win[0]*fin[130]*subvol; 
  out[32] += 0.25*win[0]*fin[140]*subvol; 
  out[33] += 0.25*win[0]*fin[153]*subvol; 
  out[34] += 0.25*win[0]*fin[154]*subvol; 
  out[35] += 0.25*win[0]*fin[186]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxvx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[5]*subvol); 
  atomicAdd(&out[6], 0.5*win[0]*fin[6]*subvol); 
  atomicAdd(&out[9], 0.5*win[0]*fin[9]*subvol); 
  atomicAdd(&out[14], 0.5*win[0]*fin[14]*subvol); 
  atomicAdd(&out[15], 0.5*win[0]*fin[15]*subvol); 
  atomicAdd(&out[18], 0.5*win[0]*fin[18]*subvol); 
  atomicAdd(&out[19], 0.5*win[0]*fin[19]*subvol); 
  atomicAdd(&out[21], 0.5*win[0]*fin[21]*subvol); 
  atomicAdd(&out[23], 0.5*win[0]*fin[23]*subvol); 
  atomicAdd(&out[24], 0.5*win[0]*fin[24]*subvol); 
  atomicAdd(&out[31], 0.5*win[0]*fin[31]*subvol); 
  atomicAdd(&out[37], 0.5*win[0]*fin[37]*subvol); 
  atomicAdd(&out[42], 0.5*win[0]*fin[42]*subvol); 
  atomicAdd(&out[43], 0.5*win[0]*fin[43]*subvol); 
  atomicAdd(&out[50], 0.5*win[0]*fin[50]*subvol); 
  atomicAdd(&out[51], 0.5*win[0]*fin[51]*subvol); 
  atomicAdd(&out[54], 0.5*win[0]*fin[54]*subvol); 
  atomicAdd(&out[56], 0.5*win[0]*fin[56]*subvol); 
  atomicAdd(&out[57], 0.5*win[0]*fin[57]*subvol); 
  atomicAdd(&out[59], 0.5*win[0]*fin[59]*subvol); 
  atomicAdd(&out[71], 0.5*win[0]*fin[71]*subvol); 
  atomicAdd(&out[83], 0.5*win[0]*fin[83]*subvol); 
  atomicAdd(&out[92], 0.5*win[0]*fin[92]*subvol); 
  atomicAdd(&out[93], 0.5*win[0]*fin[93]*subvol); 
  atomicAdd(&out[97], 0.5*win[0]*fin[97]*subvol); 
  atomicAdd(&out[102], 0.5*win[0]*fin[102]*subvol); 
  atomicAdd(&out[103], 0.5*win[0]*fin[103]*subvol); 
  atomicAdd(&out[107], 0.5*win[0]*fin[107]*subvol); 
  atomicAdd(&out[130], 0.5*win[0]*fin[130]*subvol); 
  atomicAdd(&out[140], 0.5*win[0]*fin[140]*subvol); 
  atomicAdd(&out[153], 0.5*win[0]*fin[153]*subvol); 
  atomicAdd(&out[154], 0.5*win[0]*fin[154]*subvol); 
  atomicAdd(&out[186], 0.5*win[0]*fin[186]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[2] += 0.5*win[0]*fin[2]*subvol; 
  out[3] += 0.5*win[0]*fin[3]*subvol; 
  out[5] += 0.5*win[0]*fin[5]*subvol; 
  out[6] += 0.5*win[0]*fin[6]*subvol; 
  out[9] += 0.5*win[0]*fin[9]*subvol; 
  out[14] += 0.5*win[0]*fin[14]*subvol; 
  out[15] += 0.5*win[0]*fin[15]*subvol; 
  out[18] += 0.5*win[0]*fin[18]*subvol; 
  out[19] += 0.5*win[0]*fin[19]*subvol; 
  out[21] += 0.5*win[0]*fin[21]*subvol; 
  out[23] += 0.5*win[0]*fin[23]*subvol; 
  out[24] += 0.5*win[0]*fin[24]*subvol; 
  out[31] += 0.5*win[0]*fin[31]*subvol; 
  out[37] += 0.5*win[0]*fin[37]*subvol; 
  out[42] += 0.5*win[0]*fin[42]*subvol; 
  out[43] += 0.5*win[0]*fin[43]*subvol; 
  out[50] += 0.5*win[0]*fin[50]*subvol; 
  out[51] += 0.5*win[0]*fin[51]*subvol; 
  out[54] += 0.5*win[0]*fin[54]*subvol; 
  out[56] += 0.5*win[0]*fin[56]*subvol; 
  out[57] += 0.5*win[0]*fin[57]*subvol; 
  out[59] += 0.5*win[0]*fin[59]*subvol; 
  out[71] += 0.5*win[0]*fin[71]*subvol; 
  out[83] += 0.5*win[0]*fin[83]*subvol; 
  out[92] += 0.5*win[0]*fin[92]*subvol; 
  out[93] += 0.5*win[0]*fin[93]*subvol; 
  out[97] += 0.5*win[0]*fin[97]*subvol; 
  out[102] += 0.5*win[0]*fin[102]*subvol; 
  out[103] += 0.5*win[0]*fin[103]*subvol; 
  out[107] += 0.5*win[0]*fin[107]*subvol; 
  out[130] += 0.5*win[0]*fin[130]*subvol; 
  out[140] += 0.5*win[0]*fin[140]*subvol; 
  out[153] += 0.5*win[0]*fin[153]*subvol; 
  out[154] += 0.5*win[0]*fin[154]*subvol; 
  out[186] += 0.5*win[0]*fin[186]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgyvx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.25*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.25*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.25*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.25*win[0]*fin[5]*subvol); 
  atomicAdd(&out[4], 0.25*win[0]*fin[6]*subvol); 
  atomicAdd(&out[5], 0.25*win[0]*fin[8]*subvol); 
  atomicAdd(&out[6], 0.25*win[0]*fin[13]*subvol); 
  atomicAdd(&out[7], 0.25*win[0]*fin[15]*subvol); 
  atomicAdd(&out[8], 0.25*win[0]*fin[17]*subvol); 
  atomicAdd(&out[9], 0.25*win[0]*fin[19]*subvol); 
  atomicAdd(&out[10], 0.25*win[0]*fin[21]*subvol); 
  atomicAdd(&out[11], 0.25*win[0]*fin[23]*subvol); 
  atomicAdd(&out[12], 0.25*win[0]*fin[24]*subvol); 
  atomicAdd(&out[13], 0.25*win[0]*fin[30]*subvol); 
  atomicAdd(&out[14], 0.25*win[0]*fin[36]*subvol); 
  atomicAdd(&out[15], 0.25*win[0]*fin[41]*subvol); 
  atomicAdd(&out[16], 0.25*win[0]*fin[43]*subvol); 
  atomicAdd(&out[17], 0.25*win[0]*fin[49]*subvol); 
  atomicAdd(&out[18], 0.25*win[0]*fin[51]*subvol); 
  atomicAdd(&out[19], 0.25*win[0]*fin[54]*subvol); 
  atomicAdd(&out[20], 0.25*win[0]*fin[55]*subvol); 
  atomicAdd(&out[21], 0.25*win[0]*fin[57]*subvol); 
  atomicAdd(&out[22], 0.25*win[0]*fin[59]*subvol); 
  atomicAdd(&out[23], 0.25*win[0]*fin[70]*subvol); 
  atomicAdd(&out[24], 0.25*win[0]*fin[82]*subvol); 
  atomicAdd(&out[25], 0.25*win[0]*fin[91]*subvol); 
  atomicAdd(&out[26], 0.25*win[0]*fin[93]*subvol); 
  atomicAdd(&out[27], 0.25*win[0]*fin[96]*subvol); 
  atomicAdd(&out[28], 0.25*win[0]*fin[101]*subvol); 
  atomicAdd(&out[29], 0.25*win[0]*fin[103]*subvol); 
  atomicAdd(&out[30], 0.25*win[0]*fin[107]*subvol); 
  atomicAdd(&out[31], 0.25*win[0]*fin[129]*subvol); 
  atomicAdd(&out[32], 0.25*win[0]*fin[139]*subvol); 
  atomicAdd(&out[33], 0.25*win[0]*fin[152]*subvol); 
  atomicAdd(&out[34], 0.25*win[0]*fin[154]*subvol); 
  atomicAdd(&out[35], 0.25*win[0]*fin[185]*subvol); 
#else
  out[0] += 0.25*fin[0]*win[0]*subvol; 
  out[1] += 0.25*win[0]*fin[1]*subvol; 
  out[2] += 0.25*win[0]*fin[3]*subvol; 
  out[3] += 0.25*win[0]*fin[5]*subvol; 
  out[4] += 0.25*win[0]*fin[6]*subvol; 
  out[5] += 0.25*win[0]*fin[8]*subvol; 
  out[6] += 0.25*win[0]*fin[13]*subvol; 
  out[7] += 0.25*win[0]*fin[15]*subvol; 
  out[8] += 0.25*win[0]*fin[17]*subvol; 
  out[9] += 0.25*win[0]*fin[19]*subvol; 
  out[10] += 0.25*win[0]*fin[21]*subvol; 
  out[11] += 0.25*win[0]*fin[23]*subvol; 
  out[12] += 0.25*win[0]*fin[24]*subvol; 
  out[13] += 0.25*win[0]*fin[30]*subvol; 
  out[14] += 0.25*win[0]*fin[36]*subvol; 
  out[15] += 0.25*win[0]*fin[41]*subvol; 
  out[16] += 0.25*win[0]*fin[43]*subvol; 
  out[17] += 0.25*win[0]*fin[49]*subvol; 
  out[18] += 0.25*win[0]*fin[51]*subvol; 
  out[19] += 0.25*win[0]*fin[54]*subvol; 
  out[20] += 0.25*win[0]*fin[55]*subvol; 
  out[21] += 0.25*win[0]*fin[57]*subvol; 
  out[22] += 0.25*win[0]*fin[59]*subvol; 
  out[23] += 0.25*win[0]*fin[70]*subvol; 
  out[24] += 0.25*win[0]*fin[82]*subvol; 
  out[25] += 0.25*win[0]*fin[91]*subvol; 
  out[26] += 0.25*win[0]*fin[93]*subvol; 
  out[27] += 0.25*win[0]*fin[96]*subvol; 
  out[28] += 0.25*win[0]*fin[101]*subvol; 
  out[29] += 0.25*win[0]*fin[103]*subvol; 
  out[30] += 0.25*win[0]*fin[107]*subvol; 
  out[31] += 0.25*win[0]*fin[129]*subvol; 
  out[32] += 0.25*win[0]*fin[139]*subvol; 
  out[33] += 0.25*win[0]*fin[152]*subvol; 
  out[34] += 0.25*win[0]*fin[154]*subvol; 
  out[35] += 0.25*win[0]*fin[185]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgyvx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[5]*subvol); 
  atomicAdd(&out[6], 0.5*win[0]*fin[6]*subvol); 
  atomicAdd(&out[8], 0.5*win[0]*fin[8]*subvol); 
  atomicAdd(&out[13], 0.5*win[0]*fin[13]*subvol); 
  atomicAdd(&out[15], 0.5*win[0]*fin[15]*subvol); 
  atomicAdd(&out[17], 0.5*win[0]*fin[17]*subvol); 
  atomicAdd(&out[19], 0.5*win[0]*fin[19]*subvol); 
  atomicAdd(&out[21], 0.5*win[0]*fin[21]*subvol); 
  atomicAdd(&out[23], 0.5*win[0]*fin[23]*subvol); 
  atomicAdd(&out[24], 0.5*win[0]*fin[24]*subvol); 
  atomicAdd(&out[30], 0.5*win[0]*fin[30]*subvol); 
  atomicAdd(&out[36], 0.5*win[0]*fin[36]*subvol); 
  atomicAdd(&out[41], 0.5*win[0]*fin[41]*subvol); 
  atomicAdd(&out[43], 0.5*win[0]*fin[43]*subvol); 
  atomicAdd(&out[49], 0.5*win[0]*fin[49]*subvol); 
  atomicAdd(&out[51], 0.5*win[0]*fin[51]*subvol); 
  atomicAdd(&out[54], 0.5*win[0]*fin[54]*subvol); 
  atomicAdd(&out[55], 0.5*win[0]*fin[55]*subvol); 
  atomicAdd(&out[57], 0.5*win[0]*fin[57]*subvol); 
  atomicAdd(&out[59], 0.5*win[0]*fin[59]*subvol); 
  atomicAdd(&out[70], 0.5*win[0]*fin[70]*subvol); 
  atomicAdd(&out[82], 0.5*win[0]*fin[82]*subvol); 
  atomicAdd(&out[91], 0.5*win[0]*fin[91]*subvol); 
  atomicAdd(&out[93], 0.5*win[0]*fin[93]*subvol); 
  atomicAdd(&out[96], 0.5*win[0]*fin[96]*subvol); 
  atomicAdd(&out[101], 0.5*win[0]*fin[101]*subvol); 
  atomicAdd(&out[103], 0.5*win[0]*fin[103]*subvol); 
  atomicAdd(&out[107], 0.5*win[0]*fin[107]*subvol); 
  atomicAdd(&out[129], 0.5*win[0]*fin[129]*subvol); 
  atomicAdd(&out[139], 0.5*win[0]*fin[139]*subvol); 
  atomicAdd(&out[152], 0.5*win[0]*fin[152]*subvol); 
  atomicAdd(&out[154], 0.5*win[0]*fin[154]*subvol); 
  atomicAdd(&out[185], 0.5*win[0]*fin[185]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[3] += 0.5*win[0]*fin[3]*subvol; 
  out[5] += 0.5*win[0]*fin[5]*subvol; 
  out[6] += 0.5*win[0]*fin[6]*subvol; 
  out[8] += 0.5*win[0]*fin[8]*subvol; 
  out[13] += 0.5*win[0]*fin[13]*subvol; 
  out[15] += 0.5*win[0]*fin[15]*subvol; 
  out[17] += 0.5*win[0]*fin[17]*subvol; 
  out[19] += 0.5*win[0]*fin[19]*subvol; 
  out[21] += 0.5*win[0]*fin[21]*subvol; 
  out[23] += 0.5*win[0]*fin[23]*subvol; 
  out[24] += 0.5*win[0]*fin[24]*subvol; 
  out[30] += 0.5*win[0]*fin[30]*subvol; 
  out[36] += 0.5*win[0]*fin[36]*subvol; 
  out[41] += 0.5*win[0]*fin[41]*subvol; 
  out[43] += 0.5*win[0]*fin[43]*subvol; 
  out[49] += 0.5*win[0]*fin[49]*subvol; 
  out[51] += 0.5*win[0]*fin[51]*subvol; 
  out[54] += 0.5*win[0]*fin[54]*subvol; 
  out[55] += 0.5*win[0]*fin[55]*subvol; 
  out[57] += 0.5*win[0]*fin[57]*subvol; 
  out[59] += 0.5*win[0]*fin[59]*subvol; 
  out[70] += 0.5*win[0]*fin[70]*subvol; 
  out[82] += 0.5*win[0]*fin[82]*subvol; 
  out[91] += 0.5*win[0]*fin[91]*subvol; 
  out[93] += 0.5*win[0]*fin[93]*subvol; 
  out[96] += 0.5*win[0]*fin[96]*subvol; 
  out[101] += 0.5*win[0]*fin[101]*subvol; 
  out[103] += 0.5*win[0]*fin[103]*subvol; 
  out[107] += 0.5*win[0]*fin[107]*subvol; 
  out[129] += 0.5*win[0]*fin[129]*subvol; 
  out[139] += 0.5*win[0]*fin[139]*subvol; 
  out[152] += 0.5*win[0]*fin[152]*subvol; 
  out[154] += 0.5*win[0]*fin[154]*subvol; 
  out[185] += 0.5*win[0]*fin[185]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxyvx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[3]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[5]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[6]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[15]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[19]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[21]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[23]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[24]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[43]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[51]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[54]*subvol); 
  atomicAdd(&out[12], 0.3535533905932737*win[0]*fin[57]*subvol); 
  atomicAdd(&out[13], 0.3535533905932737*win[0]*fin[59]*subvol); 
  atomicAdd(&out[14], 0.3535533905932737*win[0]*fin[93]*subvol); 
  atomicAdd(&out[15], 0.3535533905932737*win[0]*fin[103]*subvol); 
  atomicAdd(&out[16], 0.3535533905932737*win[0]*fin[107]*subvol); 
  atomicAdd(&out[17], 0.3535533905932737*win[0]*fin[154]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[3]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[5]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[6]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[15]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[19]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[21]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[23]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[24]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[43]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[51]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[54]*subvol; 
  out[12] += 0.3535533905932737*win[0]*fin[57]*subvol; 
  out[13] += 0.3535533905932737*win[0]*fin[59]*subvol; 
  out[14] += 0.3535533905932737*win[0]*fin[93]*subvol; 
  out[15] += 0.3535533905932737*win[0]*fin[103]*subvol; 
  out[16] += 0.3535533905932737*win[0]*fin[107]*subvol; 
  out[17] += 0.3535533905932737*win[0]*fin[154]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxyvx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
  atomicAdd(&out[3], win[0]*fin[3]*subvol); 
  atomicAdd(&out[5], win[0]*fin[5]*subvol); 
  atomicAdd(&out[6], win[0]*fin[6]*subvol); 
  atomicAdd(&out[15], win[0]*fin[15]*subvol); 
  atomicAdd(&out[19], win[0]*fin[19]*subvol); 
  atomicAdd(&out[21], win[0]*fin[21]*subvol); 
  atomicAdd(&out[23], win[0]*fin[23]*subvol); 
  atomicAdd(&out[24], win[0]*fin[24]*subvol); 
  atomicAdd(&out[43], win[0]*fin[43]*subvol); 
  atomicAdd(&out[51], win[0]*fin[51]*subvol); 
  atomicAdd(&out[54], win[0]*fin[54]*subvol); 
  atomicAdd(&out[57], win[0]*fin[57]*subvol); 
  atomicAdd(&out[59], win[0]*fin[59]*subvol); 
  atomicAdd(&out[93], win[0]*fin[93]*subvol); 
  atomicAdd(&out[103], win[0]*fin[103]*subvol); 
  atomicAdd(&out[107], win[0]*fin[107]*subvol); 
  atomicAdd(&out[154], win[0]*fin[154]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
  out[3] += win[0]*fin[3]*subvol; 
  out[5] += win[0]*fin[5]*subvol; 
  out[6] += win[0]*fin[6]*subvol; 
  out[15] += win[0]*fin[15]*subvol; 
  out[19] += win[0]*fin[19]*subvol; 
  out[21] += win[0]*fin[21]*subvol; 
  out[23] += win[0]*fin[23]*subvol; 
  out[24] += win[0]*fin[24]*subvol; 
  out[43] += win[0]*fin[43]*subvol; 
  out[51] += win[0]*fin[51]*subvol; 
  out[54] += win[0]*fin[54]*subvol; 
  out[57] += win[0]*fin[57]*subvol; 
  out[59] += win[0]*fin[59]*subvol; 
  out[93] += win[0]*fin[93]*subvol; 
  out[103] += win[0]*fin[103]*subvol; 
  out[107] += win[0]*fin[107]*subvol; 
  out[154] += win[0]*fin[154]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgzvx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.25*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.25*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.25*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.25*win[0]*fin[5]*subvol); 
  atomicAdd(&out[4], 0.25*win[0]*fin[6]*subvol); 
  atomicAdd(&out[5], 0.25*win[0]*fin[7]*subvol); 
  atomicAdd(&out[6], 0.25*win[0]*fin[13]*subvol); 
  atomicAdd(&out[7], 0.25*win[0]*fin[14]*subvol); 
  atomicAdd(&out[8], 0.25*win[0]*fin[17]*subvol); 
  atomicAdd(&out[9], 0.25*win[0]*fin[18]*subvol); 
  atomicAdd(&out[10], 0.25*win[0]*fin[21]*subvol); 
  atomicAdd(&out[11], 0.25*win[0]*fin[23]*subvol); 
  atomicAdd(&out[12], 0.25*win[0]*fin[24]*subvol); 
  atomicAdd(&out[13], 0.25*win[0]*fin[29]*subvol); 
  atomicAdd(&out[14], 0.25*win[0]*fin[35]*subvol); 
  atomicAdd(&out[15], 0.25*win[0]*fin[41]*subvol); 
  atomicAdd(&out[16], 0.25*win[0]*fin[42]*subvol); 
  atomicAdd(&out[17], 0.25*win[0]*fin[49]*subvol); 
  atomicAdd(&out[18], 0.25*win[0]*fin[50]*subvol); 
  atomicAdd(&out[19], 0.25*win[0]*fin[54]*subvol); 
  atomicAdd(&out[20], 0.25*win[0]*fin[55]*subvol); 
  atomicAdd(&out[21], 0.25*win[0]*fin[56]*subvol); 
  atomicAdd(&out[22], 0.25*win[0]*fin[59]*subvol); 
  atomicAdd(&out[23], 0.25*win[0]*fin[69]*subvol); 
  atomicAdd(&out[24], 0.25*win[0]*fin[81]*subvol); 
  atomicAdd(&out[25], 0.25*win[0]*fin[91]*subvol); 
  atomicAdd(&out[26], 0.25*win[0]*fin[92]*subvol); 
  atomicAdd(&out[27], 0.25*win[0]*fin[95]*subvol); 
  atomicAdd(&out[28], 0.25*win[0]*fin[101]*subvol); 
  atomicAdd(&out[29], 0.25*win[0]*fin[102]*subvol); 
  atomicAdd(&out[30], 0.25*win[0]*fin[107]*subvol); 
  atomicAdd(&out[31], 0.25*win[0]*fin[128]*subvol); 
  atomicAdd(&out[32], 0.25*win[0]*fin[138]*subvol); 
  atomicAdd(&out[33], 0.25*win[0]*fin[152]*subvol); 
  atomicAdd(&out[34], 0.25*win[0]*fin[153]*subvol); 
  atomicAdd(&out[35], 0.25*win[0]*fin[184]*subvol); 
#else
  out[0] += 0.25*fin[0]*win[0]*subvol; 
  out[1] += 0.25*win[0]*fin[1]*subvol; 
  out[2] += 0.25*win[0]*fin[2]*subvol; 
  out[3] += 0.25*win[0]*fin[5]*subvol; 
  out[4] += 0.25*win[0]*fin[6]*subvol; 
  out[5] += 0.25*win[0]*fin[7]*subvol; 
  out[6] += 0.25*win[0]*fin[13]*subvol; 
  out[7] += 0.25*win[0]*fin[14]*subvol; 
  out[8] += 0.25*win[0]*fin[17]*subvol; 
  out[9] += 0.25*win[0]*fin[18]*subvol; 
  out[10] += 0.25*win[0]*fin[21]*subvol; 
  out[11] += 0.25*win[0]*fin[23]*subvol; 
  out[12] += 0.25*win[0]*fin[24]*subvol; 
  out[13] += 0.25*win[0]*fin[29]*subvol; 
  out[14] += 0.25*win[0]*fin[35]*subvol; 
  out[15] += 0.25*win[0]*fin[41]*subvol; 
  out[16] += 0.25*win[0]*fin[42]*subvol; 
  out[17] += 0.25*win[0]*fin[49]*subvol; 
  out[18] += 0.25*win[0]*fin[50]*subvol; 
  out[19] += 0.25*win[0]*fin[54]*subvol; 
  out[20] += 0.25*win[0]*fin[55]*subvol; 
  out[21] += 0.25*win[0]*fin[56]*subvol; 
  out[22] += 0.25*win[0]*fin[59]*subvol; 
  out[23] += 0.25*win[0]*fin[69]*subvol; 
  out[24] += 0.25*win[0]*fin[81]*subvol; 
  out[25] += 0.25*win[0]*fin[91]*subvol; 
  out[26] += 0.25*win[0]*fin[92]*subvol; 
  out[27] += 0.25*win[0]*fin[95]*subvol; 
  out[28] += 0.25*win[0]*fin[101]*subvol; 
  out[29] += 0.25*win[0]*fin[102]*subvol; 
  out[30] += 0.25*win[0]*fin[107]*subvol; 
  out[31] += 0.25*win[0]*fin[128]*subvol; 
  out[32] += 0.25*win[0]*fin[138]*subvol; 
  out[33] += 0.25*win[0]*fin[152]*subvol; 
  out[34] += 0.25*win[0]*fin[153]*subvol; 
  out[35] += 0.25*win[0]*fin[184]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgzvx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[5]*subvol); 
  atomicAdd(&out[6], 0.5*win[0]*fin[6]*subvol); 
  atomicAdd(&out[7], 0.5*win[0]*fin[7]*subvol); 
  atomicAdd(&out[13], 0.5*win[0]*fin[13]*subvol); 
  atomicAdd(&out[14], 0.5*win[0]*fin[14]*subvol); 
  atomicAdd(&out[17], 0.5*win[0]*fin[17]*subvol); 
  atomicAdd(&out[18], 0.5*win[0]*fin[18]*subvol); 
  atomicAdd(&out[21], 0.5*win[0]*fin[21]*subvol); 
  atomicAdd(&out[23], 0.5*win[0]*fin[23]*subvol); 
  atomicAdd(&out[24], 0.5*win[0]*fin[24]*subvol); 
  atomicAdd(&out[29], 0.5*win[0]*fin[29]*subvol); 
  atomicAdd(&out[35], 0.5*win[0]*fin[35]*subvol); 
  atomicAdd(&out[41], 0.5*win[0]*fin[41]*subvol); 
  atomicAdd(&out[42], 0.5*win[0]*fin[42]*subvol); 
  atomicAdd(&out[49], 0.5*win[0]*fin[49]*subvol); 
  atomicAdd(&out[50], 0.5*win[0]*fin[50]*subvol); 
  atomicAdd(&out[54], 0.5*win[0]*fin[54]*subvol); 
  atomicAdd(&out[55], 0.5*win[0]*fin[55]*subvol); 
  atomicAdd(&out[56], 0.5*win[0]*fin[56]*subvol); 
  atomicAdd(&out[59], 0.5*win[0]*fin[59]*subvol); 
  atomicAdd(&out[69], 0.5*win[0]*fin[69]*subvol); 
  atomicAdd(&out[81], 0.5*win[0]*fin[81]*subvol); 
  atomicAdd(&out[91], 0.5*win[0]*fin[91]*subvol); 
  atomicAdd(&out[92], 0.5*win[0]*fin[92]*subvol); 
  atomicAdd(&out[95], 0.5*win[0]*fin[95]*subvol); 
  atomicAdd(&out[101], 0.5*win[0]*fin[101]*subvol); 
  atomicAdd(&out[102], 0.5*win[0]*fin[102]*subvol); 
  atomicAdd(&out[107], 0.5*win[0]*fin[107]*subvol); 
  atomicAdd(&out[128], 0.5*win[0]*fin[128]*subvol); 
  atomicAdd(&out[138], 0.5*win[0]*fin[138]*subvol); 
  atomicAdd(&out[152], 0.5*win[0]*fin[152]*subvol); 
  atomicAdd(&out[153], 0.5*win[0]*fin[153]*subvol); 
  atomicAdd(&out[184], 0.5*win[0]*fin[184]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[2] += 0.5*win[0]*fin[2]*subvol; 
  out[5] += 0.5*win[0]*fin[5]*subvol; 
  out[6] += 0.5*win[0]*fin[6]*subvol; 
  out[7] += 0.5*win[0]*fin[7]*subvol; 
  out[13] += 0.5*win[0]*fin[13]*subvol; 
  out[14] += 0.5*win[0]*fin[14]*subvol; 
  out[17] += 0.5*win[0]*fin[17]*subvol; 
  out[18] += 0.5*win[0]*fin[18]*subvol; 
  out[21] += 0.5*win[0]*fin[21]*subvol; 
  out[23] += 0.5*win[0]*fin[23]*subvol; 
  out[24] += 0.5*win[0]*fin[24]*subvol; 
  out[29] += 0.5*win[0]*fin[29]*subvol; 
  out[35] += 0.5*win[0]*fin[35]*subvol; 
  out[41] += 0.5*win[0]*fin[41]*subvol; 
  out[42] += 0.5*win[0]*fin[42]*subvol; 
  out[49] += 0.5*win[0]*fin[49]*subvol; 
  out[50] += 0.5*win[0]*fin[50]*subvol; 
  out[54] += 0.5*win[0]*fin[54]*subvol; 
  out[55] += 0.5*win[0]*fin[55]*subvol; 
  out[56] += 0.5*win[0]*fin[56]*subvol; 
  out[59] += 0.5*win[0]*fin[59]*subvol; 
  out[69] += 0.5*win[0]*fin[69]*subvol; 
  out[81] += 0.5*win[0]*fin[81]*subvol; 
  out[91] += 0.5*win[0]*fin[91]*subvol; 
  out[92] += 0.5*win[0]*fin[92]*subvol; 
  out[95] += 0.5*win[0]*fin[95]*subvol; 
  out[101] += 0.5*win[0]*fin[101]*subvol; 
  out[102] += 0.5*win[0]*fin[102]*subvol; 
  out[107] += 0.5*win[0]*fin[107]*subvol; 
  out[128] += 0.5*win[0]*fin[128]*subvol; 
  out[138] += 0.5*win[0]*fin[138]*subvol; 
  out[152] += 0.5*win[0]*fin[152]*subvol; 
  out[153] += 0.5*win[0]*fin[153]*subvol; 
  out[184] += 0.5*win[0]*fin[184]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxzvx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[5]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[6]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[14]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[18]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[21]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[23]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[24]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[42]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[50]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[54]*subvol); 
  atomicAdd(&out[12], 0.3535533905932737*win[0]*fin[56]*subvol); 
  atomicAdd(&out[13], 0.3535533905932737*win[0]*fin[59]*subvol); 
  atomicAdd(&out[14], 0.3535533905932737*win[0]*fin[92]*subvol); 
  atomicAdd(&out[15], 0.3535533905932737*win[0]*fin[102]*subvol); 
  atomicAdd(&out[16], 0.3535533905932737*win[0]*fin[107]*subvol); 
  atomicAdd(&out[17], 0.3535533905932737*win[0]*fin[153]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[2]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[5]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[6]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[14]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[18]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[21]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[23]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[24]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[42]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[50]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[54]*subvol; 
  out[12] += 0.3535533905932737*win[0]*fin[56]*subvol; 
  out[13] += 0.3535533905932737*win[0]*fin[59]*subvol; 
  out[14] += 0.3535533905932737*win[0]*fin[92]*subvol; 
  out[15] += 0.3535533905932737*win[0]*fin[102]*subvol; 
  out[16] += 0.3535533905932737*win[0]*fin[107]*subvol; 
  out[17] += 0.3535533905932737*win[0]*fin[153]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxzvx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], win[0]*fin[2]*subvol); 
  atomicAdd(&out[5], win[0]*fin[5]*subvol); 
  atomicAdd(&out[6], win[0]*fin[6]*subvol); 
  atomicAdd(&out[14], win[0]*fin[14]*subvol); 
  atomicAdd(&out[18], win[0]*fin[18]*subvol); 
  atomicAdd(&out[21], win[0]*fin[21]*subvol); 
  atomicAdd(&out[23], win[0]*fin[23]*subvol); 
  atomicAdd(&out[24], win[0]*fin[24]*subvol); 
  atomicAdd(&out[42], win[0]*fin[42]*subvol); 
  atomicAdd(&out[50], win[0]*fin[50]*subvol); 
  atomicAdd(&out[54], win[0]*fin[54]*subvol); 
  atomicAdd(&out[56], win[0]*fin[56]*subvol); 
  atomicAdd(&out[59], win[0]*fin[59]*subvol); 
  atomicAdd(&out[92], win[0]*fin[92]*subvol); 
  atomicAdd(&out[102], win[0]*fin[102]*subvol); 
  atomicAdd(&out[107], win[0]*fin[107]*subvol); 
  atomicAdd(&out[153], win[0]*fin[153]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
  out[2] += win[0]*fin[2]*subvol; 
  out[5] += win[0]*fin[5]*subvol; 
  out[6] += win[0]*fin[6]*subvol; 
  out[14] += win[0]*fin[14]*subvol; 
  out[18] += win[0]*fin[18]*subvol; 
  out[21] += win[0]*fin[21]*subvol; 
  out[23] += win[0]*fin[23]*subvol; 
  out[24] += win[0]*fin[24]*subvol; 
  out[42] += win[0]*fin[42]*subvol; 
  out[50] += win[0]*fin[50]*subvol; 
  out[54] += win[0]*fin[54]*subvol; 
  out[56] += win[0]*fin[56]*subvol; 
  out[59] += win[0]*fin[59]*subvol; 
  out[92] += win[0]*fin[92]*subvol; 
  out[102] += win[0]*fin[102]*subvol; 
  out[107] += win[0]*fin[107]*subvol; 
  out[153] += win[0]*fin[153]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgyzvx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[5]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[6]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[13]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[17]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[21]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[23]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[24]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[41]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[49]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[54]*subvol); 
  atomicAdd(&out[12], 0.3535533905932737*win[0]*fin[55]*subvol); 
  atomicAdd(&out[13], 0.3535533905932737*win[0]*fin[59]*subvol); 
  atomicAdd(&out[14], 0.3535533905932737*win[0]*fin[91]*subvol); 
  atomicAdd(&out[15], 0.3535533905932737*win[0]*fin[101]*subvol); 
  atomicAdd(&out[16], 0.3535533905932737*win[0]*fin[107]*subvol); 
  atomicAdd(&out[17], 0.3535533905932737*win[0]*fin[152]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[1]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[5]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[6]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[13]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[17]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[21]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[23]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[24]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[41]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[49]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[54]*subvol; 
  out[12] += 0.3535533905932737*win[0]*fin[55]*subvol; 
  out[13] += 0.3535533905932737*win[0]*fin[59]*subvol; 
  out[14] += 0.3535533905932737*win[0]*fin[91]*subvol; 
  out[15] += 0.3535533905932737*win[0]*fin[101]*subvol; 
  out[16] += 0.3535533905932737*win[0]*fin[107]*subvol; 
  out[17] += 0.3535533905932737*win[0]*fin[152]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgyzvx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], win[0]*fin[1]*subvol); 
  atomicAdd(&out[5], win[0]*fin[5]*subvol); 
  atomicAdd(&out[6], win[0]*fin[6]*subvol); 
  atomicAdd(&out[13], win[0]*fin[13]*subvol); 
  atomicAdd(&out[17], win[0]*fin[17]*subvol); 
  atomicAdd(&out[21], win[0]*fin[21]*subvol); 
  atomicAdd(&out[23], win[0]*fin[23]*subvol); 
  atomicAdd(&out[24], win[0]*fin[24]*subvol); 
  atomicAdd(&out[41], win[0]*fin[41]*subvol); 
  atomicAdd(&out[49], win[0]*fin[49]*subvol); 
  atomicAdd(&out[54], win[0]*fin[54]*subvol); 
  atomicAdd(&out[55], win[0]*fin[55]*subvol); 
  atomicAdd(&out[59], win[0]*fin[59]*subvol); 
  atomicAdd(&out[91], win[0]*fin[91]*subvol); 
  atomicAdd(&out[101], win[0]*fin[101]*subvol); 
  atomicAdd(&out[107], win[0]*fin[107]*subvol); 
  atomicAdd(&out[152], win[0]*fin[152]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
  out[1] += win[0]*fin[1]*subvol; 
  out[5] += win[0]*fin[5]*subvol; 
  out[6] += win[0]*fin[6]*subvol; 
  out[13] += win[0]*fin[13]*subvol; 
  out[17] += win[0]*fin[17]*subvol; 
  out[21] += win[0]*fin[21]*subvol; 
  out[23] += win[0]*fin[23]*subvol; 
  out[24] += win[0]*fin[24]*subvol; 
  out[41] += win[0]*fin[41]*subvol; 
  out[49] += win[0]*fin[49]*subvol; 
  out[54] += win[0]*fin[54]*subvol; 
  out[55] += win[0]*fin[55]*subvol; 
  out[59] += win[0]*fin[59]*subvol; 
  out[91] += win[0]*fin[91]*subvol; 
  out[101] += win[0]*fin[101]*subvol; 
  out[107] += win[0]*fin[107]*subvol; 
  out[152] += win[0]*fin[152]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxyzvx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[5]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[6]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[21]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[23]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[24]*subvol); 
  atomicAdd(&out[6], 0.5*win[0]*fin[54]*subvol); 
  atomicAdd(&out[7], 0.5*win[0]*fin[59]*subvol); 
  atomicAdd(&out[8], 0.5*win[0]*fin[107]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[5]*subvol; 
  out[2] += 0.5*win[0]*fin[6]*subvol; 
  out[3] += 0.5*win[0]*fin[21]*subvol; 
  out[4] += 0.5*win[0]*fin[23]*subvol; 
  out[5] += 0.5*win[0]*fin[24]*subvol; 
  out[6] += 0.5*win[0]*fin[54]*subvol; 
  out[7] += 0.5*win[0]*fin[59]*subvol; 
  out[8] += 0.5*win[0]*fin[107]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxyzvx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.0*fin[0]*win[0]*subvol); 
  atomicAdd(&out[5], 2.0*win[0]*fin[5]*subvol); 
  atomicAdd(&out[6], 2.0*win[0]*fin[6]*subvol); 
  atomicAdd(&out[21], 2.0*win[0]*fin[21]*subvol); 
  atomicAdd(&out[23], 2.0*win[0]*fin[23]*subvol); 
  atomicAdd(&out[24], 2.0*win[0]*fin[24]*subvol); 
  atomicAdd(&out[54], 2.0*win[0]*fin[54]*subvol); 
  atomicAdd(&out[59], 2.0*win[0]*fin[59]*subvol); 
  atomicAdd(&out[107], 2.0*win[0]*fin[107]*subvol); 
#else
  out[0] += 2.0*fin[0]*win[0]*subvol; 
  out[5] += 2.0*win[0]*fin[5]*subvol; 
  out[6] += 2.0*win[0]*fin[6]*subvol; 
  out[21] += 2.0*win[0]*fin[21]*subvol; 
  out[23] += 2.0*win[0]*fin[23]*subvol; 
  out[24] += 2.0*win[0]*fin[24]*subvol; 
  out[54] += 2.0*win[0]*fin[54]*subvol; 
  out[59] += 2.0*win[0]*fin[59]*subvol; 
  out[107] += 2.0*win[0]*fin[107]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.1767766952966368*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.1767766952966368*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.1767766952966368*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.1767766952966368*win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], 0.1767766952966368*win[0]*fin[4]*subvol); 
  atomicAdd(&out[5], 0.1767766952966368*win[0]*fin[6]*subvol); 
  atomicAdd(&out[6], 0.1767766952966368*win[0]*fin[7]*subvol); 
  atomicAdd(&out[7], 0.1767766952966368*win[0]*fin[8]*subvol); 
  atomicAdd(&out[8], 0.1767766952966368*win[0]*fin[9]*subvol); 
  atomicAdd(&out[9], 0.1767766952966368*win[0]*fin[10]*subvol); 
  atomicAdd(&out[10], 0.1767766952966368*win[0]*fin[11]*subvol); 
  atomicAdd(&out[11], 0.1767766952966368*win[0]*fin[12]*subvol); 
  atomicAdd(&out[12], 0.1767766952966368*win[0]*fin[17]*subvol); 
  atomicAdd(&out[13], 0.1767766952966368*win[0]*fin[18]*subvol); 
  atomicAdd(&out[14], 0.1767766952966368*win[0]*fin[19]*subvol); 
  atomicAdd(&out[15], 0.1767766952966368*win[0]*fin[20]*subvol); 
  atomicAdd(&out[16], 0.1767766952966368*win[0]*fin[22]*subvol); 
  atomicAdd(&out[17], 0.1767766952966368*win[0]*fin[24]*subvol); 
  atomicAdd(&out[18], 0.1767766952966368*win[0]*fin[25]*subvol); 
  atomicAdd(&out[19], 0.1767766952966368*win[0]*fin[26]*subvol); 
  atomicAdd(&out[20], 0.1767766952966368*win[0]*fin[27]*subvol); 
  atomicAdd(&out[21], 0.1767766952966368*win[0]*fin[28]*subvol); 
  atomicAdd(&out[22], 0.1767766952966368*win[0]*fin[35]*subvol); 
  atomicAdd(&out[23], 0.1767766952966368*win[0]*fin[36]*subvol); 
  atomicAdd(&out[24], 0.1767766952966368*win[0]*fin[37]*subvol); 
  atomicAdd(&out[25], 0.1767766952966368*win[0]*fin[38]*subvol); 
  atomicAdd(&out[26], 0.1767766952966368*win[0]*fin[39]*subvol); 
  atomicAdd(&out[27], 0.1767766952966368*win[0]*fin[40]*subvol); 
  atomicAdd(&out[28], 0.1767766952966368*win[0]*fin[45]*subvol); 
  atomicAdd(&out[29], 0.1767766952966368*win[0]*fin[46]*subvol); 
  atomicAdd(&out[30], 0.1767766952966368*win[0]*fin[47]*subvol); 
  atomicAdd(&out[31], 0.1767766952966368*win[0]*fin[53]*subvol); 
  atomicAdd(&out[32], 0.1767766952966368*win[0]*fin[55]*subvol); 
  atomicAdd(&out[33], 0.1767766952966368*win[0]*fin[56]*subvol); 
  atomicAdd(&out[34], 0.1767766952966368*win[0]*fin[57]*subvol); 
  atomicAdd(&out[35], 0.1767766952966368*win[0]*fin[58]*subvol); 
  atomicAdd(&out[36], 0.1767766952966368*win[0]*fin[60]*subvol); 
  atomicAdd(&out[37], 0.1767766952966368*win[0]*fin[65]*subvol); 
  atomicAdd(&out[38], 0.1767766952966368*win[0]*fin[66]*subvol); 
  atomicAdd(&out[39], 0.1767766952966368*win[0]*fin[67]*subvol); 
  atomicAdd(&out[40], 0.1767766952966368*win[0]*fin[68]*subvol); 
  atomicAdd(&out[41], 0.1767766952966368*win[0]*fin[75]*subvol); 
  atomicAdd(&out[42], 0.1767766952966368*win[0]*fin[76]*subvol); 
  atomicAdd(&out[43], 0.1767766952966368*win[0]*fin[77]*subvol); 
  atomicAdd(&out[44], 0.1767766952966368*win[0]*fin[87]*subvol); 
  atomicAdd(&out[45], 0.1767766952966368*win[0]*fin[88]*subvol); 
  atomicAdd(&out[46], 0.1767766952966368*win[0]*fin[89]*subvol); 
  atomicAdd(&out[47], 0.1767766952966368*win[0]*fin[95]*subvol); 
  atomicAdd(&out[48], 0.1767766952966368*win[0]*fin[96]*subvol); 
  atomicAdd(&out[49], 0.1767766952966368*win[0]*fin[97]*subvol); 
  atomicAdd(&out[50], 0.1767766952966368*win[0]*fin[98]*subvol); 
  atomicAdd(&out[51], 0.1767766952966368*win[0]*fin[99]*subvol); 
  atomicAdd(&out[52], 0.1767766952966368*win[0]*fin[100]*subvol); 
  atomicAdd(&out[53], 0.1767766952966368*win[0]*fin[106]*subvol); 
  atomicAdd(&out[54], 0.1767766952966368*win[0]*fin[109]*subvol); 
  atomicAdd(&out[55], 0.1767766952966368*win[0]*fin[114]*subvol); 
  atomicAdd(&out[56], 0.1767766952966368*win[0]*fin[122]*subvol); 
  atomicAdd(&out[57], 0.1767766952966368*win[0]*fin[123]*subvol); 
  atomicAdd(&out[58], 0.1767766952966368*win[0]*fin[124]*subvol); 
  atomicAdd(&out[59], 0.1767766952966368*win[0]*fin[134]*subvol); 
  atomicAdd(&out[60], 0.1767766952966368*win[0]*fin[135]*subvol); 
  atomicAdd(&out[61], 0.1767766952966368*win[0]*fin[136]*subvol); 
  atomicAdd(&out[62], 0.1767766952966368*win[0]*fin[137]*subvol); 
  atomicAdd(&out[63], 0.1767766952966368*win[0]*fin[148]*subvol); 
  atomicAdd(&out[64], 0.1767766952966368*win[0]*fin[149]*subvol); 
  atomicAdd(&out[65], 0.1767766952966368*win[0]*fin[150]*subvol); 
  atomicAdd(&out[66], 0.1767766952966368*win[0]*fin[159]*subvol); 
  atomicAdd(&out[67], 0.1767766952966368*win[0]*fin[167]*subvol); 
  atomicAdd(&out[68], 0.1767766952966368*win[0]*fin[178]*subvol); 
  atomicAdd(&out[69], 0.1767766952966368*win[0]*fin[179]*subvol); 
  atomicAdd(&out[70], 0.1767766952966368*win[0]*fin[180]*subvol); 
  atomicAdd(&out[71], 0.1767766952966368*win[0]*fin[198]*subvol); 
#else
  out[0] += 0.1767766952966368*fin[0]*win[0]*subvol; 
  out[1] += 0.1767766952966368*win[0]*fin[1]*subvol; 
  out[2] += 0.1767766952966368*win[0]*fin[2]*subvol; 
  out[3] += 0.1767766952966368*win[0]*fin[3]*subvol; 
  out[4] += 0.1767766952966368*win[0]*fin[4]*subvol; 
  out[5] += 0.1767766952966368*win[0]*fin[6]*subvol; 
  out[6] += 0.1767766952966368*win[0]*fin[7]*subvol; 
  out[7] += 0.1767766952966368*win[0]*fin[8]*subvol; 
  out[8] += 0.1767766952966368*win[0]*fin[9]*subvol; 
  out[9] += 0.1767766952966368*win[0]*fin[10]*subvol; 
  out[10] += 0.1767766952966368*win[0]*fin[11]*subvol; 
  out[11] += 0.1767766952966368*win[0]*fin[12]*subvol; 
  out[12] += 0.1767766952966368*win[0]*fin[17]*subvol; 
  out[13] += 0.1767766952966368*win[0]*fin[18]*subvol; 
  out[14] += 0.1767766952966368*win[0]*fin[19]*subvol; 
  out[15] += 0.1767766952966368*win[0]*fin[20]*subvol; 
  out[16] += 0.1767766952966368*win[0]*fin[22]*subvol; 
  out[17] += 0.1767766952966368*win[0]*fin[24]*subvol; 
  out[18] += 0.1767766952966368*win[0]*fin[25]*subvol; 
  out[19] += 0.1767766952966368*win[0]*fin[26]*subvol; 
  out[20] += 0.1767766952966368*win[0]*fin[27]*subvol; 
  out[21] += 0.1767766952966368*win[0]*fin[28]*subvol; 
  out[22] += 0.1767766952966368*win[0]*fin[35]*subvol; 
  out[23] += 0.1767766952966368*win[0]*fin[36]*subvol; 
  out[24] += 0.1767766952966368*win[0]*fin[37]*subvol; 
  out[25] += 0.1767766952966368*win[0]*fin[38]*subvol; 
  out[26] += 0.1767766952966368*win[0]*fin[39]*subvol; 
  out[27] += 0.1767766952966368*win[0]*fin[40]*subvol; 
  out[28] += 0.1767766952966368*win[0]*fin[45]*subvol; 
  out[29] += 0.1767766952966368*win[0]*fin[46]*subvol; 
  out[30] += 0.1767766952966368*win[0]*fin[47]*subvol; 
  out[31] += 0.1767766952966368*win[0]*fin[53]*subvol; 
  out[32] += 0.1767766952966368*win[0]*fin[55]*subvol; 
  out[33] += 0.1767766952966368*win[0]*fin[56]*subvol; 
  out[34] += 0.1767766952966368*win[0]*fin[57]*subvol; 
  out[35] += 0.1767766952966368*win[0]*fin[58]*subvol; 
  out[36] += 0.1767766952966368*win[0]*fin[60]*subvol; 
  out[37] += 0.1767766952966368*win[0]*fin[65]*subvol; 
  out[38] += 0.1767766952966368*win[0]*fin[66]*subvol; 
  out[39] += 0.1767766952966368*win[0]*fin[67]*subvol; 
  out[40] += 0.1767766952966368*win[0]*fin[68]*subvol; 
  out[41] += 0.1767766952966368*win[0]*fin[75]*subvol; 
  out[42] += 0.1767766952966368*win[0]*fin[76]*subvol; 
  out[43] += 0.1767766952966368*win[0]*fin[77]*subvol; 
  out[44] += 0.1767766952966368*win[0]*fin[87]*subvol; 
  out[45] += 0.1767766952966368*win[0]*fin[88]*subvol; 
  out[46] += 0.1767766952966368*win[0]*fin[89]*subvol; 
  out[47] += 0.1767766952966368*win[0]*fin[95]*subvol; 
  out[48] += 0.1767766952966368*win[0]*fin[96]*subvol; 
  out[49] += 0.1767766952966368*win[0]*fin[97]*subvol; 
  out[50] += 0.1767766952966368*win[0]*fin[98]*subvol; 
  out[51] += 0.1767766952966368*win[0]*fin[99]*subvol; 
  out[52] += 0.1767766952966368*win[0]*fin[100]*subvol; 
  out[53] += 0.1767766952966368*win[0]*fin[106]*subvol; 
  out[54] += 0.1767766952966368*win[0]*fin[109]*subvol; 
  out[55] += 0.1767766952966368*win[0]*fin[114]*subvol; 
  out[56] += 0.1767766952966368*win[0]*fin[122]*subvol; 
  out[57] += 0.1767766952966368*win[0]*fin[123]*subvol; 
  out[58] += 0.1767766952966368*win[0]*fin[124]*subvol; 
  out[59] += 0.1767766952966368*win[0]*fin[134]*subvol; 
  out[60] += 0.1767766952966368*win[0]*fin[135]*subvol; 
  out[61] += 0.1767766952966368*win[0]*fin[136]*subvol; 
  out[62] += 0.1767766952966368*win[0]*fin[137]*subvol; 
  out[63] += 0.1767766952966368*win[0]*fin[148]*subvol; 
  out[64] += 0.1767766952966368*win[0]*fin[149]*subvol; 
  out[65] += 0.1767766952966368*win[0]*fin[150]*subvol; 
  out[66] += 0.1767766952966368*win[0]*fin[159]*subvol; 
  out[67] += 0.1767766952966368*win[0]*fin[167]*subvol; 
  out[68] += 0.1767766952966368*win[0]*fin[178]*subvol; 
  out[69] += 0.1767766952966368*win[0]*fin[179]*subvol; 
  out[70] += 0.1767766952966368*win[0]*fin[180]*subvol; 
  out[71] += 0.1767766952966368*win[0]*fin[198]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.25*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.25*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.25*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.25*win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], 0.25*win[0]*fin[4]*subvol); 
  atomicAdd(&out[6], 0.25*win[0]*fin[6]*subvol); 
  atomicAdd(&out[7], 0.25*win[0]*fin[7]*subvol); 
  atomicAdd(&out[8], 0.25*win[0]*fin[8]*subvol); 
  atomicAdd(&out[9], 0.25*win[0]*fin[9]*subvol); 
  atomicAdd(&out[10], 0.25*win[0]*fin[10]*subvol); 
  atomicAdd(&out[11], 0.25*win[0]*fin[11]*subvol); 
  atomicAdd(&out[12], 0.25*win[0]*fin[12]*subvol); 
  atomicAdd(&out[17], 0.25*win[0]*fin[17]*subvol); 
  atomicAdd(&out[18], 0.25*win[0]*fin[18]*subvol); 
  atomicAdd(&out[19], 0.25*win[0]*fin[19]*subvol); 
  atomicAdd(&out[20], 0.25*win[0]*fin[20]*subvol); 
  atomicAdd(&out[22], 0.25*win[0]*fin[22]*subvol); 
  atomicAdd(&out[24], 0.25*win[0]*fin[24]*subvol); 
  atomicAdd(&out[25], 0.25*win[0]*fin[25]*subvol); 
  atomicAdd(&out[26], 0.25*win[0]*fin[26]*subvol); 
  atomicAdd(&out[27], 0.25*win[0]*fin[27]*subvol); 
  atomicAdd(&out[28], 0.25*win[0]*fin[28]*subvol); 
  atomicAdd(&out[35], 0.25*win[0]*fin[35]*subvol); 
  atomicAdd(&out[36], 0.25*win[0]*fin[36]*subvol); 
  atomicAdd(&out[37], 0.25*win[0]*fin[37]*subvol); 
  atomicAdd(&out[38], 0.25*win[0]*fin[38]*subvol); 
  atomicAdd(&out[39], 0.25*win[0]*fin[39]*subvol); 
  atomicAdd(&out[40], 0.25*win[0]*fin[40]*subvol); 
  atomicAdd(&out[45], 0.25*win[0]*fin[45]*subvol); 
  atomicAdd(&out[46], 0.25*win[0]*fin[46]*subvol); 
  atomicAdd(&out[47], 0.25*win[0]*fin[47]*subvol); 
  atomicAdd(&out[53], 0.25*win[0]*fin[53]*subvol); 
  atomicAdd(&out[55], 0.25*win[0]*fin[55]*subvol); 
  atomicAdd(&out[56], 0.25*win[0]*fin[56]*subvol); 
  atomicAdd(&out[57], 0.25*win[0]*fin[57]*subvol); 
  atomicAdd(&out[58], 0.25*win[0]*fin[58]*subvol); 
  atomicAdd(&out[60], 0.25*win[0]*fin[60]*subvol); 
  atomicAdd(&out[65], 0.25*win[0]*fin[65]*subvol); 
  atomicAdd(&out[66], 0.25*win[0]*fin[66]*subvol); 
  atomicAdd(&out[67], 0.25*win[0]*fin[67]*subvol); 
  atomicAdd(&out[68], 0.25*win[0]*fin[68]*subvol); 
  atomicAdd(&out[75], 0.25*win[0]*fin[75]*subvol); 
  atomicAdd(&out[76], 0.25*win[0]*fin[76]*subvol); 
  atomicAdd(&out[77], 0.25*win[0]*fin[77]*subvol); 
  atomicAdd(&out[87], 0.25*win[0]*fin[87]*subvol); 
  atomicAdd(&out[88], 0.25*win[0]*fin[88]*subvol); 
  atomicAdd(&out[89], 0.25*win[0]*fin[89]*subvol); 
  atomicAdd(&out[95], 0.25*win[0]*fin[95]*subvol); 
  atomicAdd(&out[96], 0.25*win[0]*fin[96]*subvol); 
  atomicAdd(&out[97], 0.25*win[0]*fin[97]*subvol); 
  atomicAdd(&out[98], 0.25*win[0]*fin[98]*subvol); 
  atomicAdd(&out[99], 0.25*win[0]*fin[99]*subvol); 
  atomicAdd(&out[100], 0.25*win[0]*fin[100]*subvol); 
  atomicAdd(&out[106], 0.25*win[0]*fin[106]*subvol); 
  atomicAdd(&out[109], 0.25*win[0]*fin[109]*subvol); 
  atomicAdd(&out[114], 0.25*win[0]*fin[114]*subvol); 
  atomicAdd(&out[122], 0.25*win[0]*fin[122]*subvol); 
  atomicAdd(&out[123], 0.25*win[0]*fin[123]*subvol); 
  atomicAdd(&out[124], 0.25*win[0]*fin[124]*subvol); 
  atomicAdd(&out[134], 0.25*win[0]*fin[134]*subvol); 
  atomicAdd(&out[135], 0.25*win[0]*fin[135]*subvol); 
  atomicAdd(&out[136], 0.25*win[0]*fin[136]*subvol); 
  atomicAdd(&out[137], 0.25*win[0]*fin[137]*subvol); 
  atomicAdd(&out[148], 0.25*win[0]*fin[148]*subvol); 
  atomicAdd(&out[149], 0.25*win[0]*fin[149]*subvol); 
  atomicAdd(&out[150], 0.25*win[0]*fin[150]*subvol); 
  atomicAdd(&out[159], 0.25*win[0]*fin[159]*subvol); 
  atomicAdd(&out[167], 0.25*win[0]*fin[167]*subvol); 
  atomicAdd(&out[178], 0.25*win[0]*fin[178]*subvol); 
  atomicAdd(&out[179], 0.25*win[0]*fin[179]*subvol); 
  atomicAdd(&out[180], 0.25*win[0]*fin[180]*subvol); 
  atomicAdd(&out[198], 0.25*win[0]*fin[198]*subvol); 
#else
  out[0] += 0.25*fin[0]*win[0]*subvol; 
  out[1] += 0.25*win[0]*fin[1]*subvol; 
  out[2] += 0.25*win[0]*fin[2]*subvol; 
  out[3] += 0.25*win[0]*fin[3]*subvol; 
  out[4] += 0.25*win[0]*fin[4]*subvol; 
  out[6] += 0.25*win[0]*fin[6]*subvol; 
  out[7] += 0.25*win[0]*fin[7]*subvol; 
  out[8] += 0.25*win[0]*fin[8]*subvol; 
  out[9] += 0.25*win[0]*fin[9]*subvol; 
  out[10] += 0.25*win[0]*fin[10]*subvol; 
  out[11] += 0.25*win[0]*fin[11]*subvol; 
  out[12] += 0.25*win[0]*fin[12]*subvol; 
  out[17] += 0.25*win[0]*fin[17]*subvol; 
  out[18] += 0.25*win[0]*fin[18]*subvol; 
  out[19] += 0.25*win[0]*fin[19]*subvol; 
  out[20] += 0.25*win[0]*fin[20]*subvol; 
  out[22] += 0.25*win[0]*fin[22]*subvol; 
  out[24] += 0.25*win[0]*fin[24]*subvol; 
  out[25] += 0.25*win[0]*fin[25]*subvol; 
  out[26] += 0.25*win[0]*fin[26]*subvol; 
  out[27] += 0.25*win[0]*fin[27]*subvol; 
  out[28] += 0.25*win[0]*fin[28]*subvol; 
  out[35] += 0.25*win[0]*fin[35]*subvol; 
  out[36] += 0.25*win[0]*fin[36]*subvol; 
  out[37] += 0.25*win[0]*fin[37]*subvol; 
  out[38] += 0.25*win[0]*fin[38]*subvol; 
  out[39] += 0.25*win[0]*fin[39]*subvol; 
  out[40] += 0.25*win[0]*fin[40]*subvol; 
  out[45] += 0.25*win[0]*fin[45]*subvol; 
  out[46] += 0.25*win[0]*fin[46]*subvol; 
  out[47] += 0.25*win[0]*fin[47]*subvol; 
  out[53] += 0.25*win[0]*fin[53]*subvol; 
  out[55] += 0.25*win[0]*fin[55]*subvol; 
  out[56] += 0.25*win[0]*fin[56]*subvol; 
  out[57] += 0.25*win[0]*fin[57]*subvol; 
  out[58] += 0.25*win[0]*fin[58]*subvol; 
  out[60] += 0.25*win[0]*fin[60]*subvol; 
  out[65] += 0.25*win[0]*fin[65]*subvol; 
  out[66] += 0.25*win[0]*fin[66]*subvol; 
  out[67] += 0.25*win[0]*fin[67]*subvol; 
  out[68] += 0.25*win[0]*fin[68]*subvol; 
  out[75] += 0.25*win[0]*fin[75]*subvol; 
  out[76] += 0.25*win[0]*fin[76]*subvol; 
  out[77] += 0.25*win[0]*fin[77]*subvol; 
  out[87] += 0.25*win[0]*fin[87]*subvol; 
  out[88] += 0.25*win[0]*fin[88]*subvol; 
  out[89] += 0.25*win[0]*fin[89]*subvol; 
  out[95] += 0.25*win[0]*fin[95]*subvol; 
  out[96] += 0.25*win[0]*fin[96]*subvol; 
  out[97] += 0.25*win[0]*fin[97]*subvol; 
  out[98] += 0.25*win[0]*fin[98]*subvol; 
  out[99] += 0.25*win[0]*fin[99]*subvol; 
  out[100] += 0.25*win[0]*fin[100]*subvol; 
  out[106] += 0.25*win[0]*fin[106]*subvol; 
  out[109] += 0.25*win[0]*fin[109]*subvol; 
  out[114] += 0.25*win[0]*fin[114]*subvol; 
  out[122] += 0.25*win[0]*fin[122]*subvol; 
  out[123] += 0.25*win[0]*fin[123]*subvol; 
  out[124] += 0.25*win[0]*fin[124]*subvol; 
  out[134] += 0.25*win[0]*fin[134]*subvol; 
  out[135] += 0.25*win[0]*fin[135]*subvol; 
  out[136] += 0.25*win[0]*fin[136]*subvol; 
  out[137] += 0.25*win[0]*fin[137]*subvol; 
  out[148] += 0.25*win[0]*fin[148]*subvol; 
  out[149] += 0.25*win[0]*fin[149]*subvol; 
  out[150] += 0.25*win[0]*fin[150]*subvol; 
  out[159] += 0.25*win[0]*fin[159]*subvol; 
  out[167] += 0.25*win[0]*fin[167]*subvol; 
  out[178] += 0.25*win[0]*fin[178]*subvol; 
  out[179] += 0.25*win[0]*fin[179]*subvol; 
  out[180] += 0.25*win[0]*fin[180]*subvol; 
  out[198] += 0.25*win[0]*fin[198]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.25*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.25*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.25*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.25*win[0]*fin[4]*subvol); 
  atomicAdd(&out[4], 0.25*win[0]*fin[6]*subvol); 
  atomicAdd(&out[5], 0.25*win[0]*fin[9]*subvol); 
  atomicAdd(&out[6], 0.25*win[0]*fin[11]*subvol); 
  atomicAdd(&out[7], 0.25*win[0]*fin[12]*subvol); 
  atomicAdd(&out[8], 0.25*win[0]*fin[18]*subvol); 
  atomicAdd(&out[9], 0.25*win[0]*fin[19]*subvol); 
  atomicAdd(&out[10], 0.25*win[0]*fin[20]*subvol); 
  atomicAdd(&out[11], 0.25*win[0]*fin[22]*subvol); 
  atomicAdd(&out[12], 0.25*win[0]*fin[24]*subvol); 
  atomicAdd(&out[13], 0.25*win[0]*fin[28]*subvol); 
  atomicAdd(&out[14], 0.25*win[0]*fin[37]*subvol); 
  atomicAdd(&out[15], 0.25*win[0]*fin[39]*subvol); 
  atomicAdd(&out[16], 0.25*win[0]*fin[40]*subvol); 
  atomicAdd(&out[17], 0.25*win[0]*fin[46]*subvol); 
  atomicAdd(&out[18], 0.25*win[0]*fin[47]*subvol); 
  atomicAdd(&out[19], 0.25*win[0]*fin[53]*subvol); 
  atomicAdd(&out[20], 0.25*win[0]*fin[56]*subvol); 
  atomicAdd(&out[21], 0.25*win[0]*fin[57]*subvol); 
  atomicAdd(&out[22], 0.25*win[0]*fin[58]*subvol); 
  atomicAdd(&out[23], 0.25*win[0]*fin[68]*subvol); 
  atomicAdd(&out[24], 0.25*win[0]*fin[77]*subvol); 
  atomicAdd(&out[25], 0.25*win[0]*fin[88]*subvol); 
  atomicAdd(&out[26], 0.25*win[0]*fin[89]*subvol); 
  atomicAdd(&out[27], 0.25*win[0]*fin[97]*subvol); 
  atomicAdd(&out[28], 0.25*win[0]*fin[99]*subvol); 
  atomicAdd(&out[29], 0.25*win[0]*fin[100]*subvol); 
  atomicAdd(&out[30], 0.25*win[0]*fin[106]*subvol); 
  atomicAdd(&out[31], 0.25*win[0]*fin[124]*subvol); 
  atomicAdd(&out[32], 0.25*win[0]*fin[137]*subvol); 
  atomicAdd(&out[33], 0.25*win[0]*fin[149]*subvol); 
  atomicAdd(&out[34], 0.25*win[0]*fin[150]*subvol); 
  atomicAdd(&out[35], 0.25*win[0]*fin[180]*subvol); 
#else
  out[0] += 0.25*fin[0]*win[0]*subvol; 
  out[1] += 0.25*win[0]*fin[2]*subvol; 
  out[2] += 0.25*win[0]*fin[3]*subvol; 
  out[3] += 0.25*win[0]*fin[4]*subvol; 
  out[4] += 0.25*win[0]*fin[6]*subvol; 
  out[5] += 0.25*win[0]*fin[9]*subvol; 
  out[6] += 0.25*win[0]*fin[11]*subvol; 
  out[7] += 0.25*win[0]*fin[12]*subvol; 
  out[8] += 0.25*win[0]*fin[18]*subvol; 
  out[9] += 0.25*win[0]*fin[19]*subvol; 
  out[10] += 0.25*win[0]*fin[20]*subvol; 
  out[11] += 0.25*win[0]*fin[22]*subvol; 
  out[12] += 0.25*win[0]*fin[24]*subvol; 
  out[13] += 0.25*win[0]*fin[28]*subvol; 
  out[14] += 0.25*win[0]*fin[37]*subvol; 
  out[15] += 0.25*win[0]*fin[39]*subvol; 
  out[16] += 0.25*win[0]*fin[40]*subvol; 
  out[17] += 0.25*win[0]*fin[46]*subvol; 
  out[18] += 0.25*win[0]*fin[47]*subvol; 
  out[19] += 0.25*win[0]*fin[53]*subvol; 
  out[20] += 0.25*win[0]*fin[56]*subvol; 
  out[21] += 0.25*win[0]*fin[57]*subvol; 
  out[22] += 0.25*win[0]*fin[58]*subvol; 
  out[23] += 0.25*win[0]*fin[68]*subvol; 
  out[24] += 0.25*win[0]*fin[77]*subvol; 
  out[25] += 0.25*win[0]*fin[88]*subvol; 
  out[26] += 0.25*win[0]*fin[89]*subvol; 
  out[27] += 0.25*win[0]*fin[97]*subvol; 
  out[28] += 0.25*win[0]*fin[99]*subvol; 
  out[29] += 0.25*win[0]*fin[100]*subvol; 
  out[30] += 0.25*win[0]*fin[106]*subvol; 
  out[31] += 0.25*win[0]*fin[124]*subvol; 
  out[32] += 0.25*win[0]*fin[137]*subvol; 
  out[33] += 0.25*win[0]*fin[149]*subvol; 
  out[34] += 0.25*win[0]*fin[150]*subvol; 
  out[35] += 0.25*win[0]*fin[180]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[6], 0.5*win[0]*fin[6]*subvol); 
  atomicAdd(&out[9], 0.5*win[0]*fin[9]*subvol); 
  atomicAdd(&out[11], 0.5*win[0]*fin[11]*subvol); 
  atomicAdd(&out[12], 0.5*win[0]*fin[12]*subvol); 
  atomicAdd(&out[18], 0.5*win[0]*fin[18]*subvol); 
  atomicAdd(&out[19], 0.5*win[0]*fin[19]*subvol); 
  atomicAdd(&out[20], 0.5*win[0]*fin[20]*subvol); 
  atomicAdd(&out[22], 0.5*win[0]*fin[22]*subvol); 
  atomicAdd(&out[24], 0.5*win[0]*fin[24]*subvol); 
  atomicAdd(&out[28], 0.5*win[0]*fin[28]*subvol); 
  atomicAdd(&out[37], 0.5*win[0]*fin[37]*subvol); 
  atomicAdd(&out[39], 0.5*win[0]*fin[39]*subvol); 
  atomicAdd(&out[40], 0.5*win[0]*fin[40]*subvol); 
  atomicAdd(&out[46], 0.5*win[0]*fin[46]*subvol); 
  atomicAdd(&out[47], 0.5*win[0]*fin[47]*subvol); 
  atomicAdd(&out[53], 0.5*win[0]*fin[53]*subvol); 
  atomicAdd(&out[56], 0.5*win[0]*fin[56]*subvol); 
  atomicAdd(&out[57], 0.5*win[0]*fin[57]*subvol); 
  atomicAdd(&out[58], 0.5*win[0]*fin[58]*subvol); 
  atomicAdd(&out[68], 0.5*win[0]*fin[68]*subvol); 
  atomicAdd(&out[77], 0.5*win[0]*fin[77]*subvol); 
  atomicAdd(&out[88], 0.5*win[0]*fin[88]*subvol); 
  atomicAdd(&out[89], 0.5*win[0]*fin[89]*subvol); 
  atomicAdd(&out[97], 0.5*win[0]*fin[97]*subvol); 
  atomicAdd(&out[99], 0.5*win[0]*fin[99]*subvol); 
  atomicAdd(&out[100], 0.5*win[0]*fin[100]*subvol); 
  atomicAdd(&out[106], 0.5*win[0]*fin[106]*subvol); 
  atomicAdd(&out[124], 0.5*win[0]*fin[124]*subvol); 
  atomicAdd(&out[137], 0.5*win[0]*fin[137]*subvol); 
  atomicAdd(&out[149], 0.5*win[0]*fin[149]*subvol); 
  atomicAdd(&out[150], 0.5*win[0]*fin[150]*subvol); 
  atomicAdd(&out[180], 0.5*win[0]*fin[180]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[2] += 0.5*win[0]*fin[2]*subvol; 
  out[3] += 0.5*win[0]*fin[3]*subvol; 
  out[4] += 0.5*win[0]*fin[4]*subvol; 
  out[6] += 0.5*win[0]*fin[6]*subvol; 
  out[9] += 0.5*win[0]*fin[9]*subvol; 
  out[11] += 0.5*win[0]*fin[11]*subvol; 
  out[12] += 0.5*win[0]*fin[12]*subvol; 
  out[18] += 0.5*win[0]*fin[18]*subvol; 
  out[19] += 0.5*win[0]*fin[19]*subvol; 
  out[20] += 0.5*win[0]*fin[20]*subvol; 
  out[22] += 0.5*win[0]*fin[22]*subvol; 
  out[24] += 0.5*win[0]*fin[24]*subvol; 
  out[28] += 0.5*win[0]*fin[28]*subvol; 
  out[37] += 0.5*win[0]*fin[37]*subvol; 
  out[39] += 0.5*win[0]*fin[39]*subvol; 
  out[40] += 0.5*win[0]*fin[40]*subvol; 
  out[46] += 0.5*win[0]*fin[46]*subvol; 
  out[47] += 0.5*win[0]*fin[47]*subvol; 
  out[53] += 0.5*win[0]*fin[53]*subvol; 
  out[56] += 0.5*win[0]*fin[56]*subvol; 
  out[57] += 0.5*win[0]*fin[57]*subvol; 
  out[58] += 0.5*win[0]*fin[58]*subvol; 
  out[68] += 0.5*win[0]*fin[68]*subvol; 
  out[77] += 0.5*win[0]*fin[77]*subvol; 
  out[88] += 0.5*win[0]*fin[88]*subvol; 
  out[89] += 0.5*win[0]*fin[89]*subvol; 
  out[97] += 0.5*win[0]*fin[97]*subvol; 
  out[99] += 0.5*win[0]*fin[99]*subvol; 
  out[100] += 0.5*win[0]*fin[100]*subvol; 
  out[106] += 0.5*win[0]*fin[106]*subvol; 
  out[124] += 0.5*win[0]*fin[124]*subvol; 
  out[137] += 0.5*win[0]*fin[137]*subvol; 
  out[149] += 0.5*win[0]*fin[149]*subvol; 
  out[150] += 0.5*win[0]*fin[150]*subvol; 
  out[180] += 0.5*win[0]*fin[180]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgyvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.25*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.25*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.25*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.25*win[0]*fin[4]*subvol); 
  atomicAdd(&out[4], 0.25*win[0]*fin[6]*subvol); 
  atomicAdd(&out[5], 0.25*win[0]*fin[8]*subvol); 
  atomicAdd(&out[6], 0.25*win[0]*fin[10]*subvol); 
  atomicAdd(&out[7], 0.25*win[0]*fin[12]*subvol); 
  atomicAdd(&out[8], 0.25*win[0]*fin[17]*subvol); 
  atomicAdd(&out[9], 0.25*win[0]*fin[19]*subvol); 
  atomicAdd(&out[10], 0.25*win[0]*fin[20]*subvol); 
  atomicAdd(&out[11], 0.25*win[0]*fin[22]*subvol); 
  atomicAdd(&out[12], 0.25*win[0]*fin[24]*subvol); 
  atomicAdd(&out[13], 0.25*win[0]*fin[27]*subvol); 
  atomicAdd(&out[14], 0.25*win[0]*fin[36]*subvol); 
  atomicAdd(&out[15], 0.25*win[0]*fin[38]*subvol); 
  atomicAdd(&out[16], 0.25*win[0]*fin[40]*subvol); 
  atomicAdd(&out[17], 0.25*win[0]*fin[45]*subvol); 
  atomicAdd(&out[18], 0.25*win[0]*fin[47]*subvol); 
  atomicAdd(&out[19], 0.25*win[0]*fin[53]*subvol); 
  atomicAdd(&out[20], 0.25*win[0]*fin[55]*subvol); 
  atomicAdd(&out[21], 0.25*win[0]*fin[57]*subvol); 
  atomicAdd(&out[22], 0.25*win[0]*fin[58]*subvol); 
  atomicAdd(&out[23], 0.25*win[0]*fin[67]*subvol); 
  atomicAdd(&out[24], 0.25*win[0]*fin[76]*subvol); 
  atomicAdd(&out[25], 0.25*win[0]*fin[87]*subvol); 
  atomicAdd(&out[26], 0.25*win[0]*fin[89]*subvol); 
  atomicAdd(&out[27], 0.25*win[0]*fin[96]*subvol); 
  atomicAdd(&out[28], 0.25*win[0]*fin[98]*subvol); 
  atomicAdd(&out[29], 0.25*win[0]*fin[100]*subvol); 
  atomicAdd(&out[30], 0.25*win[0]*fin[106]*subvol); 
  atomicAdd(&out[31], 0.25*win[0]*fin[123]*subvol); 
  atomicAdd(&out[32], 0.25*win[0]*fin[136]*subvol); 
  atomicAdd(&out[33], 0.25*win[0]*fin[148]*subvol); 
  atomicAdd(&out[34], 0.25*win[0]*fin[150]*subvol); 
  atomicAdd(&out[35], 0.25*win[0]*fin[179]*subvol); 
#else
  out[0] += 0.25*fin[0]*win[0]*subvol; 
  out[1] += 0.25*win[0]*fin[1]*subvol; 
  out[2] += 0.25*win[0]*fin[3]*subvol; 
  out[3] += 0.25*win[0]*fin[4]*subvol; 
  out[4] += 0.25*win[0]*fin[6]*subvol; 
  out[5] += 0.25*win[0]*fin[8]*subvol; 
  out[6] += 0.25*win[0]*fin[10]*subvol; 
  out[7] += 0.25*win[0]*fin[12]*subvol; 
  out[8] += 0.25*win[0]*fin[17]*subvol; 
  out[9] += 0.25*win[0]*fin[19]*subvol; 
  out[10] += 0.25*win[0]*fin[20]*subvol; 
  out[11] += 0.25*win[0]*fin[22]*subvol; 
  out[12] += 0.25*win[0]*fin[24]*subvol; 
  out[13] += 0.25*win[0]*fin[27]*subvol; 
  out[14] += 0.25*win[0]*fin[36]*subvol; 
  out[15] += 0.25*win[0]*fin[38]*subvol; 
  out[16] += 0.25*win[0]*fin[40]*subvol; 
  out[17] += 0.25*win[0]*fin[45]*subvol; 
  out[18] += 0.25*win[0]*fin[47]*subvol; 
  out[19] += 0.25*win[0]*fin[53]*subvol; 
  out[20] += 0.25*win[0]*fin[55]*subvol; 
  out[21] += 0.25*win[0]*fin[57]*subvol; 
  out[22] += 0.25*win[0]*fin[58]*subvol; 
  out[23] += 0.25*win[0]*fin[67]*subvol; 
  out[24] += 0.25*win[0]*fin[76]*subvol; 
  out[25] += 0.25*win[0]*fin[87]*subvol; 
  out[26] += 0.25*win[0]*fin[89]*subvol; 
  out[27] += 0.25*win[0]*fin[96]*subvol; 
  out[28] += 0.25*win[0]*fin[98]*subvol; 
  out[29] += 0.25*win[0]*fin[100]*subvol; 
  out[30] += 0.25*win[0]*fin[106]*subvol; 
  out[31] += 0.25*win[0]*fin[123]*subvol; 
  out[32] += 0.25*win[0]*fin[136]*subvol; 
  out[33] += 0.25*win[0]*fin[148]*subvol; 
  out[34] += 0.25*win[0]*fin[150]*subvol; 
  out[35] += 0.25*win[0]*fin[179]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgyvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[6], 0.5*win[0]*fin[6]*subvol); 
  atomicAdd(&out[8], 0.5*win[0]*fin[8]*subvol); 
  atomicAdd(&out[10], 0.5*win[0]*fin[10]*subvol); 
  atomicAdd(&out[12], 0.5*win[0]*fin[12]*subvol); 
  atomicAdd(&out[17], 0.5*win[0]*fin[17]*subvol); 
  atomicAdd(&out[19], 0.5*win[0]*fin[19]*subvol); 
  atomicAdd(&out[20], 0.5*win[0]*fin[20]*subvol); 
  atomicAdd(&out[22], 0.5*win[0]*fin[22]*subvol); 
  atomicAdd(&out[24], 0.5*win[0]*fin[24]*subvol); 
  atomicAdd(&out[27], 0.5*win[0]*fin[27]*subvol); 
  atomicAdd(&out[36], 0.5*win[0]*fin[36]*subvol); 
  atomicAdd(&out[38], 0.5*win[0]*fin[38]*subvol); 
  atomicAdd(&out[40], 0.5*win[0]*fin[40]*subvol); 
  atomicAdd(&out[45], 0.5*win[0]*fin[45]*subvol); 
  atomicAdd(&out[47], 0.5*win[0]*fin[47]*subvol); 
  atomicAdd(&out[53], 0.5*win[0]*fin[53]*subvol); 
  atomicAdd(&out[55], 0.5*win[0]*fin[55]*subvol); 
  atomicAdd(&out[57], 0.5*win[0]*fin[57]*subvol); 
  atomicAdd(&out[58], 0.5*win[0]*fin[58]*subvol); 
  atomicAdd(&out[67], 0.5*win[0]*fin[67]*subvol); 
  atomicAdd(&out[76], 0.5*win[0]*fin[76]*subvol); 
  atomicAdd(&out[87], 0.5*win[0]*fin[87]*subvol); 
  atomicAdd(&out[89], 0.5*win[0]*fin[89]*subvol); 
  atomicAdd(&out[96], 0.5*win[0]*fin[96]*subvol); 
  atomicAdd(&out[98], 0.5*win[0]*fin[98]*subvol); 
  atomicAdd(&out[100], 0.5*win[0]*fin[100]*subvol); 
  atomicAdd(&out[106], 0.5*win[0]*fin[106]*subvol); 
  atomicAdd(&out[123], 0.5*win[0]*fin[123]*subvol); 
  atomicAdd(&out[136], 0.5*win[0]*fin[136]*subvol); 
  atomicAdd(&out[148], 0.5*win[0]*fin[148]*subvol); 
  atomicAdd(&out[150], 0.5*win[0]*fin[150]*subvol); 
  atomicAdd(&out[179], 0.5*win[0]*fin[179]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[3] += 0.5*win[0]*fin[3]*subvol; 
  out[4] += 0.5*win[0]*fin[4]*subvol; 
  out[6] += 0.5*win[0]*fin[6]*subvol; 
  out[8] += 0.5*win[0]*fin[8]*subvol; 
  out[10] += 0.5*win[0]*fin[10]*subvol; 
  out[12] += 0.5*win[0]*fin[12]*subvol; 
  out[17] += 0.5*win[0]*fin[17]*subvol; 
  out[19] += 0.5*win[0]*fin[19]*subvol; 
  out[20] += 0.5*win[0]*fin[20]*subvol; 
  out[22] += 0.5*win[0]*fin[22]*subvol; 
  out[24] += 0.5*win[0]*fin[24]*subvol; 
  out[27] += 0.5*win[0]*fin[27]*subvol; 
  out[36] += 0.5*win[0]*fin[36]*subvol; 
  out[38] += 0.5*win[0]*fin[38]*subvol; 
  out[40] += 0.5*win[0]*fin[40]*subvol; 
  out[45] += 0.5*win[0]*fin[45]*subvol; 
  out[47] += 0.5*win[0]*fin[47]*subvol; 
  out[53] += 0.5*win[0]*fin[53]*subvol; 
  out[55] += 0.5*win[0]*fin[55]*subvol; 
  out[57] += 0.5*win[0]*fin[57]*subvol; 
  out[58] += 0.5*win[0]*fin[58]*subvol; 
  out[67] += 0.5*win[0]*fin[67]*subvol; 
  out[76] += 0.5*win[0]*fin[76]*subvol; 
  out[87] += 0.5*win[0]*fin[87]*subvol; 
  out[89] += 0.5*win[0]*fin[89]*subvol; 
  out[96] += 0.5*win[0]*fin[96]*subvol; 
  out[98] += 0.5*win[0]*fin[98]*subvol; 
  out[100] += 0.5*win[0]*fin[100]*subvol; 
  out[106] += 0.5*win[0]*fin[106]*subvol; 
  out[123] += 0.5*win[0]*fin[123]*subvol; 
  out[136] += 0.5*win[0]*fin[136]*subvol; 
  out[148] += 0.5*win[0]*fin[148]*subvol; 
  out[150] += 0.5*win[0]*fin[150]*subvol; 
  out[179] += 0.5*win[0]*fin[179]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxyvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[3]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[4]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[6]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[12]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[19]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[20]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[22]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[24]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[40]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[47]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[53]*subvol); 
  atomicAdd(&out[12], 0.3535533905932737*win[0]*fin[57]*subvol); 
  atomicAdd(&out[13], 0.3535533905932737*win[0]*fin[58]*subvol); 
  atomicAdd(&out[14], 0.3535533905932737*win[0]*fin[89]*subvol); 
  atomicAdd(&out[15], 0.3535533905932737*win[0]*fin[100]*subvol); 
  atomicAdd(&out[16], 0.3535533905932737*win[0]*fin[106]*subvol); 
  atomicAdd(&out[17], 0.3535533905932737*win[0]*fin[150]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[3]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[4]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[6]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[12]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[19]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[20]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[22]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[24]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[40]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[47]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[53]*subvol; 
  out[12] += 0.3535533905932737*win[0]*fin[57]*subvol; 
  out[13] += 0.3535533905932737*win[0]*fin[58]*subvol; 
  out[14] += 0.3535533905932737*win[0]*fin[89]*subvol; 
  out[15] += 0.3535533905932737*win[0]*fin[100]*subvol; 
  out[16] += 0.3535533905932737*win[0]*fin[106]*subvol; 
  out[17] += 0.3535533905932737*win[0]*fin[150]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxyvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
  atomicAdd(&out[3], win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], win[0]*fin[4]*subvol); 
  atomicAdd(&out[6], win[0]*fin[6]*subvol); 
  atomicAdd(&out[12], win[0]*fin[12]*subvol); 
  atomicAdd(&out[19], win[0]*fin[19]*subvol); 
  atomicAdd(&out[20], win[0]*fin[20]*subvol); 
  atomicAdd(&out[22], win[0]*fin[22]*subvol); 
  atomicAdd(&out[24], win[0]*fin[24]*subvol); 
  atomicAdd(&out[40], win[0]*fin[40]*subvol); 
  atomicAdd(&out[47], win[0]*fin[47]*subvol); 
  atomicAdd(&out[53], win[0]*fin[53]*subvol); 
  atomicAdd(&out[57], win[0]*fin[57]*subvol); 
  atomicAdd(&out[58], win[0]*fin[58]*subvol); 
  atomicAdd(&out[89], win[0]*fin[89]*subvol); 
  atomicAdd(&out[100], win[0]*fin[100]*subvol); 
  atomicAdd(&out[106], win[0]*fin[106]*subvol); 
  atomicAdd(&out[150], win[0]*fin[150]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
  out[3] += win[0]*fin[3]*subvol; 
  out[4] += win[0]*fin[4]*subvol; 
  out[6] += win[0]*fin[6]*subvol; 
  out[12] += win[0]*fin[12]*subvol; 
  out[19] += win[0]*fin[19]*subvol; 
  out[20] += win[0]*fin[20]*subvol; 
  out[22] += win[0]*fin[22]*subvol; 
  out[24] += win[0]*fin[24]*subvol; 
  out[40] += win[0]*fin[40]*subvol; 
  out[47] += win[0]*fin[47]*subvol; 
  out[53] += win[0]*fin[53]*subvol; 
  out[57] += win[0]*fin[57]*subvol; 
  out[58] += win[0]*fin[58]*subvol; 
  out[89] += win[0]*fin[89]*subvol; 
  out[100] += win[0]*fin[100]*subvol; 
  out[106] += win[0]*fin[106]*subvol; 
  out[150] += win[0]*fin[150]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgzvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.25*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.25*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.25*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.25*win[0]*fin[4]*subvol); 
  atomicAdd(&out[4], 0.25*win[0]*fin[6]*subvol); 
  atomicAdd(&out[5], 0.25*win[0]*fin[7]*subvol); 
  atomicAdd(&out[6], 0.25*win[0]*fin[10]*subvol); 
  atomicAdd(&out[7], 0.25*win[0]*fin[11]*subvol); 
  atomicAdd(&out[8], 0.25*win[0]*fin[17]*subvol); 
  atomicAdd(&out[9], 0.25*win[0]*fin[18]*subvol); 
  atomicAdd(&out[10], 0.25*win[0]*fin[20]*subvol); 
  atomicAdd(&out[11], 0.25*win[0]*fin[22]*subvol); 
  atomicAdd(&out[12], 0.25*win[0]*fin[24]*subvol); 
  atomicAdd(&out[13], 0.25*win[0]*fin[26]*subvol); 
  atomicAdd(&out[14], 0.25*win[0]*fin[35]*subvol); 
  atomicAdd(&out[15], 0.25*win[0]*fin[38]*subvol); 
  atomicAdd(&out[16], 0.25*win[0]*fin[39]*subvol); 
  atomicAdd(&out[17], 0.25*win[0]*fin[45]*subvol); 
  atomicAdd(&out[18], 0.25*win[0]*fin[46]*subvol); 
  atomicAdd(&out[19], 0.25*win[0]*fin[53]*subvol); 
  atomicAdd(&out[20], 0.25*win[0]*fin[55]*subvol); 
  atomicAdd(&out[21], 0.25*win[0]*fin[56]*subvol); 
  atomicAdd(&out[22], 0.25*win[0]*fin[58]*subvol); 
  atomicAdd(&out[23], 0.25*win[0]*fin[66]*subvol); 
  atomicAdd(&out[24], 0.25*win[0]*fin[75]*subvol); 
  atomicAdd(&out[25], 0.25*win[0]*fin[87]*subvol); 
  atomicAdd(&out[26], 0.25*win[0]*fin[88]*subvol); 
  atomicAdd(&out[27], 0.25*win[0]*fin[95]*subvol); 
  atomicAdd(&out[28], 0.25*win[0]*fin[98]*subvol); 
  atomicAdd(&out[29], 0.25*win[0]*fin[99]*subvol); 
  atomicAdd(&out[30], 0.25*win[0]*fin[106]*subvol); 
  atomicAdd(&out[31], 0.25*win[0]*fin[122]*subvol); 
  atomicAdd(&out[32], 0.25*win[0]*fin[135]*subvol); 
  atomicAdd(&out[33], 0.25*win[0]*fin[148]*subvol); 
  atomicAdd(&out[34], 0.25*win[0]*fin[149]*subvol); 
  atomicAdd(&out[35], 0.25*win[0]*fin[178]*subvol); 
#else
  out[0] += 0.25*fin[0]*win[0]*subvol; 
  out[1] += 0.25*win[0]*fin[1]*subvol; 
  out[2] += 0.25*win[0]*fin[2]*subvol; 
  out[3] += 0.25*win[0]*fin[4]*subvol; 
  out[4] += 0.25*win[0]*fin[6]*subvol; 
  out[5] += 0.25*win[0]*fin[7]*subvol; 
  out[6] += 0.25*win[0]*fin[10]*subvol; 
  out[7] += 0.25*win[0]*fin[11]*subvol; 
  out[8] += 0.25*win[0]*fin[17]*subvol; 
  out[9] += 0.25*win[0]*fin[18]*subvol; 
  out[10] += 0.25*win[0]*fin[20]*subvol; 
  out[11] += 0.25*win[0]*fin[22]*subvol; 
  out[12] += 0.25*win[0]*fin[24]*subvol; 
  out[13] += 0.25*win[0]*fin[26]*subvol; 
  out[14] += 0.25*win[0]*fin[35]*subvol; 
  out[15] += 0.25*win[0]*fin[38]*subvol; 
  out[16] += 0.25*win[0]*fin[39]*subvol; 
  out[17] += 0.25*win[0]*fin[45]*subvol; 
  out[18] += 0.25*win[0]*fin[46]*subvol; 
  out[19] += 0.25*win[0]*fin[53]*subvol; 
  out[20] += 0.25*win[0]*fin[55]*subvol; 
  out[21] += 0.25*win[0]*fin[56]*subvol; 
  out[22] += 0.25*win[0]*fin[58]*subvol; 
  out[23] += 0.25*win[0]*fin[66]*subvol; 
  out[24] += 0.25*win[0]*fin[75]*subvol; 
  out[25] += 0.25*win[0]*fin[87]*subvol; 
  out[26] += 0.25*win[0]*fin[88]*subvol; 
  out[27] += 0.25*win[0]*fin[95]*subvol; 
  out[28] += 0.25*win[0]*fin[98]*subvol; 
  out[29] += 0.25*win[0]*fin[99]*subvol; 
  out[30] += 0.25*win[0]*fin[106]*subvol; 
  out[31] += 0.25*win[0]*fin[122]*subvol; 
  out[32] += 0.25*win[0]*fin[135]*subvol; 
  out[33] += 0.25*win[0]*fin[148]*subvol; 
  out[34] += 0.25*win[0]*fin[149]*subvol; 
  out[35] += 0.25*win[0]*fin[178]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgzvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[6], 0.5*win[0]*fin[6]*subvol); 
  atomicAdd(&out[7], 0.5*win[0]*fin[7]*subvol); 
  atomicAdd(&out[10], 0.5*win[0]*fin[10]*subvol); 
  atomicAdd(&out[11], 0.5*win[0]*fin[11]*subvol); 
  atomicAdd(&out[17], 0.5*win[0]*fin[17]*subvol); 
  atomicAdd(&out[18], 0.5*win[0]*fin[18]*subvol); 
  atomicAdd(&out[20], 0.5*win[0]*fin[20]*subvol); 
  atomicAdd(&out[22], 0.5*win[0]*fin[22]*subvol); 
  atomicAdd(&out[24], 0.5*win[0]*fin[24]*subvol); 
  atomicAdd(&out[26], 0.5*win[0]*fin[26]*subvol); 
  atomicAdd(&out[35], 0.5*win[0]*fin[35]*subvol); 
  atomicAdd(&out[38], 0.5*win[0]*fin[38]*subvol); 
  atomicAdd(&out[39], 0.5*win[0]*fin[39]*subvol); 
  atomicAdd(&out[45], 0.5*win[0]*fin[45]*subvol); 
  atomicAdd(&out[46], 0.5*win[0]*fin[46]*subvol); 
  atomicAdd(&out[53], 0.5*win[0]*fin[53]*subvol); 
  atomicAdd(&out[55], 0.5*win[0]*fin[55]*subvol); 
  atomicAdd(&out[56], 0.5*win[0]*fin[56]*subvol); 
  atomicAdd(&out[58], 0.5*win[0]*fin[58]*subvol); 
  atomicAdd(&out[66], 0.5*win[0]*fin[66]*subvol); 
  atomicAdd(&out[75], 0.5*win[0]*fin[75]*subvol); 
  atomicAdd(&out[87], 0.5*win[0]*fin[87]*subvol); 
  atomicAdd(&out[88], 0.5*win[0]*fin[88]*subvol); 
  atomicAdd(&out[95], 0.5*win[0]*fin[95]*subvol); 
  atomicAdd(&out[98], 0.5*win[0]*fin[98]*subvol); 
  atomicAdd(&out[99], 0.5*win[0]*fin[99]*subvol); 
  atomicAdd(&out[106], 0.5*win[0]*fin[106]*subvol); 
  atomicAdd(&out[122], 0.5*win[0]*fin[122]*subvol); 
  atomicAdd(&out[135], 0.5*win[0]*fin[135]*subvol); 
  atomicAdd(&out[148], 0.5*win[0]*fin[148]*subvol); 
  atomicAdd(&out[149], 0.5*win[0]*fin[149]*subvol); 
  atomicAdd(&out[178], 0.5*win[0]*fin[178]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[2] += 0.5*win[0]*fin[2]*subvol; 
  out[4] += 0.5*win[0]*fin[4]*subvol; 
  out[6] += 0.5*win[0]*fin[6]*subvol; 
  out[7] += 0.5*win[0]*fin[7]*subvol; 
  out[10] += 0.5*win[0]*fin[10]*subvol; 
  out[11] += 0.5*win[0]*fin[11]*subvol; 
  out[17] += 0.5*win[0]*fin[17]*subvol; 
  out[18] += 0.5*win[0]*fin[18]*subvol; 
  out[20] += 0.5*win[0]*fin[20]*subvol; 
  out[22] += 0.5*win[0]*fin[22]*subvol; 
  out[24] += 0.5*win[0]*fin[24]*subvol; 
  out[26] += 0.5*win[0]*fin[26]*subvol; 
  out[35] += 0.5*win[0]*fin[35]*subvol; 
  out[38] += 0.5*win[0]*fin[38]*subvol; 
  out[39] += 0.5*win[0]*fin[39]*subvol; 
  out[45] += 0.5*win[0]*fin[45]*subvol; 
  out[46] += 0.5*win[0]*fin[46]*subvol; 
  out[53] += 0.5*win[0]*fin[53]*subvol; 
  out[55] += 0.5*win[0]*fin[55]*subvol; 
  out[56] += 0.5*win[0]*fin[56]*subvol; 
  out[58] += 0.5*win[0]*fin[58]*subvol; 
  out[66] += 0.5*win[0]*fin[66]*subvol; 
  out[75] += 0.5*win[0]*fin[75]*subvol; 
  out[87] += 0.5*win[0]*fin[87]*subvol; 
  out[88] += 0.5*win[0]*fin[88]*subvol; 
  out[95] += 0.5*win[0]*fin[95]*subvol; 
  out[98] += 0.5*win[0]*fin[98]*subvol; 
  out[99] += 0.5*win[0]*fin[99]*subvol; 
  out[106] += 0.5*win[0]*fin[106]*subvol; 
  out[122] += 0.5*win[0]*fin[122]*subvol; 
  out[135] += 0.5*win[0]*fin[135]*subvol; 
  out[148] += 0.5*win[0]*fin[148]*subvol; 
  out[149] += 0.5*win[0]*fin[149]*subvol; 
  out[178] += 0.5*win[0]*fin[178]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxzvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[4]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[6]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[11]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[18]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[20]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[22]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[24]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[39]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[46]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[53]*subvol); 
  atomicAdd(&out[12], 0.3535533905932737*win[0]*fin[56]*subvol); 
  atomicAdd(&out[13], 0.3535533905932737*win[0]*fin[58]*subvol); 
  atomicAdd(&out[14], 0.3535533905932737*win[0]*fin[88]*subvol); 
  atomicAdd(&out[15], 0.3535533905932737*win[0]*fin[99]*subvol); 
  atomicAdd(&out[16], 0.3535533905932737*win[0]*fin[106]*subvol); 
  atomicAdd(&out[17], 0.3535533905932737*win[0]*fin[149]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[2]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[4]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[6]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[11]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[18]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[20]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[22]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[24]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[39]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[46]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[53]*subvol; 
  out[12] += 0.3535533905932737*win[0]*fin[56]*subvol; 
  out[13] += 0.3535533905932737*win[0]*fin[58]*subvol; 
  out[14] += 0.3535533905932737*win[0]*fin[88]*subvol; 
  out[15] += 0.3535533905932737*win[0]*fin[99]*subvol; 
  out[16] += 0.3535533905932737*win[0]*fin[106]*subvol; 
  out[17] += 0.3535533905932737*win[0]*fin[149]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxzvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], win[0]*fin[2]*subvol); 
  atomicAdd(&out[4], win[0]*fin[4]*subvol); 
  atomicAdd(&out[6], win[0]*fin[6]*subvol); 
  atomicAdd(&out[11], win[0]*fin[11]*subvol); 
  atomicAdd(&out[18], win[0]*fin[18]*subvol); 
  atomicAdd(&out[20], win[0]*fin[20]*subvol); 
  atomicAdd(&out[22], win[0]*fin[22]*subvol); 
  atomicAdd(&out[24], win[0]*fin[24]*subvol); 
  atomicAdd(&out[39], win[0]*fin[39]*subvol); 
  atomicAdd(&out[46], win[0]*fin[46]*subvol); 
  atomicAdd(&out[53], win[0]*fin[53]*subvol); 
  atomicAdd(&out[56], win[0]*fin[56]*subvol); 
  atomicAdd(&out[58], win[0]*fin[58]*subvol); 
  atomicAdd(&out[88], win[0]*fin[88]*subvol); 
  atomicAdd(&out[99], win[0]*fin[99]*subvol); 
  atomicAdd(&out[106], win[0]*fin[106]*subvol); 
  atomicAdd(&out[149], win[0]*fin[149]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
  out[2] += win[0]*fin[2]*subvol; 
  out[4] += win[0]*fin[4]*subvol; 
  out[6] += win[0]*fin[6]*subvol; 
  out[11] += win[0]*fin[11]*subvol; 
  out[18] += win[0]*fin[18]*subvol; 
  out[20] += win[0]*fin[20]*subvol; 
  out[22] += win[0]*fin[22]*subvol; 
  out[24] += win[0]*fin[24]*subvol; 
  out[39] += win[0]*fin[39]*subvol; 
  out[46] += win[0]*fin[46]*subvol; 
  out[53] += win[0]*fin[53]*subvol; 
  out[56] += win[0]*fin[56]*subvol; 
  out[58] += win[0]*fin[58]*subvol; 
  out[88] += win[0]*fin[88]*subvol; 
  out[99] += win[0]*fin[99]*subvol; 
  out[106] += win[0]*fin[106]*subvol; 
  out[149] += win[0]*fin[149]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgyzvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[4]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[6]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[10]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[17]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[20]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[22]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[24]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[38]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[45]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[53]*subvol); 
  atomicAdd(&out[12], 0.3535533905932737*win[0]*fin[55]*subvol); 
  atomicAdd(&out[13], 0.3535533905932737*win[0]*fin[58]*subvol); 
  atomicAdd(&out[14], 0.3535533905932737*win[0]*fin[87]*subvol); 
  atomicAdd(&out[15], 0.3535533905932737*win[0]*fin[98]*subvol); 
  atomicAdd(&out[16], 0.3535533905932737*win[0]*fin[106]*subvol); 
  atomicAdd(&out[17], 0.3535533905932737*win[0]*fin[148]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[1]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[4]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[6]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[10]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[17]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[20]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[22]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[24]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[38]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[45]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[53]*subvol; 
  out[12] += 0.3535533905932737*win[0]*fin[55]*subvol; 
  out[13] += 0.3535533905932737*win[0]*fin[58]*subvol; 
  out[14] += 0.3535533905932737*win[0]*fin[87]*subvol; 
  out[15] += 0.3535533905932737*win[0]*fin[98]*subvol; 
  out[16] += 0.3535533905932737*win[0]*fin[106]*subvol; 
  out[17] += 0.3535533905932737*win[0]*fin[148]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgyzvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], win[0]*fin[1]*subvol); 
  atomicAdd(&out[4], win[0]*fin[4]*subvol); 
  atomicAdd(&out[6], win[0]*fin[6]*subvol); 
  atomicAdd(&out[10], win[0]*fin[10]*subvol); 
  atomicAdd(&out[17], win[0]*fin[17]*subvol); 
  atomicAdd(&out[20], win[0]*fin[20]*subvol); 
  atomicAdd(&out[22], win[0]*fin[22]*subvol); 
  atomicAdd(&out[24], win[0]*fin[24]*subvol); 
  atomicAdd(&out[38], win[0]*fin[38]*subvol); 
  atomicAdd(&out[45], win[0]*fin[45]*subvol); 
  atomicAdd(&out[53], win[0]*fin[53]*subvol); 
  atomicAdd(&out[55], win[0]*fin[55]*subvol); 
  atomicAdd(&out[58], win[0]*fin[58]*subvol); 
  atomicAdd(&out[87], win[0]*fin[87]*subvol); 
  atomicAdd(&out[98], win[0]*fin[98]*subvol); 
  atomicAdd(&out[106], win[0]*fin[106]*subvol); 
  atomicAdd(&out[148], win[0]*fin[148]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
  out[1] += win[0]*fin[1]*subvol; 
  out[4] += win[0]*fin[4]*subvol; 
  out[6] += win[0]*fin[6]*subvol; 
  out[10] += win[0]*fin[10]*subvol; 
  out[17] += win[0]*fin[17]*subvol; 
  out[20] += win[0]*fin[20]*subvol; 
  out[22] += win[0]*fin[22]*subvol; 
  out[24] += win[0]*fin[24]*subvol; 
  out[38] += win[0]*fin[38]*subvol; 
  out[45] += win[0]*fin[45]*subvol; 
  out[53] += win[0]*fin[53]*subvol; 
  out[55] += win[0]*fin[55]*subvol; 
  out[58] += win[0]*fin[58]*subvol; 
  out[87] += win[0]*fin[87]*subvol; 
  out[98] += win[0]*fin[98]*subvol; 
  out[106] += win[0]*fin[106]*subvol; 
  out[148] += win[0]*fin[148]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxyzvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[6]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[20]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[22]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[24]*subvol); 
  atomicAdd(&out[6], 0.5*win[0]*fin[53]*subvol); 
  atomicAdd(&out[7], 0.5*win[0]*fin[58]*subvol); 
  atomicAdd(&out[8], 0.5*win[0]*fin[106]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[4]*subvol; 
  out[2] += 0.5*win[0]*fin[6]*subvol; 
  out[3] += 0.5*win[0]*fin[20]*subvol; 
  out[4] += 0.5*win[0]*fin[22]*subvol; 
  out[5] += 0.5*win[0]*fin[24]*subvol; 
  out[6] += 0.5*win[0]*fin[53]*subvol; 
  out[7] += 0.5*win[0]*fin[58]*subvol; 
  out[8] += 0.5*win[0]*fin[106]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxyzvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.0*fin[0]*win[0]*subvol); 
  atomicAdd(&out[4], 2.0*win[0]*fin[4]*subvol); 
  atomicAdd(&out[6], 2.0*win[0]*fin[6]*subvol); 
  atomicAdd(&out[20], 2.0*win[0]*fin[20]*subvol); 
  atomicAdd(&out[22], 2.0*win[0]*fin[22]*subvol); 
  atomicAdd(&out[24], 2.0*win[0]*fin[24]*subvol); 
  atomicAdd(&out[53], 2.0*win[0]*fin[53]*subvol); 
  atomicAdd(&out[58], 2.0*win[0]*fin[58]*subvol); 
  atomicAdd(&out[106], 2.0*win[0]*fin[106]*subvol); 
#else
  out[0] += 2.0*fin[0]*win[0]*subvol; 
  out[4] += 2.0*win[0]*fin[4]*subvol; 
  out[6] += 2.0*win[0]*fin[6]*subvol; 
  out[20] += 2.0*win[0]*fin[20]*subvol; 
  out[22] += 2.0*win[0]*fin[22]*subvol; 
  out[24] += 2.0*win[0]*fin[24]*subvol; 
  out[53] += 2.0*win[0]*fin[53]*subvol; 
  out[58] += 2.0*win[0]*fin[58]*subvol; 
  out[106] += 2.0*win[0]*fin[106]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgvxvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.25*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.25*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.25*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.25*win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], 0.25*win[0]*fin[6]*subvol); 
  atomicAdd(&out[5], 0.25*win[0]*fin[7]*subvol); 
  atomicAdd(&out[6], 0.25*win[0]*fin[8]*subvol); 
  atomicAdd(&out[7], 0.25*win[0]*fin[9]*subvol); 
  atomicAdd(&out[8], 0.25*win[0]*fin[17]*subvol); 
  atomicAdd(&out[9], 0.25*win[0]*fin[18]*subvol); 
  atomicAdd(&out[10], 0.25*win[0]*fin[19]*subvol); 
  atomicAdd(&out[11], 0.25*win[0]*fin[24]*subvol); 
  atomicAdd(&out[12], 0.25*win[0]*fin[25]*subvol); 
  atomicAdd(&out[13], 0.25*win[0]*fin[35]*subvol); 
  atomicAdd(&out[14], 0.25*win[0]*fin[36]*subvol); 
  atomicAdd(&out[15], 0.25*win[0]*fin[37]*subvol); 
  atomicAdd(&out[16], 0.25*win[0]*fin[55]*subvol); 
  atomicAdd(&out[17], 0.25*win[0]*fin[56]*subvol); 
  atomicAdd(&out[18], 0.25*win[0]*fin[57]*subvol); 
  atomicAdd(&out[19], 0.25*win[0]*fin[65]*subvol); 
  atomicAdd(&out[20], 0.25*win[0]*fin[95]*subvol); 
  atomicAdd(&out[21], 0.25*win[0]*fin[96]*subvol); 
  atomicAdd(&out[22], 0.25*win[0]*fin[97]*subvol); 
  atomicAdd(&out[23], 0.25*win[0]*fin[134]*subvol); 
#else
  out[0] += 0.25*fin[0]*win[0]*subvol; 
  out[1] += 0.25*win[0]*fin[1]*subvol; 
  out[2] += 0.25*win[0]*fin[2]*subvol; 
  out[3] += 0.25*win[0]*fin[3]*subvol; 
  out[4] += 0.25*win[0]*fin[6]*subvol; 
  out[5] += 0.25*win[0]*fin[7]*subvol; 
  out[6] += 0.25*win[0]*fin[8]*subvol; 
  out[7] += 0.25*win[0]*fin[9]*subvol; 
  out[8] += 0.25*win[0]*fin[17]*subvol; 
  out[9] += 0.25*win[0]*fin[18]*subvol; 
  out[10] += 0.25*win[0]*fin[19]*subvol; 
  out[11] += 0.25*win[0]*fin[24]*subvol; 
  out[12] += 0.25*win[0]*fin[25]*subvol; 
  out[13] += 0.25*win[0]*fin[35]*subvol; 
  out[14] += 0.25*win[0]*fin[36]*subvol; 
  out[15] += 0.25*win[0]*fin[37]*subvol; 
  out[16] += 0.25*win[0]*fin[55]*subvol; 
  out[17] += 0.25*win[0]*fin[56]*subvol; 
  out[18] += 0.25*win[0]*fin[57]*subvol; 
  out[19] += 0.25*win[0]*fin[65]*subvol; 
  out[20] += 0.25*win[0]*fin[95]*subvol; 
  out[21] += 0.25*win[0]*fin[96]*subvol; 
  out[22] += 0.25*win[0]*fin[97]*subvol; 
  out[23] += 0.25*win[0]*fin[134]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgvxvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[6], 0.5*win[0]*fin[6]*subvol); 
  atomicAdd(&out[7], 0.5*win[0]*fin[7]*subvol); 
  atomicAdd(&out[8], 0.5*win[0]*fin[8]*subvol); 
  atomicAdd(&out[9], 0.5*win[0]*fin[9]*subvol); 
  atomicAdd(&out[17], 0.5*win[0]*fin[17]*subvol); 
  atomicAdd(&out[18], 0.5*win[0]*fin[18]*subvol); 
  atomicAdd(&out[19], 0.5*win[0]*fin[19]*subvol); 
  atomicAdd(&out[24], 0.5*win[0]*fin[24]*subvol); 
  atomicAdd(&out[25], 0.5*win[0]*fin[25]*subvol); 
  atomicAdd(&out[35], 0.5*win[0]*fin[35]*subvol); 
  atomicAdd(&out[36], 0.5*win[0]*fin[36]*subvol); 
  atomicAdd(&out[37], 0.5*win[0]*fin[37]*subvol); 
  atomicAdd(&out[55], 0.5*win[0]*fin[55]*subvol); 
  atomicAdd(&out[56], 0.5*win[0]*fin[56]*subvol); 
  atomicAdd(&out[57], 0.5*win[0]*fin[57]*subvol); 
  atomicAdd(&out[65], 0.5*win[0]*fin[65]*subvol); 
  atomicAdd(&out[95], 0.5*win[0]*fin[95]*subvol); 
  atomicAdd(&out[96], 0.5*win[0]*fin[96]*subvol); 
  atomicAdd(&out[97], 0.5*win[0]*fin[97]*subvol); 
  atomicAdd(&out[134], 0.5*win[0]*fin[134]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[2] += 0.5*win[0]*fin[2]*subvol; 
  out[3] += 0.5*win[0]*fin[3]*subvol; 
  out[6] += 0.5*win[0]*fin[6]*subvol; 
  out[7] += 0.5*win[0]*fin[7]*subvol; 
  out[8] += 0.5*win[0]*fin[8]*subvol; 
  out[9] += 0.5*win[0]*fin[9]*subvol; 
  out[17] += 0.5*win[0]*fin[17]*subvol; 
  out[18] += 0.5*win[0]*fin[18]*subvol; 
  out[19] += 0.5*win[0]*fin[19]*subvol; 
  out[24] += 0.5*win[0]*fin[24]*subvol; 
  out[25] += 0.5*win[0]*fin[25]*subvol; 
  out[35] += 0.5*win[0]*fin[35]*subvol; 
  out[36] += 0.5*win[0]*fin[36]*subvol; 
  out[37] += 0.5*win[0]*fin[37]*subvol; 
  out[55] += 0.5*win[0]*fin[55]*subvol; 
  out[56] += 0.5*win[0]*fin[56]*subvol; 
  out[57] += 0.5*win[0]*fin[57]*subvol; 
  out[65] += 0.5*win[0]*fin[65]*subvol; 
  out[95] += 0.5*win[0]*fin[95]*subvol; 
  out[96] += 0.5*win[0]*fin[96]*subvol; 
  out[97] += 0.5*win[0]*fin[97]*subvol; 
  out[134] += 0.5*win[0]*fin[134]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxvxvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[6]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[9]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[18]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[19]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[24]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[37]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[56]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[57]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[97]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[2]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[3]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[6]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[9]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[18]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[19]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[24]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[37]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[56]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[57]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[97]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxvxvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], win[0]*fin[3]*subvol); 
  atomicAdd(&out[6], win[0]*fin[6]*subvol); 
  atomicAdd(&out[9], win[0]*fin[9]*subvol); 
  atomicAdd(&out[18], win[0]*fin[18]*subvol); 
  atomicAdd(&out[19], win[0]*fin[19]*subvol); 
  atomicAdd(&out[24], win[0]*fin[24]*subvol); 
  atomicAdd(&out[37], win[0]*fin[37]*subvol); 
  atomicAdd(&out[56], win[0]*fin[56]*subvol); 
  atomicAdd(&out[57], win[0]*fin[57]*subvol); 
  atomicAdd(&out[97], win[0]*fin[97]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
  out[2] += win[0]*fin[2]*subvol; 
  out[3] += win[0]*fin[3]*subvol; 
  out[6] += win[0]*fin[6]*subvol; 
  out[9] += win[0]*fin[9]*subvol; 
  out[18] += win[0]*fin[18]*subvol; 
  out[19] += win[0]*fin[19]*subvol; 
  out[24] += win[0]*fin[24]*subvol; 
  out[37] += win[0]*fin[37]*subvol; 
  out[56] += win[0]*fin[56]*subvol; 
  out[57] += win[0]*fin[57]*subvol; 
  out[97] += win[0]*fin[97]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgyvxvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[6]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[8]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[17]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[19]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[24]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[36]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[55]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[57]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[96]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[1]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[3]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[6]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[8]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[17]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[19]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[24]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[36]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[55]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[57]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[96]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgyvxvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], win[0]*fin[1]*subvol); 
  atomicAdd(&out[3], win[0]*fin[3]*subvol); 
  atomicAdd(&out[6], win[0]*fin[6]*subvol); 
  atomicAdd(&out[8], win[0]*fin[8]*subvol); 
  atomicAdd(&out[17], win[0]*fin[17]*subvol); 
  atomicAdd(&out[19], win[0]*fin[19]*subvol); 
  atomicAdd(&out[24], win[0]*fin[24]*subvol); 
  atomicAdd(&out[36], win[0]*fin[36]*subvol); 
  atomicAdd(&out[55], win[0]*fin[55]*subvol); 
  atomicAdd(&out[57], win[0]*fin[57]*subvol); 
  atomicAdd(&out[96], win[0]*fin[96]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
  out[1] += win[0]*fin[1]*subvol; 
  out[3] += win[0]*fin[3]*subvol; 
  out[6] += win[0]*fin[6]*subvol; 
  out[8] += win[0]*fin[8]*subvol; 
  out[17] += win[0]*fin[17]*subvol; 
  out[19] += win[0]*fin[19]*subvol; 
  out[24] += win[0]*fin[24]*subvol; 
  out[36] += win[0]*fin[36]*subvol; 
  out[55] += win[0]*fin[55]*subvol; 
  out[57] += win[0]*fin[57]*subvol; 
  out[96] += win[0]*fin[96]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxyvxvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[6]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[19]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[24]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[57]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[3]*subvol; 
  out[2] += 0.5*win[0]*fin[6]*subvol; 
  out[3] += 0.5*win[0]*fin[19]*subvol; 
  out[4] += 0.5*win[0]*fin[24]*subvol; 
  out[5] += 0.5*win[0]*fin[57]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxyvxvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.0*fin[0]*win[0]*subvol); 
  atomicAdd(&out[3], 2.0*win[0]*fin[3]*subvol); 
  atomicAdd(&out[6], 2.0*win[0]*fin[6]*subvol); 
  atomicAdd(&out[19], 2.0*win[0]*fin[19]*subvol); 
  atomicAdd(&out[24], 2.0*win[0]*fin[24]*subvol); 
  atomicAdd(&out[57], 2.0*win[0]*fin[57]*subvol); 
#else
  out[0] += 2.0*fin[0]*win[0]*subvol; 
  out[3] += 2.0*win[0]*fin[3]*subvol; 
  out[6] += 2.0*win[0]*fin[6]*subvol; 
  out[19] += 2.0*win[0]*fin[19]*subvol; 
  out[24] += 2.0*win[0]*fin[24]*subvol; 
  out[57] += 2.0*win[0]*fin[57]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgzvxvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[6]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[7]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[17]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[18]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[24]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[35]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[55]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[56]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[95]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[1]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[2]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[6]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[7]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[17]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[18]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[24]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[35]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[55]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[56]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[95]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgzvxvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], win[0]*fin[2]*subvol); 
  atomicAdd(&out[6], win[0]*fin[6]*subvol); 
  atomicAdd(&out[7], win[0]*fin[7]*subvol); 
  atomicAdd(&out[17], win[0]*fin[17]*subvol); 
  atomicAdd(&out[18], win[0]*fin[18]*subvol); 
  atomicAdd(&out[24], win[0]*fin[24]*subvol); 
  atomicAdd(&out[35], win[0]*fin[35]*subvol); 
  atomicAdd(&out[55], win[0]*fin[55]*subvol); 
  atomicAdd(&out[56], win[0]*fin[56]*subvol); 
  atomicAdd(&out[95], win[0]*fin[95]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
  out[1] += win[0]*fin[1]*subvol; 
  out[2] += win[0]*fin[2]*subvol; 
  out[6] += win[0]*fin[6]*subvol; 
  out[7] += win[0]*fin[7]*subvol; 
  out[17] += win[0]*fin[17]*subvol; 
  out[18] += win[0]*fin[18]*subvol; 
  out[24] += win[0]*fin[24]*subvol; 
  out[35] += win[0]*fin[35]*subvol; 
  out[55] += win[0]*fin[55]*subvol; 
  out[56] += win[0]*fin[56]*subvol; 
  out[95] += win[0]*fin[95]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxzvxvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[6]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[18]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[24]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[56]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[2]*subvol; 
  out[2] += 0.5*win[0]*fin[6]*subvol; 
  out[3] += 0.5*win[0]*fin[18]*subvol; 
  out[4] += 0.5*win[0]*fin[24]*subvol; 
  out[5] += 0.5*win[0]*fin[56]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxzvxvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.0*fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], 2.0*win[0]*fin[2]*subvol); 
  atomicAdd(&out[6], 2.0*win[0]*fin[6]*subvol); 
  atomicAdd(&out[18], 2.0*win[0]*fin[18]*subvol); 
  atomicAdd(&out[24], 2.0*win[0]*fin[24]*subvol); 
  atomicAdd(&out[56], 2.0*win[0]*fin[56]*subvol); 
#else
  out[0] += 2.0*fin[0]*win[0]*subvol; 
  out[2] += 2.0*win[0]*fin[2]*subvol; 
  out[6] += 2.0*win[0]*fin[6]*subvol; 
  out[18] += 2.0*win[0]*fin[18]*subvol; 
  out[24] += 2.0*win[0]*fin[24]*subvol; 
  out[56] += 2.0*win[0]*fin[56]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgyzvxvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[6]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[17]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[24]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[55]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[2] += 0.5*win[0]*fin[6]*subvol; 
  out[3] += 0.5*win[0]*fin[17]*subvol; 
  out[4] += 0.5*win[0]*fin[24]*subvol; 
  out[5] += 0.5*win[0]*fin[55]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgyzvxvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.0*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 2.0*win[0]*fin[1]*subvol); 
  atomicAdd(&out[6], 2.0*win[0]*fin[6]*subvol); 
  atomicAdd(&out[17], 2.0*win[0]*fin[17]*subvol); 
  atomicAdd(&out[24], 2.0*win[0]*fin[24]*subvol); 
  atomicAdd(&out[55], 2.0*win[0]*fin[55]*subvol); 
#else
  out[0] += 2.0*fin[0]*win[0]*subvol; 
  out[1] += 2.0*win[0]*fin[1]*subvol; 
  out[6] += 2.0*win[0]*fin[6]*subvol; 
  out[17] += 2.0*win[0]*fin[17]*subvol; 
  out[24] += 2.0*win[0]*fin[24]*subvol; 
  out[55] += 2.0*win[0]*fin[55]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxyzvxvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[6]*subvol); 
  atomicAdd(&out[2], 0.7071067811865475*win[0]*fin[24]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[6]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[24]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxyzvxvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 4.0*fin[0]*win[0]*subvol); 
  atomicAdd(&out[6], 4.0*win[0]*fin[6]*subvol); 
  atomicAdd(&out[24], 4.0*win[0]*fin[24]*subvol); 
#else
  out[0] += 4.0*fin[0]*win[0]*subvol; 
  out[6] += 4.0*win[0]*fin[6]*subvol; 
  out[24] += 4.0*win[0]*fin[24]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.1767766952966368*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.1767766952966368*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.1767766952966368*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.1767766952966368*win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], 0.1767766952966368*win[0]*fin[4]*subvol); 
  atomicAdd(&out[5], 0.1767766952966368*win[0]*fin[5]*subvol); 
  atomicAdd(&out[6], 0.1767766952966368*win[0]*fin[7]*subvol); 
  atomicAdd(&out[7], 0.1767766952966368*win[0]*fin[8]*subvol); 
  atomicAdd(&out[8], 0.1767766952966368*win[0]*fin[9]*subvol); 
  atomicAdd(&out[9], 0.1767766952966368*win[0]*fin[10]*subvol); 
  atomicAdd(&out[10], 0.1767766952966368*win[0]*fin[11]*subvol); 
  atomicAdd(&out[11], 0.1767766952966368*win[0]*fin[12]*subvol); 
  atomicAdd(&out[12], 0.1767766952966368*win[0]*fin[13]*subvol); 
  atomicAdd(&out[13], 0.1767766952966368*win[0]*fin[14]*subvol); 
  atomicAdd(&out[14], 0.1767766952966368*win[0]*fin[15]*subvol); 
  atomicAdd(&out[15], 0.1767766952966368*win[0]*fin[16]*subvol); 
  atomicAdd(&out[16], 0.1767766952966368*win[0]*fin[22]*subvol); 
  atomicAdd(&out[17], 0.1767766952966368*win[0]*fin[23]*subvol); 
  atomicAdd(&out[18], 0.1767766952966368*win[0]*fin[25]*subvol); 
  atomicAdd(&out[19], 0.1767766952966368*win[0]*fin[26]*subvol); 
  atomicAdd(&out[20], 0.1767766952966368*win[0]*fin[27]*subvol); 
  atomicAdd(&out[21], 0.1767766952966368*win[0]*fin[28]*subvol); 
  atomicAdd(&out[22], 0.1767766952966368*win[0]*fin[29]*subvol); 
  atomicAdd(&out[23], 0.1767766952966368*win[0]*fin[30]*subvol); 
  atomicAdd(&out[24], 0.1767766952966368*win[0]*fin[31]*subvol); 
  atomicAdd(&out[25], 0.1767766952966368*win[0]*fin[32]*subvol); 
  atomicAdd(&out[26], 0.1767766952966368*win[0]*fin[33]*subvol); 
  atomicAdd(&out[27], 0.1767766952966368*win[0]*fin[34]*subvol); 
  atomicAdd(&out[28], 0.1767766952966368*win[0]*fin[45]*subvol); 
  atomicAdd(&out[29], 0.1767766952966368*win[0]*fin[46]*subvol); 
  atomicAdd(&out[30], 0.1767766952966368*win[0]*fin[47]*subvol); 
  atomicAdd(&out[31], 0.1767766952966368*win[0]*fin[48]*subvol); 
  atomicAdd(&out[32], 0.1767766952966368*win[0]*fin[49]*subvol); 
  atomicAdd(&out[33], 0.1767766952966368*win[0]*fin[50]*subvol); 
  atomicAdd(&out[34], 0.1767766952966368*win[0]*fin[51]*subvol); 
  atomicAdd(&out[35], 0.1767766952966368*win[0]*fin[52]*subvol); 
  atomicAdd(&out[36], 0.1767766952966368*win[0]*fin[60]*subvol); 
  atomicAdd(&out[37], 0.1767766952966368*win[0]*fin[61]*subvol); 
  atomicAdd(&out[38], 0.1767766952966368*win[0]*fin[62]*subvol); 
  atomicAdd(&out[39], 0.1767766952966368*win[0]*fin[63]*subvol); 
  atomicAdd(&out[40], 0.1767766952966368*win[0]*fin[64]*subvol); 
  atomicAdd(&out[41], 0.1767766952966368*win[0]*fin[75]*subvol); 
  atomicAdd(&out[42], 0.1767766952966368*win[0]*fin[76]*subvol); 
  atomicAdd(&out[43], 0.1767766952966368*win[0]*fin[77]*subvol); 
  atomicAdd(&out[44], 0.1767766952966368*win[0]*fin[78]*subvol); 
  atomicAdd(&out[45], 0.1767766952966368*win[0]*fin[79]*subvol); 
  atomicAdd(&out[46], 0.1767766952966368*win[0]*fin[80]*subvol); 
  atomicAdd(&out[47], 0.1767766952966368*win[0]*fin[81]*subvol); 
  atomicAdd(&out[48], 0.1767766952966368*win[0]*fin[82]*subvol); 
  atomicAdd(&out[49], 0.1767766952966368*win[0]*fin[83]*subvol); 
  atomicAdd(&out[50], 0.1767766952966368*win[0]*fin[84]*subvol); 
  atomicAdd(&out[51], 0.1767766952966368*win[0]*fin[85]*subvol); 
  atomicAdd(&out[52], 0.1767766952966368*win[0]*fin[86]*subvol); 
  atomicAdd(&out[53], 0.1767766952966368*win[0]*fin[105]*subvol); 
  atomicAdd(&out[54], 0.1767766952966368*win[0]*fin[108]*subvol); 
  atomicAdd(&out[55], 0.1767766952966368*win[0]*fin[114]*subvol); 
  atomicAdd(&out[56], 0.1767766952966368*win[0]*fin[115]*subvol); 
  atomicAdd(&out[57], 0.1767766952966368*win[0]*fin[116]*subvol); 
  atomicAdd(&out[58], 0.1767766952966368*win[0]*fin[117]*subvol); 
  atomicAdd(&out[59], 0.1767766952966368*win[0]*fin[118]*subvol); 
  atomicAdd(&out[60], 0.1767766952966368*win[0]*fin[119]*subvol); 
  atomicAdd(&out[61], 0.1767766952966368*win[0]*fin[120]*subvol); 
  atomicAdd(&out[62], 0.1767766952966368*win[0]*fin[121]*subvol); 
  atomicAdd(&out[63], 0.1767766952966368*win[0]*fin[144]*subvol); 
  atomicAdd(&out[64], 0.1767766952966368*win[0]*fin[145]*subvol); 
  atomicAdd(&out[65], 0.1767766952966368*win[0]*fin[146]*subvol); 
  atomicAdd(&out[66], 0.1767766952966368*win[0]*fin[157]*subvol); 
  atomicAdd(&out[67], 0.1767766952966368*win[0]*fin[158]*subvol); 
  atomicAdd(&out[68], 0.1767766952966368*win[0]*fin[172]*subvol); 
  atomicAdd(&out[69], 0.1767766952966368*win[0]*fin[173]*subvol); 
  atomicAdd(&out[70], 0.1767766952966368*win[0]*fin[174]*subvol); 
  atomicAdd(&out[71], 0.1767766952966368*win[0]*fin[194]*subvol); 
#else
  out[0] += 0.1767766952966368*fin[0]*win[0]*subvol; 
  out[1] += 0.1767766952966368*win[0]*fin[1]*subvol; 
  out[2] += 0.1767766952966368*win[0]*fin[2]*subvol; 
  out[3] += 0.1767766952966368*win[0]*fin[3]*subvol; 
  out[4] += 0.1767766952966368*win[0]*fin[4]*subvol; 
  out[5] += 0.1767766952966368*win[0]*fin[5]*subvol; 
  out[6] += 0.1767766952966368*win[0]*fin[7]*subvol; 
  out[7] += 0.1767766952966368*win[0]*fin[8]*subvol; 
  out[8] += 0.1767766952966368*win[0]*fin[9]*subvol; 
  out[9] += 0.1767766952966368*win[0]*fin[10]*subvol; 
  out[10] += 0.1767766952966368*win[0]*fin[11]*subvol; 
  out[11] += 0.1767766952966368*win[0]*fin[12]*subvol; 
  out[12] += 0.1767766952966368*win[0]*fin[13]*subvol; 
  out[13] += 0.1767766952966368*win[0]*fin[14]*subvol; 
  out[14] += 0.1767766952966368*win[0]*fin[15]*subvol; 
  out[15] += 0.1767766952966368*win[0]*fin[16]*subvol; 
  out[16] += 0.1767766952966368*win[0]*fin[22]*subvol; 
  out[17] += 0.1767766952966368*win[0]*fin[23]*subvol; 
  out[18] += 0.1767766952966368*win[0]*fin[25]*subvol; 
  out[19] += 0.1767766952966368*win[0]*fin[26]*subvol; 
  out[20] += 0.1767766952966368*win[0]*fin[27]*subvol; 
  out[21] += 0.1767766952966368*win[0]*fin[28]*subvol; 
  out[22] += 0.1767766952966368*win[0]*fin[29]*subvol; 
  out[23] += 0.1767766952966368*win[0]*fin[30]*subvol; 
  out[24] += 0.1767766952966368*win[0]*fin[31]*subvol; 
  out[25] += 0.1767766952966368*win[0]*fin[32]*subvol; 
  out[26] += 0.1767766952966368*win[0]*fin[33]*subvol; 
  out[27] += 0.1767766952966368*win[0]*fin[34]*subvol; 
  out[28] += 0.1767766952966368*win[0]*fin[45]*subvol; 
  out[29] += 0.1767766952966368*win[0]*fin[46]*subvol; 
  out[30] += 0.1767766952966368*win[0]*fin[47]*subvol; 
  out[31] += 0.1767766952966368*win[0]*fin[48]*subvol; 
  out[32] += 0.1767766952966368*win[0]*fin[49]*subvol; 
  out[33] += 0.1767766952966368*win[0]*fin[50]*subvol; 
  out[34] += 0.1767766952966368*win[0]*fin[51]*subvol; 
  out[35] += 0.1767766952966368*win[0]*fin[52]*subvol; 
  out[36] += 0.1767766952966368*win[0]*fin[60]*subvol; 
  out[37] += 0.1767766952966368*win[0]*fin[61]*subvol; 
  out[38] += 0.1767766952966368*win[0]*fin[62]*subvol; 
  out[39] += 0.1767766952966368*win[0]*fin[63]*subvol; 
  out[40] += 0.1767766952966368*win[0]*fin[64]*subvol; 
  out[41] += 0.1767766952966368*win[0]*fin[75]*subvol; 
  out[42] += 0.1767766952966368*win[0]*fin[76]*subvol; 
  out[43] += 0.1767766952966368*win[0]*fin[77]*subvol; 
  out[44] += 0.1767766952966368*win[0]*fin[78]*subvol; 
  out[45] += 0.1767766952966368*win[0]*fin[79]*subvol; 
  out[46] += 0.1767766952966368*win[0]*fin[80]*subvol; 
  out[47] += 0.1767766952966368*win[0]*fin[81]*subvol; 
  out[48] += 0.1767766952966368*win[0]*fin[82]*subvol; 
  out[49] += 0.1767766952966368*win[0]*fin[83]*subvol; 
  out[50] += 0.1767766952966368*win[0]*fin[84]*subvol; 
  out[51] += 0.1767766952966368*win[0]*fin[85]*subvol; 
  out[52] += 0.1767766952966368*win[0]*fin[86]*subvol; 
  out[53] += 0.1767766952966368*win[0]*fin[105]*subvol; 
  out[54] += 0.1767766952966368*win[0]*fin[108]*subvol; 
  out[55] += 0.1767766952966368*win[0]*fin[114]*subvol; 
  out[56] += 0.1767766952966368*win[0]*fin[115]*subvol; 
  out[57] += 0.1767766952966368*win[0]*fin[116]*subvol; 
  out[58] += 0.1767766952966368*win[0]*fin[117]*subvol; 
  out[59] += 0.1767766952966368*win[0]*fin[118]*subvol; 
  out[60] += 0.1767766952966368*win[0]*fin[119]*subvol; 
  out[61] += 0.1767766952966368*win[0]*fin[120]*subvol; 
  out[62] += 0.1767766952966368*win[0]*fin[121]*subvol; 
  out[63] += 0.1767766952966368*win[0]*fin[144]*subvol; 
  out[64] += 0.1767766952966368*win[0]*fin[145]*subvol; 
  out[65] += 0.1767766952966368*win[0]*fin[146]*subvol; 
  out[66] += 0.1767766952966368*win[0]*fin[157]*subvol; 
  out[67] += 0.1767766952966368*win[0]*fin[158]*subvol; 
  out[68] += 0.1767766952966368*win[0]*fin[172]*subvol; 
  out[69] += 0.1767766952966368*win[0]*fin[173]*subvol; 
  out[70] += 0.1767766952966368*win[0]*fin[174]*subvol; 
  out[71] += 0.1767766952966368*win[0]*fin[194]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.25*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.25*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.25*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.25*win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], 0.25*win[0]*fin[4]*subvol); 
  atomicAdd(&out[5], 0.25*win[0]*fin[5]*subvol); 
  atomicAdd(&out[7], 0.25*win[0]*fin[7]*subvol); 
  atomicAdd(&out[8], 0.25*win[0]*fin[8]*subvol); 
  atomicAdd(&out[9], 0.25*win[0]*fin[9]*subvol); 
  atomicAdd(&out[10], 0.25*win[0]*fin[10]*subvol); 
  atomicAdd(&out[11], 0.25*win[0]*fin[11]*subvol); 
  atomicAdd(&out[12], 0.25*win[0]*fin[12]*subvol); 
  atomicAdd(&out[13], 0.25*win[0]*fin[13]*subvol); 
  atomicAdd(&out[14], 0.25*win[0]*fin[14]*subvol); 
  atomicAdd(&out[15], 0.25*win[0]*fin[15]*subvol); 
  atomicAdd(&out[16], 0.25*win[0]*fin[16]*subvol); 
  atomicAdd(&out[22], 0.25*win[0]*fin[22]*subvol); 
  atomicAdd(&out[23], 0.25*win[0]*fin[23]*subvol); 
  atomicAdd(&out[25], 0.25*win[0]*fin[25]*subvol); 
  atomicAdd(&out[26], 0.25*win[0]*fin[26]*subvol); 
  atomicAdd(&out[27], 0.25*win[0]*fin[27]*subvol); 
  atomicAdd(&out[28], 0.25*win[0]*fin[28]*subvol); 
  atomicAdd(&out[29], 0.25*win[0]*fin[29]*subvol); 
  atomicAdd(&out[30], 0.25*win[0]*fin[30]*subvol); 
  atomicAdd(&out[31], 0.25*win[0]*fin[31]*subvol); 
  atomicAdd(&out[32], 0.25*win[0]*fin[32]*subvol); 
  atomicAdd(&out[33], 0.25*win[0]*fin[33]*subvol); 
  atomicAdd(&out[34], 0.25*win[0]*fin[34]*subvol); 
  atomicAdd(&out[45], 0.25*win[0]*fin[45]*subvol); 
  atomicAdd(&out[46], 0.25*win[0]*fin[46]*subvol); 
  atomicAdd(&out[47], 0.25*win[0]*fin[47]*subvol); 
  atomicAdd(&out[48], 0.25*win[0]*fin[48]*subvol); 
  atomicAdd(&out[49], 0.25*win[0]*fin[49]*subvol); 
  atomicAdd(&out[50], 0.25*win[0]*fin[50]*subvol); 
  atomicAdd(&out[51], 0.25*win[0]*fin[51]*subvol); 
  atomicAdd(&out[52], 0.25*win[0]*fin[52]*subvol); 
  atomicAdd(&out[60], 0.25*win[0]*fin[60]*subvol); 
  atomicAdd(&out[61], 0.25*win[0]*fin[61]*subvol); 
  atomicAdd(&out[62], 0.25*win[0]*fin[62]*subvol); 
  atomicAdd(&out[63], 0.25*win[0]*fin[63]*subvol); 
  atomicAdd(&out[64], 0.25*win[0]*fin[64]*subvol); 
  atomicAdd(&out[75], 0.25*win[0]*fin[75]*subvol); 
  atomicAdd(&out[76], 0.25*win[0]*fin[76]*subvol); 
  atomicAdd(&out[77], 0.25*win[0]*fin[77]*subvol); 
  atomicAdd(&out[78], 0.25*win[0]*fin[78]*subvol); 
  atomicAdd(&out[79], 0.25*win[0]*fin[79]*subvol); 
  atomicAdd(&out[80], 0.25*win[0]*fin[80]*subvol); 
  atomicAdd(&out[81], 0.25*win[0]*fin[81]*subvol); 
  atomicAdd(&out[82], 0.25*win[0]*fin[82]*subvol); 
  atomicAdd(&out[83], 0.25*win[0]*fin[83]*subvol); 
  atomicAdd(&out[84], 0.25*win[0]*fin[84]*subvol); 
  atomicAdd(&out[85], 0.25*win[0]*fin[85]*subvol); 
  atomicAdd(&out[86], 0.25*win[0]*fin[86]*subvol); 
  atomicAdd(&out[105], 0.25*win[0]*fin[105]*subvol); 
  atomicAdd(&out[108], 0.25*win[0]*fin[108]*subvol); 
  atomicAdd(&out[114], 0.25*win[0]*fin[114]*subvol); 
  atomicAdd(&out[115], 0.25*win[0]*fin[115]*subvol); 
  atomicAdd(&out[116], 0.25*win[0]*fin[116]*subvol); 
  atomicAdd(&out[117], 0.25*win[0]*fin[117]*subvol); 
  atomicAdd(&out[118], 0.25*win[0]*fin[118]*subvol); 
  atomicAdd(&out[119], 0.25*win[0]*fin[119]*subvol); 
  atomicAdd(&out[120], 0.25*win[0]*fin[120]*subvol); 
  atomicAdd(&out[121], 0.25*win[0]*fin[121]*subvol); 
  atomicAdd(&out[144], 0.25*win[0]*fin[144]*subvol); 
  atomicAdd(&out[145], 0.25*win[0]*fin[145]*subvol); 
  atomicAdd(&out[146], 0.25*win[0]*fin[146]*subvol); 
  atomicAdd(&out[157], 0.25*win[0]*fin[157]*subvol); 
  atomicAdd(&out[158], 0.25*win[0]*fin[158]*subvol); 
  atomicAdd(&out[172], 0.25*win[0]*fin[172]*subvol); 
  atomicAdd(&out[173], 0.25*win[0]*fin[173]*subvol); 
  atomicAdd(&out[174], 0.25*win[0]*fin[174]*subvol); 
  atomicAdd(&out[194], 0.25*win[0]*fin[194]*subvol); 
#else
  out[0] += 0.25*fin[0]*win[0]*subvol; 
  out[1] += 0.25*win[0]*fin[1]*subvol; 
  out[2] += 0.25*win[0]*fin[2]*subvol; 
  out[3] += 0.25*win[0]*fin[3]*subvol; 
  out[4] += 0.25*win[0]*fin[4]*subvol; 
  out[5] += 0.25*win[0]*fin[5]*subvol; 
  out[7] += 0.25*win[0]*fin[7]*subvol; 
  out[8] += 0.25*win[0]*fin[8]*subvol; 
  out[9] += 0.25*win[0]*fin[9]*subvol; 
  out[10] += 0.25*win[0]*fin[10]*subvol; 
  out[11] += 0.25*win[0]*fin[11]*subvol; 
  out[12] += 0.25*win[0]*fin[12]*subvol; 
  out[13] += 0.25*win[0]*fin[13]*subvol; 
  out[14] += 0.25*win[0]*fin[14]*subvol; 
  out[15] += 0.25*win[0]*fin[15]*subvol; 
  out[16] += 0.25*win[0]*fin[16]*subvol; 
  out[22] += 0.25*win[0]*fin[22]*subvol; 
  out[23] += 0.25*win[0]*fin[23]*subvol; 
  out[25] += 0.25*win[0]*fin[25]*subvol; 
  out[26] += 0.25*win[0]*fin[26]*subvol; 
  out[27] += 0.25*win[0]*fin[27]*subvol; 
  out[28] += 0.25*win[0]*fin[28]*subvol; 
  out[29] += 0.25*win[0]*fin[29]*subvol; 
  out[30] += 0.25*win[0]*fin[30]*subvol; 
  out[31] += 0.25*win[0]*fin[31]*subvol; 
  out[32] += 0.25*win[0]*fin[32]*subvol; 
  out[33] += 0.25*win[0]*fin[33]*subvol; 
  out[34] += 0.25*win[0]*fin[34]*subvol; 
  out[45] += 0.25*win[0]*fin[45]*subvol; 
  out[46] += 0.25*win[0]*fin[46]*subvol; 
  out[47] += 0.25*win[0]*fin[47]*subvol; 
  out[48] += 0.25*win[0]*fin[48]*subvol; 
  out[49] += 0.25*win[0]*fin[49]*subvol; 
  out[50] += 0.25*win[0]*fin[50]*subvol; 
  out[51] += 0.25*win[0]*fin[51]*subvol; 
  out[52] += 0.25*win[0]*fin[52]*subvol; 
  out[60] += 0.25*win[0]*fin[60]*subvol; 
  out[61] += 0.25*win[0]*fin[61]*subvol; 
  out[62] += 0.25*win[0]*fin[62]*subvol; 
  out[63] += 0.25*win[0]*fin[63]*subvol; 
  out[64] += 0.25*win[0]*fin[64]*subvol; 
  out[75] += 0.25*win[0]*fin[75]*subvol; 
  out[76] += 0.25*win[0]*fin[76]*subvol; 
  out[77] += 0.25*win[0]*fin[77]*subvol; 
  out[78] += 0.25*win[0]*fin[78]*subvol; 
  out[79] += 0.25*win[0]*fin[79]*subvol; 
  out[80] += 0.25*win[0]*fin[80]*subvol; 
  out[81] += 0.25*win[0]*fin[81]*subvol; 
  out[82] += 0.25*win[0]*fin[82]*subvol; 
  out[83] += 0.25*win[0]*fin[83]*subvol; 
  out[84] += 0.25*win[0]*fin[84]*subvol; 
  out[85] += 0.25*win[0]*fin[85]*subvol; 
  out[86] += 0.25*win[0]*fin[86]*subvol; 
  out[105] += 0.25*win[0]*fin[105]*subvol; 
  out[108] += 0.25*win[0]*fin[108]*subvol; 
  out[114] += 0.25*win[0]*fin[114]*subvol; 
  out[115] += 0.25*win[0]*fin[115]*subvol; 
  out[116] += 0.25*win[0]*fin[116]*subvol; 
  out[117] += 0.25*win[0]*fin[117]*subvol; 
  out[118] += 0.25*win[0]*fin[118]*subvol; 
  out[119] += 0.25*win[0]*fin[119]*subvol; 
  out[120] += 0.25*win[0]*fin[120]*subvol; 
  out[121] += 0.25*win[0]*fin[121]*subvol; 
  out[144] += 0.25*win[0]*fin[144]*subvol; 
  out[145] += 0.25*win[0]*fin[145]*subvol; 
  out[146] += 0.25*win[0]*fin[146]*subvol; 
  out[157] += 0.25*win[0]*fin[157]*subvol; 
  out[158] += 0.25*win[0]*fin[158]*subvol; 
  out[172] += 0.25*win[0]*fin[172]*subvol; 
  out[173] += 0.25*win[0]*fin[173]*subvol; 
  out[174] += 0.25*win[0]*fin[174]*subvol; 
  out[194] += 0.25*win[0]*fin[194]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.25*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.25*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.25*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.25*win[0]*fin[4]*subvol); 
  atomicAdd(&out[4], 0.25*win[0]*fin[5]*subvol); 
  atomicAdd(&out[5], 0.25*win[0]*fin[9]*subvol); 
  atomicAdd(&out[6], 0.25*win[0]*fin[11]*subvol); 
  atomicAdd(&out[7], 0.25*win[0]*fin[12]*subvol); 
  atomicAdd(&out[8], 0.25*win[0]*fin[14]*subvol); 
  atomicAdd(&out[9], 0.25*win[0]*fin[15]*subvol); 
  atomicAdd(&out[10], 0.25*win[0]*fin[16]*subvol); 
  atomicAdd(&out[11], 0.25*win[0]*fin[22]*subvol); 
  atomicAdd(&out[12], 0.25*win[0]*fin[23]*subvol); 
  atomicAdd(&out[13], 0.25*win[0]*fin[28]*subvol); 
  atomicAdd(&out[14], 0.25*win[0]*fin[31]*subvol); 
  atomicAdd(&out[15], 0.25*win[0]*fin[33]*subvol); 
  atomicAdd(&out[16], 0.25*win[0]*fin[34]*subvol); 
  atomicAdd(&out[17], 0.25*win[0]*fin[46]*subvol); 
  atomicAdd(&out[18], 0.25*win[0]*fin[47]*subvol); 
  atomicAdd(&out[19], 0.25*win[0]*fin[48]*subvol); 
  atomicAdd(&out[20], 0.25*win[0]*fin[50]*subvol); 
  atomicAdd(&out[21], 0.25*win[0]*fin[51]*subvol); 
  atomicAdd(&out[22], 0.25*win[0]*fin[52]*subvol); 
  atomicAdd(&out[23], 0.25*win[0]*fin[64]*subvol); 
  atomicAdd(&out[24], 0.25*win[0]*fin[77]*subvol); 
  atomicAdd(&out[25], 0.25*win[0]*fin[79]*subvol); 
  atomicAdd(&out[26], 0.25*win[0]*fin[80]*subvol); 
  atomicAdd(&out[27], 0.25*win[0]*fin[83]*subvol); 
  atomicAdd(&out[28], 0.25*win[0]*fin[85]*subvol); 
  atomicAdd(&out[29], 0.25*win[0]*fin[86]*subvol); 
  atomicAdd(&out[30], 0.25*win[0]*fin[105]*subvol); 
  atomicAdd(&out[31], 0.25*win[0]*fin[117]*subvol); 
  atomicAdd(&out[32], 0.25*win[0]*fin[121]*subvol); 
  atomicAdd(&out[33], 0.25*win[0]*fin[145]*subvol); 
  atomicAdd(&out[34], 0.25*win[0]*fin[146]*subvol); 
  atomicAdd(&out[35], 0.25*win[0]*fin[174]*subvol); 
#else
  out[0] += 0.25*fin[0]*win[0]*subvol; 
  out[1] += 0.25*win[0]*fin[2]*subvol; 
  out[2] += 0.25*win[0]*fin[3]*subvol; 
  out[3] += 0.25*win[0]*fin[4]*subvol; 
  out[4] += 0.25*win[0]*fin[5]*subvol; 
  out[5] += 0.25*win[0]*fin[9]*subvol; 
  out[6] += 0.25*win[0]*fin[11]*subvol; 
  out[7] += 0.25*win[0]*fin[12]*subvol; 
  out[8] += 0.25*win[0]*fin[14]*subvol; 
  out[9] += 0.25*win[0]*fin[15]*subvol; 
  out[10] += 0.25*win[0]*fin[16]*subvol; 
  out[11] += 0.25*win[0]*fin[22]*subvol; 
  out[12] += 0.25*win[0]*fin[23]*subvol; 
  out[13] += 0.25*win[0]*fin[28]*subvol; 
  out[14] += 0.25*win[0]*fin[31]*subvol; 
  out[15] += 0.25*win[0]*fin[33]*subvol; 
  out[16] += 0.25*win[0]*fin[34]*subvol; 
  out[17] += 0.25*win[0]*fin[46]*subvol; 
  out[18] += 0.25*win[0]*fin[47]*subvol; 
  out[19] += 0.25*win[0]*fin[48]*subvol; 
  out[20] += 0.25*win[0]*fin[50]*subvol; 
  out[21] += 0.25*win[0]*fin[51]*subvol; 
  out[22] += 0.25*win[0]*fin[52]*subvol; 
  out[23] += 0.25*win[0]*fin[64]*subvol; 
  out[24] += 0.25*win[0]*fin[77]*subvol; 
  out[25] += 0.25*win[0]*fin[79]*subvol; 
  out[26] += 0.25*win[0]*fin[80]*subvol; 
  out[27] += 0.25*win[0]*fin[83]*subvol; 
  out[28] += 0.25*win[0]*fin[85]*subvol; 
  out[29] += 0.25*win[0]*fin[86]*subvol; 
  out[30] += 0.25*win[0]*fin[105]*subvol; 
  out[31] += 0.25*win[0]*fin[117]*subvol; 
  out[32] += 0.25*win[0]*fin[121]*subvol; 
  out[33] += 0.25*win[0]*fin[145]*subvol; 
  out[34] += 0.25*win[0]*fin[146]*subvol; 
  out[35] += 0.25*win[0]*fin[174]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[5]*subvol); 
  atomicAdd(&out[9], 0.5*win[0]*fin[9]*subvol); 
  atomicAdd(&out[11], 0.5*win[0]*fin[11]*subvol); 
  atomicAdd(&out[12], 0.5*win[0]*fin[12]*subvol); 
  atomicAdd(&out[14], 0.5*win[0]*fin[14]*subvol); 
  atomicAdd(&out[15], 0.5*win[0]*fin[15]*subvol); 
  atomicAdd(&out[16], 0.5*win[0]*fin[16]*subvol); 
  atomicAdd(&out[22], 0.5*win[0]*fin[22]*subvol); 
  atomicAdd(&out[23], 0.5*win[0]*fin[23]*subvol); 
  atomicAdd(&out[28], 0.5*win[0]*fin[28]*subvol); 
  atomicAdd(&out[31], 0.5*win[0]*fin[31]*subvol); 
  atomicAdd(&out[33], 0.5*win[0]*fin[33]*subvol); 
  atomicAdd(&out[34], 0.5*win[0]*fin[34]*subvol); 
  atomicAdd(&out[46], 0.5*win[0]*fin[46]*subvol); 
  atomicAdd(&out[47], 0.5*win[0]*fin[47]*subvol); 
  atomicAdd(&out[48], 0.5*win[0]*fin[48]*subvol); 
  atomicAdd(&out[50], 0.5*win[0]*fin[50]*subvol); 
  atomicAdd(&out[51], 0.5*win[0]*fin[51]*subvol); 
  atomicAdd(&out[52], 0.5*win[0]*fin[52]*subvol); 
  atomicAdd(&out[64], 0.5*win[0]*fin[64]*subvol); 
  atomicAdd(&out[77], 0.5*win[0]*fin[77]*subvol); 
  atomicAdd(&out[79], 0.5*win[0]*fin[79]*subvol); 
  atomicAdd(&out[80], 0.5*win[0]*fin[80]*subvol); 
  atomicAdd(&out[83], 0.5*win[0]*fin[83]*subvol); 
  atomicAdd(&out[85], 0.5*win[0]*fin[85]*subvol); 
  atomicAdd(&out[86], 0.5*win[0]*fin[86]*subvol); 
  atomicAdd(&out[105], 0.5*win[0]*fin[105]*subvol); 
  atomicAdd(&out[117], 0.5*win[0]*fin[117]*subvol); 
  atomicAdd(&out[121], 0.5*win[0]*fin[121]*subvol); 
  atomicAdd(&out[145], 0.5*win[0]*fin[145]*subvol); 
  atomicAdd(&out[146], 0.5*win[0]*fin[146]*subvol); 
  atomicAdd(&out[174], 0.5*win[0]*fin[174]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[2] += 0.5*win[0]*fin[2]*subvol; 
  out[3] += 0.5*win[0]*fin[3]*subvol; 
  out[4] += 0.5*win[0]*fin[4]*subvol; 
  out[5] += 0.5*win[0]*fin[5]*subvol; 
  out[9] += 0.5*win[0]*fin[9]*subvol; 
  out[11] += 0.5*win[0]*fin[11]*subvol; 
  out[12] += 0.5*win[0]*fin[12]*subvol; 
  out[14] += 0.5*win[0]*fin[14]*subvol; 
  out[15] += 0.5*win[0]*fin[15]*subvol; 
  out[16] += 0.5*win[0]*fin[16]*subvol; 
  out[22] += 0.5*win[0]*fin[22]*subvol; 
  out[23] += 0.5*win[0]*fin[23]*subvol; 
  out[28] += 0.5*win[0]*fin[28]*subvol; 
  out[31] += 0.5*win[0]*fin[31]*subvol; 
  out[33] += 0.5*win[0]*fin[33]*subvol; 
  out[34] += 0.5*win[0]*fin[34]*subvol; 
  out[46] += 0.5*win[0]*fin[46]*subvol; 
  out[47] += 0.5*win[0]*fin[47]*subvol; 
  out[48] += 0.5*win[0]*fin[48]*subvol; 
  out[50] += 0.5*win[0]*fin[50]*subvol; 
  out[51] += 0.5*win[0]*fin[51]*subvol; 
  out[52] += 0.5*win[0]*fin[52]*subvol; 
  out[64] += 0.5*win[0]*fin[64]*subvol; 
  out[77] += 0.5*win[0]*fin[77]*subvol; 
  out[79] += 0.5*win[0]*fin[79]*subvol; 
  out[80] += 0.5*win[0]*fin[80]*subvol; 
  out[83] += 0.5*win[0]*fin[83]*subvol; 
  out[85] += 0.5*win[0]*fin[85]*subvol; 
  out[86] += 0.5*win[0]*fin[86]*subvol; 
  out[105] += 0.5*win[0]*fin[105]*subvol; 
  out[117] += 0.5*win[0]*fin[117]*subvol; 
  out[121] += 0.5*win[0]*fin[121]*subvol; 
  out[145] += 0.5*win[0]*fin[145]*subvol; 
  out[146] += 0.5*win[0]*fin[146]*subvol; 
  out[174] += 0.5*win[0]*fin[174]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.25*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.25*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.25*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.25*win[0]*fin[4]*subvol); 
  atomicAdd(&out[4], 0.25*win[0]*fin[5]*subvol); 
  atomicAdd(&out[5], 0.25*win[0]*fin[8]*subvol); 
  atomicAdd(&out[6], 0.25*win[0]*fin[10]*subvol); 
  atomicAdd(&out[7], 0.25*win[0]*fin[12]*subvol); 
  atomicAdd(&out[8], 0.25*win[0]*fin[13]*subvol); 
  atomicAdd(&out[9], 0.25*win[0]*fin[15]*subvol); 
  atomicAdd(&out[10], 0.25*win[0]*fin[16]*subvol); 
  atomicAdd(&out[11], 0.25*win[0]*fin[22]*subvol); 
  atomicAdd(&out[12], 0.25*win[0]*fin[23]*subvol); 
  atomicAdd(&out[13], 0.25*win[0]*fin[27]*subvol); 
  atomicAdd(&out[14], 0.25*win[0]*fin[30]*subvol); 
  atomicAdd(&out[15], 0.25*win[0]*fin[32]*subvol); 
  atomicAdd(&out[16], 0.25*win[0]*fin[34]*subvol); 
  atomicAdd(&out[17], 0.25*win[0]*fin[45]*subvol); 
  atomicAdd(&out[18], 0.25*win[0]*fin[47]*subvol); 
  atomicAdd(&out[19], 0.25*win[0]*fin[48]*subvol); 
  atomicAdd(&out[20], 0.25*win[0]*fin[49]*subvol); 
  atomicAdd(&out[21], 0.25*win[0]*fin[51]*subvol); 
  atomicAdd(&out[22], 0.25*win[0]*fin[52]*subvol); 
  atomicAdd(&out[23], 0.25*win[0]*fin[63]*subvol); 
  atomicAdd(&out[24], 0.25*win[0]*fin[76]*subvol); 
  atomicAdd(&out[25], 0.25*win[0]*fin[78]*subvol); 
  atomicAdd(&out[26], 0.25*win[0]*fin[80]*subvol); 
  atomicAdd(&out[27], 0.25*win[0]*fin[82]*subvol); 
  atomicAdd(&out[28], 0.25*win[0]*fin[84]*subvol); 
  atomicAdd(&out[29], 0.25*win[0]*fin[86]*subvol); 
  atomicAdd(&out[30], 0.25*win[0]*fin[105]*subvol); 
  atomicAdd(&out[31], 0.25*win[0]*fin[116]*subvol); 
  atomicAdd(&out[32], 0.25*win[0]*fin[120]*subvol); 
  atomicAdd(&out[33], 0.25*win[0]*fin[144]*subvol); 
  atomicAdd(&out[34], 0.25*win[0]*fin[146]*subvol); 
  atomicAdd(&out[35], 0.25*win[0]*fin[173]*subvol); 
#else
  out[0] += 0.25*fin[0]*win[0]*subvol; 
  out[1] += 0.25*win[0]*fin[1]*subvol; 
  out[2] += 0.25*win[0]*fin[3]*subvol; 
  out[3] += 0.25*win[0]*fin[4]*subvol; 
  out[4] += 0.25*win[0]*fin[5]*subvol; 
  out[5] += 0.25*win[0]*fin[8]*subvol; 
  out[6] += 0.25*win[0]*fin[10]*subvol; 
  out[7] += 0.25*win[0]*fin[12]*subvol; 
  out[8] += 0.25*win[0]*fin[13]*subvol; 
  out[9] += 0.25*win[0]*fin[15]*subvol; 
  out[10] += 0.25*win[0]*fin[16]*subvol; 
  out[11] += 0.25*win[0]*fin[22]*subvol; 
  out[12] += 0.25*win[0]*fin[23]*subvol; 
  out[13] += 0.25*win[0]*fin[27]*subvol; 
  out[14] += 0.25*win[0]*fin[30]*subvol; 
  out[15] += 0.25*win[0]*fin[32]*subvol; 
  out[16] += 0.25*win[0]*fin[34]*subvol; 
  out[17] += 0.25*win[0]*fin[45]*subvol; 
  out[18] += 0.25*win[0]*fin[47]*subvol; 
  out[19] += 0.25*win[0]*fin[48]*subvol; 
  out[20] += 0.25*win[0]*fin[49]*subvol; 
  out[21] += 0.25*win[0]*fin[51]*subvol; 
  out[22] += 0.25*win[0]*fin[52]*subvol; 
  out[23] += 0.25*win[0]*fin[63]*subvol; 
  out[24] += 0.25*win[0]*fin[76]*subvol; 
  out[25] += 0.25*win[0]*fin[78]*subvol; 
  out[26] += 0.25*win[0]*fin[80]*subvol; 
  out[27] += 0.25*win[0]*fin[82]*subvol; 
  out[28] += 0.25*win[0]*fin[84]*subvol; 
  out[29] += 0.25*win[0]*fin[86]*subvol; 
  out[30] += 0.25*win[0]*fin[105]*subvol; 
  out[31] += 0.25*win[0]*fin[116]*subvol; 
  out[32] += 0.25*win[0]*fin[120]*subvol; 
  out[33] += 0.25*win[0]*fin[144]*subvol; 
  out[34] += 0.25*win[0]*fin[146]*subvol; 
  out[35] += 0.25*win[0]*fin[173]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[5]*subvol); 
  atomicAdd(&out[8], 0.5*win[0]*fin[8]*subvol); 
  atomicAdd(&out[10], 0.5*win[0]*fin[10]*subvol); 
  atomicAdd(&out[12], 0.5*win[0]*fin[12]*subvol); 
  atomicAdd(&out[13], 0.5*win[0]*fin[13]*subvol); 
  atomicAdd(&out[15], 0.5*win[0]*fin[15]*subvol); 
  atomicAdd(&out[16], 0.5*win[0]*fin[16]*subvol); 
  atomicAdd(&out[22], 0.5*win[0]*fin[22]*subvol); 
  atomicAdd(&out[23], 0.5*win[0]*fin[23]*subvol); 
  atomicAdd(&out[27], 0.5*win[0]*fin[27]*subvol); 
  atomicAdd(&out[30], 0.5*win[0]*fin[30]*subvol); 
  atomicAdd(&out[32], 0.5*win[0]*fin[32]*subvol); 
  atomicAdd(&out[34], 0.5*win[0]*fin[34]*subvol); 
  atomicAdd(&out[45], 0.5*win[0]*fin[45]*subvol); 
  atomicAdd(&out[47], 0.5*win[0]*fin[47]*subvol); 
  atomicAdd(&out[48], 0.5*win[0]*fin[48]*subvol); 
  atomicAdd(&out[49], 0.5*win[0]*fin[49]*subvol); 
  atomicAdd(&out[51], 0.5*win[0]*fin[51]*subvol); 
  atomicAdd(&out[52], 0.5*win[0]*fin[52]*subvol); 
  atomicAdd(&out[63], 0.5*win[0]*fin[63]*subvol); 
  atomicAdd(&out[76], 0.5*win[0]*fin[76]*subvol); 
  atomicAdd(&out[78], 0.5*win[0]*fin[78]*subvol); 
  atomicAdd(&out[80], 0.5*win[0]*fin[80]*subvol); 
  atomicAdd(&out[82], 0.5*win[0]*fin[82]*subvol); 
  atomicAdd(&out[84], 0.5*win[0]*fin[84]*subvol); 
  atomicAdd(&out[86], 0.5*win[0]*fin[86]*subvol); 
  atomicAdd(&out[105], 0.5*win[0]*fin[105]*subvol); 
  atomicAdd(&out[116], 0.5*win[0]*fin[116]*subvol); 
  atomicAdd(&out[120], 0.5*win[0]*fin[120]*subvol); 
  atomicAdd(&out[144], 0.5*win[0]*fin[144]*subvol); 
  atomicAdd(&out[146], 0.5*win[0]*fin[146]*subvol); 
  atomicAdd(&out[173], 0.5*win[0]*fin[173]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[3] += 0.5*win[0]*fin[3]*subvol; 
  out[4] += 0.5*win[0]*fin[4]*subvol; 
  out[5] += 0.5*win[0]*fin[5]*subvol; 
  out[8] += 0.5*win[0]*fin[8]*subvol; 
  out[10] += 0.5*win[0]*fin[10]*subvol; 
  out[12] += 0.5*win[0]*fin[12]*subvol; 
  out[13] += 0.5*win[0]*fin[13]*subvol; 
  out[15] += 0.5*win[0]*fin[15]*subvol; 
  out[16] += 0.5*win[0]*fin[16]*subvol; 
  out[22] += 0.5*win[0]*fin[22]*subvol; 
  out[23] += 0.5*win[0]*fin[23]*subvol; 
  out[27] += 0.5*win[0]*fin[27]*subvol; 
  out[30] += 0.5*win[0]*fin[30]*subvol; 
  out[32] += 0.5*win[0]*fin[32]*subvol; 
  out[34] += 0.5*win[0]*fin[34]*subvol; 
  out[45] += 0.5*win[0]*fin[45]*subvol; 
  out[47] += 0.5*win[0]*fin[47]*subvol; 
  out[48] += 0.5*win[0]*fin[48]*subvol; 
  out[49] += 0.5*win[0]*fin[49]*subvol; 
  out[51] += 0.5*win[0]*fin[51]*subvol; 
  out[52] += 0.5*win[0]*fin[52]*subvol; 
  out[63] += 0.5*win[0]*fin[63]*subvol; 
  out[76] += 0.5*win[0]*fin[76]*subvol; 
  out[78] += 0.5*win[0]*fin[78]*subvol; 
  out[80] += 0.5*win[0]*fin[80]*subvol; 
  out[82] += 0.5*win[0]*fin[82]*subvol; 
  out[84] += 0.5*win[0]*fin[84]*subvol; 
  out[86] += 0.5*win[0]*fin[86]*subvol; 
  out[105] += 0.5*win[0]*fin[105]*subvol; 
  out[116] += 0.5*win[0]*fin[116]*subvol; 
  out[120] += 0.5*win[0]*fin[120]*subvol; 
  out[144] += 0.5*win[0]*fin[144]*subvol; 
  out[146] += 0.5*win[0]*fin[146]*subvol; 
  out[173] += 0.5*win[0]*fin[173]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[3]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[4]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[5]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[12]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[15]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[16]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[22]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[23]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[34]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[47]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[48]*subvol); 
  atomicAdd(&out[12], 0.3535533905932737*win[0]*fin[51]*subvol); 
  atomicAdd(&out[13], 0.3535533905932737*win[0]*fin[52]*subvol); 
  atomicAdd(&out[14], 0.3535533905932737*win[0]*fin[80]*subvol); 
  atomicAdd(&out[15], 0.3535533905932737*win[0]*fin[86]*subvol); 
  atomicAdd(&out[16], 0.3535533905932737*win[0]*fin[105]*subvol); 
  atomicAdd(&out[17], 0.3535533905932737*win[0]*fin[146]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[3]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[4]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[5]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[12]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[15]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[16]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[22]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[23]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[34]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[47]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[48]*subvol; 
  out[12] += 0.3535533905932737*win[0]*fin[51]*subvol; 
  out[13] += 0.3535533905932737*win[0]*fin[52]*subvol; 
  out[14] += 0.3535533905932737*win[0]*fin[80]*subvol; 
  out[15] += 0.3535533905932737*win[0]*fin[86]*subvol; 
  out[16] += 0.3535533905932737*win[0]*fin[105]*subvol; 
  out[17] += 0.3535533905932737*win[0]*fin[146]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
  atomicAdd(&out[3], win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], win[0]*fin[4]*subvol); 
  atomicAdd(&out[5], win[0]*fin[5]*subvol); 
  atomicAdd(&out[12], win[0]*fin[12]*subvol); 
  atomicAdd(&out[15], win[0]*fin[15]*subvol); 
  atomicAdd(&out[16], win[0]*fin[16]*subvol); 
  atomicAdd(&out[22], win[0]*fin[22]*subvol); 
  atomicAdd(&out[23], win[0]*fin[23]*subvol); 
  atomicAdd(&out[34], win[0]*fin[34]*subvol); 
  atomicAdd(&out[47], win[0]*fin[47]*subvol); 
  atomicAdd(&out[48], win[0]*fin[48]*subvol); 
  atomicAdd(&out[51], win[0]*fin[51]*subvol); 
  atomicAdd(&out[52], win[0]*fin[52]*subvol); 
  atomicAdd(&out[80], win[0]*fin[80]*subvol); 
  atomicAdd(&out[86], win[0]*fin[86]*subvol); 
  atomicAdd(&out[105], win[0]*fin[105]*subvol); 
  atomicAdd(&out[146], win[0]*fin[146]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
  out[3] += win[0]*fin[3]*subvol; 
  out[4] += win[0]*fin[4]*subvol; 
  out[5] += win[0]*fin[5]*subvol; 
  out[12] += win[0]*fin[12]*subvol; 
  out[15] += win[0]*fin[15]*subvol; 
  out[16] += win[0]*fin[16]*subvol; 
  out[22] += win[0]*fin[22]*subvol; 
  out[23] += win[0]*fin[23]*subvol; 
  out[34] += win[0]*fin[34]*subvol; 
  out[47] += win[0]*fin[47]*subvol; 
  out[48] += win[0]*fin[48]*subvol; 
  out[51] += win[0]*fin[51]*subvol; 
  out[52] += win[0]*fin[52]*subvol; 
  out[80] += win[0]*fin[80]*subvol; 
  out[86] += win[0]*fin[86]*subvol; 
  out[105] += win[0]*fin[105]*subvol; 
  out[146] += win[0]*fin[146]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgzvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.25*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.25*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.25*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.25*win[0]*fin[4]*subvol); 
  atomicAdd(&out[4], 0.25*win[0]*fin[5]*subvol); 
  atomicAdd(&out[5], 0.25*win[0]*fin[7]*subvol); 
  atomicAdd(&out[6], 0.25*win[0]*fin[10]*subvol); 
  atomicAdd(&out[7], 0.25*win[0]*fin[11]*subvol); 
  atomicAdd(&out[8], 0.25*win[0]*fin[13]*subvol); 
  atomicAdd(&out[9], 0.25*win[0]*fin[14]*subvol); 
  atomicAdd(&out[10], 0.25*win[0]*fin[16]*subvol); 
  atomicAdd(&out[11], 0.25*win[0]*fin[22]*subvol); 
  atomicAdd(&out[12], 0.25*win[0]*fin[23]*subvol); 
  atomicAdd(&out[13], 0.25*win[0]*fin[26]*subvol); 
  atomicAdd(&out[14], 0.25*win[0]*fin[29]*subvol); 
  atomicAdd(&out[15], 0.25*win[0]*fin[32]*subvol); 
  atomicAdd(&out[16], 0.25*win[0]*fin[33]*subvol); 
  atomicAdd(&out[17], 0.25*win[0]*fin[45]*subvol); 
  atomicAdd(&out[18], 0.25*win[0]*fin[46]*subvol); 
  atomicAdd(&out[19], 0.25*win[0]*fin[48]*subvol); 
  atomicAdd(&out[20], 0.25*win[0]*fin[49]*subvol); 
  atomicAdd(&out[21], 0.25*win[0]*fin[50]*subvol); 
  atomicAdd(&out[22], 0.25*win[0]*fin[52]*subvol); 
  atomicAdd(&out[23], 0.25*win[0]*fin[62]*subvol); 
  atomicAdd(&out[24], 0.25*win[0]*fin[75]*subvol); 
  atomicAdd(&out[25], 0.25*win[0]*fin[78]*subvol); 
  atomicAdd(&out[26], 0.25*win[0]*fin[79]*subvol); 
  atomicAdd(&out[27], 0.25*win[0]*fin[81]*subvol); 
  atomicAdd(&out[28], 0.25*win[0]*fin[84]*subvol); 
  atomicAdd(&out[29], 0.25*win[0]*fin[85]*subvol); 
  atomicAdd(&out[30], 0.25*win[0]*fin[105]*subvol); 
  atomicAdd(&out[31], 0.25*win[0]*fin[115]*subvol); 
  atomicAdd(&out[32], 0.25*win[0]*fin[119]*subvol); 
  atomicAdd(&out[33], 0.25*win[0]*fin[144]*subvol); 
  atomicAdd(&out[34], 0.25*win[0]*fin[145]*subvol); 
  atomicAdd(&out[35], 0.25*win[0]*fin[172]*subvol); 
#else
  out[0] += 0.25*fin[0]*win[0]*subvol; 
  out[1] += 0.25*win[0]*fin[1]*subvol; 
  out[2] += 0.25*win[0]*fin[2]*subvol; 
  out[3] += 0.25*win[0]*fin[4]*subvol; 
  out[4] += 0.25*win[0]*fin[5]*subvol; 
  out[5] += 0.25*win[0]*fin[7]*subvol; 
  out[6] += 0.25*win[0]*fin[10]*subvol; 
  out[7] += 0.25*win[0]*fin[11]*subvol; 
  out[8] += 0.25*win[0]*fin[13]*subvol; 
  out[9] += 0.25*win[0]*fin[14]*subvol; 
  out[10] += 0.25*win[0]*fin[16]*subvol; 
  out[11] += 0.25*win[0]*fin[22]*subvol; 
  out[12] += 0.25*win[0]*fin[23]*subvol; 
  out[13] += 0.25*win[0]*fin[26]*subvol; 
  out[14] += 0.25*win[0]*fin[29]*subvol; 
  out[15] += 0.25*win[0]*fin[32]*subvol; 
  out[16] += 0.25*win[0]*fin[33]*subvol; 
  out[17] += 0.25*win[0]*fin[45]*subvol; 
  out[18] += 0.25*win[0]*fin[46]*subvol; 
  out[19] += 0.25*win[0]*fin[48]*subvol; 
  out[20] += 0.25*win[0]*fin[49]*subvol; 
  out[21] += 0.25*win[0]*fin[50]*subvol; 
  out[22] += 0.25*win[0]*fin[52]*subvol; 
  out[23] += 0.25*win[0]*fin[62]*subvol; 
  out[24] += 0.25*win[0]*fin[75]*subvol; 
  out[25] += 0.25*win[0]*fin[78]*subvol; 
  out[26] += 0.25*win[0]*fin[79]*subvol; 
  out[27] += 0.25*win[0]*fin[81]*subvol; 
  out[28] += 0.25*win[0]*fin[84]*subvol; 
  out[29] += 0.25*win[0]*fin[85]*subvol; 
  out[30] += 0.25*win[0]*fin[105]*subvol; 
  out[31] += 0.25*win[0]*fin[115]*subvol; 
  out[32] += 0.25*win[0]*fin[119]*subvol; 
  out[33] += 0.25*win[0]*fin[144]*subvol; 
  out[34] += 0.25*win[0]*fin[145]*subvol; 
  out[35] += 0.25*win[0]*fin[172]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgzvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[5]*subvol); 
  atomicAdd(&out[7], 0.5*win[0]*fin[7]*subvol); 
  atomicAdd(&out[10], 0.5*win[0]*fin[10]*subvol); 
  atomicAdd(&out[11], 0.5*win[0]*fin[11]*subvol); 
  atomicAdd(&out[13], 0.5*win[0]*fin[13]*subvol); 
  atomicAdd(&out[14], 0.5*win[0]*fin[14]*subvol); 
  atomicAdd(&out[16], 0.5*win[0]*fin[16]*subvol); 
  atomicAdd(&out[22], 0.5*win[0]*fin[22]*subvol); 
  atomicAdd(&out[23], 0.5*win[0]*fin[23]*subvol); 
  atomicAdd(&out[26], 0.5*win[0]*fin[26]*subvol); 
  atomicAdd(&out[29], 0.5*win[0]*fin[29]*subvol); 
  atomicAdd(&out[32], 0.5*win[0]*fin[32]*subvol); 
  atomicAdd(&out[33], 0.5*win[0]*fin[33]*subvol); 
  atomicAdd(&out[45], 0.5*win[0]*fin[45]*subvol); 
  atomicAdd(&out[46], 0.5*win[0]*fin[46]*subvol); 
  atomicAdd(&out[48], 0.5*win[0]*fin[48]*subvol); 
  atomicAdd(&out[49], 0.5*win[0]*fin[49]*subvol); 
  atomicAdd(&out[50], 0.5*win[0]*fin[50]*subvol); 
  atomicAdd(&out[52], 0.5*win[0]*fin[52]*subvol); 
  atomicAdd(&out[62], 0.5*win[0]*fin[62]*subvol); 
  atomicAdd(&out[75], 0.5*win[0]*fin[75]*subvol); 
  atomicAdd(&out[78], 0.5*win[0]*fin[78]*subvol); 
  atomicAdd(&out[79], 0.5*win[0]*fin[79]*subvol); 
  atomicAdd(&out[81], 0.5*win[0]*fin[81]*subvol); 
  atomicAdd(&out[84], 0.5*win[0]*fin[84]*subvol); 
  atomicAdd(&out[85], 0.5*win[0]*fin[85]*subvol); 
  atomicAdd(&out[105], 0.5*win[0]*fin[105]*subvol); 
  atomicAdd(&out[115], 0.5*win[0]*fin[115]*subvol); 
  atomicAdd(&out[119], 0.5*win[0]*fin[119]*subvol); 
  atomicAdd(&out[144], 0.5*win[0]*fin[144]*subvol); 
  atomicAdd(&out[145], 0.5*win[0]*fin[145]*subvol); 
  atomicAdd(&out[172], 0.5*win[0]*fin[172]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[2] += 0.5*win[0]*fin[2]*subvol; 
  out[4] += 0.5*win[0]*fin[4]*subvol; 
  out[5] += 0.5*win[0]*fin[5]*subvol; 
  out[7] += 0.5*win[0]*fin[7]*subvol; 
  out[10] += 0.5*win[0]*fin[10]*subvol; 
  out[11] += 0.5*win[0]*fin[11]*subvol; 
  out[13] += 0.5*win[0]*fin[13]*subvol; 
  out[14] += 0.5*win[0]*fin[14]*subvol; 
  out[16] += 0.5*win[0]*fin[16]*subvol; 
  out[22] += 0.5*win[0]*fin[22]*subvol; 
  out[23] += 0.5*win[0]*fin[23]*subvol; 
  out[26] += 0.5*win[0]*fin[26]*subvol; 
  out[29] += 0.5*win[0]*fin[29]*subvol; 
  out[32] += 0.5*win[0]*fin[32]*subvol; 
  out[33] += 0.5*win[0]*fin[33]*subvol; 
  out[45] += 0.5*win[0]*fin[45]*subvol; 
  out[46] += 0.5*win[0]*fin[46]*subvol; 
  out[48] += 0.5*win[0]*fin[48]*subvol; 
  out[49] += 0.5*win[0]*fin[49]*subvol; 
  out[50] += 0.5*win[0]*fin[50]*subvol; 
  out[52] += 0.5*win[0]*fin[52]*subvol; 
  out[62] += 0.5*win[0]*fin[62]*subvol; 
  out[75] += 0.5*win[0]*fin[75]*subvol; 
  out[78] += 0.5*win[0]*fin[78]*subvol; 
  out[79] += 0.5*win[0]*fin[79]*subvol; 
  out[81] += 0.5*win[0]*fin[81]*subvol; 
  out[84] += 0.5*win[0]*fin[84]*subvol; 
  out[85] += 0.5*win[0]*fin[85]*subvol; 
  out[105] += 0.5*win[0]*fin[105]*subvol; 
  out[115] += 0.5*win[0]*fin[115]*subvol; 
  out[119] += 0.5*win[0]*fin[119]*subvol; 
  out[144] += 0.5*win[0]*fin[144]*subvol; 
  out[145] += 0.5*win[0]*fin[145]*subvol; 
  out[172] += 0.5*win[0]*fin[172]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxzvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[4]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[5]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[11]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[14]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[16]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[22]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[23]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[33]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[46]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[48]*subvol); 
  atomicAdd(&out[12], 0.3535533905932737*win[0]*fin[50]*subvol); 
  atomicAdd(&out[13], 0.3535533905932737*win[0]*fin[52]*subvol); 
  atomicAdd(&out[14], 0.3535533905932737*win[0]*fin[79]*subvol); 
  atomicAdd(&out[15], 0.3535533905932737*win[0]*fin[85]*subvol); 
  atomicAdd(&out[16], 0.3535533905932737*win[0]*fin[105]*subvol); 
  atomicAdd(&out[17], 0.3535533905932737*win[0]*fin[145]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[2]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[4]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[5]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[11]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[14]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[16]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[22]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[23]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[33]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[46]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[48]*subvol; 
  out[12] += 0.3535533905932737*win[0]*fin[50]*subvol; 
  out[13] += 0.3535533905932737*win[0]*fin[52]*subvol; 
  out[14] += 0.3535533905932737*win[0]*fin[79]*subvol; 
  out[15] += 0.3535533905932737*win[0]*fin[85]*subvol; 
  out[16] += 0.3535533905932737*win[0]*fin[105]*subvol; 
  out[17] += 0.3535533905932737*win[0]*fin[145]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxzvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], win[0]*fin[2]*subvol); 
  atomicAdd(&out[4], win[0]*fin[4]*subvol); 
  atomicAdd(&out[5], win[0]*fin[5]*subvol); 
  atomicAdd(&out[11], win[0]*fin[11]*subvol); 
  atomicAdd(&out[14], win[0]*fin[14]*subvol); 
  atomicAdd(&out[16], win[0]*fin[16]*subvol); 
  atomicAdd(&out[22], win[0]*fin[22]*subvol); 
  atomicAdd(&out[23], win[0]*fin[23]*subvol); 
  atomicAdd(&out[33], win[0]*fin[33]*subvol); 
  atomicAdd(&out[46], win[0]*fin[46]*subvol); 
  atomicAdd(&out[48], win[0]*fin[48]*subvol); 
  atomicAdd(&out[50], win[0]*fin[50]*subvol); 
  atomicAdd(&out[52], win[0]*fin[52]*subvol); 
  atomicAdd(&out[79], win[0]*fin[79]*subvol); 
  atomicAdd(&out[85], win[0]*fin[85]*subvol); 
  atomicAdd(&out[105], win[0]*fin[105]*subvol); 
  atomicAdd(&out[145], win[0]*fin[145]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
  out[2] += win[0]*fin[2]*subvol; 
  out[4] += win[0]*fin[4]*subvol; 
  out[5] += win[0]*fin[5]*subvol; 
  out[11] += win[0]*fin[11]*subvol; 
  out[14] += win[0]*fin[14]*subvol; 
  out[16] += win[0]*fin[16]*subvol; 
  out[22] += win[0]*fin[22]*subvol; 
  out[23] += win[0]*fin[23]*subvol; 
  out[33] += win[0]*fin[33]*subvol; 
  out[46] += win[0]*fin[46]*subvol; 
  out[48] += win[0]*fin[48]*subvol; 
  out[50] += win[0]*fin[50]*subvol; 
  out[52] += win[0]*fin[52]*subvol; 
  out[79] += win[0]*fin[79]*subvol; 
  out[85] += win[0]*fin[85]*subvol; 
  out[105] += win[0]*fin[105]*subvol; 
  out[145] += win[0]*fin[145]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgyzvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[4]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[5]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[10]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[13]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[16]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[22]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[23]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[32]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[45]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[48]*subvol); 
  atomicAdd(&out[12], 0.3535533905932737*win[0]*fin[49]*subvol); 
  atomicAdd(&out[13], 0.3535533905932737*win[0]*fin[52]*subvol); 
  atomicAdd(&out[14], 0.3535533905932737*win[0]*fin[78]*subvol); 
  atomicAdd(&out[15], 0.3535533905932737*win[0]*fin[84]*subvol); 
  atomicAdd(&out[16], 0.3535533905932737*win[0]*fin[105]*subvol); 
  atomicAdd(&out[17], 0.3535533905932737*win[0]*fin[144]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[1]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[4]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[5]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[10]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[13]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[16]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[22]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[23]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[32]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[45]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[48]*subvol; 
  out[12] += 0.3535533905932737*win[0]*fin[49]*subvol; 
  out[13] += 0.3535533905932737*win[0]*fin[52]*subvol; 
  out[14] += 0.3535533905932737*win[0]*fin[78]*subvol; 
  out[15] += 0.3535533905932737*win[0]*fin[84]*subvol; 
  out[16] += 0.3535533905932737*win[0]*fin[105]*subvol; 
  out[17] += 0.3535533905932737*win[0]*fin[144]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgyzvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], win[0]*fin[1]*subvol); 
  atomicAdd(&out[4], win[0]*fin[4]*subvol); 
  atomicAdd(&out[5], win[0]*fin[5]*subvol); 
  atomicAdd(&out[10], win[0]*fin[10]*subvol); 
  atomicAdd(&out[13], win[0]*fin[13]*subvol); 
  atomicAdd(&out[16], win[0]*fin[16]*subvol); 
  atomicAdd(&out[22], win[0]*fin[22]*subvol); 
  atomicAdd(&out[23], win[0]*fin[23]*subvol); 
  atomicAdd(&out[32], win[0]*fin[32]*subvol); 
  atomicAdd(&out[45], win[0]*fin[45]*subvol); 
  atomicAdd(&out[48], win[0]*fin[48]*subvol); 
  atomicAdd(&out[49], win[0]*fin[49]*subvol); 
  atomicAdd(&out[52], win[0]*fin[52]*subvol); 
  atomicAdd(&out[78], win[0]*fin[78]*subvol); 
  atomicAdd(&out[84], win[0]*fin[84]*subvol); 
  atomicAdd(&out[105], win[0]*fin[105]*subvol); 
  atomicAdd(&out[144], win[0]*fin[144]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
  out[1] += win[0]*fin[1]*subvol; 
  out[4] += win[0]*fin[4]*subvol; 
  out[5] += win[0]*fin[5]*subvol; 
  out[10] += win[0]*fin[10]*subvol; 
  out[13] += win[0]*fin[13]*subvol; 
  out[16] += win[0]*fin[16]*subvol; 
  out[22] += win[0]*fin[22]*subvol; 
  out[23] += win[0]*fin[23]*subvol; 
  out[32] += win[0]*fin[32]*subvol; 
  out[45] += win[0]*fin[45]*subvol; 
  out[48] += win[0]*fin[48]*subvol; 
  out[49] += win[0]*fin[49]*subvol; 
  out[52] += win[0]*fin[52]*subvol; 
  out[78] += win[0]*fin[78]*subvol; 
  out[84] += win[0]*fin[84]*subvol; 
  out[105] += win[0]*fin[105]*subvol; 
  out[144] += win[0]*fin[144]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxyzvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[5]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[16]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[22]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[23]*subvol); 
  atomicAdd(&out[6], 0.5*win[0]*fin[48]*subvol); 
  atomicAdd(&out[7], 0.5*win[0]*fin[52]*subvol); 
  atomicAdd(&out[8], 0.5*win[0]*fin[105]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[4]*subvol; 
  out[2] += 0.5*win[0]*fin[5]*subvol; 
  out[3] += 0.5*win[0]*fin[16]*subvol; 
  out[4] += 0.5*win[0]*fin[22]*subvol; 
  out[5] += 0.5*win[0]*fin[23]*subvol; 
  out[6] += 0.5*win[0]*fin[48]*subvol; 
  out[7] += 0.5*win[0]*fin[52]*subvol; 
  out[8] += 0.5*win[0]*fin[105]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxyzvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.0*fin[0]*win[0]*subvol); 
  atomicAdd(&out[4], 2.0*win[0]*fin[4]*subvol); 
  atomicAdd(&out[5], 2.0*win[0]*fin[5]*subvol); 
  atomicAdd(&out[16], 2.0*win[0]*fin[16]*subvol); 
  atomicAdd(&out[22], 2.0*win[0]*fin[22]*subvol); 
  atomicAdd(&out[23], 2.0*win[0]*fin[23]*subvol); 
  atomicAdd(&out[48], 2.0*win[0]*fin[48]*subvol); 
  atomicAdd(&out[52], 2.0*win[0]*fin[52]*subvol); 
  atomicAdd(&out[105], 2.0*win[0]*fin[105]*subvol); 
#else
  out[0] += 2.0*fin[0]*win[0]*subvol; 
  out[4] += 2.0*win[0]*fin[4]*subvol; 
  out[5] += 2.0*win[0]*fin[5]*subvol; 
  out[16] += 2.0*win[0]*fin[16]*subvol; 
  out[22] += 2.0*win[0]*fin[22]*subvol; 
  out[23] += 2.0*win[0]*fin[23]*subvol; 
  out[48] += 2.0*win[0]*fin[48]*subvol; 
  out[52] += 2.0*win[0]*fin[52]*subvol; 
  out[105] += 2.0*win[0]*fin[105]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgvxvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.25*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.25*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.25*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.25*win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], 0.25*win[0]*fin[5]*subvol); 
  atomicAdd(&out[5], 0.25*win[0]*fin[7]*subvol); 
  atomicAdd(&out[6], 0.25*win[0]*fin[8]*subvol); 
  atomicAdd(&out[7], 0.25*win[0]*fin[9]*subvol); 
  atomicAdd(&out[8], 0.25*win[0]*fin[13]*subvol); 
  atomicAdd(&out[9], 0.25*win[0]*fin[14]*subvol); 
  atomicAdd(&out[10], 0.25*win[0]*fin[15]*subvol); 
  atomicAdd(&out[11], 0.25*win[0]*fin[23]*subvol); 
  atomicAdd(&out[12], 0.25*win[0]*fin[25]*subvol); 
  atomicAdd(&out[13], 0.25*win[0]*fin[29]*subvol); 
  atomicAdd(&out[14], 0.25*win[0]*fin[30]*subvol); 
  atomicAdd(&out[15], 0.25*win[0]*fin[31]*subvol); 
  atomicAdd(&out[16], 0.25*win[0]*fin[49]*subvol); 
  atomicAdd(&out[17], 0.25*win[0]*fin[50]*subvol); 
  atomicAdd(&out[18], 0.25*win[0]*fin[51]*subvol); 
  atomicAdd(&out[19], 0.25*win[0]*fin[61]*subvol); 
  atomicAdd(&out[20], 0.25*win[0]*fin[81]*subvol); 
  atomicAdd(&out[21], 0.25*win[0]*fin[82]*subvol); 
  atomicAdd(&out[22], 0.25*win[0]*fin[83]*subvol); 
  atomicAdd(&out[23], 0.25*win[0]*fin[118]*subvol); 
#else
  out[0] += 0.25*fin[0]*win[0]*subvol; 
  out[1] += 0.25*win[0]*fin[1]*subvol; 
  out[2] += 0.25*win[0]*fin[2]*subvol; 
  out[3] += 0.25*win[0]*fin[3]*subvol; 
  out[4] += 0.25*win[0]*fin[5]*subvol; 
  out[5] += 0.25*win[0]*fin[7]*subvol; 
  out[6] += 0.25*win[0]*fin[8]*subvol; 
  out[7] += 0.25*win[0]*fin[9]*subvol; 
  out[8] += 0.25*win[0]*fin[13]*subvol; 
  out[9] += 0.25*win[0]*fin[14]*subvol; 
  out[10] += 0.25*win[0]*fin[15]*subvol; 
  out[11] += 0.25*win[0]*fin[23]*subvol; 
  out[12] += 0.25*win[0]*fin[25]*subvol; 
  out[13] += 0.25*win[0]*fin[29]*subvol; 
  out[14] += 0.25*win[0]*fin[30]*subvol; 
  out[15] += 0.25*win[0]*fin[31]*subvol; 
  out[16] += 0.25*win[0]*fin[49]*subvol; 
  out[17] += 0.25*win[0]*fin[50]*subvol; 
  out[18] += 0.25*win[0]*fin[51]*subvol; 
  out[19] += 0.25*win[0]*fin[61]*subvol; 
  out[20] += 0.25*win[0]*fin[81]*subvol; 
  out[21] += 0.25*win[0]*fin[82]*subvol; 
  out[22] += 0.25*win[0]*fin[83]*subvol; 
  out[23] += 0.25*win[0]*fin[118]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgvxvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[5]*subvol); 
  atomicAdd(&out[7], 0.5*win[0]*fin[7]*subvol); 
  atomicAdd(&out[8], 0.5*win[0]*fin[8]*subvol); 
  atomicAdd(&out[9], 0.5*win[0]*fin[9]*subvol); 
  atomicAdd(&out[13], 0.5*win[0]*fin[13]*subvol); 
  atomicAdd(&out[14], 0.5*win[0]*fin[14]*subvol); 
  atomicAdd(&out[15], 0.5*win[0]*fin[15]*subvol); 
  atomicAdd(&out[23], 0.5*win[0]*fin[23]*subvol); 
  atomicAdd(&out[25], 0.5*win[0]*fin[25]*subvol); 
  atomicAdd(&out[29], 0.5*win[0]*fin[29]*subvol); 
  atomicAdd(&out[30], 0.5*win[0]*fin[30]*subvol); 
  atomicAdd(&out[31], 0.5*win[0]*fin[31]*subvol); 
  atomicAdd(&out[49], 0.5*win[0]*fin[49]*subvol); 
  atomicAdd(&out[50], 0.5*win[0]*fin[50]*subvol); 
  atomicAdd(&out[51], 0.5*win[0]*fin[51]*subvol); 
  atomicAdd(&out[61], 0.5*win[0]*fin[61]*subvol); 
  atomicAdd(&out[81], 0.5*win[0]*fin[81]*subvol); 
  atomicAdd(&out[82], 0.5*win[0]*fin[82]*subvol); 
  atomicAdd(&out[83], 0.5*win[0]*fin[83]*subvol); 
  atomicAdd(&out[118], 0.5*win[0]*fin[118]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[2] += 0.5*win[0]*fin[2]*subvol; 
  out[3] += 0.5*win[0]*fin[3]*subvol; 
  out[5] += 0.5*win[0]*fin[5]*subvol; 
  out[7] += 0.5*win[0]*fin[7]*subvol; 
  out[8] += 0.5*win[0]*fin[8]*subvol; 
  out[9] += 0.5*win[0]*fin[9]*subvol; 
  out[13] += 0.5*win[0]*fin[13]*subvol; 
  out[14] += 0.5*win[0]*fin[14]*subvol; 
  out[15] += 0.5*win[0]*fin[15]*subvol; 
  out[23] += 0.5*win[0]*fin[23]*subvol; 
  out[25] += 0.5*win[0]*fin[25]*subvol; 
  out[29] += 0.5*win[0]*fin[29]*subvol; 
  out[30] += 0.5*win[0]*fin[30]*subvol; 
  out[31] += 0.5*win[0]*fin[31]*subvol; 
  out[49] += 0.5*win[0]*fin[49]*subvol; 
  out[50] += 0.5*win[0]*fin[50]*subvol; 
  out[51] += 0.5*win[0]*fin[51]*subvol; 
  out[61] += 0.5*win[0]*fin[61]*subvol; 
  out[81] += 0.5*win[0]*fin[81]*subvol; 
  out[82] += 0.5*win[0]*fin[82]*subvol; 
  out[83] += 0.5*win[0]*fin[83]*subvol; 
  out[118] += 0.5*win[0]*fin[118]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxvxvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[5]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[9]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[14]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[15]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[23]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[31]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[50]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[51]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[83]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[2]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[3]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[5]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[9]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[14]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[15]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[23]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[31]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[50]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[51]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[83]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxvxvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], win[0]*fin[3]*subvol); 
  atomicAdd(&out[5], win[0]*fin[5]*subvol); 
  atomicAdd(&out[9], win[0]*fin[9]*subvol); 
  atomicAdd(&out[14], win[0]*fin[14]*subvol); 
  atomicAdd(&out[15], win[0]*fin[15]*subvol); 
  atomicAdd(&out[23], win[0]*fin[23]*subvol); 
  atomicAdd(&out[31], win[0]*fin[31]*subvol); 
  atomicAdd(&out[50], win[0]*fin[50]*subvol); 
  atomicAdd(&out[51], win[0]*fin[51]*subvol); 
  atomicAdd(&out[83], win[0]*fin[83]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
  out[2] += win[0]*fin[2]*subvol; 
  out[3] += win[0]*fin[3]*subvol; 
  out[5] += win[0]*fin[5]*subvol; 
  out[9] += win[0]*fin[9]*subvol; 
  out[14] += win[0]*fin[14]*subvol; 
  out[15] += win[0]*fin[15]*subvol; 
  out[23] += win[0]*fin[23]*subvol; 
  out[31] += win[0]*fin[31]*subvol; 
  out[50] += win[0]*fin[50]*subvol; 
  out[51] += win[0]*fin[51]*subvol; 
  out[83] += win[0]*fin[83]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgyvxvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[5]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[8]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[13]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[15]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[23]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[30]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[49]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[51]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[82]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[1]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[3]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[5]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[8]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[13]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[15]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[23]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[30]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[49]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[51]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[82]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgyvxvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], win[0]*fin[1]*subvol); 
  atomicAdd(&out[3], win[0]*fin[3]*subvol); 
  atomicAdd(&out[5], win[0]*fin[5]*subvol); 
  atomicAdd(&out[8], win[0]*fin[8]*subvol); 
  atomicAdd(&out[13], win[0]*fin[13]*subvol); 
  atomicAdd(&out[15], win[0]*fin[15]*subvol); 
  atomicAdd(&out[23], win[0]*fin[23]*subvol); 
  atomicAdd(&out[30], win[0]*fin[30]*subvol); 
  atomicAdd(&out[49], win[0]*fin[49]*subvol); 
  atomicAdd(&out[51], win[0]*fin[51]*subvol); 
  atomicAdd(&out[82], win[0]*fin[82]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
  out[1] += win[0]*fin[1]*subvol; 
  out[3] += win[0]*fin[3]*subvol; 
  out[5] += win[0]*fin[5]*subvol; 
  out[8] += win[0]*fin[8]*subvol; 
  out[13] += win[0]*fin[13]*subvol; 
  out[15] += win[0]*fin[15]*subvol; 
  out[23] += win[0]*fin[23]*subvol; 
  out[30] += win[0]*fin[30]*subvol; 
  out[49] += win[0]*fin[49]*subvol; 
  out[51] += win[0]*fin[51]*subvol; 
  out[82] += win[0]*fin[82]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxyvxvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[5]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[15]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[23]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[51]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[3]*subvol; 
  out[2] += 0.5*win[0]*fin[5]*subvol; 
  out[3] += 0.5*win[0]*fin[15]*subvol; 
  out[4] += 0.5*win[0]*fin[23]*subvol; 
  out[5] += 0.5*win[0]*fin[51]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxyvxvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.0*fin[0]*win[0]*subvol); 
  atomicAdd(&out[3], 2.0*win[0]*fin[3]*subvol); 
  atomicAdd(&out[5], 2.0*win[0]*fin[5]*subvol); 
  atomicAdd(&out[15], 2.0*win[0]*fin[15]*subvol); 
  atomicAdd(&out[23], 2.0*win[0]*fin[23]*subvol); 
  atomicAdd(&out[51], 2.0*win[0]*fin[51]*subvol); 
#else
  out[0] += 2.0*fin[0]*win[0]*subvol; 
  out[3] += 2.0*win[0]*fin[3]*subvol; 
  out[5] += 2.0*win[0]*fin[5]*subvol; 
  out[15] += 2.0*win[0]*fin[15]*subvol; 
  out[23] += 2.0*win[0]*fin[23]*subvol; 
  out[51] += 2.0*win[0]*fin[51]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgzvxvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[5]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[7]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[13]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[14]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[23]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[29]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[49]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[50]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[81]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[1]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[2]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[5]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[7]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[13]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[14]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[23]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[29]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[49]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[50]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[81]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgzvxvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], win[0]*fin[2]*subvol); 
  atomicAdd(&out[5], win[0]*fin[5]*subvol); 
  atomicAdd(&out[7], win[0]*fin[7]*subvol); 
  atomicAdd(&out[13], win[0]*fin[13]*subvol); 
  atomicAdd(&out[14], win[0]*fin[14]*subvol); 
  atomicAdd(&out[23], win[0]*fin[23]*subvol); 
  atomicAdd(&out[29], win[0]*fin[29]*subvol); 
  atomicAdd(&out[49], win[0]*fin[49]*subvol); 
  atomicAdd(&out[50], win[0]*fin[50]*subvol); 
  atomicAdd(&out[81], win[0]*fin[81]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
  out[1] += win[0]*fin[1]*subvol; 
  out[2] += win[0]*fin[2]*subvol; 
  out[5] += win[0]*fin[5]*subvol; 
  out[7] += win[0]*fin[7]*subvol; 
  out[13] += win[0]*fin[13]*subvol; 
  out[14] += win[0]*fin[14]*subvol; 
  out[23] += win[0]*fin[23]*subvol; 
  out[29] += win[0]*fin[29]*subvol; 
  out[49] += win[0]*fin[49]*subvol; 
  out[50] += win[0]*fin[50]*subvol; 
  out[81] += win[0]*fin[81]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxzvxvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[5]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[14]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[23]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[50]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[2]*subvol; 
  out[2] += 0.5*win[0]*fin[5]*subvol; 
  out[3] += 0.5*win[0]*fin[14]*subvol; 
  out[4] += 0.5*win[0]*fin[23]*subvol; 
  out[5] += 0.5*win[0]*fin[50]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxzvxvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.0*fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], 2.0*win[0]*fin[2]*subvol); 
  atomicAdd(&out[5], 2.0*win[0]*fin[5]*subvol); 
  atomicAdd(&out[14], 2.0*win[0]*fin[14]*subvol); 
  atomicAdd(&out[23], 2.0*win[0]*fin[23]*subvol); 
  atomicAdd(&out[50], 2.0*win[0]*fin[50]*subvol); 
#else
  out[0] += 2.0*fin[0]*win[0]*subvol; 
  out[2] += 2.0*win[0]*fin[2]*subvol; 
  out[5] += 2.0*win[0]*fin[5]*subvol; 
  out[14] += 2.0*win[0]*fin[14]*subvol; 
  out[23] += 2.0*win[0]*fin[23]*subvol; 
  out[50] += 2.0*win[0]*fin[50]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgyzvxvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[5]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[13]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[23]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[49]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[2] += 0.5*win[0]*fin[5]*subvol; 
  out[3] += 0.5*win[0]*fin[13]*subvol; 
  out[4] += 0.5*win[0]*fin[23]*subvol; 
  out[5] += 0.5*win[0]*fin[49]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgyzvxvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.0*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 2.0*win[0]*fin[1]*subvol); 
  atomicAdd(&out[5], 2.0*win[0]*fin[5]*subvol); 
  atomicAdd(&out[13], 2.0*win[0]*fin[13]*subvol); 
  atomicAdd(&out[23], 2.0*win[0]*fin[23]*subvol); 
  atomicAdd(&out[49], 2.0*win[0]*fin[49]*subvol); 
#else
  out[0] += 2.0*fin[0]*win[0]*subvol; 
  out[1] += 2.0*win[0]*fin[1]*subvol; 
  out[5] += 2.0*win[0]*fin[5]*subvol; 
  out[13] += 2.0*win[0]*fin[13]*subvol; 
  out[23] += 2.0*win[0]*fin[23]*subvol; 
  out[49] += 2.0*win[0]*fin[49]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxyzvxvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[5]*subvol); 
  atomicAdd(&out[2], 0.7071067811865475*win[0]*fin[23]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[5]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[23]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxyzvxvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 4.0*fin[0]*win[0]*subvol); 
  atomicAdd(&out[5], 4.0*win[0]*fin[5]*subvol); 
  atomicAdd(&out[23], 4.0*win[0]*fin[23]*subvol); 
#else
  out[0] += 4.0*fin[0]*win[0]*subvol; 
  out[5] += 4.0*win[0]*fin[5]*subvol; 
  out[23] += 4.0*win[0]*fin[23]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgvyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.25*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.25*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.25*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.25*win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], 0.25*win[0]*fin[4]*subvol); 
  atomicAdd(&out[5], 0.25*win[0]*fin[7]*subvol); 
  atomicAdd(&out[6], 0.25*win[0]*fin[8]*subvol); 
  atomicAdd(&out[7], 0.25*win[0]*fin[9]*subvol); 
  atomicAdd(&out[8], 0.25*win[0]*fin[10]*subvol); 
  atomicAdd(&out[9], 0.25*win[0]*fin[11]*subvol); 
  atomicAdd(&out[10], 0.25*win[0]*fin[12]*subvol); 
  atomicAdd(&out[11], 0.25*win[0]*fin[22]*subvol); 
  atomicAdd(&out[12], 0.25*win[0]*fin[25]*subvol); 
  atomicAdd(&out[13], 0.25*win[0]*fin[26]*subvol); 
  atomicAdd(&out[14], 0.25*win[0]*fin[27]*subvol); 
  atomicAdd(&out[15], 0.25*win[0]*fin[28]*subvol); 
  atomicAdd(&out[16], 0.25*win[0]*fin[45]*subvol); 
  atomicAdd(&out[17], 0.25*win[0]*fin[46]*subvol); 
  atomicAdd(&out[18], 0.25*win[0]*fin[47]*subvol); 
  atomicAdd(&out[19], 0.25*win[0]*fin[60]*subvol); 
  atomicAdd(&out[20], 0.25*win[0]*fin[75]*subvol); 
  atomicAdd(&out[21], 0.25*win[0]*fin[76]*subvol); 
  atomicAdd(&out[22], 0.25*win[0]*fin[77]*subvol); 
  atomicAdd(&out[23], 0.25*win[0]*fin[114]*subvol); 
#else
  out[0] += 0.25*fin[0]*win[0]*subvol; 
  out[1] += 0.25*win[0]*fin[1]*subvol; 
  out[2] += 0.25*win[0]*fin[2]*subvol; 
  out[3] += 0.25*win[0]*fin[3]*subvol; 
  out[4] += 0.25*win[0]*fin[4]*subvol; 
  out[5] += 0.25*win[0]*fin[7]*subvol; 
  out[6] += 0.25*win[0]*fin[8]*subvol; 
  out[7] += 0.25*win[0]*fin[9]*subvol; 
  out[8] += 0.25*win[0]*fin[10]*subvol; 
  out[9] += 0.25*win[0]*fin[11]*subvol; 
  out[10] += 0.25*win[0]*fin[12]*subvol; 
  out[11] += 0.25*win[0]*fin[22]*subvol; 
  out[12] += 0.25*win[0]*fin[25]*subvol; 
  out[13] += 0.25*win[0]*fin[26]*subvol; 
  out[14] += 0.25*win[0]*fin[27]*subvol; 
  out[15] += 0.25*win[0]*fin[28]*subvol; 
  out[16] += 0.25*win[0]*fin[45]*subvol; 
  out[17] += 0.25*win[0]*fin[46]*subvol; 
  out[18] += 0.25*win[0]*fin[47]*subvol; 
  out[19] += 0.25*win[0]*fin[60]*subvol; 
  out[20] += 0.25*win[0]*fin[75]*subvol; 
  out[21] += 0.25*win[0]*fin[76]*subvol; 
  out[22] += 0.25*win[0]*fin[77]*subvol; 
  out[23] += 0.25*win[0]*fin[114]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgvyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[7], 0.5*win[0]*fin[7]*subvol); 
  atomicAdd(&out[8], 0.5*win[0]*fin[8]*subvol); 
  atomicAdd(&out[9], 0.5*win[0]*fin[9]*subvol); 
  atomicAdd(&out[10], 0.5*win[0]*fin[10]*subvol); 
  atomicAdd(&out[11], 0.5*win[0]*fin[11]*subvol); 
  atomicAdd(&out[12], 0.5*win[0]*fin[12]*subvol); 
  atomicAdd(&out[22], 0.5*win[0]*fin[22]*subvol); 
  atomicAdd(&out[25], 0.5*win[0]*fin[25]*subvol); 
  atomicAdd(&out[26], 0.5*win[0]*fin[26]*subvol); 
  atomicAdd(&out[27], 0.5*win[0]*fin[27]*subvol); 
  atomicAdd(&out[28], 0.5*win[0]*fin[28]*subvol); 
  atomicAdd(&out[45], 0.5*win[0]*fin[45]*subvol); 
  atomicAdd(&out[46], 0.5*win[0]*fin[46]*subvol); 
  atomicAdd(&out[47], 0.5*win[0]*fin[47]*subvol); 
  atomicAdd(&out[60], 0.5*win[0]*fin[60]*subvol); 
  atomicAdd(&out[75], 0.5*win[0]*fin[75]*subvol); 
  atomicAdd(&out[76], 0.5*win[0]*fin[76]*subvol); 
  atomicAdd(&out[77], 0.5*win[0]*fin[77]*subvol); 
  atomicAdd(&out[114], 0.5*win[0]*fin[114]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[2] += 0.5*win[0]*fin[2]*subvol; 
  out[3] += 0.5*win[0]*fin[3]*subvol; 
  out[4] += 0.5*win[0]*fin[4]*subvol; 
  out[7] += 0.5*win[0]*fin[7]*subvol; 
  out[8] += 0.5*win[0]*fin[8]*subvol; 
  out[9] += 0.5*win[0]*fin[9]*subvol; 
  out[10] += 0.5*win[0]*fin[10]*subvol; 
  out[11] += 0.5*win[0]*fin[11]*subvol; 
  out[12] += 0.5*win[0]*fin[12]*subvol; 
  out[22] += 0.5*win[0]*fin[22]*subvol; 
  out[25] += 0.5*win[0]*fin[25]*subvol; 
  out[26] += 0.5*win[0]*fin[26]*subvol; 
  out[27] += 0.5*win[0]*fin[27]*subvol; 
  out[28] += 0.5*win[0]*fin[28]*subvol; 
  out[45] += 0.5*win[0]*fin[45]*subvol; 
  out[46] += 0.5*win[0]*fin[46]*subvol; 
  out[47] += 0.5*win[0]*fin[47]*subvol; 
  out[60] += 0.5*win[0]*fin[60]*subvol; 
  out[75] += 0.5*win[0]*fin[75]*subvol; 
  out[76] += 0.5*win[0]*fin[76]*subvol; 
  out[77] += 0.5*win[0]*fin[77]*subvol; 
  out[114] += 0.5*win[0]*fin[114]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxvyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[4]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[9]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[11]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[12]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[22]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[28]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[46]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[47]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[77]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[2]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[3]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[4]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[9]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[11]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[12]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[22]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[28]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[46]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[47]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[77]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxvyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], win[0]*fin[4]*subvol); 
  atomicAdd(&out[9], win[0]*fin[9]*subvol); 
  atomicAdd(&out[11], win[0]*fin[11]*subvol); 
  atomicAdd(&out[12], win[0]*fin[12]*subvol); 
  atomicAdd(&out[22], win[0]*fin[22]*subvol); 
  atomicAdd(&out[28], win[0]*fin[28]*subvol); 
  atomicAdd(&out[46], win[0]*fin[46]*subvol); 
  atomicAdd(&out[47], win[0]*fin[47]*subvol); 
  atomicAdd(&out[77], win[0]*fin[77]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
  out[2] += win[0]*fin[2]*subvol; 
  out[3] += win[0]*fin[3]*subvol; 
  out[4] += win[0]*fin[4]*subvol; 
  out[9] += win[0]*fin[9]*subvol; 
  out[11] += win[0]*fin[11]*subvol; 
  out[12] += win[0]*fin[12]*subvol; 
  out[22] += win[0]*fin[22]*subvol; 
  out[28] += win[0]*fin[28]*subvol; 
  out[46] += win[0]*fin[46]*subvol; 
  out[47] += win[0]*fin[47]*subvol; 
  out[77] += win[0]*fin[77]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgyvyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[4]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[8]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[10]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[12]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[22]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[27]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[45]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[47]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[76]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[1]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[3]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[4]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[8]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[10]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[12]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[22]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[27]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[45]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[47]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[76]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgyvyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], win[0]*fin[1]*subvol); 
  atomicAdd(&out[3], win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], win[0]*fin[4]*subvol); 
  atomicAdd(&out[8], win[0]*fin[8]*subvol); 
  atomicAdd(&out[10], win[0]*fin[10]*subvol); 
  atomicAdd(&out[12], win[0]*fin[12]*subvol); 
  atomicAdd(&out[22], win[0]*fin[22]*subvol); 
  atomicAdd(&out[27], win[0]*fin[27]*subvol); 
  atomicAdd(&out[45], win[0]*fin[45]*subvol); 
  atomicAdd(&out[47], win[0]*fin[47]*subvol); 
  atomicAdd(&out[76], win[0]*fin[76]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
  out[1] += win[0]*fin[1]*subvol; 
  out[3] += win[0]*fin[3]*subvol; 
  out[4] += win[0]*fin[4]*subvol; 
  out[8] += win[0]*fin[8]*subvol; 
  out[10] += win[0]*fin[10]*subvol; 
  out[12] += win[0]*fin[12]*subvol; 
  out[22] += win[0]*fin[22]*subvol; 
  out[27] += win[0]*fin[27]*subvol; 
  out[45] += win[0]*fin[45]*subvol; 
  out[47] += win[0]*fin[47]*subvol; 
  out[76] += win[0]*fin[76]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxyvyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[12]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[22]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[47]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[3]*subvol; 
  out[2] += 0.5*win[0]*fin[4]*subvol; 
  out[3] += 0.5*win[0]*fin[12]*subvol; 
  out[4] += 0.5*win[0]*fin[22]*subvol; 
  out[5] += 0.5*win[0]*fin[47]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxyvyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.0*fin[0]*win[0]*subvol); 
  atomicAdd(&out[3], 2.0*win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], 2.0*win[0]*fin[4]*subvol); 
  atomicAdd(&out[12], 2.0*win[0]*fin[12]*subvol); 
  atomicAdd(&out[22], 2.0*win[0]*fin[22]*subvol); 
  atomicAdd(&out[47], 2.0*win[0]*fin[47]*subvol); 
#else
  out[0] += 2.0*fin[0]*win[0]*subvol; 
  out[3] += 2.0*win[0]*fin[3]*subvol; 
  out[4] += 2.0*win[0]*fin[4]*subvol; 
  out[12] += 2.0*win[0]*fin[12]*subvol; 
  out[22] += 2.0*win[0]*fin[22]*subvol; 
  out[47] += 2.0*win[0]*fin[47]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgzvyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[4]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[7]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[10]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[11]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[22]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[26]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[45]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[46]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[75]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[1]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[2]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[4]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[7]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[10]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[11]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[22]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[26]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[45]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[46]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[75]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgzvyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], win[0]*fin[2]*subvol); 
  atomicAdd(&out[4], win[0]*fin[4]*subvol); 
  atomicAdd(&out[7], win[0]*fin[7]*subvol); 
  atomicAdd(&out[10], win[0]*fin[10]*subvol); 
  atomicAdd(&out[11], win[0]*fin[11]*subvol); 
  atomicAdd(&out[22], win[0]*fin[22]*subvol); 
  atomicAdd(&out[26], win[0]*fin[26]*subvol); 
  atomicAdd(&out[45], win[0]*fin[45]*subvol); 
  atomicAdd(&out[46], win[0]*fin[46]*subvol); 
  atomicAdd(&out[75], win[0]*fin[75]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
  out[1] += win[0]*fin[1]*subvol; 
  out[2] += win[0]*fin[2]*subvol; 
  out[4] += win[0]*fin[4]*subvol; 
  out[7] += win[0]*fin[7]*subvol; 
  out[10] += win[0]*fin[10]*subvol; 
  out[11] += win[0]*fin[11]*subvol; 
  out[22] += win[0]*fin[22]*subvol; 
  out[26] += win[0]*fin[26]*subvol; 
  out[45] += win[0]*fin[45]*subvol; 
  out[46] += win[0]*fin[46]*subvol; 
  out[75] += win[0]*fin[75]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxzvyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[11]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[22]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[46]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[2]*subvol; 
  out[2] += 0.5*win[0]*fin[4]*subvol; 
  out[3] += 0.5*win[0]*fin[11]*subvol; 
  out[4] += 0.5*win[0]*fin[22]*subvol; 
  out[5] += 0.5*win[0]*fin[46]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxzvyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.0*fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], 2.0*win[0]*fin[2]*subvol); 
  atomicAdd(&out[4], 2.0*win[0]*fin[4]*subvol); 
  atomicAdd(&out[11], 2.0*win[0]*fin[11]*subvol); 
  atomicAdd(&out[22], 2.0*win[0]*fin[22]*subvol); 
  atomicAdd(&out[46], 2.0*win[0]*fin[46]*subvol); 
#else
  out[0] += 2.0*fin[0]*win[0]*subvol; 
  out[2] += 2.0*win[0]*fin[2]*subvol; 
  out[4] += 2.0*win[0]*fin[4]*subvol; 
  out[11] += 2.0*win[0]*fin[11]*subvol; 
  out[22] += 2.0*win[0]*fin[22]*subvol; 
  out[46] += 2.0*win[0]*fin[46]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgyzvyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[10]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[22]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[45]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[2] += 0.5*win[0]*fin[4]*subvol; 
  out[3] += 0.5*win[0]*fin[10]*subvol; 
  out[4] += 0.5*win[0]*fin[22]*subvol; 
  out[5] += 0.5*win[0]*fin[45]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgyzvyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.0*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 2.0*win[0]*fin[1]*subvol); 
  atomicAdd(&out[4], 2.0*win[0]*fin[4]*subvol); 
  atomicAdd(&out[10], 2.0*win[0]*fin[10]*subvol); 
  atomicAdd(&out[22], 2.0*win[0]*fin[22]*subvol); 
  atomicAdd(&out[45], 2.0*win[0]*fin[45]*subvol); 
#else
  out[0] += 2.0*fin[0]*win[0]*subvol; 
  out[1] += 2.0*win[0]*fin[1]*subvol; 
  out[4] += 2.0*win[0]*fin[4]*subvol; 
  out[10] += 2.0*win[0]*fin[10]*subvol; 
  out[22] += 2.0*win[0]*fin[22]*subvol; 
  out[45] += 2.0*win[0]*fin[45]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxyzvyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[4]*subvol); 
  atomicAdd(&out[2], 0.7071067811865475*win[0]*fin[22]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[4]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[22]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxyzvyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 4.0*fin[0]*win[0]*subvol); 
  atomicAdd(&out[4], 4.0*win[0]*fin[4]*subvol); 
  atomicAdd(&out[22], 4.0*win[0]*fin[22]*subvol); 
#else
  out[0] += 4.0*fin[0]*win[0]*subvol; 
  out[4] += 4.0*win[0]*fin[4]*subvol; 
  out[22] += 4.0*win[0]*fin[22]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgvxvyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[7]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[8]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[9]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[25]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[1]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[2]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[3]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[7]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[8]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[9]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[25]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgvxvyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], win[0]*fin[3]*subvol); 
  atomicAdd(&out[7], win[0]*fin[7]*subvol); 
  atomicAdd(&out[8], win[0]*fin[8]*subvol); 
  atomicAdd(&out[9], win[0]*fin[9]*subvol); 
  atomicAdd(&out[25], win[0]*fin[25]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
  out[1] += win[0]*fin[1]*subvol; 
  out[2] += win[0]*fin[2]*subvol; 
  out[3] += win[0]*fin[3]*subvol; 
  out[7] += win[0]*fin[7]*subvol; 
  out[8] += win[0]*fin[8]*subvol; 
  out[9] += win[0]*fin[9]*subvol; 
  out[25] += win[0]*fin[25]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxvxvyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[9]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[2]*subvol; 
  out[2] += 0.5*win[0]*fin[3]*subvol; 
  out[3] += 0.5*win[0]*fin[9]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxvxvyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.0*fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], 2.0*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 2.0*win[0]*fin[3]*subvol); 
  atomicAdd(&out[9], 2.0*win[0]*fin[9]*subvol); 
#else
  out[0] += 2.0*fin[0]*win[0]*subvol; 
  out[2] += 2.0*win[0]*fin[2]*subvol; 
  out[3] += 2.0*win[0]*fin[3]*subvol; 
  out[9] += 2.0*win[0]*fin[9]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgyvxvyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[8]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[2] += 0.5*win[0]*fin[3]*subvol; 
  out[3] += 0.5*win[0]*fin[8]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgyvxvyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.0*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 2.0*win[0]*fin[1]*subvol); 
  atomicAdd(&out[3], 2.0*win[0]*fin[3]*subvol); 
  atomicAdd(&out[8], 2.0*win[0]*fin[8]*subvol); 
#else
  out[0] += 2.0*fin[0]*win[0]*subvol; 
  out[1] += 2.0*win[0]*fin[1]*subvol; 
  out[3] += 2.0*win[0]*fin[3]*subvol; 
  out[8] += 2.0*win[0]*fin[8]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxyvxvyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[3]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[3]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxyvxvyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 4.0*fin[0]*win[0]*subvol); 
  atomicAdd(&out[3], 4.0*win[0]*fin[3]*subvol); 
#else
  out[0] += 4.0*fin[0]*win[0]*subvol; 
  out[3] += 4.0*win[0]*fin[3]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgzvxvyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[7]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[2] += 0.5*win[0]*fin[2]*subvol; 
  out[3] += 0.5*win[0]*fin[7]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgzvxvyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.0*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 2.0*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 2.0*win[0]*fin[2]*subvol); 
  atomicAdd(&out[7], 2.0*win[0]*fin[7]*subvol); 
#else
  out[0] += 2.0*fin[0]*win[0]*subvol; 
  out[1] += 2.0*win[0]*fin[1]*subvol; 
  out[2] += 2.0*win[0]*fin[2]*subvol; 
  out[7] += 2.0*win[0]*fin[7]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxzvxvyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[2]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[2]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxzvxvyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 4.0*fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], 4.0*win[0]*fin[2]*subvol); 
#else
  out[0] += 4.0*fin[0]*win[0]*subvol; 
  out[2] += 4.0*win[0]*fin[2]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgyzvxvyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[1]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[1]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgyzvxvyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 4.0*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 4.0*win[0]*fin[1]*subvol); 
#else
  out[0] += 4.0*fin[0]*win[0]*subvol; 
  out[1] += 4.0*win[0]*fin[1]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxyzvxvyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x3v_hyb_p1_avgxyzvxvyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 8.0*fin[0]*win[0]*subvol); 
#else
  out[0] += 8.0*fin[0]*win[0]*subvol; 
#endif
} 

