#include <gkyl_fpo_vlasov_kernels.h> 

GKYL_CU_DH void fpo_sgn_drag_coeff_2x3v_vy_tensor_p1_upvy(const double *drag_coeff_surf, double *sgn_drag_coeff_surf, int *const_sgn_drag_coeff_surf) {
  // drag_coeff_surf: Surface projection of drag coefficient at lower boundary.
  // sgn_drag_coeff_surf: sign(drag_coeff_surf) evaluated at quadrature points along lower surface.
  // const_sgn_drag_coeff_surf: 1 if sign(drag_coeff_surf) is constant along lower boundary, 0 otherwise. 

  const double *alpha_surf = &drag_coeff_surf[36]; 
  double *sgn_alpha_surf = &sgn_drag_coeff_surf[36]; 

  int const_sgn_alpha_surf = 1;  
  
  if (0.2*(alpha_surf[35]+alpha_surf[30])-0.19999999999999954*(alpha_surf[34]+alpha_surf[33])-0.2999999999999998*alpha_surf[32]-0.29999999999999993*(alpha_surf[31]+alpha_surf[19])+0.3*(alpha_surf[29]+alpha_surf[28]+alpha_surf[26]+alpha_surf[25])+0.22360679774997896*(alpha_surf[27]+alpha_surf[24]+alpha_surf[12]+alpha_surf[11])+0.45*(alpha_surf[23]+alpha_surf[10])-0.29999999999999977*alpha_surf[22]-0.22360679774997902*(alpha_surf[21]+alpha_surf[20]+alpha_surf[18]+alpha_surf[17])-0.45*(alpha_surf[16]+alpha_surf[15])-0.33541019662496846*(alpha_surf[14]+alpha_surf[13]+alpha_surf[4]+alpha_surf[3])+0.33541019662496846*(alpha_surf[9]+alpha_surf[8]+alpha_surf[7]+alpha_surf[6])+0.25*(alpha_surf[5]+alpha_surf[0])-0.25*(alpha_surf[2]+alpha_surf[1]) > 0.) 
    sgn_alpha_surf[0] = 1.0; 
  else  
    sgn_alpha_surf[0] = -1.0; 
  
  if (-(0.25*(alpha_surf[35]+alpha_surf[30]+alpha_surf[2]+alpha_surf[1]))+0.24999999999999994*(alpha_surf[34]+alpha_surf[33])+0.375*(alpha_surf[32]+alpha_surf[22])-0.375*(alpha_surf[29]+alpha_surf[28])-0.2795084971874737*(alpha_surf[27]+alpha_surf[12])+0.22360679774997896*(alpha_surf[24]+alpha_surf[11])+0.2795084971874738*(alpha_surf[21]+alpha_surf[20])-0.22360679774997902*(alpha_surf[18]+alpha_surf[17])-0.33541019662496846*(alpha_surf[13]+alpha_surf[3])+0.33541019662496846*(alpha_surf[7]+alpha_surf[6])+0.25*(alpha_surf[5]+alpha_surf[0]) > 0.) 
    sgn_alpha_surf[1] = 1.0; 
  else  
    sgn_alpha_surf[1] = -1.0; 
  
  if (sgn_alpha_surf[1] == sgn_alpha_surf[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2*(alpha_surf[35]+alpha_surf[30])-0.19999999999999954*(alpha_surf[34]+alpha_surf[33])-0.2999999999999998*alpha_surf[32]+0.29999999999999993*(alpha_surf[31]+alpha_surf[19])+0.3*(alpha_surf[29]+alpha_surf[28])+0.22360679774997896*(alpha_surf[27]+alpha_surf[24]+alpha_surf[12]+alpha_surf[11])-0.3*(alpha_surf[26]+alpha_surf[25])-0.45*(alpha_surf[23]+alpha_surf[10])-0.29999999999999977*alpha_surf[22]-0.22360679774997902*(alpha_surf[21]+alpha_surf[20]+alpha_surf[18]+alpha_surf[17])+0.45*(alpha_surf[16]+alpha_surf[15])+0.33541019662496846*(alpha_surf[14]+alpha_surf[7]+alpha_surf[6]+alpha_surf[4])-0.33541019662496846*(alpha_surf[13]+alpha_surf[9]+alpha_surf[8]+alpha_surf[3])+0.25*(alpha_surf[5]+alpha_surf[0])-0.25*(alpha_surf[2]+alpha_surf[1]) > 0.) 
    sgn_alpha_surf[2] = 1.0; 
  else  
    sgn_alpha_surf[2] = -1.0; 
  
  if (sgn_alpha_surf[2] == sgn_alpha_surf[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.25*(alpha_surf[35]+alpha_surf[30]+alpha_surf[2]+alpha_surf[1]))+0.25*(alpha_surf[34]+alpha_surf[33]+alpha_surf[5]+alpha_surf[0])+0.375*(alpha_surf[31]+alpha_surf[19])+0.22360679774997896*(alpha_surf[27]+alpha_surf[12])-0.375*(alpha_surf[26]+alpha_surf[25])-0.2795084971874737*(alpha_surf[24]+alpha_surf[11])-0.22360679774997902*(alpha_surf[21]+alpha_surf[20])+0.2795084971874738*(alpha_surf[18]+alpha_surf[17])-0.33541019662496846*(alpha_surf[14]+alpha_surf[4])+0.33541019662496846*(alpha_surf[9]+alpha_surf[8]) > 0.) 
    sgn_alpha_surf[3] = 1.0; 
  else  
    sgn_alpha_surf[3] = -1.0; 
  
  if (sgn_alpha_surf[3] == sgn_alpha_surf[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3125*(alpha_surf[35]+alpha_surf[30])-0.3125*(alpha_surf[34]+alpha_surf[33])-0.2795084971874737*(alpha_surf[27]+alpha_surf[24]+alpha_surf[12]+alpha_surf[11])+0.2795084971874738*(alpha_surf[21]+alpha_surf[20]+alpha_surf[18]+alpha_surf[17])+0.25*(alpha_surf[5]+alpha_surf[0])-0.25*(alpha_surf[2]+alpha_surf[1]) > 0.) 
    sgn_alpha_surf[4] = 1.0; 
  else  
    sgn_alpha_surf[4] = -1.0; 
  
  if (sgn_alpha_surf[4] == sgn_alpha_surf[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.25*(alpha_surf[35]+alpha_surf[30]+alpha_surf[2]+alpha_surf[1]))+0.25*(alpha_surf[34]+alpha_surf[33]+alpha_surf[5]+alpha_surf[0])-0.375*(alpha_surf[31]+alpha_surf[19])+0.22360679774997896*(alpha_surf[27]+alpha_surf[12])+0.375*(alpha_surf[26]+alpha_surf[25])-0.2795084971874737*(alpha_surf[24]+alpha_surf[11])-0.22360679774997902*(alpha_surf[21]+alpha_surf[20])+0.2795084971874738*(alpha_surf[18]+alpha_surf[17])+0.33541019662496846*(alpha_surf[14]+alpha_surf[4])-0.33541019662496846*(alpha_surf[9]+alpha_surf[8]) > 0.) 
    sgn_alpha_surf[5] = 1.0; 
  else  
    sgn_alpha_surf[5] = -1.0; 
  
  if (sgn_alpha_surf[5] == sgn_alpha_surf[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2*(alpha_surf[35]+alpha_surf[30])-0.19999999999999954*(alpha_surf[34]+alpha_surf[33])+0.2999999999999998*alpha_surf[32]-0.29999999999999993*(alpha_surf[31]+alpha_surf[19])-0.3*(alpha_surf[29]+alpha_surf[28])+0.22360679774997896*(alpha_surf[27]+alpha_surf[24]+alpha_surf[12]+alpha_surf[11])+0.3*(alpha_surf[26]+alpha_surf[25])-0.45*(alpha_surf[23]+alpha_surf[10])+0.29999999999999977*alpha_surf[22]-0.22360679774997902*(alpha_surf[21]+alpha_surf[20]+alpha_surf[18]+alpha_surf[17])+0.45*(alpha_surf[16]+alpha_surf[15])-0.33541019662496846*(alpha_surf[14]+alpha_surf[7]+alpha_surf[6]+alpha_surf[4])+0.33541019662496846*(alpha_surf[13]+alpha_surf[9]+alpha_surf[8]+alpha_surf[3])+0.25*(alpha_surf[5]+alpha_surf[0])-0.25*(alpha_surf[2]+alpha_surf[1]) > 0.) 
    sgn_alpha_surf[6] = 1.0; 
  else  
    sgn_alpha_surf[6] = -1.0; 
  
  if (sgn_alpha_surf[6] == sgn_alpha_surf[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.25*(alpha_surf[35]+alpha_surf[30]+alpha_surf[2]+alpha_surf[1]))+0.24999999999999994*(alpha_surf[34]+alpha_surf[33])-0.375*(alpha_surf[32]+alpha_surf[22])+0.375*(alpha_surf[29]+alpha_surf[28])-0.2795084971874737*(alpha_surf[27]+alpha_surf[12])+0.22360679774997896*(alpha_surf[24]+alpha_surf[11])+0.2795084971874738*(alpha_surf[21]+alpha_surf[20])-0.22360679774997902*(alpha_surf[18]+alpha_surf[17])+0.33541019662496846*(alpha_surf[13]+alpha_surf[3])-0.33541019662496846*(alpha_surf[7]+alpha_surf[6])+0.25*(alpha_surf[5]+alpha_surf[0]) > 0.) 
    sgn_alpha_surf[7] = 1.0; 
  else  
    sgn_alpha_surf[7] = -1.0; 
  
  if (sgn_alpha_surf[7] == sgn_alpha_surf[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2*(alpha_surf[35]+alpha_surf[30])-0.19999999999999954*(alpha_surf[34]+alpha_surf[33])+0.2999999999999998*alpha_surf[32]+0.29999999999999993*(alpha_surf[31]+alpha_surf[19])-0.3*(alpha_surf[29]+alpha_surf[28]+alpha_surf[26]+alpha_surf[25])+0.22360679774997896*(alpha_surf[27]+alpha_surf[24]+alpha_surf[12]+alpha_surf[11])+0.45*(alpha_surf[23]+alpha_surf[10])+0.29999999999999977*alpha_surf[22]-0.22360679774997902*(alpha_surf[21]+alpha_surf[20]+alpha_surf[18]+alpha_surf[17])-0.45*(alpha_surf[16]+alpha_surf[15])+0.33541019662496846*(alpha_surf[14]+alpha_surf[13]+alpha_surf[4]+alpha_surf[3])-0.33541019662496846*(alpha_surf[9]+alpha_surf[8]+alpha_surf[7]+alpha_surf[6])+0.25*(alpha_surf[5]+alpha_surf[0])-0.25*(alpha_surf[2]+alpha_surf[1]) > 0.) 
    sgn_alpha_surf[8] = 1.0; 
  else  
    sgn_alpha_surf[8] = -1.0; 
  
  if (sgn_alpha_surf[8] == sgn_alpha_surf[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.2*alpha_surf[35])+0.19999999999999954*alpha_surf[34]-0.19999999999999954*alpha_surf[33]+0.2999999999999998*alpha_surf[32]+0.29999999999999993*alpha_surf[31]+0.2*alpha_surf[30]-0.3*(alpha_surf[29]+alpha_surf[26])+0.3*(alpha_surf[28]+alpha_surf[25])-0.22360679774997896*(alpha_surf[27]+alpha_surf[24])-0.45*(alpha_surf[23]+alpha_surf[15])-0.29999999999999977*alpha_surf[22]+0.22360679774997902*(alpha_surf[21]+alpha_surf[18])-0.22360679774997902*(alpha_surf[20]+alpha_surf[17])-0.29999999999999993*alpha_surf[19]+0.45*(alpha_surf[16]+alpha_surf[10])+0.33541019662496846*(alpha_surf[14]+alpha_surf[13]+alpha_surf[8]+alpha_surf[6])+0.22360679774997896*(alpha_surf[12]+alpha_surf[11])-0.33541019662496846*(alpha_surf[9]+alpha_surf[7]+alpha_surf[4]+alpha_surf[3])-0.25*(alpha_surf[5]+alpha_surf[1])+0.25*(alpha_surf[2]+alpha_surf[0]) > 0.) 
    sgn_alpha_surf[9] = 1.0; 
  else  
    sgn_alpha_surf[9] = -1.0; 
  
  if (sgn_alpha_surf[9] == sgn_alpha_surf[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*(alpha_surf[35]+alpha_surf[2]+alpha_surf[0])-0.24999999999999994*alpha_surf[34]+0.24999999999999994*alpha_surf[33]-0.375*(alpha_surf[32]+alpha_surf[28])-0.25*(alpha_surf[30]+alpha_surf[5]+alpha_surf[1])+0.375*(alpha_surf[29]+alpha_surf[22])+0.2795084971874737*alpha_surf[27]-0.22360679774997896*alpha_surf[24]-0.2795084971874738*alpha_surf[21]+0.2795084971874738*alpha_surf[20]+0.22360679774997902*alpha_surf[18]-0.22360679774997902*alpha_surf[17]+0.33541019662496846*(alpha_surf[13]+alpha_surf[6])-0.2795084971874737*alpha_surf[12]+0.22360679774997896*alpha_surf[11]-0.33541019662496846*(alpha_surf[7]+alpha_surf[3]) > 0.) 
    sgn_alpha_surf[10] = 1.0; 
  else  
    sgn_alpha_surf[10] = -1.0; 
  
  if (sgn_alpha_surf[10] == sgn_alpha_surf[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.2*alpha_surf[35])+0.19999999999999954*alpha_surf[34]-0.19999999999999954*alpha_surf[33]+0.2999999999999998*alpha_surf[32]-0.29999999999999993*alpha_surf[31]+0.2*alpha_surf[30]-0.3*(alpha_surf[29]+alpha_surf[25])+0.3*(alpha_surf[28]+alpha_surf[26])-0.22360679774997896*(alpha_surf[27]+alpha_surf[24])+0.45*(alpha_surf[23]+alpha_surf[15])-0.29999999999999977*alpha_surf[22]+0.22360679774997902*(alpha_surf[21]+alpha_surf[18])-0.22360679774997902*(alpha_surf[20]+alpha_surf[17])+0.29999999999999993*alpha_surf[19]-0.45*(alpha_surf[16]+alpha_surf[10])-0.33541019662496846*(alpha_surf[14]+alpha_surf[8]+alpha_surf[7]+alpha_surf[3])+0.33541019662496846*(alpha_surf[13]+alpha_surf[9]+alpha_surf[6]+alpha_surf[4])+0.22360679774997896*(alpha_surf[12]+alpha_surf[11])-0.25*(alpha_surf[5]+alpha_surf[1])+0.25*(alpha_surf[2]+alpha_surf[0]) > 0.) 
    sgn_alpha_surf[11] = 1.0; 
  else  
    sgn_alpha_surf[11] = -1.0; 
  
  if (sgn_alpha_surf[11] == sgn_alpha_surf[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*(alpha_surf[35]+alpha_surf[33]+alpha_surf[2]+alpha_surf[0])-0.25*(alpha_surf[34]+alpha_surf[30]+alpha_surf[5]+alpha_surf[1])-0.375*(alpha_surf[31]+alpha_surf[25])-0.22360679774997896*alpha_surf[27]+0.375*(alpha_surf[26]+alpha_surf[19])+0.2795084971874737*alpha_surf[24]+0.22360679774997902*alpha_surf[21]-0.22360679774997902*alpha_surf[20]-0.2795084971874738*alpha_surf[18]+0.2795084971874738*alpha_surf[17]+0.33541019662496846*(alpha_surf[14]+alpha_surf[8])+0.22360679774997896*alpha_surf[12]-0.2795084971874737*alpha_surf[11]-0.33541019662496846*(alpha_surf[9]+alpha_surf[4]) > 0.) 
    sgn_alpha_surf[12] = 1.0; 
  else  
    sgn_alpha_surf[12] = -1.0; 
  
  if (sgn_alpha_surf[12] == sgn_alpha_surf[11]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.3125*(alpha_surf[35]+alpha_surf[33]))+0.3125*(alpha_surf[34]+alpha_surf[30])+0.2795084971874737*(alpha_surf[27]+alpha_surf[24])-0.2795084971874738*(alpha_surf[21]+alpha_surf[18])+0.2795084971874738*(alpha_surf[20]+alpha_surf[17])-0.2795084971874737*(alpha_surf[12]+alpha_surf[11])-0.25*(alpha_surf[5]+alpha_surf[1])+0.25*(alpha_surf[2]+alpha_surf[0]) > 0.) 
    sgn_alpha_surf[13] = 1.0; 
  else  
    sgn_alpha_surf[13] = -1.0; 
  
  if (sgn_alpha_surf[13] == sgn_alpha_surf[12]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*(alpha_surf[35]+alpha_surf[33]+alpha_surf[2]+alpha_surf[0])-0.25*(alpha_surf[34]+alpha_surf[30]+alpha_surf[5]+alpha_surf[1])+0.375*(alpha_surf[31]+alpha_surf[25])-0.22360679774997896*alpha_surf[27]-0.375*(alpha_surf[26]+alpha_surf[19])+0.2795084971874737*alpha_surf[24]+0.22360679774997902*alpha_surf[21]-0.22360679774997902*alpha_surf[20]-0.2795084971874738*alpha_surf[18]+0.2795084971874738*alpha_surf[17]-0.33541019662496846*(alpha_surf[14]+alpha_surf[8])+0.22360679774997896*alpha_surf[12]-0.2795084971874737*alpha_surf[11]+0.33541019662496846*(alpha_surf[9]+alpha_surf[4]) > 0.) 
    sgn_alpha_surf[14] = 1.0; 
  else  
    sgn_alpha_surf[14] = -1.0; 
  
  if (sgn_alpha_surf[14] == sgn_alpha_surf[13]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.2*alpha_surf[35])+0.19999999999999954*alpha_surf[34]-0.19999999999999954*alpha_surf[33]-0.2999999999999998*alpha_surf[32]+0.29999999999999993*alpha_surf[31]+0.2*alpha_surf[30]+0.3*(alpha_surf[29]+alpha_surf[25])-0.3*(alpha_surf[28]+alpha_surf[26])-0.22360679774997896*(alpha_surf[27]+alpha_surf[24])+0.45*(alpha_surf[23]+alpha_surf[15])+0.29999999999999977*alpha_surf[22]+0.22360679774997902*(alpha_surf[21]+alpha_surf[18])-0.22360679774997902*(alpha_surf[20]+alpha_surf[17])-0.29999999999999993*alpha_surf[19]-0.45*(alpha_surf[16]+alpha_surf[10])+0.33541019662496846*(alpha_surf[14]+alpha_surf[8]+alpha_surf[7]+alpha_surf[3])-0.33541019662496846*(alpha_surf[13]+alpha_surf[9]+alpha_surf[6]+alpha_surf[4])+0.22360679774997896*(alpha_surf[12]+alpha_surf[11])-0.25*(alpha_surf[5]+alpha_surf[1])+0.25*(alpha_surf[2]+alpha_surf[0]) > 0.) 
    sgn_alpha_surf[15] = 1.0; 
  else  
    sgn_alpha_surf[15] = -1.0; 
  
  if (sgn_alpha_surf[15] == sgn_alpha_surf[14]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*(alpha_surf[35]+alpha_surf[2]+alpha_surf[0])-0.24999999999999994*alpha_surf[34]+0.24999999999999994*alpha_surf[33]+0.375*(alpha_surf[32]+alpha_surf[28])-0.25*(alpha_surf[30]+alpha_surf[5]+alpha_surf[1])-0.375*(alpha_surf[29]+alpha_surf[22])+0.2795084971874737*alpha_surf[27]-0.22360679774997896*alpha_surf[24]-0.2795084971874738*alpha_surf[21]+0.2795084971874738*alpha_surf[20]+0.22360679774997902*alpha_surf[18]-0.22360679774997902*alpha_surf[17]-0.33541019662496846*(alpha_surf[13]+alpha_surf[6])-0.2795084971874737*alpha_surf[12]+0.22360679774997896*alpha_surf[11]+0.33541019662496846*(alpha_surf[7]+alpha_surf[3]) > 0.) 
    sgn_alpha_surf[16] = 1.0; 
  else  
    sgn_alpha_surf[16] = -1.0; 
  
  if (sgn_alpha_surf[16] == sgn_alpha_surf[15]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.2*alpha_surf[35])+0.19999999999999954*alpha_surf[34]-0.19999999999999954*alpha_surf[33]-0.2999999999999998*alpha_surf[32]-0.29999999999999993*alpha_surf[31]+0.2*alpha_surf[30]+0.3*(alpha_surf[29]+alpha_surf[26])-0.3*(alpha_surf[28]+alpha_surf[25])-0.22360679774997896*(alpha_surf[27]+alpha_surf[24])-0.45*(alpha_surf[23]+alpha_surf[15])+0.29999999999999977*alpha_surf[22]+0.22360679774997902*(alpha_surf[21]+alpha_surf[18])-0.22360679774997902*(alpha_surf[20]+alpha_surf[17])+0.29999999999999993*alpha_surf[19]+0.45*(alpha_surf[16]+alpha_surf[10])-0.33541019662496846*(alpha_surf[14]+alpha_surf[13]+alpha_surf[8]+alpha_surf[6])+0.22360679774997896*(alpha_surf[12]+alpha_surf[11])+0.33541019662496846*(alpha_surf[9]+alpha_surf[7]+alpha_surf[4]+alpha_surf[3])-0.25*(alpha_surf[5]+alpha_surf[1])+0.25*(alpha_surf[2]+alpha_surf[0]) > 0.) 
    sgn_alpha_surf[17] = 1.0; 
  else  
    sgn_alpha_surf[17] = -1.0; 
  
  if (sgn_alpha_surf[17] == sgn_alpha_surf[16]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.2*alpha_surf[35])-0.19999999999999954*alpha_surf[34]+0.19999999999999954*alpha_surf[33]+0.2999999999999998*alpha_surf[32]+0.29999999999999993*alpha_surf[31]+0.2*alpha_surf[30]+0.3*(alpha_surf[29]+alpha_surf[26])-0.3*(alpha_surf[28]+alpha_surf[25])-0.22360679774997896*(alpha_surf[27]+alpha_surf[24])-0.45*(alpha_surf[23]+alpha_surf[16])-0.29999999999999977*alpha_surf[22]-0.22360679774997902*(alpha_surf[21]+alpha_surf[18])+0.22360679774997902*(alpha_surf[20]+alpha_surf[17])-0.29999999999999993*alpha_surf[19]+0.45*(alpha_surf[15]+alpha_surf[10])+0.33541019662496846*(alpha_surf[14]+alpha_surf[13]+alpha_surf[9]+alpha_surf[7])+0.22360679774997896*(alpha_surf[12]+alpha_surf[11])-0.33541019662496846*(alpha_surf[8]+alpha_surf[6]+alpha_surf[4]+alpha_surf[3])-0.25*(alpha_surf[5]+alpha_surf[2])+0.25*(alpha_surf[1]+alpha_surf[0]) > 0.) 
    sgn_alpha_surf[18] = 1.0; 
  else  
    sgn_alpha_surf[18] = -1.0; 
  
  if (sgn_alpha_surf[18] == sgn_alpha_surf[17]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*(alpha_surf[35]+alpha_surf[1]+alpha_surf[0])+0.24999999999999994*alpha_surf[34]-0.24999999999999994*alpha_surf[33]-0.375*(alpha_surf[32]+alpha_surf[29])-0.25*(alpha_surf[30]+alpha_surf[5]+alpha_surf[2])+0.375*(alpha_surf[28]+alpha_surf[22])+0.2795084971874737*alpha_surf[27]-0.22360679774997896*alpha_surf[24]+0.2795084971874738*alpha_surf[21]-0.2795084971874738*alpha_surf[20]-0.22360679774997902*alpha_surf[18]+0.22360679774997902*alpha_surf[17]+0.33541019662496846*(alpha_surf[13]+alpha_surf[7])-0.2795084971874737*alpha_surf[12]+0.22360679774997896*alpha_surf[11]-0.33541019662496846*(alpha_surf[6]+alpha_surf[3]) > 0.) 
    sgn_alpha_surf[19] = 1.0; 
  else  
    sgn_alpha_surf[19] = -1.0; 
  
  if (sgn_alpha_surf[19] == sgn_alpha_surf[18]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.2*alpha_surf[35])-0.19999999999999954*alpha_surf[34]+0.19999999999999954*alpha_surf[33]+0.2999999999999998*alpha_surf[32]-0.29999999999999993*alpha_surf[31]+0.2*alpha_surf[30]+0.3*(alpha_surf[29]+alpha_surf[25])-0.3*(alpha_surf[28]+alpha_surf[26])-0.22360679774997896*(alpha_surf[27]+alpha_surf[24])+0.45*(alpha_surf[23]+alpha_surf[16])-0.29999999999999977*alpha_surf[22]-0.22360679774997902*(alpha_surf[21]+alpha_surf[18])+0.22360679774997902*(alpha_surf[20]+alpha_surf[17])+0.29999999999999993*alpha_surf[19]-0.45*(alpha_surf[15]+alpha_surf[10])-0.33541019662496846*(alpha_surf[14]+alpha_surf[9]+alpha_surf[6]+alpha_surf[3])+0.33541019662496846*(alpha_surf[13]+alpha_surf[8]+alpha_surf[7]+alpha_surf[4])+0.22360679774997896*(alpha_surf[12]+alpha_surf[11])-0.25*(alpha_surf[5]+alpha_surf[2])+0.25*(alpha_surf[1]+alpha_surf[0]) > 0.) 
    sgn_alpha_surf[20] = 1.0; 
  else  
    sgn_alpha_surf[20] = -1.0; 
  
  if (sgn_alpha_surf[20] == sgn_alpha_surf[19]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*(alpha_surf[35]+alpha_surf[34]+alpha_surf[1]+alpha_surf[0])-0.25*(alpha_surf[33]+alpha_surf[30]+alpha_surf[5]+alpha_surf[2])-0.375*(alpha_surf[31]+alpha_surf[26])-0.22360679774997896*alpha_surf[27]+0.375*(alpha_surf[25]+alpha_surf[19])+0.2795084971874737*alpha_surf[24]-0.22360679774997902*alpha_surf[21]+0.22360679774997902*alpha_surf[20]+0.2795084971874738*alpha_surf[18]-0.2795084971874738*alpha_surf[17]+0.33541019662496846*(alpha_surf[14]+alpha_surf[9])+0.22360679774997896*alpha_surf[12]-0.2795084971874737*alpha_surf[11]-0.33541019662496846*(alpha_surf[8]+alpha_surf[4]) > 0.) 
    sgn_alpha_surf[21] = 1.0; 
  else  
    sgn_alpha_surf[21] = -1.0; 
  
  if (sgn_alpha_surf[21] == sgn_alpha_surf[20]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.3125*(alpha_surf[35]+alpha_surf[34]))+0.3125*(alpha_surf[33]+alpha_surf[30])+0.2795084971874737*(alpha_surf[27]+alpha_surf[24])+0.2795084971874738*(alpha_surf[21]+alpha_surf[18])-0.2795084971874738*(alpha_surf[20]+alpha_surf[17])-0.2795084971874737*(alpha_surf[12]+alpha_surf[11])-0.25*(alpha_surf[5]+alpha_surf[2])+0.25*(alpha_surf[1]+alpha_surf[0]) > 0.) 
    sgn_alpha_surf[22] = 1.0; 
  else  
    sgn_alpha_surf[22] = -1.0; 
  
  if (sgn_alpha_surf[22] == sgn_alpha_surf[21]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*(alpha_surf[35]+alpha_surf[34]+alpha_surf[1]+alpha_surf[0])-0.25*(alpha_surf[33]+alpha_surf[30]+alpha_surf[5]+alpha_surf[2])+0.375*(alpha_surf[31]+alpha_surf[26])-0.22360679774997896*alpha_surf[27]-0.375*(alpha_surf[25]+alpha_surf[19])+0.2795084971874737*alpha_surf[24]-0.22360679774997902*alpha_surf[21]+0.22360679774997902*alpha_surf[20]+0.2795084971874738*alpha_surf[18]-0.2795084971874738*alpha_surf[17]-0.33541019662496846*(alpha_surf[14]+alpha_surf[9])+0.22360679774997896*alpha_surf[12]-0.2795084971874737*alpha_surf[11]+0.33541019662496846*(alpha_surf[8]+alpha_surf[4]) > 0.) 
    sgn_alpha_surf[23] = 1.0; 
  else  
    sgn_alpha_surf[23] = -1.0; 
  
  if (sgn_alpha_surf[23] == sgn_alpha_surf[22]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.2*alpha_surf[35])-0.19999999999999954*alpha_surf[34]+0.19999999999999954*alpha_surf[33]-0.2999999999999998*alpha_surf[32]+0.29999999999999993*alpha_surf[31]+0.2*alpha_surf[30]-0.3*(alpha_surf[29]+alpha_surf[25])+0.3*(alpha_surf[28]+alpha_surf[26])-0.22360679774997896*(alpha_surf[27]+alpha_surf[24])+0.45*(alpha_surf[23]+alpha_surf[16])+0.29999999999999977*alpha_surf[22]-0.22360679774997902*(alpha_surf[21]+alpha_surf[18])+0.22360679774997902*(alpha_surf[20]+alpha_surf[17])-0.29999999999999993*alpha_surf[19]-0.45*(alpha_surf[15]+alpha_surf[10])+0.33541019662496846*(alpha_surf[14]+alpha_surf[9]+alpha_surf[6]+alpha_surf[3])-0.33541019662496846*(alpha_surf[13]+alpha_surf[8]+alpha_surf[7]+alpha_surf[4])+0.22360679774997896*(alpha_surf[12]+alpha_surf[11])-0.25*(alpha_surf[5]+alpha_surf[2])+0.25*(alpha_surf[1]+alpha_surf[0]) > 0.) 
    sgn_alpha_surf[24] = 1.0; 
  else  
    sgn_alpha_surf[24] = -1.0; 
  
  if (sgn_alpha_surf[24] == sgn_alpha_surf[23]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*(alpha_surf[35]+alpha_surf[1]+alpha_surf[0])+0.24999999999999994*alpha_surf[34]-0.24999999999999994*alpha_surf[33]+0.375*(alpha_surf[32]+alpha_surf[29])-0.25*(alpha_surf[30]+alpha_surf[5]+alpha_surf[2])-0.375*(alpha_surf[28]+alpha_surf[22])+0.2795084971874737*alpha_surf[27]-0.22360679774997896*alpha_surf[24]+0.2795084971874738*alpha_surf[21]-0.2795084971874738*alpha_surf[20]-0.22360679774997902*alpha_surf[18]+0.22360679774997902*alpha_surf[17]-0.33541019662496846*(alpha_surf[13]+alpha_surf[7])-0.2795084971874737*alpha_surf[12]+0.22360679774997896*alpha_surf[11]+0.33541019662496846*(alpha_surf[6]+alpha_surf[3]) > 0.) 
    sgn_alpha_surf[25] = 1.0; 
  else  
    sgn_alpha_surf[25] = -1.0; 
  
  if (sgn_alpha_surf[25] == sgn_alpha_surf[24]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.2*alpha_surf[35])-0.19999999999999954*alpha_surf[34]+0.19999999999999954*alpha_surf[33]-0.2999999999999998*alpha_surf[32]-0.29999999999999993*alpha_surf[31]+0.2*alpha_surf[30]-0.3*(alpha_surf[29]+alpha_surf[26])+0.3*(alpha_surf[28]+alpha_surf[25])-0.22360679774997896*(alpha_surf[27]+alpha_surf[24])-0.45*(alpha_surf[23]+alpha_surf[16])+0.29999999999999977*alpha_surf[22]-0.22360679774997902*(alpha_surf[21]+alpha_surf[18])+0.22360679774997902*(alpha_surf[20]+alpha_surf[17])+0.29999999999999993*alpha_surf[19]+0.45*(alpha_surf[15]+alpha_surf[10])-0.33541019662496846*(alpha_surf[14]+alpha_surf[13]+alpha_surf[9]+alpha_surf[7])+0.22360679774997896*(alpha_surf[12]+alpha_surf[11])+0.33541019662496846*(alpha_surf[8]+alpha_surf[6]+alpha_surf[4]+alpha_surf[3])-0.25*(alpha_surf[5]+alpha_surf[2])+0.25*(alpha_surf[1]+alpha_surf[0]) > 0.) 
    sgn_alpha_surf[26] = 1.0; 
  else  
    sgn_alpha_surf[26] = -1.0; 
  
  if (sgn_alpha_surf[26] == sgn_alpha_surf[25]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2*(alpha_surf[35]+alpha_surf[30])+0.19999999999999954*(alpha_surf[34]+alpha_surf[33])-0.2999999999999998*alpha_surf[32]-0.29999999999999993*(alpha_surf[31]+alpha_surf[19])-0.3*(alpha_surf[29]+alpha_surf[28]+alpha_surf[26]+alpha_surf[25])+0.22360679774997896*(alpha_surf[27]+alpha_surf[24]+alpha_surf[12]+alpha_surf[11])+0.45*(alpha_surf[23]+alpha_surf[16]+alpha_surf[15]+alpha_surf[10])-0.29999999999999977*alpha_surf[22]+0.22360679774997902*(alpha_surf[21]+alpha_surf[20]+alpha_surf[18]+alpha_surf[17])-0.33541019662496846*(alpha_surf[14]+alpha_surf[13]+alpha_surf[9]+alpha_surf[8]+alpha_surf[7]+alpha_surf[6]+alpha_surf[4]+alpha_surf[3])+0.25*(alpha_surf[5]+alpha_surf[2]+alpha_surf[1]+alpha_surf[0]) > 0.) 
    sgn_alpha_surf[27] = 1.0; 
  else  
    sgn_alpha_surf[27] = -1.0; 
  
  if (sgn_alpha_surf[27] == sgn_alpha_surf[26]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.25*(alpha_surf[35]+alpha_surf[30]))-0.24999999999999994*(alpha_surf[34]+alpha_surf[33])+0.375*(alpha_surf[32]+alpha_surf[29]+alpha_surf[28]+alpha_surf[22])-0.2795084971874737*(alpha_surf[27]+alpha_surf[12])+0.22360679774997896*(alpha_surf[24]+alpha_surf[11])-0.2795084971874738*(alpha_surf[21]+alpha_surf[20])+0.22360679774997902*(alpha_surf[18]+alpha_surf[17])-0.33541019662496846*(alpha_surf[13]+alpha_surf[7]+alpha_surf[6]+alpha_surf[3])+0.25*(alpha_surf[5]+alpha_surf[2]+alpha_surf[1]+alpha_surf[0]) > 0.) 
    sgn_alpha_surf[28] = 1.0; 
  else  
    sgn_alpha_surf[28] = -1.0; 
  
  if (sgn_alpha_surf[28] == sgn_alpha_surf[27]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2*(alpha_surf[35]+alpha_surf[30])+0.19999999999999954*(alpha_surf[34]+alpha_surf[33])-0.2999999999999998*alpha_surf[32]+0.29999999999999993*(alpha_surf[31]+alpha_surf[19])-0.3*(alpha_surf[29]+alpha_surf[28])+0.22360679774997896*(alpha_surf[27]+alpha_surf[24]+alpha_surf[12]+alpha_surf[11])+0.3*(alpha_surf[26]+alpha_surf[25])-0.45*(alpha_surf[23]+alpha_surf[16]+alpha_surf[15]+alpha_surf[10])-0.29999999999999977*alpha_surf[22]+0.22360679774997902*(alpha_surf[21]+alpha_surf[20]+alpha_surf[18]+alpha_surf[17])+0.33541019662496846*(alpha_surf[14]+alpha_surf[9]+alpha_surf[8]+alpha_surf[4])-0.33541019662496846*(alpha_surf[13]+alpha_surf[7]+alpha_surf[6]+alpha_surf[3])+0.25*(alpha_surf[5]+alpha_surf[2]+alpha_surf[1]+alpha_surf[0]) > 0.) 
    sgn_alpha_surf[29] = 1.0; 
  else  
    sgn_alpha_surf[29] = -1.0; 
  
  if (sgn_alpha_surf[29] == sgn_alpha_surf[28]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.25*(alpha_surf[35]+alpha_surf[34]+alpha_surf[33]+alpha_surf[30]))+0.375*(alpha_surf[31]+alpha_surf[26]+alpha_surf[25]+alpha_surf[19])+0.22360679774997896*(alpha_surf[27]+alpha_surf[12])-0.2795084971874737*(alpha_surf[24]+alpha_surf[11])+0.22360679774997902*(alpha_surf[21]+alpha_surf[20])-0.2795084971874738*(alpha_surf[18]+alpha_surf[17])-0.33541019662496846*(alpha_surf[14]+alpha_surf[9]+alpha_surf[8]+alpha_surf[4])+0.25*(alpha_surf[5]+alpha_surf[2]+alpha_surf[1]+alpha_surf[0]) > 0.) 
    sgn_alpha_surf[30] = 1.0; 
  else  
    sgn_alpha_surf[30] = -1.0; 
  
  if (sgn_alpha_surf[30] == sgn_alpha_surf[29]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3125*(alpha_surf[35]+alpha_surf[34]+alpha_surf[33]+alpha_surf[30])-0.2795084971874737*(alpha_surf[27]+alpha_surf[24]+alpha_surf[12]+alpha_surf[11])-0.2795084971874738*(alpha_surf[21]+alpha_surf[20]+alpha_surf[18]+alpha_surf[17])+0.25*(alpha_surf[5]+alpha_surf[2]+alpha_surf[1]+alpha_surf[0]) > 0.) 
    sgn_alpha_surf[31] = 1.0; 
  else  
    sgn_alpha_surf[31] = -1.0; 
  
  if (sgn_alpha_surf[31] == sgn_alpha_surf[30]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.25*(alpha_surf[35]+alpha_surf[34]+alpha_surf[33]+alpha_surf[30]))-0.375*(alpha_surf[31]+alpha_surf[26]+alpha_surf[25]+alpha_surf[19])+0.22360679774997896*(alpha_surf[27]+alpha_surf[12])-0.2795084971874737*(alpha_surf[24]+alpha_surf[11])+0.22360679774997902*(alpha_surf[21]+alpha_surf[20])-0.2795084971874738*(alpha_surf[18]+alpha_surf[17])+0.33541019662496846*(alpha_surf[14]+alpha_surf[9]+alpha_surf[8]+alpha_surf[4])+0.25*(alpha_surf[5]+alpha_surf[2]+alpha_surf[1]+alpha_surf[0]) > 0.) 
    sgn_alpha_surf[32] = 1.0; 
  else  
    sgn_alpha_surf[32] = -1.0; 
  
  if (sgn_alpha_surf[32] == sgn_alpha_surf[31]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2*(alpha_surf[35]+alpha_surf[30])+0.19999999999999954*(alpha_surf[34]+alpha_surf[33])+0.2999999999999998*alpha_surf[32]-0.29999999999999993*(alpha_surf[31]+alpha_surf[19])+0.3*(alpha_surf[29]+alpha_surf[28])+0.22360679774997896*(alpha_surf[27]+alpha_surf[24]+alpha_surf[12]+alpha_surf[11])-0.3*(alpha_surf[26]+alpha_surf[25])-0.45*(alpha_surf[23]+alpha_surf[16]+alpha_surf[15]+alpha_surf[10])+0.29999999999999977*alpha_surf[22]+0.22360679774997902*(alpha_surf[21]+alpha_surf[20]+alpha_surf[18]+alpha_surf[17])-0.33541019662496846*(alpha_surf[14]+alpha_surf[9]+alpha_surf[8]+alpha_surf[4])+0.33541019662496846*(alpha_surf[13]+alpha_surf[7]+alpha_surf[6]+alpha_surf[3])+0.25*(alpha_surf[5]+alpha_surf[2]+alpha_surf[1]+alpha_surf[0]) > 0.) 
    sgn_alpha_surf[33] = 1.0; 
  else  
    sgn_alpha_surf[33] = -1.0; 
  
  if (sgn_alpha_surf[33] == sgn_alpha_surf[32]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.25*(alpha_surf[35]+alpha_surf[30]))-0.24999999999999994*(alpha_surf[34]+alpha_surf[33])-0.375*(alpha_surf[32]+alpha_surf[29]+alpha_surf[28]+alpha_surf[22])-0.2795084971874737*(alpha_surf[27]+alpha_surf[12])+0.22360679774997896*(alpha_surf[24]+alpha_surf[11])-0.2795084971874738*(alpha_surf[21]+alpha_surf[20])+0.22360679774997902*(alpha_surf[18]+alpha_surf[17])+0.33541019662496846*(alpha_surf[13]+alpha_surf[7]+alpha_surf[6]+alpha_surf[3])+0.25*(alpha_surf[5]+alpha_surf[2]+alpha_surf[1]+alpha_surf[0]) > 0.) 
    sgn_alpha_surf[34] = 1.0; 
  else  
    sgn_alpha_surf[34] = -1.0; 
  
  if (sgn_alpha_surf[34] == sgn_alpha_surf[33]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2*(alpha_surf[35]+alpha_surf[30])+0.19999999999999954*(alpha_surf[34]+alpha_surf[33])+0.2999999999999998*alpha_surf[32]+0.29999999999999993*(alpha_surf[31]+alpha_surf[19])+0.3*(alpha_surf[29]+alpha_surf[28]+alpha_surf[26]+alpha_surf[25])+0.22360679774997896*(alpha_surf[27]+alpha_surf[24]+alpha_surf[12]+alpha_surf[11])+0.45*(alpha_surf[23]+alpha_surf[16]+alpha_surf[15]+alpha_surf[10])+0.29999999999999977*alpha_surf[22]+0.22360679774997902*(alpha_surf[21]+alpha_surf[20]+alpha_surf[18]+alpha_surf[17])+0.33541019662496846*(alpha_surf[14]+alpha_surf[13]+alpha_surf[9]+alpha_surf[8]+alpha_surf[7]+alpha_surf[6]+alpha_surf[4]+alpha_surf[3])+0.25*(alpha_surf[5]+alpha_surf[2]+alpha_surf[1]+alpha_surf[0]) > 0.) 
    sgn_alpha_surf[35] = 1.0; 
  else  
    sgn_alpha_surf[35] = -1.0; 
  
  if (sgn_alpha_surf[35] == sgn_alpha_surf[34]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  *const_sgn_drag_coeff_surf = const_sgn_alpha_surf; 
} 
