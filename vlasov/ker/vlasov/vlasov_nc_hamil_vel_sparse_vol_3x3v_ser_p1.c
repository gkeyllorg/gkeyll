#include <gkyl_nc_hamil_gen_vol_comps_kernels.h> 
#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH GKYL_NOINLINE static void dot_pb_3x3v_ser_p1(const double *list_1, const double *list_2, double* GKYL_RESTRICT out) 
{ 
  for (int i = 0; i<64; ++i) out[0] += list_1[i]*list_2[i]; 
} 
  
GKYL_CU_DH GKYL_NOINLINE static void zero_pb_list(double *pb_list) 
{ 
  for (int i = 0; i<64; ++i) pb_list[i] = 0.0; 
} 
  
GKYL_CU_DH void vlasov_nc_hamil_vel_sparse_vol_3x3v_ser_p1(const double *w, const double *dxv, const double *vmap,
  const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f,
  const double *f_no_J, double* GKYL_RESTRICT out) 
{ 
  double alphai[64] = {0.0}; 
  double bracket_w_Hi[64] = {0.0}; 
  for (int dir = 0; dir<6; ++dir) { 
    // Compute the per-direction charateristic 
    zero_pb_list(alphai); 
    vlasov_nc_hamil_vel_sparse_vol_alpha_3x3v_ser_p1(w, dxv, dir, vmap, jacob_pos, jacob_vel, poisson_tensor_conf, hamil, alphai); 
    
    // Compute the modal basis backet coeffcients, accumulate the outputs 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp0_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[0]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp1_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[1]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp2_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[2]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp3_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[3]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp4_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[4]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp5_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[5]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp6_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[6]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp7_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[7]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp8_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[8]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp9_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[9]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp10_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[10]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp11_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[11]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp12_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[12]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp13_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[13]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp14_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[14]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp15_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[15]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp16_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[16]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp17_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[17]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp18_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[18]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp19_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[19]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp20_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[20]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp21_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[21]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp22_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[22]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp23_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[23]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp24_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[24]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp25_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[25]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp26_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[26]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp27_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[27]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp28_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[28]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp29_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[29]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp30_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[30]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp31_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[31]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp32_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[32]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp33_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[33]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp34_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[34]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp35_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[35]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp36_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[36]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp37_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[37]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp38_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[38]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp39_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[39]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp40_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[40]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp41_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[41]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp42_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[42]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp43_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[43]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp44_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[44]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp45_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[45]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp46_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[46]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp47_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[47]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp48_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[48]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp49_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[49]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp50_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[50]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp51_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[51]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp52_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[52]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp53_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[53]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp54_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[54]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp55_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[55]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp56_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[56]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp57_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[57]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp58_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[58]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp59_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[59]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp60_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[60]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp61_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[61]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp62_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[62]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp63_3x3v_ser_p1(w, dxv, dir, jacob_pos, jacob_vel, alphai, bracket_w_Hi); 
    dot_pb_3x3v_ser_p1(f, bracket_w_Hi, &out[63]); 
    
  } 
} 
