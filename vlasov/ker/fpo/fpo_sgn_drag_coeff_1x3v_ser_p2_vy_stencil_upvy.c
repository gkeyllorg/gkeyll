#include <gkyl_fpo_vlasov_kernels.h> 

GKYL_CU_DH void fpo_sgn_drag_coeff_1x3v_vy_ser_p2_upvy(const double *drag_coeff_surf, double *sgn_drag_coeff_surf, int *const_sgn_drag_coeff_surf) {
  // drag_coeff_surf: Surface projection of drag coefficient at lower boundary.
  // sgn_drag_coeff_surf: sign(drag_coeff_surf) evaluated at quadrature points along lower surface.
  // const_sgn_drag_coeff_surf: 1 if sign(drag_coeff_surf) is constant along lower boundary, 0 otherwise. 

  const double *alpha_surf = &drag_coeff_surf[20]; 
  double *sgn_alpha_surf = &sgn_drag_coeff_surf[27]; 

  int const_sgn_alpha_surf = 1;  
  
  if (0.5692099788303082*(alpha_surf[19]+alpha_surf[18]+alpha_surf[17])-0.42426406871192807*(alpha_surf[16]+alpha_surf[15]+alpha_surf[12])-0.42426406871192845*(alpha_surf[14]+alpha_surf[13]+alpha_surf[11])-0.853814968245462*alpha_surf[10]+0.3162277660168379*(alpha_surf[9]+alpha_surf[8]+alpha_surf[7])+0.6363961030678926*(alpha_surf[6]+alpha_surf[5]+alpha_surf[4])-0.4743416490252568*(alpha_surf[3]+alpha_surf[2]+alpha_surf[1])+0.3535533905932737*alpha_surf[0] > 0.) 
    sgn_alpha_surf[0] = 1.0; 
  else  
    sgn_alpha_surf[0] = -1.0; 
  
  if (-(0.711512473537885*alpha_surf[19])+0.5303300858899104*(alpha_surf[16]+alpha_surf[15])-0.42426406871192807*alpha_surf[12]-0.42426406871192845*alpha_surf[11]-0.3952847075210473*alpha_surf[9]+0.3162277660168379*(alpha_surf[8]+alpha_surf[7])+0.6363961030678926*alpha_surf[4]-0.4743416490252568*(alpha_surf[2]+alpha_surf[1])+0.3535533905932737*alpha_surf[0] > 0.) 
    sgn_alpha_surf[1] = 1.0; 
  else  
    sgn_alpha_surf[1] = -1.0; 
  
  if (sgn_alpha_surf[1] == sgn_alpha_surf[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5692099788303082*alpha_surf[19]-0.5692099788303082*(alpha_surf[18]+alpha_surf[17])-0.42426406871192807*(alpha_surf[16]+alpha_surf[15]+alpha_surf[12])+0.42426406871192845*(alpha_surf[14]+alpha_surf[13])-0.42426406871192845*alpha_surf[11]+0.853814968245462*alpha_surf[10]+0.3162277660168379*(alpha_surf[9]+alpha_surf[8]+alpha_surf[7])-0.6363961030678926*(alpha_surf[6]+alpha_surf[5])+0.6363961030678926*alpha_surf[4]+0.4743416490252568*alpha_surf[3]-0.4743416490252568*(alpha_surf[2]+alpha_surf[1])+0.3535533905932737*alpha_surf[0] > 0.) 
    sgn_alpha_surf[2] = 1.0; 
  else  
    sgn_alpha_surf[2] = -1.0; 
  
  if (sgn_alpha_surf[2] == sgn_alpha_surf[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.711512473537885*alpha_surf[18])-0.42426406871192807*alpha_surf[15]+0.5303300858899104*(alpha_surf[14]+alpha_surf[12])-0.42426406871192845*alpha_surf[13]+0.3162277660168379*(alpha_surf[9]+alpha_surf[7])-0.3952847075210473*alpha_surf[8]+0.6363961030678926*alpha_surf[5]-0.4743416490252568*(alpha_surf[3]+alpha_surf[1])+0.3535533905932737*alpha_surf[0] > 0.) 
    sgn_alpha_surf[3] = 1.0; 
  else  
    sgn_alpha_surf[3] = -1.0; 
  
  if (sgn_alpha_surf[3] == sgn_alpha_surf[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5303300858899104*(alpha_surf[15]+alpha_surf[12])-0.3952847075210473*(alpha_surf[9]+alpha_surf[8])+0.3162277660168379*alpha_surf[7]-0.4743416490252568*alpha_surf[1]+0.3535533905932737*alpha_surf[0] > 0.) 
    sgn_alpha_surf[4] = 1.0; 
  else  
    sgn_alpha_surf[4] = -1.0; 
  
  if (sgn_alpha_surf[4] == sgn_alpha_surf[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.711512473537885*alpha_surf[18]-0.42426406871192807*alpha_surf[15]-0.5303300858899104*alpha_surf[14]+0.42426406871192845*alpha_surf[13]+0.5303300858899104*alpha_surf[12]+0.3162277660168379*(alpha_surf[9]+alpha_surf[7])-0.3952847075210473*alpha_surf[8]-0.6363961030678926*alpha_surf[5]+0.4743416490252568*alpha_surf[3]-0.4743416490252568*alpha_surf[1]+0.3535533905932737*alpha_surf[0] > 0.) 
    sgn_alpha_surf[5] = 1.0; 
  else  
    sgn_alpha_surf[5] = -1.0; 
  
  if (sgn_alpha_surf[5] == sgn_alpha_surf[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.5692099788303082*(alpha_surf[19]+alpha_surf[17]))+0.5692099788303082*alpha_surf[18]+0.42426406871192807*alpha_surf[16]-0.42426406871192807*(alpha_surf[15]+alpha_surf[12])-0.42426406871192845*(alpha_surf[14]+alpha_surf[13])+0.42426406871192845*alpha_surf[11]+0.853814968245462*alpha_surf[10]+0.3162277660168379*(alpha_surf[9]+alpha_surf[8]+alpha_surf[7])-0.6363961030678926*(alpha_surf[6]+alpha_surf[4])+0.6363961030678926*alpha_surf[5]-0.4743416490252568*(alpha_surf[3]+alpha_surf[1])+0.4743416490252568*alpha_surf[2]+0.3535533905932737*alpha_surf[0] > 0.) 
    sgn_alpha_surf[6] = 1.0; 
  else  
    sgn_alpha_surf[6] = -1.0; 
  
  if (sgn_alpha_surf[6] == sgn_alpha_surf[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.711512473537885*alpha_surf[19]-0.5303300858899104*alpha_surf[16]+0.5303300858899104*alpha_surf[15]-0.42426406871192807*alpha_surf[12]+0.42426406871192845*alpha_surf[11]-0.3952847075210473*alpha_surf[9]+0.3162277660168379*(alpha_surf[8]+alpha_surf[7])-0.6363961030678926*alpha_surf[4]+0.4743416490252568*alpha_surf[2]-0.4743416490252568*alpha_surf[1]+0.3535533905932737*alpha_surf[0] > 0.) 
    sgn_alpha_surf[7] = 1.0; 
  else  
    sgn_alpha_surf[7] = -1.0; 
  
  if (sgn_alpha_surf[7] == sgn_alpha_surf[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.5692099788303082*(alpha_surf[19]+alpha_surf[18]))+0.5692099788303082*alpha_surf[17]+0.42426406871192807*alpha_surf[16]-0.42426406871192807*(alpha_surf[15]+alpha_surf[12])+0.42426406871192845*(alpha_surf[14]+alpha_surf[13]+alpha_surf[11])-0.853814968245462*alpha_surf[10]+0.3162277660168379*(alpha_surf[9]+alpha_surf[8]+alpha_surf[7])+0.6363961030678926*alpha_surf[6]-0.6363961030678926*(alpha_surf[5]+alpha_surf[4])+0.4743416490252568*(alpha_surf[3]+alpha_surf[2])-0.4743416490252568*alpha_surf[1]+0.3535533905932737*alpha_surf[0] > 0.) 
    sgn_alpha_surf[8] = 1.0; 
  else  
    sgn_alpha_surf[8] = -1.0; 
  
  if (sgn_alpha_surf[8] == sgn_alpha_surf[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.711512473537885*alpha_surf[17])-0.42426406871192807*alpha_surf[16]-0.42426406871192845*alpha_surf[14]+0.5303300858899104*(alpha_surf[13]+alpha_surf[11])+0.3162277660168379*(alpha_surf[9]+alpha_surf[8])-0.3952847075210473*alpha_surf[7]+0.6363961030678926*alpha_surf[6]-0.4743416490252568*(alpha_surf[3]+alpha_surf[2])+0.3535533905932737*alpha_surf[0] > 0.) 
    sgn_alpha_surf[9] = 1.0; 
  else  
    sgn_alpha_surf[9] = -1.0; 
  
  if (sgn_alpha_surf[9] == sgn_alpha_surf[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5303300858899104*(alpha_surf[16]+alpha_surf[11])-0.3952847075210473*(alpha_surf[9]+alpha_surf[7])+0.3162277660168379*alpha_surf[8]-0.4743416490252568*alpha_surf[2]+0.3535533905932737*alpha_surf[0] > 0.) 
    sgn_alpha_surf[10] = 1.0; 
  else  
    sgn_alpha_surf[10] = -1.0; 
  
  if (sgn_alpha_surf[10] == sgn_alpha_surf[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.711512473537885*alpha_surf[17]-0.42426406871192807*alpha_surf[16]+0.42426406871192845*alpha_surf[14]-0.5303300858899104*alpha_surf[13]+0.5303300858899104*alpha_surf[11]+0.3162277660168379*(alpha_surf[9]+alpha_surf[8])-0.3952847075210473*alpha_surf[7]-0.6363961030678926*alpha_surf[6]+0.4743416490252568*alpha_surf[3]-0.4743416490252568*alpha_surf[2]+0.3535533905932737*alpha_surf[0] > 0.) 
    sgn_alpha_surf[11] = 1.0; 
  else  
    sgn_alpha_surf[11] = -1.0; 
  
  if (sgn_alpha_surf[11] == sgn_alpha_surf[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5303300858899104*(alpha_surf[14]+alpha_surf[13])+0.3162277660168379*alpha_surf[9]-0.3952847075210473*(alpha_surf[8]+alpha_surf[7])-0.4743416490252568*alpha_surf[3]+0.3535533905932737*alpha_surf[0] > 0.) 
    sgn_alpha_surf[12] = 1.0; 
  else  
    sgn_alpha_surf[12] = -1.0; 
  
  if (sgn_alpha_surf[12] == sgn_alpha_surf[11]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932737*alpha_surf[0]-0.3952847075210473*(alpha_surf[9]+alpha_surf[8]+alpha_surf[7]) > 0.) 
    sgn_alpha_surf[13] = 1.0; 
  else  
    sgn_alpha_surf[13] = -1.0; 
  
  if (sgn_alpha_surf[13] == sgn_alpha_surf[12]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.5303300858899104*(alpha_surf[14]+alpha_surf[13]))+0.3162277660168379*alpha_surf[9]-0.3952847075210473*(alpha_surf[8]+alpha_surf[7])+0.4743416490252568*alpha_surf[3]+0.3535533905932737*alpha_surf[0] > 0.) 
    sgn_alpha_surf[14] = 1.0; 
  else  
    sgn_alpha_surf[14] = -1.0; 
  
  if (sgn_alpha_surf[14] == sgn_alpha_surf[13]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.711512473537885*alpha_surf[17]+0.42426406871192807*alpha_surf[16]-0.42426406871192845*alpha_surf[14]+0.5303300858899104*alpha_surf[13]-0.5303300858899104*alpha_surf[11]+0.3162277660168379*(alpha_surf[9]+alpha_surf[8])-0.3952847075210473*alpha_surf[7]-0.6363961030678926*alpha_surf[6]-0.4743416490252568*alpha_surf[3]+0.4743416490252568*alpha_surf[2]+0.3535533905932737*alpha_surf[0] > 0.) 
    sgn_alpha_surf[15] = 1.0; 
  else  
    sgn_alpha_surf[15] = -1.0; 
  
  if (sgn_alpha_surf[15] == sgn_alpha_surf[14]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.5303300858899104*(alpha_surf[16]+alpha_surf[11]))-0.3952847075210473*(alpha_surf[9]+alpha_surf[7])+0.3162277660168379*alpha_surf[8]+0.4743416490252568*alpha_surf[2]+0.3535533905932737*alpha_surf[0] > 0.) 
    sgn_alpha_surf[16] = 1.0; 
  else  
    sgn_alpha_surf[16] = -1.0; 
  
  if (sgn_alpha_surf[16] == sgn_alpha_surf[15]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.711512473537885*alpha_surf[17])+0.42426406871192807*alpha_surf[16]+0.42426406871192845*alpha_surf[14]-0.5303300858899104*(alpha_surf[13]+alpha_surf[11])+0.3162277660168379*(alpha_surf[9]+alpha_surf[8])-0.3952847075210473*alpha_surf[7]+0.6363961030678926*alpha_surf[6]+0.4743416490252568*(alpha_surf[3]+alpha_surf[2])+0.3535533905932737*alpha_surf[0] > 0.) 
    sgn_alpha_surf[17] = 1.0; 
  else  
    sgn_alpha_surf[17] = -1.0; 
  
  if (sgn_alpha_surf[17] == sgn_alpha_surf[16]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.5692099788303082*(alpha_surf[19]+alpha_surf[18]))+0.5692099788303082*alpha_surf[17]-0.42426406871192807*alpha_surf[16]+0.42426406871192807*(alpha_surf[15]+alpha_surf[12])-0.42426406871192845*(alpha_surf[14]+alpha_surf[13]+alpha_surf[11])+0.853814968245462*alpha_surf[10]+0.3162277660168379*(alpha_surf[9]+alpha_surf[8]+alpha_surf[7])+0.6363961030678926*alpha_surf[6]-0.6363961030678926*(alpha_surf[5]+alpha_surf[4])-0.4743416490252568*(alpha_surf[3]+alpha_surf[2])+0.4743416490252568*alpha_surf[1]+0.3535533905932737*alpha_surf[0] > 0.) 
    sgn_alpha_surf[18] = 1.0; 
  else  
    sgn_alpha_surf[18] = -1.0; 
  
  if (sgn_alpha_surf[18] == sgn_alpha_surf[17]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.711512473537885*alpha_surf[19]+0.5303300858899104*alpha_surf[16]-0.5303300858899104*alpha_surf[15]+0.42426406871192807*alpha_surf[12]-0.42426406871192845*alpha_surf[11]-0.3952847075210473*alpha_surf[9]+0.3162277660168379*(alpha_surf[8]+alpha_surf[7])-0.6363961030678926*alpha_surf[4]-0.4743416490252568*alpha_surf[2]+0.4743416490252568*alpha_surf[1]+0.3535533905932737*alpha_surf[0] > 0.) 
    sgn_alpha_surf[19] = 1.0; 
  else  
    sgn_alpha_surf[19] = -1.0; 
  
  if (sgn_alpha_surf[19] == sgn_alpha_surf[18]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.5692099788303082*(alpha_surf[19]+alpha_surf[17]))+0.5692099788303082*alpha_surf[18]-0.42426406871192807*alpha_surf[16]+0.42426406871192807*(alpha_surf[15]+alpha_surf[12])+0.42426406871192845*(alpha_surf[14]+alpha_surf[13])-0.42426406871192845*alpha_surf[11]-0.853814968245462*alpha_surf[10]+0.3162277660168379*(alpha_surf[9]+alpha_surf[8]+alpha_surf[7])-0.6363961030678926*(alpha_surf[6]+alpha_surf[4])+0.6363961030678926*alpha_surf[5]+0.4743416490252568*(alpha_surf[3]+alpha_surf[1])-0.4743416490252568*alpha_surf[2]+0.3535533905932737*alpha_surf[0] > 0.) 
    sgn_alpha_surf[20] = 1.0; 
  else  
    sgn_alpha_surf[20] = -1.0; 
  
  if (sgn_alpha_surf[20] == sgn_alpha_surf[19]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.711512473537885*alpha_surf[18]+0.42426406871192807*alpha_surf[15]+0.5303300858899104*alpha_surf[14]-0.42426406871192845*alpha_surf[13]-0.5303300858899104*alpha_surf[12]+0.3162277660168379*(alpha_surf[9]+alpha_surf[7])-0.3952847075210473*alpha_surf[8]-0.6363961030678926*alpha_surf[5]-0.4743416490252568*alpha_surf[3]+0.4743416490252568*alpha_surf[1]+0.3535533905932737*alpha_surf[0] > 0.) 
    sgn_alpha_surf[21] = 1.0; 
  else  
    sgn_alpha_surf[21] = -1.0; 
  
  if (sgn_alpha_surf[21] == sgn_alpha_surf[20]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.5303300858899104*(alpha_surf[15]+alpha_surf[12]))-0.3952847075210473*(alpha_surf[9]+alpha_surf[8])+0.3162277660168379*alpha_surf[7]+0.4743416490252568*alpha_surf[1]+0.3535533905932737*alpha_surf[0] > 0.) 
    sgn_alpha_surf[22] = 1.0; 
  else  
    sgn_alpha_surf[22] = -1.0; 
  
  if (sgn_alpha_surf[22] == sgn_alpha_surf[21]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.711512473537885*alpha_surf[18])+0.42426406871192807*alpha_surf[15]-0.5303300858899104*(alpha_surf[14]+alpha_surf[12])+0.42426406871192845*alpha_surf[13]+0.3162277660168379*(alpha_surf[9]+alpha_surf[7])-0.3952847075210473*alpha_surf[8]+0.6363961030678926*alpha_surf[5]+0.4743416490252568*(alpha_surf[3]+alpha_surf[1])+0.3535533905932737*alpha_surf[0] > 0.) 
    sgn_alpha_surf[23] = 1.0; 
  else  
    sgn_alpha_surf[23] = -1.0; 
  
  if (sgn_alpha_surf[23] == sgn_alpha_surf[22]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5692099788303082*alpha_surf[19]-0.5692099788303082*(alpha_surf[18]+alpha_surf[17])+0.42426406871192807*(alpha_surf[16]+alpha_surf[15]+alpha_surf[12])-0.42426406871192845*(alpha_surf[14]+alpha_surf[13])+0.42426406871192845*alpha_surf[11]-0.853814968245462*alpha_surf[10]+0.3162277660168379*(alpha_surf[9]+alpha_surf[8]+alpha_surf[7])-0.6363961030678926*(alpha_surf[6]+alpha_surf[5])+0.6363961030678926*alpha_surf[4]-0.4743416490252568*alpha_surf[3]+0.4743416490252568*(alpha_surf[2]+alpha_surf[1])+0.3535533905932737*alpha_surf[0] > 0.) 
    sgn_alpha_surf[24] = 1.0; 
  else  
    sgn_alpha_surf[24] = -1.0; 
  
  if (sgn_alpha_surf[24] == sgn_alpha_surf[23]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.711512473537885*alpha_surf[19])-0.5303300858899104*(alpha_surf[16]+alpha_surf[15])+0.42426406871192807*alpha_surf[12]+0.42426406871192845*alpha_surf[11]-0.3952847075210473*alpha_surf[9]+0.3162277660168379*(alpha_surf[8]+alpha_surf[7])+0.6363961030678926*alpha_surf[4]+0.4743416490252568*(alpha_surf[2]+alpha_surf[1])+0.3535533905932737*alpha_surf[0] > 0.) 
    sgn_alpha_surf[25] = 1.0; 
  else  
    sgn_alpha_surf[25] = -1.0; 
  
  if (sgn_alpha_surf[25] == sgn_alpha_surf[24]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5692099788303082*(alpha_surf[19]+alpha_surf[18]+alpha_surf[17])+0.42426406871192807*(alpha_surf[16]+alpha_surf[15]+alpha_surf[12])+0.42426406871192845*(alpha_surf[14]+alpha_surf[13]+alpha_surf[11])+0.853814968245462*alpha_surf[10]+0.3162277660168379*(alpha_surf[9]+alpha_surf[8]+alpha_surf[7])+0.6363961030678926*(alpha_surf[6]+alpha_surf[5]+alpha_surf[4])+0.4743416490252568*(alpha_surf[3]+alpha_surf[2]+alpha_surf[1])+0.3535533905932737*alpha_surf[0] > 0.) 
    sgn_alpha_surf[26] = 1.0; 
  else  
    sgn_alpha_surf[26] = -1.0; 
  
  if (sgn_alpha_surf[26] == sgn_alpha_surf[25]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  *const_sgn_drag_coeff_surf = const_sgn_alpha_surf; 
} 
