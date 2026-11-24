#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfy_2x3v_ser_p1(const double *w, const double *dxv,
  const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel,
  const double *poisson_tensor_conf, const double *hamil,
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  double dx11 = 2.0/dxv[1]; 
  double dv11 = 2.0/dxv[3]; 
  double wv = w[3]; 

  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  const double *jacob_cy_l = &jacob_pos_l[2]; 
  const double *jacob_cy_c = &jacob_pos_c[2]; 
  const double *jacob_cy_r = &jacob_pos_r[2]; 
  const double jacob_cy_l_inv = 1.0/jacob_cy_l[0]; 
  const double jacob_cy_c_inv = 1.0/jacob_cy_c[0]; 
  const double jacob_cy_r_inv = 1.0/jacob_cy_r[0]; 
  double fUpwind_r[16] = {0.0}; 
  double fUpwind_l[16] = {0.0}; 
  double Ghat_r[16] = {0.0}; 
  double Ghat_l[16] = {0.0}; 
  if (wv>0) { 

  fUpwind_r[0] = (1.224744871391589*fc[2]+0.7071067811865475*fc[0])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[1] = (1.224744871391589*fc[6]+0.7071067811865475*fc[1])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[2] = (1.224744871391589*fc[8]+0.7071067811865475*fc[3])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[3] = (1.224744871391589*fc[10]+0.7071067811865475*fc[4])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[4] = (1.224744871391589*fc[13]+0.7071067811865475*fc[5])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[5] = (1.224744871391589*fc[16]+0.7071067811865475*fc[7])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[6] = (1.224744871391589*fc[17]+0.7071067811865475*fc[9])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[7] = (1.224744871391589*fc[19]+0.7071067811865475*fc[11])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[8] = (1.224744871391589*fc[20]+0.7071067811865475*fc[12])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[9] = (1.224744871391589*fc[22]+0.7071067811865475*fc[14])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[10] = (1.224744871391589*fc[24]+0.7071067811865475*fc[15])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[11] = (1.224744871391589*fc[26]+0.7071067811865475*fc[18])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[12] = (1.224744871391589*fc[27]+0.7071067811865475*fc[21])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[13] = (1.224744871391589*fc[28]+0.7071067811865475*fc[23])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[14] = (1.224744871391589*fc[30]+0.7071067811865475*fc[25])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[15] = (1.224744871391589*fc[31]+0.7071067811865475*fc[29])*jacob_cy_c_inv*jacob_vy_inv; 

  fUpwind_l[0] = (1.224744871391589*fl[2]+0.7071067811865475*fl[0])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[1] = (1.224744871391589*fl[6]+0.7071067811865475*fl[1])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[2] = (1.224744871391589*fl[8]+0.7071067811865475*fl[3])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[3] = (1.224744871391589*fl[10]+0.7071067811865475*fl[4])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[4] = (1.224744871391589*fl[13]+0.7071067811865475*fl[5])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[5] = (1.224744871391589*fl[16]+0.7071067811865475*fl[7])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[6] = (1.224744871391589*fl[17]+0.7071067811865475*fl[9])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[7] = (1.224744871391589*fl[19]+0.7071067811865475*fl[11])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[8] = (1.224744871391589*fl[20]+0.7071067811865475*fl[12])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[9] = (1.224744871391589*fl[22]+0.7071067811865475*fl[14])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[10] = (1.224744871391589*fl[24]+0.7071067811865475*fl[15])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[11] = (1.224744871391589*fl[26]+0.7071067811865475*fl[18])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[12] = (1.224744871391589*fl[27]+0.7071067811865475*fl[21])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[13] = (1.224744871391589*fl[28]+0.7071067811865475*fl[23])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[14] = (1.224744871391589*fl[30]+0.7071067811865475*fl[25])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[15] = (1.224744871391589*fl[31]+0.7071067811865475*fl[29])*jacob_cy_l_inv*jacob_vy_inv; 

  } else { 

  fUpwind_r[0] = (0.7071067811865475*fr[0]-1.224744871391589*fr[2])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[1] = (0.7071067811865475*fr[1]-1.224744871391589*fr[6])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[2] = (0.7071067811865475*fr[3]-1.224744871391589*fr[8])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[3] = (0.7071067811865475*fr[4]-1.224744871391589*fr[10])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[4] = (0.7071067811865475*fr[5]-1.224744871391589*fr[13])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[5] = (0.7071067811865475*fr[7]-1.224744871391589*fr[16])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[6] = (0.7071067811865475*fr[9]-1.224744871391589*fr[17])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[7] = (0.7071067811865475*fr[11]-1.224744871391589*fr[19])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[8] = (0.7071067811865475*fr[12]-1.224744871391589*fr[20])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[9] = (0.7071067811865475*fr[14]-1.224744871391589*fr[22])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[10] = (0.7071067811865475*fr[15]-1.224744871391589*fr[24])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[11] = (0.7071067811865475*fr[18]-1.224744871391589*fr[26])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[12] = (0.7071067811865475*fr[21]-1.224744871391589*fr[27])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[13] = (0.7071067811865475*fr[23]-1.224744871391589*fr[28])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[14] = (0.7071067811865475*fr[25]-1.224744871391589*fr[30])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[15] = (0.7071067811865475*fr[29]-1.224744871391589*fr[31])*jacob_cy_r_inv*jacob_vy_inv; 

  fUpwind_l[0] = (0.7071067811865475*fc[0]-1.224744871391589*fc[2])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[1] = (0.7071067811865475*fc[1]-1.224744871391589*fc[6])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[2] = (0.7071067811865475*fc[3]-1.224744871391589*fc[8])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[3] = (0.7071067811865475*fc[4]-1.224744871391589*fc[10])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[4] = (0.7071067811865475*fc[5]-1.224744871391589*fc[13])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[5] = (0.7071067811865475*fc[7]-1.224744871391589*fc[16])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[6] = (0.7071067811865475*fc[9]-1.224744871391589*fc[17])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[7] = (0.7071067811865475*fc[11]-1.224744871391589*fc[19])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[8] = (0.7071067811865475*fc[12]-1.224744871391589*fc[20])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[9] = (0.7071067811865475*fc[14]-1.224744871391589*fc[22])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[10] = (0.7071067811865475*fc[15]-1.224744871391589*fc[24])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[11] = (0.7071067811865475*fc[18]-1.224744871391589*fc[26])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[12] = (0.7071067811865475*fc[21]-1.224744871391589*fc[27])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[13] = (0.7071067811865475*fc[23]-1.224744871391589*fc[28])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[14] = (0.7071067811865475*fc[25]-1.224744871391589*fc[30])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[15] = (0.7071067811865475*fc[29]-1.224744871391589*fc[31])*jacob_cy_c_inv*jacob_vy_inv; 

  } 
  Ghat_l[0] = 0.6123724356957944*(hamil[7]*fUpwind_l[9]+fUpwind_l[4]*hamil[6]+fUpwind_l[2]*hamil[4]+fUpwind_l[0]*hamil[2]); 
  Ghat_l[1] = 0.6123724356957944*(hamil[7]*fUpwind_l[12]+hamil[6]*fUpwind_l[8]+hamil[4]*fUpwind_l[5]+fUpwind_l[1]*hamil[2]); 
  Ghat_l[2] = 0.6123724356957944*(hamil[6]*fUpwind_l[9]+fUpwind_l[4]*hamil[7]+fUpwind_l[0]*hamil[4]+fUpwind_l[2]*hamil[2]); 
  Ghat_l[3] = 0.6123724356957944*(hamil[7]*fUpwind_l[14]+hamil[6]*fUpwind_l[10]+hamil[4]*fUpwind_l[7]+hamil[2]*fUpwind_l[3]); 
  Ghat_l[4] = 0.6123724356957944*(hamil[4]*fUpwind_l[9]+fUpwind_l[2]*hamil[7]+fUpwind_l[0]*hamil[6]+hamil[2]*fUpwind_l[4]); 
  Ghat_l[5] = 0.6123724356957944*(hamil[6]*fUpwind_l[12]+hamil[7]*fUpwind_l[8]+hamil[2]*fUpwind_l[5]+fUpwind_l[1]*hamil[4]); 
  Ghat_l[6] = 0.6123724356957944*(hamil[7]*fUpwind_l[15]+hamil[6]*fUpwind_l[13]+hamil[4]*fUpwind_l[11]+hamil[2]*fUpwind_l[6]); 
  Ghat_l[7] = 0.6123724356957944*(hamil[6]*fUpwind_l[14]+hamil[7]*fUpwind_l[10]+hamil[2]*fUpwind_l[7]+fUpwind_l[3]*hamil[4]); 
  Ghat_l[8] = 0.6123724356957944*(hamil[4]*fUpwind_l[12]+hamil[2]*fUpwind_l[8]+fUpwind_l[5]*hamil[7]+fUpwind_l[1]*hamil[6]); 
  Ghat_l[9] = 0.6123724356957944*(hamil[2]*fUpwind_l[9]+fUpwind_l[0]*hamil[7]+fUpwind_l[2]*hamil[6]+fUpwind_l[4]*hamil[4]); 
  Ghat_l[10] = 0.6123724356957944*(hamil[4]*fUpwind_l[14]+hamil[2]*fUpwind_l[10]+fUpwind_l[7]*hamil[7]+fUpwind_l[3]*hamil[6]); 
  Ghat_l[11] = 0.6123724356957944*(hamil[6]*fUpwind_l[15]+hamil[7]*fUpwind_l[13]+hamil[2]*fUpwind_l[11]+hamil[4]*fUpwind_l[6]); 
  Ghat_l[12] = 0.6123724356957944*(hamil[2]*fUpwind_l[12]+hamil[4]*fUpwind_l[8]+fUpwind_l[1]*hamil[7]+fUpwind_l[5]*hamil[6]); 
  Ghat_l[13] = 0.6123724356957944*(hamil[4]*fUpwind_l[15]+hamil[2]*fUpwind_l[13]+hamil[7]*fUpwind_l[11]+fUpwind_l[6]*hamil[6]); 
  Ghat_l[14] = 0.6123724356957944*(hamil[2]*fUpwind_l[14]+hamil[4]*fUpwind_l[10]+fUpwind_l[3]*hamil[7]+hamil[6]*fUpwind_l[7]); 
  Ghat_l[15] = 0.6123724356957944*(hamil[2]*fUpwind_l[15]+hamil[4]*fUpwind_l[13]+hamil[6]*fUpwind_l[11]+fUpwind_l[6]*hamil[7]); 

  Ghat_r[0] = 0.6123724356957944*(hamil[7]*fUpwind_r[9]+fUpwind_r[4]*hamil[6]+fUpwind_r[2]*hamil[4]+fUpwind_r[0]*hamil[2]); 
  Ghat_r[1] = 0.6123724356957944*(hamil[7]*fUpwind_r[12]+hamil[6]*fUpwind_r[8]+hamil[4]*fUpwind_r[5]+fUpwind_r[1]*hamil[2]); 
  Ghat_r[2] = 0.6123724356957944*(hamil[6]*fUpwind_r[9]+fUpwind_r[4]*hamil[7]+fUpwind_r[0]*hamil[4]+fUpwind_r[2]*hamil[2]); 
  Ghat_r[3] = 0.6123724356957944*(hamil[7]*fUpwind_r[14]+hamil[6]*fUpwind_r[10]+hamil[4]*fUpwind_r[7]+hamil[2]*fUpwind_r[3]); 
  Ghat_r[4] = 0.6123724356957944*(hamil[4]*fUpwind_r[9]+fUpwind_r[2]*hamil[7]+fUpwind_r[0]*hamil[6]+hamil[2]*fUpwind_r[4]); 
  Ghat_r[5] = 0.6123724356957944*(hamil[6]*fUpwind_r[12]+hamil[7]*fUpwind_r[8]+hamil[2]*fUpwind_r[5]+fUpwind_r[1]*hamil[4]); 
  Ghat_r[6] = 0.6123724356957944*(hamil[7]*fUpwind_r[15]+hamil[6]*fUpwind_r[13]+hamil[4]*fUpwind_r[11]+hamil[2]*fUpwind_r[6]); 
  Ghat_r[7] = 0.6123724356957944*(hamil[6]*fUpwind_r[14]+hamil[7]*fUpwind_r[10]+hamil[2]*fUpwind_r[7]+fUpwind_r[3]*hamil[4]); 
  Ghat_r[8] = 0.6123724356957944*(hamil[4]*fUpwind_r[12]+hamil[2]*fUpwind_r[8]+fUpwind_r[5]*hamil[7]+fUpwind_r[1]*hamil[6]); 
  Ghat_r[9] = 0.6123724356957944*(hamil[2]*fUpwind_r[9]+fUpwind_r[0]*hamil[7]+fUpwind_r[2]*hamil[6]+fUpwind_r[4]*hamil[4]); 
  Ghat_r[10] = 0.6123724356957944*(hamil[4]*fUpwind_r[14]+hamil[2]*fUpwind_r[10]+fUpwind_r[7]*hamil[7]+fUpwind_r[3]*hamil[6]); 
  Ghat_r[11] = 0.6123724356957944*(hamil[6]*fUpwind_r[15]+hamil[7]*fUpwind_r[13]+hamil[2]*fUpwind_r[11]+hamil[4]*fUpwind_r[6]); 
  Ghat_r[12] = 0.6123724356957944*(hamil[2]*fUpwind_r[12]+hamil[4]*fUpwind_r[8]+fUpwind_r[1]*hamil[7]+fUpwind_r[5]*hamil[6]); 
  Ghat_r[13] = 0.6123724356957944*(hamil[4]*fUpwind_r[15]+hamil[2]*fUpwind_r[13]+hamil[7]*fUpwind_r[11]+fUpwind_r[6]*hamil[6]); 
  Ghat_r[14] = 0.6123724356957944*(hamil[2]*fUpwind_r[14]+hamil[4]*fUpwind_r[10]+fUpwind_r[3]*hamil[7]+hamil[6]*fUpwind_r[7]); 
  Ghat_r[15] = 0.6123724356957944*(hamil[2]*fUpwind_r[15]+hamil[4]*fUpwind_r[13]+hamil[6]*fUpwind_r[11]+fUpwind_r[6]*hamil[7]); 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv11*dx11; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv11*dx11; 
  out[2] += -(1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv11*dx11); 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv11*dx11; 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv11*dx11; 
  out[5] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dv11*dx11; 
  out[6] += -(1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv11*dx11); 
  out[7] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dv11*dx11; 
  out[8] += -(1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv11*dx11); 
  out[9] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dv11*dx11; 
  out[10] += -(1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv11*dx11); 
  out[11] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dv11*dx11; 
  out[12] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dv11*dx11; 
  out[13] += -(1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv11*dx11); 
  out[14] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*dv11*dx11; 
  out[15] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*dv11*dx11; 
  out[16] += -(1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv11*dx11); 
  out[17] += -(1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dv11*dx11); 
  out[18] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*dv11*dx11; 
  out[19] += -(1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dv11*dx11); 
  out[20] += -(1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dv11*dx11); 
  out[21] += (0.7071067811865475*Ghat_l[12]-0.7071067811865475*Ghat_r[12])*dv11*dx11; 
  out[22] += -(1.224744871391589*(Ghat_r[9]+Ghat_l[9])*dv11*dx11); 
  out[23] += (0.7071067811865475*Ghat_l[13]-0.7071067811865475*Ghat_r[13])*dv11*dx11; 
  out[24] += -(1.224744871391589*(Ghat_r[10]+Ghat_l[10])*dv11*dx11); 
  out[25] += (0.7071067811865475*Ghat_l[14]-0.7071067811865475*Ghat_r[14])*dv11*dx11; 
  out[26] += -(1.224744871391589*(Ghat_r[11]+Ghat_l[11])*dv11*dx11); 
  out[27] += -(1.224744871391589*(Ghat_r[12]+Ghat_l[12])*dv11*dx11); 
  out[28] += -(1.224744871391589*(Ghat_r[13]+Ghat_l[13])*dv11*dx11); 
  out[29] += (0.7071067811865475*Ghat_l[15]-0.7071067811865475*Ghat_r[15])*dv11*dx11; 
  out[30] += -(1.224744871391589*(Ghat_r[14]+Ghat_l[14])*dv11*dx11); 
  out[31] += -(1.224744871391589*(Ghat_r[15]+Ghat_l[15])*dv11*dx11); 

  return fabs(0.5303300858899105*dv11*dx11*jacob_cy_c_inv*(1.7320508075688772*hamil[2]*jacob_vy_inv));

} 
