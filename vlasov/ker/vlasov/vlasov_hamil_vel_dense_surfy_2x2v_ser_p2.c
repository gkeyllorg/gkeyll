#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfy_2x2v_ser_p2(const double *w, const double *dxv,
  const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel,
  const double *poisson_tensor_conf, const double *hamil,
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  double dx11 = 2.0/dxv[1]; 
  double dv11 = 2.0/dxv[3]; 
  double wv = w[3]; 

  const double *jacob_vy = &jacob_vel[3]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  const double *jacob_cy_l = &jacob_pos_l[3]; 
  const double *jacob_cy_c = &jacob_pos_c[3]; 
  const double *jacob_cy_r = &jacob_pos_r[3]; 
  const double jacob_cy_l_inv = 1.0/jacob_cy_l[0]; 
  const double jacob_cy_c_inv = 1.0/jacob_cy_c[0]; 
  const double jacob_cy_r_inv = 1.0/jacob_cy_r[0]; 
  double fUpwind_r[20] = {0.0}; 
  double fUpwind_l[20] = {0.0}; 
  double Ghat_r[20] = {0.0}; 
  double Ghat_l[20] = {0.0}; 
  if (wv>0) { 

  fUpwind_r[0] = (1.5811388300841895*fc[12]+1.224744871391589*fc[2]+0.7071067811865475*fc[0])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[1] = (1.5811388300841898*fc[20]+1.224744871391589*fc[5]+0.7071067811865475*fc[1])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[2] = (1.5811388300841898*fc[22]+1.224744871391589*fc[7]+0.7071067811865475*fc[3])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[3] = (1.5811388300841898*fc[26]+1.224744871391589*fc[9]+0.7071067811865475*fc[4])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[4] = (1.5811388300841895*fc[33]+1.224744871391589*fc[15]+0.7071067811865475*fc[6])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[5] = (1.5811388300841895*fc[36]+1.224744871391589*fc[16]+0.7071067811865475*fc[8])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[6] = (1.5811388300841895*fc[38]+1.224744871391589*fc[18]+0.7071067811865475*fc[10])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[7] = (1.224744871391589*fc[19]+0.7071067811865475*fc[11])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[8] = (1.224744871391589*fc[24]+0.7071067811865475*fc[13])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[9] = (1.224744871391589*fc[29]+0.7071067811865475*fc[14])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[10] = (1.5811388300841898*fc[45]+1.224744871391589*fc[31]+0.7071067811865475*fc[17])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[11] = (1.224744871391589*fc[32]+0.7071067811865475*fc[21])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[12] = (1.224744871391589*fc[34]+0.7071067811865475*fc[23])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[13] = (1.224744871391589*fc[35]+0.7071067811865475*fc[25])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[14] = (1.224744871391589*fc[40]+0.7071067811865475*fc[27])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[15] = (1.224744871391589*fc[41]+0.7071067811865475*fc[28])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[16] = (1.224744871391589*fc[43]+0.7071067811865475*fc[30])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[17] = (1.224744871391589*fc[44]+0.7071067811865475*fc[37])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[18] = (1.224744871391589*fc[46]+0.7071067811865475*fc[39])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[19] = (1.224744871391589*fc[47]+0.7071067811865475*fc[42])*jacob_cy_c_inv*jacob_vy_inv; 

  fUpwind_l[0] = (1.5811388300841895*fl[12]+1.224744871391589*fl[2]+0.7071067811865475*fl[0])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[1] = (1.5811388300841898*fl[20]+1.224744871391589*fl[5]+0.7071067811865475*fl[1])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[2] = (1.5811388300841898*fl[22]+1.224744871391589*fl[7]+0.7071067811865475*fl[3])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[3] = (1.5811388300841898*fl[26]+1.224744871391589*fl[9]+0.7071067811865475*fl[4])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[4] = (1.5811388300841895*fl[33]+1.224744871391589*fl[15]+0.7071067811865475*fl[6])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[5] = (1.5811388300841895*fl[36]+1.224744871391589*fl[16]+0.7071067811865475*fl[8])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[6] = (1.5811388300841895*fl[38]+1.224744871391589*fl[18]+0.7071067811865475*fl[10])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[7] = (1.224744871391589*fl[19]+0.7071067811865475*fl[11])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[8] = (1.224744871391589*fl[24]+0.7071067811865475*fl[13])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[9] = (1.224744871391589*fl[29]+0.7071067811865475*fl[14])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[10] = (1.5811388300841898*fl[45]+1.224744871391589*fl[31]+0.7071067811865475*fl[17])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[11] = (1.224744871391589*fl[32]+0.7071067811865475*fl[21])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[12] = (1.224744871391589*fl[34]+0.7071067811865475*fl[23])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[13] = (1.224744871391589*fl[35]+0.7071067811865475*fl[25])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[14] = (1.224744871391589*fl[40]+0.7071067811865475*fl[27])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[15] = (1.224744871391589*fl[41]+0.7071067811865475*fl[28])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[16] = (1.224744871391589*fl[43]+0.7071067811865475*fl[30])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[17] = (1.224744871391589*fl[44]+0.7071067811865475*fl[37])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[18] = (1.224744871391589*fl[46]+0.7071067811865475*fl[39])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[19] = (1.224744871391589*fl[47]+0.7071067811865475*fl[42])*jacob_cy_l_inv*jacob_vy_inv; 

  } else { 

  fUpwind_r[0] = (1.5811388300841895*fr[12]-1.224744871391589*fr[2]+0.7071067811865475*fr[0])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[1] = (1.5811388300841898*fr[20]-1.224744871391589*fr[5]+0.7071067811865475*fr[1])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[2] = (1.5811388300841898*fr[22]-1.224744871391589*fr[7]+0.7071067811865475*fr[3])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[3] = (1.5811388300841898*fr[26]-1.224744871391589*fr[9]+0.7071067811865475*fr[4])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[4] = (1.5811388300841895*fr[33]-1.224744871391589*fr[15]+0.7071067811865475*fr[6])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[5] = (1.5811388300841895*fr[36]-1.224744871391589*fr[16]+0.7071067811865475*fr[8])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[6] = (1.5811388300841895*fr[38]-1.224744871391589*fr[18]+0.7071067811865475*fr[10])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[7] = (0.7071067811865475*fr[11]-1.224744871391589*fr[19])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[8] = (0.7071067811865475*fr[13]-1.224744871391589*fr[24])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[9] = (0.7071067811865475*fr[14]-1.224744871391589*fr[29])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[10] = (1.5811388300841898*fr[45]-1.224744871391589*fr[31]+0.7071067811865475*fr[17])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[11] = (0.7071067811865475*fr[21]-1.224744871391589*fr[32])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[12] = (0.7071067811865475*fr[23]-1.224744871391589*fr[34])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[13] = (0.7071067811865475*fr[25]-1.224744871391589*fr[35])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[14] = (0.7071067811865475*fr[27]-1.224744871391589*fr[40])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[15] = (0.7071067811865475*fr[28]-1.224744871391589*fr[41])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[16] = (0.7071067811865475*fr[30]-1.224744871391589*fr[43])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[17] = (0.7071067811865475*fr[37]-1.224744871391589*fr[44])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[18] = (0.7071067811865475*fr[39]-1.224744871391589*fr[46])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[19] = (0.7071067811865475*fr[42]-1.224744871391589*fr[47])*jacob_cy_r_inv*jacob_vy_inv; 

  fUpwind_l[0] = (1.5811388300841895*fc[12]-1.224744871391589*fc[2]+0.7071067811865475*fc[0])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[1] = (1.5811388300841898*fc[20]-1.224744871391589*fc[5]+0.7071067811865475*fc[1])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[2] = (1.5811388300841898*fc[22]-1.224744871391589*fc[7]+0.7071067811865475*fc[3])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[3] = (1.5811388300841898*fc[26]-1.224744871391589*fc[9]+0.7071067811865475*fc[4])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[4] = (1.5811388300841895*fc[33]-1.224744871391589*fc[15]+0.7071067811865475*fc[6])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[5] = (1.5811388300841895*fc[36]-1.224744871391589*fc[16]+0.7071067811865475*fc[8])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[6] = (1.5811388300841895*fc[38]-1.224744871391589*fc[18]+0.7071067811865475*fc[10])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[7] = (0.7071067811865475*fc[11]-1.224744871391589*fc[19])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[8] = (0.7071067811865475*fc[13]-1.224744871391589*fc[24])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[9] = (0.7071067811865475*fc[14]-1.224744871391589*fc[29])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[10] = (1.5811388300841898*fc[45]-1.224744871391589*fc[31]+0.7071067811865475*fc[17])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[11] = (0.7071067811865475*fc[21]-1.224744871391589*fc[32])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[12] = (0.7071067811865475*fc[23]-1.224744871391589*fc[34])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[13] = (0.7071067811865475*fc[25]-1.224744871391589*fc[35])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[14] = (0.7071067811865475*fc[27]-1.224744871391589*fc[40])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[15] = (0.7071067811865475*fc[28]-1.224744871391589*fc[41])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[16] = (0.7071067811865475*fc[30]-1.224744871391589*fc[43])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[17] = (0.7071067811865475*fc[37]-1.224744871391589*fc[44])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[18] = (0.7071067811865475*fc[39]-1.224744871391589*fc[46])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[19] = (0.7071067811865475*fc[42]-1.224744871391589*fc[47])*jacob_cy_c_inv*jacob_vy_inv; 

  } 
  Ghat_l[0] = 0.8660254037844387*hamil[6]*fUpwind_l[8]+1.9364916731037085*(fUpwind_l[6]*hamil[7]+fUpwind_l[3]*hamil[5])+0.8660254037844386*(fUpwind_l[2]*hamil[3]+fUpwind_l[0]*hamil[2]); 
  Ghat_l[1] = 0.8660254037844386*hamil[6]*fUpwind_l[12]+1.9364916731037085*(hamil[7]*fUpwind_l[10]+fUpwind_l[5]*hamil[5])+0.8660254037844386*(hamil[3]*fUpwind_l[4]+fUpwind_l[1]*hamil[2]); 
  Ghat_l[2] = 1.7320508075688772*hamil[7]*fUpwind_l[14]+0.7745966692414833*hamil[3]*fUpwind_l[8]+1.9364916731037085*fUpwind_l[3]*hamil[7]+0.7745966692414834*fUpwind_l[2]*hamil[6]+1.9364916731037085*hamil[5]*fUpwind_l[6]+0.8660254037844386*(fUpwind_l[0]*hamil[3]+fUpwind_l[2]*hamil[2]); 
  Ghat_l[3] = 1.7320508075688772*hamil[7]*fUpwind_l[16]+0.8660254037844386*hamil[6]*fUpwind_l[14]+1.7320508075688772*hamil[5]*fUpwind_l[9]+1.9364916731037085*fUpwind_l[2]*hamil[7]+0.8660254037844386*hamil[3]*fUpwind_l[6]+1.9364916731037085*fUpwind_l[0]*hamil[5]+0.8660254037844386*hamil[2]*fUpwind_l[3]; 
  Ghat_l[4] = 1.7320508075688774*hamil[7]*fUpwind_l[18]+0.7745966692414834*hamil[3]*fUpwind_l[12]+1.9364916731037085*(hamil[5]*fUpwind_l[10]+fUpwind_l[5]*hamil[7])+0.7745966692414834*fUpwind_l[4]*hamil[6]+0.8660254037844386*(hamil[2]*fUpwind_l[4]+fUpwind_l[1]*hamil[3]); 
  Ghat_l[5] = 1.7320508075688774*hamil[7]*fUpwind_l[19]+0.8660254037844387*hamil[6]*fUpwind_l[18]+1.7320508075688774*hamil[5]*fUpwind_l[15]+0.8660254037844386*hamil[3]*fUpwind_l[10]+1.9364916731037085*(fUpwind_l[4]*hamil[7]+fUpwind_l[1]*hamil[5])+0.8660254037844386*hamil[2]*fUpwind_l[5]; 
  Ghat_l[6] = 1.7320508075688774*hamil[5]*fUpwind_l[16]+0.7745966692414834*hamil[3]*fUpwind_l[14]+hamil[7]*(1.7320508075688774*(fUpwind_l[9]+fUpwind_l[8])+1.9364916731037085*fUpwind_l[0])+fUpwind_l[6]*(0.7745966692414834*hamil[6]+0.8660254037844386*hamil[2])+1.9364916731037085*fUpwind_l[2]*hamil[5]+0.8660254037844386*fUpwind_l[3]*hamil[3]; 
  Ghat_l[7] = 1.9364916731037085*(hamil[7]*fUpwind_l[17]+hamil[5]*fUpwind_l[13])+0.8660254037844387*hamil[3]*fUpwind_l[11]+0.8660254037844386*hamil[2]*fUpwind_l[7]; 
  Ghat_l[8] = 1.9364916731037085*hamil[5]*fUpwind_l[14]+(0.5532833351724881*hamil[6]+0.8660254037844386*hamil[2])*fUpwind_l[8]+1.7320508075688774*fUpwind_l[6]*hamil[7]+0.8660254037844387*fUpwind_l[0]*hamil[6]+0.7745966692414833*fUpwind_l[2]*hamil[3]; 
  Ghat_l[9] = 0.8660254037844387*hamil[3]*fUpwind_l[16]+0.8660254037844386*hamil[2]*fUpwind_l[9]+1.7320508075688774*fUpwind_l[6]*hamil[7]+1.7320508075688772*fUpwind_l[3]*hamil[5]; 
  Ghat_l[10] = 1.7320508075688772*hamil[5]*fUpwind_l[19]+0.7745966692414833*hamil[3]*fUpwind_l[18]+1.7320508075688772*hamil[7]*(fUpwind_l[15]+fUpwind_l[12])+(0.7745966692414834*hamil[6]+0.8660254037844386*hamil[2])*fUpwind_l[10]+1.9364916731037085*(fUpwind_l[1]*hamil[7]+fUpwind_l[4]*hamil[5])+0.8660254037844386*hamil[3]*fUpwind_l[5]; 
  Ghat_l[11] = 1.9364916731037085*(hamil[5]*fUpwind_l[17]+hamil[7]*fUpwind_l[13])+(0.7745966692414834*hamil[6]+0.8660254037844386*hamil[2])*fUpwind_l[11]+0.8660254037844387*hamil[3]*fUpwind_l[7]; 
  Ghat_l[12] = 1.9364916731037085*hamil[5]*fUpwind_l[18]+(0.5532833351724881*hamil[6]+0.8660254037844386*hamil[2])*fUpwind_l[12]+1.7320508075688772*hamil[7]*fUpwind_l[10]+0.8660254037844386*fUpwind_l[1]*hamil[6]+0.7745966692414834*hamil[3]*fUpwind_l[4]; 
  Ghat_l[13] = 0.8660254037844387*hamil[3]*fUpwind_l[17]+0.8660254037844386*hamil[2]*fUpwind_l[13]+1.9364916731037085*(hamil[7]*fUpwind_l[11]+hamil[5]*fUpwind_l[7]); 
  Ghat_l[14] = 1.5491933384829668*hamil[7]*fUpwind_l[16]+(0.5532833351724881*hamil[6]+0.8660254037844386*hamil[2])*fUpwind_l[14]+1.9364916731037085*hamil[5]*fUpwind_l[8]+1.7320508075688772*fUpwind_l[2]*hamil[7]+0.8660254037844386*fUpwind_l[3]*hamil[6]+0.7745966692414834*hamil[3]*fUpwind_l[6]; 
  Ghat_l[15] = 0.8660254037844387*hamil[3]*fUpwind_l[19]+0.8660254037844386*hamil[2]*fUpwind_l[15]+1.7320508075688772*hamil[7]*fUpwind_l[10]+1.7320508075688774*fUpwind_l[5]*hamil[5]; 
  Ghat_l[16] = (0.7745966692414834*hamil[6]+0.8660254037844386*hamil[2])*fUpwind_l[16]+1.5491933384829668*hamil[7]*fUpwind_l[14]+0.8660254037844387*hamil[3]*fUpwind_l[9]+1.7320508075688772*fUpwind_l[3]*hamil[7]+1.7320508075688774*hamil[5]*fUpwind_l[6]; 
  Ghat_l[17] = (0.7745966692414834*hamil[6]+0.8660254037844386*hamil[2])*fUpwind_l[17]+0.8660254037844387*hamil[3]*fUpwind_l[13]+1.9364916731037085*(hamil[5]*fUpwind_l[11]+fUpwind_l[7]*hamil[7]); 
  Ghat_l[18] = 1.5491933384829668*hamil[7]*fUpwind_l[19]+(0.5532833351724881*hamil[6]+0.8660254037844386*hamil[2])*fUpwind_l[18]+1.9364916731037085*hamil[5]*fUpwind_l[12]+0.7745966692414833*hamil[3]*fUpwind_l[10]+1.7320508075688774*fUpwind_l[4]*hamil[7]+0.8660254037844387*fUpwind_l[5]*hamil[6]; 
  Ghat_l[19] = (0.7745966692414834*hamil[6]+0.8660254037844386*hamil[2])*fUpwind_l[19]+1.5491933384829668*hamil[7]*fUpwind_l[18]+0.8660254037844387*hamil[3]*fUpwind_l[15]+1.7320508075688772*hamil[5]*fUpwind_l[10]+1.7320508075688774*fUpwind_l[5]*hamil[7]; 

  Ghat_r[0] = 0.8660254037844387*hamil[6]*fUpwind_r[8]+1.9364916731037085*(fUpwind_r[6]*hamil[7]+fUpwind_r[3]*hamil[5])+0.8660254037844386*(fUpwind_r[2]*hamil[3]+fUpwind_r[0]*hamil[2]); 
  Ghat_r[1] = 0.8660254037844386*hamil[6]*fUpwind_r[12]+1.9364916731037085*(hamil[7]*fUpwind_r[10]+fUpwind_r[5]*hamil[5])+0.8660254037844386*(hamil[3]*fUpwind_r[4]+fUpwind_r[1]*hamil[2]); 
  Ghat_r[2] = 1.7320508075688772*hamil[7]*fUpwind_r[14]+0.7745966692414833*hamil[3]*fUpwind_r[8]+1.9364916731037085*fUpwind_r[3]*hamil[7]+0.7745966692414834*fUpwind_r[2]*hamil[6]+1.9364916731037085*hamil[5]*fUpwind_r[6]+0.8660254037844386*(fUpwind_r[0]*hamil[3]+fUpwind_r[2]*hamil[2]); 
  Ghat_r[3] = 1.7320508075688772*hamil[7]*fUpwind_r[16]+0.8660254037844386*hamil[6]*fUpwind_r[14]+1.7320508075688772*hamil[5]*fUpwind_r[9]+1.9364916731037085*fUpwind_r[2]*hamil[7]+0.8660254037844386*hamil[3]*fUpwind_r[6]+1.9364916731037085*fUpwind_r[0]*hamil[5]+0.8660254037844386*hamil[2]*fUpwind_r[3]; 
  Ghat_r[4] = 1.7320508075688774*hamil[7]*fUpwind_r[18]+0.7745966692414834*hamil[3]*fUpwind_r[12]+1.9364916731037085*(hamil[5]*fUpwind_r[10]+fUpwind_r[5]*hamil[7])+0.7745966692414834*fUpwind_r[4]*hamil[6]+0.8660254037844386*(hamil[2]*fUpwind_r[4]+fUpwind_r[1]*hamil[3]); 
  Ghat_r[5] = 1.7320508075688774*hamil[7]*fUpwind_r[19]+0.8660254037844387*hamil[6]*fUpwind_r[18]+1.7320508075688774*hamil[5]*fUpwind_r[15]+0.8660254037844386*hamil[3]*fUpwind_r[10]+1.9364916731037085*(fUpwind_r[4]*hamil[7]+fUpwind_r[1]*hamil[5])+0.8660254037844386*hamil[2]*fUpwind_r[5]; 
  Ghat_r[6] = 1.7320508075688774*hamil[5]*fUpwind_r[16]+0.7745966692414834*hamil[3]*fUpwind_r[14]+hamil[7]*(1.7320508075688774*(fUpwind_r[9]+fUpwind_r[8])+1.9364916731037085*fUpwind_r[0])+fUpwind_r[6]*(0.7745966692414834*hamil[6]+0.8660254037844386*hamil[2])+1.9364916731037085*fUpwind_r[2]*hamil[5]+0.8660254037844386*fUpwind_r[3]*hamil[3]; 
  Ghat_r[7] = 1.9364916731037085*(hamil[7]*fUpwind_r[17]+hamil[5]*fUpwind_r[13])+0.8660254037844387*hamil[3]*fUpwind_r[11]+0.8660254037844386*hamil[2]*fUpwind_r[7]; 
  Ghat_r[8] = 1.9364916731037085*hamil[5]*fUpwind_r[14]+(0.5532833351724881*hamil[6]+0.8660254037844386*hamil[2])*fUpwind_r[8]+1.7320508075688774*fUpwind_r[6]*hamil[7]+0.8660254037844387*fUpwind_r[0]*hamil[6]+0.7745966692414833*fUpwind_r[2]*hamil[3]; 
  Ghat_r[9] = 0.8660254037844387*hamil[3]*fUpwind_r[16]+0.8660254037844386*hamil[2]*fUpwind_r[9]+1.7320508075688774*fUpwind_r[6]*hamil[7]+1.7320508075688772*fUpwind_r[3]*hamil[5]; 
  Ghat_r[10] = 1.7320508075688772*hamil[5]*fUpwind_r[19]+0.7745966692414833*hamil[3]*fUpwind_r[18]+1.7320508075688772*hamil[7]*(fUpwind_r[15]+fUpwind_r[12])+(0.7745966692414834*hamil[6]+0.8660254037844386*hamil[2])*fUpwind_r[10]+1.9364916731037085*(fUpwind_r[1]*hamil[7]+fUpwind_r[4]*hamil[5])+0.8660254037844386*hamil[3]*fUpwind_r[5]; 
  Ghat_r[11] = 1.9364916731037085*(hamil[5]*fUpwind_r[17]+hamil[7]*fUpwind_r[13])+(0.7745966692414834*hamil[6]+0.8660254037844386*hamil[2])*fUpwind_r[11]+0.8660254037844387*hamil[3]*fUpwind_r[7]; 
  Ghat_r[12] = 1.9364916731037085*hamil[5]*fUpwind_r[18]+(0.5532833351724881*hamil[6]+0.8660254037844386*hamil[2])*fUpwind_r[12]+1.7320508075688772*hamil[7]*fUpwind_r[10]+0.8660254037844386*fUpwind_r[1]*hamil[6]+0.7745966692414834*hamil[3]*fUpwind_r[4]; 
  Ghat_r[13] = 0.8660254037844387*hamil[3]*fUpwind_r[17]+0.8660254037844386*hamil[2]*fUpwind_r[13]+1.9364916731037085*(hamil[7]*fUpwind_r[11]+hamil[5]*fUpwind_r[7]); 
  Ghat_r[14] = 1.5491933384829668*hamil[7]*fUpwind_r[16]+(0.5532833351724881*hamil[6]+0.8660254037844386*hamil[2])*fUpwind_r[14]+1.9364916731037085*hamil[5]*fUpwind_r[8]+1.7320508075688772*fUpwind_r[2]*hamil[7]+0.8660254037844386*fUpwind_r[3]*hamil[6]+0.7745966692414834*hamil[3]*fUpwind_r[6]; 
  Ghat_r[15] = 0.8660254037844387*hamil[3]*fUpwind_r[19]+0.8660254037844386*hamil[2]*fUpwind_r[15]+1.7320508075688772*hamil[7]*fUpwind_r[10]+1.7320508075688774*fUpwind_r[5]*hamil[5]; 
  Ghat_r[16] = (0.7745966692414834*hamil[6]+0.8660254037844386*hamil[2])*fUpwind_r[16]+1.5491933384829668*hamil[7]*fUpwind_r[14]+0.8660254037844387*hamil[3]*fUpwind_r[9]+1.7320508075688772*fUpwind_r[3]*hamil[7]+1.7320508075688774*hamil[5]*fUpwind_r[6]; 
  Ghat_r[17] = (0.7745966692414834*hamil[6]+0.8660254037844386*hamil[2])*fUpwind_r[17]+0.8660254037844387*hamil[3]*fUpwind_r[13]+1.9364916731037085*(hamil[5]*fUpwind_r[11]+fUpwind_r[7]*hamil[7]); 
  Ghat_r[18] = 1.5491933384829668*hamil[7]*fUpwind_r[19]+(0.5532833351724881*hamil[6]+0.8660254037844386*hamil[2])*fUpwind_r[18]+1.9364916731037085*hamil[5]*fUpwind_r[12]+0.7745966692414833*hamil[3]*fUpwind_r[10]+1.7320508075688774*fUpwind_r[4]*hamil[7]+0.8660254037844387*fUpwind_r[5]*hamil[6]; 
  Ghat_r[19] = (0.7745966692414834*hamil[6]+0.8660254037844386*hamil[2])*fUpwind_r[19]+1.5491933384829668*hamil[7]*fUpwind_r[18]+0.8660254037844387*hamil[3]*fUpwind_r[15]+1.7320508075688772*hamil[5]*fUpwind_r[10]+1.7320508075688774*fUpwind_r[5]*hamil[7]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv11*dx11; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv11*dx11; 
  out[2] += -(1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv11*dx11); 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv11*dx11; 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv11*dx11; 
  out[5] += -(1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv11*dx11); 
  out[6] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dv11*dx11; 
  out[7] += -(1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv11*dx11); 
  out[8] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dv11*dx11; 
  out[9] += -(1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv11*dx11); 
  out[10] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dv11*dx11; 
  out[11] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dv11*dx11; 
  out[12] += (1.5811388300841895*Ghat_l[0]-1.5811388300841895*Ghat_r[0])*dv11*dx11; 
  out[13] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dv11*dx11; 
  out[14] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*dv11*dx11; 
  out[15] += -(1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv11*dx11); 
  out[16] += -(1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv11*dx11); 
  out[17] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*dv11*dx11; 
  out[18] += -(1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dv11*dx11); 
  out[19] += -(1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dv11*dx11); 
  out[20] += (1.5811388300841898*Ghat_l[1]-1.5811388300841898*Ghat_r[1])*dv11*dx11; 
  out[21] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*dv11*dx11; 
  out[22] += (1.5811388300841898*Ghat_l[2]-1.5811388300841898*Ghat_r[2])*dv11*dx11; 
  out[23] += (0.7071067811865475*Ghat_l[12]-0.7071067811865475*Ghat_r[12])*dv11*dx11; 
  out[24] += -(1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dv11*dx11); 
  out[25] += (0.7071067811865475*Ghat_l[13]-0.7071067811865475*Ghat_r[13])*dv11*dx11; 
  out[26] += (1.5811388300841898*Ghat_l[3]-1.5811388300841898*Ghat_r[3])*dv11*dx11; 
  out[27] += (0.7071067811865475*Ghat_l[14]-0.7071067811865475*Ghat_r[14])*dv11*dx11; 
  out[28] += (0.7071067811865475*Ghat_l[15]-0.7071067811865475*Ghat_r[15])*dv11*dx11; 
  out[29] += -(1.224744871391589*(Ghat_r[9]+Ghat_l[9])*dv11*dx11); 
  out[30] += (0.7071067811865475*Ghat_l[16]-0.7071067811865475*Ghat_r[16])*dv11*dx11; 
  out[31] += -(1.224744871391589*(Ghat_r[10]+Ghat_l[10])*dv11*dx11); 
  out[32] += -(1.224744871391589*(Ghat_r[11]+Ghat_l[11])*dv11*dx11); 
  out[33] += (1.5811388300841895*Ghat_l[4]-1.5811388300841895*Ghat_r[4])*dv11*dx11; 
  out[34] += -(1.224744871391589*(Ghat_r[12]+Ghat_l[12])*dv11*dx11); 
  out[35] += -(1.224744871391589*(Ghat_r[13]+Ghat_l[13])*dv11*dx11); 
  out[36] += (1.5811388300841895*Ghat_l[5]-1.5811388300841895*Ghat_r[5])*dv11*dx11; 
  out[37] += (0.7071067811865475*Ghat_l[17]-0.7071067811865475*Ghat_r[17])*dv11*dx11; 
  out[38] += (1.5811388300841895*Ghat_l[6]-1.5811388300841895*Ghat_r[6])*dv11*dx11; 
  out[39] += (0.7071067811865475*Ghat_l[18]-0.7071067811865475*Ghat_r[18])*dv11*dx11; 
  out[40] += -(1.224744871391589*(Ghat_r[14]+Ghat_l[14])*dv11*dx11); 
  out[41] += -(1.224744871391589*(Ghat_r[15]+Ghat_l[15])*dv11*dx11); 
  out[42] += (0.7071067811865475*Ghat_l[19]-0.7071067811865475*Ghat_r[19])*dv11*dx11; 
  out[43] += -(1.224744871391589*(Ghat_r[16]+Ghat_l[16])*dv11*dx11); 
  out[44] += -(1.224744871391589*(Ghat_r[17]+Ghat_l[17])*dv11*dx11); 
  out[45] += (1.5811388300841898*Ghat_l[10]-1.5811388300841898*Ghat_r[10])*dv11*dx11; 
  out[46] += -(1.224744871391589*(Ghat_r[18]+Ghat_l[18])*dv11*dx11); 
  out[47] += -(1.224744871391589*(Ghat_r[19]+Ghat_l[19])*dv11*dx11); 

  return fabs(1.25*dv11*dx11*jacob_cy_c_inv*(1.7320508075688772*hamil[2]*jacob_vy_inv));

} 
