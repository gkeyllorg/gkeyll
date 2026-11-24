#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_vel_sparse_surfx_1x2v_ser_p2(const double *w, const double *dxv,
  const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel,
  const double *poisson_tensor_conf, const double *hamil,
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[1]; 
  double wv = w[1]; 

  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  const double *jacob_cx_l = &jacob_pos_l[0]; 
  const double *jacob_cx_c = &jacob_pos_c[0]; 
  const double *jacob_cx_r = &jacob_pos_r[0]; 
  const double jacob_cx_l_inv = 1.0/jacob_cx_l[0]; 
  const double jacob_cx_c_inv = 1.0/jacob_cx_c[0]; 
  const double jacob_cx_r_inv = 1.0/jacob_cx_r[0]; 
  double fUpwind_r[8] = {0.0}; 
  double fUpwind_l[8] = {0.0}; 
  double Ghat_r[8] = {0.0}; 
  double Ghat_l[8] = {0.0}; 
  if (wv>0) { 

  fUpwind_r[0] = (1.5811388300841895*fc[7]+1.224744871391589*fc[1]+0.7071067811865475*fc[0])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[1] = (1.5811388300841898*fc[11]+1.224744871391589*fc[4]+0.7071067811865475*fc[2])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[2] = (1.5811388300841898*fc[13]+1.224744871391589*fc[5]+0.7071067811865475*fc[3])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[3] = (1.5811388300841895*fc[17]+1.224744871391589*fc[10]+0.7071067811865475*fc[6])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[4] = (1.224744871391589*fc[12]+0.7071067811865475*fc[8])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[5] = (1.224744871391589*fc[15]+0.7071067811865475*fc[9])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[6] = (1.224744871391589*fc[18]+0.7071067811865475*fc[14])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[7] = (1.224744871391589*fc[19]+0.7071067811865475*fc[16])*jacob_cx_c_inv*jacob_vx_inv; 

  fUpwind_l[0] = (1.5811388300841895*fl[7]+1.224744871391589*fl[1]+0.7071067811865475*fl[0])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[1] = (1.5811388300841898*fl[11]+1.224744871391589*fl[4]+0.7071067811865475*fl[2])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[2] = (1.5811388300841898*fl[13]+1.224744871391589*fl[5]+0.7071067811865475*fl[3])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[3] = (1.5811388300841895*fl[17]+1.224744871391589*fl[10]+0.7071067811865475*fl[6])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[4] = (1.224744871391589*fl[12]+0.7071067811865475*fl[8])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[5] = (1.224744871391589*fl[15]+0.7071067811865475*fl[9])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[6] = (1.224744871391589*fl[18]+0.7071067811865475*fl[14])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[7] = (1.224744871391589*fl[19]+0.7071067811865475*fl[16])*jacob_cx_l_inv*jacob_vx_inv; 

  } else { 

  fUpwind_r[0] = (1.5811388300841895*fr[7]-1.224744871391589*fr[1]+0.7071067811865475*fr[0])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[1] = (1.5811388300841898*fr[11]-1.224744871391589*fr[4]+0.7071067811865475*fr[2])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[2] = (1.5811388300841898*fr[13]-1.224744871391589*fr[5]+0.7071067811865475*fr[3])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[3] = (1.5811388300841895*fr[17]-1.224744871391589*fr[10]+0.7071067811865475*fr[6])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[4] = (0.7071067811865475*fr[8]-1.224744871391589*fr[12])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[5] = (0.7071067811865475*fr[9]-1.224744871391589*fr[15])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[6] = (0.7071067811865475*fr[14]-1.224744871391589*fr[18])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[7] = (0.7071067811865475*fr[16]-1.224744871391589*fr[19])*jacob_cx_r_inv*jacob_vx_inv; 

  fUpwind_l[0] = (1.5811388300841895*fc[7]-1.224744871391589*fc[1]+0.7071067811865475*fc[0])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[1] = (1.5811388300841898*fc[11]-1.224744871391589*fc[4]+0.7071067811865475*fc[2])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[2] = (1.5811388300841898*fc[13]-1.224744871391589*fc[5]+0.7071067811865475*fc[3])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[3] = (1.5811388300841895*fc[17]-1.224744871391589*fc[10]+0.7071067811865475*fc[6])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[4] = (0.7071067811865475*fc[8]-1.224744871391589*fc[12])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[5] = (0.7071067811865475*fc[9]-1.224744871391589*fc[15])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[6] = (0.7071067811865475*fc[14]-1.224744871391589*fc[18])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[7] = (0.7071067811865475*fc[16]-1.224744871391589*fc[19])*jacob_cx_c_inv*jacob_vx_inv; 

  } 
  Ghat_l[0] = 1.9364916731037085*fUpwind_l[1]*hamil[4]+0.8660254037844386*fUpwind_l[0]*hamil[1]; 
  Ghat_l[1] = (1.7320508075688772*fUpwind_l[4]+1.9364916731037085*fUpwind_l[0])*hamil[4]+0.8660254037844386*fUpwind_l[1]*hamil[1]; 
  Ghat_l[2] = 1.9364916731037085*fUpwind_l[3]*hamil[4]+0.8660254037844386*hamil[1]*fUpwind_l[2]; 
  Ghat_l[3] = hamil[4]*(1.7320508075688774*fUpwind_l[6]+1.9364916731037085*fUpwind_l[2])+0.8660254037844386*hamil[1]*fUpwind_l[3]; 
  Ghat_l[4] = 1.7320508075688772*fUpwind_l[1]*hamil[4]+0.8660254037844386*hamil[1]*fUpwind_l[4]; 
  Ghat_l[5] = 1.9364916731037085*hamil[4]*fUpwind_l[7]+0.8660254037844386*hamil[1]*fUpwind_l[5]; 
  Ghat_l[6] = 0.8660254037844386*hamil[1]*fUpwind_l[6]+1.7320508075688774*fUpwind_l[3]*hamil[4]; 
  Ghat_l[7] = 0.8660254037844386*hamil[1]*fUpwind_l[7]+1.9364916731037085*hamil[4]*fUpwind_l[5]; 

  Ghat_r[0] = 1.9364916731037085*fUpwind_r[1]*hamil[4]+0.8660254037844386*fUpwind_r[0]*hamil[1]; 
  Ghat_r[1] = (1.7320508075688772*fUpwind_r[4]+1.9364916731037085*fUpwind_r[0])*hamil[4]+0.8660254037844386*fUpwind_r[1]*hamil[1]; 
  Ghat_r[2] = 1.9364916731037085*fUpwind_r[3]*hamil[4]+0.8660254037844386*hamil[1]*fUpwind_r[2]; 
  Ghat_r[3] = hamil[4]*(1.7320508075688774*fUpwind_r[6]+1.9364916731037085*fUpwind_r[2])+0.8660254037844386*hamil[1]*fUpwind_r[3]; 
  Ghat_r[4] = 1.7320508075688772*fUpwind_r[1]*hamil[4]+0.8660254037844386*hamil[1]*fUpwind_r[4]; 
  Ghat_r[5] = 1.9364916731037085*hamil[4]*fUpwind_r[7]+0.8660254037844386*hamil[1]*fUpwind_r[5]; 
  Ghat_r[6] = 0.8660254037844386*hamil[1]*fUpwind_r[6]+1.7320508075688774*fUpwind_r[3]*hamil[4]; 
  Ghat_r[7] = 0.8660254037844386*hamil[1]*fUpwind_r[7]+1.9364916731037085*hamil[4]*fUpwind_r[5]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv10*dx10; 
  out[1] += -(1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv10*dx10); 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv10*dx10; 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv10*dx10; 
  out[4] += -(1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv10*dx10); 
  out[5] += -(1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv10*dx10); 
  out[6] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv10*dx10; 
  out[7] += (1.5811388300841895*Ghat_l[0]-1.5811388300841895*Ghat_r[0])*dv10*dx10; 
  out[8] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dv10*dx10; 
  out[9] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dv10*dx10; 
  out[10] += -(1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv10*dx10); 
  out[11] += (1.5811388300841898*Ghat_l[1]-1.5811388300841898*Ghat_r[1])*dv10*dx10; 
  out[12] += -(1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv10*dx10); 
  out[13] += (1.5811388300841898*Ghat_l[2]-1.5811388300841898*Ghat_r[2])*dv10*dx10; 
  out[14] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dv10*dx10; 
  out[15] += -(1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv10*dx10); 
  out[16] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dv10*dx10; 
  out[17] += (1.5811388300841895*Ghat_l[3]-1.5811388300841895*Ghat_r[3])*dv10*dx10; 
  out[18] += -(1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dv10*dx10); 
  out[19] += -(1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dv10*dx10); 

  return fabs(1.25*dv10*dx10*jacob_cx_c_inv*(1.7320508075688772*hamil[1]*jacob_vx_inv));

} 
