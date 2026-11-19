#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfx_2x1v_ser_p1(const double *w, const double *dxv,
  const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel,
  const double *poisson_tensor_conf, const double *hamil,
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[2]; 
  double wv = w[2]; 

  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  const double *jacob_cx_l = &jacob_pos_l[0]; 
  const double *jacob_cx_c = &jacob_pos_c[0]; 
  const double *jacob_cx_r = &jacob_pos_r[0]; 
  const double jacob_cx_l_inv = 1.0/jacob_cx_l[0]; 
  const double jacob_cx_c_inv = 1.0/jacob_cx_c[0]; 
  const double jacob_cx_r_inv = 1.0/jacob_cx_r[0]; 
  double fUpwind_r[4] = {0.0}; 
  double fUpwind_l[4] = {0.0}; 
  double Ghat_r[4] = {0.0}; 
  double Ghat_l[4] = {0.0}; 
  if (wv>0) { 

  fUpwind_r[0] = (1.224744871391589*fc[1]+0.7071067811865475*fc[0])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[1] = (1.224744871391589*fc[4]+0.7071067811865475*fc[2])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[2] = (1.224744871391589*fc[5]+0.7071067811865475*fc[3])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[3] = (1.224744871391589*fc[7]+0.7071067811865475*fc[6])*jacob_cx_c_inv*jacob_vx_inv; 

  fUpwind_l[0] = (1.224744871391589*fl[1]+0.7071067811865475*fl[0])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[1] = (1.224744871391589*fl[4]+0.7071067811865475*fl[2])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[2] = (1.224744871391589*fl[5]+0.7071067811865475*fl[3])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[3] = (1.224744871391589*fl[7]+0.7071067811865475*fl[6])*jacob_cx_l_inv*jacob_vx_inv; 

  } else { 

  fUpwind_r[0] = (0.7071067811865475*fr[0]-1.224744871391589*fr[1])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[1] = (0.7071067811865475*fr[2]-1.224744871391589*fr[4])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[2] = (0.7071067811865475*fr[3]-1.224744871391589*fr[5])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[3] = (0.7071067811865475*fr[6]-1.224744871391589*fr[7])*jacob_cx_r_inv*jacob_vx_inv; 

  fUpwind_l[0] = (0.7071067811865475*fc[0]-1.224744871391589*fc[1])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[1] = (0.7071067811865475*fc[2]-1.224744871391589*fc[4])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[2] = (0.7071067811865475*fc[3]-1.224744871391589*fc[5])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[3] = (0.7071067811865475*fc[6]-1.224744871391589*fc[7])*jacob_cx_c_inv*jacob_vx_inv; 

  } 
  Ghat_l[0] = 1.224744871391589*fUpwind_l[0]*hamil[1]; 
  Ghat_l[1] = 1.224744871391589*fUpwind_l[1]*hamil[1]; 
  Ghat_l[2] = 1.224744871391589*hamil[1]*fUpwind_l[2]; 
  Ghat_l[3] = 1.224744871391589*hamil[1]*fUpwind_l[3]; 

  Ghat_r[0] = 1.224744871391589*fUpwind_r[0]*hamil[1]; 
  Ghat_r[1] = 1.224744871391589*fUpwind_r[1]*hamil[1]; 
  Ghat_r[2] = 1.224744871391589*hamil[1]*fUpwind_r[2]; 
  Ghat_r[3] = 1.224744871391589*hamil[1]*fUpwind_r[3]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv10*dx10; 
  out[1] += -(1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv10*dx10); 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv10*dx10; 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv10*dx10; 
  out[4] += -(1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv10*dx10); 
  out[5] += -(1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv10*dx10); 
  out[6] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv10*dx10; 
  out[7] += -(1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv10*dx10); 

  return fabs(1.0606601717798212*dv10*dx10*jacob_cx_c_inv*(1.7320508075688772*hamil[1]*jacob_vx_inv));

} 
