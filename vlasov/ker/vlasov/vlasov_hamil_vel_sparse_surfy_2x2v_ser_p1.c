#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_vel_sparse_surfy_2x2v_ser_p1(const double *w, const double *dxv,
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
  double fUpwind_r[8] = {0.0}; 
  double fUpwind_l[8] = {0.0}; 
  double Ghat_r[8] = {0.0}; 
  double Ghat_l[8] = {0.0}; 
  if (wv>0) { 

  fUpwind_r[0] = (1.224744871391589*fc[2]+0.7071067811865475*fc[0])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[1] = (1.224744871391589*fc[5]+0.7071067811865475*fc[1])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[2] = (1.224744871391589*fc[7]+0.7071067811865475*fc[3])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[3] = (1.224744871391589*fc[9]+0.7071067811865475*fc[4])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[4] = (1.224744871391589*fc[11]+0.7071067811865475*fc[6])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[5] = (1.224744871391589*fc[12]+0.7071067811865475*fc[8])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[6] = (1.224744871391589*fc[14]+0.7071067811865475*fc[10])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_r[7] = (1.224744871391589*fc[15]+0.7071067811865475*fc[13])*jacob_cy_c_inv*jacob_vy_inv; 

  fUpwind_l[0] = (1.224744871391589*fl[2]+0.7071067811865475*fl[0])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[1] = (1.224744871391589*fl[5]+0.7071067811865475*fl[1])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[2] = (1.224744871391589*fl[7]+0.7071067811865475*fl[3])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[3] = (1.224744871391589*fl[9]+0.7071067811865475*fl[4])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[4] = (1.224744871391589*fl[11]+0.7071067811865475*fl[6])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[5] = (1.224744871391589*fl[12]+0.7071067811865475*fl[8])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[6] = (1.224744871391589*fl[14]+0.7071067811865475*fl[10])*jacob_cy_l_inv*jacob_vy_inv; 
  fUpwind_l[7] = (1.224744871391589*fl[15]+0.7071067811865475*fl[13])*jacob_cy_l_inv*jacob_vy_inv; 

  } else { 

  fUpwind_r[0] = (0.7071067811865475*fr[0]-1.224744871391589*fr[2])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[1] = (0.7071067811865475*fr[1]-1.224744871391589*fr[5])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[2] = (0.7071067811865475*fr[3]-1.224744871391589*fr[7])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[3] = (0.7071067811865475*fr[4]-1.224744871391589*fr[9])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[4] = (0.7071067811865475*fr[6]-1.224744871391589*fr[11])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[5] = (0.7071067811865475*fr[8]-1.224744871391589*fr[12])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[6] = (0.7071067811865475*fr[10]-1.224744871391589*fr[14])*jacob_cy_r_inv*jacob_vy_inv; 
  fUpwind_r[7] = (0.7071067811865475*fr[13]-1.224744871391589*fr[15])*jacob_cy_r_inv*jacob_vy_inv; 

  fUpwind_l[0] = (0.7071067811865475*fc[0]-1.224744871391589*fc[2])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[1] = (0.7071067811865475*fc[1]-1.224744871391589*fc[5])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[2] = (0.7071067811865475*fc[3]-1.224744871391589*fc[7])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[3] = (0.7071067811865475*fc[4]-1.224744871391589*fc[9])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[4] = (0.7071067811865475*fc[6]-1.224744871391589*fc[11])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[5] = (0.7071067811865475*fc[8]-1.224744871391589*fc[12])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[6] = (0.7071067811865475*fc[10]-1.224744871391589*fc[14])*jacob_cy_c_inv*jacob_vy_inv; 
  fUpwind_l[7] = (0.7071067811865475*fc[13]-1.224744871391589*fc[15])*jacob_cy_c_inv*jacob_vy_inv; 

  } 
  Ghat_l[0] = 0.8660254037844386*fUpwind_l[0]*hamil[2]; 
  Ghat_l[1] = 0.8660254037844386*fUpwind_l[1]*hamil[2]; 
  Ghat_l[2] = 0.8660254037844386*fUpwind_l[2]*hamil[2]; 
  Ghat_l[3] = 0.8660254037844386*hamil[2]*fUpwind_l[3]; 
  Ghat_l[4] = 0.8660254037844386*hamil[2]*fUpwind_l[4]; 
  Ghat_l[5] = 0.8660254037844386*hamil[2]*fUpwind_l[5]; 
  Ghat_l[6] = 0.8660254037844386*hamil[2]*fUpwind_l[6]; 
  Ghat_l[7] = 0.8660254037844386*hamil[2]*fUpwind_l[7]; 

  Ghat_r[0] = 0.8660254037844386*fUpwind_r[0]*hamil[2]; 
  Ghat_r[1] = 0.8660254037844386*fUpwind_r[1]*hamil[2]; 
  Ghat_r[2] = 0.8660254037844386*fUpwind_r[2]*hamil[2]; 
  Ghat_r[3] = 0.8660254037844386*hamil[2]*fUpwind_r[3]; 
  Ghat_r[4] = 0.8660254037844386*hamil[2]*fUpwind_r[4]; 
  Ghat_r[5] = 0.8660254037844386*hamil[2]*fUpwind_r[5]; 
  Ghat_r[6] = 0.8660254037844386*hamil[2]*fUpwind_r[6]; 
  Ghat_r[7] = 0.8660254037844386*hamil[2]*fUpwind_r[7]; 

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
  out[11] += -(1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv11*dx11); 
  out[12] += -(1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv11*dx11); 
  out[13] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dv11*dx11; 
  out[14] += -(1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dv11*dx11); 
  out[15] += -(1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dv11*dx11); 

  return fabs(0.75*dv11*dx11*jacob_cy_c_inv*(1.7320508075688772*hamil[2]*jacob_vy_inv));

} 
