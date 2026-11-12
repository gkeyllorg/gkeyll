#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfx_1x1v_ser_p3(const double *w, const double *dxv,
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
  double fUpwind_r[4] = {0.0}; 
  double fUpwind_l[4] = {0.0}; 
  double Ghat_r[4] = {0.0}; 
  double Ghat_l[4] = {0.0}; 
  if (wv>0) { 

  fUpwind_r[0] = (1.8708286933869707*fc[8]+1.5811388300841895*fc[4]+1.224744871391589*fc[1]+0.7071067811865475*fc[0])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[1] = (1.8708286933869707*fc[10]+1.5811388300841898*fc[6]+1.224744871391589*fc[3]+0.7071067811865475*fc[2])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[2] = (1.224744871391589*fc[7]+0.7071067811865475*fc[5])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[3] = (1.2247448713915887*fc[11]+0.7071067811865475*fc[9])*jacob_cx_c_inv*jacob_vx_inv; 

  fUpwind_l[0] = (1.8708286933869707*fl[8]+1.5811388300841895*fl[4]+1.224744871391589*fl[1]+0.7071067811865475*fl[0])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[1] = (1.8708286933869707*fl[10]+1.5811388300841898*fl[6]+1.224744871391589*fl[3]+0.7071067811865475*fl[2])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[2] = (1.224744871391589*fl[7]+0.7071067811865475*fl[5])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[3] = (1.2247448713915887*fl[11]+0.7071067811865475*fl[9])*jacob_cx_l_inv*jacob_vx_inv; 

  } else { 

  fUpwind_r[0] = (-(1.8708286933869707*fr[8])+1.5811388300841895*fr[4]-1.224744871391589*fr[1]+0.7071067811865475*fr[0])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[1] = (-(1.8708286933869707*fr[10])+1.5811388300841898*fr[6]-1.224744871391589*fr[3]+0.7071067811865475*fr[2])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[2] = (0.7071067811865475*fr[5]-1.224744871391589*fr[7])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[3] = (0.7071067811865475*fr[9]-1.2247448713915887*fr[11])*jacob_cx_r_inv*jacob_vx_inv; 

  fUpwind_l[0] = (-(1.8708286933869707*fc[8])+1.5811388300841895*fc[4]-1.224744871391589*fc[1]+0.7071067811865475*fc[0])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[1] = (-(1.8708286933869707*fc[10])+1.5811388300841898*fc[6]-1.224744871391589*fc[3]+0.7071067811865475*fc[2])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[2] = (0.7071067811865475*fc[5]-1.224744871391589*fc[7])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[3] = (0.7071067811865475*fc[9]-1.2247448713915887*fc[11])*jacob_cx_c_inv*jacob_vx_inv; 

  } 
  Ghat_l[0] = (4.183300132670378*fUpwind_l[2]+1.8708286933869707*fUpwind_l[0])*hamil[3]+2.7386127875258306*fUpwind_l[1]*hamil[2]+1.224744871391589*fUpwind_l[0]*hamil[1]; 
  Ghat_l[1] = (3.6742346141747664*fUpwind_l[3]+5.612486080160912*fUpwind_l[1])*hamil[3]+(2.4494897427831783*fUpwind_l[2]+2.7386127875258306*fUpwind_l[0])*hamil[2]+1.224744871391589*fUpwind_l[1]*hamil[1]; 
  Ghat_l[2] = (4.543441112511214*fUpwind_l[2]+4.183300132670378*fUpwind_l[0])*hamil[3]+hamil[2]*(2.405351177211819*fUpwind_l[3]+2.4494897427831783*fUpwind_l[1])+1.224744871391589*hamil[1]*fUpwind_l[2]; 
  Ghat_l[3] = (4.365266951236265*fUpwind_l[3]+3.6742346141747664*fUpwind_l[1])*hamil[3]+1.224744871391589*hamil[1]*fUpwind_l[3]+2.405351177211819*fUpwind_l[2]*hamil[2]; 

  Ghat_r[0] = (4.183300132670378*fUpwind_r[2]+1.8708286933869707*fUpwind_r[0])*hamil[3]+2.7386127875258306*fUpwind_r[1]*hamil[2]+1.224744871391589*fUpwind_r[0]*hamil[1]; 
  Ghat_r[1] = (3.6742346141747664*fUpwind_r[3]+5.612486080160912*fUpwind_r[1])*hamil[3]+(2.4494897427831783*fUpwind_r[2]+2.7386127875258306*fUpwind_r[0])*hamil[2]+1.224744871391589*fUpwind_r[1]*hamil[1]; 
  Ghat_r[2] = (4.543441112511214*fUpwind_r[2]+4.183300132670378*fUpwind_r[0])*hamil[3]+hamil[2]*(2.405351177211819*fUpwind_r[3]+2.4494897427831783*fUpwind_r[1])+1.224744871391589*hamil[1]*fUpwind_r[2]; 
  Ghat_r[3] = (4.365266951236265*fUpwind_r[3]+3.6742346141747664*fUpwind_r[1])*hamil[3]+1.224744871391589*hamil[1]*fUpwind_r[3]+2.405351177211819*fUpwind_r[2]*hamil[2]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv10*dx10; 
  out[1] += -(1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv10*dx10); 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv10*dx10; 
  out[3] += -(1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv10*dx10); 
  out[4] += (1.5811388300841895*Ghat_l[0]-1.5811388300841895*Ghat_r[0])*dv10*dx10; 
  out[5] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv10*dx10; 
  out[6] += (1.5811388300841898*Ghat_l[1]-1.5811388300841898*Ghat_r[1])*dv10*dx10; 
  out[7] += -(1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv10*dx10); 
  out[8] += -(1.8708286933869707*(Ghat_r[0]+Ghat_l[0])*dv10*dx10); 
  out[9] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv10*dx10; 
  out[10] += -(1.8708286933869707*(Ghat_r[1]+Ghat_l[1])*dv10*dx10); 
  out[11] += -(1.2247448713915887*(Ghat_r[3]+Ghat_l[3])*dv10*dx10); 

  return fabs(2.474873734152916*dv10*dx10*jacob_cx_c_inv*(2.6457513110645907*hamil[3]*jacob_vx_inv+1.7320508075688772*hamil[1]*jacob_vx_inv));

} 
