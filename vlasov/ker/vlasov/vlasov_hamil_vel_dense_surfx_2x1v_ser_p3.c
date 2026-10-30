#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_vel_dense_surfx_2x1v_ser_p3(const double *w, const double *dxv,
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
  double fUpwind_r[12] = {0.0}; 
  double fUpwind_l[12] = {0.0}; 
  double Ghat_r[12] = {0.0}; 
  double Ghat_l[12] = {0.0}; 
  if (wv>0) { 

  fUpwind_r[0] = (1.8708286933869707*fc[17]+1.5811388300841895*fc[7]+1.224744871391589*fc[1]+0.7071067811865475*fc[0])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[1] = (1.8708286933869707*fc[23]+1.5811388300841898*fc[11]+1.224744871391589*fc[4]+0.7071067811865475*fc[2])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[2] = (1.8708286933869707*fc[25]+1.5811388300841898*fc[13]+1.224744871391589*fc[5]+0.7071067811865475*fc[3])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[3] = (1.8708286933869707*fc[29]+1.5811388300841895*fc[20]+1.224744871391589*fc[10]+0.7071067811865475*fc[6])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[4] = (1.224744871391589*fc[12]+0.7071067811865475*fc[8])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[5] = (1.224744871391589*fc[15]+0.7071067811865475*fc[9])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[6] = (1.224744871391589*fc[21]+0.7071067811865475*fc[14])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[7] = (1.224744871391589*fc[22]+0.7071067811865475*fc[16])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[8] = (1.2247448713915887*fc[24]+0.7071067811865475*fc[18])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[9] = (1.2247448713915887*fc[27]+0.7071067811865475*fc[19])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[10] = (1.2247448713915887*fc[30]+0.7071067811865475*fc[26])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_r[11] = (1.2247448713915887*fc[31]+0.7071067811865475*fc[28])*jacob_cx_c_inv*jacob_vx_inv; 

  fUpwind_l[0] = (1.8708286933869707*fl[17]+1.5811388300841895*fl[7]+1.224744871391589*fl[1]+0.7071067811865475*fl[0])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[1] = (1.8708286933869707*fl[23]+1.5811388300841898*fl[11]+1.224744871391589*fl[4]+0.7071067811865475*fl[2])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[2] = (1.8708286933869707*fl[25]+1.5811388300841898*fl[13]+1.224744871391589*fl[5]+0.7071067811865475*fl[3])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[3] = (1.8708286933869707*fl[29]+1.5811388300841895*fl[20]+1.224744871391589*fl[10]+0.7071067811865475*fl[6])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[4] = (1.224744871391589*fl[12]+0.7071067811865475*fl[8])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[5] = (1.224744871391589*fl[15]+0.7071067811865475*fl[9])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[6] = (1.224744871391589*fl[21]+0.7071067811865475*fl[14])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[7] = (1.224744871391589*fl[22]+0.7071067811865475*fl[16])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[8] = (1.2247448713915887*fl[24]+0.7071067811865475*fl[18])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[9] = (1.2247448713915887*fl[27]+0.7071067811865475*fl[19])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[10] = (1.2247448713915887*fl[30]+0.7071067811865475*fl[26])*jacob_cx_l_inv*jacob_vx_inv; 
  fUpwind_l[11] = (1.2247448713915887*fl[31]+0.7071067811865475*fl[28])*jacob_cx_l_inv*jacob_vx_inv; 

  } else { 

  fUpwind_r[0] = (-(1.8708286933869707*fr[17])+1.5811388300841895*fr[7]-1.224744871391589*fr[1]+0.7071067811865475*fr[0])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[1] = (-(1.8708286933869707*fr[23])+1.5811388300841898*fr[11]-1.224744871391589*fr[4]+0.7071067811865475*fr[2])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[2] = (-(1.8708286933869707*fr[25])+1.5811388300841898*fr[13]-1.224744871391589*fr[5]+0.7071067811865475*fr[3])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[3] = (-(1.8708286933869707*fr[29])+1.5811388300841895*fr[20]-1.224744871391589*fr[10]+0.7071067811865475*fr[6])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[4] = (0.7071067811865475*fr[8]-1.224744871391589*fr[12])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[5] = (0.7071067811865475*fr[9]-1.224744871391589*fr[15])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[6] = (0.7071067811865475*fr[14]-1.224744871391589*fr[21])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[7] = (0.7071067811865475*fr[16]-1.224744871391589*fr[22])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[8] = (0.7071067811865475*fr[18]-1.2247448713915887*fr[24])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[9] = (0.7071067811865475*fr[19]-1.2247448713915887*fr[27])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[10] = (0.7071067811865475*fr[26]-1.2247448713915887*fr[30])*jacob_cx_r_inv*jacob_vx_inv; 
  fUpwind_r[11] = (0.7071067811865475*fr[28]-1.2247448713915887*fr[31])*jacob_cx_r_inv*jacob_vx_inv; 

  fUpwind_l[0] = (-(1.8708286933869707*fc[17])+1.5811388300841895*fc[7]-1.224744871391589*fc[1]+0.7071067811865475*fc[0])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[1] = (-(1.8708286933869707*fc[23])+1.5811388300841898*fc[11]-1.224744871391589*fc[4]+0.7071067811865475*fc[2])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[2] = (-(1.8708286933869707*fc[25])+1.5811388300841898*fc[13]-1.224744871391589*fc[5]+0.7071067811865475*fc[3])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[3] = (-(1.8708286933869707*fc[29])+1.5811388300841895*fc[20]-1.224744871391589*fc[10]+0.7071067811865475*fc[6])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[4] = (0.7071067811865475*fc[8]-1.224744871391589*fc[12])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[5] = (0.7071067811865475*fc[9]-1.224744871391589*fc[15])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[6] = (0.7071067811865475*fc[14]-1.224744871391589*fc[21])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[7] = (0.7071067811865475*fc[16]-1.224744871391589*fc[22])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[8] = (0.7071067811865475*fc[18]-1.2247448713915887*fc[24])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[9] = (0.7071067811865475*fc[19]-1.2247448713915887*fc[27])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[10] = (0.7071067811865475*fc[26]-1.2247448713915887*fc[30])*jacob_cx_c_inv*jacob_vx_inv; 
  fUpwind_l[11] = (0.7071067811865475*fc[28]-1.2247448713915887*fc[31])*jacob_cx_c_inv*jacob_vx_inv; 

  } 
  Ghat_l[0] = hamil[3]*(4.183300132670378*fUpwind_l[5]+1.8708286933869707*fUpwind_l[0])+2.7386127875258306*fUpwind_l[2]*hamil[2]+1.224744871391589*fUpwind_l[0]*hamil[1]; 
  Ghat_l[1] = hamil[3]*(4.183300132670379*fUpwind_l[7]+1.8708286933869707*fUpwind_l[1])+2.7386127875258306*hamil[2]*fUpwind_l[3]+1.224744871391589*fUpwind_l[1]*hamil[1]; 
  Ghat_l[2] = 3.6742346141747664*hamil[3]*fUpwind_l[9]+2.4494897427831783*hamil[2]*fUpwind_l[5]+5.612486080160912*fUpwind_l[2]*hamil[3]+2.7386127875258306*fUpwind_l[0]*hamil[2]+1.224744871391589*hamil[1]*fUpwind_l[2]; 
  Ghat_l[3] = 3.674234614174766*hamil[3]*fUpwind_l[11]+2.4494897427831783*hamil[2]*fUpwind_l[7]+fUpwind_l[3]*(5.612486080160912*hamil[3]+1.224744871391589*hamil[1])+2.7386127875258306*fUpwind_l[1]*hamil[2]; 
  Ghat_l[4] = 2.7386127875258306*hamil[2]*fUpwind_l[6]+(1.8708286933869707*hamil[3]+1.224744871391589*hamil[1])*fUpwind_l[4]; 
  Ghat_l[5] = 2.405351177211819*hamil[2]*fUpwind_l[9]+(4.543441112511214*hamil[3]+1.224744871391589*hamil[1])*fUpwind_l[5]+4.183300132670378*fUpwind_l[0]*hamil[3]+2.4494897427831783*fUpwind_l[2]*hamil[2]; 
  Ghat_l[6] = (5.612486080160912*hamil[3]+1.224744871391589*hamil[1])*fUpwind_l[6]+2.7386127875258306*hamil[2]*fUpwind_l[4]; 
  Ghat_l[7] = 2.405351177211819*hamil[2]*fUpwind_l[11]+(4.543441112511214*hamil[3]+1.224744871391589*hamil[1])*fUpwind_l[7]+4.183300132670379*fUpwind_l[1]*hamil[3]+2.4494897427831783*hamil[2]*fUpwind_l[3]; 
  Ghat_l[8] = 2.7386127875258297*hamil[2]*fUpwind_l[10]+(1.8708286933869707*hamil[3]+1.224744871391589*hamil[1])*fUpwind_l[8]; 
  Ghat_l[9] = (4.365266951236265*hamil[3]+1.224744871391589*hamil[1])*fUpwind_l[9]+2.405351177211819*hamil[2]*fUpwind_l[5]+3.6742346141747664*fUpwind_l[2]*hamil[3]; 
  Ghat_l[10] = (5.612486080160912*hamil[3]+1.224744871391589*hamil[1])*fUpwind_l[10]+2.7386127875258297*hamil[2]*fUpwind_l[8]; 
  Ghat_l[11] = (4.365266951236265*hamil[3]+1.224744871391589*hamil[1])*fUpwind_l[11]+2.405351177211819*hamil[2]*fUpwind_l[7]+3.674234614174766*fUpwind_l[3]*hamil[3]; 

  Ghat_r[0] = hamil[3]*(4.183300132670378*fUpwind_r[5]+1.8708286933869707*fUpwind_r[0])+2.7386127875258306*fUpwind_r[2]*hamil[2]+1.224744871391589*fUpwind_r[0]*hamil[1]; 
  Ghat_r[1] = hamil[3]*(4.183300132670379*fUpwind_r[7]+1.8708286933869707*fUpwind_r[1])+2.7386127875258306*hamil[2]*fUpwind_r[3]+1.224744871391589*fUpwind_r[1]*hamil[1]; 
  Ghat_r[2] = 3.6742346141747664*hamil[3]*fUpwind_r[9]+2.4494897427831783*hamil[2]*fUpwind_r[5]+5.612486080160912*fUpwind_r[2]*hamil[3]+2.7386127875258306*fUpwind_r[0]*hamil[2]+1.224744871391589*hamil[1]*fUpwind_r[2]; 
  Ghat_r[3] = 3.674234614174766*hamil[3]*fUpwind_r[11]+2.4494897427831783*hamil[2]*fUpwind_r[7]+fUpwind_r[3]*(5.612486080160912*hamil[3]+1.224744871391589*hamil[1])+2.7386127875258306*fUpwind_r[1]*hamil[2]; 
  Ghat_r[4] = 2.7386127875258306*hamil[2]*fUpwind_r[6]+(1.8708286933869707*hamil[3]+1.224744871391589*hamil[1])*fUpwind_r[4]; 
  Ghat_r[5] = 2.405351177211819*hamil[2]*fUpwind_r[9]+(4.543441112511214*hamil[3]+1.224744871391589*hamil[1])*fUpwind_r[5]+4.183300132670378*fUpwind_r[0]*hamil[3]+2.4494897427831783*fUpwind_r[2]*hamil[2]; 
  Ghat_r[6] = (5.612486080160912*hamil[3]+1.224744871391589*hamil[1])*fUpwind_r[6]+2.7386127875258306*hamil[2]*fUpwind_r[4]; 
  Ghat_r[7] = 2.405351177211819*hamil[2]*fUpwind_r[11]+(4.543441112511214*hamil[3]+1.224744871391589*hamil[1])*fUpwind_r[7]+4.183300132670379*fUpwind_r[1]*hamil[3]+2.4494897427831783*hamil[2]*fUpwind_r[3]; 
  Ghat_r[8] = 2.7386127875258297*hamil[2]*fUpwind_r[10]+(1.8708286933869707*hamil[3]+1.224744871391589*hamil[1])*fUpwind_r[8]; 
  Ghat_r[9] = (4.365266951236265*hamil[3]+1.224744871391589*hamil[1])*fUpwind_r[9]+2.405351177211819*hamil[2]*fUpwind_r[5]+3.6742346141747664*fUpwind_r[2]*hamil[3]; 
  Ghat_r[10] = (5.612486080160912*hamil[3]+1.224744871391589*hamil[1])*fUpwind_r[10]+2.7386127875258297*hamil[2]*fUpwind_r[8]; 
  Ghat_r[11] = (4.365266951236265*hamil[3]+1.224744871391589*hamil[1])*fUpwind_r[11]+2.405351177211819*hamil[2]*fUpwind_r[7]+3.674234614174766*fUpwind_r[3]*hamil[3]; 

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
  out[17] += -(1.8708286933869707*(Ghat_r[0]+Ghat_l[0])*dv10*dx10); 
  out[18] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dv10*dx10; 
  out[19] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*dv10*dx10; 
  out[20] += (1.5811388300841895*Ghat_l[3]-1.5811388300841895*Ghat_r[3])*dv10*dx10; 
  out[21] += -(1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dv10*dx10); 
  out[22] += -(1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dv10*dx10); 
  out[23] += -(1.8708286933869707*(Ghat_r[1]+Ghat_l[1])*dv10*dx10); 
  out[24] += -(1.2247448713915887*(Ghat_r[8]+Ghat_l[8])*dv10*dx10); 
  out[25] += -(1.8708286933869707*(Ghat_r[2]+Ghat_l[2])*dv10*dx10); 
  out[26] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*dv10*dx10; 
  out[27] += -(1.2247448713915887*(Ghat_r[9]+Ghat_l[9])*dv10*dx10); 
  out[28] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*dv10*dx10; 
  out[29] += -(1.8708286933869707*(Ghat_r[3]+Ghat_l[3])*dv10*dx10); 
  out[30] += -(1.2247448713915887*(Ghat_r[10]+Ghat_l[10])*dv10*dx10); 
  out[31] += -(1.2247448713915887*(Ghat_r[11]+Ghat_l[11])*dv10*dx10); 

  return fabs(2.474873734152916*dv10*dx10*jacob_cx_c_inv*(2.6457513110645907*hamil[3]*jacob_vx_inv+1.7320508075688772*hamil[1]*jacob_vx_inv));

} 
