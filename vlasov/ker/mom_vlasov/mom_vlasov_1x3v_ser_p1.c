#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void mom_vlasov_M0_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
  const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += 2.8284271247461907*f[1]*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M2ij_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
  const double *vmap, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double *vmap_vx = &vmap[0]; 
  const double *vmap_vy = &vmap[4]; 
  const double *vmap_vz = &vmap[8]; 
  const double vmap_vx0_sq = vmap_vx[0]*vmap_vx[0]; 
  const double vmap_vx1_sq = vmap_vx[1]*vmap_vx[1]; 
  const double vmap_vy0_sq = vmap_vy[0]*vmap_vy[0]; 
  const double vmap_vy1_sq = vmap_vy[1]*vmap_vy[1]; 
  const double vmap_vz0_sq = vmap_vz[0]*vmap_vz[0]; 
  const double vmap_vz1_sq = vmap_vz[1]*vmap_vz[1]; 
  out[0] += (1.4142135623730951*f[0]*vmap_vx1_sq+1.4142135623730951*f[0]*vmap_vx0_sq+2.8284271247461907*vmap_vx[0]*vmap_vx[1]*f[2])*volFact; 
  out[1] += (1.4142135623730951*f[1]*vmap_vx1_sq+1.4142135623730951*f[1]*vmap_vx0_sq+2.8284271247461907*vmap_vx[0]*vmap_vx[1]*f[5])*volFact; 
  out[2] += (1.4142135623730951*vmap_vx[1]*vmap_vy[1]*f[7]+1.4142135623730951*vmap_vx[0]*vmap_vy[1]*f[3]+1.4142135623730951*vmap_vy[0]*vmap_vx[1]*f[2]+1.4142135623730951*f[0]*vmap_vx[0]*vmap_vy[0])*volFact; 
  out[3] += (1.4142135623730951*vmap_vx[1]*vmap_vy[1]*f[11]+1.4142135623730951*vmap_vx[0]*vmap_vy[1]*f[6]+1.4142135623730951*vmap_vy[0]*vmap_vx[1]*f[5]+1.4142135623730951*vmap_vx[0]*vmap_vy[0]*f[1])*volFact; 
  out[4] += (1.4142135623730951*vmap_vx[1]*vmap_vz[1]*f[9]+1.4142135623730951*vmap_vx[0]*vmap_vz[1]*f[4]+1.4142135623730951*vmap_vz[0]*vmap_vx[1]*f[2]+1.4142135623730951*f[0]*vmap_vx[0]*vmap_vz[0])*volFact; 
  out[5] += (1.4142135623730951*vmap_vx[1]*vmap_vz[1]*f[12]+1.4142135623730951*vmap_vx[0]*vmap_vz[1]*f[8]+1.4142135623730951*vmap_vz[0]*vmap_vx[1]*f[5]+1.4142135623730951*vmap_vx[0]*vmap_vz[0]*f[1])*volFact; 
  out[6] += (1.4142135623730951*f[0]*vmap_vy1_sq+1.4142135623730951*f[0]*vmap_vy0_sq+2.8284271247461907*vmap_vy[0]*vmap_vy[1]*f[3])*volFact; 
  out[7] += (1.4142135623730951*f[1]*vmap_vy1_sq+1.4142135623730951*f[1]*vmap_vy0_sq+2.8284271247461907*vmap_vy[0]*vmap_vy[1]*f[6])*volFact; 
  out[8] += (1.4142135623730951*vmap_vy[1]*vmap_vz[1]*f[10]+1.4142135623730951*vmap_vy[0]*vmap_vz[1]*f[4]+1.4142135623730951*vmap_vz[0]*vmap_vy[1]*f[3]+1.4142135623730951*f[0]*vmap_vy[0]*vmap_vz[0])*volFact; 
  out[9] += (1.4142135623730951*vmap_vy[1]*vmap_vz[1]*f[13]+1.4142135623730951*vmap_vy[0]*vmap_vz[1]*f[8]+1.4142135623730951*vmap_vz[0]*vmap_vy[1]*f[6]+1.4142135623730951*vmap_vy[0]*vmap_vz[0]*f[1])*volFact; 
  out[10] += (1.4142135623730951*f[0]*vmap_vz1_sq+1.4142135623730951*f[0]*vmap_vz0_sq+2.8284271247461907*vmap_vz[0]*vmap_vz[1]*f[4])*volFact; 
  out[11] += (1.4142135623730951*f[1]*vmap_vz1_sq+1.4142135623730951*f[1]*vmap_vz0_sq+2.8284271247461907*vmap_vz[0]*vmap_vz[1]*f[8])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M3ijk_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
  const double *vmap, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double *vmap_vx = &vmap[0]; 
  const double *vmap_vy = &vmap[4]; 
  const double *vmap_vz = &vmap[8]; 
  const double vmap_vx0_sq = vmap_vx[0]*vmap_vx[0]; 
  const double vmap_vx0_cu = vmap_vx[0]*vmap_vx[0]*vmap_vx[0]; 
  const double vmap_vx1_sq = vmap_vx[1]*vmap_vx[1]; 
  const double vmap_vx1_cu = vmap_vx[1]*vmap_vx[1]*vmap_vx[1]; 
  const double vmap_vy0_sq = vmap_vy[0]*vmap_vy[0]; 
  const double vmap_vy0_cu = vmap_vy[0]*vmap_vy[0]*vmap_vy[0]; 
  const double vmap_vy1_sq = vmap_vy[1]*vmap_vy[1]; 
  const double vmap_vy1_cu = vmap_vy[1]*vmap_vy[1]*vmap_vy[1]; 
  const double vmap_vz0_sq = vmap_vz[0]*vmap_vz[0]; 
  const double vmap_vz0_cu = vmap_vz[0]*vmap_vz[0]*vmap_vz[0]; 
  const double vmap_vz1_sq = vmap_vz[1]*vmap_vz[1]; 
  const double vmap_vz1_cu = vmap_vz[1]*vmap_vz[1]*vmap_vz[1]; 
  out[0] += (3.0*f[0]*vmap_vx[0]*vmap_vx1_sq+1.8*f[2]*vmap_vx1_cu+3.0*vmap_vx[1]*f[2]*vmap_vx0_sq+f[0]*vmap_vx0_cu)*volFact; 
  out[1] += (3.0*vmap_vx[0]*f[1]*vmap_vx1_sq+1.8*f[5]*vmap_vx1_cu+3.0*vmap_vx[1]*f[5]*vmap_vx0_sq+f[1]*vmap_vx0_cu)*volFact; 
  out[2] += (vmap_vy[1]*f[3]*vmap_vx1_sq+f[0]*vmap_vy[0]*vmap_vx1_sq+vmap_vy[1]*f[3]*vmap_vx0_sq+f[0]*vmap_vy[0]*vmap_vx0_sq+2.0*vmap_vx[0]*vmap_vx[1]*vmap_vy[1]*f[7]+2.0*vmap_vx[0]*vmap_vy[0]*vmap_vx[1]*f[2])*volFact; 
  out[3] += (vmap_vy[1]*f[6]*vmap_vx1_sq+vmap_vy[0]*f[1]*vmap_vx1_sq+vmap_vy[1]*f[6]*vmap_vx0_sq+vmap_vy[0]*f[1]*vmap_vx0_sq+2.0*vmap_vx[0]*vmap_vx[1]*vmap_vy[1]*f[11]+2.0*vmap_vx[0]*vmap_vy[0]*vmap_vx[1]*f[5])*volFact; 
  out[4] += (vmap_vz[1]*f[4]*vmap_vx1_sq+f[0]*vmap_vz[0]*vmap_vx1_sq+vmap_vz[1]*f[4]*vmap_vx0_sq+f[0]*vmap_vz[0]*vmap_vx0_sq+2.0*vmap_vx[0]*vmap_vx[1]*vmap_vz[1]*f[9]+2.0*vmap_vx[0]*vmap_vz[0]*vmap_vx[1]*f[2])*volFact; 
  out[5] += (vmap_vz[1]*f[8]*vmap_vx1_sq+vmap_vz[0]*f[1]*vmap_vx1_sq+vmap_vz[1]*f[8]*vmap_vx0_sq+vmap_vz[0]*f[1]*vmap_vx0_sq+2.0*vmap_vx[0]*vmap_vx[1]*vmap_vz[1]*f[12]+2.0*vmap_vx[0]*vmap_vz[0]*vmap_vx[1]*f[5])*volFact; 
  out[6] += (vmap_vx[1]*f[2]*vmap_vy1_sq+f[0]*vmap_vx[0]*vmap_vy1_sq+vmap_vx[1]*f[2]*vmap_vy0_sq+f[0]*vmap_vx[0]*vmap_vy0_sq+2.0*vmap_vy[0]*vmap_vx[1]*vmap_vy[1]*f[7]+2.0*vmap_vx[0]*vmap_vy[0]*vmap_vy[1]*f[3])*volFact; 
  out[7] += (vmap_vx[1]*f[5]*vmap_vy1_sq+vmap_vx[0]*f[1]*vmap_vy1_sq+vmap_vx[1]*f[5]*vmap_vy0_sq+vmap_vx[0]*f[1]*vmap_vy0_sq+2.0*vmap_vy[0]*vmap_vx[1]*vmap_vy[1]*f[11]+2.0*vmap_vx[0]*vmap_vy[0]*vmap_vy[1]*f[6])*volFact; 
  out[8] += (vmap_vx[1]*vmap_vy[1]*vmap_vz[1]*f[14]+vmap_vx[0]*vmap_vy[1]*vmap_vz[1]*f[10]+vmap_vy[0]*vmap_vx[1]*vmap_vz[1]*f[9]+vmap_vz[0]*vmap_vx[1]*vmap_vy[1]*f[7]+vmap_vx[0]*vmap_vy[0]*vmap_vz[1]*f[4]+vmap_vx[0]*vmap_vz[0]*vmap_vy[1]*f[3]+vmap_vy[0]*vmap_vz[0]*vmap_vx[1]*f[2]+f[0]*vmap_vx[0]*vmap_vy[0]*vmap_vz[0])*volFact; 
  out[9] += (vmap_vx[1]*vmap_vy[1]*vmap_vz[1]*f[15]+vmap_vx[0]*vmap_vy[1]*vmap_vz[1]*f[13]+vmap_vy[0]*vmap_vx[1]*vmap_vz[1]*f[12]+vmap_vz[0]*vmap_vx[1]*vmap_vy[1]*f[11]+vmap_vx[0]*vmap_vy[0]*vmap_vz[1]*f[8]+vmap_vx[0]*vmap_vz[0]*vmap_vy[1]*f[6]+vmap_vy[0]*vmap_vz[0]*vmap_vx[1]*f[5]+vmap_vx[0]*vmap_vy[0]*vmap_vz[0]*f[1])*volFact; 
  out[10] += (vmap_vx[1]*f[2]*vmap_vz1_sq+f[0]*vmap_vx[0]*vmap_vz1_sq+vmap_vx[1]*f[2]*vmap_vz0_sq+f[0]*vmap_vx[0]*vmap_vz0_sq+2.0*vmap_vz[0]*vmap_vx[1]*vmap_vz[1]*f[9]+2.0*vmap_vx[0]*vmap_vz[0]*vmap_vz[1]*f[4])*volFact; 
  out[11] += (vmap_vx[1]*f[5]*vmap_vz1_sq+vmap_vx[0]*f[1]*vmap_vz1_sq+vmap_vx[1]*f[5]*vmap_vz0_sq+vmap_vx[0]*f[1]*vmap_vz0_sq+2.0*vmap_vz[0]*vmap_vx[1]*vmap_vz[1]*f[12]+2.0*vmap_vx[0]*vmap_vz[0]*vmap_vz[1]*f[8])*volFact; 
  out[12] += (3.0*f[0]*vmap_vy[0]*vmap_vy1_sq+1.8*f[3]*vmap_vy1_cu+3.0*vmap_vy[1]*f[3]*vmap_vy0_sq+f[0]*vmap_vy0_cu)*volFact; 
  out[13] += (3.0*vmap_vy[0]*f[1]*vmap_vy1_sq+1.8*f[6]*vmap_vy1_cu+3.0*vmap_vy[1]*f[6]*vmap_vy0_sq+f[1]*vmap_vy0_cu)*volFact; 
  out[14] += (vmap_vz[1]*f[4]*vmap_vy1_sq+f[0]*vmap_vz[0]*vmap_vy1_sq+vmap_vz[1]*f[4]*vmap_vy0_sq+f[0]*vmap_vz[0]*vmap_vy0_sq+2.0*vmap_vy[0]*vmap_vy[1]*vmap_vz[1]*f[10]+2.0*vmap_vy[0]*vmap_vz[0]*vmap_vy[1]*f[3])*volFact; 
  out[15] += (vmap_vz[1]*f[8]*vmap_vy1_sq+vmap_vz[0]*f[1]*vmap_vy1_sq+vmap_vz[1]*f[8]*vmap_vy0_sq+vmap_vz[0]*f[1]*vmap_vy0_sq+2.0*vmap_vy[0]*vmap_vy[1]*vmap_vz[1]*f[13]+2.0*vmap_vy[0]*vmap_vz[0]*vmap_vy[1]*f[6])*volFact; 
  out[16] += (vmap_vy[1]*f[3]*vmap_vz1_sq+f[0]*vmap_vy[0]*vmap_vz1_sq+vmap_vy[1]*f[3]*vmap_vz0_sq+f[0]*vmap_vy[0]*vmap_vz0_sq+2.0*vmap_vz[0]*vmap_vy[1]*vmap_vz[1]*f[10]+2.0*vmap_vy[0]*vmap_vz[0]*vmap_vz[1]*f[4])*volFact; 
  out[17] += (vmap_vy[1]*f[6]*vmap_vz1_sq+vmap_vy[0]*f[1]*vmap_vz1_sq+vmap_vy[1]*f[6]*vmap_vz0_sq+vmap_vy[0]*f[1]*vmap_vz0_sq+2.0*vmap_vz[0]*vmap_vy[1]*vmap_vz[1]*f[13]+2.0*vmap_vy[0]*vmap_vz[0]*vmap_vz[1]*f[8])*volFact; 
  out[18] += (3.0*f[0]*vmap_vz[0]*vmap_vz1_sq+1.8*f[4]*vmap_vz1_cu+3.0*vmap_vz[1]*f[4]*vmap_vz0_sq+f[0]*vmap_vz0_cu)*volFact; 
  out[19] += (3.0*vmap_vz[0]*f[1]*vmap_vz1_sq+1.8*f[8]*vmap_vz1_cu+3.0*vmap_vz[1]*f[8]*vmap_vz0_sq+f[1]*vmap_vz0_cu)*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_M1i_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  double dv12 = 2.0/dxv[3]; 
  const double *jacob_vz = &jacob_vel[4]; 
  const double jacob_vz_inv = 1.0/jacob_vz[0]; 
  out[0] += (1.7320508075688772*hamil[7]*f[10]+1.7320508075688772*f[4]*hamil[5]+1.7320508075688772*f[3]*hamil[4]+1.7320508075688772*f[0]*hamil[1])*dv10*jacob_vx_inv*volFact; 
  out[1] += (1.7320508075688772*hamil[7]*f[13]+1.7320508075688772*hamil[5]*f[8]+1.7320508075688772*hamil[4]*f[6]+1.7320508075688772*f[1]*hamil[1])*dv10*jacob_vx_inv*volFact; 
  out[2] += (1.7320508075688772*hamil[7]*f[9]+1.7320508075688772*f[4]*hamil[6]+1.7320508075688772*f[2]*hamil[4]+1.7320508075688772*f[0]*hamil[2])*dv11*jacob_vy_inv*volFact; 
  out[3] += (1.7320508075688772*hamil[7]*f[12]+1.7320508075688772*hamil[6]*f[8]+1.7320508075688772*hamil[4]*f[5]+1.7320508075688772*f[1]*hamil[2])*dv11*jacob_vy_inv*volFact; 
  out[4] += (1.7320508075688772*f[7]*hamil[7]+1.7320508075688772*f[3]*hamil[6]+1.7320508075688772*f[2]*hamil[5]+1.7320508075688772*f[0]*hamil[3])*dv12*jacob_vz_inv*volFact; 
  out[5] += (1.7320508075688772*hamil[7]*f[11]+1.7320508075688772*f[6]*hamil[6]+1.7320508075688772*f[5]*hamil[5]+1.7320508075688772*f[1]*hamil[3])*dv12*jacob_vz_inv*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_M2_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  out[0] += (hamil[7]*f[14]+hamil[6]*f[10]+hamil[5]*f[9]+hamil[4]*f[7]+hamil[3]*f[4]+hamil[2]*f[3]+hamil[1]*f[2]+f[0]*hamil[0])*volFact; 
  out[1] += (hamil[7]*f[15]+hamil[6]*f[13]+hamil[5]*f[12]+hamil[4]*f[11]+hamil[3]*f[8]+hamil[2]*f[6]+hamil[1]*f[5]+hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_M3i_1x3v_ser_p1(const double *w, const double *dxv, const int *idx,
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  double dv12 = 2.0/dxv[3]; 
  const double *jacob_vz = &jacob_vel[4]; 
  const double jacob_vz_inv = 1.0/jacob_vz[0]; 
  double dH_dpx[8] = {0.0}; 
  dH_dpx[0] = 1.7320508075688772*hamil[1]*dv10*jacob_vx_inv; 
  dH_dpx[2] = 1.7320508075688772*hamil[4]*dv10*jacob_vx_inv; 
  dH_dpx[3] = 1.7320508075688772*hamil[5]*dv10*jacob_vx_inv; 
  dH_dpx[6] = 1.7320508075688772*hamil[7]*dv10*jacob_vx_inv; 

  double dH_dpy[8] = {0.0}; 
  dH_dpy[0] = 1.7320508075688772*hamil[2]*dv11*jacob_vy_inv; 
  dH_dpy[1] = 1.7320508075688772*hamil[4]*dv11*jacob_vy_inv; 
  dH_dpy[3] = 1.7320508075688772*hamil[6]*dv11*jacob_vy_inv; 
  dH_dpy[5] = 1.7320508075688772*hamil[7]*dv11*jacob_vy_inv; 

  double dH_dpz[8] = {0.0}; 
  dH_dpz[0] = 1.7320508075688772*hamil[3]*dv12*jacob_vz_inv; 
  dH_dpz[1] = 1.7320508075688772*hamil[5]*dv12*jacob_vz_inv; 
  dH_dpz[2] = 1.7320508075688772*hamil[6]*dv12*jacob_vz_inv; 
  dH_dpz[4] = 1.7320508075688772*hamil[7]*dv12*jacob_vz_inv; 

  out[0] += (0.3535533905932737*dH_dpx[0]*hamil[7]*f[14]+0.3535533905932737*hamil[1]*dH_dpx[6]*f[14]+0.3535533905932737*dH_dpx[2]*hamil[5]*f[14]+0.3535533905932737*dH_dpx[3]*hamil[4]*f[14]+0.3535533905932737*dH_dpx[0]*hamil[6]*f[10]+0.3535533905932737*hamil[0]*dH_dpx[6]*f[10]+0.3535533905932737*dH_dpx[2]*hamil[3]*f[10]+0.3535533905932737*hamil[2]*dH_dpx[3]*f[10]+0.3535533905932737*dH_dpx[2]*hamil[7]*f[9]+0.3535533905932737*hamil[4]*dH_dpx[6]*f[9]+0.3535533905932737*dH_dpx[0]*hamil[5]*f[9]+0.3535533905932737*hamil[1]*dH_dpx[3]*f[9]+0.3535533905932737*dH_dpx[3]*f[7]*hamil[7]+0.3535533905932737*f[2]*dH_dpx[6]*hamil[7]+0.3535533905932737*hamil[5]*dH_dpx[6]*f[7]+0.3535533905932737*dH_dpx[0]*hamil[4]*f[7]+0.3535533905932737*hamil[1]*dH_dpx[2]*f[7]+0.3535533905932737*f[0]*dH_dpx[6]*hamil[6]+0.3535533905932737*dH_dpx[2]*f[4]*hamil[6]+0.3535533905932737*dH_dpx[3]*f[3]*hamil[6]+0.3535533905932737*hamil[2]*f[4]*dH_dpx[6]+0.3535533905932737*f[3]*hamil[3]*dH_dpx[6]+0.3535533905932737*f[2]*dH_dpx[3]*hamil[5]+0.3535533905932737*dH_dpx[2]*f[2]*hamil[4]+0.3535533905932737*dH_dpx[0]*hamil[3]*f[4]+0.3535533905932737*hamil[0]*dH_dpx[3]*f[4]+0.3535533905932737*f[0]*dH_dpx[3]*hamil[3]+0.3535533905932737*dH_dpx[0]*hamil[2]*f[3]+0.3535533905932737*hamil[0]*dH_dpx[2]*f[3]+0.3535533905932737*f[0]*dH_dpx[2]*hamil[2]+0.3535533905932737*dH_dpx[0]*hamil[1]*f[2]+0.3535533905932737*dH_dpx[0]*f[0]*hamil[0])*volFact; 
  out[1] += (0.3535533905932737*dH_dpx[0]*hamil[7]*f[15]+0.3535533905932737*hamil[1]*dH_dpx[6]*f[15]+0.3535533905932737*dH_dpx[2]*hamil[5]*f[15]+0.3535533905932737*dH_dpx[3]*hamil[4]*f[15]+0.3535533905932737*dH_dpx[0]*hamil[6]*f[13]+0.3535533905932737*hamil[0]*dH_dpx[6]*f[13]+0.3535533905932737*dH_dpx[2]*hamil[3]*f[13]+0.3535533905932737*hamil[2]*dH_dpx[3]*f[13]+0.3535533905932737*dH_dpx[2]*hamil[7]*f[12]+0.3535533905932737*hamil[4]*dH_dpx[6]*f[12]+0.3535533905932737*dH_dpx[0]*hamil[5]*f[12]+0.3535533905932737*hamil[1]*dH_dpx[3]*f[12]+0.3535533905932737*dH_dpx[3]*hamil[7]*f[11]+0.3535533905932737*hamil[5]*dH_dpx[6]*f[11]+0.3535533905932737*dH_dpx[0]*hamil[4]*f[11]+0.3535533905932737*hamil[1]*dH_dpx[2]*f[11]+0.3535533905932737*dH_dpx[2]*hamil[6]*f[8]+0.3535533905932737*hamil[2]*dH_dpx[6]*f[8]+0.3535533905932737*dH_dpx[0]*hamil[3]*f[8]+0.3535533905932737*hamil[0]*dH_dpx[3]*f[8]+0.3535533905932737*f[5]*dH_dpx[6]*hamil[7]+0.3535533905932737*dH_dpx[3]*f[6]*hamil[6]+0.3535533905932737*f[1]*dH_dpx[6]*hamil[6]+0.3535533905932737*hamil[3]*dH_dpx[6]*f[6]+0.3535533905932737*dH_dpx[0]*hamil[2]*f[6]+0.3535533905932737*hamil[0]*dH_dpx[2]*f[6]+0.3535533905932737*dH_dpx[3]*f[5]*hamil[5]+0.3535533905932737*dH_dpx[2]*hamil[4]*f[5]+0.3535533905932737*dH_dpx[0]*hamil[1]*f[5]+0.3535533905932737*f[1]*dH_dpx[3]*hamil[3]+0.3535533905932737*f[1]*dH_dpx[2]*hamil[2]+0.3535533905932737*dH_dpx[0]*hamil[0]*f[1])*volFact; 
  out[2] += (0.3535533905932737*dH_dpy[0]*hamil[7]*f[14]+0.3535533905932737*dH_dpy[1]*hamil[6]*f[14]+0.3535533905932737*hamil[2]*dH_dpy[5]*f[14]+0.3535533905932737*dH_dpy[3]*hamil[4]*f[14]+0.3535533905932737*dH_dpy[1]*hamil[7]*f[10]+0.3535533905932737*dH_dpy[0]*hamil[6]*f[10]+0.3535533905932737*hamil[4]*dH_dpy[5]*f[10]+0.3535533905932737*hamil[2]*dH_dpy[3]*f[10]+0.3535533905932737*dH_dpy[0]*hamil[5]*f[9]+0.3535533905932737*hamil[0]*dH_dpy[5]*f[9]+0.3535533905932737*dH_dpy[1]*hamil[3]*f[9]+0.3535533905932737*hamil[1]*dH_dpy[3]*f[9]+0.3535533905932737*dH_dpy[3]*f[7]*hamil[7]+0.3535533905932737*f[3]*dH_dpy[5]*hamil[7]+0.3535533905932737*dH_dpy[5]*hamil[6]*f[7]+0.3535533905932737*dH_dpy[0]*hamil[4]*f[7]+0.3535533905932737*dH_dpy[1]*hamil[2]*f[7]+0.3535533905932737*dH_dpy[3]*f[3]*hamil[6]+0.3535533905932737*f[0]*dH_dpy[5]*hamil[5]+0.3535533905932737*dH_dpy[1]*f[4]*hamil[5]+0.3535533905932737*f[2]*dH_dpy[3]*hamil[5]+0.3535533905932737*hamil[1]*f[4]*dH_dpy[5]+0.3535533905932737*f[2]*hamil[3]*dH_dpy[5]+0.3535533905932737*dH_dpy[1]*f[3]*hamil[4]+0.3535533905932737*dH_dpy[0]*hamil[3]*f[4]+0.3535533905932737*hamil[0]*dH_dpy[3]*f[4]+0.3535533905932737*f[0]*dH_dpy[3]*hamil[3]+0.3535533905932737*dH_dpy[0]*hamil[2]*f[3]+0.3535533905932737*dH_dpy[0]*hamil[1]*f[2]+0.3535533905932737*hamil[0]*dH_dpy[1]*f[2]+0.3535533905932737*f[0]*dH_dpy[1]*hamil[1]+0.3535533905932737*dH_dpy[0]*f[0]*hamil[0])*volFact; 
  out[3] += (0.3535533905932737*dH_dpy[0]*hamil[7]*f[15]+0.3535533905932737*dH_dpy[1]*hamil[6]*f[15]+0.3535533905932737*hamil[2]*dH_dpy[5]*f[15]+0.3535533905932737*dH_dpy[3]*hamil[4]*f[15]+0.3535533905932737*dH_dpy[1]*hamil[7]*f[13]+0.3535533905932737*dH_dpy[0]*hamil[6]*f[13]+0.3535533905932737*hamil[4]*dH_dpy[5]*f[13]+0.3535533905932737*hamil[2]*dH_dpy[3]*f[13]+0.3535533905932737*dH_dpy[0]*hamil[5]*f[12]+0.3535533905932737*hamil[0]*dH_dpy[5]*f[12]+0.3535533905932737*dH_dpy[1]*hamil[3]*f[12]+0.3535533905932737*hamil[1]*dH_dpy[3]*f[12]+0.3535533905932737*dH_dpy[3]*hamil[7]*f[11]+0.3535533905932737*dH_dpy[5]*hamil[6]*f[11]+0.3535533905932737*dH_dpy[0]*hamil[4]*f[11]+0.3535533905932737*dH_dpy[1]*hamil[2]*f[11]+0.3535533905932737*dH_dpy[1]*hamil[5]*f[8]+0.3535533905932737*hamil[1]*dH_dpy[5]*f[8]+0.3535533905932737*dH_dpy[0]*hamil[3]*f[8]+0.3535533905932737*hamil[0]*dH_dpy[3]*f[8]+0.3535533905932737*dH_dpy[5]*f[6]*hamil[7]+0.3535533905932737*dH_dpy[3]*f[6]*hamil[6]+0.3535533905932737*dH_dpy[1]*hamil[4]*f[6]+0.3535533905932737*dH_dpy[0]*hamil[2]*f[6]+0.3535533905932737*dH_dpy[3]*f[5]*hamil[5]+0.3535533905932737*f[1]*dH_dpy[5]*hamil[5]+0.3535533905932737*hamil[3]*dH_dpy[5]*f[5]+0.3535533905932737*dH_dpy[0]*hamil[1]*f[5]+0.3535533905932737*hamil[0]*dH_dpy[1]*f[5]+0.3535533905932737*f[1]*dH_dpy[3]*hamil[3]+0.3535533905932737*dH_dpy[1]*f[1]*hamil[1]+0.3535533905932737*dH_dpy[0]*hamil[0]*f[1])*volFact; 
  out[4] += (0.3535533905932737*dH_dpz[0]*hamil[7]*f[14]+0.3535533905932737*dH_dpz[1]*hamil[6]*f[14]+0.3535533905932737*dH_dpz[2]*hamil[5]*f[14]+0.3535533905932737*hamil[3]*dH_dpz[4]*f[14]+0.3535533905932737*dH_dpz[1]*hamil[7]*f[10]+0.3535533905932737*dH_dpz[0]*hamil[6]*f[10]+0.3535533905932737*dH_dpz[4]*hamil[5]*f[10]+0.3535533905932737*dH_dpz[2]*hamil[3]*f[10]+0.3535533905932737*dH_dpz[2]*hamil[7]*f[9]+0.3535533905932737*dH_dpz[4]*hamil[6]*f[9]+0.3535533905932737*dH_dpz[0]*hamil[5]*f[9]+0.3535533905932737*dH_dpz[1]*hamil[3]*f[9]+0.3535533905932737*dH_dpz[4]*f[4]*hamil[7]+0.3535533905932737*dH_dpz[0]*hamil[4]*f[7]+0.3535533905932737*hamil[0]*dH_dpz[4]*f[7]+0.3535533905932737*dH_dpz[1]*hamil[2]*f[7]+0.3535533905932737*hamil[1]*dH_dpz[2]*f[7]+0.3535533905932737*dH_dpz[2]*f[4]*hamil[6]+0.3535533905932737*dH_dpz[1]*f[4]*hamil[5]+0.3535533905932737*f[0]*dH_dpz[4]*hamil[4]+0.3535533905932737*dH_dpz[1]*f[3]*hamil[4]+0.3535533905932737*dH_dpz[2]*f[2]*hamil[4]+0.3535533905932737*dH_dpz[0]*hamil[3]*f[4]+0.3535533905932737*hamil[1]*f[3]*dH_dpz[4]+0.3535533905932737*f[2]*hamil[2]*dH_dpz[4]+0.3535533905932737*dH_dpz[0]*hamil[2]*f[3]+0.3535533905932737*hamil[0]*dH_dpz[2]*f[3]+0.3535533905932737*f[0]*dH_dpz[2]*hamil[2]+0.3535533905932737*dH_dpz[0]*hamil[1]*f[2]+0.3535533905932737*hamil[0]*dH_dpz[1]*f[2]+0.3535533905932737*f[0]*dH_dpz[1]*hamil[1]+0.3535533905932737*dH_dpz[0]*f[0]*hamil[0])*volFact; 
  out[5] += (0.3535533905932737*dH_dpz[0]*hamil[7]*f[15]+0.3535533905932737*dH_dpz[1]*hamil[6]*f[15]+0.3535533905932737*dH_dpz[2]*hamil[5]*f[15]+0.3535533905932737*hamil[3]*dH_dpz[4]*f[15]+0.3535533905932737*dH_dpz[1]*hamil[7]*f[13]+0.3535533905932737*dH_dpz[0]*hamil[6]*f[13]+0.3535533905932737*dH_dpz[4]*hamil[5]*f[13]+0.3535533905932737*dH_dpz[2]*hamil[3]*f[13]+0.3535533905932737*dH_dpz[2]*hamil[7]*f[12]+0.3535533905932737*dH_dpz[4]*hamil[6]*f[12]+0.3535533905932737*dH_dpz[0]*hamil[5]*f[12]+0.3535533905932737*dH_dpz[1]*hamil[3]*f[12]+0.3535533905932737*dH_dpz[0]*hamil[4]*f[11]+0.3535533905932737*hamil[0]*dH_dpz[4]*f[11]+0.3535533905932737*dH_dpz[1]*hamil[2]*f[11]+0.3535533905932737*hamil[1]*dH_dpz[2]*f[11]+0.3535533905932737*dH_dpz[4]*hamil[7]*f[8]+0.3535533905932737*dH_dpz[2]*hamil[6]*f[8]+0.3535533905932737*dH_dpz[1]*hamil[5]*f[8]+0.3535533905932737*dH_dpz[0]*hamil[3]*f[8]+0.3535533905932737*dH_dpz[1]*hamil[4]*f[6]+0.3535533905932737*hamil[1]*dH_dpz[4]*f[6]+0.3535533905932737*dH_dpz[0]*hamil[2]*f[6]+0.3535533905932737*hamil[0]*dH_dpz[2]*f[6]+0.3535533905932737*dH_dpz[2]*hamil[4]*f[5]+0.3535533905932737*hamil[2]*dH_dpz[4]*f[5]+0.3535533905932737*dH_dpz[0]*hamil[1]*f[5]+0.3535533905932737*hamil[0]*dH_dpz[1]*f[5]+0.3535533905932737*f[1]*dH_dpz[4]*hamil[4]+0.3535533905932737*f[1]*dH_dpz[2]*hamil[2]+0.3535533905932737*dH_dpz[1]*f[1]*hamil[1]+0.3535533905932737*dH_dpz[0]*hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_five_moments_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  double dv12 = 2.0/dxv[3]; 
  const double *jacob_vz = &jacob_vel[4]; 
  const double jacob_vz_inv = 1.0/jacob_vz[0]; 
  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += 2.8284271247461907*f[1]*volFact; 
  out[2] += (1.7320508075688772*hamil[7]*f[10]+1.7320508075688772*f[4]*hamil[5]+1.7320508075688772*f[3]*hamil[4]+1.7320508075688772*f[0]*hamil[1])*dv10*jacob_vx_inv*volFact; 
  out[3] += (1.7320508075688772*hamil[7]*f[13]+1.7320508075688772*hamil[5]*f[8]+1.7320508075688772*hamil[4]*f[6]+1.7320508075688772*f[1]*hamil[1])*dv10*jacob_vx_inv*volFact; 
  out[4] += (1.7320508075688772*hamil[7]*f[9]+1.7320508075688772*f[4]*hamil[6]+1.7320508075688772*f[2]*hamil[4]+1.7320508075688772*f[0]*hamil[2])*dv11*jacob_vy_inv*volFact; 
  out[5] += (1.7320508075688772*hamil[7]*f[12]+1.7320508075688772*hamil[6]*f[8]+1.7320508075688772*hamil[4]*f[5]+1.7320508075688772*f[1]*hamil[2])*dv11*jacob_vy_inv*volFact; 
  out[6] += (1.7320508075688772*f[7]*hamil[7]+1.7320508075688772*f[3]*hamil[6]+1.7320508075688772*f[2]*hamil[5]+1.7320508075688772*f[0]*hamil[3])*dv12*jacob_vz_inv*volFact; 
  out[7] += (1.7320508075688772*hamil[7]*f[11]+1.7320508075688772*f[6]*hamil[6]+1.7320508075688772*f[5]*hamil[5]+1.7320508075688772*f[1]*hamil[3])*dv12*jacob_vz_inv*volFact; 
  out[8] += (hamil[7]*f[14]+hamil[6]*f[10]+hamil[5]*f[9]+hamil[4]*f[7]+hamil[3]*f[4]+hamil[2]*f[3]+hamil[1]*f[2]+f[0]*hamil[0])*volFact; 
  out[9] += (hamil[7]*f[15]+hamil[6]*f[13]+hamil[5]*f[12]+hamil[4]*f[11]+hamil[3]*f[8]+hamil[2]*f[6]+hamil[1]*f[5]+hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_int_five_moments_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]*0.0625; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  double dv12 = 2.0/dxv[3]; 
  const double *jacob_vz = &jacob_vel[4]; 
  const double jacob_vz_inv = 1.0/jacob_vz[0]; 
  out[0] += 4.0*f[0]*volFact; 
  out[1] += (2.4494897427831783*hamil[7]*f[10]+2.4494897427831783*f[4]*hamil[5]+2.4494897427831783*f[3]*hamil[4]+2.4494897427831783*f[0]*hamil[1])*dv10*jacob_vx_inv*volFact; 
  out[2] += (2.4494897427831783*hamil[7]*f[9]+2.4494897427831783*f[4]*hamil[6]+2.4494897427831783*f[2]*hamil[4]+2.4494897427831783*f[0]*hamil[2])*dv11*jacob_vy_inv*volFact; 
  out[3] += (2.4494897427831783*f[7]*hamil[7]+2.4494897427831783*f[3]*hamil[6]+2.4494897427831783*f[2]*hamil[5]+2.4494897427831783*f[0]*hamil[3])*dv12*jacob_vz_inv*volFact; 
  out[4] += (1.4142135623730951*hamil[7]*f[14]+1.4142135623730951*hamil[6]*f[10]+1.4142135623730951*hamil[5]*f[9]+1.4142135623730951*hamil[4]*f[7]+1.4142135623730951*hamil[3]*f[4]+1.4142135623730951*hamil[2]*f[3]+1.4142135623730951*hamil[1]*f[2]+1.4142135623730951*f[0]*hamil[0])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_sparse_M1i_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  double dv12 = 2.0/dxv[3]; 
  const double *jacob_vz = &jacob_vel[4]; 
  const double jacob_vz_inv = 1.0/jacob_vz[0]; 
  out[0] += 1.7320508075688772*f[0]*hamil[1]*dv10*jacob_vx_inv*volFact; 
  out[1] += 1.7320508075688772*f[1]*hamil[1]*dv10*jacob_vx_inv*volFact; 
  out[2] += 1.7320508075688772*f[0]*hamil[2]*dv11*jacob_vy_inv*volFact; 
  out[3] += 1.7320508075688772*f[1]*hamil[2]*dv11*jacob_vy_inv*volFact; 
  out[4] += 1.7320508075688772*f[0]*hamil[3]*dv12*jacob_vz_inv*volFact; 
  out[5] += 1.7320508075688772*f[1]*hamil[3]*dv12*jacob_vz_inv*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_sparse_M2_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  out[0] += (hamil[3]*f[4]+hamil[2]*f[3]+hamil[1]*f[2]+f[0]*hamil[0])*volFact; 
  out[1] += (hamil[3]*f[8]+hamil[2]*f[6]+hamil[1]*f[5]+hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_sparse_M3i_1x3v_ser_p1(const double *w, const double *dxv, const int *idx,
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  double dv12 = 2.0/dxv[3]; 
  const double *jacob_vz = &jacob_vel[4]; 
  const double jacob_vz_inv = 1.0/jacob_vz[0]; 
  double dH_dpx[8] = {0.0}; 
  dH_dpx[0] = 1.7320508075688772*hamil[1]*dv10*jacob_vx_inv; 

  double dH_dpy[8] = {0.0}; 
  dH_dpy[0] = 1.7320508075688772*hamil[2]*dv11*jacob_vy_inv; 

  double dH_dpz[8] = {0.0}; 
  dH_dpz[0] = 1.7320508075688772*hamil[3]*dv12*jacob_vz_inv; 

  out[0] += (0.3535533905932737*dH_dpx[0]*hamil[3]*f[4]+0.3535533905932737*dH_dpx[0]*hamil[2]*f[3]+0.3535533905932737*dH_dpx[0]*hamil[1]*f[2]+0.3535533905932737*dH_dpx[0]*f[0]*hamil[0])*volFact; 
  out[1] += (0.3535533905932737*dH_dpx[0]*hamil[3]*f[8]+0.3535533905932737*dH_dpx[0]*hamil[2]*f[6]+0.3535533905932737*dH_dpx[0]*hamil[1]*f[5]+0.3535533905932737*dH_dpx[0]*hamil[0]*f[1])*volFact; 
  out[2] += (0.3535533905932737*dH_dpy[0]*hamil[3]*f[4]+0.3535533905932737*dH_dpy[0]*hamil[2]*f[3]+0.3535533905932737*dH_dpy[0]*hamil[1]*f[2]+0.3535533905932737*dH_dpy[0]*f[0]*hamil[0])*volFact; 
  out[3] += (0.3535533905932737*dH_dpy[0]*hamil[3]*f[8]+0.3535533905932737*dH_dpy[0]*hamil[2]*f[6]+0.3535533905932737*dH_dpy[0]*hamil[1]*f[5]+0.3535533905932737*dH_dpy[0]*hamil[0]*f[1])*volFact; 
  out[4] += (0.3535533905932737*dH_dpz[0]*hamil[3]*f[4]+0.3535533905932737*dH_dpz[0]*hamil[2]*f[3]+0.3535533905932737*dH_dpz[0]*hamil[1]*f[2]+0.3535533905932737*dH_dpz[0]*f[0]*hamil[0])*volFact; 
  out[5] += (0.3535533905932737*dH_dpz[0]*hamil[3]*f[8]+0.3535533905932737*dH_dpz[0]*hamil[2]*f[6]+0.3535533905932737*dH_dpz[0]*hamil[1]*f[5]+0.3535533905932737*dH_dpz[0]*hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_sparse_five_moments_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  double dv12 = 2.0/dxv[3]; 
  const double *jacob_vz = &jacob_vel[4]; 
  const double jacob_vz_inv = 1.0/jacob_vz[0]; 
  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += 2.8284271247461907*f[1]*volFact; 
  out[2] += 1.7320508075688772*f[0]*hamil[1]*dv10*jacob_vx_inv*volFact; 
  out[3] += 1.7320508075688772*f[1]*hamil[1]*dv10*jacob_vx_inv*volFact; 
  out[4] += 1.7320508075688772*f[0]*hamil[2]*dv11*jacob_vy_inv*volFact; 
  out[5] += 1.7320508075688772*f[1]*hamil[2]*dv11*jacob_vy_inv*volFact; 
  out[6] += 1.7320508075688772*f[0]*hamil[3]*dv12*jacob_vz_inv*volFact; 
  out[7] += 1.7320508075688772*f[1]*hamil[3]*dv12*jacob_vz_inv*volFact; 
  out[8] += (hamil[3]*f[4]+hamil[2]*f[3]+hamil[1]*f[2]+f[0]*hamil[0])*volFact; 
  out[9] += (hamil[3]*f[8]+hamil[2]*f[6]+hamil[1]*f[5]+hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_sparse_int_five_moments_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]*0.0625; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  double dv12 = 2.0/dxv[3]; 
  const double *jacob_vz = &jacob_vel[4]; 
  const double jacob_vz_inv = 1.0/jacob_vz[0]; 
  out[0] += 4.0*f[0]*volFact; 
  out[1] += 2.4494897427831783*f[0]*hamil[1]*dv10*jacob_vx_inv*volFact; 
  out[2] += 2.4494897427831783*f[0]*hamil[2]*dv11*jacob_vy_inv*volFact; 
  out[3] += 2.4494897427831783*f[0]*hamil[3]*dv12*jacob_vz_inv*volFact; 
  out[4] += (1.4142135623730951*hamil[3]*f[4]+1.4142135623730951*hamil[2]*f[3]+1.4142135623730951*hamil[1]*f[2]+1.4142135623730951*f[0]*hamil[0])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_phase_M1i_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  double dv12 = 2.0/dxv[3]; 
  const double *jacob_vz = &jacob_vel[4]; 
  const double jacob_vz_inv = 1.0/jacob_vz[0]; 
  out[0] += (1.224744871391589*f[13]*hamil[15]+1.224744871391589*f[10]*hamil[14]+1.224744871391589*f[8]*hamil[12]+1.224744871391589*f[6]*hamil[11]+1.224744871391589*f[4]*hamil[9]+1.224744871391589*f[3]*hamil[7]+1.224744871391589*f[1]*hamil[5]+1.224744871391589*f[0]*hamil[2])*dv10*jacob_vx_inv*volFact; 
  out[1] += (1.224744871391589*f[10]*hamil[15]+1.224744871391589*f[13]*hamil[14]+1.224744871391589*f[4]*hamil[12]+1.224744871391589*f[3]*hamil[11]+1.224744871391589*f[8]*hamil[9]+1.224744871391589*f[6]*hamil[7]+1.224744871391589*f[0]*hamil[5]+1.224744871391589*f[1]*hamil[2])*dv10*jacob_vx_inv*volFact; 
  out[2] += (1.224744871391589*f[12]*hamil[15]+1.224744871391589*f[9]*hamil[14]+1.224744871391589*f[8]*hamil[13]+1.224744871391589*f[5]*hamil[11]+1.224744871391589*f[4]*hamil[10]+1.224744871391589*f[2]*hamil[7]+1.224744871391589*f[1]*hamil[6]+1.224744871391589*f[0]*hamil[3])*dv11*jacob_vy_inv*volFact; 
  out[3] += (1.224744871391589*f[9]*hamil[15]+1.224744871391589*f[12]*hamil[14]+1.224744871391589*f[4]*hamil[13]+1.224744871391589*f[2]*hamil[11]+1.224744871391589*f[8]*hamil[10]+1.224744871391589*f[5]*hamil[7]+1.224744871391589*f[0]*hamil[6]+1.224744871391589*f[1]*hamil[3])*dv11*jacob_vy_inv*volFact; 
  out[4] += (1.224744871391589*f[11]*hamil[15]+1.224744871391589*f[7]*hamil[14]+1.224744871391589*f[6]*hamil[13]+1.224744871391589*f[5]*hamil[12]+1.224744871391589*f[3]*hamil[10]+1.224744871391589*f[2]*hamil[9]+1.224744871391589*f[1]*hamil[8]+1.224744871391589*f[0]*hamil[4])*dv12*jacob_vz_inv*volFact; 
  out[5] += (1.224744871391589*f[7]*hamil[15]+1.224744871391589*f[11]*hamil[14]+1.224744871391589*f[3]*hamil[13]+1.224744871391589*f[2]*hamil[12]+1.224744871391589*f[6]*hamil[10]+1.224744871391589*f[5]*hamil[9]+1.224744871391589*f[0]*hamil[8]+1.224744871391589*f[1]*hamil[4])*dv12*jacob_vz_inv*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_phase_M2_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  out[0] += (0.7071067811865475*f[15]*hamil[15]+0.7071067811865475*f[14]*hamil[14]+0.7071067811865475*f[13]*hamil[13]+0.7071067811865475*f[12]*hamil[12]+0.7071067811865475*f[11]*hamil[11]+0.7071067811865475*f[10]*hamil[10]+0.7071067811865475*f[9]*hamil[9]+0.7071067811865475*f[8]*hamil[8]+0.7071067811865475*f[7]*hamil[7]+0.7071067811865475*f[6]*hamil[6]+0.7071067811865475*f[5]*hamil[5]+0.7071067811865475*f[4]*hamil[4]+0.7071067811865475*f[3]*hamil[3]+0.7071067811865475*f[2]*hamil[2]+0.7071067811865475*f[1]*hamil[1]+0.7071067811865475*f[0]*hamil[0])*volFact; 
  out[1] += (0.7071067811865475*f[14]*hamil[15]+0.7071067811865475*hamil[14]*f[15]+0.7071067811865475*f[10]*hamil[13]+0.7071067811865475*hamil[10]*f[13]+0.7071067811865475*f[9]*hamil[12]+0.7071067811865475*hamil[9]*f[12]+0.7071067811865475*f[7]*hamil[11]+0.7071067811865475*hamil[7]*f[11]+0.7071067811865475*f[4]*hamil[8]+0.7071067811865475*hamil[4]*f[8]+0.7071067811865475*f[3]*hamil[6]+0.7071067811865475*hamil[3]*f[6]+0.7071067811865475*f[2]*hamil[5]+0.7071067811865475*hamil[2]*f[5]+0.7071067811865475*f[0]*hamil[1]+0.7071067811865475*hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_phase_five_moments_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  double dv12 = 2.0/dxv[3]; 
  const double *jacob_vz = &jacob_vel[4]; 
  const double jacob_vz_inv = 1.0/jacob_vz[0]; 
  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += 2.8284271247461907*f[1]*volFact; 
  out[2] += (1.224744871391589*f[13]*hamil[15]+1.224744871391589*f[10]*hamil[14]+1.224744871391589*f[8]*hamil[12]+1.224744871391589*f[6]*hamil[11]+1.224744871391589*f[4]*hamil[9]+1.224744871391589*f[3]*hamil[7]+1.224744871391589*f[1]*hamil[5]+1.224744871391589*f[0]*hamil[2])*dv10*jacob_vx_inv*volFact; 
  out[3] += (1.224744871391589*f[10]*hamil[15]+1.224744871391589*f[13]*hamil[14]+1.224744871391589*f[4]*hamil[12]+1.224744871391589*f[3]*hamil[11]+1.224744871391589*f[8]*hamil[9]+1.224744871391589*f[6]*hamil[7]+1.224744871391589*f[0]*hamil[5]+1.224744871391589*f[1]*hamil[2])*dv10*jacob_vx_inv*volFact; 
  out[4] += (1.224744871391589*f[12]*hamil[15]+1.224744871391589*f[9]*hamil[14]+1.224744871391589*f[8]*hamil[13]+1.224744871391589*f[5]*hamil[11]+1.224744871391589*f[4]*hamil[10]+1.224744871391589*f[2]*hamil[7]+1.224744871391589*f[1]*hamil[6]+1.224744871391589*f[0]*hamil[3])*dv11*jacob_vy_inv*volFact; 
  out[5] += (1.224744871391589*f[9]*hamil[15]+1.224744871391589*f[12]*hamil[14]+1.224744871391589*f[4]*hamil[13]+1.224744871391589*f[2]*hamil[11]+1.224744871391589*f[8]*hamil[10]+1.224744871391589*f[5]*hamil[7]+1.224744871391589*f[0]*hamil[6]+1.224744871391589*f[1]*hamil[3])*dv11*jacob_vy_inv*volFact; 
  out[6] += (1.224744871391589*f[11]*hamil[15]+1.224744871391589*f[7]*hamil[14]+1.224744871391589*f[6]*hamil[13]+1.224744871391589*f[5]*hamil[12]+1.224744871391589*f[3]*hamil[10]+1.224744871391589*f[2]*hamil[9]+1.224744871391589*f[1]*hamil[8]+1.224744871391589*f[0]*hamil[4])*dv12*jacob_vz_inv*volFact; 
  out[7] += (1.224744871391589*f[7]*hamil[15]+1.224744871391589*f[11]*hamil[14]+1.224744871391589*f[3]*hamil[13]+1.224744871391589*f[2]*hamil[12]+1.224744871391589*f[6]*hamil[10]+1.224744871391589*f[5]*hamil[9]+1.224744871391589*f[0]*hamil[8]+1.224744871391589*f[1]*hamil[4])*dv12*jacob_vz_inv*volFact; 
  out[8] += (0.7071067811865475*f[15]*hamil[15]+0.7071067811865475*f[14]*hamil[14]+0.7071067811865475*f[13]*hamil[13]+0.7071067811865475*f[12]*hamil[12]+0.7071067811865475*f[11]*hamil[11]+0.7071067811865475*f[10]*hamil[10]+0.7071067811865475*f[9]*hamil[9]+0.7071067811865475*f[8]*hamil[8]+0.7071067811865475*f[7]*hamil[7]+0.7071067811865475*f[6]*hamil[6]+0.7071067811865475*f[5]*hamil[5]+0.7071067811865475*f[4]*hamil[4]+0.7071067811865475*f[3]*hamil[3]+0.7071067811865475*f[2]*hamil[2]+0.7071067811865475*f[1]*hamil[1]+0.7071067811865475*f[0]*hamil[0])*volFact; 
  out[9] += (0.7071067811865475*f[14]*hamil[15]+0.7071067811865475*hamil[14]*f[15]+0.7071067811865475*f[10]*hamil[13]+0.7071067811865475*hamil[10]*f[13]+0.7071067811865475*f[9]*hamil[12]+0.7071067811865475*hamil[9]*f[12]+0.7071067811865475*f[7]*hamil[11]+0.7071067811865475*hamil[7]*f[11]+0.7071067811865475*f[4]*hamil[8]+0.7071067811865475*hamil[4]*f[8]+0.7071067811865475*f[3]*hamil[6]+0.7071067811865475*hamil[3]*f[6]+0.7071067811865475*f[2]*hamil[5]+0.7071067811865475*hamil[2]*f[5]+0.7071067811865475*f[0]*hamil[1]+0.7071067811865475*hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_phase_int_five_moments_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]*0.0625; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  double dv12 = 2.0/dxv[3]; 
  const double *jacob_vz = &jacob_vel[4]; 
  const double jacob_vz_inv = 1.0/jacob_vz[0]; 
  out[0] += 4.0*f[0]*volFact; 
  out[1] += (1.7320508075688772*f[13]*hamil[15]+1.7320508075688772*f[10]*hamil[14]+1.7320508075688772*f[8]*hamil[12]+1.7320508075688772*f[6]*hamil[11]+1.7320508075688772*f[4]*hamil[9]+1.7320508075688772*f[3]*hamil[7]+1.7320508075688772*f[1]*hamil[5]+1.7320508075688772*f[0]*hamil[2])*dv10*jacob_vx_inv*volFact; 
  out[2] += (1.7320508075688772*f[12]*hamil[15]+1.7320508075688772*f[9]*hamil[14]+1.7320508075688772*f[8]*hamil[13]+1.7320508075688772*f[5]*hamil[11]+1.7320508075688772*f[4]*hamil[10]+1.7320508075688772*f[2]*hamil[7]+1.7320508075688772*f[1]*hamil[6]+1.7320508075688772*f[0]*hamil[3])*dv11*jacob_vy_inv*volFact; 
  out[3] += (1.7320508075688772*f[11]*hamil[15]+1.7320508075688772*f[7]*hamil[14]+1.7320508075688772*f[6]*hamil[13]+1.7320508075688772*f[5]*hamil[12]+1.7320508075688772*f[3]*hamil[10]+1.7320508075688772*f[2]*hamil[9]+1.7320508075688772*f[1]*hamil[8]+1.7320508075688772*f[0]*hamil[4])*dv12*jacob_vz_inv*volFact; 
  out[4] += (f[15]*hamil[15]+f[14]*hamil[14]+f[13]*hamil[13]+f[12]*hamil[12]+f[11]*hamil[11]+f[10]*hamil[10]+f[9]*hamil[9]+f[8]*hamil[8]+f[7]*hamil[7]+f[6]*hamil[6]+f[5]*hamil[5]+f[4]*hamil[4]+f[3]*hamil[3]+f[2]*hamil[2]+f[1]*hamil[1]+f[0]*hamil[0])*volFact; 
} 
