#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void mom_vlasov_M0_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
  const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[5]*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M2ij_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
  const double *vmap, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]/4; 
  const double *vmap_vx = &vmap[0]; 
  const double *vmap_vy = &vmap[4]; 
  const double vmap_vx0_sq = vmap_vx[0]*vmap_vx[0]; 
  const double vmap_vx1_sq = vmap_vx[1]*vmap_vx[1]; 
  const double vmap_vy0_sq = vmap_vy[0]*vmap_vy[0]; 
  const double vmap_vy1_sq = vmap_vy[1]*vmap_vy[1]; 
  out[0] += (f[0]*vmap_vx1_sq+f[0]*vmap_vx0_sq+2.0*vmap_vx[0]*vmap_vx[1]*f[3])*volFact; 
  out[1] += (f[1]*vmap_vx1_sq+f[1]*vmap_vx0_sq+2.0*vmap_vx[0]*vmap_vx[1]*f[6])*volFact; 
  out[2] += (f[2]*vmap_vx1_sq+f[2]*vmap_vx0_sq+2.0*vmap_vx[0]*vmap_vx[1]*f[7])*volFact; 
  out[3] += (f[5]*vmap_vx1_sq+f[5]*vmap_vx0_sq+2.0*vmap_vx[0]*vmap_vx[1]*f[11])*volFact; 
  out[4] += (vmap_vx[1]*vmap_vy[1]*f[10]+vmap_vx[0]*vmap_vy[1]*f[4]+vmap_vy[0]*vmap_vx[1]*f[3]+f[0]*vmap_vx[0]*vmap_vy[0])*volFact; 
  out[5] += (vmap_vx[1]*vmap_vy[1]*f[13]+vmap_vx[0]*vmap_vy[1]*f[8]+vmap_vy[0]*vmap_vx[1]*f[6]+vmap_vx[0]*vmap_vy[0]*f[1])*volFact; 
  out[6] += (vmap_vx[1]*vmap_vy[1]*f[14]+vmap_vx[0]*vmap_vy[1]*f[9]+vmap_vy[0]*vmap_vx[1]*f[7]+vmap_vx[0]*vmap_vy[0]*f[2])*volFact; 
  out[7] += (vmap_vx[1]*vmap_vy[1]*f[15]+vmap_vx[0]*vmap_vy[1]*f[12]+vmap_vy[0]*vmap_vx[1]*f[11]+vmap_vx[0]*vmap_vy[0]*f[5])*volFact; 
  out[8] += (f[0]*vmap_vy1_sq+f[0]*vmap_vy0_sq+2.0*vmap_vy[0]*vmap_vy[1]*f[4])*volFact; 
  out[9] += (f[1]*vmap_vy1_sq+f[1]*vmap_vy0_sq+2.0*vmap_vy[0]*vmap_vy[1]*f[8])*volFact; 
  out[10] += (f[2]*vmap_vy1_sq+f[2]*vmap_vy0_sq+2.0*vmap_vy[0]*vmap_vy[1]*f[9])*volFact; 
  out[11] += (f[5]*vmap_vy1_sq+f[5]*vmap_vy0_sq+2.0*vmap_vy[0]*vmap_vy[1]*f[12])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M3ijk_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
  const double *vmap, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]/4; 
  const double *vmap_vx = &vmap[0]; 
  const double *vmap_vy = &vmap[4]; 
  const double vmap_vx0_sq = vmap_vx[0]*vmap_vx[0]; 
  const double vmap_vx0_cu = vmap_vx[0]*vmap_vx[0]*vmap_vx[0]; 
  const double vmap_vx1_sq = vmap_vx[1]*vmap_vx[1]; 
  const double vmap_vx1_cu = vmap_vx[1]*vmap_vx[1]*vmap_vx[1]; 
  const double vmap_vy0_sq = vmap_vy[0]*vmap_vy[0]; 
  const double vmap_vy0_cu = vmap_vy[0]*vmap_vy[0]*vmap_vy[0]; 
  const double vmap_vy1_sq = vmap_vy[1]*vmap_vy[1]; 
  const double vmap_vy1_cu = vmap_vy[1]*vmap_vy[1]*vmap_vy[1]; 
  out[0] += (2.1213203435596424*f[0]*vmap_vx[0]*vmap_vx1_sq+1.2727922061357855*f[3]*vmap_vx1_cu+2.1213203435596424*vmap_vx[1]*f[3]*vmap_vx0_sq+0.7071067811865475*f[0]*vmap_vx0_cu)*volFact; 
  out[1] += (2.1213203435596424*vmap_vx[0]*f[1]*vmap_vx1_sq+1.2727922061357855*f[6]*vmap_vx1_cu+2.1213203435596424*vmap_vx[1]*f[6]*vmap_vx0_sq+0.7071067811865475*f[1]*vmap_vx0_cu)*volFact; 
  out[2] += (2.1213203435596424*vmap_vx[0]*f[2]*vmap_vx1_sq+1.2727922061357855*f[7]*vmap_vx1_cu+2.1213203435596424*vmap_vx[1]*f[7]*vmap_vx0_sq+0.7071067811865475*f[2]*vmap_vx0_cu)*volFact; 
  out[3] += (2.1213203435596424*vmap_vx[0]*f[5]*vmap_vx1_sq+1.2727922061357855*f[11]*vmap_vx1_cu+2.1213203435596424*vmap_vx[1]*f[11]*vmap_vx0_sq+0.7071067811865475*f[5]*vmap_vx0_cu)*volFact; 
  out[4] += (0.7071067811865475*vmap_vy[1]*f[4]*vmap_vx1_sq+0.7071067811865475*f[0]*vmap_vy[0]*vmap_vx1_sq+0.7071067811865475*vmap_vy[1]*f[4]*vmap_vx0_sq+0.7071067811865475*f[0]*vmap_vy[0]*vmap_vx0_sq+1.4142135623730951*vmap_vx[0]*vmap_vx[1]*vmap_vy[1]*f[10]+1.4142135623730951*vmap_vx[0]*vmap_vy[0]*vmap_vx[1]*f[3])*volFact; 
  out[5] += (0.7071067811865475*vmap_vy[1]*f[8]*vmap_vx1_sq+0.7071067811865475*vmap_vy[0]*f[1]*vmap_vx1_sq+0.7071067811865475*vmap_vy[1]*f[8]*vmap_vx0_sq+0.7071067811865475*vmap_vy[0]*f[1]*vmap_vx0_sq+1.4142135623730951*vmap_vx[0]*vmap_vx[1]*vmap_vy[1]*f[13]+1.4142135623730951*vmap_vx[0]*vmap_vy[0]*vmap_vx[1]*f[6])*volFact; 
  out[6] += (0.7071067811865475*vmap_vy[1]*f[9]*vmap_vx1_sq+0.7071067811865475*vmap_vy[0]*f[2]*vmap_vx1_sq+0.7071067811865475*vmap_vy[1]*f[9]*vmap_vx0_sq+0.7071067811865475*vmap_vy[0]*f[2]*vmap_vx0_sq+1.4142135623730951*vmap_vx[0]*vmap_vx[1]*vmap_vy[1]*f[14]+1.4142135623730951*vmap_vx[0]*vmap_vy[0]*vmap_vx[1]*f[7])*volFact; 
  out[7] += (0.7071067811865475*vmap_vy[1]*f[12]*vmap_vx1_sq+0.7071067811865475*vmap_vy[0]*f[5]*vmap_vx1_sq+0.7071067811865475*vmap_vy[1]*f[12]*vmap_vx0_sq+0.7071067811865475*vmap_vy[0]*f[5]*vmap_vx0_sq+1.4142135623730951*vmap_vx[0]*vmap_vx[1]*vmap_vy[1]*f[15]+1.4142135623730951*vmap_vx[0]*vmap_vy[0]*vmap_vx[1]*f[11])*volFact; 
  out[8] += (0.7071067811865475*vmap_vx[1]*f[3]*vmap_vy1_sq+0.7071067811865475*f[0]*vmap_vx[0]*vmap_vy1_sq+0.7071067811865475*vmap_vx[1]*f[3]*vmap_vy0_sq+0.7071067811865475*f[0]*vmap_vx[0]*vmap_vy0_sq+1.4142135623730951*vmap_vy[0]*vmap_vx[1]*vmap_vy[1]*f[10]+1.4142135623730951*vmap_vx[0]*vmap_vy[0]*vmap_vy[1]*f[4])*volFact; 
  out[9] += (0.7071067811865475*vmap_vx[1]*f[6]*vmap_vy1_sq+0.7071067811865475*vmap_vx[0]*f[1]*vmap_vy1_sq+0.7071067811865475*vmap_vx[1]*f[6]*vmap_vy0_sq+0.7071067811865475*vmap_vx[0]*f[1]*vmap_vy0_sq+1.4142135623730951*vmap_vy[0]*vmap_vx[1]*vmap_vy[1]*f[13]+1.4142135623730951*vmap_vx[0]*vmap_vy[0]*vmap_vy[1]*f[8])*volFact; 
  out[10] += (0.7071067811865475*vmap_vx[1]*f[7]*vmap_vy1_sq+0.7071067811865475*vmap_vx[0]*f[2]*vmap_vy1_sq+0.7071067811865475*vmap_vx[1]*f[7]*vmap_vy0_sq+0.7071067811865475*vmap_vx[0]*f[2]*vmap_vy0_sq+1.4142135623730951*vmap_vy[0]*vmap_vx[1]*vmap_vy[1]*f[14]+1.4142135623730951*vmap_vx[0]*vmap_vy[0]*vmap_vy[1]*f[9])*volFact; 
  out[11] += (0.7071067811865475*vmap_vx[1]*f[11]*vmap_vy1_sq+0.7071067811865475*vmap_vx[0]*f[5]*vmap_vy1_sq+0.7071067811865475*vmap_vx[1]*f[11]*vmap_vy0_sq+0.7071067811865475*vmap_vx[0]*f[5]*vmap_vy0_sq+1.4142135623730951*vmap_vy[0]*vmap_vx[1]*vmap_vy[1]*f[15]+1.4142135623730951*vmap_vx[0]*vmap_vy[0]*vmap_vy[1]*f[12])*volFact; 
  out[12] += (2.1213203435596424*f[0]*vmap_vy[0]*vmap_vy1_sq+1.2727922061357855*f[4]*vmap_vy1_cu+2.1213203435596424*vmap_vy[1]*f[4]*vmap_vy0_sq+0.7071067811865475*f[0]*vmap_vy0_cu)*volFact; 
  out[13] += (2.1213203435596424*vmap_vy[0]*f[1]*vmap_vy1_sq+1.2727922061357855*f[8]*vmap_vy1_cu+2.1213203435596424*vmap_vy[1]*f[8]*vmap_vy0_sq+0.7071067811865475*f[1]*vmap_vy0_cu)*volFact; 
  out[14] += (2.1213203435596424*vmap_vy[0]*f[2]*vmap_vy1_sq+1.2727922061357855*f[9]*vmap_vy1_cu+2.1213203435596424*vmap_vy[1]*f[9]*vmap_vy0_sq+0.7071067811865475*f[2]*vmap_vy0_cu)*volFact; 
  out[15] += (2.1213203435596424*vmap_vy[0]*f[5]*vmap_vy1_sq+1.2727922061357855*f[12]*vmap_vy1_cu+2.1213203435596424*vmap_vy[1]*f[12]*vmap_vy0_sq+0.7071067811865475*f[5]*vmap_vy0_cu)*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_M1i_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]/4; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[3]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  out[0] += (1.7320508075688772*hamil[3]*f[4]+1.7320508075688772*f[0]*hamil[1])*dv10*jacob_vx_inv*volFact; 
  out[1] += (1.7320508075688772*hamil[3]*f[8]+1.7320508075688772*f[1]*hamil[1])*dv10*jacob_vx_inv*volFact; 
  out[2] += (1.7320508075688772*hamil[3]*f[9]+1.7320508075688772*hamil[1]*f[2])*dv10*jacob_vx_inv*volFact; 
  out[3] += (1.7320508075688772*hamil[3]*f[12]+1.7320508075688772*hamil[1]*f[5])*dv10*jacob_vx_inv*volFact; 
  out[4] += (1.7320508075688772*f[3]*hamil[3]+1.7320508075688772*f[0]*hamil[2])*dv11*jacob_vy_inv*volFact; 
  out[5] += (1.7320508075688772*hamil[3]*f[6]+1.7320508075688772*f[1]*hamil[2])*dv11*jacob_vy_inv*volFact; 
  out[6] += (1.7320508075688772*hamil[3]*f[7]+1.7320508075688772*f[2]*hamil[2])*dv11*jacob_vy_inv*volFact; 
  out[7] += (1.7320508075688772*hamil[3]*f[11]+1.7320508075688772*hamil[2]*f[5])*dv11*jacob_vy_inv*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_M2_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]/4; 
  out[0] += (hamil[3]*f[10]+hamil[2]*f[4]+hamil[1]*f[3]+f[0]*hamil[0])*volFact; 
  out[1] += (hamil[3]*f[13]+hamil[2]*f[8]+hamil[1]*f[6]+hamil[0]*f[1])*volFact; 
  out[2] += (hamil[3]*f[14]+hamil[2]*f[9]+hamil[1]*f[7]+hamil[0]*f[2])*volFact; 
  out[3] += (hamil[3]*f[15]+hamil[2]*f[12]+hamil[1]*f[11]+hamil[0]*f[5])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_M3i_2x2v_ser_p1(const double *w, const double *dxv, const int *idx,
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]/4; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[3]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  double dH_dpx[4] = {0.0}; 
  dH_dpx[0] = 1.7320508075688772*hamil[1]*dv10*jacob_vx_inv; 
  dH_dpx[2] = 1.7320508075688772*hamil[3]*dv10*jacob_vx_inv; 

  double dH_dpy[4] = {0.0}; 
  dH_dpy[0] = 1.7320508075688772*hamil[2]*dv11*jacob_vy_inv; 
  dH_dpy[1] = 1.7320508075688772*hamil[3]*dv11*jacob_vy_inv; 

  out[0] += (0.5*dH_dpx[0]*hamil[3]*f[10]+0.5*hamil[1]*dH_dpx[2]*f[10]+0.5*dH_dpx[0]*hamil[2]*f[4]+0.5*hamil[0]*dH_dpx[2]*f[4]+0.5*dH_dpx[2]*f[3]*hamil[3]+0.5*dH_dpx[0]*hamil[1]*f[3]+0.5*f[0]*dH_dpx[2]*hamil[2]+0.5*dH_dpx[0]*f[0]*hamil[0])*volFact; 
  out[1] += (0.5*dH_dpx[0]*hamil[3]*f[13]+0.5*hamil[1]*dH_dpx[2]*f[13]+0.5*dH_dpx[0]*hamil[2]*f[8]+0.5*hamil[0]*dH_dpx[2]*f[8]+0.5*dH_dpx[2]*hamil[3]*f[6]+0.5*dH_dpx[0]*hamil[1]*f[6]+0.5*f[1]*dH_dpx[2]*hamil[2]+0.5*dH_dpx[0]*hamil[0]*f[1])*volFact; 
  out[2] += (0.5*dH_dpx[0]*hamil[3]*f[14]+0.5*hamil[1]*dH_dpx[2]*f[14]+0.5*dH_dpx[0]*hamil[2]*f[9]+0.5*hamil[0]*dH_dpx[2]*f[9]+0.5*dH_dpx[2]*hamil[3]*f[7]+0.5*dH_dpx[0]*hamil[1]*f[7]+0.5*dH_dpx[2]*f[2]*hamil[2]+0.5*dH_dpx[0]*hamil[0]*f[2])*volFact; 
  out[3] += (0.5*dH_dpx[0]*hamil[3]*f[15]+0.5*hamil[1]*dH_dpx[2]*f[15]+0.5*dH_dpx[0]*hamil[2]*f[12]+0.5*hamil[0]*dH_dpx[2]*f[12]+0.5*dH_dpx[2]*hamil[3]*f[11]+0.5*dH_dpx[0]*hamil[1]*f[11]+0.5*dH_dpx[2]*hamil[2]*f[5]+0.5*dH_dpx[0]*hamil[0]*f[5])*volFact; 
  out[4] += (0.5*dH_dpy[0]*hamil[3]*f[10]+0.5*dH_dpy[1]*hamil[2]*f[10]+0.5*dH_dpy[1]*hamil[3]*f[4]+0.5*dH_dpy[0]*hamil[2]*f[4]+0.5*dH_dpy[0]*hamil[1]*f[3]+0.5*hamil[0]*dH_dpy[1]*f[3]+0.5*f[0]*dH_dpy[1]*hamil[1]+0.5*dH_dpy[0]*f[0]*hamil[0])*volFact; 
  out[5] += (0.5*dH_dpy[0]*hamil[3]*f[13]+0.5*dH_dpy[1]*hamil[2]*f[13]+0.5*dH_dpy[1]*hamil[3]*f[8]+0.5*dH_dpy[0]*hamil[2]*f[8]+0.5*dH_dpy[0]*hamil[1]*f[6]+0.5*hamil[0]*dH_dpy[1]*f[6]+0.5*dH_dpy[1]*f[1]*hamil[1]+0.5*dH_dpy[0]*hamil[0]*f[1])*volFact; 
  out[6] += (0.5*dH_dpy[0]*hamil[3]*f[14]+0.5*dH_dpy[1]*hamil[2]*f[14]+0.5*dH_dpy[1]*hamil[3]*f[9]+0.5*dH_dpy[0]*hamil[2]*f[9]+0.5*dH_dpy[0]*hamil[1]*f[7]+0.5*hamil[0]*dH_dpy[1]*f[7]+0.5*dH_dpy[1]*hamil[1]*f[2]+0.5*dH_dpy[0]*hamil[0]*f[2])*volFact; 
  out[7] += (0.5*dH_dpy[0]*hamil[3]*f[15]+0.5*dH_dpy[1]*hamil[2]*f[15]+0.5*dH_dpy[1]*hamil[3]*f[12]+0.5*dH_dpy[0]*hamil[2]*f[12]+0.5*dH_dpy[0]*hamil[1]*f[11]+0.5*hamil[0]*dH_dpy[1]*f[11]+0.5*dH_dpy[1]*hamil[1]*f[5]+0.5*dH_dpy[0]*hamil[0]*f[5])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_five_moments_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]/4; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[3]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[5]*volFact; 
  out[4] += (1.7320508075688772*hamil[3]*f[4]+1.7320508075688772*f[0]*hamil[1])*dv10*jacob_vx_inv*volFact; 
  out[5] += (1.7320508075688772*hamil[3]*f[8]+1.7320508075688772*f[1]*hamil[1])*dv10*jacob_vx_inv*volFact; 
  out[6] += (1.7320508075688772*hamil[3]*f[9]+1.7320508075688772*hamil[1]*f[2])*dv10*jacob_vx_inv*volFact; 
  out[7] += (1.7320508075688772*hamil[3]*f[12]+1.7320508075688772*hamil[1]*f[5])*dv10*jacob_vx_inv*volFact; 
  out[8] += (1.7320508075688772*f[3]*hamil[3]+1.7320508075688772*f[0]*hamil[2])*dv11*jacob_vy_inv*volFact; 
  out[9] += (1.7320508075688772*hamil[3]*f[6]+1.7320508075688772*f[1]*hamil[2])*dv11*jacob_vy_inv*volFact; 
  out[10] += (1.7320508075688772*hamil[3]*f[7]+1.7320508075688772*f[2]*hamil[2])*dv11*jacob_vy_inv*volFact; 
  out[11] += (1.7320508075688772*hamil[3]*f[11]+1.7320508075688772*hamil[2]*f[5])*dv11*jacob_vy_inv*volFact; 
  out[12] += (hamil[3]*f[10]+hamil[2]*f[4]+hamil[1]*f[3]+f[0]*hamil[0])*volFact; 
  out[13] += (hamil[3]*f[13]+hamil[2]*f[8]+hamil[1]*f[6]+hamil[0]*f[1])*volFact; 
  out[14] += (hamil[3]*f[14]+hamil[2]*f[9]+hamil[1]*f[7]+hamil[0]*f[2])*volFact; 
  out[15] += (hamil[3]*f[15]+hamil[2]*f[12]+hamil[1]*f[11]+hamil[0]*f[5])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_int_five_moments_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]*0.0625; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[3]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  out[0] += 4.0*f[0]*volFact; 
  out[1] += (3.4641016151377544*hamil[3]*f[4]+3.4641016151377544*f[0]*hamil[1])*dv10*jacob_vx_inv*volFact; 
  out[2] += (3.4641016151377544*f[3]*hamil[3]+3.4641016151377544*f[0]*hamil[2])*dv11*jacob_vy_inv*volFact; 
  out[3] += (2.0*hamil[3]*f[10]+2.0*hamil[2]*f[4]+2.0*hamil[1]*f[3]+2.0*f[0]*hamil[0])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_sparse_M1i_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]/4; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[3]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  out[0] += 1.7320508075688772*f[0]*hamil[1]*dv10*jacob_vx_inv*volFact; 
  out[1] += 1.7320508075688772*f[1]*hamil[1]*dv10*jacob_vx_inv*volFact; 
  out[2] += 1.7320508075688772*hamil[1]*f[2]*dv10*jacob_vx_inv*volFact; 
  out[3] += 1.7320508075688772*hamil[1]*f[5]*dv10*jacob_vx_inv*volFact; 
  out[4] += 1.7320508075688772*f[0]*hamil[2]*dv11*jacob_vy_inv*volFact; 
  out[5] += 1.7320508075688772*f[1]*hamil[2]*dv11*jacob_vy_inv*volFact; 
  out[6] += 1.7320508075688772*f[2]*hamil[2]*dv11*jacob_vy_inv*volFact; 
  out[7] += 1.7320508075688772*hamil[2]*f[5]*dv11*jacob_vy_inv*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_sparse_M2_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]/4; 
  out[0] += (hamil[2]*f[4]+hamil[1]*f[3]+f[0]*hamil[0])*volFact; 
  out[1] += (hamil[2]*f[8]+hamil[1]*f[6]+hamil[0]*f[1])*volFact; 
  out[2] += (hamil[2]*f[9]+hamil[1]*f[7]+hamil[0]*f[2])*volFact; 
  out[3] += (hamil[2]*f[12]+hamil[1]*f[11]+hamil[0]*f[5])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_sparse_M3i_2x2v_ser_p1(const double *w, const double *dxv, const int *idx,
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]/4; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[3]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  double dH_dpx[4] = {0.0}; 
  dH_dpx[0] = 1.7320508075688772*hamil[1]*dv10*jacob_vx_inv; 

  double dH_dpy[4] = {0.0}; 
  dH_dpy[0] = 1.7320508075688772*hamil[2]*dv11*jacob_vy_inv; 

  out[0] += (0.5*dH_dpx[0]*hamil[2]*f[4]+0.5*dH_dpx[0]*hamil[1]*f[3]+0.5*dH_dpx[0]*f[0]*hamil[0])*volFact; 
  out[1] += (0.5*dH_dpx[0]*hamil[2]*f[8]+0.5*dH_dpx[0]*hamil[1]*f[6]+0.5*dH_dpx[0]*hamil[0]*f[1])*volFact; 
  out[2] += (0.5*dH_dpx[0]*hamil[2]*f[9]+0.5*dH_dpx[0]*hamil[1]*f[7]+0.5*dH_dpx[0]*hamil[0]*f[2])*volFact; 
  out[3] += (0.5*dH_dpx[0]*hamil[2]*f[12]+0.5*dH_dpx[0]*hamil[1]*f[11]+0.5*dH_dpx[0]*hamil[0]*f[5])*volFact; 
  out[4] += (0.5*dH_dpy[0]*hamil[2]*f[4]+0.5*dH_dpy[0]*hamil[1]*f[3]+0.5*dH_dpy[0]*f[0]*hamil[0])*volFact; 
  out[5] += (0.5*dH_dpy[0]*hamil[2]*f[8]+0.5*dH_dpy[0]*hamil[1]*f[6]+0.5*dH_dpy[0]*hamil[0]*f[1])*volFact; 
  out[6] += (0.5*dH_dpy[0]*hamil[2]*f[9]+0.5*dH_dpy[0]*hamil[1]*f[7]+0.5*dH_dpy[0]*hamil[0]*f[2])*volFact; 
  out[7] += (0.5*dH_dpy[0]*hamil[2]*f[12]+0.5*dH_dpy[0]*hamil[1]*f[11]+0.5*dH_dpy[0]*hamil[0]*f[5])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_sparse_five_moments_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]/4; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[3]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[5]*volFact; 
  out[4] += 1.7320508075688772*f[0]*hamil[1]*dv10*jacob_vx_inv*volFact; 
  out[5] += 1.7320508075688772*f[1]*hamil[1]*dv10*jacob_vx_inv*volFact; 
  out[6] += 1.7320508075688772*hamil[1]*f[2]*dv10*jacob_vx_inv*volFact; 
  out[7] += 1.7320508075688772*hamil[1]*f[5]*dv10*jacob_vx_inv*volFact; 
  out[8] += 1.7320508075688772*f[0]*hamil[2]*dv11*jacob_vy_inv*volFact; 
  out[9] += 1.7320508075688772*f[1]*hamil[2]*dv11*jacob_vy_inv*volFact; 
  out[10] += 1.7320508075688772*f[2]*hamil[2]*dv11*jacob_vy_inv*volFact; 
  out[11] += 1.7320508075688772*hamil[2]*f[5]*dv11*jacob_vy_inv*volFact; 
  out[12] += (hamil[2]*f[4]+hamil[1]*f[3]+f[0]*hamil[0])*volFact; 
  out[13] += (hamil[2]*f[8]+hamil[1]*f[6]+hamil[0]*f[1])*volFact; 
  out[14] += (hamil[2]*f[9]+hamil[1]*f[7]+hamil[0]*f[2])*volFact; 
  out[15] += (hamil[2]*f[12]+hamil[1]*f[11]+hamil[0]*f[5])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_sparse_int_five_moments_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]*0.0625; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[3]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  out[0] += 4.0*f[0]*volFact; 
  out[1] += 3.4641016151377544*f[0]*hamil[1]*dv10*jacob_vx_inv*volFact; 
  out[2] += 3.4641016151377544*f[0]*hamil[2]*dv11*jacob_vy_inv*volFact; 
  out[3] += (2.0*hamil[2]*f[4]+2.0*hamil[1]*f[3]+2.0*f[0]*hamil[0])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_phase_M1i_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]/4; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[3]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  out[0] += (0.8660254037844386*f[12]*hamil[15]+0.8660254037844386*f[9]*hamil[14]+0.8660254037844386*f[8]*hamil[13]+0.8660254037844386*f[5]*hamil[11]+0.8660254037844386*f[4]*hamil[10]+0.8660254037844386*f[2]*hamil[7]+0.8660254037844386*f[1]*hamil[6]+0.8660254037844386*f[0]*hamil[3])*dv10*jacob_vx_inv*volFact; 
  out[1] += (0.8660254037844386*f[9]*hamil[15]+0.8660254037844386*f[12]*hamil[14]+0.8660254037844386*f[4]*hamil[13]+0.8660254037844386*f[2]*hamil[11]+0.8660254037844386*f[8]*hamil[10]+0.8660254037844386*f[5]*hamil[7]+0.8660254037844386*f[0]*hamil[6]+0.8660254037844386*f[1]*hamil[3])*dv10*jacob_vx_inv*volFact; 
  out[2] += (0.8660254037844386*f[8]*hamil[15]+0.8660254037844386*f[4]*hamil[14]+0.8660254037844386*f[12]*hamil[13]+0.8660254037844386*f[1]*hamil[11]+0.8660254037844386*f[9]*hamil[10]+0.8660254037844386*f[0]*hamil[7]+0.8660254037844386*f[5]*hamil[6]+0.8660254037844386*f[2]*hamil[3])*dv10*jacob_vx_inv*volFact; 
  out[3] += (0.8660254037844386*f[4]*hamil[15]+0.8660254037844386*f[8]*hamil[14]+0.8660254037844386*f[9]*hamil[13]+0.8660254037844386*hamil[10]*f[12]+0.8660254037844386*f[0]*hamil[11]+0.8660254037844386*f[1]*hamil[7]+0.8660254037844386*f[2]*hamil[6]+0.8660254037844386*hamil[3]*f[5])*dv10*jacob_vx_inv*volFact; 
  out[4] += (0.8660254037844386*f[11]*hamil[15]+0.8660254037844386*f[7]*hamil[14]+0.8660254037844386*f[6]*hamil[13]+0.8660254037844386*f[5]*hamil[12]+0.8660254037844386*f[3]*hamil[10]+0.8660254037844386*f[2]*hamil[9]+0.8660254037844386*f[1]*hamil[8]+0.8660254037844386*f[0]*hamil[4])*dv11*jacob_vy_inv*volFact; 
  out[5] += (0.8660254037844386*f[7]*hamil[15]+0.8660254037844386*f[11]*hamil[14]+0.8660254037844386*f[3]*hamil[13]+0.8660254037844386*f[2]*hamil[12]+0.8660254037844386*f[6]*hamil[10]+0.8660254037844386*f[5]*hamil[9]+0.8660254037844386*f[0]*hamil[8]+0.8660254037844386*f[1]*hamil[4])*dv11*jacob_vy_inv*volFact; 
  out[6] += (0.8660254037844386*f[6]*hamil[15]+0.8660254037844386*f[3]*hamil[14]+0.8660254037844386*f[11]*hamil[13]+0.8660254037844386*f[1]*hamil[12]+0.8660254037844386*f[7]*hamil[10]+0.8660254037844386*f[0]*hamil[9]+0.8660254037844386*f[5]*hamil[8]+0.8660254037844386*f[2]*hamil[4])*dv11*jacob_vy_inv*volFact; 
  out[7] += (0.8660254037844386*f[3]*hamil[15]+0.8660254037844386*f[6]*hamil[14]+0.8660254037844386*f[7]*hamil[13]+0.8660254037844386*f[0]*hamil[12]+0.8660254037844386*hamil[10]*f[11]+0.8660254037844386*f[1]*hamil[9]+0.8660254037844386*f[2]*hamil[8]+0.8660254037844386*hamil[4]*f[5])*dv11*jacob_vy_inv*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_phase_M2_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]/4; 
  out[0] += (0.5*f[15]*hamil[15]+0.5*f[14]*hamil[14]+0.5*f[13]*hamil[13]+0.5*f[12]*hamil[12]+0.5*f[11]*hamil[11]+0.5*f[10]*hamil[10]+0.5*f[9]*hamil[9]+0.5*f[8]*hamil[8]+0.5*f[7]*hamil[7]+0.5*f[6]*hamil[6]+0.5*f[5]*hamil[5]+0.5*f[4]*hamil[4]+0.5*f[3]*hamil[3]+0.5*f[2]*hamil[2]+0.5*f[1]*hamil[1]+0.5*f[0]*hamil[0])*volFact; 
  out[1] += (0.5*f[14]*hamil[15]+0.5*hamil[14]*f[15]+0.5*f[10]*hamil[13]+0.5*hamil[10]*f[13]+0.5*f[9]*hamil[12]+0.5*hamil[9]*f[12]+0.5*f[7]*hamil[11]+0.5*hamil[7]*f[11]+0.5*f[4]*hamil[8]+0.5*hamil[4]*f[8]+0.5*f[3]*hamil[6]+0.5*hamil[3]*f[6]+0.5*f[2]*hamil[5]+0.5*hamil[2]*f[5]+0.5*f[0]*hamil[1]+0.5*hamil[0]*f[1])*volFact; 
  out[2] += (0.5*f[13]*hamil[15]+0.5*hamil[13]*f[15]+0.5*f[10]*hamil[14]+0.5*hamil[10]*f[14]+0.5*f[8]*hamil[12]+0.5*hamil[8]*f[12]+0.5*f[6]*hamil[11]+0.5*hamil[6]*f[11]+0.5*f[4]*hamil[9]+0.5*hamil[4]*f[9]+0.5*f[3]*hamil[7]+0.5*hamil[3]*f[7]+0.5*f[1]*hamil[5]+0.5*hamil[1]*f[5]+0.5*f[0]*hamil[2]+0.5*hamil[0]*f[2])*volFact; 
  out[3] += (0.5*f[10]*hamil[15]+0.5*hamil[10]*f[15]+0.5*f[13]*hamil[14]+0.5*hamil[13]*f[14]+0.5*f[4]*hamil[12]+0.5*hamil[4]*f[12]+0.5*f[3]*hamil[11]+0.5*hamil[3]*f[11]+0.5*f[8]*hamil[9]+0.5*hamil[8]*f[9]+0.5*f[6]*hamil[7]+0.5*hamil[6]*f[7]+0.5*f[0]*hamil[5]+0.5*hamil[0]*f[5]+0.5*f[1]*hamil[2]+0.5*hamil[1]*f[2])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_phase_five_moments_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]/4; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[3]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[5]*volFact; 
  out[4] += (0.8660254037844386*f[12]*hamil[15]+0.8660254037844386*f[9]*hamil[14]+0.8660254037844386*f[8]*hamil[13]+0.8660254037844386*f[5]*hamil[11]+0.8660254037844386*f[4]*hamil[10]+0.8660254037844386*f[2]*hamil[7]+0.8660254037844386*f[1]*hamil[6]+0.8660254037844386*f[0]*hamil[3])*dv10*jacob_vx_inv*volFact; 
  out[5] += (0.8660254037844386*f[9]*hamil[15]+0.8660254037844386*f[12]*hamil[14]+0.8660254037844386*f[4]*hamil[13]+0.8660254037844386*f[2]*hamil[11]+0.8660254037844386*f[8]*hamil[10]+0.8660254037844386*f[5]*hamil[7]+0.8660254037844386*f[0]*hamil[6]+0.8660254037844386*f[1]*hamil[3])*dv10*jacob_vx_inv*volFact; 
  out[6] += (0.8660254037844386*f[8]*hamil[15]+0.8660254037844386*f[4]*hamil[14]+0.8660254037844386*f[12]*hamil[13]+0.8660254037844386*f[1]*hamil[11]+0.8660254037844386*f[9]*hamil[10]+0.8660254037844386*f[0]*hamil[7]+0.8660254037844386*f[5]*hamil[6]+0.8660254037844386*f[2]*hamil[3])*dv10*jacob_vx_inv*volFact; 
  out[7] += (0.8660254037844386*f[4]*hamil[15]+0.8660254037844386*f[8]*hamil[14]+0.8660254037844386*f[9]*hamil[13]+0.8660254037844386*hamil[10]*f[12]+0.8660254037844386*f[0]*hamil[11]+0.8660254037844386*f[1]*hamil[7]+0.8660254037844386*f[2]*hamil[6]+0.8660254037844386*hamil[3]*f[5])*dv10*jacob_vx_inv*volFact; 
  out[8] += (0.8660254037844386*f[11]*hamil[15]+0.8660254037844386*f[7]*hamil[14]+0.8660254037844386*f[6]*hamil[13]+0.8660254037844386*f[5]*hamil[12]+0.8660254037844386*f[3]*hamil[10]+0.8660254037844386*f[2]*hamil[9]+0.8660254037844386*f[1]*hamil[8]+0.8660254037844386*f[0]*hamil[4])*dv11*jacob_vy_inv*volFact; 
  out[9] += (0.8660254037844386*f[7]*hamil[15]+0.8660254037844386*f[11]*hamil[14]+0.8660254037844386*f[3]*hamil[13]+0.8660254037844386*f[2]*hamil[12]+0.8660254037844386*f[6]*hamil[10]+0.8660254037844386*f[5]*hamil[9]+0.8660254037844386*f[0]*hamil[8]+0.8660254037844386*f[1]*hamil[4])*dv11*jacob_vy_inv*volFact; 
  out[10] += (0.8660254037844386*f[6]*hamil[15]+0.8660254037844386*f[3]*hamil[14]+0.8660254037844386*f[11]*hamil[13]+0.8660254037844386*f[1]*hamil[12]+0.8660254037844386*f[7]*hamil[10]+0.8660254037844386*f[0]*hamil[9]+0.8660254037844386*f[5]*hamil[8]+0.8660254037844386*f[2]*hamil[4])*dv11*jacob_vy_inv*volFact; 
  out[11] += (0.8660254037844386*f[3]*hamil[15]+0.8660254037844386*f[6]*hamil[14]+0.8660254037844386*f[7]*hamil[13]+0.8660254037844386*f[0]*hamil[12]+0.8660254037844386*hamil[10]*f[11]+0.8660254037844386*f[1]*hamil[9]+0.8660254037844386*f[2]*hamil[8]+0.8660254037844386*hamil[4]*f[5])*dv11*jacob_vy_inv*volFact; 
  out[12] += (0.5*f[15]*hamil[15]+0.5*f[14]*hamil[14]+0.5*f[13]*hamil[13]+0.5*f[12]*hamil[12]+0.5*f[11]*hamil[11]+0.5*f[10]*hamil[10]+0.5*f[9]*hamil[9]+0.5*f[8]*hamil[8]+0.5*f[7]*hamil[7]+0.5*f[6]*hamil[6]+0.5*f[5]*hamil[5]+0.5*f[4]*hamil[4]+0.5*f[3]*hamil[3]+0.5*f[2]*hamil[2]+0.5*f[1]*hamil[1]+0.5*f[0]*hamil[0])*volFact; 
  out[13] += (0.5*f[14]*hamil[15]+0.5*hamil[14]*f[15]+0.5*f[10]*hamil[13]+0.5*hamil[10]*f[13]+0.5*f[9]*hamil[12]+0.5*hamil[9]*f[12]+0.5*f[7]*hamil[11]+0.5*hamil[7]*f[11]+0.5*f[4]*hamil[8]+0.5*hamil[4]*f[8]+0.5*f[3]*hamil[6]+0.5*hamil[3]*f[6]+0.5*f[2]*hamil[5]+0.5*hamil[2]*f[5]+0.5*f[0]*hamil[1]+0.5*hamil[0]*f[1])*volFact; 
  out[14] += (0.5*f[13]*hamil[15]+0.5*hamil[13]*f[15]+0.5*f[10]*hamil[14]+0.5*hamil[10]*f[14]+0.5*f[8]*hamil[12]+0.5*hamil[8]*f[12]+0.5*f[6]*hamil[11]+0.5*hamil[6]*f[11]+0.5*f[4]*hamil[9]+0.5*hamil[4]*f[9]+0.5*f[3]*hamil[7]+0.5*hamil[3]*f[7]+0.5*f[1]*hamil[5]+0.5*hamil[1]*f[5]+0.5*f[0]*hamil[2]+0.5*hamil[0]*f[2])*volFact; 
  out[15] += (0.5*f[10]*hamil[15]+0.5*hamil[10]*f[15]+0.5*f[13]*hamil[14]+0.5*hamil[13]*f[14]+0.5*f[4]*hamil[12]+0.5*hamil[4]*f[12]+0.5*f[3]*hamil[11]+0.5*hamil[3]*f[11]+0.5*f[8]*hamil[9]+0.5*hamil[8]*f[9]+0.5*f[6]*hamil[7]+0.5*hamil[6]*f[7]+0.5*f[0]*hamil[5]+0.5*hamil[0]*f[5]+0.5*f[1]*hamil[2]+0.5*hamil[1]*f[2])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_phase_int_five_moments_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]*0.0625; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[3]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  out[0] += 4.0*f[0]*volFact; 
  out[1] += (1.7320508075688772*f[12]*hamil[15]+1.7320508075688772*f[9]*hamil[14]+1.7320508075688772*f[8]*hamil[13]+1.7320508075688772*f[5]*hamil[11]+1.7320508075688772*f[4]*hamil[10]+1.7320508075688772*f[2]*hamil[7]+1.7320508075688772*f[1]*hamil[6]+1.7320508075688772*f[0]*hamil[3])*dv10*jacob_vx_inv*volFact; 
  out[2] += (1.7320508075688772*f[11]*hamil[15]+1.7320508075688772*f[7]*hamil[14]+1.7320508075688772*f[6]*hamil[13]+1.7320508075688772*f[5]*hamil[12]+1.7320508075688772*f[3]*hamil[10]+1.7320508075688772*f[2]*hamil[9]+1.7320508075688772*f[1]*hamil[8]+1.7320508075688772*f[0]*hamil[4])*dv11*jacob_vy_inv*volFact; 
  out[3] += (f[15]*hamil[15]+f[14]*hamil[14]+f[13]*hamil[13]+f[12]*hamil[12]+f[11]*hamil[11]+f[10]*hamil[10]+f[9]*hamil[9]+f[8]*hamil[8]+f[7]*hamil[7]+f[6]*hamil[6]+f[5]*hamil[5]+f[4]*hamil[4]+f[3]*hamil[3]+f[2]*hamil[2]+f[1]*hamil[1]+f[0]*hamil[0])*volFact; 
} 
