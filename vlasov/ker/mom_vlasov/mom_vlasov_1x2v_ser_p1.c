#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void mom_vlasov_M0_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
  const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M2ij_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
  const double *vmap, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  const double *vmap_vx = &vmap[0]; 
  const double *vmap_vy = &vmap[4]; 
  const double vmap_vx0_sq = vmap_vx[0]*vmap_vx[0]; 
  const double vmap_vx1_sq = vmap_vx[1]*vmap_vx[1]; 
  const double vmap_vy0_sq = vmap_vy[0]*vmap_vy[0]; 
  const double vmap_vy1_sq = vmap_vy[1]*vmap_vy[1]; 
  out[0] += (f[0]*vmap_vx1_sq+f[0]*vmap_vx0_sq+2.0*vmap_vx[0]*vmap_vx[1]*f[2])*volFact; 
  out[1] += (f[1]*vmap_vx1_sq+f[1]*vmap_vx0_sq+2.0*vmap_vx[0]*vmap_vx[1]*f[4])*volFact; 
  out[2] += (vmap_vx[1]*vmap_vy[1]*f[6]+vmap_vx[0]*vmap_vy[1]*f[3]+vmap_vy[0]*vmap_vx[1]*f[2]+f[0]*vmap_vx[0]*vmap_vy[0])*volFact; 
  out[3] += (vmap_vx[1]*vmap_vy[1]*f[7]+vmap_vx[0]*vmap_vy[1]*f[5]+vmap_vy[0]*vmap_vx[1]*f[4]+vmap_vx[0]*vmap_vy[0]*f[1])*volFact; 
  out[4] += (f[0]*vmap_vy1_sq+f[0]*vmap_vy0_sq+2.0*vmap_vy[0]*vmap_vy[1]*f[3])*volFact; 
  out[5] += (f[1]*vmap_vy1_sq+f[1]*vmap_vy0_sq+2.0*vmap_vy[0]*vmap_vy[1]*f[5])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M3ijk_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
  const double *vmap, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
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
  out[0] += (2.1213203435596424*f[0]*vmap_vx[0]*vmap_vx1_sq+1.2727922061357855*f[2]*vmap_vx1_cu+2.1213203435596424*vmap_vx[1]*f[2]*vmap_vx0_sq+0.7071067811865475*f[0]*vmap_vx0_cu)*volFact; 
  out[1] += (2.1213203435596424*vmap_vx[0]*f[1]*vmap_vx1_sq+1.2727922061357855*f[4]*vmap_vx1_cu+2.1213203435596424*vmap_vx[1]*f[4]*vmap_vx0_sq+0.7071067811865475*f[1]*vmap_vx0_cu)*volFact; 
  out[2] += (0.7071067811865475*vmap_vy[1]*f[3]*vmap_vx1_sq+0.7071067811865475*f[0]*vmap_vy[0]*vmap_vx1_sq+0.7071067811865475*vmap_vy[1]*f[3]*vmap_vx0_sq+0.7071067811865475*f[0]*vmap_vy[0]*vmap_vx0_sq+1.4142135623730951*vmap_vx[0]*vmap_vx[1]*vmap_vy[1]*f[6]+1.4142135623730951*vmap_vx[0]*vmap_vy[0]*vmap_vx[1]*f[2])*volFact; 
  out[3] += (0.7071067811865475*vmap_vy[1]*f[5]*vmap_vx1_sq+0.7071067811865475*vmap_vy[0]*f[1]*vmap_vx1_sq+0.7071067811865475*vmap_vy[1]*f[5]*vmap_vx0_sq+0.7071067811865475*vmap_vy[0]*f[1]*vmap_vx0_sq+1.4142135623730951*vmap_vx[0]*vmap_vx[1]*vmap_vy[1]*f[7]+1.4142135623730951*vmap_vx[0]*vmap_vy[0]*vmap_vx[1]*f[4])*volFact; 
  out[4] += (0.7071067811865475*vmap_vx[1]*f[2]*vmap_vy1_sq+0.7071067811865475*f[0]*vmap_vx[0]*vmap_vy1_sq+0.7071067811865475*vmap_vx[1]*f[2]*vmap_vy0_sq+0.7071067811865475*f[0]*vmap_vx[0]*vmap_vy0_sq+1.4142135623730951*vmap_vy[0]*vmap_vx[1]*vmap_vy[1]*f[6]+1.4142135623730951*vmap_vx[0]*vmap_vy[0]*vmap_vy[1]*f[3])*volFact; 
  out[5] += (0.7071067811865475*vmap_vx[1]*f[4]*vmap_vy1_sq+0.7071067811865475*vmap_vx[0]*f[1]*vmap_vy1_sq+0.7071067811865475*vmap_vx[1]*f[4]*vmap_vy0_sq+0.7071067811865475*vmap_vx[0]*f[1]*vmap_vy0_sq+1.4142135623730951*vmap_vy[0]*vmap_vx[1]*vmap_vy[1]*f[7]+1.4142135623730951*vmap_vx[0]*vmap_vy[0]*vmap_vy[1]*f[5])*volFact; 
  out[6] += (2.1213203435596424*f[0]*vmap_vy[0]*vmap_vy1_sq+1.2727922061357855*f[3]*vmap_vy1_cu+2.1213203435596424*vmap_vy[1]*f[3]*vmap_vy0_sq+0.7071067811865475*f[0]*vmap_vy0_cu)*volFact; 
  out[7] += (2.1213203435596424*vmap_vy[0]*f[1]*vmap_vy1_sq+1.2727922061357855*f[5]*vmap_vy1_cu+2.1213203435596424*vmap_vy[1]*f[5]*vmap_vy0_sq+0.7071067811865475*f[1]*vmap_vy0_cu)*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_M1i_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  out[0] += (1.7320508075688772*f[3]*hamil[3]+1.7320508075688772*f[0]*hamil[1])*dv10*jacob_vx_inv*volFact; 
  out[1] += (1.7320508075688772*hamil[3]*f[5]+1.7320508075688772*f[1]*hamil[1])*dv10*jacob_vx_inv*volFact; 
  out[2] += (1.7320508075688772*f[2]*hamil[3]+1.7320508075688772*f[0]*hamil[2])*dv11*jacob_vy_inv*volFact; 
  out[3] += (1.7320508075688772*hamil[3]*f[4]+1.7320508075688772*f[1]*hamil[2])*dv11*jacob_vy_inv*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_M2_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  out[0] += (hamil[3]*f[6]+hamil[2]*f[3]+hamil[1]*f[2]+f[0]*hamil[0])*volFact; 
  out[1] += (hamil[3]*f[7]+hamil[2]*f[5]+hamil[1]*f[4]+hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_M3i_1x2v_ser_p1(const double *w, const double *dxv, const int *idx,
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  double dH_dpx[4] = {0.0}; 
  dH_dpx[0] = 1.7320508075688772*hamil[1]*dv10*jacob_vx_inv; 
  dH_dpx[2] = 1.7320508075688772*hamil[3]*dv10*jacob_vx_inv; 

  double dH_dpy[4] = {0.0}; 
  dH_dpy[0] = 1.7320508075688772*hamil[2]*dv11*jacob_vy_inv; 
  dH_dpy[1] = 1.7320508075688772*hamil[3]*dv11*jacob_vy_inv; 

  out[0] += (0.5*dH_dpx[0]*hamil[3]*f[6]+0.5*hamil[1]*dH_dpx[2]*f[6]+0.5*dH_dpx[2]*f[2]*hamil[3]+0.5*dH_dpx[0]*hamil[2]*f[3]+0.5*hamil[0]*dH_dpx[2]*f[3]+0.5*f[0]*dH_dpx[2]*hamil[2]+0.5*dH_dpx[0]*hamil[1]*f[2]+0.5*dH_dpx[0]*f[0]*hamil[0])*volFact; 
  out[1] += (0.5*dH_dpx[0]*hamil[3]*f[7]+0.5*hamil[1]*dH_dpx[2]*f[7]+0.5*dH_dpx[0]*hamil[2]*f[5]+0.5*hamil[0]*dH_dpx[2]*f[5]+0.5*dH_dpx[2]*hamil[3]*f[4]+0.5*dH_dpx[0]*hamil[1]*f[4]+0.5*f[1]*dH_dpx[2]*hamil[2]+0.5*dH_dpx[0]*hamil[0]*f[1])*volFact; 
  out[2] += (0.5*dH_dpy[0]*hamil[3]*f[6]+0.5*dH_dpy[1]*hamil[2]*f[6]+0.5*dH_dpy[1]*f[3]*hamil[3]+0.5*dH_dpy[0]*hamil[2]*f[3]+0.5*dH_dpy[0]*hamil[1]*f[2]+0.5*hamil[0]*dH_dpy[1]*f[2]+0.5*f[0]*dH_dpy[1]*hamil[1]+0.5*dH_dpy[0]*f[0]*hamil[0])*volFact; 
  out[3] += (0.5*dH_dpy[0]*hamil[3]*f[7]+0.5*dH_dpy[1]*hamil[2]*f[7]+0.5*dH_dpy[1]*hamil[3]*f[5]+0.5*dH_dpy[0]*hamil[2]*f[5]+0.5*dH_dpy[0]*hamil[1]*f[4]+0.5*hamil[0]*dH_dpy[1]*f[4]+0.5*dH_dpy[1]*f[1]*hamil[1]+0.5*dH_dpy[0]*hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_five_moments_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += (1.7320508075688772*f[3]*hamil[3]+1.7320508075688772*f[0]*hamil[1])*dv10*jacob_vx_inv*volFact; 
  out[3] += (1.7320508075688772*hamil[3]*f[5]+1.7320508075688772*f[1]*hamil[1])*dv10*jacob_vx_inv*volFact; 
  out[4] += (1.7320508075688772*f[2]*hamil[3]+1.7320508075688772*f[0]*hamil[2])*dv11*jacob_vy_inv*volFact; 
  out[5] += (1.7320508075688772*hamil[3]*f[4]+1.7320508075688772*f[1]*hamil[2])*dv11*jacob_vy_inv*volFact; 
  out[6] += (hamil[3]*f[6]+hamil[2]*f[3]+hamil[1]*f[2]+f[0]*hamil[0])*volFact; 
  out[7] += (hamil[3]*f[7]+hamil[2]*f[5]+hamil[1]*f[4]+hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_int_five_moments_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[0]*dxv[1]*dxv[2]*0.125; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += (2.4494897427831783*f[3]*hamil[3]+2.4494897427831783*f[0]*hamil[1])*dv10*jacob_vx_inv*volFact; 
  out[2] += (2.4494897427831783*f[2]*hamil[3]+2.4494897427831783*f[0]*hamil[2])*dv11*jacob_vy_inv*volFact; 
  out[3] += (1.4142135623730951*hamil[3]*f[6]+1.4142135623730951*hamil[2]*f[3]+1.4142135623730951*hamil[1]*f[2]+1.4142135623730951*f[0]*hamil[0])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_sparse_M1i_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  out[0] += 1.7320508075688772*f[0]*hamil[1]*dv10*jacob_vx_inv*volFact; 
  out[1] += 1.7320508075688772*f[1]*hamil[1]*dv10*jacob_vx_inv*volFact; 
  out[2] += 1.7320508075688772*f[0]*hamil[2]*dv11*jacob_vy_inv*volFact; 
  out[3] += 1.7320508075688772*f[1]*hamil[2]*dv11*jacob_vy_inv*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_sparse_M2_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  out[0] += (hamil[2]*f[3]+hamil[1]*f[2]+f[0]*hamil[0])*volFact; 
  out[1] += (hamil[2]*f[5]+hamil[1]*f[4]+hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_sparse_M3i_1x2v_ser_p1(const double *w, const double *dxv, const int *idx,
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  double dH_dpx[4] = {0.0}; 
  dH_dpx[0] = 1.7320508075688772*hamil[1]*dv10*jacob_vx_inv; 

  double dH_dpy[4] = {0.0}; 
  dH_dpy[0] = 1.7320508075688772*hamil[2]*dv11*jacob_vy_inv; 

  out[0] += (0.5*dH_dpx[0]*hamil[2]*f[3]+0.5*dH_dpx[0]*hamil[1]*f[2]+0.5*dH_dpx[0]*f[0]*hamil[0])*volFact; 
  out[1] += (0.5*dH_dpx[0]*hamil[2]*f[5]+0.5*dH_dpx[0]*hamil[1]*f[4]+0.5*dH_dpx[0]*hamil[0]*f[1])*volFact; 
  out[2] += (0.5*dH_dpy[0]*hamil[2]*f[3]+0.5*dH_dpy[0]*hamil[1]*f[2]+0.5*dH_dpy[0]*f[0]*hamil[0])*volFact; 
  out[3] += (0.5*dH_dpy[0]*hamil[2]*f[5]+0.5*dH_dpy[0]*hamil[1]*f[4]+0.5*dH_dpy[0]*hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_sparse_five_moments_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 1.7320508075688772*f[0]*hamil[1]*dv10*jacob_vx_inv*volFact; 
  out[3] += 1.7320508075688772*f[1]*hamil[1]*dv10*jacob_vx_inv*volFact; 
  out[4] += 1.7320508075688772*f[0]*hamil[2]*dv11*jacob_vy_inv*volFact; 
  out[5] += 1.7320508075688772*f[1]*hamil[2]*dv11*jacob_vy_inv*volFact; 
  out[6] += (hamil[2]*f[3]+hamil[1]*f[2]+f[0]*hamil[0])*volFact; 
  out[7] += (hamil[2]*f[5]+hamil[1]*f[4]+hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_sparse_int_five_moments_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[0]*dxv[1]*dxv[2]*0.125; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += 2.4494897427831783*f[0]*hamil[1]*dv10*jacob_vx_inv*volFact; 
  out[2] += 2.4494897427831783*f[0]*hamil[2]*dv11*jacob_vy_inv*volFact; 
  out[3] += (1.4142135623730951*hamil[2]*f[3]+1.4142135623730951*hamil[1]*f[2]+1.4142135623730951*f[0]*hamil[0])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_phase_M1i_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  out[0] += (1.224744871391589*f[5]*hamil[7]+1.224744871391589*f[3]*hamil[6]+1.224744871391589*f[1]*hamil[4]+1.224744871391589*f[0]*hamil[2])*dv10*jacob_vx_inv*volFact; 
  out[1] += (1.224744871391589*f[3]*hamil[7]+1.224744871391589*f[5]*hamil[6]+1.224744871391589*f[0]*hamil[4]+1.224744871391589*f[1]*hamil[2])*dv10*jacob_vx_inv*volFact; 
  out[2] += (1.224744871391589*f[4]*hamil[7]+1.224744871391589*f[2]*hamil[6]+1.224744871391589*f[1]*hamil[5]+1.224744871391589*f[0]*hamil[3])*dv11*jacob_vy_inv*volFact; 
  out[3] += (1.224744871391589*f[2]*hamil[7]+1.224744871391589*f[4]*hamil[6]+1.224744871391589*f[0]*hamil[5]+1.224744871391589*f[1]*hamil[3])*dv11*jacob_vy_inv*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_phase_M2_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  out[0] += (0.7071067811865475*f[7]*hamil[7]+0.7071067811865475*f[6]*hamil[6]+0.7071067811865475*f[5]*hamil[5]+0.7071067811865475*f[4]*hamil[4]+0.7071067811865475*f[3]*hamil[3]+0.7071067811865475*f[2]*hamil[2]+0.7071067811865475*f[1]*hamil[1]+0.7071067811865475*f[0]*hamil[0])*volFact; 
  out[1] += (0.7071067811865475*f[6]*hamil[7]+0.7071067811865475*hamil[6]*f[7]+0.7071067811865475*f[3]*hamil[5]+0.7071067811865475*hamil[3]*f[5]+0.7071067811865475*f[2]*hamil[4]+0.7071067811865475*hamil[2]*f[4]+0.7071067811865475*f[0]*hamil[1]+0.7071067811865475*hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_phase_five_moments_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += (1.224744871391589*f[5]*hamil[7]+1.224744871391589*f[3]*hamil[6]+1.224744871391589*f[1]*hamil[4]+1.224744871391589*f[0]*hamil[2])*dv10*jacob_vx_inv*volFact; 
  out[3] += (1.224744871391589*f[3]*hamil[7]+1.224744871391589*f[5]*hamil[6]+1.224744871391589*f[0]*hamil[4]+1.224744871391589*f[1]*hamil[2])*dv10*jacob_vx_inv*volFact; 
  out[4] += (1.224744871391589*f[4]*hamil[7]+1.224744871391589*f[2]*hamil[6]+1.224744871391589*f[1]*hamil[5]+1.224744871391589*f[0]*hamil[3])*dv11*jacob_vy_inv*volFact; 
  out[5] += (1.224744871391589*f[2]*hamil[7]+1.224744871391589*f[4]*hamil[6]+1.224744871391589*f[0]*hamil[5]+1.224744871391589*f[1]*hamil[3])*dv11*jacob_vy_inv*volFact; 
  out[6] += (0.7071067811865475*f[7]*hamil[7]+0.7071067811865475*f[6]*hamil[6]+0.7071067811865475*f[5]*hamil[5]+0.7071067811865475*f[4]*hamil[4]+0.7071067811865475*f[3]*hamil[3]+0.7071067811865475*f[2]*hamil[2]+0.7071067811865475*f[1]*hamil[1]+0.7071067811865475*f[0]*hamil[0])*volFact; 
  out[7] += (0.7071067811865475*f[6]*hamil[7]+0.7071067811865475*hamil[6]*f[7]+0.7071067811865475*f[3]*hamil[5]+0.7071067811865475*hamil[3]*f[5]+0.7071067811865475*f[2]*hamil[4]+0.7071067811865475*hamil[2]*f[4]+0.7071067811865475*f[0]*hamil[1]+0.7071067811865475*hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_phase_int_five_moments_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[0]*dxv[1]*dxv[2]*0.125; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += (1.7320508075688772*f[5]*hamil[7]+1.7320508075688772*f[3]*hamil[6]+1.7320508075688772*f[1]*hamil[4]+1.7320508075688772*f[0]*hamil[2])*dv10*jacob_vx_inv*volFact; 
  out[2] += (1.7320508075688772*f[4]*hamil[7]+1.7320508075688772*f[2]*hamil[6]+1.7320508075688772*f[1]*hamil[5]+1.7320508075688772*f[0]*hamil[3])*dv11*jacob_vy_inv*volFact; 
  out[3] += (f[7]*hamil[7]+f[6]*hamil[6]+f[5]*hamil[5]+f[4]*hamil[4]+f[3]*hamil[3]+f[2]*hamil[2]+f[1]*hamil[1]+f[0]*hamil[0])*volFact; 
} 
