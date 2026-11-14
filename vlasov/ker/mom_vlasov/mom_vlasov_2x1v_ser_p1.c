#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void mom_vlasov_M0_2x1v_ser_p1(const double *w, const double *dxv, const int *idx, 
  const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]/2; 
  out[0] += 1.4142135623730951*f[0]*volFact; 
  out[1] += 1.4142135623730951*f[1]*volFact; 
  out[2] += 1.4142135623730951*f[2]*volFact; 
  out[3] += 1.4142135623730951*f[4]*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M2ij_2x1v_ser_p1(const double *w, const double *dxv, const int *idx, 
  const double *vmap, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]/2; 
  const double *vmap_vx = &vmap[0]; 
  const double vmap_vx0_sq = vmap_vx[0]*vmap_vx[0]; 
  const double vmap_vx1_sq = vmap_vx[1]*vmap_vx[1]; 
  out[0] += (0.7071067811865475*f[0]*vmap_vx1_sq+0.7071067811865475*f[0]*vmap_vx0_sq+1.4142135623730951*vmap_vx[0]*vmap_vx[1]*f[3])*volFact; 
  out[1] += (0.7071067811865475*f[1]*vmap_vx1_sq+0.7071067811865475*f[1]*vmap_vx0_sq+1.4142135623730951*vmap_vx[0]*vmap_vx[1]*f[5])*volFact; 
  out[2] += (0.7071067811865475*f[2]*vmap_vx1_sq+0.7071067811865475*f[2]*vmap_vx0_sq+1.4142135623730951*vmap_vx[0]*vmap_vx[1]*f[6])*volFact; 
  out[3] += (0.7071067811865475*f[4]*vmap_vx1_sq+0.7071067811865475*f[4]*vmap_vx0_sq+1.4142135623730951*vmap_vx[0]*vmap_vx[1]*f[7])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M3ijk_2x1v_ser_p1(const double *w, const double *dxv, const int *idx, 
  const double *vmap, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]/2; 
  const double *vmap_vx = &vmap[0]; 
  const double vmap_vx0_sq = vmap_vx[0]*vmap_vx[0]; 
  const double vmap_vx0_cu = vmap_vx[0]*vmap_vx[0]*vmap_vx[0]; 
  const double vmap_vx1_sq = vmap_vx[1]*vmap_vx[1]; 
  const double vmap_vx1_cu = vmap_vx[1]*vmap_vx[1]*vmap_vx[1]; 
  out[0] += (1.5*f[0]*vmap_vx[0]*vmap_vx1_sq+0.9*f[3]*vmap_vx1_cu+1.5*vmap_vx[1]*f[3]*vmap_vx0_sq+0.5*f[0]*vmap_vx0_cu)*volFact; 
  out[1] += (1.5*vmap_vx[0]*f[1]*vmap_vx1_sq+0.9*f[5]*vmap_vx1_cu+1.5*vmap_vx[1]*f[5]*vmap_vx0_sq+0.5*f[1]*vmap_vx0_cu)*volFact; 
  out[2] += (1.5*vmap_vx[0]*f[2]*vmap_vx1_sq+0.9*f[6]*vmap_vx1_cu+1.5*vmap_vx[1]*f[6]*vmap_vx0_sq+0.5*f[2]*vmap_vx0_cu)*volFact; 
  out[3] += (1.5*vmap_vx[0]*f[4]*vmap_vx1_sq+0.9*f[7]*vmap_vx1_cu+1.5*vmap_vx[1]*f[7]*vmap_vx0_sq+0.5*f[4]*vmap_vx0_cu)*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_M1i_2x1v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]/2; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  out[0] += 1.7320508075688772*f[0]*hamil[1]*dv10*jacob_vx_inv*volFact; 
  out[1] += 1.7320508075688772*f[1]*hamil[1]*dv10*jacob_vx_inv*volFact; 
  out[2] += 1.7320508075688772*hamil[1]*f[2]*dv10*jacob_vx_inv*volFact; 
  out[3] += 1.7320508075688772*hamil[1]*f[4]*dv10*jacob_vx_inv*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_M2_2x1v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]/2; 
  out[0] += (hamil[1]*f[3]+f[0]*hamil[0])*volFact; 
  out[1] += (hamil[1]*f[5]+hamil[0]*f[1])*volFact; 
  out[2] += (hamil[1]*f[6]+hamil[0]*f[2])*volFact; 
  out[3] += (hamil[1]*f[7]+hamil[0]*f[4])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_M3i_2x1v_ser_p1(const double *w, const double *dxv, const int *idx,
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]/2; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dH_dpx[2] = {0.0}; 
  dH_dpx[0] = 1.7320508075688772*hamil[1]*dv10*jacob_vx_inv; 

  out[0] += (0.7071067811865475*dH_dpx[0]*hamil[1]*f[3]+0.7071067811865475*dH_dpx[0]*f[0]*hamil[0])*volFact; 
  out[1] += (0.7071067811865475*dH_dpx[0]*hamil[1]*f[5]+0.7071067811865475*dH_dpx[0]*hamil[0]*f[1])*volFact; 
  out[2] += (0.7071067811865475*dH_dpx[0]*hamil[1]*f[6]+0.7071067811865475*dH_dpx[0]*hamil[0]*f[2])*volFact; 
  out[3] += (0.7071067811865475*dH_dpx[0]*hamil[1]*f[7]+0.7071067811865475*dH_dpx[0]*hamil[0]*f[4])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_five_moments_2x1v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]/2; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  out[0] += 1.4142135623730951*f[0]*volFact; 
  out[1] += 1.4142135623730951*f[1]*volFact; 
  out[2] += 1.4142135623730951*f[2]*volFact; 
  out[3] += 1.4142135623730951*f[4]*volFact; 
  out[4] += 1.7320508075688772*f[0]*hamil[1]*dv10*jacob_vx_inv*volFact; 
  out[5] += 1.7320508075688772*f[1]*hamil[1]*dv10*jacob_vx_inv*volFact; 
  out[6] += 1.7320508075688772*hamil[1]*f[2]*dv10*jacob_vx_inv*volFact; 
  out[7] += 1.7320508075688772*hamil[1]*f[4]*dv10*jacob_vx_inv*volFact; 
  out[8] += (hamil[1]*f[3]+f[0]*hamil[0])*volFact; 
  out[9] += (hamil[1]*f[5]+hamil[0]*f[1])*volFact; 
  out[10] += (hamil[1]*f[6]+hamil[0]*f[2])*volFact; 
  out[11] += (hamil[1]*f[7]+hamil[0]*f[4])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_int_five_moments_2x1v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[0]*dxv[1]*dxv[2]*0.125; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += 3.4641016151377544*f[0]*hamil[1]*dv10*jacob_vx_inv*volFact; 
  out[2] += (2.0*hamil[1]*f[3]+2.0*f[0]*hamil[0])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M0_upper_2x1v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *vmap, const double *jacob_vel, double v_thresh, double f_thresh, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]/2; 
  const double *px = &vmap[0]; 
  const double *jacob_vx = &jacob_vel[0]; 
  double f_nodes[8] = {0.0};
  double fquad = 0.0;
  if ((0.7071067811865475*px[0]-0.7071067811865475*px[1]) > v_thresh) { 
    fquad = -(0.3535533905932737*f[7])+0.3535533905932737*(f[6]+f[5]+f[4])-0.3535533905932737*(f[3]+f[2]+f[1])+0.3535533905932737*f[0]; 
    f_nodes[0] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
    fquad = 0.3535533905932737*f[7]-0.3535533905932737*f[6]+0.3535533905932737*f[5]-0.3535533905932737*(f[4]+f[3])+0.3535533905932737*f[2]-0.3535533905932737*f[1]+0.3535533905932737*f[0]; 
    f_nodes[2] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
    fquad = 0.3535533905932737*(f[7]+f[6])-0.3535533905932737*(f[5]+f[4]+f[3]+f[2])+0.3535533905932737*(f[1]+f[0]); 
    f_nodes[4] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
    fquad = -(0.3535533905932737*(f[7]+f[6]+f[5]))+0.3535533905932737*f[4]-0.3535533905932737*f[3]+0.3535533905932737*(f[2]+f[1]+f[0]); 
    f_nodes[6] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
  } 
  if ((0.7071067811865475*(px[1]+px[0])) > v_thresh) { 
    fquad = 0.3535533905932737*f[7]-0.3535533905932737*(f[6]+f[5])+0.3535533905932737*(f[4]+f[3])-0.3535533905932737*(f[2]+f[1])+0.3535533905932737*f[0]; 
    f_nodes[1] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
    fquad = -(0.3535533905932737*f[7])+0.3535533905932737*f[6]-0.3535533905932737*(f[5]+f[4])+0.3535533905932737*(f[3]+f[2])-0.3535533905932737*f[1]+0.3535533905932737*f[0]; 
    f_nodes[3] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
    fquad = -(0.3535533905932737*(f[7]+f[6]))+0.3535533905932737*f[5]-0.3535533905932737*f[4]+0.3535533905932737*f[3]-0.3535533905932737*f[2]+0.3535533905932737*(f[1]+f[0]); 
    f_nodes[5] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
    fquad = 0.3535533905932737*(f[7]+f[6]+f[5]+f[4]+f[3]+f[2]+f[1]+f[0]); 
    f_nodes[7] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
  } 
  out[0] += (0.5*f_nodes[7]+0.5*f_nodes[6]+0.5*f_nodes[5]+0.5*f_nodes[4]+0.5*f_nodes[3]+0.5*f_nodes[2]+0.5*f_nodes[1]+0.5*f_nodes[0])*volFact; 
  out[1] += (0.5*f_nodes[7]+0.5*f_nodes[6]+0.5*f_nodes[5]+0.5*f_nodes[4]-0.5*f_nodes[3]-0.5*f_nodes[2]-0.5*f_nodes[1]-0.5*f_nodes[0])*volFact; 
  out[2] += (0.5*f_nodes[7]+0.5*f_nodes[6]-0.5*f_nodes[5]-0.5*f_nodes[4]+0.5*f_nodes[3]+0.5*f_nodes[2]-0.5*f_nodes[1]-0.5*f_nodes[0])*volFact; 
  out[3] += (0.5*f_nodes[7]+0.5*f_nodes[6]-0.5*f_nodes[5]-0.5*f_nodes[4]-0.5*f_nodes[3]-0.5*f_nodes[2]+0.5*f_nodes[1]+0.5*f_nodes[0])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M0_lower_2x1v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *vmap, const double *jacob_vel, double v_thresh, double f_thresh, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]/2; 
  const double *px = &vmap[0]; 
  const double *jacob_vx = &jacob_vel[0]; 
  double f_nodes[8] = {0.0};
  double fquad = 0.0;
  if ((0.7071067811865475*px[0]-0.7071067811865475*px[1]) < -v_thresh) { 
    fquad = -(0.3535533905932737*f[7])+0.3535533905932737*(f[6]+f[5]+f[4])-0.3535533905932737*(f[3]+f[2]+f[1])+0.3535533905932737*f[0]; 
    f_nodes[0] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
    fquad = 0.3535533905932737*f[7]-0.3535533905932737*f[6]+0.3535533905932737*f[5]-0.3535533905932737*(f[4]+f[3])+0.3535533905932737*f[2]-0.3535533905932737*f[1]+0.3535533905932737*f[0]; 
    f_nodes[2] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
    fquad = 0.3535533905932737*(f[7]+f[6])-0.3535533905932737*(f[5]+f[4]+f[3]+f[2])+0.3535533905932737*(f[1]+f[0]); 
    f_nodes[4] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
    fquad = -(0.3535533905932737*(f[7]+f[6]+f[5]))+0.3535533905932737*f[4]-0.3535533905932737*f[3]+0.3535533905932737*(f[2]+f[1]+f[0]); 
    f_nodes[6] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
  } 
  if ((0.7071067811865475*(px[1]+px[0])) < -v_thresh) { 
    fquad = 0.3535533905932737*f[7]-0.3535533905932737*(f[6]+f[5])+0.3535533905932737*(f[4]+f[3])-0.3535533905932737*(f[2]+f[1])+0.3535533905932737*f[0]; 
    f_nodes[1] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
    fquad = -(0.3535533905932737*f[7])+0.3535533905932737*f[6]-0.3535533905932737*(f[5]+f[4])+0.3535533905932737*(f[3]+f[2])-0.3535533905932737*f[1]+0.3535533905932737*f[0]; 
    f_nodes[3] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
    fquad = -(0.3535533905932737*(f[7]+f[6]))+0.3535533905932737*f[5]-0.3535533905932737*f[4]+0.3535533905932737*f[3]-0.3535533905932737*f[2]+0.3535533905932737*(f[1]+f[0]); 
    f_nodes[5] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
    fquad = 0.3535533905932737*(f[7]+f[6]+f[5]+f[4]+f[3]+f[2]+f[1]+f[0]); 
    f_nodes[7] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
  } 
  out[0] += (0.5*f_nodes[7]+0.5*f_nodes[6]+0.5*f_nodes[5]+0.5*f_nodes[4]+0.5*f_nodes[3]+0.5*f_nodes[2]+0.5*f_nodes[1]+0.5*f_nodes[0])*volFact; 
  out[1] += (0.5*f_nodes[7]+0.5*f_nodes[6]+0.5*f_nodes[5]+0.5*f_nodes[4]-0.5*f_nodes[3]-0.5*f_nodes[2]-0.5*f_nodes[1]-0.5*f_nodes[0])*volFact; 
  out[2] += (0.5*f_nodes[7]+0.5*f_nodes[6]-0.5*f_nodes[5]-0.5*f_nodes[4]+0.5*f_nodes[3]+0.5*f_nodes[2]-0.5*f_nodes[1]-0.5*f_nodes[0])*volFact; 
  out[3] += (0.5*f_nodes[7]+0.5*f_nodes[6]-0.5*f_nodes[5]-0.5*f_nodes[4]-0.5*f_nodes[3]-0.5*f_nodes[2]+0.5*f_nodes[1]+0.5*f_nodes[0])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_phase_M1i_2x1v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]/2; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  out[0] += (0.8660254037844386*f[4]*hamil[7]+0.8660254037844386*f[2]*hamil[6]+0.8660254037844386*f[1]*hamil[5]+0.8660254037844386*f[0]*hamil[3])*dv10*jacob_vx_inv*volFact; 
  out[1] += (0.8660254037844386*f[2]*hamil[7]+0.8660254037844386*f[4]*hamil[6]+0.8660254037844386*f[0]*hamil[5]+0.8660254037844386*f[1]*hamil[3])*dv10*jacob_vx_inv*volFact; 
  out[2] += (0.8660254037844386*f[1]*hamil[7]+0.8660254037844386*f[0]*hamil[6]+0.8660254037844386*f[4]*hamil[5]+0.8660254037844386*f[2]*hamil[3])*dv10*jacob_vx_inv*volFact; 
  out[3] += (0.8660254037844386*f[0]*hamil[7]+0.8660254037844386*f[1]*hamil[6]+0.8660254037844386*f[2]*hamil[5]+0.8660254037844386*hamil[3]*f[4])*dv10*jacob_vx_inv*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_phase_M2_2x1v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]/2; 
  out[0] += (0.5*f[7]*hamil[7]+0.5*f[6]*hamil[6]+0.5*f[5]*hamil[5]+0.5*f[4]*hamil[4]+0.5*f[3]*hamil[3]+0.5*f[2]*hamil[2]+0.5*f[1]*hamil[1]+0.5*f[0]*hamil[0])*volFact; 
  out[1] += (0.5*f[6]*hamil[7]+0.5*hamil[6]*f[7]+0.5*f[3]*hamil[5]+0.5*hamil[3]*f[5]+0.5*f[2]*hamil[4]+0.5*hamil[2]*f[4]+0.5*f[0]*hamil[1]+0.5*hamil[0]*f[1])*volFact; 
  out[2] += (0.5*f[5]*hamil[7]+0.5*hamil[5]*f[7]+0.5*f[3]*hamil[6]+0.5*hamil[3]*f[6]+0.5*f[1]*hamil[4]+0.5*hamil[1]*f[4]+0.5*f[0]*hamil[2]+0.5*hamil[0]*f[2])*volFact; 
  out[3] += (0.5*f[3]*hamil[7]+0.5*hamil[3]*f[7]+0.5*f[5]*hamil[6]+0.5*hamil[5]*f[6]+0.5*f[0]*hamil[4]+0.5*hamil[0]*f[4]+0.5*f[1]*hamil[2]+0.5*hamil[1]*f[2])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_phase_five_moments_2x1v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]/2; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  out[0] += 1.4142135623730951*f[0]*volFact; 
  out[1] += 1.4142135623730951*f[1]*volFact; 
  out[2] += 1.4142135623730951*f[2]*volFact; 
  out[3] += 1.4142135623730951*f[4]*volFact; 
  out[4] += (0.8660254037844386*f[4]*hamil[7]+0.8660254037844386*f[2]*hamil[6]+0.8660254037844386*f[1]*hamil[5]+0.8660254037844386*f[0]*hamil[3])*dv10*jacob_vx_inv*volFact; 
  out[5] += (0.8660254037844386*f[2]*hamil[7]+0.8660254037844386*f[4]*hamil[6]+0.8660254037844386*f[0]*hamil[5]+0.8660254037844386*f[1]*hamil[3])*dv10*jacob_vx_inv*volFact; 
  out[6] += (0.8660254037844386*f[1]*hamil[7]+0.8660254037844386*f[0]*hamil[6]+0.8660254037844386*f[4]*hamil[5]+0.8660254037844386*f[2]*hamil[3])*dv10*jacob_vx_inv*volFact; 
  out[7] += (0.8660254037844386*f[0]*hamil[7]+0.8660254037844386*f[1]*hamil[6]+0.8660254037844386*f[2]*hamil[5]+0.8660254037844386*hamil[3]*f[4])*dv10*jacob_vx_inv*volFact; 
  out[8] += (0.5*f[7]*hamil[7]+0.5*f[6]*hamil[6]+0.5*f[5]*hamil[5]+0.5*f[4]*hamil[4]+0.5*f[3]*hamil[3]+0.5*f[2]*hamil[2]+0.5*f[1]*hamil[1]+0.5*f[0]*hamil[0])*volFact; 
  out[9] += (0.5*f[6]*hamil[7]+0.5*hamil[6]*f[7]+0.5*f[3]*hamil[5]+0.5*hamil[3]*f[5]+0.5*f[2]*hamil[4]+0.5*hamil[2]*f[4]+0.5*f[0]*hamil[1]+0.5*hamil[0]*f[1])*volFact; 
  out[10] += (0.5*f[5]*hamil[7]+0.5*hamil[5]*f[7]+0.5*f[3]*hamil[6]+0.5*hamil[3]*f[6]+0.5*f[1]*hamil[4]+0.5*hamil[1]*f[4]+0.5*f[0]*hamil[2]+0.5*hamil[0]*f[2])*volFact; 
  out[11] += (0.5*f[3]*hamil[7]+0.5*hamil[3]*f[7]+0.5*f[5]*hamil[6]+0.5*hamil[5]*f[6]+0.5*f[0]*hamil[4]+0.5*hamil[0]*f[4]+0.5*f[1]*hamil[2]+0.5*hamil[1]*f[2])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_phase_int_five_moments_2x1v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[0]*dxv[1]*dxv[2]*0.125; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += (1.7320508075688772*f[4]*hamil[7]+1.7320508075688772*f[2]*hamil[6]+1.7320508075688772*f[1]*hamil[5]+1.7320508075688772*f[0]*hamil[3])*dv10*jacob_vx_inv*volFact; 
  out[2] += (f[7]*hamil[7]+f[6]*hamil[6]+f[5]*hamil[5]+f[4]*hamil[4]+f[3]*hamil[3]+f[2]*hamil[2]+f[1]*hamil[1]+f[0]*hamil[0])*volFact; 
} 
