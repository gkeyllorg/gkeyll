#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void mom_vlasov_M0_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, 
  const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[7]*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M2ij_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, 
  const double *vmap, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  const double *vmap_vx = &vmap[0]; 
  const double *vmap_vy = &vmap[4]; 
  const double vmap_vx0_sq = vmap_vx[0]*vmap_vx[0]; 
  const double vmap_vx1_sq = vmap_vx[1]*vmap_vx[1]; 
  const double vmap_vy0_sq = vmap_vy[0]*vmap_vy[0]; 
  const double vmap_vy1_sq = vmap_vy[1]*vmap_vy[1]; 
  out[0] += (0.8944271909999159*f[8]*vmap_vx1_sq+f[0]*vmap_vx1_sq+f[0]*vmap_vx0_sq+2.0*vmap_vx[0]*vmap_vx[1]*f[2])*volFact; 
  out[1] += (0.8944271909999161*f[12]*vmap_vx1_sq+f[1]*vmap_vx1_sq+f[1]*vmap_vx0_sq+2.0*vmap_vx[0]*vmap_vx[1]*f[4])*volFact; 
  out[2] += (f[7]*vmap_vx1_sq+f[7]*vmap_vx0_sq+2.0000000000000004*vmap_vx[0]*vmap_vx[1]*f[11])*volFact; 
  out[3] += (vmap_vx[1]*vmap_vy[1]*f[6]+vmap_vx[0]*vmap_vy[1]*f[3]+vmap_vy[0]*vmap_vx[1]*f[2]+f[0]*vmap_vx[0]*vmap_vy[0])*volFact; 
  out[4] += (vmap_vx[1]*vmap_vy[1]*f[10]+vmap_vx[0]*vmap_vy[1]*f[5]+vmap_vy[0]*vmap_vx[1]*f[4]+vmap_vx[0]*vmap_vy[0]*f[1])*volFact; 
  out[5] += (vmap_vx[1]*vmap_vy[1]*f[17]+1.0000000000000002*vmap_vx[0]*vmap_vy[1]*f[13]+1.0000000000000002*vmap_vy[0]*vmap_vx[1]*f[11]+vmap_vx[0]*vmap_vy[0]*f[7])*volFact; 
  out[6] += (0.8944271909999159*f[9]*vmap_vy1_sq+f[0]*vmap_vy1_sq+f[0]*vmap_vy0_sq+2.0*vmap_vy[0]*vmap_vy[1]*f[3])*volFact; 
  out[7] += (0.8944271909999161*f[15]*vmap_vy1_sq+f[1]*vmap_vy1_sq+f[1]*vmap_vy0_sq+2.0*vmap_vy[0]*vmap_vy[1]*f[5])*volFact; 
  out[8] += (f[7]*vmap_vy1_sq+f[7]*vmap_vy0_sq+2.0000000000000004*vmap_vy[0]*vmap_vy[1]*f[13])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M3ijk_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, 
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
  out[0] += (1.8973665961010278*vmap_vx[0]*f[8]*vmap_vx1_sq+2.1213203435596424*f[0]*vmap_vx[0]*vmap_vx1_sq+1.2727922061357855*f[2]*vmap_vx1_cu+2.1213203435596424*vmap_vx[1]*f[2]*vmap_vx0_sq+0.7071067811865475*f[0]*vmap_vx0_cu)*volFact; 
  out[1] += (1.897366596101028*vmap_vx[0]*f[12]*vmap_vx1_sq+2.1213203435596424*vmap_vx[0]*f[1]*vmap_vx1_sq+1.2727922061357855*f[4]*vmap_vx1_cu+2.1213203435596424*vmap_vx[1]*f[4]*vmap_vx0_sq+0.7071067811865475*f[1]*vmap_vx0_cu)*volFact; 
  out[2] += (2.1213203435596424*vmap_vx[0]*f[7]*vmap_vx1_sq+1.2727922061357853*f[11]*vmap_vx1_cu+2.1213203435596424*vmap_vx[1]*f[11]*vmap_vx0_sq+0.7071067811865475*f[7]*vmap_vx0_cu)*volFact; 
  out[3] += (0.632455532033676*vmap_vy[1]*f[14]*vmap_vx1_sq+0.6324555320336759*vmap_vy[0]*f[8]*vmap_vx1_sq+0.7071067811865475*vmap_vy[1]*f[3]*vmap_vx1_sq+0.7071067811865475*f[0]*vmap_vy[0]*vmap_vx1_sq+0.7071067811865475*vmap_vy[1]*f[3]*vmap_vx0_sq+0.7071067811865475*f[0]*vmap_vy[0]*vmap_vx0_sq+1.4142135623730951*vmap_vx[0]*vmap_vx[1]*vmap_vy[1]*f[6]+1.4142135623730951*vmap_vx[0]*vmap_vy[0]*vmap_vx[1]*f[2])*volFact; 
  out[4] += (0.6324555320336759*vmap_vy[1]*f[18]*vmap_vx1_sq+0.632455532033676*vmap_vy[0]*f[12]*vmap_vx1_sq+0.7071067811865475*vmap_vy[1]*f[5]*vmap_vx1_sq+0.7071067811865475*vmap_vy[0]*f[1]*vmap_vx1_sq+0.7071067811865475*vmap_vy[1]*f[5]*vmap_vx0_sq+0.7071067811865475*vmap_vy[0]*f[1]*vmap_vx0_sq+1.4142135623730951*vmap_vx[0]*vmap_vx[1]*vmap_vy[1]*f[10]+1.4142135623730951*vmap_vx[0]*vmap_vy[0]*vmap_vx[1]*f[4])*volFact; 
  out[5] += (0.7071067811865475*vmap_vy[1]*f[13]*vmap_vx1_sq+0.7071067811865475*vmap_vy[0]*f[7]*vmap_vx1_sq+0.7071067811865475*vmap_vy[1]*f[13]*vmap_vx0_sq+0.7071067811865475*vmap_vy[0]*f[7]*vmap_vx0_sq+1.4142135623730951*vmap_vx[0]*vmap_vx[1]*vmap_vy[1]*f[17]+1.4142135623730951*vmap_vx[0]*vmap_vy[0]*vmap_vx[1]*f[11])*volFact; 
  out[6] += (0.632455532033676*vmap_vx[1]*f[16]*vmap_vy1_sq+0.6324555320336759*vmap_vx[0]*f[9]*vmap_vy1_sq+0.7071067811865475*vmap_vx[1]*f[2]*vmap_vy1_sq+0.7071067811865475*f[0]*vmap_vx[0]*vmap_vy1_sq+0.7071067811865475*vmap_vx[1]*f[2]*vmap_vy0_sq+0.7071067811865475*f[0]*vmap_vx[0]*vmap_vy0_sq+1.4142135623730951*vmap_vy[0]*vmap_vx[1]*vmap_vy[1]*f[6]+1.4142135623730951*vmap_vx[0]*vmap_vy[0]*vmap_vy[1]*f[3])*volFact; 
  out[7] += (0.6324555320336759*vmap_vx[1]*f[19]*vmap_vy1_sq+0.632455532033676*vmap_vx[0]*f[15]*vmap_vy1_sq+0.7071067811865475*vmap_vx[1]*f[4]*vmap_vy1_sq+0.7071067811865475*vmap_vx[0]*f[1]*vmap_vy1_sq+0.7071067811865475*vmap_vx[1]*f[4]*vmap_vy0_sq+0.7071067811865475*vmap_vx[0]*f[1]*vmap_vy0_sq+1.4142135623730951*vmap_vy[0]*vmap_vx[1]*vmap_vy[1]*f[10]+1.4142135623730951*vmap_vx[0]*vmap_vy[0]*vmap_vy[1]*f[5])*volFact; 
  out[8] += (0.7071067811865475*vmap_vx[1]*f[11]*vmap_vy1_sq+0.7071067811865475*vmap_vx[0]*f[7]*vmap_vy1_sq+0.7071067811865475*vmap_vx[1]*f[11]*vmap_vy0_sq+0.7071067811865475*vmap_vx[0]*f[7]*vmap_vy0_sq+1.4142135623730951*vmap_vy[0]*vmap_vx[1]*vmap_vy[1]*f[17]+1.4142135623730951*vmap_vx[0]*vmap_vy[0]*vmap_vy[1]*f[13])*volFact; 
  out[9] += (1.8973665961010278*vmap_vy[0]*f[9]*vmap_vy1_sq+2.1213203435596424*f[0]*vmap_vy[0]*vmap_vy1_sq+1.2727922061357855*f[3]*vmap_vy1_cu+2.1213203435596424*vmap_vy[1]*f[3]*vmap_vy0_sq+0.7071067811865475*f[0]*vmap_vy0_cu)*volFact; 
  out[10] += (1.897366596101028*vmap_vy[0]*f[15]*vmap_vy1_sq+2.1213203435596424*vmap_vy[0]*f[1]*vmap_vy1_sq+1.2727922061357855*f[5]*vmap_vy1_cu+2.1213203435596424*vmap_vy[1]*f[5]*vmap_vy0_sq+0.7071067811865475*f[1]*vmap_vy0_cu)*volFact; 
  out[11] += (2.1213203435596424*vmap_vy[0]*f[7]*vmap_vy1_sq+1.2727922061357853*f[13]*vmap_vy1_cu+2.1213203435596424*vmap_vy[1]*f[13]*vmap_vy0_sq+0.7071067811865475*f[7]*vmap_vy0_cu)*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_M1i_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[3]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  out[0] += (1.7320508075688774*hamil[7]*f[9]+3.872983346207417*f[6]*hamil[6]+3.872983346207417*f[2]*hamil[4]+1.7320508075688772*f[3]*hamil[3]+1.7320508075688772*f[0]*hamil[1])*dv10*jacob_vx_inv*volFact; 
  out[1] += (1.7320508075688772*hamil[7]*f[15]+3.872983346207417*hamil[6]*f[10]+1.7320508075688772*hamil[3]*f[5]+3.872983346207417*f[4]*hamil[4]+1.7320508075688772*f[1]*hamil[1])*dv10*jacob_vx_inv*volFact; 
  out[2] += (3.872983346207417*hamil[6]*f[17]+1.7320508075688774*hamil[3]*f[13]+3.872983346207417*hamil[4]*f[11]+1.7320508075688772*hamil[1]*f[7])*dv10*jacob_vx_inv*volFact; 
  out[3] += (1.7320508075688774*hamil[6]*f[8]+3.872983346207417*f[6]*hamil[7]+3.872983346207417*f[3]*hamil[5]+1.7320508075688772*f[2]*hamil[3]+1.7320508075688772*f[0]*hamil[2])*dv11*jacob_vy_inv*volFact; 
  out[4] += (1.7320508075688772*hamil[6]*f[12]+3.872983346207417*hamil[7]*f[10]+3.872983346207417*f[5]*hamil[5]+1.7320508075688772*hamil[3]*f[4]+1.7320508075688772*f[1]*hamil[2])*dv11*jacob_vy_inv*volFact; 
  out[5] += (3.872983346207417*hamil[7]*f[17]+3.872983346207417*hamil[5]*f[13]+1.7320508075688774*hamil[3]*f[11]+1.7320508075688772*hamil[2]*f[7])*dv11*jacob_vy_inv*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_M2_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  out[0] += (hamil[7]*f[16]+hamil[6]*f[14]+hamil[5]*f[9]+hamil[4]*f[8]+hamil[3]*f[6]+hamil[2]*f[3]+hamil[1]*f[2]+f[0]*hamil[0])*volFact; 
  out[1] += (1.0000000000000002*hamil[7]*f[19]+1.0000000000000002*hamil[6]*f[18]+1.0000000000000002*hamil[5]*f[15]+1.0000000000000002*hamil[4]*f[12]+hamil[3]*f[10]+hamil[2]*f[5]+hamil[1]*f[4]+hamil[0]*f[1])*volFact; 
  out[2] += (hamil[3]*f[17]+1.0000000000000002*hamil[2]*f[13]+1.0000000000000002*hamil[1]*f[11]+hamil[0]*f[7])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_M3i_1x2v_ser_p2(const double *w, const double *dxv, const int *idx,
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[3]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  double dH_dpx[8] = {0.0}; 
  dH_dpx[0] = 1.7320508075688772*hamil[1]*dv10*jacob_vx_inv; 
  dH_dpx[1] = 3.872983346207417*hamil[4]*dv10*jacob_vx_inv; 
  dH_dpx[2] = 1.7320508075688772*hamil[3]*dv10*jacob_vx_inv; 
  dH_dpx[3] = 3.872983346207417*hamil[6]*dv10*jacob_vx_inv; 
  dH_dpx[5] = 1.7320508075688774*hamil[7]*dv10*jacob_vx_inv; 

  double dH_dpy[8] = {0.0}; 
  dH_dpy[0] = 1.7320508075688772*hamil[2]*dv11*jacob_vy_inv; 
  dH_dpy[1] = 1.7320508075688772*hamil[3]*dv11*jacob_vy_inv; 
  dH_dpy[2] = 3.872983346207417*hamil[5]*dv11*jacob_vy_inv; 
  dH_dpy[3] = 3.872983346207417*hamil[7]*dv11*jacob_vy_inv; 
  dH_dpy[4] = 1.7320508075688772*hamil[6]*dv11*jacob_vy_inv; 

  out[0] += (0.31943828249996997*dH_dpx[5]*hamil[7]*f[16]+0.5*dH_dpx[0]*hamil[7]*f[16]+0.4*dH_dpx[3]*hamil[6]*f[16]+0.5000000000000001*dH_dpx[1]*hamil[5]*f[16]+0.5000000000000001*hamil[1]*dH_dpx[5]*f[16]+0.44721359549995804*dH_dpx[2]*hamil[3]*f[16]+0.44721359549995804*hamil[2]*dH_dpx[3]*f[16]+0.4*dH_dpx[3]*hamil[7]*f[14]+0.4472135954999579*dH_dpx[5]*hamil[6]*f[14]+0.5*dH_dpx[0]*hamil[6]*f[14]+0.5000000000000001*dH_dpx[2]*hamil[4]*f[14]+0.44721359549995804*dH_dpx[1]*hamil[3]*f[14]+0.44721359549995804*hamil[1]*dH_dpx[3]*f[14]+0.5000000000000001*dH_dpx[1]*hamil[7]*f[9]+0.31943828249996997*dH_dpx[5]*hamil[5]*f[9]+0.5*dH_dpx[0]*hamil[5]*f[9]+0.5*hamil[0]*dH_dpx[5]*f[9]+0.4472135954999579*dH_dpx[3]*hamil[3]*f[9]+0.4472135954999579*dH_dpx[2]*hamil[2]*f[9]+0.5000000000000001*dH_dpx[2]*hamil[6]*f[8]+0.5*dH_dpx[0]*hamil[4]*f[8]+0.4472135954999579*dH_dpx[3]*hamil[3]*f[8]+0.4472135954999579*dH_dpx[1]*hamil[1]*f[8]+0.44721359549995804*dH_dpx[2]*f[6]*hamil[7]+0.5000000000000001*f[2]*dH_dpx[5]*hamil[7]+0.44721359549995804*dH_dpx[3]*f[3]*hamil[7]+0.44721359549995804*dH_dpx[1]*f[6]*hamil[6]+0.44721359549995804*f[2]*dH_dpx[3]*hamil[6]+0.4472135954999579*dH_dpx[3]*hamil[5]*f[6]+0.4472135954999579*hamil[3]*dH_dpx[5]*f[6]+0.4472135954999579*dH_dpx[3]*hamil[4]*f[6]+0.5*dH_dpx[0]*hamil[3]*f[6]+0.5*hamil[0]*dH_dpx[3]*f[6]+0.5*dH_dpx[1]*hamil[2]*f[6]+0.5*hamil[1]*dH_dpx[2]*f[6]+0.5*f[0]*dH_dpx[5]*hamil[5]+0.4472135954999579*dH_dpx[2]*f[3]*hamil[5]+0.4472135954999579*hamil[2]*f[3]*dH_dpx[5]+0.4472135954999579*dH_dpx[1]*f[2]*hamil[4]+0.5*dH_dpx[1]*f[3]*hamil[3]+0.5*f[0]*dH_dpx[3]*hamil[3]+0.5*dH_dpx[2]*f[2]*hamil[3]+0.5*hamil[1]*dH_dpx[3]*f[3]+0.5*dH_dpx[0]*hamil[2]*f[3]+0.5*hamil[0]*dH_dpx[2]*f[3]+0.5*f[2]*hamil[2]*dH_dpx[3]+0.5*f[0]*dH_dpx[2]*hamil[2]+0.5*dH_dpx[0]*hamil[1]*f[2]+0.5*hamil[0]*dH_dpx[1]*f[2]+0.5*f[0]*dH_dpx[1]*hamil[1]+0.5*dH_dpx[0]*f[0]*hamil[0])*volFact; 
  out[1] += (0.31943828249997*dH_dpx[5]*hamil[7]*f[19]+0.5000000000000001*dH_dpx[0]*hamil[7]*f[19]+0.4*dH_dpx[3]*hamil[6]*f[19]+0.5*dH_dpx[1]*hamil[5]*f[19]+0.5*hamil[1]*dH_dpx[5]*f[19]+0.4472135954999579*dH_dpx[2]*hamil[3]*f[19]+0.4472135954999579*hamil[2]*dH_dpx[3]*f[19]+0.4*dH_dpx[3]*hamil[7]*f[18]+0.44721359549995804*dH_dpx[5]*hamil[6]*f[18]+0.5000000000000001*dH_dpx[0]*hamil[6]*f[18]+0.5*dH_dpx[2]*hamil[4]*f[18]+0.4472135954999579*dH_dpx[1]*hamil[3]*f[18]+0.4472135954999579*hamil[1]*dH_dpx[3]*f[18]+0.5*dH_dpx[1]*hamil[7]*f[15]+0.31943828249997*dH_dpx[5]*hamil[5]*f[15]+0.5000000000000001*dH_dpx[0]*hamil[5]*f[15]+0.5000000000000001*hamil[0]*dH_dpx[5]*f[15]+0.44721359549995804*dH_dpx[3]*hamil[3]*f[15]+0.44721359549995804*dH_dpx[2]*hamil[2]*f[15]+0.5*dH_dpx[2]*hamil[6]*f[12]+0.5000000000000001*dH_dpx[0]*hamil[4]*f[12]+0.44721359549995804*dH_dpx[3]*hamil[3]*f[12]+0.44721359549995804*dH_dpx[1]*hamil[1]*f[12]+0.44721359549995804*dH_dpx[2]*hamil[7]*f[10]+0.44721359549995804*dH_dpx[1]*hamil[6]*f[10]+0.4472135954999579*dH_dpx[3]*hamil[5]*f[10]+0.4472135954999579*hamil[3]*dH_dpx[5]*f[10]+0.4472135954999579*dH_dpx[3]*hamil[4]*f[10]+0.5*dH_dpx[0]*hamil[3]*f[10]+0.5*hamil[0]*dH_dpx[3]*f[10]+0.5*dH_dpx[1]*hamil[2]*f[10]+0.5*hamil[1]*dH_dpx[2]*f[10]+0.44721359549995804*dH_dpx[3]*f[5]*hamil[7]+0.5000000000000001*f[4]*dH_dpx[5]*hamil[7]+0.44721359549995804*dH_dpx[3]*f[4]*hamil[6]+0.4472135954999579*dH_dpx[2]*f[5]*hamil[5]+0.5*f[1]*dH_dpx[5]*hamil[5]+0.4472135954999579*hamil[2]*dH_dpx[5]*f[5]+0.5*dH_dpx[1]*hamil[3]*f[5]+0.5*hamil[1]*dH_dpx[3]*f[5]+0.5*dH_dpx[0]*hamil[2]*f[5]+0.5*hamil[0]*dH_dpx[2]*f[5]+0.4472135954999579*dH_dpx[1]*f[4]*hamil[4]+0.5*dH_dpx[2]*hamil[3]*f[4]+0.5*hamil[2]*dH_dpx[3]*f[4]+0.5*dH_dpx[0]*hamil[1]*f[4]+0.5*hamil[0]*dH_dpx[1]*f[4]+0.5*f[1]*dH_dpx[3]*hamil[3]+0.5*f[1]*dH_dpx[2]*hamil[2]+0.5*dH_dpx[1]*f[1]*hamil[1]+0.5*dH_dpx[0]*hamil[0]*f[1])*volFact; 
  out[2] += (0.44721359549995804*dH_dpx[2]*hamil[7]*f[17]+0.44721359549995804*dH_dpx[1]*hamil[6]*f[17]+0.4472135954999579*dH_dpx[3]*hamil[5]*f[17]+0.4472135954999579*hamil[3]*dH_dpx[5]*f[17]+0.4472135954999579*dH_dpx[3]*hamil[4]*f[17]+0.5*dH_dpx[0]*hamil[3]*f[17]+0.5*hamil[0]*dH_dpx[3]*f[17]+0.5*dH_dpx[1]*hamil[2]*f[17]+0.5*hamil[1]*dH_dpx[2]*f[17]+0.4472135954999579*dH_dpx[3]*hamil[7]*f[13]+0.44721359549995804*dH_dpx[2]*hamil[5]*f[13]+0.44721359549995804*hamil[2]*dH_dpx[5]*f[13]+0.5000000000000001*dH_dpx[1]*hamil[3]*f[13]+0.5000000000000001*hamil[1]*dH_dpx[3]*f[13]+0.5000000000000001*dH_dpx[0]*hamil[2]*f[13]+0.5000000000000001*hamil[0]*dH_dpx[2]*f[13]+0.5*dH_dpx[5]*hamil[7]*f[11]+0.4472135954999579*dH_dpx[3]*hamil[6]*f[11]+0.44721359549995804*dH_dpx[1]*hamil[4]*f[11]+0.5000000000000001*dH_dpx[2]*hamil[3]*f[11]+0.5000000000000001*hamil[2]*dH_dpx[3]*f[11]+0.5000000000000001*dH_dpx[0]*hamil[1]*f[11]+0.5000000000000001*hamil[0]*dH_dpx[1]*f[11]+0.5*dH_dpx[5]*hamil[5]*f[7]+0.5*dH_dpx[3]*hamil[3]*f[7]+0.5*dH_dpx[2]*hamil[2]*f[7]+0.5*dH_dpx[1]*hamil[1]*f[7]+0.5*dH_dpx[0]*hamil[0]*f[7])*volFact; 
  out[3] += (0.4472135954999579*dH_dpy[4]*hamil[7]*f[16]+0.5*dH_dpy[0]*hamil[7]*f[16]+0.4*dH_dpy[3]*hamil[6]*f[16]+0.5000000000000001*dH_dpy[1]*hamil[5]*f[16]+0.44721359549995804*dH_dpy[2]*hamil[3]*f[16]+0.44721359549995804*hamil[2]*dH_dpy[3]*f[16]+0.4*dH_dpy[3]*hamil[7]*f[14]+0.31943828249996997*dH_dpy[4]*hamil[6]*f[14]+0.5*dH_dpy[0]*hamil[6]*f[14]+0.5000000000000001*dH_dpy[2]*hamil[4]*f[14]+0.5000000000000001*hamil[2]*dH_dpy[4]*f[14]+0.44721359549995804*dH_dpy[1]*hamil[3]*f[14]+0.44721359549995804*hamil[1]*dH_dpy[3]*f[14]+0.5000000000000001*dH_dpy[1]*hamil[7]*f[9]+0.5*dH_dpy[0]*hamil[5]*f[9]+0.4472135954999579*dH_dpy[3]*hamil[3]*f[9]+0.4472135954999579*dH_dpy[2]*hamil[2]*f[9]+0.5000000000000001*dH_dpy[2]*hamil[6]*f[8]+0.31943828249996997*dH_dpy[4]*hamil[4]*f[8]+0.5*dH_dpy[0]*hamil[4]*f[8]+0.5*hamil[0]*dH_dpy[4]*f[8]+0.4472135954999579*dH_dpy[3]*hamil[3]*f[8]+0.4472135954999579*dH_dpy[1]*hamil[1]*f[8]+0.44721359549995804*dH_dpy[2]*f[6]*hamil[7]+0.44721359549995804*dH_dpy[3]*f[3]*hamil[7]+0.44721359549995804*dH_dpy[1]*f[6]*hamil[6]+0.5000000000000001*f[3]*dH_dpy[4]*hamil[6]+0.44721359549995804*f[2]*dH_dpy[3]*hamil[6]+0.4472135954999579*dH_dpy[3]*hamil[5]*f[6]+0.4472135954999579*dH_dpy[3]*hamil[4]*f[6]+0.4472135954999579*hamil[3]*dH_dpy[4]*f[6]+0.5*dH_dpy[0]*hamil[3]*f[6]+0.5*hamil[0]*dH_dpy[3]*f[6]+0.5*dH_dpy[1]*hamil[2]*f[6]+0.5*hamil[1]*dH_dpy[2]*f[6]+0.4472135954999579*dH_dpy[2]*f[3]*hamil[5]+0.5*f[0]*dH_dpy[4]*hamil[4]+0.4472135954999579*dH_dpy[1]*f[2]*hamil[4]+0.4472135954999579*hamil[1]*f[2]*dH_dpy[4]+0.5*dH_dpy[1]*f[3]*hamil[3]+0.5*f[0]*dH_dpy[3]*hamil[3]+0.5*dH_dpy[2]*f[2]*hamil[3]+0.5*hamil[1]*dH_dpy[3]*f[3]+0.5*dH_dpy[0]*hamil[2]*f[3]+0.5*hamil[0]*dH_dpy[2]*f[3]+0.5*f[2]*hamil[2]*dH_dpy[3]+0.5*f[0]*dH_dpy[2]*hamil[2]+0.5*dH_dpy[0]*hamil[1]*f[2]+0.5*hamil[0]*dH_dpy[1]*f[2]+0.5*f[0]*dH_dpy[1]*hamil[1]+0.5*dH_dpy[0]*f[0]*hamil[0])*volFact; 
  out[4] += (0.44721359549995804*dH_dpy[4]*hamil[7]*f[19]+0.5000000000000001*dH_dpy[0]*hamil[7]*f[19]+0.4*dH_dpy[3]*hamil[6]*f[19]+0.5*dH_dpy[1]*hamil[5]*f[19]+0.4472135954999579*dH_dpy[2]*hamil[3]*f[19]+0.4472135954999579*hamil[2]*dH_dpy[3]*f[19]+0.4*dH_dpy[3]*hamil[7]*f[18]+0.31943828249997*dH_dpy[4]*hamil[6]*f[18]+0.5000000000000001*dH_dpy[0]*hamil[6]*f[18]+0.5*dH_dpy[2]*hamil[4]*f[18]+0.5*hamil[2]*dH_dpy[4]*f[18]+0.4472135954999579*dH_dpy[1]*hamil[3]*f[18]+0.4472135954999579*hamil[1]*dH_dpy[3]*f[18]+0.5*dH_dpy[1]*hamil[7]*f[15]+0.5000000000000001*dH_dpy[0]*hamil[5]*f[15]+0.44721359549995804*dH_dpy[3]*hamil[3]*f[15]+0.44721359549995804*dH_dpy[2]*hamil[2]*f[15]+0.5*dH_dpy[2]*hamil[6]*f[12]+0.31943828249997*dH_dpy[4]*hamil[4]*f[12]+0.5000000000000001*dH_dpy[0]*hamil[4]*f[12]+0.5000000000000001*hamil[0]*dH_dpy[4]*f[12]+0.44721359549995804*dH_dpy[3]*hamil[3]*f[12]+0.44721359549995804*dH_dpy[1]*hamil[1]*f[12]+0.44721359549995804*dH_dpy[2]*hamil[7]*f[10]+0.44721359549995804*dH_dpy[1]*hamil[6]*f[10]+0.4472135954999579*dH_dpy[3]*hamil[5]*f[10]+0.4472135954999579*dH_dpy[3]*hamil[4]*f[10]+0.4472135954999579*hamil[3]*dH_dpy[4]*f[10]+0.5*dH_dpy[0]*hamil[3]*f[10]+0.5*hamil[0]*dH_dpy[3]*f[10]+0.5*dH_dpy[1]*hamil[2]*f[10]+0.5*hamil[1]*dH_dpy[2]*f[10]+0.44721359549995804*dH_dpy[3]*f[5]*hamil[7]+0.5000000000000001*dH_dpy[4]*f[5]*hamil[6]+0.44721359549995804*dH_dpy[3]*f[4]*hamil[6]+0.4472135954999579*dH_dpy[2]*f[5]*hamil[5]+0.5*dH_dpy[1]*hamil[3]*f[5]+0.5*hamil[1]*dH_dpy[3]*f[5]+0.5*dH_dpy[0]*hamil[2]*f[5]+0.5*hamil[0]*dH_dpy[2]*f[5]+0.4472135954999579*dH_dpy[1]*f[4]*hamil[4]+0.5*f[1]*dH_dpy[4]*hamil[4]+0.4472135954999579*hamil[1]*dH_dpy[4]*f[4]+0.5*dH_dpy[2]*hamil[3]*f[4]+0.5*hamil[2]*dH_dpy[3]*f[4]+0.5*dH_dpy[0]*hamil[1]*f[4]+0.5*hamil[0]*dH_dpy[1]*f[4]+0.5*f[1]*dH_dpy[3]*hamil[3]+0.5*f[1]*dH_dpy[2]*hamil[2]+0.5*dH_dpy[1]*f[1]*hamil[1]+0.5*dH_dpy[0]*hamil[0]*f[1])*volFact; 
  out[5] += (0.44721359549995804*dH_dpy[2]*hamil[7]*f[17]+0.44721359549995804*dH_dpy[1]*hamil[6]*f[17]+0.4472135954999579*dH_dpy[3]*hamil[5]*f[17]+0.4472135954999579*dH_dpy[3]*hamil[4]*f[17]+0.4472135954999579*hamil[3]*dH_dpy[4]*f[17]+0.5*dH_dpy[0]*hamil[3]*f[17]+0.5*hamil[0]*dH_dpy[3]*f[17]+0.5*dH_dpy[1]*hamil[2]*f[17]+0.5*hamil[1]*dH_dpy[2]*f[17]+0.4472135954999579*dH_dpy[3]*hamil[7]*f[13]+0.5*dH_dpy[4]*hamil[6]*f[13]+0.44721359549995804*dH_dpy[2]*hamil[5]*f[13]+0.5000000000000001*dH_dpy[1]*hamil[3]*f[13]+0.5000000000000001*hamil[1]*dH_dpy[3]*f[13]+0.5000000000000001*dH_dpy[0]*hamil[2]*f[13]+0.5000000000000001*hamil[0]*dH_dpy[2]*f[13]+0.4472135954999579*dH_dpy[3]*hamil[6]*f[11]+0.44721359549995804*dH_dpy[1]*hamil[4]*f[11]+0.44721359549995804*hamil[1]*dH_dpy[4]*f[11]+0.5000000000000001*dH_dpy[2]*hamil[3]*f[11]+0.5000000000000001*hamil[2]*dH_dpy[3]*f[11]+0.5000000000000001*dH_dpy[0]*hamil[1]*f[11]+0.5000000000000001*hamil[0]*dH_dpy[1]*f[11]+0.5*dH_dpy[4]*hamil[4]*f[7]+0.5*dH_dpy[3]*hamil[3]*f[7]+0.5*dH_dpy[2]*hamil[2]*f[7]+0.5*dH_dpy[1]*hamil[1]*f[7]+0.5*dH_dpy[0]*hamil[0]*f[7])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_five_moments_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[3]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[7]*volFact; 
  out[3] += (1.7320508075688774*hamil[7]*f[9]+3.872983346207417*f[6]*hamil[6]+3.872983346207417*f[2]*hamil[4]+1.7320508075688772*f[3]*hamil[3]+1.7320508075688772*f[0]*hamil[1])*dv10*jacob_vx_inv*volFact; 
  out[4] += (1.7320508075688772*hamil[7]*f[15]+3.872983346207417*hamil[6]*f[10]+1.7320508075688772*hamil[3]*f[5]+3.872983346207417*f[4]*hamil[4]+1.7320508075688772*f[1]*hamil[1])*dv10*jacob_vx_inv*volFact; 
  out[5] += (3.872983346207417*hamil[6]*f[17]+1.7320508075688774*hamil[3]*f[13]+3.872983346207417*hamil[4]*f[11]+1.7320508075688772*hamil[1]*f[7])*dv10*jacob_vx_inv*volFact; 
  out[6] += (1.7320508075688774*hamil[6]*f[8]+3.872983346207417*f[6]*hamil[7]+3.872983346207417*f[3]*hamil[5]+1.7320508075688772*f[2]*hamil[3]+1.7320508075688772*f[0]*hamil[2])*dv11*jacob_vy_inv*volFact; 
  out[7] += (1.7320508075688772*hamil[6]*f[12]+3.872983346207417*hamil[7]*f[10]+3.872983346207417*f[5]*hamil[5]+1.7320508075688772*hamil[3]*f[4]+1.7320508075688772*f[1]*hamil[2])*dv11*jacob_vy_inv*volFact; 
  out[8] += (3.872983346207417*hamil[7]*f[17]+3.872983346207417*hamil[5]*f[13]+1.7320508075688774*hamil[3]*f[11]+1.7320508075688772*hamil[2]*f[7])*dv11*jacob_vy_inv*volFact; 
  out[9] += (hamil[7]*f[16]+hamil[6]*f[14]+hamil[5]*f[9]+hamil[4]*f[8]+hamil[3]*f[6]+hamil[2]*f[3]+hamil[1]*f[2]+f[0]*hamil[0])*volFact; 
  out[10] += (1.0000000000000002*hamil[7]*f[19]+1.0000000000000002*hamil[6]*f[18]+1.0000000000000002*hamil[5]*f[15]+1.0000000000000002*hamil[4]*f[12]+hamil[3]*f[10]+hamil[2]*f[5]+hamil[1]*f[4]+hamil[0]*f[1])*volFact; 
  out[11] += (hamil[3]*f[17]+1.0000000000000002*hamil[2]*f[13]+1.0000000000000002*hamil[1]*f[11]+hamil[0]*f[7])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_int_five_moments_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[0]*dxv[1]*dxv[2]*0.125; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[3]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += (2.4494897427831783*hamil[7]*f[9]+5.477225575051662*f[6]*hamil[6]+5.477225575051662*f[2]*hamil[4]+2.4494897427831783*f[3]*hamil[3]+2.4494897427831783*f[0]*hamil[1])*dv10*jacob_vx_inv*volFact; 
  out[2] += (2.4494897427831783*hamil[6]*f[8]+5.477225575051662*f[6]*hamil[7]+5.477225575051662*f[3]*hamil[5]+2.4494897427831783*f[2]*hamil[3]+2.4494897427831783*f[0]*hamil[2])*dv11*jacob_vy_inv*volFact; 
  out[3] += (1.4142135623730951*hamil[7]*f[16]+1.4142135623730951*hamil[6]*f[14]+1.4142135623730951*hamil[5]*f[9]+1.4142135623730951*hamil[4]*f[8]+1.4142135623730951*hamil[3]*f[6]+1.4142135623730951*hamil[2]*f[3]+1.4142135623730951*hamil[1]*f[2]+1.4142135623730951*f[0]*hamil[0])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_sparse_M1i_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[3]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  out[0] += (3.872983346207417*f[2]*hamil[4]+1.7320508075688772*f[0]*hamil[1])*dv10*jacob_vx_inv*volFact; 
  out[1] += (3.872983346207417*f[4]*hamil[4]+1.7320508075688772*f[1]*hamil[1])*dv10*jacob_vx_inv*volFact; 
  out[2] += (3.872983346207417*hamil[4]*f[11]+1.7320508075688772*hamil[1]*f[7])*dv10*jacob_vx_inv*volFact; 
  out[3] += (3.872983346207417*f[3]*hamil[5]+1.7320508075688772*f[0]*hamil[2])*dv11*jacob_vy_inv*volFact; 
  out[4] += (3.872983346207417*f[5]*hamil[5]+1.7320508075688772*f[1]*hamil[2])*dv11*jacob_vy_inv*volFact; 
  out[5] += (3.872983346207417*hamil[5]*f[13]+1.7320508075688772*hamil[2]*f[7])*dv11*jacob_vy_inv*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_sparse_M2_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  out[0] += (hamil[5]*f[9]+hamil[4]*f[8]+hamil[2]*f[3]+hamil[1]*f[2]+f[0]*hamil[0])*volFact; 
  out[1] += (1.0000000000000002*hamil[5]*f[15]+1.0000000000000002*hamil[4]*f[12]+hamil[2]*f[5]+hamil[1]*f[4]+hamil[0]*f[1])*volFact; 
  out[2] += (1.0000000000000002*hamil[2]*f[13]+1.0000000000000002*hamil[1]*f[11]+hamil[0]*f[7])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_sparse_M3i_1x2v_ser_p2(const double *w, const double *dxv, const int *idx,
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[3]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  double dH_dpx[8] = {0.0}; 
  dH_dpx[0] = 1.7320508075688772*hamil[1]*dv10*jacob_vx_inv; 
  dH_dpx[1] = 3.872983346207417*hamil[4]*dv10*jacob_vx_inv; 

  double dH_dpy[8] = {0.0}; 
  dH_dpy[0] = 1.7320508075688772*hamil[2]*dv11*jacob_vy_inv; 
  dH_dpy[2] = 3.872983346207417*hamil[5]*dv11*jacob_vy_inv; 

  out[0] += (0.5000000000000001*dH_dpx[1]*hamil[5]*f[16]+0.5*dH_dpx[0]*hamil[5]*f[9]+0.5*dH_dpx[0]*hamil[4]*f[8]+0.4472135954999579*dH_dpx[1]*hamil[1]*f[8]+0.5*dH_dpx[1]*hamil[2]*f[6]+0.4472135954999579*dH_dpx[1]*f[2]*hamil[4]+0.5*dH_dpx[0]*hamil[2]*f[3]+0.5*dH_dpx[0]*hamil[1]*f[2]+0.5*hamil[0]*dH_dpx[1]*f[2]+0.5*f[0]*dH_dpx[1]*hamil[1]+0.5*dH_dpx[0]*f[0]*hamil[0])*volFact; 
  out[1] += (0.5*dH_dpx[1]*hamil[5]*f[19]+0.5000000000000001*dH_dpx[0]*hamil[5]*f[15]+0.5000000000000001*dH_dpx[0]*hamil[4]*f[12]+0.44721359549995804*dH_dpx[1]*hamil[1]*f[12]+0.5*dH_dpx[1]*hamil[2]*f[10]+0.5*dH_dpx[0]*hamil[2]*f[5]+0.4472135954999579*dH_dpx[1]*f[4]*hamil[4]+0.5*dH_dpx[0]*hamil[1]*f[4]+0.5*hamil[0]*dH_dpx[1]*f[4]+0.5*dH_dpx[1]*f[1]*hamil[1]+0.5*dH_dpx[0]*hamil[0]*f[1])*volFact; 
  out[2] += (0.5*dH_dpx[1]*hamil[2]*f[17]+0.5000000000000001*dH_dpx[0]*hamil[2]*f[13]+0.44721359549995804*dH_dpx[1]*hamil[4]*f[11]+0.5000000000000001*dH_dpx[0]*hamil[1]*f[11]+0.5000000000000001*hamil[0]*dH_dpx[1]*f[11]+0.5*dH_dpx[1]*hamil[1]*f[7]+0.5*dH_dpx[0]*hamil[0]*f[7])*volFact; 
  out[3] += (0.5000000000000001*dH_dpy[2]*hamil[4]*f[14]+0.5*dH_dpy[0]*hamil[5]*f[9]+0.4472135954999579*dH_dpy[2]*hamil[2]*f[9]+0.5*dH_dpy[0]*hamil[4]*f[8]+0.5*hamil[1]*dH_dpy[2]*f[6]+0.4472135954999579*dH_dpy[2]*f[3]*hamil[5]+0.5*dH_dpy[0]*hamil[2]*f[3]+0.5*hamil[0]*dH_dpy[2]*f[3]+0.5*f[0]*dH_dpy[2]*hamil[2]+0.5*dH_dpy[0]*hamil[1]*f[2]+0.5*dH_dpy[0]*f[0]*hamil[0])*volFact; 
  out[4] += (0.5*dH_dpy[2]*hamil[4]*f[18]+0.5000000000000001*dH_dpy[0]*hamil[5]*f[15]+0.44721359549995804*dH_dpy[2]*hamil[2]*f[15]+0.5000000000000001*dH_dpy[0]*hamil[4]*f[12]+0.5*hamil[1]*dH_dpy[2]*f[10]+0.4472135954999579*dH_dpy[2]*f[5]*hamil[5]+0.5*dH_dpy[0]*hamil[2]*f[5]+0.5*hamil[0]*dH_dpy[2]*f[5]+0.5*dH_dpy[0]*hamil[1]*f[4]+0.5*f[1]*dH_dpy[2]*hamil[2]+0.5*dH_dpy[0]*hamil[0]*f[1])*volFact; 
  out[5] += (0.5*hamil[1]*dH_dpy[2]*f[17]+0.44721359549995804*dH_dpy[2]*hamil[5]*f[13]+0.5000000000000001*dH_dpy[0]*hamil[2]*f[13]+0.5000000000000001*hamil[0]*dH_dpy[2]*f[13]+0.5000000000000001*dH_dpy[0]*hamil[1]*f[11]+0.5*dH_dpy[2]*hamil[2]*f[7]+0.5*dH_dpy[0]*hamil[0]*f[7])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_sparse_five_moments_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[3]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[7]*volFact; 
  out[3] += (3.872983346207417*f[2]*hamil[4]+1.7320508075688772*f[0]*hamil[1])*dv10*jacob_vx_inv*volFact; 
  out[4] += (3.872983346207417*f[4]*hamil[4]+1.7320508075688772*f[1]*hamil[1])*dv10*jacob_vx_inv*volFact; 
  out[5] += (3.872983346207417*hamil[4]*f[11]+1.7320508075688772*hamil[1]*f[7])*dv10*jacob_vx_inv*volFact; 
  out[6] += (3.872983346207417*f[3]*hamil[5]+1.7320508075688772*f[0]*hamil[2])*dv11*jacob_vy_inv*volFact; 
  out[7] += (3.872983346207417*f[5]*hamil[5]+1.7320508075688772*f[1]*hamil[2])*dv11*jacob_vy_inv*volFact; 
  out[8] += (3.872983346207417*hamil[5]*f[13]+1.7320508075688772*hamil[2]*f[7])*dv11*jacob_vy_inv*volFact; 
  out[9] += (hamil[5]*f[9]+hamil[4]*f[8]+hamil[2]*f[3]+hamil[1]*f[2]+f[0]*hamil[0])*volFact; 
  out[10] += (1.0000000000000002*hamil[5]*f[15]+1.0000000000000002*hamil[4]*f[12]+hamil[2]*f[5]+hamil[1]*f[4]+hamil[0]*f[1])*volFact; 
  out[11] += (1.0000000000000002*hamil[2]*f[13]+1.0000000000000002*hamil[1]*f[11]+hamil[0]*f[7])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_sparse_int_five_moments_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[0]*dxv[1]*dxv[2]*0.125; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[3]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += (5.477225575051662*f[2]*hamil[4]+2.4494897427831783*f[0]*hamil[1])*dv10*jacob_vx_inv*volFact; 
  out[2] += (5.477225575051662*f[3]*hamil[5]+2.4494897427831783*f[0]*hamil[2])*dv11*jacob_vy_inv*volFact; 
  out[3] += (1.4142135623730951*hamil[5]*f[9]+1.4142135623730951*hamil[4]*f[8]+1.4142135623730951*hamil[2]*f[3]+1.4142135623730951*hamil[1]*f[2]+1.4142135623730951*f[0]*hamil[0])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_phase_M1i_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[3]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  out[0] += (1.224744871391589*f[15]*hamil[19]+2.7386127875258306*f[10]*hamil[18]+1.224744871391589*f[13]*hamil[17]+1.224744871391589*f[9]*hamil[16]+2.7386127875258306*f[6]*hamil[14]+2.7386127875258306*f[4]*hamil[12]+1.224744871391589*f[7]*hamil[11]+1.224744871391589*f[5]*hamil[10]+2.7386127875258306*f[2]*hamil[8]+1.224744871391589*f[3]*hamil[6]+1.224744871391589*f[1]*hamil[4]+1.224744871391589*f[0]*hamil[2])*dv10*jacob_vx_inv*volFact; 
  out[1] += (1.224744871391589*f[9]*hamil[19]+2.4494897427831783*f[17]*hamil[18]+2.7386127875258306*f[6]*hamil[18]+1.0954451150103324*f[5]*hamil[17]+1.224744871391589*f[15]*hamil[16]+2.7386127875258306*f[10]*hamil[14]+1.0954451150103324*hamil[10]*f[13]+2.4494897427831783*f[11]*hamil[12]+2.7386127875258306*f[2]*hamil[12]+1.0954451150103324*f[1]*hamil[11]+1.224744871391589*f[3]*hamil[10]+2.7386127875258306*f[4]*hamil[8]+1.0954451150103324*hamil[4]*f[7]+1.224744871391589*f[5]*hamil[6]+1.224744871391589*f[0]*hamil[4]+1.224744871391589*f[1]*hamil[2])*dv10*jacob_vx_inv*volFact; 
  out[2] += (1.0954451150103324*f[15]*hamil[19]+2.4494897427831783*f[10]*hamil[18]+0.7824607964359517*f[13]*hamil[17]+1.224744871391589*f[3]*hamil[17]+2.7386127875258306*hamil[14]*f[17]+1.224744871391589*hamil[6]*f[13]+2.4494897427831783*f[4]*hamil[12]+0.7824607964359517*f[7]*hamil[11]+1.224744871391589*f[0]*hamil[11]+2.7386127875258306*hamil[8]*f[11]+1.0954451150103324*f[5]*hamil[10]+1.224744871391589*hamil[2]*f[7]+1.0954451150103324*f[1]*hamil[4])*dv10*jacob_vx_inv*volFact; 
  out[3] += (2.7386127875258306*f[10]*hamil[19]+1.224744871391589*f[12]*hamil[18]+1.224744871391589*f[11]*hamil[17]+2.7386127875258306*f[6]*hamil[16]+2.7386127875258306*f[5]*hamil[15]+1.224744871391589*f[8]*hamil[14]+1.224744871391589*f[7]*hamil[13]+1.224744871391589*f[4]*hamil[10]+2.7386127875258306*f[3]*hamil[9]+1.224744871391589*f[2]*hamil[6]+1.224744871391589*f[1]*hamil[5]+1.224744871391589*f[0]*hamil[3])*dv11*jacob_vy_inv*volFact; 
  out[4] += (2.4494897427831783*f[17]*hamil[19]+2.7386127875258306*f[6]*hamil[19]+1.224744871391589*f[8]*hamil[18]+1.0954451150103324*f[4]*hamil[17]+2.7386127875258306*f[10]*hamil[16]+2.4494897427831783*f[13]*hamil[15]+2.7386127875258306*f[3]*hamil[15]+1.224744871391589*f[12]*hamil[14]+1.0954451150103324*f[1]*hamil[13]+1.0954451150103324*hamil[10]*f[11]+1.224744871391589*f[2]*hamil[10]+2.7386127875258306*f[5]*hamil[9]+1.0954451150103324*hamil[5]*f[7]+1.224744871391589*f[4]*hamil[6]+1.224744871391589*f[0]*hamil[5]+1.224744871391589*f[1]*hamil[3])*dv11*jacob_vy_inv*volFact; 
  out[5] += (2.4494897427831783*f[10]*hamil[19]+1.0954451150103324*f[12]*hamil[18]+0.7824607964359517*f[11]*hamil[17]+1.224744871391589*f[2]*hamil[17]+2.7386127875258306*hamil[16]*f[17]+2.4494897427831783*f[5]*hamil[15]+0.7824607964359517*f[7]*hamil[13]+1.224744871391589*f[0]*hamil[13]+2.7386127875258306*hamil[9]*f[13]+1.224744871391589*hamil[6]*f[11]+1.0954451150103324*f[4]*hamil[10]+1.224744871391589*hamil[3]*f[7]+1.0954451150103324*f[1]*hamil[5])*dv11*jacob_vy_inv*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_phase_M2_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  out[0] += (0.7071067811865475*f[19]*hamil[19]+0.7071067811865475*f[18]*hamil[18]+0.7071067811865475*f[17]*hamil[17]+0.7071067811865475*f[16]*hamil[16]+0.7071067811865475*f[15]*hamil[15]+0.7071067811865475*f[14]*hamil[14]+0.7071067811865475*f[13]*hamil[13]+0.7071067811865475*f[12]*hamil[12]+0.7071067811865475*f[11]*hamil[11]+0.7071067811865475*f[10]*hamil[10]+0.7071067811865475*f[9]*hamil[9]+0.7071067811865475*f[8]*hamil[8]+0.7071067811865475*f[7]*hamil[7]+0.7071067811865475*f[6]*hamil[6]+0.7071067811865475*f[5]*hamil[5]+0.7071067811865475*f[4]*hamil[4]+0.7071067811865475*f[3]*hamil[3]+0.7071067811865475*f[2]*hamil[2]+0.7071067811865475*f[1]*hamil[1]+0.7071067811865475*f[0]*hamil[0])*volFact; 
  out[1] += (0.7071067811865475*f[16]*hamil[19]+0.7071067811865475*hamil[16]*f[19]+0.7071067811865475*f[14]*hamil[18]+0.7071067811865475*hamil[14]*f[18]+0.6324555320336759*f[10]*hamil[17]+0.6324555320336759*hamil[10]*f[17]+0.7071067811865475*f[9]*hamil[15]+0.7071067811865475*hamil[9]*f[15]+0.632455532033676*f[5]*hamil[13]+0.632455532033676*hamil[5]*f[13]+0.7071067811865475*f[8]*hamil[12]+0.7071067811865475*hamil[8]*f[12]+0.632455532033676*f[4]*hamil[11]+0.632455532033676*hamil[4]*f[11]+0.7071067811865475*f[6]*hamil[10]+0.7071067811865475*hamil[6]*f[10]+0.6324555320336759*f[1]*hamil[7]+0.6324555320336759*hamil[1]*f[7]+0.7071067811865475*f[3]*hamil[5]+0.7071067811865475*hamil[3]*f[5]+0.7071067811865475*f[2]*hamil[4]+0.7071067811865475*hamil[2]*f[4]+0.7071067811865475*f[0]*hamil[1]+0.7071067811865475*hamil[0]*f[1])*volFact; 
  out[2] += (0.6324555320336759*f[19]*hamil[19]+0.6324555320336759*f[18]*hamil[18]+0.45175395145262565*f[17]*hamil[17]+0.7071067811865475*f[6]*hamil[17]+0.7071067811865475*hamil[6]*f[17]+0.6324555320336759*f[15]*hamil[15]+0.45175395145262565*f[13]*hamil[13]+0.7071067811865475*f[3]*hamil[13]+0.7071067811865475*hamil[3]*f[13]+0.6324555320336759*f[12]*hamil[12]+0.45175395145262565*f[11]*hamil[11]+0.7071067811865475*f[2]*hamil[11]+0.7071067811865475*hamil[2]*f[11]+0.6324555320336759*f[10]*hamil[10]+0.45175395145262565*f[7]*hamil[7]+0.7071067811865475*f[0]*hamil[7]+0.7071067811865475*hamil[0]*f[7]+0.6324555320336759*f[5]*hamil[5]+0.6324555320336759*f[4]*hamil[4]+0.6324555320336759*f[1]*hamil[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_phase_five_moments_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[3]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[7]*volFact; 
  out[3] += (1.224744871391589*f[15]*hamil[19]+2.7386127875258306*f[10]*hamil[18]+1.224744871391589*f[13]*hamil[17]+1.224744871391589*f[9]*hamil[16]+2.7386127875258306*f[6]*hamil[14]+2.7386127875258306*f[4]*hamil[12]+1.224744871391589*f[7]*hamil[11]+1.224744871391589*f[5]*hamil[10]+2.7386127875258306*f[2]*hamil[8]+1.224744871391589*f[3]*hamil[6]+1.224744871391589*f[1]*hamil[4]+1.224744871391589*f[0]*hamil[2])*dv10*jacob_vx_inv*volFact; 
  out[4] += (1.224744871391589*f[9]*hamil[19]+2.4494897427831783*f[17]*hamil[18]+2.7386127875258306*f[6]*hamil[18]+1.0954451150103324*f[5]*hamil[17]+1.224744871391589*f[15]*hamil[16]+2.7386127875258306*f[10]*hamil[14]+1.0954451150103324*hamil[10]*f[13]+2.4494897427831783*f[11]*hamil[12]+2.7386127875258306*f[2]*hamil[12]+1.0954451150103324*f[1]*hamil[11]+1.224744871391589*f[3]*hamil[10]+2.7386127875258306*f[4]*hamil[8]+1.0954451150103324*hamil[4]*f[7]+1.224744871391589*f[5]*hamil[6]+1.224744871391589*f[0]*hamil[4]+1.224744871391589*f[1]*hamil[2])*dv10*jacob_vx_inv*volFact; 
  out[5] += (1.0954451150103324*f[15]*hamil[19]+2.4494897427831783*f[10]*hamil[18]+0.7824607964359517*f[13]*hamil[17]+1.224744871391589*f[3]*hamil[17]+2.7386127875258306*hamil[14]*f[17]+1.224744871391589*hamil[6]*f[13]+2.4494897427831783*f[4]*hamil[12]+0.7824607964359517*f[7]*hamil[11]+1.224744871391589*f[0]*hamil[11]+2.7386127875258306*hamil[8]*f[11]+1.0954451150103324*f[5]*hamil[10]+1.224744871391589*hamil[2]*f[7]+1.0954451150103324*f[1]*hamil[4])*dv10*jacob_vx_inv*volFact; 
  out[6] += (2.7386127875258306*f[10]*hamil[19]+1.224744871391589*f[12]*hamil[18]+1.224744871391589*f[11]*hamil[17]+2.7386127875258306*f[6]*hamil[16]+2.7386127875258306*f[5]*hamil[15]+1.224744871391589*f[8]*hamil[14]+1.224744871391589*f[7]*hamil[13]+1.224744871391589*f[4]*hamil[10]+2.7386127875258306*f[3]*hamil[9]+1.224744871391589*f[2]*hamil[6]+1.224744871391589*f[1]*hamil[5]+1.224744871391589*f[0]*hamil[3])*dv11*jacob_vy_inv*volFact; 
  out[7] += (2.4494897427831783*f[17]*hamil[19]+2.7386127875258306*f[6]*hamil[19]+1.224744871391589*f[8]*hamil[18]+1.0954451150103324*f[4]*hamil[17]+2.7386127875258306*f[10]*hamil[16]+2.4494897427831783*f[13]*hamil[15]+2.7386127875258306*f[3]*hamil[15]+1.224744871391589*f[12]*hamil[14]+1.0954451150103324*f[1]*hamil[13]+1.0954451150103324*hamil[10]*f[11]+1.224744871391589*f[2]*hamil[10]+2.7386127875258306*f[5]*hamil[9]+1.0954451150103324*hamil[5]*f[7]+1.224744871391589*f[4]*hamil[6]+1.224744871391589*f[0]*hamil[5]+1.224744871391589*f[1]*hamil[3])*dv11*jacob_vy_inv*volFact; 
  out[8] += (2.4494897427831783*f[10]*hamil[19]+1.0954451150103324*f[12]*hamil[18]+0.7824607964359517*f[11]*hamil[17]+1.224744871391589*f[2]*hamil[17]+2.7386127875258306*hamil[16]*f[17]+2.4494897427831783*f[5]*hamil[15]+0.7824607964359517*f[7]*hamil[13]+1.224744871391589*f[0]*hamil[13]+2.7386127875258306*hamil[9]*f[13]+1.224744871391589*hamil[6]*f[11]+1.0954451150103324*f[4]*hamil[10]+1.224744871391589*hamil[3]*f[7]+1.0954451150103324*f[1]*hamil[5])*dv11*jacob_vy_inv*volFact; 
  out[9] += (0.7071067811865475*f[19]*hamil[19]+0.7071067811865475*f[18]*hamil[18]+0.7071067811865475*f[17]*hamil[17]+0.7071067811865475*f[16]*hamil[16]+0.7071067811865475*f[15]*hamil[15]+0.7071067811865475*f[14]*hamil[14]+0.7071067811865475*f[13]*hamil[13]+0.7071067811865475*f[12]*hamil[12]+0.7071067811865475*f[11]*hamil[11]+0.7071067811865475*f[10]*hamil[10]+0.7071067811865475*f[9]*hamil[9]+0.7071067811865475*f[8]*hamil[8]+0.7071067811865475*f[7]*hamil[7]+0.7071067811865475*f[6]*hamil[6]+0.7071067811865475*f[5]*hamil[5]+0.7071067811865475*f[4]*hamil[4]+0.7071067811865475*f[3]*hamil[3]+0.7071067811865475*f[2]*hamil[2]+0.7071067811865475*f[1]*hamil[1]+0.7071067811865475*f[0]*hamil[0])*volFact; 
  out[10] += (0.7071067811865475*f[16]*hamil[19]+0.7071067811865475*hamil[16]*f[19]+0.7071067811865475*f[14]*hamil[18]+0.7071067811865475*hamil[14]*f[18]+0.6324555320336759*f[10]*hamil[17]+0.6324555320336759*hamil[10]*f[17]+0.7071067811865475*f[9]*hamil[15]+0.7071067811865475*hamil[9]*f[15]+0.632455532033676*f[5]*hamil[13]+0.632455532033676*hamil[5]*f[13]+0.7071067811865475*f[8]*hamil[12]+0.7071067811865475*hamil[8]*f[12]+0.632455532033676*f[4]*hamil[11]+0.632455532033676*hamil[4]*f[11]+0.7071067811865475*f[6]*hamil[10]+0.7071067811865475*hamil[6]*f[10]+0.6324555320336759*f[1]*hamil[7]+0.6324555320336759*hamil[1]*f[7]+0.7071067811865475*f[3]*hamil[5]+0.7071067811865475*hamil[3]*f[5]+0.7071067811865475*f[2]*hamil[4]+0.7071067811865475*hamil[2]*f[4]+0.7071067811865475*f[0]*hamil[1]+0.7071067811865475*hamil[0]*f[1])*volFact; 
  out[11] += (0.6324555320336759*f[19]*hamil[19]+0.6324555320336759*f[18]*hamil[18]+0.45175395145262565*f[17]*hamil[17]+0.7071067811865475*f[6]*hamil[17]+0.7071067811865475*hamil[6]*f[17]+0.6324555320336759*f[15]*hamil[15]+0.45175395145262565*f[13]*hamil[13]+0.7071067811865475*f[3]*hamil[13]+0.7071067811865475*hamil[3]*f[13]+0.6324555320336759*f[12]*hamil[12]+0.45175395145262565*f[11]*hamil[11]+0.7071067811865475*f[2]*hamil[11]+0.7071067811865475*hamil[2]*f[11]+0.6324555320336759*f[10]*hamil[10]+0.45175395145262565*f[7]*hamil[7]+0.7071067811865475*f[0]*hamil[7]+0.7071067811865475*hamil[0]*f[7]+0.6324555320336759*f[5]*hamil[5]+0.6324555320336759*f[4]*hamil[4]+0.6324555320336759*f[1]*hamil[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_phase_int_five_moments_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[0]*dxv[1]*dxv[2]*0.125; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[3]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += (1.7320508075688774*f[15]*hamil[19]+3.872983346207417*f[10]*hamil[18]+1.7320508075688774*f[13]*hamil[17]+1.7320508075688774*f[9]*hamil[16]+3.872983346207417*f[6]*hamil[14]+3.872983346207417*f[4]*hamil[12]+1.7320508075688774*f[7]*hamil[11]+1.7320508075688772*f[5]*hamil[10]+3.872983346207417*f[2]*hamil[8]+1.7320508075688772*f[3]*hamil[6]+1.7320508075688772*f[1]*hamil[4]+1.7320508075688772*f[0]*hamil[2])*dv10*jacob_vx_inv*volFact; 
  out[2] += (3.872983346207417*f[10]*hamil[19]+1.7320508075688774*f[12]*hamil[18]+1.7320508075688774*f[11]*hamil[17]+3.872983346207417*f[6]*hamil[16]+3.872983346207417*f[5]*hamil[15]+1.7320508075688774*f[8]*hamil[14]+1.7320508075688774*f[7]*hamil[13]+1.7320508075688772*f[4]*hamil[10]+3.872983346207417*f[3]*hamil[9]+1.7320508075688772*f[2]*hamil[6]+1.7320508075688772*f[1]*hamil[5]+1.7320508075688772*f[0]*hamil[3])*dv11*jacob_vy_inv*volFact; 
  out[3] += (f[19]*hamil[19]+f[18]*hamil[18]+f[17]*hamil[17]+f[16]*hamil[16]+f[15]*hamil[15]+f[14]*hamil[14]+f[13]*hamil[13]+f[12]*hamil[12]+f[11]*hamil[11]+f[10]*hamil[10]+f[9]*hamil[9]+f[8]*hamil[8]+f[7]*hamil[7]+f[6]*hamil[6]+f[5]*hamil[5]+f[4]*hamil[4]+f[3]*hamil[3]+f[2]*hamil[2]+f[1]*hamil[1]+f[0]*hamil[0])*volFact; 
} 
