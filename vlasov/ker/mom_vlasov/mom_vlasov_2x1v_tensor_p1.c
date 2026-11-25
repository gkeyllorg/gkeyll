#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void mom_vlasov_M0_2x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
  const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]/2; 
  out[0] += 1.4142135623730951*f[0]*volFact; 
  out[1] += 1.4142135623730951*f[1]*volFact; 
  out[2] += 1.4142135623730951*f[2]*volFact; 
  out[3] += 1.4142135623730951*f[4]*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M2ij_2x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
  const double *vmap, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]/2; 
  const double *vmap_vx = &vmap[0]; 
  const double vmap_vx0_sq = vmap_vx[0]*vmap_vx[0]; 
  const double vmap_vx1_sq = vmap_vx[1]*vmap_vx[1]; 
  const double vmap_vx2_sq = vmap_vx[2]*vmap_vx[2]; 
  const double vmap_vx3_sq = vmap_vx[3]*vmap_vx[3]; 
  out[0] += (0.421637021355784*f[7]*vmap_vx3_sq+0.7071067811865475*f[0]*vmap_vx3_sq+0.45175395145262565*f[7]*vmap_vx2_sq+0.7071067811865475*f[0]*vmap_vx2_sq+0.6324555320336759*f[7]*vmap_vx1_sq+0.7071067811865475*f[0]*vmap_vx1_sq+0.7071067811865475*f[0]*vmap_vx0_sq+1.2421180068162374*vmap_vx[1]*vmap_vx[3]*f[7]+1.4142135623730951*vmap_vx[0]*vmap_vx[2]*f[7]+1.2421180068162374*vmap_vx[2]*f[3]*vmap_vx[3]+1.264911064067352*vmap_vx[1]*vmap_vx[2]*f[3]+1.4142135623730951*vmap_vx[0]*vmap_vx[1]*f[3])*volFact; 
  out[1] += (0.4216370213557841*f[9]*vmap_vx3_sq+0.7071067811865475*f[1]*vmap_vx3_sq+0.4517539514526257*f[9]*vmap_vx2_sq+0.7071067811865475*f[1]*vmap_vx2_sq+0.632455532033676*f[9]*vmap_vx1_sq+0.7071067811865475*f[1]*vmap_vx1_sq+0.7071067811865475*f[1]*vmap_vx0_sq+1.2421180068162376*vmap_vx[1]*vmap_vx[3]*f[9]+1.4142135623730951*vmap_vx[0]*vmap_vx[2]*f[9]+1.2421180068162374*vmap_vx[2]*vmap_vx[3]*f[5]+1.264911064067352*vmap_vx[1]*vmap_vx[2]*f[5]+1.4142135623730951*vmap_vx[0]*vmap_vx[1]*f[5])*volFact; 
  out[2] += (0.4216370213557841*f[10]*vmap_vx3_sq+0.7071067811865475*f[2]*vmap_vx3_sq+0.4517539514526257*f[10]*vmap_vx2_sq+0.7071067811865475*f[2]*vmap_vx2_sq+0.632455532033676*f[10]*vmap_vx1_sq+0.7071067811865475*f[2]*vmap_vx1_sq+0.7071067811865475*f[2]*vmap_vx0_sq+1.2421180068162376*vmap_vx[1]*vmap_vx[3]*f[10]+1.4142135623730951*vmap_vx[0]*vmap_vx[2]*f[10]+1.2421180068162374*vmap_vx[2]*vmap_vx[3]*f[6]+1.264911064067352*vmap_vx[1]*vmap_vx[2]*f[6]+1.4142135623730951*vmap_vx[0]*vmap_vx[1]*f[6])*volFact; 
  out[3] += (0.421637021355784*f[11]*vmap_vx3_sq+0.7071067811865475*f[4]*vmap_vx3_sq+0.45175395145262565*f[11]*vmap_vx2_sq+0.7071067811865475*f[4]*vmap_vx2_sq+0.6324555320336759*f[11]*vmap_vx1_sq+0.7071067811865475*f[4]*vmap_vx1_sq+0.7071067811865475*f[4]*vmap_vx0_sq+1.2421180068162374*vmap_vx[1]*vmap_vx[3]*f[11]+1.4142135623730951*vmap_vx[0]*vmap_vx[2]*f[11]+1.2421180068162374*vmap_vx[2]*vmap_vx[3]*f[8]+1.264911064067352*vmap_vx[1]*vmap_vx[2]*f[8]+1.4142135623730951*vmap_vx[0]*vmap_vx[1]*f[8])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M3ijk_2x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
  const double *vmap, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]/2; 
  const double *vmap_vx = &vmap[0]; 
  const double vmap_vx0_sq = vmap_vx[0]*vmap_vx[0]; 
  const double vmap_vx0_cu = vmap_vx[0]*vmap_vx[0]*vmap_vx[0]; 
  const double vmap_vx1_sq = vmap_vx[1]*vmap_vx[1]; 
  const double vmap_vx1_cu = vmap_vx[1]*vmap_vx[1]*vmap_vx[1]; 
  const double vmap_vx2_sq = vmap_vx[2]*vmap_vx[2]; 
  const double vmap_vx2_cu = vmap_vx[2]*vmap_vx[2]*vmap_vx[2]; 
  const double vmap_vx3_sq = vmap_vx[3]*vmap_vx[3]; 
  const double vmap_vx3_cu = vmap_vx[3]*vmap_vx[3]*vmap_vx[3]; 
  out[0] += (2.772727272727273*vmap_vx[2]*f[7]*vmap_vx3_sq+0.8944271909999159*vmap_vx[0]*f[7]*vmap_vx3_sq+2.3*vmap_vx[1]*f[3]*vmap_vx3_sq+0.8944271909999159*f[0]*vmap_vx[2]*vmap_vx3_sq+1.5*f[0]*vmap_vx[0]*vmap_vx3_sq+0.499917348540637*f[3]*vmap_vx3_cu+0.9583148474999099*vmap_vx[0]*f[7]*vmap_vx2_sq+1.963961012123931*f[3]*vmap_vx[3]*vmap_vx2_sq+2.357142857142857*vmap_vx[1]*f[3]*vmap_vx2_sq+1.5*f[0]*vmap_vx[0]*vmap_vx2_sq+1.0714285714285714*f[7]*vmap_vx2_cu+0.31943828249996997*f[0]*vmap_vx2_cu+2.357142857142857*vmap_vx[2]*f[7]*vmap_vx1_sq+1.3416407864998738*vmap_vx[0]*f[7]*vmap_vx1_sq+1.1783766072743584*f[3]*vmap_vx[3]*vmap_vx1_sq+1.3416407864998738*f[0]*vmap_vx[2]*vmap_vx1_sq+1.5*f[0]*vmap_vx[0]*vmap_vx1_sq+0.9*f[3]*vmap_vx1_cu+1.5*vmap_vx[2]*f[7]*vmap_vx0_sq+1.5*vmap_vx[1]*f[3]*vmap_vx0_sq+0.5*f[0]*vmap_vx0_cu+3.927922024247862*vmap_vx[1]*vmap_vx[2]*vmap_vx[3]*f[7]+2.6349301969610384*vmap_vx[0]*vmap_vx[1]*vmap_vx[3]*f[7]+2.6349301969610384*vmap_vx[0]*vmap_vx[2]*f[3]*vmap_vx[3]+2.6349301969610384*f[0]*vmap_vx[1]*vmap_vx[2]*vmap_vx[3]+2.6832815729997477*vmap_vx[0]*vmap_vx[1]*vmap_vx[2]*f[3])*volFact; 
  out[1] += (2.7727272727272734*vmap_vx[2]*f[9]*vmap_vx3_sq+0.8944271909999161*vmap_vx[0]*f[9]*vmap_vx3_sq+2.3*vmap_vx[1]*f[5]*vmap_vx3_sq+0.8944271909999159*f[1]*vmap_vx[2]*vmap_vx3_sq+1.5*vmap_vx[0]*f[1]*vmap_vx3_sq+0.499917348540637*f[5]*vmap_vx3_cu+0.9583148474999098*vmap_vx[0]*f[9]*vmap_vx2_sq+1.963961012123931*vmap_vx[3]*f[5]*vmap_vx2_sq+2.357142857142857*vmap_vx[1]*f[5]*vmap_vx2_sq+1.5*vmap_vx[0]*f[1]*vmap_vx2_sq+1.0714285714285714*f[9]*vmap_vx2_cu+0.31943828249996997*f[1]*vmap_vx2_cu+2.357142857142857*vmap_vx[2]*f[9]*vmap_vx1_sq+1.3416407864998738*vmap_vx[0]*f[9]*vmap_vx1_sq+1.1783766072743584*vmap_vx[3]*f[5]*vmap_vx1_sq+1.3416407864998738*f[1]*vmap_vx[2]*vmap_vx1_sq+1.5*vmap_vx[0]*f[1]*vmap_vx1_sq+0.9*f[5]*vmap_vx1_cu+1.5*vmap_vx[2]*f[9]*vmap_vx0_sq+1.5*vmap_vx[1]*f[5]*vmap_vx0_sq+0.5*f[1]*vmap_vx0_cu+3.927922024247863*vmap_vx[1]*vmap_vx[2]*vmap_vx[3]*f[9]+2.6349301969610397*vmap_vx[0]*vmap_vx[1]*vmap_vx[3]*f[9]+2.6349301969610384*vmap_vx[0]*vmap_vx[2]*vmap_vx[3]*f[5]+2.6832815729997477*vmap_vx[0]*vmap_vx[1]*vmap_vx[2]*f[5]+2.6349301969610384*f[1]*vmap_vx[1]*vmap_vx[2]*vmap_vx[3])*volFact; 
  out[2] += (2.7727272727272734*vmap_vx[2]*f[10]*vmap_vx3_sq+0.8944271909999161*vmap_vx[0]*f[10]*vmap_vx3_sq+2.3*vmap_vx[1]*f[6]*vmap_vx3_sq+0.8944271909999159*f[2]*vmap_vx[2]*vmap_vx3_sq+1.5*vmap_vx[0]*f[2]*vmap_vx3_sq+0.499917348540637*f[6]*vmap_vx3_cu+0.9583148474999098*vmap_vx[0]*f[10]*vmap_vx2_sq+1.963961012123931*vmap_vx[3]*f[6]*vmap_vx2_sq+2.357142857142857*vmap_vx[1]*f[6]*vmap_vx2_sq+1.5*vmap_vx[0]*f[2]*vmap_vx2_sq+1.0714285714285714*f[10]*vmap_vx2_cu+0.31943828249996997*f[2]*vmap_vx2_cu+2.357142857142857*vmap_vx[2]*f[10]*vmap_vx1_sq+1.3416407864998738*vmap_vx[0]*f[10]*vmap_vx1_sq+1.1783766072743584*vmap_vx[3]*f[6]*vmap_vx1_sq+1.3416407864998738*f[2]*vmap_vx[2]*vmap_vx1_sq+1.5*vmap_vx[0]*f[2]*vmap_vx1_sq+0.9*f[6]*vmap_vx1_cu+1.5*vmap_vx[2]*f[10]*vmap_vx0_sq+1.5*vmap_vx[1]*f[6]*vmap_vx0_sq+0.5*f[2]*vmap_vx0_cu+3.927922024247863*vmap_vx[1]*vmap_vx[2]*vmap_vx[3]*f[10]+2.6349301969610397*vmap_vx[0]*vmap_vx[1]*vmap_vx[3]*f[10]+2.6349301969610384*vmap_vx[0]*vmap_vx[2]*vmap_vx[3]*f[6]+2.6832815729997477*vmap_vx[0]*vmap_vx[1]*vmap_vx[2]*f[6]+2.6349301969610384*vmap_vx[1]*f[2]*vmap_vx[2]*vmap_vx[3])*volFact; 
  out[3] += (2.772727272727273*vmap_vx[2]*f[11]*vmap_vx3_sq+0.8944271909999159*vmap_vx[0]*f[11]*vmap_vx3_sq+2.3*vmap_vx[1]*f[8]*vmap_vx3_sq+0.8944271909999159*vmap_vx[2]*f[4]*vmap_vx3_sq+1.5*vmap_vx[0]*f[4]*vmap_vx3_sq+0.499917348540637*f[8]*vmap_vx3_cu+0.9583148474999099*vmap_vx[0]*f[11]*vmap_vx2_sq+1.963961012123931*vmap_vx[3]*f[8]*vmap_vx2_sq+2.357142857142857*vmap_vx[1]*f[8]*vmap_vx2_sq+1.5*vmap_vx[0]*f[4]*vmap_vx2_sq+1.0714285714285714*f[11]*vmap_vx2_cu+0.31943828249996997*f[4]*vmap_vx2_cu+2.357142857142857*vmap_vx[2]*f[11]*vmap_vx1_sq+1.3416407864998738*vmap_vx[0]*f[11]*vmap_vx1_sq+1.1783766072743584*vmap_vx[3]*f[8]*vmap_vx1_sq+1.3416407864998738*vmap_vx[2]*f[4]*vmap_vx1_sq+1.5*vmap_vx[0]*f[4]*vmap_vx1_sq+0.9*f[8]*vmap_vx1_cu+1.5*vmap_vx[2]*f[11]*vmap_vx0_sq+1.5*vmap_vx[1]*f[8]*vmap_vx0_sq+0.5*f[4]*vmap_vx0_cu+3.927922024247862*vmap_vx[1]*vmap_vx[2]*vmap_vx[3]*f[11]+2.6349301969610384*vmap_vx[0]*vmap_vx[1]*vmap_vx[3]*f[11]+2.6349301969610384*vmap_vx[0]*vmap_vx[2]*vmap_vx[3]*f[8]+2.6832815729997477*vmap_vx[0]*vmap_vx[1]*vmap_vx[2]*f[8]+2.6349301969610384*vmap_vx[1]*vmap_vx[2]*vmap_vx[3]*f[4])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_M1i_2x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]/2; 
  const double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  double bdiv_vx[9]; 
  bdiv_vx[0] = 0.2777777777777778/jacob_vx[2]+0.4444444444444444/jacob_vx[1]+0.2777777777777778/jacob_vx[0]; 
  bdiv_vx[1] = 0.37267799624996495/jacob_vx[2]-0.37267799624996495/jacob_vx[0]; 
  bdiv_vx[2] = 0.24845199749997662/jacob_vx[2]-0.49690399499995325/jacob_vx[1]+0.24845199749997662/jacob_vx[0]; 
  bdiv_vx[3] = 0.37267799624996495/jacob_vx[2]-0.37267799624996495/jacob_vx[0]; 
  bdiv_vx[4] = 0.5/jacob_vx[2]+0.5/jacob_vx[0]; 
  bdiv_vx[5] = 0.3333333333333333/jacob_vx[2]-0.3333333333333333/jacob_vx[0]; 
  bdiv_vx[6] = 0.24845199749997662/jacob_vx[2]-0.49690399499995325/jacob_vx[1]+0.24845199749997662/jacob_vx[0]; 
  bdiv_vx[7] = 0.3333333333333333/jacob_vx[2]-0.3333333333333333/jacob_vx[0]; 
  bdiv_vx[8] = 0.2222222222222222/jacob_vx[2]+0.5555555555555556/jacob_vx[1]+0.2222222222222222/jacob_vx[0]; 
  out[0] += 3.872983346207417*hamil[2]*bdiv_vx[5]*f[7]*dv10*volFact+1.7320508075688772*hamil[1]*bdiv_vx[2]*f[7]*dv10*volFact+3.872983346207417*hamil[2]*f[3]*bdiv_vx[4]*dv10*volFact+1.7320508075688772*bdiv_vx[1]*hamil[1]*f[3]*dv10*volFact+3.872983346207417*f[0]*hamil[2]*bdiv_vx[3]*dv10*volFact+1.7320508075688772*bdiv_vx[0]*f[0]*hamil[1]*dv10*volFact; 
  out[1] += 3.872983346207417*hamil[2]*bdiv_vx[5]*f[9]*dv10*volFact+1.7320508075688772*hamil[1]*bdiv_vx[2]*f[9]*dv10*volFact+3.872983346207417*hamil[2]*bdiv_vx[4]*f[5]*dv10*volFact+1.7320508075688772*bdiv_vx[1]*hamil[1]*f[5]*dv10*volFact+3.872983346207417*f[1]*hamil[2]*bdiv_vx[3]*dv10*volFact+1.7320508075688772*bdiv_vx[0]*f[1]*hamil[1]*dv10*volFact; 
  out[2] += 3.872983346207417*hamil[2]*bdiv_vx[5]*f[10]*dv10*volFact+1.7320508075688772*hamil[1]*bdiv_vx[2]*f[10]*dv10*volFact+3.872983346207417*hamil[2]*bdiv_vx[4]*f[6]*dv10*volFact+1.7320508075688772*bdiv_vx[1]*hamil[1]*f[6]*dv10*volFact+3.872983346207417*f[2]*hamil[2]*bdiv_vx[3]*dv10*volFact+1.7320508075688772*bdiv_vx[0]*hamil[1]*f[2]*dv10*volFact; 
  out[3] += 3.872983346207417*hamil[2]*bdiv_vx[5]*f[11]*dv10*volFact+1.7320508075688772*hamil[1]*bdiv_vx[2]*f[11]*dv10*volFact+3.872983346207417*hamil[2]*bdiv_vx[4]*f[8]*dv10*volFact+1.7320508075688772*bdiv_vx[1]*hamil[1]*f[8]*dv10*volFact+3.872983346207417*hamil[2]*bdiv_vx[3]*f[4]*dv10*volFact+1.7320508075688772*bdiv_vx[0]*hamil[1]*f[4]*dv10*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_M2_2x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]/2; 
  out[0] += (hamil[2]*f[7]+hamil[1]*f[3]+f[0]*hamil[0])*volFact; 
  out[1] += (1.0000000000000002*hamil[2]*f[9]+hamil[1]*f[5]+hamil[0]*f[1])*volFact; 
  out[2] += (1.0000000000000002*hamil[2]*f[10]+hamil[1]*f[6]+hamil[0]*f[2])*volFact; 
  out[3] += (hamil[2]*f[11]+hamil[1]*f[8]+hamil[0]*f[4])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_five_moments_2x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]/2; 
  const double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  double bdiv_vx[9]; 
  bdiv_vx[0] = 0.2777777777777778/jacob_vx[2]+0.4444444444444444/jacob_vx[1]+0.2777777777777778/jacob_vx[0]; 
  bdiv_vx[1] = 0.37267799624996495/jacob_vx[2]-0.37267799624996495/jacob_vx[0]; 
  bdiv_vx[2] = 0.24845199749997662/jacob_vx[2]-0.49690399499995325/jacob_vx[1]+0.24845199749997662/jacob_vx[0]; 
  bdiv_vx[3] = 0.37267799624996495/jacob_vx[2]-0.37267799624996495/jacob_vx[0]; 
  bdiv_vx[4] = 0.5/jacob_vx[2]+0.5/jacob_vx[0]; 
  bdiv_vx[5] = 0.3333333333333333/jacob_vx[2]-0.3333333333333333/jacob_vx[0]; 
  bdiv_vx[6] = 0.24845199749997662/jacob_vx[2]-0.49690399499995325/jacob_vx[1]+0.24845199749997662/jacob_vx[0]; 
  bdiv_vx[7] = 0.3333333333333333/jacob_vx[2]-0.3333333333333333/jacob_vx[0]; 
  bdiv_vx[8] = 0.2222222222222222/jacob_vx[2]+0.5555555555555556/jacob_vx[1]+0.2222222222222222/jacob_vx[0]; 
  out[0] += 1.4142135623730951*f[0]*volFact; 
  out[4] += 3.872983346207417*hamil[2]*bdiv_vx[5]*f[7]*dv10*volFact+1.7320508075688772*hamil[1]*bdiv_vx[2]*f[7]*dv10*volFact+3.872983346207417*hamil[2]*f[3]*bdiv_vx[4]*dv10*volFact+1.7320508075688772*bdiv_vx[1]*hamil[1]*f[3]*dv10*volFact+3.872983346207417*f[0]*hamil[2]*bdiv_vx[3]*dv10*volFact+1.7320508075688772*bdiv_vx[0]*f[0]*hamil[1]*dv10*volFact; 
  out[8] += hamil[2]*f[7]*volFact+hamil[1]*f[3]*volFact+f[0]*hamil[0]*volFact; 
  out[1] += 1.4142135623730951*f[1]*volFact; 
  out[5] += 3.872983346207417*hamil[2]*bdiv_vx[5]*f[9]*dv10*volFact+1.7320508075688772*hamil[1]*bdiv_vx[2]*f[9]*dv10*volFact+3.872983346207417*hamil[2]*bdiv_vx[4]*f[5]*dv10*volFact+1.7320508075688772*bdiv_vx[1]*hamil[1]*f[5]*dv10*volFact+3.872983346207417*f[1]*hamil[2]*bdiv_vx[3]*dv10*volFact+1.7320508075688772*bdiv_vx[0]*f[1]*hamil[1]*dv10*volFact; 
  out[9] += hamil[2]*f[9]*volFact+hamil[1]*f[5]*volFact+hamil[0]*f[1]*volFact; 
  out[2] += 1.4142135623730951*f[2]*volFact; 
  out[6] += 3.872983346207417*hamil[2]*bdiv_vx[5]*f[10]*dv10*volFact+1.7320508075688772*hamil[1]*bdiv_vx[2]*f[10]*dv10*volFact+3.872983346207417*hamil[2]*bdiv_vx[4]*f[6]*dv10*volFact+1.7320508075688772*bdiv_vx[1]*hamil[1]*f[6]*dv10*volFact+3.872983346207417*f[2]*hamil[2]*bdiv_vx[3]*dv10*volFact+1.7320508075688772*bdiv_vx[0]*hamil[1]*f[2]*dv10*volFact; 
  out[10] += hamil[2]*f[10]*volFact+hamil[1]*f[6]*volFact+hamil[0]*f[2]*volFact; 
  out[3] += 1.4142135623730951*f[4]*volFact; 
  out[7] += 3.872983346207417*hamil[2]*bdiv_vx[5]*f[11]*dv10*volFact+1.7320508075688772*hamil[1]*bdiv_vx[2]*f[11]*dv10*volFact+3.872983346207417*hamil[2]*bdiv_vx[4]*f[8]*dv10*volFact+1.7320508075688772*bdiv_vx[1]*hamil[1]*f[8]*dv10*volFact+3.872983346207417*hamil[2]*bdiv_vx[3]*f[4]*dv10*volFact+1.7320508075688772*bdiv_vx[0]*hamil[1]*f[4]*dv10*volFact; 
  out[11] += hamil[2]*f[11]*volFact+hamil[1]*f[8]*volFact+hamil[0]*f[4]*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_int_five_moments_2x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[0]*dxv[1]*dxv[2]*0.125; 
  const double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  double bdiv_vx[9]; 
  bdiv_vx[0] = 0.2777777777777778/jacob_vx[2]+0.4444444444444444/jacob_vx[1]+0.2777777777777778/jacob_vx[0]; 
  bdiv_vx[1] = 0.37267799624996495/jacob_vx[2]-0.37267799624996495/jacob_vx[0]; 
  bdiv_vx[2] = 0.24845199749997662/jacob_vx[2]-0.49690399499995325/jacob_vx[1]+0.24845199749997662/jacob_vx[0]; 
  bdiv_vx[3] = 0.37267799624996495/jacob_vx[2]-0.37267799624996495/jacob_vx[0]; 
  bdiv_vx[4] = 0.5/jacob_vx[2]+0.5/jacob_vx[0]; 
  bdiv_vx[5] = 0.3333333333333333/jacob_vx[2]-0.3333333333333333/jacob_vx[0]; 
  bdiv_vx[6] = 0.24845199749997662/jacob_vx[2]-0.49690399499995325/jacob_vx[1]+0.24845199749997662/jacob_vx[0]; 
  bdiv_vx[7] = 0.3333333333333333/jacob_vx[2]-0.3333333333333333/jacob_vx[0]; 
  bdiv_vx[8] = 0.2222222222222222/jacob_vx[2]+0.5555555555555556/jacob_vx[1]+0.2222222222222222/jacob_vx[0]; 
  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += 7.745966692414834*hamil[2]*bdiv_vx[5]*f[7]*dv10*volFact+3.4641016151377544*hamil[1]*bdiv_vx[2]*f[7]*dv10*volFact+7.745966692414834*hamil[2]*f[3]*bdiv_vx[4]*dv10*volFact+3.4641016151377544*bdiv_vx[1]*hamil[1]*f[3]*dv10*volFact+7.745966692414834*f[0]*hamil[2]*bdiv_vx[3]*dv10*volFact+3.4641016151377544*bdiv_vx[0]*f[0]*hamil[1]*dv10*volFact; 
  out[2] += 2.0*hamil[2]*f[7]*volFact+2.0*hamil[1]*f[3]*volFact+2.0*f[0]*hamil[0]*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M0_upper_2x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *vmap, const double *jacob_vel, double v_thresh, double f_thresh, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]/2; 
  const double *px = &vmap[0]; 
  const double *jacob_vx = &jacob_vel[0]; 
  double f_nodes[8] = {0.0};
  double fquad = 0.0;
  if ((0.7200822998230954*px[3]-0.7071067811865475*px[1]+0.7071067811865475*px[0]) > v_thresh) { 
    fquad = -(0.3535533905932737*f[8])+0.3535533905932737*(f[6]+f[5]+f[4])-0.3535533905932737*(f[3]+f[2]+f[1])+0.3535533905932737*f[0]; 
    f_nodes[0] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
    fquad = 0.3535533905932737*f[8]-0.3535533905932737*f[6]+0.3535533905932737*f[5]-0.3535533905932737*(f[4]+f[3])+0.3535533905932737*f[2]-0.3535533905932737*f[1]+0.3535533905932737*f[0]; 
    f_nodes[2] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
    fquad = 0.3535533905932737*(f[8]+f[6])-0.3535533905932737*(f[5]+f[4]+f[3]+f[2])+0.3535533905932737*(f[1]+f[0]); 
    f_nodes[4] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
    fquad = -(0.3535533905932737*(f[8]+f[6]+f[5]))+0.3535533905932737*f[4]-0.3535533905932737*f[3]+0.3535533905932737*(f[2]+f[1]+f[0]); 
    f_nodes[6] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
  } 
  if ((0.7071067811865475*(px[1]+px[0])-0.7200822998230954*px[3]) > v_thresh) { 
    fquad = 0.3535533905932737*f[8]-0.3535533905932737*(f[6]+f[5])+0.3535533905932737*(f[4]+f[3])-0.3535533905932737*(f[2]+f[1])+0.3535533905932737*f[0]; 
    f_nodes[1] = fquad > jacob_vx[1]*f_thresh ? fquad : 0.0; 
    fquad = -(0.3535533905932737*f[8])+0.3535533905932737*f[6]-0.3535533905932737*(f[5]+f[4])+0.3535533905932737*(f[3]+f[2])-0.3535533905932737*f[1]+0.3535533905932737*f[0]; 
    f_nodes[3] = fquad > jacob_vx[1]*f_thresh ? fquad : 0.0; 
    fquad = -(0.3535533905932737*(f[8]+f[6]))+0.3535533905932737*f[5]-0.3535533905932737*f[4]+0.3535533905932737*f[3]-0.3535533905932737*f[2]+0.3535533905932737*(f[1]+f[0]); 
    f_nodes[5] = fquad > jacob_vx[1]*f_thresh ? fquad : 0.0; 
    fquad = 0.3535533905932737*(f[8]+f[6]+f[5]+f[4]+f[3]+f[2]+f[1]+f[0]); 
    f_nodes[7] = fquad > jacob_vx[1]*f_thresh ? fquad : 0.0; 
  } 
  out[0] += (0.5*f_nodes[7]+0.5*f_nodes[6]+0.5*f_nodes[5]+0.5*f_nodes[4]+0.5*f_nodes[3]+0.5*f_nodes[2]+0.5*f_nodes[1]+0.5*f_nodes[0])*volFact; 
  out[1] += (0.5*f_nodes[7]+0.5*f_nodes[6]+0.5*f_nodes[5]+0.5*f_nodes[4]-0.5*f_nodes[3]-0.5*f_nodes[2]-0.5*f_nodes[1]-0.5*f_nodes[0])*volFact; 
  out[2] += (0.5*f_nodes[7]+0.5*f_nodes[6]-0.5*f_nodes[5]-0.5*f_nodes[4]+0.5*f_nodes[3]+0.5*f_nodes[2]-0.5*f_nodes[1]-0.5*f_nodes[0])*volFact; 
  out[3] += (0.5*f_nodes[7]+0.5*f_nodes[6]-0.5*f_nodes[5]-0.5*f_nodes[4]-0.5*f_nodes[3]-0.5*f_nodes[2]+0.5*f_nodes[1]+0.5*f_nodes[0])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M0_lower_2x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *vmap, const double *jacob_vel, double v_thresh, double f_thresh, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]/2; 
  const double *px = &vmap[0]; 
  const double *jacob_vx = &jacob_vel[0]; 
  double f_nodes[8] = {0.0};
  double fquad = 0.0;
  if ((0.7200822998230954*px[3]-0.7071067811865475*px[1]+0.7071067811865475*px[0]) < -v_thresh) { 
    fquad = -(0.3535533905932737*f[8])+0.3535533905932737*(f[6]+f[5]+f[4])-0.3535533905932737*(f[3]+f[2]+f[1])+0.3535533905932737*f[0]; 
    f_nodes[0] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
    fquad = 0.3535533905932737*f[8]-0.3535533905932737*f[6]+0.3535533905932737*f[5]-0.3535533905932737*(f[4]+f[3])+0.3535533905932737*f[2]-0.3535533905932737*f[1]+0.3535533905932737*f[0]; 
    f_nodes[2] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
    fquad = 0.3535533905932737*(f[8]+f[6])-0.3535533905932737*(f[5]+f[4]+f[3]+f[2])+0.3535533905932737*(f[1]+f[0]); 
    f_nodes[4] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
    fquad = -(0.3535533905932737*(f[8]+f[6]+f[5]))+0.3535533905932737*f[4]-0.3535533905932737*f[3]+0.3535533905932737*(f[2]+f[1]+f[0]); 
    f_nodes[6] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
  } 
  if ((0.7071067811865475*(px[1]+px[0])-0.7200822998230954*px[3]) < -v_thresh) { 
    fquad = 0.3535533905932737*f[8]-0.3535533905932737*(f[6]+f[5])+0.3535533905932737*(f[4]+f[3])-0.3535533905932737*(f[2]+f[1])+0.3535533905932737*f[0]; 
    f_nodes[1] = fquad > jacob_vx[1]*f_thresh ? fquad : 0.0; 
    fquad = -(0.3535533905932737*f[8])+0.3535533905932737*f[6]-0.3535533905932737*(f[5]+f[4])+0.3535533905932737*(f[3]+f[2])-0.3535533905932737*f[1]+0.3535533905932737*f[0]; 
    f_nodes[3] = fquad > jacob_vx[1]*f_thresh ? fquad : 0.0; 
    fquad = -(0.3535533905932737*(f[8]+f[6]))+0.3535533905932737*f[5]-0.3535533905932737*f[4]+0.3535533905932737*f[3]-0.3535533905932737*f[2]+0.3535533905932737*(f[1]+f[0]); 
    f_nodes[5] = fquad > jacob_vx[1]*f_thresh ? fquad : 0.0; 
    fquad = 0.3535533905932737*(f[8]+f[6]+f[5]+f[4]+f[3]+f[2]+f[1]+f[0]); 
    f_nodes[7] = fquad > jacob_vx[1]*f_thresh ? fquad : 0.0; 
  } 
  out[0] += (0.5*f_nodes[7]+0.5*f_nodes[6]+0.5*f_nodes[5]+0.5*f_nodes[4]+0.5*f_nodes[3]+0.5*f_nodes[2]+0.5*f_nodes[1]+0.5*f_nodes[0])*volFact; 
  out[1] += (0.5*f_nodes[7]+0.5*f_nodes[6]+0.5*f_nodes[5]+0.5*f_nodes[4]-0.5*f_nodes[3]-0.5*f_nodes[2]-0.5*f_nodes[1]-0.5*f_nodes[0])*volFact; 
  out[2] += (0.5*f_nodes[7]+0.5*f_nodes[6]-0.5*f_nodes[5]-0.5*f_nodes[4]+0.5*f_nodes[3]+0.5*f_nodes[2]-0.5*f_nodes[1]-0.5*f_nodes[0])*volFact; 
  out[3] += (0.5*f_nodes[7]+0.5*f_nodes[6]-0.5*f_nodes[5]-0.5*f_nodes[4]-0.5*f_nodes[3]-0.5*f_nodes[2]+0.5*f_nodes[1]+0.5*f_nodes[0])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_phase_M1i_2x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]/2; 
  const double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  double bdiv_vx[9]; 
  bdiv_vx[0] = 0.2777777777777778/jacob_vx[2]+0.4444444444444444/jacob_vx[1]+0.2777777777777778/jacob_vx[0]; 
  bdiv_vx[1] = 0.37267799624996495/jacob_vx[2]-0.37267799624996495/jacob_vx[0]; 
  bdiv_vx[2] = 0.24845199749997662/jacob_vx[2]-0.49690399499995325/jacob_vx[1]+0.24845199749997662/jacob_vx[0]; 
  bdiv_vx[3] = 0.37267799624996495/jacob_vx[2]-0.37267799624996495/jacob_vx[0]; 
  bdiv_vx[4] = 0.5/jacob_vx[2]+0.5/jacob_vx[0]; 
  bdiv_vx[5] = 0.3333333333333333/jacob_vx[2]-0.3333333333333333/jacob_vx[0]; 
  bdiv_vx[6] = 0.24845199749997662/jacob_vx[2]-0.49690399499995325/jacob_vx[1]+0.24845199749997662/jacob_vx[0]; 
  bdiv_vx[7] = 0.3333333333333333/jacob_vx[2]-0.3333333333333333/jacob_vx[0]; 
  bdiv_vx[8] = 0.2222222222222222/jacob_vx[2]+0.5555555555555556/jacob_vx[1]+0.2222222222222222/jacob_vx[0]; 
  for (int q = 0; q < 4; ++q) { 
  const double s0 = ((q >> 1) & 1) ? 1.0 : -1.0; 
  const double s1 = ((q >> 0) & 1) ? 1.0 : -1.0; 
  double f_node[3], h_node[3]; 
  f_node[0] = 0.5*f[4]*s0*s1+0.5*f[2]*s1+0.5*f[1]*s0+0.5*f[0]; 
  f_node[1] = 0.5*f[8]*s0*s1+0.5*f[6]*s1+0.5*f[5]*s0+0.5*f[3]; 
  f_node[2] = 0.5*f[11]*s0*s1+0.5*f[10]*s1+0.5*f[9]*s0+0.5*f[7]; 
  h_node[0] = 0.5*hamil[4]*s0*s1+0.5*hamil[2]*s1+0.5*hamil[1]*s0+0.5*hamil[0]; 
  h_node[1] = 0.5*hamil[8]*s0*s1+0.5*hamil[6]*s1+0.5*hamil[5]*s0+0.5*hamil[3]; 
  h_node[2] = 0.5*hamil[11]*s0*s1+0.5*hamil[10]*s1+0.5*hamil[9]*s0+0.5*hamil[7]; 
  const double m1_x = 3.872983346207417*f_node[2]*h_node[2]*bdiv_vx[5]*dv10+3.872983346207417*f_node[1]*h_node[2]*bdiv_vx[4]*dv10+3.872983346207417*f_node[0]*h_node[2]*bdiv_vx[3]*dv10+1.7320508075688772*h_node[1]*bdiv_vx[2]*f_node[2]*dv10+1.7320508075688772*bdiv_vx[1]*f_node[1]*h_node[1]*dv10+1.7320508075688772*bdiv_vx[0]*f_node[0]*h_node[1]*dv10; 
  out[0] += 0.5*m1_x*volFact; 
  out[1] += 0.5*m1_x*s0*volFact; 
  out[2] += 0.5*m1_x*s1*volFact; 
  out[3] += 0.5*m1_x*s0*s1*volFact; 
  } 
} 
GKYL_CU_DH void mom_vlasov_hamil_phase_M2_2x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]/2; 
  out[0] += (0.5*f[11]*hamil[11]+0.5*f[10]*hamil[10]+0.5*f[9]*hamil[9]+0.5*f[8]*hamil[8]+0.5*f[7]*hamil[7]+0.5*f[6]*hamil[6]+0.5*f[5]*hamil[5]+0.5*f[4]*hamil[4]+0.5*f[3]*hamil[3]+0.5*f[2]*hamil[2]+0.5*f[1]*hamil[1]+0.5*f[0]*hamil[0])*volFact; 
  out[1] += (0.5000000000000001*f[10]*hamil[11]+0.5000000000000001*hamil[10]*f[11]+0.5000000000000001*f[7]*hamil[9]+0.5000000000000001*hamil[7]*f[9]+0.5*f[6]*hamil[8]+0.5*hamil[6]*f[8]+0.5*f[3]*hamil[5]+0.5*hamil[3]*f[5]+0.5*f[2]*hamil[4]+0.5*hamil[2]*f[4]+0.5*f[0]*hamil[1]+0.5*hamil[0]*f[1])*volFact; 
  out[2] += (0.5000000000000001*f[9]*hamil[11]+0.5000000000000001*hamil[9]*f[11]+0.5000000000000001*f[7]*hamil[10]+0.5000000000000001*hamil[7]*f[10]+0.5*f[5]*hamil[8]+0.5*hamil[5]*f[8]+0.5*f[3]*hamil[6]+0.5*hamil[3]*f[6]+0.5*f[1]*hamil[4]+0.5*hamil[1]*f[4]+0.5*f[0]*hamil[2]+0.5*hamil[0]*f[2])*volFact; 
  out[3] += (0.5*f[7]*hamil[11]+0.5*hamil[7]*f[11]+0.5*f[9]*hamil[10]+0.5*hamil[9]*f[10]+0.5*f[3]*hamil[8]+0.5*hamil[3]*f[8]+0.5*f[5]*hamil[6]+0.5*hamil[5]*f[6]+0.5*f[0]*hamil[4]+0.5*hamil[0]*f[4]+0.5*f[1]*hamil[2]+0.5*hamil[1]*f[2])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_phase_five_moments_2x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]/2; 
  const double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  double bdiv_vx[9]; 
  bdiv_vx[0] = 0.2777777777777778/jacob_vx[2]+0.4444444444444444/jacob_vx[1]+0.2777777777777778/jacob_vx[0]; 
  bdiv_vx[1] = 0.37267799624996495/jacob_vx[2]-0.37267799624996495/jacob_vx[0]; 
  bdiv_vx[2] = 0.24845199749997662/jacob_vx[2]-0.49690399499995325/jacob_vx[1]+0.24845199749997662/jacob_vx[0]; 
  bdiv_vx[3] = 0.37267799624996495/jacob_vx[2]-0.37267799624996495/jacob_vx[0]; 
  bdiv_vx[4] = 0.5/jacob_vx[2]+0.5/jacob_vx[0]; 
  bdiv_vx[5] = 0.3333333333333333/jacob_vx[2]-0.3333333333333333/jacob_vx[0]; 
  bdiv_vx[6] = 0.24845199749997662/jacob_vx[2]-0.49690399499995325/jacob_vx[1]+0.24845199749997662/jacob_vx[0]; 
  bdiv_vx[7] = 0.3333333333333333/jacob_vx[2]-0.3333333333333333/jacob_vx[0]; 
  bdiv_vx[8] = 0.2222222222222222/jacob_vx[2]+0.5555555555555556/jacob_vx[1]+0.2222222222222222/jacob_vx[0]; 
  out[0] += 1.4142135623730951*f[0]*volFact; 
  out[1] += 1.4142135623730951*f[1]*volFact; 
  out[2] += 1.4142135623730951*f[2]*volFact; 
  out[3] += 1.4142135623730951*f[4]*volFact; 
  for (int q = 0; q < 4; ++q) { 
  const double s0 = ((q >> 1) & 1) ? 1.0 : -1.0; 
  const double s1 = ((q >> 0) & 1) ? 1.0 : -1.0; 
  double f_node[3], h_node[3]; 
  f_node[0] = 0.5*f[4]*s0*s1+0.5*f[2]*s1+0.5*f[1]*s0+0.5*f[0]; 
  f_node[1] = 0.5*f[8]*s0*s1+0.5*f[6]*s1+0.5*f[5]*s0+0.5*f[3]; 
  f_node[2] = 0.5*f[11]*s0*s1+0.5*f[10]*s1+0.5*f[9]*s0+0.5*f[7]; 
  h_node[0] = 0.5*hamil[4]*s0*s1+0.5*hamil[2]*s1+0.5*hamil[1]*s0+0.5*hamil[0]; 
  h_node[1] = 0.5*hamil[8]*s0*s1+0.5*hamil[6]*s1+0.5*hamil[5]*s0+0.5*hamil[3]; 
  h_node[2] = 0.5*hamil[11]*s0*s1+0.5*hamil[10]*s1+0.5*hamil[9]*s0+0.5*hamil[7]; 
  const double m1_x = 3.872983346207417*f_node[2]*h_node[2]*bdiv_vx[5]*dv10+3.872983346207417*f_node[1]*h_node[2]*bdiv_vx[4]*dv10+3.872983346207417*f_node[0]*h_node[2]*bdiv_vx[3]*dv10+1.7320508075688772*h_node[1]*bdiv_vx[2]*f_node[2]*dv10+1.7320508075688772*bdiv_vx[1]*f_node[1]*h_node[1]*dv10+1.7320508075688772*bdiv_vx[0]*f_node[0]*h_node[1]*dv10; 
  const double m2q = f_node[2]*h_node[2]+f_node[1]*h_node[1]+f_node[0]*h_node[0]; 
  out[4] += 0.5*m1_x*volFact; 
  out[8] += 0.5*m2q*volFact; 
  out[5] += 0.5*m1_x*s0*volFact; 
  out[9] += 0.5*m2q*s0*volFact; 
  out[6] += 0.5*m1_x*s1*volFact; 
  out[10] += 0.5*m2q*s1*volFact; 
  out[7] += 0.5*m1_x*s0*s1*volFact; 
  out[11] += 0.5*m2q*s0*s1*volFact; 
  } 
} 
GKYL_CU_DH void mom_vlasov_hamil_phase_int_five_moments_2x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[0]*dxv[1]*dxv[2]*0.125; 
  const double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  double bdiv_vx[9]; 
  bdiv_vx[0] = 0.2777777777777778/jacob_vx[2]+0.4444444444444444/jacob_vx[1]+0.2777777777777778/jacob_vx[0]; 
  bdiv_vx[1] = 0.37267799624996495/jacob_vx[2]-0.37267799624996495/jacob_vx[0]; 
  bdiv_vx[2] = 0.24845199749997662/jacob_vx[2]-0.49690399499995325/jacob_vx[1]+0.24845199749997662/jacob_vx[0]; 
  bdiv_vx[3] = 0.37267799624996495/jacob_vx[2]-0.37267799624996495/jacob_vx[0]; 
  bdiv_vx[4] = 0.5/jacob_vx[2]+0.5/jacob_vx[0]; 
  bdiv_vx[5] = 0.3333333333333333/jacob_vx[2]-0.3333333333333333/jacob_vx[0]; 
  bdiv_vx[6] = 0.24845199749997662/jacob_vx[2]-0.49690399499995325/jacob_vx[1]+0.24845199749997662/jacob_vx[0]; 
  bdiv_vx[7] = 0.3333333333333333/jacob_vx[2]-0.3333333333333333/jacob_vx[0]; 
  bdiv_vx[8] = 0.2222222222222222/jacob_vx[2]+0.5555555555555556/jacob_vx[1]+0.2222222222222222/jacob_vx[0]; 
  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += 3.872983346207417*bdiv_vx[5]*f[11]*hamil[11]*dv10*volFact+3.872983346207417*bdiv_vx[4]*f[8]*hamil[11]*dv10*volFact+3.872983346207417*bdiv_vx[3]*f[4]*hamil[11]*dv10*volFact+1.7320508075688772*bdiv_vx[2]*hamil[8]*f[11]*dv10*volFact+3.872983346207417*bdiv_vx[5]*f[10]*hamil[10]*dv10*volFact+3.872983346207417*bdiv_vx[4]*f[6]*hamil[10]*dv10*volFact+3.872983346207417*f[2]*bdiv_vx[3]*hamil[10]*dv10*volFact+1.7320508075688772*bdiv_vx[2]*hamil[6]*f[10]*dv10*volFact+3.872983346207417*bdiv_vx[5]*f[9]*hamil[9]*dv10*volFact+3.872983346207417*bdiv_vx[4]*f[5]*hamil[9]*dv10*volFact+3.872983346207417*f[1]*bdiv_vx[3]*hamil[9]*dv10*volFact+1.7320508075688772*bdiv_vx[2]*hamil[5]*f[9]*dv10*volFact+1.7320508075688772*bdiv_vx[1]*f[8]*hamil[8]*dv10*volFact+1.7320508075688772*bdiv_vx[0]*f[4]*hamil[8]*dv10*volFact+3.872983346207417*bdiv_vx[5]*f[7]*hamil[7]*dv10*volFact+3.872983346207417*f[3]*bdiv_vx[4]*hamil[7]*dv10*volFact+3.872983346207417*f[0]*bdiv_vx[3]*hamil[7]*dv10*volFact+1.7320508075688772*bdiv_vx[2]*hamil[3]*f[7]*dv10*volFact+1.7320508075688772*bdiv_vx[1]*f[6]*hamil[6]*dv10*volFact+1.7320508075688772*bdiv_vx[0]*f[2]*hamil[6]*dv10*volFact+1.7320508075688772*bdiv_vx[1]*f[5]*hamil[5]*dv10*volFact+1.7320508075688772*bdiv_vx[0]*f[1]*hamil[5]*dv10*volFact+1.7320508075688772*bdiv_vx[1]*f[3]*hamil[3]*dv10*volFact+1.7320508075688772*bdiv_vx[0]*f[0]*hamil[3]*dv10*volFact; 
  out[2] += f[11]*hamil[11]*volFact+f[10]*hamil[10]*volFact+f[9]*hamil[9]*volFact+f[8]*hamil[8]*volFact+f[7]*hamil[7]*volFact+f[6]*hamil[6]*volFact+f[5]*hamil[5]*volFact+f[4]*hamil[4]*volFact+f[3]*hamil[3]*volFact+f[2]*hamil[2]*volFact+f[1]*hamil[1]*volFact+f[0]*hamil[0]*volFact; 
} 
