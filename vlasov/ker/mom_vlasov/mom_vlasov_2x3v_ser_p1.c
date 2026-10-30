#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void mom_vlasov_M0_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
  const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += 2.8284271247461907*f[1]*volFact; 
  out[2] += 2.8284271247461907*f[2]*volFact; 
  out[3] += 2.8284271247461907*f[6]*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M2ij_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
  const double *vmap, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  const double *vmap_vx = &vmap[0]; 
  const double *vmap_vy = &vmap[4]; 
  const double *vmap_vz = &vmap[8]; 
  const double vmap_vx0_sq = vmap_vx[0]*vmap_vx[0]; 
  const double vmap_vx1_sq = vmap_vx[1]*vmap_vx[1]; 
  const double vmap_vy0_sq = vmap_vy[0]*vmap_vy[0]; 
  const double vmap_vy1_sq = vmap_vy[1]*vmap_vy[1]; 
  const double vmap_vz0_sq = vmap_vz[0]*vmap_vz[0]; 
  const double vmap_vz1_sq = vmap_vz[1]*vmap_vz[1]; 
  out[0] += (1.4142135623730951*f[0]*vmap_vx1_sq+1.4142135623730951*f[0]*vmap_vx0_sq+2.8284271247461907*vmap_vx[0]*vmap_vx[1]*f[3])*volFact; 
  out[1] += (1.4142135623730951*f[1]*vmap_vx1_sq+1.4142135623730951*f[1]*vmap_vx0_sq+2.8284271247461907*vmap_vx[0]*vmap_vx[1]*f[7])*volFact; 
  out[2] += (1.4142135623730951*f[2]*vmap_vx1_sq+1.4142135623730951*f[2]*vmap_vx0_sq+2.8284271247461907*vmap_vx[0]*vmap_vx[1]*f[8])*volFact; 
  out[3] += (1.4142135623730951*f[6]*vmap_vx1_sq+1.4142135623730951*f[6]*vmap_vx0_sq+2.8284271247461907*vmap_vx[0]*vmap_vx[1]*f[16])*volFact; 
  out[4] += (1.4142135623730951*vmap_vx[1]*vmap_vy[1]*f[11]+1.4142135623730951*vmap_vx[0]*vmap_vy[1]*f[4]+1.4142135623730951*vmap_vy[0]*vmap_vx[1]*f[3]+1.4142135623730951*f[0]*vmap_vx[0]*vmap_vy[0])*volFact; 
  out[5] += (1.4142135623730951*vmap_vx[1]*vmap_vy[1]*f[18]+1.4142135623730951*vmap_vx[0]*vmap_vy[1]*f[9]+1.4142135623730951*vmap_vy[0]*vmap_vx[1]*f[7]+1.4142135623730951*vmap_vx[0]*vmap_vy[0]*f[1])*volFact; 
  out[6] += (1.4142135623730951*vmap_vx[1]*vmap_vy[1]*f[19]+1.4142135623730951*vmap_vx[0]*vmap_vy[1]*f[10]+1.4142135623730951*vmap_vy[0]*vmap_vx[1]*f[8]+1.4142135623730951*vmap_vx[0]*vmap_vy[0]*f[2])*volFact; 
  out[7] += (1.4142135623730951*vmap_vx[1]*vmap_vy[1]*f[26]+1.4142135623730951*vmap_vx[0]*vmap_vy[1]*f[17]+1.4142135623730951*vmap_vy[0]*vmap_vx[1]*f[16]+1.4142135623730951*vmap_vx[0]*vmap_vy[0]*f[6])*volFact; 
  out[8] += (1.4142135623730951*vmap_vx[1]*vmap_vz[1]*f[14]+1.4142135623730951*vmap_vx[0]*vmap_vz[1]*f[5]+1.4142135623730951*vmap_vz[0]*vmap_vx[1]*f[3]+1.4142135623730951*f[0]*vmap_vx[0]*vmap_vz[0])*volFact; 
  out[9] += (1.4142135623730951*vmap_vx[1]*vmap_vz[1]*f[21]+1.4142135623730951*vmap_vx[0]*vmap_vz[1]*f[12]+1.4142135623730951*vmap_vz[0]*vmap_vx[1]*f[7]+1.4142135623730951*vmap_vx[0]*vmap_vz[0]*f[1])*volFact; 
  out[10] += (1.4142135623730951*vmap_vx[1]*vmap_vz[1]*f[22]+1.4142135623730951*vmap_vx[0]*vmap_vz[1]*f[13]+1.4142135623730951*vmap_vz[0]*vmap_vx[1]*f[8]+1.4142135623730951*vmap_vx[0]*vmap_vz[0]*f[2])*volFact; 
  out[11] += (1.4142135623730951*vmap_vx[1]*vmap_vz[1]*f[27]+1.4142135623730951*vmap_vx[0]*vmap_vz[1]*f[20]+1.4142135623730951*vmap_vz[0]*vmap_vx[1]*f[16]+1.4142135623730951*vmap_vx[0]*vmap_vz[0]*f[6])*volFact; 
  out[12] += (1.4142135623730951*f[0]*vmap_vy1_sq+1.4142135623730951*f[0]*vmap_vy0_sq+2.8284271247461907*vmap_vy[0]*vmap_vy[1]*f[4])*volFact; 
  out[13] += (1.4142135623730951*f[1]*vmap_vy1_sq+1.4142135623730951*f[1]*vmap_vy0_sq+2.8284271247461907*vmap_vy[0]*vmap_vy[1]*f[9])*volFact; 
  out[14] += (1.4142135623730951*f[2]*vmap_vy1_sq+1.4142135623730951*f[2]*vmap_vy0_sq+2.8284271247461907*vmap_vy[0]*vmap_vy[1]*f[10])*volFact; 
  out[15] += (1.4142135623730951*f[6]*vmap_vy1_sq+1.4142135623730951*f[6]*vmap_vy0_sq+2.8284271247461907*vmap_vy[0]*vmap_vy[1]*f[17])*volFact; 
  out[16] += (1.4142135623730951*vmap_vy[1]*vmap_vz[1]*f[15]+1.4142135623730951*vmap_vy[0]*vmap_vz[1]*f[5]+1.4142135623730951*vmap_vz[0]*vmap_vy[1]*f[4]+1.4142135623730951*f[0]*vmap_vy[0]*vmap_vz[0])*volFact; 
  out[17] += (1.4142135623730951*vmap_vy[1]*vmap_vz[1]*f[23]+1.4142135623730951*vmap_vy[0]*vmap_vz[1]*f[12]+1.4142135623730951*vmap_vz[0]*vmap_vy[1]*f[9]+1.4142135623730951*vmap_vy[0]*vmap_vz[0]*f[1])*volFact; 
  out[18] += (1.4142135623730951*vmap_vy[1]*vmap_vz[1]*f[24]+1.4142135623730951*vmap_vy[0]*vmap_vz[1]*f[13]+1.4142135623730951*vmap_vz[0]*vmap_vy[1]*f[10]+1.4142135623730951*vmap_vy[0]*vmap_vz[0]*f[2])*volFact; 
  out[19] += (1.4142135623730951*vmap_vy[1]*vmap_vz[1]*f[28]+1.4142135623730951*vmap_vy[0]*vmap_vz[1]*f[20]+1.4142135623730951*vmap_vz[0]*vmap_vy[1]*f[17]+1.4142135623730951*vmap_vy[0]*vmap_vz[0]*f[6])*volFact; 
  out[20] += (1.4142135623730951*f[0]*vmap_vz1_sq+1.4142135623730951*f[0]*vmap_vz0_sq+2.8284271247461907*vmap_vz[0]*vmap_vz[1]*f[5])*volFact; 
  out[21] += (1.4142135623730951*f[1]*vmap_vz1_sq+1.4142135623730951*f[1]*vmap_vz0_sq+2.8284271247461907*vmap_vz[0]*vmap_vz[1]*f[12])*volFact; 
  out[22] += (1.4142135623730951*f[2]*vmap_vz1_sq+1.4142135623730951*f[2]*vmap_vz0_sq+2.8284271247461907*vmap_vz[0]*vmap_vz[1]*f[13])*volFact; 
  out[23] += (1.4142135623730951*f[6]*vmap_vz1_sq+1.4142135623730951*f[6]*vmap_vz0_sq+2.8284271247461907*vmap_vz[0]*vmap_vz[1]*f[20])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M3ijk_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
  const double *vmap, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
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
  out[0] += (3.0*f[0]*vmap_vx[0]*vmap_vx1_sq+1.8*f[3]*vmap_vx1_cu+3.0*vmap_vx[1]*f[3]*vmap_vx0_sq+f[0]*vmap_vx0_cu)*volFact; 
  out[1] += (3.0*vmap_vx[0]*f[1]*vmap_vx1_sq+1.8*f[7]*vmap_vx1_cu+3.0*vmap_vx[1]*f[7]*vmap_vx0_sq+f[1]*vmap_vx0_cu)*volFact; 
  out[2] += (3.0*vmap_vx[0]*f[2]*vmap_vx1_sq+1.8*f[8]*vmap_vx1_cu+3.0*vmap_vx[1]*f[8]*vmap_vx0_sq+f[2]*vmap_vx0_cu)*volFact; 
  out[3] += (3.0*vmap_vx[0]*f[6]*vmap_vx1_sq+1.8*f[16]*vmap_vx1_cu+3.0*vmap_vx[1]*f[16]*vmap_vx0_sq+f[6]*vmap_vx0_cu)*volFact; 
  out[4] += (vmap_vy[1]*f[4]*vmap_vx1_sq+f[0]*vmap_vy[0]*vmap_vx1_sq+vmap_vy[1]*f[4]*vmap_vx0_sq+f[0]*vmap_vy[0]*vmap_vx0_sq+2.0*vmap_vx[0]*vmap_vx[1]*vmap_vy[1]*f[11]+2.0*vmap_vx[0]*vmap_vy[0]*vmap_vx[1]*f[3])*volFact; 
  out[5] += (vmap_vy[1]*f[9]*vmap_vx1_sq+vmap_vy[0]*f[1]*vmap_vx1_sq+vmap_vy[1]*f[9]*vmap_vx0_sq+vmap_vy[0]*f[1]*vmap_vx0_sq+2.0*vmap_vx[0]*vmap_vx[1]*vmap_vy[1]*f[18]+2.0*vmap_vx[0]*vmap_vy[0]*vmap_vx[1]*f[7])*volFact; 
  out[6] += (vmap_vy[1]*f[10]*vmap_vx1_sq+vmap_vy[0]*f[2]*vmap_vx1_sq+vmap_vy[1]*f[10]*vmap_vx0_sq+vmap_vy[0]*f[2]*vmap_vx0_sq+2.0*vmap_vx[0]*vmap_vx[1]*vmap_vy[1]*f[19]+2.0*vmap_vx[0]*vmap_vy[0]*vmap_vx[1]*f[8])*volFact; 
  out[7] += (vmap_vy[1]*f[17]*vmap_vx1_sq+vmap_vy[0]*f[6]*vmap_vx1_sq+vmap_vy[1]*f[17]*vmap_vx0_sq+vmap_vy[0]*f[6]*vmap_vx0_sq+2.0*vmap_vx[0]*vmap_vx[1]*vmap_vy[1]*f[26]+2.0*vmap_vx[0]*vmap_vy[0]*vmap_vx[1]*f[16])*volFact; 
  out[8] += (vmap_vz[1]*f[5]*vmap_vx1_sq+f[0]*vmap_vz[0]*vmap_vx1_sq+vmap_vz[1]*f[5]*vmap_vx0_sq+f[0]*vmap_vz[0]*vmap_vx0_sq+2.0*vmap_vx[0]*vmap_vx[1]*vmap_vz[1]*f[14]+2.0*vmap_vx[0]*vmap_vz[0]*vmap_vx[1]*f[3])*volFact; 
  out[9] += (vmap_vz[1]*f[12]*vmap_vx1_sq+vmap_vz[0]*f[1]*vmap_vx1_sq+vmap_vz[1]*f[12]*vmap_vx0_sq+vmap_vz[0]*f[1]*vmap_vx0_sq+2.0*vmap_vx[0]*vmap_vx[1]*vmap_vz[1]*f[21]+2.0*vmap_vx[0]*vmap_vz[0]*vmap_vx[1]*f[7])*volFact; 
  out[10] += (vmap_vz[1]*f[13]*vmap_vx1_sq+vmap_vz[0]*f[2]*vmap_vx1_sq+vmap_vz[1]*f[13]*vmap_vx0_sq+vmap_vz[0]*f[2]*vmap_vx0_sq+2.0*vmap_vx[0]*vmap_vx[1]*vmap_vz[1]*f[22]+2.0*vmap_vx[0]*vmap_vz[0]*vmap_vx[1]*f[8])*volFact; 
  out[11] += (vmap_vz[1]*f[20]*vmap_vx1_sq+vmap_vz[0]*f[6]*vmap_vx1_sq+vmap_vz[1]*f[20]*vmap_vx0_sq+vmap_vz[0]*f[6]*vmap_vx0_sq+2.0*vmap_vx[0]*vmap_vx[1]*vmap_vz[1]*f[27]+2.0*vmap_vx[0]*vmap_vz[0]*vmap_vx[1]*f[16])*volFact; 
  out[12] += (vmap_vx[1]*f[3]*vmap_vy1_sq+f[0]*vmap_vx[0]*vmap_vy1_sq+vmap_vx[1]*f[3]*vmap_vy0_sq+f[0]*vmap_vx[0]*vmap_vy0_sq+2.0*vmap_vy[0]*vmap_vx[1]*vmap_vy[1]*f[11]+2.0*vmap_vx[0]*vmap_vy[0]*vmap_vy[1]*f[4])*volFact; 
  out[13] += (vmap_vx[1]*f[7]*vmap_vy1_sq+vmap_vx[0]*f[1]*vmap_vy1_sq+vmap_vx[1]*f[7]*vmap_vy0_sq+vmap_vx[0]*f[1]*vmap_vy0_sq+2.0*vmap_vy[0]*vmap_vx[1]*vmap_vy[1]*f[18]+2.0*vmap_vx[0]*vmap_vy[0]*vmap_vy[1]*f[9])*volFact; 
  out[14] += (vmap_vx[1]*f[8]*vmap_vy1_sq+vmap_vx[0]*f[2]*vmap_vy1_sq+vmap_vx[1]*f[8]*vmap_vy0_sq+vmap_vx[0]*f[2]*vmap_vy0_sq+2.0*vmap_vy[0]*vmap_vx[1]*vmap_vy[1]*f[19]+2.0*vmap_vx[0]*vmap_vy[0]*vmap_vy[1]*f[10])*volFact; 
  out[15] += (vmap_vx[1]*f[16]*vmap_vy1_sq+vmap_vx[0]*f[6]*vmap_vy1_sq+vmap_vx[1]*f[16]*vmap_vy0_sq+vmap_vx[0]*f[6]*vmap_vy0_sq+2.0*vmap_vy[0]*vmap_vx[1]*vmap_vy[1]*f[26]+2.0*vmap_vx[0]*vmap_vy[0]*vmap_vy[1]*f[17])*volFact; 
  out[16] += (vmap_vx[1]*vmap_vy[1]*vmap_vz[1]*f[25]+vmap_vx[0]*vmap_vy[1]*vmap_vz[1]*f[15]+vmap_vy[0]*vmap_vx[1]*vmap_vz[1]*f[14]+vmap_vz[0]*vmap_vx[1]*vmap_vy[1]*f[11]+vmap_vx[0]*vmap_vy[0]*vmap_vz[1]*f[5]+vmap_vx[0]*vmap_vz[0]*vmap_vy[1]*f[4]+vmap_vy[0]*vmap_vz[0]*vmap_vx[1]*f[3]+f[0]*vmap_vx[0]*vmap_vy[0]*vmap_vz[0])*volFact; 
  out[17] += (vmap_vx[1]*vmap_vy[1]*vmap_vz[1]*f[29]+vmap_vx[0]*vmap_vy[1]*vmap_vz[1]*f[23]+vmap_vy[0]*vmap_vx[1]*vmap_vz[1]*f[21]+vmap_vz[0]*vmap_vx[1]*vmap_vy[1]*f[18]+vmap_vx[0]*vmap_vy[0]*vmap_vz[1]*f[12]+vmap_vx[0]*vmap_vz[0]*vmap_vy[1]*f[9]+vmap_vy[0]*vmap_vz[0]*vmap_vx[1]*f[7]+vmap_vx[0]*vmap_vy[0]*vmap_vz[0]*f[1])*volFact; 
  out[18] += (vmap_vx[1]*vmap_vy[1]*vmap_vz[1]*f[30]+vmap_vx[0]*vmap_vy[1]*vmap_vz[1]*f[24]+vmap_vy[0]*vmap_vx[1]*vmap_vz[1]*f[22]+vmap_vz[0]*vmap_vx[1]*vmap_vy[1]*f[19]+vmap_vx[0]*vmap_vy[0]*vmap_vz[1]*f[13]+vmap_vx[0]*vmap_vz[0]*vmap_vy[1]*f[10]+vmap_vy[0]*vmap_vz[0]*vmap_vx[1]*f[8]+vmap_vx[0]*vmap_vy[0]*vmap_vz[0]*f[2])*volFact; 
  out[19] += (vmap_vx[1]*vmap_vy[1]*vmap_vz[1]*f[31]+vmap_vx[0]*vmap_vy[1]*vmap_vz[1]*f[28]+vmap_vy[0]*vmap_vx[1]*vmap_vz[1]*f[27]+vmap_vz[0]*vmap_vx[1]*vmap_vy[1]*f[26]+vmap_vx[0]*vmap_vy[0]*vmap_vz[1]*f[20]+vmap_vx[0]*vmap_vz[0]*vmap_vy[1]*f[17]+vmap_vy[0]*vmap_vz[0]*vmap_vx[1]*f[16]+vmap_vx[0]*vmap_vy[0]*vmap_vz[0]*f[6])*volFact; 
  out[20] += (vmap_vx[1]*f[3]*vmap_vz1_sq+f[0]*vmap_vx[0]*vmap_vz1_sq+vmap_vx[1]*f[3]*vmap_vz0_sq+f[0]*vmap_vx[0]*vmap_vz0_sq+2.0*vmap_vz[0]*vmap_vx[1]*vmap_vz[1]*f[14]+2.0*vmap_vx[0]*vmap_vz[0]*vmap_vz[1]*f[5])*volFact; 
  out[21] += (vmap_vx[1]*f[7]*vmap_vz1_sq+vmap_vx[0]*f[1]*vmap_vz1_sq+vmap_vx[1]*f[7]*vmap_vz0_sq+vmap_vx[0]*f[1]*vmap_vz0_sq+2.0*vmap_vz[0]*vmap_vx[1]*vmap_vz[1]*f[21]+2.0*vmap_vx[0]*vmap_vz[0]*vmap_vz[1]*f[12])*volFact; 
  out[22] += (vmap_vx[1]*f[8]*vmap_vz1_sq+vmap_vx[0]*f[2]*vmap_vz1_sq+vmap_vx[1]*f[8]*vmap_vz0_sq+vmap_vx[0]*f[2]*vmap_vz0_sq+2.0*vmap_vz[0]*vmap_vx[1]*vmap_vz[1]*f[22]+2.0*vmap_vx[0]*vmap_vz[0]*vmap_vz[1]*f[13])*volFact; 
  out[23] += (vmap_vx[1]*f[16]*vmap_vz1_sq+vmap_vx[0]*f[6]*vmap_vz1_sq+vmap_vx[1]*f[16]*vmap_vz0_sq+vmap_vx[0]*f[6]*vmap_vz0_sq+2.0*vmap_vz[0]*vmap_vx[1]*vmap_vz[1]*f[27]+2.0*vmap_vx[0]*vmap_vz[0]*vmap_vz[1]*f[20])*volFact; 
  out[24] += (3.0*f[0]*vmap_vy[0]*vmap_vy1_sq+1.8*f[4]*vmap_vy1_cu+3.0*vmap_vy[1]*f[4]*vmap_vy0_sq+f[0]*vmap_vy0_cu)*volFact; 
  out[25] += (3.0*vmap_vy[0]*f[1]*vmap_vy1_sq+1.8*f[9]*vmap_vy1_cu+3.0*vmap_vy[1]*f[9]*vmap_vy0_sq+f[1]*vmap_vy0_cu)*volFact; 
  out[26] += (3.0*vmap_vy[0]*f[2]*vmap_vy1_sq+1.8*f[10]*vmap_vy1_cu+3.0*vmap_vy[1]*f[10]*vmap_vy0_sq+f[2]*vmap_vy0_cu)*volFact; 
  out[27] += (3.0*vmap_vy[0]*f[6]*vmap_vy1_sq+1.8*f[17]*vmap_vy1_cu+3.0*vmap_vy[1]*f[17]*vmap_vy0_sq+f[6]*vmap_vy0_cu)*volFact; 
  out[28] += (vmap_vz[1]*f[5]*vmap_vy1_sq+f[0]*vmap_vz[0]*vmap_vy1_sq+vmap_vz[1]*f[5]*vmap_vy0_sq+f[0]*vmap_vz[0]*vmap_vy0_sq+2.0*vmap_vy[0]*vmap_vy[1]*vmap_vz[1]*f[15]+2.0*vmap_vy[0]*vmap_vz[0]*vmap_vy[1]*f[4])*volFact; 
  out[29] += (vmap_vz[1]*f[12]*vmap_vy1_sq+vmap_vz[0]*f[1]*vmap_vy1_sq+vmap_vz[1]*f[12]*vmap_vy0_sq+vmap_vz[0]*f[1]*vmap_vy0_sq+2.0*vmap_vy[0]*vmap_vy[1]*vmap_vz[1]*f[23]+2.0*vmap_vy[0]*vmap_vz[0]*vmap_vy[1]*f[9])*volFact; 
  out[30] += (vmap_vz[1]*f[13]*vmap_vy1_sq+vmap_vz[0]*f[2]*vmap_vy1_sq+vmap_vz[1]*f[13]*vmap_vy0_sq+vmap_vz[0]*f[2]*vmap_vy0_sq+2.0*vmap_vy[0]*vmap_vy[1]*vmap_vz[1]*f[24]+2.0*vmap_vy[0]*vmap_vz[0]*vmap_vy[1]*f[10])*volFact; 
  out[31] += (vmap_vz[1]*f[20]*vmap_vy1_sq+vmap_vz[0]*f[6]*vmap_vy1_sq+vmap_vz[1]*f[20]*vmap_vy0_sq+vmap_vz[0]*f[6]*vmap_vy0_sq+2.0*vmap_vy[0]*vmap_vy[1]*vmap_vz[1]*f[28]+2.0*vmap_vy[0]*vmap_vz[0]*vmap_vy[1]*f[17])*volFact; 
  out[32] += (vmap_vy[1]*f[4]*vmap_vz1_sq+f[0]*vmap_vy[0]*vmap_vz1_sq+vmap_vy[1]*f[4]*vmap_vz0_sq+f[0]*vmap_vy[0]*vmap_vz0_sq+2.0*vmap_vz[0]*vmap_vy[1]*vmap_vz[1]*f[15]+2.0*vmap_vy[0]*vmap_vz[0]*vmap_vz[1]*f[5])*volFact; 
  out[33] += (vmap_vy[1]*f[9]*vmap_vz1_sq+vmap_vy[0]*f[1]*vmap_vz1_sq+vmap_vy[1]*f[9]*vmap_vz0_sq+vmap_vy[0]*f[1]*vmap_vz0_sq+2.0*vmap_vz[0]*vmap_vy[1]*vmap_vz[1]*f[23]+2.0*vmap_vy[0]*vmap_vz[0]*vmap_vz[1]*f[12])*volFact; 
  out[34] += (vmap_vy[1]*f[10]*vmap_vz1_sq+vmap_vy[0]*f[2]*vmap_vz1_sq+vmap_vy[1]*f[10]*vmap_vz0_sq+vmap_vy[0]*f[2]*vmap_vz0_sq+2.0*vmap_vz[0]*vmap_vy[1]*vmap_vz[1]*f[24]+2.0*vmap_vy[0]*vmap_vz[0]*vmap_vz[1]*f[13])*volFact; 
  out[35] += (vmap_vy[1]*f[17]*vmap_vz1_sq+vmap_vy[0]*f[6]*vmap_vz1_sq+vmap_vy[1]*f[17]*vmap_vz0_sq+vmap_vy[0]*f[6]*vmap_vz0_sq+2.0*vmap_vz[0]*vmap_vy[1]*vmap_vz[1]*f[28]+2.0*vmap_vy[0]*vmap_vz[0]*vmap_vz[1]*f[20])*volFact; 
  out[36] += (3.0*f[0]*vmap_vz[0]*vmap_vz1_sq+1.8*f[5]*vmap_vz1_cu+3.0*vmap_vz[1]*f[5]*vmap_vz0_sq+f[0]*vmap_vz0_cu)*volFact; 
  out[37] += (3.0*vmap_vz[0]*f[1]*vmap_vz1_sq+1.8*f[12]*vmap_vz1_cu+3.0*vmap_vz[1]*f[12]*vmap_vz0_sq+f[1]*vmap_vz0_cu)*volFact; 
  out[38] += (3.0*vmap_vz[0]*f[2]*vmap_vz1_sq+1.8*f[13]*vmap_vz1_cu+3.0*vmap_vz[1]*f[13]*vmap_vz0_sq+f[2]*vmap_vz0_cu)*volFact; 
  out[39] += (3.0*vmap_vz[0]*f[6]*vmap_vz1_sq+1.8*f[20]*vmap_vz1_cu+3.0*vmap_vz[1]*f[20]*vmap_vz0_sq+f[6]*vmap_vz0_cu)*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_M1i_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[3]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  double dv12 = 2.0/dxv[4]; 
  const double *jacob_vz = &jacob_vel[4]; 
  const double jacob_vz_inv = 1.0/jacob_vz[0]; 
  out[0] += (1.7320508075688772*hamil[7]*f[15]+1.7320508075688772*f[5]*hamil[5]+1.7320508075688772*f[4]*hamil[4]+1.7320508075688772*f[0]*hamil[1])*dv10*jacob_vx_inv*volFact; 
  out[1] += (1.7320508075688772*hamil[7]*f[23]+1.7320508075688772*hamil[5]*f[12]+1.7320508075688772*hamil[4]*f[9]+1.7320508075688772*f[1]*hamil[1])*dv10*jacob_vx_inv*volFact; 
  out[2] += (1.7320508075688772*hamil[7]*f[24]+1.7320508075688772*hamil[5]*f[13]+1.7320508075688772*hamil[4]*f[10]+1.7320508075688772*hamil[1]*f[2])*dv10*jacob_vx_inv*volFact; 
  out[3] += (1.7320508075688772*hamil[7]*f[28]+1.7320508075688772*hamil[5]*f[20]+1.7320508075688772*hamil[4]*f[17]+1.7320508075688772*hamil[1]*f[6])*dv10*jacob_vx_inv*volFact; 
  out[4] += (1.7320508075688772*hamil[7]*f[14]+1.7320508075688772*f[5]*hamil[6]+1.7320508075688772*f[3]*hamil[4]+1.7320508075688772*f[0]*hamil[2])*dv11*jacob_vy_inv*volFact; 
  out[5] += (1.7320508075688772*hamil[7]*f[21]+1.7320508075688772*hamil[6]*f[12]+1.7320508075688772*hamil[4]*f[7]+1.7320508075688772*f[1]*hamil[2])*dv11*jacob_vy_inv*volFact; 
  out[6] += (1.7320508075688772*hamil[7]*f[22]+1.7320508075688772*hamil[6]*f[13]+1.7320508075688772*hamil[4]*f[8]+1.7320508075688772*f[2]*hamil[2])*dv11*jacob_vy_inv*volFact; 
  out[7] += (1.7320508075688772*hamil[7]*f[27]+1.7320508075688772*hamil[6]*f[20]+1.7320508075688772*hamil[4]*f[16]+1.7320508075688772*hamil[2]*f[6])*dv11*jacob_vy_inv*volFact; 
  out[8] += (1.7320508075688772*hamil[7]*f[11]+1.7320508075688772*f[4]*hamil[6]+1.7320508075688772*f[3]*hamil[5]+1.7320508075688772*f[0]*hamil[3])*dv12*jacob_vz_inv*volFact; 
  out[9] += (1.7320508075688772*hamil[7]*f[18]+1.7320508075688772*hamil[6]*f[9]+1.7320508075688772*hamil[5]*f[7]+1.7320508075688772*f[1]*hamil[3])*dv12*jacob_vz_inv*volFact; 
  out[10] += (1.7320508075688772*hamil[7]*f[19]+1.7320508075688772*hamil[6]*f[10]+1.7320508075688772*hamil[5]*f[8]+1.7320508075688772*f[2]*hamil[3])*dv12*jacob_vz_inv*volFact; 
  out[11] += (1.7320508075688772*hamil[7]*f[26]+1.7320508075688772*hamil[6]*f[17]+1.7320508075688772*hamil[5]*f[16]+1.7320508075688772*hamil[3]*f[6])*dv12*jacob_vz_inv*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_M2_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  out[0] += (hamil[7]*f[25]+hamil[6]*f[15]+hamil[5]*f[14]+hamil[4]*f[11]+hamil[3]*f[5]+hamil[2]*f[4]+hamil[1]*f[3]+f[0]*hamil[0])*volFact; 
  out[1] += (hamil[7]*f[29]+hamil[6]*f[23]+hamil[5]*f[21]+hamil[4]*f[18]+hamil[3]*f[12]+hamil[2]*f[9]+hamil[1]*f[7]+hamil[0]*f[1])*volFact; 
  out[2] += (hamil[7]*f[30]+hamil[6]*f[24]+hamil[5]*f[22]+hamil[4]*f[19]+hamil[3]*f[13]+hamil[2]*f[10]+hamil[1]*f[8]+hamil[0]*f[2])*volFact; 
  out[3] += (hamil[7]*f[31]+hamil[6]*f[28]+hamil[5]*f[27]+hamil[4]*f[26]+hamil[3]*f[20]+hamil[2]*f[17]+hamil[1]*f[16]+hamil[0]*f[6])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_M3i_2x3v_ser_p1(const double *w, const double *dxv, const int *idx,
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[3]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  double dv12 = 2.0/dxv[4]; 
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

  out[0] += (0.3535533905932737*dH_dpx[0]*hamil[7]*f[25]+0.3535533905932737*hamil[1]*dH_dpx[6]*f[25]+0.3535533905932737*dH_dpx[2]*hamil[5]*f[25]+0.3535533905932737*dH_dpx[3]*hamil[4]*f[25]+0.3535533905932737*dH_dpx[0]*hamil[6]*f[15]+0.3535533905932737*hamil[0]*dH_dpx[6]*f[15]+0.3535533905932737*dH_dpx[2]*hamil[3]*f[15]+0.3535533905932737*hamil[2]*dH_dpx[3]*f[15]+0.3535533905932737*dH_dpx[2]*hamil[7]*f[14]+0.3535533905932737*hamil[4]*dH_dpx[6]*f[14]+0.3535533905932737*dH_dpx[0]*hamil[5]*f[14]+0.3535533905932737*hamil[1]*dH_dpx[3]*f[14]+0.3535533905932737*dH_dpx[3]*hamil[7]*f[11]+0.3535533905932737*hamil[5]*dH_dpx[6]*f[11]+0.3535533905932737*dH_dpx[0]*hamil[4]*f[11]+0.3535533905932737*hamil[1]*dH_dpx[2]*f[11]+0.3535533905932737*f[3]*dH_dpx[6]*hamil[7]+0.3535533905932737*f[0]*dH_dpx[6]*hamil[6]+0.3535533905932737*dH_dpx[2]*f[5]*hamil[6]+0.3535533905932737*dH_dpx[3]*f[4]*hamil[6]+0.3535533905932737*hamil[2]*f[5]*dH_dpx[6]+0.3535533905932737*hamil[3]*f[4]*dH_dpx[6]+0.3535533905932737*dH_dpx[3]*f[3]*hamil[5]+0.3535533905932737*dH_dpx[0]*hamil[3]*f[5]+0.3535533905932737*hamil[0]*dH_dpx[3]*f[5]+0.3535533905932737*dH_dpx[2]*f[3]*hamil[4]+0.3535533905932737*dH_dpx[0]*hamil[2]*f[4]+0.3535533905932737*hamil[0]*dH_dpx[2]*f[4]+0.3535533905932737*f[0]*dH_dpx[3]*hamil[3]+0.3535533905932737*dH_dpx[0]*hamil[1]*f[3]+0.3535533905932737*f[0]*dH_dpx[2]*hamil[2]+0.3535533905932737*dH_dpx[0]*f[0]*hamil[0])*volFact; 
  out[1] += (0.3535533905932737*dH_dpx[0]*hamil[7]*f[29]+0.3535533905932737*hamil[1]*dH_dpx[6]*f[29]+0.3535533905932737*dH_dpx[2]*hamil[5]*f[29]+0.3535533905932737*dH_dpx[3]*hamil[4]*f[29]+0.3535533905932737*dH_dpx[0]*hamil[6]*f[23]+0.3535533905932737*hamil[0]*dH_dpx[6]*f[23]+0.3535533905932737*dH_dpx[2]*hamil[3]*f[23]+0.3535533905932737*hamil[2]*dH_dpx[3]*f[23]+0.3535533905932737*dH_dpx[2]*hamil[7]*f[21]+0.3535533905932737*hamil[4]*dH_dpx[6]*f[21]+0.3535533905932737*dH_dpx[0]*hamil[5]*f[21]+0.3535533905932737*hamil[1]*dH_dpx[3]*f[21]+0.3535533905932737*dH_dpx[3]*hamil[7]*f[18]+0.3535533905932737*hamil[5]*dH_dpx[6]*f[18]+0.3535533905932737*dH_dpx[0]*hamil[4]*f[18]+0.3535533905932737*hamil[1]*dH_dpx[2]*f[18]+0.3535533905932737*dH_dpx[2]*hamil[6]*f[12]+0.3535533905932737*hamil[2]*dH_dpx[6]*f[12]+0.3535533905932737*dH_dpx[0]*hamil[3]*f[12]+0.3535533905932737*hamil[0]*dH_dpx[3]*f[12]+0.3535533905932737*dH_dpx[3]*hamil[6]*f[9]+0.3535533905932737*hamil[3]*dH_dpx[6]*f[9]+0.3535533905932737*dH_dpx[0]*hamil[2]*f[9]+0.3535533905932737*hamil[0]*dH_dpx[2]*f[9]+0.3535533905932737*dH_dpx[6]*f[7]*hamil[7]+0.3535533905932737*dH_dpx[3]*hamil[5]*f[7]+0.3535533905932737*dH_dpx[2]*hamil[4]*f[7]+0.3535533905932737*dH_dpx[0]*hamil[1]*f[7]+0.3535533905932737*f[1]*dH_dpx[6]*hamil[6]+0.3535533905932737*f[1]*dH_dpx[3]*hamil[3]+0.3535533905932737*f[1]*dH_dpx[2]*hamil[2]+0.3535533905932737*dH_dpx[0]*hamil[0]*f[1])*volFact; 
  out[2] += (0.3535533905932737*dH_dpx[0]*hamil[7]*f[30]+0.3535533905932737*hamil[1]*dH_dpx[6]*f[30]+0.3535533905932737*dH_dpx[2]*hamil[5]*f[30]+0.3535533905932737*dH_dpx[3]*hamil[4]*f[30]+0.3535533905932737*dH_dpx[0]*hamil[6]*f[24]+0.3535533905932737*hamil[0]*dH_dpx[6]*f[24]+0.3535533905932737*dH_dpx[2]*hamil[3]*f[24]+0.3535533905932737*hamil[2]*dH_dpx[3]*f[24]+0.3535533905932737*dH_dpx[2]*hamil[7]*f[22]+0.3535533905932737*hamil[4]*dH_dpx[6]*f[22]+0.3535533905932737*dH_dpx[0]*hamil[5]*f[22]+0.3535533905932737*hamil[1]*dH_dpx[3]*f[22]+0.3535533905932737*dH_dpx[3]*hamil[7]*f[19]+0.3535533905932737*hamil[5]*dH_dpx[6]*f[19]+0.3535533905932737*dH_dpx[0]*hamil[4]*f[19]+0.3535533905932737*hamil[1]*dH_dpx[2]*f[19]+0.3535533905932737*dH_dpx[2]*hamil[6]*f[13]+0.3535533905932737*hamil[2]*dH_dpx[6]*f[13]+0.3535533905932737*dH_dpx[0]*hamil[3]*f[13]+0.3535533905932737*hamil[0]*dH_dpx[3]*f[13]+0.3535533905932737*dH_dpx[3]*hamil[6]*f[10]+0.3535533905932737*hamil[3]*dH_dpx[6]*f[10]+0.3535533905932737*dH_dpx[0]*hamil[2]*f[10]+0.3535533905932737*hamil[0]*dH_dpx[2]*f[10]+0.3535533905932737*dH_dpx[6]*hamil[7]*f[8]+0.3535533905932737*dH_dpx[3]*hamil[5]*f[8]+0.3535533905932737*dH_dpx[2]*hamil[4]*f[8]+0.3535533905932737*dH_dpx[0]*hamil[1]*f[8]+0.3535533905932737*f[2]*dH_dpx[6]*hamil[6]+0.3535533905932737*f[2]*dH_dpx[3]*hamil[3]+0.3535533905932737*dH_dpx[2]*f[2]*hamil[2]+0.3535533905932737*dH_dpx[0]*hamil[0]*f[2])*volFact; 
  out[3] += (0.3535533905932737*dH_dpx[0]*hamil[7]*f[31]+0.3535533905932737*hamil[1]*dH_dpx[6]*f[31]+0.3535533905932737*dH_dpx[2]*hamil[5]*f[31]+0.3535533905932737*dH_dpx[3]*hamil[4]*f[31]+0.3535533905932737*dH_dpx[0]*hamil[6]*f[28]+0.3535533905932737*hamil[0]*dH_dpx[6]*f[28]+0.3535533905932737*dH_dpx[2]*hamil[3]*f[28]+0.3535533905932737*hamil[2]*dH_dpx[3]*f[28]+0.3535533905932737*dH_dpx[2]*hamil[7]*f[27]+0.3535533905932737*hamil[4]*dH_dpx[6]*f[27]+0.3535533905932737*dH_dpx[0]*hamil[5]*f[27]+0.3535533905932737*hamil[1]*dH_dpx[3]*f[27]+0.3535533905932737*dH_dpx[3]*hamil[7]*f[26]+0.3535533905932737*hamil[5]*dH_dpx[6]*f[26]+0.3535533905932737*dH_dpx[0]*hamil[4]*f[26]+0.3535533905932737*hamil[1]*dH_dpx[2]*f[26]+0.3535533905932737*dH_dpx[2]*hamil[6]*f[20]+0.3535533905932737*hamil[2]*dH_dpx[6]*f[20]+0.3535533905932737*dH_dpx[0]*hamil[3]*f[20]+0.3535533905932737*hamil[0]*dH_dpx[3]*f[20]+0.3535533905932737*dH_dpx[3]*hamil[6]*f[17]+0.3535533905932737*hamil[3]*dH_dpx[6]*f[17]+0.3535533905932737*dH_dpx[0]*hamil[2]*f[17]+0.3535533905932737*hamil[0]*dH_dpx[2]*f[17]+0.3535533905932737*dH_dpx[6]*hamil[7]*f[16]+0.3535533905932737*dH_dpx[3]*hamil[5]*f[16]+0.3535533905932737*dH_dpx[2]*hamil[4]*f[16]+0.3535533905932737*dH_dpx[0]*hamil[1]*f[16]+0.3535533905932737*dH_dpx[6]*f[6]*hamil[6]+0.3535533905932737*dH_dpx[3]*hamil[3]*f[6]+0.3535533905932737*dH_dpx[2]*hamil[2]*f[6]+0.3535533905932737*dH_dpx[0]*hamil[0]*f[6])*volFact; 
  out[4] += (0.3535533905932737*dH_dpy[0]*hamil[7]*f[25]+0.3535533905932737*dH_dpy[1]*hamil[6]*f[25]+0.3535533905932737*hamil[2]*dH_dpy[5]*f[25]+0.3535533905932737*dH_dpy[3]*hamil[4]*f[25]+0.3535533905932737*dH_dpy[1]*hamil[7]*f[15]+0.3535533905932737*dH_dpy[0]*hamil[6]*f[15]+0.3535533905932737*hamil[4]*dH_dpy[5]*f[15]+0.3535533905932737*hamil[2]*dH_dpy[3]*f[15]+0.3535533905932737*dH_dpy[0]*hamil[5]*f[14]+0.3535533905932737*hamil[0]*dH_dpy[5]*f[14]+0.3535533905932737*dH_dpy[1]*hamil[3]*f[14]+0.3535533905932737*hamil[1]*dH_dpy[3]*f[14]+0.3535533905932737*dH_dpy[3]*hamil[7]*f[11]+0.3535533905932737*dH_dpy[5]*hamil[6]*f[11]+0.3535533905932737*dH_dpy[0]*hamil[4]*f[11]+0.3535533905932737*dH_dpy[1]*hamil[2]*f[11]+0.3535533905932737*f[4]*dH_dpy[5]*hamil[7]+0.3535533905932737*dH_dpy[3]*f[4]*hamil[6]+0.3535533905932737*dH_dpy[1]*f[5]*hamil[5]+0.3535533905932737*f[0]*dH_dpy[5]*hamil[5]+0.3535533905932737*dH_dpy[3]*f[3]*hamil[5]+0.3535533905932737*hamil[1]*dH_dpy[5]*f[5]+0.3535533905932737*dH_dpy[0]*hamil[3]*f[5]+0.3535533905932737*hamil[0]*dH_dpy[3]*f[5]+0.3535533905932737*f[3]*hamil[3]*dH_dpy[5]+0.3535533905932737*dH_dpy[1]*f[4]*hamil[4]+0.3535533905932737*dH_dpy[0]*hamil[2]*f[4]+0.3535533905932737*f[0]*dH_dpy[3]*hamil[3]+0.3535533905932737*dH_dpy[0]*hamil[1]*f[3]+0.3535533905932737*hamil[0]*dH_dpy[1]*f[3]+0.3535533905932737*f[0]*dH_dpy[1]*hamil[1]+0.3535533905932737*dH_dpy[0]*f[0]*hamil[0])*volFact; 
  out[5] += (0.3535533905932737*dH_dpy[0]*hamil[7]*f[29]+0.3535533905932737*dH_dpy[1]*hamil[6]*f[29]+0.3535533905932737*hamil[2]*dH_dpy[5]*f[29]+0.3535533905932737*dH_dpy[3]*hamil[4]*f[29]+0.3535533905932737*dH_dpy[1]*hamil[7]*f[23]+0.3535533905932737*dH_dpy[0]*hamil[6]*f[23]+0.3535533905932737*hamil[4]*dH_dpy[5]*f[23]+0.3535533905932737*hamil[2]*dH_dpy[3]*f[23]+0.3535533905932737*dH_dpy[0]*hamil[5]*f[21]+0.3535533905932737*hamil[0]*dH_dpy[5]*f[21]+0.3535533905932737*dH_dpy[1]*hamil[3]*f[21]+0.3535533905932737*hamil[1]*dH_dpy[3]*f[21]+0.3535533905932737*dH_dpy[3]*hamil[7]*f[18]+0.3535533905932737*dH_dpy[5]*hamil[6]*f[18]+0.3535533905932737*dH_dpy[0]*hamil[4]*f[18]+0.3535533905932737*dH_dpy[1]*hamil[2]*f[18]+0.3535533905932737*dH_dpy[1]*hamil[5]*f[12]+0.3535533905932737*hamil[1]*dH_dpy[5]*f[12]+0.3535533905932737*dH_dpy[0]*hamil[3]*f[12]+0.3535533905932737*hamil[0]*dH_dpy[3]*f[12]+0.3535533905932737*dH_dpy[5]*hamil[7]*f[9]+0.3535533905932737*dH_dpy[3]*hamil[6]*f[9]+0.3535533905932737*dH_dpy[1]*hamil[4]*f[9]+0.3535533905932737*dH_dpy[0]*hamil[2]*f[9]+0.3535533905932737*dH_dpy[3]*hamil[5]*f[7]+0.3535533905932737*hamil[3]*dH_dpy[5]*f[7]+0.3535533905932737*dH_dpy[0]*hamil[1]*f[7]+0.3535533905932737*hamil[0]*dH_dpy[1]*f[7]+0.3535533905932737*f[1]*dH_dpy[5]*hamil[5]+0.3535533905932737*f[1]*dH_dpy[3]*hamil[3]+0.3535533905932737*dH_dpy[1]*f[1]*hamil[1]+0.3535533905932737*dH_dpy[0]*hamil[0]*f[1])*volFact; 
  out[6] += (0.3535533905932737*dH_dpy[0]*hamil[7]*f[30]+0.3535533905932737*dH_dpy[1]*hamil[6]*f[30]+0.3535533905932737*hamil[2]*dH_dpy[5]*f[30]+0.3535533905932737*dH_dpy[3]*hamil[4]*f[30]+0.3535533905932737*dH_dpy[1]*hamil[7]*f[24]+0.3535533905932737*dH_dpy[0]*hamil[6]*f[24]+0.3535533905932737*hamil[4]*dH_dpy[5]*f[24]+0.3535533905932737*hamil[2]*dH_dpy[3]*f[24]+0.3535533905932737*dH_dpy[0]*hamil[5]*f[22]+0.3535533905932737*hamil[0]*dH_dpy[5]*f[22]+0.3535533905932737*dH_dpy[1]*hamil[3]*f[22]+0.3535533905932737*hamil[1]*dH_dpy[3]*f[22]+0.3535533905932737*dH_dpy[3]*hamil[7]*f[19]+0.3535533905932737*dH_dpy[5]*hamil[6]*f[19]+0.3535533905932737*dH_dpy[0]*hamil[4]*f[19]+0.3535533905932737*dH_dpy[1]*hamil[2]*f[19]+0.3535533905932737*dH_dpy[1]*hamil[5]*f[13]+0.3535533905932737*hamil[1]*dH_dpy[5]*f[13]+0.3535533905932737*dH_dpy[0]*hamil[3]*f[13]+0.3535533905932737*hamil[0]*dH_dpy[3]*f[13]+0.3535533905932737*dH_dpy[5]*hamil[7]*f[10]+0.3535533905932737*dH_dpy[3]*hamil[6]*f[10]+0.3535533905932737*dH_dpy[1]*hamil[4]*f[10]+0.3535533905932737*dH_dpy[0]*hamil[2]*f[10]+0.3535533905932737*dH_dpy[3]*hamil[5]*f[8]+0.3535533905932737*hamil[3]*dH_dpy[5]*f[8]+0.3535533905932737*dH_dpy[0]*hamil[1]*f[8]+0.3535533905932737*hamil[0]*dH_dpy[1]*f[8]+0.3535533905932737*f[2]*dH_dpy[5]*hamil[5]+0.3535533905932737*f[2]*dH_dpy[3]*hamil[3]+0.3535533905932737*dH_dpy[1]*hamil[1]*f[2]+0.3535533905932737*dH_dpy[0]*hamil[0]*f[2])*volFact; 
  out[7] += (0.3535533905932737*dH_dpy[0]*hamil[7]*f[31]+0.3535533905932737*dH_dpy[1]*hamil[6]*f[31]+0.3535533905932737*hamil[2]*dH_dpy[5]*f[31]+0.3535533905932737*dH_dpy[3]*hamil[4]*f[31]+0.3535533905932737*dH_dpy[1]*hamil[7]*f[28]+0.3535533905932737*dH_dpy[0]*hamil[6]*f[28]+0.3535533905932737*hamil[4]*dH_dpy[5]*f[28]+0.3535533905932737*hamil[2]*dH_dpy[3]*f[28]+0.3535533905932737*dH_dpy[0]*hamil[5]*f[27]+0.3535533905932737*hamil[0]*dH_dpy[5]*f[27]+0.3535533905932737*dH_dpy[1]*hamil[3]*f[27]+0.3535533905932737*hamil[1]*dH_dpy[3]*f[27]+0.3535533905932737*dH_dpy[3]*hamil[7]*f[26]+0.3535533905932737*dH_dpy[5]*hamil[6]*f[26]+0.3535533905932737*dH_dpy[0]*hamil[4]*f[26]+0.3535533905932737*dH_dpy[1]*hamil[2]*f[26]+0.3535533905932737*dH_dpy[1]*hamil[5]*f[20]+0.3535533905932737*hamil[1]*dH_dpy[5]*f[20]+0.3535533905932737*dH_dpy[0]*hamil[3]*f[20]+0.3535533905932737*hamil[0]*dH_dpy[3]*f[20]+0.3535533905932737*dH_dpy[5]*hamil[7]*f[17]+0.3535533905932737*dH_dpy[3]*hamil[6]*f[17]+0.3535533905932737*dH_dpy[1]*hamil[4]*f[17]+0.3535533905932737*dH_dpy[0]*hamil[2]*f[17]+0.3535533905932737*dH_dpy[3]*hamil[5]*f[16]+0.3535533905932737*hamil[3]*dH_dpy[5]*f[16]+0.3535533905932737*dH_dpy[0]*hamil[1]*f[16]+0.3535533905932737*hamil[0]*dH_dpy[1]*f[16]+0.3535533905932737*dH_dpy[5]*hamil[5]*f[6]+0.3535533905932737*dH_dpy[3]*hamil[3]*f[6]+0.3535533905932737*dH_dpy[1]*hamil[1]*f[6]+0.3535533905932737*dH_dpy[0]*hamil[0]*f[6])*volFact; 
  out[8] += (0.3535533905932737*dH_dpz[0]*hamil[7]*f[25]+0.3535533905932737*dH_dpz[1]*hamil[6]*f[25]+0.3535533905932737*dH_dpz[2]*hamil[5]*f[25]+0.3535533905932737*hamil[3]*dH_dpz[4]*f[25]+0.3535533905932737*dH_dpz[1]*hamil[7]*f[15]+0.3535533905932737*dH_dpz[0]*hamil[6]*f[15]+0.3535533905932737*dH_dpz[4]*hamil[5]*f[15]+0.3535533905932737*dH_dpz[2]*hamil[3]*f[15]+0.3535533905932737*dH_dpz[2]*hamil[7]*f[14]+0.3535533905932737*dH_dpz[4]*hamil[6]*f[14]+0.3535533905932737*dH_dpz[0]*hamil[5]*f[14]+0.3535533905932737*dH_dpz[1]*hamil[3]*f[14]+0.3535533905932737*dH_dpz[0]*hamil[4]*f[11]+0.3535533905932737*hamil[0]*dH_dpz[4]*f[11]+0.3535533905932737*dH_dpz[1]*hamil[2]*f[11]+0.3535533905932737*hamil[1]*dH_dpz[2]*f[11]+0.3535533905932737*dH_dpz[4]*f[5]*hamil[7]+0.3535533905932737*dH_dpz[2]*f[5]*hamil[6]+0.3535533905932737*dH_dpz[1]*f[5]*hamil[5]+0.3535533905932737*dH_dpz[0]*hamil[3]*f[5]+0.3535533905932737*dH_dpz[1]*f[4]*hamil[4]+0.3535533905932737*f[0]*dH_dpz[4]*hamil[4]+0.3535533905932737*dH_dpz[2]*f[3]*hamil[4]+0.3535533905932737*hamil[1]*dH_dpz[4]*f[4]+0.3535533905932737*dH_dpz[0]*hamil[2]*f[4]+0.3535533905932737*hamil[0]*dH_dpz[2]*f[4]+0.3535533905932737*hamil[2]*f[3]*dH_dpz[4]+0.3535533905932737*dH_dpz[0]*hamil[1]*f[3]+0.3535533905932737*hamil[0]*dH_dpz[1]*f[3]+0.3535533905932737*f[0]*dH_dpz[2]*hamil[2]+0.3535533905932737*f[0]*dH_dpz[1]*hamil[1]+0.3535533905932737*dH_dpz[0]*f[0]*hamil[0])*volFact; 
  out[9] += (0.3535533905932737*dH_dpz[0]*hamil[7]*f[29]+0.3535533905932737*dH_dpz[1]*hamil[6]*f[29]+0.3535533905932737*dH_dpz[2]*hamil[5]*f[29]+0.3535533905932737*hamil[3]*dH_dpz[4]*f[29]+0.3535533905932737*dH_dpz[1]*hamil[7]*f[23]+0.3535533905932737*dH_dpz[0]*hamil[6]*f[23]+0.3535533905932737*dH_dpz[4]*hamil[5]*f[23]+0.3535533905932737*dH_dpz[2]*hamil[3]*f[23]+0.3535533905932737*dH_dpz[2]*hamil[7]*f[21]+0.3535533905932737*dH_dpz[4]*hamil[6]*f[21]+0.3535533905932737*dH_dpz[0]*hamil[5]*f[21]+0.3535533905932737*dH_dpz[1]*hamil[3]*f[21]+0.3535533905932737*dH_dpz[0]*hamil[4]*f[18]+0.3535533905932737*hamil[0]*dH_dpz[4]*f[18]+0.3535533905932737*dH_dpz[1]*hamil[2]*f[18]+0.3535533905932737*hamil[1]*dH_dpz[2]*f[18]+0.3535533905932737*dH_dpz[4]*hamil[7]*f[12]+0.3535533905932737*dH_dpz[2]*hamil[6]*f[12]+0.3535533905932737*dH_dpz[1]*hamil[5]*f[12]+0.3535533905932737*dH_dpz[0]*hamil[3]*f[12]+0.3535533905932737*dH_dpz[1]*hamil[4]*f[9]+0.3535533905932737*hamil[1]*dH_dpz[4]*f[9]+0.3535533905932737*dH_dpz[0]*hamil[2]*f[9]+0.3535533905932737*hamil[0]*dH_dpz[2]*f[9]+0.3535533905932737*dH_dpz[2]*hamil[4]*f[7]+0.3535533905932737*hamil[2]*dH_dpz[4]*f[7]+0.3535533905932737*dH_dpz[0]*hamil[1]*f[7]+0.3535533905932737*hamil[0]*dH_dpz[1]*f[7]+0.3535533905932737*f[1]*dH_dpz[4]*hamil[4]+0.3535533905932737*f[1]*dH_dpz[2]*hamil[2]+0.3535533905932737*dH_dpz[1]*f[1]*hamil[1]+0.3535533905932737*dH_dpz[0]*hamil[0]*f[1])*volFact; 
  out[10] += (0.3535533905932737*dH_dpz[0]*hamil[7]*f[30]+0.3535533905932737*dH_dpz[1]*hamil[6]*f[30]+0.3535533905932737*dH_dpz[2]*hamil[5]*f[30]+0.3535533905932737*hamil[3]*dH_dpz[4]*f[30]+0.3535533905932737*dH_dpz[1]*hamil[7]*f[24]+0.3535533905932737*dH_dpz[0]*hamil[6]*f[24]+0.3535533905932737*dH_dpz[4]*hamil[5]*f[24]+0.3535533905932737*dH_dpz[2]*hamil[3]*f[24]+0.3535533905932737*dH_dpz[2]*hamil[7]*f[22]+0.3535533905932737*dH_dpz[4]*hamil[6]*f[22]+0.3535533905932737*dH_dpz[0]*hamil[5]*f[22]+0.3535533905932737*dH_dpz[1]*hamil[3]*f[22]+0.3535533905932737*dH_dpz[0]*hamil[4]*f[19]+0.3535533905932737*hamil[0]*dH_dpz[4]*f[19]+0.3535533905932737*dH_dpz[1]*hamil[2]*f[19]+0.3535533905932737*hamil[1]*dH_dpz[2]*f[19]+0.3535533905932737*dH_dpz[4]*hamil[7]*f[13]+0.3535533905932737*dH_dpz[2]*hamil[6]*f[13]+0.3535533905932737*dH_dpz[1]*hamil[5]*f[13]+0.3535533905932737*dH_dpz[0]*hamil[3]*f[13]+0.3535533905932737*dH_dpz[1]*hamil[4]*f[10]+0.3535533905932737*hamil[1]*dH_dpz[4]*f[10]+0.3535533905932737*dH_dpz[0]*hamil[2]*f[10]+0.3535533905932737*hamil[0]*dH_dpz[2]*f[10]+0.3535533905932737*dH_dpz[2]*hamil[4]*f[8]+0.3535533905932737*hamil[2]*dH_dpz[4]*f[8]+0.3535533905932737*dH_dpz[0]*hamil[1]*f[8]+0.3535533905932737*hamil[0]*dH_dpz[1]*f[8]+0.3535533905932737*f[2]*dH_dpz[4]*hamil[4]+0.3535533905932737*dH_dpz[2]*f[2]*hamil[2]+0.3535533905932737*dH_dpz[1]*hamil[1]*f[2]+0.3535533905932737*dH_dpz[0]*hamil[0]*f[2])*volFact; 
  out[11] += (0.3535533905932737*dH_dpz[0]*hamil[7]*f[31]+0.3535533905932737*dH_dpz[1]*hamil[6]*f[31]+0.3535533905932737*dH_dpz[2]*hamil[5]*f[31]+0.3535533905932737*hamil[3]*dH_dpz[4]*f[31]+0.3535533905932737*dH_dpz[1]*hamil[7]*f[28]+0.3535533905932737*dH_dpz[0]*hamil[6]*f[28]+0.3535533905932737*dH_dpz[4]*hamil[5]*f[28]+0.3535533905932737*dH_dpz[2]*hamil[3]*f[28]+0.3535533905932737*dH_dpz[2]*hamil[7]*f[27]+0.3535533905932737*dH_dpz[4]*hamil[6]*f[27]+0.3535533905932737*dH_dpz[0]*hamil[5]*f[27]+0.3535533905932737*dH_dpz[1]*hamil[3]*f[27]+0.3535533905932737*dH_dpz[0]*hamil[4]*f[26]+0.3535533905932737*hamil[0]*dH_dpz[4]*f[26]+0.3535533905932737*dH_dpz[1]*hamil[2]*f[26]+0.3535533905932737*hamil[1]*dH_dpz[2]*f[26]+0.3535533905932737*dH_dpz[4]*hamil[7]*f[20]+0.3535533905932737*dH_dpz[2]*hamil[6]*f[20]+0.3535533905932737*dH_dpz[1]*hamil[5]*f[20]+0.3535533905932737*dH_dpz[0]*hamil[3]*f[20]+0.3535533905932737*dH_dpz[1]*hamil[4]*f[17]+0.3535533905932737*hamil[1]*dH_dpz[4]*f[17]+0.3535533905932737*dH_dpz[0]*hamil[2]*f[17]+0.3535533905932737*hamil[0]*dH_dpz[2]*f[17]+0.3535533905932737*dH_dpz[2]*hamil[4]*f[16]+0.3535533905932737*hamil[2]*dH_dpz[4]*f[16]+0.3535533905932737*dH_dpz[0]*hamil[1]*f[16]+0.3535533905932737*hamil[0]*dH_dpz[1]*f[16]+0.3535533905932737*dH_dpz[4]*hamil[4]*f[6]+0.3535533905932737*dH_dpz[2]*hamil[2]*f[6]+0.3535533905932737*dH_dpz[1]*hamil[1]*f[6]+0.3535533905932737*dH_dpz[0]*hamil[0]*f[6])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_five_moments_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[3]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  double dv12 = 2.0/dxv[4]; 
  const double *jacob_vz = &jacob_vel[4]; 
  const double jacob_vz_inv = 1.0/jacob_vz[0]; 
  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += 2.8284271247461907*f[1]*volFact; 
  out[2] += 2.8284271247461907*f[2]*volFact; 
  out[3] += 2.8284271247461907*f[6]*volFact; 
  out[4] += (1.7320508075688772*hamil[7]*f[15]+1.7320508075688772*f[5]*hamil[5]+1.7320508075688772*f[4]*hamil[4]+1.7320508075688772*f[0]*hamil[1])*dv10*jacob_vx_inv*volFact; 
  out[5] += (1.7320508075688772*hamil[7]*f[23]+1.7320508075688772*hamil[5]*f[12]+1.7320508075688772*hamil[4]*f[9]+1.7320508075688772*f[1]*hamil[1])*dv10*jacob_vx_inv*volFact; 
  out[6] += (1.7320508075688772*hamil[7]*f[24]+1.7320508075688772*hamil[5]*f[13]+1.7320508075688772*hamil[4]*f[10]+1.7320508075688772*hamil[1]*f[2])*dv10*jacob_vx_inv*volFact; 
  out[7] += (1.7320508075688772*hamil[7]*f[28]+1.7320508075688772*hamil[5]*f[20]+1.7320508075688772*hamil[4]*f[17]+1.7320508075688772*hamil[1]*f[6])*dv10*jacob_vx_inv*volFact; 
  out[8] += (1.7320508075688772*hamil[7]*f[14]+1.7320508075688772*f[5]*hamil[6]+1.7320508075688772*f[3]*hamil[4]+1.7320508075688772*f[0]*hamil[2])*dv11*jacob_vy_inv*volFact; 
  out[9] += (1.7320508075688772*hamil[7]*f[21]+1.7320508075688772*hamil[6]*f[12]+1.7320508075688772*hamil[4]*f[7]+1.7320508075688772*f[1]*hamil[2])*dv11*jacob_vy_inv*volFact; 
  out[10] += (1.7320508075688772*hamil[7]*f[22]+1.7320508075688772*hamil[6]*f[13]+1.7320508075688772*hamil[4]*f[8]+1.7320508075688772*f[2]*hamil[2])*dv11*jacob_vy_inv*volFact; 
  out[11] += (1.7320508075688772*hamil[7]*f[27]+1.7320508075688772*hamil[6]*f[20]+1.7320508075688772*hamil[4]*f[16]+1.7320508075688772*hamil[2]*f[6])*dv11*jacob_vy_inv*volFact; 
  out[12] += (1.7320508075688772*hamil[7]*f[11]+1.7320508075688772*f[4]*hamil[6]+1.7320508075688772*f[3]*hamil[5]+1.7320508075688772*f[0]*hamil[3])*dv12*jacob_vz_inv*volFact; 
  out[13] += (1.7320508075688772*hamil[7]*f[18]+1.7320508075688772*hamil[6]*f[9]+1.7320508075688772*hamil[5]*f[7]+1.7320508075688772*f[1]*hamil[3])*dv12*jacob_vz_inv*volFact; 
  out[14] += (1.7320508075688772*hamil[7]*f[19]+1.7320508075688772*hamil[6]*f[10]+1.7320508075688772*hamil[5]*f[8]+1.7320508075688772*f[2]*hamil[3])*dv12*jacob_vz_inv*volFact; 
  out[15] += (1.7320508075688772*hamil[7]*f[26]+1.7320508075688772*hamil[6]*f[17]+1.7320508075688772*hamil[5]*f[16]+1.7320508075688772*hamil[3]*f[6])*dv12*jacob_vz_inv*volFact; 
  out[16] += (hamil[7]*f[25]+hamil[6]*f[15]+hamil[5]*f[14]+hamil[4]*f[11]+hamil[3]*f[5]+hamil[2]*f[4]+hamil[1]*f[3]+f[0]*hamil[0])*volFact; 
  out[17] += (hamil[7]*f[29]+hamil[6]*f[23]+hamil[5]*f[21]+hamil[4]*f[18]+hamil[3]*f[12]+hamil[2]*f[9]+hamil[1]*f[7]+hamil[0]*f[1])*volFact; 
  out[18] += (hamil[7]*f[30]+hamil[6]*f[24]+hamil[5]*f[22]+hamil[4]*f[19]+hamil[3]*f[13]+hamil[2]*f[10]+hamil[1]*f[8]+hamil[0]*f[2])*volFact; 
  out[19] += (hamil[7]*f[31]+hamil[6]*f[28]+hamil[5]*f[27]+hamil[4]*f[26]+hamil[3]*f[20]+hamil[2]*f[17]+hamil[1]*f[16]+hamil[0]*f[6])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_int_five_moments_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]*0.03125; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[3]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  double dv12 = 2.0/dxv[4]; 
  const double *jacob_vz = &jacob_vel[4]; 
  const double jacob_vz_inv = 1.0/jacob_vz[0]; 
  out[0] += 5.656854249492382*f[0]*volFact; 
  out[1] += (3.4641016151377544*hamil[7]*f[15]+3.4641016151377544*f[5]*hamil[5]+3.4641016151377544*f[4]*hamil[4]+3.4641016151377544*f[0]*hamil[1])*dv10*jacob_vx_inv*volFact; 
  out[2] += (3.4641016151377544*hamil[7]*f[14]+3.4641016151377544*f[5]*hamil[6]+3.4641016151377544*f[3]*hamil[4]+3.4641016151377544*f[0]*hamil[2])*dv11*jacob_vy_inv*volFact; 
  out[3] += (3.4641016151377544*hamil[7]*f[11]+3.4641016151377544*f[4]*hamil[6]+3.4641016151377544*f[3]*hamil[5]+3.4641016151377544*f[0]*hamil[3])*dv12*jacob_vz_inv*volFact; 
  out[4] += (2.0*hamil[7]*f[25]+2.0*hamil[6]*f[15]+2.0*hamil[5]*f[14]+2.0*hamil[4]*f[11]+2.0*hamil[3]*f[5]+2.0*hamil[2]*f[4]+2.0*hamil[1]*f[3]+2.0*f[0]*hamil[0])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_sparse_M1i_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[3]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  double dv12 = 2.0/dxv[4]; 
  const double *jacob_vz = &jacob_vel[4]; 
  const double jacob_vz_inv = 1.0/jacob_vz[0]; 
  out[0] += 1.7320508075688772*f[0]*hamil[1]*dv10*jacob_vx_inv*volFact; 
  out[1] += 1.7320508075688772*f[1]*hamil[1]*dv10*jacob_vx_inv*volFact; 
  out[2] += 1.7320508075688772*hamil[1]*f[2]*dv10*jacob_vx_inv*volFact; 
  out[3] += 1.7320508075688772*hamil[1]*f[6]*dv10*jacob_vx_inv*volFact; 
  out[4] += 1.7320508075688772*f[0]*hamil[2]*dv11*jacob_vy_inv*volFact; 
  out[5] += 1.7320508075688772*f[1]*hamil[2]*dv11*jacob_vy_inv*volFact; 
  out[6] += 1.7320508075688772*f[2]*hamil[2]*dv11*jacob_vy_inv*volFact; 
  out[7] += 1.7320508075688772*hamil[2]*f[6]*dv11*jacob_vy_inv*volFact; 
  out[8] += 1.7320508075688772*f[0]*hamil[3]*dv12*jacob_vz_inv*volFact; 
  out[9] += 1.7320508075688772*f[1]*hamil[3]*dv12*jacob_vz_inv*volFact; 
  out[10] += 1.7320508075688772*f[2]*hamil[3]*dv12*jacob_vz_inv*volFact; 
  out[11] += 1.7320508075688772*hamil[3]*f[6]*dv12*jacob_vz_inv*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_sparse_M2_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  out[0] += (hamil[3]*f[5]+hamil[2]*f[4]+hamil[1]*f[3]+f[0]*hamil[0])*volFact; 
  out[1] += (hamil[3]*f[12]+hamil[2]*f[9]+hamil[1]*f[7]+hamil[0]*f[1])*volFact; 
  out[2] += (hamil[3]*f[13]+hamil[2]*f[10]+hamil[1]*f[8]+hamil[0]*f[2])*volFact; 
  out[3] += (hamil[3]*f[20]+hamil[2]*f[17]+hamil[1]*f[16]+hamil[0]*f[6])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_sparse_M3i_2x3v_ser_p1(const double *w, const double *dxv, const int *idx,
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[3]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  double dv12 = 2.0/dxv[4]; 
  const double *jacob_vz = &jacob_vel[4]; 
  const double jacob_vz_inv = 1.0/jacob_vz[0]; 
  double dH_dpx[8] = {0.0}; 
  dH_dpx[0] = 1.7320508075688772*hamil[1]*dv10*jacob_vx_inv; 

  double dH_dpy[8] = {0.0}; 
  dH_dpy[0] = 1.7320508075688772*hamil[2]*dv11*jacob_vy_inv; 

  double dH_dpz[8] = {0.0}; 
  dH_dpz[0] = 1.7320508075688772*hamil[3]*dv12*jacob_vz_inv; 

  out[0] += (0.3535533905932737*dH_dpx[0]*hamil[3]*f[5]+0.3535533905932737*dH_dpx[0]*hamil[2]*f[4]+0.3535533905932737*dH_dpx[0]*hamil[1]*f[3]+0.3535533905932737*dH_dpx[0]*f[0]*hamil[0])*volFact; 
  out[1] += (0.3535533905932737*dH_dpx[0]*hamil[3]*f[12]+0.3535533905932737*dH_dpx[0]*hamil[2]*f[9]+0.3535533905932737*dH_dpx[0]*hamil[1]*f[7]+0.3535533905932737*dH_dpx[0]*hamil[0]*f[1])*volFact; 
  out[2] += (0.3535533905932737*dH_dpx[0]*hamil[3]*f[13]+0.3535533905932737*dH_dpx[0]*hamil[2]*f[10]+0.3535533905932737*dH_dpx[0]*hamil[1]*f[8]+0.3535533905932737*dH_dpx[0]*hamil[0]*f[2])*volFact; 
  out[3] += (0.3535533905932737*dH_dpx[0]*hamil[3]*f[20]+0.3535533905932737*dH_dpx[0]*hamil[2]*f[17]+0.3535533905932737*dH_dpx[0]*hamil[1]*f[16]+0.3535533905932737*dH_dpx[0]*hamil[0]*f[6])*volFact; 
  out[4] += (0.3535533905932737*dH_dpy[0]*hamil[3]*f[5]+0.3535533905932737*dH_dpy[0]*hamil[2]*f[4]+0.3535533905932737*dH_dpy[0]*hamil[1]*f[3]+0.3535533905932737*dH_dpy[0]*f[0]*hamil[0])*volFact; 
  out[5] += (0.3535533905932737*dH_dpy[0]*hamil[3]*f[12]+0.3535533905932737*dH_dpy[0]*hamil[2]*f[9]+0.3535533905932737*dH_dpy[0]*hamil[1]*f[7]+0.3535533905932737*dH_dpy[0]*hamil[0]*f[1])*volFact; 
  out[6] += (0.3535533905932737*dH_dpy[0]*hamil[3]*f[13]+0.3535533905932737*dH_dpy[0]*hamil[2]*f[10]+0.3535533905932737*dH_dpy[0]*hamil[1]*f[8]+0.3535533905932737*dH_dpy[0]*hamil[0]*f[2])*volFact; 
  out[7] += (0.3535533905932737*dH_dpy[0]*hamil[3]*f[20]+0.3535533905932737*dH_dpy[0]*hamil[2]*f[17]+0.3535533905932737*dH_dpy[0]*hamil[1]*f[16]+0.3535533905932737*dH_dpy[0]*hamil[0]*f[6])*volFact; 
  out[8] += (0.3535533905932737*dH_dpz[0]*hamil[3]*f[5]+0.3535533905932737*dH_dpz[0]*hamil[2]*f[4]+0.3535533905932737*dH_dpz[0]*hamil[1]*f[3]+0.3535533905932737*dH_dpz[0]*f[0]*hamil[0])*volFact; 
  out[9] += (0.3535533905932737*dH_dpz[0]*hamil[3]*f[12]+0.3535533905932737*dH_dpz[0]*hamil[2]*f[9]+0.3535533905932737*dH_dpz[0]*hamil[1]*f[7]+0.3535533905932737*dH_dpz[0]*hamil[0]*f[1])*volFact; 
  out[10] += (0.3535533905932737*dH_dpz[0]*hamil[3]*f[13]+0.3535533905932737*dH_dpz[0]*hamil[2]*f[10]+0.3535533905932737*dH_dpz[0]*hamil[1]*f[8]+0.3535533905932737*dH_dpz[0]*hamil[0]*f[2])*volFact; 
  out[11] += (0.3535533905932737*dH_dpz[0]*hamil[3]*f[20]+0.3535533905932737*dH_dpz[0]*hamil[2]*f[17]+0.3535533905932737*dH_dpz[0]*hamil[1]*f[16]+0.3535533905932737*dH_dpz[0]*hamil[0]*f[6])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_sparse_five_moments_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[3]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  double dv12 = 2.0/dxv[4]; 
  const double *jacob_vz = &jacob_vel[4]; 
  const double jacob_vz_inv = 1.0/jacob_vz[0]; 
  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += 2.8284271247461907*f[1]*volFact; 
  out[2] += 2.8284271247461907*f[2]*volFact; 
  out[3] += 2.8284271247461907*f[6]*volFact; 
  out[4] += 1.7320508075688772*f[0]*hamil[1]*dv10*jacob_vx_inv*volFact; 
  out[5] += 1.7320508075688772*f[1]*hamil[1]*dv10*jacob_vx_inv*volFact; 
  out[6] += 1.7320508075688772*hamil[1]*f[2]*dv10*jacob_vx_inv*volFact; 
  out[7] += 1.7320508075688772*hamil[1]*f[6]*dv10*jacob_vx_inv*volFact; 
  out[8] += 1.7320508075688772*f[0]*hamil[2]*dv11*jacob_vy_inv*volFact; 
  out[9] += 1.7320508075688772*f[1]*hamil[2]*dv11*jacob_vy_inv*volFact; 
  out[10] += 1.7320508075688772*f[2]*hamil[2]*dv11*jacob_vy_inv*volFact; 
  out[11] += 1.7320508075688772*hamil[2]*f[6]*dv11*jacob_vy_inv*volFact; 
  out[12] += 1.7320508075688772*f[0]*hamil[3]*dv12*jacob_vz_inv*volFact; 
  out[13] += 1.7320508075688772*f[1]*hamil[3]*dv12*jacob_vz_inv*volFact; 
  out[14] += 1.7320508075688772*f[2]*hamil[3]*dv12*jacob_vz_inv*volFact; 
  out[15] += 1.7320508075688772*hamil[3]*f[6]*dv12*jacob_vz_inv*volFact; 
  out[16] += (hamil[3]*f[5]+hamil[2]*f[4]+hamil[1]*f[3]+f[0]*hamil[0])*volFact; 
  out[17] += (hamil[3]*f[12]+hamil[2]*f[9]+hamil[1]*f[7]+hamil[0]*f[1])*volFact; 
  out[18] += (hamil[3]*f[13]+hamil[2]*f[10]+hamil[1]*f[8]+hamil[0]*f[2])*volFact; 
  out[19] += (hamil[3]*f[20]+hamil[2]*f[17]+hamil[1]*f[16]+hamil[0]*f[6])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_sparse_int_five_moments_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]*0.03125; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[3]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  double dv12 = 2.0/dxv[4]; 
  const double *jacob_vz = &jacob_vel[4]; 
  const double jacob_vz_inv = 1.0/jacob_vz[0]; 
  out[0] += 5.656854249492382*f[0]*volFact; 
  out[1] += 3.4641016151377544*f[0]*hamil[1]*dv10*jacob_vx_inv*volFact; 
  out[2] += 3.4641016151377544*f[0]*hamil[2]*dv11*jacob_vy_inv*volFact; 
  out[3] += 3.4641016151377544*f[0]*hamil[3]*dv12*jacob_vz_inv*volFact; 
  out[4] += (2.0*hamil[3]*f[5]+2.0*hamil[2]*f[4]+2.0*hamil[1]*f[3]+2.0*f[0]*hamil[0])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_phase_M1i_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[3]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  double dv12 = 2.0/dxv[4]; 
  const double *jacob_vz = &jacob_vel[4]; 
  const double jacob_vz_inv = 1.0/jacob_vz[0]; 
  out[0] += (0.8660254037844386*f[28]*hamil[31]+0.8660254037844386*f[24]*hamil[30]+0.8660254037844386*f[23]*hamil[29]+0.8660254037844386*f[20]*hamil[27]+0.8660254037844386*f[17]*hamil[26]+0.8660254037844386*f[15]*hamil[25]+0.8660254037844386*f[13]*hamil[22]+0.8660254037844386*f[12]*hamil[21]+0.8660254037844386*f[10]*hamil[19]+0.8660254037844386*f[9]*hamil[18]+0.8660254037844386*f[6]*hamil[16]+0.8660254037844386*f[5]*hamil[14]+0.8660254037844386*f[4]*hamil[11]+0.8660254037844386*f[2]*hamil[8]+0.8660254037844386*f[1]*hamil[7]+0.8660254037844386*f[0]*hamil[3])*dv10*jacob_vx_inv*volFact; 
  out[1] += (0.8660254037844386*f[24]*hamil[31]+0.8660254037844386*f[28]*hamil[30]+0.8660254037844386*f[15]*hamil[29]+0.8660254037844386*f[13]*hamil[27]+0.8660254037844386*f[10]*hamil[26]+0.8660254037844386*f[23]*hamil[25]+0.8660254037844386*f[20]*hamil[22]+0.8660254037844386*f[5]*hamil[21]+0.8660254037844386*f[17]*hamil[19]+0.8660254037844386*f[4]*hamil[18]+0.8660254037844386*f[2]*hamil[16]+0.8660254037844386*f[12]*hamil[14]+0.8660254037844386*f[9]*hamil[11]+0.8660254037844386*f[6]*hamil[8]+0.8660254037844386*f[0]*hamil[7]+0.8660254037844386*f[1]*hamil[3])*dv10*jacob_vx_inv*volFact; 
  out[2] += (0.8660254037844386*f[23]*hamil[31]+0.8660254037844386*f[15]*hamil[30]+0.8660254037844386*f[28]*hamil[29]+0.8660254037844386*f[12]*hamil[27]+0.8660254037844386*f[9]*hamil[26]+0.8660254037844386*f[24]*hamil[25]+0.8660254037844386*f[5]*hamil[22]+0.8660254037844386*f[20]*hamil[21]+0.8660254037844386*f[4]*hamil[19]+0.8660254037844386*f[17]*hamil[18]+0.8660254037844386*f[1]*hamil[16]+0.8660254037844386*f[13]*hamil[14]+0.8660254037844386*f[10]*hamil[11]+0.8660254037844386*f[0]*hamil[8]+0.8660254037844386*f[6]*hamil[7]+0.8660254037844386*f[2]*hamil[3])*dv10*jacob_vx_inv*volFact; 
  out[3] += (0.8660254037844386*f[15]*hamil[31]+0.8660254037844386*f[23]*hamil[30]+0.8660254037844386*f[24]*hamil[29]+0.8660254037844386*hamil[25]*f[28]+0.8660254037844386*f[5]*hamil[27]+0.8660254037844386*f[4]*hamil[26]+0.8660254037844386*f[12]*hamil[22]+0.8660254037844386*f[13]*hamil[21]+0.8660254037844386*hamil[14]*f[20]+0.8660254037844386*f[9]*hamil[19]+0.8660254037844386*f[10]*hamil[18]+0.8660254037844386*hamil[11]*f[17]+0.8660254037844386*f[0]*hamil[16]+0.8660254037844386*f[1]*hamil[8]+0.8660254037844386*f[2]*hamil[7]+0.8660254037844386*hamil[3]*f[6])*dv10*jacob_vx_inv*volFact; 
  out[4] += (0.8660254037844386*f[27]*hamil[31]+0.8660254037844386*f[22]*hamil[30]+0.8660254037844386*f[21]*hamil[29]+0.8660254037844386*f[20]*hamil[28]+0.8660254037844386*f[16]*hamil[26]+0.8660254037844386*f[14]*hamil[25]+0.8660254037844386*f[13]*hamil[24]+0.8660254037844386*f[12]*hamil[23]+0.8660254037844386*f[8]*hamil[19]+0.8660254037844386*f[7]*hamil[18]+0.8660254037844386*f[6]*hamil[17]+0.8660254037844386*f[5]*hamil[15]+0.8660254037844386*f[3]*hamil[11]+0.8660254037844386*f[2]*hamil[10]+0.8660254037844386*f[1]*hamil[9]+0.8660254037844386*f[0]*hamil[4])*dv11*jacob_vy_inv*volFact; 
  out[5] += (0.8660254037844386*f[22]*hamil[31]+0.8660254037844386*f[27]*hamil[30]+0.8660254037844386*f[14]*hamil[29]+0.8660254037844386*f[13]*hamil[28]+0.8660254037844386*f[8]*hamil[26]+0.8660254037844386*f[21]*hamil[25]+0.8660254037844386*f[20]*hamil[24]+0.8660254037844386*f[5]*hamil[23]+0.8660254037844386*f[16]*hamil[19]+0.8660254037844386*f[3]*hamil[18]+0.8660254037844386*f[2]*hamil[17]+0.8660254037844386*f[12]*hamil[15]+0.8660254037844386*f[7]*hamil[11]+0.8660254037844386*f[6]*hamil[10]+0.8660254037844386*f[0]*hamil[9]+0.8660254037844386*f[1]*hamil[4])*dv11*jacob_vy_inv*volFact; 
  out[6] += (0.8660254037844386*f[21]*hamil[31]+0.8660254037844386*f[14]*hamil[30]+0.8660254037844386*f[27]*hamil[29]+0.8660254037844386*f[12]*hamil[28]+0.8660254037844386*f[7]*hamil[26]+0.8660254037844386*f[22]*hamil[25]+0.8660254037844386*f[5]*hamil[24]+0.8660254037844386*f[20]*hamil[23]+0.8660254037844386*f[3]*hamil[19]+0.8660254037844386*f[16]*hamil[18]+0.8660254037844386*f[1]*hamil[17]+0.8660254037844386*f[13]*hamil[15]+0.8660254037844386*f[8]*hamil[11]+0.8660254037844386*f[0]*hamil[10]+0.8660254037844386*f[6]*hamil[9]+0.8660254037844386*f[2]*hamil[4])*dv11*jacob_vy_inv*volFact; 
  out[7] += (0.8660254037844386*f[14]*hamil[31]+0.8660254037844386*f[21]*hamil[30]+0.8660254037844386*f[22]*hamil[29]+0.8660254037844386*f[5]*hamil[28]+0.8660254037844386*hamil[25]*f[27]+0.8660254037844386*f[3]*hamil[26]+0.8660254037844386*f[12]*hamil[24]+0.8660254037844386*f[13]*hamil[23]+0.8660254037844386*hamil[15]*f[20]+0.8660254037844386*f[7]*hamil[19]+0.8660254037844386*f[8]*hamil[18]+0.8660254037844386*f[0]*hamil[17]+0.8660254037844386*hamil[11]*f[16]+0.8660254037844386*f[1]*hamil[10]+0.8660254037844386*f[2]*hamil[9]+0.8660254037844386*hamil[4]*f[6])*dv11*jacob_vy_inv*volFact; 
  out[8] += (0.8660254037844386*f[26]*hamil[31]+0.8660254037844386*f[19]*hamil[30]+0.8660254037844386*f[18]*hamil[29]+0.8660254037844386*f[17]*hamil[28]+0.8660254037844386*f[16]*hamil[27]+0.8660254037844386*f[11]*hamil[25]+0.8660254037844386*f[10]*hamil[24]+0.8660254037844386*f[9]*hamil[23]+0.8660254037844386*f[8]*hamil[22]+0.8660254037844386*f[7]*hamil[21]+0.8660254037844386*f[6]*hamil[20]+0.8660254037844386*f[4]*hamil[15]+0.8660254037844386*f[3]*hamil[14]+0.8660254037844386*f[2]*hamil[13]+0.8660254037844386*f[1]*hamil[12]+0.8660254037844386*f[0]*hamil[5])*dv12*jacob_vz_inv*volFact; 
  out[9] += (0.8660254037844386*f[19]*hamil[31]+0.8660254037844386*f[26]*hamil[30]+0.8660254037844386*f[11]*hamil[29]+0.8660254037844386*f[10]*hamil[28]+0.8660254037844386*f[8]*hamil[27]+0.8660254037844386*f[18]*hamil[25]+0.8660254037844386*f[17]*hamil[24]+0.8660254037844386*f[4]*hamil[23]+0.8660254037844386*f[16]*hamil[22]+0.8660254037844386*f[3]*hamil[21]+0.8660254037844386*f[2]*hamil[20]+0.8660254037844386*f[9]*hamil[15]+0.8660254037844386*f[7]*hamil[14]+0.8660254037844386*f[6]*hamil[13]+0.8660254037844386*f[0]*hamil[12]+0.8660254037844386*f[1]*hamil[5])*dv12*jacob_vz_inv*volFact; 
  out[10] += (0.8660254037844386*f[18]*hamil[31]+0.8660254037844386*f[11]*hamil[30]+0.8660254037844386*f[26]*hamil[29]+0.8660254037844386*f[9]*hamil[28]+0.8660254037844386*f[7]*hamil[27]+0.8660254037844386*f[19]*hamil[25]+0.8660254037844386*f[4]*hamil[24]+0.8660254037844386*f[17]*hamil[23]+0.8660254037844386*f[3]*hamil[22]+0.8660254037844386*f[16]*hamil[21]+0.8660254037844386*f[1]*hamil[20]+0.8660254037844386*f[10]*hamil[15]+0.8660254037844386*f[8]*hamil[14]+0.8660254037844386*f[0]*hamil[13]+0.8660254037844386*f[6]*hamil[12]+0.8660254037844386*f[2]*hamil[5])*dv12*jacob_vz_inv*volFact; 
  out[11] += (0.8660254037844386*f[11]*hamil[31]+0.8660254037844386*f[18]*hamil[30]+0.8660254037844386*f[19]*hamil[29]+0.8660254037844386*f[4]*hamil[28]+0.8660254037844386*f[3]*hamil[27]+0.8660254037844386*hamil[25]*f[26]+0.8660254037844386*f[9]*hamil[24]+0.8660254037844386*f[10]*hamil[23]+0.8660254037844386*f[7]*hamil[22]+0.8660254037844386*f[8]*hamil[21]+0.8660254037844386*f[0]*hamil[20]+0.8660254037844386*hamil[15]*f[17]+0.8660254037844386*hamil[14]*f[16]+0.8660254037844386*f[1]*hamil[13]+0.8660254037844386*f[2]*hamil[12]+0.8660254037844386*hamil[5]*f[6])*dv12*jacob_vz_inv*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_phase_M2_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  out[0] += (0.5*f[31]*hamil[31]+0.5*f[30]*hamil[30]+0.5*f[29]*hamil[29]+0.5*f[28]*hamil[28]+0.5*f[27]*hamil[27]+0.5*f[26]*hamil[26]+0.5*f[25]*hamil[25]+0.5*f[24]*hamil[24]+0.5*f[23]*hamil[23]+0.5*f[22]*hamil[22]+0.5*f[21]*hamil[21]+0.5*f[20]*hamil[20]+0.5*f[19]*hamil[19]+0.5*f[18]*hamil[18]+0.5*f[17]*hamil[17]+0.5*f[16]*hamil[16]+0.5*f[15]*hamil[15]+0.5*f[14]*hamil[14]+0.5*f[13]*hamil[13]+0.5*f[12]*hamil[12]+0.5*f[11]*hamil[11]+0.5*f[10]*hamil[10]+0.5*f[9]*hamil[9]+0.5*f[8]*hamil[8]+0.5*f[7]*hamil[7]+0.5*f[6]*hamil[6]+0.5*f[5]*hamil[5]+0.5*f[4]*hamil[4]+0.5*f[3]*hamil[3]+0.5*f[2]*hamil[2]+0.5*f[1]*hamil[1]+0.5*f[0]*hamil[0])*volFact; 
  out[1] += (0.5*f[30]*hamil[31]+0.5*hamil[30]*f[31]+0.5*f[25]*hamil[29]+0.5*hamil[25]*f[29]+0.5*f[24]*hamil[28]+0.5*hamil[24]*f[28]+0.5*f[22]*hamil[27]+0.5*hamil[22]*f[27]+0.5*f[19]*hamil[26]+0.5*hamil[19]*f[26]+0.5*f[15]*hamil[23]+0.5*hamil[15]*f[23]+0.5*f[14]*hamil[21]+0.5*hamil[14]*f[21]+0.5*f[13]*hamil[20]+0.5*hamil[13]*f[20]+0.5*f[11]*hamil[18]+0.5*hamil[11]*f[18]+0.5*f[10]*hamil[17]+0.5*hamil[10]*f[17]+0.5*f[8]*hamil[16]+0.5*hamil[8]*f[16]+0.5*f[5]*hamil[12]+0.5*hamil[5]*f[12]+0.5*f[4]*hamil[9]+0.5*hamil[4]*f[9]+0.5*f[3]*hamil[7]+0.5*hamil[3]*f[7]+0.5*f[2]*hamil[6]+0.5*hamil[2]*f[6]+0.5*f[0]*hamil[1]+0.5*hamil[0]*f[1])*volFact; 
  out[2] += (0.5*f[29]*hamil[31]+0.5*hamil[29]*f[31]+0.5*f[25]*hamil[30]+0.5*hamil[25]*f[30]+0.5*f[23]*hamil[28]+0.5*hamil[23]*f[28]+0.5*f[21]*hamil[27]+0.5*hamil[21]*f[27]+0.5*f[18]*hamil[26]+0.5*hamil[18]*f[26]+0.5*f[15]*hamil[24]+0.5*hamil[15]*f[24]+0.5*f[14]*hamil[22]+0.5*hamil[14]*f[22]+0.5*f[12]*hamil[20]+0.5*hamil[12]*f[20]+0.5*f[11]*hamil[19]+0.5*hamil[11]*f[19]+0.5*f[9]*hamil[17]+0.5*hamil[9]*f[17]+0.5*f[7]*hamil[16]+0.5*hamil[7]*f[16]+0.5*f[5]*hamil[13]+0.5*hamil[5]*f[13]+0.5*f[4]*hamil[10]+0.5*hamil[4]*f[10]+0.5*f[3]*hamil[8]+0.5*hamil[3]*f[8]+0.5*f[1]*hamil[6]+0.5*hamil[1]*f[6]+0.5*f[0]*hamil[2]+0.5*hamil[0]*f[2])*volFact; 
  out[3] += (0.5*f[25]*hamil[31]+0.5*hamil[25]*f[31]+0.5*f[29]*hamil[30]+0.5*hamil[29]*f[30]+0.5*f[15]*hamil[28]+0.5*hamil[15]*f[28]+0.5*f[14]*hamil[27]+0.5*hamil[14]*f[27]+0.5*f[11]*hamil[26]+0.5*hamil[11]*f[26]+0.5*f[23]*hamil[24]+0.5*hamil[23]*f[24]+0.5*f[21]*hamil[22]+0.5*hamil[21]*f[22]+0.5*f[5]*hamil[20]+0.5*hamil[5]*f[20]+0.5*f[18]*hamil[19]+0.5*hamil[18]*f[19]+0.5*f[4]*hamil[17]+0.5*hamil[4]*f[17]+0.5*f[3]*hamil[16]+0.5*hamil[3]*f[16]+0.5*f[12]*hamil[13]+0.5*hamil[12]*f[13]+0.5*f[9]*hamil[10]+0.5*hamil[9]*f[10]+0.5*f[7]*hamil[8]+0.5*hamil[7]*f[8]+0.5*f[0]*hamil[6]+0.5*hamil[0]*f[6]+0.5*f[1]*hamil[2]+0.5*hamil[1]*f[2])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_phase_five_moments_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[3]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  double dv12 = 2.0/dxv[4]; 
  const double *jacob_vz = &jacob_vel[4]; 
  const double jacob_vz_inv = 1.0/jacob_vz[0]; 
  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += 2.8284271247461907*f[1]*volFact; 
  out[2] += 2.8284271247461907*f[2]*volFact; 
  out[3] += 2.8284271247461907*f[6]*volFact; 
  out[4] += (0.8660254037844386*f[28]*hamil[31]+0.8660254037844386*f[24]*hamil[30]+0.8660254037844386*f[23]*hamil[29]+0.8660254037844386*f[20]*hamil[27]+0.8660254037844386*f[17]*hamil[26]+0.8660254037844386*f[15]*hamil[25]+0.8660254037844386*f[13]*hamil[22]+0.8660254037844386*f[12]*hamil[21]+0.8660254037844386*f[10]*hamil[19]+0.8660254037844386*f[9]*hamil[18]+0.8660254037844386*f[6]*hamil[16]+0.8660254037844386*f[5]*hamil[14]+0.8660254037844386*f[4]*hamil[11]+0.8660254037844386*f[2]*hamil[8]+0.8660254037844386*f[1]*hamil[7]+0.8660254037844386*f[0]*hamil[3])*dv10*jacob_vx_inv*volFact; 
  out[5] += (0.8660254037844386*f[24]*hamil[31]+0.8660254037844386*f[28]*hamil[30]+0.8660254037844386*f[15]*hamil[29]+0.8660254037844386*f[13]*hamil[27]+0.8660254037844386*f[10]*hamil[26]+0.8660254037844386*f[23]*hamil[25]+0.8660254037844386*f[20]*hamil[22]+0.8660254037844386*f[5]*hamil[21]+0.8660254037844386*f[17]*hamil[19]+0.8660254037844386*f[4]*hamil[18]+0.8660254037844386*f[2]*hamil[16]+0.8660254037844386*f[12]*hamil[14]+0.8660254037844386*f[9]*hamil[11]+0.8660254037844386*f[6]*hamil[8]+0.8660254037844386*f[0]*hamil[7]+0.8660254037844386*f[1]*hamil[3])*dv10*jacob_vx_inv*volFact; 
  out[6] += (0.8660254037844386*f[23]*hamil[31]+0.8660254037844386*f[15]*hamil[30]+0.8660254037844386*f[28]*hamil[29]+0.8660254037844386*f[12]*hamil[27]+0.8660254037844386*f[9]*hamil[26]+0.8660254037844386*f[24]*hamil[25]+0.8660254037844386*f[5]*hamil[22]+0.8660254037844386*f[20]*hamil[21]+0.8660254037844386*f[4]*hamil[19]+0.8660254037844386*f[17]*hamil[18]+0.8660254037844386*f[1]*hamil[16]+0.8660254037844386*f[13]*hamil[14]+0.8660254037844386*f[10]*hamil[11]+0.8660254037844386*f[0]*hamil[8]+0.8660254037844386*f[6]*hamil[7]+0.8660254037844386*f[2]*hamil[3])*dv10*jacob_vx_inv*volFact; 
  out[7] += (0.8660254037844386*f[15]*hamil[31]+0.8660254037844386*f[23]*hamil[30]+0.8660254037844386*f[24]*hamil[29]+0.8660254037844386*hamil[25]*f[28]+0.8660254037844386*f[5]*hamil[27]+0.8660254037844386*f[4]*hamil[26]+0.8660254037844386*f[12]*hamil[22]+0.8660254037844386*f[13]*hamil[21]+0.8660254037844386*hamil[14]*f[20]+0.8660254037844386*f[9]*hamil[19]+0.8660254037844386*f[10]*hamil[18]+0.8660254037844386*hamil[11]*f[17]+0.8660254037844386*f[0]*hamil[16]+0.8660254037844386*f[1]*hamil[8]+0.8660254037844386*f[2]*hamil[7]+0.8660254037844386*hamil[3]*f[6])*dv10*jacob_vx_inv*volFact; 
  out[8] += (0.8660254037844386*f[27]*hamil[31]+0.8660254037844386*f[22]*hamil[30]+0.8660254037844386*f[21]*hamil[29]+0.8660254037844386*f[20]*hamil[28]+0.8660254037844386*f[16]*hamil[26]+0.8660254037844386*f[14]*hamil[25]+0.8660254037844386*f[13]*hamil[24]+0.8660254037844386*f[12]*hamil[23]+0.8660254037844386*f[8]*hamil[19]+0.8660254037844386*f[7]*hamil[18]+0.8660254037844386*f[6]*hamil[17]+0.8660254037844386*f[5]*hamil[15]+0.8660254037844386*f[3]*hamil[11]+0.8660254037844386*f[2]*hamil[10]+0.8660254037844386*f[1]*hamil[9]+0.8660254037844386*f[0]*hamil[4])*dv11*jacob_vy_inv*volFact; 
  out[9] += (0.8660254037844386*f[22]*hamil[31]+0.8660254037844386*f[27]*hamil[30]+0.8660254037844386*f[14]*hamil[29]+0.8660254037844386*f[13]*hamil[28]+0.8660254037844386*f[8]*hamil[26]+0.8660254037844386*f[21]*hamil[25]+0.8660254037844386*f[20]*hamil[24]+0.8660254037844386*f[5]*hamil[23]+0.8660254037844386*f[16]*hamil[19]+0.8660254037844386*f[3]*hamil[18]+0.8660254037844386*f[2]*hamil[17]+0.8660254037844386*f[12]*hamil[15]+0.8660254037844386*f[7]*hamil[11]+0.8660254037844386*f[6]*hamil[10]+0.8660254037844386*f[0]*hamil[9]+0.8660254037844386*f[1]*hamil[4])*dv11*jacob_vy_inv*volFact; 
  out[10] += (0.8660254037844386*f[21]*hamil[31]+0.8660254037844386*f[14]*hamil[30]+0.8660254037844386*f[27]*hamil[29]+0.8660254037844386*f[12]*hamil[28]+0.8660254037844386*f[7]*hamil[26]+0.8660254037844386*f[22]*hamil[25]+0.8660254037844386*f[5]*hamil[24]+0.8660254037844386*f[20]*hamil[23]+0.8660254037844386*f[3]*hamil[19]+0.8660254037844386*f[16]*hamil[18]+0.8660254037844386*f[1]*hamil[17]+0.8660254037844386*f[13]*hamil[15]+0.8660254037844386*f[8]*hamil[11]+0.8660254037844386*f[0]*hamil[10]+0.8660254037844386*f[6]*hamil[9]+0.8660254037844386*f[2]*hamil[4])*dv11*jacob_vy_inv*volFact; 
  out[11] += (0.8660254037844386*f[14]*hamil[31]+0.8660254037844386*f[21]*hamil[30]+0.8660254037844386*f[22]*hamil[29]+0.8660254037844386*f[5]*hamil[28]+0.8660254037844386*hamil[25]*f[27]+0.8660254037844386*f[3]*hamil[26]+0.8660254037844386*f[12]*hamil[24]+0.8660254037844386*f[13]*hamil[23]+0.8660254037844386*hamil[15]*f[20]+0.8660254037844386*f[7]*hamil[19]+0.8660254037844386*f[8]*hamil[18]+0.8660254037844386*f[0]*hamil[17]+0.8660254037844386*hamil[11]*f[16]+0.8660254037844386*f[1]*hamil[10]+0.8660254037844386*f[2]*hamil[9]+0.8660254037844386*hamil[4]*f[6])*dv11*jacob_vy_inv*volFact; 
  out[12] += (0.8660254037844386*f[26]*hamil[31]+0.8660254037844386*f[19]*hamil[30]+0.8660254037844386*f[18]*hamil[29]+0.8660254037844386*f[17]*hamil[28]+0.8660254037844386*f[16]*hamil[27]+0.8660254037844386*f[11]*hamil[25]+0.8660254037844386*f[10]*hamil[24]+0.8660254037844386*f[9]*hamil[23]+0.8660254037844386*f[8]*hamil[22]+0.8660254037844386*f[7]*hamil[21]+0.8660254037844386*f[6]*hamil[20]+0.8660254037844386*f[4]*hamil[15]+0.8660254037844386*f[3]*hamil[14]+0.8660254037844386*f[2]*hamil[13]+0.8660254037844386*f[1]*hamil[12]+0.8660254037844386*f[0]*hamil[5])*dv12*jacob_vz_inv*volFact; 
  out[13] += (0.8660254037844386*f[19]*hamil[31]+0.8660254037844386*f[26]*hamil[30]+0.8660254037844386*f[11]*hamil[29]+0.8660254037844386*f[10]*hamil[28]+0.8660254037844386*f[8]*hamil[27]+0.8660254037844386*f[18]*hamil[25]+0.8660254037844386*f[17]*hamil[24]+0.8660254037844386*f[4]*hamil[23]+0.8660254037844386*f[16]*hamil[22]+0.8660254037844386*f[3]*hamil[21]+0.8660254037844386*f[2]*hamil[20]+0.8660254037844386*f[9]*hamil[15]+0.8660254037844386*f[7]*hamil[14]+0.8660254037844386*f[6]*hamil[13]+0.8660254037844386*f[0]*hamil[12]+0.8660254037844386*f[1]*hamil[5])*dv12*jacob_vz_inv*volFact; 
  out[14] += (0.8660254037844386*f[18]*hamil[31]+0.8660254037844386*f[11]*hamil[30]+0.8660254037844386*f[26]*hamil[29]+0.8660254037844386*f[9]*hamil[28]+0.8660254037844386*f[7]*hamil[27]+0.8660254037844386*f[19]*hamil[25]+0.8660254037844386*f[4]*hamil[24]+0.8660254037844386*f[17]*hamil[23]+0.8660254037844386*f[3]*hamil[22]+0.8660254037844386*f[16]*hamil[21]+0.8660254037844386*f[1]*hamil[20]+0.8660254037844386*f[10]*hamil[15]+0.8660254037844386*f[8]*hamil[14]+0.8660254037844386*f[0]*hamil[13]+0.8660254037844386*f[6]*hamil[12]+0.8660254037844386*f[2]*hamil[5])*dv12*jacob_vz_inv*volFact; 
  out[15] += (0.8660254037844386*f[11]*hamil[31]+0.8660254037844386*f[18]*hamil[30]+0.8660254037844386*f[19]*hamil[29]+0.8660254037844386*f[4]*hamil[28]+0.8660254037844386*f[3]*hamil[27]+0.8660254037844386*hamil[25]*f[26]+0.8660254037844386*f[9]*hamil[24]+0.8660254037844386*f[10]*hamil[23]+0.8660254037844386*f[7]*hamil[22]+0.8660254037844386*f[8]*hamil[21]+0.8660254037844386*f[0]*hamil[20]+0.8660254037844386*hamil[15]*f[17]+0.8660254037844386*hamil[14]*f[16]+0.8660254037844386*f[1]*hamil[13]+0.8660254037844386*f[2]*hamil[12]+0.8660254037844386*hamil[5]*f[6])*dv12*jacob_vz_inv*volFact; 
  out[16] += (0.5*f[31]*hamil[31]+0.5*f[30]*hamil[30]+0.5*f[29]*hamil[29]+0.5*f[28]*hamil[28]+0.5*f[27]*hamil[27]+0.5*f[26]*hamil[26]+0.5*f[25]*hamil[25]+0.5*f[24]*hamil[24]+0.5*f[23]*hamil[23]+0.5*f[22]*hamil[22]+0.5*f[21]*hamil[21]+0.5*f[20]*hamil[20]+0.5*f[19]*hamil[19]+0.5*f[18]*hamil[18]+0.5*f[17]*hamil[17]+0.5*f[16]*hamil[16]+0.5*f[15]*hamil[15]+0.5*f[14]*hamil[14]+0.5*f[13]*hamil[13]+0.5*f[12]*hamil[12]+0.5*f[11]*hamil[11]+0.5*f[10]*hamil[10]+0.5*f[9]*hamil[9]+0.5*f[8]*hamil[8]+0.5*f[7]*hamil[7]+0.5*f[6]*hamil[6]+0.5*f[5]*hamil[5]+0.5*f[4]*hamil[4]+0.5*f[3]*hamil[3]+0.5*f[2]*hamil[2]+0.5*f[1]*hamil[1]+0.5*f[0]*hamil[0])*volFact; 
  out[17] += (0.5*f[30]*hamil[31]+0.5*hamil[30]*f[31]+0.5*f[25]*hamil[29]+0.5*hamil[25]*f[29]+0.5*f[24]*hamil[28]+0.5*hamil[24]*f[28]+0.5*f[22]*hamil[27]+0.5*hamil[22]*f[27]+0.5*f[19]*hamil[26]+0.5*hamil[19]*f[26]+0.5*f[15]*hamil[23]+0.5*hamil[15]*f[23]+0.5*f[14]*hamil[21]+0.5*hamil[14]*f[21]+0.5*f[13]*hamil[20]+0.5*hamil[13]*f[20]+0.5*f[11]*hamil[18]+0.5*hamil[11]*f[18]+0.5*f[10]*hamil[17]+0.5*hamil[10]*f[17]+0.5*f[8]*hamil[16]+0.5*hamil[8]*f[16]+0.5*f[5]*hamil[12]+0.5*hamil[5]*f[12]+0.5*f[4]*hamil[9]+0.5*hamil[4]*f[9]+0.5*f[3]*hamil[7]+0.5*hamil[3]*f[7]+0.5*f[2]*hamil[6]+0.5*hamil[2]*f[6]+0.5*f[0]*hamil[1]+0.5*hamil[0]*f[1])*volFact; 
  out[18] += (0.5*f[29]*hamil[31]+0.5*hamil[29]*f[31]+0.5*f[25]*hamil[30]+0.5*hamil[25]*f[30]+0.5*f[23]*hamil[28]+0.5*hamil[23]*f[28]+0.5*f[21]*hamil[27]+0.5*hamil[21]*f[27]+0.5*f[18]*hamil[26]+0.5*hamil[18]*f[26]+0.5*f[15]*hamil[24]+0.5*hamil[15]*f[24]+0.5*f[14]*hamil[22]+0.5*hamil[14]*f[22]+0.5*f[12]*hamil[20]+0.5*hamil[12]*f[20]+0.5*f[11]*hamil[19]+0.5*hamil[11]*f[19]+0.5*f[9]*hamil[17]+0.5*hamil[9]*f[17]+0.5*f[7]*hamil[16]+0.5*hamil[7]*f[16]+0.5*f[5]*hamil[13]+0.5*hamil[5]*f[13]+0.5*f[4]*hamil[10]+0.5*hamil[4]*f[10]+0.5*f[3]*hamil[8]+0.5*hamil[3]*f[8]+0.5*f[1]*hamil[6]+0.5*hamil[1]*f[6]+0.5*f[0]*hamil[2]+0.5*hamil[0]*f[2])*volFact; 
  out[19] += (0.5*f[25]*hamil[31]+0.5*hamil[25]*f[31]+0.5*f[29]*hamil[30]+0.5*hamil[29]*f[30]+0.5*f[15]*hamil[28]+0.5*hamil[15]*f[28]+0.5*f[14]*hamil[27]+0.5*hamil[14]*f[27]+0.5*f[11]*hamil[26]+0.5*hamil[11]*f[26]+0.5*f[23]*hamil[24]+0.5*hamil[23]*f[24]+0.5*f[21]*hamil[22]+0.5*hamil[21]*f[22]+0.5*f[5]*hamil[20]+0.5*hamil[5]*f[20]+0.5*f[18]*hamil[19]+0.5*hamil[18]*f[19]+0.5*f[4]*hamil[17]+0.5*hamil[4]*f[17]+0.5*f[3]*hamil[16]+0.5*hamil[3]*f[16]+0.5*f[12]*hamil[13]+0.5*hamil[12]*f[13]+0.5*f[9]*hamil[10]+0.5*hamil[9]*f[10]+0.5*f[7]*hamil[8]+0.5*hamil[7]*f[8]+0.5*f[0]*hamil[6]+0.5*hamil[0]*f[6]+0.5*f[1]*hamil[2]+0.5*hamil[1]*f[2])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_phase_int_five_moments_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]*0.03125; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dv11 = 2.0/dxv[3]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  double dv12 = 2.0/dxv[4]; 
  const double *jacob_vz = &jacob_vel[4]; 
  const double jacob_vz_inv = 1.0/jacob_vz[0]; 
  out[0] += 5.656854249492382*f[0]*volFact; 
  out[1] += (1.7320508075688772*f[28]*hamil[31]+1.7320508075688772*f[24]*hamil[30]+1.7320508075688772*f[23]*hamil[29]+1.7320508075688772*f[20]*hamil[27]+1.7320508075688772*f[17]*hamil[26]+1.7320508075688772*f[15]*hamil[25]+1.7320508075688772*f[13]*hamil[22]+1.7320508075688772*f[12]*hamil[21]+1.7320508075688772*f[10]*hamil[19]+1.7320508075688772*f[9]*hamil[18]+1.7320508075688772*f[6]*hamil[16]+1.7320508075688772*f[5]*hamil[14]+1.7320508075688772*f[4]*hamil[11]+1.7320508075688772*f[2]*hamil[8]+1.7320508075688772*f[1]*hamil[7]+1.7320508075688772*f[0]*hamil[3])*dv10*jacob_vx_inv*volFact; 
  out[2] += (1.7320508075688772*f[27]*hamil[31]+1.7320508075688772*f[22]*hamil[30]+1.7320508075688772*f[21]*hamil[29]+1.7320508075688772*f[20]*hamil[28]+1.7320508075688772*f[16]*hamil[26]+1.7320508075688772*f[14]*hamil[25]+1.7320508075688772*f[13]*hamil[24]+1.7320508075688772*f[12]*hamil[23]+1.7320508075688772*f[8]*hamil[19]+1.7320508075688772*f[7]*hamil[18]+1.7320508075688772*f[6]*hamil[17]+1.7320508075688772*f[5]*hamil[15]+1.7320508075688772*f[3]*hamil[11]+1.7320508075688772*f[2]*hamil[10]+1.7320508075688772*f[1]*hamil[9]+1.7320508075688772*f[0]*hamil[4])*dv11*jacob_vy_inv*volFact; 
  out[3] += (1.7320508075688772*f[26]*hamil[31]+1.7320508075688772*f[19]*hamil[30]+1.7320508075688772*f[18]*hamil[29]+1.7320508075688772*f[17]*hamil[28]+1.7320508075688772*f[16]*hamil[27]+1.7320508075688772*f[11]*hamil[25]+1.7320508075688772*f[10]*hamil[24]+1.7320508075688772*f[9]*hamil[23]+1.7320508075688772*f[8]*hamil[22]+1.7320508075688772*f[7]*hamil[21]+1.7320508075688772*f[6]*hamil[20]+1.7320508075688772*f[4]*hamil[15]+1.7320508075688772*f[3]*hamil[14]+1.7320508075688772*f[2]*hamil[13]+1.7320508075688772*f[1]*hamil[12]+1.7320508075688772*f[0]*hamil[5])*dv12*jacob_vz_inv*volFact; 
  out[4] += (f[31]*hamil[31]+f[30]*hamil[30]+f[29]*hamil[29]+f[28]*hamil[28]+f[27]*hamil[27]+f[26]*hamil[26]+f[25]*hamil[25]+f[24]*hamil[24]+f[23]*hamil[23]+f[22]*hamil[22]+f[21]*hamil[21]+f[20]*hamil[20]+f[19]*hamil[19]+f[18]*hamil[18]+f[17]*hamil[17]+f[16]*hamil[16]+f[15]*hamil[15]+f[14]*hamil[14]+f[13]*hamil[13]+f[12]*hamil[12]+f[11]*hamil[11]+f[10]*hamil[10]+f[9]*hamil[9]+f[8]*hamil[8]+f[7]*hamil[7]+f[6]*hamil[6]+f[5]*hamil[5]+f[4]*hamil[4]+f[3]*hamil[3]+f[2]*hamil[2]+f[1]*hamil[1]+f[0]*hamil[0])*volFact; 
} 
