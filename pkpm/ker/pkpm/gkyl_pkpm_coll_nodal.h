#pragma once 
// Generated (gkylcas/pkpm/pkpm/pkpm-coll-nodal.mac): conversions between the modal expansions and 
// the configuration-nodal, velocity-modal form of the collision kernels. The nodes are the 
// 2^cdim interior Gauss-Legendre points, node q having bit d set when x_d > 0. 
#include <gkyl_util.h> 

// The p=1 configuration-space expansion c at the 2 nodes.
GKYL_CU_DH static inline void
pkpm_coll_conf_to_nodal_1x(const double *c, double* GKYL_RESTRICT n)
{
  const double t0[2] = {c[0], c[1]};
  const double t1[2] = {t0[0]-t0[1], t0[0]+t0[1]};
  n[0] = 0.7071067811865475*t1[0];
  n[1] = 0.7071067811865475*t1[1];
}

// The p=1 expansion c of the 2 node values n.
GKYL_CU_DH static inline void
pkpm_coll_nodal_to_conf_1x(const double *n, double* GKYL_RESTRICT c)
{
  const double t0[2] = {n[0], n[1]};
  const double t1[2] = {t0[0]+t0[1], t0[1]-t0[0]};
  c[0] = 0.7071067811865475*t1[0];
  c[1] = 0.7071067811865475*t1[1];
}

// The 8 modal blocks f (6 coefficients each) at the 2 nodes: g[(q*8 + b)*3 + j].
GKYL_CU_DH static inline void
pkpm_coll_phase_to_nodal_1x1v(const double *f, double* GKYL_RESTRICT g)
{
  for (int b = 0; b < 8; ++b) {
    const double *fb = &f[b*6];
    const double t0[6] = {fb[0], fb[1], fb[2], fb[3], fb[4], fb[5]};
    const double t1[6] = {t0[0]-t0[1], t0[0]+t0[1], t0[2]-t0[3], t0[2]+t0[3], t0[4]-t0[5], t0[4]+t0[5]};
    g[0 + 3*b] = 0.7071067811865475*t1[0];
    g[24 + 3*b] = 0.7071067811865475*t1[1];
    g[1 + 3*b] = 0.7071067811865475*t1[2];
    g[25 + 3*b] = 0.7071067811865475*t1[3];
    g[2 + 3*b] = 0.7071067811865475*t1[4];
    g[26 + 3*b] = 0.7071067811865475*t1[5];
  }
}

// Accumulate the nodal increment g[(q*8 + b)*3 + j] onto the 8 modal blocks out.
GKYL_CU_DH static inline void
pkpm_coll_nodal_to_phase_accumulate_1x1v(const double *g, double* GKYL_RESTRICT out)
{
  for (int b = 0; b < 8; ++b) {
    double *outb = &out[b*6];
    const double t0[6] = {g[0 + 3*b], g[24 + 3*b], g[1 + 3*b], g[25 + 3*b], g[2 + 3*b], g[26 + 3*b]};
    const double t1[6] = {t0[0]+t0[1], t0[1]-t0[0], t0[2]+t0[3], t0[3]-t0[2], t0[4]+t0[5], t0[5]-t0[4]};
    outb[0] += 0.7071067811865475*t1[0];
    outb[1] += 0.7071067811865475*t1[1];
    outb[2] += 0.7071067811865475*t1[2];
    outb[3] += 0.7071067811865475*t1[3];
    outb[4] += 0.7071067811865475*t1[4];
    outb[5] += 0.7071067811865475*t1[5];
  }
}

// The p=1 configuration-space expansion c at the 4 nodes.
GKYL_CU_DH static inline void
pkpm_coll_conf_to_nodal_2x(const double *c, double* GKYL_RESTRICT n)
{
  const double t0[4] = {c[0], c[1], c[2], c[3]};
  const double t1[4] = {t0[0]-t0[1], t0[0]+t0[1], t0[2]-t0[3], t0[2]+t0[3]};
  const double t2[4] = {t1[0]-t1[2], t1[1]-t1[3], t1[0]+t1[2], t1[1]+t1[3]};
  n[0] = 0.5*t2[0];
  n[1] = 0.5*t2[1];
  n[2] = 0.5*t2[2];
  n[3] = 0.5*t2[3];
}

// The p=1 expansion c of the 4 node values n.
GKYL_CU_DH static inline void
pkpm_coll_nodal_to_conf_2x(const double *n, double* GKYL_RESTRICT c)
{
  const double t0[4] = {n[0], n[1], n[2], n[3]};
  const double t1[4] = {t0[0]+t0[1], t0[1]-t0[0], t0[2]+t0[3], t0[3]-t0[2]};
  const double t2[4] = {t1[0]+t1[2], t1[1]+t1[3], t1[2]-t1[0], t1[3]-t1[1]};
  c[0] = 0.5*t2[0];
  c[1] = 0.5*t2[1];
  c[2] = 0.5*t2[2];
  c[3] = 0.5*t2[3];
}

// The 8 modal blocks f (12 coefficients each) at the 4 nodes: g[(q*8 + b)*3 + j].
GKYL_CU_DH static inline void
pkpm_coll_phase_to_nodal_2x1v(const double *f, double* GKYL_RESTRICT g)
{
  for (int b = 0; b < 8; ++b) {
    const double *fb = &f[b*12];
    const double t0[12] = {fb[0], fb[1], fb[2], fb[4], fb[3], fb[5], fb[6], fb[8], fb[7], fb[9], fb[10], fb[11]};
    const double t1[12] = {t0[0]-t0[1], t0[0]+t0[1], t0[2]-t0[3], t0[2]+t0[3], t0[4]-t0[5], t0[4]+t0[5], t0[6]-t0[7], t0[6]+t0[7], t0[8]-t0[9], t0[8]+t0[9], t0[10]-t0[11], t0[10]+t0[11]};
    const double t2[12] = {t1[0]-t1[2], t1[1]-t1[3], t1[0]+t1[2], t1[1]+t1[3], t1[4]-t1[6], t1[5]-t1[7], t1[4]+t1[6], t1[5]+t1[7], t1[8]-t1[10], t1[9]-t1[11], t1[8]+t1[10], t1[9]+t1[11]};
    g[0 + 3*b] = 0.5*t2[0];
    g[24 + 3*b] = 0.5*t2[1];
    g[48 + 3*b] = 0.5*t2[2];
    g[72 + 3*b] = 0.5*t2[3];
    g[1 + 3*b] = 0.5*t2[4];
    g[25 + 3*b] = 0.5*t2[5];
    g[49 + 3*b] = 0.5*t2[6];
    g[73 + 3*b] = 0.5*t2[7];
    g[2 + 3*b] = 0.5*t2[8];
    g[26 + 3*b] = 0.5*t2[9];
    g[50 + 3*b] = 0.5*t2[10];
    g[74 + 3*b] = 0.5*t2[11];
  }
}

// Accumulate the nodal increment g[(q*8 + b)*3 + j] onto the 8 modal blocks out.
GKYL_CU_DH static inline void
pkpm_coll_nodal_to_phase_accumulate_2x1v(const double *g, double* GKYL_RESTRICT out)
{
  for (int b = 0; b < 8; ++b) {
    double *outb = &out[b*12];
    const double t0[12] = {g[0 + 3*b], g[24 + 3*b], g[48 + 3*b], g[72 + 3*b], g[1 + 3*b], g[25 + 3*b], g[49 + 3*b], g[73 + 3*b], g[2 + 3*b], g[26 + 3*b], g[50 + 3*b], g[74 + 3*b]};
    const double t1[12] = {t0[0]+t0[1], t0[1]-t0[0], t0[2]+t0[3], t0[3]-t0[2], t0[4]+t0[5], t0[5]-t0[4], t0[6]+t0[7], t0[7]-t0[6], t0[8]+t0[9], t0[9]-t0[8], t0[10]+t0[11], t0[11]-t0[10]};
    const double t2[12] = {t1[0]+t1[2], t1[1]+t1[3], t1[2]-t1[0], t1[3]-t1[1], t1[4]+t1[6], t1[5]+t1[7], t1[6]-t1[4], t1[7]-t1[5], t1[8]+t1[10], t1[9]+t1[11], t1[10]-t1[8], t1[11]-t1[9]};
    outb[0] += 0.5*t2[0];
    outb[1] += 0.5*t2[1];
    outb[2] += 0.5*t2[2];
    outb[4] += 0.5*t2[3];
    outb[3] += 0.5*t2[4];
    outb[5] += 0.5*t2[5];
    outb[6] += 0.5*t2[6];
    outb[8] += 0.5*t2[7];
    outb[7] += 0.5*t2[8];
    outb[9] += 0.5*t2[9];
    outb[10] += 0.5*t2[10];
    outb[11] += 0.5*t2[11];
  }
}

// The p=1 configuration-space expansion c at the 8 nodes.
GKYL_CU_DH static inline void
pkpm_coll_conf_to_nodal_3x(const double *c, double* GKYL_RESTRICT n)
{
  const double t0[8] = {c[0], c[1], c[2], c[4], c[3], c[5], c[6], c[7]};
  const double t1[8] = {t0[0]-t0[1], t0[0]+t0[1], t0[2]-t0[3], t0[2]+t0[3], t0[4]-t0[5], t0[4]+t0[5], t0[6]-t0[7], t0[6]+t0[7]};
  const double t2[8] = {t1[0]-t1[2], t1[1]-t1[3], t1[0]+t1[2], t1[1]+t1[3], t1[4]-t1[6], t1[5]-t1[7], t1[4]+t1[6], t1[5]+t1[7]};
  const double t3[8] = {t2[0]-t2[4], t2[1]-t2[5], t2[2]-t2[6], t2[3]-t2[7], t2[0]+t2[4], t2[1]+t2[5], t2[2]+t2[6], t2[3]+t2[7]};
  n[0] = 0.3535533905932737*t3[0];
  n[1] = 0.3535533905932737*t3[1];
  n[2] = 0.3535533905932737*t3[2];
  n[3] = 0.3535533905932737*t3[3];
  n[4] = 0.3535533905932737*t3[4];
  n[5] = 0.3535533905932737*t3[5];
  n[6] = 0.3535533905932737*t3[6];
  n[7] = 0.3535533905932737*t3[7];
}

// The p=1 expansion c of the 8 node values n.
GKYL_CU_DH static inline void
pkpm_coll_nodal_to_conf_3x(const double *n, double* GKYL_RESTRICT c)
{
  const double t0[8] = {n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7]};
  const double t1[8] = {t0[0]+t0[1], t0[1]-t0[0], t0[2]+t0[3], t0[3]-t0[2], t0[4]+t0[5], t0[5]-t0[4], t0[6]+t0[7], t0[7]-t0[6]};
  const double t2[8] = {t1[0]+t1[2], t1[1]+t1[3], t1[2]-t1[0], t1[3]-t1[1], t1[4]+t1[6], t1[5]+t1[7], t1[6]-t1[4], t1[7]-t1[5]};
  const double t3[8] = {t2[0]+t2[4], t2[1]+t2[5], t2[2]+t2[6], t2[3]+t2[7], t2[4]-t2[0], t2[5]-t2[1], t2[6]-t2[2], t2[7]-t2[3]};
  c[0] = 0.3535533905932737*t3[0];
  c[1] = 0.3535533905932737*t3[1];
  c[2] = 0.3535533905932737*t3[2];
  c[4] = 0.3535533905932737*t3[3];
  c[3] = 0.3535533905932737*t3[4];
  c[5] = 0.3535533905932737*t3[5];
  c[6] = 0.3535533905932737*t3[6];
  c[7] = 0.3535533905932737*t3[7];
}

// The 8 modal blocks f (24 coefficients each) at the 8 nodes: g[(q*8 + b)*3 + j].
GKYL_CU_DH static inline void
pkpm_coll_phase_to_nodal_3x1v(const double *f, double* GKYL_RESTRICT g)
{
  for (int b = 0; b < 8; ++b) {
    const double *fb = &f[b*24];
    const double t0[24] = {fb[0], fb[1], fb[2], fb[5], fb[3], fb[6], fb[7], fb[12], fb[4], fb[8], fb[9], fb[13], fb[10], fb[14], fb[15], fb[19], fb[11], fb[16], fb[17], fb[20], fb[18], fb[21], fb[22], fb[23]};
    const double t1[24] = {t0[0]-t0[1], t0[0]+t0[1], t0[2]-t0[3], t0[2]+t0[3], t0[4]-t0[5], t0[4]+t0[5], t0[6]-t0[7], t0[6]+t0[7], t0[8]-t0[9], t0[8]+t0[9], t0[10]-t0[11], t0[10]+t0[11], t0[12]-t0[13], t0[12]+t0[13], t0[14]-t0[15], t0[14]+t0[15], t0[16]-t0[17], t0[16]+t0[17], t0[18]-t0[19], t0[18]+t0[19], t0[20]-t0[21], t0[20]+t0[21], t0[22]-t0[23], t0[22]+t0[23]};
    const double t2[24] = {t1[0]-t1[2], t1[1]-t1[3], t1[0]+t1[2], t1[1]+t1[3], t1[4]-t1[6], t1[5]-t1[7], t1[4]+t1[6], t1[5]+t1[7], t1[8]-t1[10], t1[9]-t1[11], t1[8]+t1[10], t1[9]+t1[11], t1[12]-t1[14], t1[13]-t1[15], t1[12]+t1[14], t1[13]+t1[15], t1[16]-t1[18], t1[17]-t1[19], t1[16]+t1[18], t1[17]+t1[19], t1[20]-t1[22], t1[21]-t1[23], t1[20]+t1[22], t1[21]+t1[23]};
    const double t3[24] = {t2[0]-t2[4], t2[1]-t2[5], t2[2]-t2[6], t2[3]-t2[7], t2[0]+t2[4], t2[1]+t2[5], t2[2]+t2[6], t2[3]+t2[7], t2[8]-t2[12], t2[9]-t2[13], t2[10]-t2[14], t2[11]-t2[15], t2[8]+t2[12], t2[9]+t2[13], t2[10]+t2[14], t2[11]+t2[15], t2[16]-t2[20], t2[17]-t2[21], t2[18]-t2[22], t2[19]-t2[23], t2[16]+t2[20], t2[17]+t2[21], t2[18]+t2[22], t2[19]+t2[23]};
    g[0 + 3*b] = 0.3535533905932737*t3[0];
    g[24 + 3*b] = 0.3535533905932737*t3[1];
    g[48 + 3*b] = 0.3535533905932737*t3[2];
    g[72 + 3*b] = 0.3535533905932737*t3[3];
    g[96 + 3*b] = 0.3535533905932737*t3[4];
    g[120 + 3*b] = 0.3535533905932737*t3[5];
    g[144 + 3*b] = 0.3535533905932737*t3[6];
    g[168 + 3*b] = 0.3535533905932737*t3[7];
    g[1 + 3*b] = 0.3535533905932737*t3[8];
    g[25 + 3*b] = 0.3535533905932737*t3[9];
    g[49 + 3*b] = 0.3535533905932737*t3[10];
    g[73 + 3*b] = 0.3535533905932737*t3[11];
    g[97 + 3*b] = 0.3535533905932737*t3[12];
    g[121 + 3*b] = 0.3535533905932737*t3[13];
    g[145 + 3*b] = 0.3535533905932737*t3[14];
    g[169 + 3*b] = 0.3535533905932737*t3[15];
    g[2 + 3*b] = 0.3535533905932737*t3[16];
    g[26 + 3*b] = 0.3535533905932737*t3[17];
    g[50 + 3*b] = 0.3535533905932737*t3[18];
    g[74 + 3*b] = 0.3535533905932737*t3[19];
    g[98 + 3*b] = 0.3535533905932737*t3[20];
    g[122 + 3*b] = 0.3535533905932737*t3[21];
    g[146 + 3*b] = 0.3535533905932737*t3[22];
    g[170 + 3*b] = 0.3535533905932737*t3[23];
  }
}

// Accumulate the nodal increment g[(q*8 + b)*3 + j] onto the 8 modal blocks out.
GKYL_CU_DH static inline void
pkpm_coll_nodal_to_phase_accumulate_3x1v(const double *g, double* GKYL_RESTRICT out)
{
  for (int b = 0; b < 8; ++b) {
    double *outb = &out[b*24];
    const double t0[24] = {g[0 + 3*b], g[24 + 3*b], g[48 + 3*b], g[72 + 3*b], g[96 + 3*b], g[120 + 3*b], g[144 + 3*b], g[168 + 3*b], g[1 + 3*b], g[25 + 3*b], g[49 + 3*b], g[73 + 3*b], g[97 + 3*b], g[121 + 3*b], g[145 + 3*b], g[169 + 3*b], g[2 + 3*b], g[26 + 3*b], g[50 + 3*b], g[74 + 3*b], g[98 + 3*b], g[122 + 3*b], g[146 + 3*b], g[170 + 3*b]};
    const double t1[24] = {t0[0]+t0[1], t0[1]-t0[0], t0[2]+t0[3], t0[3]-t0[2], t0[4]+t0[5], t0[5]-t0[4], t0[6]+t0[7], t0[7]-t0[6], t0[8]+t0[9], t0[9]-t0[8], t0[10]+t0[11], t0[11]-t0[10], t0[12]+t0[13], t0[13]-t0[12], t0[14]+t0[15], t0[15]-t0[14], t0[16]+t0[17], t0[17]-t0[16], t0[18]+t0[19], t0[19]-t0[18], t0[20]+t0[21], t0[21]-t0[20], t0[22]+t0[23], t0[23]-t0[22]};
    const double t2[24] = {t1[0]+t1[2], t1[1]+t1[3], t1[2]-t1[0], t1[3]-t1[1], t1[4]+t1[6], t1[5]+t1[7], t1[6]-t1[4], t1[7]-t1[5], t1[8]+t1[10], t1[9]+t1[11], t1[10]-t1[8], t1[11]-t1[9], t1[12]+t1[14], t1[13]+t1[15], t1[14]-t1[12], t1[15]-t1[13], t1[16]+t1[18], t1[17]+t1[19], t1[18]-t1[16], t1[19]-t1[17], t1[20]+t1[22], t1[21]+t1[23], t1[22]-t1[20], t1[23]-t1[21]};
    const double t3[24] = {t2[0]+t2[4], t2[1]+t2[5], t2[2]+t2[6], t2[3]+t2[7], t2[4]-t2[0], t2[5]-t2[1], t2[6]-t2[2], t2[7]-t2[3], t2[8]+t2[12], t2[9]+t2[13], t2[10]+t2[14], t2[11]+t2[15], t2[12]-t2[8], t2[13]-t2[9], t2[14]-t2[10], t2[15]-t2[11], t2[16]+t2[20], t2[17]+t2[21], t2[18]+t2[22], t2[19]+t2[23], t2[20]-t2[16], t2[21]-t2[17], t2[22]-t2[18], t2[23]-t2[19]};
    outb[0] += 0.3535533905932737*t3[0];
    outb[1] += 0.3535533905932737*t3[1];
    outb[2] += 0.3535533905932737*t3[2];
    outb[5] += 0.3535533905932737*t3[3];
    outb[3] += 0.3535533905932737*t3[4];
    outb[6] += 0.3535533905932737*t3[5];
    outb[7] += 0.3535533905932737*t3[6];
    outb[12] += 0.3535533905932737*t3[7];
    outb[4] += 0.3535533905932737*t3[8];
    outb[8] += 0.3535533905932737*t3[9];
    outb[9] += 0.3535533905932737*t3[10];
    outb[13] += 0.3535533905932737*t3[11];
    outb[10] += 0.3535533905932737*t3[12];
    outb[14] += 0.3535533905932737*t3[13];
    outb[15] += 0.3535533905932737*t3[14];
    outb[19] += 0.3535533905932737*t3[15];
    outb[11] += 0.3535533905932737*t3[16];
    outb[16] += 0.3535533905932737*t3[17];
    outb[17] += 0.3535533905932737*t3[18];
    outb[20] += 0.3535533905932737*t3[19];
    outb[18] += 0.3535533905932737*t3[20];
    outb[21] += 0.3535533905932737*t3[21];
    outb[22] += 0.3535533905932737*t3[22];
    outb[23] += 0.3535533905932737*t3[23];
  }
}

