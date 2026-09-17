#pragma once

// Private header, not for direct use in user code

#include <math.h>
#include <gkyl_array.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_wv_embed_geo.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

#include <gkyl_kann_net.h>
#include <gkyl_knutils.h>

struct wv_ten_moment {
  struct gkyl_wv_eqn eqn; // Base object.
  double k0; // Closure parameter.
  bool use_grad_closure; // Should we use gradient-based closure?
  bool use_nn_closure; // Should we use neural network-based closure?
  int poly_order; // Polynomial order of learned DG coefficients.
  struct gkyl_kann_net* ann; // Neural network architecture.
};

/**
 * Free Ten moment eqn object.
 *
 * @param ref Reference counter for Ten moment eqn
 */
void gkyl_ten_moment_free(const struct gkyl_ref_count *ref);

GKYL_CU_DH
static inline double 
sq(double x) 
{ 
  return x * x; 
}

GKYL_CU_DH
static inline void
cons_to_riem(const struct gkyl_wv_eqn *eqn,
  const double *qstate, const double *qin, double *wout)
{
  // TODO: this should use proper L matrix
  for (int i=0; i<10; ++i)
    wout[i] = qin[i];
}

GKYL_CU_DH
static inline void
riem_to_cons(const struct gkyl_wv_eqn *eqn,
  const double *qstate, const double *win, double *qout)
{
  // TODO: this should use proper L matrix
  for (int i=0; i<10; ++i)
    qout[i] = win[i];
}

/* Multiply by phi prime */
GKYL_CU_DH
static inline void 
mulByPhiPrime(double p0, double u1, double u2, double u3, const double w[10], double out[10]) 
{ 
  out[0] = w[0]; 
  out[1] = w[0]*u1+w[1]*p0; 
  out[2] = w[0]*u2+w[2]*p0; 
  out[3] = w[0]*u3+w[3]*p0; 
  out[4] = w[0]*sq(u1)+2*w[1]*p0*u1+w[4]; 
  out[5] = w[0]*u1*u2+w[1]*p0*u2+w[2]*p0*u1+w[5]; 
  out[6] = w[0]*u1*u3+w[1]*p0*u3+w[3]*p0*u1+w[6]; 
  out[7] = w[0]*sq(u2)+2*w[2]*p0*u2+w[7]; 
  out[8] = w[0]*u2*u3+w[2]*p0*u3+w[3]*p0*u2+w[8]; 
  out[9] = w[0]*sq(u3)+2*w[3]*p0*u3+w[9]; 
} 

/**
 * Computes the primitive variables given the conserved variables.
 * 
 * @param q Conserved variables
 * @param out Primitive variables
 */
GKYL_CU_DH
static inline void
gkyl_ten_moment_primitive(const double q[10], double out[10])
{
  out[0] = q[0]; 
  out[1] = q[1]/q[0]; 
  out[2] = q[2]/q[0]; 
  out[3] = q[3]/q[0]; 
  out[4] = q[4]-(q[1]*q[1])/q[0]; 
  out[5] = q[5]-(q[1]*q[2])/q[0]; 
  out[6] = q[6]-(q[1]*q[3])/q[0]; 
  out[7] = q[7]-(q[2]*q[2])/q[0]; 
  out[8] = q[8]-(q[2]*q[3])/q[0]; 
  out[9] = q[9]-(q[3]*q[3])/q[0];
}

/**
 * Computes the diagonal components of the pressure tensor
 * 
 * @param q Conserved variables
 * @param out [Pxx, Pyy, Pzz]
 */
GKYL_CU_DH
static inline void
gkyl_ten_moment_diag_pressure(const double q[10], double out[3])
{
  out[0] = q[4]-(q[1]*q[1])/q[0]; // pxx
  out[1] = q[7]-(q[2]*q[2])/q[0]; // pyy
  out[2] = q[9]-(q[3]*q[3])/q[0]; // pzz
}

/**
 * Compute maximum absolute speed.
 * 
 * @param dir Direction 
 * @param q Conserved variables
 * @return Maximum absolute speed for given q
 */
GKYL_CU_DH
static inline double
gkyl_ten_moment_max_abs_speed(const double q[10])
{
  double u = q[1]/q[0];
  double p11 = q[4] - q[0]*u*u;
  return fabs(u) + sqrt(3.0*p11/q[0]);
}

/**
 * Compute flux in direction 'dir'.
 * 
 * @param dir Direction to compute flux
 * @param Conserved variables
 * @param flux On output, the flux in direction 'dir'
 */
GKYL_CU_DH
static inline void
gkyl_ten_moment_flux(const double q[10], double flux[10])
{
  double v[10];
  gkyl_ten_moment_primitive(q, v);

  flux[0] = q[1]; // rho*u
  flux[1] = q[4]; // Pxx
  flux[2] = q[5]; // Pxy
  flux[3] = q[6]; // Pxz
  flux[4] = v[0]*v[1]*v[1]*v[1] + 3*v[1]*v[4]; // rho u^3 + 3*u*Pxx
  flux[5] = v[0]*v[1]*v[1]*v[2] + 2*v[1]*v[5] + v[2]*v[4]; // rho*u^2*v + 2*u*Pxy + v*Pxx
  flux[6] = v[0]*v[1]*v[1]*v[3] + 2*v[1]*v[6] + v[3]*v[4]; // rho*u^2*w + 2*u*Pxz + w*Pxx
  flux[7] = v[0]*v[1]*v[2]*v[2] + 2*v[2]*v[5] + v[1]*v[7]; // rho*u*v^2 + 2*v*Pxy + u*Pyy
  flux[8] = v[0]*v[1]*v[2]*v[3] + v[1]*v[8] + v[2]*v[6] + v[3]*v[5]; // rho*u*v*w + u*Pyz + v*Pxz + w*Pxy
  flux[9] = v[0]*v[1]*v[3]*v[3] + 2*v[3]*v[6] + v[1]*v[9]; // rho*u*w^2 + 2*w*Pxz + u*Pzz
}

// Ten moment perfectly reflecting wall
GKYL_CU_DH
static void
ten_moment_wall(const struct gkyl_wv_eqn* eqn, double t, int nc, const double *skin, double * GKYL_RESTRICT ghost, void *ctx)
{
  // copy density and Pxx, Pyy, and Pzz
  ghost[0] = skin[0];
  ghost[4] = skin[4];
  ghost[7] = skin[7];
  ghost[9] = skin[9];

  // zero-normal for momentum
  ghost[1] = -skin[1];
  ghost[2] = skin[2];
  ghost[3] = skin[3];

  // zero-tangent for off-diagonal components of pressure tensor
  ghost[8] = skin[8];
  ghost[6] = -skin[6];
  ghost[5] = -skin[5];
}

GKYL_CU_DH
static inline void
rot_to_local(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm,
  const double* GKYL_RESTRICT qglobal, double* GKYL_RESTRICT qlocal)
{
  // Mass density is a scalar
  qlocal[0] = qglobal[0];
  // Rotate momentum to local coordinates
  qlocal[1] = qglobal[1]*norm[0] + qglobal[2]*norm[1] + qglobal[3]*norm[2];
  qlocal[2] = qglobal[1]*tau1[0] + qglobal[2]*tau1[1] + qglobal[3]*tau1[2];
  qlocal[3] = qglobal[1]*tau2[0] + qglobal[2]*tau2[1] + qglobal[3]*tau2[2];

  // temp arrays to store rotated column vectors
  double r1[3], r2[3], r3[3];
  r1[0] = qglobal[4]*norm[0] + qglobal[5]*norm[1] + qglobal[6]*norm[2];
  r1[1] = qglobal[4]*tau1[0] + qglobal[5]*tau1[1] + qglobal[6]*tau1[2];
  r1[2] = qglobal[4]*tau2[0] + qglobal[5]*tau2[1] + qglobal[6]*tau2[2];

  r2[0] = qglobal[5]*norm[0] + qglobal[7]*norm[1] + qglobal[8]*norm[2];
  r2[1] = qglobal[5]*tau1[0] + qglobal[7]*tau1[1] + qglobal[8]*tau1[2];
  r2[2] = qglobal[5]*tau2[0] + qglobal[7]*tau2[1] + qglobal[8]*tau2[2];

  r3[0] = qglobal[6]*norm[0] + qglobal[8]*norm[1] + qglobal[9]*norm[2];
  r3[1] = qglobal[6]*tau1[0] + qglobal[8]*tau1[1] + qglobal[9]*tau1[2];
  r3[2] = qglobal[6]*tau2[0] + qglobal[8]*tau2[1] + qglobal[9]*tau2[2];

  // temp arrays to store rotated row vectors
  double v1[3], v2[3], v3[3];
  v1[0] = r1[0]*norm[0] + r2[0]*norm[1] + r3[0]*norm[2];
  v1[1] = r1[0]*tau1[0] + r2[0]*tau1[1] + r3[0]*tau1[2];
  v1[2] = r1[0]*tau2[0] + r2[0]*tau2[1] + r3[0]*tau2[2];

  v2[0] = r1[1]*norm[0] + r2[1]*norm[1] + r3[1]*norm[2];
  v2[1] = r1[1]*tau1[0] + r2[1]*tau1[1] + r3[1]*tau1[2];
  v2[2] = r1[1]*tau2[0] + r2[1]*tau2[1] + r3[1]*tau2[2]; 

  v3[0] = r1[2]*norm[0] + r2[2]*norm[1] + r3[2]*norm[2];
  v3[1] = r1[2]*tau1[0] + r2[2]*tau1[1] + r3[2]*tau1[2];
  v3[2] = r1[2]*tau2[0] + r2[2]*tau2[1] + r3[2]*tau2[2];

  qlocal[4] = v1[0];
  qlocal[5] = v1[1];
  qlocal[6] = v1[2];
  qlocal[7] = v2[1];
  qlocal[8] = v2[2];
  qlocal[9] = v3[2];
}

GKYL_CU_DH
static inline void
rot_to_global(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm,
  const double* GKYL_RESTRICT qlocal, double* GKYL_RESTRICT qglobal)
{
 
  // Mass density is a scalar
  qglobal[0] = qlocal[0];
  // Rotate momentum back to global coordinates
  qglobal[1] = qlocal[1]*norm[0] + qlocal[2]*tau1[0] + qlocal[3]*tau2[0];
  qglobal[2] = qlocal[1]*norm[1] + qlocal[2]*tau1[1] + qlocal[3]*tau2[1];
  qglobal[3] = qlocal[1]*norm[2] + qlocal[2]*tau1[2] + qlocal[3]*tau2[2];

  // temp arrays to store rotated column vectors
  double r1[3], r2[3], r3[3];
  r1[0] = qlocal[4]*norm[0] + qlocal[5]*tau1[0] + qlocal[6]*tau2[0];
  r1[1] = qlocal[4]*norm[1] + qlocal[5]*tau1[1] + qlocal[6]*tau2[1];
  r1[2] = qlocal[4]*norm[2] + qlocal[5]*tau1[2] + qlocal[6]*tau2[2];

  r2[0] = qlocal[5]*norm[0] + qlocal[7]*tau1[0] + qlocal[8]*tau2[0];
  r2[1] = qlocal[5]*norm[1] + qlocal[7]*tau1[1] + qlocal[8]*tau2[1];
  r2[2] = qlocal[5]*norm[2] + qlocal[7]*tau1[2] + qlocal[8]*tau2[2];

  r3[0] = qlocal[6]*norm[0] + qlocal[8]*tau1[0] + qlocal[9]*tau2[0];
  r3[1] = qlocal[6]*norm[1] + qlocal[8]*tau1[1] + qlocal[9]*tau2[1];
  r3[2] = qlocal[6]*norm[2] + qlocal[8]*tau1[2] + qlocal[9]*tau2[2];

  // temp arrays to store rotated row vectors
  double v1[3], v2[3], v3[3];
  v1[0] = r1[0]*norm[0] + r2[0]*tau1[0] + r3[0]*tau2[0];
  v1[1] = r1[0]*norm[1] + r2[0]*tau1[1] + r3[0]*tau2[1];
  v1[2] = r1[0]*norm[2] + r2[0]*tau1[2] + r3[0]*tau2[2];

  v2[0] = r1[1]*norm[0] + r2[1]*tau1[0] + r3[1]*tau2[0];
  v2[1] = r1[1]*norm[1] + r2[1]*tau1[1] + r3[1]*tau2[1];
  v2[2] = r1[1]*norm[2] + r2[1]*tau1[2] + r3[1]*tau2[2]; 

  v3[0] = r1[2]*norm[0] + r2[2]*tau1[0] + r3[2]*tau2[0];
  v3[1] = r1[2]*norm[1] + r2[2]*tau1[1] + r3[2]*tau2[1];
  v3[2] = r1[2]*norm[2] + r2[2]*tau1[2] + r3[2]*tau2[2];

  // Rotate pressure tensor back to local coordinates
  qglobal[4] = v1[0];
  qglobal[5] = v1[1];
  qglobal[6] = v1[2];
  qglobal[7] = v2[1];
  qglobal[8] = v2[2];
  qglobal[9] = v3[2];
}

// Waves and speeds using Roe averaging
GKYL_CU_DH
static double
wave_roe(const struct gkyl_wv_eqn *eqn,
  const double *delta, const double *ql, const double *qr, 
  double *waves, double *s)
{
  double vl[10], vr[10];
  gkyl_ten_moment_primitive(ql, vl);
  gkyl_ten_moment_primitive(qr, vr);

  // compute Roe averages
  double sqrl = sqrt(vl[0]), sqrr = sqrt(vr[0]);
  double sqr1 = 1/(sqrl+sqrr);
  
  double p0 = sqrl*sqrr;
  double p2s1 = sq(p0*sqr1);
  
  double u1 = (sqrl*vl[1] + sqrr*vr[1])*sqr1;
  double u2 = (sqrl*vl[2] + sqrr*vr[2])*sqr1;
  double u3 = (sqrl*vl[3] + sqrr*vr[3])*sqr1;
  double p11 = (sqrr*vl[4]+sqrl*vr[4])*sqr1 + 1.0/3.0*p2s1*(vr[1]-vl[1])*(vr[1]-vl[1]);
  double p12 = (sqrr*vl[5]+sqrl*vr[5])*sqr1 + 1.0/3.0*p2s1*(vr[1]-vl[1])*(vr[2]-vl[2]);
  double p13 = (sqrr*vl[6]+sqrl*vr[6])*sqr1 + 1.0/3.0*p2s1*(vr[1]-vl[1])*(vr[3]-vl[3]);
  double p22 = (sqrr*vl[7]+sqrl*vr[7])*sqr1 + 1.0/3.0*p2s1*(vr[2]-vl[2])*(vr[2]-vl[2]);
  double p23 = (sqrr*vl[8]+sqrl*vr[8])*sqr1 + 1.0/3.0*p2s1*(vr[2]-vl[2])*(vr[3]-vl[3]);
  double p33 = (sqrr*vl[9]+sqrl*vr[9])*sqr1 + 1.0/3.0*p2s1*(vr[3]-vl[3])*(vr[3]-vl[3]);

  // for multiplication by phi' we need to use unrotated values
  double v[4];
  v[1] = u1; v[2] = u2; v[3] = u3;

  double phiDelta[10];

  // pre-multiply jump (delta) by phiPrime inverse: we do this as
  // jumps are in conserved variables, while left eigenvectors used
  // below are computed from primitive variables
  phiDelta[0] = delta[0];
  phiDelta[1] = delta[1]/p0-(1.0*delta[0]*u1)/p0; 
  phiDelta[2] = delta[2]/p0-(1.0*delta[0]*u2)/p0; 
  phiDelta[3] = delta[3]/p0-(1.0*delta[0]*u3)/p0; 
  phiDelta[4] = delta[0]*sq(u1)-2.0*delta[1]*u1+delta[4]; 
  phiDelta[5] = delta[0]*u1*u2-1.0*delta[1]*u2-1.0*delta[2]*u1+delta[5]; 
  phiDelta[6] = delta[0]*u1*u3-1.0*delta[1]*u3-1.0*delta[3]*u1+delta[6]; 
  phiDelta[7] = delta[0]*sq(u2)-2.0*delta[2]*u2+delta[7]; 
  phiDelta[8] = delta[0]*u2*u3-1.0*delta[2]*u3-1.0*delta[3]*u2+delta[8]; 
  phiDelta[9] = delta[0]*sq(u3)-2.0*delta[3]*u3+delta[9];

  // predefine some constants
  double p11sq = sq(p11), p12sq = sq(p12), p13sq = sq(p13), p11th = sqrt(p11*p11*p11);
  double sqp0 = sqrt(p0), sqp11 = sqrt(p11);
  
  double leftProj[10];
  // project jumps on left eigenvectors [Gen from Maxima]
  leftProj[0] = (0.5*phiDelta[1]*sqp0*sqp11*p12)/p11sq-(0.5*phiDelta[4]*p12)/p11sq-(0.5*phiDelta[2]*sqp0)/sqp11+(0.5*phiDelta[5])/p11; 
  leftProj[1] = (0.5*phiDelta[1]*sqp0*p13)/p11th-(0.5*phiDelta[4]*p13)/p11sq-(0.5*phiDelta[3]*sqp0)/sqp11+(0.5*phiDelta[6])/p11; 
  leftProj[2] = (-(0.5*phiDelta[1]*sqp0*sqp11*p12)/p11sq)-(0.5*phiDelta[4]*p12)/p11sq+(0.5*phiDelta[2]*sqp0)/sqp11+(0.5*phiDelta[5])/p11; 
  leftProj[3] = (-(0.5*phiDelta[1]*sqp0*p13)/p11th)-(0.5*phiDelta[4]*p13)/p11sq+(0.5*phiDelta[3]*sqp0)/sqp11+(0.5*phiDelta[6])/p11; 
  leftProj[4] = (0.16666666666666666667*phiDelta[4])/p11sq-(0.2886751345948129*phiDelta[1]*sqp0)/p11th; 
  leftProj[5] = (0.2886751345948129*phiDelta[1]*sqp0)/p11th+(0.16666666666666666667*phiDelta[4])/p11sq; 
  leftProj[6] = phiDelta[0]-(0.3333333333333333333*phiDelta[4]*p0)/p11; 
  leftProj[7] = (-(0.3333333333333333333*phiDelta[4]*p11*p22)/p11sq)+(1.333333333333333333*phiDelta[4]*p12sq)/p11sq-(2.0*phiDelta[5]*p12)/p11+phiDelta[7]; 
  leftProj[8] = (-(0.3333333333333333333*phiDelta[4]*p11*p23)/p11sq)+(1.333333333333333333*phiDelta[4]*p12*p13)/p11sq-(1.0*phiDelta[5]*p13)/p11-(1.0*phiDelta[6]*p12)/p11+phiDelta[8]; 
  leftProj[9] = (-(0.3333333333333333333*phiDelta[4]*p11*p33)/p11sq)+(1.333333333333333333*phiDelta[4]*p13sq)/p11sq-(2.0*phiDelta[6]*p13)/p11+phiDelta[9]; 

  // compute waves and speeds
  double wv[10];

  // Wave 1: (ev 1 and 2 are repeated)
  s[0] = u1-sqrt(p11/p0);
  wv[0] = 0.0; 
  wv[1] = 0.0; 
  wv[2] = -(1.0*leftProj[0]*sqp11)/sqp0; 
  wv[3] = -(1.0*leftProj[1]*sqp11)/sqp0; 
  wv[4] = 0.0; 
  wv[5] = leftProj[0]*p11; 
  wv[6] = leftProj[1]*p11; 
  wv[7] = 2.0*leftProj[0]*p12; 
  wv[8] = leftProj[0]*p13+leftProj[1]*p12; 
  wv[9] = 2.0*leftProj[1]*p13;

  mulByPhiPrime(p0, v[1], v[2], v[3], wv, &waves[0]);

  // Wave 2: (ev 3 and 4 are repeated)
  s[1] = u1+sqrt(p11/p0);
  wv[0] = 0.0; 
  wv[1] = 0.0; 
  wv[2] = (leftProj[2]*sqp11)/sqp0; 
  wv[3] = (leftProj[3]*sqp11)/sqp0; 
  wv[4] = 0.0; 
  wv[5] = leftProj[2]*p11; 
  wv[6] = leftProj[3]*p11; 
  wv[7] = 2.0*leftProj[2]*p12; 
  wv[8] = leftProj[2]*p13+leftProj[3]*p12; 
  wv[9] = 2.0*leftProj[3]*p13;

  mulByPhiPrime(p0, v[1], v[2], v[3], wv, &waves[10]);
  
  // Wave 3 (ev 5)
  s[2] = u1-sqrt(3*p11/p0);
  wv[0] = leftProj[4]*p0*p11; 
  wv[1] = -(1.732050807568877*leftProj[4]*p11th)/sqp0; 
  wv[2] = -(1.732050807568877*leftProj[4]*sqp11*p12)/sqp0; 
  wv[3] = -(1.732050807568877*leftProj[4]*sqp11*p13)/sqp0; 
  wv[4] = 3.0*leftProj[4]*p11sq; 
  wv[5] = 3.0*leftProj[4]*p11*p12; 
  wv[6] = 3.0*leftProj[4]*p11*p13; 
  wv[7] = leftProj[4]*p11*p22+2.0*leftProj[4]*p12sq; 
  wv[8] = leftProj[4]*p11*p23+2.0*leftProj[4]*p12*p13; 
  wv[9] = leftProj[4]*p11*p33+2.0*leftProj[4]*p13sq;

  mulByPhiPrime(p0, v[1], v[2], v[3], wv, &waves[20]);

  // Wave 4 (ev 6)
  s[3] = u1+sqrt(3*p11/p0);
  wv[0] = leftProj[5]*p0*p11; 
  wv[1] = (1.732050807568877*leftProj[5]*p11th)/sqp0; 
  wv[2] = (1.732050807568877*leftProj[5]*sqp11*p12)/sqp0; 
  wv[3] = (1.732050807568877*leftProj[5]*sqp11*p13)/sqp0; 
  wv[4] = 3.0*leftProj[5]*p11sq; 
  wv[5] = 3.0*leftProj[5]*p11*p12; 
  wv[6] = 3.0*leftProj[5]*p11*p13; 
  wv[7] = leftProj[5]*p11*p22+2.0*leftProj[5]*p12sq; 
  wv[8] = leftProj[5]*p11*p23+2.0*leftProj[5]*p12*p13; 
  wv[9] = leftProj[5]*p11*p33+2.0*leftProj[5]*p13sq;

  mulByPhiPrime(p0, v[1], v[2], v[3], wv, &waves[30]);

  // Wave 5: (ev 7, 8, 9, 10 are repeated)
  s[4] = u1;
  wv[0] = leftProj[6]; 
  wv[1] = 0.0; 
  wv[2] = 0.0; 
  wv[3] = 0.0; 
  wv[4] = 0.0; 
  wv[5] = 0.0; 
  wv[6] = 0.0; 
  wv[7] = leftProj[7]; 
  wv[8] = leftProj[8]; 
  wv[9] = leftProj[9];

  mulByPhiPrime(p0, v[1], v[2], v[3], wv, &waves[40]);
  
  return fabs(u1)+sqrt(3*p11/p0);
}

GKYL_CU_DH
static void
qfluct_roe(const struct gkyl_wv_eqn *eqn,
  const double *ql, const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[10], *w2 = &waves[20], *w3 = &waves[30], *w4 = &waves[40];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]), s2m = fmin(0.0, s[2]), s3m = fmin(0.0, s[3]), s4m = fmin(0.0, s[4]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]), s2p = fmax(0.0, s[2]), s3p = fmax(0.0, s[3]), s4p = fmax(0.0, s[4]);

  for (int i=0; i<10; ++i) {
    amdq[i] = s0m*w0[i] + s1m*w1[i] + s2m*w2[i] + s3m*w3[i] + s4m*w4[i];
    apdq[i] = s0p*w0[i] + s1p*w1[i] + s2p*w2[i] + s3p*w3[i] + s4p*w4[i];
  }
}

GKYL_CU_DH
static void
wave_embed_absorb(const double *q, double *qphi, double *delta, void *ctx)
{
  qphi[0] = 0.0;
  qphi[1] = 0.0;
  qphi[2] = 0.0;
  qphi[3] = 0.0;
  qphi[4] = 0.0;
  qphi[5] = 0.0;
  qphi[6] = 0.0;
  qphi[7] = 0.0;
  qphi[8] = 0.0;
  qphi[9] = 0.0;
}

GKYL_CU_DH
static void
wave_embed_reflect(const double *q, double *qphi, double *delta, void *ctx)
{
  qphi[0] = q[0];
  qphi[1] = -q[1];
  qphi[2] = q[2];
  qphi[3] = q[3];
  qphi[4] = q[4];
  qphi[5] = -q[5];
  qphi[6] = -q[6];
  qphi[7] = q[7];
  qphi[8] = q[8];
  qphi[9] = q[9];
}

GKYL_CU_DH
static double
wave_embedded(const struct gkyl_wv_eqn *eqn,
  const double *delta, const double *ql, const double *qr,
  const double phil, const double phir, double *waves, double *s)
{
  double deltaphi[10] = {0.0};
  double qphi[10] = {0.0};

  double amax = 0.0;

  double fl[10], fr[10] = {0.0};
  if ((phil < 0.0) && (phir > 0.0)) {
    double sr = gkyl_ten_moment_max_abs_speed(qr);
    amax = sr;
   
    eqn->embed_geo->embed_func(qr, qphi, deltaphi, eqn->embed_geo->ctx);

    gkyl_ten_moment_flux(qphi, fl);
    gkyl_ten_moment_flux(qr, fr);

    double *w0 = &waves[0], *w1 = &waves[10];
    for (int i=0; i<10; ++i) {
      w1[i] = 0.5*((qr[i]-qphi[i]) + (fr[i]-fl[i])/amax);
    }
  }
  
  if ((phir < 0.0) && (phil > 0.0)) {
    double sl = gkyl_ten_moment_max_abs_speed(ql);
    amax = sl;
    
    eqn->embed_geo->embed_func(ql, qphi, deltaphi, eqn->embed_geo->ctx);

    gkyl_ten_moment_flux(ql, fl);
    gkyl_ten_moment_flux(qphi, fr);
    
    double *w0 = &waves[0], *w1 = &waves[10];
    for (int i=0; i<10; ++i) {
      w0[i] = 0.5*((qphi[i]-ql[i]) - (fr[i]-fl[i])/amax);
    }
  }

  s[0] = -amax;
  s[1] = amax;
  
  return s[1];
}

// Waves and speeds using Lax fluxes
GKYL_CU_DH
static double
wave_lax(const struct gkyl_wv_eqn *eqn,
  const double *delta, const double *ql, const double *qr, 
  double *waves, double *s)
{
  double sl = gkyl_ten_moment_max_abs_speed(ql);
  double sr = gkyl_ten_moment_max_abs_speed(qr);
  double amax = fmax(sl, sr);

  double fl[10], fr[10];
  gkyl_ten_moment_flux(ql, fl);
  gkyl_ten_moment_flux(qr, fr);

  double *w0 = &waves[0], *w1 = &waves[10];
  for (int i=0; i<10; ++i) {
    w0[i] = 0.5*((qr[i]-ql[i]) - (fr[i]-fl[i])/amax);
    w1[i] = 0.5*((qr[i]-ql[i]) + (fr[i]-fl[i])/amax);
  }

  s[0] = -amax;
  s[1] = amax;
  
  return s[1];
}

GKYL_CU_DH
static void
qfluct_lax(const struct gkyl_wv_eqn *eqn,
  const double *ql, const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[10];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i=0; i<10; ++i) {
    amdq[i] = s0m*w0[i] + s1m*w1[i];
    apdq[i] = s0p*w0[i] + s1p*w1[i];
  }
}

GKYL_CU_DH
static double
wave(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, 
  const double phil, const double phir, double *waves, double *s)
{
  // clear waves and wave speeds
  int mwaves = (type == GKYL_WV_HIGH_ORDER_FLUX) ? eqn->num_waves : 2;
  int meqn = eqn->num_equations;
  for (int i=0; i<mwaves; ++i) {
    double *w = &waves[i*meqn];
    for (int j=0; j<meqn; ++j)
      w[j] = 0.0;
    s[i] = 0.0;
  }
  if ((phil < 0.0) || (phir < 0.0))
    return wave_embedded(eqn, delta, ql, qr, phil, phir, waves, s);
  else {
    if (type == GKYL_WV_HIGH_ORDER_FLUX)
      return wave_roe(eqn, delta, ql, qr, waves, s);
    else
      return wave_lax(eqn, delta, ql, qr, waves, s);
  }

  return 0.0; // can't happen
}

GKYL_CU_DH
static void
qfluct(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double phil, const double phir, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX && (phil > 0.0) && (phir > 0.0))
    return qfluct_roe(eqn, ql, qr, waves, s, amdq, apdq);
  else
    return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
}

GKYL_CU_DH
static bool
check_inv(const struct gkyl_wv_eqn *eqn, const double *q)
{
  if (q[0] < 0.0) {
    return false;
  }

  double P[3];
  gkyl_ten_moment_diag_pressure(q, P);
  if (P[0] < 0.0 || P[1] < 0.0 || P[2] < 0.0) {
    return false;
  }
  
  return true;
}

GKYL_CU_DH
static double
max_speed(const struct gkyl_wv_eqn *eqn, const double *q)
{
  return gkyl_ten_moment_max_abs_speed(q);
}
