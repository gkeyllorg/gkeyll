#include <stdlib.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_mp_scheme.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <gkyl_wv_vacuum_einstein.h>
#include <gkyl_wv_vacuum_einstein_conformal.h>
#include <gkyl_wv_gr_euler.h>
#include <gkyl_wv_gr_euler_priv.h>
#include <gkyl_wave_geom.h>

#include <float.h>
#include <math.h>

// type signature for function to do recovery
typedef void (*recovery_fn_t)(int meqn,
  const double *f3m, const double *f2m, const double *fm,
  const double *fp, const double *f2p, const double *f3p,
  double *outl, double *outr);

struct gkyl_mp_scheme {
  struct gkyl_rect_grid grid; // grid object
  int ndim; // number of dimensions
  int num_up_dirs; // number of update directions
  int update_dirs[GKYL_MAX_DIM]; // directions to update
  enum gkyl_mp_recon mp_recon; // base reconstruction to use
  bool skip_mp_limiter; // should we skip MP limiter?  

  double cfl; // CFL number  
  
  const struct gkyl_wv_eqn *equation; // equation object
  struct gkyl_wave_geom *geom; // geometry object

  recovery_fn_t recovery_fn; // function to do recovery

  const struct gkyl_array *lo_mask; // optional per-cell mask (NULL = off): a face with a marked neighbour is first order
};

// Each of the recovery methods below take 6 cells, three to the left
// and three to the right, and returns the values at the left/right of
// the interface. Note that depending on the scheme, some of the input
// values may be ignored.

static inline void
c2_recovery(int meqn,
  const double *f3m, const double *f2m, const double *fm,
  const double *fp, const double *f2p, const double *f3p,
  double *outl, double *outr)
{
  // c2 is symmetric 2nd order scheme, so outl and outr are same
  for (int m=0; m<meqn; ++m)
    outr[m] = outl[m] = fm[m]/2.0 + fp[m]/2.0;
}

static inline void
c4_recovery(int meqn,
  const double *f3m, const double *f2m, const double *fm,
  const double *fp, const double *f2p, const double *f3p,
  double *outl, double *outr)
{
  // c4 is symmetric 4th order scheme, so outl and outr are same
  for (int m=0; m<meqn; ++m)
    outr[m] = outl[m] = -f2m[m]/12.0 + 7.0*fm[m]/12.0 + 7.0*fp[m]/12.0 - f2p[m]/12.0;
}

static inline void
c6_recovery(int meqn,
  const double *f3m, const double *f2m, const double *fm,
  const double *fp, const double *f2p, const double *f3p,
  double *outl, double *outr)
{
  // c6 is symmetric 6th order scheme, so outl and outr are same
  for (int m=0; m<meqn; ++m)  
    outr[m] = outl[m] =
      37.0*fp[m]/60.0+37.0*fm[m]/60.0+f3p[m]/60.0+f3m[m]/60.0-2.0*f2p[m]/15.0-2.0*f2m[m]/15.0;
}

static inline void
u1_recovery(int meqn,
  const double *f3m, const double *f2m, const double *fm,
  const double *fp, const double *f2p, const double *f3p,
  double *outl, double *outr)
{
  // u1 is upwind-biased 1st order scheme
  for (int m=0; m<meqn; ++m) {
    outl[m] = fm[m];
    outr[m] = fp[m];
  }
}

static inline void
u3_recovery(int meqn,
  const double *f3m, const double *f2m, const double *fm,
  const double *fp, const double *f2p, const double *f3p,
  double *outl, double *outr)
{
  // u3 is upwind-biased 3rd order scheme
  for (int m=0; m<meqn; ++m) {
    outl[m] = -1.0/6.0*f2m[m] + 5.0/6.0*fm[m] + 1.0/3.0*fp[m];
    outr[m] = 1.0/3.0*fm[m] + 5.0/6.0*fp[m] - 1.0/6.0*f2p[m];
  }
}

static inline void
u5_recovery(int meqn,
  const double *f3m, const double *f2m, const double *fm,
  const double *fp, const double *f2p, const double *f3p,
  double *outl, double *outr)
{
  // u5 is upwind-biased 5th order scheme
  for (int m=0; m<meqn; ++m) {
    outl[m] = 1.0/30.0*f3m[m] - 13.0/60.0*f2m[m] + 47.0/60.0*fm[m] + 9.0/20.0*fp[m] - 1.0/20.0*f2p[m];
    outr[m] = -1.0/20.0*f2m[m] + 9.0/20.0*fm[m] + 47.0/60.0*fp[m] - 13.0/60.0*f2p[m] + 1.0/30.0*f3p[m];
  }
}

static inline bool
vacuum_einstein_excision_in_stencil(const struct gkyl_wv_eqn *eqn,
  const double *qavg[6])
{
  double excision_threshold = 0.0;
  if (eqn->type == GKYL_EQN_VACUUM_EINSTEIN) {
    excision_threshold = gkyl_wv_vacuum_einstein_excision_threshold(eqn);
  }
  else if (eqn->type == GKYL_EQN_VACUUM_EINSTEIN_CONFORMAL) {
    excision_threshold = gkyl_wv_vacuum_einstein_conformal_excision_threshold(eqn);
  }
  else {
    return false;
  }

  for (int i = 0; i < 6; i++) {
    if (qavg[i][9] < excision_threshold) {
      return true;
    }
  }

  return false;
}

static inline double
minmod_2(double x, double y)
{
  if (x>0 && y>0)
    return fmin(x,y);
  if (x<0 && y<0)
    return fmax(x,y);
  return 0.0;
}

// Physical-constraint-preserving (PCP / Wu-Tang) limiter for GR-Euler. 
// gr_euler_pcp_theta(anchor, w) returns the largest theta in [0,1] such that anchor + theta*w stays in the admissible set
// {D>0, E=D+tau>0, E^2-D^2-|S|^2>0}. Applied to reconstruction only (one shared flux).
// For a static equilibrium the deviation w -> 0 so theta=1 and the WB is untouched; during collapse, anchoring at the (admissible) cell average keeps the
// face states admissible => admissibility-preserving updates (Wu & Tang 2015/2017).
static inline double
gr_euler_smallest_pos_root(double A, double B, double Cz)
{
  double scale = fmax(fabs(A), fmax(fabs(B), fabs(Cz)));
  if (!isfinite(scale)) return 0.0;
  if (scale == 0.0) return 2.0;
  A /= scale; B /= scale; Cz /= scale;
  if (fabs(A) < 1.0e-300) {
    if (B < 0.0) return -Cz / B; // single root, positive since Cz >= 0
    return 2.0; // increasing/flat from Cz >= 0: never crosses for t > 0
  }
  double disc = (B * B) - (4.0 * A * Cz);
  if (disc < 0.0) return 2.0; // no real root: same sign as the value at t=0 (>= 0)
  double sq = sqrt(disc);
  double q = -0.5 * (B + copysign(sq, B));
  double r1 = q / A;
  double r2 = q != 0.0 ? Cz / q : 0.0;
  double lo = fmin(r1, r2), hi = fmax(r1, r2);
  if (lo > 0.0) return lo;
  if (hi > 0.0) return hi;
  return 2.0;
}

// Momentum inner product gamma^{ij} u_i v_j for the covariant Valencia momentum density S_i.
static inline double
gr_euler_mom_dot(const double spatial_metric[9], const double u[3], const double vv[3])
{
  double det = (spatial_metric[0] * ((spatial_metric[4] * spatial_metric[8]) - (spatial_metric[5] * spatial_metric[7]))) -
    (spatial_metric[1] * ((spatial_metric[3] * spatial_metric[8]) - (spatial_metric[5] * spatial_metric[6]))) +
    (spatial_metric[2] * ((spatial_metric[3] * spatial_metric[7]) - (spatial_metric[4] * spatial_metric[6])));
  double inv_det = 1.0 / det;
  double inv_spatial_metric[9];
  inv_spatial_metric[0] = ((spatial_metric[4] * spatial_metric[8]) - (spatial_metric[5] * spatial_metric[7])) * inv_det;
  inv_spatial_metric[1] = ((spatial_metric[2] * spatial_metric[7]) - (spatial_metric[1] * spatial_metric[8])) * inv_det;
  inv_spatial_metric[2] = ((spatial_metric[1] * spatial_metric[5]) - (spatial_metric[2] * spatial_metric[4])) * inv_det;
  inv_spatial_metric[3] = ((spatial_metric[5] * spatial_metric[6]) - (spatial_metric[3] * spatial_metric[8])) * inv_det;
  inv_spatial_metric[4] = ((spatial_metric[0] * spatial_metric[8]) - (spatial_metric[2] * spatial_metric[6])) * inv_det;
  inv_spatial_metric[5] = ((spatial_metric[2] * spatial_metric[3]) - (spatial_metric[0] * spatial_metric[5])) * inv_det;
  inv_spatial_metric[6] = ((spatial_metric[3] * spatial_metric[7]) - (spatial_metric[4] * spatial_metric[6])) * inv_det;
  inv_spatial_metric[7] = ((spatial_metric[1] * spatial_metric[6]) - (spatial_metric[0] * spatial_metric[7])) * inv_det;
  inv_spatial_metric[8] = ((spatial_metric[0] * spatial_metric[4]) - (spatial_metric[1] * spatial_metric[3])) * inv_det;

  double s = 0.0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      s += inv_spatial_metric[(3 * i) + j] * u[i] * vv[j];
  return s;
}

// SPD test for a 3x3 spatial metric (row-major slots 9..17) via leading principal minors with a
// scale-relative floor (minors scale as scale^k). Assumes spatial_metric symmetric (caller symmetrizes first).
static inline bool
gr_euler_metric_is_spd(const double spatial_metric[9])
{
  double scale = fmax(fmax(fabs(spatial_metric[0]), fabs(spatial_metric[4])), fmax(fabs(spatial_metric[8]), 1.0e-300));
  double eps = 1.0e-12;
  // Test on the symmetrized metric (Sylvester's criterion needs symmetry); callers may pass a metric with small numerical asymmetry (e.g. the cell-average anchor).
  double spatial_metric00 = spatial_metric[0], spatial_metric11 = spatial_metric[4], spatial_metric22 = spatial_metric[8];
  double spatial_metric01 = 0.5 * (spatial_metric[1] + spatial_metric[3]), spatial_metric02 = 0.5 * (spatial_metric[2] + spatial_metric[6]), spatial_metric12 = 0.5 * (spatial_metric[5] + spatial_metric[7]);
  double m1 = spatial_metric00;
  double m2 = (spatial_metric00 * spatial_metric11) - (spatial_metric01 * spatial_metric01);
  double m3 = (spatial_metric00 * ((spatial_metric11 * spatial_metric22) - (spatial_metric12 * spatial_metric12))) - (spatial_metric01 * ((spatial_metric01 * spatial_metric22) - (spatial_metric12 * spatial_metric02))) +
    (spatial_metric02 * ((spatial_metric01 * spatial_metric12) - (spatial_metric11 * spatial_metric02)));
  return (m1 > eps * scale) && (m2 > eps * scale * scale) && (m3 > eps * scale * scale * scale);
}

// An invalid reconstructed metric falls back to its positive cell anchor.
// Maximizing a blend up to the SPD boundary creates an almost singular metric
// and unbounded coordinate wave speeds. Admissible faces remain unchanged;
// both reconstructed states still feed one shared conservative fluid flux.
static inline void
gr_euler_spd_face_limiter(double spatial_metric_face[9], const double spatial_metric_cell[9])
{
  double s;
  s = 0.5 * (spatial_metric_face[1] + spatial_metric_face[3]); spatial_metric_face[1] = s; spatial_metric_face[3] = s; // symmetrize (0,1)/(1,0)
  s = 0.5 * (spatial_metric_face[2] + spatial_metric_face[6]); spatial_metric_face[2] = s; spatial_metric_face[6] = s; // (0,2)/(2,0)
  s = 0.5 * (spatial_metric_face[5] + spatial_metric_face[7]); spatial_metric_face[5] = s; spatial_metric_face[7] = s; // (1,2)/(2,1)
  if (gr_euler_metric_is_spd(spatial_metric_face)) return;  // already SPD: no-op (linearity-preserving)

  // Do not reintroduce cell-metric asymmetry in the fallback.
  double gc[9];
  for (int m = 0; m < 9; m++) gc[m] = spatial_metric_cell[m];
  gc[1] = gc[3] = 0.5 * (spatial_metric_cell[1] + spatial_metric_cell[3]);
  gc[2] = gc[6] = 0.5 * (spatial_metric_cell[2] + spatial_metric_cell[6]);
  gc[5] = gc[7] = 0.5 * (spatial_metric_cell[5] + spatial_metric_cell[7]);
  if (!gr_euler_metric_is_spd(gc)) return; // anchor not SPD: caught upstream by the cell-average guard

  for (int m = 0; m < 9; m++) spatial_metric_face[m] = gc[m];
}

static inline double
gr_euler_surface_u1_factor(void)
{
  static double factor = -1.0;
  if (factor < 0.0) {
    factor = 1.0e5;
    const char *env = getenv("GKYL_GR_EULER_SURFACE_U1_FACTOR");
    if (env != NULL) {
      factor = atof(env);
      if (factor < 0.0) factor = 0.0;
    }
  }
  return factor;
}

static inline bool
gr_euler_low_density_cell(const struct gkyl_wv_eqn *eqn, int meqn, const double *q)
{
  if (q[27] < pow(10.0, -8.0)) return true;

  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
  double rho_atm = gr_euler->rho_atm;
  double factor = gr_euler_surface_u1_factor();
  if (!(rho_atm > 0.0) || !(factor > 0.0)) return false;

  double spatial_metric[9] = {
    q[9],  q[10], q[11],
    q[12], q[13], q[14],
    q[15], q[16], q[17],
  };
  double det = (spatial_metric[0] * ((spatial_metric[4] * spatial_metric[8]) - (spatial_metric[5] * spatial_metric[7]))) -
    (spatial_metric[1] * ((spatial_metric[3] * spatial_metric[8]) - (spatial_metric[5] * spatial_metric[6]))) +
    (spatial_metric[2] * ((spatial_metric[3] * spatial_metric[7]) - (spatial_metric[4] * spatial_metric[6])));
  if (!(isfinite(det) && det > 0.0)) return true;

  double cutoff = factor * rho_atm;
  double sqrt_det = sqrt(det);
  double rho_live = q[0] / sqrt_det;
  if (isfinite(rho_live) && rho_live <= cutoff) return true;
  if (meqn > 71 && isfinite(q[71]) && q[71] >= 0.0) {
    double rho_ref = q[71] / sqrt_det;
    if (rho_ref <= cutoff) return true;
  }
  return false;
}

static inline bool
gr_euler_surface_u1_stencil(const struct gkyl_wv_eqn *eqn, int meqn, const double *qavg[6])
{
  if (!gkyl_wv_gr_euler_atmosphere_floor() || getenv("GKYL_GR_EULER_SURFACE_U1_OFF") != NULL) return false;
  if (!gkyl_gr_euler_wb_disabled(eqn)) return false;

  for (int k = 0; k < 6; ++k)
    if (gr_euler_low_density_cell(eqn, meqn, qavg[k])) return true;
  return false;
}

bool
gkyl_mp_scheme_gr_metric_is_spd(const double spatial_metric[9])
{
  return gr_euler_metric_is_spd(spatial_metric);
}

void
gkyl_mp_scheme_gr_spd_face_limiter(double spatial_metric_face[9], const double spatial_metric_cell[9])
{
  gr_euler_spd_face_limiter(spatial_metric_face, spatial_metric_cell);
}

double
gkyl_mp_scheme_gr_euler_pcp_theta(const double qe[5], const double w[5], const double spatial_metric[9])
{
  const double frac = 1.0e-4; // keep the face D and cone at >= frac of the equilibrium's own margin.
  if (!gr_euler_metric_is_spd(spatial_metric)) return 0.0;
  for (int m=0; m<5; ++m)
    if (!isfinite(qe[m]) || !isfinite(w[m])) return 0.0;

  double D0 = qe[0];
  double E0 = qe[0] + qe[4];                                      // E_eq = D_eq + tau_eq
  double dD = w[0];
  double dE = w[0] + w[4];
  double qeS[3] = { qe[1], qe[2], qe[3] };
  double wS[3]  = { w[1], w[2], w[3] };
  double Seq2 = gr_euler_mom_dot(spatial_metric, qeS, qeS);
  double C = (E0 * E0) - (D0 * D0) - Seq2;                        // g(0) = equilibrium pressure-positive cone (> 0)
  if (!(isfinite(D0) && isfinite(E0) && isfinite(Seq2) && isfinite(C)) ||
      !(D0 > 0.0) || !(E0 > 0.0) || !(C > 0.0)) {
    return 0.0;
  }
  double Dmin = frac * D0;
  double gmin = frac * C;

  // Unlimited (theta = 1) face state.
  double D1 = D0 + dD;
  double E1 = E0 + dE;
  double s1[3] = { qe[1] + w[1], qe[2] + w[2], qe[3] + w[3] };
  double S1sq = gr_euler_mom_dot(spatial_metric, s1, s1);
  double g1 = (E1 * E1) - (D1 * D1) - S1sq;
  if (!(isfinite(D1) && isfinite(E1) && isfinite(g1))) return 0.0;
  if (isfinite(E1) && isfinite(g1) && E1 > 0.0 && D1 >= Dmin && g1 >= gmin) {
    return 1.0; // already admissible with margin: no limiting.
  }

  double theta = 1.0;

  // D(theta) = qe[0] + theta*w[0] >= Dmin (only binds when the deviation lowers D).
  if (dD < 0.0) {
    double thD = (Dmin - D0) / dD;
    if (thD < theta) theta = thD;
  }

  // g(theta) = A theta^2 + B theta + C >= gmin. value at theta=0 is C >= gmin (since gmin = frac*C).
  if (theta > 0.0) {
    double sth[3] = { qe[1] + theta * w[1], qe[2] + theta * w[2], qe[3] + theta * w[3] };
    double Dth = D0 + theta * dD;
    double Eth = E0 + theta * dE;
    double Sth2 = gr_euler_mom_dot(spatial_metric, sth, sth);
    if (Eth <= 0.0 || (Eth * Eth) - (Dth * Dth) - Sth2 < gmin) {
      double A = (dE * dE) - (dD * dD) - gr_euler_mom_dot(spatial_metric, wS, wS);
      double B = 2.0 * ((E0 * dE) - (D0 * dD) - gr_euler_mom_dot(spatial_metric, qeS, wS));
      if (!(isfinite(A) && isfinite(B))) return 0.0;
      double thg = gr_euler_smallest_pos_root(A, B, C - gmin);
      if (thg < theta) theta = thg;
    }
  }

  if (theta < 0.0) theta = 0.0;
  if (theta > 1.0) theta = 1.0;
  return theta;
}

// PCP anchor: the theta-sweep limits toward the cell average qa, but admissibility is measured with the
// reconstructed face metric spatial_metric, so the anchor itself must be admissible under it (cone
// C = E0^2 - D0^2 - gamma^{ij} S_i S_j > 0) - otherwise the theta=0 endpoint is already invalid and the PCP proof fails.
// Build a safe anchor q0 that preserves D0 and E0=D0+tau0 and reduces only the momentum by lambda so
// E0^2 - D0^2 - gamma^{ij} S0_i S0_j >= cone_floor > 0. No-op where the anchor is already admissible (smooth regions).
// Conservation-safe: only the states fed to the single shared face flux change.
void
gkyl_mp_scheme_gr_euler_pcp_safe_anchor(const double qa[5], const double spatial_metric[9], double q0[5])
{
  for (int m = 0; m < 5; m++) q0[m] = qa[m];
  double D0 = q0[0];
  double E0 = q0[0] + q0[4];
  if (!(D0 > 0.0) || !(E0 > 0.0)) return; // degenerate energy: nothing to preserve (excised/vacuum handled elsewhere)
  double S0[3] = { q0[1], q0[2], q0[3] };
  double S02 = gr_euler_mom_dot(spatial_metric, S0, S0);
  double pressure_margin = (E0 * E0) - (D0 * D0);
  double cone_floor = 1.0e-12 * fmax(E0 * E0, D0 * D0); // tiny positive cone margin (relative, scale-free)
  if (!(isfinite(S02) && isfinite(pressure_margin) && isfinite(cone_floor))) return;
  if (pressure_margin <= cone_floor) {
    // Face-only fallback for an almost pressureless/roundoff anchor: make the endpoint strictly
    // admissible without changing the cell average.
    q0[1] = 0.0;
    q0[2] = 0.0;
    q0[3] = 0.0;
    double E_safe = sqrt((D0 * D0) + cone_floor);
    q0[4] = E_safe - D0;
    return;
  }
  if (S02 > pressure_margin - cone_floor) {
    double lam = sqrt((pressure_margin - cone_floor) / S02);
    if (lam < 1.0) { q0[1] *= lam; q0[2] *= lam; q0[3] *= lam; } // reduce momentum only; D0, E0 preserved
  }
}

static inline double
minmod_4(double x, double y, double z, double w)
{
  if (x>0 && y>0 && z>0 && w>0)
    return fmin(fmin(x,y),fmin(z,w));
  if (x<0 && y<0 && z<0 && w<0)
    return fmax(fmax(x,y),fmax(z,w));
  return 0.0;
}

static inline double
median(double x, double y, double z)
{
  return x + minmod_2(y-x,z-x);
}

static inline double
min_3(double x, double y, double z)
{
  return fmin(x,fmin(y,z));
}

static inline double
max_3(double x, double y, double z)
{
  return fmax(x,fmax(y,z));
}

// MP limiter: See Eqns 3.44 - 3.57 of Peterson and Hammett SIAM
// J. Sci. Comput, vol 35, No 3 pp B576, 2013
static inline double
mp_limiter(double qe, double q2m, double q1m, double q0, double q1p, double q2p)
{
  double alpha = 4.0;
 // Suresh and Huynh recommend 1e-10, but that seems turns off the
 // limiter when the jumps are very tiny, leading to small-scale noise
 // in certain situations. Not sure what this should be. Greg Hammett
 // recommends eps = 0.0
  double eps = 0.0;
  // Eq 3.44
  double qmp = q0 + minmod_2(q1p-q0, alpha*(q0-q1m));

  //return median(qe, q0, qmp);

  // Eq 3.45
  if ((qe-q0)*(qe-qmp)<eps) return qe; // no need to apply limiter

  // Eq 3.46-3.48
  double d1m = q2m + q0 - 2*q1m;
  double d0 = q1m + q1p - 2*q0;
  double d1p = q0 + q2p - 2*q1p;
  // Eq 3.49 and 3.50
  double dp_m4 = minmod_4(4*d0-d1p, 4*d1p-d0, d0, d1p);
  double dm_m4 = minmod_4(4*d1m-d0, 4*d0-d1m, d1m, d0);

  // Eq 3.51-3.54
  double qul = q0 + alpha*(q0-q1m);
  double qavg = 0.5*(q0+q1p);
  double qmd = qavg - 0.5*dp_m4;
  double qlc = q0 + 0.5*(q0-q1m) + 4.0/3.0*dm_m4;

  // Eq 3.55 and 3.56
  double qmin = fmax(min_3(q0,q1p,qmd), min_3(q0,qul,qlc));
  double qmax = fmin(max_3(q0,q1p,qmd), max_3(q0,qul,qlc));

  // 3.57
  return median(qe, qmin, qmax);
}

static inline long
get_offset(int dir, int loc, const struct gkyl_range *range)
{
  int idx[GKYL_MAX_CDIM] = { 0, 0, 0 };
  idx[dir] = loc;
  return gkyl_range_offset(range, idx);
}

gkyl_mp_scheme*
gkyl_mp_scheme_new(const struct gkyl_mp_scheme_inp *mpinp)
{
  struct gkyl_mp_scheme *mp = gkyl_malloc(sizeof *mp);

  mp->grid = *(mpinp->grid);
  mp->ndim = mp->grid.ndim;
  
  mp->num_up_dirs = mpinp->num_up_dirs;
  for (int i=0; i<mpinp->num_up_dirs; ++i)
    mp->update_dirs[i] = mpinp->update_dirs[i];

  mp->mp_recon = mpinp->mp_recon;
  mp->skip_mp_limiter = mpinp->skip_mp_limiter;
  mp->lo_mask = NULL;
  
  mp->cfl = mpinp->cfl;
  
  mp->equation = gkyl_wv_eqn_acquire(mpinp->equation);
  mp->geom = gkyl_wave_geom_acquire(mpinp->geom);

  switch (mpinp->mp_recon) {
    case GKYL_MP_C2:
      mp->recovery_fn = c2_recovery;
      break;
    case GKYL_MP_C4:
      mp->recovery_fn = c4_recovery;
      break;
    case GKYL_MP_C6:
      mp->recovery_fn = c6_recovery;
      break;
    case GKYL_MP_U1:
      mp->recovery_fn = u1_recovery;
      break;
    case GKYL_MP_U3:
      mp->recovery_fn = u3_recovery;
      break;
    case GKYL_MP_U5:
      mp->recovery_fn = u5_recovery;
      break;
  }

  return mp;
}

void
gkyl_mp_scheme_advance(gkyl_mp_scheme *mp,
  const struct gkyl_range *update_range, const struct gkyl_array *qin,
  struct gkyl_array *qrec_l, struct gkyl_array *qrec_r,
  struct gkyl_array *amdq, struct gkyl_array *apdq,
  struct gkyl_array *cflrate, struct gkyl_array *phi, struct gkyl_array *rhs)
{
  int ndim = update_range->ndim;
  int meqn = mp->equation->num_equations;
  int mwaves = mp->equation->num_waves;

  double apdq_local[meqn], amdq_local[meqn];
  double qlocal_r[meqn], qlocal_l[meqn];
  double delta[meqn], waves[meqn*mwaves], speeds[mwaves];
  // labels for three cells to left, three cells to right of edge:
  enum { I3M, I2M, IM, IP, I2P, I3P };

  gkyl_array_clear_range(rhs, 0.0, update_range);
  
  // outer loop is over direction: the RHS is updated direction
  // by direction
  for (int d=0; d<mp->num_up_dirs; ++d) {
    int dir = mp->update_dirs[d];

    double dx = mp->grid.dx[dir];

    // compute index offsets of cells on left/right of edge
    long offsets[6];
    offsets[IP]  = get_offset(dir, 0, update_range);
    offsets[I2P] = get_offset(dir, 1, update_range);
    offsets[I3P] = get_offset(dir, 2, update_range);
    offsets[IM]  = get_offset(dir, -1, update_range);
    offsets[I2M] = get_offset(dir, -2, update_range);
    offsets[I3M] = get_offset(dir, -3, update_range);

    const double *qavg[6]; // pointers to cells attached to edge

    // create range that includes one extra layer on the upper size
    int upper[GKYL_MAX_CDIM] = { 0 };
    for (int d=0; d<update_range->ndim; ++d) upper[d] = update_range->upper[d];
    upper[dir] += 1;
    struct gkyl_range update_range_ext;
    gkyl_range_init(&update_range_ext, update_range->ndim, update_range->lower, upper);

    struct gkyl_range_iter iter;
    // loop over all cells and recover left/right edge values in
    // direction 'dir'. Note this is effectively a loop over the
    // edges, and includes upper most edge also
    gkyl_range_iter_init(&iter, &update_range_ext);
    while (gkyl_range_iter_next(&iter)) {
      // Note: edge is between cells IM and IP
      
      long loc = gkyl_range_idx(update_range, iter.idx);

      // attach pointers to cells for recovery
      for (int i=0; i<6; ++i)
        qavg[i] = gkyl_array_cfetch(qin, loc+offsets[i]);

      // qr_l is left of edge (right edge of left cell), qr_r right of
      // edge (left edge of right cell)
      double *qr_l = gkyl_array_fetch(qrec_r, loc+offsets[IM]);
      double *qr_r = gkyl_array_fetch(qrec_l, loc+offsets[IP]);

      const double *phil = gkyl_array_cfetch(phi, loc+offsets[IM]);
      const double *phir = gkyl_array_cfetch(phi, loc+offsets[IP]);

      bool use_first_order_excision_buffer = vacuum_einstein_excision_in_stencil(
        mp->equation, qavg);
      // First-order flux correction (Lemaster & Stone 2009): a face with a marked neighbour is reconstructed
      // piecewise constant. One face, shared by both cells, so the flux cancellation is preserved.
      bool lo_face = mp->lo_mask
        && (((const double *) gkyl_array_cfetch(mp->lo_mask, loc+offsets[IM]))[0] > 0.0
          || ((const double *) gkyl_array_cfetch(mp->lo_mask, loc+offsets[IP]))[0] > 0.0);
      bool gr_euler_76 = mp->equation->type == GKYL_EQN_GR_EULER && meqn == 76;
      if (use_first_order_excision_buffer || (lo_face && !gr_euler_76)) {
        // This face is shared by the cells on both sides => replacing its reconstruction once preserves the flux cancellation
        u1_recovery(meqn, qavg[I3M], qavg[I2M], qavg[IM],
          qavg[IP], qavg[I2P], qavg[I3P], qr_l, qr_r);
      }
      // GR-Euler with the static-TOV slots (meqn==76). When well-balanced, reconstruct the smooth
      // deviation w = q - q_eq and add a single-valued edge equilibrium q_eq(edge) to both sides (at a
      // static equilibrium w = 0 -> qr_l = qr_r = q_eq(edge): no jump -> exact balance). When the WB is
      // disabled (collapse stage), q_eq = 0 -> plain reconstruction of q. The PCP admissibility limiter
      // is applied afer either path (independent of the WB), so it also protects the no-WB collapse.
      else if (mp->equation->type == GKYL_EQN_GR_EULER && meqn == 76) {
        // WB on (wb=true): subtract the frozen equilibrium q_eq per stencil cell and reconstruct the smooth
        // deviation w=q-q_eq, then add back the single-valued edge equilibrium qe (which telescopes =>
        // conservative). WB off (collapse stage): q_eq=0 => plain reconstruction of q. Reference slots (>=71)
        // are carried through (they are frozen and produce no flux).
        bool wb = !gkyl_gr_euler_wb_disabled(mp->equation);
        const double *qa[6] = { qavg[I3M], qavg[I2M], qavg[IM], qavg[IP], qavg[I2P], qavg[I3P] };
        double weq[6][meqn], qeq[6][meqn];
        for (int k = 0; k < 6; ++k) {
          if (wb) {
            double qeqk[71];
            gkyl_gr_euler_equilibrium(mp->equation, qa[k], qeqk);
            for (int m = 0; m < meqn; ++m) qeq[k][m] = (m < 71) ? qeqk[m] : qa[k][m];
          }
          else {
            for (int m = 0; m < meqn; ++m) qeq[k][m] = (m < 71) ? 0.0 : qa[k][m];
          }
          for (int m = 0; m < meqn; ++m) weq[k][m] = qa[k][m] - qeq[k][m];
        }
        double wl[meqn], wr[meqn];
        bool use_surface_u1 = lo_face || gr_euler_surface_u1_stencil(mp->equation, meqn, qavg);
        if (use_surface_u1) {
          // Stellar collapse uses an artificial atmosphere. Near the surface the
          // high-order MP polynomial can overshoot the Valencia cone even with
          // admissible face limiting, so use the standard finite-volume fallback:
          // piecewise-constant reconstruction only on low-density surface stencils.
          u1_recovery(meqn, weq[0], weq[1], weq[2], weq[3], weq[4], weq[5], wl, wr);
        }
        else {
          mp->recovery_fn(meqn, weq[0], weq[1], weq[2], weq[3], weq[4], weq[5], wl, wr);
        }
        if (!use_surface_u1 && !mp->skip_mp_limiter) {
          for (int m = 0; m < meqn; ++m) {
            wr[m] = mp_limiter(wr[m], weq[5][m], weq[4][m], weq[3][m], weq[2][m], weq[1][m]);
            wl[m] = mp_limiter(wl[m], weq[0][m], weq[1][m], weq[2][m], weq[3][m], weq[4][m]);
          }
        }
        // single-valued (2nd-order centered, monotone) equilibrium at the edge, from the two adjacent cells
        for (int m = 0; m < meqn; ++m) {
          double qe = 0.5 * (qeq[2][m] + qeq[3][m]); // IM, IP (zero when WB disabled)
          qr_l[m] = qe + wl[m];
          qr_r[m] = qe + wr[m];
        }

        // SPD face-metric limiter (metric analog of the fluid PCP), applied BEFORE the fluid PCP/C2P so
        // they see a valid SPD face metric. No-op where the reconstructed metric is already SPD (smooth
        // regions => bit-identical, linearity-preserving); only repairs a non-SPD reconstructed face metric
        // by falling back to the SPD cell-average anchor. Metric slots are zero-flux for the fluid.
        gr_euler_spd_face_limiter(&qr_l[9], &qa[2][9]); // left face state, anchor = IM cell average
        gr_euler_spd_face_limiter(&qr_r[9], &qa[3][9]); // right face state, anchor = IP cell average

        // PCP (Wu-Tang) admissibility limiter: scale each reconstructed face toward its own (admissible) cell average so the face states stay in {D>0, (D+tau)^2-D^2-|S|^2>0}, making the updates
        // admissibility-preserving (Wu & Tang 2015/2017). Applied whether or not the WB is active. At a static equilibrium the deviation ~0 so theta=1 => WB reconstruction untouched.
        static int pcp_on = -1;
        if (pcp_on < 0) pcp_on = (getenv("GKYL_PCP") != NULL) ? 1 : 0;
        if (pcp_on) {
          // Anchor at a face-metric-admissible state (safe anchor) so the theta=0 endpoint is valid (strict PCP).
          double q0l[5], q0r[5];
          if (wb) {
            // Limit deviations about an admissible face equilibrium, not about
            // different neighboring cell averages. Carry that SAME equilibrium
            // into the flux subtraction, so w=0 remains an exact fixed point.
            double edge_eq[5];
            for (int m=0; m<5; ++m) edge_eq[m]=0.5*(qeq[2][m]+qeq[3][m]);
            gkyl_mp_scheme_gr_euler_pcp_safe_anchor(edge_eq,&qr_l[9],q0l);
            gkyl_mp_scheme_gr_euler_pcp_safe_anchor(edge_eq,&qr_r[9],q0r);
            for (int m=0; m<5; ++m) {
              qr_l[m]=q0l[m]+wl[m]; qr_r[m]=q0r[m]+wr[m];
            }
            qr_l[71]=q0l[0]; qr_l[72]=q0l[4];
            qr_r[71]=q0r[0]; qr_r[72]=q0r[4];
            for (int i=0; i<3; ++i) {
              qr_l[73+i]=q0l[1+i]; qr_r[73+i]=q0r[1+i];
            }
          }
          else {
            gkyl_mp_scheme_gr_euler_pcp_safe_anchor(qa[2], &qr_l[9], q0l);
            gkyl_mp_scheme_gr_euler_pcp_safe_anchor(qa[3], &qr_r[9], q0r);
          }
          double devl[5], devr[5];
          for (int m = 0; m < 5; ++m) { devl[m] = qr_l[m] - q0l[m]; devr[m] = qr_r[m] - q0r[m]; }
          double thl = gkyl_mp_scheme_gr_euler_pcp_theta(q0l, devl, &qr_l[9]);
          double thr = gkyl_mp_scheme_gr_euler_pcp_theta(q0r, devr, &qr_r[9]);
          for (int m = 0; m < 5; ++m) {
            qr_l[m] = thl == 0.0 ? q0l[m] : q0l[m] + thl * devl[m];
            qr_r[m] = thr == 0.0 ? q0r[m] : q0r[m] + thr * devr[m];
          }
        }
      }
      else {
        // recover variables at cell edge
        mp->recovery_fn(meqn, qavg[I3M], qavg[I2M], qavg[IM],
          qavg[IP], qavg[I2P], qavg[I3P],
          qr_l, qr_r);

        if (!mp->skip_mp_limiter) {
          // apply MP limiter to left and right edge recovered values
          for (int m=0; m<meqn; ++m) {
            qr_r[m] = mp_limiter(qr_r[m],
              qavg[I3P][m], qavg[I2P][m], qavg[IP][m],
              qavg[IM][m], qavg[I2M][m]);

            qr_l[m] = mp_limiter(qr_l[m],
              qavg[I3M][m], qavg[I2M][m], qavg[IM][m],
              qavg[IP][m], qavg[I2P][m]);
          }
        }
      }

      const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(mp->geom, iter.idx);

      // rotate ql and qr to local frame
      mp->equation->rotate_to_local_func(mp->equation, cg->tau1[dir], cg->tau2[dir], cg->norm[dir], qr_l, qlocal_l);
      mp->equation->rotate_to_local_func(mp->equation, cg->tau1[dir], cg->tau2[dir], cg->norm[dir], qr_r, qlocal_r);

      for (int m=0; m<meqn; ++m) delta[m] = qlocal_r[m]-qlocal_l[m];
      
      // compute waves and fluctuations
      gkyl_wv_eqn_waves(mp->equation, GKYL_WV_HIGH_ORDER_FLUX, delta,
        qlocal_l, qlocal_r, phil[0], phir[0], waves, speeds);
      gkyl_wv_eqn_qfluct(mp->equation, GKYL_WV_HIGH_ORDER_FLUX,
        qlocal_l, qlocal_r, phil[0], phir[0], waves, speeds, amdq_local, apdq_local);

      double *amdq_p = gkyl_array_fetch(amdq, loc+offsets[IM]);
      double *apdq_p = gkyl_array_fetch(apdq, loc+offsets[IP]);

      // rotate fluctuations back to global frame
      mp->equation->rotate_to_global_func(mp->equation, cg->tau1[dir], cg->tau2[dir], cg->norm[dir], amdq_local, amdq_p);
      mp->equation->rotate_to_global_func(mp->equation, cg->tau1[dir], cg->tau2[dir], cg->norm[dir], apdq_local, apdq_p);
    }

    double deltaf_local[meqn], deltaf[meqn];
    // Update RHS with contribution from flux jumps. Note this loop is
    // over interior cells
    gkyl_range_iter_init(&iter, update_range);
    while (gkyl_range_iter_next(&iter)) {
      long loc = gkyl_range_idx(update_range, iter.idx);

      // THIS NEEDS TO CHANGE WHEN DOING GENERAL GEOMETRY!!!
      const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(mp->geom, iter.idx);

      const double *qr_l = gkyl_array_cfetch(qrec_l, loc);
      const double *qr_r = gkyl_array_cfetch(qrec_r, loc);
      
      // rotate ql and qr to local frame
      mp->equation->rotate_to_local_func(mp->equation, cg->tau1[dir], cg->tau2[dir], cg->norm[dir], qr_l, qlocal_l);
      mp->equation->rotate_to_local_func(mp->equation, cg->tau1[dir], cg->tau2[dir], cg->norm[dir], qr_r, qlocal_r);

      double amax = gkyl_wv_eqn_flux_jump(mp->equation, qlocal_l, qlocal_r, deltaf_local);

      // rotate deltaf back to global frame
      mp->equation->rotate_to_global_func(mp->equation, cg->tau1[dir], cg->tau2[dir], cg->norm[dir], deltaf_local, deltaf);

      const double *amdq_p = gkyl_array_cfetch(amdq, loc);
      const double *apdq_p = gkyl_array_cfetch(apdq, loc);

      double *rhs_p = gkyl_array_fetch(rhs, loc);
      for (int m=0; m<meqn; ++m)
        rhs_p[m] += -deltaf[m]/dx - (apdq_p[m]+amdq_p[m])/dx;

      double *cflrate_p = gkyl_array_fetch(cflrate, loc);
      cflrate_p[0] += amax/dx;
    }
  }
}

double
gkyl_mp_scheme_max_dt(const gkyl_mp_scheme *mp, const struct gkyl_range *update_range,
  const struct gkyl_array *qin)
{
  double max_dt = DBL_MAX;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);
  while (gkyl_range_iter_next(&iter)) {

    for (int d=0; d<mp->num_up_dirs; ++d) {
      int dir = mp->update_dirs[d];
      double dx = mp->grid.dx[dir];

      const double *q = gkyl_array_cfetch(qin, gkyl_range_idx(update_range, iter.idx));
      double maxs = gkyl_wv_eqn_max_speed(mp->equation, q);
      max_dt = fmin(max_dt, mp->cfl*dx/maxs);
    }
    
  }


  return max_dt;  
}

void
gkyl_mp_scheme_release(gkyl_mp_scheme* mp)
{
  gkyl_wv_eqn_release(mp->equation);
  gkyl_wave_geom_release(mp->geom);
  
  gkyl_free(mp);
}

// Optional per-cell mask for first-order flux correction; NULL (the default) switches it off. The mask has the layout
// of the advanced array (one component, > 0 = marked) and must outlive its use.
void
gkyl_mp_scheme_set_lo_mask(gkyl_mp_scheme *mp, const struct gkyl_array *mask)
{
  mp->lo_mask = mask;
}
