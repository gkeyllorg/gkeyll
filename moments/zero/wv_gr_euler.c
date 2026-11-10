#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_mp_scheme.h>
#include <gkyl_wv_gr_euler.h>
#include <gkyl_wv_gr_euler_priv.h>

// Density/pressure floors for the prims recovery, they were 1.0e-8 before, now I just set a variable so gkyl_wv_gr_euler_inew can override them from inp->rho_atm / inp->p_atm when provided.
// one set of floors per process; if two GR-Euler species ever need different floors, thread them through prim_vars per equation
static double gr_euler_rho_floor = 1.0e-8;;
static double gr_euler_p_floor = 1.0e-8;;
// Atmosphere floor in the primitive recovery (gkyl_wv_gr_euler_set_atmosphere_floor). The mixed app turns it on
// for its GR-Euler species; other schemes keep the per-variable floor below.
static bool gr_euler_atm_floor = false;

void
gkyl_wv_gr_euler_set_atmosphere_floor(bool on)
{
  gr_euler_atm_floor = on;
}

bool
gkyl_wv_gr_euler_atmosphere_floor(void)
{
  return gr_euler_atm_floor;
}
#define GR_EULER_BLEND_FACTOR 1.0e5  // Low-density/surface switch for HLL/Lax blending.

// Specific enthalpy h = 1 + eps + p/rho from (rho, p), gamma-law EOS.
double
gkyl_gr_euler_specific_enthalpy(double gas_gamma, double rho, double p)
{
  return 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0)));
}

static inline double
gr_euler_spatial_metric_det(const double spatial_metric[3][3])
{
  return (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
    (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
    (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));
}

static inline void
gr_euler_spatial_metric_inv(const double spatial_metric[3][3], double spatial_det, double inv_spatial_metric[3][3])
{
  double inv_det = 1.0 / spatial_det;

  inv_spatial_metric[0][0] = ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][1])) * inv_det;
  inv_spatial_metric[0][1] = ((spatial_metric[0][2] * spatial_metric[2][1]) - (spatial_metric[0][1] * spatial_metric[2][2])) * inv_det;
  inv_spatial_metric[0][2] = ((spatial_metric[0][1] * spatial_metric[1][2]) - (spatial_metric[0][2] * spatial_metric[1][1])) * inv_det;
  inv_spatial_metric[1][0] = ((spatial_metric[1][2] * spatial_metric[2][0]) - (spatial_metric[1][0] * spatial_metric[2][2])) * inv_det;
  inv_spatial_metric[1][1] = ((spatial_metric[0][0] * spatial_metric[2][2]) - (spatial_metric[0][2] * spatial_metric[2][0])) * inv_det;
  inv_spatial_metric[1][2] = ((spatial_metric[0][2] * spatial_metric[1][0]) - (spatial_metric[0][0] * spatial_metric[1][2])) * inv_det;
  inv_spatial_metric[2][0] = ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])) * inv_det;
  inv_spatial_metric[2][1] = ((spatial_metric[0][1] * spatial_metric[2][0]) - (spatial_metric[0][0] * spatial_metric[2][1])) * inv_det;
  inv_spatial_metric[2][2] = ((spatial_metric[0][0] * spatial_metric[1][1]) - (spatial_metric[0][1] * spatial_metric[1][0])) * inv_det;
}

// Building the static-TOV equilibrium conserved state q_eq at the same point/metric as q, for the optional well-balancing.
// The frozen t=0 equilibrium is carried in the spare slots: D_eq (71), Etot_eq (72), and the covariant momentum S_eq,i (73-75, which is nonzero in Kerr-Schild where v=beta/alpha). q_eq[5..70] = the live metric (copied from q), so the equilibrium tracks the current geometry.
static void
gr_euler_tov_equilibrium(const struct wv_gr_euler *gr_euler, const double q[76], double q_eq[71])
{
  for (int i = 0; i < 71; i++) {
    q_eq[i] = q[i];
  }
  q_eq[0] = q[71]; // D_eq
  q_eq[1] = q[73]; // S_eq,x
  q_eq[2] = q[74]; // S_eq,y
  q_eq[3] = q[75]; // S_eq,z
  q_eq[4] = q[72]; // Etot_eq (slot-[4] Valencia energy E - D)
}

void
gkyl_gr_euler_equilibrium(const struct gkyl_wv_eqn *eqn, const double q[76], double q_eq[71])
{
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
  gr_euler_tov_equilibrium(gr_euler, q, q_eq);
}

bool
gkyl_gr_euler_wb_disabled(const struct gkyl_wv_eqn *eqn)
{
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
  return gr_euler->disable_well_balanced;
}

// Runtime toggle of the frozen-equilibrium well-balancing (used by the mid-run dynamical-collapse
// trigger to switch WB off at the instant the star is depleted). All WB gate sites read this field
// each step, so flipping it takes effect on the next step.
void
gkyl_gr_euler_set_wb_disabled(const struct gkyl_wv_eqn *eqn, bool disabled)
{
  struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
  gr_euler->disable_well_balanced = disabled;
}


void
gkyl_gr_euler_flux(double gas_gamma, const double q[71], double flux[71])
{
  double v[71] = { 0.0 };
  gkyl_gr_euler_prim_vars(gas_gamma, q, v);
  double rho =  v[0];
  double vx = v[1];
  double vy = v[2];
  double vz = v[3];
  double p = v[4];

  double lapse = v[5];
  double shift_x = v[6];

  double spatial_metric[3][3];
  spatial_metric[0][0] = v[9]; spatial_metric[0][1] = v[10]; spatial_metric[0][2] = v[11];
  spatial_metric[1][0] = v[12]; spatial_metric[1][1] = v[13]; spatial_metric[1][2] = v[14];
  spatial_metric[2][0] = v[15]; spatial_metric[2][1] = v[16]; spatial_metric[2][2] = v[17];

  double spatial_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
    (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
    (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));
  
  bool in_excision_region = false;
  if (v[27] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }
  
  if (!in_excision_region) {
    double vel[3];
    double v_sq = 0.0;
    vel[0] = vx; vel[1] = vy; vel[2] = vz;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        v_sq += spatial_metric[i][j] * vel[i] * vel[j];
      }
    }

    double W = 1.0 / (sqrt(1.0 - v_sq));
    if (v_sq > 1.0 - pow(10.0, -8.0)) {
      W = 1.0 / sqrt(pow(10.0, -8.0));
    }

    double h = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0)));

    double cov_vel[3] = { 0.0 };
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        cov_vel[i] += spatial_metric[i][j] * vel[j];
      }
    }

    flux[0] = (lapse * sqrt(spatial_det)) * (rho * W * (vx - (shift_x / lapse)));
    flux[1] = (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (cov_vel[0] * (vx - (shift_x / lapse))) + p);
    flux[2] = (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (cov_vel[1] * (vx - (shift_x / lapse))));
    flux[3] = (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (cov_vel[2] * (vx - (shift_x / lapse))));
    flux[4] = (lapse * sqrt(spatial_det)) * (((rho * h * (W * W)) - p - (rho * W)) * (vx - (shift_x / lapse)) + (p * vx));

    for (int i = 5; i < 71; i++) {
      flux[i] = 0.0;
    }
  }
  else {
    for (int i = 0; i < 71; i++) {
      flux[i] = 0.0;
    }
  }
}

// [BARO] Barotropic fallback in the primitive recovery (private diagnostic library, 2026-11-06; AUDIT_FLUID_STATE_WRITES,
// "BH formation and hold run", amendment). GKYL_C2P_BARO_K=K' > 0 enables it: where the gamma-law recovery has no valid
// root (the atmosphere-floor test below fails: rho or p at/below the floor trigger, v^2 >= 1 or non-finite) and D > 0,
// rho and v^i are recovered from D and S_i on the adiabat p = K' rho^Gamma instead (the "Font fix" fallback, cf.
// IllinoisGRMHD): h(D/W) sqrt(W^2 - 1) = |S|/D has a unique root W in [1, sqrt(1 + |S|^2/D^2)] for every D > 0 (the
// left side increases from 0). Used only where that density is above the atmosphere band (100 rho_floor, the
// validator's band); below it the static atmosphere applies as before. Read-only on q: the stored D, S_i and tau are
// never written. Unset: the production recovery, bit for bit.
static double
gr_euler_baro_k(void)
{
  static double k = -1.0;
  if (k < 0.0) {
    const char *e = getenv("GKYL_C2P_BARO_K");
    k = e ? atof(e) : 0.0;
    if (!(k > 0.0)) k = 0.0;
  }
  return k;
}

static bool
gr_euler_baro_recover(double gas_gamma, double K, double D, double S2, double *rho, double *p, double *Q)
{
  if (!(D > 0.0) || !isfinite(D) || !(S2 >= 0.0) || !isfinite(S2)) return false;
  double s = sqrt(S2) / D, lo = 1.0, hi = sqrt(1.0 + s * s);
  double c = gas_gamma / (gas_gamma - 1.0) * K;
  for (int it = 0; it < 200 && hi - lo > 1.0e-15 * hi; it++) {
    double W = 0.5 * (lo + hi);
    double f = (1.0 + c * pow(D / W, gas_gamma - 1.0)) * sqrt(W * W - 1.0) - s;
    if (f > 0.0) hi = W; else lo = W;
  }
  double W = 0.5 * (lo + hi);
  *rho = D / W;
  *p = K * pow(*rho, gas_gamma);
  *Q = (*rho) * (1.0 + c * pow(*rho, gas_gamma - 1.0)) * W * W; // rho h W^2
  return isfinite(*rho) && isfinite(*p) && isfinite(*Q) && *Q > 0.0;
}

// [BARO] recovery outcome per cell (gr_euler_prim_vars_impl's last argument).
#define GR_EULER_C2P_ROOT 0  // Gamma-law root (or no recovery: excision region, no atmosphere floor)
#define GR_EULER_C2P_BARO 1  // no valid Gamma-law root; barotropic fallback, live density above 100 rho_atm
#define GR_EULER_C2P_ATM 2   // no valid Gamma-law root; static atmosphere, live barotropic density <= 100 rho_atm
#define GR_EULER_C2P_FLOOR 3 // no valid Gamma-law root; static atmosphere otherwise (fallback off or failed)

static void
gr_euler_prim_vars_impl(double gas_gamma, const double q[71], double v[71], int *outcome)
{
  if (outcome) *outcome = GR_EULER_C2P_ROOT;
  double lapse = q[5];
  double shift_x = q[6];
  double shift_y = q[7];
  double shift_z = q[8];

  double spatial_metric[3][3];
  spatial_metric[0][0] = q[9]; spatial_metric[0][1] = q[10]; spatial_metric[0][2] = q[11];
  spatial_metric[1][0] = q[12]; spatial_metric[1][1] = q[13]; spatial_metric[1][2] = q[14];
  spatial_metric[2][0] = q[15]; spatial_metric[2][1] = q[16]; spatial_metric[2][2] = q[17];
  
  double extrinsic_curvature[3][3];
  extrinsic_curvature[0][0] = q[18]; extrinsic_curvature[0][1] = q[19]; extrinsic_curvature[0][2] = q[20];
  extrinsic_curvature[1][0] = q[21]; extrinsic_curvature[1][1] = q[22]; extrinsic_curvature[1][2] = q[23];
  extrinsic_curvature[2][0] = q[24]; extrinsic_curvature[2][1] = q[25]; extrinsic_curvature[2][2] = q[26];

  double lapse_der[3];
  lapse_der[0] = q[28];
  lapse_der[1] = q[29];
  lapse_der[2] = q[30];

  double shift_der[3][3];
  shift_der[0][0] = q[31]; shift_der[0][1] = q[32]; shift_der[0][2] = q[33];
  shift_der[1][0] = q[34]; shift_der[1][1] = q[35]; shift_der[1][2] = q[36];
  shift_der[2][0] = q[37]; shift_der[2][1] = q[38]; shift_der[2][2] = q[39];

  double spatial_metric_der[3][3][3];
  spatial_metric_der[0][0][0] = q[40]; spatial_metric_der[0][0][1] = q[41]; spatial_metric_der[0][0][2] = q[42];
  spatial_metric_der[0][1][0] = q[43]; spatial_metric_der[0][1][1] = q[44]; spatial_metric_der[0][1][2] = q[45];
  spatial_metric_der[0][2][0] = q[46]; spatial_metric_der[0][2][1] = q[47]; spatial_metric_der[0][2][2] = q[48];

  spatial_metric_der[1][0][0] = q[49]; spatial_metric_der[1][0][1] = q[50]; spatial_metric_der[1][0][2] = q[51];
  spatial_metric_der[1][1][0] = q[52]; spatial_metric_der[1][1][1] = q[53]; spatial_metric_der[1][1][2] = q[54];
  spatial_metric_der[1][2][0] = q[55]; spatial_metric_der[1][2][1] = q[56]; spatial_metric_der[1][2][2] = q[57];

  spatial_metric_der[2][0][0] = q[58]; spatial_metric_der[2][0][1] = q[59]; spatial_metric_der[2][0][2] = q[60];
  spatial_metric_der[2][1][0] = q[61]; spatial_metric_der[2][1][1] = q[62]; spatial_metric_der[2][1][2] = q[63];
  spatial_metric_der[2][2][0] = q[64]; spatial_metric_der[2][2][1] = q[65]; spatial_metric_der[2][2][2] = q[66];

  double evol_param = q[67];
  double x = q[68];
  double y = q[69];
  double z = q[70];

  bool in_excision_region = false;
  if (q[27] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double spatial_det = gr_euler_spatial_metric_det(spatial_metric);
    double inv_spatial_metric[3][3];
    gr_euler_spatial_metric_inv(spatial_metric, spatial_det, inv_spatial_metric);

    double D = q[0] / sqrt(spatial_det);
    double momx = q[1] / sqrt(spatial_det); // S_x
    double momy = q[2] / sqrt(spatial_det); // S_y
    double momz = q[3] / sqrt(spatial_det); // S_z
    double Etot = q[4] / sqrt(spatial_det);

    // Stored momentum is covariant Valencia S_i. Its physical norm is gamma^{ij} S_i S_j.
    double mom_cov[3] = { momx, momy, momz };
    double S2 = 0.0;
    for (int si = 0; si < 3; si++)
      for (int sj = 0; sj < 3; sj++)
        S2 += inv_spatial_metric[si][sj] * mom_cov[si] * mom_cov[sj];

    // C2P normalisation. T2 = (Etot + D)^2 - S2 = (rho h W^2)^2 - S2 scales like D^2, so it must be guarded relative to (Etot + D)^2, never against an absolute constant.
    static int c2p_legacy = -1;
    if (getenv("GKYL_C2P_LEGACY_GUARD") != NULL) {
        c2p_legacy = 1;
      }
      else {
        c2p_legacy = 0;
      }

    double T_tot = Etot + D;
    double T2 = (T_tot * T_tot) - S2;
    if (c2p_legacy) {
      if (T2 < pow(10.0, -8.0)) {
        T2 = pow(10.0, -8.0);
      }
    }
    else {
      // Relative floor: S2 <= T2 physically (v^2 < 1), so T2 -> 0 only as v -> c.
      double T2_min = 1.0e-14 * T_tot * T_tot; //a small relative margin so the sqrt stays finite without rescaling the state.
      if (!(T2 > T2_min)) {
        T2 = T2_min;
      }
    }

    double C = D / sqrt(T2);
    double C0 = (D + Etot) / sqrt(T2);

    double alpha0 = -1.0 / (gas_gamma * gas_gamma);
    double alpha1 = -2.0 * C * ((gas_gamma - 1.0) / (gas_gamma * gas_gamma));
    double alpha2 = ((gas_gamma - 2.0) / gas_gamma) * ((C0 * C0) - 1.0) + 1.0 - (C * C) * ((gas_gamma - 1.0) / gas_gamma) * ((gas_gamma - 1.0) / gas_gamma);
    double alpha4 = (C0 * C0) - 1.0;
    double eta = 2.0 * C *((gas_gamma - 1.0) / gas_gamma);

    double guess = 1.0;
    int iter = 0;
    while (iter < 100) {
      double poly = (alpha4 * (guess * guess * guess) * (guess - eta)) + (alpha2 * (guess * guess)) + (alpha1 * guess) + alpha0;
      double poly_der = alpha1 + (2.0 * alpha2 * guess) + (4.0 * alpha4 * (guess * guess * guess)) - (3.0 * eta * alpha4 * (guess * guess));

      double guess_new = guess - (poly / poly_der);
      double change = fabs(guess - guess_new);
      // Keep the converged Newton iterate, not the previous approximation.
      guess = guess_new;
      if (change < pow(10.0, -8.0)) break;
      iter += 1;
    }

    double W = 0.5 * C0 * guess * (1.0 + sqrt(1.0 + (4.0 * ((gas_gamma - 1.0) / gas_gamma) * ((1.0 - (C * guess)) / ((C0 * C0) * (guess * guess))))));
    double h = 1.0 / (C * guess);

    v[0] = D / W;
    v[4] = (v[0] * h * (W * W)) - D - Etot;
    double Q = Etot + D + v[4];
    v[1] = 0.0;
    v[2] = 0.0;
    v[3] = 0.0;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        v[1 + i] += inv_spatial_metric[i][j] * mom_cov[j] / Q;
      }
    }
    // Atmosphere floor (read-only on q, so conservation is exact: the flux and source the cell reports change, the
    // update still telescopes). If the recovered primitive state is below the floor or non-physical (rho/p <= floor,
    // superluminal, or non-finite), the fluxes and sources see a static atmosphere (rho_atm, p_atm, v = 0).
    if (gr_euler_atm_floor) {
      static double atm_tol = -1.0;
      if (atm_tol < 0.0) {
        const char *tol_env = getenv("GKYL_ATM_TOL");
        if (tol_env) {
          atm_tol = atof(tol_env);
        }
        else {
          atm_tol = 0.0;
        }
      }
      double rho_trigger = gr_euler_rho_floor * (1.0 + atm_tol);
      double p_trigger = gr_euler_p_floor * (1.0 + atm_tol);

      double vsq_chk = 0.0;
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          vsq_chk += spatial_metric[i][j] * v[1 + i] * v[1 + j];
      if (!(v[0] > rho_trigger) || !(v[4] > p_trigger) || !(vsq_chk < 1.0) || !isfinite(vsq_chk)) {
        double baro_K = gr_euler_baro_k(), baro_rho = 0.0, baro_p = 0.0, baro_Q = 0.0;
        bool baro_ok = baro_K > 0.0 && gr_euler_baro_recover(gas_gamma, baro_K, D, S2, &baro_rho, &baro_p, &baro_Q);
        if (baro_ok && baro_rho > 100.0 * gr_euler_rho_floor) {
          v[0] = baro_rho;
          v[4] = baro_p;
          for (int i = 0; i < 3; i++) {
            v[1 + i] = 0.0;
            for (int j = 0; j < 3; j++) {
              v[1 + i] += inv_spatial_metric[i][j] * mom_cov[j] / baro_Q;
            }
          }
          if (outcome) *outcome = GR_EULER_C2P_BARO;
        }
        else {
          v[0] = gr_euler_rho_floor;
          v[4] = gr_euler_p_floor;
          v[1] = 0.0; v[2] = 0.0; v[3] = 0.0;
          if (outcome) *outcome = baro_ok ? GR_EULER_C2P_ATM : GR_EULER_C2P_FLOOR;
        }
      }
    }
    else {
      // NaN-safe floor
      if (!(v[0] > gr_euler_rho_floor)) {
        v[0] = gr_euler_rho_floor;
      }
      if (!(v[4] > gr_euler_p_floor)) {
        v[4] = gr_euler_p_floor;
      }
      if (!isfinite(v[1])) { v[1] = 0.0; }
      if (!isfinite(v[2])) { v[2] = 0.0; }
      if (!isfinite(v[3])) { v[3] = 0.0; }
    }

    v[5] = lapse;
    v[6] = shift_x;
    v[7] = shift_y;
    v[8] = shift_z;

    v[9] = spatial_metric[0][0]; v[10] = spatial_metric[0][1]; v[11] = spatial_metric[0][2];
    v[12] = spatial_metric[1][0]; v[13] = spatial_metric[1][1]; v[14] = spatial_metric[1][2];
    v[15] = spatial_metric[2][0]; v[16] = spatial_metric[2][1]; v[17] = spatial_metric[2][2];

    v[18] = extrinsic_curvature[0][0]; v[19] = extrinsic_curvature[0][1]; v[20] = extrinsic_curvature[0][2];
    v[21] = extrinsic_curvature[1][0]; v[22] = extrinsic_curvature[1][1]; v[23] = extrinsic_curvature[1][2];
    v[24] = extrinsic_curvature[2][0]; v[25] = extrinsic_curvature[2][1]; v[26] = extrinsic_curvature[2][2];

    v[27] = 1.0;

    v[28] = lapse_der[0];
    v[29] = lapse_der[1];
    v[30] = lapse_der[2];

    v[31] = shift_der[0][0]; v[32] = shift_der[0][1]; v[33] = shift_der[0][2];
    v[34] = shift_der[1][0]; v[35] = shift_der[1][1]; v[36] = shift_der[1][2];
    v[37] = shift_der[2][0]; v[38] = shift_der[2][1]; v[39] = shift_der[2][2];

    v[40] = spatial_metric_der[0][0][0]; v[41] = spatial_metric_der[0][0][1]; v[42] = spatial_metric_der[0][0][2];
    v[43] = spatial_metric_der[0][1][0]; v[44] = spatial_metric_der[0][1][1]; v[45] = spatial_metric_der[0][1][2];
    v[46] = spatial_metric_der[0][2][0]; v[47] = spatial_metric_der[0][2][1]; v[48] = spatial_metric_der[0][2][2];

    v[49] = spatial_metric_der[1][0][0]; v[50] = spatial_metric_der[1][0][1]; v[51] = spatial_metric_der[1][0][2];
    v[52] = spatial_metric_der[1][1][0]; v[53] = spatial_metric_der[1][1][1]; v[54] = spatial_metric_der[1][1][2];
    v[55] = spatial_metric_der[1][2][0]; v[56] = spatial_metric_der[1][2][1]; v[57] = spatial_metric_der[1][2][2];

    v[58] = spatial_metric_der[2][0][0]; v[59] = spatial_metric_der[2][0][1]; v[60] = spatial_metric_der[2][0][2];
    v[61] = spatial_metric_der[2][1][0]; v[62] = spatial_metric_der[2][1][1]; v[63] = spatial_metric_der[2][1][2];
    v[64] = spatial_metric_der[2][2][0]; v[65] = spatial_metric_der[2][2][1]; v[66] = spatial_metric_der[2][2][2];

    v[67] = evol_param;
    v[68] = x;
    v[69] = y;
    v[70] = z;
  }
  else {
    for (int i = 0; i < 71; i++) {
      v[i] = 0.0;
    }
    
    v[27] = -1.0;
  }
}

void
gkyl_gr_euler_prim_vars(double gas_gamma, const double q[71], double v[71])
{
  gr_euler_prim_vars_impl(gas_gamma, q, v, NULL);
}

static int
gr_euler_c2p_outcome(double gas_gamma, const double q[71])
{
  double v[71];
  int outcome = GR_EULER_C2P_ROOT;
  gr_euler_prim_vars_impl(gas_gamma, q, v, &outcome);
  return outcome;
}

// [BARO] true iff GKYL_C2P_BARO_K is set and q's recovery takes the barotropic fallback (the route-1 [BARO] cell count).
bool
gkyl_gr_euler_c2p_baro_used(double gas_gamma, const double q[71])
{
  if (!(gr_euler_baro_k() > 0.0)) return false;
  return gr_euler_c2p_outcome(gas_gamma, q) == GR_EULER_C2P_BARO;
}

// [BARO] true iff GKYL_C2P_BARO_K is set, q has no valid Gamma-law root, and the barotropic branch gives the cell a
// state: the fallback (live density above 100 rho_atm) or the static atmosphere (live density at or below). The
// route-1 step check and FOFC's test call this, so their atmosphere definition is the recovery's own.
bool
gkyl_gr_euler_c2p_baro_handled(double gas_gamma, const double q[71])
{
  if (!(gr_euler_baro_k() > 0.0)) return false;
  int outcome = gr_euler_c2p_outcome(gas_gamma, q);
  return outcome == GR_EULER_C2P_BARO || outcome == GR_EULER_C2P_ATM;
}

void 
gkyl_gr_euler_inv_spatial_metric(const double q[71], double ***inv_spatial_metric)
{
  double spatial_metric[3][3];
  spatial_metric[0][0] = q[9]; spatial_metric[0][1] = q[10]; spatial_metric[0][2] = q[11];
  spatial_metric[1][0] = q[12]; spatial_metric[1][1] = q[13]; spatial_metric[1][2] = q[14];
  spatial_metric[2][0] = q[15]; spatial_metric[2][1] = q[16]; spatial_metric[2][2] = q[17];

  double spatial_det = gr_euler_spatial_metric_det(spatial_metric);
  double inv_local[3][3];
  gr_euler_spatial_metric_inv(spatial_metric, spatial_det, inv_local);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      (*inv_spatial_metric)[i][j] = inv_local[i][j];
    }
  }
}

void
gkyl_gr_euler_stress_energy_tensor(double gas_gamma, const double q[71], double ***stress_energy)
{
  double v[71] = { 0.0 };
  gkyl_gr_euler_prim_vars(gas_gamma, q, v);
  double rho = v[0];
  double vx = v[1];
  double vy = v[2];
  double vz = v[3];
  double p = v[4];

  double lapse = v[5];
  double shift_x = v[6];
  double shift_y = v[7];
  double shift_z = v[8];

  double spatial_metric[3][3];
  spatial_metric[0][0] = v[9]; spatial_metric[0][1] = v[10]; spatial_metric[0][2] = v[11];
  spatial_metric[1][0] = v[12]; spatial_metric[1][1] = v[13]; spatial_metric[1][2] = v[14];
  spatial_metric[2][0] = v[15]; spatial_metric[2][1] = v[16]; spatial_metric[2][2] = v[17];

  double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  gkyl_gr_euler_inv_spatial_metric(q, &inv_spatial_metric);

  bool in_excision_region = false;
  if (v[27] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double vel[3];
    double v_sq = 0.0;
    vel[0] = vx; vel[1] = vy; vel[2] = vz;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        v_sq += spatial_metric[i][j] * vel[i] * vel[j];
      }
    }

    double W = 1.0 / sqrt(1.0 - v_sq);
    if (v_sq > 1.0 - pow(10.0, -8.0)) {
      W = 1.0 / sqrt(pow(10.0, -8.0));
    }

    double h = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0)));

    double spacetime_vel[4];
    spacetime_vel[0] = W / lapse;
    spacetime_vel[1] = (W * vx) - (shift_x * (W / lapse));
    spacetime_vel[2] = (W * vy) - (shift_y * (W / lapse));
    spacetime_vel[3] = (W * vz) - (shift_z * (W / lapse));

    double shift[3];
    shift[0] = shift_x; shift[1] = shift_y; shift[2] = shift_z;

    double inv_spacetime_metric[4][4];
    inv_spacetime_metric[0][0] = - (1.0 / (lapse * lapse));
    for (int i = 0; i < 3; i++) {
      inv_spacetime_metric[0][i + 1] = (1.0 / (lapse * lapse)) * shift[i];
      inv_spacetime_metric[i + 1][0] = (1.0 / (lapse * lapse)) * shift[i];
    }
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        inv_spacetime_metric[i + 1][j + 1] = inv_spatial_metric[i][j] - ((1.0 / (lapse * lapse)) * shift[i] * shift[j]);
      }
    }

    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        (*stress_energy)[i][j] = (rho * h * spacetime_vel[i] * spacetime_vel[j]) + (p * inv_spacetime_metric[i][j]);
      }
    }
  }
  else {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        (*stress_energy)[i][j] = 0.0;
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(inv_spatial_metric[i]);
  }
  gkyl_free(inv_spatial_metric);
}

static inline double
gkyl_gr_euler_max_abs_speed(double gas_gamma, const double q[71])
{
  double v[71] = { 0.0 };
  gkyl_gr_euler_prim_vars(gas_gamma, q, v);
  double rho =  v[0];
  double vx = v[1];
  double vy = v[2];
  double vz = v[3];
  double p = v[4];

  double lapse = v[5];
  double shift_x = v[6];
  double shift_y = v[7];
  double shift_z = v[8];

  double spatial_metric[3][3];
  spatial_metric[0][0] = v[9]; spatial_metric[0][1] = v[10]; spatial_metric[0][2] = v[11];
  spatial_metric[1][0] = v[12]; spatial_metric[1][1] = v[13]; spatial_metric[1][2] = v[14];
  spatial_metric[2][0] = v[15]; spatial_metric[2][1] = v[16]; spatial_metric[2][2] = v[17];

  double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  gkyl_gr_euler_inv_spatial_metric(q, &inv_spatial_metric);

  double num = (gas_gamma * p) / rho;
  double den = 1.0 + ((p / rho) * (gas_gamma) / (gas_gamma - 1.0));
  double c_s = sqrt(num / den);

  bool in_excision_region = false;
  if (v[27] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  bool curved_spacetime = false;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j) {
        if (fabs(spatial_metric[i][j] - 1.0) > pow(10.0, -8.0)) {
          curved_spacetime = true;
        }
      }
      else {
        if (fabs(spatial_metric[i][j]) > pow(10.0, -8.0)) {
          curved_spacetime = true;
        }
      }
    }
  }
  if (fabs(lapse - 1.0) > pow(10.0, -8.0) || fabs(shift_x) > pow(10.0, -8.0) || fabs(shift_y) > pow(10.0, -8.0) ||
    fabs(shift_z) > pow(10.0, -8.0)) {
    curved_spacetime = true;
  }

  if (!in_excision_region) {
    if (curved_spacetime) {
      double vel[3];
      double v_sq = 0.0;
      vel[0] = vx; vel[1] = vy; vel[2] = vz;

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          v_sq += spatial_metric[i][j] * vel[i] * vel[j];
        }
      }

      double shift[3];
      shift[0] = shift_x; shift[1] = shift_y; shift[2] = shift_z;

      double material_eigs[3];
      double fast_acoustic_eigs[3];
      double slow_acoustic_eigs[3];

      for (int i = 0; i < 3; i++) {
        material_eigs[i] = (lapse * vel[i]) - shift[i];

        fast_acoustic_eigs[i] = (lapse / (1.0 - (v_sq * (c_s * c_s)))) * ((vel[i] * (1.0 - (c_s * c_s))) +
          (c_s * sqrt((1.0 - v_sq) * (inv_spatial_metric[i][i] * (1.0 - (v_sq * (c_s * c_s))) - (vel[i] * vel[i]) * (1.0 - (c_s * c_s)))))) - shift[i];
        
        slow_acoustic_eigs[i] = (lapse / (1.0 - (v_sq * (c_s * c_s)))) * ((vel[i] * (1.0 - (c_s * c_s))) -
          (c_s * sqrt((1.0 - v_sq) * (inv_spatial_metric[i][i] * (1.0 - (v_sq * (c_s * c_s))) - (vel[i] * vel[i]) * (1.0 - (c_s * c_s)))))) - shift[i];
      }

      double max_eig = 0.0;
      for (int i = 0; i < 3; i++) {
        if (fabs(material_eigs[i]) > max_eig) {
          max_eig = fabs(material_eigs[i]);
        }
        if (fabs(fast_acoustic_eigs[i]) > max_eig) {
          max_eig = fabs(fast_acoustic_eigs[i]);
        }
        if (fabs(slow_acoustic_eigs[i]) > max_eig) {
          max_eig = fabs(slow_acoustic_eigs[i]);
        }
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(inv_spatial_metric[i]);
      }
      gkyl_free(inv_spatial_metric);

      return fabs(v_sq) + max_eig;
    }
    else {
      double v_sq = sqrt((vx * vx) + (vy * vy) + (vz * vz));

      for (int i = 0; i < 3; i++) {
        gkyl_free(inv_spatial_metric[i]);
      }
      gkyl_free(inv_spatial_metric);

      return fabs(v_sq) + c_s;
    }
  }
  else {
    for (int i = 0; i < 3; i++) {
      gkyl_free(inv_spatial_metric[i]);
    }
    gkyl_free(inv_spatial_metric);

    return pow(10.0, -8.0);
  }
}

static inline void
cons_to_riem(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* qin, double* wout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < eqn->num_equations; i++) {
    wout[i] = qin[i];
  }
}

static inline void
riem_to_cons(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* win, double* qout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < eqn->num_equations; i++) {
    qout[i] = win[i];
  }
}

static void
gr_euler_wall(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 0; i < eqn->num_equations; i++) {
    ghost[i] = skin[i];
  }

  ghost[1] = -ghost[1];
}

static void
gr_euler_no_slip(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 0; i < eqn->num_equations; i++) {
    ghost[i] = skin[i];
  }

  for (int i = 1; i < 4; i++) {
    ghost[i] = -skin[i];
  }
}

static inline void
rot_to_local(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qglobal,
  double* GKYL_RESTRICT qlocal)
{
  qlocal[0] = qglobal[0];
  qlocal[1] = (qglobal[1] * norm[0]) + (qglobal[2] * norm[1]) + (qglobal[3] * norm[2]);
  qlocal[2] = (qglobal[1] * tau1[0]) + (qglobal[2] * tau1[1]) + (qglobal[3] * tau1[2]);
  qlocal[3] = (qglobal[1] * tau2[0]) + (qglobal[2] * tau2[1]) + (qglobal[3] * tau2[2]);
  qlocal[4] = qglobal[4];

  qlocal[5] = qglobal[5];
  qlocal[6] = (qglobal[6] * norm[0]) + (qglobal[7] * norm[1]) + (qglobal[8] * norm[2]);
  qlocal[7] = (qglobal[6] * tau1[0]) + (qglobal[7] * tau1[1]) + (qglobal[8] * tau1[2]);
  qlocal[8] = (qglobal[6] * tau2[0]) + (qglobal[7] * tau2[1]) + (qglobal[8] * tau2[2]);

  // Temporary arrays to store rotated column vectors.
  double r1[3], r2[3], r3[3];
  r1[0] = (qglobal[9] * norm[0]) + (qglobal[10] * norm[1]) + (qglobal[11] * norm[2]);
  r1[1] = (qglobal[9] * tau1[0]) + (qglobal[10] * tau1[1]) + (qglobal[11] * tau1[2]);
  r1[2] = (qglobal[9] * tau2[0]) + (qglobal[10] * tau2[1]) + (qglobal[11] * tau2[2]);

  r2[0] = (qglobal[12] * norm[0]) + (qglobal[13] * norm[1]) + (qglobal[14] * norm[2]);
  r2[1] = (qglobal[12] * tau1[0]) + (qglobal[13] * tau1[1]) + (qglobal[14] * tau1[2]);
  r2[2] = (qglobal[12] * tau2[0]) + (qglobal[13] * tau2[1]) + (qglobal[14] * tau2[2]);

  r3[0] = (qglobal[15] * norm[0]) + (qglobal[16] * norm[1]) + (qglobal[17] * norm[2]);
  r3[1] = (qglobal[15] * tau1[0]) + (qglobal[16] * tau1[1]) + (qglobal[17] * tau1[2]);
  r3[2] = (qglobal[15] * tau2[0]) + (qglobal[16] * tau2[1]) + (qglobal[17] * tau2[2]);

  // Temporary arrays to store rotated row vectors.
  double v1[3], v2[3], v3[3];
  v1[0] = (r1[0] * norm[0]) + (r2[0] * norm[1]) + (r3[0] * norm[2]);
  v1[1] = (r1[0] * tau1[0]) + (r2[0] * tau1[1]) + (r3[0] * tau1[2]);
  v1[2] = (r1[0] * tau2[0]) + (r2[0] * tau2[1]) + (r3[0] * tau2[2]);

  v2[0] = (r1[1] * norm[0]) + (r2[1] * norm[1]) + (r3[1] * norm[2]);
  v2[1] = (r1[1] * tau1[0]) + (r2[1] * tau1[1]) + (r3[1] * tau1[2]);
  v2[2] = (r1[1] * tau2[0]) + (r2[1] * tau2[1]) + (r3[1] * tau2[2]);

  v3[0] = (r1[2] * norm[0]) + (r2[2] * norm[1]) + (r3[2] * norm[2]);
  v3[1] = (r1[2] * tau1[0]) + (r2[2] * tau1[1]) + (r3[2] * tau1[2]);
  v3[2] = (r1[2] * tau2[0]) + (r2[2] * tau2[1]) + (r3[2] * tau2[2]);

  // Rotate spatial metric tensor to local coordinate frame.
  qlocal[9] = v1[0]; qlocal[10] = v1[1]; qlocal[11] = v1[2];
  qlocal[12] = v2[0]; qlocal[13] = v2[1]; qlocal[14] = v2[2];
  qlocal[15] = v3[0]; qlocal[16] = v3[1]; qlocal[17] = v3[2];

  // Temporary arrays to store rotated extrinsic column vectors.
  double extr_r1[3], extr_r2[3], extr_r3[3];
  extr_r1[0] = (qglobal[18] * norm[0]) + (qglobal[19] * norm[1]) + (qglobal[20] * norm[2]);
  extr_r1[1] = (qglobal[18] * tau1[0]) + (qglobal[19] * tau1[1]) + (qglobal[20] * tau1[2]);
  extr_r1[2] = (qglobal[18] * tau2[0]) + (qglobal[19] * tau2[1]) + (qglobal[20] * tau2[2]);

  extr_r2[0] = (qglobal[21] * norm[0]) + (qglobal[22] * norm[1]) + (qglobal[23] * norm[2]);
  extr_r2[1] = (qglobal[21] * tau1[0]) + (qglobal[22] * tau1[1]) + (qglobal[23] * tau1[2]);
  extr_r2[2] = (qglobal[21] * tau2[0]) + (qglobal[22] * tau2[1]) + (qglobal[23] * tau2[2]);

  extr_r3[0] = (qglobal[24] * norm[0]) + (qglobal[25] * norm[1]) + (qglobal[26] * norm[2]);
  extr_r3[1] = (qglobal[24] * tau1[0]) + (qglobal[25] * tau1[1]) + (qglobal[26] * tau1[2]);
  extr_r3[2] = (qglobal[24] * tau2[0]) + (qglobal[25] * tau2[1]) + (qglobal[26] * tau2[2]);

  // Temporary arrays to store rotated extrinsic row vectors.
  double inv_v1[3], inv_v2[3], inv_v3[3];
  inv_v1[0] = (extr_r1[0] * norm[0]) + (extr_r2[0] * norm[1]) + (extr_r3[0] * norm[2]);
  inv_v1[1] = (extr_r1[0] * tau1[0]) + (extr_r2[0] * tau1[1]) + (extr_r3[0] * tau1[2]);
  inv_v1[2] = (extr_r1[0] * tau2[0]) + (extr_r2[0] * tau2[1]) + (extr_r3[0] * tau2[2]);

  inv_v2[0] = (extr_r1[1] * norm[0]) + (extr_r2[1] * norm[1]) + (extr_r3[1] * norm[2]);
  inv_v2[1] = (extr_r1[1] * tau1[0]) + (extr_r2[1] * tau1[1]) + (extr_r3[1] * tau1[2]);
  inv_v2[2] = (extr_r1[1] * tau2[0]) + (extr_r2[1] * tau2[1]) + (extr_r3[1] * tau2[2]);

  inv_v3[0] = (extr_r1[2] * norm[0]) + (extr_r2[2] * norm[1]) + (extr_r3[2] * norm[2]);
  inv_v3[1] = (extr_r1[2] * tau1[0]) + (extr_r2[2] * tau1[1]) + (extr_r3[2] * tau1[2]);
  inv_v3[2] = (extr_r1[2] * tau2[0]) + (extr_r2[2] * tau2[1]) + (extr_r3[2] * tau2[2]);

  // Rotate extrinsic curvature tensor to local coordinate frame.
  qlocal[18] = inv_v1[0]; qlocal[19] = inv_v1[1]; qlocal[20] = inv_v1[2];
  qlocal[21] = inv_v2[0]; qlocal[22] = inv_v2[1]; qlocal[23] = inv_v2[2];
  qlocal[24] = inv_v3[0]; qlocal[25] = inv_v3[1]; qlocal[26] = inv_v3[2];

  qlocal[27] = qglobal[27];

  qlocal[28] = (qglobal[28] * norm[0]) + (qglobal[29] * norm[1]) + (qglobal[30] * norm[2]);
  qlocal[29] = (qglobal[28] * tau1[0]) + (qglobal[29] * tau1[1]) + (qglobal[30] * tau1[2]);
  qlocal[30] = (qglobal[28] * tau2[0]) + (qglobal[29] * tau2[1]) + (qglobal[30] * tau2[2]);

  // Temporary arrays to store rotated shift derivative column vectors.
  double shiftder_r1[3], shiftder_r2[3], shiftder_r3[3];
  shiftder_r1[0] = (qglobal[31] * norm[0]) + (qglobal[32] * norm[1]) + (qglobal[33] * norm[2]);
  shiftder_r1[1] = (qglobal[31] * tau1[0]) + (qglobal[32] * tau1[1]) + (qglobal[33] * tau1[2]);
  shiftder_r1[2] = (qglobal[31] * tau2[0]) + (qglobal[32] * tau2[1]) + (qglobal[33] * tau2[2]);

  shiftder_r2[0] = (qglobal[34] * norm[0]) + (qglobal[35] * norm[1]) + (qglobal[36] * norm[2]);
  shiftder_r2[1] = (qglobal[34] * tau1[0]) + (qglobal[35] * tau1[1]) + (qglobal[36] * tau1[2]);
  shiftder_r2[2] = (qglobal[34] * tau2[0]) + (qglobal[35] * tau2[1]) + (qglobal[36] * tau2[2]);

  shiftder_r3[0] = (qglobal[37] * norm[0]) + (qglobal[38] * norm[1]) + (qglobal[39] * norm[2]);
  shiftder_r3[1] = (qglobal[37] * tau1[0]) + (qglobal[38] * tau1[1]) + (qglobal[39] * tau1[2]);
  shiftder_r3[2] = (qglobal[37] * tau2[0]) + (qglobal[38] * tau2[1]) + (qglobal[39] * tau2[2]);

  // Temporary arrays to store rotated shift derivative row vectors.
  double shiftder_v1[3], shiftder_v2[3], shiftder_v3[3];
  shiftder_v1[0] = (shiftder_r1[0] * norm[0]) + (shiftder_r2[0] * norm[1]) + (shiftder_r3[0] * norm[2]);
  shiftder_v1[1] = (shiftder_r1[0] * tau1[0]) + (shiftder_r2[0] * tau1[1]) + (shiftder_r3[0] * tau1[2]);
  shiftder_v1[2] = (shiftder_r1[0] * tau2[0]) + (shiftder_r2[0] * tau2[1]) + (shiftder_r3[0] * tau2[2]);

  shiftder_v2[0] = (shiftder_r1[1] * norm[0]) + (shiftder_r2[1] * norm[1]) + (shiftder_r3[1] * norm[2]);
  shiftder_v2[1] = (shiftder_r1[1] * tau1[0]) + (shiftder_r2[1] * tau1[1]) + (shiftder_r3[1] * tau1[2]);
  shiftder_v2[2] = (shiftder_r1[1] * tau2[0]) + (shiftder_r2[1] * tau2[1]) + (shiftder_r3[1] * tau2[2]);

  shiftder_v3[0] = (shiftder_r1[2] * norm[0]) + (shiftder_r2[2] * norm[1]) + (shiftder_r3[2] * norm[2]);
  shiftder_v3[1] = (shiftder_r1[2] * tau1[0]) + (shiftder_r2[2] * tau1[1]) + (shiftder_r3[2] * tau1[2]);
  shiftder_v3[2] = (shiftder_r1[2] * tau2[0]) + (shiftder_r2[2] * tau2[1]) + (shiftder_r3[2] * tau2[2]);

  // Rotate shift vector derivative to local coordinate frame.
  qlocal[31] = shiftder_v1[0]; qlocal[32] = shiftder_v1[1]; qlocal[33] = shiftder_v1[2];
  qlocal[34] = shiftder_v2[0]; qlocal[35] = shiftder_v2[1]; qlocal[36] = shiftder_v2[2];
  qlocal[37] = shiftder_v3[0]; qlocal[38] = shiftder_v3[1]; qlocal[39] = shiftder_v3[2];

  // Temporary arrays to store rotated column vectors.
  double r11[3], r12[3], r13[3];
  double r21[3], r22[3], r23[3];
  double r31[3], r32[3], r33[3];

  r11[0] = (qglobal[40] * norm[0]) + (qglobal[41] * norm[1]) + (qglobal[42] * norm[2]);
  r11[1] = (qglobal[40] * tau1[0]) + (qglobal[41] * tau1[1]) + (qglobal[42] * tau1[2]);
  r11[2] = (qglobal[40] * tau2[0]) + (qglobal[41] * tau2[1]) + (qglobal[42] * tau2[2]);

  r12[0] = (qglobal[43] * norm[0]) + (qglobal[44] * norm[1]) + (qglobal[45] * norm[2]);
  r12[1] = (qglobal[43] * tau1[0]) + (qglobal[44] * tau1[1]) + (qglobal[45] * tau1[2]);
  r12[2] = (qglobal[43] * tau2[0]) + (qglobal[44] * tau2[1]) + (qglobal[45] * tau2[2]);

  r13[0] = (qglobal[46] * norm[0]) + (qglobal[47] * norm[1]) + (qglobal[48] * norm[2]);
  r13[1] = (qglobal[46] * tau1[0]) + (qglobal[47] * tau1[1]) + (qglobal[48] * tau1[2]);
  r13[2] = (qglobal[46] * tau2[0]) + (qglobal[47] * tau2[1]) + (qglobal[48] * tau2[2]);

  r21[0] = (qglobal[49] * norm[0]) + (qglobal[50] * norm[1]) + (qglobal[51] * norm[2]);
  r21[1] = (qglobal[49] * tau1[0]) + (qglobal[50] * tau1[1]) + (qglobal[51] * tau1[2]);
  r21[2] = (qglobal[49] * tau2[0]) + (qglobal[50] * tau2[1]) + (qglobal[51] * tau2[2]);

  r22[0] = (qglobal[52] * norm[0]) + (qglobal[53] * norm[1]) + (qglobal[54] * norm[2]);
  r22[1] = (qglobal[52] * tau1[0]) + (qglobal[53] * tau1[1]) + (qglobal[54] * tau1[2]);
  r22[2] = (qglobal[52] * tau2[0]) + (qglobal[53] * tau2[1]) + (qglobal[54] * tau2[2]);

  r23[0] = (qglobal[55] * norm[0]) + (qglobal[56] * norm[1]) + (qglobal[57] * norm[2]);
  r23[1] = (qglobal[55] * tau1[0]) + (qglobal[56] * tau1[1]) + (qglobal[57] * tau1[2]);
  r23[2] = (qglobal[55] * tau2[0]) + (qglobal[56] * tau2[1]) + (qglobal[57] * tau2[2]);

  r31[0] = (qglobal[58] * norm[0]) + (qglobal[59] * norm[1]) + (qglobal[60] * norm[2]);
  r31[1] = (qglobal[58] * tau1[0]) + (qglobal[59] * tau1[1]) + (qglobal[60] * tau1[2]);
  r31[2] = (qglobal[58] * tau2[0]) + (qglobal[59] * tau2[1]) + (qglobal[60] * tau2[2]);

  r32[0] = (qglobal[61] * norm[0]) + (qglobal[62] * norm[1]) + (qglobal[63] * norm[2]);
  r32[1] = (qglobal[61] * tau1[0]) + (qglobal[62] * tau1[1]) + (qglobal[63] * tau1[2]);
  r32[2] = (qglobal[61] * tau2[0]) + (qglobal[62] * tau2[1]) + (qglobal[63] * tau2[2]);

  r33[0] = (qglobal[64] * norm[0]) + (qglobal[65] * norm[1]) + (qglobal[66] * norm[2]);
  r33[1] = (qglobal[64] * tau1[0]) + (qglobal[65] * tau1[1]) + (qglobal[66] * tau1[2]);
  r33[2] = (qglobal[64] * tau2[0]) + (qglobal[65] * tau2[1]) + (qglobal[66] * tau2[2]);

  // Temporary arrays to store rotated row vectors.
  double s11[3], s12[3], s13[3];
  double s21[3], s22[3], s23[3];
  double s31[3], s32[3], s33[3];

  s11[0] = (r11[0] * norm[0]) + (r12[0] * norm[1]) + (r13[0] * norm[2]);
  s11[1] = (r11[1] * norm[0]) + (r12[1] * norm[1]) + (r13[1] * norm[2]);
  s11[2] = (r11[2] * norm[0]) + (r12[2] * norm[1]) + (r13[2] * norm[2]);

  s12[0] = (r11[0] * tau1[0]) + (r12[0] * tau1[1]) + (r13[0] * tau1[2]);
  s12[1] = (r11[1] * tau1[0]) + (r12[1] * tau1[1]) + (r13[1] * tau1[2]);
  s12[2] = (r11[2] * tau1[0]) + (r12[2] * tau1[1]) + (r13[2] * tau1[2]);

  s13[0] = (r11[0] * tau2[0]) + (r12[0] * tau2[1]) + (r13[0] * tau2[2]);
  s13[1] = (r11[1] * tau2[0]) + (r12[1] * tau2[1]) + (r13[1] * tau2[2]);
  s13[2] = (r11[2] * tau2[0]) + (r12[2] * tau2[1]) + (r13[2] * tau2[2]);

  s21[0] = (r21[0] * norm[0]) + (r22[0] * norm[1]) + (r23[0] * norm[2]);
  s21[1] = (r21[1] * norm[0]) + (r22[1] * norm[1]) + (r23[1] * norm[2]);
  s21[2] = (r21[2] * norm[0]) + (r22[2] * norm[1]) + (r23[2] * norm[2]);

  s22[0] = (r21[0] * tau1[0]) + (r22[0] * tau1[1]) + (r23[0] * tau1[2]);
  s22[1] = (r21[1] * tau1[0]) + (r22[1] * tau1[1]) + (r23[1] * tau1[2]);
  s22[2] = (r21[2] * tau1[0]) + (r22[2] * tau1[1]) + (r23[2] * tau1[2]);

  s23[0] = (r21[0] * tau2[0]) + (r22[0] * tau2[1]) + (r23[0] * tau2[2]);
  s23[1] = (r21[1] * tau2[0]) + (r22[1] * tau2[1]) + (r23[1] * tau2[2]);
  s23[2] = (r21[2] * tau2[0]) + (r22[2] * tau2[1]) + (r23[2] * tau2[2]);

  s31[0] = (r31[0] * norm[0]) + (r32[0] * norm[1]) + (r33[0] * norm[2]);
  s31[1] = (r31[1] * norm[0]) + (r32[1] * norm[1]) + (r33[1] * norm[2]);
  s31[2] = (r31[2] * norm[0]) + (r32[2] * norm[1]) + (r33[2] * norm[2]);

  s32[0] = (r31[0] * tau1[0]) + (r32[0] * tau1[1]) + (r33[0] * tau1[2]);
  s32[1] = (r31[1] * tau1[0]) + (r32[1] * tau1[1]) + (r33[1] * tau1[2]);
  s32[2] = (r31[2] * tau1[0]) + (r32[2] * tau1[1]) + (r33[2] * tau1[2]);

  s33[0] = (r31[0] * tau2[0]) + (r32[0] * tau2[1]) + (r33[0] * tau2[2]);
  s33[1] = (r31[1] * tau2[0]) + (r32[1] * tau2[1]) + (r33[1] * tau2[2]);
  s33[2] = (r31[2] * tau2[0]) + (r32[2] * tau2[1]) + (r33[2] * tau2[2]);
  
  // Rotate spatial metric tensor derivative to local coordinate frame.
  qlocal[40] = (s11[0] * norm[0]) + (s21[0] * norm[1]) + (s31[0] * norm[2]);
  qlocal[41] = (s11[1] * norm[0]) + (s21[1] * norm[1]) + (s31[1] * norm[2]);
  qlocal[42] = (s11[2] * norm[0]) + (s21[2] * norm[1]) + (s31[2] * norm[2]);

  qlocal[43] = (s12[0] * norm[0]) + (s22[0] * norm[1]) + (s32[0] * norm[2]);
  qlocal[44] = (s12[1] * norm[0]) + (s22[1] * norm[1]) + (s32[1] * norm[2]);
  qlocal[45] = (s12[2] * norm[0]) + (s22[2] * norm[1]) + (s32[2] * norm[2]);

  qlocal[46] = (s13[0] * norm[0]) + (s23[0] * norm[1]) + (s33[0] * norm[2]);
  qlocal[47] = (s13[1] * norm[0]) + (s23[1] * norm[1]) + (s33[1] * norm[2]);
  qlocal[48] = (s13[2] * norm[0]) + (s23[2] * norm[1]) + (s33[2] * norm[2]);

  qlocal[49] = (s11[0] * tau1[0]) + (s21[0] * tau1[1]) + (s31[0] * tau1[2]);
  qlocal[50] = (s11[1] * tau1[0]) + (s21[1] * tau1[1]) + (s31[1] * tau1[2]);
  qlocal[51] = (s11[2] * tau1[0]) + (s21[2] * tau1[1]) + (s31[2] * tau1[2]);

  qlocal[52] = (s12[0] * tau1[0]) + (s22[0] * tau1[1]) + (s32[0] * tau1[2]);
  qlocal[53] = (s12[1] * tau1[0]) + (s22[1] * tau1[1]) + (s32[1] * tau1[2]);
  qlocal[54] = (s12[2] * tau1[0]) + (s22[2] * tau1[1]) + (s32[2] * tau1[2]);

  qlocal[55] = (s13[0] * tau1[0]) + (s23[0] * tau1[1]) + (s33[0] * tau1[2]);
  qlocal[56] = (s13[1] * tau1[0]) + (s23[1] * tau1[1]) + (s33[1] * tau1[2]);
  qlocal[57] = (s13[2] * tau1[0]) + (s23[2] * tau1[1]) + (s33[2] * tau1[2]);

  qlocal[58] = (s11[0] * tau2[0]) + (s21[0] * tau2[1]) + (s31[0] * tau2[2]);
  qlocal[59] = (s11[1] * tau2[0]) + (s21[1] * tau2[1]) + (s31[1] * tau2[2]);
  qlocal[60] = (s11[2] * tau2[0]) + (s21[2] * tau2[1]) + (s31[2] * tau2[2]);

  qlocal[61] = (s12[0] * tau2[0]) + (s22[0] * tau2[1]) + (s32[0] * tau2[2]);
  qlocal[62] = (s12[1] * tau2[0]) + (s22[1] * tau2[1]) + (s32[1] * tau2[2]);
  qlocal[63] = (s12[2] * tau2[0]) + (s22[2] * tau2[1]) + (s32[2] * tau2[2]);

  qlocal[64] = (s13[0] * tau2[0]) + (s23[0] * tau2[1]) + (s33[0] * tau2[2]);
  qlocal[65] = (s13[1] * tau2[0]) + (s23[1] * tau2[1]) + (s33[1] * tau2[2]);
  qlocal[66] = (s13[2] * tau2[0]) + (s23[2] * tau2[1]) + (s33[2] * tau2[2]);

  qlocal[67] = qglobal[67];
  qlocal[68] = (qglobal[68] * norm[0]) + (qglobal[69] * norm[1]) + (qglobal[70] * norm[2]);
  qlocal[69] = (qglobal[68] * tau1[0]) + (qglobal[69] * tau1[1]) + (qglobal[70] * tau1[2]);
  qlocal[70] = (qglobal[68] * tau2[0]) + (qglobal[69] * tau2[1]) + (qglobal[70] * tau2[2]);

  // Static-TOV well-balancing carries the frozen equilibrium: q[71]=D_eq, q[72]=Etot_eq are scalars
  // (carried unrotated); q[73..75]=S_eq,i is the covariant equilibrium momentum vector, therefore we need rotation
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
  if (gr_euler->tov_eq != NULL) {
    qlocal[71] = qglobal[71];
    qlocal[72] = qglobal[72];
    qlocal[73] = (qglobal[73] * norm[0]) + (qglobal[74] * norm[1]) + (qglobal[75] * norm[2]);
    qlocal[74] = (qglobal[73] * tau1[0]) + (qglobal[74] * tau1[1]) + (qglobal[75] * tau1[2]);
    qlocal[75] = (qglobal[73] * tau2[0]) + (qglobal[74] * tau2[1]) + (qglobal[75] * tau2[2]);
  }
}

static inline void
rot_to_global(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qlocal,
  double* GKYL_RESTRICT qglobal)
{
  qglobal[0] = qlocal[0];
  qglobal[1] = (qlocal[1] * norm[0]) + (qlocal[2] * tau1[0]) + (qlocal[3] * tau2[0]);
  qglobal[2] = (qlocal[1] * norm[1]) + (qlocal[2] * tau1[1]) + (qlocal[3] * tau2[1]);
  qglobal[3] = (qlocal[1] * norm[2]) + (qlocal[2] * tau1[2]) + (qlocal[3] * tau2[2]);
  qglobal[4] = qlocal[4];

  qglobal[5] = qlocal[5];
  qglobal[6] = (qlocal[6] * norm[0]) + (qlocal[7] * tau1[0]) + (qlocal[8] * tau2[0]);
  qglobal[7] = (qlocal[6] * norm[1]) + (qlocal[7] * tau1[1]) + (qlocal[8] * tau2[1]);
  qglobal[8] = (qlocal[6] * norm[2]) + (qlocal[7] * tau1[2]) + (qlocal[8] * tau2[2]);

  // Temporary arrays to store rotated column vectors.
  double r1[3], r2[3], r3[3];
  r1[0] = (qlocal[9] * norm[0]) + (qlocal[10] * tau1[0]) + (qlocal[11] * tau2[0]);
  r1[1] = (qlocal[9] * norm[1]) + (qlocal[10] * tau1[1]) + (qlocal[11] * tau2[1]);
  r1[2] = (qlocal[9] * norm[2]) + (qlocal[10] * tau1[2]) + (qlocal[11] * tau2[2]);

  r2[0] = (qlocal[12] * norm[0]) + (qlocal[13] * tau1[0]) + (qlocal[14] * tau2[0]);
  r2[1] = (qlocal[12] * norm[1]) + (qlocal[13] * tau1[1]) + (qlocal[14] * tau2[1]);
  r2[2] = (qlocal[12] * norm[2]) + (qlocal[13] * tau1[2]) + (qlocal[14] * tau2[2]);

  r3[0] = (qlocal[15] * norm[0]) + (qlocal[16] * tau1[0]) + (qlocal[17] * tau2[0]);
  r3[1] = (qlocal[15] * norm[1]) + (qlocal[16] * tau1[1]) + (qlocal[17] * tau2[1]);
  r3[2] = (qlocal[15] * norm[2]) + (qlocal[16] * tau1[2]) + (qlocal[17] * tau2[2]);

  // Temporary arrays to store rotated row vectors.
  double v1[3], v2[3], v3[3];
  v1[0] = (r1[0] * norm[0]) + (r2[0] * tau1[0]) + (r3[0] * tau2[0]);
  v1[1] = (r1[0] * norm[1]) + (r2[0] * tau1[1]) + (r3[0] * tau2[1]);
  v1[2] = (r1[0] * norm[2]) + (r2[0] * tau1[2]) + (r3[0] * tau2[2]);

  v2[0] = (r1[1] * norm[0]) + (r2[1] * tau1[0]) + (r3[1] * tau2[0]);
  v2[1] = (r1[1] * norm[1]) + (r2[1] * tau1[1]) + (r3[1] * tau2[1]);
  v2[2] = (r1[1] * norm[2]) + (r2[1] * tau1[2]) + (r3[1] * tau2[2]);

  v3[0] = (r1[2] * norm[0]) + (r2[2] * tau1[0]) + (r3[2] * tau2[0]);
  v3[1] = (r1[2] * norm[1]) + (r2[2] * tau1[1]) + (r3[2] * tau2[1]);
  v3[2] = (r1[2] * norm[2]) + (r2[2] * tau1[2]) + (r3[2] * tau2[2]);

  // Rotate spatial metric tensor back to global coordinate frame.
  qglobal[9] = v1[0]; qglobal[10] = v1[1]; qglobal[11] = v1[2];
  qglobal[12] = v2[0]; qglobal[13] = v2[1]; qglobal[14] = v2[2];
  qglobal[15] = v3[0]; qglobal[16] = v3[1]; qglobal[17] = v3[2];

  // Temporary arrays to store rotated extrinsic column vectors.
  double extr_r1[3], extr_r2[3], extr_r3[3];
  extr_r1[0] = (qlocal[18] * norm[0]) + (qlocal[19] * tau1[0]) + (qlocal[20] * tau2[0]);
  extr_r1[1] = (qlocal[18] * norm[1]) + (qlocal[19] * tau1[1]) + (qlocal[20] * tau2[1]);
  extr_r1[2] = (qlocal[18] * norm[2]) + (qlocal[19] * tau1[2]) + (qlocal[20] * tau2[2]);

  extr_r2[0] = (qlocal[21] * norm[0]) + (qlocal[22] * tau1[0]) + (qlocal[23] * tau2[0]);
  extr_r2[1] = (qlocal[21] * norm[1]) + (qlocal[22] * tau1[1]) + (qlocal[23] * tau2[1]);
  extr_r2[2] = (qlocal[21] * norm[2]) + (qlocal[22] * tau1[2]) + (qlocal[23] * tau2[2]);

  extr_r3[0] = (qlocal[24] * norm[0]) + (qlocal[25] * tau1[0]) + (qlocal[26] * tau2[0]);
  extr_r3[1] = (qlocal[24] * norm[1]) + (qlocal[25] * tau1[1]) + (qlocal[26] * tau2[1]);
  extr_r3[2] = (qlocal[24] * norm[2]) + (qlocal[25] * tau1[2]) + (qlocal[26] * tau2[2]);

  // Temporary arrays to store rotated extrinsic row vectors.
  double inv_v1[3], inv_v2[3], inv_v3[3];
  inv_v1[0] = (extr_r1[0] * norm[0]) + (extr_r2[0] * tau1[0]) + (extr_r3[0] * tau2[0]);
  inv_v1[1] = (extr_r1[0] * norm[1]) + (extr_r2[0] * tau1[1]) + (extr_r3[0] * tau2[1]);
  inv_v1[2] = (extr_r1[0] * norm[2]) + (extr_r2[0] * tau1[2]) + (extr_r3[0] * tau2[2]);

  inv_v2[0] = (extr_r1[1] * norm[0]) + (extr_r2[1] * tau1[0]) + (extr_r3[1] * tau2[0]);
  inv_v2[1] = (extr_r1[1] * norm[1]) + (extr_r2[1] * tau1[1]) + (extr_r3[1] * tau2[1]);
  inv_v2[2] = (extr_r1[1] * norm[2]) + (extr_r2[1] * tau1[2]) + (extr_r3[1] * tau2[2]);

  inv_v3[0] = (extr_r1[2] * norm[0]) + (extr_r2[2] * tau1[0]) + (extr_r3[2] * tau2[0]);
  inv_v3[1] = (extr_r1[2] * norm[1]) + (extr_r2[2] * tau1[1]) + (extr_r3[2] * tau2[1]);
  inv_v3[2] = (extr_r1[2] * norm[2]) + (extr_r2[2] * tau1[2]) + (extr_r3[2] * tau2[2]);

  // Rotate extrinsic curvature tensor back to global coordinate frame.
  qglobal[18] = inv_v1[0]; qglobal[19] = inv_v1[1]; qglobal[20] = inv_v1[2];
  qglobal[21] = inv_v2[0]; qglobal[22] = inv_v2[1]; qglobal[23] = inv_v2[2];
  qglobal[24] = inv_v3[0]; qglobal[25] = inv_v3[1]; qglobal[26] = inv_v3[2];

  qglobal[27] = qlocal[27];

  qglobal[28] = (qlocal[28] * norm[0]) + (qlocal[29] * tau1[0]) + (qlocal[30] * tau2[0]);
  qglobal[29] = (qlocal[28] * norm[1]) + (qlocal[29] * tau1[1]) + (qlocal[30] * tau2[1]);
  qglobal[30] = (qlocal[28] * norm[2]) + (qlocal[29] * tau1[2]) + (qlocal[30] * tau2[2]);

  // Temporary arrays to store rotated shift derivative column vectors.
  double shiftder_r1[3], shiftder_r2[3], shiftder_r3[3];
  shiftder_r1[0] = (qlocal[31] * norm[0]) + (qlocal[32] * tau1[0]) + (qlocal[33] * tau2[0]);
  shiftder_r1[1] = (qlocal[31] * norm[1]) + (qlocal[32] * tau1[1]) + (qlocal[33] * tau2[1]);
  shiftder_r1[2] = (qlocal[31] * norm[2]) + (qlocal[32] * tau1[2]) + (qlocal[33] * tau2[2]);

  shiftder_r2[0] = (qlocal[34] * norm[0]) + (qlocal[35] * tau1[0]) + (qlocal[36] * tau2[0]);
  shiftder_r2[1] = (qlocal[34] * norm[1]) + (qlocal[35] * tau1[1]) + (qlocal[36] * tau2[1]);
  shiftder_r2[2] = (qlocal[34] * norm[2]) + (qlocal[35] * tau1[2]) + (qlocal[36] * tau2[2]);

  shiftder_r3[0] = (qlocal[37] * norm[0]) + (qlocal[38] * tau1[0]) + (qlocal[39] * tau2[0]);
  shiftder_r3[1] = (qlocal[37] * norm[1]) + (qlocal[38] * tau1[1]) + (qlocal[39] * tau2[1]);
  shiftder_r3[2] = (qlocal[37] * norm[2]) + (qlocal[38] * tau1[2]) + (qlocal[39] * tau2[2]);

  // Temporary arrays to store rotated shift derivative row vectors.
  double shiftder_v1[3], shiftder_v2[3], shiftder_v3[3];
  shiftder_v1[0] = (shiftder_r1[0] * norm[0]) + (shiftder_r2[0] * tau1[0]) + (shiftder_r3[0] * tau2[0]);
  shiftder_v1[1] = (shiftder_r1[0] * norm[1]) + (shiftder_r2[0] * tau1[1]) + (shiftder_r3[0] * tau2[1]);
  shiftder_v1[2] = (shiftder_r1[0] * norm[2]) + (shiftder_r2[0] * tau1[2]) + (shiftder_r3[0] * tau2[2]);

  shiftder_v2[0] = (shiftder_r1[1] * norm[0]) + (shiftder_r2[1] * tau1[0]) + (shiftder_r3[1] * tau2[0]);
  shiftder_v2[1] = (shiftder_r1[1] * norm[1]) + (shiftder_r2[1] * tau1[1]) + (shiftder_r3[1] * tau2[1]);
  shiftder_v2[2] = (shiftder_r1[1] * norm[2]) + (shiftder_r2[1] * tau1[2]) + (shiftder_r3[1] * tau2[2]);

  shiftder_v3[0] = (shiftder_r1[2] * norm[0]) + (shiftder_r2[2] * tau1[0]) + (shiftder_r3[2] * tau2[0]);
  shiftder_v3[1] = (shiftder_r1[2] * norm[1]) + (shiftder_r2[2] * tau1[1]) + (shiftder_r3[2] * tau2[1]);
  shiftder_v3[2] = (shiftder_r1[2] * norm[2]) + (shiftder_r2[2] * tau1[2]) + (shiftder_r3[2] * tau2[2]);

  // Rotate shift vector derivative back to global coordinate frame.
  qglobal[31] = shiftder_v1[0]; qglobal[32] = shiftder_v1[1]; qglobal[33] = shiftder_v1[2];
  qglobal[34] = shiftder_v2[0]; qglobal[35] = shiftder_v2[1]; qglobal[36] = shiftder_v2[2];
  qglobal[37] = shiftder_v3[0]; qglobal[38] = shiftder_v3[1]; qglobal[39] = shiftder_v3[2];

  // Temporary arrays to store rotated column vectors.
  double r11[3], r12[3], r13[3];
  double r21[3], r22[3], r23[3];
  double r31[3], r32[3], r33[3];

  r11[0] = (qlocal[40] * norm[0]) + (qlocal[41] * tau1[0]) + (qlocal[42] * tau2[0]);
  r11[1] = (qlocal[40] * norm[1]) + (qlocal[41] * tau1[1]) + (qlocal[42] * tau2[1]);
  r11[2] = (qlocal[40] * norm[2]) + (qlocal[41] * tau1[2]) + (qlocal[42] * tau2[2]);

  r12[0] = (qlocal[43] * norm[0]) + (qlocal[44] * tau1[0]) + (qlocal[45] * tau2[0]);
  r12[1] = (qlocal[43] * norm[1]) + (qlocal[44] * tau1[1]) + (qlocal[45] * tau2[1]);
  r12[2] = (qlocal[43] * norm[2]) + (qlocal[44] * tau1[2]) + (qlocal[45] * tau2[2]);

  r13[0] = (qlocal[46] * norm[0]) + (qlocal[47] * tau1[0]) + (qlocal[48] * tau2[0]);
  r13[1] = (qlocal[46] * norm[1]) + (qlocal[47] * tau1[1]) + (qlocal[48] * tau2[1]);
  r13[2] = (qlocal[46] * norm[2]) + (qlocal[47] * tau1[2]) + (qlocal[48] * tau2[2]);

  r21[0] = (qlocal[49] * norm[0]) + (qlocal[50] * tau1[0]) + (qlocal[51] * tau2[0]);
  r21[1] = (qlocal[49] * norm[1]) + (qlocal[50] * tau1[1]) + (qlocal[51] * tau2[1]);
  r21[2] = (qlocal[49] * norm[2]) + (qlocal[50] * tau1[2]) + (qlocal[51] * tau2[2]);

  r22[0] = (qlocal[52] * norm[0]) + (qlocal[53] * tau1[0]) + (qlocal[54] * tau2[0]);
  r22[1] = (qlocal[52] * norm[1]) + (qlocal[53] * tau1[1]) + (qlocal[54] * tau2[1]);
  r22[2] = (qlocal[52] * norm[2]) + (qlocal[53] * tau1[2]) + (qlocal[54] * tau2[2]);

  r23[0] = (qlocal[55] * norm[0]) + (qlocal[56] * tau1[0]) + (qlocal[57] * tau2[0]);
  r23[1] = (qlocal[55] * norm[1]) + (qlocal[56] * tau1[1]) + (qlocal[57] * tau2[1]);
  r23[2] = (qlocal[55] * norm[2]) + (qlocal[56] * tau1[2]) + (qlocal[57] * tau2[2]);

  r31[0] = (qlocal[58] * norm[0]) + (qlocal[59] * tau1[0]) + (qlocal[60] * tau2[0]);
  r31[1] = (qlocal[58] * norm[1]) + (qlocal[59] * tau1[1]) + (qlocal[60] * tau2[1]);
  r31[2] = (qlocal[58] * norm[2]) + (qlocal[59] * tau1[2]) + (qlocal[60] * tau2[2]);

  r32[0] = (qlocal[61] * norm[0]) + (qlocal[62] * tau1[0]) + (qlocal[63] * tau2[0]);
  r32[1] = (qlocal[61] * norm[1]) + (qlocal[62] * tau1[1]) + (qlocal[63] * tau2[1]);
  r32[2] = (qlocal[61] * norm[2]) + (qlocal[62] * tau1[2]) + (qlocal[63] * tau2[2]);

  r33[0] = (qlocal[64] * norm[0]) + (qlocal[65] * tau1[0]) + (qlocal[66] * tau2[0]);
  r33[1] = (qlocal[64] * norm[1]) + (qlocal[65] * tau1[1]) + (qlocal[66] * tau2[1]);
  r33[2] = (qlocal[64] * norm[2]) + (qlocal[65] * tau1[2]) + (qlocal[66] * tau2[2]);

  // Temporary arrays to store rotated row vectors.
  double s11[3], s12[3], s13[3];
  double s21[3], s22[3], s23[3];
  double s31[3], s32[3], s33[3];

  s11[0] = (r11[0] * norm[0]) + (r12[0] * tau1[0]) + (r13[0] * tau2[0]);
  s11[1] = (r11[1] * norm[0]) + (r12[1] * tau1[0]) + (r13[1] * tau2[0]);
  s11[2] = (r11[2] * norm[0]) + (r12[2] * tau1[0]) + (r13[2] * tau2[0]);

  s12[0] = (r11[0] * norm[1]) + (r12[0] * tau1[1]) + (r13[0] * tau2[1]);
  s12[1] = (r11[1] * norm[1]) + (r12[1] * tau1[1]) + (r13[1] * tau2[1]);
  s12[2] = (r11[2] * norm[1]) + (r12[2] * tau1[1]) + (r13[2] * tau2[1]);

  s13[0] = (r11[0] * norm[2]) + (r12[0] * tau1[2]) + (r13[0] * tau2[2]);
  s13[1] = (r11[1] * norm[2]) + (r12[1] * tau1[2]) + (r13[1] * tau2[2]);
  s13[2] = (r11[2] * norm[2]) + (r12[2] * tau1[2]) + (r13[2] * tau2[2]);

  s21[0] = (r21[0] * norm[0]) + (r22[0] * tau1[0]) + (r23[0] * tau2[0]);
  s21[1] = (r21[1] * norm[0]) + (r22[1] * tau1[0]) + (r23[1] * tau2[0]);
  s21[2] = (r21[2] * norm[0]) + (r22[2] * tau1[0]) + (r23[2] * tau2[0]);

  s22[0] = (r21[0] * norm[1]) + (r22[0] * tau1[1]) + (r23[0] * tau2[1]);
  s22[1] = (r21[1] * norm[1]) + (r22[1] * tau1[1]) + (r23[1] * tau2[1]);
  s22[2] = (r21[2] * norm[1]) + (r22[2] * tau1[1]) + (r23[2] * tau2[1]);

  s23[0] = (r21[0] * norm[2]) + (r22[0] * tau1[2]) + (r23[0] * tau2[2]);
  s23[1] = (r21[1] * norm[2]) + (r22[1] * tau1[2]) + (r23[1] * tau2[2]);
  s23[2] = (r21[2] * norm[2]) + (r22[2] * tau1[2]) + (r23[2] * tau2[2]);

  s31[0] = (r31[0] * norm[0]) + (r32[0] * tau1[0]) + (r33[0] * tau2[0]);
  s31[1] = (r31[1] * norm[0]) + (r32[1] * tau1[0]) + (r33[1] * tau2[0]);
  s31[2] = (r31[2] * norm[0]) + (r32[2] * tau1[0]) + (r33[2] * tau2[0]);

  s32[0] = (r31[0] * norm[1]) + (r32[0] * tau1[1]) + (r33[0] * tau2[1]);
  s32[1] = (r31[1] * norm[1]) + (r32[1] * tau1[1]) + (r33[1] * tau2[1]);
  s32[2] = (r31[2] * norm[1]) + (r32[2] * tau1[1]) + (r33[2] * tau2[1]);

  s33[0] = (r31[0] * norm[2]) + (r32[0] * tau1[2]) + (r33[0] * tau2[2]);
  s33[1] = (r31[1] * norm[2]) + (r32[1] * tau1[2]) + (r33[1] * tau2[2]);
  s33[2] = (r31[2] * norm[2]) + (r32[2] * tau1[2]) + (r33[2] * tau2[2]);

  // Rotate spatial metric tensor derivative back to global coordinate frame.
  qglobal[40] = (s11[0] * norm[0]) + (s21[0] * tau1[0]) + (s31[0] * tau2[0]);
  qglobal[41] = (s11[1] * norm[0]) + (s21[1] * tau1[0]) + (s31[1] * tau2[0]);
  qglobal[42] = (s11[2] * norm[0]) + (s21[2] * tau1[0]) + (s31[2] * tau2[0]);

  qglobal[43] = (s12[0] * norm[0]) + (s22[0] * tau1[0]) + (s32[0] * tau2[0]);
  qglobal[44] = (s12[1] * norm[0]) + (s22[1] * tau1[0]) + (s32[1] * tau2[0]);
  qglobal[45] = (s12[2] * norm[0]) + (s22[2] * tau1[0]) + (s32[2] * tau2[0]);

  qglobal[46] = (s13[0] * norm[0]) + (s23[0] * tau1[0]) + (s33[0] * tau2[0]);
  qglobal[47] = (s13[1] * norm[0]) + (s23[1] * tau1[0]) + (s33[1] * tau2[0]);
  qglobal[48] = (s13[2] * norm[0]) + (s23[2] * tau1[0]) + (s33[2] * tau2[0]);

  qglobal[49] = (s11[0] * norm[1]) + (s21[0] * tau1[1]) + (s31[0] * tau2[1]);
  qglobal[50] = (s11[1] * norm[1]) + (s21[1] * tau1[1]) + (s31[1] * tau2[1]);
  qglobal[51] = (s11[2] * norm[1]) + (s21[2] * tau1[1]) + (s31[2] * tau2[1]);

  qglobal[52] = (s12[0] * norm[1]) + (s22[0] * tau1[1]) + (s32[0] * tau2[1]);
  qglobal[53] = (s12[1] * norm[1]) + (s22[1] * tau1[1]) + (s32[1] * tau2[1]);
  qglobal[54] = (s12[2] * norm[1]) + (s22[2] * tau1[1]) + (s32[2] * tau2[1]);

  qglobal[55] = (s13[0] * norm[1]) + (s23[0] * tau1[1]) + (s33[0] * tau2[1]);
  qglobal[56] = (s13[1] * norm[1]) + (s23[1] * tau1[1]) + (s33[1] * tau2[1]);
  qglobal[57] = (s13[2] * norm[1]) + (s23[2] * tau1[1]) + (s33[2] * tau2[1]);

  qglobal[58] = (s11[0] * norm[2]) + (s21[0] * tau1[2]) + (s31[0] * tau2[2]);
  qglobal[59] = (s11[1] * norm[2]) + (s21[1] * tau1[2]) + (s31[1] * tau2[2]);
  qglobal[60] = (s11[2] * norm[2]) + (s21[2] * tau1[2]) + (s31[2] * tau2[2]);

  qglobal[61] = (s12[0] * norm[2]) + (s22[0] * tau1[2]) + (s32[0] * tau2[2]);
  qglobal[62] = (s12[1] * norm[2]) + (s22[1] * tau1[2]) + (s32[1] * tau2[2]);
  qglobal[63] = (s12[2] * norm[2]) + (s22[2] * tau1[2]) + (s32[2] * tau2[2]);

  qglobal[64] = (s13[0] * norm[2]) + (s23[0] * tau1[2]) + (s33[0] * tau2[2]);
  qglobal[65] = (s13[1] * norm[2]) + (s23[1] * tau1[2]) + (s33[1] * tau2[2]);
  qglobal[66] = (s13[2] * norm[2]) + (s23[2] * tau1[2]) + (s33[2] * tau2[2]);

  qglobal[67] = qlocal[67];
  qglobal[68] = (qlocal[68] * norm[0]) + (qlocal[69] * tau1[0]) + (qlocal[70] * tau2[0]);
  qglobal[69] = (qlocal[68] * norm[1]) + (qlocal[69] * tau1[1]) + (qlocal[70] * tau2[1]);
  qglobal[70] = (qlocal[68] * norm[2]) + (qlocal[69] * tau1[2]) + (qlocal[70] * tau2[2]);

  // Frozen-equilibrium WB slots: q[71]=D_eq, q[72]=Etot_eq are scalars (carried); q[73..75]=S_eq,i is the covariant equilibrium momentum vector, rotated back like the live momentum
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
  if (gr_euler->tov_eq != NULL) {
    qglobal[71] = qlocal[71];
    qglobal[72] = qlocal[72];
    qglobal[73] = (qlocal[73] * norm[0]) + (qlocal[74] * tau1[0]) + (qlocal[75] * tau2[0]);
    qglobal[74] = (qlocal[73] * norm[1]) + (qlocal[74] * tau1[1]) + (qlocal[75] * tau2[1]);
    qglobal[75] = (qlocal[73] * norm[2]) + (qlocal[74] * tau1[2]) + (qlocal[75] * tau2[2]);
  }
}

static double
wave_lax(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
  double gas_gamma = gr_euler->gas_gamma;

  double sl = gkyl_gr_euler_max_abs_speed(gas_gamma, ql);
  double sr = gkyl_gr_euler_max_abs_speed(gas_gamma, qr);
  double amax = fmax(sl, sr);

  double fl[71], fr[71];
  gkyl_gr_euler_flux(gas_gamma, ql, fl);
  gkyl_gr_euler_flux(gas_gamma, qr, fr);

  bool in_excision_region_l = false;
  if (ql[27] < pow(10.0, -8.0)) {
    in_excision_region_l = true;
  }

  bool in_excision_region_r = false;
  if (qr[27] < pow(10.0, -8.0)) {
    in_excision_region_r = true;
  }

  // Optional static-TOV well-balancing: feed the wave construction the equilibrium deviation dq = (qr-ql) - deq and df = (fr-fl) - dfeq, where (deq, dfeq) are the frozen TOV equilibrium jumps in the conserved variables and their fluxes.
  // amdq + apdq stays (fr-fl) - dfeq, so mass/energy (dfeq = 0 at v = 0) remain strictly conserved and the momentum keeps its (gravity-sourced) balance; a static star produces ~0 net flux.
  double deq[71] = { 0.0 }, dfeq[71] = { 0.0 };
  if (gr_euler->tov_eq && !gr_euler->disable_well_balanced && !in_excision_region_l && !in_excision_region_r) {
    double q_eq_l[71], q_eq_r[71], f_eq_l[71], f_eq_r[71];
    gr_euler_tov_equilibrium(gr_euler, ql, q_eq_l);
    gr_euler_tov_equilibrium(gr_euler, qr, q_eq_r);
    gkyl_gr_euler_flux(gas_gamma, q_eq_l, f_eq_l);
    gkyl_gr_euler_flux(gas_gamma, q_eq_r, f_eq_r);
    for (int i = 0; i < 71; i++) {
      deq[i] = q_eq_r[i] - q_eq_l[i];
      dfeq[i] = f_eq_r[i] - f_eq_l[i];
    }
  }

  int neq = eqn->num_equations;
  double *w0 = &waves[0], *w1 = &waves[neq];
  if (!in_excision_region_l && !in_excision_region_r) {
    for (int i = 0; i < 71; i++) {
      double dq = (qr[i] - ql[i]) - deq[i];
      double df = (fr[i] - fl[i]) - dfeq[i];
      w0[i] = 0.5 * (dq - (df / amax));
      w1[i] = 0.5 * (dq + (df / amax));
    }
  }
  else {
    for (int i = 0; i < 71; i++) {
      w0[i] = 0.0;
      w1[i] = 0.0;
    }
  }
  // The frozen-equilibrium slots q[71],q[72] carry zero waves, so they never evolve
  if (gr_euler->tov_eq != NULL) {
    for (int i = 71; i < neq; i++) {
      w0[i] = 0.0;
      w1[i] = 0.0;
    }
  }

  s[0] = -amax;
  s[1] = amax;

  return s[1];
}

static void
qfluct_lax(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
{
  int neq = eqn->num_equations;
  const double *w0 = &waves[0], *w1 = &waves[neq];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i = 0; i < neq; i++) { // frozen slots (i>=71) have zero waves -> amdq=apdq=0 (frozen)
    amdq[i] = (s0m * w0[i]) + (s1m * w1[i]);
    apdq[i] = (s0p * w0[i]) + (s1p * w1[i]);
  }
}

static double
wave_lax_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* delta, const double* ql, const double* qr, const double phil, const double phir, double* waves, double* s)
{
  return wave_lax(eqn, delta, ql, qr, waves, s);
}

static void
qfluct_lax_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* ql, const double* qr, const double phil, const double phir, const double* waves, const double* s,
  double* amdq, double* apdq)
{
  return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
}

static double
wave_roe(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
  double vl[71], vr[71];
  double gas_gamma = gr_euler->gas_gamma;
  
  gkyl_gr_euler_prim_vars(gas_gamma, ql, vl);
  gkyl_gr_euler_prim_vars(gas_gamma, qr, vr);

  double rho_l = vl[0];
  double vx_l = vl[1];
  double vy_l = vl[2];
  double vz_l = vl[3];
  double p_l = vl[4];

  double rho_r = vr[0];
  double vx_r = vr[1];
  double vy_r = vr[2];
  double vz_r = vr[3];
  double p_r = vr[4];

  double Etot_l = ql[4];
  double Etot_r = qr[4];

  double W_l = 1.0 / sqrt(1.0 - ((vx_l * vx_l) + (vy_l * vy_l) + (vz_l * vz_l)));
  double W_r = 1.0 / sqrt(1.0 - ((vx_r * vx_r) + (vy_r * vy_r) + (vz_r * vz_r)));

  double K_l = sqrt(Etot_l + p_l) / W_l;
  double K_r = sqrt(Etot_r + p_r) / W_r;
  double K_avg = 1.0 / (K_l + K_r);

  double v0 = ((K_l * W_l) + (K_r * W_r)) * K_avg;
  double v1 = ((K_l * W_l * vx_l) + (K_r * W_r * vx_r)) * K_avg;
  double v2 = ((K_l * W_l * vy_l) + (K_r * W_r * vy_r)) * K_avg;
  double v3 = ((K_l * W_l * vz_l) + (K_r * W_r * vz_r)) * K_avg;
  double v4 = ((p_l / K_l) + (p_r / K_r)) * K_avg;

  double c_minus = 1.0 - ((gas_gamma / (gas_gamma - 1.0)) * v4);
  double c_plus = 1.0 + ((gas_gamma / (gas_gamma - 1.0)) * v4);

  double v_alpha_sq = -(v0 * v0) + (v1 * v1) + (v2 * v2) + (v3 * v3);
  double s_sq = (0.5 * gas_gamma * v4 * (1.0 - v_alpha_sq)) - (0.5 * (gas_gamma - 1.0) * (1.0 + v_alpha_sq));
  double energy = (v0 * v0) - (v1 * v1);
  double y = sqrt(((1.0 - (gas_gamma * v4)) * energy) + s_sq);

  double k = (v0 * delta[4]) - (v1 * delta[1]);
  double v_delta = (-v0 * delta[4]) + (v1 * delta[1]) + (v2 * delta[2]) + (v3 * delta[3]);
  double a1 = -((s_sq * k) + (sqrt(s_sq) * y * ((v0 * delta[1]) - (v1 * delta[4])) + ((gas_gamma - 1.0) * energy * (delta[0] + (c_plus * v_delta))))) / (2.0 * energy * s_sq);
  double a2 = -((s_sq * k) - (sqrt(s_sq) * y * ((v0 * delta[1]) - (v1 * delta[4])) + ((gas_gamma - 1.0) * energy * (delta[0] + (c_plus * v_delta))))) / (2.0 * energy * s_sq);
  double a3 = ((2.0 * s_sq * k) + ((gas_gamma - 1.0) * energy * (delta[0] + (c_plus * v_delta)))) / (energy * s_sq);
  double a4 = delta[2] - ((k * v2) / energy);
  double a5 = delta[3] - ((k * v3) / energy);

  int neq = eqn->num_equations;
  for (int i = 0; i < neq * 3; i++) {
    waves[i] = 0;
  }

  double *wv;
  wv = &waves[0 * neq];
  wv[0] = a1 * c_minus;
  wv[1] = a1 * (v1 - ((sqrt(s_sq) * v0) / y));
  wv[2] = a1 * v2;
  wv[3] = a1 * v3;
  wv[4] = a1 * (v0 - ((sqrt(s_sq) * v1) / y));
  s[0] = (((1.0 - (gas_gamma * v4)) * v0 * v1) - (sqrt(s_sq) * y)) / (((1.0 - (gas_gamma * v4)) * v0 * v0) + s_sq);

  wv = &waves[1 * neq];
  wv[0] = (a3 * (c_minus + (s_sq / (gas_gamma - 1.0)))) - (a4 * c_plus * v2) - (a5 * c_plus * v3);
  wv[1] = a3 * v1;
  wv[2] = (a3 * v2) + a4;
  wv[3] = (a3 * v3) + a5;
  wv[4] = a3 * v0;
  s[1] = v1 / v0;

  wv = &waves[2 * neq];
  wv[0] = a2 * c_minus;
  wv[1] = a2 * (v1 + ((sqrt(s_sq) * v0) / y));
  wv[2] = a2 * v2;
  wv[3] = a2 * v3;
  wv[4] = a2 * (v0 + ((sqrt(s_sq) * v1) / y));
  s[2] = (((1.0 - (gas_gamma * v4)) * v0 * v1) + (sqrt(s_sq) * y)) / (((1.0 - (gas_gamma * v4)) * v0 * v0) + s_sq);

  return (((1.0 - (gas_gamma * v4)) * v0 * fabs(v1)) + (sqrt(s_sq) * y)) / (((1.0 - (gas_gamma * v4)) * v0 * v0) + s_sq);
}

static void
qfluct_roe(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
{
  int neq = eqn->num_equations;
  const double* w0 = &waves[0 * neq];
  const double* w1 = &waves[1 * neq];
  const double* w2 = &waves[2 * neq];

  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]), s2m = fmin(0.0, s[2]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]), s2p = fmax(0.0, s[2]);

  for (int i = 0; i < 5; i++) {
    amdq[i] = (s0m * w0[i]) + (s1m * w1[i]) + (s2m * w2[i]);
    apdq[i] = (s0p * w0[i]) + (s1p * w1[i]) + (s2p * w2[i]);
  }
  for (int i = 5; i < neq; i++) {
    amdq[i] = 0.0;
    apdq[i] = 0.0;
  }
}

static double
wave_roe_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* delta, const double* ql, const double* qr, const double phil, const double phir, double* waves, double* s)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX) {
    return wave_roe(eqn, delta, ql, qr, waves, s);
  }
  else {
    return wave_lax(eqn, delta, ql, qr, waves, s);
  }
}

static void
qfluct_roe_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* ql, const double* qr, const double phil, const double phir, const double* waves, const double* s,
  double* amdq, double* apdq)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX) {
    return qfluct_roe(eqn, ql, qr, waves, s, amdq, apdq);
  }
  else {
    return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
  }
}

static double
wave_hll(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
  double gas_gamma = gr_euler->gas_gamma;

  double vl[71], vr[71];
  gkyl_gr_euler_prim_vars(gas_gamma, ql, vl);
  gkyl_gr_euler_prim_vars(gas_gamma, qr, vr);

  double rho_l = vl[0];
  double vx_l = vl[1];
  double vy_l = vl[2];
  double vz_l = vl[3];
  double p_l = vl[4];

  double lapse_l = vl[5];
  double shift_xl = vl[6];
  double shift_yl = vl[7];
  double shift_zl = vl[8];

  double spatial_metric_l[3][3];
  spatial_metric_l[0][0] = vl[9]; spatial_metric_l[0][1] = vl[10]; spatial_metric_l[0][2] = vl[11];
  spatial_metric_l[1][0] = vl[12]; spatial_metric_l[1][1] = vl[13]; spatial_metric_l[1][2] = vl[14];
  spatial_metric_l[2][0] = vl[15]; spatial_metric_l[2][1] = vl[16]; spatial_metric_l[2][2] = vl[17];

  double **inv_spatial_metric_l = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric_l[i] = gkyl_malloc(sizeof(double[3]));
  }

  gkyl_gr_euler_inv_spatial_metric(ql, &inv_spatial_metric_l);

  bool in_excision_region_l = false;
  if (vl[27] < pow(10.0, -8.0)) {
    in_excision_region_l = true;
  }

  bool curved_spacetime_l = false;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j) {
        if (fabs(spatial_metric_l[i][j] - 1.0) > pow(10.0, -8.0)) {
          curved_spacetime_l = true;
        }
      }
      else {
        if (fabs(spatial_metric_l[i][j]) > pow(10.0, -8.0)) {
          curved_spacetime_l = true;
        }
      }
    }
  }
  if (fabs(lapse_l - 1.0) > pow(10.0, -8.0) || fabs(shift_xl) > pow(10.0, -8.0) || fabs(shift_yl) > pow(10.0, -8.0) ||
    fabs(shift_zl) > pow(10.0, -8.0)) {
    curved_spacetime_l = true;
  }

  double rho_r = vr[0];
  double vx_r = vr[1];
  double vy_r = vr[2];
  double vz_r = vr[3];
  double p_r = vr[4];

  double lapse_r = vr[5];
  double shift_xr = vr[6];
  double shift_yr = vr[7];
  double shift_zr = vr[8];

  double spatial_metric_r[3][3];
  spatial_metric_r[0][0] = vr[9]; spatial_metric_r[0][1] = vr[10]; spatial_metric_r[0][2] = vr[11];
  spatial_metric_r[1][0] = vr[12]; spatial_metric_r[1][1] = vr[13]; spatial_metric_r[1][2] = vr[14];
  spatial_metric_r[2][0] = vr[15]; spatial_metric_r[2][1] = vr[16]; spatial_metric_r[2][2] = vr[17];

  double **inv_spatial_metric_r = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric_r[i] = gkyl_malloc(sizeof(double[3]));
  }

  gkyl_gr_euler_inv_spatial_metric(qr, &inv_spatial_metric_r);

  bool in_excision_region_r = false;
  if (vr[27] < pow(10.0, -8.0)) {
    in_excision_region_r = true;
  }

  bool curved_spacetime_r = false;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j) {
        if (fabs(spatial_metric_r[i][j] - 1.0) > pow(10.0, -8.0)) {
          curved_spacetime_r = true;
        }
      }
      else {
        if (fabs(spatial_metric_r[i][j]) > pow(10.0, -8.0)) {
          curved_spacetime_r = true;
        }
      }
    }
  }
  if (fabs(lapse_r - 1.0) > pow(10.0, -8.0) || fabs(shift_xr) > pow(10.0, -8.0) || fabs(shift_yr) > pow(10.0, -8.0) ||
    fabs(shift_zr) > pow(10.0, -8.0)) {
    curved_spacetime_r = true;
  }

  double c_sl = sqrt(((gas_gamma * p_l) / rho_l) / (1.0 + ((p_l / rho_l) * (gas_gamma) / (gas_gamma - 1.0))));
  double c_sr = sqrt(((gas_gamma * p_r) / rho_r) / (1.0 + ((p_r / rho_r) * (gas_gamma) / (gas_gamma - 1.0))));

  double vx_avg = 0.5 * (vx_l + vx_r);
  double cs_avg = 0.5 * (c_sl + c_sr);

  double sl, sr;

  if (curved_spacetime_l || curved_spacetime_r) {
    double vel_l[3];
    double v_sq_l = 0.0;
    vel_l[0] = vx_l; vel_l[1] = vy_l; vel_l[2] = vz_l;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        v_sq_l += spatial_metric_l[i][j] * vel_l[i] * vel_l[j];
      }
    }

    double shift_l[3];
    shift_l[0] = shift_xl; shift_l[1] = shift_yl; shift_l[2] = shift_zl;

    double material_eigs_l[3];
    double fast_acoustic_eigs_l[3];
    double slow_acoustic_eigs_l[3];

    for (int i = 0; i < 3; i++) {
      material_eigs_l[i] = (lapse_l * vel_l[i]) - shift_l[i];

      fast_acoustic_eigs_l[i] = (lapse_l / (1.0 - (v_sq_l * (c_sl * c_sl)))) * ((vel_l[i] * (1.0 - (c_sl * c_sl))) +
        (c_sl * sqrt((1.0 - v_sq_l) * (inv_spatial_metric_l[i][i] * (1.0 - (v_sq_l * (c_sl * c_sl))) - (vel_l[i] * vel_l[i]) * (1.0 - (c_sl * c_sl)))))) - shift_l[i];
      
      slow_acoustic_eigs_l[i] = (lapse_l / (1.0 - (v_sq_l * (c_sl * c_sl)))) * ((vel_l[i] * (1.0 - (c_sl * c_sl))) -
        (c_sl * sqrt((1.0 - v_sq_l) * (inv_spatial_metric_l[i][i] * (1.0 - (v_sq_l * (c_sl * c_sl))) - (vel_l[i] * vel_l[i]) * (1.0 - (c_sl * c_sl)))))) - shift_l[i];
    }

    double max_eig_l = 0.0;
    for (int i = 0; i < 3; i++) {
      if (fabs(material_eigs_l[i]) > max_eig_l) {
        max_eig_l = fabs(material_eigs_l[i]);
      }
      if (fabs(fast_acoustic_eigs_l[i]) > max_eig_l) {
        max_eig_l = fabs(fast_acoustic_eigs_l[i]);
      }
      if (fabs(slow_acoustic_eigs_l[i]) > max_eig_l) {
        max_eig_l = fabs(slow_acoustic_eigs_l[i]);
      }
    }

    double vel_r[3];
    double v_sq_r = 0.0;
    vel_r[0] = vx_r; vel_r[1] = vy_r; vel_r[2] = vz_r;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        v_sq_r += spatial_metric_r[i][j] * vel_r[i] * vel_r[j];
      }
    }

    double shift_r[3];
    shift_r[0] = shift_xr; shift_r[1] = shift_yr; shift_r[2] = shift_zr;

    double material_eigs_r[3];
    double fast_acoustic_eigs_r[3];
    double slow_acoustic_eigs_r[3];

    for (int i = 0; i < 3; i++) {
      material_eigs_r[i] = (lapse_r * vel_r[i]) - shift_r[i];

      // (bugfix, finding 3) right acoustic eigenvalues must use the RIGHT sound speed c_sr, not c_sl.
      fast_acoustic_eigs_r[i] = (lapse_r / (1.0 - (v_sq_r * (c_sr * c_sr)))) * ((vel_r[i] * (1.0 - (c_sr * c_sr))) +
        (c_sr * sqrt((1.0 - v_sq_r) * (inv_spatial_metric_r[i][i] * (1.0 - (v_sq_r * (c_sr * c_sr))) - (vel_r[i] * vel_r[i]) * (1.0 - (c_sr * c_sr)))))) - shift_r[i];

      slow_acoustic_eigs_r[i] = (lapse_r / (1.0 - (v_sq_r * (c_sr * c_sr)))) * ((vel_r[i] * (1.0 - (c_sr * c_sr))) -
        (c_sr * sqrt((1.0 - v_sq_r) * (inv_spatial_metric_r[i][i] * (1.0 - (v_sq_r * (c_sr * c_sr))) - (vel_r[i] * vel_r[i]) * (1.0 - (c_sr * c_sr)))))) - shift_r[i];
    }

    double max_eig_r = 0.0;
    for (int i = 0; i < 3; i++) {
      if (fabs(material_eigs_r[i]) > max_eig_r) {
        max_eig_r = fabs(material_eigs_r[i]);
      }
      if (fabs(fast_acoustic_eigs_r[i]) > max_eig_r) {
        max_eig_r = fabs(fast_acoustic_eigs_r[i]);
      }
      if (fabs(slow_acoustic_eigs_r[i]) > max_eig_r) {
        max_eig_r = fabs(slow_acoustic_eigs_r[i]);
      }
    }

    // Warning for later: these HLL signal speeds are NOT a correct bound on the characteristic cone. max_eig_{l,r} are already coordinate eigenvalues (lapse*v - shift +/- acoustic),
    // but they are then averaged and pushed through the relativistic velocity-addition formula a second time - double-counting the boost. The correct HLL bounds are the endpoint extrema of the coordinate
    // eigenvalues: sl = min(lambda^-_L, lambda^-_R), sr = max(lambda^+_L, lambda^+_R). 
    // The GR-Euler driver uses the LAX solver, so this path is unused for the TOV/collapse work;
    // HLL (and the flat branch below, which uses averaged states rather than endpoint min/max) must be reworked and validated before use with the curved 76-var system.
    double max_eig_avg = 0.5 * (max_eig_l + max_eig_r);
    
    sl = (vx_avg - max_eig_avg) / (1.0 - (vx_avg * max_eig_avg));
    sr = (vx_avg + max_eig_avg) / (1.0 + (vx_avg * max_eig_avg));
  }
  else {
    sl = (vx_avg - cs_avg) / (1.0 - (vx_avg * cs_avg));
    sr = (vx_avg + cs_avg) / (1.0 + (vx_avg * cs_avg));
  }

  double fl[71], fr[71];
  gkyl_gr_euler_flux(gas_gamma, ql, fl);
  gkyl_gr_euler_flux(gas_gamma, qr, fr);

  // Optional static-TOV well-balancing (identical to the Lax path): feed the HLL wave construction the equilibrium DEVIATION (qr-ql)-deq and (fr-fl)-dfeq. The standard HLL waves w0 = qm-ql, w1 = qr-qm are algebraically w0 = (sr*dq - df)/(sr-sl),
  // w1 = (df - sl*dq)/(sr-sl) with dq = qr-ql, df = fr-fl; so subtracting (deq, dfeq) keeps sl*w0 + sr*w1 = df = (fr-fl)-dfeq, i.e. amdq+apdq is unchanged -> strictly flux- and fluctuation-conservative (mass/energy: dfeq=0; momentum keeps its gravity-sourced balance).
  double deq[71] = { 0.0 }, dfeq[71] = { 0.0 };
  if (gr_euler->tov_eq && !gr_euler->disable_well_balanced && !in_excision_region_l && !in_excision_region_r) {
    double q_eq_l[71], q_eq_r[71], f_eq_l[71], f_eq_r[71];
    gr_euler_tov_equilibrium(gr_euler, ql, q_eq_l);
    gr_euler_tov_equilibrium(gr_euler, qr, q_eq_r);
    gkyl_gr_euler_flux(gas_gamma, q_eq_l, f_eq_l);
    gkyl_gr_euler_flux(gas_gamma, q_eq_r, f_eq_r);
    for (int i = 0; i < 71; i++) {
      deq[i] = q_eq_r[i] - q_eq_l[i];
      dfeq[i] = f_eq_r[i] - f_eq_l[i];
    }
  }

  double inv = 1.0 / (sr - sl);
  int neq = eqn->num_equations;
  double *w0 = &waves[0], *w1 = &waves[neq];
  if (!in_excision_region_l && !in_excision_region_r) {
    for (int i = 0; i < 71; i++) {
      double dq = (qr[i] - ql[i]) - deq[i];
      double df = (fr[i] - fl[i]) - dfeq[i];
      w0[i] = ((sr * dq) - df) * inv;
      w1[i] = (df - (sl * dq)) * inv;
    }
  }
  else {
    for (int i = 0; i < 71; i++) {
      w0[i] = 0.0;
      w1[i] = 0.0;
    }
  }
  // Static-TOV well-balancing only: the frozen-equilibrium slots q[71],q[72] carry zero waves, so they never evolve
  if (gr_euler->tov_eq != NULL) {
    for (int i = 71; i < neq; i++) {
      w0[i] = 0.0;
      w1[i] = 0.0;
    }
  }

  s[0] = sl;
  s[1] = sr;

  return fmax(fabs(sl), fabs(sr));
}

static void
qfluct_hll(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
{
  int neq = eqn->num_equations;
  const double *w0 = &waves[0], *w1 = &waves[neq];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i = 0; i < neq; i++) { // frozen slots (i>=71) have zero waves -> amdq=apdq=0 (frozen)
    amdq[i] = (s0m * w0[i]) + (s1m * w1[i]);
    apdq[i] = (s0p * w0[i]) + (s1p * w1[i]);
  }
}

static double
wave_hll_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* delta, const double* ql, const double* qr, const double phil, const double phir, double* waves, double* s)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX) {
    return wave_hll(eqn, delta, ql, qr, waves, s);
  }
  else {
    return wave_lax(eqn, delta, ql, qr, waves, s);
  }

  return 0.0; // Unreachable code.
}

static void
qfluct_hll_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* ql, const double* qr, const double phil, const double phir, const double* waves, const double* s,
  double* amdq, double* apdq)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX) {
    return qfluct_hll(eqn, ql, qr, waves, s, amdq, apdq);
  }
  else {
    return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
  }
}

// HLL/Lax: HLL in the bulk, Lax near the atmosphere/surface. 
// Both wave constructions keep amdq+apdq=df (the qfluct is identical), so the blend is strictly flux- and fluctuation-conservative. 
static double
wave_hll_lax(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
  double vl[71], vr[71];
  gkyl_gr_euler_prim_vars(gr_euler->gas_gamma, ql, vl);
  gkyl_gr_euler_prim_vars(gr_euler->gas_gamma, qr, vr);

  double rho_min = fmin(vl[0], vr[0]);
  if (rho_min < GR_EULER_BLEND_FACTOR * gr_euler_rho_floor) {
    return wave_lax(eqn, delta, ql, qr, waves, s);
  }
  return wave_hll(eqn, delta, ql, qr, waves, s);
}

static double
wave_hll_lax_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* delta, const double* ql, const double* qr, const double phil, const double phir, double* waves, double* s)
{
  return wave_hll_lax(eqn, delta, ql, qr, waves, s);
}

static double
flux_jump(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, double* flux_jump)
{
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
  double gas_gamma = gr_euler->gas_gamma;

  double fr[71], fl[71];
  gkyl_gr_euler_flux(gas_gamma, ql, fl);
  gkyl_gr_euler_flux(gas_gamma, qr, fr);

  bool in_excision_region_l = false;
  if (ql[27] < pow(10.0, -8.0)) {
    in_excision_region_l = true;
  }

  bool in_excision_region_r = false;
  if (qr[27] < pow(10.0, -8.0)) {
    in_excision_region_r = true;
  }

  if (!in_excision_region_l && !in_excision_region_r) {
    for (int m = 0; m < 71; m++) {
      flux_jump[m] = fr[m] - fl[m];
    }
    // subtract the equilibrium flux jump so this within-cell flux difference becomes the deviation flux jump
    if (gr_euler->tov_eq != NULL && !gr_euler->disable_well_balanced) { // called by the MP scheme
      double q_eq_l[71], q_eq_r[71], f_eq_l[71], f_eq_r[71];
      gr_euler_tov_equilibrium(gr_euler, ql, q_eq_l);
      gr_euler_tov_equilibrium(gr_euler, qr, q_eq_r);
      gkyl_gr_euler_flux(gas_gamma, q_eq_l, f_eq_l);
      gkyl_gr_euler_flux(gas_gamma, q_eq_r, f_eq_r);
      for (int m = 0; m < 71; m++) {
        flux_jump[m] -= (f_eq_r[m] - f_eq_l[m]);
      }
    }
  }
  else {
    for (int m = 0; m < 71; m++) {
      flux_jump[m] = 0.0;
    }
  }
  // Frozen WB slots (>= 71) carry no flux
  for (int m = 71; m < eqn->num_equations; m++) {
    flux_jump[m] = 0.0;
  }

  double amaxl = gkyl_gr_euler_max_abs_speed(gas_gamma, ql);
  double amaxr = gkyl_gr_euler_max_abs_speed(gas_gamma, qr);

  return fmax(amaxl, amaxr);
}

static bool
check_inv(const struct gkyl_wv_eqn* eqn, const double* q)
{
  // The forward Valencia cone also requires E=D+tau>0; squaring alone
  // admits its unphysical negative-energy branch.
  // sqrt(det) > 0 cancels in the sign test, so the densitized slots can be used directly.
  if (!isfinite(q[27])) return false;
  if (q[27] < pow(10.0, -8.0)) {
    return true; // excised cell: not an evolved fluid cell, nothing to check
  }
  for (int m=0; m<eqn->num_equations; ++m)
    if (!isfinite(q[m])) return false;
  double DE = q[0] + q[4]; // sqrt(gamma) (D + tau)
  if (!(isfinite(DE) && DE > 0.0 && q[0] > 0.0 && q[5] > 0.0))
    return false;
  if (!gkyl_mp_scheme_gr_metric_is_spd(&q[9])) return false;
  double spatial_metric[3][3];
  spatial_metric[0][0] = q[9];  spatial_metric[0][1] = q[10]; spatial_metric[0][2] = q[11];
  spatial_metric[1][0] = q[12]; spatial_metric[1][1] = q[13]; spatial_metric[1][2] = q[14];
  spatial_metric[2][0] = q[15]; spatial_metric[2][1] = q[16]; spatial_metric[2][2] = q[17];
  double spatial_det = gr_euler_spatial_metric_det(spatial_metric);
  if (!(isfinite(spatial_det) && spatial_det > 0.0)) {
    return false;
  }
  double inv_spatial_metric[3][3];
  gr_euler_spatial_metric_inv(spatial_metric, spatial_det, inv_spatial_metric);

  // Scaling keeps the cone test finite for large, otherwise valid states.
  double scale = fmax(DE, fmax(q[0], fmax(fabs(q[1]), fmax(fabs(q[2]), fabs(q[3])))));
  double D = q[0]/scale, E = DE/scale;
  double S[3] = {q[1]/scale, q[2]/scale, q[3]/scale};
  double S2 = 0.0;
  for (int si = 0; si < 3; si++)
    for (int sj = 0; sj < 3; sj++)
      S2 += inv_spatial_metric[si][sj] * S[si] * S[sj];
  return isfinite(S2) && S2 >= 0.0 && (E * E) - (D * D) - S2 > 0.0;
}

static double
max_speed(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
  double gas_gamma = gr_euler->gas_gamma;

  return gkyl_gr_euler_max_abs_speed(gas_gamma, q);
}

static inline void
gr_euler_cons_to_diag(const struct gkyl_wv_eqn* eqn, const double* qin, double* diag)
{
  for (int i = 0; i < 5; i++) {
    diag[i] = qin[i];
  }
}

// DEAD CODE for all problems (function below)
// Raw geometric/gravitational source (no well-balancing); factored so it can also be evaluated at the frozen TOV equilibrium state for the WB subtraction below.
static inline void
gr_euler_raw_source(double gas_gamma, const double* qin, double* sout)
{
  double v[71] = { 0.0 };
  gkyl_gr_euler_prim_vars(gas_gamma, qin, v);
  double rho = v[0];
  double vx = v[1];
  double vy = v[2];
  double vz = v[3];
  double p = v[4];

  double lapse = v[5];
  double shift_x = v[6];
  double shift_y = v[7];
  double shift_z = v[8];

  double spatial_metric[3][3];
  spatial_metric[0][0] = v[9]; spatial_metric[0][1] = v[10]; spatial_metric[0][2] = v[11];
  spatial_metric[1][0] = v[12]; spatial_metric[1][1] = v[13]; spatial_metric[1][2] = v[14];
  spatial_metric[2][0] = v[15]; spatial_metric[2][1] = v[16]; spatial_metric[2][2] = v[17];

  double extrinsic_curvature[3][3];
  extrinsic_curvature[0][0] = v[18]; extrinsic_curvature[0][1] = v[19]; extrinsic_curvature[0][2] = v[20];
  extrinsic_curvature[1][0] = v[21]; extrinsic_curvature[1][1] = v[22]; extrinsic_curvature[1][2] = v[23];
  extrinsic_curvature[2][0] = v[24]; extrinsic_curvature[2][1] = v[25]; extrinsic_curvature[2][2] = v[26];

  double lapse_der[3];
  lapse_der[0] = v[28];
  lapse_der[1] = v[29];
  lapse_der[2] = v[30];

  double shift_der[3][3];
  shift_der[0][0] = v[31]; shift_der[0][1] = v[32]; shift_der[0][2] = v[33];
  shift_der[1][0] = v[34]; shift_der[1][1] = v[35]; shift_der[1][2] = v[36];
  shift_der[2][0] = v[37]; shift_der[2][1] = v[38]; shift_der[2][2] = v[39];

  double spatial_metric_der[3][3][3];
  spatial_metric_der[0][0][0] = v[40]; spatial_metric_der[0][0][1] = v[41]; spatial_metric_der[0][0][2] = v[42];
  spatial_metric_der[0][1][0] = v[43]; spatial_metric_der[0][1][1] = v[44]; spatial_metric_der[0][1][2] = v[45];
  spatial_metric_der[0][2][0] = v[46]; spatial_metric_der[0][2][1] = v[47]; spatial_metric_der[0][2][2] = v[48];

  spatial_metric_der[1][0][0] = v[49]; spatial_metric_der[1][0][1] = v[50]; spatial_metric_der[1][0][2] = v[51];
  spatial_metric_der[1][1][0] = v[52]; spatial_metric_der[1][1][1] = v[53]; spatial_metric_der[1][1][2] = v[54];
  spatial_metric_der[1][2][0] = v[55]; spatial_metric_der[1][2][1] = v[56]; spatial_metric_der[1][2][2] = v[57];

  spatial_metric_der[2][0][0] = v[58]; spatial_metric_der[2][0][1] = v[59]; spatial_metric_der[2][0][2] = v[60];
  spatial_metric_der[2][1][0] = v[61]; spatial_metric_der[2][1][1] = v[62]; spatial_metric_der[2][1][2] = v[63];
  spatial_metric_der[2][2][0] = v[64]; spatial_metric_der[2][2][1] = v[65]; spatial_metric_der[2][2][2] = v[66];

  double **stress_energy = gkyl_malloc(sizeof(double*[4]));
  for (int i = 0; i < 4; i++) {
    stress_energy[i] = gkyl_malloc(sizeof(double[4]));
  }

  gkyl_gr_euler_stress_energy_tensor(gas_gamma, qin, &stress_energy);

  bool in_excision_region = false;
  if (v[27] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double shift[3];
    shift[0] = shift_x; shift[1] = shift_y; shift[2] = shift_z;

    double vel[3];
    double v_sq = 0.0;
    vel[0] = vx; vel[1] = vy; vel[2] = vz;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        v_sq += spatial_metric[i][j] * vel[i] * vel[j];
      }
    }

    double W = 1.0 / (sqrt(1.0 - v_sq));
    if (v_sq > 1.0 - pow(10.0, -8.0)) {
      W = 1.0 / sqrt(1.0 - pow(10.0, -8.0));
    }

    double h = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0)));
    
    double cov_vel[3] = { 0.0 };
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        cov_vel[i] += spatial_metric[i][j] * vel[j];
      }
    }

    double mom[3];
    mom[0] = (rho * h) * (W * W) * cov_vel[0];
    mom[1] = (rho * h) * (W * W) * cov_vel[1];
    mom[2] = (rho * h) * (W * W) * cov_vel[2];

    // Energy density source.
    sout[4] = 0.0;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        sout[4] += stress_energy[0][0] * shift[i] * shift[j] * extrinsic_curvature[i][j];
        sout[4] += 2.0 * stress_energy[0][i + 1] * shift[j] * extrinsic_curvature[i][j];
        sout[4] += stress_energy[i + 1][j + 1] * extrinsic_curvature[i][j];
      }

      sout[4] += stress_energy[0][0] * shift[i] * lapse_der[i];
      sout[4] -= stress_energy[0][i + 1] * lapse_der[i];
    }

    // Momentum density sources.
    for (int j = 0; j < 3; j++) {
      sout[1 + j] = -stress_energy[0][0] * lapse * lapse_der[j];

      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          sout[1 + j] += 0.5 * stress_energy[0][0] * shift[k] * shift[l] * spatial_metric_der[j][k][l];
          sout[1 + j] += 0.5 * stress_energy[k + 1][l + 1] * spatial_metric_der[j][k][l];
        }

        sout[1 + j] += (mom[k] / lapse) * shift_der[j][k];

        for (int i = 0; i < 3; i++) {
          sout[1 + j] += stress_energy[0][i + 1] * shift[k] * spatial_metric_der[j][i][k];
        }
      }
    }
  }
  else {
    for (int i = 0; i < 71; i++) {
      sout[i] = 0.0;
    }
  }

  for (int i = 0; i < 4; i++) {
    gkyl_free(stress_energy[i]);
  }
  gkyl_free(stress_energy);
}

static void
gr_euler_source(const struct gkyl_wv_eqn* eqn, const double* qin, double* sout)
{
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);

  gr_euler_raw_source(gr_euler->gas_gamma, qin, sout);

  // Optional static-TOV well-balancing: subtract the equilibrium geometric source sigma_eq = source(q_eq), so a static star leaves the (gravity-sourced) momentum unmoved.
  // Only the fluid components (1..4) are touched; mass (0) and the metric (5..70) carry no source.
  if (gr_euler->tov_eq && !gr_euler->disable_well_balanced) {
    double q_eq[71];
    gr_euler_tov_equilibrium(gr_euler, qin, q_eq);
    double sout_eq[71] = { 0.0 };
    gr_euler_raw_source(gr_euler->gas_gamma, q_eq, sout_eq);
    for (int i = 0; i < 71; i++) {
      sout[i] -= sout_eq[i];
    }
  }
  for (int i = 71; i < eqn->num_equations; i++) {
    sout[i] = 0.0;
  }
}

void
gkyl_gr_euler_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_wv_eqn* base = container_of(ref, struct gkyl_wv_eqn, ref_count);

  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct wv_gr_euler *gr_euler = container_of(base->on_dev, struct wv_gr_euler, eqn);
    gkyl_cu_free(gr_euler);
  }

  struct wv_gr_euler *gr_euler = container_of(base, struct wv_gr_euler, eqn);
  gkyl_free(gr_euler);
}

struct gkyl_wv_eqn*
gkyl_wv_gr_euler_new(double gas_gamma, enum gkyl_spacetime_gauge spacetime_gauge, int reinit_freq, struct gkyl_gr_spacetime* spacetime, bool use_gpu)
{
  return gkyl_wv_gr_euler_inew(&(struct gkyl_wv_gr_euler_inp) {
      .gas_gamma = gas_gamma,
      .spacetime_gauge = spacetime_gauge,
      .reinit_freq = reinit_freq,
      .spacetime = spacetime,
      .rp_type = WV_GR_EULER_RP_HLL,
      .use_gpu = use_gpu,
    }
  );
}

struct gkyl_wv_eqn*
gkyl_wv_gr_euler_inew(const struct gkyl_wv_gr_euler_inp* inp)
{
  struct wv_gr_euler *gr_euler = gkyl_malloc(sizeof(struct wv_gr_euler));

  gr_euler->eqn.type = GKYL_EQN_GR_EULER;
  if (inp->tov_eq != NULL) {
    // Well-balancing state: physical+metric (0..70) + frozen conserved equilibrium D_eq/Etot_eq/S_eq (71..75).
    gr_euler->eqn.num_equations = 76;
  }
  else {
    gr_euler->eqn.num_equations = 71;
  }
  gr_euler->eqn.num_diag = 5;

  gr_euler->gas_gamma = inp->gas_gamma;
  gr_euler->spacetime = inp->spacetime;
  gr_euler->spacetime_gauge = inp->spacetime_gauge;
  gr_euler->reinit_freq = inp->reinit_freq;
  gr_euler->tov_eq = inp->tov_eq; // Optional static-TOV well-balancing (NULL -> inactive) (the equilibrium comes from the discrete slots, so no EOS constants needed)

  gr_euler->disable_well_balanced = inp->disable_well_balanced;

  if (inp->rho_atm > 0.0) gr_euler_rho_floor = inp->rho_atm;
  if (inp->p_atm > 0.0) gr_euler_p_floor = inp->p_atm;
  gr_euler->rho_atm = gr_euler_rho_floor;
  gr_euler->p_atm = gr_euler_p_floor;

  if (inp->rp_type == WV_GR_EULER_RP_LAX) {
    gr_euler->eqn.num_waves = 2;
    gr_euler->eqn.waves_func = wave_lax_l;
    gr_euler->eqn.qfluct_func = qfluct_lax_l;
  }
  else if (inp->rp_type == WV_GR_EULER_RP_ROE) {
    gr_euler->eqn.num_waves = 3;
    gr_euler->eqn.waves_func = wave_roe_l;
    gr_euler->eqn.qfluct_func = qfluct_roe_l;
  }
  else if (inp->rp_type == WV_GR_EULER_RP_HLL) {
    gr_euler->eqn.num_waves = 2;
    gr_euler->eqn.waves_func = wave_hll_l;
    gr_euler->eqn.qfluct_func = qfluct_hll_l;
  }
  else if (inp->rp_type == WV_GR_EULER_RP_HLL_LAX) {
    gr_euler->eqn.num_waves = 2;
    gr_euler->eqn.waves_func = wave_hll_lax_l;
    gr_euler->eqn.qfluct_func = qfluct_hll_l;  // identical to qfluct_lax;
  }

  gr_euler->eqn.flux_jump = flux_jump;
  gr_euler->eqn.check_inv_func = check_inv;
  gr_euler->eqn.max_speed_func = max_speed;
  gr_euler->eqn.rotate_to_local_func = rot_to_local;
  gr_euler->eqn.rotate_to_global_func = rot_to_global;

  gr_euler->eqn.wall_bc_func = gr_euler_wall;
  gr_euler->eqn.no_slip_bc_func = gr_euler_no_slip;

  gr_euler->eqn.cons_to_riem = cons_to_riem;
  gr_euler->eqn.riem_to_cons = riem_to_cons;

  gr_euler->eqn.cons_to_diag = gr_euler_cons_to_diag;

  gr_euler->eqn.source_func = gr_euler_source;

  gr_euler->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(gr_euler->eqn.flags);
  gr_euler->eqn.ref_count = gkyl_ref_count_init(gkyl_gr_euler_free);
  gr_euler->eqn.on_dev = &gr_euler->eqn; // On the CPU, the equation object points to itself.

  gr_euler->eqn.embed_geo = NULL;

  return &gr_euler->eqn;
}

double
gkyl_wv_gr_euler_gas_gamma(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
  double gas_gamma = gr_euler->gas_gamma;

  return gas_gamma;
}

enum gkyl_spacetime_gauge
gkyl_wv_gr_euler_spacetime_gauge(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
  enum gkyl_spacetime_gauge spacetime_gauge = gr_euler->spacetime_gauge;

  return spacetime_gauge;
}

int
gkyl_wv_gr_euler_reinit_freq(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
  int reinit_freq = gr_euler->reinit_freq;

  return reinit_freq;
}

struct gkyl_gr_spacetime*
gkyl_wv_gr_euler_spacetime(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
  struct gkyl_gr_spacetime *spacetime = gr_euler->spacetime;

  return spacetime;
}
