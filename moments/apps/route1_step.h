// Route 1: method-of-lines SSP-RK3 step for the conformal spacetime (diagnostic build only).
// Included into a build-time copy of mom_update_mixed.c (see build.py); GKYL_ROUTE1 unset leaves the
// production step untouched.
//
// One SSP-RK3 step of dq/dt = L(q) - L(e) replaces production's driver/source/transport/source/driver
// composition. e is the spacetime well-balancing reference production already uses (equilibrium state,
// equilibrium source S_eq), so e stays an exact fixed point exactly as in production.
//   rows gt, alpha, K, D, A, V: central face flux 1/2 (F(q_i) + F(q_{i+1})) (single-valued, conservative)
//                               + vacuum source + matter source (frozen midpoint fluid) - S_eq
//   alpha: evolved as ln alpha (keeps C_A = A - D0 ln alpha exact under RK)
//   chi: pure source S_pure = S_64 + beta.P + 2 chi B_kk (route-1 design, identity_test.c)
//   P_k: D0_k S_pure (P flux = -S_pure);  H = sym D0 L(P);  B = 1/2 D0 L(beta)
//   beta: Gamma-driver rate mu (Gamma - Gamma_ref) - eta (beta - beta_ref) (production driver, gated)
//   KO: matched conservative face flux sigma/(16h) delta^4 on q - e (ln alpha for alpha), all evolved
//       rows; G = 0 on faces whose 4-cell stencil leaves the interior.
// Env: GKYL_ROUTE1=1, GKYL_ROUTE1_KO (sigma, default 0.3), GKYL_ROUTE1_KO_GHOSTS (KO up to the boundary,
// default 0), GKYL_ROUTE1_KO_SBP (boundary-modified MSN dissipation instead of the cut KO; see
// route1_ko_sbp), GKYL_ROUTE1_CHECK (report reduction-identity drift every n steps, default 20),
// GKYL_ROUTE1_NEWRAD (radiative ghost cells evolved with the stages, default 0; see route1_newrad.h),
// GKYL_ROUTE1_BAND (upwind band width in cells, default 0 = off; GKYL_ROUTE1_BAND_S overrides its speed;
// GKYL_ROUTE1_BAND_DERIVED=1 bands only the parent rows and gives D, A, V the increments their identities
// imply), GKYL_ROUTE1_MONITOR (boundary monitor interval in time units, default 0 = off; with
// GKYL_ROUTE1_IDB=1 it also reports C_D, C_A, C_V by cell depth), GKYL_ROUTE1_PEEL (the peeled boundary:
// NewRad on the parents, derived rows on the ghosts, KO to the boundary, centred B; refuses band and buffer;
// see route1_newrad.h), GKYL_ROUTE1_DRIVER_METRIC (the driver's Gammatilde from 1/2 D0 gt instead of D; see
// route1_driver_metric), GKYL_ROUTE1_CHARBC (the characteristic boundary flux: B and P in flux form, F(e + Pi (q - e))
// on the outer faces, no ghost reads; see route1_charbc.h). Every coupled step of this build also keeps the rest-mass ledger (route1_ledger_*); the
// fluid switches GKYL_ROUTE1_BOOST and GKYL_ROUTE1_FLOOR_ONLY are described above route1_boost_enabled.

#include <gkyl_wv_vacuum_einstein_conformal_priv.h>

// [FX] Fluid-only excision inside the apparent horizon (this build only; GKYL_FLUID_EXCISION=f with 0 < f < 1;
// unset = off and the step is the frozen build's). Hydro excision as in Hawke, Loffler & Nerozzi, PRD 71, 104006
// (2005): the spacetime keeps evolving everywhere, only the fluid is cut out.
// - Mask: once the driver's spherical AH certificate has reported a valid root R_AH in 3 consecutive scans
//   (route1_fx_ah), the cells with r < f R_AH. The mask only grows, and only when the legality check below holds on
//   its new boundary.
// - Excised cells keep their stored state frozen (no direct write; their RHS is dropped), so their matter still
//   sources gravity exactly as before they were excised.
// - Each forward-Euler stage evaluates the RHS on a scratch copy whose excised cells hold the full row of the live
//   cell just outside the mask on the same ray (radial outflow fill). A face between a live and an excised cell
//   then carries the live-side flux, and the live cells' update is the flux-conservative one.
// - The excised cells' D-row RHS sums (faces inside the mask telescope) to the D flux through the mask surface. It
//   goes into the black-hole reservoir D_BH, combined with the SSP-RK3 weights and committed with the accepted step,
//   so sum D (grid) + D_BH changes only by the outer-boundary flux ([FX] lines).
// - Legality: on every live cell next to the mask, the outgoing radial light speed
//   -beta^i xi_i + alpha sqrt(gamma^ij xi_i xi_j) (xi = x/r) must be negative. It bounds every fluid
//   characteristic, so nothing leaves the mask. A violation is logged as [FX-ILLEGAL].
// GKYL_FX_DRIVER_REF_SAVE / _LOAD (restart support): write / read the Gamma-driver reference captured at activation.
static double route1_fx_f = -2.0, route1_fx_r = 0.0, route1_fx_req = 0.0;
static int route1_fx_nvalid = 0;
static long route1_fx_n = 0, route1_fx_illegal = 0, route1_fx_nstep = 0;
static struct gkyl_array *route1_fx_mask = NULL, *route1_fx_src = NULL, *route1_fx_scr = NULL;
static double route1_fx_res = 0.0, route1_fx_stage[3], route1_fx_pend = 0.0, route1_fx_bh = 0.0;
static double route1_fx_lam = -INFINITY, route1_fx_lam_r = 0.0;

// [FXR] GKYL_FX_REMOVE=1 (with GKYL_FLUID_EXCISION): matter removal at masking, the published practice for collapse to
// a black hole (matter deep inside the horizon reset to atmosphere, as in the ETK spherical-coordinates code; the hydro
// variables excised in the 1D collapse of Bernuzzi et al.). Each cell that joins the mask has its rows D, S_i, tau
// moved into the black-hole reservoir (D_BH, S_BH, tau_BH) and set to the static atmosphere of the primitive floor
// (rho_atm, p_atm, v = 0: D = sqrt(gamma) rho_atm, S_i = 0, tau = sqrt(gamma) p_atm / (Gamma - 1)). This happens once,
// right after the accepted step in which the cell joined the mask, after that step's ledger sums: the ledger books it
// in "between", and sum D (grid) + D_BH does not change ([FX-REMOVE]). The cell then stays frozen (RHS dropped) as
// before, and the D flux into the mask keeps going into D_BH. Unset: the frozen fx build.
static int route1_fxr = -1;
static struct gkyl_array *route1_fx_new = NULL;
static double route1_fx_sbh[3] = { 0.0, 0.0, 0.0 }, route1_fx_tbh = 0.0;

static bool
route1_fx_remove_on(void)
{
  if (route1_fxr < 0) {
    const char *v = getenv("GKYL_FX_REMOVE");
    route1_fxr = v && atoi(v) != 0;
  }
  return route1_fxr;
}

static double
route1_fx_frac(void)
{
  if (route1_fx_f < -1.0) {
    const char *v = getenv("GKYL_FLUID_EXCISION");
    route1_fx_f = v ? atof(v) : 0.0;
    if (!(route1_fx_f > 0.0 && route1_fx_f < 1.0)) route1_fx_f = 0.0;
  }
  return route1_fx_f;
}

// [HKO] KO inside the apparent horizon (this build only; GKYL_ROUTE1_KO_AH=sigma_in, unset = off and the step is the
// fx build's). As in Baiotti & Rezzolla, PRL 97, 141101 (2006), the KO coefficient may vary in space: sigma_in on
// every face that touches a cell with r < f_ko R_AH (GKYL_ROUTE1_KO_AH_F, default 1), the global sigma elsewhere.
// The radius comes from the driver's spherical AH certificate after 3 consecutive valid roots and only grows.
// The coefficient is a face value (max of its two cells), so the KO term stays a face difference.
static double route1_hko_sig = -1.0, route1_hko_f = 1.0, route1_hko_r = 0.0;
static int route1_hko_nv = 0;

// [RUS] GKYL_ROUTE1_RUS_AH=1 (this build only): on the same faces a Rusanov term 1/2 s (v_R - v_L) replaces the KO
// term (v = the KO variables: q - e, ln alpha for alpha; the KO rows). s = the larger of the two cells' fastest speed:
// the library's max |characteristic speed| (shift and 1+log included) or |beta^d| + sqrt(4 mu/3 gt^dd), the
// longitudinal Gamma-driver speed (mu local: kappa alpha^2 chi^2 when weighted). One value per face: a face difference.
static int route1_hko_rus = -1;

static bool
route1_rus(void)
{
  if (route1_hko_rus < 0) {
    const char *v = getenv("GKYL_ROUTE1_RUS_AH");
    route1_hko_rus = v && atoi(v) != 0;
  }
  return route1_hko_rus;
}

static double
route1_hko_sigma(void)
{
  if (route1_hko_sig < 0.0) {
    const char *v = getenv("GKYL_ROUTE1_KO_AH"), *f = getenv("GKYL_ROUTE1_KO_AH_F");
    route1_hko_sig = v ? atof(v) : 0.0;
    if (!(route1_hko_sig > 0.0)) route1_hko_sig = 0.0;
    route1_hko_f = f ? atof(f) : 1.0;
  }
  return route1_hko_sig;
}

static void
route1_hko_ah(double t, int valid, double R)
{
  if (!(route1_hko_sigma() > 0.0) && !route1_rus()) return;
  route1_hko_nv = (valid && R > 0.0) ? route1_hko_nv + 1 : 0;
  if (route1_hko_nv >= 3 && route1_hko_f*R > route1_hko_r) {
    route1_hko_r = route1_hko_f*R;
    if (route1_rus())
      fprintf(stderr, "[RUS] t=%.6g Rusanov 1/2 s (v_R - v_L) instead of KO on faces touching r < %.4f (%g x R_AH)\n",
        t, route1_hko_r, route1_hko_f);
    else
      fprintf(stderr, "[HKO] t=%.6g KO sigma %g on faces touching r < %.4f (%g x R_AH)\n", t, route1_hko_sig,
        route1_hko_r, route1_hko_f);
  }
}

// [RUS] fastest |speed| of a cell (any direction): library max speed, or |beta^d| + sqrt(4 mu/3 gt^dd).
static double
route1_rus_speed(const struct moment_species *esp, const double *q, double mu)
{
  double s = gkyl_wv_eqn_max_speed(esp->equation, q);
  if (!isfinite(s)) s = 0.0;
  const double *g = q; // conformal metric, slots 0..8 (row-major 3x3)
  double c00 = g[4]*g[8] - g[5]*g[7], c11 = g[0]*g[8] - g[2]*g[6], c22 = g[0]*g[4] - g[1]*g[3];
  double det = g[0]*c00 - g[1]*(g[3]*g[8] - g[5]*g[6]) + g[2]*(g[3]*g[7] - g[4]*g[6]);
  double gii[3] = { c00/det, c11/det, c22/det };
  for (int d = 0; d < 3; ++d) {
    double v = fabs(q[52 + d]) + sqrt(fmax(0.0, 4.0*mu/3.0*gii[d]));
    if (isfinite(v) && v > s) s = v;
  }
  return s;
}

// Called by the driver after each [AHCERT-ROOT] line.
void
route1_fx_ah(double t, int valid, double R)
{
  route1_hko_ah(t, valid, R);
  if (!(route1_fx_frac() > 0.0)) return;
  route1_fx_nvalid = (valid && R > 0.0) ? route1_fx_nvalid + 1 : 0;
  if (route1_fx_nvalid >= 3 && route1_fx_f*R > route1_fx_req) route1_fx_req = route1_fx_f*R;
}

static bool
route1_fx_masked(long l)
{
  return route1_fx_mask && ((const double *) gkyl_array_cfetch(route1_fx_mask, l))[0] > 0.0;
}

static double
route1_fx_rad(const gkyl_moment_app *app, const int *idx, double *xc)
{
  xc[0] = xc[1] = xc[2] = 0.0;
  gkyl_rect_grid_cell_center(&app->grid, idx, xc);
  double r = 0.0;
  for (int d = 0; d < app->ndim; ++d) r += xc[d]*xc[d];
  return sqrt(r);
}

// Outgoing radial light speed of fluid row q (lapse 5, shift 6-8, metric 9-17) at xc, r > 0.
static double
route1_fx_lambda(const double *q, const double *xc, double r)
{
  double g[3][3] = { { q[9], q[10], q[11] }, { q[12], q[13], q[14] }, { q[15], q[16], q[17] } }, gi[3][3];
  sym3_inverse(g, gi);
  double xi[3] = { xc[0]/r, xc[1]/r, xc[2]/r }, nn = 0.0, bn = 0.0;
  for (int i = 0; i < 3; ++i) {
    bn += q[6 + i]*xi[i];
    for (int j = 0; j < 3; ++j) nn += gi[i][j]*xi[i]*xi[j];
  }
  return -bn + q[5]*sqrt(fmax(nn, 0.0));
}

// Max outgoing light speed over the live cells (r >= rr) with a face neighbour inside rr.
static double
route1_fx_legal(const gkyl_moment_app *app, const struct gkyl_array *f, double rr, double *rat)
{
  double lmax = -INFINITY;
  *rat = 0.0;
  struct gkyl_range_iter it;
  gkyl_range_iter_init(&it, &app->local);
  while (gkyl_range_iter_next(&it)) {
    double xc[3], xn[3];
    double r = route1_fx_rad(app, it.idx, xc);
    if (r < rr) continue;
    bool adj = false;
    for (int d = 0; d < app->ndim && !adj; ++d)
      for (int s = -1; s <= 1 && !adj; s += 2) {
        int j[3] = { it.idx[0], it.idx[1], it.idx[2] };
        j[d] += s;
        if (j[d] < app->local.lower[d] || j[d] > app->local.upper[d]) continue;
        adj = route1_fx_rad(app, j, xn) < rr;
      }
    if (!adj) continue;
    double lam = route1_fx_lambda(gkyl_array_cfetch(f, gkyl_range_idx(&app->local_ext, it.idx)), xc, r);
    if (lam > lmax) { lmax = lam; *rat = r; }
  }
  return lmax;
}

// Grow the mask to the requested radius when legal; check the current boundary every call.
static void
route1_fx_update(gkyl_moment_app *app, struct moment_species *fsp, double t)
{
  if (!(route1_fx_frac() > 0.0)) return;
  double rat;
  if (route1_fx_req > route1_fx_r) {
    double lam = route1_fx_legal(app, fsp->fcurr, route1_fx_req, &rat);
    if (!(lam < 0.0)) {
      gkyl_moment_app_cout(app, stderr, "[FX] t=%.6g growth to r=%.4f refused: outgoing light speed %+.4e at r=%.4f\n",
        t, route1_fx_req, lam, rat);
      route1_fx_req = route1_fx_r;
    }
    else {
      long vol = app->local_ext.volume;
      if (!route1_fx_mask) {
        route1_fx_mask = gkyl_array_new(GKYL_DOUBLE, 1, vol);
        route1_fx_src = gkyl_array_new(GKYL_DOUBLE, 1, vol);
        route1_fx_scr = gkyl_array_new(GKYL_DOUBLE, fsp->num_equations, vol);
        gkyl_array_clear(route1_fx_mask, 0.0);
        gkyl_array_clear(route1_fx_src, -1.0);
        if (route1_fx_remove_on()) {
          route1_fx_new = gkyl_array_new(GKYL_DOUBLE, 1, vol);
          gkyl_array_clear(route1_fx_new, 0.0);
        }
      }
      route1_fx_r = route1_fx_req;
      double dxm = app->grid.dx[0];
      long n = 0;
      struct gkyl_range_iter it;
      gkyl_range_iter_init(&it, &app->local);
      while (gkyl_range_iter_next(&it)) {
        double xc[3];
        double r = route1_fx_rad(app, it.idx, xc);
        if (r < route1_fx_r) {
          long lm = gkyl_range_idx(&app->local_ext, it.idx);
          double *mk = gkyl_array_fetch(route1_fx_mask, lm);
          if (route1_fx_new && !(mk[0] > 0.0)) ((double *) gkyl_array_fetch(route1_fx_new, lm))[0] = 1.0; // [FXR]
          mk[0] = 1.0;
          n += 1;
        }
      }
      gkyl_range_iter_init(&it, &app->local);
      while (gkyl_range_iter_next(&it)) {
        long l = gkyl_range_idx(&app->local_ext, it.idx);
        if (!route1_fx_masked(l)) continue;
        double xc[3], r = route1_fx_rad(app, it.idx, xc), u[3] = { 1.0, 0.0, 0.0 };
        if (r > 0.0) for (int d = 0; d < 3; ++d) u[d] = xc[d]/r;
        long src = -1;
        for (int k = 0; k < 40 && src < 0; ++k) {
          double rs = route1_fx_r + (0.75 + 0.5*k)*dxm;
          int j[3] = { 0, 0, 0 };
          for (int d = 0; d < app->ndim; ++d) {
            int v = app->local.lower[d] + (int) floor((rs*u[d] - app->grid.lower[d])/app->grid.dx[d]);
            j[d] = v < app->local.lower[d] ? app->local.lower[d] : v > app->local.upper[d] ? app->local.upper[d] : v;
          }
          long lj = gkyl_range_idx(&app->local_ext, j);
          if (!route1_fx_masked(lj)) src = lj;
        }
        ((double *) gkyl_array_fetch(route1_fx_src, l))[0] = (double) src;
      }
      route1_fx_n = n;
      double lam = route1_fx_legal(app, fsp->fcurr, route1_fx_r, &rat);
      gkyl_moment_app_cout(app, stderr, "[FX] t=%.6g mask r < %.4f (f = %g x R_AH): %ld cells; outgoing light speed"
        " at its boundary <= %+.4e (r=%.4f); D_BH=%.10e\n", t, route1_fx_r, route1_fx_f, n, lam, rat, route1_fx_bh);
    }
  }
  if (route1_fx_r > 0.0) {
    route1_fx_lam = route1_fx_legal(app, fsp->fcurr, route1_fx_r, &route1_fx_lam_r);
    if (!(route1_fx_lam < 0.0)) {
      route1_fx_illegal += 1;
      if (route1_fx_illegal <= 20 || route1_fx_illegal % 100 == 0)
        gkyl_moment_app_cout(app, stderr, "[FX-ILLEGAL] t=%.6g outgoing light speed %+.4e at r=%.4f on the mask"
          " boundary r=%.4f (count %ld)\n", t, route1_fx_lam, route1_fx_lam_r, route1_fx_r, route1_fx_illegal);
    }
  }
}

// RHS input with the outflow fill (the scratch copy), or fin itself when no mask is active.
static const struct gkyl_array *
route1_fx_fill(const gkyl_moment_app *app, const struct gkyl_array *fin)
{
  if (!(route1_fx_r > 0.0)) return fin;
  gkyl_array_copy(route1_fx_scr, fin);
  size_t esz = fin->esznc;
  struct gkyl_range_iter it;
  gkyl_range_iter_init(&it, &app->local);
  while (gkyl_range_iter_next(&it)) {
    long l = gkyl_range_idx(&app->local_ext, it.idx);
    if (!route1_fx_masked(l)) continue;
    long s = (long) ((const double *) gkyl_array_cfetch(route1_fx_src, l))[0];
    if (s >= 0) memcpy(gkyl_array_fetch(route1_fx_scr, l), gkyl_array_cfetch(fin, s), esz);
  }
  return route1_fx_scr;
}

// Drop the excised cells' RHS; route1_fx_res = h x (their D-row RHS) x cell volume (the D flux into the mask).
static void
route1_fx_drop(const gkyl_moment_app *app, struct gkyl_array *rhs, double h)
{
  route1_fx_res = 0.0;
  if (!(route1_fx_r > 0.0)) return;
  double dv = 1.0, sum = 0.0;
  for (int d = 0; d < app->ndim; ++d) dv *= app->grid.dx[d];
  struct gkyl_range_iter it;
  gkyl_range_iter_init(&it, &app->local);
  while (gkyl_range_iter_next(&it)) {
    long l = gkyl_range_idx(&app->local_ext, it.idx);
    if (!route1_fx_masked(l)) continue;
    double *q = gkyl_array_fetch(rhs, l);
    sum += q[0];
    for (size_t c = 0; c < rhs->ncomp; ++c) q[c] = 0.0;
  }
  route1_fx_res = h*sum*dv;
}

// Excised cells of f set from fsrc (identical by construction; removes the round-off of the RK3 combination).
static void
route1_fx_keep(const gkyl_moment_app *app, struct gkyl_array *f, const struct gkyl_array *fsrc)
{
  if (!(route1_fx_r > 0.0)) return;
  struct gkyl_range_iter it;
  gkyl_range_iter_init(&it, &app->local);
  while (gkyl_range_iter_next(&it)) {
    long l = gkyl_range_idx(&app->local_ext, it.idx);
    if (route1_fx_masked(l)) memcpy(gkyl_array_fetch(f, l), gkyl_array_cfetch(fsrc, l), f->esznc);
  }
}

// [FXR] Move the newly masked cells' rows 0-4 into the reservoir and set them to the static atmosphere.
static void
route1_fx_remove(gkyl_moment_app *app, struct moment_species *fsp, double t)
{
  if (!route1_fx_new) return;
  const struct wv_gr_euler *ge = container_of(fsp->equation, struct wv_gr_euler, eqn);
  double dv = 1.0, rd = 0.0, rs[3] = { 0.0, 0.0, 0.0 }, rt = 0.0;
  long n = 0;
  for (int d = 0; d < app->ndim; ++d) dv *= app->grid.dx[d];
  struct gkyl_range_iter it;
  gkyl_range_iter_init(&it, &app->local);
  while (gkyl_range_iter_next(&it)) {
    long l = gkyl_range_idx(&app->local_ext, it.idx);
    double *pn = gkyl_array_fetch(route1_fx_new, l);
    if (!(pn[0] > 0.0)) continue;
    pn[0] = 0.0;
    double *q = gkyl_array_fetch(fsp->fcurr, l);
    double g[3][3] = { { q[9], q[10], q[11] }, { q[12], q[13], q[14] }, { q[15], q[16], q[17] } }, gi[3][3];
    double sq = sqrt(fmax(sym3_inverse(g, gi), 0.0));
    double atm[5] = { sq*ge->rho_atm, 0.0, 0.0, 0.0, sq*ge->p_atm/(ge->gas_gamma - 1.0) };
    rd += (q[0] - atm[0])*dv;
    for (int i = 0; i < 3; ++i) rs[i] += q[1 + i]*dv;
    rt += (q[4] - atm[4])*dv;
    for (int c = 0; c < 5; ++c) q[c] = atm[c];
    n += 1;
  }
  if (n == 0) return;
  route1_fx_bh += rd;
  for (int i = 0; i < 3; ++i) route1_fx_sbh[i] += rs[i];
  route1_fx_tbh += rt;
  gkyl_moment_app_cout(app, stderr, "[FX-REMOVE] t=%.6g cells=%ld removed D=%.14e S=%.6e,%.6e,%.6e tau=%.6e;"
    " reservoir D_BH=%.14e S_BH=%.6e,%.6e,%.6e tau_BH=%.6e\n", t, n, rd, rs[0], rs[1], rs[2], rt, route1_fx_bh,
    route1_fx_sbh[0], route1_fx_sbh[1], route1_fx_sbh[2], route1_fx_tbh);
}

static void
route1_fx_report(const gkyl_moment_app *app, const struct gkyl_array *f, double t)
{
  if (!(route1_fx_frac() > 0.0)) return;
  route1_fx_nstep += 1;
  if (route1_fx_nstep % 10 != 0) return;
  double dv = 1.0, live = 0.0, mask = 0.0;
  for (int d = 0; d < app->ndim; ++d) dv *= app->grid.dx[d];
  struct gkyl_range_iter it;
  gkyl_range_iter_init(&it, &app->local);
  while (gkyl_range_iter_next(&it)) {
    long l = gkyl_range_idx(&app->local_ext, it.idx);
    double D = ((const double *) gkyl_array_cfetch(f, l))[0]*dv;
    if (route1_fx_masked(l)) mask += D; else live += D;
  }
  gkyl_moment_app_cout(app, stderr, "[FX] t=%.6g r_mask=%.4f n=%ld D_live=%.14e D_mask=%.14e D_BH=%.14e"
    " sum=%.14e lambda_max=%+.3e illegal=%ld\n", t, route1_fx_r, route1_fx_n, live, mask, route1_fx_bh,
    live + mask + route1_fx_bh, route1_fx_lam, route1_fx_illegal);
}

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static bool
route1_enabled(void)
{
  const char *v = getenv("GKYL_ROUTE1");
  return v && atoi(v) != 0;
}

static double
route1_sigma(void)
{
  const char *v = getenv("GKYL_ROUTE1_KO");
  return v ? atof(v) : 0.3;
}

static const double route1_frames[3][3][3] = { // tau1, tau2, norm (as identity_test.c / split.c)
  { {0,1,0}, {0,0,1}, {1,0,0} },
  { {0,0,1}, {1,0,0}, {0,1,0} },
  { {1,0,0}, {0,1,0}, {0,0,1} } };

static void
route1_flux(const struct moment_species *esp, int d, const double *q, double *f)
{
  double ql[77], fl[77];
  gkyl_wv_eqn_rotate_to_local(esp->equation, route1_frames[d][0], route1_frames[d][1], route1_frames[d][2], q, ql);
  gkyl_vacuum_einstein_conformal_flux(esp->vacuum_einstein_conformal_excision_threshold,
    esp->vacuum_einstein_conformal_spacetime_slicing, esp->vacuum_einstein_conformal_spacetime_evolution, ql, fl);
  gkyl_wv_eqn_rotate_to_global(esp->equation, route1_frames[d][0], route1_frames[d][1], route1_frames[d][2], fl, f);
}

// Gammatilde^i reconstructed from D_kij, exactly as gamma_driver_shift_solve.
static void
route1_gamma(const double *q, double Gam[3])
{
  double gt[3][3], gti[3][3];
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      gt[i][j] = q[3*i + j];
  sym3_inverse(gt, gti);
  for (int i = 0; i < 3; ++i) {
    double acc = 0.0;
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 3; ++k) {
        double gjk = gti[j][k];
        if (gjk == 0.0) continue;
        for (int l = 0; l < 3; ++l)
          acc += gjk * gti[i][l] * (q[19 + 9*j + 3*l + k] + q[19 + 9*k + 3*l + j] - q[19 + 9*l + 3*j + k]);
      }
    Gam[i] = isfinite(acc) ? acc : 0.0;
  }
}

static double
route1_spure(const struct moment_species *esp, const double *q)
{
  double s[77];
  gkyl_wv_eqn_source(esp->equation, q, s);
  return s[64] + q[52]*q[65] + q[53]*q[66] + q[54]*q[67] + 2*q[64]*(q[55] + q[59] + q[63]);
}

static long
route1_nb(const gkyl_moment_app *app, const int *idx, int d, int off)
{
  int j[GKYL_MAX_DIM];
  for (int k = 0; k < app->ndim; ++k) j[k] = idx[k];
  j[d] += off;
  return gkyl_range_idx(&app->local_ext, j);
}

static bool
route1_in_local(const gkyl_moment_app *app, const int *idx, int d, int off)
{
  int v = idx[d] + off;
  return v >= app->local.lower[d] && v <= app->local.upper[d];
}

// GKYL_ROUTE1_KO_GHOSTS=1: KO stencils may reach into the BC-filled ghost layers, so the face flux of
// the dissipation also acts on the boundary cells (as ETK applies dissipation up to the boundary). The KO
// flux through a boundary face is then one more boundary flux: interior sums still change only through
// boundary fluxes. Default 0: G = 0 on any face whose stencil leaves the interior.
static bool
route1_ko_ghosts(void)
{
  const char *v = getenv("GKYL_ROUTE1_KO_GHOSTS");
  return v && atoi(v) != 0;
}

// GKYL_ROUTE1_KO_SBP=1: boundary-modified dissipation of Mattsson, Svard & Nordstrom (J. Sci. Comput. 21,
// 2004) in place of the cut KO. Per direction the rate is -sigma/(16h) D^T D w on the interior cells, where D
// is the undivided second difference s_c = w_{c-1} - 2 w_c + w_{c+1}, with the boundary cell's row taken
// one-sided (equal to its neighbour's). In face-flux form G_{k+1/2} = m_{k+1} s_{k+1} - m_k s_k, with m = 1
// inside, m = 2 on the second cell from each boundary (it also carries the boundary cell's row) and m = 0 on
// the boundary cells, so G = 0 on the boundary faces. This equals the cut KO away from the boundary, is
// conservative, is negative semidefinite (the cut KO is not: zeroing faces leaves an indefinite boundary term
// in its energy), reaches every cell including edges and corners, and reads no ghost cell.
static bool
route1_ko_sbp(void)
{
  const char *v = getenv("GKYL_ROUTE1_KO_SBP");
  return v && atoi(v) != 0;
}

static double
route1_ko_sbp_weight(const gkyl_moment_app *app, int d, int c)
{
  int lo = app->local.lower[d], up = app->local.upper[d];
  if (c <= lo || c >= up) return 0.0;
  return (c == lo + 1 || c == up - 1) ? 2.0 : 1.0;
}

static bool
route1_in_ext(const gkyl_moment_app *app, const int *idx, int d, int off)
{
  int v = idx[d] + off;
  return v >= app->local_ext.lower[d] && v <= app->local_ext.upper[d];
}

static bool
route1_global_periodic(const gkyl_moment_app *app, int d)
{
  for (int k = 0; k < app->num_periodic_dir; ++k)
    if (app->periodic_dirs[k] == d) return true;
  return false;
}

// GKYL_ROUTE1_B_CENTRED=1 (NewRad only): B = 1/2 D0 (beta rate) with the centred stencil on the boundary
// cells too, reading the radiated beta rate of the ghosts, so one stencil ties B to beta across the faces (as
// P and H already are). Default 0: one-sided inward at a non-periodic global face (as production).
static bool
route1_b_centred(void)
{
  const char *v = getenv("GKYL_ROUTE1_B_CENTRED");
  return v && atoi(v) != 0;
}

// Periodic directions: copy the ghosts of a rate-like array (rate, lp) from across the seam, so stencils that
// read them see the same values as anywhere else. Fills edge and corner ghosts too. No-op without periodic
// directions.
static void
route1_periodic_sync(const gkyl_moment_app *app, struct gkyl_array *a)
{
  if (app->num_periodic_dir == 0) return;
  struct gkyl_range_iter it;
  gkyl_range_iter_init(&it, &app->local_ext);
  while (gkyl_range_iter_next(&it)) {
    int j[GKYL_MAX_DIM];
    bool wrap = false;
    for (int d = 0; d < app->ndim; ++d) {
      j[d] = it.idx[d];
      if (!route1_global_periodic(app, d)) continue;
      int n = app->local.upper[d] - app->local.lower[d] + 1;
      if (j[d] < app->local.lower[d]) { j[d] += n; wrap = true; }
      else if (j[d] > app->local.upper[d]) { j[d] -= n; wrap = true; }
    }
    if (wrap)
      memcpy(gkyl_array_fetch(a, gkyl_range_idx(&app->local_ext, it.idx)),
        gkyl_array_cfetch(a, gkyl_range_idx(&app->local_ext, j)), a->esznc);
  }
}

// GKYL_ROUTE1_DRIVER_BUFFER=W (> 0): gauge buffer of width W before the outer faces (a slab along each
// non-periodic global face, by distance to the nearest one). Across it the driver's mu goes smoothly (C2) to
// GKYL_ROUTE1_DRIVER_BUFFER_MU times its value at the face (default 0), and eta to
// GKYL_ROUTE1_DRIVER_BUFFER_ETA (default: unchanged). B = 1/2 D0 (beta rate) carries the gradient of the
// profile, so the reduction identity is untouched. Default off: mu and eta are returned unchanged.
static void
route1_driver_buffer(const gkyl_moment_app *app, const int *idx, double *mu, double *eta)
{
  static int init = 0;
  static double width = 0.0, mu_face = 0.0, eta_face = -1.0;
  if (!init) {
    const char *v = getenv("GKYL_ROUTE1_DRIVER_BUFFER");
    width = v ? atof(v) : 0.0;
    v = getenv("GKYL_ROUTE1_DRIVER_BUFFER_MU");
    mu_face = v ? atof(v) : 0.0;
    v = getenv("GKYL_ROUTE1_DRIVER_BUFFER_ETA");
    eta_face = v ? atof(v) : -1.0;
    init = 1;
  }
  if (width <= 0) return;
  double x[3] = { 0.0, 0.0, 0.0 }, dist = 1e300;
  gkyl_rect_grid_cell_center(&app->grid, idx, x);
  for (int d = 0; d < app->ndim; ++d)
    if (!route1_global_periodic(app, d))
      dist = fmin(dist, fmin(x[d] - app->grid.lower[d], app->grid.upper[d] - x[d]));
  double s = 1.0 - dist/width; // 0 at the inner edge of the buffer, 1 at the face
  if (s <= 0) return;
  double b = s >= 1 ? 1.0 : s*s*s*(10 - 15*s + 6*s*s);
  if (eta_face >= 0) *eta += b*(eta_face - *eta);
  *mu *= 1 - b*(1 - mu_face);
}

// Depth in cells below the nearest non-periodic global face: 0 on the outermost interior cell, negative on
// ghosts, large when every direction is periodic.
static int
route1_depth(const gkyl_moment_app *app, const int *idx)
{
  int dep = 1 << 30;
  for (int d = 0; d < app->ndim; ++d) {
    if (route1_global_periodic(app, d)) continue;
    int a = idx[d] - app->global.lower[d], b = app->global.upper[d] - idx[d];
    if (a < dep) dep = a;
    if (b < dep) dep = b;
  }
  return dep;
}

// GKYL_ROUTE1_BAND=W (> 0): upwind band of W cells before the non-periodic outer faces (see route1_rhs).
static int
route1_band_width(void)
{
  const char *v = getenv("GKYL_ROUTE1_BAND");
  return v ? atoi(v) : 0;
}

// Ramp of the band: 1 on the boundary face (outer cell depth -1), cos^2 down to 0 at depth W.
static double
route1_band_weight(int d, int W)
{
  if (d < 0) return 1.0;
  if (d >= W) return 0.0;
  double c = cos(0.5*M_PI*(d + 1.0)/(W + 1.0));
  return c*c;
}

// GKYL_ROUTE1_BAND_DERIVED=1: the band's face terms act on the parent rows only (gt, ln alpha, K, beta, chi);
// D, A and V take the increments implied by their parents (route1_band_derive), as B follows beta and P
// follows chi. Default 0: every evolved row gets its own face term.
static bool
route1_band_derived(void)
{
  const char *v = getenv("GKYL_ROUTE1_BAND_DERIVED");
  return v && atoi(v) != 0;
}

// The band's one global speed: GKYL_ROUTE1_BAND_S, or the largest characteristic speed over the box and the
// three directions in the state of the first call. Fixed for the run, so the band stays linear.
static double
route1_band_speed(gkyl_moment_app *app, const struct moment_species *esp, const struct gkyl_array *u)
{
  static double s = -1.0;
  if (s >= 0) return s;
  const char *v = getenv("GKYL_ROUTE1_BAND_S");
  if (v) s = atof(v);
  else {
    s = 0.0;
    struct gkyl_range_iter it;
    gkyl_range_iter_init(&it, &app->local);
    while (gkyl_range_iter_next(&it)) {
      const double *q = gkyl_array_cfetch(u, gkyl_range_idx(&app->local_ext, it.idx));
      for (int d = 0; d < app->ndim; ++d) {
        double ql[77];
        gkyl_wv_eqn_rotate_to_local(esp->equation, route1_frames[d][0], route1_frames[d][1], route1_frames[d][2], q, ql);
        double c = gkyl_wv_eqn_max_speed(esp->equation, ql);
        if (isfinite(c) && c > s) s = c;
      }
    }
  }
  gkyl_moment_app_cout(app, stderr, "[ROUTE1] upwind band: W=%d cells, s=%g (%s)\n", route1_band_width(), s,
    v ? "GKYL_ROUTE1_BAND_S" : "largest characteristic speed at the first step");
  if (route1_band_derived())
    gkyl_moment_app_cout(app, stderr, "[ROUTE1] band derived: face terms on gt, ln alpha, K, beta, chi; D, A, V"
      " follow (GKYL_ROUTE1_BAND_DERIVED)\n");
  return s;
}

struct route1_work {
  struct gkyl_array *L, *G, *sp, *lp, *u1, *u2, *rate;
  struct gkyl_array *band;  // GKYL_ROUTE1_BAND: the band's contribution to the rate
  struct gkyl_array *bandv; // GKYL_ROUTE1_BAND_DERIVED: gt^kl R_kl, gt^kl R_ki (route1_band_derive)
  struct gkyl_array *idb0;  // GKYL_ROUTE1_IDB: C_D 27, C_A 3, C_V 3 on every local cell at route-1 start
  struct gkyl_array *Fe[3], *spe; // F_d(e), S_pure(e): constant, computed once
  struct gkyl_array *id0;   // reduction residuals at route-1 start (C_D 27, C_A 3, C_P 3, C_B 9, C_H 9)
  struct gkyl_array *gh;    // GKYL_ROUTE1_NEWRAD: evolved ghost states, kept between steps
  long steps;
  bool lalpha_checked;
};
static struct route1_work route1_w;

#include "route1_newrad.h"
#include "route1_charbc.h"
#include "route1_symbol.h"

// Band increments of D, A and V implied by the parent band rates R in band (gt rows 0..8, ln alpha row 9):
//   D_kij: 1/2 D0_k R_ij        A_k: D0_k R_9
//   V_i:   1/2 D0_i (gt^kl R_kl) - 1/2 D0_l (gt^kl R_ki)     (Vtilde_i = gt^kl (D_ikl - D_lki))
// On a boundary cell D0 reads the ghost's share of the band. The NewRad ghost rate carries (r_p/r)^N times
// the actual rate of its inward cell p (in alpha, not ln alpha), so the band moves the ghost's gt by
// (r_p/r)^N R_ij(p) and its ln alpha by (r_p/r)^N (alpha_p/alpha) R_9(p); without NewRad the share is 0.
// Each increment is a centred difference of one value per cell, i.e. a difference of single-valued face
// values, so interior sums change only through the boundary faces. C_D and C_A are kept exactly, boundary
// cells included. C_V is kept up to the terms R d(gt^kl) and (d_t gt^kl)(D_ikl - D_lki), second order in the
// deviation, which the flux form leaves out. aux (10 components per cell) holds gt^kl R_kl and gt^kl R_ki
// (index 1 + 3l + i), on the band cells and the face ghosts.
static void
route1_band_derive(gkyl_moment_app *app, const struct gkyl_array *u, struct gkyl_array *band,
  struct gkyl_array *aux, int band_w)
{
  int ndim = app->ndim;
  bool nr = route1_newrad_enabled();
  double npow = route1_newrad_power();
  struct gkyl_range_iter it;
  gkyl_array_clear(aux, 0.0);
  gkyl_range_iter_init(&it, &app->local_ext);
  while (gkyl_range_iter_next(&it)) {
    int s[3];
    int layer = route1_ghost_dir(app, it.idx, s);
    long l = gkyl_range_idx(&app->local_ext, it.idx);
    const double *q = gkyl_array_cfetch(u, l);
    double *R = gkyl_array_fetch(band, l);
    if (layer == 0) {
      if (route1_depth(app, it.idx) > band_w) continue;
    }
    else {
      // face ghosts (one boundary direction) are the only ghosts a centred D0 on local reads
      if (!nr || layer != 1 || abs(s[0]) + abs(s[1]) + abs(s[2]) != 1) continue;
      int ip[GKYL_MAX_DIM];
      for (int d = 0; d < ndim; ++d) ip[d] = it.idx[d] - s[d];
      long lp = gkyl_range_idx(&app->local_ext, ip);
      const double *Rp = gkyl_array_cfetch(band, lp), *qp = gkyl_array_cfetch(u, lp);
      double scale = npow >= 0 ? pow(route1_radius(app, ip)/route1_radius(app, it.idx), npow) : 0.0;
      for (int m = 0; m < 9; ++m) R[m] = scale*Rp[m];
      R[9] = scale*qp[9]*Rp[9]/q[9];
    }
    double *a = gkyl_array_fetch(aux, l), gt[3][3], gti[3][3];
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        gt[i][j] = q[3*i + j];
    sym3_inverse(gt, gti);
    for (int k = 0; k < 3; ++k)
      for (int m = 0; m < 3; ++m) {
        a[0] += gti[k][m]*R[3*k + m];
        for (int i = 0; i < 3; ++i) a[1 + 3*m + i] += gti[k][m]*R[3*k + i];
      }
  }
  gkyl_range_iter_init(&it, &app->local);
  while (gkyl_range_iter_next(&it)) {
    if (route1_depth(app, it.idx) > band_w + 1) continue;
    double *bd = gkyl_array_fetch(band, gkyl_range_idx(&app->local_ext, it.idx));
    for (int m = 19; m < 52; ++m) bd[m] = 0.0;
    for (int k = 0; k < ndim; ++k) {
      double h = app->grid.dx[k];
      const double *rp = gkyl_array_cfetch(band, route1_nb(app, it.idx, k, 1));
      const double *rm = gkyl_array_cfetch(band, route1_nb(app, it.idx, k, -1));
      const double *ap = gkyl_array_cfetch(aux, route1_nb(app, it.idx, k, 1));
      const double *am = gkyl_array_cfetch(aux, route1_nb(app, it.idx, k, -1));
      for (int ij = 0; ij < 9; ++ij) bd[19 + 9*k + ij] = 0.5*(rp[ij] - rm[ij])/(2*h);
      bd[46 + k] = (rp[9] - rm[9])/(2*h);
      bd[49 + k] += 0.5*(ap[0] - am[0])/(2*h);
      for (int i = 0; i < 3; ++i) bd[49 + i] -= 0.5*(ap[1 + 3*k + i] - am[1 + 3*k + i])/(2*h);
    }
  }
}

// Reduction residuals of state u at interior cell idx (central D0, ln alpha).
static void
route1_residuals(const gkyl_moment_app *app, const struct gkyl_array *u, const int *idx, double *c)
{
  const double *q = gkyl_array_cfetch(u, gkyl_range_idx(&app->local_ext, idx));
  int n = 0;
  for (int k = 0; k < 3; ++k) {
    const double *p = gkyl_array_cfetch(u, route1_nb(app, idx, k, 1)), *m = gkyl_array_cfetch(u, route1_nb(app, idx, k, -1));
    double h = app->grid.dx[k];
    for (int ij = 0; ij < 9; ++ij) c[n++] = q[19+9*k+ij] - 0.5*(p[ij]-m[ij])/(2*h);
  }
  for (int k = 0; k < 3; ++k) {
    const double *p = gkyl_array_cfetch(u, route1_nb(app, idx, k, 1)), *m = gkyl_array_cfetch(u, route1_nb(app, idx, k, -1));
    double h = app->grid.dx[k];
    c[n++] = q[46+k] - (log(p[9])-log(m[9]))/(2*h);
    c[n++] = q[65+k] - (p[64]-m[64])/(2*h);
    for (int j = 0; j < 3; ++j) c[n++] = q[55+3*k+j] - 0.5*(p[52+j]-m[52+j])/(2*h);
    for (int j = 0; j < 3; ++j) c[n++] = q[68+3*k+j] - 0.5*((p[65+j]-m[65+j])/(2*h));
  }
  // H: complete the symmetric part with D0_j P_k
  n = 27;
  for (int k = 0; k < 3; ++k) {
    n += 2 + 3;
    for (int j = 0; j < 3; ++j) {
      const double *p = gkyl_array_cfetch(u, route1_nb(app, idx, j, 1)), *m = gkyl_array_cfetch(u, route1_nb(app, idx, j, -1));
      c[n++] -= 0.5*(p[65+k]-m[65+k])/(2*app->grid.dx[j]);
    }
  }
}
#define ROUTE1_NRES 51

static bool
route1_check_interior(const gkyl_moment_app *app, const int *idx)
{
  for (int d = 0; d < app->ndim; ++d)
    if (idx[d] - app->local.lower[d] < 3 || app->local.upper[d] - idx[d] < 3) return false;
  return true;
}

static void
route1_report(gkyl_moment_app *app, const struct gkyl_array *u, double t)
{
  const char *names[5] = { "C_D", "C_A", "C_P", "C_B", "C_H" };
  double worst[5] = { 0 }, c[ROUTE1_NRES];
  bool init = route1_w.id0 == NULL;
  if (init) route1_w.id0 = gkyl_array_new(GKYL_DOUBLE, ROUTE1_NRES, app->local_ext.volume);
  struct gkyl_range_iter it;
  gkyl_range_iter_init(&it, &app->local);
  while (gkyl_range_iter_next(&it)) {
    if (!route1_check_interior(app, it.idx)) continue;
    double *c0 = gkyl_array_fetch(route1_w.id0, gkyl_range_idx(&app->local_ext, it.idx));
    route1_residuals(app, u, it.idx, c);
    if (init) { memcpy(c0, c, sizeof c); continue; }
    for (int m = 0; m < ROUTE1_NRES; ++m) {
      int g = m < 27 ? 0 : ((m-27) % 8 == 0 ? 1 : (m-27) % 8 == 1 ? 2 : (m-27) % 8 < 5 ? 3 : 4);
      worst[g] = fmax(worst[g], fabs(c[m] - c0[m]));
    }
  }
  if (init) {
    gkyl_moment_app_cout(app, stderr, "[ROUTE1] t=%.6g sigma=%g ko_ghosts=%d: reduction residuals captured\n", t,
      route1_sigma(), route1_ko_ghosts() || route1_peel());
    if (route1_peel())
      gkyl_moment_app_cout(app, stderr, "[ROUTE1] peel: NewRad on gt, ln alpha, K, V, beta, chi; D, A, P, B, H derived"
        " on the ghosts after every stage; KO to the boundary in every direction; B centred\n");
    for (int s = 0; s < app->num_species && route1_newrad_enabled(); ++s)
      if (app->species[s].has_vacuum_einstein_conformal)
        gkyl_moment_app_cout(app, stderr, "[ROUTE1] newrad ghosts: N=%g v_g=%g oldspeed=%d dpow=%g\n", route1_newrad_power(),
          route1_newrad_vgauge(app->species[s].vacuum_einstein_conformal_spacetime_slicing), route1_newrad_oldspeed(),
          route1_newrad_dpow());
    if (route1_newrad_enabled() && (route1_b_centred() || route1_peel() || route1_newrad_beta_mu() > 0))
      gkyl_moment_app_cout(app, stderr, "[ROUTE1] shift boundary: B centred on boundary cells=%d, beta split mu=%g\n",
        route1_b_centred() || route1_peel(), route1_newrad_beta_mu());
    if (getenv("GKYL_ROUTE1_DRIVER_BUFFER"))
      gkyl_moment_app_cout(app, stderr, "[ROUTE1] driver buffer: width=%s mu_face=%s eta_face=%s\n",
        getenv("GKYL_ROUTE1_DRIVER_BUFFER"),
        getenv("GKYL_ROUTE1_DRIVER_BUFFER_MU") ? getenv("GKYL_ROUTE1_DRIVER_BUFFER_MU") : "0",
        getenv("GKYL_ROUTE1_DRIVER_BUFFER_ETA") ? getenv("GKYL_ROUTE1_DRIVER_BUFFER_ETA") : "unchanged");
    return;
  }
  gkyl_moment_app_cout(app, stderr, "[ROUTE1-ID] t=%.6g step=%ld max|C-C0| (interior, 3 cells in):", t, route1_w.steps);
  for (int g = 0; g < 5; ++g) gkyl_moment_app_cout(app, stderr, " %s=%.2e", names[g], worst[g]);
  gkyl_moment_app_cout(app, stderr, "\n");
}

static bool
route1_idb_enabled(void)
{
  const char *v = getenv("GKYL_ROUTE1_IDB");
  return v && atoi(v) != 0;
}

// GKYL_ROUTE1_IDB=1 (with the monitor): one line [ROUTE1-IDB] per monitor line with max|C - C0| of C_D (27,
// D_kij - 1/2 D0_k gt_ij), C_A (3, A_k - D0_k ln alpha) and C_V (3, V_i - gt^kl (D_ikl - D_lki)) by depth of
// the cell below the nearest non-periodic face: 0, 1, 2, 3, 4-5, >= 6. D0 on the boundary cells reads the
// ghosts. The first call (valid ghosts at route-1 start) captures C0. Diagnostic; default off.
static void
route1_idb(gkyl_moment_app *app, const struct gkyl_array *u, double t)
{
  bool init = route1_w.idb0 == NULL;
  if (init) route1_w.idb0 = gkyl_array_new(GKYL_DOUBLE, 33, app->local_ext.volume);
  double worst[3][6] = { { 0 } };
  struct gkyl_range_iter it;
  gkyl_range_iter_init(&it, &app->local);
  while (gkyl_range_iter_next(&it)) {
    long l = gkyl_range_idx(&app->local_ext, it.idx);
    const double *q = gkyl_array_cfetch(u, l);
    double c[33] = { 0 }, gt[3][3], gti[3][3];
    for (int k = 0; k < app->ndim; ++k) {
      const double *p = gkyl_array_cfetch(u, route1_nb(app, it.idx, k, 1)), *m = gkyl_array_cfetch(u, route1_nb(app, it.idx, k, -1));
      double h = app->grid.dx[k];
      for (int ij = 0; ij < 9; ++ij) c[9*k + ij] = q[19 + 9*k + ij] - 0.5*(p[ij] - m[ij])/(2*h);
      c[27 + k] = q[46 + k] - (log(p[9]) - log(m[9]))/(2*h);
    }
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        gt[i][j] = q[3*i + j];
    sym3_inverse(gt, gti);
    for (int i = 0; i < 3; ++i) {
      c[30 + i] = q[49 + i];
      for (int k = 0; k < 3; ++k)
        for (int m = 0; m < 3; ++m)
          c[30 + i] -= gti[k][m]*(q[19 + 9*i + 3*k + m] - q[19 + 9*m + 3*k + i]);
    }
    double *c0 = gkyl_array_fetch(route1_w.idb0, l);
    if (init) { memcpy(c0, c, sizeof c); continue; }
    int d = route1_depth(app, it.idx), b = d <= 3 ? d : (d <= 5 ? 4 : 5);
    for (int m = 0; m < 33; ++m) {
      int g = m < 27 ? 0 : (m < 30 ? 1 : 2);
      worst[g][b] = fmax(worst[g][b], fabs(c[m] - c0[m]));
    }
  }
  if (init) {
    gkyl_moment_app_cout(app, stderr, "[ROUTE1] t=%.6g identities by depth (GKYL_ROUTE1_IDB): C0 captured\n", t);
    return;
  }
  const char *names[3] = { "C_D", "C_A", "C_V" };
  gkyl_moment_app_cout(app, stderr, "[ROUTE1-IDB] t=%.6g max|C-C0| at depth 0,1,2,3,4-5,>=6:", t);
  for (int g = 0; g < 3; ++g)
    gkyl_moment_app_cout(app, stderr, " %s=%.2e,%.2e,%.2e,%.2e,%.2e,%.2e", names[g], worst[g][0], worst[g][1],
      worst[g][2], worst[g][3], worst[g][4], worst[g][5]);
  gkyl_moment_app_cout(app, stderr, "\n");
}

// GKYL_ROUTE1_MONITOR=T (> 0): every T time units one line [ROUTE1-MON] with max|alpha - 1|, max|trK|,
// max|beta|, the maximum over the 27 D and over the 3 A components (q - e; flat values without e), over the
// outer 3 layers (depth < 3) and over the box. Values of this rank only. Diagnostic; default off.
static void
route1_monitor(gkyl_moment_app *app, const struct gkyl_array *u, const struct gkyl_array *e, double t)
{
  static long last = -1;
  const char *ev = getenv("GKYL_ROUTE1_MONITOR");
  double every = ev ? atof(ev) : 0.0;
  if (every <= 0) return;
  long slot = (long) floor(t/every + 1e-6);
  if (slot == last) return;
  last = slot;
  double v[2][5] = { { 0 } };
  struct gkyl_range_iter it;
  gkyl_range_iter_init(&it, &app->local);
  while (gkyl_range_iter_next(&it)) {
    long l = gkyl_range_idx(&app->local_ext, it.idx);
    const double *q = gkyl_array_cfetch(u, l), *qe = e ? gkyl_array_cfetch(e, l) : NULL;
    double c[5] = { 0 }, b2 = 0.0;
    c[0] = fabs(q[9] - (qe ? qe[9] : 1.0));
    c[1] = fabs(route1_trk(q) - (qe ? route1_trk(qe) : 0.0));
    for (int j = 0; j < 3; ++j) {
      double d = q[52 + j] - (qe ? qe[52 + j] : 0.0);
      b2 += d*d;
    }
    c[2] = sqrt(b2);
    for (int m = 19; m < 46; ++m) c[3] = fmax(c[3], fabs(q[m] - (qe ? qe[m] : 0.0)));
    for (int m = 46; m < 49; ++m) c[4] = fmax(c[4], fabs(q[m] - (qe ? qe[m] : 0.0)));
    bool outer = route1_depth(app, it.idx) < 3;
    for (int k = 0; k < 5; ++k) {
      v[1][k] = fmax(v[1][k], c[k]);
      if (outer) v[0][k] = fmax(v[0][k], c[k]);
    }
  }
  gkyl_moment_app_cout(app, stderr, "[ROUTE1-MON] t=%.6g outer3 a=%.4e k=%.4e b=%.4e D=%.4e A=%.4e"
    " box a=%.4e k=%.4e b=%.4e D=%.4e A=%.4e\n", t, v[0][0], v[0][1], v[0][2], v[0][3], v[0][4],
    v[1][0], v[1][1], v[1][2], v[1][3], v[1][4]);
  if (route1_idb_enabled() && route1_w.idb0) route1_idb(app, u, t);
}

// GKYL_ROUTE1_DRIVER_METRIC=1: the driver's Gammatilde^i from the parent metric: D_kij replaced by 1/2 D_k gt_ij
// with route1_dr's stencil (centred on every interior cell; the boundary cells read the ghosts), and the reference
// captured the same way at activation. It equals production's Gammatilde when C_D = 0, so the continuum gauge is
// unchanged; the reduction-constraint error no longer enters the shift. Default 0: from the evolved D (production).
static bool
route1_driver_metric(void)
{
  const char *v = getenv("GKYL_ROUTE1_DRIVER_METRIC");
  return v && atoi(v) != 0;
}

static void
route1_gamma_metric(const gkyl_moment_app *app, const struct gkyl_array *u, const int *idx, double Gam[3])
{
  double q[46] = { 0.0 };   // route1_gamma reads gt (0-8) and D (19-45)
  const double *q0 = gkyl_array_cfetch(u, gkyl_range_idx(&app->local_ext, idx));
  for (int m = 0; m < 9; ++m) q[m] = q0[m];
  for (int k = 0; k < app->ndim; ++k)
    for (int ij = 0; ij < 9; ++ij) q[19 + 9*k + ij] = 0.5*route1_dr(app, u, ij, idx, k);
  route1_gamma(q, Gam);
}

// Gammatilde part of the driver reference on the cells the driver acts on, with route1_gamma_metric (u with valid
// ghosts); the beta part stays as captured.
static void
route1_driver_metric_reference(gkyl_moment_app *app, struct moment_species *esp, const struct gkyl_array *u)
{
  double alpha_floor = gamma_driver_alpha_floor(esp);
  struct gkyl_range_iter it;
  gkyl_range_iter_init(&it, &app->local);
  while (gkyl_range_iter_next(&it)) {
    long l = gkyl_range_idx(&app->local_ext, it.idx);
    if (gamma_driver_cell_active(gkyl_array_cfetch(u, l), alpha_floor, true))
      route1_gamma_metric(app, u, it.idx, gkyl_array_fetch(esp->gamma_driver_reference, l));
  }
}

// [ZD] Damping of Z_i = -C_V_i in the V rows (this build only; GKYL_ROUTE1_ZDAMP=kappa_1, unset = off).
// d_t V_i gets -kappa_1 (C_V_i(u) - C_V_i(e)), with C_V_i = V_i - gt^kl (D_ikl - D_lki) (the route1_idb identity) and
// e the equilibrium when the WB step is active (the term then vanishes on it). This is the Z_i damping of Z4
// (Gundlach et al., gr-qc/0504114: d_t Z_i contains -kappa_1 alpha Z_i) with the black-hole prescription
// kappa_1 -> kappa_1 / alpha of Alic, Kastaun & Rezzolla, PRD 88, 064049 (2013), Eq. (27): the product is kappa_1, so
// the damping does not vanish where the lapse collapses. With Theta = 0 this is the only damping term.
static double
route1_zdamp(void)
{
  static double k = -1.0;
  if (k < 0.0) {
    const char *v = getenv("GKYL_ROUTE1_ZDAMP");
    k = v ? atof(v) : 0.0;
    if (!(k > 0.0)) k = 0.0;
    if (k > 0.0) fprintf(stderr, "[ZD] d_t V_i -= %g (C_V_i(u) - C_V_i(e)) (kappa_1/alpha prescription)\n", k);
  }
  return k;
}

static void
route1_cv(const double *q, double cv[3])
{
  double gt[3][3], gti[3][3];
  for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) gt[i][j] = q[3*i + j];
  sym3_inverse(gt, gti);
  for (int i = 0; i < 3; ++i) {
    cv[i] = q[49 + i];
    for (int k = 0; k < 3; ++k)
      for (int m = 0; m < 3; ++m)
        cv[i] -= gti[k][m]*(q[19 + 9*i + 3*k + m] - q[19 + 9*m + 3*k + i]);
  }
}

// [FC] Face-flux ordering-constraint term (this build only; GKYL_ROUTE1_FACECONS=c, unset = off), the finite-volume
// analogue of the gamma_2 term of first-order reductions (Lindblom et al., CQG 23, S447 (2006); Owen, PRD 76, 044019
// (2007)). For a primary p with derivative variable d_k = (1/f) d_k p, the compact face constraint at i+1/2 along k is
//   C_{i+1/2} = (p_{i+1} - p_i)/dx - f (d_{k,i} + d_{k,i+1})/2,
// which sees the two-cell modes of p that the centred identity d_k - D0_k p / f cannot (D0 is blind to them). One
// flux per face, used by both cells: G_{i+1/2} = nu C_{i+1/2}, nu = c dx, rate(p) += (G_{i+1/2} - G_{i-1/2})/dx, on
// w = u - e as KO, no alpha factor, G = 0 on a face with a cell outside local (non-periodic). Pairs (p, d, f):
// ln alpha / A_k (1), gt_ij / D_kij (2), beta^j / B_kj (2). GKYL_ROUTE1_FACECONS_FIELDS=a: ln alpha only;
// default agb (all three). Spectral bound of the term: 4c/dx per direction (12c/dx in 3D).
static double
route1_fc(int *fields)
{
  static double c = -1.0;
  static int f = 0;
  if (c < 0.0) {
    const char *v = getenv("GKYL_ROUTE1_FACECONS");
    c = v ? atof(v) : 0.0;
    if (!(c > 0.0)) c = 0.0;
    const char *fs = getenv("GKYL_ROUTE1_FACECONS_FIELDS");
    f = fs && strcmp(fs, "a") == 0 ? 1 : 7; // 1 ln alpha, 2 gt, 4 beta
    if (c > 0.0) fprintf(stderr, "[FC] face-flux constraint term nu = %g dx on%s%s%s\n", c,
      f & 1 ? " ln alpha/A" : "", f & 2 ? " gt/D" : "", f & 4 ? " beta/B" : "");
  }
  if (fields) *fields = f;
  return c;
}

// rate = L(u) - L(e). u must have valid ghosts. Slot 9 of rate is d(ln alpha)/dt.
static void
route1_rhs(gkyl_moment_app *app, struct moment_species *esp, const struct gkyl_array *u,
  const struct gkyl_array *e, const struct gkyl_array *seq, const struct gkyl_array *fluid_f,
  double gas_gamma, double t, struct gkyl_array *rate)
{
  struct route1_work *w = &route1_w;
  int ndim = app->ndim, nc = esp->num_equations;
  double sigma = route1_sigma();
  enum gkyl_spacetime_evolution evolution = esp->vacuum_einstein_conformal_spacetime_evolution;
  gkyl_array_clear(rate, 0.0);
  struct gkyl_range_iter it;

  // 1. Central flux divergence of F(u) - F(e) on rows 0..51 (alpha carries no flux; checked once).
  //    GKYL_ROUTE1_CHARBC: done after step 3 by route1_charbc_divergence, with B and P rows and the boundary flux.
  bool cb = route1_charbc();
  for (int d = 0; d < ndim && !cb; ++d) {
    gkyl_range_iter_init(&it, &app->local_ext);
    while (gkyl_range_iter_next(&it)) {
      long l = gkyl_range_idx(&app->local_ext, it.idx);
      double *g = gkyl_array_fetch(w->G, l);
      route1_flux(esp, d, gkyl_array_cfetch(u, l), g);
      if (e) {
        if (!w->Fe[d]) w->Fe[d] = gkyl_array_new(GKYL_DOUBLE, nc, app->local_ext.volume);
        double *fe = gkyl_array_fetch(w->Fe[d], l);
        if (!w->spe) route1_flux(esp, d, gkyl_array_cfetch(e, l), fe);
        for (int m = 0; m < nc; ++m) g[m] -= fe[m];
      }
      if (!w->lalpha_checked && fabs(g[9]) > 0) {
        gkyl_moment_app_cout(app, stderr, "[ROUTE1] WARNING alpha row carries flux %.3e; ln-alpha pairing assumes none\n", g[9]);
        w->lalpha_checked = true;
      }
    }
    double h = app->grid.dx[d];
    gkyl_range_iter_init(&it, &app->local);
    while (gkyl_range_iter_next(&it)) {
      double *r = gkyl_array_fetch(rate, gkyl_range_idx(&app->local_ext, it.idx));
      const double *gp = gkyl_array_cfetch(w->G, route1_nb(app, it.idx, d, 1));
      const double *gm = gkyl_array_cfetch(w->G, route1_nb(app, it.idx, d, -1));
      for (int m = 0; m < 52; ++m) if (m != 9) r[m] -= (gp[m] - gm[m])/(2*h);
    }
  }

  // 2. Pointwise source on rows 0..51: vacuum + matter (frozen fluid, live metric) - S_eq.
  //    Alpha as ln alpha: S_9(u)/alpha_u - S_9(e)/alpha_e keeps C_A exact.
  gkyl_range_iter_init(&it, &app->local);
  while (gkyl_range_iter_next(&it)) {
    long l = gkyl_range_idx(&app->local_ext, it.idx);
    const double *q = gkyl_array_cfetch(u, l);
    double *r = gkyl_array_fetch(rate, l), s[77], mrhs[77] = { 0 };
    gkyl_wv_eqn_source(esp->equation, q, s);
    if (fluid_f) {
      double qf[71];
      memcpy(qf, gkyl_array_cfetch(fluid_f, l), sizeof qf);
      sync_fluid_metric_from_einstein(q, qf, 0.0, true);
      add_einstein_matter_source(q, qf, gas_gamma, 0.0, true, evolution, mrhs);
    }
    for (int m = 0; m < 52; ++m) if (m != 9) r[m] += s[m] + mrhs[m];
    r[9] = (s[9] + mrhs[9])/q[9];
    if (e) {
      const double *qe = gkyl_array_cfetch(e, l), *se = gkyl_array_cfetch(seq, l);
      for (int m = 0; m < 52; ++m) if (m != 9) r[m] -= se[m];
      r[9] -= se[9]/qe[9];
    }
  }

  // 3. chi: S_pure(u) - S_pure(e) on all of local_ext (ghost states are BC-filled).
  //    The first call also fills the S_pure(e) cache, which marks F_d(e) as cached too.
  bool fill = e && !w->spe;
  if (fill) w->spe = gkyl_array_new(GKYL_DOUBLE, 1, app->local_ext.volume);
  gkyl_range_iter_init(&it, &app->local_ext);
  while (gkyl_range_iter_next(&it)) {
    long l = gkyl_range_idx(&app->local_ext, it.idx);
    double *sp = gkyl_array_fetch(w->sp, l);
    double se = 0.0;
    if (e) {
      double *c = gkyl_array_fetch(w->spe, l);
      if (fill) c[0] = route1_spure(esp, gkyl_array_cfetch(e, l));
      se = c[0];
    }
    sp[0] = route1_spure(esp, gkyl_array_cfetch(u, l)) - se;
  }
  if (cb) {
    struct route1_drv drv = { .on = ndim >= 2 && mixed_experimental_flag_enabled("GKYL_GAMMA_DRIVER")
      && !mixed_dyn_collapse_hold_phase(app) && esp->gamma_driver_reference, .mu = 0.75, .eta = 1.0,
      .alpha_floor = gamma_driver_alpha_floor(esp) };
    if (drv.on) {
      drv.mu = atof(getenv("GKYL_GAMMA_DRIVER")); if (drv.mu == 0.0) drv.mu = 0.75;
      if (getenv("GKYL_GAMMA_DRIVER_ETA")) drv.eta = atof(getenv("GKYL_GAMMA_DRIVER_ETA"));
    }
    route1_charbc_setup(app, esp, &drv);
    route1_charbc_divergence(app, esp, u, e, e ? w->spe : NULL, &drv, w->G, rate);
  }

  // 3b. GKYL_ROUTE1_BAND=W: upwind band before the non-periodic outer faces. Each face flux gets the Rusanov
  //     term -1/2 s w (v_R - v_L) on the evolved rows (gt, ln alpha, K, D, A, V, beta, chi; v = q - e), with one
  //     global speed s and the ramp w of the depth of the face's outer cell. One value per face, used by both
  //     cells: conservative. The chi part goes into S_pure, so P = D0 (chi rate) and H carry it; B follows the
  //     beta rate (step 5). GKYL_ROUTE1_BAND_DERIVED: no face terms on D, A, V, which then take the parents'
  //     derived increments (route1_band_derive).
  int band_w = route1_band_width();
  if (band_w > 0) {
    double s = route1_band_speed(app, esp, u);
    bool derived = route1_band_derived();
    if (!w->band) w->band = gkyl_array_new(GKYL_DOUBLE, nc, app->local_ext.volume);
    if (derived && !w->bandv) w->bandv = gkyl_array_new(GKYL_DOUBLE, 10, app->local_ext.volume);
    gkyl_array_clear(w->band, 0.0);
    gkyl_range_iter_init(&it, &app->local);
    while (gkyl_range_iter_next(&it)) {
      int dc = route1_depth(app, it.idx);
      if (dc > band_w) continue;
      long l = gkyl_range_idx(&app->local_ext, it.idx);
      double *bd = gkyl_array_fetch(w->band, l);
      const double *q = gkyl_array_cfetch(u, l), *qe = e ? gkyl_array_cfetch(e, l) : NULL;
      for (int d = 0; d < ndim; ++d) {
        for (int side = -1; side <= 1; side += 2) {
          int jn[GKYL_MAX_DIM];
          for (int k = 0; k < ndim; ++k) jn[k] = it.idx[k];
          jn[d] += side;
          int dn = route1_depth(app, jn);
          double wf = route1_band_weight(dc < dn ? dc : dn, band_w);
          if (wf == 0.0) continue;
          long ln = gkyl_range_idx(&app->local_ext, jn);
          const double *qn = gkyl_array_cfetch(u, ln), *qne = e ? gkyl_array_cfetch(e, ln) : NULL;
          double c = 0.5*s*wf/app->grid.dx[d];
          for (int m = 0; m < 65; ++m) {
            if ((m >= 55 && m < 64) || (derived && m >= 19 && m < 52)) continue;
            double a = m == 9 ? log(q[9]) : q[m], b = m == 9 ? log(qn[9]) : qn[m];
            if (e) {
              a -= m == 9 ? log(qe[9]) : qe[m];
              b -= m == 9 ? log(qne[9]) : qne[m];
            }
            bd[m] += c*(b - a);
          }
        }
      }
      ((double *) gkyl_array_fetch(w->sp, l))[0] += bd[64];
    }
    if (derived) route1_band_derive(app, u, w->band, w->bandv, band_w);
  }

  // 4. Gamma-driver rate for beta (production gating and reference).
  bool driver = ndim >= 2 && mixed_experimental_flag_enabled("GKYL_GAMMA_DRIVER")
    && !mixed_dyn_collapse_hold_phase(app) && esp->gamma_driver_reference;
  double mu = 0.75, eta = 1.0;
  if (driver) {
    mu = atof(getenv("GKYL_GAMMA_DRIVER")); if (mu == 0.0) mu = 0.75;
    if (getenv("GKYL_GAMMA_DRIVER_ETA")) eta = atof(getenv("GKYL_GAMMA_DRIVER_ETA"));
  }
  // GKYL_GAMMA_DRIVER_WEIGHT=kappa: mu = kappa alpha^2 chi^2 per cell (as stage1_review_patches/gauge_weight_2026_10_31);
  // the shift speeds are then sqrt(kappa) and sqrt(4 kappa/3) of the light speed, below it and below the 1+log lapse
  // speed for kappa < 3/4 and alpha < 3/(2 kappa). Default off: the constant mu above.
  double kappa = getenv("GKYL_GAMMA_DRIVER_WEIGHT") ? atof(getenv("GKYL_GAMMA_DRIVER_WEIGHT")) : 0.0;
  static bool kappa_said = false;
  if (driver && kappa > 0.0 && !kappa_said) {
    gkyl_moment_app_cout(app, stderr, "[ROUTE1] driver mu = %g alpha^2 chi^2 per cell (GKYL_GAMMA_DRIVER_WEIGHT)\n", kappa);
    kappa_said = true;
  }
  double alpha_floor = gamma_driver_alpha_floor(esp);
  bool ko_ext = route1_ko_ghosts() || route1_peel();
  bool ko_sbp = route1_ko_sbp();
  if (ko_sbp && ko_ext) {
    fprintf(stderr, "[ROUTE1] GKYL_ROUTE1_KO_SBP reads no ghost cell; it refuses KO_GHOSTS and PEEL\n"); exit(2);
  }

  // 5. KO on w = u - e (ln alpha for alpha), rows 0..54 and 64..67; G = 0 if the stencil leaves local
  //    (or local_ext with GKYL_ROUTE1_KO_GHOSTS). Along a periodic direction it reads the ghosts (no cut).
  //    Then the grid-defined rows: P (local + one ghost layer), beta, B, H.
  gkyl_range_iter_init(&it, &app->local);
  while (gkyl_range_iter_next(&it)) {
    long l = gkyl_range_idx(&app->local_ext, it.idx);
    double *r = gkyl_array_fetch(rate, l);
    const double *q = gkyl_array_cfetch(u, l);
    r[64] = ((const double *)gkyl_array_cfetch(w->sp, l))[0];
    if (driver && gamma_driver_cell_active(q, alpha_floor, true)) {
      const double *ref = gkyl_array_cfetch(esp->gamma_driver_reference, l);
      double Gam[3], mub = kappa > 0.0 ? kappa*q[9]*q[9]*q[64]*q[64] : mu, etab = eta;
      route1_driver_buffer(app, it.idx, &mub, &etab);
      if (route1_driver_metric()) route1_gamma_metric(app, u, it.idx, Gam);
      else route1_gamma(q, Gam);
      for (int i = 0; i < ndim; ++i) r[52+i] = mub*(Gam[i] - ref[i]) - etab*(q[52+i] - ref[3+i]);
    }
    if (route1_zdamp() > 0.0) { // [ZD]
      double cv[3];
      route1_cv(q, cv);
      if (e) {
        double ce[3];
        route1_cv(gkyl_array_cfetch(e, l), ce);
        for (int i = 0; i < 3; ++i) cv[i] -= ce[i];
      }
      for (int i = 0; i < 3; ++i) r[49 + i] -= route1_zdamp()*cv[i];
    }
    if (band_w > 0) { // chi's part is already in S_pure
      const double *bd = gkyl_array_cfetch(w->band, l);
      for (int m = 0; m < 55; ++m) r[m] += bd[m];
    }
    if (sigma > 0) {
      for (int d = 0; d < ndim; ++d) {
        double c = sigma/(16*app->grid.dx[d]);
        // faces i-1/2 (stencil i-2..i+1) and i+1/2 (stencil i-1..i+2)
        #define ROUTE1_KO_IN(o) (ko_ext || route1_global_periodic(app, d) ? route1_in_ext(app, it.idx, d, o) \
          : route1_in_local(app, it.idx, d, o))
        bool fm = ROUTE1_KO_IN(-2) && ROUTE1_KO_IN(1);
        bool fp = ROUTE1_KO_IN(-1) && ROUTE1_KO_IN(2);
        if (!ko_sbp && !fm && !fp) continue;
        double mw[5]; // GKYL_ROUTE1_KO_SBP weights of cells i-2..i+2
        for (int o = -2; o <= 2; ++o) mw[o+2] = route1_ko_sbp_weight(app, d, it.idx[d] + o);
        const double *s[5];
        const double *se[5] = { 0 };
        for (int o = -2; o <= 2; ++o) {
          bool ok = ROUTE1_KO_IN(o) || o == 0;
          long lo = ok ? route1_nb(app, it.idx, d, o) : l;
          s[o+2] = gkyl_array_cfetch(u, lo);
          if (e) se[o+2] = gkyl_array_cfetch(e, lo);
        }
        double cp = c, cm = c; // [HKO] face coefficients at i+1/2 and i-1/2
        double rs_p = 0.0, rs_m = 0.0; // [RUS] face speeds at i+1/2 and i-1/2 (0 = no Rusanov term)
        if (route1_hko_r > 0.0) {
          double rr[3];
          for (int o = -1; o <= 1; ++o) {
            int jdx[GKYL_MAX_DIM];
            for (int dd = 0; dd < ndim; ++dd) jdx[dd] = it.idx[dd];
            jdx[d] += o;
            double xc[3] = { 0, 0, 0 };
            gkyl_rect_grid_cell_center(&app->grid, jdx, xc);
            rr[o+1] = sqrt(xc[0]*xc[0] + xc[1]*xc[1] + xc[2]*xc[2]);
          }
          double cin = route1_rus() ? 0.0 : fmax(c, route1_hko_sig/(16*app->grid.dx[d]));
          bool mp = rr[1] < route1_hko_r || rr[2] < route1_hko_r, mm = rr[1] < route1_hko_r || rr[0] < route1_hko_r;
          if (mp) cp = cin;
          if (mm) cm = cin;
          if (route1_rus()) { // [RUS] face speeds; the KO term is off on these faces (cin = 0)
            double sc[3];
            for (int o = -1; o <= 1; ++o) {
              const double *qo = s[o+2];
              double muc = kappa > 0.0 ? kappa*qo[9]*qo[9]*qo[64]*qo[64] : mu;
              sc[o+1] = route1_rus_speed(esp, qo, muc);
            }
            rs_p = mp && ROUTE1_KO_IN(1) ? fmax(sc[1], sc[2]) : 0.0;
            rs_m = mm && ROUTE1_KO_IN(-1) ? fmax(sc[1], sc[0]) : 0.0;
          }
        }
        for (int m = 0; m < 68; ++m) {
          if (m >= 55 && m < 64 && !cb) continue; // B: KO of its own only when it is a flux variable
          double v[5];
          for (int o = 0; o < 5; ++o) {
            double a = m == 9 ? log(s[o][9]) : s[o][m];
            double b = e ? (m == 9 ? log(se[o][9]) : se[o][m]) : 0.0;
            v[o] = a - b;
          }
          double gp, gm;
          if (ko_sbp) { // G_{i+1/2} = m_{i+1} s_{i+1} - m_i s_i; s only where its weight is nonzero
            double sm = mw[1] ? v[0] - 2*v[1] + v[2] : 0.0, s0 = mw[2] ? v[1] - 2*v[2] + v[3] : 0.0,
                   sp = mw[3] ? v[2] - 2*v[3] + v[4] : 0.0;
            gp = mw[3]*sp - mw[2]*s0;
            gm = mw[2]*s0 - mw[1]*sm;
          }
          else {
            gp = fp ? v[4] - 3*v[3] + 3*v[2] - v[1] : 0.0;
            gm = fm ? v[3] - 3*v[2] + 3*v[1] - v[0] : 0.0;
          }
          r[m] -= cp == cm ? cp*(gp - gm) : cp*gp - cm*gm; // equal faces: the frozen expression, bit for bit
          if (rs_p > 0.0 || rs_m > 0.0) // [RUS] -(F_{i+1/2} - F_{i-1/2})/dx of the face terms -1/2 s (v_R - v_L)
            r[m] += 0.5*(rs_p*(v[3] - v[2]) - rs_m*(v[2] - v[1]))/app->grid.dx[d];
        }
        #undef ROUTE1_KO_IN
      }
    }
    int fc_fields;
    if (route1_fc(&fc_fields) > 0.0) { // [FC]
      for (int d = 0; d < ndim; ++d) {
        double h = app->grid.dx[d], nu = route1_fc(NULL)*h;
        bool per = route1_global_periodic(app, d);
        bool in[3];
        const double *s[3], *se[3] = { 0 };
        for (int o = -1; o <= 1; ++o) {
          in[o+1] = o == 0 || (per ? route1_in_ext(app, it.idx, d, o) : route1_in_local(app, it.idx, d, o));
          long lo = in[o+1] ? route1_nb(app, it.idx, d, o) : l;
          s[o+1] = gkyl_array_cfetch(u, lo);
          if (e) se[o+1] = gkyl_array_cfetch(e, lo);
        }
        for (int m = 0; m < 13; ++m) { // 0: ln alpha; 1..9: gt_ij; 10..12: beta^j
          if (m == 0 && !(fc_fields & 1)) continue;
          if (m >= 1 && m <= 9 && !(fc_fields & 2)) continue;
          if (m >= 10 && !(fc_fields & 4)) continue;
          int pr = m == 0 ? 9 : m <= 9 ? m - 1 : 52 + (m - 10);
          int dr = m == 0 ? 46 + d : m <= 9 ? 19 + 9*d + (m - 1) : 55 + 3*d + (m - 10);
          double f = m == 0 ? 1.0 : 2.0, pv[3], dv[3];
          for (int o = 0; o < 3; ++o) {
            pv[o] = (m == 0 ? log(s[o][9]) : s[o][pr]) - (e ? (m == 0 ? log(se[o][9]) : se[o][pr]) : 0.0);
            dv[o] = s[o][dr] - (e ? se[o][dr] : 0.0);
          }
          double gp = in[2] ? nu*((pv[2] - pv[1])/h - 0.5*f*(dv[1] + dv[2])) : 0.0;
          double gm = in[0] ? nu*((pv[1] - pv[0])/h - 0.5*f*(dv[0] + dv[1])) : 0.0;
          r[pr] += (gp - gm)/h;
        }
      }
    }
  }
  // NewRad: radiative ghost rates, then P on every cell: D_k (chi rate) (+ KO already added on local).
  bool nr = route1_newrad_enabled();
  if (nr) {
    route1_newrad_rates(app, esp, u, e, rate, w->sp);
    route1_periodic_sync(app, w->sp); // slab: periodic ghosts, and edge ghosts from their seam copies
    gkyl_range_iter_init(&it, &app->local_ext);
    while (gkyl_range_iter_next(&it)) {
      int s[3];
      long l = gkyl_range_idx(&app->local_ext, it.idx);
      double *lp = gkyl_array_fetch(w->lp, l);
      double *r = route1_ghost_dir(app, it.idx, s) == 0 ? gkyl_array_fetch(rate, l) : NULL;
      for (int k = 0; k < 3; ++k) {
        lp[k] = (k < ndim ? route1_dr(app, w->sp, 0, it.idx, k) : 0.0) + (r ? r[65+k] : 0.0);
        if (r) r[65+k] = lp[k];
      }
    }
  }
  // GKYL_ROUTE1_CHARBC: P's rate (flux divergence + KO) is complete on local; lp = it.
  gkyl_range_iter_init(&it, &app->local);
  while (cb && gkyl_range_iter_next(&it)) {
    long l = gkyl_range_idx(&app->local_ext, it.idx);
    const double *r = gkyl_array_cfetch(rate, l);
    double *lp = gkyl_array_fetch(w->lp, l);
    for (int k = 0; k < 3; ++k) lp[k] = r[65+k];
  }
  // P on local plus the face-neighbour ghosts: D0_k S_pure (+ KO already added on local).
  gkyl_range_iter_init(&it, &app->local_ext);
  while (!nr && !cb && gkyl_range_iter_next(&it)) {
    int out = 0;
    for (int d = 0; d < ndim; ++d)
      out += it.idx[d] < app->local.lower[d] || it.idx[d] > app->local.upper[d]
        ? (it.idx[d] < app->local.lower[d]-1 || it.idx[d] > app->local.upper[d]+1 ? 2 : 1) : 0;
    if (out > 1) continue;
    long l = gkyl_range_idx(&app->local_ext, it.idx);
    double *lp = gkyl_array_fetch(w->lp, l);
    double *r = out == 0 ? gkyl_array_fetch(rate, l) : NULL;
    for (int k = 0; k < 3; ++k) {
      double v = 0.0;
      if (k < ndim) {
        const double *p = gkyl_array_cfetch(w->sp, route1_nb(app, it.idx, k, 1));
        const double *m = gkyl_array_cfetch(w->sp, route1_nb(app, it.idx, k, -1));
        v = (p[0] - m[0])/(2*app->grid.dx[k]);
      }
      lp[k] = v + (r ? r[65+k] : 0.0);
      if (r) r[65+k] = lp[k];
    }
  }
  // Periodic directions: the stencils below read lp and the beta rate across the seam.
  route1_periodic_sync(app, w->lp);
  route1_periodic_sync(app, rate);
  // H = sym D0 L(P); B = 1/2 D0 L(beta), one-sided inward at a non-periodic global face (as production)
  // unless GKYL_ROUTE1_B_CENTRED with NewRad.
  bool b_centred = nr && (route1_b_centred() || route1_peel());
  // GKYL_ROUTE1_CHARBC: H = sym d L(P) with route1_charbc_d (one-sided on the boundary cells); B already has its rate.
  gkyl_range_iter_init(&it, &app->local);
  while (cb && gkyl_range_iter_next(&it)) {
    double *r = gkyl_array_fetch(rate, gkyl_range_idx(&app->local_ext, it.idx));
    double dlp[3][3] = { { 0 } };
    for (int j = 0; j < ndim; ++j)
      for (int k = 0; k < 3; ++k) dlp[j][k] = route1_charbc_d(app, w->lp, k, it.idx, j);
    for (int k = 0; k < 3; ++k)
      for (int j = 0; j < 3; ++j) r[68+3*k+j] = 0.5*(dlp[k][j] + dlp[j][k]);
  }
  gkyl_range_iter_init(&it, &app->local);
  while (!cb && gkyl_range_iter_next(&it)) {
    long l = gkyl_range_idx(&app->local_ext, it.idx);
    double *r = gkyl_array_fetch(rate, l);
    double dlp[3][3] = { { 0 } };
    for (int j = 0; j < ndim; ++j) {
      const double *p = gkyl_array_cfetch(w->lp, route1_nb(app, it.idx, j, 1));
      const double *m = gkyl_array_cfetch(w->lp, route1_nb(app, it.idx, j, -1));
      for (int k = 0; k < 3; ++k) dlp[j][k] = (p[k] - m[k])/(2*app->grid.dx[j]);
    }
    for (int k = 0; k < 3; ++k)
      for (int j = 0; j < 3; ++j) r[68+3*k+j] = 0.5*(dlp[k][j] + dlp[j][k]);
    for (int i = 0; i < ndim; ++i) {
      double h = app->grid.dx[i];
      bool per = route1_global_periodic(app, i);
      int inward = 0;
      if (!per && !b_centred && it.idx[i] == app->global.lower[i]) inward = 1;
      if (!per && !b_centred && it.idx[i] == app->global.upper[i]) inward = -1;
      for (int j = 0; j < 3; ++j) {
        double der;
        if (inward) {
          const double *r1 = gkyl_array_cfetch(rate, route1_nb(app, it.idx, i, inward));
          const double *r2 = gkyl_array_cfetch(rate, route1_nb(app, it.idx, i, 2*inward));
          der = inward*(-3*r[52+j] + 4*r1[52+j] - r2[52+j])/(2*h);
        }
        else {
          const double *rp = gkyl_array_cfetch(rate, route1_nb(app, it.idx, i, 1));
          const double *rm = gkyl_array_cfetch(rate, route1_nb(app, it.idx, i, -1));
          der = (rp[52+j] - rm[52+j])/(2*h);
        }
        r[55+3*i+j] = 0.5*der;
      }
    }
  }
  if (nr) route1_newrad_derived(app, w->lp, rate);
}

// out = a*x + b*(y + dt*rate) on rng, slot 9 geometric (ln alpha). a = 0 means out = y + dt*rate.
static void
route1_combine(gkyl_moment_app *app, const struct gkyl_range *rng, int nc, double a, const struct gkyl_array *x,
  double b, const struct gkyl_array *y, double dt, const struct gkyl_array *rate, struct gkyl_array *out)
{
  struct gkyl_range_iter it;
  gkyl_range_iter_init(&it, rng);
  while (gkyl_range_iter_next(&it)) {
    long l = gkyl_range_idx(&app->local_ext, it.idx);
    const double *xs = a != 0 ? gkyl_array_cfetch(x, l) : NULL, *ys = gkyl_array_cfetch(y, l);
    const double *r = gkyl_array_cfetch(rate, l);
    double *o = gkyl_array_fetch(out, l);
    for (int m = 0; m < nc; ++m) {
      if (m == 9) {
        double ly = log(ys[9]) + dt*r[9];
        o[9] = a != 0 ? exp(a*log(xs[9]) + b*ly) : exp(ly);
      }
      else {
        double yv = ys[m] + dt*r[m];
        o[m] = a != 0 ? a*xs[m] + b*yv : yv;
      }
    }
  }
}

// One SSP-RK3 step of f[0] in place. Returns false on an invalid stage (alpha or chi <= 0, NaN).
static bool
route1_step(gkyl_moment_app *app, struct moment_species *esp, const struct gkyl_array *fluid_f,
  double gas_gamma, double tcurr, double dt, const struct gkyl_array *e, const struct gkyl_array *seq)
{
  struct route1_work *w = &route1_w;
  int nc = esp->num_equations;
  if (esp->vacuum_einstein_conformal_excision_threshold != 0) {
    gkyl_moment_app_cout(app, stderr, "[ROUTE1] excision is not supported; refusing the step\n");
    return false;
  }
  if (!w->L) {
    long vol = app->local_ext.volume;
    w->L = gkyl_array_new(GKYL_DOUBLE, nc, vol);
    w->G = gkyl_array_new(GKYL_DOUBLE, nc, vol);
    w->sp = gkyl_array_new(GKYL_DOUBLE, 1, vol);
    w->lp = gkyl_array_new(GKYL_DOUBLE, 3, vol);
    w->u1 = gkyl_array_new(GKYL_DOUBLE, nc, vol);
    w->u2 = gkyl_array_new(GKYL_DOUBLE, nc, vol);
    gkyl_array_clear(w->sp, 0.0); gkyl_array_clear(w->lp, 0.0);
    route1_report(app, esp->f[0], tcurr);
  }
  // Driver reference: captured at activation exactly as production (dt = 0 call); with
  // GKYL_ROUTE1_DRIVER_METRIC its Gammatilde part is captured again below, once the ghosts are valid.
  bool ref_new = false;
  if (mixed_experimental_flag_enabled("GKYL_GAMMA_DRIVER") && !mixed_dyn_collapse_hold_phase(app)
      && !esp->gamma_driver_reference) {
    const char *ld = getenv("GKYL_FX_DRIVER_REF_LOAD"), *sv = getenv("GKYL_FX_DRIVER_REF_SAVE");
    if (ld) { // [FX] restart: the reference captured at activation in the original run
      esp->gamma_driver_reference = gkyl_array_new(GKYL_DOUBLE, 6, app->local_ext.volume);
      FILE *fp = fopen(ld, "rb");
      size_t want = esp->gamma_driver_reference->size*6, got = fp ? fread(esp->gamma_driver_reference->data, sizeof(double), want, fp) : 0;
      if (fp) fclose(fp);
      if (got != want) { gkyl_moment_app_cout(app, stderr, "[FX] driver reference %s: read %zu of %zu doubles\n", ld, got, want); exit(3); }
      gkyl_moment_app_cout(app, stderr, "[FX] t=%.17g driver reference loaded from %s\n", tcurr, ld);
    }
    else gamma_driver_shift_rk2(app, esp, tcurr, 0.0);
    ref_new = esp->gamma_driver_reference != NULL;
    if (sv && esp->gamma_driver_reference && !ld) {
      FILE *fp = fopen(sv, "wb");
      size_t n = fp ? fwrite(esp->gamma_driver_reference->data, sizeof(double), esp->gamma_driver_reference->size*6, fp) : 0;
      if (fp) fclose(fp);
      gkyl_moment_app_cout(app, stderr, "[FX] t=%.17g driver reference saved to %s (%zu doubles)\n", tcurr, sv, n);
    }
  }

  struct gkyl_array *u0 = esp->f[0], *u1 = w->u1, *u2 = w->u2, *L = w->L;
  if (seq == NULL) e = NULL;
  if (route1_symbol()) { // diagnostic: [SYMBOL] lines (=1) or the [SYM2] variant scan (=2), then stop
    if (atoi(getenv("GKYL_ROUTE1_SYMBOL")) == 2) route1_symbol_scan(app, esp, u0); else route1_symbol_lines(app, esp, u0);
    exit(0);
  }
  // NewRad: the ghosts are evolved with the stages (no BC fill inside the step).
  bool nr = route1_newrad_enabled();
  if (nr && app->num_periodic_dir >= app->ndim) {
    gkyl_moment_app_cout(app, stderr, "[ROUTE1] newrad needs at least one non-periodic direction; refusing the step\n");
    return false;
  }
  bool peel = route1_peel();
  if (route1_charbc()) {
    bool same = true;
    for (int d = 0; d < app->ndim; ++d)
      same = same && app->local.lower[d] == app->global.lower[d] && app->local.upper[d] == app->global.upper[d];
    if (nr || peel || route1_band_width() > 0 || route1_ko_ghosts() || route1_ko_sbp()
        || getenv("GKYL_ROUTE1_DRIVER_BUFFER") || getenv("GKYL_GAMMA_DRIVER_WEIGHT") || !same) {
      gkyl_moment_app_cout(app, stderr, "[ROUTE1] charbc refuses NEWRAD, PEEL, BAND, KO_GHOSTS, KO_SBP, DRIVER_BUFFER, GAMMA_DRIVER_WEIGHT and"
        " domain decomposition; refusing the step\n");
      return false;
    }
  }
  if (peel && (route1_band_width() > 0 || getenv("GKYL_ROUTE1_DRIVER_BUFFER"))) {
    gkyl_moment_app_cout(app, stderr, "[ROUTE1] peel runs with no band and no buffer; refusing the step\n");
    return false;
  }
  const struct gkyl_range *rng = nr ? &app->local_ext : &app->local;
  if (!nr || !w->gh) moment_species_apply_bc(app, tcurr, esp, u0);
  if (nr) { route1_newrad_ghosts(app, u0, false); route1_periodic_sync(app, u0); }
  if (ref_new && route1_driver_metric()) route1_driver_metric_reference(app, esp, u0);
  if (route1_idb_enabled() && !w->idb0) route1_idb(app, u0, tcurr);
  route1_rhs(app, esp, u0, e, seq, fluid_f, gas_gamma, tcurr, L);
  gkyl_array_copy_range(u1, u0, &app->local_ext);
  route1_combine(app, rng, nc, 0.0, NULL, 1.0, u0, dt, L, u1);
  if (nr) route1_periodic_sync(app, u1);
  if (peel) route1_peel_derive(app, u1);
  if (!nr) moment_species_apply_bc(app, tcurr+dt, esp, u1);
  route1_rhs(app, esp, u1, e, seq, fluid_f, gas_gamma, tcurr+dt, L);
  gkyl_array_copy_range(u2, u0, &app->local_ext);
  route1_combine(app, rng, nc, 0.75, u0, 0.25, u1, dt, L, u2);
  if (nr) route1_periodic_sync(app, u2);
  if (peel) route1_peel_derive(app, u2);
  if (!nr) moment_species_apply_bc(app, tcurr+0.5*dt, esp, u2);
  route1_rhs(app, esp, u2, e, seq, fluid_f, gas_gamma, tcurr+0.5*dt, L);
  route1_combine(app, rng, nc, 1.0/3.0, u0, 2.0/3.0, u2, dt, L, u1);
  if (nr) route1_periodic_sync(app, u1);
  if (peel) route1_peel_derive(app, u1);
  bool ok = true;
  struct gkyl_range_iter it;
  gkyl_range_iter_init(&it, rng);
  while (gkyl_range_iter_next(&it)) {
    const double *q = gkyl_array_cfetch(u1, gkyl_range_idx(&app->local_ext, it.idx));
    if (!(q[9] > 0 && q[64] > 0 && isfinite(q[10]) && isfinite(q[19]))) { ok = false; break; }
  }
  if (!ok) return false;
  gkyl_array_copy_range(u0, u1, rng);
  if (nr) route1_newrad_ghosts(app, u0, true);
  else moment_species_apply_bc(app, tcurr+dt, esp, u0);
  w->steps++;
  const char *ck = getenv("GKYL_ROUTE1_CHECK");
  long every = ck ? atol(ck) : 20;
  if (every > 0 && w->steps % every == 0) route1_report(app, u0, tcurr+dt);
  route1_monitor(app, u0, e, tcurr+dt);
  return true;
}

// Ledger of the fluid's stored conserved rows over every accepted coupled step (this build; reads only, this
// rank). Rows D, S_x, S_y, S_z, tau (fluid rows 0-4), each summed times the cell volume over the local cells, are
// taken at the step's start, after each fluid MP SSP-RK3 update, after each excision freeze and after each
// atmosphere repair. Buckets: flux (the updates: flux difference, plus the geometric source for S_i and tau),
// excis (mixed_restore_excised_fluid_cells: excised cells set back to their pre-update value), reset
// (mixed_gr_euler_conserved_atmosphere_reset) and between (from one accepted step's end to the next one's start:
// the collapse trigger, driver-level edits). Only flux and source may change the stored state; excis, reset and
// between are direct writes. The repair is also counted cell by cell (route1_ledger_cell, hooked where it writes
// rows 0-4), by mechanism: proj (in-band projection toward the atmosphere, every row), lift (out-of-band energy
// lift, tau) and cap (out-of-band momentum cap, S_i), with cell counts, signed sums and the absolute sums |dS|
// (over i) and |dtau|; check is the largest per-step difference between those per-cell sums and the reset bucket
// (round-off); sites counts the repaired cells by call: after fluid update 1, after the spacetime step and the
// metric copy, after fluid update 2 (which part of the step leaves cells inadmissible). With no D source a flux
// update changes Sum D only through the outer faces: face integrates dt times sum over the outer faces of
// max |D| (alpha + |beta|) dA of the boundary cell and its first ghost, a bound on that flux for |v| < 1, and
// excess is the largest flux-update change of one step beyond its bound (> 0: the update changes Sum D inside the
// box). Every 0.5 time units, cumulative since the first step: [LEDGER] for D (resid = dSumD - flux - reset -
// between - excis, round-off) and [LEDGER-NF] for every direct write, per row.
// floor n counts, per call site, the cells the GKYL_ROUTE1_FLOOR_ONLY validator accepted outside the cone (accepted
// steps); short counts boosted fluid updates the CFL shortened (each redone at the shorter length). Under
// GKYL_ROUTE1_FLOOR_ONLY the ledger asserts reset = excis = 0 exactly in every accepted step and stops the run
// otherwise.
enum { R1LG_FLUX, R1LG_RESET, R1LG_EXCISION };   // buckets of route1_ledger_mark
enum { R1LG_PROJ, R1LG_LIFT, R1LG_CAP };         // repair mechanisms of route1_ledger_cell

static struct {
  bool init;
  long last_slot;
  double sum0[5], end[5], cur[5], between[5], face, excess, check;
  double tb[3][5], b[3][5];   // change per bucket and row: this step, accepted steps
  double tm[3][7], m[3][7];   // repair per mechanism: signed rows 0-4, |dS|, |dtau|
  long tn[3], n[3];           // repaired cells per mechanism
  long ts[3], s[3];           // repaired cells per call site
  long ta[3], a[3];           // cells accepted outside the cone per call site (GKYL_ROUTE1_FLOOR_ONLY)
  long nshort;                // boosted fluid updates shortened by the CFL
  int nflux, nreset;          // marks so far in this step
} route1_lg;

// Fluid source and atmosphere (this build). The source's lapse-gradient part, dtau/dt = -S^i d_i alpha and
// dS_j/dt = -E d_j alpha (E = D + tau; rows densitized; sources_explicit.c), is a boost generator on (E, S_j): its
// exact flow keeps D and E^2 - gamma^ij S_i S_j, hence the admissible cone E^2 - D^2 - S^2 > 0, for any step. A
// forward-Euler stage lowers E^2 - S^2 by g^2 h^2 (E^2 - S_par^2) (g = |d alpha|), which takes cold matter
// (eps < (g h)^2 / 2) out of the cone (AUDIT_FLUID_STATE_WRITES_2026_11_04.md).
// GKYL_ROUTE1_BOOST=1 (while the fluid's well-balancing is off): each fluid MP SSP-RK3 update of length h becomes
// exact boost h/2, RK3 with the rest of the source, exact boost h/2 (Strang, second order). The RK3 sees zero in the
// fluid's d alpha rows 28-30, which enter only those two terms; they are restored after it. Frozen geometry, as in
// the RK3. Default 0: the whole source inside the RK3 (production).
// GKYL_ROUTE1_BOOST=2: the integrating-factor (Lawson) form of the same RK3 instead, with the exact boost as the
// factor (route1_fluid_lawson): the boosts sit between the stages, so every flux sees a stage state rather than one
// kicked by h/2 of gravity, and the update is third order like the RK3.
// GKYL_ROUTE1_FLOOR_ONLY=1: the allowed atmosphere only. The conserved repair
// (mixed_gr_euler_conserved_atmosphere_reset) is skipped, and the validator accepts a cell outside the cone, or with
// D or E <= 0, whose density marker is at most GKYL_ATM_CONS_RESET_FACTOR (default 100) rho_atm, the band the repair
// projected: the C2P floor (GKYL_ATM_RESET) evaluates its fluxes and sources at the static atmosphere, and its
// stored state changes only by flux and source. Denser cells must stay in the cone (the step is rejected otherwise).
static bool
route1_boost_enabled(void)
{
  const char *v = getenv("GKYL_ROUTE1_BOOST");
  return v && atoi(v) != 0;
}

static bool
route1_boost_lawson(void)
{
  const char *v = getenv("GKYL_ROUTE1_BOOST");
  return v && atoi(v) == 2;
}

static bool
route1_floor_only(void)
{
  const char *v = getenv("GKYL_ROUTE1_FLOOR_ONLY");
  return v && atoi(v) != 0;
}

static void
route1_ledger_sums(const gkyl_moment_app *app, const struct gkyl_array *f, double s[5])
{
  for (int c = 0; c < 5; ++c) s[c] = 0.0;
  struct gkyl_range_iter it;
  gkyl_range_iter_init(&it, &app->local);
  while (gkyl_range_iter_next(&it)) {
    const double *q = gkyl_array_cfetch(f, gkyl_range_idx(&app->local_ext, it.idx));
    for (int c = 0; c < 5; ++c) s[c] += q[c];
  }
  for (int c = 0; c < 5; ++c) s[c] *= app->grid.cellVolume;
}

static double
route1_ledger_face(const gkyl_moment_app *app, const struct gkyl_array *f)
{
  double out = 0.0;
  struct gkyl_range_iter it;
  gkyl_range_iter_init(&it, &app->local);
  while (gkyl_range_iter_next(&it)) {
    for (int d = 0; d < app->ndim; ++d) {
      if (route1_global_periodic(app, d)) continue;
      int side = it.idx[d] == app->global.lower[d] ? -1 : (it.idx[d] == app->global.upper[d] ? 1 : 0);
      if (!side) continue;
      double m = 0.0;
      for (int o = 0; o <= 1; ++o) {
        const double *q = gkyl_array_cfetch(f, route1_nb(app, it.idx, d, o*side));
        m = fmax(m, fabs(q[0])*(fabs(q[5]) + sqrt(q[6]*q[6] + q[7]*q[7] + q[8]*q[8])));
      }
      out += m*app->grid.cellVolume/app->grid.dx[d];
    }
  }
  return out;
}

static void
route1_ledger_begin(gkyl_moment_app *app, int fidx)
{
  route1_fx_pend = 0.0; // [FX] reservoir pending for this step attempt
  if (fidx < 0) return;
  double s[5];
  route1_ledger_sums(app, app->species[fidx].fcurr, s);
  if (!route1_lg.init) {
    route1_lg.init = true;
    memcpy(route1_lg.sum0, s, sizeof s);
    memcpy(route1_lg.end, s, sizeof s);
    route1_lg.last_slot = -1;
  }
  for (int c = 0; c < 5; ++c) {
    route1_lg.between[c] += s[c] - route1_lg.end[c];
    route1_lg.end[c] = route1_lg.cur[c] = s[c];
  }
  memset(route1_lg.tb, 0, sizeof route1_lg.tb);
  memset(route1_lg.tm, 0, sizeof route1_lg.tm);
  memset(route1_lg.tn, 0, sizeof route1_lg.tn);
  memset(route1_lg.ts, 0, sizeof route1_lg.ts);
  memset(route1_lg.ta, 0, sizeof route1_lg.ta);
  route1_lg.nflux = route1_lg.nreset = 0;
}

static void
route1_ledger_mark(gkyl_moment_app *app, int fidx, int bucket)
{
  if (fidx < 0 || !route1_lg.init) return;
  double s[5];
  route1_ledger_sums(app, app->species[fidx].fcurr, s);
  for (int c = 0; c < 5; ++c) {
    route1_lg.tb[bucket][c] += s[c] - route1_lg.cur[c];
    route1_lg.cur[c] = s[c];
  }
  route1_lg.nflux += bucket == R1LG_FLUX;
  route1_lg.nreset += bucket == R1LG_RESET;
}

// One cell written by the atmosphere repair: its rows 0-4 before and after.
static void
route1_ledger_cell(const gkyl_moment_app *app, bool out_of_band, const double qold[5], const double qnew[5])
{
  if (!route1_lg.init) return;
  int k = !out_of_band ? R1LG_PROJ : (qnew[4] != qold[4] ? R1LG_LIFT : R1LG_CAP);
  double v = app->grid.cellVolume, *a = route1_lg.tm[k];
  for (int c = 0; c < 5; ++c) a[c] += (qnew[c] - qold[c])*v;
  for (int c = 1; c < 4; ++c) a[5] += fabs(qnew[c] - qold[c])*v;
  a[6] += fabs(qnew[4] - qold[4])*v;
  route1_lg.tn[k] += 1;
  route1_lg.ts[route1_lg.nflux >= 2 ? 2 : (route1_lg.nreset >= 1 ? 1 : 0)] += 1;
}

#ifdef ROUTE1_BARO
// [BARO] barotropic fallback in the primitive recovery (private library scan_tmp/route1_baro_2026_11_06, env
// GKYL_C2P_BARO_K; AUDIT_FLUID_STATE_WRITES, "BH formation and hold run", amendment). Built only with R1_BARO=1.
extern bool gkyl_gr_euler_c2p_baro_used(double gas_gamma, const double q[71]);
extern bool gkyl_gr_euler_c2p_baro_handled(double gas_gamma, const double q[71]);

// One [BARO] line per accepted step while any cell's recovery takes the fallback: count, cells that entered / left
// it since the previous accepted step (flip-flop check), radius range, D / sqrt(gamma) range. Reads only.
static void
route1_baro_log(gkyl_moment_app *app, int fidx, double t)
{
  static char *was = NULL;
  if (fidx < 0 || !getenv("GKYL_C2P_BARO_K")) return;
  const struct moment_species *fsp = &app->species[fidx];
  const struct gkyl_array *f = fsp->fcurr;
  if (!was) was = calloc(app->local.volume, 1);
  int n = 0, in = 0, out = 0;
  double rmin = INFINITY, rmax = 0.0, dmin = INFINITY, dmax = 0.0;
  struct gkyl_range_iter it;
  gkyl_range_iter_init(&it, &app->local);
  for (long k = 0; gkyl_range_iter_next(&it); ++k) {
    const double *q = gkyl_array_cfetch(f, gkyl_range_idx(&app->local_ext, it.idx));
    bool used = q[27] >= pow(10.0, -8.0) && gkyl_gr_euler_c2p_baro_used(fsp->gr_euler_gas_gamma, q);
    if (used != (bool) was[k]) { if (used) ++in; else ++out; was[k] = used; }
    if (!used) continue;
    ++n;
    double xc[3] = { 0, 0, 0 }, r = 0.0;
    gkyl_rect_grid_cell_center(&app->grid, it.idx, xc);
    for (int d = 0; d < app->ndim; ++d) r += xc[d]*xc[d];
    r = sqrt(r);
    double g[3][3] = { { q[9], q[10], q[11] }, { q[12], q[13], q[14] }, { q[15], q[16], q[17] } }, gi[3][3];
    double dn = q[0]/sqrt(sym3_inverse(g, gi));
    rmin = fmin(rmin, r); rmax = fmax(rmax, r); dmin = fmin(dmin, dn); dmax = fmax(dmax, dn);
  }
  if (n || out)
    gkyl_moment_app_cout(app, stderr, "[BARO] t=%.8g cells=%d entered=%d left=%d r=[%.3f,%.3f] D/sqrtg=[%.3e,%.3e]\n",
      t, n, in, out, n ? rmin : 0.0, rmax, n ? dmin : 0.0, dmax);
}
#endif

static void
route1_ledger_end(gkyl_moment_app *app, int fidx, double tcurr, double dt)
{
#ifdef ROUTE1_BARO
  route1_baro_log(app, fidx, tcurr + dt);
#endif
  route1_fx_bh += route1_fx_pend; // [FX] accepted step: commit the D that crossed into the mask
  route1_fx_pend = 0.0;
  if (fidx >= 0) route1_fx_report(app, app->species[fidx].fcurr, tcurr + dt);
  if (fidx < 0 || !route1_lg.init) return;
  const struct gkyl_array *f = app->species[fidx].fcurr;
  double s[5], fb = dt*route1_ledger_face(app, f);
  route1_ledger_sums(app, f, s);
  for (int k = 0; k < 3; ++k) {
    for (int c = 0; c < 5; ++c) route1_lg.b[k][c] += route1_lg.tb[k][c];
    for (int c = 0; c < 7; ++c) route1_lg.m[k][c] += route1_lg.tm[k][c];
    route1_lg.n[k] += route1_lg.tn[k];
    route1_lg.s[k] += route1_lg.ts[k];
    route1_lg.a[k] += route1_lg.ta[k];
  }
  for (int c = 0; c < 5 && route1_floor_only(); ++c)
    if (route1_lg.tb[R1LG_RESET][c] != 0.0 || route1_lg.tb[R1LG_EXCISION][c] != 0.0) {
      gkyl_moment_app_cout(app, stderr, "[LEDGER-ASSERT] t=%.6g row %d: reset %.6e excis %.6e; only flux and source"
        " may change the stored state. Stopping.\n", tcurr + dt, c, route1_lg.tb[R1LG_RESET][c],
        route1_lg.tb[R1LG_EXCISION][c]);
      exit(3);
    }
  for (int c = 0; c < 5; ++c)
    route1_lg.check = fmax(route1_lg.check, fabs(route1_lg.tm[R1LG_PROJ][c] + route1_lg.tm[R1LG_LIFT][c]
      + route1_lg.tm[R1LG_CAP][c] - route1_lg.tb[R1LG_RESET][c]));
  route1_lg.face += fb;
  route1_lg.excess = fmax(route1_lg.excess, fabs(route1_lg.tb[R1LG_FLUX][0]) - fb);
  memcpy(route1_lg.end, s, sizeof s);
  route1_fx_remove(app, &app->species[fidx], tcurr + dt); // [FXR] after this step's sums: booked in "between"
  long slot = (long) floor((tcurr + dt)/0.5 + 1e-6);
  if (slot == route1_lg.last_slot) return;
  route1_lg.last_slot = slot;
  double (*b)[5] = route1_lg.b, (*m)[7] = route1_lg.m, *w = route1_lg.between, *x = b[R1LG_EXCISION];
  double d = s[0] - route1_lg.sum0[0];
  gkyl_moment_app_cout(app, stderr, "[LEDGER] t=%.6g SumD=%.14e dSumD=%.6e flux=%.6e reset=%.6e between=%.6e"
    " face<=%.3e excess=%.3e resid=%.3e excis=%.6e\n", tcurr + dt, s[0], d, b[R1LG_FLUX][0], b[R1LG_RESET][0], w[0],
    route1_lg.face, route1_lg.excess, d - b[R1LG_FLUX][0] - b[R1LG_RESET][0] - w[0] - x[0], x[0]);
  gkyl_moment_app_cout(app, stderr, "[LEDGER-NF] t=%.6g proj n=%ld D=%.3e S=%.3e,%.3e,%.3e tau=%.3e |dS|=%.3e"
    " |dtau|=%.3e lift n=%ld tau=%.3e cap n=%ld |dS|=%.3e sites n=%ld,%ld,%ld excis D=%.3e S=%.3e,%.3e,%.3e"
    " tau=%.3e between D=%.3e S=%.3e,%.3e,%.3e tau=%.3e check=%.1e floor n=%ld,%ld,%ld short=%ld\n", tcurr + dt,
    route1_lg.n[R1LG_PROJ],
    m[R1LG_PROJ][0], m[R1LG_PROJ][1], m[R1LG_PROJ][2], m[R1LG_PROJ][3], m[R1LG_PROJ][4], m[R1LG_PROJ][5], m[R1LG_PROJ][6],
    route1_lg.n[R1LG_LIFT], m[R1LG_LIFT][4], route1_lg.n[R1LG_CAP], m[R1LG_CAP][5], route1_lg.s[0],
    route1_lg.s[1], route1_lg.s[2], x[0], x[1], x[2], x[3], x[4],
    w[0], w[1], w[2], w[3], w[4], route1_lg.check, route1_lg.a[0], route1_lg.a[1], route1_lg.a[2], route1_lg.nshort);
}

// Exact boost over time s on the local cells of f (excised cells skipped), d alpha from da (3 rows) when given,
// else from f's rows 28-30: with g = |d alpha|, n_j = d_j alpha / g and S_par = n^j S_j,
// E' = E cosh(gs) - S_par sinh(gs), S_par' = S_par cosh(gs) - E sinh(gs); S_j across n unchanged.
// Written as increments with cosh - 1 = 2 sinh^2(gs/2), so tau = E - D keeps its digits in cold matter.
static void
route1_boost_array(gkyl_moment_app *app, struct gkyl_array *f, const struct gkyl_array *da, double s)
{
  struct gkyl_range_iter it;
  gkyl_range_iter_init(&it, &app->local);
  while (gkyl_range_iter_next(&it)) {
    long l = gkyl_range_idx(&app->local_ext, it.idx);
    double *q = gkyl_array_fetch(f, l);
    if (q[27] < 1e-8) continue;
    double g[3][3] = { { q[9], q[10], q[11] }, { q[12], q[13], q[14] }, { q[15], q[16], q[17] } }, gi[3][3];
    if (!(sym3_inverse(g, gi) > 0.0)) continue;
    const double *a = da ? gkyl_array_cfetch(da, l) : q + 28;
    double au[3] = { 0.0, 0.0, 0.0 }, g2 = 0.0;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) au[i] += gi[i][j]*a[j];
      g2 += au[i]*a[i];
    }
    if (!(g2 > 0.0)) continue;
    double gn = sqrt(g2), sp = (au[0]*q[1] + au[1]*q[2] + au[2]*q[3])/gn, E = q[0] + q[4];
    double sh = sinh(gn*s), hh = sinh(0.5*gn*s), cm = 2.0*hh*hh;
    double dsp = sp*cm - E*sh;
    q[4] += E*cm - sp*sh;
    for (int i = 0; i < 3; ++i) q[1 + i] += dsp*a[i]/gn;
  }
}

static void
route1_fluid_boost(gkyl_moment_app *app, struct moment_species *fsp, double s)
{
  route1_boost_array(app, fsp->fcurr, NULL, s);
}

// One forward-Euler stage of the fluid alone, as production's forward_euler (mom_update_ssp_rk.c) takes it for the
// GR fluid of the mixed step, where its coupling and excision hooks do nothing: RHS, CFL, accumulate, BCs.
static void
route1_fluid_fe(gkyl_moment_app *app, struct moment_species *fsp, double t, double dt,
  const struct gkyl_array *fin, struct gkyl_array *fout, struct gkyl_update_status *st)
{
  app->stat.nfeuler += 1;
  double dtmin = moment_species_rhs(app, fsp, route1_fx_fill(app, fin), fout); // [FX] outflow fill (off: fin)
  double rel = (dt - dtmin)/dt;
  if (rel > 0 && rel < 0.01) dtmin = dt;
  double dta = st->dt_actual = dt < dtmin ? dt : dtmin;
  st->dt_suggested = dtmin;
  route1_fx_drop(app, fout, dta); // [FX] excised cells: RHS dropped, D flux into the mask recorded
  gkyl_array_accumulate_range(gkyl_array_scale_range(fout, dta, &app->local), 1.0, fin, &app->local);
  moment_species_apply_bc(app, t + dta, fsp, fout);
}

// GKYL_ROUTE1_BOOST=2: production's SSP-RK3 in integrating-factor (Lawson) form, the factor being the exact boost
// e^{sL}, and N the RHS without the lapse-gradient terms (the d alpha rows zero):
//   u1 = e^{hL} (u + h N(u)),  u2 = 3/4 e^{hL/2} u + 1/4 e^{-hL/2} (u1 + h N(u1)),
//   u' = 1/3 e^{hL} u + 2/3 e^{hL/2} (u2 + h N(u2)).
// Each stage is production's flux stage followed by boosts (of either sign, all of which keep the cone), combined
// with production's weights; the fluxes are evaluated at u, u1, u2. Production's CFL restarts. The boosts touch
// the local cells, so the BCs are applied again after each combination. Only for a mixed step whose fluid is the
// only species the RK3 updates (no field), as production's forward_euler then reduces to route1_fluid_fe.
static struct gkyl_update_status
route1_fluid_lawson(gkyl_moment_app *app, int fidx, double t, double h)
{
  struct moment_species *fsp = &app->species[fidx];
  static struct gkyl_array *fb, *da;
  long vol = app->local_ext.volume;
  if (!fb) {
    fb = gkyl_array_new(GKYL_DOUBLE, fsp->num_equations, vol);
    da = gkyl_array_new(GKYL_DOUBLE, 3, vol);
  }
  struct gkyl_array *f0 = fsp->f0, *f1 = fsp->f1, *fn = fsp->fnew;
  const struct gkyl_range *r = &app->local_ext;
  for (long l = 0; l < vol; ++l) {
    double *q = gkyl_array_fetch(f0, l), *z = gkyl_array_fetch(da, l);
    for (int c = 0; c < 3; ++c) { z[c] = q[28 + c]; q[28 + c] = 0.0; }
  }
  struct gkyl_update_status st = { .success = true };
  for (;;) {
    route1_fluid_fe(app, fsp, t, h, f0, f1, &st);
    h = st.dt_actual;
    route1_boost_array(app, f1, da, h);
    moment_species_apply_bc(app, t + h, fsp, f1);
    route1_fluid_fe(app, fsp, t + h, h, f1, fn, &st);
    if (st.dt_actual < h) { h = st.dt_actual; route1_lg.nshort += 1; continue; }
    route1_boost_array(app, fn, da, -0.5*h);
    gkyl_array_copy_range(fb, f0, r);
    route1_boost_array(app, fb, da, 0.5*h);
    gkyl_array_accumulate_range(gkyl_array_set_range(f1, 0.75, fb, r), 0.25, fn, r);
    moment_species_apply_bc(app, t + 0.5*h, fsp, f1);
    route1_fluid_fe(app, fsp, t + 0.5*h, h, f1, fn, &st);
    if (st.dt_actual < h) { h = st.dt_actual; route1_lg.nshort += 1; continue; }
    route1_boost_array(app, fn, da, 0.5*h);
    route1_boost_array(app, fb, da, 0.5*h);
    gkyl_array_accumulate_range(gkyl_array_set_range(f1, 1.0/3.0, fb, r), 2.0/3.0, fn, r);
    gkyl_array_copy_range(f0, f1, r);
    break;
  }
  for (long l = 0; l < vol; ++l) {
    double *q = gkyl_array_fetch(f0, l);
    const double *z = gkyl_array_cfetch(da, l);
    for (int c = 0; c < 3; ++c) q[28 + c] = z[c];
  }
  moment_species_apply_bc(app, t + h, fsp, f0);
  st.dt_actual = h;
  return st;
}

// GKYL_ROUTE1_FOFC=1: first-order flux correction (Lemaster & Stone 2009; Athena++'s FOFC) inside production's fluid
// SSP-RK3. After each forward-Euler stage the local cells are checked with the validator's per-cell test
// (mixed_gr_euler_state_valid). The failing cells are marked and the stage is redone from its input, with every face
// that touches a marked cell reconstructed piecewise constant. That face is shared by both cells, so the flux still
// telescopes. The mask goes through gkyl_mp_scheme_set_lo_mask, which exists only in the private library copy
// scan_tmp/route1_fofc_2026_11_05. Marks accumulate over at most 3 redos of a stage and are cleared before the next
// stage. Each stage stays a forward-Euler step of production's RHS, combined with production's weights and CFL
// restarts. When no cell fails, the update is production's. Also while well-balanced: the face then reconstructs the
// deviation from equilibrium piecewise constant (mp_scheme.c, meqn 76), so a static equilibrium stays exact. Same
// scope as route1_fluid_lawson: the fluid is the only species the RK3 updates, and there is no field.
extern void gkyl_mp_scheme_set_lo_mask(gkyl_mp_scheme *mp, const struct gkyl_array *mask);

static bool
route1_fofc(void)
{
  const char *v = getenv("GKYL_ROUTE1_FOFC");
  return v && atoi(v) != 0;
}

// True when cell q fails the validator's test (mixed_gr_euler_state_valid, same order of operations); marker = its
// density marker (NAN when not reached).
static bool
route1_fofc_bad(const struct moment_species *fsp, const double *q, double rho_atm, double *marker)
{
  *marker = NAN;
  for (int m = 0; m < fsp->num_equations; ++m) if (!isfinite(q[m])) return true;
  if (!(q[27] >= pow(10.0, -8.0))) return false;
  if (!(q[0] > 0.0)) return true;
  double g[3][3] = { { q[9], q[10], q[11] }, { q[12], q[13], q[14] }, { q[15], q[16], q[17] } }, gi[3][3];
  double det = sym3_inverse(g, gi);
  if (!(isfinite(det) && det > 0.0)) return true;
  double S[3] = { q[1], q[2], q[3] }, de = q[0] + q[4], s2 = 0.0;
  for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) s2 += gi[i][j]*S[i]*S[j];
  double cone = mixed_gr_euler_momentum_capacity(q) - s2;
  *marker = mixed_gr_euler_density_marker(fsp, q, sqrt(det));
  bool in_band = *marker <= mixed_atm_band*rho_atm;
#ifdef ROUTE1_BARO
  if (!(cone > 0.0 || in_band)) in_band = gkyl_gr_euler_c2p_baro_handled(fsp->gr_euler_gas_gamma, q); // [BARO] the recovery's own test
#endif
  return !(isfinite(de) && isfinite(s2) && isfinite(cone) && de > 0.0 && (cone > 0.0 || in_band));
}

// Marks (mask = 1) the local cells of f that fail the test; returns how many were newly marked. Logs one [FOFC] line
// when any cell fails: failing count, newly marked, radius range, density marker / rho_atm range.
static int
route1_fofc_mark(gkyl_moment_app *app, const struct moment_species *fsp, const struct gkyl_array *f,
  struct gkyl_array *mask, double t, int stage, int pass)
{
  double rho_atm = container_of(fsp->equation, struct wv_gr_euler, eqn)->rho_atm;
  int nbad = 0, nnew = 0;
  double rmin = INFINITY, rmax = 0.0, mmin = INFINITY, mmax = 0.0;
  struct gkyl_range_iter it;
  gkyl_range_iter_init(&it, &app->local);
  while (gkyl_range_iter_next(&it)) {
    long l = gkyl_range_idx(&app->local_ext, it.idx);
    double mk;
    if (route1_fx_masked(l)) continue; // [FX] excised: not an evolved fluid cell
    if (!route1_fofc_bad(fsp, gkyl_array_cfetch(f, l), rho_atm, &mk)) continue;
    double *w = gkyl_array_fetch(mask, l), xc[3] = { 0, 0, 0 }, r = 0.0;
    nbad += 1;
    if (!(w[0] > 0.0)) { w[0] = 1.0; nnew += 1; }
    gkyl_rect_grid_cell_center(&app->grid, it.idx, xc);
    for (int d = 0; d < app->ndim; ++d) r += xc[d]*xc[d];
    r = sqrt(r);
    rmin = fmin(rmin, r); rmax = fmax(rmax, r);
    if (isfinite(mk)) { mmin = fmin(mmin, mk/rho_atm); mmax = fmax(mmax, mk/rho_atm); }
  }
  if (nbad)
    gkyl_moment_app_cout(app, stderr, "[FOFC] t=%.8g stage=%d pass=%d bad=%d new=%d r=[%.3f,%.3f]"
      " marker/rho_atm=[%.3e,%.3e]\n", t, stage, pass, nbad, nnew, rmin, rmax, mmin, mmax);
  return nnew;
}

// One forward-Euler stage with the correction. Stage 1 keeps a CFL-shortened length (as production); for stages 2
// and 3 a shortened length returns at once, as production then restarts the RK3. Cells that still fail after the
// redos show as a [FOFC] line with new=0 (the validator then rejects the step, as in production).
static void
route1_fofc_fe(gkyl_moment_app *app, struct moment_species *fsp, double t, double dt,
  const struct gkyl_array *fin, struct gkyl_array *fout, struct gkyl_update_status *st, struct gkyl_array *mask,
  int stage)
{
  route1_fluid_fe(app, fsp, t, dt, fin, fout, st);
  bool marked = false;
  for (int pass = 0; ; ++pass) {
    if (stage > 1 && st->dt_actual < dt) break;
    if (route1_fofc_mark(app, fsp, fout, mask, t, stage, pass) == 0 || pass == 3) break;
    marked = true;
    double h = st->dt_actual;
    gkyl_mp_scheme_set_lo_mask(fsp->mp_slvr, mask);
    route1_fluid_fe(app, fsp, t, h, fin, fout, st);
    gkyl_mp_scheme_set_lo_mask(fsp->mp_slvr, NULL);
    if (stage == 1) dt = h;
  }
  if (marked) gkyl_array_clear(mask, 0.0);
}

static struct gkyl_update_status
route1_fluid_fofc(gkyl_moment_app *app, int fidx, double t, double h)
{
  struct moment_species *fsp = &app->species[fidx];
  static struct gkyl_array *mask;
  if (!mask) {
    mask = gkyl_array_new(GKYL_DOUBLE, 1, app->local_ext.volume);
    gkyl_array_clear(mask, 0.0);
    gkyl_moment_app_cout(app, stderr, "[FOFC] on: fluid meqn = %d (first-order faces through the surface-U1 path when"
      " 76)\n", fsp->num_equations);
  }
  struct gkyl_array *f0 = fsp->f0, *f1 = fsp->f1, *fn = fsp->fnew;
  const struct gkyl_range *r = &app->local_ext;
  struct gkyl_update_status st = { .success = true };
  route1_fx_update(app, fsp, t); // [FX] mask growth (legal only) and boundary legality check
  for (;;) {
    route1_fofc_fe(app, fsp, t, h, f0, f1, &st, mask, 1);
    route1_fx_stage[0] = route1_fx_res;
    h = st.dt_actual;
    route1_fofc_fe(app, fsp, t + h, h, f1, fn, &st, mask, 2);
    route1_fx_stage[1] = route1_fx_res;
    if (st.dt_actual < h) { h = st.dt_actual; continue; }
    gkyl_array_accumulate_range(gkyl_array_set_range(f1, 0.75, f0, r), 0.25, fn, r);
    route1_fx_keep(app, f1, f0);
    route1_fofc_fe(app, fsp, t + 0.5*h, h, f1, fn, &st, mask, 3);
    route1_fx_stage[2] = route1_fx_res;
    if (st.dt_actual < h) { h = st.dt_actual; continue; }
    gkyl_array_accumulate_range(gkyl_array_set_range(f1, 1.0/3.0, f0, r), 2.0/3.0, fn, r);
    route1_fx_keep(app, f1, f0);
    gkyl_array_copy_range(f0, f1, r);
    // [FX] D into the mask over the RK3: u0 + h (L0/6 + L1/6 + 2 L2/3)
    route1_fx_pend += route1_fx_stage[0]/6.0 + route1_fx_stage[1]/6.0 + 2.0*route1_fx_stage[2]/3.0;
    break;
  }
  st.dt_actual = h;
  return st;
}

// One fluid MP SSP-RK3 update of length h at t (a half-step of the mixed step), Strang-split with the exact boost
// under GKYL_ROUTE1_BOOST (see route1_boost_enabled); otherwise production's update. When the CFL shortens the
// update (the RK3 restarts from its input at the shorter length), it is redone from the pre-boost state at the length
// taken, so the two boost halves always match the RK3; a last shortening (after 3 redos) is closed by boosting the
// remainder, so the total boost always equals the length taken.
static struct gkyl_update_status
route1_fluid_rk3(gkyl_moment_app *app, int fidx, double t, double h)
{
  struct moment_species *fsp = &app->species[fidx];
  bool wb = fsp->num_equations > 71 && !gkyl_gr_euler_wb_disabled(fsp->equation);
  if (route1_fofc() && fsp->has_gr_euler && fsp->scheme_type == GKYL_MOMENT_MP) { // WB too: U1 of the deviation
    int n = 0;
    for (int i = 0; i < app->num_species; ++i) n += app->species[i].scheme_type != GKYL_MOMENT_WAVE_PROP;
    if (n == 1 && !app->has_field && !route1_boost_enabled()) return route1_fluid_fofc(app, fidx, t, h);
  }
  if (!route1_boost_enabled() || !fsp->has_gr_euler || wb) return moment_update_ssp_rk3_at(app, t, h);
  if (route1_boost_lawson()) {
    int n = 0;
    for (int i = 0; i < app->num_species; ++i) n += app->species[i].scheme_type != GKYL_MOMENT_WAVE_PROP;
    if (n == 1 && !app->has_field) return route1_fluid_lawson(app, fidx, t, h);
  }
  static struct gkyl_array *f0, *da;
  long vol = app->local_ext.volume;
  if (!f0) {
    f0 = gkyl_array_new(GKYL_DOUBLE, fsp->num_equations, vol);
    da = gkyl_array_new(GKYL_DOUBLE, 3, vol);
  }
  gkyl_array_copy_range(f0, fsp->fcurr, &app->local_ext);
  struct gkyl_update_status st;
  for (int k = 0; ; ++k) {
    if (k) gkyl_array_copy_range(fsp->fcurr, f0, &app->local_ext);
    route1_fluid_boost(app, fsp, 0.5*h);
    moment_species_apply_bc(app, t, fsp, fsp->fcurr);
    for (long l = 0; l < vol; ++l) {
      double *q = gkyl_array_fetch(fsp->fcurr, l), *z = gkyl_array_fetch(da, l);
      for (int c = 0; c < 3; ++c) { z[c] = q[28 + c]; q[28 + c] = 0.0; }
    }
    st = moment_update_ssp_rk3_at(app, t, h);
    for (long l = 0; l < vol; ++l) {
      double *q = gkyl_array_fetch(fsp->fcurr, l);
      const double *z = gkyl_array_cfetch(da, l);
      for (int c = 0; c < 3; ++c) q[28 + c] = z[c];
    }
    if (!st.success || st.dt_actual == h || k == 3) break;
    route1_lg.nshort += 1;
    h = st.dt_actual;
  }
  if (st.success) route1_fluid_boost(app, fsp, st.dt_actual - 0.5*h);
  return st;
}

// The validator under GKYL_ROUTE1_FLOOR_ONLY (see route1_boost_enabled); counts the cells it accepts outside the
// cone by call site (after fluid update 1, after the metric copy, after fluid update 2).
static bool
route1_state_valid_floor(gkyl_moment_app *app, const struct moment_species *fsp)
{
  if (!fsp->has_gr_euler) return true;
  const struct wv_gr_euler *gr = container_of(fsp->equation, struct wv_gr_euler, eqn);
  const char *fe = getenv("GKYL_ATM_CONS_RESET_FACTOR");
  double band = (fe && atof(fe) > 1.0 ? atof(fe) : 100.0)*gr->rho_atm;
  int bad = 0;
  long acc = 0;
  struct gkyl_range_iter it;
  gkyl_range_iter_init(&it, &app->local);
  while (!bad && gkyl_range_iter_next(&it)) {
    const double *q = gkyl_array_cfetch(fsp->fcurr, gkyl_range_idx(&app->local_ext, it.idx));
    for (int m = 0; m < fsp->num_equations; ++m) bad = bad || !isfinite(q[m]);
    if (bad || q[27] < 1e-8) continue;
    double g[3][3] = { { q[9], q[10], q[11] }, { q[12], q[13], q[14] }, { q[15], q[16], q[17] } }, gi[3][3];
    double det = sym3_inverse(g, gi);
    if (!(isfinite(det) && det > 0.0)) { bad = 1; continue; }
    double s2 = 0.0;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) s2 += gi[i][j]*q[1 + i]*q[1 + j];
    if (q[0] > 0.0 && q[0] + q[4] > 0.0 && mixed_gr_euler_momentum_capacity(q) - s2 > 0.0) continue;
    if (mixed_gr_euler_density_marker(fsp, q, sqrt(det)) <= band) acc += 1;
    else bad = 1;
  }
  int gbad = bad;
  gkyl_comm_allreduce(app->comm, GKYL_INT, GKYL_MAX, 1, &bad, &gbad);
  if (route1_lg.init) route1_lg.ta[route1_lg.nflux >= 2 ? 2 : (route1_lg.nreset >= 2 ? 1 : 0)] += acc;
  return gbad == 0;
}
