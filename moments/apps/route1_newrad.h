// Route-1 radiative outer boundary (GKYL_ROUTE1_NEWRAD=1). Included by route1_step.h.
//
// The ghost cells are the boundary points of ETK NewRad: they are evolved in the same SSP-RK3 stages as
// the interior, with a radiative rate in place of the PDE. The interior update stays the central
// face-flux difference, so the boundary enters only through the single-valued boundary-face flux and the
// update stays conservative. Production NewRad instead overwrites the outermost interior layer after the
// step by one forward-Euler update over dt.
//
// Per component, with w = q - e (held equilibrium; flat values when there is none), n = x/r
// (Alcubierre et al. gr-qc/0206072 Eqs. 108-113; ETK NewRad newrad.cc):
//   d_t w = -v n.grad w - v w/r + aux (r_p/r)^N,  aux = actual rate at p + v n.grad w(p) + v w(p)/r_p,
//   p = one cell inward along every boundary direction (diagonal at edges and corners), processed after
//   the ghost layer it belongs to. Derivatives: one-sided second order inward in boundary directions,
//   centered otherwise, at the boundary point and at p alike (newrad.cc).
// Speeds (ETK ML_BSSN: v0 = sqrt(harmonicF) for alpha, trK and the lapse driver, 1 otherwise):
//   alpha, A_k: v_g = sqrt(f(alpha = 1)) (1+log: sqrt 2; harmonic, geodesic: 1)
//   K_ij: the trace-free part at 1 and trK = chi^2 gt^ij K_ij at v_g
//   gt, D, V, beta, chi: 1
// Grid-defined rows follow their parents as in the interior: P = D (chi rate), H = sym D (P rate),
// B = 1/2 D (beta rate), D centered where both neighbours exist, one-sided inward at the outer layer.
// Env: GKYL_ROUTE1_NEWRAD (0/1), GKYL_ROUTE1_NEWRAD_N (N, default 3 as the paper; < 0: no correction),
//   GKYL_ROUTE1_NEWRAD_VG (override v_g), GKYL_ROUTE1_NEWRAD_OLDSPEED (1: production speeds, every K_ij
//   at 1 and alpha, A_k at sqrt(2/alpha); for the speed A/B test only), GKYL_ROUTE1_NEWRAD_DPOW (falloff
//   power p of the derivative variables D, A, V in w = u(r - vt)/r^p, i.e. d_t w = -v n.grad w - p v w/r;
//   default 1 = the paper's law for every field; 2 = one power faster than the fields they differentiate).
// GKYL_ROUTE1_PEEL=1 (implies NewRad): radiate the parents, derive the derivatives. NewRad evolves only gt,
//   ln alpha, K, V, beta and chi on the ghosts; after every stage D, A, P, B and H there are set from the ghost
//   states with the interior's identities (route1_peel_derive). KO reaches into the ghosts in every direction
//   and B is centred on the boundary cells (route1_step.h). No band, no buffer.

static bool
route1_peel(void)
{
  const char *v = getenv("GKYL_ROUTE1_PEEL");
  return v && atoi(v) != 0;
}

static bool
route1_newrad_enabled(void)
{
  const char *v = getenv("GKYL_ROUTE1_NEWRAD");
  return (v && atoi(v) != 0) || route1_peel();
}

static double
route1_newrad_power(void)
{
  const char *v = getenv("GKYL_ROUTE1_NEWRAD_N");
  return v ? atof(v) : 3.0;
}

static bool
route1_newrad_oldspeed(void)
{
  const char *v = getenv("GKYL_ROUTE1_NEWRAD_OLDSPEED");
  return v && atoi(v) != 0;
}

static double
route1_newrad_dpow(void)
{
  const char *v = getenv("GKYL_ROUTE1_NEWRAD_DPOW");
  return v ? atof(v) : 1.0;
}

// GKYL_ROUTE1_NEWRAD_BETA_SPLIT=1: the shift leaves at the integrated driver's own speeds, sqrt(4 mu/3) along
// n = x/r and sqrt(mu) across it (mu from GKYL_GAMMA_DRIVER, 0.75 if unset or 0). Returns mu, or 0 when off
// (every beta component then leaves at 1).
static double
route1_newrad_beta_mu(void)
{
  const char *v = getenv("GKYL_ROUTE1_NEWRAD_BETA_SPLIT");
  if (!(v && atoi(v) != 0)) return 0.0;
  const char *m = getenv("GKYL_GAMMA_DRIVER");
  double mu = m ? atof(m) : 0.0;
  return mu > 0.0 ? mu : 0.75;
}

static double
route1_newrad_vgauge(enum gkyl_spacetime_slicing slicing)
{
  const char *v = getenv("GKYL_ROUTE1_NEWRAD_VG");
  if (v) return atof(v);
  return slicing == GKYL_1PLUSLOG_SLICING ? sqrt(2.0) : 1.0;
}

// Boundary directions s (+1 above upper, -1 below lower, 0 inside) and ghost layer (0 = interior).
// Periodic directions are not boundaries: their ghosts count as interior here (filled across the seam).
static int
route1_ghost_dir(const gkyl_moment_app *app, const int *idx, int s[3])
{
  int layer = 0;
  for (int d = 0; d < 3; ++d) {
    s[d] = 0;
    if (d >= app->ndim || route1_global_periodic(app, d)) continue;
    int o = 0;
    if (idx[d] > app->local.upper[d]) { s[d] = 1; o = idx[d] - app->local.upper[d]; }
    else if (idx[d] < app->local.lower[d]) { s[d] = -1; o = app->local.lower[d] - idx[d]; }
    if (o > layer) layer = o;
  }
  return layer;
}

static double
route1_nr_flat(int m)
{
  return (m == 0 || m == 4 || m == 8 || m == 9 || m == 64) ? 1.0 : 0.0;
}

static double
route1_trk(const double *q)
{
  double gt[3][3], gti[3][3], tr = 0.0;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) gt[i][j] = q[3*i + j];
  sym3_inverse(gt, gti);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) tr += gti[i][j]*q[10 + 3*i + j];
  return q[64]*q[64]*tr;
}

// d/dt trK from the rates of chi, gt and K at one cell (rate[9] is not used).
static double
route1_trk_rate(const double *q, const double *rt)
{
  double gt[3][3], gti[3][3], dgt[3][3] = { { 0 } };
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) gt[i][j] = q[3*i + j];
  sym3_inverse(gt, gti);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b) dgt[i][j] -= gti[i][a]*rt[3*a + b]*gti[b][j];
  double chi = q[64], trk = 0.0, out = 0.0;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      trk += gti[i][j]*q[10 + 3*i + j];
      out += chi*chi*(dgt[i][j]*q[10 + 3*i + j] + gti[i][j]*rt[10 + 3*i + j]);
    }
  return out + 2*chi*rt[64]*trk;
}

// Scalar field radiated by NewRad: component m of q - e, or trK(q) - trK(e) for m < 0.
struct route1_nr_field { const struct gkyl_array *u, *e; int m; };

static double
route1_nr_val(const struct route1_nr_field *f, long l)
{
  const double *q = gkyl_array_cfetch(f->u, l);
  const double *qe = f->e ? gkyl_array_cfetch(f->e, l) : NULL;
  if (f->m < 0) return route1_trk(q) - (qe ? route1_trk(qe) : 0.0);
  return q[f->m] - (qe ? qe[f->m] : route1_nr_flat(f->m));
}

static double
route1_radius(const gkyl_moment_app *app, const int *idx)
{
  double x[3] = { 0.0, 0.0, 0.0 };
  gkyl_rect_grid_cell_center(&app->grid, idx, x);
  return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}

// v n.grad w + pw v w/r at cell p, with the stencil orientation s of the boundary point (pw: falloff power).
static double
route1_nr_sommerfeld(const gkyl_moment_app *app, const struct route1_nr_field *f, const int *p,
  const int *s, double v, double pw)
{
  double x[3] = { 0.0, 0.0, 0.0 };
  gkyl_rect_grid_cell_center(&app->grid, p, x);
  double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
  double w0 = route1_nr_val(f, gkyl_range_idx(&app->local_ext, p)), adv = 0.0;
  for (int d = 0; d < app->ndim; ++d) {
    double h = app->grid.dx[d], der;
    if (s[d]) {
      double w1 = route1_nr_val(f, route1_nb(app, p, d, -s[d]));
      double w2 = route1_nr_val(f, route1_nb(app, p, d, -2*s[d]));
      der = s[d]*(3*w0 - 4*w1 + w2)/(2*h);
    }
    else {
      der = (route1_nr_val(f, route1_nb(app, p, d, 1)) - route1_nr_val(f, route1_nb(app, p, d, -1)))/(2*h);
    }
    adv += x[d]/r*der;
  }
  return v*adv + pw*v*w0/r;
}

// Derivative of component c of a along d: centered when both neighbours lie in local_ext, otherwise
// one-sided second order inward.
static double
route1_dr(const gkyl_moment_app *app, const struct gkyl_array *a, int c, const int *idx, int d)
{
  double h = app->grid.dx[d];
  const double *f0 = gkyl_array_cfetch(a, route1_nb(app, idx, d, 0));
  bool lo = route1_in_ext(app, idx, d, -1), hi = route1_in_ext(app, idx, d, 1);
  if (lo && hi) {
    const double *fp = gkyl_array_cfetch(a, route1_nb(app, idx, d, 1));
    const double *fm = gkyl_array_cfetch(a, route1_nb(app, idx, d, -1));
    return (fp[c] - fm[c])/(2*h);
  }
  int t = hi ? -1 : 1;
  const double *f1 = gkyl_array_cfetch(a, route1_nb(app, idx, d, -t));
  const double *f2 = gkyl_array_cfetch(a, route1_nb(app, idx, d, -2*t));
  return t*(3*f0[c] - 4*f1[c] + f2[c])/(2*h);
}

// Derivative of component m of a state along d with route1_dr's stencil; m = 9 differentiates ln alpha.
static double
route1_peel_dr(const gkyl_moment_app *app, const struct gkyl_array *u, int m, const int *idx, int d)
{
  #define ROUTE1_PV(o) (m == 9 ? log(((const double *) gkyl_array_cfetch(u, route1_nb(app, idx, d, o)))[9]) \
    : ((const double *) gkyl_array_cfetch(u, route1_nb(app, idx, d, o)))[m])
  double h = app->grid.dx[d];
  bool lo = route1_in_ext(app, idx, d, -1), hi = route1_in_ext(app, idx, d, 1);
  double out;
  if (lo && hi) out = (ROUTE1_PV(1) - ROUTE1_PV(-1))/(2*h);
  else {
    int t = hi ? -1 : 1;
    out = t*(3*ROUTE1_PV(0) - 4*ROUTE1_PV(-t) + ROUTE1_PV(-2*t))/(2*h);
  }
  #undef ROUTE1_PV
  return out;
}

// GKYL_ROUTE1_PEEL: the grid-defined rows on every ghost cell (faces, edges, corners) from the ghost states, with
// the identities the interior keeps (route1_residuals): D_kij = 1/2 D_k gt_ij, A_k = D_k ln alpha, P_k = D_k chi,
// B_kj = 1/2 D_k beta^j, then H_kj = 1/2 (D_k P_j + D_j P_k), which reads P on the neighbouring ghosts.
static void
route1_peel_derive(const gkyl_moment_app *app, struct gkyl_array *u)
{
  int ndim = app->ndim;
  struct gkyl_range_iter it;
  for (int pass = 0; pass < 2; ++pass) {
    gkyl_range_iter_init(&it, &app->local_ext);
    while (gkyl_range_iter_next(&it)) {
      int s[3];
      if (route1_ghost_dir(app, it.idx, s) == 0) continue;
      double *q = gkyl_array_fetch(u, gkyl_range_idx(&app->local_ext, it.idx));
      for (int k = 0; k < 3; ++k) {
        bool on = k < ndim;
        if (pass == 0) {
          for (int ij = 0; ij < 9; ++ij) q[19 + 9*k + ij] = on ? 0.5*route1_peel_dr(app, u, ij, it.idx, k) : 0.0;
          q[46 + k] = on ? route1_peel_dr(app, u, 9, it.idx, k) : 0.0;
          q[65 + k] = on ? route1_peel_dr(app, u, 64, it.idx, k) : 0.0;
          for (int j = 0; j < 3; ++j) q[55 + 3*k + j] = on ? 0.5*route1_peel_dr(app, u, 52 + j, it.idx, k) : 0.0;
        }
        else {
          for (int j = 0; j < 3; ++j)
            q[68 + 3*k + j] = 0.5*((on ? route1_peel_dr(app, u, 65 + j, it.idx, k) : 0.0)
              + (j < ndim ? route1_peel_dr(app, u, 65 + k, it.idx, j) : 0.0));
        }
      }
    }
  }
}

// Radiative rates on every ghost cell (rows 0..54 except B, and chi). Interior rates must be complete
// for rows 0..54 and 64. Slot 9 of rate is d(ln alpha)/dt, as in the interior. Also sets sp on the ghosts
// to their chi rate, so P = D (chi rate) holds across the boundary.
static void
route1_newrad_rates(gkyl_moment_app *app, struct moment_species *esp, const struct gkyl_array *u,
  const struct gkyl_array *e, struct gkyl_array *rate, struct gkyl_array *sp)
{
  enum gkyl_spacetime_slicing slicing = esp->vacuum_einstein_conformal_spacetime_slicing;
  double npow = route1_newrad_power(), vg = route1_newrad_vgauge(slicing), dpow = route1_newrad_dpow();
  double bmu = route1_newrad_beta_mu();
  bool old = route1_newrad_oldspeed(), peel = route1_peel();
  struct gkyl_range_iter it;
  for (int layer = 1; layer <= 2; ++layer) {
    gkyl_range_iter_init(&it, &app->local_ext);
    while (gkyl_range_iter_next(&it)) {
      int s[3], ip[GKYL_MAX_DIM];
      if (route1_ghost_dir(app, it.idx, s) != layer) continue;
      for (int d = 0; d < app->ndim; ++d) ip[d] = it.idx[d] - s[d];
      double r = route1_radius(app, it.idx), rp = route1_radius(app, ip);
      double scale = npow >= 0 ? pow(rp/r, npow) : 0.0;
      long l = gkyl_range_idx(&app->local_ext, it.idx), lpi = gkyl_range_idx(&app->local_ext, ip);
      const double *q = gkyl_array_cfetch(u, l), *qp = gkyl_array_cfetch(u, lpi);
      const double *rtp = gkyl_array_cfetch(rate, lpi);
      double *rt = gkyl_array_fetch(rate, l);
      for (int m = 0; m < 65; ++m) {
        if (m >= 55 && m < 64) continue;
        if (peel && m >= 19 && m <= 48) continue; // D, A: derived after the stage (route1_peel_derive)
        if (bmu > 0 && m >= 52 && m <= 54) continue; // split below
        bool gauge = m == 9 || (m >= 46 && m <= 48);
        double v = gauge ? (old ? newrad_speed(m, q[9], slicing) : vg) : 1.0;
        double pw = (m >= 19 && m <= 51) ? dpow : 1.0; // D, A, V
        struct route1_nr_field f = { u, e, m };
        double act = m == 9 ? qp[9]*rtp[9] : rtp[m];
        double aux = act + route1_nr_sommerfeld(app, &f, ip, s, v, pw);
        double rr = -route1_nr_sommerfeld(app, &f, it.idx, s, v, pw) + aux*scale;
        rt[m] = m == 9 ? rr/q[9] : rr;
      }
      if (bmu > 0) {
        // beta: the radiative rate at the longitudinal and the transverse speed, projected along and across n.
        double x[3] = { 0.0, 0.0, 0.0 }, n[3], rl[3], rtr[3], nl = 0.0, nt = 0.0;
        double vl = sqrt(4.0*bmu/3.0), vt = sqrt(bmu);
        gkyl_rect_grid_cell_center(&app->grid, it.idx, x);
        for (int j = 0; j < 3; ++j) {
          n[j] = x[j]/r;
          struct route1_nr_field f = { u, e, 52 + j };
          rl[j] = -route1_nr_sommerfeld(app, &f, it.idx, s, vl, 1.0)
            + (rtp[52 + j] + route1_nr_sommerfeld(app, &f, ip, s, vl, 1.0))*scale;
          rtr[j] = -route1_nr_sommerfeld(app, &f, it.idx, s, vt, 1.0)
            + (rtp[52 + j] + route1_nr_sommerfeld(app, &f, ip, s, vt, 1.0))*scale;
        }
        for (int j = 0; j < 3; ++j) { nl += n[j]*rl[j]; nt += n[j]*rtr[j]; }
        for (int j = 0; j < 3; ++j) rt[52 + j] = n[j]*nl + rtr[j] - n[j]*nt;
      }
      if (!old) {
        // Replace the trace of the speed-1 K rates by the radiated trK rate.
        struct route1_nr_field f = { u, e, -1 };
        double aux = route1_trk_rate(qp, rtp) + route1_nr_sommerfeld(app, &f, ip, s, vg, 1.0);
        double rk = -route1_nr_sommerfeld(app, &f, it.idx, s, vg, 1.0) + aux*scale;
        double gt[3][3], gti[3][3], chi2 = q[64]*q[64], tr1 = 0.0;
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j) gt[i][j] = q[3*i + j];
        sym3_inverse(gt, gti);
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j) tr1 += chi2*gti[i][j]*rt[10 + 3*i + j];
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j) rt[10 + 3*i + j] += gt[i][j]/chi2*(rk - tr1)/3.0;
      }
      ((double *) gkyl_array_fetch(sp, l))[0] = rt[64];
    }
  }
}

// Grid-defined rows on the ghosts, after route1_newrad_rates and the P-rate pass (lp on all of local_ext).
static void
route1_newrad_derived(gkyl_moment_app *app, const struct gkyl_array *lp, struct gkyl_array *rate)
{
  struct gkyl_range_iter it;
  gkyl_range_iter_init(&it, &app->local_ext);
  while (gkyl_range_iter_next(&it)) {
    int s[3];
    if (route1_ghost_dir(app, it.idx, s) == 0) continue;
    double *rt = gkyl_array_fetch(rate, gkyl_range_idx(&app->local_ext, it.idx));
    const double *lpc = gkyl_array_cfetch(lp, gkyl_range_idx(&app->local_ext, it.idx));
    double dlp[3][3] = { { 0 } };
    for (int k = 0; k < 3; ++k) rt[65 + k] = lpc[k];
    for (int j = 0; j < app->ndim; ++j)
      for (int k = 0; k < 3; ++k) dlp[j][k] = route1_dr(app, lp, k, it.idx, j);
    for (int k = 0; k < 3; ++k)
      for (int j = 0; j < 3; ++j) rt[68 + 3*k + j] = 0.5*(dlp[k][j] + dlp[j][k]);
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        rt[55 + 3*i + j] = i < app->ndim ? 0.5*route1_dr(app, rate, 52 + j, it.idx, i) : 0.0;
  }
}

// Ghost states persist across steps (other code may refill ghosts between steps). The first call
// captures the BC-filled ghosts as the boundary points' initial data.
static void
route1_newrad_ghosts(gkyl_moment_app *app, struct gkyl_array *u, bool save)
{
  struct route1_work *w = &route1_w;
  if (!w->gh) {
    if (route1_peel()) route1_peel_derive(app, u);
    w->gh = gkyl_array_new(GKYL_DOUBLE, u->ncomp, app->local_ext.volume);
    gkyl_array_copy_range(w->gh, u, &app->local_ext);
    return;
  }
  struct gkyl_range_iter it;
  gkyl_range_iter_init(&it, &app->local_ext);
  while (gkyl_range_iter_next(&it)) {
    int s[3];
    if (route1_ghost_dir(app, it.idx, s) == 0) continue;
    long l = gkyl_range_idx(&app->local_ext, it.idx);
    if (save) memcpy(gkyl_array_fetch(w->gh, l), gkyl_array_cfetch(u, l), u->esznc);
    else memcpy(gkyl_array_fetch(u, l), gkyl_array_cfetch(w->gh, l), u->esznc);
  }
}
