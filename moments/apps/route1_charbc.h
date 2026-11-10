// Route-1 characteristic boundary flux (GKYL_ROUTE1_CHARBC=1). Included by route1_step.h.
//
// The outer boundary is a flux. Every derivative variable is a flux variable, every primary is pointwise:
//   face-flux rows: K, D, A, V (library flux), B_dj = -1/2 R_beta^j (the Gamma-driver rate of beta^j), P_d = -S_pure;
//   pointwise rows: gt, ln alpha, beta, chi (sources only);  H: grid-defined, sym D0 L(P), one-sided on the boundary cells.
// Interior faces carry the central flux 1/2 (F(q_i) + F(q_i+1)), so in the interior B's rate is 1/2 D0 R_beta and P's is
// D0 S_pure, as before; KO now also acts on B (in the interior KO(1/2 D0 beta) = 1/2 D0 KO(beta), so nothing changes there).
// On a non-periodic global face the flux is F(u_b) with u_b = e + Pi (q_c - e): q_c the boundary cell, e the held
// equilibrium (flat when there is none), Pi the projector onto the outgoing and zero-speed characteristic fields of the
// Jacobian of the face-flux vector at flat space, one per axis and side, computed once (LAPACK dgeev). Incoming fields are
// zero relative to e. ln alpha stands for alpha in q - e. Each boundary face of an edge or corner cell gets its own flux.
// No spacetime stencil reads a ghost cell (periodic directions excepted, where the ghosts are copies across the seam).
// Kidder et al., PRD 71 064020 (2005); Lindblom et al., CQG 23 S447 (2006).
// Env: GKYL_ROUTE1_CHARBC (0/1), GKYL_ROUTE1_CHARBC_DUMP (path: write the three 77x77 Jacobians, row-major, then go on).

extern int dgeev_(char *jobvl, char *jobvr, int *n, double *a, int *lda, double *wr, double *wi, double *vl, int *ldvl,
  double *vr, int *ldvr, double *work, int *lwork, int *info);
extern int dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
extern int dgetri_(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);

static bool
route1_charbc(void)
{
  const char *v = getenv("GKYL_ROUTE1_CHARBC");
  return v && atoi(v) != 0;
}

// Gamma-driver parameters of route1_rhs (production gating).
struct route1_drv {
  bool on;
  double mu, eta, alpha_floor;
};

// Face-flux vector along d: library rows 0-51 (alpha's row 9 is zero), B_d (55+3d+j) = -1/2 R_beta^j, P_d (65+d) =
// -(S_pure(q) - spe); every other row 0. ref: the cell's driver reference (Gammatilde 0-2, beta 3-5), NULL = zero.
static void
route1_face_flux(const struct moment_species *esp, int d, const double *q, const double *ref, double spe,
  const struct route1_drv *drv, double *f)
{
  double fl[77];
  route1_flux(esp, d, q, fl);
  for (int m = 0; m < 77; ++m) f[m] = (m < 52 && m != 9) ? fl[m] : 0.0;
  if (drv->on && gamma_driver_cell_active(q, drv->alpha_floor, true)) {
    double Gam[3];
    route1_gamma(q, Gam);
    for (int j = 0; j < 3; ++j)
      f[55 + 3*d + j] = -0.5*(drv->mu*(Gam[j] - (ref ? ref[j] : 0.0)) - drv->eta*(q[52 + j] - (ref ? ref[3 + j] : 0.0)));
  }
  f[65 + d] = -(route1_spure(esp, q) - spe);
}

static struct {
  bool ready;
  double pi[3][2][77*77]; // [axis][0 lower face, 1 upper face], row-major: (Pi x)_m = sum_n pi[m*77 + n] x_n
} route1_cb;

// Projectors of every axis from the flat-space Jacobian (once). Reports the speeds and the eigenvector condition number
// (check K1, NEWRAD_AUDIT.md 2026-11-05 13:30); refuses (exit 2) when the Jacobian is not diagonalizable.
static void
route1_charbc_setup(gkyl_moment_app *app, const struct moment_species *esp, const struct route1_drv *drv)
{
  if (route1_cb.ready) return;
  const int n = 77;
  double q0[77] = { 0.0 };
  q0[0] = q0[4] = q0[8] = 1.0; q0[9] = 1.0; q0[64] = 1.0;
  const char *dump = getenv("GKYL_ROUTE1_CHARBC_DUMP");
  FILE *fd = dump ? fopen(dump, "wb") : NULL;
  for (int d = 0; d < app->ndim; ++d) {
    static double J[77*77], A[77*77], VR[77*77], VI[77*77], wr[77], wi[77], work[77*64];
    int ipiv[77];
    for (int c = 0; c < n; ++c) { // column c: d F / d q_c (q_9 as ln alpha: alpha = 1 there, so the same)
      double qp[77], qm[77], fp[77], fm[77], eps = 1e-6;
      memcpy(qp, q0, sizeof qp); memcpy(qm, q0, sizeof qm);
      qp[c] += eps; qm[c] -= eps;
      route1_face_flux(esp, d, qp, NULL, route1_spure(esp, q0), drv, fp);
      route1_face_flux(esp, d, qm, NULL, route1_spure(esp, q0), drv, fm);
      for (int m = 0; m < n; ++m) J[m*n + c] = (fp[m] - fm[m])/(2*eps);
    }
    if (fd) fwrite(J, sizeof(double), n*n, fd);
    for (int m = 0; m < n; ++m) for (int c = 0; c < n; ++c) A[c*n + m] = J[m*n + c]; // column-major
    int N = n, lwork = 77*64, info = 0;
    char nv = 'N', vv = 'V';
    dgeev_(&nv, &vv, &N, A, &N, wr, wi, NULL, &N, VR, &N, work, &lwork, &info);
    double wimax = 0.0, vnorm = 0.0, inorm = 0.0;
    for (int m = 0; m < n; ++m) wimax = fmax(wimax, fabs(wi[m]));
    memcpy(VI, VR, sizeof VI);
    for (int m = 0; m < n*n; ++m) vnorm += VR[m]*VR[m];
    int linfo = 0;
    if (info == 0) { dgetrf_(&N, &N, VI, &N, ipiv, &linfo); if (!linfo) dgetri_(&N, VI, &N, ipiv, work, &lwork, &linfo); }
    for (int m = 0; m < n*n; ++m) inorm += VI[m]*VI[m];
    double cond = (info || linfo) ? INFINITY : sqrt(vnorm*inorm);
    // speeds, grouped
    double sp[77]; int ns = 0;
    for (int m = 0; m < n; ++m) {
      int k = 0;
      while (k < ns && fabs(sp[k] - wr[m]) > 1e-6) ++k;
      if (k == ns) sp[ns++] = wr[m];
    }
    char buf[2048]; int off = 0;
    for (int k = 0; k < ns; ++k) {
      int mult = 0;
      for (int m = 0; m < n; ++m) mult += fabs(wr[m] - sp[k]) <= 1e-6;
      off += snprintf(buf + off, sizeof buf - off, " %+.6f(x%d)", sp[k], mult);
    }
    gkyl_moment_app_cout(app, stderr, "[ROUTE1] charbc axis %d: speeds%s; max|Im| %.1e; eigenvector condition %.3e\n",
      d, buf, wimax, cond);
    if (!(cond < 1e8) || wimax > 1e-8) {
      gkyl_moment_app_cout(app, stderr, "[ROUTE1] charbc axis %d: Jacobian not diagonalizable (K1 fails); refusing\n", d);
      exit(2);
    }
    // Pi_side = VR diag(keep) VR^-1; lower face (outward -d) keeps lambda < tol, upper face keeps lambda > -tol.
    for (int side = 0; side < 2; ++side) {
      double *P = route1_cb.pi[d][side];
      for (int m = 0; m < n; ++m)
        for (int c = 0; c < n; ++c) {
          double s = 0.0;
          for (int k = 0; k < n; ++k) {
            bool keep = side ? wr[k] > -1e-8 : wr[k] < 1e-8;
            if (keep) s += VR[k*n + m]*VI[c*n + k]; // VR(m,k) * VR^-1(k,c), column-major
          }
          P[m*n + c] = s;
        }
    }
  }
  if (fd) { fclose(fd); gkyl_moment_app_cout(app, stderr, "[ROUTE1] charbc Jacobians written to %s\n", dump); }
  gkyl_moment_app_cout(app, stderr, "[ROUTE1] charbc: boundary faces carry F(e + Pi (q - e)); KO cut at the boundary,"
    " no ghost reads; B and P in flux form, H one-sided on the boundary cells\n");
  route1_cb.ready = true;
}

// u_b = e + Pi (q - e) on the side's face of axis d (q_9 - e_9 taken in ln alpha).
static void
route1_charbc_state(int d, int side, const double *q, const double *e, double *ub)
{
  static const double flat[77] = { [0] = 1.0, [4] = 1.0, [8] = 1.0, [9] = 1.0, [64] = 1.0 };
  const double *b = e ? e : flat;
  double w[77];
  for (int m = 0; m < 77; ++m) w[m] = m == 9 ? log(q[9]/b[9]) : q[m] - b[m];
  const double *P = route1_cb.pi[d][side];
  for (int m = 0; m < 77; ++m) {
    double s = 0.0;
    for (int c = 0; c < 77; ++c) s += P[m*77 + c]*w[c];
    ub[m] = m == 9 ? b[9]*exp(s) : b[m] + s;
  }
}

// Step 1 under the switch: rate -= (G_i+1/2 - G_i-1/2)/h on the face-flux rows, G the face-flux vector minus its value at e,
// central inside, F(u_b) - F(e) on non-periodic global faces. sp: S_pure(u) - S_pure(e) on local (and periodic ghosts).
static void
route1_charbc_divergence(gkyl_moment_app *app, const struct moment_species *esp, const struct gkyl_array *u,
  const struct gkyl_array *e, const struct gkyl_array *spe, const struct route1_drv *drv, struct gkyl_array *G,
  struct gkyl_array *rate)
{
  struct route1_work *w = &route1_w;
  int ndim = app->ndim, nc = esp->num_equations;
  const struct gkyl_array *ref = drv->on ? esp->gamma_driver_reference : NULL;
  struct gkyl_range_iter it;
  for (int d = 0; d < ndim; ++d) {
    bool per = route1_global_periodic(app, d);
    // cell flux G = F_d(q) - F_d(e) on local (and on the periodic ghosts along d)
    gkyl_range_iter_init(&it, &app->local_ext);
    while (gkyl_range_iter_next(&it)) {
      bool inside = true;
      for (int k = 0; k < ndim; ++k) {
        bool out = it.idx[k] < app->local.lower[k] || it.idx[k] > app->local.upper[k];
        if (out && !(k == d && per)) inside = false;
      }
      if (!inside) continue;
      long l = gkyl_range_idx(&app->local_ext, it.idx);
      double *g = gkyl_array_fetch(G, l), fe[77];
      const double *se = spe ? gkyl_array_cfetch(spe, l) : NULL;
      route1_face_flux(esp, d, gkyl_array_cfetch(u, l), ref ? gkyl_array_cfetch(ref, l) : NULL, se ? se[0] : 0.0, drv, g);
      if (e) {
        route1_face_flux(esp, d, gkyl_array_cfetch(e, l), ref ? gkyl_array_cfetch(ref, l) : NULL, se ? se[0] : 0.0, drv, fe);
        for (int m = 0; m < nc && m < 77; ++m) g[m] -= fe[m];
      }
    }
    double h = app->grid.dx[d];
    gkyl_range_iter_init(&it, &app->local);
    while (gkyl_range_iter_next(&it)) {
      long l = gkyl_range_idx(&app->local_ext, it.idx);
      double *r = gkyl_array_fetch(rate, l);
      const double *g0 = gkyl_array_cfetch(G, l);
      double fp[77], fm[77];
      for (int side = 0; side < 2; ++side) {
        double *ff = side ? fp : fm;
        bool wall = !per && it.idx[d] == (side ? app->global.upper[d] : app->global.lower[d]);
        if (wall) {
          const double *q = gkyl_array_cfetch(u, l), *qe = e ? gkyl_array_cfetch(e, l) : NULL;
          const double *rf = ref ? gkyl_array_cfetch(ref, l) : NULL;
          const double *se = spe ? gkyl_array_cfetch(spe, l) : NULL;
          double ub[77], fe[77];
          route1_charbc_state(d, side, q, qe, ub);
          route1_face_flux(esp, d, ub, rf, se ? se[0] : 0.0, drv, ff);
          if (qe) {
            route1_face_flux(esp, d, qe, rf, se ? se[0] : 0.0, drv, fe);
            for (int m = 0; m < 77; ++m) ff[m] -= fe[m];
          }
        }
        else {
          const double *gn = gkyl_array_cfetch(G, route1_nb(app, it.idx, d, side ? 1 : -1));
          for (int m = 0; m < 77; ++m) ff[m] = 0.5*(g0[m] + gn[m]);
        }
      }
      for (int m = 0; m < 77; ++m) {
        bool row = (m < 52 && m != 9) || (m >= 55 && m < 64) || (m >= 65 && m < 68);
        if (row) r[m] -= (fp[m] - fm[m])/h;
      }
    }
  }
  (void) w;
}

// d/dx_d of component c of a over local: centred inside (and across a periodic seam), one-sided second order inward on a
// non-periodic global boundary cell. Reads no non-periodic ghost.
static double
route1_charbc_d(const gkyl_moment_app *app, const struct gkyl_array *a, int c, const int *idx, int d)
{
  double h = app->grid.dx[d];
  bool per = route1_global_periodic(app, d);
  int t = 0;
  if (!per && idx[d] == app->global.lower[d]) t = 1;
  if (!per && idx[d] == app->global.upper[d]) t = -1;
  const double *f0 = gkyl_array_cfetch(a, route1_nb(app, idx, d, 0));
  if (!t) {
    const double *fp = gkyl_array_cfetch(a, route1_nb(app, idx, d, 1)), *fm = gkyl_array_cfetch(a, route1_nb(app, idx, d, -1));
    return (fp[c] - fm[c])/(2*h);
  }
  const double *f1 = gkyl_array_cfetch(a, route1_nb(app, idx, d, t)), *f2 = gkyl_array_cfetch(a, route1_nb(app, idx, d, 2*t));
  return t*(-3*f0[c] + 4*f1[c] - f2[c])/(2*h);
}
