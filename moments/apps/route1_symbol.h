// [SYMBOL] principal-symbol check along radial lines (diagnostic; GKYL_ROUTE1_SYMBOL=1 prints at the first route-1
// step, then exits 0). Included by route1_step.h after route1_charbc.h.
//
// Per cell and axis d: the Jacobian of the face-flux vector (route1_face_flux: library rows, B rows -1/2 R_beta, P rows
// -S_pure) at the cell's state, restricted to the derivative variables (K, D, A, V, B, P); primaries (gt, alpha, beta,
// chi) enter only at lower order. Plus the principal part the K source carries through the chi Hessian: H = sym d(P),
// so d_t K_ij gets sum_ab s_ij,ab H_ab with s = d S_K / d H (library source, numerically), i.e. a flux entry
// -1/2 (s_ij,db + s_ij,bd) in column P_b along d (NEWRAD_AUDIT.md, "Result, K1", item 2).
// Driver variants: mu = 0.75, mu = 0.6, mu = 0.5 alpha^2 chi^2 (weighted). Reported: eigenvector condition number of the
// 54x54 symbol (Frobenius norms, LAPACK unit columns), max |Im lambda|, and the positive speeds along the line's axis.
// Lines: +x from the centre (radial axis 0) and the diagonal i = j = k (axes 0, 1, 2 all oblique).

static bool
route1_symbol(void)
{
  const char *v = getenv("GKYL_ROUTE1_SYMBOL");
  return v && atoi(v) != 0;
}

static const int route1_sym_var[54] = { 10, 11, 12, 13, 14, 15, 16, 17, 18,
  19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45,
  46, 47, 48, 49, 50, 51, 55, 56, 57, 58, 59, 60, 61, 62, 63, 65, 66, 67 };

// [SYM2] variants (GKYL_ROUTE1_SYMBOL=2): mode 0 as built; 1 adds the shift advection beta^k d_k beta^j = 2 beta^k B_kj
// (B_dj flux -beta^k B_kj, B following the beta rate); 2 writes the advection of the gradient families D_kij, A_k, B_kj,
// P_k in advective form beta^l d_l X_k (flux -beta^l X_k along l) in place of the flux form d_k(beta^l X_l), i.e. adds the
// curl ordering constraint times beta, as first-order systems that are hyperbolic for any shift do (cf. Lindblom et al.
// 2006, FO generalized harmonic); mode 2 includes the shift advection. Out: the smallest |lambda| above 1e-7 and the
// count of |lambda| <= 1e-7 (zero-speed modes) of the last call.
static int route1_sym_mode = 0, route1_sym_verify = 0, route1_sym_nzero;
static double route1_sym_minabs;

// cond, max|Im|, positive speeds (sorted, distinct to 1e-4) of the symbol of axis d at q for driver strength mu.
static double
route1_symbol_cell(const struct moment_species *esp, int d, const double *q, double mu, double *wimax, double *spd,
  int *nspd)
{
  const int n = 54;
  struct route1_drv drv = { .on = true, .mu = mu, .eta = 1.0, .alpha_floor = gamma_driver_alpha_floor(esp) };
  static double J[54*54], A[54*54], VR[54*54], VI[54*54], wr[54], wi[54], work[54*64];
  int ipiv[54];
  double spe = route1_spure(esp, q);
  for (int cc = 0; cc < n; ++cc) {
    int c = route1_sym_var[cc];
    double qp[77], qm[77], fp[77], fm[77], eps = 1e-6*fmax(1.0, fabs(q[c]));
    memcpy(qp, q, sizeof qp); memcpy(qm, q, sizeof qm);
    qp[c] += eps; qm[c] -= eps;
    route1_face_flux(esp, d, qp, NULL, spe, &drv, fp);
    route1_face_flux(esp, d, qm, NULL, spe, &drv, fm);
    for (int mm = 0; mm < n; ++mm) J[mm*n + cc] = (fp[route1_sym_var[mm]] - fm[route1_sym_var[mm]])/(2*eps);
  }
  // chi-Hessian part of the K source: s[ij][ab] = d S_{K_ij} / d H_ab
  double s[9][9];
  for (int ab = 0; ab < 9; ++ab) {
    double qp[77], qm[77], sp[77], sm[77], eps = 1e-6;
    memcpy(qp, q, sizeof qp); memcpy(qm, q, sizeof qm);
    qp[68 + ab] += eps; qm[68 + ab] -= eps;
    gkyl_wv_eqn_source(esp->equation, qp, sp); gkyl_wv_eqn_source(esp->equation, qm, sm);
    for (int ij = 0; ij < 9; ++ij) s[ij][ab] = (sp[10 + ij] - sm[10 + ij])/(2*eps);
  }
  for (int ij = 0; ij < 9; ++ij) // K rows are the first 9 of the list; P_b columns are the last 3
    for (int b = 0; b < 3; ++b)
      J[ij*n + 51 + b] -= 0.5*(s[ij][3*d + b] + s[ij][3*b + d]);
  if (route1_sym_verify) { // rotation check: rows X_d against transverse columns X_l, expected -beta^l
    fprintf(stderr, "[SYM2] axis %d: A_d/A_l", d);
    for (int l = 0; l < 3; ++l) fprintf(stderr, " %+.6f(%+.6f)", J[(36 + d)*n + 36 + l], -q[52 + l]);
    fprintf(stderr, " | D_d00/D_l00");
    for (int l = 0; l < 3; ++l) fprintf(stderr, " %+.6f(%+.6f)", J[(9 + 9*d)*n + 9 + 9*l], -q[52 + l]);
    fprintf(stderr, " | P_d/P_l");
    for (int l = 0; l < 3; ++l) fprintf(stderr, " %+.6f(%+.6f)", J[(51 + d)*n + 51 + l], -q[52 + l]);
    fprintf(stderr, " | B_d0/B_l0");
    for (int l = 0; l < 3; ++l) fprintf(stderr, " %+.6f", J[(42 + 3*d)*n + 42 + 3*l]);
    fprintf(stderr, "\n");
  }
  if (route1_sym_mode == 1) {
    for (int j = 0; j < 3; ++j) for (int k = 0; k < 3; ++k) J[(42 + 3*d + j)*n + 42 + 3*k + j] -= q[52 + k];
  }
  else if (route1_sym_mode == 2) { // families: list base, stride over the derivative index, sub-components
    const int fb[4] = { 9, 36, 42, 51 }, fs[4] = { 9, 1, 3, 1 }, fn[4] = { 9, 1, 3, 1 };
    for (int f = 0; f < 4; ++f)
      for (int s2 = 0; s2 < fn[f]; ++s2) {
        for (int l = 0; l < 3; ++l) // drop the flux-form cross terms -beta^l X_l (l != d); B has none
          if (l != d && f != 2) J[(fb[f] + fs[f]*d + s2)*n + fb[f] + fs[f]*l + s2] += q[52 + l];
        for (int k = 0; k < 3; ++k) // advective term -beta^d X_k (X_d already has it, except B)
          if (k != d || f == 2) J[(fb[f] + fs[f]*k + s2)*n + fb[f] + fs[f]*k + s2] -= q[52 + d];
      }
  }
  for (int mm = 0; mm < n; ++mm) for (int cc = 0; cc < n; ++cc) A[cc*n + mm] = J[mm*n + cc];
  int N = n, lwork = 54*64, info = 0, linfo = 0;
  char nv = 'N', vv = 'V';
  dgeev_(&nv, &vv, &N, A, &N, wr, wi, NULL, &N, VR, &N, work, &lwork, &info);
  double vnorm = 0.0, inorm = 0.0;
  *wimax = 0.0;
  for (int m = 0; m < n; ++m) *wimax = fmax(*wimax, fabs(wi[m]));
  route1_sym_minabs = INFINITY; route1_sym_nzero = 0;
  for (int m = 0; m < n; ++m) {
    double a = hypot(wr[m], wi[m]);
    if (a <= 1e-7) ++route1_sym_nzero; else route1_sym_minabs = fmin(route1_sym_minabs, a);
  }
  memcpy(VI, VR, sizeof VI);
  for (int m = 0; m < n*n; ++m) vnorm += VR[m]*VR[m];
  if (info == 0) { dgetrf_(&N, &N, VI, &N, ipiv, &linfo); if (!linfo) dgetri_(&N, VI, &N, ipiv, work, &lwork, &linfo); }
  for (int m = 0; m < n*n; ++m) inorm += VI[m]*VI[m];
  *nspd = 0;
  for (int m = 0; m < n; ++m) {
    if (wr[m] <= 1e-6) continue;
    int k = 0;
    while (k < *nspd && fabs(spd[k] - wr[m]) > 1e-4) ++k;
    if (k == *nspd && *nspd < 12) spd[(*nspd)++] = wr[m];
  }
  for (int a = 1; a < *nspd; ++a) for (int b = a; b > 0 && spd[b] < spd[b-1]; --b) { double t = spd[b]; spd[b] = spd[b-1]; spd[b-1] = t; }
  return (info || linfo) ? INFINITY : sqrt(vnorm*inorm);
}

static void
route1_symbol_lines(gkyl_moment_app *app, const struct moment_species *esp, const struct gkyl_array *u)
{
  int c[3], ndim = app->ndim;
  for (int d = 0; d < ndim; ++d) c[d] = (app->global.lower[d] + app->global.upper[d] + 1)/2; // first cell with x >= 0
  gkyl_moment_app_cout(app, stderr, "[SYMBOL] columns: line r alpha chi | per variant (mu=0.75, mu=0.6, mu=0.5 a^2 chi^2):"
    " max cond over axes, max|Im|, speeds along axis 0\n");
  for (int line = 0; line < 2; ++line) {
    for (int s = 0; ; ++s) {
      int idx[3];
      for (int d = 0; d < ndim; ++d) idx[d] = c[d] + ((line == 1 || d == 0) ? s : 0);
      bool in = true;
      for (int d = 0; d < ndim; ++d) in = in && idx[d] <= app->local.upper[d];
      if (!in) break;
      const double *q = gkyl_array_cfetch(u, gkyl_range_idx(&app->local_ext, idx));
      double xc[3] = { 0, 0, 0 }, r = 0.0;
      gkyl_rect_grid_cell_center(&app->grid, idx, xc);
      for (int d = 0; d < ndim; ++d) r += xc[d]*xc[d];
      r = sqrt(r);
      char buf[1024]; int off = 0;
      for (int v = 0; v < 3; ++v) {
        double mu = v == 0 ? 0.75 : v == 1 ? 0.6 : 0.5*q[9]*q[9]*q[64]*q[64];
        double cmax = 0.0, wim = 0.0, spd0[12]; int ns0 = 0;
        for (int d = 0; d < ndim; ++d) {
          double wi, spd[12]; int ns;
          double cd = route1_symbol_cell(esp, d, q, mu, &wi, spd, &ns);
          cmax = fmax(cmax, cd); wim = fmax(wim, wi);
          if (d == 0) { ns0 = ns; memcpy(spd0, spd, sizeof spd0); }
        }
        off += snprintf(buf + off, sizeof buf - off, " | %8.1e %7.1e", cmax, wim);
        for (int k = 0; k < ns0; ++k) off += snprintf(buf + off, sizeof buf - off, " %.4f", spd0[k]);
      }
      gkyl_moment_app_cout(app, stderr, "[SYMBOL] %s r=%6.3f alpha=%.4f chi=%.4f%s\n", line ? "diag" : "x", r, q[9], q[64], buf);
    }
  }
  gkyl_moment_app_cout(app, stderr, "[SYMBOL] done\n");
}

// [SYM2] GKYL_ROUTE1_SYMBOL=2: rotation check at the first cell, then the +x line and the diagonal out to r = 6, at the
// cells and at three states linearly interpolated between neighbours, for modes 0, 1, 2 (mu = 0.75): max cond over axes,
// max |Im|, and along axis 0 the smallest nonzero |lambda| and the number of zero-speed modes.
static void
route1_symbol_scan(gkyl_moment_app *app, const struct moment_species *esp, const struct gkyl_array *u)
{
  int c[3], ndim = app->ndim;
  for (int d = 0; d < ndim; ++d) c[d] = (app->global.lower[d] + app->global.upper[d] + 1)/2;
  fprintf(stderr, "[SYM2] columns: line r alpha | mode 0 (as built), 1 (+shift advection), 2 (+advective gradient families):"
    " max cond, max|Im|, min nonzero |lambda| axis 0, zero modes axis 0\n");
  for (int line = 0; line < 2; ++line) {
    for (int s = 0; ; ++s) {
      int i0[3], i1[3];
      bool in = true;
      for (int d = 0; d < ndim; ++d) {
        i0[d] = c[d] + ((line == 1 || d == 0) ? s : 0);
        i1[d] = c[d] + ((line == 1 || d == 0) ? s + 1 : 0);
        in = in && i1[d] <= app->local.upper[d];
      }
      if (!in) break;
      const double *q0 = gkyl_array_cfetch(u, gkyl_range_idx(&app->local_ext, i0));
      const double *q1 = gkyl_array_cfetch(u, gkyl_range_idx(&app->local_ext, i1));
      double x0[3] = { 0, 0, 0 }, x1[3] = { 0, 0, 0 };
      gkyl_rect_grid_cell_center(&app->grid, i0, x0);
      gkyl_rect_grid_cell_center(&app->grid, i1, x1);
      double r0 = sqrt(x0[0]*x0[0] + x0[1]*x0[1] + x0[2]*x0[2]);
      if (r0 > 6.0) break;
      for (int fq = 0; fq < 4; ++fq) {
        double f = 0.25*fq, q[77], r = 0.0;
        for (int m = 0; m < 77; ++m) q[m] = (1 - f)*q0[m] + f*q1[m];
        for (int d = 0; d < 3; ++d) { double x = (1 - f)*x0[d] + f*x1[d]; r += x*x; }
        r = sqrt(r);
        char buf[1024]; int off = 0;
        for (int mode = 0; mode < 3; ++mode) {
          route1_sym_mode = mode;
          double cmax = 0.0, wim = 0.0, mina = 0.0; int nz = 0;
          for (int d = 0; d < ndim; ++d) {
            double wi, spd[12]; int ns;
            route1_sym_verify = (line == 0 && s == 0 && fq == 0 && mode == 0);
            double cd = route1_symbol_cell(esp, d, q, 0.75, &wi, spd, &ns);
            cmax = fmax(cmax, cd); wim = fmax(wim, wi);
            if (d == 0) { mina = route1_sym_minabs; nz = route1_sym_nzero; }
          }
          off += snprintf(buf + off, sizeof buf - off, " | %8.1e %7.1e %.4f %2d", cmax, wim, mina, nz);
        }
        route1_sym_mode = 0; route1_sym_verify = 0;
        fprintf(stderr, "[SYM2] %s r=%6.3f alpha=%.4f%s\n", line ? "diag" : "x", r, q[9], buf);
      }
    }
  }
  fprintf(stderr, "[SYM2] done\n");
}
