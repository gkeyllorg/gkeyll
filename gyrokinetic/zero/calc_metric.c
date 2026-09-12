#include <gkyl_calc_metric.h>
#include <gkyl_calc_metric_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_nodal_ops.h>
#include <gkyl_array_ops_priv.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

// test commit for slack channel 5 (after new clone and commits and pushes included on front end)

gkyl_calc_metric*
gkyl_calc_metric_new(const struct gkyl_basis *cbasis, const struct gkyl_rect_grid *grid,
  const struct gkyl_range *global, const struct gkyl_range *global_ext,
  const struct gkyl_range *local, const struct gkyl_range *local_ext,
  bool exit_at_checks, bool use_gpu)
{
  gkyl_calc_metric *up = gkyl_malloc(sizeof(gkyl_calc_metric));

  up->exit_at_checks = exit_at_checks;
  up->cbasis = cbasis;
  up->cdim = cbasis->ndim;
  up->cnum_basis = cbasis->num_basis;
  up->poly_order = cbasis->poly_order;
  up->grid = grid;
  up->use_gpu = use_gpu;
  up->num_cells = up->grid->cells;
  up->n2m = gkyl_nodal_ops_new(up->cbasis, up->grid, up->use_gpu);

  up->global = *global;
  up->global_ext = *global_ext;

  up->local = *local;
  up->local_ext = *local_ext;
  up->rz_expected_jacobian_sign = 0;
  up->rz_jacobian_valid = true;
  up->rz_jacobian_guard_nonfatal = false;

  return up;
}

void
gkyl_calc_metric_set_signed_jacobian_guard_nonfatal(gkyl_calc_metric *up,
  bool nonfatal)
{
  up->rz_jacobian_guard_nonfatal = nonfatal;
}

bool
gkyl_calc_metric_signed_jacobian_valid(const gkyl_calc_metric *up)
{
  return up->rz_jacobian_valid;
}

int
gkyl_calc_metric_signed_jacobian_sign(const gkyl_calc_metric *up)
{
  return up->rz_expected_jacobian_sign;
}

static inline double calc_metric(double dxdz[3][3], int i, int j) 
{
  double sum = 0;
  for (int k = 0; k < 3; ++k)
    sum += dxdz[k][i - 1] * dxdz[k][j - 1];
  return sum;
} 

// Calculates e^1 = e_2 x e_3 /J
static inline void
calc_dual(double J, const double e_2[3], const double e_3[3], double e1[3])
{
  e1[0] = (e_2[1]*e_3[2] - e_2[2]*e_3[1] )/J;
  e1[1] = -(e_2[0]*e_3[2] - e_2[2]*e_3[0] )/J;
  e1[2] = (e_2[0]*e_3[1] - e_2[1]*e_3[0] )/J;
}

static inline void
matTvec(double M[3][3], double v[3], double result[3]) {
  for (int i = 0; i < 3; i++) {
    result[i] = 0.0;
    for (int j = 0; j < 3; j++) {
      result[i] += M[j][i] * v[j];
    }
  }
}

static inline double dot(const double a[3], const double b[3]) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static inline double
signed_jacobian_rz(double R, double dRdpsi, double dZdpsi, double dRdtheta, double dZdtheta)
{
  return R*(dRdpsi*dZdtheta - dRdtheta*dZdpsi);
}

struct signed_jacobian_guard_state {
  int expected_sign;
  bool valid;
  bool nonfatal;
};

static bool
signed_jacobian_guard_enabled(void)
{
  static int enabled = -1;
  if (enabled < 0) {
    const char *value = getenv("GKYL_MAP_JACOBIAN_SIGN_GUARD");
    enabled = value && value[0] != '\0' && value[0] != '0';
    if (enabled) {
      fprintf(stderr, "GKYL_SIGNED_JACOBIAN_GUARD enabled version=1\n");
      fflush(stderr);
    }
  }
  return enabled != 0;
}

static void
signed_jacobian_guard_check(struct signed_jacobian_guard_state *state,
  double J, bool enforce_sign, const char *location, int dir,
  const int cidx[3])
{
  bool enabled = signed_jacobian_guard_enabled();
  if (!isfinite(J) || J == 0.0) {
    state->valid = false;
    if (!enabled || state->nonfatal)
      return;
    fprintf(stderr,
      "GKYL_SIGNED_JACOBIAN_GUARD failure version=1 location=%s dir=%d index=%d,%d,%d J=%.17g reason=%s\n",
      location, dir, cidx[0], cidx[1], cidx[2], J,
      isfinite(J) ? "zero" : "nonfinite");
    abort();
  }
  if (!enforce_sign)
    return;
  int sign = J < 0.0 ? -1 : 1;
  if (state->expected_sign == 0)
    state->expected_sign = sign;
  else if (sign != state->expected_sign) {
    state->valid = false;
    if (!enabled || state->nonfatal)
      return;
    fprintf(stderr,
      "GKYL_SIGNED_JACOBIAN_GUARD failure version=1 location=%s dir=%d index=%d,%d,%d J=%.17g expected_sign=%d reason=sign_reversal\n",
      location, dir, cidx[0], cidx[1], cidx[2], J,
      state->expected_sign);
    abort();
  }
}

static inline void cross(const double a[3], const double b[3], double c[3]) {
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

static inline void
check_orthonormality(const double tan[9], const double dual[9], bool exit_at_check) {
  // Check that the coordinate system has tangent/dual vectors
  // satisfying orthonormality.

  // Calculate e_i \dot e^j.
  double prod[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      prod[i][j] = tan[3*i +0]*dual[3*j +0] + tan[3*i +1]*dual[3*j +1] + tan[3*i +2]*dual[3*j +2];
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if ( i==j && prod[i][j] < 0 ) {
        fprintf(stderr, "calc_metric.c: Orthonormality violated : e_%d . e^%d = %.6e\n", i+1, j+1, prod[i][j]);
        fprintf(stderr, "calc_metric.c: diagonal products = [%.6e, %.6e, %.6e]\n",
          prod[0][0], prod[1][1], prod[2][2]);
        assert(!exit_at_check);
      }
    }
  }

  // The loop above tests only the SIGN of the diagonal, so a metric wrong by a
  // scale factor, or carrying non-zero off-diagonal products, passes it
  // silently.  Measure the full identity |e_i . e^j - delta_ij| as well.
  //
  // This REPORTS by default and is fatal only under
  // GKYL_METRIC_STRICT_ORTHONORMALITY=1, following the same opt-in pattern as
  // the seam-participation guard.  The coordinate system is genuinely singular
  // at X-point and divertor-plate corners, so a check that aborted by default
  // would reject geometry the plan explicitly expects to build; what matters is
  // WHERE the residual is large, not merely that it is non-zero somewhere.
  //
  // The tolerance is deliberately loose by default and tunable, because the
  // representative residual has not yet been measured across devices; pick a
  // threshold from that measurement rather than asserting one up front.
  {
    static double tol = -1.0;
    static int strict = -1;
    static long reported = 0;
    const long report_cap = 20;
    if (tol < 0.0) {
      const char *s = getenv("GKYL_METRIC_ORTHONORMALITY_TOL");
      tol = (s && s[0]) ? atof(s) : 1.0e-6;
      if (!(tol > 0.0)) tol = 1.0e-6;
    }
    if (strict < 0) {
      const char *s = getenv("GKYL_METRIC_STRICT_ORTHONORMALITY");
      strict = (s && s[0] && s[0] != '0') ? 1 : 0;
    }
    double worst = 0.0; int wi = 0, wj = 0;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        double d = fabs(prod[i][j] - (i == j ? 1.0 : 0.0));
        if (d > worst) { worst = d; wi = i; wj = j; }
      }
    }
    if (worst > tol) {
      if (reported < report_cap) {
        fprintf(stderr,
          "TOK_METRIC_ORTHONORMALITY residual=%.9e at e_%d.e^%d tol=%.3e "
          "diag=[%.9e,%.9e,%.9e]\n",
          worst, wi+1, wj+1, tol, prod[0][0], prod[1][1], prod[2][2]);
        if (++reported == report_cap)
          fprintf(stderr, "TOK_METRIC_ORTHONORMALITY further reports suppressed "
                          "after %ld\n", report_cap);
      }
      if (strict) {
        fprintf(stderr, "calc_metric.c: orthonormality residual %.9e exceeds "
                        "tolerance %.3e under GKYL_METRIC_STRICT_ORTHONORMALITY\n",
                        worst, tol);
        assert(!exit_at_check);
      }
    }
  }
}

static inline void
check_right_handed(const double tan[9], const double dual[9], bool exit_at_check) {
  // Check that the coordinate system is right handed.
  const double *e1 = &tan[0];
  const double *e2 = &tan[3];
  const double *e3 = &tan[6];

  double cross_e2e3[3];
  cross(e2, e3, cross_e2e3);

  double J = dot(e1, cross_e2e3);

  if (J < 0.0) {
    fprintf(stderr, "calc_metric.c: Left-handed coordinate system, J = e_1 . (e_2 x e_3) = %.6e < 0.\n", J);
    assert(!exit_at_check);
  }
  else if (J == 0.0) {
    fprintf(stderr, "calc_metric.c: Degenerate coordinate system, J = %.6e\n", J);
    assert(!exit_at_check);
  }
}

static inline void
check_parallel(double *v1, double *v2, bool exit_at_check) {
  // Check v1 and v2 are parallel by checking that:
  //   |v1 x v2 | < eps 
  //   |v1 . v2 - 1| < eps
  const double eps = 1e-3;

  double cx = v1[1]*v2[2] - v1[2]*v2[1];
  double cy = v1[2]*v2[0] - v1[0]*v2[2];
  double cz = v1[0]*v2[1] - v1[1]*v2[0];

  double c_mag = sqrt(cx*cx + cy*cy + cz*cz);

  double dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
  if (fabs(c_mag) < eps && fabs(dot-1.0) < eps)
    return;
  else {
    fprintf(stderr, "calc_metric.c: inconsistent B & mapc2p (hat{b} not parallel to e_3; |b . e_3|=%.6e, |b x e_3|=%.6e).\n",dot,c_mag);
    assert(!exit_at_check);
  }
}

static inline void
check_axisymmetric(struct gkyl_array* arr, struct gkyl_range *range, bool exit_at_check) {
  const double rel_tol = 1e-6;
  const double abs_tol = 1e-11;
  double reldiff;
  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  int cidx[3];
  int cidx_prev[3];

  for (int ip=range->lower[PSI_IDX]; ip<=range->upper[PSI_IDX]; ++ip) {
    for (int it=range->lower[TH_IDX]; it<=range->upper[TH_IDX]; ++it) {

      double g_ij_avg[6] = {0.0};
      int num_al = range->upper[AL_IDX] - range->lower[AL_IDX] + 1;
      for (int ia=range->lower[AL_IDX]; ia<=range->upper[AL_IDX]; ++ia) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        const double *g_ij_n = gkyl_array_cfetch(arr, gkyl_range_idx(range, cidx));
        for (int k = 0; k < 6; ++k) {
            g_ij_avg[k] += g_ij_n[k] / num_al;
        }
      }

      for (int ia=range->lower[AL_IDX] + 1; ia<=range->upper[AL_IDX]; ++ia) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        const double *g_ij_n = gkyl_array_cfetch(arr, gkyl_range_idx(range, cidx));

        cidx_prev[PSI_IDX] = ip;
        cidx_prev[AL_IDX] = ia - 1;
        cidx_prev[TH_IDX] = it;
        const double *g_ij_prev = gkyl_array_cfetch(arr, gkyl_range_idx(range, cidx_prev));

        for (int k = 0; k < 6; ++k) {
          reldiff = fabs(g_ij_n[k] - g_ij_avg[k])/(rel_tol*fabs(g_ij_avg[k]) + abs_tol);
          if (reldiff < 1) {
            return;
          }
          else {
            fprintf(stderr, "calc_metric.c: Axisymmetry violated at ip=%d, it=%d, ia=%d. g_ij component %d variation %.6e exceeds tolerance\n", ip, it, ia, k, reldiff);
            assert(!exit_at_check);
          }
        }

      }
    }
  }
}

static int metric_diag_block_id = -1;

struct metric_diag_stats {
  long nodes;
  long bad_jacobgeo;
  long bad_bmag;
  long bad_gdiag;
  double min_jacobgeo, max_jacobgeo, min_abs_jacobgeo;
  double min_bmag, max_bmag;
  double min_g11, min_g22, min_g33;
  double max_g11, max_g22, max_g33;
  int min_abs_jacobgeo_idx[3];
  int max_g11_idx[3];
  int max_g22_idx[3];
  int max_g33_idx[3];
};

static void
metric_diag_stats_init(struct metric_diag_stats *s)
{
  s->nodes = 0;
  s->bad_jacobgeo = 0;
  s->bad_bmag = 0;
  s->bad_gdiag = 0;
  s->min_jacobgeo = DBL_MAX;
  s->max_jacobgeo = -DBL_MAX;
  s->min_abs_jacobgeo = DBL_MAX;
  s->min_bmag = DBL_MAX;
  s->max_bmag = -DBL_MAX;
  s->min_g11 = DBL_MAX;
  s->min_g22 = DBL_MAX;
  s->min_g33 = DBL_MAX;
  s->max_g11 = -DBL_MAX;
  s->max_g22 = -DBL_MAX;
  s->max_g33 = -DBL_MAX;
  for (int d=0; d<3; ++d) {
    s->min_abs_jacobgeo_idx[d] = 0;
    s->max_g11_idx[d] = 0;
    s->max_g22_idx[d] = 0;
    s->max_g33_idx[d] = 0;
  }
}

static void
metric_diag_update_minmax(double v, double *minv, double *maxv)
{
  if (v < *minv) *minv = v;
  if (v > *maxv) *maxv = v;
}

static void
metric_diag_scan(struct metric_diag_stats *s, const struct gkyl_range *range,
  const struct gkyl_array *jacobgeo, const struct gkyl_array *bmag,
  const struct gkyl_array *g_ij)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(range, iter.idx);
    const double *j_n = gkyl_array_cfetch(jacobgeo, linidx);
    const double *b_n = gkyl_array_cfetch(bmag, linidx);
    const double *g_n = gkyl_array_cfetch(g_ij, linidx);

    s->nodes++;
    if (!isfinite(j_n[0]) || j_n[0] <= 0.0)
      s->bad_jacobgeo++;
    else {
      metric_diag_update_minmax(j_n[0], &s->min_jacobgeo, &s->max_jacobgeo);
      double abs_j = fabs(j_n[0]);
      if (abs_j < s->min_abs_jacobgeo) {
        s->min_abs_jacobgeo = abs_j;
        for (int d=0; d<3; ++d)
          s->min_abs_jacobgeo_idx[d] = iter.idx[d];
      }
    }

    if (!isfinite(b_n[0]) || b_n[0] <= 0.0)
      s->bad_bmag++;
    else
      metric_diag_update_minmax(b_n[0], &s->min_bmag, &s->max_bmag);

    if (!isfinite(g_n[0]) || !isfinite(g_n[3]) || !isfinite(g_n[5]) ||
        g_n[0] <= 0.0 || g_n[3] <= 0.0 || g_n[5] <= 0.0) {
      s->bad_gdiag++;
    }
    else {
      metric_diag_update_minmax(g_n[0], &s->min_g11, &s->max_g11);
      metric_diag_update_minmax(g_n[3], &s->min_g22, &s->max_g22);
      metric_diag_update_minmax(g_n[5], &s->min_g33, &s->max_g33);
      if (g_n[0] >= s->max_g11) {
        for (int d=0; d<3; ++d)
          s->max_g11_idx[d] = iter.idx[d];
      }
      if (g_n[3] >= s->max_g22) {
        for (int d=0; d<3; ++d)
          s->max_g22_idx[d] = iter.idx[d];
      }
      if (g_n[5] >= s->max_g33) {
        for (int d=0; d<3; ++d)
          s->max_g33_idx[d] = iter.idx[d];
      }
    }
  }
}

static void
metric_diag_write(const char *location, int dir, const struct metric_diag_stats *s)
{
  const char *csv_path = getenv("GKYL_METRIC_DIAG_CSV");
  if (!csv_path || !csv_path[0])
    return;
  const char *case_name = getenv("GKYL_METRIC_DIAG_CASE");
  if (!case_name || !case_name[0])
    case_name = "unknown";

  fprintf(stderr,
    "METRIC_DIAG case=%s block=%d location=%s dir=%d nodes=%ld bad_jacobgeo=%ld "
    "bad_bmag=%ld bad_gdiag=%ld min_jacobgeo=%.16e max_jacobgeo=%.16e "
    "min_abs_jacobgeo=%.16e min_abs_jacobgeo_idx=(%d,%d,%d) "
    "min_bmag=%.16e max_bmag=%.16e min_g11=%.16e min_g22=%.16e min_g33=%.16e "
    "max_g11=%.16e max_g11_idx=(%d,%d,%d) max_g22=%.16e max_g22_idx=(%d,%d,%d) "
    "max_g33=%.16e max_g33_idx=(%d,%d,%d)\n",
    case_name, metric_diag_block_id, location, dir, s->nodes, s->bad_jacobgeo, s->bad_bmag,
    s->bad_gdiag, s->min_jacobgeo, s->max_jacobgeo, s->min_abs_jacobgeo,
    s->min_abs_jacobgeo_idx[0], s->min_abs_jacobgeo_idx[1], s->min_abs_jacobgeo_idx[2],
    s->min_bmag, s->max_bmag, s->min_g11, s->min_g22, s->min_g33,
    s->max_g11, s->max_g11_idx[0], s->max_g11_idx[1], s->max_g11_idx[2],
    s->max_g22, s->max_g22_idx[0], s->max_g22_idx[1], s->max_g22_idx[2],
    s->max_g33, s->max_g33_idx[0], s->max_g33_idx[1], s->max_g33_idx[2]);

  FILE *csv = fopen(csv_path, "a");
  if (!csv) {
    fprintf(stderr, "METRIC_DIAG failed_to_open_csv=%s\n", csv_path);
    return;
  }

  fprintf(csv,
    "%s,%d,%s,%d,%ld,%ld,%ld,%ld,%.16e,%.16e,%.16e,%d,%d,%d,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%d,%d,%d,%.16e,%d,%d,%d,%.16e,%d,%d,%d\n",
    case_name, metric_diag_block_id, location, dir, s->nodes, s->bad_jacobgeo,
    s->bad_bmag, s->bad_gdiag, s->min_jacobgeo, s->max_jacobgeo,
    s->min_abs_jacobgeo, s->min_abs_jacobgeo_idx[0], s->min_abs_jacobgeo_idx[1],
    s->min_abs_jacobgeo_idx[2], s->min_bmag, s->max_bmag, s->min_g11,
    s->min_g22, s->min_g33, s->max_g11, s->max_g11_idx[0], s->max_g11_idx[1],
    s->max_g11_idx[2], s->max_g22, s->max_g22_idx[0], s->max_g22_idx[1],
    s->max_g22_idx[2], s->max_g33, s->max_g33_idx[0], s->max_g33_idx[1],
    s->max_g33_idx[2]);
  fclose(csv);
}

static void
metric_diag_write_interior(const struct gk_geometry *gk_geom)
{
  metric_diag_block_id++;
  struct metric_diag_stats s;
  metric_diag_stats_init(&s);
  metric_diag_scan(&s, &gk_geom->nrange_int, gk_geom->geo_int.jacobgeo_nodal,
    gk_geom->geo_int.bmag_nodal, gk_geom->geo_int.g_ij_nodal);
  metric_diag_write("interior", -1, &s);
}

static void
metric_diag_write_surface(int dir, const struct gk_geometry *gk_geom)
{
  struct metric_diag_stats s;
  metric_diag_stats_init(&s);
  metric_diag_scan(&s, &gk_geom->nrange_surf[dir],
    gk_geom->geo_surf[dir].jacobgeo_nodal, gk_geom->geo_surf[dir].bmag_nodal,
    gk_geom->geo_surf[dir].g_ij_nodal);
  metric_diag_write("surface", dir, &s);
}

void gkyl_calc_metric_advance_rz( gkyl_calc_metric *up, struct gkyl_range *nrange,
  struct gkyl_array *mc2p_nodal_fd, struct gkyl_array *ddtheta_nodal,
  struct gkyl_array *bmag_nodal, double *dzc, struct gkyl_array *gFld,
  struct gkyl_array *tanvecFld, struct gkyl_array *dualFld,
  struct gkyl_array *dualmagFld, struct gkyl_array *normFld,
  struct gkyl_array *jFld, struct gkyl_array* bcartFld, const struct gkyl_range *update_range)
{
  struct gkyl_array* gFld_nodal = gkyl_array_new(GKYL_DOUBLE, 6, nrange->volume);
  struct gkyl_array* jFld_nodal = gkyl_array_new(GKYL_DOUBLE, 1, nrange->volume);
  struct gkyl_array* bcartFld_nodal = gkyl_array_new(GKYL_DOUBLE, 3, nrange->volume);
  struct gkyl_array* tanvecFld_nodal = gkyl_array_new(GKYL_DOUBLE, 9, nrange->volume);
  struct gkyl_array* dualFld_nodal = gkyl_array_new(GKYL_DOUBLE, 9, nrange->volume);
  struct gkyl_array* dualmagFld_nodal = gkyl_array_new(GKYL_DOUBLE, 3, nrange->volume);
  struct gkyl_array* normFld_nodal = gkyl_array_new(GKYL_DOUBLE, 9, nrange->volume);
  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { R_IDX, Z_IDX, PHI_IDX }; // arrangement of cartesian coordinates
  int cidx[3];
  struct signed_jacobian_guard_state jac_guard = { 0 };
  for(int ia=nrange->lower[AL_IDX]; ia<=nrange->upper[AL_IDX]; ++ia){
    for (int ip=nrange->lower[PSI_IDX]; ip<=nrange->upper[PSI_IDX]; ++ip) {
      for (int it=nrange->lower[TH_IDX]; it<=nrange->upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        const double *mc2p_n = gkyl_array_cfetch(mc2p_nodal_fd, gkyl_range_idx(nrange, cidx));
        double dxdz[3][3];

        if((ip == nrange->lower[PSI_IDX]) && (up->local.lower[PSI_IDX]== up->global.lower[PSI_IDX]) ) {
          dxdz[0][0] = (-3*mc2p_n[R_IDX] + 4*mc2p_n[6+R_IDX] - mc2p_n[12+R_IDX] )/dzc[0]/2;
          dxdz[1][0] = (-3*mc2p_n[Z_IDX] + 4*mc2p_n[6+Z_IDX] - mc2p_n[12+Z_IDX] )/dzc[0]/2;
          dxdz[2][0] = (-3*mc2p_n[PHI_IDX] + 4*mc2p_n[6+PHI_IDX] - mc2p_n[12+PHI_IDX] );
        }
        else if((ip == nrange->upper[PSI_IDX]) && (up->local.upper[PSI_IDX]== up->global.upper[PSI_IDX])) {
          dxdz[0][0] = (3*mc2p_n[R_IDX] - 4*mc2p_n[3+R_IDX] + mc2p_n[9+R_IDX] )/dzc[0]/2;
          dxdz[1][0] = (3*mc2p_n[Z_IDX] - 4*mc2p_n[3+Z_IDX] + mc2p_n[9+Z_IDX] )/dzc[0]/2;
          dxdz[2][0] = (3*mc2p_n[PHI_IDX] - 4*mc2p_n[3+PHI_IDX] + mc2p_n[9+PHI_IDX] );
        }
        else {
          dxdz[0][0] = -(mc2p_n[3 +R_IDX] -   mc2p_n[6+R_IDX])/2/dzc[0];
          dxdz[1][0] = -(mc2p_n[3 +Z_IDX] -   mc2p_n[6+Z_IDX])/2/dzc[0];
          dxdz[2][0] = -(mc2p_n[3 +PHI_IDX] -   mc2p_n[6+PHI_IDX]);
        }
        // Take into account wrapping of cyclic coordinate phi
        if (dxdz[2][0] < -M_PI) {
          dxdz[2][0] += 2*M_PI;
        } else if (dxdz[2][0] > M_PI) {
          dxdz[2][0] -= 2*M_PI;
        }
        dxdz[2][0] = dxdz[2][0]/2.0/dzc[0];

        if((ia == nrange->lower[AL_IDX]) && (up->local.lower[AL_IDX]== up->global.lower[AL_IDX]) ) {
          dxdz[0][1] = (-3*mc2p_n[R_IDX] +  4*mc2p_n[18+R_IDX] -  mc2p_n[24+R_IDX])/dzc[1]/2;
          dxdz[1][1] = (-3*mc2p_n[Z_IDX] +  4*mc2p_n[18+Z_IDX] -  mc2p_n[24+Z_IDX])/dzc[1]/2;
          dxdz[2][1] = (-3*mc2p_n[PHI_IDX] +  4*mc2p_n[18+PHI_IDX] -  mc2p_n[24+PHI_IDX])/dzc[1]/2;
        }
        else if((ia == nrange->upper[AL_IDX])  && (up->local.upper[AL_IDX]== up->global.upper[AL_IDX])){
          dxdz[0][1] = (3*mc2p_n[R_IDX] -  4*mc2p_n[15+R_IDX] +  mc2p_n[21+R_IDX] )/dzc[1]/2;
          dxdz[1][1] = (3*mc2p_n[Z_IDX] -  4*mc2p_n[15+Z_IDX] +  mc2p_n[21+Z_IDX] )/dzc[1]/2;
          dxdz[2][1] = (3*mc2p_n[PHI_IDX] -  4*mc2p_n[15+PHI_IDX] +  mc2p_n[21+PHI_IDX] )/dzc[1]/2;
        }
        else {
          dxdz[0][1] = -(mc2p_n[15 +R_IDX] -  mc2p_n[18 +R_IDX])/2/dzc[1];
          dxdz[1][1] = -(mc2p_n[15 +Z_IDX] -  mc2p_n[18 +Z_IDX])/2/dzc[1];
          dxdz[2][1] = -(mc2p_n[15 +PHI_IDX] -  mc2p_n[18 +PHI_IDX])/2/dzc[1];
        }

        if((it == nrange->lower[TH_IDX]) && (up->local.lower[TH_IDX]== up->global.lower[TH_IDX])){
          dxdz[0][2] = (-3*mc2p_n[R_IDX] + 4*mc2p_n[30+R_IDX] - mc2p_n[36+R_IDX])/dzc[2]/2;
          dxdz[1][2] = (-3*mc2p_n[Z_IDX] + 4*mc2p_n[30+Z_IDX] - mc2p_n[36+Z_IDX])/dzc[2]/2;
          dxdz[2][2] = (-3*mc2p_n[PHI_IDX] + 4*mc2p_n[30+PHI_IDX] - mc2p_n[36+PHI_IDX])/dzc[2]/2;
        }
        else if((it == nrange->upper[TH_IDX]) && (up->local.upper[TH_IDX]== up->global.upper[TH_IDX])){
          dxdz[0][2] = (3*mc2p_n[R_IDX] - 4*mc2p_n[27+R_IDX] + mc2p_n[33+R_IDX] )/dzc[2]/2;
          dxdz[1][2] = (3*mc2p_n[Z_IDX] - 4*mc2p_n[27+Z_IDX] + mc2p_n[33+Z_IDX] )/dzc[2]/2;
          dxdz[2][2] = (3*mc2p_n[PHI_IDX] - 4*mc2p_n[27+PHI_IDX] + mc2p_n[33+PHI_IDX] )/dzc[2]/2;
        }
        else {
          dxdz[0][2] = -(mc2p_n[27 +R_IDX] - mc2p_n[30 +R_IDX])/2/dzc[2];
          dxdz[1][2] = -(mc2p_n[27 +Z_IDX] - mc2p_n[30 +Z_IDX])/2/dzc[2];
          dxdz[2][2] = -(mc2p_n[27 +PHI_IDX] - mc2p_n[30 +PHI_IDX])/2/dzc[2];
        }

        // Take into account wrapping of cyclic coordinate phi
        if (dxdz[2][0] < -M_PI) {
          dxdz[2][0] += 2*M_PI;
        } else if (dxdz[2][0] > M_PI) {
          dxdz[2][0] -= 2*M_PI;
        }
        dxdz[2][0] = dxdz[2][0]/2.0/dzc[0];

        // Use exact expressions for dR/dtheta and dZ/dtheta
        double *ddtheta_n = gkyl_array_fetch(ddtheta_nodal, gkyl_range_idx(nrange, cidx));
        dxdz[0][2] = ddtheta_n[0]; // dR/dtheta
        dxdz[1][2] = ddtheta_n[1]; // dZ/dtheta

        // use exact expressions for d/dalpha
        dxdz[0][1] = 0.0; // dR/dalpha
        dxdz[1][1] = 0.0; // dZ/dalpha
        dxdz[2][1] = -1.0; // dphi/dalpha

        // dxdz is in cylindrical coords, calculate J as
        // J = R(dR/dpsi*dZ/dtheta - dR/dtheta*dZ/dpsi)
        double *jFld_n= gkyl_array_fetch(jFld_nodal, gkyl_range_idx(nrange, cidx));
        double R = mc2p_n[R_IDX];
        double signed_J = signed_jacobian_rz(R, dxdz[0][0], dxdz[1][0],
          dxdz[0][2], dxdz[1][2]);
        signed_jacobian_guard_check(&jac_guard, signed_J, true,
          "rz", -1, cidx);
        jFld_n[0] = fabs(signed_J);

        // Calculate dphi/dtheta based on the divergence free condition
        // on B: 1 = J*B/sqrt(g_33)
        double *bmag_n = gkyl_array_fetch(bmag_nodal, gkyl_range_idx(nrange, cidx));
        double dphidtheta = (jFld_n[0]*jFld_n[0]*bmag_n[0]*bmag_n[0] - dxdz[0][2]*dxdz[0][2] - dxdz[1][2]*dxdz[1][2])/R/R;
        dphidtheta = sqrt(dphidtheta);
        // Recover sign from exact dphidtheta = F(psi)/R/\grad(psi).
        if (ddtheta_n[2] < 0) {
          dphidtheta = -dphidtheta;
        }

        double *gFld_n= gkyl_array_fetch(gFld_nodal, gkyl_range_idx(nrange, cidx));
        gFld_n[0] = dxdz[0][0]*dxdz[0][0] + R*R*dxdz[2][0]*dxdz[2][0] + dxdz[1][0]*dxdz[1][0]; 
        gFld_n[1] = R*R*dxdz[2][0]; 
        gFld_n[2] = dxdz[0][0]*dxdz[0][2] + R*R*dxdz[2][0]*dphidtheta + dxdz[1][0]*dxdz[1][2];
        gFld_n[3] = R*R; 
        gFld_n[4] = R*R*dphidtheta;
        gFld_n[5] = dxdz[0][2]*dxdz[0][2] + R*R*dphidtheta*dphidtheta + dxdz[1][2]*dxdz[1][2]; 

        // Calculate cartesian components of bhat
        double *bcartFld_n= gkyl_array_fetch(bcartFld_nodal, gkyl_range_idx(nrange, cidx));
        double phi = mc2p_n[PHI_IDX];
        double b3 = 1/sqrt(gFld_n[5]);
        bcartFld_n[0] = b3*(dxdz[0][2]*cos(phi) - R*sin(phi)*dphidtheta);
        bcartFld_n[1] = b3*(dxdz[0][2]*sin(phi) + R*cos(phi)*dphidtheta);
        bcartFld_n[2] = b3*(dxdz[1][2]);

        // Set cartesian components of tangents and duals
        double Z = mc2p_n[Z_IDX];
        double J = jFld_n[0];
        double *tanvecFld_n= gkyl_array_fetch(tanvecFld_nodal, gkyl_range_idx(nrange, cidx));
        tanvecFld_n[0] = dxdz[0][0]*cos(phi) - R*sin(phi)*dxdz[2][0]; 
        tanvecFld_n[1] = dxdz[0][0]*sin(phi)  + R*cos(phi)*dxdz[2][0]; 
        tanvecFld_n[2] = dxdz[1][0];

        tanvecFld_n[3] = +R*sin(phi); 
        tanvecFld_n[4] = -R*cos(phi); 
        tanvecFld_n[5] = 0.0; 

        tanvecFld_n[6] = dxdz[0][2]*cos(phi) - R*sin(phi)*dphidtheta; 
        tanvecFld_n[7] = dxdz[0][2]*sin(phi)  + R*cos(phi)*dphidtheta; 
        tanvecFld_n[8] = dxdz[1][2];

        double *dualFld_n= gkyl_array_fetch(dualFld_nodal, gkyl_range_idx(nrange, cidx));
        dualFld_n[0] = -R/J*cos(phi)*dxdz[1][2];
        dualFld_n[1] = -R/J*sin(phi)*dxdz[1][2];
        dualFld_n[2] = +R/J*dxdz[0][2];

        dualFld_n[3] =  1/J * ( dxdz[1][0]*dxdz[0][2]*sin(phi) + dxdz[1][0]*R*cos(phi)*dphidtheta
                               -dxdz[1][2]*dxdz[0][0]*sin(phi) - dxdz[1][2]*R*cos(phi)*dxdz[2][0] );
        dualFld_n[4] = -1/J * ( dxdz[1][0]*dxdz[0][2]*cos(phi) + dxdz[1][0]*R*sin(phi)*dphidtheta
                               -dxdz[1][2]*dxdz[0][0]*cos(phi) - dxdz[1][2]*R*sin(phi)*dxdz[2][0] );
        dualFld_n[5] =  R/J * ( dxdz[0][2]*dxdz[2][0] - dxdz[0][0]*dphidtheta);

        dualFld_n[6] = +R/J*cos(phi)*dxdz[1][0];
        dualFld_n[7] = +R/J*sin(phi)*dxdz[1][0];
        dualFld_n[8] = -R/J*dxdz[0][0];

        double norm1 = sqrt(dualFld_n[0]*dualFld_n[0] + dualFld_n[1]*dualFld_n[1] + dualFld_n[2]*dualFld_n[2]);
        double norm2 = sqrt(dualFld_n[3]*dualFld_n[3] + dualFld_n[4]*dualFld_n[4] + dualFld_n[5]*dualFld_n[5]);
        double norm3 = sqrt(dualFld_n[6]*dualFld_n[6] + dualFld_n[7]*dualFld_n[7] + dualFld_n[8]*dualFld_n[8]);

        double *dualmagFld_n = gkyl_array_fetch(dualmagFld_nodal, gkyl_range_idx(nrange, cidx));
        dualmagFld_n[0] = norm1;
        dualmagFld_n[1] = norm2;
        dualmagFld_n[2] = norm3;
        
        // Set normal vectors
        double *normFld_n = gkyl_array_fetch(normFld_nodal, gkyl_range_idx(nrange, cidx));
        normFld_n[0] = dualFld_n[0]/norm1;
        normFld_n[1] = dualFld_n[1]/norm1;
        normFld_n[2] = dualFld_n[2]/norm1;

        normFld_n[3] = dualFld_n[3]/norm2;
        normFld_n[4] = dualFld_n[4]/norm2;
        normFld_n[5] = dualFld_n[5]/norm2;

        normFld_n[6] = dualFld_n[6]/norm3;
        normFld_n[7] = dualFld_n[7]/norm3;
        normFld_n[8] = dualFld_n[8]/norm3;
      }
    }
  }
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 6, gFld_nodal, gFld, false);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 1, jFld_nodal, jFld, false);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 3, bcartFld_nodal, bcartFld, false);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 9, tanvecFld_nodal, tanvecFld, false);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 9, dualFld_nodal, dualFld, false);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 3, dualmagFld_nodal, dualmagFld, false);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 9, normFld_nodal, normFld, false);
  gkyl_array_release(gFld_nodal);
  gkyl_array_release(jFld_nodal);
  gkyl_array_release(bcartFld_nodal);
  gkyl_array_release(tanvecFld_nodal);
  gkyl_array_release(dualFld_nodal);
  gkyl_array_release(dualmagFld_nodal);
  gkyl_array_release(normFld_nodal);
}

void 
gkyl_calc_metric_advance_rz_interior(gkyl_calc_metric *up, struct gk_geometry *gk_geom)
{
  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { R_IDX, Z_IDX, PHI_IDX }; // arrangement of cartesian coordinates
  int cidx[3];
  struct signed_jacobian_guard_state jac_guard = {
    .valid = true,
    .nonfatal = up->rz_jacobian_guard_nonfatal,
  };
  for(int ia=gk_geom->nrange_int.lower[AL_IDX]; ia<=gk_geom->nrange_int.upper[AL_IDX]; ++ia){
    for (int ip=gk_geom->nrange_int.lower[PSI_IDX]; ip<=gk_geom->nrange_int.upper[PSI_IDX]; ++ip) {
      for (int it=gk_geom->nrange_int.lower[TH_IDX]; it<=gk_geom->nrange_int.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        const double *mc2p_n = gkyl_array_cfetch(gk_geom->geo_int.mc2p_nodal_fd, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        double dxdz[3][3];

        dxdz[0][0] = -(mc2p_n[3 +R_IDX] -   mc2p_n[6+R_IDX])/2/gk_geom->dzc[0];
        dxdz[1][0] = -(mc2p_n[3 +Z_IDX] -   mc2p_n[6+Z_IDX])/2/gk_geom->dzc[0];
        dxdz[2][0] = -(mc2p_n[3 +PHI_IDX] -   mc2p_n[6+PHI_IDX]);

        // Take into account wrapping of cyclic coordinate phi
        if (dxdz[2][0] < -M_PI) {
          dxdz[2][0] += 2*M_PI;
        } else if (dxdz[2][0] > M_PI) {
          dxdz[2][0] -= 2*M_PI;
        }
        dxdz[2][0] = dxdz[2][0]/2.0/gk_geom->dzc[0];

        // Use exact expressions for dR/dtheta and dZ/dtheta
        double *ddtheta_n = gkyl_array_fetch(gk_geom->geo_int.ddtheta_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        dxdz[0][2] = ddtheta_n[0];
        dxdz[1][2] = ddtheta_n[1];

        // Get position map deriv for psi
        double *ddpsi_n = gkyl_array_fetch(gk_geom->geo_int.ddpsi_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));

        // dxdz is in cylindrical coords, calculate J as
        // J = R(dR/dpsi*dZ/dtheta - dR/dtheta*dZ/dpsi)
        double *jFld_n= gkyl_array_fetch(gk_geom->geo_int.jacobgeo_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        double R = mc2p_n[R_IDX];
        double J = signed_jacobian_rz(R, dxdz[0][0], dxdz[1][0], dxdz[0][2], dxdz[1][2]);
        signed_jacobian_guard_check(&jac_guard, J, true,
          "interior", -1, cidx);
        jFld_n[0] = fabs(J);
        double alpha_sign = J < 0.0 ? -1.0 : 1.0;

        // Calculate dphi/dtheta based on the divergence free condition
        // on B: 1 = J*B/sqrt(g_33)
        double *bmag_n = gkyl_array_fetch(gk_geom->geo_int.bmag_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        double dphidtheta = (jFld_n[0]*jFld_n[0]*bmag_n[0]*bmag_n[0]/ddpsi_n[0]/ddpsi_n[0] - dxdz[0][2]*dxdz[0][2] - dxdz[1][2]*dxdz[1][2])/R/R;
        // Argument is >= 0 analytically; clamp away roundoff so sqrt does not return NaN.
        dphidtheta = sqrt(fmax(0.0, dphidtheta));
        // Recover sign from exact dphidtheta = F(psi)/R/\grad(psi).
        if (ddtheta_n[2] < 0) {
          dphidtheta = -dphidtheta;
        }

        // AS 2/22/25 It seems that now that we are using interior points,
        // cmag comes out fine without directtly enforcing the condition
        //double dphidtheta = ddtheta_n[2];

        double *gFld_n= gkyl_array_fetch(gk_geom->geo_int.g_ij_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        gFld_n[0] = dxdz[0][0]*dxdz[0][0] + R*R*dxdz[2][0]*dxdz[2][0] + dxdz[1][0]*dxdz[1][0]; 
        gFld_n[1] = alpha_sign*R*R*dxdz[2][0];
        gFld_n[2] = dxdz[0][0]*dxdz[0][2] + R*R*dxdz[2][0]*dphidtheta + dxdz[1][0]*dxdz[1][2];
        gFld_n[3] = R*R; 
        gFld_n[4] = alpha_sign*R*R*dphidtheta;
        gFld_n[5] = dxdz[0][2]*dxdz[0][2] + R*R*dphidtheta*dphidtheta + dxdz[1][2]*dxdz[1][2]; 

        // Calculate cartesian components of bhat
        double *bcartFld_n= gkyl_array_fetch(gk_geom->geo_int.bcart_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        double phi = mc2p_n[PHI_IDX];
        double b3 = 1/sqrt(gFld_n[5]);
        bcartFld_n[0] = b3*(dxdz[0][2]*cos(phi) - R*sin(phi)*dphidtheta);
        bcartFld_n[1] = b3*(dxdz[0][2]*sin(phi) + R*cos(phi)*dphidtheta);
        bcartFld_n[2] = b3*(dxdz[1][2]);

        // Set cartesian components of tangents and duals
        double Z = mc2p_n[Z_IDX];
        double *tanvecFld_n= gkyl_array_fetch(gk_geom->geo_int.dxdz_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        tanvecFld_n[0] = dxdz[0][0]*cos(phi) - R*sin(phi)*dxdz[2][0]; 
        tanvecFld_n[1] = dxdz[0][0]*sin(phi)  + R*cos(phi)*dxdz[2][0]; 
        tanvecFld_n[2] = dxdz[1][0];

        tanvecFld_n[3] = alpha_sign*(-R*sin(phi));
        tanvecFld_n[4] = alpha_sign*(+R*cos(phi));
        tanvecFld_n[5] = 0.0; 

        tanvecFld_n[6] = dxdz[0][2]*cos(phi) - R*sin(phi)*dphidtheta; 
        tanvecFld_n[7] = dxdz[0][2]*sin(phi)  + R*cos(phi)*dphidtheta; 
        tanvecFld_n[8] = dxdz[1][2];

        double *dualFld_n= gkyl_array_fetch(gk_geom->geo_int.dzdx_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        const double e1[3] = { tanvecFld_n[0], tanvecFld_n[1], tanvecFld_n[2] };
        const double e2[3] = { tanvecFld_n[3], tanvecFld_n[4], tanvecFld_n[5] };
        const double e3[3] = { tanvecFld_n[6], tanvecFld_n[7], tanvecFld_n[8] };
        calc_dual(jFld_n[0], e2, e3, &dualFld_n[0]);
        calc_dual(jFld_n[0], e3, e1, &dualFld_n[3]);
        calc_dual(jFld_n[0], e1, e2, &dualFld_n[6]);

        // Check that the coordinate system has tangent/dual vectors
        // satisfying orthonormality, and that it's right handed.
        check_orthonormality(tanvecFld_n, dualFld_n, up->exit_at_checks);
        check_right_handed(tanvecFld_n, dualFld_n, up->exit_at_checks);

        double norm1 = sqrt(dualFld_n[0]*dualFld_n[0] + dualFld_n[1]*dualFld_n[1] + dualFld_n[2]*dualFld_n[2]);
        double norm2 = sqrt(dualFld_n[3]*dualFld_n[3] + dualFld_n[4]*dualFld_n[4] + dualFld_n[5]*dualFld_n[5]);
        double norm3 = sqrt(dualFld_n[6]*dualFld_n[6] + dualFld_n[7]*dualFld_n[7] + dualFld_n[8]*dualFld_n[8]);

        double *dualmagFld_n = gkyl_array_fetch(gk_geom->geo_int.dualmag_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        dualmagFld_n[0] = norm1;
        dualmagFld_n[1] = norm2;
        dualmagFld_n[2] = norm3;
        
        // Set normal vectors
        double *normFld_n = gkyl_array_fetch(gk_geom->geo_int.normals_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        normFld_n[0] = dualFld_n[0]/norm1;
        normFld_n[1] = dualFld_n[1]/norm1;
        normFld_n[2] = dualFld_n[2]/norm1;

        normFld_n[3] = dualFld_n[3]/norm2;
        normFld_n[4] = dualFld_n[4]/norm2;
        normFld_n[5] = dualFld_n[5]/norm2;

        normFld_n[6] = dualFld_n[6]/norm3;
        normFld_n[7] = dualFld_n[7]/norm3;
        normFld_n[8] = dualFld_n[8]/norm3;

        // Set e^m \dot curl(bhat) 
        double *curlbhat_n = gkyl_array_fetch(gk_geom->geo_int.curlbhat_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        double *dualcurlbhat_n = gkyl_array_fetch(gk_geom->geo_int.dualcurlbhat_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        dualcurlbhat_n[0] = dualFld_n[0]*curlbhat_n[0] +  dualFld_n[1]*curlbhat_n[1] + dualFld_n[2]*curlbhat_n[2];
        dualcurlbhat_n[1] = dualFld_n[3]*curlbhat_n[0] +  dualFld_n[4]*curlbhat_n[1] + dualFld_n[5]*curlbhat_n[2];
        dualcurlbhat_n[2] = dualFld_n[6]*curlbhat_n[0] +  dualFld_n[7]*curlbhat_n[1] + dualFld_n[8]*curlbhat_n[2];


        // Set e^3 \dot B 
        double *B3_n = gkyl_array_fetch(gk_geom->geo_int.B3_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        B3_n[0] = bmag_n[0]/sqrt(gFld_n[5]);

        // set e^3 \dot B /|B|
        double *dualcurlbhatoverB_n = gkyl_array_fetch(gk_geom->geo_int.dualcurlbhatoverB_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        dualcurlbhatoverB_n[0] = dualcurlbhat_n[0]/bmag_n[0];
        dualcurlbhatoverB_n[1] = dualcurlbhat_n[1]/bmag_n[0];
        dualcurlbhatoverB_n[2] = dualcurlbhat_n[2]/bmag_n[0];

        // set B^3/B = 1/sqrt(g_33)
        double *rtg33inv_n = gkyl_array_fetch(gk_geom->geo_int.rtg33inv_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        rtg33inv_n[0] = 1.0/sqrt(gFld_n[5]);

        // set b_i/JB
        double* bioverJB_n = gkyl_array_fetch(gk_geom->geo_int.bioverJB_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        bioverJB_n[0] = gFld_n[2]/sqrt(gFld_n[5])/J/bmag_n[0];
        bioverJB_n[1] = gFld_n[4]/sqrt(gFld_n[5])/J/bmag_n[0];
        bioverJB_n[2] = gFld_n[5]/sqrt(gFld_n[5])/J/bmag_n[0];

      }
    }
  }
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &gk_geom->local, 6, gk_geom->geo_int.g_ij_nodal, gk_geom->geo_int.g_ij, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &gk_geom->local, 1, gk_geom->geo_int.jacobgeo_nodal, gk_geom->geo_int.jacobgeo, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &gk_geom->local, 3, gk_geom->geo_int.bcart_nodal, gk_geom->geo_int.bcart, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &gk_geom->local, 9, gk_geom->geo_int.dxdz_nodal, gk_geom->geo_int.dxdz, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &gk_geom->local, 9, gk_geom->geo_int.dzdx_nodal, gk_geom->geo_int.dzdx, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &gk_geom->local, 3, gk_geom->geo_int.dualmag_nodal, gk_geom->geo_int.dualmag, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &gk_geom->local, 9, gk_geom->geo_int.normals_nodal, gk_geom->geo_int.normals, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &gk_geom->local, 3, gk_geom->geo_int.dualcurlbhat_nodal, gk_geom->geo_int.dualcurlbhat, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &gk_geom->local, 3, gk_geom->geo_int.dualcurlbhatoverB_nodal, gk_geom->geo_int.dualcurlbhatoverB, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &gk_geom->local, 1, gk_geom->geo_int.rtg33inv_nodal, gk_geom->geo_int.rtg33inv, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &gk_geom->local, 3, gk_geom->geo_int.bioverJB_nodal, gk_geom->geo_int.bioverJB, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &gk_geom->local, 1, gk_geom->geo_int.B3_nodal, gk_geom->geo_int.B3, true);
  up->rz_expected_jacobian_sign = jac_guard.expected_sign;
  up->rz_jacobian_valid = jac_guard.valid;
  metric_diag_write_interior(gk_geom);
}

void gkyl_calc_metric_advance_rz_surface(gkyl_calc_metric *up, int dir, struct gk_geometry *gk_geom)
{
  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { R_IDX, Z_IDX, PHI_IDX }; // arrangement of cartesian coordinates
  int cidx[3];
  struct signed_jacobian_guard_state jac_guard = {
    .expected_sign = up->rz_expected_jacobian_sign,
    .valid = true,
    .nonfatal = up->rz_jacobian_guard_nonfatal,
  };
  for(int ia=gk_geom->nrange_surf[dir].lower[AL_IDX]; ia<=gk_geom->nrange_surf[dir].upper[AL_IDX]; ++ia){
    for (int ip=gk_geom->nrange_surf[dir].lower[PSI_IDX]; ip<=gk_geom->nrange_surf[dir].upper[PSI_IDX]; ++ip) {
      for (int it=gk_geom->nrange_surf[dir].lower[TH_IDX]; it<=gk_geom->nrange_surf[dir].upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        const double *mc2p_n = gkyl_array_cfetch(gk_geom->geo_surf[dir].mc2p_nodal_fd, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        double dxdz[3][3];

        if((ip == gk_geom->nrange_surf[dir].lower[PSI_IDX]) && (up->local.lower[PSI_IDX]== up->global.lower[PSI_IDX]) && dir==0) {
          dxdz[0][0] = (-3*mc2p_n[R_IDX] + 4*mc2p_n[6+R_IDX] - mc2p_n[12+R_IDX] )/gk_geom->dzc[0]/2;
          dxdz[1][0] = (-3*mc2p_n[Z_IDX] + 4*mc2p_n[6+Z_IDX] - mc2p_n[12+Z_IDX] )/gk_geom->dzc[0]/2;
          dxdz[2][0] = (-3*mc2p_n[PHI_IDX] + 4*mc2p_n[6+PHI_IDX] - mc2p_n[12+PHI_IDX] );
        }
        else if((ip == gk_geom->nrange_surf[dir].upper[PSI_IDX]) && (up->local.upper[PSI_IDX]== up->global.upper[PSI_IDX]) && dir==0) {
          dxdz[0][0] = (3*mc2p_n[R_IDX] - 4*mc2p_n[3+R_IDX] + mc2p_n[9+R_IDX] )/gk_geom->dzc[0]/2;
          dxdz[1][0] = (3*mc2p_n[Z_IDX] - 4*mc2p_n[3+Z_IDX] + mc2p_n[9+Z_IDX] )/gk_geom->dzc[0]/2;
          dxdz[2][0] = (3*mc2p_n[PHI_IDX] - 4*mc2p_n[3+PHI_IDX] + mc2p_n[9+PHI_IDX] );
        }
        else {
          dxdz[0][0] = -(mc2p_n[3 +R_IDX] -   mc2p_n[6+R_IDX])/2/gk_geom->dzc[0];
          dxdz[1][0] = -(mc2p_n[3 +Z_IDX] -   mc2p_n[6+Z_IDX])/2/gk_geom->dzc[0];
          dxdz[2][0] = -(mc2p_n[3 +PHI_IDX] -   mc2p_n[6+PHI_IDX]);
        }
        // Take into account wrapping of cyclic coordinate phi
        if (dxdz[2][0] < -M_PI) {
          dxdz[2][0] += 2*M_PI;
        } else if (dxdz[2][0] > M_PI) {
          dxdz[2][0] -= 2*M_PI;
        }
        dxdz[2][0] = dxdz[2][0]/2.0/gk_geom->dzc[0];

        // Use exact expressions for dR/dtheta and dZ/dtheta
        double *ddtheta_n = gkyl_array_fetch(gk_geom->geo_surf[dir].ddtheta_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        dxdz[0][2] = ddtheta_n[0]; // dR/dtheta
        dxdz[1][2] = ddtheta_n[1]; // dZ/dtheta

        // Get position map deriv for psi
        double *ddpsi_n = gkyl_array_fetch(gk_geom->geo_surf[dir].ddpsi_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));

        double *bmag_n = gkyl_array_fetch(gk_geom->geo_surf[dir].bmag_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        double R = mc2p_n[R_IDX];

        // At a physical psi boundary, the one-sided finite-difference
        // stencil can acquire a large, non-smooth component tangent to a
        // separatrix contour.  Retain that tangential component (it is set
        // by the theta coordinate), but obtain the normal component from
        // the independently evaluated EQDSK quantities.  This enforces
        //
        //   J |B| / (dPsi/dpsi) = sqrt(g_33)
        //
        // without changing the requested psi/rho extent.
        bool psi_boundary = dir == PSI_IDX &&
          ((ip == gk_geom->nrange_surf[dir].lower[PSI_IDX] &&
              up->local.lower[PSI_IDX] == up->global.lower[PSI_IDX]) ||
           (ip == gk_geom->nrange_surf[dir].upper[PSI_IDX] &&
              up->local.upper[PSI_IDX] == up->global.upper[PSI_IDX]));
        if (psi_boundary) {
          double tangent_sq = dxdz[0][2]*dxdz[0][2] + dxdz[1][2]*dxdz[1][2];
          double g33_exact = tangent_sq + R*R*ddtheta_n[2]*ddtheta_n[2];
          double cross_rz = dxdz[0][0]*dxdz[1][2] - dxdz[0][2]*dxdz[1][0];
          // Preserve the orientation established by the block's interior
          // mapping.  The raw one-sided radial stencil supplies only the
          // tangential projection here; its sign can be noisy at a
          // separatrix and must not flip surface metrics or normals.
          double cross_sign = up->rz_expected_jacobian_sign != 0
            ? up->rz_expected_jacobian_sign
            : (cross_rz < 0.0 ? -1.0 : 1.0);
          double cross_exact = cross_sign*fabs(ddpsi_n[0])*sqrt(g33_exact)/(R*bmag_n[0]);
          double tangent_projection =
            (dxdz[0][0]*dxdz[0][2] + dxdz[1][0]*dxdz[1][2])/tangent_sq;
          dxdz[0][0] = tangent_projection*dxdz[0][2] + cross_exact*dxdz[1][2]/tangent_sq;
          dxdz[1][0] = tangent_projection*dxdz[1][2] - cross_exact*dxdz[0][2]/tangent_sq;
        }

        // dxdz is in cylindrical coords, calculate J as
        // J = R(dR/dpsi*dZ/dtheta - dR/dtheta*dZ/dpsi)
        long surf_idx =
          gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx);

        double *jFld_n = gkyl_array_fetch(
          gk_geom->geo_surf[dir].jacobgeo_nodal,
          surf_idx
        );

        double *jSignedFld_n = gkyl_array_fetch(
          gk_geom->geo_surf[dir].jacobgeo_signed_nodal,
          surf_idx
        );

        double J = signed_jacobian_rz(
          R,
          dxdz[0][0], dxdz[1][0],
          dxdz[0][2], dxdz[1][2]
        );

        signed_jacobian_guard_check(
          &jac_guard, J, true, "surface", dir, cidx
        );

        // Diagnostic orientation-preserving value.
        jSignedFld_n[0] = J;

        // Existing positive physical volume factor.
        jFld_n[0] = fabs(J);
        double alpha_sign = J < 0.0 ? -1.0 : 1.0;

        // Calculate dphi/dtheta based on the divergence free condition
        // on B: 1 = J*B/sqrt(g_33)
        // Argument is >= 0 analytically; clamp away roundoff so sqrt does not return NaN.
        double dphidtheta = psi_boundary ? ddtheta_n[2] :
          sqrt(fmax(0.0, (jFld_n[0]*jFld_n[0]*bmag_n[0]*bmag_n[0]/ddpsi_n[0]/ddpsi_n[0]
            - dxdz[0][2]*dxdz[0][2] - dxdz[1][2]*dxdz[1][2])/R/R));
        // Recover sign from exact dphidtheta = F(psi)/R/\grad(psi).
        if (!psi_boundary && ddtheta_n[2] < 0)
          dphidtheta = -dphidtheta;

        double *gFld_n= gkyl_array_fetch(gk_geom->geo_surf[dir].g_ij_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        gFld_n[0] = dxdz[0][0]*dxdz[0][0] + R*R*dxdz[2][0]*dxdz[2][0] + dxdz[1][0]*dxdz[1][0]; 
        gFld_n[1] = alpha_sign*R*R*dxdz[2][0];
        gFld_n[2] = dxdz[0][0]*dxdz[0][2] + R*R*dxdz[2][0]*dphidtheta + dxdz[1][0]*dxdz[1][2];
        gFld_n[3] = R*R; 
        gFld_n[4] = alpha_sign*R*R*dphidtheta;
        gFld_n[5] = dxdz[0][2]*dxdz[0][2] + R*R*dphidtheta*dphidtheta + dxdz[1][2]*dxdz[1][2]; 

        // Calculate cmag, bi, and jtot_inv
        double *biFld_n= gkyl_array_fetch(gk_geom->geo_surf[dir].b_i_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        biFld_n[0] = gFld_n[2]/sqrt(gFld_n[5]);
        biFld_n[1] = gFld_n[4]/sqrt(gFld_n[5]);
        biFld_n[2] = gFld_n[5]/sqrt(gFld_n[5]);

        double *cmagFld_n= gkyl_array_fetch(gk_geom->geo_surf[dir].cmag_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        cmagFld_n[0] = jFld_n[0]*bmag_n[0]/sqrt(gFld_n[5]);
        double *jtotinvFld_n= gkyl_array_fetch(gk_geom->geo_surf[dir].jacobtot_inv_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        jtotinvFld_n[0] = 1.0/(jFld_n[0]*bmag_n[0]);

        // Calculate cartesian components of bhat
        double *bcartFld_n= gkyl_array_fetch(gk_geom->geo_surf[dir].bcart_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        double phi = mc2p_n[PHI_IDX];
        double b3 = 1/sqrt(gFld_n[5]);
        bcartFld_n[0] = b3*(dxdz[0][2]*cos(phi) - R*sin(phi)*dphidtheta);
        bcartFld_n[1] = b3*(dxdz[0][2]*sin(phi) + R*cos(phi)*dphidtheta);
        bcartFld_n[2] = b3*(dxdz[1][2]);

        // Set cartesian components of tangents and duals
        double Z = mc2p_n[Z_IDX];
        double *tanvecFld_n= gkyl_array_fetch(gk_geom->geo_surf[dir].dxdz_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        tanvecFld_n[0] = dxdz[0][0]*cos(phi) - R*sin(phi)*dxdz[2][0]; 
        tanvecFld_n[1] = dxdz[0][0]*sin(phi)  + R*cos(phi)*dxdz[2][0]; 
        tanvecFld_n[2] = dxdz[1][0];

        tanvecFld_n[3] = alpha_sign*(-R*sin(phi));
        tanvecFld_n[4] = alpha_sign*(+R*cos(phi));
        tanvecFld_n[5] = 0.0; 

        tanvecFld_n[6] = dxdz[0][2]*cos(phi) - R*sin(phi)*dphidtheta; 
        tanvecFld_n[7] = dxdz[0][2]*sin(phi)  + R*cos(phi)*dphidtheta; 
        tanvecFld_n[8] = dxdz[1][2];

        double *dualFld_n= gkyl_array_fetch(gk_geom->geo_surf[dir].dzdx_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        const double e1[3] = { tanvecFld_n[0], tanvecFld_n[1], tanvecFld_n[2] };
        const double e2[3] = { tanvecFld_n[3], tanvecFld_n[4], tanvecFld_n[5] };
        const double e3[3] = { tanvecFld_n[6], tanvecFld_n[7], tanvecFld_n[8] };
        calc_dual(jFld_n[0], e2, e3, &dualFld_n[0]);
        calc_dual(jFld_n[0], e3, e1, &dualFld_n[3]);
        calc_dual(jFld_n[0], e1, e2, &dualFld_n[6]);

        // Check that the coordinate system has tangent/dual vectors
        // satisfying orthonormality, and that it's right handed.
        check_orthonormality(tanvecFld_n, dualFld_n, up->exit_at_checks);
        check_right_handed(tanvecFld_n, dualFld_n, up->exit_at_checks);

        double norm1 = sqrt(dualFld_n[0]*dualFld_n[0] + dualFld_n[1]*dualFld_n[1] + dualFld_n[2]*dualFld_n[2]);
        double norm2 = sqrt(dualFld_n[3]*dualFld_n[3] + dualFld_n[4]*dualFld_n[4] + dualFld_n[5]*dualFld_n[5]);
        double norm3 = sqrt(dualFld_n[6]*dualFld_n[6] + dualFld_n[7]*dualFld_n[7] + dualFld_n[8]*dualFld_n[8]);

        double *dualmagFld_n = gkyl_array_fetch(gk_geom->geo_surf[dir].dualmag_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        dualmagFld_n[0] = norm1;
        dualmagFld_n[1] = norm2;
        dualmagFld_n[2] = norm3;
        
        // Set normal vectors
        double *normFld_n = gkyl_array_fetch(gk_geom->geo_surf[dir].normals_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        normFld_n[0] = dualFld_n[0]/norm1;
        normFld_n[1] = dualFld_n[1]/norm1;
        normFld_n[2] = dualFld_n[2]/norm1;

        normFld_n[3] = dualFld_n[3]/norm2;
        normFld_n[4] = dualFld_n[4]/norm2;
        normFld_n[5] = dualFld_n[5]/norm2;

        normFld_n[6] = dualFld_n[6]/norm3;
        normFld_n[7] = dualFld_n[7]/norm3;
        normFld_n[8] = dualFld_n[8]/norm3;

        // Set lenr
        double *lenr_n = gkyl_array_fetch(gk_geom->geo_surf[dir].lenr_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        lenr_n[0] = J*dualmagFld_n[dir];

        // Set n^3 \dot B 
        double *B3_n = gkyl_array_fetch(gk_geom->geo_surf[dir].B3_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        B3_n[0] = bmag_n[0]/sqrt(gFld_n[5])/norm3;

        // Set n^m \dot curl(bhat) 
        double *curlbhat_n = gkyl_array_fetch(gk_geom->geo_surf[dir].curlbhat_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        double *normcurlbhat_n = gkyl_array_fetch(gk_geom->geo_surf[dir].normcurlbhat_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        normcurlbhat_n[0] = normFld_n[3*dir+0]*curlbhat_n[0] +  normFld_n[3*dir+1]*curlbhat_n[1] + normFld_n[3*dir+2]*curlbhat_n[2];

        // set bimpactangle = arcsin(1/sqrt(g_33 * g^33))
        double *bimpactangle_n = gkyl_array_fetch(gk_geom->geo_surf[dir].bimpactangle_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        bimpactangle_n[0] = asin(1.0/(sqrt(gFld_n[5]) * norm3));

      }
    }
  }
  up->rz_jacobian_valid = up->rz_jacobian_valid && jac_guard.valid;
  metric_diag_write_surface(dir, gk_geom);
}



void
gkyl_calc_metric_advance_rz_neut_interior( gkyl_calc_metric *up, struct gk_geometry* gk_geom)
{
  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { R_IDX, Z_IDX, PHI_IDX }; // arrangement of cartesian coordinates
  int cidx[3];
  for (int ia=gk_geom->nrange_int.lower[AL_IDX]; ia<=gk_geom->nrange_int.upper[AL_IDX]; ++ia){
    for (int ip=gk_geom->nrange_int.lower[PSI_IDX]; ip<=gk_geom->nrange_int.upper[PSI_IDX]; ++ip) {
      for (int it=gk_geom->nrange_int.lower[TH_IDX]; it<=gk_geom->nrange_int.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        const double *mc2p_n = gkyl_array_cfetch(gk_geom->geo_int.mc2p_nodal_fd, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        double dxdz[3][3];
        dxdz[0][0] = -(mc2p_n[3 +R_IDX] -   mc2p_n[6+R_IDX])/2/gk_geom->dzc[0];
        dxdz[1][0] = -(mc2p_n[3 +Z_IDX] -   mc2p_n[6+Z_IDX])/2/gk_geom->dzc[0];


        // dphi/dpsi =0
        dxdz[2][0] = 0.0;
        // Use exact expressions for dR/dtheta and dZ/dtheta, dphi/dtheta
        double *ddtheta_n = gkyl_array_fetch(gk_geom->geo_int.ddtheta_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        dxdz[0][2] = ddtheta_n[0];
        dxdz[1][2] = ddtheta_n[1];
        dxdz[2][2] = 0.0;
        double dphidtheta = dxdz[2][2];

        // dxdz is in cylindrical coords, calculate J as
        // J = R(dR/dpsi*dZ/dtheta - dR/dtheta*dZ/dpsi)
        double R = mc2p_n[R_IDX];
        double jac = sqrt(R*R*(   dxdz[0][0]*dxdz[0][0]*dxdz[1][2]*dxdz[1][2]
                               +  dxdz[0][2]*dxdz[0][2]*dxdz[1][0]*dxdz[1][0]
                               -2*dxdz[0][0]*dxdz[0][2]*dxdz[1][0]*dxdz[1][2] )) ;

        double *gFld_n= gkyl_array_fetch(gk_geom->geo_int.g_ij_neut_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        gFld_n[0] = dxdz[0][0]*dxdz[0][0] + R*R*dxdz[2][0]*dxdz[2][0] + dxdz[1][0]*dxdz[1][0]; 
        gFld_n[1] = R*R*dxdz[2][0]; 
        gFld_n[2] = dxdz[0][0]*dxdz[0][2] + R*R*dxdz[2][0]*dphidtheta + dxdz[1][0]*dxdz[1][2];
        gFld_n[3] = R*R; 
        gFld_n[4] = R*R*dphidtheta;
        gFld_n[5] = dxdz[0][2]*dxdz[0][2] + R*R*dphidtheta*dphidtheta + dxdz[1][2]*dxdz[1][2]; 

        double *grFld_n= gkyl_array_fetch(gk_geom->geo_int.gij_neut_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        grFld_n[0] = R*R/jac/jac*(dxdz[1][2]*dxdz[1][2] + dxdz[0][2]*dxdz[0][2] );
        grFld_n[1] = 0.0;
        grFld_n[2] = -R*R/jac/jac*(dxdz[0][0]*dxdz[0][2] + dxdz[1][0]*dxdz[1][2] );
        grFld_n[3] = 1/R/R;
        grFld_n[4] = 0.0;
        grFld_n[5] = R*R/jac/jac*(dxdz[0][0]*dxdz[0][0] + dxdz[1][0]*dxdz[1][0] );

      }
    }
  }
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &up->local, 6, gk_geom->geo_int.g_ij_neut_nodal, gk_geom->geo_int.g_ij_neut, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &up->local, 6, gk_geom->geo_int.gij_neut_nodal, gk_geom->geo_int.gij_neut, true);
}

void gkyl_calc_metric_advance_mirror(
  gkyl_calc_metric *up, struct gkyl_range *nrange,
  struct gkyl_array *mc2p_nodal_fd, struct gkyl_array *ddtheta_nodal,
  struct gkyl_array *bmag_nodal, double *dzc, struct gkyl_array *gFld,
  struct gkyl_array *tanvecFld,
  struct gkyl_array *dualFld,
  struct gkyl_array *dualmagFld,
  struct gkyl_array *normFld,
  struct gkyl_array *jFld, struct gkyl_array* bcartFld, const struct gkyl_range *update_range)
{
  struct gkyl_array* gFld_nodal = gkyl_array_new(GKYL_DOUBLE, 6, nrange->volume);
  struct gkyl_array* jFld_nodal = gkyl_array_new(GKYL_DOUBLE, 1, nrange->volume);
  struct gkyl_array* bcartFld_nodal = gkyl_array_new(GKYL_DOUBLE, 3, nrange->volume);
  struct gkyl_array* tanvecFld_nodal = gkyl_array_new(GKYL_DOUBLE, 9, nrange->volume);
  struct gkyl_array* dualFld_nodal = gkyl_array_new(GKYL_DOUBLE, 9, nrange->volume);
  struct gkyl_array* dualmagFld_nodal = gkyl_array_new(GKYL_DOUBLE, 3, nrange->volume);
  struct gkyl_array* normFld_nodal = gkyl_array_new(GKYL_DOUBLE, 9, nrange->volume);
  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { R_IDX, Z_IDX, PHI_IDX }; // arrangement of cartesian coordinates
  int cidx[3];
  for (int ia=nrange->lower[AL_IDX]; ia<=nrange->upper[AL_IDX]; ++ia){
    for (int ip=nrange->lower[PSI_IDX]; ip<=nrange->upper[PSI_IDX]; ++ip) {
      for (int it=nrange->lower[TH_IDX]; it<=nrange->upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        const double *mc2p_n = gkyl_array_cfetch(mc2p_nodal_fd, gkyl_range_idx(nrange, cidx));
        double dxdz[3][3];

        if((ip == nrange->lower[PSI_IDX]) && (up->local.lower[PSI_IDX]== up->global.lower[PSI_IDX]) ) {
          dxdz[0][0] = (-3*mc2p_n[R_IDX] + 4*mc2p_n[6+R_IDX] - mc2p_n[12+R_IDX] )/dzc[0]/2;
          dxdz[1][0] = (-3*mc2p_n[Z_IDX] + 4*mc2p_n[6+Z_IDX] - mc2p_n[12+Z_IDX] )/dzc[0]/2;
          dxdz[2][0] = (-3*mc2p_n[PHI_IDX] + 4*mc2p_n[6+PHI_IDX] - mc2p_n[12+PHI_IDX] )/dzc[0]/2;
        }
        else if((ip == nrange->upper[PSI_IDX]) && (up->local.upper[PSI_IDX]== up->global.upper[PSI_IDX])) {
          dxdz[0][0] = (3*mc2p_n[R_IDX] - 4*mc2p_n[3+R_IDX] + mc2p_n[9+R_IDX] )/dzc[0]/2;
          dxdz[1][0] = (3*mc2p_n[Z_IDX] - 4*mc2p_n[3+Z_IDX] + mc2p_n[9+Z_IDX] )/dzc[0]/2;
          dxdz[2][0] = (3*mc2p_n[PHI_IDX] - 4*mc2p_n[3+PHI_IDX] + mc2p_n[9+PHI_IDX] )/dzc[0]/2;
        }
        else {
          dxdz[0][0] = -(mc2p_n[3 +R_IDX] -   mc2p_n[6+R_IDX])/2/dzc[0];
          dxdz[1][0] = -(mc2p_n[3 +Z_IDX] -   mc2p_n[6+Z_IDX])/2/dzc[0];
          dxdz[2][0] = -(mc2p_n[3 +PHI_IDX] -   mc2p_n[6+PHI_IDX])/2/dzc[0];
        }

        if((ia == nrange->lower[AL_IDX]) && (up->local.lower[AL_IDX]== up->global.lower[AL_IDX]) ) {
          dxdz[0][1] = (-3*mc2p_n[R_IDX] +  4*mc2p_n[18+R_IDX] -  mc2p_n[24+R_IDX])/dzc[1]/2;
          dxdz[1][1] = (-3*mc2p_n[Z_IDX] +  4*mc2p_n[18+Z_IDX] -  mc2p_n[24+Z_IDX])/dzc[1]/2;
          dxdz[2][1] = (-3*mc2p_n[PHI_IDX] +  4*mc2p_n[18+PHI_IDX] -  mc2p_n[24+PHI_IDX])/dzc[1]/2;
        }
        else if((ia == nrange->upper[AL_IDX])  && (up->local.upper[AL_IDX]== up->global.upper[AL_IDX])){
          dxdz[0][1] = (3*mc2p_n[R_IDX] -  4*mc2p_n[15+R_IDX] +  mc2p_n[21+R_IDX] )/dzc[1]/2;
          dxdz[1][1] = (3*mc2p_n[Z_IDX] -  4*mc2p_n[15+Z_IDX] +  mc2p_n[21+Z_IDX] )/dzc[1]/2;
          dxdz[2][1] = (3*mc2p_n[PHI_IDX] -  4*mc2p_n[15+PHI_IDX] +  mc2p_n[21+PHI_IDX] )/dzc[1]/2;
        }
        else {
          dxdz[0][1] = -(mc2p_n[15 +R_IDX] -  mc2p_n[18 +R_IDX])/2/dzc[1];
          dxdz[1][1] = -(mc2p_n[15 +Z_IDX] -  mc2p_n[18 +Z_IDX])/2/dzc[1];
          dxdz[2][1] = -(mc2p_n[15 +PHI_IDX] -  mc2p_n[18 +PHI_IDX])/2/dzc[1];
        }

        if((it == nrange->lower[TH_IDX]) && (up->local.lower[TH_IDX]== up->global.lower[TH_IDX])){
          dxdz[0][2] = (-3*mc2p_n[R_IDX] + 4*mc2p_n[30+R_IDX] - mc2p_n[36+R_IDX])/dzc[2]/2;
          dxdz[1][2] = (-3*mc2p_n[Z_IDX] + 4*mc2p_n[30+Z_IDX] - mc2p_n[36+Z_IDX])/dzc[2]/2;
          dxdz[2][2] = (-3*mc2p_n[PHI_IDX] + 4*mc2p_n[30+PHI_IDX] - mc2p_n[36+PHI_IDX])/dzc[2]/2;
        }
        else if((it == nrange->upper[TH_IDX]) && (up->local.upper[TH_IDX]== up->global.upper[TH_IDX])){
          dxdz[0][2] = (3*mc2p_n[R_IDX] - 4*mc2p_n[27+R_IDX] + mc2p_n[33+R_IDX] )/dzc[2]/2;
          dxdz[1][2] = (3*mc2p_n[Z_IDX] - 4*mc2p_n[27+Z_IDX] + mc2p_n[33+Z_IDX] )/dzc[2]/2;
          dxdz[2][2] = (3*mc2p_n[PHI_IDX] - 4*mc2p_n[27+PHI_IDX] + mc2p_n[33+PHI_IDX] )/dzc[2]/2;
        }
        else {
          dxdz[0][2] = -(mc2p_n[27 +R_IDX] - mc2p_n[30 +R_IDX])/2/dzc[2];
          dxdz[1][2] = -(mc2p_n[27 +Z_IDX] - mc2p_n[30 +Z_IDX])/2/dzc[2];
          dxdz[2][2] = -(mc2p_n[27 +PHI_IDX] - mc2p_n[30 +PHI_IDX])/2/dzc[2];
        }

        // Use exact expressions for dphidtheta, dR/dtheta, and dZ/dtheta
        double *ddtheta_n = gkyl_array_fetch(ddtheta_nodal, gkyl_range_idx(nrange, cidx));
        dxdz[0][2] = ddtheta_n[0]; // dR/dtheta
        dxdz[1][2] = ddtheta_n[1]; // dZ/dtheta
        dxdz[2][2] = ddtheta_n[2]; // dphi/dtheta

        double R = mc2p_n[R_IDX];
        double *bmag_n = gkyl_array_fetch(bmag_nodal, gkyl_range_idx(nrange, cidx));

        double *gFld_n= gkyl_array_fetch(gFld_nodal, gkyl_range_idx(nrange, cidx));
        gFld_n[5] = dxdz[0][2]*dxdz[0][2] + dxdz[1][2]*dxdz[1][2];

        // dxdz is in cylindrical coordinates. Caculate dR/dpsi as
        // dR/dpsi = (1/(dZ/dtheta) ) * [sqrt(g_33)/RB + dR/dtheta*dZ/dtheta]
        // because this is more reliable near R=0 than using finite differences
        double dRdpsi = 1/dxdz[1][2]*(sqrt(gFld_n[5])/bmag_n[0]/R + dxdz[0][2]*dxdz[1][0]);
        // Calculate J as J = R(dR/dpsi*dZ/dtheta - dR/dtheta*dZ/dpsi)
        double *jFld_n= gkyl_array_fetch(jFld_nodal, gkyl_range_idx(nrange, cidx));
        jFld_n[0] = sqrt(R*R*( dRdpsi*dRdpsi*dxdz[1][2]*dxdz[1][2]
                              +dxdz[0][2]*dxdz[0][2]*dxdz[1][0]*dxdz[1][0]
                              -2*dRdpsi*dxdz[0][2]*dxdz[1][0]*dxdz[1][2] ));

        gFld_n[0] = dRdpsi*dRdpsi + dxdz[1][0]*dxdz[1][0]; 
        gFld_n[1] = 0.0; 
        gFld_n[2] = dRdpsi*dxdz[0][2] + dxdz[1][0]*dxdz[1][2];
        gFld_n[3] = R*R; 
        gFld_n[4] = 0.0; 

        // Now do bcart
        double *bcartFld_n= gkyl_array_fetch(bcartFld_nodal, gkyl_range_idx(nrange, cidx));
        double phi = mc2p_n[PHI_IDX];
        double b3 = 1/sqrt(gFld_n[5]);
        bcartFld_n[0] = b3*dxdz[0][2]*cos(phi);
        bcartFld_n[1] = b3*dxdz[0][2]*sin(phi);
        bcartFld_n[2] = b3*dxdz[1][2];

        // Set cartesian components of tangents and duals
        double Z = mc2p_n[Z_IDX];
        double J = jFld_n[0];
        double *tanvecFld_n= gkyl_array_fetch(tanvecFld_nodal, gkyl_range_idx(nrange, cidx));
        tanvecFld_n[0] = dxdz[0][0]*cos(phi) - R*sin(phi)*dxdz[2][0]; 
        tanvecFld_n[1] = dxdz[0][0]*sin(phi)  + R*cos(phi)*dxdz[2][0]; 
        tanvecFld_n[2] = dxdz[1][0];

        tanvecFld_n[3] = +R*sin(phi); 
        tanvecFld_n[4] = -R*cos(phi); 
        tanvecFld_n[5] = 0.0; 

        tanvecFld_n[6] = dxdz[0][2]*cos(phi); 
        tanvecFld_n[7] = dxdz[0][2]*sin(phi); 
        tanvecFld_n[8] = dxdz[1][2];

        double *dualFld_n= gkyl_array_fetch(dualFld_nodal, gkyl_range_idx(nrange, cidx));
        dualFld_n[0] = -R/J*cos(phi)*dxdz[1][2];
        dualFld_n[1] = -R/J*sin(phi)*dxdz[1][2];
        dualFld_n[2] = +R/J*dxdz[0][2];

        dualFld_n[3] =  1/J * (dxdz[1][0]*dxdz[0][2]*sin(phi) - dxdz[1][2]*dxdz[0][0]*sin(phi) - dxdz[1][2]*R*cos(phi)*dxdz[2][0] );
        dualFld_n[4] = -1/J * (dxdz[1][0]*dxdz[0][2]*cos(phi) - dxdz[1][2]*dxdz[0][0]*cos(phi) + dxdz[1][2]*R*sin(phi)*dxdz[2][0] );
        dualFld_n[5] =  R/J * dxdz[0][2]*dxdz[2][0];

        dualFld_n[6] = +R/J*cos(phi)*dxdz[1][0];
        dualFld_n[7] = +R/J*sin(phi)*dxdz[1][0];
        dualFld_n[8] = -R/J*dxdz[0][0];

        double norm1 = sqrt(dualFld_n[0]*dualFld_n[0] + dualFld_n[1]*dualFld_n[1] + dualFld_n[2]*dualFld_n[2]);
        double norm2 = sqrt(dualFld_n[3]*dualFld_n[3] + dualFld_n[4]*dualFld_n[4] + dualFld_n[5]*dualFld_n[5]);
        double norm3 = sqrt(dualFld_n[6]*dualFld_n[6] + dualFld_n[7]*dualFld_n[7] + dualFld_n[8]*dualFld_n[8]);

        double *dualmagFld_n = gkyl_array_fetch(dualmagFld_nodal, gkyl_range_idx(nrange, cidx));
        dualmagFld_n[0] = norm1;
        dualmagFld_n[1] = norm2;
        dualmagFld_n[2] = norm3;
        
        // Set normal vectors
        double *normFld_n = gkyl_array_fetch(normFld_nodal, gkyl_range_idx(nrange, cidx));
        normFld_n[0] = dualFld_n[0]/norm1;
        normFld_n[1] = dualFld_n[1]/norm1;
        normFld_n[2] = dualFld_n[2]/norm1;

        normFld_n[3] = dualFld_n[3]/norm2;
        normFld_n[4] = dualFld_n[4]/norm2;
        normFld_n[5] = dualFld_n[5]/norm2;

        normFld_n[6] = dualFld_n[6]/norm3;
        normFld_n[7] = dualFld_n[7]/norm3;
        normFld_n[8] = dualFld_n[8]/norm3;
      }
    }
  }
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 6, gFld_nodal, gFld, false);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 1, jFld_nodal, jFld, false);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 3, bcartFld_nodal, bcartFld, false);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 9, tanvecFld_nodal, tanvecFld, false);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 9, dualFld_nodal, dualFld, false);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 3, dualmagFld_nodal, dualmagFld, false);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 9, normFld_nodal, normFld, false);
  gkyl_array_release(gFld_nodal);
  gkyl_array_release(jFld_nodal);
  gkyl_array_release(bcartFld_nodal);
  gkyl_array_release(tanvecFld_nodal);
  gkyl_array_release(dualFld_nodal);
  gkyl_array_release(dualmagFld_nodal);
  gkyl_array_release(normFld_nodal);
}

void gkyl_calc_metric_advance_mirror_interior(
  gkyl_calc_metric *up, struct gkyl_range *nrange,
  struct gkyl_array *mc2p_nodal_fd, struct gkyl_array *ddtheta_nodal,
  struct gkyl_array *bmag_nodal, double *dzc, struct gkyl_array *gFld,
  struct gkyl_array *tanvecFld,
  struct gkyl_array *dualFld,
  struct gkyl_array *dualmagFld,
  struct gkyl_array *normFld,
  struct gkyl_array *jFld, struct gkyl_array* bcartFld, const struct gkyl_range *update_range)
{
  struct gkyl_array* gFld_nodal = gkyl_array_new(GKYL_DOUBLE, 6, nrange->volume);
  struct gkyl_array* jFld_nodal = gkyl_array_new(GKYL_DOUBLE, 1, nrange->volume);
  struct gkyl_array* bcartFld_nodal = gkyl_array_new(GKYL_DOUBLE, 3, nrange->volume);
  struct gkyl_array* tanvecFld_nodal = gkyl_array_new(GKYL_DOUBLE, 9, nrange->volume);
  struct gkyl_array* dualFld_nodal = gkyl_array_new(GKYL_DOUBLE, 9, nrange->volume);
  struct gkyl_array* dualmagFld_nodal = gkyl_array_new(GKYL_DOUBLE, 3, nrange->volume);
  struct gkyl_array* normFld_nodal = gkyl_array_new(GKYL_DOUBLE, 9, nrange->volume);
  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { R_IDX, Z_IDX, PHI_IDX }; // arrangement of cartesian coordinates
  int cidx[3];
  for (int ia=nrange->lower[AL_IDX]; ia<=nrange->upper[AL_IDX]; ++ia){
    for (int ip=nrange->lower[PSI_IDX]; ip<=nrange->upper[PSI_IDX]; ++ip) {
      for (int it=nrange->lower[TH_IDX]; it<=nrange->upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        const double *mc2p_n = gkyl_array_cfetch(mc2p_nodal_fd, gkyl_range_idx(nrange, cidx));
        double dxdz[3][3];

        dxdz[0][0] = -(mc2p_n[3 +R_IDX] -   mc2p_n[6+R_IDX])/2/dzc[0];
        dxdz[1][0] = -(mc2p_n[3 +Z_IDX] -   mc2p_n[6+Z_IDX])/2/dzc[0];
        dxdz[2][0] = -(mc2p_n[3 +PHI_IDX] -   mc2p_n[6+PHI_IDX])/2/dzc[0];

        // Use exact expressions for dR/dtheta, dZ/dtheta, and dphidtheta 
        double *ddtheta_n = gkyl_array_fetch(ddtheta_nodal, gkyl_range_idx(nrange, cidx));
        dxdz[0][2] = ddtheta_n[0]; // dR/dtheta
        dxdz[1][2] = ddtheta_n[1]; // dZ/dtheta
        dxdz[2][2] = ddtheta_n[2]; // dphi/dtheta

        double R = mc2p_n[R_IDX];
        double *bmag_n = gkyl_array_fetch(bmag_nodal, gkyl_range_idx(nrange, cidx));

        double *gFld_n= gkyl_array_fetch(gFld_nodal, gkyl_range_idx(nrange, cidx));
        gFld_n[5] = dxdz[0][2]*dxdz[0][2] + dxdz[1][2]*dxdz[1][2];

        // dxdz is in cylindrical coordinates. Caculate dR/dpsi as
        // dR/dpsi = (1/(dZ/dtheta) ) * [sqrt(g_33)/RB + dR/dtheta*dZ/dtheta]
        // because this is more reliable near R=0 than using finite differences
        double dRdpsi = 1/dxdz[1][2]*(sqrt(gFld_n[5])/bmag_n[0]/R + dxdz[0][2]*dxdz[1][0]);
        // Calculate J as J = R(dR/dpsi*dZ/dtheta - dR/dtheta*dZ/dpsi)
        double *jFld_n= gkyl_array_fetch(jFld_nodal, gkyl_range_idx(nrange, cidx));
        jFld_n[0] = sqrt(R*R*( dRdpsi*dRdpsi*dxdz[1][2]*dxdz[1][2]
                              +dxdz[0][2]*dxdz[0][2]*dxdz[1][0]*dxdz[1][0]
                              -2*dRdpsi*dxdz[0][2]*dxdz[1][0]*dxdz[1][2] ));

        gFld_n[0] = dRdpsi*dRdpsi + dxdz[1][0]*dxdz[1][0]; 
        gFld_n[1] = 0.0; 
        gFld_n[2] = dRdpsi*dxdz[0][2] + dxdz[1][0]*dxdz[1][2];
        gFld_n[3] = R*R; 
        gFld_n[4] = 0.0; 

        // Now do bcart
        double *bcartFld_n= gkyl_array_fetch(bcartFld_nodal, gkyl_range_idx(nrange, cidx));
        double phi = mc2p_n[PHI_IDX];
        double b3 = 1/sqrt(gFld_n[5]);
        bcartFld_n[0] = b3*dxdz[0][2]*cos(phi);
        bcartFld_n[1] = b3*dxdz[0][2]*sin(phi);
        bcartFld_n[2] = b3*dxdz[1][2];

        // Set cartesian components of tangents and duals
        double Z = mc2p_n[Z_IDX];
        double J = jFld_n[0];
        double *tanvecFld_n= gkyl_array_fetch(tanvecFld_nodal, gkyl_range_idx(nrange, cidx));
        tanvecFld_n[0] = dxdz[0][0]*cos(phi) - R*sin(phi)*dxdz[2][0]; 
        tanvecFld_n[1] = dxdz[0][0]*sin(phi)  + R*cos(phi)*dxdz[2][0]; 
        tanvecFld_n[2] = dxdz[1][0];

        tanvecFld_n[3] = +R*sin(phi); 
        tanvecFld_n[4] = -R*cos(phi); 
        tanvecFld_n[5] = 0.0; 

        tanvecFld_n[6] = dxdz[0][2]*cos(phi); 
        tanvecFld_n[7] = dxdz[0][2]*sin(phi); 
        tanvecFld_n[8] = dxdz[1][2];

        double *dualFld_n= gkyl_array_fetch(dualFld_nodal, gkyl_range_idx(nrange, cidx));
        dualFld_n[0] = -R/J*cos(phi)*dxdz[1][2];
        dualFld_n[1] = -R/J*sin(phi)*dxdz[1][2];
        dualFld_n[2] = +R/J*dxdz[0][2];

        dualFld_n[3] =  1/J * (dxdz[1][0]*dxdz[0][2]*sin(phi) - dxdz[1][2]*dxdz[0][0]*sin(phi) - dxdz[1][2]*R*cos(phi)*dxdz[2][0] );
        dualFld_n[4] = -1/J * (dxdz[1][0]*dxdz[0][2]*cos(phi) - dxdz[1][2]*dxdz[0][0]*cos(phi) - dxdz[1][2]*R*sin(phi)*dxdz[2][0] );
        dualFld_n[5] =  R/J * dxdz[0][2]*dxdz[2][0];

        dualFld_n[6] = +R/J*cos(phi)*dxdz[1][0];
        dualFld_n[7] = +R/J*sin(phi)*dxdz[1][0];
        dualFld_n[8] = -R/J*dxdz[0][0];

        double norm1 = sqrt(dualFld_n[0]*dualFld_n[0] + dualFld_n[1]*dualFld_n[1] + dualFld_n[2]*dualFld_n[2]);
        double norm2 = sqrt(dualFld_n[3]*dualFld_n[3] + dualFld_n[4]*dualFld_n[4] + dualFld_n[5]*dualFld_n[5]);
        double norm3 = sqrt(dualFld_n[6]*dualFld_n[6] + dualFld_n[7]*dualFld_n[7] + dualFld_n[8]*dualFld_n[8]);

        double *dualmagFld_n = gkyl_array_fetch(dualmagFld_nodal, gkyl_range_idx(nrange, cidx));
        dualmagFld_n[0] = norm1;
        dualmagFld_n[1] = norm2;
        dualmagFld_n[2] = norm3;
        
        // Set normal vectors
        double *normFld_n = gkyl_array_fetch(normFld_nodal, gkyl_range_idx(nrange, cidx));
        normFld_n[0] = dualFld_n[0]/norm1;
        normFld_n[1] = dualFld_n[1]/norm1;
        normFld_n[2] = dualFld_n[2]/norm1;

        normFld_n[3] = dualFld_n[3]/norm2;
        normFld_n[4] = dualFld_n[4]/norm2;
        normFld_n[5] = dualFld_n[5]/norm2;

        normFld_n[6] = dualFld_n[6]/norm3;
        normFld_n[7] = dualFld_n[7]/norm3;
        normFld_n[8] = dualFld_n[8]/norm3;
      }
    }
  }
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 6, gFld_nodal, gFld, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 1, jFld_nodal, jFld, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 3, bcartFld_nodal, bcartFld, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 9, tanvecFld_nodal, tanvecFld, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 9, dualFld_nodal, dualFld, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 3, dualmagFld_nodal, dualmagFld, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 9, normFld_nodal, normFld, true);
  gkyl_array_release(gFld_nodal);
  gkyl_array_release(jFld_nodal);
  gkyl_array_release(bcartFld_nodal);
  gkyl_array_release(tanvecFld_nodal);
  gkyl_array_release(dualFld_nodal);
  gkyl_array_release(dualmagFld_nodal);
  gkyl_array_release(normFld_nodal);
}

void gkyl_calc_metric_advance_mirror_surface(
  gkyl_calc_metric *up, int dir, struct gkyl_range *nrange,
  struct gkyl_array *mc2p_nodal_fd, struct gkyl_array *ddtheta_nodal,
  struct gkyl_array *bmag_nodal, double *dzc,
  struct gkyl_array *jFld_nodal,
  struct gkyl_array *biFld_nodal,
  struct gkyl_array *cmagFld_nodal,
  struct gkyl_array *jtotinvFld_nodal,
  const struct gkyl_range *update_range)
{
  struct gkyl_array* gFld_nodal = gkyl_array_new(GKYL_DOUBLE, 6, nrange->volume);
  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { R_IDX, Z_IDX, PHI_IDX }; // arrangement of cartesian coordinates
  int cidx[3];
  for(int ia=nrange->lower[AL_IDX]; ia<=nrange->upper[AL_IDX]; ++ia){
    for (int ip=nrange->lower[PSI_IDX]; ip<=nrange->upper[PSI_IDX]; ++ip) {
      for (int it=nrange->lower[TH_IDX]; it<=nrange->upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        const double *mc2p_n = gkyl_array_cfetch(mc2p_nodal_fd, gkyl_range_idx(nrange, cidx));
        double dxdz[3][3];

        if((ip == nrange->lower[PSI_IDX]) && (up->local.lower[PSI_IDX]== up->global.lower[PSI_IDX]) && dir==0) {
          dxdz[0][0] = (-3*mc2p_n[R_IDX] + 4*mc2p_n[6+R_IDX] - mc2p_n[12+R_IDX] )/dzc[0]/2;
          dxdz[1][0] = (-3*mc2p_n[Z_IDX] + 4*mc2p_n[6+Z_IDX] - mc2p_n[12+Z_IDX] )/dzc[0]/2;
          dxdz[2][0] = (-3*mc2p_n[PHI_IDX] + 4*mc2p_n[6+PHI_IDX] - mc2p_n[12+PHI_IDX] )/dzc[0]/2;
        }
        else if((ip == nrange->upper[PSI_IDX]) && (up->local.upper[PSI_IDX]== up->global.upper[PSI_IDX]) && dir==0) {
          dxdz[0][0] = (3*mc2p_n[R_IDX] - 4*mc2p_n[3+R_IDX] + mc2p_n[9+R_IDX] )/dzc[0]/2;
          dxdz[1][0] = (3*mc2p_n[Z_IDX] - 4*mc2p_n[3+Z_IDX] + mc2p_n[9+Z_IDX] )/dzc[0]/2;
          dxdz[2][0] = (3*mc2p_n[PHI_IDX] - 4*mc2p_n[3+PHI_IDX] + mc2p_n[9+PHI_IDX] )/dzc[0]/2;
        }
        else {
          dxdz[0][0] = -(mc2p_n[3 +R_IDX] -   mc2p_n[6+R_IDX])/2/dzc[0];
          dxdz[1][0] = -(mc2p_n[3 +Z_IDX] -   mc2p_n[6+Z_IDX])/2/dzc[0];
          dxdz[2][0] = -(mc2p_n[3 +PHI_IDX] -   mc2p_n[6+PHI_IDX])/2/dzc[0];
        }

        // Use exact expressions for dR/dtheta, dZ/dtheta, and dphidtheta 
        double *ddtheta_n = gkyl_array_fetch(ddtheta_nodal, gkyl_range_idx(nrange, cidx));
        dxdz[0][2] = ddtheta_n[0]; // dR/dtheta
        dxdz[1][2] = ddtheta_n[1]; // dZ/dtheta
        dxdz[2][2] = ddtheta_n[2]; // dphi/dtheta

        double R = mc2p_n[R_IDX];
        double *bmag_n = gkyl_array_fetch(bmag_nodal, gkyl_range_idx(nrange, cidx));

        double *gFld_n= gkyl_array_fetch(gFld_nodal, gkyl_range_idx(nrange, cidx));
        gFld_n[5] = dxdz[0][2]*dxdz[0][2] + dxdz[1][2]*dxdz[1][2];

        // dxdz is in cylindrical coordinates. Caculate dR/dpsi as
        // dR/dpsi = (1/(dZ/dtheta) ) * [sqrt(g_33)/RB + dR/dtheta*dZ/dtheta]
        // because this is more reliable near R=0 than using finite differences
        double dRdpsi = 1/dxdz[1][2]*(sqrt(gFld_n[5])/bmag_n[0]/R + dxdz[0][2]*dxdz[1][0]);
        // Calculate J as J = R(dR/dpsi*dZ/dtheta - dR/dtheta*dZ/dpsi)
        double *jFld_n= gkyl_array_fetch(jFld_nodal, gkyl_range_idx(nrange, cidx));
        jFld_n[0] = sqrt(R*R*( dRdpsi*dRdpsi*dxdz[1][2]*dxdz[1][2]
                              +dxdz[0][2]*dxdz[0][2]*dxdz[1][0]*dxdz[1][0]
                              -2*dRdpsi*dxdz[0][2]*dxdz[1][0]*dxdz[1][2] ));

        gFld_n[0] = dRdpsi*dRdpsi + dxdz[1][0]*dxdz[1][0]; 
        gFld_n[1] = 0.0; 
        gFld_n[2] = dRdpsi*dxdz[0][2] + dxdz[1][0]*dxdz[1][2];
        gFld_n[3] = R*R; 
        gFld_n[4] = 0.0; 

        // Calculate cmag, bi, and jtot_inv
        double *biFld_n= gkyl_array_fetch(biFld_nodal, gkyl_range_idx(nrange, cidx));
        biFld_n[0] = gFld_n[2]/sqrt(gFld_n[5]);
        biFld_n[1] = gFld_n[4]/sqrt(gFld_n[5]);
        biFld_n[2] = gFld_n[5]/sqrt(gFld_n[5]);

        double *cmagFld_n= gkyl_array_fetch(cmagFld_nodal, gkyl_range_idx(nrange, cidx));
        cmagFld_n[0] = jFld_n[0]*bmag_n[0]/sqrt(gFld_n[5]);
        double *jtotinvFld_n= gkyl_array_fetch(jtotinvFld_nodal, gkyl_range_idx(nrange, cidx));
        jtotinvFld_n[0] = 1.0/(jFld_n[0]*bmag_n[0]);
      }
    }
  }
  gkyl_array_release(gFld_nodal);
}

void gkyl_calc_metric_advance(gkyl_calc_metric *up, struct gkyl_range *nrange,
  struct gkyl_array *mc2p_nodal_fd, double *dzc,
  struct gkyl_array *gFld,
  struct gkyl_array *tanvecFld,
  struct gkyl_array *dualFld,
  struct gkyl_array *dualmagFld,
  struct gkyl_array *normFld,
  const struct gkyl_range *update_range)
{
  struct gkyl_array* gFld_nodal = gkyl_array_new(GKYL_DOUBLE, 6, nrange->volume);
  struct gkyl_array* tanvecFld_nodal = gkyl_array_new(GKYL_DOUBLE, 9, nrange->volume);
  struct gkyl_array* dualFld_nodal = gkyl_array_new(GKYL_DOUBLE, 9, nrange->volume);
  struct gkyl_array* dualmagFld_nodal = gkyl_array_new(GKYL_DOUBLE, 3, nrange->volume);
  struct gkyl_array* normFld_nodal = gkyl_array_new(GKYL_DOUBLE, 9, nrange->volume);
  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates
  int cidx[3];
  
  for (int ia=nrange->lower[AL_IDX]; ia<=nrange->upper[AL_IDX]; ++ia) {
    for (int ip=nrange->lower[PSI_IDX]; ip<=nrange->upper[PSI_IDX]; ++ip) {
      for (int it=nrange->lower[TH_IDX]; it<=nrange->upper[TH_IDX]; ++it) {
          cidx[PSI_IDX] = ip;
          cidx[AL_IDX] = ia;
          cidx[TH_IDX] = it;
          const double *mc2p_n = gkyl_array_cfetch(mc2p_nodal_fd, gkyl_range_idx(nrange, cidx));
          double dxdz[3][3]; // tan vecs at node
          double dzdx[3][3]; // duals at node

          if((ip == nrange->lower[PSI_IDX]) && (up->local.lower[PSI_IDX]== up->global.lower[PSI_IDX]) ) {
            dxdz[0][0] = (-3*mc2p_n[X_IDX] + 4*mc2p_n[6+X_IDX] - mc2p_n[12+X_IDX] )/dzc[0]/2;
            dxdz[1][0] = (-3*mc2p_n[Y_IDX] + 4*mc2p_n[6+Y_IDX] - mc2p_n[12+Y_IDX] )/dzc[0]/2;
            dxdz[2][0] = (-3*mc2p_n[Z_IDX] + 4*mc2p_n[6+Z_IDX] - mc2p_n[12+Z_IDX] )/dzc[0]/2;
          }
          else if((ip == nrange->upper[PSI_IDX]) && (up->local.upper[PSI_IDX]== up->global.upper[PSI_IDX])) {
            dxdz[0][0] = (3*mc2p_n[X_IDX] - 4*mc2p_n[3+X_IDX] + mc2p_n[9+X_IDX] )/dzc[0]/2;
            dxdz[1][0] = (3*mc2p_n[Y_IDX] - 4*mc2p_n[3+Y_IDX] + mc2p_n[9+Y_IDX] )/dzc[0]/2;
            dxdz[2][0] = (3*mc2p_n[Z_IDX] - 4*mc2p_n[3+Z_IDX] + mc2p_n[9+Z_IDX] )/dzc[0]/2;
          }
          else{
            dxdz[0][0] = -(mc2p_n[3 +X_IDX] -   mc2p_n[6+X_IDX])/2/dzc[0];
            dxdz[1][0] = -(mc2p_n[3 +Y_IDX] -   mc2p_n[6+Y_IDX])/2/dzc[0];
            dxdz[2][0] = -(mc2p_n[3 +Z_IDX] -   mc2p_n[6+Z_IDX])/2/dzc[0];
          }


          if((ia == nrange->lower[AL_IDX]) && (up->local.lower[AL_IDX]== up->global.lower[AL_IDX]) ) {
            dxdz[0][1] = (-3*mc2p_n[X_IDX] +  4*mc2p_n[18+X_IDX] -  mc2p_n[24+X_IDX])/dzc[1]/2;
            dxdz[1][1] = (-3*mc2p_n[Y_IDX] +  4*mc2p_n[18+Y_IDX] -  mc2p_n[24+Y_IDX])/dzc[1]/2;
            dxdz[2][1] = (-3*mc2p_n[Z_IDX] +  4*mc2p_n[18+Z_IDX] -  mc2p_n[24+Z_IDX])/dzc[1]/2;
          }
          else if((ia == nrange->upper[AL_IDX])  && (up->local.upper[AL_IDX]== up->global.upper[AL_IDX])) {
            dxdz[0][1] = (3*mc2p_n[X_IDX] -  4*mc2p_n[15+X_IDX] +  mc2p_n[21+X_IDX] )/dzc[1]/2;
            dxdz[1][1] = (3*mc2p_n[Y_IDX] -  4*mc2p_n[15+Y_IDX] +  mc2p_n[21+Y_IDX] )/dzc[1]/2;
            dxdz[2][1] = (3*mc2p_n[Z_IDX] -  4*mc2p_n[15+Z_IDX] +  mc2p_n[21+Z_IDX] )/dzc[1]/2;
          }
          else {
            dxdz[0][1] = -(mc2p_n[15 +X_IDX] -  mc2p_n[18 +X_IDX])/2/dzc[1];
            dxdz[1][1] = -(mc2p_n[15 +Y_IDX] -  mc2p_n[18 +Y_IDX])/2/dzc[1];
            dxdz[2][1] = -(mc2p_n[15 +Z_IDX] -  mc2p_n[18 +Z_IDX])/2/dzc[1];
          }

          if((it == nrange->lower[TH_IDX]) && (up->local.lower[TH_IDX]== up->global.lower[TH_IDX])) {
            dxdz[0][2] = (-3*mc2p_n[X_IDX] + 4*mc2p_n[30+X_IDX] - mc2p_n[36+X_IDX])/dzc[2]/2;
            dxdz[1][2] = (-3*mc2p_n[Y_IDX] + 4*mc2p_n[30+Y_IDX] - mc2p_n[36+Y_IDX])/dzc[2]/2;
            dxdz[2][2] = (-3*mc2p_n[Z_IDX] + 4*mc2p_n[30+Z_IDX] - mc2p_n[36+Z_IDX])/dzc[2]/2;
          }
          else if((it == nrange->upper[TH_IDX]) && (up->local.upper[TH_IDX]== up->global.upper[TH_IDX])) {
            dxdz[0][2] = (3*mc2p_n[X_IDX] - 4*mc2p_n[27+X_IDX] + mc2p_n[33+X_IDX] )/dzc[2]/2;
            dxdz[1][2] = (3*mc2p_n[Y_IDX] - 4*mc2p_n[27+Y_IDX] + mc2p_n[33+Y_IDX] )/dzc[2]/2;
            dxdz[2][2] = (3*mc2p_n[Z_IDX] - 4*mc2p_n[27+Z_IDX] + mc2p_n[33+Z_IDX] )/dzc[2]/2;
          }
          else{
            dxdz[0][2] = -(mc2p_n[27 +X_IDX] - mc2p_n[30 +X_IDX])/2/dzc[2];
            dxdz[1][2] = -(mc2p_n[27 +Y_IDX] - mc2p_n[30 +Y_IDX])/2/dzc[2];
            dxdz[2][2] = -(mc2p_n[27 +Z_IDX] - mc2p_n[30 +Z_IDX])/2/dzc[2];
          }

          double *gFld_n= gkyl_array_fetch(gFld_nodal, gkyl_range_idx(nrange, cidx));
          gFld_n[0] = calc_metric(dxdz, 1, 1); 
          gFld_n[1] = calc_metric(dxdz, 1, 2); 
          gFld_n[2] = calc_metric(dxdz, 1, 3); 
          gFld_n[3] = calc_metric(dxdz, 2, 2); 
          gFld_n[4] = calc_metric(dxdz, 2, 3); 
          gFld_n[5] = calc_metric(dxdz, 3, 3); 

          double Jsq = gFld_n[0]*( gFld_n[3]*gFld_n[5] - gFld_n[4]*gFld_n[4] )
                      -gFld_n[1]*( gFld_n[1]*gFld_n[5] - gFld_n[4]*gFld_n[2] )
                      +gFld_n[2]*( gFld_n[1]*gFld_n[4] - gFld_n[3]*gFld_n[2] );
          double J = sqrt(Jsq);
          double e_1[3], e_2[3], e_3[3];
          e_1[0] = dxdz[0][0]; e_1[1] = dxdz[1][0]; e_1[2] = dxdz[2][0];
          e_2[0] = dxdz[0][1]; e_2[1] = dxdz[1][1]; e_2[2] = dxdz[2][1];
          e_3[0] = dxdz[0][2]; e_3[1] = dxdz[1][2]; e_3[2] = dxdz[2][2];
          calc_dual(J, e_2, e_3, dzdx[0]);
          calc_dual(J, e_3, e_1, dzdx[1]);
          calc_dual(J, e_1, e_2, dzdx[2]);

          double *dualFld_n= gkyl_array_fetch(dualFld_nodal, gkyl_range_idx(nrange, cidx));
          dualFld_n[0] = dzdx[0][0];
          dualFld_n[1] = dzdx[0][1];
          dualFld_n[2] = dzdx[0][2];
          dualFld_n[3] = dzdx[1][0];
          dualFld_n[4] = dzdx[1][1];
          dualFld_n[5] = dzdx[1][2];
          dualFld_n[6] = dzdx[2][0];
          dualFld_n[7] = dzdx[2][1];
          dualFld_n[8] = dzdx[2][2];

          double *tanvecFld_n= gkyl_array_fetch(tanvecFld_nodal, gkyl_range_idx(nrange, cidx));
          tanvecFld_n[0] = dxdz[0][0]; 
          tanvecFld_n[1] = dxdz[1][0]; 
          tanvecFld_n[2] = dxdz[2][0]; 
          tanvecFld_n[3] = dxdz[0][1]; 
          tanvecFld_n[4] = dxdz[1][1]; 
          tanvecFld_n[5] = dxdz[2][1]; 
          tanvecFld_n[6] = dxdz[0][2]; 
          tanvecFld_n[7] = dxdz[1][2]; 
          tanvecFld_n[8] = dxdz[2][2]; 

          double norm1 = sqrt(dualFld_n[0]*dualFld_n[0] + dualFld_n[1]*dualFld_n[1] + dualFld_n[2]*dualFld_n[2]);
          double norm2 = sqrt(dualFld_n[3]*dualFld_n[3] + dualFld_n[4]*dualFld_n[4] + dualFld_n[5]*dualFld_n[5]);
          double norm3 = sqrt(dualFld_n[6]*dualFld_n[6] + dualFld_n[7]*dualFld_n[7] + dualFld_n[8]*dualFld_n[8]);

          double *dualmagFld_n = gkyl_array_fetch(dualmagFld_nodal, gkyl_range_idx(nrange, cidx));
          dualmagFld_n[0] = norm1;
          dualmagFld_n[1] = norm2;
          dualmagFld_n[2] = norm3;
          
          // Set normal vectors
          double *normFld_n = gkyl_array_fetch(normFld_nodal, gkyl_range_idx(nrange, cidx));
          normFld_n[0] = dualFld_n[0]/norm1;
          normFld_n[1] = dualFld_n[1]/norm1;
          normFld_n[2] = dualFld_n[2]/norm1;

          normFld_n[3] = dualFld_n[3]/norm2;
          normFld_n[4] = dualFld_n[4]/norm2;
          normFld_n[5] = dualFld_n[5]/norm2;

          normFld_n[6] = dualFld_n[6]/norm3;
          normFld_n[7] = dualFld_n[7]/norm3;
          normFld_n[8] = dualFld_n[8]/norm3;
      }
    }
  }
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 6, gFld_nodal, gFld, false);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 9, tanvecFld_nodal, tanvecFld, false);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 9, dualFld_nodal, dualFld, false);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 3, dualmagFld_nodal, dualmagFld, false);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 9, normFld_nodal, normFld, false);
  gkyl_array_release(gFld_nodal);
  gkyl_array_release(tanvecFld_nodal);
  gkyl_array_release(dualFld_nodal);
  gkyl_array_release(dualmagFld_nodal);
  gkyl_array_release(normFld_nodal);
}

void gkyl_calc_metric_advance_interior(gkyl_calc_metric *up, struct gk_geometry *gk_geom)
{
  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates
  int cidx[3];
  
  for(int ia=gk_geom->nrange_int.lower[AL_IDX]; ia<=gk_geom->nrange_int.upper[AL_IDX]; ++ia){
    for (int ip=gk_geom->nrange_int.lower[PSI_IDX]; ip<=gk_geom->nrange_int.upper[PSI_IDX]; ++ip) {
      for (int it=gk_geom->nrange_int.lower[TH_IDX]; it<=gk_geom->nrange_int.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        const double *mc2p_n = gkyl_array_cfetch(gk_geom->geo_int.mc2p_nodal_fd, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        double dxdz[3][3]; // tan vecs at node
        double dzdx[3][3]; // duals at node

        dxdz[0][0] = -(mc2p_n[3 +X_IDX] - mc2p_n[6+X_IDX])/2/gk_geom->dzc[0];
        dxdz[1][0] = -(mc2p_n[3 +Y_IDX] - mc2p_n[6+Y_IDX])/2/gk_geom->dzc[0];
        dxdz[2][0] = -(mc2p_n[3 +Z_IDX] - mc2p_n[6+Z_IDX])/2/gk_geom->dzc[0];


        dxdz[0][1] = -(mc2p_n[15 +X_IDX] - mc2p_n[18 +X_IDX])/2/gk_geom->dzc[1];
        dxdz[1][1] = -(mc2p_n[15 +Y_IDX] - mc2p_n[18 +Y_IDX])/2/gk_geom->dzc[1];
        dxdz[2][1] = -(mc2p_n[15 +Z_IDX] - mc2p_n[18 +Z_IDX])/2/gk_geom->dzc[1];

        dxdz[0][2] = -(mc2p_n[27 +X_IDX] - mc2p_n[30 +X_IDX])/2/gk_geom->dzc[2];
        dxdz[1][2] = -(mc2p_n[27 +Y_IDX] - mc2p_n[30 +Y_IDX])/2/gk_geom->dzc[2];
        dxdz[2][2] = -(mc2p_n[27 +Z_IDX] - mc2p_n[30 +Z_IDX])/2/gk_geom->dzc[2];


        const double *bhat_n = gkyl_array_cfetch(gk_geom->geo_int.b_i_nodal_fd, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        double dbhatdz[3][3]; // tan vecs at node

        dbhatdz[0][0] = -(bhat_n[3 +X_IDX] - bhat_n[6+X_IDX])/2/gk_geom->dzc[0];
        dbhatdz[1][0] = -(bhat_n[3 +Y_IDX] - bhat_n[6+Y_IDX])/2/gk_geom->dzc[0];
        dbhatdz[2][0] = -(bhat_n[3 +Z_IDX] - bhat_n[6+Z_IDX])/2/gk_geom->dzc[0];

        dbhatdz[0][1] = -(bhat_n[15 +X_IDX] - bhat_n[18 +X_IDX])/2/gk_geom->dzc[1];
        dbhatdz[1][1] = -(bhat_n[15 +Y_IDX] - bhat_n[18 +Y_IDX])/2/gk_geom->dzc[1];
        dbhatdz[2][1] = -(bhat_n[15 +Z_IDX] - bhat_n[18 +Z_IDX])/2/gk_geom->dzc[1];

        dbhatdz[0][2] = -(bhat_n[27 +X_IDX] - bhat_n[30 +X_IDX])/2/gk_geom->dzc[2];
        dbhatdz[1][2] = -(bhat_n[27 +Y_IDX] - bhat_n[30 +Y_IDX])/2/gk_geom->dzc[2];
        dbhatdz[2][2] = -(bhat_n[27 +Z_IDX] - bhat_n[30 +Z_IDX])/2/gk_geom->dzc[2];

        double *gFld_n= gkyl_array_fetch(gk_geom->geo_int.g_ij_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        gFld_n[0] = calc_metric(dxdz, 1, 1); 
        gFld_n[1] = calc_metric(dxdz, 1, 2); 
        gFld_n[2] = calc_metric(dxdz, 1, 3); 
        gFld_n[3] = calc_metric(dxdz, 2, 2); 
        gFld_n[4] = calc_metric(dxdz, 2, 3); 
        gFld_n[5] = calc_metric(dxdz, 3, 3); 

        double Jsq = gFld_n[0]*( gFld_n[3]*gFld_n[5] - gFld_n[4]*gFld_n[4] )
                    -gFld_n[1]*( gFld_n[1]*gFld_n[5] - gFld_n[4]*gFld_n[2] )
                    +gFld_n[2]*( gFld_n[1]*gFld_n[4] - gFld_n[3]*gFld_n[2] );
        double J = sqrt(Jsq);
        double e_1[3], e_2[3], e_3[3];
        e_1[0] = dxdz[0][0]; e_1[1] = dxdz[1][0]; e_1[2] = dxdz[2][0];
        e_2[0] = dxdz[0][1]; e_2[1] = dxdz[1][1]; e_2[2] = dxdz[2][1];
        e_3[0] = dxdz[0][2]; e_3[1] = dxdz[1][2]; e_3[2] = dxdz[2][2];
        calc_dual(J, e_2, e_3, dzdx[0]);
        calc_dual(J, e_3, e_1, dzdx[1]);
        calc_dual(J, e_1, e_2, dzdx[2]);

        double *dualFld_n= gkyl_array_fetch(gk_geom->geo_int.dzdx_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        dualFld_n[0] = dzdx[0][0];
        dualFld_n[1] = dzdx[0][1];
        dualFld_n[2] = dzdx[0][2];
        dualFld_n[3] = dzdx[1][0];
        dualFld_n[4] = dzdx[1][1];
        dualFld_n[5] = dzdx[1][2];
        dualFld_n[6] = dzdx[2][0];
        dualFld_n[7] = dzdx[2][1];
        dualFld_n[8] = dzdx[2][2];

        double *tanvecFld_n= gkyl_array_fetch(gk_geom->geo_int.dxdz_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        tanvecFld_n[0] = dxdz[0][0]; 
        tanvecFld_n[1] = dxdz[1][0]; 
        tanvecFld_n[2] = dxdz[2][0]; 
        tanvecFld_n[3] = dxdz[0][1]; 
        tanvecFld_n[4] = dxdz[1][1]; 
        tanvecFld_n[5] = dxdz[2][1]; 
        tanvecFld_n[6] = dxdz[0][2]; 
        tanvecFld_n[7] = dxdz[1][2]; 
        tanvecFld_n[8] = dxdz[2][2]; 

        // Check that the coordinate system has tangent/dual vectors
        // satisfying orthonormality, and that it's right handed.
        check_orthonormality(tanvecFld_n, dualFld_n, up->exit_at_checks);
        check_right_handed(tanvecFld_n, dualFld_n, up->exit_at_checks);

        // Check if bhat and e_3 are parallel.
        double bhat_vec[3] = {bhat_n[X_IDX], bhat_n[Y_IDX], bhat_n[Z_IDX]};
        double e_3_norm[3] = {tanvecFld_n[6]/sqrt(gFld_n[5]),
                              tanvecFld_n[7]/sqrt(gFld_n[5]),
                              tanvecFld_n[8]/sqrt(gFld_n[5])};
        check_parallel(bhat_vec, e_3_norm, up->exit_at_checks);

        double norm1 = sqrt(dualFld_n[0]*dualFld_n[0] + dualFld_n[1]*dualFld_n[1] + dualFld_n[2]*dualFld_n[2]);
        double norm2 = sqrt(dualFld_n[3]*dualFld_n[3] + dualFld_n[4]*dualFld_n[4] + dualFld_n[5]*dualFld_n[5]);
        double norm3 = sqrt(dualFld_n[6]*dualFld_n[6] + dualFld_n[7]*dualFld_n[7] + dualFld_n[8]*dualFld_n[8]);

        double *dualmagFld_n = gkyl_array_fetch(gk_geom->geo_int.dualmag_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        dualmagFld_n[0] = norm1;
        dualmagFld_n[1] = norm2;
        dualmagFld_n[2] = norm3;
        
        // Set normal vectors
        double *normFld_n = gkyl_array_fetch(gk_geom->geo_int.normals_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        normFld_n[0] = dualFld_n[0]/norm1;
        normFld_n[1] = dualFld_n[1]/norm1;
        normFld_n[2] = dualFld_n[2]/norm1;

        normFld_n[3] = dualFld_n[3]/norm2;
        normFld_n[4] = dualFld_n[4]/norm2;
        normFld_n[5] = dualFld_n[5]/norm2;

        normFld_n[6] = dualFld_n[6]/norm3;
        normFld_n[7] = dualFld_n[7]/norm3;
        normFld_n[8] = dualFld_n[8]/norm3;

        double *bmag_n = gkyl_array_fetch(gk_geom->geo_int.bmag_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        // Set e^m \dot curl(bhat) 
        double *curlbhat_n = gkyl_array_fetch(gk_geom->geo_int.curlbhat_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        // I first need the derivatives of B_X,Y,Z wrt XYZ
        double dbhatdX[3][3];
        matTvec(dzdx, dbhatdz[0], dbhatdX[0]);
        matTvec(dzdx, dbhatdz[1], dbhatdX[1]);
        matTvec(dzdx, dbhatdz[2], dbhatdX[2]);
        curlbhat_n[0] = (dbhatdX[2][1] - dbhatdX[1][2]);
        curlbhat_n[1] = (dbhatdX[0][2] - dbhatdX[2][0]);
        curlbhat_n[2] = (dbhatdX[1][0] - dbhatdX[0][1]);
        double *dualcurlbhat_n = gkyl_array_fetch(gk_geom->geo_int.dualcurlbhat_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        dualcurlbhat_n[0] = dualFld_n[0]*curlbhat_n[0] + dualFld_n[1]*curlbhat_n[1] + dualFld_n[2]*curlbhat_n[2];
        dualcurlbhat_n[1] = dualFld_n[3]*curlbhat_n[0] + dualFld_n[4]*curlbhat_n[1] + dualFld_n[5]*curlbhat_n[2];
        dualcurlbhat_n[2] = dualFld_n[6]*curlbhat_n[0] + dualFld_n[7]*curlbhat_n[1] + dualFld_n[8]*curlbhat_n[2];

        // Set e^3 \dot B 
        double *B3_n = gkyl_array_fetch(gk_geom->geo_int.B3_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        B3_n[0] = bmag_n[0]/sqrt(gFld_n[5]);

        // set e^3 \dot B /|B|
        double *dualcurlbhatoverB_n = gkyl_array_fetch(gk_geom->geo_int.dualcurlbhatoverB_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        dualcurlbhatoverB_n[0] = dualcurlbhat_n[0]/bmag_n[0];
        dualcurlbhatoverB_n[1] = dualcurlbhat_n[1]/bmag_n[0];
        dualcurlbhatoverB_n[2] = dualcurlbhat_n[2]/bmag_n[0];

        // set B^3/B = 1/sqrt(g_33)
        double *rtg33inv_n = gkyl_array_fetch(gk_geom->geo_int.rtg33inv_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        rtg33inv_n[0] = 1.0/sqrt(gFld_n[5]);

        // set b_i/JB
        double* bioverJB_n = gkyl_array_fetch(gk_geom->geo_int.bioverJB_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        bioverJB_n[0] = gFld_n[2]/sqrt(gFld_n[5])/J/bmag_n[0];
        bioverJB_n[1] = gFld_n[4]/sqrt(gFld_n[5])/J/bmag_n[0];
        bioverJB_n[2] = gFld_n[5]/sqrt(gFld_n[5])/J/bmag_n[0];
      }
    }
  }

  check_axisymmetric(gk_geom->geo_int.g_ij_nodal, &gk_geom->nrange_int, up->exit_at_checks);

  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &gk_geom->local, 6, gk_geom->geo_int.g_ij_nodal, gk_geom->geo_int.g_ij, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &gk_geom->local, 9, gk_geom->geo_int.dxdz_nodal, gk_geom->geo_int.dxdz, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &gk_geom->local, 9, gk_geom->geo_int.dzdx_nodal, gk_geom->geo_int.dzdx, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &gk_geom->local, 3, gk_geom->geo_int.dualmag_nodal, gk_geom->geo_int.dualmag, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &gk_geom->local, 9, gk_geom->geo_int.normals_nodal, gk_geom->geo_int.normals, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &gk_geom->local, 3, gk_geom->geo_int.dualcurlbhat_nodal, gk_geom->geo_int.dualcurlbhat, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &gk_geom->local, 3, gk_geom->geo_int.dualcurlbhatoverB_nodal, gk_geom->geo_int.dualcurlbhatoverB, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &gk_geom->local, 1, gk_geom->geo_int.rtg33inv_nodal, gk_geom->geo_int.rtg33inv, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &gk_geom->local, 3, gk_geom->geo_int.bioverJB_nodal, gk_geom->geo_int.bioverJB, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &gk_geom->local, 1, gk_geom->geo_int.B3_nodal, gk_geom->geo_int.B3, true);
}

void gkyl_calc_metric_advance_surface(gkyl_calc_metric *up, int dir, struct gk_geometry *gk_geom)
{
  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates
  int cidx[3];
  
  for(int ia=gk_geom->nrange_surf[dir].lower[AL_IDX]; ia<=gk_geom->nrange_surf[dir].upper[AL_IDX]; ++ia){
    for (int ip=gk_geom->nrange_surf[dir].lower[PSI_IDX]; ip<=gk_geom->nrange_surf[dir].upper[PSI_IDX]; ++ip) {
      for (int it=gk_geom->nrange_surf[dir].lower[TH_IDX]; it<=gk_geom->nrange_surf[dir].upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        const double *mc2p_n = gkyl_array_cfetch(gk_geom->geo_surf[dir].mc2p_nodal_fd, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        double dxdz[3][3]; // tan vecs at node
        double dzdx[3][3]; // duals at node

        const double *bhat_n = gkyl_array_cfetch(gk_geom->geo_surf[dir].b_i_nodal_fd, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        double dbhatdz[3][3]; // tan vecs at node

        if((ip == gk_geom->nrange_surf[dir].lower[PSI_IDX]) && (up->local.lower[PSI_IDX]== up->global.lower[PSI_IDX]) && dir==0) {
          dxdz[0][0] = (-3*mc2p_n[X_IDX] + 4*mc2p_n[6+X_IDX] - mc2p_n[12+X_IDX] )/gk_geom->dzc[0]/2;
          dxdz[1][0] = (-3*mc2p_n[Y_IDX] + 4*mc2p_n[6+Y_IDX] - mc2p_n[12+Y_IDX] )/gk_geom->dzc[0]/2;
          dxdz[2][0] = (-3*mc2p_n[Z_IDX] + 4*mc2p_n[6+Z_IDX] - mc2p_n[12+Z_IDX] )/gk_geom->dzc[0]/2;

          dbhatdz[0][0] = (-3*bhat_n[X_IDX] + 4*bhat_n[6+X_IDX] - bhat_n[12+X_IDX] )/gk_geom->dzc[0]/2;
          dbhatdz[1][0] = (-3*bhat_n[Y_IDX] + 4*bhat_n[6+Y_IDX] - bhat_n[12+Y_IDX] )/gk_geom->dzc[0]/2;
          dbhatdz[2][0] = (-3*bhat_n[Z_IDX] + 4*bhat_n[6+Z_IDX] - bhat_n[12+Z_IDX] )/gk_geom->dzc[0]/2;
        }
        else if((ip == gk_geom->nrange_surf[dir].upper[PSI_IDX]) && (up->local.upper[PSI_IDX]== up->global.upper[PSI_IDX]) && dir==0) {
          dxdz[0][0] = (3*mc2p_n[X_IDX] - 4*mc2p_n[3+X_IDX] + mc2p_n[9+X_IDX] )/gk_geom->dzc[0]/2;
          dxdz[1][0] = (3*mc2p_n[Y_IDX] - 4*mc2p_n[3+Y_IDX] + mc2p_n[9+Y_IDX] )/gk_geom->dzc[0]/2;
          dxdz[2][0] = (3*mc2p_n[Z_IDX] - 4*mc2p_n[3+Z_IDX] + mc2p_n[9+Z_IDX] )/gk_geom->dzc[0]/2;

          dbhatdz[0][0] = (3*bhat_n[X_IDX] - 4*bhat_n[3+X_IDX] + bhat_n[9+X_IDX] )/gk_geom->dzc[0]/2;
          dbhatdz[1][0] = (3*bhat_n[Y_IDX] - 4*bhat_n[3+Y_IDX] + bhat_n[9+Y_IDX] )/gk_geom->dzc[0]/2;
          dbhatdz[2][0] = (3*bhat_n[Z_IDX] - 4*bhat_n[3+Z_IDX] + bhat_n[9+Z_IDX] )/gk_geom->dzc[0]/2;
        }
        else{
          dxdz[0][0] = -(mc2p_n[3 +X_IDX] -   mc2p_n[6+X_IDX])/2/gk_geom->dzc[0];
          dxdz[1][0] = -(mc2p_n[3 +Y_IDX] -   mc2p_n[6+Y_IDX])/2/gk_geom->dzc[0];
          dxdz[2][0] = -(mc2p_n[3 +Z_IDX] -   mc2p_n[6+Z_IDX])/2/gk_geom->dzc[0];

          dbhatdz[0][0] = -(bhat_n[3 +X_IDX] -   bhat_n[6+X_IDX])/2/gk_geom->dzc[0];
          dbhatdz[1][0] = -(bhat_n[3 +Y_IDX] -   bhat_n[6+Y_IDX])/2/gk_geom->dzc[0];
          dbhatdz[2][0] = -(bhat_n[3 +Z_IDX] -   bhat_n[6+Z_IDX])/2/gk_geom->dzc[0];
        }


        if((ia == gk_geom->nrange_surf[dir].lower[AL_IDX]) && (up->local.lower[AL_IDX]== up->global.lower[AL_IDX]) && dir==1 ) {
          dxdz[0][1] = (-3*mc2p_n[X_IDX] +  4*mc2p_n[18+X_IDX] -  mc2p_n[24+X_IDX])/gk_geom->dzc[1]/2;
          dxdz[1][1] = (-3*mc2p_n[Y_IDX] +  4*mc2p_n[18+Y_IDX] -  mc2p_n[24+Y_IDX])/gk_geom->dzc[1]/2;
          dxdz[2][1] = (-3*mc2p_n[Z_IDX] +  4*mc2p_n[18+Z_IDX] -  mc2p_n[24+Z_IDX])/gk_geom->dzc[1]/2;

          dbhatdz[0][1] = (-3*bhat_n[X_IDX] +  4*bhat_n[18+X_IDX] -  bhat_n[24+X_IDX])/gk_geom->dzc[1]/2;
          dbhatdz[1][1] = (-3*bhat_n[Y_IDX] +  4*bhat_n[18+Y_IDX] -  bhat_n[24+Y_IDX])/gk_geom->dzc[1]/2;
          dbhatdz[2][1] = (-3*bhat_n[Z_IDX] +  4*bhat_n[18+Z_IDX] -  bhat_n[24+Z_IDX])/gk_geom->dzc[1]/2;
        }
        else if((ia == gk_geom->nrange_surf[dir].upper[AL_IDX])  && (up->local.upper[AL_IDX]== up->global.upper[AL_IDX]) && dir==1 ) {
          dxdz[0][1] = (3*mc2p_n[X_IDX] -  4*mc2p_n[15+X_IDX] +  mc2p_n[21+X_IDX] )/gk_geom->dzc[1]/2;
          dxdz[1][1] = (3*mc2p_n[Y_IDX] -  4*mc2p_n[15+Y_IDX] +  mc2p_n[21+Y_IDX] )/gk_geom->dzc[1]/2;
          dxdz[2][1] = (3*mc2p_n[Z_IDX] -  4*mc2p_n[15+Z_IDX] +  mc2p_n[21+Z_IDX] )/gk_geom->dzc[1]/2;

          dbhatdz[0][1] = (3*bhat_n[X_IDX] -  4*bhat_n[15+X_IDX] +  bhat_n[21+X_IDX] )/gk_geom->dzc[1]/2;
          dbhatdz[1][1] = (3*bhat_n[Y_IDX] -  4*bhat_n[15+Y_IDX] +  bhat_n[21+Y_IDX] )/gk_geom->dzc[1]/2;
          dbhatdz[2][1] = (3*bhat_n[Z_IDX] -  4*bhat_n[15+Z_IDX] +  bhat_n[21+Z_IDX] )/gk_geom->dzc[1]/2;
        }
        else {
          dxdz[0][1] = -(mc2p_n[15 +X_IDX] -  mc2p_n[18 +X_IDX])/2/gk_geom->dzc[1];
          dxdz[1][1] = -(mc2p_n[15 +Y_IDX] -  mc2p_n[18 +Y_IDX])/2/gk_geom->dzc[1];
          dxdz[2][1] = -(mc2p_n[15 +Z_IDX] -  mc2p_n[18 +Z_IDX])/2/gk_geom->dzc[1];

          dbhatdz[0][1] = -(bhat_n[15 +X_IDX] -  bhat_n[18 +X_IDX])/2/gk_geom->dzc[1];
          dbhatdz[1][1] = -(bhat_n[15 +Y_IDX] -  bhat_n[18 +Y_IDX])/2/gk_geom->dzc[1];
          dbhatdz[2][1] = -(bhat_n[15 +Z_IDX] -  bhat_n[18 +Z_IDX])/2/gk_geom->dzc[1];
        }

        if((it == gk_geom->nrange_surf[dir].lower[TH_IDX]) && (up->local.lower[TH_IDX]== up->global.lower[TH_IDX]) && dir==2 ) {
          dxdz[0][2] = (-3*mc2p_n[X_IDX] + 4*mc2p_n[30+X_IDX] - mc2p_n[36+X_IDX])/gk_geom->dzc[2]/2;
          dxdz[1][2] = (-3*mc2p_n[Y_IDX] + 4*mc2p_n[30+Y_IDX] - mc2p_n[36+Y_IDX])/gk_geom->dzc[2]/2;
          dxdz[2][2] = (-3*mc2p_n[Z_IDX] + 4*mc2p_n[30+Z_IDX] - mc2p_n[36+Z_IDX])/gk_geom->dzc[2]/2;

          dbhatdz[0][2] = (-3*bhat_n[X_IDX] + 4*bhat_n[30+X_IDX] - bhat_n[36+X_IDX])/gk_geom->dzc[2]/2;
          dbhatdz[1][2] = (-3*bhat_n[Y_IDX] + 4*bhat_n[30+Y_IDX] - bhat_n[36+Y_IDX])/gk_geom->dzc[2]/2;
          dbhatdz[2][2] = (-3*bhat_n[Z_IDX] + 4*bhat_n[30+Z_IDX] - bhat_n[36+Z_IDX])/gk_geom->dzc[2]/2;
        }
        else if((it == gk_geom->nrange_surf[dir].upper[TH_IDX]) && (up->local.upper[TH_IDX]== up->global.upper[TH_IDX]) && dir==2 ) {
          dxdz[0][2] = (3*mc2p_n[X_IDX] - 4*mc2p_n[27+X_IDX] + mc2p_n[33+X_IDX] )/gk_geom->dzc[2]/2;
          dxdz[1][2] = (3*mc2p_n[Y_IDX] - 4*mc2p_n[27+Y_IDX] + mc2p_n[33+Y_IDX] )/gk_geom->dzc[2]/2;
          dxdz[2][2] = (3*mc2p_n[Z_IDX] - 4*mc2p_n[27+Z_IDX] + mc2p_n[33+Z_IDX] )/gk_geom->dzc[2]/2;

          dbhatdz[0][2] = (3*bhat_n[X_IDX] - 4*bhat_n[27+X_IDX] + bhat_n[33+X_IDX] )/gk_geom->dzc[2]/2;
          dbhatdz[1][2] = (3*bhat_n[Y_IDX] - 4*bhat_n[27+Y_IDX] + bhat_n[33+Y_IDX] )/gk_geom->dzc[2]/2;
          dbhatdz[2][2] = (3*bhat_n[Z_IDX] - 4*bhat_n[27+Z_IDX] + bhat_n[33+Z_IDX] )/gk_geom->dzc[2]/2;
        }
        else{
          dxdz[0][2] = -(mc2p_n[27 +X_IDX] - mc2p_n[30 +X_IDX])/2/gk_geom->dzc[2];
          dxdz[1][2] = -(mc2p_n[27 +Y_IDX] - mc2p_n[30 +Y_IDX])/2/gk_geom->dzc[2];
          dxdz[2][2] = -(mc2p_n[27 +Z_IDX] - mc2p_n[30 +Z_IDX])/2/gk_geom->dzc[2];

          dbhatdz[0][2] = -(bhat_n[27 +X_IDX] - bhat_n[30 +X_IDX])/2/gk_geom->dzc[2];
          dbhatdz[1][2] = -(bhat_n[27 +Y_IDX] - bhat_n[30 +Y_IDX])/2/gk_geom->dzc[2];
          dbhatdz[2][2] = -(bhat_n[27 +Z_IDX] - bhat_n[30 +Z_IDX])/2/gk_geom->dzc[2];
        }

        double *gFld_n= gkyl_array_fetch(gk_geom->geo_surf[dir].g_ij_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        gFld_n[0] = calc_metric(dxdz, 1, 1); 
        gFld_n[1] = calc_metric(dxdz, 1, 2); 
        gFld_n[2] = calc_metric(dxdz, 1, 3); 
        gFld_n[3] = calc_metric(dxdz, 2, 2); 
        gFld_n[4] = calc_metric(dxdz, 2, 3); 
        gFld_n[5] = calc_metric(dxdz, 3, 3); 

        double Jsq = gFld_n[0]*( gFld_n[3]*gFld_n[5] - gFld_n[4]*gFld_n[4] )
                    -gFld_n[1]*( gFld_n[1]*gFld_n[5] - gFld_n[4]*gFld_n[2] )
                    +gFld_n[2]*( gFld_n[1]*gFld_n[4] - gFld_n[3]*gFld_n[2] );
        double J = sqrt(Jsq);
        double *jFld_n= gkyl_array_fetch(gk_geom->geo_surf[dir].jacobgeo_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        jFld_n[0] = J;
        // Calculate cmag, bi, and jtot_inv
        double *biFld_n= gkyl_array_fetch(gk_geom->geo_surf[dir].b_i_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        biFld_n[0] = gFld_n[2]/sqrt(gFld_n[5]);
        biFld_n[1] = gFld_n[4]/sqrt(gFld_n[5]);
        biFld_n[2] = gFld_n[5]/sqrt(gFld_n[5]);

        double *bmag_n = gkyl_array_fetch(gk_geom->geo_surf[dir].bmag_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        double *cmagFld_n= gkyl_array_fetch(gk_geom->geo_surf[dir].cmag_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        cmagFld_n[0] = jFld_n[0]*bmag_n[0]/sqrt(gFld_n[5]);
        double *jtotinvFld_n= gkyl_array_fetch(gk_geom->geo_surf[dir].jacobtot_inv_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        jtotinvFld_n[0] = 1.0/(jFld_n[0]*bmag_n[0]);

        double e_1[3], e_2[3], e_3[3];
        e_1[0] = dxdz[0][0]; e_1[1] = dxdz[1][0]; e_1[2] = dxdz[2][0];
        e_2[0] = dxdz[0][1]; e_2[1] = dxdz[1][1]; e_2[2] = dxdz[2][1];
        e_3[0] = dxdz[0][2]; e_3[1] = dxdz[1][2]; e_3[2] = dxdz[2][2];
        calc_dual(J, e_2, e_3, dzdx[0]);
        calc_dual(J, e_3, e_1, dzdx[1]);
        calc_dual(J, e_1, e_2, dzdx[2]);

        double *dualFld_n= gkyl_array_fetch(gk_geom->geo_surf[dir].dzdx_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        dualFld_n[0] = dzdx[0][0];
        dualFld_n[1] = dzdx[0][1];
        dualFld_n[2] = dzdx[0][2];
        dualFld_n[3] = dzdx[1][0];
        dualFld_n[4] = dzdx[1][1];
        dualFld_n[5] = dzdx[1][2];
        dualFld_n[6] = dzdx[2][0];
        dualFld_n[7] = dzdx[2][1];
        dualFld_n[8] = dzdx[2][2];

        double *tanvecFld_n= gkyl_array_fetch(gk_geom->geo_surf[dir].dxdz_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        tanvecFld_n[0] = dxdz[0][0]; 
        tanvecFld_n[1] = dxdz[1][0]; 
        tanvecFld_n[2] = dxdz[2][0]; 
        tanvecFld_n[3] = dxdz[0][1]; 
        tanvecFld_n[4] = dxdz[1][1]; 
        tanvecFld_n[5] = dxdz[2][1]; 
        tanvecFld_n[6] = dxdz[0][2]; 
        tanvecFld_n[7] = dxdz[1][2]; 
        tanvecFld_n[8] = dxdz[2][2]; 

        // Check that the coordinate system has tangent/dual vectors
        // satisfying orthonormality, and that it's right handed.
        check_orthonormality(tanvecFld_n, dualFld_n, up->exit_at_checks);
        check_right_handed(tanvecFld_n, dualFld_n, up->exit_at_checks);

        double norm1 = sqrt(dualFld_n[0]*dualFld_n[0] + dualFld_n[1]*dualFld_n[1] + dualFld_n[2]*dualFld_n[2]);
        double norm2 = sqrt(dualFld_n[3]*dualFld_n[3] + dualFld_n[4]*dualFld_n[4] + dualFld_n[5]*dualFld_n[5]);
        double norm3 = sqrt(dualFld_n[6]*dualFld_n[6] + dualFld_n[7]*dualFld_n[7] + dualFld_n[8]*dualFld_n[8]);

        double *dualmagFld_n = gkyl_array_fetch(gk_geom->geo_surf[dir].dualmag_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        dualmagFld_n[0] = norm1;
        dualmagFld_n[1] = norm2;
        dualmagFld_n[2] = norm3;
        
        // Set normal vectors
        double *normFld_n = gkyl_array_fetch(gk_geom->geo_surf[dir].normals_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        normFld_n[0] = dualFld_n[0]/norm1;
        normFld_n[1] = dualFld_n[1]/norm1;
        normFld_n[2] = dualFld_n[2]/norm1;

        normFld_n[3] = dualFld_n[3]/norm2;
        normFld_n[4] = dualFld_n[4]/norm2;
        normFld_n[5] = dualFld_n[5]/norm2;

        normFld_n[6] = dualFld_n[6]/norm3;
        normFld_n[7] = dualFld_n[7]/norm3;
        normFld_n[8] = dualFld_n[8]/norm3;

        // Set lenr
        double *lenr_n = gkyl_array_fetch(gk_geom->geo_surf[dir].lenr_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        lenr_n[0] = J*dualmagFld_n[dir];

        // Set n^3 \dot B 
        double *B3_n = gkyl_array_fetch(gk_geom->geo_surf[dir].B3_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        B3_n[0] = bmag_n[0]/sqrt(gFld_n[5])/norm3;

        // Set n^m \dot curl(bhat) 
        double *curlbhat_n = gkyl_array_fetch(gk_geom->geo_surf[dir].curlbhat_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        // I first need the derivatives of B_X,Y,Z wrt XYZ
        double dbhatdX[3][3];
        matTvec(dzdx, dbhatdz[0], dbhatdX[0]);
        matTvec(dzdx, dbhatdz[1], dbhatdX[1]);
        matTvec(dzdx, dbhatdz[2], dbhatdX[2]);
        curlbhat_n[0] = (dbhatdX[2][1] - dbhatdX[1][2]);
        curlbhat_n[1] = (dbhatdX[0][2] - dbhatdX[2][0]);
        curlbhat_n[2] = (dbhatdX[1][0] - dbhatdX[0][1]);
        double *normcurlbhat_n = gkyl_array_fetch(gk_geom->geo_surf[dir].normcurlbhat_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        normcurlbhat_n[0] = normFld_n[3*dir+0]*curlbhat_n[0] +  normFld_n[3*dir+1]*curlbhat_n[1] + normFld_n[3*dir+2]*curlbhat_n[2];

         // set bimpactangle = arcsin(1/sqrt(g_33 * g^33))
         double *bimpactangle_n = gkyl_array_fetch(gk_geom->geo_surf[dir].bimpactangle_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
         bimpactangle_n[0] = asin(1.0/(sqrt(gFld_n[5]) * norm3));

      }
    }
  }

  check_axisymmetric(gk_geom->geo_surf[dir].g_ij_nodal, &gk_geom->nrange_surf[dir], up->exit_at_checks);
}

void gkyl_calc_metric_advance_bcart(gkyl_calc_metric *up, struct gkyl_range *nrange,
  struct gkyl_array *biFld, struct gkyl_array *dualFld, struct gkyl_array *bcartFld,
  const struct gkyl_range *update_range)
{
  struct gkyl_array* bcartFld_nodal = gkyl_array_new(GKYL_DOUBLE, 3, nrange->volume);
  struct gkyl_array* biFld_nodal = gkyl_array_new(GKYL_DOUBLE, 3, nrange->volume);
  struct gkyl_array* dualFld_nodal = gkyl_array_new(GKYL_DOUBLE, 9, nrange->volume);
  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates
  int cidx[3];

  // Fill the inputs
  gkyl_nodal_ops_m2n(up->n2m, up->cbasis, up->grid, nrange, update_range, 9, dualFld_nodal, dualFld, true);
  gkyl_nodal_ops_m2n(up->n2m, up->cbasis, up->grid, nrange, update_range, 3, biFld_nodal, biFld, true);
  
  for(int ia=nrange->lower[AL_IDX]; ia<=nrange->upper[AL_IDX]; ++ia) {
    for (int ip=nrange->lower[PSI_IDX]; ip<=nrange->upper[PSI_IDX]; ++ip) {
      for (int it=nrange->lower[TH_IDX]; it<=nrange->upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;

        double *biFld_n= gkyl_array_fetch(biFld_nodal, gkyl_range_idx(nrange, cidx));
        double *dualFld_n= gkyl_array_fetch(dualFld_nodal, gkyl_range_idx(nrange, cidx));
        double dzdx[3][3]; // duals at node
        dzdx[0][0] = dualFld_n[0];
        dzdx[0][1] = dualFld_n[1];
        dzdx[0][2] = dualFld_n[2];
        dzdx[1][0] = dualFld_n[3];
        dzdx[1][1] = dualFld_n[4];
        dzdx[1][2] = dualFld_n[5];
        dzdx[2][0] = dualFld_n[6];
        dzdx[2][1] = dualFld_n[7];
        dzdx[2][2] = dualFld_n[8];
        double *bcartFld_n= gkyl_array_fetch(bcartFld_nodal, gkyl_range_idx(nrange, cidx));
        bcartFld_n[0] = dzdx[0][0]*biFld_n[0] + dzdx[1][0]*biFld_n[1] + dzdx[2][0]*biFld_n[2];
        bcartFld_n[1] = dzdx[0][1]*biFld_n[0] + dzdx[1][1]*biFld_n[1] + dzdx[2][1]*biFld_n[2];
        bcartFld_n[2] = dzdx[0][2]*biFld_n[0] + dzdx[1][2]*biFld_n[1] + dzdx[2][2]*biFld_n[2];
      }
    }
  }
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 3, bcartFld_nodal, bcartFld, false);
  gkyl_array_release(bcartFld_nodal);
  gkyl_array_release(biFld_nodal);
  gkyl_array_release(dualFld_nodal);
}

void
gkyl_calc_metric_release(gkyl_calc_metric* up)
{
  gkyl_nodal_ops_release(up->n2m);
  gkyl_free(up);
}
