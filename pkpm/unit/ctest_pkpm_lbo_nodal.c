// The modal PKPM LBO kernels against a per-node 1V operator: at p=1 in configuration space every
// weak product is a pointwise product at the Gauss-Legendre nodes, so the kernel of any cdim,
// converted to the nodes, must equal the 1x1v kernel applied at each node to x-independent data
// holding that node's values, on all 8 blocks. The batched-LU primitive-moment solve must equal
// the 2x2 solve at each node.
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <acutest.h>
#include <gkyl_lbo_pkpm_kernels.h>
#include <gkyl_mat.h>
#include <gkyl_pkpm_coll_nodal.h>
#include <gkyl_prim_lbo_pkpm_kernels.h>
#include <gkyl_util.h>

#define NB 8
#define NV 3

typedef double (*vol_t)(const double *, const double *, const double *, const double *, const double *, double *);
typedef double (*surf_t)(const double *, const double *, const double *, const double *, const double *, const double *, const double *, double *);
typedef double (*bdry_t)(const double *, const double *, const double *, const double *, int, const double *, const double *, double *);
typedef void (*map_t)(const double *in, double *GKYL_RESTRICT out);
typedef void (*prim_t)(struct gkyl_mat *, struct gkyl_mat *, const double *, const double *, const double *);

struct lbo_kernels {
  vol_t drag_vol, diff_vol;
  surf_t drag_surf, diff_surf;
  bdry_t drag_bdry, diff_bdry;
  map_t conf_to_nodal, phase_to_nodal;
  prim_t prim;
};

#define LBO_KERNELS(c)                              \
  {lbo_pkpm_drag_vol_##c##x1v_ser_p1,               \
   lbo_pkpm_diff_vol_##c##x1v_ser_p1,               \
   lbo_pkpm_drag_surfvpar_##c##x1v_ser_p1,          \
   lbo_pkpm_diff_surfvpar_##c##x1v_ser_p1,          \
   lbo_pkpm_drag_boundary_surfvpar_##c##x1v_ser_p1, \
   lbo_pkpm_diff_boundary_surfvpar_##c##x1v_ser_p1, \
   pkpm_coll_conf_to_nodal_##c##x,                  \
   pkpm_coll_phase_to_nodal_##c##x1v,               \
   pkpm_self_prim_moments_##c##x1v_ser_p1}

static const struct lbo_kernels kern[3] = {LBO_KERNELS(1), LBO_KERNELS(2), LBO_KERNELS(3)};

// Agreement relative to the larger of 1 and the magnitudes: the values can be near zero.
static bool
close_to(double a, double b, double tol)
{
  return fabs(a - b) <= tol * fmax(1.0, fmax(fabs(a), fabs(b)));
}

static double
rnd(void)
{
  return 2.0 * rand() / (double)RAND_MAX - 1.0;
}

// The 1x1v modal cell whose x-independent data are the 8 blocks of node values gq[24].
static void
node_to_1x(const double *gq, double *f1)
{
  double g1[2 * NB * NV];
  memcpy(&g1[0], gq, sizeof(double[NB * NV]));
  memcpy(&g1[NB * NV], gq, sizeof(double[NB * NV]));
  memset(f1, 0, sizeof(double[NB * 6]));
  pkpm_coll_nodal_to_phase_accumulate_1x1v(g1, f1);
}

enum kernel_kind { DRAG_VOL, DIFF_VOL, DRAG_SURF, DIFF_SURF, DRAG_BDRY, DIFF_BDRY, NUM_KINDS };

static void
apply(
  const struct lbo_kernels *k, enum kernel_kind kind, int edge, const double *w, const double *dxv,
  const double *nu, const double *nuP, const double *fl, const double *fc, const double *fr,
  double *out
)
{
  switch (kind) {
    case DRAG_VOL:
      k->drag_vol(w, dxv, nu, nuP, fc, out);
      break;
    case DIFF_VOL:
      k->diff_vol(w, dxv, nu, nuP, fc, out);
      break;
    case DRAG_SURF:
      k->drag_surf(w, dxv, nu, nuP, fl, fc, fr, out);
      break;
    case DIFF_SURF:
      k->diff_surf(w, dxv, nu, nuP, fl, fc, fr, out);
      break;
    case DRAG_BDRY:
      k->drag_bdry(w, dxv, nu, nuP, edge, fc, fr, out);
      break;
    case DIFF_BDRY:
      k->diff_bdry(w, dxv, nu, nuP, edge, fc, fr, out);
      break;
    default:
      break;
  }
}

static void
test_kernels(int cdim)
{
  const struct lbo_kernels *k = &kern[cdim - 1], *k1 = &kern[0];
  int nc = 1 << cdim, np = NV * nc;
  double nu[8], nuP[16], nu_n[8], nuU_n[8], nuT_n[8];
  double fl[NB * 24], fc[NB * 24], fr[NB * 24], out[NB * 24];
  double gl[NB * 8 * NV], gc[NB * 8 * NV], gr[NB * 8 * NV], gout[NB * 8 * NV];
  double f1l[NB * 6], f1c[NB * 6], f1r[NB * 6], o1[NB * 6], g1out[2 * NB * NV];
  double w[4], dxv[4];

  for (int t = 0; t < 100; ++t) {
    double scale = pow(sqrt(2.0), cdim);
    for (int i = 0; i < nc; ++i) {
      nu[i] = (i == 0 ? 2.0 * scale : 0.0) + 0.5 * rnd();
      nuP[i] = rnd();
      nuP[nc + i] = (i == 0 ? 1.5 * scale : 0.0) + 0.5 * rnd();
    }
    for (int i = 0; i < NB * np; ++i) {
      fl[i] = rnd();
      fc[i] = rnd();
      fr[i] = rnd();
    }
    for (int d = 0; d < cdim; ++d) {
      w[d] = rnd();
      dxv[d] = 0.5 + fabs(rnd());
    }
    w[cdim] = 3.0 * rnd();
    dxv[cdim] = 0.2 + fabs(rnd());
    double w1[2] = {0.0, w[cdim]}, d1[2] = {1.0, dxv[cdim]};

    k->conf_to_nodal(nu, nu_n);
    k->conf_to_nodal(nuP, nuU_n);
    k->conf_to_nodal(&nuP[nc], nuT_n);
    k->phase_to_nodal(fl, gl);
    k->phase_to_nodal(fc, gc);
    k->phase_to_nodal(fr, gr);

    for (int kind = 0; kind < NUM_KINDS; ++kind) {
      int num_edge = (kind == DRAG_BDRY || kind == DIFF_BDRY) ? 2 : 1;
      for (int e = 0; e < num_edge; ++e) {
        int edge = 2 * e - 1;
        memset(out, 0, sizeof(out));
        apply(k, kind, edge, w, dxv, nu, nuP, fl, fc, fr, out);
        k->phase_to_nodal(out, gout);

        for (int q = 0; q < nc; ++q) {
          double sq = sqrt(2.0);
          double nu1[2] = {sq * nu_n[q], 0.0};
          double P1[4] = {sq * nuU_n[q], 0.0, sq * nuT_n[q], 0.0};
          node_to_1x(&gl[q * NB * NV], f1l);
          node_to_1x(&gc[q * NB * NV], f1c);
          node_to_1x(&gr[q * NB * NV], f1r);
          memset(o1, 0, sizeof(o1));
          apply(k1, kind, edge, w1, d1, nu1, P1, f1l, f1c, f1r, o1);
          k1->phase_to_nodal(o1, g1out);
          for (int i = 0; i < NB * NV; ++i) {
            // the 1x1v output is x-independent, and equals the node values of the full kernel
            TEST_CHECK(close_to(g1out[i], g1out[NB * NV + i], 1e-12));
            TEST_CHECK(close_to(g1out[i], gout[q * NB * NV + i], 1e-11));
          }
        }
      }
    }
  }
}

// Solve the dense system M x = b in place by Gaussian elimination with pivoting.
static void
solve_dense(int n, double M[16][17], double *x)
{
  for (int p = 0; p < n; ++p) {
    int pv = p;
    for (int r = p + 1; r < n; ++r) {
      if (fabs(M[r][p]) > fabs(M[pv][p])) {
        pv = r;
      }
    }
    for (int j = 0; j <= n; ++j) {
      double tmp = M[p][j];
      M[p][j] = M[pv][j];
      M[pv][j] = tmp;
    }
    for (int r = 0; r < n; ++r) {
      if (r != p) {
        double fac = M[r][p] / M[p][p];
        for (int j = p; j <= n; ++j) {
          M[r][j] -= fac * M[p][j];
        }
      }
    }
  }
  for (int i = 0; i < n; ++i) {
    x[i] = M[i][n] / M[i][i];
  }
}

static void
test_prim(int cdim)
{
  const struct lbo_kernels *k = &kern[cdim - 1];
  int nc = 1 << cdim, n = 2 * nc;
  double scale = pow(sqrt(2.0), cdim), nu[8] = {1.0};
  struct gkyl_mat *A = gkyl_mat_new(n, n, 0.0), *b = gkyl_mat_new(n, 1, 0.0);

  for (int t = 0; t < 100; ++t) {
    double moms[9 * 8] = {0.0}, bc[2 * 8] = {0.0}, M[16][17], sol[16];
    for (int i = 0; i < nc; ++i) {
      for (int c = 0; c < 3; ++c) {
        moms[c * nc + i] = (i == 0 ? 2.0 * scale : 0.0) + 0.3 * rnd() / nc;
      }
      moms[3 * nc + i] = 0.3 * rnd();
      bc[i] = 0.01 * rnd();
      bc[nc + i] = 0.01 * rnd();
    }
    k->prim(A, b, moms, bc, nu);
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        M[i][j] = gkyl_mat_get(A, i, j);
      }
      M[i][n] = gkyl_mat_get(b, i, 0);
    }
    solve_dense(n, M, sol);

    double m0[8], m2[8], mG[8], m1[8], cM[8], cE[8], u[8], vt[8];
    k->conf_to_nodal(&moms[0], m0);
    k->conf_to_nodal(&moms[nc], m2);
    k->conf_to_nodal(&moms[2 * nc], mG);
    k->conf_to_nodal(&moms[3 * nc], m1);
    k->conf_to_nodal(&bc[0], cM);
    k->conf_to_nodal(&bc[nc], cE);
    k->conf_to_nodal(&sol[0], u);
    k->conf_to_nodal(&sol[nc], vt);
    for (int q = 0; q < nc; ++q) {
      double a12 = -cM[q], a22 = 3.0 * m0[q] - cE[q], r2 = 2.0 * mG[q] + m2[q];
      double det = m0[q] * a22 - a12 * m1[q];
      TEST_CHECK(close_to(u[q], (m1[q] * a22 - a12 * r2) / det, 1e-13));
      TEST_CHECK(close_to(vt[q], (m0[q] * r2 - m1[q] * m1[q]) / det, 1e-13));
    }
  }
  gkyl_mat_release(A);
  gkyl_mat_release(b);
}

void
test_kernels_1x()
{
  test_kernels(1);
}

void
test_kernels_2x()
{
  test_kernels(2);
}

void
test_kernels_3x()
{
  test_kernels(3);
}

void
test_prim_1x()
{
  test_prim(1);
}

void
test_prim_2x()
{
  test_prim(2);
}

void
test_prim_3x()
{
  test_prim(3);
}

TEST_LIST = {
  {"kernels_1x", test_kernels_1x},
  {"kernels_2x", test_kernels_2x},
  {"kernels_3x", test_kernels_3x},
  {"prim_1x", test_prim_1x},
  {"prim_2x", test_prim_2x},
  {"prim_3x", test_prim_3x},
  {NULL, NULL}
};
