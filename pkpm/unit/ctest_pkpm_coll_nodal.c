// Tests for the configuration-nodal, velocity-modal conversions of the PKPM collision kernels: they
// must agree with the basis functions, invert each other, and turn p=1 weak products into pointwise
// products at the nodes.
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <acutest.h>
#include <gkyl_basis.h>
#include <gkyl_binop_mul_ser.h>
#include <gkyl_pkpm_coll_nodal.h>
#include <gkyl_util.h>

#define NUM_BLOCKS 8
#define NUM_VMODES 3

typedef void (*conf_map_t)(const double *in, double *GKYL_RESTRICT out);
typedef void (*mul_t)(const double *f, const double *g, double *fg);

static const conf_map_t conf_to_nodal[] = {
  pkpm_coll_conf_to_nodal_1x, pkpm_coll_conf_to_nodal_2x, pkpm_coll_conf_to_nodal_3x
};
static const conf_map_t nodal_to_conf[] = {
  pkpm_coll_nodal_to_conf_1x, pkpm_coll_nodal_to_conf_2x, pkpm_coll_nodal_to_conf_3x
};
static const conf_map_t phase_to_nodal[] = {
  pkpm_coll_phase_to_nodal_1x1v, pkpm_coll_phase_to_nodal_2x1v, pkpm_coll_phase_to_nodal_3x1v
};
static const conf_map_t nodal_to_phase_acc[] = {
  pkpm_coll_nodal_to_phase_accumulate_1x1v, pkpm_coll_nodal_to_phase_accumulate_2x1v,
  pkpm_coll_nodal_to_phase_accumulate_3x1v
};
static const mul_t conf_mul[] = {binop_mul_1d_ser_p1, binop_mul_2d_ser_p1, binop_mul_3d_ser_p1};

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

static void
node_coord(int cdim, int q, double *x)
{
  for (int d = 0; d < cdim; ++d) {
    x[d] = ((q >> d) & 1) ? 1.0 / sqrt(3.0) : -1.0 / sqrt(3.0);
  }
}

static void
test_conf(int cdim)
{
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, 1);
  int nc = cbasis.num_basis;
  double c[8], n[8], back[8], x[3];
  for (int t = 0; t < 50; ++t) {
    for (int i = 0; i < nc; ++i) {
      c[i] = rnd();
    }
    conf_to_nodal[cdim - 1](c, n);
    for (int q = 0; q < nc; ++q) {
      node_coord(cdim, q, x);
      TEST_CHECK(close_to(n[q], cbasis.eval_expand(x, c), 1e-14));
    }
    nodal_to_conf[cdim - 1](n, back);
    for (int i = 0; i < nc; ++i) {
      TEST_CHECK(close_to(back[i], c[i], 1e-14));
    }
  }
}

// The weak product of two p=1 expansions takes the pointwise product at the nodes.
static void
test_product(int cdim)
{
  int nc = 1 << cdim;
  double a[8], b[8], ab[8], na[8], nb[8], nab[8];
  for (int t = 0; t < 50; ++t) {
    for (int i = 0; i < nc; ++i) {
      a[i] = rnd();
      b[i] = rnd();
    }
    conf_mul[cdim - 1](a, b, ab);
    conf_to_nodal[cdim - 1](a, na);
    conf_to_nodal[cdim - 1](b, nb);
    conf_to_nodal[cdim - 1](ab, nab);
    for (int q = 0; q < nc; ++q) {
      TEST_CHECK(close_to(nab[q], na[q] * nb[q], 1e-13));
    }
  }
}

static void
test_phase(int cdim)
{
  struct gkyl_basis pbasis, vbasis;
  gkyl_cart_modal_hybrid(&pbasis, cdim, 1);
  gkyl_cart_modal_serendip(&vbasis, 1, 2);
  int np = pbasis.num_basis, nc = 1 << cdim;
  TEST_CHECK(np == nc * NUM_VMODES);

  static const double gv[3] = {-0.7745966692414834, 0.0, 0.7745966692414834};
  static const double wv[3] = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};

  double f[NUM_BLOCKS * 24], g[NUM_BLOCKS * 8 * NUM_VMODES], back[NUM_BLOCKS * 24];
  for (int t = 0; t < 20; ++t) {
    for (int i = 0; i < NUM_BLOCKS * np; ++i) {
      f[i] = rnd();
    }
    phase_to_nodal[cdim - 1](f, g);

    for (int b = 0; b < NUM_BLOCKS; ++b) {
      for (int q = 0; q < nc; ++q) {
        double eta[4], psi_expand[NUM_VMODES] = {0.0};
        node_coord(cdim, q, eta);
        for (int m = 0; m < 3; ++m) {
          double pv[NUM_VMODES];
          eta[cdim] = gv[m];
          vbasis.eval(&gv[m], pv);
          double value = pbasis.eval_expand(eta, &f[b * np]);
          for (int j = 0; j < NUM_VMODES; ++j) {
            psi_expand[j] += wv[m] * value * pv[j];
          }
        }
        for (int j = 0; j < NUM_VMODES; ++j) {
          TEST_CHECK(close_to(g[(q * NUM_BLOCKS + b) * NUM_VMODES + j], psi_expand[j], 1e-13));
        }
      }
    }

    memset(back, 0, sizeof(back));
    nodal_to_phase_acc[cdim - 1](g, back);
    for (int i = 0; i < NUM_BLOCKS * np; ++i) {
      TEST_CHECK(close_to(back[i], f[i], 1e-13));
    }
  }
}

void
test_conf_1x()
{
  test_conf(1);
}

void
test_conf_2x()
{
  test_conf(2);
}

void
test_conf_3x()
{
  test_conf(3);
}

void
test_product_1x()
{
  test_product(1);
}

void
test_product_2x()
{
  test_product(2);
}

void
test_product_3x()
{
  test_product(3);
}

void
test_phase_1x()
{
  test_phase(1);
}

void
test_phase_2x()
{
  test_phase(2);
}

void
test_phase_3x()
{
  test_phase(3);
}

TEST_LIST = {{"conf_1x", test_conf_1x},       {"conf_2x", test_conf_2x},
             {"conf_3x", test_conf_3x},       {"product_1x", test_product_1x},
             {"product_2x", test_product_2x}, {"product_3x", test_product_3x},
             {"phase_1x", test_phase_1x},     {"phase_2x", test_phase_2x},
             {"phase_3x", test_phase_3x},     {NULL, NULL}};
