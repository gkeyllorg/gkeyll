// The conversions between modal expansions and the configuration-nodal, velocity-modal form:
// they must agree with the basis functions, invert each other, act block by block, and turn
// p=1 weak products into pointwise products at the nodes.
#include <math.h>
#include <stdlib.h>

#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_conf_nodal.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

#define NUM_BLOCKS 3 // expansions per cell
#define NUM_CELLS 4 // cells per direction
#define NUM_VMODES 3 // velocity modes of the 1v hybrid basis

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
fill_random(struct gkyl_array *arr)
{
  double *arr_d = arr->data;
  for (size_t i = 0; i < arr->size * arr->ncomp; ++i) {
    arr_d[i] = rnd();
  }
}

static void
test_conf(int cdim)
{
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, 1);
  int nc = cbasis.num_basis;
  int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  for (int d = 0; d < cdim; ++d) {
    lower[d] = 1;
    upper[d] = NUM_CELLS;
  }
  struct gkyl_range crange;
  gkyl_range_init(&crange, cdim, lower, upper);

  struct gkyl_array *fmodal = gkyl_array_new(GKYL_DOUBLE, NUM_BLOCKS * nc, crange.volume);
  struct gkyl_array *gmodal = gkyl_array_new(GKYL_DOUBLE, NUM_BLOCKS * nc, crange.volume);
  struct gkyl_array *fg = gkyl_array_new(GKYL_DOUBLE, NUM_BLOCKS * nc, crange.volume);
  struct gkyl_array *fnodal = gkyl_array_new(GKYL_DOUBLE, NUM_BLOCKS * nc, crange.volume);
  struct gkyl_array *gnodal = gkyl_array_new(GKYL_DOUBLE, NUM_BLOCKS * nc, crange.volume);
  struct gkyl_array *fgnodal = gkyl_array_new(GKYL_DOUBLE, NUM_BLOCKS * nc, crange.volume);
  struct gkyl_array *back = gkyl_array_new(GKYL_DOUBLE, NUM_BLOCKS * nc, crange.volume);
  fill_random(fmodal);
  fill_random(gmodal);
  for (int b = 0; b < NUM_BLOCKS; ++b) {
    gkyl_dg_mul_op_range(&cbasis, b, fg, b, fmodal, b, gmodal, &crange);
  }

  struct gkyl_dg_conf_nodal *up = gkyl_dg_conf_nodal_new(&cbasis, 0, false);
  gkyl_dg_conf_nodal_conf_to_nodal(up, &crange, fmodal, fnodal);
  gkyl_dg_conf_nodal_conf_to_nodal(up, &crange, gmodal, gnodal);
  gkyl_dg_conf_nodal_conf_to_nodal(up, &crange, fg, fgnodal);
  gkyl_dg_conf_nodal_conf_to_modal(up, &crange, fnodal, back);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &crange);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&crange, iter.idx);
    const double *fmodal_d = gkyl_array_cfetch(fmodal, loc);
    const double *fnodal_d = gkyl_array_cfetch(fnodal, loc);
    const double *gnodal_d = gkyl_array_cfetch(gnodal, loc);
    const double *fgnodal_d = gkyl_array_cfetch(fgnodal, loc);
    const double *back_d = gkyl_array_cfetch(back, loc);
    for (int b = 0; b < NUM_BLOCKS; ++b) {
      for (int q = 0; q < nc; ++q) {
        double x[3];
        node_coord(cdim, q, x);
        // node values of the expansion, and the weak product as the product of node values
        TEST_CHECK(close_to(fnodal_d[b * nc + q], cbasis.eval_expand(x, &fmodal_d[b * nc]), 1e-14));
        TEST_CHECK(
          close_to(fgnodal_d[b * nc + q], fnodal_d[b * nc + q] * gnodal_d[b * nc + q], 1e-13)
        );
      }
      for (int i = 0; i < nc; ++i) {
        TEST_CHECK(close_to(back_d[b * nc + i], fmodal_d[b * nc + i], 1e-14));
      }
    }
  }

  gkyl_dg_conf_nodal_release(up);
  gkyl_array_release(fmodal);
  gkyl_array_release(gmodal);
  gkyl_array_release(fg);
  gkyl_array_release(fnodal);
  gkyl_array_release(gnodal);
  gkyl_array_release(fgnodal);
  gkyl_array_release(back);
}

static void
test_phase(int cdim)
{
  struct gkyl_basis cbasis, pbasis, vbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, 1);
  gkyl_cart_modal_hybrid(&pbasis, cdim, 1);
  gkyl_cart_modal_serendip(&vbasis, 1, 2);
  int np = pbasis.num_basis, nc = cbasis.num_basis;
  TEST_CHECK(np == nc * NUM_VMODES);

  int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  for (int d = 0; d < cdim + 1; ++d) {
    lower[d] = 1;
    upper[d] = NUM_CELLS;
  }
  struct gkyl_range prange;
  gkyl_range_init(&prange, cdim + 1, lower, upper);

  struct gkyl_array *fmodal = gkyl_array_new(GKYL_DOUBLE, NUM_BLOCKS * np, prange.volume);
  struct gkyl_array *fnodal = gkyl_array_new(GKYL_DOUBLE, NUM_BLOCKS * np, prange.volume);
  struct gkyl_array *back = gkyl_array_new(GKYL_DOUBLE, NUM_BLOCKS * np, prange.volume);
  fill_random(fmodal);
  gkyl_array_clear(back, 0.0);

  struct gkyl_dg_conf_nodal *up = gkyl_dg_conf_nodal_new(&cbasis, &pbasis, false);
  gkyl_dg_conf_nodal_phase_to_nodal(up, &prange, fmodal, fnodal);
  gkyl_dg_conf_nodal_phase_to_modal_accumulate(up, &prange, fnodal, back);

  // 3-point Gauss-Legendre quadrature in the velocity, exact for the products of p=2 functions.
  static const double gv[3] = {-0.7745966692414834, 0.0, 0.7745966692414834};
  static const double wv[3] = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &prange);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&prange, iter.idx);
    const double *fmodal_d = gkyl_array_cfetch(fmodal, loc);
    const double *fnodal_d = gkyl_array_cfetch(fnodal, loc);
    const double *back_d = gkyl_array_cfetch(back, loc);
    for (int b = 0; b < NUM_BLOCKS; ++b) {
      for (int q = 0; q < nc; ++q) {
        // the velocity modes of the expansion at node q
        double eta[4], vel_modes[NUM_VMODES] = {0.0};
        node_coord(cdim, q, eta);
        for (int m = 0; m < 3; ++m) {
          double vel_basis[NUM_VMODES];
          eta[cdim] = gv[m];
          vbasis.eval(&gv[m], vel_basis);
          double value = pbasis.eval_expand(eta, &fmodal_d[b * np]);
          for (int j = 0; j < NUM_VMODES; ++j) {
            vel_modes[j] += wv[m] * value * vel_basis[j];
          }
        }
        for (int j = 0; j < NUM_VMODES; ++j) {
          TEST_CHECK(close_to(fnodal_d[b * np + q * NUM_VMODES + j], vel_modes[j], 1e-13));
        }
      }
      for (int i = 0; i < np; ++i) {
        TEST_CHECK(close_to(back_d[b * np + i], fmodal_d[b * np + i], 1e-13));
      }
    }
  }

  gkyl_dg_conf_nodal_release(up);
  gkyl_array_release(fmodal);
  gkyl_array_release(fnodal);
  gkyl_array_release(back);
}

static void
test_conf_1x(void)
{
  test_conf(1);
}

static void
test_conf_2x(void)
{
  test_conf(2);
}

static void
test_conf_3x(void)
{
  test_conf(3);
}

static void
test_phase_1x1v(void)
{
  test_phase(1);
}

static void
test_phase_2x1v(void)
{
  test_phase(2);
}

static void
test_phase_3x1v(void)
{
  test_phase(3);
}

TEST_LIST = {
  {"conf_1x", test_conf_1x},
  {"conf_2x", test_conf_2x},
  {"conf_3x", test_conf_3x},
  {"phase_1x1v", test_phase_1x1v},
  {"phase_2x1v", test_phase_2x1v},
  {"phase_3x1v", test_phase_3x1v},
  {NULL, NULL}
};
