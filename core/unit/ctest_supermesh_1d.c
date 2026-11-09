#include <acutest.h>

#include <math.h>
#include <stdlib.h>

#include <gkyl_supermesh_1d.h>

#define NMAX 64

static void
uniform_edges(int num_cells, double lower, double dx, double *edges)
{
  for (int i = 0; i <= num_cells; ++i) {
    edges[i] = lower + i * dx;
  }
}

static void
test_identical(void)
{
  double edges[NMAX + 1];
  struct gkyl_supermesh_1d_interval iv[3 * NMAX];
  uniform_edges(8, -1.0, 0.25, edges);
  int n = gkyl_supermesh_1d_merge(9, edges, 1.0, 0.0, 9, edges, 1.0, 0.0, 0, 0, 1e-13, iv);
  TEST_CHECK(n == 8);
  for (int i = 0; i < n; ++i) {
    TEST_CHECK(iv[i].k_l == i);
    TEST_CHECK(iv[i].k_r == i);
    TEST_CHECK(gkyl_compare(iv[i].lo, edges[i], 1e-14));
    TEST_CHECK(gkyl_compare(iv[i].up, edges[i + 1], 1e-14));
  }
}

static void
test_shifted(void)
{
  double edges[NMAX + 1];
  struct gkyl_supermesh_1d_interval iv[3 * NMAX];
  uniform_edges(4, 0.0, 1.0, edges);
  // x_r = x + 0.25, so the right grid covers x in [-0.25, 3.75]
  int n = gkyl_supermesh_1d_merge(5, edges, 1.0, 0.0, 5, edges, 1.0, 0.25, 0, 0, 1e-13, iv);
  TEST_CHECK(n == 9);
  TEST_CHECK(gkyl_compare(iv[0].lo, -0.25, 1e-14));
  TEST_CHECK(iv[0].k_l == -1);
  TEST_CHECK(iv[0].k_r == 0);
  TEST_CHECK(gkyl_compare(iv[1].lo, 0.0, 1e-14));
  TEST_CHECK(gkyl_compare(iv[1].up, 0.75, 1e-14));
  TEST_CHECK(iv[1].k_l == 0);
  TEST_CHECK(iv[1].k_r == 0);
  TEST_CHECK(iv[2].k_l == 0);
  TEST_CHECK(iv[2].k_r == 1);
  TEST_CHECK(iv[8].k_l == 3);
  TEST_CHECK(iv[8].k_r == -1);
  TEST_CHECK(gkyl_compare(iv[8].up, 4.0, 1e-14));
}

static void
test_dilated(void)
{
  double edges[NMAX + 1];
  struct gkyl_supermesh_1d_interval iv[3 * NMAX];
  uniform_edges(2, -1.0, 1.0, edges);
  // x_r = 2 x: the right grid covers x in [-0.5, 0.5]
  int n = gkyl_supermesh_1d_merge(3, edges, 1.0, 0.0, 3, edges, 2.0, 0.0, 0, 0, 1e-13, iv);
  TEST_CHECK(n == 4);
  TEST_CHECK(iv[0].k_l == 0 && iv[0].k_r == -1);
  TEST_CHECK(gkyl_compare(iv[0].up, -0.5, 1e-14));
  TEST_CHECK(iv[1].k_l == 0 && iv[1].k_r == 0);
  TEST_CHECK(iv[2].k_l == 1 && iv[2].k_r == 1);
  TEST_CHECK(iv[3].k_l == 1 && iv[3].k_r == -1);
}

static void
test_disjoint(void)
{
  double edges[NMAX + 1];
  struct gkyl_supermesh_1d_interval iv[3 * NMAX];
  uniform_edges(2, 0.0, 0.5, edges);
  // the right grid covers [2, 3]: the gap (1, 2) between the two is dropped
  int n = gkyl_supermesh_1d_merge(3, edges, 1.0, 0.0, 3, edges, 1.0, -2.0, 0, 0, 1e-13, iv);
  TEST_CHECK(n == 4);
  TEST_CHECK(iv[1].k_l == 1 && iv[1].k_r == -1);
  TEST_CHECK(gkyl_compare(iv[1].up, 1.0, 1e-14));
  TEST_CHECK(iv[2].k_l == -1 && iv[2].k_r == 0);
  TEST_CHECK(gkyl_compare(iv[2].lo, 2.0, 1e-14));
  // the same with the grids swapped
  n = gkyl_supermesh_1d_merge(3, edges, 1.0, -2.0, 3, edges, 1.0, 0.0, 0, 0, 1e-13, iv);
  TEST_CHECK(n == 4);
  TEST_CHECK(iv[0].k_l == -1 && iv[0].k_r == 0);
  TEST_CHECK(iv[3].k_l == 1 && iv[3].k_r == -1);
}

static void
test_kinks(void)
{
  double edges[NMAX + 1];
  struct gkyl_supermesh_1d_interval iv[3 * NMAX];
  uniform_edges(2, 0.0, 1.0, edges);
  // outside (-5, 9), on an edge (1) and inside a cell (0.5)
  double kinks[4] = {-5.0, 0.5, 1.0, 9.0};
  int n = gkyl_supermesh_1d_merge(3, edges, 1.0, 0.0, 3, edges, 1.0, 0.0, 4, kinks, 1e-13, iv);
  TEST_CHECK(n == 3);
  TEST_CHECK(gkyl_compare(iv[0].up, 0.5, 1e-14));
  TEST_CHECK(iv[0].k_l == 0 && iv[1].k_l == 0 && iv[2].k_l == 1);
  TEST_CHECK(gkyl_compare(iv[1].up, 1.0, 1e-14));
}

// The measure of the union of the two mapped grids (each a single connected interval).
static double
union_measure(double lo_l, double up_l, double lo_r, double up_r)
{
  double overlap = fmin(up_l, up_r) - fmax(lo_l, lo_r);
  return (up_l - lo_l) + (up_r - lo_r) - fmax(overlap, 0.0);
}

static void
test_random(void)
{
  double edges_l[NMAX + 1], edges_r[NMAX + 1];
  struct gkyl_supermesh_1d_interval iv[3 * NMAX];
  srand(4242);
  for (int trial = 0; trial < 2000; ++trial) {
    int nl = 1 + rand() % NMAX, nr = 1 + rand() % NMAX, num_kinks = rand() % 4;
    uniform_edges(nl, -1.0 - 3.0 * rand() / RAND_MAX, 0.05 + 0.2 * rand() / RAND_MAX, edges_l);
    uniform_edges(nr, -1.0 - 3.0 * rand() / RAND_MAX, 0.05 + 0.2 * rand() / RAND_MAX, edges_r);
    double r_l = exp(4.0 * (double)rand() / RAND_MAX - 2.0),
           r_r = exp(4.0 * (double)rand() / RAND_MAX - 2.0);
    double c_l = 6.0 * rand() / RAND_MAX - 3.0, c_r = 6.0 * rand() / RAND_MAX - 3.0;
    double kinks[3];
    for (int k = 0; k < num_kinks; ++k) {
      kinks[k] = -4.0 + 8.0 * (k + (double)rand() / RAND_MAX) / num_kinks;
    }
    int n = gkyl_supermesh_1d_merge(
      nl + 1, edges_l, r_l, c_l, nr + 1, edges_r, r_r, c_r, num_kinks, kinks, 1e-13, iv
    );
    TEST_CHECK(n <= gkyl_supermesh_1d_max_intervals(nl + 1, nr + 1, num_kinks));
    double len = 0.0;
    for (int i = 0; i < n; ++i) {
      TEST_CHECK(iv[i].up > iv[i].lo);
      TEST_CHECK(iv[i].k_l >= 0 || iv[i].k_r >= 0);
      if (i > 0) {
        TEST_CHECK(iv[i].lo >= iv[i - 1].up - 1e-12);
      }
      len += iv[i].up - iv[i].lo;
      // each end of the interval lies in the cell the interval claims, in each grid
      if (iv[i].k_l >= 0) {
        TEST_CHECK(r_l * iv[i].lo + c_l >= edges_l[iv[i].k_l] - 1e-10);
        TEST_CHECK(r_l * iv[i].up + c_l <= edges_l[iv[i].k_l + 1] + 1e-10);
      }
      if (iv[i].k_r >= 0) {
        TEST_CHECK(r_r * iv[i].lo + c_r >= edges_r[iv[i].k_r] - 1e-10);
        TEST_CHECK(r_r * iv[i].up + c_r <= edges_r[iv[i].k_r + 1] + 1e-10);
      }
    }
    double measure = union_measure(
      (edges_l[0] - c_l) / r_l, (edges_l[nl] - c_l) / r_l, (edges_r[0] - c_r) / r_r,
      (edges_r[nr] - c_r) / r_r
    );
    TEST_CHECK(fabs(len - measure) < 1e-9);
  }
}

TEST_LIST = {
  {"test_identical", test_identical},
  {"test_shifted", test_shifted},
  {"test_dilated", test_dilated},
  {"test_disjoint", test_disjoint},
  {"test_kinks", test_kinks},
  {"test_random", test_random},
  {NULL, NULL}
};
