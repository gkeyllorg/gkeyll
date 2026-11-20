// Tests for the 1V Maxwellian of the PKPM collision operators: projection, moments and correction.
#include <math.h>

#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_pkpm_lte.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

struct lte_setup {
  struct gkyl_rect_grid grid;
  struct gkyl_range range;
  struct gkyl_array *f;
  gkyl_pkpm_lte *up;
};

static struct lte_setup
lte_setup_new(double w_max, int num_cells)
{
  struct lte_setup s;
  double lower[] = {-w_max}, upper[] = {w_max};
  int cells[] = {num_cells}, lo[] = {1}, up[] = {num_cells};
  gkyl_rect_grid_init(&s.grid, 1, lower, upper, cells);
  gkyl_range_init(&s.range, 1, lo, up);
  s.f = gkyl_array_new(GKYL_DOUBLE, 3, s.range.volume);
  s.up = gkyl_pkpm_lte_inew(&(struct gkyl_pkpm_lte_inp){
    .vel_grid = &s.grid,
    .vel_range = &s.range,
    .num_quad = 4,
    .eps = 1.0e-14,
    .max_iter = 100,
  });
  return s;
}

static void
lte_setup_release(struct lte_setup *s)
{
  gkyl_array_release(s->f);
  gkyl_pkpm_lte_release(s->up);
}

// A fine grid and a wide domain: the projection has the moments of the Maxwellian.
void
test_proj_moments()
{
  struct lte_setup s = lte_setup_new(14.0, 256);
  double target[3] = {2.5, 0.3, 0.7}, m[3];
  gkyl_pkpm_lte_proj(s.up, target, s.f);
  gkyl_pkpm_lte_moments(s.up, s.f, m);
  TEST_CHECK(fabs(m[0] - 2.5) < 1e-12);
  TEST_CHECK(fabs(m[1] - 2.5 * 0.3) < 1e-12);
  TEST_CHECK(fabs(m[2] - 2.5 * (0.3 * 0.3 + 0.7)) < 1e-12);
  TEST_MSG("moments %.15g %.15g %.15g", m[0], m[1], m[2]);
  lte_setup_release(&s);
}

// The moments integrate the p=2 expansion exactly: f = 1 on the domain has M2 = 2 w_max^3/3.
void
test_moments_exact()
{
  double w_max = 3.0;
  struct lte_setup s = lte_setup_new(w_max, 7);
  gkyl_array_clear(s.f, 0.0);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(s.f, gkyl_range_idx(&s.range, iter.idx));
    d[0] = sqrt(2.0);
  }
  double m[3];
  gkyl_pkpm_lte_moments(s.up, s.f, m);
  TEST_CHECK(fabs(m[0] - 2.0 * w_max) < 1e-13);
  TEST_CHECK(fabs(m[1]) < 1e-13);
  TEST_CHECK(fabs(m[2] - 2.0 * w_max * w_max * w_max / 3.0) < 1e-13);
  lte_setup_release(&s);
}

// The correction reaches the target moments on truncated and coarse grids.
static void
check_correct(double w_max, int num_cells, double n, double u, double vtsq)
{
  struct lte_setup s = lte_setup_new(w_max, num_cells);
  double target[3] = {n, u, vtsq}, m[3];
  struct gkyl_pkpm_lte_status status = gkyl_pkpm_lte_correct(s.up, target, s.f);
  TEST_CHECK(status.converged);
  TEST_MSG(
    "w_max %g, Nw %d: error %g after %d iterations", w_max, num_cells, status.error, status.num_iter
  );
  gkyl_pkpm_lte_moments(s.up, s.f, m);
  TEST_CHECK(fabs(m[0] - n) < 2e-14 * n);
  TEST_CHECK(fabs(m[1] - n * u) < 2e-14 * n);
  TEST_CHECK(fabs(m[2] - n * (u * u + vtsq)) < 2e-14 * n * (u * u + vtsq));
  lte_setup_release(&s);
}

void
test_correct_unit()
{
  check_correct(5.0, 16, 1.0, 0.0, 1.0);
  check_correct(6.0, 8, 1.0, 0.0, 1.0);
  check_correct(6.0, 12, 1.0, 0.0, 1.0);
  check_correct(6.0, 64, 1.0, 0.0, 1.0);
  check_correct(8.0, 24, 1.0, 0.0, 1.0);
}

void
test_correct_general()
{
  check_correct(6.0, 32, 1.7, 0.4, 0.8);
  check_correct(8.0, 48, 0.3, -0.6, 1.5);
}

TEST_LIST = {
  {"proj_moments", test_proj_moments},
  {"moments_exact", test_moments_exact},
  {"correct_unit", test_correct_unit},
  {"correct_general", test_correct_general},
  {NULL, NULL}
};
