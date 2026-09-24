#include <gkyl_array.h>
#include <gkyl_util.h>
#include <acutest.h>
#include <gkyl_array_rio.h>
#include <gkyl_proj_powsqrt_on_basis.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>

#include "math.h"

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

void eval_fun_1x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double Lx = 2.*M_PI;
  fout[0] = 0.5*(1.+cos(0.5*2.*M_PI*x/Lx));
}
void eval_powsqrt_fun_1x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double Lx = 2.*M_PI;
  double exponent = 3;
  double fun[1];
  eval_fun_1x(t, xn, fun, ctx);
  fout[0] = pow( sqrt(fun[0]), exponent);
}

void
test_1x(int poly_order, bool use_gpu)
{
  double Lx = 2.*M_PI;
  double lower[] = {-Lx/2.}, upper[] = {Lx/2.};
  int cells[] = {32};

  int ndim = sizeof(lower)/sizeof(lower[0]);

  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  int ghost[] = { 1 };
  struct gkyl_range local, local_ext; // Local, local-ext phase-space ranges.
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create f, g=pow(sqrt(f), q) arrays.
  struct gkyl_array *f, *g;
  f = mkarr(basis.num_basis, local_ext.volume);
  g = mkarr(basis.num_basis, local_ext.volume);
  struct gkyl_array *f_cu, *g_cu;
  if (use_gpu) { // Create device copies
    f_cu  = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
    g_cu  = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  }

  gkyl_proj_on_basis *proj_f = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_fun_1x, NULL);

  gkyl_proj_on_basis_advance(proj_f, 0.0, &local, f);

  if (use_gpu) {
    // Copy host array to device.
    gkyl_array_copy(f_cu, f);
  }

  // Create pow(sqrt( ), ) updater.
  gkyl_proj_powsqrt_on_basis *proj_up = gkyl_proj_powsqrt_on_basis_new(&basis, poly_order+1, use_gpu);

  if (use_gpu) {
    gkyl_proj_powsqrt_on_basis_advance(proj_up, &local, 3, f_cu, g_cu);
    gkyl_array_copy(g, g_cu);
  } else {
    gkyl_proj_powsqrt_on_basis_advance(proj_up, &local, 3, f, g);
  }

  // Project expected g.
  struct gkyl_array *gA;
  gA = mkarr(basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis *proj_g = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_powsqrt_fun_1x, NULL);
  gkyl_proj_on_basis_advance(proj_g, 0.0, &local, gA);

  for (int k=0; k<cells[0]; k++) {
    int idx[] = {k+1};
    long linidx = gkyl_range_idx(&local, idx);
    const double *g_p = gkyl_array_cfetch(g, linidx);
    const double *gA_p = gkyl_array_cfetch(gA, linidx);
    for (int m=0; m<basis.num_basis; m++) {
      TEST_CHECK( gkyl_compare(gA_p[m], g_p[m], 1e-10) );
      TEST_MSG("Expected: %.13e in cell (%d)", gA_p[m], idx[0]);
      TEST_MSG("Produced: %.13e", g_p[m]);
    }
  }

  gkyl_array_release(f); gkyl_array_release(g); gkyl_array_release(gA);
  if (use_gpu) {
    gkyl_array_release(f_cu); gkyl_array_release(g_cu);
  }

  gkyl_proj_on_basis_release(proj_f);
  gkyl_proj_on_basis_release(proj_g);

  gkyl_proj_powsqrt_on_basis_release(proj_up);
}

void eval_fun_2x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double Lx = 2.*M_PI, Ly = 8.*M_PI;
  fout[0] = 0.5*(1.+cos(0.5*2.*M_PI*x/Lx))*(y+Ly);
}
void eval_powsqrt_fun_2x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double Lx = 2.*M_PI, Ly = 8.*M_PI;
  double exponent = 3;
  double fun[1];
  eval_fun_2x(t, xn, fun, ctx);
  fout[0] = pow( sqrt(fun[0]), exponent);
}

void
test_2x(int poly_order, bool use_gpu)
{
  double Lx = 2.*M_PI, Ly = 8.*M_PI;
  double lower[] = {-Lx/2., -Ly/2.}, upper[] = {Lx/2., Ly/2.};
  int cells[] = {16, 64};

  int ndim = sizeof(lower)/sizeof(lower[0]);

  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  int ghost[] = { 1, 1 };
  struct gkyl_range local, local_ext; // Local, local-ext phase-space ranges.
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create f, g=pow(sqrt(f), q) arrays.
  struct gkyl_array *f, *g;
  f = mkarr(basis.num_basis, local_ext.volume);
  g = mkarr(basis.num_basis, local_ext.volume);
  struct gkyl_array *f_cu, *g_cu;
  if (use_gpu) { // Create device copies
    f_cu  = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
    g_cu  = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  }

  gkyl_proj_on_basis *proj_f = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_fun_2x, NULL);

  gkyl_proj_on_basis_advance(proj_f, 0.0, &local, f);

  if (use_gpu) {
    // Copy host array to device.
    gkyl_array_copy(f_cu, f);
  }

  // Create pow(sqrt( ), ) updater.
  gkyl_proj_powsqrt_on_basis *proj_up = gkyl_proj_powsqrt_on_basis_new(&basis, poly_order+1, use_gpu);

  if (use_gpu) {
    gkyl_proj_powsqrt_on_basis_advance(proj_up, &local, 3, f_cu, g_cu);
    gkyl_array_copy(g, g_cu);
  } else {
    gkyl_proj_powsqrt_on_basis_advance(proj_up, &local, 3, f, g);
  }

  // Project expected g.
  struct gkyl_array *gA;
  gA = mkarr(basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis *proj_g = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_powsqrt_fun_2x, NULL);
  gkyl_proj_on_basis_advance(proj_g, 0.0, &local, gA);

  for (int j=0; j<cells[0]; j++) {
    for (int k=0; k<cells[0]; k++) {
      int idx[] = {j+1,k+1};
      long linidx = gkyl_range_idx(&local, idx);
      const double *g_p = gkyl_array_cfetch(g, linidx);
      const double *gA_p = gkyl_array_cfetch(gA, linidx);
      for (int m=0; m<basis.num_basis; m++) {
        TEST_CHECK( gkyl_compare(gA_p[m], g_p[m], 4e-9) );
        TEST_MSG("Expected: %.13e in cell (%d)", gA_p[m], idx[0]);
        TEST_MSG("Produced: %.13e", g_p[m]);
      }
    }
  }

  gkyl_array_release(f); gkyl_array_release(g); gkyl_array_release(gA);
  if (use_gpu) {
    gkyl_array_release(f_cu); gkyl_array_release(g_cu);
  }

  gkyl_proj_on_basis_release(proj_f);
  gkyl_proj_on_basis_release(proj_g);

  gkyl_proj_powsqrt_on_basis_release(proj_up);
}

void eval_fun_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double Lx = 2.*M_PI, Ly = 8.*M_PI, Lz = 4.*M_PI;
  fout[0] = 0.5*(1.+cos(0.5*2.*M_PI*x/Lx))*(y+Ly)*(Lz/2.-0.5*fabs(z));
}
void eval_powsqrt_fun_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double Lx = 2.*M_PI, Ly = 8.*M_PI, Lz = 4.*M_PI;
  double exponent = 3;
  double fun[1];
  eval_fun_3x(t, xn, fun, ctx);
  fout[0] = pow( sqrt(fun[0]), exponent);
}

void
test_3x(int poly_order, bool use_gpu)
{
  double Lx = 2.*M_PI, Ly = 8.*M_PI, Lz = 4.*M_PI;
  double lower[] = {-Lx/2., -Ly/2., -Lz/2.}, upper[] = {Lx/2., Ly/2., Lz/2.};
  int cells[] = {16, 64, 32};

  int ndim = sizeof(lower)/sizeof(lower[0]);

  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  int ghost[] = { 1, 1, 1 };
  struct gkyl_range local, local_ext; // Local, local-ext phase-space ranges.
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create f, g=pow(sqrt(f), q) arrays.
  struct gkyl_array *f, *g;
  f = mkarr(basis.num_basis, local_ext.volume);
  g = mkarr(basis.num_basis, local_ext.volume);
  struct gkyl_array *f_cu, *g_cu;
  if (use_gpu) { // Create device copies
    f_cu  = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
    g_cu  = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  }

  gkyl_proj_on_basis *proj_f = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_fun_3x, NULL);

  gkyl_proj_on_basis_advance(proj_f, 0.0, &local, f);

  if (use_gpu) {
    // Copy host array to device.
    gkyl_array_copy(f_cu, f);
  }

  // Create pow(sqrt( ), ) updater.
  gkyl_proj_powsqrt_on_basis *proj_up = gkyl_proj_powsqrt_on_basis_new(&basis, poly_order+1, use_gpu);

  if (use_gpu) {
    gkyl_proj_powsqrt_on_basis_advance(proj_up, &local, 3, f_cu, g_cu);
    gkyl_array_copy(g, g_cu);
  } else {
    gkyl_proj_powsqrt_on_basis_advance(proj_up, &local, 3, f, g);
  }

  // Project expected g.
  struct gkyl_array *gA;
  gA = mkarr(basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis *proj_g = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_powsqrt_fun_3x, NULL);
  gkyl_proj_on_basis_advance(proj_g, 0.0, &local, gA);

  for (int i=0; i<cells[0]; i++) {
    for (int j=0; j<cells[0]; j++) {
      for (int k=0; k<cells[0]; k++) {
        int idx[] = {i+1,j+1,k+1};
        long linidx = gkyl_range_idx(&local, idx);
        const double *g_p = gkyl_array_cfetch(g, linidx);
        const double *gA_p = gkyl_array_cfetch(gA, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(gA_p[m], g_p[m], 1e-7) );
          TEST_MSG("Expected: %.13e in cell (%d)", gA_p[m], idx[0]);
          TEST_MSG("Produced: %.13e", g_p[m]);
        }
      }
    }
  }

  gkyl_array_release(f); gkyl_array_release(g); gkyl_array_release(gA);
  if (use_gpu) {
    gkyl_array_release(f_cu); gkyl_array_release(g_cu);
  }

  gkyl_proj_on_basis_release(proj_f);
  gkyl_proj_on_basis_release(proj_g);

  gkyl_proj_powsqrt_on_basis_release(proj_up);
}

void test_1x_p1() { test_1x(1, false); }
void test_1x_p2() { test_1x(2, false); }

void test_2x_p1() { test_2x(1, false); }
void test_2x_p2() { test_2x(2, false); }

void test_3x_p1() { test_3x(1, false); }
void test_3x_p2() { test_3x(2, false); }

#ifdef GKYL_HAVE_CUDA
void test_1x_p1_gpu() { test_1x(1, true); }
void test_1x_p2_gpu() { test_1x(2, true); }

void test_2x_p1_gpu() { test_2x(1, true); }
void test_2x_p2_gpu() { test_2x(2, true); }

void test_3x_p1_gpu() { test_3x(1, true); }
void test_3x_p2_gpu() { test_3x(2, true); }
#endif

TEST_LIST = {
  { "test_1x_p1", test_1x_p1 },
  { "test_1x_p2", test_1x_p2 },

  { "test_2x_p1", test_2x_p1 },
  { "test_2x_p2", test_2x_p2 },

  { "test_3x_p1", test_3x_p1 },
  { "test_3x_p2", test_3x_p2 },
#ifdef GKYL_HAVE_CUDA
  { "test_1x_p1_gpu", test_1x_p1_gpu },
  { "test_1x_p2_gpu", test_1x_p2_gpu },

  { "test_2x_p1_gpu", test_2x_p1_gpu },
  { "test_2x_p2_gpu", test_2x_p2_gpu },

  { "test_3x_p1_gpu", test_3x_p1_gpu },
  { "test_3x_p2_gpu", test_3x_p2_gpu },
#endif
  { NULL, NULL },
};
