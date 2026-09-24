#include <gkyl_array.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_ops.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_util.h>
#include <gkyl_array_rio.h>
#include <acutest.h>
#include <gkyl_gyrokinetic_pol_density.h>

#include <stdio.h>
#include <stdlib.h>

// Allocate array (filled with zeros).
static struct gkyl_array*
mkarr(bool on_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (on_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}


void evalFunc1x_1(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = 1.0;
}

void evalFunc1x_quad(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = -pow(xn[0]-0.5, 2) + 1.0;
}

void evalFunc2x_quad(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = (-pow(xn[0]-0.5, 2) + 1.0) * (-pow(xn[1]-0.5, 2) + 1.0);
}

void evalFunc3x_quad(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = (-pow(xn[0]-0.5, 2) + 1.0) * (-pow(xn[1]-0.5, 2) + 1.0) * (-pow(xn[2]-0.5, 2) + 1.0);
}

void test_1x_flat( bool use_gpu ) 
{  
  int cells[] = {8};
  int poly_order = 1;
  double lower[] = {0.0}, upper[] = {1.0};
  int dim = sizeof(lower)/sizeof(lower[0]);

  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  int ghost[] = { 1, 1 };
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);

  struct gkyl_array *npol = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  struct gkyl_array *npol_ho = use_gpu? mkarr(false, npol->ncomp, npol->size)
                                      : gkyl_array_acquire(npol);
  struct gkyl_array *epsilon = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  struct gkyl_array *epsilon_ho = use_gpu? mkarr(false, epsilon->ncomp, epsilon->size)
                                         : gkyl_array_acquire(epsilon);

  struct gkyl_eval_on_nodes *epsilon_proj = gkyl_eval_on_nodes_new(&grid, &basis,
    1, evalFunc1x_1, NULL);
  gkyl_eval_on_nodes_advance(epsilon_proj, 0.0, &localRange, epsilon_ho);
  gkyl_eval_on_nodes_release(epsilon_proj);
  gkyl_array_copy(epsilon, epsilon_ho);

  struct gkyl_basis phi_pol_basis;
  gkyl_cart_modal_tensor(&phi_pol_basis, 1, poly_order+1);

  struct gkyl_array *phi_pol = mkarr(use_gpu, phi_pol_basis.num_basis, localRange_ext.volume);
  struct gkyl_array *phi_pol_ho = use_gpu? mkarr(false, phi_pol->ncomp, phi_pol->size)
                                         : gkyl_array_acquire(phi_pol);

  struct gkyl_eval_on_nodes *phi_pol_proj = gkyl_eval_on_nodes_new(&grid, &phi_pol_basis,
    1, evalFunc1x_1, NULL);
  gkyl_eval_on_nodes_advance(phi_pol_proj, 0.0, &localRange, phi_pol_ho);
  gkyl_eval_on_nodes_release(phi_pol_proj);

  gkyl_array_copy(phi_pol, phi_pol_ho);

  struct gkyl_gyrokinetic_pol_density* npol_op = gkyl_gyrokinetic_pol_density_new(basis, grid, use_gpu);
  gkyl_gyrokinetic_pol_density_advance(npol_op, &localRange, epsilon, phi_pol, npol);

  gkyl_array_copy(npol_ho, npol);

  // read the components of npol
  struct gkyl_range_iter conf_iter;
  gkyl_range_iter_init(&conf_iter, &localRange);
  while (gkyl_range_iter_next(&conf_iter)) {
    long linidx = gkyl_range_idx(&localRange, conf_iter.idx);
    double *npol_d = gkyl_array_fetch(npol_ho, linidx);
    TEST_CHECK( gkyl_compare(npol_d[0], sqrt(2), 1e-14) );
    TEST_CHECK( gkyl_compare(npol_d[1], 0.0, 1e-12) );
  }
  gkyl_array_release(npol_ho);
  gkyl_array_release(npol);
  gkyl_array_release(epsilon_ho);
  gkyl_array_release(epsilon);
  gkyl_array_release(phi_pol_ho);
  gkyl_array_release(phi_pol);
  gkyl_gyrokinetic_pol_density_release(npol_op);
}

void test_1x_flat_cpu() 
{
  test_1x_flat(false);
}

void test_1x_flat_gpu() 
{
  test_1x_flat(true);
}

void test_1x_quad( bool use_gpu ) 
{  
  int cells[] = {8};
  int poly_order = 1;
  double lower[] = {0.0}, upper[] = {1.0};
  int dim = sizeof(lower)/sizeof(lower[0]);
  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  int ghost[] = { 1, 1 };
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);

  struct gkyl_array *npol = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  struct gkyl_array *npol_ho = use_gpu? mkarr(false, npol->ncomp, npol->size)
                                      : gkyl_array_acquire(npol);
  struct gkyl_array *epsilon = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  struct gkyl_array *epsilon_ho = use_gpu? mkarr(false, epsilon->ncomp, epsilon->size)
                                         : gkyl_array_acquire(epsilon);

  struct gkyl_eval_on_nodes *epsilon_proj = gkyl_eval_on_nodes_new(&grid, &basis,
    1, evalFunc1x_1, NULL);
  gkyl_eval_on_nodes_advance(epsilon_proj, 0.0, &localRange, epsilon_ho);
  gkyl_eval_on_nodes_release(epsilon_proj);
  gkyl_array_copy(epsilon, epsilon_ho);

  struct gkyl_basis phi_pol_basis;
  gkyl_cart_modal_tensor(&phi_pol_basis, 1, poly_order+1);

  struct gkyl_array *phi_pol = mkarr(use_gpu, phi_pol_basis.num_basis, localRange_ext.volume);
  struct gkyl_array *phi_pol_ho = use_gpu? mkarr(false, phi_pol->ncomp, phi_pol->size)
                                         : gkyl_array_acquire(phi_pol);

  struct gkyl_eval_on_nodes *phi_pol_proj = gkyl_eval_on_nodes_new(&grid, &phi_pol_basis,
    1, evalFunc1x_quad, NULL);
  gkyl_eval_on_nodes_advance(phi_pol_proj, 0.0, &localRange, phi_pol_ho);
  gkyl_eval_on_nodes_release(phi_pol_proj);
  gkyl_array_copy(phi_pol, phi_pol_ho);

  struct gkyl_gyrokinetic_pol_density* npol_op = gkyl_gyrokinetic_pol_density_new(basis, grid, use_gpu);
  gkyl_gyrokinetic_pol_density_advance(npol_op, &localRange, epsilon, phi_pol, npol);

  gkyl_array_copy(npol_ho, npol);

  // read the components of npol
  struct gkyl_range_iter conf_iter;
  gkyl_range_iter_init(&conf_iter, &localRange);
  while (gkyl_range_iter_next(&conf_iter)) {
    long linidx = gkyl_range_idx(&localRange, conf_iter.idx);
    double *npol_d = gkyl_array_fetch(npol_ho, linidx);
    double tol = 1e-5;
    if (conf_iter.idx[0] == 1) {
      TEST_CHECK( gkyl_compare(npol_d[0], 1.141683, tol) );
      TEST_MSG("npol_d[0] = %g", npol_d[0]);
      TEST_CHECK( gkyl_compare(npol_d[1], 0.044652, tol) );
      TEST_MSG("npol_d[1] = %g", npol_d[1]);
    } else if (conf_iter.idx[0] == 2) {
      TEST_CHECK( gkyl_compare(npol_d[0], 1.274265, tol) );
      TEST_MSG("npol_d[0] = %g", npol_d[0]);
      TEST_CHECK( gkyl_compare(npol_d[1], 0.031894, tol) );
      TEST_MSG("npol_d[1] = %g", npol_d[1]);
    } else if (conf_iter.idx[0] == 3) {
      TEST_CHECK( gkyl_compare(npol_d[0], 1.362654, tol) );
      TEST_MSG("npol_d[0] = %g", npol_d[0]);
      TEST_CHECK( gkyl_compare(npol_d[1], 0.019137, tol) );
      TEST_MSG("npol_d[1] = %g", npol_d[1]);
    } else if (conf_iter.idx[0] == 4) {
      TEST_CHECK( gkyl_compare(npol_d[0], 1.406848, tol) );
      TEST_MSG("npol_d[0] = %g", npol_d[0]);
      TEST_CHECK( gkyl_compare(npol_d[1], 0.006379, tol) );
      TEST_MSG("npol_d[1] = %g", npol_d[1]);
    } else if (conf_iter.idx[0] == 5) {
      TEST_CHECK( gkyl_compare(npol_d[0], 1.406848, tol) );
      TEST_MSG("npol_d[0] = %g", npol_d[0]);
      TEST_CHECK( gkyl_compare(npol_d[1], -0.006379, tol) );
      TEST_MSG("npol_d[1] = %g", npol_d[1]);
    } else if (conf_iter.idx[0] == 6) {
      TEST_CHECK( gkyl_compare(npol_d[0], 1.362654, tol) );
      TEST_MSG("npol_d[0] = %g", npol_d[0]);
      TEST_CHECK( gkyl_compare(npol_d[1], -0.019137, tol) );
      TEST_MSG("npol_d[1] = %g", npol_d[1]);
    } else if (conf_iter.idx[0] == 7) {
      TEST_CHECK( gkyl_compare(npol_d[0], 1.274265, tol) );
      TEST_MSG("npol_d[0] = %g", npol_d[0]);
      TEST_CHECK( gkyl_compare(npol_d[1], -0.031894, tol) );
      TEST_MSG("npol_d[1] = %g", npol_d[1]);
    } else if (conf_iter.idx[0] == 8) {
      TEST_CHECK( gkyl_compare(npol_d[0], 1.141683, tol) );
      TEST_MSG("npol_d[0] = %g", npol_d[0]);
      TEST_CHECK( gkyl_compare(npol_d[1], -0.044652, tol) );
      TEST_MSG("npol_d[1] = %g", npol_d[1]);
    }
  }

  gkyl_array_release(npol_ho);
  gkyl_array_release(npol);
  gkyl_array_release(epsilon_ho);
  gkyl_array_release(epsilon);
  gkyl_array_release(phi_pol_ho);
  gkyl_array_release(phi_pol);
  gkyl_gyrokinetic_pol_density_release(npol_op);
}

void test_1x_quad_cpu() 
{
  test_1x_quad(false);
}

void test_1x_quad_gpu() 
{
  test_1x_quad(true);
}

void
test_2x_quad( bool use_gpu ) 
{
  int cells[] = {7, 7};
  int poly_order = 1;
  double lower[] = {0.0, 0.0}, upper[] = {1.0, 1.0};
  double time = 0.0;
  int dim = sizeof(lower)/sizeof(lower[0]);
  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  int ghost[] = { 1, 1 };
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);

  struct gkyl_array *npol = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  struct gkyl_array *npol_ho = use_gpu? mkarr(false, npol->ncomp, npol->size)
                                      : gkyl_array_acquire(npol);
  struct gkyl_array *epsilon = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  struct gkyl_array *epsilon_ho = use_gpu? mkarr(false, epsilon->ncomp, epsilon->size)
                                         : gkyl_array_acquire(epsilon);

  struct gkyl_eval_on_nodes *epsilon_proj = gkyl_eval_on_nodes_new(&grid, &basis,
    1, evalFunc1x_1, NULL);
  gkyl_eval_on_nodes_advance(epsilon_proj, time, &localRange, epsilon_ho);
  gkyl_eval_on_nodes_release(epsilon_proj);
  gkyl_array_copy(epsilon, epsilon_ho);

  struct gkyl_basis phi_pol_basis;
  gkyl_cart_modal_tensor(&phi_pol_basis, dim, poly_order+1);

  struct gkyl_array *phi_pol = mkarr(use_gpu, phi_pol_basis.num_basis, localRange_ext.volume);
  struct gkyl_array *phi_pol_ho = use_gpu? mkarr(false, phi_pol->ncomp, phi_pol->size)
                                         : gkyl_array_acquire(phi_pol);

  struct gkyl_eval_on_nodes *phi_pol_proj = gkyl_eval_on_nodes_new(&grid, &phi_pol_basis,
    1, evalFunc2x_quad, NULL);
  gkyl_eval_on_nodes_advance(phi_pol_proj, time, &localRange, phi_pol_ho);
  gkyl_eval_on_nodes_release(phi_pol_proj);
  gkyl_array_copy(phi_pol, phi_pol_ho);

  struct gkyl_gyrokinetic_pol_density* npol_op = gkyl_gyrokinetic_pol_density_new(basis, grid, use_gpu);
  gkyl_gyrokinetic_pol_density_advance(npol_op, &localRange, epsilon, phi_pol, npol);

  gkyl_array_copy(npol_ho, npol);

  // read the components of npol
  struct gkyl_range_iter conf_iter;
  gkyl_range_iter_init(&conf_iter, &localRange);
  while (gkyl_range_iter_next(&conf_iter)) {
    long linidx = gkyl_range_idx(&localRange, conf_iter.idx);
    double *npol_d = gkyl_array_fetch(npol_ho, linidx);
    double tol = 1e-12;
    if (conf_iter.idx[1] == 1) {
      TEST_CHECK( gkyl_compare(npol_d[0], 3.2585034013604615, tol) );
      TEST_CHECK( gkyl_compare(npol_d[1], 0.0, tol) );
      TEST_CHECK( gkyl_compare(npol_d[2], 0.1413919026586975, tol) );
      TEST_CHECK( gkyl_compare(npol_d[3], 0.0, tol) );
    } else if (conf_iter.idx[1] == 2) {
      TEST_CHECK( gkyl_compare(npol_d[0], 3.6666666666666150, tol) );
      TEST_CHECK( gkyl_compare(npol_d[1], 0.0, tol) );
      TEST_CHECK( gkyl_compare(npol_d[2], 0.0942612684391317, tol) );
      TEST_CHECK( gkyl_compare(npol_d[3], 0.0, tol) );
    } else if (conf_iter.idx[1] == 3) {
      TEST_CHECK( gkyl_compare(npol_d[0], 3.9115646258503931, tol) );
      TEST_CHECK( gkyl_compare(npol_d[1], 0.0, tol) );
      TEST_CHECK( gkyl_compare(npol_d[2], 0.0471306342195649, tol) );
      TEST_CHECK( gkyl_compare(npol_d[3], 0.0, tol) );
    } else if (conf_iter.idx[1] == 4) {
      TEST_CHECK( gkyl_compare(npol_d[0], 3.9931972789115648, tol) );
      TEST_CHECK( gkyl_compare(npol_d[1], 0.0, tol) );
      TEST_CHECK( gkyl_compare(npol_d[2], 0.0, tol) );
      TEST_CHECK( gkyl_compare(npol_d[3], 0.0, tol) );
    } else if (conf_iter.idx[1] == 5) {
      TEST_CHECK( gkyl_compare(npol_d[0], 3.9115646258503931, tol) );
      TEST_CHECK( gkyl_compare(npol_d[1], 0.0, tol) );
      TEST_CHECK( gkyl_compare(npol_d[2], -0.0471306342195649, tol) );
      TEST_CHECK( gkyl_compare(npol_d[3], 0.0, tol) );
    } else if (conf_iter.idx[1] == 6) {
      TEST_CHECK( gkyl_compare(npol_d[0], 3.6666666666666150, tol) );
      TEST_CHECK( gkyl_compare(npol_d[1], 0.0, tol) );
      TEST_CHECK( gkyl_compare(npol_d[2], -0.0942612684391317, tol) );
      TEST_CHECK( gkyl_compare(npol_d[3], 0.0, tol) );
    } else if (conf_iter.idx[1] == 7) {
      TEST_CHECK( gkyl_compare(npol_d[0], 3.2585034013604615, tol) );
      TEST_CHECK( gkyl_compare(npol_d[1], 0.0, tol) );
      TEST_CHECK( gkyl_compare(npol_d[2], -0.1413919026586975, tol) );
      TEST_CHECK( gkyl_compare(npol_d[3], 0.0, tol) );
    }
  }

  gkyl_array_release(npol_ho);
  gkyl_array_release(npol);
  gkyl_array_release(epsilon_ho);
  gkyl_array_release(epsilon);
  gkyl_array_release(phi_pol_ho);
  gkyl_array_release(phi_pol);
  gkyl_gyrokinetic_pol_density_release(npol_op);
}

void test_2x_quad_cpu() 
{
  test_2x_quad(false);
}

void test_2x_quad_gpu() 
{
  test_2x_quad(true);
}

void
test_3x_flat( bool use_gpu ) 
{
  int cells[] = {7, 7, 7};
  int poly_order = 1;
  double lower[] = {0.0, 0.0, 0.0}, upper[] = {1.0, 1.0, 1.0};
  double time = 0.0;
  int dim = sizeof(lower)/sizeof(lower[0]);
  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  int ghost[] = { 1, 1 };
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);

  struct gkyl_array *npol = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  struct gkyl_array *npol_ho = use_gpu? mkarr(false, npol->ncomp, npol->size)
                                      : gkyl_array_acquire(npol);
  struct gkyl_array *epsilon = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  struct gkyl_array *epsilon_ho = use_gpu? mkarr(false, epsilon->ncomp, epsilon->size)
                                         : gkyl_array_acquire(epsilon);

  struct gkyl_eval_on_nodes *epsilon_proj = gkyl_eval_on_nodes_new(&grid, &basis,
    1, evalFunc1x_1, NULL);
  gkyl_eval_on_nodes_advance(epsilon_proj, time, &localRange, epsilon_ho);
  gkyl_eval_on_nodes_release(epsilon_proj);
  gkyl_array_copy(epsilon, epsilon_ho);

  struct gkyl_basis phi_pol_basis;
  gkyl_cart_modal_tensor(&phi_pol_basis, dim, poly_order+1);

  struct gkyl_array *phi_pol = mkarr(use_gpu, phi_pol_basis.num_basis, localRange_ext.volume);
  struct gkyl_array *phi_pol_ho = use_gpu? mkarr(false, phi_pol->ncomp, phi_pol->size)
                                         : gkyl_array_acquire(phi_pol);

  struct gkyl_eval_on_nodes *phi_pol_proj = gkyl_eval_on_nodes_new(&grid, &phi_pol_basis,
    1, evalFunc1x_quad, NULL);
  gkyl_eval_on_nodes_advance(phi_pol_proj, time, &localRange, phi_pol_ho);
  gkyl_eval_on_nodes_release(phi_pol_proj);
  gkyl_array_copy(phi_pol, phi_pol_ho);

  struct gkyl_gyrokinetic_pol_density* npol_op = gkyl_gyrokinetic_pol_density_new(basis, grid, use_gpu);
  gkyl_gyrokinetic_pol_density_advance(npol_op, &localRange, epsilon, phi_pol, npol);

  gkyl_array_copy(npol_ho, npol);

  // read the components of npol
  struct gkyl_range_iter conf_iter;
  gkyl_range_iter_init(&conf_iter, &localRange);
  while (gkyl_range_iter_next(&conf_iter)) {
    long linidx = gkyl_range_idx(&localRange, conf_iter.idx);
    double *npol_d = gkyl_array_fetch(npol_ho, linidx);
    double tol = 1e-12;
    TEST_CHECK( gkyl_compare(npol_d[0], 5.6568542494927714, tol) );
    TEST_CHECK( gkyl_compare(npol_d[1], 0.0, tol) );
    TEST_CHECK( gkyl_compare(npol_d[2], 0.0, tol) );
    TEST_CHECK( gkyl_compare(npol_d[3], 0.0, tol) );
    TEST_CHECK( gkyl_compare(npol_d[4], 0.0, tol) );
    TEST_CHECK( gkyl_compare(npol_d[5], 0.0, tol) );
    TEST_CHECK( gkyl_compare(npol_d[6], 0.0, tol) );
    TEST_CHECK( gkyl_compare(npol_d[7], 0.0, tol) );
  }

  gkyl_array_release(npol_ho);
  gkyl_array_release(npol);
  gkyl_array_release(epsilon_ho);
  gkyl_array_release(epsilon);
  gkyl_array_release(phi_pol_ho);
  gkyl_array_release(phi_pol);
  gkyl_gyrokinetic_pol_density_release(npol_op);
}

void
test_3x_flat_cpu() 
{
  test_3x_flat(false);
}

void
test_3x_flat_gpu() 
{
  test_3x_flat(true);
}

TEST_LIST = {
  { "test_1x_flat_cpu", test_1x_flat_cpu },
  { "test_1x_quad_cpu", test_1x_quad_cpu },
  { "test_2x_quad_cpu", test_2x_quad_cpu },
  { "test_3x_flat_cpu", test_3x_flat_cpu },
#ifdef GKYL_HAVE_CUDA
  { "test_1x_flat_gpu", test_1x_flat_gpu },
  { "test_1x_quad_gpu", test_1x_quad_gpu },
  { "test_2x_quad_gpu", test_2x_quad_gpu },
  { "test_3x_flat_gpu", test_3x_flat_gpu },
#endif
  { NULL, NULL },
};
