#include "gkyl_array.h"
#include <acutest.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_dg_calc_gk_neut_hamil.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_ops.h>
#include <stdio.h>

// Function to allocate a gkyl array, zero-initialized, on CPU or GPU
static struct gkyl_array* mkarr(bool on_gpu, long nc, long size) {
  return on_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
                : gkyl_array_new(GKYL_DOUBLE, nc, size);
}

void eval_gxx(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0;
}
void eval_gyy(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 2.0;
}
void eval_gzz(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 3.0;
}

void
test_hamil(int cdim, bool use_gpu)
{
  // construct phase grid
  // construct phase basis
  // construct conf and phase ranges

  double xmax = 1.0;
  double vmax = 1.0;
  int poly_order = 1;
  int vdim = 3;
  int pdim = cdim + vdim;

  double lower[pdim], upper[pdim];
  int cells[pdim], ghost[pdim];
  int Nx = 4;
  int Nv = 16;

  for(int i=0; i<pdim; ++i) {
    if (i<cdim) {
      lower[i] = -xmax;
      upper[i] =  xmax;
      cells[i] = Nx;
    }
    else {
      lower[i] = -vmax;
      upper[i] =  vmax;
      cells[i] = Nv;
    }
    ghost[i] = 0;
  }

  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&confGrid, ghost, &confRange_ext, &confRange);

  struct gkyl_rect_grid grid;
  struct gkyl_range range, range_ext;
  gkyl_rect_grid_init(&grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&grid, ghost, &range_ext, &range);

  // basis functions
  struct gkyl_basis pbasis, cbasis; // phase-space, conf-space basis

  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  gkyl_cart_modal_tensor(&pbasis, pdim, poly_order);

  // construct gji gkyl_array
  gkyl_proj_on_basis *proj_gxx = gkyl_proj_on_basis_new(&confGrid, &cbasis,
    poly_order+1, 1, eval_gxx, NULL);
  gkyl_proj_on_basis *proj_gyy = gkyl_proj_on_basis_new(&confGrid, &cbasis,
    poly_order+1, 1, eval_gyy, NULL);
  gkyl_proj_on_basis *proj_gzz = gkyl_proj_on_basis_new(&confGrid, &cbasis,
    poly_order+1, 1, eval_gzz, NULL);
  
  struct gkyl_array *gij = mkarr(use_gpu, 6*cbasis.num_basis, confRange.volume);
  struct gkyl_array *gxx = mkarr(use_gpu, cbasis.num_basis, confRange.volume);
  struct gkyl_array *gxy = mkarr(use_gpu, cbasis.num_basis, confRange.volume);
  struct gkyl_array *gxz = mkarr(use_gpu, cbasis.num_basis, confRange.volume);
  struct gkyl_array *gyy = mkarr(use_gpu, cbasis.num_basis, confRange.volume);
  struct gkyl_array *gyz = mkarr(use_gpu, cbasis.num_basis, confRange.volume);
  struct gkyl_array *gzz = mkarr(use_gpu, cbasis.num_basis, confRange.volume);

  struct gkyl_array *hamil, *hamil_ho;
  hamil = mkarr(use_gpu, pbasis.num_basis, range_ext.volume);
  hamil_ho = use_gpu ? mkarr(false, hamil->ncomp, hamil->size)
                     : gkyl_array_acquire(hamil);

  struct gkyl_array *gxx_ho, *gyy_ho, *gzz_ho;
  gxx_ho = use_gpu ? mkarr(false, gxx->ncomp, gxx->size)
                     : gkyl_array_acquire(gxx);
  gyy_ho = use_gpu ? mkarr(false, gyy->ncomp, gyy->size)
                     : gkyl_array_acquire(gyy);
  gzz_ho = use_gpu ? mkarr(false, gzz->ncomp, gzz->size)
                     : gkyl_array_acquire(gzz);

  gkyl_proj_on_basis_advance(proj_gxx, 1.0, &confRange, gxx_ho);
  gkyl_proj_on_basis_advance(proj_gyy, 1.0, &confRange, gyy_ho);
  gkyl_proj_on_basis_advance(proj_gzz, 1.0, &confRange, gzz_ho);

  if (use_gpu) {
    gkyl_array_copy(gxx, gxx_ho);
    gkyl_array_copy(gyy, gyy_ho);
    gkyl_array_copy(gzz, gzz_ho);
  }

  gkyl_array_set_offset(gij, 1.0, gxx, 0);
  gkyl_array_set_offset(gij, 1.0, gxy, cbasis.num_basis);
  gkyl_array_set_offset(gij, 1.0, gxz, 2*cbasis.num_basis);
  gkyl_array_set_offset(gij, 1.0, gyy, 3*cbasis.num_basis);
  gkyl_array_set_offset(gij, 1.0, gyz, 4*cbasis.num_basis);
  gkyl_array_set_offset(gij, 1.0, gzz, 5*cbasis.num_basis);

  struct gkyl_dg_calc_gk_neut_hamil* hamil_calc = gkyl_dg_calc_gk_neut_hamil_new(&grid, &pbasis, cdim, use_gpu);
  gkyl_dg_calc_gk_neut_hamil_calc(hamil_calc, &confRange, &range, gij, hamil);

  if (use_gpu) {
    gkyl_array_copy(hamil_ho, hamil);
  }

  /* char fname[1024]; */
  /* sprintf(fname, "ctest_gk_neut_hamil_%dx.gkyl", cdim); */
  /* gkyl_grid_sub_array_write(&grid, &range, 0, hamil, fname); */
  
  // test against predicted value
  if (cdim==3) {
    const double *fv = gkyl_array_cfetch(hamil_ho, gkyl_range_idx(&range_ext, (int[6]){1, 1, 1, 1, 1, 1}));
    double p1_vals[] = { 2.1125000000000000e+01, -2.2204460492503131e-16,  3.3306690738754696e-16,
			 2.2204460492503131e-16, -2.7063293868263760e-01, -5.4126587736527476e-01,
			 -8.1189881604791214e-01, -8.3266726846886741e-17, -1.6653345369377348e-16,
			 -1.3877787807814457e-16,  0.0000000000000000e+00,  2.7755575615628914e-16,
			 -8.3266726846886741e-17,  0.0000000000000000e+00,  2.7755575615628914e-17,
			 0.0000000000000000e+00,  0.0000000000000000e+00,  5.2735593669694936e-16,
			 -1.9428902930940239e-16, -2.2204460492503131e-16,  0.0000000000000000e+00,
			 2.2204460492503131e-16,  1.3877787807814457e-16, -2.7755575615628914e-17,
			 2.7755575615628914e-17, -2.7755575615628914e-17, -2.7755575615628914e-17,
			 -1.1102230246251565e-16,  0.0000000000000000e+00, -5.5511151231257827e-17,
			 0.0000000000000000e+00, -1.1102230246251565e-16,  0.0000000000000000e+00,
			 0.0000000000000000e+00,  0.0000000000000000e+00,  5.5511151231257827e-17,
			 -2.2204460492503131e-16,  2.2204460492503131e-16, -2.7755575615628914e-17,
			 1.6653345369377348e-16, -2.2204460492503131e-16, -5.5511151231257827e-17,
			 -6.9388939039072284e-17, -4.8572257327350599e-17,  3.4694469519536142e-17,
			 -4.8572257327350599e-17, -1.1102230246251565e-16, -1.1102230246251565e-16,
			 -5.5511151231257827e-17,  0.0000000000000000e+00,  5.5511151231257827e-17,
			 -1.1102230246251565e-16, -1.3877787807814457e-16,  5.5511151231257827e-17,
			 4.1633363423443370e-17,  5.5511151231257827e-17, -2.7755575615628914e-17,
			 -6.9388939039072284e-18,  1.7347234759768071e-17,  2.0816681711721685e-17,
			 1.3877787807814457e-17,  2.7755575615628914e-17, -2.0816681711721685e-17,
			 -8.7378663975128048e-18};
    for (int i=0; i<pbasis.num_basis; ++i) {
      TEST_CHECK( gkyl_compare_double(p1_vals[i], fv[i], 1e-12) ); 
    }
  }

  gkyl_array_release(hamil); 
  gkyl_array_release(gij); 
  gkyl_array_release(gxx); 
  gkyl_array_release(gxy); 
  gkyl_array_release(gxz); 
  gkyl_array_release(gyy); 
  gkyl_array_release(gyz); 
  gkyl_array_release(gzz); 
  gkyl_array_release(hamil_ho);
  gkyl_array_release(gxx_ho);
  gkyl_array_release(gyy_ho);
  gkyl_array_release(gzz_ho);
  gkyl_proj_on_basis_release(proj_gxx); 
  gkyl_proj_on_basis_release(proj_gyy);
  gkyl_proj_on_basis_release(proj_gzz);
  gkyl_dg_calc_gk_neut_hamil_release(hamil_calc);
}

void test_gkneut_hamil_3x_ho() { test_hamil(3, false); }

#ifdef GKYL_HAVE_CUDA
void test_gkneut_hamil_3x_dev() { test_hamil(3, true); }
#endif

TEST_LIST = {
  { "test_gkneut_hamil_3x_ho", test_gkneut_hamil_3x_ho },
#ifdef GKYL_HAVE_CUDA
  { "test_gkneut_hamil_3x_dev", test_gkneut_hamil_3x_dev },
#endif
  { NULL, NULL },
};
