#include <acutest.h>
#include <math.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>

void f_1d(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 2 + x;
}
void g_1d(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 2*x*x + 8;
}

void fv2_1d(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 2. + x;
  fout[1] = 1. - x;
}
void gv2_1d(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 2*x*x + 8;
  fout[1] = x*x - 8;
}

void fv3_1d(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 2. + x;
  fout[1] = 1. - x;
  fout[2] = 0.5 + x;
}
void gv3_1d(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 2*x*x + 8;
  fout[1] = x*x - 8;
  fout[2] = 0.5*x*x + 4;
}

void check_dot_product_1d(const double *fv1_d, const double *gv1_d, const double *fvdgv1_d,
                          const double *fv2_d, const double *gv2_d, const double *fvdgv2_d,
                          const double *fv3_d, const double *gv3_d, const double *fvdgv3_d,
                          int poly_order) {
  if (poly_order == 1) {
    TEST_CHECK( gkyl_compare(fvdgv1_d[0], 0.7071067811865475*(fv1_d[1]*gv1_d[1]+fv1_d[0]*gv1_d[0]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv1_d[1], 0.7071067811865475*(fv1_d[0]*gv1_d[1]+gv1_d[0]*fv1_d[1]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv2_d[0], 0.7071067811865475*(fv2_d[3]*gv2_d[3]+fv2_d[2]*gv2_d[2]+fv2_d[1]*gv2_d[1]+fv2_d[0]*gv2_d[0]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv2_d[1], 0.7071067811865475*(fv2_d[2]*gv2_d[3]+gv2_d[2]*fv2_d[3]+fv2_d[0]*gv2_d[1]+gv2_d[0]*fv2_d[1]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv3_d[0], 0.7071067811865475*(fv3_d[5]*gv3_d[5]+fv3_d[4]*gv3_d[4]+fv3_d[3]*gv3_d[3]+fv3_d[2]*gv3_d[2]+fv3_d[1]*gv3_d[1]+fv3_d[0]*gv3_d[0]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv3_d[1], 0.7071067811865475*(fv3_d[4]*gv3_d[5]+gv3_d[4]*fv3_d[5]+fv3_d[2]*gv3_d[3]+gv3_d[2]*fv3_d[3]+fv3_d[0]*gv3_d[1]+gv3_d[0]*fv3_d[1]), 1e-12) );
  } else if (poly_order == 2) {
    TEST_CHECK( gkyl_compare(fvdgv1_d[0], 0.7071067811865475*(fv1_d[2]*gv1_d[2]+fv1_d[1]*gv1_d[1]+fv1_d[0]*gv1_d[0]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv1_d[1], 0.1414213562373095*(4.47213595499958*fv1_d[1]*gv1_d[2]+4.47213595499958*gv1_d[1]*fv1_d[2]+5.0*fv1_d[0]*gv1_d[1]+5.0*gv1_d[0]*fv1_d[1]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv1_d[2], 0.02020305089104421*((22.3606797749979*fv1_d[2]+35.0*fv1_d[0])*gv1_d[2]+35.0*gv1_d[0]*fv1_d[2]+31.30495168499706*fv1_d[1]*gv1_d[1]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv2_d[0], 0.7071067811865475*(fv2_d[5]*gv2_d[5]+fv2_d[4]*gv2_d[4]+fv2_d[3]*gv2_d[3]+fv2_d[2]*gv2_d[2]+fv2_d[1]*gv2_d[1]+fv2_d[0]*gv2_d[0]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv2_d[1], 0.1414213562373095*(4.47213595499958*fv2_d[4]*gv2_d[5]+4.47213595499958*gv2_d[4]*fv2_d[5]+5.0*fv2_d[3]*gv2_d[4]+5.0*gv2_d[3]*fv2_d[4]+4.47213595499958*fv2_d[1]*gv2_d[2]+4.47213595499958*gv2_d[1]*fv2_d[2]+5.0*fv2_d[0]*gv2_d[1]+5.0*gv2_d[0]*fv2_d[1]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv2_d[2], 0.02020305089104421*((22.3606797749979*fv2_d[5]+35.0*fv2_d[3])*gv2_d[5]+35.0*gv2_d[3]*fv2_d[5]+31.30495168499706*fv2_d[4]*gv2_d[4]+(22.3606797749979*fv2_d[2]+35.0*fv2_d[0])*gv2_d[2]+35.0*gv2_d[0]*fv2_d[2]+31.30495168499706*fv2_d[1]*gv2_d[1]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv3_d[0], 0.7071067811865475*(fv3_d[8]*gv3_d[8]+fv3_d[7]*gv3_d[7]+fv3_d[6]*gv3_d[6]+fv3_d[5]*gv3_d[5]+fv3_d[4]*gv3_d[4]+fv3_d[3]*gv3_d[3]+fv3_d[2]*gv3_d[2]+fv3_d[1]*gv3_d[1]+fv3_d[0]*gv3_d[0]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv3_d[1], 0.1414213562373095*(4.47213595499958*fv3_d[7]*gv3_d[8]+4.47213595499958*gv3_d[7]*fv3_d[8]+5.0*fv3_d[6]*gv3_d[7]+5.0*gv3_d[6]*fv3_d[7]+4.47213595499958*fv3_d[4]*gv3_d[5]+4.47213595499958*gv3_d[4]*fv3_d[5]+5.0*fv3_d[3]*gv3_d[4]+5.0*gv3_d[3]*fv3_d[4]+4.47213595499958*fv3_d[1]*gv3_d[2]+4.47213595499958*gv3_d[1]*fv3_d[2]+5.0*fv3_d[0]*gv3_d[1]+5.0*gv3_d[0]*fv3_d[1]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv3_d[2], 0.02020305089104421*((22.3606797749979*fv3_d[8]+35.0*fv3_d[6])*gv3_d[8]+35.0*gv3_d[6]*fv3_d[8]+31.30495168499706*fv3_d[7]*gv3_d[7]+(22.3606797749979*fv3_d[5]+35.0*fv3_d[3])*gv3_d[5]+35.0*gv3_d[3]*fv3_d[5]+31.30495168499706*fv3_d[4]*gv3_d[4]+(22.3606797749979*fv3_d[2]+35.0*fv3_d[0])*gv3_d[2]+35.0*gv3_d[0]*fv3_d[2]+31.30495168499706*fv3_d[1]*gv3_d[1]), 1e-12) );
  }
  return;
}

void
test_1d(int poly_order, bool use_gpu)
{
  double lower[] = {0.0}, upper[] = {1.0};
  int cells[] = {2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // projection updater for dist-function
  struct gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1, f_1d, NULL);
  struct gkyl_proj_on_basis *projDistg = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1, g_1d, NULL);
  // projection updaters for vector fields.
  struct gkyl_proj_on_basis *projfv2 = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 2, fv2_1d, NULL);
  struct gkyl_proj_on_basis *projgv2 = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 2, gv2_1d, NULL);
  struct gkyl_proj_on_basis *projfv3 = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 3, fv3_1d, NULL);
  struct gkyl_proj_on_basis *projgv3 = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 3, gv3_1d, NULL);

  // create array range: no ghost-cells in velocity space
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);  

  // create distribution function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *distg = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *distf_cu, *distg_cu;
  if (use_gpu) {
    distf_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    distg_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  }

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &arr_range, distf);
  gkyl_proj_on_basis_advance(projDistg, 0.0, &arr_range, distg);

  if (use_gpu) {
    // copy host array to device
    gkyl_array_copy(distf_cu, distf);
    gkyl_array_copy(distg_cu, distg);
  }

  struct gkyl_array *f_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *g_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(f_bar, 0.0);
  gkyl_array_clear(g_bar, 0.0);

  // Vector fields for dot product.
  struct gkyl_array *fv1 = gkyl_array_new(GKYL_DOUBLE, 1*basis.num_basis, arr_range.volume);
  struct gkyl_array *gv1 = gkyl_array_new(GKYL_DOUBLE, 1*basis.num_basis, arr_range.volume);
  struct gkyl_array *fv2 = gkyl_array_new(GKYL_DOUBLE, 2*basis.num_basis, arr_range.volume);
  struct gkyl_array *gv2 = gkyl_array_new(GKYL_DOUBLE, 2*basis.num_basis, arr_range.volume);
  struct gkyl_array *fv3 = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, arr_range.volume);
  struct gkyl_array *gv3 = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, arr_range.volume);
  struct gkyl_array *fv1_cu, *gv1_cu, *fv2_cu, *gv2_cu, *fv3_cu, *gv3_cu;
  if (use_gpu) {
    fv1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1*basis.num_basis, arr_range.volume);
    gv1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1*basis.num_basis, arr_range.volume);
    fv2_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 2*basis.num_basis, arr_range.volume);
    gv2_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 2*basis.num_basis, arr_range.volume);
    fv3_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, arr_range.volume);
    gv3_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, arr_range.volume);
  }

  // project vector fields on basis.
  gkyl_proj_on_basis_advance(projDistf, 0.0, &arr_range, fv1);
  gkyl_proj_on_basis_advance(projDistg, 0.0, &arr_range, gv1);
  gkyl_proj_on_basis_advance(projfv2, 0.0, &arr_range, fv2);
  gkyl_proj_on_basis_advance(projgv2, 0.0, &arr_range, gv2);
  gkyl_proj_on_basis_advance(projfv3, 0.0, &arr_range, fv3);
  gkyl_proj_on_basis_advance(projgv3, 0.0, &arr_range, gv3);
  if (use_gpu) {
    // copy host array to device
    gkyl_array_copy(fv1_cu, fv1);
    gkyl_array_copy(gv1_cu, gv1);
    gkyl_array_copy(fv2_cu, fv2);
    gkyl_array_copy(gv2_cu, gv2);
    gkyl_array_copy(fv3_cu, fv3);
    gkyl_array_copy(gv3_cu, gv3);
  }

  struct gkyl_array *h, *f_bar_cu, *g_bar_cu, *h_cu;
  struct gkyl_array *fvdgv1, *fvdgv2, *fvdgv3, *fvdgv1_cu, *fvdgv2_cu, *fvdgv3_cu;
  if (use_gpu) {
    f_bar_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    g_bar_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    gkyl_array_clear(f_bar_cu, 0.0);
    gkyl_array_clear(g_bar_cu, 0.0);

    // Product array only needs to be initialized on GPU.
    h_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    gkyl_array_clear(h_cu, 0.0);

    // Dot product of fv and gv on the GPU.
    fvdgv1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    fvdgv2_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    fvdgv3_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    gkyl_array_clear(fvdgv1_cu, 0.0);
    gkyl_array_clear(fvdgv2_cu, 0.0);
    gkyl_array_clear(fvdgv3_cu, 0.0);
  } else {
    h = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    gkyl_array_clear(h, 0.0);
  }

  // Dot product of fv and gv.
  fvdgv1 = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  fvdgv2 = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  fvdgv3 = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(fvdgv1, 0.0);
  gkyl_array_clear(fvdgv2, 0.0);
  gkyl_array_clear(fvdgv3, 0.0);

  gkyl_dg_bin_op_mem *mem;
  if (use_gpu) {
    // allocate memory
    mem = gkyl_dg_bin_op_mem_cu_dev_new(f_bar->size, basis.num_basis);
  
    // h = f*g
    gkyl_dg_mul_op(&basis, 0, h_cu, 0, distf_cu, 0, distg_cu);
    // f_bar = h/g = f
    gkyl_dg_div_op(mem, &basis, 0, f_bar_cu, 0, h_cu, 0, distg_cu);
    // g_bar = h/f = g
    gkyl_dg_div_op(mem, &basis, 0, g_bar_cu, 0, h_cu, 0, distf_cu);

    // fvdgv = fv . gv
    gkyl_dg_dot_product_op(&basis, fvdgv1_cu, fv1_cu, gv1_cu);
    gkyl_dg_dot_product_op(&basis, fvdgv2_cu, fv2_cu, gv2_cu);
    gkyl_dg_dot_product_op(&basis, fvdgv3_cu, fv3_cu, gv3_cu);
  
    // copy from device and check if things are ok
    gkyl_array_copy(f_bar, f_bar_cu);
    gkyl_array_copy(g_bar, g_bar_cu);
    gkyl_array_copy(g_bar, g_bar_cu);
    gkyl_array_copy(fvdgv1, fvdgv1_cu);
    gkyl_array_copy(fvdgv2, fvdgv2_cu);
    gkyl_array_copy(fvdgv3, fvdgv3_cu);
  } else {
    // allocate memory
    mem = gkyl_dg_bin_op_mem_new(f_bar->size, basis.num_basis);

    // h = f*g
    gkyl_dg_mul_op(&basis, 0, h, 0, distf, 0, distg);
    // f_bar = h/g = f
    gkyl_dg_div_op(mem, &basis, 0, f_bar, 0, h, 0, distg);
    // g_bar = h/f = g
    gkyl_dg_div_op(mem, &basis, 0, g_bar, 0, h, 0, distf);

    // fvdgv = fv . gv
    gkyl_dg_dot_product_op(&basis, fvdgv1, fv1, gv1);
    gkyl_dg_dot_product_op(&basis, fvdgv2, fv2, gv2);
    gkyl_dg_dot_product_op(&basis, fvdgv3, fv3, gv3);
  }

  for (size_t i=0; i<arr_range.volume; ++i) {
    const double *f_d = gkyl_array_cfetch(distf, i);
    const double *fbar_d = gkyl_array_cfetch(f_bar, i);
    const double *g_d = gkyl_array_cfetch(distg, i);
    const double *gbar_d = gkyl_array_cfetch(g_bar, i);
    for (int k=0; k<basis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare(f_d[k], fbar_d[k], 1e-12) );
      TEST_CHECK( gkyl_compare(g_d[k], gbar_d[k], 1e-12) );
    }
    const double *fv1_d = gkyl_array_cfetch(fv1, i);
    const double *gv1_d = gkyl_array_cfetch(gv1, i);
    const double *fvdgv1_d = gkyl_array_cfetch(fvdgv1, i);
    const double *fv2_d = gkyl_array_cfetch(fv2, i);
    const double *gv2_d = gkyl_array_cfetch(gv2, i);
    const double *fvdgv2_d = gkyl_array_cfetch(fvdgv2, i);
    const double *fv3_d = gkyl_array_cfetch(fv3, i);
    const double *gv3_d = gkyl_array_cfetch(gv3, i);
    const double *fvdgv3_d = gkyl_array_cfetch(fvdgv3, i);
    check_dot_product_1d(fv1_d, gv1_d, fvdgv1_d,
                         fv2_d, gv2_d, fvdgv2_d,
                         fv3_d, gv3_d, fvdgv3_d, poly_order);
  }

  // Test range methods
  if (use_gpu) {
    gkyl_array_clear(f_bar_cu, 0.0);
    gkyl_array_clear(g_bar_cu, 0.0);
    gkyl_array_clear(h_cu, 0.0);
    // h = f*g
    gkyl_dg_mul_op_range(&basis, 0, h_cu, 0, distf_cu, 0, distg_cu, &arr_range);
    // f_bar = h/g = f
    gkyl_dg_div_op_range(mem, &basis, 0, f_bar_cu, 0, h_cu, 0, distg_cu, &arr_range);
    // g_bar = h/f = g
    gkyl_dg_div_op_range(mem, &basis, 0, g_bar_cu, 0, h_cu, 0, distf_cu, &arr_range);

    // fvdgv = fv . gv
    gkyl_array_clear(fvdgv1_cu, 0.0);
    gkyl_array_clear(fvdgv2_cu, 0.0);
    gkyl_array_clear(fvdgv3_cu, 0.0);
    gkyl_dg_dot_product_op_range(&basis, fvdgv1_cu, fv1_cu, gv1_cu, &arr_range);
    gkyl_dg_dot_product_op_range(&basis, fvdgv2_cu, fv2_cu, gv2_cu, &arr_range);
    gkyl_dg_dot_product_op_range(&basis, fvdgv3_cu, fv3_cu, gv3_cu, &arr_range);
  
    // copy from device and check if things are ok
    gkyl_array_copy(f_bar, f_bar_cu);
    gkyl_array_copy(g_bar, g_bar_cu);
    gkyl_array_copy(fvdgv1, fvdgv1_cu);
    gkyl_array_copy(fvdgv2, fvdgv2_cu);
    gkyl_array_copy(fvdgv3, fvdgv3_cu);
  } else {
    gkyl_array_clear(f_bar, 0.0);
    gkyl_array_clear(g_bar, 0.0);
    gkyl_array_clear(h, 0.0);
    // h = f*g
    gkyl_dg_mul_op_range(&basis, 0, h, 0, distf, 0, distg, &arr_range);
    // f_bar = h/g = f
    gkyl_dg_div_op_range(mem, &basis, 0, f_bar, 0, h, 0, distg, &arr_range);
    // g_bar = h/f = g
    gkyl_dg_div_op_range(mem, &basis, 0, g_bar, 0, h, 0, distf, &arr_range);

    // fvdgv = fv . gv
    gkyl_array_clear(fvdgv1, 0.0);
    gkyl_array_clear(fvdgv2, 0.0);
    gkyl_array_clear(fvdgv3, 0.0);
    gkyl_dg_dot_product_op_range(&basis, fvdgv1, fv1, gv1, &arr_range);
    gkyl_dg_dot_product_op_range(&basis, fvdgv2, fv2, gv2, &arr_range);
    gkyl_dg_dot_product_op_range(&basis, fvdgv3, fv3, gv3, &arr_range);
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &arr_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&arr_range, iter.idx);
    const double *f_d = gkyl_array_cfetch(distf, loc);
    const double *fbar_d = gkyl_array_cfetch(f_bar, loc);
    const double *g_d = gkyl_array_cfetch(distg, loc);
    const double *gbar_d = gkyl_array_cfetch(g_bar, loc);
    for (int k=0; k<basis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare(f_d[k], fbar_d[k], 1e-12) );
      TEST_CHECK( gkyl_compare(g_d[k], gbar_d[k], 1e-12) );
    }
    const double *fv1_d = gkyl_array_cfetch(fv1, loc);
    const double *gv1_d = gkyl_array_cfetch(gv1, loc);
    const double *fvdgv1_d = gkyl_array_cfetch(fvdgv1, loc);
    const double *fv2_d = gkyl_array_cfetch(fv2, loc);
    const double *gv2_d = gkyl_array_cfetch(gv2, loc);
    const double *fvdgv2_d = gkyl_array_cfetch(fvdgv2, loc);
    const double *fv3_d = gkyl_array_cfetch(fv3, loc);
    const double *gv3_d = gkyl_array_cfetch(gv3, loc);
    const double *fvdgv3_d = gkyl_array_cfetch(fvdgv3, loc);
    check_dot_product_1d(fv1_d, gv1_d, fvdgv1_d,
                         fv2_d, gv2_d, fvdgv2_d,
                         fv3_d, gv3_d, fvdgv3_d, poly_order);
  }

  // mean ops
  double al2[2];
  if (use_gpu) {
    struct gkyl_array *mvals_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 2, arr_range.volume);
    gkyl_array_clear(mvals_cu, 0.0);

    // means are stored in h[0]
    gkyl_dg_calc_average_range(&basis, 0, mvals_cu, 0, distf_cu, arr_range);
    // L2 are stored in h[1]
    gkyl_dg_calc_l2_range(&basis, 1, mvals_cu, 0, distf_cu, arr_range);

    double* al2_cu = (double*) gkyl_cu_malloc(sizeof(double[2]));
    gkyl_array_reduce_range(al2_cu, mvals_cu, GKYL_SUM, &arr_range);
  
    gkyl_cu_memcpy(al2, al2_cu, sizeof(double[2]), GKYL_CU_MEMCPY_D2H);
    gkyl_cu_free(al2_cu);
    gkyl_array_release(mvals_cu);
  } else {
    struct gkyl_array *mvals = gkyl_array_new(GKYL_DOUBLE, 2, arr_range.volume);
    gkyl_array_clear(mvals, 0.0);

    // means are stored in h[0]
    gkyl_dg_calc_average_range(&basis, 0, mvals, 0, distf, arr_range);
    // L2 are stored in h[1]
    gkyl_dg_calc_l2_range(&basis, 1, mvals, 0, distf, arr_range);

    gkyl_array_reduce_range(al2, mvals, GKYL_SUM, &arr_range);
    gkyl_array_release(mvals);
  }


  double vol = grid.cellVolume;
  TEST_CHECK( gkyl_compare(al2[0]*vol, 2.5, 1e-14) );
  TEST_CHECK( gkyl_compare(al2[1]*vol, 19.0/3.0, 1e-14) );
  
  gkyl_proj_on_basis_release(projDistf);
  gkyl_proj_on_basis_release(projDistg);
  gkyl_proj_on_basis_release(projfv2);
  gkyl_proj_on_basis_release(projgv2);
  gkyl_proj_on_basis_release(projfv3);
  gkyl_proj_on_basis_release(projgv3);
  gkyl_array_release(distf);
  gkyl_array_release(distg);
  gkyl_array_release(f_bar);
  gkyl_array_release(g_bar);
  gkyl_array_release(fv1);
  gkyl_array_release(gv1);
  gkyl_array_release(fv2);
  gkyl_array_release(gv2);
  gkyl_array_release(fv3);
  gkyl_array_release(gv3);
  gkyl_array_release(fvdgv1);
  gkyl_array_release(fvdgv2);
  gkyl_array_release(fvdgv3);
  gkyl_dg_bin_op_mem_release(mem);
  if (use_gpu) {
    gkyl_array_release(distf_cu);
    gkyl_array_release(distg_cu);
    gkyl_array_release(f_bar_cu);
    gkyl_array_release(g_bar_cu);
    gkyl_array_release(h_cu);
    gkyl_array_release(fvdgv1_cu);
    gkyl_array_release(fvdgv2_cu);
    gkyl_array_release(fvdgv3_cu);
  } else {
    gkyl_array_release(h);
  }
}

void
test_inv_1d(int poly_order, bool use_gpu)
{
  double lower[] = {0.0}, upper[] = {1.0};
  int cells[] = {2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // Project fields.
  struct gkyl_proj_on_basis *projf = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1, f_1d, NULL);

  // Create array range.
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);  

  // Create field arrays.
  struct gkyl_array *ffld, *ffld_inv;
  struct gkyl_array *ffld_ho, *ffld_inv_ho;
  ffld_ho = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local.volume);
  ffld_inv_ho = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local.volume);
  if (use_gpu) {
    ffld = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local.volume);
    ffld_inv = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local.volume);
  } else {
    ffld = ffld_ho;
    ffld_inv = ffld_inv_ho;
  }

  // Project the field onto basis.
  gkyl_proj_on_basis_advance(projf, 0.0, &local, ffld_ho);
  gkyl_array_copy(ffld, ffld_ho);

  // Invert the field and check its results.
  gkyl_dg_inv_op(&basis, 0, ffld_inv, 0, ffld);
  gkyl_array_copy(ffld_inv_ho, ffld_inv);

  for (size_t i=0; i<local.volume; ++i) {
    const double *A = gkyl_array_cfetch(ffld_ho, i);
    const double *A_inv = gkyl_array_cfetch(ffld_inv_ho, i);

    const double A0R2 = pow(A[0],2);
    const double A1R2 = pow(A[1],2);
    double det = -0.5*(A1R2-1.0*A0R2);
  
    double A_inv_expected[basis.num_basis];
    A_inv_expected[0] = A[0]/det;
    A_inv_expected[1] = -(1.0*A[1])/det;

    for (int k=0; k<basis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare(A_inv_expected[k], A_inv[k], 1e-12) );
    }
  }

  // Test the range method.
  gkyl_array_clear(ffld_inv, 0.0);
  gkyl_dg_inv_op_range(&basis, 0, ffld_inv, 0, ffld, &local);
  gkyl_array_copy(ffld_inv_ho, ffld_inv);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&local, iter.idx);

    const double *A = gkyl_array_cfetch(ffld_ho, loc);
    const double *A_inv = gkyl_array_cfetch(ffld_inv_ho, loc);

    const double A0R2 = pow(A[0],2);
    const double A1R2 = pow(A[1],2);
    double det = -0.5*(A1R2-1.0*A0R2);
  
    double A_inv_expected[basis.num_basis];
    A_inv_expected[0] = A[0]/det;
    A_inv_expected[1] = -(1.0*A[1])/det;

    for (int k=0; k<basis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare(A_inv_expected[k], A_inv[k], 1e-12) );
    }
  }

  gkyl_proj_on_basis_release(projf);
  gkyl_array_release(ffld);
  gkyl_array_release(ffld_inv);
  if (use_gpu) {
    gkyl_array_release(ffld_ho);
    gkyl_array_release(ffld_inv_ho);
  }
}

void f_1d2d(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1 + x;
}

void f_2d(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double y = xn[1];
  fout[0] = 2 + x + y;
}
void g_2d(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double y = xn[1];
  fout[0] = 2*x*y + 8;
}

void fv2_2d(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double y = xn[1];
  fout[0] = 2 + x + y;
  fout[1] = 1 - x + y;
}
void gv2_2d(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double y = xn[1];
  fout[0] = 2*x*y + 8;
  fout[1] = x*y - 4;
}

void fv3_2d(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double y = xn[1];
  fout[0] = 2 + x + y;
  fout[1] = 1 - x + y;
  fout[2] = 0.5 + x - y;
}
void gv3_2d(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double y = xn[1];
  fout[0] = 2*x*y + 8;
  fout[1] = x*y - 4;
  fout[2] = 0.5*x*y + 2;
}

void check_dot_product_2d(const double *fv2_d, const double *gv2_d, const double *fvdgv2_d,
                          const double *fv3_d, const double *gv3_d, const double *fvdgv3_d,
                          int poly_order) {
  if (poly_order == 1) {
    TEST_CHECK( gkyl_compare(fvdgv2_d[0], 0.5*(fv2_d[7]*gv2_d[7]+fv2_d[6]*gv2_d[6]+fv2_d[5]*gv2_d[5]+fv2_d[4]*gv2_d[4]+fv2_d[3]*gv2_d[3]+fv2_d[2]*gv2_d[2]+fv2_d[1]*gv2_d[1]+fv2_d[0]*gv2_d[0]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv2_d[1], 0.5*(fv2_d[6]*gv2_d[7]+gv2_d[6]*fv2_d[7]+fv2_d[4]*gv2_d[5]+gv2_d[4]*fv2_d[5]+fv2_d[2]*gv2_d[3]+gv2_d[2]*fv2_d[3]+fv2_d[0]*gv2_d[1]+gv2_d[0]*fv2_d[1]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv2_d[2], 0.5*(fv2_d[5]*gv2_d[7]+gv2_d[5]*fv2_d[7]+fv2_d[4]*gv2_d[6]+gv2_d[4]*fv2_d[6]+fv2_d[1]*gv2_d[3]+gv2_d[1]*fv2_d[3]+fv2_d[0]*gv2_d[2]+gv2_d[0]*fv2_d[2]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv2_d[3], 0.5*(fv2_d[4]*gv2_d[7]+gv2_d[4]*fv2_d[7]+fv2_d[5]*gv2_d[6]+gv2_d[5]*fv2_d[6]+fv2_d[0]*gv2_d[3]+gv2_d[0]*fv2_d[3]+fv2_d[1]*gv2_d[2]+gv2_d[1]*fv2_d[2]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv3_d[0], 0.5*(fv3_d[11]*gv3_d[11]+fv3_d[10]*gv3_d[10]+fv3_d[9]*gv3_d[9]+fv3_d[8]*gv3_d[8]+fv3_d[7]*gv3_d[7]+fv3_d[6]*gv3_d[6]+fv3_d[5]*gv3_d[5]+fv3_d[4]*gv3_d[4]+fv3_d[3]*gv3_d[3]+fv3_d[2]*gv3_d[2]+fv3_d[1]*gv3_d[1]+fv3_d[0]*gv3_d[0]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv3_d[1], 0.5*(fv3_d[10]*gv3_d[11]+gv3_d[10]*fv3_d[11]+fv3_d[8]*gv3_d[9]+gv3_d[8]*fv3_d[9]+fv3_d[6]*gv3_d[7]+gv3_d[6]*fv3_d[7]+fv3_d[4]*gv3_d[5]+gv3_d[4]*fv3_d[5]+fv3_d[2]*gv3_d[3]+gv3_d[2]*fv3_d[3]+fv3_d[0]*gv3_d[1]+gv3_d[0]*fv3_d[1]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv3_d[2], 0.5*(fv3_d[9]*gv3_d[11]+gv3_d[9]*fv3_d[11]+fv3_d[8]*gv3_d[10]+gv3_d[8]*fv3_d[10]+fv3_d[5]*gv3_d[7]+gv3_d[5]*fv3_d[7]+fv3_d[4]*gv3_d[6]+gv3_d[4]*fv3_d[6]+fv3_d[1]*gv3_d[3]+gv3_d[1]*fv3_d[3]+fv3_d[0]*gv3_d[2]+gv3_d[0]*fv3_d[2]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv3_d[3], 0.5*(fv3_d[8]*gv3_d[11]+gv3_d[8]*fv3_d[11]+fv3_d[9]*gv3_d[10]+gv3_d[9]*fv3_d[10]+fv3_d[4]*gv3_d[7]+gv3_d[4]*fv3_d[7]+fv3_d[5]*gv3_d[6]+gv3_d[5]*fv3_d[6]+fv3_d[0]*gv3_d[3]+gv3_d[0]*fv3_d[3]+fv3_d[1]*gv3_d[2]+gv3_d[1]*fv3_d[2]), 1e-12) );
  } else if (poly_order == 2) {
    TEST_CHECK( gkyl_compare(fvdgv2_d[0], 0.5*(fv2_d[15]*gv2_d[15]+fv2_d[14]*gv2_d[14]+fv2_d[13]*gv2_d[13]+fv2_d[12]*gv2_d[12]+fv2_d[11]*gv2_d[11]+fv2_d[10]*gv2_d[10]+fv2_d[9]*gv2_d[9]+fv2_d[8]*gv2_d[8]+fv2_d[7]*gv2_d[7]+fv2_d[6]*gv2_d[6]+fv2_d[5]*gv2_d[5]+fv2_d[4]*gv2_d[4]+fv2_d[3]*gv2_d[3]+fv2_d[2]*gv2_d[2]+fv2_d[1]*gv2_d[1]+fv2_d[0]*gv2_d[0]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv2_d[1], 0.03333333333333333*(15.0*fv2_d[13]*gv2_d[15]+15.0*gv2_d[13]*fv2_d[15]+13.41640786499874*fv2_d[11]*gv2_d[14]+13.41640786499874*gv2_d[11]*fv2_d[14]+13.41640786499874*fv2_d[9]*gv2_d[12]+13.41640786499874*gv2_d[9]*fv2_d[12]+15.0*fv2_d[10]*gv2_d[11]+15.0*gv2_d[10]*fv2_d[11]+15.0*fv2_d[8]*gv2_d[9]+15.0*gv2_d[8]*fv2_d[9]+15.0*fv2_d[5]*gv2_d[7]+15.0*gv2_d[5]*fv2_d[7]+13.41640786499874*fv2_d[3]*gv2_d[6]+13.41640786499874*gv2_d[3]*fv2_d[6]+13.41640786499874*fv2_d[1]*gv2_d[4]+13.41640786499874*gv2_d[1]*fv2_d[4]+15.0*fv2_d[2]*gv2_d[3]+15.0*gv2_d[2]*fv2_d[3]+15.0*fv2_d[0]*gv2_d[1]+15.0*gv2_d[0]*fv2_d[1]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv2_d[2], 0.03333333333333333*(13.41640786499874*fv2_d[11]*gv2_d[15]+13.41640786499874*gv2_d[11]*fv2_d[15]+15.0*fv2_d[12]*gv2_d[14]+15.0*gv2_d[12]*fv2_d[14]+13.41640786499874*fv2_d[10]*gv2_d[13]+13.41640786499874*gv2_d[10]*fv2_d[13]+15.0*fv2_d[9]*gv2_d[11]+15.0*gv2_d[9]*fv2_d[11]+15.0*fv2_d[8]*gv2_d[10]+15.0*gv2_d[8]*fv2_d[10]+13.41640786499874*fv2_d[3]*gv2_d[7]+13.41640786499874*gv2_d[3]*fv2_d[7]+15.0*fv2_d[4]*gv2_d[6]+15.0*gv2_d[4]*fv2_d[6]+13.41640786499874*fv2_d[2]*gv2_d[5]+13.41640786499874*gv2_d[2]*fv2_d[5]+15.0*fv2_d[1]*gv2_d[3]+15.0*gv2_d[1]*fv2_d[3]+15.0*fv2_d[0]*gv2_d[2]+15.0*gv2_d[0]*fv2_d[2]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv2_d[3], 0.03333333333333333*((12.0*fv2_d[14]+13.41640786499874*fv2_d[10])*gv2_d[15]+(12.0*gv2_d[14]+13.41640786499874*gv2_d[10])*fv2_d[15]+13.41640786499874*fv2_d[9]*gv2_d[14]+13.41640786499874*gv2_d[9]*fv2_d[14]+13.41640786499874*fv2_d[11]*gv2_d[13]+13.41640786499874*gv2_d[11]*fv2_d[13]+13.41640786499874*fv2_d[11]*gv2_d[12]+13.41640786499874*gv2_d[11]*fv2_d[12]+15.0*fv2_d[8]*gv2_d[11]+15.0*gv2_d[8]*fv2_d[11]+15.0*fv2_d[9]*gv2_d[10]+15.0*gv2_d[9]*fv2_d[10]+(12.0*fv2_d[6]+13.41640786499874*fv2_d[2])*gv2_d[7]+(12.0*gv2_d[6]+13.41640786499874*gv2_d[2])*fv2_d[7]+13.41640786499874*fv2_d[1]*gv2_d[6]+13.41640786499874*gv2_d[1]*fv2_d[6]+13.41640786499874*fv2_d[3]*gv2_d[5]+13.41640786499874*gv2_d[3]*fv2_d[5]+13.41640786499874*fv2_d[3]*gv2_d[4]+13.41640786499874*gv2_d[3]*fv2_d[4]+15.0*fv2_d[0]*gv2_d[3]+15.0*gv2_d[0]*fv2_d[3]+15.0*fv2_d[1]*gv2_d[2]+15.0*gv2_d[1]*fv2_d[2]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv2_d[4], 0.004761904761904762*(93.91485505499116*fv2_d[15]*gv2_d[15]+(67.0820393249937*fv2_d[14]+105.0*fv2_d[10])*gv2_d[14]+105.0*gv2_d[10]*fv2_d[14]+(67.0820393249937*fv2_d[12]+105.0*fv2_d[8])*gv2_d[12]+105.0*gv2_d[8]*fv2_d[12]+93.91485505499116*fv2_d[11]*gv2_d[11]+93.91485505499116*fv2_d[9]*gv2_d[9]+93.91485505499116*fv2_d[7]*gv2_d[7]+(67.0820393249937*fv2_d[6]+105.0*fv2_d[2])*gv2_d[6]+105.0*gv2_d[2]*fv2_d[6]+(67.0820393249937*fv2_d[4]+105.0*fv2_d[0])*gv2_d[4]+105.0*gv2_d[0]*fv2_d[4]+93.91485505499116*fv2_d[3]*gv2_d[3]+93.91485505499116*fv2_d[1]*gv2_d[1]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv2_d[5], 0.004761904761904762*((67.0820393249937*fv2_d[15]+105.0*fv2_d[9])*gv2_d[15]+105.0*gv2_d[9]*fv2_d[15]+93.91485505499116*fv2_d[14]*gv2_d[14]+(67.0820393249937*fv2_d[13]+105.0*fv2_d[8])*gv2_d[13]+105.0*gv2_d[8]*fv2_d[13]+93.91485505499116*fv2_d[11]*gv2_d[11]+93.91485505499116*fv2_d[10]*gv2_d[10]+(67.0820393249937*fv2_d[7]+105.0*fv2_d[1])*gv2_d[7]+105.0*gv2_d[1]*fv2_d[7]+93.91485505499116*fv2_d[6]*gv2_d[6]+(67.0820393249937*fv2_d[5]+105.0*fv2_d[0])*gv2_d[5]+105.0*gv2_d[0]*fv2_d[5]+93.91485505499116*fv2_d[3]*gv2_d[3]+93.91485505499116*fv2_d[2]*gv2_d[2]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv2_d[6], 0.004761904761904762*(84.0*fv2_d[11]*gv2_d[15]+84.0*gv2_d[11]*fv2_d[15]+(93.91485505499116*fv2_d[13]+67.0820393249937*fv2_d[12]+105.0*fv2_d[8])*gv2_d[14]+(93.91485505499116*gv2_d[13]+67.0820393249937*gv2_d[12]+105.0*gv2_d[8])*fv2_d[14]+105.0*fv2_d[10]*gv2_d[12]+105.0*gv2_d[10]*fv2_d[12]+93.91485505499116*fv2_d[9]*gv2_d[11]+93.91485505499116*gv2_d[9]*fv2_d[11]+84.0*fv2_d[3]*gv2_d[7]+84.0*gv2_d[3]*fv2_d[7]+(93.91485505499116*fv2_d[5]+67.0820393249937*fv2_d[4]+105.0*fv2_d[0])*gv2_d[6]+(93.91485505499116*gv2_d[5]+67.0820393249937*gv2_d[4]+105.0*gv2_d[0])*fv2_d[6]+105.0*fv2_d[2]*gv2_d[4]+105.0*gv2_d[2]*fv2_d[4]+93.91485505499116*fv2_d[1]*gv2_d[3]+93.91485505499116*gv2_d[1]*fv2_d[3]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv2_d[7], 0.004761904761904762*((67.0820393249937*fv2_d[13]+93.91485505499116*fv2_d[12]+105.0*fv2_d[8])*gv2_d[15]+(67.0820393249937*gv2_d[13]+93.91485505499116*gv2_d[12]+105.0*gv2_d[8])*fv2_d[15]+84.0*fv2_d[11]*gv2_d[14]+84.0*gv2_d[11]*fv2_d[14]+105.0*fv2_d[9]*gv2_d[13]+105.0*gv2_d[9]*fv2_d[13]+93.91485505499116*fv2_d[10]*gv2_d[11]+93.91485505499116*gv2_d[10]*fv2_d[11]+(67.0820393249937*fv2_d[5]+93.91485505499116*fv2_d[4]+105.0*fv2_d[0])*gv2_d[7]+(67.0820393249937*gv2_d[5]+93.91485505499116*gv2_d[4]+105.0*gv2_d[0])*fv2_d[7]+84.0*fv2_d[3]*gv2_d[6]+84.0*gv2_d[3]*fv2_d[6]+105.0*fv2_d[1]*gv2_d[5]+105.0*gv2_d[1]*fv2_d[5]+93.91485505499116*fv2_d[2]*gv2_d[3]+93.91485505499116*gv2_d[2]*fv2_d[3]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv3_d[0], 0.5*(fv3_d[23]*gv3_d[23]+fv3_d[22]*gv3_d[22]+fv3_d[21]*gv3_d[21]+fv3_d[20]*gv3_d[20]+fv3_d[19]*gv3_d[19]+fv3_d[18]*gv3_d[18]+fv3_d[17]*gv3_d[17]+fv3_d[16]*gv3_d[16]+fv3_d[15]*gv3_d[15]+fv3_d[14]*gv3_d[14]+fv3_d[13]*gv3_d[13]+fv3_d[12]*gv3_d[12]+fv3_d[11]*gv3_d[11]+fv3_d[10]*gv3_d[10]+fv3_d[9]*gv3_d[9]+fv3_d[8]*gv3_d[8]+fv3_d[7]*gv3_d[7]+fv3_d[6]*gv3_d[6]+fv3_d[5]*gv3_d[5]+fv3_d[4]*gv3_d[4]+fv3_d[3]*gv3_d[3]+fv3_d[2]*gv3_d[2]+fv3_d[1]*gv3_d[1]+fv3_d[0]*gv3_d[0]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv3_d[1], 0.03333333333333333*(15.0*fv3_d[21]*gv3_d[23]+15.0*gv3_d[21]*fv3_d[23]+13.41640786499874*fv3_d[19]*gv3_d[22]+13.41640786499874*gv3_d[19]*fv3_d[22]+13.41640786499874*fv3_d[17]*gv3_d[20]+13.41640786499874*gv3_d[17]*fv3_d[20]+15.0*fv3_d[18]*gv3_d[19]+15.0*gv3_d[18]*fv3_d[19]+15.0*fv3_d[16]*gv3_d[17]+15.0*gv3_d[16]*fv3_d[17]+15.0*fv3_d[13]*gv3_d[15]+15.0*gv3_d[13]*fv3_d[15]+13.41640786499874*fv3_d[11]*gv3_d[14]+13.41640786499874*gv3_d[11]*fv3_d[14]+13.41640786499874*fv3_d[9]*gv3_d[12]+13.41640786499874*gv3_d[9]*fv3_d[12]+15.0*fv3_d[10]*gv3_d[11]+15.0*gv3_d[10]*fv3_d[11]+15.0*fv3_d[8]*gv3_d[9]+15.0*gv3_d[8]*fv3_d[9]+15.0*fv3_d[5]*gv3_d[7]+15.0*gv3_d[5]*fv3_d[7]+13.41640786499874*fv3_d[3]*gv3_d[6]+13.41640786499874*gv3_d[3]*fv3_d[6]+13.41640786499874*fv3_d[1]*gv3_d[4]+13.41640786499874*gv3_d[1]*fv3_d[4]+15.0*fv3_d[2]*gv3_d[3]+15.0*gv3_d[2]*fv3_d[3]+15.0*fv3_d[0]*gv3_d[1]+15.0*gv3_d[0]*fv3_d[1]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv3_d[2], 0.03333333333333333*(13.41640786499874*fv3_d[19]*gv3_d[23]+13.41640786499874*gv3_d[19]*fv3_d[23]+15.0*fv3_d[20]*gv3_d[22]+15.0*gv3_d[20]*fv3_d[22]+13.41640786499874*fv3_d[18]*gv3_d[21]+13.41640786499874*gv3_d[18]*fv3_d[21]+15.0*fv3_d[17]*gv3_d[19]+15.0*gv3_d[17]*fv3_d[19]+15.0*fv3_d[16]*gv3_d[18]+15.0*gv3_d[16]*fv3_d[18]+13.41640786499874*fv3_d[11]*gv3_d[15]+13.41640786499874*gv3_d[11]*fv3_d[15]+15.0*fv3_d[12]*gv3_d[14]+15.0*gv3_d[12]*fv3_d[14]+13.41640786499874*fv3_d[10]*gv3_d[13]+13.41640786499874*gv3_d[10]*fv3_d[13]+15.0*fv3_d[9]*gv3_d[11]+15.0*gv3_d[9]*fv3_d[11]+15.0*fv3_d[8]*gv3_d[10]+15.0*gv3_d[8]*fv3_d[10]+13.41640786499874*fv3_d[3]*gv3_d[7]+13.41640786499874*gv3_d[3]*fv3_d[7]+15.0*fv3_d[4]*gv3_d[6]+15.0*gv3_d[4]*fv3_d[6]+13.41640786499874*fv3_d[2]*gv3_d[5]+13.41640786499874*gv3_d[2]*fv3_d[5]+15.0*fv3_d[1]*gv3_d[3]+15.0*gv3_d[1]*fv3_d[3]+15.0*fv3_d[0]*gv3_d[2]+15.0*gv3_d[0]*fv3_d[2]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv3_d[3], 0.03333333333333333*((12.0*fv3_d[22]+13.41640786499874*fv3_d[18])*gv3_d[23]+(12.0*gv3_d[22]+13.41640786499874*gv3_d[18])*fv3_d[23]+13.41640786499874*fv3_d[17]*gv3_d[22]+13.41640786499874*gv3_d[17]*fv3_d[22]+13.41640786499874*fv3_d[19]*gv3_d[21]+13.41640786499874*gv3_d[19]*fv3_d[21]+13.41640786499874*fv3_d[19]*gv3_d[20]+13.41640786499874*gv3_d[19]*fv3_d[20]+15.0*fv3_d[16]*gv3_d[19]+15.0*gv3_d[16]*fv3_d[19]+15.0*fv3_d[17]*gv3_d[18]+15.0*gv3_d[17]*fv3_d[18]+(12.0*fv3_d[14]+13.41640786499874*fv3_d[10])*gv3_d[15]+(12.0*gv3_d[14]+13.41640786499874*gv3_d[10])*fv3_d[15]+13.41640786499874*fv3_d[9]*gv3_d[14]+13.41640786499874*gv3_d[9]*fv3_d[14]+13.41640786499874*fv3_d[11]*gv3_d[13]+13.41640786499874*gv3_d[11]*fv3_d[13]+13.41640786499874*fv3_d[11]*gv3_d[12]+13.41640786499874*gv3_d[11]*fv3_d[12]+15.0*fv3_d[8]*gv3_d[11]+15.0*gv3_d[8]*fv3_d[11]+15.0*fv3_d[9]*gv3_d[10]+15.0*gv3_d[9]*fv3_d[10]+(12.0*fv3_d[6]+13.41640786499874*fv3_d[2])*gv3_d[7]+(12.0*gv3_d[6]+13.41640786499874*gv3_d[2])*fv3_d[7]+13.41640786499874*fv3_d[1]*gv3_d[6]+13.41640786499874*gv3_d[1]*fv3_d[6]+13.41640786499874*fv3_d[3]*gv3_d[5]+13.41640786499874*gv3_d[3]*fv3_d[5]+13.41640786499874*fv3_d[3]*gv3_d[4]+13.41640786499874*gv3_d[3]*fv3_d[4]+15.0*fv3_d[0]*gv3_d[3]+15.0*gv3_d[0]*fv3_d[3]+15.0*fv3_d[1]*gv3_d[2]+15.0*gv3_d[1]*fv3_d[2]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv3_d[4], 0.004761904761904762*(93.91485505499116*fv3_d[23]*gv3_d[23]+(67.0820393249937*fv3_d[22]+105.0*fv3_d[18])*gv3_d[22]+105.0*gv3_d[18]*fv3_d[22]+(67.0820393249937*fv3_d[20]+105.0*fv3_d[16])*gv3_d[20]+105.0*gv3_d[16]*fv3_d[20]+93.91485505499116*fv3_d[19]*gv3_d[19]+93.91485505499116*fv3_d[17]*gv3_d[17]+93.91485505499116*fv3_d[15]*gv3_d[15]+(67.0820393249937*fv3_d[14]+105.0*fv3_d[10])*gv3_d[14]+105.0*gv3_d[10]*fv3_d[14]+(67.0820393249937*fv3_d[12]+105.0*fv3_d[8])*gv3_d[12]+105.0*gv3_d[8]*fv3_d[12]+93.91485505499116*fv3_d[11]*gv3_d[11]+93.91485505499116*fv3_d[9]*gv3_d[9]+93.91485505499116*fv3_d[7]*gv3_d[7]+(67.0820393249937*fv3_d[6]+105.0*fv3_d[2])*gv3_d[6]+105.0*gv3_d[2]*fv3_d[6]+(67.0820393249937*fv3_d[4]+105.0*fv3_d[0])*gv3_d[4]+105.0*gv3_d[0]*fv3_d[4]+93.91485505499116*fv3_d[3]*gv3_d[3]+93.91485505499116*fv3_d[1]*gv3_d[1]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv3_d[5], 0.004761904761904762*((67.0820393249937*fv3_d[23]+105.0*fv3_d[17])*gv3_d[23]+105.0*gv3_d[17]*fv3_d[23]+93.91485505499116*fv3_d[22]*gv3_d[22]+(67.0820393249937*fv3_d[21]+105.0*fv3_d[16])*gv3_d[21]+105.0*gv3_d[16]*fv3_d[21]+93.91485505499116*fv3_d[19]*gv3_d[19]+93.91485505499116*fv3_d[18]*gv3_d[18]+(67.0820393249937*fv3_d[15]+105.0*fv3_d[9])*gv3_d[15]+105.0*gv3_d[9]*fv3_d[15]+93.91485505499116*fv3_d[14]*gv3_d[14]+(67.0820393249937*fv3_d[13]+105.0*fv3_d[8])*gv3_d[13]+105.0*gv3_d[8]*fv3_d[13]+93.91485505499116*fv3_d[11]*gv3_d[11]+93.91485505499116*fv3_d[10]*gv3_d[10]+(67.0820393249937*fv3_d[7]+105.0*fv3_d[1])*gv3_d[7]+105.0*gv3_d[1]*fv3_d[7]+93.91485505499116*fv3_d[6]*gv3_d[6]+(67.0820393249937*fv3_d[5]+105.0*fv3_d[0])*gv3_d[5]+105.0*gv3_d[0]*fv3_d[5]+93.91485505499116*fv3_d[3]*gv3_d[3]+93.91485505499116*fv3_d[2]*gv3_d[2]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv3_d[6], 0.004761904761904762*(84.0*fv3_d[19]*gv3_d[23]+84.0*gv3_d[19]*fv3_d[23]+(93.91485505499116*fv3_d[21]+67.0820393249937*fv3_d[20]+105.0*fv3_d[16])*gv3_d[22]+(93.91485505499116*gv3_d[21]+67.0820393249937*gv3_d[20]+105.0*gv3_d[16])*fv3_d[22]+105.0*fv3_d[18]*gv3_d[20]+105.0*gv3_d[18]*fv3_d[20]+93.91485505499116*fv3_d[17]*gv3_d[19]+93.91485505499116*gv3_d[17]*fv3_d[19]+84.0*fv3_d[11]*gv3_d[15]+84.0*gv3_d[11]*fv3_d[15]+(93.91485505499116*fv3_d[13]+67.0820393249937*fv3_d[12]+105.0*fv3_d[8])*gv3_d[14]+(93.91485505499116*gv3_d[13]+67.0820393249937*gv3_d[12]+105.0*gv3_d[8])*fv3_d[14]+105.0*fv3_d[10]*gv3_d[12]+105.0*gv3_d[10]*fv3_d[12]+93.91485505499116*fv3_d[9]*gv3_d[11]+93.91485505499116*gv3_d[9]*fv3_d[11]+84.0*fv3_d[3]*gv3_d[7]+84.0*gv3_d[3]*fv3_d[7]+(93.91485505499116*fv3_d[5]+67.0820393249937*fv3_d[4]+105.0*fv3_d[0])*gv3_d[6]+(93.91485505499116*gv3_d[5]+67.0820393249937*gv3_d[4]+105.0*gv3_d[0])*fv3_d[6]+105.0*fv3_d[2]*gv3_d[4]+105.0*gv3_d[2]*fv3_d[4]+93.91485505499116*fv3_d[1]*gv3_d[3]+93.91485505499116*gv3_d[1]*fv3_d[3]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv3_d[7], 0.004761904761904762*((67.0820393249937*fv3_d[21]+93.91485505499116*fv3_d[20]+105.0*fv3_d[16])*gv3_d[23]+(67.0820393249937*gv3_d[21]+93.91485505499116*gv3_d[20]+105.0*gv3_d[16])*fv3_d[23]+84.0*fv3_d[19]*gv3_d[22]+84.0*gv3_d[19]*fv3_d[22]+105.0*fv3_d[17]*gv3_d[21]+105.0*gv3_d[17]*fv3_d[21]+93.91485505499116*fv3_d[18]*gv3_d[19]+93.91485505499116*gv3_d[18]*fv3_d[19]+(67.0820393249937*fv3_d[13]+93.91485505499116*fv3_d[12]+105.0*fv3_d[8])*gv3_d[15]+(67.0820393249937*gv3_d[13]+93.91485505499116*gv3_d[12]+105.0*gv3_d[8])*fv3_d[15]+84.0*fv3_d[11]*gv3_d[14]+84.0*gv3_d[11]*fv3_d[14]+105.0*fv3_d[9]*gv3_d[13]+105.0*gv3_d[9]*fv3_d[13]+93.91485505499116*fv3_d[10]*gv3_d[11]+93.91485505499116*gv3_d[10]*fv3_d[11]+(67.0820393249937*fv3_d[5]+93.91485505499116*fv3_d[4]+105.0*fv3_d[0])*gv3_d[7]+(67.0820393249937*gv3_d[5]+93.91485505499116*gv3_d[4]+105.0*gv3_d[0])*fv3_d[7]+84.0*fv3_d[3]*gv3_d[6]+84.0*gv3_d[3]*fv3_d[6]+105.0*fv3_d[1]*gv3_d[5]+105.0*gv3_d[1]*fv3_d[5]+93.91485505499116*fv3_d[2]*gv3_d[3]+93.91485505499116*gv3_d[2]*fv3_d[3]), 1e-12) );
  }
  return;
}

void
test_2d(int poly_order, bool use_gpu)
{
  double lower[] = {0.0, 0.0}, upper[] = {1.0, 1.0};
  int cells[] = {2, 2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // projection updater for dist-function
  struct gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1, f_2d, NULL);
  struct gkyl_proj_on_basis *projDistg = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1, g_2d, NULL);
  // projection updaters for vector fields.
  struct gkyl_proj_on_basis *projfv2 = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 2, fv2_2d, NULL);
  struct gkyl_proj_on_basis *projgv2 = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 2, gv2_2d, NULL);
  struct gkyl_proj_on_basis *projfv3 = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 3, fv3_2d, NULL);
  struct gkyl_proj_on_basis *projgv3 = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 3, gv3_2d, NULL);

  // create array range: no ghost-cells in velocity space
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  // create distribution function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *distg = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *distf_cu, *distg_cu;
  if (use_gpu) {
    distf_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    distg_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  }

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &arr_range, distf);
  gkyl_proj_on_basis_advance(projDistg, 0.0, &arr_range, distg);

  if (use_gpu) {
    // copy host array to device
    gkyl_array_copy(distf_cu, distf);
    gkyl_array_copy(distg_cu, distg);
  }

  // create conf-space grid:
  double clower[] = {lower[0]}, cupper[] = {upper[0]};
  int ccells[] = {cells[0]};
  int cdim = sizeof(clower)/sizeof(clower[0]);
  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, cdim, clower, cupper, ccells);
  struct gkyl_range arr_crange, arr_ext_crange;
  gkyl_create_grid_ranges(&cgrid, nghost, &arr_ext_crange, &arr_crange);
  // conf-space basis functions
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  // create a conf-space factor
  struct gkyl_array *cfield = gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, arr_crange.volume);
  // project conf-space function onto basis.
  struct gkyl_proj_on_basis *proj_cfield = gkyl_proj_on_basis_new(&cgrid, &cbasis, poly_order+1, 1, f_1d2d, NULL);
  gkyl_proj_on_basis_advance(proj_cfield, 0.0, &arr_crange, cfield);
  struct gkyl_array *cfield_cu;
  if (use_gpu) {
    cfield_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, cbasis.num_basis, arr_crange.volume);
    gkyl_array_copy(cfield_cu, cfield);
  }

  struct gkyl_array *f_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *g_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *w_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(w_bar, 0.0);
  gkyl_array_clear(f_bar, 0.0);
  gkyl_array_clear(g_bar, 0.0);

  // Vector fields for dot product.
  struct gkyl_array *fv2 = gkyl_array_new(GKYL_DOUBLE, 2*basis.num_basis, arr_range.volume);
  struct gkyl_array *gv2 = gkyl_array_new(GKYL_DOUBLE, 2*basis.num_basis, arr_range.volume);
  struct gkyl_array *fv3 = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, arr_range.volume);
  struct gkyl_array *gv3 = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, arr_range.volume);
  struct gkyl_array *fv2_cu, *gv2_cu, *fv3_cu, *gv3_cu;
  if (use_gpu) {
    fv2_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 2*basis.num_basis, arr_range.volume);
    gv2_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 2*basis.num_basis, arr_range.volume);
    fv3_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, arr_range.volume);
    gv3_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, arr_range.volume);
  }

  // project vector fields on basis.
  gkyl_proj_on_basis_advance(projfv2, 0.0, &arr_range, fv2);
  gkyl_proj_on_basis_advance(projgv2, 0.0, &arr_range, gv2);
  gkyl_proj_on_basis_advance(projfv3, 0.0, &arr_range, fv3);
  gkyl_proj_on_basis_advance(projgv3, 0.0, &arr_range, gv3);
  if (use_gpu) {
    // copy host array to device
    gkyl_array_copy(fv2_cu, fv2);
    gkyl_array_copy(gv2_cu, gv2);
    gkyl_array_copy(fv3_cu, fv3);
    gkyl_array_copy(gv3_cu, gv3);
  }

  struct gkyl_array *h, *f_bar_cu, *g_bar_cu, *w_bar_cu, *h_cu;
  struct gkyl_array *fvdgv2, *fvdgv3, *fvdgv2_cu, *fvdgv3_cu;
  if (use_gpu) {
    f_bar_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    g_bar_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    w_bar_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    gkyl_array_clear(f_bar_cu, 0.0);
    gkyl_array_clear(g_bar_cu, 0.0);
    gkyl_array_clear(w_bar_cu, 0.0);

    // Product array only needs to be initialized on GPU.
    h_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    gkyl_array_clear(h_cu, 0.0);

    // Dot product of fv and gv on the GPU.
    fvdgv2_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    fvdgv3_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    gkyl_array_clear(fvdgv2_cu, 0.0);
    gkyl_array_clear(fvdgv3_cu, 0.0);
  } else {
    h = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    gkyl_array_clear(h, 0.0);
  }
  
  // Dot product of fv and gv.
  fvdgv2 = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  fvdgv3 = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(fvdgv2, 0.0);
  gkyl_array_clear(fvdgv3, 0.0);

  gkyl_dg_bin_op_mem *mem;
  if (use_gpu) {
    // allocate memory
    mem = gkyl_dg_bin_op_mem_cu_dev_new(f_bar->size, basis.num_basis);

    // h = f*g
    gkyl_dg_mul_op(&basis, 0, h_cu, 0, distf_cu, 0, distg_cu);
    // f_bar = h/g = f
    gkyl_dg_div_op(mem, &basis, 0, f_bar_cu, 0, h_cu, 0, distg_cu);
    // g_bar = h/f = g
    gkyl_dg_div_op(mem, &basis, 0, g_bar_cu, 0, h_cu, 0, distf_cu);

    // fvdgv = fv . gv
    gkyl_dg_dot_product_op(&basis, fvdgv2_cu, fv2_cu, gv2_cu);
    gkyl_dg_dot_product_op(&basis, fvdgv3_cu, fv3_cu, gv3_cu);
  
    // copy from device and check if things are ok
    gkyl_array_copy(f_bar, f_bar_cu);
    gkyl_array_copy(g_bar, g_bar_cu);
    gkyl_array_copy(fvdgv2, fvdgv2_cu);
    gkyl_array_copy(fvdgv3, fvdgv3_cu);
  } else {
    // allocate memory
    mem = gkyl_dg_bin_op_mem_new(f_bar->size, basis.num_basis);

    // h = f*g
    gkyl_dg_mul_op(&basis, 0, h, 0, distf, 0, distg);
    // f_bar = h/g = f
    gkyl_dg_div_op(mem, &basis, 0, f_bar, 0, h, 0, distg);
    // g_bar = h/f = g
    gkyl_dg_div_op(mem, &basis, 0, g_bar, 0, h, 0, distf);

    // fvdgv = fv . gv
    gkyl_dg_dot_product_op(&basis, fvdgv2, fv2, gv2);
    gkyl_dg_dot_product_op(&basis, fvdgv3, fv3, gv3);
  }

  for (size_t i=0; i<arr_range.volume; ++i) {
    const double *f_d = gkyl_array_cfetch(distf, i);
    const double *fbar_d = gkyl_array_cfetch(f_bar, i);
    const double *g_d = gkyl_array_cfetch(distg, i);
    const double *gbar_d = gkyl_array_cfetch(g_bar, i);
    for (int k=0; k<basis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare(f_d[k], fbar_d[k], 1e-12) );
      TEST_CHECK( gkyl_compare(g_d[k], gbar_d[k], 1e-12) );
    }
    const double *fv2_d = gkyl_array_cfetch(fv2, i);
    const double *gv2_d = gkyl_array_cfetch(gv2, i);
    const double *fvdgv2_d = gkyl_array_cfetch(fvdgv2, i);
    const double *fv3_d = gkyl_array_cfetch(fv3, i);
    const double *gv3_d = gkyl_array_cfetch(gv3, i);
    const double *fvdgv3_d = gkyl_array_cfetch(fvdgv3, i);
    check_dot_product_2d(fv2_d, gv2_d, fvdgv2_d,
                         fv3_d, gv3_d, fvdgv3_d, poly_order);
  }

  // Test range methods
  if (use_gpu) {
    gkyl_array_clear(f_bar_cu, 0.0);
    gkyl_array_clear(g_bar_cu, 0.0);
    gkyl_array_clear(w_bar_cu, 0.0);
    gkyl_array_clear(h_cu, 0.0);
    // h = f*g
    gkyl_dg_mul_op_range(&basis, 0, h_cu, 0, distf_cu, 0, distg_cu, &arr_range);
    // f_bar = h/g = f
    gkyl_dg_div_op_range(mem, &basis, 0, f_bar_cu, 0, h_cu, 0, distg_cu, &arr_range);
    // g_bar = h/f = g
    gkyl_dg_div_op_range(mem, &basis, 0, g_bar_cu, 0, h_cu, 0, distf_cu, &arr_range);
    // w = cfield*f
    gkyl_dg_mul_conf_phase_op_range(&cbasis, &basis, w_bar_cu, cfield_cu, distf_cu, &arr_crange, &arr_range);

    // fvdgv = fv . gv
    gkyl_array_clear(fvdgv2_cu, 0.0);
    gkyl_array_clear(fvdgv3_cu, 0.0);
    gkyl_dg_dot_product_op_range(&basis, fvdgv2_cu, fv2_cu, gv2_cu, &arr_range);
    gkyl_dg_dot_product_op_range(&basis, fvdgv3_cu, fv3_cu, gv3_cu, &arr_range);
  
    // copy from device and check if things are ok
    gkyl_array_copy(f_bar, f_bar_cu);
    gkyl_array_copy(g_bar, g_bar_cu);
    gkyl_array_copy(w_bar, w_bar_cu);
    gkyl_array_copy(fvdgv2, fvdgv2_cu);
    gkyl_array_copy(fvdgv3, fvdgv3_cu);
  } else {
    gkyl_array_clear(f_bar, 0.0);
    gkyl_array_clear(g_bar, 0.0);
    gkyl_array_clear(w_bar, 0.0);
    gkyl_array_clear(h, 0.0);
    // h = f*g
    gkyl_dg_mul_op_range(&basis, 0, h, 0, distf, 0, distg, &arr_range);
    // f_bar = h/g = f
    gkyl_dg_div_op_range(mem, &basis, 0, f_bar, 0, h, 0, distg, &arr_range);
    // g_bar = h/f = g
    gkyl_dg_div_op_range(mem, &basis, 0, g_bar, 0, h, 0, distf, &arr_range);
    // w = cfield*f
    gkyl_dg_mul_conf_phase_op_range(&cbasis, &basis, w_bar, cfield, distf, &arr_crange, &arr_range);

    // fvdgv = fv . gv
    gkyl_array_clear(fvdgv2, 0.0);
    gkyl_array_clear(fvdgv3, 0.0);
    gkyl_dg_dot_product_op_range(&basis, fvdgv2, fv2, gv2, &arr_range);
    gkyl_dg_dot_product_op_range(&basis, fvdgv3, fv3, gv3, &arr_range);
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &arr_range);

  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&arr_range, iter.idx);
    const double *f_d = gkyl_array_cfetch(distf, loc);
    const double *fbar_d = gkyl_array_cfetch(f_bar, loc);
    const double *g_d = gkyl_array_cfetch(distg, loc);
    const double *gbar_d = gkyl_array_cfetch(g_bar, loc);
    for (int k=0; k<basis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare(f_d[k], fbar_d[k], 1e-12) );
      TEST_CHECK( gkyl_compare(g_d[k], gbar_d[k], 1e-12) );
    }

    const double *wbar_d = gkyl_array_cfetch(w_bar, loc);
    int cidx[cdim];
    for (int d=0; d<cdim; d++) cidx[d] = iter.idx[d];
    long cloc = gkyl_range_idx(&arr_crange, cidx);
    const double *cf_d = gkyl_array_cfetch(cfield, cloc);
    if (poly_order==1) {
      TEST_CHECK( gkyl_compare(0.7071067811865475*(cf_d[0]*f_d[0]+cf_d[1]*f_d[1]), wbar_d[0], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*(cf_d[0]*f_d[1]+cf_d[1]*f_d[0]), wbar_d[1], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*(cf_d[1]*f_d[3]+cf_d[0]*f_d[2]), wbar_d[2], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*(cf_d[0]*f_d[3]+cf_d[1]*f_d[2]), wbar_d[3], 1e-12) );
    } else if (poly_order==2) {
      TEST_CHECK( gkyl_compare(0.7071067811865475*(cf_d[0]*f_d[0]+cf_d[1]*f_d[1]+cf_d[2]*f_d[4]), wbar_d[0], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6324555320336759*(cf_d[1]*f_d[4]+cf_d[2]*f_d[1])
                              +0.7071067811865475*(cf_d[0]*f_d[1]+cf_d[1]*f_d[0]), wbar_d[1], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*(cf_d[0]*f_d[2]+cf_d[1]*f_d[3]+cf_d[2]*f_d[6]), wbar_d[2], 1e-12) );
      TEST_CHECK( gkyl_compare(0.632455532033676*(cf_d[1]*f_d[6]+cf_d[2]*f_d[3])
                              +0.7071067811865475*(cf_d[0]*f_d[3]+cf_d[1]*f_d[2]), wbar_d[3], 1e-12) );
      TEST_CHECK( gkyl_compare(0.4517539514526256*cf_d[2]*f_d[4]
                              +0.7071067811865475*(cf_d[0]*f_d[4]+cf_d[2]*f_d[0])
                              +0.6324555320336759*cf_d[1]*f_d[1], wbar_d[4], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*(cf_d[1]*f_d[7]+cf_d[0]*f_d[5]), wbar_d[5], 1e-12) );
      TEST_CHECK( gkyl_compare(0.4517539514526256*cf_d[2]*f_d[6]
                              +0.7071067811865475*(cf_d[0]*f_d[6]+cf_d[2]*f_d[2])
                              +0.632455532033676*cf_d[1]*f_d[3], wbar_d[6], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6324555320336759*cf_d[2]*f_d[7]
                              +0.7071067811865475*(cf_d[0]*f_d[7]+cf_d[1]*f_d[5]), wbar_d[7], 1e-12) );
    } else if (poly_order==3) {
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[3]*f_d[8]+0.7071067811865475*cf_d[2]*f_d[4]+0.7071067811865475*cf_d[1]*f_d[1]+0.7071067811865475*cf_d[0]*f_d[0], wbar_d[0], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6210590034081186*cf_d[2]*f_d[8]+0.6210590034081186*cf_d[3]*f_d[4]+0.6324555320336759*cf_d[1]*f_d[4]+0.6324555320336759*f_d[1]*cf_d[2]+0.7071067811865475*cf_d[0]*f_d[1]+0.7071067811865475*f_d[0]*cf_d[1], wbar_d[1], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865474*cf_d[3]*f_d[10]+0.7071067811865475*cf_d[2]*f_d[6]+0.7071067811865475*cf_d[1]*f_d[3]+0.7071067811865475*cf_d[0]*f_d[2], wbar_d[2], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6210590034081187*cf_d[2]*f_d[10]+0.6210590034081187*cf_d[3]*f_d[6]+0.632455532033676*cf_d[1]*f_d[6]+0.6324555320336759*cf_d[2]*f_d[3]+0.7071067811865475*cf_d[0]*f_d[3]+0.7071067811865475*cf_d[1]*f_d[2], wbar_d[3], 1e-12) );
      TEST_CHECK( gkyl_compare(0.421637021355784*cf_d[3]*f_d[8]+0.6210590034081186*cf_d[1]*f_d[8]+0.4517539514526256*cf_d[2]*f_d[4]+0.7071067811865475*cf_d[0]*f_d[4]+0.6210590034081186*f_d[1]*cf_d[3]+0.7071067811865475*f_d[0]*cf_d[2]+0.6324555320336759*cf_d[1]*f_d[1], wbar_d[4], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[1]*f_d[7]+0.7071067811865475*cf_d[0]*f_d[5], wbar_d[5], 1e-12) );
      TEST_CHECK( gkyl_compare(0.4216370213557839*cf_d[3]*f_d[10]+0.6210590034081187*cf_d[1]*f_d[10]+0.4517539514526256*cf_d[2]*f_d[6]+0.7071067811865475*cf_d[0]*f_d[6]+0.6210590034081187*cf_d[3]*f_d[3]+0.632455532033676*cf_d[1]*f_d[3]+0.7071067811865475*cf_d[2]*f_d[2], wbar_d[6], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6324555320336759*cf_d[2]*f_d[7]+0.7071067811865475*cf_d[0]*f_d[7]+0.7071067811865475*cf_d[1]*f_d[5], wbar_d[7], 1e-12) );
      TEST_CHECK( gkyl_compare(0.421637021355784*cf_d[2]*f_d[8]+0.7071067811865475*cf_d[0]*f_d[8]+0.421637021355784*cf_d[3]*f_d[4]+0.6210590034081186*cf_d[1]*f_d[4]+0.7071067811865475*f_d[0]*cf_d[3]+0.6210590034081186*f_d[1]*cf_d[2], wbar_d[8], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865474*cf_d[1]*f_d[11]+0.7071067811865475*cf_d[0]*f_d[9], wbar_d[9], 1e-12) );
      TEST_CHECK( gkyl_compare(0.421637021355784*cf_d[2]*f_d[10]+0.7071067811865475*cf_d[0]*f_d[10]+0.4216370213557839*cf_d[3]*f_d[6]+0.6210590034081187*cf_d[1]*f_d[6]+0.6210590034081187*cf_d[2]*f_d[3]+0.7071067811865474*f_d[2]*cf_d[3], wbar_d[10], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6324555320336759*cf_d[2]*f_d[11]+0.7071067811865475*cf_d[0]*f_d[11]+0.7071067811865474*cf_d[1]*f_d[9], wbar_d[11], 1e-12) );
    }

    const double *fv2_d = gkyl_array_cfetch(fv2, loc);
    const double *gv2_d = gkyl_array_cfetch(gv2, loc);
    const double *fvdgv2_d = gkyl_array_cfetch(fvdgv2, loc);
    const double *fv3_d = gkyl_array_cfetch(fv3, loc);
    const double *gv3_d = gkyl_array_cfetch(gv3, loc);
    const double *fvdgv3_d = gkyl_array_cfetch(fvdgv3, loc);
    check_dot_product_2d(fv2_d, gv2_d, fvdgv2_d,
                         fv3_d, gv3_d, fvdgv3_d, poly_order);
  }

  // mean ops
  double al2[2];
  if (use_gpu) {
    struct gkyl_array *mvals_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 2, arr_range.volume);
    gkyl_array_clear(mvals_cu, 0.0);

    // means are stored in h[0]
    gkyl_dg_calc_average_range(&basis, 0, mvals_cu, 0, distf_cu, arr_range);
    // L2 are stored in h[1]
    gkyl_dg_calc_l2_range(&basis, 1, mvals_cu, 0, distf_cu, arr_range);

    double* al2_cu = (double*) gkyl_cu_malloc(sizeof(double[2]));
    gkyl_array_reduce_range(al2_cu, mvals_cu, GKYL_SUM, &arr_range);

    gkyl_cu_memcpy(al2, al2_cu, sizeof(double[2]), GKYL_CU_MEMCPY_D2H);
    gkyl_cu_free(al2_cu);
    gkyl_array_release(mvals_cu);
  } else {
    struct gkyl_array *mvals = gkyl_array_new(GKYL_DOUBLE, 2, arr_range.volume);
    gkyl_array_clear(mvals, 0.0);

    // means are stored in h[0]
    gkyl_dg_calc_average_range(&basis, 0, mvals, 0, distf, arr_range);
    // L2 are stored in h[1]
    gkyl_dg_calc_l2_range(&basis, 1, mvals, 0, distf, arr_range);

    gkyl_array_reduce_range(al2, mvals, GKYL_SUM, &arr_range);
    gkyl_array_release(mvals);
  }

  double vol = grid.cellVolume;
  TEST_CHECK( gkyl_compare(al2[0]*vol, 3.0, 1e-14) );
  TEST_CHECK( gkyl_compare(al2[1]*vol, 55.0/6.0, 1e-14) );  
  
  gkyl_proj_on_basis_release(proj_cfield);
  gkyl_proj_on_basis_release(projDistf);
  gkyl_proj_on_basis_release(projDistg);
  gkyl_proj_on_basis_release(projfv2);
  gkyl_proj_on_basis_release(projgv2);
  gkyl_proj_on_basis_release(projfv3);
  gkyl_proj_on_basis_release(projgv3);
  gkyl_array_release(distf);
  gkyl_array_release(distg);
  gkyl_array_release(f_bar);
  gkyl_array_release(g_bar);
  gkyl_array_release(cfield);
  gkyl_array_release(w_bar);
  gkyl_array_release(fv2);
  gkyl_array_release(gv2);
  gkyl_array_release(fv3);
  gkyl_array_release(gv3);
  gkyl_array_release(fvdgv2);
  gkyl_array_release(fvdgv3);
  gkyl_dg_bin_op_mem_release(mem);  
  if (use_gpu) {
    gkyl_array_release(distf_cu);
    gkyl_array_release(distg_cu);
    gkyl_array_release(cfield_cu);
    gkyl_array_release(f_bar_cu);
    gkyl_array_release(g_bar_cu);
    gkyl_array_release(w_bar_cu);
    gkyl_array_release(h_cu);
    gkyl_array_release(fvdgv2_cu);
    gkyl_array_release(fvdgv3_cu);
  } else {
    gkyl_array_release(h);
  }
}

void
test_inv_2d(int poly_order, bool use_gpu)
{
  double lower[] = {0.0, 0.0}, upper[] = {1.0, 1.0};
  int cells[] = {2, 2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // Project fields.
  struct gkyl_proj_on_basis *projf = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1, f_2d, NULL);

  // Create array range.
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);  

  // Create field arrays.
  struct gkyl_array *ffld, *ffld_inv, *iden;
  struct gkyl_array *ffld_ho, *ffld_inv_ho, *iden_ho;
  ffld_ho = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local.volume);
  ffld_inv_ho = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local.volume);
  iden_ho = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local.volume);
  if (use_gpu) {
    ffld = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local.volume);
    ffld_inv = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local.volume);
    iden = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local.volume);
  } else {
    ffld = ffld_ho;
    ffld_inv = ffld_inv_ho;
    iden = iden_ho; 
  }

  // Project the field onto basis.
  gkyl_proj_on_basis_advance(projf, 0.0, &local, ffld_ho);
  gkyl_array_copy(ffld, ffld_ho);

  // Invert the field and check its results.
  gkyl_dg_inv_op(&basis, 0, ffld_inv, 0, ffld);
  gkyl_array_copy(ffld_inv_ho, ffld_inv);

  for (size_t i=0; i<local.volume; ++i) {
    const double *A = gkyl_array_cfetch(ffld_ho, i);
    const double *A_inv = gkyl_array_cfetch(ffld_inv_ho, i);

    double A_inv_expected[basis.num_basis];

    const double A0R2 = pow(A[0],2);
    const double A0R3 = pow(A[0],3);
    const double A0R4 = pow(A[0],4);
    const double A1R2 = pow(A[1],2);
    const double A1R3 = pow(A[1],3);
    const double A1R4 = pow(A[1],4);
    const double A2R2 = pow(A[2],2);
    const double A2R3 = pow(A[2],3);
    const double A2R4 = pow(A[2],4);
    const double A3R2 = pow(A[3],2);
    const double A3R3 = pow(A[3],3);
    const double A3R4 = pow(A[3],4);
  
    double det = 0.0625*(A3R4+((-2.0*A2R2)-2.0*A1R2-2.0*A0R2)*A3R2+8.0*A[0]*A[1]*A[2]*A[3]+A2R4+((-2.0*A1R2)-2.0*A0R2)*A2R2+A1R4-2.0*A0R2*A1R2+A0R4);
  
    A_inv_expected[0] = -(0.25*(A[0]*A3R2-2.0*A[1]*A[2]*A[3]+A[0]*A2R2+A[0]*A1R2-1.0*A0R3))/det;
    A_inv_expected[1] = -(0.25*(A[1]*A3R2-2.0*A[0]*A[2]*A[3]+A[1]*A2R2-1.0*A1R3+A0R2*A[1]))/det;
    A_inv_expected[2] = -(0.25*(A[2]*A3R2-2.0*A[0]*A[1]*A[3]-1.0*A2R3+(A1R2+A0R2)*A[2]))/det;
    A_inv_expected[3] = (0.25*(A3R3+((-1.0*A2R2)-1.0*A1R2-1.0*A0R2)*A[3]+2.0*A[0]*A[1]*A[2]))/det;

    for (int k=0; k<basis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare(A_inv_expected[k], A_inv[k], 1e-12) );
    }
  }

  // Test the range method.
  gkyl_array_clear(ffld_inv, 0.0);
  gkyl_dg_inv_op_range(&basis, 0, ffld_inv, 0, ffld, &local);
  gkyl_array_copy(ffld_inv_ho, ffld_inv);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&local, iter.idx);

    const double *A = gkyl_array_cfetch(ffld_ho, loc);
    const double *A_inv = gkyl_array_cfetch(ffld_inv_ho, loc);

    double A_inv_expected[basis.num_basis];

    const double A0R2 = pow(A[0],2);
    const double A0R3 = pow(A[0],3);
    const double A0R4 = pow(A[0],4);
    const double A1R2 = pow(A[1],2);
    const double A1R3 = pow(A[1],3);
    const double A1R4 = pow(A[1],4);
    const double A2R2 = pow(A[2],2);
    const double A2R3 = pow(A[2],3);
    const double A2R4 = pow(A[2],4);
    const double A3R2 = pow(A[3],2);
    const double A3R3 = pow(A[3],3);
    const double A3R4 = pow(A[3],4);
  
    double det = 0.0625*(A3R4+((-2.0*A2R2)-2.0*A1R2-2.0*A0R2)*A3R2+8.0*A[0]*A[1]*A[2]*A[3]+A2R4+((-2.0*A1R2)-2.0*A0R2)*A2R2+A1R4-2.0*A0R2*A1R2+A0R4);
  
    A_inv_expected[0] = -(0.25*(A[0]*A3R2-2.0*A[1]*A[2]*A[3]+A[0]*A2R2+A[0]*A1R2-1.0*A0R3))/det;
    A_inv_expected[1] = -(0.25*(A[1]*A3R2-2.0*A[0]*A[2]*A[3]+A[1]*A2R2-1.0*A1R3+A0R2*A[1]))/det;
    A_inv_expected[2] = -(0.25*(A[2]*A3R2-2.0*A[0]*A[1]*A[3]-1.0*A2R3+(A1R2+A0R2)*A[2]))/det;
    A_inv_expected[3] = (0.25*(A3R3+((-1.0*A2R2)-1.0*A1R2-1.0*A0R2)*A[3]+2.0*A[0]*A[1]*A[2]))/det;

    for (int k=0; k<basis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare(A_inv_expected[k], A_inv[k], 1e-12) );
    }
  }

  // Check if A.A_inv = 1.
  gkyl_dg_mul_op_range(&basis, 0, iden, 0, ffld, 0, ffld_inv, &local);
  gkyl_array_copy(iden_ho, iden);
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&local, iter.idx);

    const double *iden_c = gkyl_array_cfetch(iden_ho, loc);
    TEST_CHECK( gkyl_compare(pow(sqrt(2.0),ndim), iden_c[0], 1e-12) );
    for (int k=1; k<basis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare(0.0, iden_c[k], 1e-12) );
    }
  }

  gkyl_proj_on_basis_release(projf);
  gkyl_array_release(ffld);
  gkyl_array_release(ffld_inv);
  gkyl_array_release(iden);
  if (use_gpu) {
    gkyl_array_release(ffld_ho);
    gkyl_array_release(ffld_inv_ho);
    gkyl_array_release(iden_ho);
  }
}

void f_1d3d(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 5 + x;
}

void f_3d(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double y = xn[1];
  double z = xn[2];
  fout[0] = 5 + x + y + z;
}
void g_3d(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double y = xn[1];
  double z = xn[2];
  fout[0] = ((x*y*z + 8)*(x*y*z + 8) + 100);
}

void fg_3d(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double y = xn[1];
  double z = xn[2];

  double f[1], g[1];
  f_3d(t, xn, f, ctx);
  g_3d(t, xn, g, ctx);
  
  fout[0] = f[0]*g[0];
}

void fv3_3d(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double y = xn[1];
  double z = xn[2];
  fout[0] = 5 + x + y + z;
  fout[1] = 2.5 - x + y - z;
  fout[2] = 12.5 + x - y + z;
}
void gv3_3d(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double y = xn[1];
  double z = xn[2];
  fout[0] = ((x*y*z + 8)*(x*y*z + 8) + 100);
  fout[1] = ((x*y*z - 4)*(x*y*z + 4) - 50);
  fout[2] = ((x*y*z + 2)*(x*y*z - 2) + 25);
}

void check_dot_product_3d(const double *fv3_d, const double *gv3_d, const double *fvdgv3_d,
                          int poly_order) {
  if (poly_order == 1) {
    TEST_CHECK( gkyl_compare(fvdgv3_d[0], 0.3535533905932737*(fv3_d[23]*gv3_d[23]+fv3_d[22]*gv3_d[22]+fv3_d[21]*gv3_d[21]+fv3_d[20]*gv3_d[20]+fv3_d[19]*gv3_d[19]+fv3_d[18]*gv3_d[18]+fv3_d[17]*gv3_d[17]+fv3_d[16]*gv3_d[16]+fv3_d[15]*gv3_d[15]+fv3_d[14]*gv3_d[14]+fv3_d[13]*gv3_d[13]+fv3_d[12]*gv3_d[12]+fv3_d[11]*gv3_d[11]+fv3_d[10]*gv3_d[10]+fv3_d[9]*gv3_d[9]+fv3_d[8]*gv3_d[8]+fv3_d[7]*gv3_d[7]+fv3_d[6]*gv3_d[6]+fv3_d[5]*gv3_d[5]+fv3_d[4]*gv3_d[4]+fv3_d[3]*gv3_d[3]+fv3_d[2]*gv3_d[2]+fv3_d[1]*gv3_d[1]+fv3_d[0]*gv3_d[0]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv3_d[1], 0.3535533905932737*(fv3_d[22]*gv3_d[23]+gv3_d[22]*fv3_d[23]+fv3_d[19]*gv3_d[21]+gv3_d[19]*fv3_d[21]+fv3_d[18]*gv3_d[20]+gv3_d[18]*fv3_d[20]+fv3_d[16]*gv3_d[17]+gv3_d[16]*fv3_d[17]+fv3_d[14]*gv3_d[15]+gv3_d[14]*fv3_d[15]+fv3_d[11]*gv3_d[13]+gv3_d[11]*fv3_d[13]+fv3_d[10]*gv3_d[12]+gv3_d[10]*fv3_d[12]+fv3_d[8]*gv3_d[9]+gv3_d[8]*fv3_d[9]+fv3_d[6]*gv3_d[7]+gv3_d[6]*fv3_d[7]+fv3_d[3]*gv3_d[5]+gv3_d[3]*fv3_d[5]+fv3_d[2]*gv3_d[4]+gv3_d[2]*fv3_d[4]+fv3_d[0]*gv3_d[1]+gv3_d[0]*fv3_d[1]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv3_d[2], 0.3535533905932737*(fv3_d[21]*gv3_d[23]+gv3_d[21]*fv3_d[23]+fv3_d[19]*gv3_d[22]+gv3_d[19]*fv3_d[22]+fv3_d[17]*gv3_d[20]+gv3_d[17]*fv3_d[20]+fv3_d[16]*gv3_d[18]+gv3_d[16]*fv3_d[18]+fv3_d[13]*gv3_d[15]+gv3_d[13]*fv3_d[15]+fv3_d[11]*gv3_d[14]+gv3_d[11]*fv3_d[14]+fv3_d[9]*gv3_d[12]+gv3_d[9]*fv3_d[12]+fv3_d[8]*gv3_d[10]+gv3_d[8]*fv3_d[10]+fv3_d[5]*gv3_d[7]+gv3_d[5]*fv3_d[7]+fv3_d[3]*gv3_d[6]+gv3_d[3]*fv3_d[6]+fv3_d[1]*gv3_d[4]+gv3_d[1]*fv3_d[4]+fv3_d[0]*gv3_d[2]+gv3_d[0]*fv3_d[2]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv3_d[3], 0.3535533905932737*(fv3_d[20]*gv3_d[23]+gv3_d[20]*fv3_d[23]+fv3_d[18]*gv3_d[22]+gv3_d[18]*fv3_d[22]+fv3_d[17]*gv3_d[21]+gv3_d[17]*fv3_d[21]+fv3_d[16]*gv3_d[19]+gv3_d[16]*fv3_d[19]+fv3_d[12]*gv3_d[15]+gv3_d[12]*fv3_d[15]+fv3_d[10]*gv3_d[14]+gv3_d[10]*fv3_d[14]+fv3_d[9]*gv3_d[13]+gv3_d[9]*fv3_d[13]+fv3_d[8]*gv3_d[11]+gv3_d[8]*fv3_d[11]+fv3_d[4]*gv3_d[7]+gv3_d[4]*fv3_d[7]+fv3_d[2]*gv3_d[6]+gv3_d[2]*fv3_d[6]+fv3_d[1]*gv3_d[5]+gv3_d[1]*fv3_d[5]+fv3_d[0]*gv3_d[3]+gv3_d[0]*fv3_d[3]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv3_d[4], 0.3535533905932737*(fv3_d[19]*gv3_d[23]+gv3_d[19]*fv3_d[23]+fv3_d[21]*gv3_d[22]+gv3_d[21]*fv3_d[22]+fv3_d[16]*gv3_d[20]+gv3_d[16]*fv3_d[20]+fv3_d[17]*gv3_d[18]+gv3_d[17]*fv3_d[18]+fv3_d[11]*gv3_d[15]+gv3_d[11]*fv3_d[15]+fv3_d[13]*gv3_d[14]+gv3_d[13]*fv3_d[14]+fv3_d[8]*gv3_d[12]+gv3_d[8]*fv3_d[12]+fv3_d[9]*gv3_d[10]+gv3_d[9]*fv3_d[10]+fv3_d[3]*gv3_d[7]+gv3_d[3]*fv3_d[7]+fv3_d[5]*gv3_d[6]+gv3_d[5]*fv3_d[6]+fv3_d[0]*gv3_d[4]+gv3_d[0]*fv3_d[4]+fv3_d[1]*gv3_d[2]+gv3_d[1]*fv3_d[2]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv3_d[5], 0.3535533905932737*(fv3_d[18]*gv3_d[23]+gv3_d[18]*fv3_d[23]+fv3_d[20]*gv3_d[22]+gv3_d[20]*fv3_d[22]+fv3_d[16]*gv3_d[21]+gv3_d[16]*fv3_d[21]+fv3_d[17]*gv3_d[19]+gv3_d[17]*fv3_d[19]+fv3_d[10]*gv3_d[15]+gv3_d[10]*fv3_d[15]+fv3_d[12]*gv3_d[14]+gv3_d[12]*fv3_d[14]+fv3_d[8]*gv3_d[13]+gv3_d[8]*fv3_d[13]+fv3_d[9]*gv3_d[11]+gv3_d[9]*fv3_d[11]+fv3_d[2]*gv3_d[7]+gv3_d[2]*fv3_d[7]+fv3_d[4]*gv3_d[6]+gv3_d[4]*fv3_d[6]+fv3_d[0]*gv3_d[5]+gv3_d[0]*fv3_d[5]+fv3_d[1]*gv3_d[3]+gv3_d[1]*fv3_d[3]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv3_d[6], 0.3535533905932737*(fv3_d[17]*gv3_d[23]+gv3_d[17]*fv3_d[23]+fv3_d[16]*gv3_d[22]+gv3_d[16]*fv3_d[22]+fv3_d[20]*gv3_d[21]+gv3_d[20]*fv3_d[21]+fv3_d[18]*gv3_d[19]+gv3_d[18]*fv3_d[19]+fv3_d[9]*gv3_d[15]+gv3_d[9]*fv3_d[15]+fv3_d[8]*gv3_d[14]+gv3_d[8]*fv3_d[14]+fv3_d[12]*gv3_d[13]+gv3_d[12]*fv3_d[13]+fv3_d[10]*gv3_d[11]+gv3_d[10]*fv3_d[11]+fv3_d[1]*gv3_d[7]+gv3_d[1]*fv3_d[7]+fv3_d[0]*gv3_d[6]+gv3_d[0]*fv3_d[6]+fv3_d[4]*gv3_d[5]+gv3_d[4]*fv3_d[5]+fv3_d[2]*gv3_d[3]+gv3_d[2]*fv3_d[3]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv3_d[7], 0.3535533905932737*(fv3_d[16]*gv3_d[23]+gv3_d[16]*fv3_d[23]+fv3_d[17]*gv3_d[22]+gv3_d[17]*fv3_d[22]+fv3_d[18]*gv3_d[21]+gv3_d[18]*fv3_d[21]+fv3_d[19]*gv3_d[20]+gv3_d[19]*fv3_d[20]+fv3_d[8]*gv3_d[15]+gv3_d[8]*fv3_d[15]+fv3_d[9]*gv3_d[14]+gv3_d[9]*fv3_d[14]+fv3_d[10]*gv3_d[13]+gv3_d[10]*fv3_d[13]+fv3_d[11]*gv3_d[12]+gv3_d[11]*fv3_d[12]+fv3_d[0]*gv3_d[7]+gv3_d[0]*fv3_d[7]+fv3_d[1]*gv3_d[6]+gv3_d[1]*fv3_d[6]+fv3_d[2]*gv3_d[5]+gv3_d[2]*fv3_d[5]+fv3_d[3]*gv3_d[4]+gv3_d[3]*fv3_d[4]), 1e-12) );
  } else if (poly_order == 2) {
    // Not checked/implemented because the kernels are too big.
  }
  return;
}

void
test_3d(int poly_order, bool use_gpu)
{
  double lower[] = {0.0, 0.0, 0.0}, upper[] = {1.0, 1.0, 1.0};
  int cells[] = {2, 2, 2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // projection updater for dist-function
  struct gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1, f_3d, NULL);
  struct gkyl_proj_on_basis *projDistg = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1, g_3d, NULL);
  // projection updaters for vector fields.
  struct gkyl_proj_on_basis *projfv3 = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 3, fv3_2d, NULL);
  struct gkyl_proj_on_basis *projgv3 = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 3, gv3_2d, NULL);

  // create array range: no ghost-cells in velocity space
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  // create distribution function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *distg = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *distf_cu, *distg_cu;
  if (use_gpu) {
    distf_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    distg_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  }

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &arr_range, distf);
  gkyl_proj_on_basis_advance(projDistg, 0.0, &arr_range, distg);

  if (use_gpu) {
    // copy host array to device
    gkyl_array_copy(distf_cu, distf);
    gkyl_array_copy(distg_cu, distg);
  }

  // create conf-space grid:
  double clower[] = {lower[0]}, cupper[] = {upper[0]};
  int ccells[] = {cells[0]};
  int cdim = sizeof(clower)/sizeof(clower[0]);
  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, cdim, clower, cupper, ccells);
  struct gkyl_range arr_crange, arr_ext_crange;
  gkyl_create_grid_ranges(&cgrid, nghost, &arr_ext_crange, &arr_crange);
  // conf-space basis functions
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  // create a conf-space factor
  struct gkyl_array *cfield = gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, arr_crange.volume);
  // project conf-space function onto basis.
  struct gkyl_proj_on_basis *proj_cfield = gkyl_proj_on_basis_new(&cgrid, &cbasis, poly_order+1, 1, f_1d3d, NULL);
  gkyl_proj_on_basis_advance(proj_cfield, 0.0, &arr_crange, cfield);
  struct gkyl_array *cfield_cu;
  if (use_gpu) {
    cfield_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, cbasis.num_basis, arr_crange.volume);
    gkyl_array_copy(cfield_cu, cfield);
  }

  struct gkyl_array *f_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *g_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *w_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(f_bar, 0.0);
  gkyl_array_clear(g_bar, 0.0);
  gkyl_array_clear(w_bar, 0.0);

  // Vector fields for dot product.
  struct gkyl_array *fv3 = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, arr_range.volume);
  struct gkyl_array *gv3 = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, arr_range.volume);
  struct gkyl_array *fv3_cu, *gv3_cu;
  if (use_gpu) {
    fv3_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, arr_range.volume);
    gv3_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, arr_range.volume);
  }

  // project vector fields on basis.
  gkyl_proj_on_basis_advance(projfv3, 0.0, &arr_range, fv3);
  gkyl_proj_on_basis_advance(projgv3, 0.0, &arr_range, gv3);
  if (use_gpu) {
    // copy host array to device
    gkyl_array_copy(fv3_cu, fv3);
    gkyl_array_copy(gv3_cu, gv3);
  }

  struct gkyl_array *h, *f_bar_cu, *g_bar_cu, *w_bar_cu, *h_cu;
  struct gkyl_array *fvdgv3, *fvdgv3_cu;
  if (use_gpu) {
    f_bar_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    g_bar_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    w_bar_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    gkyl_array_clear(f_bar_cu, 0.0);
    gkyl_array_clear(g_bar_cu, 0.0);
    gkyl_array_clear(w_bar_cu, 0.0);

    // Product array only needs to be initialized on GPU
    h_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    gkyl_array_clear(h_cu, 0.0);

    // Dot product of fv and gv on the GPU.
    fvdgv3_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    gkyl_array_clear(fvdgv3_cu, 0.0);
  } else {
    h = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    gkyl_array_clear(h, 0.0);
  }
  
  // Dot product of fv and gv.
  fvdgv3 = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(fvdgv3, 0.0);

  gkyl_dg_bin_op_mem *mem;
  if (use_gpu) {
    // allocate memory
    mem = gkyl_dg_bin_op_mem_cu_dev_new(f_bar->size, basis.num_basis);

    // h = f*g
    gkyl_dg_mul_op(&basis, 0, h_cu, 0, distf_cu, 0, distg_cu);
    // f_bar = h/g = f
    gkyl_dg_div_op(mem, &basis, 0, f_bar_cu, 0, h_cu, 0, distg_cu);
    // g_bar = h/f = g
    gkyl_dg_div_op(mem, &basis, 0, g_bar_cu, 0, h_cu, 0, distf_cu);

    // fvdgv = fv . gv
    gkyl_dg_dot_product_op(&basis, fvdgv3_cu, fv3_cu, gv3_cu);
  
    // copy from device and check if things are ok
    gkyl_array_copy(f_bar, f_bar_cu);
    gkyl_array_copy(g_bar, g_bar_cu);
    gkyl_array_copy(fvdgv3, fvdgv3_cu);
  } else {
    // allocate memory
    mem = gkyl_dg_bin_op_mem_new(f_bar->size, basis.num_basis);

    // h = f*g
    gkyl_dg_mul_op(&basis, 0, h, 0, distf, 0, distg);
    // f_bar = h/g = f
    gkyl_dg_div_op(mem, &basis, 0, f_bar, 0, h, 0, distg);
    // g_bar = h/f = g
    gkyl_dg_div_op(mem, &basis, 0, g_bar, 0, h, 0, distf);

    // fvdgv = fv . gv
    gkyl_dg_dot_product_op(&basis, fvdgv3, fv3, gv3);
  }

  for (size_t i=0; i<arr_range.volume; ++i) {
    const double *f_d = gkyl_array_cfetch(distf, i);
    const double *fbar_d = gkyl_array_cfetch(f_bar, i);
    const double *g_d = gkyl_array_cfetch(distg, i);
    const double *gbar_d = gkyl_array_cfetch(g_bar, i);
    for (int k=0; k<basis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare(f_d[k], fbar_d[k], 1e-10) );
      TEST_CHECK( gkyl_compare(g_d[k], gbar_d[k], 1e-10) );
    }
    const double *fv3_d = gkyl_array_cfetch(fv3, i);
    const double *gv3_d = gkyl_array_cfetch(gv3, i);
    const double *fvdgv3_d = gkyl_array_cfetch(fvdgv3, i);
    check_dot_product_3d(fv3_d, gv3_d, fvdgv3_d, poly_order);
  }

  // Test range methods
  if (use_gpu) {
    gkyl_array_clear(f_bar_cu, 0.0);
    gkyl_array_clear(g_bar_cu, 0.0);
    gkyl_array_clear(w_bar_cu, 0.0);
    gkyl_array_clear(h_cu, 0.0);
    // h = f*g
    gkyl_dg_mul_op_range(&basis, 0, h_cu, 0, distf_cu, 0, distg_cu, &arr_range);
    // f_bar = h/g = f
    gkyl_dg_div_op_range(mem, &basis, 0, f_bar_cu, 0, h_cu, 0, distg_cu, &arr_range);
    // g_bar = h/f = g
    gkyl_dg_div_op_range(mem, &basis, 0, g_bar_cu, 0, h_cu, 0, distf_cu, &arr_range);
    // w = cfield*f
    gkyl_dg_mul_conf_phase_op_range(&cbasis, &basis, w_bar_cu, cfield_cu, distf_cu, &arr_crange, &arr_range);

    // fvdgv = fv . gv
    gkyl_array_clear(fvdgv3_cu, 0.0);
    gkyl_dg_dot_product_op_range(&basis, fvdgv3_cu, fv3_cu, gv3_cu, &arr_range);
  
    // copy from device and check if things are ok
    gkyl_array_copy(f_bar, f_bar_cu);
    gkyl_array_copy(g_bar, g_bar_cu);
    gkyl_array_copy(w_bar, w_bar_cu);
    gkyl_array_copy(fvdgv3, fvdgv3_cu);
  } else {
    gkyl_array_clear(f_bar, 0.0);
    gkyl_array_clear(g_bar, 0.0);
    gkyl_array_clear(w_bar, 0.0);
    gkyl_array_clear(h, 0.0);
    // h = f*g
    gkyl_dg_mul_op_range(&basis, 0, h, 0, distf, 0, distg, &arr_range);
    // f_bar = h/g = f
    gkyl_dg_div_op_range(mem, &basis, 0, f_bar, 0, h, 0, distg, &arr_range);
    // g_bar = h/f = g
    gkyl_dg_div_op_range(mem, &basis, 0, g_bar, 0, h, 0, distf, &arr_range);
    // w = cfield*f
    gkyl_dg_mul_conf_phase_op_range(&cbasis, &basis, w_bar, cfield, distf, &arr_crange, &arr_range);

    // fvdgv = fv . gv
    gkyl_array_clear(fvdgv3, 0.0);
    gkyl_dg_dot_product_op_range(&basis, fvdgv3, fv3, gv3, &arr_range);
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &arr_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&arr_range, iter.idx);
    const double *f_d = gkyl_array_cfetch(distf, loc);
    const double *fbar_d = gkyl_array_cfetch(f_bar, loc);
    const double *g_d = gkyl_array_cfetch(distg, loc);
    const double *gbar_d = gkyl_array_cfetch(g_bar, loc);
    for (int k=0; k<basis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare(f_d[k], fbar_d[k], 1e-10) );
      TEST_CHECK( gkyl_compare(g_d[k], gbar_d[k], 1e-10) );
    }

    const double *wbar_d = gkyl_array_cfetch(w_bar, loc);
    int cidx[cdim];
    for (int d=0; d<cdim; d++) cidx[d] = iter.idx[d];
    long cloc = gkyl_range_idx(&arr_crange, cidx);
    const double *cf_d = gkyl_array_cfetch(cfield, cloc);
    if (poly_order==1) {
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[1]*f_d[1]+0.7071067811865475*cf_d[0]*f_d[0], wbar_d[0], 1e-8) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[0]*f_d[1]+0.7071067811865475*f_d[0]*cf_d[1], wbar_d[1], 1e-8) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[1]*f_d[4]+0.7071067811865475*cf_d[0]*f_d[2], wbar_d[2], 1e-8) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[1]*f_d[5]+0.7071067811865475*cf_d[0]*f_d[3], wbar_d[3], 1e-8) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[0]*f_d[4]+0.7071067811865475*cf_d[1]*f_d[2], wbar_d[4], 1e-8) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[0]*f_d[5]+0.7071067811865475*cf_d[1]*f_d[3], wbar_d[5], 1e-8) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[1]*f_d[7]+0.7071067811865475*cf_d[0]*f_d[6], wbar_d[6], 1e-8) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[0]*f_d[7]+0.7071067811865475*cf_d[1]*f_d[6], wbar_d[7], 1e-8) );
    } else if (poly_order==2) {
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[2]*f_d[7]+0.7071067811865475*cf_d[1]*f_d[1]+0.7071067811865475*cf_d[0]*f_d[0], wbar_d[0], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6324555320336759*cf_d[1]*f_d[7]+0.6324555320336759*f_d[1]*cf_d[2]+0.7071067811865475*cf_d[0]*f_d[1]+0.7071067811865475*f_d[0]*cf_d[1], wbar_d[1], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[2]*f_d[11]+0.7071067811865475*cf_d[1]*f_d[4]+0.7071067811865475*cf_d[0]*f_d[2], wbar_d[2], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[2]*f_d[13]+0.7071067811865475*cf_d[1]*f_d[5]+0.7071067811865475*cf_d[0]*f_d[3], wbar_d[3], 1e-12) );
      TEST_CHECK( gkyl_compare(0.632455532033676*cf_d[1]*f_d[11]+0.6324555320336759*cf_d[2]*f_d[4]+0.7071067811865475*cf_d[0]*f_d[4]+0.7071067811865475*cf_d[1]*f_d[2], wbar_d[4], 1e-12) );
      TEST_CHECK( gkyl_compare(0.632455532033676*cf_d[1]*f_d[13]+0.6324555320336759*cf_d[2]*f_d[5]+0.7071067811865475*cf_d[0]*f_d[5]+0.7071067811865475*cf_d[1]*f_d[3], wbar_d[5], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[2]*f_d[17]+0.7071067811865475*cf_d[1]*f_d[10]+0.7071067811865475*cf_d[0]*f_d[6], wbar_d[6], 1e-12) );
      TEST_CHECK( gkyl_compare(0.4517539514526256*cf_d[2]*f_d[7]+0.7071067811865475*cf_d[0]*f_d[7]+0.7071067811865475*f_d[0]*cf_d[2]+0.6324555320336759*cf_d[1]*f_d[1], wbar_d[7], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[1]*f_d[12]+0.7071067811865475*cf_d[0]*f_d[8], wbar_d[8], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[1]*f_d[15]+0.7071067811865475*cf_d[0]*f_d[9], wbar_d[9], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6324555320336759*cf_d[1]*f_d[17]+0.6324555320336759*cf_d[2]*f_d[10]+0.7071067811865475*cf_d[0]*f_d[10]+0.7071067811865475*cf_d[1]*f_d[6], wbar_d[10], 1e-12) );
      TEST_CHECK( gkyl_compare(0.4517539514526256*cf_d[2]*f_d[11]+0.7071067811865475*cf_d[0]*f_d[11]+0.632455532033676*cf_d[1]*f_d[4]+0.7071067811865475*cf_d[2]*f_d[2], wbar_d[11], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6324555320336759*cf_d[2]*f_d[12]+0.7071067811865475*cf_d[0]*f_d[12]+0.7071067811865475*cf_d[1]*f_d[8], wbar_d[12], 1e-12) );
      TEST_CHECK( gkyl_compare(0.4517539514526256*cf_d[2]*f_d[13]+0.7071067811865475*cf_d[0]*f_d[13]+0.632455532033676*cf_d[1]*f_d[5]+0.7071067811865475*cf_d[2]*f_d[3], wbar_d[13], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[1]*f_d[18]+0.7071067811865475*cf_d[0]*f_d[14], wbar_d[14], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6324555320336759*cf_d[2]*f_d[15]+0.7071067811865475*cf_d[0]*f_d[15]+0.7071067811865475*cf_d[1]*f_d[9], wbar_d[15], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[1]*f_d[19]+0.7071067811865475*cf_d[0]*f_d[16], wbar_d[16], 1e-12) );
      TEST_CHECK( gkyl_compare(0.4517539514526256*cf_d[2]*f_d[17]+0.7071067811865475*cf_d[0]*f_d[17]+0.6324555320336759*cf_d[1]*f_d[10]+0.7071067811865475*cf_d[2]*f_d[6], wbar_d[17], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6324555320336759*cf_d[2]*f_d[18]+0.7071067811865475*cf_d[0]*f_d[18]+0.7071067811865475*cf_d[1]*f_d[14], wbar_d[18], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6324555320336759*cf_d[2]*f_d[19]+0.7071067811865475*cf_d[0]*f_d[19]+0.7071067811865475*cf_d[1]*f_d[16], wbar_d[19], 1e-12) );
    } else if (poly_order==3) {
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[3]*f_d[17]+0.7071067811865475*cf_d[2]*f_d[7]+0.7071067811865475*cf_d[1]*f_d[1]+0.7071067811865475*cf_d[0]*f_d[0], wbar_d[0], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6210590034081186*cf_d[2]*f_d[17]+0.6210590034081186*cf_d[3]*f_d[7]+0.6324555320336759*cf_d[1]*f_d[7]+0.6324555320336759*f_d[1]*cf_d[2]+0.7071067811865475*cf_d[0]*f_d[1]+0.7071067811865475*f_d[0]*cf_d[1], wbar_d[1], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865474*cf_d[3]*f_d[23]+0.7071067811865475*cf_d[2]*f_d[11]+0.7071067811865475*cf_d[1]*f_d[4]+0.7071067811865475*cf_d[0]*f_d[2], wbar_d[2], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865474*cf_d[3]*f_d[25]+0.7071067811865475*cf_d[2]*f_d[13]+0.7071067811865475*cf_d[1]*f_d[5]+0.7071067811865475*cf_d[0]*f_d[3], wbar_d[3], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6210590034081187*cf_d[2]*f_d[23]+0.6210590034081187*cf_d[3]*f_d[11]+0.632455532033676*cf_d[1]*f_d[11]+0.6324555320336759*cf_d[2]*f_d[4]+0.7071067811865475*cf_d[0]*f_d[4]+0.7071067811865475*cf_d[1]*f_d[2], wbar_d[4], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6210590034081187*cf_d[2]*f_d[25]+0.6210590034081187*cf_d[3]*f_d[13]+0.632455532033676*cf_d[1]*f_d[13]+0.6324555320336759*cf_d[2]*f_d[5]+0.7071067811865475*cf_d[0]*f_d[5]+0.7071067811865475*cf_d[1]*f_d[3], wbar_d[5], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[3]*f_d[29]+0.7071067811865475*cf_d[2]*f_d[20]+0.7071067811865475*cf_d[1]*f_d[10]+0.7071067811865475*cf_d[0]*f_d[6], wbar_d[6], 1e-12) );
      TEST_CHECK( gkyl_compare(0.421637021355784*cf_d[3]*f_d[17]+0.6210590034081186*cf_d[1]*f_d[17]+0.4517539514526256*cf_d[2]*f_d[7]+0.7071067811865475*cf_d[0]*f_d[7]+0.6210590034081186*f_d[1]*cf_d[3]+0.7071067811865475*f_d[0]*cf_d[2]+0.6324555320336759*cf_d[1]*f_d[1], wbar_d[7], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[1]*f_d[12]+0.7071067811865475*cf_d[0]*f_d[8], wbar_d[8], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[1]*f_d[15]+0.7071067811865475*cf_d[0]*f_d[9], wbar_d[9], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6210590034081186*cf_d[2]*f_d[29]+0.6210590034081186*cf_d[3]*f_d[20]+0.6324555320336759*cf_d[1]*f_d[20]+0.6324555320336759*cf_d[2]*f_d[10]+0.7071067811865475*cf_d[0]*f_d[10]+0.7071067811865475*cf_d[1]*f_d[6], wbar_d[10], 1e-12) );
      TEST_CHECK( gkyl_compare(0.4216370213557839*cf_d[3]*f_d[23]+0.6210590034081187*cf_d[1]*f_d[23]+0.4517539514526256*cf_d[2]*f_d[11]+0.7071067811865475*cf_d[0]*f_d[11]+0.6210590034081187*cf_d[3]*f_d[4]+0.632455532033676*cf_d[1]*f_d[4]+0.7071067811865475*cf_d[2]*f_d[2], wbar_d[11], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6324555320336759*cf_d[2]*f_d[12]+0.7071067811865475*cf_d[0]*f_d[12]+0.7071067811865475*cf_d[1]*f_d[8], wbar_d[12], 1e-12) );
      TEST_CHECK( gkyl_compare(0.4216370213557839*cf_d[3]*f_d[25]+0.6210590034081187*cf_d[1]*f_d[25]+0.4517539514526256*cf_d[2]*f_d[13]+0.7071067811865475*cf_d[0]*f_d[13]+0.6210590034081187*cf_d[3]*f_d[5]+0.632455532033676*cf_d[1]*f_d[5]+0.7071067811865475*cf_d[2]*f_d[3], wbar_d[13], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[1]*f_d[21]+0.7071067811865475*cf_d[0]*f_d[14], wbar_d[14], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6324555320336759*cf_d[2]*f_d[15]+0.7071067811865475*cf_d[0]*f_d[15]+0.7071067811865475*cf_d[1]*f_d[9], wbar_d[15], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[1]*f_d[22]+0.7071067811865475*cf_d[0]*f_d[16], wbar_d[16], 1e-12) );
      TEST_CHECK( gkyl_compare(0.421637021355784*cf_d[2]*f_d[17]+0.7071067811865475*cf_d[0]*f_d[17]+0.421637021355784*cf_d[3]*f_d[7]+0.6210590034081186*cf_d[1]*f_d[7]+0.7071067811865475*f_d[0]*cf_d[3]+0.6210590034081186*f_d[1]*cf_d[2], wbar_d[17], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865474*cf_d[1]*f_d[24]+0.7071067811865475*cf_d[0]*f_d[18], wbar_d[18], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865474*cf_d[1]*f_d[27]+0.7071067811865475*cf_d[0]*f_d[19], wbar_d[19], 1e-12) );
      TEST_CHECK( gkyl_compare(0.421637021355784*cf_d[3]*f_d[29]+0.6210590034081186*cf_d[1]*f_d[29]+0.4517539514526256*cf_d[2]*f_d[20]+0.7071067811865475*cf_d[0]*f_d[20]+0.6210590034081186*cf_d[3]*f_d[10]+0.6324555320336759*cf_d[1]*f_d[10]+0.7071067811865475*cf_d[2]*f_d[6], wbar_d[20], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6324555320336759*cf_d[2]*f_d[21]+0.7071067811865475*cf_d[0]*f_d[21]+0.7071067811865475*cf_d[1]*f_d[14], wbar_d[21], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6324555320336759*cf_d[2]*f_d[22]+0.7071067811865475*cf_d[0]*f_d[22]+0.7071067811865475*cf_d[1]*f_d[16], wbar_d[22], 1e-12) );
      TEST_CHECK( gkyl_compare(0.421637021355784*cf_d[2]*f_d[23]+0.7071067811865475*cf_d[0]*f_d[23]+0.4216370213557839*cf_d[3]*f_d[11]+0.6210590034081187*cf_d[1]*f_d[11]+0.6210590034081187*cf_d[2]*f_d[4]+0.7071067811865474*f_d[2]*cf_d[3], wbar_d[23], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6324555320336759*cf_d[2]*f_d[24]+0.7071067811865475*cf_d[0]*f_d[24]+0.7071067811865474*cf_d[1]*f_d[18], wbar_d[24], 1e-12) );
      TEST_CHECK( gkyl_compare(0.421637021355784*cf_d[2]*f_d[25]+0.7071067811865475*cf_d[0]*f_d[25]+0.4216370213557839*cf_d[3]*f_d[13]+0.6210590034081187*cf_d[1]*f_d[13]+0.6210590034081187*cf_d[2]*f_d[5]+0.7071067811865474*cf_d[3]*f_d[3], wbar_d[25], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865474*cf_d[1]*f_d[30]+0.7071067811865475*cf_d[0]*f_d[26], wbar_d[26], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6324555320336759*cf_d[2]*f_d[27]+0.7071067811865475*cf_d[0]*f_d[27]+0.7071067811865474*cf_d[1]*f_d[19], wbar_d[27], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865474*cf_d[1]*f_d[31]+0.7071067811865475*cf_d[0]*f_d[28], wbar_d[28], 1e-12) );
      TEST_CHECK( gkyl_compare(0.421637021355784*cf_d[2]*f_d[29]+0.7071067811865475*cf_d[0]*f_d[29]+0.421637021355784*cf_d[3]*f_d[20]+0.6210590034081186*cf_d[1]*f_d[20]+0.6210590034081186*cf_d[2]*f_d[10]+0.7071067811865475*cf_d[3]*f_d[6], wbar_d[29], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6324555320336759*cf_d[2]*f_d[30]+0.7071067811865475*cf_d[0]*f_d[30]+0.7071067811865474*cf_d[1]*f_d[26], wbar_d[30], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6324555320336759*cf_d[2]*f_d[31]+0.7071067811865475*cf_d[0]*f_d[31]+0.7071067811865474*cf_d[1]*f_d[28], wbar_d[31], 1e-12) );
    }

    const double *fv3_d = gkyl_array_cfetch(fv3, loc);
    const double *gv3_d = gkyl_array_cfetch(gv3, loc);
    const double *fvdgv3_d = gkyl_array_cfetch(fvdgv3, loc);
    check_dot_product_3d(fv3_d, gv3_d, fvdgv3_d, poly_order);
  }

  // mean ops
  double al2[2];
  if (use_gpu) {
    struct gkyl_array *mvals_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 2, arr_range.volume);
    gkyl_array_clear(mvals_cu, 0.0);

    // means are stored in h[0]
    gkyl_dg_calc_average_range(&basis, 0, mvals_cu, 0, distf_cu, arr_range);
    // L2 are stored in h[1]
    gkyl_dg_calc_l2_range(&basis, 1, mvals_cu, 0, distf_cu, arr_range);

    double* al2_cu = (double*) gkyl_cu_malloc(sizeof(double[2]));
    gkyl_array_reduce_range(al2_cu, mvals_cu, GKYL_SUM, &arr_range);

    gkyl_cu_memcpy(al2, al2_cu, sizeof(double[2]), GKYL_CU_MEMCPY_D2H);
    gkyl_cu_free(al2_cu);
    gkyl_array_release(mvals_cu);
  } else {
    struct gkyl_array *mvals = gkyl_array_new(GKYL_DOUBLE, 2, arr_range.volume);
    gkyl_array_clear(mvals, 0.0);

    // means are stored in h[0]
    gkyl_dg_calc_average_range(&basis, 0, mvals, 0, distf, arr_range);
    // L2 are stored in h[1]
    gkyl_dg_calc_l2_range(&basis, 1, mvals, 0, distf, arr_range);

    gkyl_array_reduce_range(al2, mvals, GKYL_SUM, &arr_range);
    gkyl_array_release(mvals);
  }

  double vol = grid.cellVolume;
  TEST_CHECK( gkyl_compare(al2[0]*vol, 6.5, 1e-14) );
  TEST_CHECK( gkyl_compare(al2[1]*vol, 85.0/2.0, 1e-14) );

  gkyl_proj_on_basis_release(proj_cfield);
  gkyl_proj_on_basis_release(projDistf);
  gkyl_proj_on_basis_release(projDistg);
  gkyl_proj_on_basis_release(projfv3);
  gkyl_proj_on_basis_release(projgv3);
  gkyl_array_release(distf);
  gkyl_array_release(distg);
  gkyl_array_release(f_bar);
  gkyl_array_release(g_bar);
  gkyl_array_release(cfield);
  gkyl_array_release(w_bar);
  gkyl_array_release(fv3);
  gkyl_array_release(gv3);
  gkyl_array_release(fvdgv3);
  gkyl_dg_bin_op_mem_release(mem);
  if (use_gpu) {
    gkyl_array_release(distf_cu);
    gkyl_array_release(distg_cu);
    gkyl_array_release(cfield_cu);
    gkyl_array_release(f_bar_cu);
    gkyl_array_release(g_bar_cu);
    gkyl_array_release(w_bar_cu);
    gkyl_array_release(h_cu);
    gkyl_array_release(fvdgv3_cu);
  } else {
    gkyl_array_release(h);
  }
}

void f_2d2d(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  fout[0] = 5 + x*y;
}

void f_4d(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double y = xn[1];
  double vx = xn[2];
  double vy = xn[3];
  fout[0] = 5 + (x + y)*(vx+vy);
}

void
test_4d(int poly_order, bool use_gpu)
{
  double lower[] = {0.0, 0.0, 0.0, 0.0}, upper[] = {1.0, 1.0, 1.0, 1.0};
  int cells[] = {2, 2, 2, 2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // projection updater for dist-function
  struct gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1, f_4d, NULL);

  // create array range: no ghost-cells in velocity space
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  // create distribution function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *distf_cu;
  if (use_gpu) {
    distf_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  }

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &arr_range, distf);

  if (use_gpu) {
    // copy host array to device
    gkyl_array_copy(distf_cu, distf);
  }

  // create conf-space grid:
  double clower[] = {lower[0],lower[1]}, cupper[] = {upper[0],upper[1]};
  int ccells[] = {cells[0],cells[1]};
  int cdim = sizeof(clower)/sizeof(clower[0]);
  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, cdim, clower, cupper, ccells);
  struct gkyl_range arr_crange, arr_ext_crange;
  gkyl_create_grid_ranges(&cgrid, nghost, &arr_ext_crange, &arr_crange);
  // conf-space basis functions
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  // create a conf-space factor
  struct gkyl_array *cfield = gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, arr_crange.volume);
  // project conf-space function onto basis.
  struct gkyl_proj_on_basis *proj_cfield = gkyl_proj_on_basis_new(&cgrid, &cbasis, poly_order+1, 1, f_2d2d, NULL);
  gkyl_proj_on_basis_advance(proj_cfield, 0.0, &arr_crange, cfield);
  struct gkyl_array *cfield_cu;
  if (use_gpu) {
    cfield_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, cbasis.num_basis, arr_crange.volume);
    gkyl_array_copy(cfield_cu, cfield);
  }

  struct gkyl_array *w_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(w_bar, 0.0);

  struct gkyl_array *w_bar_cu;
  if (use_gpu) {
    w_bar_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    gkyl_array_clear(w_bar_cu, 0.0);
  }

  // Test range methods
  if (use_gpu) {
    gkyl_array_clear(w_bar_cu, 0.0);
    // w = cfield*f
    gkyl_dg_mul_conf_phase_op_range(&cbasis, &basis, w_bar_cu, cfield_cu, distf_cu, &arr_crange, &arr_range);

    // copy from device and check if things are ok
    gkyl_array_copy(w_bar, w_bar_cu);
  } else {
    gkyl_array_clear(w_bar, 0.0);
    // w = cfield*f
    gkyl_dg_mul_conf_phase_op_range(&cbasis, &basis, w_bar, cfield, distf, &arr_crange, &arr_range);
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &arr_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&arr_range, iter.idx);
    const double *f_d = gkyl_array_cfetch(distf, loc);

    const double *wbar_d = gkyl_array_cfetch(w_bar, loc);
    int cidx[cdim];
    for (int d=0; d<cdim; d++) cidx[d] = iter.idx[d];
    long cloc = gkyl_range_idx(&arr_crange, cidx);
    const double *cf_d = gkyl_array_cfetch(cfield, cloc);
    if (poly_order==1) {
      TEST_CHECK( gkyl_compare(0.5*(cf_d[3]*f_d[5]+cf_d[2]*f_d[2]+cf_d[1]*f_d[1]+cf_d[0]*f_d[0]),     wbar_d[0], 1e-8) );
      TEST_CHECK( gkyl_compare(0.5*(cf_d[2]*f_d[5]+f_d[2]*cf_d[3]+cf_d[0]*f_d[1]+f_d[0]*cf_d[1]),     wbar_d[1], 1e-8) );
      TEST_CHECK( gkyl_compare(0.5*(cf_d[1]*f_d[5]+f_d[1]*cf_d[3]+cf_d[0]*f_d[2]+f_d[0]*cf_d[2]),     wbar_d[2], 1e-8) );
      TEST_CHECK( gkyl_compare(0.5*(cf_d[3]*f_d[11]+cf_d[2]*f_d[7]+cf_d[1]*f_d[6]+cf_d[0]*f_d[3]),    wbar_d[3], 1e-8) );
      TEST_CHECK( gkyl_compare(0.5*(cf_d[3]*f_d[12]+cf_d[2]*f_d[9]+cf_d[1]*f_d[8]+cf_d[0]*f_d[4]),    wbar_d[4], 1e-8) );
      TEST_CHECK( gkyl_compare(0.5*(cf_d[0]*f_d[5]+f_d[0]*cf_d[3]+cf_d[1]*f_d[2]+f_d[1]*cf_d[2]),     wbar_d[5], 1e-8) );
      TEST_CHECK( gkyl_compare(0.5*(cf_d[2]*f_d[11]+cf_d[3]*f_d[7]+cf_d[0]*f_d[6]+cf_d[1]*f_d[3]),    wbar_d[6], 1e-8) );
      TEST_CHECK( gkyl_compare(0.5*(cf_d[1]*f_d[11]+cf_d[0]*f_d[7]+cf_d[3]*f_d[6]+cf_d[2]*f_d[3]),    wbar_d[7], 1e-8) );
      TEST_CHECK( gkyl_compare(0.5*(cf_d[2]*f_d[12]+cf_d[3]*f_d[9]+cf_d[0]*f_d[8]+cf_d[1]*f_d[4]),    wbar_d[8], 1e-8) );
      TEST_CHECK( gkyl_compare(0.5*(cf_d[1]*f_d[12]+cf_d[0]*f_d[9]+cf_d[3]*f_d[8]+cf_d[2]*f_d[4]),    wbar_d[9], 1e-8) );
      TEST_CHECK( gkyl_compare(0.5*(cf_d[3]*f_d[15]+cf_d[2]*f_d[14]+cf_d[1]*f_d[13]+cf_d[0]*f_d[10]), wbar_d[10], 1e-8) );
      TEST_CHECK( gkyl_compare(0.5*(cf_d[0]*f_d[11]+cf_d[1]*f_d[7]+cf_d[2]*f_d[6]+cf_d[3]*f_d[3]),    wbar_d[11], 1e-8) );
      TEST_CHECK( gkyl_compare(0.5*(cf_d[0]*f_d[12]+cf_d[1]*f_d[9]+cf_d[2]*f_d[8]+cf_d[3]*f_d[4]),    wbar_d[12], 1e-8) );
      TEST_CHECK( gkyl_compare(0.5*(cf_d[2]*f_d[15]+cf_d[3]*f_d[14]+cf_d[0]*f_d[13]+cf_d[1]*f_d[10]), wbar_d[13], 1e-8) );
      TEST_CHECK( gkyl_compare(0.5*(cf_d[1]*f_d[15]+cf_d[0]*f_d[14]+cf_d[3]*f_d[13]+cf_d[2]*f_d[10]), wbar_d[14], 1e-8) );
      TEST_CHECK( gkyl_compare(0.5*(cf_d[0]*f_d[15]+cf_d[1]*f_d[14]+cf_d[2]*f_d[13]+cf_d[3]*f_d[10]), wbar_d[15], 1e-8) );
    } else if (poly_order==2) {
      TEST_CHECK( gkyl_compare(0.5*(cf_d[7]*f_d[20]+cf_d[6]*f_d[19]+cf_d[5]*f_d[12]+cf_d[4]*f_d[11]+cf_d[3]*f_d[5]+cf_d[2]*f_d[2]+cf_d[1]*f_d[1]+cf_d[0]*f_d[0]), wbar_d[0], 1e-8) );
      TEST_CHECK( gkyl_compare(0.03333333333333333*(15.0*cf_d[5]*f_d[20]+13.41640786499874*cf_d[3]*f_d[19]+15.0*cf_d[7]*f_d[12]+13.41640786499874*cf_d[1]*f_d[11]+13.41640786499874*f_d[5]*cf_d[6]+15.0*cf_d[2]*f_d[5]+13.41640786499874*f_d[1]*cf_d[4]+15.0*f_d[2]*cf_d[3]+15.0*cf_d[0]*f_d[1]+15.0*f_d[0]*cf_d[1]), wbar_d[1], 1e-8) );
      TEST_CHECK( gkyl_compare(0.03333333333333333*(13.41640786499874*cf_d[3]*f_d[20]+15.0*cf_d[4]*f_d[19]+13.41640786499874*cf_d[2]*f_d[12]+15.0*cf_d[6]*f_d[11]+13.41640786499874*f_d[5]*cf_d[7]+15.0*cf_d[1]*f_d[5]+13.41640786499874*f_d[2]*cf_d[5]+15.0*f_d[1]*cf_d[3]+15.0*cf_d[0]*f_d[2]+15.0*f_d[0]*cf_d[2]), wbar_d[2], 1e-8) );
      TEST_CHECK( gkyl_compare(0.03333333333333333*(15.0*cf_d[7]*f_d[33]+15.0*cf_d[6]*f_d[32]+15.0*cf_d[5]*f_d[22]+15.0*cf_d[4]*f_d[21]+15.0*cf_d[3]*f_d[15]+15.0*cf_d[2]*f_d[7]+15.0*cf_d[1]*f_d[6]+15.0*cf_d[0]*f_d[3]), wbar_d[3], 1e-8) );
      TEST_CHECK( gkyl_compare(0.03333333333333333*(15.0*cf_d[7]*f_d[36]+15.0*cf_d[6]*f_d[35]+15.0*cf_d[5]*f_d[26]+15.0*cf_d[4]*f_d[25]+15.0*cf_d[3]*f_d[16]+15.0*cf_d[2]*f_d[9]+15.0*cf_d[1]*f_d[8]+15.0*cf_d[0]*f_d[4]), wbar_d[4], 1e-8) );
      TEST_CHECK( gkyl_compare(0.03333333333333333*((12.0*cf_d[6]+13.41640786499874*cf_d[2])*f_d[20]+(12.0*cf_d[7]+13.41640786499874*cf_d[1])*f_d[19]+13.41640786499874*cf_d[3]*f_d[12]+13.41640786499874*cf_d[3]*f_d[11]+13.41640786499874*f_d[2]*cf_d[7]+13.41640786499874*f_d[1]*cf_d[6]+(13.41640786499874*cf_d[5]+13.41640786499874*cf_d[4]+15.0*cf_d[0])*f_d[5]+15.0*f_d[0]*cf_d[3]+15.0*cf_d[1]*f_d[2]+15.0*f_d[1]*cf_d[2]), wbar_d[5], 1e-8) );
      TEST_CHECK( gkyl_compare(0.03333333333333333*(15.0*cf_d[5]*f_d[33]+13.41640786499874*cf_d[3]*f_d[32]+15.0*cf_d[7]*f_d[22]+13.41640786499874*cf_d[1]*f_d[21]+(13.41640786499874*cf_d[6]+15.0*cf_d[2])*f_d[15]+15.0*cf_d[3]*f_d[7]+(13.41640786499874*cf_d[4]+15.0*cf_d[0])*f_d[6]+15.0*cf_d[1]*f_d[3]), wbar_d[6], 1e-8) );
      TEST_CHECK( gkyl_compare(0.03333333333333333*(13.41640786499874*cf_d[3]*f_d[33]+15.0*cf_d[4]*f_d[32]+13.41640786499874*cf_d[2]*f_d[22]+15.0*cf_d[6]*f_d[21]+(13.41640786499874*cf_d[7]+15.0*cf_d[1])*f_d[15]+(13.41640786499874*cf_d[5]+15.0*cf_d[0])*f_d[7]+15.0*cf_d[3]*f_d[6]+15.0*cf_d[2]*f_d[3]), wbar_d[7], 1e-8) );
      TEST_CHECK( gkyl_compare(0.03333333333333333*(15.0*cf_d[5]*f_d[36]+13.41640786499874*cf_d[3]*f_d[35]+15.0*cf_d[7]*f_d[26]+13.41640786499874*cf_d[1]*f_d[25]+(13.41640786499874*cf_d[6]+15.0*cf_d[2])*f_d[16]+15.0*cf_d[3]*f_d[9]+(13.41640786499874*cf_d[4]+15.0*cf_d[0])*f_d[8]+15.0*cf_d[1]*f_d[4]), wbar_d[8], 1e-8) );
      TEST_CHECK( gkyl_compare(0.03333333333333333*(13.41640786499874*cf_d[3]*f_d[36]+15.0*cf_d[4]*f_d[35]+13.41640786499874*cf_d[2]*f_d[26]+15.0*cf_d[6]*f_d[25]+(13.41640786499874*cf_d[7]+15.0*cf_d[1])*f_d[16]+(13.41640786499874*cf_d[5]+15.0*cf_d[0])*f_d[9]+15.0*cf_d[3]*f_d[8]+15.0*cf_d[2]*f_d[4]), wbar_d[9], 1e-8) );
      TEST_CHECK( gkyl_compare(0.5*(cf_d[7]*f_d[45]+cf_d[6]*f_d[44]+cf_d[5]*f_d[38]+cf_d[4]*f_d[37]+cf_d[3]*f_d[31]+cf_d[2]*f_d[18]+cf_d[1]*f_d[17]+cf_d[0]*f_d[10]), wbar_d[10], 1e-8) );
      TEST_CHECK( gkyl_compare(0.004761904761904762*(93.91485505499116*cf_d[7]*f_d[20]+(67.0820393249937*cf_d[6]+105.0*cf_d[2])*f_d[19]+(67.0820393249937*cf_d[4]+105.0*cf_d[0])*f_d[11]+105.0*f_d[2]*cf_d[6]+93.91485505499116*cf_d[3]*f_d[5]+105.0*f_d[0]*cf_d[4]+93.91485505499116*cf_d[1]*f_d[1]), wbar_d[11], 1e-8) );
      TEST_CHECK( gkyl_compare(0.004761904761904762*((67.0820393249937*cf_d[7]+105.0*cf_d[1])*f_d[20]+93.91485505499116*cf_d[6]*f_d[19]+(67.0820393249937*cf_d[5]+105.0*cf_d[0])*f_d[12]+105.0*f_d[1]*cf_d[7]+93.91485505499116*cf_d[3]*f_d[5]+105.0*f_d[0]*cf_d[5]+93.91485505499116*cf_d[2]*f_d[2]), wbar_d[12], 1e-8) );
      TEST_CHECK( gkyl_compare(0.03333333333333333*(15.0*cf_d[3]*f_d[34]+15.0*cf_d[2]*f_d[24]+15.0*cf_d[1]*f_d[23]+15.0*cf_d[0]*f_d[13]), wbar_d[13], 1e-8) );
      TEST_CHECK( gkyl_compare(0.03333333333333333*(15.0*cf_d[3]*f_d[41]+15.0*cf_d[2]*f_d[29]+15.0*cf_d[1]*f_d[28]+15.0*cf_d[0]*f_d[14]), wbar_d[14], 1e-8) );
      TEST_CHECK( gkyl_compare(0.006666666666666667*((60.00000000000001*cf_d[6]+67.0820393249937*cf_d[2])*f_d[33]+(60.00000000000001*cf_d[7]+67.0820393249937*cf_d[1])*f_d[32]+67.08203932499369*cf_d[3]*f_d[22]+67.08203932499369*cf_d[3]*f_d[21]+(67.0820393249937*cf_d[5]+67.0820393249937*cf_d[4]+75.0*cf_d[0])*f_d[15]+(67.08203932499369*cf_d[7]+75.0*cf_d[1])*f_d[7]+(67.08203932499369*cf_d[6]+75.0*cf_d[2])*f_d[6]+75.0*cf_d[3]*f_d[3]), wbar_d[15], 1e-8) );
      TEST_CHECK( gkyl_compare(0.006666666666666667*((60.00000000000001*cf_d[6]+67.0820393249937*cf_d[2])*f_d[36]+(60.00000000000001*cf_d[7]+67.0820393249937*cf_d[1])*f_d[35]+67.08203932499369*cf_d[3]*f_d[26]+67.08203932499369*cf_d[3]*f_d[25]+(67.0820393249937*cf_d[5]+67.0820393249937*cf_d[4]+75.0*cf_d[0])*f_d[16]+(67.08203932499369*cf_d[7]+75.0*cf_d[1])*f_d[9]+(67.08203932499369*cf_d[6]+75.0*cf_d[2])*f_d[8]+75.0*cf_d[3]*f_d[4]), wbar_d[16], 1e-8) );
      TEST_CHECK( gkyl_compare(0.03333333333333333*(15.0*cf_d[5]*f_d[45]+13.41640786499874*cf_d[3]*f_d[44]+15.0*cf_d[7]*f_d[38]+13.41640786499874*cf_d[1]*f_d[37]+(13.41640786499874*cf_d[6]+15.0*cf_d[2])*f_d[31]+15.0*cf_d[3]*f_d[18]+(13.41640786499874*cf_d[4]+15.0*cf_d[0])*f_d[17]+15.0*cf_d[1]*f_d[10]), wbar_d[17], 1e-8) );
      TEST_CHECK( gkyl_compare(0.03333333333333333*(13.41640786499874*cf_d[3]*f_d[45]+15.0*cf_d[4]*f_d[44]+13.41640786499874*cf_d[2]*f_d[38]+15.0*cf_d[6]*f_d[37]+(13.41640786499874*cf_d[7]+15.0*cf_d[1])*f_d[31]+(13.41640786499874*cf_d[5]+15.0*cf_d[0])*f_d[18]+15.0*cf_d[3]*f_d[17]+15.0*cf_d[2]*f_d[10]), wbar_d[18], 1e-8) );
      TEST_CHECK( gkyl_compare(0.004761904761904762*(84.0*cf_d[3]*f_d[20]+(93.91485505499116*cf_d[5]+67.0820393249937*cf_d[4]+105.0*cf_d[0])*f_d[19]+93.91485505499116*cf_d[6]*f_d[12]+(67.0820393249937*cf_d[6]+105.0*cf_d[2])*f_d[11]+84.0*f_d[5]*cf_d[7]+105.0*f_d[0]*cf_d[6]+93.91485505499116*cf_d[1]*f_d[5]+105.0*f_d[2]*cf_d[4]+93.91485505499116*f_d[1]*cf_d[3]), wbar_d[19], 1e-8) );
      TEST_CHECK( gkyl_compare(0.004761904761904762*((67.0820393249937*cf_d[5]+93.91485505499116*cf_d[4]+105.0*cf_d[0])*f_d[20]+84.0*cf_d[3]*f_d[19]+(67.0820393249937*cf_d[7]+105.0*cf_d[1])*f_d[12]+93.91485505499116*cf_d[7]*f_d[11]+105.0*f_d[0]*cf_d[7]+84.0*f_d[5]*cf_d[6]+93.91485505499116*cf_d[2]*f_d[5]+105.0*f_d[1]*cf_d[5]+93.91485505499116*f_d[2]*cf_d[3]), wbar_d[20], 1e-8) );
      TEST_CHECK( gkyl_compare(0.004761904761904762*(93.91485505499116*cf_d[7]*f_d[33]+(67.0820393249937*cf_d[6]+105.0*cf_d[2])*f_d[32]+(67.0820393249937*cf_d[4]+105.0*cf_d[0])*f_d[21]+93.91485505499116*cf_d[3]*f_d[15]+105.0*cf_d[6]*f_d[7]+93.91485505499116*cf_d[1]*f_d[6]+105.0*f_d[3]*cf_d[4]), wbar_d[21], 1e-8) );
      TEST_CHECK( gkyl_compare(0.004761904761904762*((67.0820393249937*cf_d[7]+105.0*cf_d[1])*f_d[33]+93.91485505499116*cf_d[6]*f_d[32]+(67.0820393249937*cf_d[5]+105.0*cf_d[0])*f_d[22]+93.91485505499116*cf_d[3]*f_d[15]+93.91485505499116*cf_d[2]*f_d[7]+105.0*f_d[6]*cf_d[7]+105.0*f_d[3]*cf_d[5]), wbar_d[22], 1e-8) );
      TEST_CHECK( gkyl_compare(0.03333333333333333*((13.41640786499874*cf_d[6]+15.0*cf_d[2])*f_d[34]+15.0*cf_d[3]*f_d[24]+(13.41640786499874*cf_d[4]+15.0*cf_d[0])*f_d[23]+15.0*cf_d[1]*f_d[13]), wbar_d[23], 1e-8) );
      TEST_CHECK( gkyl_compare(0.03333333333333333*((13.41640786499874*cf_d[7]+15.0*cf_d[1])*f_d[34]+(13.41640786499874*cf_d[5]+15.0*cf_d[0])*f_d[24]+15.0*cf_d[3]*f_d[23]+15.0*cf_d[2]*f_d[13]), wbar_d[24], 1e-8) );
      TEST_CHECK( gkyl_compare(0.004761904761904762*(93.91485505499116*cf_d[7]*f_d[36]+(67.0820393249937*cf_d[6]+105.0*cf_d[2])*f_d[35]+(67.0820393249937*cf_d[4]+105.0*cf_d[0])*f_d[25]+93.91485505499116*cf_d[3]*f_d[16]+105.0*cf_d[6]*f_d[9]+93.91485505499116*cf_d[1]*f_d[8]+105.0*cf_d[4]*f_d[4]), wbar_d[25], 1e-8) );
      TEST_CHECK( gkyl_compare(0.004761904761904762*((67.0820393249937*cf_d[7]+105.0*cf_d[1])*f_d[36]+93.91485505499116*cf_d[6]*f_d[35]+(67.0820393249937*cf_d[5]+105.0*cf_d[0])*f_d[26]+93.91485505499116*cf_d[3]*f_d[16]+93.91485505499116*cf_d[2]*f_d[9]+105.0*cf_d[7]*f_d[8]+105.0*f_d[4]*cf_d[5]), wbar_d[26], 1e-8) );
      TEST_CHECK( gkyl_compare(0.03333333333333333*(15.0*cf_d[3]*f_d[46]+15.0*cf_d[2]*f_d[40]+15.0*cf_d[1]*f_d[39]+15.0*cf_d[0]*f_d[27]), wbar_d[27], 1e-8) );
      TEST_CHECK( gkyl_compare(0.03333333333333333*((13.41640786499874*cf_d[6]+15.0*cf_d[2])*f_d[41]+15.0*cf_d[3]*f_d[29]+(13.41640786499874*cf_d[4]+15.0*cf_d[0])*f_d[28]+15.0*cf_d[1]*f_d[14]), wbar_d[28], 1e-8) );
      TEST_CHECK( gkyl_compare(0.03333333333333333*((13.41640786499874*cf_d[7]+15.0*cf_d[1])*f_d[41]+(13.41640786499874*cf_d[5]+15.0*cf_d[0])*f_d[29]+15.0*cf_d[3]*f_d[28]+15.0*cf_d[2]*f_d[14]), wbar_d[29], 1e-8) );
      TEST_CHECK( gkyl_compare(0.03333333333333333*(15.0*cf_d[3]*f_d[47]+15.0*cf_d[2]*f_d[43]+15.0*cf_d[1]*f_d[42]+15.0*cf_d[0]*f_d[30]), wbar_d[30], 1e-8) );
      TEST_CHECK( gkyl_compare(0.03333333333333333*((12.0*cf_d[6]+13.41640786499874*cf_d[2])*f_d[45]+(12.0*cf_d[7]+13.41640786499874*cf_d[1])*f_d[44]+13.41640786499874*cf_d[3]*f_d[38]+13.41640786499874*cf_d[3]*f_d[37]+(13.41640786499874*cf_d[5]+13.41640786499874*cf_d[4]+15.0*cf_d[0])*f_d[31]+(13.41640786499874*cf_d[7]+15.0*cf_d[1])*f_d[18]+(13.41640786499874*cf_d[6]+15.0*cf_d[2])*f_d[17]+15.0*cf_d[3]*f_d[10]), wbar_d[31], 1e-8) );
      TEST_CHECK( gkyl_compare(9.523809523809524e-4*(420.0*cf_d[3]*f_d[33]+(469.5742752749559*cf_d[5]+335.4101966249685*cf_d[4]+525.0*cf_d[0])*f_d[32]+469.5742752749559*cf_d[6]*f_d[22]+(335.4101966249685*cf_d[6]+525.0000000000001*cf_d[2])*f_d[21]+(420.0000000000001*cf_d[7]+469.5742752749559*cf_d[1])*f_d[15]+525.0*cf_d[4]*f_d[7]+469.5742752749559*cf_d[3]*f_d[6]+525.0000000000001*f_d[3]*cf_d[6]), wbar_d[32], 1e-8) );
      TEST_CHECK( gkyl_compare(9.523809523809524e-4*((335.4101966249685*cf_d[5]+469.5742752749559*cf_d[4]+525.0*cf_d[0])*f_d[33]+420.0*cf_d[3]*f_d[32]+(335.4101966249685*cf_d[7]+525.0000000000001*cf_d[1])*f_d[22]+469.5742752749559*cf_d[7]*f_d[21]+(420.0000000000001*cf_d[6]+469.5742752749559*cf_d[2])*f_d[15]+469.5742752749559*cf_d[3]*f_d[7]+525.0000000000001*f_d[3]*cf_d[7]+525.0*cf_d[5]*f_d[6]), wbar_d[33], 1e-8) );
      TEST_CHECK( gkyl_compare(0.03333333333333333*((13.41640786499874*cf_d[5]+13.41640786499874*cf_d[4]+15.0*cf_d[0])*f_d[34]+(13.41640786499874*cf_d[7]+15.0*cf_d[1])*f_d[24]+(13.41640786499874*cf_d[6]+15.0*cf_d[2])*f_d[23]+15.0*cf_d[3]*f_d[13]), wbar_d[34], 1e-8) );
      TEST_CHECK( gkyl_compare(9.523809523809524e-4*(420.0*cf_d[3]*f_d[36]+(469.5742752749559*cf_d[5]+335.4101966249685*cf_d[4]+525.0*cf_d[0])*f_d[35]+469.5742752749559*cf_d[6]*f_d[26]+(335.4101966249685*cf_d[6]+525.0000000000001*cf_d[2])*f_d[25]+(420.0000000000001*cf_d[7]+469.5742752749559*cf_d[1])*f_d[16]+525.0*cf_d[4]*f_d[9]+469.5742752749559*cf_d[3]*f_d[8]+525.0000000000001*f_d[4]*cf_d[6]), wbar_d[35], 1e-8) );
      TEST_CHECK( gkyl_compare(9.523809523809524e-4*((335.4101966249685*cf_d[5]+469.5742752749559*cf_d[4]+525.0*cf_d[0])*f_d[36]+420.0*cf_d[3]*f_d[35]+(335.4101966249685*cf_d[7]+525.0000000000001*cf_d[1])*f_d[26]+469.5742752749559*cf_d[7]*f_d[25]+(420.0000000000001*cf_d[6]+469.5742752749559*cf_d[2])*f_d[16]+469.5742752749559*cf_d[3]*f_d[9]+525.0*cf_d[5]*f_d[8]+525.0000000000001*f_d[4]*cf_d[7]), wbar_d[36], 1e-8) );
      TEST_CHECK( gkyl_compare(0.004761904761904762*(93.91485505499116*cf_d[7]*f_d[45]+(67.0820393249937*cf_d[6]+105.0*cf_d[2])*f_d[44]+(67.0820393249937*cf_d[4]+105.0*cf_d[0])*f_d[37]+93.91485505499116*cf_d[3]*f_d[31]+105.0*cf_d[6]*f_d[18]+93.91485505499116*cf_d[1]*f_d[17]+105.0*cf_d[4]*f_d[10]), wbar_d[37], 1e-8) );
      TEST_CHECK( gkyl_compare(0.004761904761904762*((67.0820393249937*cf_d[7]+105.0*cf_d[1])*f_d[45]+93.91485505499116*cf_d[6]*f_d[44]+(67.0820393249937*cf_d[5]+105.0*cf_d[0])*f_d[38]+93.91485505499116*cf_d[3]*f_d[31]+93.91485505499116*cf_d[2]*f_d[18]+105.0*cf_d[7]*f_d[17]+105.0*cf_d[5]*f_d[10]), wbar_d[38], 1e-8) );
      TEST_CHECK( gkyl_compare(0.03333333333333333*((13.41640786499874*cf_d[6]+15.0*cf_d[2])*f_d[46]+15.0*cf_d[3]*f_d[40]+(13.41640786499874*cf_d[4]+15.0*cf_d[0])*f_d[39]+15.0*cf_d[1]*f_d[27]), wbar_d[39], 1e-8) );
      TEST_CHECK( gkyl_compare(0.03333333333333333*((13.41640786499874*cf_d[7]+15.0*cf_d[1])*f_d[46]+(13.41640786499874*cf_d[5]+15.0*cf_d[0])*f_d[40]+15.0*cf_d[3]*f_d[39]+15.0*cf_d[2]*f_d[27]), wbar_d[40], 1e-8) );
      TEST_CHECK( gkyl_compare(0.03333333333333333*((13.41640786499874*cf_d[5]+13.41640786499874*cf_d[4]+15.0*cf_d[0])*f_d[41]+(13.41640786499874*cf_d[7]+15.0*cf_d[1])*f_d[29]+(13.41640786499874*cf_d[6]+15.0*cf_d[2])*f_d[28]+15.0*cf_d[3]*f_d[14]), wbar_d[41], 1e-8) );
      TEST_CHECK( gkyl_compare(0.03333333333333333*((13.41640786499874*cf_d[6]+15.0*cf_d[2])*f_d[47]+15.0*cf_d[3]*f_d[43]+(13.41640786499874*cf_d[4]+15.0*cf_d[0])*f_d[42]+15.0*cf_d[1]*f_d[30]), wbar_d[42], 1e-8) );
      TEST_CHECK( gkyl_compare(0.03333333333333333*((13.41640786499874*cf_d[7]+15.0*cf_d[1])*f_d[47]+(13.41640786499874*cf_d[5]+15.0*cf_d[0])*f_d[43]+15.0*cf_d[3]*f_d[42]+15.0*cf_d[2]*f_d[30]), wbar_d[43], 1e-8) );
      TEST_CHECK( gkyl_compare(0.004761904761904762*(84.0*cf_d[3]*f_d[45]+(93.91485505499116*cf_d[5]+67.0820393249937*cf_d[4]+105.0*cf_d[0])*f_d[44]+93.91485505499116*cf_d[6]*f_d[38]+(67.0820393249937*cf_d[6]+105.0*cf_d[2])*f_d[37]+(84.0*cf_d[7]+93.91485505499116*cf_d[1])*f_d[31]+105.0*cf_d[4]*f_d[18]+93.91485505499116*cf_d[3]*f_d[17]+105.0*cf_d[6]*f_d[10]), wbar_d[44], 1e-8) );
      TEST_CHECK( gkyl_compare(0.004761904761904762*((67.0820393249937*cf_d[5]+93.91485505499116*cf_d[4]+105.0*cf_d[0])*f_d[45]+84.0*cf_d[3]*f_d[44]+(67.0820393249937*cf_d[7]+105.0*cf_d[1])*f_d[38]+93.91485505499116*cf_d[7]*f_d[37]+(84.0*cf_d[6]+93.91485505499116*cf_d[2])*f_d[31]+93.91485505499116*cf_d[3]*f_d[18]+105.0*cf_d[5]*f_d[17]+105.0*cf_d[7]*f_d[10]), wbar_d[45], 1e-8) );
      TEST_CHECK( gkyl_compare(0.03333333333333333*((13.41640786499874*cf_d[5]+13.41640786499874*cf_d[4]+15.0*cf_d[0])*f_d[46]+(13.41640786499874*cf_d[7]+15.0*cf_d[1])*f_d[40]+(13.41640786499874*cf_d[6]+15.0*cf_d[2])*f_d[39]+15.0*cf_d[3]*f_d[27]), wbar_d[46], 1e-8) );
      TEST_CHECK( gkyl_compare(0.03333333333333333*((13.41640786499874*cf_d[5]+13.41640786499874*cf_d[4]+15.0*cf_d[0])*f_d[47]+(13.41640786499874*cf_d[7]+15.0*cf_d[1])*f_d[43]+(13.41640786499874*cf_d[6]+15.0*cf_d[2])*f_d[42]+15.0*cf_d[3]*f_d[30]), wbar_d[47], 1e-8) );
    }
  }

  gkyl_proj_on_basis_release(proj_cfield);
  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
  gkyl_array_release(cfield);
  gkyl_array_release(w_bar);
  if (use_gpu) {
    gkyl_array_release(distf_cu);
    gkyl_array_release(cfield_cu);
    gkyl_array_release(w_bar_cu);
  }
}

void test_bin_ops_1d_p1_ho(){ test_1d(1, false); }
void test_bin_ops_1d_p2_ho(){ test_1d(2, false); }
void test_bin_ops_1d_p3_ho(){ test_1d(3, false); }

void test_bin_ops_inv_1d_p1_ho(){ test_inv_1d(1, false); }

void test_bin_ops_2d_p1_ho(){ test_2d(1, false); }
void test_bin_ops_2d_p2_ho(){ test_2d(2, false); }
void test_bin_ops_2d_p3_ho(){ test_2d(3, false); }

void test_bin_ops_inv_2d_p1_ho(){ test_inv_2d(1, false); }

void test_bin_ops_3d_p1_ho(){ test_3d(1, false); }
void test_bin_ops_3d_p2_ho(){ test_3d(2, false); }

void test_bin_ops_4d_p1_ho(){ test_4d(1, false); }
void test_bin_ops_4d_p2_ho(){ test_4d(2, false); }

void f_3d_p3(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double y = xn[1];
  double z = xn[2];
  fout[0] = ((5 + x)*(5 + y)*(5 + z)*(5 + x)*(5 + y)*(5 + z)*(5 + x)*(5 + y)*(5 + z) + 100);
}

void g_3d_p3(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double y = xn[1];
  double z = xn[2];
  fout[0] = ((x*y*z + 8)*(x*y*z + 8)*(x*y*z + 8) + 100);
}

void
test_bin_ops_3d_p3_ho()
{
  int poly_order = 3;
  double lower[] = {0.0, 0.0, 0.0}, upper[] = {1.0, 1.0, 1.0};
  int cells[] = {2, 2, 2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // projection updater for dist-function
  struct gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis, 5, 1, f_3d_p3, NULL);
  struct gkyl_proj_on_basis *projDistg = gkyl_proj_on_basis_new(&grid, &basis, 5, 1, g_3d_p3, NULL);

  // create array range: no ghost-cells in velocity space
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  // create distribution function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *distg = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  
  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &arr_range, distf);
  gkyl_proj_on_basis_advance(projDistg, 0.0, &arr_range, distg);

  struct gkyl_array *h = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(h, 0.0);

  struct gkyl_array *f_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(f_bar, 0.0);

  struct gkyl_array *g_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(g_bar, 0.0);

  // allocate memory
  gkyl_dg_bin_op_mem *mem = gkyl_dg_bin_op_mem_new(f_bar->size, basis.num_basis);  

  // h = f*g
  gkyl_dg_mul_op(&basis, 0, h, 0, distf, 0, distg);
  // f_bar = h/g = f
  gkyl_dg_div_op(mem, &basis, 0, f_bar, 0, h, 0, distg);
  // g_bar = h/f = g
  gkyl_dg_div_op(mem, &basis, 0, g_bar, 0, h, 0, distf);

  for (size_t i=0; i<arr_range.volume; ++i) {
    const double *f_d = gkyl_array_cfetch(distf, i);
    const double *fbar_d = gkyl_array_cfetch(f_bar, i);
    const double *g_d = gkyl_array_cfetch(distg, i);
    const double *gbar_d = gkyl_array_cfetch(g_bar, i);
    for (int k=0; k<basis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare(f_d[k], fbar_d[k], 1e-10) );
      TEST_CHECK( gkyl_compare(g_d[k], gbar_d[k], 1e-10) );
    }
  }

  // Test range methods
  gkyl_array_clear(h, 0.0);
  gkyl_array_clear(f_bar, 0.0);
  gkyl_array_clear(g_bar, 0.0);
  // h = f*g
  gkyl_dg_mul_op_range(&basis, 0, h, 0, distf, 0, distg, &arr_range);
  // f_bar = h/g = f
  gkyl_dg_div_op_range(mem, &basis, 0, f_bar, 0, h, 0, distg, &arr_range);
  // g_bar = h/f = g
  gkyl_dg_div_op_range(mem, &basis, 0, g_bar, 0, h, 0, distf, &arr_range);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &arr_range);

  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&arr_range, iter.idx);
    const double *f_d = gkyl_array_cfetch(distf, loc);
    const double *fbar_d = gkyl_array_cfetch(f_bar, loc);
    const double *g_d = gkyl_array_cfetch(distg, loc);
    const double *gbar_d = gkyl_array_cfetch(g_bar, loc);
    for (int k=0; k<basis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare(f_d[k], fbar_d[k], 1e-10) );
      TEST_CHECK( gkyl_compare(g_d[k], gbar_d[k], 1e-10) );
    }
  }
  
  gkyl_proj_on_basis_release(projDistf);
  gkyl_proj_on_basis_release(projDistg);
  gkyl_array_release(distf);
  gkyl_array_release(distg);
  gkyl_array_release(f_bar);
  gkyl_array_release(g_bar);
  gkyl_array_release(h);
  gkyl_dg_bin_op_mem_release(mem);  
}

// conf-phase accumulate subtraction test: pout = g(x,vx,vy) - f(x)*ones
// where f(x)=x, g(x,vx,vy)=sin(x)*cos(vx)+vy, ones=1
// verifies gkyl_dg_mul_conf_phase_op_accumulate_range with a=-1 and pop=constant 1.

void f_sub_eval(double t, const double *xn, double* restrict fout, void *ctx) {
  fout[0] = xn[0]; // f(x) = x
}

void g_eval(double t, const double *xn, double* restrict fout, void *ctx) {
  double x = xn[0], vx = xn[1], vy = xn[2];
  fout[0] = sin(x)*cos(vx) + vy; // g(x,vx,vy) = sin(x)*cos(vx) + vy
}

void ones_eval(double t, const double *xn, double* restrict fout, void *ctx) {
  fout[0] = 1.0;
}

void h_sol_eval(double t, const double *xn, double* restrict fout, void *ctx) {
  double x = xn[0], vx = xn[1], vy = xn[2];
  fout[0] = sin(x)*cos(vx) + vy - x; // h = g - f
}

void
test_subspace_accumulate(int poly_order, bool use_gpu)
{
  // Main dimension
  int dim = 3;
  double lower[] = {0.0, -M_PI/2.0, -M_PI/2.0};
  double upper[] = {M_PI, M_PI/2.0, M_PI/2.0};
  int cells[] = {8, 4, 4};
  // Sub dimension: 1D, x in [0, pi]
  int subdim = 1;
  double sublower[] = {0.0}, subupper[] = {M_PI};
  int ccells[] = {8};

  // grids
  struct gkyl_rect_grid fullgrid, subgrid;
  gkyl_rect_grid_init(&fullgrid, dim, lower, upper, cells);
  gkyl_rect_grid_init(&subgrid, subdim, sublower, subupper, ccells);

  // basis
  struct gkyl_basis basis, subbasis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);
  gkyl_cart_modal_serendip(&subbasis, subdim, poly_order);

  // ranges (no ghost cells)
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&fullgrid, nghost, &local_ext, &local);

  struct gkyl_range sub_local, sub_local_ext;
  gkyl_create_grid_ranges(&subgrid, nghost, &sub_local_ext, &sub_local);

  // projection updaters
  struct gkyl_proj_on_basis *proj_f = gkyl_proj_on_basis_new(&subgrid, &subbasis, poly_order+1, 1, f_sub_eval, NULL);
  struct gkyl_proj_on_basis *proj_g = gkyl_proj_on_basis_new(&fullgrid, &basis, poly_order+1, 1, g_eval, NULL);
  struct gkyl_proj_on_basis *proj_ones = gkyl_proj_on_basis_new(&fullgrid, &basis, poly_order+1, 1, ones_eval, NULL);
  struct gkyl_proj_on_basis *proj_h = gkyl_proj_on_basis_new(&fullgrid, &basis, poly_order+1, 1, h_sol_eval, NULL);

  // allocate arrays
  struct gkyl_array *f_sub_ho = gkyl_array_new(GKYL_DOUBLE, subbasis.num_basis, sub_local.volume);
  struct gkyl_array *g_full_ho = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local.volume);
  struct gkyl_array *ones_ho = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local.volume);
  struct gkyl_array *h_sol_ho = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local.volume);
  struct gkyl_array *pout_ho = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local.volume);

  // project onto basis with host arrays
  gkyl_proj_on_basis_advance(proj_f, 0.0, &sub_local, f_sub_ho);
  gkyl_proj_on_basis_advance(proj_g, 0.0, &local, g_full_ho);
  gkyl_proj_on_basis_advance(proj_ones, 0.0, &local, ones_ho);
  gkyl_proj_on_basis_advance(proj_h, 0.0, &local, h_sol_ho);

  // If GPU is enabled, make device copies of arrays and copy data to device
  struct gkyl_array *f_sub = use_gpu ? 
    gkyl_array_cu_dev_new(GKYL_DOUBLE, subbasis.num_basis, sub_local.volume) : gkyl_array_acquire(f_sub_ho);
  struct gkyl_array *g_full = use_gpu ? 
    gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local.volume) : gkyl_array_acquire(g_full_ho);
  struct gkyl_array *ones_phase = use_gpu ? 
    gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local.volume) : gkyl_array_acquire(ones_ho);
  struct gkyl_array *h_expect = use_gpu ? 
    gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local.volume) : gkyl_array_acquire(h_sol_ho);
  struct gkyl_array *pout = use_gpu ? 
    gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local.volume) : gkyl_array_acquire(pout_ho);

  // copy data to device if using GPU
  gkyl_array_copy(f_sub, f_sub_ho);
  gkyl_array_copy(g_full, g_full_ho);
  gkyl_array_copy(ones_phase, ones_ho);

  // pout = g, then pout += -1 * f_sub * ones_phase  =>  pout = g - f
  gkyl_array_copy(pout, g_full);
  gkyl_dg_mul_conf_phase_op_accumulate_range(&subbasis, &basis, pout, -1.0,
    f_sub, ones_phase, &sub_local, &local);

  // Copy pout back to host if using GPU
  gkyl_array_copy(pout_ho, pout);

  // check: pout should equal h_expect coefficient-by-coefficient
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&local, iter.idx);
    const double *pout_d = gkyl_array_cfetch(pout_ho,     loc);
    const double *h_d    = gkyl_array_cfetch(h_sol_ho, loc);
    for (int k = 0; k < basis.num_basis; ++k)
      TEST_CHECK( gkyl_compare(h_d[k], pout_d[k], 1e-11) );
  }

  // cleanup
  gkyl_proj_on_basis_release(proj_f);
  gkyl_proj_on_basis_release(proj_g);
  gkyl_proj_on_basis_release(proj_ones);
  gkyl_proj_on_basis_release(proj_h);
  gkyl_array_release(f_sub_ho);
  gkyl_array_release(g_full_ho);
  gkyl_array_release(ones_ho);
  gkyl_array_release(h_sol_ho);
  gkyl_array_release(pout_ho);
}

void test_bin_ops_conf_phase_accumulate_subtract_p1_ho() { test_subspace_accumulate(1, false); }
void test_bin_ops_conf_phase_accumulate_subtract_p2_ho() { test_subspace_accumulate(2, false); }

// Cuda specific tests
#ifdef GKYL_HAVE_CUDA

void test_bin_ops_1d_p1_dev(){ test_1d(1, true); }
void test_bin_ops_1d_p2_dev(){ test_1d(2, true); }
void test_bin_ops_1d_p3_dev(){ test_1d(3, true); }

void test_bin_ops_inv_1d_p1_dev(){ test_inv_1d(1, true); }

void test_bin_ops_2d_p1_dev(){ test_2d(1, true); }
void test_bin_ops_2d_p2_dev(){ test_2d(2, true); }
void test_bin_ops_2d_p3_dev(){ test_2d(3, true); }

void test_bin_ops_inv_2d_p1_dev(){ test_inv_2d(1, true); }

void test_bin_ops_3d_p1_dev(){ test_3d(1, true); }
void test_bin_ops_3d_p2_dev(){ test_3d(2, true); }

void
test_bin_ops_3d_p3_dev()
{
  int poly_order = 3;
  double lower[] = {0.0, 0.0, 0.0}, upper[] = {1.0, 1.0, 1.0};
  int cells[] = {2, 2, 2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // projection updater for dist-function
  struct gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis, 5, 1, f_3d_p3, NULL);
  struct gkyl_proj_on_basis *projDistg = gkyl_proj_on_basis_new(&grid, &basis, 5, 1, g_3d_p3, NULL);

  // create array range: no ghost-cells in velocity space
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  // create distribution function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *distg = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);

  // make device copies of arrays
  struct gkyl_array *distf_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *distg_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);  

  // project distribution function on basis on CPU
  gkyl_proj_on_basis_advance(projDistf, 0.0, &arr_range, distf);
  gkyl_proj_on_basis_advance(projDistg, 0.0, &arr_range, distg);

  // copy host array to device
  gkyl_array_copy(distf_cu, distf);
  gkyl_array_copy(distg_cu, distg);

  // Product array only needs to be initialized on GPU
  struct gkyl_array *h_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(h_cu, 0.0);

  struct gkyl_array *f_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *f_bar_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(f_bar, 0.0);
  gkyl_array_clear(f_bar_cu, 0.0);

  struct gkyl_array *g_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *g_bar_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(g_bar, 0.0);
  gkyl_array_clear(g_bar_cu, 0.0);

  gkyl_dg_bin_op_mem *mem = gkyl_dg_bin_op_mem_cu_dev_new(f_bar->size, basis.num_basis);

  // h = f*g
  gkyl_dg_mul_op(&basis, 0, h_cu, 0, distf_cu, 0, distg_cu);
  // f_bar = h/g = f
  gkyl_dg_div_op(mem, &basis, 0, f_bar_cu, 0, h_cu, 0, distg_cu);
  // g_bar = h/f = g
  gkyl_dg_div_op(mem, &basis, 0, g_bar_cu, 0, h_cu, 0, distf_cu);

  // copy from device and check if things are ok
  gkyl_array_copy(f_bar, f_bar_cu);
  gkyl_array_copy(g_bar, g_bar_cu);
  for (size_t i=0; i<arr_range.volume; ++i) {
    const double *f_d = gkyl_array_cfetch(distf, i);
    const double *fbar_d = gkyl_array_cfetch(f_bar, i);
    const double *g_d = gkyl_array_cfetch(distg, i);
    const double *gbar_d = gkyl_array_cfetch(g_bar, i);
    for (int k=0; k<basis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare(f_d[k], fbar_d[k], 1e-10) );
      TEST_CHECK( gkyl_compare(g_d[k], gbar_d[k], 1e-10) );
    }
  }

  // Test range methods
  gkyl_array_clear(h_cu, 0.0);
  gkyl_array_clear(f_bar_cu, 0.0);
  gkyl_array_clear(g_bar_cu, 0.0);
  // h = f*g
  gkyl_dg_mul_op_range(&basis, 0, h_cu, 0, distf_cu, 0, distg_cu, &arr_range);
  // f_bar = h/g = f
  gkyl_dg_div_op_range(mem, &basis, 0, f_bar_cu, 0, h_cu, 0, distg_cu, &arr_range);
  // g_bar = h/f = g
  gkyl_dg_div_op_range(mem, &basis, 0, g_bar_cu, 0, h_cu, 0, distf_cu, &arr_range);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &arr_range);
  // copy from device and check if things are ok
  gkyl_array_copy(f_bar, f_bar_cu);
  gkyl_array_copy(g_bar, g_bar_cu);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&arr_range, iter.idx);
    const double *f_d = gkyl_array_cfetch(distf, loc);
    const double *fbar_d = gkyl_array_cfetch(f_bar, loc);
    const double *g_d = gkyl_array_cfetch(distg, loc);
    const double *gbar_d = gkyl_array_cfetch(g_bar, loc);
    for (int k=0; k<basis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare(f_d[k], fbar_d[k], 1e-10) );
      TEST_CHECK( gkyl_compare(g_d[k], gbar_d[k], 1e-10) );
    }
  }
  
  gkyl_proj_on_basis_release(projDistf);
  gkyl_proj_on_basis_release(projDistg);
  gkyl_array_release(distf);
  gkyl_array_release(distg);
  gkyl_array_release(f_bar);
  gkyl_array_release(g_bar);
  gkyl_array_release(distf_cu);
  gkyl_array_release(distg_cu);
  gkyl_array_release(f_bar_cu);
  gkyl_array_release(g_bar_cu);
  gkyl_array_release(h_cu);
  gkyl_dg_bin_op_mem_release(mem);
}

void test_bin_ops_conf_phase_accumulate_subtract_p1_dev() { test_subspace_accumulate(1, true); }
void test_bin_ops_conf_phase_accumulate_subtract_p2_dev() { test_subspace_accumulate(2, true); }

#endif

TEST_LIST = {
  { "test_bin_ops_1d_p1_ho", test_bin_ops_1d_p1_ho },
  { "test_bin_ops_inv_1d_p1_ho", test_bin_ops_inv_1d_p1_ho },
  { "test_bin_ops_1d_p2_ho", test_bin_ops_1d_p2_ho },
  { "test_bin_ops_1d_p3_ho", test_bin_ops_1d_p3_ho },
  { "test_bin_ops_2d_p1_ho", test_bin_ops_2d_p1_ho },
  { "test_bin_ops_inv_2d_p1_ho", test_bin_ops_inv_2d_p1_ho },
  { "test_bin_ops_2d_p2_ho", test_bin_ops_2d_p2_ho },
  { "test_bin_ops_2d_p3_ho", test_bin_ops_2d_p3_ho },
  { "test_bin_ops_3d_p1_ho", test_bin_ops_3d_p1_ho },
  { "test_bin_ops_3d_p2_ho", test_bin_ops_3d_p2_ho },
  { "test_bin_ops_3d_p3_ho", test_bin_ops_3d_p3_ho },
  { "test_bin_ops_4d_p1_ho", test_bin_ops_4d_p1_ho },
  { "test_bin_ops_4d_p2_ho", test_bin_ops_4d_p2_ho },
  { "test_bin_ops_conf_phase_accumulate_subtract_p1_ho", test_bin_ops_conf_phase_accumulate_subtract_p1_ho },
  { "test_bin_ops_conf_phase_accumulate_subtract_p2_ho", test_bin_ops_conf_phase_accumulate_subtract_p2_ho },
#ifdef GKYL_HAVE_CUDA
  { "test_bin_ops_1d_p1_dev", test_bin_ops_1d_p1_dev },
  { "test_bin_ops_inv_1d_p1_dev", test_bin_ops_inv_1d_p1_dev },
  { "test_bin_ops_1d_p2_dev", test_bin_ops_1d_p2_dev },
  { "test_bin_ops_1d_p3_dev", test_bin_ops_1d_p3_dev },
  { "test_bin_ops_2d_p1_dev", test_bin_ops_2d_p1_dev },
  { "test_bin_ops_inv_2d_p1_dev", test_bin_ops_inv_2d_p1_dev },
  { "test_bin_ops_2d_p2_dev", test_bin_ops_2d_p2_dev },
  { "test_bin_ops_2d_p3_dev", test_bin_ops_2d_p3_dev },
  { "test_bin_ops_3d_p1_dev", test_bin_ops_3d_p1_dev },
  { "test_bin_ops_3d_p2_dev", test_bin_ops_3d_p2_dev },
  { "test_bin_ops_3d_p3_dev", test_bin_ops_3d_p3_dev },
  { "test_bin_ops_conf_phase_accumulate_subtract_p1_dev", test_bin_ops_conf_phase_accumulate_subtract_p1_dev },
  { "test_bin_ops_conf_phase_accumulate_subtract_p2_dev", test_bin_ops_conf_phase_accumulate_subtract_p2_dev },
#endif
  { NULL, NULL },
};