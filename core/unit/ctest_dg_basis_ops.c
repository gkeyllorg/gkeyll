#include <acutest.h>

#include <gkyl_basis.h>
#include <gkyl_dg_basis_ops.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_range.h>
#include <gkyl_proj_on_basis.h>


static inline double sq(double x) { return x*x; }
static inline double cub(double x) { return x*x*x; }

// Allocate array (filled with zeros).
static struct gkyl_array*
mkarr(bool use_gpu, long nc, long size)
{
  struct gkyl_array* a = use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
	                        : gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

static void
eval_array_at_coord_1d_func(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xc[0];
  double Lx = 5.0;
  fout[0] = 3.3*0.5*(1.0+cos((2*M_PI/Lx)*x));
}

void
test_eval_array_at_coord_1d_p_hodev(int poly_order, bool use_gpu)
{
  double lower[] = { 0.0 }, upper[] = { 5.0 };
  int cells[] = { 8 };
  double eval_coord_ho[] = {3.0}; // Result = 3.3*0.5*(1+sin((2*pi/5)*3.0)) = 0.315121959

  int ndim = sizeof(lower)/sizeof(lower[0]);

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  int nghost[GKYL_MAX_DIM];
  for (int d=0; d<ndim; d++)
    nghost[d] = 1;

  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);

  struct gkyl_basis basis;
  gkyl_cart_modal_tensor(&basis, ndim, poly_order);
  struct gkyl_basis *basis_on_dev = use_gpu? gkyl_cart_modal_tensor_cu_dev_new(ndim, poly_order)
                                           : gkyl_cart_modal_tensor_new(ndim, poly_order);

  struct gkyl_array *fld = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *fld_ho = use_gpu? mkarr(false, basis.num_basis, local_ext.volume)
	                             : gkyl_array_acquire(fld);

  gkyl_proj_on_basis *proj_op = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_array_at_coord_1d_func, NULL);
  gkyl_proj_on_basis_advance(proj_op, 0.0, &local, fld_ho);
  gkyl_array_copy(fld, fld_ho);

  double *eval_coord, *fld_at_coord;
  if (use_gpu) {
    eval_coord = gkyl_cu_malloc(sizeof(double));
    fld_at_coord = gkyl_cu_malloc(sizeof(double));
    gkyl_cu_memcpy(eval_coord, eval_coord_ho, ndim*sizeof(double), GKYL_CU_MEMCPY_H2D);
  }
  else {
    eval_coord = gkyl_malloc(sizeof(double));
    fld_at_coord = gkyl_malloc(sizeof(double));
    memcpy(eval_coord, eval_coord_ho, ndim*sizeof(double));
  }

  gkyl_dg_basis_ops_eval_array_at_coord_comp(fld, eval_coord, basis_on_dev, &grid, &local, fld_at_coord);

  double fld_at_coord_ho[1];
  if (use_gpu)
    gkyl_cu_memcpy(fld_at_coord_ho, fld_at_coord, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(fld_at_coord_ho, fld_at_coord, sizeof(double));

  double accepted_val = 3.121168140e-01;
  TEST_CHECK( gkyl_compare_double(fld_at_coord_ho[0], accepted_val, 1.0e-8) );
  TEST_MSG( "Got: %.9e | Expected: %.9e\n",fld_at_coord_ho[0], accepted_val );

  if (use_gpu) {
    gkyl_cart_modal_basis_release_cu(basis_on_dev);
    gkyl_cu_free(eval_coord);
    gkyl_cu_free(fld_at_coord);
  }
  else {
    gkyl_cart_modal_basis_release(basis_on_dev);
    gkyl_free(eval_coord);
    gkyl_free(fld_at_coord);
  }
  gkyl_proj_on_basis_release(proj_op);
  gkyl_array_release(fld_ho);
  gkyl_array_release(fld);
}

void
test_cubic_1d(void)
{
  double val[2] = { 1.0, 2.0 };
  double grad[2] = { -1.0, -2.0 };
  double coeff[4] = { 0.0 };
  
  gkyl_dg_calc_cubic_1d(val, grad, coeff);

  struct gkyl_basis b3;
  gkyl_cart_modal_serendip(&b3, 1, 3);

  TEST_CHECK( gkyl_compare_double(val[0], b3.eval_expand((double[1]) { -1.0 }, coeff), 1.0e-15) );
  TEST_CHECK( gkyl_compare_double(val[1], b3.eval_expand((double[1]) { 1.0 }, coeff), 1.0e-15) );

  TEST_CHECK( gkyl_compare_double(grad[0], b3.eval_grad_expand(0, (double[1]) { -1.0 }, coeff), 1.0e-15) );
  TEST_CHECK( gkyl_compare_double(grad[1], b3.eval_grad_expand(0, (double[1]) { 1.0 }, coeff), 1.0e-15) );
}

void
test_cubic_2d(void)
{
  double val[4] = { 1.0, 2.0, 3.0, 4.0 };
  double gradx[4] = { -1.0, -2.0, -3.0, -4.0 };
  double grady[4] = { 1.0, 2.0, 3.0, 4.0 };
  double gradxy[4] = { 1.5, 2.5, 3.5, 4.5 };
  double coeff[16] = { 0.0 };
  
  gkyl_dg_calc_cubic_2d(val, gradx, grady, gradxy, coeff);

  struct gkyl_basis b3;
  gkyl_cart_modal_tensor(&b3, 2, 3);

  TEST_CHECK( gkyl_compare_double(val[0], b3.eval_expand((double[2]) { -1.0, -1.0 }, coeff), 1.0e-15) );
  TEST_CHECK( gkyl_compare_double(val[1], b3.eval_expand((double[2]) { -1.0, 1.0 }, coeff), 1.0e-15) );
  TEST_CHECK( gkyl_compare_double(val[2], b3.eval_expand((double[2]) { 1.0, -1.0 }, coeff), 1.0e-15) );
  TEST_CHECK( gkyl_compare_double(val[3], b3.eval_expand((double[2]) { 1.0, 1.0 }, coeff), 1.0e-15) );

  TEST_CHECK( gkyl_compare_double(gradx[0], b3.eval_grad_expand(0, (double[2]) { -1.0, -1.0 }, coeff), 1.0e-15) );
  TEST_CHECK( gkyl_compare_double(gradx[1], b3.eval_grad_expand(0, (double[2]) { -1.0, 1.0 }, coeff), 1.0e-15) );
  TEST_CHECK( gkyl_compare_double(gradx[2], b3.eval_grad_expand(0, (double[2]) { 1.0, -1.0 }, coeff), 1.0e-15) );
  TEST_CHECK( gkyl_compare_double(gradx[3], b3.eval_grad_expand(0, (double[2]) { 1.0, 1.0 }, coeff), 1.0e-15) );

  TEST_CHECK( gkyl_compare_double(grady[0], b3.eval_grad_expand(1, (double[2]) { -1.0, -1.0 }, coeff), 1.0e-15) );
  TEST_CHECK( gkyl_compare_double(grady[1], b3.eval_grad_expand(1, (double[2]) { -1.0, 1.0 }, coeff), 1.0e-15) );
  TEST_CHECK( gkyl_compare_double(grady[2], b3.eval_grad_expand(1, (double[2]) { 1.0, -1.0 }, coeff), 1.0e-15) );
  TEST_CHECK( gkyl_compare_double(grady[3], b3.eval_grad_expand(1, (double[2]) { 1.0, 1.0 }, coeff), 1.0e-15) );
}

void
test_cubic_evalf_2d(void)
{
  double lower[] = { 0.0, 0.0 }, upper[] = { 5.0, 5.0 };
  int cells[] = { 8, 8 };

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  // nodal grid used in IO so we can plot things
  double nc_lower[] = { lower[0] - 0.5*grid.dx[0], lower[1] - 0.5*grid.dx[1] };
  double nc_upper[] = { upper[0] + 0.5*grid.dx[0], upper[1] + 0.5*grid.dx[1] };
  int nc_cells[] = { cells[0] + 1, cells[1] + 1 };
  struct gkyl_rect_grid nc_grid;
  gkyl_rect_grid_init(&nc_grid, 2, nc_lower, nc_upper, nc_cells);

  struct gkyl_range local, local_ext;
  int nghost[GKYL_MAX_CDIM] = { 0, 0 };  
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);

  struct gkyl_range nc_local, nc_local_ext;
  gkyl_create_grid_ranges(&nc_grid, nghost, &nc_local_ext, &nc_local);

  struct gkyl_basis basis;
  gkyl_cart_modal_tensor(&basis, 2, 3);

  struct gkyl_array *psi_nodal = gkyl_array_new(GKYL_DOUBLE, 1, (cells[0]+1)*(cells[1]+1));
  struct gkyl_array *psi_cubic = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_dg_basis_op_mem *mem = gkyl_dg_alloc_cubic_2d(cells);

  double xn[2];
  
  // initialize 2D nodal values f = x^y^2 which we can represent exactly
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &nc_local);
  while (gkyl_range_iter_next(&iter)) {
    long nidx = gkyl_range_idx(&nc_local, iter.idx);
    
    gkyl_rect_grid_ll_node(&grid, iter.idx, xn);
    
    double *pn = gkyl_array_fetch(psi_nodal, nidx);
    pn[0] = sq(xn[0])*sq(xn[1]);
  }
  // compute cubic expansion
  gkyl_dg_calc_cubic_2d_from_nodal_vals(mem, cells, grid.dx,
    psi_nodal, psi_cubic);
    
  // compute evalf function from nodal values
  struct gkyl_basis_ops_evalf *evf = gkyl_dg_basis_ops_evalf_new(&grid, psi_nodal);

  // Test wgrad function: f, df/dx, and df/dy
  xn[0] = 2.5;
  xn[1] = 2.5;
  double fout[4];
  evf->eval_cubic_wgrad(0.0, xn, fout, evf->ctx);
  TEST_CHECK( gkyl_compare_double(fout[0], sq(xn[0])*sq(xn[1]), 1.0e-14) );
  TEST_CHECK( gkyl_compare_double(fout[1], 2*xn[0]*sq(xn[1]), 1.0e-14) );
  TEST_CHECK( gkyl_compare_double(fout[2], 2*sq(xn[0])*xn[1], 1.0e-14) );

  // Test wgrad2 function: f, d2f/dx2, d2f/dy2, and d2f/dxdy
  evf->eval_cubic_wgrad2(0.0, xn, fout, evf->ctx);
  TEST_CHECK( gkyl_compare_double(fout[0], sq(xn[0])*sq(xn[1]), 1.0e-14) );
  TEST_CHECK( gkyl_compare_double(fout[1], 2*sq(xn[1]), 1.0e-14) );
  TEST_CHECK( gkyl_compare_double(fout[2], 2*sq(xn[0]), 1.0e-14) );
  TEST_CHECK( gkyl_compare_double(fout[3], 4*xn[0]*xn[1], 1.0e-14) );

  gkyl_array_release(psi_nodal);
  gkyl_array_release(psi_cubic);
  
  gkyl_dg_basis_op_mem_release(mem);
  gkyl_dg_basis_ops_evalf_release(evf);
}

void
test_eval_array_at_coord_1d_ho() {
  // p = 1
  test_eval_array_at_coord_1d_p_hodev(1, false);
}

#ifdef GKYL_HAVE_CUDA
void
test_eval_array_at_coord_1d_dev() {
  // p = 1
  test_eval_array_at_coord_1d_p_hodev(1, true);
}
#endif

TEST_LIST = {
  { "test_eval_array_at_coord_1d_ho", test_eval_array_at_coord_1d_ho },
  { "cubic_1d", test_cubic_1d },
  { "cubic_2d", test_cubic_2d },
  { "cubic_evalf_2d", test_cubic_evalf_2d },
#ifdef GKYL_HAVE_CUDA
  { "test_eval_array_at_coord_1d_dev", test_eval_array_at_coord_1d_dev },
#endif
  { NULL, NULL },
};

