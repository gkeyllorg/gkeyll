#include <acutest.h>
#include <math.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_ops.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <gkyl_basis.h>
#include <gkyl_deflate_zsurf.h>
#include <gkyl_nodal_ops.h>

#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_fem_poisson.h>

#include <gkyl_fem_parproj.h>

#include <gkyl_deflated_dg_bin_ops.h>


void
proj_jac(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0];
  double z = xn[1];
  fout[0] = 1.0 + z*z*cos(x);
}

void
proj_rho(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0];
  double z = xn[1];
  fout[0] = z*sin((2.*M_PI/(2.*M_PI))*x);
}

// Check continuity along last dim in 2x
void check_continuity_2x(struct gkyl_rect_grid grid, struct gkyl_range range, struct gkyl_basis basis, struct gkyl_array *field)
{
  struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, grid.ndim, basis.num_basis);
  basis.node_list(gkyl_array_fetch(nodes, 0));
  // Check continuity
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  int nidx;
  long lin_nidx;
  int idx[3];
  const double *node_i  ;

  while (gkyl_range_iter_next(&iter)) {
    if (iter.idx[1] != range.upper[1]) {
      long lidx = gkyl_range_idx(&range, iter.idx);
      idx[0] = iter.idx[0];
      idx[1] = iter.idx[1] + 1;
      long lidx_up = gkyl_range_idx(&range, idx);
      double *arr = gkyl_array_fetch(field, lidx);
      double *arr_up = gkyl_array_fetch(field, lidx_up);
      node_i  = gkyl_array_cfetch(nodes, 2);
      double temp1 = basis.eval_expand(node_i, arr);
      node_i  = gkyl_array_cfetch(nodes, 3);
      double temp2 = basis.eval_expand(node_i, arr);
      node_i  = gkyl_array_cfetch(nodes, 0);
      double temp_up1 = basis.eval_expand(node_i, arr_up);
      node_i  = gkyl_array_cfetch(nodes, 1);
      double temp_up2 = basis.eval_expand(node_i, arr_up);
      TEST_CHECK( gkyl_compare(temp1, temp_up1, 1e-12) );
      TEST_CHECK( gkyl_compare(temp2, temp_up2, 1e-12) );
    }
  }
  gkyl_array_release(nodes);
}

void check_same(struct gkyl_range range, struct gkyl_basis basis, struct gkyl_array *field1, struct gkyl_array* field2)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(&range, iter.idx);
    const double *f1 = gkyl_array_cfetch(field1, lidx);
    const double *f2 = gkyl_array_cfetch(field2, lidx);
    for(int i = 0; i< basis.num_basis; i++)
      TEST_CHECK( gkyl_compare(f1[i], f2[i], 1e-10) );
  }
}



void
test_bop(bool use_gpu)
{
  int c_lop = 0;
  int c_rop = 0;
  int c_oop = 0;
  // create the 2d field
  // create xz grid
  double lower[] = { -M_PI, 0.0 }, upper[] = { M_PI, 1.0 };
  int cells[] = { 12, 8 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  //ranges
  struct gkyl_range local, local_ext;
  int nghost[GKYL_MAX_CDIM] = { 1, 1 };
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);
  struct gkyl_range global_sub_range = local;

  // basis function
  int poly_order = 1;
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 2, poly_order);

  struct gkyl_basis *basis_on_dev = use_gpu? gkyl_cart_modal_serendip_cu_dev_new(2, poly_order)
                                             : gkyl_cart_modal_serendip_new(2, poly_order);

  // project initial functions on 2d field
  struct gkyl_array *rho = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&grid, &basis, 2, 1, &proj_rho, 0);
  gkyl_proj_on_basis_advance(proj, 0.0, &local, rho);
  gkyl_proj_on_basis_release(proj);

  struct gkyl_array *jac= gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_eval_on_nodes *eon = gkyl_eval_on_nodes_new(&grid, &basis, 1, &proj_jac, 0);
  gkyl_eval_on_nodes_advance(eon, 0.0, &local, jac);
  gkyl_eval_on_nodes_release(eon);

  //gkyl_grid_sub_array_write(&grid, &local, 0, rho, "rho.gkyl");
  //gkyl_grid_sub_array_write(&grid, &local, 0, jac, "jac.gkyl");

  // Get C = rho.J
  struct gkyl_array *Cxz = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_dg_mul_op(&basis, 0, Cxz, 0,rho, 0, jac);
  //gkyl_grid_sub_array_write(&grid, &local, 0, Cxz, "Cxz.gkyl");

  struct gkyl_array *Cxz_dev =  use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume) : gkyl_array_acquire(Cxz);
  gkyl_array_copy(Cxz_dev, Cxz);
  struct gkyl_array *jac_dev = use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume) : gkyl_array_acquire(jac);
  gkyl_array_copy(jac_dev, jac);

  // Allocate Exz which will store C/J, the smoothed version of it, and
  // Fxz = Exz*J which will hold the final result
  struct gkyl_array *Exz = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  struct gkyl_array *Exz_smooth = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  struct gkyl_array *Fxz = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);

  struct gkyl_array *Exz_dev = use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume) : gkyl_array_acquire(Exz);
  struct gkyl_array *Exz_smooth_dev = use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume) : gkyl_array_acquire(Exz_smooth);
  struct gkyl_array *Fxz_dev = use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume) : gkyl_array_acquire(Fxz);

  // Make deflated array operator
  struct gkyl_deflated_dg_bin_ops* operator = gkyl_deflated_dg_bin_ops_new(grid, basis_on_dev, basis, local, use_gpu);

  // Divide
  gkyl_deflated_dg_bin_ops_div(operator, c_oop, Exz_dev, c_lop, Cxz_dev, c_rop, jac_dev);

  gkyl_array_copy(Exz, Exz_dev);
  //gkyl_grid_sub_array_write(&grid, &local, 0, Exz, "Exz.gkyl");

  // Smooth E
  struct gkyl_fem_parproj *fem_parproj = gkyl_fem_parproj_new(&local, &grid, &basis, GKYL_FEM_PARPROJ_NONE, 0, 0, 0, use_gpu);
  gkyl_fem_parproj_set_rhs(fem_parproj, Exz_dev, Exz_dev);
  gkyl_fem_parproj_solve(fem_parproj, Exz_smooth_dev);
  gkyl_fem_parproj_release(fem_parproj);

  gkyl_array_copy(Exz_smooth, Exz_smooth_dev);
  //gkyl_grid_sub_array_write(&grid, &local, 0, Exz_smooth, "Exz_smooth.gkyl");

  // Multiply
  gkyl_deflated_dg_bin_ops_mul(operator, c_oop, Fxz_dev, c_lop, Exz_smooth_dev, c_rop, jac_dev);
  
  gkyl_array_copy(Fxz, Fxz_dev);
  //gkyl_grid_sub_array_write(&grid, &local, 0, Fxz, "Fxz.gkyl");

  check_continuity_2x(grid, local, basis, Fxz);
  
  if (use_gpu)
    gkyl_cart_modal_basis_release_cu(basis_on_dev);
  else
    gkyl_cart_modal_basis_release(basis_on_dev);

  gkyl_array_release(rho);
  gkyl_array_release(jac);
  gkyl_array_release(Cxz);
  gkyl_array_release(Exz);
  gkyl_array_release(Exz_smooth);
  gkyl_array_release(Fxz);
  gkyl_array_release(Cxz_dev);
  gkyl_array_release(jac_dev);
  gkyl_array_release(Exz_dev);
  gkyl_array_release(Exz_smooth_dev);
  gkyl_array_release(Fxz_dev);
  


  gkyl_deflated_dg_bin_ops_release(operator);

}

void test_deflated_bop_ho(void) { test_bop(false); }
void test_deflated_bop_dev(void) { test_bop(true); }

TEST_LIST = {
  { "test_deflated_bop_ho", test_deflated_bop_ho},
#ifdef GKYL_HAVE_CUDA
  { "test_deflated_bop_dev", test_deflated_bop_dev},
#endif
  { NULL, NULL },
};

