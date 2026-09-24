#include <acutest.h>
#include <math.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_ops.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <gkyl_basis.h>
#include <gkyl_nodal_ops.h>
#include <gkyl_deflate_zsurf.h>


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
proj_func(double t, const double *xn, double *fout, void *ctx)
{
  fout[0] = cos(xn[0])*sin(xn[1]);
}

void
proj_func3d(double t, const double *xn, double *fout, void *ctx)
{
  fout[0] = cos(2*xn[0])*sin(xn[1])*xn[2]*xn[2]*xn[2];
}

void
test_p1_2x(){
  // create  grid, ranges, basis
  double lower[] = { 0.0, -1.5 }, upper[] = { 1.5, 1.5 };
  int cells[] = { 8, 16 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);
  struct gkyl_range local, local_ext;
  int nghost[GKYL_MAX_CDIM] = { 1, 1 };
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);
  int poly_order = 1;
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 2, poly_order);
  bool use_gpu = false;
#ifdef GKYL_HAVE_CUDA
  use_gpu = true;
  struct gkyl_basis *basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
  gkyl_cart_modal_serendip_cu_dev(basis_on_dev, 2, poly_order);
#else
  struct gkyl_basis *basis_on_dev = &basis;
#endif
  
  // Project initial function
  struct gkyl_array *funcdg = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  struct gkyl_eval_on_nodes *eon = gkyl_eval_on_nodes_new(&grid, &basis, 1, &proj_func, 0);
  gkyl_eval_on_nodes_advance(eon, 0.0, &local, funcdg);
  gkyl_grid_sub_array_write(&grid, &local, 0, funcdg, "proj_func.gkyl");
#ifdef GKYL_HAVE_CUDA
  struct gkyl_array *funcdg_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_array_copy(funcdg_dev, funcdg);
#else
  struct gkyl_array *funcdg_dev = funcdg;
#endif

  // Construct nrange and nodal field
  int nodes[3] = { 1, 1, 1 };
  for (int d=0; d<grid.ndim; ++d)
    nodes[d] = grid.cells[d] + 1;
  struct gkyl_range nrange;
  gkyl_range_init_from_shape(&nrange, grid.ndim, nodes);
  struct gkyl_array* nodal_fld = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
#ifdef GKYL_HAVE_CUDA
  struct gkyl_array *nodal_fld_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
#else
  struct gkyl_array *nodal_fld_dev = nodal_fld;
#endif


  // Create output field
  struct gkyl_array *funcdg2 = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
#ifdef GKYL_HAVE_CUDA
  struct gkyl_array *funcdg2_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
#else
  struct gkyl_array *funcdg2_dev = funcdg2;
#endif

  // Trnasform forward and back
  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&basis, &grid, use_gpu);
  gkyl_nodal_ops_m2n(n2m, basis_on_dev, &grid, &nrange, &local, 1, nodal_fld_dev, funcdg_dev, false);
  gkyl_nodal_ops_n2m(n2m, basis_on_dev, &grid, &nrange, &local, 1, nodal_fld_dev, funcdg2_dev, false);

#ifdef GKYL_HAVE_CUDA
  gkyl_array_copy(funcdg2, funcdg2_dev);
#endif
  check_same(local, basis, funcdg, funcdg2);
  gkyl_grid_sub_array_write(&grid, &local, 0, funcdg2, "proj_func2.gkyl");



#ifdef GKYL_HAVE_CUDA
  gkyl_cu_free(basis_on_dev);
  gkyl_array_release(funcdg_dev);
  gkyl_array_release(funcdg2_dev);
  gkyl_array_release(nodal_fld_dev);
#endif
  gkyl_array_release(funcdg);
  gkyl_array_release(funcdg2);
  gkyl_array_release(nodal_fld);
  gkyl_eval_on_nodes_release(eon);
  gkyl_nodal_ops_release(n2m);

}

void
test_p1_3x(){
  // create  grid, ranges, basis
  double lower[] = { 0.0, -1.5, -1.0}, upper[] = { 1.5, 1.5, 1.0 };
  int cells[] = { 2, 4, 2 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 3, lower, upper, cells);
  struct gkyl_range local, local_ext;
  int nghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);
  int poly_order = 1;
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 3, poly_order);
  bool use_gpu = false;
#ifdef GKYL_HAVE_CUDA
  use_gpu = true;
  struct gkyl_basis *basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
  gkyl_cart_modal_serendip_cu_dev(basis_on_dev, 3, poly_order);
#else
  struct gkyl_basis *basis_on_dev = &basis;
#endif
  
  // Project initial function
  struct gkyl_array *funcdg = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  struct gkyl_eval_on_nodes *eon = gkyl_eval_on_nodes_new(&grid, &basis, 1, &proj_func3d, 0);
  gkyl_eval_on_nodes_advance(eon, 0.0, &local, funcdg);
  gkyl_grid_sub_array_write(&grid, &local, 0, funcdg, "proj_func3d.gkyl");
#ifdef GKYL_HAVE_CUDA
  struct gkyl_array *funcdg_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_array_copy(funcdg_dev, funcdg);
#else
  struct gkyl_array *funcdg_dev = funcdg;
#endif

  // Construct nrange and nodal field
  int nodes[3] = { 1, 1, 1 };
  for (int d=0; d<grid.ndim; ++d)
    nodes[d] = grid.cells[d] + 1;
  struct gkyl_range nrange;
  gkyl_range_init_from_shape(&nrange, grid.ndim, nodes);
  struct gkyl_array* nodal_fld = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
#ifdef GKYL_HAVE_CUDA
  struct gkyl_array *nodal_fld_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
#else
  struct gkyl_array *nodal_fld_dev = nodal_fld;
#endif


  // Create output field
  struct gkyl_array *funcdg2 = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
#ifdef GKYL_HAVE_CUDA
  struct gkyl_array *funcdg2_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
#else
  struct gkyl_array *funcdg2_dev = funcdg2;
#endif

  // Trnasform forward and back
  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&basis, &grid, use_gpu);
  gkyl_nodal_ops_m2n(n2m, basis_on_dev, &grid, &nrange, &local, 1, nodal_fld_dev, funcdg_dev, false);
  gkyl_nodal_ops_n2m(n2m, basis_on_dev, &grid, &nrange, &local, 1, nodal_fld_dev, funcdg2_dev, false);

#ifdef GKYL_HAVE_CUDA
  gkyl_array_copy(funcdg2, funcdg2_dev);
#endif
  check_same(local, basis, funcdg, funcdg2);
  gkyl_grid_sub_array_write(&grid, &local, 0, funcdg2, "proj_func3d_2.gkyl");



#ifdef GKYL_HAVE_CUDA
  gkyl_cu_free(basis_on_dev);
  gkyl_array_release(funcdg_dev);
  gkyl_array_release(funcdg2_dev);
  gkyl_array_release(nodal_fld_dev);
#endif
  gkyl_array_release(funcdg);
  gkyl_array_release(funcdg2);
  gkyl_array_release(nodal_fld);
  gkyl_eval_on_nodes_release(eon);
  gkyl_nodal_ops_release(n2m);

}



void
test_p1_interior_2x(){
  // create  grid, ranges, basis
  double lower[] = { 0.0, -1.5 }, upper[] = { 1.5, 1.5 };
  int cells[] = { 2, 4 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);
  struct gkyl_range local, local_ext;
  int nghost[GKYL_MAX_CDIM] = { 1, 1 };
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);
  int poly_order = 1;
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 2, poly_order);
  bool use_gpu = false;
  struct gkyl_basis *basis_on_dev = &basis;
  
  // Project initial function
  struct gkyl_array *funcdg = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  struct gkyl_eval_on_nodes *eon = gkyl_eval_on_nodes_new(&grid, &basis, 1, &proj_func, 0);
  gkyl_eval_on_nodes_advance(eon, 0.0, &local, funcdg);
  gkyl_grid_sub_array_write(&grid, &local, 0, funcdg, "proj_func.gkyl");
  struct gkyl_array *funcdg_dev = funcdg;

  // Construct nrange and nodal field
  int nodes[3] = { 1, 1, 1 };
  for (int d=0; d<grid.ndim; ++d)
    nodes[d] = 2*grid.cells[d];
  struct gkyl_range nrange;
  gkyl_range_init_from_shape(&nrange, grid.ndim, nodes);
  struct gkyl_array* nodal_fld = gkyl_array_new(GKYL_DOUBLE, 1, nrange.volume);
  struct gkyl_array *nodal_fld_dev = nodal_fld;


  // Create output field
  struct gkyl_array *funcdg2 = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  struct gkyl_array *funcdg2_dev = funcdg2;

  // Trnasform forward and back
  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&basis, &grid, use_gpu);
  gkyl_nodal_ops_m2n(n2m, basis_on_dev, &grid, &nrange, &local, 1, nodal_fld_dev, funcdg_dev, true);
  gkyl_nodal_ops_n2m(n2m, basis_on_dev, &grid, &nrange, &local, 1, nodal_fld_dev, funcdg2_dev, true);

  check_same(local, basis, funcdg, funcdg2);
  gkyl_grid_sub_array_write(&grid, &local, 0, funcdg2, "proj_func2.gkyl");



  gkyl_array_release(funcdg);
  gkyl_array_release(funcdg2);
  gkyl_array_release(nodal_fld);
  gkyl_eval_on_nodes_release(eon);
  gkyl_nodal_ops_release(n2m);

}

void
test_p1_interior_3x(){
  // create  grid, ranges, basis
  double lower[] = { 0.0, -1.5, -1.0}, upper[] = { 1.5, 1.5, 1.0 };
  int cells[] = { 2, 4, 2 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 3, lower, upper, cells);
  struct gkyl_range local, local_ext;
  int nghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);
  int poly_order = 1;
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 3, poly_order);
  bool use_gpu = false;
  struct gkyl_basis *basis_on_dev = &basis;
  
  // Project initial function
  struct gkyl_array *funcdg = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  struct gkyl_eval_on_nodes *eon = gkyl_eval_on_nodes_new(&grid, &basis, 1, &proj_func3d, 0);
  gkyl_eval_on_nodes_advance(eon, 0.0, &local, funcdg);
  gkyl_grid_sub_array_write(&grid, &local, 0, funcdg, "proj_func3d.gkyl");
  struct gkyl_array *funcdg_dev = funcdg;

  // Construct nrange and nodal field
  int nodes[3] = { 1, 1, 1 };
  for (int d=0; d<grid.ndim; ++d)
    nodes[d] = grid.cells[d]*2;
  struct gkyl_range nrange;
  gkyl_range_init_from_shape(&nrange, grid.ndim, nodes);
  struct gkyl_array* nodal_fld = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
  struct gkyl_array *nodal_fld_dev = nodal_fld;


  // Create output field
  struct gkyl_array *funcdg2 = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  struct gkyl_array *funcdg2_dev = funcdg2;

  // Trnasform forward and back
  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&basis, &grid, use_gpu);
  gkyl_nodal_ops_m2n(n2m, basis_on_dev, &grid, &nrange, &local, 1, nodal_fld_dev, funcdg_dev, true);
  gkyl_nodal_ops_n2m(n2m, basis_on_dev, &grid, &nrange, &local, 1, nodal_fld_dev, funcdg2_dev, true);

  check_same(local, basis, funcdg, funcdg2);
  gkyl_grid_sub_array_write(&grid, &local, 0, funcdg2, "proj_func3d_2.gkyl");



  gkyl_array_release(funcdg);
  gkyl_array_release(funcdg2);
  gkyl_array_release(nodal_fld);
  gkyl_eval_on_nodes_release(eon);
  gkyl_nodal_ops_release(n2m);

}




void
test_p1_deflated(){
  // create  grid, ranges, basis
  double lower[] = { 0.0, -1.5 }, upper[] = { 1.5, 1.5 };
  int cells[] = { 8, 16 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);
  struct gkyl_range local, local_ext;
  int nghost[GKYL_MAX_CDIM] = { 1, 1 };
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);
  int poly_order = 1;
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 2, poly_order);
  bool use_gpu = false;
#ifdef GKYL_HAVE_CUDA
  use_gpu = true;
  struct gkyl_basis *basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
  gkyl_cart_modal_serendip_cu_dev(basis_on_dev, 2, poly_order);
#else
  struct gkyl_basis *basis_on_dev = &basis;
#endif

  // Project Initial Function
  struct gkyl_array *funcdg = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_eval_on_nodes *eon = gkyl_eval_on_nodes_new(&grid, &basis, 1, &proj_func, 0);
  gkyl_eval_on_nodes_advance(eon, 0.0, &local, funcdg);
  gkyl_grid_sub_array_write(&grid, &local, 0, funcdg, "proj_func.gkyl");
#ifdef GKYL_HAVE_CUDA
  struct gkyl_array *funcdg_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_array_copy(funcdg_dev, funcdg);
#else
  struct gkyl_array *funcdg_dev = funcdg;
#endif

  // Construct Nodal range, field, and nodal ops object
  int nodes[3] = { 1, 1, 1 };
  for (int d=0; d<grid.ndim; ++d)
    nodes[d] = grid.cells[d] + 1;
  struct gkyl_range nrange;
  gkyl_range_init_from_shape(&nrange, grid.ndim, nodes);
  struct gkyl_array* nodal_fld = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
#ifdef GKYL_HAVE_CUDA
  struct gkyl_array *nodal_fld_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
#else
  struct gkyl_array *nodal_fld_dev = nodal_fld;
#endif
  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&basis, &grid, use_gpu);

  // Create deflated 1d grid, ranges, basis, nodal range
  struct gkyl_rect_grid deflated_grid;
  struct gkyl_basis deflated_basis;
  struct gkyl_range deflated_local, deflated_local_ext;
  double deflated_lower[1] = { grid.lower[0]}, deflated_upper[1] = { grid.upper[0]};
  int deflated_cells[1] = { grid.cells[0] };
  gkyl_rect_grid_init(&deflated_grid, 1, deflated_lower, deflated_upper, deflated_cells);
  int deflated_nghost[GKYL_MAX_CDIM] = { 1 };
  gkyl_create_grid_ranges(&deflated_grid, deflated_nghost, &deflated_local_ext, &deflated_local);
  gkyl_cart_modal_serendip(&deflated_basis, 1, poly_order);
#ifdef GKYL_HAVE_CUDA
  use_gpu = true;
  struct gkyl_basis *deflated_basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
  gkyl_cart_modal_serendip_cu_dev(deflated_basis_on_dev, deflated_grid.ndim, poly_order);
#else
  struct gkyl_basis *deflated_basis_on_dev = &deflated_basis;
#endif
  int deflated_nodes[1];
  for (int d=0; d<deflated_grid.ndim; ++d)
    deflated_nodes[d] = deflated_grid.cells[d] + 1;
  struct gkyl_range deflated_nrange;
  gkyl_range_init_from_shape(&deflated_nrange, deflated_grid.ndim, deflated_nodes);

  // Create deflated modal field 
  struct gkyl_array *deflated_funcdg = gkyl_array_new(GKYL_DOUBLE, deflated_basis.num_basis, local_ext.volume);
#ifdef GKYL_HAVE_CUDA
  struct gkyl_array *deflated_funcdg_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, deflated_basis.num_basis, local_ext.volume);
#else
  struct gkyl_array *deflated_funcdg_dev = deflated_funcdg;
#endif

  // Created deflator objects and 1d nodal ops object
  struct gkyl_deflate_zsurf *deflator_lo = gkyl_deflate_zsurf_new(&basis, &deflated_basis, 0, use_gpu);
  struct gkyl_deflate_zsurf *deflator_up = gkyl_deflate_zsurf_new(&basis, &deflated_basis, 1, use_gpu);
  struct gkyl_nodal_ops *n2m_1d = gkyl_nodal_ops_new(&deflated_basis, &deflated_grid, use_gpu);

  // Loop along last dim (z) and call the 1d deflator to populate the 2d nodal field
  for(int zidx = local.lower[1]; zidx <= local.upper[1]; zidx++){
    gkyl_deflate_zsurf_advance(deflator_lo, zidx, &local, &deflated_local, funcdg_dev, deflated_funcdg_dev, 1);
    gkyl_nodal_ops_m2n_deflated(n2m_1d, deflated_basis_on_dev, &deflated_grid, &nrange, &deflated_nrange, &deflated_local, 1, nodal_fld_dev, deflated_funcdg_dev, zidx-1);
    if (zidx == local.upper[1]) {
      gkyl_deflate_zsurf_advance(deflator_up, zidx, &local, &deflated_local, funcdg_dev, deflated_funcdg_dev, 1);
      gkyl_nodal_ops_m2n_deflated(n2m_1d, deflated_basis_on_dev, &deflated_grid, &nrange, &deflated_nrange, &deflated_local, 1, nodal_fld_dev, deflated_funcdg_dev, zidx);
    }
  }

  // Create output field
  struct gkyl_array *funcdg2 = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
#ifdef GKYL_HAVE_CUDA
  struct gkyl_array *funcdg2_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
#else
  struct gkyl_array *funcdg2_dev = funcdg2;
#endif

  // Transform back to modal
  gkyl_nodal_ops_n2m(n2m, basis_on_dev, &grid, &nrange, &local, 1, nodal_fld_dev, funcdg2_dev, false);
#ifdef GKYL_HAVE_CUDA
  gkyl_array_copy(funcdg2, funcdg2_dev);
#endif
  check_same(local, basis, funcdg, funcdg2);
  gkyl_grid_sub_array_write(&grid, &local, 0, funcdg2, "proj_func2.gkyl");

#ifdef GKYL_HAVE_CUDA
  gkyl_cu_free(basis_on_dev);
  gkyl_array_release(funcdg_dev);
  gkyl_array_release(deflated_funcdg_dev);
  gkyl_array_release(funcdg2_dev);
  gkyl_array_release(nodal_fld_dev);
#endif
  gkyl_array_release(funcdg);
  gkyl_array_release(deflated_funcdg);
  gkyl_array_release(funcdg2);
  gkyl_array_release(nodal_fld);

  gkyl_eval_on_nodes_release(eon);
  gkyl_nodal_ops_release(n2m);
  gkyl_nodal_ops_release(n2m_1d);
  gkyl_deflate_zsurf_release(deflator_lo);
  gkyl_deflate_zsurf_release(deflator_up);

}

void test_p1_deflated_3d(){
  // create  grid, ranges, basis
  double lower[] = { 0.0, -1.5, -2 }, upper[] = { 1.5, 1.5,2 };
  int cells[] = { 8, 16 , 12};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 3, lower, upper, cells);
  struct gkyl_range local, local_ext;
  int nghost[GKYL_MAX_CDIM] = { 1, 1 ,1};
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);
  int poly_order = 1;
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 3, poly_order);
  bool use_gpu = false;
#ifdef GKYL_HAVE_CUDA
  use_gpu = true;
  struct gkyl_basis *basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
  gkyl_cart_modal_serendip_cu_dev(basis_on_dev, 3, poly_order);
#else
  struct gkyl_basis *basis_on_dev = &basis;
#endif

  // Project Initial Function
  struct gkyl_array *funcdg = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_eval_on_nodes *eon = gkyl_eval_on_nodes_new(&grid, &basis, 1, &proj_func3d, 0);
  gkyl_eval_on_nodes_advance(eon, 0.0, &local, funcdg);
  gkyl_grid_sub_array_write(&grid, &local, 0, funcdg, "proj_func.gkyl");
#ifdef GKYL_HAVE_CUDA
  struct gkyl_array *funcdg_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_array_copy(funcdg_dev, funcdg);
#else
  struct gkyl_array *funcdg_dev = funcdg;
#endif

  // Construct Nodal range, field, and nodal ops object
  int nodes[3] = { 1, 1, 1 };
  for (int d=0; d<grid.ndim; ++d)
    nodes[d] = grid.cells[d] + 1;
  struct gkyl_range nrange;
  gkyl_range_init_from_shape(&nrange, grid.ndim, nodes);
  struct gkyl_array* nodal_fld = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
#ifdef GKYL_HAVE_CUDA
  struct gkyl_array *nodal_fld_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
#else
  struct gkyl_array *nodal_fld_dev = nodal_fld;
#endif
  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&basis, &grid, use_gpu);

  // Create deflated 2d grid, ranges, basis, nodal range
  struct gkyl_rect_grid deflated_grid;
  struct gkyl_basis deflated_basis;
  struct gkyl_range deflated_local, deflated_local_ext;
  double deflated_lower[2] = { grid.lower[0], grid.lower[1] };
  double deflated_upper[2] = { grid.upper[0], grid.upper[1] };
  int deflated_cells[2] = { grid.cells[0], grid.cells[1] };
  gkyl_rect_grid_init(&deflated_grid, 2, deflated_lower, deflated_upper, deflated_cells);
  int deflated_nghost[GKYL_MAX_CDIM] = { 1, 1 };
  gkyl_create_grid_ranges(&deflated_grid, deflated_nghost, &deflated_local_ext, &deflated_local);
  gkyl_cart_modal_serendip(&deflated_basis, 2, poly_order);
#ifdef GKYL_HAVE_CUDA
  use_gpu = true;
  struct gkyl_basis *deflated_basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
  gkyl_cart_modal_serendip_cu_dev(deflated_basis_on_dev, deflated_grid.ndim, poly_order);
#else
  struct gkyl_basis *deflated_basis_on_dev = &deflated_basis;
#endif
  int deflated_nodes[2];
  for (int d=0; d<deflated_grid.ndim; ++d)
    deflated_nodes[d] = deflated_grid.cells[d] + 1;
  struct gkyl_range deflated_nrange;
  gkyl_range_init_from_shape(&deflated_nrange, deflated_grid.ndim, deflated_nodes);

  // Create deflated modal field 
  struct gkyl_array *deflated_funcdg = gkyl_array_new(GKYL_DOUBLE, deflated_basis.num_basis, local_ext.volume);
#ifdef GKYL_HAVE_CUDA
  struct gkyl_array *deflated_funcdg_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, deflated_basis.num_basis, local_ext.volume);
#else
  struct gkyl_array *deflated_funcdg_dev = deflated_funcdg;
#endif

  // Created deflator objects and 1d nodal ops object
  struct gkyl_deflate_zsurf *deflator_lo = gkyl_deflate_zsurf_new(&basis, &deflated_basis, 0, use_gpu);
  struct gkyl_deflate_zsurf *deflator_up = gkyl_deflate_zsurf_new(&basis, &deflated_basis, 1, use_gpu);
  struct gkyl_nodal_ops *n2m_2d = gkyl_nodal_ops_new(&deflated_basis, &deflated_grid, use_gpu);

  // Loop along last dim (z) and call the 1d deflator to populate the 2d nodal field
  for(int zidx = local.lower[2]; zidx <= local.upper[2]; zidx++){
    gkyl_deflate_zsurf_advance(deflator_lo, zidx, &local, &deflated_local, funcdg_dev, deflated_funcdg_dev, 1);
    gkyl_nodal_ops_m2n_deflated(n2m_2d, deflated_basis_on_dev, &deflated_grid, &nrange, &deflated_nrange, &deflated_local, 1, nodal_fld_dev, deflated_funcdg_dev, zidx-1);
    if (zidx == local.upper[2]) {
      gkyl_deflate_zsurf_advance(deflator_up, zidx, &local, &deflated_local, funcdg_dev, deflated_funcdg_dev, 1);
      gkyl_nodal_ops_m2n_deflated(n2m_2d, deflated_basis_on_dev, &deflated_grid, &nrange, &deflated_nrange, &deflated_local, 1, nodal_fld_dev, deflated_funcdg_dev, zidx);
    }
  }

  // Create output field
  struct gkyl_array *funcdg2 = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
#ifdef GKYL_HAVE_CUDA
  struct gkyl_array *funcdg2_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
#else
  struct gkyl_array *funcdg2_dev = funcdg2;
#endif

  // Transform back to modal
  gkyl_nodal_ops_n2m(n2m, basis_on_dev, &grid, &nrange, &local, 1, nodal_fld_dev, funcdg2_dev, false);
#ifdef GKYL_HAVE_CUDA
  gkyl_array_copy(funcdg2, funcdg2_dev);
#endif
  check_same(local, basis, funcdg, funcdg2);
  gkyl_grid_sub_array_write(&grid, &local, 0, funcdg2, "proj_func2.gkyl");

#ifdef GKYL_HAVE_CUDA
  gkyl_cu_free(basis_on_dev);
  gkyl_array_release(funcdg_dev);
  gkyl_array_release(deflated_funcdg_dev);
  gkyl_array_release(funcdg2_dev);
  gkyl_array_release(nodal_fld_dev);
#endif
  gkyl_array_release(funcdg);
  gkyl_array_release(deflated_funcdg);
  gkyl_array_release(funcdg2);
  gkyl_array_release(nodal_fld);

  gkyl_eval_on_nodes_release(eon);
  gkyl_nodal_ops_release(n2m);
  gkyl_nodal_ops_release(n2m_2d);
  gkyl_deflate_zsurf_release(deflator_lo);
  gkyl_deflate_zsurf_release(deflator_up);

}


void
test_p2_btype(enum gkyl_basis_type basis_type){
  // create  grid
  double lower[] = { 0.0, -1.5 }, upper[] = { 1.5, 1.5 };
  // as ellipitical surfaces are exact, we only need 1 cell in each
  // direction
  int cells[] = { 8, 16 };


  int poly_order = 2;
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  //  ranges
  struct gkyl_range local, local_ext;
  int nghost[GKYL_MAX_CDIM] = { 1, 1 };
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);

  //  basis function
  struct gkyl_basis basis;
  if(basis_type == GKYL_BASIS_MODAL_SERENDIPITY)
    gkyl_cart_modal_serendip(&basis, 2, poly_order);
  if(basis_type == GKYL_BASIS_MODAL_TENSOR)
    gkyl_cart_modal_tensor(&basis, 2, poly_order);

  struct gkyl_array *funcdg = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  
  gkyl_eval_on_nodes *eon = gkyl_eval_on_nodes_new(&grid, &basis, 1, &proj_func, 0);
  gkyl_eval_on_nodes_advance(eon, 0.0, &local, funcdg);
  gkyl_eval_on_nodes_release(eon);
  gkyl_grid_sub_array_write(&grid, &local, 0, funcdg, "proj_func.gkyl");

  int nodes[3] = { 1, 1, 1 };
                   
  for (int d=0; d<grid.ndim; ++d)
    nodes[d] = 2*(grid.cells[d]) + 1;

  struct gkyl_range nrange;
  gkyl_range_init_from_shape(&nrange, grid.ndim, nodes);


  struct gkyl_array* nodal_fld = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&basis, &grid, false);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange, &local, 1, nodal_fld, funcdg, false);

  struct gkyl_array *funcdg2 = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_nodal_ops_n2m(n2m, &basis, &grid, &nrange, &local, 1, nodal_fld, funcdg2, false);
  gkyl_nodal_ops_release(n2m);
  check_same(local, basis, funcdg, funcdg2);
  gkyl_grid_sub_array_write(&grid, &local, 0, funcdg2, "proj_func2.gkyl");

  gkyl_array_release(funcdg);
  gkyl_array_release(funcdg2);
  gkyl_array_release(nodal_fld);
  
}

void
test_p2_ser()
{
  return test_p2_btype(GKYL_BASIS_MODAL_SERENDIPITY);
}

void
test_p2_tensor()
{
  return test_p2_btype(GKYL_BASIS_MODAL_TENSOR);
}


TEST_LIST = {
  { "test_p1_interior_2x", test_p1_interior_2x},
  { "test_p1_interior_3x", test_p1_interior_3x},
  {"test_p1_2x", test_p1_2x},
  { "test_p1_3x", test_p1_3x},
  { "test_p1_deflated", test_p1_deflated},
  { "test_p1_deflated_3d", test_p1_deflated_3d},
  { "test_p2_ser", test_p2_ser},
  { "test_p2_tensor", test_p2_tensor},
  { NULL, NULL },
};
