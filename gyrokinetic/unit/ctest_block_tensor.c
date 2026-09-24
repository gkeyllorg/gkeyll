
#include <acutest.h>
#include <math.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_ops.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <gkyl_basis.h>
#include <gkyl_eval_on_nodes.h>

#include <gkyl_bc_block_tensor.h>
#include <gkyl_bc_block_tensor_priv.h>

void
proj_one(double t, const double *xn, double *fout, void *ctx)
{
  fout[0] = 1;
  fout[1] = 0;
  fout[2] = 0;
  fout[3] = 0;
  fout[4] = 1;
  fout[5] = 0;
  fout[6] = 0;
  fout[7] = 0;
  fout[8] = 1;
}

void
proj_cyl(double t, const double *xn, double *fout, void *ctx)
{
  double r = xn[0];
  double theta = xn[1];

  fout[0] = cos(theta);
  fout[1] = 0;
  fout[2] = sin(theta);

  fout[3] = 0;
  fout[4] = 1;
  fout[5] = 0;

  fout[6] = -r*sin(theta);
  fout[7] = 0;
  fout[8] = r*cos(theta);
}



void test_cartesian_2x_onecell()
{
  double lower[] = { -1, -1 }, upper[] = { 1, 1 };
  int cells[] = { 8, 16 };
  int dim = sizeof(lower)/sizeof(lower[0]);
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  //ranges
  struct gkyl_range local, local_ext;
  int nghost[GKYL_MAX_CDIM] = { 1, 1 };
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);
  // basis function
  int poly_order = 1;
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  // Create cartesian tangent vectors and duals
  struct gkyl_array *dxdz= gkyl_array_new(GKYL_DOUBLE, 9*basis.num_basis, local_ext.volume);
  struct gkyl_array *dzdx= gkyl_array_new(GKYL_DOUBLE, 9*basis.num_basis, local_ext.volume);
  gkyl_eval_on_nodes *proj = gkyl_eval_on_nodes_new(&grid, &basis, 9, &proj_one, 0);
  gkyl_eval_on_nodes_advance(proj, 0.0, &local, dxdz);
  gkyl_eval_on_nodes_advance(proj, 0.0, &local, dzdx);
  gkyl_eval_on_nodes_release(proj);
  gkyl_grid_sub_array_write(&grid, &local, 0, dxdz, "dxdz.gkyl");
  gkyl_grid_sub_array_write(&grid, &local, 0, dzdx, "dzdx.gkyl");

  struct bc_block_tensor *bt = gkyl_bc_block_tensor_new(&grid, &local, &local_ext, &basis, false);

  int edge1 = 1; //upper edge
  int edge2 = 0; //lower edge
  int dir = 1; // second direction
  int idx1[GKYL_MAX_DIM] = {4,16};
  int idx2[GKYL_MAX_DIM] = {4,1};
  long loc1 = gkyl_range_idx(&local, idx1);
  long loc2 = gkyl_range_idx(&local, idx2);
  const double *dxdz_i = gkyl_array_fetch(dxdz, loc1);
  const double *dzdx_j = gkyl_array_fetch(dzdx, loc2);

  int num_nodes = 2;
  double tj_i[dim*dim*num_nodes];

  calc_tensor(bt, dir, edge1, edge2, dzdx_j, dxdz_i, tj_i);
  for(int n = 0; n<num_nodes; n++)
    for (int i = 0; i < 4; i++){
      printf("tj_i[%d][%d] = %g\n", i,n, tj_i[i*num_nodes + n]);
    }

  printf("\n OR \n");
  for (int i = 0; i < 8; i++)
    printf("tj_i[%d] = %g\n", i, tj_i[i]);

  gkyl_array_release(dxdz);
  gkyl_array_release(dzdx);
  gkyl_bc_block_tensor_release(bt);
}

void test_cartesian_2x_z()
{
  // Block 2 grid
  double lower2[] = { -1.0, 0.0 }, upper2[] = { 1.0, 1.0 };
  int cells2[] = { 8, 16 };
  int dim = sizeof(lower2)/sizeof(lower2[0]);
  struct gkyl_rect_grid grid2;
  gkyl_rect_grid_init(&grid2, dim, lower2, upper2, cells2);
  struct gkyl_range local2, local_ext2;
  int nghost[GKYL_MAX_CDIM] = { 1, 1 };
  gkyl_create_grid_ranges(&grid2, nghost, &local_ext2, &local2);

  // Block 1 grid
  double lower1[] = { -1.0, -1 }, upper1[] = { 1.0, 0.0 };
  int cells1[] = { 8, 16 };
  struct gkyl_rect_grid grid1;
  gkyl_rect_grid_init(&grid1, dim, lower1, upper1, cells1);
  struct gkyl_range local1, local_ext1;
  gkyl_create_grid_ranges(&grid1, nghost, &local_ext1, &local1);



  // Common basis
  int poly_order = 1;
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  // Block 2 duals
  struct gkyl_array *dzdx2 = gkyl_array_new(GKYL_DOUBLE, 9*basis.num_basis, local_ext2.volume);
  gkyl_eval_on_nodes *proj2 = gkyl_eval_on_nodes_new(&grid2, &basis, 9, &proj_one, 0);
  gkyl_eval_on_nodes_advance(proj2, 0.0, &local2, dzdx2);
  gkyl_eval_on_nodes_release(proj2);
  gkyl_grid_sub_array_write(&grid2, &local2, 0, dzdx2, "dzdx2.gkyl");

  // Block 1 tangents
  struct gkyl_array *dxdz1 = gkyl_array_new(GKYL_DOUBLE, 9*basis.num_basis, local_ext1.volume);
  gkyl_eval_on_nodes *proj1 = gkyl_eval_on_nodes_new(&grid1, &basis, 9, &proj_one, 0);
  gkyl_eval_on_nodes_advance(proj1, 0.0, &local1, dxdz1);
  gkyl_eval_on_nodes_release(proj1);
  gkyl_grid_sub_array_write(&grid1, &local1, 0, dxdz1, "dxdz1.gkyl");



  struct bc_block_tensor *bt = gkyl_bc_block_tensor_new(&grid2, &local2, &local_ext2, &basis, false);
  int edge1 = 1; //upper edge
  int edge2 = 0; //lower edge
  int dir = 1; // second direction
  gkyl_bc_block_tensor_advance(bt, dir, edge1, edge2, dxdz1, dzdx2, &local1, &local2);
  gkyl_grid_sub_array_write(&grid2, &local2, 0,bt->tensor, "tji.gkyl");

  gkyl_array_release(dxdz1);
  gkyl_array_release(dzdx2);
  gkyl_bc_block_tensor_release(bt);
}

void test_cartesian_2x_x()
{
  // Block 2 grid
  double lower2[] = { -1.0, 0.0 }, upper2[] = { 1.0, 1.0 };
  int cells2[] = { 8, 16 };
  int dim = sizeof(lower2)/sizeof(lower2[0]);
  struct gkyl_rect_grid grid2;
  gkyl_rect_grid_init(&grid2, dim, lower2, upper2, cells2);
  struct gkyl_range local2, local_ext2;
  int nghost[GKYL_MAX_CDIM] = { 1, 1 };
  gkyl_create_grid_ranges(&grid2, nghost, &local_ext2, &local2);

  // Block 1 grid
  double lower1[] = { -1.0, -1 }, upper1[] = { 0.0, 1.0 };
  int cells1[] = { 8, 16 };
  struct gkyl_rect_grid grid1;
  gkyl_rect_grid_init(&grid1, dim, lower1, upper1, cells1);
  struct gkyl_range local1, local_ext1;
  gkyl_create_grid_ranges(&grid1, nghost, &local_ext1, &local1);



  // Common basis
  int poly_order = 1;
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  // Block 2 duals
  struct gkyl_array *dzdx2 = gkyl_array_new(GKYL_DOUBLE, 9*basis.num_basis, local_ext2.volume);
  gkyl_eval_on_nodes *proj2 = gkyl_eval_on_nodes_new(&grid2, &basis, 9, &proj_one, 0);
  gkyl_eval_on_nodes_advance(proj2, 0.0, &local2, dzdx2);
  gkyl_eval_on_nodes_release(proj2);
  gkyl_grid_sub_array_write(&grid2, &local2, 0, dzdx2, "dzdx2.gkyl");

  // Block 1 tangents
  struct gkyl_array *dxdz1 = gkyl_array_new(GKYL_DOUBLE, 9*basis.num_basis, local_ext1.volume);
  gkyl_eval_on_nodes *proj1 = gkyl_eval_on_nodes_new(&grid1, &basis, 9, &proj_one, 0);
  gkyl_eval_on_nodes_advance(proj1, 0.0, &local1, dxdz1);
  gkyl_eval_on_nodes_release(proj1);
  gkyl_grid_sub_array_write(&grid1, &local1, 0, dxdz1, "dxdz1.gkyl");



  struct bc_block_tensor *bt = gkyl_bc_block_tensor_new(&grid2, &local2, &local_ext2, &basis, false);
  int edge1 = 1; //upper edge
  int edge2 = 0; //lower edge
  int dir = 0; // first direction
  gkyl_bc_block_tensor_advance(bt, dir, edge1, edge2, dxdz1, dzdx2, &local1, &local2);
  gkyl_grid_sub_array_write(&grid2, &local2, 0, bt->tensor, "tji.gkyl");

  gkyl_array_release(dxdz1);
  gkyl_array_release(dzdx2);
  gkyl_bc_block_tensor_release(bt);
}

// Block 2 cartesian and block 1 is cylindrical
// Make it such that the blocks line up
void test_cyl_cart_2x_z()
{
  // Block 2 grid
  double lower2[] = { 0.5, -1.0 }, upper2[] = { 1.0, 0.0 };
  int cells2[] = { 8, 16 };
  int dim = sizeof(lower2)/sizeof(lower2[0]);
  struct gkyl_rect_grid grid2;
  gkyl_rect_grid_init(&grid2, dim, lower2, upper2, cells2);
  struct gkyl_range local2, local_ext2;
  int nghost[GKYL_MAX_CDIM] = { 1, 1 };
  gkyl_create_grid_ranges(&grid2, nghost, &local_ext2, &local2);

  // Block 1 grid
  double lower1[] = { 0.5, 0.0 }, upper1[] = { 1.0, M_PI/2 };
  int cells1[] = { 8, 16 };
  struct gkyl_rect_grid grid1;
  gkyl_rect_grid_init(&grid1, dim, lower1, upper1, cells1);
  struct gkyl_range local1, local_ext1;
  gkyl_create_grid_ranges(&grid1, nghost, &local_ext1, &local1);



  // Common basis
  int poly_order = 1;
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  // Block 2 duals
  struct gkyl_array *dzdx2 = gkyl_array_new(GKYL_DOUBLE, 9*basis.num_basis, local_ext2.volume);
  gkyl_eval_on_nodes *proj2 = gkyl_eval_on_nodes_new(&grid2, &basis, 9, &proj_one, 0);
  gkyl_eval_on_nodes_advance(proj2, 0.0, &local2, dzdx2);
  gkyl_eval_on_nodes_release(proj2);
  gkyl_grid_sub_array_write(&grid2, &local2, 0, dzdx2, "dzdx2.gkyl");

  // Block 1 tangents
  struct gkyl_array *dxdz1 = gkyl_array_new(GKYL_DOUBLE, 9*basis.num_basis, local_ext1.volume);
  gkyl_eval_on_nodes *proj1 = gkyl_eval_on_nodes_new(&grid1, &basis, 9, &proj_cyl, 0);
  gkyl_eval_on_nodes_advance(proj1, 0.0, &local1, dxdz1);
  gkyl_eval_on_nodes_release(proj1);
  gkyl_grid_sub_array_write(&grid1, &local1, 0, dxdz1, "dxdz1.gkyl");



  struct bc_block_tensor *bt = gkyl_bc_block_tensor_new(&grid2, &local2, &local_ext2, &basis, false);
  int edge1 = 0; // lower edge
  int edge2 = 1; // upper edge
  int dir = 1; // second direction
  gkyl_bc_block_tensor_advance(bt, dir, edge1, edge2, dxdz1, dzdx2, &local1, &local2);
  gkyl_grid_sub_array_write(&grid2, &local2, 0, bt->tensor, "tji.gkyl");

  gkyl_array_release(dxdz1);
  gkyl_array_release(dzdx2);
  gkyl_bc_block_tensor_release(bt);
}


void test_cartesian_3x_onecell()
{
  double lower[] = { -1, -1 ,-1}, upper[] = { 1, 1 ,1};
  int cells[] = { 8, 16,4 };
  int dim = sizeof(lower)/sizeof(lower[0]);
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  //ranges
  struct gkyl_range local, local_ext;
  int nghost[GKYL_MAX_CDIM] = { 1, 1 ,1};
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);
  // basis function
  int poly_order = 1;
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  // Create cartesian tangent vectors and duals
  struct gkyl_array *dxdz= gkyl_array_new(GKYL_DOUBLE, 9*basis.num_basis, local_ext.volume);
  struct gkyl_array *dzdx= gkyl_array_new(GKYL_DOUBLE, 9*basis.num_basis, local_ext.volume);
  gkyl_eval_on_nodes *proj = gkyl_eval_on_nodes_new(&grid, &basis, 9, &proj_one, 0);
  gkyl_eval_on_nodes_advance(proj, 0.0, &local, dxdz);
  gkyl_eval_on_nodes_advance(proj, 0.0, &local, dzdx);
  gkyl_eval_on_nodes_release(proj);
  gkyl_grid_sub_array_write(&grid, &local, 0, dxdz, "dxdz.gkyl");
  gkyl_grid_sub_array_write(&grid, &local, 0, dzdx, "dzdx.gkyl");

  struct bc_block_tensor *bt = gkyl_bc_block_tensor_new(&grid, &local, &local_ext, &basis, false);

  int edge1 = 1; //upper edge
  int edge2 = 0; //lower edge
  int dir = 2; // third direction
  int idx1[3] = {4,8,4};
  int idx2[3] = {4,8,4};
  long loc1 = gkyl_range_idx(&local, idx1);
  long loc2 = gkyl_range_idx(&local, idx2);
  const double *dxdz_i = gkyl_array_fetch(dxdz, loc1);
  const double *dzdx_j = gkyl_array_fetch(dzdx, loc2);

  int num_nodes = 4;
  double tj_i[dim*dim*num_nodes];

  calc_tensor(bt, dir, edge1, edge2, dzdx_j, dxdz_i, tj_i);
  for(int n = 0; n<num_nodes; n++)
    for (int i = 0; i < 9; i++){
      printf("tj_i[%d][%d] = %g\n", i,n, tj_i[i*num_nodes + n]);
    }

  //printf("\n OR \n");
  //for (int i = 0; i < 36; i++)
  //  printf("tj_i[%d] = %g\n", i, tj_i[i]);
  gkyl_array_release(dxdz);
  gkyl_array_release(dzdx);
  gkyl_bc_block_tensor_release(bt);


}



TEST_LIST = {
  //{ "test_cartesian_2x_z", test_cartesian_2x_z},
  //{ "test_cartesian_2x_x", test_cartesian_2x_x},
  { "test_cyl_cart_2x_z", test_cyl_cart_2x_z},
  { NULL, NULL },
};
