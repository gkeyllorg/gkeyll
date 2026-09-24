/**
 * This code tests the skin_surf_from_ghost (ssfg) updater. It projects a field on a basis, sets ghost cell values, and then
 * copies the field values from ghost cells to skin cells. The test checks that the skin cell value meets the ghost cell value.
 */

#include "gkyl_array.h"
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_rio.h>
#include <acutest.h>
#include <gkyl_skin_surf_from_ghost.h>  

// Evaluate the projection of the modal representation inside a cell from 1x to 3x
double eval_f(const double *phi, double x, double y, double z, int cdim);

// Function to allocate a gkyl array, zero-initialized, on CPU or GPU
static struct gkyl_array* mkarr(bool on_gpu, long nc, long size) {
  return on_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
                : gkyl_array_new(GKYL_DOUBLE, nc, size);
}

// Analytical field evaluation function, setting field to zero
void eval_field(double t, const double *xn, double* restrict fout, void *ctx) {
  fout[0] = 0.0;
}

// Function to set up and test the ghost-to-skin surface copy updater
void test_ssfg(int cdim, int poly_order, bool use_gpu, enum gkyl_edge_loc edge, int dir, bool control) { 
  double lower[cdim], upper[cdim];
  int cells[cdim];
  switch(cdim) {
    case 3:
      lower[0] = 0.0; lower[1] = 0.0; lower[2] = 0.0;
      upper[0] = 1.0; upper[1] = 1.0; upper[2] = 1.0;
      cells[0] = 2;   cells[1] = 4;   cells[2] = 6;
      break;
    case 2:
      lower[0] = 0.0; lower[1] = 0.0;
      upper[0] = 1.0; upper[1] = 1.0;
      cells[0] = 4;   cells[1] = 6;
      break;
    case 1:
      lower[0] = 0.0;
      upper[0] = 1.0;
      cells[0] = 6;
      break;
    default:
      fprintf(stderr, "Invalid cdim value: %d\n", cdim);
      exit(1);
  }
  const int ndim = sizeof(cells)/sizeof(cells[0]);

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);
  double dx = grid.dx[dir];

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, cdim, poly_order);

  int ghost[GKYL_MAX_CDIM] = {1, 1, 1};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct gkyl_range skin_r, ghost_r;
  gkyl_skin_ghost_ranges(&skin_r, &ghost_r, dir, edge, &local_ext, ghost);

  struct gkyl_array *field_ho, *field;
  field = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  field_ho = use_gpu ? mkarr(false, field->ncomp, field->size)
                     : gkyl_array_acquire(field);

  gkyl_proj_on_basis *proj_field = gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, eval_field, NULL);
  gkyl_proj_on_basis_advance(proj_field, 0.0, &local, field_ho);
  gkyl_proj_on_basis_release(proj_field);
  
  // Set ghost cell modal values
  double fg0 = sqrt(pow(2,cdim));
  double fg1 = 1.0;
  double fg2 = 2.0;
  double fg3 = 3.0;
  double fg4 = 4.0;
  double fg5 = 5.0;
  double fg6 = 6.0;
  double fg7 = 7.0;
  struct gkyl_range_iter iter_ghost;
  gkyl_range_iter_init(&iter_ghost, &ghost_r);
  while (gkyl_range_iter_next(&iter_ghost)) {
    long linidx = gkyl_range_idx(&ghost_r, iter_ghost.idx);
    double *f_i = gkyl_array_fetch(field_ho, linidx);
    f_i[0] = fg0;
    f_i[1] = fg1;
    if (cdim > 1){
      f_i[2] = fg2;
      f_i[3] = fg3;
    }
    if (cdim > 2){
      f_i[4] = fg4;
      f_i[5] = fg5;
      f_i[6] = fg6;
      f_i[7] = fg7;
    }
  }

  // Copy field values to the GPU if necessary
  gkyl_array_copy(field, field_ho);

  // Initialize the skin-surf updater and call it if control is false
  gkyl_skin_surf_from_ghost* up = gkyl_skin_surf_from_ghost_new(dir, edge, basis, &skin_r, &ghost_r, use_gpu);
  if (!control) // to test identity operation
    gkyl_skin_surf_from_ghost_advance(up, field);
  gkyl_skin_surf_from_ghost_release(up);

  // Copy field values back to host for checking
  gkyl_array_copy(field_ho, field);

  // Check that field values in the skin cells meet the expected result
  int gidx[GKYL_MAX_DIM];
  struct gkyl_range_iter iter_skin;
  gkyl_range_iter_init(&iter_skin, &skin_r);
  while (gkyl_range_iter_next(&iter_skin)) {
    // get skin cell modal values
    long skin_linidx = gkyl_range_idx(&skin_r, iter_skin.idx);
    const double *fskin = gkyl_array_cfetch(field_ho, skin_linidx);
    // Get ghost cell corresponding to skin cell
    gkyl_copy_int_arr(ndim, iter_skin.idx, gidx);
    gidx[dir] = edge == GKYL_LOWER_EDGE? iter_skin.idx[dir]-1 : iter_skin.idx[dir]+1; 
    long ghost_linidx = gkyl_range_idx(&ghost_r, gidx);
    const double *fghost = (const double*) gkyl_array_cfetch(field_ho, ghost_linidx);

    double xs = 0.0, ys = 0.0, zs = 0.0, xg = 0.0, yg = 0.0, zg = 0.0;
    char dir_c;
    double edge_skin = edge == GKYL_UPPER_EDGE? 1.0 : -1.0;
    switch(dir) {
      case 0:
        xs = edge_skin; 
        xg = -edge_skin;
        dir_c = 'x';
        break;
      case 1:
        ys = edge_skin; 
        yg = -edge_skin; 
        dir_c = 'y';
        break;
      case 2:
        zs = edge_skin; 
        zg = -edge_skin; 
        dir_c = 'z';
        break;
    }
    double vskin  = eval_f(fskin,xs,ys,zs,cdim);
    double vghost = eval_f(fghost,xg,yg,zg,cdim);
    double check_val;

    // Test values based on control mode
    check_val = control? 0.0 : vghost;
    TEST_CHECK(gkyl_compare(vskin, check_val, 1e-14));
  }
  gkyl_array_release(field);
  gkyl_array_release(field_ho);
}

// Tests for 3D case on CPU
void test_ssfg_ho() {
  bool use_gpu    = false;
  int  poly_order = 1;
  for (int cdim = 1; cdim <= 3; cdim++) {
    for (int control = 0; control <= 0; control++) {
      for (int dir = 0; dir <= cdim-1; dir++) {
        test_ssfg(cdim, poly_order, use_gpu, GKYL_UPPER_EDGE, dir, control == 1);
        test_ssfg(cdim, poly_order, use_gpu, GKYL_LOWER_EDGE, dir, control == 1);
      }
    }
  }
}

// Tests for 3D case on GPU (if available)
void test_ssfg_dev() {
  int  poly_order = 1;
  bool use_gpu    = true;
  for (int cdim = 1; cdim <= 3; cdim++) {
    for (int control = 0; control <= 0; control++) {
      for (int dir = 0; dir <= cdim-1; dir++) {
        test_ssfg(cdim, poly_order, use_gpu, GKYL_UPPER_EDGE, dir, control==1);
        test_ssfg(cdim, poly_order, use_gpu, GKYL_LOWER_EDGE, dir, control==1);
      }
    }
  }
}

// Evaluate the projection of the modal representation inside a cell from 1x to 3x
double eval_f(const double *phi, const double x, const double y, const double z, const int cdim) {
  double sqrt2 = sqrt(2.0);
  double sqrt3 = sqrt(3.0);
  double denom = pow(2.0, 1.5);
  switch (cdim) {
    case 1:
      return (sqrt3 * phi[1] * x) / sqrt2 + phi[0] / sqrt2;
    case 2:
      return (3 * phi[3] * x * y) / 2.0 + (sqrt3 * phi[2] * y) / 2.0
             + (sqrt3 * phi[1] * x) / 2.0 + phi[0] / 2.0;
    case 3:
      return (pow(3.0, 1.5) * phi[7] * x * y * z) / denom + (3 * phi[6] * y * z) / denom
             + (3 * phi[5] * x * z) / denom + (sqrt3 * phi[3] * z) / denom
             + (3 * phi[4] * x * y) / denom + (sqrt3 * phi[2] * y) / denom
             + (sqrt3 * phi[1] * x) / denom + phi[0] / denom;
    default:
      fprintf(stderr, "Invalid cdim value: %d. Must be 0, 1, or 2.\n", cdim);
      return 0.0;
  }
}

// List of tests for the test framework
TEST_LIST = {
  { "test_ssfg_ho", test_ssfg_ho },
#ifdef GKYL_HAVE_CUDA
  { "test_ssfg_dev", test_ssfg_dev },
#endif
  { NULL, NULL },
};