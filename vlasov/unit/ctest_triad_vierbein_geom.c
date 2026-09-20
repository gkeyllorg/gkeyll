#include <acutest.h>

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov_triad_geom.h>
#include <gkyl_vlasov_triad_geom_priv.h>
#include <gkyl_util.h>

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

void
test_triad_math_1v()
{
  // Test: Can the correct bracket be constructed for:
  // 1D streched coordinates: x(z^1) = (z^1)^2 \sigma_x
  int vdim = 1;

  // Accuracy of the tests
  double eps = 1e-12;
  double vierbein[1];

  // Randomly choosen coordinates
  double z1 = 3.2;

  // Assign covaraint tangent basis 
  vierbein[0] = 2*z1;

  // 1. Test kernel_vierbein_inv_1v()
  double vierbein_inv[1];
  kernel_vierbein_inv_1v(vierbein, vierbein_inv);

  // vierbein_inv only has symmetric entries as well
  TEST_CHECK(gkyl_compare_double(vierbein_inv[0], 1.0 / (2.0 * z1), eps));

  // 2. Test kernel_conf_poisson_tensor_1v()
  double conf_poisson_tensor[1];

  // Assign the triad basis gradient
  double vierbein_gradient[1];

  // gradient in z1
  vierbein_gradient[0] = 2.0;

  kernel_conf_poisson_tensor_vierbein_1v(vierbein_inv,
    vierbein_gradient, conf_poisson_tensor);

  // Pi_{xx} block
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[0], 1.0 / (2*z1), eps));

}

void eval_vierbein_1v(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = 1.0;
}

void eval_vierbein_gradient_1v(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = 0.0;
}

void
test_triad_1x1v_flat_conf(int poly_order)
{

  double lower[] = {0.1, -1.0}, upper[] = {1.0, 1.0};
  int cells[] = {2, 2};
  int vdim = 1, cdim = 1;
  int ndim = cdim+vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};
  double velLower[] = {lower[1]}, velUpper[] = {upper[1]};
  int velCells[] = {cells[1]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

    struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  // velocity range
  int velGhost[] = {0, 0, 0};
  struct gkyl_range velLocal, velLocal_ext; 
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);

  int confGhost[] = { 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = { confGhost[0], 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Construct the input map
  struct gkyl_vlasov_triad_geom_inp inp_triad_geom = { 0 };
  inp_triad_geom.use_vierbein = true;
  inp_triad_geom.use_preset_geom = false;
  inp_triad_geom.eval_vierbein = eval_vierbein_1v; 
  inp_triad_geom.eval_vierbein_gradient = eval_vierbein_gradient_1v; 
  inp_triad_geom.eval_vierbein_ctx = 0; 
  inp_triad_geom.eval_vierbein_gradient_ctx = 0; 

  // Make the memory for arrays (modal)
  struct gkyl_array *conf_poisson_tensor; // Configuration space Poisson tensor representation

  // Size of the PT
  int num_pt_indices[3] = { 1 , 6, 18 }; 

  // Allocate arrays for covariant tangent basis 
  conf_poisson_tensor = mkarr(confBasis.num_basis*num_pt_indices[vdim-1], confLocal_ext.volume);

  // Construct the Geometry for this configuration
  gkyl_vlasov_triad_geom_new(&confGrid, &confLocal, confBasis, 
    &grid, &local, basis, inp_triad_geom, conf_poisson_tensor);

  // Iterate over the grid, conf space, checking output
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &confLocal);
  while (gkyl_range_iter_next(&iter)) {
    const double *conf_poisson_tensor_d = gkyl_array_cfetch(conf_poisson_tensor, gkyl_range_idx(&confLocal, iter.idx));
    
    int NC = confBasis.num_basis;
    TEST_CHECK( gkyl_compare_double(conf_poisson_tensor_d[0], sqrt(2.0), 1e-12) );
    TEST_CHECK( gkyl_compare_double(conf_poisson_tensor_d[1], sqrt(0.0), 1e-12) );
    TEST_CHECK( gkyl_compare_double(conf_poisson_tensor_d[2], sqrt(0.0), 1e-12) );
  }


  // Release the memory
  gkyl_array_release(conf_poisson_tensor);
}




void eval_vierbein_2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = 1.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = 1.0;
}

void eval_vierbein_gradient_2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = 0.0;

  fout[4] = 0.0;
  fout[5] = 0.0;
  fout[6] = 0.0;
  fout[7] = 0.0;
}

void
test_triad_1x2v_flat_conf(int poly_order)
{

  double lower[] = {0.1, -1.0, -1.0}, upper[] = {1.0, 1.0, 1.0};
  int cells[] = {2, 2, 2};
  int vdim = 2, cdim = 1;
  int ndim = cdim+vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};
  double velLower[] = {lower[1], lower[2]}, velUpper[] = {upper[1], upper[2]};
  int velCells[] = {cells[1], cells[2]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

    struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  // velocity range
  int velGhost[] = { 0, 0 };
  struct gkyl_range velLocal, velLocal_ext; 
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);

  int confGhost[] = { 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = { confGhost[0], 0, 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Construct the input map
  struct gkyl_vlasov_triad_geom_inp inp_triad_geom = { 0 };
  inp_triad_geom.use_vierbein = true;
  inp_triad_geom.use_preset_geom = false;
  inp_triad_geom.eval_vierbein = eval_vierbein_2v; 
  inp_triad_geom.eval_vierbein_gradient = eval_vierbein_gradient_2v; 
  inp_triad_geom.eval_vierbein_ctx = 0; 
  inp_triad_geom.eval_vierbein_gradient_ctx = 0; 

  // Make the memory for arrays (modal)
  struct gkyl_array *conf_poisson_tensor; // Configuration space Poisson tensor representation

  // Size of the PT
  int num_pt_indices[3] = { 1 , 6, 18 }; 

  // Allocate arrays for covariant tangent basis 
  conf_poisson_tensor = mkarr(confBasis.num_basis*num_pt_indices[vdim-1], confLocal_ext.volume);

  // Construct the Geometry for this configuration
  gkyl_vlasov_triad_geom_new(&confGrid, &confLocal, confBasis, 
    &grid, &local, basis, inp_triad_geom, conf_poisson_tensor);

  // Iterate over the grid, conf space, checking output
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &confLocal);
  while (gkyl_range_iter_next(&iter)) {
    const double *conf_poisson_tensor_d = gkyl_array_cfetch(conf_poisson_tensor, gkyl_range_idx(&confLocal, iter.idx));
    
    int NC = confBasis.num_basis;

    // index size of the symmetric vierbein
    int n_sym = vdim*(vdim+1)/2; 

    /* Shape of conf_poisson_tensor non-zero elements.
    // xx - components are not used by 1x2v, but might still be saved/ non-zero
    // -- - compoentns are zero/implied by symmetry
    // Used - components actualled used by the kernels 
    // 
    +----+----+----+----+     +----+----+----+----+     
    | -- | -- | 00 | 01 |     | -- | xx |Used|Used|     
    +----+----+----+----+     +----+----+----+----+     
    | -- | -- | 02 | 03 |     | xx | xx | xx | xx |     
    +----+----+----+----+     +----+----+----+----+     
    | -- | -- | -- |04+j|     | -- | xx | -- |Used|      
    +----+----+----+----+     +----+----+----+----+     
    | -- | -- | -- | -- |     | -- | xx | -- | -- |     
    +----+----+----+----+     +----+----+----+----+     
    */

    // pt at the first and second points of conflocal
    double conf_pt_vals[18] = { 1.4142135623730951e+00, 0.0, 0.0,    // PT-Comp 00
      0.0, 0.0, 0.0, // PT-Comp 01
      0.0, 0.0, 0.0, // PT-Comp 02
      1.4142135623730951e+00, 0.0, 0.0, // PT-Comp 03
      0.0, 0.0, 0.0, // PT-Comp 04 + 0
      0.0, 0.0, 0.0  }; // PT-Comp 05

    for (int k = 0; k<num_pt_indices[vdim-1]; ++k) {
      for (int m = 0; m<NC; ++m) {
        int test_idx = k*NC + m;
        double expected = conf_pt_vals[test_idx];
        TEST_CHECK( gkyl_compare_double(conf_poisson_tensor_d[test_idx], expected, 1e-12) );
      }
    }
  }


  // Release the memory
  gkyl_array_release(conf_poisson_tensor);
}

void eval_vierbein_annulus_2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double q_r = xn[0];
  double q_theta = 0.0;
  fout[0] = 1.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = q_r;
}

void eval_vierbein_gradient_annulus_2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double q_r = xn[0];
  double q_theta = 0.0;
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = 1.0;

  fout[4] = 0.0;
  fout[5] = 0.0;
  fout[6] = 0.0;
  fout[7] = 0.0;
}

void
test_triad_1x2v_annulus_conf(int poly_order)
{

  double lower[] = {0.1, -1.0, -1.0}, upper[] = {1.0, 1.0, 1.0};
  int cells[] = {2, 2, 2};
  int vdim = 2, cdim = 1;
  int ndim = cdim+vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};
  double velLower[] = {lower[1], lower[2]}, velUpper[] = {upper[1], upper[2]};
  int velCells[] = {cells[1], cells[2]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

    struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  // velocity range
  int velGhost[] = { 0, 0 };
  struct gkyl_range velLocal, velLocal_ext; 
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);

  int confGhost[] = { 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = { confGhost[0], 0, 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Construct the input map
  struct gkyl_vlasov_triad_geom_inp inp_triad_geom = { 0 };
  inp_triad_geom.use_vierbein = true;
  inp_triad_geom.use_preset_geom = false;
  inp_triad_geom.eval_vierbein = eval_vierbein_annulus_2v; 
  inp_triad_geom.eval_vierbein_gradient = eval_vierbein_gradient_annulus_2v; 
  inp_triad_geom.eval_vierbein_ctx = 0; 
  inp_triad_geom.eval_vierbein_gradient_ctx = 0;  

  // Make the memory for arrays (modal)
  struct gkyl_array *conf_poisson_tensor; // Configuration space Poisson tensor representation

  // Size of the PT
  int num_pt_indices[3] = { 1 , 6, 18 }; 

  // Allocate arrays for covariant tangent basis 
  conf_poisson_tensor = mkarr(confBasis.num_basis*num_pt_indices[vdim-1], confLocal_ext.volume);

  // Construct the Geometry for this configuration
  gkyl_vlasov_triad_geom_new(&confGrid, &confLocal, confBasis, 
    &grid, &local, basis, inp_triad_geom, conf_poisson_tensor);

  // Iterate over the grid, conf space, checking output
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &confLocal);
  while (gkyl_range_iter_next(&iter)) {
    const double *conf_poisson_tensor_d = gkyl_array_cfetch(conf_poisson_tensor, gkyl_range_idx(&confLocal, iter.idx));
    
    int NC = confBasis.num_basis;

    // index size of the symmetric vierbein
    int n_sym = vdim*(vdim+1)/2; 

    // Double precision values of the expected output
    // pt at the first and second points of conflocal
    double conf_pt_pnt1_vals[18] = { 1.4142135623730951e+00, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      5.6865230654862211e+00, -3.3402132856134248e+00, 1.1941468087349127e+00, 0.0, 0.0, 0.0,
      5.6865230654862211e+00, -3.3402132856134248e+00, 1.1941468087349127e+00  };
    double conf_pt_pnt2_vals[18] = { 1.4142135623730951e+00, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      1.8807796203407585e+00, -3.3402132856134242e-01, 5.0077124237270627e-02, 0.0, 0.0, 0.0,
      1.8807796203407585e+00, -3.3402132856134242e-01, 5.0077124237270627e-02  };

    for (int k = 0; k<num_pt_indices[vdim-1]; ++k) {
      for (int m = 0; m<NC; ++m) {
        int test_idx = k*NC + m;
        double expected = 0;
        if (iter.idx[0] == 1) expected = conf_pt_pnt1_vals[test_idx];
        if (iter.idx[0] == 2) expected = conf_pt_pnt2_vals[test_idx];
        //printf("conf_poisson_tensor_d[%d]: %1.16e\n",test_idx,conf_poisson_tensor_d[test_idx]);
        TEST_CHECK( gkyl_compare_double(conf_poisson_tensor_d[test_idx], expected, 1e-12) );
      }
    }
  }

  // Release the memory
  gkyl_array_release(conf_poisson_tensor);
}

void
test_triad_math_2v()
{
  // Test: Can the correct bracket be constructed for:
  // 2D Annular Disk Coordiantes
  int vdim = 2;

  // Accuracy of the tests
  double eps = 1e-12;
  double vierbein[4];

  // Randomly choosen coordinates
  double r = 3.2;
  double theta = 1.4;

  // Assign covaraint tangent basis 
  vierbein[0] = 1.0;
  vierbein[1] = 0.0;
  vierbein[2] = 0.0;
  vierbein[3] = r;

  // 1. Test kernel_vierbein_inv_2v()
  double vierbein_inv[4];
  kernel_vierbein_inv_2v(vierbein, vierbein_inv);

  // vierbein_inv only has symmetric entries as well
  TEST_CHECK(gkyl_compare_double(vierbein_inv[0], 1.0, eps));
  TEST_CHECK(gkyl_compare_double(vierbein_inv[1], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(vierbein_inv[2], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(vierbein_inv[3], 1.0/r, eps));

  // 2. Test kernel_conf_poisson_tensor_2v()
  double conf_poisson_tensor[6];

  // Assign the triad basis gradient
  double vierbein_gradient[8];

  // gradient in r
  vierbein_gradient[0] = 0.0;
  vierbein_gradient[1] = 0.0;

  vierbein_gradient[2] = 0.0;
  vierbein_gradient[3] = 1.0;

  // gradient in theta
  vierbein_gradient[4] = 0.0;
  vierbein_gradient[5] = 0.0;

  vierbein_gradient[6] = 0.0;
  vierbein_gradient[7] = 0.0;

  kernel_conf_poisson_tensor_vierbein_2v(vierbein_inv, vierbein_gradient, conf_poisson_tensor);

  // Pi_{xx} block
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[0], 1.0, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[1], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[2], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[3], 1/r, eps));

  // Pi_{pp} block, p_r coeffs 
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[4], 0.0, eps));

  // Pi_{pp} block, p_theta coeffs 
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[5], 1.0/r, eps));
}

void
test_triad_math_3v()
{

  // Test: Can the correct bracket be constructed for:
  // 3D Spherical Coordiantes
  int vdim = 3;

  // Accuracy of the tests
  double eps = 1e-12;
  double vierbein[9];

  // Randomly choosen coordinates
  double r = 3.2;
  double theta = 1.4;
  double phi = 0.8;

  // Assign covaraint tangent basis 
  vierbein[0] = 1.0;
  vierbein[1] = 0.0;
  vierbein[2] = 0.0;

  vierbein[3] = 0.0;
  vierbein[4] = r;
  vierbein[5] = 0.0;

  vierbein[6] = 0.0;
  vierbein[7] = 0.0;
  vierbein[8] = r * sin(theta);

  // 1. Test kernel_vierbein_inv_3v()
  double vierbein_inv[9];
  kernel_vierbein_inv_3v(vierbein, vierbein_inv);

  // vierbein_inv only has symmetric entries as well
  TEST_CHECK(gkyl_compare_double(vierbein_inv[0], 1.0, eps));
  TEST_CHECK(gkyl_compare_double(vierbein_inv[1], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(vierbein_inv[2], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(vierbein_inv[3], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(vierbein_inv[4], 1.0/r, eps));
  TEST_CHECK(gkyl_compare_double(vierbein_inv[5], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(vierbein_inv[6], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(vierbein_inv[7], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(vierbein_inv[8], 1.0/(r*sin(theta)), eps));

  // 5. Test kernel_conf_poisson_tensor_3v()
  double conf_poisson_tensor[18];

  // Assign the triad basis gradient
  double vierbein_gradient[27];

  // gradient in r
  vierbein_gradient[0] = 0.0;
  vierbein_gradient[1] = 0.0;
  vierbein_gradient[2] = 0.0;

  vierbein_gradient[3] = 0.0;
  vierbein_gradient[4] = 1.0;
  vierbein_gradient[5] = 0.0;

  vierbein_gradient[6] = 0.0;
  vierbein_gradient[7] = 0.0;
  vierbein_gradient[8] = sin(theta);

  // gradient in theta
  vierbein_gradient[9] = 0.0;
  vierbein_gradient[10] = 0.0;
  vierbein_gradient[11] = 0.0;

  vierbein_gradient[12] = 0.0;
  vierbein_gradient[13] = 0.0;
  vierbein_gradient[14] = 0.0;

  vierbein_gradient[15] = 0.0;
  vierbein_gradient[16] = 0.0;
  vierbein_gradient[17] = r*cos(theta);

  // gradient in phi
  vierbein_gradient[18] = 0.0;
  vierbein_gradient[19] = 0.0;
  vierbein_gradient[20] = 0.0;

  vierbein_gradient[21] = 0.0;
  vierbein_gradient[22] = 0.0;
  vierbein_gradient[23] = 0.0;

  vierbein_gradient[24] = 0.0;
  vierbein_gradient[25] = 0.0;
  vierbein_gradient[26] = 0.0;

  kernel_conf_poisson_tensor_vierbein_3v(vierbein_inv, vierbein_gradient, conf_poisson_tensor);

  // Pi_{xx} block
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[0], 1.0, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[1], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[2], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[3], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[4], 1.0/r, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[5], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[6], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[7], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[8], 1.0/(r*sin(theta)), eps));

  // Pi_{pp} block, p_r coeffs 
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[9], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[10], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[11], 0.0, eps));

  // Pi_{pp} block, p_theta coeffs 
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[12], 1.0/r, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[13], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[14], 0.0, eps));

  // Pi_{pp} block, p_phi coeffs 
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[15], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[16], 1.0/r, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[17], cos(theta)/(r*sin(theta)), eps));

}

void
eval_ks_rphi_vierbein_2v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Parameters
  double a = 0.2;
  double M = 0.7;

  // Coordinates
  double r = xn[0];
  double phi = xn[1];

  // Intermediate Variables
  double rho_sq = r * r ;

  // Metric spatial covariant components
  double h_rr = ( 1.0 + 2.0 * M * r / rho_sq);
  double h_pp = ( rho_sq + a * a * ( 1.0 + 2.0 * M * r / rho_sq )   );
  double h_thth = rho_sq;
  double h_rp = - a * ( 1.0 + 2.0 * M * r / rho_sq );

  // Vierbein: e_i^a = g_i . sigma^a
  fout[0] = sqrt( h_rr );
  fout[1] = 0.0;
  fout[2] = h_rp / sqrt( h_rr );
  fout[3] = sqrt( h_pp - h_rp * h_rp / h_rr );
}

void
eval_ks_rphi_vierbein_gradient_2v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  double a = 0.2;
  double M = 0.7;

  double r = xn[0];
  double phi = xn[1];

  // Intermediate Variables
  double rho_sq = r * r ;

  // Gradient w.r.t. r: d(e_i^a)/dr
  fout[0] = (M * ( - r * r)) / (pow(rho_sq, 1.5) * sqrt(rho_sq + 2.0 * M * r));
  fout[1] = 0.0;
  fout[2] = -(M * a * ( - r * r)) / (pow(rho_sq, 1.5) * sqrt(rho_sq + 2.0 * M * r));
  fout[3] = r / sqrt(rho_sq);

  // Gradient w.r.t. phi: d(e_i^a)/dphi
  fout[4] = 0.0;
  fout[5] = 0.0;
  fout[6] = 0.0;
  fout[7] = 0.0;

}

void
test_triad_2x2v_rphi_ks_pnt(int poly_order)
{

  double lower[] = {1.2, 1.0, -1.0, -1.0}, upper[] = {1.201, 1.001, 1.0, 1.0};
  int cells[] = {2, 2, 2, 2};
  int vdim = 2, cdim = 2;
  int ndim = cdim+vdim;

  double confLower[] = {lower[0], lower[1]}, confUpper[] = {upper[0], upper[1]};
  int confCells[] = {cells[0], cells[1]};
  double velLower[] = { lower[2], lower[3]}, velUpper[] = { upper[2], upper[3]};
  int velCells[] = {cells[2], cells[3]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

    struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  // velocity range
  int velGhost[] = { 0, 0 };
  struct gkyl_range velLocal, velLocal_ext; 
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);

  int confGhost[] = { 1, 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = { confGhost[0], confGhost[1], 0, 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Construct the input map
  struct gkyl_vlasov_triad_geom_inp inp_triad_geom = { 0 };
  inp_triad_geom.use_vierbein = true;
  inp_triad_geom.use_preset_geom = false;
  inp_triad_geom.eval_vierbein = eval_ks_rphi_vierbein_2v; 
  inp_triad_geom.eval_vierbein_gradient = eval_ks_rphi_vierbein_gradient_2v; 
  inp_triad_geom.eval_vierbein_ctx = 0; 
  inp_triad_geom.eval_vierbein_gradient_ctx = 0;  

  // Make the memory for arrays (modal)
  struct gkyl_array *conf_poisson_tensor; // Configuration space Poisson tensor representation

  // Size of the PT
  int num_pt_indices[3] = { 1 , 6, 18 }; 

  // Allocate arrays for covariant tangent basis 
  conf_poisson_tensor = mkarr(confBasis.num_basis*num_pt_indices[vdim-1], confLocal_ext.volume);

  // Construct the Geometry for this configuration
  gkyl_vlasov_triad_geom_new(&confGrid, &confLocal, confBasis, 
    &grid, &local, basis, inp_triad_geom, conf_poisson_tensor);

  // Iterate over the grid, conf space, checking output
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &confLocal);
  while (gkyl_range_iter_next(&iter)) {
    const double *conf_poisson_tensor_d = gkyl_array_cfetch(conf_poisson_tensor, gkyl_range_idx(&confLocal, iter.idx));
    
    int NC = confBasis.num_basis;

    // index size of the symmetric vierbein
    int n_sym = vdim*(vdim+1)/2; 

    // Double precision values of the expected output
    // pt at the first and second points of conflocal
    double conf_pt_pnt1_vals[48] = {
      1.3588086393098029e+00,  4.3989488197243091e-05,  0.0000000000000000e+00,  0.0000000000000000e+00,
     -2.8209032673875084e-09,  1.5884597076836006e-16,  0.0000000000000000e+00, -4.4892812507152318e-21,

      3.3326390817298623e-01, -4.0077069914578266e-05,  0.0000000000000000e+00,  4.6259292692714869e-18,
      4.3107081734112447e-09,  4.0163355468972541e-17,  0.0000000000000000e+00,  3.3710617091548934e-21,

      0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
      0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,

      1.6663195408649309e+00, -2.0038534957297145e-04,  0.0000000000000000e+00,  0.0000000000000000e+00,
      2.1553540958340271e-08,  1.8375096125123659e-16,  0.0000000000000000e+00,  1.7028157575878893e-20,

      7.4747990717562870e-02, -2.2127022055885313e-05,  0.0000000000000000e+00,  0.0000000000000000e+00,
      3.9991364151008178e-09,  7.6698926571441654e-18,  0.0000000000000000e+00,  1.9701838207745808e-21,

      1.1321046895815821e+00, -9.9492398515537507e-05,  0.0000000000000000e+00, -1.8503717077085947e-17,
      8.3511852151311742e-09,  1.3368141403477171e-16,  0.0000000000000000e+00,  9.9229411160645424e-21
    };

    double conf_pt_pnt2_vals[48] = {
      1.3588086393098029e+00,  4.3989488197243091e-05,  0.0000000000000000e+00,  0.0000000000000000e+00,
     -2.8209032673875084e-09,  1.5884597076836006e-16,  0.0000000000000000e+00, -4.4892812507152318e-21,

      3.3326390817298623e-01, -4.0077069914578266e-05,  0.0000000000000000e+00,  4.6259292692714869e-18,
      4.3107081734112447e-09,  4.0163355468972541e-17,  0.0000000000000000e+00,  3.3710617091548934e-21,

      0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
      0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,

      1.6663195408649309e+00, -2.0038534957297145e-04,  0.0000000000000000e+00,  0.0000000000000000e+00,
      2.1553540958340271e-08,  1.8375096125123659e-16,  0.0000000000000000e+00,  1.7028157575878893e-20,

      7.4747990717562870e-02, -2.2127022055885313e-05,  0.0000000000000000e+00,  0.0000000000000000e+00,
      3.9991364151008178e-09,  7.6698926571441654e-18,  0.0000000000000000e+00,  1.9701838207745808e-21,

      1.1321046895815821e+00, -9.9492398515537507e-05,  0.0000000000000000e+00, -1.8503717077085947e-17,
      8.3511852151311742e-09,  1.3368141403477171e-16,  0.0000000000000000e+00,  9.9229411160645424e-21
    };

    double conf_pt_pnt3_vals[48] = {
      1.3589609855304008e+00,  4.3967646125575739e-05,  0.0000000000000000e+00,  0.0000000000000000e+00,
     -2.8186954269208275e-09,  1.2945958434809000e-16,  0.0000000000000000e+00, -3.5482011072808203e-21,

      3.3312513494050400e-01, -4.0043700164429930e-05,  0.0000000000000000e+00,  4.6259292692714869e-18,
      4.3053254227575864e-09,  3.0521154308908952e-17,  0.0000000000000000e+00,  4.6283179257499702e-21,

      0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
      0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,

      1.6656256747025202e+00, -2.0021850082239002e-04,  0.0000000000000000e+00, -1.8503717077085947e-17,
      2.1526627144731539e-08,  2.0492889448999093e-16,  0.0000000000000000e+00,  2.1966070824554527e-20,

      7.4671394086855941e-02, -2.2096072509741311e-05,  0.0000000000000000e+00,  0.0000000000000000e+00,
      3.9920039782016186e-09,  8.0867145517377461e-18,  0.0000000000000000e+00,  2.4846788383995600e-21,

      1.1317601498076699e+00, -9.9427743788891502e-05,  0.0000000000000000e+00,  1.8503717077085947e-17,
      8.3425950860921977e-09,  9.8703706655908211e-17,  0.0000000000000000e+00,  9.1561817091435301e-21
    };

    double conf_pt_pnt4_vals[48] = {
      1.3589609855304008e+00,  4.3967646125575739e-05,  0.0000000000000000e+00,  0.0000000000000000e+00,
     -2.8186954269208275e-09,  1.2945958434809000e-16,  0.0000000000000000e+00, -3.5482011072808203e-21,

      3.3312513494050400e-01, -4.0043700164429930e-05,  0.0000000000000000e+00,  4.6259292692714869e-18,
      4.3053254227575864e-09,  3.0521154308908952e-17,  0.0000000000000000e+00,  4.6283179257499702e-21,

      0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
      0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,

      1.6656256747025202e+00, -2.0021850082239002e-04,  0.0000000000000000e+00, -1.8503717077085947e-17,
      2.1526627144731539e-08,  2.0492889448999093e-16,  0.0000000000000000e+00,  2.1966070824554527e-20,

      7.4671394086855941e-02, -2.2096072509741311e-05,  0.0000000000000000e+00,  0.0000000000000000e+00,
      3.9920039782016186e-09,  8.0867145517377461e-18,  0.0000000000000000e+00,  2.4846788383995600e-21,

      1.1317601498076699e+00, -9.9427743788891502e-05,  0.0000000000000000e+00,  1.8503717077085947e-17,
      8.3425950860921977e-09,  9.8703706655908211e-17,  0.0000000000000000e+00,  9.1561817091435301e-21
    };


    for (int k = 0; k<num_pt_indices[vdim-1]; ++k) {
      for (int m = 0; m<NC; ++m) {
        int test_idx = k*NC + m;
        double expected = 0;
        if (iter.idx[0] == 1 && iter.idx[1] == 1) expected = conf_pt_pnt1_vals[test_idx];
        if (iter.idx[0] == 1 && iter.idx[1] == 2) expected = conf_pt_pnt2_vals[test_idx];
        if (iter.idx[0] == 2 && iter.idx[1] == 1) expected = conf_pt_pnt3_vals[test_idx];
        if (iter.idx[0] == 2 && iter.idx[1] == 2) expected = conf_pt_pnt4_vals[test_idx];
        //printf("conf_poisson_tensor_d[%d]: %1.16e\n",test_idx,conf_poisson_tensor_d[test_idx]);
        // GJ (6/9/26) The test is low because this appears okay on GPU but fails even at 1e-10
        TEST_CHECK( gkyl_compare_double(conf_poisson_tensor_d[test_idx], expected, 1e-8) );
      }
    }
  }

  // Release the memory
  gkyl_array_release(conf_poisson_tensor);
}

void eval_vierbein_obl_sph_3v(double t, const double *xn, double* restrict fout, void *ctx)
{
  // Parameters
  double a = 0.2;
  double M = 0.7;

  // Coordinates
  double r = xn[0];
  double theta = 1.0;
  double phi = 1.0;

  // Intermediate Variables
  double rho_sq = r * r + a * a * cos(theta) * cos(theta);

  // Metric spatial covariant components
  double h_rr = ( 1.0 + 2.0 * M * r / rho_sq);
  double h_pp = sin(theta) * sin(theta) * ( rho_sq + a * a * ( 1.0 + 2.0 * M * r / rho_sq ) * sin(theta) * sin(theta)  );
  double h_thth = rho_sq;
  double h_rp = - a * ( 1.0 + 2.0 * M * r / rho_sq ) * sin(theta) * sin(theta);

  // Vierbein: e_i^a = g_i . sigma^a
  fout[0] = sqrt( h_rr );
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = 0.0;
  fout[4] = sqrt( rho_sq );
  fout[5] = 0.0;
  fout[6] = h_rp / sqrt( h_rr );
  fout[7] = 0.0;
  fout[8] = sqrt( h_pp - h_rp * h_rp / h_rr );
}

void eval_vierbein_gradient_obl_sph_3v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double a = 0.2;
  double M = 0.7;

  double r = xn[0];
  double theta = 1.0;
  double phi = 1.0;

  // Intermediate Variables
  double rho_sq = r * r + a * a * cos(theta) * cos(theta);

  // Gradient w.r.t. r: d(e_i^a)/dr
  fout[0] = (M * (a * a * cos(theta) * cos(theta) - r * r)) / (pow(rho_sq, 1.5) * sqrt(rho_sq + 2.0 * M * r));
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = 0.0;
  fout[4] = r / sqrt(rho_sq);
  fout[5] = 0.0;
  fout[6] = -(M * a * sin(theta) * sin(theta) * (a * a * cos(theta) * cos(theta) - r * r)) / (pow(rho_sq, 1.5) * sqrt(rho_sq + 2.0 * M * r));
  fout[7] = 0.0;
  fout[8] = (r * sin(theta)) / sqrt(rho_sq);

  // Gradient w.r.t. theta: d(e_i^a)/dtheta
  fout[9]  = (M * a * a * r * sin(2.0 * theta)) / (pow(rho_sq, 1.5) * sqrt(rho_sq + 2.0 * M * r));
  fout[10] = 0.0;
  fout[11] = 0.0;
  fout[12] = 0.0;
  fout[13] = -(a * a * sin(2.0 * theta)) / (2.0 * sqrt(rho_sq));
  fout[14] = 0.0;
  fout[15] = -(2.0 * a * cos(theta) * sin(theta) * (a * a * a * a * cos(theta) * cos(theta) * cos(theta) * cos(theta) + 2.0 * M * r * r * r + r * r * r * r + 2.0 * a * a * r * r * cos(theta) * cos(theta) + 2.0 * M * a * a * r * cos(theta) * cos(theta) + M * a * a * r * sin(theta) * sin(theta))) / (pow(rho_sq, 1.5) * sqrt(rho_sq + 2.0 * M * r));
  fout[16] = 0.0;
  fout[17] = (sin(4.0 * theta) * a * a + 2.0 * sin(2.0 * theta) * r * r) / (2.0 * sqrt(1.0 - cos(2.0 * theta)) * sqrt(a * a * cos(2.0 * theta) + a * a + 2.0 * r * r));

  // Gradient w.r.t. theta: d(e_i^a)/dphi
  fout[18] = 0.0;
  fout[19] = 0.0;
  fout[20] = 0.0;
  fout[21] = 0.0;
  fout[22] = 0.0;
  fout[23] = 0.0;
  fout[24] = 0.0;
  fout[25] = 0.0;
  fout[26] = 0.0;
}

void
test_triad_3v_obl_sph_conf(int poly_order)
{

  double lower[] = {0.5, -1.0, -1.0, -1.0}, upper[] = {1.5, 1.0, 1.0, 1.0};
  int cells[] = {2, 2, 2, 2};
  int vdim = 3, cdim = 1;
  int ndim = cdim+vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};
  double velLower[] = {lower[1], lower[2], lower[3]}, velUpper[] = {upper[1], upper[2], upper[3]};
  int velCells[] = {cells[1], cells[2], cells[3]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

    struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  // velocity range
  int velGhost[] = { 0, 0, 0 };
  struct gkyl_range velLocal, velLocal_ext; 
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);

  int confGhost[] = { 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = { confGhost[0], 0, 0, 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Construct the input map
  struct gkyl_vlasov_triad_geom_inp inp_triad_geom = { 0 };
  inp_triad_geom.use_vierbein = true;
  inp_triad_geom.use_preset_geom = false;
  inp_triad_geom.eval_vierbein = eval_vierbein_obl_sph_3v; 
  inp_triad_geom.eval_vierbein_gradient = eval_vierbein_gradient_obl_sph_3v; 
  inp_triad_geom.eval_vierbein_ctx = 0; 
  inp_triad_geom.eval_vierbein_gradient_ctx = 0;  

  // Make the memory for arrays (modal)
  struct gkyl_array *conf_poisson_tensor; // Configuration space Poisson tensor representation

  // Size of the PT
  int num_pt_indices[3] = { 1 , 6, 18 }; 

  // Allocate arrays for covariant tangent basis 
  conf_poisson_tensor = mkarr(confBasis.num_basis*num_pt_indices[vdim-1], confLocal_ext.volume);

  // Construct the Geometry for this configuration
  gkyl_vlasov_triad_geom_new(&confGrid, &confLocal, confBasis, 
    &grid, &local, basis, inp_triad_geom, conf_poisson_tensor);

  // Iterate over the grid, conf space, checking output
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &confLocal);
  while (gkyl_range_iter_next(&iter)) {
    const double *conf_poisson_tensor_d = gkyl_array_cfetch(conf_poisson_tensor, gkyl_range_idx(&confLocal, iter.idx));
    
    int NC = confBasis.num_basis;

    // index size of the symmetric vierbein
    int n_sym = vdim*(vdim+1)/2; 

    // Double precision values of the expected output
    // pt at the first and second points of conflocal
    double conf_pt_pnt1_vals[54] = {
      8.3618107409643117e-01,  5.1457473130182417e-02, -4.1806945016435162e-03,
     -0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
      3.2637882582504923e-01, -6.6002566736818952e-02,  1.0986860997434837e-02,
      0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
      1.9393349962004920e+00, -3.9218563639414727e-01,  6.5283659186080201e-02,
     -0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
      0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
      0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
      2.3046962179488957e+00, -4.6607149084728239e-01,  7.7582781064024775e-02,
     -5.1231097570368596e-02,  3.6136517546674948e-02, -1.1641616453888226e-02,
      1.4916600138443334e-01, -6.9669818422142496e-02,  1.7263355736463301e-02,
     -1.0477839575729653e+00,  5.1172555580643986e-01, -1.3640211096772709e-01,
      1.1180531478575433e+00, -1.4554580599782768e-01,  1.7632304440892480e-02,
      0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
     -4.6066819422732763e-01,  1.8911477290676781e-01, -4.5071700813560885e-02,
      0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
      1.1180531478575431e+00, -1.4554580599782768e-01,  1.7632304440892591e-02,
      1.1695894129651985e+00, -2.0365162950883728e-01,  2.7132936474682574e-02
    };

    double conf_pt_pnt2_vals[54] = {
      9.7118391687045558e-01,  2.9561911060611016e-02, -1.7984560831458061e-03,
     -0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
      1.9227801821038229e-01, -2.2622673560699226e-02,  2.3097493850544964e-03,
      0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
      1.1425112789496739e+00, -1.3442337269575541e-01,  1.3724474323863870e-02,
     -0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
      0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
      0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
      1.3577548122000940e+00, -1.5974807821382403e-01,  1.6310098115857333e-02,
     -7.6223193723774160e-03,  3.1845707381018599e-03, -6.8607093083360800e-04,
      4.1607608429150726e-02, -1.2112759983408558e-02,  2.0164185412722604e-03,
     -2.9711762483787008e-01,  8.0478190612048625e-02, -1.3391251162867419e-02,
      7.7865238138877158e-01, -6.6390243442267571e-02,  5.2519835403373391e-03,
      0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
     -1.5689662645424968e-01,  3.7219985505675783e-02, -5.6043258481635787e-03,
      0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
      7.7865238138877146e-01, -6.6390243442267585e-02,  5.2519835403373946e-03,
      7.1944415265280803e-01, -8.1198844969437509e-02,  7.7997880735233126e-03
    };


    for (int k = 0; k<num_pt_indices[vdim-1]; ++k) {
      for (int m = 0; m<NC; ++m) {
        int test_idx = k*NC + m;
        double expected = 0;
        if (iter.idx[0] == 1) expected = conf_pt_pnt1_vals[test_idx];
        if (iter.idx[0] == 2) expected = conf_pt_pnt2_vals[test_idx];
        //printf("conf_poisson_tensor_d[%d]: %1.16e\n",test_idx,conf_poisson_tensor_d[test_idx]);
        TEST_CHECK( gkyl_compare_double(conf_poisson_tensor_d[test_idx], expected, 1e-12) );
      }
    }
  }

  // Release the memory
  gkyl_array_release(conf_poisson_tensor);
}


void
test_triad_3v_obl_sph_conf_pnt(int poly_order)
{

  // The small radial domain is to make the PT computable point-like
  double lower[] = {1.2, -1.0, -1.0, -1.0}, upper[] = {1.201, 1.0, 1.0, 1.0};
  int cells[] = {2, 2, 2, 2};
  int vdim = 3, cdim = 1;
  int ndim = cdim+vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};
  double velLower[] = {lower[1], lower[2], lower[3]}, velUpper[] = {upper[1], upper[2], upper[3]};
  int velCells[] = {cells[1], cells[2], cells[3]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

    struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  // velocity range
  int velGhost[] = { 0, 0, 0 };
  struct gkyl_range velLocal, velLocal_ext; 
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);

  int confGhost[] = { 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = { confGhost[0], 0, 0, 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Construct the input map
  struct gkyl_vlasov_triad_geom_inp inp_triad_geom = { 0 };
  inp_triad_geom.use_vierbein = true;
  inp_triad_geom.use_preset_geom = false;
  inp_triad_geom.eval_vierbein = eval_vierbein_obl_sph_3v; 
  inp_triad_geom.eval_vierbein_gradient = eval_vierbein_gradient_obl_sph_3v; 
  inp_triad_geom.eval_vierbein_ctx = 0; 
  inp_triad_geom.eval_vierbein_gradient_ctx = 0;  

  // Make the memory for arrays (modal)
  struct gkyl_array *conf_poisson_tensor; // Configuration space Poisson tensor representation

  // Size of the PT
  int num_pt_indices[3] = { 1 , 6, 18 }; 

  // Allocate arrays for covariant tangent basis 
  conf_poisson_tensor = mkarr(confBasis.num_basis*num_pt_indices[vdim-1], confLocal_ext.volume);

  // Construct the Geometry for this configuration
  gkyl_vlasov_triad_geom_new(&confGrid, &confLocal, confBasis, 
    &grid, &local, basis, inp_triad_geom, conf_poisson_tensor);

  // Iterate over the grid, conf space, checking output
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &confLocal);
  while (gkyl_range_iter_next(&iter)) {
    const double *conf_poisson_tensor_d = gkyl_array_cfetch(conf_poisson_tensor, gkyl_range_idx(&confLocal, iter.idx));
    
    int NC = confBasis.num_basis;

    // index size of the symmetric vierbein
    int n_sym = vdim*(vdim+1)/2; 

    // Double precision values of the expected output
    // pt at the first and second points of conflocal
    double conf_pt_pnt1_vals[54] = {
      9.6290932667188200e-01,  3.0557192303253710e-05, -1.9060309741679256e-09,
     -0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
      1.9749649683737491e-01, -2.3559225423118542e-05,  2.5034786937286079e-09,
      0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
      1.1735193512493030e+00, -1.3998834094392952e-04,  1.4875609027276800e-08,
     -0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
      0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
      0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
      1.3946046535605872e+00, -1.6636145924360524e-04,  1.7678100999715272e-08,
     -7.8825693985710692e-03,  3.2533094023435649e-06, -7.5997620951496514e-10,
      4.3421825783067178e-02, -1.2614050461930800e-05,  2.2226953151201167e-09,
     -3.0909497157018945e-01,  8.3552360881739072e-05, -1.4753654335173865e-08,
      7.9580674413801245e-01, -6.8907451802267957e-05,  5.6345035264203492e-09,
      0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
     -1.6322310167933296e-01,  3.8783726313331308e-05, -6.1565631082591592e-09,
      0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
      7.9580674413801233e-01, -6.8907451802086651e-05,  5.6345035819315004e-09,
      7.3881286793473677e-01, -8.4626522972194968e-05,  8.4315775465537861e-09
    };

    double conf_pt_pnt2_vals[54] = {
      9.6301515432534812e-01,  3.0542433744200089e-05, -1.9046126736654344e-09,
     -0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
      1.9741491886046661e-01, -2.3539845475652491e-05,  2.5004028958501497e-09,
      0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
      1.1730346169068171e+00, -1.3987318577010459e-04,  1.4857332885849471e-08,
     -0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
      0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
      0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
      1.3940285976403748e+00, -1.6622460939880196e-04,  1.7656381705803026e-08,
     -7.8713097929162569e-03,  3.2474291031404502e-06, -7.5831105403301365e-10,
      4.3378159233906165e-02, -1.2596848405148490e-05,  2.2188576674860956e-09,
     -3.0880573552501606e-01,  8.3438181690105028e-05, -1.4727293548494199e-08,
      7.9556811729279164e-01, -6.8863828717794862e-05,  5.6289295621407699e-09,
      0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
     -1.6308883346440881e-01,  3.8736076679782285e-05, -6.1465219895234848e-09,
      0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
      7.9556811729279164e-01, -6.8863828717749529e-05,  5.6289296089518828e-09,
      7.3851982613910860e-01, -8.4561248795774556e-05,  8.4221442670652718e-09
    };


    for (int k = 0; k<num_pt_indices[vdim-1]; ++k) {
      for (int m = 0; m<NC; ++m) {
        int test_idx = k*NC + m;
        double expected = 0;
        if (iter.idx[0] == 1) expected = conf_pt_pnt1_vals[test_idx];
        if (iter.idx[0] == 2) expected = conf_pt_pnt2_vals[test_idx];
        //printf("conf_poisson_tensor_d[%d]: %1.16e\n",test_idx,conf_poisson_tensor_d[test_idx]);
        //printf("conf_poisson_tensor_d[%d]: %1.16e (expected: %1.16e)\n",test_idx,conf_poisson_tensor_d[test_idx], expected);
        TEST_CHECK( gkyl_compare_double(conf_poisson_tensor_d[test_idx], expected, 1e-8) );
      }
    }
  }

  // Release the memory
  gkyl_array_release(conf_poisson_tensor);
}

void test_triad_1v() { test_triad_math_1v(); }
void test_triad_1x1v_flat() { test_triad_1x1v_flat_conf(2); }
void test_triad_1x2v_flat() { test_triad_1x2v_flat_conf(2); }
void test_triad_1x2v_annulus() { test_triad_1x2v_annulus_conf(2); }
void test_triad_2v() { test_triad_math_2v(); }
void test_triad_3v() { test_triad_math_3v(); }
void test_triad_2x2v_rphi_ks() { test_triad_2x2v_rphi_ks_pnt(2); }
void test_triad_3v_obl_sph() { test_triad_3v_obl_sph_conf(2); }
void test_triad_3v_obl_sph_pnt() { test_triad_3v_obl_sph_conf_pnt(2); }

TEST_LIST = {
  { "test_triad_1v", test_triad_1v}, 
  { "test_triad_1x1v_flat", test_triad_1x1v_flat}, 
  { "test_triad_1x2v_flat", test_triad_1x2v_flat}, 
  { "test_triad_1x2v_annulus", test_triad_1x2v_annulus}, 
  { "test_triad_2v", test_triad_2v},
  { "test_triad_3v", test_triad_3v},
  { "test_triad_2x2v_rphi_ks", test_triad_2x2v_rphi_ks},
  { "test_triad_3v_obl_sph", test_triad_3v_obl_sph},
  { "test_triad_3v_obl_sph_pnt", test_triad_3v_obl_sph_pnt},
  {NULL, NULL}
};
