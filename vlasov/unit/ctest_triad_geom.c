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

  // 1. Test kernel_metric_1v()
  double cov_tangent_basis[1];
  double h_ij[1];

  // Randomly choosen coordinates
  double z1 = 3.2;

  // Assign covaraint tangent basis 
  cov_tangent_basis[0] = 2*z1;

  kernel_metric_1v(cov_tangent_basis, h_ij);

  // h_ij only has symmetric entries
  TEST_CHECK(gkyl_compare_double(h_ij[0], 4.0 * z1 * z1, eps));

  // 2. Test kernel_metric_inv_1v()
  double h_ij_inv[1];
  kernel_metric_inv_1v(h_ij, h_ij_inv);

  // h_ij_inv only has symmetric entries as well
  TEST_CHECK(gkyl_compare_double(h_ij_inv[0], 1.0 / (4.0 * z1 * z1), eps));

  // 3. Test kernel_metric_det_1v()
  double det_h[1];
  kernel_metric_det_1v(h_ij, det_h);

  TEST_CHECK(gkyl_compare_double(det_h[0], 2.0 * z1, eps));

  // 4. Test compute_nu_inv_1v()

  double triad_basis[1];
  double nu_inv[1];

  // Assign the triad basis 
  triad_basis[0] = 1;

  compute_nu_inv_1v(h_ij_inv, triad_basis, cov_tangent_basis, nu_inv);

  TEST_CHECK(gkyl_compare_double(nu_inv[0], 1.0 / (2.0 * z1), eps));

  // 5. Test kernel_conf_poisson_tensor_1v()
  double conf_poisson_tensor[1];

  // Assign the triad basis gradient
  double triad_basis_gradient[1];

  // gradient in r
  triad_basis_gradient[0] = 0.0;

  kernel_conf_poisson_tensor_1v(h_ij_inv, triad_basis, cov_tangent_basis,
    triad_basis_gradient, conf_poisson_tensor);

  // Pi_{xx} block
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[0], 1.0 / (2*z1), eps));

}

void eval_triad_basis_1v(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = 1;
}

void eval_cov_tangent_basis_1v(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = 1;
}

void eval_triad_basis_gradient_1v(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = 0;
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
  struct gkyl_vlasov_triad_geom_inp inp_basis_vectors = { 0 };
  inp_basis_vectors.use_vierbein = false;
  inp_basis_vectors.use_preset_geom = false;
  inp_basis_vectors.eval_cov_tangent_basis = eval_cov_tangent_basis_1v; 
  inp_basis_vectors.eval_triad_basis = eval_triad_basis_1v; 
  inp_basis_vectors.eval_triad_basis_gradient = eval_triad_basis_gradient_1v; 
  inp_basis_vectors.eval_cov_tangent_basis_ctx = 0; 
  inp_basis_vectors.eval_triad_basis_ctx = 0; 
  inp_basis_vectors.eval_triad_basis_gradient_ctx = 0; 

  // Make the memory for arrays (modal)
  struct gkyl_array *cov_tangent_basis; // Covariant tangent basis
  struct gkyl_array *triad_basis; // Triad basis
  struct gkyl_array *conf_poisson_tensor; // Configuration space Poisson tensor representation

  // Size of the PT
  int num_pt_indices[3] = { 1 , 6, 18 }; 

  // Allocate arrays for covariant tangent basis 
  cov_tangent_basis = mkarr(basis.num_basis*vdim*vdim, local_ext.volume);
  triad_basis = mkarr(basis.num_basis*vdim*vdim, local_ext.volume);
  conf_poisson_tensor = mkarr(basis.num_basis*num_pt_indices[vdim-1], local_ext.volume);

  // Construct the Geometry for this configuration
  gkyl_vlasov_triad_geom_new(&confGrid, &confLocal, confBasis, 
    &grid, &local, basis, inp_basis_vectors, conf_poisson_tensor);

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
  gkyl_array_release(cov_tangent_basis);
  gkyl_array_release(triad_basis);
  gkyl_array_release(conf_poisson_tensor);
}




void eval_triad_basis_flat_2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = 1.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = 1.0;
}

void eval_cov_tangent_basis_flat_2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = 1.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = 1.0;
}

void eval_triad_basis_gradient_flat_2v(double t, const double *xn, double* restrict fout, void *ctx)
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
  struct gkyl_vlasov_triad_geom_inp inp_basis_vectors = { 0 };
  inp_basis_vectors.use_vierbein = false;
  inp_basis_vectors.use_preset_geom = false;
  inp_basis_vectors.eval_cov_tangent_basis = eval_cov_tangent_basis_flat_2v; 
  inp_basis_vectors.eval_triad_basis = eval_triad_basis_flat_2v; 
  inp_basis_vectors.eval_triad_basis_gradient = eval_triad_basis_gradient_flat_2v; 
  inp_basis_vectors.eval_cov_tangent_basis_ctx = 0; 
  inp_basis_vectors.eval_triad_basis_ctx = 0; 
  inp_basis_vectors.eval_triad_basis_gradient_ctx = 0; 

  // Make the memory for arrays (modal)
  struct gkyl_array *cov_tangent_basis; // Covariant tangent basis
  struct gkyl_array *triad_basis; // Triad basis
  struct gkyl_array *conf_poisson_tensor; // Configuration space Poisson tensor representation

  // Size of the PT
  int num_pt_indices[3] = { 1 , 6, 18 }; 

  // Allocate arrays for covariant tangent basis 
  cov_tangent_basis = mkarr(basis.num_basis*vdim*vdim, local_ext.volume);
  triad_basis = mkarr(basis.num_basis*vdim*vdim, local_ext.volume);
  conf_poisson_tensor = mkarr(basis.num_basis*num_pt_indices[vdim-1], local_ext.volume);

  // Construct the Geometry for this configuration
  gkyl_vlasov_triad_geom_new(&confGrid, &confLocal, confBasis, 
    &grid, &local, basis, inp_basis_vectors, conf_poisson_tensor);

  // Iterate over the grid, conf space, checking output
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &confLocal);
  while (gkyl_range_iter_next(&iter)) {
    const double *conf_poisson_tensor_d = gkyl_array_cfetch(conf_poisson_tensor, gkyl_range_idx(&confLocal, iter.idx));
    
    int NC = confBasis.num_basis;

    // index size of the symmetric h_ij
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
  gkyl_array_release(cov_tangent_basis);
  gkyl_array_release(triad_basis);
  gkyl_array_release(conf_poisson_tensor);
}

void eval_cov_tangent_basis_annulus_2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double q_r = xn[0];
  double q_theta = 0.0;
  fout[0] = cos(q_theta);
  fout[1] = sin(q_theta);
  fout[2] = - q_r * sin(q_theta);
  fout[3] = q_r * cos(q_theta);
}

void eval_triad_basis_annulus_2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double q_r = xn[0];
  double q_theta = 0.0;
  fout[0] = cos(q_theta);
  fout[1] = sin(q_theta);
  fout[2] = - sin(q_theta);
  fout[3] = cos(q_theta);
}

void eval_triad_basis_gradient_annulus_2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double q_r = xn[0];
  double q_theta = 0.0;
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = 0.0;

  fout[4] = sin(q_theta);
  fout[5] = cos(q_theta);
  fout[6] = - cos(q_theta);
  fout[7] = sin(q_theta);
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
  struct gkyl_vlasov_triad_geom_inp inp_basis_vectors = { 0 };
  inp_basis_vectors.use_vierbein = false;
  inp_basis_vectors.use_preset_geom = false;
  inp_basis_vectors.eval_cov_tangent_basis = eval_cov_tangent_basis_annulus_2v; 
  inp_basis_vectors.eval_triad_basis = eval_triad_basis_annulus_2v; 
  inp_basis_vectors.eval_triad_basis_gradient = eval_triad_basis_gradient_annulus_2v; 
  inp_basis_vectors.eval_cov_tangent_basis_ctx = 0; 
  inp_basis_vectors.eval_triad_basis_ctx = 0; 
  inp_basis_vectors.eval_triad_basis_gradient_ctx = 0; 

  // Make the memory for arrays (modal)
  struct gkyl_array *cov_tangent_basis; // Covariant tangent basis
  struct gkyl_array *triad_basis; // Triad basis
  struct gkyl_array *conf_poisson_tensor; // Configuration space Poisson tensor representation

  // Size of the PT
  int num_pt_indices[3] = { 1 , 6, 18 }; 

  // Allocate arrays for covariant tangent basis 
  cov_tangent_basis = mkarr(basis.num_basis*vdim*vdim, local_ext.volume);
  triad_basis = mkarr(basis.num_basis*vdim*vdim, local_ext.volume);
  conf_poisson_tensor = mkarr(basis.num_basis*num_pt_indices[vdim-1], local_ext.volume);

  // Construct the Geometry for this configuration
  gkyl_vlasov_triad_geom_new(&confGrid, &confLocal, confBasis, 
    &grid, &local, basis, inp_basis_vectors, conf_poisson_tensor);

  // Iterate over the grid, conf space, checking output
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &confLocal);
  while (gkyl_range_iter_next(&iter)) {
    const double *conf_poisson_tensor_d = gkyl_array_cfetch(conf_poisson_tensor, gkyl_range_idx(&confLocal, iter.idx));
    
    int NC = confBasis.num_basis;

    // index size of the symmetric h_ij
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
  gkyl_array_release(cov_tangent_basis);
  gkyl_array_release(triad_basis);
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

  // 1. Test kernel_metric_2v()
  double cov_tangent_basis[4];
  double h_ij[3];

  // Randomly choosen coordinates
  double r = 3.2;
  double theta = 1.4;

  // Assign covaraint tangent basis 
  cov_tangent_basis[0] = cos(theta);
  cov_tangent_basis[1] = sin(theta);
  cov_tangent_basis[2] = - r * sin(theta);
  cov_tangent_basis[3] = r * cos(theta);

  kernel_metric_2v(cov_tangent_basis, h_ij);

  // h_ij only has symmetric entries
  TEST_CHECK(gkyl_compare_double(h_ij[0], 1.0, eps));
  TEST_CHECK(gkyl_compare_double(h_ij[1], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(h_ij[2], r*r, eps));

  // 2. Test kernel_metric_inv_2v()
  double h_ij_inv[3];
  kernel_metric_inv_2v(h_ij, h_ij_inv);

  // h_ij_inv only has symmetric entries as well
  TEST_CHECK(gkyl_compare_double(h_ij_inv[0], 1.0, eps));
  TEST_CHECK(gkyl_compare_double(h_ij_inv[1], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(h_ij_inv[2], 1.0/(r*r), eps));

  // 3. Test kernel_metric_det_2v()
  double det_h[1];
  kernel_metric_det_2v(h_ij, det_h);

  TEST_CHECK(gkyl_compare_double(det_h[0], r, eps));

  // 4. Test compute_nu_inv_2v()

  double triad_basis[4];
  double nu_inv[4];

  // Assign the triad basis 
  triad_basis[0] = cos(theta);
  triad_basis[1] = sin(theta);
  triad_basis[2] = - sin(theta);
  triad_basis[3] = cos(theta);

  compute_nu_inv_2v(h_ij_inv, triad_basis, cov_tangent_basis, nu_inv);

  TEST_CHECK(gkyl_compare_double(nu_inv[0], 1.0, eps));
  TEST_CHECK(gkyl_compare_double(nu_inv[1], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(nu_inv[2], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(nu_inv[3], 1/r, eps));

  // 5. Test kernel_conf_poisson_tensor_2v()
  double conf_poisson_tensor[6];

  // Assign the triad basis gradient
  double triad_basis_gradient[8];

  // gradient in r
  triad_basis_gradient[0] = 0.0;
  triad_basis_gradient[1] = 0.0;

  triad_basis_gradient[2] = 0.0;
  triad_basis_gradient[3] = 0.0;

  // gradient in theta
  triad_basis_gradient[4] = -sin(theta);
  triad_basis_gradient[5] = cos(theta);

  triad_basis_gradient[6] = -cos(theta);
  triad_basis_gradient[7] = -sin(theta);

  kernel_conf_poisson_tensor_2v(h_ij_inv, triad_basis, cov_tangent_basis,
    triad_basis_gradient, conf_poisson_tensor);

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

  // 1. Test kernel_metric_3v()
  double cov_tangent_basis[9];
  double h_ij[6];

  // Randomly choosen coordinates
  double r = 3.2;
  double theta = 1.4;
  double phi = 0.8;

  // Assign covaraint tangent basis 
  cov_tangent_basis[0] = sin(theta) * cos(phi);
  cov_tangent_basis[1] = sin(theta) * sin(phi);
  cov_tangent_basis[2] = cos(theta);

  cov_tangent_basis[3] = r * cos(theta) * cos(phi);
  cov_tangent_basis[4] = r * cos(theta) * sin(phi);
  cov_tangent_basis[5] = -r * sin(theta);

  cov_tangent_basis[6] = - r * sin(theta) * sin(phi);
  cov_tangent_basis[7] = r * sin(theta) * cos(phi);
  cov_tangent_basis[8] = 0.0;

  kernel_metric_3v(cov_tangent_basis, h_ij);

  // h_ij only has symmetric entries
  TEST_CHECK(gkyl_compare_double(h_ij[0], 1.0, eps));
  TEST_CHECK(gkyl_compare_double(h_ij[1], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(h_ij[2], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(h_ij[3], r*r, eps));
  TEST_CHECK(gkyl_compare_double(h_ij[4], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(h_ij[5], r*r*sin(theta)*sin(theta), eps));

  // 2. Test kernel_metric_inv_3v()
  double h_ij_inv[6];
  kernel_metric_inv_3v(h_ij, h_ij_inv);

  // h_ij_inv only has symmetric entries as well
  TEST_CHECK(gkyl_compare_double(h_ij_inv[0], 1.0, eps));
  TEST_CHECK(gkyl_compare_double(h_ij_inv[1], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(h_ij_inv[2], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(h_ij_inv[3], 1.0/(r*r), eps));
  TEST_CHECK(gkyl_compare_double(h_ij_inv[4], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(h_ij_inv[5], 1.0/(r*r*sin(theta)*sin(theta)), eps));

  // 3. Test kernel_metric_det_3v()
  double det_h[1];
  kernel_metric_det_3v(h_ij, det_h);

  TEST_CHECK(gkyl_compare_double(det_h[0], r*r*sin(theta), eps));

  // 4. Test compute_nu_inv_3v()

  double triad_basis[9];
  double nu_inv[9];

  // Assign the triad basis 
  triad_basis[0] = sin(theta) * cos(phi);
  triad_basis[1] = sin(theta) * sin(phi);
  triad_basis[2] = cos(theta);

  triad_basis[3] = cos(theta) * cos(phi);
  triad_basis[4] = cos(theta) * sin(phi);
  triad_basis[5] = - sin(theta);

  triad_basis[6] = - sin(phi);
  triad_basis[7] = cos(phi);
  triad_basis[8] = 0.0;

  compute_nu_inv_3v(h_ij_inv, triad_basis, cov_tangent_basis, nu_inv);

  TEST_CHECK(gkyl_compare_double(nu_inv[0], 1.0, eps));
  TEST_CHECK(gkyl_compare_double(nu_inv[1], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(nu_inv[2], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(nu_inv[3], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(nu_inv[4], 1/r, eps));
  TEST_CHECK(gkyl_compare_double(nu_inv[5], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(nu_inv[6], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(nu_inv[7], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(nu_inv[8], 1/(r*sin(theta)), eps));

  // 5. Test kernel_conf_poisson_tensor_3v()
  double conf_poisson_tensor[18];

  // Assign the triad basis gradient
  double triad_basis_gradient[27];

  // gradient in r
  triad_basis_gradient[0] = 0.0;
  triad_basis_gradient[1] = 0.0;
  triad_basis_gradient[2] = 0.0;

  triad_basis_gradient[3] = 0.0;
  triad_basis_gradient[4] = 0.0;
  triad_basis_gradient[5] = 0.0;

  triad_basis_gradient[6] = 0.0;
  triad_basis_gradient[7] = 0.0;
  triad_basis_gradient[8] = 0.0;

  // gradient in theta
  triad_basis_gradient[9] = cos(phi)*cos(theta);
  triad_basis_gradient[10] = cos(theta)*sin(phi);
  triad_basis_gradient[11] = -sin(theta);

  triad_basis_gradient[12] = -cos(phi)*sin(theta);
  triad_basis_gradient[13] = -sin(phi)*sin(theta);
  triad_basis_gradient[14] = -cos(theta);

  triad_basis_gradient[15] = 0.0;
  triad_basis_gradient[16] = 0.0;
  triad_basis_gradient[17] = 0.0;

  // gradient in phi
  triad_basis_gradient[18] = -sin(phi)*sin(theta);
  triad_basis_gradient[19] = cos(phi)*sin(theta);
  triad_basis_gradient[20] = 0.0;

  triad_basis_gradient[21] = -cos(theta)*sin(phi);
  triad_basis_gradient[22] = cos(phi)*cos(theta);
  triad_basis_gradient[23] = 0.0;

  triad_basis_gradient[24] = -cos(phi);
  triad_basis_gradient[25] = -sin(phi);
  triad_basis_gradient[26] = 0.0;

  kernel_conf_poisson_tensor_3v(h_ij_inv, triad_basis, cov_tangent_basis,
    triad_basis_gradient, conf_poisson_tensor);

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

// Fill per-cell, node-major nodal arrays by sampling the evalf callbacks at the
// conf-basis nodes — the same sampling gkyl_vlasov_triad_geom_from_basis does
// internally, so from_nodal fed these arrays must reproduce it exactly.
static void
fill_nodal_from_evalf(const struct gkyl_rect_grid *cgrid, const struct gkyl_range *crange,
  const struct gkyl_basis *cbasis, int vdim,
  evalf_t eval_cov_tangent_basis, evalf_t eval_triad_basis, evalf_t eval_triad_basis_gradient,
  struct gkyl_array *cov_tangent_basis_nodal, struct gkyl_array *triad_basis_nodal,
  struct gkyl_array *triad_basis_gradient_nodal)
{
  int num_basis = cbasis->num_basis;
  struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, cgrid->ndim, num_basis);
  cbasis->node_list(gkyl_array_fetch(nodes, 0));

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, crange);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_rect_grid_cell_center(cgrid, iter.idx, xc);
    long lidx = gkyl_range_idx(crange, iter.idx);
    double *ctb = gkyl_array_fetch(cov_tangent_basis_nodal, lidx);
    double *tb = gkyl_array_fetch(triad_basis_nodal, lidx);
    double *tbg = gkyl_array_fetch(triad_basis_gradient_nodal, lidx);
    for (int i=0; i<num_basis; ++i) {
      log_to_comp(cgrid->ndim, gkyl_array_cfetch(nodes, i), cgrid->dx, xc, xmu);
      eval_cov_tangent_basis(0.0, xmu, &ctb[i*vdim*vdim], NULL);
      eval_triad_basis(0.0, xmu, &tb[i*vdim*vdim], NULL);
      eval_triad_basis_gradient(0.0, xmu, &tbg[i*vdim*vdim*vdim], NULL);
    }
  }
  gkyl_array_release(nodes);
}

// A/B the from_nodal constructor against the evalf-callback constructor on an
// analytic geometry: identical nodal inputs must give an identical Poisson tensor.
static void
test_triad_from_nodal(int cdim, int vdim, int poly_order,
  const double *confLower, const double *confUpper, const int *confCells,
  evalf_t eval_cov_tangent_basis, evalf_t eval_triad_basis, evalf_t eval_triad_basis_gradient)
{
  int ndim = cdim+vdim;
  int num_pt_indices[3] = { 1 , 6, 18 };
  int num_pt = num_pt_indices[vdim-1];

  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  int cells[GKYL_MAX_DIM];
  for (int d=0; d<cdim; ++d) {
    lower[d] = confLower[d]; upper[d] = confUpper[d]; cells[d] = confCells[d];
  }
  for (int d=cdim; d<ndim; ++d) {
    lower[d] = -1.0; upper[d] = 1.0; cells[d] = 2;
  }

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  struct gkyl_basis basis, confBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[GKYL_MAX_DIM] = { 1, 1, 1 };
  struct gkyl_range confLocal, confLocal_ext;
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<cdim; ++d) ghost[d] = confGhost[d];
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  int NC = confBasis.num_basis;

  // Reference: the evalf-callback constructor
  struct gkyl_vlasov_triad_geom_inp inp_basis_vectors = { 0 };
  inp_basis_vectors.eval_cov_tangent_basis = eval_cov_tangent_basis;
  inp_basis_vectors.eval_triad_basis = eval_triad_basis;
  inp_basis_vectors.eval_triad_basis_gradient = eval_triad_basis_gradient;

  struct gkyl_array *pt_ref = mkarr(NC*num_pt, confLocal_ext.volume);
  gkyl_vlasov_triad_geom_new(&confGrid, &confLocal, confBasis,
    &grid, &local, basis, inp_basis_vectors, pt_ref);

  // from_nodal fed the same callbacks sampled at the conf-basis nodes
  struct gkyl_array *cov_tangent_basis_nodal = mkarr(NC*vdim*vdim, confLocal_ext.volume);
  struct gkyl_array *triad_basis_nodal = mkarr(NC*vdim*vdim, confLocal_ext.volume);
  struct gkyl_array *triad_basis_gradient_nodal = mkarr(NC*vdim*vdim*vdim, confLocal_ext.volume);

  fill_nodal_from_evalf(&confGrid, &confLocal, &confBasis, vdim,
    eval_cov_tangent_basis, eval_triad_basis, eval_triad_basis_gradient,
    cov_tangent_basis_nodal, triad_basis_nodal, triad_basis_gradient_nodal);

  struct gkyl_array *pt_nodal = mkarr(NC*num_pt, confLocal_ext.volume);
  gkyl_vlasov_triad_geom_from_nodal(&confGrid, &confLocal, confBasis,
    &grid, &local, basis,
    cov_tangent_basis_nodal, triad_basis_nodal, triad_basis_gradient_nodal, pt_nodal);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &confLocal);
  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(&confLocal, iter.idx);
    const double *ref_d = gkyl_array_cfetch(pt_ref, lidx);
    const double *nod_d = gkyl_array_cfetch(pt_nodal, lidx);
    for (int k=0; k<NC*num_pt; ++k) {
      // 1e-14, not exact: log_to_comp is a static inline instantiated in both this
      // TU and vlasov_triad_geom.c, and fast-math may round the node coordinate
      // differently between the two before it reaches the trig calls.
      TEST_CHECK( gkyl_compare_double(nod_d[k], ref_d[k], 1e-14) );
      TEST_MSG("cell (%d), coeff %d: from_nodal %1.16e vs from_basis %1.16e",
        iter.idx[0], k, nod_d[k], ref_d[k]);
    }
  }

  gkyl_array_release(pt_ref);
  gkyl_array_release(pt_nodal);
  gkyl_array_release(cov_tangent_basis_nodal);
  gkyl_array_release(triad_basis_nodal);
  gkyl_array_release(triad_basis_gradient_nodal);
}

// Spherical coordinates at fixed (theta, phi) as a function of the single conf
// coordinate r — exercises the vdim=3, 18-component Poisson tensor path.
static const double sph_theta_c = 1.4, sph_phi_c = 0.8;

void eval_cov_tangent_basis_sph_3v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double r = xn[0], theta = sph_theta_c, phi = sph_phi_c;
  fout[0] = sin(theta)*cos(phi);
  fout[1] = sin(theta)*sin(phi);
  fout[2] = cos(theta);

  fout[3] = r*cos(theta)*cos(phi);
  fout[4] = r*cos(theta)*sin(phi);
  fout[5] = -r*sin(theta);

  fout[6] = -r*sin(theta)*sin(phi);
  fout[7] = r*sin(theta)*cos(phi);
  fout[8] = 0.0;
}

void eval_triad_basis_sph_3v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double theta = sph_theta_c, phi = sph_phi_c;
  fout[0] = sin(theta)*cos(phi);
  fout[1] = sin(theta)*sin(phi);
  fout[2] = cos(theta);

  fout[3] = cos(theta)*cos(phi);
  fout[4] = cos(theta)*sin(phi);
  fout[5] = -sin(theta);

  fout[6] = -sin(phi);
  fout[7] = cos(phi);
  fout[8] = 0.0;
}

void eval_triad_basis_gradient_sph_3v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double theta = sph_theta_c, phi = sph_phi_c;
  // gradient in r
  for (int k=0; k<9; ++k) fout[k] = 0.0;

  // gradient in theta
  fout[9] = cos(phi)*cos(theta);
  fout[10] = cos(theta)*sin(phi);
  fout[11] = -sin(theta);

  fout[12] = -cos(phi)*sin(theta);
  fout[13] = -sin(phi)*sin(theta);
  fout[14] = -cos(theta);

  fout[15] = 0.0;
  fout[16] = 0.0;
  fout[17] = 0.0;

  // gradient in phi
  fout[18] = -sin(phi)*sin(theta);
  fout[19] = cos(phi)*sin(theta);
  fout[20] = 0.0;

  fout[21] = -cos(theta)*sin(phi);
  fout[22] = cos(phi)*cos(theta);
  fout[23] = 0.0;

  fout[24] = -cos(phi);
  fout[25] = -sin(phi);
  fout[26] = 0.0;
}

// A/B the interior-node (Gauss-Legendre) from_nodal variant: fill the nodal
// arrays at the GL node coordinates, then check the modal Poisson tensor
// against the pointwise kernels + quad_nodal_to_modal computed directly here.
// This pins down the global-nodal-range gather (node->cell mapping and node
// ordering), which is the part the gyrokinetic bridge relies on.
static void
test_triad_from_nodal_interior(int cdim, int vdim,
  const double *confLower, const double *confUpper, const int *confCells,
  evalf_t eval_cov_tangent_basis, evalf_t eval_triad_basis, evalf_t eval_triad_basis_gradient)
{
  int poly_order = 1;
  int ndim = cdim+vdim;
  int num_pt_indices[3] = { 1 , 6, 18 };
  int num_pt = num_pt_indices[vdim-1];

  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  int cells[GKYL_MAX_DIM];
  for (int d=0; d<cdim; ++d) {
    lower[d] = confLower[d]; upper[d] = confUpper[d]; cells[d] = confCells[d];
  }
  for (int d=cdim; d<ndim; ++d) {
    lower[d] = -1.0; upper[d] = 1.0; cells[d] = 2;
  }

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  struct gkyl_basis basis, confBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  struct gkyl_range confLocal, confLocal_ext;
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<cdim; ++d) ghost[d] = confGhost[d];
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  int NC = confBasis.num_basis; // == number of GL interior nodes per cell at p1

  // Global interior-node range: 2 nodes per direction per cell, zero-based.
  int nlower[GKYL_MAX_CDIM] = { 0 }, nupper[GKYL_MAX_CDIM] = { 0 };
  for (int d=0; d<cdim; ++d) nupper[d] = 2*confCells[d]-1;
  struct gkyl_range nrange;
  gkyl_range_init(&nrange, cdim, nlower, nupper);

  struct gkyl_array *ctb_n = mkarr(vdim*vdim, nrange.volume);
  struct gkyl_array *tb_n = mkarr(vdim*vdim, nrange.volume);
  struct gkyl_array *tbg_n = mkarr(vdim*vdim*vdim, nrange.volume);

  double glpt = 1.0/sqrt(3.0);
  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];
  int nidx[GKYL_MAX_CDIM];
  struct gkyl_range_iter citer;
  gkyl_range_iter_init(&citer, &confLocal);
  while (gkyl_range_iter_next(&citer)) {
    gkyl_rect_grid_cell_center(&confGrid, citer.idx, xc);
    for (int i=0; i<NC; ++i) {
      for (int j=0; j<cdim; ++j) {
        int b = (i >> (cdim-1-j)) & 1;
        nidx[j] = (citer.idx[j]-confLocal.lower[j])*2 + b;
        xmu[j] = xc[j] + (2*b-1)*0.5*confGrid.dx[j]*glpt;
      }
      long ln = gkyl_range_idx(&nrange, nidx);
      eval_cov_tangent_basis(0.0, xmu, gkyl_array_fetch(ctb_n, ln), NULL);
      eval_triad_basis(0.0, xmu, gkyl_array_fetch(tb_n, ln), NULL);
      eval_triad_basis_gradient(0.0, xmu, gkyl_array_fetch(tbg_n, ln), NULL);
    }
  }

  struct gkyl_array *pt_out = mkarr(NC*num_pt, confLocal_ext.volume);
  gkyl_vlasov_triad_geom_from_nodal_interior(&confGrid, &confLocal, confBasis,
    &grid, &local, basis, &nrange, ctb_n, tb_n, tbg_n, pt_out);

  // Reference: pointwise kernels at the GL nodes + quad_nodal_to_modal, straight
  // from the callbacks with no global-nodal-range in between.
  metric_t kern_metric = choose_metric_kern(vdim);
  metric_inv_t kern_metric_inv = choose_metric_inv_kern(vdim);
  conf_poisson_tensor_t kern_pt = choose_conf_poisson_tensor_kern(vdim);

  double ctb[9], tb[9], tbg[27], h_ij[6], h_ij_inv[6];
  double pt_at_nodes[8*18], fnodal[8], pt_expected[8*18];

  gkyl_range_iter_init(&citer, &confLocal);
  while (gkyl_range_iter_next(&citer)) {
    gkyl_rect_grid_cell_center(&confGrid, citer.idx, xc);
    for (int i=0; i<NC; ++i) {
      for (int j=0; j<cdim; ++j) {
        int b = (i >> (cdim-1-j)) & 1;
        xmu[j] = xc[j] + (2*b-1)*0.5*confGrid.dx[j]*glpt;
      }
      eval_cov_tangent_basis(0.0, xmu, ctb, NULL);
      eval_triad_basis(0.0, xmu, tb, NULL);
      eval_triad_basis_gradient(0.0, xmu, tbg, NULL);
      kern_metric(ctb, h_ij);
      kern_metric_inv(h_ij, h_ij_inv);
      kern_pt(h_ij_inv, tb, ctb, tbg, &pt_at_nodes[i*num_pt]);
    }
    for (int c=0; c<num_pt; ++c) {
      for (int i=0; i<NC; ++i)
        fnodal[i] = pt_at_nodes[i*num_pt + c];
      for (int k=0; k<NC; ++k)
        confBasis.quad_nodal_to_modal(fnodal, &pt_expected[c*NC], k);
    }

    const double *out_d = gkyl_array_cfetch(pt_out, gkyl_range_idx(&confLocal, citer.idx));
    for (int k=0; k<NC*num_pt; ++k) {
      // 1e-13, not exact: fast-math compiles the callback/kernel chains of the two
      // loops differently, and the analytic phi-cancellations in the vdim=3 Pi
      // contractions then differ at a few-ulp level between them.
      TEST_CHECK( gkyl_compare_double(out_d[k], pt_expected[k], 1e-13) );
      TEST_MSG("cell (%d,%d,%d), coeff %d: got %1.16e expected %1.16e",
        citer.idx[0], cdim>1 ? citer.idx[1] : 0, cdim>2 ? citer.idx[2] : 0,
        k, out_d[k], pt_expected[k]);
    }
  }

  gkyl_array_release(ctb_n);
  gkyl_array_release(tb_n);
  gkyl_array_release(tbg_n);
  gkyl_array_release(pt_out);
}

// Spherical-form fields varying in ALL conf directions (r=x, theta(y), phi(z)) —
// any transposed node gather in the interior bridge shows up as a mismatch.
void eval_cov_tangent_basis_sph_3x3v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double r = xn[0], theta = 1.0 + 0.2*xn[1], phi = 0.5 + 0.3*xn[2];
  fout[0] = sin(theta)*cos(phi);
  fout[1] = sin(theta)*sin(phi);
  fout[2] = cos(theta);

  fout[3] = r*cos(theta)*cos(phi);
  fout[4] = r*cos(theta)*sin(phi);
  fout[5] = -r*sin(theta);

  fout[6] = -r*sin(theta)*sin(phi);
  fout[7] = r*sin(theta)*cos(phi);
  fout[8] = 0.0;
}

void eval_triad_basis_sph_3x3v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double theta = 1.0 + 0.2*xn[1], phi = 0.5 + 0.3*xn[2];
  fout[0] = sin(theta)*cos(phi);
  fout[1] = sin(theta)*sin(phi);
  fout[2] = cos(theta);

  fout[3] = cos(theta)*cos(phi);
  fout[4] = cos(theta)*sin(phi);
  fout[5] = -sin(theta);

  fout[6] = -sin(phi);
  fout[7] = cos(phi);
  fout[8] = 0.0;
}

void eval_triad_basis_gradient_sph_3x3v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double theta = 1.0 + 0.2*xn[1], phi = 0.5 + 0.3*xn[2];
  // gradient in r
  for (int k=0; k<9; ++k) fout[k] = 0.0;

  // gradient in theta
  fout[9] = cos(phi)*cos(theta);
  fout[10] = cos(theta)*sin(phi);
  fout[11] = -sin(theta);

  fout[12] = -cos(phi)*sin(theta);
  fout[13] = -sin(phi)*sin(theta);
  fout[14] = -cos(theta);

  fout[15] = 0.0;
  fout[16] = 0.0;
  fout[17] = 0.0;

  // gradient in phi
  fout[18] = -sin(phi)*sin(theta);
  fout[19] = cos(phi)*sin(theta);
  fout[20] = 0.0;

  fout[21] = -cos(theta)*sin(phi);
  fout[22] = cos(phi)*cos(theta);
  fout[23] = 0.0;

  fout[24] = -cos(phi);
  fout[25] = -sin(phi);
  fout[26] = 0.0;
}

void
test_triad_from_nodal_interior_1x2v_annulus_p1()
{
  double confLower[] = {0.1}, confUpper[] = {1.0};
  int confCells[] = {4};
  test_triad_from_nodal_interior(1, 2, confLower, confUpper, confCells,
    eval_cov_tangent_basis_annulus_2v, eval_triad_basis_annulus_2v, eval_triad_basis_gradient_annulus_2v);
}

void
test_triad_from_nodal_interior_3x3v_spherical_p1()
{
  double confLower[] = {0.1, 0.0, 0.0}, confUpper[] = {1.0, 1.0, 1.0};
  int confCells[] = {3, 2, 2};
  test_triad_from_nodal_interior(3, 3, confLower, confUpper, confCells,
    eval_cov_tangent_basis_sph_3x3v, eval_triad_basis_sph_3x3v, eval_triad_basis_gradient_sph_3x3v);
}

// --- from_tangents_interior (b-aligned Gram-Schmidt + FD gradients) tests ---

typedef void (*nodal_fill_t)(const double *xmu, double *out);

// Cylindrical field-aligned coordinates (r, alpha, z), field along z:
// tangents e_r, e_alpha, e_z; the b-aligned GS triad is exactly (rhat, alphahat, zhat).
static void fill_cyl_tangents(const double *xmu, double *tan)
{
  double r = xmu[0], a = xmu[1];
  tan[0] = cos(a);    tan[1] = sin(a);   tan[2] = 0.0;
  tan[3] = -r*sin(a); tan[4] = r*cos(a); tan[5] = 0.0;
  tan[6] = 0.0;       tan[7] = 0.0;      tan[8] = 1.0;
}

static void fill_cyl_triad(const double *xmu, double *tri)
{
  double a = xmu[1];
  tri[0] = cos(a);  tri[1] = sin(a); tri[2] = 0.0;
  tri[3] = -sin(a); tri[4] = cos(a); tri[5] = 0.0;
  tri[6] = 0.0;     tri[7] = 0.0;    tri[8] = 1.0;
}

static void fill_cyl_triad_grad(const double *xmu, double *grad)
{
  double a = xmu[1];
  for (int k=0; k<27; ++k) grad[k] = 0.0;
  // d/dalpha of (rhat, alphahat, zhat)
  grad[9] = -sin(a);  grad[10] = cos(a);
  grad[12] = -cos(a); grad[13] = -sin(a);
}

static void fill_cyl_bhat(const double *xmu, double *b)
{
  b[0] = 0.0; b[1] = 0.0; b[2] = 1.0;
}

static void fill_flat_tangents(const double *xmu, double *tan)
{
  for (int k=0; k<9; ++k) tan[k] = 0.0;
  tan[0] = 1.0; tan[4] = 1.0; tan[8] = 1.0;
}

static void fill_flat_triad_grad(const double *xmu, double *grad)
{
  for (int k=0; k<27; ++k) grad[k] = 0.0;
}

// A/B the b-aligned construction against the interior constructor fed the
// ANALYTIC triad and gradients: the difference is the GS reconstruction (exact
// for these geometries) plus the finite-difference gradient error.
static void
test_triad_from_tangents_interior(nodal_fill_t fill_tan, nodal_fill_t fill_tri,
  nodal_fill_t fill_grad, nodal_fill_t fill_bhat,
  const double *confLower, const double *confUpper, const int *confCells, double tol)
{
  int cdim = 3, vdim = 3, poly_order = 1;
  int ndim = cdim+vdim;
  int num_pt = 18;

  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  int cells[GKYL_MAX_DIM];
  for (int d=0; d<cdim; ++d) {
    lower[d] = confLower[d]; upper[d] = confUpper[d]; cells[d] = confCells[d];
  }
  for (int d=cdim; d<ndim; ++d) {
    lower[d] = -1.0; upper[d] = 1.0; cells[d] = 2;
  }

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  struct gkyl_basis basis, confBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  struct gkyl_range confLocal, confLocal_ext;
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<cdim; ++d) ghost[d] = confGhost[d];
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  int NC = confBasis.num_basis;

  int nlower[GKYL_MAX_CDIM] = { 0 }, nupper[GKYL_MAX_CDIM] = { 0 };
  for (int d=0; d<cdim; ++d) nupper[d] = 2*confCells[d]-1;
  struct gkyl_range nrange;
  gkyl_range_init(&nrange, cdim, nlower, nupper);

  struct gkyl_array *tan_n = mkarr(9, nrange.volume);
  struct gkyl_array *bhat_n = mkarr(3, nrange.volume);
  struct gkyl_array *tri_n = mkarr(9, nrange.volume);
  struct gkyl_array *grad_n = mkarr(27, nrange.volume);

  double glpt = 1.0/sqrt(3.0);
  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];
  int nidx[GKYL_MAX_CDIM];
  struct gkyl_range_iter citer;
  gkyl_range_iter_init(&citer, &confLocal);
  while (gkyl_range_iter_next(&citer)) {
    gkyl_rect_grid_cell_center(&confGrid, citer.idx, xc);
    for (int i=0; i<NC; ++i) {
      for (int j=0; j<cdim; ++j) {
        int b = (i >> (cdim-1-j)) & 1;
        nidx[j] = (citer.idx[j]-confLocal.lower[j])*2 + b;
        xmu[j] = xc[j] + (2*b-1)*0.5*confGrid.dx[j]*glpt;
      }
      long ln = gkyl_range_idx(&nrange, nidx);
      fill_tan(xmu, gkyl_array_fetch(tan_n, ln));
      fill_bhat(xmu, gkyl_array_fetch(bhat_n, ln));
      fill_tri(xmu, gkyl_array_fetch(tri_n, ln));
      fill_grad(xmu, gkyl_array_fetch(grad_n, ln));
    }
  }

  struct gkyl_array *pt_a = mkarr(NC*num_pt, confLocal_ext.volume);
  gkyl_vlasov_triad_geom_from_tangents_interior(&confGrid, &confLocal, confBasis,
    &grid, &local, basis, &nrange, tan_n, bhat_n, true, pt_a);

  struct gkyl_array *pt_b = mkarr(NC*num_pt, confLocal_ext.volume);
  gkyl_vlasov_triad_geom_from_nodal_interior(&confGrid, &confLocal, confBasis,
    &grid, &local, basis, &nrange, tan_n, tri_n, grad_n, pt_b);

  gkyl_range_iter_init(&citer, &confLocal);
  while (gkyl_range_iter_next(&citer)) {
    long lidx = gkyl_range_idx(&confLocal, citer.idx);
    const double *a_d = gkyl_array_cfetch(pt_a, lidx);
    const double *b_d = gkyl_array_cfetch(pt_b, lidx);
    for (int k=0; k<NC*num_pt; ++k) {
      TEST_CHECK( gkyl_compare_double(a_d[k], b_d[k], tol) );
      TEST_MSG("cell (%d,%d,%d), coeff %d: from_tangents %1.16e vs analytic %1.16e",
        citer.idx[0], citer.idx[1], citer.idx[2], k, a_d[k], b_d[k]);
    }
  }

  gkyl_array_release(tan_n);
  gkyl_array_release(bhat_n);
  gkyl_array_release(tri_n);
  gkyl_array_release(grad_n);
  gkyl_array_release(pt_a);
  gkyl_array_release(pt_b);
}

void
test_triad_from_tangents_interior_3x3v_flat()
{
  double confLower[] = {1.0, 0.0, 0.0}, confUpper[] = {2.0, 0.4, 1.0};
  int confCells[] = {2, 2, 2};
  // Constant tangents: FD gradients are exactly zero, so the match is exact.
  test_triad_from_tangents_interior(fill_flat_tangents, fill_flat_tangents,
    fill_flat_triad_grad, fill_cyl_bhat, confLower, confUpper, confCells, 1e-14);
}

void
test_triad_from_tangents_interior_3x3v_cylindrical()
{
  double confLower[] = {1.0, 0.0, 0.0}, confUpper[] = {2.0, 0.4, 1.0};
  int confCells[] = {4, 8, 4};
  // The triad varies as trig in alpha, so the second-order FD gradients carry
  // an O(h^2) truncation error relative to the analytic connection terms.
  test_triad_from_tangents_interior(fill_cyl_tangents, fill_cyl_triad,
    fill_cyl_triad_grad, fill_cyl_bhat, confLower, confUpper, confCells, 1e-3);
}

void
test_triad_from_nodal_1x2v_annulus_p1()
{
  double confLower[] = {0.1}, confUpper[] = {1.0};
  int confCells[] = {4};
  test_triad_from_nodal(1, 2, 1, confLower, confUpper, confCells,
    eval_cov_tangent_basis_annulus_2v, eval_triad_basis_annulus_2v, eval_triad_basis_gradient_annulus_2v);
}

void
test_triad_from_nodal_1x2v_annulus_p2()
{
  double confLower[] = {0.1}, confUpper[] = {1.0};
  int confCells[] = {4};
  test_triad_from_nodal(1, 2, 2, confLower, confUpper, confCells,
    eval_cov_tangent_basis_annulus_2v, eval_triad_basis_annulus_2v, eval_triad_basis_gradient_annulus_2v);
}

void
test_triad_from_nodal_1x3v_spherical_p1()
{
  double confLower[] = {0.1}, confUpper[] = {1.0};
  int confCells[] = {4};
  test_triad_from_nodal(1, 3, 1, confLower, confUpper, confCells,
    eval_cov_tangent_basis_sph_3v, eval_triad_basis_sph_3v, eval_triad_basis_gradient_sph_3v);
}

void test_triad_1v() { test_triad_math_1v(); }
void test_triad_1x1v_flat() { test_triad_1x1v_flat_conf(2); }
void test_triad_1x2v_flat() { test_triad_1x2v_flat_conf(2); }
void test_triad_1x2v_annulus() { test_triad_1x2v_annulus_conf(2); }
void test_triad_2v() { test_triad_math_2v(); }
void test_triad_3v() { test_triad_math_3v(); }

TEST_LIST = {
  { "test_triad_1v", test_triad_1v}, 
  { "test_triad_1x1v_flat", test_triad_1x1v_flat}, 
  { "test_triad_1x2v_flat", test_triad_1x2v_flat}, 
  { "test_triad_1x2v_annulus", test_triad_1x2v_annulus}, 
  { "test_triad_from_nodal_1x2v_annulus_p1", test_triad_from_nodal_1x2v_annulus_p1},
  { "test_triad_from_nodal_1x2v_annulus_p2", test_triad_from_nodal_1x2v_annulus_p2},
  { "test_triad_from_nodal_1x3v_spherical_p1", test_triad_from_nodal_1x3v_spherical_p1},
  { "test_triad_from_nodal_interior_1x2v_annulus_p1", test_triad_from_nodal_interior_1x2v_annulus_p1},
  { "test_triad_from_nodal_interior_3x3v_spherical_p1", test_triad_from_nodal_interior_3x3v_spherical_p1},
  { "test_triad_from_tangents_interior_3x3v_flat", test_triad_from_tangents_interior_3x3v_flat},
  { "test_triad_from_tangents_interior_3x3v_cylindrical", test_triad_from_tangents_interior_3x3v_cylindrical},
  { "test_triad_2v", test_triad_2v},
  { "test_triad_3v", test_triad_3v},
  {NULL, NULL}
};
