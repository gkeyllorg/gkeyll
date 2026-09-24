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


void
proj_func(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0];
  double z = xn[1];
  fout[0] = z*cos(x);
}

void
proj_func2(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0];
  double z = xn[1];
  fout[0] = sin((2.*M_PI/(2.*M_PI))*x);
}

void evalFunc1x_neumannx_dirichletx(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double a = 5.0;
  double c0 = 0.;
  double c1 = a/12. - 1./2.;
  fout[0] = -(1.-a*pow(x,2));
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
test_deflate_inflate(bool use_gpu)
{
  // Create the 2d field.
  // Create xz grid.
  double lower[] = { -M_PI, 0.0 }, upper[] = { M_PI, 1.0 };
  int cells[] = { 12, 8 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  // Ranges.
  struct gkyl_range local, local_ext;
  int nghost[GKYL_MAX_CDIM] = { 1, 1 };
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);

  // Basis function.
  int poly_order = 1;
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 2, poly_order);
  struct gkyl_basis *basis_on_dev = use_gpu? gkyl_cart_modal_serendip_cu_dev_new(2, poly_order)
                                           : gkyl_cart_modal_serendip_new(2, poly_order);

  // Project initial function on 2d field.
  struct gkyl_array *field = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_eval_on_nodes *proj = gkyl_eval_on_nodes_new(&grid, &basis, 1, &proj_func, 0);
  gkyl_eval_on_nodes_advance(proj, 0.0, &local, field);
  gkyl_eval_on_nodes_release(proj);
  //gkyl_grid_sub_array_write(&grid, &local, 0, field, "in_field.gkyl");
  struct gkyl_array *field_dev = use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume)
	                                 : gkyl_array_acquire(field);
  gkyl_array_copy(field_dev, field);

  // Create deflated 1d grid, ranges, basis, and field.
  double deflated_lower[] = { -M_PI}, deflated_upper[] = { M_PI};
  int deflated_cells[] = { 12};
  struct gkyl_rect_grid deflated_grid;
  gkyl_rect_grid_init(&deflated_grid, 1, deflated_lower, deflated_upper, deflated_cells);

  // Ranges.
  struct gkyl_range deflated_local, deflated_local_ext;
  int deflated_nghost[GKYL_MAX_CDIM] = { 1 };
  gkyl_create_grid_ranges(&deflated_grid, deflated_nghost, &deflated_local_ext, &deflated_local);

  // Deflated basis function.
  int deflated_poly_order = 1;
  struct gkyl_basis deflated_basis;
  gkyl_cart_modal_serendip(&deflated_basis, 1, deflated_poly_order);
  struct gkyl_basis *deflated_basis_on_dev = use_gpu? gkyl_cart_modal_serendip_cu_dev_new(1, poly_order)
                                                    : gkyl_cart_modal_serendip_new(1, poly_order);
  
  // Deflated field.
  struct gkyl_array *deflated_field = gkyl_array_new(GKYL_DOUBLE, deflated_basis.num_basis, deflated_local_ext.volume);
  struct gkyl_array *deflated_field_dev = use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, deflated_basis.num_basis, deflated_local_ext.volume)
	                                          : gkyl_array_acquire(deflated_field);

  // Create nrange and the 2d nodal array to be populated.
  int nodes[2] = { 1, 1 };
  if (poly_order == 1){
    for (int d=0; d<grid.ndim; ++d)
      nodes[d] = grid.cells[d] + 1;
  }
  if (poly_order == 2){
    for (int d=0; d<grid.ndim; ++d)
      nodes[d] = 2*(grid.cells[d]) + 1;
  }
  struct gkyl_range nrange;
  gkyl_range_init_from_shape(&nrange, grid.ndim, nodes);
  struct gkyl_array* nodal_fld = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
  struct gkyl_array* nodal_fld_dev = use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, grid.ndim, nrange.volume)
	                                     : gkyl_array_acquire(nodal_fld);

  // Create the deflated nodal range and 1d array that will be used as an intermediate.
  int deflated_nodes[3] = { 1, 1, 1 };
  if (poly_order == 1){
    for (int d=0; d<deflated_grid.ndim; ++d)
      deflated_nodes[d] = deflated_grid.cells[d] + 1;
  }
  if (poly_order == 2){
    for (int d=0; d<deflated_grid.ndim; ++d)
      deflated_nodes[d] = 2*(deflated_grid.cells[d]) + 1;
  }
  struct gkyl_range deflated_nrange;
  gkyl_range_init_from_shape(&deflated_nrange, deflated_grid.ndim, deflated_nodes);

  // Now for each z slice we want to deflate and fill the 2d nodal array.
  gkyl_deflate_zsurf *deflator_lo = gkyl_deflate_zsurf_new(&basis, &deflated_basis, 0, use_gpu);
  gkyl_deflate_zsurf *deflator_up = gkyl_deflate_zsurf_new(&basis, &deflated_basis, 1, use_gpu);
  
  struct gkyl_nodal_ops *n2m_1d = gkyl_nodal_ops_new(&deflated_basis, &deflated_grid, use_gpu);

  int ctr = 0;
  for(int zidx = local.lower[1]; zidx <= local.upper[1]; zidx++){
    // First deflate.
    gkyl_deflate_zsurf_advance(deflator_lo, zidx, &local, &deflated_local, field_dev, deflated_field_dev, 1);
    // Modal to Nodal in 1d -> Store the result in the 2d nodal field.
    gkyl_nodal_ops_m2n_deflated(n2m_1d, deflated_basis_on_dev, 
      &deflated_grid, &nrange, &deflated_nrange, &deflated_local, 1, 
      nodal_fld_dev, deflated_field_dev, ctr);
    ctr += 1;
    if (zidx == local.upper[1]){
      gkyl_deflate_zsurf_advance(deflator_up, zidx, &local, &deflated_local, field_dev, deflated_field_dev, 1);
      gkyl_nodal_ops_m2n_deflated(n2m_1d, deflated_basis_on_dev, 
        &deflated_grid, &nrange, &deflated_nrange, &deflated_local, 1, 
        nodal_fld_dev, deflated_field_dev, ctr);
    }
  }

  struct gkyl_array *out_field = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  struct gkyl_array *out_field_dev = use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume) : gkyl_array_acquire(out_field);

  // Convert back to modal and do a check.
  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&basis, &grid, use_gpu);
  gkyl_nodal_ops_n2m(n2m, basis_on_dev, &grid, &nrange, &local, 1, nodal_fld_dev, out_field_dev, false);
  gkyl_array_copy(out_field, out_field_dev);
  //gkyl_grid_sub_array_write(&grid, &local, 0, out_field, "out_field.gkyl");

  check_same(local, basis, field, out_field);

  if (use_gpu) {
    gkyl_cart_modal_basis_release_cu(basis_on_dev);
    gkyl_cart_modal_basis_release_cu(deflated_basis_on_dev);
  }
  else {
    gkyl_cart_modal_basis_release(basis_on_dev);
    gkyl_cart_modal_basis_release(deflated_basis_on_dev);
  }

  gkyl_array_release(field);
  gkyl_array_release(out_field);
  gkyl_array_release(deflated_field);
  gkyl_array_release(nodal_fld);
  gkyl_array_release(field_dev);
  gkyl_array_release(out_field_dev);
  gkyl_array_release(deflated_field_dev);
  gkyl_array_release(nodal_fld_dev);
  
  gkyl_nodal_ops_release(n2m_1d);
  gkyl_nodal_ops_release(n2m);
  gkyl_deflate_zsurf_release(deflator_lo);
  gkyl_deflate_zsurf_release(deflator_up);
}

void
test_poisson_slices()
{
  // Create the 2d field.
  // Create xz grid.
  double lower[] = { -M_PI, 0.0 }, upper[] = { 3*M_PI/2, 1.0 };
  int cells[] = { 12, 8 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  // Ranges.
  struct gkyl_range local, local_ext;
  int nghost[GKYL_MAX_CDIM] = { 1, 1 };
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);

  // Basis function.
  int poly_order = 1;
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 2, poly_order);

  // Project initial function on 2d field.
  struct gkyl_array *field = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_eval_on_nodes *proj = gkyl_eval_on_nodes_new(&grid, &basis, 1, &proj_func, 0);
  gkyl_eval_on_nodes_advance(proj, 0.0, &local, field);
  gkyl_eval_on_nodes_release(proj);
  //gkyl_grid_sub_array_write(&grid, &local, 0, field, "in_field.gkyl");

  // Create deflated 1d grid, ranges, basis, and field.
  // Create xz grid.
  double deflated_lower[] = { -M_PI}, deflated_upper[] = { M_PI};
  int deflated_cells[] = { 12};
  struct gkyl_rect_grid deflated_grid;
  gkyl_rect_grid_init(&deflated_grid, 1, deflated_lower, deflated_upper, deflated_cells);

  // Ranges,
  struct gkyl_range deflated_local, deflated_local_ext;
  int deflated_nghost[GKYL_MAX_CDIM] = { 1 };
  gkyl_create_grid_ranges(&deflated_grid, deflated_nghost, &deflated_local_ext, &deflated_local);

  // Basis function.
  int deflated_poly_order = 1;
  struct gkyl_basis deflated_basis;
  gkyl_cart_modal_serendip(&deflated_basis, 1, deflated_poly_order);

  // Field.
  struct gkyl_array *deflated_field = gkyl_array_new(GKYL_DOUBLE, deflated_basis.num_basis, deflated_local_ext.volume);
  struct gkyl_array *deflated_phi = gkyl_array_new(GKYL_DOUBLE, deflated_basis.num_basis, deflated_local_ext.volume);

  // Create nrange and the 2d nodal array to be populated.
  int nodes[2] = { 1, 1 };
  if (poly_order == 1){
    for (int d=0; d<grid.ndim; ++d)
      nodes[d] = grid.cells[d] + 1;
  }

  if (poly_order == 2){
    for (int d=0; d<grid.ndim; ++d)
      nodes[d] = 2*(grid.cells[d]) + 1;
  }

  struct gkyl_range nrange;
  gkyl_range_init_from_shape(&nrange, grid.ndim, nodes);
  struct gkyl_array* nodal_fld = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange.volume);

  // Create the deflated nodal range and 1d array that will be used as an intermediate.
  int deflated_nodes[3] = { 1, 1, 1 };
  if (poly_order == 1){
    for (int d=0; d<deflated_grid.ndim; ++d)
      deflated_nodes[d] = deflated_grid.cells[d] + 1;
  }

  if (poly_order == 2){
    for (int d=0; d<deflated_grid.ndim; ++d)
      deflated_nodes[d] = 2*(deflated_grid.cells[d]) + 1;
  }

  struct gkyl_range deflated_nrange;
  gkyl_range_init_from_shape(&deflated_nrange, deflated_grid.ndim, deflated_nodes);
  struct gkyl_array* deflated_nodal_fld = gkyl_array_new(GKYL_DOUBLE, deflated_grid.ndim, deflated_nrange.volume);

  // Now for each z slice we want to deflate
  // and then to a m2n to give us a nodal array at that z slice
  // Then fill the correct nodes in the 2d nodal array
  gkyl_deflate_zsurf *deflator_lo = gkyl_deflate_zsurf_new(&basis, &deflated_basis, 0, false);
  gkyl_deflate_zsurf *deflator_up = gkyl_deflate_zsurf_new(&basis, &deflated_basis, 1, false);

  struct gkyl_poisson_bc poisson_bc;
  poisson_bc.lo_type[0] = GKYL_POISSON_NEUMANN;
  poisson_bc.up_type[0] = GKYL_POISSON_DIRICHLET;
  poisson_bc.lo_value[0].v[0] = 0.;
  poisson_bc.up_value[0].v[0] = 0.;

  struct gkyl_array *epsilon = gkyl_array_new(GKYL_DOUBLE, deflated_basis.num_basis, deflated_local_ext.volume);
  gkyl_array_shiftc(epsilon, sqrt(2.0), 0); // Sets weight=1.
  struct gkyl_fem_poisson *fem_poisson = gkyl_fem_poisson_new(&deflated_local, &deflated_grid,
    deflated_basis, &poisson_bc, NULL, epsilon, 0, false, false);

  struct gkyl_nodal_ops *n2m_1d = gkyl_nodal_ops_new(&basis, &grid, false);
  int nidx[2];
  for(int zidx = local.lower[1]; zidx <= local.upper[1]; zidx++){
    // first deflate
    gkyl_deflate_zsurf_advance(deflator_lo, zidx, &local, &deflated_local, field, deflated_field, 1);
    // do the poisson solve
    gkyl_fem_poisson_set_rhs(fem_poisson, deflated_field, NULL);
    gkyl_fem_poisson_solve(fem_poisson, deflated_phi);
    // then nodal to modal
    gkyl_nodal_ops_m2n(n2m_1d, &deflated_basis, &deflated_grid, &deflated_nrange, &deflated_local, 1, deflated_nodal_fld, deflated_phi, false);
    // now loop through the 2d nrange and populate
    nidx[1] = zidx-1;
    for(int ix = 0; ix <=nrange.upper[0]; ix++){
      nidx[0] = ix;
      long lin_nidx = gkyl_range_idx(&nrange, nidx);
      long lin_nidx_deflated = gkyl_range_idx(&deflated_nrange, &ix);
      const double* input = gkyl_array_fetch(deflated_nodal_fld, lin_nidx_deflated);
      double* output = gkyl_array_fetch(nodal_fld, lin_nidx);
      output[0] = input[0];
    }
    if (zidx == local.upper[1]){
      // first deflate
      gkyl_deflate_zsurf_advance(deflator_up, zidx, &local, &deflated_local, field, deflated_field, 1);
      // do the poisson solve
      gkyl_fem_poisson_set_rhs(fem_poisson, deflated_field, NULL);
      gkyl_fem_poisson_solve(fem_poisson, deflated_phi);
      // then nodal to modal
      gkyl_nodal_ops_m2n(n2m_1d, &deflated_basis, &deflated_grid, &deflated_nrange, &deflated_local, 1, deflated_nodal_fld, deflated_phi, false);
      // now loop through the 2d nrange and populate
      nidx[1] = zidx;
      for(int ix = 0; ix <= nrange.upper[0]; ix++){
        nidx[0] = ix;
        long lin_nidx = gkyl_range_idx(&nrange, nidx);
        long lin_nidx_deflated = gkyl_range_idx(&deflated_nrange, &ix);
        const double* input = gkyl_array_fetch(deflated_nodal_fld, lin_nidx_deflated);
        double* output = gkyl_array_fetch(nodal_fld, lin_nidx);
        output[0] = input[0];
      }
    }

  }

  struct gkyl_array *out_field = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&basis, &grid, false);
  gkyl_nodal_ops_n2m(n2m, &basis, &grid, &nrange, &local, 1, nodal_fld, out_field, false);
  gkyl_nodal_ops_release(n2m);
  //gkyl_grid_sub_array_write(&grid, &local, 0, out_field, "out_field.gkyl");

  // Release allocated resources
  gkyl_array_release(field);
  gkyl_array_release(deflated_field);
  gkyl_array_release(deflated_phi);
  gkyl_array_release(nodal_fld);
  gkyl_array_release(deflated_nodal_fld);
  gkyl_array_release(out_field);
  gkyl_array_release(epsilon);
  
  gkyl_nodal_ops_release(n2m_1d);
  gkyl_deflate_zsurf_release(deflator_lo);
  gkyl_deflate_zsurf_release(deflator_up);
  gkyl_fem_poisson_release(fem_poisson);
}

void test_deflate_inflate_ho(void) { test_deflate_inflate(false); }
void test_deflate_inflate_dev(void) { test_deflate_inflate(true); }

TEST_LIST = {
  { "test_deflate_inflate_ho", test_deflate_inflate_ho},
  { "test_poisson_slices", test_poisson_slices},
#ifdef GKYL_HAVE_CUDA
  { "test_deflate_inflate_dev", test_deflate_inflate_dev},
#endif
  { NULL, NULL },
};