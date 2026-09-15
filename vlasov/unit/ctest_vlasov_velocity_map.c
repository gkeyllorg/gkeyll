#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_velocity_map.h>

// Different polynomial mappings in each velocity dimension so we can
// non-trivially check that the per-direction data is stored in the correct
// order. All are monotone on [-2, 2].
static double map_quad(double v) { return v + 0.1*v*v; }
static double map_quad_prime(double v) { return 1.0 + 0.2*v; }
static double map_cubic(double v) { return v + 0.05*v*v*v; }
static double map_quartic(double v) { return v + 0.01*v*v*v*v; }

static void
eval_map_quad(double t, const double *vc, double *vp, void *ctx)
{
  vp[0] = map_quad(vc[0]);
}
static void
eval_map_cubic(double t, const double *vc, double *vp, void *ctx)
{
  vp[0] = map_cubic(vc[0]);
}
static void
eval_map_quartic(double t, const double *vc, double *vp, void *ctx)
{
  vp[0] = map_quartic(vc[0]);
}

typedef double (*map_t)(double v);
static const map_t test_maps[] = { map_quad, map_cubic, map_quartic };
static const evalf_t test_eval_maps[] = { eval_map_quad, eval_map_cubic, eval_map_quartic };

// Location of the four 1D cubic coefficients of direction d inside
// vmap_pgkyl, and the scale factor relating them to the 1D coefficients.
static const int pgkyl_idx[3][3][4] = {
  { {0, 1, 2, 3}, {0}, {0} }, // 1V
  { {0, 1, 4, 8}, {0, 2, 5, 9}, {0} }, // 2V
  { {0, 1, 7, 17}, {0, 2, 8, 18}, {0, 3, 9, 19} }, // 3V
};
static const double pgkyl_scale[3] = { 1.0, 1.4142135623730951, 2.0 };

// Nodal values and the nodal derivatives used by the C^1 cubic construction:
// the cubic interpolates the map at the cell vertices and interpolates
// second-order finite-difference derivatives there (centered in the interior,
// one-sided at the domain boundaries); see dg_basis_ops.c. For maps of degree
// <= 2 these finite differences are exact, so the whole cubic is exact.
static void
calc_nodal_vals_and_grads(map_t fn, double lower, double dx, int ncells,
  double *node_f, double *node_g)
{
  for (int j=0; j<ncells+1; ++j)
    node_f[j] = fn(lower + j*dx);

  for (int j=1; j<ncells; ++j)
    node_g[j] = (node_f[j+1]-node_f[j-1])/(2.0*dx);
  node_g[0] = -(node_f[2]-4.0*node_f[1]+3.0*node_f[0])/(2.0*dx);
  node_g[ncells] = (3.0*node_f[ncells]-4.0*node_f[ncells-1]+node_f[ncells-2])/(2.0*dx);
}

// Test the construction of the velocity map object: check that the C^1 cubic
// in each direction matches the map's value and derivative at the cell
// vertices, is C^1 across interior interfaces, and that all the Jacobian
// arrays are the cubic's derivative evaluated at the appropriate
// Gauss-Legendre points.
static void
test_vmap(int vdim, int poly_order, bool use_gpu)
{
  double lower[] = { -2.0, -2.0, -2.0 }, upper[] = { 2.0, 2.0, 2.0 };
  int cells[] = { 8, 6, 4 };
  int nghost[] = { 0, 0, 0 };

  struct gkyl_rect_grid vgrid;
  gkyl_rect_grid_init(&vgrid, vdim, lower, upper, cells);

  struct gkyl_range local_vel, local_ext_vel;
  gkyl_create_grid_ranges(&vgrid, nghost, &local_ext_vel, &local_vel);

  struct gkyl_basis basis_vel;
  gkyl_cart_modal_tensor(&basis_vel, vdim, poly_order);

  struct gkyl_vlasov_velocity_map_inp inp_vmap[GKYL_MAX_CDIM] = { 0 };
  for (int d=0; d<vdim; ++d)
    inp_vmap[d].eval_vmap = test_eval_maps[d];

  struct gkyl_vlasov_velocity_map *vvm = gkyl_vlasov_velocity_map_new(&vgrid,
    &local_vel, &basis_vel, inp_vmap, false, use_gpu);

  TEST_CHECK( !vvm->is_identity );
  TEST_CHECK( vvm->is_mapped );
  TEST_CHECK( vvm->rep == GKYL_VLASOV_VMAP_C1_CUBIC );
  TEST_CHECK( gkyl_vlasov_velocity_map_is_cu_dev(vvm) == use_gpu );

  // Exercise the reference counting.
  struct gkyl_vlasov_velocity_map *vvm2 = gkyl_vlasov_velocity_map_acquire(vvm);
  gkyl_vlasov_velocity_map_release(vvm2);

  // 1D cubic basis for evaluating the per-direction map.
  struct gkyl_basis b1;
  gkyl_cart_modal_tensor(&b1, 1, 3);

  // Expected nodal values and derivatives per direction.
  double node_f[3][16], node_g[3][16];
  for (int d=0; d<vdim; ++d)
    calc_nodal_vals_and_grads(test_maps[d], lower[d], vgrid.dx[d], cells[d],
      node_f[d], node_g[d]);

  // Quadrature range for indexing jacob_vel_gauss, matching the ordering used
  // by the projection and divide/rescale kernels.
  int qshape[GKYL_MAX_DIM];
  for (int d=0; d<vdim; ++d) qshape[d] = poly_order+1;
  struct gkyl_range qrange;
  gkyl_range_init_from_shape(&qrange, vdim, qshape);

  struct gkyl_basis vmap_pgkyl_basis;
  gkyl_cart_modal_serendip(&vmap_pgkyl_basis, vdim, 3);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local_vel);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&local_vel, iter.idx);
    const double *vmap_d = gkyl_array_cfetch(vvm->vmap_host, loc);
    const double *jacob_vel_d = gkyl_array_cfetch(vvm->jacob_vel_host, loc);
    const double *jacob_vel_surf_d = gkyl_array_cfetch(vvm->jacob_vel_surf_host, loc);
    const double *jacob_vel_gauss_d = gkyl_array_cfetch(vvm->jacob_vel_gauss_host, loc);
    const double *vmap_pgkyl_d = gkyl_array_cfetch(vvm->vmap_pgkyl_host, loc);
    const double *vmap_avg_d = gkyl_array_cfetch(vvm->vmap_avg_pgkyl_host, loc);

    for (int d=0; d<vdim; ++d) {
      const double *c_d = &vmap_d[4*d];
      int jl = iter.idx[d]-1, jr = iter.idx[d]; // vertex node indices of this cell
      double dzdv = 2.0/vgrid.dx[d];

      // 1. Vertex values: the cubic interpolates the map at cell vertices.
      double zlo[] = { -1.0 }, zup[] = { 1.0 };
      TEST_CHECK( gkyl_compare_double(b1.eval_expand(zlo, c_d), node_f[d][jl], 1e-12) );
      TEST_MSG("vdim=%d dir=%d cell=%d: lower vertex value %.15e expected %.15e",
        vdim, d, iter.idx[d], b1.eval_expand(zlo, c_d), node_f[d][jl]);
      TEST_CHECK( gkyl_compare_double(b1.eval_expand(zup, c_d), node_f[d][jr], 1e-12) );

      // 2. Vertex derivatives: the cubic interpolates the (finite-difference)
      // nodal derivatives; for the quadratic direction these equal the
      // analytic derivative.
      double grad_lo = dzdv*b1.eval_grad_expand(0, zlo, c_d);
      double grad_up = dzdv*b1.eval_grad_expand(0, zup, c_d);
      TEST_CHECK( gkyl_compare_double(grad_lo, node_g[d][jl], 1e-12) );
      TEST_MSG("vdim=%d dir=%d cell=%d: lower vertex deriv %.15e expected %.15e",
        vdim, d, iter.idx[d], grad_lo, node_g[d][jl]);
      TEST_CHECK( gkyl_compare_double(grad_up, node_g[d][jr], 1e-12) );
      if (d == 0) {
        // Quadratic map: finite differences (and hence the whole cubic) are exact.
        double vlo = lower[d] + jl*vgrid.dx[d], vup = lower[d] + jr*vgrid.dx[d];
        TEST_CHECK( gkyl_compare_double(grad_lo, map_quad_prime(vlo), 1e-12) );
        TEST_CHECK( gkyl_compare_double(grad_up, map_quad_prime(vup), 1e-12) );
      }

      // 3. jacob_vel/jacob_vel_surf: cubic's derivative at the volume and
      // surface Gauss-Legendre points.
      const double *ord_vol = gkyl_gauss_ordinates[poly_order+1];
      for (int i=0; i<poly_order+1; ++i) {
        double z[] = { ord_vol[i] };
        double expected = dzdv*b1.eval_grad_expand(0, z, c_d);
        TEST_CHECK( gkyl_compare_double(jacob_vel_d[d*(poly_order+1)+i], expected, 1e-12) );
        TEST_MSG("vdim=%d dir=%d cell=%d quad=%d: jacob_vel %.15e expected %.15e",
          vdim, d, iter.idx[d], i, jacob_vel_d[d*(poly_order+1)+i], expected);
        if (d == 0) {
          // Quadratic map: also check against the analytic Jacobian at the
          // mapped quadrature point.
          double vc = lower[d] + (jl+0.5)*vgrid.dx[d] + 0.5*vgrid.dx[d]*ord_vol[i];
          TEST_CHECK( gkyl_compare_double(jacob_vel_d[d*(poly_order+1)+i],
            map_quad_prime(vc), 1e-12) );
        }
      }
      const double *ord_surf = gkyl_gauss_ordinates[poly_order+2];
      for (int i=0; i<poly_order+2; ++i) {
        double z[] = { ord_surf[i] };
        double expected = dzdv*b1.eval_grad_expand(0, z, c_d);
        TEST_CHECK( gkyl_compare_double(jacob_vel_surf_d[d*(poly_order+2)+i], expected, 1e-12) );
      }

      // 6. I/O arrays: cell average and coefficient placement in the
      // full-dimensional serendipity representation.
      TEST_CHECK( gkyl_compare_double(vmap_avg_d[d], c_d[0]/sqrt(2.0), 1e-12) );
      const int *pidx = pgkyl_idx[vdim-1][d];
      const double *pgkyl_dir = &vmap_pgkyl_d[d*vmap_pgkyl_basis.num_basis];
      for (int i=0; i<4; ++i) {
        TEST_CHECK( gkyl_compare_double(pgkyl_dir[pidx[i]], pgkyl_scale[vdim-1]*c_d[i], 1e-12) );
        TEST_MSG("vdim=%d dir=%d cell=%d coeff=%d: vmap_pgkyl %.15e expected %.15e",
          vdim, d, iter.idx[d], i, pgkyl_dir[pidx[i]], pgkyl_scale[vdim-1]*c_d[i]);
      }
    }

    // 5. jacob_vel_gauss: product over directions of the per-direction
    // Jacobian, with the last velocity index varying fastest.
    struct gkyl_range_iter qiter;
    gkyl_range_iter_init(&qiter, &qrange);
    while (gkyl_range_iter_next(&qiter)) {
      long qloc = gkyl_range_idx(&qrange, qiter.idx);
      double expected = 1.0;
      for (int d=0; d<vdim; ++d) {
        double z[] = { gkyl_gauss_ordinates[poly_order+1][qiter.idx[d]] };
        expected *= (2.0/vgrid.dx[d])*b1.eval_grad_expand(0, z, &vmap_d[4*d]);
      }
      TEST_CHECK( gkyl_compare_double(jacob_vel_gauss_d[qloc], expected, 1e-12) );
      TEST_MSG("vdim=%d cell=(%d): jacob_vel_gauss[%ld] %.15e expected %.15e",
        vdim, iter.idx[0], qloc, jacob_vel_gauss_d[qloc], expected);
    }
  }

  // 4. C^1 continuity: value and derivative of the per-direction cubic match
  // across every interior interface. Check along the first row of cells in
  // each direction.
  for (int d=0; d<vdim; ++d) {
    for (int i=local_vel.lower[d]; i<local_vel.upper[d]; ++i) {
      int idxl[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];
      for (int k=0; k<vdim; ++k) idxl[k] = idxr[k] = local_vel.lower[k];
      idxl[d] = i; idxr[d] = i+1;

      const double *cl = (const double*) gkyl_array_cfetch(vvm->vmap_host,
        gkyl_range_idx(&local_vel, idxl)) + 4*d;
      const double *cr = (const double*) gkyl_array_cfetch(vvm->vmap_host,
        gkyl_range_idx(&local_vel, idxr)) + 4*d;

      double zup[] = { 1.0 }, zlo[] = { -1.0 };
      TEST_CHECK( gkyl_compare_double(b1.eval_expand(zup, cl), b1.eval_expand(zlo, cr), 1e-12) );
      TEST_CHECK( gkyl_compare_double(b1.eval_grad_expand(0, zup, cl),
        b1.eval_grad_expand(0, zlo, cr), 1e-12) );
      TEST_MSG("vdim=%d dir=%d interface=%d: C^1 continuity", vdim, d, i);
    }
  }

  gkyl_vlasov_velocity_map_release(vvm);
}

// Test the construction of the C^0 linear velocity map used with Serendipity
// velocity bases: the linear function in each direction interpolates the map
// at the cell vertices (the Gauss-Lobatto nodes), the quadratic and cubic
// coefficients of the degenerate-cubic storage vanish identically, every
// Jacobian array holds the cell's constant Jacobian (the ratio of physical
// to computational cell width) at all of its quadrature slots, and the map
// is continuous across interior interfaces.
static void
test_vmap_ser(int vdim, int poly_order, bool use_gpu)
{
  double lower[] = { -2.0, -2.0, -2.0 }, upper[] = { 2.0, 2.0, 2.0 };
  int cells[] = { 8, 6, 4 };
  int nghost[] = { 0, 0, 0 };

  struct gkyl_rect_grid vgrid;
  gkyl_rect_grid_init(&vgrid, vdim, lower, upper, cells);

  struct gkyl_range local_vel, local_ext_vel;
  gkyl_create_grid_ranges(&vgrid, nghost, &local_ext_vel, &local_vel);

  struct gkyl_basis basis_vel;
  gkyl_cart_modal_serendip(&basis_vel, vdim, poly_order);

  struct gkyl_vlasov_velocity_map_inp inp_vmap[GKYL_MAX_CDIM] = { 0 };
  for (int d=0; d<vdim; ++d)
    inp_vmap[d].eval_vmap = test_eval_maps[d];

  struct gkyl_vlasov_velocity_map *vvm = gkyl_vlasov_velocity_map_new(&vgrid,
    &local_vel, &basis_vel, inp_vmap, false, use_gpu);

  TEST_CHECK( !vvm->is_identity );
  TEST_CHECK( vvm->is_mapped );
  TEST_CHECK( vvm->rep == GKYL_VLASOV_VMAP_C0_LINEAR );
  TEST_CHECK( gkyl_vlasov_velocity_map_is_cu_dev(vvm) == use_gpu );

  // 1D cubic basis for evaluating the per-direction (degenerate cubic) map.
  struct gkyl_basis b1;
  gkyl_cart_modal_tensor(&b1, 1, 3);

  struct gkyl_basis vmap_pgkyl_basis;
  gkyl_cart_modal_serendip(&vmap_pgkyl_basis, vdim, 3);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local_vel);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&local_vel, iter.idx);
    const double *vmap_d = gkyl_array_cfetch(vvm->vmap_host, loc);
    const double *jacob_vel_d = gkyl_array_cfetch(vvm->jacob_vel_host, loc);
    const double *jacob_vel_surf_d = gkyl_array_cfetch(vvm->jacob_vel_surf_host, loc);
    const double *jacob_vel_gauss_d = gkyl_array_cfetch(vvm->jacob_vel_gauss_host, loc);
    const double *vmap_pgkyl_d = gkyl_array_cfetch(vvm->vmap_pgkyl_host, loc);
    const double *vmap_avg_d = gkyl_array_cfetch(vvm->vmap_avg_pgkyl_host, loc);

    double jac_tot = 1.0;
    for (int d=0; d<vdim; ++d) {
      const double *c_d = &vmap_d[4*d];
      double vlo = lower[d] + (iter.idx[d]-1)*vgrid.dx[d];
      double vup = lower[d] + iter.idx[d]*vgrid.dx[d];

      // 1. The linear map interpolates the mapping at the cell vertices, and
      // the quadratic and cubic coefficients of the degenerate cubic vanish.
      double zlo[] = { -1.0 }, zup[] = { 1.0 };
      TEST_CHECK( gkyl_compare_double(b1.eval_expand(zlo, c_d), test_maps[d](vlo), 1e-14) );
      TEST_MSG("vdim=%d dir=%d cell=%d: lower vertex value %.15e expected %.15e",
        vdim, d, iter.idx[d], b1.eval_expand(zlo, c_d), test_maps[d](vlo));
      TEST_CHECK( gkyl_compare_double(b1.eval_expand(zup, c_d), test_maps[d](vup), 1e-14) );
      TEST_CHECK( c_d[2] == 0.0 );
      TEST_CHECK( c_d[3] == 0.0 );

      // 2. The constant Jacobian (physical over computational cell width) is
      // replicated at every volume and surface quadrature slot.
      double jac = (test_maps[d](vup)-test_maps[d](vlo))/vgrid.dx[d];
      jac_tot *= jac;
      for (int i=0; i<poly_order+1; ++i) {
        TEST_CHECK( gkyl_compare_double(jacob_vel_d[d*(poly_order+1)+i], jac, 1e-14) );
        TEST_MSG("vdim=%d dir=%d cell=%d quad=%d: jacob_vel %.15e expected %.15e",
          vdim, d, iter.idx[d], i, jacob_vel_d[d*(poly_order+1)+i], jac);
      }
      for (int i=0; i<poly_order+2; ++i)
        TEST_CHECK( gkyl_compare_double(jacob_vel_surf_d[d*(poly_order+2)+i], jac, 1e-14) );

      // 3. I/O arrays: cell average and coefficient placement (zeros in the
      // quadratic and cubic slots).
      TEST_CHECK( gkyl_compare_double(vmap_avg_d[d], c_d[0]/sqrt(2.0), 1e-14) );
      const int *pidx = pgkyl_idx[vdim-1][d];
      const double *pgkyl_dir = &vmap_pgkyl_d[d*vmap_pgkyl_basis.num_basis];
      for (int i=0; i<4; ++i)
        TEST_CHECK( gkyl_compare_double(pgkyl_dir[pidx[i]], pgkyl_scale[vdim-1]*c_d[i], 1e-14) );
    }

    // 4. Total Jacobian: product over directions of the per-direction
    // constants, at every Gauss-Legendre point.
    for (int i=0; i<vvm->jacob_vel_gauss_host->ncomp; ++i)
      TEST_CHECK( gkyl_compare_double(jacob_vel_gauss_d[i], jac_tot, 1e-14) );
  }

  // 5. C^0 continuity: the map's value matches across every interior
  // interface (to roundoff; the nodes are evaluated per cell from the cell
  // center). The Jacobian is deliberately discontinuous there.
  for (int d=0; d<vdim; ++d) {
    for (int i=local_vel.lower[d]; i<local_vel.upper[d]; ++i) {
      int idxl[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];
      for (int k=0; k<vdim; ++k) idxl[k] = idxr[k] = local_vel.lower[k];
      idxl[d] = i; idxr[d] = i+1;

      const double *cl = (const double*) gkyl_array_cfetch(vvm->vmap_host,
        gkyl_range_idx(&local_vel, idxl)) + 4*d;
      const double *cr = (const double*) gkyl_array_cfetch(vvm->vmap_host,
        gkyl_range_idx(&local_vel, idxr)) + 4*d;

      double zup[] = { 1.0 }, zlo[] = { -1.0 };
      TEST_CHECK( gkyl_compare_double(b1.eval_expand(zup, cl), b1.eval_expand(zlo, cr), 1e-14) );
      TEST_MSG("vdim=%d dir=%d interface=%d: C^0 continuity %.15e vs %.15e",
        vdim, d, i, b1.eval_expand(zup, cl), b1.eval_expand(zlo, cr));
    }
  }

  gkyl_vlasov_velocity_map_release(vvm);
}

// Identity (uniform grid) default: all input mappings NULL. Valid for both
// representations: the C^1 cubic and the C^0 linear reproduce the identity
// map exactly, with unit Jacobians everywhere.
static void
test_vmap_identity(int vdim, int poly_order, bool use_gpu, bool ser_basis)
{
  double lower[] = { -2.0, -2.0, -2.0 }, upper[] = { 2.0, 2.0, 2.0 };
  int cells[] = { 8, 6, 4 };
  int nghost[] = { 0, 0, 0 };

  struct gkyl_rect_grid vgrid;
  gkyl_rect_grid_init(&vgrid, vdim, lower, upper, cells);

  struct gkyl_range local_vel, local_ext_vel;
  gkyl_create_grid_ranges(&vgrid, nghost, &local_ext_vel, &local_vel);

  struct gkyl_basis basis_vel;
  if (ser_basis)
    gkyl_cart_modal_serendip(&basis_vel, vdim, poly_order);
  else
    gkyl_cart_modal_tensor(&basis_vel, vdim, poly_order);

  struct gkyl_vlasov_velocity_map_inp inp_vmap[GKYL_MAX_CDIM] = { 0 };

  struct gkyl_vlasov_velocity_map *vvm = gkyl_vlasov_velocity_map_new(&vgrid,
    &local_vel, &basis_vel, inp_vmap, false, use_gpu);

  TEST_CHECK( vvm->is_identity );
  TEST_CHECK( vvm->is_mapped ); // always true: single code path.
  TEST_CHECK( vvm->rep == (ser_basis || poly_order == 1 ?
    GKYL_VLASOV_VMAP_C0_LINEAR : GKYL_VLASOV_VMAP_C1_CUBIC) );

  struct gkyl_basis b1;
  gkyl_cart_modal_tensor(&b1, 1, 3);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local_vel);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&local_vel, iter.idx);
    const double *vmap_d = gkyl_array_cfetch(vvm->vmap_host, loc);
    const double *jacob_vel_d = gkyl_array_cfetch(vvm->jacob_vel_host, loc);
    const double *jacob_vel_surf_d = gkyl_array_cfetch(vvm->jacob_vel_surf_host, loc);
    const double *jacob_vel_gauss_d = gkyl_array_cfetch(vvm->jacob_vel_gauss_host, loc);

    for (int d=0; d<vdim; ++d) {
      // Identity map: vmap reproduces v at the vertices, all Jacobians are 1.
      double zlo[] = { -1.0 }, zup[] = { 1.0 };
      double vlo = lower[d] + (iter.idx[d]-1)*vgrid.dx[d];
      double vup = lower[d] + iter.idx[d]*vgrid.dx[d];
      TEST_CHECK( gkyl_compare_double(b1.eval_expand(zlo, &vmap_d[4*d]), vlo, 1e-13) );
      TEST_CHECK( gkyl_compare_double(b1.eval_expand(zup, &vmap_d[4*d]), vup, 1e-13) );
      for (int i=0; i<poly_order+1; ++i)
        TEST_CHECK( gkyl_compare_double(jacob_vel_d[d*(poly_order+1)+i], 1.0, 1e-13) );
      for (int i=0; i<poly_order+2; ++i)
        TEST_CHECK( gkyl_compare_double(jacob_vel_surf_d[d*(poly_order+2)+i], 1.0, 1e-13) );
    }
    for (int i=0; i<vvm->jacob_vel_gauss_host->ncomp; ++i)
      TEST_CHECK( gkyl_compare_double(jacob_vel_gauss_d[i], 1.0, 1e-13) );
  }

  gkyl_vlasov_velocity_map_release(vvm);
}

static void test_vmap_1v_p2(void) { test_vmap(1, 2, false); }
static void test_vmap_1v_p3(void) { test_vmap(1, 3, false); }
static void test_vmap_2v_p2(void) { test_vmap(2, 2, false); }
static void test_vmap_3v_p2(void) { test_vmap(3, 2, false); }
static void test_vmap_identity_1v_p2(void) { test_vmap_identity(1, 2, false, false); }
static void test_vmap_identity_2v_p2(void) { test_vmap_identity(2, 2, false, false); }
static void test_vmap_identity_3v_p2(void) { test_vmap_identity(3, 2, false, false); }

static void test_vmap_ser_1v_p1(void) { test_vmap_ser(1, 1, false); }
static void test_vmap_ser_1v_p2(void) { test_vmap_ser(1, 2, false); }
static void test_vmap_ser_1v_p3(void) { test_vmap_ser(1, 3, false); }
static void test_vmap_ser_2v_p1(void) { test_vmap_ser(2, 1, false); }
static void test_vmap_ser_2v_p2(void) { test_vmap_ser(2, 2, false); }
static void test_vmap_ser_3v_p1(void) { test_vmap_ser(3, 1, false); }
static void test_vmap_ser_3v_p2(void) { test_vmap_ser(3, 2, false); }
static void test_vmap_ser_identity_1v_p1(void) { test_vmap_identity(1, 1, false, true); }
static void test_vmap_ser_identity_2v_p2(void) { test_vmap_identity(2, 2, false, true); }
static void test_vmap_ser_identity_3v_p1(void) { test_vmap_identity(3, 1, false, true); }

#ifdef GKYL_HAVE_CUDA
// On GPUs the solver arrays are device-resident: copy them back to the host
// and compare against a host-built map. The map itself is constructed on the
// CPU and copied over, so this checks the device copy of every solver array
// for both representations, for both user mappings and the identity
// (uniform grid) default.
static void
test_vmap_gpu(int vdim, int poly_order, bool identity, bool ser_basis)
{
  double lower[] = { -2.0, -2.0, -2.0 }, upper[] = { 2.0, 2.0, 2.0 };
  int cells[] = { 8, 6, 4 };
  int nghost[] = { 0, 0, 0 };

  struct gkyl_rect_grid vgrid;
  gkyl_rect_grid_init(&vgrid, vdim, lower, upper, cells);

  struct gkyl_range local_vel, local_ext_vel;
  gkyl_create_grid_ranges(&vgrid, nghost, &local_ext_vel, &local_vel);

  struct gkyl_basis basis_vel;
  if (ser_basis)
    gkyl_cart_modal_serendip(&basis_vel, vdim, poly_order);
  else
    gkyl_cart_modal_tensor(&basis_vel, vdim, poly_order);

  struct gkyl_vlasov_velocity_map_inp inp_vmap[GKYL_MAX_CDIM] = { 0 };
  if (!identity)
    for (int d=0; d<vdim; ++d)
      inp_vmap[d].eval_vmap = test_eval_maps[d];

  struct gkyl_vlasov_velocity_map *vvm_ho = gkyl_vlasov_velocity_map_new(&vgrid,
    &local_vel, &basis_vel, inp_vmap, false, false);
  struct gkyl_vlasov_velocity_map *vvm = gkyl_vlasov_velocity_map_new(&vgrid,
    &local_vel, &basis_vel, inp_vmap, false, true);

  TEST_CHECK( gkyl_vlasov_velocity_map_is_cu_dev(vvm) );
  TEST_CHECK( vvm->is_identity == identity );

  struct gkyl_array *dev_arrays[] = { vvm->vmap, vvm->jacob_vel, vvm->jacob_vel_surf, vvm->jacob_vel_gauss };
  struct gkyl_array *ho_arrays[] = { vvm_ho->vmap, vvm_ho->jacob_vel, vvm_ho->jacob_vel_surf, vvm_ho->jacob_vel_gauss };
  for (int n=0; n<4; ++n) {
    struct gkyl_array *buf = gkyl_array_new(GKYL_DOUBLE, dev_arrays[n]->ncomp, dev_arrays[n]->size);
    gkyl_array_copy(buf, dev_arrays[n]);
    for (long i=0; i<buf->size; ++i) {
      const double *a = gkyl_array_cfetch(buf, i);
      const double *b = gkyl_array_cfetch(ho_arrays[n], i);
      for (int c=0; c<buf->ncomp; ++c)
        TEST_CHECK( gkyl_compare_double(a[c], b[c], 1e-15) );
    }
    gkyl_array_release(buf);
  }

  gkyl_vlasov_velocity_map_release(vvm_ho);
  gkyl_vlasov_velocity_map_release(vvm);
}

static void test_vmap_1v_p2_gpu(void) { test_vmap_gpu(1, 2, false, false); }
static void test_vmap_1v_p3_gpu(void) { test_vmap_gpu(1, 3, false, false); }
static void test_vmap_2v_p2_gpu(void) { test_vmap_gpu(2, 2, false, false); }
static void test_vmap_3v_p2_gpu(void) { test_vmap_gpu(3, 2, false, false); }
static void test_vmap_identity_2v_p2_gpu(void) { test_vmap_gpu(2, 2, true, false); }
static void test_vmap_ser_1v_p1_gpu(void) { test_vmap_gpu(1, 1, false, true); }
static void test_vmap_ser_2v_p2_gpu(void) { test_vmap_gpu(2, 2, false, true); }
static void test_vmap_ser_3v_p1_gpu(void) { test_vmap_gpu(3, 1, false, true); }
static void test_vmap_ser_identity_2v_p1_gpu(void) { test_vmap_gpu(2, 1, true, true); }
#endif

TEST_LIST = {
  { "vmap_1v_p2", test_vmap_1v_p2 },
  { "vmap_1v_p3", test_vmap_1v_p3 },
  { "vmap_2v_p2", test_vmap_2v_p2 },
  { "vmap_3v_p2", test_vmap_3v_p2 },
  { "vmap_identity_1v_p2", test_vmap_identity_1v_p2 },
  { "vmap_identity_2v_p2", test_vmap_identity_2v_p2 },
  { "vmap_identity_3v_p2", test_vmap_identity_3v_p2 },
  { "vmap_ser_1v_p1", test_vmap_ser_1v_p1 },
  { "vmap_ser_1v_p2", test_vmap_ser_1v_p2 },
  { "vmap_ser_1v_p3", test_vmap_ser_1v_p3 },
  { "vmap_ser_2v_p1", test_vmap_ser_2v_p1 },
  { "vmap_ser_2v_p2", test_vmap_ser_2v_p2 },
  { "vmap_ser_3v_p1", test_vmap_ser_3v_p1 },
  { "vmap_ser_3v_p2", test_vmap_ser_3v_p2 },
  { "vmap_ser_identity_1v_p1", test_vmap_ser_identity_1v_p1 },
  { "vmap_ser_identity_2v_p2", test_vmap_ser_identity_2v_p2 },
  { "vmap_ser_identity_3v_p1", test_vmap_ser_identity_3v_p1 },
#ifdef GKYL_HAVE_CUDA
  { "vmap_1v_p2_gpu", test_vmap_1v_p2_gpu },
  { "vmap_1v_p3_gpu", test_vmap_1v_p3_gpu },
  { "vmap_2v_p2_gpu", test_vmap_2v_p2_gpu },
  { "vmap_3v_p2_gpu", test_vmap_3v_p2_gpu },
  { "vmap_identity_2v_p2_gpu", test_vmap_identity_2v_p2_gpu },
  { "vmap_ser_1v_p1_gpu", test_vmap_ser_1v_p1_gpu },
  { "vmap_ser_2v_p2_gpu", test_vmap_ser_2v_p2_gpu },
  { "vmap_ser_3v_p1_gpu", test_vmap_ser_3v_p1_gpu },
  { "vmap_ser_identity_2v_p1_gpu", test_vmap_ser_identity_2v_p1_gpu },
#endif
  { NULL, NULL },
};
