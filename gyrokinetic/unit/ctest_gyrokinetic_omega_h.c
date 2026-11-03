#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_gyrokinetic_omega_h.h>

/**
 * Evaluate analytic density and geometry fields representable in the tested basis.
 *
 * @param kind Manufactured test case.
 * @param ndim Configuration dimension.
 * @param poly_order Polynomial order.
 * @param cell_tag First (1), last (2), or other (0) radial cell.
 * @param xc Coordinates in the reference cell.
 * @param density Output Jacobian-weighted density.
 * @param geometry Output geometry factor.
 */
static void
sample_fields(
  int kind, int ndim, int poly_order, int cell_tag, const double *xc, double *density,
  double *geometry
)
{
  *density = 9.0;
  *geometry = 2.0;
  if (kind == 1) {
    *density = cell_tag == 1 ? 100.0 : 1.0;
    *geometry = cell_tag == 1 ? 1.0 : 10.0;
  } else if (kind == 2) {
    *density = 4.0;
    *geometry = 2.0;
    for (int dir = 0; dir < ndim; ++dir) {
      *density += 0.2 * (dir + 1) * xc[dir];
      *geometry += 0.1 * (dir + 1) * xc[dir];
      if (poly_order == 2) {
        *density -= 0.1 * xc[dir] * xc[dir];
        *geometry -= 0.2 * xc[dir] * xc[dir];
      }
    }
  } else if (kind == 3) {
    // For p1 the product peaks at the center, above both Gauss-node values.
    *density = 1.0 + xc[0];
    *geometry = 1.0 - 0.5 * xc[0];
  } else if (kind == 5) {
    *density = cell_tag == 2 ? 121.0 : 1.0;
    *geometry = 1.0;
  } else if (kind == 6) {
    *density = 9.0;
    *geometry = -2.0;
  } else if (kind == 4) {
    *density = -2.0;
    *geometry = 3.0;
  }
}

/**
 * Check analytic products, ghost exclusion, partitioning and repeated reductions.
 *
 * @param ndim Configuration dimension.
 * @param poly_order Polynomial order.
 * @param tensor Whether to use a tensor basis instead of serendipity.
 * @param use_gpu Whether to exercise the CUDA path.
 */
static void
check_estimator(int ndim, int poly_order, bool tensor, bool use_gpu)
{
  struct gkyl_basis basis;
  if (tensor) {
    gkyl_cart_modal_tensor(&basis, ndim, poly_order);
  } else {
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  }
  int lower[3] = {1, 1, 1}, upper[3] = {515, 1, 1};
  int ext_lower[3] = {0, 0, 0}, ext_upper[3] = {516, 2, 2};
  if (ndim > 1) {
    upper[0] = 9;
    ext_upper[0] = 10;
    upper[1] = 3;
    ext_upper[1] = 4;
    if (ndim == 3) {
      upper[2] = 2;
      ext_upper[2] = 3;
    }
  }
  struct gkyl_range local_ext, local;
  gkyl_range_init(&local_ext, ndim, ext_lower, ext_upper);
  gkyl_sub_range_init(&local, &local_ext, lower, upper);
  struct gkyl_array *density = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  struct gkyl_array *geometry = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  struct gkyl_array *density_dev =
    use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume) :
              gkyl_array_acquire(density);
  struct gkyl_array *geometry_dev =
    use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume) :
              gkyl_array_acquire(geometry);
  double *nodes = gkyl_malloc(ndim * basis.num_basis * sizeof(double));
  double *density_nodal = gkyl_malloc(basis.num_basis * sizeof(double));
  double *geometry_nodal = gkyl_malloc(basis.num_basis * sizeof(double));
  basis.node_list(nodes);

  for (int kind = 0; kind < 7; ++kind) {
    // Large ghost coefficients must never enter the interior maximum.
    gkyl_array_clear(density, 1e8);
    gkyl_array_clear(geometry, 1e8);
    double expected = 0.0;
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &local);
    while (gkyl_range_iter_next(&iter)) {
      int cell_tag = iter.idx[0] == lower[0] ? 1 : (iter.idx[0] == upper[0] ? 2 : 0);
      long loc = gkyl_range_idx(&local, iter.idx);
      for (int node = 0; node < basis.num_basis; ++node) {
        sample_fields(
          kind, ndim, poly_order, cell_tag, &nodes[node * ndim], &density_nodal[node],
          &geometry_nodal[node]
        );
      }
      basis.nodal_to_modal(density_nodal, gkyl_array_fetch(density, loc));
      basis.nodal_to_modal(geometry_nodal, gkyl_array_fetch(geometry, loc));
      int num_quad = 1;
      for (int dir = 0; dir < ndim; ++dir) {
        num_quad *= poly_order + 1;
      }
      for (int node = 0; node <= num_quad; ++node) {
        double xc[3] = {0.0}, density_ref, geometry_ref;
        int ordinal = node;
        if (node < num_quad) {
          for (int dir = 0; dir < ndim; ++dir) {
            xc[dir] = gkyl_gauss_ordinates[poly_order + 1][ordinal % (poly_order + 1)];
            ordinal /= poly_order + 1;
          }
        }
        sample_fields(kind, ndim, poly_order, cell_tag, xc, &density_ref, &geometry_ref);
        expected = fmax(expected, fmax(0.0, density_ref) * geometry_ref * geometry_ref);
      }
    }
    gkyl_array_copy(density_dev, density);
    gkyl_array_copy(geometry_dev, geometry);
    struct gkyl_gyrokinetic_omega_h *up =
      gkyl_gyrokinetic_omega_h_new(&basis, &local, geometry_dev, use_gpu);
    double result = gkyl_gyrokinetic_omega_h_advance(up, density_dev);
    TEST_CHECK(gkyl_compare_double(result, expected, 2e-13));
    TEST_MSG(
      "dim=%d p=%d kind=%d actual=%.17g expected=%.17g", ndim, poly_order, kind, result, expected
    );
    if (kind == 0) {
      TEST_CHECK(gkyl_compare_double(result, 36.0, 2e-13));
    }
    if (kind == 1) {
      TEST_CHECK(gkyl_compare_double(result, 100.0, 2e-13));
    }
    if (kind == 3) {
      TEST_CHECK(gkyl_compare_double(result, 1.0, 2e-13));
    }

    // Splitting the same parent range must preserve the global maximum.
    double split_max = 0.0;
    for (int part = 0; part < 2; ++part) {
      int split_lo[3], split_up[3];
      for (int dir = 0; dir < ndim; ++dir) {
        split_lo[dir] = lower[dir];
        split_up[dir] = upper[dir];
      }
      if (part == 0) {
        split_up[0] = upper[0] / 2;
      } else {
        split_lo[0] = upper[0] / 2 + 1;
      }
      struct gkyl_range split;
      gkyl_sub_range_init(&split, &local_ext, split_lo, split_up);
      struct gkyl_gyrokinetic_omega_h *part_up =
        gkyl_gyrokinetic_omega_h_new(&basis, &split, geometry_dev, use_gpu);
      split_max = fmax(split_max, gkyl_gyrokinetic_omega_h_advance(part_up, density_dev));
      gkyl_gyrokinetic_omega_h_release(part_up);
    }
    TEST_CHECK(gkyl_compare_double(split_max, result, 2e-13));

    // A smaller subsequent result must not retain a stale GPU atomic maximum.
    gkyl_array_scale(density_dev, 0.25);
    TEST_CHECK(
      gkyl_compare_double(gkyl_gyrokinetic_omega_h_advance(up, density_dev), 0.25 * expected, 2e-13)
    );
    gkyl_array_clear(density_dev, 0.0);
    TEST_CHECK(gkyl_gyrokinetic_omega_h_advance(up, density_dev) == 0.0);
    gkyl_gyrokinetic_omega_h_release(up);
  }
  int empty_upper[3] = {0, 1, 1};
  struct gkyl_range empty;
  gkyl_sub_range_init(&empty, &local_ext, lower, empty_upper);
  struct gkyl_gyrokinetic_omega_h *empty_up =
    gkyl_gyrokinetic_omega_h_new(&basis, &empty, geometry_dev, use_gpu);
  TEST_CHECK(gkyl_gyrokinetic_omega_h_advance(empty_up, density_dev) == 0.0);
  gkyl_gyrokinetic_omega_h_release(empty_up);
  gkyl_gyrokinetic_omega_h_release(NULL);
  gkyl_free(nodes);
  gkyl_free(density_nodal);
  gkyl_free(geometry_nodal);
  gkyl_array_release(density);
  gkyl_array_release(geometry);
  gkyl_array_release(density_dev);
  gkyl_array_release(geometry_dev);
}

/**
 * Exercise configuration dimensions 1--3 and orders 1--2 on CPU.
 *
 */
static void
test_cpu(void)
{
  for (int ndim = 1; ndim <= 3; ++ndim) {
    for (int order = 1; order <= 2; ++order) {
      for (int tensor = 0; tensor <= 1; ++tensor) {
        check_estimator(ndim, order, tensor, false);
      }
    }
  }
}

#ifdef GKYL_HAVE_CUDA
/**
 * Exercise the same analytic cases and CUDA block reduction sizes on the GPU.
 *
 */
static void
test_gpu(void)
{
  for (int ndim = 1; ndim <= 3; ++ndim) {
    for (int order = 1; order <= 2; ++order) {
      for (int tensor = 0; tensor <= 1; ++tensor) {
        check_estimator(ndim, order, tensor, true);
      }
    }
  }
}
#endif

TEST_LIST = {
  {"omega_h_cpu", test_cpu},
#ifdef GKYL_HAVE_CUDA
  {"omega_h_gpu", test_gpu},
#endif
  {NULL, NULL}
};
