#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_dg_reduce.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_util.h>
#include <time.h>

void
test_reduce_dg(bool use_gpu)
{
  int poly_order = 1;
  int ncomp = 3;
  double lower[] = {-M_PI}, upper[] = {M_PI};
  int cells[] = {20};

  int ndim = sizeof(lower) / sizeof(lower[0]);

  struct gkyl_basis *basis;
  if (use_gpu) {
    basis = gkyl_cu_malloc(sizeof(struct gkyl_basis));
    gkyl_cart_modal_serendip_cu_dev(basis, ndim, poly_order);
  } else {
    basis = gkyl_malloc(sizeof(struct gkyl_basis));
    gkyl_cart_modal_serendip(basis, ndim, poly_order);
  }
  struct gkyl_basis basis_ho;
  gkyl_cart_modal_serendip(&basis_ho, ndim, poly_order);

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  int ghost[ndim];
  for (int d = 0; d < ndim; d++) {
    ghost[d] = 1;
  }
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct gkyl_array *arr =
    use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, ncomp * basis_ho.num_basis, local_ext.volume) :
              gkyl_array_new(GKYL_DOUBLE, ncomp * basis_ho.num_basis, local_ext.volume);
  struct gkyl_array *arr_ho = use_gpu ? gkyl_array_new(GKYL_DOUBLE, arr->ncomp, arr->size) :
                                        gkyl_array_acquire(arr);

  // Load 1D Gauss-Legendre nodes.
  int num_quad = poly_order + 1;
  double ordinates1[num_quad];
  memcpy(ordinates1, gkyl_gauss_ordinates[num_quad], sizeof(double[num_quad]));

  // Create range to loop over nodes.
  int qshape[GKYL_MAX_DIM];
  for (int i = 0; i < ndim; ++i) {
    qshape[i] = num_quad;
  }
  struct gkyl_range qrange;
  gkyl_range_init_from_shape(&qrange, ndim, qshape);

  // Create nodes array.
  int num_nodes = qrange.volume;
  struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, ndim, num_nodes);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &qrange);
  while (gkyl_range_iter_next(&iter)) {
    long linc = gkyl_range_idx(&qrange, iter.idx);
    double *nod = gkyl_array_fetch(nodes, linc);
    for (int i = 0; i < ndim; ++i) {
      nod[i] = ordinates1[iter.idx[i] - qrange.lower[i]];
    }
  }

  // Populate arr with a function evaluated at nodes (transformed to modal).
  // Keep track of nodal max, min and sum.
  double arr_max[ncomp], arr_min[ncomp], arr_sum[ncomp];
  for (size_t i = 0; i < ncomp; ++i) {
    arr_max[i] = -1e20;
    arr_min[i] = 1e20;
    arr_sum[i] = 0.0;
  }

  for (size_t i = 0; i < arr->size; ++i) {
    double *arr_c = gkyl_array_fetch(arr_ho, i);

    int idx[GKYL_MAX_DIM];
    gkyl_range_inv_idx(&local_ext, i, idx);

    double xc[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(&grid, idx, xc);

    for (int ci = 0; ci < ncomp; ci++) {
      double arr_nodal[num_nodes];
      for (size_t k = 0; k < num_nodes; k++) {
        const double *nod = gkyl_array_cfetch(nodes, k);
        double x[GKYL_MAX_DIM];
        for (int d = 0; d < ndim; d++) {
          x[d] = xc[d] + 0.5 * grid.dx[0] * nod[d];
        }

        arr_nodal[k] = (ci + 1);
        for (int d = 0; d < ndim; d++) {
          arr_nodal[k] *= sin(((d + 1) * 2.0 * M_PI / (upper[d] - lower[d])) * x[d]);
        }

        arr_max[ci] = GKYL_MAX2(arr_max[ci], arr_nodal[k]);
        arr_min[ci] = GKYL_MIN2(arr_min[ci], arr_nodal[k]);
        arr_sum[ci] += arr_nodal[k];
      }

      for (size_t k = 0; k < basis_ho.num_basis; k++) {
        basis_ho.quad_nodal_to_modal(arr_nodal, &arr_c[ci * basis_ho.num_basis], k);
      }
    }
  }
  gkyl_array_copy(arr, arr_ho);

  // Reduce each vector component.
  double *amax, *amin, *asum;
  if (use_gpu) {
    amax = gkyl_cu_malloc(ncomp * sizeof(double));
    amin = gkyl_cu_malloc(ncomp * sizeof(double));
    asum = gkyl_cu_malloc(ncomp * sizeof(double));
  } else {
    amax = gkyl_malloc(ncomp * sizeof(double));
    amin = gkyl_malloc(ncomp * sizeof(double));
    asum = gkyl_malloc(ncomp * sizeof(double));
  }

  for (int ci = 0; ci < ncomp; ci++) {
    gkyl_array_dg_reducec(&amax[ci], arr, ci, GKYL_MAX, basis);
    gkyl_array_dg_reducec(&amin[ci], arr, ci, GKYL_MIN, basis);
    gkyl_array_dg_reducec(&asum[ci], arr, ci, GKYL_SUM, basis);
  }

  // Check results.
  double amax_ho[ncomp], amin_ho[ncomp], asum_ho[ncomp];
  if (use_gpu) {
    gkyl_cu_memcpy(amax_ho, amax, ncomp * sizeof(double), GKYL_CU_MEMCPY_D2H);
    gkyl_cu_memcpy(amin_ho, amin, ncomp * sizeof(double), GKYL_CU_MEMCPY_D2H);
    gkyl_cu_memcpy(asum_ho, asum, ncomp * sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    memcpy(amax_ho, amax, ncomp * sizeof(double));
    memcpy(amin_ho, amin, ncomp * sizeof(double));
    memcpy(asum_ho, asum, ncomp * sizeof(double));
  }

  for (int c = 0; c < ncomp; ++c) {
    TEST_CHECK(gkyl_compare(arr_max[c], amax_ho[c], 1e-14));
    TEST_MSG("%d MAX: Expected: %g | Got: %g\n", c, arr_max[c], amax_ho[c]);
    TEST_CHECK(gkyl_compare(arr_min[c], amin_ho[c], 1e-14));
    TEST_MSG("%d MIN: Expected: %g | Got: %g\n", c, arr_min[c], amin_ho[c]);
    TEST_CHECK(gkyl_compare(arr_sum[c], asum_ho[c], 1e-14));
    TEST_MSG("%d SUM: Expected: %g | Got: %g\n", c, arr_sum[c], asum_ho[c]);
  }

  gkyl_array_release(arr_ho);
  gkyl_array_release(arr);
  gkyl_array_release(nodes);
  if (use_gpu) {
    gkyl_cart_modal_basis_release_cu(basis);
    gkyl_cu_free(amax);
    gkyl_cu_free(amin);
    gkyl_cu_free(asum);
  } else {
    gkyl_cart_modal_basis_release(basis);
    gkyl_free(amax);
    gkyl_free(amin);
    gkyl_free(asum);
  }
}

void
test_reduce_dg_range(bool use_gpu)
{
  int poly_order = 1;
  int ncomp = 3;
  double lower[] = {-M_PI}, upper[] = {M_PI};
  int cells[] = {20};

  int ndim = sizeof(lower) / sizeof(lower[0]);

  struct gkyl_basis *basis;
  if (use_gpu) {
    basis = gkyl_cu_malloc(sizeof(struct gkyl_basis));
    gkyl_cart_modal_serendip_cu_dev(basis, ndim, poly_order);
  } else {
    basis = gkyl_malloc(sizeof(struct gkyl_basis));
    gkyl_cart_modal_serendip(basis, ndim, poly_order);
  }
  struct gkyl_basis basis_ho;
  gkyl_cart_modal_serendip(&basis_ho, ndim, poly_order);

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  int ghost[ndim];
  for (int d = 0; d < ndim; d++) {
    ghost[d] = 1;
  }
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct gkyl_array *arr =
    use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, ncomp * basis_ho.num_basis, local_ext.volume) :
              gkyl_array_new(GKYL_DOUBLE, ncomp * basis_ho.num_basis, local_ext.volume);
  struct gkyl_array *arr_ho = use_gpu ? gkyl_array_new(GKYL_DOUBLE, arr->ncomp, arr->size) :
                                        gkyl_array_acquire(arr);

  // Load 1D Gauss-Legendre nodes.
  int num_quad = poly_order + 1;
  double ordinates1[num_quad];
  memcpy(ordinates1, gkyl_gauss_ordinates[num_quad], sizeof(double[num_quad]));

  // Create range to loop over nodes.
  int qshape[GKYL_MAX_DIM];
  for (int i = 0; i < ndim; ++i) {
    qshape[i] = num_quad;
  }
  struct gkyl_range qrange;
  gkyl_range_init_from_shape(&qrange, ndim, qshape);

  // Create nodes array.
  int num_nodes = qrange.volume;
  struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, ndim, num_nodes);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &qrange);
  while (gkyl_range_iter_next(&iter)) {
    long linc = gkyl_range_idx(&qrange, iter.idx);
    double *nod = gkyl_array_fetch(nodes, linc);
    for (int i = 0; i < ndim; ++i) {
      nod[i] = ordinates1[iter.idx[i] - qrange.lower[i]];
    }
  }

  // Populate arr with a function evaluated at nodes (transformed to modal).
  double arr_max[ncomp], arr_min[ncomp], arr_sum[ncomp];
  for (size_t i = 0; i < ncomp; ++i) {
    arr_max[i] = -1e20;
    arr_min[i] = 1e20;
    arr_sum[i] = 0.0;
  }

  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&local, iter.idx);

    double *arr_c = gkyl_array_fetch(arr_ho, linidx);

    double xc[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(&grid, iter.idx, xc);

    for (int ci = 0; ci < ncomp; ci++) {
      double arr_nodal[num_nodes];
      for (size_t k = 0; k < num_nodes; k++) {
        const double *nod = gkyl_array_cfetch(nodes, k);
        double x[GKYL_MAX_DIM];
        for (int d = 0; d < ndim; d++) {
          x[d] = xc[d] + 0.5 * grid.dx[0] * nod[d];
        }

        arr_nodal[k] = (ci + 1);
        for (int d = 0; d < ndim; d++) {
          arr_nodal[k] *= sin(((d + 1) * 2.0 * M_PI / (upper[d] - lower[d])) * x[d]);
        }

        arr_max[ci] = GKYL_MAX2(arr_max[ci], arr_nodal[k]);
        arr_min[ci] = GKYL_MIN2(arr_min[ci], arr_nodal[k]);
        arr_sum[ci] += arr_nodal[k];
      }

      for (size_t k = 0; k < basis_ho.num_basis; k++) {
        basis_ho.quad_nodal_to_modal(arr_nodal, &arr_c[ci * basis_ho.num_basis], k);
      }
    }
  }
  gkyl_array_copy(arr, arr_ho);

  // Reduce each vector component.
  double *amax, *amin, *asum;
  if (use_gpu) {
    amax = gkyl_cu_malloc(ncomp * sizeof(double));
    amin = gkyl_cu_malloc(ncomp * sizeof(double));
    asum = gkyl_cu_malloc(ncomp * sizeof(double));
  } else {
    amax = gkyl_malloc(ncomp * sizeof(double));
    amin = gkyl_malloc(ncomp * sizeof(double));
    asum = gkyl_malloc(ncomp * sizeof(double));
  }
  for (int ci = 0; ci < ncomp; ci++) {
    gkyl_array_dg_reducec_range(&amax[ci], arr, ci, GKYL_MAX, basis, &local);
    gkyl_array_dg_reducec_range(&amin[ci], arr, ci, GKYL_MIN, basis, &local);
    gkyl_array_dg_reducec_range(&asum[ci], arr, ci, GKYL_SUM, basis, &local);
  }

  // Check results.
  double amax_ho[ncomp], amin_ho[ncomp], asum_ho[ncomp];
  if (use_gpu) {
    gkyl_cu_memcpy(amax_ho, amax, ncomp * sizeof(double), GKYL_CU_MEMCPY_D2H);
    gkyl_cu_memcpy(amin_ho, amin, ncomp * sizeof(double), GKYL_CU_MEMCPY_D2H);
    gkyl_cu_memcpy(asum_ho, asum, ncomp * sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    memcpy(amax_ho, amax, ncomp * sizeof(double));
    memcpy(amin_ho, amin, ncomp * sizeof(double));
    memcpy(asum_ho, asum, ncomp * sizeof(double));
  }

  for (int c = 0; c < ncomp; ++c) {
    TEST_CHECK(gkyl_compare(arr_max[c], amax_ho[c], 1e-14));
    TEST_MSG("%d MAX: Expected: %g | Got: %g\n", c, arr_max[c], amax_ho[c]);
    TEST_CHECK(gkyl_compare(arr_min[c], amin_ho[c], 1e-14));
    TEST_MSG("%d MIN: Expected: %g | Got: %g\n", c, arr_min[c], amin_ho[c]);
    TEST_CHECK(gkyl_compare(arr_sum[c], asum_ho[c], 1e-14));
    TEST_MSG("%d SUM: Expected: %g | Got: %g\n", c, arr_sum[c], asum_ho[c]);
  }

  gkyl_array_release(arr_ho);
  gkyl_array_release(arr);
  gkyl_array_release(nodes);
  if (use_gpu) {
    gkyl_cart_modal_basis_release_cu(basis);
    gkyl_cu_free(amax);
    gkyl_cu_free(amin);
    gkyl_cu_free(asum);
  } else {
    gkyl_cart_modal_basis_release(basis);
    gkyl_free(amax);
    gkyl_free(amin);
    gkyl_free(asum);
  }
}

void
test_reduce_dg_ho()
{
  test_reduce_dg(false);
}

void
test_reduce_dg_range_ho()
{
  test_reduce_dg_range(false);
}

static void
test_reduce_dg_edges(bool use_gpu, bool pinned)
{
  // Check DG reductions across multiple blocks, partial blocks and parent-strided
  // subranges, with ghosts, repeated output reuse and empty selections.
  struct gkyl_basis basis_ho;
  gkyl_cart_modal_serendip(&basis_ho, 2, 1);

  struct gkyl_basis *basis = &basis_ho;
  if (use_gpu) {
    basis = gkyl_cu_malloc(sizeof(*basis));
    gkyl_cart_modal_serendip_cu_dev(basis, 2, 1);
  }

  int upper = 2 * GKYL_DEFAULT_NUM_THREADS + 2;

  struct gkyl_range parent, local, empty;
  gkyl_range_init(&parent, 2, (int[]){-2, 4}, (int[]){upper + 2, 8});
  gkyl_sub_range_init(&local, &parent, (int[]){0, 5}, (int[]){upper, 7});
  gkyl_sub_range_init(&empty, &parent, (int[]){0, 5}, (int[]){-1, 7});

  struct gkyl_array *arr_ho = gkyl_array_new(GKYL_DOUBLE, 2 * basis_ho.num_basis, parent.volume);
  struct gkyl_array *arr = use_gpu ?
                             gkyl_array_cu_dev_new(GKYL_DOUBLE, arr_ho->ncomp, arr_ho->size) :
                             gkyl_array_acquire(arr_ho);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &parent);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&parent, iter.idx);
    double *coeff = gkyl_array_fetch(arr_ho, loc);
    bool interior = gkyl_range_contains_idx(&local, iter.idx);
    for (int comp = 0; comp < 2; ++comp) {
      double base = interior ? 1 + 0.01 * loc : (loc % 2 ? -1e4 : 1e4);
      coeff[4 * comp] = 2 * (comp + 1) * base;
      coeff[4 * comp + 1] = 0.3 * (comp + 1);
      coeff[4 * comp + 2] = -0.2 * (comp + 1);
      coeff[4 * comp + 3] = 0.1;
    }
  }

  gkyl_array_copy(arr, arr_ho);
  double result_ho;
  double *result =
    use_gpu ? (pinned ? gkyl_cu_malloc_host(sizeof(double)) : gkyl_cu_malloc(sizeof(double))) :
              &result_ho;

  const enum gkyl_array_op ops[] = {GKYL_MIN, GKYL_MAX, GKYL_SUM};
  for (int selection = 0; selection < 4; ++selection) {
    const struct gkyl_range *range = selection == 0 ? &parent : selection == 1 ? &local : &empty;

    // Only the host size is consulted for empty whole-array reductions.
    struct gkyl_array empty_arr = *arr;
    empty_arr.size = 0;

    for (int comp = 0; comp < 2; ++comp) {
      double expected[] = {DBL_MAX, -DBL_MAX, 0};
      gkyl_range_iter_init(&iter, range);
      while (gkyl_range_iter_next(&iter)) {
        const double *coeff = gkyl_array_cfetch(arr_ho, gkyl_range_idx(range, iter.idx));
        for (int node = 0; node < 4; ++node) {
          double coord[] = {gkyl_gauss_ordinates[2][node % 2], gkyl_gauss_ordinates[2][node / 2]};
          double value = basis_ho.eval_expand(coord, coeff + comp * basis_ho.num_basis);
          expected[0] = fmin(expected[0], value);
          expected[1] = fmax(expected[1], value);
          expected[2] += value;
        }
      }

      for (int repeat = 0; repeat < 3; ++repeat) {
        for (int operation = 0; operation < 3; ++operation) {
          result_ho = 123.0 + repeat;
          if (use_gpu && !pinned) {
            gkyl_cu_memcpy(result, &result_ho, sizeof(double), GKYL_CU_MEMCPY_H2D);
          } else {
            *result = result_ho;
          }
          if (selection == 0 || selection == 3) {
            gkyl_array_dg_reducec(
              result, selection == 3 ? &empty_arr : arr, comp, ops[operation], basis
            );
          } else {
            gkyl_array_dg_reducec_range(result, arr, comp, ops[operation], basis, range);
          }
          if (use_gpu && !pinned) {
            gkyl_cu_memcpy(&result_ho, result, sizeof(double), GKYL_CU_MEMCPY_D2H);
          } else {
            result_ho = *result;
          }
          double target = expected[operation];
          TEST_CHECK(
            selection >= 2 ? result_ho == target :
                             fabs(result_ho - target) <= 2e-12 * fmax(1.0, fabs(target))
          );
          TEST_MSG(
            "selection=%d component=%d op=%d repeat=%d pinned=%d: %.17g != %.17g", selection, comp,
            operation, repeat, pinned, result_ho, target
          );
        }
      }
    }
  }

  gkyl_array_release(arr);
  gkyl_array_release(arr_ho);
  if (use_gpu) {
    gkyl_cart_modal_basis_release_cu(basis);
    if (pinned) {
      gkyl_cu_free_host(result);
    } else {
      gkyl_cu_free(result);
    }
  }
}

static void
test_reduce_dg_edges_ho(void)
{
  test_reduce_dg_edges(false, false);
}

// CUDA specific tests
#ifdef GKYL_HAVE_CUDA

void
test_reduce_dg_dev()
{
  test_reduce_dg(true);
}

void
test_reduce_dg_range_dev()
{
  test_reduce_dg_range(true);
}

static void
test_reduce_dg_edges_dev()
{
  test_reduce_dg_edges(true, false);
  test_reduce_dg_edges(true, true);
}
#endif

TEST_LIST = {
  {"array_reduce_dg_edges_ho", test_reduce_dg_edges_ho},
  {"array_reduce_dg_ho", test_reduce_dg_ho},
  {"array_reduce_dg_range_ho", test_reduce_dg_range_ho},
#ifdef GKYL_HAVE_CUDA
  {"array_reduce_dg_edges_dev", test_reduce_dg_edges_dev},
  {"array_reduce_dg_dev", test_reduce_dg_dev},
  {"array_reduce_dg_range_dev", test_reduce_dg_range_dev},
#endif
  {NULL, NULL}
};
