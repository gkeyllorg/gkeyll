#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_util.h>
#include <time.h>
#include <string.h>

/**
 * Allocate a double array on the selected backend.
 *
 * @param use_gpu Whether to allocate on the GPU.
 * @param ncomp Number of components per cell.
 * @param size Number of cells.
 * @return New array; the caller must release it.
 */
static struct gkyl_array *
mkarr(bool use_gpu, int ncomp, long size)
{
  return use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, ncomp, size) :
                   gkyl_array_new(GKYL_DOUBLE, ncomp, size);
}

void
test_array_reduce_dummy_ho()
{
}

void
test_reduce_ho()
{
  int ncomp = 3, ncells = 200;
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, ncomp, ncells);

  for (size_t i = 0; i < arr->size; ++i) {
    double *d = gkyl_array_fetch(arr, i);
    for (size_t c = 0; c < ncomp; ++c) {
      d[c] = 0.5 * i + 0.1 * c;
    }
  }

  double amin[ncomp], amax[ncomp];
  gkyl_array_reduce(&amin[0], arr, GKYL_MIN);
  gkyl_array_reduce(&amax[0], arr, GKYL_MAX);

  for (size_t c = 0; c < ncomp; ++c) {
    TEST_CHECK(amin[c] == 0.1 * c);
    TEST_CHECK(amax[c] == 0.5 * (ncells - 1) + 0.1 * c);
  }

  gkyl_array_release(arr);
}

void
test_reduce_range_ho()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);

  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 3, range.volume);

  int count = -1000;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *d = gkyl_array_fetch(arr, loc);
    d[0] = count + 0.5;
    d[1] = count + 1.5;
    d[2] = count + 2.5;

    count += 1;
  }

  double amin[3], amax[3];
  gkyl_array_reduce_range(amin, arr, GKYL_MIN, &range);

  TEST_CHECK(amin[0] == -999.5);
  TEST_CHECK(amin[1] == -998.5);
  TEST_CHECK(amin[2] == -997.5);

  gkyl_array_reduce_range(amax, arr, GKYL_MAX, &range);

  TEST_CHECK(amax[0] == -800.5);
  TEST_CHECK(amax[1] == -799.5);
  TEST_CHECK(amax[2] == -798.5);

  gkyl_array_release(arr);
}

void
test_sum_reduce_range_ho()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);

  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 3, range.volume);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *d = gkyl_array_fetch(arr, loc);
    d[0] = 0.5;
    d[1] = 1.5;
    d[2] = 2.5;
  }

  double asum[3];
  gkyl_array_reduce_range(asum, arr, GKYL_SUM, &range);

  TEST_CHECK(asum[0] == 0.5 * range.volume);
  TEST_CHECK(asum[1] == 1.5 * range.volume);
  TEST_CHECK(asum[2] == 2.5 * range.volume);

  gkyl_array_release(arr);
}

static const enum gkyl_array_op extended_ops[] = {GKYL_MIN,     GKYL_MAX,    GKYL_SUM,
                                                  GKYL_ABS_MAX, GKYL_SQ_SUM, GKYL_CELL_MEAN_SQ_SUM};

/**
 * Check supported reductions against expected results, repeating GPU calls to expose races.
 *
 * @param use_gpu Whether to run on the GPU.
 * @param arr_ho Host input fixture with three components per cell.
 * @param wgt_ho Host weight fixture matching arr_ho.
 * @param range Selected cells, or NULL for the whole array.
 * @param expected Results indexed by extended_ops and component.
 * @param weighted Whether to use the weighted reduction API.
 * @param pinned Whether GPU results use pinned host storage instead of device storage.
 */
static void
check_reductions(
  bool use_gpu, const struct gkyl_array *arr_ho, const struct gkyl_array *wgt_ho,
  const struct gkyl_range *range, const double expected[6][3], bool weighted, bool pinned
)
{
  struct gkyl_array *arr = mkarr(use_gpu, 3, arr_ho->size);
  struct gkyl_array *wgt = mkarr(use_gpu, 3, wgt_ho->size);

  gkyl_array_copy(arr, arr_ho);
  gkyl_array_copy(wgt, wgt_ho);

  double result_ho[3];
#ifdef GKYL_DETERMINISTIC_REDUCTIONS_GPU
  double first_sum[3];
#endif
  double *result = use_gpu ? (pinned ? gkyl_cu_malloc_host(sizeof(result_ho)) :
                                       gkyl_cu_malloc(sizeof(result_ho))) :
                             result_ho;

  int repeats = use_gpu ? 5 : 1;
  for (int rep = 0; rep < repeats; ++rep) {
    for (int operation = 0; operation < (weighted ? 6 : 5); ++operation) {
      enum gkyl_array_op op = extended_ops[operation];
      if (weighted) {
        if (range) {
          gkyl_array_reduce_weighted_range(result, arr, wgt, op, range);
        } else {
          gkyl_array_reduce_weighted(result, arr, wgt, op);
        }
      } else {
        if (range) {
          gkyl_array_reduce_range(result, arr, op, range);
        } else {
          gkyl_array_reduce(result, arr, op);
        }
      }

      if (use_gpu) {
        if (pinned) {
          memcpy(result_ho, result, sizeof(result_ho));
        } else {
          gkyl_cu_memcpy(result_ho, result, sizeof(result_ho), GKYL_CU_MEMCPY_D2H);
        }
      }

      for (int comp = 0; comp < 3; ++comp) {
#ifdef GKYL_DETERMINISTIC_REDUCTIONS_GPU
        if (use_gpu && !weighted && op == GKYL_SUM) {
          if (rep == 0) {
            first_sum[comp] = result_ho[comp];
          } else {
            TEST_CHECK(memcmp(&result_ho[comp], &first_sum[comp], sizeof(double)) == 0);
            TEST_MSG(
              "sum changed on repeat %d, component %d: %.17g vs %.17g", rep, comp, result_ho[comp],
              first_sum[comp]
            );
          }
        }
#endif
        // Parallel sum ordering can differ, but extrema/empty identities are exact.
        double target = expected[operation][comp];
        if (operation == 0 || operation == 1 || operation == 3) {
          TEST_CHECK(result_ho[comp] == target);
        } else {
          TEST_CHECK(fabs(result_ho[comp] - target) <= 2e-12 * fmax(1.0, fabs(target)));
        }
        TEST_MSG(
          "weighted=%d op=%d comp=%d rep=%d: %.17g != %.17g", weighted, op, comp, rep,
          result_ho[comp], target
        );
      }
    }
  }

  if (use_gpu) {
    if (pinned) {
      gkyl_cu_free_host(result);
    } else {
      gkyl_cu_free(result);
    }
  }
  gkyl_array_release(arr);
  gkyl_array_release(wgt);
}

static void
array_reduce_extended_test(bool use_gpu)
{
  // Test reduction formulas, block-boundary sizes, subrange exclusion, and empty identities.

  // Hand-calculated results distinguish wx^2 from (wx)^2, and mean square from RMS.
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 3, 2);
  struct gkyl_array *wgt = gkyl_array_new(GKYL_DOUBLE, 3, 2);

  const double values[] = {1, -2, 3, -4, 5, -6};
  const double weights[] = {2, -3, 0.5, -1, 4, 2};
  memcpy(arr->data, values, sizeof(values));
  memcpy(wgt->data, weights, sizeof(weights));

  const double plain[6][3] = {{-4, -2, -6}, {1, 5, 3}, {-3, 3, -3}, {4, 5, 6}, {17, 29, 45}, {0}};
  const double weighted[6][3] = {{2, 6, -12},     {4, 20, 1.5},
                                 {6, 26, -10.5},  {4, 20, 12},
                                 {-14, 88, 76.5}, {-49.0 / 3.0, 266.0 / 3.0, 161.0 / 3.0}};

  struct gkyl_range full;
  gkyl_range_init(&full, 1, (int[]){1}, (int[]){2});

  for (int region = 0; region < 2; ++region) {
    check_reductions(use_gpu, arr, wgt, region ? &full : NULL, plain, false, false);
    check_reductions(use_gpu, arr, wgt, region ? &full : NULL, weighted, true, false);
  }

  if (use_gpu) {
    check_reductions(true, arr, wgt, NULL, plain, false, true);
    check_reductions(true, arr, wgt, &full, weighted, true, true);
  }

  gkyl_array_release(arr);
  gkyl_array_release(wgt);

  const int block_size = 256; // Match the default CUDA block size in CPU tests too.
  const int counts[] = {1, block_size - 1, block_size, block_size + 1, 5 * block_size + 13};
  for (int shape = 0; shape < 7; ++shape) {
    struct gkyl_range parent, selection;
    if (shape < 5) {
      gkyl_range_init(&parent, 1, (int[]){-2}, (int[]){counts[shape] - 3});
      selection = parent;
    } else {
      gkyl_range_init(&parent, 3, (int[]){-2, 1, 4}, (int[]){2, 5, 8});
      gkyl_sub_range_init(
        &selection, &parent, (int[]){-1, 2, 5}, shape == 5 ? (int[]){1, 4, 7} : (int[]){-2, 4, 7}
      );
    }

    struct gkyl_array *input = gkyl_array_new(GKYL_DOUBLE, 3, parent.volume);
    struct gkyl_array *weight = gkyl_array_new(GKYL_DOUBLE, 3, parent.volume);

    double expected[2][6][3] = {0};
    for (int kind = 0; kind < 2; ++kind) {
      for (int comp = 0; comp < 3; ++comp) {
        expected[kind][0][comp] = DBL_MAX;
        expected[kind][1][comp] = -DBL_MAX;
      }
    }

    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &parent);
    while (gkyl_range_iter_next(&iter)) {
      long cell = gkyl_range_idx(&parent, iter.idx);
      bool selected = selection.volume > 0 && gkyl_range_contains_idx(&selection, iter.idx);

      double *data = gkyl_array_fetch(input, cell);
      double *wdata = gkyl_array_fetch(weight, cell);

      double sum_sq = 0.0;
      for (int comp = 0; comp < 3; ++comp) {
        data[comp] = (cell % 2 ? -1.0 : 1.0) * (1.0 + (cell % 19) / 8.0 + comp / 4.0);
        // Extrema in different CUDA blocks, plus large excluded-cell sentinels.
        if (cell == 0 && comp == 0) {
          data[comp] = -8192;
        }
        if (cell == parent.volume - 1 && comp == 1) {
          data[comp] = 16384;
        }
        if (!selected) {
          data[comp] = 1e8;
        }
        wdata[comp] = (comp == 1 ? -1.0 : 1.0) * (0.25 + (cell % 7) / 16.0 + comp / 8.0);
        sum_sq += data[comp] * data[comp];
      }

      if (selected) {
        for (int kind = 0; kind < 2; ++kind) {
          for (int comp = 0; comp < 3; ++comp) {
            double factor = kind ? wdata[comp] : 1.0;
            double product = factor * data[comp];
            expected[kind][0][comp] = fmin(expected[kind][0][comp], product);
            expected[kind][1][comp] = fmax(expected[kind][1][comp], product);
            expected[kind][2][comp] += product;
            expected[kind][3][comp] = fmax(expected[kind][3][comp], fabs(product));
            expected[kind][4][comp] += factor * data[comp] * data[comp];
            expected[kind][5][comp] += factor * (sum_sq / 3.0);
          }
        }
      }
    }

    check_reductions(use_gpu, input, weight, &selection, expected[0], false, false);
    check_reductions(use_gpu, input, weight, &selection, expected[1], true, false);

    if (shape < 5) {
      check_reductions(use_gpu, input, weight, NULL, expected[0], false, false);
      check_reductions(use_gpu, input, weight, NULL, expected[1], true, false);
    }

    gkyl_array_release(input);
    gkyl_array_release(weight);
  }
}

static void
test_array_reduce_extended_ho(void)
{
  array_reduce_extended_test(false);
}

// Fast-math permits reassociation, so evaluation-order checks require a strict build.
#ifndef __FAST_MATH__
/**
 * Check weighted squared sums whose unweighted squares overflow or underflow.
 *
 * @param use_gpu Whether to run on the GPU.
 */
static void
array_reduce_weighted_sq_sum_extreme_test(bool use_gpu)
{
  struct gkyl_array *arr_ho = gkyl_array_new(GKYL_DOUBLE, 2, 3);
  struct gkyl_array *wgt_ho = gkyl_array_new(GKYL_DOUBLE, 2, 3);
  const double values[] = {0, 0, -1e200, 1e-200, 0, 0};
  const double weights[] = {1, 1, 1e-200, -1e200, 1, 1};
  const double expected[] = {1e200, -1e-200};
  memcpy(arr_ho->data, values, sizeof(values));
  memcpy(wgt_ho->data, weights, sizeof(weights));

  struct gkyl_array *arr = mkarr(use_gpu, 2, 3);
  struct gkyl_array *wgt = mkarr(use_gpu, 2, 3);
  gkyl_array_copy(arr, arr_ho);
  gkyl_array_copy(wgt, wgt_ho);
  struct gkyl_range parent, selection;
  gkyl_range_init(&parent, 1, (int[]){-1}, (int[]){1});
  gkyl_sub_range_init(&selection, &parent, (int[]){0}, (int[]){0});

  double result_ho[2];
  double *result = use_gpu ? gkyl_cu_malloc(sizeof(result_ho)) : result_ho;
  for (int region = 0; region < 2; ++region) {
    if (region) {
      gkyl_array_reduce_weighted_range(result, arr, wgt, GKYL_SQ_SUM, &selection);
    } else {
      gkyl_array_reduce_weighted(result, arr, wgt, GKYL_SQ_SUM);
    }
    if (use_gpu) {
      gkyl_cu_memcpy(result_ho, result, sizeof(result_ho), GKYL_CU_MEMCPY_D2H);
    }
    for (int comp = 0; comp < 2; ++comp) {
      // A purely relative check must reject zero for the tiny nonzero result.
      TEST_CHECK(fabs((result_ho[comp] - expected[comp]) / expected[comp]) < 1e-14);
      TEST_MSG("region=%d comp=%d: %.17g != %.17g", region, comp, result_ho[comp], expected[comp]);
    }
  }
  if (use_gpu) {
    gkyl_cu_free(result);
  }
  gkyl_array_release(arr);
  gkyl_array_release(wgt);
  gkyl_array_release(arr_ho);
  gkyl_array_release(wgt_ho);
}

/**
 * Run the extreme weighted squared-sum checks on the CPU in strict builds.
 */
static void
test_array_reduce_weighted_sq_sum_extreme_ho(void)
{
  array_reduce_weighted_sq_sum_extreme_test(false);
}
#endif

// CUDA specific tests
#ifdef GKYL_HAVE_CUDA

void
test_array_reduce_max_dev()
{
  unsigned long numComp = 3, numCells = 10;
  // Create host array and device copy
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, numComp, numCells);
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, numComp, numCells);
  // Initialize data
  double *a1_d = a1->data;
  for (unsigned i = 0; i < numCells; ++i) {
    for (unsigned k = 0; k < numComp; ++k) {
      a1_d[i * numComp + k] = (double)i + (double)k * 0.10;
    }
  }
  gkyl_array_copy(a1_cu, a1);
  // create a device and host arrays to store reduction
  double *a1max = gkyl_malloc(numComp * sizeof(double));
  double *a1max_cu = (double *)gkyl_cu_malloc(numComp * sizeof(double));

  // Component-wise reduce array.
  gkyl_array_reduce(a1max_cu, a1_cu, GKYL_MAX);

  // Copy to host and check values.
  gkyl_cu_memcpy(a1max, a1max_cu, numComp * sizeof(double), GKYL_CU_MEMCPY_D2H);
  for (unsigned k = 0; k < numComp; ++k) {
    TEST_CHECK(gkyl_compare(a1max[k], (double)(numCells - 1) + (double)k * 0.10, 1e-14));
  }

  gkyl_free(a1max);
  gkyl_cu_free(a1max_cu);
  gkyl_array_release(a1);
  gkyl_array_release(a1_cu);
}

void
test_array_reduce_max_big_dev()
{
  unsigned long numComp = 12, numCells = 1000;
  // Create host array and device copy
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, numComp, numCells);
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, numComp, numCells);
  // Initialize data
  double *a1_d = a1->data;
  for (unsigned i = 0; i < numCells; ++i) {
    for (unsigned k = 0; k < numComp; ++k) {
      a1_d[i * numComp + k] = (double)i + (double)k * 0.10;
    }
  }
  gkyl_array_copy(a1_cu, a1);
  // create a device and host arrays to store reduction
  double *a1max = gkyl_malloc(numComp * sizeof(double));
  double *a1max_cu = (double *)gkyl_cu_malloc(numComp * sizeof(double));

  // Component-wise reduce array.
  gkyl_array_reduce(a1max_cu, a1_cu, GKYL_MAX);

  // Copy to host and check values.
  gkyl_cu_memcpy(a1max, a1max_cu, numComp * sizeof(double), GKYL_CU_MEMCPY_D2H);
  for (unsigned k = 0; k < numComp; ++k) {
    // Make print statements to manually check the comparison
    // printf("a1max[%d] = %g\n", k, a1max[k]);
    // printf("numCells-1+(double)k*0.10 = %g\n", (double)(numCells-1)+(double)k*0.10);
    TEST_CHECK(gkyl_compare(a1max[k], (double)(numCells - 1) + (double)k * 0.10, 1e-14));
  }

  gkyl_free(a1max);
  gkyl_cu_free(a1max_cu);
  gkyl_array_release(a1);
  gkyl_array_release(a1_cu);
}

void
test_array_reduce_range_1d_max_dev()
{
  unsigned long numComp = 1, numCells = 10;
  // create host array and device copy
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, numComp, numCells);
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, numComp, numCells);
  // initialize data
  double *a1_d = a1->data;
  for (unsigned i = 0; i < numCells; ++i) {
    for (unsigned k = 0; k < numComp; ++k) {
      a1_d[i * numComp + k] = (double)i + (double)k * 0.10;
    }
  }
  gkyl_array_copy(a1_cu, a1);
  // create a device and host arrays to store reduction
  double *a1max = gkyl_malloc(numComp * sizeof(double));
  double *a1max_cu = (double *)gkyl_cu_malloc(numComp * sizeof(double));

  int lower[] = {1}, upper[] = {numCells};
  int sublower[] = {3}, subupper[] = {7};
  struct gkyl_range range, subrange;
  gkyl_range_init(&range, 1, lower, upper);
  gkyl_sub_range_init(&subrange, &range, sublower, subupper);

  // Component-wise reduce array on range
  gkyl_array_reduce_range(a1max_cu, a1_cu, GKYL_MAX, &range);
  // Copy to host and check values.
  gkyl_cu_memcpy(a1max, a1max_cu, numComp * sizeof(double), GKYL_CU_MEMCPY_D2H);
  for (unsigned k = 0; k < numComp; ++k) {
    TEST_CHECK(gkyl_compare(a1max[k], (double)(upper[0] - 1) + (double)k * 0.10, 1e-14));
  }

  // Component-wise reduce array on subrange
  gkyl_cu_memset(a1max_cu, 0, sizeof(double) * numComp);
  gkyl_array_reduce_range(a1max_cu, a1_cu, GKYL_MAX, &subrange);

  // Copy to host and check values.
  gkyl_cu_memcpy(a1max, a1max_cu, numComp * sizeof(double), GKYL_CU_MEMCPY_D2H);
  for (unsigned k = 0; k < numComp; ++k) {
    TEST_CHECK(gkyl_compare(a1max[k], (double)(subupper[0] - 1) + (double)k * 0.10, 1e-14));
  }

  gkyl_free(a1max);
  gkyl_cu_free(a1max_cu);
  gkyl_array_release(a1);
  gkyl_array_release(a1_cu);
}

void
test_array_reduce_range_2d_max_dev()
{
  int cells[] = {8, 10};
  int ghost[] = {1, 0};
  double lower[] = {0., -1.};
  double upper[] = {1., 1.};

  struct gkyl_rect_grid grid;
  struct gkyl_range range, range_ext;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);
  gkyl_create_grid_ranges(&grid, ghost, &range_ext, &range);

  unsigned long numComp = 1;
  unsigned long numCells = range_ext.volume;
  // create host array and device copy
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, numComp, numCells);
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, numComp, numCells);
  struct gkyl_array *b1 = gkyl_array_new(GKYL_DOUBLE, 1, numCells);
  struct gkyl_array *b1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, numCells);
  // initialize data
  double *a1_d = a1->data;
  double *b1_d = b1->data;
  for (unsigned i = 0; i < numCells; ++i) {
    b1_d[i] = (double)i + 1000;
    for (unsigned k = 0; k < numComp; ++k) {
      a1_d[i * numComp + k] = (double)numCells - 1 - i + (double)k * 0.10;
    }
  }
  gkyl_array_copy(a1_cu, a1);
  gkyl_array_copy(b1_cu, b1);
  // create a device and host arrays to store reduction
  double *a1max = gkyl_malloc(numComp * sizeof(double));
  double *a1max_correct = gkyl_malloc(numComp * sizeof(double));
  double *a1max_cu = (double *)gkyl_cu_malloc(numComp * sizeof(double));
  double *b1max = gkyl_malloc(1 * sizeof(double));
  double *b1max_cu = (double *)gkyl_cu_malloc(1 * sizeof(double));
  gkyl_cu_memset(a1max_cu, 0, sizeof(double) * numComp);
  gkyl_cu_memset(b1max_cu, 0, sizeof(double) * 1);

  // Component-wise reduce array on range
  gkyl_array_reduce_range(b1max_cu, b1_cu, GKYL_MAX, &range_ext);
  // Copy to host and check values.
  gkyl_cu_memcpy(b1max, b1max_cu, sizeof(double), GKYL_CU_MEMCPY_D2H);
  TEST_CHECK(gkyl_compare(b1max[0], (double)(1000 + numCells - 1), 1e-14));

  // Component-wise reduce array on subrange
  gkyl_cu_memset(a1max_cu, 0, sizeof(double) * numComp);
  gkyl_array_reduce_range(a1max_cu, a1_cu, GKYL_MAX, &range);

  gkyl_array_reduce_range(a1max_correct, a1, GKYL_MAX, &range);

  // Copy to host and check values.
  gkyl_cu_memcpy(a1max, a1max_cu, numComp * sizeof(double), GKYL_CU_MEMCPY_D2H);
  for (unsigned k = 0; k < numComp; ++k) {
    TEST_CHECK(gkyl_compare(a1max[k], a1max_correct[k], 1e-14));
  }

  gkyl_free(a1max);
  gkyl_cu_free(a1max_cu);
  gkyl_array_release(a1);
  gkyl_array_release(a1_cu);

  gkyl_cu_free(b1max);
  gkyl_cu_free(b1max_cu);
  gkyl_array_release(b1);
  gkyl_array_release(b1_cu);
}

void
test_cu_array_reduce_range_max_timer(int NX, int NY, int VX, int VY)
{
  int cells[] = {NX, NY, VX, VY};
  int ghost[] = {1, 1, 0, 0};
  double lower[] = {0.0, 0.0, 0.0, 0.0};
  double upper[] = {1.0, 1.0, 1.0, 1.0};

  int ndim = 4;
  struct gkyl_rect_grid grid;
  struct gkyl_range range, range_ext;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  gkyl_create_grid_ranges(&grid, ghost, &range_ext, &range);

  unsigned long numComp = 1;
  unsigned long numCells = range_ext.volume;
  // create host array and device copy
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, numComp, numCells);
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, numComp, numCells);
  struct gkyl_array *b1 = gkyl_array_new(GKYL_DOUBLE, 1, numCells);
  struct gkyl_array *b1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, numCells);
  // initialize data
  double *a1_d = a1->data;
  double *b1_d = b1->data;
  for (unsigned i = 0; i < numCells; ++i) {
    b1_d[i] = (double)i + 1000;
    for (unsigned k = 0; k < numComp; ++k) {
      a1_d[i * numComp + k] = (double)numCells - 1 - i + (double)k * 0.10;
    }
  }
  gkyl_array_copy(a1_cu, a1);
  gkyl_array_copy(b1_cu, b1);
  // create a device and host arrays to store reduction
  double *a1max = gkyl_malloc(numComp * sizeof(double));
  double *a1max_correct = gkyl_malloc(numComp * sizeof(double));
  double *a1max_cu = (double *)gkyl_cu_malloc(numComp * sizeof(double));

  double *b1max = gkyl_malloc(1 * sizeof(double));
  double *b1max_cu = (double *)gkyl_cu_malloc(1 * sizeof(double));

  gkyl_cu_memset(a1max_cu, 0, sizeof(double) * numComp);
  gkyl_cu_memset(b1max_cu, 0, sizeof(double) * 1);

  // Component-wise reduce array on subrange
  gkyl_cu_memset(a1max_cu, 0, sizeof(double) * numComp);

  struct timespec tm = gkyl_wall_clock();

  for (int i = 0; i < 100; ++i) {
    gkyl_array_reduce_range(a1max_cu, a1_cu, GKYL_MAX, &range);
  }
  double red_tm = gkyl_time_diff_now_sec(tm);

  // printf("100 reductions on (%d,%d,%d,%d) took %g sec\n", NX, NY, VX, VY,
  // red_tm);

  gkyl_free(a1max);
  gkyl_cu_free(a1max_cu);
  gkyl_array_release(a1);
  gkyl_array_release(a1_cu);

  gkyl_cu_free(b1max);
  gkyl_cu_free(b1max_cu);
  gkyl_array_release(b1);
  gkyl_array_release(b1_cu);
}

void
test_array_reduce_range_max_timer_32x32x40x40_dev()
{
  test_cu_array_reduce_range_max_timer(32, 32, 40, 40);
}

void
test_array_reduce_range_max_timer_32x32x32x32_dev()
{
  test_cu_array_reduce_range_max_timer(32, 32, 32, 32);
}

static void
test_array_reduce_extended_dev(void)
{
  array_reduce_extended_test(true);
}
#ifdef GKYL_DETERMINISTIC_REDUCTIONS_GPU
/**
 * Check that multiblock sums of interleaved components and a strided subrange repeat exactly.
 */
static void
test_array_reduce_sum_deterministic_dev(void)
{
  struct gkyl_range parent, selection;
  gkyl_range_init(&parent, 2, (int[]){1, 1}, (int[]){37, 37});
  gkyl_sub_range_init(&selection, &parent, (int[]){2, 2}, (int[]){36, 36});

  struct gkyl_array *host = gkyl_array_new(GKYL_DOUBLE, 3, parent.volume);
  for (long cell = 0; cell < parent.volume; ++cell) {
    double *value = gkyl_array_fetch(host, cell);
    for (int comp = 0; comp < 3; ++comp) {
      value[comp] = ((cell % 17) - 8) * 0.1 / (comp + 1);
    }
  }
  struct gkyl_array *device = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3, parent.volume);
  gkyl_array_copy(device, host);
  double *result_dev = gkyl_cu_malloc(3 * sizeof(double));

  for (int use_range = 0; use_range < 2; ++use_range) {
    double expected[3], first[3], result[3];
    if (use_range) {
      gkyl_array_reduce_range(expected, host, GKYL_SUM, &selection);
    } else {
      gkyl_array_reduce(expected, host, GKYL_SUM);
    }
    for (int rep = 0; rep < 5; ++rep) {
      if (use_range) {
        gkyl_array_reduce_range(result_dev, device, GKYL_SUM, &selection);
      } else {
        gkyl_array_reduce(result_dev, device, GKYL_SUM);
      }
      gkyl_cu_memcpy(result, result_dev, sizeof(result), GKYL_CU_MEMCPY_D2H);
      for (int comp = 0; comp < 3; ++comp) {
        TEST_CHECK(fabs(result[comp] - expected[comp]) < 1e-11);
        if (rep == 0) {
          first[comp] = result[comp];
        } else {
          TEST_CHECK(memcmp(&result[comp], &first[comp], sizeof(double)) == 0);
        }
      }
    }
  }

  gkyl_cu_free(result_dev);
  gkyl_array_release(device);
  gkyl_array_release(host);
}
#endif
#ifndef __FAST_MATH__
/**
 * Run the extreme weighted squared-sum checks on the GPU in strict builds.
 */
static void
test_array_reduce_weighted_sq_sum_extreme_dev(void)
{
  array_reduce_weighted_sq_sum_extreme_test(true);
}
#endif
#endif

TEST_LIST = {
  {"array_reduce_dummy_ho", test_array_reduce_dummy_ho},
  {"array_reduce_ho", test_reduce_ho},
  {"array_reduce_range_ho", test_reduce_range_ho},
  {"array_reduce_sum_range_ho", test_sum_reduce_range_ho},
  {"array_reduce_extended_ho", test_array_reduce_extended_ho},
#ifndef __FAST_MATH__
  {"array_reduce_weighted_sq_sum_extreme_ho", test_array_reduce_weighted_sq_sum_extreme_ho},
#endif
#ifdef GKYL_HAVE_CUDA
  {"array_reduce_max_dev", test_array_reduce_max_dev},
  {"array_reduce_max_big_dev", test_array_reduce_max_big_dev},
  {"array_reduce_range_1d_max_dev", test_array_reduce_range_1d_max_dev},
  {"array_reduce_range_2d_max_dev", test_array_reduce_range_2d_max_dev},
  {"array_reduce_range_max_timer_32x32x40x40_dev", test_array_reduce_range_max_timer_32x32x40x40_dev
  },
  {"array_reduce_range_max_timer_32x32x32x32_dev", test_array_reduce_range_max_timer_32x32x32x32_dev
  },
  {"array_reduce_extended_dev", test_array_reduce_extended_dev},
#ifdef GKYL_DETERMINISTIC_REDUCTIONS_GPU
  {"array_reduce_sum_deterministic_dev", test_array_reduce_sum_deterministic_dev},
#endif
#ifndef __FAST_MATH__
  {"array_reduce_weighted_sq_sum_extreme_dev", test_array_reduce_weighted_sq_sum_extreme_dev},
#endif
#endif
  {NULL, NULL}
};
