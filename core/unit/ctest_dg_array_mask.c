#include <acutest.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_dg_array_mask.h>
#include <gkyl_dg_array_mask_priv.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

#include <math.h>

// Helper function to create test arrays
static struct gkyl_array *
mkarr(bool use_gpu, long nc, long size)
{
  struct gkyl_array *a = use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
                                 : gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

// Test basic mask creation and initialization
void test_mask_new(bool use_gpu)
{
  int shape[] = { 10, 20 };
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);

  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS,
    .threshold = 1e-10,
    .conf_rng = &range,
    .conf_rng_ext = &range,
    .use_gpu = use_gpu
  };

  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);

  TEST_CHECK(mask != NULL);
  TEST_CHECK(mask->type == GKYL_DG_ARRAY_MASK_C0_LESS);
  TEST_CHECK(mask->use_gpu == use_gpu);
  TEST_CHECK(mask->mask_arr != NULL);
  TEST_CHECK(mask->mask_arr->size == range.volume);
  TEST_CHECK(mask->mask_arr->ncomp == 1);

  // Check that mask is initialized to -1.0 (false)
  const struct gkyl_array *mask_arr = gkyl_dg_array_mask_get_mask(mask);
  struct gkyl_array *mask_ho = mkarr(false, mask_arr->ncomp, mask_arr->size);
  gkyl_array_copy(mask_ho, mask_arr);

  for (unsigned i = 0; i < range.volume; ++i) {
    const double *mask_d = gkyl_array_cfetch(mask_ho, i);
    TEST_CHECK(gkyl_compare(mask_d[0], -1.0, 1e-14) );
  }

  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

// Test mask with NONE type
void test_mask_none_type(bool use_gpu)
{
  struct gkyl_dg_array_mask_inp mask_inp = {
    .use_gpu = use_gpu
  };
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);

  TEST_CHECK(mask->type == GKYL_DG_ARRAY_MASK_NONE);
  TEST_CHECK(mask->use_gpu == use_gpu);
  TEST_CHECK(mask->threshold == 0.0);

  gkyl_dg_array_mask_release(mask);
}

// Test mask advance with threshold
void test_mask_advance_threshold(bool use_gpu)
{
  int shape[] = { 10 };
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 1, shape);

  double threshold = 0.5;
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS,
    .threshold = threshold,
    .conf_rng = &range,
    .conf_rng_ext = &range,
    .use_gpu = use_gpu
  };

  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);

  // Create test array on host with varying values
  struct gkyl_array *arr = mkarr(use_gpu, 1, range.volume);
  struct gkyl_array *arr_ho = mkarr(false, arr->ncomp, arr->size);
  for (unsigned i = 0; i < range.volume; ++i) {
    // Values alternate above and below threshold
    double *arr_d = gkyl_array_fetch(arr_ho, i);
    arr_d[0] = (i % 2 == 0) ? 0.0 : 10.0;
  }

  // Copy to device if needed and run advance
  gkyl_array_copy(arr, arr_ho);

  double global_max = 0.0;
  gkyl_array_reduce(&global_max, arr_ho, GKYL_MAX);
  gkyl_dg_array_mask_advance_threshold(mask, global_max);
  gkyl_dg_array_mask_advance(mask, arr);

  // Copy mask back to host for verification
  const struct gkyl_array *mask_arr = gkyl_dg_array_mask_get_mask(mask);
  struct gkyl_array *mask_ho = mkarr(false, mask_arr->ncomp, mask_arr->size);
  gkyl_array_copy(mask_ho, mask_arr);

  // Check mask values
  double expected_threshold = threshold * sqrt(2.0); // 1D scaling
  for (unsigned i = 0; i < range.volume; ++i) {
    const double *mask_d = gkyl_array_cfetch(mask_ho, i);
    const double *arr_d = gkyl_array_cfetch(arr_ho, i);
    if (fabs(arr_d[0]) < expected_threshold) {
      TEST_CHECK(gkyl_compare(mask_d[0], 1.0, 1e-14) );  // True (masked)
    }
    else {
      TEST_CHECK(gkyl_compare(mask_d[0], -1.0, 1e-14) );  // False (not masked)
    }
  }

  gkyl_array_release(arr_ho);
  gkyl_array_release(arr);
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

// Test mask advance with all values below threshold
void test_mask_advance_all_below(bool use_gpu)
{
  int shape[] = { 5, 5 };
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);

  double threshold = 10.0;
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS,
    .threshold = threshold,
    .conf_rng = &range,
    .conf_rng_ext = &range,
    .use_gpu = use_gpu
  };

  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);

  // Create test array with all values below threshold
  struct gkyl_array *arr = mkarr(use_gpu, 1, range.volume);
  struct gkyl_array *arr_ho = mkarr(false, arr->ncomp, arr->size);
  gkyl_array_clear(arr_ho, 1e-12);
  gkyl_array_copy(arr, arr_ho);
  
  double global_max = 0.0;
  gkyl_array_reduce(&global_max, arr_ho, GKYL_MAX);
  gkyl_dg_array_mask_advance_threshold(mask, global_max);
  gkyl_dg_array_mask_advance(mask, arr);

  // Copy mask back to host for verification
  const struct gkyl_array *mask_arr = gkyl_dg_array_mask_get_mask(mask);
  struct gkyl_array *mask_ho = mkarr(false, mask_arr->ncomp, mask_arr->size);
  gkyl_array_copy(mask_ho, mask_arr);

  // All cells should be masked (1.0)
  for (unsigned i = 0; i < range.volume; ++i) {
    const double *mask_d = gkyl_array_cfetch(mask_ho, i);
    TEST_CHECK(gkyl_compare(mask_d[0], 1.0, 1e-14) );
  }

  gkyl_array_release(arr_ho);
  gkyl_array_release(arr);
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

// Test mask advance with all values above threshold
void test_mask_advance_all_above(bool use_gpu)
{
  int shape[] = { 5, 5 };
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);

  double threshold = 1e-15;
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS,
    .threshold = threshold,
    .conf_rng = &range,
    .conf_rng_ext = &range,
    .use_gpu = use_gpu
  };

  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);

  // Create test array with all values above threshold
  struct gkyl_array *arr = mkarr(use_gpu, 1, range.volume);
  struct gkyl_array *arr_ho = mkarr(false, arr->ncomp, arr->size);
  gkyl_array_clear(arr_ho, 100.0);
  gkyl_array_copy(arr, arr_ho);
  
  double global_max = 0.0;
  gkyl_array_reduce(&global_max, arr_ho, GKYL_MAX);
  gkyl_dg_array_mask_advance_threshold(mask, global_max);
  gkyl_dg_array_mask_advance(mask, arr);

  // Copy mask back to host for verification
  const struct gkyl_array *mask_arr = gkyl_dg_array_mask_get_mask(mask);
  struct gkyl_array *mask_ho = mkarr(false, mask_arr->ncomp, mask_arr->size);
  gkyl_array_copy(mask_ho, mask_arr);

  // No cells should be masked (-1.0)
  for (unsigned i = 0; i < range.volume; ++i) {
    const double *mask_d = gkyl_array_cfetch(mask_ho, i);
    TEST_CHECK(gkyl_compare(mask_d[0], -1.0, 1e-14) );
  }

  gkyl_array_release(arr_ho);
  gkyl_array_release(arr);
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

// Test mask advance with negative values
void test_mask_advance_negative_values(bool use_gpu)
{
  int shape[] = { 8 };
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 1, shape);

  double threshold = 0.5;
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS,
    .threshold = threshold,
    .conf_rng = &range,
    .conf_rng_ext = &range,
    .use_gpu = use_gpu
  };

  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);

  // Create test array with negative and positive values
  struct gkyl_array *arr = mkarr(use_gpu, 1, range.volume);
  struct gkyl_array *arr_ho = mkarr(false, arr->ncomp, arr->size);
  double vals[] = { 0.1, -0.1, 1.0, -1.0, 0.3, -0.3, 2.0, -2.0 };
  for (unsigned i = 0; i < range.volume; ++i) {
    double *arr_d = gkyl_array_fetch(arr_ho, i);
    arr_d[0] = vals[i];
  }
  gkyl_array_copy(arr, arr_ho);
  
  double global_max = 0.0;
  gkyl_array_reduce(&global_max, arr_ho, GKYL_MAX);
  gkyl_dg_array_mask_advance_threshold(mask, global_max);
  gkyl_dg_array_mask_advance(mask, arr);

  // Copy mask back to host for verification
  const struct gkyl_array *mask_arr = gkyl_dg_array_mask_get_mask(mask);
  struct gkyl_array *mask_ho = mkarr(false, mask_arr->ncomp, mask_arr->size);
  gkyl_array_copy(mask_ho, mask_arr);

  // Check mask values (uses fabs, so negative values treated as positive)
  double expected_threshold = threshold * sqrt(2.0);
  for (unsigned i = 0; i < range.volume; ++i) {
    const double *mask_d = gkyl_array_cfetch(mask_ho, i);
    const double *arr_d = gkyl_array_cfetch(arr_ho, i);
    if (fabs(arr_d[0]) < expected_threshold) {
      TEST_CHECK(gkyl_compare(mask_d[0], 1.0, 1e-14) );
    }
    else {
      TEST_CHECK(gkyl_compare(mask_d[0], -1.0, 1e-14) );
    }
  }

  gkyl_array_release(arr_ho);
  gkyl_array_release(arr);
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

// Test mask advance with GREATER_THAN threshold
void test_mask_advance_greater_than_threshold(bool use_gpu)
{
  int shape[] = { 10 };
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 1, shape);

  double threshold = 0.5;
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_GREATER,
    .threshold = threshold,
    .conf_rng = &range,
    .conf_rng_ext = &range,
    .use_gpu = use_gpu
  };

  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);

  // Create test array with varying values
  struct gkyl_array *arr = mkarr(use_gpu, 1, range.volume);
  struct gkyl_array *arr_ho = mkarr(false, arr->ncomp, arr->size);
  for (unsigned i = 0; i < range.volume; ++i) {
    // Values alternate above and below threshold
    double *arr_d = gkyl_array_fetch(arr_ho, i);
    arr_d[0] = (i % 2 == 0) ? 0.1 : 1.0;
  }
  gkyl_array_copy(arr, arr_ho);
  
  double global_max = 0.0;
  gkyl_array_reduce(&global_max, arr_ho, GKYL_MAX);
  gkyl_dg_array_mask_advance_threshold(mask, global_max);
  gkyl_dg_array_mask_advance(mask, arr);

  // Copy mask back to host for verification
  const struct gkyl_array *mask_arr = gkyl_dg_array_mask_get_mask(mask);
  struct gkyl_array *mask_ho = mkarr(false, mask_arr->ncomp, mask_arr->size);
  gkyl_array_copy(mask_ho, mask_arr);

  // Check mask values (opposite of LESS_THAN)
  double expected_threshold = threshold * sqrt(2.0); // 1D scaling
  for (unsigned i = 0; i < range.volume; ++i) {
    const double *mask_d = gkyl_array_cfetch(mask_ho, i);
    const double *arr_d = gkyl_array_cfetch(arr_ho, i);
    if (fabs(arr_d[0]) > expected_threshold) {
      TEST_CHECK(gkyl_compare(mask_d[0], 1.0, 1e-14) );  // True (masked)
    }
    else {
      TEST_CHECK(gkyl_compare(mask_d[0], -1.0, 1e-14) );  // False (not masked)
    }
  }

  gkyl_array_release(arr_ho);
  gkyl_array_release(arr);
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

// Test GREATER_THAN mask with all values above threshold
void test_mask_advance_greater_than_all_above(bool use_gpu)
{
  int shape[] = { 5, 5 };
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);

  double threshold = 1e-15;
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_GREATER,
    .threshold = threshold,
    .conf_rng = &range,
    .conf_rng_ext = &range,
    .use_gpu = use_gpu
  };

  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);

  // Create test array with all values above threshold
  struct gkyl_array *arr = mkarr(use_gpu, 1, range.volume);
  struct gkyl_array *arr_ho = mkarr(false, arr->ncomp, arr->size);
  gkyl_array_clear(arr_ho, 100.0);
  gkyl_array_copy(arr, arr_ho);

  double global_max = 0.0;
  gkyl_array_reduce(&global_max, arr_ho, GKYL_MAX);
  gkyl_dg_array_mask_advance_threshold(mask, global_max);
  gkyl_dg_array_mask_advance(mask, arr);

  // Copy mask back to host for verification
  const struct gkyl_array *mask_arr = gkyl_dg_array_mask_get_mask(mask);
  struct gkyl_array *mask_ho = mkarr(false, mask_arr->ncomp, mask_arr->size);
  gkyl_array_copy(mask_ho, mask_arr);

  // All cells should be masked (1.0) since values > threshold
  for (unsigned i = 0; i < range.volume; ++i) {
    const double *mask_d = gkyl_array_cfetch(mask_ho, i);
    TEST_CHECK(gkyl_compare(mask_d[0], 1.0, 1e-14) );
  }

  gkyl_array_release(arr_ho);
  gkyl_array_release(arr);
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

// Test GREATER_THAN mask with all values below threshold
void test_mask_advance_greater_than_all_below(bool use_gpu)
{
  int shape[] = { 5, 5 };
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);

  double threshold = 10.0;
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_GREATER,
    .threshold = threshold,
    .conf_rng = &range,
    .conf_rng_ext = &range,
    .use_gpu = use_gpu
  };

  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);

  // Create test array with all values below threshold
  struct gkyl_array *arr = mkarr(use_gpu, 1, range.volume);
  struct gkyl_array *arr_ho = mkarr(false, arr->ncomp, arr->size);
  gkyl_array_clear(arr_ho, 1e-12);
  gkyl_array_copy(arr, arr_ho);
  
  double global_max = 0.0;
  gkyl_array_reduce(&global_max, arr_ho, GKYL_MAX);
  gkyl_dg_array_mask_advance_threshold(mask, global_max);
  gkyl_dg_array_mask_advance(mask, arr);

  // Copy mask back to host for verification
  const struct gkyl_array *mask_arr = gkyl_dg_array_mask_get_mask(mask);
  struct gkyl_array *mask_ho = mkarr(false, mask_arr->ncomp, mask_arr->size);
  gkyl_array_copy(mask_ho, mask_arr);

  // No cells should be masked (-1.0) since values < threshold
  for (unsigned i = 0; i < range.volume; ++i) {
    const double *mask_d = gkyl_array_cfetch(mask_ho, i);
    TEST_CHECK(gkyl_compare(mask_d[0], -1.0, 1e-14) );
  }

  gkyl_array_release(arr_ho);
  gkyl_array_release(arr);
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

// Test GREATER_THAN mask with negative values
void test_mask_advance_greater_than_negative_values(bool use_gpu)
{
  int shape[] = { 8 };
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 1, shape);

  double threshold = 0.5;
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_GREATER,
    .threshold = threshold,
    .conf_rng = &range,
    .conf_rng_ext = &range,
    .use_gpu = use_gpu
  };

  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);

  // Create test array with negative and positive values
  struct gkyl_array *arr = mkarr(use_gpu, 1, range.volume);
  struct gkyl_array *arr_ho = mkarr(false, arr->ncomp, arr->size);
  double vals[] = { 0.1, -0.1, 1.0, -1.0, 0.3, -0.3, 2.0, -2.0 };
  for (unsigned i = 0; i < range.volume; ++i) {
    double *arr_d = gkyl_array_fetch(arr_ho, i);
    arr_d[0] = vals[i];
  }
  gkyl_array_copy(arr, arr_ho);
  
  double global_max = 0.0;
  gkyl_array_reduce(&global_max, arr_ho, GKYL_MAX);
  gkyl_dg_array_mask_advance_threshold(mask, global_max);
  gkyl_dg_array_mask_advance(mask, arr);

  // Copy mask back to host for verification
  const struct gkyl_array *mask_arr = gkyl_dg_array_mask_get_mask(mask);
  struct gkyl_array *mask_ho = mkarr(false, mask_arr->ncomp, mask_arr->size);
  gkyl_array_copy(mask_ho, mask_arr);

  // Check mask values (uses fabs, so negative values treated as positive)
  double expected_threshold = threshold * sqrt(2.0);
  for (unsigned i = 0; i < range.volume; ++i) {
    const double *arr_d = gkyl_array_cfetch(arr_ho, i);
    const double *mask_d = gkyl_array_cfetch(mask_ho, i);
    if (fabs(arr_d[0]) > expected_threshold) {
      TEST_CHECK(gkyl_compare(mask_d[0], 1.0, 1e-14) );
    }
    else {
      TEST_CHECK(gkyl_compare(mask_d[0], -1.0, 1e-14) );
    }
  }

  gkyl_array_release(arr_ho);
  gkyl_array_release(arr);
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

// Test mask eval inline function
void test_mask_eval(bool use_gpu)
{
  int shape[] = { 6 };
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 1, shape);

  double threshold = 0.5;
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS,
    .threshold = threshold,
    .conf_rng = &range,
    .conf_rng_ext = &range,
    .use_gpu = use_gpu
  };

  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);

  // Use advance to set mask values via a test array
  struct gkyl_array *arr = mkarr(use_gpu, 1, range.volume);
  struct gkyl_array *arr_ho = mkarr(false, arr->ncomp, arr->size);
  // Set up values so that mask will be: true, false, true, false, true, false
  double *arr_d;
  double vals[] = { 0.1, 1.0, 0.2, 2.0, 0.3, 5.0 };
  for (int i = 0; i < shape[0]; i++) {
    arr_d = gkyl_array_fetch(arr_ho, i);
    arr_d[0] = vals[i];
  }
  gkyl_array_copy(arr, arr_ho);
  
  double global_max = 0.0;
  gkyl_array_reduce(&global_max, arr_ho, GKYL_MAX);
  gkyl_dg_array_mask_advance_threshold(mask, global_max);
  gkyl_dg_array_mask_advance(mask, arr);

  // Test eval function
  // This part does not work on the GPU because the boolean is evaluated and
  // returned on the GPU, but needs a host evaluation for the TEST_CHECK
  bool *mask_cond_ho, *mask_cond;
  mask_cond_ho = gkyl_malloc(sizeof(bool));
  if (use_gpu)
    mask_cond = gkyl_cu_malloc(sizeof(bool));
  else
    mask_cond = gkyl_malloc(sizeof(bool));

  bool mask_cond_ref[] = {true, false, true, false, true, false};
  for (int i=0; i<shape[0]; i++) {
    gkyl_dg_array_mask_eval_idx(mask, &i, mask_cond);
    if (use_gpu)
      gkyl_cu_memcpy(mask_cond_ho, mask_cond, sizeof(bool), GKYL_CU_MEMCPY_D2H);
    else
      memcpy(mask_cond_ho, mask_cond, sizeof(bool));

    TEST_CHECK(mask_cond_ho[0] == mask_cond_ref[i]);
    TEST_MSG("Got %d at i=%d. Expected %d\n",mask_cond_ho[0], i, mask_cond_ref[i]);
  }

  gkyl_free(mask_cond_ho);
  if (use_gpu)
    gkyl_cu_free(mask_cond);
  else
    gkyl_free(mask_cond);

  gkyl_array_release(arr_ho);
  gkyl_array_release(arr);
  gkyl_dg_array_mask_release(mask);
}

// Test mask eval with NONE type
void test_mask_eval_none_type(bool use_gpu)
{
  struct gkyl_dg_array_mask_inp mask_inp = {
    .use_gpu = use_gpu
  };
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);

  // For NONE type, eval should always return false regardless of index
  bool *mask_cond_ho, *mask_cond;
  mask_cond_ho = gkyl_malloc(sizeof(bool));
  if (use_gpu)
    mask_cond = gkyl_cu_malloc(sizeof(bool));
  else
    mask_cond = gkyl_malloc(sizeof(bool));

  bool mask_cond_ref[] = {false, false, false, false, false};
  int idx[] = {0, 1, 2, 3, 100};
  for (int i=0; i<5; i++) {
    gkyl_dg_array_mask_eval_idx(mask, &idx[i], mask_cond);
    if (use_gpu)
      gkyl_cu_memcpy(mask_cond_ho, mask_cond, sizeof(bool), GKYL_CU_MEMCPY_D2H);
    else
      memcpy(mask_cond_ho, mask_cond, sizeof(bool));

    TEST_CHECK(mask_cond_ho[0] == mask_cond_ref[i]);
    TEST_MSG("Got %d at idx=%d. Expected %d\n",mask_cond_ho[0], idx[i], mask_cond_ref[i]);
  }

  gkyl_free(mask_cond_ho);
  if (use_gpu)
    gkyl_cu_free(mask_cond);
  else
    gkyl_free(mask_cond);

  gkyl_dg_array_mask_release(mask);
}

// Test mask scale_by_cell
void test_mask_scale_by_cell(bool use_gpu)
{
  int shape[] = { 5 };
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 1, shape);

  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS,
    .threshold = 0.5,
    .conf_rng = &range,
    .conf_rng_ext = &range,
    .use_gpu = use_gpu
  };

  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);

  // Set mask values
  // Create an array with alternating 1.0 and -1.0, then use mask_advance to set mask values
  struct gkyl_array *arr = mkarr(use_gpu, 1, range.volume);
  struct gkyl_array *arr_ho = mkarr(false, arr->ncomp, arr->size);
  for (unsigned i = 0; i < range.volume; ++i) {
    double *arr_d = gkyl_array_fetch(arr_ho, i);
    arr_d[0] = (i % 2 == 0) ? 0.0 : 1.0; // 0.0 will be masked (mask=1.0), 1.0 will not (mask=-1.0)
  }
  gkyl_array_copy(arr, arr_ho);
  
  double global_max = 0.0;
  gkyl_array_reduce(&global_max, arr_ho, GKYL_MAX);
  gkyl_dg_array_mask_advance_threshold(mask, global_max);
  gkyl_dg_array_mask_advance(mask, arr);

  // Save original mask values
  const struct gkyl_array *mask_arr = gkyl_dg_array_mask_get_mask(mask);
  struct gkyl_array *mask_orig_ho = mkarr(false, mask_arr->ncomp, mask_arr->size);
  gkyl_array_copy(mask_orig_ho, mask_arr);
  double mask_orig[range.volume];
  for (unsigned i = 0; i < range.volume; ++i) {
    const double *mask_d = gkyl_array_cfetch(mask_orig_ho, i);
    mask_orig[i] = mask_d[0];
  }

  // Create array to multiply (1 component per cell)
  gkyl_array_clear(arr_ho, 0.0);
  for (unsigned i = 0; i < range.volume; ++i) {
    double *arr_d = gkyl_array_fetch(arr_ho, i);
    arr_d[0] = 2.0 + i;
  }

  gkyl_array_copy(arr, arr_ho);

  // Scale mask by array
  gkyl_dg_array_mask_scale_by_cell(mask, arr);

  // Copy mask back to host for verification
  struct gkyl_array *mask_ho = mkarr(false, mask_arr->ncomp, mask_arr->size);
  gkyl_array_copy(mask_ho, mask_arr);

  // Check results: mask should be multiplied by arr values
  for (unsigned cell = 0; cell < range.volume; ++cell) {
    const double *arr_val = gkyl_array_cfetch(arr_ho, cell);
    double expected = mask_orig[cell] * arr_val[0];
    const double *mask_d = gkyl_array_cfetch(mask_ho, cell);
    TEST_CHECK(gkyl_compare(mask_d[0], expected, 1e-14) );
  }

  gkyl_array_release(arr_ho);
  gkyl_array_release(arr);
  gkyl_array_release(mask_orig_ho);
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

// Test mask acquire and release
void test_mask_acquire_release(bool use_gpu)
{
  int shape[] = { 3, 3 };
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);

  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS,
    .threshold = 1.0,
    .conf_rng = &range,
    .conf_rng_ext = &range,
    .use_gpu = use_gpu
  };

  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);

  // Acquire reference
  struct gkyl_dg_array_mask *mask2 = gkyl_dg_array_mask_acquire(mask);
  TEST_CHECK(mask2 == mask);

  // Release both references
  gkyl_dg_array_mask_release(mask2);
  gkyl_dg_array_mask_release(mask);
}

// Test threshold scaling with different dimensions
void test_mask_threshold_scaling(bool use_gpu)
{
  // Test 1D
  int shape1d[] = { 10 };
  struct gkyl_range range1d;
  gkyl_range_init_from_shape(&range1d, 1, shape1d);

  double threshold = 1.0;
  struct gkyl_dg_array_mask_inp mask_inp1d = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS,
    .threshold = threshold,
    .conf_rng = &range1d,
    .conf_rng_ext = &range1d,
    .use_gpu = use_gpu
  };

  struct gkyl_dg_array_mask *mask1d = gkyl_dg_array_mask_new(mask_inp1d);
  double expected1d = threshold * sqrt(2.0);
  TEST_CHECK(gkyl_compare(mask1d->threshold, expected1d, 1e-14) );
  gkyl_dg_array_mask_release(mask1d);

  // Test 2D
  int shape2d[] = { 5, 5 };
  struct gkyl_range range2d;
  gkyl_range_init_from_shape(&range2d, 2, shape2d);

  struct gkyl_dg_array_mask_inp mask_inp2d = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS,
    .threshold = threshold,
    .conf_rng = &range2d,
    .conf_rng_ext = &range2d,
    .use_gpu = use_gpu
  };

  struct gkyl_dg_array_mask *mask2d = gkyl_dg_array_mask_new(mask_inp2d);
  double expected2d = threshold * pow(sqrt(2.0), 2);
  TEST_CHECK(gkyl_compare(mask2d->threshold, expected2d, 1e-14) );
  gkyl_dg_array_mask_release(mask2d);

  // Test 3D
  int shape3d[] = { 3, 3, 3 };
  struct gkyl_range range3d;
  gkyl_range_init_from_shape(&range3d, 3, shape3d);

  struct gkyl_dg_array_mask_inp mask_inp3d = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS,
    .threshold = threshold,
    .conf_rng = &range3d,
    .conf_rng_ext = &range3d,
    .use_gpu = use_gpu
  };

  struct gkyl_dg_array_mask *mask3d = gkyl_dg_array_mask_new(mask_inp3d);
  double expected3d = threshold * pow(sqrt(2.0), 3);
  TEST_CHECK(gkyl_compare(mask3d->threshold, expected3d, 1e-14) );
  gkyl_dg_array_mask_release(mask3d);
}

// Test fractional threshold mask (global)
void test_mask_advance_frac_threshold(bool use_gpu)
{
  int shape[] = { 10 };
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 1, shape);

  double frac_threshold = 0.5; // 50% of max value
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS_FRAC,
    .threshold = frac_threshold,
    .conf_rng = &range,
    .conf_rng_ext = &range,
    .use_gpu = use_gpu
  };

  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);

  // Create test array with known max value
  struct gkyl_array *arr = mkarr(use_gpu, 1, range.volume);
  struct gkyl_array *arr_ho = mkarr(false, arr->ncomp, arr->size);
  double max_val = 10.0;
  for (unsigned i = 0; i < range.volume; ++i) {
    double *arr_d = gkyl_array_fetch(arr_ho, i);
    arr_d[0] = (i + 1) * 1.0; // Values from 1.0 to 10.0
  }
  gkyl_array_copy(arr, arr_ho);
  double global_max = 0.0;
  gkyl_array_reduce(&global_max, arr_ho, GKYL_MAX);
  gkyl_dg_array_mask_advance_threshold(mask, global_max);
  gkyl_dg_array_mask_advance(mask, arr);

  // Copy mask back to host for verification`
  const struct gkyl_array *mask_arr = gkyl_dg_array_mask_get_mask(mask);
  struct gkyl_array *mask_ho = mkarr(false, mask_arr->ncomp, mask_arr->size);
  gkyl_array_copy(mask_ho, mask_arr);

  // Expected threshold is 50% of max (10.0) = 5.0
  double expected_threshold = frac_threshold * max_val;
  for (unsigned i = 0; i < range.volume; ++i) {
    const double *mask_d = gkyl_array_cfetch(mask_ho, i);
    const double *arr_d = gkyl_array_cfetch(arr_ho, i);
    if (fabs(arr_d[0]) < expected_threshold) {
      TEST_CHECK(gkyl_compare(mask_d[0], 1.0, 1e-14) );  // True (masked)
    }
    else {
      TEST_CHECK(gkyl_compare(mask_d[0], -1.0, 1e-14) );  // False (not masked)
    }
  }

  gkyl_array_release(arr_ho);
  gkyl_array_release(arr);
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

// Test fractional threshold mask GREATER_THAN (global)
void test_mask_advance_frac_threshold_greater(bool use_gpu)
{
  int shape[] = { 8 };
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 1, shape);

  double frac_threshold = 0.3; // 30% of max value
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_GREATER_FRAC,
    .threshold = frac_threshold,
    .conf_rng = &range,
    .conf_rng_ext = &range,
    .use_gpu = use_gpu
  };

  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);

  // Create test array with varying values
  struct gkyl_array *arr = mkarr(use_gpu, 1, range.volume);
  struct gkyl_array *arr_ho = mkarr(false, arr->ncomp, arr->size);
  double vals[] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 };
  double max_val = 8.0;
  for (unsigned i = 0; i < range.volume; ++i) {
    double *arr_d = gkyl_array_fetch(arr_ho, i);
    arr_d[0] = vals[i];
  }
  gkyl_array_copy(arr, arr_ho);
  
  double global_max = 0.0;
  gkyl_array_reduce(&global_max, arr_ho, GKYL_MAX);
  gkyl_dg_array_mask_advance_threshold(mask, global_max);
  gkyl_dg_array_mask_advance(mask, arr);

  // Copy mask back to host for verification
  const struct gkyl_array *mask_arr = gkyl_dg_array_mask_get_mask(mask);
  struct gkyl_array *mask_ho = mkarr(false, mask_arr->ncomp, mask_arr->size);
  gkyl_array_copy(mask_ho, mask_arr);

  // Expected threshold is 30% of max (8.0) = 2.4
  double expected_threshold = frac_threshold * max_val;
  for (unsigned i = 0; i < range.volume; ++i) {
    const double *mask_d = gkyl_array_cfetch(mask_ho, i);
    const double *arr_d = gkyl_array_cfetch(arr_ho, i);
    if (fabs(arr_d[0]) > expected_threshold) {
      TEST_CHECK(gkyl_compare(mask_d[0], 1.0, 1e-14) );  // True (masked)
    }
    else {
      TEST_CHECK(gkyl_compare(mask_d[0], -1.0, 1e-14) );  // False (not masked)
    }
  }

  gkyl_array_release(arr_ho);
  gkyl_array_release(arr);
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

// Test spatial fractional threshold mask (velocity-space dependent)
void test_mask_advance_frac_threshold_spatial(bool use_gpu)
{
  // Set up 2x2 configuration space, 3x3 velocity space
  int conf_shape[] = { 2, 2 };
  int vel_shape[] = { 3, 3 };
  int phase_shape[] = { 2, 2, 3, 3 };

  struct gkyl_range conf_rng, vel_rng, phase_rng;
  gkyl_range_init_from_shape(&conf_rng, 2, conf_shape);
  gkyl_range_init_from_shape(&vel_rng, 2, vel_shape);
  gkyl_range_init_from_shape(&phase_rng, 4, phase_shape);

  double frac_threshold = 0.5; // 50% of local max
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS_FRAC_CONF,
    .threshold = frac_threshold,
    .phase_rng = &phase_rng,
    .phase_rng_ext = &phase_rng,
    .conf_rng = &conf_rng,
    .conf_rng_ext = &conf_rng,
    .vel_rng = &vel_rng,
    .use_gpu = use_gpu
  };

  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);

  // Create test array with spatially-varying max values
  // For each config cell, create a different max value in velocity space
  struct gkyl_array *arr = mkarr(use_gpu, 1, phase_rng.volume);
  struct gkyl_array *arr_ho = mkarr(false, arr->ncomp, arr->size);

  struct gkyl_range_iter iter_conf, iter_vel;
  gkyl_range_iter_init(&iter_conf, &conf_rng);

  while (gkyl_range_iter_next(&iter_conf)) {
    long conf_idx = gkyl_range_idx(&conf_rng, iter_conf.idx);
    // Set different max for each config cell: 10.0, 20.0, 30.0, 40.0
    double local_max = (conf_idx + 1) * 10.0;

    gkyl_range_iter_init(&iter_vel, &vel_rng);
    while (gkyl_range_iter_next(&iter_vel)) {
      int pidx[4];
      for (int d = 0; d < 2; d++) {
        pidx[d] = iter_conf.idx[d];
      }
      for (int d = 0; d < 2; d++) {
        pidx[2 + d] = iter_vel.idx[d];
      }

      long vel_idx = gkyl_range_idx(&vel_rng, iter_vel.idx);
      long phase_idx = gkyl_range_idx(&phase_rng, pidx);

      double *arr_d = gkyl_array_fetch(arr_ho, phase_idx);
      // Create gradient in velocity space: max at vel_idx=8, min at vel_idx=0
      arr_d[0] = local_max * (vel_idx + 1) / 9.0;
    }
  }

  gkyl_array_copy(arr, arr_ho);
  
  double global_max = 0.0;
  gkyl_array_reduce(&global_max, arr_ho, GKYL_MAX);
  gkyl_dg_array_mask_advance_threshold(mask, global_max);
  gkyl_dg_array_mask_advance(mask, arr);

  // Copy mask back to host for verification
  const struct gkyl_array *mask_arr = gkyl_dg_array_mask_get_mask(mask);
  struct gkyl_array *mask_ho = mkarr(false, mask_arr->ncomp, mask_arr->size);
  gkyl_array_copy(mask_ho, mask_arr);

  // Verify: each config cell should have threshold = 50% of its local max
  gkyl_range_iter_init(&iter_conf, &conf_rng);
  while (gkyl_range_iter_next(&iter_conf)) {
    long conf_idx = gkyl_range_idx(&conf_rng, iter_conf.idx);
    double local_max = (conf_idx + 1) * 10.0;
    double expected_threshold = frac_threshold * local_max;

    gkyl_range_iter_init(&iter_vel, &vel_rng);
    while (gkyl_range_iter_next(&iter_vel)) {
      int pidx[4];
      for (int d = 0; d < 2; d++) {
        pidx[d] = iter_conf.idx[d];
      }
      for (int d = 0; d < 2; d++) {
        pidx[2 + d] = iter_vel.idx[d];
      }
      long phase_idx = gkyl_range_idx(&phase_rng, pidx);

      const double *mask_d = gkyl_array_cfetch(mask_ho, phase_idx);
      const double *arr_d = gkyl_array_cfetch(arr_ho, phase_idx);

      if (fabs(arr_d[0]) < expected_threshold) {
        TEST_CHECK(gkyl_compare(mask_d[0], 1.0, 1e-14) );  // True (masked)
      }
      else {
        TEST_CHECK(gkyl_compare(mask_d[0], -1.0, 1e-14) );  // False (not masked)
      }
    }
  }

  gkyl_array_release(arr_ho);
  gkyl_array_release(arr);
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

// Test spatial fractional threshold GREATER_THAN
void test_mask_advance_frac_threshold_spatial_greater(bool use_gpu)
{
  // Simpler test: 2 config cells, 4 velocity cells
  int conf_shape[] = { 2 };
  int vel_shape[] = { 4 };
  int phase_shape[] = { 2, 4 };

  struct gkyl_range conf_rng, vel_rng, phase_rng;
  gkyl_range_init_from_shape(&conf_rng, 1, conf_shape);
  gkyl_range_init_from_shape(&vel_rng, 1, vel_shape);
  gkyl_range_init_from_shape(&phase_rng, 2, phase_shape);

  double frac_threshold = 0.6; // 60% of local max
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_GREATER_FRAC_CONF,
    .threshold = frac_threshold,
    .phase_rng = &phase_rng,
    .phase_rng_ext = &phase_rng,
    .conf_rng = &conf_rng,
    .conf_rng_ext = &conf_rng,
    .vel_rng = &vel_rng,
    .use_gpu = use_gpu
  };

  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);

  // Create test array: config cell 0 has max=10, config cell 1 has max=20
  struct gkyl_array *arr = mkarr(use_gpu, 1, phase_rng.volume);
  struct gkyl_array *arr_ho = mkarr(false, arr->ncomp, arr->size);

  for (int ic = 0; ic < 2; ic++) {
    double local_max = (ic + 1) * 10.0;
    for (int iv = 0; iv < 4; iv++) {
      int pidx[2] = { ic, iv };
      long phase_idx = gkyl_range_idx(&phase_rng, pidx);
      double *arr_d = gkyl_array_fetch(arr_ho, phase_idx);
      // Values: 2.5, 5.0, 7.5, 10.0 for ic=0 and 5.0, 10.0, 15.0, 20.0 for ic=1
      arr_d[0] = local_max * (iv + 1) / 4.0;
    }
  }

  gkyl_array_copy(arr, arr_ho);
  
  double global_max = 0.0;
  gkyl_array_reduce(&global_max, arr_ho, GKYL_MAX);
  gkyl_dg_array_mask_advance_threshold(mask, global_max);
  gkyl_dg_array_mask_advance(mask, arr);

  // Copy mask back to host for verification
  const struct gkyl_array *mask_arr = gkyl_dg_array_mask_get_mask(mask);
  struct gkyl_array *mask_ho = mkarr(false, mask_arr->ncomp, mask_arr->size);
  gkyl_array_copy(mask_ho, mask_arr);

  // Verify GREATER_THAN logic with spatially-dependent thresholds
  for (int ic = 0; ic < 2; ic++) {
    double local_max = (ic + 1) * 10.0;
    double expected_threshold = frac_threshold * local_max; // 6.0 and 12.0

    for (int iv = 0; iv < 4; iv++) {
      int pidx[2] = { ic, iv };
      long phase_idx = gkyl_range_idx(&phase_rng, pidx);

      const double *mask_d = gkyl_array_cfetch(mask_ho, phase_idx);
      const double *arr_d = gkyl_array_cfetch(arr_ho, phase_idx);

      if (fabs(arr_d[0]) > expected_threshold) {
        TEST_CHECK(gkyl_compare(mask_d[0], 1.0, 1e-14) );  // True (masked)
      }
      else {
        TEST_CHECK(gkyl_compare(mask_d[0], -1.0, 1e-14) );  // False (not masked)
      }
    }
  }

  gkyl_array_release(arr_ho);
  gkyl_array_release(arr);
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

// Create some tests with a separate extended range

// Test mask advance with value threshold and extended range (ghost cells)
void test_mask_advance_threshold_ext_range(bool use_gpu, int ncell, int nghost_cell)
{
  double lower[] = { 0.0 };
  double upper[] = { 8.0 };
  int cells[] = { ncell };
  int ndim = sizeof(lower) / sizeof(lower[0]);

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  int nghost[GKYL_MAX_DIM];
  for (int d = 0; d < ndim; d++) {
    nghost[d] = nghost_cell;
  }

  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);

  double threshold = 0.5;
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS,
    .threshold = threshold,
    .conf_rng = &local,
    .conf_rng_ext = &local_ext,
    .use_gpu = use_gpu
  };

  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);

  // Verify mask array has extended range size
  TEST_CHECK(mask->mask_arr->size == local_ext.volume);

  // Create test array on extended range with varying values
  struct gkyl_array *arr = mkarr(use_gpu, 1, local_ext.volume);
  struct gkyl_array *arr_ho = mkarr(false, arr->ncomp, arr->size);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local_ext);
  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(&local_ext, iter.idx);
    double *arr_d = gkyl_array_fetch(arr_ho, lidx);
    // Alternate values: ghost cells and interior cells
    arr_d[0] = (iter.idx[0] % 2 == 0) ? 0.0 : 10.0;
  }

  gkyl_array_copy(arr, arr_ho);
  
  double global_max = 0.0;
  gkyl_array_reduce(&global_max, arr_ho, GKYL_MAX);
  gkyl_dg_array_mask_advance_threshold(mask, global_max);
  gkyl_dg_array_mask_advance(mask, arr);

  // Copy mask back to host for verification
  const struct gkyl_array *mask_arr = gkyl_dg_array_mask_get_mask(mask);
  struct gkyl_array *mask_ho = mkarr(false, mask_arr->ncomp, mask_arr->size);
  gkyl_array_copy(mask_ho, mask_arr);

  // Check mask values only on interior (local) range
  double expected_threshold = threshold * sqrt(2.0); // 1D scaling
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(&local_ext, iter.idx);
    const double *mask_d = gkyl_array_cfetch(mask_ho, lidx);
    const double *arr_d = gkyl_array_cfetch(arr_ho, lidx);
    if (fabs(arr_d[0]) < expected_threshold) {
      TEST_CHECK(gkyl_compare(mask_d[0], 1.0, 1e-14) );  // True (masked)
    }
    else {
      TEST_CHECK(gkyl_compare(mask_d[0], -1.0, 1e-14) );  // False (not masked)
    }
  }

  // Test gkyl_dg_array_mask_eval and gkyl_dg_array_mask_eval_idx
  // These are important because eval_idx uses the phase range to compute lidx
  if (!use_gpu) {
    // When the mask is on the GPU, the eval function cannot be called from the host side.
    gkyl_range_iter_init(&iter, &local);
    while (gkyl_range_iter_next(&iter)) {
      long lidx = gkyl_range_idx(&local_ext, iter.idx);
      const double *arr_d = gkyl_array_cfetch(arr_ho, lidx);
      bool expected_mask = fabs(arr_d[0]) < expected_threshold;

      // Test eval_idx with multi-dimensional index
      TEST_CHECK(gkyl_dg_array_mask_eval_idx_ker(mask, iter.idx) == expected_mask);
    }
  }

  gkyl_array_release(arr_ho);
  gkyl_array_release(arr);
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

// Test fractional threshold mask with extended range (ghost cells)
void test_mask_advance_frac_threshold_ext_range(bool use_gpu)
{
  double lower[] = { 0.0 };
  double upper[] = { 10.0 };
  int cells[] = { 10 };
  int ndim = sizeof(lower) / sizeof(lower[0]);

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  int nghost[GKYL_MAX_DIM];
  for (int d = 0; d < ndim; d++) {
    nghost[d] = 1;
  }

  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);

  double frac_threshold = 0.5; // 50% of max value
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS_FRAC,
    .threshold = frac_threshold,
    .conf_rng = &local,
    .conf_rng_ext = &local_ext,
    .use_gpu = use_gpu
  };

  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);

  // Verify mask array has extended range size
  TEST_CHECK(mask->mask_arr->size == local_ext.volume);

  // Create test array with known max value
  // Put max in interior, smaller values in ghost cells
  struct gkyl_array *arr = mkarr(use_gpu, 1, local_ext.volume);
  struct gkyl_array *arr_ho = mkarr(false, arr->ncomp, arr->size);
  double max_val = 10.0;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local_ext);
  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(&local_ext, iter.idx);
    double *arr_d = gkyl_array_fetch(arr_ho, lidx);
    // Ghost cells (idx 0 and 11) get small values, interior cells get gradient
    if (gkyl_range_contains_idx(&local, iter.idx)) {
      // Interior: values from 1.0 to 10.0 based on position in local range
      int local_pos = iter.idx[0] - local.lower[0];
      arr_d[0] = (local_pos + 1) * 1.0;
    }
    else {
      // Ghost cells: small value
      arr_d[0] = 0.5;
    }
  }

  gkyl_array_copy(arr, arr_ho);

  double global_max = 0.0;
  gkyl_array_reduce(&global_max, arr_ho, GKYL_MAX);
  gkyl_dg_array_mask_advance_threshold(mask, global_max);
  gkyl_dg_array_mask_advance(mask, arr);

  // Copy mask back to host for verification
  const struct gkyl_array *mask_arr = gkyl_dg_array_mask_get_mask(mask);
  struct gkyl_array *mask_ho = mkarr(false, mask_arr->ncomp, mask_arr->size);
  gkyl_array_copy(mask_ho, mask_arr);

  // Expected threshold is 50% of max (10.0) = 5.0
  double expected_threshold = frac_threshold * max_val;

  // Check mask values on interior (local) range
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(&local_ext, iter.idx);
    const double *mask_d = gkyl_array_cfetch(mask_ho, lidx);
    const double *arr_d = gkyl_array_cfetch(arr_ho, lidx);
    if (fabs(arr_d[0]) < expected_threshold) {
      TEST_CHECK(gkyl_compare(mask_d[0], 1.0, 1e-14) );  // True (masked)
    }
    else {
      TEST_CHECK(gkyl_compare(mask_d[0], -1.0, 1e-14) );  // False (not masked)
    }
  }

  // Test gkyl_dg_array_mask_eval and gkyl_dg_array_mask_eval_idx
  // These are important because eval_idx uses the phase range to compute lidx
  if (!use_gpu) {
    // When the mask is on the GPU, the eval function cannot be called from the host side.
    gkyl_range_iter_init(&iter, &local);
    while (gkyl_range_iter_next(&iter)) {
      long lidx = gkyl_range_idx(&local_ext, iter.idx);
      const double *arr_d = gkyl_array_cfetch(arr_ho, lidx);
      bool expected_mask = fabs(arr_d[0]) < expected_threshold;

      // Test eval_idx with multi-dimensional index
      TEST_CHECK(gkyl_dg_array_mask_eval_idx_ker(mask, iter.idx) == expected_mask);
    }
  }

  gkyl_array_release(arr_ho);
  gkyl_array_release(arr);
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

// Test spatial fractional threshold mask with extended range (ghost cells)
void test_mask_advance_frac_threshold_spatial_ext_range(bool use_gpu)
{
  // Set up 1D config space with ghost cells, 1D velocity space
  // Using gkyl_create_grid_ranges for both config and phase space
  // to ensure consistent indexing

  double conf_lower[] = { 0.0 };
  double conf_upper[] = { 4.0 };
  int conf_cells[] = { 4 };
  int conf_ndim = 1;

  struct gkyl_rect_grid conf_grid;
  gkyl_rect_grid_init(&conf_grid, conf_ndim, conf_lower, conf_upper, conf_cells);

  int conf_nghost[GKYL_MAX_DIM];
  for (int d = 0; d < conf_ndim; d++) {
    conf_nghost[d] = 1;
  }

  struct gkyl_range conf_local, conf_local_ext;
  gkyl_create_grid_ranges(&conf_grid, conf_nghost, &conf_local_ext, &conf_local);

  // Velocity space using grid ranges as well for consistency
  double vel_lower[] = { 0.0 };
  double vel_upper[] = { 4.0 };
  int vel_cells[] = { 4 };
  int vel_ndim = 1;

  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vel_ndim, vel_lower, vel_upper, vel_cells);

  int vel_nghost[GKYL_MAX_DIM];
  vel_nghost[0] = 0; // No ghost cells in velocity

  struct gkyl_range vel_local, vel_local_ext;
  gkyl_create_grid_ranges(&vel_grid, vel_nghost, &vel_local_ext, &vel_local);

  // Phase space = config x velocity, with ghost cells in config direction only
  double phase_lower[] = { 0.0, 0.0 };
  double phase_upper[] = { 4.0, 4.0 };
  int phase_cells[] = { 4, 4 };
  int phase_ndim = 2;

  struct gkyl_rect_grid phase_grid;
  gkyl_rect_grid_init(&phase_grid, phase_ndim, phase_lower, phase_upper, phase_cells);

  int phase_nghost[GKYL_MAX_DIM];
  phase_nghost[0] = 1; // Ghost cells in config direction
  phase_nghost[1] = 0; // No ghost cells in velocity direction

  struct gkyl_range phase_local, phase_local_ext;
  gkyl_create_grid_ranges(&phase_grid, phase_nghost, &phase_local_ext, &phase_local);

  double frac_threshold = 0.5; // 50% of local max
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS_FRAC_CONF,
    .threshold = frac_threshold,
    .phase_rng = &phase_local,
    .phase_rng_ext = &phase_local_ext,
    .conf_rng = &conf_local,
    .conf_rng_ext = &conf_local_ext,
    .vel_rng = &vel_local,
    .use_gpu = use_gpu
  };

  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);

  // Verify mask array has extended range size
  TEST_CHECK(mask->mask_arr->size == phase_local_ext.volume);

  // Create test array with spatially-varying max values
  struct gkyl_array *arr = mkarr(use_gpu, 1, phase_local_ext.volume);
  struct gkyl_array *arr_ho = mkarr(false, arr->ncomp, arr->size);

  struct gkyl_range_iter iter_conf, iter_vel;
  gkyl_range_iter_init(&iter_conf, &conf_local_ext);
  while (gkyl_range_iter_next(&iter_conf)) {
    // Set different max for each config cell (including ghost cells)
    int conf_pos = iter_conf.idx[0] - conf_local_ext.lower[0];
    double local_max = (conf_pos + 1) * 10.0; // 10, 20, 30, 40, 50, 60 for 6 cells

    gkyl_range_iter_init(&iter_vel, &vel_local);
    while (gkyl_range_iter_next(&iter_vel)) {
      int pidx[2] = { iter_conf.idx[0], iter_vel.idx[0] };
      long phase_idx = gkyl_range_idx(&phase_local_ext, pidx);

      int vel_pos = iter_vel.idx[0] - vel_local.lower[0];
      double *arr_d = gkyl_array_fetch(arr_ho, phase_idx);
      // Create gradient in velocity space
      arr_d[0] = local_max * (vel_pos + 1) / 4.0;
    }
  }

  gkyl_array_copy(arr, arr_ho);
  
  double global_max = 0.0;
  gkyl_array_reduce(&global_max, arr_ho, GKYL_MAX);
  gkyl_dg_array_mask_advance_threshold(mask, global_max);
  gkyl_dg_array_mask_advance(mask, arr);

  // Copy mask back to host for verification
  const struct gkyl_array *mask_arr = gkyl_dg_array_mask_get_mask(mask);
  struct gkyl_array *mask_ho = mkarr(false, mask_arr->ncomp, mask_arr->size);
  gkyl_array_copy(mask_ho, mask_arr);

  // Verify: each config cell should have threshold = 50% of its local max
  // Only check interior cells
  gkyl_range_iter_init(&iter_conf, &conf_local);
  while (gkyl_range_iter_next(&iter_conf)) {
    int conf_pos = iter_conf.idx[0] - conf_local_ext.lower[0];
    double local_max = (conf_pos + 1) * 10.0;
    double expected_threshold = frac_threshold * local_max;

    gkyl_range_iter_init(&iter_vel, &vel_local);
    while (gkyl_range_iter_next(&iter_vel)) {
      int pidx[2] = { iter_conf.idx[0], iter_vel.idx[0] };
      long phase_idx = gkyl_range_idx(&phase_local_ext, pidx);

      const double *mask_d = gkyl_array_cfetch(mask_ho, phase_idx);
      const double *arr_d = gkyl_array_cfetch(arr_ho, phase_idx);

      if (fabs(arr_d[0]) < expected_threshold) {
        TEST_CHECK(gkyl_compare(mask_d[0], 1.0, 1e-14) );  // True (masked)
      }
      else {
        TEST_CHECK(gkyl_compare(mask_d[0], -1.0, 1e-14) );  // False (not masked)
      }
    }
  }

  // Test gkyl_dg_array_mask_eval and gkyl_dg_array_mask_eval_idx
  // These are important because eval_idx uses the phase range to compute lidx
  if (!use_gpu) {
    // When the mask is on the GPU, the eval function cannot be called from the host side.
    gkyl_range_iter_init(&iter_conf, &conf_local);
    while (gkyl_range_iter_next(&iter_conf)) {
      int conf_pos = iter_conf.idx[0] - conf_local_ext.lower[0];
      double local_max = (conf_pos + 1) * 10.0;
      double expected_threshold = frac_threshold * local_max;

      gkyl_range_iter_init(&iter_vel, &vel_local);
      while (gkyl_range_iter_next(&iter_vel)) {
        int pidx[2] = { iter_conf.idx[0], iter_vel.idx[0] };
        long phase_idx = gkyl_range_idx(&phase_local_ext, pidx);

        const double *arr_d = gkyl_array_cfetch(arr_ho, phase_idx);
        bool expected_mask = fabs(arr_d[0]) < expected_threshold;

        // Test eval_idx with multi-dimensional index
        TEST_CHECK(gkyl_dg_array_mask_eval_idx_ker(mask, pidx) == expected_mask);
      }
    }
  }

  gkyl_array_release(arr_ho);
  gkyl_array_release(arr);
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

// CPU test wrappers
void test_mask_new_ho()
{
  test_mask_new(false);
}

void test_mask_none_type_ho()
{
  test_mask_none_type(false);
}

void test_mask_advance_threshold_ho()
{
  test_mask_advance_threshold(false);
}

void test_mask_advance_all_below_ho()
{
  test_mask_advance_all_below(false);
}

void test_mask_advance_all_above_ho()
{
  test_mask_advance_all_above(false);
}

void test_mask_advance_negative_values_ho()
{
  test_mask_advance_negative_values(false);
}

void test_mask_advance_greater_than_threshold_ho()
{
  test_mask_advance_greater_than_threshold(false);
}

void test_mask_advance_greater_than_all_above_ho()
{
  test_mask_advance_greater_than_all_above(false);
}

void test_mask_advance_greater_than_all_below_ho()
{
  test_mask_advance_greater_than_all_below(false);
}

void test_mask_advance_greater_than_negative_values_ho()
{
  test_mask_advance_greater_than_negative_values(false);
}

void test_mask_eval_ho()
{
  test_mask_eval(false);
}

void test_mask_eval_none_type_ho()
{
  test_mask_eval_none_type(false);
}

void test_mask_scale_by_cell_ho()
{
  test_mask_scale_by_cell(false);
}

void test_mask_acquire_release_ho()
{
  test_mask_acquire_release(false);
}

void test_mask_threshold_scaling_ho()
{
  test_mask_threshold_scaling(false);
}

void test_mask_advance_frac_threshold_ho()
{
  test_mask_advance_frac_threshold(false);
}

void test_mask_advance_frac_threshold_greater_ho()
{
  test_mask_advance_frac_threshold_greater(false);
}

void test_mask_advance_frac_threshold_spatial_ho()
{
  test_mask_advance_frac_threshold_spatial(false);
}

void test_mask_advance_frac_threshold_spatial_greater_ho()
{
  test_mask_advance_frac_threshold_spatial_greater(false);
}

void test_mask_advance_threshold_ext_range_ho_1()
{
  test_mask_advance_threshold_ext_range(false, 9, 3);
}

void test_mask_advance_threshold_ext_range_ho_2()
{
  test_mask_advance_threshold_ext_range(false, 9, 2);
}

void test_mask_advance_threshold_ext_range_ho_3()
{
  test_mask_advance_threshold_ext_range(false, 10, 3);
}

void test_mask_advance_threshold_ext_range_ho_4()
{
  test_mask_advance_threshold_ext_range(false, 10, 2);
}

void test_mask_advance_frac_threshold_ext_range_ho()
{
  test_mask_advance_frac_threshold_ext_range(false);
}

void test_mask_advance_frac_threshold_spatial_ext_range_ho()
{
  test_mask_advance_frac_threshold_spatial_ext_range(false);
}

#ifdef GKYL_HAVE_CUDA

// GPU test wrappers
void test_mask_new_dev()
{
  test_mask_new(true);
}

void test_mask_none_type_dev()
{
  test_mask_none_type(true);
}

void test_mask_advance_threshold_dev()
{
  test_mask_advance_threshold(true);
}

void test_mask_advance_all_below_dev()
{
  test_mask_advance_all_below(true);
}

void test_mask_advance_all_above_dev()
{
  test_mask_advance_all_above(true);
}

void test_mask_advance_negative_values_dev()
{
  test_mask_advance_negative_values(true);
}

void test_mask_advance_greater_than_threshold_dev()
{
  test_mask_advance_greater_than_threshold(true);
}

void test_mask_advance_greater_than_all_above_dev()
{
  test_mask_advance_greater_than_all_above(true);
}

void test_mask_advance_greater_than_all_below_dev()
{
  test_mask_advance_greater_than_all_below(true);
}

void test_mask_advance_greater_than_negative_values_dev()
{
  test_mask_advance_greater_than_negative_values(true);
}

void test_mask_eval_none_type_dev()
{
  test_mask_eval_none_type(true);
}

void test_mask_scale_by_cell_dev()
{
  test_mask_scale_by_cell(true);
}

void test_mask_acquire_release_dev()
{
  test_mask_acquire_release(true);
}

void test_mask_threshold_scaling_dev()
{
  test_mask_threshold_scaling(true);
}

void test_mask_advance_frac_threshold_dev()
{
  test_mask_advance_frac_threshold(true);
}

void test_mask_advance_frac_threshold_greater_dev()
{
  test_mask_advance_frac_threshold_greater(true);
}

void test_mask_advance_frac_threshold_spatial_dev()
{
  test_mask_advance_frac_threshold_spatial(true);
}

void test_mask_advance_frac_threshold_spatial_greater_dev()
{
  test_mask_advance_frac_threshold_spatial_greater(true);
}

void test_mask_advance_threshold_ext_range_dev_1()
{
  test_mask_advance_threshold_ext_range(true, 9, 3);
}

void test_mask_advance_threshold_ext_range_dev_2()
{
  test_mask_advance_threshold_ext_range(true, 9, 2);
}

void test_mask_advance_threshold_ext_range_dev_3()
{
  test_mask_advance_threshold_ext_range(true, 10, 3);
}

void test_mask_advance_threshold_ext_range_dev_4()
{
  test_mask_advance_threshold_ext_range(true, 10, 2);
}

void test_mask_advance_frac_threshold_ext_range_dev()
{
  test_mask_advance_frac_threshold_ext_range(true);
}

void test_mask_advance_frac_threshold_spatial_ext_range_dev()
{
  test_mask_advance_frac_threshold_spatial_ext_range(true);
}

#endif

TEST_LIST = {
  { "mask_new", test_mask_new_ho },
  { "mask_none_type", test_mask_none_type_ho },
  { "mask_advance_threshold", test_mask_advance_threshold_ho },
  { "mask_advance_all_below", test_mask_advance_all_below_ho },
  { "mask_advance_all_above", test_mask_advance_all_above_ho },
  { "mask_advance_negative_values", test_mask_advance_negative_values_ho },
  { "mask_advance_greater_than_threshold", test_mask_advance_greater_than_threshold_ho },
  { "mask_advance_greater_than_all_above", test_mask_advance_greater_than_all_above_ho },
  { "mask_advance_greater_than_all_below", test_mask_advance_greater_than_all_below_ho },
  { "mask_advance_greater_than_neg_vals", test_mask_advance_greater_than_negative_values_ho },
  { "mask_eval", test_mask_eval_ho },
  { "mask_eval_none_type", test_mask_eval_none_type_ho },
  { "mask_scale_by_cell", test_mask_scale_by_cell_ho },
  { "mask_acquire_release", test_mask_acquire_release_ho },
  { "mask_threshold_scaling", test_mask_threshold_scaling_ho },
  { "mask_advance_frac_threshold", test_mask_advance_frac_threshold_ho },
  { "mask_advance_frac_threshold_greater", test_mask_advance_frac_threshold_greater_ho },
  { "mask_advance_frac_threshold_spatial", test_mask_advance_frac_threshold_spatial_ho },
  { "mask_advance_frac_threshold_spatial_greater",
    test_mask_advance_frac_threshold_spatial_greater_ho },
  { "mask_advance_threshold_ext_range_1", test_mask_advance_threshold_ext_range_ho_1 },
  { "mask_advance_threshold_ext_range_2", test_mask_advance_threshold_ext_range_ho_2 },
  { "mask_advance_threshold_ext_range_3", test_mask_advance_threshold_ext_range_ho_3 },
  { "mask_advance_threshold_ext_range_4", test_mask_advance_threshold_ext_range_ho_4 },
  { "mask_advance_frac_threshold_ext_range", test_mask_advance_frac_threshold_ext_range_ho },
  { "mask_advance_frac_threshold_spatial_ext_range",
    test_mask_advance_frac_threshold_spatial_ext_range_ho },
#ifdef GKYL_HAVE_CUDA
  { "mask_new_gpu", test_mask_new_dev },
  { "mask_none_type_gpu", test_mask_none_type_dev },
  { "mask_advance_threshold_gpu", test_mask_advance_threshold_dev },
  { "mask_advance_all_below_gpu", test_mask_advance_all_below_dev },
  { "mask_advance_all_above_gpu", test_mask_advance_all_above_dev },
  { "mask_advance_negative_values_gpu", test_mask_advance_negative_values_dev },
  { "mask_advance_greater_than_threshold_gpu", test_mask_advance_greater_than_threshold_dev },
  { "mask_advance_greater_than_all_above_gpu", test_mask_advance_greater_than_all_above_dev },
  { "mask_advance_greater_than_all_below_gpu", test_mask_advance_greater_than_all_below_dev },
  { "mask_advance_greater_than_neg_vals_gpu", test_mask_advance_greater_than_negative_values_dev },
  { "mask_eval_none_type_gpu", test_mask_eval_none_type_dev },
  { "mask_scale_by_cell_gpu", test_mask_scale_by_cell_dev },
  { "mask_acquire_release_gpu", test_mask_acquire_release_dev },
  { "mask_threshold_scaling_gpu", test_mask_threshold_scaling_dev },
  { "mask_advance_frac_threshold_gpu", test_mask_advance_frac_threshold_dev },
  { "mask_advance_frac_threshold_greater_gpu", test_mask_advance_frac_threshold_greater_dev },
  { "mask_advance_frac_threshold_spatial_gpu", test_mask_advance_frac_threshold_spatial_dev },
  { "mask_advance_frac_threshold_spatial_greater_gpu",
    test_mask_advance_frac_threshold_spatial_greater_dev },
  { "mask_advance_threshold_ext_range_1_gpu", test_mask_advance_threshold_ext_range_dev_1 },
  { "mask_advance_threshold_ext_range_2_gpu", test_mask_advance_threshold_ext_range_dev_2 },
  { "mask_advance_threshold_ext_range_3_gpu", test_mask_advance_threshold_ext_range_dev_3 },
  { "mask_advance_threshold_ext_range_4_gpu", test_mask_advance_threshold_ext_range_dev_4 },
  { "mask_advance_frac_threshold_ext_range_gpu", test_mask_advance_frac_threshold_ext_range_dev },
  { "mask_advance_frac_threshold_spatial_ext_range_gpu",
    test_mask_advance_frac_threshold_spatial_ext_range_dev },
#endif
  { NULL, NULL },
};
