#include <math.h>

#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_gpython.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

static void
check_quad_roundtrip(gpython_basis *basis, int expected_num_quad)
{
  TEST_ASSERT(basis != NULL);
  int num_basis = gpython_basis_num_basis(basis);
  int num_quad = gpython_basis_num_quad(basis);
  TEST_ASSERT(num_quad == expected_num_quad);

  double *modal = gkyl_calloc(num_basis, sizeof(double));
  double *quad = gkyl_malloc(num_quad * sizeof(double));
  double *restored = gkyl_malloc(num_basis * sizeof(double));

  // Every mode, including the highest velocity modes, must survive projection.
  for (int mode = 0; mode < num_basis; ++mode) {
    modal[mode] = 1.0;
    for (int node = 0; node < num_quad; ++node) {
      quad[node] = NAN;
    }
    for (int index = 0; index < num_basis; ++index) {
      restored[index] = NAN;
    }
    TEST_CHECK(gpython_basis_modal_to_quad(basis, modal, quad) == 0);
    TEST_CHECK(gpython_basis_quad_to_modal(basis, quad, restored) == 0);
    for (int index = 0; index < num_basis; ++index) {
      TEST_CHECK(gkyl_compare_double(restored[index], modal[index], 1e-12));
    }
    modal[mode] = 0.0;
  }

  gkyl_free(restored);
  gkyl_free(quad);
  gkyl_free(modal);
  gpython_basis_release(basis);
}

void
test_gpython_quad_cartesian()
{
  const char *types[] = {"serendipity", "tensor"};
  TEST_CHECK(gpython_api_version() == GPYTHON_API_VERSION);
  for (int type = 0; type < 2; ++type) {
    for (int order = 1; order <= 2; ++order) {
      int num_quad = 1;
      for (int ndim = 1; ndim <= 3; ++ndim) {
        num_quad *= order + 1;
        check_quad_roundtrip(gpython_basis_new(types[type], ndim, order), num_quad);
      }
    }
  }
}

void
test_gpython_quad_hybrid()
{
  for (int cdim = 1; cdim <= 3; ++cdim) {
    int num_quad = 1 << cdim;
    for (int vdim = 1; vdim <= 3; ++vdim) {
      num_quad *= 3;
      check_quad_roundtrip(gpython_basis_new_hybrid("hybrid", cdim, vdim), num_quad);
    }
    for (int vdim = 1; vdim <= 2; ++vdim) {
      check_quad_roundtrip(
        gpython_basis_new_hybrid("gkhybrid", cdim, vdim), 3 * (1 << (cdim + vdim - 1))
      );
    }
  }
}

void
test_gpython_quad_unavailable()
{
  const char *types[] = {"serendipity", "tensor"};
  for (int type = 0; type < 2; ++type) {
    for (int order = 0; order <= 3; order += 3) {
      gpython_basis *basis = gpython_basis_new(types[type], 1, order);
      TEST_ASSERT(basis != NULL);
      double modal[4] = {0.0}, quad[4] = {0.0};
      TEST_CHECK(gpython_basis_num_quad(basis) == 0);
      TEST_CHECK(gpython_basis_modal_to_quad(basis, modal, quad) == -1);
      TEST_CHECK(gpython_basis_quad_to_modal(basis, quad, modal) == -1);
      gpython_basis_release(basis);
    }
  }
}

static void
eval_average_quadratic(double tm, const double *xn, double *out, void *ctx)
{
  double x = xn[0], y = xn[1];
  out[0] = x * x + x + 3.0 * y + 2.0;
}

static void
eval_average_weight(double tm, const double *xn, double *out, void *ctx)
{
  out[0] = 2.0 + xn[0];
}

void
test_gpython_proj_on_basis_invalid()
{
  const double lower[] = {-1.0, -1.0}, upper[] = {1.0, 1.0};
  const int cells[] = {2, 3}, wrong_cells[] = {2, 4}, empty_cells[] = {0, 3};
  struct gkyl_range range;
  gkyl_range_init_from_shape1(&range, 2, cells);
  gpython_basis *basis = gpython_basis_new("serendipity", range.ndim, 2);
  gpython_array *field = gpython_array_new(gpython_basis_num_basis(basis), range.volume);

  TEST_CHECK(
    gpython_proj_on_basis(
      1, lower, upper, cells, basis, 3, 1, eval_average_quadratic, NULL, 0.0, field
    ) != 0
  );
  TEST_CHECK(
    gpython_proj_on_basis(
      2, lower, upper, cells, basis, -1, 1, eval_average_quadratic, NULL, 0.0, field
    ) != 0
  );
  TEST_CHECK(
    gpython_proj_on_basis(
      2, lower, upper, cells, basis, 3, 0, eval_average_quadratic, NULL, 0.0, field
    ) != 0
  );
  TEST_CHECK(
    gpython_proj_on_basis(
      2, lower, upper, cells, basis, 3, 2, eval_average_quadratic, NULL, 0.0, field
    ) != 0
  );
  TEST_CHECK(
    gpython_proj_on_basis(2, lower, upper, cells, basis, 3, 1, NULL, NULL, 0.0, field) != 0
  );
  TEST_CHECK(
    gpython_proj_on_basis(
      2, lower, upper, wrong_cells, basis, 3, 1, eval_average_quadratic, NULL, 0.0, field
    ) != 0
  );
  TEST_CHECK(
    gpython_proj_on_basis(
      2, lower, upper, empty_cells, basis, 3, 1, eval_average_quadratic, NULL, 0.0, field
    ) != 0
  );
  TEST_CHECK(
    gpython_proj_on_basis(
      2, upper, lower, cells, basis, 3, 1, eval_average_quadratic, NULL, 0.0, field
    ) != 0
  );

  gpython_array_release(field);
  gpython_basis_release(basis);
}

static void
check_average_quadratic(const int *cells)
{
  const double lower[] = {-1.0, -1.0}, upper[] = {1.0, 1.0};
  const int avg_dim[] = {1, 0}, ghost[] = {0, 0};
  struct gkyl_rect_grid grid, target_grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);
  gkyl_rect_grid_init(&target_grid, 1, &lower[1], &upper[1], &cells[1]);
  struct gkyl_range local, local_ext, target_local, target_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  gkyl_create_grid_ranges(&target_grid, ghost, &target_ext, &target_local);

  gpython_basis *basis = gpython_basis_new("serendipity", 2, 2);
  gpython_basis *target_basis = gpython_basis_new("serendipity", 1, 2);
  int num_basis = gpython_basis_num_basis(basis);
  int target_num_basis = gpython_basis_num_basis(target_basis);
  gpython_array *field = gpython_array_new(num_basis, local_ext.volume);
  gpython_array *weight = gpython_array_new(num_basis, local_ext.volume);
  gpython_array *result = gpython_array_new(target_num_basis, target_ext.volume);
  int num_quad = gpython_basis_poly_order(basis) + 1;
  TEST_ASSERT(
    gpython_proj_on_basis(
      grid.ndim, lower, upper, cells, basis, num_quad, 1, eval_average_quadratic, NULL, 0.0, field
    ) == 0
  );
  TEST_ASSERT(
    gpython_proj_on_basis(
      grid.ndim, lower, upper, cells, basis, num_quad, 1, eval_average_weight, NULL, 0.0, weight
    ) == 0
  );

  for (int weighted = 0; weighted < 2; ++weighted) {
    TEST_CHECK(
      gpython_array_average(
        grid.ndim, lower, upper, cells, basis, target_basis, target_grid.ndim, target_grid.cells,
        avg_dim, weighted ? weight : NULL, field, result
      ) == 0
    );
    // Integrating x gives 7/3+3*y, or 5/2+3*y with the weight 2+x.
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &target_local);
    while (gkyl_range_iter_next(&iter)) {
      double center[1];
      gkyl_rect_grid_cell_center(&target_grid, iter.idx, center);
      const double *coeff = gpython_array_cfetch(result, gkyl_range_idx(&target_local, iter.idx));
      double expected[] = {
        sqrt(2.0) * ((weighted ? 2.5 : 7.0 / 3.0) + 3.0 * center[0]),
        sqrt(3.0 / 2.0) * target_grid.dx[0], 0.0
      };
      for (int mode = 0; mode < target_num_basis; ++mode) {
        TEST_CHECK(gkyl_compare_double(coeff[mode], expected[mode], 2e-13));
      }
    }
  }
  gpython_array_release(result);
  gpython_array_release(weight);
  gpython_array_release(field);
  gpython_basis_release(target_basis);
  gpython_basis_release(basis);
}

void
test_gpython_average_quadratic()
{
  check_average_quadratic((int[]){1, 1});
  check_average_quadratic((int[]){2, 3});
}

struct projection_fields_ctx {
  double constant[2];
  double slope[2][2];
};

static void
eval_projection_fields(double tm, const double *xn, double *out, void *ctx)
{
  const struct projection_fields_ctx *params = ctx;
  for (int component = 0; component < 2; ++component) {
    out[component] = params->constant[component] + params->slope[component][0] * xn[0] +
                     params->slope[component][1] * xn[1] + tm;
  }
}

void
test_gpython_projection_multiple_fields()
{
  gpython_basis *basis = gpython_basis_new("serendipity", 2, 1);
  const double lower[] = {-1.0, 0.25}, upper[] = {2.0, 2.25};
  const int cells[] = {2, 3}, ghost[] = {0, 0};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct projection_fields_ctx ctx = {.constant = {1.0, 3.0}, .slope = {{1.0, 2.0}, {-2.0, 4.0}}};
  const double tm = 0.25;
  int num_basis = gpython_basis_num_basis(basis);
  int nfields = sizeof(ctx.constant) / sizeof(ctx.constant[0]);
  gpython_array *field = gpython_array_new(nfields * num_basis, local_ext.volume);

  // Start with nonzero entries to check that projection overwrites every mode.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    double *coeff = gpython_array_fetch(field, gkyl_range_idx(&local, iter.idx));
    for (int mode = 0; mode < nfields * num_basis; ++mode) {
      coeff[mode] = -99.0;
    }
  }
  TEST_ASSERT(
    gpython_proj_on_basis(
      grid.ndim, lower, upper, cells, basis, 0, nfields, eval_projection_fields, &ctx, tm, field
    ) == 0
  );

  // Check the projected coefficients independently before reducing dimensions.
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    double center[2];
    gkyl_rect_grid_cell_center(&grid, iter.idx, center);
    const double *coeff = gpython_array_cfetch(field, gkyl_range_idx(&local, iter.idx));
    for (int component = 0; component < nfields; ++component) {
      double expected[] = {
        2.0 * (ctx.constant[component] + ctx.slope[component][0] * center[0] +
               ctx.slope[component][1] * center[1] + tm),
        ctx.slope[component][0] * grid.dx[0] / sqrt(3.0),
        ctx.slope[component][1] * grid.dx[1] / sqrt(3.0), 0.0
      };
      for (int mode = 0; mode < num_basis; ++mode) {
        TEST_CHECK(gkyl_compare_double(coeff[component * num_basis + mode], expected[mode], 2e-13));
      }
    }
  }
  for (int direction = 0; direction < 2; ++direction) {
    int kept = 1 - direction;
    struct gkyl_rect_grid target_grid;
    gkyl_rect_grid_init(&target_grid, 1, &lower[kept], &upper[kept], &cells[kept]);
    struct gkyl_range target_local, target_ext;
    gkyl_create_grid_ranges(&target_grid, ghost, &target_ext, &target_local);
    gpython_basis *target_basis =
      gpython_basis_new("serendipity", target_grid.ndim, gpython_basis_poly_order(basis));
    int target_num_basis = gpython_basis_num_basis(target_basis);
    for (int point = 0; point < 3; ++point) {
      double coordinate = lower[direction] + 0.5 * point * (upper[direction] - lower[direction]);
      int btype, order, cdim, vdim;
      gpython_array *result = gpython_eval_at_coord_proj(
        basis, 2, 2, lower, upper, cells, 1, &direction, &coordinate, 1, &cells[kept], field,
        &btype, &order, &cdim, &vdim
      );
      TEST_ASSERT(result != NULL);
      TEST_CHECK(order == gpython_basis_poly_order(target_basis));
      TEST_CHECK(gpython_array_ncomp(result) == nfields * target_num_basis);
      TEST_CHECK(gpython_array_size(result) == target_ext.volume);
      gkyl_range_iter_init(&iter, &target_local);
      while (gkyl_range_iter_next(&iter)) {
        double center[1];
        gkyl_rect_grid_cell_center(&target_grid, iter.idx, center);
        const double *coeff = gpython_array_cfetch(result, gkyl_range_idx(&target_local, iter.idx));
        for (int component = 0; component < nfields; ++component) {
          int offset = component * target_num_basis;
          double mean = ctx.constant[component] + ctx.slope[component][direction] * coordinate +
                        ctx.slope[component][kept] * center[0] + tm;
          TEST_CHECK(gkyl_compare_double(coeff[offset], sqrt(2.0) * mean, 2e-13));
          TEST_CHECK(gkyl_compare_double(
            coeff[offset + 1], target_grid.dx[0] * ctx.slope[component][kept] / sqrt(6.0), 2e-13
          ));
        }
      }
      gpython_array_release(result);
    }
    gpython_basis_release(target_basis);
  }
  int directions[] = {0, 1}, cells_tar[] = {1};
  double coordinates[] = {0.7, 1.6};
  int btype, order, cdim, vdim;
  gpython_array *result = gpython_eval_at_coord_proj(
    basis, 2, 2, lower, upper, cells, 2, directions, coordinates, 1, cells_tar, field, &btype,
    &order, &cdim, &vdim
  );
  TEST_ASSERT(result != NULL);
  TEST_CHECK(order == 0 && cdim == 1 && vdim == 0);
  TEST_CHECK(gpython_array_ncomp(result) == nfields);
  const double *coeff = gpython_array_cfetch(result, 0);
  for (int component = 0; component < nfields; ++component) {
    double expected = ctx.constant[component] + ctx.slope[component][0] * coordinates[0] +
                      ctx.slope[component][1] * coordinates[1] + tm;
    TEST_CHECK(gkyl_compare_double(coeff[component], sqrt(2.0) * expected, 2e-13));
  }
  gpython_array_release(result);
  gpython_array_release(field);
  gpython_basis_release(basis);
}

TEST_LIST = {
  {"quad_cartesian", test_gpython_quad_cartesian},
  {"quad_hybrid", test_gpython_quad_hybrid},
  {"quad_unavailable", test_gpython_quad_unavailable},
  {"proj_on_basis_invalid", test_gpython_proj_on_basis_invalid},
  {"average_quadratic", test_gpython_average_quadratic},
  {"projection_multiple_fields", test_gpython_projection_multiple_fields},
  {NULL, NULL}
};
