#include <stdio.h>
#include <stdlib.h>
#include <acutest.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_position_map.h>
#include <gkyl_position_map_priv.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_calc_bmag.h>
#include <float.h>
#include <math.h>

void
test_nonuniform_position_map(
  double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx
)
{
  double poly_order = 2;
  double z = xn[0];
  double left = 0.25;
  double right = 0.75;
  if (z < -left) {
    fout[0] = z;
  } else if (z < right) {
    fout[0] = -pow(z - right, poly_order) / fabs(pow(left - right, poly_order - 1)) + right;
  } else {
    fout[0] = z;
  }
}

void
test_nonuniform_position_map_slope(
  double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx
)
{
  double poly_order = 2;
  double z = xn[0];
  double left = 0.25;
  double right = 0.75;
  if (z < -left) {
    fout[0] = 1.0;
  } else if (z < right) {
    fout[0] =
      -poly_order * pow(z - right, poly_order - 1) / fabs(pow(left - right, poly_order - 1));
  } else {
    fout[0] = 1.0;
  }
}

void
test_identity_position_map(
  double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx
)
{
  fout[0] = xn[0];
}

void
test_nonuniform_position_map_3x(
  double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx
)
{
  double poly_order = 2;
  double left = 0.25;
  double right = 0.75;
  for (int i = 0; i < 3; i++) {
    double z = xn[i];
    if (z < -left) {
      fout[i] = z;
    } else if (z < right) {
      fout[i] = -pow(z - right, poly_order) / fabs(pow(left - right, poly_order - 1)) + right;
    } else {
      fout[i] = z;
    }
  }
}

void
bmag_func(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0];
  double s = 0.6;
  double c = 0.;
  // double B = (4*pow(s*(x-c),2) - 0.3*pow(s*(x-c),4) + 1)*exp(-pow(s*(x-c),2));
  double B = 1 / (1 + 100 * pow(x - M_PI / 2, 2)) + 1 / (1 + 100 * pow(x + M_PI / 2, 2));
  fout[0] = B;
}

void
test_position_map_init_1x_ho()
{
  int cells[] = {32};
  int poly_order = 1;
  double lower[] = {0.0}, upper[] = {1.0};
  int dim = sizeof(lower) / sizeof(lower[0]);
  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);
  // Ranges
  int ghost[] = {1, 1};
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  struct gkyl_position_map_inp pos_map_inp = {
    .maps =
      {test_nonuniform_position_map, test_nonuniform_position_map, test_nonuniform_position_map},
    .ctxs = {NULL, NULL, NULL},
  };

  struct gkyl_position_map *pos_map = gkyl_position_map_new(
    pos_map_inp, grid, localRange, localRange_ext, localRange, localRange_ext, basis
  );

  TEST_CHECK(pos_map->to_optimize == 0);
  TEST_CHECK(pos_map->grid.ndim == 1);
  TEST_CHECK(pos_map->local.ndim == 1);
  TEST_CHECK(pos_map->local_ext.ndim == 1);
  TEST_CHECK(pos_map->basis.ndim == 1);
  TEST_CHECK(pos_map->basis.poly_order == 1);

  gkyl_position_map_release(pos_map);
}

void
test_position_map_init_1x_null_ho()
{
  int cells[] = {8};
  int poly_order = 1;
  double lower[] = {0.0}, upper[] = {1.0};
  int dim = sizeof(lower) / sizeof(lower[0]);
  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);
  // Ranges
  int ghost[] = {1, 1};
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  struct gkyl_position_map_inp pos_map_inp = {};

  struct gkyl_position_map *pos_map = gkyl_position_map_new(
    pos_map_inp, grid, localRange, localRange_ext, localRange, localRange_ext, basis
  );

  TEST_CHECK(pos_map->id == GKYL_PMAP_USER_INPUT);
  for (double i = 0; i < 1; i = i + 0.1) {
    double x[1] = {i};
    double y[1];
    pos_map->maps[0](0.0, x, y, pos_map->ctxs[0]);
    TEST_CHECK(y[0] == x[0]);

    pos_map->maps[1](0.0, x, y, pos_map->ctxs[1]);
    TEST_CHECK(y[0] == x[0]);

    pos_map->maps[2](0.0, x, y, pos_map->ctxs[2]);
    TEST_CHECK(y[0] == x[0]);
  }
  TEST_CHECK(pos_map->to_optimize == 0);
  TEST_CHECK(pos_map->grid.ndim == 1);
  TEST_CHECK(pos_map->local.ndim == 1);
  TEST_CHECK(pos_map->local_ext.ndim == 1);
  TEST_CHECK(pos_map->basis.ndim == 1);
  TEST_CHECK(pos_map->basis.poly_order == 1);

  gkyl_position_map_release(pos_map);
}

void
test_position_map_init_2x_ho()
{
  int cells[] = {8, 8};
  int poly_order = 1;
  double lower[] = {0.0, 0.0}, upper[] = {1.0, 1.0};
  int dim = sizeof(lower) / sizeof(lower[0]);
  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);
  // Ranges
  int ghost[] = {1, 1};
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  struct gkyl_position_map_inp pos_map_inp = {
    .maps =
      {test_nonuniform_position_map, test_nonuniform_position_map, test_nonuniform_position_map},
    .ctxs = {0, 0, 0},
  };

  struct gkyl_position_map *pos_map = gkyl_position_map_new(
    pos_map_inp, grid, localRange, localRange_ext, localRange, localRange_ext, basis
  );

  TEST_CHECK(pos_map->to_optimize == 0);
  TEST_CHECK(pos_map->grid.ndim == 2);
  TEST_CHECK(pos_map->local.ndim == 2);
  TEST_CHECK(pos_map->local_ext.ndim == 2);
  TEST_CHECK(pos_map->basis.ndim == 2);
  TEST_CHECK(pos_map->basis.poly_order == 1);

  gkyl_position_map_release(pos_map);
}

void
test_position_map_init_3x_ho()
{
  int cells[] = {8, 8, 8};
  int poly_order = 1;
  double lower[] = {0.0, 0.0, 0.0}, upper[] = {1.0, 1.0, 1.0};
  int dim = sizeof(lower) / sizeof(lower[0]);
  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);
  // Ranges
  int ghost[] = {1, 1, 1};
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  struct gkyl_position_map_inp pos_map_inp = {
    .maps =
      {test_nonuniform_position_map, test_nonuniform_position_map, test_nonuniform_position_map},
    .ctxs = {0, 0, 0},
  };

  struct gkyl_position_map *pos_map = gkyl_position_map_new(
    pos_map_inp, grid, localRange, localRange_ext, localRange, localRange_ext, basis
  );

  TEST_CHECK(pos_map->to_optimize == 0);
  TEST_CHECK(pos_map->grid.ndim == 3);
  TEST_CHECK(pos_map->local.ndim == 3);
  TEST_CHECK(pos_map->local_ext.ndim == 3);
  TEST_CHECK(pos_map->basis.ndim == 3);
  TEST_CHECK(pos_map->basis.poly_order == 1);

  gkyl_position_map_release(pos_map);
}

void
test_position_map_set_ho()
{
  int cells[] = {8, 8, 8};
  int poly_order = 1;
  double lower[] = {0.0, 0.0, 0.0}, upper[] = {1.0, 1.0, 1.0};
  int dim = sizeof(lower) / sizeof(lower[0]);
  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);
  // Ranges
  int ghost[] = {1, 1, 1};
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  struct gkyl_position_map_inp pos_map_inp = {
    .maps =
      {test_nonuniform_position_map, test_nonuniform_position_map, test_nonuniform_position_map},
    .ctxs = {0, 0, 0},
  };

  struct gkyl_position_map *pos_map = gkyl_position_map_new(
    pos_map_inp, grid, localRange, localRange_ext, localRange, localRange_ext, basis
  );

  struct gkyl_array *pmap_arr_set =
    gkyl_array_new(GKYL_DOUBLE, 3 * pos_map->basis.num_basis, pos_map->local_ext.volume);
  gkyl_array_clear(pmap_arr_set, 1.0);

  gkyl_position_map_set_mc2nu(pos_map, pmap_arr_set);

  double *pos_map_i = pos_map->mc2nu->data;
  for (unsigned i = 0; i < pos_map->mc2nu->size; ++i) {
    TEST_CHECK(gkyl_compare(pos_map_i[i], 1.0, 1e-14));
  }

  gkyl_array_release(pmap_arr_set);
  gkyl_position_map_release(pos_map);
}

void
test_position_map_eval_mc2nu_ho()
{
  int cells[] = {8, 8, 8};
  int poly_order = 2;
  double lower[] = {0.0, 0.0, 0.0}, upper[] = {1.0, 1.0, 1.0};
  int dim = sizeof(lower) / sizeof(lower[0]);
  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);
  // Ranges
  int ghost[] = {1, 1, 1};
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  struct gkyl_position_map_inp pos_map_inp = {
    .maps =
      {test_nonuniform_position_map, test_nonuniform_position_map, test_nonuniform_position_map},
    .ctxs = {0, 0, 0},
  };

  struct gkyl_position_map *pos_map = gkyl_position_map_new(
    pos_map_inp, grid, localRange, localRange_ext, localRange, localRange_ext, basis
  );

  struct gkyl_array *pmap_arr_set =
    gkyl_array_new(GKYL_DOUBLE, 3 * pos_map->basis.num_basis, pos_map->local_ext.volume);

  gkyl_proj_on_basis *projDistf =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 3, test_nonuniform_position_map_3x, 0);
  gkyl_proj_on_basis_advance(projDistf, 0.0, &localRange, pmap_arr_set);
  gkyl_proj_on_basis_release(projDistf);

  gkyl_position_map_set_mc2nu(pos_map, pmap_arr_set);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 5; k++) {
        double x[3] = {i / 10.0, j / 10.0, k / 10.0};
        double x_fa[3];
        gkyl_position_map_eval_mc2nu(pos_map, x, x_fa);
        double x_analytic[3];
        test_nonuniform_position_map(0.0, &x[0], &x_analytic[0], 0);
        test_nonuniform_position_map(0.0, &x[1], &x_analytic[1], 0);
        test_nonuniform_position_map(0.0, &x[2], &x_analytic[2], 0);
        for (int d = 0; d < 3; ++d) {
          TEST_CHECK(gkyl_compare(x_fa[d], x_analytic[d], 1e-12));
        }
      }
    }
  }

  gkyl_array_release(pmap_arr_set);
  gkyl_position_map_release(pos_map);
}

void
test_position_map_slope_ho()
{
  int cells[] = {8, 8, 8};
  int poly_order = 2;
  double lower[] = {0.0, 0.0, 0.0}, upper[] = {1.0, 1.0, 1.0};
  int dim = sizeof(lower) / sizeof(lower[0]);
  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);
  // Ranges
  int ghost[] = {1, 1, 1};
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  struct gkyl_position_map_inp pos_map_inp = {
    .maps =
      {test_nonuniform_position_map, test_nonuniform_position_map, test_nonuniform_position_map},
    .ctxs = {0, 0, 0},
  };

  struct gkyl_position_map *pos_map = gkyl_position_map_new(
    pos_map_inp, grid, localRange, localRange_ext, localRange, localRange_ext, basis
  );

  struct gkyl_array *pmap_arr_set =
    gkyl_array_new(GKYL_DOUBLE, 3 * pos_map->basis.num_basis, pos_map->local_ext.volume);

  gkyl_proj_on_basis *projDistf =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 3, test_nonuniform_position_map_3x, 0);
  gkyl_proj_on_basis_advance(projDistf, 0.0, &localRange, pmap_arr_set);
  gkyl_proj_on_basis_release(projDistf);

  gkyl_position_map_set_mc2nu(pos_map, pmap_arr_set);

  for (int i = 0; i < 8; i++) {
    for (int j = 0; j < 8; j++) {
      for (int k = 0; k < 8; k++) {
        double x[3] = {i / 8.0, j / 8.0, k / 8.0};
        if (x[0] == 0.25 || x[0] == 0.75) {
          continue;
        }
        if (x[1] == 0.25 || x[1] == 0.75) {
          continue;
        }
        if (x[2] == 0.25 || x[2] == 0.75) {
          continue;
        }
        double x_analytic[3];
        test_nonuniform_position_map_slope(0.0, &x[0], &x_analytic[0], 0);
        test_nonuniform_position_map_slope(0.0, &x[1], &x_analytic[1], 0);
        test_nonuniform_position_map_slope(0.0, &x[2], &x_analytic[2], 0);
        double slope[3], finite_diff_slope;
        slope[0] = gkyl_position_map_slope(pos_map, 0, x[0], 1e-6, grid.lower[0], grid.upper[0]);
        slope[1] = gkyl_position_map_slope(pos_map, 1, x[1], 1e-6, grid.lower[1], grid.upper[1]);
        slope[2] = gkyl_position_map_slope(pos_map, 2, x[2], 1e-6, grid.lower[2], grid.upper[2]);
        for (int d = 0; d < 3; ++d) {
          TEST_CHECK(gkyl_compare(slope[d], x_analytic[d], 1e-6));
        }
      }
    }
  }
  gkyl_array_release(pmap_arr_set);
  gkyl_position_map_release(pos_map);
}

void
test_position_polynomial_map_optimize_1x_ho()
{
  int cells[] = {64};
  int poly_order = 1;
  double lower[] = {-M_PI + 1e-2}, upper[] = {M_PI - 1e-2};
  int dim = sizeof(lower) / sizeof(lower[0]);
  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);
  // Ranges
  int ghost[] = {1, 1};
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  struct gkyl_position_map_inp pos_map_inp = {
    .id = GKYL_PMAP_CONSTANT_DB_POLYNOMIAL,
    .map_strength = 1.0,
  };

  struct gkyl_position_map *pos_map = gkyl_position_map_new(
    pos_map_inp, grid, localRange, localRange_ext, localRange, localRange_ext, basis
  );

  // Project bmag_func onto bmag_global
  struct gkyl_array *bmag_global =
    gkyl_array_new(GKYL_DOUBLE, basis.num_basis, localRange_ext.volume);
  gkyl_proj_on_basis *projB =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, bmag_func, 0);
  gkyl_proj_on_basis_advance(projB, 0.0, &localRange, bmag_global);
  gkyl_proj_on_basis_release(projB);

  struct gkyl_rect_grid grid3D;
  double lower3D[] = {0.4, -0.1, lower[0]}, upper3D[] = {0.6, 0.1, upper[0]};
  int cells3D[] = {1, 1, cells[0]};
  gkyl_rect_grid_init(&grid3D, 3, lower3D, upper3D, cells3D);
  int ghost3D[] = {1, 1, 1};
  struct gkyl_range localRange3D, localRange3D_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid3D, ghost3D, &localRange3D_ext, &localRange3D);

  gkyl_position_map_optimize(pos_map, grid3D, localRange3D);
  gkyl_position_map_set_bmag(pos_map, NULL, bmag_global);
  gkyl_position_map_optimize(pos_map, grid3D, localRange3D);

  TEST_CHECK(pos_map->to_optimize == true);
  TEST_CHECK(gkyl_compare(pos_map->constB_ctx->theta_throat, 1.565796, 1e-6));
  TEST_CHECK(gkyl_compare(pos_map->constB_ctx->Bmag_throat, 1.093613, 1e-6));
  TEST_CHECK(gkyl_compare(pos_map->constB_ctx->psi, 0.5, 1e-6));
  TEST_CHECK(gkyl_compare(pos_map->constB_ctx->alpha, 0.0, 1e-6));
  TEST_CHECK(gkyl_compare(pos_map->constB_ctx->map_strength, 1.0, 1e-6));
  TEST_CHECK(pos_map->constB_ctx->map_order_center == 2);
  TEST_CHECK(pos_map->constB_ctx->map_order_expander == 3);
  TEST_CHECK(pos_map->constB_ctx->N_theta_boundaries == 65);

  gkyl_position_map_release(pos_map);
  gkyl_array_release(bmag_global);
}

void
test_position_map_numeric_optimize_1x_ho()
{
  int cells[] = {64};
  int poly_order = 1;
  double lower[] = {-M_PI + 1e-2}, upper[] = {M_PI - 1e-2};
  int dim = sizeof(lower) / sizeof(lower[0]);
  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);
  // Ranges
  int ghost[] = {1, 1};
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  struct gkyl_position_map_inp pos_map_inp = {
    .id = GKYL_PMAP_CONSTANT_DB_NUMERIC,
    .map_strength = 1.0,
  };

  struct gkyl_position_map *pos_map = gkyl_position_map_new(
    pos_map_inp, grid, localRange, localRange_ext, localRange, localRange_ext, basis
  );

  // Project bmag_func onto bmag_global
  struct gkyl_array *bmag_global =
    gkyl_array_new(GKYL_DOUBLE, basis.num_basis, localRange_ext.volume);
  gkyl_proj_on_basis *projB =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, bmag_func, 0);
  gkyl_proj_on_basis_advance(projB, 0.0, &localRange, bmag_global);
  gkyl_proj_on_basis_release(projB);

  struct gkyl_rect_grid grid3D;
  double lower3D[] = {0.4, -0.1, lower[0]}, upper3D[] = {0.6, 0.1, upper[0]};
  int cells3D[] = {1, 1, cells[0]};
  gkyl_rect_grid_init(&grid3D, 3, lower3D, upper3D, cells3D);
  int ghost3D[] = {1, 1, 1};
  struct gkyl_range localRange3D, localRange3D_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid3D, ghost3D, &localRange3D_ext, &localRange3D);

  gkyl_position_map_optimize(pos_map, grid3D, localRange3D);
  gkyl_position_map_set_bmag(pos_map, NULL, bmag_global);
  gkyl_position_map_optimize(pos_map, grid3D, localRange3D);

  double theta_extrema_analytic[5] = {lower[0], lower[0] / 2, 0.0, upper[0] / 2, upper[0]};

  TEST_CHECK(pos_map->constB_ctx->num_extrema == 5);
  TEST_CHECK(gkyl_compare(pos_map->constB_ctx->theta_extrema[0], theta_extrema_analytic[0], 1e-15));
  TEST_CHECK(gkyl_compare(pos_map->constB_ctx->theta_extrema[1], theta_extrema_analytic[1], 1e-15));
  TEST_CHECK(gkyl_compare(pos_map->constB_ctx->theta_extrema[2], theta_extrema_analytic[2], 1e-15));
  TEST_CHECK(gkyl_compare(pos_map->constB_ctx->theta_extrema[3], theta_extrema_analytic[3], 1e-15));
  TEST_CHECK(gkyl_compare(pos_map->constB_ctx->theta_extrema[4], theta_extrema_analytic[4], 1e-15));

  gkyl_position_map_release(pos_map);
  gkyl_array_release(bmag_global);
}

void
test_position_map_numeric_calculate_1x_ho()
{
  int cells[] = {64};
  int poly_order = 1;
  double lower[] = {-M_PI + 1e-2}, upper[] = {M_PI - 1e-2};
  int dim = sizeof(lower) / sizeof(lower[0]);
  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);
  // Ranges
  int ghost[] = {1, 1};
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  struct gkyl_position_map_inp pos_map_inp = {
    .id = GKYL_PMAP_CONSTANT_DB_NUMERIC,
    .map_strength = 1.0,
  };

  struct gkyl_position_map *pos_map = gkyl_position_map_new(
    pos_map_inp, grid, localRange, localRange_ext, localRange, localRange_ext, basis
  );

  // Project bmag_func onto bmag_global
  struct gkyl_array *bmag_global =
    gkyl_array_new(GKYL_DOUBLE, basis.num_basis, localRange_ext.volume);
  gkyl_proj_on_basis *projB =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, bmag_func, 0);
  gkyl_proj_on_basis_advance(projB, 0.0, &localRange, bmag_global);
  gkyl_proj_on_basis_release(projB);

  struct gkyl_rect_grid grid3D;
  double lower3D[] = {0.4, -0.1, lower[0]}, upper3D[] = {0.6, 0.1, upper[0]};
  int cells3D[] = {1, 1, cells[0]};
  gkyl_rect_grid_init(&grid3D, 3, lower3D, upper3D, cells3D);
  int ghost3D[] = {1, 1, 1};
  struct gkyl_range localRange3D, localRange3D_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid3D, ghost3D, &localRange3D_ext, &localRange3D);

  gkyl_position_map_optimize(pos_map, grid3D, localRange3D);
  gkyl_position_map_set_bmag(pos_map, NULL, bmag_global);
  gkyl_position_map_optimize(pos_map, grid3D, localRange3D);

  double theta_map = 1.0;
  pos_map->maps[2](0.0, &theta_map, &theta_map, pos_map->ctxs[2]);
  TEST_CHECK(gkyl_compare(theta_map, 1.505924, 1e-5));

  gkyl_position_map_release(pos_map);
  gkyl_array_release(bmag_global);
}

void
test_null_map_is_identity()
{
  struct gkyl_position_map *pmap = gkyl_position_map_null_new();

  TEST_CHECK(pmap != NULL);
  TEST_CHECK(pmap->id == GKYL_PMAP_USER_INPUT);
  TEST_CHECK(pmap->to_optimize == false);

  // Each of the three maps must act as the identity, and each derivative as 1.
  for (double z = -1.0; z <= 1.0; z += 0.25) {
    for (int i = 0; i < 3; i++) {
      double x[1] = {z}, y[1] = {0.0};
      pmap->maps[i](0.0, x, y, pmap->ctxs[i]);
      TEST_CHECK(gkyl_compare(y[0], z, 1e-15));

      double dy[1] = {0.0};
      pmap->map_derivs[i](0.0, x, dy, pmap->ctxs[i]);
      TEST_CHECK(gkyl_compare(dy[0], 1.0, 1e-15));
    }
  }

  gkyl_position_map_release(pmap);
}

static void
nonuniform_map_1d(double t, const double *xn, double *fout, void *ctx)
{
  // A simple smooth monotone map z -> z + 0.1*sin(z) (identity-like near 0).
  fout[0] = xn[0] + 0.1 * sin(xn[0]);
}

void
test_inew_constructor()
{
  int cells[] = {16};
  int poly_order = 1;
  double lower[] = {-1.0}, upper[] = {1.0};
  int dim = 1;

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);
  int ghost[] = {1};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  struct gkyl_position_map_inew_inp inp = {
    .pmap_info =
      {
        .maps = {nonuniform_map_1d, nonuniform_map_1d, nonuniform_map_1d},
        .ctxs = {NULL, NULL, NULL},
      },
    .grid = grid,
    .local = local,
    .local_ext = local_ext,
    .global = local,
    .global_ext = local_ext,
    .basis = basis,
  };

  struct gkyl_position_map *pmap = gkyl_position_map_inew(inp);

  TEST_CHECK(pmap != NULL);
  TEST_CHECK(pmap->grid.ndim == 1);
  TEST_CHECK(pmap->basis.poly_order == 1);
  TEST_CHECK(pmap->id == GKYL_PMAP_USER_INPUT);

  // The user map must be wired in: check a couple of values.
  double x[1] = {0.3}, y[1];
  pmap->maps[0](0.0, x, y, pmap->ctxs[0]);
  TEST_CHECK(gkyl_compare(y[0], 0.3 + 0.1 * sin(0.3), 1e-14));

  gkyl_position_map_release(pmap);
}

void
test_acquire_refcount()
{
  struct gkyl_position_map *pmap = gkyl_position_map_null_new();

  struct gkyl_position_map *pmap2 = gkyl_position_map_acquire(pmap);
  TEST_CHECK(pmap2 == pmap);

  // Drop the first reference; the object must survive for use through pmap2.
  gkyl_position_map_release(pmap);

  double x[1] = {0.5}, y[1];
  pmap2->maps[0](0.0, x, y, pmap2->ctxs[0]);
  TEST_CHECK(gkyl_compare(y[0], 0.5, 1e-15));

  gkyl_position_map_release(pmap2);
}

void
test_set_compression()
{
  struct gkyl_position_map *pmap = gkyl_position_map_null_new();

  // With both compression factors zero the maps stay at the identity backups,
  // but the geometric parameters must be recorded.
  pmap->xpt_ctx->compression_factor = 0.0;
  pmap->xpt_ctx->radial_compression_factor = 0.0;

  double zcut = 1.25, zcenter = 0.1, w = 0.05, psisep = 0.7;
  gkyl_position_map_set_compression(pmap, zcut, zcenter, w, psisep);

  TEST_CHECK(gkyl_compare(pmap->xpt_ctx->zcut, zcut, 1e-15));
  TEST_CHECK(gkyl_compare(pmap->xpt_ctx->zcenter, zcenter, 1e-15));
  TEST_CHECK(gkyl_compare(pmap->xpt_ctx->w, w, 1e-15));
  TEST_CHECK(gkyl_compare(pmap->xpt_ctx->psisep, psisep, 1e-15));

  // Maps 0 and 1 fall back to identity backups when compression is disabled.
  double x[1] = {0.4}, y[1];
  pmap->maps[0](0.0, x, y, pmap->ctxs[0]);
  TEST_CHECK(gkyl_compare(y[0], 0.4, 1e-15));
  pmap->maps[1](0.0, x, y, pmap->ctxs[1]);
  TEST_CHECK(gkyl_compare(y[0], 0.4, 1e-15));

  gkyl_position_map_release(pmap);
}

struct bounded_map_ctx {
  int calls, outside;
  bool cubic;
};

static void bounded_map(double t, const double *xn, double *out, void *ctx)
{
  struct bounded_map_ctx *map = ctx;
  map->calls++;
  if (xn[0] < 0.0 || xn[0] > 1.0) {
    map->outside++;
  }
  out[0] = 0.1 + 0.7 * xn[0] + 0.2 * xn[0] * xn[0] * (map->cubic ? xn[0] : 1.0);
}

static void bounded_map_derivative(double t, const double *xn, double *out, void *ctx)
{
  struct bounded_map_ctx *map = ctx;
  out[0] = 0.7 + (map->cubic ? 0.6 * xn[0] * xn[0] : 0.4 * xn[0]);
}

static void test_slope_global_bounds(void)
{
  // The map index remains three-dimensional on 1D and 2D simulation grids.
  for (int cdim = 1; cdim <= 3; ++cdim) {
    struct gkyl_rect_grid grid;
    gkyl_rect_grid_init(&grid, cdim, (double[]){0, 0, 0}, (double[]){1, 1, 1}, (int[]){8, 8, 8});
    struct gkyl_range local, ext;
    gkyl_create_grid_ranges(&grid, (int[]){1, 1, 1}, &ext, &local);
    struct gkyl_basis basis;
    gkyl_cart_modal_serendip(&basis, cdim, 1);
    struct bounded_map_ctx map = {0};
    struct gkyl_position_map *pmap = gkyl_position_map_new(
      (struct gkyl_position_map_inp
      ){.maps = {bounded_map, bounded_map, bounded_map}, .ctxs = {&map, &map, &map}},
      grid, local, ext, local, ext, basis
    );
    double points[] = {0.0, 1e-14, 0.025, 0.37, 0.975, 1.0 - 1e-14, 1.0};
    for (int dim = 0; dim < 3; ++dim) {
      for (int point = 0; point < sizeof(points) / sizeof(points[0]); ++point) {
        for (int big_step = 0; big_step < 2; ++big_step) {
          double actual =
            gkyl_position_map_slope(pmap, dim, points[point], big_step ? 2.0 : 0.1, 0.0, 1.0);
          TEST_CHECK(gkyl_compare(actual, 0.7 + 0.4 * points[point], 2e-13));
        }
      }
    }
    TEST_CHECK(map.outside == 0);
    // Second-order convergence for a cubic, including one-sided endpoints.
    map.cubic = true;
    for (int point = 0; point < 3; ++point) {
      double x = point == 0 ? 0.0 : point == 1 ? 0.37 : 1.0;
      double exact;
      bounded_map_derivative(0, &x, &exact, &map);
      double previous = 0.0;
      for (int level = 0; level < 5; ++level) {
        double actual = gkyl_position_map_slope(pmap, 2, x, 0.1 / pow(2, level), 0.0, 1.0);
        double error = fabs(actual - exact);
        if (level) {
          TEST_CHECK(gkyl_compare(previous / error, 4.0, 1e-8));
        }
        previous = error;
      }
    }
    // Analytic derivatives must bypass the mapping entirely.
    pmap->use_map_derivs = true;
    pmap->map_derivs[2] = bounded_map_derivative;
    map.calls = 0;
    double x = 0.37, exact;
    bounded_map_derivative(0, &x, &exact, &map);
    TEST_CHECK(gkyl_position_map_slope(pmap, 2, x, 0.1, 0, 1) == exact);
    TEST_CHECK(map.calls == 0);
    gkyl_position_map_release(pmap);
  }
}

static void test_compression_parameters(void)
{
  struct gkyl_position_map *pmap = gkyl_position_map_null_new();
  pmap->xpt_ctx->compression_factor = 0.6;
  pmap->xpt_ctx->radial_compression_factor = 0.25;
  gkyl_position_map_set_compression(pmap, 1.0, 0.0, 0.5, 0.5);
  TEST_CHECK(gkyl_compare(gkyl_position_map_slope(pmap, 0, 0.5, 0.01, 0, 1), 0.25, 1e-14));
  TEST_CHECK(gkyl_compare(gkyl_position_map_slope(pmap, 2, 1.0, 0.01, -1, 1), 0.6, 1e-14));
  double point = 0.75, mapped;
  pmap->maps[0](0, &point, &mapped, pmap->ctxs[0]);
  TEST_CHECK(gkyl_compare(mapped, point - 0.75 * 0.5 / M_PI, 1e-14));
  pmap->xpt_ctx->compression_factor = pmap->xpt_ctx->radial_compression_factor = 1.0;
  gkyl_position_map_set_compression(pmap, 1.0, 0.0, 0.5, 0.5);
  for (int dim = 0; dim < 3; ++dim) {
    pmap->maps[dim](0, &point, &mapped, pmap->ctxs[dim]);
    TEST_CHECK(gkyl_compare(mapped, point, 1e-14));
    TEST_CHECK(gkyl_compare(gkyl_position_map_slope(pmap, dim, point, 0.01, 0, 1), 1.0, 1e-14));
  }
  // Disabling compression restores the original map and its derivative.
  struct bounded_map_ctx map = {0};
  pmap->xpt_ctx->compression_factor = pmap->xpt_ctx->radial_compression_factor = 0.0;
  pmap->xpt_ctx->maps_backup[0] = bounded_map;
  pmap->xpt_ctx->map_derivs_backup[0] = 0; // exercise numerical fallback
  pmap->xpt_ctx->ctxs_backup[0] = &map;
  gkyl_position_map_set_compression(pmap, 1.0, 0.0, 0.5, 0.5);
  TEST_CHECK(
    gkyl_compare(gkyl_position_map_slope(pmap, 0, point, 0.01, 0, 1), 0.7 + 0.4 * point, 1e-13)
  );
  gkyl_position_map_release(pmap);
}

static void asymmetric_bmag(double t, const double *xn, double *out, void *ctx)
{
  int profile = *(int *)ctx;
  if (profile == 0) {
    out[0] = xn[0] < 0.25 ? 2.0 - xn[0] : 1.75 + 2.0 * (xn[0] - 0.25);
  } else {
    // Include a constant field with roundoff-sized oscillations.
    out[0] = 2.0 + (profile == 2 ? 8.0 * DBL_EPSILON * (cos(4.0 * M_PI * xn[0]) + xn[0]) : 0.0);
  }
}

static void test_numeric_degenerate_and_strength(void)
{
  for (int profile = 0; profile < 3; ++profile) {
    struct gkyl_rect_grid grid, grid3;
    gkyl_rect_grid_init(&grid, 1, (double[]){-1}, (double[]){1}, (int[]){16});
    gkyl_rect_grid_init(&grid3, 3, (double[]){-2, -1, -1}, (double[]){-1, 1, 1}, (int[]){1, 1, 16});
    struct gkyl_range local, ext, global3, ext3;
    gkyl_create_grid_ranges(&grid, (int[]){1}, &ext, &local);
    gkyl_create_grid_ranges(&grid3, (int[]){1, 1, 1}, &ext3, &global3);
    struct gkyl_basis basis;
    gkyl_cart_modal_serendip(&basis, 1, 1);
    struct gkyl_position_map *pmap = gkyl_position_map_new(
      (struct gkyl_position_map_inp){.id = GKYL_PMAP_CONSTANT_DB_NUMERIC, .map_strength = 1.0},
      grid, local, ext, local, ext, basis
    );
    struct gkyl_array *bmag = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, ext.volume);
    bool uniform = profile != 0;
    struct gkyl_proj_on_basis *proj =
      gkyl_proj_on_basis_new(&grid, &basis, 2, 1, asymmetric_bmag, &profile);
    gkyl_proj_on_basis_advance(proj, 0.0, &local, bmag);
    gkyl_proj_on_basis_release(proj);
    gkyl_position_map_optimize(pmap, grid3, global3);
    gkyl_position_map_set_bmag(pmap, NULL, bmag);
    gkyl_position_map_optimize(pmap, grid3, global3);
    TEST_CHECK(gkyl_compare(pmap->constB_ctx->psi, -1.5, 1e-14));
    if (uniform) {
      TEST_CHECK(pmap->constB_ctx->num_extrema == 2);
      TEST_CHECK(pmap->constB_ctx->theta_extrema[0] == -1.0);
      TEST_CHECK(pmap->constB_ctx->theta_extrema[1] == 1.0);
      TEST_CHECK(pmap->constB_ctx->dB_cell == 0.0);
    }
    for (int strength = 0; strength < 3; ++strength) {
      pmap->constB_ctx->map_strength = 0.5 * strength;
      for (int point = 0; point <= 32; ++point) {
        // Include the exact computational location of the asymmetric minimum.
        double x = point == 16 ? -1.0 + 2.0 * 1.25 / 2.75 : -1.0 + point / 16.0;
        // Invert the cumulative |dB| of the piecewise-linear field analytically.
        double change = 2.75 * (x + 1.0) / 2.0;
        double full = change <= 1.25 ? -1.0 + change : 0.25 + (change - 1.25) / 2.0;
        double actual;
        pmap->maps[2](0, &x, &actual, pmap->ctxs[2]);
        double expected = uniform ? x : (1.0 - 0.5 * strength) * x + 0.5 * strength * full;
        TEST_CHECK(gkyl_compare(actual, expected, 2e-12));
        TEST_MSG(
          "profile %d strength %g x %g actual %g expected %g", profile, 0.5 * strength, x, actual,
          expected
        );
      }
    }
    if (!uniform) {
      pmap->constB_ctx->enable_maximum_slope_limits_at_min_B = true;
      pmap->constB_ctx->enable_maximum_slope_limits_at_max_B = true;
      pmap->constB_ctx->maximum_slope_at_min_B = 2.0;
      pmap->constB_ctx->maximum_slope_at_max_B = 2.0;
      double x_min = -1.0 + 2.0 * 1.25 / 2.75;
      for (int strength = 1; strength <= 2; ++strength) {
        pmap->constB_ctx->map_strength = strength * 0.5;
        double left = x_min - 1e-8, right = x_min + 1e-8, mapped_left, mapped_right;
        pmap->maps[2](0, &left, &mapped_left, pmap->ctxs[2]);
        pmap->maps[2](0, &right, &mapped_right, pmap->ctxs[2]);
        TEST_CHECK(mapped_right >= mapped_left);
        TEST_CHECK(mapped_right - mapped_left < 4e-8);
        double previous = -1.0;
        for (int i = 1; i <= 256; ++i) {
          double x = -1.0 + i / 128.0, mapped;
          pmap->maps[2](0, &x, &mapped, pmap->ctxs[2]);
          TEST_CHECK(mapped > previous);
          TEST_CHECK(mapped - previous <= 2.0 / 128.0 + 1e-12);
          previous = mapped;
        }
      }
    }
    // Replacing B must release the previous optimization storage.
    gkyl_position_map_set_bmag(pmap, NULL, bmag);
    gkyl_position_map_optimize(pmap, grid3, global3);
    if (uniform) {
      TEST_CHECK(pmap->constB_ctx->num_extrema == 2);
      TEST_CHECK(pmap->constB_ctx->dB_cell == 0.0);
    }
    gkyl_array_release(bmag);
    gkyl_position_map_release(pmap);
  }
}

TEST_LIST = {
  {"numeric_degenerate_and_strength", test_numeric_degenerate_and_strength},
  {"slope_global_bounds", test_slope_global_bounds},
  {"compression_parameters", test_compression_parameters},
  {"test_position_map_init_1x_ho", test_position_map_init_1x_ho},
  {"test_position_map_init_1x_null_ho", test_position_map_init_1x_null_ho},
  {"test_position_map_init_2x_ho", test_position_map_init_2x_ho},
  {"test_position_map_init_3x_ho", test_position_map_init_3x_ho},
  {"test_position_map_set_ho", test_position_map_set_ho},
  {"test_position_map_eval_mc2nu_ho", test_position_map_eval_mc2nu_ho},
  {"test_position_map_slope_ho", test_position_map_slope_ho},
  {"test_position_polynomial_map_optimize_1x_ho", test_position_polynomial_map_optimize_1x_ho},
  {"test_position_map_numeric_optimize_1x_ho", test_position_map_numeric_optimize_1x_ho},
  {"test_position_map_numeric_calculate_1x_ho", test_position_map_numeric_calculate_1x_ho},
  {"test_position_map_null_map_is_identity", test_null_map_is_identity},
  {"test_position_map_inew_constructor", test_inew_constructor},
  {"test_position_map_acquire_refcount", test_acquire_refcount},
  {"test_position_map_set_compression", test_set_compression},
  {NULL, NULL}
};
