#include <acutest.h>

#include <math.h>

#include <gkyl_gyrokinetic_priv.h>

struct wall_profile_ctx {
  double offset;
  double time_slope;
  double space_slope;
};

static void
wall_profile(double t, const double *xn, double *fout, void *ctx)
{
  const struct wall_profile_ctx *profile = ctx;
  fout[0] = profile->offset + profile->time_slope*t + profile->space_slope*xn[0];
}

static gkyl_gyrokinetic_app
make_app(void)
{
  gkyl_gyrokinetic_app app = {
    .cdim = 1,
    .use_gpu = false,
  };
  double lower[] = { -1.0 }, upper[] = { 1.0 };
  int cells[] = { 4 }, ghost[] = { 1 };

  gkyl_rect_grid_init(&app.grid, app.cdim, lower, upper, cells);
  gkyl_cart_modal_serendip(&app.basis, app.cdim, 1);
  gkyl_create_grid_ranges(&app.grid, ghost, &app.local_ext, &app.local);

  return app;
}

static void
check_profile(const gkyl_gyrokinetic_app *app, const struct gkyl_array *phi,
  const struct wall_profile_ctx *profile, double tm)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local_ext);
  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(&app->local_ext, iter.idx);
    const double *phi_c = gkyl_array_cfetch(phi, lidx);
    double xc[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(&app->grid, iter.idx, xc);

    const double eta[] = { -1.0, 0.0, 1.0 };
    for (int i=0; i<3; ++i) {
      double x = xc[0] + 0.5*app->grid.dx[0]*eta[i];
      double expected = profile->offset + profile->time_slope*tm + profile->space_slope*x;
      TEST_CHECK(fabs(app->basis.eval_expand(&eta[i], phi_c)-expected) < 1e-12);
    }
  }
}

static void
test_species_wall_profiles(void)
{
  // A linear profile lies exactly in the P1 basis. Compare its DG expansion
  // at both cell endpoints and the center with the closed-form profile, first
  // at t=0 and then after advancing each wall independently to t=0.625.
  gkyl_gyrokinetic_app app = make_app();
  struct wall_profile_ctx lower_ctx = {
    .offset = 1.5,
    .time_slope = 2.0,
    .space_slope = -0.25,
  };
  struct wall_profile_ctx upper_ctx = {
    .offset = -0.5,
    .time_slope = -1.0,
    .space_slope = 0.75,
  };
  struct gk_species species = {
    .lower_bc[0] = {
      .type = GKYL_BC_GK_SPECIES_SHEATH,
      .aux_profile = wall_profile,
      .aux_ctx = &lower_ctx,
    },
    .upper_bc[0] = {
      .type = GKYL_BC_GK_SPECIES_SHEATH,
      .aux_profile = wall_profile,
      .aux_ctx = &upper_ctx,
    },
  };

  gk_species_phi_wall_init(&app, &species.lower_bc[0], &species.phi_wall_lo);
  gk_species_phi_wall_init(&app, &species.upper_bc[0], &species.phi_wall_up);
  TEST_ASSERT(species.phi_wall_lo.phi != 0);
  TEST_ASSERT(species.phi_wall_up.phi != 0);
  TEST_ASSERT(species.phi_wall_lo.projector != 0);
  TEST_ASSERT(species.phi_wall_up.projector != 0);
  TEST_ASSERT(species.phi_wall_lo.advance_func != 0);
  TEST_ASSERT(species.phi_wall_up.advance_func != 0);
  check_profile(&app, species.phi_wall_lo.phi, &lower_ctx, 0.0);
  check_profile(&app, species.phi_wall_up.phi, &upper_ctx, 0.0);

  double tm = 0.625;
  gk_species_phi_wall_advance(&app, &species.phi_wall_lo, tm);
  check_profile(&app, species.phi_wall_lo.phi, &lower_ctx, tm);
  check_profile(&app, species.phi_wall_up.phi, &upper_ctx, 0.0);
  gk_species_phi_wall_advance(&app, &species.phi_wall_up, tm);
  check_profile(&app, species.phi_wall_up.phi, &upper_ctx, tm);

  gk_species_phi_wall_release(&app, &species.phi_wall_lo);
  gk_species_phi_wall_release(&app, &species.phi_wall_up);
}

static void
test_grounded_and_non_sheath_walls(void)
{
  // A sheath with no auxiliary profile represents a grounded wall, so every
  // DG coefficient is exactly zero. A non-sheath boundary needs no wall-
  // potential array or projector at all.
  gkyl_gyrokinetic_app app = make_app();
  struct gk_species species = {
    .lower_bc[0] = {
      .type = GKYL_BC_GK_SPECIES_SHEATH,
    },
    .upper_bc[0] = {
      .type = GKYL_BC_GK_SPECIES_REFLECT,
    },
  };

  gk_species_phi_wall_init(&app, &species.lower_bc[0], &species.phi_wall_lo);
  gk_species_phi_wall_init(&app, &species.upper_bc[0], &species.phi_wall_up);
  TEST_ASSERT(species.phi_wall_lo.phi != 0);
  TEST_CHECK(species.phi_wall_lo.phi_host == species.phi_wall_lo.phi);
  TEST_CHECK(species.phi_wall_lo.projector == 0);
  TEST_ASSERT(species.phi_wall_lo.advance_func != 0);
  TEST_CHECK(species.phi_wall_up.phi == 0);
  TEST_CHECK(species.phi_wall_up.phi_host == 0);
  TEST_CHECK(species.phi_wall_up.projector == 0);
  TEST_ASSERT(species.phi_wall_up.advance_func != 0);

  gk_species_phi_wall_advance(&app, &species.phi_wall_lo, 1.0);
  gk_species_phi_wall_advance(&app, &species.phi_wall_up, 1.0);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app.local_ext);
  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(&app.local_ext, iter.idx);
    const double *phi_c = gkyl_array_cfetch(species.phi_wall_lo.phi, lidx);
    for (int k=0; k<app.basis.num_basis; ++k)
      TEST_CHECK(phi_c[k] == 0.0);
  }

  gk_species_phi_wall_release(&app, &species.phi_wall_lo);
  gk_species_phi_wall_release(&app, &species.phi_wall_up);
}

TEST_LIST = {
  { "test_species_wall_profiles", test_species_wall_profiles },
  { "test_grounded_and_non_sheath_walls", test_grounded_and_non_sheath_walls },
  { 0 },
};
