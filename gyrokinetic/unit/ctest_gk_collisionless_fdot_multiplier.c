#include <acutest.h>

#include <gkyl_gyrokinetic_priv.h>

static void
mapc2p(double t, const double *xc, double *xp, void *ctx)
{
  for (int d = 0; d < 3; ++d) {
    xp[d] = xc[d];
  }
}

static void
bfield(double t, const double *xc, double *out, void *ctx)
{
  out[0] = out[1] = 0.0;
  out[2] = 1.0 + 3.0 * pow(sin(xc[2]), 2);
}

static void
initial_dist(double t, const double *xn, double *out, void *ctx)
{
  double bmag = 1.0 + 3.0 * pow(sin(xn[0]), 2);
  out[0] = (1.0 + 0.2 * cos(2.0 * xn[0])) * (1.0 + 0.1 * pow(xn[1], 4)) *
           exp(-0.5 * xn[1] * xn[1] - 2.0 * xn[2] * bmag);
}

static void
source_dist(double t, const double *xn, double *out, void *ctx)
{
  out[0] = 0.01 * exp(-0.5 * xn[1] * xn[1] - xn[2]);
}

static void
collision_frequency(double t, const double *xn, double *out, void *ctx)
{
  out[0] = 0.03;
}

static struct gkyl_gyrokinetic_collisionless
collisionless_input(enum gkyl_gk_collisionless_type type, bool enable)
{
  return (struct gkyl_gyrokinetic_collisionless){
    .type = type,
    .scale_factor = 0.01,
    .time_rate_multiplier =
      {
        .num_multipliers = enable ? 1 : 0,
        .multiplier[0] =
          {
            .type = GKYL_GK_FDOT_MULTIPLIER_FIXED_FACTOR_TIMES_OMEGA_MAX,
            .cellwise_const = true,
            .cfl_factor_times_omega_max = 0.1,
            .write_diagnostics = true,
          },
      },
  };
}

static gkyl_gyrokinetic_app *
new_app(bool use_gpu, enum gkyl_gk_collisionless_type type, bool collisions)
{
  struct gkyl_gk inp = {
    .name = "ctest_gk_collisionless_fdot_multiplier",
    .cdim = 1,
    .lower = {-M_PI / 2.0},
    .upper = {M_PI / 2.0},
    .cells = {8},
    .poly_order = 1,
    .basis_type = GKYL_BASIS_MODAL_SERENDIPITY,
    .geometry = {.geometry_id = GKYL_GEOMETRY_MAPC2P, .mapc2p = mapc2p, .bfield_func = bfield},
    .num_periodic_dir = 1,
    .periodic_dirs = {0},
    .num_species = 1,
    .species = {{
      .name = "ion",
      .charge = 1.0,
      .mass = 1.0,
      .vdim = 2,
      .lower = {-6.0, 0.0},
      .upper = {6.0, 6.0},
      .cells = {12, 8},
      .polarization_density = 1.0,
      .projection = {.proj_id = GKYL_PROJ_FUNC, .func = initial_dist},
      .collisions =
        {.collision_id = collisions ? GKYL_LBO_COLLISIONS : 0, .self_nu = collision_frequency},
      .source =
        {
          .source_id = GKYL_PROJ_SOURCE,
          .num_sources = 1,
          .projection[0] = {.proj_id = GKYL_PROJ_FUNC, .func = source_dist},
        },
    }},
    .field =
      {
        .gkfield_id = GKYL_GK_FIELD_BOLTZMANN,
        .electron_mass = 1.0,
        .electron_charge = -1.0,
        .electron_temp = 1.0,
        .zero_init_field = true,
        .is_static = true,
      },
    .parallelism = {.use_gpu = use_gpu, .cuts = {1}},
  };
  inp.species[0].collisionless = collisionless_input(type, true);
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(false, use_gpu, stderr);
  inp.parallelism.comm = comm;
  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&inp);
  gkyl_gyrokinetic_comms_release(comm);
  gkyl_gyrokinetic_app_apply_ic(app, 0.0);
  return app;
}

static struct gkyl_array *
copy_to_host(const struct gkyl_array *array)
{
  struct gkyl_array *host = gkyl_array_new(GKYL_DOUBLE, array->ncomp, array->size);
  gkyl_array_copy(host, array);
  return host;
}

static double
compute_rhs(gkyl_gyrokinetic_app *app)
{
  struct gk_species *species = &app->species[0];
  const struct gkyl_array *fin[] = {species->f};
  struct gkyl_array *fout[] = {species->f1};
  struct gkyl_array **bflux[] = {species->bflux.f1};
  struct gkyl_update_status status = {.success = true};
  gyrokinetic_rhs(app, 0.0, DBL_MAX, fin, fout, bflux, NULL, NULL, NULL, &status);
  return status.dt_suggested;
}

// Compare the full app RHS against independently assembled collisionless and other terms.
// The latter include a nonzero source and, optionally, LBO collisions.
static void
check_rhs(
  gkyl_gyrokinetic_app *app, const struct gkyl_array *rhs_base, const struct gkyl_array *cfl_base,
  const struct gkyl_array *rhs_cls, const struct gkyl_array *cfl_cls, double multiplier_scale,
  bool loss_cone, double dt
)
{
  struct gk_species *species = &app->species[0];
  struct gkyl_array *rhs = copy_to_host(species->f1);
  struct gkyl_array *cfl = copy_to_host(species->cflrate);
  struct gkyl_array *loss = loss_cone ? copy_to_host(species->fdot_mult.multiplier) : NULL;
  struct gkyl_array *screen =
    multiplier_scale >= 0.0 ? copy_to_host(species->collisionless.fdot_mult.multiplier) : NULL;
  double omega_max = 0.0;
  gkyl_array_reduce_range(&omega_max, cfl_cls, GKYL_MAX, &species->local);

  double expected_max_cfl = 0.0, max_rhs_error = 0.0, max_rhs = 0.0;
  double max_cfl_error = 0.0, max_other_rhs = 0.0;
  int trapped_cells = 0, passing_cells = 0, screened_cells = 0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &species->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&species->local, iter.idx);
    const double *base = gkyl_array_cfetch(rhs_base, loc);
    const double *cls = gkyl_array_cfetch(rhs_cls, loc);
    const double *actual = gkyl_array_cfetch(rhs, loc);
    double cls_cfl = ((const double *)gkyl_array_cfetch(cfl_cls, loc))[0];
    double base_cfl = ((const double *)gkyl_array_cfetch(cfl_base, loc))[0];
    double actual_cfl = ((const double *)gkyl_array_cfetch(cfl, loc))[0];
    double mask = loss ? ((const double *)gkyl_array_cfetch(loss, loc))[0] : 1.0;
    double mult = multiplier_scale < 0.0 ?
                    1.0 :
                    multiplier_scale * (cls_cfl > 0.0 ? fmin(1.0, 0.1 * omega_max / cls_cfl) : 1.0);
    if (screen) {
      TEST_CHECK(gkyl_compare(((const double *)gkyl_array_cfetch(screen, loc))[0], mult, 2e-12));
    }
    screened_cells += mult < 1.0;
    trapped_cells += mask == 0.0;
    passing_cells += mask == 1.0;
    for (int k = 0; k < species->basis.num_basis; ++k) {
      double other = base[k] - cls[k];
      double expected = mask * (mult * cls[k] + other);
      max_rhs_error = fmax(max_rhs_error, fabs(actual[k] - expected));
      max_rhs = fmax(max_rhs, fabs(expected));
      max_other_rhs = fmax(max_other_rhs, fabs(other));
    }
    double expected_cfl = mask * (mult * cls_cfl + base_cfl - cls_cfl);
    max_cfl_error = fmax(max_cfl_error, fabs(actual_cfl - expected_cfl));
    expected_max_cfl = fmax(expected_max_cfl, expected_cfl);
  }
  TEST_CHECK(max_rhs_error < 2e-12 * max_rhs);
  TEST_MSG("RHS error %.17g, norm %.17g", max_rhs_error, max_rhs);
  TEST_CHECK(max_other_rhs > 1e-8);
  TEST_CHECK(max_cfl_error < 2e-12 * expected_max_cfl);
  TEST_CHECK(gkyl_compare(dt, app->cfl / expected_max_cfl, 2e-12));
  TEST_MSG("dt %.17g, expected %.17g", dt, app->cfl / expected_max_cfl);
  if (multiplier_scale >= 0.0 && omega_max > 0.0) {
    TEST_CHECK(screened_cells > 0);
  }
  if (loss_cone) {
    TEST_CHECK(trapped_cells > 0);
    TEST_CHECK(passing_cells > 0);
    gkyl_array_release(loss);
  }
  gkyl_array_release(rhs);
  gkyl_array_release(cfl);
  if (screen) {
    gkyl_array_release(screen);
  }
}

static void
run_collisionless_multiplier(bool use_gpu, enum gkyl_gk_collisionless_type type, bool collisions)
{
  gkyl_gyrokinetic_app *app = new_app(use_gpu, type, collisions);
  struct gk_species *species = &app->species[0];
  double dt_initial = compute_rhs(app);
  struct gkyl_array *rhs_initial = copy_to_host(species->f1);

  struct gkyl_gyrokinetic_collisionless cls_inp = collisionless_input(type, false);
  gkyl_gyrokinetic_app_reset_species_collisionless(app, 0.0, "ion", cls_inp);
  double dt_base = compute_rhs(app);
  struct gkyl_array *rhs_base = copy_to_host(species->f1);
  struct gkyl_array *cfl_base = copy_to_host(species->cflrate);

  gkyl_array_clear(species->f1, 0.0);
  gkyl_array_clear(species->cflrate, 0.0);
  gk_species_collisionless_rhs(app, species, &species->collisionless, species->f, species->f1);
  struct gkyl_array *rhs_cls = copy_to_host(species->f1);
  struct gkyl_array *cfl_cls = copy_to_host(species->cflrate);

  // Enable through the public reset API, as at the start of an OAP.
  cls_inp = collisionless_input(type, true);
  gkyl_gyrokinetic_app_reset_species_collisionless(app, 0.0, "ion", cls_inp);
  double dt_screened = compute_rhs(app);
  check_rhs(app, rhs_base, cfl_base, rhs_cls, cfl_cls, 1.0, false, dt_screened);
  TEST_CHECK(gkyl_compare(dt_screened, dt_initial, 2e-12));
  struct gkyl_array *rhs_reset = copy_to_host(species->f1);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &species->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&species->local, iter.idx);
    const double *initial = gkyl_array_cfetch(rhs_initial, loc);
    const double *reset = gkyl_array_cfetch(rhs_reset, loc);
    for (int k = 0; k < species->basis.num_basis; ++k) {
      TEST_CHECK(gkyl_compare(initial[k], reset[k], 2e-12));
    }
  }
  if (type == GKYL_GK_COLLISIONLESS_NONE) {
    TEST_CHECK(gkyl_compare(dt_screened, dt_base, 2e-12));
  } else if (collisions) {
    TEST_CHECK(dt_screened > dt_base);
    TEST_CHECK(dt_screened < 10.0 * dt_base);
  } else {
    TEST_CHECK(gkyl_compare(dt_screened, 10.0 * dt_base, 2e-12));
  }

  // Multiple collisionless components compose without rescaling collisions or sources.
  cls_inp.time_rate_multiplier.num_multipliers = 2;
  cls_inp.time_rate_multiplier.multiplier[1] = (struct gkyl_gyrokinetic_fdot_multiplier_comp){
    .type = GKYL_GK_FDOT_MULTIPLIER_CONSTANT,
    .cellwise_const = true,
    .time_dilation_scale_const = 0.5,
  };
  gkyl_gyrokinetic_app_reset_species_collisionless(app, 0.0, "ion", cls_inp);
  double dt_chain = compute_rhs(app);
  check_rhs(app, rhs_base, cfl_base, rhs_cls, cfl_cls, 0.5, false, dt_chain);

  // Reusing a slot as NONE must discard the previous constant (also used by omega_H).
  cls_inp.time_rate_multiplier.multiplier[1].type = GKYL_GK_FDOT_MULTIPLIER_NONE;
  gkyl_gyrokinetic_app_reset_species_collisionless(app, 0.0, "ion", cls_inp);
  double dt_none = compute_rhs(app);
  check_rhs(app, rhs_base, cfl_base, rhs_cls, cfl_cls, 1.0, false, dt_none);
  TEST_CHECK(
    gk_fdot_multiplier_get_time_dilation_scale_const(app, &species->collisionless.fdot_mult) == 1.0
  );
  cls_inp.time_rate_multiplier.multiplier[1].type = GKYL_GK_FDOT_MULTIPLIER_CONSTANT;
  gkyl_gyrokinetic_app_reset_species_collisionless(app, 0.0, "ion", cls_inp);

  // A species-wide loss cone still screens the total RHS during OAP.
  struct gkyl_gyrokinetic_fdot_multiplier loss_inp = {
    .num_multipliers = 1,
    .multiplier[0] =
      {.type = GKYL_GK_FDOT_MULTIPLIER_LOSS_CONE, .cellwise_const = true, .write_diagnostics = true
      },
  };
  gkyl_gyrokinetic_app_reset_species_fdot_multiplier(app, 0.0, "ion", loss_inp);
  double dt_loss = compute_rhs(app);
  check_rhs(app, rhs_base, cfl_base, rhs_cls, cfl_cls, 0.5, true, dt_loss);

  // Diagnostics from the two chains have distinct names and independent enable flags.
  gk_species_fdot_multiplier_write(app, species, &species->fdot_mult, 0.0, 0);
  gk_species_collisionless_write_diags(app, species, &species->collisionless, 0.0, 0);
  const char *diag_files[] = {
    "ctest_gk_collisionless_fdot_multiplier-ion_fdot_multiplier_0.gkyl",
    "ctest_gk_collisionless_fdot_multiplier-ion_collisionless_fdot_multiplier_0.gkyl"
  };
  const struct gk_fdot_multiplier *diag_mult[] = {
    &species->fdot_mult, &species->collisionless.fdot_mult
  };
  for (int i = 0; i < 2; ++i) {
    struct gkyl_rect_grid grid;
    struct gkyl_array *written = gkyl_array_new(GKYL_DOUBLE, 1, species->local_ext.volume);
    TEST_CHECK(
      gkyl_grid_sub_array_read(&grid, &species->local, written, diag_files[i]) ==
      GKYL_ARRAY_RIO_SUCCESS
    );
    struct gkyl_array *expected = copy_to_host(diag_mult[i]->multiplier);
    gkyl_range_iter_init(&iter, &species->local);
    while (gkyl_range_iter_next(&iter)) {
      long loc = gkyl_range_idx(&species->local, iter.idx);
      TEST_CHECK(
        ((const double *)gkyl_array_cfetch(written, loc))[0] ==
        ((const double *)gkyl_array_cfetch(expected, loc))[0]
      );
    }
    gkyl_array_release(written);
    gkyl_array_release(expected);
    remove(diag_files[i]);
  }

  // Disabling the collisionless chain must leave the total loss-cone chain intact.
  cls_inp.time_rate_multiplier.num_multipliers = 0;
  gkyl_gyrokinetic_app_reset_species_collisionless(app, 0.0, "ion", cls_inp);
  double dt_loss_only = compute_rhs(app);
  check_rhs(app, rhs_base, cfl_base, rhs_cls, cfl_cls, -1.0, true, dt_loss_only);

  // FDP reset restores the unscreened RHS and CFL, and repeated resets are safe.
  loss_inp.num_multipliers = 0;
  gkyl_gyrokinetic_app_reset_species_fdot_multiplier(app, 0.0, "ion", loss_inp);
  gkyl_gyrokinetic_app_reset_species_collisionless(app, 0.0, "ion", cls_inp);
  double dt_restored = compute_rhs(app);
  check_rhs(app, rhs_base, cfl_base, rhs_cls, cfl_cls, -1.0, false, dt_restored);
  TEST_CHECK(gkyl_compare(dt_restored, dt_base, 2e-12));

  gkyl_array_release(rhs_initial);
  gkyl_array_release(rhs_reset);
  gkyl_array_release(rhs_base);
  gkyl_array_release(cfl_base);
  gkyl_array_release(rhs_cls);
  gkyl_array_release(cfl_cls);
  gkyl_gyrokinetic_app_release(app);
}

static void
test_es(void)
{
  run_collisionless_multiplier(false, GKYL_GK_COLLISIONLESS_ES, false);
}

static void
test_es_lbo(void)
{
  run_collisionless_multiplier(false, GKYL_GK_COLLISIONLESS_ES, true);
}

static void
test_zero_cfl(void)
{
  run_collisionless_multiplier(false, GKYL_GK_COLLISIONLESS_NONE, true);
}

#ifdef GKYL_HAVE_CUDA
static void
test_es_cu(void)
{
  run_collisionless_multiplier(true, GKYL_GK_COLLISIONLESS_ES, false);
}

static void
test_es_lbo_cu(void)
{
  run_collisionless_multiplier(true, GKYL_GK_COLLISIONLESS_ES, true);
}

static void
test_zero_cfl_cu(void)
{
  run_collisionless_multiplier(true, GKYL_GK_COLLISIONLESS_NONE, true);
}
#endif

TEST_LIST = {
  {"collisionless_fdot_es", test_es},
  {"collisionless_fdot_es_lbo", test_es_lbo},
  {"collisionless_fdot_zero_cfl", test_zero_cfl},
#ifdef GKYL_HAVE_CUDA
  {"collisionless_fdot_es_cu", test_es_cu},
  {"collisionless_fdot_es_lbo_cu", test_es_lbo_cu},
  {"collisionless_fdot_zero_cfl_cu", test_zero_cfl_cu},
#endif
  {NULL, NULL}
};
