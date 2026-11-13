// Verify that LBO exchange is invariant under configuration-coordinate rescaling.
#include <gkyl_gyrokinetic_priv.h>
#include <math.h>
#include <stdio.h>

struct test_ctx {
  double scale, curvature, density, flow, temperature;
};

static void
mapc2p(double t, const double *zc, double *xp, void *ctx)
{
  struct test_ctx *tc = ctx;
  xp[0] = zc[0];
  xp[1] = zc[1];
  xp[2] = tc->scale * (zc[2] + tc->curvature * zc[2] * zc[2]);
}

static void
bfield(double t, const double *zc, double *out, void *ctx)
{
  out[0] = out[1] = 0.0;
  out[2] = 1.0;
}

static void
density(double t, const double *xn, double *out, void *ctx)
{
  out[0] = ((struct test_ctx *)ctx)->density;
}

static void
flow(double t, const double *xn, double *out, void *ctx)
{
  out[0] = ((struct test_ctx *)ctx)->flow;
}

static void
temperature(double t, const double *xn, double *out, void *ctx)
{
  out[0] = ((struct test_ctx *)ctx)->temperature;
}

static void
frequency(double t, const double *xn, double *out, void *ctx)
{
  out[0] = 0.1;
}

static gkyl_gyrokinetic_app *
make_app(struct test_ctx *ctx, bool prescribed)
{
  struct gkyl_gk inp = {
    .name = "gk_lbo_jacobian",
    .cdim = 1,
    .lower = {-0.5},
    .upper = {0.5},
    .cells = {2},
    .poly_order = 1,
    .basis_type = GKYL_BASIS_MODAL_SERENDIPITY,
    .geometry =
      {
        .geometry_id = GKYL_GEOMETRY_MAPC2P,
        .mapc2p = mapc2p,
        .c2p_ctx = ctx,
        .bfield_func = bfield,
      },
    .num_periodic_dir = 1,
    .periodic_dirs = {0},
    .num_species = 2,
    .field = {.zero_init_field = true, .is_static = true},
  };
  for (int s = 0; s < 2; ++s) {
    inp.species[s] = (struct gkyl_gyrokinetic_species){
      .mass = s == 0 ? 1.0 : 4.0,
      .charge = s == 0 ? -1.0 : 1.0,
      .vdim = 2,
      .lower = {-10.0, 0.0},
      .upper = {10.0, 40.0},
      .cells = {32, 24},
      .polarization_density = 1.0,
      .projection =
        {
          .proj_id = GKYL_PROJ_BIMAXWELLIAN,
          .density = density,
          .ctx_density = &ctx[s],
          .upar = flow,
          .ctx_upar = &ctx[s],
          .temppar = temperature,
          .ctx_temppar = &ctx[s],
          .tempperp = temperature,
          .ctx_tempperp = &ctx[s],
          .correct_all_moms = true,
        },
      .collisionless = {.type = GKYL_GK_COLLISIONLESS_ES_NO_BY},
      .collisions =
        {
          .collision_id = GKYL_LBO_COLLISIONS,
          .den_ref = 1.0,
          .temp_ref = ctx[s].temperature,
          .eps0 = 1.0,
          .hbar = 0.001,
          .bmag_ref = 1.0,
          .num_cross_collisions = 1,
          .self_nu = prescribed ? frequency : NULL,
          .cross_nu = {prescribed ? frequency : NULL},
        },
    };
    snprintf(inp.species[s].name, sizeof(inp.species[s].name), "s%d", s);
    snprintf(
      inp.species[s].collisions.collide_with[0], sizeof(inp.species[s].collisions.collide_with[0]),
      "s%d", 1 - s
    );
  }
  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&inp);
  gkyl_gyrokinetic_app_apply_ic(app, 0.0);
  for (int s = 0; s < 2; ++s) {
    struct gk_species *sp = &app->species[s];
    gk_species_lbo_moms(app, sp, &sp->lbo, sp->f);
  }
  for (int s = 0; s < 2; ++s) {
    struct gk_species *sp = &app->species[s];
    gk_species_lbo_cross_moms(app, sp, &sp->lbo, sp->f);
    gkyl_array_clear(sp->f1, 0.0);
    gk_species_lbo_rhs(app, sp, &sp->lbo, sp->f, sp->f1);
  }
  return app;
}

static int
compare(
  const char *name, const struct gkyl_array *ref, const struct gkyl_array *other,
  const struct gkyl_range *range, double scale
)
{
  double err = 0.0, norm = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);
    const double *a = gkyl_array_cfetch(ref, loc), *b = gkyl_array_cfetch(other, loc);
    for (int k = 0; k < ref->ncomp; ++k) {
      if (!isfinite(a[k]) || !isfinite(b[k])) {
        return 1;
      }
      err = fmax(err, fabs(a[k] - b[k] / scale));
      norm = fmax(norm, fabs(a[k]));
    }
  }
  if (err > 1.e-8 * fmax(norm, 1.e-12)) {
    fprintf(stderr, "%s: relative error %.12g\n", name, err / norm);
    return 1;
  }
  return 0;
}

// Wide velocity bounds make boundary corrections negligible in this rate check.
static int
check_exchange(gkyl_gyrokinetic_app *app)
{
  int failures = 0;
  for (int s = 0; s < 2; ++s) {
    struct gk_species *sp = &app->species[s];
    double mass = sp->info.mass, other_mass = app->species[1 - s].info.mass;
    double temp = s == 0 ? 1.0 : 2.0, other_temp = s == 0 ? 2.0 : 1.0;
    double expected_flow = s == 0 ? -0.2 : 0.3;
    double expected_vtsq =
      temp / mass + (other_temp - temp + other_mass * 0.5 * 0.5 / 3.0) / (mass + other_mass);
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &app->local);
    while (gkyl_range_iter_next(&iter)) {
      long loc = gkyl_range_idx(&app->local, iter.idx);
      const double *self = gkyl_array_cfetch(sp->lbo.prim_moms, loc);
      const double *total = gkyl_array_cfetch(sp->lbo.nu_prim_moms, loc);
      const double *nu_self = gkyl_array_cfetch(sp->lbo.self_nu, loc);
      const double *nu_cross = gkyl_array_cfetch(sp->lbo.cross_nu[0], loc);
      double cross_flow = (total[0] - nu_self[0] * self[0] / sqrt(2.0)) / nu_cross[0];
      int nb = app->basis.num_basis;
      double cross_vtsq = (total[nb] - nu_self[0] * self[nb] / sqrt(2.0)) / nu_cross[0];
      if (fabs(cross_flow - expected_flow) > 1.e-6 || fabs(cross_vtsq - expected_vtsq) > 1.e-6) {
        fprintf(
          stderr, "species %d: cross (u, vtSq) = (%.12g, %.12g), expected (%.12g, %.12g)\n", s,
          cross_flow, cross_vtsq, expected_flow, expected_vtsq
        );
        ++failures;
      }
    }
  }
  return failures;
}

int
main(void)
{
  int failures = 0;
  for (int prescribed = 0; prescribed < 2; ++prescribed) {
    for (int varying = 0; varying < 2; ++varying) {
      struct test_ctx ctx[2] = {
        {
          .scale = 1.0,
          .curvature = varying ? 0.2 : 0.0,
          .density = 1.0,
          .flow = 0.3,
          .temperature = 1.0,
        },
        {.density = 1.0, .flow = -0.2, .temperature = 2.0},
      };
      gkyl_gyrokinetic_app *ref = make_app(ctx, prescribed);
      ctx[0].scale = 4.0;
      gkyl_gyrokinetic_app *scaled = make_app(ctx, prescribed);
      if (!prescribed && !varying) {
        failures += check_exchange(ref);
        failures += check_exchange(scaled);
      }
      for (int s = 0; s < 2; ++s) {
        struct gk_species *a = &ref->species[s], *b = &scaled->species[s];
        failures +=
          compare("self primitives", a->lbo.prim_moms, b->lbo.prim_moms, &ref->local, 1.0);
        failures +=
          compare("cross frequency", a->lbo.cross_nu[0], b->lbo.cross_nu[0], &ref->local, 1.0);
        failures +=
          compare("exchange coefficient", a->lbo.alpha_E, b->lbo.alpha_E, &ref->local, 4.0);
        failures += compare(
          "assembled primitives", a->lbo.nu_prim_moms, b->lbo.nu_prim_moms, &ref->local, 1.0
        );
        failures += compare("collision RHS", a->f1, b->f1, &a->local, 4.0);
      }
      gkyl_gyrokinetic_app_release(ref);
      gkyl_gyrokinetic_app_release(scaled);
    }
  }
  printf("LBO Jacobian regression: %s (%d failures)\n", failures ? "FAIL" : "PASS", failures);
  return failures ? 1 : 0;
}
