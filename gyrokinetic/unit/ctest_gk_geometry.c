#include <acutest.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_array_ops.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>

struct query_ctx {
  int ndim, order, field;
};

/**
 * Supply independently known polynomial fields for modal-query tests.
 *
 * @param t Unused time.
 * @param xc Computational coordinates.
 * @param out Field values.
 * @param vctx Dimension, order and field selector.
 */
static void
query_polynomial(double t, const double *xc, double *out, void *vctx)
{
  struct query_ctx *ctx = vctx;
  double value = 0;
  for (int d = 0; d < ctx->ndim; ++d) {
    value += (d + 1) * xc[d] + (ctx->order == 2 ? 0.1 * xc[d] * xc[d] : 0);
  }
  if (ctx->field == 0) {
    out[0] = 2 + 0.1 * value;
  } else if (ctx->field == 1) {
    out[0] = 0.2 + 0.01 * value;
    out[1] = -0.7;
    out[2] = 0.4;
  } else {
    out[0] = 2 + 0.05 * value;
    out[1] = 0.3;
    out[2] = 7.2;
  }
}

/**
 * Allocate just the arrays used by queries, with a ghosted index layout.
 *
 * @param ndim Configuration dimension.
 * @param order Polynomial order.
 * @return Manufactured host geometry.
 */
static struct gk_geometry
query_geometry(int ndim, int order)
{
  struct gk_geometry geo = {.geometry_id = GKYL_GEOMETRY_TOKAMAK};
  gkyl_rect_grid_init(&geo.grid, ndim, (double[]){0, 0, 0}, (double[]){1, 1, 1}, (int[]){4, 4, 4});
  gkyl_cart_modal_serendip(&geo.basis, ndim, order);
  gkyl_create_grid_ranges(&geo.grid, (int[]){1, 1, 1}, &geo.local_ext, &geo.local);
  geo.geo_int.bmag = gkyl_array_new(GKYL_DOUBLE, geo.basis.num_basis, geo.local_ext.volume);
  geo.geo_int.bcart = gkyl_array_new(GKYL_DOUBLE, 3 * geo.basis.num_basis, geo.local_ext.volume);
  geo.geo_int.mc2p = gkyl_array_new(GKYL_DOUBLE, 3 * geo.basis.num_basis, geo.local_ext.volume);
  struct gkyl_array *arrays[] = {geo.geo_int.bmag, geo.geo_int.bcart, geo.geo_int.mc2p};
  for (int field = 0; field < 3; ++field) {
    struct query_ctx ctx = {ndim, order, field};
    gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(
      &geo.grid, &geo.basis, order + 1, field == 0 ? 1 : 3, query_polynomial, &ctx
    );
    gkyl_proj_on_basis_advance(proj, 0, &geo.local, arrays[field]);
    gkyl_proj_on_basis_release(proj);
  }
  return geo;
}

/**
 * Free the manufactured geometry's three arrays.
 *
 * @param geo Manufactured geometry.
 */
static void
query_release(struct gk_geometry *geo)
{
  gkyl_array_release(geo->geo_int.bmag);
  gkyl_array_release(geo->geo_int.bcart);
  gkyl_array_release(geo->geo_int.mc2p);
}

/**
 * Compare polynomial field queries in 1x, 2x and 3x at orders one and two.
 *
 */
static void
test_modal_queries(void)
{
  for (int ndim = 1; ndim <= 3; ++ndim) {
    for (int order = 1; order <= 2; ++order) {
      struct gk_geometry geo = query_geometry(ndim, order);
      for (int point = 0; point <= 10; ++point) {
        double xc[] = {0.1 * point, 0.1 * point, 0.1 * point};
        double bmag, expected_bmag, bhat[3], map[3], actual[3], expected[3];
        struct query_ctx ctx = {ndim, order, 0};
        query_polynomial(0, xc, &expected_bmag, &ctx);
        ctx.field = 1;
        query_polynomial(0, xc, bhat, &ctx);
        ctx.field = 2;
        query_polynomial(0, xc, map, &ctx);
        gkyl_gk_geometry_geo_int_bmag_at_comp_coord(&geo, xc, &bmag);
        TEST_CHECK(gkyl_compare(bmag, expected_bmag, 2e-13));
        gkyl_gk_geometry_geo_int_bfieldcart_at_comp_coord(&geo, xc, actual);
        for (int d = 0; d < 3; ++d) {
          TEST_CHECK(gkyl_compare(actual[d], expected_bmag * bhat[d], 2e-13));
        }
        gkyl_gk_geometry_geo_int_comp2cyl_coords(&geo, xc, actual);
        TEST_CHECK(gkyl_compare(actual[0], map[0], 2e-13));
        TEST_CHECK(gkyl_compare(actual[1], map[2], 2e-13)); // Unwrapped phi > 2pi.
        TEST_CHECK(gkyl_compare(actual[2], map[1], 2e-13));
        expected[0] = bmag * (cos(map[2]) * bhat[0] + sin(map[2]) * bhat[1]);
        expected[1] = bmag * (-sin(map[2]) * bhat[0] + cos(map[2]) * bhat[1]);
        expected[2] = bmag * bhat[2];
        gkyl_gk_geometry_geo_int_bfieldcyl_at_comp_coord(&geo, xc, actual);
        for (int d = 0; d < 3; ++d) {
          TEST_CHECK(gkyl_compare(actual[d], expected[d], 2e-13));
        }
        gkyl_gk_geometry_geo_int_bfieldtor_at_comp_coord(&geo, xc, actual);
        TEST_CHECK(gkyl_compare(actual[0], expected[1], 2e-13));
        TEST_CHECK(gkyl_compare(actual[1], hypot(expected[0], expected[2]), 2e-13));
      }
      query_release(&geo);
    }
  }
}

/**
 * Set a vector or scalar constant in every modal cell.
 *
 * @param geo Manufactured geometry.
 * @param array Array to set.
 * @param values Constant components.
 */
static void
query_constant(struct gk_geometry *geo, struct gkyl_array *array, const double *values)
{
  gkyl_array_clear(array, 0);
  int nb = geo->basis.num_basis;
  for (long cell = 0; cell < array->size; ++cell) {
    double *coeff = gkyl_array_fetch(array, cell);
    for (int comp = 0; comp < array->ncomp / nb; ++comp) {
      coeff[comp * nb] = values[comp] * pow(2, 0.5 * geo->basis.ndim);
    }
  }
}

/**
 * Check spatial conversions, shaping, field reversal, and zero poloidal field.
 *
 */
static void
test_coordinate_queries(void)
{
  struct gk_geometry geo = query_geometry(2, 1);
  double xc[] = {0.3, 0.8}, out[3];
  struct gkyl_gen_miller_eqdsk_inp miller = {.Rmajor_axis = 2, .Zvert_axis = 0.2, .kappa = 1};
  for (int shape = 0; shape < 3; ++shape) {
    miller.kappa = shape ? 1.5 : 0; // Zero selects the generator default.
    miller.delta = shape ? 0.35 : 0;
    miller.a_shift = shape == 2 ? 0.4 : 0;
    for (int angle = -15; angle <= 15; ++angle) {
      double radius = 0.4, theta = angle * 0.2, phi = -0.6;
      double Rmajor = 2 - miller.a_shift * radius * radius / 4 +
                      radius * cos(theta + asin(miller.delta) * sin(theta));
      double Zvert = 0.2 + (shape ? 1.5 : 1) * radius * sin(theta);
      for (int cartesian = 0; cartesian < 2; ++cartesian) {
        geo.geometry_id = cartesian ? GKYL_GEOMETRY_MAPC2P : GKYL_GEOMETRY_TOKAMAK;
        double mapped[3] = {Rmajor, Zvert, phi};
        if (cartesian) {
          mapped[0] = Rmajor * cos(phi);
          mapped[1] = Rmajor * sin(phi);
          mapped[2] = Zvert;
        }
        query_constant(&geo, geo.geo_int.mc2p, mapped);
        gkyl_gk_geometry_geo_int_comp2miller_coords(&geo, xc, &miller, out);
        TEST_CHECK(gkyl_compare(out[0], radius, 2e-12));
        TEST_CHECK(gkyl_compare(out[1], theta, 2e-12));
        TEST_CHECK(gkyl_compare(out[2], phi, 2e-12));
        gkyl_gk_geometry_geo_int_comp2tor_coords(&geo, xc, 2, 0.2, out);
        TEST_CHECK(gkyl_compare(out[0], hypot(Rmajor - 2, Zvert - 0.2), 2e-12));
        TEST_CHECK(gkyl_compare(out[1], atan2(Zvert - 0.2, Rmajor - 2), 2e-12));
        for (int sign = -1; sign <= 1; sign += 2) {
          query_constant(&geo, geo.geo_int.bmag, (double[]){2});
          query_constant(&geo, geo.geo_int.bcart, (double[]){-sign * sin(phi), sign * cos(phi), 0});
          gkyl_gk_geometry_geo_int_bfieldtor_at_comp_coord(&geo, xc, out);
          TEST_CHECK(gkyl_compare(out[0], 2 * sign, 2e-13));
          TEST_CHECK(fabs(out[1]) < 1e-13);
        }
      }
    }
  }
  // A valid root near the shift-induced fold must not be rejected because
  // the initial bracketing step overshoots the fold.
  geo.geometry_id = GKYL_GEOMETRY_TOKAMAK;
  miller = (struct gkyl_gen_miller_eqdsk_inp){.Rmajor_axis = 2, .kappa = 1, .a_shift = 2};
  double radius = 0.9999;
  query_constant(&geo, geo.geo_int.mc2p, (double[]){2 + radius - 0.5 * radius * radius, 0, 0});
  gkyl_gk_geometry_geo_int_comp2miller_coords(&geo, xc, &miller, out);
  TEST_CHECK(gkyl_compare(out[0], radius, 1e-10));
  query_release(&geo);
}

/**
 * Check deterministic traces, including a local MPI-like subdomain boundary.
 *
 */
static void
test_query_traces(void)
{
  struct gk_geometry geo = query_geometry(1, 1);
  query_constant(&geo, geo.geo_int.bmag, (double[]){1});
  double *cell = gkyl_array_fetch(geo.geo_int.bmag, gkyl_range_idx(&geo.local, (int[]){3}));
  cell[0] = 3 * sqrt(2.0);
  double out;
  gkyl_gk_geometry_geo_int_bmag_at_comp_coord(&geo, (double[]){0.5}, &out);
  TEST_CHECK(gkyl_compare(out, 3, 1e-14));
  gkyl_sub_range_init(&geo.local, &geo.local_ext, (int[]){1}, (int[]){2});
  gkyl_gk_geometry_geo_int_bmag_at_comp_coord(&geo, (double[]){0.5}, &out);
  TEST_CHECK(gkyl_compare(out, 1, 1e-14));
  query_release(&geo);
}

/**
 * Verify fail-fast behavior in child processes without terminating the suite.
 *
 */
static void
test_query_failures(void)
{
  struct gk_geometry geo = query_geometry(1, 1);
  for (int scenario = 0; scenario < 9; ++scenario) {
    int pipefd[2];
    TEST_ASSERT(pipe(pipefd) == 0);
    fflush(NULL);
    pid_t pid = fork();
    TEST_ASSERT(pid >= 0);
    if (pid == 0) {
      close(pipefd[0]);
      dup2(pipefd[1], STDERR_FILENO);
      close(pipefd[1]);
      double xc[] = {0.3}, out[3];
      struct gkyl_gen_miller_eqdsk_inp miller = {.Rmajor_axis = 2, .kappa = 1};
      switch (scenario) {
        case 0:
          gkyl_gk_geometry_geo_int_bmag_at_comp_coord(NULL, xc, out);
          break;
        case 1:
          xc[0] = 1.1;
          gkyl_gk_geometry_geo_int_bmag_at_comp_coord(&geo, xc, out);
          break;
        case 2:
          xc[0] = strtod("nan", NULL);
          gkyl_gk_geometry_geo_int_bmag_at_comp_coord(&geo, xc, out);
          break;
        case 3:
          geo.geometry_id = GKYL_GEOMETRY_FROMFILE;
          gkyl_gk_geometry_geo_int_comp2cyl_coords(&geo, xc, out);
          break;
        case 4:
          query_constant(&geo, geo.geo_int.mc2p, (double[]){2, 0, 0});
          gkyl_gk_geometry_geo_int_comp2tor_coords(&geo, xc, 2, 0, out);
          break;
        case 5:
          miller.kappa = -1;
          gkyl_gk_geometry_geo_int_comp2miller_coords(&geo, xc, &miller, out);
          break;
        case 6:
          query_constant(&geo, geo.geo_int.mc2p, (double[]){10, 0, 0});
          miller.a_shift = 2;
          gkyl_gk_geometry_geo_int_comp2miller_coords(&geo, xc, &miller, out);
          break;
        case 7:
          query_constant(&geo, geo.geo_int.mc2p, (double[]){0, 0, 0});
          gkyl_gk_geometry_geo_int_comp2cyl_coords(&geo, xc, out);
          break;
        case 8:
          gkyl_sub_range_init(&geo.local, &geo.local_ext, (int[]){3}, (int[]){4});
          gkyl_gk_geometry_geo_int_bmag_at_comp_coord(&geo, xc, out);
          break;
      }
      _exit(0);
    }
    close(pipefd[1]);
    char message[512] = {0};
    ssize_t count = read(pipefd[0], message, sizeof message - 1);
    close(pipefd[0]);
    int status;
    TEST_ASSERT(waitpid(pid, &status, 0) == pid);
    TEST_CHECK(WIFEXITED(status) && WEXITSTATUS(status) == EXIT_FAILURE);
    TEST_MSG("scenario=%d status=%d message=%s", scenario, status, message);
    TEST_CHECK(count > 0 && strstr(message, "geo_int query:") != NULL);
  }
  query_release(&geo);
}

TEST_LIST = {
  {"modal_queries", test_modal_queries},
  {"coordinate_queries", test_coordinate_queries},
  {"query_traces", test_query_traces},
  {"query_failures", test_query_failures},
  {NULL, NULL}
};
