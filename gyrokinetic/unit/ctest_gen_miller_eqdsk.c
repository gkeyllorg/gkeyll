#include <acutest.h>
#include <gkyl_gen_miller_eqdsk_priv.h>
#include <gkyl_dg_basis_ops.h>
#include <gkyl_efit.h>
#include <gkyl_gk_geometry_tok.h>
#include <gkyl_position_map.h>
#include <gkyl_tok_geo.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h> // mkstemp(): create a uniquely named temporary file.
#include <unistd.h> // unlink(): remove a file's directory entry.

/**
 * Return a constant positive safety-factor magnitude.
 *
 * @param R_omp Unused outboard major radius [m].
 * @param ctx Unused context.
 * @return Safety-factor magnitude.
 */
static double
constant_q(double R_omp, void *ctx)
{
  return 2;
}

/**
 * Evaluate the LTX safety-factor magnitude.
 *
 * @param R_omp Outboard major radius [m].
 * @param ctx Unused context.
 * @return Safety-factor magnitude.
 */
static double
ltx_q(double R_omp, void *ctx)
{
  double radius = R_omp - 0.406051632, a = 0.61183 - 0.406051632;
  return 4.3153848 * exp(0.5 * 2.6899871 * (radius * radius / (a * a) - 1));
}

/**
 * Evaluate the DIII-D safety-factor magnitude.
 *
 * @param R_omp Outboard major radius [m].
 * @param ctx Unused context.
 * @return Safety-factor magnitude.
 */
static double
d3d_q(double R_omp, void *ctx)
{
  return ((49.46395467479657 * R_omp - 260.79513158768754) * R_omp + 458.42618139184754) * R_omp -
         267.63441353752336;
}

/**
 * Evaluate the TCV safety-factor magnitude.
 *
 * @param R_omp Outboard major radius [m].
 * @param ctx Unused context.
 * @return Safety-factor magnitude.
 */
static double
tcv_q(double R_omp, void *ctx)
{
  return ((497.3420166252413 * R_omp - 1408.736172826569) * R_omp + 1331.4134861681464) * R_omp -
         419.00692601227627;
}

/**
 * Return a nonfinite safety factor to exercise validation.
 *
 * @param R_omp Unused outboard major radius [m].
 * @param ctx Unused context.
 * @return NaN.
 */
static double
invalid_q(double R_omp, void *ctx)
{
  // Parse at runtime so -ffast-math cannot fold a constant NaN return.
  return strtod("nan", NULL);
}

/**
 * Check analytic circular flux, both field signs, and invalid inputs.
 *
 */
static void
test_circle(void)
{
  struct gkyl_gen_miller_eqdsk_inp inp = {
    .Rmajor_axis = 2,
    .bmag_axis = 1.5,
    .Rmajor_LCFS = 2.5,
    .qsafety = constant_q,
  };

  gkyl_gen_miller_eqdsk *up = gkyl_gen_miller_eqdsk_inew(&inp);
  TEST_ASSERT(up != NULL);

  const struct gkyl_gen_miller_eqdsk_info *info = gkyl_gen_miller_eqdsk_get_info(up);

  for (int j = 0; j <= 100; ++j) {
    double radius = info->rminor_guard * j / 100;
    double exact = 1.5 * radius * radius / (2 + sqrt(4 - radius * radius));
    double flux = gkyl_gen_miller_eqdsk_psi_at_R_omp(up, 2 + radius);
    TEST_CHECK(fabs(flux - exact) < 1e-8 * info->psi_boundary);
    TEST_CHECK(fabs(gkyl_gen_miller_eqdsk_R_omp_at_psi(up, flux) - (2 + radius)) < 1e-8);
  }

  TEST_CHECK(gkyl_gen_miller_eqdsk_advance(up));

  gkyl_gen_miller_eqdsk_release(up);

  inp.bmag_axis = -1.5;
  up = gkyl_gen_miller_eqdsk_inew(&inp);
  TEST_ASSERT(up != NULL);
  double rminor = 0.3;
  double exact = 1.5 * rminor * rminor / (2 + sqrt(4 - rminor * rminor));
  TEST_CHECK(fabs(gkyl_gen_miller_eqdsk_psi_at_r(up, rminor) - exact) < 1e-10);
  gkyl_gen_miller_eqdsk_release(up);
  double invalid_fields[] = {0, strtod("nan", NULL), strtod("inf", NULL)};
  for (int j = 0; j < 3; ++j) {
    inp.bmag_axis = invalid_fields[j];
    TEST_CHECK(gkyl_gen_miller_eqdsk_inew(&inp) == NULL);
    TEST_MSG("Invalid bmag_axis case %d: %.17g", j, inp.bmag_axis);
  }
  inp.bmag_axis = 1.5;
  inp.qsafety = invalid_q;
  TEST_CHECK(gkyl_gen_miller_eqdsk_inew(&inp) == NULL);

  inp.qsafety = constant_q;
  inp.a_shift = 10;
  TEST_CHECK(gkyl_gen_miller_eqdsk_inew(&inp) == NULL);

  inp.a_shift = 0;
  inp.num_nodes_Rmajor = 64;
  TEST_CHECK(gkyl_gen_miller_eqdsk_inew(&inp) == NULL);
}

/**
 * Verify GEQDSK profiles and reconstructed numerical tokamak geometry.
 *
 * @param filename Generated GEQDSK path.
 * @param inp Generator parameters.
 * @param generator Generator supplying reference flux and coverage.
 */
static void
check_reader(
  const char *filename, const struct gkyl_gen_miller_eqdsk_inp *inp,
  const gkyl_gen_miller_eqdsk *generator
)
{
  struct gkyl_efit_inp efit_inp = {.rz_poly_order = 2, .flux_poly_order = 1, .reflect = false};
  snprintf(efit_inp.filepath, sizeof efit_inp.filepath, "%s", filename);

  struct gkyl_efit *efit = gkyl_efit_new(&efit_inp);

  const struct gkyl_gen_miller_eqdsk_info *info = gkyl_gen_miller_eqdsk_get_info(generator);

  TEST_CHECK(fabs(efit->sibry / info->psi_boundary - 1) < 1e-8);
  TEST_CHECK(fabs(efit->zmaxis - inp->Zvert_axis) < 1e-9);

  // The axis is a degenerate critical point. Check any reported X-points lie
  // outside the physical simulation annulus, rather than rejecting the axis.
  for (int j = 0; j < efit->num_xpts_cubic; ++j) {
    double point[2] = {efit->Rxpt_cubic[j], efit->Zxpt_cubic[j]}, value[3];
    efit->evf->eval_cubic_wgrad(0, point, value, efit->evf->ctx);
    double lower = gkyl_gen_miller_eqdsk_psi_at_R_omp(generator, inp->Rmajor_omp_bounds[0]);
    double upper = gkyl_gen_miller_eqdsk_psi_at_R_omp(generator, inp->Rmajor_omp_bounds[1]);
    TEST_CHECK(value[0] < lower || value[0] > upper);
  }

  double max_bp_error = 0, max_flux_error = 0, max_q_error = 0;
  for (int surface = 0; surface < 5; ++surface) {
    double R_omp = inp->Rmajor_omp_bounds[0] +
                   (inp->Rmajor_omp_bounds[1] - inp->Rmajor_omp_bounds[0]) * surface / 4;
    double flux = gkyl_gen_miller_eqdsk_psi_at_R_omp(generator, R_omp);
    double offset = R_omp - inp->Rmajor_axis;
    double radius = 2 * offset / (1 + sqrt(1 - 2 * inp->a_shift * offset / inp->Rmajor_axis));
    double h = radius * 1e-5;
    double omp_lo = inp->Rmajor_axis -
                    inp->a_shift * (radius - h) * (radius - h) / (2 * inp->Rmajor_axis) + radius -
                    h;
    double omp_hi = inp->Rmajor_axis -
                    inp->a_shift * (radius + h) * (radius + h) / (2 * inp->Rmajor_axis) + radius +
                    h;
    double dpsi = (gkyl_gen_miller_eqdsk_psi_at_R_omp(generator, omp_hi) -
                   gkyl_gen_miller_eqdsk_psi_at_R_omp(generator, omp_lo)) /
                  (2 * h);
    double q_integral = 0;
    for (int j = 0; j < 1024; ++j) {
      double theta = 2 * M_PI * j / 1024, rz[2], value[3];
      gkyl_gen_miller_eqdsk_RZ_at_rtheta(generator, radius, theta, rz);
      efit->evf->eval_cubic_wgrad(0, rz, value, efit->evf->ctx);
      double eta = asin(inp->delta), angle = theta + eta * sin(theta);
      double Rr = cos(angle) - inp->a_shift * radius / inp->Rmajor_axis;
      double Rt = -radius * (1 + eta * cos(theta)) * sin(angle);
      double Zr = inp->kappa * sin(theta), Zt = inp->kappa * radius * cos(theta);
      double dl = hypot(Rt, Zt), bp = hypot(value[1], value[2]) / rz[0];
      double expected_bp = dpsi * dl / (rz[0] * (Rr * Zt - Rt * Zr));
      max_bp_error = fmax(max_bp_error, fabs(bp / expected_bp - 1));
      max_flux_error = fmax(max_flux_error, fabs(value[0] - flux) / info->psi_boundary);
      q_integral += inp->Rmajor_axis * inp->bmag_axis * dl / (rz[0] * rz[0] * bp) / 1024;
    }
    max_q_error = fmax(
      max_q_error,
      fabs(q_integral / copysign(inp->qsafety(R_omp, inp->qsafety_ctx), inp->bmag_axis) - 1)
    );
  }

  TEST_CHECK(max_bp_error < 1e-3);
  TEST_MSG("Bp error = %.9g", max_bp_error);
  TEST_CHECK(max_flux_error < 1e-5);
  TEST_MSG("flux error = %.9g", max_flux_error);
  TEST_CHECK(max_q_error < 1e-3);
  TEST_MSG("q error = %.9g", max_q_error);

  TEST_CHECK(fabs(efit->bcentr - inp->bmag_axis) < 1e-8);
  struct gkyl_range_iter fit;
  gkyl_range_iter_init(&fit, &efit->fluxlocal);
  while (gkyl_range_iter_next(&fit)) {
    long loc = gkyl_range_idx(&efit->fluxlocal, fit.idx);
    const double *F = gkyl_array_cfetch(efit->fpolflux, loc);
    const double *q = gkyl_array_cfetch(efit->qflux, loc);
    TEST_CHECK(fabs(F[0] / sqrt(2) - inp->Rmajor_axis * inp->bmag_axis) < 1e-8);
    TEST_CHECK(q[0] * inp->bmag_axis > 0);
  }
  gkyl_efit_release(efit);

  // Check roots of the actual quadratic and cubic EFIT reconstructions.
  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_GEOMETRY_TOKAMAK_IWL,
    .rleft = info->RmajorZvert_bounds[0],
    .rright = info->RmajorZvert_bounds[1],
    .rmin = info->RmajorZvert_bounds[0],
    .rmax = info->RmajorZvert_bounds[1],
    .zmin = info->RmajorZvert_bounds[2],
    .zmax = info->RmajorZvert_bounds[3],
  };

  for (int cubic = 0; cubic < 2; ++cubic) {
    ginp.use_cubics = cubic;
    ginp.use_hyperbolic_numbers = true;
    struct gkyl_tok_geo *geo = gkyl_tok_geo_new(&efit_inp, &ginp);
    geo->rmin = ginp.rmin;
    geo->rmax = ginp.rmax;
    double R[4], slope[4], dR[4], dZ[4];
    int roots = gkyl_tok_geo_R_psiZ(geo, info->psi_boundary, inp->Zvert_axis, 4, R, slope, dR, dZ);
    TEST_CHECK(roots == 2);
    double error = 1e9;
    for (int j = 0; j < roots; ++j) {
      error = fmin(error, fabs(R[j] - inp->Rmajor_LCFS));
    }
    TEST_CHECK(error < 1e-4 * info->rminor_LCFS);
    gkyl_tok_geo_release(geo);
  }

  ginp.use_cubics = false;
  struct gkyl_gk_geometry_inp geometry = {
    .geometry_id = GKYL_GEOMETRY_TOKAMAK,
    .efit_info = efit_inp,
    .tok_grid_info = ginp,
    .has_LCFS = true,
    .x_LCFS = info->psi_boundary,
  };
  double lower[] = {0.95 * info->psi_boundary, -0.01, -M_PI + 1e-10};
  double upper[] = {1.05 * info->psi_boundary, 0.01, M_PI - 1e-10};
  int cells[] = {2, 1, 2}, ghosts[] = {1, 1, 1};

  gkyl_rect_grid_init(&geometry.grid, 3, lower, upper, cells);
  gkyl_cart_modal_serendip(&geometry.basis, 3, 1);
  gkyl_create_grid_ranges(&geometry.grid, ghosts, &geometry.local_ext, &geometry.local);
  geometry.global = geometry.local;
  geometry.global_ext = geometry.local_ext;
  geometry.geo_grid = geometry.grid;
  geometry.geo_basis = geometry.basis;
  geometry.geo_local = geometry.local;
  geometry.geo_local_ext = geometry.local_ext;
  geometry.geo_global = geometry.global;
  geometry.geo_global_ext = geometry.global_ext;
  geometry.position_map = gkyl_position_map_null_new();

  struct gk_geometry *gk_geo = gkyl_gk_geometry_tok_new(&geometry);
  TEST_CHECK(gkyl_gk_geometry_reduce_bmag(gk_geo, GKYL_MIN) > 0);

  struct gkyl_range_iter nodes;
  gkyl_range_iter_init(&nodes, &gk_geo->nrange_int);
  while (gkyl_range_iter_next(&nodes)) {
    long loc = gkyl_range_idx(&gk_geo->nrange_int, nodes.idx);
    const double *position = gkyl_array_cfetch(gk_geo->geo_int.mc2p_nodal, loc);
    const double *bhat = gkyl_array_cfetch(gk_geo->geo_int.bcart_nodal, loc);
    const double *bmag = gkyl_array_cfetch(gk_geo->geo_int.bmag_nodal, loc);
    const double *tangent = gkyl_array_cfetch(gk_geo->geo_int.dxdz_nodal, loc);
    double phi = position[2], Rmajor = position[0];
    double Bphi = bmag[0] * (-bhat[0] * sin(phi) + bhat[1] * cos(phi));
    TEST_CHECK(fabs(Bphi * Rmajor / (inp->Rmajor_axis * inp->bmag_axis) - 1) < 1e-7);
    double normal[3] = {
      tangent[4] * tangent[8] - tangent[5] * tangent[7],
      tangent[5] * tangent[6] - tangent[3] * tangent[8],
      tangent[3] * tangent[7] - tangent[4] * tangent[6]
    };
    TEST_CHECK(tangent[0] * normal[0] + tangent[1] * normal[1] + tangent[2] * normal[2] > 0);
    // The field follows increasing poloidal coordinate for either toroidal sign.
    TEST_CHECK(bhat[2] * tangent[8] >= 0);
  }
  gkyl_gk_geometry_release(gk_geo);
  gkyl_position_map_release(geometry.position_map);
}

/**
 * Compute unit-radius Miller arc length independently for reference quadrature.
 *
 * @param theta Miller angle [rad].
 * @param ctx Generator input containing elongation and triangularity.
 * @return Arc length per unit radius and angle.
 */
static double
reference_arc_speed(double theta, void *ctx)
{
  const struct gkyl_gen_miller_eqdsk_inp *inp = ctx;
  double eta = asin(inp->delta);
  double Rt = -(1 + eta * cos(theta)) * sin(theta + eta * sin(theta));
  return sqrt(Rt * Rt + pow(inp->kappa * cos(theta), 2));
}

/**
 * Integrate reference arc length with a separately resolved composite midpoint rule.
 *
 * @param inp Miller shape parameters.
 * @param theta Upper Miller angle [rad]; the lower limit is zero.
 * @return Unit-radius poloidal arc length.
 */
static double
reference_arc(struct gkyl_gen_miller_eqdsk_inp *inp, double theta)
{
  int count = 131072;
  double step = theta / count, integral = 0;
  for (int j = 0; j < count; ++j) {
    integral += reference_arc_speed((j + 0.5) * step, inp);
  }
  return step * integral;
}

/**
 * Check flux inversion and arc-coordinate conversion for circular and shifted TCV shapes.
 *
 */
static void
test_coordinates(void)
{
  struct gkyl_gen_miller_eqdsk_inp inp = {
    .Rmajor_axis = 2,
    .bmag_axis = 1.5,
    .Rmajor_LCFS = 2.5,
    .kappa = 1,
    .qsafety = constant_q,
  };
  for (int shape = 0; shape < 2; ++shape) {
    if (shape) {
      inp.kappa = 1.4;
      inp.delta = -0.38;
      inp.a_shift = 0.5;
    }
    gkyl_gen_miller_eqdsk *up = gkyl_gen_miller_eqdsk_inew(&inp);
    TEST_ASSERT(up != NULL);
    const struct gkyl_gen_miller_eqdsk_info *info = gkyl_gen_miller_eqdsk_get_info(up);
    for (int j = 0; j <= 10; ++j) {
      double radius = info->rminor_guard * j / 10;
      double psi = gkyl_gen_miller_eqdsk_psi_at_r(up, radius);
      TEST_CHECK(fabs(gkyl_gen_miller_eqdsk_r_at_psi(up, psi) - radius) < 1e-9);
      double rz[2];
      gkyl_gen_miller_eqdsk_RZ_at_rtheta(up, radius, 0, rz);
      TEST_CHECK(fabs(gkyl_gen_miller_eqdsk_psi_at_R_omp(up, rz[0]) - psi) < 1e-12);
    }
    double half_arc = reference_arc(&inp, M_PI);
    for (int j = -10; j <= 10; ++j) {
      double theta_tok = M_PI * j / 10;
      double theta = gkyl_gen_miller_eqdsk_theta_at_theta_tok(up, theta_tok);
      double arc = theta == 0 ? 0 : reference_arc(&inp, fabs(theta));
      TEST_CHECK(fabs(M_PI * arc / half_arc - fabs(theta_tok)) < 1e-10);
      TEST_MSG(
        "shape=%d, theta_tok=%.17g, theta=%.17g, arc residual=%.17g", shape, theta_tok, theta,
        M_PI * arc / half_arc - fabs(theta_tok)
      );
      TEST_CHECK(theta * theta_tok >= 0);
      if (!shape) {
        TEST_CHECK(fabs(theta - theta_tok) < 1e-12);
      }
    }

    gkyl_gen_miller_eqdsk_release(up);
  }
}

/**
 * Check both orientations for LTX, DIII-D, and TCV geometries.
 *
 */
static void
test_devices(void)
{
  const char *names[] = {"ltx", "d3d", "tcv"};
  struct gkyl_gen_miller_eqdsk_inp inputs[] = {
    {
      .Rmajor_axis = 0.406051632,
      .bmag_axis = 0.240606108,
      .Rmajor_LCFS = 0.61183,
      .kappa = 1.3,
      .delta = 0.4,
      .qsafety = ltx_q,
      .has_Rmajor_omp_bounds = true,
      .Rmajor_omp_bounds = {0.58, 0.7},
    },
    {
      .Rmajor_axis = 1.6,
      .Zvert_axis = 0.013055028,
      .bmag_axis = 2 * 1.6486461 / 1.6,
      .Rmajor_LCFS = 2.17,
      .kappa = 1.35,
      .delta = 0.4,
      .qsafety = d3d_q,
      .has_Rmajor_omp_bounds = true,
      .Rmajor_omp_bounds = {2.17 - 5 * 0.15 / 8, 2.17 + 3 * 0.15 / 8},
    },
    {
      .Rmajor_axis = 0.8867856264,
      .Zvert_axis = 0.1414361745,
      .bmag_axis = 1.4,
      .Rmajor_LCFS = 1.0870056099999,
      .kappa = 1.4,
      .delta = -0.38,
      .a_shift = 0.5,
      .qsafety = tcv_q,
      .has_Rmajor_omp_bounds = true,
      .Rmajor_omp_bounds = {1.0513586729299085, 1.1555930803275622},
    }
  };

  // Independent Python values at the same tolerance and physical parameters.
  double expected_boundary[] = {0.003518612110703164, 0.28091394499783906, 0.03483576483233173};
  for (int orientation = 0; orientation < 2; ++orientation) {
    for (int device = 0; device < 3; ++device) {
      inputs[device].bmag_axis = (orientation ? -1 : 1) * fabs(inputs[device].bmag_axis);
      char filename[] = "/tmp/gkyl_miller_XXXXXX";
      int descriptor = mkstemp(filename);
      TEST_ASSERT(descriptor >= 0);
      close(descriptor);
      gkyl_gen_miller_eqdsk *up = gkyl_gen_miller_eqdsk_inew(&inputs[device]);

      TEST_ASSERT(up != NULL);
      TEST_CHECK(
        fabs(gkyl_gen_miller_eqdsk_get_info(up)->psi_boundary / expected_boundary[device] - 1) <
        1e-8
      );
      TEST_ASSERT(gkyl_gen_miller_eqdsk_advance(up));

      TEST_ASSERT(gkyl_gen_miller_eqdsk_write(up, filename));
      check_reader(filename, &inputs[device], up);
      unlink(filename);
      // Optional Python-generated files exercise the same native reader checks.
      const char *references = getenv("GKYL_MILLER_REFERENCE_DIR");
      if (references && orientation) {
        char path[1024];
        snprintf(path, sizeof path, "%s/miller_%s.geqdsk", references, names[device]);
        check_reader(path, &inputs[device], up);
      }
      gkyl_gen_miller_eqdsk_release(up);
    }
  }
}

/**
 * Check NaN bits even when the unit test is compiled with fast-math.
 *
 * @param value Return value to classify.
 * @return Whether the value is an IEEE NaN.
 */
static bool
miller_is_nan(double value)
{
  uint64_t bits;
  memcpy(&bits, &value, sizeof bits);
  volatile uint64_t exponent = bits & UINT64_C(0x7ff0000000000000);
  volatile uint64_t mantissa = bits & UINT64_C(0x000fffffffffffff);
  return exponent == UINT64_C(0x7ff0000000000000) && mantissa != 0;
}

/**
 * Check the parameter-only forward angle map against independent arc quadrature
 * and the existing generator inverse, including the small endpoint exclusion.
 *
 */
static void
test_forward_angle(void)
{
  struct gkyl_gen_miller_eqdsk_inp inp = {0};
  for (int shift = -1; shift <= 1; ++shift) {
    inp.a_shift = shift * 0.4;
    for (int kappa_default = 0; kappa_default <= 1; ++kappa_default) {
      inp.kappa = kappa_default;
      for (int j = -16; j <= 16; ++j) {
        double angle = M_PI * j / 16;
        TEST_CHECK(gkyl_gen_miller_eqdsk_theta_tok_at_theta(&inp, angle) == angle);
      }
    }
  }
  for (int shape = 0; shape < 3; ++shape) {
    inp = (struct gkyl_gen_miller_eqdsk_inp){
      .Rmajor_axis = 2,
      .Rmajor_LCFS = 2.6,
      .bmag_axis = -1.5,
      .qsafety = constant_q,
      .kappa = 1.6,
      .delta = (shape - 1) * 0.35,
      .a_shift = 0.4,
    };
    gkyl_gen_miller_eqdsk *up = gkyl_gen_miller_eqdsk_inew(&inp);
    TEST_ASSERT(up != NULL);
    double half_arc = reference_arc(&inp, M_PI), previous = -4;
    for (int j = -16; j <= 16; ++j) {
      double angle = M_PI * j / 16;
      double theta = gkyl_gen_miller_eqdsk_theta_tok_at_theta(&inp, angle);
      double expected = copysign(M_PI * reference_arc(&inp, fabs(angle)) / half_arc, angle);
      TEST_CHECK(!miller_is_nan(theta) && theta > previous);
      TEST_CHECK(fabs(theta - expected) < 1e-10);
      TEST_MSG("shape=%d angle=%.17g theta=%.17g reference=%.17g", shape, angle, theta, expected);
      TEST_CHECK(fabs(gkyl_gen_miller_eqdsk_theta_at_theta_tok(up, theta) - angle) < 1e-10);
      TEST_CHECK(theta == -gkyl_gen_miller_eqdsk_theta_tok_at_theta(&inp, -angle));
      struct gkyl_gen_miller_eqdsk_inp other = inp;
      other.a_shift = -0.4;
      other.Zvert_axis = 0.2;
      other.Rmajor_axis = 3;
      TEST_CHECK(theta == gkyl_gen_miller_eqdsk_theta_tok_at_theta(&other, angle));
      previous = theta;
    }
    TEST_CHECK(gkyl_gen_miller_eqdsk_theta_tok_at_theta(&inp, 0) == 0);
    TEST_CHECK(gkyl_gen_miller_eqdsk_theta_tok_at_theta(&inp, M_PI) == M_PI);
    TEST_CHECK(gkyl_gen_miller_eqdsk_theta_tok_at_theta(&inp, -M_PI) == -M_PI);
    double angle = M_PI - 5e-11;
    double theta = gkyl_gen_miller_eqdsk_theta_tok_at_theta(&inp, angle);
    TEST_CHECK(theta < M_PI);
    TEST_CHECK(fabs((M_PI - theta) / (M_PI - angle) - M_PI * inp.kappa / half_arc) < 2e-5);
    TEST_CHECK(
      fabs(
        gkyl_gen_miller_eqdsk_theta_tok_at_theta(&inp, 1e-15) / 1e-15 - M_PI * inp.kappa / half_arc
      ) < 1e-10
    );
    gkyl_gen_miller_eqdsk_release(up);
  }
  TEST_CHECK(miller_is_nan(gkyl_gen_miller_eqdsk_theta_tok_at_theta(NULL, 0)));
  double invalid_angles[] = {-4, 4, strtod("nan", NULL), strtod("inf", NULL)};
  for (int j = 0; j < 4; ++j) {
    TEST_CHECK(miller_is_nan(gkyl_gen_miller_eqdsk_theta_tok_at_theta(&inp, invalid_angles[j])));
  }
  double invalid_kappa[] = {-1, strtod("nan", NULL), strtod("inf", NULL)};
  for (int j = 0; j < 3; ++j) {
    inp.kappa = invalid_kappa[j];
    TEST_CHECK(miller_is_nan(gkyl_gen_miller_eqdsk_theta_tok_at_theta(&inp, 0.3)));
  }
  inp.kappa = 1;
  double invalid_delta[] = {-1, 1, strtod("nan", NULL), strtod("inf", NULL)};
  for (int j = 0; j < 4; ++j) {
    inp.delta = invalid_delta[j];
    TEST_CHECK(miller_is_nan(gkyl_gen_miller_eqdsk_theta_tok_at_theta(&inp, 0.3)));
  }
}

/**
 * Check parameter-only outboard maps, shared surface geometry, and invalid branches.
 *
 */
static void
test_outboard_maps(void)
{
  // These helpers need neither a field nor a safety-factor callback.
  struct gkyl_gen_miller_eqdsk_inp inp = {.Rmajor_axis = 2};
  double shifts[] = {-0.4, -1e-14, 0, 1e-14, 0.4};
  for (int j = 0; j < 5; ++j) {
    inp.a_shift = shifts[j];
    TEST_CHECK(gkyl_gen_miller_eqdsk_R_omp_at_rminor(&inp, 0) == 2);
    TEST_CHECK(gkyl_gen_miller_eqdsk_rminor_at_R_omp(&inp, 2) == 0);
    for (int k = 1; k <= 20; ++k) {
      double rminor = 0.03 * k;
      double R_omp = gkyl_gen_miller_eqdsk_R_omp_at_rminor(&inp, rminor);
      double expected = 2 + rminor - shifts[j] * rminor * rminor / 4;
      TEST_CHECK(fabs(R_omp - expected) < 1e-14);
      TEST_CHECK(fabs(gkyl_gen_miller_eqdsk_rminor_at_R_omp(&inp, R_omp) - rminor) < 1e-14);
    }
  }

  inp.a_shift = 0.4;
  inp.kappa = 1.5;
  inp.delta = -0.3;
  inp.bmag_axis = -1.5;
  inp.qsafety = constant_q;
  inp.Rmajor_LCFS = gkyl_gen_miller_eqdsk_R_omp_at_rminor(&inp, 0.5);
  gkyl_gen_miller_eqdsk *up = gkyl_gen_miller_eqdsk_inew(&inp);
  TEST_ASSERT(up != NULL);
  TEST_CHECK(fabs(gkyl_gen_miller_eqdsk_get_info(up)->rminor_LCFS - 0.5) < 1e-14);
  for (int k = 0; k <= 10; ++k) {
    double rminor = 0.05 * k, rz[2];
    gkyl_gen_miller_eqdsk_RZ_at_rtheta(up, rminor, 0, rz);
    TEST_CHECK(fabs(gkyl_gen_miller_eqdsk_R_omp_at_rminor(&inp, rminor) - rz[0]) < 1e-14);
  }
  gkyl_gen_miller_eqdsk_release(up);

  TEST_CHECK(miller_is_nan(gkyl_gen_miller_eqdsk_R_omp_at_rminor(NULL, 0.5)));
  TEST_CHECK(miller_is_nan(gkyl_gen_miller_eqdsk_rminor_at_R_omp(NULL, 2.5)));
  double bad_rminor[] = {-0.1, 5, 6, strtod("nan", NULL), strtod("inf", NULL)};
  double bad_R_omp[] = {1.9, 4.5, 5, strtod("nan", NULL), strtod("inf", NULL)};
  for (int j = 0; j < 5; ++j) {
    TEST_CHECK(miller_is_nan(gkyl_gen_miller_eqdsk_R_omp_at_rminor(&inp, bad_rminor[j])));
    TEST_CHECK(miller_is_nan(gkyl_gen_miller_eqdsk_rminor_at_R_omp(&inp, bad_R_omp[j])));
  }
  // Approach, but do not include, the increasing branch's turning point.
  double rnear = 4.99;
  double Rnear = gkyl_gen_miller_eqdsk_R_omp_at_rminor(&inp, rnear);
  TEST_CHECK(fabs(gkyl_gen_miller_eqdsk_rminor_at_R_omp(&inp, Rnear) - rnear) < 1e-10);

  double bad_axis[] = {0, -1, strtod("nan", NULL), strtod("inf", NULL)};
  for (int j = 0; j < 4; ++j) {
    inp.Rmajor_axis = bad_axis[j];
    TEST_CHECK(miller_is_nan(gkyl_gen_miller_eqdsk_R_omp_at_rminor(&inp, 0.5)));
    TEST_CHECK(miller_is_nan(gkyl_gen_miller_eqdsk_rminor_at_R_omp(&inp, 2.5)));
  }
  inp.Rmajor_axis = 2;
  double bad_shift[] = {strtod("nan", NULL), strtod("inf", NULL)};
  for (int j = 0; j < 2; ++j) {
    inp.a_shift = bad_shift[j];
    TEST_CHECK(miller_is_nan(gkyl_gen_miller_eqdsk_R_omp_at_rminor(&inp, 0.5)));
    TEST_CHECK(miller_is_nan(gkyl_gen_miller_eqdsk_rminor_at_R_omp(&inp, 2.5)));
  }
  inp.a_shift = -1e308;
  TEST_CHECK(miller_is_nan(gkyl_gen_miller_eqdsk_R_omp_at_rminor(&inp, 10)));
  TEST_CHECK(miller_is_nan(gkyl_gen_miller_eqdsk_rminor_at_R_omp(&inp, 10)));
}

/**
 * Check signed flux normalization independently of the Miller generator's orientation.
 *
 */
static void
test_flux_extent_conversion(void)
{
  double rho[2] = {0.5, 1.2}, psi[2];
  TEST_CHECK(gen_miller_flux_bounds(GKYL_GK_GEOMETRY_RADIAL_COORD_RHO_POL, rho, 3, 7, psi));
  TEST_CHECK(fabs(psi[0] - 4) < 1e-14 && fabs(psi[1] - 8.76) < 1e-14);
  TEST_CHECK(gen_miller_flux_bounds(GKYL_GK_GEOMETRY_RADIAL_COORD_RHO_POL, rho, 3, -1, psi));
  TEST_CHECK(fabs(psi[0] + 2.76) < 1e-14 && fabs(psi[1] - 2) < 1e-14);
  double flux[2] = {-2.76, 2};
  TEST_CHECK(gen_miller_flux_bounds(GKYL_GK_GEOMETRY_RADIAL_COORD_PSI_POL, flux, 3, -1, psi));
  TEST_CHECK(psi[0] == flux[0] && psi[1] == flux[1]);
  rho[1] = 1;
  TEST_CHECK(gen_miller_flux_bounds(GKYL_GK_GEOMETRY_RADIAL_COORD_RHO_POL, rho, 3, -1, psi));
  TEST_CHECK(psi[0] == -1);
  TEST_CHECK(!gen_miller_flux_bounds(GKYL_GK_GEOMETRY_RADIAL_COORD_RHO_POL, rho, 3, 3, psi));
  flux[1] = 3.1; // Wrong side of the axis for decreasing-outward flux.
  TEST_CHECK(!gen_miller_flux_bounds(GKYL_GK_GEOMETRY_RADIAL_COORD_PSI_POL, flux, 3, -1, psi));
  double invalid[][2] = {
    {0, 1},
    {-1, 1},
    {1, 1},
    {2, 1},
    {0.1, 1e308},
    {0.1, strtod("inf", NULL)},
    {strtod("nan", NULL), 1}
  };
  for (int j = 0; j < sizeof invalid / sizeof invalid[0]; ++j) {
    TEST_CHECK(!gen_miller_flux_bounds(GKYL_GK_GEOMETRY_RADIAL_COORD_RHO_POL, invalid[j], 3, 7, psi)
    );
  }
}

/**
 * Compare all extent representations for core, SOL, and LCFS-crossing domains.
 *
 */
static void
test_domain_extents(void)
{
  for (int shape = 0; shape < 2; ++shape) {
    for (int region = 0; region < 3; ++region) {
      struct gkyl_gen_miller_eqdsk_inp inp = {
        .Rmajor_axis = 2,
        .bmag_axis = shape ? -1.5 : 1.5,
        .kappa = shape ? 1.4 : 1,
        .delta = shape ? -0.3 : 0,
        .a_shift = shape ? 0.3 : 0,
        .qsafety = constant_q,
        .guard_fraction = 0.1,
        .rtol = 1e-9,
      };
      inp.Rmajor_LCFS = gkyl_gen_miller_eqdsk_R_omp_at_rminor(&inp, 0.5);
      double radii[3][2] = {{0.2, 0.4}, {0.6, 0.65}, {0.4, 0.6}};
      inp.has_Rmajor_omp_bounds = true;
      for (int j = 0; j < 2; ++j) {
        inp.Rmajor_omp_bounds[j] = gkyl_gen_miller_eqdsk_R_omp_at_rminor(&inp, radii[region][j]);
      }
      gkyl_gen_miller_eqdsk *reference = gkyl_gen_miller_eqdsk_inew(&inp);
      TEST_ASSERT(reference != NULL);
      double expected[2], rho[2];
      char boundary[32];
      snprintf(
        boundary, sizeof boundary, "%.9E", gkyl_gen_miller_eqdsk_get_info(reference)->psi_boundary
      );
      for (int j = 0; j < 2; ++j) {
        expected[j] = gkyl_gen_miller_eqdsk_psi_at_R_omp(reference, inp.Rmajor_omp_bounds[j]);
        rho[j] = sqrt(expected[j] / strtod(boundary, NULL));
      }
      gkyl_gen_miller_eqdsk_release(reference);
      for (int rep = 0; rep < 3; ++rep) {
        inp.x_extent_type = rep;
        inp.has_Rmajor_omp_bounds = false; // Exercise automatic coverage.
        const double *bounds = rep == 0 ? expected : rep == 1 ? rho : inp.Rmajor_omp_bounds;
        double snapshot[2] = {bounds[0], bounds[1]}, resolved[2];
        gkyl_gen_miller_eqdsk *up = gkyl_gen_miller_eqdsk_inew_for_domain(&inp, bounds, resolved);
        TEST_ASSERT(up != NULL);
        const struct gkyl_gen_miller_eqdsk_info *info = gkyl_gen_miller_eqdsk_get_info(up);
        TEST_CHECK(bounds[0] == snapshot[0] && bounds[1] == snapshot[1]);
        TEST_CHECK(!inp.has_Rmajor_omp_bounds);
        for (int j = 0; j < 2; ++j) {
          TEST_CHECK(fabs(resolved[j] - expected[j]) < 1e-8 * info->psi_boundary);
          if (rep == 0) {
            TEST_CHECK(resolved[j] == expected[j]);
          }
          TEST_CHECK(
            fabs(gkyl_gen_miller_eqdsk_r_at_psi(up, resolved[j]) - radii[region][j]) < 1e-8
          );
        }
        TEST_CHECK(info->rminor_required >= radii[region][1] - 1e-8);
        TEST_CHECK(info->psi_guard > resolved[1]);
        if (shape == 0 && region == 1 && rep == 1) {
          TEST_CHECK(gkyl_gen_miller_eqdsk_advance(up));
        }
        gkyl_gen_miller_eqdsk_release(up);
      }
    }
  }
}

/**
 * Check explicit coverage, serialized LCFS endpoints, and domain-constructor failures.
 *
 */
static void
test_domain_invalid(void)
{
  struct gkyl_gen_miller_eqdsk_inp inp = {
    .Rmajor_axis = 2,
    .bmag_axis = 1.5,
    .Rmajor_LCFS = 2.5,
    .qsafety = constant_q,
    .x_extent_type = GKYL_GK_GEOMETRY_RADIAL_COORD_RHO_POL,
    .has_Rmajor_omp_bounds = true,
    .Rmajor_omp_bounds = {2.1, 2.5},
  };
  double bounds[2] = {0.5, 1}, resolved[2];
  gkyl_gen_miller_eqdsk *up = gkyl_gen_miller_eqdsk_inew_for_domain(&inp, bounds, resolved);
  TEST_ASSERT(up != NULL);
  char boundary[32];
  snprintf(boundary, sizeof boundary, "%.9E", gkyl_gen_miller_eqdsk_get_info(up)->psi_boundary);
  TEST_CHECK(resolved[1] == strtod(boundary, NULL));
  // Check the real writer and reader agree with rho=1 exactly.
  char filename[] = "/tmp/gkyl_miller_extent_XXXXXX";
  int fd = mkstemp(filename);
  TEST_ASSERT(fd >= 0);
  close(fd);
  TEST_ASSERT(gkyl_gen_miller_eqdsk_advance(up));
  TEST_ASSERT(gkyl_gen_miller_eqdsk_write(up, filename));
  struct gkyl_efit_inp einp = {.rz_poly_order = 2, .flux_poly_order = 1};
  snprintf(einp.filepath, sizeof einp.filepath, "%s", filename);
  struct gkyl_efit *efit = gkyl_efit_new(&einp);
  TEST_CHECK(resolved[1] == efit->sibry);
  gkyl_efit_release(efit);
  unlink(filename);
  gkyl_gen_miller_eqdsk_release(up);
  bounds[1] = 1.01; // Explicit coverage must not be silently expanded into the guard.
  TEST_CHECK(gkyl_gen_miller_eqdsk_inew_for_domain(&inp, bounds, resolved) == NULL);
  bounds[0] = 0.01;
  bounds[1] = 0.5;
  TEST_CHECK(gkyl_gen_miller_eqdsk_inew_for_domain(&inp, bounds, resolved) == NULL);
  inp.has_Rmajor_omp_bounds = false;
  double invalid[][2] = {
    {0, 1},      {-0.1, 1}, {1, 1}, {2, 1}, {strtod("nan", NULL), 1}, {0.5, strtod("inf", NULL)},
    {0.5, 1e308}
  };
  for (int j = 0; j < sizeof invalid / sizeof invalid[0]; ++j) {
    TEST_CHECK(gkyl_gen_miller_eqdsk_inew_for_domain(&inp, invalid[j], resolved) == NULL);
  }
  bounds[0] = 0.5;
  bounds[1] = 100;
  TEST_CHECK(gkyl_gen_miller_eqdsk_inew_for_domain(&inp, bounds, resolved) == NULL);
  bounds[1] = 0.9;
  TEST_CHECK(gkyl_gen_miller_eqdsk_inew_for_domain(NULL, bounds, resolved) == NULL);
  TEST_CHECK(gkyl_gen_miller_eqdsk_inew_for_domain(&inp, NULL, resolved) == NULL);
  TEST_CHECK(gkyl_gen_miller_eqdsk_inew_for_domain(&inp, bounds, NULL) == NULL);
  inp.x_extent_type = 100;
  TEST_CHECK(gkyl_gen_miller_eqdsk_inew_for_domain(&inp, bounds, resolved) == NULL);
  inp.x_extent_type = GKYL_GK_GEOMETRY_RADIAL_COORD_R_OMP;
  bounds[0] = 2;
  bounds[1] = 2.4;
  TEST_CHECK(gkyl_gen_miller_eqdsk_inew_for_domain(&inp, bounds, resolved) == NULL);
  inp.x_extent_type = GKYL_GK_GEOMETRY_RADIAL_COORD_PSI_POL;
  bounds[0] = -0.1;
  bounds[1] = 0.1;
  TEST_CHECK(gkyl_gen_miller_eqdsk_inew_for_domain(&inp, bounds, resolved) == NULL);
}

TEST_LIST = {
  {"miller_flux_extent_conversion", test_flux_extent_conversion},
  {"miller_domain_extents", test_domain_extents},
  {"miller_domain_invalid", test_domain_invalid},
  {"miller_circle_and_invalid_inputs", test_circle},
  {"miller_devices_and_reader", test_devices},
  {"miller_coordinate_conversions", test_coordinates},
  {"miller_forward_angle", test_forward_angle},
  {"miller_outboard_maps", test_outboard_maps},
  {NULL, NULL}
};
