#include <acutest.h>
#include <gkyl_gk_geometry_tok.h>
#include <gkyl_gk_geometry_priv.h>
#include <gkyl_position_map.h>
#include <gkyl_tok_geo_priv.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

struct tok_field_ctx {
  double psi_sign, tor_sign, shear;
};

/**
 * Check finiteness without relying on floating-point predicates disabled by fast-math.
 *
 * @param value Floating-point value to classify.
 * @return Whether the value is finite.
 */
static bool
tok_is_finite(double value)
{
  uint64_t bits;
  memcpy(&bits, &value, sizeof bits);
  volatile uint64_t exponent = bits & UINT64_C(0x7ff0000000000000);
  return exponent != UINT64_C(0x7ff0000000000000);
}

/**
 * Write nested analytic flux surfaces and a linear F profile, constant in the SOL.
 *
 * @param stream Temporary GEQDSK stream owned by the caller.
 * @param ctx Flux orientation, toroidal orientation, and F slope.
 */
static void
tok_write_equilibrium(FILE *stream, const struct tok_field_ctx *ctx)
{
  const int count = 65;
  double boundary = ctx->psi_sign * 0.125;
  double header[] = {2, 2, 2, 1, 0, 2, 0, 0,        boundary, ctx->tor_sign * 1.5,
                     0, 0, 0, 2, 0, 0, 0, boundary, 0,        0};
  fprintf(stream, "%-48s%4d%4d%4d\n", "Tokamak field and curl reference", 0, count, count);
  for (int i = 0; i < 20; ++i) {
    fprintf(stream, "%16.9E%s", header[i], i % 5 == 4 ? "\n" : "");
  }
  // FPOL, pressure, FF', and pressure derivative, in axis-to-boundary order.
  for (int profile = 0; profile < 4; ++profile) {
    for (int i = 0; i < count; ++i) {
      double psi = boundary * i / (count - 1);
      double fpol = ctx->tor_sign * (3 + ctx->shear * psi);
      double value = profile == 0 ? fpol : profile == 2 ? fpol * ctx->tor_sign * ctx->shear : 0;
      fprintf(stream, "%16.9E%s", value, i % 5 == 4 ? "\n" : "");
    }
  }
  for (int iz = 0; iz < count; ++iz) {
    for (int ir = 0; ir < count; ++ir) {
      double Rmajor = 1 + 2.0 * ir / (count - 1), Zvert = -1 + 2.0 * iz / (count - 1);
      double offset = Rmajor - 2;
      double psi =
        ctx->psi_sign * (0.5 * (offset * offset + Zvert * Zvert) + 0.05 * offset * offset * offset);
      fprintf(stream, "%16.9E%s", psi, (iz * count + ir) % 5 == 4 ? "\n" : "");
    }
  }
  fprintf(stream, "\n");
  for (int i = 0; i < count; ++i) {
    fprintf(stream, "%16.9E%s", 2.0, i % 5 == 4 ? "\n" : "");
  }
  fprintf(stream, "\n");
  fprintf(stream, "%5d%5d\n", count, 0);
  for (int i = 0; i < count; ++i) {
    double theta = 2 * M_PI * i / (count - 1);
    double lo = 0, hi = 1;
    for (int j = 0; j < 52; ++j) {
      double rminor = (lo + hi) / 2, offset = rminor * cos(theta);
      if (0.5 * rminor * rminor + 0.05 * offset * offset * offset < 0.125) {
        lo = rminor;
      } else {
        hi = rminor;
      }
    }
    double rminor = (lo + hi) / 2;
    fprintf(stream, "%16.9E%16.9E\n", 2 + rminor * cos(theta), rminor * sin(theta));
  }
}

/**
 * Evaluate the analytic Cartesian unit field using B = grad(psi) cross grad(phi) + F grad(phi).
 *
 * @param ctx Signed equilibrium parameters.
 * @param position Cartesian position [m].
 * @param bhat Output Cartesian unit field.
 */
static void
tok_reference_bhat(const struct tok_field_ctx *ctx, const double position[3], double bhat[3])
{
  double Rmajor = hypot(position[0], position[1]), Zvert = position[2];
  double offset = Rmajor - 2;
  double flux = 0.5 * (offset * offset + Zvert * Zvert) + 0.05 * offset * offset * offset;
  double psi = ctx->psi_sign * fmin(flux, 0.125);
  double fpol = ctx->tor_sign * (3 + ctx->shear * psi);
  double Br = -ctx->psi_sign * Zvert / Rmajor;
  double Bz = ctx->psi_sign * (offset + 0.15 * offset * offset) / Rmajor;
  double Bphi = fpol / Rmajor, bmag = sqrt(Br * Br + Bphi * Bphi + Bz * Bz);
  bhat[0] = (Br * position[0] - Bphi * position[1]) / (Rmajor * bmag);
  bhat[1] = (Br * position[1] + Bphi * position[0]) / (Rmajor * bmag);
  bhat[2] = Bz / bmag;
}

/**
 * Compare nodal unit fields and their curl with independent Cartesian finite differences.
 *
 * @param ctx Signed equilibrium parameters.
 * @param range Nodal range for interior or face geometry.
 * @param position Stored cylindrical positions in (R,Z,phi) order.
 * @param bhat Stored Cartesian unit field.
 * @param curl Stored Cartesian curl of the unit field.
 * @param errors Output running maximum errors in bhat, core curl, and SOL curl.
 */
static void
tok_check_field(
  const struct tok_field_ctx *ctx, const struct gkyl_range *range,
  const struct gkyl_array *position, const struct gkyl_array *bhat, const struct gkyl_array *curl,
  double errors[3]
)
{
  struct gkyl_range_iter iter;
  int samples[2] = {0};
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);
    const double *rzphi = gkyl_array_cfetch(position, loc);
    const double *actual_b = gkyl_array_cfetch(bhat, loc);
    const double *actual_curl = gkyl_array_cfetch(curl, loc);
    for (int comp = 0; comp < 3; ++comp) {
      TEST_ASSERT(tok_is_finite(rzphi[comp]));
      TEST_ASSERT(tok_is_finite(actual_b[comp]));
      TEST_ASSERT(tok_is_finite(actual_curl[comp]));
    }
    TEST_ASSERT(rzphi[0] > 0);
    double xyz[] = {rzphi[0] * cos(rzphi[2]), rzphi[0] * sin(rzphi[2]), rzphi[1]};
    double expected_b[3], grad_b[3][3];
    tok_reference_bhat(ctx, xyz, expected_b);
    double offset = rzphi[0] - 2;
    double flux = 0.5 * (offset * offset + rzphi[1] * rzphi[1]) + 0.05 * offset * offset * offset;
    // The piecewise F profile has a derivative jump at the LCFS; omit that trace.
    if (fabs(flux - 0.125) < 1e-5) {
      continue;
    }
    samples[flux < 0.125 ? 0 : 1] += 1;
    const double step = 1e-5;
    for (int dir = 0; dir < 3; ++dir) {
      double plus[] = {xyz[0], xyz[1], xyz[2]}, minus[] = {xyz[0], xyz[1], xyz[2]};
      double bplus[3], bminus[3];
      plus[dir] += step;
      minus[dir] -= step;
      tok_reference_bhat(ctx, plus, bplus);
      tok_reference_bhat(ctx, minus, bminus);
      for (int comp = 0; comp < 3; ++comp) {
        grad_b[dir][comp] = (bplus[comp] - bminus[comp]) / (2 * step);
      }
    }
    double expected_curl[] = {
      grad_b[1][2] - grad_b[2][1], grad_b[2][0] - grad_b[0][2], grad_b[0][1] - grad_b[1][0]
    };
    for (int comp = 0; comp < 3; ++comp) {
      errors[0] = fmax(errors[0], fabs(actual_b[comp] - expected_b[comp]));
      int region = flux < 0.125 ? 1 : 2;
      errors[region] = fmax(errors[region], fabs(actual_curl[comp] - expected_curl[comp]));
    }
  }
  TEST_CHECK(samples[0] > 0 && samples[1] > 0);
}

// Independent polar-angle quadrature for the synthetic cubic flux (not a traced contour).
static double
tok_reference_q(const struct tok_field_ctx *ctx, double psi)
{
  double fpol = ctx->tor_sign * (3 + ctx->shear * ctx->psi_sign * fmin(ctx->psi_sign * psi, 0.125));
  double integral = 0;
  const int nodes = 4096;
  for (int i = 0; i < nodes; ++i) {
    double angle = 2 * M_PI * (i + 0.5) / nodes, cosine = cos(angle),
           cube = cosine * cosine * cosine;
    double radius = sqrt(2 * ctx->psi_sign * psi);
    for (int j = 0; j < 8; ++j) {
      radius -=
        (0.5 * radius * radius + 0.05 * radius * radius * radius * cube - ctx->psi_sign * psi) /
        (radius + 0.15 * radius * radius * cube);
    }
    integral += 1 / ((2 + radius * cosine) * (1 + 0.15 * radius * cube));
  }
  return ctx->psi_sign * fpol * integral / nodes;
}

// Check the signed full-surface q despite the partial simulation theta interval.
static void
tok_check_q(const struct tok_field_ctx *ctx, const struct gk_geometry *geo)
{
  struct gkyl_range_iter iter;
  double error = 0;
  gkyl_range_iter_init(&iter, &geo->local);
  while (gkyl_range_iter_next(&iter)) {
    double center[3];
    gkyl_rect_grid_cell_center(&geo->grid, iter.idx, center);
    const double *modal =
      gkyl_array_cfetch(geo->geo_int.qprofile, gkyl_range_idx(&geo->local, iter.idx));
    for (int side = -1; side <= 1; side += 2) {
      double eta[] = {side / sqrt(3.0), 0, 0};
      double psi = center[0] + eta[0] * geo->grid.dx[0] / 2;
      double expected = tok_reference_q(ctx, psi);
      double actual = geo->basis.eval_expand(eta, modal);
      // Without a limiter cut the full IWL connection covers the whole contour,
      // including outside the nominal LCFS of this analytic equilibrium.
      TEST_CHECK(tok_is_finite(actual));
      error = fmax(error, fabs(actual - expected));
    }
  }
  TEST_CHECK(error < 1e-3);
  TEST_MSG("signed q error=%.9g psi_sign=%g F_sign=%g", error, ctx->psi_sign, ctx->tor_sign);
}

// Shift a nonuniform parallel map to differentiate the actual geometry map.
static void
tok_shift_theta(double t, const double *xn, double *out, void *ctx)
{
  double *params = ctx;
  out[0] = xn[0] + params[0] + params[1] * sin(xn[0]);
}

// Compare independently shifted position maps with the supplied contour derivatives.
static void
tok_check_map(struct gk_geometry *geometry, struct gkyl_gk_geometry_inp *inp)
{
  inp->tok_grid_info.cgrid = inp->grid;
  inp->tok_grid_info.cbasis = inp->basis;
  struct gkyl_tok_geo *tok = gkyl_tok_geo_new(&inp->efit_info, &inp->tok_grid_info);
  struct gkyl_array *deriv = gkyl_array_clone(geometry->geo_int.ddtheta_nodal);
  struct gkyl_array *plus = gkyl_array_clone(geometry->geo_int.mc2p_nodal);
  bool original_derivs = inp->position_map->use_map_derivs;
  inp->position_map->use_map_derivs = false;
  mc2nu_t original = inp->position_map->maps[2];
  void *original_ctx = inp->position_map->ctxs[2];
  double params[] = {0, 0.07};
  inp->position_map->maps[2] = tok_shift_theta;
  inp->position_map->ctxs[2] = params;
  gkyl_tok_geo_calc_interior(
    geometry, &geometry->nrange_int, geometry->dzc, tok, &inp->tok_grid_info, inp->position_map
  );
  gkyl_array_copy(deriv, geometry->geo_int.ddtheta_nodal);
  const double step = 1e-4;
  params[0] = step;
  gkyl_tok_geo_calc_interior(
    geometry, &geometry->nrange_int, geometry->dzc, tok, &inp->tok_grid_info, inp->position_map
  );
  gkyl_array_copy(plus, geometry->geo_int.mc2p_nodal);
  params[0] = -step;
  gkyl_tok_geo_calc_interior(
    geometry, &geometry->nrange_int, geometry->dzc, tok, &inp->tok_grid_info, inp->position_map
  );
  double error[3] = {0};
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &geometry->nrange_int);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&geometry->nrange_int, iter.idx);
    const double *upper = gkyl_array_cfetch(plus, loc);
    const double *lower = gkyl_array_cfetch(geometry->geo_int.mc2p_nodal, loc);
    const double *expected = gkyl_array_cfetch(deriv, loc);
    double theta = gk_geometry_node_coord(
      &geometry->grid, &geometry->local, &geometry->global, &geometry->nrange_int, 2, iter.idx[2],
      geometry->basis.poly_order + 1, true
    );
    double slope = 1 + params[1] * cos(theta);
    for (int c = 0; c < 3; ++c) {
      error[c] = fmax(error[c], fabs((upper[c] - lower[c]) / (2 * step) * slope - expected[c]));
    }
  }
  for (int c = 0; c < 3; ++c) {
    TEST_CHECK(error[c] < 4e-4);
    TEST_MSG(
      "map derivative comp=%d flux orientation=%d error=%.9g", c, tok->theta_orientation, error[c]
    );
  }
  inp->position_map->use_map_derivs = original_derivs;
  inp->position_map->maps[2] = original;
  inp->position_map->ctxs[2] = original_ctx;
  gkyl_array_release(plus);
  gkyl_array_release(deriv);
  gkyl_tok_geo_release(tok);
}

/**
 * Check both field signs and flux orientations at interior and face nodes, in core and SOL.
 *
 * @param cubic Whether to use the cubic rather than quadratic EQDSK reconstruction.
 */
static void
tok_check_curl(bool cubic)
{
  for (int psi_sign = 1; psi_sign >= -1; psi_sign -= 2) {
    for (int tor_sign = -1; tor_sign <= 1; tor_sign += 2) {
      for (int varying = 0; varying < 2; ++varying) {
        struct tok_field_ctx ctx = {psi_sign, tor_sign, varying ? 0.4 : 0};
        char filename[] = "/tmp/gkyl_tok_curl_XXXXXX";
        int descriptor = mkstemp(filename);
        TEST_ASSERT(descriptor >= 0);
        FILE *stream = fdopen(descriptor, "w");
        TEST_ASSERT(stream != NULL);
        tok_write_equilibrium(stream, &ctx);
        TEST_ASSERT(fclose(stream) == 0);

        struct gkyl_gk_geometry_inp inp = {
          .geometry_id = GKYL_GEOMETRY_TOKAMAK,
          .efit_info = {.rz_poly_order = 2, .flux_poly_order = 1},
          .tok_grid_info =
            {
              .ftype = GKYL_GEOMETRY_TOKAMAK_IWL,
              .rleft = 1,
              .rright = 3,
              .rmin = 1,
              .rmax = 3,
              .rclose = 2.5,
              .zmin = -1,
              .zmax = 1,
              .use_cubics = cubic,
              .use_hyperbolic_numbers = true,
            },
          .position_map = gkyl_position_map_null_new(),
        };
        snprintf(inp.efit_info.filepath, sizeof inp.efit_info.filepath, "%s", filename);
        double lower[] = {psi_sign > 0 ? 0.08 : -0.17, -0.2, -2.8};
        double upper[] = {psi_sign > 0 ? 0.17 : -0.08, 0.2, 2.8};
        gkyl_rect_grid_init(&inp.grid, 3, lower, upper, (int[]){2, 1, 2});
        gkyl_cart_modal_serendip(&inp.basis, 3, 1);
        gkyl_create_grid_ranges(&inp.grid, (int[]){1, 1, 1}, &inp.local_ext, &inp.local);
        inp.global = inp.geo_local = inp.geo_global = inp.local;
        inp.global_ext = inp.geo_local_ext = inp.geo_global_ext = inp.local_ext;
        inp.geo_grid = inp.grid;
        inp.geo_basis = inp.basis;
        struct gk_geometry *geo = gkyl_gk_geometry_tok_new(&inp);
        tok_check_q(&ctx, geo);
        double errors[3] = {0};
        tok_check_field(
          &ctx, &geo->nrange_int, geo->geo_int.mc2p_nodal, geo->geo_int.bcart_nodal,
          geo->geo_int.curlbhat_nodal, errors
        );
        for (int dir = 0; dir < 3; ++dir) {
          tok_check_field(
            &ctx, &geo->nrange_surf[dir], geo->geo_surf[dir].mc2p_nodal_fd,
            geo->geo_surf[dir].bcart_nodal, geo->geo_surf[dir].curlbhat_nodal, errors
          );
        }
        for (int quantity = 0; quantity < 3; ++quantity) {
          // Both reconstructions have second-derivative error from sampled flux;
          // doubling the R-Z resolution halves the observed curl error.
          TEST_CHECK(errors[quantity] < (quantity == 0 ? 5e-5 : 4e-3));
          TEST_MSG(
            "cubic=%d psi_sign=%d tor_sign=%d varying_F=%d quantity=%d error=%.9g", cubic, psi_sign,
            tor_sign, varying, quantity, errors[quantity]
          );
        }
        if (!varying) {
          tok_check_map(geo, &inp);
        }
        gkyl_gk_geometry_release(geo);
        gkyl_position_map_release(inp.position_map);
        TEST_CHECK(unlink(filename) == 0);
      }
    }
  }
}

/**
 * Check the quadratic EQDSK field and curl reconstruction.
 *
 */
static void
test_tok_curl_quadratic(void)
{
  tok_check_curl(false);
}

/**
 * Check the cubic EQDSK field and curl reconstruction.
 *
 */
static void
test_tok_curl_cubic(void)
{
  tok_check_curl(true);
}

/**
 * Check that the quadratic limit of the cubic root solver returns signed contour tangents.
 *
 */
static void
test_tok_cubic_quadratic_limit(void)
{
  struct gkyl_basis basis;
  gkyl_cart_modal_tensor(&basis, 2, 3);
  double nodes[] = {-sqrt(3.0 / 5), 0, sqrt(3.0 / 5)}, weights[] = {5.0 / 9, 8.0 / 9, 5.0 / 9};
  for (int sign = -1; sign <= 1; sign += 2) {
    double modal[16] = {0};
    // Three-point Gauss quadrature exactly projects this quadratic onto a p3 basis.
    for (int ir = 0; ir < 3; ++ir) {
      for (int iz = 0; iz < 3; ++iz) {
        double point[] = {nodes[ir], nodes[iz]}, values[16];
        basis.eval(point, values);
        double flux = sign * 0.5 * ((point[0] - 0.1) * (point[0] - 0.1) + point[1] * point[1]);
        for (int mode = 0; mode < 16; ++mode) {
          modal[mode] += weights[ir] * weights[iz] * flux * values[mode];
        }
      }
    }
    double center[] = {2, 0}, spacing[] = {2, 2}, Zvert = 0.2;
    struct RdRdZ_sol algorithms[] = {
      calc_RdR_p3_hyperbolic(modal, sign * 0.125, Zvert, center, spacing),
      calc_RdR_p3(modal, sign * 0.125, Zvert, center, spacing)
    };
    for (int method = 0; method < 2; ++method) {
      struct RdRdZ_sol roots = algorithms[method];
      TEST_ASSERT(roots.nsol == 2);
      for (int i = 0; i < roots.nsol; ++i) {
        double offset = roots.R[i] - 2.1;
        TEST_CHECK(fabs(offset * offset + Zvert * Zvert - 0.25) < 1e-12);
        TEST_CHECK(fabs(roots.dR[i] + 2 * sign * Zvert) < 1e-12);
        TEST_CHECK(fabs(roots.dZ[i] - 2 * sign * offset) < 1e-12);
      }
    }
  }
}

// Compare core q with EQDSK and check finite, identical full-connection SOL
// profiles for the unsplit contour and each of its three block types.
static void
test_tok_qprofile_branches(void)
{
  struct gkyl_efit_inp einp = {
    .filepath = "gyrokinetic/data/eqdsk/asdex.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
  };

  struct gkyl_tok_geo_grid_inp inp = {
    .ftype = GKYL_GEOMETRY_TOKAMAK_CORE,
    .rmin = 0,
    .rmax = 5,
    .rclose = 2.5,
    .rright = 2.5,
    .rleft = 0.7,
    .zmin = -1.3,
    .zmax = 1.0,
    .zmin_left = -1.2,
    .zmin_right = -1.0,
  };
  struct gkyl_tok_geo *geo = gkyl_tok_geo_new(&einp, &inp);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 1, 1);

  double lower[] = {0.01}, upper[] = {0.08};
  int cells[] = {4};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);

  struct gkyl_array *profile = gkyl_array_new(GKYL_DOUBLE, 2, cells[0]);
  struct gkyl_array *split = gkyl_array_new(GKYL_DOUBLE, 2, cells[0]);

  gkyl_tok_geo_calc_qprofile(geo, &inp, &grid, &basis, profile);

  double max_error = 0;
  for (int cell = 0; cell < cells[0]; ++cell) {
    for (int side = -1; side <= 1; side += 2) {
      double eta = side / sqrt(3.0);
      double psi = grid.lower[0] + (cell + 0.5 * (1 + eta)) * grid.dx[0];
      int idx = geo->frange.lower[0] + (int)((psi - geo->fgrid.lower[0]) / geo->fgrid.dx[0]);

      double center;
      gkyl_rect_grid_cell_center(&geo->fgrid, &idx, &center);

      double efeta = 2 * (psi - center) / geo->fgrid.dx[0];
      double expected = geo->fbasis.eval_expand(
        &efeta, gkyl_array_cfetch(geo->qdg, gkyl_range_idx(&geo->frange, &idx))
      );

      double actual = basis.eval_expand(&eta, gkyl_array_cfetch(profile, cell));
      TEST_CHECK(tok_is_finite(actual));
      max_error = fmax(max_error, fabs(actual / expected - 1));
    }
  }
  TEST_CHECK(max_error < 0.02);
  TEST_MSG("core q relative error against EQDSK = %.9g", max_error);

  // A cell face lies at the separatrix, while the Gauss nodes remain in the SOL.
  lower[0] = geo->psisep;
  upper[0] = lower[0] + 0.01;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);

  inp.ftype = GKYL_GEOMETRY_TOKAMAK_LSN_SOL;
  gkyl_tok_geo_calc_qprofile(geo, &inp, &grid, &basis, profile);

  enum gkyl_tok_geo_type types[] = {
    GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO, GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID,
    GKYL_GEOMETRY_TOKAMAK_LSN_SOL_UP
  };
  for (int part = 0; part < 3; ++part) {
    inp.ftype = types[part];
    gkyl_tok_geo_calc_qprofile(geo, &inp, &grid, &basis, split);
    for (int cell = 0; cell < cells[0]; ++cell) {
      const double *full = gkyl_array_cfetch(profile, cell);
      const double *piece = gkyl_array_cfetch(split, cell);
      for (int mode = 0; mode < 2; ++mode) {
        TEST_CHECK(tok_is_finite(piece[mode]));
        TEST_CHECK(fabs(piece[mode] - full[mode]) < 1e-12);
      }
    }
  }

  gkyl_array_release(split);
  gkyl_array_release(profile);
  gkyl_tok_geo_release(geo);
}

TEST_LIST = {
  {"tok_qprofile_branches", test_tok_qprofile_branches},
  {"tok_cubic_quadratic_limit", test_tok_cubic_quadratic_limit},
  {"tok_curl_quadratic", test_tok_curl_quadratic},
  {"tok_curl_cubic", test_tok_curl_cubic},
  {NULL, NULL}
};
