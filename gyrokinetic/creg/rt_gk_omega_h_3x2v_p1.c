#include <math.h>
#include <stdio.h>

#include <gkyl_array_ops.h>
#include <gkyl_gyrokinetic_priv.h>

/**
 * Constant shear map: J=1, g^xx=1+s^2, g^xy=-s, g^yy=g_zz=1.
 *
 * @param t Time (unused).
 * @param xc Computational coordinates.
 * @param xp Output Cartesian coordinates.
 * @param ctx Pointer to the constant shear.
 */
static void
mapc2p(double t, const double *xc, double *xp, void *ctx)
{
  double shear = *(double *)ctx;
  xp[0] = xc[0] + shear * xc[1];
  xp[1] = xc[1];
  xp[2] = xc[2];
}

/**
 * Uniform Cartesian magnetic field parallel to the third coordinate.
 *
 * @param t Time (unused).
 * @param xc Computational coordinates (unused).
 * @param out Output Cartesian magnetic-field components.
 * @param ctx Context (unused).
 */
static void
bfield(double t, const double *xc, double *out, void *ctx)
{
  out[0] = out[1] = 0.0;
  out[2] = 2.0;
}

/**
 * Constant density and temperature for the species projection.
 *
 * @param t Time (unused).
 * @param xc Computational coordinates (unused).
 * @param out Output unit density or temperature.
 * @param ctx Context (unused).
 */
static void
unit_profile(double t, const double *xc, double *out, void *ctx)
{
  out[0] = 1.0;
}

/**
 * Zero parallel flow for the species projection.
 *
 * @param t Time (unused).
 * @param xc Computational coordinates (unused).
 * @param out Output parallel velocity.
 * @param ctx Context (unused).
 */
static void
zero_profile(double t, const double *xc, double *out, void *ctx)
{
  out[0] = 0.0;
}

/**
 * Check the app's cached factor against the analytic symmetric contraction.
 *
 * @param shear Constant coordinate shear.
 * @return Zero on success, one on failure.
 */
static int
check_shear(double shear)
{
  struct gkyl_gk inp = {
    .name = "gk_omega_h_3x2v_p1",
    .cdim = 3,
    .lower = {-1.0, -1.0, -1.0},
    .upper = {1.0, 1.0, 1.0},
    .cells = {2, 2, 2},
    .poly_order = 1,
    .basis_type = GKYL_BASIS_MODAL_SERENDIPITY,
    .geometry =
      {
        .geometry_id = GKYL_GEOMETRY_MAPC2P,
        .mapc2p = mapc2p,
        .c2p_ctx = &shear,
        .bfield_func = bfield,
      },
    .num_periodic_dir = 2,
    .periodic_dirs = {1, 2},
    .num_species = 1,
    .species = {{
      .name = "ion",
      .charge = 1.0,
      .mass = 1.0,
      .vdim = 2,
      .lower = {-4.0, 0.0},
      .upper = {4.0, 4.0},
      .cells = {4, 4},
      .polarization_density = 1.0,
      .projection =
        {
          .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
          .density = unit_profile,
          .temp = unit_profile,
          .upar = zero_profile,
        },
      .collisionless = {.type = GKYL_GK_COLLISIONLESS_ES},
      .bcs =
        {{.dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX},
         {.dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX}},
    }},
    .field =
      {
        .polarization_bmag = 2.0,
        .poisson_bcs =
          {{.dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET},
           {.dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET},
           {.dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_PERIODIC},
           {.dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_PERIODIC}},
      },
  };
  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&inp);
  struct gkyl_array *density =
    gkyl_array_new(GKYL_DOUBLE, app->basis.num_basis, app->local_ext.volume);
  // Jn=1; the normalized constant basis function in 3x is 1/sqrt(8).
  gkyl_array_shiftc(density, sqrt(8.0), 0);
  double measured = gkyl_gyrokinetic_omega_h_advance(app->omegaH, density);

  // C/(JB)=1, polarization weight=n_pol*m/B_pol^2=1/4, kx=ky=pi/2.
  double kperp = M_PI / 2.0, kpar = 2.0 * M_PI;
  double expected =
    kpar * kpar / (0.25 * kperp * kperp * ((1.0 + shear * shear) - 2.0 * shear + 1.0));
  double rel_error = fabs(measured / expected - 1.0);
  printf(
    "shear=%g: omegaH_geo_fac^2=%g, expected=%g, relative error=%.3g\n", shear, measured, expected,
    rel_error
  );
  gkyl_array_release(density);
  gkyl_gyrokinetic_app_release(app);
  return !(rel_error < 1e-9);
}

/**
 * Cover both signs of the mixed term and the zero-shear limit on the CPU.
 *
 * @return Zero if all checks pass, one otherwise.
 */
int
main(void)
{
  int failures = check_shear(0.0);
  failures += check_shear(0.75);
  failures += check_shear(-0.75);
  return failures ? 1 : 0;
}
