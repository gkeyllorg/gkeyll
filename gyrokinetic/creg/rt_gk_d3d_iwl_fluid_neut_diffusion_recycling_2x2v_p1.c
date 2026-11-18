// DIII-D inner-wall-limited 2x2v regression for a locally recycled,
// diffusive fluid neutral.
//
// The configuration coordinates are x (radial) and theta (parallel). The
// domain contains only the open-field-line region, so both parallel
// boundaries are material targets. The mapped neutral equation is
//
//   d_t(J n_n) = d_i[J D_n g_n^{ij} d_j n_n] - J nu_iz n_n,
//   D_n = v_ti^2/(n_e <sigma v>_cx).
//
// At each target the diffusive neutral flux is replaced by the recycled flux,
// so the neutral flux entering through the target is locally
//
//   Gamma_n = R Gamma_i,
//
// where R is the recycling fraction and Gamma_i is the incident ion flux.
// The test exercises reaction-rate diffusion, ionization loss, recycling,
// mapped geometry, and the matrix-free backward-Euler diffusion solve.

#include <math.h>
#include <stdio.h>

#include <gkyl_const.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_gyrokinetic_run.h>
#include <gkyl_math.h>
#include <rt_arg_parse.h>

struct d3d_neut_ctx {
  // Miller-like magnetic-geometry parameters.
  double a_shift, Z_axis, R_axis, R0, a_mid, x_inner, B0;
  double kappa, delta, q0, Cy, x_LCFS;

  // Species parameters.
  double me, mi, qe, n0, Te0, Ti0, Tn;

  // Configuration-space domain and recycling coefficient.
  double x_min, x_max, z_min, z_max, Lz, rec_frac;

  // Reference thermal speeds.
  double vte0, vti0;
};

static double
r_x(double x, const struct d3d_neut_ctx *ctx)
{
  return x + ctx->a_mid - ctx->x_inner;
}

static double
qprofile(double R)
{
  const double a[] = {
    154.51071835546747, -921.8584472748003, 1842.1077075366113, -1231.619813170522
  };
  return ((a[0] * R + a[1]) * R + a[2]) * R + a[3];
}

static double
R_rtheta(double r, double theta, void *vctx)
{
  struct d3d_neut_ctx *ctx = vctx;
  return ctx->R_axis - ctx->a_shift * r * r / (2.0 * ctx->R_axis) +
         r * cos(theta + asin(ctx->delta) * sin(theta));
}

static double
Z_rtheta(double r, double theta, void *vctx)
{
  struct d3d_neut_ctx *ctx = vctx;
  return ctx->Z_axis + ctx->kappa * r * sin(theta);
}

static double
dRdr(double r, double theta, void *vctx)
{
  struct d3d_neut_ctx *ctx = vctx;
  return -ctx->a_shift * r / ctx->R_axis + cos(theta + asin(ctx->delta) * sin(theta));
}

static double
dRdtheta(double r, double theta, void *vctx)
{
  struct d3d_neut_ctx *ctx = vctx;
  return -r * sin(theta + asin(ctx->delta) * sin(theta)) * (1.0 + asin(ctx->delta) * cos(theta));
}

static double
dZdr(double r, double theta, void *vctx)
{
  struct d3d_neut_ctx *ctx = vctx;
  return ctx->kappa * sin(theta);
}

static double
dZdtheta(double r, double theta, void *vctx)
{
  struct d3d_neut_ctx *ctx = vctx;
  return ctx->kappa * r * cos(theta);
}

static double
Jr(double r, double theta, void *ctx)
{
  return R_rtheta(r, theta, ctx) * (dRdr(r, theta, ctx) * dZdtheta(r, theta, ctx) -
                                    dRdtheta(r, theta, ctx) * dZdr(r, theta, ctx));
}

struct integrand_ctx {
  struct d3d_neut_ctx *app;
  double r;
};

static double
integrand(double theta, void *vctx)
{
  struct integrand_ctx *ctx = vctx;
  return Jr(ctx->r, theta, ctx->app) / pow(R_rtheta(ctx->r, theta, ctx->app), 2);
}

static double
dPsidr(double r, double theta, void *vctx)
{
  struct d3d_neut_ctx *ctx = vctx;
  struct integrand_ctx int_ctx = {.app = ctx, .r = r};
  struct gkyl_qr_res integral = gkyl_dbl_exp(integrand, &int_ctx, 0.0, 2.0 * M_PI, 7, 1.0e-10);
  const double R_omp = R_rtheta(r, 0.0, ctx);
  return ctx->B0 * ctx->R0 * integral.res / (2.0 * M_PI * qprofile(R_omp));
}

static double
alpha(double r, double theta, double phi, void *vctx)
{
  struct d3d_neut_ctx *ctx = vctx;
  double twrap = theta;
  while (twrap < -M_PI) {
    twrap += 2.0 * M_PI;
  }
  while (twrap > M_PI) {
    twrap -= 2.0 * M_PI;
  }
  struct integrand_ctx int_ctx = {.app = ctx, .r = r};
  struct gkyl_qr_res integral;
  if (twrap > 0.0) {
    integral = gkyl_dbl_exp(integrand, &int_ctx, 0.0, twrap, 7, 1.0e-10);
  } else {
    integral = gkyl_dbl_exp(integrand, &int_ctx, twrap, 0.0, 7, 1.0e-10);
    integral.res = -integral.res;
  }
  return phi - ctx->B0 * ctx->R0 * integral.res / dPsidr(r, theta, ctx);
}

static double
Bphi(double R, void *vctx)
{
  struct d3d_neut_ctx *ctx = vctx;
  return ctx->B0 * ctx->R0 / R;
}

static double
gradr(double r, double theta, void *ctx)
{
  return R_rtheta(r, theta, ctx) / Jr(r, theta, ctx) *
         sqrt(pow(dRdtheta(r, theta, ctx), 2) + pow(dZdtheta(r, theta, ctx), 2));
}

static double
radial_shape(double x, const struct d3d_neut_ctx *ctx)
{
  const double s = fmax(0.0, (x - ctx->x_LCFS) / (ctx->x_max - ctx->x_LCFS));
  return 0.25 + 0.75 * exp(-s * s / 0.36);
}

static void
eval_density(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *vctx)
{
  struct d3d_neut_ctx *ctx = vctx;
  fout[0] = ctx->n0 * radial_shape(xn[0], ctx);
}

static void
eval_elc_temp(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *vctx)
{
  struct d3d_neut_ctx *ctx = vctx;
  fout[0] = ctx->Te0 * (0.25 + 0.75 * radial_shape(xn[0], ctx));
}

static void
eval_ion_temp(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *vctx)
{
  struct d3d_neut_ctx *ctx = vctx;
  fout[0] = ctx->Ti0 * (0.25 + 0.75 * radial_shape(xn[0], ctx));
}

static void
eval_outward_upar(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *vctx)
{
  struct d3d_neut_ctx *ctx = vctx;
  const double cs = sqrt(ctx->Te0 / ctx->mi);
  fout[0] = 0.6 * cs * tanh(4.0 * xn[1] / ctx->Lz);
}

static void
eval_neut_density(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *vctx)
{
  fout[0] = 1.0e15;
}

static void
eval_zero_udrift(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = fout[1] = fout[2] = 0.0;
}

static void
eval_neut_temp(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *vctx)
{
  struct d3d_neut_ctx *ctx = vctx;
  fout[0] = ctx->Tn;
}

// For 2x mapc2p geometry Gkeyll inflates (x,theta) to (x,0,theta)
// before calling these functions, matching the established DIII-D 2x input.
static void
mapc2p(double t, const double *GKYL_RESTRICT xc, double *GKYL_RESTRICT xp, void *vctx)
{
  struct d3d_neut_ctx *ctx = vctx;
  const double x = xc[0], y = xc[1], theta = xc[2];
  const double r = r_x(x, ctx);
  const double R = R_rtheta(r, theta, ctx);
  const double phi = y / ctx->Cy + alpha(r, theta, 0.0, ctx);
  xp[0] = R * cos(phi);
  xp[1] = R * sin(phi);
  xp[2] = Z_rtheta(r, theta, ctx);
}

static void
bfield_func(double t, const double *GKYL_RESTRICT xc, double *GKYL_RESTRICT fout, void *vctx)
{
  struct d3d_neut_ctx *ctx = vctx;
  const double x = xc[0], y = xc[1], theta = xc[2];
  const double r = r_x(x, ctx);
  const double R = R_rtheta(r, theta, ctx);
  const double Bt = Bphi(R, ctx);
  const double Bp = dPsidr(r, theta, ctx) / R * gradr(r, theta, ctx);
  const double drdtheta = dRdtheta(r, theta, ctx);
  const double dzdtheta = dZdtheta(r, theta, ctx);
  const double den = sqrt(drdtheta * drdtheta + dzdtheta * dzdtheta);
  const double Br = Bp * drdtheta / den;
  const double Bz = Bp * dzdtheta / den;
  const double phi = y / ctx->Cy + alpha(r, theta, 0.0, ctx);
  // Set Cartesian components of the magnetic field.
  fout[0] = Br * cos(phi) + Bt * sin(phi);
  fout[1] = Br * sin(phi) - Bt * cos(phi);
  fout[2] = Bz;
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Init(&argc, &argv);
  }
#endif

  const double R_axis = 1.7074685, R_LCFSmid = 2.17;
  const double a_shift = 1.0, x_inner = 0.10, x_outer = 0.05;
  const double Rmid_min = R_LCFSmid - x_inner;
  const double Rmid_max = R_LCFSmid + x_outer;
  const double R0 = 0.5 * (Rmid_min + Rmid_max);
  const double a_mid =
    R_axis / a_shift -
    sqrt(R_axis * (R_axis - 2.0 * a_shift * R_LCFSmid + 2.0 * a_shift * R_axis)) / a_shift;
  const double B_axis = 2.0;
  struct d3d_neut_ctx ctx = {
    .a_shift = a_shift,
    .Z_axis = -0.0014645315,
    .R_axis = R_axis,
    .R0 = R0,
    .a_mid = a_mid,
    .x_inner = x_inner,
    .B0 = B_axis * R_axis / R0,
    .kappa = 1.35,
    .delta = -0.4,
    .x_LCFS = R_LCFSmid - Rmid_min,
    .me = GKYL_ELECTRON_MASS,
    .mi = 2.014 * GKYL_PROTON_MASS,
    .qe = GKYL_ELEMENTARY_CHARGE,
    .n0 = 2.0e19,
    .Te0 = 100.0 * GKYL_ELEMENTARY_CHARGE,
    .Ti0 = 100.0 * GKYL_ELEMENTARY_CHARGE,
    .Tn = 2.0 * GKYL_ELEMENTARY_CHARGE,
    .x_min = R_LCFSmid - Rmid_min,
    .x_max = Rmid_max - Rmid_min,
    .Lz = 2.0 * M_PI - 1.0e-10,
    .rec_frac = 0.95,
  };
  ctx.z_min = -0.5 * ctx.Lz;
  ctx.z_max = 0.5 * ctx.Lz;
  ctx.q0 = qprofile(R0);
  ctx.Cy = (R0 - R_axis) / ctx.q0;
  ctx.vte0 = sqrt(ctx.Te0 / ctx.me);
  ctx.vti0 = sqrt(ctx.Ti0 / ctx.mi);

  const int nx = APP_ARGS_CHOOSE(app_args.xcells[0], 6);
  const int nz = APP_ARGS_CHOOSE(app_args.xcells[1], 8);
  const int nv = APP_ARGS_CHOOSE(app_args.vcells[0], 12);
  const int nmu = APP_ARGS_CHOOSE(app_args.vcells[1], 8);
  double t_end = 5.0e-7;
  int num_frames = 2;
  if (app_args.opt_args[0] != '\0') {
    sscanf(app_args.opt_args, "t_end=%lf,num_frames=%d", &t_end, &num_frames);
  }

  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  // Electron species.
  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = -ctx.qe,
    .mass = ctx.me,
    .vdim = 2,
    .lower = {-6.0 * ctx.vte0, 0.0},
    .upper = {6.0 * ctx.vte0, ctx.me * pow(6.0 * ctx.vte0, 2) / (2.0 * ctx.B0)},
    .cells = {nv, nmu},
    .polarization_density = ctx.n0,
    .projection =
      {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
        .density = eval_density,
        .ctx_density = &ctx,
        .upar = eval_outward_upar,
        .ctx_upar = &ctx,
        .temp = eval_elc_temp,
        .ctx_temp = &ctx,
      },
    .correct =
      {.correct_all_moms = true, .use_last_converged = true, .iter_eps = 1.0e-10, .max_iter = 20},
    .collisionless = {.type = GKYL_GK_COLLISIONLESS_ES},
    .bcs =
      {{.dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX},
       {.dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX},
       {.dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
       {.dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB}},
  };

  // Ion species.
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.qe,
    .mass = ctx.mi,
    .vdim = 2,
    .lower = {-6.0 * ctx.vti0, 0.0},
    .upper = {6.0 * ctx.vti0, ctx.mi * pow(6.0 * ctx.vti0, 2) / (2.0 * ctx.B0)},
    .cells = {nv, nmu},
    .polarization_density = ctx.n0,
    .projection =
      {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
        .density = eval_density,
        .ctx_density = &ctx,
        .upar = eval_outward_upar,
        .ctx_upar = &ctx,
        .temp = eval_ion_temp,
        .ctx_temp = &ctx,
      },
    .correct =
      {.correct_all_moms = true, .use_last_converged = true, .iter_eps = 1.0e-10, .max_iter = 20},
    .collisionless = {.type = GKYL_GK_COLLISIONLESS_ES},
    .bcs =
      {{.dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX},
       {.dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX},
       {.dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
       {.dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB}},
    .boundary_flux_diagnostics = {.num_diag_moments = 1, .diag_moments = {GKYL_F_MOMENT_M0}},
    .num_diag_moments = 1,
    .diag_moments = {GKYL_F_MOMENT_M0},
  };

  // Recycling input shared by the lower and upper parallel targets.
  struct gkyl_gyrokinetic_emission_inp recycling = {
    .num_species = 1,
    .in_species = {"ion"},
    .recycling_frac = ctx.rec_frac,
    .emission_temp = ctx.Tn,
  };

  // Fluid-neutral species.
  struct gkyl_gyrokinetic_neut_species neut = {
    .name = "D0",
    .mass = ctx.mi,
    .vdim = 0,
    .gas_gamma = 5.0 / 3.0,
    .projection =
      {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
        .density = eval_neut_density,
        .ctx_density = &ctx,
        .udrift = eval_zero_udrift,
        .ctx_udrift = &ctx,
        .temp = eval_neut_temp,
        .ctx_temp = &ctx,
      },
    .diffusion =
      {
        // Construct D_n from the production CX rate and advance diffusion with
        // the matrix-free backward-Euler solver.
        .use_reaction_rates = true,
        .write_diagnostics = true,
        .is_implicit = true,
        .implicit_tol = 1.0e-10,
        .implicit_max_iter = 200,
        .min_collision_frequency = 1.0e4,
      },
    .react_neut =
      {
        // CX supplies D_n; ionization depletes the neutral density.
        .write_diagnostics = true,
        .num_react = 2,
        .react_type =
          {{
             .react_id = GKYL_REACT_CX,
             .type_self = GKYL_SELF_PARTNER,
             .ion_id = GKYL_ION_D,
             .elc_nm = "elc",
             .ion_nm = "ion",
             .partner_nm = "D0",
             .ion_mass = ctx.mi,
             .partner_mass = ctx.mi,
           },
           {
             .react_id = GKYL_REACT_IZ,
             .type_self = GKYL_SELF_DONOR,
             .ion_id = GKYL_ION_H,
             .elc_nm = "elc",
             .ion_nm = "ion",
             .donor_nm = "D0",
             .charge_state = 0,
             .ion_mass = ctx.mi,
             .elc_mass = ctx.me,
           }},
      },
    .num_diag_moments = 3,
    .diag_moments = {GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2},
    .bcs =
      {// Neutral flux entering through each target is R times the incident ion flux.
       {
         .dir = 1,
         .edge = GKYL_LOWER_EDGE,
         .type = GKYL_BC_GK_SPECIES_RECYCLE,
         .emission = recycling,
       },
       {
         .dir = 1,
         .edge = GKYL_UPPER_EDGE,
         .type = GKYL_BC_GK_SPECIES_RECYCLE,
         .emission = recycling,
       }
      },
  };

  // Gyrokinetic app input.
  struct gkyl_gk app_inp = {
    .cdim = 2,
    .lower = {ctx.x_min, ctx.z_min},
    .upper = {ctx.x_max, ctx.z_max},
    .cells = {nx, nz},
    .poly_order = 1,
    .basis_type = app_args.basis_type,
    .cfl_frac = 0.4,
    .geometry =
      {
        .geometry_id = GKYL_GEOMETRY_MAPC2P,
        .world = {0.0},
        .mapc2p = mapc2p,
        .c2p_ctx = &ctx,
        .bfield_func = bfield_func,
        .bfield_ctx = &ctx,
        .has_LCFS = false,
        .x_LCFS = ctx.x_LCFS,
      },
    .num_species = 2,
    .species = {elc, ion},
    .num_neut_species = 1,
    .neut_species = {neut},
    .field = {.is_static = true, .zero_init_field = true},
    .parallelism =
      {.use_gpu = app_args.use_gpu, .cuts = {app_args.cuts[0], app_args.cuts[1]}, .comm = comm},
  };
  snprintf(app_inp.name, sizeof(app_inp.name), "%s", app_args.app_name);
  struct gkyl_gyrokinetic_run_inp run_inp = {
    .app_inp = app_inp,
    .time_stepping =
      {
        .t_end = t_end,
        .num_frames = num_frames,
        .write_phase_freq = 1.0,
        .int_diag_calc_num = 20,
        .dt_failure_tol = 1.0e-4,
        .num_failures_max = 20,
        .is_restart = app_args.is_restart,
        .restart_frame = app_args.restart_frame,
        .num_steps = app_args.num_steps,
      },
  };
  gkyl_gyrokinetic_run_simulation(&run_inp);
  gkyl_gyrokinetic_comms_release(comm);
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif
  return 0;
}
