#include <math.h>
#include <stdio.h>

#include <gkyl_const.h>
#include <gkyl_gen_miller_eqdsk.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_gyrokinetic_run.h>
#include <gkyl_math.h>
#include <rt_arg_parse.h>

struct d3d_neut_ctx {
  // Miller geometry parameters.
  double a_shift, Zvert_axis, Rmajor_axis, R0, rminor_LCFS, x_inner, B0;
  double kappa, delta, x_LCFS;
  struct gkyl_gen_miller_eqdsk_inp miller;
  gkyl_gen_miller_eqdsk *equilibrium; // Flux table for profile callbacks; owned by main.

  // Species parameters.
  double me, mi, qe, n0, Te0, Ti0;

  // Neutral temperature: T_n = w*Tn_cold + (1-w)*Tn_ion_frac*T_i, with the
  // cold weight w decaying from the SOL targets over Tn_cold_width (in z).
  double Tn_cold, Tn_ion_frac, Tn_cold_width;

  // Configuration-space domain (x in outboard major radius before the app
  // converts it to psi) and recycling coefficients of the targets and the
  // outer wall.
  double x_min, x_max, z_min, z_max, Lx, Lz, rec_frac_target, rec_frac_wall;

  // Adaptive core source parameters.
  bool adapt_energy_srcCORE, adapt_particle_srcCORE;
  double center_srcCORE[2], sigma_srcCORE[2];
  double energy_srcCORE, particle_srcCORE, floor_srcCORE;

  // Reference thermal speeds.
  double vte0, vti0;
};

static void
miller_kappa(double t, const double *xc, double *GKYL_RESTRICT out, void *vctx)
{
  const struct d3d_neut_ctx *ctx = vctx;
  out[0] = ctx->kappa;
}

static void
miller_delta(double t, const double *xc, double *GKYL_RESTRICT out, void *vctx)
{
  const struct d3d_neut_ctx *ctx = vctx;
  out[0] = ctx->delta;
}

static void
miller_a_shift(double t, const double *xc, double *GKYL_RESTRICT out, void *vctx)
{
  const struct d3d_neut_ctx *ctx = vctx;
  out[0] = ctx->a_shift;
}

static double
qsafety(double R, void *ctx)
{
  // Cubic fit to the DIII-D q profile of the analytic input, used for
  // R >= R_match, which covers the simulated domain. The fit turns negative
  // near the axis, so inside R_match q continues as the exponential with the
  // same value and slope; it only sets psi offsets inside the domain.
  const double a[] = {
    154.51071835546747, -921.8584472748003, 1842.1077075366113, -1231.619813170522
  };
  const double R_match = 2.0;
  const double Rq = fmax(R, R_match);
  const double q = ((a[0] * Rq + a[1]) * Rq + a[2]) * Rq + a[3];
  if (R >= R_match) {
    return q;
  }
  const double dqdR = (3.0 * a[0] * Rq + 2.0 * a[1]) * Rq + a[2];
  return q * exp(dqdR / q * (R - R_match));
}

/**
 * Convert the analytic input's radial coordinate x = r - r_LCFS + x_inner to psi.
 *
 * @param ctx Simulation context with an initialized flux table.
 * @param x Radial coordinate [m].
 * @return Poloidal flux divided by 2*pi [Wb/rad].
 */
static double
psi_at_x(const struct d3d_neut_ctx *ctx, double x)
{
  return gkyl_gen_miller_eqdsk_psi_at_r(ctx->equilibrium, x + ctx->rminor_LCFS - ctx->x_inner);
}

/**
 * Recover the analytic input's radial coordinate x = r - r_LCFS + x_inner from psi.
 *
 * @param ctx Simulation context with an initialized flux table.
 * @param psi Poloidal flux divided by 2*pi [Wb/rad].
 * @return Radial coordinate [m].
 */
static double
x_at_psi(const struct d3d_neut_ctx *ctx, double psi)
{
  return gkyl_gen_miller_eqdsk_r_at_psi(ctx->equilibrium, psi) - ctx->rminor_LCFS + ctx->x_inner;
}

// Miller angle at the normalized poloidal arc length theta_tok.
static double
theta_at_z(const struct d3d_neut_ctx *ctx, double z)
{
  return gkyl_gen_miller_eqdsk_theta_at_theta_tok(ctx->equilibrium, z);
}

struct source_profile_ctx {
  const struct d3d_neut_ctx *app;
  double center[2], sigma[2], floor;
};

// Core source Gaussian in the analytic input's radial coordinate and Miller angle.
static void
source_profile(double t, const double *xn, double *fout, void *vctx)
{
  const struct source_profile_ctx *src = vctx;
  const double dx = (x_at_psi(src->app, xn[0]) - src->center[0]) / src->sigma[0];
  const double dtheta = (theta_at_z(src->app, xn[1]) - src->center[1]) / src->sigma[1];
  fout[0] = exp(-0.5 * (dx * dx + dtheta * dtheta)) + src->floor;
}

// Radial tanh profile, a*tanh(b*(x - c)) + d, fitted to experimental data.
static double
tanh_profile(double x, double a, double b, double c, double d)
{
  return a * tanh(b * (x - c)) + d;
}

static void
eval_density(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *vctx)
{
  const double x = x_at_psi(vctx, xn[0]);
  fout[0] = tanh_profile(x, -8.8e18, 34.0, 0.093, 1.1e19);
}

static void
eval_elc_temp(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *vctx)
{
  const double eV = GKYL_ELEMENTARY_CHARGE;
  const double x = x_at_psi(vctx, xn[0]);
  fout[0] = tanh_profile(x, -240.0 * eV, 32.5, 0.08, 260.0 * eV);
}

static void
eval_ion_temp(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *vctx)
{
  const double eV = GKYL_ELEMENTARY_CHARGE;
  const double x = x_at_psi(vctx, xn[0]);
  fout[0] = tanh_profile(x, -1.6e3 * eV, 3.2, 0.16, 14.0 * eV);
}

static void
eval_zero_upar(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *vctx)
{
  fout[0] = 0.0;
}

// Ion parallel flow: zero inside the LCFS and, in the SOL, rising to
// -c_s at the lower target and +c_s at the upper target within about 0.5 rad of
// each target. c_s = sqrt(Te0/m_i) uses the reference electron temperature.
static void
eval_ion_upar(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *vctx)
{
  struct d3d_neut_ctx *ctx = vctx;
  const double z = xn[1];
  const double cs = sqrt(ctx->Te0 / ctx->mi);
  const double target_factor =
    1.0 / cosh((z - 0.5 * ctx->Lz) / 0.5) - 1.0 / cosh((z + 0.5 * ctx->Lz) / 0.5);
  fout[0] = xn[0] < ctx->x_LCFS ? 0.0 : cs * target_factor;
}

// Weight of the cold (recycled) neutral temperature: 1 at the SOL
// targets, decaying over Tn_cold_width along z, and 0 in the core.
static void
eval_neut_cold_weight(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *vctx)
{
  struct d3d_neut_ctx *ctx = vctx;
  const double x = xn[0], z = xn[1];
  if (x < ctx->x_LCFS) {
    fout[0] = 0.0;
    return;
  }
  const double w = exp(-(z - ctx->z_min) / ctx->Tn_cold_width) +
                   exp(-(ctx->z_max - z) / ctx->Tn_cold_width);
  fout[0] = fmin(w, 1.0);
}

// Neutral density peaked symmetrically at both targets and suppressed inside the LCFS.
static void
eval_neut_density(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *vctx)
{
  struct d3d_neut_ctx *ctx = vctx;
  const double x = x_at_psi(ctx, xn[0]), theta = theta_at_z(ctx, xn[1]);
  const double x_LCFS = ctx->x_inner; // LCFS position in the analytic radial coordinate.
  const double alpha_neut = 0.1 * ctx->n0;
  const double radial_factor = 1.0 - tanh((x_LCFS - x) / 0.05);
  // Symmetric in theta: the same peak at both targets.
  const double dist_to_target = 0.5 * ctx->Lz - fabs(theta);
  fout[0] = alpha_neut * radial_factor * (pow(1.0 / cosh(dist_to_target / 0.5), 2) + 1.0e-6);
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
  const double x = x_at_psi(ctx, xn[0]);
  fout[0] = 0.6 * ctx->Ti0 * (0.5 * (1.0 + tanh(2.0 - 25.0 * x)) + 0.1);
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

  const double Rmajor_axis = 1.7074685, R_LCFSmid = 2.17;
  const double a_shift = 1.0, x_inner = 0.10, x_outer = 0.05;
  const double R0 = R_LCFSmid + 0.5 * (x_outer - x_inner); // Outboard reference radius [m].
  // B_phi clockwise viewed from above, so the ion grad-B drift is downward.
  // The Miller generator fixes the plasma current along -phi, so here the
  // current is parallel to B_phi.
  const double bmag_axis = -2.0;

  // The Miller shape callbacks need the final context address.
  struct d3d_neut_ctx ctx = {
    .a_shift = a_shift,
    .Zvert_axis = -0.0014645315,
    .Rmajor_axis = Rmajor_axis,
    .R0 = R0,
    .x_inner = x_inner,
    .kappa = 1.35,
    .delta = -0.4,
    .me = GKYL_ELECTRON_MASS,
    .mi = 2.01410177811 * GKYL_PROTON_MASS,
    .qe = GKYL_ELEMENTARY_CHARGE,
    .n0 = 2.0e19,
    .Te0 = 300.0 * GKYL_ELEMENTARY_CHARGE,
    .Ti0 = 300.0 * GKYL_ELEMENTARY_CHARGE,
    .Lx = x_inner + x_outer,
    .Lz = 2.0 * M_PI - 1.0e-10,
    .rec_frac_target = 0.98,
    .rec_frac_wall = 0.9,
    .Tn_cold = 3.0 * GKYL_ELEMENTARY_CHARGE, // Franck-Condon recycled atoms.
    .Tn_ion_frac = 0.5,
    .Tn_cold_width = 0.1, // About a 3 eV atom's CX mean free path [rad].
  };
  ctx.miller = (struct gkyl_gen_miller_eqdsk_inp){
    .input_coord_limits_type = {GKYL_GK_GEOMETRY_COORD_RADIAL_R_OMP, 0, 0},
    .Rmajor_axis = Rmajor_axis,
    .Zvert_axis = ctx.Zvert_axis,
    .bmag_axis = bmag_axis,
    .Rmajor_LCFS = R_LCFSmid,
    .kappa = miller_kappa,
    .kappa_ctx = &ctx,
    .delta = miller_delta,
    .delta_ctx = &ctx,
    .a_shift = miller_a_shift,
    .a_shift_ctx = &ctx,
    .qsafety = qsafety,
    .has_Rmajor_omp_bounds = true,
  };
  const double rminor_LCFS_guess = gkyl_gen_miller_eqdsk_rminor_at_R_omp(&ctx.miller, R_LCFSmid);
  ctx.x_max = gkyl_gen_miller_eqdsk_R_omp_at_rminor(&ctx.miller, rminor_LCFS_guess + x_outer);
  ctx.miller.Rmajor_omp_bounds[0] = Rmajor_axis;
  ctx.miller.Rmajor_omp_bounds[1] = ctx.x_max;
  ctx.equilibrium = gkyl_gen_miller_eqdsk_inew(&ctx.miller);
  int failed = ctx.equilibrium == NULL;
  if (!failed) {
    const struct gkyl_gen_miller_eqdsk_info *info = gkyl_gen_miller_eqdsk_get_info(ctx.equilibrium);
    ctx.rminor_LCFS = info->rminor_LCFS;
    // Put the LCFS at 2/3 of the psi grid, keeping the analytic outer surface.
    ctx.x_LCFS = gkyl_round_to_e16_9(info->psi_boundary);
    const double psi_upper = psi_at_x(&ctx, ctx.Lx);
    const double psi_lower = (ctx.x_LCFS - (2.0 / 3.0) * psi_upper) / (1.0 / 3.0);
    ctx.x_min = gkyl_gen_miller_eqdsk_R_omp_at_psi(ctx.equilibrium, psi_lower);
    ctx.B0 = gkyl_gen_miller_eqdsk_bmag_at_RZ(ctx.equilibrium, R0, ctx.Zvert_axis);
    failed = !(psi_lower > info->psi_axis) || !isfinite(ctx.x_min) || !isfinite(ctx.B0) ||
             ctx.B0 <= 0.0;
  }
  ctx.z_min = -0.5 * ctx.Lz;
  ctx.z_max = 0.5 * ctx.Lz;
  ctx.vte0 = sqrt(ctx.Te0 / ctx.me);
  ctx.vti0 = sqrt(ctx.Ti0 / ctx.mi);

  // Adaptive core source: injects 0.75 MW per species and replaces the
  // particles lost through the inner radial boundary. The binormal coordinate
  // is the field-line label alpha, whose full torus spans 2 pi, so a 2x volume
  // integral is per radian of alpha and the equivalent power is 0.75 MW/(2 pi).
  ctx.adapt_energy_srcCORE = true;
  ctx.adapt_particle_srcCORE = true;
  ctx.energy_srcCORE = 0.75e6 / (2.0 * M_PI);
  ctx.particle_srcCORE = 0.0;
  ctx.center_srcCORE[1] = 0.25 * ctx.Lz; // Ion source; the electron source is at -Lz/4.
  ctx.sigma_srcCORE[0] = 0.03 * ctx.Lx;
  ctx.sigma_srcCORE[1] = ctx.Lz / 6.0;
  ctx.floor_srcCORE = 1.0e-10;

  // The LCFS sits at 2/3 of the radial domain, so nx must be a multiple of 3.
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

  failed = failed || nx < 3 || nx % 3 != 0;
  int any_failed;
  gkyl_comm_allreduce_host(comm, GKYL_INT, GKYL_MAX, 1, &failed, &any_failed);
  if (any_failed) {
    fprintf(stderr, "DIII-D Miller input requires a valid equilibrium and nx divisible by 3.\n");
    gkyl_gen_miller_eqdsk_release(ctx.equilibrium);
    gkyl_gyrokinetic_comms_release(comm);
#ifdef GKYL_HAVE_MPI
    if (app_args.use_mpi) {
      MPI_Finalize();
    }
#endif
    return 1;
  }

  // The core source is centered on the inner radial boundary.
  ctx.center_srcCORE[0] =
    x_at_psi(&ctx, gkyl_gen_miller_eqdsk_psi_at_R_omp(ctx.equilibrium, ctx.x_min));
  struct source_profile_ctx core_elc = {
    .app = &ctx,
    .center = {ctx.center_srcCORE[0], -ctx.center_srcCORE[1]},
    .sigma = {ctx.sigma_srcCORE[0], ctx.sigma_srcCORE[1]},
    .floor = ctx.floor_srcCORE,
  };
  struct source_profile_ctx core_ion = core_elc;
  core_ion.center[1] = ctx.center_srcCORE[1];

  // Core sources, adapting to the losses through the inner radial boundary.
  struct gkyl_gyrokinetic_projection proj_srcCORE_e = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_GAUSSIAN,
    .configuration_profile = source_profile,
    .configuration_profile_ctx = &core_elc,
    .total_num_particles = ctx.particle_srcCORE,
    .total_kin_energy = ctx.energy_srcCORE,
    .temp_min = 0.5 * ctx.Te0,
    .temp_max = 5.0 * ctx.Te0,
    .f_floor = ctx.floor_srcCORE,
  };
  struct gkyl_gyrokinetic_projection proj_srcCORE_i = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_GAUSSIAN,
    .configuration_profile = source_profile,
    .configuration_profile_ctx = &core_ion,
    .total_num_particles = ctx.particle_srcCORE,
    .total_kin_energy = ctx.energy_srcCORE,
    .temp_max = 5.0 * ctx.Te0,
    .f_floor = ctx.floor_srcCORE,
  };
  struct gkyl_gyrokinetic_adapt_source adapt_srcCORE_e = {
    .adapt_to_species = "elc",
    .adapt_particle = ctx.adapt_particle_srcCORE,
    .adapt_energy = ctx.adapt_energy_srcCORE,
    .num_boundaries = 1, // Only the inner radial boundary.
    .dir = {0},
    .edge = {GKYL_LOWER_EDGE},
  };
  struct gkyl_gyrokinetic_adapt_source adapt_srcCORE_i = {
    .adapt_to_species = "ion",
    .adapt_particle = ctx.adapt_particle_srcCORE,
    .adapt_energy = ctx.adapt_energy_srcCORE,
    .num_boundaries = 1, // Only the inner radial boundary.
    .dir = {0},
    .edge = {GKYL_LOWER_EDGE},
  };

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
        .upar = eval_zero_upar,
        .ctx_upar = &ctx,
        .temp = eval_elc_temp,
        .ctx_temp = &ctx,
      },
    .collisionless = {.type = GKYL_GK_COLLISIONLESS_ES},
    .react_neut =
      {
        .num_react = 2,
        .write_diagnostics = true,
        .react_type =
          {{
             .react_id = GKYL_REACT_IZ,
             .type_self = GKYL_SELF_ELC,
             .ion_id = GKYL_ION_H,
             .elc_nm = "elc",
             .ion_nm = "ion",
             .donor_nm = "D0",
             .charge_state = 0,
             .ion_mass = ctx.mi,
             .elc_mass = ctx.me,
           },
           {
             .react_id = GKYL_REACT_RECOMB,
             .type_self = GKYL_SELF_ELC,
             .ion_id = GKYL_ION_H,
             .elc_nm = "elc",
             .ion_nm = "ion",
             .recvr_nm = "D0",
             .charge_state = 0,
             .ion_mass = ctx.mi,
             .elc_mass = ctx.me,
           }},
      },
    .source =
      {
        .source_id = GKYL_PROJ_SOURCE,
        .num_sources = 1,
        .num_adapt_sources = 1,
        .projection[0] = proj_srcCORE_e,
        .adapt[0] = adapt_srcCORE_e,
        .diagnostics =
          {
            .num_diag_moments = 1,
            .diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
            .num_integrated_diag_moments = 1,
            .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
          },
      },
    .bcs =
      {{.dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
       {.dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
       {.dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH_CONDUCTING},
       {.dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH_CONDUCTING}},
    .num_diag_moments = 9,
    .diag_moments =
      {GKYL_F_MOMENT_HAMILTONIAN, GKYL_F_MOMENT_BIMAXWELLIAN, GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1,
       GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M3PAR,
       GKYL_F_MOMENT_M3PERP},
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
    .boundary_flux_diagnostics =
      {.num_integrated_diag_moments = 1, .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN}},
    .time_rate_diagnostics = true,
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
        .upar = eval_ion_upar,
        .ctx_upar = &ctx,
        .temp = eval_ion_temp,
        .ctx_temp = &ctx,
      },
    .collisionless = {.type = GKYL_GK_COLLISIONLESS_ES},
    .react_neut =
      {
        .num_react = 3,
        .write_diagnostics = true,
        .react_type =
          {{
             .react_id = GKYL_REACT_IZ,
             .type_self = GKYL_SELF_ION,
             .ion_id = GKYL_ION_H,
             .elc_nm = "elc",
             .ion_nm = "ion",
             .donor_nm = "D0",
             .charge_state = 0,
             .ion_mass = ctx.mi,
             .elc_mass = ctx.me,
           },
           {
             .react_id = GKYL_REACT_RECOMB,
             .type_self = GKYL_SELF_ION,
             .ion_id = GKYL_ION_H,
             .elc_nm = "elc",
             .ion_nm = "ion",
             .recvr_nm = "D0",
             .charge_state = 0,
             .ion_mass = ctx.mi,
             .elc_mass = ctx.me,
           },
           {
             .react_id = GKYL_REACT_CX,
             .type_self = GKYL_SELF_ION,
             .ion_id = GKYL_ION_H,
             .ion_nm = "ion",
             .partner_nm = "D0",
             .ion_mass = ctx.mi,
             .partner_mass = ctx.mi,
           }},
      },
    .source =
      {
        .source_id = GKYL_PROJ_SOURCE,
        .num_sources = 1,
        .num_adapt_sources = 1,
        .projection[0] = proj_srcCORE_i,
        .adapt[0] = adapt_srcCORE_i,
        .diagnostics =
          {
            .num_diag_moments = 1,
            .diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
            .num_integrated_diag_moments = 1,
            .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
          },
      },
    .bcs =
      {{.dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
       {.dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
       {.dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH_CONDUCTING},
       {.dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH_CONDUCTING}},
    .num_diag_moments = 9,
    .diag_moments =
      {GKYL_F_MOMENT_HAMILTONIAN, GKYL_F_MOMENT_BIMAXWELLIAN, GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1,
       GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M3PAR,
       GKYL_F_MOMENT_M3PERP},
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
    // M0 boundary fluxes are kept to check the recycled particle balance.
    .boundary_flux_diagnostics =
      {
        .num_diag_moments = 1,
        .diag_moments = {GKYL_F_MOMENT_M0},
        .num_integrated_diag_moments = 2,
        .integrated_diag_moments = {GKYL_F_MOMENT_M0, GKYL_F_MOMENT_HAMILTONIAN},
      },
    .time_rate_diagnostics = true,
  };

  // Recycling of the ion flux at the parallel targets and at the outer wall.
  struct gkyl_gyrokinetic_emission_inp recycling_target = {
    .num_species = 1,
    .in_species = {"ion"},
    .recycling_frac = ctx.rec_frac_target,
  };
  struct gkyl_gyrokinetic_emission_inp recycling_wall = {
    .num_species = 1,
    .in_species = {"ion"},
    .recycling_frac = ctx.rec_frac_wall,
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
        // Construct D_n = (T_n/m_n)/(n_e <sigma v>_cx) from the neutral
        // temperature and the production CX rate, and advance diffusion with
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
        // CX supplies D_n; ionization depletes and recombination replenishes
        // the neutral density.
        .write_diagnostics = true,
        .num_react = 3,
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
           },
           {
             .react_id = GKYL_REACT_RECOMB,
             .type_self = GKYL_SELF_RECVR,
             .ion_id = GKYL_ION_H,
             .elc_nm = "elc",
             .ion_nm = "ion",
             .recvr_nm = "D0",
             .charge_state = 0,
             .ion_mass = ctx.mi,
             .elc_mass = ctx.me,
           }},
      },
    // Recycled atoms are cold near the targets and heat toward a fraction of
    // T_i by charge exchange away from them.
    .temp_model =
      {
        .temp_cold = ctx.Tn_cold,
        .ion_temp_frac = ctx.Tn_ion_frac,
        .ion_nm = "ion",
        .cold_weight = eval_neut_cold_weight,
        .cold_weight_ctx = &ctx,
      },
    .num_diag_moments = 3,
    .diag_moments = {GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2},
    .bcs =
      {
        // Neutrals leaving through the inner (core) boundary are absorbed:
        // they would ionize deeper in the core.
        {.dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
        // The neutral flux entering through the outer wall and through each
        // SOL target is R times the incident ion flux there.
        {
          .dir = 0,
          .edge = GKYL_UPPER_EDGE,
          .type = GKYL_BC_GK_SPECIES_RECYCLE,
          .emission = recycling_wall,
        },
        {
          .dir = 1,
          .edge = GKYL_LOWER_EDGE,
          .type = GKYL_BC_GK_SPECIES_RECYCLE,
          .emission = recycling_target,
        },
        {
          .dir = 1,
          .edge = GKYL_UPPER_EDGE,
          .type = GKYL_BC_GK_SPECIES_RECYCLE,
          .emission = recycling_target,
        },
      },
  };

  // Zero potential where the LCFS meets the targets, as in the TCV IWL tests.
  struct gkyl_poisson_bias_line target_corner_bcs[] = {
    {.perp_dirs = {0, 1}, .perp_coords = {ctx.x_LCFS, ctx.z_min}, .val = 0.0},
    {.perp_dirs = {0, 1}, .perp_coords = {ctx.x_LCFS, ctx.z_max}, .val = 0.0},
  };
  struct gkyl_poisson_bias_line_list bias_line_list = {.num_bias_line = 2, .bl = target_corner_bcs};

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
        // The app writes <app-name>-miller.geqdsk, builds tokamak geometry from it
        // and derives the LCFS split from that equilibrium.
        .geometry_id = GKYL_GEOMETRY_MILLER,
        .miller_grid_info = ctx.miller,
      },
    .num_species = 2,
    .species = {elc, ion},
    .num_neut_species = 1,
    .neut_species = {neut},
    .field =
      {
        .gkfield_id = GKYL_GK_FIELD_ES,
        .polarization_bmag = ctx.B0,
        .poisson_bcs =
          {
            {.dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}},
            {.dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}},
          },
        .bias_line_list = &bias_line_list,
        .time_rate_diagnostics = true,
      },
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
  gkyl_gen_miller_eqdsk_release(ctx.equilibrium);
  gkyl_gyrokinetic_comms_release(comm);
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif
  return 0;
}
