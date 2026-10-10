// Pulsar polar-cap discharge with an adaptive pair source (special-relativistic Vlasov-Maxwell,
// 1x1v), a scaled 1D version of the setup of Chernoglazov, Philippov and Timokhin (2024, ApJL
// 974, L32) without photons: a positron atmosphere 100 exp(-x) held by gravity (acceleration -T
// where the atmosphere density exceeds 0.5) over a magnetospheric pair plasma, electrons with the
// Goldreich-Julian charge density -(1 + 0.8 x/L) and a flow that carries the applied
// super-Goldreich-Julian current J0 = 2 on top of the positron flux, so that E = 0 satisfies
// Gauss's law for the perturbation field and the plasma carries J0 at t = 0. Fixed-function
// boundaries keep the star-side and magnetospheric states in the ghost cells. As the atmosphere
// settles the electron flux can no longer carry J0 and a gap opens above the atmosphere,
// dE/dt = -(j + J0); once particles pass p > 25 (a proxy for curvature-photon emission) each
// species injects pairs at 0.5 (n_e + n_p above the threshold) per unit time, those bred by
// particles moving away from the star drifting at p = 2 and those bred by particles moving toward
// the star at p = -2, with T = 0.2 mc^2 and the rescaled densities smoothed by ten passes of the
// Gaussian filter; the cascade screens the gap into a sustained discharge. The momentum grid is
// linear (cell 0.3) at the origin and exponential out to p = 200.
// Figures of merit (serendipity p2, 256 x 192 cells): the gap field above the atmosphere reaches
// -0.73, -1.05 and -1.38 at t = 10, 15 and 20 (x = 6 to 9; field energy 1.51, 5.89 and 17.0), the
// cascade starts at t = 29.8, 368 pairs per unit area are created by t = 60, when the discharge
// injects 23.9 pairs per unit time with the electron density above the threshold at 1.83 and the
// gap field at -2.28 (x = 10); tensor p1 and tensor p2 give the same picture to 4% and 0.5%, and
// the C driver and its Lua twin agree to round-off. The run takes about five minutes.
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>
#include <gkyl_util.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#ifdef GKYL_HAVE_NCCL
#include <gkyl_nccl_comm.h>
#endif
#endif

#include <rt_arg_parse.h>

struct pulsar_sr_ctx {
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using normalized code units).
  double epsilon0; // Permittivity of free space.
  double mu0; // Permeability of free space.
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.
  double mass_pos; // Positron mass.
  double charge_pos; // Positron charge.

  double T; // Temperature of the plasma and of the injected pairs (units of mc^2).
  double J0; // Applied current (x-direction).
  double grav; // Strength of gravity (inverse scale height of the atmosphere).
  double n_atm; // Atmosphere density at the star.
  double n_mag; // Positron density of the magnetospheric plasma.
  double v_mag; // Flow velocity of the magnetospheric positrons (units of c).
  double rate_src; // Pair injection rate per unit density above the threshold (each species).
  double p_thresh; // Momentum above which particles breed pairs.
  double f_thresh; // Distribution function below which the threshold moment is not accumulated.
  double u_src; // Drift (four-) velocity of the injected pairs (x-direction).
  int num_filters; // Number of Gaussian filter passes on the rescaled source density.

  // Simulation parameters.
  int Nx; // Cell count (configuration space: x-direction).
  int Npx; // Cell count (momentum space: px-direction).
  double Lx; // Domain size (configuration space: x-direction).
  double px_max; // Domain boundary (momentum space: px-direction).
  double px_lin; // Momentum map: cell size of the linear part of the map at the origin.
  int poly_order; // Polynomial order.
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  int field_energy_calcs; // Number of times to calculate field energy.
  int integrated_mom_calcs; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct pulsar_sr_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  double epsilon0 = 1.0; // Permittivity of free space.
  double mu0 = 1.0; // Permeability of free space.
  double mass_elc = 1.0; // Electron mass.
  double charge_elc = -1.0; // Electron charge.
  double mass_pos = 1.0; // Positron mass.
  double charge_pos = 1.0; // Positron charge.

  double T = 0.2; // Temperature of the plasma and of the injected pairs (units of mc^2).
  double J0 = 2.0; // Applied current (x-direction).
  double grav = 1.0; // Strength of gravity (inverse scale height of the atmosphere).
  double n_atm = 100.0; // Atmosphere density at the star.
  double n_mag = 5.0; // Positron density of the magnetospheric plasma.
  double v_mag = 0.8; // Flow velocity of the magnetospheric positrons (units of c).
  double rate_src = 0.5; // Pair injection rate per unit density above the threshold (each species).
  double p_thresh = 25.0; // Momentum above which particles breed pairs.
  double f_thresh =
    1.0e-3; // Distribution function below which the threshold moment is not accumulated.
  double u_src = 2.0; // Drift (four-) velocity of the injected pairs (x-direction).
  int num_filters = 10; // Number of Gaussian filter passes on the rescaled source density.

  // Simulation parameters.
  int Nx = 256; // Cell count (configuration space: x-direction).
  int Npx = 192; // Cell count (momentum space: px-direction).
  double Lx = 50.0; // Domain size (configuration space: x-direction).
  double px_max = 200.0; // Domain boundary (momentum space: px-direction).
  double px_lin = 0.14; // Momentum map: cell size of the linear part of the map at the origin.
  int poly_order = 2; // Polynomial order.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 60.0; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs =
    INT_MAX; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct pulsar_sr_ctx ctx = {
    .pi = pi,
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .mass_pos = mass_pos,
    .charge_pos = charge_pos,
    .T = T,
    .J0 = J0,
    .grav = grav,
    .n_atm = n_atm,
    .n_mag = n_mag,
    .v_mag = v_mag,
    .rate_src = rate_src,
    .p_thresh = p_thresh,
    .f_thresh = f_thresh,
    .u_src = u_src,
    .num_filters = num_filters,
    .Nx = Nx,
    .Npx = Npx,
    .Lx = Lx,
    .px_max = px_max,
    .px_lin = px_lin,
    .poly_order = poly_order,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .num_frames = num_frames,
    .field_energy_calcs = field_energy_calcs,
    .integrated_mom_calcs = integrated_mom_calcs,
    .integrated_L2_f_calcs = integrated_L2_f_calcs,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };

  return ctx;
}

// Lab-frame positron density: the atmosphere over the magnetospheric plasma.
// The volatile intermediates keep clang from fusing the products and sums into FMAs, so the
// C driver and its Lua twin evaluate the same floating-point operations.
static double
pos_density(const struct pulsar_sr_ctx *app, double x)
{
  volatile double n_atm = app->n_atm * exp(-app->grav * x);
  return n_atm + app->n_mag;
}

// Lab-frame electron density: the positron density plus the Goldreich-Julian charge density
// -(1 + 0.8 x/L), so that Gauss's law holds for the perturbation field E = 0 at t = 0.
static double
elc_density(const struct pulsar_sr_ctx *app, double x)
{
  volatile double gj = (0.8 * x) / app->Lx;
  volatile double n_gj = 1.0 + gj;
  return n_gj + pos_density(app, x);
}

void
evalElcDensityInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct pulsar_sr_ctx *app = ctx;
  double x = xn[0];

  // The electron flux carries the applied current on top of the positron flux.
  double n = elc_density(app, x);
  volatile double flux_pos = app->n_mag * app->v_mag;
  volatile double beta = (app->J0 + flux_pos) / n;
  volatile double beta_sq = beta * beta;

  // Rest-frame electron density for the LTE projection.
  fout[0] = n * sqrt(1.0 - beta_sq);
}

void
evalElcVDriftInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct pulsar_sr_ctx *app = ctx;
  double x = xn[0];

  // The electron flux carries the applied current on top of the positron flux.
  volatile double flux_pos = app->n_mag * app->v_mag;
  volatile double beta = (app->J0 + flux_pos) / elc_density(app, x);
  volatile double beta_sq = beta * beta;

  // Electron drift (four-) velocity.
  fout[0] = beta / sqrt(1.0 - beta_sq);
}

void
evalPosDensityInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct pulsar_sr_ctx *app = ctx;
  double x = xn[0];

  // The magnetospheric positron flux is shared with the atmosphere.
  double n = pos_density(app, x);
  volatile double beta = (app->n_mag * app->v_mag) / n;
  volatile double beta_sq = beta * beta;

  // Rest-frame positron density for the LTE projection.
  fout[0] = n * sqrt(1.0 - beta_sq);
}

void
evalPosVDriftInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct pulsar_sr_ctx *app = ctx;
  double x = xn[0];

  // The magnetospheric positron flux is shared with the atmosphere.
  volatile double beta = (app->n_mag * app->v_mag) / pos_density(app, x);
  volatile double beta_sq = beta * beta;

  // Positron drift (four-) velocity.
  fout[0] = beta / sqrt(1.0 - beta_sq);
}

void
evalTempInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct pulsar_sr_ctx *app = ctx;

  // Temperature.
  fout[0] = app->T;
}

void
evalAppliedAccel(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct pulsar_sr_ctx *app = ctx;
  double x = xn[0];

  double n_atm = app->n_atm * exp(-app->grav * x);

  // Gravity acts where the atmosphere is denser than 0.5.
  double accel = 0.0;
  if (n_atm > 0.5) {
    accel = -app->grav * app->T;
  }

  fout[0] = accel, fout[1] = 0.0, fout[2] = 0.0;
}

void
evalDensitySource(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct pulsar_sr_ctx *app = ctx;

  // Pair injection rate per unit density above the threshold. The LTE projection takes the
  // rest-frame density, so divide by the Lorentz factor of the drift to inject rate_src pairs
  // per unit time in the lab frame.
  volatile double u_sq = app->u_src * app->u_src;
  fout[0] = app->rate_src / sqrt(1.0 + u_sq);
}

void
evalVDriftSourceUp(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct pulsar_sr_ctx *app = ctx;

  // Drift (four-) velocity of the pairs bred by particles moving away from the star.
  fout[0] = app->u_src;
}

void
evalVDriftSourceLo(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct pulsar_sr_ctx *app = ctx;

  // Drift (four-) velocity of the pairs bred by particles moving toward the star.
  fout[0] = -app->u_src;
}

void
evalFieldInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  // Set electric field.
  fout[0] = 0.0, fout[1] = 0.0, fout[2] = 0.0;
  // Set magnetic field.
  fout[3] = 0.0, fout[4] = 0.0, fout[5] = 0.0;
  // Set auxiliary fields.
  fout[6] = 0.0, fout[7] = 0.0;
}

void
evalAppliedCurrent(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct pulsar_sr_ctx *app = ctx;

  // Applied current (x-direction).
  fout[0] = app->J0, fout[1] = 0.0, fout[2] = 0.0;
}

void
mapc2p_px(double t, const double *GKYL_RESTRICT vc, double *GKYL_RESTRICT vp, void *ctx)
{
  struct pulsar_sr_ctx *app = ctx;
  double px_c = vc[0];

  double px_max = app->px_max;
  double px_lin = app->px_lin;
  int Npx = app->Npx;

  // Linear-to-exponential momentum map: cells of size px_lin at the origin, growing
  // exponentially to the domain boundary px_max.
  volatile double lin = px_lin * px_c * Npx;
  if (px_c < 0.0) {
    volatile double ex = exp(-px_c * log(px_max));
    vp[0] = lin - ex + 1.0;
  } else {
    volatile double ex = exp(px_c * log(px_max));
    vp[0] = lin + ex - 1.0;
  }
}

// Snap a trigger whose next time exceeds t_end only by round-off back to t_end.
static void
snap_trigger_to_t_end(struct gkyl_tm_trigger *trig, double t_end)
{
  if (trig->tcurr > t_end && trig->tcurr <= t_end * (1.0 + 1.0e-10)) {
    trig->tcurr = t_end;
  }
}

void
write_data(struct gkyl_tm_trigger *iot, gkyl_vlasov_app *app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_write) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_vlasov_app_write(app, t_curr, frame);
    gkyl_vlasov_app_write_field_energy(app);
    gkyl_vlasov_app_write_integrated_mom(app);
    gkyl_vlasov_app_write_integrated_L2_f(app);
    gkyl_vlasov_app_write_mom(app, t_curr, frame);
  }
}

void
calc_field_energy(struct gkyl_tm_trigger *fet, gkyl_vlasov_app *app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(fet, t_curr) || force_calc) {
    gkyl_vlasov_app_calc_field_energy(app, t_curr);
  }
}

void
calc_integrated_mom(
  struct gkyl_tm_trigger *imt, gkyl_vlasov_app *app, double t_curr, bool force_calc
)
{
  if (gkyl_tm_trigger_check_and_bump(imt, t_curr) || force_calc) {
    gkyl_vlasov_app_calc_integrated_mom(app, t_curr);
  }
}

void
calc_integrated_L2_f(
  struct gkyl_tm_trigger *l2t, gkyl_vlasov_app *app, double t_curr, bool force_calc
)
{
  if (gkyl_tm_trigger_check_and_bump(l2t, t_curr) || force_calc) {
    gkyl_vlasov_app_calc_integrated_L2_f(app, t_curr);
  }
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

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct pulsar_sr_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int NPX = APP_ARGS_CHOOSE(app_args.vcells[0], ctx.Npx);
  ctx.Npx = NPX; // The momentum map uses the cell count.

  int nrank = 1; // Number of processors in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif

  int ccells[] = {NX};
  int cdim = sizeof(ccells) / sizeof(ccells[0]);

  int cuts[cdim];
#ifdef GKYL_HAVE_MPI
  for (int d = 0; d < cdim; d++) {
    if (app_args.use_mpi) {
      cuts[d] = app_args.cuts[d];
    } else {
      cuts[d] = 1;
    }
  }
#else
  for (int d = 0; d < cdim; d++) {
    cuts[d] = 1;
  }
#endif

  // Construct communicator for use in app.
  struct gkyl_comm *comm;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_gpu && app_args.use_mpi) {
#ifdef GKYL_HAVE_NCCL
    comm = gkyl_nccl_comm_new(&(struct gkyl_nccl_comm_inp){.mpi_comm = MPI_COMM_WORLD});
#else
    printf(" Using -g and -M together requires NCCL.\n");
    assert(0 == 1);
#endif
  } else if (app_args.use_mpi) {
    comm = gkyl_mpi_comm_new(&(struct gkyl_mpi_comm_inp){.mpi_comm = MPI_COMM_WORLD});
  } else {
    comm = gkyl_null_comm_inew(&(struct gkyl_null_comm_inp){.use_gpu = app_args.use_gpu});
  }
#else
  comm = gkyl_null_comm_inew(&(struct gkyl_null_comm_inp){.use_gpu = app_args.use_gpu});
#endif

  int my_rank;
  gkyl_comm_get_rank(comm, &my_rank);
  int comm_size;
  gkyl_comm_get_size(comm, &comm_size);

  int ncuts = 1;
  for (int d = 0; d < cdim; d++) {
    ncuts *= cuts[d];
  }

  if (ncuts != comm_size) {
    if (my_rank == 0) {
      fprintf(stderr, "*** Number of ranks, %d, does not match total cuts, %d!\n", comm_size, ncuts);
    }
    goto mpifinalize;
  }

  // Electrons.
  struct gkyl_vlasov_kinetic_species elc = {
    .model_id = GKYL_MODEL_SR,
    .bcx = {.lower = {.type = GKYL_SPECIES_FIXED_FUNC}, .upper = {.type = GKYL_SPECIES_FIXED_FUNC}},
    .app_accel = evalAppliedAccel,
    .app_accel_ctx = &ctx,
    .lower = {-1.0},
    .upper = {1.0},
    .cells = {NPX},
    .mapc2p_vel = {{.mapc2p_vel_func = mapc2p_px, .mapc2p_vel_ctx = &ctx}},

    .num_init = 1,
    .projection[0] =
      {
        .proj_id = GKYL_PROJ_VLASOV_LTE,
        .density = evalElcDensityInit,
        .ctx_density = &ctx,
        .temp = evalTempInit,
        .ctx_temp = &ctx,
        .V_drift = evalElcVDriftInit,
        .ctx_V_drift = &ctx,
        .correct_all_moms = true,
        .use_last_converged = true,
      },

    // Adaptive pair injection, rescaled by the density above the threshold in both species.
    .source =
      {
        .source_id = GKYL_PROJ_ADAPT_DENSITY_SOURCE,
        // Particles moving away from the star breed outgoing pairs, particles moving
        // toward the star breed pairs that bombard the surface.
        .num_cross_source = 2,
        .source_with = {"pos", "pos"},
        .source_with_v_thresh = {ctx.p_thresh, ctx.p_thresh},
        .source_with_f_thresh = {ctx.f_thresh, ctx.f_thresh},
        .source_with_upper_half = {true, false},
        .source_with_proj = {0, 1},
        .write_source = true,
        .filter = true,
        .num_filters = ctx.num_filters,
        .num_sources = 2,
        .projection[0] =
          {
            .proj_id = GKYL_PROJ_VLASOV_LTE,
            .density = evalDensitySource,
            .ctx_density = &ctx,
            .temp = evalTempInit,
            .ctx_temp = &ctx,
            .V_drift = evalVDriftSourceUp,
            .ctx_V_drift = &ctx,
            .correct_all_moms = true,
            .use_last_converged = true,
          },
        .projection[1] =
          {
            .proj_id = GKYL_PROJ_VLASOV_LTE,
            .density = evalDensitySource,
            .ctx_density = &ctx,
            .temp = evalTempInit,
            .ctx_temp = &ctx,
            .V_drift = evalVDriftSourceLo,
            .ctx_V_drift = &ctx,
            .correct_all_moms = true,
            .use_last_converged = true,
          },
      },

    // The LTE moments are left out: the source diagnostics would divide by the zero
    // source before any particle has passed the threshold.
    .num_diag_moments = 2,
    .diag_moments = {GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1},
  };

  // Positrons.
  struct gkyl_vlasov_kinetic_species pos = {
    .model_id = GKYL_MODEL_SR,
    .bcx = {.lower = {.type = GKYL_SPECIES_FIXED_FUNC}, .upper = {.type = GKYL_SPECIES_FIXED_FUNC}},
    .app_accel = evalAppliedAccel,
    .app_accel_ctx = &ctx,
    .lower = {-1.0},
    .upper = {1.0},
    .cells = {NPX},
    .mapc2p_vel = {{.mapc2p_vel_func = mapc2p_px, .mapc2p_vel_ctx = &ctx}},

    .num_init = 1,
    .projection[0] =
      {
        .proj_id = GKYL_PROJ_VLASOV_LTE,
        .density = evalPosDensityInit,
        .ctx_density = &ctx,
        .temp = evalTempInit,
        .ctx_temp = &ctx,
        .V_drift = evalPosVDriftInit,
        .ctx_V_drift = &ctx,
        .correct_all_moms = true,
        .use_last_converged = true,
      },

    // Adaptive pair injection, rescaled by the density above the threshold in both species.
    .source =
      {
        .source_id = GKYL_PROJ_ADAPT_DENSITY_SOURCE,
        // Particles moving away from the star breed outgoing pairs, particles moving
        // toward the star breed pairs that bombard the surface.
        .num_cross_source = 2,
        .source_with = {"elc", "elc"},
        .source_with_v_thresh = {ctx.p_thresh, ctx.p_thresh},
        .source_with_f_thresh = {ctx.f_thresh, ctx.f_thresh},
        .source_with_upper_half = {true, false},
        .source_with_proj = {0, 1},
        .write_source = true,
        .filter = true,
        .num_filters = ctx.num_filters,
        .num_sources = 2,
        .projection[0] =
          {
            .proj_id = GKYL_PROJ_VLASOV_LTE,
            .density = evalDensitySource,
            .ctx_density = &ctx,
            .temp = evalTempInit,
            .ctx_temp = &ctx,
            .V_drift = evalVDriftSourceUp,
            .ctx_V_drift = &ctx,
            .correct_all_moms = true,
            .use_last_converged = true,
          },
        .projection[1] =
          {
            .proj_id = GKYL_PROJ_VLASOV_LTE,
            .density = evalDensitySource,
            .ctx_density = &ctx,
            .temp = evalTempInit,
            .ctx_temp = &ctx,
            .V_drift = evalVDriftSourceLo,
            .ctx_V_drift = &ctx,
            .correct_all_moms = true,
            .use_last_converged = true,
          },
      },

    // The LTE moments are left out: the source diagnostics would divide by the zero
    // source before any particle has passed the threshold.
    .num_diag_moments = 2,
    .diag_moments = {GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1},
  };

  // Field.
  struct gkyl_vlasov_field field = {
    .epsilon0 = ctx.epsilon0,
    .mu0 = ctx.mu0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,

    .init = evalFieldInit,
    .ctx = &ctx,

    .app_current = evalAppliedCurrent,
    .app_current_ctx = &ctx,
  };

  // Vlasov-Maxwell app.
  struct gkyl_vm app_inp = {
    .cdim = 1,
    .vdim = 1,
    .lower = {0.0},
    .upper = {ctx.Lx},
    .cells = {NX},

    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,
    .cfl_frac = ctx.cfl_frac,

    .num_periodic_dir = 0,
    .periodic_dirs = {},

    .num_species = 2,
    .species =
      {
        {
          .name = "elc",
          .charge = ctx.charge_elc,
          .mass = ctx.mass_elc,
          .type = GKYL_SPECIES_VLASOV,
          .kinetic = elc,
        },
        {
          .name = "pos",
          .charge = ctx.charge_pos,
          .mass = ctx.mass_pos,
          .type = GKYL_SPECIES_VLASOV,
          .kinetic = pos,
        },
      },

    .field = field,

    .parallelism = {.use_gpu = app_args.use_gpu, .cuts = {app_args.cuts[0]}, .comm = comm},
  };

  // Create app object.
  // Set app output name from the executable name (argv[0]).
  snprintf(app_inp.name, sizeof(app_inp.name), "%s", app_args.app_name);
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&app_inp);

  // Initial and final simulation times.
  double t_curr = 0.0, t_end = ctx.t_end;

  // Initialize simulation.
  int frame_curr = 0;
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status =
      gkyl_vlasov_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_vlasov_app_cout(
        app, stderr, "*** Failed to read restart file! (%s)\n",
        gkyl_array_rio_status_msg(status.io_status)
      );
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_vlasov_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_vlasov_app_cout(app, stdout, " at time = %g\n", t_curr);
  } else {
    gkyl_vlasov_app_apply_ic(app, t_curr);
  }

  // Create trigger for field energy.
  int field_energy_calcs = ctx.field_energy_calcs;
  struct gkyl_tm_trigger fe_trig = {
    .dt = t_end / field_energy_calcs,
    .tcurr = t_curr,
    .curr = frame_curr,
  };

  calc_field_energy(&fe_trig, app, t_curr, false);

  // Create trigger for integrated moments.
  int integrated_mom_calcs = ctx.integrated_mom_calcs;
  struct gkyl_tm_trigger im_trig = {
    .dt = t_end / integrated_mom_calcs,
    .tcurr = t_curr,
    .curr = frame_curr,
  };

  calc_integrated_mom(&im_trig, app, t_curr, false);

  // Create trigger for integrated L2 norm of the distribution function.
  int integrated_L2_f_calcs = ctx.integrated_L2_f_calcs;
  struct gkyl_tm_trigger l2f_trig = {
    .dt = t_end / integrated_L2_f_calcs,
    .tcurr = t_curr,
    .curr = frame_curr,
  };

  calc_integrated_L2_f(&l2f_trig, app, t_curr, false);

  // Create trigger for IO.
  int num_frames = ctx.num_frames;
  struct gkyl_tm_trigger io_trig = {
    .dt = t_end / num_frames,
    .tcurr = frame_curr * (t_end / num_frames),
    .curr = frame_curr,
  };

  write_data(&io_trig, app, t_curr, false);

  // Compute initial guess of maximum stable time-step.
  double dt = t_end - t_curr;

  // The requested time-step is shortened near the end of the simulation so that
  // the final step lands exactly on t_end.
  bool is_dt_clipped = false; // Was the requested dt shortened below the stable dt?
  bool is_last_step = true; // Does the requested dt reach t_end?

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = ctx.dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx.num_failures_max;

  long step = 1;
  while ((t_curr < t_end) && (step <= app_args.num_steps)) {
    gkyl_vlasov_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = gkyl_vlasov_update(app, dt);
    gkyl_vlasov_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

    if (!status.success) {
      gkyl_vlasov_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }

    // Only a step that took the full requested dt counts as shortened/final.
    bool took_requested_dt = status.dt_actual == dt;
    bool was_dt_clipped = is_dt_clipped && took_requested_dt;
    if (is_last_step && took_requested_dt) {
      // Avoid round-off leaving t_curr just short of t_end.
      t_curr = t_end;
      // Trigger times are accumulated sums and can exceed t_end by round-off.
      // Snap them so the final frame and diagnostics are still produced.
      snap_trigger_to_t_end(&fe_trig, t_end);
      snap_trigger_to_t_end(&im_trig, t_end);
      snap_trigger_to_t_end(&l2f_trig, t_end);
      snap_trigger_to_t_end(&io_trig, t_end);
    } else {
      t_curr += status.dt_actual;
    }

    // Request the next time-step. If the remaining time fits in one stable step,
    // take exactly the remaining time; if it fits in less than two, split it into
    // two equal steps so the final step is never a sliver of the stable dt.
    double t_left = t_end - t_curr;
    dt = status.dt_suggested;
    is_dt_clipped = false;
    is_last_step = false;
    if (t_left <= dt) {
      dt = t_left;
      is_dt_clipped = true;
      is_last_step = true;
    } else if (t_left < 2.0 * dt) {
      dt = 0.5 * t_left;
      is_dt_clipped = true;
    }

    calc_field_energy(&fe_trig, app, t_curr, false);
    calc_integrated_mom(&im_trig, app, t_curr, false);
    calc_integrated_L2_f(&l2f_trig, app, t_curr, false);
    write_data(&io_trig, app, t_curr, false);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    } else if (!was_dt_clipped && status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_vlasov_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_vlasov_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_vlasov_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_vlasov_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_vlasov_app_cout(
          app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max
        );

        calc_field_energy(&fe_trig, app, t_curr, true);
        calc_integrated_mom(&im_trig, app, t_curr, true);
        calc_integrated_L2_f(&l2f_trig, app, t_curr, true);
        write_data(&io_trig, app, t_curr, true);

        break;
      }
    } else {
      num_failures = 0;
    }

    step += 1;
  }

  calc_field_energy(&fe_trig, app, t_curr, false);
  calc_integrated_mom(&im_trig, app, t_curr, false);
  calc_integrated_L2_f(&l2f_trig, app, t_curr, false);
  write_data(&io_trig, app, t_curr, false);
  gkyl_vlasov_app_stat_write(app);

  struct gkyl_vlasov_stat stat = gkyl_vlasov_app_stat(app);

  gkyl_vlasov_app_cout(app, stdout, "\n");
  gkyl_vlasov_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_vlasov_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_vlasov_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_vlasov_app_cout(
      app, stdout, "  Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]
    );
    gkyl_vlasov_app_cout(
      app, stdout, "  Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]
    );
  }
  gkyl_vlasov_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_vlasov_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_vlasov_app_cout(
    app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm
  );
  gkyl_vlasov_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_vlasov_app_cout(
    app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm
  );
  gkyl_vlasov_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

  gkyl_vlasov_app_cout(app, stdout, "Number of write calls %ld\n", stat.n_io);
  double io_tm =
    stat.field_io_tm + stat.species_io_tm + stat.field_diag_io_tm + stat.species_diag_io_tm;
  gkyl_vlasov_app_cout(app, stdout, "IO time took %g secs \n", io_tm);

freeresources:
  // Free resources after simulation completion.
  gkyl_comm_release(comm);
  gkyl_vlasov_app_release(app);

mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif

  return 0;
}
