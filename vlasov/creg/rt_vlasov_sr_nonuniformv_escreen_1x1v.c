// Electric-field screening by pair creation (special-relativistic Vlasov-Maxwell, 1x1v): the
// setup of Tolman, Philippov and Timokhin (2022, ApJL 933, L37) as adapted in the Gkeyll
// relativistic Vlasov paper. An electron-positron plasma (n = 1 per species, T = mc^2) sits in a
// uniform electric field E0 = 50 while cold pairs (T = 0.1 mc^2) drifting at u = 5 are injected at
// 0.5 per species per unit time. The injected pairs screen the field, which then oscillates at the
// relativistic plasma frequency with an amplitude damped by the fresh pairs spun up at every zero
// crossing. The momentum grid is linear (cell 0.4) near the origin and exponential out to p = 400
// so the particles accelerated during screening (p ~ 2 sqrt(xi)/3 = 330, xi = E0^3/S = 2.5e5)
// stay on the grid. Reference: a Lagrangian k = 0 solution of the same equations (scratch
// kzero.py). Figures of merit (16 x 96 cells, 7214 steps): screening at t = 8.42 (reference 8.42,
// cold-beam estimate (sqrt(n0^2 + S E0) - n0)/S = 8.2), peak field 20.2, 12.4 and 9.4 at
// t = 12.2, 19.5 and 24.8 (reference 20.1, 12.1 and 9.1), envelope decay rate 0.0613 over the
// first three peaks against 0.0639 for the reference and 0.0624 for the PIC fit of Tolman et al.
// (a = 0.43 xi^-0.47 per t0 = mc/(e E0)); the rms field tracks the reference to 2.3% of E0
// (0.6% in L1); the injected mass is exact.
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

struct escreen_sr_ctx {
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using normalized code units).
  double epsilon0; // Permittivity of free space.
  double mu0; // Permeability of free space.
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.
  double mass_pos; // Positron mass.
  double charge_pos; // Positron charge.

  double n0; // Initial number density of each species.
  double T0; // Initial temperature (units of mc^2).
  double E0; // Initial electric field.
  double n_src; // Pair injection rate (number density per unit time, each species).
  double T_src; // Temperature of the injected pairs (units of mc^2).
  double u_src; // Drift (four-) velocity of the injected pairs (x-direction).

  double perturb; // Relative amplitude of the seeded density modes.
  int num_modes; // Number of seeded density modes.

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

struct escreen_sr_ctx
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

  double n0 = 1.0; // Initial number density of each species.
  double T0 = 1.0; // Initial temperature (units of mc^2).
  double E0 = 50.0; // Initial electric field.
  double n_src = 0.5; // Pair injection rate (number density per unit time, each species).
  double T_src = 0.1; // Temperature of the injected pairs (units of mc^2).
  double u_src = 5.0; // Drift (four-) velocity of the injected pairs (x-direction).

  double perturb = 1.0e-3; // Relative amplitude of the seeded density modes.
  int num_modes = 4; // Number of seeded density modes.

  // Simulation parameters.
  int Nx = 16; // Cell count (configuration space: x-direction).
  int Npx = 96; // Cell count (momentum space: px-direction).
  double Lx = 1.0; // Domain size (configuration space: x-direction).
  double px_max = 400.0; // Domain boundary (momentum space: px-direction).
  double px_lin = 0.14; // Momentum map: cell size of the linear part of the map at the origin.
  int poly_order = 2; // Polynomial order.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 30.0; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs =
    INT_MAX; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct escreen_sr_ctx ctx = {
    .pi = pi,
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .mass_pos = mass_pos,
    .charge_pos = charge_pos,
    .n0 = n0,
    .T0 = T0,
    .E0 = E0,
    .n_src = n_src,
    .T_src = T_src,
    .u_src = u_src,
    .perturb = perturb,
    .num_modes = num_modes,
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

// Seeded density modes: sum_j (perturb / j) cos(j k x + 2 pi j 0.3), j = 1 .. num_modes.
static double
density_modes(const struct escreen_sr_ctx *app, double x)
{
  double pi = app->pi;
  double kx = 2.0 * pi / app->Lx;
  double sum = 0.0;
  for (int j = 1; j <= app->num_modes; j++) {
    sum += (app->perturb / j) * cos((j * kx * x) + (2.0 * pi * j * 0.3));
  }
  return sum;
}

// The electric field that balances the seeded charge density through Gauss's law.
static double
field_modes(const struct escreen_sr_ctx *app, double x)
{
  double pi = app->pi;
  double kx = 2.0 * pi / app->Lx;
  double sum = 0.0;
  for (int j = 1; j <= app->num_modes; j++) {
    sum += (app->perturb / j) * sin((j * kx * x) + (2.0 * pi * j * 0.3)) / (j * kx);
  }
  return 2.0 * app->n0 * sum;
}

void
evalElcDensityInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct escreen_sr_ctx *app = ctx;
  double x = xn[0];

  // Electron number density with the seeded modes removed.
  fout[0] = app->n0 * (1.0 - density_modes(app, x));
}

void
evalPosDensityInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct escreen_sr_ctx *app = ctx;
  double x = xn[0];

  // Positron number density with the seeded modes added.
  fout[0] = app->n0 * (1.0 + density_modes(app, x));
}

void
evalTempInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct escreen_sr_ctx *app = ctx;

  // Initial temperature.
  fout[0] = app->T0;
}

void
evalVDriftInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  // The initial plasma is at rest.
  fout[0] = 0.0;
}

void
evalDensitySource(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct escreen_sr_ctx *app = ctx;

  // Pair injection rate. The LTE projection takes the rest-frame density, so divide by the
  // Lorentz factor of the drift to inject n_src pairs per unit time in the lab frame.
  fout[0] = app->n_src / sqrt(1.0 + (app->u_src * app->u_src));
}

void
evalTempSource(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct escreen_sr_ctx *app = ctx;

  // Temperature of the injected pairs.
  fout[0] = app->T_src;
}

void
evalVDriftSource(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct escreen_sr_ctx *app = ctx;

  // Drift (four-) velocity of the injected pairs.
  fout[0] = app->u_src;
}

void
evalFieldInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct escreen_sr_ctx *app = ctx;
  double x = xn[0];

  double Ex = app->E0 + field_modes(app, x); // Total electric field (x-direction).

  // Set electric field.
  fout[0] = Ex, fout[1] = 0.0, fout[2] = 0.0;
  // Set magnetic field.
  fout[3] = 0.0, fout[4] = 0.0, fout[5] = 0.0;
  // Set auxiliary fields.
  fout[6] = 0.0, fout[7] = 0.0;
}

void
mapc2p_px(double t, const double *GKYL_RESTRICT vc, double *GKYL_RESTRICT vp, void *ctx)
{
  struct escreen_sr_ctx *app = ctx;
  double px_c = vc[0];

  double px_max = app->px_max;
  double px_lin = app->px_lin;
  int Npx = app->Npx;

  // Linear-to-exponential momentum map: cells of size px_lin at the origin, growing
  // exponentially to the domain boundary px_max.
  if (px_c < 0.0) {
    vp[0] = (px_lin * px_c * Npx) - exp(-px_c * log(px_max)) + 1.0;
  } else {
    vp[0] = (px_lin * px_c * Npx) + exp(px_c * log(px_max)) - 1.0;
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

  struct escreen_sr_ctx ctx = create_ctx(); // Context for initialization functions.

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
        .V_drift = evalVDriftInit,
        .ctx_V_drift = &ctx,
        .correct_all_moms = true,
        .use_last_converged = true,
      },

    // Pair injection.
    .source =
      {
        .source_id = GKYL_PROJ_SOURCE,
        .num_sources = 1,
        .projection[0] =
          {
            .proj_id = GKYL_PROJ_VLASOV_LTE,
            .density = evalDensitySource,
            .ctx_density = &ctx,
            .temp = evalTempSource,
            .ctx_temp = &ctx,
            .V_drift = evalVDriftSource,
            .ctx_V_drift = &ctx,
            .correct_all_moms = true,
            .use_last_converged = true,
          },
      },

    .num_diag_moments = 3,
    .diag_moments = {GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_LTE},
  };

  // Positrons.
  struct gkyl_vlasov_kinetic_species pos = {
    .model_id = GKYL_MODEL_SR,
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
        .V_drift = evalVDriftInit,
        .ctx_V_drift = &ctx,
        .correct_all_moms = true,
        .use_last_converged = true,
      },

    // Pair injection.
    .source =
      {
        .source_id = GKYL_PROJ_SOURCE,
        .num_sources = 1,
        .projection[0] =
          {
            .proj_id = GKYL_PROJ_VLASOV_LTE,
            .density = evalDensitySource,
            .ctx_density = &ctx,
            .temp = evalTempSource,
            .ctx_temp = &ctx,
            .V_drift = evalVDriftSource,
            .ctx_V_drift = &ctx,
            .correct_all_moms = true,
            .use_last_converged = true,
          },
      },

    .num_diag_moments = 3,
    .diag_moments = {GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_LTE},
  };

  // Field.
  struct gkyl_vlasov_field field = {
    .epsilon0 = ctx.epsilon0,
    .mu0 = ctx.mu0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,

    .init = evalFieldInit,
    .ctx = &ctx,
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

    .num_periodic_dir = 1,
    .periodic_dirs = {0},

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
