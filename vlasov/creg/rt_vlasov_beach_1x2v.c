// Propagation into a plasma wave beach with kinetic electrons and the Maxwell field (1x2v), the
// kinetic counterpart of vlasov/creg/rt_dg_5m_mom_beach_1x.c, about a held background density.
// The setup follows entry JE8 of Ammar's Simulation Journal
// (https://ammar-hakim.org/sj/je/je8/je8-plasmabeach.html), adapted from Section III. A. of D. N.
// Smithe (2007), "Finite-difference time-domain simulation of fusion plasmas at radiofrequency time
// scales", Physics of Plasmas 14 (5): 056104: the plasma frequency falls as omega_p deltaT = 25 (1
// - x)^5 with deltaT = 0.01 / c, here on a uniform background with omega_p deltaT = 0.025, and an
// antenna current J_y = -J_0 sin(omega_drive t) on 1 - 1/400 < x < 1 drives a wave with
// omega_drive deltaT = pi / 10 toward the cutoff omega_p = omega_drive at x = 0.583. The upper
// boundary of the field, next to the antenna, is absorbing (GKYL_FIELD_ABSORB). The electron
// current is explicit in the SSP-RK3 stages, so the time step resolves the plasma frequency of
// every cell (omega_p dt = 0.9 in the cell at x = 0), where the 5-moment beach, with its implicit
// source step, takes a quarter of the steps. The Maxwellian electrons (1 eV, velocity grid +-6
// v_te) are held about their initial state
// (First_Fourier_Harmonic_Notes/background_density_gradient.tex): the app derives the body force
// that holds them and subtracts what the discretization leaves. Without the hold the pressure
// gradient pushes the electrons until their space charge stops them, and E_x rings about the
// ambipolar field -(T_e/e) d ln n/dx at the local plasma frequency, up to 54 V/m. The response is
// linear in the antenna current: J_0 = 1 A/m^3 (1e-12 in the 5-moment beach) puts the electron
// drift at 2.4e-5 v_te, above the round-off of the current of a distribution function.
//
// Figures of merit at t = 5e-9 (100 x 8 x 8 cells, serendipity p2, 4062 steps), from the cell
// averages: on 0.3 < x < 0.95, E_y is that of the 5-moment beach, per unit antenna current, times
// 1.1105 to 3.9e-3 in L2. The factor is the 10/9 of its current that the 3-point quadrature of
// the applied current gives the antenna, which covers a quarter of the last cell. The maximum of
// |E_y| is 0.992 V/m at x = 0.625, 2.11 times the amplitude mu_0 c J_0 / 800 the antenna radiates
// to each side in free space (1.90 times that of the projected current, as in the 5-moment beach),
// and |E_y| stays below 0.3% of the maximum for x < 0.485, beyond the cutoff. |E_x| stays below
// 1e-7 V/m, the density at its initial value to 6e-11 and |u_x| below 3e-11 v_te. 12 x 12
// velocity cells give the same E_y to four digits.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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

struct beach_ctx {
  // Mathematical constants (dimensionless).
  double pi; // Pi.

  // Physical constants (using non-normalized physical units).
  double epsilon0; // Permittivity of free space.
  double mu0; // Permeability of free space.
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.

  double Te; // Electron temperature.
  double J0; // Reference current density (Amps / m^3).
  double omega_p_floor; // Plasma frequency of the uniform background (times deltaT).

  // Derived physical quantities (using non-normalized physical units).
  double light_speed; // Speed of light.
  double vte; // Electron thermal velocity.

  // Simulation parameters.
  int Nx; // Cell count (configuration space: x-direction).
  int Nvx; // Cell count (velocity space: vx-direction).
  int Nvy; // Cell count (velocity space: vy-direction).
  double Lx; // Domain size (configuration space: x-direction).
  double vx_max; // Domain boundary (velocity space: vx-direction).
  double vy_max; // Domain boundary (velocity space: vy-direction).
  double Lx100; // Domain size over 100 (x-direction).
  double
    x_last_edge; // Lower edge of the antenna, which reaches the upper boundary (low density side).
  int poly_order; // Polynomial order.
  double cfl_frac; // CFL coefficient.

  double deltaT; // Arbitrary constant, with units of time.
  double factor; // Numerical factor for calculation of electron number density.
  double omega_drive; // Drive current angular frequency.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  int field_energy_calcs; // Number of times to calculate field energy.
  int integrated_mom_calcs; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct beach_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI; // Pi.

  // Physical constants (using non-normalized physical units).
  double epsilon0 = 8.854187817620389850536563031710750260608e-12; // Permittivity of free space.
  double mu0 = 12.56637061435917295385057353311801153679e-7; // Permeability of free space.
  double mass_elc = 9.10938215e-31; // Electron mass.
  double charge_elc = -1.602176487e-19; // Electron charge.

  double Te = 1.0 * (-charge_elc); // Electron temperature (1 eV).
  double J0 = 1.0; // Reference current density (Amps / m^3).
  double omega_p_floor = 0.025; // Plasma frequency of the uniform background (times deltaT).

  // Derived physical quantities (using non-normalized physical units).
  double light_speed = 1.0 / sqrt(mu0 * epsilon0); // Speed of light.
  double vte = sqrt(Te / mass_elc); // Electron thermal velocity.

  // Simulation parameters.
  int Nx = 100; // Cell count (configuration space: x-direction).
  int Nvx = 8; // Cell count (velocity space: vx-direction).
  int Nvy = 8; // Cell count (velocity space: vy-direction).
  double Lx = 1.0; // Domain size (configuration space: x-direction).
  double vx_max = 6.0 * vte; // Domain boundary (velocity space: vx-direction).
  double vy_max = 6.0 * vte; // Domain boundary (velocity space: vy-direction).
  double Lx100 = Lx / 100.0; // Domain size over 100 (x-direction).
  double x_last_edge =
    Lx -
    Lx / 400.0; // Lower edge of the antenna, which reaches the upper boundary (low density side).
  int poly_order = 2; // Polynomial order.
  double cfl_frac = 0.9; // CFL coefficient.

  double deltaT = Lx100 / light_speed; // Arbitrary constant, with units of time.
  double factor =
    deltaT * deltaT * charge_elc * charge_elc /
    (mass_elc * epsilon0); // Numerical factor for calculation of electron number density.
  double omega_drive = pi / 10.0 / deltaT; // Drive current angular frequency.

  double t_end = 5.0e-9; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs =
    INT_MAX; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct beach_ctx ctx = {
    .pi = pi,
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .Te = Te,
    .J0 = J0,
    .omega_p_floor = omega_p_floor,
    .light_speed = light_speed,
    .vte = vte,
    .Nx = Nx,
    .Nvx = Nvx,
    .Nvy = Nvy,
    .Lx = Lx,
    .vx_max = vx_max,
    .vy_max = vy_max,
    .Lx100 = Lx100,
    .x_last_edge = x_last_edge,
    .poly_order = poly_order,
    .cfl_frac = cfl_frac,
    .deltaT = deltaT,
    .factor = factor,
    .omega_drive = omega_drive,
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

// Electron distribution function: a Maxwellian at rest with the plasma frequency profile (times
// deltaT) on the uniform background. The initial condition and the held background.
void
evalElcInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct beach_ctx *app = ctx;
  double x = xn[0], vx = xn[1], vy = xn[2];

  double pi = app->pi;
  double factor = app->factor;
  double omega_p_floor = app->omega_p_floor;
  double vte = app->vte;

  // Plasma frequency profile (times deltaT), on the uniform background.
  double omegaPdt = 25.0 * (1.0 - x) * (1.0 - x) * (1.0 - x) * (1.0 - x) * (1.0 - x);
  double ne =
    (omegaPdt * omegaPdt + omega_p_floor * omega_p_floor) / factor; // Electron number density.

  // Set electron distribution function.
  fout[0] = ne / (2.0 * pi * vte * vte) * exp(-(vx * vx + vy * vy) / (2.0 * vte * vte));
}

void
evalFieldInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  // Set electric field, magnetic field and correction potentials.
  for (int i = 0; i < 8; ++i) {
    fout[i] = 0.0;
  }
}

void
evalAppCurrent(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct beach_ctx *app = ctx;
  double x = xn[0];

  double J0 = app->J0;
  double x_last_edge = app->x_last_edge;
  double omega_drive = app->omega_drive;

  // Applied current (y-direction): the antenna next to the upper boundary.
  double app_y = 0.0;
  if (x > x_last_edge) {
    app_y = -J0 * sin(omega_drive * t);
  }

  // Set applied current.
  fout[0] = 0.0;
  fout[1] = app_y;
  fout[2] = 0.0;
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

  struct beach_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int NVX = APP_ARGS_CHOOSE(app_args.vcells[0], ctx.Nvx);
  int NVY = APP_ARGS_CHOOSE(app_args.vcells[1], ctx.Nvy);

  // Field.
  struct gkyl_vlasov_field field = {
    .epsilon0 = ctx.epsilon0,
    .mu0 = ctx.mu0,

    .init = evalFieldInit,
    .ctx = &ctx,

    .app_current = evalAppCurrent,
    .app_current_ctx = &ctx,
    .app_current_evolve = true,

    .bcx = {GKYL_FIELD_COPY, GKYL_FIELD_ABSORB},
  };

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

  // Electrons.
  struct gkyl_vlasov_kinetic_species elc = {
    .lower = {-ctx.vx_max, -ctx.vy_max},
    .upper = {ctx.vx_max, ctx.vy_max},
    .cells = {NVX, NVY},

    .num_init = 1,
    .projection[0] = {.proj_id = GKYL_PROJ_FUNC, .func = evalElcInit, .ctx_func = &ctx},
    // The initial state is held as the background.
    .background = {.proj_id = GKYL_PROJ_FUNC, .func = evalElcInit, .ctx_func = &ctx},

    .bcx = {.lower = {.type = GKYL_SPECIES_COPY}, .upper = {.type = GKYL_SPECIES_COPY}},

    .num_diag_moments = 3,
    .diag_moments = {GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2},
  };

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

  // Vlasov-Maxwell app.
  struct gkyl_vm app_inp = {

    .cdim = 1,
    .vdim = 2,
    .lower = {0.0},
    .upper = {ctx.Lx},
    .cells = {NX},

    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,
    .cfl_frac = ctx.cfl_frac,

    .num_periodic_dir = 0,
    .periodic_dirs = {},

    .num_species = 1,
    .species =
      {
        {
          .name = "elc",
          .type = GKYL_SPECIES_VLASOV,
          .charge = ctx.charge_elc,
          .mass = ctx.mass_elc,
          .kinetic = elc,
        },
      },

    .field = field,

    .parallelism = {.use_gpu = app_args.use_gpu, .cuts = {cuts[0]}, .comm = comm},
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
    .tcurr = gkyl_tm_trigger_anchor(t_end / num_frames, t_curr),
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
