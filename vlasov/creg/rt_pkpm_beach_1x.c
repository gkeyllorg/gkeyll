// Propagation into a plasma wave beach with PKPM electrons in a magnetic field along the direction
// of propagation (1x, normalized parallel velocity w = v_par/xi), the magnetized counterpart of
// vlasov/creg/rt_dg_5m_mom_beach_1x.c. The plasma frequency falls as omega_p deltaT = 25 (1 - x)^5
// with deltaT = 0.01 / c, on a uniform background with omega_p deltaT = 0.025, and an antenna
// current J_y = -J_0 sin(omega_drive t) on 1 - 1/400 < x < 1 drives a wave with omega_drive deltaT
// = pi / 10 into the density gradient. The external field B_x gives Omega_ce = 0.5 omega_drive: the
// wave splits into its two circular polarizations, with the cutoff of the right-handed one at x =
// 0.611 (omega_p^2 = omega (omega - Omega_ce)) and of the left-handed one at x = 0.566 (omega_p^2 =
// omega (omega + Omega_ce)), and its plane of polarization rotates. J_0 is such that the wave the
// antenna radiates to each side in free space has a magnetic field of 0.02 B_x, so the magnetic
// unit vector changes in time up to the edge. The upper boundary of the field, next to the antenna,
// is absorbing (GKYL_FIELD_ABSORB). The fluid-field coupling is the staggered operator split with
// the implicit fluid-Maxwell source step, which takes the applied current. The electrons (1 eV)
// are held about their initial state (the fluid state evaluated at the nodes;
// First_Fourier_Harmonic_Notes/background_density_gradient.tex): the app derives the body force
// along b that holds the parallel pressure gradient, v_te^2 d ln n/dx away from the boundary
// cells, and subtracts what the discretization leaves. Without the hold the pressure gradient
// pushes the electrons until their space charge stops them, and E_x rings about the ambipolar
// field -(T_e/e) d ln n/dx at the local plasma frequency, up to 56 V/m.
//
// Figures of merit at t = 5e-9 (200 x 8 cells, serendipity p1, 1000 steps), from the cell averages:
// the maximum of |E_perp| = sqrt(E_y^2 + E_z^2) is 2.79e5 V/m at x = 0.6275, 1.74 times the
// amplitude mu_0 c J_0 / 800 the antenna radiates to each side in free space; the maxima of |E_y|
// and |E_z| are 1.33 and 1.17 times that amplitude. |E_perp| stays below 0.7% of its maximum for x
// < 0.466, beyond both cutoffs, where |E_x| stays below 0.02 V/m (22 V/m without the hold) and the
// density changes by less than 1e-8; in the wave E_x reaches 3.2e3 V/m. The norm of the
// first-harmonic blocks of the distribution is at most 0.42 of that of F (at x = 0.9075), and 0.24
// in the cell at the absorbing edge.

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
  double pi;

  // Physical constants (using non-normalized physical units).
  double epsilon0; // Permittivity of free space.
  double mu0; // Permeability of free space.
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.

  double Te; // Electron temperature (energy units).
  double omega_c_frac; // Electron cyclotron frequency over the drive frequency.
  double dB_B0; // Amplitude of the wave the antenna radiates in free space, over B0.
  double omega_p_floor; // Plasma frequency of the uniform background (times deltaT).

  // Derived physical quantities (using non-normalized physical units).
  double light_speed; // Speed of light.
  double B0; // External magnetic field (x-direction).
  double J0; // Antenna current density.

  // Simulation parameters.
  int Nx; // Cell count (configuration space: x-direction).
  int Nvx; // Cell count (velocity space: normalized parallel velocity w).
  double Lx; // Domain size (configuration space: x-direction).
  double Lx100; // Domain size over 100 (x-direction).
  double x_last_edge; // Lower edge of the antenna, which reaches the upper boundary.
  double vx_max; // Domain boundary (velocity space: w = v_par/xi).
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
  double pi = M_PI;

  // Physical constants (using non-normalized physical units).
  double epsilon0 = 8.854187817620389850536563031710750260608e-12; // Permittivity of free space.
  double mu0 = 12.56637061435917295385057353311801153679e-7; // Permeability of free space.
  double mass_elc = 9.10938215e-31; // Electron mass.
  double charge_elc = -1.602176487e-19; // Electron charge.

  double Te = 1.602176487e-19; // Electron temperature (1 eV).
  double omega_c_frac = 0.5; // Electron cyclotron frequency over the drive frequency.
  double dB_B0 = 0.02; // Amplitude of the wave the antenna radiates in free space, over B0.
  double omega_p_floor = 0.025; // Plasma frequency of the uniform background (times deltaT).

  // Derived physical quantities (using non-normalized physical units).
  double light_speed = 1.0 / sqrt(mu0 * epsilon0); // Speed of light.

  // Simulation parameters.
  int Nx = 200; // Cell count (configuration space: x-direction).
  int Nvx = 8; // Cell count (velocity space: normalized parallel velocity w).
  double Lx = 1.0; // Domain size (configuration space: x-direction).
  double Lx100 = Lx / 100.0; // Domain size over 100 (x-direction).
  double x_last_edge =
    Lx - Lx / 400.0; // Lower edge of the antenna, which reaches the upper boundary.
  double vx_max = 8.0; // Domain boundary (velocity space: w = v_par/xi).
  int poly_order = 1; // Polynomial order.
  double cfl_frac = 0.9; // CFL coefficient.

  double deltaT = Lx100 / light_speed; // Arbitrary constant, with units of time.
  double factor =
    deltaT * deltaT * charge_elc * charge_elc /
    (mass_elc * epsilon0); // Numerical factor for calculation of electron number density.
  double omega_drive = pi / 10.0 / deltaT; // Drive current angular frequency.

  double B0 = omega_c_frac * omega_drive * mass_elc / fabs(charge_elc); // External magnetic field.
  // A current sheet of width Lx - x_last_edge radiates B = mu0 J0 (Lx - x_last_edge) / 2 to each
  // side in free space.
  double J0 = 2.0 * dB_B0 * B0 / (mu0 * (Lx - x_last_edge)); // Antenna current density.

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
    .omega_c_frac = omega_c_frac,
    .dB_B0 = dB_B0,
    .omega_p_floor = omega_p_floor,
    .light_speed = light_speed,
    .B0 = B0,
    .J0 = J0,
    .Nx = Nx,
    .Nvx = Nvx,
    .Lx = Lx,
    .Lx100 = Lx100,
    .x_last_edge = x_last_edge,
    .vx_max = vx_max,
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

void
evalElcDistInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct beach_ctx *app = ctx;
  double w = xn[1];

  double pi = app->pi;

  double F = exp(-0.5 * w * w) / sqrt(2.0 * pi); // Normalized Maxwellian in w = v_par/xi.

  // Set normalized electron distribution functions [F, G] (G = (T_perp/T) F = F).
  fout[0] = F;
  fout[1] = F;

  // No first harmonic: the blocks [M0_perp (x, y, z), N_perp (x, y, z)] vanish.
  fout[2] = 0.0;
  fout[3] = 0.0;
  fout[4] = 0.0;
  fout[5] = 0.0;
  fout[6] = 0.0;
  fout[7] = 0.0;
}

void
evalElcFluidInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct beach_ctx *app = ctx;
  double x = xn[0];

  double mass_elc = app->mass_elc;
  double Te = app->Te;
  double omega_p_floor = app->omega_p_floor;
  double factor = app->factor;

  // Plasma frequency profile (times deltaT), on the uniform background.
  double omegaPdt = 25.0 * (1.0 - x) * (1.0 - x) * (1.0 - x) * (1.0 - x) * (1.0 - x);
  double ne =
    (omegaPdt * omegaPdt + omega_p_floor * omega_p_floor) / factor; // Electron number density.

  // Set electron fluid state [rho, rho ux, rho uy, rho uz, E].
  fout[0] = mass_elc * ne;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = 0.0;
  fout[4] = 1.5 * ne * Te;
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
evalExternalFieldInit(
  double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx
)
{
  struct beach_ctx *app = ctx;

  double B0 = app->B0;

  // Set external electric and magnetic fields.
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = B0;
  fout[4] = 0.0;
  fout[5] = 0.0;
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
  struct gkyl_vlasov_kinetic_species elc_kinetic = {
    .lower = {-ctx.vx_max},
    .upper = {ctx.vx_max},
    .cells = {NVX},

    .num_init = 1,
    .projection[0] = {.proj_id = GKYL_PROJ_FUNC, .func = evalElcDistInit, .ctx_func = &ctx},
    // The initial state is held as the background.
    .background = {.proj_id = GKYL_PROJ_FUNC, .func = evalElcDistInit, .ctx_func = &ctx},

    .bcx = {.lower = {.type = GKYL_SPECIES_COPY}, .upper = {.type = GKYL_SPECIES_COPY}},
  };

  struct gkyl_vlasov_fluid_species elc_fluid = {
    .init = evalElcFluidInit,
    .ctx = &ctx,
    // The initial state is held as the background, evaluated at the nodes (continuous).
    .background = evalElcFluidInit,
    .background_ctx = &ctx,
    .eval_on_nodes = true,

    .bcx = {GKYL_SPECIES_COPY, GKYL_SPECIES_COPY},
  };

  // Field.
  struct gkyl_vlasov_field field = {
    .epsilon0 = ctx.epsilon0,
    .mu0 = ctx.mu0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,

    .init = evalFieldInit,
    .ctx = &ctx,

    .ext_em = evalExternalFieldInit,
    .ext_em_ctx = &ctx,

    .app_current = evalAppCurrent,
    .app_current_ctx = &ctx,
    .app_current_evolve = true,

    .bcx = {GKYL_FIELD_COPY, GKYL_FIELD_ABSORB},
  };

  // PKPM app.
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

    .num_species = 1,
    .species = {{
      .name = "elc",
      .charge = ctx.charge_elc,
      .mass = ctx.mass_elc,
      .type = GKYL_SPECIES_PKPM,
      .kinetic = elc_kinetic,
      .fluid = elc_fluid,
    }},

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
  gkyl_vlasov_app_cout(
    app, stdout, "Fluid species RHS calc took %g secs\n", stat.fluid_species_rhs_tm
  );
  gkyl_vlasov_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_vlasov_app_cout(app, stdout, "Species PKPM vars took %g secs\n", stat.species_pkpm_vars_tm);
  gkyl_vlasov_app_cout(
    app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm
  );
  gkyl_vlasov_app_cout(
    app, stdout, "EM variables (bvar) calc took %g secs\n", stat.field_em_vars_tm
  );
  gkyl_vlasov_app_cout(
    app, stdout, "Current evaluation and accumulate took %g secs\n", stat.current_tm
  );
  gkyl_vlasov_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

  gkyl_vlasov_app_cout(app, stdout, "Number of write calls %ld\n", stat.n_io);
  gkyl_vlasov_app_cout(app, stdout, "IO time took %g secs \n", stat.species_io_tm);

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
