// Rayleigh-Taylor instability for a neutral gas with LBO collisions with the PKPM system of
// equations (2x, normalized parallel velocity w = v_par/xi).
// A gas of density 2 above one of density 1 in a gravitational acceleration 0.1 along -y, in the
// box [0, 1/4] x [0, 1/2] with reflecting walls at the bottom and the top and periodic in x, to t
// = 3.5. The interface is flat at y = 1/4, the pressure is hydrostatic with 0.01 at the top wall,
// and the perturbation is the velocity of the incompressible interface mode with one wavelength. A
// static uniform magnetic field along y, the direction of gravity, sets the field direction b.
// Gravity is an applied acceleration (explicit sources).
// Moderate collisions advanced explicitly, nu = 20.0: the mean free path of the hottest gas is 0.7
// cells. This is the collision frequency at which the diffusion step of the LBO operator is about
// the transport step; rt_pkpm_neut_bgk_rt_2x has the same collision frequency and grid.
// Figures of merit (12 x 24 x 8 cells) at t = 3.5, with nu dt = 0.017 on average: the integral of
// rho u_y^2 is 2.77e-06 at the end of the run (7.46e-07 at t = 0); mass is conserved to 3e-13.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_pkpm.h>
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

struct rt_ctx {
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using normalized code units).
  double epsilon0; // Permittivity of free space.
  double mu0; // Permeability of free space.
  double mass; // Neutral mass.
  double charge; // Neutral charge.

  double grav; // Gravitational acceleration (along -y).
  double rho_top; // Top mass density.
  double rho_bot; // Bottom mass density.
  double p_ref; // Pressure at the top wall.
  double pert_max; // Amplitude of the velocity perturbation.
  double nu; // Collision frequency.
  double B0; // Magnetic field along y (sets the field direction b, along gravity).

  // Simulation parameters.
  int Nx; // Cell count (configuration space: x-direction).
  int Ny; // Cell count (configuration space: y-direction).
  int Nvx; // Cell count (velocity space: normalized parallel velocity w).
  double Lx; // Domain size (configuration space: x-direction).
  double Ly; // Domain size (configuration space: y-direction).
  double vx_max; // Domain boundary (velocity space: w = v_par/xi).
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

struct rt_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  double epsilon0 = 1.0; // Permittivity of free space.
  double mu0 = 1.0; // Permeability of free space.
  double mass = 1.0; // Neutral mass.
  double charge = 0.0; // Neutral charge.

  double grav = 0.1; // Gravitational acceleration (along -y).
  double rho_top = 2.0; // Top mass density.
  double rho_bot = 1.0; // Bottom mass density.
  double p_ref = 0.01; // Pressure at the top wall.
  double pert_max = 0.01; // Amplitude of the velocity perturbation.
  double nu = 20.0; // Collision frequency.
  double B0 = 1.0; // Magnetic field along y (sets the field direction b, along gravity).

  // Simulation parameters.
  int Nx = 12; // Cell count (configuration space: x-direction).
  int Ny = 24; // Cell count (configuration space: y-direction).
  int Nvx = 8; // Cell count (velocity space: normalized parallel velocity w).
  double Lx = 0.25; // Domain size (configuration space: x-direction).
  double Ly = 0.5; // Domain size (configuration space: y-direction).
  double vx_max = 6.0; // Domain boundary (velocity space: w = v_par/xi).
  int poly_order = 1; // Polynomial order.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 3.5; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs =
    INT_MAX; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct rt_ctx ctx = {
    .pi = pi,
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .mass = mass,
    .charge = charge,
    .grav = grav,
    .rho_top = rho_top,
    .rho_bot = rho_bot,
    .p_ref = p_ref,
    .pert_max = pert_max,
    .nu = nu,
    .B0 = B0,
    .Nx = Nx,
    .Ny = Ny,
    .Nvx = Nvx,
    .Lx = Lx,
    .Ly = Ly,
    .vx_max = vx_max,
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

void
evalNeutDistInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct rt_ctx *app = ctx;
  double w = xn[2];

  double pi = app->pi;

  double F = exp(-0.5 * w * w) / sqrt(2.0 * pi); // Normalized Maxwellian in w = v_par/xi.

  // Set normalized neutral distribution functions [F, G] (G = (T_perp/T) F = F).
  fout[0] = F;
  fout[1] = F;
}

void
evalNeutFluidInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct rt_ctx *app = ctx;
  double x = xn[0], y = xn[1];

  double pi = app->pi;
  double grav = app->grav;
  double rho_top = app->rho_top;
  double rho_bot = app->rho_bot;
  double p_ref = app->p_ref;
  double pert_max = app->pert_max;
  double Lx = app->Lx;
  double Ly = app->Ly;

  double kx = 2.0 * pi / Lx; // Wavenumber of the perturbation.
  double yloc = 0.5 * Ly; // Interface.

  double rho = 0.0;
  double pres = 0.0;
  if (y > yloc) {
    rho = rho_top; // Mass density (top).
    pres = p_ref + rho_top * grav * (Ly - y); // Hydrostatic pressure (top).
  } else {
    rho = rho_bot; // Mass density (bottom).
    pres = p_ref + rho_top * grav * (Ly - yloc) +
           rho_bot * grav * (yloc - y); // Hydrostatic pressure (bottom).
  }
  // Velocity of the incompressible interface mode.
  double decay = exp(-kx * fabs(y - yloc));
  double vx = (y > yloc ? 1.0 : -1.0) * pert_max * sin(kx * x) * decay;
  double vy = pert_max * cos(kx * x) * decay;
  double E = 0.5 * rho * (vx * vx + vy * vy) + 1.5 * pres; // Total energy density.

  // Set neutral fluid state [rho, rho ux, rho uy, rho uz, E].
  fout[0] = rho;
  fout[1] = rho * vx;
  fout[2] = rho * vy;
  fout[3] = 0.0;
  fout[4] = E;
}

void
evalAppAccel(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct rt_ctx *app = ctx;

  double grav = app->grav;

  // Set applied acceleration.
  fout[0] = 0.0;
  fout[1] = -grav;
  fout[2] = 0.0;
}

void
evalNu(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct rt_ctx *app = ctx;

  double nu = app->nu;

  // Set collision frequency.
  fout[0] = nu;
}

void
evalFieldInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct rt_ctx *app = ctx;

  double B0 = app->B0;

  // Set electric field, magnetic field and correction potentials.
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = 0.0;
  fout[4] = B0;
  fout[5] = 0.0;
  fout[6] = 0.0;
  fout[7] = 0.0;
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
write_data(struct gkyl_tm_trigger *iot, gkyl_pkpm_app *app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_write) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_pkpm_app_write(app, t_curr, frame);
    gkyl_pkpm_app_write_field_energy(app);
    gkyl_pkpm_app_write_integrated_mom(app);
    gkyl_pkpm_app_write_integrated_L2_f(app);
  }
}

void
calc_field_energy(struct gkyl_tm_trigger *fet, gkyl_pkpm_app *app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(fet, t_curr) || force_calc) {
    gkyl_pkpm_app_calc_field_energy(app, t_curr);
  }
}

void
calc_integrated_mom(struct gkyl_tm_trigger *imt, gkyl_pkpm_app *app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(imt, t_curr) || force_calc) {
    gkyl_pkpm_app_calc_integrated_mom(app, t_curr);
  }
}

void
calc_integrated_L2_f(struct gkyl_tm_trigger *l2t, gkyl_pkpm_app *app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(l2t, t_curr) || force_calc) {
    gkyl_pkpm_app_calc_integrated_L2_f(app, t_curr);
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

  struct rt_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], ctx.Ny);
  int NVX = APP_ARGS_CHOOSE(app_args.vcells[0], ctx.Nvx);

  int nrank = 1; // Number of processors in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif

  int ccells[] = {NX, NY};
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

  // Neutrals.
  struct gkyl_pkpm_species neut = {
    .name = "neut",
    .charge = ctx.charge,
    .mass = ctx.mass,
    .lower = {-ctx.vx_max},
    .upper = {ctx.vx_max},
    .cells = {NVX},

    .init_dist = evalNeutDistInit,
    .ctx_dist = &ctx,
    .init_fluid = evalNeutFluidInit,
    .ctx_fluid = &ctx,

    .app_accel = evalAppAccel,
    .app_accel_ctx = &ctx,

    .bcy = {GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT},

    .collisions = {.collision_id = GKYL_LBO_COLLISIONS, .self_nu = evalNu, .ctx = &ctx},
  };

  // Field.
  struct gkyl_pkpm_field field = {
    .epsilon0 = ctx.epsilon0,
    .mu0 = ctx.mu0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,

    .init = evalFieldInit,
    .ctx = &ctx,

    .is_static = true,
  };

  // PKPM app.
  struct gkyl_pkpm app_inp = {

    .cdim = 2,
    .vdim = 1,
    .lower = {0.0, 0.0},
    .upper = {ctx.Lx, ctx.Ly},
    .cells = {NX, NY},

    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,
    .cfl_frac = ctx.cfl_frac,

    .num_periodic_dir = 1,
    .periodic_dirs = {0},

    .num_species = 1,
    .species = {neut},

    .field = field,

    .use_explicit_source = true,

    .parallelism =
      {.use_gpu = app_args.use_gpu, .cuts = {app_args.cuts[0], app_args.cuts[1]}, .comm = comm},
  };
  // Create app object.
  // Set app output name from the executable name (argv[0]).
  snprintf(app_inp.name, sizeof(app_inp.name), "%s", app_args.app_name);
  gkyl_pkpm_app *app = gkyl_pkpm_app_new(&app_inp);

  // Initial and final simulation times.
  double t_curr = 0.0, t_end = ctx.t_end;

  // Initialize simulation.
  int frame_curr = 0;
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status =
      gkyl_pkpm_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_pkpm_app_cout(
        app, stderr, "*** Failed to read restart file! (%s)\n",
        gkyl_array_rio_status_msg(status.io_status)
      );
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_pkpm_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_pkpm_app_cout(app, stdout, " at time = %g\n", t_curr);
  } else {
    gkyl_pkpm_app_apply_ic(app, t_curr);
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
    gkyl_pkpm_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = gkyl_pkpm_update(app, dt);
    gkyl_pkpm_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

    if (!status.success) {
      gkyl_pkpm_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
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

      gkyl_pkpm_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_pkpm_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_pkpm_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_pkpm_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_pkpm_app_cout(
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
  gkyl_pkpm_app_stat_write(app);

  struct gkyl_pkpm_stat stat = gkyl_pkpm_app_stat(app);

  gkyl_pkpm_app_cout(app, stdout, "\n");
  gkyl_pkpm_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_pkpm_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_pkpm_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_pkpm_app_cout(
      app, stdout, "  Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]
    );
    gkyl_pkpm_app_cout(
      app, stdout, "  Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]
    );
  }
  gkyl_pkpm_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_pkpm_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_pkpm_app_cout(
    app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm
  );
  gkyl_pkpm_app_cout(
    app, stdout, "Fluid species RHS calc took %g secs\n", stat.fluid_species_rhs_tm
  );
  gkyl_pkpm_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_pkpm_app_cout(app, stdout, "Species PKPM vars took %g secs\n", stat.species_pkpm_vars_tm);
  gkyl_pkpm_app_cout(
    app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm
  );
  gkyl_pkpm_app_cout(app, stdout, "EM variables (bvar) calc took %g secs\n", stat.field_em_vars_tm);
  gkyl_pkpm_app_cout(
    app, stdout, "Current evaluation and accumulate took %g secs\n", stat.current_tm
  );
  gkyl_pkpm_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

  gkyl_pkpm_app_cout(app, stdout, "Number of write calls %ld\n", stat.n_io);
  gkyl_pkpm_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

freeresources:
  // Free resources after simulation completion.
  gkyl_comm_release(comm);
  gkyl_pkpm_app_release(app);

mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif

  return 0;
}
