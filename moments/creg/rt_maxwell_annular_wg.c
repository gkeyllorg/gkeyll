#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_euler.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

#include <rt_arg_parse.h>

struct annular_wg_ctx
{
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using normalized code units).
  double epsilon0; // Permittivity of free space.
  double mu0; // Permeability of free space.

  double w_mode; // Frequency of wave mode.
  int bessel_order; // Spherical Bessel function order.
  double a_coeff; // Spherical Bessel function coefficient (first kind).
  double b_coeff; // Spherical Bessel function coefficient (second kind).

  // Derived physical quantities (using normalized code units).
  double t_period; // Time period of wave mode.

  // Simulation parameters.
  int Nr; // Cell count (radial direction).
  int Ntheta; // Cell count (angular direction).
  double Lr; // Domain size (radial direction).
  double Ltheta; // Domain size (angular direction).
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  int field_energy_writes; // Number of times to output field energy.
  int integrated_mom_writes; // Number of times to output integrated moments.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct annular_wg_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  double epsilon0 = 1.0; // Permittivity of free space.
  double mu0 = 1.0; // Permeability of free space.

  double w_mode = 1.19318673737701; // Frequency of wave mode.
  int bessel_order = 2; // Spherical Bessel function order.
  double a_coeff = 1.0; // Spherical Bessel function coefficient (first kind).
  double b_coeff = 0.9904672582498093; // Spherical Bessel function coefficient (second kind).

  // Derived physical quantities (using normalized code units).
  double t_period = 2.0 * pi / w_mode; // Time period of wave mode.

  // Simulation parameters.
  int Nr = 32; // Cell count (radial direction).
  int Ntheta = 32 * 6; // Cell count (angular direction).
  double Lr = 3.0; // Domain size (radial direction).
  double Ltheta = 2.0 * pi; // Domain size (angular direction).
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 2.0 * t_period; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int field_energy_writes = 1; // Number of times to output field energy.
  int integrated_mom_writes = 1; // Number of times to output integrated moments.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct annular_wg_ctx ctx = {
    .pi = pi,
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .w_mode = w_mode,
    .bessel_order = bessel_order,
    .a_coeff = a_coeff,
    .b_coeff = b_coeff,
    .t_period = t_period,
    .Nr = Nr,
    .Ntheta = Ntheta,
    .Lr = Lr,
    .Ltheta = Ltheta,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .num_frames = num_frames,
    .field_energy_writes = field_energy_writes,
    .integrated_mom_writes = integrated_mom_writes,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };

  return ctx;
}

void
evalFieldInit(double t, const double* GKYL_RESTRICT zc, double* GKYL_RESTRICT fout, void *ctx)
{
  double r = zc[0], theta = zc[1];
  struct annular_wg_ctx *app = ctx;

  double w_mode = app->w_mode;
  int bessel_order = app->bessel_order;
  double a_coeff = app->a_coeff;
  double b_coeff = app->b_coeff;

  double Ez_radial = (a_coeff * jn(bessel_order, r * w_mode)) + (b_coeff * yn(bessel_order, r * w_mode));

  double Ex = 0.0; // Total electric field (x-direction).
  double Ey = 0.0; // Total electric field (y-direction).
  double Ez = Ez_radial * cos(2.0 * theta); // Total electric field (z-direction).
  
  double Bx = 0.0; // Total magnetic field (x-direction).
  double By = 0.0; // Total magnetic field (y-direction).
  double Bz = 0.0; // Total magnetic field (z-direction).

  // Set electric field.
  fout[0] = Ex, fout[1] = Ey; fout[2] = Ez;
  // Set magnetic field.
  fout[3] = Bx, fout[4] = By; fout[5] = Bz;
  // Set correction potentials.
  fout[6] = 0.0; fout[7] = 0.0;
}

static inline void
mapc2p(double t, const double* GKYL_RESTRICT zc, double* GKYL_RESTRICT xp, void* ctx)
{
  double r = zc[0], theta = zc[1];

  // Set physical coordinates (x, y) from computational coordinates (r, theta).
  xp[0] = r * cos(theta);
  xp[1] = r * sin(theta);
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_moment_app* app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr)) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_moment_app_write(app, t_curr, frame);
  }
}

void
write_field_energy(struct gkyl_tm_trigger* fet, gkyl_moment_app* app, double t_curr)
{
  if (gkyl_tm_trigger_check_and_bump(fet, t_curr)) {
    gkyl_moment_app_calc_field_energy(app, t_curr);
    gkyl_moment_app_write_field_energy(app);
  }
}

void
write_integrated_mom(struct gkyl_tm_trigger* imt, gkyl_moment_app* app, double t_curr)
{
  if (gkyl_tm_trigger_check_and_bump(imt, t_curr)) {
    gkyl_moment_app_calc_integrated_mom(app, t_curr);
    gkyl_moment_app_write_integrated_mom(app);
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

  struct annular_wg_ctx ctx = create_ctx(); // Context for initialization functions.

  int NR = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nr);
  int NTHETA = APP_ARGS_CHOOSE(app_args.xcells[1], ctx.Ntheta);

  // Field.
  struct gkyl_moment_field field = {
    .epsilon0 = ctx.epsilon0, .mu0 = ctx.mu0,
    
    .limiter = GKYL_NO_LIMITER,
    .init = evalFieldInit,
    .ctx = &ctx,

    .bcx = { GKYL_FIELD_PEC_WALL, GKYL_FIELD_PEC_WALL },
  };

  int nrank = 1; // Number of processes in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif

  // Create global range.
  int cells[] = { NR, NTHETA };
  int dim = sizeof(cells) / sizeof(cells[0]);

  int cuts[dim];
#ifdef GKYL_HAVE_MPI
  for (int d = 0; d < dim; d++) {
    if (app_args.use_mpi) {
      cuts[d] = app_args.cuts[d];
    }
    else {
      cuts[d] = 1;
    }
  }
#else
  for (int d = 0; d < dim; d++) {
    cuts[d] = 1;
  }
#endif

  // Construct communicator for use in app.
  struct gkyl_comm *comm;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
      }
    );
  }
  else {
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .use_gpu = app_args.use_gpu
      }
    );
  }
#else
  comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .use_gpu = app_args.use_gpu
    }
  );
#endif

  int my_rank;
  gkyl_comm_get_rank(comm, &my_rank);
  int comm_size;
  gkyl_comm_get_size(comm, &comm_size);

  int ncuts = 1;
  for (int d = 0; d < dim; d++) {
    ncuts *= cuts[d];
  }

  if (ncuts != comm_size) {
    if (my_rank == 0) {
      fprintf(stderr, "*** Number of ranks, %d, does not match total cuts, %d!\n", comm_size, ncuts);
    }
    goto mpifinalize;
  }

  // Moment app.
  struct gkyl_moment app_inp = {

    .ndim = 2,
    .lower = { 2.0, 0.0 },
    .upper = { 2.0 + ctx.Lr, ctx.Ltheta },
    .cells = { NR, NTHETA },

    .mapc2p = mapc2p,

    .num_periodic_dir = 1,
    .periodic_dirs = { 1 },

    .cfl_frac = ctx.cfl_frac,

    .field = field,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0], app_args.cuts[1] },
      .comm = comm,
    },
  };

  // Create app object.
  // Set app output name from the executable name (argv[0]).
  snprintf(app_inp.name, sizeof(app_inp.name), "%s", app_args.app_name);
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  // Initial and final simulation times.
  double t_curr = 0.0, t_end = ctx.t_end;

  // Initialize simulation.
  int frame_curr = 0;
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_moment_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_moment_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n", gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_moment_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_moment_app_cout(app, stdout, " at time = %g\n", t_curr);
  }
  else {
    gkyl_moment_app_apply_ic(app, t_curr);
  }

  // Create trigger for IO.
  int num_frames = ctx.num_frames;
  struct gkyl_tm_trigger io_trig = { .dt = t_end / num_frames, .tcurr = t_curr, .curr = frame_curr };

  write_data(&io_trig, app, t_curr, false);

  // Create trigger for field energy.
  int field_energy_writes = ctx.field_energy_writes;
  struct gkyl_tm_trigger fe_trig = { .dt = t_end / field_energy_writes, .tcurr = t_curr, .curr = frame_curr };

  write_field_energy(&fe_trig, app, t_curr);

  // Create trigger for integrated moments.
  int integrated_mom_writes = ctx.integrated_mom_writes;
  struct gkyl_tm_trigger im_trig = { .dt = t_end / integrated_mom_writes, .tcurr = t_curr, .curr = frame_curr };

  write_integrated_mom(&im_trig, app, t_curr);

  // Compute initial guess of maximum stable time-step.
  double dt = t_end - t_curr;

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = ctx.dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx.num_failures_max;

  long step = 1;
  while ((t_curr < t_end) && (step <= app_args.num_steps)) {
    gkyl_moment_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = gkyl_moment_update(app, dt);
    gkyl_moment_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      gkyl_moment_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    write_data(&io_trig, app, t_curr, false);
    write_field_energy(&fe_trig, app, t_curr);
    write_integrated_mom(&im_trig, app, t_curr);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_moment_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_moment_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_moment_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_moment_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_moment_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);
        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }

  write_data(&io_trig, app, t_curr, false);
  write_field_energy(&fe_trig, app, t_curr);
  write_integrated_mom(&im_trig, app, t_curr);
  gkyl_moment_app_stat_write(app);

  struct gkyl_moment_stat stat = gkyl_moment_app_stat(app);

  gkyl_moment_app_cout(app, stdout, "\n");
  gkyl_moment_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_moment_app_cout(app, stdout, "Number of failed time-steps %ld\n", stat.nfail);
  gkyl_moment_app_cout(app, stdout, "Species updates took %g secs\n", stat.species_tm);
  gkyl_moment_app_cout(app, stdout, "Field updates took %g secs\n", stat.field_tm);
  gkyl_moment_app_cout(app, stdout, "Source updates took %g secs\n", stat.sources_tm);
  gkyl_moment_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

freeresources:
  // Free resources after simulation completion.
  gkyl_comm_release(comm);
  gkyl_moment_app_release(app);  
  
mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif
  
  return 0;
}