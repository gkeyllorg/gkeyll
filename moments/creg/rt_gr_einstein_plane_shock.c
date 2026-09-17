// Shock tube test on a plane-symmetric Gowdy spacetime for the coupled fluid-Einstein equations, assuming an ultra-relativistic equation of state.
// Input parameters taken from the initial conditions in Section 9 (Riemann problem 1), from the article:
// A. P. Barnes, P. G. Lefloch, B. G. Schmidt and J. M. Stewart (2004), "The Glimm scheme for perfect fluids on plane-symmetric Gowdy spacetimes",
// Classical and Quantum Gravity, Volume 21 (22): 5043.
// https://iopscience.iop.org/article/10.1088/0264-9381/21/22/003

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_gr_medium.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

#include <rt_arg_parse.h>

struct gr_einstein_plane_shock_ctx
{
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using normalized code units).
  double gas_gamma; // Adiabatic index.
  double kappa; // Stress-energy prefactor in the Einstein field equations.

  double exp_2a; // Exponential appearing in dt and dx metric terms.

  // Derived physical quantities (using normalized code units).
  double rhol; // Left fluid mass density.
  double rhor; // Right fluid mass density.

  double Etot_l; // Left fluid total energy density.
  double Etot_r; // Right fluid total energy density.

  // Simulation parameters.
  int Nx; // Cell count (x-direction).
  double Lx; // Domain size (x-direction).
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  int field_energy_calcs; // Number of times to calculate field energy.
  int integrated_mom_calcs; // Number of times to calculate integrated moments.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct gr_einstein_plane_shock_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  double gas_gamma = 4.0 / 3.0; // Adiabatic index.
  double kappa = 8.0 * pi; // Stress-energy prefactor in the Einstein field equations.

  double exp_2a = exp(-4.0); // Exponential appearing in dt and dx metric terms.

  // Derived physical quantities (using normalized code units).
  double rhol = 100.0 / kappa; // Left fluid mass density.
  double rhor = 1.0 / kappa; // Right fluid mass density.

  double Etot_l = rhol; // Left fluid total energy density.
  double Etot_r = rhor; // Right fluid total energy density.

  // Simulation parameters.
  int Nx = 4096; // Cell count (x-direction).
  double Lx = 2.0; // Domain size (x-direction).
  double cfl_frac = 0.95; // CFL coefficient.

  double t_end = 0.5; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gr_einstein_plane_shock_ctx ctx = {
    .pi = pi,
    .gas_gamma = gas_gamma,
    .kappa = kappa,
    .exp_2a = exp_2a,
    .rhol = rhol,
    .rhor = rhor,
    .Etot_l = Etot_l,
    .Etot_r = Etot_r,
    .Nx = Nx,
    .Lx = Lx,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .num_frames = num_frames,
    .field_energy_calcs = field_energy_calcs,
    .integrated_mom_calcs = integrated_mom_calcs,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };

  return ctx;
}

void
evalGRMediumInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0];
  struct gr_einstein_plane_shock_ctx *app = ctx;

  double gas_gamma = app->gas_gamma;
  double kappa = app->kappa;

  double exp_2a = app->exp_2a;

  double Etot_l = app->Etot_l;
  double Etot_r = app->Etot_r;

  double Etot = 0.0;

  if (x < 0.0) {
    Etot = Etot_l; // Fluid total energy density (left).
  }
  else {
    Etot = Etot_r; // Fluid total energy density (right).
  }

  double a_dt = 0.0; // Time derivative of metric term a.
  double a_dx = 0.0; // Space derivative of metric term a.
  double b_dt = 0.0; // Time derivative of metric term b.
  double b_dx = -sqrt((kappa * exp_2a * Etot) / 3.0) * tan((0.5 * x * sqrt(3.0 * kappa * exp_2a * Etot))); // Space derivative of metric term b.
  double c_dt = 0.0; // Time derivative of metric term c.
  double c_dx = 0.0; // Space derivative of metric term c.

  double b_dx_plus = -sqrt((kappa * exp_2a * Etot) / 3.0) * tan((0.5 * (x + (0.5 * pow(10.0, -8.0))) * sqrt(3.0 * kappa * exp_2a * Etot)));
  double b_dx_minus = -sqrt((kappa * exp_2a * Etot) / 3.0) * tan((0.5 * (x - (0.5 * pow(10.0, -8.0))) * sqrt(3.0 * kappa * exp_2a * Etot)));

  double a_dt_dx = 0.0; // Mixed space-time derivative of metric term a.
  double a_dx_dx = 0.0; // Second space derivative of metric term a.
  double b_dt_dx = 0.0; // Mixed space-time derivative of metric term b.
  double b_dx_dx = (b_dx_plus - b_dx_minus) / pow(10.0, -8.0); // Second space derivative of metric term b.
  double c_dt_dx = 0.0; // Mixed space-time derivative of metric term c.
  double c_dx_dx = 0.0; // Second space derivative of metric term c.

  double mom_x = 0.0; // Fluid momentum (x-direction).

  // Set exponential appearing in dt and dx metric terms.
  fout[0] = exp_2a;
  // Set first time and space derivatives of metric terms.
  fout[1] = a_dt; fout[2] = a_dx;
  fout[3] = b_dt; fout[4] = b_dx;
  fout[5] = c_dt; fout[6] = c_dx;
  // Set second time and space derivatives of metric terms.
  fout[7] = a_dt_dx; fout[8] = a_dx_dx;
  fout[9] = b_dt_dx; fout[10] = b_dx_dx;
  fout[11] = c_dt_dx; fout[12] = c_dx_dx;
  // Set fluid total energy density.
  fout[13] = Etot;
  // Set fluid momentum density.
  fout[14] = mom_x;
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_moment_app* app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_write) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_moment_app_write(app, t_curr, frame);
    gkyl_moment_app_write_field_energy(app);
    gkyl_moment_app_write_integrated_mom(app);
  }
}

void
calc_field_energy(struct gkyl_tm_trigger* fet, gkyl_moment_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(fet, t_curr) || force_calc) {
    gkyl_moment_app_calc_field_energy(app, t_curr);
  }
}

void
calc_integrated_mom(struct gkyl_tm_trigger* imt, gkyl_moment_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(imt, t_curr) || force_calc) {
    gkyl_moment_app_calc_integrated_mom(app, t_curr);
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

  struct gr_einstein_plane_shock_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);

  // Fluid equations.
  struct gkyl_wv_eqn *gr_medium = gkyl_wv_gr_medium_new(ctx.gas_gamma, ctx.kappa, app_args.use_gpu);

  struct gkyl_moment_species fluid = {
    .name = "gr_medium",
    .equation = gr_medium,
    
    .init = evalGRMediumInit,
    .ctx = &ctx,

    .has_einstein_medium = true,
    .medium_gas_gamma = ctx.gas_gamma,
    .medium_kappa = ctx.kappa,

    .bcx = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
  };

  int nrank = 1; // Number of processes in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif

  // Create global range.
  int cells[] = { NX };
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

    .ndim = 1,
    .lower = { -0.5 * ctx.Lx },
    .upper = { 0.5 * ctx.Lx }, 
    .cells = { NX },

    .cfl_frac = ctx.cfl_frac,

    .num_species = 1,
    .species = { fluid },

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0] },
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

  // Create trigger for field energy.
  int field_energy_calcs = ctx.field_energy_calcs;
  struct gkyl_tm_trigger fe_trig = { .dt = t_end / field_energy_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_field_energy(&fe_trig, app, t_curr, false);

  // Create trigger for integrated moments.
  int integrated_mom_calcs = ctx.integrated_mom_calcs;
  struct gkyl_tm_trigger im_trig = { .dt = t_end / integrated_mom_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_integrated_mom(&im_trig, app, t_curr, false);

  // Create trigger for IO.
  int num_frames = ctx.num_frames;
  struct gkyl_tm_trigger io_trig = { .dt = t_end / num_frames, .tcurr = frame_curr * (t_end / num_frames), .curr = frame_curr };

  write_data(&io_trig, app, t_curr, false);

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

    calc_field_energy(&fe_trig, app, t_curr, false);
    calc_integrated_mom(&im_trig, app, t_curr, false);
    write_data(&io_trig, app, t_curr, false);

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

        calc_field_energy(&fe_trig, app, t_curr, true);
        calc_integrated_mom(&im_trig, app, t_curr, true);
        write_data(&io_trig, app, t_curr, true);

        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }

  calc_field_energy(&fe_trig, app, t_curr, false);
  calc_integrated_mom(&im_trig, app, t_curr, false);
  write_data(&io_trig, app, t_curr, false);
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
  gkyl_wv_eqn_release(gr_medium);
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
