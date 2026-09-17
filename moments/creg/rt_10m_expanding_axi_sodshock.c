#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_ten_moment.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

#include <rt_arg_parse.h>
#include <gkyl_kann_net.h>
#include <gkyl_knutils.h>

struct expanding_axi_sodshock_ctx
{
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using normalized code units).
  double gas_gamma; // Adiabatic index.
  double U0; // (Initial) comoving plasma velocity.
  double R0; // (Initial) radial distance from expansion/contraction center.

  double rhol; // Left/inner fluid mass density.
  double ul; // Left/inner fluid velocity.
  double pl; // Left/inner fluid pressure.

  double rhor; // Right/outer fluid mass density.
  double ur; // Right/outer fluid velocity.
  double pr; // Right/outer fluid pressure.

  // Simulation parameters.
  int Nr; // Cell count (radial direction).
  int Ntheta; // Cell count (angular direction).
  double Lr; // Domain size (radial direction).
  double Ltheta; // Domain size (angular direction).
  double k0; // Closure parameter.
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  int field_energy_calcs; // Number of times to calculate field energy.
  int integrated_mom_calcs; // Number of times to calculate integrated moments.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.

  double rloc; // Fluid boundary (radial coordinate).

  // Neural network parameters.
  bool use_nn_closure; // Use neural network-based closure?
  int poly_order; // Polynomial order of learned DG coefficients.
  const char* nn_closure_file; // File path of neural network to use.
};

struct expanding_axi_sodshock_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  double gas_gamma = 5.0 / 3.0; // Adiabatic index.
  double U0 = 1.0; // (Initial) comoving plasma velocity.
  double R0 = 1.0; // (Initial) radial distance from expansion/contraction center.

  double rhol = 3.0; // Left/inner fluid mass density.
  double ul = 0.0; // Left/inner fluid velocity.
  double pl = 3.0; // Left/inner fluid pressure.

  double rhor = 1.0; // Right/outer fluid mass density.
  double ur = 0.0; // Right/outer fluid velocity.
  double pr = 1.0; // Right/outer fluid pressure.

  // Simulation parameters.
  int Nr = 128; // Cell count (radial direction).
  int Ntheta = 128 * 6; // Cell count (angular direction).
  double Lr = 1.0; // Domain size (radial direction).
  double Ltheta = 2.0 * pi; // Domain size (angular direction).
  double k0 = 0.1; // Closure parameter.
  double cfl_frac = 0.9; // CFL coefficient.

  double t_end = 0.2; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  double rloc = 0.5 * (0.25 + 1.25); // Fluid boundary (radial coordinate).

  // Neural network parameters.
  bool use_nn_closure = false; // Use neural network-based closure?
  int poly_order = 1; // Polynomial order of learned DG coefficients.
  const char* nn_closure_file = "moments/data/neural_nets/pkpm_periodic_es_shock_p1_moms_nn_1"; // File path of neural network to use.

  struct expanding_axi_sodshock_ctx ctx = {
    .pi = pi,
    .gas_gamma = gas_gamma,
    .U0 = U0,
    .R0 = R0,
    .rhol = rhol,
    .ul = ul,
    .pl = pl,
    .rhor = rhor,
    .ur = ur,
    .pr = pr,
    .Nr = Nr,
    .Ntheta = Ntheta,
    .Lr = Lr,
    .Ltheta = Ltheta,
    .k0 = k0,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .num_frames = num_frames,
    .field_energy_calcs = field_energy_calcs,
    .integrated_mom_calcs = integrated_mom_calcs,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
    .rloc = rloc,
    .use_nn_closure = use_nn_closure,
    .poly_order = poly_order,
    .nn_closure_file = nn_closure_file,
  };

  return ctx;
}

void
eval10mInit(double t, const double* GKYL_RESTRICT zc, double* GKYL_RESTRICT fout, void* ctx)
{
  double r = zc[0];
  struct expanding_axi_sodshock_ctx *app = ctx;

  double gas_gamma = app->gas_gamma;

  double rhol = app->rhol;
  double ul = app->ul;
  double pl = app->pl;

  double rhor = app->rhor;
  double ur = app->ur;
  double pr = app->pr;

  double rloc = app->rloc;

  double rho = 0.0;
  double u = 0.0;
  double p = 0.0;

  if (r < rloc) {
    rho = rhol; // Fluid mass density (left/inner).
    u = ul; // Fluid velocity (left/inner).
    p = pl; // Fluid pressure (left/inner).
  }
  else {
    rho = rhor; // Fluid mass density (right/outer).
    u = ur; // Fluid velocity (right/outer).
    p = pr; // Fluid pressure (right/outer).
  }

  double mom_x = rho * u; // Fluid momentum density (x-direction).
  double mom_y = 0.0; // Fluid momentum density (y-direction).
  double mom_z = 0.0; // Fluid momentum density (z-direction).

  double pr_xx = p + (0.5 * rho * u * u); // Fluid pressure tensor (xx-component).
  double pr_xy = 0.0; // Fluid pressure tensor (xy-component).
  double pr_xz = 0.0; // Fluid pressure tensor (xz-component).
  double pr_yy = p; // Fluid pressure tensor (yy-component).
  double pr_yz = 0.0; // Fluid pressure tensor (yz-component).
  double pr_zz = p; // Fluid pressure tensor (zz-component).
  
  // Set fluid mass density.
  fout[0] = rho;
  // Set fluid momentum density.
  fout[1] = mom_x; fout[2] = mom_y; fout[3] = mom_z;
  // Set fluid pressure tensor.
  fout[4] = pr_xx; fout[5] = pr_xy; fout[6] = pr_xz;
  fout[7] = pr_yy; fout[8] = pr_yz; fout[9] = pr_zz;
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

  struct expanding_axi_sodshock_ctx ctx = create_ctx(); // Context for initialization functions.

  int NR = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nr);
  int NTHETA = APP_ARGS_CHOOSE(app_args.xcells[1], ctx.Ntheta);

  struct gkyl_kann_net **ann = gkyl_malloc(sizeof(struct gkyl_kann_net*) * 1);
  if (ctx.use_nn_closure) {
    const char *fmt_10m = "%s-%s.dat";
    int sz_10m = gkyl_calc_strlen(fmt_10m, ctx.nn_closure_file, "10m");
    char fileNm_10m[sz_10m + 1];
    snprintf(fileNm_10m, sizeof fileNm_10m, fmt_10m, ctx.nn_closure_file, "10m");
    FILE *file_10m = fopen(fileNm_10m, "r");
    if (file_10m != NULL) {
      ann[0] = gkyl_kann_net_load(fileNm_10m, app_args.use_gpu);
      fclose(file_10m);
    }
    else {
      ann[0] = 0;
      ctx.use_nn_closure = false;
      fprintf(stderr, "Neural network for 10m species not found! Disabling NN-based closure.\n");
    }
  }

  // Fluid equations.
  struct gkyl_wv_eqn *ten_moment = gkyl_wv_ten_moment_new(ctx.k0, false, ctx.use_nn_closure, ctx.poly_order, ann[0], app_args.use_gpu);

  struct gkyl_moment_species fluid = {
    .name = "10m",
    .equation = ten_moment,
    
    .init = eval10mInit,
    .ctx = &ctx,

    .has_volume_sources = true,
    .volume_gas_gamma = ctx.gas_gamma,
    .volume_U0 = ctx.U0,
    .volume_R0 = ctx.R0,

    .bcx = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
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
    .lower = { 0.25, 0.0 },
    .upper = { 0.25 + ctx.Lr, ctx.Ltheta },
    .cells = { NR, NTHETA },

    .mapc2p = mapc2p,

    .num_periodic_dir = 1,
    .periodic_dirs = { 1 },

    .cfl_frac = ctx.cfl_frac,

    .num_species = 1,
    .species = { fluid },

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
  if (ctx.use_nn_closure) {
    for (int i = 0; i < app_inp.num_species; i++) {
      gkyl_kann_net_release(ann[i]);
    }
  }
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
  gkyl_wv_eqn_release(ten_moment);
  gkyl_comm_release(comm);
  gkyl_moment_app_release(app);
  gkyl_free(ann);
  
mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif
  
  return 0;
}
