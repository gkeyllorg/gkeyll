// 2D Kelvin-Helmholtz instability test for the Euler equations using the DG/Vlasov solver.
// Input parameters match the initial conditions found in Section 4 of the article:
// M. Sementilli, R. Zangeneh and J. Chen (2024), "Influence of Cross Perturbations on Turbulent Kelvin-Helmholtz Instability",
// Fluids, Volume 9 (3): 52.
// https://www.mdpi.com/2311-5521/9/3/52

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>
#include <gkyl_util.h>
#include <gkyl_wv_euler.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#ifdef GKYL_HAVE_NCCL
#include <gkyl_nccl_comm.h>
#endif
#endif

#include <rt_arg_parse.h>

struct kh_2d_ctx
{
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using normalized code units).
  double gas_gamma; // Adiabatic idex.

  double rhol; // Left/inner fluid mass density.
  double ul; // Left/inner fluid velocity (x-direction).
  double pl; // Left/inner fluid pressure.

  double rhor; // Right/outer fluid mass density.
  double ur; // Right/outer fluid velocity (x-direction).
  double pr; // Right/outer fluid pressure.

  double yloc; // Fluid boundary (y-coordinate).

  // Simulation parameters.
  int Nx; // Cell count (configuration space: x-direction).
  int Ny; // Cell count (configuration space: y-direction).
  double Lx; // Domain size (configuration space: x-direction).
  double Ly; // Domain size (configuration space: y-direction).
  int poly_order; // Polynomial order.
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct kh_2d_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  double gas_gamma = 1.4; // Adiabatic index.

  double rhol = 2.0; // Left/inner fluid mass density.
  double ul = -0.5; // Left/inner fluid velocity (x-direction).
  double pl = 2.5; // Left/inner fluid pressure.

  double rhor = 1.0; // Right/outer fluid mass density.
  double ur = 0.5; // Right/outer fluid velocity (x-direction).
  double pr = 2.5; // Right/outer fluid pressure.

  double yloc = 0.25; // Fluid boundary (y-coordinate).

  // Simulation parameters.
  int Nx = 64; // Cell count (configuration space: x-direction).
  int Ny = 64; // Cell count (configuration space: y-direction).
  double Lx = 1.0; // Domain size (configuration space: x-direction).
  double Ly = 1.0; // Domain size (configuration space: y-direction).
  int poly_order = 1; // Polynomial order.
  double cfl_frac = 0.9; // CFL coefficient.

  double t_end = 10.0; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct kh_2d_ctx ctx = {
    .pi = pi,
    .gas_gamma = gas_gamma,
    .rhol = rhol,
    .ul = ul,
    .pl = pl,
    .rhor = rhor,
    .ur = ur,
    .pr = pr,
    .yloc = yloc,
    .Nx = Nx,
    .Ny = Ny,
    .Lx = Lx,
    .Ly = Ly,
    .poly_order = poly_order,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .num_frames = num_frames,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };

  return ctx;
}

void
evalEulerInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct kh_2d_ctx *app = ctx;

  double pi = app -> pi;

  double gas_gamma = app -> gas_gamma;
  
  double rhol = app -> rhol;
  double ul = app -> ul;
  double pl = app -> pl;

  double rhor = app -> rhor;
  double ur = app -> ur;
  double pr = app -> pr;

  double yloc = app -> yloc;

  double rho = 0.0;
  double vx = 0.0;
  double vy = 0.0;
  double p = 0.0;

  if (fabs(y) < yloc) {
    rho = rhol; // Fluid mass density (left/inner).
    vx = ul; // Fluid x-velocity (left/inner).
    p = pl; // Fluid pressure (left/inner).
  }
  else {
    rho = rhor; // Fluid mass density (right/outer).
    vx = ur; // Fluid x-velocity (right/outer).
    p = pr; // Fluid pressure (right/outer).
  }

  pcg64_random_t rng = gkyl_pcg64_init(0); // Random number generator;

  double alpha = 1.0e-2;
  double k = 2.0 * pi;

  for (int i = 0; i < 16; i++) {
    for (int j = 0; j < 16; j++) {
      vx += alpha * gkyl_pcg64_rand_double(&rng) * sin(i * k * x + j * k * y + 2.0 * pi * gkyl_pcg64_rand_double(&rng));
      vy += alpha * gkyl_pcg64_rand_double(&rng) * sin(i * k * x + j * k * y + 2.0 * pi * gkyl_pcg64_rand_double(&rng));
    }
  }

  // Set fluid mass density.
  fout[0] = rho;
  // Set fluid momentum density.
  fout[1] = rho * vx; fout[2] = rho * vy; fout[3] = 0.0;
  // Set fluid total energy density.
  fout[4] = p / (gas_gamma - 1.0) + 0.5 * rho * (vx * vx + vy * vy);
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_vlasov_app* app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr)) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_vlasov_app_write(app, t_curr, iot->curr - 1);
    gkyl_vlasov_app_write_mom(app, t_curr, iot->curr - 1);
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

  struct kh_2d_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], ctx.Ny);

  // Fluid equations.
  struct gkyl_wv_eqn *euler = gkyl_wv_euler_new(ctx.gas_gamma, app_args.use_gpu);

  struct gkyl_vlasov_fluid_species fluid = {
    .equation = euler,
    .init = evalEulerInit,
    .ctx = &ctx,
  };

  int nrank = 1; // Number of processors in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif  

  int ccells[] = { NX, NY };
  int cdim = sizeof(ccells) / sizeof(ccells[0]);

  int cuts[cdim];
#ifdef GKYL_HAVE_MPI  
  for (int d = 0; d < cdim; d++) {
    if (app_args.use_mpi) {
      cuts[d] = app_args.cuts[d];
    }
    else {
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
    comm = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
      }
    );
#else
    printf(" Using -g and -M together requires NCCL.\n");
    assert(0 == 1);
#endif
  }
  else if (app_args.use_mpi) {
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

    .cdim = 2, .vdim = 0,
    .lower = { -0.5 * ctx.Lx, -0.5 * ctx.Ly },
    .upper = { 0.5 * ctx.Lx, 0.5 * ctx.Ly }, 
    .cells = { NX, NY },

    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,
    .cfl_frac = ctx.cfl_frac,

    .num_periodic_dir = 2,
    .periodic_dirs = { 0, 1 },   

    .num_species = 1,
    .species = {
      { .name = "euler",
        .type = GKYL_SPECIES_FLUID, .fluid = fluid },
    },


    .skip_field = true,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0], app_args.cuts[1] },
      .comm = comm,
    },
  };

  // Create app object.
  // Set app output name from the executable name (argv[0]).
  snprintf(app_inp.name, sizeof(app_inp.name), "%s", app_args.app_name);
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&app_inp);

  // Initial and final simulation times.
  double t_curr = 0.0, t_end = ctx.t_end;

  // Create trigger for IO.
  int num_frames = ctx.num_frames;
  struct gkyl_tm_trigger io_trig = { .dt = t_end / num_frames };

  // Initialize simulation.
  gkyl_vlasov_app_apply_ic(app, t_curr);
  write_data(&io_trig, app, t_curr, false);

  // Compute initial guess of maximum stable time-step.
  double dt = t_end - t_curr;

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

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    write_data(&io_trig, app, t_curr, false);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_vlasov_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_vlasov_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_vlasov_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_vlasov_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_vlasov_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);
        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }

  write_data(&io_trig, app, t_curr, false);
  gkyl_vlasov_app_stat_write(app);

  struct gkyl_vlasov_stat stat = gkyl_vlasov_app_stat(app);

  gkyl_vlasov_app_cout(app, stdout, "\n");
  gkyl_vlasov_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_vlasov_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_vlasov_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_vlasov_app_cout(app, stdout, "  Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_vlasov_app_cout(app, stdout, "  Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  gkyl_vlasov_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_vlasov_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_vlasov_app_cout(app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm);
  gkyl_vlasov_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_vlasov_app_cout(app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm);
  gkyl_vlasov_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

  gkyl_vlasov_app_cout(app, stdout, "Number of write calls %ld\n", stat.n_io);
  double io_tm =  stat.field_io_tm + stat.species_io_tm + stat.field_diag_io_tm + stat.species_diag_io_tm;
  gkyl_vlasov_app_cout(app, stdout, "IO time took %g secs \n", io_tm);

  // Free resources after simulation completion.
  gkyl_wv_eqn_release(euler);
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
