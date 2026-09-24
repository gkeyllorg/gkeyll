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

struct elc_heat_flux_ctx
{
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using normalized code units).
  double gas_gamma; // Adiabatic index.
  double epsilon0; // Permittivity of free space.
  double mu0; // Permeability of free space.
  double mass_ion; // Ion mass.
  double charge_ion; // Ion charge.
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.

  double n0; // Reference number density.
  double coll_fac; // Collision factor.

  // Derived physical quantities (using normalized code units).
  double T_high; // High reference temperature (eV).
  double T_low; // Low reference temperature (eV).
  double vte; // Electron thermal velocity.

  double rho_elc; // Electron mass density.
  double rho_ion; // Ion mass density.

  double E_ion; // Ion total energy density.
  double E_elc_lower; // Electron lower boundary total energy density.
  double E_elc_upper; // Electron upper boundary total energy density.

  double tau; // Collision time.
  double lambda; // Collision wavelength.

  // Simulation parameters.
  int Nx; // Cell count (x-direction).
  double Lx; // Domain size (x-direction).
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct elc_heat_flux_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  double gas_gamma = 5.0 / 3.0; // Adiabatic index.
  double epsilon0 = 8.854187817620389850536563031710750260608e-12; // Permittivity of free space.
  double mu0 = 12.56637061435917295385057353311801153679e-7; // Permeability of free space.
  double mass_ion = 1.672621637e-27; // Ion mass.
  double charge_ion = 1.602176487e-19; // Ion charge.
  double mass_elc = 9.10938215e-31; // Electron mass.
  double charge_elc = -1.602176487e-19; // Electron charge.

  double n0 = 1.0e18; // Reference number density.
  double coll_fac = 1.0; // Collision factor.

  // Derived physical quantities (using normalized code units).
  double T_high = 100.0 * charge_ion; // High reference temperature (eV).
  double T_low = 1.0 * charge_ion; // Low reference temperature (eV).
  double vte = sqrt(T_low / mass_elc); // Electron thermal velocity.

  double rho_elc = n0 * mass_elc; // Electron mass density.
  double rho_ion = n0 * mass_ion; // Ion mass density.

  double E_ion = n0 * T_low / (gas_gamma - 1.0); // Ion total energy density.
  double E_elc_lower = n0 * T_high / (gas_gamma - 1.0); // Electron lower boundary total energy density.
  double E_elc_upper = n0 * T_low / (gas_gamma - 1.0); // Electron upper boundary total energy density.

  double tau = 6.0 * sqrt(2.0 * pi * mass_elc * T_low * pi * T_low * pi * T_low) * epsilon0 * epsilon0 /
    (charge_ion * charge_ion * charge_ion * charge_ion * n0); // Collision time.
  double lambda = vte * tau; // Collision wavelength.

  // Simulation parameters.
  int Nx = 64; // Cell count (x-direction).
  double Lx = 1000.0 * lambda; // Domain size (x-direction).
  double cfl_frac = 0.05; // CFL coefficient.

  double t_end = 1.0 * tau; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct elc_heat_flux_ctx ctx = {
    .pi = pi,
    .gas_gamma = gas_gamma,
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .mass_ion = mass_ion,
    .charge_ion = charge_ion,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .n0 = n0,
    .coll_fac = coll_fac,
    .T_high = T_high,
    .T_low = T_low,
    .vte = vte,
    .rho_elc = rho_elc,
    .rho_ion = rho_ion,
    .E_ion = E_ion,
    .E_elc_lower = E_elc_lower,
    .E_elc_upper = E_elc_upper,
    .tau = tau,
    .lambda = lambda,
    .Nx = Nx,
    .Lx = Lx,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .num_frames = num_frames,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };

  return ctx;
}

void
evalElcInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0];
  struct elc_heat_flux_ctx *app = ctx;

  double gas_gamma = app->gas_gamma;

  double T_high = app->T_high;
  double T_low = app->T_low;

  double n0 = app->n0;
  double rho_elc = app->rho_elc;

  double Lx = app->Lx;

  double T = T_high + (T_low - T_high) * x / Lx;
  double E_elc = n0 * T / (gas_gamma - 1.0);

  // Set electron mass density.
  fout[0] = rho_elc;
  // Set electron momentum density.
  fout[1] = 0.0; fout[2] = 0.0; fout[3] = 0.0;
  // Set electron total energy density.
  fout[4] = E_elc;
}

void
evalIonInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct elc_heat_flux_ctx *app = ctx;

  double rho_ion = app->rho_ion;
  double E_ion = app->E_ion;

  // Set ion mass density.
  fout[0] = rho_ion;
  // Set ion momentum density.
  fout[1] = 0.0; fout[2] = 0.0; fout[3] = 0.0;
  // Set ion total energy density.
  fout[4] = E_ion;
}

void
evalFieldInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Set electric field.
  fout[0] = 0.0, fout[1] = 0.0; fout[2] = 0.0;
  // Set magnetic field.
  fout[3] = 0.0, fout[4] = 0.0; fout[5] = 0.0;
  // Set correction potentials.
  fout[6] = 0.0; fout[7] = 0.0;
}

void
evalElcLowerBC(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  struct elc_heat_flux_ctx *app = ctx;

  double rho_elc = app->rho_elc;
  double E_elc_lower = app->E_elc_lower;

  // Set electron lower boundary mass density.
  ghost[0] = rho_elc;
  // Set electron lower boundary momentum density.
  ghost[1] = 0.0; ghost[2] = 0.0; ghost[3] = 0.0;
  // Set electron lower boundary total energy density.
  ghost[4] = E_elc_lower;
}

void
evalElcUpperBC(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  struct elc_heat_flux_ctx *app = ctx;

  double rho_elc = app->rho_elc;
  double E_elc_upper = app->E_elc_upper;

  // Set electron upper boundary mass density.
  ghost[0] = rho_elc;
  // Set electron upper boundary momentum density.
  ghost[1] = 0.0; ghost[2] = 0.0; ghost[3] = 0.0;
  // Set electron upper boundary total energy density.
  ghost[4] = E_elc_upper;
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

  struct elc_heat_flux_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);

  // Electron/ion equations.
  struct gkyl_wv_eqn *elc_euler = gkyl_wv_euler_new(ctx.gas_gamma, app_args.use_gpu);
  struct gkyl_wv_eqn *ion_euler = gkyl_wv_euler_new(ctx.gas_gamma, app_args.use_gpu);

  struct gkyl_moment_species elc = {
    .name = "elc",
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .equation = elc_euler,
    
    .init = evalElcInit,
    .ctx = &ctx,

    .type_brag = GKYL_BRAG_UNMAG_FULL,

    .bcx = { GKYL_SPECIES_FUNC, GKYL_SPECIES_FUNC },  
    .bcx_func = { evalElcLowerBC, evalElcUpperBC}, 
  };

  struct gkyl_moment_species ion = {
    .name = "ion",
    .charge = ctx.charge_ion, .mass = ctx.mass_ion,
    .equation = ion_euler,
    
    .init = evalIonInit,
    .ctx = &ctx,  

    .bcx = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
  }; 

  // Field.
  struct gkyl_moment_field field = {
    .epsilon0 = ctx.epsilon0, .mu0 = ctx.mu0,
    .mag_error_speed_fact = 1.0,
    
    .is_static = false,
    .init = evalFieldInit,
    .ctx = &ctx,
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
  struct gkyl_range global_r;
  gkyl_create_global_range(dim, cells, &global_r);

  // Create decomposition.
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

  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(dim, cuts, &global_r);

  // Construct communicator for use in app.
  struct gkyl_comm *comm;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
        .decomp = decomp
      }
    );
  }
  else {
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .decomp = decomp,
        .use_gpu = app_args.use_gpu
      }
    );
  }
#else
  comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .decomp = decomp,
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
    .lower = { 0.0 },
    .upper = { ctx.Lx },
    .cells = { NX },

    .num_periodic_dir = 0,
    .periodic_dirs = { },
    .cfl_frac = ctx.cfl_frac,

    .has_braginskii = true,
    .coll_fac = ctx.coll_fac,

    .num_species = 2,
    .species = { elc, ion },

    .field = field,

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

  // Create trigger for IO.
  int num_frames = ctx.num_frames;
  struct gkyl_tm_trigger io_trig = { .dt = t_end / num_frames };

  // Initialize simulation.
  gkyl_moment_app_apply_ic(app, t_curr);
  write_data(&io_trig, app, t_curr, false);

  // Compute estimate of maximum stable time-step.
  double dt = gkyl_moment_app_max_dt(app);

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
  gkyl_moment_app_stat_write(app);

  struct gkyl_moment_stat stat = gkyl_moment_app_stat(app);

  gkyl_moment_app_cout(app, stdout, "\n");
  gkyl_moment_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_moment_app_cout(app, stdout, "Number of failed time-steps %ld\n", stat.nfail);
  gkyl_moment_app_cout(app, stdout, "Species updates took %g secs\n", stat.species_tm);
  gkyl_moment_app_cout(app, stdout, "Field updates took %g secs\n", stat.field_tm);
  gkyl_moment_app_cout(app, stdout, "Source updates took %g secs\n", stat.sources_tm);
  gkyl_moment_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

  // Free resources after simulation completion.
  gkyl_wv_eqn_release(elc_euler);
  gkyl_wv_eqn_release(ion_euler);
  gkyl_rect_decomp_release(decomp);
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
