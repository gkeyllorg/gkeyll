// Geospace Environmental Modeling (GEM) magnetic reconnection test for the 5-moment equations.
// Input parameters match the equilibrium and initial conditions in Section 2, from the article:
// J. Birn et al. (2001), "Geospace Environmental Modeling (GEM) Magnetic Reconnection Challenge",
// Journal of Geophysical Research: Space Physics, Volume 106 (A3): 3715-3719.
// https://agupubs.onlinelibrary.wiley.com/doi/10.1029/1999JA900449

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

struct gem_ctx
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

  double di; // Ion skin depth.
  double n0; // Reference in-plane number density.
  double nb_over_n0; // Background number density / reference number density.
  double B0; // Reference magnetic field strength.
  double Ti_over_Te; // Ion temperature / electron temperature.
  double beta; // Total plasma beta.
  
  // Derived physical quantities (using normalized code units).
  double lambda; // Current sheet width.
  double psi0; // Reference magnetic scalar potential.

  double Ti_frac; // Fraction of total temperature from ions.
  double Te_frac; // Fraction of total temperature from electrons.
  double T_tot; // Total temperature.

  double omega_ci; // Ion cyclotron frequency.

  // Simulation parameters.
  int Nx; // Cell count (x-direction).
  int Ny; // Cell count (y-direction).
  double Lx; // Domain size (x-direction).
  double Ly; // Domain size (y-direction).
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  int field_energy_calcs; // Number of times to calculate field energy.
  int integrated_mom_calcs; // Number of times to calculate integrated moments.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct gem_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  double gas_gamma = 5.0 / 3.0; // Adiabatic index.
  double epsilon0 = 1.0; // Permittivity of free space.
  double mu0 = 1.0; // Permeability of free space.
  double mass_ion = 1.0; // Ion mass.
  double charge_ion = 1.0; // Ion charge.
  double mass_elc = 1.0 / 25.0; // Electron mass.
  double charge_elc = -1.0; // Electron charge.

  double di = 1.0; // Ion skin depth.
  double n0 = 1.0; // Reference in-plane number density.
  double nb_over_n0 = 0.2; // Background number density / reference number density.
  double B0 = 0.1; // Reference magnetic field strength.
  double Ti_over_Te = 5.0; // Ion temperature / electron temperature.
  double beta = 1.0; // Total plasma beta.
  
  // Derived physical quantities (using normalized code units).
  double lambda = 0.5 * di; // Current sheet width.
  double psi0 = 0.1 * B0 * di; // Reference magnetic scalar potential.

  double Ti_frac = Ti_over_Te / (1.0 + Ti_over_Te); // Fraction of total temperature from ions.
  double Te_frac = 1.0 / (1.0 + Ti_over_Te); // Fraction of total temperature from electrons.
  double T_tot = beta * (B0 * B0) / 2.0 / n0; // Total temperature.

  double omega_ci = fabs(charge_ion * B0 / mass_ion); // Ion cyclotron frequency.

  // Simulation parameters.
  int Nx = 128; // Cell count (x-direction).
  int Ny = 64; // Cell count (y-direction).
  double Lx = 25.6 * di; // Domain size (x-direction).
  double Ly = 12.8 * di; // Domain size (y-direction).
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 25.0 / omega_ci; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gem_ctx ctx = {
    .pi = pi,
    .gas_gamma = gas_gamma,
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .mass_ion = mass_ion,
    .charge_ion = charge_ion,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .di = di,
    .n0 = n0,
    .nb_over_n0 = nb_over_n0,
    .B0 = B0,
    .Ti_over_Te = Ti_over_Te,
    .beta = beta,
    .lambda = lambda,
    .psi0 = psi0,
    .Ti_frac = Ti_frac,
    .Te_frac = Te_frac,
    .T_tot = T_tot,
    .omega_ci = omega_ci,
    .Nx = Nx,
    .Ny = Ny,
    .Lx = Lx,
    .Ly = Ly,
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
evalElcInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct gem_ctx *app = ctx;

  double gas_gamma = app->gas_gamma;
  double mass_elc = app->mass_elc;
  double charge_elc = app->charge_elc;
  double lambda = app->lambda;
  double n0 = app->n0;
  double nb_over_n0 = app->nb_over_n0;
  double B0 = app->B0;
  double beta = app->beta;

  double Te_frac = app->Te_frac;
  double T_tot = app->T_tot;

  double sech_sq = (1.0 / cosh(y / lambda)) * (1.0 / cosh(y / lambda)); // Hyperbolic secant squared.

  double n = n0 * (sech_sq + nb_over_n0); // Total number density.
  double Jz = -(B0 / lambda) * sech_sq; // Total current density (z-direction).

  double rhoe = n * mass_elc; // Electron mass density.
  double mome_x = 0.0; // Electron momentum density (x-direction).
  double mome_y = 0.0; // Electron momentum density (y-direction).
  double mome_z = (mass_elc / charge_elc) * Jz * Te_frac; // Electron momentum density (z-direction).
  double Ee_tot = n * T_tot * Te_frac / (gas_gamma - 1.0) + 0.5 * mome_z * mome_z / rhoe; // Electron total energy density.

  // Set electron mass density.
  fout[0] = rhoe;
  // Set electron momentum density.
  fout[1] = mome_x; fout[2] = mome_y; fout[3] = mome_z;
  // Set electron total energy density.
  fout[4] = Ee_tot;
}

void
evalIonInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct gem_ctx *app = ctx;

  double gas_gamma = app->gas_gamma;
  double mass_ion = app->mass_ion;
  double charge_ion = app->charge_ion;
  double lambda = app->lambda;
  double n0 = app->n0;
  double nb_over_n0 = app->nb_over_n0;
  double B0 = app->B0;
  double beta = app->beta;

  double Ti_frac = app->Ti_frac;
  double T_tot = app->T_tot;

  double sech_sq = (1.0 / cosh(y / lambda)) * (1.0 / cosh(y / lambda)); // Hyperbolic secant squared.

  double n = n0 * (sech_sq + nb_over_n0); // Total number density.
  double Jz = -(B0 / lambda) * sech_sq; // Total current density (z-direction).

  double rhoi = n * mass_ion; // Ion mass density.
  double momi_x = 0.0; // Ion momentum density (x-direction).
  double momi_y = 0.0; // Ion momentum density (y-direction).
  double momi_z = (mass_ion / charge_ion) * Jz * Ti_frac; // Ion momentum density (z-direction).
  double Ei_tot = n * T_tot * Ti_frac / (gas_gamma - 1.0) + 0.5 * momi_z * momi_z / rhoi; // Ion total energy density.

  // Set ion mass density.
  fout[0] = rhoi;
  // Set ion momentum density.
  fout[1] = momi_x; fout[2] = momi_y; fout[3] = momi_z;
  // Set ion total energy density.
  fout[4] = Ei_tot;
}

void
evalFieldInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct gem_ctx *app = ctx;

  double pi = app->pi;

  double lambda = app->lambda;
  double B0 = app->B0;

  double psi0 = app->psi0;

  double Lx = app->Lx;
  double Ly = app->Ly;

  double Bxb = B0 * tanh(y / lambda); // Total magnetic field strength.

  double Ex = 0.0; // Total electric field (x-direction).
  double Ey = 0.0; // Total electric field (y-direction).
  double Ez = 0.0; // Total electric field (z-direction).

  double Bx = Bxb - psi0 * (pi / Ly) * cos(2.0 * pi * x / Lx) * sin(pi * y / Ly); // Total magnetic field (x-direction).
  double By = psi0 * (2.0 * pi / Lx) * sin(2.0 * pi * x / Lx) * cos(pi * y / Ly); // Total magnetic field (y-direction).
  double Bz = 0.0; // Total magnetic field (z-direction).

  // Set electric field.
  fout[0] = Ex, fout[1] = Ey; fout[2] = Ez;
  // Set magnetic field.
  fout[3] = Bx, fout[4] = By; fout[5] = Bz;
  // Set correction potentials.
  fout[6] = 0.0; fout[7] = 0.0;
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

  struct gem_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], ctx.Ny);

  // Electron/ion equations.
  struct gkyl_wv_eqn *elc_euler = gkyl_wv_euler_new(ctx.gas_gamma, app_args.use_gpu);
  struct gkyl_wv_eqn *ion_euler = gkyl_wv_euler_new(ctx.gas_gamma, app_args.use_gpu);

  struct gkyl_moment_species elc = {
    .name = "elc",
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .equation = elc_euler,
    
    .init = evalElcInit,
    .ctx = &ctx,

    .bcy = { GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT },
  };

  struct gkyl_moment_species ion = {
    .name = "ion",
    .charge = ctx.charge_ion, .mass = ctx.mass_ion,
    .equation = ion_euler,
    
    .init = evalIonInit,
    .ctx = &ctx,

    .bcy = { GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT },    
  };

  // Field.
  struct gkyl_moment_field field = {
    .epsilon0 = ctx.epsilon0, .mu0 = ctx.mu0,
    .mag_error_speed_fact = 1.0,
    
    .init = evalFieldInit,
    .ctx = &ctx,
    
    .bcy = { GKYL_FIELD_PEC_WALL, GKYL_FIELD_PEC_WALL },
  };

  int nrank = 1; // Number of processes in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif

  int cells[] = { NX, NY };
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
    .lower = { -0.5 * ctx.Lx, -0.5 * ctx.Ly },
    .upper = { 0.5 * ctx.Lx, 0.5 * ctx.Ly }, 
    .cells = { NX, NY },

    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },
    .cfl_frac = ctx.cfl_frac,

    .num_species = 2,
    .species = { elc, ion },

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
  struct gkyl_tm_trigger io_trig = { .dt = t_end / num_frames, .tcurr = t_curr, .curr = frame_curr };

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
  gkyl_wv_eqn_release(elc_euler);
  gkyl_wv_eqn_release(ion_euler);
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
