// Propagation into a plasma wave beach test for the PKPM system of equations.
// Input parameters match the initial conditions found in entry JE8 of Ammar's Simulation Journal (https://ammar-hakim.org/sj/je/je8/je8-plasmabeach.html), adapted from Section III. A. of the article:
// D. N. Smithe (2007), "Finite-difference time-domain simulation of fusion plasmas at radiofrequency time scales",
// Physics of Plasmas, Volume 14 (5): 056104.
// https://pubs.aip.org/aip/pop/article/14/5/056104/929539/Finite-difference-time-domain-simulation-of-fusion
// NOTE: This simulation is a magnetized variant of JE8 and the Smithe Physics of Plasmas paper
// There is an initial Bx that causes the driven Jy to rotate between Jy and Jz. 

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

struct mom_beach_ctx
{
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using non-normalized physical units).
  double gas_gamma; // Adiabatic index.
  double epsilon0; // Permittivity of free space.
  double mu0; // Permeability of free space.
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.
  double vt_elc; // Electron thermal velocity.
  double nu_elc; // Electron-electron collision frequency. 
  double density_floor; // Electron density floor to avoid negative density in low-density region. 

  double J0; // Reference current density (Amps / m^3).
  
  // Derived physical quantities (using non-normalized physical units).
  double light_speed; // Speed of light.

  // Simulation parameters.
  int Nx; // Cell count (x-direction).
  double Lx; // Domain size (x-direction).
  double Lx100; // Domain size over 100 (x-direction).
  double x_last_edge; // Location of center of last cell.
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  double init_dt; // Initial time step guess so first step does not generate NaN
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.

  double deltaT; // Arbitrary constant, with units of time.
  double factor; // Numerical factor for calculation of electron number density.
  double omega_drive; // Drive current angular frequency.
};

struct mom_beach_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using non-normalized physical units).
  double gas_gamma = 5.0 / 3.0; // Adiabatic index.
  double epsilon0 = 8.854187817620389850536563031710750260608e-12; // Permittivity of free space.
  double mu0 = 12.56637061435917295385057353311801153679e-7; // Permeability of free space.
  double mass_elc = 9.10938215e-31; // Electron mass.
  double charge_elc = -1.602176487e-19; // Electron charge.
  double T_elc = -charge_elc; // Electron temperature. 
  double vt_elc = sqrt(T_elc/mass_elc); // Electron thermal velocity. 
  double density_floor = 1.0e1; // Electron density floor to avoid negative density in low-density region. 

  double J0 = 1.0e-12; // Reference current density (Amps / m^3).
  
  // Derived physical quantities (using non-normalized physical units).
  double light_speed = 1.0 / sqrt(mu0 * epsilon0); // Speed of light.

  // Simulation parameters.
  int Nx = 200; // Cell count (x-direction).
  double Lx = 1.0; // Domain size (x-direction).
  double Lx100 = Lx / 100.0; // Domain size over 100 (x-direction).
  double x_last_edge = Lx - Lx / Nx; // Location of center of last upper cell (low density side).
  double cfl_frac = 1.0; // CFL coefficient.
  double t_end = 5.0e-9; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  double deltaT = Lx100 / light_speed; // Arbitrary constant, with units of time.
  double factor = deltaT * deltaT * charge_elc * charge_elc / (mass_elc * epsilon0); // Numerical factor for calculation of electron number density.
  double omega_drive = pi / 10.0 / deltaT; // Drive current angular frequency.
  // initial dt guess so first step does not generate NaN
  double init_dt = ((Lx/Nx)/light_speed)/(5.0);

  // Coulomb logarithms.
  double n0 = 1.0/factor; // reference density
  double log_lambda_elc = 6.6 - 0.5 * log(n0 / 1.0e20) + 1.5 * log(T_elc / (-charge_elc));
  // Collision frequencies.
  double nu_elc = log_lambda_elc * pow(charge_elc,4) * n0 /
    (6.0 * sqrt(2.0) * pow(pi,3.0/2.0) * pow(epsilon0,2) * sqrt(mass_elc) * pow(T_elc,3.0/2.0));

  struct mom_beach_ctx ctx = {
    .pi = pi,
    .gas_gamma = gas_gamma,
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .density_floor = density_floor, 
    .J0 = J0,
    .vt_elc = vt_elc, 
    .nu_elc = nu_elc, 
    .light_speed = light_speed,
    .Nx = Nx,
    .Lx = Lx,
    .Lx100 = Lx100,
    .x_last_edge = x_last_edge,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .num_frames = num_frames,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
    .deltaT = deltaT,
    .factor = factor,
    .omega_drive = omega_drive,
    .init_dt = init_dt, 
  };

  return ctx;
}

static inline double
maxwellian(double n, double v, double vth)
{
  double v2 = v*v;
  return n/sqrt(2*M_PI*vth*vth)*exp(-v2/(2*vth*vth));
}

void
evalDistFuncElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], v = xn[1];
  struct mom_beach_ctx *app = ctx;

  double gas_gamma = app->gas_gamma;
  double epsilon0 = app->epsilon0;
  double mass_elc = app->mass_elc;
  double charge_elc = app->charge_elc;
  double density_floor = app->density_floor; 

  double light_speed = app->light_speed;

  double Lx100 = app->Lx100;

  double factor = app ->factor;

  double omegaPdt = 25.0 * (1.0 - x) * (1.0 - x) * (1.0 - x) * (1.0 - x) * (1.0 - x); // Plasma frequency profile.
  double ne = omegaPdt * omegaPdt / factor + density_floor; // Electron number density (with a floor for avoiding negative density).
  double vt_elc = app->vt_elc; // Electron thermal velocity. 

  fout[0] = maxwellian(ne, v, vt_elc);
  fout[1] = vt_elc*vt_elc*maxwellian(ne, v, vt_elc);
}

void
evalFluidElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0];
  // no initial flow
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
}

void
evalFieldFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Set electric field.
  fout[0] = 0.0, fout[1] = 0.0; fout[2] = 0.0;
  // Set magnetic field.
  fout[3] = 0.0, fout[4] = 0.0; fout[5] = 0.0;
  // Set correction potentials.
  fout[6] = 0.0; fout[7] = 0.0;
}

void
evalExtEmFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct mom_beach_ctx *app = ctx;
  double x = xn[0];
  double B_x = 1.0;
  
  fout[0] = 0.0; fout[1] = 0.0, fout[2] = 0.0;
  fout[3] = B_x; fout[4] = 0.0; fout[5] = 0.0;
}

void
evalAppCurrent(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0];
  struct mom_beach_ctx *app = ctx;

  double pi = app->pi;

  double J0 = app->J0;

  double light_speed = app->light_speed;

  double Lx = app->Lx;
  double Nx = app->Nx;
  double Lx100 = app->Lx100;
  double x_last_edge = app->x_last_edge;

  double omega_drive = app->omega_drive;
  fout[0] = 0.0;
  fout[2] = 0.0;
  if (x > x_last_edge) {
    // Set applied current.
    fout[1] = -J0 * sin(omega_drive * t);
  }
}

void
evalNuElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct mom_beach_ctx *app = ctx;
  fout[0] = app->nu_elc;
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_pkpm_app* app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr)) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_pkpm_app_write(app, t_curr, iot->curr-1);
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

  struct mom_beach_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int NV = APP_ARGS_CHOOSE(app_args.vcells[0], 32);

  // electrons
  struct gkyl_pkpm_species elc = {
    .name = "elc",
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .lower = { -6.0 * ctx.vt_elc},
    .upper = { 6.0 * ctx.vt_elc}, 
    .cells = { NV },

    .ctx_dist = &ctx,
    .ctx_fluid = &ctx,
    .init_dist = evalDistFuncElc,
    .init_fluid = evalFluidElc,

    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,

      .ctx = &ctx,
      .self_nu = evalNuElc,
    }, 
    .bcx = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
  };

  // field
  struct gkyl_pkpm_field field = {
    .epsilon0 = ctx.epsilon0, .mu0 = ctx.mu0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,

    .ctx = &ctx,
    .init = evalFieldFunc,

    .ext_em = evalExtEmFunc,
    .ext_em_ctx = &ctx,
    .app_current = evalAppCurrent, 
    .app_current_ctx = &ctx,     
    .app_current_evolve = true,
    .bcx = { GKYL_FIELD_COPY, GKYL_FIELD_COPY }, 
  };

  int nrank = 1; // Number of processes in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
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

  int ccells[] = { NX };
  int cdim = sizeof(ccells) / sizeof(ccells[0]);
  int ncuts = 1;
  for (int d = 0; d < cdim; d++) {
    ncuts *= app_args.cuts[d];
  }

  if (ncuts != comm_size) {
    if (my_rank == 0) {
      fprintf(stderr, "*** Number of ranks, %d, does not match total cuts, %d!\n", comm_size, ncuts);
    }
    goto mpifinalize;
  }

  // PKPM app
  struct gkyl_pkpm app_inp = {

    .cdim = 1, .vdim = 1,
    .lower = { 0.0 },
    .upper = { ctx.Lx }, 
    .cells = { NX },
    .poly_order = 2,
    .basis_type = app_args.basis_type,
    .cfl_frac = ctx.cfl_frac,

    .num_periodic_dir = 0,
    .periodic_dirs = { },

    .num_species = 1,
    .species = { elc },
    .field = field,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0] },
      .comm = comm,
    },
  };

  // create app object
  // Set app output name from the executable name (argv[0]).
  snprintf(app_inp.name, sizeof(app_inp.name), "%s", app_args.app_name);
  gkyl_pkpm_app *app = gkyl_pkpm_app_new(&app_inp);

  // start, end and initial time-step
  double t_curr = 0.0, t_end = ctx.t_end;
  double dt = ctx.init_dt;

  // Create trigger for IO.
  int num_frames = ctx.num_frames;
  struct gkyl_tm_trigger io_trig = { .dt = t_end / num_frames };

  // initialize simulation
  gkyl_pkpm_app_apply_ic(app, t_curr);
  write_data(&io_trig, app, t_curr, false);  

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

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    write_data(&io_trig, app, t_curr, false);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_pkpm_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_pkpm_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_pkpm_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_pkpm_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_pkpm_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);
        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }
  write_data(&io_trig, app, t_curr, false);
  gkyl_pkpm_app_stat_write(app);

  // fetch simulation statistics
  struct gkyl_pkpm_stat stat = gkyl_pkpm_app_stat(app);

  gkyl_pkpm_app_cout(app, stdout, "\n");
  gkyl_pkpm_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_pkpm_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_pkpm_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_pkpm_app_cout(app, stdout, "Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_pkpm_app_cout(app, stdout, "Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  gkyl_pkpm_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_pkpm_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_pkpm_app_cout(app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm);
  gkyl_pkpm_app_cout(app, stdout, "Fluid Species RHS calc took %g secs\n", stat.fluid_species_rhs_tm);
  gkyl_pkpm_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_pkpm_app_cout(app, stdout, "Species PKPM Vars took %g secs\n", stat.species_pkpm_vars_tm);
  gkyl_pkpm_app_cout(app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm);
  gkyl_pkpm_app_cout(app, stdout, "EM Variables (bvar) calculation took %g secs\n", stat.field_em_vars_tm);
  gkyl_pkpm_app_cout(app, stdout, "Current evaluation and accumulate took %g secs\n", stat.current_tm);

  gkyl_pkpm_app_cout(app, stdout, "Species BCs took %g secs\n", stat.species_bc_tm);
  gkyl_pkpm_app_cout(app, stdout, "Fluid Species BCs took %g secs\n", stat.fluid_species_bc_tm);
  gkyl_pkpm_app_cout(app, stdout, "Field BCs took %g secs\n", stat.field_bc_tm);
  
  gkyl_pkpm_app_cout(app, stdout, "Updates took %g secs\n", stat.total_tm);
  
  gkyl_pkpm_app_cout(app, stdout, "Number of write calls %ld,\n", stat.n_io);
  gkyl_pkpm_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

  gkyl_comm_release(comm);

  // simulation complete, free app
  gkyl_pkpm_app_release(app);

  mpifinalize:
  ;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif  
  
  return 0;
}
