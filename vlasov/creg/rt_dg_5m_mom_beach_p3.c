// Propagation into a plasma wave beach test for the 5-moment equations using the DG/Vlasov solver, with third-order polynomial reconstruction.
// Input parameters match the initial conditions found in entry JE8 of Ammar's Simulation Journal (https://ammar-hakim.org/sj/je/je8/je8-plasmabeach.html), adapted from Section III. A. of the article:
// D. N. Smithe (2007), "Finite-difference time-domain simulation of fusion plasmas at radiofrequency time scales",
// Physics of Plasmas, Volume 14 (5): 056104.
// https://pubs.aip.org/aip/pop/article/14/5/056104/929539/Finite-difference-time-domain-simulation-of-fusion

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

  double J0; // Reference current density (Amps / m^3).
  
  // Derived physical quantities (using non-normalized physical units).
  double light_speed; // Speed of light.

  // Simulation parameters.
  int Nx; // Cell count (x-direction).
  double Lx; // Domain size (x-direction).
  int poly_order; // Polynomial order.
  double Lx100; // Domain size over 100 (x-direction).
  double x_last_edge; // Location of center of last cell.
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  int field_energy_calcs; // Number of times to calculate field energy.
  int integrated_mom_calcs; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs; // Number of times to calculate integrated L2 norm of distribution function.
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

  double J0 = 1.0e-12; // Reference current density (Amps / m^3).
  
  // Derived physical quantities (using non-normalized physical units).
  double light_speed = 1.0 / sqrt(mu0 * epsilon0); // Speed of light.

  // Simulation parameters.
  int Nx = 200; // Cell count (x-direction).
  double Lx = 1.0; // Domain size (x-direction).
  int poly_order = 3; // Polynomial order.
  double Lx100 = Lx / 100.0; // Domain size over 100 (x-direction).
  double x_last_edge = Lx - Lx / Nx; // Location of center of last upper cell (low density side).
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 5.0e-9; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs = INT_MAX; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  double deltaT = Lx100 / light_speed; // Arbitrary constant, with units of time.
  double factor = deltaT * deltaT * charge_elc * charge_elc / (mass_elc * epsilon0); // Numerical factor for calculation of electron number density.
  double omega_drive = pi / 10.0 / deltaT; // Drive current angular frequency.

  struct mom_beach_ctx ctx = {
    .pi = pi,
    .gas_gamma = gas_gamma,
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .J0 = J0,
    .light_speed = light_speed,
    .Nx = Nx,
    .Lx = Lx,
    .poly_order = poly_order,
    .Lx100 = Lx100,
    .x_last_edge = x_last_edge,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .num_frames = num_frames,
    .field_energy_calcs = field_energy_calcs,
    .integrated_mom_calcs = integrated_mom_calcs,
    .integrated_L2_f_calcs = integrated_L2_f_calcs,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
    .deltaT = deltaT,
    .factor = factor,
    .omega_drive = omega_drive,
  };

  return ctx;
}

void
evalElcInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0];
  struct mom_beach_ctx *app = ctx;

  double gas_gamma = app->gas_gamma;
  double mass_elc = app->mass_elc;
  double charge_elc = app->charge_elc;

  double factor = app ->factor;

  double omegaPdt = 25.0 * (1.0 - x) * (1.0 - x) * (1.0 - x) * (1.0 - x) * (1.0 - x); // Plasma frequency profile.
  double ne = omegaPdt * omegaPdt / factor + 1.0; // Electron number density (with density floor).

  double rhoe = mass_elc * ne; // Electron mass density.
  double mome_x = 0.0; // Electron momentum density (x-direction).
  double mome_y = 0.0; // Electron momentum density (y-direction).
  double mome_z = 0.0; // Electron momentum density (z-direction).
  double Ee_tot = ne * (-charge_elc) / (gas_gamma - 1.0); // Electron total energy density.

  // Set electron mass density.
  fout[0] = rhoe;
  // Set electron momentum density.
  fout[1] = mome_x; fout[2] = mome_y; fout[3] = mome_z;
  // Set electron total energy density.
  fout[4] = Ee_tot;
}

void
evalFieldInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double Ex = 0.0; // Total electric field (x-direction).
  double Ey = 0.0; // Total electric field (y-direction).
  double Ez = 0.0; // Total electric field (z-direction).

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

void
evalAppCurrent(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0];
  struct mom_beach_ctx *app = ctx;

  double J0 = app->J0;
  double x_last_edge = app->x_last_edge;
  double omega_drive = app->omega_drive;

  double app_x = 0.0; // Applied current (x-direction).
  double app_y = 0.0;
  double app_z = 0.0; // Applied current (z-direction).

  if (x > x_last_edge) {
    app_y = -J0 * sin(omega_drive * t); // Applied current (y-direction, right).
  }
  else {
    app_y = 0.0; // Applied current (y-direction, left).
  }
  
  // Set applied current.
  fout[0] = app_x; fout[1] = app_y; fout[2] = app_z;
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_vlasov_app* app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr)) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_vlasov_app_write(app, t_curr, frame);
    gkyl_vlasov_app_write_field_energy(app);
    gkyl_vlasov_app_write_integrated_mom(app);
    gkyl_vlasov_app_write_integrated_L2_f(app);

    gkyl_vlasov_app_calc_mom(app);
    gkyl_vlasov_app_write_mom(app, t_curr, frame);
  }
}

void
calc_field_energy(struct gkyl_tm_trigger* fet, gkyl_vlasov_app* app, double t_curr)
{
  if (gkyl_tm_trigger_check_and_bump(fet, t_curr)) {
    gkyl_vlasov_app_calc_field_energy(app, t_curr);
  }
}

void
calc_integrated_mom(struct gkyl_tm_trigger* imt, gkyl_vlasov_app* app, double t_curr)
{
  if (gkyl_tm_trigger_check_and_bump(imt, t_curr)) {
    gkyl_vlasov_app_calc_integrated_mom(app, t_curr);
  }
}

void
calc_integrated_L2_f(struct gkyl_tm_trigger* l2t, gkyl_vlasov_app* app, double t_curr)
{
  if (gkyl_tm_trigger_check_and_bump(l2t, t_curr)) {
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

  struct mom_beach_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);

  // Electron equations.
  struct gkyl_wv_eqn *elc_euler = gkyl_wv_euler_new(ctx.gas_gamma, app_args.use_gpu);

  struct gkyl_vlasov_fluid_species elc = {
    .name = "elc",
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .equation = elc_euler,
    .init = evalElcInit,
    .ctx = &ctx,

    .bcx = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
  };

  // Field.
  struct gkyl_vlasov_field field = {
    .epsilon0 = ctx.epsilon0, .mu0 = ctx.mu0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,
    .limit_em = true,
    
    .init = evalFieldInit,
    .ctx = &ctx,
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

  // Vlasov-Maxwell app.
  struct gkyl_vm app_inp = {

    .cdim = 1, .vdim = 0,
    .lower = { 0.0 },
    .upper = { ctx.Lx },
    .cells = { NX },

    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,
    .cfl_frac = ctx.cfl_frac,

    .num_periodic_dir = 0,
    .periodic_dirs = { },

    .num_species = 0,
    .species = { },

    .num_fluid_species = 1,
    .fluid_species = { elc },

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
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&app_inp);

  // Initial and final simulation times.
  double t_curr = 0.0, t_end = ctx.t_end;

  // Initialize simulation.
  int frame_curr = 0;
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_vlasov_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_vlasov_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n", gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_vlasov_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_vlasov_app_cout(app, stdout, " at time = %g\n", t_curr);
  }
  else {
    gkyl_vlasov_app_apply_ic(app, t_curr);
  }

  // Create trigger for field energy.
  int field_energy_calcs = ctx.field_energy_calcs;
  struct gkyl_tm_trigger fe_trig = { .dt = t_end / field_energy_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_field_energy(&fe_trig, app, t_curr);

  // Create trigger for integrated moments.
  int integrated_mom_calcs = ctx.integrated_mom_calcs;
  struct gkyl_tm_trigger im_trig = { .dt = t_end / integrated_mom_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_integrated_mom(&im_trig, app, t_curr);

  // Create trigger for integrated L2 norm of the distribution function.
  int integrated_L2_f_calcs = ctx.integrated_L2_f_calcs;
  struct gkyl_tm_trigger l2f_trig = { .dt = t_end / integrated_L2_f_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_integrated_L2_f(&l2f_trig, app, t_curr);

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
    gkyl_vlasov_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = gkyl_vlasov_update(app, dt);
    gkyl_vlasov_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      gkyl_vlasov_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_field_energy(&fe_trig, app, t_curr);
    calc_integrated_mom(&im_trig, app, t_curr);
    calc_integrated_L2_f(&l2f_trig, app, t_curr);
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

  calc_field_energy(&fe_trig, app, t_curr);
  calc_integrated_mom(&im_trig, app, t_curr);
  calc_integrated_L2_f(&l2f_trig, app, t_curr);
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
  gkyl_vlasov_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);
  
freeresources:
  // Free resources after simulation completion.
  gkyl_wv_eqn_release(elc_euler);
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
