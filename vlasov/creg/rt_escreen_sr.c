#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#ifdef GKYL_HAVE_NCCL
#include <gkyl_nccl_comm.h>
#endif
#endif

#include <rt_arg_parse.h>

struct escreen_ctx 
{
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using normalized code units).
  double epsilon0; // Permittivity of free space.
  double mu0; // Permeability of free space.
  double mass_ion; // Positron mass.
  double charge_ion; // Positron charge.
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.

  double n0; // Reference density.
  double T; // Temperature (units of mc^2).
  double E0; // Electric field amplitude.
  double noise_amp; // Noise level for perturbation.
  int mode_init; // Initial wave mode to perturb with noise.
  int mode_final; // Final wave mode to perturb with noise.

  // Simulation parameters.
  int Nx; // Cell count (configuration space: x-direction).
  int Npx; // Cell count (momentum space: px-direction).
  double Lx; // Domain size (configuration space: x-direction).
  double px_max; // Domain boundary (momentum space: px-direction).
  int poly_order; // Polynomial order.
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
   int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct escreen_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  double epsilon0 = 1.0; // Permittivity of free space.
  double mu0 = 1.0; // Permeability of free space.
  double mass_ion = 1.0; // Positron mass.
  double charge_ion = 1.0; // Positron charge.
  double mass_elc = 1.0; // Electron mass.
  double charge_elc = -1.0; // Electron charge.

  double n0 = 1.0; // Reference density.
  double T = 1.0; // Reference temperature (in units of mc^2).

  double E0 = 100.0; // Electric field amplitude
  double noise_amp = 1.0e-6 * E0; // Noise level for perturbation.
  int mode_init = 1; // Initial wave mode to perturb with noise.
  int mode_final = 4; // Final wave mode to perturb with noise.

  // Simulation parameters.
  int Nx = 16; // Cell count (configuration space: x-direction).
  int Npx = 512; // Cell count (velocity space: vx-direction).
  double Lx = 1.0; // Domain size (configuration space: x-direction).
  double px_max = 768.0; // Domain boundary (velocity space: vx-direction).
  int poly_order = 2; // Polynomial order.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 10.0; // Final simulation time.
  int num_frames = 2; // Number of output frames.
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  return (struct escreen_ctx) {
    .pi = pi,
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .mass_ion = mass_ion,
    .charge_ion = charge_ion,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .n0 = n0, 
    .T = T,
    .E0 = E0, 
    .noise_amp = noise_amp,
    .mode_init = mode_init,
    .mode_final = mode_final,
    .Nx = Nx,
    .Npx = Npx,
    .Lx = Lx,
    .px_max = px_max,
    .poly_order = poly_order,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .num_frames = num_frames,
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };
}

void
evalDensityInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct escreen_ctx *app = ctx;

  double n = app->n0;

  // Set density (no perturbation).
  fout[0] = n;
}

void
evalDensityPerturbInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct escreen_ctx *app = ctx;
  double x = xn[0];

  // Perturbation to initial density to satisfy div(E) = rho_c
  double n = app->n0;
  double pi = app->pi;  
  double noise_amp = app->noise_amp;
  double mode_init = app->mode_init;
  double mode_final = app->mode_final;

  double E_x = app->E0;
  double alpha = noise_amp * E_x; // Applied amplitude.
  double Lx = app->Lx;
  double kx = 2.0 * pi / Lx; // Wave number (x-direction).
  pcg64_random_t rng = gkyl_pcg64_init(0); // Random number generator.
  for (int i = mode_init; i < mode_final; i++) {
    n -= alpha * gkyl_pcg64_rand_double(&rng) * cos(i * kx * x + 2.0 * pi * gkyl_pcg64_rand_double(&rng)); 
  }

  // Set density (no perturbation).
  fout[0] = n;
}

void
evalTempInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct escreen_ctx *app = ctx;

  double T = app->T;

  // Set temperature.
  fout[0] = T;
}

void
evalVDriftInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct escreen_ctx *app = ctx;

  // Set drift (four-) velocity.
  fout[0] = 0.0;
}

void
evalFieldFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct escreen_ctx *app = ctx;
  double x = xn[0];
  double pi = app->pi;

  double noise_amp = app->noise_amp;
  double mode_init = app->mode_init;
  double mode_final = app->mode_final;

  double E_x = app->E0;
  double alpha = noise_amp * E_x; // Applied amplitude.
  double Lx = app->Lx;
  double kx = 2.0 * pi / Lx; // Wave number (x-direction).
  pcg64_random_t rng = gkyl_pcg64_init(0); // Random number generator.

  // Perturbation to initial electric field
  for (int i = mode_init; i < mode_final; i++) {
    E_x += alpha * gkyl_pcg64_rand_double(&rng) * sin(i * kx * x + 2.0 * pi * gkyl_pcg64_rand_double(&rng)); 
  }

  fout[0] = E_x; fout[1] = 0.0, fout[2] = 0.0;
  fout[3] = 0.0; fout[4] = 0.0; fout[5] = 0.0;
  fout[6] = 0.0; fout[7] = 0.0;
}

void

calc_integrated_diagnostics(struct gkyl_tm_trigger* iot, gkyl_vlasov_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_calc) {
    gkyl_vlasov_app_calc_field_energy(app, t_curr);
  }
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_vlasov_app* app, double t_curr, bool force_write)
{
  bool trig_now = gkyl_tm_trigger_check_and_bump(iot, t_curr);
  if (trig_now || force_write) {
    int frame = (!trig_now) && force_write? iot->curr : iot->curr-1;

    gkyl_vlasov_app_write(app, t_curr, frame);
    gkyl_vlasov_app_write_mom(app, t_curr, frame);

    gkyl_vlasov_app_calc_field_energy(app, t_curr);
    gkyl_vlasov_app_write_field_energy(app);
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

  struct escreen_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int NPX = APP_ARGS_CHOOSE(app_args.vcells[0], ctx.Npx);

  int nrank = 1; // Number of processors in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif  

  int ccells[] = { NX };
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

  // electrons
  struct gkyl_vlasov_kinetic_species elc = {
    .model_id = GKYL_MODEL_SR,
    .lower = { -ctx.px_max },
    .upper = { ctx.px_max }, 
    .cells = { NPX },

    .num_init = 1, 
    .projection[0] = {
      .proj_id = GKYL_PROJ_VLASOV_LTE,
      .density = evalDensityPerturbInit,
      .ctx_density = &ctx,
      .temp = evalTempInit,
      .ctx_temp = &ctx,
      .V_drift = evalVDriftInit,
      .ctx_V_drift = &ctx,
      .correct_all_moms = true,
      .use_last_converged = true, 
    },

    // Source is the same as initial condition without perturbation.
    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1, 
      .projection[0] = {
        .proj_id = GKYL_PROJ_VLASOV_LTE,
        .density = evalDensityInit,
        .ctx_density = &ctx,
        .temp = evalTempInit,
        .ctx_temp = &ctx,
        .V_drift = evalVDriftInit,
        .ctx_V_drift = &ctx,
        .correct_all_moms = true,
        .use_last_converged = true, 
      },
    },

    .num_diag_moments = 2,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1 },
  };

  // positrons
  struct gkyl_vlasov_kinetic_species pos = {
    .model_id = GKYL_MODEL_SR,
    .lower = { -ctx.px_max },
    .upper = { ctx.px_max }, 
    .cells = { NPX },

    .num_init = 1, 
    // No perturbation in positron initial conditions. 
    .projection[0] = {
      .proj_id = GKYL_PROJ_VLASOV_LTE,
      .density = evalDensityInit,
      .ctx_density = &ctx,
      .temp = evalTempInit,
      .ctx_temp = &ctx,
      .V_drift = evalVDriftInit,
      .ctx_V_drift = &ctx,
      .correct_all_moms = true,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1, 
      .projection[0] = {
        .proj_id = GKYL_PROJ_VLASOV_LTE,
        .density = evalDensityInit,
        .ctx_density = &ctx,
        .temp = evalTempInit,
        .ctx_temp = &ctx,
        .V_drift = evalVDriftInit,
        .ctx_V_drift = &ctx,
        .correct_all_moms = true,
      },
    },

    .num_diag_moments = 2,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1 },
  };

  // field
  struct gkyl_vlasov_field field = {
    .epsilon0 = ctx.epsilon0, .mu0 = ctx.mu0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,

    .ctx = &ctx,
    .init = evalFieldFunc,
  };

  // Vlasov-Maxwell app.
  struct gkyl_vm app_inp = {
    
    .cdim = 1, .vdim = 1,
    .lower = { -0.5 * ctx.Lx },
    .upper = { 0.5 * ctx.Lx },
    .cells = { NX },

    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,
    .cfl_frac = ctx.cfl_frac,

    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },

    .num_species = 2,
    .species = {
      { .name = "elc", .charge = ctx.charge_elc, .mass = ctx.mass_elc,
        .type = GKYL_SPECIES_VLASOV, .kinetic = elc },
      { .name = "pos", .charge = ctx.charge_ion, .mass = ctx.mass_ion,
        .type = GKYL_SPECIES_VLASOV, .kinetic = pos },
    },

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
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&app_inp);

  // Initial and final simulation times.
  int frame_curr = 0;
  double t_curr = 0.0, t_end = ctx.t_end;
  // Initialize simulation.
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_vlasov_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_vlasov_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n",
        gkyl_array_rio_status_msg(status.io_status));
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
  
  // Create triggers for IO.
  int num_frames = ctx.num_frames, num_int_diag_calc = ctx.int_diag_calc_num;
  struct gkyl_tm_trigger trig_write = { .dt = t_end/num_frames, .tcurr = t_curr, .curr = frame_curr };
  struct gkyl_tm_trigger trig_calc_intdiag = { .dt = t_end/GKYL_MAX2(num_frames, num_int_diag_calc),
    .tcurr = t_curr, .curr = frame_curr };

  // Write out ICs (if restart, it overwrites the restart frame).
  calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, false);
  write_data(&trig_write, app, t_curr, false);

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

    calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, t_curr > t_end);
    write_data(&trig_write, app, t_curr, t_curr > t_end);

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
        calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, true);
        write_data(&trig_write, app, t_curr, true);
        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }

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
