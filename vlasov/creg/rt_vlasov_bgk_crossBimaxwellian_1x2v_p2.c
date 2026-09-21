#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_vlasov.h>
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

struct bgk_cross_ctx
{
  // Mathematical constants (dimensionless).
  double pi;

  int cdim, vdim; // Dimensionality.

  // Physical constants (using physical units).
  double epsilon0; // Permittivity of free space.
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.
  double mass_ion; // Proton mass.
  double charge_ion; // Proton charge.

  double n0; // Reference number density (1 / m^3).

  double T_par_elc; // Parallel electron temperature.
  double T_par_ion; // Parallel ion temperature.
  double alpha; // Ratio of perpendicular to parallel temperatures.

  // Derived physical quantities (using physical units).
  double T_perp_elc; // Perpendicular electron temperature.
  double T_perp_ion; // Perpendicular ion temperature.

  double vte_par; // Parallel electron thermal velocity.
  double vti_par; // Parallel ion thermal velocity.
  double vte_perp; // Perpendicular electron thermal velocity.
  double vti_perp; // Perpendicular ion thermal velocity.

  double log_lambda_elc; // Electron Coulomb logarithm.
  double log_lambda_ion; // Ion Coulomb logarithm.
  double nu_elc; // Electron collision frequency.
  double nu_ion; // Ion collision frequency.

  double Te; // Electron temperature.
  double Ti; // Ion temperature.
  double vte; // Electron thermal velocity.
  double vti; // Ion thermal velocity.
  double ue_par; // Electron bulk velocity (x-direction).
  double ui_par; // Ion bulk velocity (x-direction).

  // Simulation parameters.
  int Nx; // Cell count (configuration space: x-direction).
  int Nvx; // Cell count (velocity space: vx-direction).
  int Nvy; // Cell count (velocity space: vy-direction).
  double Lx; // Domain size (configuration space: x-direction).
  double vx_max_elc; // Electron domain boundary (velocity space: vx-direction).
  double vx_max_ion; // Ion domain boundary (velocity space: vx-direction).
  double vy_max_elc; // Electron domain boundary (velocity space: vy-direction).
  double vy_max_ion; // Ion domain boundary (velocity space: vy-direction).
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

struct bgk_cross_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = GKYL_PI;

  int cdim = 1, vdim = 2; // Dimensionality.

  // Physical constants (using physical units).
  double epsilon0 = GKYL_EPSILON0; // Permittivity of free space.
  double mass_elc = GKYL_ELECTRON_MASS; // Electron mass.
  double mass_ion = 2.014 * GKYL_PROTON_MASS; // Proton mass.
  double charge_elc = -GKYL_ELEMENTARY_CHARGE; // Electron charge.
  double charge_ion = GKYL_ELEMENTARY_CHARGE; // Proton charge.

  double T_par_elc = 300.0 * GKYL_ELEMENTARY_CHARGE; // Parallel electron temperature.
  double T_par_ion = 200.0 * GKYL_ELEMENTARY_CHARGE; // Parallel ion temperature.
  double alpha = 1.3; // Ratio of perpendicular to parallel temperatures.

  double n0 = 7.0e19; //  Reference number density (1 / m^3).

  // Derived physical quantities (using physical units).
  double T_perp_elc = alpha * T_par_elc; // Perpendicular electron temperature.
  double T_perp_ion = alpha * T_par_ion; // Perpendicular ion temperature.

  double vte_par = sqrt(T_par_elc / mass_elc); // Parallel electron thermal velocity.
  double vti_par = sqrt(T_par_ion / mass_ion); // Parallel ion thermal velocity.
  double vte_perp = sqrt(T_perp_elc / mass_elc); // Perpendicular electron thermal velocity.
  double vti_perp = sqrt(T_perp_ion / mass_ion); // Perpendicular ion thermal velocity.

  double log_lambda_elc = 6.6 - 0.5 * log(n0 / 1.0e20) + 1.5 * log(T_par_elc / charge_ion); // Electron Coulomb logarithm.
  double log_lambda_ion = 6.6 - 0.5 * log(n0 / 1.0e20) + 1.5 * log(T_par_ion / charge_ion); // Ion Coulomb logarithm.
  double nu_elc = log_lambda_elc * pow(charge_ion, 4.0) * n0 /
    (12.0 * pow(GKYL_PI, 3.0 / 2.0) * pow(epsilon0, 2.0) * sqrt(mass_elc) * pow(T_par_elc, 3.0 / 2.0)); // Electron collision frequency.
  double nu_ion = log_lambda_ion * pow(charge_ion, 4.0) * n0 /
    (12.0 * pow(GKYL_PI, 3.0 / 2.0) * pow(epsilon0, 2.0) * sqrt(mass_ion) * pow(T_par_ion, 3.0 / 2.0)); // Ion collision frequency.
  
  double Te = (T_par_elc + (2.0 * T_perp_elc)) / 3.0; // Electron temperature.
  double Ti = (T_par_ion + (2.0 * T_perp_ion)) / 3.0; // Ion temperature.
  double vte = sqrt(Te / mass_elc); // Electron thermal velocity.
  double vti = sqrt(Ti / mass_ion); // Ion thermal velocity.
  double ue_par = 0.5 * sqrt(mass_elc / mass_ion) * vte; // Parallel electron velocity.
  double ui_par = 50.0 * (mass_elc / mass_ion) * vti; // Parallel ion velocity.

  // Simulation parameters.
  int Nx = 1; // Cell count (configuration space: x-direction).
  int Nvx = 16; // Cell count (velocity space: vx-direction).
  int Nvy = 16; // Cell count (velocity space: vy-direction).
  double Lx = 4.0; // Domain size (configuration space: x-direction).
  double vx_max_elc = 6.0 * vte_par; // First electron domain boundary (velocity space: vx-direction).
  double vx_max_ion = 6.0 * vti_par; // Second ion domain boundary (velocity space: vx-direction).
  double vy_max_elc = 6.0 * vte_par; // First electron domain boundary (velocity space: vy-direction).
  double vy_max_ion = 6.0 * vti_par; // Second ion domain boundary (velocity space: vy-direction).
  int poly_order = 2; // Polynomial order.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 200 / nu_ion; // Final simulation time.
  int num_frames = 6000; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs = INT_MAX; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.
  
  struct bgk_cross_ctx ctx = {
    .pi = pi,
    .cdim = cdim,
    .vdim = vdim,
    .epsilon0 = epsilon0,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .mass_ion = mass_ion,
    .charge_ion = charge_ion,
    .T_par_elc = T_par_elc,
    .T_par_ion = T_par_ion,
    .n0 = n0,
    .T_perp_elc = T_perp_elc,
    .T_perp_ion = T_perp_ion,
    .vte_par = vte_par,
    .vti_par = vti_par,
    .vte_perp = vte_perp,
    .vti_perp = vti_perp,
    .log_lambda_elc = log_lambda_elc,
    .nu_elc = nu_elc,
    .log_lambda_ion = log_lambda_ion,
    .nu_ion = nu_ion,
    .Te = Te,
    .Ti = Ti,
    .vte = vte,
    .vti = vti,
    .ue_par = ue_par,
    .ui_par = ui_par,
    .Nx = Nx,
    .Nvx = Nvx,
    .Nvy = Nvy,
    .Lx = Lx,
    .vx_max_elc = vx_max_elc,
    .vx_max_ion = vx_max_ion,
    .vy_max_elc = vy_max_elc,
    .vy_max_ion = vy_max_ion,
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
evalElcInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct bgk_cross_ctx *app = ctx;

  double vx = xn[1], vy = xn[2];

  double pi = app->pi;
  double n0 = app->n0;
  double vte_par = app->vte_par;
  double vte_perp = app->vte_perp;
  double ue_par = app->ue_par;

  double vsq_par = pow((vx - ue_par), 2.0);
  double vsq_perp = pow(vy, 2.0);

  double n = (n0 / (2.0 * pow(pi, 3.0 / 2.0) * vte_par * pow(vte_perp, 2.0)) ) *
    exp( (-vsq_par/pow(vte_par, 2.0)) + (-vsq_perp/pow(vte_perp, 2.0)) ); // Bimaxwellian distribution function.

  // Set distribution function.
  fout[0] = n;
}

void
evalIonInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct bgk_cross_ctx *app = ctx;

  double vx = xn[1], vy = xn[2];

  double pi = app->pi;
  double n0 = app->n0;
  double vti_par = app->vti_par;
  double vti_perp = app->vti_perp;
  double ui_par = app->ui_par;

  double vsq_par = pow((vx - ui_par), 2.0);
  double vsq_perp = pow(vy, 2.0);

  double n = (n0 / (2.0 * pow(pi, 3.0 / 2.0) * vti_par * pow(vti_perp, 2.0)) ) *
    exp( (-vsq_par/pow(vti_par, 2.0)) + (-vsq_perp/pow(vti_perp, 2.0)) ); // Bimaxwellian distribution function.

  // Set distribution function.
  fout[0] = n;
}

void
evalNuEE(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct bgk_cross_ctx *app = ctx;

  double nu_elc = app->nu_elc;

  // Set collision frequency.
  fout[0] = nu_elc;
}

void evalNuEI(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct bgk_cross_ctx *app = ctx;

  double nu_elc = app->nu_elc;

  // Set collision frequency.
  fout[0] = sqrt(2.0) * nu_elc;
}

void
evalNuII(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct bgk_cross_ctx *app = ctx;

  double nu_ion = app->nu_ion;

  // Set collision frequency.
  fout[0] = nu_ion;
}

void
evalNuIE(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct bgk_cross_ctx *app = ctx;

  double mass_elc = app->mass_elc;
  double mass_ion = app->mass_ion;
  double nu_ion = app->nu_ion;

  // Set collision frequency.
  fout[0] = sqrt(2.0) * sqrt(mass_elc / mass_ion) * nu_ion;
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_vlasov_app* app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_write) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_vlasov_app_write(app, t_curr, frame);
    gkyl_vlasov_app_write_field_energy(app);
    gkyl_vlasov_app_write_integrated_mom(app);
    gkyl_vlasov_app_write_integrated_L2_f(app);
    gkyl_vlasov_app_write_mom(app, t_curr, frame);
  }
}

void
calc_field_energy(struct gkyl_tm_trigger* fet, gkyl_vlasov_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(fet, t_curr) || force_calc) {
    gkyl_vlasov_app_calc_field_energy(app, t_curr);
  }
}

void
calc_integrated_mom(struct gkyl_tm_trigger* imt, gkyl_vlasov_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(imt, t_curr) || force_calc) {
    gkyl_vlasov_app_calc_integrated_mom(app, t_curr);
  }
}

void
calc_integrated_L2_f(struct gkyl_tm_trigger* l2t, gkyl_vlasov_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(l2t, t_curr) || force_calc) {
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

  struct bgk_cross_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int NVX = APP_ARGS_CHOOSE(app_args.vcells[0], ctx.Nvx);
  int NVY = APP_ARGS_CHOOSE(app_args.vcells[1], ctx.Nvy);

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

  // Electrons.
  struct gkyl_vlasov_kinetic_species elc = {
    .lower = { -ctx.vx_max_elc, -ctx.vy_max_elc },
    .upper = { ctx.vx_max_elc, ctx.vy_max_elc }, 
    .cells = { NVX, NVY },
    .num_init = 1, 
    .projection[0] = {
      .proj_id = GKYL_PROJ_FUNC,
      .func = evalElcInit,
      .ctx_func = &ctx,
    },

    .correct = {
      .correct_all_moms = true,
      .use_last_converged = true,
      .iter_eps = 1e-12,
      .max_iter = 10,
    },

    .collisions =  {
      .collision_id = GKYL_BGK_COLLISIONS,
      .self_nu = evalNuEE,
      .self_nu_ctx = &ctx,
      .num_cross_collisions = 1,
      .collide_with = { "ion" },
      .cross_nu = { evalNuEI },
      .cross_nu_ctx = { &ctx },
      // Reference values for log Lambda
      .den_ref = ctx.n0,
      .temp_ref = ctx.Te,
    },
    
    .num_diag_moments = 3,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2IJ },
  };

  // Second neutral species.
  struct gkyl_vlasov_kinetic_species ion = {
    .lower = { -ctx.vx_max_ion, -ctx.vy_max_ion },
    .upper = { ctx.vx_max_ion, ctx.vy_max_ion }, 
    .cells = { NVX, NVY },

    .num_init = 1, 
    .projection[0] = {
      .proj_id = GKYL_PROJ_FUNC,
      .func = evalIonInit,
      .ctx_func = &ctx,
    },

    .correct = {
      .correct_all_moms = true,
      .use_last_converged = true,
      .iter_eps = 1e-12,
      .max_iter = 10,
    },

    .collisions =  {
      .collision_id = GKYL_BGK_COLLISIONS,
      .self_nu = evalNuII,
      .self_nu_ctx = &ctx,
      .num_cross_collisions = 1,
      .collide_with = { "elc" },
      .cross_nu = { evalNuIE },
      .cross_nu_ctx = { &ctx },
      // Reference values for log Lambda
      .den_ref = ctx.n0,
      .temp_ref = ctx.Ti,
    },
    
    .num_diag_moments = 3,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2IJ },
  };

    // Vlasov-Maxwell app.
  struct gkyl_vm app_inp = {
    .name = "vlasov_bgk_crossBimax_1x2v_p2",

    .cdim = 1, .vdim = 2,
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
      { .name = "ion", .charge = ctx.charge_ion, .mass = ctx.mass_ion,
        .type = GKYL_SPECIES_VLASOV, .kinetic = ion },
    },

    .skip_field = true,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0] },
      .comm = comm,
    },
  };

  // Create app object.
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

  calc_field_energy(&fe_trig, app, t_curr, false);

  // Create trigger for integrated moments.
  int integrated_mom_calcs = ctx.integrated_mom_calcs;
  struct gkyl_tm_trigger im_trig = { .dt = t_end / integrated_mom_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_integrated_mom(&im_trig, app, t_curr, false);

  // Create trigger for integrated L2 norm of the distribution function.
  int integrated_L2_f_calcs = ctx.integrated_L2_f_calcs;
  struct gkyl_tm_trigger l2f_trig = { .dt = t_end / integrated_L2_f_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_integrated_L2_f(&l2f_trig, app, t_curr, false);

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

    calc_field_energy(&fe_trig, app, t_curr, false);
    calc_integrated_mom(&im_trig, app, t_curr, false);
    calc_integrated_L2_f(&l2f_trig, app, t_curr, false);
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

        calc_field_energy(&fe_trig, app, t_curr, true);
        calc_integrated_mom(&im_trig, app, t_curr, true);
        calc_integrated_L2_f(&l2f_trig, app, t_curr, true);
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
  calc_integrated_L2_f(&l2f_trig, app, t_curr, false);
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