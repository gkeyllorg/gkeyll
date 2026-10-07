// Cross-species relaxation of an electron-ion plasma under LBO collisions (Vlasov, 1x2v).
// Homogeneous electrons (T = 1, drifting at half the ion thermal speed) and ions
// (mass 100, T = 0.5, at rest) with constant collision frequencies nu_ee = nu_ei = 1, nu_ii = 0.1
// and nu_ie = (m_e / m_i) nu_ei, so that the pair exchanges momentum and energy conservatively. Run
// for 20 electron collision times. Particle number, total momentum and total energy are conserved
// to 1e-12. The drift difference decays at nu_ei (measured 1.01) and the temperature difference at
// 2 (nu_ei m_e + nu_ie m_i) / (m_e + m_i) = 0.0396 (measured 0.0396).

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
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

struct cross_relax_ctx {
  // Physical constants (using normalized code units).
  double mass_elc; // Electron mass.
  double mass_ion; // Ion mass.
  double charge_elc; // Electron charge.
  double charge_ion; // Ion charge.

  double n0; // Reference number density.
  double Te; // Electron temperature.
  double Ti; // Ion temperature.
  double ue; // Electron drift velocity (x-direction).
  double ui; // Ion drift velocity (x-direction).
  double nu_ee; // Electron-electron collision frequency.
  double nu_ei; // Electron-ion collision frequency.
  double nu_ii; // Ion-ion collision frequency.
  double nu_ie; // Ion-electron collision frequency.

  // Derived physical quantities (using normalized code units).
  double vte; // Electron thermal velocity.
  double vti; // Ion thermal velocity.

  // Simulation parameters.
  int Nx; // Cell count (configuration space: x-direction).
  int Nvx; // Cell count (velocity space: vx-direction).
  int Nvy; // Cell count (velocity space: vy-direction).
  double Lx; // Domain size (configuration space: x-direction).
  double vx_max_elc; // Electron domain boundary (velocity space: vx-direction).
  double vy_max_elc; // Electron domain boundary (velocity space: vy-direction).
  double vx_max_ion; // Ion domain boundary (velocity space: vx-direction).
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

struct cross_relax_ctx
create_ctx(void)
{
  // Physical constants (using normalized code units).
  double mass_elc = 1.0; // Electron mass.
  double mass_ion = 100.0; // Ion mass.
  double charge_elc = -1.0; // Electron charge.
  double charge_ion = 1.0; // Ion charge.

  double n0 = 1.0; // Reference number density.
  double Te = 1.0; // Electron temperature.
  double Ti = 0.5; // Ion temperature.
  double nu_ee = 1.0; // Electron-electron collision frequency.
  double nu_ei = 1.0; // Electron-ion collision frequency.
  double nu_ii = 0.1; // Ion-ion collision frequency.
  double nu_ie = (mass_elc / mass_ion) * nu_ei; // Ion-electron collision frequency.

  // Derived physical quantities (using normalized code units).
  double vte = sqrt(Te / mass_elc); // Electron thermal velocity.
  double vti = sqrt(Ti / mass_ion); // Ion thermal velocity.
  double ue = 0.5 * vti; // Electron drift velocity (x-direction).
  double ui = 0.0; // Ion drift velocity (x-direction).

  // Simulation parameters.
  int Nx = 2; // Cell count (configuration space: x-direction).
  int Nvx = 16; // Cell count (velocity space: vx-direction).
  int Nvy = 16; // Cell count (velocity space: vy-direction).
  double Lx = 5.0; // Domain size (configuration space: x-direction).
  double vx_max_elc = 5.0 * vte; // Electron domain boundary (velocity space: vx-direction).
  double vy_max_elc = 5.0 * vte; // Electron domain boundary (velocity space: vy-direction).
  double vx_max_ion = 5.0 * vti; // Ion domain boundary (velocity space: vx-direction).
  double vy_max_ion = 5.0 * vti; // Ion domain boundary (velocity space: vy-direction).
  int poly_order = 2; // Polynomial order.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 20.0; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs =
    INT_MAX; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct cross_relax_ctx ctx = {
    .mass_elc = mass_elc,
    .mass_ion = mass_ion,
    .charge_elc = charge_elc,
    .charge_ion = charge_ion,
    .n0 = n0,
    .Te = Te,
    .Ti = Ti,
    .ue = ue,
    .ui = ui,
    .nu_ee = nu_ee,
    .nu_ei = nu_ei,
    .nu_ii = nu_ii,
    .nu_ie = nu_ie,
    .vte = vte,
    .vti = vti,
    .Nx = Nx,
    .Nvx = Nvx,
    .Nvy = Nvy,
    .Lx = Lx,
    .vx_max_elc = vx_max_elc,
    .vy_max_elc = vy_max_elc,
    .vx_max_ion = vx_max_ion,
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
evalElcDensity(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct cross_relax_ctx *app = ctx;

  // Set number density.
  fout[0] = app->n0;
}

void
evalElcTemp(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct cross_relax_ctx *app = ctx;

  // Set isotropic temperature.
  fout[0] = app->Te;
}

void
evalElcVDrift(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct cross_relax_ctx *app = ctx;

  // Set drift velocity.
  fout[0] = app->ue;
  fout[1] = 0.0;
}

void
evalIonDensity(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct cross_relax_ctx *app = ctx;

  // Set number density.
  fout[0] = app->n0;
}

void
evalIonTemp(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct cross_relax_ctx *app = ctx;

  // Set isotropic temperature.
  fout[0] = app->Ti;
}

void
evalIonVDrift(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct cross_relax_ctx *app = ctx;

  // Set drift velocity.
  fout[0] = app->ui;
  fout[1] = 0.0;
}

void
evalNuEE(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct cross_relax_ctx *app = ctx;

  // Set collision frequency.
  fout[0] = app->nu_ee;
}

void
evalNuEI(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct cross_relax_ctx *app = ctx;

  // Set collision frequency.
  fout[0] = app->nu_ei;
}

void
evalNuII(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct cross_relax_ctx *app = ctx;

  // Set collision frequency.
  fout[0] = app->nu_ii;
}

void
evalNuIE(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct cross_relax_ctx *app = ctx;

  // Set collision frequency.
  fout[0] = app->nu_ie;
}

// Snap a trigger whose next time exceeds t_end only by round-off back to t_end.
static void
snap_trigger_to_t_end(struct gkyl_tm_trigger *trig, double t_end)
{
  if (trig->tcurr > t_end && trig->tcurr <= t_end * (1.0 + 1.0e-10)) {
    trig->tcurr = t_end;
  }
}

void
write_data(struct gkyl_tm_trigger *iot, gkyl_vlasov_app *app, double t_curr, bool force_write)
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
calc_field_energy(struct gkyl_tm_trigger *fet, gkyl_vlasov_app *app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(fet, t_curr) || force_calc) {
    gkyl_vlasov_app_calc_field_energy(app, t_curr);
  }
}

void
calc_integrated_mom(
  struct gkyl_tm_trigger *imt, gkyl_vlasov_app *app, double t_curr, bool force_calc
)
{
  if (gkyl_tm_trigger_check_and_bump(imt, t_curr) || force_calc) {
    gkyl_vlasov_app_calc_integrated_mom(app, t_curr);
  }
}

void
calc_integrated_L2_f(
  struct gkyl_tm_trigger *l2t, gkyl_vlasov_app *app, double t_curr, bool force_calc
)
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

  struct cross_relax_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int NVX = APP_ARGS_CHOOSE(app_args.vcells[0], ctx.Nvx);
  int NVY = APP_ARGS_CHOOSE(app_args.vcells[1], ctx.Nvy);

  int nrank = 1; // Number of processors in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif

  int ccells[] = {NX};
  int cdim = sizeof(ccells) / sizeof(ccells[0]);

  int cuts[cdim];
#ifdef GKYL_HAVE_MPI
  for (int d = 0; d < cdim; d++) {
    if (app_args.use_mpi) {
      cuts[d] = app_args.cuts[d];
    } else {
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
    comm = gkyl_nccl_comm_new(&(struct gkyl_nccl_comm_inp){.mpi_comm = MPI_COMM_WORLD});
#else
    printf(" Using -g and -M together requires NCCL.\n");
    assert(0 == 1);
#endif
  } else if (app_args.use_mpi) {
    comm = gkyl_mpi_comm_new(&(struct gkyl_mpi_comm_inp){.mpi_comm = MPI_COMM_WORLD});
  } else {
    comm = gkyl_null_comm_inew(&(struct gkyl_null_comm_inp){.use_gpu = app_args.use_gpu});
  }
#else
  comm = gkyl_null_comm_inew(&(struct gkyl_null_comm_inp){.use_gpu = app_args.use_gpu});
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
    .model_id = GKYL_MODEL_DEFAULT,
    .lower = {-ctx.vx_max_elc, -ctx.vy_max_elc},
    .upper = {ctx.vx_max_elc, ctx.vy_max_elc},
    .cells = {NVX, NVY},

    .num_init = 1,
    .projection[0] =
      {
        .proj_id = GKYL_PROJ_VLASOV_LTE,
        .density = evalElcDensity,
        .ctx_density = &ctx,
        .temp = evalElcTemp,
        .ctx_temp = &ctx,
        .V_drift = evalElcVDrift,
        .ctx_V_drift = &ctx,
        .correct_all_moms = true,
        .use_last_converged = true,
      },
    .collisions =
      {
        .collision_id = GKYL_LBO_COLLISIONS,
        .self_nu = evalNuEE,
        .self_nu_ctx = &ctx,
        .num_cross_collisions = 1,
        .collide_with = {"ion"},
        .cross_nu = {evalNuEI},
        .cross_nu_ctx = {&ctx},
      },

    .num_diag_moments = 4,
    .diag_moments = {GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_LTE},
  };

  // Ions.
  struct gkyl_vlasov_kinetic_species ion = {
    .model_id = GKYL_MODEL_DEFAULT,
    .lower = {-ctx.vx_max_ion, -ctx.vy_max_ion},
    .upper = {ctx.vx_max_ion, ctx.vy_max_ion},
    .cells = {NVX, NVY},

    .num_init = 1,
    .projection[0] =
      {
        .proj_id = GKYL_PROJ_VLASOV_LTE,
        .density = evalIonDensity,
        .ctx_density = &ctx,
        .temp = evalIonTemp,
        .ctx_temp = &ctx,
        .V_drift = evalIonVDrift,
        .ctx_V_drift = &ctx,
        .correct_all_moms = true,
        .use_last_converged = true,
      },
    .collisions =
      {
        .collision_id = GKYL_LBO_COLLISIONS,
        .self_nu = evalNuII,
        .self_nu_ctx = &ctx,
        .num_cross_collisions = 1,
        .collide_with = {"elc"},
        .cross_nu = {evalNuIE},
        .cross_nu_ctx = {&ctx},
      },

    .num_diag_moments = 4,
    .diag_moments = {GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_LTE},
  };

  // Vlasov app (no field).
  struct gkyl_vm app_inp = {

    .cdim = 1,
    .vdim = 2,
    .lower = {0.0},
    .upper = {ctx.Lx},
    .cells = {NX},

    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,
    .cfl_frac = ctx.cfl_frac,

    .num_periodic_dir = 1,
    .periodic_dirs = {0},

    .num_species = 2,
    .species =
      {
        {
          .name = "elc",
          .charge = ctx.charge_elc,
          .mass = ctx.mass_elc,
          .type = GKYL_SPECIES_VLASOV,
          .kinetic = elc,
        },
        {
          .name = "ion",
          .charge = ctx.charge_ion,
          .mass = ctx.mass_ion,
          .type = GKYL_SPECIES_VLASOV,
          .kinetic = ion,
        },
      },

    .skip_field = true,

    .parallelism = {.use_gpu = app_args.use_gpu, .cuts = {app_args.cuts[0]}, .comm = comm},
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
    struct gkyl_app_restart_status status =
      gkyl_vlasov_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_vlasov_app_cout(
        app, stderr, "*** Failed to read restart file! (%s)\n",
        gkyl_array_rio_status_msg(status.io_status)
      );
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_vlasov_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_vlasov_app_cout(app, stdout, " at time = %g\n", t_curr);
  } else {
    gkyl_vlasov_app_apply_ic(app, t_curr);
  }

  // Create trigger for field energy.
  int field_energy_calcs = ctx.field_energy_calcs;
  struct gkyl_tm_trigger fe_trig = {
    .dt = t_end / field_energy_calcs,
    .tcurr = t_curr,
    .curr = frame_curr,
  };

  calc_field_energy(&fe_trig, app, t_curr, false);

  // Create trigger for integrated moments.
  int integrated_mom_calcs = ctx.integrated_mom_calcs;
  struct gkyl_tm_trigger im_trig = {
    .dt = t_end / integrated_mom_calcs,
    .tcurr = t_curr,
    .curr = frame_curr,
  };

  calc_integrated_mom(&im_trig, app, t_curr, false);

  // Create trigger for integrated L2 norm of the distribution function.
  int integrated_L2_f_calcs = ctx.integrated_L2_f_calcs;
  struct gkyl_tm_trigger l2f_trig = {
    .dt = t_end / integrated_L2_f_calcs,
    .tcurr = t_curr,
    .curr = frame_curr,
  };

  calc_integrated_L2_f(&l2f_trig, app, t_curr, false);

  // Create trigger for IO.
  int num_frames = ctx.num_frames;
  struct gkyl_tm_trigger io_trig = {
    .dt = t_end / num_frames,
    .tcurr = frame_curr * (t_end / num_frames),
    .curr = frame_curr,
  };

  write_data(&io_trig, app, t_curr, false);

  // Compute initial guess of maximum stable time-step.
  double dt = t_end - t_curr;

  // The requested time-step is shortened near the end of the simulation so that
  // the final step lands exactly on t_end.
  bool is_dt_clipped = false; // Was the requested dt shortened below the stable dt?
  bool is_last_step = true; // Does the requested dt reach t_end?

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

    // Only a step that took the full requested dt counts as shortened/final.
    bool took_requested_dt = status.dt_actual == dt;
    bool was_dt_clipped = is_dt_clipped && took_requested_dt;
    if (is_last_step && took_requested_dt) {
      // Avoid round-off leaving t_curr just short of t_end.
      t_curr = t_end;
      // Trigger times are accumulated sums and can exceed t_end by round-off.
      // Snap them so the final frame and diagnostics are still produced.
      snap_trigger_to_t_end(&fe_trig, t_end);
      snap_trigger_to_t_end(&im_trig, t_end);
      snap_trigger_to_t_end(&l2f_trig, t_end);
      snap_trigger_to_t_end(&io_trig, t_end);
    } else {
      t_curr += status.dt_actual;
    }

    // Request the next time-step. If the remaining time fits in one stable step,
    // take exactly the remaining time; if it fits in less than two, split it into
    // two equal steps so the final step is never a sliver of the stable dt.
    double t_left = t_end - t_curr;
    dt = status.dt_suggested;
    is_dt_clipped = false;
    is_last_step = false;
    if (t_left <= dt) {
      dt = t_left;
      is_dt_clipped = true;
      is_last_step = true;
    } else if (t_left < 2.0 * dt) {
      dt = 0.5 * t_left;
      is_dt_clipped = true;
    }

    calc_field_energy(&fe_trig, app, t_curr, false);
    calc_integrated_mom(&im_trig, app, t_curr, false);
    calc_integrated_L2_f(&l2f_trig, app, t_curr, false);
    write_data(&io_trig, app, t_curr, false);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    } else if (!was_dt_clipped && status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_vlasov_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_vlasov_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_vlasov_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_vlasov_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_vlasov_app_cout(
          app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max
        );

        calc_field_energy(&fe_trig, app, t_curr, true);
        calc_integrated_mom(&im_trig, app, t_curr, true);
        calc_integrated_L2_f(&l2f_trig, app, t_curr, true);
        write_data(&io_trig, app, t_curr, true);

        break;
      }
    } else {
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
    gkyl_vlasov_app_cout(
      app, stdout, "  Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]
    );
    gkyl_vlasov_app_cout(
      app, stdout, "  Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]
    );
  }
  gkyl_vlasov_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_vlasov_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_vlasov_app_cout(
    app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm
  );
  gkyl_vlasov_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_vlasov_app_cout(
    app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm
  );
  gkyl_vlasov_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

  gkyl_vlasov_app_cout(app, stdout, "Number of write calls %ld\n", stat.n_io);
  double io_tm =
    stat.field_io_tm + stat.species_io_tm + stat.field_diag_io_tm + stat.species_diag_io_tm;
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
