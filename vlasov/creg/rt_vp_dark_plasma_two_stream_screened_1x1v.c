// Two-stream instability of a charged dark-matter pair plasma with the Vlasov-Poisson system
// including self-gravity and a screened (massive dark photon) electrostatic interaction (1x1v). The
// configuration follows rt_vp_dark_plasma_two_stream_1x1v: electrons and positrons (m = 1, |q| = 1,
// vt = 1) each consist of two counter-streaming Maxwellian beams at +-5 vt with density n0/2, the
// electron beams seeded with random-amplitude, random-phase density noise of amplitude 1e-4 in
// modes 1-32, the positrons unperturbed. Here the electrostatic potential is screened, -nabla^2 phi
// + mu_sq phi = rho_c/epsilon0 with mu_sq = 0.04 (screening length 5 lambda_D), and the
// gravitational coupling is alpha_g = 1e-3 (lambda_J = sqrt(1000) lambda_D = 31.6 lambda_D, omega_J
// = 0.0316 omega_pe) so that lambda_D < 1/mu < lambda_J are separated. The box is 10 pi lambda_J
// long (fundamental k lambda_J = 0.2) and the run lasts 10 Jeans times.
// Kinetic linear theory for the screened counter-streaming beams (dielectric 1 + mu_sq/k^2 + beam
// terms): box modes 2-36 are electrostatically unstable, the fastest being mode 24 (k lambda_D =
// 0.15) with gamma = 0.231 omega_pe (half the unscreened 0.473), so the electrostatic energy grows
// at up to 2 gamma = 0.463.
// Figures of merit (serendipity p2, 128 x 64 cells): the electrostatic energy grows at 0.448 over t
// = 20-40 and peaks at 5.25e3 at t = 49.5; the gravitational energy peaks at 2.03e4 at t = 282;
// electron and positron numbers conserved to 5e-13 and total energy (kinetic plus epsilon0/2 of the
// electrostatic energy file, int |grad phi|^2 + mu_sq phi^2, minus half the gravitational energy
// file, int |grad phi_g|^2/alpha_g) to 8e-6.

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

struct dark_plasma_ctx {
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using normalized code units).
  double epsilon0; // Permittivity of free space.
  double alpha_g; // Gravitational coupling 4 pi G epsilon0 m^2/q^2.
  double mu_sq; // Inverse screening length squared (dark photon mass).
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.
  double mass_pos; // Positron mass.
  double charge_pos; // Positron charge.

  double n0; // Reference number density.
  double T; // Temperature.
  double Vx_drift; // Beam drift velocity (x-direction).

  double delta_n; // Applied perturbation amplitude.
  int mode_init; // Initial wave mode to perturb with noise.
  int mode_final; // Final wave mode to perturb with noise.

  // Derived physical quantities (using normalized code units).
  double vte; // Electron thermal velocity.
  double omega_pe; // Electron plasma frequency.
  double lambda_D; // Electron Debye length.
  double lambda_J; // Jeans length.
  double omega_J; // Jeans frequency.

  // Simulation parameters.
  int Nx; // Cell count (configuration space: x-direction).
  int Nvx; // Cell count (velocity space: vx-direction).
  double Lx; // Domain size (configuration space: x-direction).
  double kx; // Perturbed wave number (x-direction).
  double vx_max; // Domain boundary (velocity space: vx-direction).
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

struct dark_plasma_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  double epsilon0 = 1.0; // Permittivity of free space.
  double alpha_g = 1.0e-3; // Gravitational coupling 4 pi G epsilon0 m^2/q^2.
  double mu_sq = 4.0e-2; // Inverse screening length squared (dark photon mass).
  double mass_elc = 1.0; // Electron mass.
  double charge_elc = -1.0; // Electron charge.
  double mass_pos = 1.0; // Positron mass.
  double charge_pos = 1.0; // Positron charge.

  double n0 = 1.0; // Reference number density.
  double T = 1.0; // Temperature.
  double Vx_drift = 5.0; // Beam drift velocity (x-direction).

  double delta_n = 1.0e-4; // Applied perturbation amplitude.
  int mode_init = 1; // Initial wave mode to perturb with noise.
  int mode_final = 32; // Final wave mode to perturb with noise.

  // Derived physical quantities (using normalized code units).
  double vte = sqrt(T / mass_elc); // Electron thermal velocity.
  double omega_pe =
    sqrt((charge_elc * charge_elc) * n0 / (epsilon0 * mass_elc)); // Electron plasma frequency.
  double lambda_D = vte / omega_pe; // Electron Debye length.
  double lambda_J = lambda_D / sqrt(alpha_g); // Jeans length.
  double omega_J = vte / lambda_J; // Jeans frequency.

  // Simulation parameters.
  int Nx = 128; // Cell count (configuration space: x-direction).
  int Nvx = 64; // Cell count (velocity space: vx-direction).
  double kx =
    0.2 / lambda_J; // Perturbed wave number (x-direction); the fundamental mode of the box.
  double Lx = 2.0 * pi / kx; // Domain size (configuration space: x-direction).
  double vx_max = 32.0 * vte; // Domain boundary (velocity space: vx-direction).
  int poly_order = 2; // Polynomial order.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 10.0 / omega_J; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs =
    INT_MAX; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct dark_plasma_ctx ctx = {
    .pi = pi,
    .epsilon0 = epsilon0,
    .alpha_g = alpha_g,
    .mu_sq = mu_sq,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .mass_pos = mass_pos,
    .charge_pos = charge_pos,
    .n0 = n0,
    .T = T,
    .Vx_drift = Vx_drift,
    .delta_n = delta_n,
    .mode_init = mode_init,
    .mode_final = mode_final,
    .vte = vte,
    .omega_pe = omega_pe,
    .lambda_D = lambda_D,
    .lambda_J = lambda_J,
    .omega_J = omega_J,
    .Nx = Nx,
    .Nvx = Nvx,
    .Lx = Lx,
    .kx = kx,
    .vx_max = vx_max,
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
evalElcDensityInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct dark_plasma_ctx *app = ctx;
  double x = xn[0];

  double pi = app->pi;
  double n0 = app->n0;
  double delta_n = app->delta_n;
  double kx = app->kx;

  // Random-amplitude, random-phase noise in modes mode_init-mode_final. The
  // generator is re-seeded at every point so every point sees the same noise.
  srand(0);
  double perturb = 0.0;
  for (int i = app->mode_init; i <= app->mode_final; i++) {
    double amplitude = (double)rand() / RAND_MAX;
    double phase = 2.0 * pi * (double)rand() / RAND_MAX;
    perturb += delta_n * amplitude * cos(i * kx * x + phase);
  }

  // Set electron beam number density.
  fout[0] = 0.5 * (1.0 + perturb) * n0;
}

void
evalTempInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct dark_plasma_ctx *app = ctx;

  // Set isotropic temperature.
  fout[0] = app->T;
}

void
evalVDriftLeftInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct dark_plasma_ctx *app = ctx;

  // Set left-going beam drift velocity.
  fout[0] = -app->Vx_drift;
}

void
evalVDriftRightInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct dark_plasma_ctx *app = ctx;

  // Set right-going beam drift velocity.
  fout[0] = app->Vx_drift;
}

void
evalPosDensityInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct dark_plasma_ctx *app = ctx;

  // Set positron beam number density.
  fout[0] = 0.5 * app->n0;
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

  struct dark_plasma_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int NVX = APP_ARGS_CHOOSE(app_args.vcells[0], ctx.Nvx);

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

  // Electrons: two counter-streaming beams with density noise.
  struct gkyl_vlasov_kinetic_species elc = {
    .lower = {-ctx.vx_max},
    .upper = {ctx.vx_max},
    .cells = {NVX},

    .num_init = 2,
    .projection[0] =
      {
        .proj_id = GKYL_PROJ_VLASOV_LTE,
        .density = evalElcDensityInit,
        .ctx_density = &ctx,
        .temp = evalTempInit,
        .ctx_temp = &ctx,
        .V_drift = evalVDriftRightInit,
        .ctx_V_drift = &ctx,
        .correct_all_moms = true,
      },
    .projection[1] =
      {
        .proj_id = GKYL_PROJ_VLASOV_LTE,
        .density = evalElcDensityInit,
        .ctx_density = &ctx,
        .temp = evalTempInit,
        .ctx_temp = &ctx,
        .V_drift = evalVDriftLeftInit,
        .ctx_V_drift = &ctx,
        .correct_all_moms = true,
      },

    .num_diag_moments = 3,
    .diag_moments = {GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2},
  };

  // Positrons: two unperturbed counter-streaming beams.
  struct gkyl_vlasov_kinetic_species pos = {
    .lower = {-ctx.vx_max},
    .upper = {ctx.vx_max},
    .cells = {NVX},

    .num_init = 2,
    .projection[0] =
      {
        .proj_id = GKYL_PROJ_VLASOV_LTE,
        .density = evalPosDensityInit,
        .ctx_density = &ctx,
        .temp = evalTempInit,
        .ctx_temp = &ctx,
        .V_drift = evalVDriftRightInit,
        .ctx_V_drift = &ctx,
        .correct_all_moms = true,
      },
    .projection[1] =
      {
        .proj_id = GKYL_PROJ_VLASOV_LTE,
        .density = evalPosDensityInit,
        .ctx_density = &ctx,
        .temp = evalTempInit,
        .ctx_temp = &ctx,
        .V_drift = evalVDriftLeftInit,
        .ctx_V_drift = &ctx,
        .correct_all_moms = true,
      },

    .num_diag_moments = 3,
    .diag_moments = {GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2},
  };

  // Field: screened electrostatics and self-gravity.
  struct gkyl_vlasov_field field = {
    .epsilon0 = ctx.epsilon0,
    .alpha_g = ctx.alpha_g,
    .mu_sq = ctx.mu_sq,

    .poisson_bcs = {.lo_type = {GKYL_POISSON_PERIODIC}, .up_type = {GKYL_POISSON_PERIODIC}},
  };

  // Vlasov-Poisson app.
  struct gkyl_vm app_inp = {

    .cdim = 1,
    .vdim = 1,
    .lower = {-0.5 * ctx.Lx},
    .upper = {0.5 * ctx.Lx},
    .cells = {NX},

    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,
    .cfl_frac = ctx.cfl_frac,

    .num_periodic_dir = 1,
    .periodic_dirs = {0},

    .num_species = 2,
    .species =
      {{
         .name = "elc",
         .charge = ctx.charge_elc,
         .mass = ctx.mass_elc,
         .type = GKYL_SPECIES_VLASOV,
         .kinetic = elc,
       },
       {
         .name = "pos",
         .charge = ctx.charge_pos,
         .mass = ctx.mass_pos,
         .type = GKYL_SPECIES_VLASOV,
         .kinetic = pos,
       }},

    .field = field,
    .is_electrostatic = true,

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
