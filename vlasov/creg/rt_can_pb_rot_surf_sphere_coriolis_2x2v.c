// Coriolis deflection of a bump launched along a meridian on a rotating unit sphere, rotating-frame
// Hamiltonian form H = p^2/2 - omega p_phi (2x2v, collisionless). A Gaussian bump of temperature
// 0.01 just below the equator moves poleward with unit momentum while co-rotating with the sphere;
// in the inertial frame its mean momentum follows the great circle through (theta_0, phi_0) =
// (1.4, pi) with physical velocity (-1, omega sin theta_0), x(t) = x_0 cos(s t) + (v/s) sin(s t)
// with s = |v| = 1.405, and in the rotating frame phi_rot = phi_inertial - omega t, so the centroid
// must reach (theta, phi) = (1.113, 3.167) at t = 0.3
// (thermal corrections are of order T t^2 = 1e-3); the drift in phi relative to the meridian is the
// Coriolis effect. Figures of merit (serendipity p2, 10 x 20 x 11 x 8 cells, 78 steps): at t = 0.3
// the centroid is at (theta, phi) = (1.1234, 3.1690) against the exact (1.1129, 3.1673), errors
// +0.011 and +0.002 rad with cells of 0.039 and 0.105 rad; mass changes by 4e-3.

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

struct sphere_coriolis_ctx {
  double mass; // Neutral mass.
  double charge; // Neutral charge.
  double n0; // Peak number density of the bump.
  double T0; // Temperature of the bump.
  double vt; // Thermal velocity.
  double omega; // Rotation rate of the sphere.
  double theta_0; // Polar angle of the bump.
  double phi_0; // Azimuth of the bump.
  double L_bump; // Squared width of the bump.
  double V_theta_0; // Polar momentum of the bump (poleward).

  int Ntheta; // Cell count (configuration space: polar direction).
  int Nphi; // Cell count (configuration space: azimuthal direction).
  int Nvtheta; // Cell count (velocity space: polar momentum).
  int Nvphi; // Cell count (velocity space: azimuthal momentum).
  double Ltheta; // Domain size (configuration space: polar direction).
  double Lphi; // Domain size (configuration space: azimuthal direction).
  double ptheta_lo; // Lower boundary of the polar-momentum domain.
  double
    ptheta_hi; // Upper boundary of the polar-momentum domain (the Coriolis force raises p_theta).
  double pphi_lo; // Lower boundary of the azimuthal-momentum domain (around omega sin^2 theta_0).
  double pphi_hi; // Upper boundary of the azimuthal-momentum domain.
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

struct sphere_coriolis_ctx
create_ctx(void)
{
  double mass = 1.0; // Neutral mass.
  double charge = 0.0; // Neutral charge.
  double n0 = 1.0; // Peak number density of the bump.
  double T0 = 0.01; // Temperature of the bump.
  double vt = sqrt(T0); // Thermal velocity.
  double omega = 1.0; // Rotation rate of the sphere.
  double theta_0 = 1.4; // Polar angle of the bump.
  double phi_0 = M_PI; // Azimuth of the bump.
  double L_bump = 0.1; // Squared width of the bump.
  double V_theta_0 = -1.0; // Polar momentum of the bump (poleward).

  int Ntheta = 10; // Cell count (configuration space: polar direction).
  int Nphi = 20; // Cell count (configuration space: azimuthal direction).
  int Nvtheta = 11; // Cell count (velocity space: polar momentum).
  int Nvphi = 8; // Cell count (velocity space: azimuthal momentum).
  double Ltheta = (3.0 * M_PI) / 8.0; // Domain size (configuration space: polar direction).
  double Lphi = 2.0 * M_PI; // Domain size (configuration space: azimuthal direction).
  double ptheta_lo = -1.5; // Lower boundary of the polar-momentum domain.
  double ptheta_hi =
    -0.4; // Upper boundary of the polar-momentum domain (the Coriolis force raises p_theta).
  double pphi_lo =
    0.55; // Lower boundary of the azimuthal-momentum domain (around omega sin^2 theta_0).
  double pphi_hi = 1.35; // Upper boundary of the azimuthal-momentum domain.
  int poly_order = 2; // Polynomial order.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 0.3; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs =
    INT_MAX; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct sphere_coriolis_ctx ctx = {
    .mass = mass,
    .charge = charge,
    .n0 = n0,
    .T0 = T0,
    .vt = vt,
    .omega = omega,
    .theta_0 = theta_0,
    .phi_0 = phi_0,
    .L_bump = L_bump,
    .V_theta_0 = V_theta_0,
    .Ntheta = Ntheta,
    .Nphi = Nphi,
    .Nvtheta = Nvtheta,
    .Nvphi = Nvphi,
    .Ltheta = Ltheta,
    .Lphi = Lphi,
    .ptheta_lo = ptheta_lo,
    .ptheta_hi = ptheta_hi,
    .pphi_lo = pphi_lo,
    .pphi_hi = pphi_hi,
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
evalInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct sphere_coriolis_ctx *app = ctx;
  double theta = xn[0], phi = xn[1], p_theta_dot = xn[2], p_phi_dot = xn[3];
  double s2 = sin(theta) * sin(theta);
  double dth = theta - app->theta_0, dph = (phi - app->phi_0) * sin(theta);

  // Gaussian bump at (theta_0, phi_0) carrying a Maxwellian centred on the canonical momenta
  // (V_theta_0, omega sin^2 theta) of a parcel moving along the meridian and co-rotating.
  double n = app->n0 * exp(-((dth * dth) + (dph * dph)) / (0.5 * app->L_bump));
  double dpt = p_theta_dot - app->V_theta_0, dpp = p_phi_dot - (app->omega * s2);
  double efact = ((dpt * dpt) + (dpp * dpp / s2)) / (2.0 * app->T0);

  // Metric determinant times n times the normalized Maxwellian (its measure is 2 pi T sin theta).
  fout[0] = n * exp(-efact) / (2.0 * M_PI * app->T0);
}

void
evalHamiltonian(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct sphere_coriolis_ctx *app = ctx;
  double q_theta = xn[0];
  double p_theta_dot = xn[2], p_phi_dot = xn[3];
  double s2 = sin(q_theta) * sin(q_theta);

  // Rotating-frame Hamiltonian H = p_theta^2/2 + p_phi^2/(2 sin^2 theta) - omega p_phi.
  fout[0] = (0.5 * p_theta_dot * p_theta_dot) + (0.5 * p_phi_dot * p_phi_dot / s2) -
            (app->omega * p_phi_dot);
}

void
evalInvMetric(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct sphere_coriolis_ctx *app = ctx;
  double q_theta = xn[0];
  // Set inverse metric tensor (aa, ab, bb components).
  fout[0] = 1.0;
  fout[1] = 0.0;
  fout[2] = 1.0 / (sin(q_theta) * sin(q_theta));
}

void
evalMetric(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct sphere_coriolis_ctx *app = ctx;
  double q_theta = xn[0];
  // Set metric tensor (aa, ab, bb components).
  fout[0] = 1.0;
  fout[1] = 0.0;
  fout[2] = sin(q_theta) * sin(q_theta);
}

void
evalMetricDet(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct sphere_coriolis_ctx *app = ctx;
  double q_theta = xn[0];
  // Set metric tensor determinant (square root of det g).
  fout[0] = sin(q_theta);
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

  struct sphere_coriolis_ctx ctx = create_ctx(); // Context for initialization functions.

  int NTHETA = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Ntheta);
  int NPHI = APP_ARGS_CHOOSE(app_args.xcells[1], ctx.Nphi);
  int NVTHETA = APP_ARGS_CHOOSE(app_args.vcells[0], ctx.Nvtheta);
  int NVPHI = APP_ARGS_CHOOSE(app_args.vcells[1], ctx.Nvphi);

  int nrank = 1; // Number of processors in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif

  int ccells[] = {NTHETA, NPHI};
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

  // Neutral species.
  struct gkyl_vlasov_kinetic_species neut = {
    .model_id = GKYL_MODEL_CANONICAL_PB,
    .lower = {ctx.ptheta_lo, ctx.pphi_lo},
    .upper = {ctx.ptheta_hi, ctx.pphi_hi},
    .cells = {NVTHETA, NVPHI},

    .hamil = evalHamiltonian,
    .hamil_ctx = &ctx,
    .h_ij = evalMetric,
    .h_ij_ctx = &ctx,
    .h_ij_inv = evalInvMetric,
    .h_ij_inv_ctx = &ctx,
    .det_h = evalMetricDet,
    .det_h_ctx = &ctx,

    .num_init = 1,
    .projection[0] = {.proj_id = GKYL_PROJ_FUNC, .func = evalInit, .ctx_func = &ctx},

    .bcx = {.lower = {.type = GKYL_SPECIES_REFLECT}, .upper = {.type = GKYL_SPECIES_REFLECT}},

    .num_diag_moments = 4,
    .diag_moments = {GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_LTE, GKYL_F_MOMENT_ENERGY},
  };

  // Vlasov app (no field).
  struct gkyl_vm app_inp = {

    .cdim = 2,
    .vdim = 2,
    .lower = {M_PI / 4.0, 0.0},
    .upper = {(M_PI / 4.0) + ctx.Ltheta, ctx.Lphi},
    .cells = {NTHETA, NPHI},

    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,
    .cfl_frac = ctx.cfl_frac,

    .num_periodic_dir = 1,
    .periodic_dirs = {1},

    .num_species = 1,
    .species = {{
      .name = "neut",
      .charge = ctx.charge,
      .mass = ctx.mass,
      .type = GKYL_SPECIES_VLASOV,
      .kinetic = neut,
    }},

    .skip_field = true,

    .parallelism =
      {.use_gpu = app_args.use_gpu, .cuts = {app_args.cuts[0], app_args.cuts[1]}, .comm = comm},
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
