// Circular geodesic of a phase-space blob around a Kerr black hole with spin a = 0.8 M
// (prograde orbit) in ingoing Kerr-Schild coordinates (equatorial r-phi preset geometry) through
// the GR triad Vlasov solver (2x2v, collisionless). The blob starts at r0 = 8 M with the tetrad
// momenta of the prograde circular orbit (Boyer-Lindquist E and L converted to the Kerr-Schild
// covariant momenta, which have a radial component because g_tr is nonzero), so its density
// centroid must stay at r0 and advance in phi at the orbital frequency Omega = sqrt(M) / (r0^(3/2)
// + a sqrt(M)) over t = 25 M; the Kerr and Schwarzschild frequencies differ by 3.5 percent. Figures
// of merit (serendipity p2, 12 x 24 x 6 x 6 cells, 137 steps): at t = 25 M the density centroid is
// at (r, phi) = (7.980, 4.215) against the exact (8, 4.209), i.e. within 0.06 cells in r and 0.07
// cells in phi after Omega t = 1.067 rad; the centroid lags the Schwarzschild driver's by 0.041 rad
// against the exact frame-dragging difference of 0.038 rad; mass changes by 3.5e-4 through the
// absorbing walls.

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

struct geodesic_ctx {
  double mass; // Neutral mass.
  double charge; // Neutral charge.
  double mass_bh; // Black hole mass M.
  double spin_bh; // Black hole spin a = J/M.
  double r0; // Orbit radius.
  double phi0; // Initial azimuth of the blob.

  double v_c; // Orbital velocity parameter sqrt(M / r0).
  double E_c; // Specific energy of the prograde circular orbit (Boyer-Lindquist).
  double L_c; // Specific angular momentum of the circular orbit.
  double g_tt; // Boyer-Lindquist metric component g_tt at r0 (equator).
  double g_tp; // Boyer-Lindquist metric component g_tphi at r0.
  double g_pp; // Boyer-Lindquist metric component g_phiphi at r0.
  double det_g; // Determinant of the (t, phi) metric block.
  double p_t_up; // Contravariant p^t of the orbit.
  double p_phi_up; // Contravariant p^phi of the orbit.
  double h_rr; // Kerr-Schild spatial metric h_rr at r0.
  double h_rp; // Kerr-Schild spatial metric h_rphi at r0.
  double h_pp; // Kerr-Schild spatial metric h_phiphi at r0.
  double p_r_cov; // Covariant Kerr-Schild p_r of the orbit (nonzero since g_tr != 0).
  double pr0; // Tetrad momentum p_r of the orbit.
  double pphi0; // Tetrad momentum p_phi of the orbit.
  double omega_orb; // Orbital angular frequency dphi/dt of the circular orbit.

  int Nr; // Cell count (configuration space: radial direction).
  int Nphi; // Cell count (configuration space: azimuthal direction).
  int Nvr; // Cell count (velocity space: radial momentum).
  int Nvphi; // Cell count (velocity space: azimuthal momentum).
  double r_min; // Inner radius of the domain.
  double r_max; // Outer radius of the domain.
  double pr_lo; // Lower boundary of the radial momentum domain.
  double pr_hi; // Upper boundary of the radial momentum domain.
  double pphi_lo; // Lower boundary of the azimuthal momentum domain.
  double pphi_hi; // Upper boundary of the azimuthal momentum domain.
  double sigma_r; // Blob width in r (half a cell).
  double sigma_phi; // Blob width in phi (half a cell).
  double sigma_pr; // Blob width in p_r (a quarter cell).
  double sigma_pphi; // Blob width in p_phi (a quarter cell).
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

struct geodesic_ctx
create_ctx(void)
{
  double mass = 1.0; // Neutral mass.
  double charge = 0.0; // Neutral charge.
  double mass_bh = 1.0; // Black hole mass M.
  double spin_bh = 0.8; // Black hole spin a = J/M.
  double r0 = 8.0; // Orbit radius.
  double phi0 = M_PI; // Initial azimuth of the blob.

  double v_c = sqrt(mass_bh / r0); // Orbital velocity parameter sqrt(M / r0).
  double E_c = (1.0 - (2.0 * v_c * v_c) + (spin_bh * v_c * v_c * v_c)) /
               sqrt(
                 1.0 - (3.0 * v_c * v_c) + (2.0 * spin_bh * v_c * v_c * v_c)
               ); // Specific energy of the prograde circular orbit (Boyer-Lindquist).
  double L_c =
    sqrt(mass_bh * r0) *
    (1.0 - (2.0 * spin_bh * v_c * v_c * v_c) + (spin_bh * spin_bh * v_c * v_c * v_c * v_c)) /
    sqrt(
      1.0 - (3.0 * v_c * v_c) + (2.0 * spin_bh * v_c * v_c * v_c)
    ); // Specific angular momentum of the circular orbit.
  double g_tt =
    -(1.0 - (2.0 * mass_bh / r0)); // Boyer-Lindquist metric component g_tt at r0 (equator).
  double g_tp = -2.0 * mass_bh * spin_bh / r0; // Boyer-Lindquist metric component g_tphi at r0.
  double g_pp =
    (r0 * r0) + (spin_bh * spin_bh) +
    (2.0 * mass_bh * spin_bh * spin_bh / r0); // Boyer-Lindquist metric component g_phiphi at r0.
  double det_g = (g_tt * g_pp) - (g_tp * g_tp); // Determinant of the (t, phi) metric block.
  double p_t_up = ((g_pp * (-E_c)) - (g_tp * L_c)) / det_g; // Contravariant p^t of the orbit.
  double p_phi_up = ((-g_tp * (-E_c)) + (g_tt * L_c)) / det_g; // Contravariant p^phi of the orbit.
  double h_rr = 1.0 + (2.0 * mass_bh / r0); // Kerr-Schild spatial metric h_rr at r0.
  double h_rp = -spin_bh * (1.0 + (2.0 * mass_bh / r0)); // Kerr-Schild spatial metric h_rphi at r0.
  double h_pp = (r0 * r0) + (spin_bh * spin_bh * (1.0 + (2.0 * mass_bh / r0))
                            ); // Kerr-Schild spatial metric h_phiphi at r0.
  double p_r_cov =
    ((2.0 * mass_bh / r0) * p_t_up) +
    (h_rp * p_phi_up); // Covariant Kerr-Schild p_r of the orbit (nonzero since g_tr != 0).
  double pr0 = p_r_cov / sqrt(h_rr); // Tetrad momentum p_r of the orbit.
  double pphi0 = (L_c - ((h_rp / sqrt(h_rr)) * pr0)) /
                 sqrt(h_pp - (h_rp * h_rp / h_rr)); // Tetrad momentum p_phi of the orbit.
  double omega_orb = sqrt(mass_bh) / ((r0 * sqrt(r0)) + (spin_bh * sqrt(mass_bh))
                                     ); // Orbital angular frequency dphi/dt of the circular orbit.

  int Nr = 12; // Cell count (configuration space: radial direction).
  int Nphi = 24; // Cell count (configuration space: azimuthal direction).
  int Nvr = 6; // Cell count (velocity space: radial momentum).
  int Nvphi = 6; // Cell count (velocity space: azimuthal momentum).
  double r_min = 6.0; // Inner radius of the domain.
  double r_max = 10.0; // Outer radius of the domain.
  double pr_lo = 0.1; // Lower boundary of the radial momentum domain.
  double pr_hi = 0.5; // Upper boundary of the radial momentum domain.
  double pphi_lo = 0.3; // Lower boundary of the azimuthal momentum domain.
  double pphi_hi = 0.6; // Upper boundary of the azimuthal momentum domain.
  double sigma_r = (r_max - r_min) / (2.0 * Nr); // Blob width in r (half a cell).
  double sigma_phi = 2.0 * M_PI / (2.0 * Nphi); // Blob width in phi (half a cell).
  double sigma_pr = (pr_hi - pr_lo) / (4.0 * Nvr); // Blob width in p_r (a quarter cell).
  double sigma_pphi = (pphi_hi - pphi_lo) / (4.0 * Nvphi); // Blob width in p_phi (a quarter cell).
  int poly_order = 2; // Polynomial order.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 25.0; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs =
    INT_MAX; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct geodesic_ctx ctx = {
    .mass = mass,
    .charge = charge,
    .mass_bh = mass_bh,
    .spin_bh = spin_bh,
    .r0 = r0,
    .phi0 = phi0,
    .v_c = v_c,
    .E_c = E_c,
    .L_c = L_c,
    .g_tt = g_tt,
    .g_tp = g_tp,
    .g_pp = g_pp,
    .det_g = det_g,
    .p_t_up = p_t_up,
    .p_phi_up = p_phi_up,
    .h_rr = h_rr,
    .h_rp = h_rp,
    .h_pp = h_pp,
    .p_r_cov = p_r_cov,
    .pr0 = pr0,
    .pphi0 = pphi0,
    .omega_orb = omega_orb,
    .Nr = Nr,
    .Nphi = Nphi,
    .Nvr = Nvr,
    .Nvphi = Nvphi,
    .r_min = r_min,
    .r_max = r_max,
    .pr_lo = pr_lo,
    .pr_hi = pr_hi,
    .pphi_lo = pphi_lo,
    .pphi_hi = pphi_hi,
    .sigma_r = sigma_r,
    .sigma_phi = sigma_phi,
    .sigma_pr = sigma_pr,
    .sigma_pphi = sigma_pphi,
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
  struct geodesic_ctx *app = ctx;
  double r = xn[0], phi = xn[1], p_r = xn[2], p_phi = xn[3];

  // Gaussian blob in phase space centred on the circular orbit (r0, phi0, pr0, pphi0).
  double dr = (r - app->r0) / app->sigma_r, dphi = (phi - app->phi0) / app->sigma_phi;
  double dpr = (p_r - app->pr0) / app->sigma_pr, dpphi = (p_phi - app->pphi0) / app->sigma_pphi;
  fout[0] = exp(-0.5 * ((dr * dr) + (dphi * dphi) + (dpr * dpr) + (dpphi * dpphi)));
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

  struct geodesic_ctx ctx = create_ctx(); // Context for initialization functions.

  int NR = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nr);
  int NPHI = APP_ARGS_CHOOSE(app_args.xcells[1], ctx.Nphi);
  int NVR = APP_ARGS_CHOOSE(app_args.vcells[0], ctx.Nvr);
  int NVPHI = APP_ARGS_CHOOSE(app_args.vcells[1], ctx.Nvphi);

  int nrank = 1; // Number of processors in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif

  int ccells[] = {NR, NPHI};
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
    .model_id = GKYL_MODEL_TRIAD_GR,
    .lower = {ctx.pr_lo, ctx.pphi_lo},
    .upper = {ctx.pr_hi, ctx.pphi_hi},
    .cells = {NVR, NVPHI},

    .num_init = 1,
    .projection[0] = {.proj_id = GKYL_PROJ_FUNC, .func = evalInit, .ctx_func = &ctx},

    .bcx = {.lower = {.type = GKYL_SPECIES_ABSORB}, .upper = {.type = GKYL_SPECIES_ABSORB}},

    .num_diag_moments = 2,
    .diag_moments = {GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1},
  };

  // Preset geometry.
  struct gkyl_vlasov_geom geom = {
    .use_preset_geom = true,
    .triad_preset_geom_type = GKYL_TRIAD_GR_KERR_SCHILD_RPHI,
    .mass_bh = ctx.mass_bh,
    .spin_bh = ctx.spin_bh,
  };

  // Vlasov app (no field).
  struct gkyl_vm app_inp = {
    .geom = geom,

    .cdim = 2,
    .vdim = 2,
    .lower = {ctx.r_min, 0.0},
    .upper = {ctx.r_max, 2.0 * M_PI},
    .cells = {NR, NPHI},

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
