// Sod shock tube of a neutral gas with implicit BGK collisions through the triad Vlasov solver on a
// cylinder in (r, z) with a single periodic axial cell and the preset cylindrical geometry, jump at
// r = 1 on 0.5 < r < 1.5, on cubic velocity maps (2x3v). The states (n, T) = (1, 1) and
// (0.125, 0.8) are at rest on either side of the jump; with collision frequency NU the solution is
// the Euler one with gamma = 5/3. The driver runs the tensor p = 1 hybrid basis
// (p = 2 in velocity space), the configuration used for kinetic neutrals in the gyrokinetic solver,
// and must reproduce the 1x3v driver. The reference is a quasi-1D finite-volume Euler solution with
// the metric determinant as area factor (reflecting walls), and the exact planar Riemann solution
// is quoted for contrast. Figures of merit (tensor p1 hybrid, 16 x 1 x 6 x 6 x 6 cells on cubic
// velocity maps, 52 steps): density L1 error 0.0373 against the quasi-1D reference (0.0364 against
// the planar exact solution); shock at 1.2031 (reference 1.1821) within one cell; cell-averaged
// density within 0.027 (relative L1 0.012) of the serendipity p2 1x3v driver at 16 radial cells;
// mass and energy conserved to round-off.

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

struct cylindrical_sodshock_ctx {
  double mass; // Neutral mass.
  double charge; // Neutral charge.
  double nl; // Left/inner number density.
  double Tl; // Left/inner temperature.
  double nr; // Right/outer number density.
  double Tr; // Right/outer temperature.
  double vt; // Thermal velocity.
  double nu; // Collision frequency.
  double midplane; // Radial location of the jump.

  int Nr; // Cell count (configuration space: radial direction).
  int Nz; // Cell count (configuration space: axial direction).
  int Nvr; // Cell count (velocity space: radial direction).
  int Nvz; // Cell count (velocity space: axial direction).
  int Nvtheta; // Cell count (velocity space: angular direction).
  double Lr; // Domain size (configuration space: radial direction).
  double Lz; // Domain size (configuration space: axial direction).
  double vr_max; // Domain boundary (velocity space: radial direction).
  double vz_max; // Domain boundary (velocity space: axial direction).
  double vtheta_max; // Domain boundary (velocity space: angular direction).
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

struct cylindrical_sodshock_ctx
create_ctx(void)
{
  double mass = 1.0; // Neutral mass.
  double charge = 0.0; // Neutral charge.
  double nl = 1.0; // Left/inner number density.
  double Tl = 1.0; // Left/inner temperature.
  double nr = 0.125; // Right/outer number density.
  double Tr = 0.8; // Right/outer temperature.
  double vt = 1.0; // Thermal velocity.
  double nu = 1000.0; // Collision frequency.
  double midplane = 1.0; // Radial location of the jump.

  int Nr = 16; // Cell count (configuration space: radial direction).
  int Nz = 1; // Cell count (configuration space: axial direction).
  int Nvr = 6; // Cell count (velocity space: radial direction).
  int Nvz = 6; // Cell count (velocity space: axial direction).
  int Nvtheta = 6; // Cell count (velocity space: angular direction).
  double Lr = 1.0; // Domain size (configuration space: radial direction).
  double Lz = 1.0; // Domain size (configuration space: axial direction).
  double vr_max = 8.0 * vt; // Domain boundary (velocity space: radial direction).
  double vz_max = 8.0 * vt; // Domain boundary (velocity space: axial direction).
  double vtheta_max = 8.0 * vt; // Domain boundary (velocity space: angular direction).
  int poly_order = 1; // Polynomial order.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 0.1; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs =
    INT_MAX; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct cylindrical_sodshock_ctx ctx = {
    .mass = mass,
    .charge = charge,
    .nl = nl,
    .Tl = Tl,
    .nr = nr,
    .Tr = Tr,
    .vt = vt,
    .nu = nu,
    .midplane = midplane,
    .Nr = Nr,
    .Nz = Nz,
    .Nvr = Nvr,
    .Nvz = Nvz,
    .Nvtheta = Nvtheta,
    .Lr = Lr,
    .Lz = Lz,
    .vr_max = vr_max,
    .vz_max = vz_max,
    .vtheta_max = vtheta_max,
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
evalDensityInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct cylindrical_sodshock_ctx *app = ctx;
  double r = xn[0];

  // Left and right states, times the metric determinant.
  double n = (r < app->midplane) ? app->nl : app->nr;
  fout[0] = (r)*n;
}

void
evalTempInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct cylindrical_sodshock_ctx *app = ctx;
  double r = xn[0];

  // Set isotropic temperature.
  fout[0] = (r < app->midplane) ? app->Tl : app->Tr;
}

void
evalVDriftInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct cylindrical_sodshock_ctx *app = ctx;
  // Set drift velocity in the orthonormal frame (the gas is at rest).
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
}

void
evalNu(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct cylindrical_sodshock_ctx *app = ctx;
  // Set collision frequency.
  fout[0] = app->nu;
}

struct vel_map_ctx {
  double vmax; // Physical velocity-space extent.
  double s; // Cubic stretching of the map: v = vc (1 + s (vc/vmax)^2) / (1 + s).
};

void
mapc2p_vel_dir(double t, const double *xc, double *xp, void *ctx)
{
  struct vel_map_ctx *vctx = ctx;
  double vc = xc[0], vm = vctx->vmax, s = vctx->s;
  xp[0] = vc * (1.0 + s * (vc / vm) * (vc / vm)) / (1.0 + s);
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

  struct cylindrical_sodshock_ctx ctx = create_ctx(); // Context for initialization functions.

  int NR = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nr);
  int NZ = APP_ARGS_CHOOSE(app_args.xcells[1], ctx.Nz);
  int NVR = APP_ARGS_CHOOSE(app_args.vcells[0], ctx.Nvr);
  int NVZ = APP_ARGS_CHOOSE(app_args.vcells[1], ctx.Nvz);
  int NVTHETA = APP_ARGS_CHOOSE(app_args.vcells[2], ctx.Nvtheta);

  int nrank = 1; // Number of processors in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif

  int ccells[] = {NR, NZ};
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

  // Cubic velocity maps: linear core at the origin, stretched tails.
  struct vel_map_ctx vmap_ctx_vr = {.vmax = ctx.vr_max, .s = 1.5};
  struct vel_map_ctx vmap_ctx_vz = {.vmax = ctx.vz_max, .s = 1.5};
  struct vel_map_ctx vmap_ctx_vtheta = {.vmax = ctx.vtheta_max, .s = 1.5};

  // Neutral species.
  struct gkyl_vlasov_kinetic_species neut = {
    .model_id = GKYL_MODEL_TRIAD,
    .lower = {-ctx.vr_max, -ctx.vz_max, -ctx.vtheta_max},
    .upper = {ctx.vr_max, ctx.vz_max, ctx.vtheta_max},
    .cells = {NVR, NVZ, NVTHETA},

    .mapc2p_vel =
      {{.mapc2p_vel_func = mapc2p_vel_dir, .mapc2p_vel_ctx = &vmap_ctx_vr},
       {.mapc2p_vel_func = mapc2p_vel_dir, .mapc2p_vel_ctx = &vmap_ctx_vz},
       {.mapc2p_vel_func = mapc2p_vel_dir, .mapc2p_vel_ctx = &vmap_ctx_vtheta}},

    .num_init = 1,
    .projection[0] =
      {
        .proj_id = GKYL_PROJ_VLASOV_LTE,
        .density = evalDensityInit,
        .ctx_density = &ctx,
        .temp = evalTempInit,
        .ctx_temp = &ctx,
        .V_drift = evalVDriftInit,
        .ctx_V_drift = &ctx,
        .correct_all_moms = true,
        .iter_eps = 0.0,
        .max_iter = 0,
        .use_last_converged = false,
      },

    .collisions =
      {
        .collision_id = GKYL_BGK_COLLISIONS,
        .self_nu = evalNu,
        .self_nu_ctx = &ctx,
        .is_implicit = true,
      },

    .correct =
      {.correct_all_moms = true, .iter_eps = 1.0e-12, .max_iter = 100, .use_last_converged = false},

    .bcx = {.lower = {.type = GKYL_SPECIES_REFLECT}, .upper = {.type = GKYL_SPECIES_REFLECT}},

    .num_diag_moments = 4,
    .diag_moments = {GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_LTE, GKYL_F_MOMENT_ENERGY},
  };

  // Preset geometry.
  struct gkyl_vlasov_geom geom = {
    .use_preset_geom = true,
    .triad_preset_geom_type = GKYL_TRIAD_CYLINDRICAL_RZ,
    .mass_bh = 0.0,
    .spin_bh = 0.0,
  };

  // Vlasov app (no field).
  struct gkyl_vm app_inp = {
    .geom = geom,

    .cdim = 2,
    .vdim = 3,
    .lower = {0.5, 0.0},
    .upper = {0.5 + ctx.Lr, ctx.Lz},
    .cells = {NR, NZ},

    .poly_order = ctx.poly_order,
    .basis_type = GKYL_BASIS_MODAL_TENSOR,
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
