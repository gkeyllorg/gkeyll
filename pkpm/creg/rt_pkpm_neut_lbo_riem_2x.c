// Two-dimensional Riemann problem for a neutral gas with LBO collisions with the PKPM system of
// equations (2x, normalized parallel velocity w = v_par/xi).
// Four constant states meeting at (0.5, 0.5) in the unit square (configuration 3 of Lax & Liu),
// with (density, x-velocity, y-velocity, pressure) = (1.5, 0, 0, 1.5) upper right, (0.5323, 1.206,
// 0, 0.3) upper left, (0.138, 1.206, 1.206, 0.029) lower left and (0.5323, 0, 1.206, 0.3) lower
// right, to t = 0.3. A static uniform magnetic field along the diagonal x = y sets the field
// direction b, at 45 degrees to the grid.
// Moderate collisions advanced explicitly, nu = 100.0: the mean free path of the hottest gas is
// 0.2 cells. This is the collision frequency at which the diffusion step of the LBO operator is
// about the transport step; rt_pkpm_neut_bgk_riem_2x has the same collision frequency and grid.
// Figures of merit (20 x 20 x 8 cells) at t = 0.3, with nu dt = 0.012 on average: against the
// Euler solution of the moments app (400 x 400 cells) averaged over the cells of this grid the L1
// differences of the cell averages are 3.2% (density) and 3.8% (pressure); the density stays
// between 0.1379 and 1.5726, and the solution is symmetric about the diagonal to round-off.

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

struct riem_ctx {
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using normalized code units).
  double epsilon0; // Permittivity of free space.
  double mu0; // Permeability of free space.
  double mass; // Neutral mass.
  double charge; // Neutral charge.

  double rho_ul; // Upper left mass density.
  double u_ul; // Upper left x-velocity.
  double v_ul; // Upper left y-velocity.
  double p_ul; // Upper left pressure.

  double rho_ur; // Upper right mass density.
  double u_ur; // Upper right x-velocity.
  double v_ur; // Upper right y-velocity.
  double p_ur; // Upper right pressure.

  double rho_ll; // Lower left mass density.
  double u_ll; // Lower left x-velocity.
  double v_ll; // Lower left y-velocity.
  double p_ll; // Lower left pressure.

  double rho_lr; // Lower right mass density.
  double u_lr; // Lower right x-velocity.
  double v_lr; // Lower right y-velocity.
  double p_lr; // Lower right pressure.

  double loc; // Corner of the four states (both x and y coordinates).
  double nu; // Collision frequency.
  double B0; // Magnetic field along each of x and y (sets the field direction b).

  // Simulation parameters.
  int Nx; // Cell count (configuration space: x-direction).
  int Ny; // Cell count (configuration space: y-direction).
  int Nvx; // Cell count (velocity space: normalized parallel velocity w).
  double Lx; // Domain size (configuration space: x-direction).
  double Ly; // Domain size (configuration space: y-direction).
  double vx_max; // Domain boundary (velocity space: w = v_par/xi).
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

struct riem_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  double epsilon0 = 1.0; // Permittivity of free space.
  double mu0 = 1.0; // Permeability of free space.
  double mass = 1.0; // Neutral mass.
  double charge = 0.0; // Neutral charge.

  double rho_ul = 0.5323; // Upper left mass density.
  double u_ul = 1.206; // Upper left x-velocity.
  double v_ul = 0.0; // Upper left y-velocity.
  double p_ul = 0.3; // Upper left pressure.

  double rho_ur = 1.5; // Upper right mass density.
  double u_ur = 0.0; // Upper right x-velocity.
  double v_ur = 0.0; // Upper right y-velocity.
  double p_ur = 1.5; // Upper right pressure.

  double rho_ll = 0.138; // Lower left mass density.
  double u_ll = 1.206; // Lower left x-velocity.
  double v_ll = 1.206; // Lower left y-velocity.
  double p_ll = 0.029; // Lower left pressure.

  double rho_lr = 0.5323; // Lower right mass density.
  double u_lr = 0.0; // Lower right x-velocity.
  double v_lr = 1.206; // Lower right y-velocity.
  double p_lr = 0.3; // Lower right pressure.

  double loc = 0.5; // Corner of the four states (both x and y coordinates).
  double nu = 100.0; // Collision frequency.
  double B0 = 1.0 / sqrt(2.0); // Magnetic field along each of x and y (sets the field direction b).

  // Simulation parameters.
  int Nx = 20; // Cell count (configuration space: x-direction).
  int Ny = 20; // Cell count (configuration space: y-direction).
  int Nvx = 8; // Cell count (velocity space: normalized parallel velocity w).
  double Lx = 1.0; // Domain size (configuration space: x-direction).
  double Ly = 1.0; // Domain size (configuration space: y-direction).
  double vx_max = 6.0; // Domain boundary (velocity space: w = v_par/xi).
  int poly_order = 1; // Polynomial order.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 0.3; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs =
    INT_MAX; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct riem_ctx ctx = {
    .pi = pi,
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .mass = mass,
    .charge = charge,
    .rho_ul = rho_ul,
    .u_ul = u_ul,
    .v_ul = v_ul,
    .p_ul = p_ul,
    .rho_ur = rho_ur,
    .u_ur = u_ur,
    .v_ur = v_ur,
    .p_ur = p_ur,
    .rho_ll = rho_ll,
    .u_ll = u_ll,
    .v_ll = v_ll,
    .p_ll = p_ll,
    .rho_lr = rho_lr,
    .u_lr = u_lr,
    .v_lr = v_lr,
    .p_lr = p_lr,
    .loc = loc,
    .nu = nu,
    .B0 = B0,
    .Nx = Nx,
    .Ny = Ny,
    .Nvx = Nvx,
    .Lx = Lx,
    .Ly = Ly,
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
evalNeutDistInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct riem_ctx *app = ctx;
  double w = xn[2];

  double pi = app->pi;

  double F = exp(-0.5 * w * w) / sqrt(2.0 * pi); // Normalized Maxwellian in w = v_par/xi.

  // Set normalized neutral distribution functions [F, G] (G = (T_perp/T) F = F).
  fout[0] = F;
  fout[1] = F;
}

void
evalNeutFluidInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct riem_ctx *app = ctx;
  double x = xn[0], y = xn[1];

  double loc = app->loc;

  double rho = 0.0;
  double vx = 0.0;
  double vy = 0.0;
  double pres = 0.0;

  if (y > loc) {
    if (x < loc) {
      rho = app->rho_ul; // Mass density (upper left).
      vx = app->u_ul; // x-velocity (upper left).
      vy = app->v_ul; // y-velocity (upper left).
      pres = app->p_ul; // Pressure (upper left).
    } else {
      rho = app->rho_ur; // Mass density (upper right).
      vx = app->u_ur; // x-velocity (upper right).
      vy = app->v_ur; // y-velocity (upper right).
      pres = app->p_ur; // Pressure (upper right).
    }
  } else {
    if (x < loc) {
      rho = app->rho_ll; // Mass density (lower left).
      vx = app->u_ll; // x-velocity (lower left).
      vy = app->v_ll; // y-velocity (lower left).
      pres = app->p_ll; // Pressure (lower left).
    } else {
      rho = app->rho_lr; // Mass density (lower right).
      vx = app->u_lr; // x-velocity (lower right).
      vy = app->v_lr; // y-velocity (lower right).
      pres = app->p_lr; // Pressure (lower right).
    }
  }
  double E = 0.5 * rho * (vx * vx + vy * vy) + 1.5 * pres; // Total energy density.

  // Set neutral fluid state [rho, rho ux, rho uy, rho uz, E].
  fout[0] = rho;
  fout[1] = rho * vx;
  fout[2] = rho * vy;
  fout[3] = 0.0;
  fout[4] = E;
}

void
evalNu(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct riem_ctx *app = ctx;

  double nu = app->nu;

  // Set collision frequency.
  fout[0] = nu;
}

void
evalFieldInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct riem_ctx *app = ctx;

  double B0 = app->B0;

  // Set electric field, magnetic field and correction potentials.
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = B0;
  fout[4] = B0;
  fout[5] = 0.0;
  fout[6] = 0.0;
  fout[7] = 0.0;
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
write_data(struct gkyl_tm_trigger *iot, gkyl_pkpm_app *app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_write) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_pkpm_app_write(app, t_curr, frame);
    gkyl_pkpm_app_write_field_energy(app);
    gkyl_pkpm_app_write_integrated_mom(app);
    gkyl_pkpm_app_write_integrated_L2_f(app);
  }
}

void
calc_field_energy(struct gkyl_tm_trigger *fet, gkyl_pkpm_app *app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(fet, t_curr) || force_calc) {
    gkyl_pkpm_app_calc_field_energy(app, t_curr);
  }
}

void
calc_integrated_mom(struct gkyl_tm_trigger *imt, gkyl_pkpm_app *app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(imt, t_curr) || force_calc) {
    gkyl_pkpm_app_calc_integrated_mom(app, t_curr);
  }
}

void
calc_integrated_L2_f(struct gkyl_tm_trigger *l2t, gkyl_pkpm_app *app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(l2t, t_curr) || force_calc) {
    gkyl_pkpm_app_calc_integrated_L2_f(app, t_curr);
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

  struct riem_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], ctx.Ny);
  int NVX = APP_ARGS_CHOOSE(app_args.vcells[0], ctx.Nvx);

  int nrank = 1; // Number of processors in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif

  int ccells[] = {NX, NY};
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

  // Neutrals.
  struct gkyl_pkpm_species neut = {
    .name = "neut",
    .charge = ctx.charge,
    .mass = ctx.mass,
    .lower = {-ctx.vx_max},
    .upper = {ctx.vx_max},
    .cells = {NVX},

    .init_dist = evalNeutDistInit,
    .ctx_dist = &ctx,
    .init_fluid = evalNeutFluidInit,
    .ctx_fluid = &ctx,

    .collisions = {.collision_id = GKYL_LBO_COLLISIONS, .self_nu = evalNu, .ctx = &ctx},
  };

  // Field.
  struct gkyl_pkpm_field field = {
    .epsilon0 = ctx.epsilon0,
    .mu0 = ctx.mu0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,

    .init = evalFieldInit,
    .ctx = &ctx,

    .is_static = true,
  };

  // PKPM app.
  struct gkyl_pkpm app_inp = {

    .cdim = 2,
    .vdim = 1,
    .lower = {0.0, 0.0},
    .upper = {ctx.Lx, ctx.Ly},
    .cells = {NX, NY},

    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,
    .cfl_frac = ctx.cfl_frac,

    .num_periodic_dir = 0,
    .periodic_dirs = {},

    .num_species = 1,
    .species = {neut},

    .field = field,

    .parallelism =
      {.use_gpu = app_args.use_gpu, .cuts = {app_args.cuts[0], app_args.cuts[1]}, .comm = comm},
  };
  // Create app object.
  // Set app output name from the executable name (argv[0]).
  snprintf(app_inp.name, sizeof(app_inp.name), "%s", app_args.app_name);
  gkyl_pkpm_app *app = gkyl_pkpm_app_new(&app_inp);

  // Initial and final simulation times.
  double t_curr = 0.0, t_end = ctx.t_end;

  // Initialize simulation.
  int frame_curr = 0;
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status =
      gkyl_pkpm_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_pkpm_app_cout(
        app, stderr, "*** Failed to read restart file! (%s)\n",
        gkyl_array_rio_status_msg(status.io_status)
      );
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_pkpm_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_pkpm_app_cout(app, stdout, " at time = %g\n", t_curr);
  } else {
    gkyl_pkpm_app_apply_ic(app, t_curr);
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
    gkyl_pkpm_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = gkyl_pkpm_update(app, dt);
    gkyl_pkpm_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

    if (!status.success) {
      gkyl_pkpm_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
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

      gkyl_pkpm_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_pkpm_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_pkpm_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_pkpm_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_pkpm_app_cout(
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
  gkyl_pkpm_app_stat_write(app);

  struct gkyl_pkpm_stat stat = gkyl_pkpm_app_stat(app);

  gkyl_pkpm_app_cout(app, stdout, "\n");
  gkyl_pkpm_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_pkpm_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_pkpm_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_pkpm_app_cout(
      app, stdout, "  Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]
    );
    gkyl_pkpm_app_cout(
      app, stdout, "  Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]
    );
  }
  gkyl_pkpm_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_pkpm_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_pkpm_app_cout(
    app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm
  );
  gkyl_pkpm_app_cout(
    app, stdout, "Fluid species RHS calc took %g secs\n", stat.fluid_species_rhs_tm
  );
  gkyl_pkpm_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_pkpm_app_cout(app, stdout, "Species PKPM vars took %g secs\n", stat.species_pkpm_vars_tm);
  gkyl_pkpm_app_cout(
    app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm
  );
  gkyl_pkpm_app_cout(app, stdout, "EM variables (bvar) calc took %g secs\n", stat.field_em_vars_tm);
  gkyl_pkpm_app_cout(
    app, stdout, "Current evaluation and accumulate took %g secs\n", stat.current_tm
  );
  gkyl_pkpm_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

  gkyl_pkpm_app_cout(app, stdout, "Number of write calls %ld\n", stat.n_io);
  gkyl_pkpm_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

freeresources:
  // Free resources after simulation completion.
  gkyl_comm_release(comm);
  gkyl_pkpm_app_release(app);

mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif

  return 0;
}
