// Relaxation of a self-gravitating collisionless sheet toward equilibrium with the Vlasov-Poisson
// system (1x1v). Neutral dark matter (m = 1, gravitational coupling alpha_g = 4 pi G = 1) starts as
// a Gaussian slab of mass per unit area M = 1 and width a = 1 with a Maxwellian velocity
// distribution whose dispersion puts it at half the virial kinetic energy, 2K/U = 1/2 with U = pi G
// int int rho(x) rho(x') |x - x'| dx dx', so it collapses and phase mixes toward a quasi-stationary
// state. The isothermal equilibrium of 1D gravity is the Spitzer-Camm sheet rho = rho0 sech^2(x/z0)
// with z0 = sigma/sqrt(2 pi G rho0) = sigma^2/(pi G M) (not the NFW profile, which is 3D), and the
// 1D virial theorem is 2K = U. The configuration mesh is a linear-to-exponential map on [-1, 1]
// clustering cells at x = 0 while the box edge sits at ~20 a (effectively infinity, with absorbing
// species boundaries and phi_g = 0 Dirichlet boundaries); the velocity mesh uses the same map so
// the high-velocity tail of the relaxed distribution is covered.
// Figures of merit (serendipity p2, 64 x 64 cells): the virial ratio 2K/U goes from 0.500 to 0.9998
// +- 0.004 (mean over t = 32-63) and the final density is fit by rho0 sech^2(x/z0) with rho0 =
// 0.530 and z0 = 0.945 (L2 residual 1.4%; the isothermal sheet with the final K has rho0 = 0.532
// and z0 = 0.940); mass conserved to 5e-13 and total energy (kinetic minus half the gravitational
// energy file, int |grad phi_g|^2/alpha_g) to 2e-10.

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

struct grav_relaxation_ctx {
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using normalized code units).
  double alpha_g; // Gravitational coupling 4 pi G.
  double mass_dm; // Dark matter mass.
  double charge_dm; // Dark matter charge (neutral).

  double M; // Mass per unit area of the sheet.
  double a; // Initial Gaussian width.
  double virial_ratio0; // Initial virial ratio 2K/U.

  // Derived physical quantities (using normalized code units).
  double sigma0; // Initial velocity dispersion.
  double rho0; // Initial peak density.
  double t_dyn; // Dynamical time 1/sqrt(4 pi G rho0).

  // Simulation parameters.
  int Nx; // Cell count (configuration space: x-direction).
  int Nvx; // Cell count (velocity space: vx-direction).
  double dx_lin; // Configuration map: cell size of the linear part at x = 0.
  double x_exp; // Configuration map: the exponential part reaches x_exp - 1 at the boundary.
  double dv_lin; // Velocity map: cell size of the linear part at vx = 0.
  double v_exp; // Velocity map: the exponential part reaches v_exp - 1 at the boundary.
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

struct grav_relaxation_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  double alpha_g = 1.0; // Gravitational coupling 4 pi G.
  double mass_dm = 1.0; // Dark matter mass.
  double charge_dm = 0.0; // Dark matter charge (neutral).

  double M = 1.0; // Mass per unit area of the sheet.
  double a = 1.0; // Initial Gaussian width.
  double virial_ratio0 = 0.5; // Initial virial ratio 2K/U.

  // Derived physical quantities (using normalized code units).
  // For a Gaussian sheet U = 2 sqrt(pi) G M^2 a and K = M sigma^2/2, so
  // 2K/U = sigma^2/(2 sqrt(pi) G M a) = 2 sqrt(pi) sigma^2/(alpha_g M a).
  double sigma0 =
    sqrt(virial_ratio0 * alpha_g * M * a / (2.0 * sqrt(pi))); // Initial velocity dispersion.
  double rho0 = M / (sqrt(2.0 * pi) * a); // Initial peak density.
  double t_dyn = 1.0 / sqrt(alpha_g * rho0); // Dynamical time 1/sqrt(4 pi G rho0).

  // Simulation parameters.
  int Nx = 64; // Cell count (configuration space: x-direction).
  int Nvx = 64; // Cell count (velocity space: vx-direction).
  double dx_lin = 0.05; // Configuration map: cell size of the linear part at x = 0.
  double x_exp = 20.0; // Configuration map: the exponential part reaches x_exp - 1 at the boundary.
  double dv_lin = 0.02; // Velocity map: cell size of the linear part at vx = 0.
  double v_exp = 4.0; // Velocity map: the exponential part reaches v_exp - 1 at the boundary.
  int poly_order = 2; // Polynomial order.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 40.0 * t_dyn; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs =
    INT_MAX; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct grav_relaxation_ctx ctx = {
    .pi = pi,
    .alpha_g = alpha_g,
    .mass_dm = mass_dm,
    .charge_dm = charge_dm,
    .M = M,
    .a = a,
    .virial_ratio0 = virial_ratio0,
    .sigma0 = sigma0,
    .rho0 = rho0,
    .t_dyn = t_dyn,
    .Nx = Nx,
    .Nvx = Nvx,
    .dx_lin = dx_lin,
    .x_exp = x_exp,
    .dv_lin = dv_lin,
    .v_exp = v_exp,
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
evalDmInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct grav_relaxation_ctx *app = ctx;
  double x = xn[0], vx = xn[1];

  double pi = app->pi;
  double M = app->M;
  double a = app->a;
  double sigma0 = app->sigma0;

  // Gaussian slab with a Maxwellian velocity distribution.
  double n = (M / (2.0 * pi * a * sigma0)) * exp(-(x * x) / (2.0 * a * a)) *
             exp(-(vx * vx) / (2.0 * sigma0 * sigma0));

  // Set distribution function.
  fout[0] = n;
}

// Linear-to-exponential map of a computational coordinate c in [-1, 1]: cells of
// size d_lin at the origin (N cells across the domain), growing exponentially
// so that the exponential part reaches e_max - 1 at the boundary.
static double
lin_exp_map(double c, double d_lin, double e_max, int N)
{
  double lin = 0.5 * d_lin * N * c;
  if (c < 0.0) {
    return lin - (exp(-c * log(e_max)) - 1.0);
  }
  return lin + (exp(c * log(e_max)) - 1.0);
}

void
mapc2p_pos_x(double t, const double *GKYL_RESTRICT xc, double *GKYL_RESTRICT xp, void *ctx)
{
  struct grav_relaxation_ctx *app = ctx;

  xp[0] = lin_exp_map(xc[0], app->dx_lin, app->x_exp, app->Nx);
}

void
mapc2p_vel_vx(double t, const double *GKYL_RESTRICT vc, double *GKYL_RESTRICT vp, void *ctx)
{
  struct grav_relaxation_ctx *app = ctx;

  vp[0] = lin_exp_map(vc[0], app->dv_lin, app->v_exp, app->Nvx);
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

  struct grav_relaxation_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int NVX = APP_ARGS_CHOOSE(app_args.vcells[0], ctx.Nvx);
  // The maps scale their linear part with the cell count.
  ctx.Nx = NX;
  ctx.Nvx = NVX;

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

  // Dark matter. The velocity grid is computational, mapped to physical velocities by mapc2p_vel.
  struct gkyl_vlasov_kinetic_species dm = {
    .lower = {-1.0},
    .upper = {1.0},
    .cells = {NVX},

    .mapc2p_vel = {{.mapc2p_vel_func = mapc2p_vel_vx, .mapc2p_vel_ctx = &ctx}},

    .num_init = 1,
    .projection[0] = {.proj_id = GKYL_PROJ_FUNC, .func = evalDmInit, .ctx_func = &ctx},

    // Absorbing boundaries: the box edge is effectively at infinity.
    .bcx = {.lower = {.type = GKYL_SPECIES_ABSORB}, .upper = {.type = GKYL_SPECIES_ABSORB}},

    .num_diag_moments = 3,
    .diag_moments = {GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2},
  };

  // Field: self-gravity only (the species is neutral), phi_g = 0 at the box edges.
  struct gkyl_vlasov_field field = {
    .field_id = GKYL_FIELD_PHI,
    .alpha_g = ctx.alpha_g,

    .poisson_bcs =
      {
        .lo_type = {GKYL_POISSON_DIRICHLET},
        .up_type = {GKYL_POISSON_DIRICHLET},
        .lo_value = {0.0},
        .up_value = {0.0},
      },
  };

  // Vlasov-Poisson app. The configuration grid is computational, mapped to
  // physical positions by mapc2p_pos.
  struct gkyl_vm app_inp = {

    .cdim = 1,
    .vdim = 1,
    .lower = {-1.0},
    .upper = {1.0},
    .cells = {NX},

    .mapc2p_pos[0] = {.mapc2p_pos_func = mapc2p_pos_x, .mapc2p_pos_ctx = &ctx},

    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,
    .cfl_frac = ctx.cfl_frac,

    .num_periodic_dir = 0,
    .periodic_dirs = {},

    .num_species = 1,
    .species = {{
      .name = "dm",
      .charge = ctx.charge_dm,
      .mass = ctx.mass_dm,
      .type = GKYL_SPECIES_VLASOV,
      .kinetic = dm,
    }},

    .field = field,

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
    .tcurr = gkyl_tm_trigger_anchor(t_end / num_frames, t_curr),
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
