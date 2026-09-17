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

struct sphere_khi_ctx
{
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using normalized code units).
  double mass; // Neutral mass.
  double charge; // Neutral charge.

  double nl; // Left/inner number density.
  double Pl; // Left/inner pressure.
  double V_theta_drift_l; // Left/inner drift velocity (polar angular direction).
  double V_phi_drift_l; // Left/inner drift velocity (azimuthal angular direction).

  double nr; // Right/outer number density.
  double Pr; // Right/outer pressure.
  double V_theta_drift_r; // Right/outer drift velocity (polar angular direction).
  double V_phi_drift_r; // Right/outer drift velocity (azimuthal angular direction).

  double vt; // Thermal velocity.
  double nu; // Collision frequency.

  // Simulation parameters.
  int Ntheta; // Cell count (configuration space: polar angular direction).
  int Nphi; // Cell count (configuration space: azimuthal angular direction).
  int Nvtheta; // Cell count (velocity space: polar angular direction).
  int Nvphi; // Cell count (velocity space: azimuthal angular direction).
  double Ltheta; // Domain size (configuration space: polar angular direction).
  double Lphi; // Domain size (configuration space: azimuthal angular direction).
  double vtheta_max; // Domain boundary (velocity space: polar angular direction).
  double vphi_max; // Domain boundary (velocity space: azimuthal angular direction).
  int poly_order; // Polynomial order.
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  int field_energy_calcs; // Number of times to calculate field energy.
  int integrated_mom_calcs; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.

  double R; // Radius of the sphere.
  double midplane; // Polar angular midplane location.
  double theta_loc; // Polar angular boundary location designating jump in quantities.
};

struct sphere_khi_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  double mass = 1.0; // Neutral mass.
  double charge = 0.0; // Neutral charge.

  double nl = 2.0; // Left/inner number density.
  double Pl = 2.5; // Left/inner pressure.
  double V_theta_drift_l = 0.0; // Left/inner drift velocity (polar angular direction).
  double V_phi_drift_l = -0.5; // Left/inner drift velocity (azimuthal angular direction).

  double nr = 1.0; // Right/outer number density.
  double Pr = 2.5; // Right/outer pressure.
  double V_theta_drift_r = 0.0; // Right/outer drift velocity (polar angular direction).
  double V_phi_drift_r = 0.5; // Right/outer drift velocity (azimuthal angular direction).

  double vt = 1.0; // Thermal velocity.
  double nu = 15000.0; // Collision frequency.

  // Simulation parameters.
  int Ntheta = 32; // Cell count (configuration space: polar angular direction).
  int Nphi = 32; // Cell count (configuration space: azimuthal angular direction).
  int Nvtheta = 8; // Cell count (velocity space: polar angular direction).
  int Nvphi = 8; // Cell count (velocity space: azimuthal angular direction).
  double Ltheta = pi / 2.0; // Domain size (configuration space: polar angular direction).
  double Lphi = 2.0 * pi; // Domain size (configuration space: azimuthal angular direction).
  double vtheta_max = 8.0 * vt; // Domain boundary (velocity space: polar angular direction).
  double vphi_max = 8.0 * vt; // Domain boundary (velocity space: azimuthal angular direction).
  int poly_order = 2; // Polynomial order.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 0.01; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs = INT_MAX; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  double R = 1.0; // Radius of the sphere.
  double midplane = pi / 2.0; // Polar angular midplane location.
  double theta_loc = pi / 8.0; // Polar angular boundary location designating jump in quantities.

  struct sphere_khi_ctx ctx = {
    .pi = pi,
    .mass = mass,
    .charge = charge,
    .nl = nl,
    .Pl = Pl,
    .V_theta_drift_l = V_theta_drift_l,
    .V_phi_drift_l = V_phi_drift_l,
    .nr = nr,
    .Pr = Pr,
    .V_theta_drift_r = V_theta_drift_r,
    .V_phi_drift_r = V_phi_drift_r,
    .vt = vt,
    .nu = nu,
    .Ntheta = Ntheta,
    .Nphi = Nphi,
    .Nvtheta = Nvtheta,
    .Nvphi = Nvphi,
    .Ltheta = Ltheta,
    .Lphi = Lphi,
    .vtheta_max = vtheta_max,
    .vphi_max = vphi_max,
    .poly_order = poly_order,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .num_frames = num_frames,
    .field_energy_calcs = field_energy_calcs,
    .integrated_mom_calcs = integrated_mom_calcs,
    .integrated_L2_f_calcs = integrated_L2_f_calcs,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
    .R = R,
    .midplane = midplane,
    .theta_loc = theta_loc,
  };

  return ctx;
}

void
evalDensityInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sphere_khi_ctx *app = ctx;
  double theta = xn[0];

  double nl = app->nl;
  double nr = app->nr;

  double R = app->R;
  double midplane = app->midplane;
  double theta_loc = app->theta_loc;

  double n = 0.0;

  if (fabs(theta - midplane) < theta_loc) {
    n = nl; // Total number density (left/inner).
  }
  else {
    n = nr; // Total number density (right/outer).
  }

  double metric_det = (R * R) * sin(theta);

  // Set total number density.
  fout[0] = metric_det * n;
}

void
evalTempInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sphere_khi_ctx *app = ctx;
  double theta = xn[0];

  double nl = app->nl;
  double nr = app->nr;
  double Pl = app->Pl;
  double Pr = app->Pr;

  double midplane = app->midplane;
  double theta_loc = app->theta_loc;

  double T = 0.0;

  if (fabs(theta - midplane) < theta_loc) {
    T = Pl/nl; // Isotropic temperature (left/inner).
  }
  else {
    T = Pr/nr; // Isotropic temperature (right/outer).
  }

  // Set isotropic temperature.
  fout[0] = T;
}

void
evalVDriftInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sphere_khi_ctx *app = ctx;
  double theta = xn[0];
  double phi = xn[1];
  double pi = app->pi;

  double V_theta_drift_l = app->V_theta_drift_l;
  double V_phi_drift_l = app->V_phi_drift_l;

  double V_theta_drift_r = app->V_theta_drift_r;
  double V_phi_drift_r = app->V_phi_drift_r;

  double midplane = app->midplane;
  double theta_loc = app->theta_loc;

  double V_theta_drift = 0.0;
  double V_phi_drift = 0.0;

  if (fabs(theta - midplane) < theta_loc) {
    V_phi_drift = V_phi_drift_l; // Azimuthal angular drift velocity (left/inner).
  }
  else {
    V_phi_drift = V_phi_drift_r; // Azimuthal angular drift velocity (right/outer).
  }

  pcg64_random_t rng = gkyl_pcg64_init(0); // Random number generator;

  // Initalize noise
  double alpha = 1.0e-2;
  double k_theta = 2.0 * pi;
  double k_phi = 2.0 * pi/app->Lphi;

  for (int i = 0; i < 16; i++) {
    for (int j = 0; j < 16; j++) {
      V_theta_drift += alpha * gkyl_pcg64_rand_double(&rng) * sin(i * k_theta * theta + j * k_phi * phi + 2.0 * pi * gkyl_pcg64_rand_double(&rng));
      V_phi_drift += alpha * gkyl_pcg64_rand_double(&rng) * sin(i * k_theta * theta + j * k_phi * phi + 2.0 * pi * gkyl_pcg64_rand_double(&rng));
    }
  }

  // Set total drift velocity.
  fout[0] = V_theta_drift; fout[1] = V_phi_drift;
}

void
evalNu(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sphere_khi_ctx *app = ctx;

  double nu = app->nu;

  // Set collision frequency.
  fout[0] = nu;
}

void
evalHamiltonian(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sphere_khi_ctx *app = ctx;
  double q_theta = xn[0], p_theta_dot = xn[2], p_phi_dot = xn[3];
  
  double R = app->R;

  double inv_metric_theta_theta = 1.0 / (R * R);
  double inv_metric_theta_phi = 0.0;
  double inv_metric_phi_phi = 1.0 / ((R * sin(q_theta)) * (R * sin(q_theta)));

  double hamiltonian = (0.5 * inv_metric_theta_theta * p_theta_dot * p_theta_dot) + (0.5 * (2.0 * inv_metric_theta_phi * p_theta_dot * p_phi_dot)) +
    (0.5 * inv_metric_phi_phi * p_phi_dot * p_phi_dot); // Canonical Hamiltonian.
  
  // Set canonical Hamiltonian.
  fout[0] = hamiltonian;
}

void
evalInvMetric(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sphere_khi_ctx *app = ctx;
  double q_theta = xn[0];

  double R = app->R;

  double inv_metric_theta_theta = 1.0 / (R * R); // Inverse metric tensor (polar-polar component).
  double inv_metric_theta_phi = 0.0; // Inverse metric tensor (polar-azimuthal component).
  double inv_metric_phi_phi = 1.0 / ((R * sin(q_theta)) * (R * sin(q_theta))); // Inverse metric tensor (azimuthal-azimuthal component).
  
  // Set inverse metric tensor.
  fout[0] = inv_metric_theta_theta; fout[1] = inv_metric_theta_phi; fout[2] = inv_metric_phi_phi;
}

void
evalMetric(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sphere_khi_ctx *app = ctx;
  double q_theta = xn[0];

  double R = app->R;

  double metric_theta_theta = R * R; // Metric tensor (polar-polar component).
  double metric_theta_phi = 0.0; // Metric tensor (polar-azimuthal component).
  double metric_phi_phi = (R * sin(q_theta)) * (R * sin(q_theta)); // Metric tensor (azimuthal-azimuthal component).
  
  // Set metric tensor.
  fout[0] = metric_theta_theta; fout[1] = metric_theta_phi; fout[2] = metric_phi_phi;
}

void
evalMetricDet(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sphere_khi_ctx *app = ctx;
  double q_theta = xn[0];

  double R = app->R;

  double metric_det = (R * R) * sin(q_theta); // Metric tensor determinant.
  
  // Set metric tensor determinant.
  fout[0] = metric_det;
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

    gkyl_vlasov_app_calc_mom(app);
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

  struct sphere_khi_ctx ctx = create_ctx(); // Context for initialization functions.

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

  int ccells[] = { NTHETA, NPHI };
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

  // Neutral species.
  struct gkyl_vlasov_species neut = {
    .name = "neut",
    .model_id = GKYL_MODEL_CANONICAL_PB,
    .charge = ctx.charge, .mass = ctx.mass,
    .lower = { -ctx.vtheta_max, -ctx.vphi_max },
    .upper = { ctx.vtheta_max, ctx.vphi_max },
    .cells = { NVTHETA, NVPHI },

    .hamil = evalHamiltonian,
    .hamil_ctx = &ctx,
    .h_ij = evalMetric,
    .h_ij_ctx = &ctx,
    .h_ij_inv = evalInvMetric,
    .h_ij_inv_ctx = &ctx,
    .det_h = evalMetricDet,
    .det_h_ctx = &ctx,
    .output_f_lte = true,

    .num_init = 1, 
    .projection[0] = {
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
    .collisions =  {
      .collision_id = GKYL_BGK_COLLISIONS,
      .self_nu = evalNu,
      .ctx = &ctx,
      .has_implicit_coll_scheme = true,
      .correct_all_moms = true,
      .iter_eps = 0.0,
      .max_iter = 0,
      .use_last_converged = false,
    },

    .bcx = {
      .lower = { .type = GKYL_SPECIES_REFLECT, },
      .upper = { .type = GKYL_SPECIES_REFLECT, },
    },
    
    .num_diag_moments = 3,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_LTE },
  };

  // Vlasov-Maxwell app.
  struct gkyl_vm app_inp = {

   .cdim = 2, .vdim = 2, 
   .lower = { ctx.pi / 4.0, 0.0 },
   .upper = { (ctx.pi / 4.0) + ctx.Ltheta, ctx.Lphi },
   .cells = { NTHETA, NPHI },

   .poly_order = ctx.poly_order,
   .basis_type = app_args.basis_type,
   .cfl_frac = ctx.cfl_frac,

   .num_periodic_dir = 1,
   .periodic_dirs = { 1 },

   .num_species = 1,
   .species = { neut },

   .skip_field = true,

   .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0], app_args.cuts[1] },
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
  gkyl_vlasov_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

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
