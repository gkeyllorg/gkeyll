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

struct cylindrical_sodshock_ctx
{
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using normalized code units).
  double mass; // Neutral mass.
  double charge; // Neutral charge.

  double nl; // Left/inner number density.
  double Tl; // Left/inner temperature.
  double V_r_drift_l; // Left/inner drift velocity (radial direction).
  double V_theta_drift_l; // Left/inner drift velocity (angular direction).
  double V_z_drift_l; // Left/inner drift velocity (z-direction).

  double nr; // Right/outer number density.
  double Tr; // Right/outer temperature.
  double V_r_drift_r; // Right/outer drift velocity (radial direction).
  double V_theta_drift_r; // Right/outer drift velocity (angular direction).
  double V_z_drift_r; // Right/outer drift velocity (z-direction).

  double vt; // Thermal velocity.
  double nu; // Collision frequency.

  // Simulation parameters.
  int Nr; // Cell count (configuration space: radial direction).
  int Ntheta; // Cell count (configuration space: angular direction).
  int Nz; // Cell count (configuration space: z-direction).
  int Nvr; // Cell count (velocity space: radial direction).
  int Nvtheta; // Cell count (velocity space: angular direction).
  int Nvz; // Cell count (velocity space: z-direction).
  double Lr; // Domain size (configuration space: radial direction).
  double Ltheta; // Domain size (configuration space: angular direction).
  double Lz; // Domain size (configuration space: z-direction).
  double vr_max; // Domain boundary (velocity space: radial direction).
  double vtheta_max; // Domain boundary (velocity space: angular direction).
  double vz_max; // Domain boundary (velocity space: z-direction).
  int poly_order; // Polynomial order.
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  int field_energy_calcs; // Number of times to calculate field energy.
  int integrated_mom_calcs; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.

  double midplane; // Radial midplane location designating jump in quantities.
};

struct cylindrical_sodshock_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  double mass = 1.0; // Neutral mass.
  double charge = 0.0; // Neutral charge.

  double nl = 1.0; // Left/inner number density.
  double Tl = 1.0; // Left/inner temperature.
  double V_r_drift_l = 0.0; // Left/inner drift velocity (radial direction).
  double V_theta_drift_l = 0.0; // Left/inner drift velocity (angular direction).
  double V_z_drift_l = 0.0; // Left/inner drift velocity (z-direction).

  double nr = 0.125; // Right/outer number density.
  double Tr = sqrt(0.1 / 0.125); // Right/outer temperature.
  double V_r_drift_r = 0.0; // Right/outer drift velocity (radial direction).
  double V_theta_drift_r = 0.0; // Right/outer drift velocity (angular direction).
  double V_z_drift_r = 0.0; // Right/outer drift velocity (z-direction).

  double vt = 1.0; // Thermal velocity.
  double nu = 15000.0; // Collision frequency.

  // Simulation parameters.
  int Nr = 32; // Cell count (configuration space: radial direction).
  int Ntheta = 1; // Cell count (configuration space: angular direction).
  int Nz = 1; // Cell count (configuration space: z-direction).
  int Nvr = 8; // Cell count (velocity space: radial direction).
  int Nvtheta = 8; // Cell count (velocity space: angular direction).
  int Nvz = 8; // Cell count (velocity space: z-direction).
  double Lr = 1.0; // Domain size (configuration space: radial direction).
  double Ltheta = 2.0 * pi; // Domain size (configuration space: angular direction).
  double Lz = 2.0; // Domain size (configuration space: z-direction).
  double vr_max = 8.0 * vt; // Domain boundary (velocity space: radial direction).
  double vtheta_max = 8.0 * vt; // Domain boundary (velocity space: angular direction).
  double vz_max = 8.0 * vt; // Domain boundary (velocity space: z-direction).
  int poly_order = 1; // Polynomial order.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 0.1; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs = INT_MAX; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  double midplane = 1.0; // Radial midplane location designating jump in quantities.

  struct cylindrical_sodshock_ctx ctx = {
    .pi = pi,
    .mass = mass,
    .charge = charge,
    .nl = nl,
    .Tl = Tl,
    .V_r_drift_l = V_r_drift_l,
    .V_theta_drift_l = V_theta_drift_l,
    .V_z_drift_l = V_z_drift_l,
    .nr = nr,
    .Tr = Tr,
    .V_r_drift_r = V_r_drift_r,
    .V_theta_drift_r = V_theta_drift_r,
    .V_z_drift_r = V_z_drift_r,
    .vt = vt,
    .nu = nu,
    .Nr = Nr,
    .Ntheta = Ntheta,
    .Nz = Nz,
    .Nvr = Nvr,
    .Nvtheta = Nvtheta,
    .Nvz = Nvz,
    .Lr = Lr,
    .Ltheta = Ltheta,
    .Lz = Lz,
    .vr_max = vr_max,
    .vtheta_max = vtheta_max,
    .vz_max = vz_max,
    .poly_order = poly_order,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .num_frames = num_frames,
    .field_energy_calcs = field_energy_calcs,
    .integrated_mom_calcs = integrated_mom_calcs,
    .integrated_L2_f_calcs = integrated_L2_f_calcs,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
    .midplane = midplane,
  };

  return ctx;
}

void
evalDensityInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct cylindrical_sodshock_ctx *app = ctx;
  double r = xn[0];

  double nl = app->nl;
  double nr = app->nr;
  double midplane = app->midplane;

  double n = 0.0;

  if (r < midplane) {
    n = nl; // Total number density (left/inner).
  }
  else {
    n = nr; // Total number density (right/outer).
  }

  double metric_det = r;

  // Set total number density.
  fout[0] = metric_det * n;
}

void
evalTempInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct cylindrical_sodshock_ctx *app = ctx;
  double r = xn[0];

  double Tl = app->Tl;
  double Tr = app->Tr;
  double midplane = app->midplane;

  double T = 0.0;

  if (r < midplane) {
    T = Tl; // Isotropic temperature (left/inner).
  }
  else {
    T = Tr; // Isotropic temperature (right/outer).
  }

  // Set isotropic temperature.
  fout[0] = T;
}

void
evalVDriftInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct cylindrical_sodshock_ctx *app = ctx;
  double r = xn[0];

  double V_r_drift_l = app->V_r_drift_l;
  double V_theta_drift_l = app->V_theta_drift_l;
  double V_z_drift_l = app->V_z_drift_l;

  double V_r_drift_r = app->V_r_drift_r;
  double V_theta_drift_r = app->V_theta_drift_r;
  double V_z_drift_r = app->V_z_drift_r;

  double midplane = app->midplane;

  double V_r_drift = 0.0;
  double V_theta_drift = 0.0;
  double V_z_drift = 0.0;

  if (r < midplane) {
    V_r_drift = V_r_drift_l; // Radial drift velocity (left/inner).
    V_theta_drift = V_theta_drift_l; // Angular drift velocity (left/inner).
    V_z_drift = V_z_drift_l; // Z drift velocity (left/inner).
  }
  else {
    V_r_drift = V_r_drift_r; // Radial drift velocity (right/outer).
    V_theta_drift = V_theta_drift_r; // Angular drift velocity (right/outer).
    V_z_drift = V_z_drift_r; // Z drift velocity (right/outer).
  }

  // Set total drift velocity.
  fout[0] = V_r_drift; fout[1] = V_theta_drift; fout[2] = V_z_drift;
}

void
evalNu(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct cylindrical_sodshock_ctx *app = ctx;

  double nu = app->nu;

  // Set collision frequency.
  fout[0] = nu;
}

void
evalHamiltonian(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double q_r = xn[0], p_r_dot = xn[3], p_theta_dot = xn[4], p_z_dot = xn[5];

  double inv_metric_r_r = 1.0;
  double inv_metric_r_theta = 0.0;
  double inv_metric_r_z = 0.0;
  double inv_metric_theta_theta = 1.0 / (q_r * q_r);
  double inv_metric_theta_z = 0.0;
  double inv_metric_z_z = 1.0;

  double hamiltonian = (0.5 * inv_metric_r_r * p_r_dot * p_r_dot) + (0.5 * (2.0 * inv_metric_r_theta * p_r_dot * p_theta_dot)) +
    (0.5 * (2.0 * inv_metric_r_z * p_r_dot * p_z_dot)) + (0.5 * inv_metric_theta_theta * p_theta_dot * p_theta_dot) +
    (0.5 * (2.0 * inv_metric_theta_z * p_theta_dot * p_z_dot)) + (0.5 * inv_metric_z_z * p_z_dot * p_z_dot); // Canonical Hamiltonian.
  
  // Set canonical Hamiltonian.
  fout[0] = hamiltonian;
}

void
evalInvMetric(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double q_r = xn[0];

  double inv_metric_r_r = 1.0; // Inverse metric tensor (radial-radial component).
  double inv_metric_r_theta = 0.0; // Inverse metric tensor (radial-angular component).
  double inv_metric_r_z = 0.0; // Inverse metric tensor (radial-z component).
  double inv_metric_theta_theta = 1.0 / (q_r * q_r); // Inverse metric tensor (angular-angular component).
  double inv_metric_theta_z = 0.0; // Inverse metric tensor (angular-z component).
  double inv_metric_z_z = 1.0; // Inverse metric tensor (z-z component).
  
  // Set inverse metric tensor.
  fout[0] = inv_metric_r_r; fout[1] = inv_metric_r_theta; fout[2] = inv_metric_r_z;
  fout[3] = inv_metric_theta_theta; fout[4] = inv_metric_theta_z; fout[5] = inv_metric_z_z;
}

void
evalMetric(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double q_r = xn[0];

  double metric_r_r = 1.0; // Metric tensor (radial-radial component).
  double metric_r_theta = 0.0; // Metric tensor (radial-angular component).
  double metric_r_z = 0.0; // Metric tensor (radial-z component).
  double metric_theta_theta = q_r * q_r; // Metric tensor (angular-angular component).
  double metric_theta_z = 0.0; // Metric tensor (angular-z component).
  double metric_z_z = 1.0; // Metric tensor (z-z component).
  
  // Set metric tensor.
  fout[0] = metric_r_r; fout[1] = metric_r_theta; fout[2] = metric_r_z;
  fout[3] = metric_theta_theta; fout[4] = metric_theta_z; fout[5] = metric_z_z;
}

void
evalMetricDet(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double q_r = xn[0];

  double metric_det = q_r; // Metric tensor determinant.
  
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

  struct cylindrical_sodshock_ctx ctx = create_ctx(); // Context for initialization functions.

  int NR = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nr);
  int NTHETA = APP_ARGS_CHOOSE(app_args.xcells[1], ctx.Ntheta);
  int NZ = APP_ARGS_CHOOSE(app_args.xcells[2], ctx.Nz);
  int NVR = APP_ARGS_CHOOSE(app_args.vcells[0], ctx.Nvr);
  int NVTHETA = APP_ARGS_CHOOSE(app_args.vcells[1], ctx.Nvtheta);
  int NVZ = APP_ARGS_CHOOSE(app_args.vcells[2], ctx.Nvz);

  int nrank = 1; // Number of processors in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif  

  int ccells[] = { NR, NTHETA, NZ };
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
  struct gkyl_vlasov_kinetic_species neut = {
    .model_id = GKYL_MODEL_CANONICAL_PB,
    .lower = { -ctx.vr_max, -ctx.vtheta_max, -ctx.vz_max },
    .upper = { ctx.vr_max, ctx.vtheta_max, ctx.vz_max },
    .cells = { NVR, NVTHETA, NVZ },

    .hamil = evalHamiltonian,
    .hamil_ctx = &ctx,
    .h_ij = evalMetric,
    .h_ij_ctx = &ctx,
    .h_ij_inv = evalInvMetric,
    .h_ij_inv_ctx = &ctx,
    .det_h = evalMetricDet,
    .det_h_ctx = &ctx,

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
      .self_nu_ctx = &ctx,
      .is_implicit = true,
    },

    .correct = {
      .correct_all_moms = true,
      .iter_eps = 1.0e-12,
      .max_iter = 100,
      .use_last_converged = false,
    },

    .bcx = {
      .lower = { .type = GKYL_SPECIES_REFLECT, },
      .upper = { .type = GKYL_SPECIES_REFLECT, },
    },
    
    .num_diag_moments = 4,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_LTE, GKYL_F_MOMENT_ENERGY },
  };

  // Vlasov-Maxwell app.
  struct gkyl_vm app_inp = {

   .cdim = 3, .vdim = 3,
   .lower = { 0.5, 0.0, -0.5 * ctx.Lz },
   .upper = { 0.5 + ctx.Lr, ctx.Ltheta, 0.5 * ctx.Lz },
   .cells = { NR, NTHETA, NZ },

   .poly_order = ctx.poly_order,
   .basis_type = app_args.basis_type,
   .cfl_frac = ctx.cfl_frac,

   .num_periodic_dir = 2,
   .periodic_dirs = { 1, 2 },

   .num_species = 1,
    .species = {
      { .name = "neut", .charge = ctx.charge, .mass = ctx.mass,
        .type = GKYL_SPECIES_VLASOV, .kinetic = neut },
    },

   .skip_field = true,

   .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0], app_args.cuts[1], app_args.cuts[2] },
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
