// 2D ring-accretion problem onto a static (Schwarzschild) black hole, for the general relativistic can-pb.
// Input parameters describe an asymmetrical ring of cold relativistic gas accreting onto a non-rotating black hole.

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

struct star_newtonian_static_ctx
{
  // Mathematical constants (dimensionless).
  double pi;

  // Spacetime parameters (using geometric units).
  double star_mass; // Mass of the black hole.

  // Physical constants (using normalized code units) For particles.
  double mass; // Neutral mass.
  double charge; // Neutral charge.
  double nu; // Collision frequency.
  double vt; // Thermal velocity.

  // Simulation parameters.
  int Nr; // Cell count (configuration space: radial direction).
  int Ntheta; // Cell count (configuration space: azimuthal angular direction).
  int Nvr; // Cell count (velocity space: radial direction).
  int Nvtheta; // Cell count (velocity space: azimuthal angular direction).
  double Lr_min; // Domain size maximum radius (configuration space: radial direction).
  double Lr_max; // Domain size maximum radius (configuration space: radial direction).
  double Ltheta_min; // Domain size min (configuration space: azimuthal angular direction).
  double Ltheta_max; // Domain size max (configuration space: azimuthal angular direction).
  double v_r_max; // Domain boundary (velocity space: radial direction).
  double v_theta_max; // Domain boundary (velocity space: azimuthal angular direction).
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

struct star_newtonian_static_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  double mass = 1.0; // Neutral mass.
  double charge = 0.0; // Neutral charge.
  double nu = 2*pi; // Collision frequency.
    
  // Spacetime parameters (using geometric units).
  double star_mass = 1.0; // Mass of the black hole.
  double bh_spin = 0.0; // Spin of the black hole.

  double bh_pos_x = 0.0; // Position of the black hole (x-direction).
  double bh_pos_y = 0.0; // Position of the black hole (y-direction).
  double bh_pos_z = 0.0; // Position of the black hole (z-direction).

  double vt = 1.0; // Thermal velocity.

  // Simulation parameters.
  int Nr = 16; // Cell count (configuration space: radial direction).
  int Ntheta = 32; // Cell count (configuration space: azimuthal angular direction).
  int Nvr = 8; // Cell count (velocity space: radial direction).
  int Nvtheta = 32; // Cell count (velocity space: azimuthal angular direction).
  double Lr_min = 0.5; // Domain size radius min (configuration space: radial direction).
  double Lr_max = 2.0; // Domain size radius max (configuration space: radial direction).
  double Ltheta_min = 0.0; // Domain size minimum (configuration space: azimuthal angular direction).
  double Ltheta_max = 2.0 * pi; // Domain size maximum (configuration space: azimuthal angular direction).
  double v_r_max = 0.5 * vt; // Domain boundary (velocity space: radial direction).
  double v_theta_max = 3.0 * vt; // Domain boundary (velocity space: azimuthal angular direction).
  int poly_order = 2; // Polynomial order.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 0.001; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs = INT_MAX; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct star_newtonian_static_ctx ctx = {
    .pi = pi,
    .mass = mass,
    .charge = charge,
    .star_mass = star_mass,
    .vt = vt,
    .nu = nu,
    .Nr = Nr,
    .Ntheta = Ntheta,
    .Nvr = Nvr,
    .Nvtheta = Nvtheta,
    .Lr_min = Lr_min,
    .Lr_max = Lr_max,
    .Ltheta_min = Ltheta_min,
    .Ltheta_max = Ltheta_max,
    .v_r_max = v_r_max,
    .v_theta_max = v_theta_max,
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
evalNu(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct star_newtonian_static_ctx *app = ctx;

  double nu = app->nu;

  // Set collision frequency.
  fout[0] = nu;
}

void
evalInitialf(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct star_newtonian_static_ctx *app = ctx;
  double q_r = xn[0], q_theta = xn[1], p_r_dot = xn[2], p_theta_dot = xn[3];
  double metric_det = q_r; // Metric tensor determinant.

  double dr = (app->Lr_max - app->Lr_min)/app->Nr;
  double dtheta = (app->Ltheta_max - app->Ltheta_min)/app->Ntheta;
  double dpr = (2.0*app->v_r_max)/app->Nvr;
  double dptheta = (2.0*app->v_theta_max)/app->Nvtheta;

  double pi = app->pi;
  double f = 0.0; // background distribution function
 
  // Initalize f's positisional amplitude in space  
  // Init with a sqaure function in momentum space
  // around orbital speed p_theta = r (for circular orbits)
  // and r_dot = 0
  const double sigma = 0.15;
  if ((p_r_dot > - 0.1) && (p_r_dot < 0.1)){
    if ((p_theta_dot > sqrt(q_r - 0.1)) && (p_theta_dot < sqrt(q_r + 0.1))){
      if (q_r < 1.0) {
          f = 0.1 + 0.9 * exp(-((q_r - 1.0) * (q_r - 1.0)) / (2.0 * sigma * sigma));
      } 
      else if (q_r > 1.5) {
          f =  0.1 + 0.9 * exp(-((q_r - 1.5) * (q_r - 1.5)) / (2.0 * sigma * sigma));
      } 
      else {
          f = 1.0;
      }
    }
  }

  //printf("dr: %1.3e, dtheta = %1.3e, dpr = %1.3e, dptheta = %1.3e \n",dr, dtheta, dpr, dptheta);

  // Set distribution function, multiply by jacobian for the correct density
  double jac = q_r;
  fout[0] = jac*f;
}

void 
evalHamiltonian(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct star_newtonian_static_ctx *app = ctx;

  double q_r = xn[0], q_theta = xn[1], p_r_dot = xn[2], p_theta_dot = xn[3];

  double inv_metric_r_r = 1.0; // Inverse metric tensor (radial-radial component).
  double inv_metric_r_theta = 0.0; // Inverse metric tensor (radial-angular component).
  double inv_metric_theta_theta = 1.0 / (q_r * q_r); // Inverse metric tensor (angular-angular component).
  
  // H = \alpha \gamma - \beta \cdot p
  double hamiltonian = (0.5 * inv_metric_r_r * p_r_dot * p_r_dot) + (0.5 * (2.0 * inv_metric_r_theta * p_r_dot * p_theta_dot)) +
    (0.5 * inv_metric_theta_theta * p_theta_dot * p_theta_dot) - 1.0/q_r; // Canonical Hamiltonian, with gravity (Newtonian)
  
  // Set canonical Hamiltonian.
  fout[0] = hamiltonian;
}

void
evalInvMetric(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct star_newtonian_static_ctx *app = ctx;
  double q_r = xn[0];

  double inv_metric_r_r = 1.0; // Inverse metric tensor (radial-radial component).
  double inv_metric_r_theta = 0.0; // Inverse metric tensor (radial-angular component).
  double inv_metric_theta_theta = 1.0 / (q_r * q_r); // Inverse metric tensor (angular-angular component).
  
  // Set inverse metric tensor.
  fout[0] = inv_metric_r_r; fout[1] = inv_metric_r_theta; fout[2] = inv_metric_theta_theta;
}

void
evalMetric(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct star_newtonian_static_ctx *app = ctx;
  double q_r = xn[0];

  double metric_r_r = 1.0; //Metric tensor (radial-radial component).
  double metric_r_theta = 0.0; //Metric tensor (radial-angular component).
  double metric_theta_theta = q_r * q_r; //Metric tensor (angular-angular component).
  
  // Set metric tensor.
  fout[0] = metric_r_r; fout[1] = metric_r_theta; fout[2] = metric_theta_theta;
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

  struct star_newtonian_static_ctx ctx = create_ctx(); // Context for initialization functions.

  int Nr = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nr);
  int Ntheta = APP_ARGS_CHOOSE(app_args.xcells[1], ctx.Ntheta);
  int Nvr = APP_ARGS_CHOOSE(app_args.vcells[0], ctx.Nvr);
  int Nvtheta = APP_ARGS_CHOOSE(app_args.vcells[1], ctx.Nvtheta);

  int nrank = 1; // Number of processors in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif  

  int ccells[] = { Nr, Ntheta };
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
    .lower = { -ctx.v_r_max, -1.0 },
    .upper = { ctx.v_r_max, ctx.v_theta_max},
    .cells = { Nvr, Nvtheta },

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
      .proj_id = GKYL_PROJ_FUNC,
      .func = evalInitialf,
      .ctx_func = &ctx,
    },
    .collisions =  {
      .collision_id = GKYL_BGK_COLLISIONS,
      .self_nu = evalNu,
      .self_nu_ctx = &ctx,
      .is_implicit = true,
    },

    .correct = {
      .correct_all_moms = false,
      .iter_eps = 0.0,
      .max_iter = 0,
      .use_last_converged = false,
      .output_f_lte = false,
    },

    .bcx = {
      .lower = { .type = GKYL_SPECIES_ABSORB, },
      .upper = { .type = GKYL_SPECIES_ABSORB, },
    },
    
    .num_diag_moments = 2,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1_FROM_H },
  };

  // Vlasov-Maxwell app.
  struct gkyl_vm app_inp = {

   .cdim = 2, .vdim = 2, 
   .lower = { ctx.Lr_min, ctx.Ltheta_min },
   .upper = { ctx.Lr_max, ctx.Ltheta_max },
   .cells = { Nr, Ntheta },

   .poly_order = ctx.poly_order,
   .basis_type = app_args.basis_type,
   .cfl_frac = ctx.cfl_frac,

   .num_periodic_dir = 1,
   .periodic_dirs = { 1 },

   .num_species = 1,
    .species = {
      { .name = "neut", .charge = ctx.charge, .mass = ctx.mass,
        .type = GKYL_SPECIES_VLASOV, .kinetic = neut },
    },

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
