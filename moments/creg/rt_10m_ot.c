#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_ten_moment.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

#include <rt_arg_parse.h>
#include <gkyl_kann_net.h>
#include <gkyl_knutils.h>

struct ot_ctx
{
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using normalized code units).
  double epsilon0; // Permittivity of free space.
  double mu0; // Permeability of free space.
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.
  double mass_ion; // Ion mass.
  double charge_ion; // Ion charge.

  double n0; // Reference number density.
  double vAe; // Electron Alfven velocity.
  double beta; // Plasma beta.

  // Derived physical quantities (using normalized code units).
  double B0; // Reference magnetic field strength.
  double vAi; // Ion Alfven velocity.

  double vte; // Electron thermal velocity.
  double vti; // Ion thermal velocity.

  double omega_ci; // Ion cyclotron frequency.
  double d_i; // Ion sound inertial length.

  double delta_B0; // Reference magnetic field strength perturbation.
  double delta_u0; // Reference fluid velocity perturbation.

  // Simulation parameters.
  int Nx; // Cell count (x-direction).
  int Ny; // Cell count (y-direction).
  double Lx; // Domain size (x-direction).
  double Ly; // Domain size (y-direction).
  double k0; // Closure parameter.
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  int field_energy_calcs; // Number of times to calculate field energy.
  int integrated_mom_calcs; // Number of times to calculate integrated moments.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.

  // Neural network parameters.
  bool use_nn_closure; // Use neural network-based closure?
  int poly_order; // Polynomial order of learned DG coefficients.
  const char* nn_closure_file; // File path of neural network to use.
};

struct ot_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  double epsilon0 = 1.0; // Permittivity of free space.
  double mu0 = 1.0; // Permeability of free space.
  double mass_elc = 1.0; // Electron mass.
  double charge_elc = -1.0; // Electron charge.
  double mass_ion = 25.0; // Ion mass.
  double charge_ion = 1.0; // Ion charge.

  double n0 = 1.0; // Reference number density.
  double vAe = 0.5; // Electron Alfven velocity.
  double beta = 0.08; // Plasma beta.

  // Derived physical quantities (using normalized code units).
  double B0 = vAe * sqrt(mu0 * n0 * mass_elc); // Reference magnetic field strength.
  double vAi = vAe / sqrt(mass_ion); // Ion Alfven velocity.

  double vte = vAe * sqrt(beta / 2.0); // Electron thermal velocity.
  double vti = vte / sqrt(mass_ion); // Ion thermal velocity.

  double omega_ci = charge_ion * B0 / mass_ion; // Ion cyclotron frequency.
  double d_i = vAi / omega_ci; // Ion sound inertial length.

  double delta_B0 = 0.2 * B0; // Reference magnetic field strength perturbation.
  double delta_u0 = 0.2 * vAi; // Reference fluid velocity perturbation.

  // Simulation parameters.
  int Nx = 128; // Cell count (x-direction).
  int Ny = 128; // Cell count (y-direction).
  double Lx = 20.48 * d_i; // Domain size (x-direction).
  double Ly = 20.48 * d_i; // Domain size (y-direction).
  double k0 = 5.0; // Closure parameter.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 75.0 / omega_ci; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  // Neural network parameters.
  bool use_nn_closure = false; // Use neural network-based closure?
  int poly_order = 1; // Polynomial order of learned DG coefficients.
  const char* nn_closure_file = "moments/data/neural_nets/pkpm_ot_p1_moms_nn_1"; // File path of neural network to use.

  struct ot_ctx ctx = {
    .pi = pi,
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .mass_ion = mass_ion,
    .charge_ion = charge_ion,
    .n0 = n0,
    .vAe = vAe,
    .beta = beta,
    .B0 = B0,
    .vAi = vAi,
    .vte = vte,
    .vti = vti,
    .omega_ci = omega_ci,
    .d_i = d_i,
    .delta_B0 = delta_B0,
    .delta_u0 = delta_u0,
    .Nx = Nx,
    .Ny = Ny,
    .Lx = Lx,
    .Ly = Ly,
    .k0 = k0,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .num_frames = num_frames,
    .field_energy_calcs = field_energy_calcs,
    .integrated_mom_calcs = integrated_mom_calcs,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
    .use_nn_closure = use_nn_closure,
    .poly_order = poly_order,
    .nn_closure_file = nn_closure_file,
  };

  return ctx;
}

void
evalElcInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct ot_ctx *app = ctx;
  double x = xn[0], y = xn[1];

  double pi = app->pi;
  double mu0 = app->mu0;
  double mass_elc = app->mass_elc;
  double charge_ion = app->charge_ion;

  double n0 = app->n0;
  double vte = app->vte;

  double delta_B0 = app->delta_B0;
  double delta_u0 = app->delta_u0;

  double Lx = app->Lx;
  double Ly = app->Ly;

  double Jz = (delta_B0 * (4.0 * pi / Lx) * cos(4.0 * pi * x / Lx) + delta_B0 * (2.0 * pi / Ly) * cos(2.0 * pi * y / Ly)) / mu0;

  double vx_drift = -delta_u0 * sin(2.0 * pi * y / Ly); // Electron drift velocity (x-direction).
  double vy_drift = delta_u0 * sin(2.0 * pi * x / Lx); // Electron drift velocity (y-direction).
  double vz_drift = -Jz / charge_ion; // Electron drift velocity (z-direction).

  double rhoe = n0 * mass_elc; // Electron mass density.
  double mome_x = mass_elc * vx_drift; // Electron total momentum density (x-direction).
  double mome_y = mass_elc * vy_drift; // Electron total momentum density (y-direction).
  double mome_z = mass_elc * vz_drift; // Electron total momentum density (z-direction).
  double pre = vte * vte * rhoe; // Electron pressure (scalar).

  double pre_xx = pre + (mome_x * mome_x) / rhoe; // Electron pressure tensor (xx-component).
  double pre_xy = 0.0; // Electron pressure tensor (xy-component).
  double pre_xz = 0.0; // Electron pressure tensor (xz-component).
  double pre_yy = pre + (mome_y * mome_y) / rhoe; // Electron pressure tensor (yy-component).
  double pre_yz = 0.0; // Electron pressure tensor (yz-component).
  double pre_zz = pre + (mome_z * mome_z) / rhoe; // Electron pressure tensor (zz-component).

  // Set electron mass density.
  fout[0] = rhoe;
  // Set electron momentum density.
  fout[1] = mome_x; fout[2] = mome_y; fout[3] = mome_z;
  // Set electron pressure tensor.
  fout[4] = pre_xx; fout[5] = pre_xy; fout[6] = pre_xz;
  fout[7] = pre_yy; fout[8] = pre_yz; fout[9] = pre_zz;
}

void
evalIonInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct ot_ctx *app = ctx;
  double x = xn[0], y = xn[1];

  double pi = app->pi;
  double mass_ion = app->mass_ion;

  double n0 = app->n0;
  double vti = app->vti;

  double delta_u0 = app->delta_u0;

  double Lx = app->Lx;
  double Ly = app->Ly;

  double vx_drift = -delta_u0 * sin(2.0 * pi * y / Ly); // Ion drift velocity (x-direction).
  double vy_drift = delta_u0 * sin(2.0 * pi * x / Lx); // Ion drift velocity (y-direction).
  double vz_drift = 0.0; // Ion drift velocity (z-direction).

  double rhoi = n0 * mass_ion; // Ion mass density.
  double momi_x = mass_ion * vx_drift; // Ion total momentum density (x-direction).
  double momi_y = mass_ion * vy_drift; // Ion total momentum density (y-direction).
  double momi_z = mass_ion * vz_drift; // Ion total momentum density (z-direction).
  double pri = vti * vti * rhoi; // Ion pressure (scalar).

  double pri_xx = pri + (momi_x * momi_x) / rhoi; // Ion pressure tensor (xx-component).
  double pri_xy = 0.0; // Ion pressure tensor (xy-component).
  double pri_xz = 0.0; // Ion pressure tensor (xz-component).
  double pri_yy = pri + (momi_y * momi_y) / rhoi; // Ion pressure tensor (yy-component).
  double pri_yz = 0.0; // Ion pressure tensor (yz-component).
  double pri_zz = pri + (momi_z * momi_z) / rhoi; // Ion pressure tensor (zz-component).

  // Set ion mass density.
  fout[0] = rhoi;
  // Set ion momentum density.
  fout[1] = momi_x; fout[2] = momi_y; fout[3] = momi_z;
  // Set ion pressure tensor.
  fout[4] = pri_xx; fout[5] = pri_xy; fout[6] = pri_xz;
  fout[7] = pri_yy; fout[8] = pri_yz; fout[9] = pri_zz;
}

void
evalFieldInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct ot_ctx *app = ctx;
  double x = xn[0], y = xn[1];

  double pi = app->pi;
  double mu0 = app->mu0;
  double mass_elc = app->mass_elc;
  double charge_ion = app->charge_ion;

  double B0 = app->B0;
  double delta_B0 = app->delta_B0;
  double delta_u0 = app->delta_u0;

  double Lx = app->Lx;
  double Ly = app->Ly;

  double Jz = (delta_B0 * (4.0 * pi / Lx) * cos(4.0 * pi * x / Lx) + delta_B0 * (2.0 * pi / Ly) * cos(2.0 * pi * y / Ly)) / mu0;

  double Bx = -delta_B0 * sin(2.0 * pi * y / Ly); // Total magnetic field (x-direction).
  double By = delta_B0 * sin(4.0 * pi * x / Lx); // Total magnetic field (y-direction).
  double Bz = B0; // Total magnetic field (z-direction).

  double vx_drift = -delta_u0 * sin(2.0 * pi * y / Ly); // Electron drift velocity (x-direction).
  double vy_drift = delta_u0 * sin(2.0 * pi * x / Lx); // Electron drift velocity (y-direction).
  double vz_drift = -Jz / charge_ion; // Electron drift velocity (z-direction).

  double Ex = -((vy_drift * Bz) - (vz_drift * By)); // Total electric field (x-direction).
  double Ey = -((vz_drift * Bx) - (vx_drift * Bz)); // Total electric field (y-direction).
  double Ez = -((vx_drift * By) - (vy_drift * Bx)); // Total electric field (z-direction).
  
  // Set electric field.
  fout[0] = Ex; fout[1] = Ey, fout[2] = Ez;
  // Set magnetic field.
  fout[3] = Bx; fout[4] = By; fout[5] = Bz;
  // Set correction potentials.
  fout[6] = 0.0; fout[7] = 0.0;
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_moment_app* app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_write) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_moment_app_write(app, t_curr, frame);
    gkyl_moment_app_write_field_energy(app);
    gkyl_moment_app_write_integrated_mom(app);
  }
}

void
calc_field_energy(struct gkyl_tm_trigger* fet, gkyl_moment_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(fet, t_curr) || force_calc) {
    gkyl_moment_app_calc_field_energy(app, t_curr);
  }
}

void
calc_integrated_mom(struct gkyl_tm_trigger* imt, gkyl_moment_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(imt, t_curr) || force_calc) {
    gkyl_moment_app_calc_integrated_mom(app, t_curr);
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

  struct ot_ctx ctx = create_ctx(); // Context for initialization functions.
    
  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], ctx.Ny);

  struct gkyl_kann_net **ann = gkyl_malloc(sizeof(struct gkyl_kann_net*) * 2);
  if (ctx.use_nn_closure) {
    const char *fmt_elc = "%s-%s.dat";
    int sz_elc = gkyl_calc_strlen(fmt_elc, ctx.nn_closure_file, "elc");
    char fileNm_elc[sz_elc + 1];
    snprintf(fileNm_elc, sizeof fileNm_elc, fmt_elc, ctx.nn_closure_file, "elc");
    FILE *file_elc = fopen(fileNm_elc, "r");
    if (file_elc != NULL) {
      ann[0] = gkyl_kann_net_load(fileNm_elc, app_args.use_gpu);
      fclose(file_elc);
    }
    else {
      ann[0] = 0;
      ctx.use_nn_closure = false;
      fprintf(stderr, "Neural network for elc species not found! Disabling NN-based closure.\n");
    }

    const char *fmt_ion = "%s-%s.dat";
    int sz_ion = gkyl_calc_strlen(fmt_ion, ctx.nn_closure_file, "ion");
    char fileNm_ion[sz_ion + 1];
    snprintf(fileNm_ion, sizeof fileNm_ion, fmt_ion, ctx.nn_closure_file, "ion");
    FILE *file_ion = fopen(fileNm_ion, "r");
    if (file_ion != NULL) {
      ann[1] = gkyl_kann_net_load(fileNm_ion, app_args.use_gpu);
      fclose(file_ion);
    }
    else {
      ann[1] = 0;
      ctx.use_nn_closure = false;
      fprintf(stderr, "Neural network for ion species not found! Disabling NN-based closure.\n");
    }
  }
    
  // Electron/ion equations.
  struct gkyl_wv_eqn *elc_ten_moment = gkyl_wv_ten_moment_new(ctx.k0, false, ctx.use_nn_closure, ctx.poly_order, ann[0], app_args.use_gpu);
  struct gkyl_wv_eqn *ion_ten_moment = gkyl_wv_ten_moment_new(ctx.k0, false, ctx.use_nn_closure, ctx.poly_order, ann[1], app_args.use_gpu);

  struct gkyl_moment_species elc = {
    .name = "elc",
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .equation = elc_ten_moment,

    .init = evalElcInit,
    .ctx = &ctx,
  };
    
  struct gkyl_moment_species ion = {
    .name = "ion",
    .charge = ctx.charge_ion, .mass = ctx.mass_ion,
    .equation = ion_ten_moment,

    .init = evalIonInit,
    .ctx = &ctx,
  };

  // Field.
  struct gkyl_moment_field field = {
    .epsilon0 = ctx.epsilon0, .mu0 = ctx.mu0,
    .mag_error_speed_fact = 1.0,
        
    .init = evalFieldInit,
    .ctx = &ctx,
  };
  int nrank = 1; // Number of processes in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif

  // Create global range.
  int cells[] = { NX, NY };
  int dim = sizeof(cells) / sizeof(cells[0]);

  int cuts[dim];
#ifdef GKYL_HAVE_MPI
  for (int d = 0; d < dim; d++) {
    if (app_args.use_mpi) {
      cuts[d] = app_args.cuts[d];
    }
    else {
      cuts[d] = 1;
    }
  }
#else
  for (int d = 0; d < dim; d++) {
    cuts[d] = 1;
  }
#endif

  // Construct communicator for use in app.
  struct gkyl_comm *comm;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
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
  for (int d = 0; d < dim; d++) {
    ncuts *= cuts[d];
  }

  if (ncuts != comm_size) {
    if (my_rank == 0) {
      fprintf(stderr, "*** Number of ranks, %d, does not match total cuts, %d!\n", comm_size, ncuts);
    }
    goto mpifinalize;
  }

  // Moment app.
  struct gkyl_moment app_inp = {

    .ndim = 2,
    .lower = { 0.0, 0.0 },
    .upper = { ctx.Lx, ctx.Ly },
    .cells = { NX, NY },

    .num_periodic_dir = 2,
    .periodic_dirs = { 0, 1 },
    .cfl_frac = ctx.cfl_frac,

    .num_species = 2,
    .species = { elc, ion },

    .field = field,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0], app_args.cuts[1] },
      .comm = comm,
    },
  };

  // Create app object.
  // Set app output name from the executable name (argv[0]).
  snprintf(app_inp.name, sizeof(app_inp.name), "%s", app_args.app_name);
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  // Initial and final simulation times.
  double t_curr = 0.0, t_end = ctx.t_end;

  // Initialize simulation.
  int frame_curr = 0;
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_moment_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_moment_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n", gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_moment_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_moment_app_cout(app, stdout, " at time = %g\n", t_curr);
  }
  else {
    gkyl_moment_app_apply_ic(app, t_curr);
  }

  // Create trigger for field energy.
  int field_energy_calcs = ctx.field_energy_calcs;
  struct gkyl_tm_trigger fe_trig = { .dt = t_end / field_energy_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_field_energy(&fe_trig, app, t_curr, false);

  // Create trigger for integrated moments.
  int integrated_mom_calcs = ctx.integrated_mom_calcs;
  struct gkyl_tm_trigger im_trig = { .dt = t_end / integrated_mom_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_integrated_mom(&im_trig, app, t_curr, false);

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
    gkyl_moment_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = gkyl_moment_update(app, dt);
    gkyl_moment_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      gkyl_moment_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_field_energy(&fe_trig, app, t_curr, false);
    calc_integrated_mom(&im_trig, app, t_curr, false);
    write_data(&io_trig, app, t_curr, false);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_moment_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_moment_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_moment_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_moment_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_moment_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);

        calc_field_energy(&fe_trig, app, t_curr, true);
        calc_integrated_mom(&im_trig, app, t_curr, true);
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
  write_data(&io_trig, app, t_curr, false);
  if (ctx.use_nn_closure) {
    for (int i = 0; i < app_inp.num_species; i++) {
      gkyl_kann_net_release(ann[i]);
    }
  }
  gkyl_moment_app_stat_write(app);

  struct gkyl_moment_stat stat = gkyl_moment_app_stat(app);

  gkyl_moment_app_cout(app, stdout, "\n");
  gkyl_moment_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_moment_app_cout(app, stdout, "Number of failed time-steps %ld\n", stat.nfail);
  gkyl_moment_app_cout(app, stdout, "Species updates took %g secs\n", stat.species_tm);
  gkyl_moment_app_cout(app, stdout, "Field updates took %g secs\n", stat.field_tm);
  gkyl_moment_app_cout(app, stdout, "Source updates took %g secs\n", stat.sources_tm);
  gkyl_moment_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

freeresources:
  // Free resources after simulation completion.
  gkyl_wv_eqn_release(elc_ten_moment);
  gkyl_wv_eqn_release(ion_ten_moment);
  gkyl_comm_release(comm);
  gkyl_moment_app_release(app);
  gkyl_free(ann);
  
mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif
  
  return 0;
}