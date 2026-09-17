// Lower-Hybrid Drift Instability (LHDI) test for the 10-moment equations.
// Input parameters match the initial conditions in Section 5.2, from the article:
// J. Ng, A. Hakim, J. Juno and A. Bhattacharjee (2019), "Drift Instabilities in Thin Current Sheets Using a Two-Fluid Model With Pressure Tensor Effects",
// Journal of Geophysical Research: Space Physics, Volume 124 (5): 3331-3346.
// https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018JA026313

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

struct lhdi_ctx
{
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using normalized code units).
  double epsilon0; // Permittivity of free space.
  double mu0; // Permeability of free space.
  double mass_ion; // Ion mass.
  double charge_ion; // Ion charge.
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.
  double Te_over_Ti; // Electron temperature / ion temperature.

  double n0; // Reference number density.
  double nb_over_n0; // Background number density / reference number density.
  double vte; // Electron thermal velocity.

  double beta; // Electron plasma beta.

  double noise_amp; // Noise level for perturbation.
  double mode; // Wave mode to perturb with noise.

  // Derived physical quantities (using normalized code units).
  double Te; // Electron temperature.
  double Ti; // Ion temperature.
  double vti; // Ion thermal velocity.

  double vAe; // Electron Alfven velocity.
  double B0; // Reference magnetic field strength (derived from normalization of mass_elc and n0).
  double vAi; // Ion Alfven velocity.

  double omega_ci; // Ion cyclotron frequency.
  double omega_ce; // Electron cyclotron frequency.

  double larmor_ion; // Ion Larmor radius.
  double larmor_elc; // Electron Larmor radius.

  double l; // Current sheet width.
  
  double nb; // Background number density.
  double Te_frac; // Fraction of total temperature from electrons.
  double Ti_frac; // Fraction of total temperature from ions.
  double ix; // Current (x-direction).
  double iy; // Current (y-direction).

  // Simulation parameters.
  int Nx; // Cell count (x-direction).
  int Ny; // Cell count (y-direction).
  double Lx; // Domain size (x-direction).
  double Ly; // Domain size (y-direction).
  double k0_elc; // Closure parameter for electrons.
  double k0_ion; // Closure parameter for ions.
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

struct lhdi_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  double epsilon0 = 1.0; // Permittivity of free space.
  double mu0 = 1.0; // Permeability of free space.
  double mass_ion = 36.0; // Ion mass.
  double charge_ion = 1.0; // Ion charge.
  double mass_elc = 1.0; // Electron mass.
  double charge_elc = -1.0; // Electron charge.
  double Te_over_Ti = 0.1; // Electron temperature / ion temperature.

  double n0 = 1.0; // Reference number density.
  double nb_over_n0 = 0.001; // Background number density / reference number density.
  double vte = 0.06; // Electron thermal velocity.

  double beta = 1.0 / 11.0; // Electron plasma beta.

  double noise_amp = 0.0001; // Noise level for perturbation.
  double mode = 8.0; // Wave mode to perturb with noise.

  // Derived physical quantities (using normalized code units).
  double Te = vte * vte * mass_elc / 2.0; // Electron temperature.
  double Ti = Te / Te_over_Ti; // Ion temperature.
  double vti = sqrt(2.0 * Ti / mass_ion); // Ion thermal velocity.

  double vAe = vte / sqrt(beta); // Electron Alfven velocity.
  double B0 = vAe; // Reference magnetic field strength (derived from normalization of mass_elc and n0).
  double vAi = vAe / sqrt(mass_ion); // Ion Alfven velocity.

  double omega_ci = charge_ion * B0 / mass_ion; // Ion cyclotron frequency.
  double omega_ce = charge_ion * B0 / mass_elc; // Electron cyclotron frequency.

  double larmor_ion = vti / omega_ci; // Ion Larmor radius.
  double larmor_elc = vte / omega_ce; // Electron Larmor radius.

  double l = larmor_ion; // Current sheet width.
  
  double nb = n0 * nb_over_n0; // Background number density.
  double Te_frac = Te / (Te + Ti); // Fraction of total temperature from electrons.
  double Ti_frac = 1.0 - Te_frac; // Fraction of total temperature from ions.
  double ix = mode; // Current (x-direction).
  double iy = 1.0; // Current (y-direction).

  // Simulation parameters.
  int Nx = 64; // Cell count (x-direction).
  int Ny = 128; // Cell count (y-direction).
  double Lx = 6.4 * l; // Domain size (x-direction).
  double Ly = 12.8 * l; // Domain size (y-direction).
  double k0_elc = 1.0; // Closure parameter for electrons.
  double k0_ion = 1.0 / 6.0; // Closure parameter for ions.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 1100.0; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  // Neural network parameters.
  bool use_nn_closure = false; // Use neural network-based closure?
  int poly_order = 1; // Polynomial order of learned DG coefficients.
  const char* nn_closure_file = "moments/data/neural_nets/pkpm_ot_p1_moms_nn_1"; // File path of neural network to use.

  struct lhdi_ctx ctx = {
    .pi = pi,
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .mass_ion = mass_ion,
    .charge_ion = charge_ion,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .Te_over_Ti = Te_over_Ti,
    .n0 = n0,
    .nb_over_n0 = nb_over_n0,
    .vte = vte,
    .beta = beta,
    .noise_amp = noise_amp,
    .mode = mode,
    .Te = Te,
    .Ti = Ti,
    .vti = vti,
    .vAe = vAe,
    .B0 = B0,
    .vAi = vAi,
    .omega_ci = omega_ci,
    .omega_ce = omega_ce,
    .larmor_ion = larmor_ion,
    .larmor_elc = larmor_elc,
    .l = l,
    .nb = n0 * nb_over_n0,
    .Te_frac = Te_frac,
    .Ti_frac = Ti_frac,
    .ix = ix,
    .iy = iy,
    .Nx = Nx,
    .Ny = Ny,
    .Lx = Lx,
    .Ly = Ly,
    .k0_elc = k0_elc,
    .k0_ion = k0_ion,
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
  double x = xn[0], y = xn[1];
  struct lhdi_ctx *app = ctx;

  double pi = app->pi;

  double mass_elc = app->mass_elc;
  double charge_elc = app->charge_elc;

  double n0 = app->n0;

  double noise_amp = app->noise_amp;
  double mode = app->mode;
  
  double Te = app->Te;

  double B0 = app->B0;
  
  double l = app->l;

  double Te_frac = app->Te_frac;
  double ix = app->ix;
  double iy = app->iy;

  double Lx = app->Lx;
  double Ly = app->Ly;

  double sech_sq = (1.0 / cosh(y / l)) * (1.0 / cosh(y / l)); // Hyperbolic secant squared.

  double n = n0 * sech_sq; // Total number density.
  double Jx_noise = -noise_amp * (iy * pi / Ly) * sin(iy * pi * y / Ly) * sin(ix * 2.0 * pi* x / Lx) / mode; // Current density noise (x-direction).
  double Jy_noise = -noise_amp * (ix * 2.0 * pi / Lx) * cos(iy * pi * y / Ly) * cos(ix * 2.0 * pi * x / Lx) / mode; // Current density noise (y-direction).

  double Jx = (B0 / l) * (-sech_sq) + Jx_noise; // Total current density, with noise (x-direction).
  double Jy = Jy_noise; // Total current density, with noise (y-direction).

  double rhoe = n * mass_elc; // Electron mass density.
  double mome_x = (mass_elc / charge_elc) * Jx * Te_frac; // Electron momentum density (x-direction).
  double mome_y = (mass_elc / charge_elc) * Jy * Te_frac; // Electron momentum density (y-direction).
  double mome_z = 0.0; // Electron momentum density (z-direction).
  double pre = n * Te; // Electron pressure (scalar).

  double pre_xx = pre + (mome_x * mome_x) / rhoe; // Electron pressure tensor (xx-component).
  double pre_xy = (mome_x * mome_y) / rhoe; // Electron pressure tensor (xy-component).
  double pre_xz = 0.0; // Electron pressure tensor (xz-component).
  double pre_yy = pre + (mome_y * mome_y) / rhoe; // Electron pressure tensor (yy-component).
  double pre_yz = 0.0; // Electron pressure tensor (yz-component).
  double pre_zz = pre; // Electron pressure tensor (zz-component).

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
  double x = xn[0], y = xn[1];
  struct lhdi_ctx *app = ctx;

  double pi = app->pi;

  double mass_ion = app->mass_ion;
  double charge_ion = app->charge_ion;

  double n0 = app->n0;

  double noise_amp = app->noise_amp;
  double mode = app->mode;
  
  double Ti = app->Ti;

  double B0 = app->B0;
  
  double l = app->l;

  double Ti_frac = app->Ti_frac;
  double ix = app->ix;
  double iy = app->iy;

  double Lx = app->Lx;
  double Ly = app->Ly;

  double sech_sq = (1.0 / cosh(y / l)) * (1.0 / cosh(y / l)); // Hyperbolic secant squared.

  double n = n0 * sech_sq; // Total number density.
  double Jx_noise = -noise_amp * (iy * pi / Ly) * sin(iy * pi * y / Ly) * sin(ix * 2.0 * pi* x / Lx) / mode; // Current density noise (x-direction).
  double Jy_noise = -noise_amp * (ix * 2.0 * pi / Lx) * cos(iy * pi * y / Ly) * cos(ix * 2.0 * pi * x / Lx) / mode; // Current density noise (y-direction).

  double Jx = (B0 / l) * (-sech_sq) + Jx_noise; // Total current density, with noise (x-direction).
  double Jy = Jy_noise; // Total current density, with noise (y-direction).

  double rhoi = n * mass_ion; // Ion mass density.
  double momi_x = (mass_ion / charge_ion) * Jx * Ti_frac; // Ion momentum density (x-direction).
  double momi_y = (mass_ion / charge_ion) * Jy * Ti_frac; // Ion momentum density (y-direction).
  double momi_z = 0.0; // Ion momentum density (z-direction).
  double pri = n * Ti; // Ion pressure (scalar).

  double pri_xx = pri + (momi_x * momi_x) / rhoi; // Ion pressure tensor (xx-component).
  double pri_xy = (momi_x * momi_y) / rhoi; // Ion pressure tensor (xy-component).
  double pri_xz = 0.0; // Ion pressure tensor (xz-component).
  double pri_yy = pri + (momi_y * momi_y) / rhoi; // Ion pressure tensor (yy-component).
  double pri_yz = 0.0; // Ion pressure tensor (yz-component).
  double pri_zz = pri; // Ion pressure tensor (zz-component).

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
  double x = xn[0], y = xn[1];
  struct lhdi_ctx *app = ctx;

  double pi = app->pi;

  double noise_amp = app->noise_amp;
  double mode = app->mode;

  double B0 = app->B0;

  double l = app->l;

  double ix = app->ix;
  double iy = app->iy;
  
  double Lx = app->Lx;
  double Ly = app->Ly;

  double Bz_noise = noise_amp * cos(iy * pi * y / Ly) * sin(ix * 2.0 * pi * x / Lx) / mode;
  
  double Ex = 0.0; // Total electric field (x-direction).
  double Ey = 0.0; // Total electric field (y-direction).
  double Ez = 0.0; // Total electric field (z-direction).

  double Bx = 0.0; // Total magnetic field (x-direction).
  double By = 0.0; // Total magnetic field (y-direction).
  double Bz = -B0 * tanh(y / l) + Bz_noise; // Total magnetic field (z-direction).

  // Set electric field.
  fout[0] = Ex; fout[1] = Ey; fout[2] = Ez;
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

  struct lhdi_ctx ctx = create_ctx(); // Context for initialization functions.

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
  struct gkyl_wv_eqn *elc_ten_moment = gkyl_wv_ten_moment_new(ctx.k0_elc, false, ctx.use_nn_closure, ctx.poly_order, ann[0], app_args.use_gpu);
  struct gkyl_wv_eqn *ion_ten_moment = gkyl_wv_ten_moment_new(ctx.k0_ion, false, ctx.use_nn_closure, ctx.poly_order, ann[1], app_args.use_gpu);

  struct gkyl_moment_species elc = {
    .name = "elc",
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .equation = elc_ten_moment,
    
    .init = evalElcInit,
    .ctx = &ctx,

    .bcy = { GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT },
  };

  struct gkyl_moment_species ion = {
    .name = "ion",
    .charge = ctx.charge_ion, .mass = ctx.mass_ion,
    .equation = ion_ten_moment,
    
    .init = evalIonInit,
    .ctx = &ctx,

    .bcy = { GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT },    
  };

  // Field.
  struct gkyl_moment_field field = {
    .epsilon0 = ctx.epsilon0, .mu0 = ctx.mu0,
    .mag_error_speed_fact = 1.0,
    
    .init = evalFieldInit,
    .ctx = &ctx,
    
    .bcy = { GKYL_FIELD_PEC_WALL, GKYL_FIELD_PEC_WALL },
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
    .lower = { -0.5 * ctx.Lx, -0.5 * ctx.Ly},
    .upper = { 0.5 * ctx.Lx, 0.5 * ctx.Ly},
    .cells = { NX, NY },

    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },
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
