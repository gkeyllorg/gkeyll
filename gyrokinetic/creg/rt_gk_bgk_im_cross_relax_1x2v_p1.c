#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_gyrokinetic_run.h>
#include <gkyl_util.h>

#include <rt_arg_parse.h>

struct bgk_cross_relax_ctx
{ 
  int cdim, vdim; // Dimensionality.

  // Physical constants (using non-normalized physical units).
  double epsilon0; // Permittivity of free space.
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.
  double mass_ion; // Proton mass.
  double charge_ion; // Proton charge.

  double T_par_elc; // Parallel electron temperature.
  double T_par_ion; // Parallel ion temperature.
  double alpha; // Ratio of perpendicular to parallel temperatures.

  double n0; // Reference number density (1 / m^3).
  double B0; // Reference magnetic field strength (Tesla).

  // Derived physical quantities (using non-normalized physical units).
  double T_perp_elc; // Perpendicular electron temperature.
  double T_perp_ion; // Perpendicular ion temperature.

  double Te; // Electron temperature.
  double Ti; // Ion temperature.

  double log_lambda_elc; // Electron Coulomb logarithm.
  double log_lambda_ion; // Ion Coulomb logarithm.
  double nu_elc; // Electron collision frequency.
  double nu_ion; // Ion collision frequency.

  double vte_par; // Parallel electron thermal velocity.
  double vti_par; // Parallel ion thermal velocity.
  double vte; // Electron thermal velocity.
  double vti; // Ion thermal velocity.

  double upar_elc; // Parallel electron velocity.
  double upar_ion; // Parallel ion velocity.

  // Simulation parameters.
  int Nz; // Cell count (configuration space: z-direction).
  int Nvpar; // Cell count (velocity space: parallel velocity direction).
  int Nmu; // Cell count (velocity space: magnetic moment direction).
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  double Lz; // Domain size (configuration space: z-direction).
  double vpar_max_elc; // Domain boundary (electron velocity space: parallel velocity direction).
  double mu_max_elc; // Domain boundary (electron velocity space: magnetic moment direction).
  double vpar_max_ion; // Domain boundary (ion velocity space: parallel velocity direction).
  double mu_max_ion; // Domain boundary (ion velocity space: magnetic moment direction).
  int poly_order; // Polynomial order.
  double cfl_frac; // CFL coefficient.

  double t_end; // End time.
  int num_frames; // Number of output frames.
  double write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct bgk_cross_relax_ctx
create_ctx(void)
{
  int cdim = 1, vdim = 2; // Dimensionality.

  // Physical constants (using non-normalized physical units).
  double epsilon0 = GKYL_EPSILON0; // Permittivity of free space.
  double mass_elc = GKYL_ELECTRON_MASS; // Electron mass.
  double mass_ion = 2.014 * GKYL_PROTON_MASS; // Proton mass.
  double charge_elc = -GKYL_ELEMENTARY_CHARGE; // Electron charge.
  double charge_ion = GKYL_ELEMENTARY_CHARGE; // Proton charge.

  double T_par_elc = 300.0 * GKYL_ELEMENTARY_CHARGE; // Parallel electron temperature.
  double T_par_ion = 200.0 * GKYL_ELEMENTARY_CHARGE; // Parallel ion temperature.
  double alpha = 1.3; // Ration of perpendicular to parallel temperatures.

  double n0 = 7.0e19; //  Reference number density (1 / m^3).
  double B0 = 1.0; // Reference magnetic field strength (Tesla).

  // Derived physical quantities (using non-normalized physical units).
  double T_perp_elc = alpha * T_par_elc; // Perpendicular electron temperature.
  double T_perp_ion = alpha * T_par_ion; // Perpendicular ion temperature.

  double Te = (T_par_elc + (2.0 * T_perp_elc)) / 3.0; // Electron temperature.
  double Ti = (T_par_ion + (2.0 * T_perp_ion)) / 3.0; // Ion temperature.

  double log_lambda_elc = 6.6 - 0.5 * log(n0 / 1.0e20) + 1.5 * log(T_par_elc / charge_ion); // Electron Coulomb logarithm.
  double log_lambda_ion = 6.6 - 0.5 * log(n0 / 1.0e20) + 1.5 * log(T_par_ion / charge_ion); // Ion Coulomb logarithm.
  double nu_elc = log_lambda_elc * pow(charge_ion, 4.0) * n0 /
    (6.0 * sqrt(2.0) * pow(M_PI, 3.0 / 2.0) * pow(epsilon0, 2.0) * sqrt(mass_elc) * pow(T_par_elc, 3.0 / 2.0)); // Electron collision frequency.
  double nu_ion = log_lambda_ion * pow(charge_ion, 4.0) * n0 /
    (12.0 * pow(M_PI, 3.0 / 2.0) * pow(epsilon0, 2.0) * sqrt(mass_ion) * pow(T_par_ion, 3.0 / 2.0)); // Ion collision frequency.
  
  double vte_par = sqrt(T_par_elc / mass_elc); // Parallel electron thermal velocity.
  double vti_par = sqrt(T_par_ion / mass_ion); // Parallel ion thermal velocity.
  double vte = sqrt(Te / mass_elc); // Electron thermal velocity.
  double vti = sqrt(Ti / mass_ion); // Ion thermal velocity.

  double upar_elc = 0.5 * sqrt(mass_elc / mass_ion) * vte; // Parallel electron velocity.
  double upar_ion = 50.0 * (mass_elc / mass_ion) * vti; // Parallel ion velocity.

  // Simulation parameters.
  int Nz = 1; // Cell count (configuration space: z-direction).
  int Nvpar = 16; // Cell count (velocity space: parallel velocity direction).
  int Nmu = 16; // Cell count (velocity space: magnetic moment direction).
  double Lz = 4.0; // Domain size (configuration space: z-direction).
  double vpar_max_elc = 5.0 * vte_par; // Domain boundary (electron velocity space: parallel velocity direction).
  double mu_max_elc = mass_elc * pow(5.0 * vte_par, 2.0) / (2.0 * B0); // Domain boundary (electron velocity space: magnetic moment direction).
  double vpar_max_ion = 5.0 * vti_par; // Domain boundary (ion velocity space: parallel velocity direction).
  double mu_max_ion = mass_ion * pow(5.0 * vti_par, 2.0) / (2.0 * B0); // Domain boundary (ion velocity space: magnetic moment direction).
  int poly_order = 1; // Polynomial order.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 0.1 / nu_ion; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  double write_phase_freq = 0.2; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.
  
  struct bgk_cross_relax_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .epsilon0 = epsilon0,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .mass_ion = mass_ion,
    .charge_ion = charge_ion,
    .T_par_elc = T_par_elc,
    .T_par_ion = T_par_ion,
    .n0 = n0,
    .B0 = B0,
    .T_perp_elc = T_perp_elc,
    .T_perp_ion = T_perp_ion,
    .Te = Te,
    .Ti = Ti,
    .log_lambda_elc = log_lambda_elc,
    .nu_elc = nu_elc,
    .log_lambda_ion = log_lambda_ion,
    .nu_ion = nu_ion,
    .vte_par = vte_par,
    .vti_par = vti_par,
    .vte = vte,
    .vti = vti,
    .upar_elc = upar_elc,
    .upar_ion = upar_ion,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nz, Nvpar, Nmu},
    .Lz = Lz,
    .vpar_max_elc = vpar_max_elc,
    .mu_max_elc = mu_max_elc,
    .vpar_max_ion = vpar_max_ion,
    .mu_max_ion = mu_max_ion,
    .poly_order = poly_order,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .num_frames = num_frames,
    .write_phase_freq = write_phase_freq,
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };

  return ctx;
}

void
evalElcDensityInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct bgk_cross_relax_ctx *app = ctx;

  double n0 = app->n0;

  // Set electron total number density.
  fout[0] = n0;
}

void
evalElcTparInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct bgk_cross_relax_ctx *app = ctx;

  double T_par_elc = app->T_par_elc;

  // Set electron parallel temperature.
  fout[0] = T_par_elc;
}

void
evalElcTperpInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct bgk_cross_relax_ctx *app = ctx;

  double T_perp_elc = app->T_perp_elc;

  // Set electron perpendicular temperature.
  fout[0] = T_perp_elc;
}

void
evalElcUparInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct bgk_cross_relax_ctx *app = ctx;

  double upar_elc = app->upar_elc;

  // Set electron parallel velocity.
  fout[0] = upar_elc;
}

void
evalIonDensityInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct bgk_cross_relax_ctx *app = ctx;

  double n0 = app->n0;

  // Set ion total number density.
  fout[0] = n0;
}

void
evalIonTparInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct bgk_cross_relax_ctx *app = ctx;

  double T_par_ion = app->T_par_ion;

  // Set ion parallel temperature.
  fout[0] = T_par_ion;
}

void
evalIonTperpInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct bgk_cross_relax_ctx *app = ctx;

  double T_perp_ion = app->T_perp_ion;

  // Set ion perpendicular temperature.
  fout[0] = T_perp_ion;
}

void
evalIonUparInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct bgk_cross_relax_ctx *app = ctx;
  
  double upar_ion = app->upar_ion;

  // Set ion parallel velocity.
  fout[0] = upar_ion;
}

static inline void
mapc2p(double t, const double* GKYL_RESTRICT zc, double* GKYL_RESTRICT xp, void* ctx)
{
  // Set physical coordinates (X, Y, Z) from computational coordinates (x, y, z).
  xp[0] = zc[0]; xp[1] = zc[1]; xp[2] = zc[2];
}

void
bfield_func(double t, const double* GKYL_RESTRICT zc, double* GKYL_RESTRICT fout, void* ctx)
{
  struct bgk_cross_relax_ctx *app = ctx;

  double B0 = app->B0;

  // zc are computational coords. 
  // Set Cartesian components of magnetic field.
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = app->B0;
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) MPI_Init(&argc, &argv);
#endif

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct bgk_cross_relax_ctx ctx = create_ctx(); // Context for init functions.

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  // Electrons.
  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .vdim = ctx.vdim,
    .lower = { -ctx.vpar_max_elc, 0.0 },
    .upper = { ctx.vpar_max_elc, ctx.mu_max_elc },
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_BIMAXWELLIAN,
      .density = evalElcDensityInit,
      .ctx_density = &ctx,
      .temppar = evalElcTparInit,
      .ctx_temppar = &ctx,
      .tempperp = evalElcTperpInit,
      .ctx_tempperp = &ctx,
      .upar = evalElcUparInit,
      .ctx_upar = &ctx,
      .correct_all_moms = true,
    },

    .correct = {
      .correct_all_moms = true, 
      .use_last_converged = true, 
      .iter_eps = 1e-12,
      .max_iter = 10,
    },   

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES_NO_BY,
    },

    .collisions =  {
      .collision_id = GKYL_BGK_COLLISIONS,
      .den_ref = ctx.n0,
      .temp_ref = ctx.Te,
      .num_cross_collisions = 1,
      .collide_with = { "ion" },
      .is_implicit = true,
    },

    .num_diag_moments = 6,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_BIMAXWELLIAN },
  };

  // Ions.
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.charge_ion, .mass = ctx.mass_ion,
    .vdim = ctx.vdim,
    .lower = { -ctx.vpar_max_ion, 0.0 },
    .upper = { ctx.vpar_max_ion, ctx.mu_max_ion },
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_BIMAXWELLIAN,
      .density = evalIonDensityInit,
      .ctx_density = &ctx,
      .temppar = evalIonTparInit,
      .ctx_temppar = &ctx,
      .tempperp = evalIonTperpInit,
      .ctx_tempperp = &ctx,
      .upar = evalIonUparInit,
      .ctx_upar = &ctx,
      .correct_all_moms = true,
    },

    .correct = {
      .correct_all_moms = true, 
      .use_last_converged = true, 
      .iter_eps = 1e-12,
      .max_iter = 10,
    },       

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES_NO_BY,
    },

    .collisions =  {
      .collision_id = GKYL_BGK_COLLISIONS,
      .den_ref = ctx.n0,
      .temp_ref = ctx.Ti,
      .num_cross_collisions = 1,
      .collide_with = { "elc" },
      .is_implicit = true,
    },
    
    .num_diag_moments = 6,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_BIMAXWELLIAN },
  };

  // Field.
  struct gkyl_gyrokinetic_field field = {
    .zero_init_field = true, // Don't compute the field at t = 0.
    .is_static = true, // Don't evolve the field in time.
  };

  // Gyrokinetic app.
  struct gkyl_gk app_inp = {

    .cdim = ctx.cdim,
    .lower = { -0.5 * ctx.Lz },
    .upper = { 0.5 * ctx.Lz },
    .cells = { cells_x[0] },

    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,
    .cfl_frac = ctx.cfl_frac,

    .geometry = {
      .geometry_id = GKYL_GEOMETRY_MAPC2P,
      .world = { 0.0, 0.0 },

      .mapc2p = mapc2p,
      .c2p_ctx = &ctx,
      .bfield_func = bfield_func,
      .bfield_ctx = &ctx
    },

    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },

    .num_species = 2,
    .species = { elc, ion },

    .field = field,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0] },
      .comm = comm,
    },
  };


  // Set app output name from the executable name (argv[0]).
  snprintf(app_inp.name, sizeof(app_inp.name), "%s", app_args.app_name);
  struct gkyl_gyrokinetic_run_inp run_inp = {
    .app_inp = app_inp,
    .time_stepping = {
      .t_end = ctx.t_end,
      .num_frames = ctx.num_frames,
      .write_phase_freq = ctx.write_phase_freq,
      .int_diag_calc_num = ctx.int_diag_calc_num,
      .dt_failure_tol = ctx.dt_failure_tol,
      .num_failures_max = ctx.num_failures_max,
      .is_restart = app_args.is_restart,
      .restart_frame = app_args.restart_frame,
      .num_steps = app_args.num_steps,
    },
  };

  gkyl_gyrokinetic_run_simulation(&run_inp);

  gkyl_gyrokinetic_comms_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif

  return 0;
}
