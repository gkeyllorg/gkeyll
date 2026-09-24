#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_gyrokinetic_run.h>

#include <rt_arg_parse.h>

struct rad_ctx
{
  int cdim, vdim; // Dimensionality.
  
  // Physical constants (using non-normalized physical units).
  double epsilon0; // Permittivity of free space.
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.
  double mass_ion; // Proton mass.
  double charge_ion; // Proton charge.

  double Te; // Electron temperature.
  double Ti; // Ion temperature.
  double B0; // Reference magnetic field strength (Tesla).
  double n0; // Reference number density (1 / m^3).

  double nu_frac; // Collision frequency fraction.

  double k_perp_rho_s; // Product of perpendicular wavenumber and ion-sound gyroradius.

  double log_lambda_elc; // Electron Coulomb logarithm.
  double log_lambda_ion; // Ion Coulomb logarithm.
  double nu_elc; // Electron collision frequency.
  double nu_ion; // Ion collision frequency.
  double nu_elc_ion; // Electron-ion collision frequency.
  double nu_ion_elc; // Ion-electron collision frequency.

  double c_s; // Sound speed.
  double vte; // Electron thermal velocity.
  double vti; // Ion thermal velocity.
  double omega_ci; // Ion cyclotron frequency.
  double rho_s; // Ion-sound gyroradius.

  double k_perp; // Perpendicular wavenumber (for Poisson solver).

  // Simulation parameters.
  int Nz; // Cell count (configuration space: z-direction).
  int Nvpar; // Cell count (velocity space: parallel velocity direction).
  int Nmu; // Cell count (velocity space: magnetic moment direction).
  int Nvneut; // Cell count (velocity space: each direction for neutral)
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  double Lz; // Domain size (configuration space: z-direction).
  double vpar_max_elc; // Domain boundary (electron velocity space: parallel velocity direction).
  double mu_max_elc; // Domain boundary (electron velocity space: magnetic moment direction).
  double vpar_max_ion; // Domain boundary (ion velocity space: parallel velocity direction).
  double mu_max_ion; // Domain boundary (ion velocity space: magnetic moment direction).

  double t_end; // End time.
  int num_frames; // Number of output frames.
  double write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct rad_ctx
create_ctx(void)
{
  int cdim = 1, vdim = 2; // Dimensionality.

  // Physical constants (using non-normalized physical units).
  double epsilon0 = GKYL_EPSILON0; // Permittivity of free space.
  double mass_elc = GKYL_ELECTRON_MASS; // Electron mass.
  double charge_elc = -GKYL_ELEMENTARY_CHARGE; // Electron charge.
  double mass_ion = GKYL_PROTON_MASS; // Proton mass.
  double charge_ion = GKYL_ELEMENTARY_CHARGE; // Proton charge.

  double Te = 0.72 * GKYL_ELEMENTARY_CHARGE; // Electron temperature.
  double Ti = 30.0 * GKYL_ELEMENTARY_CHARGE; // Ion temperature.
  double B0 = 1.0; // Reference magnetic field strength (Tesla).
  double n0 = 1.0e19; //  Reference number density (1 / m^3).

  double nu_frac = 0.1; // Collision frequency fraction.

  double k_perp_rho_s = 0.1; // Product of perpendicular wavenumber and ion-sound gyroradius.

  // Coulomb logarithms.
  double log_lambda_elc = 6.6 - 0.5 * log(n0 / 1.0e20) + 1.5 * log(Te / charge_ion);
  double log_lambda_ion = 6.6 - 0.5 * log(n0 / 1.0e20) + 1.5 * log(Ti / charge_ion);

  // Collision frequencies.
  double nu_elc = nu_frac * log_lambda_elc * pow(charge_ion,4) * n0 /
    (6.0 * sqrt(2.0) * pow(M_PI,3.0/2.0) * pow(epsilon0,2) * sqrt(mass_elc) * pow(Te,3.0/2.0));
  double nu_ion = nu_frac * log_lambda_ion * pow(charge_ion,4) * n0 /
    (12.0 * pow(M_PI,3.0/2.0) * pow(epsilon0,2) * sqrt(mass_ion) * pow(Ti,3.0/2.0));
  double nu_elc_ion = nu_elc*sqrt(2.0);
  double nu_ion_elc = nu_elc_ion*(mass_elc/mass_ion);
  
  double c_s = sqrt(Te / mass_ion); // Sound speed.
  double vte = sqrt(Te / mass_elc); // Electron thermal velocity.
  double vti = sqrt(Ti / mass_ion); // Ion thermal velocity.
  double omega_ci = fabs(charge_ion * B0 / mass_ion); // Ion cyclotron frequency.
  double rho_s = c_s / omega_ci; // Ion-sound gyroradius.

  double k_perp = k_perp_rho_s / rho_s; // Perpendicular wavenumber (for Poisson solver).

  // Simulation parameters.
  int Nz = 4; // Cell count (configuration space: z-direction).
  int Nvpar = 64; // Cell count (velocity space: parallel velocity direction).
  int Nmu = 32; // Cell count (velocity space: magnetic moment direction).
  int Nvneut = 8; // Cell count (velocity space: each direction for neutral)
  double Lz = 2 * 5.596343e-02; // Domain size (configuration space: z-direction).
  double vpar_max_elc = 8.0 * vte; // Domain size (electron velocity space: parallel velocity direction).
  double mu_max_elc = 0.75 * mass_elc * (4.0 * vte) * (4.0 * vte) / (2.0 * B0); // Domain size (electron velocity space: magnetic moment direction).
  double vpar_max_ion = 8.0 * vti; // Domain size (ion velocity space: parallel velocity direction).
  double mu_max_ion = 0.75 * mass_ion * (4.0 * vti) * (4.0 * vti) / (2.0 * B0); // Domain size (ion velocity space: magnetic moment direction).

  double t_end = 1.0e-10; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  double write_phase_freq = 0.2; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.
  
  struct rad_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .epsilon0 = epsilon0,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .mass_ion = mass_ion,
    .charge_ion = charge_ion,
    .Te = Te,
    .Ti = Ti,
    .B0 = B0,
    .n0 = n0,
    .nu_frac = nu_frac,
    .k_perp_rho_s = k_perp_rho_s,
    .log_lambda_elc = log_lambda_elc,
    .nu_elc = nu_elc,
    .log_lambda_ion = log_lambda_ion,
    .nu_ion = nu_ion,
    .nu_elc_ion = nu_elc_ion,
    .nu_ion_elc = nu_ion_elc,
    .c_s = c_s,
    .vte = vte,
    .vti = vti,
    .omega_ci = omega_ci,
    .rho_s = rho_s,
    .k_perp = k_perp,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .Nvneut = Nvneut,
    .cells = {Nz, Nvpar, Nmu, Nvneut},
    .Lz = Lz,
    .vpar_max_elc = vpar_max_elc,
    .mu_max_elc = mu_max_elc,
    .vpar_max_ion = vpar_max_ion,
    .mu_max_ion = mu_max_ion,
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
evalDensityInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct rad_ctx *app = ctx;

  double n0 = app -> n0;
  double x = xn[0];
  double Lz = app -> Lz;

  // Set number density.
  fout[0] = fmax( n0 *  pow(x+1.0,120.0), n0*1e-6);
}

void
evalElcDensityInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct rad_ctx *app = ctx;

  double n0 = app -> n0;
  double x = xn[0];
  double Lz = app -> Lz;

  // Set number density.
  fout[0] = fmax( n0 *  pow(x+1.0,120.0), n0*1e-6);
}

void
evalUparInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Set parallel velocity.
  fout[0] = 0.0;
}

void
evalTempElcInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct rad_ctx *app = ctx;

  double Te = app->Te;

  // Set electron temperature.
  fout[0] = Te;
}

void
evalTempIonInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct rad_ctx *app = ctx;

  double Ti = app->Ti;

  // Set ion temperature.
  fout[0] = Ti;
}

void
evalNuElc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct rad_ctx *app = ctx;

  double nu_elc = app->nu_elc;

  // Set electron collision frequency.
  fout[0] = nu_elc;
}

void
evalNuIon(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct rad_ctx *app = ctx;

  double nu_ion = app->nu_ion;

  // Set ion collision frequency.
  fout[0] = app->nu_ion;
}

void
evalNuElcIon(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct rad_ctx *app = ctx;

  double nu_elc_ion = app->nu_elc_ion;

  // Set electron-ion collision frequency.
  fout[0] = nu_elc_ion;
}

void
evalNuIonElc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct rad_ctx *app = ctx;

  double nu_ion_elc = app->nu_ion_elc;

  // Set ion-electron collision frequency.
  fout[0] = app->nu_ion_elc;
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
  struct rad_ctx *app = ctx;

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

  struct rad_ctx ctx = create_ctx(); // Context for init functions.

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  // Electron species.
  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .vdim = ctx.vdim,
    .lower = { -0.5 * ctx.vpar_max_elc, 0.0 },
    .upper = { 0.5 * ctx.vpar_max_elc, ctx.mu_max_elc },
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0/2,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .density = evalElcDensityInit,
      .ctx_density = &ctx,
      .upar = evalUparInit,
      .ctx_upar = &ctx,
      .temp = evalTempElcInit,
      .ctx_temp = &ctx,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .self_nu = evalNuElc,
      .self_nu_ctx = &ctx,
      .num_cross_collisions = 2,
      .collide_with = { "ion", "elc2" },
      .cross_nu = { evalNuElcIon, evalNuElc, },
      .cross_nu_ctx = &ctx,
      .den_ref = ctx.n0/2,
      .temp_ref = ctx.Te,
    },

    //Emissivity for all three of these should be 0
    .radiation = {
      .radiation_id = GKYL_GK_RADIATION,
      .num_cross_collisions = 3, 
      .collide_with = { "ion", "test_sp_1", "test_sp_2"},
      .atomic_Z = {1, 3, 3},
      .charge_state = {0, 0, 0},
      .num_of_densities = {1, 1, 5},
      .reference_ne = ctx.n0/10,
      .min_ne = ctx.n0/100,
      .max_ne = ctx.n0,
      .te_min_model = GKYL_CONST_TE,
      .Te_min = 12.0 * GKYL_ELEMENTARY_CHARGE,
    },
    
    .num_diag_moments = 7,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP },
  };

  // 2nd Electron species.
  struct gkyl_gyrokinetic_species elc2 = {
    .name = "elc2",
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .vdim = ctx.vdim,
    .lower = { -0.5 * ctx.vpar_max_elc, 0.0 },
    .upper = { 0.5 * ctx.vpar_max_elc, ctx.mu_max_elc },
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0/2,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .density = evalElcDensityInit,
      .ctx_density = &ctx,
      .upar = evalUparInit,
      .ctx_upar = &ctx,
      .temp = evalTempElcInit,
      .ctx_temp = &ctx,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .self_nu = evalNuElc,
      .self_nu_ctx = &ctx,
      .num_cross_collisions = 2,
      .collide_with = { "ion", "elc" },
      .cross_nu = { evalNuElcIon, evalNuElc, },
      .cross_nu_ctx = &ctx,
      .den_ref = ctx.n0/2,
      .temp_ref = ctx.Te,
    },

    //Emissivity of ion and test_sp_1 should be 0
    .radiation = {
      .radiation_id = GKYL_GK_RADIATION, 
      .num_cross_collisions = 5, 
      .collide_with = { "ion", "test_sp_1", "test_sp_2", "test_sp_3", "test_sp_4"},
      .atomic_Z = {1, 3, 3, 3, 3},
      .charge_state = {0, 0, 0, 0, 0},
      .num_of_densities = {1, 1, 3, 8, 15},
      .reference_ne = ctx.n0,
      .min_ne = ctx.n0*1.01,
      .max_ne = ctx.n0*10000,
      },
    
    .num_diag_moments = 7,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP },
  };
  
  // Ion species.
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.charge_ion, .mass = ctx.mass_ion,
    .vdim = ctx.vdim,
    .lower = { -ctx.vpar_max_ion, 0.0 },
    .upper = { ctx.vpar_max_ion, ctx.mu_max_ion },
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .density = evalDensityInit,
      .ctx_density = &ctx,
      .upar= evalUparInit,
      .ctx_upar = &ctx,
      .temp = evalTempIonInit,
      .ctx_temp = &ctx,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .self_nu = evalNuElc,
      .self_nu_ctx = &ctx,
      .num_cross_collisions = 2,
      .collide_with = { "elc", "elc2", },
      .cross_nu = { evalNuIonElc, evalNuIonElc, },
      .cross_nu_ctx = &ctx,
      .den_ref = ctx.n0,
      .temp_ref = ctx.Ti,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 0,
    },

    .num_diag_moments = 7,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP },
  };

  // D0
  struct gkyl_gyrokinetic_neut_species test_sp_1 = {
    .name = "test_sp_1", .mass = ctx.mass_ion,
    .vdim = ctx.vdim+1,
    .lower = { -ctx.vpar_max_ion, -ctx.vpar_max_ion, -ctx.vpar_max_ion},
    .upper = { ctx.vpar_max_ion, ctx.vpar_max_ion, ctx.vpar_max_ion },
    .cells = { ctx.Nvneut, ctx.Nvneut, ctx.Nvneut},
    .is_static = true,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = evalDensityInit,
      .ctx_upar = &ctx,
      .udrift = evalUparInit,
      .ctx_temp = &ctx,
      .temp = evalTempIonInit,
    },

    .num_diag_moments = 3,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2},
  };

  // Second D0
  struct gkyl_gyrokinetic_neut_species test_sp_2 = {
    .name = "test_sp_2", .mass = ctx.mass_ion,
    .vdim = ctx.vdim+1,
    .lower = { -ctx.vpar_max_ion, -ctx.vpar_max_ion, -ctx.vpar_max_ion},
    .upper = { ctx.vpar_max_ion, ctx.vpar_max_ion, ctx.vpar_max_ion },
    .cells = { ctx.Nvneut, ctx.Nvneut, ctx.Nvneut},
    .is_static = true,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = evalDensityInit,
      .ctx_upar = &ctx,
      .udrift = evalUparInit,
      .ctx_temp = &ctx,
      .temp = evalTempIonInit,
    },

    .num_diag_moments = 3,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2},
  };
  // Third D0
  struct gkyl_gyrokinetic_neut_species test_sp_3 = {
    .name = "test_sp_3", .mass = ctx.mass_ion,
    .vdim = ctx.vdim+1,
    .lower = { -ctx.vpar_max_ion, -ctx.vpar_max_ion, -ctx.vpar_max_ion},
    .upper = { ctx.vpar_max_ion, ctx.vpar_max_ion, ctx.vpar_max_ion },
    .cells = { ctx.Nvneut, ctx.Nvneut, ctx.Nvneut},
    .is_static = true,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = evalDensityInit,
      .ctx_upar = &ctx,
      .udrift = evalUparInit,
      .ctx_temp = &ctx,
      .temp = evalTempIonInit,
    },

    .num_diag_moments = 3,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2},
  };
  // Fourth D0
  struct gkyl_gyrokinetic_neut_species test_sp_4 = {
    .name = "test_sp_4", .mass = ctx.mass_ion,
    .vdim = ctx.vdim+1,
    .lower = { -ctx.vpar_max_ion, -ctx.vpar_max_ion, -ctx.vpar_max_ion},
    .upper = { ctx.vpar_max_ion, ctx.vpar_max_ion, ctx.vpar_max_ion },
    .cells = { ctx.Nvneut, ctx.Nvneut, ctx.Nvneut},
    .is_static = true,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = evalDensityInit,
      .ctx_upar = &ctx,
      .udrift = evalUparInit,
      .ctx_temp = &ctx,
      .temp = evalTempIonInit,
    },

    .num_diag_moments = 3,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2},
  };
  
  // Field.
  struct gkyl_gyrokinetic_field field = {
    .kperpSq = ctx.k_perp * ctx.k_perp,
    .zero_init_field = true, // Don't compute the field at t=0.
    .is_static = true, // Don't update the field in time.
  };

  // GK app.
  struct gkyl_gk app_inp = {

    .cdim = ctx.cdim,
    .lower = { -0.5 * ctx.Lz },
    .upper = { 0.5 * ctx.Lz },
    .cells = { cells_x[0] },
    .poly_order = 1,
    .basis_type = app_args.basis_type,

    .geometry = {
      .geometry_id = GKYL_GEOMETRY_MAPC2P,
      .world = { 0.0, 0.0 },
      .mapc2p = mapc2p,
      .c2p_ctx = &ctx,
      .bfield_func = bfield_func,
      .bfield_ctx = &ctx,
    },

    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },

    .num_species = 3,
    .species = { elc, elc2, ion },

    .num_neut_species = 4,
    .neut_species = {test_sp_1, test_sp_2, test_sp_3, test_sp_4},

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
