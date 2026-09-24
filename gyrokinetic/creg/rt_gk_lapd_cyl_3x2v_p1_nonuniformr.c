#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_gyrokinetic_run.h>

#include <rt_arg_parse.h>

struct lapd_cyl_ctx
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

  double Te_src; // Source electron temperature.
  double r_src; // Source radial extent.
  double L_src; // Source length.
  double S0; // Source reference number density.
  double floor_src; // Minimum source intensity;

  // Simulation parameters.
  int Nr; // Cell count (configuration space: radial direction).
  int Ntheta; // Cell count (configuration space: angular direction).
  int Nz; // Cell count (configuration space: z-direction).
  int Nvpar; // Cell count (velocity space: parallel velocity direction).
  int Nmu; // Cell count (velocity space: magnetic moment direction).
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  double Lr; // Domain size (configuration space: radial direction).
  double Ltheta; // Domain size (configuration space: angular direction).
  double Lz; // Domain size (configuration space: z-direction).
  double L_perp; // Perpendicular length of domain.
  double vpar_max_elc; // Domain boundary (electron velocity space: parallel velocity direction).
  double mu_max_elc; // Domain boundary (electron velocity space: magnetic moment direction).
  double vpar_max_ion; // Domain boundary (ion velocity space: parallel velocity direction).
  double mu_max_ion; // Domain boundary (ion velocity space: magnetic moment direction).

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  double write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct lapd_cyl_ctx
create_ctx(void)
{
  int cdim = 3, vdim = 2; // Dimensionality.

  // Physical constants (using non-normalized physical units).
  double epsilon0 = GKYL_EPSILON0; // Permittivity of free space.
  double mass_elc = GKYL_PROTON_MASS * 3.973 / 400.0; // Electron mass.
  double charge_elc = -GKYL_ELEMENTARY_CHARGE; // Electron charge.
  double mass_ion = 3.973 * GKYL_PROTON_MASS; // Proton mass.
  double charge_ion = GKYL_ELEMENTARY_CHARGE; // Proton charge.

  double Te = 6.0 * GKYL_ELEMENTARY_CHARGE; // Electron temperature.
  double Ti = 1.0 * GKYL_ELEMENTARY_CHARGE; // Ion temperature.
  double B0 = 0.0398; // Reference magnetic field strength (Tesla).
  double n0 = 2.0e18; //  Reference number density (1 / m^3).

  double nu_frac = 0.1; // Collision frequency fraction.

  // Coulomb logarithms.
  double log_lambda_elc = 6.6 - 0.5 * log(n0 / 1.0e20) + 1.5 * log(Te / charge_ion);
  double log_lambda_ion = 6.6 - 0.5 * log(n0 / 1.0e20) + 1.5 * log(Ti / charge_ion);

  // Collision frequencies.
  double nu_elc = nu_frac * log_lambda_elc * pow(charge_ion,4) * n0 /
    (6.0 * sqrt(2.0) * pow(M_PI,3.0/2.0) * pow(epsilon0,2) * sqrt(mass_elc) * pow(Te,3.0/2.0));
  double nu_ion = nu_frac * log_lambda_ion * pow(charge_ion,4) * n0 /
    (12.0 * pow(M_PI,3.0/2.0) * pow(epsilon0,2) * sqrt(mass_ion) * pow(Ti,3.0/2.0));
  double nu_elc_ion = nu_elc*sqrt(2.0);
  double nu_ion_elc = nu_ion_elc*(mass_elc/mass_ion);
  
  double c_s = sqrt(Te / mass_ion); // Sound speed.
  double vte = sqrt(Te / mass_elc); // Electron thermal velocity.
  double vti = sqrt(Ti / mass_ion); // Ion thermal velocity.
  double omega_ci = fabs(charge_ion * B0 / mass_ion); // Ion cyclotron frequency.
  double rho_s = c_s / omega_ci; // Ion-sound gyroradius.

  double Te_src = 6.8 * GKYL_ELEMENTARY_CHARGE; // Source electron temperature.
  double r_src = 20.0 * rho_s; // Source radial extent.
  double L_src = 0.5 * rho_s; // Source length.
  double S0 = 1.08 * n0 * c_s / (36.0 * 40.0 * rho_s); // Source reference number density.
  double floor_src = 0.01; // Minimum source intensity.

  // Simulation parameters.
  int Nr = 10; // Cell count (configuration space: radial direction).
  int Ntheta = 10; // Cell count (configuration space: angular direction).
  int Nz = 8; // Cell count (configuration space: z-direction).
  int Nvpar = 8; // Cell count (velocity space: parallel velocity direction).
  int Nmu = 4; // Cell count (velocity space: magnetic moment direction).
  double Lr = 47.5 * rho_s; // Domain size (configuration space: radial direction).
  double Ltheta = 2.0 * M_PI; // Domain size (configuration space: angular direction).
  double Lz = 36.0 * 40.0 * rho_s; // Domain size (configuration space: z-direction).
  double L_perp = 100.0 * rho_s; // Perpendicular length of domain.
  double vpar_max_elc = 4.0 * vte; // Domain boundary (electron velocity space: parallel velocity direction).
  double mu_max_elc = (3.0 / 2.0) * 0.5 * mass_elc * pow(4.0 * vte,2) / (2.0 * B0); // Domain boundary (electron velocity space: magnetic moment direction).
  double vpar_max_ion = 4.0 * vti; // Domain boundary (ion velocity space: parallel velocity direction).
  double mu_max_ion = (3.0 / 2.0) * 0.5 * mass_ion * pow(4.0 * vti,2) / (2.0 * B0); // Domain boundary (ion velocity space: magnetic moment direction).

  double t_end = 1.0e-6; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  double write_phase_freq = 0.2; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.
  
  struct lapd_cyl_ctx ctx = {
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
    .Te_src = Te_src,
    .r_src = r_src,
    .L_src = L_src,
    .S0 = S0,
    .floor_src = floor_src,
    .Nr = Nr,
    .Ntheta = Ntheta,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nr, Ntheta, Nz, Nvpar, Nmu},
    .Lr = Lr,
    .Ltheta = Ltheta,
    .Lz = Lz,
    .L_perp = L_perp,
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
  struct lapd_cyl_ctx *app = ctx;
  double r = xn[0];

  double n0 = app->n0;
  double L_perp = app->L_perp;

  pcg32_random_t rng = gkyl_pcg32_init(0);
  double perturb = 2.0e-3 * (1.0 - 0.5 * gkyl_pcg32_rand_double(&rng));

  double n = 0.0;

  if (r < 0.5 * L_perp) {
    n = ((1.0 - (1.0 / 20.0)) * pow(1.0 - (r / (0.5 * L_perp)) * (r / (0.5 * L_perp)), 3.0) + (1.0 / 20.0)) * n0 * (1.0 + perturb);
  }
  else {
    n = (1.0 / 20.0) * n0 * (1.0 + perturb);
  }

  // Set number density.
  fout[0] = n;
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
  struct lapd_cyl_ctx *app = ctx;
  double r = xn[0];

  double Te = app->Te;
  double L_perp = app->L_perp;

  double T = 0.0;

  if (r < 0.5 * L_perp) {
    T = ((1.0 - (1.0 / 5.0)) * pow(1.0 - (r / (0.5 * L_perp)) * (r / (0.5 * L_perp)), 3.0) + (1.0 / 5.0)) * Te;
  }
  else {
    T = (1.0 / 5.0) * Te;
  }

  // Set electron temperature.
  fout[0] = T;
}

void
evalTempIonInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct lapd_cyl_ctx *app = ctx;

  double Ti = app->Ti;

  // Set ion temperature.
  fout[0] = Ti;
}

void
evalSourceDensityInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct lapd_cyl_ctx *app = ctx;
  double r = xn[0];

  double r_src = app->r_src;
  double L_src = app->L_src;
  double S0 = app->S0;
  double floor_src = app->floor_src;

  // Set source number density.
  fout[0] = S0 * (floor_src + (1.0 - floor_src) * 0.5 * (1.0 - tanh((r - r_src) / L_src)));
}

void
evalSourceUparInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Set source parallel velocity.
  fout[0] = 0.0;
}

void
evalSourceTempElcInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct lapd_cyl_ctx *app = ctx;
  double r = xn[0];

  double Te_src = app->Te_src;
  double L_perp = app->L_perp;

  double T = 0.0;

  if (r < 0.5 * L_perp) {
    T = (1.0 - (1.0 / 2.5) * pow(1.0 - (r / (0.5 * L_perp)) * (r / (0.5 * L_perp)), 3.0) + (1.0 / 2.5)) * Te_src;
  }
  else {
    T = (1.0 / 2.5) * Te_src;
  }

  // Set electron source temperature.
  fout[0] = T;
}

void
evalSourceTempIonInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct lapd_cyl_ctx *app = ctx;

  double Ti = app->Ti;

  // Set ion source temperature.
  fout[0] = Ti;
}

void
evalNuElc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct lapd_cyl_ctx *app = ctx;

  double nu_elc = app->nu_elc;

  // Set electron collision frequency.
  fout[0] = nu_elc;
}

void
evalNuIon(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct lapd_cyl_ctx *app = ctx;

  double nu_ion = app->nu_ion;

  // Set ion collision frequency.
  fout[0] = nu_ion;
}

void
evalNuElcIon(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct lapd_cyl_ctx *app = ctx;

  double nu_elc_ion = app->nu_elc_ion;

  // Set electron-ion collision frequency.
  fout[0] = nu_elc_ion;
}

void
evalNuIonElc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct lapd_cyl_ctx *app = ctx;

  double nu_ion_elc = app->nu_ion_elc;

  // Set electron-ion collision frequency.
  fout[0] = nu_ion_elc;
}

void
mc2nu_r(double t, const double* GKYL_RESTRICT xc, double* GKYL_RESTRICT xnu, void* ctx)
{
  struct lapd_cyl_ctx *app = ctx;
  double r = xc[0];
  double r_min = 2.5 * app->rho_s;
  double r_max = 2.5 * app->rho_s + app->Lr;
  double poly_order = 1.4;
  xnu[0] = pow(r - r_min, poly_order)*pow(r_max - r_min, 1 - poly_order) + r_min;
}

static inline void
mapc2p(double t, const double* GKYL_RESTRICT zc, double* GKYL_RESTRICT xp, void* ctx)
{
  double r = zc[0], theta = zc[1], z = zc[2];

  // Set physical coordinates (X, Y, Z) from computational coordinates (r, theta, z).
  xp[0] = r * cos(theta);
  xp[1] = r * sin(theta);
  xp[2] = z;
}

void
bfield_func(double t, const double* GKYL_RESTRICT zc, double* GKYL_RESTRICT fout, void* ctx)
{
  struct lapd_cyl_ctx *app = ctx;

  double B0 = app->B0;

  // zc are computational coords. 
  // Set Cartesian components of magnetic field.
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = B0;
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

  struct lapd_cyl_ctx ctx = create_ctx(); // Context for init functions.

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  // Electron species.
  struct gkyl_gyrokinetic_projection elc_ic = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
    .density = evalDensityInit,
    .ctx_density = &ctx,
    .upar = evalUparInit,
    .ctx_upar = &ctx,
    .temp = evalTempElcInit,
    .ctx_temp = &ctx,
  };

  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .vdim = ctx.vdim,
    .lower = { -ctx.vpar_max_elc, 0.0 },
    .upper = { ctx.vpar_max_elc, ctx.mu_max_elc },
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    .projection = elc_ic,

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .self_nu = evalNuElc,
      .self_nu_ctx = &ctx,
      .num_cross_collisions = 1,
      .collide_with = { "ion" },
      .cross_nu = { evalNuElcIon, },
      .cross_nu_ctx = &ctx,
      .den_ref = ctx.n0,
      .temp_ref = ctx.Te, 
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .density = evalSourceDensityInit,
        .ctx_density = &ctx,
        .upar = evalSourceUparInit,
        .ctx_upar = &ctx,
        .temp = evalSourceTempElcInit,
        .ctx_temp = &ctx,
      }, 
    },
    
    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_FIXED_FUNC, .projection = elc_ic, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX, },
      { .dir = 2, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH, },
      { .dir = 2, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH, },
    },

    .num_diag_moments = 7,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP },
  };

  // Ion species.
  struct gkyl_gyrokinetic_projection ion_ic = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
    .density = evalDensityInit,
    .ctx_density = &ctx,
    .upar = evalUparInit,
    .ctx_upar = &ctx,
    .temp = evalTempIonInit,
    .ctx_temp = &ctx,
  };

  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.charge_ion, .mass = ctx.mass_ion,
    .vdim = ctx.vdim,
    .lower = { -ctx.vpar_max_ion, 0.0},
    .upper = { ctx.vpar_max_ion, ctx.mu_max_ion},
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    .projection = ion_ic,

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .self_nu = evalNuIon,
      .self_nu_ctx = &ctx,
      .num_cross_collisions = 1,
      .collide_with = { "elc" },
      .cross_nu = { evalNuIonElc, },
      .cross_nu_ctx = &ctx,
      .den_ref = ctx.n0,
      .temp_ref = ctx.Ti, 
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
        .density = evalSourceDensityInit,
        .ctx_density = &ctx,
        .upar = evalSourceUparInit,
        .ctx_upar = &ctx,
        .temp = evalSourceTempIonInit,
        .ctx_temp = &ctx,
      }, 
    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_FIXED_FUNC, .projection = ion_ic, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX, },
      { .dir = 2, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH, },
      { .dir = 2, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH, },
    },

    .num_diag_moments = 7,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP },
  };

  // Field.
  struct gkyl_gyrokinetic_field field = {
    .poisson_bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_NEUMANN, .value = {0.0} },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
      { .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_PERIODIC, },
      { .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_PERIODIC, },
    },
  };

  // GK app.
  struct gkyl_gk app_inp = {

    .cdim = ctx.cdim,
    .lower = { (2.5 * ctx.rho_s), -0.5 * ctx.Ltheta, -0.5 * ctx.Lz },
    .upper = { (2.5 * ctx.rho_s) + ctx.Lr, 0.5 * ctx.Ltheta, 0.5 * ctx.Lz },
    .cells = { cells_x[0], cells_x[1], cells_x[2] },
    .poly_order = 1,
    .basis_type = app_args.basis_type,

    .geometry = {
      .geometry_id = GKYL_GEOMETRY_MAPC2P,
      .mapc2p = mapc2p,
      .c2p_ctx = &ctx,
      .bfield_func = bfield_func,
      .bfield_ctx = &ctx,
      .position_map_info = {
        .maps[0] = mc2nu_r,
        .ctxs[0] = &ctx,
      },
    },

    .num_periodic_dir = 1,
    .periodic_dirs = { 1 },

    .num_species = 2,
    .species = { elc, ion },
    .field = field,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0], app_args.cuts[1], app_args.cuts[2] },
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
