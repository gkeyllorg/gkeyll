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

struct ar_react_ctx
{ 
  int cdim, vdim; // Dimensionality.

  // Physical constants (using non-normalized physical units).
  double epsilon0; // Permittivity of free space.
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.
  double mass_ion; // Proton mass.
  double charge_ion; // Proton charge.
  double mass_Ar1; // Ar1+ ion mass.
  double charge_Ar1; // Ar1+ ion charge.
  double mass_Ar2; // Ar2+ ion mass.
  double charge_Ar2; // Ar2+ ion charge.

  double Te; // Electron temperature.
  double Ti; // Ion temperature.
  double TAr1; // Ar1+ temperature.
  double TAr2; // Ar2+ temperature.
  double n0_elc; // Electron reference number density (1 / m^3).

  double B_axis; // Magnetic field axis (simple toroidal coordinates).
  double R0; // Major radius (simple toroidal coordinates).
  double a0; // Minor axis (simple toroidal coordinates).

  double nu_frac; // Collision frequency fraction.

  double k_perp_rho_s; // Product of perpendicular wavenumber and ion-sound gyroradius.

  // Derived physical quantities (using non-normalized physical units).
  double n0_Ar1; // Ar1+ reference number density (1 / m^3).
  double n0_Ar2; // Ar2+ reference number density (1 / m^3).
  double n0_ion; // Ion reference number density (1 / m^3.)

  double R; // Radial coordinate (simple toroidal coordinates).
  double B0; // Reference magnetic field strength (Tesla).

  double log_lambda_elc; // Electron Coulomb logarithm.
  double log_lambda_ion; // Ion Coulomb logarithm.
  double nu_elc; // Electron collision frequency.
  double nu_ion; // Ion collision frequency.

  double c_s; // Sound speed.
  double vte; // Electron thermal velocity.
  double vti; // Ion thermal velocity.
  double vtAr1; // Ar1+ thermal velocity.
  double vtAr2; // Ar2+ thermal velocity.
  double omega_ci; // Ion cyclotron frequency.
  double rho_s; // Ion-sound gyroradius.

  double k_perp; // Perpendicular wavenumber (for Poisson solver).

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
  double vpar_max_Ar1; // Domain boundary (Ar1+ velocity space: parallel velocity direction).
  double mu_max_Ar1; // Domain boundary (Ar1+ velocity space: magnetic moment direction).
  double vpar_max_Ar2; // Domain boundary (Ar2+ velocity space: parallel velocity direction).
  double mu_max_Ar2; // Domain boundary (Ar2+ velocity space: magnetic moment direction).
  int poly_order; // Polynomial order.
  double cfl_frac; // CFL coefficient.

  double t_end; // End time.
  int num_frames; // Number of output frames.
  double write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct ar_react_ctx
create_ctx(void)
{
  int cdim = 1, vdim = 2; // Dimensionality.

  // Physical constants (using non-normalized physical units).
  double epsilon0 = GKYL_EPSILON0; // Permittivity of free space.
  double mass_elc = GKYL_ELECTRON_MASS; // Electron mass.
  double mass_ion = 2.014 * GKYL_PROTON_MASS; // Proton mass.
  double mass_Ar1 = 39.95 * GKYL_PROTON_MASS; // Ar1+ ion mass.
  double mass_Ar2 = 39.95 * GKYL_PROTON_MASS; // Ar2+ ion mass.
  double charge_elc = -GKYL_ELEMENTARY_CHARGE; // Electron charge.
  double charge_ion = GKYL_ELEMENTARY_CHARGE; // Proton charge.
  double charge_Ar1 = GKYL_ELEMENTARY_CHARGE; // Ar1+ ion charge.
  double charge_Ar2 = 2.0 * GKYL_ELEMENTARY_CHARGE; // Ar2+ ion charge.

  double Te = 200.0 * GKYL_ELEMENTARY_CHARGE; // Electron temperature.
  double Ti = 200.0 * GKYL_ELEMENTARY_CHARGE; // Ion temperature.
  double TAr1 = 200.0 * GKYL_ELEMENTARY_CHARGE; // Ar1+ temperature.
  double TAr2 = 200.0 * GKYL_ELEMENTARY_CHARGE; // Ar2+ temperature.
  double n0_elc = 1.0e21; //  Electron reference number density (1 / m^3).

  double B_axis = 0.5; // Magnetic field axis (simple toroidal coordinates).
  double R0 = 0.85; // Major radius (simple toroidal coordinates).
  double a0 = 0.15; // Minor axis (simple toroidal coordinates).

  double nu_frac = 0.25; // Collision frequency fraction.

  double k_perp_rho_s = 0.3; // Product of perpendicular wavenumber and ion-sound gyroradius.

  // Derived physical quantities (using non-normalized physical units).
  double n0_Ar1 = n0_elc * 1.0e-1; // Ar1+ reference number density (1 / m^3).
  double n0_Ar2 = n0_elc * 1.0e-2; // Ar2+ reference number density (1 / m^3).
  double n0_ion = n0_elc - n0_Ar1 - (2.0 * n0_Ar2); // Ion reference number density (1 / m^3).

  double R = R0 + a0; // Radial coordinate (simple toroidal coordinates).
  double B0 = B_axis * (R0 / R); // Reference magnetic field strength (Tesla).

  double log_lambda_elc = 6.6 - 0.5 * log(n0_elc / 1.0e20) + 1.5 * log(Te / charge_ion); // Electron Coulomb logarithm.
  double log_lambda_ion = 6.6 - 0.5 * log(n0_elc / 1.0e20) + 1.5 * log(Ti / charge_ion); // Ion Coulomb logarithm.
  double nu_elc = nu_frac * log_lambda_elc * pow(charge_ion, 4.0) * n0_elc /
    (6.0 * sqrt(2.0) * pow(M_PI, 3.0 / 2.0) * pow(epsilon0, 2.0) * sqrt(mass_elc) * pow(Te, 3.0 / 2.0)); // Electron collision frequency.
  double nu_ion = nu_frac * log_lambda_ion * pow(charge_ion, 4.0) * n0_elc /
    (12.0 * pow(M_PI, 3.0 / 2.0) * pow(epsilon0, 2.0) * sqrt(mass_ion) * pow(Ti, 3.0 / 2.0)); // Ion collision frequency.
  
  double c_s = sqrt(Te / mass_ion); // Sound speed.
  double vte = sqrt(Te / mass_elc); // Electron thermal velocity.
  double vti = sqrt(Ti / mass_ion); // Ion thermal velocity.
  double vtAr1 = sqrt(TAr1 / mass_Ar1); // Ar1+ thermal velocity.
  double vtAr2 = sqrt(TAr2 / mass_Ar2); // Ar2+ thermal velocity.
  double omega_ci = fabs(charge_ion * B0 / mass_ion); // Ion cyclotron frequency.
  double rho_s = c_s / omega_ci; // Ion-sound gyroradius.

  double k_perp = k_perp_rho_s / rho_s; // Perpendicular wavenumber (for Poisson solver).

  // Simulation parameters.
  int Nz = 4; // Cell count (configuration space: z-direction).
  int Nvpar = 6; // Cell count (velocity space: parallel velocity direction).
  int Nmu = 4; // Cell count (velocity space: magnetic moment direction).
  double Lz = 4.0; // Domain size (configuration space: z-direction).
  double vpar_max_elc = 6.0 * vte; // Domain boundary (electron velocity space: parallel velocity direction).
  double mu_max_elc = (3.0 / 2.0) * 0.5 * mass_elc * pow(4.0 * vte, 2.0) / (2.0 * B0); // Domain boundary (electron velocity space: magnetic moment direction).
  double vpar_max_ion = 6.0 * vti; // Domain boundary (ion velocity space: parallel velocity direction).
  double mu_max_ion = (3.0 / 2.0) * 0.5 * mass_ion * pow(4.0 * vti, 2.0) / (2.0 * B0); // Domain boundary (ion velocity space: magnetic moment direction).
  double vpar_max_Ar1 = 6.0 * vtAr1; // Domain boundary (Ar1+ velocity space: parallel velocity direction).
  double mu_max_Ar1 = (3.0 / 2.0) * 0.5 * mass_Ar1 * pow(4.0 * vtAr1, 2.0) / (2.0 * B0); // Domain boundary (Ar1+ velocity space: magnetic moment direction).
  double vpar_max_Ar2 = 6.0 * vtAr2; // Domain boundary (Ar2+ velocity space: parallel velocity direction).
  double mu_max_Ar2 = (3.0 / 2.0) * 0.5 * mass_Ar2 * pow(4.0 * vtAr2, 2.0) / (2.0 * B0); // Domain boundary (Ar2+ velocity space: magnetic moment direction).
  int poly_order = 1; // Polynomial order.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 1.0e-7; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  double write_phase_freq = 0.2; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.
  
  struct ar_react_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .epsilon0 = epsilon0,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .mass_ion = mass_ion,
    .charge_ion = charge_ion,
    .mass_Ar1 = mass_Ar1,
    .charge_Ar1 = charge_Ar1,
    .mass_Ar2 = mass_Ar2,
    .charge_Ar2 = charge_Ar2,
    .Te = Te,
    .Ti = Ti,
    .TAr1 = TAr1,
    .TAr2 = TAr2,
    .n0_elc = n0_elc,
    .B_axis = B_axis,
    .R0 = R0,
    .a0 = a0,
    .nu_frac = nu_frac,
    .k_perp_rho_s = k_perp_rho_s,
    .n0_Ar1 = n0_Ar1,
    .n0_Ar2 = n0_Ar2,
    .n0_ion = n0_ion,
    .R = R,
    .B0 = B0,
    .log_lambda_elc = log_lambda_elc,
    .nu_elc = nu_elc,
    .log_lambda_ion = log_lambda_ion,
    .nu_ion = nu_ion,
    .c_s = c_s,
    .vte = vte,
    .vti = vti,
    .vtAr1 = vtAr1,
    .vtAr2 = vtAr2,
    .omega_ci = omega_ci,
    .rho_s = rho_s,
    .k_perp = k_perp,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nz, Nvpar, Nmu},
    .Lz = Lz,
    .vpar_max_elc = vpar_max_elc,
    .mu_max_elc = mu_max_elc,
    .vpar_max_ion = vpar_max_ion,
    .mu_max_ion = mu_max_ion,
    .vpar_max_Ar1 = vpar_max_Ar1,
    .mu_max_Ar1 = mu_max_Ar1,
    .vpar_max_Ar2 = vpar_max_Ar2,
    .mu_max_Ar2 = mu_max_Ar2,
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
  struct ar_react_ctx *app = ctx;

  double n0_elc = app->n0_elc;

  // Set electron total number density.
  fout[0] = n0_elc;
}

void
evalElcTempInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct ar_react_ctx *app = ctx;

  double Te = app->Te;

  // Set electron isotropic temperature.
  fout[0] = Te;
}

void
evalElcUparInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Set electron parallel velocity.
  fout[0] = 0.0;
}

void
evalIonDensityInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct ar_react_ctx *app = ctx;

  double n0_ion = app->n0_ion;

  // Set ion total number density.
  fout[0] = n0_ion;
}

void
evalIonTempInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct ar_react_ctx *app = ctx;

  double Ti = app->Ti;

  // Set ion isotropic temperature.
  fout[0] = Ti;
}

void
evalIonUparInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Set ion parallel velocity.
  fout[0] = 0.0;
}

void
evalAr1DensityInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct ar_react_ctx *app = ctx;

  double n0_Ar1 = app->n0_Ar1;

  // Set Ar1+ ion total number density.
  fout[0] = n0_Ar1;
}

void
evalAr1TempInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct ar_react_ctx *app = ctx;

  double TAr1 = app->TAr1;

  // Set Ar1+ ion isotropic temperature.
  fout[0] = TAr1;
}

void
evalAr1UparInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Set Ar1+ ion parallel velocity.
  fout[0] = 0.0;
}

void
evalAr2DensityInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct ar_react_ctx *app = ctx;

  double n0_Ar2 = app->n0_Ar2;

  // Set Ar2+ ion total number density.
  fout[0] = n0_Ar2;
}

void
evalAr2TempInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct ar_react_ctx *app = ctx;

  double TAr2 = app->TAr2;

  // Set Ar2+ ion isotropic temperature.
  fout[0] = TAr2;
}

void
evalAr2UparInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Set Ar2+ ion parallel velocity.
  fout[0] = 0.0;
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
  struct ar_react_ctx *app = ctx;

  double B0 = app->B0;

  // zc are computational coords. 
  // Set Cartesian components of magnetic field.
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = app->B0;
}

static inline void
mapc2p_vel_elc(double t, const double* GKYL_RESTRICT vc, double* GKYL_RESTRICT vp, void* ctx)
{
  struct ar_react_ctx *app = ctx;
  double cvpar = vc[0], cmu = vc[1];

  double vpar_max_elc = app->vpar_max_elc;
  double mu_max_elc = app->mu_max_elc;

  double vpar = 0.0;
  double mu = 0.0;

  if (cvpar < 0.0) {
    vpar = -vpar_max_elc * (cvpar * cvpar);
  }
  else {
    vpar = vpar_max_elc * (cvpar * cvpar);
  }
  mu = mu_max_elc * (cmu * cmu);

  // Set rescaled electron velocity space coordinates (vpar, mu) from old velocity space coordinates (cvpar, cmu):
  vp[0] = vpar; vp[1] = mu;
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

  struct ar_react_ctx ctx = create_ctx(); // Context for init functions.

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
    .lower = { -1.0, 0.0 },
    .upper = { 1.0, 1.0 },
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0_elc,

    .mapc2p = {
      .mapping = mapc2p_vel_elc,
      .ctx = &ctx,
    },

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .density = evalElcDensityInit,
      .ctx_density = &ctx,
      .temp = evalElcTempInit,
      .ctx_temp = &ctx,
      .upar = evalElcUparInit,
      .ctx_upar = &ctx,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .den_ref = ctx.n0_elc,
      .temp_ref = ctx.Te,
      .num_cross_collisions = 3,
      .collide_with = { "ion", "Ar1", "Ar2" },
    },

    .react = {
      .num_react = 2,
      .react_type = {
        {
          .react_id = GKYL_REACT_IZ, 
          .type_self = GKYL_SELF_ELC, 
          .ion_id = GKYL_ION_AR, 
          .elc_nm = "elc", 
          .ion_nm = "Ar2", 
          .donor_nm = "Ar1", 
          .charge_state = 1, 
          .ion_mass = ctx.mass_Ar2,
          .elc_mass = ctx.mass_elc, 
        }, 
        {
          .react_id = GKYL_REACT_RECOMB,
          .type_self = GKYL_SELF_ELC,
          .ion_id = GKYL_ION_AR,
          .elc_nm = "elc",
          .ion_nm = "Ar2",
          .recvr_nm = "Ar1",
          .charge_state = 1,
          .ion_mass = ctx.mass_Ar2,
          .elc_mass = ctx.mass_elc,
        },
      },
    },
    
    .num_diag_moments = 5,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP },
  };

  // Ions.
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.charge_ion, .mass = ctx.mass_ion,
    .vdim = ctx.vdim,
    .lower = { -ctx.vpar_max_ion, 0.0 },
    .upper = { ctx.vpar_max_ion, ctx.mu_max_ion },
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0_ion,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .density = evalIonDensityInit,
      .ctx_density = &ctx,
      .temp = evalIonTempInit,
      .ctx_temp = &ctx,
      .upar = evalIonUparInit,
      .ctx_upar = &ctx,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .den_ref = ctx.n0_elc,
      .temp_ref = ctx.Ti,
      .num_cross_collisions = 3,
      .collide_with = { "elc", "Ar1", "Ar2" },
    },
    
    .num_diag_moments = 5,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP },
  };

  // Ar1+ ions.
  struct gkyl_gyrokinetic_species Ar1 = {
    .name = "Ar1",
    .vdim = ctx.vdim,
    .charge = ctx.charge_Ar1, .mass = ctx.mass_Ar1,
    .lower = { -ctx.vpar_max_Ar1, 0.0 },
    .upper = { ctx.vpar_max_Ar1, ctx.mu_max_Ar1 },
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0_Ar1,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .density = evalAr1DensityInit,
      .ctx_density = &ctx,
      .temp = evalAr1TempInit,
      .ctx_temp = &ctx,
      .upar = evalAr1UparInit,
      .ctx_upar = &ctx,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .den_ref = ctx.n0_Ar1,
      .temp_ref = ctx.TAr1,
      .num_cross_collisions = 3,
      .collide_with = { "elc", "ion", "Ar2" },
    },

    .react = {
      .num_react = 2,
      .react_type = {
        {
          .react_id = GKYL_REACT_IZ, 
          .type_self = GKYL_SELF_DONOR, 
          .ion_id = GKYL_ION_AR, 
          .elc_nm = "elc", 
          .ion_nm = "Ar2", 
          .donor_nm = "Ar1",
          .charge_state = 1, 
          .ion_mass = ctx.mass_Ar2,
          .elc_mass = ctx.mass_elc, 
        }, 
        {
          .react_id = GKYL_REACT_RECOMB,
          .type_self = GKYL_SELF_RECVR,
          .ion_id = GKYL_ION_AR,
          .elc_nm = "elc",
          .ion_nm = "Ar2",
          .recvr_nm = "Ar1",
          .charge_state = 1,
          .ion_mass = ctx.mass_Ar2,
          .elc_mass = ctx.mass_elc,
        },
      },
    },
    
    .num_diag_moments = 5,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP },
  };

  // Ar2+ ions.
  struct gkyl_gyrokinetic_species Ar2 = {
    .name = "Ar2",
    .charge = ctx.charge_Ar2, .mass = ctx.mass_Ar2,
    .vdim = ctx.vdim,
    .lower = { -ctx.vpar_max_Ar2, 0.0 },
    .upper = { ctx.vpar_max_Ar2, ctx.mu_max_Ar2 },
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0_Ar2,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .density = evalAr2DensityInit,
      .ctx_density = &ctx,
      .temp = evalAr2TempInit,
      .ctx_temp = &ctx,
      .upar = evalAr2UparInit,
      .ctx_upar = &ctx,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .den_ref = ctx.n0_Ar2,
      .temp_ref = ctx.TAr2,
      .num_cross_collisions = 3,
      .collide_with = { "elc", "ion", "Ar1" },
    },

    .react = {
      .num_react = 2,
      .react_type = {
        {
          .react_id = GKYL_REACT_IZ, 
          .type_self = GKYL_SELF_ION,
          .ion_id = GKYL_ION_AR, 
          .elc_nm = "elc", 
          .ion_nm = "Ar2", 
          .donor_nm = "Ar1",
          .charge_state = 1, 
          .ion_mass = ctx.mass_Ar2,
          .elc_mass = ctx.mass_elc, 
        }, 
        {
          .react_id = GKYL_REACT_RECOMB,
          .type_self = GKYL_SELF_ION,
          .ion_id = GKYL_ION_AR,
          .elc_nm = "elc",
          .ion_nm = "Ar2",
          .recvr_nm = "Ar1",
          .charge_state = 1,
          .ion_mass = ctx.mass_Ar2,
          .elc_mass = ctx.mass_elc,
        },
      },
    },
    
    .num_diag_moments = 5,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP },
  };

  // Field.
  struct gkyl_gyrokinetic_field field = {
    .kperpSq = ctx.k_perp * ctx.k_perp,

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

    .num_species = 4,
    .species = { elc, ion, Ar1, Ar2 },

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
