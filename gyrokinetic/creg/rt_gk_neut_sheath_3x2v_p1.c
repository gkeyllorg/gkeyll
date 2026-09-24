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

struct sheath_ctx
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
  double T0; // Neutrals. 
  double n0; // Reference number density (1 / m^3).

  double B_axis; // Magnetic field axis (simple toroidal coordinates).
  double R0; // Major radius (simple toroidal coordinates).
  double a0; // Minor axis (simple toroidal coordinates).

  double nu_frac; // Collision frequency fraction.
  double rec_frac; 

  // Derived physical quantities (using non-normalized physical units).
  double R; // Radial coordinate (simple toroidal coordinates).
  double B0; // Reference magnetic field strength (Tesla).
  
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

  double n_src; // Source number density.
  double T_src; // Source temperature.
  double xmu_src; // Source mean position (x-direction).
  double xsigma_src; // Source standard deviation (x-direction).
  double floor_src; // Minimum source intensity.

  // Simulation parameters.
  int Nx; // Number of cells along x.
  int Ny; // Number of cells along y.
  int Nz; // Number of cells along z.
  int Nvpar; // Number of cells along vpar.
  int Nmu; // Number of cells along mu.
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  double Lx; // Domain length along x.
  double Ly; // Domain length along y.
  double Lz; // Domain length along z.
  double vpar_max_elc; // Maximum parallel electron velocity.
  double mu_max_elc;   // Maximum electron magnetic moment.
  double vpar_max_ion; // Maximum parallel ion velocity.
  double mu_max_ion;   // Maximum ion magnetic moment.
  double vmax_neut; // Neutrals.
  int poly_order; // Polynomial order.
  double cfl_frac; // CFL coefficient.

  double t_end; // End time.
  int num_frames; // Number of output frames.
  double write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct sheath_ctx
create_ctx(void)
{
  int cdim = 3, vdim = 2; // Dimensionality.

  // Physical constants (using non-normalized physical units).
  double epsilon0 = GKYL_EPSILON0; // Permittivity of free space.
  double mass_elc = GKYL_ELECTRON_MASS; // Electron mass.
  double mass_ion = 2.014 * GKYL_PROTON_MASS; // Proton mass.
  double charge_elc = -GKYL_ELEMENTARY_CHARGE; // Electron charge.
  double charge_ion = GKYL_ELEMENTARY_CHARGE; // Proton charge.

  double Te = 40.0 * GKYL_ELEMENTARY_CHARGE; // Electron temperature.
  double Ti = 40.0 * GKYL_ELEMENTARY_CHARGE; // Ion temperature.
  double T0 = 10.0 * GKYL_ELEMENTARY_CHARGE;
  double n0 = 7.0e18; //  Reference number density (1 / m^3).

  double B_axis = 0.5; // Magnetic field axis (simple toroidal coordinates).
  double R0 = 0.85; // Major radius (simple toroidal coordinates).
  double a0 = 0.15; // Minor axis (simple toroidal coordinates).

  double nu_frac = 0.1; // Collision frequency fraction.
  double rec_frac = 1.0;
  
  // Derived physical quantities (using non-normalized physical units).
  double R = R0 + a0; // Radial coordinate (simple toroidal coordinates).
  double B0 = B_axis * (R0 / R); // Reference magnetic field strength (Tesla).

  double log_lambda_elc = 6.6 - 0.5 * log(n0 / 1.0e20) + 1.5 * log(Te / charge_ion); // Electron Coulomb logarithm.
  double log_lambda_ion = 6.6 - 0.5 * log(n0 / 1.0e20) + 1.5 * log(Ti / charge_ion); // Ion Coulomb logarithm.
  double nu_elc = nu_frac * log_lambda_elc * pow(charge_ion, 4.0) * n0 /
    (6.0 * sqrt(2.0) * pow(M_PI, 3.0 / 2.0) * pow(epsilon0, 2.0) * sqrt(mass_elc) * pow(Te, 3.0 / 2.0)); // Electron collision frequency.
  double nu_ion = nu_frac * log_lambda_ion * pow(charge_ion, 4.0) * n0 /
    (12.0 * pow(M_PI, 3.0 / 2.0) * pow(epsilon0, 2.0) * sqrt(mass_ion) * pow(Ti, 3.0 / 2.0)); // Ion collision frequency.
  double nu_elc_ion = nu_elc*sqrt(2.0);
  double nu_ion_elc = nu_elc_ion*(mass_elc/mass_ion);
  
  double c_s = sqrt(Te / mass_ion); // Sound speed.
  double vte = sqrt(Te / mass_elc); // Electron thermal velocity.
  double vti = sqrt(Ti / mass_ion); // Ion thermal velocity.
  double vtn = sqrt(T0 / mass_ion);
  double omega_ci = fabs(charge_ion * B0 / mass_ion); // Ion cyclotron frequency.
  double rho_s = c_s / omega_ci; // Ion-sound gyroradius.

  double n_src = 1.4690539 * 3.612270e23; // Source number density.
  double T_src = 2.0 * Te; // Source temperature.
  double xmu_src = R; // Source mean position (x-direction).
  double xsigma_src = 0.005; // Source standard deviation (x-direction).
  double floor_src = 0.1; // Minimum source intensity.

  // Grid parameters.
  int Nx = 4; // Number of cells along x.
  int Ny = 2; // Number of cells along y.
  int Nz = 8; // Number of cells along z.
  int Nvpar = 6; // Number of cells along vpar.
  int Nmu = 4; // Number of cells along mu.
  double Lx = 50.0*rho_s; // Domain length along x.
  double Ly = 100.0*rho_s; // Domain length along 1.
  double Lz = 4.0; // Domain length along z.
  double vpar_max_elc = 4.0*vte; // Maximum parallel electron velocity.
  double mu_max_elc = (3.0/2.0)*0.5*mass_elc*pow(4.0*vte,2.0)/(2.0*B0); // Maximum electron magnetic moment.
  double vpar_max_ion = 4.0*vti; // Maximum parallel ion velocity.
  double mu_max_ion = (3.0/2.0)*0.5*mass_ion*pow(4.0*vti,2.0)/(2.0*B0); // Maximum ion magnetic moment.
  double vmax_neut = 4.0*vtn; 
  int poly_order = 1; // Polynomial order.
  double cfl_frac = 0.50; // CFL coefficient.

  double t_end = 2.0e-7; // Final simulation time.
  int num_frames = 2; // Number of output frames.
  double write_phase_freq = 0.2; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct sheath_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .epsilon0 = epsilon0,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .mass_ion = mass_ion,
    .charge_ion = charge_ion,
    .Te = Te,
    .Ti = Ti,
    .T0 = T0,
    .n0 = n0,
    .B_axis = B_axis,
    .R0 = R0,
    .a0 = a0,
    .nu_frac = nu_frac,
    .rec_frac = rec_frac,
    .R = R,
    .B0 = B0,
    .log_lambda_elc = log_lambda_elc,
    .log_lambda_ion = log_lambda_ion,
    .nu_elc = nu_elc,
    .nu_ion = nu_ion,
    .nu_elc_ion = nu_elc_ion,
    .nu_ion_elc = nu_ion_elc,
    .c_s = c_s,
    .vte = vte,
    .vti = vti,
    .omega_ci = omega_ci,
    .rho_s = rho_s,
    .n_src = n_src,
    .T_src = T_src,
    .xmu_src = xmu_src,
    .xsigma_src = xsigma_src,
    .floor_src = floor_src,
    .Nx = Nx,
    .Ny = Ny,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nx,Ny,Nz,Nvpar,Nmu},
    .Lx = Lx,
    .Ly = Ly,
    .Lz = Lz,
    .vpar_max_elc = vpar_max_elc,
    .mu_max_elc = mu_max_elc,
    .vpar_max_ion = vpar_max_ion,
    .mu_max_ion = mu_max_ion,
    .vmax_neut = vmax_neut,
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
  struct sheath_ctx *app = ctx;
  double x = xn[0], z = xn[2];

  double mass_ion = app->mass_ion;

  double n_src = app->n_src;
  double T_src = app->T_src;
  double xmu_src = app->xmu_src;
  double xsigma_src = app->xsigma_src;
  double floor_src = app->floor_src;

  double Lz = app->Lz;

  double src_density = GKYL_MAX2(exp(-((x - xmu_src) * (x - xmu_src)) / ((2.0 * xsigma_src) * (2.0 * xsigma_src))), floor_src) * n_src;
  double src_temp = 0.0;
  double n = 0;

  if (x < xmu_src + 3.0 * xsigma_src) {
    src_temp = T_src;
  }
  else {
    src_temp = (3.0 / 8.0) * T_src;
  }

  double c_s_src = sqrt((5.0 / 3.0) * src_temp / mass_ion);
  double n_peak = 4.0 * sqrt(5.0) / 3.0 / c_s_src * (0.125 * Lz) * src_density;

  if (fabs(z) <= 0.25 * Lz) {
    n = 0.5 * n_peak * (1.0 + sqrt(1.0 - (z / (0.25 * Lz)) * (z / (0.25 * Lz)))); // Electron total number density (left).
  }
  else {
    n = 0.5 * n_peak; // Electron total number density (right).
  }

  // Set electron total number density.
  fout[0] = n;
}

void
evalElcTempInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sheath_ctx *app = ctx;
  double x = xn[0];

  double Te = app->Te;

  double xmu_src = app->xmu_src;
  double xsigma_src = app->xsigma_src;

  double T = 0.0;

  if (x < xmu_src + 3.0 * xsigma_src) {
    T = (5.0 / 4.0) * Te; // Electron isotropic temperature (left).
  }
  else {
    T = 0.5 * Te; // Electron isotropic temperature (right).
  }

  // Set electron isotropic temperature.
  fout[0] = T;
}

void
evalElcUparInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Set electron parallel velocity.
  fout[0] = 0.0;
}

void
evalElcSourceDensityInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sheath_ctx *app = ctx;
  double x = xn[0], z = xn[2];

  double n_src = app->n_src;
  double xmu_src = app->xmu_src;
  double xsigma_src = app->xsigma_src;
  double floor_src = app->floor_src;

  double Lz = app->Lz;

  double n = 0.0;

  if (fabs(z) < 0.25 * Lz) {
    n = GKYL_MAX2(exp(-((x - xmu_src) * (x - xmu_src)) / ((2.0 * xsigma_src) * (2.0 * xsigma_src))),
      floor_src) * n_src; // Electron source total number density (left).
  }
  else {
    n = 1.0e-40 * n_src; // Electron source total number density (right).
  }

  // Set electron source total number density.
  fout[0] = n;
}

void
evalElcSourceTempInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sheath_ctx *app = ctx;
  double x = xn[0];

  double T_src = app->T_src;
  double xmu_src = app->xmu_src;
  double xsigma_src = app->xsigma_src;

  double T = 0.0;

  if (x < xmu_src + 3.0 * xsigma_src) {
    T = T_src; // Electron source isotropic temperature (left).
  }
  else {
    T = (3.0 / 8.0) * T_src; // Electron source isotropic temperature (right).
  }

  // Set electron source isotropic temperature.
  fout[0] = T;
}

void
evalElcSourceUparInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Set electron source parallel velocity.
  fout[0] = 0.0;
}

void
evalIonDensityInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sheath_ctx *app = ctx;
  double x = xn[0], z = xn[2];

  double mass_ion = app->mass_ion;

  double n_src = app->n_src;
  double T_src = app->T_src;
  double xmu_src = app->xmu_src;
  double xsigma_src = app->xsigma_src;
  double floor_src = app->floor_src;

  double Lz = app->Lz;

  double src_density = GKYL_MAX2(exp(-((x - xmu_src) * (x - xmu_src)) / ((2.0 * xsigma_src) * (2.0 * xsigma_src))), floor_src) * n_src;
  double src_temp = 0.0;
  double n = 0;

  if (x < xmu_src + 3.0 * xsigma_src) {
    src_temp = T_src;
  }
  else {
    src_temp = (3.0 / 8.0) * T_src;
  }

  double c_s_src = sqrt((5.0 / 3.0) * src_temp / mass_ion);
  double n_peak = 4.0 * sqrt(5.0) / 3.0 / c_s_src * (0.125 * Lz) * src_density;

  if (fabs(z) <= 0.25 * Lz) {
    n = 0.5 * n_peak * (1.0 + sqrt(1.0 - (z / (0.25 * Lz)) * (z / (0.25 * Lz)))); // Ion total number density (left).
  }
  else {
    n = 0.5 * n_peak; // Ion total number density (right).
  }

  // Set ion total number density.
  fout[0] = n;
}

void
evalIonTempInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sheath_ctx *app = ctx;
  double x = xn[0];

  double Ti = app->Ti;

  double xmu_src = app->xmu_src;
  double xsigma_src = app->xsigma_src;

  double T = 0.0;

  if (x < xmu_src + 3.0 * xsigma_src) {
    T = (5.0 / 4.0) * Ti; // Ion isotropic temperature (left).
  }
  else {
    T = 0.5 * Ti; // Ion isotropic temperature (right).
  }

  // Set ion isotropic temperature.
  fout[0] = T;
}

void
evalIonUparInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Set ion parallel velocity.
  fout[0] = 0.0;
}

// Neutral profiles.
void neut_density_init(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sheath_ctx *app = ctx;
  double z = xn[2];
  double n0 = app->n0;
  double Lz = app->Lz;
  double n = 0.0;

  //Set number density.
  if (z <= 0) {
    n = n0*(pow(1.0/cosh(-(Lz/2. + z)/0.2),2.0) + 1.e-6);
  }
  else {
    n = n0*(pow(1.0/cosh((-Lz/2. + z)/0.2),2.0) + 1.e-6);
  }
  fout[0] = n; 
}

void unit_density(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sheath_ctx *app = ctx;
  fout[0] = 1.0; 
}

void temp_neut(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct sheath_ctx *app = ctx;
  double T = app->T0;
  fout[0] = T;
}

void udrift(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0; 
  fout[1] = 0.0;
  fout[2] = 0.0;
}


void
evalIonSourceDensityInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sheath_ctx *app = ctx;
  double x = xn[0], z = xn[2];

  double n_src = app->n_src;
  double xmu_src = app->xmu_src;
  double xsigma_src = app->xsigma_src;
  double floor_src = app->floor_src;

  double Lz = app->Lz;

  double n = 0.0;

  if (fabs(z) < 0.25 * Lz) {
    n = GKYL_MAX2(exp(-((x - xmu_src) * (x - xmu_src)) / ((2.0 * xsigma_src) * (2.0 * xsigma_src))),
      floor_src) * n_src; // Ion source total number density (left).
  }
  else {
    n = 1.0e-40 * n_src; // Ion source total number density (right).
  }

  // Set ion source total number density.
  fout[0] = n;
}

void
evalIonSourceTempInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sheath_ctx *app = ctx;
  double x = xn[0];

  double T_src = app->T_src;
  double xmu_src = app->xmu_src;
  double xsigma_src = app->xsigma_src;

  double T = 0.0;

  if (x < xmu_src + 3.0 * xsigma_src) {
    T = T_src; // Ion source isotropic temperature (left).
  }
  else {
    T = (3.0 / 8.0) * T_src; // Ion source isotropic temperature (right).
  }

  // Set ion source isotropic temperature.
  fout[0] = T;
}

void
evalIonSourceUparInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Set ion source parallel velocity.
  fout[0] = 0.0;
}

void
evalNuElc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sheath_ctx *app = ctx;

  double nu_elc = app->nu_elc;

  // Set electron collision frequency.
  fout[0] = nu_elc;
}

void
evalNuIon(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sheath_ctx *app = ctx;

  double nu_ion = app->nu_ion;

  // Set ion collision frequency.
  fout[0] = nu_ion;
}

void
evalNuElcIon(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sheath_ctx *app = ctx;

  double nu_elc_ion = app->nu_elc_ion;

  // Set electron-ion collision frequency.
  fout[0] = nu_elc_ion;
}

void
evalNuIonElc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sheath_ctx *app = ctx;

  double nu_ion_elc = app->nu_ion_elc;

  // Set ion-electron collision frequency.
  fout[0] = nu_ion_elc;
}

static inline void
mapc2p(double t, const double* GKYL_RESTRICT zc, double* GKYL_RESTRICT xp, void* ctx)
{
  struct sheath_ctx *app = ctx;
  double x = zc[0], y = zc[1], z = zc[2];

  double R0 = app->R0;
  double a0 = app->a0;

  double R = x;
  double phi = z / (R0 + a0);
  double X = R * cos(phi);
  double Y = R * sin(phi);
  double Z = y;

  // Set physical coordinates (X, Y, Z) from computational coordinates (x, y, z).
//  xp[0] = X; xp[1] = Y; xp[2] = Z;
  xp[0] = x; xp[1] = y; xp[2] = z;
}

void
bfield_func(double t, const double* GKYL_RESTRICT zc, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sheath_ctx *app = ctx;
  double x = zc[0];

  double B0 = app->B0;
  double R = app->R;

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

  struct sheath_ctx ctx = create_ctx(); // Context for init functions.

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  struct gkyl_gyrokinetic_emission_inp neut_bc = {
    .num_species = 1,
    .in_species = { "ion" },
    .recycling_frac = ctx.rec_frac,
    .emission_temp = ctx.T0,
  };

  struct gkyl_gyrokinetic_projection recyc_proj = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
    .ctx_density = &ctx,
    .density = unit_density,
    .ctx_upar = &ctx,
    .udrift= udrift,
    .ctx_temp = &ctx,
    .temp = temp_neut,
  };
  
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
        .density = evalElcSourceDensityInit,
        .ctx_density = &ctx,
        .temp = evalElcSourceTempInit,
        .ctx_temp = &ctx,
        .upar = evalElcSourceUparInit,
        .ctx_upar = &ctx,
      },
      .diagnostics = {
        .num_diag_moments = 1,
        .diag_moments = { GKYL_F_MOMENT_M0M1M2PARM2PERP },
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
      }
    },

    .react_neut = {
      .num_react = 2,
      .react_type = {
        { .react_id = GKYL_REACT_IZ,
          .type_self = GKYL_SELF_ELC,
          .ion_id = GKYL_ION_H,
	  .elc_nm = "elc",
          .ion_nm = "ion",
          .donor_nm = "D0",
	  .charge_state = 0,
          .ion_mass = ctx.mass_ion,
          .elc_mass = ctx.mass_elc,
        },
	{ .react_id = GKYL_REACT_RECOMB,
          .type_self = GKYL_SELF_ELC,
          .ion_id = GKYL_ION_H,
	  .elc_nm = "elc",
          .ion_nm = "ion",
          .recvr_nm = "D0",
	  .charge_state = 0,
          .ion_mass = ctx.mass_ion,
          .elc_mass = ctx.mass_elc,
        },
      },
    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX, },
      { .dir = 2, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH, },
      { .dir = 2, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH, },
    },

    .num_diag_moments = 2,
    .diag_moments = { GKYL_F_MOMENT_M0M1M2PARM2PERP, GKYL_F_MOMENT_BIMAXWELLIAN },
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    .time_rate_diagnostics = true,

    .boundary_flux_diagnostics = {
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    },
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
      .self_nu = evalNuIon,
      .self_nu_ctx = &ctx,
      .num_cross_collisions = 1,
      .collide_with = { "elc" },
      .cross_nu = { evalNuIonElc, },
      .cross_nu_ctx = &ctx,
      .den_ref = ctx.n0,
      .temp_ref = ctx.Ti,
    },

    .react_neut = {
      .num_react = 3,
      .react_type = {
        { .react_id = GKYL_REACT_IZ,
          .type_self = GKYL_SELF_ION,
          .ion_id = GKYL_ION_H,
	  .elc_nm = "elc",
          .ion_nm = "ion",
          .donor_nm = "D0",
	  .charge_state = 0,
          .ion_mass = ctx.mass_ion,
          .elc_mass = ctx.mass_elc,
        },
	{ .react_id = GKYL_REACT_RECOMB,
          .type_self = GKYL_SELF_ION,
          .ion_id = GKYL_ION_H,
	  .elc_nm = "elc",
          .ion_nm = "ion",
          .recvr_nm = "D0",
	  .charge_state = 0,
          .ion_mass = ctx.mass_ion,
          .elc_mass = ctx.mass_elc,
        },
	{ .react_id = GKYL_REACT_CX,
          .type_self = GKYL_SELF_ION,
          .ion_id = GKYL_ION_H,
          .ion_nm = "ion",
          .partner_nm = "D0",
          .ion_mass = ctx.mass_ion,
          .partner_mass = ctx.mass_ion,
        },
      },
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,

      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
        .density = evalIonSourceDensityInit,
        .ctx_density = &ctx,
        .temp = evalIonSourceTempInit,
        .ctx_temp = &ctx,
        .upar = evalIonSourceUparInit,
        .ctx_upar = &ctx,
      },
      .diagnostics = {
        .num_diag_moments = 1,
        .diag_moments = { GKYL_F_MOMENT_M0M1M2PARM2PERP },
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
//        .time_integrated = true,
      }
    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX, },
      { .dir = 2, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH, },
      { .dir = 2, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH, },
    },
    
    .num_diag_moments = 2,
    .diag_moments = { GKYL_F_MOMENT_M0M1M2PARM2PERP, GKYL_F_MOMENT_BIMAXWELLIAN },
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    .time_rate_diagnostics = true,

    .boundary_flux_diagnostics = {
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    },
  };

  struct gkyl_gyrokinetic_neut_species D0 = {
    .name = "D0",
    .mass = ctx.mass_ion,
    .vdim = ctx.vdim+1,
    .lower = { -ctx.vmax_neut, -ctx.vmax_neut, -ctx.vmax_neut},
    .upper = { ctx.vmax_neut, ctx.vmax_neut, ctx.vmax_neut },
    .cells = { cells_v[1], cells_v[1], cells_v[1] },
    
    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = neut_density_init,
      .ctx_upar = &ctx,
      .udrift= udrift,
      .ctx_temp = &ctx,
      .temp = temp_neut,      
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_NEUTRAL,
    },

    .react_neut = {
      .num_react = 3,
      .react_type = {
        { .react_id = GKYL_REACT_IZ,
          .type_self = GKYL_SELF_DONOR,
          .ion_id = GKYL_ION_H,
    	  .elc_nm = "elc",
          .ion_nm = "ion",
          .donor_nm = "D0",
    	  .charge_state = 0,
          .ion_mass = ctx.mass_ion,
          .elc_mass = ctx.mass_elc,
        },
    	{ .react_id = GKYL_REACT_RECOMB,
          .type_self = GKYL_SELF_RECVR,
          .ion_id = GKYL_ION_H,
    	  .elc_nm = "elc",
          .ion_nm = "ion",
          .donor_nm = "D0",
    	  .charge_state = 0,
          .ion_mass = ctx.mass_ion,
          .elc_mass = ctx.mass_elc,
        },
    	{ .react_id = GKYL_REACT_CX,
          .type_self = GKYL_SELF_PARTNER,
          .ion_id = GKYL_ION_H,
          .ion_nm = "ion",
          .partner_nm = "D0",
          .ion_mass = ctx.mass_ion,
          .partner_mass = ctx.mass_ion,
        },
      },
    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB, },
      { .dir = 2, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_RECYCLE, .emission = neut_bc, .write_diagnostics = true, },
      { .dir = 2, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_RECYCLE, .emission = neut_bc, .write_diagnostics = true, },
    },
   
    .num_diag_moments = 3,
    .diag_moments = { GKYL_F_MOMENT_M1_FROM_H, GKYL_F_MOMENT_ENERGY, GKYL_F_MOMENT_LTE},
  };

  // Field.
  struct gkyl_gyrokinetic_field field = {
    .poisson_bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
      { .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_PERIODIC },
      { .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_PERIODIC },
    },
    .time_rate_diagnostics = true,
  };

  // Gyrokinetic app.
  struct gkyl_gk app_inp = {

    .cdim = ctx.cdim,
    .lower = { ctx.R - (0.5 * ctx.Lx), -0.5 * ctx.Ly, -0.5 * ctx.Lz},
    .upper = { ctx.R + (0.5 * ctx.Lx),  0.5 * ctx.Ly,  0.5 * ctx.Lz},
    .cells = { cells_x[0], cells_x[1], cells_x[2] },

    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,
    .cfl_frac = ctx.cfl_frac,

    .geometry = {
      .geometry_id = GKYL_GEOMETRY_MAPC2P,
      .mapc2p = mapc2p,
      .c2p_ctx = &ctx,
      .bfield_func = bfield_func,
      .bfield_ctx = &ctx
    },

    .num_periodic_dir = 1,
    .periodic_dirs = { 1 },

    .num_species = 2,
    .species = { elc, ion },
    .num_neut_species = 1,
    .neut_species = { D0 },

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
