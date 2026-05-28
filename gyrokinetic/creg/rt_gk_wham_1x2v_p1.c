#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_eqn_type.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_gyrokinetic_run.h>
#include <gkyl_math.h>

#include <rt_arg_parse.h>

// Define the context of the simulation. This is basically all the globals
struct gk_mirror_ctx
{
  int cdim, vdim; // Dimensionality.

  // Plasma parameters
  double mi;
  double qi;
  double me;
  double qe;
  double Te0;
  double n0;
  double B_p;
  double beta;
  double tau;
  double Ti0;
  double kperpRhos;
  // Electron-electron collision freq.
  double logLambdaElc;
  double nuElc;
  // Ion-ion collision freq.
  double logLambdaIon;
  double nuIon;
  // Thermal speeds.
  double vti;
  double vte;
  double c_s;
  // Gyrofrequencies and gyroradii.
  double omega_ci;
  double rho_s;
  double kperp; // Perpendicular wavenumber in SI units.
  // Axial coordinate Z extents. Endure that Z=0 is not on
  double z_min;
  double z_max;
  double psi_eval;
  double z_m;
  // Physics parameters at mirror throat
  double n_m;
  double Te_m;
  double Ti_m;
  double Ti_perp0;
  double Ti_par0;
  double Ti_perp_m;
  double Ti_par_m;
  double Te_perp0;
  double Te_par0;
  double Te_perp_m;
  double Te_par_m;
  double cs_m;
  // Source parameters
  double NSrcIon;
  double lineLengthSrcIon;
  double sigSrcIon;
  double NSrcFloorIon;
  double TSrc0Ion;
  double TSrcFloorIon;
  double NSrcElc;
  double lineLengthSrcElc;
  double sigSrcElc;
  double NSrcFloorElc;
  double TSrc0Elc;
  double TSrcFloorElc;
  // Grid parameters
  double vpar_max_ion;
  double vpar_max_elc;
  double mu_max_ion;
  double mu_max_elc;
  int Nz;
  int Nvpar;
  int Nmu;
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  int poly_order;

  double t_end; // End time.
  int num_frames; // Number of output frames.
  double write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};


// -- Source functions.
void
eval_density_elc_source(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double z = xn[0];
  double NSrc = app->NSrcElc;
  double zSrc = app->lineLengthSrcElc;
  double sigSrc = app->sigSrcElc;
  double NSrcFloor = app->NSrcFloorElc;
  if (fabs(z) <= app->z_m)
  {
    fout[0] = fmax(NSrcFloor, (NSrc / sqrt(2.0 * M_PI * pow(sigSrc, 2.))) *
                                  exp(-1 * pow((z - zSrc), 2) / (2.0 * pow(sigSrc, 2.))));
  }
  else
  {
    fout[0] = 1e-16;
  }
}

void
eval_upar_elc_source(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
eval_temp_elc_source(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double z = xn[0];
  double sigSrc = app->sigSrcElc;
  double TSrc0 = app->TSrc0Elc;
  double Tfloor = app->TSrcFloorElc;
  if (fabs(z) <= 2.0 * sigSrc)
  {
    fout[0] = TSrc0;
  }
  else
  {
    fout[0] = Tfloor;
  }
}

void
eval_density_ion_source(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double z = xn[0];
  double NSrc = app->NSrcIon;
  double zSrc = app->lineLengthSrcIon;
  double sigSrc = app->sigSrcIon;
  double NSrcFloor = app->NSrcFloorIon;
  if (fabs(z) <= app->z_m)
  {
    fout[0] = fmax(NSrcFloor, (NSrc / sqrt(2.0 * M_PI * pow(sigSrc, 2))) *
                                  exp(-1 * pow((z - zSrc), 2) / (2.0 * pow(sigSrc, 2))));
  }
  else
  {
    fout[0] = 1e-16;
  }
}

void
eval_upar_ion_source(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
eval_temp_ion_source(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double z = xn[0];
  double sigSrc = app->sigSrcIon;
  double TSrc0 = app->TSrc0Ion;
  double Tfloor = app->TSrcFloorIon;
  if (fabs(z) <= 2.0 * sigSrc)
  {
    fout[0] = TSrc0;
  }
  else
  {
    fout[0] = Tfloor;
  }
}

// Electrons initial conditions
void
eval_density_elc(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double z = xn[0];
  double z_m = app->z_m;
  double sigma = 0.9*z_m;
  if (fabs(z) <= sigma)
  {
    fout[0] = 0.5*app->n0*(1. + tanh(10. * sigma * fabs(sigma - fabs(z))));
  }
  else
  {
    fout[0] = 0.5*app->n0*exp(-5 * (fabs(sigma - fabs(z))));
  }
}

void
eval_upar_elc(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double z = xn[0];
  double cs_m = app->cs_m;
  double z_m = app->z_m;
  double z_max = app->z_max;
  if (fabs(z) <= z_m)
  {
    fout[0] = 0.0;
  }
  else
  {
    fout[0] = fabs(z) / z * cs_m * tanh(3 * (z_max - z_m) * fabs(fabs(z) - z_m)); // Maybe put a 5 here
  }
}

void
eval_temp_par_elc(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double z = xn[0];
  double z_m = app->z_m;
  double Te_par0 = app->Te_par0;
  double Te_par_m = app->Te_par_m;
  if (fabs(z) <= z_m)
  {
    fout[0] = Te_par_m+(Te_par0-Te_par_m)*tanh(4 * fabs(z_m - fabs(z)));
  }
  else
  {
    fout[0] = Te_par_m;
  }
}

void
eval_temp_perp_elc(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double z = xn[0];
  double z_m = app->z_m;
  double Te_perp0 = app->Te_perp0;
  double Te_perp_m = app->Te_perp_m;
  if (fabs(z) <= z_m)
  {
   fout[0] = Te_perp_m - Te_perp0*tanh(3.*fabs(z_m-fabs(z)));
  }
  else
  {
    fout[0] = Te_perp_m * GKYL_MAX2(1.e-3, exp(-5. * (fabs(z_m - fabs(z)))));
  }
}

// Ion initial conditions
void
eval_density_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double z = xn[0];
  double z_m = app->z_m;
  double sigma = 0.9*z_m;
  if (fabs(z) <= sigma)
  {
    fout[0] = 0.5*app->n0*(1. + tanh(10. * sigma * fabs(sigma - fabs(z))));
  }
  else
  {
    fout[0] = 0.5*app->n0* exp(-5 * (fabs(sigma - fabs(z))));
  }
}

void
eval_upar_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double z = xn[0];
  double cs_m = app->cs_m;
  double z_m = app->z_m;
  double z_max = app->z_max;
  if (fabs(z) <= z_m)
  {
    fout[0] = 0.0;
  }
  else
  {
    fout[0] = fabs(z) / z * cs_m * tanh(3 * (z_max - z_m) * fabs(fabs(z) - z_m)); // Maybe put a 5 here
  }
}
void
eval_temp_par_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double z = xn[0];
  double z_m = app->z_m;
  double Ti_par0 = app->Ti_par0;
  double Ti_par_m = app->Ti_par_m;
  if (fabs(z) <= z_m)
  {
    fout[0] = Ti_par_m + (Ti_par0 - Ti_par_m) * tanh(4 * fabs(z_m - fabs(z)));
  }
  else
  {
    fout[0] = Ti_par_m * GKYL_MAX2(1.e-2, 4 * log(fabs(fabs(z) - z_m) + 1));
  }
}

void
eval_temp_perp_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double z = xn[0];
  double z_m = app->z_m;
  double Ti_perp0 = app->Ti_perp0;
  double Ti_perp_m = app->Ti_perp_m;
  if (fabs(z) <= z_m)
  {
    fout[0] = Ti_perp_m + (Ti_perp0 - Ti_perp_m) * tanh(3. * fabs(z_m - fabs(z)));
  }
  else
  {
    fout[0] = Ti_perp_m * GKYL_MAX2(1.e-3, exp(-5. * (fabs(z_m - fabs(z)))));
  }
}

void mapc2p_vel_ion(double t, const double *vc, double* GKYL_RESTRICT vp, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double vpar_max_ion = app->vpar_max_ion;
  double mu_max_ion = app->mu_max_ion;

  double cvpar = vc[0], cmu = vc[1];
  double b = 1.45;
  double linear_velocity_threshold = 1./3.;
  double frac_linear = 1/b*atan(linear_velocity_threshold*tan(b));
  if (fabs(cvpar) < frac_linear) {
    double func_frac = tan(frac_linear*b) / tan(b);
    vp[0] = vpar_max_ion*func_frac*cvpar/frac_linear;
  }
  else {
    vp[0] = vpar_max_ion*tan(cvpar*b)/tan(b);
  }
  // Quadratic map in mu.
  vp[1] = mu_max_ion*pow(cmu,2);
}

void mapc2p_vel_elc(double t, const double *vc, double* GKYL_RESTRICT vp, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double vpar_max_elc = app->vpar_max_elc;
  double mu_max_elc = app->mu_max_elc;

  double cvpar = vc[0], cmu = vc[1];
  double b = 1.45;
  double linear_velocity_threshold = 1./6.;
  double frac_linear = 1/b*atan(linear_velocity_threshold*tan(b));
  if (fabs(cvpar) < frac_linear) {
    double func_frac = tan(frac_linear*b) / tan(b);
    vp[0] = vpar_max_elc*func_frac*cvpar/frac_linear;
  }
  else {
    vp[0] = vpar_max_elc*tan(cvpar*b)/tan(b);
  }
  // Quadratic map in mu.
  vp[1] = mu_max_elc*pow(cmu,2);
}

struct gk_mirror_ctx
create_ctx(void)
{
  int cdim = 1, vdim = 2; // Dimensionality.

  // Universal constant parameters.
  double eps0 = GKYL_EPSILON0;
  double mu0 = GKYL_MU0; // Not sure if this is right
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mp = GKYL_PROTON_MASS; // ion mass
  double me = GKYL_ELECTRON_MASS;
  double qi = eV;  // ion charge
  double qe = -eV; // electron charge

  // Plasma parameters.
  double mi = 2.014 * mp;
  double Te0 = 940 * eV;
  double n0 = 3e19;
  double B_p = 0.53;
  double beta = 0.4;
  double tau = pow(B_p, 2.) * beta / (2.0 * mu0 * n0 * Te0) - 1.;
  double Ti0 = tau * Te0;
  double kperpRhos = 0.1;

  double nuFrac = 1.0;
  // Electron-electron collision freq.
  double logLambdaElc = 6.6 - 0.5 * log(n0 / 1e20) + 1.5 * log(Te0 / eV);
  double nuElc = nuFrac * logLambdaElc * pow(eV, 4.) * n0 /
                 (6. * sqrt(2.) * pow(M_PI, 3. / 2.) * pow(eps0, 2.) * sqrt(me) * pow(Te0, 3. / 2.));
  // Ion-ion collision freq.
  double logLambdaIon = 6.6 - 0.5 * log(n0 / 1e20) + 1.5 * log(Ti0 / eV);
  double nuIon = nuFrac * logLambdaIon * pow(eV, 4.) * n0 /
                 (12 * pow(M_PI, 3. / 2.) * pow(eps0, 2.) * sqrt(mi) * pow(Ti0, 3. / 2.));

  // Thermal speeds.
  double vti = sqrt(Ti0 / mi);
  double vte = sqrt(Te0 / me);
  double c_s = sqrt(Te0 / mi);

  // Gyrofrequencies and gyroradii.
  double omega_ci = eV * B_p / mi;
  double rho_s = c_s / omega_ci;

  // Perpendicular wavenumber in SI units:
  double kperp = kperpRhos / rho_s;

  // Geometry parameters.
  // Axial coordinate Z extents. Endure that Z=0 is not on
  // the boundary of a cell (due to AD errors).
  double z_min = -2.0;
  double z_max =  2.0;
  double psi_eval = 2.5e-3;

  double z_m = 1.0;

  // Source parameters
  double NSrcIon = 3.1715e23 / 8.0;
  double lineLengthSrcIon = 0.0;
  double sigSrcIon = z_m / 4.0;
  double NSrcFloorIon = 0.05 * NSrcIon;
  double TSrc0Ion = Ti0 * 1.25;
  double TSrcFloorIon = TSrc0Ion / 8.0;
  double NSrcElc = NSrcIon;
  double lineLengthSrcElc = lineLengthSrcIon;
  double sigSrcElc = sigSrcIon;
  double NSrcFloorElc = NSrcFloorIon;
  double TSrc0Elc = TSrc0Ion / tau;
  double TSrcFloorElc = TSrcFloorIon / tau;

  // Initial conditions parameters
  double Ti_perp0 = 10000 * eV;
  double Ti_perp_m = 15000 * eV;
  double Ti_par0 = 7500 * eV;
  double Ti_par_m = 1000 * eV;

  double Te_par0 = 1800 * eV;
  double Te_par_m = 300 * eV;
  double Te_perp0 = 2000 * eV;
  double Te_perp_m = 3000 * eV;

  // Physics parameters at mirror throat
  double n_m = 1.105617e19;
  double Te_m = 346.426583 * eV;
  double Ti_m = 3081.437703 * eV;
  double cs_m = 4.037740e5;

  // Grid parameters
  double vpar_max_elc = 20 * vte;
  double mu_max_elc = me * pow(3. * vte, 2.) / (2. * B_p);
  double vpar_max_ion = 20 * vti;
  double mu_max_ion = mi * pow(3. * vti, 2.) / (2. * B_p);
  int Nz = 32;
  int Nvpar = 8; // Number of cells in the paralell velocity direction 96
  int Nmu = 8;  // Number of cells in the mu direction 192
  int poly_order = 1;

  double t_end = 1e-8;
  int num_frames = 1;
  double write_phase_freq = 0.2; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gk_mirror_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .mi = mi,
    .qi = qi,
    .me = me,
    .qe = qe,
    .Te0 = Te0,
    .n0 = n0,
    .B_p = B_p,
    .beta = beta,
    .tau = tau,
    .Ti0 = Ti0,
    .kperpRhos = kperpRhos,
    .logLambdaElc = logLambdaElc,
    .nuElc = nuElc,
    .logLambdaIon = logLambdaIon,
    .nuIon = nuIon,
    .vti = vti,
    .vte = vte,
    .c_s = c_s,
    .omega_ci = omega_ci,
    .rho_s = rho_s,
    .kperp = kperp,
    .z_min = z_min,
    .z_max = z_max,
    .psi_eval = psi_eval,
    .z_m = z_m,
    .n_m = n_m,
    .Te_m = Te_m,
    .Ti_m = Ti_m,
    .Ti_perp0 = Ti_perp0,
    .Ti_par0 = Ti_par0,
    .Ti_perp_m = Ti_perp_m,
    .Ti_par_m = Ti_par_m,
    .Te_par0 = Te_par0,
    .Te_par_m = Te_par_m,
    .Te_perp0 = Te_perp0,
    .Te_perp_m = Te_perp_m,
    .cs_m = cs_m,
    .NSrcIon = NSrcIon,
    .lineLengthSrcIon = lineLengthSrcIon,
    .sigSrcIon = sigSrcIon,
    .NSrcFloorIon = NSrcFloorIon,
    .TSrc0Ion = TSrc0Ion,
    .TSrcFloorIon = TSrcFloorIon,
    .NSrcElc = NSrcElc,
    .lineLengthSrcElc = lineLengthSrcElc,
    .sigSrcElc = sigSrcElc,
    .NSrcFloorElc = NSrcFloorElc,
    .TSrc0Elc = TSrc0Elc,
    .TSrcFloorElc = TSrcFloorElc,
    .vpar_max_ion = vpar_max_ion,
    .vpar_max_elc = vpar_max_elc,
    .mu_max_ion = mu_max_ion,
    .mu_max_elc = mu_max_elc,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nz, Nvpar, Nmu},
    .poly_order = poly_order,
    .t_end = t_end,
    .num_frames = num_frames,
    .write_phase_freq = write_phase_freq,
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };
  return ctx;
}


int main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) MPI_Init(&argc, &argv);
#endif

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct gk_mirror_ctx ctx = create_ctx(); // Context for init functions.

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.qe,
    .mass = ctx.me,
    .vdim = ctx.vdim,
    .lower = {-1.0, 0.0},
    .upper = { 1.0, 1.0},
    .cells = { cells_v[0], cells_v[1] },

    .polarization_density = ctx.n0,

    .mapc2p = {
      .mapping = mapc2p_vel_elc,
      .ctx = &ctx,
    },

    .projection = {
      .proj_id = GKYL_PROJ_BIMAXWELLIAN,
      .ctx_density = &ctx,
      .density = eval_density_elc,
      .ctx_upar = &ctx,
      .upar= eval_upar_elc,
      .ctx_temppar = &ctx,
      .temppar = eval_temp_par_elc,
      .ctx_tempperp = &ctx,
      .tempperp = eval_temp_perp_elc,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .den_ref = ctx.n0,
      .temp_ref = ctx.Te0,
      .num_cross_collisions = 1,
      .collide_with = { "ion" },
      .write_diagnostics = true,
    },

    .damping = {
      .type = GKYL_GK_DAMPING_LOW_PASS_FILTER,
      .rate_const = 1e-6,
      .write_rate = true,
      .write_fbar = true,
      .cellwise_const = true,
    },
    
    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
        .ctx_density = &ctx,
        .density = eval_density_elc_source,
        .ctx_upar = &ctx,
        .upar= eval_upar_elc_source,
        .ctx_temp = &ctx,
        .temp = eval_temp_elc_source,
      },
      .diagnostics = {
        .num_diag_moments = 5,
        .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP },
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
      }
    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH, },
    },

    .write_omega_cfl = true,
    .num_diag_moments = 8,
    .diag_moments = {GKYL_F_MOMENT_BIMAXWELLIAN, GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP },
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    .time_rate_diagnostics = true,

    .boundary_flux_diagnostics = {
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    },
  };

  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.qi,
    .mass = ctx.mi,
    .vdim = ctx.vdim,
    .lower = {-1.0, 0.0},
    .upper = { 1.0, 1.0},
    .cells = { cells_v[0], cells_v[1] },

    .polarization_density = ctx.n0,

    .mapc2p = {
      .mapping = mapc2p_vel_ion,
      .ctx = &ctx,
    },

    .projection = {
      .proj_id = GKYL_PROJ_BIMAXWELLIAN,
      .ctx_density = &ctx,
      .density = eval_density_ion,
      .ctx_upar = &ctx,
      .upar= eval_upar_ion,
      .ctx_temppar = &ctx,
      .temppar = eval_temp_par_ion,
      .ctx_tempperp = &ctx,
      .tempperp = eval_temp_perp_ion,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .den_ref = ctx.n0,
      .temp_ref = ctx.Ti0,
      .num_cross_collisions = 1,
      .collide_with = { "elc" },
      .write_diagnostics = true,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
        .ctx_density = &ctx,
        .density = eval_density_ion_source,
        .ctx_upar = &ctx,
        .upar= eval_upar_ion_source,
        .ctx_temp = &ctx,
        .temp = eval_temp_ion_source,
      },
      .diagnostics = {
        .num_diag_moments = 5,
        .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP },
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
      }
    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH, },
    },

    .write_omega_cfl = true,
    .num_diag_moments = 8,
    .diag_moments = {GKYL_F_MOMENT_BIMAXWELLIAN, GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP },
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    .time_rate_diagnostics = true,

    .boundary_flux_diagnostics = {
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    },
  };

  struct gkyl_gyrokinetic_field field = {
    .polarization_bmag = ctx.B_p,
    .kperpSq = pow(ctx.kperp, 2.),
    .time_rate_diagnostics = true,
  };

  struct gkyl_mirror_geo_grid_inp grid_inp = {
    .filename_psi = "gyrokinetic/data/unit/wham_hires.geqdsk_psi.gkyl", // psi file to use
    .rclose = 0.2, // closest R to region of interest
    .zmin = -2.0,  // Z of lower boundary
    .zmax =  2.0,  // Z of upper boundary
    .include_axis = false, // Include R=0 axis in grid
    .fl_coord = GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z, // coordinate system for psi grid
  };

  // GK app
  struct gkyl_gk app_inp = {
    .cdim = ctx.cdim,
    .lower = {ctx.z_min},
    .upper = {ctx.z_max},
    .cells = { cells_x[0] },
    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,

    .geometry = {
      .geometry_id = GKYL_GEOMETRY_MIRROR,
      .world = {ctx.psi_eval, 0.0},
      .mirror_grid_info = grid_inp,
    },

    .num_periodic_dir = 0,
    .periodic_dirs = {},

    .num_species = 2,
    .species = {elc, ion},

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
