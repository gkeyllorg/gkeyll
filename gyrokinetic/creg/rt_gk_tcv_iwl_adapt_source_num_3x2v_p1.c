#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_eqn_type.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_gen_miller_eqdsk.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_gyrokinetic_run.h>
#include <gkyl_math.h>

#include <rt_arg_parse.h>

// Numerical-equilibrium counterpart of rt_gk_tcv_iwl_adapt_source_3x2v_p1.c.
// Computational coordinates are (psi, alpha, normalized poloidal arc length).
// Keep the original outer surface and put the LCFS at 2/5 of the uniform psi grid;
// this slightly changes the inner surface. Nx must therefore be divisible by 5.
// Density, temperature, and Gaussian envelopes are evaluated in the original minor-radius
// coordinate and Miller angle using generator helpers. These profiles are independent of
// alpha; only the binormal domain width needs conversion by Cy.

// Define the context of the simulation. This stores global parameters.
struct gk_app_ctx {
  int cdim, vdim;
  // Geometry and magnetic field parameters
  double a_shift, Zvert_axis, Rmajor_axis, R0, a_mid, x_inner, r0, B0, kappa, delta, q0, Cy, Bref,
    x_LCFS;
  struct gkyl_gen_miller_eqdsk_inp miller;
  gkyl_gen_miller_eqdsk *equilibrium; // Flux table for profile callbacks; owned by main.
  // Plasma parameters
  int num_species;
  double me, qe, mi, qi, n0, Te0, Ti0;
  // Collision parameters
  double nuFrac, nuElc, nuIon;
  // Source parameters
  int num_sources;
  bool adapt_energy_srcCORE, adapt_particle_srcCORE;
  double center_srcCORE[3], sigma_srcCORE[3];
  double energy_srcCORE, particle_srcCORE;
  double floor_srcCORE;
  bool adapt_energy_srcRECY, adapt_particle_srcRECY;
  double center_srcRECY[3], sigma_srcRECY[3];
  double energy_srcRECY, particle_srcRECY;
  double floor_srcRECY;
  // Grid parameters
  double Lx, Ly, Lz;
  double x_min, x_max, y_min, y_max, z_min, z_max;
  int Nx, Ny, Nz, Nvpar, Nmu;
  int cells[GKYL_MAX_DIM], poly_order;
  double vpar_max_elc, mu_max_elc, vpar_max_ion, mu_max_ion;
  // Simulation control parameters
  double t_end, write_phase_freq;
  int num_frames, int_diag_calc_num, num_failures_max;
  double dt_failure_tol;
};

static double
qsafety(double R, void *ctx)
{
  // Evaluate the cubic fit to the TCV negative-triangularity q profile (discharge 65130).
  return ((497.3420166252413 * R - 1408.736172826569) * R + 1331.4134861681464) * R -
         419.00692601227627;
}

/**
 * Convert the original input's radial coordinate to poloidal flux.
 *
 * @param app Simulation context with an initialized flux table.
 * @param x Original radial coordinate, r - a_mid + x_inner [m].
 * @return Poloidal flux divided by 2*pi [Wb/rad].
 */
static double
psi_at_x(const struct gk_app_ctx *app, double x)
{
  return gkyl_gen_miller_eqdsk_psi_at_r(app->equilibrium, x + app->a_mid - app->x_inner);
}

/**
 * Recover the original radial coordinate for the density and temperature fits.
 *
 * @param app Simulation context with an initialized flux table.
 * @param psi Poloidal flux divided by 2*pi [Wb/rad].
 * @return Original radial coordinate, r - a_mid + x_inner [m].
 */
static double
x_at_psi(const struct gk_app_ctx *app, double psi)
{
  return gkyl_gen_miller_eqdsk_r_at_psi(app->equilibrium, psi) - app->a_mid + app->x_inner;
}

struct source_profile_ctx {
  const struct gk_app_ctx *app;
  double center[3], sigma[3], floor;
};

static void
source_profile(double t, const double *xn, double *fout, void *ctx)
{
  // A source Gaussian in the analytic input's radial coordinate and Miller angle.
  const struct source_profile_ctx *src = ctx;
  double x = x_at_psi(src->app, xn[0]);
  double theta = gkyl_gen_miller_eqdsk_theta_at_theta_tok(src->app->equilibrium, xn[2]);
  double dx = (x - src->center[0]) / src->sigma[0];
  // The analytic input wraps parallel Gaussian distances in the IWL domain.
  double dtheta = remainder(theta - src->center[2], src->app->Lz) / src->sigma[2];
  fout[0] = exp(-0.5 * (dx * dx + dtheta * dtheta)) + src->floor;
}

static void
zero_func(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

static void
density_init(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  // Density fit at the minor radius corresponding to psi.
  struct gk_app_ctx *app = ctx;
  double x = x_at_psi(app, xn[0]);
  double n0 = 5e19;
  double x0 = -0.03;
  double c1 = 0.5;
  double c2 = 8.0;
  double c3 = 0.005;
  fout[0] = n0 * (c1 * (1. + tanh(c2 * (-10 * (x + x0)))) + c3);
}

static void
temp_elc(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  // Electron temperature fit at the minor radius corresponding to psi.
  struct gk_app_ctx *app = ctx;
  double x = x_at_psi(app, xn[0]);
  double T0 = 200 * GKYL_ELEMENTARY_CHARGE;
  double x0 = -0.03; // position of the transition region
  double c0 = 1.3; // multiplicative factor
  double c1 = 0.5; // control the temperature at the _core
  double c2 = 8.0; // control the width of the transition region
  double c3 = 0.1; // control the temperature at the SOL
  fout[0] = c0 * T0 * (c1 * (1. + tanh(c2 * (-10 * (x + x0)))) + c3);
}

static void
temp_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  // Ion temperature fit at the minor radius corresponding to psi.
  struct gk_app_ctx *app = ctx;
  double x = x_at_psi(app, xn[0]);
  double T0 = 200 * GKYL_ELEMENTARY_CHARGE;
  double x0 = -0.04; // position of the transition region
  double c0 = 1.0; // multiplicative factor
  double c1 = 0.5; // control the temperature at the _core
  double c2 = 3.0; // control the width of the transition region
  double c3 = 0.2; // control the temperature at the SOL
  fout[0] = c0 * T0 * (c1 * (1. + tanh(c2 * (-10 * (x + x0)))) + c3);
}

static void
mapc2p_vel_elc(double t, const double *vc, double *GKYL_RESTRICT vp, void *ctx)
{
  // Map electron velocity coordinates to physical parallel velocity and magnetic moment.
  struct gk_app_ctx *app = ctx;
  double vpar_max_elc = app->vpar_max_elc;
  double mu_max_elc = app->mu_max_elc;
  double cvpar = vc[0], cmu = vc[1];
  // Linear map up to vpar_max/2, then quadratic.
  if (fabs(cvpar) <= 0.5) {
    vp[0] = vpar_max_elc * cvpar;
  } else if (cvpar < -0.5) {
    vp[0] = -vpar_max_elc * 2.0 * pow(cvpar, 2);
  } else {
    vp[0] = vpar_max_elc * 2.0 * pow(cvpar, 2);
  }
  // Quadratic map in mu.
  vp[1] = mu_max_elc * pow(cmu, 2);
}

static void
mapc2p_vel_ion(double t, const double *vc, double *GKYL_RESTRICT vp, void *ctx)
{
  // Map ion velocity coordinates to physical parallel velocity and magnetic moment.
  struct gk_app_ctx *app = ctx;
  double vpar_max_ion = app->vpar_max_ion;
  double mu_max_ion = app->mu_max_ion;
  double cvpar = vc[0], cmu = vc[1];
  // Linear map up to vpar_max/2, then quadratic.
  if (fabs(cvpar) <= 0.5) {
    vp[0] = vpar_max_ion * cvpar;
  } else if (cvpar < -0.5) {
    vp[0] = -vpar_max_ion * 2.0 * pow(cvpar, 2);
  } else {
    vp[0] = vpar_max_ion * 2.0 * pow(cvpar, 2);
  }
  // Quadratic map in mu.
  vp[1] = mu_max_ion * pow(cmu, 2);
}

static struct gk_app_ctx
create_ctx(void)
{
  int cdim = 3, vdim = 2; // Dimensionality.
  // Universal constant parameters.
  double eps0 = GKYL_EPSILON0, eV = GKYL_ELEMENTARY_CHARGE;
  double mp = GKYL_PROTON_MASS, me = GKYL_ELECTRON_MASS;
  double qi = eV; // ion charge
  double qe = -eV; // electron charge

  // Geometry and magnetic field.
  double a_shift = 0.5; // Parameter in Shafranov shift.
  double Zvert_axis = 0.1414361745; // Magnetic axis height [m], passed explicitly to the generator.
  double Rmajor_axis = 0.8867856264; // Magnetic axis major radius [m].
  double bmag_axis = 1.4; // Magnetic field at the magnetic axis [T].
  double R_LCFSmid = 1.0870056099999; // Major radius of the LCFS at the outboard midplane [m].
  double x_inner = 0.04; // Radial extent inside LCFS
  double x_outer = 0.08; // Radial extent outside LCFS
  double Rmid_min = R_LCFSmid - x_inner; // Minimum midplane major radius of simulation box [m].
  double Rmid_max = R_LCFSmid + x_outer; // Maximum midplane major radius of simulation box [m].
  double R0 = 0.5 * (Rmid_min + Rmid_max); // Major radius of the simulation box [m].

  // Minor radius at outboard midplane [m]. Redefine it with
  // Shafranov shift, to ensure LCFS radial location.
  double a_mid =
    fabs(a_shift) < 1e-13 ?
      R_LCFSmid - Rmajor_axis :
      Rmajor_axis / a_shift -
        sqrt(Rmajor_axis * (Rmajor_axis - 2 * a_shift * R_LCFSmid + 2 * a_shift * Rmajor_axis)) /
          a_shift;

  double r0 = R0 - Rmajor_axis; // Minor radius of the simulation box [m].
  double B0 = bmag_axis * (Rmajor_axis / R0); // Magnetic field magnitude in the simulation box [T].
  double kappa = 1.4; // Elongation (=1 for no elongation).
  double delta = -0.38; // Triangularity (=0 for no triangularity).

  // Plasma parameters. Chosen based on the value of a cubic sline
  // between the last TS data inside the LCFS and the probe data in
  // in the far SOL, near R=0.475 m.
  int num_species = 2;
  double AMU = 2.01410177811;
  double mi = mp * AMU; // Deuterium ions.
  double Te0 = 100 * eV; // Electron reference temperature [J].
  double Ti0 = 100 * eV; // Ion reference temperature [J].
  double n0 = 2.0e19; // Reference density [1/m^3].
  double Bref = 1.129; // Reference magnetic field [T].
  double vte = sqrt(Te0 / me), vti = sqrt(Ti0 / mi);
  double c_s = sqrt(Te0 / mi);
  double omega_ci = fabs(qi * B0 / mi);
  double rho_s = c_s / omega_ci;

  // Configuration domain parameters
  double Lx = Rmid_max - Rmid_min; // Domain size along x.
  double x_min = 0.;
  double x_max = Lx;
  double x_LCFS = R_LCFSmid - Rmid_min; // Radial location of the last closed flux surface.
  double q0 = qsafety(R0, 0); // Safety factor in the center of domain.
  double Cy = r0 / q0; // Normalization in binormal coordinate.

  double Ly = 150 * rho_s; // Domain size along y.
  // Adjust the domain size along y to have integer toroidal mode number.
  // We need: 2*pi*Cy/Ly = integer.
  Ly = 2. * M_PI * Cy / round(2. * M_PI * Cy / Ly);
  double y_min = -Ly / 2.;
  double y_max = Ly / 2.;

  double Lz = 2. * M_PI - 1e-10; // Domain size along magnetic field.
  double z_min = -Lz / 2.;
  double z_max = Lz / 2.;

  // Collision frequencies
  double nuFrac = 0.1;
  // Electron-electron collision freq.
  double logLambdaElc = 6.6 - 0.5 * log(n0 / 1e20) + 1.5 * log(Ti0 / eV);
  double nuElc = nuFrac * logLambdaElc * pow(eV, 4) * n0 /
                 (6 * sqrt(2.) * pow(M_PI, 3. / 2.) * pow(eps0, 2) * sqrt(me) * pow(Te0, 3. / 2.));
  // Ion-ion collision freq.
  double logLambdaIon = 6.6 - 0.5 * log(n0 / 1e20) + 1.5 * log(Ti0 / eV);
  double nuIon = nuFrac * logLambdaIon * pow(eV, 4) * n0 /
                 (12 * pow(M_PI, 3. / 2.) * pow(eps0, 2) * sqrt(mi) * pow(Ti0, 3. / 2.));

  // Source parameters
  int num_sources = 2;
  double P_exp = 0.34e6; // P_sol measured [W]
  double vol_frac = 1.0 / (2. * M_PI * Cy / Ly); // Toroidal fraction of the simulation box.
  double P_inj =
    P_exp * vol_frac /
    num_species; // Injection power normalized to the volume fraction and per species [W]
  // Core source:
  // - Injects the prescribed power in the core region.
  // - Replaces particles lost through the inner radial boundary.
  bool adapt_energy_srcCORE =
    true; // The source will compensate the losses in energy according to given boundaries.
  bool adapt_particle_srcCORE =
    true; // The source will compensate the losses in particle according to given boundaries.
  double energy_srcCORE = P_inj; // What the source must inject in energy [W]
  double particle_srcCORE = 0.0; // What the source must inject in particle [1/s]
  double center_srcCORE[3] = {x_min, 0.0, -Lz / 4}; // This is the position of the ion source,
  double sigma_srcCORE[3] = {0.03 * Lx, 0.0, Lz / 6}; // the electron source will be at +Lz/4.
  double floor_srcCORE = 1e-10;
  // Recycling source:
  // - Reinjects particles that are absorbed by the wall.
  // - Energy is free to leave the system.
  bool adapt_energy_srcRECY = false;
  bool adapt_particle_srcRECY = true;
  double energy_srcRECY = 0.0; // [W]
  double particle_srcRECY = 0.0; // [1/s]
  double center_srcRECY[3] = {0.5 * x_LCFS, 0.0, M_PI};
  double sigma_srcRECY[3] = {0.25 * x_LCFS, 0.0, 0.05 * Lz};
  double floor_srcRECY = 1e-10;

  // Grid parameters (reduced resolution for the regression test, minimal recommended values in comments)
  int Nx = 15; // LCFS at 2/5 of the psi domain; must be divisible by 5.
  int Ny = 8; // (16)
  int Nz = 8; // (12)
  int Nvpar = 4; // (12)
  int Nmu = 4; // (8)
  int poly_order = 1;
  // Velocity box dimensions
  double vpar_max_elc = 5. * vte;
  double mu_max_elc = 1. * me * pow(4 * vte, 2) / (2 * B0);
  double vpar_max_ion = 5. * vti;
  double mu_max_ion = 1. * mi * pow(4 * vti, 2) / (2 * B0);
  double t_end = 1.e-6;
  int num_frames = 1;
  double write_phase_freq = 1.0;
  int int_diag_calc_num = num_frames * 100;
  double dt_failure_tol = 1.0e-3; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gk_app_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .a_shift = a_shift,
    .Rmajor_axis = Rmajor_axis,
    .Zvert_axis = Zvert_axis,
    .R0 = R0,
    .a_mid = a_mid,
    .x_inner = x_inner,
    .r0 = r0,
    .B0 = B0,
    .kappa = kappa,
    .delta = delta,
    .q0 = q0,
    .Cy = Cy,
    .Lx = Lx,
    .Ly = Ly,
    .Lz = Lz,
    .x_min = x_min,
    .x_max = x_max,
    .y_min = y_min,
    .y_max = y_max,
    .z_min = z_min,
    .z_max = z_max,
    .Bref = Bref,
    .x_LCFS = x_LCFS,
    .num_species = num_species,
    .me = me,
    .qe = qe,
    .mi = mi,
    .qi = qi,
    .n0 = n0,
    .Te0 = Te0,
    .Ti0 = Ti0,
    .nuFrac = nuFrac,
    .nuElc = nuElc,
    .nuIon = nuIon,
    .num_sources = num_sources,
    .adapt_energy_srcCORE = adapt_energy_srcCORE,
    .adapt_particle_srcCORE = adapt_particle_srcCORE,
    .center_srcCORE = {center_srcCORE[0], center_srcCORE[1], center_srcCORE[2]},
    .sigma_srcCORE = {sigma_srcCORE[0], sigma_srcCORE[1], sigma_srcCORE[2]},
    .energy_srcCORE = energy_srcCORE,
    .particle_srcCORE = particle_srcCORE,
    .floor_srcCORE = floor_srcCORE,
    .adapt_energy_srcRECY = adapt_energy_srcRECY,
    .adapt_particle_srcRECY = adapt_particle_srcRECY,
    .center_srcRECY = {center_srcRECY[0], center_srcRECY[1], center_srcRECY[2]},
    .sigma_srcRECY = {sigma_srcRECY[0], sigma_srcRECY[1], sigma_srcRECY[2]},
    .energy_srcRECY = energy_srcRECY,
    .particle_srcRECY = particle_srcRECY,
    .floor_srcRECY = floor_srcRECY,
    .Nx = Nx,
    .Ny = Ny,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nx, Ny, Nz, Nvpar, Nmu},
    .poly_order = poly_order,
    .vpar_max_elc = vpar_max_elc,
    .mu_max_elc = mu_max_elc,
    .vpar_max_ion = vpar_max_ion,
    .mu_max_ion = mu_max_ion,
    .write_phase_freq = write_phase_freq,
    .t_end = t_end,
    .num_frames = num_frames,
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };
  // Use a fixed coverage envelope so both flux tables (here and in the app) agree.
  // R_omp includes the quadratic shift; Rmid_min/max above are nominal reference radii.
  double outer_r = a_mid + x_outer;
  double outer_R = Rmajor_axis - a_shift * outer_r * outer_r / (2 * Rmajor_axis) + outer_r;
  ctx.miller = (struct gkyl_gen_miller_eqdsk_inp){
    .x_extent_type = GKYL_GK_GEOMETRY_RADIAL_COORD_R_OMP,
    .Rmajor_axis = Rmajor_axis,
    .Zvert_axis = Zvert_axis,
    .bmag_axis = -fabs(bmag_axis), // Match the analytic input's negative toroidal field.
    .Rmajor_LCFS = R_LCFSmid,
    .kappa = kappa,
    .delta = delta,
    .a_shift = a_shift,
    .qsafety = qsafety,
    .has_Rmajor_omp_bounds = true,
    .Rmajor_omp_bounds = {Rmajor_axis, outer_R},
  };
  ctx.equilibrium = gkyl_gen_miller_eqdsk_inew(&ctx.miller);
  if (!ctx.equilibrium) {
    return ctx;
  }

  // Match the LCFS value written with E16.9 precision in the generated GEQDSK.
  char boundary_text[32];
  snprintf(
    boundary_text, sizeof boundary_text, "%.9E",
    gkyl_gen_miller_eqdsk_get_info(ctx.equilibrium)->psi_boundary
  );
  ctx.x_LCFS = strtod(boundary_text, 0); // Bias-line coordinates must be psi, not R_omp.
  double psi_upper = psi_at_x(&ctx, x_max);
  double psi_lower = (ctx.x_LCFS - 0.4 * psi_upper) / 0.6;
  ctx.x_min = gkyl_gen_miller_eqdsk_R_omp_at_psi(ctx.equilibrium, psi_lower);
  ctx.x_max = outer_R; // Miller app inputs use R_omp; the app converts them to psi.

  // Numerical tokamak y is alpha [rad], whereas the analytic input used y = Cy*alpha [m].
  // Power normalization above already uses this same toroidal fraction.
  ctx.Ly /= Cy;
  ctx.y_min /= Cy;
  ctx.y_max /= Cy;
  return ctx;
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

  struct gk_app_ctx ctx = create_ctx(); // context for init functions

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d = 0; d < ctx.cdim; d++) {
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  }
  for (int d = 0; d < ctx.vdim; d++) {
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim + d]);
  }

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  int failed = ctx.equilibrium == 0 || cells_x[0] < 5 || cells_x[0] % 5 != 0;
  int any_failed;
  gkyl_comm_allreduce_host(comm, GKYL_INT, GKYL_MAX, 1, &failed, &any_failed);
  if (any_failed) {
    fprintf(stderr, "TCV Miller input requires a valid equilibrium and Nx divisible by 5.\n");
    gkyl_gen_miller_eqdsk_release(ctx.equilibrium);
    gkyl_gyrokinetic_comms_release(comm);
#ifdef GKYL_HAVE_MPI
    if (app_args.use_mpi) {
      MPI_Finalize();
    }
#endif
    return 1;
  }

  int rank;
  gkyl_comm_get_rank(comm, &rank);
  if (rank == 0) {
    printf(
      "TCV Miller: R_omp=[%.12g, %.12g] m, psi_LCFS=%.12g Wb/rad (2/5 of radial grid).\n",
      ctx.x_min, ctx.x_max, ctx.x_LCFS
    );
  }

  struct source_profile_ctx core_elc = {
    .app = &ctx,
    .center = {ctx.center_srcCORE[0], 0, -ctx.center_srcCORE[2]},
    .sigma = {ctx.sigma_srcCORE[0], 0, ctx.sigma_srcCORE[2]},
    .floor = ctx.floor_srcCORE,
  };
  struct source_profile_ctx core_ion = core_elc;
  core_ion.center[2] = ctx.center_srcCORE[2];
  struct source_profile_ctx recycling = {
    .app = &ctx,
    .center = {ctx.center_srcRECY[0], 0, ctx.center_srcRECY[2]},
    .sigma = {ctx.sigma_srcRECY[0], 0, ctx.sigma_srcRECY[2]},
    .floor = ctx.floor_srcRECY,
  };

  // Definition of the sources projections.
  struct gkyl_gyrokinetic_projection proj_srcCORE_e = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_GAUSSIAN,
    .configuration_profile = source_profile,
    .configuration_profile_ctx = &core_elc,
    .total_num_particles = ctx.particle_srcCORE,
    .total_kin_energy = ctx.energy_srcCORE,
    .temp_max = 5.0 * ctx.Te0,
    .temp_min = 0.1 * ctx.Te0,
    .f_floor = ctx.floor_srcCORE,
  };
  struct gkyl_gyrokinetic_projection proj_srcCORE_i = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_GAUSSIAN,
    .configuration_profile = source_profile,
    .configuration_profile_ctx = &core_ion,
    .total_num_particles = ctx.particle_srcCORE,
    .total_kin_energy = ctx.energy_srcCORE,
    .temp_max = 5.0 * ctx.Ti0,
    .temp_min = 0.1 * ctx.Ti0,
    .f_floor = ctx.floor_srcCORE,
  };
  struct gkyl_gyrokinetic_projection proj_srcRECY_e = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_GAUSSIAN,
    .configuration_profile = source_profile,
    .configuration_profile_ctx = &recycling,
    .total_num_particles = ctx.particle_srcRECY,
    .total_kin_energy = ctx.energy_srcRECY,
    .temp_max = 5.0 * ctx.Te0,
    .temp_min = 0.1 * ctx.Te0,
    .f_floor = ctx.floor_srcRECY,
  };
  struct gkyl_gyrokinetic_projection proj_srcRECY_i = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_GAUSSIAN,
    .configuration_profile = source_profile,
    .configuration_profile_ctx = &recycling,
    .total_num_particles = ctx.particle_srcRECY,
    .total_kin_energy = ctx.energy_srcRECY,
    .temp_max = 5.0 * ctx.Ti0,
    .temp_min = 0.1 * ctx.Ti0,
    .f_floor = ctx.floor_srcRECY,
  };

  // Definition of the source adaptation structures.
  struct gkyl_gyrokinetic_adapt_source adapt_srcCORE_e = {
    .adapt_to_species = "elc",
    .adapt_particle = ctx.adapt_particle_srcCORE,
    .adapt_energy = ctx.adapt_energy_srcCORE,
    .num_boundaries = 1, // Only the inner radial boundary.
    .dir = {0},
    .edge = {GKYL_LOWER_EDGE},
  };
  struct gkyl_gyrokinetic_adapt_source adapt_srcCORE_i = {
    .adapt_to_species = "ion",
    .adapt_particle = ctx.adapt_particle_srcCORE,
    .adapt_energy = ctx.adapt_energy_srcCORE,
    .num_boundaries = 1, // Only the inner radial boundary.
    .dir = {0},
    .edge = {GKYL_LOWER_EDGE},
  };
  struct gkyl_gyrokinetic_adapt_source adapt_srcRECY_e = {
    .adapt_to_species = "ion", // Adapted to ions to maintain neutrality.
    .adapt_particle = ctx.adapt_particle_srcRECY,
    .adapt_energy = ctx.adapt_energy_srcRECY,
    .num_boundaries = 3, // outer radial boundary + both z boundaries.
    .dir = {0, 2, 2},
    .edge = {GKYL_UPPER_EDGE, GKYL_LOWER_EDGE, GKYL_UPPER_EDGE},
  };
  struct gkyl_gyrokinetic_adapt_source adapt_srcRECY_i = {
    .adapt_to_species = "ion",
    .adapt_particle = ctx.adapt_particle_srcRECY,
    .adapt_energy = ctx.adapt_energy_srcRECY,
    .num_boundaries = 3, // outer radial boundary + both z boundaries.
    .dir = {0, 2, 2},
    .edge = {GKYL_UPPER_EDGE, GKYL_LOWER_EDGE, GKYL_UPPER_EDGE},
  };

  // Electrons
  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.qe,
    .mass = ctx.me,
    .vdim = ctx.vdim,
    .lower = {-1.0 / sqrt(2.0), 0.0},
    .upper = {1.0 / sqrt(2.0), 1.0},
    .cells = {cells_v[0], cells_v[1]},
    .polarization_density = ctx.n0,

    .mapc2p = {.mapping = mapc2p_vel_elc, .ctx = &ctx},

    // .init_from_file = {
    //    .type = GKYL_IC_IMPORT_F,
    //    .file_name = "restart-elc.gkyl",
    //    .jacobtot_inv_file_name = "restart-jacobtot_inv.gkyl",
    // },

    .projection =
      {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
        .ctx_density = &ctx,
        .ctx_upar = &ctx,
        .ctx_temp = &ctx,
        .density = density_init,
        .upar = zero_func,
        .temp = temp_elc,
      },

    .collisionless = {.type = GKYL_GK_COLLISIONLESS_ES},

    .collisions =
      {
        .collision_id = GKYL_LBO_COLLISIONS,
        .nu_frac = ctx.nuFrac,
        .den_ref = ctx.n0, // Density used to calculate coulomb logarithm
        .temp_ref = ctx.Te0, // Temperature used to calculate coulomb logarithm
        .num_cross_collisions = 1,
        .collide_with = {"ion"},
      },

    .source =
      {
        .source_id = GKYL_PROJ_SOURCE,
        .num_sources = ctx.num_sources,
        .num_adapt_sources = ctx.num_sources,
        .projection[0] = proj_srcCORE_e,
        .adapt[0] = adapt_srcCORE_e,
        .projection[1] = proj_srcRECY_e,
        .adapt[1] = adapt_srcRECY_e,
        .diagnostics =
          {
            .num_diag_moments = 1,
            .diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
            .num_integrated_diag_moments = 1,
            .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
          },
      },

    .bcs =
      {{.dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
       {.dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
       {.dir = 2, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH_CONDUCTING},
       {.dir = 2, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH_CONDUCTING}},

    .num_diag_moments = 9,
    .diag_moments =
      {GKYL_F_MOMENT_HAMILTONIAN, GKYL_F_MOMENT_BIMAXWELLIAN, GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1,
       GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M3PAR,
       GKYL_F_MOMENT_M3PERP},

    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},

    .boundary_flux_diagnostics =
      {
        .num_diag_moments = 1,
        .diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
      },

    .time_rate_diagnostics = true,
  };

  // Ions
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.qi,
    .mass = ctx.mi,
    .vdim = ctx.vdim,
    .lower = {-1.0 / sqrt(2.0), 0.0},
    .upper = {1.0 / sqrt(2.0), 1.0},
    .cells = {cells_v[0], cells_v[1]},
    .polarization_density = ctx.n0,

    .mapc2p = {.mapping = mapc2p_vel_ion, .ctx = &ctx},

    // .init_from_file = {
    //    .type = GKYL_IC_IMPORT_F,
    //    .file_name = "restart-ion.gkyl",
    //    .jacobtot_inv_file_name = "restart-jacobtot_inv.gkyl",
    // },

    .projection =
      {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
        .ctx_density = &ctx,
        .ctx_upar = &ctx,
        .ctx_temp = &ctx,
        .density = density_init,
        .upar = zero_func,
        .temp = temp_ion,
      },

    .collisionless = {.type = GKYL_GK_COLLISIONLESS_ES},

    .collisions =
      {
        .collision_id = GKYL_LBO_COLLISIONS,
        .nu_frac = ctx.nuFrac,
        .den_ref = ctx.n0, // Density used to calculate coulomb logarithm
        .temp_ref = ctx.Ti0, // Temperature used to calculate coulomb logarithm
        .num_cross_collisions = 1,
        .collide_with = {"elc"},
      },

    .source =
      {
        .source_id = GKYL_PROJ_SOURCE,
        .num_sources = ctx.num_sources,
        .num_adapt_sources = ctx.num_sources,
        .projection[0] = proj_srcCORE_i,
        .adapt[0] = adapt_srcCORE_i,
        .projection[1] = proj_srcRECY_i,
        .adapt[1] = adapt_srcRECY_i,
        .diagnostics =
          {
            .num_diag_moments = 1,
            .diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
            .num_integrated_diag_moments = 1,
            .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
          },
      },

    .bcs =
      {{.dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
       {.dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
       {.dir = 2, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH_CONDUCTING},
       {.dir = 2, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH_CONDUCTING}},

    .num_diag_moments = 9,
    .diag_moments =
      {GKYL_F_MOMENT_HAMILTONIAN, GKYL_F_MOMENT_BIMAXWELLIAN, GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1,
       GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M3PAR,
       GKYL_F_MOMENT_M3PERP},

    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},

    .boundary_flux_diagnostics =
      {
        .num_diag_moments = 1,
        .diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
      },

    .time_rate_diagnostics = true,
  };

  struct gkyl_poisson_bias_line target_corner_bcs[] = {
    {
      .perp_dirs = {0, 2}, // Directions perpendicular to line.
      .perp_coords = {ctx.x_LCFS, ctx.z_min}, // Coordinates of the line in perpendicular directions.
      .val = 0.0, // Biasing value.
    },
    {
      .perp_dirs = {0, 2}, // Directions perpendicular to line.
      .perp_coords = {ctx.x_LCFS, ctx.z_max}, // Coordinates of the line in perpendicular directions.
      .val = 0.0, // Biasing value.
    }
  };

  struct gkyl_poisson_bias_line_list bias_line_list = {.num_bias_line = 2, .bl = target_corner_bcs};

  // field
  struct gkyl_gyrokinetic_field field = {
    .gkfield_id = GKYL_GK_FIELD_ES,
    .polarization_bmag = ctx.Bref,
    .poisson_bcs =
      {{.dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}},
       {.dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}}},
    .bias_line_list = &bias_line_list,
    .time_rate_diagnostics = true,
  };

  // The app writes <app-name>-miller.geqdsk and initializes numerical tokamak geometry.
  // It also derives the LCFS split and core twist shifts from that equilibrium.
  struct gkyl_gyrokinetic_geometry geometry = {
    .geometry_id = GKYL_GEOMETRY_MILLER,
    .miller_grid_info = ctx.miller,
  };

  // Parallelism
  struct gkyl_app_parallelism_inp parallelism = {
    .comm = comm,
    .cuts = {app_args.cuts[0], app_args.cuts[1], app_args.cuts[2]},
    .use_gpu = app_args.use_gpu,
  };

  // GK app
  struct gkyl_gk app_inp = {

    .cfl_frac_omegaH = 1.0e9,
    .cfl_frac = 1.0,

    .cdim = ctx.cdim,
    .lower = {ctx.x_min, ctx.y_min, ctx.z_min},
    .upper = {ctx.x_max, ctx.y_max, ctx.z_max},
    .cells = {cells_x[0], cells_x[1], cells_x[2]},
    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,

    .geometry = geometry,

    .num_periodic_dir = 1,
    .periodic_dirs = {1},

    .num_species = ctx.num_species,
    .species = {elc, ion},

    .field = field,

    .parallelism = parallelism,
  };

  // Set app output name from the executable name (argv[0]).
  snprintf(app_inp.name, sizeof(app_inp.name), "%s", app_args.app_name);

  struct gkyl_gyrokinetic_run_inp run_inp = {
    .app_inp = app_inp,
    .time_stepping =
      {
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

  gkyl_gen_miller_eqdsk_release(ctx.equilibrium);

  gkyl_gyrokinetic_comms_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif

  return 0;
}
