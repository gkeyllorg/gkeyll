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
struct gk_mirror_ctx {
  int cdim, vdim; // Dimensionality.
  // Plasma parameters
  double mi;
  double qi;
  double me;
  double qe;
  double Te0;
  double n0;
  double B_p; // Hardcoded magnetic field at midplane 0.53 for double lorentzian
  double Bmag_midp; // Magnetic field magnitude at midplane (Z=0)
  double beta;
  double tau;
  double Ti0;
  double nuFrac;
  // Ion-ion collision freq.
  double logLambdaIon;
  double nuIon;
  double vti, vte;
  double RatZeq0; // Radius of the field line at Z=0.
  double kperp; // Perpendicular wavenumber for ES gyrokinetics.
  // Axial coordinate Z extents. Endure that Z=0 is not on
  double z_min;
  double z_max;
  double psi_eval;
  double psi_max;
  double psi_min;
  double theta_eval;
  double theta_min;
  double theta_max;
  // Physics parameters at mirror throat
  double vpar_max_ion;
  double mu_max_ion;
  double vpar_max_elc;
  double mu_max_elc;

  int Npsi;
  int Ntheta;
  int Nz;
  int Nvpar;
  int Nmu;
  int Nvpar_elc;
  int Nmu_elc;
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  int poly_order;

  // Source parameters
  double ion_source_amplitude;
  double ion_source_sigma;
  double ion_source_temp;

  double t_end; // End time.
  int num_frames; // Number of output frames.
  int num_phases; // Number of phases.
  double
    write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  double
    int_diag_calc_freq; // Frequency of calculating integrated diagnostics (as a factor of num_frames).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.

  // Geometry parameters for Lorentzian mirror
  double mcB; // Magnetic field parameter
  double mcB_inner; // Magnetic field parameter for inner mirror
  double mcB_outer; // Magnetic field parameter for outer mirror
  double gamma; // Width parameter for Lorentzian profile
  double gamma_inner; // Width parameter for Lorentzian profile for inner mirror
  double gamma_outer; // Width parameter for Lorentzian profile for outer mirror
  double Z_m; // Mirror throat location
  double L_center; // Length of central region
  double L_plugs; // Length of plugs
  double Z_min; // Minimum Z coordinate
  double Z_max; // Maximum Z coordinate
  double psi_in; // Working variable for psi integration
  double z_in; // Working variable for z integration
};

// Evaluate initial conditions
// I think ICs should be constructed to match the expander rather than the center
// as the center is quickly filled in the OAP
void initial_density(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  fout[0] = 1e17;
}

void initial_upar(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double z = xn[1];
  double c_s = 7 * sqrt(app->Te0 / app->mi);
  if (fabs(z) <= app->Z_m) {
    fout[0] = 0.0;
  } else {
    fout[0] = fabs(z) / z * c_s * tanh(4 * (app->Z_max - app->Z_m) * fabs(fabs(z) - app->Z_m));
  }
}

void eval_zero(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void initial_temp_ion(
  double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx
)
{
  struct gk_mirror_ctx *app = ctx;
  fout[0] = app->Ti0 / 10.0;
}

void initial_temp_elc(
  double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx
)
{
  struct gk_mirror_ctx *app = ctx;
  fout[0] = app->Te0;
}

void mapc2p_vel_ion(double t, const double *vc, double *GKYL_RESTRICT vp, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double vpar_max_ion = app->vpar_max_ion;
  double mu_max_ion = app->mu_max_ion;

  double cvpar = vc[0], cmu = vc[1];
  double b = 1.4;
  vp[0] = vpar_max_ion * tan(cvpar * b) / tan(b);
  vp[1] = mu_max_ion * pow(cmu, 2); // Cubic map in mu.
}

void mapc2p_vel_elc(double t, const double *vc, double *GKYL_RESTRICT vp, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double mu_max_elc = app->mu_max_elc;
  double vpar_max_elc = app->vpar_max_elc;
  double cvpar = vc[0], cmu = vc[1];
  vp[0] = vpar_max_elc * cvpar;
  vp[1] = mu_max_elc * pow(cmu, 4); // Cubic map in mu.
}

struct gk_mirror_ctx create_ctx(void)
{
  int cdim = 2, vdim = 2; // Dimensionality.
  int poly_order = 1;

  // Universal constant parameters.
  double eps0 = GKYL_EPSILON0;
  double mu0 = GKYL_MU0; // Not sure if this is right
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mp = GKYL_PROTON_MASS; // ion mass
  double me = GKYL_ELECTRON_MASS;
  double qi = eV; // ion charge
  double qe = -eV; // electron charge

  // Plasma parameters.
  double mi = 2.014 * mp;
  double Te0 = 940 * eV;
  double n0 = 3e19;
  double B_p = 0.53; // Bmag at z=0
  double beta = 0.4;
  double tau = pow(B_p, 2.) * beta / (2.0 * mu0 * n0 * Te0) - 1.;
  double Ti0 = tau * Te0;

  // Ion-ion collision freq.
  double nuFrac = 1.0;
  double logLambdaIon = 6.6 - 0.5 * log(n0 / 1e20) + 1.5 * log(Ti0 / eV);
  double nuIon = nuFrac * logLambdaIon * pow(eV, 4.) * n0 /
                 (12 * pow(M_PI, 3. / 2.) * pow(eps0, 2.) * sqrt(mi) * pow(Ti0, 3. / 2.));

  // Thermal speeds.
  double vti = sqrt(Ti0 / mi);
  double vte = sqrt(Te0 / me);

  // Grid parameters
  double vpar_max_ion = 16 * vti;
  double mu_max_ion = mi * pow(3. * vti, 2.) / (2. * B_p);
  double vpar_max_elc = 4 * vte;
  double mu_max_elc = me * pow(4. * vte, 2.) / (2. * B_p);

  int Nz = 32;
  int Npsi = 2;
  int Nvpar = 16;
  int Nmu = 8;
  int Nvpar_elc = 8;
  int Nmu_elc = 8;

  // Geometry parameters.
  double RatZeq0 = 0.10; // Radius of the field line at Z=0.
  double Z_min = -2.5;
  double Z_max = 2.5;
  double mcB = 3.691260;
  double gamma = 0.226381;
  double Z_m = 0.98;

  // Calculate phase structure
  double write_phase_freq =
    1; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  double int_diag_calc_freq =
    100; // Frequency of calculating integrated diagnostics (as a factor of num_frames).
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gk_mirror_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .poly_order = poly_order,

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

    .nuFrac = nuFrac,
    .logLambdaIon = logLambdaIon,
    .nuIon = nuIon,
    .vti = vti,
    .vte = vte,
    .RatZeq0 = RatZeq0,
    .vpar_max_ion = vpar_max_ion,
    .mu_max_ion = mu_max_ion,
    .vpar_max_elc = vpar_max_elc,
    .mu_max_elc = mu_max_elc,

    .Npsi = Npsi,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .Nvpar_elc = Nvpar_elc,
    .Nmu_elc = Nmu_elc,
    .cells = {Npsi, Nz, Nvpar, Nmu},
    .Z_min = Z_min,
    .Z_max = Z_max,

    .write_phase_freq = write_phase_freq,
    .int_diag_calc_freq = int_diag_calc_freq,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,

    .mcB = mcB,
    .gamma = gamma,
    .Z_m = Z_m
  };

  // Populate a couple more values in the context.
  ctx.psi_max = 3e-3;
  ctx.psi_min = 1e-5;
  ctx.psi_eval = (ctx.psi_max + ctx.psi_min) / 2.0;

  ctx.z_min = -2.5;
  ctx.z_max = 2.5;

  return ctx;
}

int main(int argc, char **argv)
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

  struct gk_mirror_ctx ctx = create_ctx(); // Context for init functions.

  int rank = 0;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  }
#endif
  // if (rank == 0)
  //   print_ctx(&ctx);

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d = 0; d < ctx.cdim; d++) {
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  }
  for (int d = 0; d < ctx.vdim; d++) {
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim + d]);
  }

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.qi,
    .mass = ctx.mi,
    .vdim = ctx.vdim,
    .lower = {-1.0, 0.0},
    .upper = {1.0, 1.0},
    .cells = {cells_v[0], cells_v[1]},
    .polarization_density = ctx.n0,

    .projection =
      {.proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
       .density = initial_density,
       .ctx_density = &ctx,
       .upar = initial_upar,
       .ctx_upar = &ctx,
       .temp = initial_temp_ion,
       .ctx_temp = &ctx},

    .mapc2p = {.mapping = mapc2p_vel_ion, .ctx = &ctx},

    .collisionless = {.type = GKYL_GK_COLLISIONLESS_ES, .write_diagnostics = true},

    .collisions =
      {.collision_id = GKYL_LBO_COLLISIONS,
       .den_ref = ctx.n0,
       .temp_ref = ctx.Ti0,
       .num_cross_collisions = 1,
       .collide_with = {"elc"},
       .write_diagnostics = true},

    .write_omega_cfl = true,
    .num_diag_moments = 8,
    .diag_moments =
      {GKYL_F_MOMENT_BIMAXWELLIAN, GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2,
       GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP},
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = {GKYL_F_MOMENT_M0M1M2PARM2PERP},
    .time_rate_diagnostics = true,

    .boundary_flux_diagnostics =
      {.num_integrated_diag_moments = 1, .integrated_diag_moments = {GKYL_F_MOMENT_M0M1M2PARM2PERP}}
  };

  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.qe,
    .mass = ctx.me,
    .vdim = ctx.vdim,
    .lower = {-1.0, 0.0},
    .upper = {1.0, 1.0},
    .cells = {ctx.Nvpar_elc, ctx.Nmu_elc},

    .polarization_density = ctx.n0,

    .mapc2p = {.mapping = mapc2p_vel_elc, .ctx = &ctx},

    .projection =
      {.proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
       .density = initial_density,
       .ctx_density = &ctx,
       .upar = eval_zero,
       .ctx_upar = &ctx,
       .temp = initial_temp_elc,
       .ctx_temp = &ctx,
       .correct_all_moms = true},

    .collisions =
      {.collision_id = GKYL_LBO_COLLISIONS,
       .den_ref = ctx.n0,
       .temp_ref = ctx.Te0,
       .num_cross_collisions = 1,
       .collide_with = {"ion"},
       .write_diagnostics = true,
       .not_in_dfdt = true},

    .num_diag_moments = 1,
    .diag_moments = {GKYL_F_MOMENT_MAXWELLIAN}
  };

  struct gkyl_gyrokinetic_field field = {
    .gkfield_id = GKYL_GK_FIELD_BOLTZMANN,
    .electron_mass = ctx.me,
    .electron_charge = ctx.qe,
    .electron_temp = ctx.Te0,
    .is_static = false
  };

  struct gkyl_mirror_geo_grid_inp grid_inp = {
    .filename_psi = "gyrokinetic/data/eqdsk/lorentzian_R10.geqdsk_psi.gkyl", // psi file to use
    .rclose = 0.2, // closest R to region of interest
    .zmin = -2.5, // Z of lower boundary
    .zmax = 2.5, // Z of upper boundary
    .include_axis = false, // Include R=0 axis in grid
    .fl_coord = GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z // coordinate system for psi grid
  };

  struct gkyl_gk app_inp = {
    // GK app
    .name = "gk_lorentzian_mirror",
    .cdim = ctx.cdim,
    .upper = {ctx.psi_max, ctx.Z_max},
    .lower = {ctx.psi_min, ctx.Z_min},
    .cells = {cells_x[0], cells_x[1]},
    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,

    .geometry = {.geometry_id = GKYL_GEOMETRY_MIRROR, .world = {0.0}, .mirror_grid_info = grid_inp},

    .num_periodic_dir = 0,
    .periodic_dirs = {},

    .num_species = 2,
    .species = {ion, elc},

    .field = field,

    .parallelism =
      {.use_gpu = app_args.use_gpu, .cuts = {app_args.cuts[0], app_args.cuts[1]}, .comm = comm}
  };
  struct gkyl_gyrokinetic_run_inp run_inp = {
    .app_inp = app_inp,
    .time_stepping =
      {.t_end = 5e-9,
       .num_frames = 50,
       .write_phase_freq = 1,
       .int_diag_calc_num = 50 * 100,
       .dt_failure_tol = ctx.dt_failure_tol,
       .num_failures_max = ctx.num_failures_max,
       .is_restart = false,
       .restart_frame = 0,
       .num_steps = app_args.num_steps},
    .print_verbosity =
      {.disable_timings = true, .estimate_completion_time = true, .enabled = true, .frequency = 1}
  };

  gkyl_gyrokinetic_run_simulation(&run_inp);

  gkyl_gyrokinetic_comms_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif
  return 0;
}
