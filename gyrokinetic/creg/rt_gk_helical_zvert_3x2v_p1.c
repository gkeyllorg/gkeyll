/*
 * This regression test is based on Section 5.2, "Helical SOL configuration 
 * including magnetic shear", Thesis of Noah Mandell, 2021.
*/
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
struct gk_app_ctx {
  int cdim, vdim; // Dimensionality.
  
  // Geometry and magnetic field.
  double R_axis;
  double a0;
  double B_axis;
  double Bvx0;
  double x0;
  double R0;
  double H;
  int n; // Exponent for vertical magnetic field profile.

  // Plasma parameters.
  double me; double qe;
  double mi; double qi;
  double n0; double Te0; double Ti0; 

  // Collisions.
  double nu_frac;

  // Source parameters.
  double Psrc;
  double S0;
  double xSource;
  double lambdaSource;

  // Grid parameters.
  double Lx; // Domain size in radial direction.
  double Ly; // Domain size in binormal direction.
  double Lz; // Domain size along magnetic field.
  double x_min; double x_max;
  double y_min; double y_max;
  double z_min; double z_max;
  int Nx;
  int Ny;
  int Nz;
  int Nvpar;
  int Nmu;
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  int poly_order;
  double vpar_max_elc; double mu_max_elc;
  double vpar_max_ion; double mu_max_ion;

  double t_end; // End time.
  int num_frames; // Number of output frames.
  double write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
  double max_run_time; // Maximum run time in seconds, 0 means no limit.
};

// Common source density profiles.
double sourceDensity(double t, const double * GKYL_RESTRICT xn, void *ctx)
{
  double x = xn[0], z = xn[2];
  double sourceFloor = 0.1;

  struct gk_app_ctx *app = ctx;
  double S0 = app->S0;
  double lambdaSource = app->lambdaSource;
  double xSource = app->xSource;
  double Lz = app->Lz;

  if (x < xSource) {
    sourceFloor = 0.5;
  }
  if (fabs(z) < Lz/4) {
    return 0.90625*S0*fmax(exp(-0.5*pow((x-xSource)/lambdaSource,2)),sourceFloor);
  } else {
    return 1e-10;
  }
}

// Common source temperature profile.
double sourceTemperature(double t, const double * GKYL_RESTRICT xn, void *ctx)
{
  double x = xn[0], z = xn[2];
  struct gk_app_ctx *app = ctx;
  double xSource = app->xSource;
  double lambdaSource = app->lambdaSource;
  if (x < xSource + 3 * lambdaSource) {
    return 80*GKYL_ELEMENTARY_CHARGE;
  } else {
    return 30*GKYL_ELEMENTARY_CHARGE;
  }
}

// Initial density.
double densityInit(double t, const double * GKYL_RESTRICT xn, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  struct gk_app_ctx *app = ctx;
  double Ls = app->Lz/4;
  double xSource[3] = {x, y, 0};
  double effectiveSource = sourceDensity(t, xSource, ctx);
  double c_ss = sqrt(5/3*sourceTemperature(t, xSource, ctx)/app->mi);
  double nPeak = 4*sqrt(5)/3/c_ss*Ls*effectiveSource/2;
  pcg64_random_t rng = gkyl_pcg64_init(0);
  double perturb = 1e-3*(gkyl_pcg64_rand_double(&rng) - 0.5)*2.0;
  if (fabs(z) <= Ls) {
    return nPeak * (1 + sqrt(1-pow(z/Ls,2)))/2 * (1+perturb);
  } else {
    return nPeak/2 * (1+perturb);
  }
}

// Initial temperature.
double temperatureInit(double t, const double * GKYL_RESTRICT xn, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  struct gk_app_ctx *app = ctx;
  double xSource = app->xSource;
  double lambdaSource = app->lambdaSource;
  if (x < xSource + 3*lambdaSource) {
    return 80*GKYL_ELEMENTARY_CHARGE;
  } else {
    return 20*GKYL_ELEMENTARY_CHARGE;
  }
}

// Initial ion drift speed.
double driftSpeed(const double * GKYL_RESTRICT xn, void *ctx){
  double x = xn[0], y = xn[1], z = xn[2];
  struct gk_app_ctx *app = ctx;
  double xSource = app->xSource;
  double lambdaSource = app->lambdaSource;
  double Lz = app->Lz;
  double mi = app->mi;
  double Te;
  if (x < xSource + 3*lambdaSource) {
    Te = 50*GKYL_ELEMENTARY_CHARGE;
  } else {
    Te = 20*GKYL_ELEMENTARY_CHARGE;
  }
  double Ls = Lz/4;
  if (fabs(z) <= Ls) {
    return z/Ls*sqrt(Te/mi);
  } else {
    return (z > 0 ? 1.0 : -1.0)*sqrt(Te/mi);
  }
}

// Mapping
double Rx(const double *xc, void *ctx)
{
  return xc[0];
}
double Zx(const double *xc, void *ctx)
{
  return xc[2];
}

// Magnetic field functions.
double Bphi(const double *xc, void *ctx)
{
  double x = xc[0];
  struct gk_app_ctx *app = ctx;
  double B_axis = app->B_axis;
  double R_axis = app->R_axis;
  double R = Rx(xc, ctx);
  return B_axis*R_axis/R;
}

double Bvert(const double *xc, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double Bvx0 = app->Bvx0;
  int n = app->n;
  double R0 = app->R0;
  double R = Rx(xc, ctx);
  return Bvx0 * pow(R/R0, n);
}

double Bmag(const double *xc, void *ctx)
{
  double Bt = Bphi(xc, ctx);
  double Bv = Bvert(xc, ctx);
  return Bv * sqrt(1 + pow(Bt/Bv,2));
}

// Toroidal angle coordinate.
double phix(const double *xc, void *ctx)
{
  double x = xc[0], y = xc[1], z = xc[2];
  struct gk_app_ctx *app = ctx;
  double R_axis = app->R_axis;
  double Lz = app->Lz;
  double Bt = Bphi(xc, ctx);
  double Bv = Bvert(xc, ctx);
  return y/R_axis + (Bt * z)/(Bv * x);
}

double qprofile(const double *xc, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double H = app->H;
  double R = Rx(xc, ctx);
  double Bt = Bphi(xc, ctx);
  double Bv = Bvert(xc, ctx);

  return (H * Bt)/(2*GKYL_PI*R*Bv);
}

// Interface function calls.
void source_density(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{ fout[0] = sourceDensity(t, xn, ctx); }
void source_temperature(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{ fout[0] = sourceTemperature(t, xn, ctx); }
void density_init(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{ fout[0] = densityInit(t, xn, ctx); }
void temp_init(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{ fout[0] = temperatureInit(t, xn, ctx); }
void upar_ion_init(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{ fout[0] = driftSpeed(xn, ctx); }
void bfield_func(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{ 
  // Cartesian formulation of the magnetic field.
  double Bt = Bphi(xc, ctx);
  double Bv = Bvert(xc, ctx);
  double phi = phix(xc, ctx);
  fout[0] = -Bt*sin(phi);
  fout[1] = Bt*cos(phi);
  fout[2] = Bv;
}
void zero_func(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{ fout[0] = 0.0; }


// Geometry evaluation functions for the gk app
void mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  // Map to cylindrical (R, Z, phi) coordinates.
  double R   = Rx(xc, ctx);
  double Z   = Zx(xc, ctx);
  double phi = phix(xc, ctx);
  // Map to Cartesian (X, Y, Z) coordinates.
  double X = R*cos(phi);
  double Y = R*sin(phi);

  xp[0] = X; xp[1] = Y; xp[2] = Z;
}

struct gk_app_ctx
create_ctx(void)
{
  int cdim = 3, vdim = 2; // Dimensionality.

  // Universal constant parameters.
  double eps0 = GKYL_EPSILON0, eV = GKYL_ELEMENTARY_CHARGE;
  double mp = GKYL_PROTON_MASS, me = GKYL_ELECTRON_MASS;
  double qi = eV; // ion charge
  double qe = -eV; // electron charge

  // Plasma parameters. Chosen based on the value of a cubic sline
  // between the last TS data inside the LCFS and the probe data in
  // in the far SOL, near R=0.475 m.
  double AMU = 2.01410177811;
  double mi = mp*AMU; // Deuterium ions.
  double Te0 = 40*eV;
  double Ti0 = 40*eV;
  double n0 = 7e18*10; // [1/m^3]

  // Geometry and magnetic field.
  double B_axis = 0.5;
  double R_axis = 0.85;
  double a0 = 0.5;
  double H = 2.4; // Poloidal length (used to get the field line pitch).
  double Lcx0 = 8.0; // Connection length in the middle of the domain.
  int shear = -2;

  double x0 = R_axis + a0;
  double R0 = x0;
  
  const double vec0[3] = {x0, 0.0, 0.0};
  double B0 = Bphi(vec0, &(struct gk_app_ctx){.B_axis=B_axis, .R_axis=R_axis});
  double Bvx0 = H * B0 / Lcx0;
  double n = shear + 2;

  // Source parameters.
  // double P_SOL = 0.62*1e6; // Power crossing the separatrix, in Watts.
  // double Psrc = P_SOL * Ly / (2*M_PI*Rc);
  double S0 = 5.7691e23*10; // Taken from the input file.
  double xSource = x0 - 0.05;
  double lambdaSource = 0.005;

  // Collisions;
  double nuFrac = 0.1;

  // Derived parameters.
  double vte = sqrt(Te0/me), vti = sqrt(Ti0/mi); // Thermal speeds.
  double c_s = sqrt(Te0/mi); // Sound speed.
  double omega_ci = fabs(qi*B_axis/mi); // Ion cyclotron frequency.
  double rho_s = c_s/omega_ci; // Ion sound gyroradius.

  // Box size.
  double Lx = 56*rho_s;
  double Ly = 100*rho_s*H/Lcx0;
  double Lz = H; // [m]

  double x_min = x0 - Lx/2;
  double x_max = x0 + Lx/2;
  double y_min = -Ly/2;
  double y_max = Ly/2;
  double z_min = -Lz/2;
  double z_max = Lz/2;

  // Grid parameters
  int Nx = 8;
  int Ny = 4;
  int Nz = 4;
  int Nvpar = 4;
  int Nmu = 2;
  int poly_order = 1;

  double vpar_max_elc = 4.*vte;
  double mu_max_elc = 12*me*pow(vte,2)/(2*B_axis);
  double vpar_max_ion = 4.*vti;
  double mu_max_ion = 12*mi*pow(vti,2)/(2*B_axis);

  double t_end = 1.e-6; // End time, should terminate in 43 steps.
  int num_frames = 1;
  double write_phase_freq = 1.0; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gk_app_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,

    .B_axis = B_axis,
    .Bvx0 = Bvx0,
    .R_axis = R_axis,
    .R0 = R0,
    .a0 = a0,
    .H = H,
    .x0 = x0,
    .n = n,

    .Lx = Lx,
    .Ly = Ly,
    .Lz = Lz,
    .x_min = x_min, .x_max = x_max,
    .y_min = y_min, .y_max = y_max,
    .z_min = z_min, .z_max = z_max,

    .me = me, .qe = qe,
    .mi = mi, .qi = qi,
    .n0 = n0, .Te0 = Te0, .Ti0 = Ti0,
  
    .nu_frac = nuFrac,
    
    .S0 = S0,
    .xSource = xSource,
    .lambdaSource = lambdaSource,
  
    .Nx = Nx,
    .Ny = Ny,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nx, Ny, Nz, Nvpar, Nmu},
    .poly_order = poly_order,
    .vpar_max_elc = vpar_max_elc, .mu_max_elc = mu_max_elc,
    .vpar_max_ion = vpar_max_ion, .mu_max_ion = mu_max_ion,

    .t_end = t_end, .num_frames = num_frames,
    .write_phase_freq = write_phase_freq,
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };
  return ctx;
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

  struct gk_app_ctx ctx = create_ctx(); // Context for init functions.

  // Extract variables from command line arguments.
  sscanf(app_args.opt_args, "max_run_time=%lf",&ctx.max_run_time);

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  // electrons
  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.qe, .mass = ctx.me,
    .vdim = ctx.vdim,
    .lower = {-ctx.vpar_max_elc, 0.0},
    .upper = { ctx.vpar_max_elc, ctx.mu_max_elc},
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .ctx_upar = &ctx,
      .ctx_temp = &ctx,
      .density = density_init,
      .upar = zero_func,
      .temp = temp_init,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_EM,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .num_cross_collisions = 1,
      .collide_with = { "elc" },
      .den_ref = ctx.n0,
      .temp_ref = ctx.Te0,
      .nu_frac = ctx.nu_frac,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
        .ctx_density = &ctx,
        .ctx_upar = &ctx,
        .ctx_temp = &ctx,
        .density = source_density,
        .upar = zero_func,
        .temp = source_temperature,
      },
      .diagnostics = {
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
      },
    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX },
      { .dir = 2, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH, },
      { .dir = 2, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH, },
    },

    .num_diag_moments = 9,
    .diag_moments = {GKYL_F_MOMENT_HAMILTONIAN, GKYL_F_MOMENT_BIMAXWELLIAN, 
      GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, 
      GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP},

    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },

    .boundary_flux_diagnostics = {
      .num_diag_moments = 1,
      .diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    },
    
    .time_rate_diagnostics = true,
  };

  // ions
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.qi, .mass = ctx.mi,
    .vdim = ctx.vdim,
    .lower = { -ctx.vpar_max_ion, 0.0},
    .upper = { ctx.vpar_max_ion, ctx.mu_max_ion},
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .ctx_upar = &ctx,
      .ctx_temp = &ctx,
      .density = density_init,
      .upar = upar_ion_init,
      .temp = temp_init,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_EM,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .num_cross_collisions = 1,
      .collide_with = { "ion" },
      .den_ref = ctx.n0,
      .temp_ref = ctx.Ti0, 
      .nu_frac = ctx.nu_frac,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
        .ctx_density = &ctx,
        .ctx_upar = &ctx,
        .ctx_temp = &ctx,
        .density = source_density,
        .upar = zero_func,
        .temp = source_temperature,
      },
      .diagnostics = {
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
      },
    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX },
      { .dir = 2, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH },
      { .dir = 2, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH },
    },

    .num_diag_moments = 9,
    .diag_moments = {GKYL_F_MOMENT_HAMILTONIAN, GKYL_F_MOMENT_BIMAXWELLIAN, 
      GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, 
      GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP},

    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },

    .boundary_flux_diagnostics = {
      .num_diag_moments = 1,
      .diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    },

    .time_rate_diagnostics = true,
  };

  // field
  struct gkyl_gyrokinetic_field field = {
    .poisson_bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    },
    .ampere_bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    },
    .mu0 = GKYL_MU0,
    .time_rate_diagnostics = true,
  };

  // GK app
  struct gkyl_gk app_inp = {

    .cfl_frac_omegaH = 1.0,
    .cfl_frac = 0.9,

    .cdim = ctx.cdim,
    .lower = { ctx.x_min, ctx.y_min, ctx.z_min },
    .upper = { ctx.x_max, ctx.y_max, ctx.z_max },
    .cells = { cells_x[0], cells_x[1], cells_x[2] },
    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,

    .geometry = {
      .geometry_id = GKYL_GEOMETRY_MAPC2P,
      .mapc2p = mapc2p, // mapping of computational to physical space
      .c2p_ctx = &ctx,
      .bfield_func = bfield_func, // magnetic field
      .bfield_ctx = &ctx,
    },

    .num_periodic_dir = 1,
    .periodic_dirs = { 1 },

    .num_species = 2,
    .species = { elc, ion },
    .field = field,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .comm = comm,
      .cuts = { app_args.cuts[0], app_args.cuts[1], app_args.cuts[2] },
    }
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
    // .print_verbosity = {
    //   .enabled = true,
    //   .frequency = 1.0,
    // }
  };

  gkyl_gyrokinetic_run_simulation(&run_inp);
  
  gkyl_gyrokinetic_comms_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif

  return 0;
}