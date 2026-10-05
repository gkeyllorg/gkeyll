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
  double B_axis;
  double R0;
  double Rc;
  double Bv0;

  // Plasma parameters.
  double me;
  double qe;
  double mi;
  double qi;
  double n0;
  double Te0;
  double Ti0;
  double Ntot;
  double center_ic[3], sigma_ic[3];

  // Collisions.
  double nu_frac;

  // Adaptive source parameters.
  double P_SOL;
  int num_sources;
  double center_srcCORE[3], sigma_srcCORE[3];
  double center_srcRECYlo[3], sigma_srcRECYlo[3];
  double center_srcRECYup[3], sigma_srcRECYup[3];
  double floor_src;

  // Grid parameters.
  double Lx; // Domain size in radial direction.
  double Ly; // Domain size in binormal direction.
  double Lz; // Domain size along magnetic field.
  double x_min;
  double x_max;
  double y_min;
  double y_max;
  double z_min;
  double z_max;
  int Nx;
  int Ny;
  int Nz;
  int Nvpar;
  int Nmu;
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  int poly_order;
  double vpar_max_elc;
  double mu_max_elc;
  double vpar_max_ion;
  double mu_max_ion;

  double t_end; // End time.
  int num_frames; // Number of output frames.
  double
    write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

// Radial coordinate mapping.
double
Rx(const double *xc, void *ctx)
{
  return xc[0];
}

// Toroidal magnetic field.
double
Bphi(const double *xc, void *ctx)
{
  double x = xc[0];
  struct gk_app_ctx *app = ctx;
  double B_axis = app->B_axis;
  double R0 = app->R0;
  double R = Rx(xc, ctx);
  return B_axis * R0 / R;
}

// Vertical magnetic field.
double
Bvert(const double *xc, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double Bv0 = app->Bv0;
  return Bv0;
}

// Magnetic field magnitude.
double
bmag(const double *xc, void *ctx)
{
  double Bt = Bphi(xc, ctx);
  double Bv = Bvert(xc, ctx);
  return sqrt(pow(Bt, 2) + pow(Bv, 2));
}

// Field line pitch.
double
thetax(const double *xc, void *ctx)
{
  return asin(Bvert(xc, ctx) / bmag(xc, ctx));
}

// Parallel coordinate mapping.
double
Zx(const double *xc, void *ctx)
{
  double theta = thetax(xc, ctx);
  return xc[2] * sin(theta);
}

double
phix(const double *xc, void *ctx)
{
  double x = xc[0], y = xc[1], z = xc[2];
  struct gk_app_ctx *app = ctx;
  double Rc = app->Rc;
  double Bt = Bphi(xc, ctx);
  double Bv = Bvert(xc, ctx);

  // Helical sheared mapping. Passes right hand check and conserves particle.
  double theta = thetax(xc, ctx);
  return y / Rc + (Bt * z * sin(theta)) / (Bv * x);
}

void
bfield_func(double t, const double *xc, double *GKYL_RESTRICT fout, void *ctx)
{
  // Cartesian formulation of the magnetic field.
  double Bt = Bphi(xc, ctx);
  double Bv = Bvert(xc, ctx);
  double phi = phix(xc, ctx);
  fout[0] = -Bt * sin(phi);
  fout[1] = Bt * cos(phi);
  fout[2] = Bv;
}

void
zero_func(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

// Geometry evaluation functions for the gk app
void
mapc2p(double t, const double *xc, double *GKYL_RESTRICT xp, void *ctx)
{
  double x = xc[0], y = xc[1], z = xc[2];

  // Map to cylindrical (R, Z, phi) coordinates.
  double R = Rx(xc, ctx);
  double Z = Zx(xc, ctx);
  double phi = phix(xc, ctx);
  // Map to Cartesian (X, Y, Z) coordinates.
  double X = R * cos(phi);
  double Y = R * sin(phi);

  xp[0] = X;
  xp[1] = Y;
  xp[2] = Z;
}

// Taken from rt gk d3d 3x2c, is this the non uniform v grid mapping?
void
mapc2p_vel_elc(double t, const double *vc, double *GKYL_RESTRICT vp, void *ctx)
{
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

void
mapc2p_vel_ion(double t, const double *vc, double *GKYL_RESTRICT vp, void *ctx)
{
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
  double mi = mp * AMU; // Deuterium ions.
  double Te0 = 40 * eV;
  double Ti0 = 40 * eV;
  double n0 = 7e19; // [1/m^3] (10 times NSTX ref.)
  double Ntot = 6.0e18; // Total number of particles in the simulation box.

  // Geometry and magnetic field.
  double B_axis = 0.5;
  double R0 = 0.85;
  double a0 = 0.5;
  double Lp = 2.4; // Poloidal length at x0.
  double Lt = 8.0; // Toroidal length at x0.
  double Rc = R0 + a0;
  double B0 = B_axis * R0 / Rc;

  double sintheta = Lp / Lt;
  double Bv0 = B0 * sintheta;

  // Collisions;
  double nuFrac = 0.1;

  // Derived parameters.
  double vte = sqrt(Te0 / me), vti = sqrt(Ti0 / mi); // Thermal speeds.
  double c_s = sqrt(Te0 / mi); // Sound speed.
  double omega_ci = fabs(qi * B0 / mi); // Ion cyclotron frequency.
  double rho_s = c_s / omega_ci; // Ion sound gyroradius.

  // Box size.
  double Lx = 50 * rho_s;
  double Ly = 50 * rho_s;
  double Lz = Lt; // [m]

  double x_min = Rc - Lx / 2;
  double x_max = Rc + Lx / 2;
  double y_min = -Ly / 2;
  double y_max = Ly / 2;
  double z_min = -Lz / 2;
  double z_max = Lz / 2;
  double center_ic[3] = {x_min, 0.0, 0.0}; // Initial condition centered at z=0.
  double sigma_ic[3] = {Lx, 0.0, Lz}; // Initial condition width along x and z.

  // 1+2 adaptive source setup:
  // - CORE: at the left boundary (x=xmin), injects a fixed power P_SOL
  //   and compensate for energy loss through inner radial boundary.
  // - RECY: neutral recycling model, put back particles lost through wall
  //   at 5 eV, no energy compensation. one for up, one for low.
  double P_SOL = 8.1e5; // Power injected in the SOL [W].
  int num_sources = 3;
  double center_srcCORE[3] = {x_min, 0.0, 0.0}; // Source centered at x=xmin, z=0.
  double sigma_srcCORE[3] = {Lx / 16, 0.0, Lz / 8}; // Source width along x and z.
  double center_srcRECYup[3] = {x_max, 0.0, Lz / 2}; // Source centered at x=xmax, z=Lz/2 (periodic).
  double sigma_srcRECYup[3] = {Lx, 0.0, Lz / 16}; // Source width along x and z.
  double center_srcRECYlo[3] = {
    x_max, 0.0, -Lz / 2
  }; // Source centered at x=xmax, z=-Lz/2 (periodic).
  double sigma_srcRECYlo[3] = {Lx, 0.0, Lz / 16}; // Source width along x and z.
  double floor_src = 1e-10;

  // Grid parameters
  int Nx = 4;
  int Ny = 4;
  int Nz = 4;
  int Nvpar = 4;
  int Nmu = 2;
  int poly_order = 1;

  double vpar_max_elc = 4. * vte;
  double mu_max_elc = 12 * me * pow(vte, 2) / (2 * B0);
  double vpar_max_ion = 4. * vti;
  double mu_max_ion = 12 * mi * pow(vti, 2) / (2 * B0);

  double t_end = 0.1e-6; // End time, should terminate in 43 steps.
  int num_frames = 1;
  double write_phase_freq =
    1.0; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num = num_frames * 100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gk_app_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .B_axis = B_axis,
    .R0 = R0,
    .Rc = Rc,
    .Bv0 = Bv0,
    .Lx = Lx,
    .Ly = Ly,
    .Lz = Lz,
    .x_min = x_min,
    .x_max = x_max,
    .y_min = y_min,
    .y_max = y_max,
    .z_min = z_min,
    .z_max = z_max,

    .me = me,
    .qe = qe,
    .mi = mi,
    .qi = qi,
    .n0 = n0,
    .Te0 = Te0,
    .Ti0 = Ti0,
    .Ntot = Ntot,
    .center_ic = {center_ic[0], center_ic[1], center_ic[2]},
    .sigma_ic = {sigma_ic[0], sigma_ic[1], sigma_ic[2]},

    .nu_frac = nuFrac,

    .P_SOL = P_SOL,
    .num_sources = num_sources,
    .center_srcCORE = {center_srcCORE[0], center_srcCORE[1], center_srcCORE[2]},
    .sigma_srcCORE = {sigma_srcCORE[0], sigma_srcCORE[1], sigma_srcCORE[2]},
    .center_srcRECYlo = {center_srcRECYlo[0], center_srcRECYlo[1], center_srcRECYlo[2]},
    .sigma_srcRECYlo = {sigma_srcRECYlo[0], sigma_srcRECYlo[1], sigma_srcRECYlo[2]},
    .center_srcRECYup = {center_srcRECYup[0], center_srcRECYup[1], center_srcRECYup[2]},
    .sigma_srcRECYup = {sigma_srcRECYup[0], sigma_srcRECYup[1], sigma_srcRECYup[2]},
    .floor_src = floor_src,

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

    .t_end = t_end,
    .num_frames = num_frames,
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
  if (app_args.use_mpi) {
    MPI_Init(&argc, &argv);
  }
#endif

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct gk_app_ctx ctx = create_ctx(); // Context for init functions.

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d = 0; d < ctx.cdim; d++) {
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  }
  for (int d = 0; d < ctx.vdim; d++) {
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim + d]);
  }

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  // Projection parameters for the sources.
  struct gkyl_gyrokinetic_projection proj_srcCORE_e = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_GAUSSIAN,
    .gaussian_mean = {ctx.center_srcCORE[0], ctx.center_srcCORE[1], ctx.center_srcCORE[2]},
    .gaussian_std_dev = {ctx.sigma_srcCORE[0], ctx.sigma_srcCORE[1], ctx.sigma_srcCORE[2]},
    .total_num_particles = 0.0,
    .total_kin_energy = ctx.P_SOL / 2,
    .temp_max = 5.0 * ctx.Te0,
    .temp_min = 0.1 * ctx.Te0,
    .f_floor = ctx.floor_src,
  };
  struct gkyl_gyrokinetic_adapt_source adapt_srcCORE_e = {
    .adapt_to_species = "elc",
    .adapt_particle = true,
    .adapt_energy = true,
    .num_boundaries = 1,
    .dir = {0},
    .edge = {GKYL_LOWER_EDGE},
  };

  struct gkyl_gyrokinetic_projection proj_srcRECYlo_e = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_GAUSSIAN,
    .gaussian_mean = {ctx.center_srcRECYlo[0], ctx.center_srcRECYlo[1], ctx.center_srcRECYlo[2]},
    .gaussian_std_dev = {ctx.sigma_srcRECYlo[0], ctx.sigma_srcRECYlo[1], ctx.sigma_srcRECYlo[2]},
    .total_num_particles = 0.0,
    .total_kin_energy = 0.0,
    .temp_max = 5.0 * ctx.Te0,
    .temp_min = 5.0 * 1.602e-19, // 5 eV
    .f_floor = ctx.floor_src,
  };
  struct gkyl_gyrokinetic_adapt_source adapt_srcRECYlo_e = {
    .adapt_to_species = "ion", // ambipolar
    .adapt_particle = true,
    .adapt_energy = false,
    .num_boundaries = 1,
    .dir = {2},
    .edge = {GKYL_LOWER_EDGE},
  };

  struct gkyl_gyrokinetic_projection proj_srcRECYup_e = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_GAUSSIAN,
    .gaussian_mean = {ctx.center_srcRECYup[0], ctx.center_srcRECYup[1], ctx.center_srcRECYup[2]},
    .gaussian_std_dev = {ctx.sigma_srcRECYup[0], ctx.sigma_srcRECYup[1], ctx.sigma_srcRECYup[2]},
    .total_num_particles = 0.0,
    .total_kin_energy = 0.0,
    .temp_max = 5.0 * ctx.Te0,
    .temp_min = 5.0 * 1.602e-19, // 5 eV
    .f_floor = ctx.floor_src,
  };
  struct gkyl_gyrokinetic_adapt_source adapt_srcRECYup_e = {
    .adapt_to_species = "ion", // ambipolar
    .adapt_particle = true,
    .adapt_energy = false,
    .num_boundaries = 1,
    .dir = {2},
    .edge = {GKYL_UPPER_EDGE},
  };

  // electrons
  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.qe,
    .mass = ctx.me,
    .vdim = ctx.vdim,
    .cells = {cells_v[0], cells_v[1]},
    .polarization_density = ctx.n0,

    .lower = {-1.0 / sqrt(2.0), 0.0},
    .upper = {1.0 / sqrt(2.0), 1.0},
    .mapc2p = {.mapping = mapc2p_vel_elc, .ctx = &ctx},

    .projection =
      {
        .proj_id = GKYL_PROJ_MAXWELLIAN_GAUSSIAN,
        .gaussian_mean = {ctx.center_ic[0], ctx.center_ic[1], ctx.center_ic[2]},
        .gaussian_std_dev = {ctx.sigma_ic[0], ctx.sigma_ic[1], ctx.sigma_ic[2]},
        .total_num_particles = ctx.Ntot,
        .total_kin_energy = ctx.Te0 * ctx.Ntot,
        .temp_max = 5.0 * ctx.Te0,
        .temp_min = 0.1 * ctx.Te0,
        .f_floor = 1e-10,
      },

    .collisionless = {.type = GKYL_GK_COLLISIONLESS_ES},

    .collisions =
      {
        .collision_id = GKYL_LBO_COLLISIONS,
        .num_cross_collisions = 1,
        .collide_with = {"elc"},
        .den_ref = ctx.n0,
        .temp_ref = ctx.Te0,
        .nu_frac = ctx.nu_frac,
      },

    .source =
      {
        .source_id = GKYL_PROJ_SOURCE,
        .num_sources = ctx.num_sources,
        .num_adapt_sources = ctx.num_sources,
        .projection[0] = proj_srcCORE_e,
        .adapt[0] = adapt_srcCORE_e,
        .projection[1] = proj_srcRECYlo_e,
        .adapt[1] = adapt_srcRECYlo_e,
        .projection[2] = proj_srcRECYup_e,
        .adapt[2] = adapt_srcRECYup_e,
        .diagnostics =
          {.num_integrated_diag_moments = 1, .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN}
          },
      },

    .bcs =
      {{.dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
       {.dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
       {
         .dir = 2,
         .edge = GKYL_LOWER_EDGE,
         .type = GKYL_BC_GK_SPECIES_SHEATH_SURROGATE,
         .aux_str = "gyrokinetic/data/nn_model/nn_model_sheath_bc_conv_MPE.kann",
       },
       {
         .dir = 2,
         .edge = GKYL_UPPER_EDGE,
         .type = GKYL_BC_GK_SPECIES_SHEATH_SURROGATE,
         .aux_str = "gyrokinetic/data/nn_model/nn_model_sheath_bc_conv_MPE.kann",
       }},

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

  // Projection parameters for the sources.
  struct gkyl_gyrokinetic_projection proj_srcCORE_i = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_GAUSSIAN,
    .gaussian_mean = {ctx.center_srcCORE[0], ctx.center_srcCORE[1], ctx.center_srcCORE[2]},
    .gaussian_std_dev = {ctx.sigma_srcCORE[0], ctx.sigma_srcCORE[1], ctx.sigma_srcCORE[2]},
    .total_num_particles = 0.0,
    .total_kin_energy = ctx.P_SOL / 2,
    .temp_max = 5.0 * ctx.Ti0,
    .temp_min = 0.1 * ctx.Ti0,
    .f_floor = ctx.floor_src,
  };
  struct gkyl_gyrokinetic_adapt_source adapt_srcCORE_i = {
    .adapt_to_species = "ion",
    .adapt_particle = true,
    .adapt_energy = true,
    .num_boundaries = 1,
    .dir = {0},
    .edge = {GKYL_LOWER_EDGE},
  };

  struct gkyl_gyrokinetic_projection proj_srcRECYlo_i = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_GAUSSIAN,
    .gaussian_mean = {ctx.center_srcRECYlo[0], ctx.center_srcRECYlo[1], ctx.center_srcRECYlo[2]},
    .gaussian_std_dev = {ctx.sigma_srcRECYlo[0], ctx.sigma_srcRECYlo[1], ctx.sigma_srcRECYlo[2]},
    .total_num_particles = 0.0,
    .total_kin_energy = 0.0,
    .temp_max = 5.0 * ctx.Ti0,
    .temp_min = 5.0 * 1.602e-19, // 5 eV
    .f_floor = ctx.floor_src,
  };
  struct gkyl_gyrokinetic_adapt_source adapt_srcRECYlo_i = {
    .adapt_to_species = "ion", // ambipolar
    .adapt_particle = true,
    .adapt_energy = false,
    .num_boundaries = 1,
    .dir = {2},
    .edge = {GKYL_LOWER_EDGE},
  };

  struct gkyl_gyrokinetic_projection proj_srcRECYup_i = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_GAUSSIAN,
    .gaussian_mean = {ctx.center_srcRECYup[0], ctx.center_srcRECYup[1], ctx.center_srcRECYup[2]},
    .gaussian_std_dev = {ctx.sigma_srcRECYup[0], ctx.sigma_srcRECYup[1], ctx.sigma_srcRECYup[2]},
    .total_num_particles = 0.0,
    .total_kin_energy = 0.0,
    .temp_max = 5.0 * ctx.Ti0,
    .temp_min = 5.0 * 1.602e-19, // 5 eV
    .f_floor = ctx.floor_src,
  };
  struct gkyl_gyrokinetic_adapt_source adapt_srcRECYup_i = {
    .adapt_to_species = "ion", // ambipolar
    .adapt_particle = true,
    .adapt_energy = false,
    .num_boundaries = 1,
    .dir = {2},
    .edge = {GKYL_UPPER_EDGE},
  };

  // ions
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.qi,
    .mass = ctx.mi,
    .vdim = ctx.vdim,
    .cells = {cells_v[0], cells_v[1]},
    .polarization_density = ctx.n0,

    .lower = {-1.0 / sqrt(2.0), 0.0},
    .upper = {1.0 / sqrt(2.0), 1.0},
    .mapc2p = {.mapping = mapc2p_vel_ion, .ctx = &ctx},

    .projection =
      {
        .proj_id = GKYL_PROJ_MAXWELLIAN_GAUSSIAN,
        .gaussian_mean = {ctx.center_ic[0], ctx.center_ic[1], ctx.center_ic[2]},
        .gaussian_std_dev = {ctx.sigma_ic[0], ctx.sigma_ic[1], ctx.sigma_ic[2]},
        .total_num_particles = ctx.Ntot,
        .total_kin_energy = ctx.Ti0 * ctx.Ntot,
        .temp_max = 5.0 * ctx.Ti0,
        .temp_min = 0.1 * ctx.Ti0,
        .f_floor = 1e-10,
      },

    .collisionless = {.type = GKYL_GK_COLLISIONLESS_ES},

    .collisions =
      {
        .collision_id = GKYL_LBO_COLLISIONS,
        .num_cross_collisions = 1,
        .collide_with = {"elc"},
        .den_ref = ctx.n0,
        .temp_ref = ctx.Te0,
        .nu_frac = ctx.nu_frac,
      },

    .source =
      {
        .source_id = GKYL_PROJ_SOURCE,
        .num_sources = ctx.num_sources,
        .num_adapt_sources = ctx.num_sources,
        .projection[0] = proj_srcCORE_i,
        .adapt[0] = adapt_srcCORE_i,
        .projection[1] = proj_srcRECYlo_i,
        .adapt[1] = adapt_srcRECYlo_i,
        .projection[2] = proj_srcRECYup_i,
        .adapt[2] = adapt_srcRECYup_i,
        .diagnostics =
          {.num_integrated_diag_moments = 1, .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN}
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

  // field
  struct gkyl_gyrokinetic_field field = {
    .poisson_bcs =
      {{.dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}},
       {.dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}}},
    .time_rate_diagnostics = true,
  };

  // GK app
  struct gkyl_gk app_inp = {

    .cfl_frac_omegaH = 1.0,
    .cfl_frac = 1.0,

    .cdim = ctx.cdim,
    .lower = {ctx.x_min, ctx.y_min, ctx.z_min},
    .upper = {ctx.x_max, ctx.y_max, ctx.z_max},
    .cells = {cells_x[0], cells_x[1], cells_x[2]},
    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,

    .geometry =
      {
        .geometry_id = GKYL_GEOMETRY_MAPC2P,
        .mapc2p = mapc2p, // mapping of computational to physical space
        .c2p_ctx = &ctx,
        .bfield_func = bfield_func, // magnetic field
        .bfield_ctx = &ctx,
      },

    .num_periodic_dir = 1,
    .periodic_dirs = {1},

    .num_species = 2,
    .species = {elc, ion},
    .field = field,

    .parallelism =
      {
        .use_gpu = app_args.use_gpu,
        .comm = comm,
        .cuts = {app_args.cuts[0], app_args.cuts[1], app_args.cuts[2]},
      },
  };

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
    .print_verbosity = {.enabled = true},
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