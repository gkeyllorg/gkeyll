#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_efit.h>
#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_gyrokinetic_run.h>
#include <gkyl_tok_geo.h>

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
  double n0; // Reference number density (1 / m^3).

  double B_axis; // Magnetic field axis (simple toroidal coordinates).
  double R0; // Major radius (simple toroidal coordinates).
  double a0; // Minor axis (simple toroidal coordinates).

  double nu_frac; // Collision frequency fraction.

  double k_perp_rho_s; // Product of perpendicular wavenumber and ion-sound gyroradius.

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

  double k_perp; // Perpendicular wavenumber (for Poisson solver).

  double n_src; // Source number density.
  double T_src; // Source temperature.

  double c_s_src; // Source sound speed.
  double n_peak; // Peak number density.

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

struct sheath_ctx
create_ctx(void)
{
  int cdim = 1, vdim = 2; // Dimensionality.

  // Physical constants (using non-normalized physical units).
  double epsilon0 = GKYL_EPSILON0; // Permittivity of free space.
  double mass_elc = GKYL_ELECTRON_MASS; // Electron mass.
  double mass_ion = 2.014 * GKYL_PROTON_MASS; // Proton mass.
  double charge_elc = -GKYL_ELEMENTARY_CHARGE; // Electron charge.
  double charge_ion = GKYL_ELEMENTARY_CHARGE; // Proton charge.

  double Te = 40.0 * GKYL_ELEMENTARY_CHARGE; // Electron temperature.
  double Ti = 40.0 * GKYL_ELEMENTARY_CHARGE; // Ion temperature.
  double n0 = 7.0e18; //  Reference number density (1 / m^3).

  double B_axis = 0.5; // Magnetic field axis (simple toroidal coordinates).
  double R0 = 0.85; // Major radius (simple toroidal coordinates).
  double a0 = 0.15; // Minor axis (simple toroidal coordinates).

  double nu_frac = 0.1; // Collision frequency fraction.

  double k_perp_rho_s = 0.3; // Product of perpendicular wavenumber and ion-sound gyroradius.

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
  double omega_ci = fabs(charge_ion * B0 / mass_ion); // Ion cyclotron frequency.
  double rho_s = c_s / omega_ci; // Ion-sound gyroradius.

  double k_perp = k_perp_rho_s / rho_s; // Perpendicular wavenumber (for Poisson solver).

  double n_src = 2.870523e21; // Source number density.
  double T_src = 2.0 * Te; // Source temperature.

  double c_s_src = sqrt((5.0 / 3.0) * T_src / mass_ion); // Source sound speed.
  double n_peak = 4.0 * sqrt(5.0) / 3.0 / c_s_src * 0.5 * n_src; // Peak number density.

  // Simulation parameters.
  int Nz = 8; // Cell count (configuration space: z-direction).
  int Nvpar = 6; // Cell count (velocity space: parallel velocity direction).
  int Nmu = 4; // Cell count (velocity space: magnetic moment direction).
  double Lz = 4.0; // Domain size (configuration space: z-direction).
  double vpar_max_elc = 4.0 * vte; // Domain boundary (electron velocity space: parallel velocity direction).
  double mu_max_elc = (3.0 / 2.0) * 0.5 * mass_elc * pow(4.0 * vte, 2.0) / (2.0 * B0); // Domain boundary (electron velocity space: magnetic moment direction).
  double vpar_max_ion = 4.0 * vti; // Domain boundary (ion velocity space: parallel velocity direction).
  double mu_max_ion = (3.0 / 2.0) * 0.5 * mass_ion * pow(4.0 * vti, 2.0) / (2.0 * B0); // Domain boundary (ion velocity space: magnetic moment direction).
  int poly_order = 1; // Polynomial order.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 6.0e-6; // Final simulation time.
  int num_frames = 1; // Number of output frames.
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
    .n0 = n0,
    .B_axis = B_axis,
    .R0 = R0,
    .a0 = a0,
    .nu_frac = nu_frac,
    .k_perp_rho_s = k_perp_rho_s,
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
    .k_perp = k_perp,
    .n_src = n_src,
    .T_src = T_src,
    .c_s_src = c_s_src,
    .n_peak = n_peak,
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
  struct sheath_ctx *app = ctx;
  double z = xn[0];

  double n_peak = app->n_peak;
  double Lz = app->Lz;

  double n = 0.0;

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
evalElcSourceDensityInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sheath_ctx *app = ctx;
  double z = xn[0];

  double n_src = app->n_src;
  double Lz = app->Lz;

  double n = 0.0;

  if (fabs(z) < 0.25 * Lz) {
    n = n_src; // Electron source total number density (left).
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

  double T_src = app->T_src;

  // Set electron source isotropic temperature.
  fout[0] = T_src;
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
  double z = xn[0];

  double n_peak = app->n_peak;
  double Lz = app->Lz;

  double n = 0.0;

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
evalIonSourceDensityInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sheath_ctx *app = ctx;
  double z = xn[0];

  double n_src = app->n_src;
  double Lz = app->Lz;

  double n = 0.0;

  if (fabs(z) < 0.25 * Lz) {
    n = n_src; // Ion source total number density (left).
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

  double T_src = app->T_src;

  // Set ion source isotropic temperature.
  fout[0] = T_src;
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
  // Set physical coordinates (X, Y, Z) from computational coordinates (x, y, z).
  xp[0] = zc[0]; xp[1] = zc[1]; xp[2] = zc[2];
}

void
bfield_func(double t, const double* GKYL_RESTRICT zc, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sheath_ctx *app = ctx;

  double B0 = app->B0;

  // zc are computational coords. 
  // Set Cartesian components of magnetic field.
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[0] = B0;
}

struct gkyl_gk_block_geom*
create_gk_block_geom(void *ctx)
{

  struct sheath_ctx *app = ctx;

  struct gkyl_gk_block_geom *bgeom = gkyl_gk_block_geom_new(1, 3);

  /* Block layout and coordinates

   x  
   ^  
   |
   1  +------------------+------------------+------------------+
   |  |b0                |b1                |b2                |
   |  |lower SOL         |middle            |upper             |
   |  |                  |                  |                  |
   0  +------------------+------------------+------------------+
      0 -----------------1------------------2-------------------3> z

      Edges that touch coincide are physically connected unless
      otherwise indicated by a special symbol. Edges with a special
      symbol such as o,x,%, or % are instead connected to the other
      edge with the same symbol. Edges that do not coincide with
      another edge are a physical boundary.
  */  


  int nz=app->Nz;
  double Lz = app->Lz;

  // block 0. Lower SOL.
  gkyl_gk_block_geom_set_block(bgeom, 0, &(struct gkyl_gk_block_geom_info) {
      .lower = { -Lz/2.0 },
      .upper = { -Lz/4.0 },
      .cells = { nz/4},
      .cuts = { 1 },
      .geometry = {
        .geometry_id = GKYL_GEOMETRY_MAPC2P,
        .world = { 0.0, 0.0 },
        .mapc2p = mapc2p,
        .c2p_ctx = app,
        .bfield_func = bfield_func,
        .bfield_ctx = app 
      },

      .connections[0] = { // z-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 1, .dir = 0, .edge = GKYL_LOWER_POSITIVE},
      }
    }
  );

  // block 1. Middle SOL.
  gkyl_gk_block_geom_set_block(bgeom, 1, &(struct gkyl_gk_block_geom_info) {
      .lower = { -Lz/4.0},
      .upper = { Lz/4.0},
      .cells = { nz/2},
      .cuts = { 1 },
      .geometry = {
        .geometry_id = GKYL_GEOMETRY_MAPC2P,
        .world = { 0.0, 0.0 },
        .mapc2p = mapc2p,
        .c2p_ctx = app,
        .bfield_func = bfield_func,
        .bfield_ctx = app
      },

      .connections[0] = { // z-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 2, .dir = 0, .edge = GKYL_LOWER_POSITIVE},
      }

    }
  );

  // block 2. Upper SOL.
  gkyl_gk_block_geom_set_block(bgeom, 2, &(struct gkyl_gk_block_geom_info) {
      .lower = { Lz/4.0},
      .upper = { Lz/2.0},
      .cells = { nz/4},
      .cuts = { 1 },
      .geometry = {
        .geometry_id = GKYL_GEOMETRY_MAPC2P,
        .world = { 0.0, 0.0 },
        .mapc2p = mapc2p,
        .c2p_ctx = app,
        .bfield_func = bfield_func,
        .bfield_ctx = app
      },

      .connections[0] = { // z-direction connections
        { .bid = 1, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
      }
    }
  );

  return bgeom;
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

  // Construct block geometry.
  struct gkyl_gk_block_geom *bgeom = create_gk_block_geom(&ctx);
  int nblocks = gkyl_gk_block_geom_num_blocks(bgeom);

  // Elc Species
  // all data is common across blocks
  struct gkyl_gyrokinetic_multib_species_pb elc_blocks[1];
  elc_blocks[0] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 0,
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
        .num_diag_moments = 5,
        .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP },
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
//        .time_integrated = true,
      }
    },


  };


  struct gkyl_gyrokinetic_bc elc_phys_bcs[] = {
    // block 1 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    { .bidx = 2, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH },
  };

  struct gkyl_gyrokinetic_multib_species elc = {
    .name = "elc",
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .vdim = ctx.vdim,
    .lower = { -ctx.vpar_max_elc, 0.0},
    .upper = {  ctx.vpar_max_elc, ctx.mu_max_elc}, 
    .cells = { cells_v[0], cells_v[1] },

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

    .num_diag_moments = 7,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP },
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    .time_rate_diagnostics = true,
    .boundary_flux_diagnostics = {
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
//      .time_integrated = true,
    },

    .duplicate_across_blocks = true,
    .blocks = elc_blocks,
    .num_physical_bcs = 2,
    .bcs = elc_phys_bcs,
  };


  // Ion Species
  // all data is common across blocks
  struct gkyl_gyrokinetic_multib_species_pb ion_blocks[1];
  ion_blocks[0] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 0,

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
        .num_diag_moments = 5,
        .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP },
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
//        .time_integrated = true,
      }
    },

  };


  struct gkyl_gyrokinetic_bc ion_phys_bcs[] = {
    { .bidx = 0, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    { .bidx = 2, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH },
  };

  struct gkyl_gyrokinetic_multib_species ion = {
    .name = "ion",
    .charge = ctx.charge_ion, .mass = ctx.mass_ion,
    .vdim = ctx.vdim,
    .lower = { -ctx.vpar_max_ion, 0.0},
    .upper = {  ctx.vpar_max_ion, ctx.mu_max_ion}, 
    .cells = { cells_v[0], cells_v[1] },

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

    .num_diag_moments = 7,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP },
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    .time_rate_diagnostics = true,
    .boundary_flux_diagnostics = {
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
//      .time_integrated = true,
    },

    .duplicate_across_blocks = true,
    .blocks = ion_blocks,
    .num_physical_bcs = 2,
    .bcs = ion_phys_bcs,
  };

 

  // Field object
  struct gkyl_gyrokinetic_multib_field_pb field_blocks[1];
  field_blocks[0] = (struct gkyl_gyrokinetic_multib_field_pb) {
    .polarization_bmag = ctx.B0,
    .kperpSq = ctx.k_perp * ctx.k_perp,
    .time_rate_diagnostics = true,
  };

  struct gkyl_gyrokinetic_bc field_phys_bcs[] = {
  };

  struct gkyl_gyrokinetic_multib_field field = {
    .duplicate_across_blocks = true,
    .blocks = field_blocks, 
    .num_physical_bcs = 0, 
    .bcs = field_phys_bcs,
  };

  struct gkyl_gyrokinetic_multib *app_inp = gkyl_calloc(1,sizeof(struct gkyl_gyrokinetic_multib));
  app_inp->cdim = ctx.cdim;
  app_inp->poly_order = 1;
  app_inp->basis_type = app_args.basis_type;
  app_inp->use_gpu = app_args.use_gpu;
  app_inp->cfl_frac = 1.0;
  app_inp->gk_block_geom = bgeom;
  app_inp->num_species = 2;
  app_inp->species[0] = elc;
  app_inp->species[1] = ion;
  app_inp->field = field;
  app_inp->comm = comm;


  // Set app output name from the executable name (argv[0]).
  snprintf(app_inp->name, sizeof(app_inp->name), "%s", app_args.app_name);
  struct gkyl_gyrokinetic_run_inp run_inp = {
    .app_type = GKYL_GK_MULTIB,
    .multib_app_inp = *app_inp,
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

  gkyl_gk_block_geom_release(bgeom);
  gkyl_gyrokinetic_comms_release(comm);
  gkyl_free(app_inp);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif

  return 0;
}
