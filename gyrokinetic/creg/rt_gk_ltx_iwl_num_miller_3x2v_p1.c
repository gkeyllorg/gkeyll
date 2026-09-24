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
  double Lz;        // Domain size along magnetic field.
  double Ly;        // Domain size in binormal direction
  double z_min;  double z_max;
  double psi_min;  double psi_max;
  double psi_LCFS;    // Radial location of the last closed flux surface.
  double Lx_core; // Radial extent of core in psi
  double Lx; // Total radial extent in psi
  // Plasma parameters.
  double me;  double qe;
  double mi;  double qi;
  double n0;  double Te0;  double Ti0; 

  // Collisions.
  double nuFrac;
  double nuElc;  double nuIon;
  double nuElcIon;  double nuIonElc;

  // Source parameters.
  double n_src;  double Te_src;  double Ti_src;

  // Grid parameters.
  int Nz;
  int Ny;
  int Nvpar;
  int Nmu;
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  int poly_order;
  double vpar_max_elc;  double mu_max_elc;
  double vpar_max_ion;  double mu_max_ion;

  double t_end; // End time.
  int num_frames; // Number of output frames.
  double write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

void pfunc_upper(double s, double* RZ){
  RZ[0] = 0.14047;
  RZ[1] = -(s-0.061)*0.6;
}

void pfunc_lower(double s, double* RZ){
  RZ[0] = 0.14047;
  RZ[1] = (s-0.061)*0.6;
}

// Electron source profiles.
void density_src(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double Lz = app->Lz;
  double n_src = app->n_src;
  double x = xn[0];
  double z = xn[2];

  if ( app->psi_max - x < app->Lx_core/6.0 ) {
    fout[0] = app->n_src;
  } else {
    fout[0] = app->n_src*1.0e-5;
  }
}
void upar_elc_src(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}
void temp_elc_src(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double Te_src = app->Te_src;
  fout[0] = Te_src;
}

// Ion source profiles.
void upar_ion_src(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}
void temp_ion_src(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double Ti_src = app->Ti_src;
  fout[0] = Ti_src;
}

// Ion initial conditions
void density_init(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  fout[0] = app->n0;
}
void upar_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  fout[0] = 0.0;
}
void temp_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double Ti0 = app->Ti0;
  fout[0] = Ti0;
}

// Electron initial conditions
void upar_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  fout[0] = 0.0;
}
void temp_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double Te0 = app->Te0;
  fout[0] = Te0;
}

// Collision frequencies.
void evalNuElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  fout[0] = app->nuElc;
}
void evalNuIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  fout[0] = app->nuIon;
}
void evalNuElcIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  fout[0] = app->nuElcIon;
}
void evalNuIonElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  fout[0] = app->nuIonElc;
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

  // Geometry and magnetic field.
  double Lz = 2.0*(M_PI-1e-14);    // Domain size along magnetic field.
  double Ly = 2.0*M_PI/8.0;
  double B0 = 0.24;
  double psi_LCFS = 0.003172759514520552; // psi at LCFS. Taken from efit
  double psi_min = psi_LCFS - 0.0004; // inner flux surface of domain
  double psi_max = psi_LCFS + 0.0004; // outer flux surface of domain
  double Lx = psi_max - psi_min;
  double Lx_core = psi_LCFS - psi_min;

  // Plasma parameters. Chosen based on the value of a cubic sline
  // between the last TS data inside the LCFS and the probe data in
  // in the far SOL, near R=0.475 m.
  double mi  = mp;   // Hydrogen ions.
  double Te0 = 178*eV;
  double Ti0 = 70*eV;
  double n0  = 1.78e18;   // [1/m^3]

  double nuFrac = 1.0;
  // Electron-electron collision freq.
  double logLambdaElc = 6.6 - 0.5 * log(n0/1e20) + 1.5 * log(Ti0/eV);
  double nuElc = nuFrac * logLambdaElc * pow(eV, 4) * n0 /
    (6*sqrt(2.) * pow(M_PI,3./2.) * pow(eps0,2) * sqrt(me) * pow(Te0,3./2.));
  // Ion-ion collision freq.
  double logLambdaIon = 6.6 - 0.5 * log(n0/1e20) + 1.5 * log(Ti0/eV);
  double nuIon = nuFrac * logLambdaIon * pow(eV, 4) * n0 /
    (12 * pow(M_PI,3./2.) * pow(eps0,2) * sqrt(mi) * pow(Ti0,3./2.));
  double nuElcIon = nuElc*sqrt(2.0);
  double nuIonElc = nuElcIon*(me/mi);

  double vte = sqrt(Te0/me), vti = sqrt(Ti0/mi); // Thermal speeds.

  double c_s = sqrt(Te0/mi);
  double omega_ci = fabs(qi*B0/mi);
  double rho_s = c_s/omega_ci;

  // Source parameters
  double n_src = 1.95e22;
  double Te_src = 410*eV;
  double Ti_src = 40*eV;

  // Grid parameters
  int Nx = 4;
  int Nz = 32;
  int Ny = 2;
  int Nvpar = 6;
  int Nmu = 4;
  int poly_order = 1;

  double vpar_max_elc = 4.*vte;
  double mu_max_elc = me*pow(1.5*4*vte,2)/(2*B0);
  double vpar_max_ion = 4.*vti;
  double mu_max_ion = mi*pow(1.5*4*vti,2)/(2*B0);

  double t_end = 6.0e-8;
  int num_frames = 1;
  double write_phase_freq = 0.2; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gk_app_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .Lz     = Lz    ,
    .Ly     = Ly    ,
    .z_min = -Lz/2.,  .z_max = Lz/2.,
    .psi_min = psi_min,  .psi_max = psi_max,
    .psi_LCFS = psi_LCFS,
    .Lx = Lx,
    .Lx_core = Lx_core,
  
    .me = me,  .qe = qe,
    .mi = mi,  .qi = qi,
    .n0 = n0,  .Te0 = Te0,  .Ti0 = Ti0,
  
    .nuFrac = nuFrac,
    .nuElc = nuElc,  .nuIon = nuIon,
    .nuElcIon = nuElcIon,  .nuIonElc = nuIonElc,
  
    .n_src = n_src,  .Te_src = Te_src,  .Ti_src = Ti_src,
  
    .Nz = Nz,
    .Ny = Ny,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nx, Ny, Nz, Nvpar, Nmu},
    .poly_order = poly_order,
    .vpar_max_elc = vpar_max_elc,  .mu_max_elc = mu_max_elc,
    .vpar_max_ion = vpar_max_ion,  .mu_max_ion = mu_max_ion,

    .t_end = t_end,  .num_frames = num_frames,
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

  struct gk_app_ctx ctx = create_ctx(); // Context for init functions.

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
    .lower = { -ctx.vpar_max_elc, 0.0},
    .upper = {  ctx.vpar_max_elc, ctx.mu_max_elc}, 
    .cells = { cells_v[0], cells_v[1] },

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = density_init,
      .ctx_upar = &ctx,
      .upar= upar_elc,
      .ctx_temp = &ctx,
      .temp = temp_elc,      
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
      .temp_ref = ctx.Te0,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = density_src,
        .ctx_upar = &ctx,
        .upar= upar_elc_src,
        .ctx_temp = &ctx,
        .temp = temp_elc_src,  
      }, 
    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB, },
      { .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH, },
      { .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH, },
    },

    .num_diag_moments = 5,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP },
  };

  // ions
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.qi, .mass = ctx.mi,
    .vdim = ctx.vdim,
    .lower = { -ctx.vpar_max_ion, 0.0},
    .upper = {  ctx.vpar_max_ion, ctx.mu_max_ion}, 
    .cells = { cells_v[0], cells_v[1] },

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = density_init,
      .ctx_upar = &ctx,
      .upar= upar_ion,
      .ctx_temp = &ctx,
      .temp = temp_ion,      
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
      .temp_ref = ctx.Ti0,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = density_src,
        .ctx_upar = &ctx,
        .upar= upar_ion_src,
        .ctx_temp = &ctx,
        .temp = temp_ion_src,  
      }, 
    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB, },
      { .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH, },
      { .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH, },
    },
    
    .num_diag_moments = 5,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP },
  };

  // field
  struct gkyl_gyrokinetic_field field = {
    .gkfield_id = GKYL_GK_FIELD_ES,
    .poisson_bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}, },
      { .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_PERIODIC, .value = {0.0}, },
      { .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_PERIODIC, .value = {0.0}, },
    },
  };

  struct gkyl_efit_inp efit_inp = {
    // psiRZ and related inputs
    .filepath = "gyrokinetic/data/eqdsk/ltx_miller.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect = true,
  };

  struct gkyl_tok_geo_grid_inp grid_inp = {
    .ftype = GKYL_GEOMETRY_TOKAMAK_IWL,
    .rleft= 0.1,
    .rright = 0.7,
    .rmin = 0.1,
    .rmax = 0.7,
    .zmin = -0.35,
    .zmax = 0.35,
    .plate_spec = false,
  }; 

  // GK app
  struct gkyl_gk app_inp = {
    .cdim = ctx.cdim,
    .lower = { ctx.psi_min, -ctx.Ly/2.0, ctx.z_min },
    .upper = { ctx.psi_max, ctx.Ly/2.0, ctx.z_max },
    .cells = { cells_x[0], cells_x[1], cells_x[2] },
    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,

    .geometry = {
      .world = {0.0},
      .geometry_id = GKYL_GEOMETRY_TOKAMAK,
      .efit_info = efit_inp,
      .tok_grid_info = grid_inp,
      .has_LCFS = true,
      .x_LCFS = ctx.psi_LCFS, // Location of last closed flux surface.
    },

    .num_periodic_dir = 1,
    .periodic_dirs = {1},

    .num_species = 2,
    .species = { elc, ion },
    .field = field,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0], app_args.cuts[1] , app_args.cuts[2]},
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
