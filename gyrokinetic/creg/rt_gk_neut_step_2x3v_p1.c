#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_fem_parproj.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_gyrokinetic_run.h>
#include <gkyl_tok_geo.h>

#include <rt_arg_parse.h>

struct gk_step_ctx {
  int cdim, vdim; // Dimensionality.
  double chargeElc; // electron charge
  double massElc; // electron mass
  double chargeIon; // ion charge
  double massIon; // ion mass
  double Te; // electron temperature
  double Ti; // ion temperature
  double vtIon;
  double vtElc;
  double vtD0; // D0 thermal velocity.
  double nD0;
  double TD0;
  double nuElc; // electron collision frequency
  double nuIon; // ion collision frequency
  double nuD0;
  double nuFrac; // Factor to multiply collision frequencies
  double B0; // reference magnetic field
  double n0; // reference density
  double nsource;
  // Source parameters
  double T_source; // Source electron temperature
  double cx;
  double cz;
  // Simulation parameters
  int Nx; // Cell count (configuration space: x-direction).
  int Nz; // Cell count (configuration space: z-direction).
  int Nvpar; // Cell count (velocity space: parallel velocity direction).
  int Nmu; // Cell count (velocity space: magnetic moment direction).
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  double lower_x, upper_x; // Limits of the domain along x.
  double Lx; // Domain size (configuration space: x-direction).
  double Lz; // Domain size (configuration space: z-direction).
  double vpar_max_elc; // Velocity space extents in vparallel for electrons
  double mu_max_elc; // Velocity space extents in mu for electrons
  double vpar_max_ion; // Velocity space extents in vparallel for ions
  double mu_max_ion; // Velocity space extents in mu for ions
  double vpar_max_D0; // Domain boundary (D0 velocity space).
  double t_end; // end time
  int num_frames; // number of output frames
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

void
eval_density(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double x = xn[0], z = xn[1];
  double n0 = app->n0;
  double cx = app->cx;
  double cz = app->cz;
  double xcenter = 1.2014;
  double n = n0*exp(-(x-xcenter)*(x-xcenter)/2/cx/cx) * exp(-z*z/2/cz/cz);
  if (n/n0 < 1e-1)
    n = n0*1e-1;
  fout[0] = n;
}

void unit_density(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct gk_app_ctx *app = ctx;
  fout[0] = 1.0; 
}

void
eval_upar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
evalDensityInitD0(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct gk_step_ctx *app = ctx;
  double z = xn[0];

  double n0 = app->nD0;
  fout[0] = n0;
}

void
evalTempD0Init(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct gk_step_ctx *app = ctx;

  double TD0 = app->TD0;

  // Set electron temperature.
  fout[0] = TD0;
}

void
eval_udrift(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0; 
  fout[1] = 0.0;
  fout[2] = 0.0;
}



void
eval_temp_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double T = app->Te;
  fout[0] = T;
}

void
eval_temp_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double T = app->Ti;
  fout[0] = T;
}

void
eval_density_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double x = xn[0], z = xn[1];
  double nsource = app->nsource;
  fout[0] = nsource;
  //double cx = app->cx;
  //double cz = app->cz;
  //double xcenter = 1.2014;
  if (fabs(z) > app->Lz/20.0*(3.0/4.0) )
    fout[0] = nsource;
  else
    fout[0] = nsource*1.e-3;
}

void
eval_upar_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
eval_temp_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double n0 = app->n0;
  double T = app->T_source;
  fout[0] = T;
}

void
evalNuElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  fout[0] = app->nuElc;
}

void
evalNuIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  fout[0] = app->nuIon;
}

void
evalNuD0(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  fout[0] = app->nuD0;
}

double plasma_frequency(double n, double m)
{
  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  return sqrt(n*eV*eV/m/eps0);
}

struct gk_step_ctx
create_ctx(void)
{
  int cdim = 2, vdim = 2; // Dimensionality.

  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mi = 2.014*GKYL_PROTON_MASS; // ion mass
  double me = GKYL_ELECTRON_MASS;
  double qi = eV; // ion charge
  double qe = -eV; // electron charge

  double B0 = 2.51; // Magnetic field magnitude in Tesla
  double Te = 25.0 * GKYL_ELEMENTARY_CHARGE; // Electron temperature.
  double Ti = 25.0 * GKYL_ELEMENTARY_CHARGE; // Ion temperature.
  double n0 = 1.0e19; //  Reference number density (1 / m^3).
                             
  double TD0 = 10.0 * GKYL_ELEMENTARY_CHARGE;
  double nD0 = 1.0*n0;

  // Derived parameters.
  double vtIon = sqrt(Ti/mi);
  double vtElc = sqrt(Te/me);
  double vtD0 = sqrt(TD0/mi);

  // Source parameters.
  double nsource = 2.870523e25; // Source number density.
  double T_source = 2.0*Te;
  double cx = 0.0065612*9;
  double cz = 0.4916200;

  // Collision parameters.
  double nuFrac = 0.25;
  double logLambdaElc = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Te/eV);
  double nuElc = nuFrac*logLambdaElc*pow(eV, 4.0)*n0/(6.0*sqrt(2.0)*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(me)*(Te*sqrt(Te)));  // collision freq

  double logLambdaIon = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Ti/eV);
  double nuIon = nuFrac*logLambdaIon*pow(eV, 4.0)*n0/(12.0*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(mi)*(Ti*sqrt(Ti)));
  double nuD0 = 1e8;

  // Simulation box size (m).
  double lower_x = 0.934;
  double upper_x = 1.334;
  double Lx = upper_x - lower_x;
  double Lz = (M_PI-1e-14)*2.0;

  double vpar_max_elc = 4.0*vtElc;
  double mu_max_elc = 12*me*vtElc*vtElc/(2.0*B0);

  double vpar_max_ion = 4.0*vtIon;
  double mu_max_ion = 12*mi*vtIon*vtIon/(2.0*B0);

  double vpar_max_D0 = vtD0;


  // Number of cells.
  int Nx = 4;
  int Nz = 8;
  int Nvpar = 12;
  int Nmu = 8;

  double t_end = 1e-6; 
  double num_frames = 1;
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gk_step_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .chargeElc = qe, 
    .massElc = me, 
    .chargeIon = qi, 
    .massIon = mi,
    .Te = Te, 
    .Ti = Ti, 
    .TD0 = TD0,
    .vtIon = vtIon,
    .vtElc = vtElc,
    .vtD0 = vtD0,
    .nuElc = nuElc, 
    .nuIon = nuIon, 
    .nuD0 = nuD0,
    .nuFrac = nuFrac,
    .B0 = B0, 
    .n0 = n0, 
    .nD0 = nD0,
    .T_source = T_source, 
    .nsource = nsource,
    .cx = cx,
    .cz = cz,
    .lower_x = lower_x,
    .upper_x = upper_x,
    .Lx = Lx,
    .Lz = Lz, 
    .vpar_max_elc = vpar_max_elc, 
    .mu_max_elc = mu_max_elc, 
    .vpar_max_ion = vpar_max_ion, 
    .mu_max_ion = mu_max_ion, 
    .vpar_max_D0= vpar_max_D0,
    .Nx = Nx,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nx, Nz, Nvpar, Nmu},
    .t_end = t_end, 
    .num_frames = num_frames, 
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

  struct gk_step_ctx ctx = create_ctx(); // Context for init functions.

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  // neutral Deuterium
  struct gkyl_gyrokinetic_neut_species D0 = {
    .name = "D0", .mass = ctx.massIon,
    .vdim = ctx.vdim+1,
    .lower = { -ctx.vpar_max_D0*2.3, -ctx.vpar_max_D0*25.1, -ctx.vpar_max_D0*9.0},
    .upper = { ctx.vpar_max_D0*2.3, ctx.vpar_max_D0*25.1, ctx.vpar_max_D0*9.0 },
    .cells = { 24,24,24 },
    .is_static = false,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = evalDensityInitD0,
      .ctx_upar = &ctx,
      .udrift= eval_udrift,
      .ctx_temp = &ctx,
      .temp = evalTempD0Init,      
      //.correct_all_moms = true,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_NEUTRAL,
    },

    .source = { 
       .source_id = GKYL_PROJ_SOURCE, 
       .num_sources = 1, 
       .projection[0] = { 
         .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
         .density = eval_density_source, 
         .ctx_density = &ctx, 
         .udrift = eval_udrift, 
         .ctx_upar = &ctx, 
         .temp = evalTempD0Init, 
         .ctx_temp = &ctx,
       }, 
    }, 

//    .collisions =  { 
//      .collision_id = GKYL_BGK_COLLISIONS, 
//      .self_nu = evalNuD0, 
//      .ctx = &ctx, 
//      //.is_implicit = true, 
//    }, 

//    .positivity = {
//      .type = GKYL_GK_POSITIVITY_SHIFT,
//      .write_diagnostics = true,
//    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB, },
      { .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB, },
      { .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB, },
    },
    
    .num_diag_moments = 4,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1_FROM_H, GKYL_F_MOMENT_ENERGY, GKYL_F_MOMENT_LTE},
  };

  // Field.
  struct gkyl_gyrokinetic_field field = {
    .is_static = true,
    .zero_init_field = true,
    .poisson_bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}, },
    },
    .time_rate_diagnostics = false,
  };

  struct gkyl_efit_inp efit_inp = {
    // psiRZ and related inputs
    .filepath = "gyrokinetic/data/eqdsk/step.geqdsk",   // equilibrium to use
    .rz_poly_order = 2,                       // polynomial order for psi(R,Z) used for field line tracing
    .flux_poly_order = 1,                     // polynomial order for fpol(psi)
    .reflect = true,                          // Reflect lower half of psi(R,Z) for up-down symmetry
  };

  struct gkyl_tok_geo_grid_inp grid_inp = {
    .ftype = GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT, // type of geometry
    .rclose = 6.2,            // closest R to region of interest
    .rright= 6.2,             // Closest R to outboard SOL
    .rleft= 2.0,              // closest R to inboard SOL
    .rmin = 1.1,              // smallest R in machine
    .rmax = 6.2,              // largest R in machine
    .zmin = -8.3,         // Z of lower divertor plate
    .zmax = 8.3,          // Z of upper divertor plate
  };

  // GK app
  struct gkyl_gk app_inp = {

    .cdim = ctx.cdim,
    .lower = { ctx.lower_x, -ctx.Lz/2.0 },
    .upper = { ctx.upper_x,  ctx.Lz/2.0 },
    .cells = { cells_x[0], cells_x[1] },
    .poly_order = 1,
    .basis_type = app_args.basis_type,
    .cfl_frac = 1.0,
    .cfl_frac_omegaH = 1e10,

    .geometry = {
      .world = {0.0},
      .geometry_id = GKYL_GEOMETRY_TOKAMAK,
      .efit_info = efit_inp,
      .tok_grid_info = grid_inp,
    },

    .num_periodic_dir = 0,
    .periodic_dirs = {  },

    .num_species = 0,
    .species = { },

    .num_neut_species = 1,
    .neut_species = { D0 },

    .field = field,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0], app_args.cuts[1] },
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
      .write_phase_freq = 1,
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
