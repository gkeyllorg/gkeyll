#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_gyrokinetic_run.h>
#include <gkyl_tok_geo.h>

#include <rt_arg_parse.h>

struct gk_solovev_ctx {
  int cdim, vdim; // Dimensionality.

  double chargeElc; // Electron charge.
  double massElc; // Electron mass.
  double chargeIon; // Ion charge.
  double massIon; // Ion mass.
  double Te; // Electron temperature.
  double Ti; // Ion temperature.
  double c_s; // Sound speed.
  double nuElc; // Electron collision frequency.
  double nuIon; // Ion collision frequency.
  double nuElcIon; // Electron-ion collision frequency.
  double nuIonElc; // Ion-electron collision frequency.
  double B0; // Reference magnetic field.
  double n0; // Reference density.
  // Source parameters.
  double lambda_source;
  double x_source;
  // Domain parameters.            
  int Nx, Ny, Nz; // Number of cells in x,y,z.
  int Nvpar, Nmu; // Number of cells in vpar,mu.
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  double Lx; // Box size in x
  double Ly; // Box size in y
  double Lz; // Box size in z
  double vpar_max_elc; // Velocity space extents in vparallel for electrons
  double mu_max_elc; // Velocity space extents in mu for electrons
  double vpar_max_ion; // Velocity space extents in vparallel for ions
  double mu_max_ion; // Velocity space extents in mu for ions

  double t_end; // end time
  int num_frames; // number of output frames
  double write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

void
eval_density(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_solovev_ctx *app = ctx;
  double x = xn[0], y = xn[1], z = xn[2];
  double lambda_source = app->lambda_source;
  double x_source = app->x_source;
  double Lz = app->Lz;
  double Ls = Lz/4;
  double floor = 0.1;


  // find source density at z = 0
  double source_density = 0;
  double source_floor = 1e-10;
  if (x > x_source)
     source_floor = 1e-2; // higher floor to left of source peak
  source_density = fmax(exp(-(x-x_source)*(x-x_source)/((2*lambda_source)*(2*lambda_source))), source_floor);

  // find source temp at z = 0
  double source_temp = 0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  if (x > x_source - 3*lambda_source)
     source_temp = 80*eV;
  else
     source_temp = 30*eV;

  // now compute initial desity
  double effective_source = 3.800419e+23*fmax(source_density, floor);
  double c_ss = sqrt(5.0/3.0*source_temp/app->massIon);
  double n_peak = 4*sqrt(5)/3/c_ss*Ls*effective_source/2;
  double perturb = 0;
  if (fabs(z) <= Ls)
     fout[0]  = n_peak*(1+sqrt(1-(z/Ls)*(z/Ls)))/2*(1+perturb);
  else
     fout[0] = n_peak/2*(1+perturb);
}

void
eval_upar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
eval_temp_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_solovev_ctx *app = ctx;
  double x = xn[0], y = xn[1], z = xn[2];
  double lambda_source = app->lambda_source;
  double x_source = app->x_source;
  double eV = GKYL_ELEMENTARY_CHARGE;
  if (x > x_source - 3*lambda_source)
     fout[0] = 50*eV;
  else
     fout[0] = 20*eV;
}

void
eval_temp_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_solovev_ctx *app = ctx;
  double x = xn[0], y = xn[1], z = xn[2];
  double lambda_source = app->lambda_source;
  double x_source = app->x_source;
  double eV = GKYL_ELEMENTARY_CHARGE;
  if (x > x_source - 3*lambda_source)
     fout[0] = 50*eV;
  else
     fout[0] = 20*eV;
}

void
eval_density_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_solovev_ctx *app = ctx;
  double x = xn[0], y = xn[1], z = xn[2];
  double lambda_source = app->lambda_source;
  double x_source = app->x_source;
  double Lz = app->Lz;
  double source_floor = 1e-10;
  if (x > x_source)
     source_floor = 1e-2; // higher floor to left of source peak
  if (fabs(z) < Lz/4)
     fout[0] = 3.800419e+23*fmax(exp(-(x-x_source)*(x-x_source)/((2*lambda_source)*(2*lambda_source))), source_floor);
  else
     fout[0] = 3.800419e+23*1e-40;
}

void
eval_upar_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
eval_temp_elc_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_solovev_ctx *app = ctx;
  double x = xn[0], y = xn[1], z = xn[2];
  double lambda_source = app->lambda_source;
  double x_source = app->x_source;
  double eV = GKYL_ELEMENTARY_CHARGE;
  if (x > x_source - 3*lambda_source)
     fout[0] = 80*eV;
  else
     fout[0] = 30*eV;
}

void
eval_temp_ion_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_solovev_ctx *app = ctx;
  double x = xn[0], y = xn[1], z = xn[2];
  double lambda_source = app->lambda_source;
  double x_source = app->x_source;
  double eV = GKYL_ELEMENTARY_CHARGE;
  if (x > x_source - 3*lambda_source)
     fout[0] = 80*eV;
  else
     fout[0] = 30*eV;
}

void
evalNuElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_solovev_ctx *app = ctx;
  fout[0] = app->nuElc;
}

void
evalNuIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_solovev_ctx *app = ctx;
  fout[0] = app->nuIon;
}

void
evalNuElcIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_solovev_ctx *app = ctx;
  fout[0] = app->nuElcIon;
}

void
evalNuIonElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_solovev_ctx *app = ctx;
  fout[0] = app->nuIonElc;
}

struct gk_solovev_ctx
create_ctx(void)
{
  int cdim = 3, vdim = 2; // Dimensionality.

  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mi = 2.014*GKYL_PROTON_MASS; // ion mass
  double me = GKYL_ELECTRON_MASS;
  double qi = eV; // ion charge
  double qe = -eV; // electron charge

  double Te = 40.0*eV;
  double Ti = 40.0*eV;
  double B0 = 0.55; // Magnetic field magnitude in Tesla
  double n0 = 7.0e18; // Particle density in 1/m^3

  // Derived parameters.
  double vtIon = sqrt(Ti/mi);
  double vtElc = sqrt(Te/me);
  double c_s = sqrt(Te/mi);
  double omega_ci = fabs(qi*B0/mi);
  double rho_s = c_s/omega_ci;


  // Collision parameters.
  double nuFrac = 0.1;
  double logLambdaElc = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Te/eV);
  double nuElc = nuFrac*logLambdaElc*pow(eV, 4.0)*n0/(6.0*sqrt(2.0)*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(me)*(Te*sqrt(Te)));

  double logLambdaIon = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Ti/eV);
  double nuIon = nuFrac*logLambdaIon*pow(eV, 4.0)*n0/(12.0*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(mi)*(Ti*sqrt(Ti)));
  double nuElcIon = nuElc*sqrt(2.0);
  double nuIonElc = nuElcIon*(me/mi);

  // Simulation box size (psi, alpha, theta).
  double q0 = 2.0;
  double r0 = 0.40705706492831;
  double Lx = 0.02;
  double Ly = 50 * rho_s * q0/r0 ; // should be 0.107890816895
  double Lz = (M_PI-1e-14)*2.0 ; // Domain size (configuration space: z-direction).

  // Source parameters.
  double x_source = -0.07;
  double lambda_source = 0.03*Lx;

  // Velocity Grid
  double vpar_max_elc = 4.0*vtElc;
  double mu_max_elc = 0.75*me*(4.0*vtElc)*(4.0*vtElc)/(2.0*B0);

  double vpar_max_ion = 4.0*vtIon;
  double mu_max_ion = 0.75*mi*(4.0*vtIon)*(4.0*vtIon)/(2.0*B0);

  int Nx = 4; // Number of cells in x.
  int Ny = 2; // Number of cells in y.
  int Nz = 8; // Number of cells in z.
  int Nvpar = 6; // Number of cells in vpar.
  int Nmu = 4; // Number of cells in mu.

  double t_end = 4.0e-7; 
  double num_frames = 1;
  double write_phase_freq = 0.2; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gk_solovev_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .chargeElc = qe, 
    .massElc = me, 
    .chargeIon = qi, 
    .massIon = mi,
    .Te = Te, 
    .Ti = Ti, 
    .c_s = c_s, 
    .nuElc = nuElc, 
    .nuIon = nuIon, 
    .nuElcIon = nuElcIon, 
    .nuIonElc = nuIonElc, 
    .B0 = B0, 
    .n0 = n0, 
    .Lx = Lx, 
    .Ly = Ly,  
    .Lz = Lz, 
    .lambda_source = lambda_source,
    .x_source = x_source,
    .vpar_max_elc = vpar_max_elc, 
    .mu_max_elc = mu_max_elc, 
    .vpar_max_ion = vpar_max_ion, 
    .mu_max_ion = mu_max_ion, 
    .Nx = Nx,
    .Ny = Ny,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nx, Ny, Nz, Nvpar, Nmu},
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
  if (app_args.use_mpi) MPI_Init(&argc, &argv);
#endif

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct gk_solovev_ctx ctx = create_ctx(); // Context for init functions.

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
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .vdim = ctx.vdim,
    .lower = { -ctx.vpar_max_elc, 0.0},
    .upper = {  ctx.vpar_max_elc, ctx.mu_max_elc}, 
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_elc,      
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
        .ctx_density = &ctx,
        .density = eval_density_source,
        .ctx_upar = &ctx,
        .upar= eval_upar_source,
        .ctx_temp = &ctx,
        .temp = eval_temp_elc_source,      
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
        
    .num_diag_moments = 7,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP },
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    .time_rate_diagnostics = true,

    .boundary_flux_diagnostics = {
      .num_diag_moments = 1,
      .diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
//      .time_integrated = true,
    },
  };

  // ions
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .vdim = ctx.vdim,
    .lower = { -ctx.vpar_max_ion, 0.0},
    .upper = {  ctx.vpar_max_ion, ctx.mu_max_ion}, 
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_ion,      
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

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = eval_density_source,
        .ctx_upar = &ctx,
        .upar= eval_upar_source,
        .ctx_temp = &ctx,
        .temp = eval_temp_ion_source,      
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
    
    .num_diag_moments = 7,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP },
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    .time_rate_diagnostics = true,

    .boundary_flux_diagnostics = {
      .num_diag_moments = 1,
      .diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
//      .time_integrated = true,
    },
  };

  // field
  struct gkyl_gyrokinetic_field field = {
    .poisson_bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
      { .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_PERIODIC },
      { .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_PERIODIC },
    },
    .time_rate_diagnostics = true,
  };

  struct gkyl_efit_inp efit_inp = {
    // psiRZ and related inputs
    .filepath = "gyrokinetic/data/eqdsk/solovev.geqdsk", // equilibrium to use
    .rz_poly_order = 2,                        // polynomial order for psi(R,Z) used for field line tracing
    .flux_poly_order = 1,                      // polynomial order for fpol(psi)
  };
  
  struct gkyl_tok_geo_grid_inp grid_inp = {
    .ftype = GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT,    // type of geometry
    .rclose = 3.0,               // closest R to region of interest
    .rright = 3.0,               // Closest R to outboard SOL
    .rleft = 0.1,                // closest R to inboard SOL
    .rmin = 0.1,                 // smallest R in machine
    .rmax = 3.5,                 // largest R in machine
    .zmin = -1.5,                // Z of lower divertor plate
    .zmax = 1.5,                 // Z of upper divertor plate
  }; 

  // GK app
  struct gkyl_gk app_inp = {

    .cdim = ctx.cdim,
    .lower = { -0.08, -ctx.Ly/2.0, -ctx.Lz/2.0 },
    .upper = { -0.06, ctx.Ly/2.0, ctx.Lz/2.0 },
    .cells = { cells_x[0], cells_x[1], cells_x[2] },
    .poly_order = 1,
    .basis_type = app_args.basis_type,

    .geometry = {
      .geometry_id = GKYL_GEOMETRY_TOKAMAK,
      .efit_info = efit_inp,
      .tok_grid_info = grid_inp,
    },

    .num_periodic_dir = 1,
    .periodic_dirs = { 1 },

    .num_species = 2,
    .species = { elc, ion },
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
