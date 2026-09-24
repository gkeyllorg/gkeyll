#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_efit.h>
#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_gyrokinetic_run.h>
#include <gkyl_tok_geo.h>

#include <rt_arg_parse.h>

void shaped_pfunc_lower_outer(double s, double* RZ){
  RZ[0] = 3.5+2.0*s;
  RZ[1] = -8.29;
}

void shaped_pfunc_upper_outer(double s, double* RZ){
  RZ[0] = 3.5+2.0*s;
  RZ[1] = 8.29;
}

void shaped_pfunc_upper_inner(double s, double* RZ){
  RZ[0] = 1.651 + (1.8 - 1.651)*s;
  RZ[1] = 6.331 + (6.777 - 6.331)*s;
}

void shaped_pfunc_lower_inner(double s, double* RZ){
  RZ[0] = 1.65 + (1.8 - 1.65)*s;
  RZ[1] = -(6.33 + (6.777 - 6.33)*s);
}

struct gkyl_gk_block_geom*
create_gk_block_geom(void)
{
  struct gkyl_gk_block_geom *bgeom = gkyl_gk_block_geom_new(1, 3);

  /* Block layout and coordinates

   x  
   ^  
   |
   1  +------------------+------------------+------------------+
   |  |b0                |b1                |b2                |
   |  |lower outer SOL   |middle outer sol  |upper outer sol   |
   |  |                  |                  |                  |
   0  +------------------+------------------+------------------+
      0 -----------1------------2------------3 -> z

      Edges that touch coincide are physically connected unless
      otherwise indicated by a special symbol. Edges with a special
      symbol such as o,x,%, or % are instead connected to the other
      edge with the same symbol. Edges that do not coincide with
      another edge are a physical boundary.
  */  



  struct gkyl_efit_inp efit_inp = {
    // psiRZ and related inputs
    .filepath = "gyrokinetic/data/eqdsk/step.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect = true,
  };

  struct gkyl_efit *efit = gkyl_efit_new(&efit_inp);
  double psisep = efit->psisep;
  psisep = 1.4688;
  gkyl_efit_release(efit);
  double psi_up_core = 1.8;
  double psi_up_pf = 1.8;
  double psi_lo_outer_sol = 0.934;
  double psi_lo_inner_sol = 1.45;

  int npsi_outer_sol = 10;

  double ntheta_lower  = 4;
  double ntheta_middle = 8;
  double ntheta_upper  = 4;

  // Note that for tokamak multi-block simulations, 
  // these theta limits are just placeholders and will be 
  // reset in the multi-block app.
  double theta_lo = -M_PI + 1e-14, theta_up = M_PI - 1e-14;

  // block 0. Lower outer SOL.
  gkyl_gk_block_geom_set_block(bgeom, 0, &(struct gkyl_gk_block_geom_info) {
      .lower = { theta_lo },
      .upper = { theta_up },
      .cells = { ntheta_lower},
      .cuts = { 1 },
      .geometry = {
        .world = {1.2014, 0.0},
        .geometry_id = GKYL_GEOMETRY_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO,
          .rright = 6.2,
          .rleft = 1.1,
          .rmin = 2.1,
          .rmax = 6.2,
          .zmin = -8.29,
          .zmax = 8.29,
          .plate_spec = false,
          .plate_func_lower = shaped_pfunc_lower_outer,
          .plate_func_upper = shaped_pfunc_upper_outer,
        }
      },
      
      .connections[0] = { // z-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 1, .dir = 0, .edge = GKYL_LOWER_POSITIVE}
      }
    }
  );

  // block 1. Middle outer SOL.
  gkyl_gk_block_geom_set_block(bgeom, 1, &(struct gkyl_gk_block_geom_info) {
      .lower = { theta_lo },
      .upper = { theta_up },
      .cells = { ntheta_middle},
      .cuts = { 1 },
      .geometry = {
        .world = {1.2014, 0.0},
        .geometry_id = GKYL_GEOMETRY_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID,
          .rright = 6.2,
          .rleft = 1.1,
          .rmin = 2.1,
          .rmax = 6.2,
          .plate_spec = false,
          .plate_func_lower = shaped_pfunc_lower_outer,
          .plate_func_upper = shaped_pfunc_upper_outer,
          .zmin = -8.29,
          .zmax = 8.29,
        }
      },
      
      .connections[0] = { // z-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 2, .dir = 0, .edge = GKYL_LOWER_POSITIVE}
      }
    }
  );

  // block 2. Upper outer SOL.
  gkyl_gk_block_geom_set_block(bgeom, 2, &(struct gkyl_gk_block_geom_info) {
      .lower = { theta_lo },
      .upper = { theta_up },
      .cells = { ntheta_upper},
      .cuts = { 1 },
      .geometry = {
        .world = {1.2014, 0.0},
        .geometry_id = GKYL_GEOMETRY_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP,
          .rright = 6.2,
          .rleft = 1.1,
          .rmin = 2.1,
          .rmax = 6.2,
          .plate_spec = false,
          .plate_func_lower = shaped_pfunc_lower_outer,
          .plate_func_upper = shaped_pfunc_upper_outer,
          .zmin = -8.29,
          .zmax = 8.29,
        }
      },
      
      .connections[0] = { // z-direction connections
        { .bid = 1, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 2, .dir = 0, .edge = GKYL_PHYSICAL} // physical boundary
      }
    }
  );

  return bgeom;
}

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
  double nuElc; // electron collision frequency
  double nuIon; // ion collision frequency
  double nuFrac; // Factor to multiply collision frequencies
  double B0; // reference magnetic field
  double n0; // reference density
  double nsource;
  // Source parameters
  double T_source; // Source electron temperature
  double cx;
  double cz;
  double k_perp; // Perpendicular wavenumber (for Poisson solver).
  // Simulation parameters
  int Nz; // Cell count (configuration space: z-direction).
  int Nvpar; // Cell count (velocity space: parallel velocity direction).
  int Nmu; // Cell count (velocity space: magnetic moment direction).
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  double vpar_max_elc; // Velocity space extents in vparallel for electrons
  double mu_max_elc; // Velocity space extents in mu for electrons
  double vpar_max_ion; // Velocity space extents in vparallel for ions
  double mu_max_ion; // Velocity space extents in mu for ions

  double t_end; // End time.
  int num_frames; // Number of output frames.
  double write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct gk_step_ctx
create_ctx(void)
{
  int cdim = 1, vdim = 2; // Dimensionality.

  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mi = 2.014*GKYL_PROTON_MASS; // ion mass
  double me = GKYL_ELECTRON_MASS;
  double qi = eV; // ion charge
  double qe = -eV; // electron charge

  double Te = 364*eV;
  double Ti = 534*eV;
  double B0 = 2.51; // Magnetic field magnitude in Tesla
  double n0 = 3.0e19; // Particle density in 1/m^3
                             
  // Derived parameters.
  double vtIon = sqrt(Ti/mi);
  double vtElc = sqrt(Te/me);

  double k_perp_rho_s = 0.3; // Product of perpendicular wavenumber and ion-sound gyroradius.
  double c_s = sqrt(Te / mi); // Sound speed.
  double omega_ci = fabs(qi* B0 / mi); // Ion cyclotron frequency.
  double rho_s = c_s / omega_ci; // Ion-sound gyroradius.

  double k_perp = k_perp_rho_s / rho_s; // Perpendicular wavenumber (for Poisson solver).

  // Source parameters.
  double nsource = 3.9e23; // peak source rate in particles/m^3/s 
  double T_source = 1037.0*eV;
  double cz = 0.8;
  double cx = 0.0065612*4;

  // Collision parameters.
  double nuFrac = 0.25;
  double logLambdaElc = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Te/eV);
  double nuElc = nuFrac*logLambdaElc*pow(eV, 4.0)*n0/(6.0*sqrt(2.0)*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(me)*(Te*sqrt(Te)));  // collision freq

  double logLambdaIon = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Ti/eV);
  double nuIon = nuFrac*logLambdaIon*pow(eV, 4.0)*n0/(12.0*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(mi)*(Ti*sqrt(Ti)));

  // Simulation box size (m).
  double vpar_max_elc = 6.0*vtElc;
  double mu_max_elc = 12*me*vtElc*vtElc/(2.0*B0);

  double vpar_max_ion = 6.0*vtIon;
  double mu_max_ion = 12*mi*vtIon*vtIon/(2.0*B0);

  // Number of cells.
  int Nz = 4;
  int Nvpar = 12;
  int Nmu = 8;

  double t_end = 1.0e-4; 
  double num_frames = 10;
  double write_phase_freq = 0.2; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
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
    .vtIon = vtIon,
    .vtElc = vtElc,
    .nuElc = nuElc, 
    .nuIon = nuIon, 
    .nuFrac = nuFrac,
    .B0 = B0, 
    .n0 = n0, 
    .k_perp = k_perp,
    .T_source = T_source, 
    .nsource = nsource,
    .cx = cx,
    .cz = cz,
    .vpar_max_elc = vpar_max_elc, 
    .mu_max_elc = mu_max_elc, 
    .vpar_max_ion = vpar_max_ion, 
    .mu_max_ion = mu_max_ion, 
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nz, Nvpar, Nmu},
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
initDensity(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = 1.2014, z = xn[0];

  struct gk_step_ctx *app = ctx;
  double n0 = app->n0;
  double cx = app->cx;
  double cz = app->cz;
  double xcenter = 1.2014;
  double n = n0*exp(-(x-xcenter)*(x-xcenter)/2/cx/cx) * exp(-z*z/2/cz/cz);
  if (n/n0 < 1e-1)
    n = n0*1e-1;
  fout[0] = n;
}

void
sourceDensity(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = 1.2014, z = xn[0];

  struct gk_step_ctx *app = ctx;
  double nsource = app->nsource;
  double cx = app->cx;
  double cz = app->cz;
  double xcenter = 1.2014;
  double n = nsource*exp(-(x-xcenter)*(x-xcenter)/2/cx/cx) * exp(-z*z/2/cz/cz);
  if (n/nsource < 1e-3)
    n = nsource*1e-3;
  fout[0] = n;
}

void
initTempElc(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  double T = input->Te;
  fout[0] = T;
}

void
initTempIon(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  double T = input->Ti;
  fout[0] = T;
}

void
initUpar(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = 0.0;
}

void
sourceTemp(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  double T = input->T_source;
  fout[0] = T;
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

  struct gk_step_ctx ctx = create_ctx(); // Context for init functions.

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  // construct block geometry
  struct gkyl_gk_block_geom *bgeom = create_gk_block_geom();
  int nblocks = gkyl_gk_block_geom_num_blocks(bgeom);

  // Elc Species
  // all data is common across blocks
  struct gkyl_gyrokinetic_multib_species_pb elc_blocks[1];
  elc_blocks[0] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = initDensity,
      .ctx_upar = &ctx,
      .upar = initUpar,
      .ctx_temp = &ctx,
      .temp = initTempElc,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = sourceDensity,
        .ctx_upar = &ctx,
        .upar = initUpar,
        .ctx_temp = &ctx,
        .temp = sourceTemp,      
      }, 
    },

  };


  struct gkyl_gyrokinetic_bc elc_phys_bcs[] = {
    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    { .bidx = 2, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH }
  };

  struct gkyl_gyrokinetic_multib_species elc = {
    .name = "elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .vdim = ctx.vdim,
    .lower = { -ctx.vpar_max_elc, 0.0},
    .upper = {  ctx.vpar_max_elc, ctx.mu_max_elc}, 
    .cells = { cells_v[0], cells_v[1] },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES_NO_BY,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .nu_frac = ctx.nuFrac,
      .den_ref = ctx.n0, // Density used to calculate coulomb logarithm
      .temp_ref = ctx.Te, // Temperature used to calculate coulomb logarithm
      .bmag_ref =  2.51,
      .num_cross_collisions = 1,
      .collide_with = { "ion" },
    },

    .num_diag_moments = 7,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP },

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
      .ctx_density = &ctx,
      .density = initDensity,
      .ctx_upar = &ctx,
      .upar = initUpar,
      .ctx_temp = &ctx,
      .temp = initTempIon,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = sourceDensity,
        .ctx_upar = &ctx,
        .upar = initUpar,
        .ctx_temp = &ctx,
        .temp = sourceTemp,      
      }, 
    },

  };

  struct gkyl_gyrokinetic_bc ion_phys_bcs[] = {
    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    { .bidx = 2, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH }
  };

  struct gkyl_gyrokinetic_multib_species ion = {
    .name = "ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .vdim = ctx.vdim,
    .lower = { -ctx.vpar_max_ion, 0.0},
    .upper = {  ctx.vpar_max_ion, ctx.mu_max_ion}, 
    .cells = { cells_v[0], cells_v[1] },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES_NO_BY,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .nu_frac = ctx.nuFrac,
      .den_ref = ctx.n0, // Density used to calculate coulomb logarithm
      .temp_ref = ctx.Ti, // Temperature used to calculate coulomb logarithm
      .bmag_ref =  2.51,
      .num_cross_collisions = 1,
      .collide_with = { "elc"},
    },

    .num_diag_moments = 7,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP },

    .duplicate_across_blocks = true,
    .blocks = ion_blocks,
    .num_physical_bcs = 2,
    .bcs = ion_phys_bcs,
  };

  // Field object
  struct gkyl_gyrokinetic_multib_field_pb field_blocks[1];
  field_blocks[0] = (struct gkyl_gyrokinetic_multib_field_pb) {
    .polarization_bmag = 2.51,
    .kperpSq = ctx.k_perp*ctx.k_perp,
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
