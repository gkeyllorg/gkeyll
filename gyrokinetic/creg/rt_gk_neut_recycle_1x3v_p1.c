#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_fem_parproj.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_gyrokinetic_run.h>

#include <rt_arg_parse.h>

struct sheath_ctx
{
  int cdim, vdim; // Dimensionality.
  
  // Physical constants (using non-normalized physical units).
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.
  double mass_ion; // Proton mass.
  double charge_ion; // Proton charge.
  double mass_neut; // Neutral mass.

  double Te; // Electron temperature.
  double Ti; // Ion temperature.
  double Tn; // Neutral temperature.
  double n0; // Reference number density (1 / m^3).

  double rec_frac; // Recycling fraction for neutral BCs.
  double Tn_recycle; // Neutral recycling temperature.

  double k_perp_rho_s; // Product of perpendicular wavenumber and ion-sound gyroradius.
 
  double B0; // Reference magnetic field strength (Tesla).

  double c_s; // Sound speed.
  double vte; // Electron thermal velocity.
  double vti; // Ion thermal velocity.
  double vtn; // Neutral thermal velocity.
  double omega_ci; // Ion cyclotron frequency.
  double rho_s; // Ion-sound gyroradius.

  double k_perp; // Perpendicular wavenumber (for Poisson solver).

  double n_src; // Source number density.
  double T_src; // Source temperature.

  int Nz; // Number of cells along magnetic field.
  int Nvpar; // Number of cells in vpar.
  int Nmu; // Number of cells in mu.
  int Nv; // Number of cells in neutral v.
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  double Lz; // Domain size along z.
  double vpar_max_elc; // Maximum electron vpar.
  double mu_max_elc; // Maximum electron mu.
  double vpar_max_ion; // Maximum ion vpar.
  double mu_max_ion; // Maximum ion mu.
  double v_max_neut; // Maximum neutral velocity.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct sheath_ctx
create_ctx(void)
{
  int cdim = 1, vdim = 2; // Dimensionality.

  // Physical constants (using non-normalized physical units).
  double eV = GKYL_ELEMENTARY_CHARGE; // Elementary charge.
  double mass_elc = GKYL_ELECTRON_MASS; // Electron mass.
  double mass_ion = GKYL_PROTON_MASS; // Proton mass.
  double mass_neut = mass_ion; // Neutral mass.
  double charge_elc = -eV; // Electron charge.
  double charge_ion =  eV; // Proton charge.

  double Te = 30.0 * eV; // Electron temperature.
  double Ti = 60.0 * eV; // Ion temperature.
  double Tn = 60.0 * eV; // Neutral temperature
  double n0 = 5.0e18; //  Reference number density (1 / m^3).

  double rec_frac = 1.0; // Recycling coefficient for neutral BCs.
  double Tn_recycle = 10.0 * eV; // Neutral recycling temperature

  double k_perp_rho_s = 0.2; // Product of perpendicular wavenumber and ion-sound gyroradius.

  // Derived physical quantities (using non-normalized physical units).
  double B0 = 0.5; // Reference magnetic field strength (Tesla).

  double vte = sqrt(Te / mass_elc); // Electron thermal velocity.
  double vti = sqrt(Ti / mass_ion); // Ion thermal velocity.
  double vtn = sqrt(Tn / mass_neut); // Neutral thermal velocity.

  double c_s = sqrt(Te / mass_ion); // Sound speed.
  double omega_ci = fabs(charge_ion / mass_ion); // Ion cyclotron frequency.
  double rho_s = c_s / omega_ci; // Ion-sound gyroradius.

  double k_perp = k_perp_rho_s / rho_s; // Perpendicular wavenumber (for Poisson solver).

  double n_src = pow(n0,2.0)*8e-21; // Source number density.
  double T_src = 10.0*eV; // Source temperature.

  int Nz = 224; // Number of cells along magnetic field.
  int Nvpar = 16; // Number of cells in vpar.
  int Nmu = 4; // Number of cells in mu.
  int Nv = 16; // Number of cells in neutral v.

  double Lz = 40.0; // Domain size along z.
  double vpar_max_elc = 4.0*vte; // Maximum electron vpar.
  double mu_max_elc = 12.0*mass_elc*pow(vte,2)/(2.0*B0); // Maximum electron mu.
  double vpar_max_ion = 4.0 * vti; // Maximum ion vpar.
  double mu_max_ion = 12.0*mass_ion*pow(vti,2)/(2.0*B0); // Maximum ion mu.
  double v_max_neut = 4.0 * vtn; // Maximum neutral velocity.
 
  double t_end = .00625*Lz/c_s; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.
  
  struct sheath_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .mass_ion = mass_ion,
    .charge_ion = charge_ion,
    .mass_neut = mass_neut,
    .Te = Te,
    .Ti = Ti,
    .Tn = Tn,
    .n0 = n0,
    .rec_frac = rec_frac,
    .Tn_recycle = Tn_recycle,
    .k_perp_rho_s = k_perp_rho_s,
    .B0 = B0,
    .vte = vte,
    .vti = vti,
    .vtn = vtn,
    .c_s = c_s,
    .omega_ci = omega_ci,
    .rho_s = rho_s,
    .k_perp = k_perp,
    .n_src = n_src,
    .T_src = T_src,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .Nv = Nv,
    .cells = {Nz, Nvpar, Nmu},
    .Lz = Lz,
    .vpar_max_elc = vpar_max_elc,
    .mu_max_elc = mu_max_elc,
    .vpar_max_ion = vpar_max_ion,
    .mu_max_ion = mu_max_ion,
    .v_max_neut = v_max_neut,
    .t_end = t_end,
    .num_frames = num_frames,
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };

  return ctx;
}

void
evalSourceDensityInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Neutral source rate.
  struct sheath_ctx *app = ctx;
  double z = xn[0];

  double n_src = app->n_src;

  fout[0] = n_src;
}

void
evalSourceUdriftInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Neutral source drift velocity.
  fout[0] = 0.0; 
  fout[1] = 0.0;
  fout[2] = 0.0;
}

void
evalSourceTempInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Neutral source temperature.
  struct sheath_ctx *app = ctx;

  double T_src = app->T_src;

  fout[0] = T_src;
}

void
evalDensityInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Electron/ion density.
  struct sheath_ctx *app = ctx;
  double z = xn[0];

  double n = app->n0;
  fout[0] = n;
}

static double
sech(double x)
{
  // Hyperbolic secant.
  return 1.0/cosh(x);
}

void
evalDensityNeutInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Neutral density.
  struct sheath_ctx *app = ctx;
  double z = xn[0];
  double n0 = app->n0;
  double Lz = app->Lz;
  double den = 0.0;

  double z0 = 0.2;
  double flr = 1e-6;

  // Set number density.
  if (z <= 0) {
    den = n0*(pow(sech((-Lz/2-z)/z0),2) + flr);
  }
  else {
    den = n0*(pow(sech((Lz/2-z)/z0),2) + flr);
  }
  fout[0] = den;
}

void
evalUparInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Electron/ion parallel drift speed.
  struct sheath_ctx *app = ctx;
  double z = xn[0]; 
  fout[0] = app->c_s*z/(app->Lz/2.);
}

void
evalUdriftInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  // Neutral drift velocity.
  fout[0] = 0.0; 
  fout[1] = 0.0;
  fout[2] = 0.0;
}

void
evalTempElcInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Electron temperature.
  struct sheath_ctx *app = ctx;

  double Te = app->Te;

  fout[0] = Te;
}

void
evalTempIonInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Ion temperature.
  struct sheath_ctx *app = ctx;

  double Ti = app->Ti;

  fout[0] = Ti;
}

void
evalTempNeutInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  // Neutral temperature.
  struct sheath_ctx *app = ctx;
  double Tn = app->Tn;
  fout[0] = Tn/6.0;
}

static inline void
mapc2p(double t, const double* GKYL_RESTRICT zc, double* GKYL_RESTRICT xp, void* ctx)
{
  // Set physical coordinates (X, Y, Z) from computational coordinates (x, y, z).
  double x = zc[0], y = zc[1], z = zc[2];

  xp[0] = x; xp[1] = y; xp[2] = z;
}

void
bfield_func(double t, const double* GKYL_RESTRICT zc, double* GKYL_RESTRICT fout, void *ctx)
{
  struct sheath_ctx *app = ctx;

  double B0 = app->B0;

  // Set Cartesian components of magnetic field.
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = app->B0;
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

  struct sheath_ctx ctx = create_ctx(); // Context for initialization functions.

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  int my_rank = 0;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    gkyl_comm_get_rank(comm, &my_rank);
#endif

  struct gkyl_gyrokinetic_emission_inp neut_bc = {
    .num_species = 1,
    .in_species = { "ion" },
    .recycling_frac = ctx.rec_frac,
    .emission_temp = ctx.Tn_recycle,
  };

  // Electron species.
  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .vdim = ctx.vdim,
    .lower = { -ctx.vpar_max_elc, 0.0 },
    .upper = {  ctx.vpar_max_elc, ctx.mu_max_elc },
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    .is_static = true,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .density = evalDensityInit,
      .ctx_density = &ctx,
      .upar = evalUparInit,
      .ctx_upar = &ctx,
      .temp = evalTempElcInit,
      .ctx_temp = &ctx,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
    },

    .num_diag_moments = 4,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP },
  };

  // Ion species.
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.charge_ion, .mass = ctx.mass_ion,
    .vdim = ctx.vdim,
    .lower = { -ctx.vpar_max_ion, 0.0 },
    .upper = {  ctx.vpar_max_ion, ctx.mu_max_ion },
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    .is_static = true,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .density = evalDensityInit,
      .ctx_density = &ctx,
      .upar = evalUparInit,
      .ctx_upar = &ctx,
      .temp = evalTempIonInit,
      .ctx_temp = &ctx,      
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
    },

    .num_diag_moments = 4,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP },
  };

  struct gkyl_gyrokinetic_neut_species neut = {
    .name = "neut", .mass = ctx.mass_neut,
    .vdim = ctx.vdim+1,
    .lower = { -ctx.v_max_neut, -ctx.v_max_neut, -ctx.v_max_neut },
    .upper = {  ctx.v_max_neut,  ctx.v_max_neut,  ctx.v_max_neut },
    .cells = { ctx.Nv, ctx.Nv, ctx.Nv, },

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = evalDensityNeutInit,
      .ctx_upar = &ctx,
      .udrift= evalUdriftInit,
      .ctx_temp = &ctx,
      .temp = evalTempNeutInit,      
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_NEUTRAL,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
        .density = evalSourceDensityInit,
        .ctx_density = &ctx,
        .udrift = evalSourceUdriftInit,
        .ctx_upar = &ctx,
        .temp = evalSourceTempInit,
        .ctx_temp = &ctx,
      },
    },

    .react_neut = {
      .num_react = 2,
      .react_type = {
        { .react_id = GKYL_REACT_IZ,
          .type_self = GKYL_SELF_DONOR,
          .ion_id = GKYL_ION_H,
          .elc_nm = "elc",
          .ion_nm = "ion",
          .donor_nm = "neut",
          .charge_state = 0,
          .ion_mass = ctx.mass_ion,
          .elc_mass = ctx.mass_elc,
        },
        { .react_id = GKYL_REACT_CX,
          .type_self = GKYL_SELF_PARTNER,
          .ion_id = GKYL_ION_H,
          .ion_nm = "ion",
          .partner_nm = "neut",
          .ion_mass = ctx.mass_ion,
          .partner_mass = ctx.mass_ion,
        },
      },
    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_RECYCLE, .emission = neut_bc, .write_diagnostics = true, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_RECYCLE, .emission = neut_bc, .write_diagnostics = true, },
    },
    
    .num_diag_moments = 3,
    .diag_moments = { GKYL_F_MOMENT_M1_FROM_H, GKYL_F_MOMENT_ENERGY, GKYL_F_MOMENT_LTE},
  };

  // Field.
  struct gkyl_gyrokinetic_field field = {
    .kperpSq = ctx.k_perp * ctx.k_perp,
    .is_static = true,
    .zero_init_field = true,
  };

  // Gyrokinetic app.
  struct gkyl_gk app_inp = {

    .cdim = ctx.cdim,
    .lower = { -0.5 * ctx.Lz},
    .upper = {  0.5 * ctx.Lz},
    .cells = { cells_x[0] },

    .poly_order = 1,
    .basis_type = app_args.basis_type,
    .cfl_frac = 1.0,

    .geometry = {
      .geometry_id = GKYL_GEOMETRY_MAPC2P,
      .world = {0.0, 0.0},
      .mapc2p = mapc2p,
      .c2p_ctx = &ctx,
      .bfield_func = bfield_func,
      .bfield_ctx = &ctx
    },

    .num_periodic_dir = 0,
    .periodic_dirs = { },

    .num_species = 2,
    .species = { elc, ion },
    .num_neut_species = 1,
    .neut_species = { neut },

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
