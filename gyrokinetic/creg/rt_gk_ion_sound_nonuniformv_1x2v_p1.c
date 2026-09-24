#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_gyrokinetic_run.h>
#include <gkyl_util.h>

#include <rt_arg_parse.h>

struct ion_sound_ctx
{ 
  int cdim, vdim; // Dimensionality.

  // Physical constants (using non-normalized physical units).
  double mass_ion; // Proton mass.
  double charge_ion; // Proton charge.
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.

  double Te; // Electron temperature.
  double Ti; // Ion temperature.
  double n0; // Reference number density.
  double B0; // Reference magnetic field strength.

  double k_perp; // Perpendicular wave number (for Poisson solver).

  double alpha; // Applied perturbation amplitude.
  double kz; // Perturbed wave number (z-direction).

  // Derived physical quantities (using non-normalized physical units).
  double nu_ion; // Ion collision frequency.
  double nu_elc; // Electron collision frequency.

  double vte; // Electron thermal velocity.
  double vti; // Ion thermal velocity.

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

struct ion_sound_ctx
create_ctx(void)
{
  int cdim = 1, vdim = 2; // Dimensionality.

  // Physical constants (using normalized code units).
  double mass_ion = 1.0; // Ion mass.
  double charge_ion = 1.0; // Ion charge.
  double mass_elc = 1.0 / 1836.16; // Electron mass.
  double charge_elc = -1.0; // Electron charge.

  double Te = 1.0; // Electron temperature.
  double Ti = 1.0; // Ion temperature.
  double n0 = 1.0; // Reference number density.
  double B0 = 1.0; // Reference magnetic field strength.

  double k_perp = 0.1; // Perpendicular wave number (for Poisson solver).

  double alpha = 0.01; // Applied perturbation amplitude.
  double kz = 0.5; // Perturbed wave number (z-direction).

  // Derived physical quantities (using non-normalized physical units).
  double nu_ion = 2.0; // Ion collision frequency.
  double nu_elc = nu_ion * sqrt(mass_ion / mass_elc); // Electron collision frequency.
  
  double vte = sqrt(Te / mass_elc); // Electron thermal velocity.
  double vti = sqrt(Ti / mass_ion); // Ion thermal velocity.

  // Simulation parameters.
  int Nz = 8; // Cell count (configuration space: z-direction).
  int Nvpar = 64; // Cell count (velocity space: parallel velocity direction).
  int Nmu = 12; // Cell count (velocity space: magnetic moment direction).
  double Lz = 2.0 * M_PI / kz; // Domain size (configuration space: z-direction).
  double vpar_max_elc = 6.0 * vte; // Domain boundary (electron velocity space: parallel velocity direction).
  double mu_max_elc = mass_elc * pow(6.0 * vte, 2.0) / (2.0 * B0); // Domain boundary (electron velocity space: magnetic moment direction).
  double vpar_max_ion = 6.0 * vti; // Domain boundary (ion velocity space: parallel velocity direction).
  double mu_max_ion = mass_ion * pow(6.0 * vti, 2.0) / (2.0 * B0); // Domain boundary (ion velocity space: magnetic moment direction).
  int poly_order = 1; // Polynomial order.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 0.25; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  double write_phase_freq = 0.2; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.
  
  struct ion_sound_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .mass_ion = mass_ion,
    .charge_ion = charge_ion,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .Te = Te,
    .Ti = Ti,
    .n0 = n0,
    .B0 = B0,
    .k_perp = k_perp,
    .alpha = alpha,
    .kz = kz,
    .nu_ion = nu_ion,
    .nu_elc = nu_elc,
    .vte = vte,
    .vti = vti,
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
  struct ion_sound_ctx *app = ctx;

  double n0 = app->n0;

  // Set electron total number density.
  fout[0] = n0;
}

void
evalElcTempInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct ion_sound_ctx *app = ctx;

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
evalIonDensityInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct ion_sound_ctx *app = ctx;
  double z = xn[0];

  double n0 = app->n0;
  double alpha = app->alpha;
  double kz = app->kz;

  double n = (1.0 + alpha * cos(kz * z)) * n0; // Ion total number density.

  // Set ion total number density.
  fout[0] = n;
}

void
evalIonTempInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct ion_sound_ctx *app = ctx;

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
evalElcNu(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct ion_sound_ctx *app = ctx;

  double nu_elc = app->nu_elc;

  // Set electron collision frequency.
  fout[0] = nu_elc;
}

void
evalIonNu(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct ion_sound_ctx *app = ctx;

  double nu_ion = app->nu_ion;

  // Set ion collision frequency.
  fout[0] = nu_ion;
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
  struct ion_sound_ctx *app = ctx;

  double B0 = app->B0;

  // zc are computational coords. 
  // Set Cartesian components of magnetic field.
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = app->B0;
}

static inline void
mapc2p_vel_elc(double t, const double* GKYL_RESTRICT vc, double* GKYL_RESTRICT vp, void* ctx)
{
  struct ion_sound_ctx *app = ctx;
  double cvpar = vc[0], cmu = vc[1];

  double vpar_max_elc = app->vpar_max_elc;
  double mu_max_elc = app->mu_max_elc;

  double vpar = 0.0;
  double mu = 0.0;

  if (cvpar < 0.0) {
    vpar = -vpar_max_elc * (cvpar * cvpar);
  }
  else {
    vpar = vpar_max_elc * (cvpar * cvpar);
  }
  mu = mu_max_elc * (cmu * cmu);

  // Set rescaled electron velocity space coordinates (vpar, mu) from old velocity space coordinates (cvpar, cmu):
  vp[0] = vpar; vp[1] = mu;
}

static inline void
mapc2p_vel_ion(double t, const double* GKYL_RESTRICT vc, double* GKYL_RESTRICT vp, void* ctx)
{
  struct ion_sound_ctx *app = ctx;
  double cvpar = vc[0], cmu = vc[1];

  double vpar_max_ion = app->vpar_max_ion;
  double mu_max_ion = app->mu_max_ion;

  double vpar = 0.0;
  double mu = 0.0;

  if (cvpar < 0.0) {
    vpar = -vpar_max_ion * (cvpar * cvpar);
  }
  else {
    vpar = vpar_max_ion * (cvpar * cvpar);
  }
  mu = mu_max_ion * (cmu * cmu);

  // Set rescaled ion velocity space coordinates (vpar, mu) from old velocity space coordinates (cvpar, cmu):
  vp[0] = vpar; vp[1] = mu;
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

  struct ion_sound_ctx ctx = create_ctx(); // Context for init functions.

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  // Electrons.
  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .vdim = ctx.vdim,
    .lower = { -1.0, 0.0 },
    .upper = { 1.0, 1.0 },
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    .mapc2p = {
      .mapping = mapc2p_vel_elc,
      .ctx = &ctx,
    },

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .density = evalElcDensityInit,
      .ctx_density = &ctx,
      .temp = evalElcTempInit,
      .ctx_temp = &ctx,
      .upar = evalElcUparInit,
      .ctx_upar = &ctx,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
    },

    .num_diag_moments = 5,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP },
  };

  // Ions.
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.charge_ion, .mass = ctx.mass_ion,
    .vdim = ctx.vdim,
    .lower = { -1.0, 0.0 },
    .upper = { 1.0, 1.0 },
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    .mapc2p = {
      .mapping = mapc2p_vel_ion,
      .ctx = &ctx,
    },

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .density = evalIonDensityInit,
      .ctx_density = &ctx,
      .temp = evalIonTempInit,
      .ctx_temp = &ctx,
      .upar = evalIonUparInit,
      .ctx_upar = &ctx,
    },
    
    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
    },

    .num_diag_moments = 5,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP },
  };

  // Field.
  struct gkyl_gyrokinetic_field field = {
    .kperpSq = ctx.k_perp * ctx.k_perp,
  };

  // Gyrokinetic app.
  struct gkyl_gk app_inp = {

    .cdim = ctx.cdim,
    .lower = { -0.5 * ctx.Lz },
    .upper = { 0.5 * ctx.Lz },
    .cells = { cells_x[0] },

    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,
    .cfl_frac = ctx.cfl_frac,

    .geometry = {
      .geometry_id = GKYL_GEOMETRY_MAPC2P,
      .world = { 0.0, 0.0 },

      .mapc2p = mapc2p,
      .c2p_ctx = &ctx,
      .bfield_func = bfield_func,
      .bfield_ctx = &ctx
    },

    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },

    .num_species = 2,
    .species = { elc, ion },

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
