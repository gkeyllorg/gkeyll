#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_gyrokinetic_run.h>
#include <gkyl_util.h>

#include <rt_arg_parse.h>

struct lbo_relax_ctx
{
  int cdim, vdim; // Dimensionality.

  // Physical constants (using normalized code units).
  double mass; // Top hat/bump mass.
  double charge; // Top hat/bump charge.

  double B0; // Reference magnetic field strength.
  double n0; // Reference number density.
  double u0; // Reference velocity.
  double vt; // Top hat Maxwellian thermal velocity.
  double nu; // Collision frequency.

  double ab; // Bump Maxwellian amplitude.
  double sb; // Bump Maxwellian softening factor, to avoid divergence.
  double vtb; // Bump Maxwellian thermal velocity.

  // Derived physical quantities (using normalized code units).
  double ub; // Bump location (in velocity space).

  // Simulation parameters.
  int Nz; // Cell count (configuration space: z-direction).
  int Nvpar; // Cell count (velocity space: parallel velocity direction).
  int Nmu; // Cell count (velocity space: magnetic moment direction).
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  double Lz; // Domain size (configuration space: z-direction).
  double vpar_max; // Domain boundary (velocity space: parallel velocity direction).
  double mu_max; // Domain boundary (velocity space: magnetic moment direction).
  int poly_order; // Polynomial order.
  double cfl_frac; // CFL coefficient.

  double t_end; // End time.
  int num_frames; // Number of output frames.
  double write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct lbo_relax_ctx
create_ctx(void)
{
  int cdim = 1, vdim = 2; // Dimensionality.

  // Physical constants (using normalized code units).
  double mass = 1.0; // Top hat/bump mass.
  double charge = 1.0; // Top hat/bump charge.

  double B0 = 1.0; // Reference magnetic field strength.
  double n0 = 1.0; // Reference number density.
  double u0 = 0.0; // Reference velocity.
  double vt = 1.0 / 3.0; // Top hat Maxwellian thermal velocity.
  double nu = 0.01; // Collision frequency.

  double ab = sqrt(0.1); // Bump Maxwellian amplitude.
  double sb = 0.12; // Bump Maxwellian softening factor, to avoid divergence.
  double vtb = 1.0; // Bump Maxwellian thermal velocity.

  // Derived physical quantities (using normalized code units).
  double ub = 4.0 * sqrt((pow(3.0 * vt / 2.0, 2.0)) / 3.0); // Bump location (in velocity space).

  // Simulation parameters.
  int Nz = 4; // Cell count (configuration space: z-direction).
  int Nvpar = 32; // Cell count (velocity space: parallel velocity direction).
  int Nmu = 16; // Cell count (velocity space: magnetic moment direction).
  double Lz = 1.0; // Domain size (configuration space: x-direction).
  double vpar_max = 8.0 * vt; // Domain boundary (velocity space: parallel velocity direction).
  double mu_max = 0.5 * mass * pow(3.5 * vt, 2.0) / B0; // Domain boundary (velocity space: magnetic moment direction).
  int poly_order = 1; // Polynomial order.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 0.01 / nu; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  double write_phase_freq = 0.2; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct lbo_relax_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .mass = mass,
    .charge = charge,
    .B0 = B0,
    .n0 = n0,
    .u0 = u0,
    .vt = vt,
    .nu = nu,
    .ab = ab,
    .sb = sb,
    .vtb = vtb,
    .ub = ub,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nz, Nvpar, Nmu},
    .Lz = Lz,
    .vpar_max = vpar_max,
    .mu_max = mu_max,
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
evalTopHatInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct lbo_relax_ctx *app = ctx;
  double vpar = xn[1];

  double n0 = app->n0;
  double vt = app->vt;

  double v0 = sqrt(3.0) * vt;

  double n = 0.0;

  if (fabs(vpar) < v0) {
    n = n0 / 2.0 / v0; // Distribution function (low velocity).
  }
  else {
    n = 0.0; // Distribution function (high velocity).
  }

  // Set distribution function.
  fout[0] = n;
}

void
evalBumpInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct lbo_relax_ctx *app = ctx;
  double vpar = xn[1], mu = xn[2];

  double B0 = app->B0;
  double n0 = app->n0;
  double u0 = app->u0;
  double vt = app->vt;

  double ab = app->ab;
  double sb = app->sb;
  double vtb = app->vtb;

  double ub = app->ub;

  double v_sq = ((vpar - u0) / (sqrt(2.0) * vt)) * ((vpar - u0) / (sqrt(2.0) * vt)) + mu * B0;
  double vb_sq = ((vpar - u0) / (sqrt(2.0) * vtb)) * ((vpar - u0) / (sqrt(2.0) * vtb)) + mu * B0;
  
  double n = (n0 / sqrt(2.0 * M_PI * vt)) * exp(-v_sq) + (n0 / sqrt(2.0 * M_PI * vtb)) *
    exp(-vb_sq) * (ab * ab) / ((vpar - ub) * (vpar - ub) + sb * sb); // Distribution function.

  // Set distribution function.
  fout[0] = n;
}

void
evalTopHatNu(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct lbo_relax_ctx *app = ctx;

  double nu = app->nu;

  // Set collision frequency.
  fout[0] = nu;
}

void
evalBumpNu(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct lbo_relax_ctx *app = ctx;

  double nu = app->nu;

  // Set collision frequency.
  fout[0] = nu;
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
  struct lbo_relax_ctx *app = ctx;
  
  double B0 = app->B0;

  // zc are computational coords. 
  // Set Cartesian components of magnetic field.
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = app->B0;
}

static inline void
mapc2p_vel_square(double t, const double* GKYL_RESTRICT vc, double* GKYL_RESTRICT vp, void* ctx)
{
  struct lbo_relax_ctx *app = ctx;
  double cvpar = vc[0], cmu = vc[1];

  double vpar_max = app->vpar_max;
  double mu_max = app->mu_max;

  double vpar = 0.0;
  double mu = 0.0;

  if (cvpar < 0.0) {
    vpar = -vpar_max * (cvpar * cvpar);
  }
  else {
    vpar = vpar_max * (cvpar * cvpar);
  }
  mu = mu_max * (cmu * cmu);

  // Set rescaled top hat velocity space coordinates (vpar, mu) from old velocity space coordinates (cvpar, cmu):
  vp[0] = vpar; vp[1] = mu;
}

static inline void
mapc2p_vel_bump(double t, const double* GKYL_RESTRICT vc, double* GKYL_RESTRICT vp, void* ctx)
{
  struct lbo_relax_ctx *app = ctx;
  double cvpar = vc[0], cmu = vc[1];

  double vpar_max = app->vpar_max;
  double mu_max = app->mu_max;

  double vpar = 0.0;
  double mu = 0.0;

  if (cvpar < 0.0) {
    vpar = -vpar_max * (cvpar * cvpar);
  }
  else {
    vpar = vpar_max * (cvpar * cvpar);
  }
  mu = mu_max * (cmu * cmu);

  // Set rescaled bump velocity space coordinates (vpar, mu) from old velocity space coordinates (cvpar, cmu):
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

  struct lbo_relax_ctx ctx = create_ctx(); // Context for init functions.

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  // Top hat species.
  struct gkyl_gyrokinetic_species square = {
    .name = "square",
    .charge = ctx.charge, .mass = ctx.mass,
    .vdim = ctx.vdim,
    .lower = { -1.0, 0.0 },
    .upper = { 1.0, 1.0 },
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    .mapc2p = {
      .mapping = mapc2p_vel_square,
      .ctx = &ctx,
    },

    .projection = {
      .proj_id = GKYL_PROJ_FUNC,
      .func = evalTopHatInit,
      .ctx_func = &ctx,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .self_nu = evalTopHatNu,
      .self_nu_ctx = &ctx,
    },
    
    .num_diag_moments = 7,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP },
  };

  // Bump species.
  struct gkyl_gyrokinetic_species bump = {
    .name = "bump",
    .charge = ctx.charge, .mass = ctx.mass,
    .vdim = ctx.vdim,
    .lower = { -1.0, 0.0 },
    .upper = { 1.0, 1.0 },
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    .mapc2p = {
      .mapping = mapc2p_vel_bump,
      .ctx = &ctx,
    },

    .projection = {
      .proj_id = GKYL_PROJ_FUNC,
      .func = evalBumpInit,
      .ctx_func = &ctx,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .self_nu = evalBumpNu,
      .self_nu_ctx = &ctx,
    },
    
    .num_diag_moments = 7,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP },
  };

  // Field.
  struct gkyl_gyrokinetic_field field = {
    .gkfield_id = GKYL_GK_FIELD_BOLTZMANN,

    .electron_mass = ctx.mass,
    .electron_charge = ctx.charge,
    .electron_temp = ctx.vt * ctx.vt * ctx.mass,

    .zero_init_field = true, // Don't compute the field at t = 0.
    .is_static = true, // Don't evolve the field in time.
  };

  // Gyrokinetic app.
  struct gkyl_gk app_inp = {

    .cdim = ctx.cdim,
    .lower = { 0.0 },
    .upper = { ctx.Lz },
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
    .species = { square, bump },

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
