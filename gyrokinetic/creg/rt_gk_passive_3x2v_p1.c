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

struct passive_ctx
{
  int cdim, vdim; // Dimensionality.

  // Physical constants (using non-normalized physical units).
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.

  double Te; // Electron temperature.
  double n0; // Reference number density (1 / m^3).
  double B0; // Reference magnetic field strength (Tesla).
  
  double ux, uy, uz; // Passive advection velocity.
  double f_amplitude; // Amplitude of the distribution.
  double f_floor; // Floor of the distribution.

  // Simulation parameters.
  int Nx; // Cell count (configuration space: x-direction).
  int Ny; // Cell count (configuration space: y-direction).
  int Nz; // Cell count (configuration space: z-direction).
  int Nvpar; // Cell count (velocity space: parallel velocity direction).
  int Nmu; // Cell count (velocity space: magnetic moment direction).
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  double Lx; // Domain size (configuration space: x-direction).
  double Ly; // Domain size (configuration space: y-direction).
  double Lz; // Domain size (configuration space: z-direction).
  double vpar_max_elc; // Domain boundary (electron velocity space: parallel velocity direction).
  double mu_max_elc; // Domain boundary (electron velocity space: magnetic moment direction).
  int poly_order; // Polynomial order.
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  double write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct passive_ctx
create_ctx(void)
{
  int cdim = 3, vdim = 2; // Dimensionality.

  // Physical constants (using non-normalized physical units).
  double mass_elc = GKYL_ELECTRON_MASS; // Electron mass.
  double charge_elc = -GKYL_ELEMENTARY_CHARGE; // Electron charge.

  double Te = 1.0 * GKYL_ELEMENTARY_CHARGE; // Electron temperature.
  double n0 = 1.0e18; //  Reference number density (1 / m^3).
  double B0 = 1.0; // Reference magnetic field.

  double ux = 0.0, uy = 0.0, uz = 1.0; // Passive advection velocity.
  double f_amplitude = 1.0; // Amplitude of the distribution.
  double f_floor = 1.0e-10; // Floor of the distribution.

  double vte = sqrt(Te/mass_elc); // Electron thermal velocity.

  // Simulation parameters.
  int Nx = 16; // Cell count (configuration space: x-direction).
  int Ny = 16; // Cell count (configuration space: y-direction).
  int Nz = 16; // Cell count (configuration space: z-direction).
  int Nvpar = 2; // Cell count (velocity space: parallel velocity direction).
  int Nmu = 2; // Cell count (velocity space: magnetic moment direction).
  double Lx = 1.0; // Domain size (configuration space: x-direction).
  double Ly = 1.0; // Domain size (configuration space: y-direction).
  double Lz = 1.0; // Domain size (configuration space: z-direction).
  double vpar_max_elc = 4.0*vte; // Domain boundary (electron velocity space: parallel velocity direction).
  double mu_max_elc = mass_elc*pow(vpar_max_elc,2.0)/(2.0*B0); // Domain boundary (electron velocity space: magnetic moment direction).
  int poly_order = 1; // Polynomial order.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 1.; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  double write_phase_freq = 1.0; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.
  
  struct passive_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .Te = Te,
    .n0 = n0,
    .B0 = B0,
    .ux = ux, .uy = uy, uz = uz,
    .f_amplitude = f_amplitude,
    .f_floor = f_floor,
    .Nx = Nx,
    .Ny = Ny,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nx, Ny, Nz, Nvpar, Nmu},
    .Lx = Lx,
    .Ly = Ly,
    .Lz = Lz,
    .vpar_max_elc = vpar_max_elc,
    .mu_max_elc = mu_max_elc,
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
distf_elc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1], z = xn[2], vpar = xn[3], mu = xn[4];

  struct passive_ctx *app = ctx;
  double Lx = app->Lx;
  double Ly = app->Ly;
  double Lz = app->Lz;
  double f_amplitude = app->f_amplitude;
  double f_floor = app->f_floor;

  // Cube
  double rx2 = pow(x-Lx/2,2);
  double ry2 = pow(y-Ly/2,2);
  double rz2 = pow(z-Lz/2,2);

  if (rx2 < pow(Lx/4,2) && ry2 < pow(Ly/4,2) && rz2 < pow(Lz/4,2))
    fout[0] = f_amplitude;
  else
    fout[0] = f_floor;
}

void
passive_velocity_elc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];

  struct passive_ctx *app = ctx;
  double ux = app->ux;
  double uy = app->uy;
  double uz = app->uz;

  fout[0] = ux;
  fout[1] = uy;
  fout[2] = uz;
}

static inline void
mapc2p(double t, const double* GKYL_RESTRICT zc, double* GKYL_RESTRICT xp, void* ctx)
{
  double x = zc[0], y = zc[1], z = zc[2];

  struct passive_ctx *app = ctx;

  xp[0] = x; xp[1] = y; xp[2] = z;
}

void
bfield_func(double t, const double* GKYL_RESTRICT zc, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = zc[0], y = zc[1], z = zc[2];

  struct passive_ctx *app = ctx;
  double B0 = app->B0;

  // zc are computational coords. 
  // Set Cartesian components of magnetic field.
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = B0;
}

void bc_shift_func_lo(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xc[0], y = xc[1], z = xc[2];

  struct gk_app_ctx *app = ctx;

  fout[0] = -(x-0.5);
}

void bc_shift_func_up(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  bc_shift_func_lo(t,xc,fout,ctx);
  fout[0] *= -1.0;
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

  struct passive_ctx ctx = create_ctx(); // Context for init functions.

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
    .lower = { -ctx.vpar_max_elc, 0.0 },
    .upper = {  ctx.vpar_max_elc, ctx.mu_max_elc },
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_FUNC,
      .func = distf_elc,
      .ctx_func = &ctx,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_PASSIVE,
      .passive_speeds = passive_velocity_elc, 
      .passive_speeds_ctx = &ctx,
      .write_diagnostics = true,
    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_COPY, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_COPY, },
      { .dir = 2, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_TWISTSHIFT, },
      { .dir = 2, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_TWISTSHIFT, },
    },

    .num_diag_moments = 3,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, },
//    .num_integrated_diag_moments = 1,
//    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
//    .time_rate_diagnostics = true,
//
//    .boundary_flux_diagnostics = {
//      .num_diag_moments = 1,
//      .diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
//      .num_integrated_diag_moments = 1,
//      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
////      .time_integrated = true,
//    },
  };


  // Field.
  struct gkyl_gyrokinetic_field field = {
    .gkfield_id = GKYL_GK_FIELD_BOLTZMANN,
    .zero_init_field = true,
    .is_static = true,
  };

  // Gyrokinetic app.
  struct gkyl_gk app_inp = {

    .cdim = ctx.cdim,
    .lower = { 0.0   , 0.0   , 0.0    },
    .upper = { ctx.Lx, ctx.Ly, ctx.Lz },
    .cells = { cells_x[0], cells_x[1], cells_x[2] },

    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,
    .cfl_frac = ctx.cfl_frac,

    .geometry = {
      .geometry_id = GKYL_GEOMETRY_MAPC2P,
      .world = { },
      .mapc2p = mapc2p,
      .c2p_ctx = &ctx,
      .bfield_func = bfield_func,
      .bfield_ctx = &ctx,
      .parallel_lower_bc_shift_func = bc_shift_func_lo,
      .parallel_upper_bc_shift_func = bc_shift_func_up,
      .parallel_lower_bc_shift_ctx = &ctx,
      .parallel_upper_bc_shift_ctx = &ctx,
    },

    .num_periodic_dir = 1,
    .periodic_dirs = { 1, },

    .num_species = 1,
    .species = { elc, },

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
//    .print_verbosity = {
//      .enabled = true,
//      .disable_timings = true,
//    },
  };

  gkyl_gyrokinetic_run_simulation(&run_inp);
  
  gkyl_gyrokinetic_comms_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif

  return 0;
}
