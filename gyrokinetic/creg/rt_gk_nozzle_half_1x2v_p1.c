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
struct gk_nozzle_ctx
{
  int cdim, vdim; // Dimensionality.

  // Plasma parameters
  double me, mi;
  double qe, qi;
  double n_init;
  double Te_init;
  double Ti_init;
  double nu_ion;
  double B_p;
  // Thermal speeds.
  double vti;
  // Gyrofrequencies and gyroradii.
  double z_min;
  double z_max;
  double psi_eval;
  // Grid parameters
  double vpar_max_ion;
  double mu_max_ion;
  int Nz;
  int Nvpar;
  int Nmu;
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  int poly_order;

  double t_end; // End time.
  int num_frames; // Number of output frames.
  double write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

void
eval_density_ion_init(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_nozzle_ctx *app = ctx;
  double z = xn[0];
  if (fabs(z) < 0.2){
    fout[0] = app->n_init;
  } else {
    fout[0] = 1e-6 * app->n_init;
  }
}

void
eval_upar_ion_init(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
eval_temp_ion_init(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_nozzle_ctx *app = ctx;
  double z = xn[0];
  fout[0] = app->Ti_init;
}

void
eval_nu_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_nozzle_ctx *app = ctx;
  fout[0] = app->nu_ion;
}

void mapc2p_vel_ion(double t, const double *vc, double* GKYL_RESTRICT vp, void *ctx)
{
  struct gk_nozzle_ctx *app = ctx;
  double vpar_max_ion = app->vpar_max_ion;
  double mu_max_ion = app->mu_max_ion;

  double cvpar = vc[0], cmu = vc[1];
  // double b = 1.45;
  // double linear_velocity_threshold = 1./6.;
  // double frac_linear = 1/b*atan(linear_velocity_threshold*tan(b));
  // if (fabs(cvpar) < frac_linear) {
  //   double func_frac = tan(frac_linear*b) / tan(b);
  //   vp[0] = vpar_max_ion*func_frac*cvpar/frac_linear;
  // }
  // else {
  //   vp[0] = vpar_max_ion*tan(cvpar*b)/tan(b);
  // }
  // Quadratic map in mu.
  // vp[1] = mu_max_ion*pow(cmu,2);

  vp[0] = vpar_max_ion*cvpar;
  vp[1] = mu_max_ion*cmu;
}

struct gk_nozzle_ctx
create_ctx(void)
{
  int cdim = 1, vdim = 2; // Dimensionality.

  // Universal constant parameters.
  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double me = GKYL_ELECTRON_MASS; // Electron mass.
  double mi = GKYL_PROTON_MASS; // Ion mass.
  double qe = -eV;  // Electron charge.
  double qi = eV;  // Ion charge.

  double z_min =  0.0;
  double z_max =  1.0;
  double psi_eval = 1e-5;
  double B_p = 0.008;

  // Plasma parameters
  double n_init = 3e19;
  double Te_init = 1e3 * eV;
  double Ti_init = 10e3 * eV;
  double vti = sqrt(Ti_init / mi);
  printf("vti = %g\n", vti);

  // Grid parameters
  double vpar_max_ion = 6 * vti;
  double mu_max_ion = mi * pow(8. * vti, 2.) / (2. * B_p);
  printf("mu_max_ion = %g\n", mu_max_ion);
  printf("vpar_max_ion = %g\n", vpar_max_ion);
  int Nz = 16;
  int Nvpar = 32; // Number of cells in the paralell velocity direction 96
  int Nmu = 32;  // Number of cells in the mu direction 192
  int poly_order = 1;

  // double loglambda_ion = 6.6 - 0.5 * log(n_init / 1e20) + 1.5 * log(Ti_init / eV);
  // double nu_frac = 100.0;
  // double nu_ion = nu_frac * loglambda_ion * pow(eV, 4.) * n_init /
  //                (12 * pow(M_PI, 3. / 2.) * pow(eps0, 2.) * sqrt(mi) * pow(Ti_init, 3. / 2.));
  
  double nu_ion = 1 / (5e-3);
  printf("nu_ion = %g\n", nu_ion);
  printf("1/nu_ion = %g\n", 1.0/nu_ion);

  double t_end = 10e-9;
  int num_frames = 100;
  double write_phase_freq = 0.2; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.


  struct gk_nozzle_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .me = me,
    .mi = mi,
    .qe = qe,
    .qi = qi,
    .nu_ion = nu_ion,
    .B_p = B_p,
    .z_min = z_min,
    .z_max = z_max,
    .psi_eval = psi_eval,
    .vpar_max_ion = vpar_max_ion,
    .mu_max_ion = mu_max_ion,
    .n_init = n_init,
    .Te_init = Te_init,
    .Ti_init = Ti_init,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nz, Nvpar, Nmu},
    .poly_order = poly_order,
    .t_end = t_end,
    .num_frames = num_frames,
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

  struct gk_nozzle_ctx ctx = create_ctx(); // Context for init functions.

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  struct gkyl_gyrokinetic_projection proj_init = {
    .proj_id = GKYL_PROJ_BIMAXWELLIAN,
    .density = eval_density_ion_init,
    .ctx_density = &ctx,
    .upar = eval_upar_ion_init,
    .ctx_upar = &ctx,
    .temppar = eval_temp_ion_init,
    .ctx_temppar = &ctx,
    .tempperp = eval_temp_ion_init,
    .ctx_tempperp = &ctx,
    .correct_all_moms = true,
  };

  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.qi,
    .mass = ctx.mi,
    .vdim = ctx.vdim,
    .lower = {-ctx.vpar_max_ion, 0.0},
    .upper = { ctx.vpar_max_ion, ctx.mu_max_ion},
    .cells = { cells_v[0], cells_v[1] },

    .projection = proj_init,

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .self_nu = eval_nu_ion,
      .self_nu_ctx = &ctx,
    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_FIXED_FUNC, .projection = proj_init, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB, },
    },

    .num_diag_moments = 6,
    .diag_moments = {GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_BIMAXWELLIAN},
  };

  struct gkyl_mirror_geo_grid_inp grid_inp = {
    .filename_psi = "gyrokinetic/data/unit/single_coil.geqdsk_psi.gkyl", // psi file to use
    .rclose = 0.2, // closest R to region of interest
    .zmin = -1.0,  // Z of lower boundary
    .zmax =  1.0,  // Z of upper boundary 
    .include_axis = false, // Include R=0 axis in grid
    .fl_coord = GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z, // coordinate system for psi grid
  };

  struct gkyl_gyrokinetic_field field = {
    .gkfield_id = GKYL_GK_FIELD_BOLTZMANN,

    .electron_mass = ctx.me,
    .electron_charge = ctx.qe,
    .electron_temp = ctx.Te_init,
    .polarization_bmag = ctx.B_p, // Issue here. B0 from soloviev, so not sure what to do. Ours is not constant

    .zero_init_field = true, // Don't compute the field at t=0.
    .is_static = true, // Don't update the field in time.
  };

  // GK app
  struct gkyl_gk app_inp = {
    .cdim = ctx.cdim,
    .lower = {ctx.z_min},
    .upper = {ctx.z_max},
    .cells = { cells_x[0] },
    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,

    .geometry = {
      .geometry_id = GKYL_GEOMETRY_MIRROR,
      .world = {ctx.psi_eval, 0.0},
      .mirror_grid_info = grid_inp,
    },

    .field = field,

    .num_periodic_dir = 0,
    .periodic_dirs = {},

    .num_species = 1,
    .species = {ion},

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
