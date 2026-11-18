#include <math.h>
#include <stdio.h>

#include <gkyl_gyrokinetic.h>
#include <gkyl_gyrokinetic_run.h>

#include <rt_arg_parse.h>

struct fluid_neut_diffusion_ctx {
  double n0; // Reference neutral number density.
  double amplitude; // Relative density-perturbation amplitude.
  double temperature; // Neutral temperature.
  double diffusion_coeff; // Constant neutral diffusion coefficient.
  double map_amp; // Amplitude of the nonuniform-x map.
  double shear; // Shear of X with respect to z.
};

static void
eval_neut_density(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *vctx)
{
  struct fluid_neut_diffusion_ctx *ctx = vctx;
  const double z = xn[2];
  double X = xn[0] + ctx->map_amp * sin(xn[0]) + ctx->shear * sin(z);
  fout[0] = ctx->n0 * (1.0 + ctx->amplitude * cos(X) * cos(z));
}

static void
eval_neut_udrift(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *vctx)
{
  fout[0] = 0.2;
  fout[1] = 0.0;
  fout[2] = 0.0;
}

static void
eval_neut_temp(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *vctx)
{
  struct fluid_neut_diffusion_ctx *ctx = vctx;
  fout[0] = ctx->temperature;
}

static void
mapc2p(double t, const double *GKYL_RESTRICT xc, double *GKYL_RESTRICT xp, void *vctx)
{
  struct fluid_neut_diffusion_ctx *ctx = vctx;
  // Set physical coordinates (X,Y,Z) from computational coordinates.
  xp[0] = xc[0] + ctx->map_amp * sin(xc[0]) + ctx->shear * sin(xc[2]);
  xp[1] = xc[1];
  xp[2] = xc[2];
}

static void
bfield_func(double t, const double *GKYL_RESTRICT xc, double *GKYL_RESTRICT fout, void *vctx)
{
  struct fluid_neut_diffusion_ctx *ctx = vctx;
  // Set a field tangent to the sheared z coordinate lines.
  fout[0] = ctx->shear * cos(xc[2]);
  fout[1] = 0.0;
  fout[2] = 1.0;
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

  struct fluid_neut_diffusion_ctx ctx = {
    .n0 = 1.0,
    .amplitude = 0.1,
    .temperature = 1.0,
    .diffusion_coeff = 0.1,
    .map_amp = 0.2,
    .shear = 0.25,
  };

  int cells_x = APP_ARGS_CHOOSE(app_args.xcells[0], 24);
  int cells_y = APP_ARGS_CHOOSE(app_args.xcells[1], 8);
  int cells_z = APP_ARGS_CHOOSE(app_args.xcells[2], 24);
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  // Fluid-neutral species.
  struct gkyl_gyrokinetic_neut_species neut = {
    .name = "n",
    .mass = 1.0,
    .vdim = 0,
    .gas_gamma = 5.0 / 3.0,

    .projection =
      {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
        .density = eval_neut_density,
        .ctx_density = &ctx,
        .udrift = eval_neut_udrift,
        .ctx_udrift = &ctx,
        .temp = eval_neut_temp,
        .ctx_temp = &ctx,
      },

    // Configuration-space density diffusion.
    .diffusion =
      {
        .D = ctx.diffusion_coeff,
        .is_implicit = true,
        .implicit_tol = 1.0e-12,
        .implicit_max_iter = 200,
        .write_diagnostics = true,
      },

    .num_diag_moments = 3,
    .diag_moments = {GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2},
  };

  // Gyrokinetic app input.
  struct gkyl_gk app_inp = {
    .cdim = 3,
    .lower = {-M_PI, -M_PI, -M_PI},
    .upper = {M_PI, M_PI, M_PI},
    .cells = {cells_x, cells_y, cells_z},
    .poly_order = 1,
    .basis_type = app_args.basis_type,
    .cfl_frac = 0.9,

    .geometry =
      {
        .geometry_id = GKYL_GEOMETRY_MAPC2P,
        .mapc2p = mapc2p,
        .c2p_ctx = &ctx,
        .bfield_func = bfield_func,
        .bfield_ctx = &ctx,
      },

    .num_periodic_dir = 3,
    .periodic_dirs = {0, 1, 2},

    .num_species = 0,
    .num_neut_species = 1,
    .neut_species = {neut},

    .field = {.is_static = true, .zero_init_field = true},

    .parallelism =
      {
        .use_gpu = app_args.use_gpu,
        .cuts = {app_args.cuts[0], app_args.cuts[1], app_args.cuts[2]},
        .comm = comm,
      },
  };

  snprintf(app_inp.name, sizeof(app_inp.name), "%s", app_args.app_name);
  struct gkyl_gyrokinetic_run_inp run_inp = {
    .app_inp = app_inp,
    .time_stepping =
      {
        .t_end = 0.2,
        .num_frames = 1,
        .write_phase_freq = 1.0,
        .int_diag_calc_num = 10,
        .dt_failure_tol = 1.0e-4,
        .num_failures_max = 20,
        .is_restart = app_args.is_restart,
        .restart_frame = app_args.restart_frame,
        .num_steps = app_args.num_steps,
      },
  };

  gkyl_gyrokinetic_run_simulation(&run_inp);
  gkyl_gyrokinetic_comms_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif

  return 0;
}
