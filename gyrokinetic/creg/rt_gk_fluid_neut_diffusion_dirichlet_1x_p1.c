#include <math.h>
#include <stdio.h>

#include <gkyl_gyrokinetic.h>
#include <gkyl_gyrokinetic_run.h>
#include <rt_arg_parse.h>

struct fluid_neut_diffusion_ctx {
  double zlo, zup; // Configuration-space boundaries.
  double nlo, nup; // Dirichlet neutral number densities.
  double amplitude; // Initial sine-mode amplitude.
  double temperature; // Neutral temperature.
  double diffusion; // Constant neutral diffusion coefficient.
};

static void
eval_neut_density(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *vctx)
{
  struct fluid_neut_diffusion_ctx *ctx = vctx;
  const double L = ctx->zup - ctx->zlo;
  const double s = (xn[0] - ctx->zlo) / L;
  const double nss = ctx->nlo + (ctx->nup - ctx->nlo) * s;
  fout[0] = nss + ctx->amplitude * exp(-ctx->diffusion * pow(M_PI / L, 2) * t) * sin(M_PI * s);
}

static void
eval_neut_udrift(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = fout[1] = fout[2] = 0.0;
}

static void
eval_neut_temp(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *vctx)
{
  struct fluid_neut_diffusion_ctx *ctx = vctx;
  fout[0] = ctx->temperature;
}

static void
mapc2p(double t, const double *GKYL_RESTRICT xc, double *GKYL_RESTRICT xp, void *ctx)
{
  // Set physical coordinates (X,Y,Z) from computational coordinates.
  xp[0] = xc[0];
  xp[1] = xc[1];
  xp[2] = xc[2];
}

static void
bfield_func(double t, const double *GKYL_RESTRICT xc, double *GKYL_RESTRICT fout, void *ctx)
{
  // Set Cartesian components of the magnetic field.
  fout[0] = 0.0;
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
    .zlo = -0.5,
    .zup = 0.5,
    .nlo = 1.0,
    .nup = 1.7,
    .amplitude = 0.3,
    .temperature = 1.0,
    .diffusion = 0.4,
  };
  const int cells = APP_ARGS_CHOOSE(app_args.xcells[0], 32);
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  // Fluid-neutral species.
  struct gkyl_gyrokinetic_neut_species neut = {
    .name = "n",
    .mass = 2.0,
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
        .D = ctx.diffusion,
        .is_implicit = true,
        .implicit_tol = 1.0e-12,
        .implicit_max_iter = 200,
        .write_diagnostics = true,
        .lower_bc_type = GKYL_NEUT_FLUID_DIFFUSION_DIRICHLET,
        .upper_bc_type = GKYL_NEUT_FLUID_DIFFUSION_DIRICHLET,
        .lower_bc_density = ctx.nlo,
        .upper_bc_density = ctx.nup,
      },
    .num_diag_moments = 1,
    .diag_moments = {GKYL_F_MOMENT_M0},
  };

  // Gyrokinetic app input.
  struct gkyl_gk app_inp = {
    .cdim = 1,
    .lower = {ctx.zlo},
    .upper = {ctx.zup},
    .cells = {cells},
    .poly_order = 1,
    .basis_type = app_args.basis_type,
    .cfl_frac = 0.8,
    .geometry = {.geometry_id = GKYL_GEOMETRY_MAPC2P, .mapc2p = mapc2p, .bfield_func = bfield_func},
    .num_species = 0,
    .num_neut_species = 1,
    .neut_species = {neut},
    .field = {.is_static = true, .zero_init_field = true},
    .parallelism = {.use_gpu = app_args.use_gpu, .cuts = {app_args.cuts[0]}, .comm = comm},
  };

  snprintf(app_inp.name, sizeof(app_inp.name), "%s", app_args.app_name);
  struct gkyl_gyrokinetic_run_inp run_inp = {
    .app_inp = app_inp,
    .time_stepping =
      {
        .t_end = 0.5,
        .num_frames = 5,
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
