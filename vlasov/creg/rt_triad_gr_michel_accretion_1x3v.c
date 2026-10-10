// Michel accretion of a Juttner gas onto a Schwarzschild black hole in ingoing Kerr-Schild
// coordinates through the GR triad Vlasov solver with implicit BGK collisions (1x3v, radial preset
// geometry with cubic velocity maps). The exact transonic Michel solution with sonic radius 8 M
// (temperature 0.037 at r = 15 M rising to 0.15 inside the horizon, radial velocity -0.14 to -0.97)
// is tabulated and projected as the initial state; the outer boundary holds the Michel state and
// the inner boundary at r = 1.8 M, inside the horizon, absorbs the inflow. The accretion rate r^2 n
// U^r computed from the rest-frame LTE moments must equal the exact -16 M^2 at every radius, and
// the profile must stay stationary over the run. Figures of merit (serendipity p2, 16 x 6 x 6 x 6
// cells on cubic velocity maps, 332 steps): the accretion rate r^2 n U^r from the rest-frame LTE
// moments is -16.007 at t = 0 (exact -16, constant in r to 3.4e-3) and -16.17 at t = 20 M
// (constant to 2.8e-2); over the run the density, temperature and drift profiles change by at most
// 1.7e-2, 2.0e-2 and 1.5e-2 inside r = 12 M; mass changes by -1.9e-2 through the horizon boundary
// as the inner cells settle.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>
#include <gkyl_util.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#ifdef GKYL_HAVE_NCCL
#include <gkyl_nccl_comm.h>
#endif
#endif

#include <rt_arg_parse.h>

struct michel_ctx {
  double mass; // Neutral mass.
  double charge; // Neutral charge.
  double mass_bh; // Black hole mass M.
  double nu; // Collision frequency.
  double vt; // Momentum scale (units of m c).

  int Nr; // Cell count (configuration space: radial direction).
  int Nvr; // Cell count (velocity space: radial momentum).
  int Nvtheta; // Cell count (velocity space: polar momentum).
  int Nvphi; // Cell count (velocity space: azimuthal momentum).
  double r_min; // Inner radius (inside the horizon; absorbing).
  double r_max; // Outer radius (fixed to the Michel state).
  double vr_max; // Domain boundary (velocity space: radial momentum).
  double vtheta_max; // Domain boundary (velocity space: polar momentum).
  double vphi_max; // Domain boundary (velocity space: azimuthal momentum).
  int poly_order; // Polynomial order.
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  int field_energy_calcs; // Number of times to calculate field energy.
  int integrated_mom_calcs; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct michel_ctx
create_ctx(void)
{
  double mass = 1.0; // Neutral mass.
  double charge = 0.0; // Neutral charge.
  double mass_bh = 1.0; // Black hole mass M.
  double nu = 15000.0; // Collision frequency.
  double vt = 1.0; // Momentum scale (units of m c).

  int Nr = 16; // Cell count (configuration space: radial direction).
  int Nvr = 6; // Cell count (velocity space: radial momentum).
  int Nvtheta = 6; // Cell count (velocity space: polar momentum).
  int Nvphi = 6; // Cell count (velocity space: azimuthal momentum).
  double r_min = 1.8; // Inner radius (inside the horizon; absorbing).
  double r_max = 15.0; // Outer radius (fixed to the Michel state).
  double vr_max = 1.5 * vt; // Domain boundary (velocity space: radial momentum).
  double vtheta_max = 1.5 * vt; // Domain boundary (velocity space: polar momentum).
  double vphi_max = 1.5 * vt; // Domain boundary (velocity space: azimuthal momentum).
  int poly_order = 2; // Polynomial order.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 20.0; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs =
    INT_MAX; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct michel_ctx ctx = {
    .mass = mass,
    .charge = charge,
    .mass_bh = mass_bh,
    .nu = nu,
    .vt = vt,
    .Nr = Nr,
    .Nvr = Nvr,
    .Nvtheta = Nvtheta,
    .Nvphi = Nvphi,
    .r_min = r_min,
    .r_max = r_max,
    .vr_max = vr_max,
    .vtheta_max = vtheta_max,
    .vphi_max = vphi_max,
    .poly_order = poly_order,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .num_frames = num_frames,
    .field_energy_calcs = field_energy_calcs,
    .integrated_mom_calcs = integrated_mom_calcs,
    .integrated_L2_f_calcs = integrated_L2_f_calcs,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };

  return ctx;
}

// Michel accretion profile for a Juttner gas onto a Schwarzschild black hole (M = 1), sonic radius
// 8 M: rest-frame density, temperature and the radial spatial four-velocity relative to the
// Kerr-Schild normal observer, tabulated on 67 radii in [1.8, 15] (see the header).
#define MICHEL_NPTS 67
static const double michel_r[MICHEL_NPTS] = {
  1.800000,  2.000000,  2.200000,  2.400000,  2.600000,  2.800000,  3.000000,  3.200000,  3.400000,
  3.600000,  3.800000,  4.000000,  4.200000,  4.400000,  4.600000,  4.800000,  5.000000,  5.200000,
  5.400000,  5.600000,  5.800000,  6.000000,  6.200000,  6.400000,  6.600000,  6.800000,  7.000000,
  7.200000,  7.400000,  7.600000,  7.800000,  8.000000,  8.200000,  8.400000,  8.600000,  8.800000,
  9.000000,  9.200000,  9.400000,  9.600000,  9.800000,  10.000000, 10.200000, 10.400000, 10.600000,
  10.800000, 11.000000, 11.200000, 11.400000, 11.600000, 11.800000, 12.000000, 12.200000, 12.400000,
  12.600000, 12.800000, 13.000000, 13.200000, 13.400000, 13.600000, 13.800000, 14.000000, 14.200000,
  14.400000, 14.600000, 14.800000, 15.000000,
};
static const double michel_n[MICHEL_NPTS] = {
  6.164246e+00, 5.373017e+00, 4.749860e+00, 4.247912e+00, 3.835976e+00, 3.492549e+00, 3.202359e+00,
  2.954288e+00, 2.740060e+00, 2.553398e+00, 2.389463e+00, 2.244466e+00, 2.115403e+00, 1.999860e+00,
  1.895884e+00, 1.801870e+00, 1.716494e+00, 1.638653e+00, 1.567421e+00, 1.502014e+00, 1.441767e+00,
  1.386110e+00, 1.334552e+00, 1.286670e+00, 1.242095e+00, 1.200505e+00, 1.161620e+00, 1.125190e+00,
  1.090997e+00, 1.058846e+00, 1.028565e+00, 1.000000e+00, 9.730129e-01, 9.474801e-01, 9.232903e-01,
  9.003430e-01, 8.785475e-01, 8.578215e-01, 8.380903e-01, 8.192859e-01, 8.013460e-01, 7.842141e-01,
  7.678381e-01, 7.521704e-01, 7.371672e-01, 7.227883e-01, 7.089963e-01, 6.957572e-01, 6.830390e-01,
  6.708124e-01, 6.590503e-01, 6.477272e-01, 6.368197e-01, 6.263058e-01, 6.161653e-01, 6.063789e-01,
  5.969291e-01, 5.877990e-01, 5.789732e-01, 5.704371e-01, 5.621770e-01, 5.541801e-01, 5.464342e-01,
  5.389280e-01, 5.316509e-01, 5.245928e-01, 5.177441e-01,
};
static const double michel_T[MICHEL_NPTS] = {
  1.512292e-01, 1.409009e-01, 1.321331e-01, 1.245834e-01, 1.180048e-01, 1.122143e-01, 1.070730e-01,
  1.024735e-01, 9.833132e-02, 9.457910e-02, 9.116228e-02, 8.803625e-02, 8.516412e-02, 8.251511e-02,
  8.006333e-02, 7.778680e-02, 7.566677e-02, 7.368711e-02, 7.183388e-02, 7.009497e-02, 6.845977e-02,
  6.691900e-02, 6.546444e-02, 6.408886e-02, 6.278578e-02, 6.154947e-02, 6.037475e-02, 5.925702e-02,
  5.819210e-02, 5.717623e-02, 5.620601e-02, 5.527834e-02, 5.439041e-02, 5.353964e-02, 5.272370e-02,
  5.194042e-02, 5.118784e-02, 5.046414e-02, 4.976765e-02, 4.909682e-02, 4.845022e-02, 4.782653e-02,
  4.722454e-02, 4.664308e-02, 4.608112e-02, 4.553766e-02, 4.501179e-02, 4.450264e-02, 4.400940e-02,
  4.353134e-02, 4.306774e-02, 4.261795e-02, 4.218134e-02, 4.175733e-02, 4.134537e-02, 4.094494e-02,
  4.055556e-02, 4.017677e-02, 3.980813e-02, 3.944923e-02, 3.909968e-02, 3.875912e-02, 3.842720e-02,
  3.810358e-02, 3.778796e-02, 3.748003e-02, 3.717951e-02,
};
static const double michel_u[MICHEL_NPTS] = {
  -5.141481e-02, -5.150117e-02, -5.133975e-02, -5.099059e-02, -5.049884e-02, -4.989886e-02,
  -4.921711e-02, -4.847410e-02, -4.768588e-02, -4.686505e-02, -4.602158e-02, -4.516340e-02,
  -4.429679e-02, -4.342677e-02, -4.255735e-02, -4.169170e-02, -4.083238e-02, -3.998138e-02,
  -3.914028e-02, -3.831033e-02, -3.749247e-02, -3.668743e-02, -3.589573e-02, -3.511776e-02,
  -3.435376e-02, -3.360388e-02, -3.286819e-02, -3.214667e-02, -3.143927e-02, -3.074587e-02,
  -3.006633e-02, -2.940047e-02, -2.874810e-02, -2.810899e-02, -2.748291e-02, -2.686962e-02,
  -2.626888e-02, -2.568042e-02, -2.510399e-02, -2.453932e-02, -2.398616e-02, -2.344425e-02,
  -2.291333e-02, -2.239314e-02, -2.188344e-02, -2.138398e-02, -2.089451e-02, -2.041480e-02,
  -1.994460e-02, -1.948371e-02, -1.903188e-02, -1.858890e-02, -1.815457e-02, -1.772866e-02,
  -1.731098e-02, -1.690134e-02, -1.649953e-02, -1.610537e-02, -1.571868e-02, -1.533928e-02,
  -1.496700e-02, -1.460167e-02, -1.424313e-02, -1.389122e-02, -1.354578e-02, -1.320666e-02,
  -1.287373e-02,
};

// Linear interpolation of a tabulated profile in r.
static double
michel_interp(const double *tab, double r)
{
  if (r <= michel_r[0]) {
    return tab[0];
  }
  if (r >= michel_r[MICHEL_NPTS - 1]) {
    return tab[MICHEL_NPTS - 1];
  }
  int i = 0;
  while (michel_r[i + 1] < r) {
    i = i + 1;
  }
  double w = (r - michel_r[i]) / (michel_r[i + 1] - michel_r[i]);
  return tab[i] * (1.0 - w) + tab[i + 1] * w;
}

void
evalDensityInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct michel_ctx *app = ctx;
  double r = xn[0];

  // Rest-frame density of the Michel profile times the metric determinant sqrt(h).
  fout[0] = (r * sqrt((2.0 * app->mass_bh * r) + (r * r))) * michel_interp(michel_n, r);
}

void
evalTempInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct michel_ctx *app = ctx;
  double r = xn[0];

  // Rest-frame temperature of the Michel profile.
  fout[0] = michel_interp(michel_T, r);
}

void
evalVDriftInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct michel_ctx *app = ctx;
  double r = xn[0];

  // Radial spatial four-velocity relative to the Kerr-Schild normal observer; no angular flow.
  fout[0] = michel_interp(michel_u, r);
  fout[1] = 0.0;
  fout[2] = 0.0;
}

void
evalNu(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct michel_ctx *app = ctx;
  // Set collision frequency.
  fout[0] = app->nu;
}

struct vel_map_ctx {
  double vmax; // Physical velocity-space extent.
  double s; // Cubic stretching of the map: v = vc (1 + s (vc/vmax)^2) / (1 + s).
};

void
mapc2p_vel_dir(double t, const double *xc, double *xp, void *ctx)
{
  struct vel_map_ctx *vctx = ctx;
  double vc = xc[0], vm = vctx->vmax, s = vctx->s;
  xp[0] = vc * (1.0 + s * (vc / vm) * (vc / vm)) / (1.0 + s);
}

// Snap a trigger whose next time exceeds t_end only by round-off back to t_end.
static void
snap_trigger_to_t_end(struct gkyl_tm_trigger *trig, double t_end)
{
  if (trig->tcurr > t_end && trig->tcurr <= t_end * (1.0 + 1.0e-10)) {
    trig->tcurr = t_end;
  }
}

void
write_data(struct gkyl_tm_trigger *iot, gkyl_vlasov_app *app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_write) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_vlasov_app_write(app, t_curr, frame);
    gkyl_vlasov_app_write_field_energy(app);
    gkyl_vlasov_app_write_integrated_mom(app);
    gkyl_vlasov_app_write_integrated_L2_f(app);
    gkyl_vlasov_app_write_mom(app, t_curr, frame);
  }
}

void
calc_field_energy(struct gkyl_tm_trigger *fet, gkyl_vlasov_app *app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(fet, t_curr) || force_calc) {
    gkyl_vlasov_app_calc_field_energy(app, t_curr);
  }
}

void
calc_integrated_mom(
  struct gkyl_tm_trigger *imt, gkyl_vlasov_app *app, double t_curr, bool force_calc
)
{
  if (gkyl_tm_trigger_check_and_bump(imt, t_curr) || force_calc) {
    gkyl_vlasov_app_calc_integrated_mom(app, t_curr);
  }
}

void
calc_integrated_L2_f(
  struct gkyl_tm_trigger *l2t, gkyl_vlasov_app *app, double t_curr, bool force_calc
)
{
  if (gkyl_tm_trigger_check_and_bump(l2t, t_curr) || force_calc) {
    gkyl_vlasov_app_calc_integrated_L2_f(app, t_curr);
  }
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

  struct michel_ctx ctx = create_ctx(); // Context for initialization functions.

  int NR = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nr);
  int NVR = APP_ARGS_CHOOSE(app_args.vcells[0], ctx.Nvr);
  int NVTHETA = APP_ARGS_CHOOSE(app_args.vcells[1], ctx.Nvtheta);
  int NVPHI = APP_ARGS_CHOOSE(app_args.vcells[2], ctx.Nvphi);

  int nrank = 1; // Number of processors in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif

  int ccells[] = {NR};
  int cdim = sizeof(ccells) / sizeof(ccells[0]);

  int cuts[cdim];
#ifdef GKYL_HAVE_MPI
  for (int d = 0; d < cdim; d++) {
    if (app_args.use_mpi) {
      cuts[d] = app_args.cuts[d];
    } else {
      cuts[d] = 1;
    }
  }
#else
  for (int d = 0; d < cdim; d++) {
    cuts[d] = 1;
  }
#endif

  // Construct communicator for use in app.
  struct gkyl_comm *comm;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_gpu && app_args.use_mpi) {
#ifdef GKYL_HAVE_NCCL
    comm = gkyl_nccl_comm_new(&(struct gkyl_nccl_comm_inp){.mpi_comm = MPI_COMM_WORLD});
#else
    printf(" Using -g and -M together requires NCCL.\n");
    assert(0 == 1);
#endif
  } else if (app_args.use_mpi) {
    comm = gkyl_mpi_comm_new(&(struct gkyl_mpi_comm_inp){.mpi_comm = MPI_COMM_WORLD});
  } else {
    comm = gkyl_null_comm_inew(&(struct gkyl_null_comm_inp){.use_gpu = app_args.use_gpu});
  }
#else
  comm = gkyl_null_comm_inew(&(struct gkyl_null_comm_inp){.use_gpu = app_args.use_gpu});
#endif

  int my_rank;
  gkyl_comm_get_rank(comm, &my_rank);
  int comm_size;
  gkyl_comm_get_size(comm, &comm_size);

  int ncuts = 1;
  for (int d = 0; d < cdim; d++) {
    ncuts *= cuts[d];
  }

  if (ncuts != comm_size) {
    if (my_rank == 0) {
      fprintf(stderr, "*** Number of ranks, %d, does not match total cuts, %d!\n", comm_size, ncuts);
    }
    goto mpifinalize;
  }

  // Cubic velocity maps: linear core at the origin, stretched tails.
  struct vel_map_ctx vmap_ctx_vr = {.vmax = ctx.vr_max, .s = 1.5};
  struct vel_map_ctx vmap_ctx_vtheta = {.vmax = ctx.vtheta_max, .s = 1.5};
  struct vel_map_ctx vmap_ctx_vphi = {.vmax = ctx.vphi_max, .s = 1.5};

  // Neutral species.
  struct gkyl_vlasov_kinetic_species neut = {
    .model_id = GKYL_MODEL_TRIAD_GR,
    .lower = {-ctx.vr_max, -ctx.vtheta_max, -ctx.vphi_max},
    .upper = {ctx.vr_max, ctx.vtheta_max, ctx.vphi_max},
    .cells = {NVR, NVTHETA, NVPHI},

    .mapc2p_vel =
      {{.mapc2p_vel_func = mapc2p_vel_dir, .mapc2p_vel_ctx = &vmap_ctx_vr},
       {.mapc2p_vel_func = mapc2p_vel_dir, .mapc2p_vel_ctx = &vmap_ctx_vtheta},
       {.mapc2p_vel_func = mapc2p_vel_dir, .mapc2p_vel_ctx = &vmap_ctx_vphi}},

    .num_init = 1,
    .projection[0] =
      {
        .proj_id = GKYL_PROJ_VLASOV_LTE,
        .density = evalDensityInit,
        .ctx_density = &ctx,
        .temp = evalTempInit,
        .ctx_temp = &ctx,
        .V_drift = evalVDriftInit,
        .ctx_V_drift = &ctx,
        .correct_all_moms = true,
        .iter_eps = 0.0,
        .max_iter = 0,
        .use_last_converged = false,
      },

    .collisions =
      {
        .collision_id = GKYL_BGK_COLLISIONS,
        .self_nu = evalNu,
        .self_nu_ctx = &ctx,
        .is_implicit = true,
      },

    .correct =
      {.correct_all_moms = true, .iter_eps = 1.0e-12, .max_iter = 100, .use_last_converged = false},

    .bcx = {.lower = {.type = GKYL_SPECIES_ABSORB}, .upper = {.type = GKYL_SPECIES_FIXED_FUNC}},

    .num_diag_moments = 4,
    .diag_moments = {GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_LTE, GKYL_F_MOMENT_ENERGY},
  };

  // Preset geometry.
  struct gkyl_vlasov_geom geom = {
    .use_preset_geom = true,
    .triad_preset_geom_type = GKYL_TRIAD_GR_KERR_SCHILD_R,
    .mass_bh = ctx.mass_bh,
    .spin_bh = 0.0,
  };

  // Vlasov app (no field).
  struct gkyl_vm app_inp = {
    .geom = geom,

    .cdim = 1,
    .vdim = 3,
    .lower = {ctx.r_min},
    .upper = {ctx.r_max},
    .cells = {NR},

    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,
    .cfl_frac = ctx.cfl_frac,

    .num_periodic_dir = 0,
    .periodic_dirs = {},

    .num_species = 1,
    .species = {{
      .name = "neut",
      .charge = ctx.charge,
      .mass = ctx.mass,
      .type = GKYL_SPECIES_VLASOV,
      .kinetic = neut,
    }},

    .skip_field = true,

    .parallelism = {.use_gpu = app_args.use_gpu, .cuts = {app_args.cuts[0]}, .comm = comm},
  };

  // Create app object.
  // Set app output name from the executable name (argv[0]).
  snprintf(app_inp.name, sizeof(app_inp.name), "%s", app_args.app_name);
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&app_inp);

  // Initial and final simulation times.
  double t_curr = 0.0, t_end = ctx.t_end;

  // Initialize simulation.
  int frame_curr = 0;
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status =
      gkyl_vlasov_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_vlasov_app_cout(
        app, stderr, "*** Failed to read restart file! (%s)\n",
        gkyl_array_rio_status_msg(status.io_status)
      );
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_vlasov_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_vlasov_app_cout(app, stdout, " at time = %g\n", t_curr);
  } else {
    gkyl_vlasov_app_apply_ic(app, t_curr);
  }

  // Create trigger for field energy.
  int field_energy_calcs = ctx.field_energy_calcs;
  struct gkyl_tm_trigger fe_trig = {
    .dt = t_end / field_energy_calcs,
    .tcurr = t_curr,
    .curr = frame_curr,
  };

  calc_field_energy(&fe_trig, app, t_curr, false);

  // Create trigger for integrated moments.
  int integrated_mom_calcs = ctx.integrated_mom_calcs;
  struct gkyl_tm_trigger im_trig = {
    .dt = t_end / integrated_mom_calcs,
    .tcurr = t_curr,
    .curr = frame_curr,
  };

  calc_integrated_mom(&im_trig, app, t_curr, false);

  // Create trigger for integrated L2 norm of the distribution function.
  int integrated_L2_f_calcs = ctx.integrated_L2_f_calcs;
  struct gkyl_tm_trigger l2f_trig = {
    .dt = t_end / integrated_L2_f_calcs,
    .tcurr = t_curr,
    .curr = frame_curr,
  };

  calc_integrated_L2_f(&l2f_trig, app, t_curr, false);

  // Create trigger for IO.
  int num_frames = ctx.num_frames;
  struct gkyl_tm_trigger io_trig = {
    .dt = t_end / num_frames,
    .tcurr = frame_curr * (t_end / num_frames),
    .curr = frame_curr,
  };

  write_data(&io_trig, app, t_curr, false);

  // Compute initial guess of maximum stable time-step.
  double dt = t_end - t_curr;

  // The requested time-step is shortened near the end of the simulation so that
  // the final step lands exactly on t_end.
  bool is_dt_clipped = false; // Was the requested dt shortened below the stable dt?
  bool is_last_step = true; // Does the requested dt reach t_end?

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = ctx.dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx.num_failures_max;

  long step = 1;
  while ((t_curr < t_end) && (step <= app_args.num_steps)) {
    gkyl_vlasov_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = gkyl_vlasov_update(app, dt);
    gkyl_vlasov_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

    if (!status.success) {
      gkyl_vlasov_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }

    // Only a step that took the full requested dt counts as shortened/final.
    bool took_requested_dt = status.dt_actual == dt;
    bool was_dt_clipped = is_dt_clipped && took_requested_dt;
    if (is_last_step && took_requested_dt) {
      // Avoid round-off leaving t_curr just short of t_end.
      t_curr = t_end;
      // Trigger times are accumulated sums and can exceed t_end by round-off.
      // Snap them so the final frame and diagnostics are still produced.
      snap_trigger_to_t_end(&fe_trig, t_end);
      snap_trigger_to_t_end(&im_trig, t_end);
      snap_trigger_to_t_end(&l2f_trig, t_end);
      snap_trigger_to_t_end(&io_trig, t_end);
    } else {
      t_curr += status.dt_actual;
    }

    // Request the next time-step. If the remaining time fits in one stable step,
    // take exactly the remaining time; if it fits in less than two, split it into
    // two equal steps so the final step is never a sliver of the stable dt.
    double t_left = t_end - t_curr;
    dt = status.dt_suggested;
    is_dt_clipped = false;
    is_last_step = false;
    if (t_left <= dt) {
      dt = t_left;
      is_dt_clipped = true;
      is_last_step = true;
    } else if (t_left < 2.0 * dt) {
      dt = 0.5 * t_left;
      is_dt_clipped = true;
    }

    calc_field_energy(&fe_trig, app, t_curr, false);
    calc_integrated_mom(&im_trig, app, t_curr, false);
    calc_integrated_L2_f(&l2f_trig, app, t_curr, false);
    write_data(&io_trig, app, t_curr, false);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    } else if (!was_dt_clipped && status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_vlasov_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_vlasov_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_vlasov_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_vlasov_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_vlasov_app_cout(
          app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max
        );

        calc_field_energy(&fe_trig, app, t_curr, true);
        calc_integrated_mom(&im_trig, app, t_curr, true);
        calc_integrated_L2_f(&l2f_trig, app, t_curr, true);
        write_data(&io_trig, app, t_curr, true);

        break;
      }
    } else {
      num_failures = 0;
    }

    step += 1;
  }

  calc_field_energy(&fe_trig, app, t_curr, false);
  calc_integrated_mom(&im_trig, app, t_curr, false);
  calc_integrated_L2_f(&l2f_trig, app, t_curr, false);
  write_data(&io_trig, app, t_curr, false);
  gkyl_vlasov_app_stat_write(app);

  struct gkyl_vlasov_stat stat = gkyl_vlasov_app_stat(app);

  gkyl_vlasov_app_cout(app, stdout, "\n");
  gkyl_vlasov_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_vlasov_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_vlasov_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_vlasov_app_cout(
      app, stdout, "  Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]
    );
    gkyl_vlasov_app_cout(
      app, stdout, "  Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]
    );
  }
  gkyl_vlasov_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_vlasov_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_vlasov_app_cout(
    app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm
  );
  gkyl_vlasov_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_vlasov_app_cout(
    app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm
  );
  gkyl_vlasov_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

  gkyl_vlasov_app_cout(app, stdout, "Number of write calls %ld\n", stat.n_io);
  double io_tm =
    stat.field_io_tm + stat.species_io_tm + stat.field_diag_io_tm + stat.species_diag_io_tm;
  gkyl_vlasov_app_cout(app, stdout, "IO time took %g secs \n", io_tm);

freeresources:
  // Free resources after simulation completion.
  gkyl_comm_release(comm);
  gkyl_vlasov_app_release(app);

mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif

  return 0;
}
