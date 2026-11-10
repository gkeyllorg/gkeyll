// LAPD in cylindrical coordinates, uniform radial grid and nonuniform velocity grid.
// Reproduces the setup of Sec. 4.3.3 of Francisquez et al., JCP 558 (2026) 114852
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_gyrokinetic_run.h>

#include <rt_arg_parse.h>

struct vel_map_ctx {
  double vpar_max; // Physical parallel velocity extent.
  double mu_max; // Physical magnetic moment extent.
};

struct lapd_cyl_ctx {
  int cdim, vdim; // Dimensionality.

  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.
  double mass_ion; // Ion mass.
  double charge_ion; // Ion charge.

  double n_max; // Peak density.
  double Te_max; // Peak electron temperature.
  double Ti_max; // Peak ion temperature.
  double n_ref; // Reference density.
  double Te_ref; // Reference electron temperature.
  double Ti_ref; // Reference ion temperature.
  double B0; // Magnetic field strength (Tesla).
  double nu_frac; // Collision frequency fraction.

  double S0_src; // Source amplitude.
  double Te_src; // Source electron temperature.
  double Ti_src; // Source ion temperature.
  double r_src; // Source radial extent.
  double L_src; // Source radial transition length.
  double floor_src; // Minimum source intensity.

  double L_perp; // Radial scale length of the profiles.
  double r_min; // Inner radius.
  double r_max; // Outer radius.
  double Lz; // Domain size in z.

  struct vel_map_ctx vmap_elc; // Electron velocity map parameters.
  struct vel_map_ctx vmap_ion; // Ion velocity map parameters.
  double vpar_max_c; // Computational parallel velocity extent.
  double mu_max_c; // Computational magnetic moment extent.

  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  double
    write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct lapd_cyl_ctx
create_ctx(void)
{
  int cdim = 3, vdim = 2; // Dimensionality.

  double eV = GKYL_ELEMENTARY_CHARGE;
  double mass_ion = 3.973 * GKYL_PROTON_MASS; // Ion mass.
  double mass_elc = mass_ion / 400.0; // Electron mass.
  double charge_ion = eV; // Ion charge.
  double charge_elc = -eV; // Electron charge.

  double n_max = 2.0e18; // Peak density.
  double Te_max = 4.0 * eV; // Peak electron temperature.
  double Ti_max = 1.0 * eV; // Peak ion temperature.
  double n_ref = 0.5 * n_max; // Reference density.
  double Te_ref = 0.5 * Te_max; // Reference electron temperature.
  double Ti_ref = 0.5 * Ti_max; // Reference ion temperature.
  double B0 = 0.0398; // Magnetic field strength (Tesla).

  // Collision frequency fraction.
  // Equiv. to 0.1 before the Coulomb logarithm change.
  double nu_frac = 0.0415;

  double diam = 1.2; // Domain diameter.
  double Lz = 18.0; // Domain size in z.
  double r_max = 0.5 * diam; // Outer radius.
  double r_min = 0.05 * r_max; // Inner radius.
  double L_perp = 0.8 * diam; // Radial scale length of the profiles.

  // Physical velocity space extents.
  double vt_elc_ref = sqrt(Te_ref / mass_elc);
  double vt_ion_ref = sqrt(Ti_ref / mass_ion);
  struct vel_map_ctx vmap_elc = {
    .vpar_max = 6.0 * vt_elc_ref,
    .mu_max = 0.5 * mass_elc * pow(6.0 * vt_elc_ref, 2) / B0,
  };
  struct vel_map_ctx vmap_ion = {
    .vpar_max = 6.0 * vt_ion_ref,
    .mu_max = 0.5 * mass_ion * pow(6.0 * vt_ion_ref, 2) / B0,
  };

  // Source parameters.
  double c_s_max = sqrt(Te_max / mass_ion); // Peak sound speed.
  double S0_src = 1.08 * n_max * c_s_max / Lz; // Source amplitude.
  double Te_src = 12.0 * eV; // Source electron temperature.
  double Ti_src = 1.0 * eV; // Source ion temperature.
  double r_src = 0.25; // Source radial extent.
  double L_src = 2.0 * 0.625 / 100.0; // Source radial transition length.
  double floor_src = 0.01; // Minimum source intensity.

  // Number of cells (remove denominator for same resolution as JCP).
  int Nr = 36 / 9, Ntheta = 36 / 9, Nz = 8 / 4;
  int Nvpar = 8 / 2, Nmu = 4 / 2;

  double t_end = 6.0e-3 / 100.0; // Final simulation time.
  int num_frames = 600 / 600; // Number of output frames.
  double write_phase_freq =
    0.2 / 0.2; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num = num_frames * 100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct lapd_cyl_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .mass_ion = mass_ion,
    .charge_ion = charge_ion,
    .n_max = n_max,
    .Te_max = Te_max,
    .Ti_max = Ti_max,
    .n_ref = n_ref,
    .Te_ref = Te_ref,
    .Ti_ref = Ti_ref,
    .B0 = B0,
    .nu_frac = nu_frac,
    .S0_src = S0_src,
    .Te_src = Te_src,
    .Ti_src = Ti_src,
    .r_src = r_src,
    .L_src = L_src,
    .floor_src = floor_src,
    .L_perp = L_perp,
    .r_min = r_min,
    .r_max = r_max,
    .Lz = Lz,
    .vmap_elc = vmap_elc,
    .vmap_ion = vmap_ion,
    .vpar_max_c = 1.0 / sqrt(2.0),
    .mu_max_c = 1.0,
    .cells = {Nr, Ntheta, Nz, Nvpar, Nmu},
    .t_end = t_end,
    .num_frames = num_frames,
    .write_phase_freq = write_phase_freq,
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };

  return ctx;
}

// (1-c_edge)*(1-(2r/L_perp)^2)^3 + c_edge inside r < L_perp/2, c_edge outside.
static inline double
radial_profile(double r, double c_edge, double L_perp)
{
  if (r < 0.5 * L_perp) {
    return (1.0 - c_edge) * pow(1.0 - pow(r / (0.5 * L_perp), 2), 3) + c_edge;
  }
  return c_edge;
}

void
evalDensityInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct lapd_cyl_ctx *app = ctx;
  double r = xn[0];

  // Reseeded on every call, so this is a uniform factor (as in the original input).
  pcg32_random_t rng = gkyl_pcg32_init(0);
  double perturb = 2.0e-3 * (1.0 - 0.5 * gkyl_pcg32_rand_double(&rng));

  fout[0] = app->n_max * radial_profile(r, 1.0 / 20.0, app->L_perp) * (1.0 + perturb);
}

void
evalUparInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
evalTempElcInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct lapd_cyl_ctx *app = ctx;
  fout[0] = app->Te_max * radial_profile(xn[0], 1.0 / 5.0, app->L_perp);
}

void
evalTempIonInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct lapd_cyl_ctx *app = ctx;
  fout[0] = app->Ti_max;
}

void
evalSourceDensityInit(
  double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx
)
{
  struct lapd_cyl_ctx *app = ctx;
  double r = xn[0];

  double floor_src = app->floor_src;
  fout[0] = app->S0_src *
            (floor_src + (1.0 - floor_src) * 0.5 * (1.0 - tanh((r - app->r_src) / app->L_src)));
}

void
evalSourceUparInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
evalSourceTempElcInit(
  double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx
)
{
  struct lapd_cyl_ctx *app = ctx;
  fout[0] = app->Te_src * radial_profile(xn[0], 1.0 / 2.5, app->L_perp);
}

void
evalSourceTempIonInit(
  double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx
)
{
  struct lapd_cyl_ctx *app = ctx;
  fout[0] = app->Ti_src;
}

static inline void
mapc2p(double t, const double *GKYL_RESTRICT zc, double *GKYL_RESTRICT xp, void *ctx)
{
  double r = zc[0], theta = zc[1], z = zc[2];

  // Set physical coordinates (X, Y, Z) from computational coordinates (r, theta, z).
  xp[0] = r * cos(theta);
  xp[1] = r * sin(theta);
  xp[2] = z;
}

void
bfield_func(double t, const double *GKYL_RESTRICT zc, double *GKYL_RESTRICT fout, void *ctx)
{
  struct lapd_cyl_ctx *app = ctx;

  // Set Cartesian components of magnetic field.
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = app->B0;
}

void
mapc2p_vel(double t, const double *GKYL_RESTRICT vc, double *GKYL_RESTRICT vp, void *ctx)
{
  struct vel_map_ctx *vmap = ctx;
  double cvpar = vc[0], cmu = vc[1];

  // Linear in vpar up to vpar_max/2, quadratic beyond.
  if (fabs(cvpar) <= 0.5) {
    vp[0] = vmap->vpar_max * cvpar;
  } else {
    vp[0] = copysign(2.0 * vmap->vpar_max * cvpar * cvpar, cvpar);
  }

  // Quadratic in mu.
  vp[1] = vmap->mu_max * cmu * cmu;
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

  struct lapd_cyl_ctx ctx = create_ctx(); // Context for init functions.

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d = 0; d < ctx.cdim; d++) {
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  }
  for (int d = 0; d < ctx.vdim; d++) {
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim + d]);
  }

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  // Electron species.
  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.charge_elc,
    .mass = ctx.mass_elc,
    .vdim = ctx.vdim,
    .lower = {-ctx.vpar_max_c, 0.0},
    .upper = {ctx.vpar_max_c, ctx.mu_max_c},
    .cells = {cells_v[0], cells_v[1]},
    .polarization_density = ctx.n_ref,

    .mapc2p = {.mapping = mapc2p_vel, .ctx = &ctx.vmap_elc},

    .projection =
      {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
        .density = evalDensityInit,
        .ctx_density = &ctx,
        .upar = evalUparInit,
        .ctx_upar = &ctx,
        .temp = evalTempElcInit,
        .ctx_temp = &ctx,
      },

    .collisionless = {.type = GKYL_GK_COLLISIONLESS_ES},

    .collisions =
      {
        .collision_id = GKYL_LBO_COLLISIONS,
        .den_ref = ctx.n_ref,
        .temp_ref = ctx.Te_ref,
        .nu_frac = ctx.nu_frac,
        .num_cross_collisions = 1,
        .collide_with = {"ion"},
      },

    .source =
      {
        .source_id = GKYL_PROJ_SOURCE,
        .num_sources = 1,
        .projection[0] =
          {
            .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
            .density = evalSourceDensityInit,
            .ctx_density = &ctx,
            .upar = evalSourceUparInit,
            .ctx_upar = &ctx,
            .temp = evalSourceTempElcInit,
            .ctx_temp = &ctx,
          },
        .diagnostics =
          {.num_integrated_diag_moments = 1, .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN}
          },
      },

    .bcs =
      {{.dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
       {.dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX},
       {.dir = 2, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH_CONDUCTING},
       {.dir = 2, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH_CONDUCTING}},

    .num_diag_moments = 4,
    .diag_moments =
      {GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_BIMAXWELLIAN},
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
    .time_rate_diagnostics = true,
    .num_balance_diags = 2,
    .balance_diags = {GKYL_BALANCE_GK_PARTICLE, GKYL_BALANCE_GK_ENERGY},

    .boundary_flux_diagnostics =
      {.num_integrated_diag_moments = 1, .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN}},
  };

  // Ion species.
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.charge_ion,
    .mass = ctx.mass_ion,
    .vdim = ctx.vdim,
    .lower = {-ctx.vpar_max_c, 0.0},
    .upper = {ctx.vpar_max_c, ctx.mu_max_c},
    .cells = {cells_v[0], cells_v[1]},
    .polarization_density = ctx.n_ref,

    .mapc2p = {.mapping = mapc2p_vel, .ctx = &ctx.vmap_ion},

    .projection =
      {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
        .density = evalDensityInit,
        .ctx_density = &ctx,
        .upar = evalUparInit,
        .ctx_upar = &ctx,
        .temp = evalTempIonInit,
        .ctx_temp = &ctx,
      },

    .collisionless = {.type = GKYL_GK_COLLISIONLESS_ES},

    .collisions =
      {
        .collision_id = GKYL_LBO_COLLISIONS,
        .den_ref = ctx.n_ref,
        .temp_ref = ctx.Ti_ref,
        .nu_frac = ctx.nu_frac,
        .num_cross_collisions = 1,
        .collide_with = {"elc"},
      },

    .source =
      {
        .source_id = GKYL_PROJ_SOURCE,
        .num_sources = 1,
        .projection[0] =
          {
            .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
            .density = evalSourceDensityInit,
            .ctx_density = &ctx,
            .upar = evalSourceUparInit,
            .ctx_upar = &ctx,
            .temp = evalSourceTempIonInit,
            .ctx_temp = &ctx,
          },
        .diagnostics =
          {.num_integrated_diag_moments = 1, .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN}
          },
      },

    .bcs =
      {{.dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
       {.dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX},
       {.dir = 2, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH_CONDUCTING},
       {.dir = 2, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH_CONDUCTING}},

    .num_diag_moments = 4,
    .diag_moments =
      {GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_BIMAXWELLIAN},
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
    .time_rate_diagnostics = true,
    .num_balance_diags = 2,
    .balance_diags = {GKYL_BALANCE_GK_PARTICLE, GKYL_BALANCE_GK_ENERGY},

    .boundary_flux_diagnostics =
      {.num_integrated_diag_moments = 1, .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN}},
  };

  // Field.
  struct gkyl_gyrokinetic_field field = {
    .poisson_bcs =
      {{.dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}},
       {.dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}},
       {.dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_PERIODIC},
       {.dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_PERIODIC}},
    .time_rate_diagnostics = true,
  };

  // GK app.
  struct gkyl_gk app_inp = {

    .cdim = ctx.cdim,
    .lower = {ctx.r_min, -M_PI, -0.5 * ctx.Lz},
    .upper = {ctx.r_max, M_PI, 0.5 * ctx.Lz},
    .cells = {cells_x[0], cells_x[1], cells_x[2]},
    .poly_order = 1,
    .basis_type = app_args.basis_type,

    .geometry =
      {
        .geometry_id = GKYL_GEOMETRY_MAPC2P,
        .mapc2p = mapc2p,
        .c2p_ctx = &ctx,
        .bfield_func = bfield_func,
        .bfield_ctx = &ctx,
      },

    .num_periodic_dir = 1,
    .periodic_dirs = {1},

    .diagnostics =
      {.num_balance_diags = 2, .balance_diags = {GKYL_BALANCE_GK_PARTICLE, GKYL_BALANCE_GK_ENERGY}},

    .num_species = 2,
    .species = {elc, ion},
    .field = field,

    .parallelism =
      {
        .use_gpu = app_args.use_gpu,
        .cuts = {app_args.cuts[0], app_args.cuts[1], app_args.cuts[2]},
        .comm = comm,
      },
  };

  // Set app output name from the executable name (argv[0]).
  snprintf(app_inp.name, sizeof(app_inp.name), "%s", app_args.app_name);
  struct gkyl_gyrokinetic_run_inp run_inp = {
    .app_inp = app_inp,
    .time_stepping =
      {
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
    .print_verbosity = {.enabled = false, .frequency = 0.01},
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
