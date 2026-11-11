// Sheath formation with the PKPM system of equations (1x, normalized parallel velocity
// w = v_par/xi), the counterpart of vlasov/creg/rt_vlasov_sheath_1x1v.c.
// Half domain: symmetry boundary at x = 0 (reflecting species, symmetry-wall field),
// absorbing conducting wall at x = Lx. A boundary-flux source replaces the particles lost
// to the wall, and a BGK-form temperature source at the source temperature acts upstream
// (x < L_nu) to hold the pre-sheath plasma at that temperature; the region next to
// the wall is collisionless. The position map clusters the cells at the wall.
// Run to t = 1000/omega_pe, when the wall loss is 0.6% below its steady value.
// Figures of merit at t = 1000 (Vlasov 1x1v in parentheses): midplane-to-wall potential drop
// e*dphi/Te = 1.929 (1.925), ion wall flux / (n_mid c_s) = 0.506 (0.509), ion speed at the
// wall / c_s = 1.635 (1.629), with c_s = sqrt((Te + Ti)/mi). Over the presheath (the 13 cells
// upstream of the three sheath cells) the cell averages of the density and the flow velocity
// of both species are within 1.1% of Vlasov (rms 0.4%), the largest difference being in the
// cell at the sheath edge; the wall cell has the electron density 2.2% low and the electron
// velocity 2.0% high.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_pkpm.h>
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

struct sheath_ctx {
  // Physical constants (using normalized code units).
  double epsilon0; // Permittivity of free space.
  double mu0; // Permeability of free space.
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.
  double mass_ion; // Ion mass.
  double charge_ion; // Ion charge.

  double n0; // Reference number density.
  double Te_over_Ti; // Electron temperature / ion temperature.
  double B0; // Reference magnetic field strength.

  // Derived physical quantities (using normalized code units).
  double vte; // Electron thermal velocity.
  double vti; // Ion thermal velocity.
  double lambda_D; // Electron Debye length.
  double omega_pe; // Electron plasma frequency.
  double nu_Q_elc; // Rate of the electron temperature source (upstream).
  double nu_Q_ion; // Rate of the ion temperature source (upstream).

  // Simulation parameters.
  int Nx; // Cell count (configuration space: x-direction).
  int Nvx; // Cell count (velocity space: normalized parallel velocity w).
  double Lx; // Domain size (configuration space: x-direction).
  double Ls; // Extent of the source region.
  double L_nu; // Extent of the temperature-source region.
  double map_A; // Amplitude of the position map clustering the cells at the wall.
  double vx_max_elc; // Domain boundary (electron velocity space: w = v_par/xi).
  double vx_max_ion; // Domain boundary (ion velocity space: w = v_par/xi).
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

struct sheath_ctx
create_ctx(void)
{
  // Physical constants (using normalized code units).
  double epsilon0 = 1.0; // Permittivity of free space.
  double mu0 = 1.0 / 25.0; // Permeability of free space (speed of light = 5 vte).
  double mass_elc = 1.0; // Electron mass.
  double charge_elc = -1.0; // Electron charge.
  double mass_ion = 100.0; // Ion mass.
  double charge_ion = 1.0; // Ion charge.

  double n0 = 1.0; // Reference number density.
  double Te_over_Ti = 1.0; // Electron temperature / ion temperature.
  double B0 = 1.0; // Reference magnetic field strength.

  // Derived physical quantities (using normalized code units).
  double vte = 1.0; // Electron thermal velocity.
  double vti = vte / sqrt(Te_over_Ti * mass_ion); // Ion thermal velocity.

  double lambda_D = sqrt(
    epsilon0 * vte * vte * mass_elc / (n0 * charge_ion * charge_ion)
  ); // Electron Debye length.
  double omega_pe =
    sqrt(n0 * charge_ion * charge_ion / (epsilon0 * mass_elc)); // Electron plasma frequency.
  double nu_Q_elc = vte / (5.0 * lambda_D); // Rate of the electron temperature source (upstream).
  double nu_Q_ion = vti / (5.0 * lambda_D); // Rate of the ion temperature source (upstream).

  // Simulation parameters.
  int Nx = 16; // Cell count (configuration space: x-direction).
  int Nvx = 16; // Cell count (velocity space: normalized parallel velocity w).
  double Lx = 64.0 * lambda_D; // Domain size (configuration space: x-direction).
  double Ls = 40.0 * lambda_D; // Extent of the source region.
  double L_nu = 32.0 * lambda_D; // Extent of the temperature-source region.
  double map_A = 0.8; // Amplitude of the position map clustering the cells at the wall.
  double vx_max_elc = 5.0; // Domain boundary (electron velocity space: w = v_par/xi).
  double vx_max_ion = 5.0; // Domain boundary (ion velocity space: w = v_par/xi).
  int poly_order = 1; // Polynomial order.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 1000.0 / omega_pe; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs =
    INT_MAX; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct sheath_ctx ctx = {
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .mass_ion = mass_ion,
    .charge_ion = charge_ion,
    .n0 = n0,
    .Te_over_Ti = Te_over_Ti,
    .B0 = B0,
    .vte = vte,
    .vti = vti,
    .lambda_D = lambda_D,
    .omega_pe = omega_pe,
    .nu_Q_elc = nu_Q_elc,
    .nu_Q_ion = nu_Q_ion,
    .Nx = Nx,
    .Nvx = Nvx,
    .Lx = Lx,
    .Ls = Ls,
    .L_nu = L_nu,
    .map_A = map_A,
    .vx_max_elc = vx_max_elc,
    .vx_max_ion = vx_max_ion,
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

void
evalElcDistInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  double w = xn[1]; // Normalized parallel velocity w = v_par/xi.

  // Normalized isotropic Maxwellian: F = exp(-w^2/2)/sqrt(2 pi) (M0 = 1, M2 = 1) and
  // G = (T_perp/T) F = F; the density and temperature enter through the fluid state.
  double F = exp(-0.5 * w * w) / sqrt(2.0 * M_PI);

  // Set normalized electron distribution functions.
  fout[0] = F;
  fout[1] = F;
}

void
evalElcFluidInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct sheath_ctx *app = ctx;

  double rho = app->mass_elc * app->n0; // Electron mass density.
  double E = 1.5 * rho * app->vte * app->vte; // Electron total energy density (no flow).

  // Set electron fluid state [rho, rho ux, rho uy, rho uz, E].
  fout[0] = rho;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = 0.0;
  fout[4] = E;
}

void
evalIonDistInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  double w = xn[1]; // Normalized parallel velocity w = v_par/xi.

  // Normalized isotropic Maxwellian: F = exp(-w^2/2)/sqrt(2 pi) (M0 = 1, M2 = 1) and
  // G = (T_perp/T) F = F; the density and temperature enter through the fluid state.
  double F = exp(-0.5 * w * w) / sqrt(2.0 * M_PI);

  // Set normalized ion distribution functions.
  fout[0] = F;
  fout[1] = F;
}

void
evalIonFluidInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct sheath_ctx *app = ctx;

  double rho = app->mass_ion * app->n0; // Ion mass density.
  double E = 1.5 * rho * app->vti * app->vti; // Ion total energy density (no flow).

  // Set ion fluid state [rho, rho ux, rho uy, rho uz, E].
  fout[0] = rho;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = 0.0;
  fout[4] = E;
}

void
evalFieldInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct sheath_ctx *app = ctx;

  // Set electric field.
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  // Set magnetic field (the guide field is external, see evalExternalFieldInit).
  fout[3] = 0.0;
  fout[4] = 0.0;
  fout[5] = 0.0;
  // Set correction potentials.
  fout[6] = 0.0;
  fout[7] = 0.0;
}

// Static external guide field along x (b = x_hat): not subject to the field boundary
// conditions, so the wall's image field does not reverse it.
void
evalExternalFieldInit(
  double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx
)
{
  struct sheath_ctx *app = ctx;

  // Set external electric field.
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  // Set external magnetic field.
  fout[3] = app->B0;
  fout[4] = 0.0;
  fout[5] = 0.0;
}

// Rates of the temperature source: nu_Q_elc and nu_Q_ion upstream (x < L_nu), vanishing over
// 6 lambda_D towards the collisionless region next to the wall.
void
evalElcSourceRate(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct sheath_ctx *app = ctx;
  double x = xn[0];

  // Set electron temperature-source rate.
  fout[0] = app->nu_Q_elc / (1.0 + exp((x - app->L_nu) / (6.0 * app->lambda_D)));
}

void
evalIonSourceRate(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct sheath_ctx *app = ctx;
  double x = xn[0];

  // Set ion temperature-source rate.
  fout[0] = app->nu_Q_ion / (1.0 + exp((x - app->L_nu) / (6.0 * app->lambda_D)));
}

// Source profile: a density 2 (Ls - x)/Ls inside the source region (it integrates to Ls,
// so the boundary-flux source replaces the wall losses exactly), at the plasma
// temperature and without drift. The temperature source relaxes to the same temperature.
void
evalSourceDensity(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct sheath_ctx *app = ctx;
  double x = xn[0], Ls = app->Ls;

  // Set source number density.
  fout[0] = x < Ls ? 2.0 * (Ls - x) / Ls : 0.0;
}

void
evalElcSourceTemp(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct sheath_ctx *app = ctx;

  // Set electron source temperature.
  fout[0] = app->mass_elc * app->vte * app->vte;
}

void
evalIonSourceTemp(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct sheath_ctx *app = ctx;

  // Set ion source temperature.
  fout[0] = app->mass_ion * app->vti * app->vti;
}

void
evalSourceVDrift(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  // Set source drift velocity.
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
}

// Position map clustering the cells at the wall (x = Lx): the cell size goes from
// (1 + A) Lx/Nx at the symmetry boundary to (1 - A) Lx/Nx at the wall.
void
mapc2p_pos_x(double t, const double *GKYL_RESTRICT xc, double *GKYL_RESTRICT xp, void *ctx)
{
  struct sheath_ctx *app = ctx;
  double Lx = app->Lx, A = app->map_A;
  double x_c = xc[0];

  xp[0] = (1.0 + A) * x_c - A * (x_c * x_c) / Lx;
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
write_data(struct gkyl_tm_trigger *iot, gkyl_pkpm_app *app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_write) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_pkpm_app_write(app, t_curr, frame);
    gkyl_pkpm_app_write_field_energy(app);
    gkyl_pkpm_app_write_integrated_mom(app);
    gkyl_pkpm_app_write_integrated_L2_f(app);
  }
}

void
calc_field_energy(struct gkyl_tm_trigger *fet, gkyl_pkpm_app *app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(fet, t_curr) || force_calc) {
    gkyl_pkpm_app_calc_field_energy(app, t_curr);
  }
}

void
calc_integrated_mom(struct gkyl_tm_trigger *imt, gkyl_pkpm_app *app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(imt, t_curr) || force_calc) {
    gkyl_pkpm_app_calc_integrated_mom(app, t_curr);
  }
}

void
calc_integrated_L2_f(struct gkyl_tm_trigger *l2t, gkyl_pkpm_app *app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(l2t, t_curr) || force_calc) {
    gkyl_pkpm_app_calc_integrated_L2_f(app, t_curr);
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

  struct sheath_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int NVX = APP_ARGS_CHOOSE(app_args.vcells[0], ctx.Nvx);

  int nrank = 1; // Number of processors in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif

  int ccells[] = {NX};
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

  // Electrons.
  struct gkyl_pkpm_species elc = {
    .name = "elc",
    .charge = ctx.charge_elc,
    .mass = ctx.mass_elc,
    .lower = {-ctx.vx_max_elc},
    .upper = {ctx.vx_max_elc},
    .cells = {NVX},

    .init_dist = evalElcDistInit,
    .ctx_dist = &ctx,
    .init_fluid = evalElcFluidInit,
    .ctx_fluid = &ctx,

    .source_bgk =
      {
        .source_bgk_id = GKYL_SOURCE_BGK_HEATING,
        .rate_profile = evalElcSourceRate,
        .rate_profile_ctx = &ctx,
        .temp_shape = evalElcSourceTemp,
        .temp_shape_ctx = &ctx,
      },

    .source =
      {
        .source_id = GKYL_BFLUX_SOURCE,
        .source_length = ctx.Ls,
        .source_species = "ion",
        .write_source = true,
        .density = evalSourceDensity,
        .ctx_density = &ctx,
        .temp = evalElcSourceTemp,
        .ctx_temp = &ctx,
        .V_drift = evalSourceVDrift,
        .ctx_V_drift = &ctx,
      },

    .bcx = {GKYL_SPECIES_REFLECT, GKYL_SPECIES_ABSORB},
  };

  // Ions.
  struct gkyl_pkpm_species ion = {
    .name = "ion",
    .charge = ctx.charge_ion,
    .mass = ctx.mass_ion,
    .lower = {-ctx.vx_max_ion},
    .upper = {ctx.vx_max_ion},
    .cells = {NVX},

    .init_dist = evalIonDistInit,
    .ctx_dist = &ctx,
    .init_fluid = evalIonFluidInit,
    .ctx_fluid = &ctx,

    .source_bgk =
      {
        .source_bgk_id = GKYL_SOURCE_BGK_HEATING,
        .rate_profile = evalIonSourceRate,
        .rate_profile_ctx = &ctx,
        .temp_shape = evalIonSourceTemp,
        .temp_shape_ctx = &ctx,
      },

    .source =
      {
        .source_id = GKYL_BFLUX_SOURCE,
        .source_length = ctx.Ls,
        .source_species = "ion",
        .write_source = true,
        .density = evalSourceDensity,
        .ctx_density = &ctx,
        .temp = evalIonSourceTemp,
        .ctx_temp = &ctx,
        .V_drift = evalSourceVDrift,
        .ctx_V_drift = &ctx,
      },

    .bcx = {GKYL_SPECIES_REFLECT, GKYL_SPECIES_ABSORB},
  };

  // Field.
  struct gkyl_pkpm_field field = {
    .epsilon0 = ctx.epsilon0,
    .mu0 = ctx.mu0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,

    .init = evalFieldInit,
    .ctx = &ctx,
    .ext_em = evalExternalFieldInit,
    .ext_em_ctx = &ctx,

    .bcx = {GKYL_FIELD_SYM_WALL, GKYL_FIELD_PEC_WALL},
  };

  // PKPM app.
  struct gkyl_pkpm app_inp = {

    .cdim = 1,
    .vdim = 1,
    .lower = {0.0},
    .upper = {ctx.Lx},
    .cells = {NX},
    .mapc2p_pos[0] = {.mapc2p_pos_func = mapc2p_pos_x, .mapc2p_pos_ctx = &ctx},

    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,
    .cfl_frac = ctx.cfl_frac,

    .num_periodic_dir = 0,
    .periodic_dirs = {},

    .num_species = 2,
    .species = {elc, ion},

    .field = field,

    .parallelism = {.use_gpu = app_args.use_gpu, .cuts = {app_args.cuts[0]}, .comm = comm},
  };

  // Create app object.
  // Set app output name from the executable name (argv[0]).
  snprintf(app_inp.name, sizeof(app_inp.name), "%s", app_args.app_name);
  gkyl_pkpm_app *app = gkyl_pkpm_app_new(&app_inp);

  // Initial and final simulation times.
  double t_curr = 0.0, t_end = ctx.t_end;

  // Initialize simulation.
  int frame_curr = 0;
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status =
      gkyl_pkpm_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_pkpm_app_cout(
        app, stderr, "*** Failed to read restart file! (%s)\n",
        gkyl_array_rio_status_msg(status.io_status)
      );
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_pkpm_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_pkpm_app_cout(app, stdout, " at time = %g\n", t_curr);
  } else {
    gkyl_pkpm_app_apply_ic(app, t_curr);
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
    gkyl_pkpm_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = gkyl_pkpm_update(app, dt);
    gkyl_pkpm_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

    if (!status.success) {
      gkyl_pkpm_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
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

      gkyl_pkpm_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_pkpm_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_pkpm_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_pkpm_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_pkpm_app_cout(
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
  gkyl_pkpm_app_stat_write(app);

  struct gkyl_pkpm_stat stat = gkyl_pkpm_app_stat(app);

  gkyl_pkpm_app_cout(app, stdout, "\n");
  gkyl_pkpm_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_pkpm_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_pkpm_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_pkpm_app_cout(
      app, stdout, "  Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]
    );
    gkyl_pkpm_app_cout(
      app, stdout, "  Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]
    );
  }
  gkyl_pkpm_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_pkpm_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_pkpm_app_cout(
    app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm
  );
  gkyl_pkpm_app_cout(
    app, stdout, "Fluid species RHS calc took %g secs\n", stat.fluid_species_rhs_tm
  );
  gkyl_pkpm_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_pkpm_app_cout(app, stdout, "Species PKPM vars took %g secs\n", stat.species_pkpm_vars_tm);
  gkyl_pkpm_app_cout(
    app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm
  );
  gkyl_pkpm_app_cout(app, stdout, "EM variables (bvar) calc took %g secs\n", stat.field_em_vars_tm);
  gkyl_pkpm_app_cout(
    app, stdout, "Current evaluation and accumulate took %g secs\n", stat.current_tm
  );
  gkyl_pkpm_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

  gkyl_pkpm_app_cout(app, stdout, "Number of write calls %ld\n", stat.n_io);
  gkyl_pkpm_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

freeresources:
  // Free resources after simulation completion.
  gkyl_comm_release(comm);
  gkyl_pkpm_app_release(app);

mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif

  return 0;
}
