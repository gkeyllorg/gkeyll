// Sheath formation at an electron-emitting wall with the Vlasov-Poisson system of equations, run
// to a steady state. The setup of the sheath test (symmetry boundary at x = 0, conducting wall
// at x = Lx, a boundary-flux source replacing the particles lost to the wall, BGK collisions
// holding the pre-sheath at its source temperature, the position map clustering the cells at
// the wall) with the real ion mass and secondary electron emission from the wall: the energy-
// dependent emission model of Bradshaw and Srinivasan (2024, PSST 33, 035008) with the copper
// fits of Furman and Pivi (true secondaries with a Gaussian spectrum, elastic backscattering)
// scaled to Te = 10 eV.
// Figures of merit at steady state (serendipity p2, 16 x 24 cells, t = 4000/omega_pe): the flux-
// weighted true-secondary yield is 0.188 and the total yield (emitted over incoming electron flux,
// from the net electron wall flux) 0.59, below the Hobbs-Wesson critical yield 1 - 8.3 sqrt(me/mi)
// = 0.81, so the sheath stays monotonic; the midplane-to-wall potential drop e*dphi/Te = 2.30
// against 3.23 without emission, a reduction of 0.93 where the Hobbs-Wesson floating-potential
// shift -ln(1 - delta) gives 0.89; ion wall flux / (n_mid c_s) = 0.72, ion speed at the wall / c_s
// = 2.34; the net electron and ion wall fluxes balance to 1e-4. The emitted distribution and the
// yield diagnostic (<species>_bc_up_yield.gkyl: yield and incoming flux) are written at every
// frame.
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_bc_emission.h>
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

struct sheath_ctx {
  // Physical constants (using normalized code units).
  double epsilon0; // Permittivity of free space.
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.
  double mass_ion; // Ion mass.
  double charge_ion; // Ion charge.

  double n0; // Reference number density.
  double Vx_drift_elc; // Electron drift velocity (x-direction).
  double Vx_drift_ion; // Ion drift velocity (x-direction).

  // Derived physical quantities (using normalized code units).
  double Te; // Electron temperature.
  double Ti; // Ion temperature.

  double vte; // Electron thermal velocity.
  double vti; // Ion thermal velocity.

  double lambda_D; // Electron Debye length.
  double omega_pe; // Electron plasma frequency.

  double nu_ee; // Electron-electron collision frequency (upstream).
  double nu_ii; // Ion-ion collision frequency (upstream).

  // Secondary electron emission from the wall: the copper fits of Furman and Pivi (2002) with
  // every energy expressed in units of the electron temperature (the fits are in eV; here
  // Te = 10 eV). The emitted spectrum is the Gaussian in log-energy, the true-secondary yield
  // is the Furman-Pivi fit and the elastic backscattering is the Furman-Pivi low-energy fit.
  double Te_eV; // Electron temperature the emission fits are scaled with (eV).
  double E_0; // Spectrum: peak energy (units of Te).
  double tau; // Spectrum: width in log-energy.
  double deltahat_ts; // Yield: peak true-secondary yield.
  double Ehat_ts; // Yield: energy of the peak yield (units of Te).
  double t1, t2, t3, t4; // Yield: angular-dependence parameters (unused in 1V).
  double s; // Yield: shape parameter.
  double P1_inf; // Elastic: backscattering probability at high energy.
  double P1_hat; // Elastic: backscattering probability at zero energy.
  double E_hat; // Elastic: energy of the peak backscattering (units of Te).
  double W; // Elastic: decay energy of the backscattering (units of Te).
  double p; // Elastic: shape parameter.
  double t_bound; // Time over which the emission is ramped up from zero.

  // Simulation parameters.
  int Nx; // Cell count (configuration space: x-direction).
  int Nvx; // Cell count (velocity space: vx-direction).
  double Lx; // Domain size (configuration space: x-direction).
  double Ls; // Domain size (source).
  double L_nu; // Extent of the collisional (fixed-temperature) region.
  double vx_max_elc; // Domain boundary (electron velocity space: vx-direction).
  double vx_max_ion; // Domain boundary (ion velocity space: vx-direction).
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
  double mass_elc = 1.0; // Electron mass.
  double charge_elc = -1.0; // Electron charge.
  double mass_ion = 1836.153; // Ion mass.
  double charge_ion = 1.0; // Ion charge.

  double n0 = 1.0; // Reference number density.
  double Vx_drift_elc = 0.0; // Electron drift velocity (x-direction).
  double Vx_drift_ion = 0.0; // Ion drift velocity (x-direction).

  // Derived physical quantities (using normalized code units).
  double Te = 1.0 * charge_ion; // Electron temperature.
  double Ti = 1.0 * charge_ion; // Ion temperature.

  double vte = sqrt(Te / mass_elc); // Electron thermal velocity.
  double vti = sqrt(Ti / mass_ion); // Ion thermal velocity.

  double lambda_D = sqrt(epsilon0 * Te / (n0 * charge_ion * charge_ion)); // Electron Debye length.
  double omega_pe =
    sqrt(n0 * charge_ion * charge_ion / (epsilon0 * mass_elc)); // Electron plasma frequency.

  double nu_ee = vte / (5.0 * lambda_D); // Electron-electron collision frequency (upstream).
  double nu_ii = vti / (5.0 * lambda_D); // Ion-ion collision frequency (upstream).

  // Secondary electron emission (copper, Furman and Pivi 2002), energies in units of Te.
  double Te_eV = 10.0; // Electron temperature the emission fits are scaled with (eV).
  double E_0 = 1.97 / Te_eV; // Spectrum: peak energy (units of Te).
  double tau = 0.88; // Spectrum: width in log-energy.
  double deltahat_ts = 1.885; // Yield: peak true-secondary yield.
  double Ehat_ts = 276.8 / Te_eV; // Yield: energy of the peak yield (units of Te).
  double t1 = 0.66, t2 = 0.8, t3 = 0.7, t4 = 1.0; // Yield: angular-dependence parameters.
  double s = 1.54; // Yield: shape parameter.
  double P1_inf = 0.02; // Elastic: backscattering probability at high energy.
  double P1_hat = 0.496; // Elastic: backscattering probability at zero energy.
  double E_hat = 1.0e-6 / Te_eV; // Elastic: energy of the peak backscattering (units of Te).
  double W = 60.86 / Te_eV; // Elastic: decay energy of the backscattering (units of Te).
  double p = 1.0; // Elastic: shape parameter.

  // Simulation parameters.
  int Nx = 16; // Cell count (configuration space: x-direction).
  int Nvx = 24; // Cell count (velocity space: vx-direction).
  double Lx = 64.0 * lambda_D; // Domain size (configuration space: x-direction).
  double Ls = 40.0 * lambda_D; // Domain size (source).
  double L_nu = 32.0 * lambda_D; // Extent of the collisional (fixed-temperature) region.
  double vx_max_elc = 5.0 * vte; // Domain boundary (electron velocity space: vx-direction).
  double vx_max_ion = 8.0 * vti; // Domain boundary (ion velocity space: vx-direction).
  int poly_order = 2; // Polynomial order.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 4000.0 / omega_pe; // Final simulation time.
  double t_bound = 500.0 / omega_pe; // Time over which the emission is ramped up from zero.
  int num_frames = 1; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs =
    INT_MAX; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct sheath_ctx ctx = {
    .epsilon0 = epsilon0,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .mass_ion = mass_ion,
    .charge_ion = charge_ion,
    .n0 = n0,
    .Vx_drift_elc = Vx_drift_elc,
    .Vx_drift_ion = Vx_drift_ion,
    .Te = Te,
    .Ti = Ti,
    .vte = vte,
    .vti = vti,
    .lambda_D = lambda_D,
    .omega_pe = omega_pe,
    .nu_ee = nu_ee,
    .nu_ii = nu_ii,
    .Te_eV = Te_eV,
    .E_0 = E_0,
    .tau = tau,
    .deltahat_ts = deltahat_ts,
    .Ehat_ts = Ehat_ts,
    .t1 = t1,
    .t2 = t2,
    .t3 = t3,
    .t4 = t4,
    .s = s,
    .P1_inf = P1_inf,
    .P1_hat = P1_hat,
    .E_hat = E_hat,
    .W = W,
    .p = p,
    .t_bound = t_bound,
    .Nx = Nx,
    .Nvx = Nvx,
    .Lx = Lx,
    .Ls = Ls,
    .L_nu = L_nu,
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
evalElcDensityInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct sheath_ctx *app = ctx;

  double n0 = app->n0;

  // Set electron total number density.
  fout[0] = n0;
}

void
evalElcTempInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct sheath_ctx *app = ctx;

  double Te = app->Te;

  // Set electron isotropic temperature..
  fout[0] = Te;
}

void
evalElcVDriftInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct sheath_ctx *app = ctx;

  double Vx_drift_elc = app->Vx_drift_elc;

  // Set electron drift velocity.
  fout[0] = Vx_drift_elc;
}

void
evalElcSourceDensityInit(
  double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx
)
{
  struct sheath_ctx *app = ctx;
  double x = xn[0];

  double Ls = app->Ls;

  double n = 0.0;

  if (x < Ls) {
    n = 2.0 * (Ls - x) / Ls; // Electron source total number density (inside the source region).
  } else {
    n = 0.0; // Electron source total number density (outside the source region).
  }

  // Set electron source total number density.
  fout[0] = n;
}

void
evalElcSourceTempInit(
  double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx
)
{
  struct sheath_ctx *app = ctx;

  double Te = app->Te;

  // Set electron source isotropic temperature.
  fout[0] = Te;
}

void
evalElcSourceVDriftInit(
  double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx
)
{
  struct sheath_ctx *app = ctx;

  double Vx_drift_elc = app->Vx_drift_elc;

  // Set electron source drift velocity.
  fout[0] = Vx_drift_elc;
}

void
evalIonDensityInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct sheath_ctx *app = ctx;

  double n0 = app->n0;

  // Set ion total number density.
  fout[0] = n0;
}

void
evalIonTempInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct sheath_ctx *app = ctx;

  double Ti = app->Ti;

  // Set ion isotropic temperature..
  fout[0] = Ti;
}

void
evalIonVDriftInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct sheath_ctx *app = ctx;

  double Vx_drift_ion = app->Vx_drift_ion;

  // Set ion drift velocity.
  fout[0] = Vx_drift_ion;
}

void
evalIonSourceDensityInit(
  double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx
)
{
  struct sheath_ctx *app = ctx;
  double x = xn[0];

  double Ls = app->Ls;

  double n = 0.0;

  if (x < Ls) {
    n = 2.0 * (Ls - x) / Ls; // Ion source total number density (inside the source region).
  } else {
    n = 0.0; // Ion source total number density (outside the source region).
  }

  // Set ion source total number density.
  fout[0] = n;
}

void
evalIonSourceTempInit(
  double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx
)
{
  struct sheath_ctx *app = ctx;

  double Ti = app->Ti;

  // Set ion source isotropic temperature..
  fout[0] = Ti;
}

void
evalIonSourceVDriftInit(
  double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx
)
{
  struct sheath_ctx *app = ctx;

  double Vx_drift_ion = app->Vx_drift_ion;

  // Set ion source drift velocity.
  fout[0] = Vx_drift_ion;
}

// Position map clustering the cells at the wall (x = Lx): the cell size goes from
// (1 + A) * Lx / Nx at the symmetry boundary to (1 - A) * Lx / Nx at the wall.
void
mapc2p_pos(double t, const double *GKYL_RESTRICT xc, double *GKYL_RESTRICT xp, void *ctx)
{
  struct sheath_ctx *app = ctx;
  double x_c = xc[0];

  double Lx = app->Lx;
  double A = 0.8;

  xp[0] = (1.0 + A) * x_c - A * (x_c * x_c) / Lx;
}

void
evalElcNu(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct sheath_ctx *app = ctx;

  double x = xn[0];

  double nu_ee = app->nu_ee;
  double lambda_D = app->lambda_D;
  double L_nu = app->L_nu;

  double nu = nu_ee / (1.0 + exp((x - L_nu) / (6.0 * lambda_D))); // Electron collision frequency.

  // Set electron collision frequency.
  fout[0] = nu;
}

void
evalIonNu(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct sheath_ctx *app = ctx;

  double x = xn[0];

  double nu_ii = app->nu_ii;
  double lambda_D = app->lambda_D;
  double L_nu = app->L_nu;

  double nu = nu_ii / (1.0 + exp((x - L_nu) / (6.0 * lambda_D))); // Ion collision frequency.

  // Set ion collision frequency.
  fout[0] = nu;
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

  // Secondary electron emission at the wall, driven by the electrons hitting it. The models
  // take the unit charge that converts m v^2 / 2 to the energy unit of the fits; with the fits
  // scaled to Te and Te = 1 in code units, that charge is 1.
  char impact_species[1][128] = {"elc"};
  struct gkyl_emission_spectrum_model *spectrum_model[1] = {
    gkyl_emission_spectrum_gaussian_new(1.0, ctx.E_0, ctx.tau, app_args.use_gpu),
  };
  struct gkyl_emission_yield_model *yield_model[1] = {
    gkyl_emission_yield_furman_pivi_new(
      1.0, ctx.deltahat_ts, ctx.Ehat_ts, ctx.t1, ctx.t2, ctx.t3, ctx.t4, ctx.s, app_args.use_gpu
    ),
  };
  struct gkyl_emission_elastic_model *elastic_model = gkyl_emission_elastic_furman_pivi_new(
    1.0, ctx.P1_inf, ctx.P1_hat, ctx.E_hat, ctx.W, ctx.p, app_args.use_gpu
  );
  struct gkyl_bc_emission_ctx *emission = gkyl_bc_emission_new(
    1, ctx.t_bound, true, spectrum_model, yield_model, elastic_model, impact_species
  );

  // Electrons.
  struct gkyl_vlasov_kinetic_species elc = {
    .lower = {-ctx.vx_max_elc},
    .upper = {ctx.vx_max_elc},
    .cells = {NVX},

    .num_init = 1,
    .projection[0] =
      {
        .proj_id = GKYL_PROJ_VLASOV_LTE,
        .density = evalElcDensityInit,
        .ctx_density = &ctx,
        .temp = evalElcTempInit,
        .ctx_temp = &ctx,
        .V_drift = evalElcVDriftInit,
        .ctx_V_drift = &ctx,
        .correct_all_moms = true,
      },
    .collisions =
      {
        .collision_id = GKYL_BGK_COLLISIONS,
        .self_nu = evalElcNu,
        .self_nu_ctx = &ctx,
        .fixed_temp_relax = true,
      },

    .source =
      {
        .source_id = GKYL_BFLUX_SOURCE,
        .source_length = ctx.Ls,
        .source_species = "ion",

        .num_sources = 1,
        .projection[0] =
          {
            .proj_id = GKYL_PROJ_VLASOV_LTE,
            .density = evalElcSourceDensityInit,
            .ctx_density = &ctx,
            .temp = evalElcSourceTempInit,
            .ctx_temp = &ctx,
            .V_drift = evalElcSourceVDriftInit,
            .ctx_V_drift = &ctx,
          },
      },

    .bcx =
      {
        .lower = {.type = GKYL_SPECIES_REFLECT},
        .upper = {.type = GKYL_SPECIES_EMISSION, .aux_ctx = emission},
      },

    .num_diag_moments = 3,
    .diag_moments = {GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2},
  };

  // Ions.
  struct gkyl_vlasov_kinetic_species ion = {
    .lower = {-ctx.vx_max_ion},
    .upper = {ctx.vx_max_ion},
    .cells = {NVX},

    .num_init = 1,
    .projection[0] =
      {
        .proj_id = GKYL_PROJ_VLASOV_LTE,
        .density = evalIonDensityInit,
        .ctx_density = &ctx,
        .temp = evalIonTempInit,
        .ctx_temp = &ctx,
        .V_drift = evalIonVDriftInit,
        .ctx_V_drift = &ctx,
        .correct_all_moms = true,
      },
    .collisions =
      {
        .collision_id = GKYL_BGK_COLLISIONS,
        .self_nu = evalIonNu,
        .self_nu_ctx = &ctx,
        .fixed_temp_relax = true,
      },

    .source =
      {
        .source_id = GKYL_BFLUX_SOURCE,
        .source_length = ctx.Ls,
        .source_species = "ion",

        .num_sources = 1,
        .projection[0] =
          {
            .proj_id = GKYL_PROJ_VLASOV_LTE,
            .density = evalIonSourceDensityInit,
            .ctx_density = &ctx,
            .temp = evalIonSourceTempInit,
            .ctx_temp = &ctx,
            .V_drift = evalIonSourceVDriftInit,
            .ctx_V_drift = &ctx,
          },
      },

    .bcx = {.lower = {.type = GKYL_SPECIES_REFLECT}, .upper = {.type = GKYL_SPECIES_ABSORB}},

    .num_diag_moments = 3,
    .diag_moments = {GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2},
  };

  // Field.
  struct gkyl_vlasov_field field = {
    .epsilon0 = ctx.epsilon0,

    .poisson_bcs =
      {
        .lo_type = {GKYL_POISSON_NEUMANN},
        .up_type = {GKYL_POISSON_DIRICHLET},

        .lo_value = {0.0},
        .up_value = {0.0},
      },
  };

  // Vlasov-Poisson app.
  struct gkyl_vm app_inp = {

    .cdim = 1,
    .vdim = 1,
    .lower = {0.0},
    .upper = {ctx.Lx},
    .cells = {NX},
    .mapc2p_pos[0] = {.mapc2p_pos_func = mapc2p_pos, .mapc2p_pos_ctx = &ctx},

    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,
    .cfl_frac = ctx.cfl_frac,

    .num_periodic_dir = 0,
    .periodic_dirs = {},

    .num_species = 2,
    .species =
      {{
         .name = "elc",
         .charge = ctx.charge_elc,
         .mass = ctx.mass_elc,
         .type = GKYL_SPECIES_VLASOV,
         .kinetic = elc,
       },
       {
         .name = "ion",
         .charge = ctx.charge_ion,
         .mass = ctx.mass_ion,
         .type = GKYL_SPECIES_VLASOV,
         .kinetic = ion,
       }},

    .field = field,
    .is_electrostatic = true,

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
  gkyl_bc_emission_release(emission);
  gkyl_emission_spectrum_model_release(spectrum_model[0]);
  gkyl_emission_yield_model_release(yield_model[0]);
  gkyl_emission_elastic_model_release(elastic_model);

mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif

  return 0;
}
