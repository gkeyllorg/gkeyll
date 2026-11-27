// Buneman instability with PKPM electrons and a DG Euler fluid for the ions (1x), the mixed
// counterpart of vlasov/creg/rt_pkpm_buneman_explicit_1x.c (both species PKPM) and
// vlasov/creg/rt_vlasov_hybrid_buneman_1x1v.c (kinetic electrons in one velocity dimension with the
// same ion fluid and the same electron self-collisions; collisions between a kinetic and a fluid
// species are not supported, so that test has no electron-ion exchange). Parameters follow entry
// JE33 of Ammar's Simulation Journal (https://ammar-hakim.org/sj/je/je33/je33-buneman.html):
// m_i/m_e = 25, electron drift 0.159 so that k v_d = omega_pe for the box mode k = 2 pi, vte =
// 0.02, vti = 0.001, density perturbation 1e-5 of the electrons with the electric field from
// Gauss's law. The electrons are a PKPM species: the drift is their fluid velocity and F is a
// Maxwellian in w = v_par/xi; a uniform magnetic field along x sets the field direction. The ions
// are a fluid (gamma = 5/3, HLLC flux) with the temperature of the kinetic test. The ghost current
// of each species cancels its mean current, so that Ampere's law equals the periodic Poisson solve.
// LBO self-collisions act on the electrons at nu_ee = 1e-4 omega_pe; the electrons and the ion
// fluid exchange momentum and energy at nu_ei = nu_ee and nu_ie = (m_e/m_i) nu_ei. The friction
// slows the drift by 0.4% up to t = 40 (rate nu_ei + nu_ie) and raises the electron temperature by
// 16% (dT_e/dt = (2/3) (m_i/(m_e + m_i)) nu_ei m_e v_d^2).
// The fluid-field coupling is advanced explicitly within the SSP-RK3 step (use_explicit_sources):
// the Lorentz force is in the right-hand side of both species and their currents are on the Ampere
// right-hand side of every stage; rt_pkpm_hybrid_buneman_1x.c is the same problem with the
// staggered operator split.
// Collisionless linear theory with kinetic electrons and the adiabatic ion fluid at k = 2 pi: gamma
// = 0.19764 omega_pe, omega_r = 0.14251 omega_pe (0.19758 and 0.14254 with kinetic ions).
// Figures of merit (16 cells, 64 in w): the field-energy history int E_x^2 dx grows at gamma =
// 0.1875 over t = 20-65 (0.1895 without collisions; both species PKPM: 0.1874; the test with
// kinetic electrons and the same ion fluid, which has no friction and no frictional heating:
// 0.1946). It peaks at 3.5e-3 at t = 72.3 (both PKPM: 3.4e-3 at t = 72.3) and is 1.8e-3 at t = 80.
// Electron and ion number are conserved to 1.2e-12. The total energy (field and the fluid energies
// of both species from the histories) changes by at most 2e-13 of itself up to t = 40 and by at
// most 6e-11 up to t = 80. The mean of E_x is at round-off (below 1e-15) at t = 80, where the
// largest cell average of E_x is 0.085.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>
#include <gkyl_wv_euler.h>
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

struct buneman_ctx {
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using normalized code units).
  double epsilon0; // Permittivity of free space.
  double mu0; // Permeability of free space.
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.
  double mass_ion; // Ion mass.
  double charge_ion; // Ion charge.

  double n0; // Reference number density.
  double Vx_drift_elc; // Electron drift velocity (x-direction).
  double Vx_drift_ion; // Ion drift velocity (x-direction).
  double vte; // Electron thermal velocity.
  double vti; // Ion thermal velocity.
  double gas_gamma; // Adiabatic index of the ion fluid.
  double alpha; // Applied perturbation amplitude.
  double B0; // Magnetic field along x (sets the field direction b).

  double Lx; // Domain size (configuration space: x-direction).
  double kx; // Perturbed wave number (x-direction).
  double omega_pe; // Electron plasma frequency.
  double nu_ee; // Electron-electron collision frequency.
  double nu_ei; // Electron-ion collision frequency.
  double nu_ie; // Ion-electron collision frequency.

  // Simulation parameters.
  int Nx; // Cell count (configuration space: x-direction).
  int Nvx; // Cell count (electron velocity space: normalized parallel velocity w).
  double vx_max_elc; // Domain boundary (electron velocity space: w = v_par/xi).
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

struct buneman_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  double epsilon0 = 1.0; // Permittivity of free space.
  double mu0 = 1.0; // Permeability of free space.
  double mass_elc = 1.0; // Electron mass.
  double charge_elc = -1.0; // Electron charge.
  double mass_ion = 25.0; // Ion mass.
  double charge_ion = 1.0; // Ion charge.

  double n0 = 1.0; // Reference number density.
  double Vx_drift_elc = 0.159; // Electron drift velocity (x-direction).
  double Vx_drift_ion = 0.0; // Ion drift velocity (x-direction).
  double vte = 0.02; // Electron thermal velocity.
  double vti = 0.001; // Ion thermal velocity.
  double gas_gamma = 5.0 / 3.0; // Adiabatic index of the ion fluid.
  double alpha = 1.0e-5; // Applied perturbation amplitude.
  double B0 = 1.0; // Magnetic field along x (sets the field direction b).

  double Lx = 1.0; // Domain size (configuration space: x-direction).
  double kx = 2.0 * pi / Lx; // Perturbed wave number (x-direction).
  double omega_pe =
    sqrt(n0 * charge_ion * charge_ion / (epsilon0 * mass_elc)); // Electron plasma frequency.
  double nu_ee = 1.0e-4 * omega_pe; // Electron-electron collision frequency.
  double nu_ei = nu_ee; // Electron-ion collision frequency.
  double nu_ie = (mass_elc / mass_ion) * nu_ei; // Ion-electron collision frequency.

  // Simulation parameters.
  int Nx = 16; // Cell count (configuration space: x-direction).
  int Nvx = 64; // Cell count (electron velocity space: normalized parallel velocity w).
  double vx_max_elc = 8.0; // Domain boundary (electron velocity space: w = v_par/xi).
  int poly_order = 1; // Polynomial order.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 80.0 / omega_pe; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs =
    INT_MAX; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct buneman_ctx ctx = {
    .pi = pi,
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .mass_ion = mass_ion,
    .charge_ion = charge_ion,
    .n0 = n0,
    .Vx_drift_elc = Vx_drift_elc,
    .Vx_drift_ion = Vx_drift_ion,
    .vte = vte,
    .vti = vti,
    .gas_gamma = gas_gamma,
    .alpha = alpha,
    .B0 = B0,
    .Lx = Lx,
    .kx = kx,
    .omega_pe = omega_pe,
    .nu_ee = nu_ee,
    .nu_ei = nu_ei,
    .nu_ie = nu_ie,
    .Nx = Nx,
    .Nvx = Nvx,
    .vx_max_elc = vx_max_elc,
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
  struct buneman_ctx *app = ctx;
  double w = xn[1];

  double pi = app->pi;

  double F = exp(-0.5 * w * w) / sqrt(2.0 * pi); // Normalized Maxwellian in w = v_par/xi.

  // Set normalized electron distribution functions [F, G] (G = (T_perp/T) F = F).
  fout[0] = F;
  fout[1] = F;

  // No first harmonic: the blocks [M0_perp (x, y, z), N_perp (x, y, z)] vanish.
  fout[2] = 0.0;
  fout[3] = 0.0;
  fout[4] = 0.0;
  fout[5] = 0.0;
  fout[6] = 0.0;
  fout[7] = 0.0;
}

void
evalElcFluidInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct buneman_ctx *app = ctx;
  double x = xn[0];

  double mass_elc = app->mass_elc;
  double n0 = app->n0;
  double alpha = app->alpha;
  double kx = app->kx;
  double Vx_drift_elc = app->Vx_drift_elc;
  double vte = app->vte;

  double rho = mass_elc * n0 * (1.0 + alpha * cos(kx * x)); // Electron mass density.
  double ux = Vx_drift_elc; // Electron drift.
  double E = 0.5 * rho * ux * ux + 1.5 * rho * vte * vte; // Electron total energy density.

  // Set electron fluid state [rho, rho ux, rho uy, rho uz, E].
  fout[0] = rho;
  fout[1] = rho * ux;
  fout[2] = 0.0;
  fout[3] = 0.0;
  fout[4] = E;
}

void
evalIonFluidInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct buneman_ctx *app = ctx;

  double mass_ion = app->mass_ion;
  double n0 = app->n0;
  double Vx_drift_ion = app->Vx_drift_ion;
  double vti = app->vti;

  double rho = mass_ion * n0; // Ion mass density.
  double ux = Vx_drift_ion; // Ion drift.
  // Ion total energy density.
  double E = 0.5 * rho * ux * ux + rho * vti * vti / (app->gas_gamma - 1.0);

  // Set ion fluid state [rho, rho ux, rho uy, rho uz, E].
  fout[0] = rho;
  fout[1] = rho * ux;
  fout[2] = 0.0;
  fout[3] = 0.0;
  fout[4] = E;
}

void
evalNuEE(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct buneman_ctx *app = ctx;

  // Set electron-electron collision frequency.
  fout[0] = app->nu_ee;
}

void
evalNuEI(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct buneman_ctx *app = ctx;

  // Set electron-ion collision frequency.
  fout[0] = app->nu_ei;
}

void
evalNuIE(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct buneman_ctx *app = ctx;

  // Set ion-electron collision frequency.
  fout[0] = app->nu_ie;
}

void
evalFieldInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct buneman_ctx *app = ctx;
  double x = xn[0];

  double alpha = app->alpha;
  double kx = app->kx;

  double Ex = -alpha * sin(kx * x) / kx; // Electric field from Gauss's law.

  // Set electric field, magnetic field and correction potentials.
  fout[0] = Ex;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = 0.0;
  fout[4] = 0.0;
  fout[5] = 0.0;
  fout[6] = 0.0;
  fout[7] = 0.0;
}

void
evalExternalFieldInit(
  double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx
)
{
  struct buneman_ctx *app = ctx;

  double B0 = app->B0;

  // Set external electric and magnetic fields.
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = B0;
  fout[4] = 0.0;
  fout[5] = 0.0;
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

  struct buneman_ctx ctx = create_ctx(); // Context for initialization functions.

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
  struct gkyl_vlasov_kinetic_species elc_kinetic = {
    .lower = {-ctx.vx_max_elc},
    .upper = {ctx.vx_max_elc},
    .cells = {NVX},

    .num_init = 1,
    .projection[0] = {.proj_id = GKYL_PROJ_FUNC, .func = evalElcDistInit, .ctx_func = &ctx},

    .collisions =
      {
        .collision_id = GKYL_LBO_COLLISIONS,
        .self_nu = evalNuEE,
        .self_nu_ctx = &ctx,
        .num_cross_collisions = 1,
        .collide_with = {"ion"},
        .cross_nu = {evalNuEI},
        .cross_nu_ctx = {&ctx},
      },
  };

  struct gkyl_vlasov_fluid_species elc_fluid = {.init = evalElcFluidInit, .ctx = &ctx};

  // Ions (adiabatic fluid, HLLC fluxes): they collide with the electrons through the exchange.
  struct gkyl_wv_euler_inp euler_inp = {
    .gas_gamma = ctx.gas_gamma,
    .rp_type = WV_EULER_RP_HLLC,
    .use_gpu = app_args.use_gpu,
  };
  struct gkyl_wv_eqn *ion_euler = gkyl_wv_euler_inew(&euler_inp);
  struct gkyl_vlasov_fluid_species ion_fluid = {
    .equation = ion_euler,
    .init = evalIonFluidInit,
    .ctx = &ctx,

    .collisions =
      {
        .num_cross_collisions = 1,
        .collide_with = {"elc"},
        .cross_nu = {evalNuIE},
        .cross_nu_ctx = {&ctx},
      },
  };

  // Field.
  struct gkyl_vlasov_field field = {
    .epsilon0 = ctx.epsilon0,
    .mu0 = ctx.mu0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,

    .init = evalFieldInit,
    .ctx = &ctx,

    .ext_em = evalExternalFieldInit,
    .ext_em_ctx = &ctx,

    .use_ghost_current = true,
  };

  // Vlasov app with a PKPM species (electrons) and a fluid species (ions).
  struct gkyl_vm app_inp = {

    .cdim = 1,
    .vdim = 1,
    .lower = {0.0},
    .upper = {ctx.Lx},
    .cells = {NX},

    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,
    .cfl_frac = ctx.cfl_frac,

    .use_explicit_sources = true,

    .num_periodic_dir = 1,
    .periodic_dirs = {0},

    .num_species = 2,
    .species =
      {{
         .name = "elc",
         .charge = ctx.charge_elc,
         .mass = ctx.mass_elc,
         .type = GKYL_SPECIES_PKPM,
         .kinetic = elc_kinetic,
         .fluid = elc_fluid,
       },
       {
         .name = "ion",
         .charge = ctx.charge_ion,
         .mass = ctx.mass_ion,
         .type = GKYL_SPECIES_FLUID,
         .fluid = ion_fluid,
       }},

    .field = field,

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
    .tcurr = gkyl_tm_trigger_anchor(t_end / num_frames, t_curr),
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
  gkyl_vlasov_app_cout(
    app, stdout, "Fluid species RHS calc took %g secs\n", stat.fluid_species_rhs_tm
  );
  gkyl_vlasov_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_vlasov_app_cout(app, stdout, "Species PKPM vars took %g secs\n", stat.species_pkpm_vars_tm);
  gkyl_vlasov_app_cout(
    app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm
  );
  gkyl_vlasov_app_cout(
    app, stdout, "EM variables (bvar) calc took %g secs\n", stat.field_em_vars_tm
  );
  gkyl_vlasov_app_cout(
    app, stdout, "Current evaluation and accumulate took %g secs\n", stat.current_tm
  );
  gkyl_vlasov_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

  gkyl_vlasov_app_cout(app, stdout, "Number of write calls %ld\n", stat.n_io);
  gkyl_vlasov_app_cout(app, stdout, "IO time took %g secs \n", stat.species_io_tm);

freeresources:
  // Free resources after simulation completion.
  gkyl_comm_release(comm);
  gkyl_vlasov_app_release(app);
  gkyl_wv_eqn_release(ion_euler);

mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif

  return 0;
}
