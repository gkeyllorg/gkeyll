// Parallel ion cyclotron wave with kinetic ions and fluid electrons (hybrid Vlasov-Maxwell, 1x3v).
// A left-hand circularly polarized wave propagates along the background magnetic field (x) at k d_i
// = 1 in a plasma with beta_i = beta_e = 0.1, vA / c = 0.1 and mi / me = 100, and is damped by the
// cyclotron resonance of the Maxwellian ions. The fluid-EM coupling is the implicit staggered
// split: the source step takes the ion current as known and solves the electron momentum with E,
// the transport step moves the ions with that E.
//
// Linear theory (Maxwell's equations with the displacement current, kinetic ions, and cold
// electrons, since the pressure does not enter a transverse wave along B): omega^2 - c^2 k^2 +
// omega_pi^2 (omega / (sqrt(2) k v_ti)) Z((omega - Omega_ci) / (sqrt(2) k v_ti)) - omega_pe^2 omega
// / (omega + |Omega_ce|) = 0 gives omega = 0.51572 Omega_ci and a damping rate of 0.05099 Omega_ci.
// Kinetic electrons give 0.51575 and 0.05101.
//
// The wave starts with a magnetic field of 1e-4 B0, Maxwellian ions and an electron fluid with
// transverse velocities u = -1.5953 vA dB / B0 and -0.5992 vA dB / B0, and E = 0.6031 vA (dB / B0)
// x B0. These amplitudes solve the linear initial value problem for which the three other normal
// modes at this wave number and polarization are not excited: the whistler that propagates the
// opposite way (omega = -1.566 Omega_ci) and two high-frequency electron modes. The amplitude
// relations of the wave itself start the whistler, which is not damped, at 4.5% of the field. The
// Maxwellian is not the kinetic eigenfunction: the rest of the perturbation phase-mixes during the
// first period, which leaves the wave 0.985 of the initial field, so the figures of merit are taken
// after it.
//
// Figures of merit (tensor p1 with the hybrid phase basis, 16 x 16 x 6 x 6 cells, 2327 steps,
// omega_pi dt = 0.131, omega_pe dt = |Omega_ce| dt = 1.31, set by the light wave): between t = 1
// and 2.5 wave periods the k d_i = 1 mode of B_y + i B_z decays at 0.05153 Omega_ci (+1.1%) with
// omega = 0.51658 Omega_ci (+0.2%); the history of the transverse magnetic energy over the same
// interval gives 0.05145 Omega_ci (+0.9%). The total energy changes by -5e-13 (7e-5 of the initial
// magnetic energy of the wave).

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

struct hybrid_ion_cyclotron_ctx {
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
  double B0; // Background magnetic field (along x).
  double beta_ion; // Ion plasma beta.
  double Te; // Electron temperature.
  double Ti; // Ion temperature.
  double gas_gamma; // Adiabatic index of the electron fluid.
  double vti; // Ion thermal velocity.
  double vA; // Alfven speed.
  double omega_ci; // Ion cyclotron frequency.

  double alpha; // Amplitude of the magnetic field of the wave, relative to B0.
  double k0; // Wave number.
  double omega_r; // Real frequency of the wave (linear theory).
  double ion_vel_fac; // Ion velocity of the initial condition: u = -ion_vel_fac vA dB / B0.
  double elc_vel_fac; // Electron velocity of the initial condition: u = -elc_vel_fac vA dB / B0.
  double E_fac; // Electric field of the initial condition: E = E_fac vA (dB / B0) x B0.

  // Simulation parameters.
  int Nx; // Cell count (configuration space: x-direction).
  int Nvx; // Cell count (velocity space: vx-direction).
  int Nvy; // Cell count (velocity space: vy-direction).
  int Nvz; // Cell count (velocity space: vz-direction).
  double Lx; // Domain size (configuration space: x-direction).
  double vx_max; // Domain boundary (velocity space: vx-direction).
  double vy_max; // Domain boundary (velocity space: vy-direction).
  double vz_max; // Domain boundary (velocity space: vz-direction).
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

struct hybrid_ion_cyclotron_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  double epsilon0 = 1.0; // Permittivity of free space.
  double mu0 = 1.0; // Permeability of free space.
  double mass_elc = 1.0; // Electron mass.
  double charge_elc = -1.0; // Electron charge.
  double mass_ion = 100.0; // Ion mass.
  double charge_ion = 1.0; // Ion charge.

  double n0 = 1.0; // Reference number density.
  double B0 = 1.0; // Background magnetic field (Alfven speed 0.1 c).
  double beta_ion = 0.1; // Ion plasma beta.
  double Ti = 0.5 * beta_ion * B0 * B0 / (mu0 * n0); // Ion temperature.
  double Te = Ti; // Electron temperature.
  double gas_gamma = 5.0 / 3.0; // Adiabatic index of the electron fluid.
  double vti = sqrt(Ti / mass_ion); // Ion thermal velocity.
  double vA = B0 / sqrt(mu0 * n0 * mass_ion); // Alfven speed.
  double omega_ci = charge_ion * B0 / mass_ion; // Ion cyclotron frequency.
  double omega_pi =
    sqrt(n0 * charge_ion * charge_ion / (epsilon0 * mass_ion)); // Ion plasma frequency.
  double di = 1.0 / (sqrt(epsilon0 * mu0) * omega_pi); // Ion inertial length.

  double alpha = 1.0e-4; // Amplitude of the magnetic field of the wave, relative to B0.
  double k0 = 1.0 / di; // Wave number (k d_i = 1).
  // Root of the dispersion relation of kinetic ions and cold fluid electrons (see the header).
  double omega_r = 0.515720033889 * omega_ci; // Real frequency of the wave.
  // Amplitudes of the initial condition that starts the ion cyclotron wave alone (see the header).
  double ion_vel_fac = 1.595254032963; // Ion velocity: u = -ion_vel_fac vA dB / B0.
  double elc_vel_fac = 0.599201362932; // Electron velocity: u = -elc_vel_fac vA dB / B0.
  double E_fac = 0.603120213183; // Electric field: E = E_fac vA (dB / B0) x B0.

  // Simulation parameters.
  int Nx = 16; // Cell count (configuration space: x-direction).
  int Nvx = 16; // Cell count (velocity space: vx-direction).
  int Nvy = 6; // Cell count (velocity space: vy-direction).
  int Nvz = 6; // Cell count (velocity space: vz-direction).
  double Lx = 2.0 * pi / k0; // Domain size (configuration space: x-direction).
  double vx_max = 5.0 * vti; // Domain boundary (velocity space: vx-direction).
  double vy_max = 4.5 * vti; // Domain boundary (velocity space: vy-direction).
  double vz_max = 4.5 * vti; // Domain boundary (velocity space: vz-direction).
  int poly_order = 1; // Polynomial order.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 2.5 * (2.0 * pi / omega_r); // Final simulation time (2.5 wave periods).
  int num_frames = 5; // Number of output frames (two per wave period).
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs =
    INT_MAX; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct hybrid_ion_cyclotron_ctx ctx = {
    .pi = pi,
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .mass_ion = mass_ion,
    .charge_ion = charge_ion,
    .n0 = n0,
    .B0 = B0,
    .beta_ion = beta_ion,
    .Te = Te,
    .Ti = Ti,
    .gas_gamma = gas_gamma,
    .vti = vti,
    .vA = vA,
    .omega_ci = omega_ci,
    .alpha = alpha,
    .k0 = k0,
    .omega_r = omega_r,
    .ion_vel_fac = ion_vel_fac,
    .elc_vel_fac = elc_vel_fac,
    .E_fac = E_fac,
    .Nx = Nx,
    .Nvx = Nvx,
    .Nvy = Nvy,
    .Nvz = Nvz,
    .Lx = Lx,
    .vx_max = vx_max,
    .vy_max = vy_max,
    .vz_max = vz_max,
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

// Magnetic field of the wave relative to B0, circularly polarized in the sense the ions gyrate.
static void
wave_dB(const struct hybrid_ion_cyclotron_ctx *app, double x, double *dBy, double *dBz)
{
  *dBy = -app->alpha * sin(app->k0 * x);
  *dBz = app->alpha * cos(app->k0 * x);
}

void
evalIonInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct hybrid_ion_cyclotron_ctx *app = ctx;
  double x = xn[0], vx = xn[1], vy = xn[2], vz = xn[3];

  double pi = app->pi;
  double vti = app->vti;
  double n0 = app->n0;

  // Transverse velocity of the initial condition. The ions start as a Maxwellian with this
  // velocity, which is not the kinetic eigenfunction of the damped wave: the rest of the
  // perturbation phase-mixes during the first period.
  double dBy, dBz;
  wave_dB(app, x, &dBy, &dBz);
  double uy = -app->ion_vel_fac * app->vA * dBy;
  double uz = -app->ion_vel_fac * app->vA * dBz;

  // Set distribution function.
  double v_sq = vx * vx + (vy - uy) * (vy - uy) + (vz - uz) * (vz - uz);
  fout[0] = n0 / pow(2.0 * pi * vti * vti, 1.5) * exp(-v_sq / (2.0 * vti * vti));
}

void
evalElcInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct hybrid_ion_cyclotron_ctx *app = ctx;
  double x = xn[0];

  double gas_gamma = app->gas_gamma;
  double n0 = app->n0;

  // Transverse velocity of the initial condition (the density and the pressure are uniform).
  double dBy, dBz;
  wave_dB(app, x, &dBy, &dBz);
  double uy = -app->elc_vel_fac * app->vA * dBy;
  double uz = -app->elc_vel_fac * app->vA * dBz;

  // Set mass density, momentum density and total energy density.
  double rho = app->mass_elc * n0;
  fout[0] = rho;
  fout[1] = 0.0;
  fout[2] = rho * uy;
  fout[3] = rho * uz;
  fout[4] = n0 * app->Te / (gas_gamma - 1.0) + 0.5 * rho * (uy * uy + uz * uz);
}

void
evalFieldInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct hybrid_ion_cyclotron_ctx *app = ctx;
  double x = xn[0];

  double B0 = app->B0;
  double dBy, dBz;
  wave_dB(app, x, &dBy, &dBz);

  // Set electric field: E = E_fac vA (dB / B0) x B0, with B0 along x.
  fout[0] = 0.0;
  fout[1] = app->E_fac * app->vA * B0 * dBz;
  fout[2] = -app->E_fac * app->vA * B0 * dBy;
  // Set magnetic field.
  fout[3] = B0;
  fout[4] = B0 * dBy;
  fout[5] = B0 * dBz;
  // Set correction potentials.
  fout[6] = 0.0;
  fout[7] = 0.0;
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

  struct hybrid_ion_cyclotron_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int NVX = APP_ARGS_CHOOSE(app_args.vcells[0], ctx.Nvx);
  int NVY = APP_ARGS_CHOOSE(app_args.vcells[1], ctx.Nvy);
  int NVZ = APP_ARGS_CHOOSE(app_args.vcells[2], ctx.Nvz);

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

  // Ions (kinetic).
  struct gkyl_vlasov_kinetic_species ion = {
    .lower = {-ctx.vx_max, -ctx.vy_max, -ctx.vz_max},
    .upper = {ctx.vx_max, ctx.vy_max, ctx.vz_max},
    .cells = {NVX, NVY, NVZ},

    .num_init = 1,
    .projection[0] = {.proj_id = GKYL_PROJ_FUNC, .func = evalIonInit, .ctx_func = &ctx},

    .num_diag_moments = 3,
    .diag_moments = {GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2},
  };

  // Electrons (fluid, HLLC fluxes).
  struct gkyl_wv_euler_inp euler_inp = {
    .gas_gamma = ctx.gas_gamma,
    .rp_type = WV_EULER_RP_HLLC,
    .use_gpu = app_args.use_gpu,
  };
  struct gkyl_wv_eqn *elc_euler = gkyl_wv_euler_inew(&euler_inp);
  struct gkyl_vlasov_fluid_species elc = {.equation = elc_euler, .init = evalElcInit, .ctx = &ctx};

  // Field.
  struct gkyl_vlasov_field field = {
    .epsilon0 = ctx.epsilon0,
    .mu0 = ctx.mu0,

    .init = evalFieldInit,
    .ctx = &ctx,
  };

  // Vlasov-Maxwell app. The tensor basis at polynomial order 1 gives the kinetic species the
  // hybrid phase-space basis (order 1 in x, order 2 in v).
  struct gkyl_vm app_inp = {

    .cdim = 1,
    .vdim = 3,
    .lower = {0.0},
    .upper = {ctx.Lx},
    .cells = {NX},

    .poly_order = ctx.poly_order,
    .basis_type = GKYL_BASIS_MODAL_TENSOR,
    .cfl_frac = ctx.cfl_frac,

    .num_periodic_dir = 1,
    .periodic_dirs = {0},

    .num_species = 2,
    .species =
      {{
         .name = "ion",
         .charge = ctx.charge_ion,
         .mass = ctx.mass_ion,
         .type = GKYL_SPECIES_VLASOV,
         .kinetic = ion,
       },
       {
         .name = "elc",
         .charge = ctx.charge_elc,
         .mass = ctx.mass_elc,
         .type = GKYL_SPECIES_FLUID,
         .fluid = elc,
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
  gkyl_wv_eqn_release(elc_euler);
  gkyl_vlasov_app_release(app);

mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif

  return 0;
}