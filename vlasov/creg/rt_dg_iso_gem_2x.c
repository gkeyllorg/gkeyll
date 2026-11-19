// Geospace Environment Modeling (GEM) magnetic reconnection challenge with two DG isothermal Euler
// fluids and the Maxwell field of the Vlasov app (2x), the counterpart of
// moments/creg/rt_iso_gem.c. The initial conditions follow J. Birn et al. (2001), "Geospace
// Environmental Modeling (GEM) Magnetic Reconnection Challenge", Journal of Geophysical Research:
// Space Physics 106 (A3): 3715-3719: a Harris sheet of half-width 0.5 d_i with background density
// 0.2, T_i / T_e = 5, m_i / m_e = 25, c = 10 v_A, and a flux perturbation of 0.1 B_0 d_i; periodic
// in x, with reflecting walls for the fluids and perfectly conducting walls for the field in y.
// Each fluid keeps the temperature of its initial condition (p = v_t^2 rho). The fluid-field
// coupling is the staggered operator split with the implicit fluid-Maxwell source step, in which
// the isothermal fluids carry no energy. The fluids use the Lax flux of the isothermal wave
// equation object (rp_type) with the vertex positivity limiter of the density, and the field its
// slope limiter (limit_em) and divergence cleaning of the magnetic field (mgnErrorSpeedFactor = 1).
//
// Figures of merit at t = 25 / Omega_ci (64 x 32 cells, serendipity p1, 4167 steps): the
// reconnected flux int |B_y(x, y = 0)| dx / 2 grows from 0.0200 to 0.1561 (the moments finite-
// volume solution on the same cells: 0.0200 to 0.1265). The out-of-plane field reaches |B_z| =
// 0.027 and the electron speed 0.20 (cell averages); the cell averages of the correction potential
// stay below 2.2e-3, the densities between 0.12 and 2.4, and the field energy falls from 1.511 to
// 1.052. The run amplifies round-off, so the late-time figures agree between runs only to a few
// percent. B_y changes sign three times along the midplane (once in the moments solution). The
// mirror symmetry about the midplane y = 0, kept to round-off by the moments solution, is lost to
// 1.9% of the maximum of the electron z-momentum.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>
#include <gkyl_util.h>
#include <gkyl_wv_iso_euler.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#ifdef GKYL_HAVE_NCCL
#include <gkyl_nccl_comm.h>
#endif
#endif

#include <rt_arg_parse.h>

struct gem_ctx {
  // Mathematical constants (dimensionless).
  double pi; // Pi.

  // Physical constants (using normalized code units).
  double epsilon0; // Permittivity of free space.
  double mu0; // Permeability of free space.
  double mass_ion; // Ion mass.
  double charge_ion; // Ion charge.
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.

  double di; // Ion skin depth.
  double n0; // Reference in-plane number density.
  double nb_over_n0; // Background number density / reference number density.
  double B0; // Reference magnetic field strength.
  double Ti_over_Te; // Ion temperature / electron temperature.
  double beta; // Total plasma beta.

  double lambda; // Current sheet width.
  double psi0; // Reference magnetic scalar potential.

  // Derived physical quantities (using normalized code units).
  double Ti_frac; // Fraction of total temperature from ions.
  double Te_frac; // Fraction of total temperature from electrons.
  double T_tot; // Total temperature.
  double omega_ci; // Ion cyclotron frequency.

  // Simulation parameters.
  int Nx; // Cell count (x-direction).
  int Ny; // Cell count (y-direction).
  double Lx; // Domain size (x-direction).
  double Ly; // Domain size (y-direction).
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

struct gem_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI; // Pi.

  // Physical constants (using normalized code units).
  double epsilon0 = 1.0; // Permittivity of free space.
  double mu0 = 1.0; // Permeability of free space.
  double mass_ion = 1.0; // Ion mass.
  double charge_ion = 1.0; // Ion charge.
  double mass_elc = 1.0 / 25.0; // Electron mass.
  double charge_elc = -1.0; // Electron charge.

  double di = 1.0; // Ion skin depth.
  double n0 = 1.0; // Reference in-plane number density.
  double nb_over_n0 = 0.2; // Background number density / reference number density.
  double B0 = 0.1; // Reference magnetic field strength.
  double Ti_over_Te = 5.0; // Ion temperature / electron temperature.
  double beta = 1.0; // Total plasma beta.

  double lambda = 0.5 * di; // Current sheet width.
  double psi0 = 0.1 * B0 * di; // Reference magnetic scalar potential.

  // Derived physical quantities (using normalized code units).
  double Ti_frac = Ti_over_Te / (1.0 + Ti_over_Te); // Fraction of total temperature from ions.
  double Te_frac = 1.0 / (1.0 + Ti_over_Te); // Fraction of total temperature from electrons.
  double T_tot = beta * (B0 * B0) / 2.0 / n0; // Total temperature.
  double omega_ci = fabs(charge_ion * B0 / mass_ion); // Ion cyclotron frequency.

  // Simulation parameters.
  int Nx = 64; // Cell count (x-direction).
  int Ny = 32; // Cell count (y-direction).
  double Lx = 25.6 * di; // Domain size (x-direction).
  double Ly = 12.8 * di; // Domain size (y-direction).
  int poly_order = 1; // Polynomial order.
  double cfl_frac = 0.9; // CFL coefficient.

  double t_end = 25.0 / omega_ci; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int field_energy_calcs = 100; // Number of times to calculate field energy.
  int integrated_mom_calcs = 100; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs =
    INT_MAX; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gem_ctx ctx = {
    .pi = pi,
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .mass_ion = mass_ion,
    .charge_ion = charge_ion,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .di = di,
    .n0 = n0,
    .nb_over_n0 = nb_over_n0,
    .B0 = B0,
    .Ti_over_Te = Ti_over_Te,
    .beta = beta,
    .lambda = lambda,
    .psi0 = psi0,
    .Ti_frac = Ti_frac,
    .Te_frac = Te_frac,
    .T_tot = T_tot,
    .omega_ci = omega_ci,
    .Nx = Nx,
    .Ny = Ny,
    .Lx = Lx,
    .Ly = Ly,
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
evalElcInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gem_ctx *app = ctx;
  double y = xn[1];

  double mass = app->mass_elc;
  double charge = app->charge_elc;
  double n0 = app->n0;
  double nb_over_n0 = app->nb_over_n0;
  double B0 = app->B0;
  double lambda = app->lambda;
  double T_frac = app->Te_frac;

  double sech_sq =
    (1.0 / cosh(y / lambda)) * (1.0 / cosh(y / lambda)); // Hyperbolic secant squared.

  double n = n0 * (sech_sq + nb_over_n0); // Total number density.
  double Jz = -(B0 / lambda) * sech_sq; // Total current density (z-direction).

  double rho = n * mass; // Mass density.
  double mom_z = (mass / charge) * Jz * T_frac; // Momentum density (z-direction).

  // Set mass density and momentum density.
  fout[0] = rho;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = mom_z;
}

void
evalIonInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gem_ctx *app = ctx;
  double y = xn[1];

  double mass = app->mass_ion;
  double charge = app->charge_ion;
  double n0 = app->n0;
  double nb_over_n0 = app->nb_over_n0;
  double B0 = app->B0;
  double lambda = app->lambda;
  double T_frac = app->Ti_frac;

  double sech_sq =
    (1.0 / cosh(y / lambda)) * (1.0 / cosh(y / lambda)); // Hyperbolic secant squared.

  double n = n0 * (sech_sq + nb_over_n0); // Total number density.
  double Jz = -(B0 / lambda) * sech_sq; // Total current density (z-direction).

  double rho = n * mass; // Mass density.
  double mom_z = (mass / charge) * Jz * T_frac; // Momentum density (z-direction).

  // Set mass density and momentum density.
  fout[0] = rho;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = mom_z;
}

void
evalFieldInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gem_ctx *app = ctx;
  double x = xn[0], y = xn[1];

  double pi = app->pi;
  double B0 = app->B0;
  double lambda = app->lambda;
  double psi0 = app->psi0;
  double Lx = app->Lx;
  double Ly = app->Ly;

  double Bxb = B0 * tanh(y / lambda); // Magnetic field of the current sheet.

  // Set electric field.
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  // Set magnetic field (current sheet plus the perturbation).
  fout[3] = Bxb - psi0 * (pi / Ly) * cos(2.0 * pi * x / Lx) * sin(pi * y / Ly);
  fout[4] = psi0 * (2.0 * pi / Lx) * sin(2.0 * pi * x / Lx) * cos(pi * y / Ly);
  fout[5] = 0.0;
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

  struct gem_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], ctx.Ny);

  // Electron/ion equations (Lax fluxes): isothermal, with the thermal velocities of the initial
  // temperatures.
  struct gkyl_wv_eqn *elc_iso_euler = gkyl_wv_iso_euler_inew(&(struct gkyl_wv_iso_euler_inp){
    .vt = sqrt(ctx.T_tot * ctx.Te_frac / ctx.mass_elc),
    .rp_type = WV_ISO_EULER_RP_LAX,
    .use_gpu = app_args.use_gpu,
  });
  struct gkyl_wv_eqn *ion_iso_euler = gkyl_wv_iso_euler_inew(&(struct gkyl_wv_iso_euler_inp){
    .vt = sqrt(ctx.T_tot * ctx.Ti_frac / ctx.mass_ion),
    .rp_type = WV_ISO_EULER_RP_LAX,
    .use_gpu = app_args.use_gpu,
  });

  struct gkyl_vlasov_fluid_species elc = {
    .equation = elc_iso_euler,
    .init = evalElcInit,
    .ctx = &ctx,

    .bcy = {GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT},
  };

  struct gkyl_vlasov_fluid_species ion = {
    .equation = ion_iso_euler,
    .init = evalIonInit,
    .ctx = &ctx,

    .bcy = {GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT},
  };

  // Field.
  struct gkyl_vlasov_field field = {
    .epsilon0 = ctx.epsilon0,
    .mu0 = ctx.mu0,
    .mgnErrorSpeedFactor = 1.0,
    .limit_em = true,

    .init = evalFieldInit,
    .ctx = &ctx,

    .bcy = {GKYL_FIELD_PEC_WALL, GKYL_FIELD_PEC_WALL},
  };

  int nrank = 1; // Number of processors in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif

  int ccells[] = {NX, NY};
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

  // Vlasov app (fluids and field).
  struct gkyl_vm app_inp = {

    .cdim = 2,
    .vdim = 0,
    .lower = {-0.5 * ctx.Lx, -0.5 * ctx.Ly},
    .upper = {0.5 * ctx.Lx, 0.5 * ctx.Ly},
    .cells = {NX, NY},

    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,
    .cfl_frac = ctx.cfl_frac,

    .num_periodic_dir = 1,
    .periodic_dirs = {0},

    .num_species = 2,
    .species =
      {{
         .name = "elc",
         .type = GKYL_SPECIES_FLUID,
         .charge = ctx.charge_elc,
         .mass = ctx.mass_elc,
         .fluid = elc,
       },
       {
         .name = "ion",
         .type = GKYL_SPECIES_FLUID,
         .charge = ctx.charge_ion,
         .mass = ctx.mass_ion,
         .fluid = ion,
       }},

    .field = field,

    .parallelism = {.use_gpu = app_args.use_gpu, .cuts = {cuts[0], cuts[1]}, .comm = comm},
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
  gkyl_wv_eqn_release(elc_iso_euler);
  gkyl_wv_eqn_release(ion_iso_euler);
  gkyl_vlasov_app_release(app);

mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif

  return 0;
}
