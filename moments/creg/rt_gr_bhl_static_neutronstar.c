// 2D Bondi-Hoyle-Lyttleton accretion problem onto a static neutron star, for the general relativistic Euler equations.
// Input parameters describe wind accretion of a cold relativistic gas onto a non-rotating neutron star.
// Based on the analytical solution for stiff relativistic fluids presented in the article:
// L. I. Petrich, S. L. Shapiro and S. A. Teukolsky (1988), "Accretion onto a moving black hole: An exact solution",
// Physical Review Letters, Volume 60 (18): 1781-1784.
// https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.60.1781

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_gr_euler.h>
#include <gkyl_gr_minkowski.h>
#include <gkyl_gr_neutronstar.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

#include <rt_arg_parse.h>

struct bhl_static_ctx
{
  // Physical constants (using normalized code units).
  double gas_gamma; // Adiabatic index.

  double rhol; // Left fluid mass density.
  double ul; // Left fluid velocity.
  double pl; // Left fluid pressure.

  double rhor; // Right fluid mass density.
  double ur; // Right fluid velocity.
  double pr; // Right fluid pressure.

  // Spacetime parameters (using geometric units).
  double mass; // Mass of the neutron star.
  double spin; // Spin of the neutron star.

  double pos_x; // Position of the neutron star (x-direction).
  double pos_y; // Position of the neutron star (y-direction).
  double pos_z; // Position of the neutron star (z-direction).

  // Derived spacetime parameters (using geometric units).
  double ang_mom; // Angular momentum of the neutron star
  double spin_dimensionless; // Dimensionless spin of the neutron star.

  double alpha; // Neutron star alpha parameter (mass quadrupole moment coefficient).
  double beta; // Neutron star beta parameter (spin octupole moment coefficient).
  double gamma; // Neutron star gamma parameter (mass hexadecapole moment coefficient).

  double mass_quadrupole; // Neutron star mass quadrupole.
  double spin_octupole; // Neutron star spin octupole.
  double mass_hexadecapole; // Neutron star mass hexadecapole;

  // Pointer to spacetime metric.
  struct gkyl_gr_spacetime *spacetime;

  // Simulation parameters.
  int Nx; // Cell count (x-direction).
  int Ny; // Cell count (y-direction).
  double Lx; // Domain size (x-direction).
  double Ly; // Domain size (y-direction).
  double cfl_frac; // CFL coefficient.

  enum gkyl_spacetime_gauge spacetime_gauge; // Spacetime gauge choice.
  int reinit_freq; // Spacetime reinitialization frequency.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  int field_energy_calcs; // Number of times to calculate field energy.
  int integrated_mom_calcs; // Number of times to calculate integrated moments.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.

  double x_loc; // Shock location (x-direction).
};

struct bhl_static_ctx
create_ctx(void)
{
  // Physical constants (using normalized code units).
  double gas_gamma = 5.0 / 3.0; // Adiabatic index.

  double rhol = 3.0; // Left fluid mass density.
  double ul = 0.3; // Left fluid velocity.
  double pl = 0.05; // Left fluid pressure.

  double rhor = 0.01; // Right fluid mass density.
  double ur = 0.0; // Right fluid velocity.
  double pr = 0.01; // Right fluid pressure.

  // Spacetime parameters (using geometric units).
  double mass = 0.3; // Mass of the neutron star.
  double spin = 0.0; // Spin of the neutron star.

  double pos_x = 2.5; // Position of the neutron star (x-direction).
  double pos_y = 2.5; // Position of the neutron star (y-direction).
  double pos_z = 0.0; // Position of the neutron star (z-direction).

  // Derived spacetime parameters.
  double ang_mom = spin * mass; // Angular momentum of the neutron star.
  double spin_dimensionless = ang_mom / (mass * mass); // Dimensionless spin of the neutron star.

  double alpha = 5.0; // Neutron star alpha parameter (mass quadrupole moment coefficient).
  double beta = pow(-0.36 + (1.48 * pow(sqrt(alpha), 0.65)), 3.0); // Neutron star beta parameter (spin octupole moment coefficient).
  double gamma = pow(-4.749 + (0.27613 * pow(sqrt(alpha), 1.5146)) + (5.5168 * pow(sqrt(alpha), 0.22229)), 4.0); // Neutron star gamma parameter (mass hexadecapole moment coefficient).

  double mass_quadrupole = -alpha * (spin_dimensionless * spin_dimensionless) * (mass * mass * mass); // Neutron star mass quadrupole.
  double spin_octupole = -beta * (spin_dimensionless * spin_dimensionless * spin_dimensionless) * (mass * mass * mass * mass); // Neutron star spin octupole.
  double mass_hexadecapole = gamma * (spin_dimensionless * spin_dimensionless * spin_dimensionless * spin_dimensionless) * (mass * mass * mass * mass * mass); // Neutron star mass hexadecapole.

  // Pointer to spacetime metric.
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_neutronstar_new(false, mass, spin, mass_quadrupole, spin_octupole, mass_hexadecapole, pos_x, pos_y, pos_z);

  // Simulation parameters.
  int Nx = 256; // Cell count (x-direction).
  int Ny = 256; // Cell count (y-direction).
  double Lx = 5.0; // Domain size (x-direction).
  double Ly = 5.0; // Domain size (y-direction).
  double cfl_frac = 0.95; // CFL coefficient.

  enum gkyl_spacetime_gauge spacetime_gauge = GKYL_STATIC_GAUGE; // Spacetime gauge choice.
  int reinit_freq = 100; // Spacetime reinitialization frequency.

  double t_end = 15.0; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  double x_loc = 1.0; // Shock location (x-direction).

  struct bhl_static_ctx ctx = {
    .gas_gamma = gas_gamma,
    .rhol = rhol,
    .ul = ul,
    .pl = pl,
    .rhor = rhor,
    .ur = ur,
    .pr = pr,
    .mass = mass,
    .spin = spin,
    .pos_x = pos_x,
    .pos_y = pos_y,
    .pos_z = pos_z,
    .ang_mom = ang_mom,
    .spin_dimensionless = spin_dimensionless,
    .alpha = alpha,
    .beta = beta,
    .gamma = gamma,
    .mass_quadrupole = mass_quadrupole,
    .spin_octupole = spin_octupole,
    .mass_hexadecapole = mass_hexadecapole,
    .spacetime = spacetime,
    .Nx = Nx,
    .Ny = Ny,
    .Lx = Lx,
    .Ly = Ly,
    .cfl_frac = cfl_frac,
    .spacetime_gauge = spacetime_gauge,
    .reinit_freq = reinit_freq,
    .t_end = t_end,
    .num_frames = num_frames,
    .field_energy_calcs = field_energy_calcs,
    .integrated_mom_calcs = integrated_mom_calcs,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
    .x_loc = x_loc,
  };

  return ctx;
}

void
evalGREulerInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct bhl_static_ctx *app = ctx;

  double gas_gamma = app->gas_gamma;

  double rhol = app->rhol;
  double ul = app->ul;
  double pl = app->pl;

  double rhor = app->rhor;
  double ur = app->ur;
  double pr = app->pr;

  struct gkyl_gr_spacetime *spacetime = app->spacetime;

  double x_loc = app->x_loc;

  double Lx = app->Lx;
  double Ly = app->Ly;

  double rho = 0.0;
  double u = 0.0;
  double p = 0.0;

  if (x < x_loc) {
    rho = rhol; // Fluid mass density (left).
    u = ul; // Fluid velocity (left).
    p = pl; // Fluid pressure (left).
  }
  else {
    rho = rhor; // Fluid mass density (right).
    u = ur; // Fluid velocity (right).
    p = pr; // Fluid pressure (right).
  }
  
  double spatial_det, lapse;
  double *shift = gkyl_malloc(sizeof(double[3]));
  bool in_excision_region;

  double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
  }

  double *lapse_der = gkyl_malloc(sizeof(double[3]));
  double **shift_der = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    shift_der[i] = gkyl_malloc(sizeof(double[3]));
  }

  double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));

    for (int j = 0; j < 3; j++) {
      spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
    }
  }

  spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, 0.0, &spatial_det);
  spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse);
  spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &shift);
  spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);
  
  spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spatial_metric);
  spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature);

  spacetime->lapse_function_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der);
  spacetime->shift_vector_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der);
  spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der);

  double *vel = gkyl_malloc(sizeof(double[3]));
  double v_sq = 0.0;
  vel[0] = u; vel[1] = 0.0; vel[2] = 0.0;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      v_sq += spatial_metric[i][j] * vel[i] * vel[j];
    }
  }

  double W = 1.0 / (sqrt(1.0 - v_sq));
  if (v_sq > 1.0 - pow(10.0, -8.0)) {
    W = 1.0 / sqrt(1.0 - pow(10.0, -8.0));
  }

  double h = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0)));
  
  double rho_rel = sqrt(spatial_det) * rho * W; // Fluid relativistic mass density.
  double mom_x = sqrt(spatial_det) * rho * h * (W * W) * u; // Fluid momentum density (x-direction).
  double mom_y = 0.0; // Fluid momentum density (y-direction).
  double mom_z = 0.0; // Fluid momentum density (z-direction).
  double Etot = sqrt(spatial_det) * ((rho * h * (W * W)) - p - (rho * W)); // Fluid total energy density.

  // Set fluid relativistic mass density.
  fout[0] = rho_rel;
  // Set fluid momentum density.
  fout[1] = mom_x; fout[2] = mom_y; fout[3] = mom_z;
  // Set fluid total energy density.
  fout[4] = Etot;

  // Set lapse gauge variable.
  fout[5] = lapse;
  // Set shift gauge variables.
  fout[6] = shift[0]; fout[7] = shift[1]; fout[8] = shift[2];

  // Set spatial metric tensor.
  fout[9] = spatial_metric[0][0]; fout[10] = spatial_metric[0][1]; fout[11] = spatial_metric[0][2];
  fout[12] = spatial_metric[1][0]; fout[13] = spatial_metric[1][1]; fout[14] = spatial_metric[1][2];
  fout[15] = spatial_metric[2][0]; fout[16] = spatial_metric[2][1]; fout[17] = spatial_metric[2][2];

  // Set extrinsic curvature tensor.
  fout[18] = extrinsic_curvature[0][0]; fout[19] = extrinsic_curvature[0][1]; fout[20] = extrinsic_curvature[0][2];
  fout[21] = extrinsic_curvature[1][0]; fout[22] = extrinsic_curvature[1][1]; fout[23] = extrinsic_curvature[1][2];
  fout[24] = extrinsic_curvature[2][0]; fout[25] = extrinsic_curvature[2][1]; fout[26] = extrinsic_curvature[2][2];

  // Set excision boundary conditions.
  if (in_excision_region) {
    fout[27] = -1.0;
  }
  else {
    fout[27] = 1.0;
  }

  // Set lapse function derivatives.
  fout[28] = lapse_der[0]; fout[29] = lapse_der[1]; fout[30] = lapse_der[2];
  // Set shift vector derivatives.
  fout[31] = shift_der[0][0]; fout[32] = shift_der[0][1]; fout[33] = shift_der[0][2];
  fout[34] = shift_der[1][0]; fout[35] = shift_der[1][1]; fout[36] = shift_der[1][2];
  fout[37] = shift_der[2][0]; fout[38] = shift_der[2][1]; fout[39] = shift_der[2][2];

  // Set spatial metric tensor derivatives.
  fout[40] = spatial_metric_der[0][0][0]; fout[41] = spatial_metric_der[0][0][1]; fout[42] = spatial_metric_der[0][0][2];
  fout[43] = spatial_metric_der[0][1][0]; fout[44] = spatial_metric_der[0][1][1]; fout[45] = spatial_metric_der[0][1][2];
  fout[46] = spatial_metric_der[0][2][0]; fout[47] = spatial_metric_der[0][2][1]; fout[48] = spatial_metric_der[0][2][2];

  fout[49] = spatial_metric_der[1][0][0]; fout[50] = spatial_metric_der[1][0][1]; fout[51] = spatial_metric_der[1][0][2];
  fout[52] = spatial_metric_der[1][1][0]; fout[53] = spatial_metric_der[1][1][1]; fout[54] = spatial_metric_der[1][1][2];
  fout[55] = spatial_metric_der[1][2][0]; fout[56] = spatial_metric_der[1][2][1]; fout[57] = spatial_metric_der[1][2][2];

  fout[58] = spatial_metric_der[2][0][0]; fout[59] = spatial_metric_der[2][0][1]; fout[60] = spatial_metric_der[2][0][2];
  fout[61] = spatial_metric_der[2][1][0]; fout[62] = spatial_metric_der[2][1][1]; fout[63] = spatial_metric_der[2][1][2];
  fout[64] = spatial_metric_der[2][2][0]; fout[65] = spatial_metric_der[2][2][1]; fout[66] = spatial_metric_der[2][2][2];

  // Set evolution parameter.
  fout[67] = 0.0;

  // Set spatial coordinates.
  fout[68] = x; fout[69] = y; fout[70] = 0.0;

  if (in_excision_region) {
    for (int i = 0; i < 67; i++) {
      fout[i] = 0.0;
    }
    fout[27] = -1.0;
  }

  // Free all tensorial quantities.
  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_metric[i]);
    gkyl_free(extrinsic_curvature[i]);
    gkyl_free(shift_der[i]);

    for (int j = 0; j < 3; j++) {
      gkyl_free(spatial_metric_der[i][j]);
    }
    gkyl_free(spatial_metric_der[i]);
  }
  gkyl_free(spatial_metric);
  gkyl_free(extrinsic_curvature);
  gkyl_free(shift);
  gkyl_free(vel);
  gkyl_free(lapse_der);
  gkyl_free(shift_der);
  gkyl_free(spatial_metric_der);
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_moment_app* app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_write) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_moment_app_write(app, t_curr, frame);
    gkyl_moment_app_write_field_energy(app);
    gkyl_moment_app_write_integrated_mom(app);
  }
}

void
calc_field_energy(struct gkyl_tm_trigger* fet, gkyl_moment_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(fet, t_curr) || force_calc) {
    gkyl_moment_app_calc_field_energy(app, t_curr);
  }
}

void
calc_integrated_mom(struct gkyl_tm_trigger* imt, gkyl_moment_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(imt, t_curr) || force_calc) {
    gkyl_moment_app_calc_integrated_mom(app, t_curr);
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

  struct bhl_static_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], ctx.Ny);

  // Fluid equations.
  struct gkyl_wv_eqn *gr_euler = gkyl_wv_gr_euler_new(ctx.gas_gamma, ctx.spacetime_gauge, ctx.reinit_freq, ctx.spacetime, app_args.use_gpu);

  struct gkyl_moment_species fluid = {
    .name = "gr_euler",
    .equation = gr_euler,
    
    .init = evalGREulerInit,
    .force_low_order_flux = true, // Use Lax fluxes.
    .ctx = &ctx,

    .has_gr_euler = true,
    .gr_euler_gas_gamma = ctx.gas_gamma,

    .bcx = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
    .bcy = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
  };

  int nrank = 1; // Number of processes in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif

  // Create global range.
  int cells[] = { NX, NY };
  int dim = sizeof(cells) / sizeof(cells[0]);

  int cuts[dim];
#ifdef GKYL_HAVE_MPI
  for (int d = 0; d < dim; d++) {
    if (app_args.use_mpi) {
      cuts[d] = app_args.cuts[d];
    }
    else {
      cuts[d] = 1;
    }
  }
#else
  for (int d = 0; d < dim; d++) {
    cuts[d] = 1;
  }
#endif

  // Construct communicator for use in app.
  struct gkyl_comm *comm;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
      }
    );
  }
  else {
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .use_gpu = app_args.use_gpu
      }
    );
  }
#else
  comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .use_gpu = app_args.use_gpu
    }
  );
#endif

  int my_rank;
  gkyl_comm_get_rank(comm, &my_rank);
  int comm_size;
  gkyl_comm_get_size(comm, &comm_size);

  int ncuts = 1;
  for (int d = 0; d < dim; d++) {
    ncuts *= cuts[d];
  }

  if (ncuts != comm_size) {
    if (my_rank == 0) {
      fprintf(stderr, "*** Number of ranks, %d, does not match total cuts, %d!\n", comm_size, ncuts);
    }
    goto mpifinalize;
  }

  // Moment app.
  struct gkyl_moment app_inp = {

    .ndim = 2,
    .lower = { 0.0, 0.0 },
    .upper = { ctx.Lx, ctx.Ly },
    .cells = { NX, NY },

    .scheme_type = GKYL_MOMENT_WAVE_PROP,
    .mp_recon = app_args.mp_recon,

    .cfl_frac = ctx.cfl_frac,

    .num_species = 1,
    .species = { fluid },

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0], app_args.cuts[1] },
      .comm = comm,
    },
  };

  // Create app object.
  // Set app output name from the executable name (argv[0]).
  snprintf(app_inp.name, sizeof(app_inp.name), "%s", app_args.app_name);
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  // Initial and final simulation times.
  double t_curr = 0.0, t_end = ctx.t_end;

  // Initialize simulation.
  int frame_curr = 0;
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_moment_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_moment_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n", gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_moment_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_moment_app_cout(app, stdout, " at time = %g\n", t_curr);
  }
  else {
    gkyl_moment_app_apply_ic(app, t_curr);
  }

  // Create trigger for field energy.
  int field_energy_calcs = ctx.field_energy_calcs;
  struct gkyl_tm_trigger fe_trig = { .dt = t_end / field_energy_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_field_energy(&fe_trig, app, t_curr, false);

  // Create trigger for integrated moments.
  int integrated_mom_calcs = ctx.integrated_mom_calcs;
  struct gkyl_tm_trigger im_trig = { .dt = t_end / integrated_mom_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_integrated_mom(&im_trig, app, t_curr, false);

  // Create trigger for IO.
  int num_frames = ctx.num_frames;
  struct gkyl_tm_trigger io_trig = { .dt = t_end / num_frames, .tcurr = frame_curr * (t_end / num_frames), .curr = frame_curr };

  write_data(&io_trig, app, t_curr, false);

  // Compute initial guess of maximum stable time-step.
  double dt = t_end - t_curr;

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = ctx.dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx.num_failures_max;

  long step = 1;
  while ((t_curr < t_end) && (step <= app_args.num_steps)) {
    gkyl_moment_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = gkyl_moment_update(app, dt);
    gkyl_moment_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      gkyl_moment_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_field_energy(&fe_trig, app, t_curr, false);
    calc_integrated_mom(&im_trig, app, t_curr, false);
    write_data(&io_trig, app, t_curr, false);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_moment_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_moment_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_moment_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_moment_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_moment_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);

        calc_field_energy(&fe_trig, app, t_curr, true);
        calc_integrated_mom(&im_trig, app, t_curr, true);
        write_data(&io_trig, app, t_curr, true);

        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }

  calc_field_energy(&fe_trig, app, t_curr, false);
  calc_integrated_mom(&im_trig, app, t_curr, false);
  write_data(&io_trig, app, t_curr, false);
  gkyl_moment_app_stat_write(app);

  struct gkyl_moment_stat stat = gkyl_moment_app_stat(app);

  gkyl_moment_app_cout(app, stdout, "\n");
  gkyl_moment_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_moment_app_cout(app, stdout, "Number of failed time-steps %ld\n", stat.nfail);
  gkyl_moment_app_cout(app, stdout, "Species updates took %g secs\n", stat.species_tm);
  gkyl_moment_app_cout(app, stdout, "Field updates took %g secs\n", stat.field_tm);
  gkyl_moment_app_cout(app, stdout, "Source updates took %g secs\n", stat.sources_tm);
  gkyl_moment_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

freeresources:
  // Free resources after simulation completion.
  gkyl_wv_eqn_release(gr_euler);
  gkyl_gr_spacetime_release(ctx.spacetime);
  gkyl_comm_release(comm);
  gkyl_moment_app_release(app);  
  
mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif
  
  return 0;
}
