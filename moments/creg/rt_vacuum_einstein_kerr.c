#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_vacuum_einstein.h>
#include <gkyl_gr_minkowski.h>
#include <gkyl_gr_blackhole.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

#include <rt_arg_parse.h>

struct einstein_kerr_ctx
{
  // Spacetime parameters (using geometric units).
  double mass; // Mass of the black hole.
  double spin; // Spin of the black hole.

  double pos_x; // Position of the black hole (x-direction).
  double pos_y; // Position of the black hole (y-direction).
  double pos_z; // Position of the black hole (z-direction).

  // Pointer to spacetime metric.
  struct gkyl_gr_spacetime *spacetime;

  // Evolution parameters.
  double excision_threshold; // Excision threshold (lapse).
  enum gkyl_spacetime_slicing spacetime_slicing; // Spacetime slicing condition.
  enum gkyl_spacetime_evolution spacetime_evolution; // Spacetime evolution system.

  // Simulation parameters.
  int Nx; // Cell count (x-direction).
  int Ny; // Cell count (y-direction).
  double Lx; // Domain size (x-direction).
  double Ly; // Domain size (y-direction).
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  int field_energy_calcs; // Number of times to calculate field energy.
  int integrated_mom_calcs; // Number of times to calculate integrated moments.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct einstein_kerr_ctx
create_ctx(void)
{
  // Spacetime parameters (using geometric units).
  double mass = 0.5; // Mass of the black hole.
  double spin = -0.3; // Spin of the black hole.

  double pos_x = 5.0; // Position of the black hole (x-direction).
  double pos_y = 5.0; // Position of the black hole (y-direction).
  double pos_z = 0.0; // Position of the black hole (z-direction).

  // Pointer to spacetime metric.
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_new(false, mass, spin, pos_x, pos_y, pos_z);

  // Evolution parameters.
  double excision_threshold = 0.3; // Excision threshold (lapse).
  enum gkyl_spacetime_slicing spacetime_slicing = GKYL_1PLUSLOG_SLICING; // Spacetime slicing condition.
  enum gkyl_spacetime_evolution spacetime_evolution = GKYL_EINSTEIN_EVOLUTION; // Spacetime evolution system.

  // Simulation parameters.
  int Nx = 256; // Cell count (x-direction).
  int Ny = 256; // Cell count (y-direction).
  double Lx = 10.0; // Domain size (x-direction).
  double Ly = 10.0; // Domain size (y-direction).
  double cfl_frac = 0.8; // CFL coefficient.

  double t_end = 10.0; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct einstein_kerr_ctx ctx = {
    .mass = mass,
    .spin = spin,
    .pos_x = pos_x,
    .pos_y = pos_y,
    .pos_z = pos_z,
    .spacetime = spacetime,
    .excision_threshold = excision_threshold,
    .spacetime_slicing = spacetime_slicing,
    .spacetime_evolution = spacetime_evolution,
    .Nx = Nx,
    .Ny = Ny,
    .Lx = Lx,
    .Ly = Ly,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .num_frames = num_frames,
    .field_energy_calcs = field_energy_calcs,
    .integrated_mom_calcs = integrated_mom_calcs,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };

  return ctx;
}

void
evalVacuumEinsteinInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct einstein_kerr_ctx *app = ctx;

  struct gkyl_gr_spacetime *spacetime = app->spacetime;

  double spatial_det, lapse;
  double *shift = gkyl_malloc(sizeof(double[3]));
  bool in_excision_region;

  double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
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
  spacetime->spatial_inv_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &inv_spatial_metric);
  spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature);

  spacetime->lapse_function_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der);
  spacetime->shift_vector_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der);
  spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        spatial_metric_der[i][j][k] = 0.5 * spatial_metric_der[i][j][k];
      }

      shift_der[i][j] = 0.5 * shift_der[i][j];
    }
  }

  for (int i = 0; i < 3; i++) {
    lapse_der[i] = lapse_der[i] / lapse;
  }

  double spatial_metric_der_raised1[3][3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        spatial_metric_der_raised1[k][i][j] = 0.0;
        
        for (int l = 0; l < 3; l++) {
          spatial_metric_der_raised1[k][i][j] += inv_spatial_metric[k][l] * spatial_metric_der[l][i][j];
        }
      }
    }
  }
  
  double spatial_metric_der_raised3[3][3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        spatial_metric_der_raised3[i][j][k] = 0.0;

        for (int l = 0; l < 3; l++) {
          spatial_metric_der_raised3[i][j][k] += inv_spatial_metric[l][k] * spatial_metric_der[i][j][l];
        }
      }
    }
  }

  double aux_vect[3];
  for (int i = 0; i < 3; i++) {
    aux_vect[i] = 0.0;

    for (int s = 0; s < 3; s++) {
      aux_vect[i] += spatial_metric_der_raised3[i][s][s];
      aux_vect[i] -= spatial_metric_der_raised1[s][s][i];  
    }
  }

  // Set spatial metric tensor.
  fout[0] = spatial_metric[0][0]; fout[1] = spatial_metric[0][1]; fout[2] = spatial_metric[0][2];
  fout[3] = spatial_metric[1][0]; fout[4] = spatial_metric[1][1]; fout[5] = spatial_metric[1][2];
  fout[6] = spatial_metric[2][0]; fout[7] = spatial_metric[2][1]; fout[8] = spatial_metric[2][2];

  // Set lapse gauge variable.
  fout[9] = lapse;

  // Set extrinsic curvature tensor.
  fout[10] = extrinsic_curvature[0][0]; fout[11] = extrinsic_curvature[0][1]; fout[12] = extrinsic_curvature[0][2];
  fout[13] = extrinsic_curvature[1][0]; fout[14] = extrinsic_curvature[1][1]; fout[15] = extrinsic_curvature[1][2];
  fout[16] = extrinsic_curvature[2][0]; fout[17] = extrinsic_curvature[2][1]; fout[18] = extrinsic_curvature[2][2];

  // Set spatial metric tensor derivatives.
  fout[19] = spatial_metric_der[0][0][0]; fout[20] = spatial_metric_der[0][0][1]; fout[21] = spatial_metric_der[0][0][2];
  fout[22] = spatial_metric_der[0][1][0]; fout[23] = spatial_metric_der[0][1][1]; fout[24] = spatial_metric_der[0][1][2];
  fout[25] = spatial_metric_der[0][2][0]; fout[26] = spatial_metric_der[0][2][1]; fout[27] = spatial_metric_der[0][2][2];

  fout[28] = spatial_metric_der[1][0][0]; fout[29] = spatial_metric_der[1][0][1]; fout[30] = spatial_metric_der[1][0][2];
  fout[31] = spatial_metric_der[1][1][0]; fout[32] = spatial_metric_der[1][1][1]; fout[33] = spatial_metric_der[1][1][2];
  fout[34] = spatial_metric_der[1][2][0]; fout[35] = spatial_metric_der[1][2][1]; fout[36] = spatial_metric_der[1][2][2];

  fout[37] = spatial_metric_der[2][0][0]; fout[38] = spatial_metric_der[2][0][1]; fout[39] = spatial_metric_der[2][0][2];
  fout[40] = spatial_metric_der[2][1][0]; fout[41] = spatial_metric_der[2][1][1]; fout[42] = spatial_metric_der[2][1][2];
  fout[43] = spatial_metric_der[2][2][0]; fout[44] = spatial_metric_der[2][2][1]; fout[45] = spatial_metric_der[2][2][2];

  // Set lapse function derivatives.
  fout[46] = lapse_der[0]; fout[47] = lapse_der[1]; fout[48] = lapse_der[2];

  // Set auxiliary vector.
  fout[49] = aux_vect[0]; fout[50] = aux_vect[1]; fout[51] = aux_vect[2];

  // Set shift gauge variables.
  fout[52] = shift[0]; fout[53] = shift[1]; fout[54] = shift[2];

  // Set shift vector derivatives.
  fout[55] = shift_der[0][0]; fout[56] = shift_der[0][1]; fout[57] = shift_der[0][2];
  fout[58] = shift_der[1][0]; fout[59] = shift_der[1][1]; fout[60] = shift_der[1][2];
  fout[61] = shift_der[2][0]; fout[62] = shift_der[2][1]; fout[63] = shift_der[2][2];

  if (in_excision_region) {
    for (int i = 0; i < 64; i++) {
      fout[i] = 0.0;
    }
  }

  // Free all tensorial quantities.
  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_metric[i]);
    gkyl_free(inv_spatial_metric[i]);
    gkyl_free(extrinsic_curvature[i]);
    gkyl_free(shift_der[i]);

    for (int j = 0; j < 3; j++) {
      gkyl_free(spatial_metric_der[i][j]);
    }
    gkyl_free(spatial_metric_der[i]);
  }
  gkyl_free(spatial_metric);
  gkyl_free(inv_spatial_metric);
  gkyl_free(extrinsic_curvature);
  gkyl_free(shift);
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

  struct einstein_kerr_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], ctx.Ny);

  // Einstein equations.
  struct gkyl_wv_eqn *vacuum_einstein = gkyl_wv_vacuum_einstein_new(ctx.excision_threshold, ctx.spacetime_slicing, ctx.spacetime_evolution, app_args.use_gpu);

  struct gkyl_moment_species einstein = {
    .name = "vacuum_einstein",
    .equation = vacuum_einstein,
    
    .init = evalVacuumEinsteinInit,
    .force_low_order_flux = false, // Use HLL fluxes.
    .ctx = &ctx,

    .has_vacuum_einstein = true,
    .vacuum_einstein_excision_threshold = ctx.excision_threshold,
    .vacuum_einstein_spacetime_slicing = ctx.spacetime_slicing,
    .vacuum_einstein_spacetime_evolution = ctx.spacetime_evolution,

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
    .name = "vacuum_einstein_kerr",

    .ndim = 2,
    .lower = { 0.0, 0.0 },
    .upper = { ctx.Lx, ctx.Ly },
    .cells = { NX, NY },

    .scheme_type = GKYL_MOMENT_WAVE_PROP,
    .mp_recon = app_args.mp_recon,

    .cfl_frac = ctx.cfl_frac,

    .num_species = 1,
    .species = { einstein },

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0], app_args.cuts[1] },
      .comm = comm,
    },
  };

  // Create app object.
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
  gkyl_wv_eqn_release(vacuum_einstein);
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
