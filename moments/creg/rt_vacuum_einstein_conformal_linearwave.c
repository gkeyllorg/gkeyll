#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_vacuum_einstein_conformal.h>
#include <gkyl_gr_minkowski.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

#include <rt_arg_parse.h>

struct einstein_conformal_linearwave_ctx
{
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using normalized code units).
  double amp; // Wave amplitude.

  // Pointer to spacetime metric.
  struct gkyl_gr_spacetime *spacetime;

  // Evolution parameters.
  double excision_threshold; // Excision threshold (lapse).
  enum gkyl_spacetime_slicing spacetime_slicing; // Spacetime slicing condition.
  enum gkyl_spacetime_evolution spacetime_evolution; // Spacetime evolution system.

  // Simulation parameters.
  int Nx; // Cell count (x-direction).
  double Lx; // Domain size (x-direction).
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  int field_energy_calcs; // Number of times to calculate field energy.
  int integrated_mom_calcs; // Number of times to calculate integrated moments.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct einstein_conformal_linearwave_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  double amp = pow(10.0, -8.0); // Wave amplitude.

  // Pointer to spacetime metric.
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);

  // Evolution parameters.
  double excision_threshold = 0.3; // Excision threshold (lapse).
  enum gkyl_spacetime_slicing spacetime_slicing = GKYL_HARMONIC_SLICING; // Spacetime slicing condition.
  enum gkyl_spacetime_evolution spacetime_evolution = GKYL_EINSTEIN_EVOLUTION; // Spacetime evolution system.

  // Simulation parameters.
  int Nx = 50; // Cell count (x-direction).
  double Lx = 1.0; // Domain size (x-direction).
  double cfl_frac = 0.95; // CFL coefficient.

  double t_end = 1000.0; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct einstein_conformal_linearwave_ctx ctx = {
    .pi = pi,
    .amp = amp,
    .spacetime = spacetime,
    .excision_threshold = excision_threshold,
    .spacetime_slicing = spacetime_slicing,
    .spacetime_evolution = spacetime_evolution,
    .Nx = Nx,
    .Lx = Lx,
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
evalVacuumEinsteinConformalInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0];
  struct einstein_conformal_linearwave_ctx *app = ctx;

  double pi = app->pi;
  double amp = app->amp;
  struct gkyl_gr_spacetime *spacetime = app->spacetime;

  double conformal_spatial_det, conformal_lapse;
  double *conformal_shift = gkyl_malloc(sizeof(double[3]));
  bool in_excision_region;

  double **conformal_spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    conformal_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  double **inv_conformal_spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_conformal_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  double **conformal_extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    conformal_extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
  }

  double *conformal_lapse_der = gkyl_malloc(sizeof(double[3]));
  double **conformal_shift_der = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    conformal_shift_der[i] = gkyl_malloc(sizeof(double[3]));
  }

  double ***conformal_spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
  for (int i = 0; i < 3; i++) {
    conformal_spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));

    for (int j = 0; j < 3; j++) {
      conformal_spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
    }
  }

  double conformal_fact, bssn_conformal_fact;
  double *conformal_fact_der = gkyl_malloc(sizeof(double[3]));
  double *bssn_conformal_fact_der = gkyl_malloc(sizeof(double[3]));
  double **bssn_conformal_fact_der2 = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    bssn_conformal_fact_der2[i] = gkyl_malloc(sizeof(double[3]));
  }

  spacetime->spatial_metric_det_func(spacetime, 0.0, x, 0.0, 0.0, &conformal_spatial_det);
  spacetime->lapse_function_func(spacetime, 0.0, x, 0.0, 0.0, &conformal_lapse);
  spacetime->shift_vector_func(spacetime, 0.0, x, 0.0, 0.0, &conformal_shift);
  spacetime->excision_region_func(spacetime, 0.0, x, 0.0, 0.0, &in_excision_region);
  
  spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, 0.0, 0.0, &conformal_spatial_metric);
  spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, 0.0, 0.0, 1.0, 1.0, 1.0, &conformal_extrinsic_curvature);

  spacetime->conformal_factor_func(spacetime, 0.0, x, 0.0, 0.0, &conformal_fact);
  spacetime->bssn_conformal_factor_func(spacetime, 0.0, x, 0.0, 0.0, &bssn_conformal_fact);

  spacetime->conformal_factor_der_func(spacetime, 0.0, x, 0.0, 0.0, 1.0, 1.0, 1.0, &conformal_fact_der);
  spacetime->bssn_conformal_factor_der_func(spacetime, 0.0, x, 0.0, 0.0, 1.0, 1.0, 1.0, &bssn_conformal_fact_der);
  spacetime->bssn_conformal_factor_der2_func(spacetime, 0.0, x, 0.0, 0.0, 1.0, 1.0, 1.0, &bssn_conformal_fact_der2);

  spacetime->lapse_function_der_func(spacetime, 0.0, x, 0.0, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &conformal_lapse_der);
  spacetime->shift_vector_der_func(spacetime, 0.0, x, 0.0, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &conformal_shift_der);
  spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, 0.0, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &conformal_spatial_metric_der);

  double b = amp * sin(2.0 * pi * x);
  conformal_spatial_det = 1.0 - (b * b);
  conformal_fact = pow(conformal_spatial_det, 1.0 / 12.0);
  bssn_conformal_fact = 1.0 / (conformal_fact * conformal_fact);

  conformal_spatial_metric[1][1] = (1.0 + b) / (conformal_fact * conformal_fact * conformal_fact * conformal_fact);
  conformal_spatial_metric[2][2] = (1.0 - b) / (conformal_fact * conformal_fact * conformal_fact * conformal_fact);

  conformal_extrinsic_curvature[1][1] = -amp * pi * cos(2.0 * pi * x);
  conformal_extrinsic_curvature[2][2] = amp * pi * cos(2.0 * pi * x);

  conformal_spatial_metric_der[0][1][1] = 2.0 * amp * pi * cos(2.0 * pi * x);
  conformal_spatial_metric_der[0][2][2] = -2.0 * amp * pi * cos(2.0 * pi * x);

  inv_conformal_spatial_metric[1][1] = (1.0 / (1.0 + b)) * (conformal_fact * conformal_fact * conformal_fact * conformal_fact);;
  inv_conformal_spatial_metric[2][2] = (1.0 / (1.0 + b)) * (conformal_fact * conformal_fact * conformal_fact * conformal_fact);;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        conformal_spatial_metric_der[i][j][k] = 0.5 * conformal_spatial_metric_der[i][j][k];
        conformal_spatial_metric_der[i][j][k] /= (conformal_fact * conformal_fact * conformal_fact * conformal_fact);
        conformal_spatial_metric_der[i][j][k] -= 2.0 * conformal_fact_der[i] * conformal_spatial_metric[j][k];
      }

      conformal_shift_der[i][j] = 0.5 * conformal_shift_der[i][j];
    }
  }

  for (int i = 0; i < 3; i++) {
    conformal_lapse_der[i] = conformal_lapse_der[i] / conformal_lapse;
  }

  double conformal_spatial_metric_der_raised1[3][3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        conformal_spatial_metric_der_raised1[k][i][j] = 0.0;
        
        for (int l = 0; l < 3; l++) {
          conformal_spatial_metric_der_raised1[k][i][j] += inv_conformal_spatial_metric[k][l] * conformal_spatial_metric_der[l][i][j];
        }
      }
    }
  }
  
  double conformal_spatial_metric_der_raised3[3][3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        conformal_spatial_metric_der_raised3[i][j][k] = 0.0;

        for (int l = 0; l < 3; l++) {
          conformal_spatial_metric_der_raised3[i][j][k] += inv_conformal_spatial_metric[l][k] * conformal_spatial_metric_der[i][j][l];
        }
      }
    }
  }

  double conformal_aux_vect[3];
  for (int i = 0; i < 3; i++) {
    conformal_aux_vect[i] = 0.0;

    for (int s = 0; s < 3; s++) {
      conformal_aux_vect[i] += conformal_spatial_metric_der_raised3[i][s][s];
      conformal_aux_vect[i] -= conformal_spatial_metric_der_raised1[s][s][i];  
    }

    conformal_aux_vect[i] -= 4.0 * conformal_fact_der[i];
  }

  // Set conformal spatial metric tensor.
  fout[0] = conformal_spatial_metric[0][0]; fout[1] = conformal_spatial_metric[0][1]; fout[2] = conformal_spatial_metric[0][2];
  fout[3] = conformal_spatial_metric[1][0]; fout[4] = conformal_spatial_metric[1][1]; fout[5] = conformal_spatial_metric[1][2];
  fout[6] = conformal_spatial_metric[2][0]; fout[7] = conformal_spatial_metric[2][1]; fout[8] = conformal_spatial_metric[2][2];

  // Set conformal lapse gauge variable.
  fout[9] = conformal_lapse;

  // Set conformal extrinsic curvature tensor.
  fout[10] = conformal_extrinsic_curvature[0][0]; fout[11] = conformal_extrinsic_curvature[0][1]; fout[12] = conformal_extrinsic_curvature[0][2];
  fout[13] = conformal_extrinsic_curvature[1][0]; fout[14] = conformal_extrinsic_curvature[1][1]; fout[15] = conformal_extrinsic_curvature[1][2];
  fout[16] = conformal_extrinsic_curvature[2][0]; fout[17] = conformal_extrinsic_curvature[2][1]; fout[18] = conformal_extrinsic_curvature[2][2];

  // Set conformal spatial metric tensor derivatives.
  fout[19] = conformal_spatial_metric_der[0][0][0]; fout[20] = conformal_spatial_metric_der[0][0][1]; fout[21] = conformal_spatial_metric_der[0][0][2];
  fout[22] = conformal_spatial_metric_der[0][1][0]; fout[23] = conformal_spatial_metric_der[0][1][1]; fout[24] = conformal_spatial_metric_der[0][1][2];
  fout[25] = conformal_spatial_metric_der[0][2][0]; fout[26] = conformal_spatial_metric_der[0][2][1]; fout[27] = conformal_spatial_metric_der[0][2][2];

  fout[28] = conformal_spatial_metric_der[1][0][0]; fout[29] = conformal_spatial_metric_der[1][0][1]; fout[30] = conformal_spatial_metric_der[1][0][2];
  fout[31] = conformal_spatial_metric_der[1][1][0]; fout[32] = conformal_spatial_metric_der[1][1][1]; fout[33] = conformal_spatial_metric_der[1][1][2];
  fout[34] = conformal_spatial_metric_der[1][2][0]; fout[35] = conformal_spatial_metric_der[1][2][1]; fout[36] = conformal_spatial_metric_der[1][2][2];

  fout[37] = conformal_spatial_metric_der[2][0][0]; fout[38] = conformal_spatial_metric_der[2][0][1]; fout[39] = conformal_spatial_metric_der[2][0][2];
  fout[40] = conformal_spatial_metric_der[2][1][0]; fout[41] = conformal_spatial_metric_der[2][1][1]; fout[42] = conformal_spatial_metric_der[2][1][2];
  fout[43] = conformal_spatial_metric_der[2][2][0]; fout[44] = conformal_spatial_metric_der[2][2][1]; fout[45] = conformal_spatial_metric_der[2][2][2];

  // Set conformal lapse function derivatives.
  fout[46] = conformal_lapse_der[0]; fout[47] = conformal_lapse_der[1]; fout[48] = conformal_lapse_der[2];

  // Set conformal auxiliary vector.
  fout[49] = conformal_aux_vect[0]; fout[50] = conformal_aux_vect[1]; fout[51] = conformal_aux_vect[2];

  // Set conformal shift gauge variables.
  fout[52] = conformal_shift[0]; fout[53] = conformal_shift[1]; fout[54] = conformal_shift[2];

  // Set conformal shift vector derivatives.
  fout[55] = conformal_shift_der[0][0]; fout[56] = conformal_shift_der[0][1]; fout[57] = conformal_shift_der[0][2];
  fout[58] = conformal_shift_der[1][0]; fout[59] = conformal_shift_der[1][1]; fout[60] = conformal_shift_der[1][2];
  fout[61] = conformal_shift_der[2][0]; fout[62] = conformal_shift_der[2][1]; fout[63] = conformal_shift_der[2][2];

  // Set BSSN conformal factor.
  fout[64] = bssn_conformal_fact;

  // Set BSSN conformal factor first derivatives.
  fout[65] = bssn_conformal_fact_der[0]; fout[66] = bssn_conformal_fact_der[1]; fout[67] = bssn_conformal_fact_der[2];

  // Set BSSN conformal factor second derivatives.
  fout[68] = bssn_conformal_fact_der2[0][0]; fout[69] = bssn_conformal_fact_der2[0][1]; fout[70] = bssn_conformal_fact_der2[0][2];
  fout[71] = bssn_conformal_fact_der2[1][0]; fout[72] = bssn_conformal_fact_der2[1][1]; fout[73] = bssn_conformal_fact_der2[1][2];
  fout[74] = bssn_conformal_fact_der2[2][0]; fout[75] = bssn_conformal_fact_der2[2][1]; fout[76] = bssn_conformal_fact_der2[2][2];

  if (in_excision_region) {
    for (int i = 0; i < 77; i++) {
      fout[i] = 0.0;
    }
  }

  // Free all tensorial quantities.
  for (int i = 0; i < 3; i++) {
    gkyl_free(conformal_spatial_metric[i]);
    gkyl_free(inv_conformal_spatial_metric[i]);
    gkyl_free(conformal_extrinsic_curvature[i]);
    gkyl_free(conformal_shift_der[i]);
    gkyl_free(bssn_conformal_fact_der2[i]);

    for (int j = 0; j < 3; j++) {
      gkyl_free(conformal_spatial_metric_der[i][j]);
    }
    gkyl_free(conformal_spatial_metric_der[i]);
  }
  gkyl_free(conformal_spatial_metric);
  gkyl_free(inv_conformal_spatial_metric);
  gkyl_free(conformal_extrinsic_curvature);
  gkyl_free(conformal_shift);
  gkyl_free(conformal_lapse_der);
  gkyl_free(conformal_shift_der);
  gkyl_free(conformal_spatial_metric_der);
  gkyl_free(conformal_fact_der);
  gkyl_free(bssn_conformal_fact_der);
  gkyl_free(bssn_conformal_fact_der2);
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

  struct einstein_conformal_linearwave_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);

  // Conformal Einstein equations.
  struct gkyl_wv_eqn *vacuum_einstein_conformal = gkyl_wv_vacuum_einstein_conformal_new(ctx.excision_threshold, ctx.spacetime_slicing, ctx.spacetime_evolution, app_args.use_gpu);

  struct gkyl_moment_species einstein_conformal = {
    .name = "vacuum_einstein_conformal",
    .equation = vacuum_einstein_conformal,
    
    .init = evalVacuumEinsteinConformalInit,
    .force_low_order_flux = false, // Use HLL fluxes.
    .ctx = &ctx,

    .has_vacuum_einstein_conformal = true,
    .vacuum_einstein_conformal_excision_threshold = ctx.excision_threshold,
    .vacuum_einstein_conformal_spacetime_slicing = ctx.spacetime_slicing,
    .vacuum_einstein_conformal_spacetime_evolution = ctx.spacetime_evolution,
  };

  int nrank = 1; // Number of processes in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif

  // Create global range.
  int cells[] = { NX };
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
    .name = "vacuum_einstein_conformal_linearwave",

    .ndim = 1,
    .lower = { -0.5 * ctx.Lx },
    .upper = { 0.5 * ctx.Lx },
    .cells = { NX },

    .scheme_type = GKYL_MOMENT_WAVE_PROP,
    .mp_recon = app_args.mp_recon,

    .cfl_frac = ctx.cfl_frac,

    .num_species = 1,
    .species = { einstein_conformal },

    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0] },
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
  gkyl_wv_eqn_release(vacuum_einstein_conformal);
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
