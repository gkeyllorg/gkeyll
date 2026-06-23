#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>
#include <gkyl_bc_emission.h>
#include <rt_arg_parse.h>

struct sheath_ctx {
  int cdim, vdim; // Dimensionality.
  double epsilon0;
  double mu0;
  double q0;
  double chargeElc; // electron charge
  double massElc; // electron mass
  double chargeIon; // ion charge
  double massIon; // ion mass
  double n0;
  double Te; // electron to ion temperature ratio
  double Ti;
  double vte; // electron thermal velocity
  double vti; // ion thermal velocity
  double lambda_D;
  double Lx; // size of the box
  double Ls;
  double omega_pe;
  double phi;
  double deltahat_ts;
  double Ehat_ts;
  double t1;
  double t2;
  double t3;
  double t4;
  double s;
  double P1_inf;
  double P1_hat;
  double E_hat;
  double W;
  double p;
  int Nx;
  int Nvx;
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  int num_emission_species;
  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  int int_diag_calc_num; // Number of times to compute integrated diagnostics.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

static inline double sq(double x) { return x*x; }

void
evalDistFuncElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct sheath_ctx *app = ctx;
  double x = xn[0], v = xn[1];
  double vt = app->vte;
  double n = app->n0;
  double fv = n/sqrt(2.0*M_PI*sq(vt))*(exp(-sq(v)/(2*sq(vt))));
  fout[0] = fv;
}

void
evalDistFuncElcSource(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct sheath_ctx *app = ctx;
  double x = xn[0], v = xn[1];
  double vt = app->vte;
  double Ls = app->Ls;
  double fv = 1.0/sqrt(2.0*M_PI*sq(vt))*(exp(-sq(v)/(2*sq(vt))));
  if(fabs(x) < Ls) {
    fout[0] = (Ls - fabs(x))/Ls*fv;
  } else {
    fout[0] = 0.0;
  }
}

void
evalDistFuncIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct sheath_ctx *app = ctx;
  double x = xn[0], v = xn[1];
  double vt = app->vti;
  double n = app->n0;
  double fv = n/sqrt(2.0*M_PI*sq(vt))*(exp(-sq(v)/(2*sq(vt))));
  fout[0] = fv;
}

void
evalDistFuncIonSource(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct sheath_ctx *app = ctx;
  double x = xn[0], v = xn[1];
  double vt = app->vti;
  double Ls = app->Ls;
  double fv = 1.0/sqrt(2.0*M_PI*sq(vt))*(exp(-sq(v)/(2*sq(vt))));
  if(fabs(x) < Ls) {
    fout[0] = (Ls - fabs(x))/Ls*fv;
  } else {
    fout[0] = 0.0;
  }
}

void
evalFieldFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct sheath_ctx *app = ctx;
  double x = xn[0];
  
  fout[0] = 0.0; fout[1] = 0.0, fout[2] = 0.0;
  fout[3] = 0.0; fout[4] = 0.0; fout[5] = 0.0;
  fout[6] = 0.0; fout[7] = 0.0;
}

struct sheath_ctx
create_ctx(void)
{
  int cdim = 1, vdim = 1; // Dimensionality.

  double epsilon0 = 8.854e-12;
  double mu0 = 1.257e-6;

  double massElc = 9.109e-31;
  double massIon = 1836.153*massElc;
  double q0 = 1.602e-19;

  double n0 = 1.0e17;
  double Te0 = 10.0*q0;
  double Ti0 = 10.0*q0;

  double vte0 = sqrt(Te0/massElc);
  double vti0 = sqrt(Ti0/massIon);

  double omega_pe = sqrt(n0*q0*q0/(epsilon0*massElc));
  double lambda_D = sqrt(epsilon0*Te0/(n0*q0*q0));

  double Lx = 128.0*lambda_D;
  double Ls = 100.0*lambda_D;

  int Nx = 128;
  int Nvx = 32;
  
  // SEE parameters
  double phi = 4.68;
  double deltahat_ts = 1.885;
  double Ehat_ts = 276.8;
  double t1 = 0.66;
  double t2 = 0.8;
  double t3 = 0.7;
  double t4 = 1.0;
  double s = 1.54;
  double P1_inf = 0.02;
  double P1_hat = 0.496;
  double E_hat = 1.0e-6;
  double W = 60.86;
  double p = 1.0;

  double t_end = 10.0/omega_pe; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct sheath_ctx ctx = {
    .cdim = cdim,  .vdim = vdim,
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .q0 = q0,
    .chargeElc = -q0,
    .massElc = massElc,
    .chargeIon = q0,
    .massIon = massIon,
    .n0 = n0,
    .Te = Te0,
    .Ti = Ti0,
    .vte = vte0,
    .vti = vti0,
    .lambda_D = lambda_D,
    .Lx = Lx,
    .Ls = Ls,
    .omega_pe = omega_pe,
    .phi = phi,
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
    .Nx = Nx,
    .Nvx = Nvx, 
    .cells = {Nx, Nvx},
    .num_emission_species = 1,
    .t_end = t_end,
    .num_frames = num_frames,
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };
  return ctx;
}

void
calc_integrated_diagnostics(struct gkyl_tm_trigger* iot, gkyl_vlasov_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_calc) {
    gkyl_vlasov_app_calc_field_energy(app, t_curr);
    gkyl_vlasov_app_calc_integrated_mom(app, t_curr);
  }
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_vlasov_app* app, double t_curr, bool force_write)
{
  bool trig_now = gkyl_tm_trigger_check_and_bump(iot, t_curr);
  if (trig_now || force_write) {
    int frame = (!trig_now) && force_write? iot->curr : iot->curr-1;

    gkyl_vlasov_app_write(app, t_curr, frame);

    gkyl_vlasov_app_calc_mom(app);
    gkyl_vlasov_app_write_mom(app, t_curr, frame);

    gkyl_vlasov_app_calc_field_energy(app, t_curr);
    gkyl_vlasov_app_write_field_energy(app);

    gkyl_vlasov_app_calc_integrated_mom(app, t_curr);
    gkyl_vlasov_app_write_integrated_mom(app);
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

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_vlasov_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  char in_species[1][128] = { "elc" };
  struct gkyl_bc_emission_ctx *bc_ctx = gkyl_bc_emission_secondary_electron_copper_new(ctx.num_emission_species, 0.0, in_species, app_args.use_gpu);
  /* struct gkyl_spectrum_model *spectrum_model[1]; */
  /* spectrum_model[0] = gkyl_spectrum_chung_everhart_new(ctx.q0, ctx.phi, app_args.use_gpu); */
  /* struct gkyl_yield_model *yield_model[1]; */
  /* yield_model[0] = gkyl_yield_furman_pivi_new(ctx.q0, ctx.deltahat_ts, ctx.Ehat_ts, ctx.t1, ctx.t2, ctx.t3, ctx.t4, ctx.s, app_args.use_gpu); */
  /* struct gkyl_elastic_model *elastic_model = gkyl_elastic_furman_pivi_new(ctx.q0, ctx.P1_inf, ctx.P1_hat, ctx.E_hat, ctx.W, ctx.p, app_args.use_gpu); */
  /* struct gkyl_bc_emission_ctx *bc_ctx = gkyl_bc_emission_new(ctx.num_emission_species, 0.0, true, spectrum_model, yield_model, elastic_model, in_species); */

  // electrons
  struct gkyl_vlasov_species elc = {
    .name = "elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -4.0*ctx.vte},
    .upper = { 4.0*ctx.vte}, 
    .cells = { cells_v[0] },

    .num_init = 1,
    .projection[0] = {
      .proj_id = GKYL_PROJ_FUNC,
      .func = evalDistFuncElc,
      .ctx_func = &ctx,
    },

    .source = {
      .source_id = GKYL_BFLUX_SOURCE,
      .source_length = ctx.Ls,
      .source_species = "ion",
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_FUNC,
        .func = evalDistFuncElcSource,
        .ctx_func = &ctx,
      },
    },

    .bcx = {
      .lower = { .type = GKYL_SPECIES_REFLECT, },
      .upper = { .type = GKYL_SPECIES_EMISSION,
                 .aux_ctx = bc_ctx, },
    },
    
    .num_diag_moments = 3,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2 },
  };

  // ions
  struct gkyl_vlasov_species ion = {
    .name = "ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { -4.0*ctx.vti},
    .upper = { 4.0*ctx.vti}, 
    .cells = { cells_v[0] },

    .num_init = 1,
    .projection[0] = {
      .proj_id = GKYL_PROJ_FUNC,
      .func = evalDistFuncIon,
      .ctx_func = &ctx,
    },

    .source = {
      .source_id = GKYL_BFLUX_SOURCE,
      .source_length = ctx.Ls,
      .source_species = "ion",
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_FUNC,
        .func = evalDistFuncIonSource,
        .ctx_func = &ctx,
      },
    },

    .bcx = {
      .lower = { .type = GKYL_SPECIES_REFLECT, },
      .upper = { .type = GKYL_SPECIES_ABSORB, },
    },
    
    .num_diag_moments = 3,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2 },
  };

  // Field.
  struct gkyl_vlasov_field field = {
    .epsilon0 = ctx.epsilon0,
    .poisson_bcs = {
      .lo_type = { GKYL_POISSON_NEUMANN },
      .up_type = { GKYL_POISSON_DIRICHLET },
      .lo_value = { 0.0 }, .up_value = { 0.0 }
    },
  };

  // VM app
  struct gkyl_vm app_inp = {

    .cdim = ctx.cdim, .vdim = ctx.vdim,
    .lower = { -ctx.Lx },
    .upper = { ctx.Lx },
    .cells = { cells_x[0] },
    .poly_order = 2,
    .basis_type = app_args.basis_type,

    .num_periodic_dir = 0,
    .periodic_dirs = { },

    .num_species = 2,
    .species = { elc, ion },
    .field = field,
    .is_electrostatic = true,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0] },
      .comm = comm,
    }

  };

  // Create app object.
  // Set app output name from the executable name (argv[0]).
  snprintf(app_inp.name, sizeof(app_inp.name), "%s", app_args.app_name);
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&app_inp);

  // Initial & final simulation times.
  int frame_curr = 0;
  double t_curr = 0.0, t_end = ctx.t_end;
  // Initialize simulation.
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_vlasov_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_vlasov_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n",
        gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_vlasov_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_vlasov_app_cout(app, stdout, " at time = %g\n", t_curr);
  }
  else {
    gkyl_vlasov_app_apply_ic(app, t_curr);
  }

  // Create triggers for IO.
  int num_frames = ctx.num_frames, num_int_diag_calc = ctx.int_diag_calc_num;
  struct gkyl_tm_trigger trig_write = { .dt = t_end/num_frames, .tcurr = t_curr, .curr = frame_curr };
  struct gkyl_tm_trigger trig_calc_intdiag = { .dt = t_end/GKYL_MAX2(num_frames, num_int_diag_calc),
    .tcurr = t_curr, .curr = frame_curr };

  // Write out ICs (if restart, it overwrites the restart frame).
  calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, false);
  write_data(&trig_write, app, t_curr, false);

  // Compute initial guess of maximum stable time-step.
  double dt = t_end - t_curr;

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

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, t_curr > t_end);
    write_data(&trig_write, app, t_curr, t_curr > t_end);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_vlasov_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_vlasov_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_vlasov_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_vlasov_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_vlasov_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);
	calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, true);
        write_data(&trig_write, app, t_curr, true);
        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }

  gkyl_vlasov_app_stat_write(app);

  struct gkyl_vlasov_stat stat = gkyl_vlasov_app_stat(app);

  gkyl_vlasov_app_cout(app, stdout, "\n");
  gkyl_vlasov_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_vlasov_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_vlasov_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_vlasov_app_cout(app, stdout, "  Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_vlasov_app_cout(app, stdout, "  Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  gkyl_vlasov_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_vlasov_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_vlasov_app_cout(app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm);
  gkyl_vlasov_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_vlasov_app_cout(app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm);
  gkyl_vlasov_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

  gkyl_vlasov_app_cout(app, stdout, "Number of write calls %ld\n", stat.n_io);
  gkyl_vlasov_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

  freeresources:
  // Free resources after simulation completion.
  gkyl_vlasov_app_release(app);
  gkyl_vlasov_comms_release(comm);
  /* for (int i=0; i<ctx.num_emission_species; ++i) { */
  /*   gkyl_spectrum_model_release(spectrum_model[i]); */
  /*   gkyl_yield_model_release(yield_model[i]); */
  /* } */
  /* gkyl_elastic_model_release(elastic_model); */
  gkyl_bc_emission_release(bc_ctx);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif

  return 0;
}
