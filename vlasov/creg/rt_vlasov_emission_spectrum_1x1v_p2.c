#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_bc_emission.h>
#include <gkyl_vlasov.h>
#include <rt_arg_parse.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#ifdef GKYL_HAVE_NCCL
#include <gkyl_nccl_comm.h>
#endif
#endif

struct sheath_ctx {
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
  int Nv;
  int num_emission_species;
  double t_end;
  int num_frames;
  double dt_failure_tol;
  int num_failures_max;
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
    fout[0] = 2*(Ls - fabs(x))/Ls*fv;
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
    fout[0] = 2*(Ls - fabs(x))/Ls*fv;
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
  double massElc = 9.109e-31;
  double q0 = 1.602e-19;

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
  struct sheath_ctx ctx = {
    .epsilon0 = 8.854e-12,
    .mu0 = 1.257e-6,
    .q0 = q0,
    .chargeElc = -q0,
    .massElc = massElc,
    .chargeIon = q0,
    .massIon = 1836.153*massElc,
    .n0 = 1.0e17,
    .Te = 10.0*q0,
    .Ti = 10.0*q0,
    .vte = sqrt(ctx.Te/massElc),
    .vti = sqrt(ctx.Ti/ctx.massIon),
    .lambda_D = sqrt(ctx.epsilon0*ctx.Te/(ctx.n0*q0*q0)),
    .Lx = 128.0*ctx.lambda_D,
    .Ls = 100.0*ctx.lambda_D,
    .omega_pe = sqrt(ctx.n0*q0*q0/(ctx.epsilon0*massElc)),
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
    .Nx = 128,
    .Nv = 32,
    .num_emission_species = 1,
    .t_end = 10.0/ctx.omega_pe,
    .num_frames = 1,
    .dt_failure_tol = 1.0e-4,
    .num_failures_max = 20,
  };
  return ctx;
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_vlasov_app* app, double t_curr, bool force_write)
{
  gkyl_vlasov_app_calc_integrated_mom(app, t_curr);
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr)) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_vlasov_app_write(app, t_curr, iot->curr - 1);

    gkyl_vlasov_app_calc_mom(app);
    gkyl_vlasov_app_write_mom(app, t_curr, iot->curr - 1);
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
  int NV = APP_ARGS_CHOOSE(app_args.vcells[0], ctx.Nv);

  int nrank = 1; // Number of processors in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif

  int ccells[] = { NX };
  int cdim = sizeof(ccells) / sizeof(ccells[0]);

  int cuts[cdim];
#ifdef GKYL_HAVE_MPI  
  for (int d = 0; d < cdim; d++) {
    if (app_args.use_mpi) {
      cuts[d] = app_args.cuts[d];
    }
    else {
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
    comm = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
      }
    );
#else
    printf(" Using -g and -M together requires NCCL.\n");
    assert(0 == 1);
#endif
  }
  else if (app_args.use_mpi) {
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
  for (int d = 0; d < cdim; d++) {
    ncuts *= cuts[d];
  }

  if (ncuts != comm_size) {
    if (my_rank == 0) {
      fprintf(stderr, "*** Number of ranks, %d, does not match total cuts, %d!\n", comm_size, ncuts);
    }
    goto mpifinalize;
  }

  char in_species[1][128] = { "elc" };

  // Copper preset for emission BC where models and parameters are chosen automatically
  struct gkyl_bc_emission_ctx *bc_ctx =
    gkyl_bc_emission_secondary_electron_copper_new(ctx.num_emission_species, 0.0,
    in_species, app_args.use_gpu);

  // Full specification of the BC with user-defined models and parameters
  // Uncomment and use this if no material preset is available
  /* struct gkyl_emission_spectrum_model *spectrum_model[1]; */
  /* spectrum_model[0] = gkyl_emission_spectrum_chung_everhart_new(ctx.q0, ctx.phi, app_args.use_gpu); */
  /* struct gkyl_emission_yield_model *yield_model[1]; */
  /* yield_model[0] = gkyl_emission_yield_furman_pivi_new(ctx.q0, ctx.deltahat_ts, */
  /*   ctx.Ehat_ts, ctx.t1, ctx.t2, ctx.t3, ctx.t4, ctx.s, app_args.use_gpu); */
  /* struct gkyl_emission_elastic_model *elastic_model = */
  /*   gkyl_emission_elastic_furman_pivi_new(ctx.q0, ctx.P1_inf, ctx.P1_hat, ctx.E_hat, */
  /*   ctx.W, ctx.p, app_args.use_gpu); */
  /* struct gkyl_bc_emission_ctx *bc_ctx = gkyl_bc_emission_new(ctx.num_emission_species, */
  /*   0.0, true, spectrum_model, yield_model, elastic_model, in_species); */

  // electrons
  struct gkyl_vlasov_species elc = {
    .name = "elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -4.0*ctx.vte},
    .upper = { 4.0*ctx.vte}, 
    .cells = { NV },

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
    .cells = { NV },

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

  // field
  struct gkyl_vlasov_field field = {
    .epsilon0 = ctx.epsilon0, .mu0 = ctx.mu0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,

    .ctx = &ctx,
    .init = evalFieldFunc,

    .bcx = { GKYL_FIELD_SYM_WALL, GKYL_FIELD_PEC_WALL }
  };

  // VM app
  struct gkyl_vm app_inp = {

    .cdim = 1, .vdim = 1,
    .lower = { 0.0 },
    .upper = { ctx.Lx },
    .cells = { NX },
    .poly_order = 2,
    .basis_type = app_args.basis_type,

    .num_periodic_dir = 0,
    .periodic_dirs = { },

    .num_species = 2,
    .species = { elc, ion },
    .field = field,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0] },
      .comm = comm,
    },
  };

  // Create app object.
  // Set app output name from the executable name (argv[0]).
  snprintf(app_inp.name, sizeof(app_inp.name), "%s", app_args.app_name);
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&app_inp);

  // Initial and final simulation times.
  double t_curr = 0.0, t_end = ctx.t_end;

  // Create trigger for IO.
  int num_frames = ctx.num_frames;
  struct gkyl_tm_trigger io_trig = { .dt = t_end / num_frames };

  // Initialize simulation.
  gkyl_vlasov_app_apply_ic(app, t_curr);
  write_data(&io_trig, app, t_curr, false);

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

    write_data(&io_trig, app, t_curr, false);

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
        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }

  write_data(&io_trig, app, t_curr, false);
  gkyl_vlasov_app_write_integrated_mom(app);
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

  // Free resources after simulation completion.
  gkyl_comm_release(comm);
  gkyl_vlasov_app_release(app);
  /* for (int i=0; i<ctx.num_emission_species; ++i) { */
  /*   gkyl_spectrum_model_release(spectrum_model[i]); */
  /*   gkyl_yield_model_release(yield_model[i]); */
  /* } */
  /* gkyl_elastic_model_release(elastic_model); */
  gkyl_bc_emission_release(bc_ctx);

mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif

  return 0;
}
