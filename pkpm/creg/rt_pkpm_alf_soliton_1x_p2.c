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

struct pkpm_alf_ctx {
  double epsilon0;
  double mu0;
  double chargeElc; // electron charge
  double massElc; // electron mass
  double chargeIon; // ion charge
  double massIon; // ion mass
  double Te_Ti; // electron to ion temperature ratio
  double n0;
  double vAe;
  double B0;
  double beta_elc;
  double beta_ion;
  double vt_elc;
  double vt_ion;
  double nuElc;
  double nuIon;
  double di; 
  double a;
  double delta_B0;
  double Lx; // Domain size (x-direction).
  int Nx; // Cell count (x-direction).
  double cfl_frac; // CFL coefficient.
  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  double init_dt; // Initial time step guess so first step does not generate NaN
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
  bool use_gpu;
};

static inline double
maxwellian(double n, double v, double vth)
{
  double v2 = v*v;
  return n/sqrt(2*M_PI*vth*vth)*exp(-v2/(2*vth*vth));
}

static inline double
sech(double x)
{
  return 1.0/(cosh(x));
}

static inline double
sech2(double x)
{
  return 1.0/(cosh(x)*cosh(x));
}

static inline double
tanh2(double x)
{
  return tanh(x)*tanh(x);
}

void
evalDistFuncElc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_alf_ctx *app = ctx;
  
  double x = xn[0], vx = xn[1];

  double me = app->massElc;
  double B0 = app->B0;
  double n0 = app->n0;
  double beta_elc = app->beta_elc; 
  double B0perp = app->delta_B0;
  double a = app->a;
  double di = app->di; 
  double n = 2.0*a*sech(a*x/di) + n0;

  double arg = 0.5*a*x/di;
  double phi = 4.0*atan(tanh(arg));
  double B_x = B0;
  double B_y = B0perp*sin(phi);
  double B_z = B0perp*cos(phi);

  double magB2 = 0.5*(B_x*B_x + B_y*B_y + B_z*B_z);
  double Te = magB2*beta_elc/n;
  double vt_elc = sqrt(Te/me);
  
  double fv = maxwellian(n, vx, vt_elc);
    
  fout[0] = fv;
  fout[1] = vt_elc*vt_elc*fv;
}
void
evalDistFuncIon(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_alf_ctx *app = ctx;
  
  double x = xn[0], vx = xn[1];

  double mi = app->massIon;
  double B0 = app->B0;
  double n0 = app->n0;
  double beta_ion = app->beta_ion; 
  double B0perp = app->delta_B0;
  double a = app->a;
  double di = app->di; 
  double n = 2.0*a*sech(a*x/di) + n0;

  double arg = 0.5*a*x/di;
  double phi = 4.0*atan(tanh(arg));
  double B_x = B0;
  double B_y = B0perp*sin(phi);
  double B_z = B0perp*cos(phi);

  double magB2 = 0.5*(B_x*B_x + B_y*B_y + B_z*B_z);
  double Ti = magB2*beta_ion/n;
  double vt_ion = sqrt(Ti/mi);

  double fv = maxwellian(n, vx, vt_ion);
    
  fout[0] = fv;
  fout[1] = vt_ion*vt_ion*fv;
}

void
evalFluidElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_alf_ctx *app = ctx;
  
  double x = xn[0];

  double qe = app->chargeElc;
  double qi = app->chargeIon;
  double me = app->massElc;
  double mi = app->massIon;
  double B0perp = app->delta_B0;

  double a = app->a;
  double di = app->di; 
  double arg = 0.5*a*x/di;
  double phi = 4.0*atan(tanh(arg));
  double Jy = 2.0*B0perp*(a*sech2(arg)*sin(phi)/(di*tanh2(arg) + di));
  double Jz = 2.0*B0perp*(a*sech2(arg)*cos(phi)/(di*tanh2(arg) + di));

  double vdrift_x = 0.0;
  double vdrift_y = Jy/qe;
  double vdrift_z = Jz/qe;

  fout[0] = me*vdrift_x;
  fout[1] = me*vdrift_y;
  fout[2] = me*vdrift_z;
}

void
evalFluidIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_alf_ctx *app = ctx;
  
  double x = xn[0];

  double qe = app->chargeElc;
  double qi = app->chargeIon;
  double me = app->massElc;
  double mi = app->massIon;
  double B0perp = app->delta_B0;

  double vdrift_x = 0.0;
  double vdrift_y = 0.0;
  double vdrift_z = 0.0;

  fout[0] = mi*vdrift_x;
  fout[1] = mi*vdrift_y;
  fout[2] = mi*vdrift_z;
}

void
evalFieldFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_alf_ctx *app = ctx;

  double x = xn[0];

  double qe = app->chargeElc;
  double qi = app->chargeIon;
  double L = app->Lx;
  double B0 = app->B0;
  double B0perp = app->delta_B0;

  double a = app->a;
  double di = app->di; 
  double arg = 0.5*a*x/di;
  double phi = 4.0*atan(tanh(arg));

  double B_x = B0;
  double B_y = B0perp*sin(phi);
  double B_z = B0perp*cos(phi);

  // Assumes qi = abs(qe)
  double Jy = 2.0*B0perp*(a*sech2(arg)*sin(phi)/(di*tanh2(arg) + di));
  double Jz = 2.0*B0perp*(a*sech2(arg)*cos(phi)/(di*tanh2(arg) + di));
  double n = 2.0*a*sech(a*x/di);
  double u_xe = 0.0;
  double u_ye = Jy/(n*qe);
  double u_ze = Jz/(n*qe);

  // E = - v_e x B ~  (J - u) x B
  double E_x = - (u_ye*B_z - u_ze*B_y);
  double E_y = - (u_ze*B_x - u_xe*B_z);
  double E_z = - (u_xe*B_y - u_ye*B_x);
  
  fout[0] = E_x; fout[1] = E_y, fout[2] = E_z;
  fout[3] = B_x; fout[4] = B_y; fout[5] = B_z;
  fout[6] = 0.0; fout[7] = 0.0;
}

void
evalNuElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_alf_ctx *app = ctx;
  fout[0] = app->nuElc;
}

void
evalNuIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_alf_ctx *app = ctx;
  fout[0] = app->nuIon;
}

struct pkpm_alf_ctx
create_ctx(void)
{
  double epsilon0 = 1.0; // permittivity of free space
  double mu0 = 1.0; // pemiability of free space

  double massElc = 1.0; // electron mass
  double chargeElc = -1.0; // electron charge
  double massIon = 100.0; // ion mass
  double chargeIon = 1.0; // ion charge

  double Te_Ti = 1.0; // ratio of electron to ion temperature
  // initial conditions
  double a = 0.01;
  double n0 = 1.0; // initial number density 
  double vAe = 0.25;
  double beta_elc = 0.5;

  double B0 = vAe*sqrt(mu0*n0*massElc);
  double delta_B0 = a*B0;
  double vt_elc = vAe*sqrt(beta_elc/2.0);
  // ion velocities
  double vAi = vAe/sqrt(massIon);
  double vt_ion = vt_elc/sqrt(massIon*Te_Ti); //Ti/Te = 1.0
  double beta_ion = beta_elc/Te_Ti;

  // ion cyclotron frequency and gyroradius
  double omegaCi = chargeIon*B0/massIon;
  double di = vAi/omegaCi;
  double rhoi = sqrt(2.0)*vt_ion/omegaCi;

  // collision frequencies
  double nuElc = 0.0001*omegaCi;
  double nuIon = 0.0001*omegaCi/sqrt(massIon);

  double Lx = 10.0*(di/a);
  int Nx = 16; 
  double dx = Lx/Nx;
  double cfl_frac = 1.0; // CFL coefficient.
  double t_end = 10000.0/omegaCi; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.
  double init_dt = (Lx/Nx)/(5.0);

  struct pkpm_alf_ctx ctx = {
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .chargeElc = chargeElc,
    .massElc = massElc,
    .chargeIon = chargeIon,
    .massIon = massIon,
    .Te_Ti = Te_Ti,
    .n0 = n0,
    .vAe = vAe,
    .B0 = B0,
    .beta_elc = beta_elc,
    .beta_ion = beta_ion,
    .vt_elc = vt_elc,
    .vt_ion = vt_ion,
    .nuElc = nuElc,
    .nuIon = nuIon,
    .di = di, 
    .a = a, 
    .delta_B0 = delta_B0,
    .Lx = Lx,
    .Nx = Nx, 
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .num_frames = num_frames,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
    .init_dt = init_dt, 
  };
  return ctx;
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_pkpm_app* app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr)) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_pkpm_app_write(app, t_curr, iot->curr-1);
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

  struct pkpm_alf_ctx ctx = create_ctx(); // context for init functions

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int VX = APP_ARGS_CHOOSE(app_args.vcells[0], 32);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
  
  // electrons
  struct gkyl_pkpm_species elc = {
    .name = "elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -6.0 * ctx.vt_elc},
    .upper = { 6.0 * ctx.vt_elc}, 
    .cells = { VX },

    .ctx_dist = &ctx,
    .ctx_fluid = &ctx,
    .init_dist = evalDistFuncElc,
    .init_fluid = evalFluidElc,

    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,

      .ctx = &ctx,
      .self_nu = evalNuElc,
    },    
  };
  
  // ions
  struct gkyl_pkpm_species ion = {
    .name = "ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { -6.0 * ctx.vt_ion},
    .upper = { 6.0 * ctx.vt_ion}, 
    .cells = { VX },

    .ctx_dist = &ctx,
    .ctx_fluid = &ctx,
    .init_dist = evalDistFuncIon,
    .init_fluid = evalFluidIon,

    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,

      .ctx = &ctx,
      .self_nu = evalNuIon,
    },    
  };

  // field
  struct gkyl_pkpm_field field = {
    .epsilon0 = 1.0, .mu0 = 1.0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,

    .ctx = &ctx,
    .init = evalFieldFunc
  };

  int nrank = 1; // Number of processes in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
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

  int ccells[] = { NX };
  int cdim = sizeof(ccells) / sizeof(ccells[0]);
  int ncuts = 1;
  for (int d = 0; d < cdim; d++) {
    ncuts *= app_args.cuts[d];
  }

  if (ncuts != comm_size) {
    if (my_rank == 0) {
      fprintf(stderr, "*** Number of ranks, %d, does not match total cuts, %d!\n", comm_size, ncuts);
    }
    goto mpifinalize;
  }

  // pkpm app
  struct gkyl_pkpm app_inp = {

    .cdim = 1, .vdim = 1,
    .lower = { -ctx.Lx },
    .upper = { ctx.Lx },
    .cells = { NX},
    .poly_order = 2,
    .basis_type = app_args.basis_type,
    .cfl_frac = ctx.cfl_frac,

    // .use_explicit_source = true, 
    
    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },

    .num_species = 2,
    .species = { elc, ion },
    .field = field,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0] },
      .comm = comm,
    },
  };

  // create app object
  // Set app output name from the executable name (argv[0]).
  snprintf(app_inp.name, sizeof(app_inp.name), "%s", app_args.app_name);
  gkyl_pkpm_app *app = gkyl_pkpm_app_new(&app_inp);

  // start, end and initial time-step
  double t_curr = 0.0, t_end = ctx.t_end;
  double dt = ctx.init_dt;

  // Create trigger for IO.
  int num_frames = ctx.num_frames;
  struct gkyl_tm_trigger io_trig = { .dt = t_end / num_frames };

  // initialize simulation
  gkyl_pkpm_app_apply_ic(app, t_curr);
  write_data(&io_trig, app, t_curr, false);  

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

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    write_data(&io_trig, app, t_curr, false);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_pkpm_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_pkpm_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_pkpm_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_pkpm_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_pkpm_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);
        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }
  write_data(&io_trig, app, t_curr, false);
  gkyl_pkpm_app_stat_write(app);

  // fetch simulation statistics
  struct gkyl_pkpm_stat stat = gkyl_pkpm_app_stat(app);

  gkyl_pkpm_app_cout(app, stdout, "\n");
  gkyl_pkpm_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_pkpm_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_pkpm_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_pkpm_app_cout(app, stdout, "Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_pkpm_app_cout(app, stdout, "Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  gkyl_pkpm_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_pkpm_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_pkpm_app_cout(app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm);
  gkyl_pkpm_app_cout(app, stdout, "Fluid Species RHS calc took %g secs\n", stat.fluid_species_rhs_tm);
  gkyl_pkpm_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_pkpm_app_cout(app, stdout, "Species PKPM Vars took %g secs\n", stat.species_pkpm_vars_tm);
  gkyl_pkpm_app_cout(app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm);
  gkyl_pkpm_app_cout(app, stdout, "EM Variables (bvar) calculation took %g secs\n", stat.field_em_vars_tm);
  gkyl_pkpm_app_cout(app, stdout, "Current evaluation and accumulate took %g secs\n", stat.current_tm);

  gkyl_pkpm_app_cout(app, stdout, "Species BCs took %g secs\n", stat.species_bc_tm);
  gkyl_pkpm_app_cout(app, stdout, "Fluid Species BCs took %g secs\n", stat.fluid_species_bc_tm);
  gkyl_pkpm_app_cout(app, stdout, "Field BCs took %g secs\n", stat.field_bc_tm);
  
  gkyl_pkpm_app_cout(app, stdout, "Updates took %g secs\n", stat.total_tm);
  
  gkyl_pkpm_app_cout(app, stdout, "Number of write calls %ld,\n", stat.n_io);
  gkyl_pkpm_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

  gkyl_comm_release(comm);

  // simulation complete, free app
  gkyl_pkpm_app_release(app);

  mpifinalize:
  ;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif  
  
  return 0;
}
