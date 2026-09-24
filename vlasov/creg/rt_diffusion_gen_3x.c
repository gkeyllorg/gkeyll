#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>
#include <gkyl_util.h>
#include <gkyl_wv_advect.h>
#include <rt_arg_parse.h>

struct sim_ctx {
  double a; // advection velocity
  double Lx; // size of the box
};

void
evalInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct sim_ctx *app = ctx;
  double x = xn[0], y = xn[1], z = xn[2];
  if (x > -1 & x < 1 & y > -1 & y < 1 & z > -1 & z < 1)
    fout[0] = 1.0;
  else
    fout[0] = 0.0;
}

void
D(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct sim_ctx *app = ctx;
  double x = xn[0], y = xn[1], z = xn[2];
  fout[0] = 1.0;//x + 2.0;
  fout[1] = 1.0;
  fout[2] = 1.0;
  fout[3] = 1.0;//y + 2.0;
  fout[4] = 1.0;
  fout[5] = 1.0;//z + 2.0;
}

void
eval_advect_vel(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct sim_ctx *app = ctx;
  double x = xn[0], y = xn[1], z = xn[2];
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
}

struct sim_ctx
create_ctx(void)
{
  struct sim_ctx ctx = {
    .Lx = 2,
  };
  return ctx;
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
  struct sim_ctx ctx = create_ctx(); // context for init functions
  
  // Equation object for getting equation type, 
  // advection velocity is set by eval_advect_vel function. 
  double c = 1.0;
  struct gkyl_wv_eqn *advect = gkyl_wv_advect_new(c, false);
  
  struct gkyl_vlasov_fluid_species f = {


    .ctx = &ctx,
    .init = evalInit,
    .equation = advect,
    .advection = {
      .velocity = eval_advect_vel,
      .velocity_ctx = &ctx,
    },
    .diffusion = {.Dij = D, .Dij_ctx = 0},
  };  

  // VM app
  struct gkyl_vm vm = {
    .name = "rt-diffusion-gen-3x",

    .cdim = 3, .vdim = 0,
    .lower = { -2, -2, -2 },
    .upper = { 2, 2, 2 },
    .cells = { 16, 16, 16 },
    .poly_order = 2,
    .basis_type = app_args.basis_type,

    .num_periodic_dir = 3,
    .periodic_dirs = { 0, 1, 2 },

    .num_species = 1,
    .species = {
      { .name = "f", .charge = 0.0, .mass = 1.0,
        .type = GKYL_SPECIES_FLUID, .fluid = f },
    },
    .cfl_frac = 0.9,

    .skip_field = true,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
    },
  };

  // create app object
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&vm);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 0.01;
  double dt = tend-tcurr;

  // initialize simulation
  gkyl_vlasov_app_apply_ic(app, tcurr);
  
  gkyl_vlasov_app_write(app, tcurr, 0);
  gkyl_vlasov_app_write_mom(app, tcurr, 0);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    printf("Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_vlasov_update(app, dt);
    printf(" dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      printf("** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;
    step += 1;
  }

  gkyl_vlasov_app_write(app, tcurr, 1);
  gkyl_vlasov_app_write_mom(app, tcurr, 1);
  gkyl_vlasov_app_stat_write(app);

  // fetch simulation statistics
  struct gkyl_vlasov_stat stat = gkyl_vlasov_app_stat(app);

  // simulation complete, free app
  gkyl_wv_eqn_release(advect);
  gkyl_vlasov_app_release(app);

  printf("\n");
  printf("Number of update calls %ld\n", stat.nup);
  printf("Number of forward-Euler calls %ld\n", stat.nfeuler);
  printf("Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    printf("Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    printf("Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  printf("Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  printf("Species RHS calc took %g secs\n", stat.species_rhs_tm);
  //printf("Field RHS calc took %g secs\n", stat.field_rhs_tm);
  //printf("Current evaluation and accumulate took %g secs\n", stat.current_tm);
  printf("Updates took %g secs\n", stat.total_tm);
  
  return 0;
}
