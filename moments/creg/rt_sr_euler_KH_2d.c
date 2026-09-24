#include <math.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_sr_euler.h>
#include <rt_arg_parse.h>

struct sr_euler_ctx {
  double gas_gamma; // gas constant
};

void
evalSREulerInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct sr_euler_ctx *app = ctx;
  double gas_gamma = app->gas_gamma;
  
  double x = xn[0], y = xn[1];
  //ICs from KH test from sec 4.3.2 Stone Athena++ 2020
  
  double rho, rhou = 1.5, rhol = 0.5, p = 20.;
  double pi = 3.141592653589793238462643383279502884;
  
  double u = 0.25*tanh(100*y);
  double v = sin(2*pi*x)*exp(-100*y*y)/400;
  
  rho = rhol;
  if (y > 0.) {
    rho = rhou;
  }
  
  double gamma = 1 / sqrt(1 - u*u - v*v);
  double rhoh = gas_gamma * p / (gas_gamma - 1)  + rho;
  
  fout[0] = gamma*rho;
  fout[1] = gamma*gamma*rhoh - p;
  fout[2] = gamma*gamma*rhoh*u;
  fout[3] = gamma*gamma*rhoh*v;
  fout[4] = 0.;
}

struct sr_euler_ctx
sr_euler_ctx(void)
{
  return (struct sr_euler_ctx) { .gas_gamma = 4./3. };
}

void
write_data(struct gkyl_tm_trigger *iot, const gkyl_moment_app *app, double tcurr)
{
  if (gkyl_tm_trigger_check_and_bump(iot, tcurr))
    gkyl_moment_app_write(app, tcurr, iot->curr-1);
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 128);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], 64);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
  struct sr_euler_ctx ctx = sr_euler_ctx(); // context for init functions

  // equation object
  struct gkyl_wv_eqn *sr_euler = gkyl_wv_sr_euler_new(ctx.gas_gamma);

  struct gkyl_moment_species fluid = {
    .name = "sr_euler",

    .equation = sr_euler,

    .ctx = &ctx,
    .init = evalSREulerInit,
  };

  // VM app
  struct gkyl_moment app_inp = {

    .ndim = 2,
    .lower = { 0.0, -0.25 },
    .upper = { 1.0, 0.25 }, 
    .cells = { NX, NY },

    .cfl_frac = 0.9,

    .num_species = 1,
    .species = { fluid },
    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },
  };

  // create app object
  // Set app output name from the executable name (argv[0]).
  snprintf(app_inp.name, sizeof(app_inp.name), "%s", app_args.app_name);
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 5.;
  int nframe = 1;
  
  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = tend/nframe };
  
  // initialize simulation
  gkyl_moment_app_apply_ic(app, tcurr);
  write_data(&io_trig, app, tcurr);
  //gkyl_moment_app_write(app, tcurr, 0);

  // compute estimate of maximum stable time-step
  double dt = gkyl_moment_app_max_dt(app);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    printf("Taking time-step %ld at t = %g ...", step, tcurr);
    struct gkyl_update_status status = gkyl_moment_update(app, dt);
    printf(" dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      printf("** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;

    write_data(&io_trig, app, tcurr);
    
    step += 1;
  }

  //gkyl_moment_app_write(app, tcurr, 1);
  gkyl_moment_app_stat_write(app);

  struct gkyl_moment_stat stat = gkyl_moment_app_stat(app);

  // simulation complete, free resources
  gkyl_wv_eqn_release(sr_euler);
  gkyl_moment_app_release(app);

  printf("\n");
  printf("Number of update calls %ld\n", stat.nup);
  printf("Number of failed time-steps %ld\n", stat.nfail);
  printf("Species updates took %g secs\n", stat.species_tm);
  printf("Field updates took %g secs\n", stat.field_tm);
  printf("Total updates took %g secs\n", stat.total_tm);
  
  return 0;
}
