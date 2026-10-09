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
evalSREulerInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct sr_euler_ctx *app = ctx;
  double gas_gamma = app->gas_gamma;

  double x = xn[0];
  //ICs from test 4, shock tube, in Eulderink 1995

  double rhol = 10.0, ul = 0.0, pl = 40. / 3.;
  double rhor = 1.0, ur = 0.0, pr = 2. / (3.e7);

  double rho = rhor, u = ur, p = pr;
  if (x < 45.) {
    rho = rhol;
    u = ul;
    p = pl;
  }

  double gamma = 1 / sqrt(1 - u * u);
  double rhoh = gas_gamma * p / (gas_gamma - 1) + rho;

  fout[0] = gamma * rho;
  fout[1] = gamma * gamma * rhoh - p;
  fout[2] = gamma * gamma * rhoh * u;
  fout[3] = 0.;
  fout[4] = 0.;
}

struct sr_euler_ctx
sr_euler_ctx(void)
{
  return (struct sr_euler_ctx){.gas_gamma = 5. / 3.};
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 100);

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

    .bcx = {GKYL_SPECIES_COPY, GKYL_SPECIES_COPY},
  };

  // VM app
  struct gkyl_moment app_inp = {

    .ndim = 1,
    .lower = {0.0},
    .upper = {100.0},
    .cells = {NX},

    .cfl_frac = 0.9,

    .num_species = 1,
    .species = {fluid},
  };

  // create app object
  // Set app output name from the executable name (argv[0]).
  snprintf(app_inp.name, sizeof(app_inp.name), "%s", app_args.app_name);
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 50.;

  // initialize simulation
  gkyl_moment_app_apply_ic(app, tcurr);
  gkyl_moment_app_write(app, tcurr, 0);

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

    step += 1;
  }

  gkyl_moment_app_write(app, tcurr, 1);
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
