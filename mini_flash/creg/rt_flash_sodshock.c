#include <math.h>

#include <gkyl_const.h>
#include <gkyl_flash.h>
#include <gkyl_util.h>

#include <rt_arg_parse.h>

struct flash_sod_ctx {
  double mass;
  double rhol, uxl, pl;
  double rhor, uxr, pr;
  double vthl;
  double Bx0;
  double gas_gamma;
  double nu;
  int NX_fluid;
  int NV;
  int num_frames;
  double t_end;
};

struct flash_sod_ctx create_ctx() {
  double mass = 1.0;

  // Sod shock tube
  double nl = 1.0;
  double uxl = 0.0;
  double pl = 1.0;

  double nr = 0.125;
  double uxr = 0.0;
  double pr = 0.1;

  double rhol = mass*nl;
  double rhor = mass*nr;
  double vthl = sqrt(pl/rhol);

  double gas_gamma = 5.0/3.0;
  double nu = 100.0;

  struct flash_sod_ctx ctx = {
    .mass = mass,
    .rhol = rhol, 
    .uxl = uxl, 
    .pl = pl,
    .rhor = rhor, 
    .uxr = uxr, 
    .pr = pr,
    .vthl = vthl,
    .Bx0 = 1.0,
    .gas_gamma = gas_gamma,
    .nu = nu,
    .NX_fluid = 64,
    .NV = 16,
    .num_frames = 10,
    .t_end = 0.3,
  };

  return ctx;
}

void init_dist(double t, const double *xn, double *fout, void *ctx) {
  struct flash_sod_ctx *app = (struct flash_sod_ctx *)ctx;
  double x = xn[0], vx = xn[1];

  double n, vth;
  if (x <= 0.0) {
    n = app->rhol/app->mass;
    vth = sqrt(app->pl/app->rhol);
  }
  else {
    n = app->rhor/app->mass;
    vth = sqrt(app->pr/app->rhor);
  }

  double F0 = n/sqrt(2.0*GKYL_PI*vth*vth) * exp(-(vx*vx)/(2.0*vth*vth));

  double G = vth * vth * F0;
  fout[0] = F0;
  fout[1] = G;
}

void init_fluid(double t, const double *xn, double *fout, void *ctx) {
  struct flash_sod_ctx *app = (struct flash_sod_ctx *)ctx;
  double gas_gamma = app->gas_gamma;
  double x = xn[0];

  double rho = 0.0;
  double ux = 0.0;
  double press = 0.0;
  if (x <= 0.0) {
    rho = app->rhol;
    ux = app->uxl;
    press = app->pl;
  }
  else {
    rho = app->rhor;
    ux = app->uxr;
    press = app->pr;
  }

  double uy = 0.0;
  double uz = 0.0;
  double q_sq = ux*ux + uy*uy + uz*uz;
  double ener = press/(gas_gamma-1.0) + 0.5*rho*q_sq;

  // Set pressure, momentum density, energy density.
  fout[0] = rho;
  fout[1] = rho*ux;
  fout[2] = 0.0;
  fout[3] = rho*uz;
  fout[4] = ener;
}

void init_nu(double t, const double *xn, double *fout, void *ctx) {
  struct flash_sod_ctx *app = ctx;
  fout[0] = app->nu;
}

void init_field(double t, const double *xn, double *fout, void *ctx) {
  struct flash_sod_ctx *app = ctx;

  double Bx0 = app->Bx0;

  double Ex = 0.0; // Total electric field (x-direction).
  double Ey = 0.0; // Total electric field (y-direction).
  double Ez = 0.0; // Total electric field (z-direction).

  double Bx = Bx0; // Total magnetic field (x-direction).
  double By = 0.0; // Total magnetic field (y-direction).
  double Bz = 0.0; // Total magnetic field (z-direction).

  // Electric field.
  fout[0] = Ex; fout[1] = Ey; fout[2] = Ez;
  // Magnetic field.
  fout[3] = Bx; fout[4] = By; fout[5] = Bz;
  // Correction potentials.
  fout[6] = 0.0; fout[7] = 0.0;
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_flash_app* app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_write) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
      iot->curr += 1;
    }

  gkyl_flash_app_write(app, t_curr, frame);
  }
}

int main(int argc, char **argv) {
  struct gkyl_app_args args = parse_app_args(argc, argv);
  struct flash_sod_ctx ctx = create_ctx();

  int NX = APP_ARGS_CHOOSE(args.xcells[0], ctx.NX_fluid);
  int NV = APP_ARGS_CHOOSE(args.vcells[0], ctx.NV);

  // char fmt_kinetic[128] = "kinetic-%d";
  // char fmt_fluid[128] = "fluid-%d";
  
  struct gkyl_flash_species spec = {
    .mass = ctx.mass,
    .name_kinetic = "kinetic",
    .name_fluid = "fluid",
    .cells = { NV },
    .lower = { -5.0*ctx.vthl },
    .upper = { 5.0*ctx.vthl },


    .collisions = {
      .num_cross_collisions = 0,
      .ctx = &ctx,
      .self_nu = init_nu,
    },

    .init_dist = init_dist,
    .init_fluid = init_fluid,
    .ctx = &ctx,

    .gas_gamma = ctx.gas_gamma,
    .limiter = GKYL_MIN_MOD,

    .bcx = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },

    .write_pkpm_diags = true,
  };
  // snprintf(spec.name_kinetic, sizeof spec.name_kinetic, fmt_kinetic, (int)NX);
  // snprintf(spec.name_fluid, sizeof spec.name_fluid, fmt_fluid, (int)NX);

  struct gkyl_flash_field field = {
    .mu0 = 1.0, .epsilon0 = 1.0,
    .is_static = true,
    .ctx = &ctx,
    .init = init_field
  };

  struct gkyl_flash_inp inp = {
    .name = "flash_sodshock",
    .cdim = 1, .vdim = 1,
    .lower_conf = { -1.0 },
    .upper_conf = { 1.0 },

    .cells_fluid = { NX },

    .cfl_frac = 1.0,

    .num_species = 1,
    .species = { spec },

    .field = field,

    .num_periodic_dir = 0,
    .periodic_dirs = { 0 }
  };

  gkyl_flash_app *app = gkyl_flash_app_new(&inp);
  gkyl_flash_app_apply_ic(app, 0.0);

  double dt = 1.0e-6;
  double tcurr = 0.0;
  struct gkyl_update_status status;

  struct gkyl_tm_trigger io_trig = { .curr = 0, .dt = ctx.t_end/ctx.num_frames, .tcurr = tcurr };
  write_data(&io_trig, app, tcurr, false);

  double tend = ctx.t_end;
  long step = 1;
  struct timespec tm = gkyl_wall_clock();

  while ((tcurr < tend) && (step <= args.num_steps)) {
    gkyl_flash_app_cout(app, stdout, "Taking time step %ld at t = %g ... ", step, tcurr);
    status = gkyl_flash_update(app, dt);
    gkyl_flash_app_cout(app, stdout, "dt = %g\n", status.dt_actual);

    tcurr += status.dt_actual;
    dt = status.dt_suggested;

    step += 1;

    if (args.num_steps < INT_MAX) {
      int frame_step = (ctx.num_frames>args.num_steps) ? 1 : (int)args.num_steps/ctx.num_frames;
      if ((step%frame_step == 0) && (step != 0)) {
        gkyl_flash_app_write(app, tcurr, io_trig.curr);
        io_trig.curr += 1;
      }
    }
    else {
      write_data(&io_trig, app, tcurr, false);
    }
  }

  double runtime = gkyl_time_diff_now_sec(tm);

  struct gkyl_flash_stat stat = gkyl_flash_app_stat(app);

  gkyl_flash_app_cout(app, stdout, "\n");
  gkyl_flash_app_cout(app, stdout, "Took %ld steps to final time %g.\n\n", step-1, tcurr);
  gkyl_flash_app_cout(app, stdout, "Total time spent initializing app: %5.5f\n", stat.init_tm);
  gkyl_flash_app_cout(app, stdout, "Time spent initializing all species: %5.5f\n", stat.init_species_tm);
  gkyl_flash_app_cout(app, stdout, "Number of forward Euler calls: %d\n", stat.nfeuler);
  gkyl_flash_app_cout(app, stdout, "Number of failed update calls: %d\n", stat.nfail);
  gkyl_flash_app_cout(app, stdout, "Species updates took %5.5f sec\n", stat.species_rhs_tm);
  gkyl_flash_app_cout(app, stdout, "Fluid updates took %5.5f sec\n", stat.fluid_rhs_tm);
  gkyl_flash_app_cout(app, stdout, "I/O took %5.5f sec\n", stat.io_tm);
  gkyl_flash_app_cout(app, stdout, "Total time: %5.5f\n", runtime);

  gkyl_flash_app_release(app);

  return 0;
}
