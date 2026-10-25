#include <math.h>

#include <gkyl_const.h>
#include <gkyl_flash_block_geom.h>
#include <gkyl_flash_multib.h>
#include <gkyl_null_comm.h>
#include <gkyl_util.h>

#ifdef GKYL_HAVE_MPI
#include <gkyl_mpi_comm.h>
#include <mpi.h>
#endif

#include <rt_arg_parse.h>

struct flash_sod_ctx {
  double mass; // Species mass.
  double rhol, uxl, pl; // Left side density, velocity, and pressure.
  double rhor, uxr, pr; // Right side density, velocity, and pressure.
  double vthl; // Left side thermal velocity.
  double Bx0; // Magnetic field in kinetic direction, ALWAYS set to 1.0.
  double gas_gamma; // Adiabatic index.
  double nu; // Collision frequency
  double Lx; // Domain length in x.
  double xmin; // Starting value of grid in x-direction.
  int num_blocks; // Number of blocks of equal size.
  int fluid_cells_per_block; // Number of fluid cells per block.
  int NX_fluid; // Total number of fluid cells (summed over blocks).
  int NV; // Number of velocity space cells.
  int num_frames; // Number of output frames.
  double t_end; // Simulation end time.
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
  double nu = 1000.0;

  double Lx = 1.0;
  double xmin = 0.0;

  int num_blocks = 3;
  int fluid_cells_per_block = 32;
  int NX_fluid = num_blocks*fluid_cells_per_block;
  int NV = 16;

  int num_frames = 10;
  double t_end = 0.01;

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
    .Lx = Lx,
    .xmin = xmin,
    .num_blocks = num_blocks,
    .fluid_cells_per_block = fluid_cells_per_block,
    .NX_fluid = NX_fluid,
    .NV = NV,
    .num_frames = num_frames,
    .t_end = t_end,
  };

  return ctx;
}

void init_dist(double t, const double *xn, double *fout, void *ctx) {
  struct flash_sod_ctx *app = (struct flash_sod_ctx *)ctx;
  double x = xn[0], vx = xn[1];

  double n, vth;
  if (fabs(x) <= 0.5) {
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
  if (fabs(x) <= 0.5) {
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
  fout[2] = rho*uy;
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

  // Electric field.
  fout[0] = 0.0; fout[1] = 0.0; fout[2] = 0.0;
  // Magnetic field.
  fout[3] = Bx0; fout[4] = 0.0; fout[5] = 0.0;
  // Correction potentials.
  fout[6] = 0.0; fout[7] = 0.0;
}

// Split the domain into num_blocks blocks of equal size.
//
//   +-----------+-----------+-----------+
//   |b0         |b1         |b2         |
//   +-----------+-----------+-----------+
//  xmin                              xmin+Lx
//
// The lower edge of b0 and the upper edge of the final block are physical boundaries.
struct gkyl_flash_block_geom*
create_block_geom(const struct flash_sod_ctx *ctx, int NX)
{
  int nb = ctx->num_blocks;
  struct gkyl_flash_block_geom *bgeom = gkyl_flash_block_geom_new(1, nb);

  // Divide the domain into blocks of equal size.
  double dx_block = ctx->Lx/nb;
  for (int b=0; b<nb; ++b) {
    struct gkyl_target_edge lo_conn = b == 0 ?
      (struct gkyl_target_edge) { .bid = b, .dir = 0, .edge = GKYL_PHYSICAL } :
      (struct gkyl_target_edge) { .bid = b-1, .dir = 0, .edge = GKYL_UPPER_POSITIVE };
    struct gkyl_target_edge up_conn = b == nb-1 ?
      (struct gkyl_target_edge) { .bid = b, .dir = 0, .edge = GKYL_PHYSICAL } :
      (struct gkyl_target_edge) { .bid = b+1, .dir = 0, .edge = GKYL_LOWER_POSITIVE };

    gkyl_flash_block_geom_set_block(bgeom, b, &(struct gkyl_flash_block_geom_info) {
        .block_id = b,
        .lower = { ctx->xmin + b*dx_block },
        .upper = { ctx->xmin + (b+1)*dx_block },
        .cells = { ctx->fluid_cells_per_block },
        .cuts = { 1 },
        .connections[0] = { lo_conn, up_conn },
      }
    );
  }

  return bgeom;
}

int main(int argc, char **argv) {
  struct gkyl_app_args args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
  if (args.use_mpi) {
    MPI_Init(&argc, &argv);
  }
#endif

  struct flash_sod_ctx ctx = create_ctx();

  int NX = APP_ARGS_CHOOSE(args.xcells[0], ctx.NX_fluid);
  int NV = APP_ARGS_CHOOSE(args.vcells[0], ctx.NV);

  // Each block needs an even number of fluid cells (kinetic grid is 2x coarser).
  if (NX % (2*ctx.num_blocks) != 0) {
    fprintf(stderr, "Number of fluid cells (%d) must be divisible by %d.\n", NX, 2*ctx.num_blocks);
    return 1;
  }

  // Construct communicator for use in app.
  struct gkyl_comm *comm = 0;
#ifdef GKYL_HAVE_MPI
  if (args.use_mpi)
    comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) { .mpi_comm = MPI_COMM_WORLD } );
  else
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) { .use_gpu = false } );
#else
  comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) { .use_gpu = false } );
#endif

  // Construct block geometry.
  struct gkyl_flash_block_geom *bgeom = create_block_geom(&ctx, NX);
  int nblocks = gkyl_flash_block_geom_num_blocks(bgeom);

  // Initial conditions are the same function of global x in every block.
  struct gkyl_flash_multib_species_inp_per_block spec_blocks[1] = {
    {
      .block_id = 0,
      .ctx = &ctx,
      .init_dist = init_dist,
      .init_fluid = init_fluid,
    }
  };

  // Set domain boundary conditions for exterior blocks.
  struct gkyl_flash_block_physical_bcs spec_bcs[] = {
    { .block_id = 0, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_COPY },
    { .block_id = nblocks-1, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_COPY },
  };

  struct gkyl_flash_multib_species spec = {
    .name_kin = "kinetic",
    .name_flu = "fluid",
    .mass = ctx.mass,
    .cells = { NV },
    .lower = { -5.0*ctx.vthl },
    .upper = { 5.0*ctx.vthl },

    .gas_gamma = ctx.gas_gamma,
    .limiter = GKYL_MIN_MOD,

    // Domain boundary conditions
    .num_physical_bcs = 2,
    .bcs = spec_bcs,

    .write_pkpm_diags = true,

    // Duplicate initial condition function across all blocks?
    .duplicate_ics = true,
    .blocks = spec_blocks,

    .collisions = {
      .num_cross_collisions = 0,
      .ctx = &ctx,
      .self_nu = init_nu,
    },
  };

  // EM field initial conditions. 
  // Only specifying one because it will be duplicated across blocks.
  struct gkyl_flash_multib_field_inp_per_block field_blocks[1] = {
    {
      .block_id = 0,
      .ctx = &ctx,
      .init = init_field,
    }
  };

  // Initialize multiblock field.
  struct gkyl_flash_multib_field field = {
    .mu0 = 1.0, .epsilon0 = 1.0,
    .is_static = true,
    .duplicate_ics = true,
    .blocks = field_blocks,
  };

  // Initialize multiblock app input.
  struct gkyl_flash_multib_inp inp = {
    .name = "multib_sodshock",
    .cdim = 1, .vdim = 1,
    .num_blocks_global = nblocks,
    .block_geom = bgeom,

    .cfl_frac = 1.0,

    .num_periodic_dir = 0,
    .periodic_dirs = { 0 },

    .num_species = 1,
    .species = { spec },

    .field = field,

    .comm = comm,
  };

  gkyl_flash_multib_app *app = gkyl_flash_multib_app_new(&inp);
  gkyl_flash_multib_app_apply_ic(app, 0.0);

  double tcurr = 0.0;
  double tend = ctx.t_end;
  double frame_dt = tend/ctx.num_frames;
  int frame = 0;

  gkyl_flash_multib_app_write(app, tcurr, frame++);

  double dt = tend;
  long step = 1;
  struct timespec tm = gkyl_wall_clock();

  while ((tcurr < tend) && (step <= args.num_steps)) {
    printf("Taking time step %ld at t = %g ... ", step, tcurr);
    struct gkyl_update_status status = gkyl_flash_multib_update(app, fmin(dt, tend - tcurr));
    printf("dt = %g\n", status.dt_actual);

    tcurr += status.dt_actual;
    dt = status.dt_suggested;
    step += 1;

    if (tcurr >= frame*frame_dt || tcurr >= tend) {
      gkyl_flash_multib_app_write(app, tcurr, frame++);
    }
  }

  double runtime = gkyl_time_diff_now_sec(tm);

  printf("\n");
  printf("Took %ld steps to final time %g.\n", step-1, tcurr);
  printf("Total time: %5.5f\n", runtime);

  gkyl_flash_multib_app_release(app);
  gkyl_flash_block_geom_release(bgeom);
  gkyl_comm_release(comm);

#ifdef GKYL_HAVE_MPI
  if (args.use_mpi) {
    MPI_Finalize();
  }
#endif

  return 0;
}
