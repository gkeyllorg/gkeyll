#include <gkyl_alloc.h>
#include <gkyl_moment_multib.h>
#include <gkyl_null_comm.h>
#include <gkyl_wv_euler.h>

#include <rt_arg_parse.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

// Gas constant
static const double gas_gamma = 1.4;

struct gkyl_block_geom*
create_block_geom(void)
{
  struct gkyl_block_geom *bgeom = gkyl_block_geom_new(2, 3);

  /* Block layout and coordinates

   Y  
   ^  
   |  
   2  +------+
   |  |b0    |
   |  |      |
   1  +------+-----+
   |  |b1    |b2   |
   |  |      |     |
   0  +------+-----+

      0 -----1-----2 -> X
  */  

  // block 0
  gkyl_block_geom_set_block(bgeom, 0, &(struct gkyl_block_geom_info) {
      .lower = { 0, 1 },
      .upper = { 1, 2 },
      .cells = { 128, 128 },
      .cuts = { 1, 1 },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }, // physical boundary
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }  // physical boundary
      },
      .connections[1] = { // y-direction connections
        { .bid = 1, .dir = 1, .edge = GKYL_UPPER_POSITIVE },
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL } // physical boundary
      }
    }
  );
  
  // block 1
  gkyl_block_geom_set_block(bgeom, 1, &(struct gkyl_block_geom_info) {
      .lower = { 0, 0 },
      .upper = { 1, 1 },
      .cells = { 128, 128 },
      .cuts = { 1, 1 },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }, // physical boundary
        { .bid = 2, .dir = 0, .edge = GKYL_LOWER_POSITIVE }
      },
      .connections[1] = { // y-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL }, // physical boundary
        { .bid = 0, .dir = 1, .edge = GKYL_LOWER_POSITIVE }
      }
    }
  );

  // block 2
  gkyl_block_geom_set_block(bgeom, 2, &(struct gkyl_block_geom_info) {
      .lower = { 1, 0 },
      .upper = { 2, 1 },
      .cells = { 128, 128 },
      .cuts = { 1, 1 },
      
      .connections[0] = { // x-direction connections
        { .bid = 1, .dir = 0, .edge = GKYL_UPPER_POSITIVE },
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL } // physical boundary
      },
      .connections[1] = { // y-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL }, // physical boundary
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL } // physical boundary
      }
    }
  );

  return bgeom;
}

static void
write_data(struct gkyl_tm_trigger* iot, gkyl_moment_multib_app* app,
  double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr)) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }
    gkyl_moment_multib_app_write(app, t_curr, frame);
  }
}

void
initFluidSod(double t, const double *xn, double* restrict fout, void *ctx)
{
  double xsloc = 1.25, ysloc = 1.5;
  double x = xn[0], y = xn[1];

  double rho = 0.125, pr = 0.1;
  if (y>ysloc || x>xsloc) {
    rho = 1.0;
    pr = 1.0;
  }
  
  fout[0] = rho;
  fout[1] = 0.0; fout[2] = 0.0; fout[3] = 0.0;
  fout[4] = pr/(gas_gamma-1);
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  struct gkyl_comm *comm = 0;
  if (app_args.use_mpi) {
#ifdef GKYL_HAVE_MPI
    MPI_Init(&argc, &argv);
    comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
        .sync_corners = true
      }
    );
#endif
  }
  if (comm == 0)
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) { } );

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
  
  // construct block geometry
  struct gkyl_block_geom *bgeom = create_block_geom();
  int nblocks = gkyl_block_geom_num_blocks(bgeom);

  struct gkyl_wv_eqn *euler_eqn = gkyl_wv_euler_inew( &(struct gkyl_wv_euler_inp) {
      .gas_gamma = gas_gamma
    }
  );

  // all data is common across blocks
  struct gkyl_moment_multib_species_pb euler_blocks[1];
  euler_blocks[0] = (struct gkyl_moment_multib_species_pb) {
    .init = initFluidSod,
  };

  struct gkyl_block_physical_bcs euler_phys_bcs[] = {
    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_REFLECT },
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_REFLECT },
    { .bidx = 0, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_COPY },
    // block 1 BCs
    { .bidx = 1, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_REFLECT },
    { .bidx = 1, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_REFLECT },
    // block 2 BCs
    { .bidx = 2, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_COPY },
    { .bidx = 2, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_REFLECT },
    { .bidx = 2, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_REFLECT },
  };

  struct gkyl_moment_multib_species euler = {
    .name = "euler",
    .charge = 0.0,
    .mass = 1.0,
    .equation = euler_eqn,

    .duplicate_across_blocks = true,
    .blocks = euler_blocks,

    .num_physical_bcs = 8,
    .bcs = euler_phys_bcs,
  };

  struct gkyl_moment_multib app_inp = {

    .block_geom = bgeom,
    .cfl_frac = 0.9,

    .num_species = 1,
    .species = { euler },

    .comm = comm
  };

  // Set app output name from the executable name (argv[0]).
  snprintf(app_inp.name, sizeof(app_inp.name), "%s", app_args.app_name);
  struct gkyl_moment_multib_app *app = gkyl_moment_multib_app_new(&app_inp);
  if (app == 0) {
    fprintf(stderr, "Unable to construct multiblock App! (Probably too few or too many ranks)\n");
    goto finish;
  }
  
  gkyl_moment_multib_app_write_topo(app);

  // Initial and final simulation times.
  double t_curr = 0.0, t_end = 0.6;

  // Create trigger for IO.
  int num_frames = 4;
  struct gkyl_tm_trigger io_trig = { .dt = t_end / num_frames };

  // Initialize simulation.
  gkyl_moment_multib_app_apply_ic(app, t_curr);
  write_data(&io_trig, app, t_curr, false);

  gkyl_comm_release(comm);
  gkyl_block_geom_release(bgeom);
  gkyl_wv_eqn_release(euler_eqn);
  gkyl_moment_multib_app_release(app);

  finish:

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif
  
  return 0;
}
