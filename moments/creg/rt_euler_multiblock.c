#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_block_topo.h>
#include <gkyl_fv_proj.h>
#include <gkyl_moment.h>
#include <gkyl_null_pool.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_thread_pool.h>
#include <gkyl_util.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wave_prop.h>
#include <gkyl_wv_apply_bc.h>
#include <gkyl_wv_euler.h>
#include <gkyl_wv_maxwell.h>
#include <rt_arg_parse.h>

#include <thpool.h>

// Gas constant
static const double gas_gamma = 1.4;

// ranges for use in BCs
struct skin_ghost_ranges {
  struct gkyl_range lower_skin[2];
  struct gkyl_range lower_ghost[2];

  struct gkyl_range upper_skin[2];
  struct gkyl_range upper_ghost[2];
};

void
skin_ghost_ranges_init(struct skin_ghost_ranges *sgr,
  const struct gkyl_range *parent, const int *ghost)
{
  int ndim = parent->ndim;
  
  for (int d=0; d<ndim; ++d) {
    gkyl_skin_ghost_ranges(&sgr->lower_skin[d], &sgr->lower_ghost[d],
      d, GKYL_LOWER_EDGE, parent, ghost);
    gkyl_skin_ghost_ranges(&sgr->upper_skin[d], &sgr->upper_ghost[d],
      d, GKYL_UPPER_EDGE, parent, ghost);
  }
}

struct gkyl_block_topo*
create_block_topo()
{
  struct gkyl_block_topo *btopo = gkyl_block_topo_new(2, 3);

  /* Block layout

     +------+
     |0     |
     |      |
     +------+-----+
     |1     |2    |
     |      |     |
     +------+-----+
    
  */  

  // block 0
  btopo->conn[0] = (struct gkyl_block_connections) {
    .connections[0] = { // x-direction connections
      { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }, // physical boundary
      { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }  // physical boundary
    },
    .connections[1] = { // y-direction connections
      { .bid = 1, .dir = 1, .edge = GKYL_UPPER_POSITIVE },
      { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL } // physical boundary
    }
  };
  // block 1
  btopo->conn[1] = (struct gkyl_block_connections) {
    .connections[0] = { // x-direction connections
      { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }, // physical boundary
      { .bid = 2, .dir = 0, .edge = GKYL_LOWER_POSITIVE }
    },
    .connections[1] = { // y-direction connections
      { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL }, // physical boundary
      { .bid = 0, .dir = 1, .edge = GKYL_LOWER_POSITIVE }
    }
  };
  // block 2
  btopo->conn[2] = (struct gkyl_block_connections) {
    .connections[0] = { // x-direction connections
      { .bid = 1, .dir = 0, .edge = GKYL_UPPER_POSITIVE },
      { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL } // physical boundary
    },
    .connections[1] = { // y-direction connections
      { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL }, // physical boundary
      { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL } // physical boundary
    }
  };

  return btopo;
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

// Block data and methods on individual blocks
struct block_data {
  struct gkyl_rect_grid grid; // grid
  gkyl_fv_proj *fv_proj; // project initial conditions

  // extended and update ranges
  struct gkyl_range ext_range, range;
  // arrays for solution
  struct gkyl_array *fdup, *f[3];

  // geometry object
  struct gkyl_wave_geom *geom;

  // equation system
  struct gkyl_wv_eqn *euler;
  // updaters
  gkyl_wave_prop *slvr[2]; // solver in each direction

  struct skin_ghost_ranges skin_ghost; // conf-space skin/ghost
  struct gkyl_array *bc_buffer; // buffer for use in block BCs

  // boundary conditions on lower/upper edges in each direction
  gkyl_wv_apply_bc *lower_bc[2], *upper_bc[2];
};

void
block_bc_updaters_init(struct block_data *bdata, const struct gkyl_block_connections *conn)
{
  int nghost[] = { 2, 2, 2 };

  // create updaters for physical boundaries (at present, all assumed
  // to be solid walls)
  for (int d=0; d<2; ++d) {
    
    bdata->lower_bc[d] = bdata->upper_bc[d] = 0;

    // create BC updater in dir 'd' on lower edge
    if (conn->connections[d][0].edge == GKYL_PHYSICAL)
      bdata->lower_bc[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->euler, bdata->geom,
        d, GKYL_LOWER_EDGE, nghost, bdata->euler->wall_bc_func, 0);

    // create BC updater in dir 'd' on upper edge
    if (conn->connections[d][1].edge == GKYL_PHYSICAL)
      bdata->upper_bc[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->euler, bdata->geom,
        d, GKYL_UPPER_EDGE, nghost, bdata->euler->wall_bc_func, 0);
  }

  // create skin/ghost region
  skin_ghost_ranges_init(&bdata->skin_ghost, &bdata->ext_range, nghost);
  // allocate buffer for inter-block BCs
  long buff_sz = 0;
  for (int d=0; d<2; ++d) {
    long vol = bdata->skin_ghost.lower_skin[d].volume;
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  bdata->bc_buffer = gkyl_array_new(GKYL_DOUBLE, 5, buff_sz);
}

void
block_bc_updaters_release(struct block_data *bdata)
{
  for (int d=0; d<2; ++d) {
    if (bdata->lower_bc[d])
      gkyl_wv_apply_bc_release(bdata->lower_bc[d]);
    if (bdata->upper_bc[d])
      gkyl_wv_apply_bc_release(bdata->upper_bc[d]);
  }
  gkyl_array_release(bdata->bc_buffer);
}

void
block_bc_updaters_apply(const struct block_data *bdata, double tm, struct gkyl_array *fld)
{
  for (int d=0; d<2; ++d) {
    if (bdata->lower_bc[d])
      gkyl_wv_apply_bc_advance(bdata->lower_bc[d], tm, &bdata->range, fld);
    if (bdata->upper_bc[d])
      gkyl_wv_apply_bc_advance(bdata->upper_bc[d], tm, &bdata->range, fld);
  }
}

void
sync_blocks(const struct gkyl_block_topo *btopo, const struct block_data bdata[],
  struct gkyl_array *fld[])
{
  for (int i=0; i<btopo->num_blocks; ++i) {
    
    for (int d=0; d<btopo->ndim; ++d) {
      const struct gkyl_target_edge *te = btopo->conn[i].connections[d];

      // lower-edge
      if (te[0].edge != GKYL_PHYSICAL) {

        struct gkyl_array *bc_buffer = bdata[i].bc_buffer;
        
        // copy skin-cell data to buffer
        gkyl_array_copy_to_buffer(
          bc_buffer->data, fld[i], &(bdata[i].skin_ghost.lower_skin[d]));

        int tbid = te[0].bid, tdir = te[0].dir;
          
        // copy buffer to ghost-cells of target block        
        switch (te[0].edge) {
          case GKYL_LOWER_POSITIVE:
          case GKYL_LOWER_NEGATIVE:
            gkyl_array_copy_from_buffer(fld[tbid], bc_buffer->data, &(bdata[tbid].skin_ghost.lower_ghost[tdir]));
            break;

          case GKYL_UPPER_POSITIVE:
          case GKYL_UPPER_NEGATIVE:
            gkyl_array_copy_from_buffer(fld[tbid], bc_buffer->data, &(bdata[tbid].skin_ghost.upper_ghost[tdir]));
            break;

          default:
            ;
        }
      }

      // upper-edge
      if (te[1].edge != GKYL_PHYSICAL) {

        struct gkyl_array *bc_buffer = bdata[i].bc_buffer;
        
        // copy skin-cell data to buffer
        gkyl_array_copy_to_buffer(
          bc_buffer->data, fld[i], &(bdata[i].skin_ghost.upper_skin[d]));

        int tbid = te[1].bid, tdir = te[1].dir;
          
        // copy buffer to ghost-cells of target block        
        switch (te[1].edge) {
          case GKYL_LOWER_POSITIVE:
          case GKYL_LOWER_NEGATIVE:
            gkyl_array_copy_from_buffer(fld[tbid], bc_buffer->data, &(bdata[tbid].skin_ghost.lower_ghost[tdir]));
            break;

          case GKYL_UPPER_POSITIVE:
          case GKYL_UPPER_NEGATIVE:
            gkyl_array_copy_from_buffer(fld[tbid], bc_buffer->data, &(bdata[tbid].skin_ghost.upper_ghost[tdir]));
            break;

          default:
            ;
        }
      }      
    }
  }
}

void
block_data_write(const char *fileNm, const struct block_data *bdata)
{ 
  gkyl_grid_sub_array_write(&bdata->grid, &bdata->range, 0, bdata->f[0], fileNm);
}

double
block_data_max_dt(const struct block_data *bdata)
{
  double dt = DBL_MAX;
  for (int d=0; d<2; ++d)
    dt = fmin(dt, gkyl_wave_prop_max_dt(bdata->slvr[d], &bdata->range, bdata->f[0]));
  return dt;
}

// job pool info for updating blocks
struct update_block_ctx {
  const struct block_data *bdata; // block data
  int bidx; // block index (for debuging)
  int dir; // direction
  double tcurr, dt; // current time and time-step
  struct gkyl_wave_prop_status stat; // status of wave propagation (on output)
};

void
update_block_job_func(void *ctx)
{
  struct update_block_ctx *ubctx = ctx;
  int d = ubctx->dir;
  double tcurr = ubctx->tcurr, dt = ubctx->dt;
  const struct block_data *bdata = ubctx->bdata;

  // run wave-prop updater
  ubctx->stat = gkyl_wave_prop_advance(bdata->slvr[d], tcurr, dt,
    &bdata->range, NULL, bdata->f[d], bdata->f[d+1]);

  // apply block-local boundary conditions
  block_bc_updaters_apply(bdata, tcurr, bdata->f[d+1]);
}

struct gkyl_update_status
update_all_blocks(const struct gkyl_job_pool *job_pool,
  const struct gkyl_block_topo *btopo, const struct block_data bdata[], double tcurr, double dt)
{
  int num_blocks = btopo->num_blocks;
  double dt_suggested = DBL_MAX;


  for (int d=0; d<2; ++d) {

    // initialize block ctx data
    struct update_block_ctx block_ctx[num_blocks];
    for (int i=0; i<num_blocks; ++i)
      block_ctx[i] = (struct update_block_ctx) {
        .bdata = &bdata[i],
        .tcurr = tcurr,
        .dir = d,
        .dt = dt,
        .bidx = i,
    };
    
    for (int i=0; i<num_blocks; ++i)
      gkyl_job_pool_add_work(job_pool, update_block_job_func, &block_ctx[i]);
    gkyl_job_pool_wait(job_pool);
  
    struct gkyl_array *fld[num_blocks];
    for (int i=0; i<num_blocks; ++i) {

      // return immediately if a block failed
      if (block_ctx[i].stat.success == 0)
        return (struct gkyl_update_status) {
          .success = 0,
          .dt_suggested = block_ctx[i].stat.dt_suggested
        };

      dt_suggested = fmin(dt_suggested, block_ctx[i].stat.dt_suggested);
      fld[i] = bdata[i].f[d+1]; // for use in block-boundary sync
    }

    // sync all block boundaries
    sync_blocks(btopo, bdata, fld);
  }
  
  return (struct gkyl_update_status) {
    .success = 1,
    .dt_suggested = dt_suggested
  };  
}

struct sim_stats {
  int nfail;
};

// context and functions for various tasks
void
init_job_func(void *ctx)
{
  struct block_data *bdata = ctx;
  gkyl_fv_proj_advance(bdata->fv_proj, 0.0, &bdata->ext_range, bdata->f[0]);
}

struct copy_job_ctx {
  int bidx;
  const struct gkyl_array *inp;
  struct gkyl_array *out;
};
void
copy_job_func(void *ctx)
{
  struct copy_job_ctx *jctx = ctx;
  //printf("Inside copy_job_func with index %d\n", jctx->bidx);
  gkyl_array_copy(jctx->out, jctx->inp);
}


// function that takes a time-step
struct gkyl_update_status
update(const struct gkyl_job_pool *job_pool,
  const struct gkyl_block_topo *btopo, const struct block_data bdata[],
  double tcurr, double dt0, struct sim_stats *stats)
{
  int num_blocks = btopo->num_blocks;
  double dt_suggested = DBL_MAX;
  
  // time-stepper states
  enum {
    UPDATE_DONE = 0,
    PRE_UPDATE,
    POST_UPDATE,
    FLUID_UPDATE,
    UPDATE_REDO,
  } state = PRE_UPDATE;

  struct copy_job_ctx copy_ctx[num_blocks];

  double dt = dt0;
  while (state != UPDATE_DONE) {
    switch (state) {
      case PRE_UPDATE:
        state = FLUID_UPDATE; // next state
          
        // copy old solution in case we need to redo this step

        // create context objects ...
        for (int i=0; i<num_blocks; ++i)
          copy_ctx[i] = (struct copy_job_ctx) {
            .bidx = i,
            .inp = bdata[i].f[0],
            .out = bdata[i].fdup
          };

        // ... run jobs
        for (int i=0; i<num_blocks; ++i)
          gkyl_job_pool_add_work(job_pool, copy_job_func, &copy_ctx[i]);
        gkyl_job_pool_wait(job_pool);

        break;
          
      case FLUID_UPDATE:
        state = POST_UPDATE; // next state
          
        struct gkyl_update_status s = update_all_blocks(job_pool, btopo, bdata, tcurr, dt);
        if (!s.success) {
          stats->nfail += 1;
          dt = s.dt_suggested;
          state = UPDATE_REDO;
          break;
        }
        dt_suggested = fmin(dt_suggested, s.dt_suggested);

        break;

      case POST_UPDATE:
        state = UPDATE_DONE;

        // copy solution in prep for next time-step

        for (int i=0; i<num_blocks; ++i)
          copy_ctx[i] = (struct copy_job_ctx) {
            .bidx = i,
            .inp = bdata[i].f[2],
            .out = bdata[i].f[0]
          };
        
        for (int i=0; i<num_blocks; ++i)
          gkyl_job_pool_add_work(job_pool, copy_job_func, &copy_ctx[i]);
        gkyl_job_pool_wait(job_pool);
          
        break;

      case UPDATE_REDO:
        state = PRE_UPDATE; // start all-over again

        // restore solution and retake step

        for (int i=0; i<num_blocks; ++i)
          copy_ctx[i] = (struct copy_job_ctx) {
            .bidx = i,
            .inp = bdata[i].fdup,
            .out = bdata[i].f[0]
          };
        
        for (int i=0; i<num_blocks; ++i)
          gkyl_job_pool_add_work(job_pool, copy_job_func, &copy_ctx[i]);
        gkyl_job_pool_wait(job_pool);
          
        break;

      case UPDATE_DONE: // unreachable code! (suppresses warning)
        break;
    }
  }

  return (struct gkyl_update_status) {
    .success = 1,
    .dt_actual = dt,
    .dt_suggested = dt_suggested,
  };
}

void
write_sol(const char *fbase, int num_blocks, const struct block_data bdata[])
{ 
  for (int i=0; i<num_blocks; ++i) {
    const char *fmt = "%s_b%d.gkyl";
    int sz = snprintf(0, 0, fmt, fbase, i);
    char fileNm[sz+1]; // ensures no buffer overflow  
    snprintf(fileNm, sizeof fileNm, fmt, fbase, i);
    block_data_write(fileNm, &bdata[i]);
  }
}

double
max_dt(int num_blocks, const struct block_data bdata[])
{
  double dt = DBL_MAX;
  for (int i=0; i<num_blocks; ++i)
    dt = fmin(dt, block_data_max_dt(&bdata[i]));
  return dt;
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
  
  int num_blocks = 3, nx = 128, ny = 128;
  struct block_data bdata[num_blocks];

  struct gkyl_job_pool *job_pool = gkyl_thread_pool_new(app_args.num_threads);

  // construct grid for each block
  gkyl_rect_grid_init(&bdata[0].grid, 2,
    (double []) { 0, 1 },
    (double []) { 1, 2 },
    (int []) { nx, ny }
  );
  gkyl_rect_grid_init(&bdata[1].grid, 2,
    (double []) { 0, 0 },
    (double []) { 1, 1 },
    (int []) { nx, ny }
  );
  gkyl_rect_grid_init(&bdata[2].grid, 2,
    (double []) { 1, 0 },
    (double []) { 2, 1 },
    (int []) { nx, ny }
  );

  // create projection updaters
  for (int i=0; i<num_blocks; ++i)
    bdata[i].fv_proj = gkyl_fv_proj_new(&bdata[i].grid, 2, 5, initFluidSod, 0);

  // create ranges and geometry
  for (int i=0; i<num_blocks; ++i) {
    gkyl_create_grid_ranges(&bdata[i].grid, (int []) { 2, 2 }, &bdata[i].ext_range, &bdata[i].range);
    bdata[i].geom = gkyl_wave_geom_new(&bdata[i].grid, &bdata[i].ext_range,
      0, 0, false);
  }

  // create FV updaters for dimensional sweeps
  for (int i=0; i<num_blocks; ++i) {
    bdata[i].euler = gkyl_wv_euler_new(1.4, false);

    for (int d=0; d<2; ++d)
      bdata[i].slvr[d] = gkyl_wave_prop_new( &(struct gkyl_wave_prop_inp) {
          .grid = &bdata[i].grid,
          .equation = bdata[i].euler,
          .limiter = GKYL_MONOTONIZED_CENTERED,
          .num_up_dirs = 1,
          .update_dirs = { d },
          .cfl = 0.95,
          .geom = bdata[i].geom
        }
      );
  }

  struct gkyl_block_topo *btopo = create_block_topo();  

  // create BC updaters
  for (int i=0; i<num_blocks; ++i)
    block_bc_updaters_init(&bdata[i], &btopo->conn[i]);

  // allocate fields
  for (int i=0; i<num_blocks; ++i) {
    bdata[i].fdup = gkyl_array_new(GKYL_DOUBLE, 5, bdata[i].ext_range.volume);
    for (int d=0; d<3; ++d) 
      bdata[i].f[d] = gkyl_array_new(GKYL_DOUBLE, 5, bdata[i].ext_range.volume);
  }

  // apply initial conditions
  for (int i=0; i<num_blocks; ++i)
    gkyl_job_pool_add_work(job_pool, init_job_func, &bdata[i]);
  gkyl_job_pool_wait(job_pool);

  // write initial conditions to file
  write_sol("euler_multiblock_0", num_blocks, bdata);

  // run simulation
  double tcurr = 0.0, tend = 0.6;
  double dt = max_dt(num_blocks, bdata);

  struct sim_stats stats = { };

  struct timespec tm_start = gkyl_wall_clock();

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    printf("Taking time-step %ld at t = %g ...", step, tcurr);
    struct gkyl_update_status status = update(job_pool, btopo, bdata, tcurr, dt, &stats);
    printf(" dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      printf("** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;

    step += 1;
  }

  double tm_total_sec = gkyl_time_diff_now_sec(tm_start);

  write_sol("euler_multiblock_1", num_blocks, bdata);

  printf("Total run-time: %g. failed steps: %d\n", tm_total_sec, stats.nfail);

  // free data
  for (int i=0; i<num_blocks; ++i) {
    gkyl_fv_proj_release(bdata[i].fv_proj);

    gkyl_wv_eqn_release(bdata[i].euler);
    block_bc_updaters_release(&bdata[i]);

    gkyl_wave_geom_release(bdata[i].geom);

    for (int d=0; d<2; ++d)
      gkyl_wave_prop_release(bdata[i].slvr[d]);
    
    gkyl_array_release(bdata[i].fdup);
    for (int d=0; d<3; ++d) 
      gkyl_array_release(bdata[i].f[d]);
  }
  
  gkyl_block_topo_release(btopo);
  gkyl_job_pool_release(job_pool);
  
  return 0;
}
