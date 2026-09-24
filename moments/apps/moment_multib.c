#include <gkyl_array_rio_priv.h>
#include <gkyl_elem_type_priv.h>
#include <gkyl_moment_multib.h>
#include <gkyl_moment_multib_priv.h>
#include <gkyl_rrobin_decomp.h>

#include <mpack.h>

// compute total number of ranges specified by cuts
static inline int
calc_cuts(int ndim, const int *cuts)
{
  int tc = 1;
  for (int d=0; d<ndim; ++d) tc *= cuts[d];
  return tc;
}

// simple linear search to check if val occurs in lst
static bool
has_int(int n, int val, const int *lst)
{
  for (int i=0; i<n; ++i)
    if (val == lst[i])
      return true;
  return false;
}

// compute total and maximum number of cuts
static void
calc_tot_and_max_cuts(const struct gkyl_block_geom *block_geom, int tot_max[2])
{
  int ndim = gkyl_block_geom_ndim(block_geom);
  int num_blocks = gkyl_block_geom_num_blocks(block_geom);

  int max_cuts = 0, tot_cuts = 0;
  for (int i=0; i<num_blocks; ++i) {
    const struct gkyl_block_geom_info *bgi = gkyl_block_geom_get_block(block_geom, i);
    int ncuts = calc_cuts(ndim, bgi->cuts);
    max_cuts = ncuts > max_cuts ? ncuts : max_cuts;
    tot_cuts += ncuts;
  }
  tot_max[0] = tot_cuts;
  tot_max[1] = max_cuts;
}

// construct the mpack meta-data for multi-block data files
static struct gkyl_msgpack_data *
moment_multib_meta(struct moment_multib_output_meta meta)
{
  struct gkyl_msgpack_data *mt = gkyl_malloc(sizeof *mt);

  mt->meta_sz = 0;
  mpack_writer_t writer;
  mpack_writer_init_growable(&writer, &mt->meta, &mt->meta_sz);

  // add some data to mpack
  mpack_build_map(&writer);
  
  mpack_write_cstr(&writer, "time");
  mpack_write_double(&writer, meta.stime);

  mpack_write_cstr(&writer, "frame");
  mpack_write_i64(&writer, meta.frame);

  mpack_write_cstr(&writer, "topo_file");
  mpack_write_cstr(&writer, meta.topo_file_name);

  mpack_write_cstr(&writer, "app_name");
  mpack_write_cstr(&writer, meta.app_name);

  mpack_complete_map(&writer);

  int status = mpack_writer_destroy(&writer);

  if (status != mpack_ok) {
    free(mt->meta); // we need to use free here as mpack does its own malloc
    gkyl_free(mt);
    mt = 0;
  }

  return mt;  
}

// write out multi-block data files
static int
moment_multib_data_write(const char *fname, struct moment_multib_output_meta meta)
{
  enum gkyl_array_rio_status status = GKYL_ARRAY_RIO_FOPEN_FAILED;
  FILE *fp = 0;
  int err;
  with_file (fp, fname, "w") {
    struct gkyl_msgpack_data *amet = moment_multib_meta(meta);
    if (amet) {
      status = gkyl_header_meta_write_fp( &(struct gkyl_array_header_info) {
          .file_type = gkyl_file_type_int[GKYL_MULTI_BLOCK_DATA_FILE],
          .meta_size = amet->meta_sz,
          .meta = amet->meta
        },
        fp
      );
      MPACK_FREE(amet->meta);
      gkyl_free(amet);
    }
    else {
      status = GKYL_ARRAY_RIO_META_FAILED;
    }
  }
  return status;
}

// construct single-block App for given block ID
static struct gkyl_moment_app *
singleb_app_new(const struct gkyl_moment_multib *mbinp, int bid,
  const struct gkyl_moment_multib_app *mbapp)
{
  int ndim = gkyl_block_geom_ndim(mbapp->block_geom);
  int num_blocks = gkyl_block_geom_num_blocks(mbapp->block_geom);

  const struct gkyl_block_geom_info *bgi =
    gkyl_block_geom_get_block(mbapp->block_geom, bid);

  // construct top-level single-block input struct
  struct gkyl_moment app_inp = { };

  strcpy(app_inp.name, mbinp->name);
  if (num_blocks > 1) {
    cstr app_name = cstr_from_fmt("%s_b%d", mbinp->name, bid);
    strcpy(app_inp.name, app_name.str);
    cstr_drop(&app_name);
  }

  app_inp.ndim = ndim;
  for (int i=0; i<ndim; ++i) {
    app_inp.lower[i] = bgi->lower[i];
    app_inp.upper[i] = bgi->upper[i];
    app_inp.cells[i] = bgi->cells[i];
  }

  app_inp.cfl_frac = mbinp->cfl_frac;

  app_inp.scheme_type = mbinp->scheme_type;
  app_inp.mp_recon = mbinp->mp_recon;
  app_inp.skip_mp_limiter = mbinp->skip_mp_limiter;
  app_inp.use_hybrid_flux_kep = mbinp->use_hybrid_flux_kep;

  app_inp.num_skip_dirs = mbinp->num_skip_dirs;
  for (int d=0; d<3; ++d)
    app_inp.skip_dirs[d] = mbinp->skip_dirs[d];

  int num_species = app_inp.num_species = mbinp->num_species;
  // construct each species input block
  for (int i=0; i<num_species; ++i) {
    const struct gkyl_moment_multib_species *sp = &mbinp->species[i];
    
    struct gkyl_moment_species species_inp = { };
    strcpy(species_inp.name, sp->name);

    species_inp.charge = sp->charge;
    species_inp.mass = sp->mass;
    
    species_inp.equation = sp->equation;
    species_inp.limiter = sp->limiter;
    species_inp.split_type = sp->split_type;

    species_inp.is_static = sp->is_static;
    species_inp.force_low_order_flux = sp->force_low_order_flux;

    // choose proper block-specific species input
    const struct gkyl_moment_multib_species_pb *sp_pb = &sp->blocks[0];
    if (!sp->duplicate_across_blocks) {
      for (int i=0; i<num_blocks; ++i)
        if (bid == sp->blocks[i].block_id) {
          sp_pb = &sp->blocks[i];
          break;
        }
    }

    species_inp.ctx = sp_pb->ctx;
    species_inp.init = sp_pb->init;

    species_inp.app_accel_ctx = sp_pb->app_accel_ctx;
    species_inp.app_accel = sp_pb->app_accel;
    species_inp.app_accel_evolve = sp_pb->app_accel_evolve;

    species_inp.nT_source_ctx = sp_pb->nT_source_ctx;
    species_inp.nT_source_func = sp_pb->nT_source_func;
    species_inp.nT_source_set_only_once = sp_pb->nT_source_set_only_once;

    // by default, skip BCs altogether
    for (int e=0; e<2; ++e) {
      species_inp.bcx[e] = GKYL_SPECIES_SKIP;
      species_inp.bcy[e] = GKYL_SPECIES_SKIP;
      species_inp.bcz[e] = GKYL_SPECIES_SKIP;
    }

    // set species physical BCs: we need to search through the list of
    // physical BCs and set the appropriate input to single-block
    // species inp
    for (int i=0; i<sp->num_physical_bcs; ++i) {
      if (bid == sp->bcs[i].bidx) {

        int e = sp->bcs[i].edge;
        if (sp->bcs[i].dir == 0)
          species_inp.bcx[e] = sp->bcs[i].bc_type;
        else if (sp->bcs[i].dir == 1)
          species_inp.bcy[e] = sp->bcs[i].bc_type;
        else
          species_inp.bcz[e] = sp->bcs[i].bc_type;
      }
    }

    // copy species input into app input
    memcpy(&app_inp.species[i], &species_inp, sizeof(struct gkyl_moment_species));
  }


  // construct field input
  if (mbinp->field.blocks) { // assumption is that blocks = 0 when there is no field
    struct gkyl_moment_field field_inp = { };
    const struct gkyl_moment_multib_field *fld = &mbinp->field;

    field_inp.epsilon0 = fld->epsilon0;
    field_inp.mu0 = fld->mu0;
    field_inp.elc_error_speed_fact = fld->elc_error_speed_fact;
    field_inp.mag_error_speed_fact = fld->mag_error_speed_fact;

    field_inp.limiter = fld->limiter;
    field_inp.is_static = fld->is_static;

    // choose proper block-specific field input
    const struct gkyl_moment_multib_field_pb *fld_pb = &fld->blocks[0];
    if (!fld->duplicate_across_blocks) {
      for (int i=0; i<num_blocks; ++i)
        if (bid == fld->blocks[i].block_id) {
          fld_pb = &fld->blocks[i];
          break;
        }
    }

    field_inp.ctx = fld_pb->ctx;
    field_inp.init = fld_pb->init;

    field_inp.app_current_ctx = fld_pb->app_current_ctx;
    field_inp.app_current = fld_pb->app_current;
    field_inp.app_current_evolve = fld_pb->app_current_evolve;
    field_inp.t_ramp_curr = fld_pb->t_ramp_curr;

    field_inp.ext_em_ctx = fld_pb->ext_em_ctx;
    field_inp.ext_em = fld_pb->ext_em;
    field_inp.ext_em_evolve = fld_pb->ext_em_evolve;
    field_inp.t_ramp_E = fld_pb->t_ramp_E;

    field_inp.use_explicit_em_coupling = fld_pb->use_explicit_em_coupling;

    // by default, skip BCs altogether
    for (int e=0; e<2; ++e) {
      field_inp.bcx[e] = GKYL_FIELD_SKIP;
      field_inp.bcy[e] = GKYL_FIELD_SKIP;
      field_inp.bcz[e] = GKYL_FIELD_SKIP;
    }

    // set field physical BCs: we need to search through the list of
    // physical BCs and set the appropriate input to single-block
    // field inp
    for (int i=0; i<fld->num_physical_bcs; ++i) {
      if (bid == fld->bcs[i].bidx) {

        int e = fld->bcs[i].edge;
        if (fld->bcs[i].dir == 0)
          field_inp.bcx[e] = fld->bcs[i].bc_type;
        else if (fld->bcs[i].dir == 1)
          field_inp.bcy[e] = fld->bcs[i].bc_type;
        else
          field_inp.bcz[e] = fld->bcs[i].bc_type;
      }
    }    

    // copy field input into app input
    memcpy(&app_inp.field, &field_inp, sizeof(struct gkyl_moment_field));
  }

  struct gkyl_comm *comm = mbapp->block_comms[bid];

  struct gkyl_app_parallelism_inp parallel_inp = {};
//  parallel_inp.use_gpu = mbinp->use_gpu; // MF 2024/09/18: moments don't have use_gpu.
  for (int d=0; d<ndim; ++d) parallel_inp.cuts[d] = bgi->cuts[d];
  parallel_inp.comm = comm;
  // copy parallelism input into app input.
  memcpy(&app_inp.parallelism, &parallel_inp, sizeof(struct gkyl_app_parallelism_inp));

  return gkyl_moment_app_new(&app_inp);
}


struct gkyl_moment_multib_app *
gkyl_moment_multib_app_new(const struct gkyl_moment_multib *mbinp)
{
  int my_rank;
  gkyl_comm_get_rank(mbinp->comm, &my_rank);
  int num_ranks;
  gkyl_comm_get_size(mbinp->comm, &num_ranks);

  int tot_max[2];
  calc_tot_and_max_cuts(mbinp->block_geom, tot_max);
  if ((num_ranks > tot_max[0]) || (num_ranks < tot_max[1]))
    return 0;

  struct gkyl_moment_multib_app *mbapp = gkyl_malloc(sizeof(*mbapp));
  strcpy(mbapp->name, mbinp->name);
  mbapp->comm = gkyl_comm_acquire(mbinp->comm);  
  
  mbapp->block_geom = gkyl_block_geom_acquire(mbinp->block_geom);
  mbapp->block_topo = gkyl_block_geom_topo(mbinp->block_geom);
  
  int ndim = gkyl_block_geom_ndim(mbapp->block_geom);
  int num_blocks = gkyl_block_geom_num_blocks(mbapp->block_geom);

  // construct round-robin decomposition
  int *branks = gkyl_malloc(sizeof(int[num_blocks]));
  for (int i=0; i<num_blocks; ++i) {
    const struct gkyl_block_geom_info *bgi = gkyl_block_geom_get_block(mbapp->block_geom, i);
    branks[i] = calc_cuts(ndim, bgi->cuts);
  }
  mbapp->round_robin = gkyl_rrobin_decomp_new(num_ranks, num_blocks, branks);

  int num_local_blocks = 0;
  mbapp->local_blocks = gkyl_malloc(sizeof(int[num_blocks]));

  int lidx = 0;
  int *rank_list = gkyl_malloc(sizeof(int[num_ranks])); // this is larger than needed

  mbapp->decomp = gkyl_malloc(num_blocks*sizeof(struct gkyl_rect_decomp*));
    
  // construct list of block communicators: there are as many
  // communicators as blocks. Not all communicators are valid on each
  // rank. The total number of valid communicators is
  // num_local_blocks.
  mbapp->block_comms = gkyl_malloc(num_blocks*sizeof(struct gkyl_comm *));
  for (int i=0; i<num_blocks; ++i) {
    gkyl_rrobin_decomp_getranks(mbapp->round_robin, i, rank_list);

    bool is_my_rank_in_decomp = has_int(branks[i], my_rank, rank_list);

    if (is_my_rank_in_decomp) {
      mbapp->local_blocks[lidx++] = i;
      num_local_blocks += 1;      
    }

    const struct gkyl_block_geom_info *bgi = gkyl_block_geom_get_block(mbapp->block_geom, i);
    struct gkyl_range block_global_range;
    gkyl_create_global_range(ndim, bgi->cells, &block_global_range);

    mbapp->decomp[i] = gkyl_rect_decomp_new_from_cuts(
      ndim, bgi->cuts, &block_global_range);

    bool status;
    mbapp->block_comms[i] = gkyl_comm_create_comm_from_ranks(mbinp->comm,
      branks[i], rank_list, mbapp->decomp[i], &status);
  }
  gkyl_free(rank_list);
  mbapp->num_local_blocks = num_local_blocks;  

  printf("Rank %d handles %d Apps\n", my_rank, num_local_blocks);
  for (int i=0; i<num_local_blocks; ++i)
    printf("  Rank %d handles block %d\n", my_rank, mbapp->local_blocks[i]);

  mbapp->num_species = 0;
  mbapp->singleb_apps = 0;
  
  if (num_local_blocks > 0) {
    mbapp->num_species = mbinp->num_species;
    mbapp->singleb_apps = gkyl_malloc(num_local_blocks*sizeof(struct gkyl_moment_app*));
  }

  if (!mbinp->field.blocks)
    mbapp->has_field = false;

  for (int i=0; i<mbinp->num_species; ++i)
    strcpy(mbapp->species_name[i], mbinp->species[i].name);

  // create single-block App for all local blocks this rank handles
  for (int i=0; i<num_local_blocks; ++i)
    mbapp->singleb_apps[i] = singleb_app_new(mbinp, mbapp->local_blocks[i], mbapp);

  // construct send/recv lists for each block this rank handles
  mbapp->send_conn = gkyl_malloc(sizeof(struct gkyl_multib_comm_conn *[num_local_blocks]));
  mbapp->recv_conn = gkyl_malloc(sizeof(struct gkyl_multib_comm_conn *[num_local_blocks]));

  mbapp->stat = (struct gkyl_moment_stat) {
  };

  gkyl_free(branks);
  
  return mbapp;
}

double
gkyl_moment_multib_app_max_dt(gkyl_moment_multib_app *app)
{
  // TODO
  return 0;
}


void
gkyl_moment_multib_app_apply_ic(gkyl_moment_multib_app* app, double t0)
{
  app->tcurr = t0;
  gkyl_moment_multib_app_apply_ic_field(app, t0);
  for (int i=0; i<app->num_species; ++i)
    gkyl_moment_multib_app_apply_ic_species(app, i, t0);
}

void
gkyl_moment_multib_app_apply_ic_field(gkyl_moment_multib_app* app, double t0)
{
  app->tcurr = t0;
  for (int i=0; i<app->num_local_blocks; ++i)
    gkyl_moment_app_apply_ic_field(app->singleb_apps[i], t0);
  gkyl_comm_barrier(app->comm);
}

void
gkyl_moment_multib_app_apply_ic_species(gkyl_moment_multib_app* app, int sidx, double t0)
{
  app->tcurr = t0;
  for (int i=0; i<app->num_local_blocks; ++i)
    gkyl_moment_app_apply_ic_species(app->singleb_apps[i], sidx, t0);
  gkyl_comm_barrier(app->comm);
}

struct gkyl_app_restart_status
gkyl_moment_multib_app_from_frame_field(gkyl_moment_multib_app *app,
  int frame)
{
  // TODO
  return (struct gkyl_app_restart_status) { };
}

struct gkyl_app_restart_status
gkyl_moment_multib_app_from_frame_species(gkyl_moment_multib_app *app,
  int sidx, int frame)
{
  // TODO
  return (struct gkyl_app_restart_status) { };  
}

// private function to handle variable argument list for printing
static void
v_moment_app_cout(const gkyl_moment_multib_app* app, FILE *fp, const char *fmt, va_list argp)
{
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if ((rank == 0) && fp)
    vfprintf(fp, fmt, argp);
}

void
gkyl_moment_multib_app_cout(const gkyl_moment_multib_app* app, FILE *fp, const char *fmt, ...)
{
  va_list argp;
  va_start(argp, fmt);
  v_moment_app_cout(app, fp, fmt, argp);
  va_end(argp);
}

void
gkyl_moment_multib_app_write_topo(const gkyl_moment_multib_app* app)
{
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (0 == rank) {
    cstr file_name = cstr_from_fmt("%s_btopo.gkyl", app->name);
    gkyl_block_topo_write(app->block_topo, file_name.str);
    cstr_drop(&file_name);
  }
}

void
gkyl_moment_multib_app_write(const gkyl_moment_multib_app* app, double tm, int frame)
{
  gkyl_moment_multib_app_write_field(app, tm, frame);
  for (int i=0; i<app->num_species; ++i)
    gkyl_moment_multib_app_write_species(app, i, tm, frame);
}

void
gkyl_moment_multib_app_write_field(const gkyl_moment_multib_app *app, double tm, int frame)
{
  for (int i=0; i<app->num_local_blocks; ++i)
    gkyl_moment_app_write_field(app->singleb_apps[i], tm, frame);

  if (app->has_field) {
    int rank;
    gkyl_comm_get_rank(app->comm, &rank);
    if (0 == rank) {
      cstr file_name = cstr_from_fmt("%s-%s_%d.gkyl", app->name, "field", frame);
      cstr topo_file_name = cstr_from_fmt("%s_btopo.gkyl", app->name);
      
      moment_multib_data_write(file_name.str, (struct moment_multib_output_meta) {
          .frame = frame,
          .stime = tm,
          .topo_file_name = topo_file_name.str,
          .app_name = app->name
        }
      );
      
      cstr_drop(&topo_file_name);
      cstr_drop(&file_name);
    }
  }
  
  gkyl_comm_barrier(app->comm);
}

void
gkyl_moment_multib_app_write_species(const gkyl_moment_multib_app* app, int sidx, double tm, int frame)
{
  for (int i=0; i<app->num_local_blocks; ++i)
    gkyl_moment_app_write_species(app->singleb_apps[i], sidx, tm, frame);

  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (0 == rank) {
    cstr file_name = cstr_from_fmt("%s-%s_%d.gkyl", app->name, app->species_name[sidx], frame);
    cstr topo_file_name = cstr_from_fmt("%s_btopo.gkyl", app->name);
      
    moment_multib_data_write(file_name.str, (struct moment_multib_output_meta) {
        .frame = frame,
        .stime = tm,
        .topo_file_name = topo_file_name.str,
        .app_name = app->name
      }
    );
    
    cstr_drop(&topo_file_name);
    cstr_drop(&file_name);
  }  
  
  gkyl_comm_barrier(app->comm);
}

void
gkyl_moment_multib_app_write_field_energy(gkyl_moment_multib_app *app)
{
  // TODO
}

void
gkyl_moment_multib_app_write_integrated_mom(gkyl_moment_multib_app *app)
{
  // TODO
}

void
gkyl_moment_multib_app_stat_write(const gkyl_moment_multib_app *app)
{
  // TODO
}

struct gkyl_update_status
gkyl_moment_multib_update(gkyl_moment_multib_app *app, double dt)
{
  // TODO
  return (struct gkyl_update_status) { };
}

void
gkyl_moment_multib_app_calc_field_energy(gkyl_moment_multib_app *app, double tm)
{
  // TODO  
}

void
gkyl_moment_multib_app_get_field_energy(gkyl_moment_multib_app *app, double *vals)
{
  // TODO
}

void
gkyl_moment_multib_app_calc_integrated_mom(gkyl_moment_multib_app *app, double tm)
{
  // TODO
}

struct gkyl_moment_stat
gkyl_moment_multib_app_stat(gkyl_moment_multib_app *app)
{
  return app->stat;
}

void
gkyl_moment_multib_app_release(gkyl_moment_multib_app* mbapp)
{
  if (mbapp->singleb_apps) {
    for (int i=0; i<mbapp->num_local_blocks; ++i)
      gkyl_moment_app_release(mbapp->singleb_apps[i]);
    gkyl_free(mbapp->singleb_apps);

    for (int i=0; i<mbapp->num_local_blocks; ++i) {
      //gkyl_multib_comm_conn_release(mbapp->send_conn[i]);
      //gkyl_multib_comm_conn_release(mbapp->recv_conn[i]);
    }

    gkyl_free(mbapp->send_conn);
    gkyl_free(mbapp->recv_conn);
  }


  int num_blocks = gkyl_block_geom_num_blocks(mbapp->block_geom);

  for (int i=0; i<num_blocks; ++i)
    gkyl_rect_decomp_release(mbapp->decomp[i]);
  gkyl_free(mbapp->decomp);

  for (int i=0; i<num_blocks; ++i)
    gkyl_comm_release(mbapp->block_comms[i]);
  gkyl_free(mbapp->block_comms);
  gkyl_comm_release(mbapp->comm);

  gkyl_rrobin_decomp_release(mbapp->round_robin);
  
  gkyl_block_geom_release(mbapp->block_geom);
  gkyl_block_topo_release(mbapp->block_topo);
  gkyl_free(mbapp->local_blocks);    
  
  gkyl_free(mbapp);
}
