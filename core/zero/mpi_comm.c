#ifdef GKYL_HAVE_MPI

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_rio_format_desc.h>
#include <gkyl_array_rio_priv.h>
#include <gkyl_comm_priv.h>
#include <gkyl_elem_type_priv.h>
#include <gkyl_mpi_comm.h>
#include <gkyl_mpi_comm_priv.h>

#include <assert.h>
#include <errno.h>
#include <string.h>

// Mapping of Gkeyll type to MPI_Datatype
static MPI_Datatype g2_mpi_datatype[] = {
  [GKYL_INT] = MPI_INT,
  [GKYL_INT_64] = MPI_INT64_T,
  [GKYL_FLOAT] = MPI_FLOAT,
  [GKYL_DOUBLE] = MPI_DOUBLE
};

// Mapping of Gkeyll ops to MPI_Op
static MPI_Op g2_mpi_op[] = {
  [GKYL_MIN] = MPI_MIN,
  [GKYL_MAX] = MPI_MAX,
  [GKYL_SUM] = MPI_SUM
};

struct gkyl_comm_state {
  MPI_Request req;
  MPI_Status stat;
};

struct extra_mpi_comm_inp {
  bool is_comm_allocated; // is MPI_Comm allocated?
};

// Internal method to create a new MPI communicator
static struct gkyl_comm* mpi_comm_new(
  const struct gkyl_mpi_comm_inp *inp, const struct extra_mpi_comm_inp *extra_inp);

static void
comm_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_comm *comm = container_of(ref, struct gkyl_comm, ref_count);
  struct mpi_comm *mpi = container_of(comm, struct mpi_comm, priv_comm.pub_comm);

  int ndim = mpi->decomp->ndim;
  gkyl_rect_decomp_release(mpi->decomp);

  gkyl_rect_decomp_neigh_release(mpi->neigh);
  for (int d=0; d<ndim; ++d)
    gkyl_rect_decomp_neigh_release(mpi->per_neigh[d]);
  
  for (int i=0; i<MAX_RECV_NEIGH; ++i)
    gkyl_mem_buff_release(mpi->recv[i].buff);
  
  for (int i=0; i<MAX_RECV_NEIGH; ++i)
    gkyl_mem_buff_release(mpi->send[i].buff);
  
  gkyl_mem_buff_release(mpi->allgather_buff_local.buff);
  gkyl_mem_buff_release(mpi->allgather_buff_global.buff);

  if (mpi->is_mcomm_allocated)
    MPI_Comm_free(&mpi->mcomm);

  gkyl_free(mpi);
}

static int
get_rank(struct gkyl_comm *comm, int *rank)
{
  struct mpi_comm *mpi = container_of(comm, struct mpi_comm, priv_comm.pub_comm);
  MPI_Comm_rank(mpi->mcomm, rank);
  return 0;
}

static int
get_size(struct gkyl_comm *comm, int *sz)
{
  struct mpi_comm *mpi = container_of(comm, struct mpi_comm, priv_comm.pub_comm);
  MPI_Comm_size(mpi->mcomm, sz);
  return 0;
}

static int
array_send(struct gkyl_array *array, int dest, int tag, struct gkyl_comm *comm)
{
  size_t vol = array->esznc*array->size;
  struct mpi_comm *mpi = container_of(comm, struct mpi_comm, priv_comm.pub_comm);  
  int ret = MPI_Send(array->data, vol, MPI_CHAR, dest, tag, mpi->mcomm); 
  return ret == MPI_SUCCESS ? 0 : 1;
}

static int
array_isend(struct gkyl_array *array, int dest, int tag, struct gkyl_comm *comm, struct gkyl_comm_state *state)
{
  size_t vol = array->esznc*array->size;
  struct mpi_comm *mpi = container_of(comm, struct mpi_comm, priv_comm.pub_comm);  
  int ret = MPI_Isend(array->data, vol, MPI_CHAR, dest, tag, mpi->mcomm, &state->req); 
  return ret == MPI_SUCCESS ? 0 : 1;
}

static int
array_recv(struct gkyl_array *array, int src, int tag, struct gkyl_comm *comm)
{
  size_t vol = array->esznc*array->size;
  struct mpi_comm *mpi = container_of(comm, struct mpi_comm, priv_comm.pub_comm);  
  MPI_Status stat;
  int ret = MPI_Recv(array->data, vol, MPI_CHAR, src, tag, mpi->mcomm, &stat); 
  return ret == MPI_SUCCESS ? 0 : 1;
}

static int
array_irecv(struct gkyl_array *array, int src, int tag, struct gkyl_comm *comm, struct gkyl_comm_state *state)
{
  size_t vol = array->esznc*array->size;
  struct mpi_comm *mpi = container_of(comm, struct mpi_comm, priv_comm.pub_comm);  
  int ret = MPI_Irecv(array->data, vol, MPI_CHAR, src, tag, mpi->mcomm, &state->req); 
  return ret == MPI_SUCCESS ? 0 : 1;
}

static int
allreduce(struct gkyl_comm *comm, enum gkyl_elem_type type,
  enum gkyl_array_op op, int nelem, const void *inp,
  void *out)
{
  struct mpi_comm *mpi = container_of(comm, struct mpi_comm, priv_comm.pub_comm);  
  int ret =
    MPI_Allreduce(inp, out, nelem, g2_mpi_datatype[type], g2_mpi_op[op], mpi->mcomm);
  return ret == MPI_SUCCESS ? 0 : 1;
}

static int
array_allgather(struct gkyl_comm *comm,
  const struct gkyl_range *local, const struct gkyl_range *global, 
  const struct gkyl_array *array_local, struct gkyl_array *array_global)
{
  assert(array_global->esznc == array_local->esznc);

  struct mpi_comm *mpi = container_of(comm, struct mpi_comm, priv_comm.pub_comm);

  struct gkyl_range gather_range;

  int rank;
  MPI_Comm_rank(mpi->mcomm, &rank);

  assert(local->volume == mpi->decomp->ranges[rank].volume);
  assert(global->volume == mpi->decomp->parent_range.volume);

  // potentially re-size local buffer volume
  size_t send_vol = array_local->esznc*mpi->decomp->ranges[rank].volume;
  if (gkyl_mem_buff_size(mpi->allgather_buff_local.buff) < send_vol)
    gkyl_mem_buff_resize(mpi->allgather_buff_local.buff, send_vol);

  // potentially re-size global buffer volume
  size_t buff_global_vol = array_local->esznc*mpi->decomp->parent_range.volume;
  if (gkyl_mem_buff_size(mpi->allgather_buff_global.buff) < buff_global_vol)
    gkyl_mem_buff_resize(mpi->allgather_buff_global.buff, buff_global_vol);

  // copy data to local buffer
  gkyl_array_copy_to_buffer(gkyl_mem_buff_data(mpi->allgather_buff_local.buff), 
    array_local, local);

  size_t nelem = array_local->esznc*mpi->decomp->ranges[rank].volume;
  // gather data into global buffer
  int ret = 
    MPI_Allgather(gkyl_mem_buff_data(mpi->allgather_buff_local.buff), nelem, MPI_CHAR, 
      gkyl_mem_buff_data(mpi->allgather_buff_global.buff), nelem, MPI_CHAR, mpi->mcomm);

  // copy data to global array
  int idx = 0;
  for (int r=0; r<mpi->decomp->ndecomp; ++r) {
    int isrecv = gkyl_sub_range_intersect(
      &gather_range, global, &mpi->decomp->ranges[r]);
    gkyl_array_copy_from_buffer(array_global, 
      gkyl_mem_buff_data(mpi->allgather_buff_global.buff) + idx, &gather_range);
    idx += array_local->esznc*gather_range.volume;
  }

  return 0;
}

static int
array_bcast(struct gkyl_comm *comm, const struct gkyl_array *asend,
  struct gkyl_array *arecv, int root)
{
  assert(asend->esznc == arecv->esznc);
  assert(asend->size == arecv->size);

  struct mpi_comm *mpi = container_of(comm, struct mpi_comm, priv_comm.pub_comm);

  size_t nelem = asend->esznc*asend->size;
  int ret = 
    MPI_Bcast(asend->data, nelem, MPI_CHAR, root, mpi->mcomm);

  return 0;
}

static int
sync(struct gkyl_comm *comm,
  const struct gkyl_range *local, const struct gkyl_range *local_ext,
  struct gkyl_array *array, bool use_corners)
{
  struct mpi_comm *mpi = container_of(comm, struct mpi_comm, priv_comm.pub_comm);

  int rank; gkyl_comm_get_rank(comm, &rank);

  int nghost[GKYL_MAX_DIM] = { 0 };
  for (int i=0; i<mpi->decomp->ndim; ++i)
    nghost[i] = local_ext->upper[i]-local->upper[i];

  int nridx = 0;
  int tag = MPI_BASE_TAG;

  // post nonblocking recv to get data into ghost-cells  
  for (int n=0; n<mpi->neigh->num_neigh; ++n) {
    int nid = mpi->neigh->neigh[n];
    int n_dir = mpi->neigh->dir[n];
    int n_edge = mpi->neigh->edge[n];
    
    struct gkyl_range skin;
    if (use_corners)
      gkyl_skin_ghost_with_corners_ranges(&skin, &mpi->recv[nridx].range, n_dir, n_edge,
        local_ext, nghost);
    else
      gkyl_skin_ghost_ranges(&skin, &mpi->recv[nridx].range, n_dir, n_edge,
        local_ext, nghost);
    
    size_t recv_vol = array->esznc*mpi->recv[nridx].range.volume;

    if (recv_vol>0) {
      if (gkyl_mem_buff_size(mpi->recv[nridx].buff) < recv_vol)
        gkyl_mem_buff_resize(mpi->recv[nridx].buff, recv_vol);

      MPI_Irecv(gkyl_mem_buff_data(mpi->recv[nridx].buff),
        recv_vol, MPI_CHAR, nid, tag, mpi->mcomm, &mpi->recv[nridx].status);

      nridx += 1;
    }
  }

  int nsidx = 0;
  
  // post non-blocking sends of skin-cell data to neighbors
  for (int n=0; n<mpi->neigh->num_neigh; ++n) {
    int nid = mpi->neigh->neigh[n];
    int n_dir = mpi->neigh->dir[n];
    int n_edge = mpi->neigh->edge[n];    
    
    struct gkyl_range ghost;
    if (use_corners)
      gkyl_skin_ghost_with_corners_ranges(&mpi->send[nsidx].range, &ghost, n_dir, n_edge,
        local_ext, nghost);
    else
      gkyl_skin_ghost_ranges(&mpi->send[nsidx].range, &ghost, n_dir, n_edge,
        local_ext, nghost);
    
    size_t send_vol = array->esznc*mpi->send[nsidx].range.volume;

    if (send_vol>0) {
      if (gkyl_mem_buff_size(mpi->send[nsidx].buff) < send_vol)
        gkyl_mem_buff_resize(mpi->send[nsidx].buff, send_vol);
      
      gkyl_array_copy_to_buffer(gkyl_mem_buff_data(mpi->send[nsidx].buff),
        array, &(mpi->send[nsidx].range));

      MPI_Isend(gkyl_mem_buff_data(mpi->send[nsidx].buff),
        send_vol, MPI_CHAR, nid, tag, mpi->mcomm, &mpi->send[nsidx].status);

      nsidx += 1;
    }
  }

  // complete send
  for (int s=0; s<nsidx; ++s) {
    int issend = mpi->send[s].range.volume;
    if (issend)
      MPI_Wait(&mpi->send[s].status, MPI_STATUS_IGNORE);
  }

  // complete recv, copying data into ghost-cells
  for (int r=0; r<nridx; ++r) {
    int isrecv = mpi->recv[r].range.volume;
    if (isrecv) {
      MPI_Wait(&mpi->recv[r].status, MPI_STATUS_IGNORE);
      
      gkyl_array_copy_from_buffer(array,
        gkyl_mem_buff_data(mpi->recv[r].buff),
        &(mpi->recv[r].range)
      );
    }
  }
  
  return 0;
}

static int
array_sync(struct gkyl_comm *comm,
  const struct gkyl_range *local, const struct gkyl_range *local_ext,
  struct gkyl_array *array)
{
  struct mpi_comm *mpi = container_of(comm, struct mpi_comm, priv_comm.pub_comm);  
  sync(comm, local, local_ext, array, false);
  if (mpi->sync_corners) {
    for (int i=1; i<mpi->decomp->ndim; ++i)
      sync(comm, local, local_ext, array, true);
  }
  
  return 0;
}

static int
per_send_tag(const struct gkyl_range *dir_edge,
  int dir, int e)
{
  int base_tag = MPI_BASE_PER_TAG;
  return base_tag + gkyl_range_idx(dir_edge, (int[]) { dir, e });
}
static int
per_recv_tag(const struct gkyl_range *dir_edge,
  int dir, int e)
{
  int base_tag = MPI_BASE_PER_TAG;
  return base_tag + gkyl_range_idx(dir_edge, (int[]) { dir, (e+1)%2 });
}

static int
per_sync(struct gkyl_comm *comm, const struct gkyl_range *local,
  const struct gkyl_range *local_ext,
  int nper_dirs, const int *per_dirs, struct gkyl_array *array, bool use_corners)
{
  struct mpi_comm *mpi = container_of(comm, struct mpi_comm, priv_comm.pub_comm);

  if (!mpi->touches_any_edge) return 0; // nothing to sync

  int nghost[GKYL_MAX_DIM];
  for (int i=0; i<mpi->decomp->ndim; ++i)
    nghost[i] = local_ext->upper[i]-local->upper[i];

  int nridx = 0;
  int edge_type[] = { GKYL_LOWER_EDGE, GKYL_UPPER_EDGE };

  // post nonblocking recv to get data into ghost-cells
  for (int i=0; i<nper_dirs; ++i) {
    int dir = per_dirs[i];

    for (int e=0; e<2; ++e) {
      if (mpi->is_on_edge[e][dir]) {

        for (int pn=0; pn<mpi->per_neigh[dir]->num_neigh; ++pn) {
          int nid = mpi->per_neigh[dir]->neigh[pn];

          struct gkyl_range skin;
          if (use_corners)
            gkyl_skin_ghost_with_corners_ranges(&skin, &mpi->recv[nridx].range, dir, edge_type[e],
              local_ext, nghost);
          else
            gkyl_skin_ghost_ranges(&skin, &mpi->recv[nridx].range, dir, edge_type[e],
              local_ext, nghost);

          size_t recv_vol = array->esznc*mpi->recv[nridx].range.volume;

          if (recv_vol>0) {
            if (gkyl_mem_buff_size(mpi->recv[nridx].buff) < recv_vol)
              gkyl_mem_buff_resize(mpi->recv[nridx].buff, recv_vol);

            int tag = per_recv_tag(&mpi->dir_edge, dir, e);
            MPI_Irecv(gkyl_mem_buff_data(mpi->recv[nridx].buff),
              recv_vol, MPI_CHAR, nid, tag, mpi->mcomm, &mpi->recv[nridx].status);

            nridx += 1;
          }
        }
      }
    }
  }

  int nsidx = 0;  

  // post non-blocking send skin-cell data to neighbors
  for (int i=0; i<nper_dirs; ++i) {
    int dir = per_dirs[i];

    for (int e=0; e<2; ++e) {
      if (mpi->is_on_edge[e][dir]) {

        for (int pn=0; pn<mpi->per_neigh[dir]->num_neigh; ++pn) {
          int nid = mpi->per_neigh[dir]->neigh[pn];

          struct gkyl_range ghost;
          if (use_corners)
            gkyl_skin_ghost_with_corners_ranges(&mpi->send[nsidx].range, &ghost, dir, edge_type[e],
              local_ext, nghost);
          else
            gkyl_skin_ghost_ranges(&mpi->send[nsidx].range, &ghost, dir, edge_type[e],
              local_ext, nghost);

          size_t send_vol = array->esznc*mpi->send[nsidx].range.volume;

          if (send_vol>0) {
            if (gkyl_mem_buff_size(mpi->send[nsidx].buff) < send_vol)
              gkyl_mem_buff_resize(mpi->send[nsidx].buff, send_vol);
            
            gkyl_array_copy_to_buffer(gkyl_mem_buff_data(mpi->send[nsidx].buff),
              array, &(mpi->send[nsidx].range));

            int tag = per_send_tag(&mpi->dir_edge, dir, e);
            MPI_Isend(gkyl_mem_buff_data(mpi->send[nsidx].buff),
              send_vol, MPI_CHAR, nid, tag, mpi->mcomm, &mpi->send[nsidx].status);

            nsidx += 1;
          }
        }
      }
    }
  }

  // complete send
  for (int s=0; s<nsidx; ++s) {
    int issend = mpi->send[s].range.volume;
    if (issend)
      MPI_Wait(&mpi->send[s].status, MPI_STATUS_IGNORE);
  }

  // complete recv, copying data into ghost-cells
  for (int r=0; r<nridx; ++r) {
    int isrecv = mpi->recv[r].range.volume;
    if (isrecv) {
      MPI_Wait(&mpi->recv[r].status, MPI_STATUS_IGNORE);
      
      gkyl_array_copy_from_buffer(array,
        gkyl_mem_buff_data(mpi->recv[r].buff),
        &(mpi->recv[r].range)
      );
    }
  }

  mpi->nrecv = nridx > mpi->nrecv ? nridx : mpi->nrecv;
  
  return 0;
}

static int
array_per_sync(struct gkyl_comm *comm, const struct gkyl_range *local,
  const struct gkyl_range *local_ext,
  int nper_dirs, const int *per_dirs, struct gkyl_array *array)
{
  struct mpi_comm *mpi = container_of(comm, struct mpi_comm, priv_comm.pub_comm);  
  per_sync(comm, local, local_ext, nper_dirs, per_dirs, array, false);
  if (mpi->sync_corners) {
    for (int i=1; i<mpi->decomp->ndim; ++i)
      per_sync(comm, local, local_ext, nper_dirs, per_dirs, array, true);
  }
  
  return 0;
}

static int
barrier(struct gkyl_comm *comm)
{
  struct mpi_comm *mpi = container_of(comm, struct mpi_comm, priv_comm.pub_comm);
  MPI_Barrier(mpi->mcomm);
  return 0;
}

// set of functions to help with parallel array output using MPI-IO
static void
sub_array_decomp_write(struct mpi_comm *comm,
  const struct gkyl_rect_decomp *decomp,
  const struct gkyl_range *range,
  const struct gkyl_msgpack_data *meta,
  const struct gkyl_array *arr, MPI_File fp)
{
#define _F(loc) gkyl_array_cfetch(arr, loc)

  int rank;
  MPI_Comm_rank(comm->mcomm, &rank);

  // seek to appropriate place in the file, depending on rank
  size_t hdr_sz = gkyl_base_hdr_size(meta->meta_sz) + gkyl_file_type_3_hrd_size(range->ndim);
  size_t file_loc = hdr_sz +
    arr->esznc*comm->local_range_offset +
    rank*gkyl_file_type_3_range_hrd_size(range->ndim);

  MPI_Offset fp_offset = file_loc;
  MPI_File_seek(fp, fp_offset, MPI_SEEK_SET);

  do {
    char *buff; size_t buff_sz;
    FILE *fbuff = open_memstream(&buff, &buff_sz);
  
    uint64_t loidx[GKYL_MAX_DIM] = {0}, upidx[GKYL_MAX_DIM] = {0};
    for (int d = 0; d < range->ndim; ++d) {
      loidx[d] = range->lower[d];
      upidx[d] = range->upper[d];
    }
    
    fwrite(loidx, sizeof(uint64_t), range->ndim, fbuff);
    fwrite(upidx, sizeof(uint64_t), range->ndim, fbuff);
    uint64_t sz = range->volume;
    fwrite(&sz, sizeof(uint64_t), 1, fbuff);

    fclose(fbuff);

    MPI_Status status;
    MPI_File_write(fp, buff, buff_sz, MPI_CHAR, &status);

    free(buff);
    
  } while (0);

  // construct skip iterator to allow writing (potentially) in chunks
  // rather than element by element or requiring a copy of data. Note:
  // We must use "range" here and not decomp->ranges[rank] as the
  // latter is not a sub-range of the local extended
  // range.
  struct gkyl_range_skip_iter skip;
  gkyl_range_skip_iter_init(&skip, range);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &skip.range);

  MPI_Status status;
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&skip.range, iter.idx);
    MPI_File_write(fp, _F(start), arr->esznc*skip.delta, MPI_CHAR, &status);
  }

#undef _F
}

static int
grid_sub_array_decomp_write_fp(struct mpi_comm *comm,
  const struct gkyl_rect_grid *grid,
  const struct gkyl_rect_decomp *decomp,
  const struct gkyl_range *range,
  const struct gkyl_msgpack_data *meta,
  const struct gkyl_array *arr, MPI_File fp)
{
  char *buff; size_t buff_sz;
  FILE *fbuff = open_memstream(&buff, &buff_sz);

  // write header to a char buffer
  gkyl_grid_sub_array_header_write_fp(grid,
    &(struct gkyl_array_header_info) {
      .file_type = gkyl_file_type_int[GKYL_MULTI_RANGE_DATA_FILE],
      .etype = arr->type,
      .esznc = arr->esznc,
      .tot_cells = decomp->parent_range.volume,
      .meta_size = meta ? meta->meta_sz : 0,
      .meta = meta ? meta->meta : 0
    },
    fbuff
  );
  uint64_t nrange = decomp->ndecomp;
  fwrite(&nrange, sizeof(uint64_t), 1, fbuff);
  fclose(fbuff);

  // actual header IO is only done by rank 0
  int rank;
  MPI_Comm_rank(comm->mcomm, &rank);
  if (rank == 0) {
    MPI_Status status;
    MPI_File_write(fp, buff, buff_sz, MPI_CHAR, &status);
  }
  free(buff);

  struct gkyl_msgpack_data zero_meta = (struct gkyl_msgpack_data) {
    .meta_sz = 0,
    .meta = 0
  };
  
  // write data in array
  sub_array_decomp_write(comm, decomp, range,
    meta ? meta : &zero_meta,
    arr, fp);
  return errno;
}

static int array_write(struct gkyl_comm *comm,
  const struct gkyl_rect_grid *grid,
  const struct gkyl_range *range,
  const struct gkyl_msgpack_data *meta,
  const struct gkyl_array *arr, const char *fname)
{
  struct mpi_comm *mpi = container_of(comm, struct mpi_comm, priv_comm.pub_comm);
  MPI_File fp;
  int err =
    MPI_File_open(mpi->mcomm, fname, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
  if (err != MPI_SUCCESS)
    return err;
  err = grid_sub_array_decomp_write_fp(mpi, grid, mpi->decomp, range, meta, arr, fp);
  MPI_File_close(&fp);
  return err;
}

static int
array_read(struct gkyl_comm *comm,
  const struct gkyl_rect_grid *grid, const struct gkyl_range *range,
  struct gkyl_array *arr, const char *fname)
{
  struct gkyl_rect_grid fgrid;
  int status = gkyl_grid_sub_array_read(&fgrid, range, arr, fname);
  if (status == 0) {
    if (!gkyl_rect_grid_cmp(grid, &fgrid))
      status = 1;
  }
  return status;
}

static struct gkyl_comm*
extend_comm(const struct gkyl_comm *comm, const struct gkyl_range *erange)
{
  struct mpi_comm *mpi = container_of(comm, struct mpi_comm, priv_comm.pub_comm);

  // extend internal decomp object and create a new communicator
  struct gkyl_rect_decomp *ext_decomp = gkyl_rect_decomp_extended_new(erange, mpi->decomp);
  struct gkyl_comm *ext_comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
      .mpi_comm = mpi->mcomm,
      .decomp = ext_decomp
    }
  );
  gkyl_rect_decomp_release(ext_decomp);
  return ext_comm;
}

static struct gkyl_comm*
split_comm(const struct gkyl_comm *comm, int color, struct gkyl_rect_decomp *new_decomp)
{
  struct mpi_comm *mpi = container_of(comm, struct mpi_comm, priv_comm.pub_comm);
  int rank;
  MPI_Comm_rank(mpi->mcomm, &rank);
  MPI_Comm new_mcomm;
  int ret = MPI_Comm_split(mpi->mcomm, color, rank, &new_mcomm);
  assert(ret == MPI_SUCCESS);

  return gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
      .mpi_comm = new_mcomm,
      .sync_corners = mpi->sync_corners,
      .decomp = new_decomp,
    }
  );
}

static struct gkyl_comm*
create_comm_from_ranks(const struct gkyl_comm *comm,
  int nranks, const int *ranks, struct gkyl_rect_decomp *new_decomp,
  bool *is_valid)
{
  struct mpi_comm *mpi = container_of(comm, struct mpi_comm, priv_comm.pub_comm);

  MPI_Group group;
  MPI_Comm_group(mpi->mcomm, &group);

  MPI_Group new_group;
  MPI_Group_incl(group, nranks, ranks, &new_group);

  MPI_Comm new_mcomm;
  MPI_Comm_create_group(mpi->mcomm, new_group, 0, &new_mcomm);

  *is_valid = false;
  struct gkyl_comm *new_comm = 0;
  if (MPI_COMM_NULL != new_mcomm) {
    *is_valid = true;

    new_comm = mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
        .mpi_comm = new_mcomm,
        .sync_corners = mpi->sync_corners,
        .decomp = new_decomp,
      },
      &(struct extra_mpi_comm_inp) {
        .is_comm_allocated = true
      }
    );
  }

  MPI_Group_free(&group);
  MPI_Group_free(&new_group);
  
  return new_comm;
}

static struct gkyl_comm_state *
comm_state_new(struct gkyl_comm *comm)
{
  struct gkyl_comm_state *state = gkyl_malloc(sizeof *state);
  return state;
}

static void
comm_state_release(struct gkyl_comm_state *state)
{
  gkyl_free(state);
}

static void
comm_state_wait(struct gkyl_comm_state *state)
{
  MPI_Wait(&state->req, &state->stat);
}

static struct gkyl_comm*
mpi_comm_new(const struct gkyl_mpi_comm_inp *inp,
  const struct extra_mpi_comm_inp *extra_inp)
{
  struct mpi_comm *mpi = gkyl_malloc(sizeof *mpi);
  strcpy(mpi->priv_comm.pub_comm.id, "mpi_comm");

  mpi->is_mcomm_allocated = extra_inp->is_comm_allocated;
  mpi->mcomm = inp->mpi_comm;
  mpi->sync_corners = inp->sync_corners;

  int comm_size;
  MPI_Comm_size(inp->mpi_comm, &comm_size);

  mpi->priv_comm.pub_comm.has_decomp = true;
  if (0 == inp->decomp) {
    mpi->priv_comm.pub_comm.has_decomp = false;
    
    // construct a dummy decomposition
    mpi->decomp =
      gkyl_rect_decomp_new_from_cuts_and_cells(1,
        (int[]) { comm_size }, (int[]) { comm_size });
  }
  else {
    mpi->decomp = gkyl_rect_decomp_acquire(inp->decomp);
  }

  int rank;
  MPI_Comm_rank(inp->mpi_comm, &rank);

  mpi->local_range_offset = gkyl_rect_decomp_calc_offset(mpi->decomp, rank);

  // NOTE: we are not computing corner neighbors as the corner syncs
  // are handled by two calls to sync method instead
  mpi->neigh = gkyl_rect_decomp_calc_neigh(mpi->decomp, false, rank);
    
  for (int d=0; d<mpi->decomp->ndim; ++d)
    // NOTE: we are not computing corner periodic neighbors as the
    // corner syncs are handled by two calls to periodic sync method
    // instead
    mpi->per_neigh[d] =
      gkyl_rect_decomp_calc_periodic_neigh(mpi->decomp, d, false, rank);
  
  gkyl_range_init(&mpi->dir_edge, 2, (int[]) { 0, 0 }, (int[]) { GKYL_MAX_DIM, 2 });
  
  int num_touches = 0;
  for (int d=0; d<mpi->decomp->ndim; ++d) {
    mpi->is_on_edge[0][d] = gkyl_range_is_on_lower_edge(
      d, &mpi->decomp->ranges[rank], &mpi->decomp->parent_range);
    mpi->is_on_edge[1][d] = gkyl_range_is_on_upper_edge(
      d, &mpi->decomp->ranges[rank], &mpi->decomp->parent_range);
    num_touches += mpi->is_on_edge[0][d] + mpi->is_on_edge[1][d];
  }
  mpi->touches_any_edge = num_touches > 0 ? true : false;
  
  mpi->nrecv = 0;
  for (int i=0; i<MAX_RECV_NEIGH; ++i)
    mpi->recv[i].buff = gkyl_mem_buff_new(16);
  
  mpi->nsend = 0;
  for (int i=0; i<MAX_RECV_NEIGH; ++i)
    mpi->send[i].buff = gkyl_mem_buff_new(16);
  
  mpi->allgather_buff_local.buff = gkyl_mem_buff_new(16);
  mpi->allgather_buff_global.buff = gkyl_mem_buff_new(16);
  
  mpi->priv_comm.gkyl_array_sync = array_sync;
  mpi->priv_comm.gkyl_array_per_sync = array_per_sync;
  mpi->priv_comm.gkyl_array_write = array_write;
  mpi->priv_comm.gkyl_array_read = array_read;
  mpi->priv_comm.gkyl_array_allgather = array_allgather;
  mpi->priv_comm.gkyl_array_allgather_host = array_allgather;
  
  mpi->priv_comm.get_rank = get_rank;
  mpi->priv_comm.get_size = get_size;
  mpi->priv_comm.barrier = barrier;
  mpi->priv_comm.gkyl_array_bcast = array_bcast;
  mpi->priv_comm.gkyl_array_bcast_host = array_bcast;
  mpi->priv_comm.allreduce = allreduce;
  mpi->priv_comm.allreduce_host = allreduce;
  mpi->priv_comm.extend_comm = extend_comm;
  mpi->priv_comm.split_comm = split_comm;
  mpi->priv_comm.create_comm_from_ranks = create_comm_from_ranks;
  
  mpi->priv_comm.pub_comm.ref_count = gkyl_ref_count_init(comm_free);

  return &mpi->priv_comm.pub_comm;
}

struct gkyl_comm*
gkyl_mpi_comm_new(const struct gkyl_mpi_comm_inp *inp)
{
  return mpi_comm_new(inp, &(struct extra_mpi_comm_inp) {
      .is_comm_allocated = false
    }
  );
}

#endif
