#include <gkyl_cuda_mpi_comm_priv.h>
#include <gkyl_multib_comm_conn_priv.h>

#if defined(GKYL_HAVE_MPI) && defined(GKYL_HAVE_CUDA)

#include <gkyl_array_ops.h>
#include <gkyl_comm_io.h>
#include <gkyl_mpi_comm.h>
#include <assert.h>
#include <limits.h>
#include <string.h>

static MPI_Datatype cuda_mpi_types[] = {
  [GKYL_INT] = MPI_INT,
  [GKYL_INT_64] = MPI_INT64_T,
  [GKYL_FLOAT] = MPI_FLOAT,
  [GKYL_DOUBLE] = MPI_DOUBLE,
  [GKYL_LONG] = MPI_LONG
};
static MPI_Op cuda_mpi_ops[] = {[GKYL_MIN] = MPI_MIN, [GKYL_MAX] = MPI_MAX, [GKYL_SUM] = MPI_SUM};

static struct cuda_mpi_comm *cuda_mpi_from_comm(const struct gkyl_comm *comm)
{
  return container_of(comm, struct cuda_mpi_comm, priv_comm.pub_comm);
}

static void cuda_mpi_check(struct cuda_mpi_comm *mpi, int status)
{
  if (status != MPI_SUCCESS) {
    char message[MPI_MAX_ERROR_STRING];
    int len;
    MPI_Error_string(status, message, &len);
    fprintf(stderr, "CUDA MPI rank %d: %s\n", mpi->rank, message);
    MPI_Abort(mpi->mcomm, status);
  }
}

static int cuda_mpi_count(struct cuda_mpi_comm *mpi, size_t count)
{
  // Do not silently truncate MPI-3 byte counts or strand peers on an error.
  if (count > INT_MAX) {
    cuda_mpi_check(mpi, MPI_ERR_COUNT);
  }
  return (int)count;
}

static void cuda_mpi_reserve(struct cuda_mpi_comm *mpi, int capacity)
{
  if (capacity <= mpi->capacity) {
    return;
  }
  mpi->send = gkyl_realloc(mpi->send, capacity * sizeof *mpi->send);
  mpi->recv = gkyl_realloc(mpi->recv, capacity * sizeof *mpi->recv);
  mpi->requests = gkyl_realloc(mpi->requests, 2 * capacity * sizeof *mpi->requests);
  for (int idx = mpi->capacity; idx < capacity; ++idx) {
    mpi->send[idx].buffer = gkyl_mem_buff_cu_new(16);
    mpi->recv[idx].buffer = gkyl_mem_buff_cu_new(16);
  }
  mpi->capacity = capacity;
}

static void cuda_mpi_message(
  struct cuda_mpi_comm *mpi, struct cuda_mpi_msg *msg, struct gkyl_array *array,
  const struct gkyl_range *range, int peer, int tag
)
{
  assert(gkyl_array_is_cu_dev(array));
  // Gkeyll's CUDA range-packing kernels operate on double components.
  if (array->type != GKYL_DOUBLE) {
    cuda_mpi_check(mpi, MPI_ERR_TYPE);
  }
  msg->array = array;
  msg->range = *range;
  msg->peer = peer;
  msg->tag = tag;
  msg->count = cuda_mpi_count(mpi, array->esznc * range->volume);
  if (gkyl_mem_buff_size(msg->buffer) < msg->count) {
    gkyl_mem_buff_resize(msg->buffer, msg->count);
  }
}

// All packing precedes MPI. MPI completion precedes unpacking. Reusable device
// buffers never pass through application-owned host staging memory. MPI may
// still stage internally, according to its selected transport and thresholds.
static int cuda_mpi_exchange(struct cuda_mpi_comm *mpi, int nsend, int nrecv)
{
  checkCuda(cudaDeviceSynchronize());
  for (int idx = 0; idx < nsend; ++idx) {
    struct cuda_mpi_msg *msg = &mpi->send[idx];
    if (msg->count > 0) {
      gkyl_array_copy_to_buffer(gkyl_mem_buff_data(msg->buffer), msg->array, &msg->range);
    }
  }
  checkCuda(cudaDeviceSynchronize());
  for (int idx = 0; idx < nrecv; ++idx) {
    struct cuda_mpi_msg *msg = &mpi->recv[idx];
    cuda_mpi_check(
      mpi, MPI_Irecv(
             gkyl_mem_buff_data(msg->buffer), msg->count, MPI_BYTE, msg->peer, msg->tag, mpi->mcomm,
             &mpi->requests[idx]
           )
    );
  }
  for (int idx = 0; idx < nsend; ++idx) {
    struct cuda_mpi_msg *msg = &mpi->send[idx];
    cuda_mpi_check(
      mpi, MPI_Isend(
             gkyl_mem_buff_data(msg->buffer), msg->count, MPI_BYTE, msg->peer, msg->tag, mpi->mcomm,
             &mpi->requests[nrecv + idx]
           )
    );
  }
  cuda_mpi_check(mpi, MPI_Waitall(nsend + nrecv, mpi->requests, MPI_STATUSES_IGNORE));
  for (int idx = 0; idx < nrecv; ++idx) {
    struct cuda_mpi_msg *msg = &mpi->recv[idx];
    if (msg->count > 0) {
      gkyl_array_copy_from_buffer(msg->array, gkyl_mem_buff_data(msg->buffer), &msg->range);
    }
  }
  checkCuda(cudaDeviceSynchronize());
  return 0;
}

static int cuda_mpi_get_rank(struct gkyl_comm *comm, int *rank)
{
  *rank = cuda_mpi_from_comm(comm)->rank;
  return 0;
}

static int cuda_mpi_get_size(struct gkyl_comm *comm, int *size)
{
  *size = cuda_mpi_from_comm(comm)->size;
  return 0;
}

static int cuda_mpi_barrier(struct gkyl_comm *comm)
{
  checkCuda(cudaDeviceSynchronize());
  return MPI_Barrier(cuda_mpi_from_comm(comm)->mcomm);
}

static void cuda_mpi_group(void)
{
}

static int cuda_mpi_allreduce_host(
  struct gkyl_comm *comm, enum gkyl_elem_type type, enum gkyl_array_op op, int count,
  const void *inp, void *out
)
{
  struct cuda_mpi_comm *mpi = cuda_mpi_from_comm(comm);
  if (type < GKYL_INT || type > GKYL_LONG) {
    return MPI_ERR_TYPE;
  }
  if (op < GKYL_MIN || op > GKYL_SUM) {
    return MPI_ERR_OP;
  }
  if (count < 0) {
    return MPI_ERR_COUNT;
  }
  return MPI_Allreduce(
    inp == out ? MPI_IN_PLACE : inp, out, count, cuda_mpi_types[type], cuda_mpi_ops[op], mpi->mcomm
  );
}

static int cuda_mpi_allreduce(
  struct gkyl_comm *comm, enum gkyl_elem_type type, enum gkyl_array_op op, int count,
  const void *inp, void *out
)
{
  checkCuda(cudaDeviceSynchronize());
  return cuda_mpi_allreduce_host(comm, type, op, count, inp, out);
}

static int cuda_mpi_bcast_host(
  struct gkyl_comm *comm, const struct gkyl_array *send, struct gkyl_array *recv, int root
)
{
  struct cuda_mpi_comm *mpi = cuda_mpi_from_comm(comm);
  // Only root needs a send array; the result is always placed in recv.
  if (mpi->rank == root && send != recv) {
    assert(send->esznc == recv->esznc && send->size == recv->size);
    gkyl_array_copy(recv, send);
  }
  return MPI_Bcast(
    recv->data, cuda_mpi_count(mpi, recv->esznc * recv->size), MPI_BYTE, root, mpi->mcomm
  );
}

static int cuda_mpi_bcast(
  struct gkyl_comm *comm, const struct gkyl_array *send, struct gkyl_array *recv, int root
)
{
  struct cuda_mpi_comm *mpi = cuda_mpi_from_comm(comm);
  checkCuda(cudaDeviceSynchronize());
  if (mpi->rank == root && send != recv) {
    assert(send->esznc == recv->esznc && send->size == recv->size);
    gkyl_array_copy(recv, send);
    checkCuda(cudaDeviceSynchronize());
  }
  return MPI_Bcast(
    recv->data, cuda_mpi_count(mpi, recv->esznc * recv->size), MPI_BYTE, root, mpi->mcomm
  );
}

static int cuda_mpi_gather(
  struct gkyl_comm *comm, const struct gkyl_range *local, const struct gkyl_range *global,
  const struct gkyl_array *src, struct gkyl_array *dst, bool device
)
{
  struct cuda_mpi_comm *mpi = cuda_mpi_from_comm(comm);
  assert(src->esznc == dst->esznc);
  if (device && (src->type != GKYL_DOUBLE || dst->type != GKYL_DOUBLE)) {
    cuda_mpi_check(mpi, MPI_ERR_TYPE);
  }
  assert(local->volume == mpi->decomp->ranges[mpi->rank].volume);
  assert(global->volume == mpi->decomp->parent_range.volume);
  size_t total = 0;
  for (int rank = 0; rank < mpi->size; ++rank) {
    mpi->offsets[rank] = cuda_mpi_count(mpi, total);
    mpi->counts[rank] = cuda_mpi_count(mpi, src->esznc * mpi->decomp->ranges[rank].volume);
    total += mpi->counts[rank];
  }
  cuda_mpi_count(mpi, total);
  gkyl_mem_buff local_buff = mpi->gather_local[device], global_buff = mpi->gather_global[device];
  if (gkyl_mem_buff_size(local_buff) < mpi->counts[mpi->rank]) {
    gkyl_mem_buff_resize(local_buff, mpi->counts[mpi->rank]);
  }
  if (gkyl_mem_buff_size(global_buff) < total) {
    gkyl_mem_buff_resize(global_buff, total);
  }
  if (device) {
    checkCuda(cudaDeviceSynchronize());
  }
  if (local->volume > 0) {
    gkyl_array_copy_to_buffer(gkyl_mem_buff_data(local_buff), src, local);
  }
  if (device) {
    checkCuda(cudaDeviceSynchronize());
  }
  int status = MPI_Allgatherv(
    gkyl_mem_buff_data(local_buff), mpi->counts[mpi->rank], MPI_BYTE,
    gkyl_mem_buff_data(global_buff), mpi->counts, mpi->offsets, MPI_BYTE, mpi->mcomm
  );
  if (status != MPI_SUCCESS) {
    return status;
  }
  for (int rank = 0; rank < mpi->size; ++rank) {
    struct gkyl_range range;
    gkyl_sub_range_intersect(&range, global, &mpi->decomp->ranges[rank]);
    if (range.volume > 0) {
      gkyl_array_copy_from_buffer(dst, gkyl_mem_buff_data(global_buff) + mpi->offsets[rank], &range);
    }
  }
  if (device) {
    checkCuda(cudaDeviceSynchronize());
  }
  return 0;
}

static int cuda_mpi_allgather(
  struct gkyl_comm *comm, const struct gkyl_range *local, const struct gkyl_range *global,
  const struct gkyl_array *src, struct gkyl_array *dst
)
{
  return cuda_mpi_gather(comm, local, global, src, dst, true);
}

static int cuda_mpi_allgather_host(
  struct gkyl_comm *comm, const struct gkyl_range *local, const struct gkyl_range *global,
  const struct gkyl_array *src, struct gkyl_array *dst
)
{
  return cuda_mpi_gather(comm, local, global, src, dst, false);
}

static int cuda_mpi_sync(
  struct gkyl_comm *comm, const struct gkyl_range *local, const struct gkyl_range *ext,
  struct gkyl_array *array
)
{
  struct cuda_mpi_comm *mpi = cuda_mpi_from_comm(comm);
  int ghost[GKYL_MAX_DIM] = {0};
  for (int dir = 0; dir < local->ndim; ++dir) {
    ghost[dir] = ext->upper[dir] - local->upper[dir];
  }
  int count = mpi->neigh->num_neigh;
  cuda_mpi_reserve(mpi, count);
  int passes = mpi->sync_corners ? local->ndim : 1;
  for (int pass = 0; pass < passes; ++pass) {
    for (int idx = 0; idx < count; ++idx) {
      int dir = mpi->neigh->dir[idx], edge = mpi->neigh->edge[idx];
      struct gkyl_range skin, halo;
      if (pass == 0) {
        gkyl_skin_ghost_ranges(&skin, &halo, dir, edge, ext, ghost);
      } else {
        gkyl_skin_ghost_with_corners_ranges(&skin, &halo, dir, edge, ext, ghost);
      }
      int peer = mpi->neigh->neigh[idx];
      cuda_mpi_message(mpi, &mpi->send[idx], array, &skin, peer, 100 + 2 * dir + edge);
      cuda_mpi_message(mpi, &mpi->recv[idx], array, &halo, peer, 100 + 2 * dir + (1 - edge));
    }
    cuda_mpi_exchange(mpi, count, count);
  }
  return 0;
}

static int cuda_mpi_per_sync(
  struct gkyl_comm *comm, const struct gkyl_range *local, const struct gkyl_range *ext, int ndirs,
  const int *dirs, struct gkyl_array *array
)
{
  struct cuda_mpi_comm *mpi = cuda_mpi_from_comm(comm);
  int ghost[GKYL_MAX_DIM] = {0};
  for (int dir = 0; dir < local->ndim; ++dir) {
    ghost[dir] = ext->upper[dir] - local->upper[dir];
  }
  cuda_mpi_reserve(mpi, 2 * local->ndim);
  int passes = mpi->sync_corners ? local->ndim : 1;
  for (int pass = 0; pass < passes; ++pass) {
    int count = 0;
    for (int idx = 0; idx < ndirs; ++idx) {
      int dir = dirs[idx];
      for (int edge = 0; edge < 2; ++edge) {
        bool on_edge = edge == 0 ? local->lower[dir] == mpi->decomp->parent_range.lower[dir] :
                                   local->upper[dir] == mpi->decomp->parent_range.upper[dir];
        if (!on_edge) {
          continue;
        }
        for (int nb = 0; nb < mpi->periodic[dir]->num_neigh; ++nb) {
          struct gkyl_range skin, halo;
          if (pass == 0) {
            gkyl_skin_ghost_ranges(&skin, &halo, dir, edge, ext, ghost);
          } else {
            gkyl_skin_ghost_with_corners_ranges(&skin, &halo, dir, edge, ext, ghost);
          }
          int peer = mpi->periodic[dir]->neigh[nb];
          cuda_mpi_reserve(mpi, count + 1);
          cuda_mpi_message(mpi, &mpi->send[count], array, &skin, peer, 200 + 2 * dir + edge);
          cuda_mpi_message(mpi, &mpi->recv[count], array, &halo, peer, 200 + 2 * dir + (1 - edge));
          ++count;
        }
      }
    }
    cuda_mpi_exchange(mpi, count, count);
  }
  return 0;
}

static int cuda_mpi_write(
  struct gkyl_comm *comm, const struct gkyl_rect_grid *grid, const struct gkyl_range *range,
  const struct gkyl_msgpack_data *meta, const struct gkyl_array *array, const char *fname
)
{
  assert(!gkyl_array_is_cu_dev(array));
  return gkyl_comm_array_write(cuda_mpi_from_comm(comm)->host, grid, range, meta, array, fname);
}

static int cuda_mpi_read(
  struct gkyl_comm *comm, const struct gkyl_rect_grid *grid, const struct gkyl_range *range,
  struct gkyl_array *array, const char *fname
)
{
  assert(!gkyl_array_is_cu_dev(array));
  return gkyl_comm_array_read(cuda_mpi_from_comm(comm)->host, grid, range, array, fname);
}

static struct gkyl_comm *
cuda_mpi_extend(const struct gkyl_comm *comm, const struct gkyl_range *range)
{
  struct cuda_mpi_comm *mpi = cuda_mpi_from_comm(comm);
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_extended_new(range, mpi->decomp);
  struct gkyl_comm *result = gkyl_cuda_mpi_comm_new(&(struct gkyl_cuda_mpi_comm_inp){
    .mpi_comm = mpi->mcomm, .decomp = decomp, .sync_corners = mpi->sync_corners, .device_set = true
  });
  gkyl_rect_decomp_release(decomp);
  return result;
}

static struct gkyl_comm *
cuda_mpi_split(const struct gkyl_comm *comm, int color, struct gkyl_rect_decomp *decomp)
{
  struct cuda_mpi_comm *mpi = cuda_mpi_from_comm(comm);
  MPI_Comm split;
  cuda_mpi_check(mpi, MPI_Comm_split(mpi->mcomm, color, mpi->rank, &split));
  if (split == MPI_COMM_NULL) {
    return 0;
  }
  struct gkyl_comm *result = gkyl_cuda_mpi_comm_new(&(
    struct gkyl_cuda_mpi_comm_inp
  ){.mpi_comm = split, .decomp = decomp, .sync_corners = mpi->sync_corners, .device_set = true});
  MPI_Comm_free(&split);
  return result;
}

static struct gkyl_comm *cuda_mpi_from_ranks(
  const struct gkyl_comm *comm, int nranks, const int *ranks, struct gkyl_rect_decomp *decomp,
  bool *valid
)
{
  struct cuda_mpi_comm *mpi = cuda_mpi_from_comm(comm);
  MPI_Group parent, group;
  MPI_Comm subset = MPI_COMM_NULL;
  cuda_mpi_check(mpi, MPI_Comm_group(mpi->mcomm, &parent));
  cuda_mpi_check(mpi, MPI_Group_incl(parent, nranks, ranks, &group));
  int group_rank;
  MPI_Group_rank(group, &group_rank);
  if (group_rank != MPI_UNDEFINED) {
    cuda_mpi_check(mpi, MPI_Comm_create_group(mpi->mcomm, group, 0, &subset));
  }
  MPI_Group_free(&group);
  MPI_Group_free(&parent);
  *valid = subset != MPI_COMM_NULL;
  if (!*valid) {
    return 0;
  }
  struct gkyl_comm *result = gkyl_cuda_mpi_comm_new(&(
    struct gkyl_cuda_mpi_comm_inp
  ){.mpi_comm = subset, .decomp = decomp, .sync_corners = mpi->sync_corners, .device_set = true});
  MPI_Comm_free(&subset);
  return result;
}

static void cuda_mpi_free(const struct gkyl_ref_count *ref)
{
  struct cuda_mpi_comm *mpi = cuda_mpi_from_comm(container_of(ref, struct gkyl_comm, ref_count));
  for (int idx = 0; idx < mpi->capacity; ++idx) {
    gkyl_mem_buff_release(mpi->send[idx].buffer);
    gkyl_mem_buff_release(mpi->recv[idx].buffer);
  }
  gkyl_free(mpi->send);
  gkyl_free(mpi->recv);
  gkyl_free(mpi->requests);
  gkyl_free(mpi->counts);
  gkyl_free(mpi->offsets);
  for (int idx = 0; idx < 2; ++idx) {
    gkyl_mem_buff_release(mpi->gather_local[idx]);
    gkyl_mem_buff_release(mpi->gather_global[idx]);
  }
  gkyl_rect_decomp_neigh_release(mpi->neigh);
  for (int dir = 0; dir < mpi->decomp->ndim; ++dir) {
    gkyl_rect_decomp_neigh_release(mpi->periodic[dir]);
  }
  gkyl_rect_decomp_release(mpi->decomp);
  gkyl_comm_release(mpi->host);
  MPI_Comm_free(&mpi->mcomm);
  gkyl_free(mpi);
}

struct gkyl_comm *gkyl_cuda_mpi_comm_new(const struct gkyl_cuda_mpi_comm_inp *inp)
{
  if (inp->mpi_comm == MPI_COMM_NULL) {
    return 0;
  }
  struct cuda_mpi_comm *mpi = gkyl_calloc(1, sizeof *mpi);
  mpi->mcomm = inp->mpi_comm;
  cuda_mpi_check(mpi, MPI_Comm_dup(inp->mpi_comm, &mpi->mcomm));
  MPI_Comm_rank(mpi->mcomm, &mpi->rank);
  MPI_Comm_size(mpi->mcomm, &mpi->size);
  if (!inp->device_set) {
    MPI_Comm shared;
    cuda_mpi_check(
      mpi, MPI_Comm_split_type(mpi->mcomm, MPI_COMM_TYPE_SHARED, mpi->rank, MPI_INFO_NULL, &shared)
    );
    int local_rank, devices;
    MPI_Comm_rank(shared, &local_rank);
    MPI_Comm_free(&shared);
    checkCuda(cudaGetDeviceCount(&devices));
    if (devices == 0) {
      cuda_mpi_check(mpi, MPI_ERR_OTHER);
    }
    checkCuda(cudaSetDevice(local_rank % devices));
  }
  mpi->decomp =
    inp->decomp ?
      gkyl_rect_decomp_acquire(inp->decomp) :
      gkyl_rect_decomp_new_from_cuts_and_cells(1, (int[]){mpi->size}, (int[]){mpi->size});
  if (mpi->decomp->ndecomp != mpi->size) {
    cuda_mpi_check(mpi, MPI_ERR_DIMS);
  }
  mpi->host = gkyl_mpi_comm_new(&(
    struct gkyl_mpi_comm_inp
  ){.mpi_comm = mpi->mcomm, .decomp = inp->decomp, .sync_corners = inp->sync_corners});
  mpi->sync_corners = inp->sync_corners;
  mpi->neigh = gkyl_rect_decomp_calc_neigh(mpi->decomp, false, mpi->rank);
  for (int dir = 0; dir < mpi->decomp->ndim; ++dir) {
    mpi->periodic[dir] = gkyl_rect_decomp_calc_periodic_neigh(mpi->decomp, dir, false, mpi->rank);
  }
  cuda_mpi_reserve(mpi, 2 * mpi->decomp->ndim);
  mpi->counts = gkyl_malloc(mpi->size * sizeof *mpi->counts);
  mpi->offsets = gkyl_malloc(mpi->size * sizeof *mpi->offsets);
  mpi->gather_local[0] = gkyl_mem_buff_new(16);
  mpi->gather_global[0] = gkyl_mem_buff_new(16);
  mpi->gather_local[1] = gkyl_mem_buff_cu_new(16);
  mpi->gather_global[1] = gkyl_mem_buff_cu_new(16);
  struct gkyl_comm_priv *priv = &mpi->priv_comm;
  strcpy(priv->pub_comm.id, "cuda_mpi_comm");
  priv->pub_comm.has_decomp = inp->decomp != 0;
  priv->pub_comm.ref_count = gkyl_ref_count_init(cuda_mpi_free);
  priv->get_rank = cuda_mpi_get_rank;
  priv->get_size = cuda_mpi_get_size;
  priv->barrier = cuda_mpi_barrier;
  priv->allreduce = cuda_mpi_allreduce;
  priv->allreduce_host = cuda_mpi_allreduce_host;
  priv->gkyl_array_bcast = cuda_mpi_bcast;
  priv->gkyl_array_bcast_host = cuda_mpi_bcast_host;
  priv->gkyl_array_allgather = cuda_mpi_allgather;
  priv->gkyl_array_allgather_host = cuda_mpi_allgather_host;
  priv->gkyl_array_sync = cuda_mpi_sync;
  priv->gkyl_array_per_sync = cuda_mpi_per_sync;
  priv->gkyl_array_write = cuda_mpi_write;
  priv->gkyl_array_read = cuda_mpi_read;
  priv->comm_group_call_start = cuda_mpi_group;
  priv->comm_group_call_end = cuda_mpi_group;
  priv->extend_comm = cuda_mpi_extend;
  priv->split_comm = cuda_mpi_split;
  priv->create_comm_from_ranks = cuda_mpi_from_ranks;
  return &priv->pub_comm;
}

int gkyl_multib_comm_conn_array_transfer_cuda_mpi(
  struct gkyl_comm *comm, int nblocks, const int *blocks, struct gkyl_multib_comm_conn **send,
  struct gkyl_multib_comm_conn **recv, struct gkyl_array **src, struct gkyl_array **dst
)
{
  struct cuda_mpi_comm *mpi = cuda_mpi_from_comm(comm);
  int nsend = 0, nrecv = 0;
  for (int block = 0; block < nblocks; ++block) {
    nsend += send[block]->num_comm_conn;
    nrecv += recv[block]->num_comm_conn;
  }
  cuda_mpi_reserve(mpi, nsend > nrecv ? nsend : nrecv);
  nsend = nrecv = 0;
  int *tag_max, flag;
  MPI_Comm_get_attr(mpi->mcomm, MPI_TAG_UB, &tag_max, &flag);
  for (int block = 0; block < nblocks; ++block) {
    for (int idx = 0; idx < send[block]->num_comm_conn; ++idx) {
      const struct gkyl_comm_conn *conn = &send[block]->comm_conn[idx];
      long tag = 4242L + 1000L * mpi->rank + 100L * conn->src_edge + blocks[block];
      if (tag > (flag ? *tag_max : 32767)) {
        cuda_mpi_check(mpi, MPI_ERR_TAG);
      }
      cuda_mpi_message(mpi, &mpi->send[nsend++], src[block], &conn->range, conn->rank, tag);
    }
    for (int idx = 0; idx < recv[block]->num_comm_conn; ++idx) {
      const struct gkyl_comm_conn *conn = &recv[block]->comm_conn[idx];
      long tag = 4242L + 1000L * conn->rank + 100L * conn->tar_edge + conn->block_id;
      if (tag > (flag ? *tag_max : 32767)) {
        cuda_mpi_check(mpi, MPI_ERR_TAG);
      }
      cuda_mpi_message(mpi, &mpi->recv[nrecv++], dst[block], &conn->range, conn->rank, tag);
    }
  }
  return cuda_mpi_exchange(mpi, nsend, nrecv);
}

#else

int gkyl_multib_comm_conn_array_transfer_cuda_mpi(
  struct gkyl_comm *comm, int nblocks, const int *blocks, struct gkyl_multib_comm_conn **send,
  struct gkyl_multib_comm_conn **recv, struct gkyl_array **src, struct gkyl_array **dst
)
{
  return 1;
}

#endif
