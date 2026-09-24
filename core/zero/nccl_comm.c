#ifdef GKYL_HAVE_NCCL

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_rio_format_desc.h>
#include <gkyl_array_rio_priv.h>
#include <gkyl_comm_priv.h>
#include <gkyl_elem_type_priv.h>
#include <gkyl_mpi_comm.h>
#include <gkyl_nccl_comm.h>
#include <gkyl_nccl_comm_priv.h>
#include <gkyl_comm_io.h>
#include <gkyl_util.h>

#include <assert.h>
#include <errno.h>
#include <string.h>

// Mapping of Gkeyll type to ncclDataType_t
static ncclDataType_t g2_nccl_datatype[] = {
  [GKYL_INT] = ncclInt,
  [GKYL_INT_64] = ncclInt64,
  [GKYL_FLOAT] = ncclFloat,
  [GKYL_DOUBLE] = ncclDouble,
};

// Mapping of Gkeyll ops to ncclRedOp_t.
static ncclRedOp_t g2_nccl_op[] = {
  [GKYL_MIN] = ncclMin,
  [GKYL_MAX] = ncclMax,
  [GKYL_SUM] = ncclSum,
};

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

struct extra_nccl_comm_inp {
  bool is_comm_allocated; // is MPI_Comm allocated?
};

// Internal method to create a new NCCL communicator
static struct gkyl_comm* nccl_comm_new(
  const struct gkyl_nccl_comm_inp *inp, const struct extra_nccl_comm_inp *extra_inp);

struct gkyl_comm_state {
  ncclComm_t *ncomm;
  int tag;
  cudaStream_t *custream;
  int peer;
};

static struct gkyl_comm_state *
comm_state_new(struct gkyl_comm *comm)
{
  struct gkyl_comm_state *state = gkyl_malloc(sizeof *state);

  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, priv_comm.pub_comm);
  state->ncomm = &nccl->ncomm;
  state->custream = &nccl->custream;

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
  checkCuda(cudaStreamSynchronize(*(state->custream)));
}

// Ensure completion of NCCL collectives that are expected to behave synchronously.
static inline void
wait_for_nccl_collective(struct nccl_comm *nccl, ncclResult_t status)
{
  checkNCCL(status);
  checkCuda(cudaStreamSynchronize(nccl->custream));
}

static void
comm_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_comm *comm = container_of(ref, struct gkyl_comm, ref_count);
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, priv_comm.pub_comm);

  int ndim = nccl->decomp->ndim;
  gkyl_rect_decomp_release(nccl->decomp);

  gkyl_rect_decomp_neigh_release(nccl->neigh);
  for (int d=0; d<ndim; ++d)
    gkyl_rect_decomp_neigh_release(nccl->per_neigh[d]);

  for (int i=0; i<MAX_RECV_NEIGH; ++i)
    gkyl_mem_buff_release(nccl->recv[i].buff);

  for (int i=0; i<MAX_RECV_NEIGH; ++i)
    gkyl_mem_buff_release(nccl->send[i].buff);

  gkyl_mem_buff_release(nccl->allgather_buff_local.buff);
  gkyl_mem_buff_release(nccl->allgather_buff_global.buff);

  // Finalize NCCL comm.
  checkCuda(cudaStreamSynchronize(nccl->custream));
  checkCuda(cudaDeviceSynchronize());
  ncclCommDestroy(nccl->ncomm);

  if (nccl->is_mcomm_allocated)
    MPI_Comm_free(&nccl->mcomm);

  gkyl_comm_release(nccl->mpi_comm);
  gkyl_free(nccl);
}

static int
get_rank(struct gkyl_comm *comm, int *rank)
{
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, priv_comm.pub_comm);
  *rank = nccl->rank;
  return 0;
}

static int
get_size(struct gkyl_comm *comm, int *sz)
{
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, priv_comm.pub_comm);
  *sz = nccl->size;
  return 0;
}

static int
barrier(struct gkyl_comm *comm)
{
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, priv_comm.pub_comm);
  checkCuda(cudaStreamSynchronize(nccl->custream));
  MPI_Barrier(nccl->mcomm);
  return 0;
}

static int
array_write(struct gkyl_comm *comm,
  const struct gkyl_rect_grid *grid,
  const struct gkyl_range *range,
  const struct gkyl_msgpack_data *meta,
  const struct gkyl_array *arr, const char *fname)
{
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, priv_comm.pub_comm);
  return gkyl_comm_array_write(nccl->mpi_comm, grid, range, meta, arr, fname);
}

static int
array_read(struct gkyl_comm *comm,
  const struct gkyl_rect_grid *grid, const struct gkyl_range *range,
  struct gkyl_array *arr, const char *fname)
{
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, priv_comm.pub_comm);
  return gkyl_comm_array_read(nccl->mpi_comm, grid, range, arr, fname);
}

static int
array_send(struct gkyl_array *array, int dest, int tag, struct gkyl_comm *comm)
{
  size_t vol = array->ncomp*array->size;
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, priv_comm.pub_comm);
  checkNCCL(ncclSend(array->data, vol, g2_nccl_datatype[array->type], dest, nccl->ncomm, nccl->custream));
  checkCuda(cudaStreamSynchronize(nccl->custream));
  return 0;
}

static int
array_recv(struct gkyl_array *array, int src, int tag, struct gkyl_comm *comm)
{
  size_t vol = array->ncomp*array->size;
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, priv_comm.pub_comm);
  checkNCCL(ncclRecv(array->data, vol, g2_nccl_datatype[array->type], src, nccl->ncomm, nccl->custream));
  checkCuda(cudaStreamSynchronize(nccl->custream));
  return 0;
}

static int
array_isend(struct gkyl_array *array, int dest, int tag, struct gkyl_comm *comm, struct gkyl_comm_state *state)
{
  size_t vol = array->ncomp*array->size;
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, priv_comm.pub_comm);
  checkNCCL(ncclSend(array->data, vol, g2_nccl_datatype[array->type], dest, nccl->ncomm, nccl->custream));
  state->tag = tag;
  state->peer = dest;
  return 0;
}

static int
array_irecv(struct gkyl_array *array, int src, int tag, struct gkyl_comm *comm, struct gkyl_comm_state *state)
{
  size_t vol = array->ncomp*array->size;
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, priv_comm.pub_comm);
  checkNCCL(ncclRecv(array->data, vol, g2_nccl_datatype[array->type], src, nccl->ncomm, nccl->custream));
  state->tag = tag;
  state->peer = src;
  return 0;
}

static int
allreduce(struct gkyl_comm *comm, enum gkyl_elem_type type,
  enum gkyl_array_op op, int nelem, const void *inp,
  void *out)
{
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, priv_comm.pub_comm);
  ncclResult_t status = ncclAllReduce(inp, out, nelem,
    g2_nccl_datatype[type], g2_nccl_op[op], nccl->ncomm, nccl->custream);
  wait_for_nccl_collective(nccl, status);
  return 0;
}

static int
allreduce_host(struct gkyl_comm *comm, enum gkyl_elem_type type,
  enum gkyl_array_op op, int nelem, const void *inp,
  void *out)
{
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, priv_comm.pub_comm);  
  return gkyl_comm_allreduce(nccl->mpi_comm, type, op, nelem, inp, out);
}

static int
array_allgather(struct gkyl_comm *comm,
  const struct gkyl_range *local, const struct gkyl_range *global, 
  const struct gkyl_array *array_local, struct gkyl_array *array_global)
{
  assert(array_global->esznc == array_local->esznc);

  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, priv_comm.pub_comm);

  struct gkyl_range gather_range;

  int rank = nccl->rank;

  assert(local->volume == nccl->decomp->ranges[rank].volume);
  assert(global->volume == nccl->decomp->parent_range.volume);

  // potentially re-size local buffer volume
  size_t send_vol = array_local->esznc*nccl->decomp->ranges[rank].volume;
  if (gkyl_mem_buff_size(nccl->allgather_buff_local.buff) < send_vol)
    gkyl_mem_buff_resize(nccl->allgather_buff_local.buff, send_vol);

  // potentially re-size global buffer volume
  size_t buff_global_vol = array_local->esznc*nccl->decomp->parent_range.volume;
  if (gkyl_mem_buff_size(nccl->allgather_buff_global.buff) < buff_global_vol)
    gkyl_mem_buff_resize(nccl->allgather_buff_global.buff, buff_global_vol);

  // copy data to local buffer
  gkyl_array_copy_to_buffer(gkyl_mem_buff_data(nccl->allgather_buff_local.buff), 
    array_local, local);

  size_t nelem = array_local->esznc*nccl->decomp->ranges[rank].volume;
  // gather data into global buffer
  ncclResult_t status = ncclAllGather(
    gkyl_mem_buff_data(nccl->allgather_buff_local.buff),
    gkyl_mem_buff_data(nccl->allgather_buff_global.buff),
    nelem, ncclChar, nccl->ncomm, nccl->custream);
  wait_for_nccl_collective(nccl, status);

  // copy data to global array
  int idx = 0;
  for (int r=0; r<nccl->decomp->ndecomp; ++r) {
    int isrecv = gkyl_sub_range_intersect(
      &gather_range, global, &nccl->decomp->ranges[r]);
    gkyl_array_copy_from_buffer(array_global, 
      gkyl_mem_buff_data(nccl->allgather_buff_global.buff) + idx, &gather_range);
    idx += array_local->esznc*gather_range.volume;
  }

  return 0;
}

static int
array_allgather_host(struct gkyl_comm *comm,
  const struct gkyl_range *local, const struct gkyl_range *global, 
  const struct gkyl_array *array_local, struct gkyl_array *array_global)
{
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, priv_comm.pub_comm);
  return gkyl_comm_array_allgather_host(nccl->mpi_comm, local, global, array_local, array_global);
}

static int
array_bcast(struct gkyl_comm *comm, const struct gkyl_array *asend,
  struct gkyl_array *arecv, int root)
{
  assert(asend->esznc == arecv->esznc);
  assert(asend->size == arecv->size);

  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, priv_comm.pub_comm);

  size_t nelem = asend->ncomp*asend->size;

  ncclResult_t status = ncclBroadcast(asend->data, arecv->data, nelem,
    g2_nccl_datatype[asend->type], root, nccl->ncomm, nccl->custream);
  wait_for_nccl_collective(nccl, status);

  return 0;
}

static int
array_bcast_host(struct gkyl_comm *comm, const struct gkyl_array *asend,
  struct gkyl_array *arecv, int root)
{
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, priv_comm.pub_comm);
  return gkyl_comm_array_bcast_host(nccl->mpi_comm, asend, arecv, root);
}

static void
group_call_start()
{
  checkNCCL(ncclGroupStart());
}

static void
group_call_end()
{
  checkNCCL(ncclGroupEnd());
}

static int
array_sync(struct gkyl_comm *comm, const struct gkyl_range *local,
  const struct gkyl_range *local_ext, struct gkyl_array *array)
{
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, priv_comm.pub_comm);

  int elo[GKYL_MAX_DIM], eup[GKYL_MAX_DIM];
  for (int i=0; i<nccl->decomp->ndim; ++i)
    elo[i] = eup[i] = local_ext->upper[i]-local->upper[i];

  // Phase 1: Prepare recv ranges and buffers before the NCCL group.
  int recv_nids[MAX_RECV_NEIGH];
  int nridx = 0;
  for (int n=0; n<nccl->neigh->num_neigh; ++n) {
    int nid = nccl->neigh->neigh[n];

    int isrecv = gkyl_sub_range_intersect(
      &nccl->recv[nridx].range, local_ext, &nccl->decomp->ranges[nid]);
    size_t recv_vol = array->esznc*nccl->recv[nridx].range.volume;

    if (isrecv) {
      if (gkyl_mem_buff_size(nccl->recv[nridx].buff) < recv_vol)
        gkyl_mem_buff_resize(nccl->recv[nridx].buff, recv_vol);

      recv_nids[nridx] = nid;
      nridx += 1;
    }
  }

  // Phase 2: Prepare send ranges, resize buffers, and copy skin-cell data
  // into send buffers before the NCCL group.
  int send_nids[MAX_RECV_NEIGH];
  int nsidx = 0;
  for (int n=0; n<nccl->neigh->num_neigh; ++n) {
    int nid = nccl->neigh->neigh[n];

    struct gkyl_range neigh_ext;
    gkyl_range_extend(&neigh_ext, &nccl->decomp->ranges[nid], elo, eup);

    int issend = gkyl_sub_range_intersect(
      &nccl->send[nsidx].range, local, &neigh_ext);
    size_t send_vol = array->esznc*nccl->send[nsidx].range.volume;

    if (issend) {
      if (gkyl_mem_buff_size(nccl->send[nsidx].buff) < send_vol)
        gkyl_mem_buff_resize(nccl->send[nsidx].buff, send_vol);

      gkyl_array_copy_to_buffer(gkyl_mem_buff_data(nccl->send[nsidx].buff),
        array, &(nccl->send[nsidx].range));

      send_nids[nsidx] = nid;
      nsidx += 1;
    }
  }

  bool has_nccl_ops = (nridx > 0) || (nsidx > 0);
  if (has_nccl_ops) {
    // Phase 3: NCCL group with only ncclRecv/ncclSend calls inside.
    // No CUDA memory operations or kernel launches between group start/end.
    checkNCCL(ncclGroupStart());

    for (int r=0; r<nridx; ++r) {
      size_t recv_vol = array->esznc*nccl->recv[r].range.volume;
      checkNCCL(ncclRecv(gkyl_mem_buff_data(nccl->recv[r].buff),
        recv_vol, ncclChar, recv_nids[r], nccl->ncomm, nccl->custream));
    }

    for (int s=0; s<nsidx; ++s) {
      size_t send_vol = array->esznc*nccl->send[s].range.volume;
      checkNCCL(ncclSend(gkyl_mem_buff_data(nccl->send[s].buff),
        send_vol, ncclChar, send_nids[s], nccl->ncomm, nccl->custream));
    }

    checkNCCL(ncclGroupEnd());

    // Phase 4: Complete sends and recvs.
    checkCuda(cudaStreamSynchronize(nccl->custream));
  }

  // Phase 5: Copy received data into ghost-cells.
  for (int r=0; r<nridx; ++r) {
    gkyl_array_copy_from_buffer(array,
      gkyl_mem_buff_data(nccl->recv[r].buff),
      &(nccl->recv[r].range)
    );
  }

  return 0;
}

static int
array_per_sync(struct gkyl_comm *comm, const struct gkyl_range *local,
  const struct gkyl_range *local_ext,
  int nper_dirs, const int *per_dirs, struct gkyl_array *array)
{
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, priv_comm.pub_comm);

  if (!nccl->touches_any_edge) return 0; // nothing to sync

  int elo[GKYL_MAX_DIM], eup[GKYL_MAX_DIM];
  for (int i=0; i<nccl->decomp->ndim; ++i)
    elo[i] = eup[i] = local_ext->upper[i]-local->upper[i];

  int nridx = 0;
  int nsidx = 0;

  int shift_sign[] = { -1, 1 };

  // Handle self-periodic cases first (local copies, no NCCL communication).
  for (int i=0; i<nper_dirs; ++i) {
    int dir = per_dirs[i];

    for (int e=0; e<2; ++e) {
      if (nccl->is_on_edge[e][dir]) {
        int nid = nccl->per_neigh[dir]->neigh[0];
        if (nid == nccl->rank) {
          int delta[GKYL_MAX_DIM] = { 0 };
          delta[dir] = shift_sign[e]*gkyl_range_shape(&nccl->decomp->parent_range, dir);

          struct gkyl_range neigh_shift, neigh_shift_ext;
          gkyl_range_shift(&neigh_shift, &nccl->decomp->ranges[nid], delta);

          int isrecv = gkyl_sub_range_intersect(
            &nccl->recv[nridx].range, local_ext, &neigh_shift);

          delta[dir] *= -1;
          gkyl_range_shift(&neigh_shift, &nccl->decomp->ranges[nid], delta);
          gkyl_range_extend(&neigh_shift_ext, &neigh_shift, elo, eup);
          int issend = gkyl_sub_range_intersect(
            &nccl->send[nsidx].range, local, &neigh_shift_ext);

          size_t recv_vol = array->esznc*nccl->recv[nridx].range.volume;
          if (gkyl_mem_buff_size(nccl->recv[nridx].buff) < recv_vol)
            gkyl_mem_buff_resize(nccl->recv[nridx].buff, recv_vol);

          gkyl_array_copy_to_buffer(gkyl_mem_buff_data(nccl->recv[nridx].buff), array, &(nccl->send[nsidx].range));
          gkyl_array_copy_from_buffer(array, gkyl_mem_buff_data(nccl->recv[nridx].buff), &(nccl->recv[nridx].range));

          nridx += 1;
          nsidx += 1;
        }
      }
    }
  }

  // Phase 1: Prepare all multi-rank periodic recv/send ranges, resize buffers,
  // and copy skin-cell data into send buffers BEFORE the NCCL group.
  int nridx_nccl_start = nridx;
  int nsidx_nccl_start = nsidx;
  int recv_nids[MAX_RECV_NEIGH];
  int send_nids[MAX_RECV_NEIGH];

  for (int i=0; i<nper_dirs; ++i) {
    int dir = per_dirs[i];

    for (int e=0; e<2; ++e) {
      if (nccl->is_on_edge[e][dir]) {
        int nid = nccl->per_neigh[dir]->neigh[0];
        if (nid != nccl->rank) {
          int delta[GKYL_MAX_DIM] = { 0 };
          delta[dir] = shift_sign[e]*gkyl_range_shape(&nccl->decomp->parent_range, dir);

          if (nccl->per_neigh[dir]->num_neigh == 1) {

            struct gkyl_range neigh_shift;
            gkyl_range_shift(&neigh_shift, &nccl->decomp->ranges[nid], delta);

            int isrecv = gkyl_sub_range_intersect(
              &nccl->recv[nridx].range, local_ext, &neigh_shift);
            size_t recv_vol = array->esznc*nccl->recv[nridx].range.volume;

            if (isrecv) {
              if (gkyl_mem_buff_size(nccl->recv[nridx].buff) < recv_vol)
                gkyl_mem_buff_resize(nccl->recv[nridx].buff, recv_vol);

              recv_nids[nridx] = nid;
              nridx += 1;
            }

            struct gkyl_range neigh_shift_ext;
            gkyl_range_shift(&neigh_shift, &nccl->decomp->ranges[nid], delta);
            gkyl_range_extend(&neigh_shift_ext, &neigh_shift, elo, eup);

            int issend = gkyl_sub_range_intersect(
              &nccl->send[nsidx].range, local, &neigh_shift_ext);
            size_t send_vol = array->esznc*nccl->send[nsidx].range.volume;

            if (issend) {
              if (gkyl_mem_buff_size(nccl->send[nsidx].buff) < send_vol)
                gkyl_mem_buff_resize(nccl->send[nsidx].buff, send_vol);

              gkyl_array_copy_to_buffer(gkyl_mem_buff_data(nccl->send[nsidx].buff),
                array, &(nccl->send[nsidx].range));

              send_nids[nsidx] = nid;
              nsidx += 1;
            }
          }
        }
      }
    }
  }

  // Phase 2: NCCL group with only ncclRecv/ncclSend calls inside.
  // Only enter the group if this rank has actual operations to post.
  // Ranks that touch an edge in a non-periodic direction (e.g. velocity
  // in a phase-space decomposition) but not in any periodic direction
  // would otherwise post an empty ncclGroupStart/ncclGroupEnd, which
  // triggers a segfault in some NCCL versions (e.g. 2.25).
  bool has_nccl_ops = (nridx > nridx_nccl_start) || (nsidx > nsidx_nccl_start);

  if (has_nccl_ops) {
    checkNCCL(ncclGroupStart());

    for (int r=nridx_nccl_start; r<nridx; ++r) {
      size_t recv_vol = array->esznc*nccl->recv[r].range.volume;
      checkNCCL(ncclRecv(gkyl_mem_buff_data(nccl->recv[r].buff),
        recv_vol, ncclChar, recv_nids[r], nccl->ncomm, nccl->custream));
    }

    for (int s=nsidx_nccl_start; s<nsidx; ++s) {
      size_t send_vol = array->esznc*nccl->send[s].range.volume;
      checkNCCL(ncclSend(gkyl_mem_buff_data(nccl->send[s].buff),
        send_vol, ncclChar, send_nids[s], nccl->ncomm, nccl->custream));
    }

    checkNCCL(ncclGroupEnd());

    // Phase 3: Complete all sends and recvs.
    checkCuda(cudaStreamSynchronize(nccl->custream));
  }

  // Phase 4: Copy received data into ghost-cells.
  for (int r=nridx_nccl_start; r<nridx; ++r) {
    gkyl_array_copy_from_buffer(array,
      gkyl_mem_buff_data(nccl->recv[r].buff),
      &(nccl->recv[r].range)
    );
  }

  nccl->nrecv = nridx > nccl->nrecv ? nridx : nccl->nrecv;

  return 0;
}

static struct gkyl_comm*
extend_comm(const struct gkyl_comm *comm, const struct gkyl_range *erange)
{
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, priv_comm.pub_comm);

  // extend internal decomp object and create a new communicator
  struct gkyl_rect_decomp *ext_decomp = gkyl_rect_decomp_extended_new(erange, nccl->decomp);
  struct gkyl_comm *ext_comm = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
      .mpi_comm = nccl->mcomm,
      .decomp = ext_decomp,
      .sync_corners = nccl->sync_corners,
      .device_set = 1,
      .custream = nccl->custream,
    }
  );
  gkyl_rect_decomp_release(ext_decomp);
  return ext_comm;
}

static struct gkyl_comm*
split_comm(const struct gkyl_comm *comm, int color, struct gkyl_rect_decomp *new_decomp)
{
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, priv_comm.pub_comm);
  MPI_Comm new_mcomm;
  int ret = MPI_Comm_split(nccl->mcomm, color, nccl->rank, &new_mcomm);
  assert(ret == MPI_SUCCESS);

  struct gkyl_comm *newcomm = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
      .mpi_comm = new_mcomm,
      .decomp = new_decomp,
      .device_set = 1,
      .custream = nccl->custream,
    }
  );
  return newcomm;
}

static struct gkyl_comm*
create_comm_from_ranks(const struct gkyl_comm *comm,
  int nranks, const int *ranks, struct gkyl_rect_decomp *new_decomp,
  bool *is_valid)
{
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, priv_comm.pub_comm);

  MPI_Group group;
  MPI_Comm_group(nccl->mcomm, &group);

  MPI_Group new_group;
  MPI_Group_incl(group, nranks, ranks, &new_group);

  MPI_Comm new_mcomm;
  MPI_Comm_create_group(nccl->mcomm, new_group, 0, &new_mcomm);

  *is_valid = false;
  struct gkyl_comm *new_comm = 0;
  if (MPI_COMM_NULL != new_mcomm) {
    *is_valid = true;

    new_comm = nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
        .mpi_comm = new_mcomm,
        .sync_corners = nccl->sync_corners,
        .decomp = new_decomp,
        .device_set = 1,
        .custream = nccl->custream,
      },
      &(struct extra_nccl_comm_inp) {
        .is_comm_allocated = true
      }
    );
  }

  MPI_Group_free(&group);
  MPI_Group_free(&new_group);

  return new_comm;
}

struct gkyl_comm*
nccl_comm_new(const struct gkyl_nccl_comm_inp *inp,
  const struct extra_nccl_comm_inp *extra_inp)
{
  struct nccl_comm *nccl = gkyl_malloc(sizeof *nccl);
  strcpy(nccl->priv_comm.pub_comm.id, "nccl_comm");
  
  nccl->is_mcomm_allocated = extra_inp->is_comm_allocated;
  nccl->mcomm = inp->mpi_comm;

  nccl->mpi_comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
      .mpi_comm = nccl->mcomm,
      .decomp = inp->decomp, 
      .sync_corners = inp->sync_corners, 
    }
  );
  MPI_Comm_rank(nccl->mcomm, &nccl->rank);
  MPI_Comm_size(nccl->mcomm, &nccl->size);

  if (inp->device_set == 0) {
    int num_devices[1];
    checkCuda(cudaGetDeviceCount(num_devices));
  
    int local_rank = nccl->rank % num_devices[0];
    checkCuda(cudaSetDevice(local_rank));
  }

  ncclUniqueId nId;
  if (nccl->rank == 0) ncclGetUniqueId(&nId);
  MPI_Bcast((void *)&nId, sizeof(nId), MPI_BYTE, 0, nccl->mcomm);

  if (inp->custream == 0)
    checkCuda(cudaStreamCreate(&nccl->custream));
  else
    nccl->custream = inp->custream;

  // Initialize NCCL comm in blocking mode. Nonblocking only defers host-side
  // enqueue (forcing ncclCommGetAsyncError polling) and does NOT enable
  // compute/comm overlap; that requires separate CUDA streams plus events.
  ncclConfig_t config = NCCL_CONFIG_INITIALIZER;
  config.blocking = 1;
  checkNCCL(ncclCommInitRankConfig(&nccl->ncomm, nccl->size, nId, nccl->rank, &config));
  checkCuda(cudaStreamSynchronize(nccl->custream));

  nccl->sync_corners = inp->sync_corners;

  nccl->priv_comm.pub_comm.has_decomp = true;
  if (0 == inp->decomp) {
    nccl->priv_comm.pub_comm.has_decomp = false;

    // Construct a dummy decomposition.
    nccl->decomp = gkyl_rect_decomp_new_from_cuts_and_cells(1,
      (int[]) { nccl->size }, (int[]) { nccl->size });
  }
  else {
    nccl->decomp = gkyl_rect_decomp_acquire(inp->decomp);
  }

  nccl->neigh = gkyl_rect_decomp_calc_neigh(nccl->decomp, inp->sync_corners, nccl->rank);
  for (int d=0; d<nccl->decomp->ndim; ++d)
    nccl->per_neigh[d] =
      gkyl_rect_decomp_calc_periodic_neigh(nccl->decomp, d, false, nccl->rank);
  
  nccl->nrecv = 0;
  for (int i=0; i<MAX_RECV_NEIGH; ++i)
    nccl->recv[i].buff = gkyl_mem_buff_cu_new(16);
  
  nccl->nsend = 0;
  for (int i=0; i<MAX_RECV_NEIGH; ++i)
    nccl->send[i].buff = gkyl_mem_buff_cu_new(16);

  nccl->allgather_buff_local.buff = gkyl_mem_buff_cu_new(16);
  nccl->allgather_buff_global.buff = gkyl_mem_buff_cu_new(16);

  gkyl_range_init(&nccl->dir_edge, 2, (int[]) { 0, 0 }, (int[]) { GKYL_MAX_DIM, 2 });

  int num_touches = 0;
  for (int d=0; d<nccl->decomp->ndim; ++d) {
    nccl->is_on_edge[0][d] = gkyl_range_is_on_lower_edge(
      d, &nccl->decomp->ranges[nccl->rank], &nccl->decomp->parent_range);
    nccl->is_on_edge[1][d] = gkyl_range_is_on_upper_edge(
      d, &nccl->decomp->ranges[nccl->rank], &nccl->decomp->parent_range);
    num_touches += nccl->is_on_edge[0][d] + nccl->is_on_edge[1][d];
  }
  nccl->touches_any_edge = num_touches > 0 ? true : false;
  
  nccl->local_range_offset = gkyl_rect_decomp_calc_offset(nccl->decomp, nccl->rank);

  
  nccl->priv_comm.gkyl_array_sync = array_sync;
  nccl->priv_comm.gkyl_array_per_sync = array_per_sync;
  nccl->priv_comm.gkyl_array_write = array_write;
  nccl->priv_comm.gkyl_array_read = array_read;
  
  nccl->priv_comm.get_rank = get_rank;
  nccl->priv_comm.get_size = get_size;
  nccl->priv_comm.barrier = barrier;
  nccl->priv_comm.allreduce = allreduce;
  nccl->priv_comm.allreduce_host = allreduce_host;
// MF 2024/09/12: disable these for now per 498b7d1569eaa9285ae59581bd22dab124672f7b.
//  nccl->priv_comm.gkyl_array_send = array_send;
//  nccl->priv_comm.gkyl_array_isend = array_isend;
//  nccl->priv_comm.gkyl_array_recv = array_recv;
//  nccl->priv_comm.gkyl_array_irecv = array_irecv;
//  nccl->priv_comm.comm_state_new = comm_state_new;
//  nccl->priv_comm.comm_state_release = comm_state_release;
//  nccl->priv_comm.comm_state_wait = comm_state_wait;
  nccl->priv_comm.gkyl_array_allgather = array_allgather;
  nccl->priv_comm.gkyl_array_allgather_host = array_allgather_host;
  nccl->priv_comm.gkyl_array_bcast = array_bcast;
  nccl->priv_comm.gkyl_array_bcast_host = array_bcast_host;
  nccl->priv_comm.comm_group_call_start = group_call_start;
  nccl->priv_comm.comm_group_call_end = group_call_end;
  nccl->priv_comm.extend_comm = extend_comm;
  nccl->priv_comm.split_comm = split_comm;
  nccl->priv_comm.create_comm_from_ranks = create_comm_from_ranks;
  nccl->priv_comm.pub_comm.ref_count = gkyl_ref_count_init(comm_free);

  return &nccl->priv_comm.pub_comm;
}

struct gkyl_comm*
gkyl_nccl_comm_new(const struct gkyl_nccl_comm_inp *inp)
{
  return nccl_comm_new(inp, &(struct extra_nccl_comm_inp) {
      .is_comm_allocated = false
    }
  );
}

#endif
