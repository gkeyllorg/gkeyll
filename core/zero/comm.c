#include <gkyl_comm_priv.h>

struct gkyl_comm*
gkyl_comm_acquire(const struct gkyl_comm *comm)
{
  gkyl_ref_count_inc(&comm->ref_count);
  return (struct gkyl_comm*) comm;
}

void
gkyl_comm_release(const struct gkyl_comm *comm)
{
  if (comm)
    gkyl_ref_count_dec(&comm->ref_count);
}

int
gkyl_comm_get_rank(struct gkyl_comm *pcomm, int *rank)
{
  struct gkyl_comm_priv *comm = container_of(pcomm, struct gkyl_comm_priv, pub_comm);
  return comm->get_rank(pcomm, rank);
}

int
gkyl_comm_get_size(struct gkyl_comm *pcomm, int *sz)
{
  struct gkyl_comm_priv *comm = container_of(pcomm, struct gkyl_comm_priv, pub_comm);
  return comm->get_size(pcomm, sz);
}

int
gkyl_comm_allreduce(struct gkyl_comm *pcomm, enum gkyl_elem_type type,
  enum gkyl_array_op op, int nelem, const void *inp, void *out)
{
  struct gkyl_comm_priv *comm = container_of(pcomm, struct gkyl_comm_priv, pub_comm);  
  return comm->allreduce(pcomm, type, op, nelem, inp, out);
}

int
gkyl_comm_allreduce_host(struct gkyl_comm *pcomm, enum gkyl_elem_type type,
  enum gkyl_array_op op, int nelem, const void *inp, void *out)
{
  struct gkyl_comm_priv *comm = container_of(pcomm, struct gkyl_comm_priv, pub_comm);  
  return comm->allreduce_host(pcomm, type, op, nelem, inp, out);
}

int
gkyl_comm_array_allgather(struct gkyl_comm *pcomm, 
  const struct gkyl_range *local, const struct gkyl_range *global,
  const struct gkyl_array *array_local, struct gkyl_array *array_global)
{
  struct gkyl_comm_priv *comm = container_of(pcomm, struct gkyl_comm_priv, pub_comm);  
  return comm->gkyl_array_allgather(pcomm, local, global, array_local, array_global);
}

int
gkyl_comm_array_allgather_host(struct gkyl_comm *pcomm, 
  const struct gkyl_range *local, const struct gkyl_range *global,
  const struct gkyl_array *array_local, struct gkyl_array *array_global)
{
  struct gkyl_comm_priv *comm = container_of(pcomm, struct gkyl_comm_priv, pub_comm);  
  return comm->gkyl_array_allgather_host(pcomm, local, global, array_local, array_global);
}

int
gkyl_comm_array_bcast(struct gkyl_comm *pcomm, 
  const struct gkyl_array *array_send, struct gkyl_array *array_recv, int root)
{
  struct gkyl_comm_priv *comm = container_of(pcomm, struct gkyl_comm_priv, pub_comm);  
  return comm->gkyl_array_bcast(pcomm, array_send, array_recv, root);
}

int
gkyl_comm_array_bcast_host(struct gkyl_comm *pcomm, 
  const struct gkyl_array *array_send, struct gkyl_array *array_recv, int root)
{
  struct gkyl_comm_priv *comm = container_of(pcomm, struct gkyl_comm_priv, pub_comm);
  return comm->gkyl_array_bcast_host(pcomm, array_send, array_recv, root);
}

int
gkyl_comm_array_sync(struct gkyl_comm *pcomm,
  const struct gkyl_range *local,
  const struct gkyl_range *local_ext,
  struct gkyl_array *array)
{
  struct gkyl_comm_priv *comm = container_of(pcomm, struct gkyl_comm_priv, pub_comm);  
  comm->barrier(pcomm);
  return comm->gkyl_array_sync(pcomm, local, local_ext, array);
}

int
gkyl_comm_array_per_sync(struct gkyl_comm *pcomm,
  const struct gkyl_range *local,
  const struct gkyl_range *local_ext,
  int nper_dirs, const int *per_dirs,
  struct gkyl_array *array)
{
  struct gkyl_comm_priv *comm = container_of(pcomm, struct gkyl_comm_priv, pub_comm);  
  comm->barrier(pcomm);
  return comm->gkyl_array_per_sync(pcomm, local, local_ext,
    nper_dirs, per_dirs, array);
}

int
gkyl_comm_barrier(struct gkyl_comm *pcomm)
{
  struct gkyl_comm_priv *comm = container_of(pcomm, struct gkyl_comm_priv, pub_comm);  
  return comm->barrier(pcomm);
}

void
gkyl_comm_group_call_start(struct gkyl_comm *pcomm)
{
  struct gkyl_comm_priv *comm = container_of(pcomm, struct gkyl_comm_priv, pub_comm);  
  comm->comm_group_call_start();
}

void
gkyl_comm_group_call_end(struct gkyl_comm *pcomm)
{
  struct gkyl_comm_priv *comm = container_of(pcomm, struct gkyl_comm_priv, pub_comm);  
  comm->comm_group_call_end();
}

int
gkyl_comm_array_write(struct gkyl_comm *pcomm,
  const struct gkyl_rect_grid *grid,
  const struct gkyl_range *range,
  const struct gkyl_msgpack_data *meta,
  const struct gkyl_array *arr, const char *fname)
{
  struct gkyl_comm_priv *comm = container_of(pcomm, struct gkyl_comm_priv, pub_comm);  
  int status = comm->gkyl_array_write(pcomm, grid, range, meta, arr, fname);
  gkyl_comm_barrier(pcomm);
  return status;
}

int
gkyl_comm_array_read(struct gkyl_comm *pcomm,
  const struct gkyl_rect_grid *grid, const struct gkyl_range *range,
  struct gkyl_array *arr, const char *fname)
{
  struct gkyl_comm_priv *comm = container_of(pcomm, struct gkyl_comm_priv, pub_comm);  
  int status = comm->gkyl_array_read(pcomm, grid, range, arr, fname);
  gkyl_comm_barrier(pcomm);
  return status;
}

struct gkyl_comm*
gkyl_comm_extend_comm(const struct gkyl_comm *pcomm,
  const struct gkyl_range *erange)
{
  struct gkyl_comm_priv *comm = container_of(pcomm, struct gkyl_comm_priv, pub_comm);
  return comm->extend_comm(pcomm, erange);
}

struct gkyl_comm*
gkyl_comm_split_comm(const struct gkyl_comm *pcomm, int color,
  struct gkyl_rect_decomp *new_decomp)
{
  struct gkyl_comm_priv *comm = container_of(pcomm, struct gkyl_comm_priv, pub_comm);  
  return comm->split_comm(pcomm, color, new_decomp);
}

struct gkyl_comm *
gkyl_comm_create_comm_from_ranks(const struct gkyl_comm *pcomm, int nranks,
  const int *ranks, struct gkyl_rect_decomp *new_decomp,
  bool *is_valid)
{
  struct gkyl_comm_priv *comm = container_of(pcomm, struct gkyl_comm_priv, pub_comm);  
  return comm->create_comm_from_ranks(pcomm, nranks, ranks, new_decomp, is_valid);
}
