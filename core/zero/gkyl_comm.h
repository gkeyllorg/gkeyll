#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_elem_type.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_ref_count.h>

// Structure holding data and function pointers to communicate various
// Gkeyll objects across multi-region or multi-block domains
struct gkyl_comm {
  char id[128]; // string ID for communcator.
  bool has_decomp; // flag to indicate if comm has an associated decomp.

  struct gkyl_ref_count ref_count; // reference count.
};

/**
 * Get rank of communicator.
 *
 * @param comm Communicator.
 * @param rank On output, the rank.
 * @return error code: 0 for success.
 */
int gkyl_comm_get_rank(struct gkyl_comm *comm, int *rank);

/**
 * Get number of ranks in communicator
 *
 * @param comm Communicator.
 * @param rank On output, the rank.
 * @return error code: 0 for success.
 */
int gkyl_comm_get_size(struct gkyl_comm *comm, int *sz);

/**
 * All reduce values across domains.
 *
 * @param comm Communicator.
 * @param type Data-type of element.
 * @param op Operator to use in reduction.
 * @param nelem Number of elemets in inp and out.
 * @param inp Local values on domain.
 * @param out Reduced values.
 * @return error code: 0 for success.
 */
int gkyl_comm_allreduce(struct gkyl_comm *comm, enum gkyl_elem_type type,
  enum gkyl_array_op op, int nelem, const void *inp, void *out);

/**
 * All reduce values across domains on the host/MPI communicator.
 *
 * @param comm Communicator.
 * @param type Data-type of element.
 * @param op Operator to use in reduction.
 * @param nelem Number of elemets in inp and out.
 * @param inp Local values on domain.
 * @param out Reduced values.
 * @return error code: 0 for success.
 */
int gkyl_comm_allreduce_host(struct gkyl_comm *comm, enum gkyl_elem_type type,
  enum gkyl_array_op op, int nelem, const void *inp, void *out);

/**
 * Gather all local data into a global array on each process.
 *
 * @param comm Communicator.
 * @param local Local range for array.
 * @param global Global range for array.
 * @param array_local Local array.
 * @param array_global Global array.
 * @return error code: 0 for success.
 */
int gkyl_comm_array_allgather(struct gkyl_comm *comm, 
  const struct gkyl_range *local, const struct gkyl_range *global,
  const struct gkyl_array *array_local, struct gkyl_array *array_global);

/**
 * Gather all local data on host into a global array on each process.
 *
 * @param comm Communicator.
 * @param local Local range for array.
 * @param global Global range for array.
 * @param array_local Local array.
 * @param array_global Global array.
 * @return error code: 0 for success.
 */
int gkyl_comm_array_allgather_host(struct gkyl_comm *comm, 
  const struct gkyl_range *local, const struct gkyl_range *global,
  const struct gkyl_array *array_local, struct gkyl_array *array_global);

/**
 * Broadcast an array to other processes.
 *
 * @param comm Communicator.
 * @param array_send Array to send (only in rank 'root').
 * @param array_recv Receive buffer array.
 * @param root Broadcasting process.
 * @return error code: 0 for success.
 */
int gkyl_comm_array_bcast(struct gkyl_comm *comm, 
  const struct gkyl_array *array_send, struct gkyl_array *array_recv, int root);

/**
 * Broadcast a host side array to other processes.
 *
 * @param comm Communicator.
 * @param array_send Array to send (only in rank 'root').
 * @param array_recv Receive buffer array.
 * @param root Broadcasting process.
 * @return error code: 0 for success.
 */
int gkyl_comm_array_bcast_host(struct gkyl_comm *comm, 
  const struct gkyl_array *array_send, struct gkyl_array *array_recv, int root);

/**
 * Synchronize array across domain.
 *
 * @param comm Communicator.
 * @param local Local range for array: sub-range of local_ext.
 * @param local_ext Extended range, i.e. range over which array is defined.
 * @param array Array to synchronize.
 * @return error code: 0 for success.
 */
int gkyl_comm_array_sync(struct gkyl_comm *comm,
  const struct gkyl_range *local,
  const struct gkyl_range *local_ext,
  struct gkyl_array *array);

/**
 * Synchronize array across domain in periodic directions.
 *
 * @param comm Communicator.
 * @param local Local range for array: sub-range of local_ext.
 * @param local_ext Extended range, i.e. range over which array is defined.
 * @param nper_dirs Number of periodic directions.
 * @param per_dirs Directions that are periodic.
 * @param array Array to synchronize.
 * @return error code: 0 for success.
 */
int gkyl_comm_array_per_sync(struct gkyl_comm *comm,
  const struct gkyl_range *local,
  const struct gkyl_range *local_ext,
  int nper_dirs, const int *per_dirs,
  struct gkyl_array *array);

/**
 * Barrier across domains
 *
 * @param comm Communicator.
 * @return error code: 0 for success.
 */
int gkyl_comm_barrier(struct gkyl_comm *comm);


/**
 * Start and end a group call
 * 
 * @param comm Communicator.
 */
void gkyl_comm_group_call_start(struct gkyl_comm *comm);
void gkyl_comm_group_call_end(struct gkyl_comm *comm);

/**
 * Create a new communcator that extends the communcator to work on a
 * extended domain specified by erange. (Each range handled by the
 * communicator is extended by a tensor-product with erange). The
 * returned communicator must be freed by calling gkyl_comm_release.
 *
 * @param comm Communicator.
 * @param erange Range to extend by.
 * @return Newly created communicator.
 */
struct gkyl_comm* gkyl_comm_extend_comm(const struct gkyl_comm *comm,
  const struct gkyl_range *erange);

/**
 * Split a communicator into a new communcator based on color. All
 * ranks with the same color will form the new communcator. In the input @a
 * new_decomp can be NULL.
 *
 * @param comm Communicator.
 * @param color All ranks of same color will share a communicator.
 * @param new_decomp Decomp object to associate new communicator. Can be NULL.
 * @return Newly created communicator.
 */
struct gkyl_comm* gkyl_comm_split_comm(const struct gkyl_comm *comm, int color,
  struct gkyl_rect_decomp *new_decomp);

/**
 * Create a new communicator that incudes a subset of ranks in @a
 * comm. This call can return a NULL if the communicator is not valid
 * on the parent calling rank. In this case the is_valid flag is also
 * set to false.
 *
 * @param comm Communicator.
 * @param nrank Number of ranks to include.
 * @param ranks List of ranks to include.
 * @param new_decomp Decomp object to associate new communicator. Can be NULL.
 * @param is_valid On output, true if comm is usable, false otherwise.
 * @return Newly created communicator.
 */
struct gkyl_comm* gkyl_comm_create_comm_from_ranks(const struct gkyl_comm *comm, int nranks,
  const int *ranks, struct gkyl_rect_decomp *new_decomp,
  bool *is_valid);

/**
 * Acquire pointer to communicator
 *
 * @param comm Communicator to to get acquire.
 * @return Acquired comm obj pointer.
 */
struct gkyl_comm* gkyl_comm_acquire(const struct gkyl_comm *comm);

/**
 * Release communicator memory.
 *
 * @param comm Communicator to release.
 */
void gkyl_comm_release(const struct gkyl_comm *comm);
