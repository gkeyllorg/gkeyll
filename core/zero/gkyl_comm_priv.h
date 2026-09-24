#pragma once

#include <gkyl_comm.h>

// Get local "rank"
typedef int (*get_rank_t)(struct gkyl_comm *comm, int *rank);

// Get number of ranks
typedef int (*get_size_t)(struct gkyl_comm *comm, int *sz);

// "Reduce" all elements of @a type in array @a data and store output in @a out
typedef int (*allreduce_t)(struct gkyl_comm *comm, enum gkyl_elem_type type,
  enum gkyl_array_op op, int nelem, const void *inp, void *out);

// Gather local arrays into global array on each process.
typedef int (*gkyl_array_allgather_t)(struct gkyl_comm *comm,
  const struct gkyl_range *local, const struct gkyl_range *global,
  const struct gkyl_array *array_local, struct gkyl_array *array_global);

// Broadcast array to other processes.
typedef int (*gkyl_array_bcast_t)(struct gkyl_comm *comm,
  const struct gkyl_array *array_send, struct gkyl_array *array_recv, int root);

// "Synchronize" @a array across the regions or blocks.
typedef int (*gkyl_array_sync_t)(struct gkyl_comm *comm,
  const struct gkyl_range *local, const struct gkyl_range *local_ext,
  struct gkyl_array *array);

// "Synchronize" @a array across the periodic directions
typedef int (*gkyl_array_per_sync_t)(struct gkyl_comm *comm,
  const struct gkyl_range *local,
  const struct gkyl_range *local_ext,
  int nper_dirs, const int *per_dirs,
  struct gkyl_array *array);

// Write array to specified file
typedef int (*gkyl_array_write_t)(struct gkyl_comm *comm,
  const struct gkyl_rect_grid *grid,
  const struct gkyl_range *range,
  const struct gkyl_msgpack_data *meta,                                
  const struct gkyl_array *arr, const char *fname);

// Read array from specified file
typedef int (*gkyl_array_read_t)(struct gkyl_comm *comm,
  const struct gkyl_rect_grid *grid, const struct gkyl_range *range,
  struct gkyl_array *arr, const char *fname);

// Create a new communicator that extends the communcator to work on a
// extended domain specified by erange
typedef struct gkyl_comm* (*extend_comm_t)(const struct gkyl_comm *comm,
  const struct gkyl_range *erange);

// Create a new communicator by splitting a comm, and choosing members
// of new communicator according to the color rank. It can be used with
// a new decomp object, or the same one used for the parent comm, depending
// of the use case.
typedef struct gkyl_comm* (*split_comm_t)(const struct gkyl_comm *comm,
  int color, struct gkyl_rect_decomp *new_decomp);

// Create a new communicator from the input comm that takes a list of
// ranks to include in it.
typedef struct gkyl_comm *(*create_comm_from_ranks_t)(
  const struct gkyl_comm *comm, int nranks, const int *ranks,
  struct gkyl_rect_decomp *new_decomp,
  bool *is_valid
);

// Barrier
typedef int (*barrier_t)(struct gkyl_comm *comm);

// Start and end a group call (e.g. in NCCL).
typedef void (*comm_group_call_start_t)();
typedef void (*comm_group_call_end_t)();

// Private structure: not available to public-facing API but to
// specific comm types
struct gkyl_comm_priv {
  struct gkyl_comm pub_comm; // public facing communicator

  // FOLLOWING DO NOT NEED A DECOMP
  
  get_rank_t get_rank; // get local rank function.
  get_size_t get_size; // get number of ranks.
  barrier_t barrier; // barrier.
  allreduce_t allreduce; // all reduce function
  allreduce_t allreduce_host; // all reduce using the host (MPI) communicator
  
  extend_comm_t extend_comm; // extend communcator
  split_comm_t split_comm;   // split communicator.
  create_comm_from_ranks_t  create_comm_from_ranks; // communictor from ranks

  comm_group_call_start_t comm_group_call_start; // start a group call
  comm_group_call_end_t comm_group_call_end; // end a group call

  // FOLLOWING NEED A DECOMP 

  gkyl_array_allgather_t gkyl_array_allgather; // gather local arrays to global array
  gkyl_array_allgather_t gkyl_array_allgather_host; // gather local arrays to global array on host

  gkyl_array_bcast_t gkyl_array_bcast; // broadcast array to other processes
  gkyl_array_bcast_t gkyl_array_bcast_host; // broadcast host side array to other processes

  gkyl_array_sync_t gkyl_array_sync; // sync array
  gkyl_array_per_sync_t gkyl_array_per_sync; // sync array in periodic dirs

  gkyl_array_write_t gkyl_array_write; // array output
  gkyl_array_read_t gkyl_array_read; // array input
};
