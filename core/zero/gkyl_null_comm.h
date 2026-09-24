#pragma once

#include <gkyl_comm.h>

#include <gkyl_rect_decomp.h>

// input to create new MPI communicator
struct gkyl_null_comm_inp {
  bool use_gpu; // flag to use if this communicator is on GPUs  
  const struct gkyl_rect_decomp *decomp; // pre-computed decomposition
  bool sync_corners; // should we sync corners?
};

/**
 * Return a new "null" communicator, i.e. a communicator for a single
 * core calculation.
 *
 * @param inp Input struct to use for initialization
 * @return New communicator
 */
struct gkyl_comm *gkyl_null_comm_inew(const struct gkyl_null_comm_inp *inp);

