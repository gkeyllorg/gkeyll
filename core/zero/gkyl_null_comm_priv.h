#pragma once

// Private header for mpi_comm. Do not include in user-facing header files.

#include <gkyl_null_comm.h>
#include <gkyl_alloc.h>

// ranges for use in BCs
struct skin_ghost_ranges {
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];

  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];

  long max_vol; // maximum vol of send/recv region
};

// define long -> skin_ghost_ranges ...
#define i_key long
#define i_val struct skin_ghost_ranges
#define i_tag l2sgr
#include <stc/cmap.h>
// ... done with map definition

// Private struct
struct null_comm {
  struct gkyl_comm_priv priv_comm; // base communicator
  struct gkyl_rect_decomp *decomp; // pre-computed decomposition

  bool use_gpu; // flag to use if this communicator is on GPUs
  bool sync_corners; // should we sync corners?
  
  struct gkyl_range grange; // range to "hash" ghost layout

  cmap_l2sgr l2sgr; // map from long -> skin_ghost_ranges
  cmap_l2sgr l2sgr_wc; // map from long -> skin_ghost_ranges with corners
  
  gkyl_mem_buff pbuff; // CUDA buffer for periodic BCs
};

