#pragma once

#include <gkyl_orbit_average_gyrokinetic.h>

enum { GKYL_ORBIT_NUM_BASIS = 12, GKYL_ORBIT_NUM_QUAD = 32, GKYL_ORBIT_NUM_POS = 40 };

/** A coupled magnetic-moment cell; modal cells are contiguous in its list. */
struct gkyl_orbit_average_gyrokinetic_block {
  long cell_offset; // Offset in cell_indices.
  int num_cells;
  int num_modes; // Number of modal coefficients = num_cells*num_basis.
  int rank;
  long basis_offset; // Offset in q and dual; both have rank*num_modes entries.
  long reduce_offset; // Offset in reduced coefficients.
};

struct gkyl_orbit_average_gyrokinetic {
  bool use_gpu, positivity;
  int num_basis, num_pos, num_blocks;
  long num_cells, num_rejected, num_reduced;
  long basis_size;
  struct gkyl_orbit_average_gyrokinetic_block *blocks;
  long *cell_indices;
  double *q, *dual; // Column-major Q and M Q for M-orthonormal invariant columns.
  double *eval; // Row-major positivity-node phase basis; num_pos*num_basis.
  double *coeff, *correction, *blend; // Preallocated update work arrays.
  void *gpu; // CUDA-owned host descriptor with device buffers.
};

#ifdef GKYL_HAVE_CUDA
void gkyl_orbit_average_gyrokinetic_cu_init(gkyl_orbit_average_gyrokinetic *up);
void gkyl_orbit_average_gyrokinetic_cu_release(gkyl_orbit_average_gyrokinetic *up);
struct gkyl_orbit_average_gyrokinetic_result gkyl_orbit_average_gyrokinetic_advance_cu(
  gkyl_orbit_average_gyrokinetic *up, double dt, double frequency, struct gkyl_array *f
);
#endif
