#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_fem_gyroaverage.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_superlu_ops.h>

// Maximum number of perpendicular dimensions and of nodes in a cell.
#define FEM_GYROAVG_MAX_DIM 2
#define FEM_GYROAVG_MAX_NODES 4

// Maximum arc length of a gyroangle quadrature interval.
#define FEM_GYROAVG_MAX_ARC (M_PI / 4.0)

// Sparse matrix in compressed sparse row format.
struct fem_gyroavg_csr {
  long num_rows; // Number of rows.
  long num_nonzeros; // Number of nonzero entries.
  long *row_ptr; // Index of the first entry of each row (num_rows+1 values).
  long *col_idx; // Column of each entry.
  double *vals; // Value of each entry.
};

// Batch of updaters on the device, used to apply the updater of each of num_fld fields on the
// GPU. The matrices of all the fields and magnetic moments are stacked so that each step of
// the operator is a single kernel launch. Nodal vectors are stacked likewise: the one of the
// s-th field and n-th magnetic moment starts at (s*num_mu+n)*num_nodes.
struct gkyl_fem_gyroaverage_batch {
  int ndim; // Number of perpendicular dimensions.
  int num_basis; // Number of basis functions (and of nodes in a cell).
  int num_mu; // Number of magnetic moments.
  int num_fld; // Number of fields.
  long num_nodes; // Number of nodes of a field.
  long num_cells; // Number of cells of a field.
  int num_nodes_dir[FEM_GYROAVG_MAX_DIM]; // Number of nodes in each direction.
  // Load matrix, load_i = int phi_i*f dX: row i is the node i, and column cell*num_basis+k is
  // the k-th DG coefficient of f in the cell-th cell (the same for all fields).
  struct fem_gyroavg_csr load_mat;
  // Gyroaverage matrices B stacked: row (s*num_mu+n)*num_nodes+i (s-th field, n-th magnetic
  // moment, node i), column s*num_nodes+j.
  struct fem_gyroavg_csr gyro_mat;
  // Transposes of the gyroaverage matrices stacked: row s*num_nodes+j, column
  // (s*num_mu+n)*num_nodes+i.
  struct fem_gyroavg_csr gyro_mat_trans;
  // Inverse (row-major) of the mass matrix of the 1D nodal basis in each direction. The mass
  // matrix of the plane is their tensor product.
  double *mass_inv[FEM_GYROAVG_MAX_DIM];
  long node_stride[FEM_GYROAVG_MAX_DIM]; // Stride of each direction in the node index.
  long *cell_nodes; // Global index of the nodes of each cell (num_basis per cell).
  long *cell_linidx; // Linear index of each cell into the array of a field.
  double *nod2mod; // Nodal to modal transformation (FEM_GYROAVG_MAX_NODES^2, row-major).
  // Buffers for the nodal vectors of all fields and magnetic moments.
  double *nodal_in, *nodal_buff, *nodal_acc, *nodal_mass_buff;
};

// Updater type.
struct gkyl_fem_gyroaverage {
  int ndim; // Number of perpendicular dimensions.
  int num_basis; // Number of basis functions (and of nodes in a cell).
  struct gkyl_rect_grid grid; // Grid of the perpendicular plane.
  struct gkyl_basis basis; // Basis on the perpendicular plane.
  struct gkyl_range local; // Range the operator acts on.
  bool is_periodic[FEM_GYROAVG_MAX_DIM]; // Whether each direction is periodic.
  bool reflect_bnd; // Whether the field is reflected about non-periodic boundaries.
  int num_cells[FEM_GYROAVG_MAX_DIM]; // Number of cells in each direction.
  int num_nodes_dir[FEM_GYROAVG_MAX_DIM]; // Number of nodes in each direction.
  long num_nodes; // Total number of nodes.
  int num_mu; // Number of magnetic moments.

  // Offset (0 or 1) in each direction of each node of a cell from its lower corner.
  int node_offset[FEM_GYROAVG_MAX_NODES][FEM_GYROAVG_MAX_DIM];
  // Modal coefficients (first index) of each nodal basis function (second index).
  double nod2mod[FEM_GYROAVG_MAX_NODES][FEM_GYROAVG_MAX_NODES];
  double load_fac; // Cell volume over 2^ndim.
  long *cell_nodes; // Global index of the nodes of each cell (num_basis per cell).

  struct fem_gyroavg_csr *gyro_mat; // Gyroaverage matrix B of each magnetic moment.
  struct gkyl_superlu_prob *mass_prob; // Mass matrix solver.

  // Buffers for nodal vectors.
  double *nodal_in, *nodal_buff, *nodal_acc;

  bool use_gpu; // Whether to run on the GPU.
  struct gkyl_fem_gyroaverage_batch *batch; // Batch of this updater alone (if use_gpu).
};

#ifdef GKYL_HAVE_CUDA
/**
 * Gyroaverage each field of a batch on the GPU. See gkyl_fem_gyroaverage_batch_advance for
 * documentation.
 */
void gkyl_fem_gyroaverage_batch_advance_cu(
  struct gkyl_fem_gyroaverage_batch *batch, const struct gkyl_array *fin, struct gkyl_array *fout
);

/**
 * Apply the adjoint of the gyroaverage to each field of a batch on the GPU. See
 * gkyl_fem_gyroaverage_batch_advance_adjoint for documentation.
 */
void gkyl_fem_gyroaverage_batch_advance_adjoint_cu(
  struct gkyl_fem_gyroaverage_batch *batch, const struct gkyl_array *fin, struct gkyl_array *fout
);
#endif
