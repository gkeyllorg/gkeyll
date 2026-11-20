/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_fem_gyroaverage.h>
#include <gkyl_fem_gyroaverage_priv.h>
#include <gkyl_util.h>
}

static int
fem_gyroavg_num_blocks(long num_threads)
{
  // Number of blocks of GKYL_DEFAULT_NUM_THREADS threads needed for num_threads threads.
  return (int)((num_threads + GKYL_DEFAULT_NUM_THREADS - 1) / GKYL_DEFAULT_NUM_THREADS);
}

__global__ static void
fem_gyroavg_load_cu_kernel(
  int num_fld, int num_comp, long fld_size, long num_nodes, int num_basis, const long *row_ptr,
  const long *col_idx, const double *vals, const long *cell_linidx, const struct gkyl_array *fdg,
  double *load
)
{
  // Compute the load vector int phi_i*f dX, phi_i being the nodal basis, of each of the
  // num_comp DG components f of each of the num_fld fields in fdg. The s-th field is in the
  // cells of fdg starting at s*fld_size, and the coefficients of its n-th component start at
  // n*num_basis. Its load starts at (s*num_comp+n)*num_nodes.
  long num_tot = num_fld * num_comp * num_nodes;
  for (unsigned long tid = threadIdx.x + blockIdx.x * blockDim.x; tid < num_tot;
       tid += gridDim.x * blockDim.x) {
    long vec = tid / num_nodes, node = tid % num_nodes;
    long fld = vec / num_comp, comp = vec % num_comp;
    double row_sum = 0.0;
    for (long k = row_ptr[node]; k < row_ptr[node + 1]; k++) {
      const double *fdg_d = (const double *)gkyl_array_cfetch(
        fdg, fld * fld_size + cell_linidx[col_idx[k] / num_basis]
      );
      row_sum += vals[k] * fdg_d[comp * num_basis + col_idx[k] % num_basis];
    }
    load[tid] = row_sum;
  }
}

__global__ static void
fem_gyroavg_csr_mat_vec_cu_kernel(
  long num_rows, const long *row_ptr, const long *col_idx, const double *vals, const double *xin,
  double *yout
)
{
  // Compute yout = mat*xin for a matrix in compressed sparse row format.
  for (unsigned long row = threadIdx.x + blockIdx.x * blockDim.x; row < num_rows;
       row += gridDim.x * blockDim.x) {
    double row_sum = 0.0;
    for (long k = row_ptr[row]; k < row_ptr[row + 1]; k++) {
      row_sum += vals[k] * xin[col_idx[k]];
    }
    yout[row] = row_sum;
  }
}

__global__ static void
fem_gyroavg_mass_inv_cu_kernel(
  long num_vec, long num_nodes, int num_nodes_dir, long node_stride, const double *mass_inv,
  const double *xin, double *yout
)
{
  // Apply the inverse of the 1D mass matrix along one direction (with num_nodes_dir nodes and
  // stride node_stride in the node index) to each of the num_vec nodal vectors in xin.
  long num_tot = num_vec * num_nodes;
  for (unsigned long tid = threadIdx.x + blockIdx.x * blockDim.x; tid < num_tot;
       tid += gridDim.x * blockDim.x) {
    long node = tid % num_nodes;
    long node_dir = (node / node_stride) % num_nodes_dir;
    // First node of the line of nodes along this direction that contains this node.
    const double *xin_line = &xin[tid - node_dir * node_stride];
    const double *mass_inv_row = &mass_inv[node_dir * num_nodes_dir];
    double row_sum = 0.0;
    for (int j = 0; j < num_nodes_dir; j++) {
      row_sum += mass_inv_row[j] * xin_line[j * node_stride];
    }
    yout[tid] = row_sum;
  }
}

__global__ static void
fem_gyroavg_nodal_to_dg_cu_kernel(
  int num_fld, int num_comp, long fld_size, long num_cells, long num_nodes, int num_basis,
  const long *cell_nodes, const long *cell_linidx, const double *nod2mod, const double *nodal,
  struct gkyl_array *fdg
)
{
  // Compute the DG coefficients of each of the num_comp continuous components of each of the
  // num_fld fields, with nodal values in nodal (those of the n-th component of the s-th field
  // start at (s*num_comp+n)*num_nodes). The fields are stored in fdg as in
  // fem_gyroavg_load_cu_kernel.
  long num_tot = num_fld * num_comp * num_cells;
  for (unsigned long tid = threadIdx.x + blockIdx.x * blockDim.x; tid < num_tot;
       tid += gridDim.x * blockDim.x) {
    long vec = tid / num_cells, cell = tid % num_cells;
    long fld = vec / num_comp, comp = vec % num_comp;
    const long *nodes = &cell_nodes[cell * num_basis];
    const double *nodal_vec = &nodal[vec * num_nodes];
    double *fdg_d = (double *)gkyl_array_fetch(fdg, fld * fld_size + cell_linidx[cell]);
    for (int k = 0; k < num_basis; k++) {
      double coeff = 0.0;
      for (int a = 0; a < num_basis; a++) {
        coeff += nod2mod[k * FEM_GYROAVG_MAX_NODES + a] * nodal_vec[nodes[a]];
      }
      fdg_d[comp * num_basis + k] = coeff;
    }
  }
}

static void
fem_gyroavg_mass_solve_cu(
  const struct gkyl_fem_gyroaverage_batch *batch, long num_vec, const double *rhs, double *sol
)
{
  // Solve M*sol = rhs for num_vec right-hand sides, M being the mass matrix of the nodal basis,
  // by applying the inverse of the 1D mass matrix of each direction in turn. The arrays rhs and
  // sol must differ.
  int num_blocks = fem_gyroavg_num_blocks(num_vec * batch->num_nodes);
  const double *xin = rhs;
  for (int d = batch->ndim - 1; d >= 0; d--) {
    double *yout = d == 0 ? sol : batch->nodal_mass_buff;
    fem_gyroavg_mass_inv_cu_kernel<<<num_blocks, GKYL_DEFAULT_NUM_THREADS>>>(
      num_vec, batch->num_nodes, batch->num_nodes_dir[d], batch->node_stride[d], batch->mass_inv[d],
      xin, yout
    );
    xin = yout;
  }
}

void
gkyl_fem_gyroaverage_batch_advance_cu(
  struct gkyl_fem_gyroaverage_batch *batch, const struct gkyl_array *fin, struct gkyl_array *fout
)
{
  int num_fld = batch->num_fld, num_mu = batch->num_mu, num_basis = batch->num_basis;
  long num_nodes = batch->num_nodes, num_cells = batch->num_cells;
  long fld_size = fin->size / num_fld;
  long num_vec_mu = (long)num_fld * num_mu;

  // Nodal values of the projection of the inputs onto the continuous space.
  fem_gyroavg_load_cu_kernel<<<
    fem_gyroavg_num_blocks(num_fld * num_nodes), GKYL_DEFAULT_NUM_THREADS>>>(
    num_fld, 1, fld_size, num_nodes, num_basis, batch->load_mat.row_ptr, batch->load_mat.col_idx,
    batch->load_mat.vals, batch->cell_linidx, fin->on_dev, batch->nodal_buff
  );
  fem_gyroavg_mass_solve_cu(batch, num_fld, batch->nodal_buff, batch->nodal_in);

  // Gyroaverage at all the magnetic moments.
  fem_gyroavg_csr_mat_vec_cu_kernel<<<
    fem_gyroavg_num_blocks(num_vec_mu * num_nodes), GKYL_DEFAULT_NUM_THREADS>>>(
    batch->gyro_mat.num_rows, batch->gyro_mat.row_ptr, batch->gyro_mat.col_idx,
    batch->gyro_mat.vals, batch->nodal_in, batch->nodal_buff
  );
  fem_gyroavg_mass_solve_cu(batch, num_vec_mu, batch->nodal_buff, batch->nodal_acc);
  fem_gyroavg_nodal_to_dg_cu_kernel<<<
    fem_gyroavg_num_blocks(num_vec_mu * num_cells), GKYL_DEFAULT_NUM_THREADS>>>(
    num_fld, num_mu, fld_size, num_cells, num_nodes, num_basis, batch->cell_nodes,
    batch->cell_linidx, batch->nod2mod, batch->nodal_acc, fout->on_dev
  );
}

void
gkyl_fem_gyroaverage_batch_advance_adjoint_cu(
  struct gkyl_fem_gyroaverage_batch *batch, const struct gkyl_array *fin, struct gkyl_array *fout
)
{
  int num_fld = batch->num_fld, num_mu = batch->num_mu, num_basis = batch->num_basis;
  long num_nodes = batch->num_nodes, num_cells = batch->num_cells;
  long fld_size = fin->size / num_fld;
  long num_vec_mu = (long)num_fld * num_mu;

  // Nodal values of the projection of the inputs at each magnetic moment.
  fem_gyroavg_load_cu_kernel<<<
    fem_gyroavg_num_blocks(num_vec_mu * num_nodes), GKYL_DEFAULT_NUM_THREADS>>>(
    num_fld, num_mu, fld_size, num_nodes, num_basis, batch->load_mat.row_ptr,
    batch->load_mat.col_idx, batch->load_mat.vals, batch->cell_linidx, fin->on_dev,
    batch->nodal_buff
  );
  fem_gyroavg_mass_solve_cu(batch, num_vec_mu, batch->nodal_buff, batch->nodal_in);

  // Sum over the magnetic moments of the transposed gyroaverage matrices.
  fem_gyroavg_csr_mat_vec_cu_kernel<<<
    fem_gyroavg_num_blocks(num_fld * num_nodes), GKYL_DEFAULT_NUM_THREADS>>>(
    batch->gyro_mat_trans.num_rows, batch->gyro_mat_trans.row_ptr, batch->gyro_mat_trans.col_idx,
    batch->gyro_mat_trans.vals, batch->nodal_in, batch->nodal_acc
  );
  fem_gyroavg_mass_solve_cu(batch, num_fld, batch->nodal_acc, batch->nodal_buff);
  fem_gyroavg_nodal_to_dg_cu_kernel<<<
    fem_gyroavg_num_blocks(num_fld * num_cells), GKYL_DEFAULT_NUM_THREADS>>>(
    num_fld, 1, fld_size, num_cells, num_nodes, num_basis, batch->cell_nodes, batch->cell_linidx,
    batch->nod2mod, batch->nodal_buff, fout->on_dev
  );
}
