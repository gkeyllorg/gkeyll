/* -*- c++ -*- */

#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_vlasov_conf_flux_surf.h>
#include <gkyl_dg_vlasov_conf_flux_surf_priv.h>
#include <gkyl_util.h> 
}

#include <cassert>

static void
gkyl_parallelize_components_kernel_launch_dims(dim3* dimGrid, dim3* dimBlock, gkyl_range range, int ncomp)
{
  // Create a 2D thread grid so we launch ncomp*range.volume number of threads
  // so we can parallelize over components too
  dimBlock->y = ncomp; // ncomp *must* be less than 256
  dimGrid->y = 1;
  dimBlock->x = GKYL_DEFAULT_NUM_THREADS/ncomp;
  dimGrid->x = gkyl_int_div_up(range.volume, dimBlock->x);
}

__global__ void
gkyl_dg_vlasov_conf_flux_surf_advance_cu_kernel(struct gkyl_dg_vlasov_conf_flux_surf *up,
  struct gkyl_range conf_range, struct gkyl_range phase_range, struct gkyl_range phase_range_ext,
  const struct gkyl_array *vmap, const struct gkyl_array *jacob_pos, const struct gkyl_array *jacob_vel_surf,
  const struct gkyl_array *poisson_tensor_conf, const struct gkyl_array *hamil,
  const struct gkyl_array *fin, struct gkyl_array *cflrate, struct gkyl_array *conf_flux_surf)
{
  // Per-node |alpha| values for the block, reduced to alpha_max per (cell, dir)
  // by the threadIdx.y == 0 thread of each cell.
  __shared__ double alpha_smem[GKYL_DEFAULT_NUM_THREADS];

  int pdim = up->pdim;
  int cdim = up->cdim;

  // 2D thread grid: linc2 indexes the surface node (i_node: transverse
  // configuration-space nodes, m_node: velocity-space nodes, matching the CPU
  // dispatch's i-major/m-minor node order); threads in x index the phase-space
  // cell. blockDim.y == num_nodes_conf*num_nodes_vel by construction.
  int num_nodes_vel = up->num_nodes_vel;
  int num_nodes = up->num_nodes_conf*num_nodes_vel;
  int linc2 = threadIdx.y;
  int i_node = linc2 / num_nodes_vel;
  int m_node = linc2 % num_nodes_vel;

  // No grid-stride loop: the launch covers phase_range.volume exactly so the
  // __syncthreads() barriers below are reached uniformly by every thread in
  // the block (threads past the end of the range only skip the guarded work).
  unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
  bool valid = linc1 < phase_range.volume;

  int idx[GKYL_MAX_DIM], idx_l[GKYL_MAX_DIM], idx_r[GKYL_MAX_DIM], idx_vel[GKYL_MAX_DIM];
  int idx_hamil[GKYL_MAX_DIM];
  double xcC[GKYL_MAX_DIM];
  const double *poisson_tensor_conf_d = 0, *hamil_d = 0, *f_c = 0;
  const double *vmap_d = 0, *jacob_vel_surf_d = 0, *jacob_pos_c = 0;
  double *cflrate_d = 0, *flux = 0;
  if (valid) {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&phase_range, linc1, idx);
    long cidx = gkyl_range_idx(&conf_range, idx);
    long pidx = gkyl_range_idx(&phase_range, idx);

    for (int i=0; i<up->hamil_dim; ++i) {
      idx_hamil[i] = idx[up->hamil_offset+i];
    }
    long hidx = gkyl_range_idx(&up->hamil_range, idx_hamil);

    // Grab the cell center location for NC bracket calculation
    gkyl_rect_grid_cell_center(&up->phase_grid, idx, xcC);

    f_c = (const double*) gkyl_array_cfetch(fin, pidx);
    cflrate_d = (double*) gkyl_array_fetch(cflrate, pidx);
    poisson_tensor_conf_d = (const double*) gkyl_array_cfetch(poisson_tensor_conf, cidx);
    hamil_d = (const double*) gkyl_array_cfetch(hamil, hidx);
    flux = (double*) gkyl_array_fetch(conf_flux_surf, pidx);
    for (int i=0; i<pdim-cdim; ++i) {
      idx_vel[i] = idx[cdim+i];
    }
    long vidx = gkyl_range_idx(&up->vel_range, idx_vel);
    vmap_d = (const double*) gkyl_array_cfetch(vmap, vidx);
    jacob_vel_surf_d = (const double*) gkyl_array_cfetch(jacob_vel_surf, vidx);
    jacob_pos_c = (const double*) gkyl_array_cfetch(jacob_pos, cidx);
  }

  // Each cell owns *lower* fluxes in each configuration-space direction.
  // So we need the distribution function in our current cell, and the cell
  // one index lower in each direction. If we are at the lower configuration-space
  // edge, we call ghost cells
  for (int dir = 0; dir<cdim; ++dir) {
    // Lower face owned by this cell (the left neighbor may be a ghost cell).
    // hamil_pt_edge = -1: evaluate the Hamiltonian/PT on this cell's lower face.
    const double *f_l = 0;
    double node_alpha = 0.0;
    const double *jacob_pos_l = 0;
    if (valid) {
      gkyl_copy_int_arr(pdim, idx, idx_l);
      idx_l[dir] = idx_l[dir]-1;
      long pidx_l = gkyl_range_idx(&phase_range, idx_l);
      f_l = (const double*) gkyl_array_cfetch(fin, pidx_l);
      long cidx_l = gkyl_range_idx(&conf_range, idx_l);
      jacob_pos_l = (const double*) gkyl_array_cfetch(jacob_pos, cidx_l);
      double alpha = up->hamil_alpha_quad[dir](i_node, m_node, -1, xcC, up->phase_grid.dx,
        vmap_d, jacob_pos_c, jacob_vel_surf_d, poisson_tensor_conf_d, hamil_d);
      node_alpha = up->lax_flux_nodal[dir](i_node, m_node, jacob_pos_l, jacob_pos_c, alpha, f_l, f_c, flux);
    }
    alpha_smem[threadIdx.x + blockDim.x*threadIdx.y] = node_alpha;
    __syncthreads();

    if (valid && threadIdx.y == 0) {
      // Reduce alpha_max in the CPU dispatch's node order so the fmax chain,
      // and hence the CFL estimate, matches the CPU loop exactly.
      double alpha_max = 0.0;
      for (int n=0; n<num_nodes; ++n) {
        alpha_max = fmax(alpha_max, alpha_smem[threadIdx.x + blockDim.x*n]);
      }
      double cfl = up->lax_cfl[dir](up->phase_grid.dx, jacob_pos_l, jacob_pos_c, alpha_max);
      // Always compute the flux, but if we are below threshold, ignore the stable time step estimate.
      if (fabs(f_l[0]) < up->skip_cell_thresh &&
          fabs(f_c[0]) < up->skip_cell_thresh) {
        cfl = 0.0;
      }
      cflrate_d[0] += cfl;
    }
    // alpha_smem is reused by the boundary pass and the next direction.
    __syncthreads();

    // If at the right boundary compute flux owned by the point in the ghost cell
    bool at_upper = valid && (idx[dir] == phase_range.upper[dir]);
    const double *f_r = 0;
    double *flux_r = 0, *cflrate_d_r = 0;
    const double *jacob_pos_r = 0;
    node_alpha = 0.0;
    if (at_upper) {
      // Index the right cell (ghost cell)
      gkyl_copy_int_arr(pdim, idx, idx_r);
      idx_r[dir] = idx_r[dir]+1;
      long pidx_r = gkyl_range_idx(&phase_range_ext, idx_r);

      f_r = (const double*) gkyl_array_cfetch(fin, pidx_r);
      flux_r = (double*) gkyl_array_fetch(conf_flux_surf, pidx_r);
      cflrate_d_r = (double*) gkyl_array_fetch(cflrate, pidx_r);

      /* As a concequence of not having ghost cells for PT/Hamil, they are shifted here
        and evaluated in the kernels at the upper boundary +1. This is allowed by continuity of hamil/pt */
      double xcR[GKYL_MAX_DIM];
      gkyl_rect_grid_cell_center(&up->phase_grid, idx_r, xcR);
      // Ghost-owned flux: the current cell is the l side, the ghost the c side
      // (ghost jacob_pos = skin value by the extended-range convention).
      long cidx_r = gkyl_range_idx(&conf_range, idx_r);
      jacob_pos_r = (const double*) gkyl_array_cfetch(jacob_pos, cidx_r);
      double alpha = up->hamil_alpha_quad[dir](i_node, m_node, 1, xcR, up->phase_grid.dx,
        vmap_d, jacob_pos_r, jacob_vel_surf_d, poisson_tensor_conf_d, hamil_d);
      node_alpha = up->lax_flux_nodal[dir](i_node, m_node, jacob_pos_c, jacob_pos_r, alpha, f_c, f_r, flux_r);
    }
    alpha_smem[threadIdx.x + blockDim.x*threadIdx.y] = node_alpha;
    __syncthreads();

    if (at_upper && threadIdx.y == 0) {
      double alpha_max = 0.0;
      for (int n=0; n<num_nodes; ++n) {
        alpha_max = fmax(alpha_max, alpha_smem[threadIdx.x + blockDim.x*n]);
      }
      double cfl = up->lax_cfl[dir](up->phase_grid.dx, jacob_pos_c, jacob_pos_r, alpha_max);
      if (fabs(f_c[0]) < up->skip_cell_thresh &&
          fabs(f_r[0]) < up->skip_cell_thresh) {
        cfl = 0.0;
      }
      cflrate_d_r[0] += cfl;
    }
    __syncthreads();
  }
}

void
gkyl_dg_vlasov_conf_flux_surf_advance_cu(struct gkyl_dg_vlasov_conf_flux_surf *up,
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, const struct gkyl_range *phase_range_ext,
  const struct gkyl_array *poisson_tensor_conf, const struct gkyl_array *hamil,
  const struct gkyl_array *fin, struct gkyl_array *cflrate, struct gkyl_array *conf_flux_surf)
{
  // 2D thread grid: parallelize over surface nodes as well as phase-space cells.
  int num_nodes = up->num_nodes_conf*up->num_nodes_vel;
  assert(num_nodes <= GKYL_DEFAULT_NUM_THREADS);
  dim3 dimGrid, dimBlock;
  gkyl_parallelize_components_kernel_launch_dims(&dimGrid, &dimBlock, *phase_range, num_nodes);
  gkyl_dg_vlasov_conf_flux_surf_advance_cu_kernel<<<dimGrid, dimBlock>>>(up->on_dev,
    *conf_range, *phase_range, *phase_range_ext,
    up->vmap->on_dev, up->jacob_pos->on_dev, up->jacob_vel_surf->on_dev,
    poisson_tensor_conf->on_dev,
    hamil->on_dev, fin->on_dev, cflrate->on_dev, conf_flux_surf->on_dev);
}

// CUDA kernel to set device pointers to canonical pb vars kernel functions
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
gkyl_dg_vlasov_conf_flux_surf_set_cu_dev_ptrs(struct gkyl_dg_vlasov_conf_flux_surf *up,
  enum gkyl_basis_type b_type, int cdim, int vdim, int poly_order, 
  enum gkyl_model_id model_id, enum gkyl_hamil_id hamil_id, bool use_lo)
{
  // Sparse (separable) vs. dense velocity-space Hamiltonian kernel selection.
  bool hamil_sparse = (hamil_id == GKYL_HAMIL_VEL_SPARSE);

  int kernel_index = cv_index[cdim].vdim[vdim];   
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      if ( use_lo ) {
        up->lax_flux_nodal[0] = ser_lax_flux_nodal_x_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[0] = ser_lax_flux_nodal_x_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[1] = ser_lax_flux_nodal_y_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[1] = ser_lax_flux_nodal_y_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[2] = ser_lax_flux_nodal_z_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[2] = ser_lax_flux_nodal_z_cfl_kernels[kernel_index].kernels[poly_order];
      }
      else {
        up->lax_flux_nodal[0] = ser_ho_lax_flux_nodal_x_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[0] = ser_ho_lax_flux_nodal_x_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[1] = ser_ho_lax_flux_nodal_y_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[1] = ser_ho_lax_flux_nodal_y_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[2] = ser_ho_lax_flux_nodal_z_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[2] = ser_ho_lax_flux_nodal_z_cfl_kernels[kernel_index].kernels[poly_order];
      }

      // Only have Hamiltonian forces in general geometry.
      if (model_id == GKYL_MODEL_TRIAD) {
        if ( use_lo ) {
          up->hamil_alpha_quad[0] = hamil_sparse ?
            ser_hamil_vel_sparse_alpha_quad_x_kernels[kernel_index].kernels[poly_order] :
            ser_hamil_vel_dense_alpha_quad_x_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[1] = hamil_sparse ?
            ser_hamil_vel_sparse_alpha_quad_y_kernels[kernel_index].kernels[poly_order] :
            ser_hamil_vel_dense_alpha_quad_y_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[2] = hamil_sparse ?
            ser_hamil_vel_sparse_alpha_quad_z_kernels[kernel_index].kernels[poly_order] :
            ser_hamil_vel_dense_alpha_quad_z_kernels[kernel_index].kernels[poly_order];
        } 
        else {
          up->hamil_alpha_quad[0] = hamil_sparse ?
            ser_hamil_vel_sparse_ho_alpha_quad_x_kernels[kernel_index].kernels[poly_order] :
            ser_hamil_vel_dense_ho_alpha_quad_x_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[1] = hamil_sparse ?
            ser_hamil_vel_sparse_ho_alpha_quad_y_kernels[kernel_index].kernels[poly_order] :
            ser_hamil_vel_dense_ho_alpha_quad_y_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[2] = hamil_sparse ?
            ser_hamil_vel_sparse_ho_alpha_quad_z_kernels[kernel_index].kernels[poly_order] :
            ser_hamil_vel_dense_ho_alpha_quad_z_kernels[kernel_index].kernels[poly_order];
        }
      }
      else if (model_id == GKYL_MODEL_TRIAD_GR
        || model_id == GKYL_MODEL_CANONICAL_PB || model_id == GKYL_MODEL_CANONICAL_PB_GR) {
        // Full phase-space Hamiltonian: alpha_dir = P . grad_v H evaluated at
        // the surface nodes. Canonical-PB models supply the identity Poisson
        // tensor, reducing this to the canonical streaming speed dH/dv_dir.
        if ( use_lo ) {
          up->hamil_alpha_quad[0] = ser_hamil_phase_alpha_quad_x_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[1] = ser_hamil_phase_alpha_quad_y_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[2] = ser_hamil_phase_alpha_quad_z_kernels[kernel_index].kernels[poly_order];
        } 
        else {
          up->hamil_alpha_quad[0] = ser_hamil_phase_ho_alpha_quad_x_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[1] = ser_hamil_phase_ho_alpha_quad_y_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[2] = ser_hamil_phase_ho_alpha_quad_z_kernels[kernel_index].kernels[poly_order];
        }
      }

      break;

    case GKYL_BASIS_MODAL_TENSOR:
      if (use_lo) {
        up->lax_flux_nodal[0] = tensor_lax_flux_nodal_x_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[0] = tensor_lax_flux_nodal_x_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[1] = tensor_lax_flux_nodal_y_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[1] = tensor_lax_flux_nodal_y_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[2] = tensor_lax_flux_nodal_z_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[2] = tensor_lax_flux_nodal_z_cfl_kernels[kernel_index].kernels[poly_order];
      } 
      else {
        up->lax_flux_nodal[0] = tensor_ho_lax_flux_nodal_x_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[0] = tensor_ho_lax_flux_nodal_x_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[1] = tensor_ho_lax_flux_nodal_y_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[1] = tensor_ho_lax_flux_nodal_y_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[2] = tensor_ho_lax_flux_nodal_z_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[2] = tensor_ho_lax_flux_nodal_z_cfl_kernels[kernel_index].kernels[poly_order];
      }
      
      if (model_id == GKYL_MODEL_TRIAD_GR
        || model_id == GKYL_MODEL_CANONICAL_PB || model_id == GKYL_MODEL_CANONICAL_PB_GR) {
        // Full phase-space Hamiltonian: alpha_dir = P . grad_v H evaluated at
        // the surface nodes. Canonical-PB models supply the identity Poisson
        // tensor, reducing this to the canonical streaming speed dH/dv_dir.
        // Only the p=1 tensor hybrid has a phase-space Hamiltonian representation.
        if ( use_lo ) {
          up->hamil_alpha_quad[0] = tensor_hamil_phase_alpha_quad_x_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[1] = tensor_hamil_phase_alpha_quad_y_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[2] = tensor_hamil_phase_alpha_quad_z_kernels[kernel_index].kernels[poly_order];
        }
        else {
          up->hamil_alpha_quad[0] = tensor_hamil_phase_ho_alpha_quad_x_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[1] = tensor_hamil_phase_ho_alpha_quad_y_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[2] = tensor_hamil_phase_ho_alpha_quad_z_kernels[kernel_index].kernels[poly_order];
        }
      }
      else if (model_id == GKYL_MODEL_TRIAD) {
        // Triad bracket on the tensor p=1 hybrid (sparse or dense velocity-space
        // Hamiltonian): per-node inverse velocity-map Jacobians of the C^1 cubic map.
        if ( use_lo ) {
          up->hamil_alpha_quad[0] = hamil_sparse ?
            tensor_hamil_vel_sparse_alpha_quad_x_kernels[kernel_index].kernels[poly_order] :
            tensor_hamil_vel_dense_alpha_quad_x_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[1] = hamil_sparse ?
            tensor_hamil_vel_sparse_alpha_quad_y_kernels[kernel_index].kernels[poly_order] :
            tensor_hamil_vel_dense_alpha_quad_y_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[2] = hamil_sparse ?
            tensor_hamil_vel_sparse_alpha_quad_z_kernels[kernel_index].kernels[poly_order] :
            tensor_hamil_vel_dense_alpha_quad_z_kernels[kernel_index].kernels[poly_order];
        }
        else {
          up->hamil_alpha_quad[0] = hamil_sparse ?
            tensor_hamil_vel_sparse_ho_alpha_quad_x_kernels[kernel_index].kernels[poly_order] :
            tensor_hamil_vel_dense_ho_alpha_quad_x_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[1] = hamil_sparse ?
            tensor_hamil_vel_sparse_ho_alpha_quad_y_kernels[kernel_index].kernels[poly_order] :
            tensor_hamil_vel_dense_ho_alpha_quad_y_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[2] = hamil_sparse ?
            tensor_hamil_vel_sparse_ho_alpha_quad_z_kernels[kernel_index].kernels[poly_order] :
            tensor_hamil_vel_dense_ho_alpha_quad_z_kernels[kernel_index].kernels[poly_order];
        }
      }

      break;

    default:
      assert(false);
      break;    
  } 
  // Set assembly functions for computing fluxes. 
  up->conf_flux_surf = conf_flux_surf_nodes;
  // Surface node counts for the per-node dispatch: (p+1) points per direction,
  // p+2 for the higher-order (anti-aliasing) and tensor (cubic-map) kernels.
  // The tensor p=1 hybrid (p=1 conf x p=2 vel) is anisotropic: 2 nodes per
  // configuration direction, 3 (lo) or 4 (ho) per velocity direction.
  int nq_conf = poly_order + 1, nq_vel = poly_order + 1;
  if ((poly_order > 1) && !use_lo) { nq_conf = poly_order + 2; nq_vel = poly_order + 2; }
  if (b_type == GKYL_BASIS_MODAL_TENSOR) {
    if (poly_order == 1) {
      nq_conf = 2;
      nq_vel = use_lo ? 3 : 4;
    }
    else {
      nq_conf = poly_order + 2;
      nq_vel = poly_order + 2;
    }
  }
  up->num_nodes_conf = 1;
  for (int d=0; d<cdim-1; ++d) up->num_nodes_conf *= nq_conf;
  up->num_nodes_vel = 1;
  for (int d=0; d<vdim; ++d) up->num_nodes_vel *= nq_vel;
}

gkyl_dg_vlasov_conf_flux_surf*
gkyl_dg_vlasov_conf_flux_surf_cu_dev_inew(const struct gkyl_dg_vlasov_conf_flux_surf_inp *inp)
{
  struct gkyl_dg_vlasov_conf_flux_surf *up = (struct gkyl_dg_vlasov_conf_flux_surf*) gkyl_malloc(sizeof(*up));

  int cdim = inp->conf_basis->ndim, pdim = inp->phase_basis->ndim, vdim = pdim-cdim;
  int poly_order = inp->conf_basis->poly_order;

  up->phase_grid = *inp->phase_grid;
  up->cdim = cdim;
  up->pdim = pdim;
  up->use_gpu = true; 

  // Are we skipping cells with small phase space density?
  if (inp->skip_cell_thresh > 0.0) {
    up->skip_cell_thresh = inp->skip_cell_thresh * pow(sqrt(2.0), pdim);
  }
  else {
    up->skip_cell_thresh = -1.0;
  }  

  // Determine Hamiltonian dimensionality and index offset for indexing Hamiltonian
  // from an input phase space index. 
  up->hamil_range = *inp->hamil_range; 
  if (inp->model_id == GKYL_MODEL_DEFAULT || inp->model_id == GKYL_MODEL_SR || inp->model_id == GKYL_MODEL_TRIAD) {
    up->hamil_dim = vdim; 
    up->hamil_offset = cdim; 
  }
  else {
    up->hamil_dim = pdim; 
    up->hamil_offset = 0; 
  }
  up->vel_range = *inp->vel_range;
  // Mesh maps: acquired for lifetime safety; raw device pointers unpacked at
  // launch time from the host-side borrowed array objects.
  assert(inp->vel_map);
  up->vel_map = gkyl_vlasov_velocity_map_acquire(inp->vel_map);
  up->vmap = inp->vel_map->vmap;
  up->jacob_vel_surf = inp->vel_map->jacob_vel_surf;
  assert(inp->pos_map);
  up->pos_map = gkyl_vlasov_position_map_acquire(inp->pos_map);
  up->jacob_pos = inp->pos_map->jacob_pos;

  // Host mirror of the surface node counts set on the device struct by
  // set_cu_dev_ptrs; the advance wrapper needs them to size the 2D
  // (cells x nodes) kernel launch.
  int nq_conf = poly_order + 1, nq_vel = poly_order + 1;
  if ((poly_order > 1) && !inp->use_lo) { nq_conf = poly_order + 2; nq_vel = poly_order + 2; }
  if (inp->conf_basis->b_type == GKYL_BASIS_MODAL_TENSOR) {
    if (poly_order == 1) {
      nq_conf = 2;
      nq_vel = inp->use_lo ? 3 : 4;
    }
    else {
      nq_conf = poly_order + 2;
      nq_vel = poly_order + 2;
    }
  }
  up->num_nodes_conf = 1;
  for (int d=0; d<cdim-1; ++d) up->num_nodes_conf *= nq_conf;
  up->num_nodes_vel = 1;
  for (int d=0; d<vdim; ++d) up->num_nodes_vel *= nq_vel;

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_dg_vlasov_conf_flux_surf *up_cu = (struct gkyl_dg_vlasov_conf_flux_surf*) gkyl_cu_malloc(sizeof(*up_cu));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_dg_vlasov_conf_flux_surf), GKYL_CU_MEMCPY_H2D);

  gkyl_dg_vlasov_conf_flux_surf_set_cu_dev_ptrs<<<1,1>>>(up_cu, inp->conf_basis->b_type, 
    cdim, vdim, poly_order, inp->model_id, inp->hamil_id, inp->use_lo);  

  // set parent on_dev pointer
  up->on_dev = up_cu;
 
  return up;
}
