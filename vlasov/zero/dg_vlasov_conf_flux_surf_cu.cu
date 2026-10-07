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
gkyl_parallelize_components_kernel_launch_dims(
  dim3 *dimGrid, dim3 *dimBlock, gkyl_range range, int ncomp
)
{
  // Create a 2D thread grid so we launch ncomp*range.volume number of threads
  // so we can parallelize over components too
  dimBlock->y = ncomp; // ncomp *must* be less than 256
  dimGrid->y = 1;
  dimBlock->x = GKYL_DEFAULT_NUM_THREADS / ncomp;
  dimGrid->x = gkyl_int_div_up(range.volume, dimBlock->x);
}

__global__ void
gkyl_dg_vlasov_conf_flux_surf_advance_cu_kernel(
  struct gkyl_dg_vlasov_conf_flux_surf *up, struct gkyl_range conf_range,
  struct gkyl_range phase_range, struct gkyl_range phase_range_ext, const struct gkyl_array *vmap,
  const struct gkyl_array *jacob_pos, const struct gkyl_array *jacob_vel_surf,
  const struct gkyl_array *poisson_tensor_conf, const struct gkyl_array *hamil,
  const struct gkyl_array *fin, struct gkyl_array *cflrate, struct gkyl_array *conf_flux_surf
)
{
  // Per-node |alpha| values for the block, reduced to alpha_max per (cell, face)
  // by the last node thread of each cell (the projection threads are the first
  // num_surf_basis node threads, so the two run in different warps and overlap).
  __shared__ double alpha_smem[GKYL_DEFAULT_NUM_THREADS];
  // Stage-1 arrays of the sum-factorized nodal f evaluation, G[m*NA + a] per
  // side, filled cooperatively (one thread per (inner node, outer shape) item;
  // there are at most as many items as surface nodes); each cell owns a block
  // of blockDim.y entries per side.
  __shared__ double G_smem[2 * GKYL_DEFAULT_NUM_THREADS];
  // This face's nodal Lax flux, one entry per surface node of each cell of the
  // block; projected onto the surface modal basis by the lax_prj stage (one
  // thread per surface mode) and only then written to the flux array.
  __shared__ double F_smem[GKYL_DEFAULT_NUM_THREADS];
  // Shared factors of the force producer (dynamic shared memory, sized by the
  // launcher): per cell, the outer factors O[t*NO + i] of every term t, then
  // the inner factors I[t*NI + m]; the force at a node is their dot product.
  extern __shared__ double alpha_factors_smem[];

  int pdim = up->pdim;
  int cdim = up->cdim;

  // 2D thread grid: linc2 indexes the surface node (i_node: transverse
  // configuration-space nodes, m_node: velocity-space nodes, matching the CPU
  // dispatch's i-major/m-minor node order); threads in x index the phase-space
  // cell. blockDim.y == num_nodes_conf*num_nodes_vel by construction.
  int num_nodes_vel = up->num_nodes_vel;
  int num_nodes = up->num_nodes_conf * num_nodes_vel;
  int linc2 = threadIdx.y;
  int i_node = linc2 / num_nodes_vel;
  int m_node = linc2 % num_nodes_vel;
  double *G_a = &G_smem[2 * blockDim.y * threadIdx.x];
  double *G_b = G_a + blockDim.y;
  double *F_cell = &F_smem[blockDim.y * threadIdx.x];
  const int NO = up->num_nodes_conf, NI = up->num_nodes_vel;
  double *O = &alpha_factors_smem[threadIdx.x * up->alpha_nterms_max * (NO + NI)];
  double *I = O + up->alpha_nterms_max * NO;

  // No grid-stride loop: the launch covers phase_range.volume exactly so the
  // __syncthreads() barriers below are reached uniformly by every thread in
  // the block (threads past the end of the range only skip the guarded work).
  unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x;
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

    for (int i = 0; i < up->hamil_dim; ++i) {
      idx_hamil[i] = idx[up->hamil_offset + i];
    }
    long hidx = gkyl_range_idx(&up->hamil_range, idx_hamil);

    // Grab the cell center location for NC bracket calculation
    gkyl_rect_grid_cell_center(&up->phase_grid, idx, xcC);

    f_c = (const double *)gkyl_array_cfetch(fin, pidx);
    cflrate_d = (double *)gkyl_array_fetch(cflrate, pidx);
    poisson_tensor_conf_d = (const double *)gkyl_array_cfetch(poisson_tensor_conf, cidx);
    hamil_d = (const double *)gkyl_array_cfetch(hamil, hidx);
    flux = (double *)gkyl_array_fetch(conf_flux_surf, pidx);
    for (int i = 0; i < pdim - cdim; ++i) {
      idx_vel[i] = idx[cdim + i];
    }
    long vidx = gkyl_range_idx(&up->vel_range, idx_vel);
    vmap_d = (const double *)gkyl_array_cfetch(vmap, vidx);
    jacob_vel_surf_d = (const double *)gkyl_array_cfetch(jacob_vel_surf, vidx);
    jacob_pos_c = (const double *)gkyl_array_cfetch(jacob_pos, cidx);
  }

  // Each cell owns *lower* fluxes in each configuration-space direction.
  // So we need the distribution function in our current cell, and the cell
  // one index lower in each direction. If we are at the lower configuration-space
  // edge, we call ghost cells
  for (int dir = 0; dir < cdim; ++dir) {
    // Lower face owned by this cell (the left neighbor may be a ghost cell).
    // hamil_pt_edge = -1: evaluate the Hamiltonian/PT on this cell's lower face.
    const double *f_l = 0;
    double node_alpha = 0.0;
    const double *jacob_pos_l = 0;
    int nt = 0;
    if (valid) {
      gkyl_copy_int_arr(pdim, idx, idx_l);
      idx_l[dir] = idx_l[dir] - 1;
      long pidx_l = gkyl_range_idx(&phase_range, idx_l);
      f_l = (const double *)gkyl_array_cfetch(fin, pidx_l);
      long cidx_l = gkyl_range_idx(&conf_range, idx_l);
      jacob_pos_l = (const double *)gkyl_array_cfetch(jacob_pos, cidx_l);
      // Shared work of the cell over its node threads: the force producer's
      // outer/inner factors on this cell's lower face (hamil_pt_edge = -1) and
      // stage 1 of the nodal f evaluation (l, c sides).
      nt = up->hamil_alpha_shared[dir](
        linc2, blockDim.y, xcC, up->phase_grid.dx, -1, vmap_d, jacob_pos_c, jacob_vel_surf_d,
        poisson_tensor_conf_d, hamil_d, O, I
      );
      up->lax_g[dir](linc2, f_l, f_c, G_a, G_b);
    }
    __syncthreads();
    if (valid) {
      // Per-node work: the force at this node is the dot product of the shared
      // factors over all terms.
      double alpha = 0.0;
      for (int t = 0; t < nt; ++t) {
        alpha += O[t * NO + i_node] * I[t * NI + m_node];
      }
      node_alpha =
        up->lax_flux_nodal[dir](i_node, m_node, jacob_pos_l, jacob_pos_c, alpha, G_a, G_b, F_cell);
    }
    alpha_smem[threadIdx.x + blockDim.x * threadIdx.y] = node_alpha;
    __syncthreads();

    if (valid) {
      // Stage 3: project this thread's surface mode of the cell's nodal flux
      // onto the surface modal basis (a no-op for linc2 >= num_surf_basis).
      up->lax_prj[dir](linc2, F_cell, flux);
    }
    if (valid && threadIdx.y == blockDim.y - 1) {
      // Reduce alpha_max in the CPU dispatch's node order so the fmax chain,
      // and hence the CFL estimate, matches the CPU loop exactly.
      double alpha_max = 0.0;
      for (int n = 0; n < num_nodes; ++n) {
        alpha_max = fmax(alpha_max, alpha_smem[threadIdx.x + blockDim.x * n]);
      }
      double cfl = up->lax_cfl[dir](up->phase_grid.dx, jacob_pos_l, jacob_pos_c, alpha_max);
      // Always compute the flux, but if we are below threshold, ignore the stable time step estimate.
      if (fabs(f_l[0]) < up->skip_cell_thresh && fabs(f_c[0]) < up->skip_cell_thresh) {
        cfl = 0.0;
      }
      cflrate_d[0] += cfl;
    }
    // alpha_smem, G_smem and F_smem are reused by the boundary pass and the next direction.
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
      idx_r[dir] = idx_r[dir] + 1;
      long pidx_r = gkyl_range_idx(&phase_range_ext, idx_r);

      f_r = (const double *)gkyl_array_cfetch(fin, pidx_r);
      flux_r = (double *)gkyl_array_fetch(conf_flux_surf, pidx_r);
      cflrate_d_r = (double *)gkyl_array_fetch(cflrate, pidx_r);

      /* As a concequence of not having ghost cells for PT/Hamil, they are shifted here
        and evaluated in the kernels at the upper boundary +1. This is allowed by continuity of hamil/pt */
      // Ghost-owned flux: the current cell is the l side, the ghost the c side
      // (ghost jacob_pos = skin value by the extended-range convention).
      long cidx_r = gkyl_range_idx(&conf_range, idx_r);
      jacob_pos_r = (const double *)gkyl_array_cfetch(jacob_pos, cidx_r);
      // Shared work for the ghost-owned upper face (hamil_pt_edge = +1): the
      // producer's factors and stage 1 of the nodal f evaluation (c, r sides).
      double xcR[GKYL_MAX_DIM];
      gkyl_rect_grid_cell_center(&up->phase_grid, idx_r, xcR);
      nt = up->hamil_alpha_shared[dir](
        linc2, blockDim.y, xcR, up->phase_grid.dx, 1, vmap_d, jacob_pos_r, jacob_vel_surf_d,
        poisson_tensor_conf_d, hamil_d, O, I
      );
      up->lax_g[dir](linc2, f_c, f_r, G_a, G_b);
    }
    __syncthreads();
    if (at_upper) {
      double alpha = 0.0;
      for (int t = 0; t < nt; ++t) {
        alpha += O[t * NO + i_node] * I[t * NI + m_node];
      }
      node_alpha =
        up->lax_flux_nodal[dir](i_node, m_node, jacob_pos_c, jacob_pos_r, alpha, G_a, G_b, F_cell);
    }
    alpha_smem[threadIdx.x + blockDim.x * threadIdx.y] = node_alpha;
    __syncthreads();

    if (at_upper) {
      up->lax_prj[dir](linc2, F_cell, flux_r);
    }
    if (at_upper && threadIdx.y == blockDim.y - 1) {
      double alpha_max = 0.0;
      for (int n = 0; n < num_nodes; ++n) {
        alpha_max = fmax(alpha_max, alpha_smem[threadIdx.x + blockDim.x * n]);
      }
      double cfl = up->lax_cfl[dir](up->phase_grid.dx, jacob_pos_c, jacob_pos_r, alpha_max);
      if (fabs(f_c[0]) < up->skip_cell_thresh && fabs(f_r[0]) < up->skip_cell_thresh) {
        cfl = 0.0;
      }
      cflrate_d_r[0] += cfl;
    }
    __syncthreads();
  }
}

void
gkyl_dg_vlasov_conf_flux_surf_advance_cu(
  struct gkyl_dg_vlasov_conf_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_range *phase_range_ext,
  const struct gkyl_array *poisson_tensor_conf, const struct gkyl_array *hamil,
  const struct gkyl_array *fin, struct gkyl_array *cflrate, struct gkyl_array *conf_flux_surf
)
{
  // 2D thread grid: parallelize over surface nodes as well as phase-space cells.
  int num_nodes = up->num_nodes_conf * up->num_nodes_vel;
  assert(num_nodes <= GKYL_DEFAULT_NUM_THREADS);
  dim3 dimGrid, dimBlock;
  gkyl_parallelize_components_kernel_launch_dims(&dimGrid, &dimBlock, *phase_range, num_nodes);
  // Dynamic shared memory for the force producer's shared factors.
  size_t alpha_smem =
    sizeof(double) * dimBlock.x * up->alpha_nterms_max * (up->num_nodes_conf + up->num_nodes_vel);
  gkyl_dg_vlasov_conf_flux_surf_advance_cu_kernel<<<dimGrid, dimBlock, alpha_smem>>>(
    up->on_dev, *conf_range, *phase_range, *phase_range_ext, up->vmap->on_dev,
    up->jacob_pos->on_dev, up->jacob_vel_surf->on_dev, poisson_tensor_conf->on_dev, hamil->on_dev,
    fin->on_dev, cflrate->on_dev, conf_flux_surf->on_dev
  );
}

// CUDA kernel to set device pointers to canonical pb vars kernel functions
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void
gkyl_dg_vlasov_conf_flux_surf_set_cu_dev_ptrs(
  struct gkyl_dg_vlasov_conf_flux_surf *up, enum gkyl_basis_type b_type, int cdim, int vdim,
  int poly_order, enum gkyl_model_id model_id, enum gkyl_hamil_id hamil_id, bool use_lo
)
{
  // Sparse (separable) vs. dense velocity-space Hamiltonian kernel selection.
  bool hamil_sparse = (hamil_id == GKYL_HAMIL_VEL_SPARSE);
  // By default, no configuration-space force (no terms).
  for (int d = 0; d < cdim; ++d) {
    up->hamil_alpha_shared[d] = no_hamil_alpha_shared;
  }

  int kernel_index = cv_index[cdim].vdim[vdim];
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      if (use_lo) {
        up->lax_g[0] = ser_lax_flux_nodal_x_g_kernels[kernel_index].kernels[poly_order];
        up->lax_prj[0] = ser_lax_flux_nodal_x_prj_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[0] = ser_lax_flux_nodal_x_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[0] = ser_lax_flux_nodal_x_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_g[1] = ser_lax_flux_nodal_y_g_kernels[kernel_index].kernels[poly_order];
        up->lax_prj[1] = ser_lax_flux_nodal_y_prj_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[1] = ser_lax_flux_nodal_y_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[1] = ser_lax_flux_nodal_y_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_g[2] = ser_lax_flux_nodal_z_g_kernels[kernel_index].kernels[poly_order];
        up->lax_prj[2] = ser_lax_flux_nodal_z_prj_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[2] = ser_lax_flux_nodal_z_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[2] = ser_lax_flux_nodal_z_cfl_kernels[kernel_index].kernels[poly_order];
      } else {
        up->lax_g[0] = ser_ho_lax_flux_nodal_x_g_kernels[kernel_index].kernels[poly_order];
        up->lax_prj[0] = ser_ho_lax_flux_nodal_x_prj_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[0] = ser_ho_lax_flux_nodal_x_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[0] = ser_ho_lax_flux_nodal_x_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_g[1] = ser_ho_lax_flux_nodal_y_g_kernels[kernel_index].kernels[poly_order];
        up->lax_prj[1] = ser_ho_lax_flux_nodal_y_prj_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[1] = ser_ho_lax_flux_nodal_y_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[1] = ser_ho_lax_flux_nodal_y_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_g[2] = ser_ho_lax_flux_nodal_z_g_kernels[kernel_index].kernels[poly_order];
        up->lax_prj[2] = ser_ho_lax_flux_nodal_z_prj_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[2] = ser_ho_lax_flux_nodal_z_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[2] = ser_ho_lax_flux_nodal_z_cfl_kernels[kernel_index].kernels[poly_order];
      }

      // Only have Hamiltonian forces in general geometry.
      if (model_id == GKYL_MODEL_TRIAD) {
        if (use_lo) {
          up->hamil_alpha_shared[0] =
            hamil_sparse ?
              ser_hamil_vel_sparse_alpha_quad_x_kernels[kernel_index].kernels[poly_order] :
              ser_hamil_vel_dense_alpha_quad_x_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[1] =
            hamil_sparse ?
              ser_hamil_vel_sparse_alpha_quad_y_kernels[kernel_index].kernels[poly_order] :
              ser_hamil_vel_dense_alpha_quad_y_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[2] =
            hamil_sparse ?
              ser_hamil_vel_sparse_alpha_quad_z_kernels[kernel_index].kernels[poly_order] :
              ser_hamil_vel_dense_alpha_quad_z_kernels[kernel_index].kernels[poly_order];
        } else {
          up->hamil_alpha_shared[0] =
            hamil_sparse ?
              ser_hamil_vel_sparse_ho_alpha_quad_x_kernels[kernel_index].kernels[poly_order] :
              ser_hamil_vel_dense_ho_alpha_quad_x_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[1] =
            hamil_sparse ?
              ser_hamil_vel_sparse_ho_alpha_quad_y_kernels[kernel_index].kernels[poly_order] :
              ser_hamil_vel_dense_ho_alpha_quad_y_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[2] =
            hamil_sparse ?
              ser_hamil_vel_sparse_ho_alpha_quad_z_kernels[kernel_index].kernels[poly_order] :
              ser_hamil_vel_dense_ho_alpha_quad_z_kernels[kernel_index].kernels[poly_order];
        }
      } else if (model_id == GKYL_MODEL_TRIAD_GR || model_id == GKYL_MODEL_CANONICAL_PB ||
                 model_id == GKYL_MODEL_CANONICAL_PB_GR) {
        // Full phase-space Hamiltonian: alpha_dir = P . grad_v H evaluated at
        // the surface nodes. Canonical-PB models supply the identity Poisson
        // tensor, reducing this to the canonical streaming speed dH/dv_dir.
        if (use_lo) {
          up->hamil_alpha_shared[0] =
            ser_hamil_phase_alpha_quad_x_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[1] =
            ser_hamil_phase_alpha_quad_y_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[2] =
            ser_hamil_phase_alpha_quad_z_kernels[kernel_index].kernels[poly_order];
        } else {
          up->hamil_alpha_shared[0] =
            ser_hamil_phase_ho_alpha_quad_x_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[1] =
            ser_hamil_phase_ho_alpha_quad_y_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[2] =
            ser_hamil_phase_ho_alpha_quad_z_kernels[kernel_index].kernels[poly_order];
        }
      }

      break;

    case GKYL_BASIS_MODAL_TENSOR:
      if (use_lo) {
        up->lax_g[0] = tensor_lax_flux_nodal_x_g_kernels[kernel_index].kernels[poly_order];
        up->lax_prj[0] = tensor_lax_flux_nodal_x_prj_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[0] = tensor_lax_flux_nodal_x_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[0] = tensor_lax_flux_nodal_x_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_g[1] = tensor_lax_flux_nodal_y_g_kernels[kernel_index].kernels[poly_order];
        up->lax_prj[1] = tensor_lax_flux_nodal_y_prj_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[1] = tensor_lax_flux_nodal_y_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[1] = tensor_lax_flux_nodal_y_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_g[2] = tensor_lax_flux_nodal_z_g_kernels[kernel_index].kernels[poly_order];
        up->lax_prj[2] = tensor_lax_flux_nodal_z_prj_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[2] = tensor_lax_flux_nodal_z_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[2] = tensor_lax_flux_nodal_z_cfl_kernels[kernel_index].kernels[poly_order];
      } else {
        up->lax_g[0] = tensor_ho_lax_flux_nodal_x_g_kernels[kernel_index].kernels[poly_order];
        up->lax_prj[0] = tensor_ho_lax_flux_nodal_x_prj_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[0] =
          tensor_ho_lax_flux_nodal_x_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[0] = tensor_ho_lax_flux_nodal_x_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_g[1] = tensor_ho_lax_flux_nodal_y_g_kernels[kernel_index].kernels[poly_order];
        up->lax_prj[1] = tensor_ho_lax_flux_nodal_y_prj_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[1] =
          tensor_ho_lax_flux_nodal_y_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[1] = tensor_ho_lax_flux_nodal_y_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_g[2] = tensor_ho_lax_flux_nodal_z_g_kernels[kernel_index].kernels[poly_order];
        up->lax_prj[2] = tensor_ho_lax_flux_nodal_z_prj_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[2] =
          tensor_ho_lax_flux_nodal_z_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[2] = tensor_ho_lax_flux_nodal_z_cfl_kernels[kernel_index].kernels[poly_order];
      }

      if (model_id == GKYL_MODEL_TRIAD_GR || model_id == GKYL_MODEL_CANONICAL_PB ||
          model_id == GKYL_MODEL_CANONICAL_PB_GR) {
        // Full phase-space Hamiltonian: alpha_dir = P . grad_v H evaluated at
        // the surface nodes. Canonical-PB models supply the identity Poisson
        // tensor, reducing this to the canonical streaming speed dH/dv_dir.
        // Only the p=1 tensor hybrid has a phase-space Hamiltonian representation.
        if (use_lo) {
          up->hamil_alpha_shared[0] =
            tensor_hamil_phase_alpha_quad_x_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[1] =
            tensor_hamil_phase_alpha_quad_y_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[2] =
            tensor_hamil_phase_alpha_quad_z_kernels[kernel_index].kernels[poly_order];
        } else {
          up->hamil_alpha_shared[0] =
            tensor_hamil_phase_ho_alpha_quad_x_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[1] =
            tensor_hamil_phase_ho_alpha_quad_y_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[2] =
            tensor_hamil_phase_ho_alpha_quad_z_kernels[kernel_index].kernels[poly_order];
        }
      } else if (model_id == GKYL_MODEL_TRIAD) {
        // Triad bracket on the tensor p=1 hybrid (sparse or dense velocity-space
        // Hamiltonian): per-node inverse velocity-map Jacobians of the C^1 cubic map.
        if (use_lo) {
          up->hamil_alpha_shared[0] =
            hamil_sparse ?
              tensor_hamil_vel_sparse_alpha_quad_x_kernels[kernel_index].kernels[poly_order] :
              tensor_hamil_vel_dense_alpha_quad_x_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[1] =
            hamil_sparse ?
              tensor_hamil_vel_sparse_alpha_quad_y_kernels[kernel_index].kernels[poly_order] :
              tensor_hamil_vel_dense_alpha_quad_y_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[2] =
            hamil_sparse ?
              tensor_hamil_vel_sparse_alpha_quad_z_kernels[kernel_index].kernels[poly_order] :
              tensor_hamil_vel_dense_alpha_quad_z_kernels[kernel_index].kernels[poly_order];
        } else {
          up->hamil_alpha_shared[0] =
            hamil_sparse ?
              tensor_hamil_vel_sparse_ho_alpha_quad_x_kernels[kernel_index].kernels[poly_order] :
              tensor_hamil_vel_dense_ho_alpha_quad_x_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[1] =
            hamil_sparse ?
              tensor_hamil_vel_sparse_ho_alpha_quad_y_kernels[kernel_index].kernels[poly_order] :
              tensor_hamil_vel_dense_ho_alpha_quad_y_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[2] =
            hamil_sparse ?
              tensor_hamil_vel_sparse_ho_alpha_quad_z_kernels[kernel_index].kernels[poly_order] :
              tensor_hamil_vel_dense_ho_alpha_quad_z_kernels[kernel_index].kernels[poly_order];
        }
      }

      break;

    default:
      assert(false);
      break;
  }
  // The device kernel inlines the surface assembly; the whole-surface CPU
  // dispatch pointer is not used on the device.
  up->conf_flux_surf = 0;
  // Size of the force-factor buffers: the largest term count over the
  // directions (copied back to the host, which sizes the launch's shared memory).
  up->alpha_nterms_max = 0;
  for (int d = 0; d < cdim; ++d) {
    up->alpha_nterms_max = GKYL_MAX2(up->alpha_nterms_max, conf_flux_surf_alpha_nterms(up, d));
  }
}

gkyl_dg_vlasov_conf_flux_surf *
gkyl_dg_vlasov_conf_flux_surf_cu_dev_inew(const struct gkyl_dg_vlasov_conf_flux_surf_inp *inp)
{
  struct gkyl_dg_vlasov_conf_flux_surf *up =
    (struct gkyl_dg_vlasov_conf_flux_surf *)gkyl_malloc(sizeof(*up));

  int cdim = inp->conf_basis->ndim, pdim = inp->phase_basis->ndim, vdim = pdim - cdim;
  int poly_order = inp->conf_basis->poly_order;

  up->phase_grid = *inp->phase_grid;
  up->cdim = cdim;
  up->pdim = pdim;
  up->use_gpu = true;

  // Are we skipping cells with small phase space density?
  if (inp->skip_cell_thresh > 0.0) {
    up->skip_cell_thresh = inp->skip_cell_thresh * pow(sqrt(2.0), pdim);
  } else {
    up->skip_cell_thresh = -1.0;
  }

  // Determine Hamiltonian dimensionality and index offset for indexing Hamiltonian
  // from an input phase space index.
  up->hamil_range = *inp->hamil_range;
  if (inp->model_id == GKYL_MODEL_DEFAULT || inp->model_id == GKYL_MODEL_SR ||
      inp->model_id == GKYL_MODEL_TRIAD) {
    up->hamil_dim = vdim;
    up->hamil_offset = cdim;
  } else {
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

  // Surface node counts and modal size of the stored flux (the advance wrapper
  // sizes the 2D (cells x nodes) kernel launch from them).
  conf_flux_surf_num_nodes(
    gkyl_basis_phase_kernel_type(inp->conf_basis, inp->phase_basis), cdim, vdim, poly_order,
    inp->use_lo, &up->num_nodes_conf, &up->num_nodes_vel
  );
  up->num_surf_basis = conf_flux_surf_num_surf_basis(
    gkyl_basis_phase_kernel_type(inp->conf_basis, inp->phase_basis), cdim, vdim, poly_order
  );
  assert(up->num_surf_basis <= up->num_nodes_conf * up->num_nodes_vel);

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_dg_vlasov_conf_flux_surf *up_cu =
    (struct gkyl_dg_vlasov_conf_flux_surf *)gkyl_cu_malloc(sizeof(*up_cu));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_dg_vlasov_conf_flux_surf), GKYL_CU_MEMCPY_H2D);

  gkyl_dg_vlasov_conf_flux_surf_set_cu_dev_ptrs<<<1, 1>>>(
    up_cu, gkyl_basis_phase_kernel_type(inp->conf_basis, inp->phase_basis), cdim, vdim, poly_order,
    inp->model_id, inp->hamil_id, inp->use_lo
  );

  // The device selection knows the producer's term count; the host needs its
  // maximum to size the launch's dynamic shared memory.
  gkyl_cu_memcpy(&up->alpha_nterms_max, &up_cu->alpha_nterms_max, sizeof(int), GKYL_CU_MEMCPY_D2H);
  assert(
    up->alpha_nterms_max * (up->num_nodes_conf + up->num_nodes_vel) <=
    GKYL_VLASOV_CONF_FLUX_SURF_MAX_ALPHA_FACTORS
  );

  // set parent on_dev pointer
  up->on_dev = up_cu;

  return up;
}
