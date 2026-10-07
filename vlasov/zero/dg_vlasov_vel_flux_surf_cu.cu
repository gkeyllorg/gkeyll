/* -*- c++ -*- */

#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_vlasov_vel_flux_surf.h>
#include <gkyl_dg_vlasov_vel_flux_surf_priv.h>
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
gkyl_dg_vlasov_vel_flux_surf_advance_cu_kernel(
  struct gkyl_dg_vlasov_vel_flux_surf *up, struct gkyl_range conf_range,
  struct gkyl_range phase_range, const struct gkyl_array *vmap, const struct gkyl_array *jacob_pos,
  const struct gkyl_array *jacob_vel_surf, const struct gkyl_array *poisson_tensor_conf,
  const struct gkyl_array *hamil, const struct gkyl_array *qmem, const struct gkyl_array *pot_tot,
  const struct gkyl_array *rad, const struct gkyl_array *fin, struct gkyl_array *cflrate,
  struct gkyl_array *vel_flux_surf
)
{
  // Per-node |alpha| values for the block, reduced to alpha_max per (cell, dir)
  // by the last node thread of each cell: the projection threads are the first
  // num_surf_basis node threads, so the reduction runs in a different warp and
  // overlaps with them instead of serializing behind them.
  __shared__ double alpha_smem[GKYL_DEFAULT_NUM_THREADS];
  // Stage-1 arrays of the sum-factorized nodal f evaluation, G[j*NA + a] per
  // side, filled cooperatively: each cell owns a block of blockDim.y entries
  // per side (there are at most as many (inner node, outer shape) items as
  // surface nodes), so one thread per item fills them and every node thread
  // of the cell then reads its stage-2 dot product from shared memory.
  __shared__ double G_smem[2 * GKYL_DEFAULT_NUM_THREADS];
  // This direction's nodal Lax flux, one entry per surface node of each cell
  // of the block; projected onto the surface modal basis by the lax_prj stage
  // (one thread per surface mode) and only then written to the flux array.
  __shared__ double F_smem[GKYL_DEFAULT_NUM_THREADS];
  // Shared factors of the force producers (dynamic shared memory, sized by the
  // launcher): per cell, the outer factors O[t*NO + i] of every term t of every
  // producer, then the inner factors I[t*NI + j]. The force at a node is the
  // dot product of the two over all terms.
  extern __shared__ double alpha_factors_smem[];

  int pdim = up->pdim;
  int cdim = up->cdim;
  int vdim = pdim - cdim;

  // 2D thread grid: linc2 indexes the surface node (i_node: configuration-space
  // nodes, j_node: transverse velocity-space nodes, matching the CPU dispatch's
  // i-major/j-minor node order); threads in x index the phase-space cell.
  // blockDim.y == num_nodes_conf*num_nodes_vel by construction of the launch.
  int num_nodes_vel = up->num_nodes_vel;
  int num_nodes = up->num_nodes_conf * num_nodes_vel;
  int linc2 = threadIdx.y;
  int i_node = linc2 / num_nodes_vel;
  int j_node = linc2 % num_nodes_vel;
  double *G_l = &G_smem[2 * blockDim.y * threadIdx.x];
  double *G_c = G_l + blockDim.y;
  double *F_cell = &F_smem[blockDim.y * threadIdx.x];
  const int NO = up->num_nodes_conf, NI = up->num_nodes_vel;
  double *O = &alpha_factors_smem[threadIdx.x * up->alpha_nterms_max * (NO + NI)];
  double *I = O + up->alpha_nterms_max * NO;

  // No grid-stride loop: the launch covers phase_range.volume exactly so the
  // __syncthreads() barriers below are reached uniformly by every thread in
  // the block (threads past the end of the range only skip the guarded work).
  unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x;
  bool valid = linc1 < phase_range.volume;

  int idx[GKYL_MAX_DIM], idx_l[GKYL_MAX_DIM], idx_vel[GKYL_MAX_DIM], idx_hamil[GKYL_MAX_DIM];
  double xcC[GKYL_MAX_DIM];
  long vidx = 0;
  const double *poisson_tensor_conf_d = 0, *hamil_d = 0, *qmem_d = 0, *pot_tot_d = 0, *rad_d = 0;
  const double *jacob_pos_d = 0, *vmap_d = 0, *f_c = 0;
  double *cflrate_d = 0, *flux = 0;
  if (valid) {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&phase_range, linc1, idx);
    long cidx = gkyl_range_idx(&conf_range, idx);
    long pidx = gkyl_range_idx(&phase_range, idx);

    for (int i = 0; i < vdim; ++i) {
      idx_vel[i] = idx[cdim + i];
    }
    vidx = gkyl_range_idx(&up->vel_range, idx_vel);

    for (int i = 0; i < up->hamil_dim; ++i) {
      idx_hamil[i] = idx[up->hamil_offset + i];
    }
    long hidx = gkyl_range_idx(&up->hamil_range, idx_hamil);

    // Grab the cell center location for NC bracket calculation
    gkyl_rect_grid_cell_center(&up->phase_grid, idx, xcC);

    poisson_tensor_conf_d = (const double *)gkyl_array_cfetch(poisson_tensor_conf, cidx);
    hamil_d = (const double *)gkyl_array_cfetch(hamil, hidx);
    qmem_d = qmem ? (const double *)gkyl_array_cfetch(qmem, cidx) : 0;
    pot_tot_d = pot_tot ? (const double *)gkyl_array_cfetch(pot_tot, cidx) : 0;
    rad_d = rad ? (const double *)gkyl_array_cfetch(rad, vidx) : 0;
    jacob_pos_d = jacob_pos ? (const double *)gkyl_array_cfetch(jacob_pos, cidx) : 0;
    vmap_d = (const double *)gkyl_array_cfetch(vmap, vidx);
    f_c = (const double *)gkyl_array_cfetch(fin, pidx);
    cflrate_d = (double *)gkyl_array_fetch(cflrate, pidx);
    flux = (double *)gkyl_array_fetch(vel_flux_surf, pidx);
  }

  // Each cell owns *lower* fluxes in each velocity-space direction.
  // So we need the distribution function in our current cell, and the cell
  // one index lower in each direction. If we are at the lower velocity-space
  // edge, we call special kernels because we have *no* ghost cells in velocity
  // space, so we cannot index the distribution function in the lower direction
  // along that velocity-space edge.
  for (int dir = 0; dir < vdim; ++dir) {
    bool edge = valid && (idx[cdim + dir] == phase_range.lower[cdim + dir]);
    const double *f_l = 0, *jacob_vel_l_d = 0, *jacob_vel_d = 0;
    double node_alpha = 0.0;
    int nt = 0;
    if (valid && !edge) {
      gkyl_copy_int_arr(pdim, idx, idx_l);
      idx_l[cdim + dir] = idx_l[cdim + dir] - 1;
      long pidx_l = gkyl_range_idx(&phase_range, idx_l);
      f_l = (const double *)gkyl_array_cfetch(fin, pidx_l);
      // Velocity-space Jacobian of the lower neighbor in this direction, for
      // the minimum-Jacobian time-step estimate of the C^0 linear map.
      int idx_vel_l[GKYL_MAX_DIM];
      for (int i = 0; i < vdim; ++i) {
        idx_vel_l[i] = idx_l[cdim + i];
      }
      long vidx_l = gkyl_range_idx(&up->vel_range, idx_vel_l);
      jacob_vel_l_d = (const double *)gkyl_array_cfetch(jacob_vel_surf, vidx_l);
      jacob_vel_d = (const double *)gkyl_array_cfetch(jacob_vel_surf, vidx);

      // Shared work of the cell, cooperatively over its node threads: the
      // outer/inner factors of every force producer (one item per outer or
      // inner node; each producer's terms follow the previous one's) and stage
      // 1 of the nodal f evaluation (one item per (inner node, outer shape)).
      const double *dx = up->phase_grid.dx;
      nt += up->hamil_alpha_shared[dir](
        linc2, blockDim.y, xcC, dx, vmap_d, jacob_pos_d, jacob_vel_d, poisson_tensor_conf_d,
        hamil_d, O + nt * NO, I + nt * NI
      );
      nt += up->E_alpha_shared[dir](linc2, blockDim.y, dx, qmem_d, O + nt * NO, I + nt * NI);
      nt += up->phi_alpha_shared[dir](
        linc2, blockDim.y, dx, jacob_pos_d, pot_tot_d, O + nt * NO, I + nt * NI
      );
      nt += up->B_alpha_shared[dir](
        linc2, blockDim.y, dx, jacob_vel_d, hamil_d, qmem_d, O + nt * NO, I + nt * NI
      );
      nt += up->rad_alpha_shared[dir](linc2, blockDim.y, dx, rad_d, O + nt * NO, I + nt * NI);
      up->lax_g[dir](linc2, f_l, f_c, G_l, G_c);
    }
    __syncthreads();

    if (valid && !edge) {
      // Per-node work: the force at this thread's surface node is the dot
      // product of the shared outer and inner factors over all terms; nodal
      // Lax flux from the stage-1 arrays into the nodal buffer, keeping |alpha|
      // for the CFL reduction below.
      double alpha = 0.0;
      for (int t = 0; t < nt; ++t) {
        alpha += O[t * NO + i_node] * I[t * NI + j_node];
      }
      node_alpha = up->lax_flux_nodal[dir](i_node, j_node, jacob_vel_d, alpha, G_l, G_c, F_cell);
    }
    alpha_smem[threadIdx.x + blockDim.x * threadIdx.y] = node_alpha;
    __syncthreads();

    if (valid && !edge) {
      // Stage 3: project this thread's surface mode of the cell's nodal flux
      // onto the surface modal basis (a no-op for linc2 >= num_surf_basis).
      up->lax_prj[dir](linc2, F_cell, flux);
    }
    if (valid && threadIdx.y == blockDim.y - 1) {
      if (edge) {
        cflrate_d[0] += up->vel_flux_surf_edge(
          up, dir, xcC, up->phase_grid.dx,
          jacob_vel_surf ? (const double *)gkyl_array_cfetch(jacob_vel_surf, vidx) : 0,
          poisson_tensor_conf_d, hamil_d, qmem_d, pot_tot_d, rad_d, f_c, flux
        );
      } else {
        // Reduce alpha_max in the CPU dispatch's node order so the fmax chain,
        // and hence the CFL estimate, matches the CPU loop exactly.
        double alpha_max = 0.0;
        for (int n = 0; n < num_nodes; ++n) {
          alpha_max = fmax(alpha_max, alpha_smem[threadIdx.x + blockDim.x * n]);
        }
        double cfl = up->lax_cfl[dir](up->phase_grid.dx, jacob_vel_l_d, jacob_vel_d, alpha_max);
        // Always compute the flux, but if we are below threshold, ignore the stable time step estimate.
        if (fabs(f_l[0]) < up->skip_cell_thresh && fabs(f_c[0]) < up->skip_cell_thresh) {
          cfl = 0.0;
        }
        cflrate_d[0] += cfl;
      }
    }
    // alpha_smem, G_smem and F_smem are reused by the next direction.
    __syncthreads();
  }
}

void
gkyl_dg_vlasov_vel_flux_surf_advance_cu(
  struct gkyl_dg_vlasov_vel_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *poisson_tensor_conf,
  const struct gkyl_array *hamil, const struct gkyl_array *qmem, const struct gkyl_array *pot_tot,
  const struct gkyl_array *rad, const struct gkyl_array *fin, struct gkyl_array *cflrate,
  struct gkyl_array *vel_flux_surf
)
{
  // 2D thread grid: parallelize over surface nodes as well as phase-space cells.
  int num_nodes = up->num_nodes_conf * up->num_nodes_vel;
  assert(num_nodes <= GKYL_DEFAULT_NUM_THREADS);
  dim3 dimGrid, dimBlock;
  gkyl_parallelize_components_kernel_launch_dims(&dimGrid, &dimBlock, *phase_range, num_nodes);
  // Dynamic shared memory for the force producers' shared factors: per cell of
  // the block, alpha_nterms_max outer + inner factor rows.
  size_t alpha_smem =
    sizeof(double) * dimBlock.x * up->alpha_nterms_max * (up->num_nodes_conf + up->num_nodes_vel);
  gkyl_dg_vlasov_vel_flux_surf_advance_cu_kernel<<<dimGrid, dimBlock, alpha_smem>>>(
    up->on_dev, *conf_range, *phase_range, up->vmap->on_dev, up->jacob_pos->on_dev,
    up->jacob_vel_surf->on_dev, poisson_tensor_conf->on_dev, hamil->on_dev, qmem ? qmem->on_dev : 0,
    pot_tot ? pot_tot->on_dev : 0, rad ? rad->on_dev : 0, fin->on_dev, cflrate->on_dev,
    vel_flux_surf->on_dev
  );
}

// CUDA kernel to set device pointers to canonical pb vars kernel functions
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void
gkyl_dg_vlasov_vel_flux_surf_set_cu_dev_ptrs(
  struct gkyl_dg_vlasov_vel_flux_surf *up, enum gkyl_basis_type b_type, int cdim, int vdim,
  int poly_order, enum gkyl_model_id model_id, enum gkyl_hamil_id hamil_id, bool has_E,
  bool has_phi, bool has_B, bool has_rad, bool use_lo
)
{
  // Sparse (separable) vs. dense velocity-space Hamiltonian kernel selection.
  bool hamil_sparse = (hamil_id == GKYL_HAMIL_VEL_SPARSE);
  // By default, we have no forces from Hamiltonian, E, B, or phi.
  for (int d = 0; d < vdim; ++d) {
    up->hamil_alpha_shared[d] = no_hamil_alpha_shared;
    up->E_alpha_shared[d] = no_E_alpha_shared;
    up->phi_alpha_shared[d] = no_phi_alpha_shared;
    up->B_alpha_shared[d] = no_B_alpha_shared;
    up->rad_alpha_shared[d] = no_rad_alpha_shared;
  }

  int kernel_index = cv_index[cdim].vdim[vdim];
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      if (use_lo) {
        up->lax_g[0] = ser_lax_flux_nodal_vx_g_kernels[kernel_index].kernels[poly_order];
        up->lax_prj[0] = ser_lax_flux_nodal_vx_prj_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[0] = ser_lax_flux_nodal_vx_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[0] = ser_lax_flux_nodal_vx_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_g[1] = ser_lax_flux_nodal_vy_g_kernels[kernel_index].kernels[poly_order];
        up->lax_prj[1] = ser_lax_flux_nodal_vy_prj_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[1] = ser_lax_flux_nodal_vy_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[1] = ser_lax_flux_nodal_vy_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_g[2] = ser_lax_flux_nodal_vz_g_kernels[kernel_index].kernels[poly_order];
        up->lax_prj[2] = ser_lax_flux_nodal_vz_prj_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[2] = ser_lax_flux_nodal_vz_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[2] = ser_lax_flux_nodal_vz_cfl_kernels[kernel_index].kernels[poly_order];
      } else {
        up->lax_g[0] = ser_ho_lax_flux_nodal_vx_g_kernels[kernel_index].kernels[poly_order];
        up->lax_prj[0] = ser_ho_lax_flux_nodal_vx_prj_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[0] = ser_ho_lax_flux_nodal_vx_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[0] = ser_ho_lax_flux_nodal_vx_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_g[1] = ser_ho_lax_flux_nodal_vy_g_kernels[kernel_index].kernels[poly_order];
        up->lax_prj[1] = ser_ho_lax_flux_nodal_vy_prj_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[1] = ser_ho_lax_flux_nodal_vy_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[1] = ser_ho_lax_flux_nodal_vy_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_g[2] = ser_ho_lax_flux_nodal_vz_g_kernels[kernel_index].kernels[poly_order];
        up->lax_prj[2] = ser_ho_lax_flux_nodal_vz_prj_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[2] = ser_ho_lax_flux_nodal_vz_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[2] = ser_ho_lax_flux_nodal_vz_cfl_kernels[kernel_index].kernels[poly_order];
      }

      // Only have Hamiltonian forces in general geometry.
      if (model_id == GKYL_MODEL_CANONICAL_PB || model_id == GKYL_MODEL_CANONICAL_PB_GR) {
        if (use_lo) {
          up->hamil_alpha_shared[0] =
            ser_hamil_phase_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[1] =
            ser_hamil_phase_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[2] =
            ser_hamil_phase_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        } else {
          up->hamil_alpha_shared[0] =
            ser_hamil_phase_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[1] =
            ser_hamil_phase_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[2] =
            ser_hamil_phase_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }
      } else if (model_id == GKYL_MODEL_TRIAD) {
        if (use_lo) {
          up->hamil_alpha_shared[0] =
            hamil_sparse ?
              ser_nc_hamil_vel_sparse_alpha_quad_vx_kernels[kernel_index].kernels[poly_order] :
              ser_nc_hamil_vel_dense_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[1] =
            hamil_sparse ?
              ser_nc_hamil_vel_sparse_alpha_quad_vy_kernels[kernel_index].kernels[poly_order] :
              ser_nc_hamil_vel_dense_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[2] =
            hamil_sparse ?
              ser_nc_hamil_vel_sparse_alpha_quad_vz_kernels[kernel_index].kernels[poly_order] :
              ser_nc_hamil_vel_dense_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        } else {
          up->hamil_alpha_shared[0] =
            hamil_sparse ?
              ser_nc_hamil_vel_sparse_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order] :
              ser_nc_hamil_vel_dense_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[1] =
            hamil_sparse ?
              ser_nc_hamil_vel_sparse_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order] :
              ser_nc_hamil_vel_dense_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[2] =
            hamil_sparse ?
              ser_nc_hamil_vel_sparse_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order] :
              ser_nc_hamil_vel_dense_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }
      } else if (model_id == GKYL_MODEL_TRIAD_GR) {
        if (use_lo) {
          up->hamil_alpha_shared[0] =
            ser_nc_hamil_phase_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[1] =
            ser_nc_hamil_phase_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[2] =
            ser_nc_hamil_phase_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        } else {
          up->hamil_alpha_shared[0] =
            ser_nc_hamil_phase_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[1] =
            ser_nc_hamil_phase_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[2] =
            ser_nc_hamil_phase_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }
      }

      if (use_lo) {
        if (has_E) {
          up->E_alpha_shared[0] = ser_E_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->E_alpha_shared[1] = ser_E_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->E_alpha_shared[2] = ser_E_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }

        if (has_phi) {
          up->phi_alpha_shared[0] = ser_phi_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->phi_alpha_shared[1] = ser_phi_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->phi_alpha_shared[2] = ser_phi_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }

        if (has_B) {
          if (hamil_id == GKYL_HAMIL_PHASE) {
            up->B_alpha_shared[0] =
              ser_B_hamil_phase_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
            up->B_alpha_shared[1] =
              ser_B_hamil_phase_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
            up->B_alpha_shared[2] =
              ser_B_hamil_phase_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
          } else {
            up->B_alpha_shared[0] =
              hamil_sparse ?
                ser_B_hamil_vel_sparse_alpha_quad_vx_kernels[kernel_index].kernels[poly_order] :
                ser_B_hamil_vel_dense_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
            up->B_alpha_shared[1] =
              hamil_sparse ?
                ser_B_hamil_vel_sparse_alpha_quad_vy_kernels[kernel_index].kernels[poly_order] :
                ser_B_hamil_vel_dense_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
            up->B_alpha_shared[2] =
              hamil_sparse ?
                ser_B_hamil_vel_sparse_alpha_quad_vz_kernels[kernel_index].kernels[poly_order] :
                ser_B_hamil_vel_dense_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
          }
        }

        if (has_rad) {
          up->rad_alpha_shared[0] = ser_rad_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->rad_alpha_shared[1] = ser_rad_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->rad_alpha_shared[2] = ser_rad_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }
      } else {
        if (has_E) {
          up->E_alpha_shared[0] = ser_E_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->E_alpha_shared[1] = ser_E_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->E_alpha_shared[2] = ser_E_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }

        if (has_phi) {
          up->phi_alpha_shared[0] =
            ser_phi_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->phi_alpha_shared[1] =
            ser_phi_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->phi_alpha_shared[2] =
            ser_phi_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }

        if (has_B) {
          if (hamil_id == GKYL_HAMIL_PHASE) {
            up->B_alpha_shared[0] =
              ser_B_ho_hamil_phase_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
            up->B_alpha_shared[1] =
              ser_B_ho_hamil_phase_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
            up->B_alpha_shared[2] =
              ser_B_ho_hamil_phase_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
          } else {
            up->B_alpha_shared[0] =
              hamil_sparse ?
                ser_B_ho_hamil_vel_sparse_alpha_quad_vx_kernels[kernel_index].kernels[poly_order] :
                ser_B_ho_hamil_vel_dense_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
            up->B_alpha_shared[1] =
              hamil_sparse ?
                ser_B_ho_hamil_vel_sparse_alpha_quad_vy_kernels[kernel_index].kernels[poly_order] :
                ser_B_ho_hamil_vel_dense_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
            up->B_alpha_shared[2] =
              hamil_sparse ?
                ser_B_ho_hamil_vel_sparse_alpha_quad_vz_kernels[kernel_index].kernels[poly_order] :
                ser_B_ho_hamil_vel_dense_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
          }
        }

        if (has_rad) {
          up->rad_alpha_shared[0] =
            ser_rad_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->rad_alpha_shared[1] =
            ser_rad_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->rad_alpha_shared[2] =
            ser_rad_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }
      }

      break;

    case GKYL_BASIS_MODAL_TENSOR:
      // Only the tensor p=1 hybrid has distinct lo/ho surface variants; the
      // plain and ho lists share the (high-order by design) kernels at p>1.
      if (use_lo) {
        up->lax_g[0] = tensor_lax_flux_nodal_vx_g_kernels[kernel_index].kernels[poly_order];
        up->lax_prj[0] = tensor_lax_flux_nodal_vx_prj_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[0] = tensor_lax_flux_nodal_vx_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[0] = tensor_lax_flux_nodal_vx_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_g[1] = tensor_lax_flux_nodal_vy_g_kernels[kernel_index].kernels[poly_order];
        up->lax_prj[1] = tensor_lax_flux_nodal_vy_prj_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[1] = tensor_lax_flux_nodal_vy_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[1] = tensor_lax_flux_nodal_vy_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_g[2] = tensor_lax_flux_nodal_vz_g_kernels[kernel_index].kernels[poly_order];
        up->lax_prj[2] = tensor_lax_flux_nodal_vz_prj_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[2] = tensor_lax_flux_nodal_vz_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[2] = tensor_lax_flux_nodal_vz_cfl_kernels[kernel_index].kernels[poly_order];
      } else {
        up->lax_g[0] = tensor_ho_lax_flux_nodal_vx_g_kernels[kernel_index].kernels[poly_order];
        up->lax_prj[0] = tensor_ho_lax_flux_nodal_vx_prj_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[0] =
          tensor_ho_lax_flux_nodal_vx_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[0] = tensor_ho_lax_flux_nodal_vx_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_g[1] = tensor_ho_lax_flux_nodal_vy_g_kernels[kernel_index].kernels[poly_order];
        up->lax_prj[1] = tensor_ho_lax_flux_nodal_vy_prj_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[1] =
          tensor_ho_lax_flux_nodal_vy_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[1] = tensor_ho_lax_flux_nodal_vy_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_g[2] = tensor_ho_lax_flux_nodal_vz_g_kernels[kernel_index].kernels[poly_order];
        up->lax_prj[2] = tensor_ho_lax_flux_nodal_vz_prj_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[2] =
          tensor_ho_lax_flux_nodal_vz_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[2] = tensor_ho_lax_flux_nodal_vz_cfl_kernels[kernel_index].kernels[poly_order];
      }

      // Only have Hamiltonian forces in general geometry: the p=1 tensor
      // hybrid is the only tensor basis with a phase-space Hamiltonian
      // representation.
      if (model_id == GKYL_MODEL_CANONICAL_PB || model_id == GKYL_MODEL_CANONICAL_PB_GR) {
        if (use_lo) {
          up->hamil_alpha_shared[0] =
            tensor_hamil_phase_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[1] =
            tensor_hamil_phase_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[2] =
            tensor_hamil_phase_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        } else {
          up->hamil_alpha_shared[0] =
            tensor_hamil_phase_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[1] =
            tensor_hamil_phase_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[2] =
            tensor_hamil_phase_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }
      } else if (model_id == GKYL_MODEL_TRIAD) {
        // Triad bracket on the tensor p=1 hybrid (sparse or dense velocity-space
        // Hamiltonian): per-node inverse velocity-map Jacobians of the C^1 cubic
        // map and the cubic vmap in the omega = v.Pi momentum factor.
        if (use_lo) {
          up->hamil_alpha_shared[0] =
            hamil_sparse ?
              tensor_nc_hamil_vel_sparse_alpha_quad_vx_kernels[kernel_index].kernels[poly_order] :
              tensor_nc_hamil_vel_dense_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[1] =
            hamil_sparse ?
              tensor_nc_hamil_vel_sparse_alpha_quad_vy_kernels[kernel_index].kernels[poly_order] :
              tensor_nc_hamil_vel_dense_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[2] =
            hamil_sparse ?
              tensor_nc_hamil_vel_sparse_alpha_quad_vz_kernels[kernel_index].kernels[poly_order] :
              tensor_nc_hamil_vel_dense_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        } else {
          up->hamil_alpha_shared[0] =
            hamil_sparse ?
              tensor_nc_hamil_vel_sparse_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order] :
              tensor_nc_hamil_vel_dense_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[1] =
            hamil_sparse ?
              tensor_nc_hamil_vel_sparse_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order] :
              tensor_nc_hamil_vel_dense_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[2] =
            hamil_sparse ?
              tensor_nc_hamil_vel_sparse_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order] :
              tensor_nc_hamil_vel_dense_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }
      } else if (model_id == GKYL_MODEL_TRIAD_GR) {
        if (use_lo) {
          up->hamil_alpha_shared[0] =
            tensor_nc_hamil_phase_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[1] =
            tensor_nc_hamil_phase_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[2] =
            tensor_nc_hamil_phase_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        } else {
          up->hamil_alpha_shared[0] =
            tensor_nc_hamil_phase_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[1] =
            tensor_nc_hamil_phase_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_shared[2] =
            tensor_nc_hamil_phase_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }
      }

      if (use_lo) {
        if (has_E) {
          up->E_alpha_shared[0] = tensor_E_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->E_alpha_shared[1] = tensor_E_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->E_alpha_shared[2] = tensor_E_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }

        if (has_phi) {
          up->phi_alpha_shared[0] =
            tensor_phi_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->phi_alpha_shared[1] =
            tensor_phi_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->phi_alpha_shared[2] =
            tensor_phi_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }

        if (has_B) {
          // No phase-space Hamiltonian magnetic-force kernels for tensor p>1;
          // phase runs keep the no-op defaults set above. (The tensor p=1
          // hybrid does have B_hamil_phase kernels; hook them in with the
          // phase-Hamiltonian hybrid support.)
          if (hamil_id == GKYL_HAMIL_PHASE) {
            up->B_alpha_shared[0] =
              tensor_B_hamil_phase_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
            up->B_alpha_shared[1] =
              tensor_B_hamil_phase_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
            up->B_alpha_shared[2] =
              tensor_B_hamil_phase_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
          } else {
            up->B_alpha_shared[0] =
              hamil_sparse ?
                tensor_B_hamil_vel_sparse_alpha_quad_vx_kernels[kernel_index].kernels[poly_order] :
                tensor_B_hamil_vel_dense_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
            up->B_alpha_shared[1] =
              hamil_sparse ?
                tensor_B_hamil_vel_sparse_alpha_quad_vy_kernels[kernel_index].kernels[poly_order] :
                tensor_B_hamil_vel_dense_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
            up->B_alpha_shared[2] =
              hamil_sparse ?
                tensor_B_hamil_vel_sparse_alpha_quad_vz_kernels[kernel_index].kernels[poly_order] :
                tensor_B_hamil_vel_dense_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
          }
        }

        if (has_rad) {
          up->rad_alpha_shared[0] =
            tensor_rad_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->rad_alpha_shared[1] =
            tensor_rad_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->rad_alpha_shared[2] =
            tensor_rad_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }
      } else {
        if (has_E) {
          up->E_alpha_shared[0] =
            tensor_E_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->E_alpha_shared[1] =
            tensor_E_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->E_alpha_shared[2] =
            tensor_E_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }

        if (has_phi) {
          up->phi_alpha_shared[0] =
            tensor_phi_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->phi_alpha_shared[1] =
            tensor_phi_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->phi_alpha_shared[2] =
            tensor_phi_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }

        if (has_B) {
          // No phase-space Hamiltonian magnetic-force kernels for tensor p>1;
          // phase runs keep the no-op defaults set above. (The tensor p=1
          // hybrid does have B_hamil_phase kernels; hook them in with the
          // phase-Hamiltonian hybrid support.)
          if (hamil_id == GKYL_HAMIL_PHASE) {
            up->B_alpha_shared[0] =
              tensor_B_ho_hamil_phase_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
            up->B_alpha_shared[1] =
              tensor_B_ho_hamil_phase_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
            up->B_alpha_shared[2] =
              tensor_B_ho_hamil_phase_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
          } else {
            up->B_alpha_shared[0] =
              hamil_sparse ?
                tensor_B_ho_hamil_vel_sparse_alpha_quad_vx_kernels[kernel_index]
                  .kernels[poly_order] :
                tensor_B_ho_hamil_vel_dense_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
            up->B_alpha_shared[1] =
              hamil_sparse ?
                tensor_B_ho_hamil_vel_sparse_alpha_quad_vy_kernels[kernel_index]
                  .kernels[poly_order] :
                tensor_B_ho_hamil_vel_dense_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
            up->B_alpha_shared[2] =
              hamil_sparse ?
                tensor_B_ho_hamil_vel_sparse_alpha_quad_vz_kernels[kernel_index]
                  .kernels[poly_order] :
                tensor_B_ho_hamil_vel_dense_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
          }
        }

        if (has_rad) {
          up->rad_alpha_shared[0] =
            tensor_rad_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->rad_alpha_shared[1] =
            tensor_rad_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->rad_alpha_shared[2] =
            tensor_rad_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }
      }

      break;

    default:
      assert(false);
      break;
  }
  // The device kernel inlines the surface assembly; the whole-surface CPU
  // dispatch pointer is not used on the device.
  up->vel_flux_surf = 0;
  // Currently only support zero-flux boundary, so edge velocity flux an empty function (flux = 0.0).
  up->vel_flux_surf_edge = no_vel_flux_surf_edge;
  // Size of the force-factor buffers: the largest total term count over the
  // directions (copied back to the host, which sizes the launch's shared memory).
  up->alpha_nterms_max = 0;
  for (int d = 0; d < vdim; ++d) {
    up->alpha_nterms_max = GKYL_MAX2(up->alpha_nterms_max, vel_flux_surf_alpha_nterms(up, d));
  }
}

gkyl_dg_vlasov_vel_flux_surf *
gkyl_dg_vlasov_vel_flux_surf_cu_dev_inew(const struct gkyl_dg_vlasov_vel_flux_surf_inp *inp)
{
  struct gkyl_dg_vlasov_vel_flux_surf *up =
    (struct gkyl_dg_vlasov_vel_flux_surf *)gkyl_malloc(sizeof(*up));

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
  if (inp->hamil_id == GKYL_HAMIL_PHASE) {
    up->hamil_dim = pdim;
    up->hamil_offset = 0;
  } else {
    up->hamil_dim = vdim;
    up->hamil_offset = cdim;
  }
  // The velocity map is required: it provides the velocity-space Jacobian at
  // surface quadrature points and the velocity-space range used to index it.
  // The host pointers below are not dereferenced on device (the advance wrapper
  // passes the raw device array pointer to the kernel as an argument).
  assert(inp->vel_map);
  assert(inp->pos_map);
  up->vel_range = inp->vel_map->local_vel;
  up->vel_map = 0;
  up->jacob_vel_surf = 0;
  up->vmap = 0;
  up->pos_map = 0;
  up->jacob_pos = 0;

  // Surface node counts and modal size of the stored flux (the advance wrapper
  // sizes the 2D (cells x nodes) kernel launch from them).
  vel_flux_surf_num_nodes(
    gkyl_basis_phase_kernel_type(inp->conf_basis, inp->phase_basis), cdim, vdim, poly_order,
    inp->use_lo, &up->num_nodes_conf, &up->num_nodes_vel
  );
  up->num_surf_basis = vel_flux_surf_num_surf_basis(
    gkyl_basis_phase_kernel_type(inp->conf_basis, inp->phase_basis), cdim, vdim, poly_order
  );
  assert(up->num_surf_basis <= up->num_nodes_conf * up->num_nodes_vel);

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_dg_vlasov_vel_flux_surf *up_cu =
    (struct gkyl_dg_vlasov_vel_flux_surf *)gkyl_cu_malloc(sizeof(*up_cu));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_dg_vlasov_vel_flux_surf), GKYL_CU_MEMCPY_H2D);

  gkyl_dg_vlasov_vel_flux_surf_set_cu_dev_ptrs<<<1, 1>>>(
    up_cu, gkyl_basis_phase_kernel_type(inp->conf_basis, inp->phase_basis), cdim, vdim, poly_order,
    inp->model_id, inp->hamil_id, inp->has_E, inp->has_phi, inp->has_B, inp->has_rad, inp->use_lo
  );

  // The device selection knows the producers' term counts; the host needs
  // their maximum to size the launch's dynamic shared memory.
  gkyl_cu_memcpy(&up->alpha_nterms_max, &up_cu->alpha_nterms_max, sizeof(int), GKYL_CU_MEMCPY_D2H);
  assert(
    up->alpha_nterms_max * (up->num_nodes_conf + up->num_nodes_vel) <=
    GKYL_VLASOV_VEL_FLUX_SURF_MAX_ALPHA_FACTORS
  );

  // set parent on_dev pointer
  up->on_dev = up_cu;

  // Host-side updater stores the acquired map and host array pointers.
  up->vel_map = gkyl_vlasov_velocity_map_acquire(inp->vel_map);
  up->jacob_vel_surf = inp->vel_map->jacob_vel_surf;
  up->vmap = inp->vel_map->vmap;
  up->pos_map = gkyl_vlasov_position_map_acquire(inp->pos_map);
  up->jacob_pos = inp->pos_map->jacob_pos;

  return up;
}
