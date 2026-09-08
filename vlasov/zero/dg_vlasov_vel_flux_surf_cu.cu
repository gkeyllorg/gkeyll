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
gkyl_dg_vlasov_vel_flux_surf_advance_cu_kernel(struct gkyl_dg_vlasov_vel_flux_surf *up,
  struct gkyl_range conf_range, struct gkyl_range phase_range,
  const struct gkyl_array *vmap, const struct gkyl_array *jacob_pos, const struct gkyl_array *jacob_vel_surf, const struct gkyl_array *poisson_tensor_conf, const struct gkyl_array *hamil,
  const struct gkyl_array *qmem, const struct gkyl_array *pot_tot, const struct gkyl_array *rad,
  const struct gkyl_array *fin, struct gkyl_array *cflrate, struct gkyl_array *vel_flux_surf)
{
  // Per-node |alpha| values for the block, reduced to alpha_max per (cell, dir)
  // by the threadIdx.y == 0 thread of each cell.
  __shared__ double alpha_smem[GKYL_DEFAULT_NUM_THREADS];

  int pdim = up->pdim;
  int cdim = up->cdim;
  int vdim = pdim - cdim;

  // 2D thread grid: linc2 indexes the surface node (i_node: configuration-space
  // nodes, j_node: transverse velocity-space nodes, matching the CPU dispatch's
  // i-major/j-minor node order); threads in x index the phase-space cell.
  // blockDim.y == num_nodes_conf*num_nodes_vel by construction of the launch.
  int num_nodes_vel = up->num_nodes_vel;
  int num_nodes = up->num_nodes_conf*num_nodes_vel;
  int linc2 = threadIdx.y;
  int i_node = linc2 / num_nodes_vel;
  int j_node = linc2 % num_nodes_vel;

  // No grid-stride loop: the launch covers phase_range.volume exactly so the
  // __syncthreads() barriers below are reached uniformly by every thread in
  // the block (threads past the end of the range only skip the guarded work).
  unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
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

    for (int i=0; i<vdim; ++i) {
      idx_vel[i] = idx[cdim+i];
    }
    vidx = gkyl_range_idx(&up->vel_range, idx_vel);

    for (int i=0; i<up->hamil_dim; ++i) {
      idx_hamil[i] = idx[up->hamil_offset+i];
    }
    long hidx = gkyl_range_idx(&up->hamil_range, idx_hamil);

    // Grab the cell center location for NC bracket calculation
    gkyl_rect_grid_cell_center(&up->phase_grid, idx, xcC);

    poisson_tensor_conf_d = (const double*) gkyl_array_cfetch(poisson_tensor_conf, cidx);
    hamil_d = (const double*) gkyl_array_cfetch(hamil, hidx);
    qmem_d = qmem ? (const double*) gkyl_array_cfetch(qmem, cidx) : 0;
    pot_tot_d = pot_tot ? (const double*) gkyl_array_cfetch(pot_tot, cidx) : 0;
    rad_d = rad ? (const double*) gkyl_array_cfetch(rad, vidx) : 0;
    jacob_pos_d = jacob_pos ? (const double*) gkyl_array_cfetch(jacob_pos, cidx) : 0;
    vmap_d = (const double*) gkyl_array_cfetch(vmap, vidx);
    f_c = (const double*) gkyl_array_cfetch(fin, pidx);
    cflrate_d = (double*) gkyl_array_fetch(cflrate, pidx);
    flux = (double*) gkyl_array_fetch(vel_flux_surf, pidx);
  }

  // Each cell owns *lower* fluxes in each velocity-space direction.
  // So we need the distribution function in our current cell, and the cell
  // one index lower in each direction. If we are at the lower velocity-space
  // edge, we call special kernels because we have *no* ghost cells in velocity
  // space, so we cannot index the distribution function in the lower direction
  // along that velocity-space edge.
  for (int dir = 0; dir<vdim; ++dir) {
    bool edge = valid && (idx[cdim+dir] == phase_range.lower[cdim+dir]);
    const double *f_l = 0, *jacob_vel_l_d = 0, *jacob_vel_d = 0;
    double node_alpha = 0.0;
    if (valid && !edge) {
      gkyl_copy_int_arr(pdim, idx, idx_l);
      idx_l[cdim+dir] = idx_l[cdim+dir]-1;
      long pidx_l = gkyl_range_idx(&phase_range, idx_l);
      f_l = (const double*) gkyl_array_cfetch(fin, pidx_l);
      // Velocity-space Jacobian of the lower neighbor in this direction, for
      // the minimum-Jacobian time-step estimate of the C^0 linear map.
      int idx_vel_l[GKYL_MAX_DIM];
      for (int i=0; i<vdim; ++i) {
        idx_vel_l[i] = idx_l[cdim+i];
      }
      long vidx_l = gkyl_range_idx(&up->vel_range, idx_vel_l);
      jacob_vel_l_d = (const double*) gkyl_array_cfetch(jacob_vel_surf, vidx_l);
      jacob_vel_d = (const double*) gkyl_array_cfetch(jacob_vel_surf, vidx);

      // Compose the force producers into a register alpha at this thread's
      // surface node and apply the per-node Lax flux (disjoint nodal writes),
      // keeping |alpha| for the CFL reduction below.
      double alpha = up->hamil_alpha_quad[dir](i_node, j_node, xcC, up->phase_grid.dx, vmap_d, jacob_pos_d, jacob_vel_d, poisson_tensor_conf_d, hamil_d)
        + up->E_alpha_quad[dir](i_node, j_node, up->phase_grid.dx, qmem_d)
        + up->phi_alpha_quad[dir](i_node, j_node, up->phase_grid.dx, jacob_pos_d, pot_tot_d)
        + up->B_alpha_quad[dir](i_node, j_node, up->phase_grid.dx, jacob_vel_d, hamil_d, qmem_d)
        + up->rad_alpha_quad[dir](i_node, j_node, up->phase_grid.dx, rad_d);
      node_alpha = up->lax_flux_nodal[dir](i_node, j_node, jacob_vel_d, alpha, f_l, f_c, flux);
    }
    alpha_smem[threadIdx.x + blockDim.x*threadIdx.y] = node_alpha;
    __syncthreads();

    if (valid && threadIdx.y == 0) {
      if (edge) {
        cflrate_d[0] += up->vel_flux_surf_edge(up, dir, xcC, up->phase_grid.dx,
          jacob_vel_surf ? (const double*) gkyl_array_cfetch(jacob_vel_surf, vidx) : 0, poisson_tensor_conf_d,
          hamil_d, qmem_d, pot_tot_d, rad_d, f_c, flux);
      }
      else {
        // Reduce alpha_max in the CPU dispatch's node order so the fmax chain,
        // and hence the CFL estimate, matches the CPU loop exactly.
        double alpha_max = 0.0;
        for (int n=0; n<num_nodes; ++n) {
          alpha_max = fmax(alpha_max, alpha_smem[threadIdx.x + blockDim.x*n]);
        }
        double cfl = up->lax_cfl[dir](up->phase_grid.dx, jacob_vel_l_d, jacob_vel_d, alpha_max);
        // Always compute the flux, but if we are below threshold, ignore the stable time step estimate.
        if (fabs(f_l[0]) < up->skip_cell_thresh &&
            fabs(f_c[0]) < up->skip_cell_thresh) {
          cfl = 0.0;
        }
        cflrate_d[0] += cfl;
      }
    }
    // alpha_smem is reused by the next direction.
    __syncthreads();
  }
}

void
gkyl_dg_vlasov_vel_flux_surf_advance_cu(struct gkyl_dg_vlasov_vel_flux_surf *up,
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const struct gkyl_array *poisson_tensor_conf, const struct gkyl_array *hamil,
  const struct gkyl_array *qmem, const struct gkyl_array *pot_tot, const struct gkyl_array *rad,
  const struct gkyl_array *fin, struct gkyl_array *cflrate, struct gkyl_array *vel_flux_surf)
{
  // 2D thread grid: parallelize over surface nodes as well as phase-space cells.
  int num_nodes = up->num_nodes_conf*up->num_nodes_vel;
  assert(num_nodes <= GKYL_DEFAULT_NUM_THREADS);
  dim3 dimGrid, dimBlock;
  gkyl_parallelize_components_kernel_launch_dims(&dimGrid, &dimBlock, *phase_range, num_nodes);
  gkyl_dg_vlasov_vel_flux_surf_advance_cu_kernel<<<dimGrid, dimBlock>>>(up->on_dev,
    *conf_range, *phase_range, up->vmap->on_dev, up->jacob_pos->on_dev, up->jacob_vel_surf->on_dev, poisson_tensor_conf->on_dev,
    hamil->on_dev, qmem ? qmem->on_dev : 0, pot_tot ? pot_tot->on_dev : 0, rad ? rad->on_dev : 0,
    fin->on_dev, cflrate->on_dev, vel_flux_surf->on_dev);
}

// CUDA kernel to set device pointers to canonical pb vars kernel functions
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
gkyl_dg_vlasov_vel_flux_surf_set_cu_dev_ptrs(struct gkyl_dg_vlasov_vel_flux_surf *up,
  enum gkyl_basis_type b_type, int cdim, int vdim, int poly_order, 
  enum gkyl_model_id model_id, enum gkyl_hamil_id hamil_id, bool has_E, bool has_phi, bool has_B, bool has_rad, bool use_lo)
{
  // Sparse (separable) vs. dense velocity-space Hamiltonian kernel selection.
  bool hamil_sparse = (hamil_id == GKYL_HAMIL_VEL_SPARSE);
  // By default, we have no forces from Hamiltonian, E, B, or phi. 
  for (int d=0; d<vdim; ++d) {
    up->hamil_alpha_quad[d] = no_hamil_alpha_quad; 
    up->E_alpha_quad[d] = no_E_alpha_quad;
    up->phi_alpha_quad[d] = no_phi_alpha_quad; 
    up->B_alpha_quad[d] = no_B_alpha_quad;
    up->rad_alpha_quad[d] = no_rad_alpha_quad; 
  } 

  int kernel_index = cv_index[cdim].vdim[vdim];   
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      if ( use_lo ) {
        up->lax_flux_nodal[0] = ser_lax_flux_nodal_vx_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[0] = ser_lax_flux_nodal_vx_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[1] = ser_lax_flux_nodal_vy_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[1] = ser_lax_flux_nodal_vy_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[2] = ser_lax_flux_nodal_vz_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[2] = ser_lax_flux_nodal_vz_cfl_kernels[kernel_index].kernels[poly_order];
      }
      else {
        up->lax_flux_nodal[0] = ser_ho_lax_flux_nodal_vx_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[0] = ser_ho_lax_flux_nodal_vx_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[1] = ser_ho_lax_flux_nodal_vy_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[1] = ser_ho_lax_flux_nodal_vy_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[2] = ser_ho_lax_flux_nodal_vz_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[2] = ser_ho_lax_flux_nodal_vz_cfl_kernels[kernel_index].kernels[poly_order];
      }

      // Only have Hamiltonian forces in general geometry. 
      if (model_id == GKYL_MODEL_CANONICAL_PB || model_id == GKYL_MODEL_CANONICAL_PB_GR) {
        if ( use_lo ) {
          up->hamil_alpha_quad[0] = ser_hamil_phase_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[1] = ser_hamil_phase_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[2] = ser_hamil_phase_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }
        else {
          up->hamil_alpha_quad[0] = ser_hamil_phase_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[1] = ser_hamil_phase_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[2] = ser_hamil_phase_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }
      }
      else if (model_id == GKYL_MODEL_TRIAD) {
        if ( use_lo ) {
          up->hamil_alpha_quad[0] = hamil_sparse ?
            ser_nc_hamil_vel_sparse_alpha_quad_vx_kernels[kernel_index].kernels[poly_order] :
            ser_nc_hamil_vel_dense_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[1] = hamil_sparse ?
            ser_nc_hamil_vel_sparse_alpha_quad_vy_kernels[kernel_index].kernels[poly_order] :
            ser_nc_hamil_vel_dense_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[2] = hamil_sparse ?
            ser_nc_hamil_vel_sparse_alpha_quad_vz_kernels[kernel_index].kernels[poly_order] :
            ser_nc_hamil_vel_dense_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        } 
        else {
          up->hamil_alpha_quad[0] = hamil_sparse ?
            ser_nc_hamil_vel_sparse_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order] :
            ser_nc_hamil_vel_dense_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[1] = hamil_sparse ?
            ser_nc_hamil_vel_sparse_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order] :
            ser_nc_hamil_vel_dense_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[2] = hamil_sparse ?
            ser_nc_hamil_vel_sparse_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order] :
            ser_nc_hamil_vel_dense_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }
      }
      else if (model_id == GKYL_MODEL_TRIAD_GR) {
        if ( use_lo ) {
          up->hamil_alpha_quad[0] = ser_nc_hamil_phase_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[1] = ser_nc_hamil_phase_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[2] = ser_nc_hamil_phase_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        } 
        else {
          up->hamil_alpha_quad[0] = ser_nc_hamil_phase_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[1] = ser_nc_hamil_phase_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[2] = ser_nc_hamil_phase_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }
      }

      if ( use_lo ) {
        if (has_E) {
          up->E_alpha_quad[0] = ser_E_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->E_alpha_quad[1] = ser_E_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->E_alpha_quad[2] = ser_E_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }

        if (has_phi) {
          up->phi_alpha_quad[0] = ser_phi_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->phi_alpha_quad[1] = ser_phi_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->phi_alpha_quad[2] = ser_phi_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }

        if (has_B) {
          if (hamil_id == GKYL_HAMIL_PHASE) {
            up->B_alpha_quad[0] = ser_B_hamil_phase_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
            up->B_alpha_quad[1] = ser_B_hamil_phase_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
            up->B_alpha_quad[2] = ser_B_hamil_phase_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
          }
          else {
          up->B_alpha_quad[0] = hamil_sparse ?
            ser_B_hamil_vel_sparse_alpha_quad_vx_kernels[kernel_index].kernels[poly_order] :
            ser_B_hamil_vel_dense_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->B_alpha_quad[1] = hamil_sparse ?
            ser_B_hamil_vel_sparse_alpha_quad_vy_kernels[kernel_index].kernels[poly_order] :
            ser_B_hamil_vel_dense_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->B_alpha_quad[2] = hamil_sparse ?
            ser_B_hamil_vel_sparse_alpha_quad_vz_kernels[kernel_index].kernels[poly_order] :
            ser_B_hamil_vel_dense_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
          }
        }

        if (has_rad) {
          up->rad_alpha_quad[0] = ser_rad_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->rad_alpha_quad[1] = ser_rad_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->rad_alpha_quad[2] = ser_rad_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }    
      } 
      else {
        if (has_E) {
          up->E_alpha_quad[0] = ser_E_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->E_alpha_quad[1] = ser_E_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->E_alpha_quad[2] = ser_E_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }

        if (has_phi) {
          up->phi_alpha_quad[0] = ser_phi_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->phi_alpha_quad[1] = ser_phi_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->phi_alpha_quad[2] = ser_phi_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }

        if (has_B) {
          if (hamil_id == GKYL_HAMIL_PHASE) {
            up->B_alpha_quad[0] = ser_B_ho_hamil_phase_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
            up->B_alpha_quad[1] = ser_B_ho_hamil_phase_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
            up->B_alpha_quad[2] = ser_B_ho_hamil_phase_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
          }
          else {
          up->B_alpha_quad[0] = hamil_sparse ?
            ser_B_ho_hamil_vel_sparse_alpha_quad_vx_kernels[kernel_index].kernels[poly_order] :
            ser_B_ho_hamil_vel_dense_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->B_alpha_quad[1] = hamil_sparse ?
            ser_B_ho_hamil_vel_sparse_alpha_quad_vy_kernels[kernel_index].kernels[poly_order] :
            ser_B_ho_hamil_vel_dense_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->B_alpha_quad[2] = hamil_sparse ?
            ser_B_ho_hamil_vel_sparse_alpha_quad_vz_kernels[kernel_index].kernels[poly_order] :
            ser_B_ho_hamil_vel_dense_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
          }
        }

        if (has_rad) {
          up->rad_alpha_quad[0] = ser_rad_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->rad_alpha_quad[1] = ser_rad_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->rad_alpha_quad[2] = ser_rad_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }    
      }  

      break;

    case GKYL_BASIS_MODAL_TENSOR:
      // Only the tensor p=1 hybrid has distinct lo/ho surface variants; the
      // plain and ho lists share the (high-order by design) kernels at p>1.
      if ( use_lo ) {
        up->lax_flux_nodal[0] = tensor_lax_flux_nodal_vx_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[0] = tensor_lax_flux_nodal_vx_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[1] = tensor_lax_flux_nodal_vy_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[1] = tensor_lax_flux_nodal_vy_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[2] = tensor_lax_flux_nodal_vz_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[2] = tensor_lax_flux_nodal_vz_cfl_kernels[kernel_index].kernels[poly_order];
      }
      else {
        up->lax_flux_nodal[0] = tensor_ho_lax_flux_nodal_vx_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[0] = tensor_ho_lax_flux_nodal_vx_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[1] = tensor_ho_lax_flux_nodal_vy_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[1] = tensor_ho_lax_flux_nodal_vy_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[2] = tensor_ho_lax_flux_nodal_vz_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[2] = tensor_ho_lax_flux_nodal_vz_cfl_kernels[kernel_index].kernels[poly_order];
      }

      // Only have Hamiltonian forces in general geometry: the p=1 tensor
      // hybrid is the only tensor basis with a phase-space Hamiltonian
      // representation.
      if (model_id == GKYL_MODEL_CANONICAL_PB || model_id == GKYL_MODEL_CANONICAL_PB_GR) {
        if ( use_lo ) {
          up->hamil_alpha_quad[0] = tensor_hamil_phase_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[1] = tensor_hamil_phase_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[2] = tensor_hamil_phase_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }
        else {
          up->hamil_alpha_quad[0] = tensor_hamil_phase_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[1] = tensor_hamil_phase_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[2] = tensor_hamil_phase_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }
      }
      else if (model_id == GKYL_MODEL_TRIAD) {
        // Triad bracket on the tensor p=1 hybrid (sparse or dense velocity-space
        // Hamiltonian): per-node inverse velocity-map Jacobians of the C^1 cubic
        // map and the cubic vmap in the omega = v.Pi momentum factor.
        if ( use_lo ) {
          up->hamil_alpha_quad[0] = hamil_sparse ?
            tensor_nc_hamil_vel_sparse_alpha_quad_vx_kernels[kernel_index].kernels[poly_order] :
            tensor_nc_hamil_vel_dense_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[1] = hamil_sparse ?
            tensor_nc_hamil_vel_sparse_alpha_quad_vy_kernels[kernel_index].kernels[poly_order] :
            tensor_nc_hamil_vel_dense_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[2] = hamil_sparse ?
            tensor_nc_hamil_vel_sparse_alpha_quad_vz_kernels[kernel_index].kernels[poly_order] :
            tensor_nc_hamil_vel_dense_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }
        else {
          up->hamil_alpha_quad[0] = hamil_sparse ?
            tensor_nc_hamil_vel_sparse_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order] :
            tensor_nc_hamil_vel_dense_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[1] = hamil_sparse ?
            tensor_nc_hamil_vel_sparse_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order] :
            tensor_nc_hamil_vel_dense_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[2] = hamil_sparse ?
            tensor_nc_hamil_vel_sparse_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order] :
            tensor_nc_hamil_vel_dense_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }
      }
      else if (model_id == GKYL_MODEL_TRIAD_GR) {
        if ( use_lo ) {
          up->hamil_alpha_quad[0] = tensor_nc_hamil_phase_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[1] = tensor_nc_hamil_phase_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[2] = tensor_nc_hamil_phase_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }
        else {
          up->hamil_alpha_quad[0] = tensor_nc_hamil_phase_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[1] = tensor_nc_hamil_phase_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[2] = tensor_nc_hamil_phase_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }
      }

      if ( use_lo ) {
        if (has_E) {
          up->E_alpha_quad[0] = tensor_E_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->E_alpha_quad[1] = tensor_E_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->E_alpha_quad[2] = tensor_E_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }

        if (has_phi) {
          up->phi_alpha_quad[0] = tensor_phi_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->phi_alpha_quad[1] = tensor_phi_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->phi_alpha_quad[2] = tensor_phi_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }

        if (has_B) {
          // No phase-space Hamiltonian magnetic-force kernels for tensor p>1;
          // phase runs keep the no-op defaults set above. (The tensor p=1
          // hybrid does have B_hamil_phase kernels; hook them in with the
          // phase-Hamiltonian hybrid support.)
          if (hamil_id == GKYL_HAMIL_PHASE) {
            up->B_alpha_quad[0] = tensor_B_hamil_phase_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
            up->B_alpha_quad[1] = tensor_B_hamil_phase_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
            up->B_alpha_quad[2] = tensor_B_hamil_phase_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
          }
          else {
          up->B_alpha_quad[0] = hamil_sparse ?
            tensor_B_hamil_vel_sparse_alpha_quad_vx_kernels[kernel_index].kernels[poly_order] :
            tensor_B_hamil_vel_dense_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->B_alpha_quad[1] = hamil_sparse ?
            tensor_B_hamil_vel_sparse_alpha_quad_vy_kernels[kernel_index].kernels[poly_order] :
            tensor_B_hamil_vel_dense_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->B_alpha_quad[2] = hamil_sparse ?
            tensor_B_hamil_vel_sparse_alpha_quad_vz_kernels[kernel_index].kernels[poly_order] :
            tensor_B_hamil_vel_dense_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
          }
        }

        if (has_rad) {
          up->rad_alpha_quad[0] = tensor_rad_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->rad_alpha_quad[1] = tensor_rad_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->rad_alpha_quad[2] = tensor_rad_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }
      }
      else {
        if (has_E) {
          up->E_alpha_quad[0] = tensor_E_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->E_alpha_quad[1] = tensor_E_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->E_alpha_quad[2] = tensor_E_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }

        if (has_phi) {
          up->phi_alpha_quad[0] = tensor_phi_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->phi_alpha_quad[1] = tensor_phi_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->phi_alpha_quad[2] = tensor_phi_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }

        if (has_B) {
          // No phase-space Hamiltonian magnetic-force kernels for tensor p>1;
          // phase runs keep the no-op defaults set above. (The tensor p=1
          // hybrid does have B_hamil_phase kernels; hook them in with the
          // phase-Hamiltonian hybrid support.)
          if (hamil_id == GKYL_HAMIL_PHASE) {
            up->B_alpha_quad[0] = tensor_B_ho_hamil_phase_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
            up->B_alpha_quad[1] = tensor_B_ho_hamil_phase_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
            up->B_alpha_quad[2] = tensor_B_ho_hamil_phase_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
          }
          else {
          up->B_alpha_quad[0] = hamil_sparse ?
            tensor_B_ho_hamil_vel_sparse_alpha_quad_vx_kernels[kernel_index].kernels[poly_order] :
            tensor_B_ho_hamil_vel_dense_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->B_alpha_quad[1] = hamil_sparse ?
            tensor_B_ho_hamil_vel_sparse_alpha_quad_vy_kernels[kernel_index].kernels[poly_order] :
            tensor_B_ho_hamil_vel_dense_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->B_alpha_quad[2] = hamil_sparse ?
            tensor_B_ho_hamil_vel_sparse_alpha_quad_vz_kernels[kernel_index].kernels[poly_order] :
            tensor_B_ho_hamil_vel_dense_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
          }
        }

        if (has_rad) {
          up->rad_alpha_quad[0] = tensor_rad_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->rad_alpha_quad[1] = tensor_rad_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->rad_alpha_quad[2] = tensor_rad_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        }
      }

      break;

    default:
      assert(false);
      break;    
  } 
  // Set assembly functions for computing fluxes. 
  up->vel_flux_surf = vel_flux_surf_nodes;
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
  for (int d=0; d<cdim; ++d) up->num_nodes_conf *= nq_conf;
  up->num_nodes_vel = 1;
  for (int d=0; d<vdim-1; ++d) up->num_nodes_vel *= nq_vel;
  // Currently only support zero-flux boundary, so edge velocity flux an empty function (flux = 0.0). 
  up->vel_flux_surf_edge = no_vel_flux_surf_edge;   
}

gkyl_dg_vlasov_vel_flux_surf*
gkyl_dg_vlasov_vel_flux_surf_cu_dev_inew(const struct gkyl_dg_vlasov_vel_flux_surf_inp *inp)
{
  struct gkyl_dg_vlasov_vel_flux_surf *up = (struct gkyl_dg_vlasov_vel_flux_surf*) gkyl_malloc(sizeof(*up));

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
  if (inp->hamil_id == GKYL_HAMIL_PHASE) {
    up->hamil_dim = pdim;
    up->hamil_offset = 0;
  }
  else {
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

  // Host mirror of the surface node counts set on the device struct by
  // set_cu_dev_ptrs; the advance wrapper needs them to size the 2D
  // (cells x nodes) kernel launch. Keep in sync with set_cu_dev_ptrs: the
  // tensor p=1 hybrid is anisotropic (2 conf nodes, 3/4 vel nodes by use_lo).
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
  for (int d=0; d<cdim; ++d) up->num_nodes_conf *= nq_conf;
  up->num_nodes_vel = 1;
  for (int d=0; d<vdim-1; ++d) up->num_nodes_vel *= nq_vel;

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_dg_vlasov_vel_flux_surf *up_cu = (struct gkyl_dg_vlasov_vel_flux_surf*) gkyl_cu_malloc(sizeof(*up_cu));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_dg_vlasov_vel_flux_surf), GKYL_CU_MEMCPY_H2D);

  gkyl_dg_vlasov_vel_flux_surf_set_cu_dev_ptrs<<<1,1>>>(up_cu, inp->conf_basis->b_type,
    cdim, vdim, poly_order, inp->model_id, inp->hamil_id, inp->has_E, inp->has_phi, inp->has_B, inp->has_rad, inp->use_lo);

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
