/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_vlasov.h>    
#include <gkyl_dg_vlasov_priv.h>
}

#include <cassert>

// CUDA kernel to set device pointers to range object and vlasov kernel function
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
dg_vlasov_set_cu_dev_ptrs(struct dg_vlasov *vlasov, enum gkyl_basis_type b_type,
  int cdim, int vdim, int poly_order, enum gkyl_model_id model_id, enum gkyl_hamil_id hamil_id,
  bool has_E, bool has_phi, bool has_B, bool has_rad, bool use_lo)
{
  vlasov->eqn.vol_term = vlasov_vol;
  vlasov->eqn.surf_term = surf;
  vlasov->eqn.boundary_surf_term = boundary_surf;
  vlasov->eqn.boundary_diag_term = boundary_diag;

  // By default, we have no forces from E, phi, B, or radiation. 
  vlasov->E_vol = no_E_vol; 
  vlasov->phi_vol = no_phi_vol;
  vlasov->Bx_vol = no_B_vol; 
  vlasov->By_vol = no_B_vol; 
  vlasov->Bz_vol = no_B_vol;  
  vlasov->rad_vol = no_rad_vol; 

  const gkyl_dg_vlasov_stream_surf_kern_list *stream_surf_x_kernels, 
    *stream_surf_y_kernels, 
    *stream_surf_z_kernels;

  const gkyl_dg_vlasov_stream_boundary_surf_kern_list *stream_boundary_surf_x_kernels, 
    *stream_boundary_surf_y_kernels,
    *stream_boundary_surf_z_kernels;

  const gkyl_dg_vlasov_stream_surf_from_flux_kern_list *stream_surf_from_flux_x_kernels, 
    *stream_surf_from_flux_y_kernels, 
    *stream_surf_from_flux_z_kernels;

  const gkyl_dg_vlasov_stream_boundary_surf_from_flux_kern_list *stream_boundary_surf_from_flux_x_kernels, 
    *stream_boundary_surf_from_flux_y_kernels,
    *stream_boundary_surf_from_flux_z_kernels;

  const gkyl_dg_vlasov_accel_surf_kern_list *accel_surf_vx_kernels, 
    *accel_surf_vy_kernels, 
    *accel_surf_vz_kernels;

  const gkyl_dg_vlasov_accel_boundary_surf_kern_list *accel_boundary_surf_vx_kernels, 
    *accel_boundary_surf_vy_kernels,
    *accel_boundary_surf_vz_kernels;
  
  int kernel_index = cv_index[cdim].vdim[vdim]; 
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      // Set function pointers for individual pieces of the volume update.    
      if (model_id == GKYL_MODEL_DEFAULT || model_id == GKYL_MODEL_SR) {
        if (hamil_id == GKYL_HAMIL_VEL_SPARSE) {
          vlasov->hamil_vol = ser_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order];

          if (has_B) {
            vlasov->Bx_vol = ser_Bx_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order];
            vlasov->By_vol = ser_By_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order];
            vlasov->Bz_vol = ser_Bz_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order];
          }

          stream_surf_x_kernels = ser_stream_hamil_vel_sparse_surf_x_kernels;
          stream_surf_y_kernels = ser_stream_hamil_vel_sparse_surf_y_kernels;
          stream_surf_z_kernels = ser_stream_hamil_vel_sparse_surf_z_kernels;
          stream_boundary_surf_x_kernels = ser_stream_hamil_vel_sparse_boundary_surf_x_kernels;
          stream_boundary_surf_y_kernels = ser_stream_hamil_vel_sparse_boundary_surf_y_kernels;
          stream_boundary_surf_z_kernels = ser_stream_hamil_vel_sparse_boundary_surf_z_kernels;
        }
        else {
          vlasov->hamil_vol = ser_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];

          if (has_B) {
            vlasov->Bx_vol = ser_Bx_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];
            vlasov->By_vol = ser_By_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];
            vlasov->Bz_vol = ser_Bz_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];
          }

          stream_surf_x_kernels = ser_stream_hamil_vel_dense_surf_x_kernels;
          stream_surf_y_kernels = ser_stream_hamil_vel_dense_surf_y_kernels;
          stream_surf_z_kernels = ser_stream_hamil_vel_dense_surf_z_kernels;
          stream_boundary_surf_x_kernels = ser_stream_hamil_vel_dense_boundary_surf_x_kernels;
          stream_boundary_surf_y_kernels = ser_stream_hamil_vel_dense_boundary_surf_y_kernels;
          stream_boundary_surf_z_kernels = ser_stream_hamil_vel_dense_boundary_surf_z_kernels;
        }
        if (has_rad) vlasov->rad_vol = ser_rad_vol_kernels[kernel_index].kernels[poly_order];
      }
      else if (model_id == GKYL_MODEL_TRIAD || model_id == GKYL_MODEL_TRIAD_GR) {
        if (model_id == GKYL_MODEL_TRIAD) {
          vlasov->hamil_vol = (hamil_id == GKYL_HAMIL_VEL_SPARSE) ?
            ser_nc_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order] :
            ser_nc_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];
        }
        else if (model_id == GKYL_MODEL_TRIAD_GR) {
          // GR triads use the full phase-space Hamiltonian
          vlasov->hamil_vol = ser_nc_hamil_phase_vol_kernels[kernel_index].kernels[poly_order];
        }

        if ( use_lo ) {
          stream_surf_from_flux_x_kernels = ser_stream_surf_x_kernels;
          stream_surf_from_flux_y_kernels = ser_stream_surf_y_kernels;
          stream_surf_from_flux_z_kernels = ser_stream_surf_z_kernels;
          
          stream_boundary_surf_from_flux_x_kernels = ser_stream_boundary_surf_x_kernels;
          stream_boundary_surf_from_flux_y_kernels = ser_stream_boundary_surf_y_kernels;
          stream_boundary_surf_from_flux_z_kernels = ser_stream_boundary_surf_z_kernels;
        } 
        else {
          stream_surf_from_flux_x_kernels = ser_stream_ho_surf_x_kernels;
          stream_surf_from_flux_y_kernels = ser_stream_ho_surf_y_kernels;
          stream_surf_from_flux_z_kernels = ser_stream_ho_surf_z_kernels;
          
          stream_boundary_surf_from_flux_x_kernels = ser_stream_boundary_ho_surf_x_kernels;
          stream_boundary_surf_from_flux_y_kernels = ser_stream_boundary_ho_surf_y_kernels;
          stream_boundary_surf_from_flux_z_kernels = ser_stream_boundary_ho_surf_z_kernels;   
        }
      }
      else {
        // Canonical-PB models: volume term keeps the inline phase-Hamiltonian
        // kernels; streaming goes through the precomputed configuration-space
        // fluxes and the Hamiltonian-agnostic from-flux surface consumers,
        // exactly like the triad models.
        vlasov->hamil_vol = ser_hamil_phase_vol_kernels[kernel_index].kernels[poly_order];

        if ( use_lo ) {
          stream_surf_from_flux_x_kernels = ser_stream_surf_x_kernels;
          stream_surf_from_flux_y_kernels = ser_stream_surf_y_kernels;
          stream_surf_from_flux_z_kernels = ser_stream_surf_z_kernels;

          stream_boundary_surf_from_flux_x_kernels = ser_stream_boundary_surf_x_kernels;
          stream_boundary_surf_from_flux_y_kernels = ser_stream_boundary_surf_y_kernels;
          stream_boundary_surf_from_flux_z_kernels = ser_stream_boundary_surf_z_kernels;
        }
        else {
          stream_surf_from_flux_x_kernels = ser_stream_ho_surf_x_kernels;
          stream_surf_from_flux_y_kernels = ser_stream_ho_surf_y_kernels;
          stream_surf_from_flux_z_kernels = ser_stream_ho_surf_z_kernels;

          stream_boundary_surf_from_flux_x_kernels = ser_stream_boundary_ho_surf_x_kernels;
          stream_boundary_surf_from_flux_y_kernels = ser_stream_boundary_ho_surf_y_kernels;
          stream_boundary_surf_from_flux_z_kernels = ser_stream_boundary_ho_surf_z_kernels;
        }
      }
      if (has_E) vlasov->E_vol = ser_E_vol_kernels[kernel_index].kernels[poly_order];
      if (has_B) {
        if (model_id == GKYL_MODEL_TRIAD_GR) {
          vlasov->Bx_vol = ser_Bx_hamil_phase_vol_kernels[kernel_index].kernels[poly_order];
          vlasov->By_vol = ser_By_hamil_phase_vol_kernels[kernel_index].kernels[poly_order];
          vlasov->Bz_vol = ser_Bz_hamil_phase_vol_kernels[kernel_index].kernels[poly_order];
        }
        else {
          vlasov->Bx_vol = (hamil_id == GKYL_HAMIL_VEL_SPARSE) ?
            ser_Bx_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order] :
            ser_Bx_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];
          vlasov->By_vol = (hamil_id == GKYL_HAMIL_VEL_SPARSE) ?
            ser_By_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order] :
            ser_By_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];
          vlasov->Bz_vol = (hamil_id == GKYL_HAMIL_VEL_SPARSE) ?
            ser_Bz_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order] :
            ser_Bz_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];
        }
      }
      if (has_phi) vlasov->phi_vol = ser_phi_vol_kernels[kernel_index].kernels[poly_order];

      if (use_lo) {
        accel_surf_vx_kernels = ser_accel_surf_vx_kernels;
        accel_surf_vy_kernels = ser_accel_surf_vy_kernels;
        accel_surf_vz_kernels = ser_accel_surf_vz_kernels;

        accel_boundary_surf_vx_kernels = ser_accel_boundary_surf_vx_kernels;
        accel_boundary_surf_vy_kernels = ser_accel_boundary_surf_vy_kernels;
        accel_boundary_surf_vz_kernels = ser_accel_boundary_surf_vz_kernels;
      }
      else {
        accel_surf_vx_kernels = ser_accel_ho_surf_vx_kernels;
        accel_surf_vy_kernels = ser_accel_ho_surf_vy_kernels;
        accel_surf_vz_kernels = ser_accel_ho_surf_vz_kernels;

        accel_boundary_surf_vx_kernels = ser_accel_boundary_ho_surf_vx_kernels;
        accel_boundary_surf_vy_kernels = ser_accel_boundary_ho_surf_vy_kernels;
        accel_boundary_surf_vz_kernels = ser_accel_boundary_ho_surf_vz_kernels;
      }
      
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      // Set function pointers for individual pieces of the volume update. 
      if (model_id == GKYL_MODEL_DEFAULT || model_id == GKYL_MODEL_SR) {
        if (hamil_id == GKYL_HAMIL_VEL_SPARSE) {
          vlasov->hamil_vol = tensor_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order];

          if (has_B) {
            vlasov->Bx_vol = tensor_Bx_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order];
            vlasov->By_vol = tensor_By_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order];
            vlasov->Bz_vol = tensor_Bz_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order];
          }

          stream_surf_x_kernels = tensor_stream_hamil_vel_sparse_surf_x_kernels;
          stream_surf_y_kernels = tensor_stream_hamil_vel_sparse_surf_y_kernels;
          stream_surf_z_kernels = tensor_stream_hamil_vel_sparse_surf_z_kernels;

          stream_boundary_surf_x_kernels = tensor_stream_hamil_vel_sparse_boundary_surf_x_kernels;
          stream_boundary_surf_y_kernels = tensor_stream_hamil_vel_sparse_boundary_surf_y_kernels;
          stream_boundary_surf_z_kernels = tensor_stream_hamil_vel_sparse_boundary_surf_z_kernels;
        }
        else {
          vlasov->hamil_vol = tensor_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];

          if (has_B) {
            vlasov->Bx_vol = tensor_Bx_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];
            vlasov->By_vol = tensor_By_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];
            vlasov->Bz_vol = tensor_Bz_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];
          }

          stream_surf_x_kernels = tensor_stream_hamil_vel_dense_surf_x_kernels;
          stream_surf_y_kernels = tensor_stream_hamil_vel_dense_surf_y_kernels;
          stream_surf_z_kernels = tensor_stream_hamil_vel_dense_surf_z_kernels;

          stream_boundary_surf_x_kernels = tensor_stream_hamil_vel_dense_boundary_surf_x_kernels;
          stream_boundary_surf_y_kernels = tensor_stream_hamil_vel_dense_boundary_surf_y_kernels;
          stream_boundary_surf_z_kernels = tensor_stream_hamil_vel_dense_boundary_surf_z_kernels;
        }
        if (has_rad) vlasov->rad_vol = tensor_rad_vol_kernels[kernel_index].kernels[poly_order];
      }
      else if (model_id == GKYL_MODEL_TRIAD || model_id == GKYL_MODEL_TRIAD_GR) {
        if (model_id == GKYL_MODEL_TRIAD) {
          vlasov->hamil_vol = (hamil_id == GKYL_HAMIL_VEL_SPARSE) ?
            tensor_nc_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order] :
            tensor_nc_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];
        }
        else if (model_id == GKYL_MODEL_TRIAD_GR) {
          // GR triads use the full phase-space Hamiltonian
          vlasov->hamil_vol = tensor_nc_hamil_phase_vol_kernels[kernel_index].kernels[poly_order];
        }

        if ( use_lo ) {
          stream_surf_from_flux_x_kernels = tensor_stream_surf_x_kernels;
          stream_surf_from_flux_y_kernels = tensor_stream_surf_y_kernels;
          stream_surf_from_flux_z_kernels = tensor_stream_surf_z_kernels;

          stream_boundary_surf_from_flux_x_kernels = tensor_stream_boundary_surf_x_kernels;
          stream_boundary_surf_from_flux_y_kernels = tensor_stream_boundary_surf_y_kernels;
          stream_boundary_surf_from_flux_z_kernels = tensor_stream_boundary_surf_z_kernels;
        }
        else {
          stream_surf_from_flux_x_kernels = tensor_stream_ho_surf_x_kernels;
          stream_surf_from_flux_y_kernels = tensor_stream_ho_surf_y_kernels;
          stream_surf_from_flux_z_kernels = tensor_stream_ho_surf_z_kernels;

          stream_boundary_surf_from_flux_x_kernels = tensor_stream_boundary_ho_surf_x_kernels;
          stream_boundary_surf_from_flux_y_kernels = tensor_stream_boundary_ho_surf_y_kernels;
          stream_boundary_surf_from_flux_z_kernels = tensor_stream_boundary_ho_surf_z_kernels;
        }
      }
      else {
        // Canonical-PB models: only the p=1 tensor hybrid has a phase-space
        // Hamiltonian representation; volume keeps the inline phase-Hamiltonian
        // kernel and streaming goes through the from-flux surface consumers.
        vlasov->hamil_vol = tensor_hamil_phase_vol_kernels[kernel_index].kernels[poly_order];

        if ( use_lo ) {
          stream_surf_from_flux_x_kernels = tensor_stream_surf_x_kernels;
          stream_surf_from_flux_y_kernels = tensor_stream_surf_y_kernels;
          stream_surf_from_flux_z_kernels = tensor_stream_surf_z_kernels;

          stream_boundary_surf_from_flux_x_kernels = tensor_stream_boundary_surf_x_kernels;
          stream_boundary_surf_from_flux_y_kernels = tensor_stream_boundary_surf_y_kernels;
          stream_boundary_surf_from_flux_z_kernels = tensor_stream_boundary_surf_z_kernels;
        }
        else {
          stream_surf_from_flux_x_kernels = tensor_stream_ho_surf_x_kernels;
          stream_surf_from_flux_y_kernels = tensor_stream_ho_surf_y_kernels;
          stream_surf_from_flux_z_kernels = tensor_stream_ho_surf_z_kernels;

          stream_boundary_surf_from_flux_x_kernels = tensor_stream_boundary_ho_surf_x_kernels;
          stream_boundary_surf_from_flux_y_kernels = tensor_stream_boundary_ho_surf_y_kernels;
          stream_boundary_surf_from_flux_z_kernels = tensor_stream_boundary_ho_surf_z_kernels;
        }
      }
      if (has_E) vlasov->E_vol = tensor_E_vol_kernels[kernel_index].kernels[poly_order];
      if (has_B) {
        if (model_id == GKYL_MODEL_TRIAD_GR) {
          vlasov->Bx_vol = tensor_Bx_hamil_phase_vol_kernels[kernel_index].kernels[poly_order];
          vlasov->By_vol = tensor_By_hamil_phase_vol_kernels[kernel_index].kernels[poly_order];
          vlasov->Bz_vol = tensor_Bz_hamil_phase_vol_kernels[kernel_index].kernels[poly_order];
        }
        else {
          vlasov->Bx_vol = (hamil_id == GKYL_HAMIL_VEL_SPARSE) ?
            tensor_Bx_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order] :
            tensor_Bx_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];
          vlasov->By_vol = (hamil_id == GKYL_HAMIL_VEL_SPARSE) ?
            tensor_By_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order] :
            tensor_By_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];
          vlasov->Bz_vol = (hamil_id == GKYL_HAMIL_VEL_SPARSE) ?
            tensor_Bz_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order] :
            tensor_Bz_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];
        }
      }
      if (has_phi) vlasov->phi_vol = tensor_phi_vol_kernels[kernel_index].kernels[poly_order];

      // Only the tensor p=1 hybrid has distinct lo/ho surface variants; the
      // plain and ho lists share the (high-order by design) kernels at p>1.
      if ( use_lo ) {
        accel_surf_vx_kernels = tensor_accel_surf_vx_kernels;
        accel_surf_vy_kernels = tensor_accel_surf_vy_kernels;
        accel_surf_vz_kernels = tensor_accel_surf_vz_kernels;

        accel_boundary_surf_vx_kernels = tensor_accel_boundary_surf_vx_kernels;
        accel_boundary_surf_vy_kernels = tensor_accel_boundary_surf_vy_kernels;
        accel_boundary_surf_vz_kernels = tensor_accel_boundary_surf_vz_kernels;
      }
      else {
        accel_surf_vx_kernels = tensor_accel_ho_surf_vx_kernels;
        accel_surf_vy_kernels = tensor_accel_ho_surf_vy_kernels;
        accel_surf_vz_kernels = tensor_accel_ho_surf_vz_kernels;

        accel_boundary_surf_vx_kernels = tensor_accel_boundary_ho_surf_vx_kernels;
        accel_boundary_surf_vy_kernels = tensor_accel_boundary_ho_surf_vy_kernels;
        accel_boundary_surf_vz_kernels = tensor_accel_boundary_ho_surf_vz_kernels;
      }

      break;

    default:
      assert(false);
      break;    
  }

  if (model_id == GKYL_MODEL_TRIAD || hamil_id == GKYL_HAMIL_PHASE) {
    vlasov->stream_surf_from_flux[0] = stream_surf_from_flux_x_kernels[kernel_index].kernels[poly_order];
    vlasov->stream_surf_from_flux[1] = stream_surf_from_flux_y_kernels[kernel_index].kernels[poly_order];
    vlasov->stream_surf_from_flux[2] = stream_surf_from_flux_z_kernels[kernel_index].kernels[poly_order];

    vlasov->stream_boundary_surf_from_flux[0] = stream_boundary_surf_from_flux_x_kernels[kernel_index].kernels[poly_order];
    vlasov->stream_boundary_surf_from_flux[1] = stream_boundary_surf_from_flux_y_kernels[kernel_index].kernels[poly_order];
    vlasov->stream_boundary_surf_from_flux[2] = stream_boundary_surf_from_flux_z_kernels[kernel_index].kernels[poly_order]; 
  }
  else {
    vlasov->stream_surf[0] = stream_surf_x_kernels[kernel_index].kernels[poly_order];
    vlasov->stream_surf[1] = stream_surf_y_kernels[kernel_index].kernels[poly_order];
    vlasov->stream_surf[2] = stream_surf_z_kernels[kernel_index].kernels[poly_order];

    vlasov->stream_boundary_surf[0] = stream_boundary_surf_x_kernels[kernel_index].kernels[poly_order];
    vlasov->stream_boundary_surf[1] = stream_boundary_surf_y_kernels[kernel_index].kernels[poly_order];
    vlasov->stream_boundary_surf[2] = stream_boundary_surf_z_kernels[kernel_index].kernels[poly_order]; 
  }

  vlasov->accel_surf[0] = accel_surf_vx_kernels[kernel_index].kernels[poly_order];
  vlasov->accel_surf[1] = accel_surf_vy_kernels[kernel_index].kernels[poly_order];
  vlasov->accel_surf[2] = accel_surf_vz_kernels[kernel_index].kernels[poly_order];

  vlasov->accel_boundary_surf[0] = accel_boundary_surf_vx_kernels[kernel_index].kernels[poly_order];
  vlasov->accel_boundary_surf[1] = accel_boundary_surf_vy_kernels[kernel_index].kernels[poly_order];
  vlasov->accel_boundary_surf[2] = accel_boundary_surf_vz_kernels[kernel_index].kernels[poly_order];  
}

struct gkyl_dg_eqn*
gkyl_dg_vlasov_cu_dev_inew(const struct gkyl_dg_vlasov_inp *inp)
{
  struct dg_vlasov *vlasov = (struct dg_vlasov*) gkyl_malloc(sizeof(*vlasov));

  int cdim = inp->conf_basis->ndim, pdim = inp->phase_basis->ndim, vdim = pdim-cdim;
  int poly_order = inp->conf_basis->poly_order;

  vlasov->cdim = cdim;
  vlasov->pdim = pdim;  

  // Are we skipping cells with small phase space density?
  if (inp->skip_cell_thresh > 0.0) {
    vlasov->skip_cell_thresh = inp->skip_cell_thresh * pow(sqrt(2.0), pdim);
  }
  else {
    vlasov->skip_cell_thresh = -1.0;
  }

  // Determine Hamiltonian dimensionality and index offset for indexing Hamiltonian
  // from an input phase space index.
  if (inp->hamil_id == GKYL_HAMIL_PHASE) {
    vlasov->hamil_dim = pdim;
    vlasov->hamil_offset = 0;
  }
  else {
    vlasov->hamil_dim = vdim;
    vlasov->hamil_offset = cdim;
  }
  vlasov->hamil_range = *inp->hamil_range;
  vlasov->conf_range = *inp->conf_range;
  vlasov->phase_range = *inp->phase_range;

  // The velocity map is required: it provides the velocity-space range used
  // to index per-velocity-cell quantities (Jacobian, radiation drag).
  assert(inp->vel_map);
  vlasov->vel_map = 0; // Host pointer set after the device memcpy below.
  vlasov->vel_range = inp->vel_map->local_vel;
  vlasov->use_vmap = inp->vel_map->is_mapped;
  vlasov->jacob_vel = 0;
  vlasov->vmap = 0;
  if (vlasov->use_vmap) {
    // Unpack the raw device pointers for use inside kernels; lifetime is
    // guaranteed by acquiring the vel_map object on the host side below.
    vlasov->jacob_vel = inp->vel_map->jacob_vel->on_dev;
    vlasov->vmap = inp->vel_map->vmap->on_dev;
  }
  // Position map is required; unpack the raw device Jacobian pointer (lifetime
  // guaranteed by acquiring the pos_map object on the host side below).
  assert(inp->pos_map);
  vlasov->pos_map = 0; // Host pointer set after the device memcpy below.
  vlasov->jacob_pos = inp->pos_map->jacob_pos->on_dev;
  struct gkyl_array *poisson_tensor_conf_ho = gkyl_array_acquire(inp->poisson_tensor_conf);
  struct gkyl_array *hamil_ho = gkyl_array_acquire(inp->hamil); 
  // Optional force/drag inputs (absent e.g. for gyrokinetic neutral species).
  struct gkyl_array *qmem_ho = inp->qmem ? gkyl_array_acquire(inp->qmem) : 0; 
  struct gkyl_array *pot_tot_ho = inp->pot_tot ? gkyl_array_acquire(inp->pot_tot) : 0; 
  struct gkyl_array *rad_ho = inp->rad ? gkyl_array_acquire(inp->rad) : 0;
  vlasov->conf_flux_surf = 0;
  struct gkyl_array *conf_flux_surf_ho = 0;
  vlasov->use_conf_flux_surf = false;
  if (inp->model_id == GKYL_MODEL_TRIAD || inp->hamil_id == GKYL_HAMIL_PHASE) {
    vlasov->use_conf_flux_surf = true;
    conf_flux_surf_ho = gkyl_array_acquire(inp->conf_flux_surf);
    vlasov->conf_flux_surf = conf_flux_surf_ho->on_dev; 
  }
  struct gkyl_array *vel_flux_surf_ho = gkyl_array_acquire(inp->vel_flux_surf); 
  struct gkyl_array *f_no_J_ho = gkyl_array_acquire(inp->f_no_J); 
  // store pointers to on_dev for copying over to device. 
  vlasov->poisson_tensor_conf = poisson_tensor_conf_ho->on_dev; 
  vlasov->hamil = hamil_ho->on_dev; 
  vlasov->qmem = qmem_ho ? qmem_ho->on_dev : 0; 
  vlasov->pot_tot = pot_tot_ho ? pot_tot_ho->on_dev : 0; 
  vlasov->rad = rad_ho ? rad_ho->on_dev : 0; 
  vlasov->vel_flux_surf = vel_flux_surf_ho->on_dev; 
  vlasov->f_no_J =f_no_J_ho->on_dev; 
  vlasov->eqn.num_equations = 1;

  vlasov->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(vlasov->eqn.flags);
  vlasov->eqn.ref_count = gkyl_ref_count_init(gkyl_vlasov_free);

  // copy the host struct to device struct
  struct dg_vlasov *vlasov_cu = (struct dg_vlasov*) gkyl_cu_malloc(sizeof(struct dg_vlasov));
  gkyl_cu_memcpy(vlasov_cu, vlasov, sizeof(struct dg_vlasov), GKYL_CU_MEMCPY_H2D);

  dg_vlasov_set_cu_dev_ptrs<<<1,1>>>(vlasov_cu, inp->conf_basis->b_type, 
    cdim, vdim, poly_order, inp->model_id, inp->hamil_id,
    inp->has_E, inp->has_phi, inp->has_B, inp->has_rad, inp->use_lo);

  // set parent on_dev pointer
  vlasov->eqn.on_dev = &vlasov_cu->eqn;

  // Host-side equation object should store host pointers.
  vlasov->vel_map = gkyl_vlasov_velocity_map_acquire(inp->vel_map);
  vlasov->jacob_vel = vlasov->use_vmap ? inp->vel_map->jacob_vel : 0;
  vlasov->pos_map = gkyl_vlasov_position_map_acquire(inp->pos_map);
  vlasov->jacob_pos = inp->pos_map->jacob_pos;
  vlasov->poisson_tensor_conf = poisson_tensor_conf_ho;
  vlasov->hamil = hamil_ho; 
  vlasov->qmem = qmem_ho; 
  vlasov->pot_tot = pot_tot_ho; 
  vlasov->rad = rad_ho; 
  vlasov->conf_flux_surf = conf_flux_surf_ho; 
  vlasov->vel_flux_surf = vel_flux_surf_ho; 
  vlasov->f_no_J = f_no_J_ho; 

  return &vlasov->eqn;
}
