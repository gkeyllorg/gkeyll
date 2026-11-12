#include "gkyl_dg_eqn.h"
#include <assert.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_dg_vlasov.h>
#include <gkyl_dg_vlasov_priv.h>
#include <gkyl_util.h>

void
gkyl_vlasov_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_eqn *base = container_of(ref, struct gkyl_dg_eqn, ref_count);

  if (gkyl_dg_eqn_is_cu_dev(base)) {
    // free inner on_dev object
    struct dg_vlasov *vlasov = container_of(base->on_dev, struct dg_vlasov, eqn);
    gkyl_cu_free(vlasov);
  }

  struct dg_vlasov *vlasov = container_of(base, struct dg_vlasov, eqn);
  gkyl_array_release(vlasov->poisson_tensor_conf);
  gkyl_array_release(vlasov->hamil);
  if (vlasov->qmem) {
    gkyl_array_release(vlasov->qmem);
  }
  if (vlasov->pot_tot) {
    gkyl_array_release(vlasov->pot_tot);
  }
  if (vlasov->rad) {
    gkyl_array_release(vlasov->rad);
  }
  if (vlasov->vel_map) {
    gkyl_vlasov_velocity_map_release(vlasov->vel_map);
  }
  if (vlasov->pos_map) {
    gkyl_vlasov_position_map_release(vlasov->pos_map);
  }
  if (vlasov->use_conf_flux_surf) {
    gkyl_array_release(vlasov->conf_flux_surf);
  }
  gkyl_array_release(vlasov->vel_flux_surf);
  gkyl_array_release(vlasov->f_no_J);
  gkyl_free(vlasov);
}

// Host-side check that every volume kernel the requested combination of basis,
// model and Hamiltonian representation needs actually exists in the tables
// (some rows are NULL, e.g. the 3x3v phase-space B kernels). Runs before the
// device constructor too, so an unsupported combination fails with an assert
// instead of a NULL device function pointer.
static void
dg_vlasov_check_vol_kernels(const struct gkyl_dg_vlasov_inp *inp)
{
  int cdim = inp->conf_basis->ndim, pdim = inp->phase_basis->ndim, vdim = pdim - cdim;
  int po = inp->conf_basis->poly_order;
  int ki = cv_index[cdim].vdim[vdim];
  assert(ki != -1);
  bool ten =
    (gkyl_basis_phase_kernel_type(inp->conf_basis, inp->phase_basis) == GKYL_BASIS_MODAL_TENSOR);
  bool sparse = (inp->hamil_id == GKYL_HAMIL_VEL_SPARSE);
  bool phase = (inp->hamil_id == GKYL_HAMIL_PHASE);
  bool triad = (inp->model_id == GKYL_MODEL_TRIAD || inp->model_id == GKYL_MODEL_TRIAD_GR);

  const gkyl_dg_vlasov_hamil_vol_kern_list *hl;
  if (triad) {
    hl = phase ?
           (ten ? tensor_nc_hamil_phase_vol_kernels : ser_nc_hamil_phase_vol_kernels) :
         sparse ?
           (ten ? tensor_nc_hamil_vel_sparse_vol_kernels : ser_nc_hamil_vel_sparse_vol_kernels) :
           (ten ? tensor_nc_hamil_vel_dense_vol_kernels : ser_nc_hamil_vel_dense_vol_kernels);
  } else {
    hl = phase  ? (ten ? tensor_hamil_phase_vol_kernels : ser_hamil_phase_vol_kernels) :
         sparse ? (ten ? tensor_hamil_vel_sparse_vol_kernels : ser_hamil_vel_sparse_vol_kernels) :
                  (ten ? tensor_hamil_vel_dense_vol_kernels : ser_hamil_vel_dense_vol_kernels);
  }
  assert(hl[ki].kernels[po]);

  if (inp->has_B) {
    const gkyl_dg_vlasov_B_vol_kern_list *bx, *by, *bz;
    if (phase) {
      bx = ten ? tensor_Bx_hamil_phase_vol_kernels : ser_Bx_hamil_phase_vol_kernels;
      by = ten ? tensor_By_hamil_phase_vol_kernels : ser_By_hamil_phase_vol_kernels;
      bz = ten ? tensor_Bz_hamil_phase_vol_kernels : ser_Bz_hamil_phase_vol_kernels;
    } else if (sparse) {
      bx = ten ? tensor_Bx_hamil_vel_sparse_vol_kernels : ser_Bx_hamil_vel_sparse_vol_kernels;
      by = ten ? tensor_By_hamil_vel_sparse_vol_kernels : ser_By_hamil_vel_sparse_vol_kernels;
      bz = ten ? tensor_Bz_hamil_vel_sparse_vol_kernels : ser_Bz_hamil_vel_sparse_vol_kernels;
    } else {
      bx = ten ? tensor_Bx_hamil_vel_dense_vol_kernels : ser_Bx_hamil_vel_dense_vol_kernels;
      by = ten ? tensor_By_hamil_vel_dense_vol_kernels : ser_By_hamil_vel_dense_vol_kernels;
      bz = ten ? tensor_Bz_hamil_vel_dense_vol_kernels : ser_Bz_hamil_vel_dense_vol_kernels;
    }
    assert(bx[ki].kernels[po]);
    assert(by[ki].kernels[po]);
    assert(bz[ki].kernels[po]);
  }
  if (inp->has_E) {
    assert((ten ? tensor_E_vol_kernels : ser_E_vol_kernels)[ki].kernels[po]);
  }
  if (inp->has_phi) {
    assert((ten ? tensor_phi_vol_kernels : ser_phi_vol_kernels)[ki].kernels[po]);
  }
  if (inp->has_rad) {
    assert((ten ? tensor_rad_vol_kernels : ser_rad_vol_kernels)[ki].kernels[po]);
  }
}

struct gkyl_dg_eqn *
gkyl_dg_vlasov_inew(const struct gkyl_dg_vlasov_inp *inp)
{
  dg_vlasov_check_vol_kernels(inp);

#ifdef GKYL_HAVE_CUDA
  if (inp->use_gpu) {
    return gkyl_dg_vlasov_cu_dev_inew(inp);
  }
#endif
  struct dg_vlasov *vlasov = gkyl_malloc(sizeof(*vlasov));

  int cdim = inp->conf_basis->ndim, pdim = inp->phase_basis->ndim, vdim = pdim - cdim;
  int poly_order = inp->conf_basis->poly_order;

  vlasov->cdim = cdim;
  vlasov->pdim = pdim;

  // Are we skipping cells with small phase space density?
  if (inp->skip_cell_thresh > 0.0) {
    vlasov->skip_cell_thresh = inp->skip_cell_thresh * pow(sqrt(2.0), pdim);
  } else {
    vlasov->skip_cell_thresh = -1.0;
  }

  // Determine Hamiltonian dimensionality and index offset for indexing Hamiltonian
  // from an input phase space index.
  if (inp->hamil_id == GKYL_HAMIL_PHASE) {
    vlasov->hamil_dim = pdim;
    vlasov->hamil_offset = 0;
  } else {
    vlasov->hamil_dim = vdim;
    vlasov->hamil_offset = cdim;
  }
  vlasov->hamil_range = *inp->hamil_range;
  vlasov->conf_range = *inp->conf_range;
  vlasov->phase_range = *inp->phase_range;

  // The velocity map is required: it provides the velocity-space range used
  // to index per-velocity-cell quantities (Jacobian, radiation drag).
  assert(inp->vel_map);
  vlasov->vel_map = gkyl_vlasov_velocity_map_acquire(inp->vel_map);
  vlasov->vel_range = inp->vel_map->local_vel;
  vlasov->use_vmap = inp->vel_map->is_mapped;
  vlasov->jacob_vel = 0;
  vlasov->vmap = 0;
  if (vlasov->use_vmap) {
    // Borrowed pointers; kept alive by the acquired vel_map.
    vlasov->jacob_vel = inp->vel_map->jacob_vel;
    vlasov->vmap = inp->vel_map->vmap;
  }
  // The position map is required: it provides the (per-conf-cell constant)
  // configuration-space Jacobian used to transform the streaming term.
  assert(inp->pos_map);
  vlasov->pos_map = gkyl_vlasov_position_map_acquire(inp->pos_map);
  vlasov->jacob_pos = inp->pos_map->jacob_pos; // Borrowed; kept alive by the acquired pos_map.
  vlasov->poisson_tensor_conf = gkyl_array_acquire(inp->poisson_tensor_conf);
  vlasov->hamil = gkyl_array_acquire(inp->hamil);
  vlasov->qmem = inp->qmem ? gkyl_array_acquire(inp->qmem) : 0;
  vlasov->pot_tot = inp->pot_tot ? gkyl_array_acquire(inp->pot_tot) : 0;
  vlasov->rad = inp->rad ? gkyl_array_acquire(inp->rad) : 0;
  vlasov->use_conf_flux_surf = false;
  if (inp->model_id == GKYL_MODEL_TRIAD || inp->hamil_id == GKYL_HAMIL_PHASE) {
    vlasov->use_conf_flux_surf = true;
    vlasov->conf_flux_surf = gkyl_array_acquire(inp->conf_flux_surf);
  }
  vlasov->vel_flux_surf = gkyl_array_acquire(inp->vel_flux_surf);
  vlasov->f_no_J = gkyl_array_acquire(inp->f_no_J);

  vlasov->eqn.num_equations = 1;
  vlasov->eqn.vol_term = vlasov_vol;
  vlasov->eqn.surf_term = surf;
  vlasov->eqn.boundary_surf_term = boundary_surf;
  vlasov->eqn.boundary_diag_term = boundary_diag;

  // By default, we have no forces from E, B, phi, or radiation.
  vlasov->E_vol = no_E_vol;
  vlasov->Bx_vol = no_B_vol;
  vlasov->By_vol = no_B_vol;
  vlasov->Bz_vol = no_B_vol;
  vlasov->phi_vol = no_phi_vol;
  vlasov->rad_vol = no_rad_vol;

  const gkyl_dg_vlasov_stream_surf_kern_list *stream_surf_x_kernels, *stream_surf_y_kernels,
    *stream_surf_z_kernels;

  const gkyl_dg_vlasov_stream_boundary_surf_kern_list *stream_boundary_surf_x_kernels,
    *stream_boundary_surf_y_kernels, *stream_boundary_surf_z_kernels;

  const gkyl_dg_vlasov_stream_surf_from_flux_kern_list *stream_surf_from_flux_x_kernels,
    *stream_surf_from_flux_y_kernels, *stream_surf_from_flux_z_kernels;

  const gkyl_dg_vlasov_stream_boundary_surf_from_flux_kern_list
    *stream_boundary_surf_from_flux_x_kernels,
    *stream_boundary_surf_from_flux_y_kernels, *stream_boundary_surf_from_flux_z_kernels;

  const gkyl_dg_vlasov_accel_surf_kern_list *accel_surf_vx_kernels, *accel_surf_vy_kernels,
    *accel_surf_vz_kernels;

  const gkyl_dg_vlasov_accel_boundary_surf_kern_list *accel_boundary_surf_vx_kernels,
    *accel_boundary_surf_vy_kernels, *accel_boundary_surf_vz_kernels;

  int kernel_index = cv_index[cdim].vdim[vdim];
  switch (gkyl_basis_phase_kernel_type(inp->conf_basis, inp->phase_basis)) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      // Set function pointers for individual pieces of the volume update.
      if (inp->model_id == GKYL_MODEL_DEFAULT || inp->model_id == GKYL_MODEL_SR) {
        if (inp->hamil_id == GKYL_HAMIL_VEL_SPARSE) {
          vlasov->hamil_vol = ser_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order];

          if (inp->has_B) {
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
        } else {
          vlasov->hamil_vol = ser_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];

          if (inp->has_B) {
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
        if (inp->has_rad) {
          vlasov->rad_vol = ser_rad_vol_kernels[kernel_index].kernels[poly_order];
        }

      } else if (inp->model_id == GKYL_MODEL_TRIAD || inp->model_id == GKYL_MODEL_TRIAD_GR) {
        if (inp->model_id == GKYL_MODEL_TRIAD) {
          vlasov->hamil_vol =
            (inp->hamil_id == GKYL_HAMIL_VEL_SPARSE) ?
              ser_nc_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order] :
              ser_nc_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];
        } else if (inp->model_id == GKYL_MODEL_TRIAD_GR) {
          vlasov->hamil_vol = ser_nc_hamil_phase_vol_kernels[kernel_index].kernels[poly_order];
        }

        if (inp->use_lo) {
          stream_surf_from_flux_x_kernels = ser_stream_surf_x_kernels;
          stream_surf_from_flux_y_kernels = ser_stream_surf_y_kernels;
          stream_surf_from_flux_z_kernels = ser_stream_surf_z_kernels;

          stream_boundary_surf_from_flux_x_kernels = ser_stream_boundary_surf_x_kernels;
          stream_boundary_surf_from_flux_y_kernels = ser_stream_boundary_surf_y_kernels;
          stream_boundary_surf_from_flux_z_kernels = ser_stream_boundary_surf_z_kernels;
        } else {
          stream_surf_from_flux_x_kernels = ser_stream_ho_surf_x_kernels;
          stream_surf_from_flux_y_kernels = ser_stream_ho_surf_y_kernels;
          stream_surf_from_flux_z_kernels = ser_stream_ho_surf_z_kernels;

          stream_boundary_surf_from_flux_x_kernels = ser_stream_boundary_ho_surf_x_kernels;
          stream_boundary_surf_from_flux_y_kernels = ser_stream_boundary_ho_surf_y_kernels;
          stream_boundary_surf_from_flux_z_kernels = ser_stream_boundary_ho_surf_z_kernels;
        }
      } else {
        // Canonical-PB models: volume term keeps the inline phase-Hamiltonian
        // kernels; streaming goes through the precomputed configuration-space
        // fluxes (conf_flux_surf updater with the identity Poisson tensor) and
        // the Hamiltonian-agnostic from-flux surface consumers, exactly like
        // the triad models.
        vlasov->hamil_vol = ser_hamil_phase_vol_kernels[kernel_index].kernels[poly_order];

        if (inp->use_lo) {
          stream_surf_from_flux_x_kernels = ser_stream_surf_x_kernels;
          stream_surf_from_flux_y_kernels = ser_stream_surf_y_kernels;
          stream_surf_from_flux_z_kernels = ser_stream_surf_z_kernels;

          stream_boundary_surf_from_flux_x_kernels = ser_stream_boundary_surf_x_kernels;
          stream_boundary_surf_from_flux_y_kernels = ser_stream_boundary_surf_y_kernels;
          stream_boundary_surf_from_flux_z_kernels = ser_stream_boundary_surf_z_kernels;
        } else {
          stream_surf_from_flux_x_kernels = ser_stream_ho_surf_x_kernels;
          stream_surf_from_flux_y_kernels = ser_stream_ho_surf_y_kernels;
          stream_surf_from_flux_z_kernels = ser_stream_ho_surf_z_kernels;

          stream_boundary_surf_from_flux_x_kernels = ser_stream_boundary_ho_surf_x_kernels;
          stream_boundary_surf_from_flux_y_kernels = ser_stream_boundary_ho_surf_y_kernels;
          stream_boundary_surf_from_flux_z_kernels = ser_stream_boundary_ho_surf_z_kernels;
        }
      }
      if (inp->has_E) {
        vlasov->E_vol = ser_E_vol_kernels[kernel_index].kernels[poly_order];
      }
      if (inp->has_B) {
        // Phase-space Hamiltonians (triad-GR and canonical-PB models) use the
        // phase-B kernels, matching the velocity-flux updater which also keys on
        // hamil_id; velocity-space Hamiltonians use the sparse/dense families.
        if (inp->hamil_id == GKYL_HAMIL_PHASE) {
          vlasov->Bx_vol = ser_Bx_hamil_phase_vol_kernels[kernel_index].kernels[poly_order];
          vlasov->By_vol = ser_By_hamil_phase_vol_kernels[kernel_index].kernels[poly_order];
          vlasov->Bz_vol = ser_Bz_hamil_phase_vol_kernels[kernel_index].kernels[poly_order];
        } else {
          vlasov->Bx_vol = (inp->hamil_id == GKYL_HAMIL_VEL_SPARSE) ?
                             ser_Bx_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order] :
                             ser_Bx_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];
          vlasov->By_vol = (inp->hamil_id == GKYL_HAMIL_VEL_SPARSE) ?
                             ser_By_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order] :
                             ser_By_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];
          vlasov->Bz_vol = (inp->hamil_id == GKYL_HAMIL_VEL_SPARSE) ?
                             ser_Bz_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order] :
                             ser_Bz_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];
        }
      }
      if (inp->has_phi) {
        vlasov->phi_vol = ser_phi_vol_kernels[kernel_index].kernels[poly_order];
      }

      if (inp->use_lo) {
        accel_surf_vx_kernels = ser_accel_surf_vx_kernels;
        accel_surf_vy_kernels = ser_accel_surf_vy_kernels;
        accel_surf_vz_kernels = ser_accel_surf_vz_kernels;
        accel_boundary_surf_vx_kernels = ser_accel_boundary_surf_vx_kernels;
        accel_boundary_surf_vy_kernels = ser_accel_boundary_surf_vy_kernels;
        accel_boundary_surf_vz_kernels = ser_accel_boundary_surf_vz_kernels;
      } else {
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
      if (inp->model_id == GKYL_MODEL_DEFAULT || inp->model_id == GKYL_MODEL_SR) {
        if (inp->hamil_id == GKYL_HAMIL_VEL_SPARSE) {
          vlasov->hamil_vol = tensor_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order];

          if (inp->has_B) {
            vlasov->Bx_vol =
              tensor_Bx_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order];
            vlasov->By_vol =
              tensor_By_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order];
            vlasov->Bz_vol =
              tensor_Bz_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order];
          }

          stream_surf_x_kernels = tensor_stream_hamil_vel_sparse_surf_x_kernels;
          stream_surf_y_kernels = tensor_stream_hamil_vel_sparse_surf_y_kernels;
          stream_surf_z_kernels = tensor_stream_hamil_vel_sparse_surf_z_kernels;

          stream_boundary_surf_x_kernels = tensor_stream_hamil_vel_sparse_boundary_surf_x_kernels;
          stream_boundary_surf_y_kernels = tensor_stream_hamil_vel_sparse_boundary_surf_y_kernels;
          stream_boundary_surf_z_kernels = tensor_stream_hamil_vel_sparse_boundary_surf_z_kernels;
        } else {
          vlasov->hamil_vol = tensor_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];

          if (inp->has_B) {
            vlasov->Bx_vol =
              tensor_Bx_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];
            vlasov->By_vol =
              tensor_By_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];
            vlasov->Bz_vol =
              tensor_Bz_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];
          }

          stream_surf_x_kernels = tensor_stream_hamil_vel_dense_surf_x_kernels;
          stream_surf_y_kernels = tensor_stream_hamil_vel_dense_surf_y_kernels;
          stream_surf_z_kernels = tensor_stream_hamil_vel_dense_surf_z_kernels;

          stream_boundary_surf_x_kernels = tensor_stream_hamil_vel_dense_boundary_surf_x_kernels;
          stream_boundary_surf_y_kernels = tensor_stream_hamil_vel_dense_boundary_surf_y_kernels;
          stream_boundary_surf_z_kernels = tensor_stream_hamil_vel_dense_boundary_surf_z_kernels;
        }
        if (inp->has_rad) {
          vlasov->rad_vol = tensor_rad_vol_kernels[kernel_index].kernels[poly_order];
        }
      } else if (inp->model_id == GKYL_MODEL_TRIAD || inp->model_id == GKYL_MODEL_TRIAD_GR) {
        if (inp->model_id == GKYL_MODEL_TRIAD) {
          vlasov->hamil_vol =
            (inp->hamil_id == GKYL_HAMIL_VEL_SPARSE) ?
              tensor_nc_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order] :
              tensor_nc_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];
        } else if (inp->model_id == GKYL_MODEL_TRIAD_GR) {
          vlasov->hamil_vol = tensor_nc_hamil_phase_vol_kernels[kernel_index].kernels[poly_order];
        }

        if (inp->use_lo) {
          stream_surf_from_flux_x_kernels = tensor_stream_surf_x_kernels;
          stream_surf_from_flux_y_kernels = tensor_stream_surf_y_kernels;
          stream_surf_from_flux_z_kernels = tensor_stream_surf_z_kernels;

          stream_boundary_surf_from_flux_x_kernels = tensor_stream_boundary_surf_x_kernels;
          stream_boundary_surf_from_flux_y_kernels = tensor_stream_boundary_surf_y_kernels;
          stream_boundary_surf_from_flux_z_kernels = tensor_stream_boundary_surf_z_kernels;
        } else {
          stream_surf_from_flux_x_kernels = tensor_stream_ho_surf_x_kernels;
          stream_surf_from_flux_y_kernels = tensor_stream_ho_surf_y_kernels;
          stream_surf_from_flux_z_kernels = tensor_stream_ho_surf_z_kernels;

          stream_boundary_surf_from_flux_x_kernels = tensor_stream_boundary_ho_surf_x_kernels;
          stream_boundary_surf_from_flux_y_kernels = tensor_stream_boundary_ho_surf_y_kernels;
          stream_boundary_surf_from_flux_z_kernels = tensor_stream_boundary_ho_surf_z_kernels;
        }
      } else {
        // Canonical-PB models: only the p=1 tensor hybrid has a phase-space
        // Hamiltonian representation. The volume term keeps the inline
        // phase-Hamiltonian kernel; streaming goes through the precomputed
        // configuration-space fluxes and the Hamiltonian-agnostic from-flux
        // surface consumers, exactly like the Serendipity path.
        vlasov->hamil_vol = tensor_hamil_phase_vol_kernels[kernel_index].kernels[poly_order];

        if (inp->use_lo) {
          stream_surf_from_flux_x_kernels = tensor_stream_surf_x_kernels;
          stream_surf_from_flux_y_kernels = tensor_stream_surf_y_kernels;
          stream_surf_from_flux_z_kernels = tensor_stream_surf_z_kernels;

          stream_boundary_surf_from_flux_x_kernels = tensor_stream_boundary_surf_x_kernels;
          stream_boundary_surf_from_flux_y_kernels = tensor_stream_boundary_surf_y_kernels;
          stream_boundary_surf_from_flux_z_kernels = tensor_stream_boundary_surf_z_kernels;
        } else {
          stream_surf_from_flux_x_kernels = tensor_stream_ho_surf_x_kernels;
          stream_surf_from_flux_y_kernels = tensor_stream_ho_surf_y_kernels;
          stream_surf_from_flux_z_kernels = tensor_stream_ho_surf_z_kernels;

          stream_boundary_surf_from_flux_x_kernels = tensor_stream_boundary_ho_surf_x_kernels;
          stream_boundary_surf_from_flux_y_kernels = tensor_stream_boundary_ho_surf_y_kernels;
          stream_boundary_surf_from_flux_z_kernels = tensor_stream_boundary_ho_surf_z_kernels;
        }
      }
      if (inp->has_E) {
        vlasov->E_vol = tensor_E_vol_kernels[kernel_index].kernels[poly_order];
      }
      if (inp->has_B) {
        // Phase-space Hamiltonians (triad-GR and canonical-PB models) use the
        // phase-B kernels, matching the velocity-flux updater which also keys on
        // hamil_id; velocity-space Hamiltonians use the sparse/dense families.
        if (inp->hamil_id == GKYL_HAMIL_PHASE) {
          vlasov->Bx_vol = tensor_Bx_hamil_phase_vol_kernels[kernel_index].kernels[poly_order];
          vlasov->By_vol = tensor_By_hamil_phase_vol_kernels[kernel_index].kernels[poly_order];
          vlasov->Bz_vol = tensor_Bz_hamil_phase_vol_kernels[kernel_index].kernels[poly_order];
        } else {
          vlasov->Bx_vol =
            (inp->hamil_id == GKYL_HAMIL_VEL_SPARSE) ?
              tensor_Bx_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order] :
              tensor_Bx_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];
          vlasov->By_vol =
            (inp->hamil_id == GKYL_HAMIL_VEL_SPARSE) ?
              tensor_By_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order] :
              tensor_By_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];
          vlasov->Bz_vol =
            (inp->hamil_id == GKYL_HAMIL_VEL_SPARSE) ?
              tensor_Bz_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order] :
              tensor_Bz_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];
        }
      }
      if (inp->has_phi) {
        vlasov->phi_vol = tensor_phi_vol_kernels[kernel_index].kernels[poly_order];
      }

      // Only the tensor p=1 hybrid has distinct lo/ho surface variants; the
      // plain and ho lists share the (high-order by design) kernels at p>1.
      if (inp->use_lo) {
        accel_surf_vx_kernels = tensor_accel_surf_vx_kernels;
        accel_surf_vy_kernels = tensor_accel_surf_vy_kernels;
        accel_surf_vz_kernels = tensor_accel_surf_vz_kernels;
        accel_boundary_surf_vx_kernels = tensor_accel_boundary_surf_vx_kernels;
        accel_boundary_surf_vy_kernels = tensor_accel_boundary_surf_vy_kernels;
        accel_boundary_surf_vz_kernels = tensor_accel_boundary_surf_vz_kernels;
      } else {
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

  if (inp->model_id == GKYL_MODEL_TRIAD || inp->hamil_id == GKYL_HAMIL_PHASE) {
    vlasov->stream_surf_from_flux[0] =
      stream_surf_from_flux_x_kernels[kernel_index].kernels[poly_order];
    vlasov->stream_surf_from_flux[1] =
      stream_surf_from_flux_y_kernels[kernel_index].kernels[poly_order];
    vlasov->stream_surf_from_flux[2] =
      stream_surf_from_flux_z_kernels[kernel_index].kernels[poly_order];

    vlasov->stream_boundary_surf_from_flux[0] =
      stream_boundary_surf_from_flux_x_kernels[kernel_index].kernels[poly_order];
    vlasov->stream_boundary_surf_from_flux[1] =
      stream_boundary_surf_from_flux_y_kernels[kernel_index].kernels[poly_order];
    vlasov->stream_boundary_surf_from_flux[2] =
      stream_boundary_surf_from_flux_z_kernels[kernel_index].kernels[poly_order];
  } else {
    vlasov->stream_surf[0] = stream_surf_x_kernels[kernel_index].kernels[poly_order];
    vlasov->stream_surf[1] = stream_surf_y_kernels[kernel_index].kernels[poly_order];
    vlasov->stream_surf[2] = stream_surf_z_kernels[kernel_index].kernels[poly_order];

    vlasov->stream_boundary_surf[0] =
      stream_boundary_surf_x_kernels[kernel_index].kernels[poly_order];
    vlasov->stream_boundary_surf[1] =
      stream_boundary_surf_y_kernels[kernel_index].kernels[poly_order];
    vlasov->stream_boundary_surf[2] =
      stream_boundary_surf_z_kernels[kernel_index].kernels[poly_order];
  }

  vlasov->accel_surf[0] = accel_surf_vx_kernels[kernel_index].kernels[poly_order];
  vlasov->accel_surf[1] = accel_surf_vy_kernels[kernel_index].kernels[poly_order];
  vlasov->accel_surf[2] = accel_surf_vz_kernels[kernel_index].kernels[poly_order];

  vlasov->accel_boundary_surf[0] = accel_boundary_surf_vx_kernels[kernel_index].kernels[poly_order];
  vlasov->accel_boundary_surf[1] = accel_boundary_surf_vy_kernels[kernel_index].kernels[poly_order];
  vlasov->accel_boundary_surf[2] = accel_boundary_surf_vz_kernels[kernel_index].kernels[poly_order];

  // ensure non-NULL pointers
  assert(vlasov->hamil_vol);
  assert(vlasov->E_vol && vlasov->phi_vol && vlasov->rad_vol);
  assert(vlasov->Bx_vol && vlasov->By_vol && vlasov->Bz_vol);
  for (int i = 0; i < cdim; ++i) {
    if (inp->model_id == GKYL_MODEL_TRIAD || inp->hamil_id == GKYL_HAMIL_PHASE) {
      assert(vlasov->stream_surf_from_flux[i]);
      assert(vlasov->stream_boundary_surf_from_flux[i]);
    } else {
      assert(vlasov->stream_surf[i]);
      assert(vlasov->stream_boundary_surf[i]);
    }
  }
  for (int i = 0; i < vdim; ++i) {
    assert(vlasov->accel_surf[i]);
    assert(vlasov->accel_boundary_surf[i]);
  }

  vlasov->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(vlasov->eqn.flags);

  vlasov->eqn.ref_count = gkyl_ref_count_init(gkyl_vlasov_free);
  vlasov->eqn.on_dev = &vlasov->eqn; // CPU eqn obj points to itself

  return &vlasov->eqn;
}
