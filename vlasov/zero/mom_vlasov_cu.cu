/* -*- c++ -*- */

#include <assert.h>
#include <stdlib.h>
#include <string.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_mom_vlasov.h>
#include <gkyl_mom_vlasov_priv.h>
#include <gkyl_util.h>
}

static int
v_num_mom(int vdim, enum gkyl_distribution_moments mom_type)
{
  int m3ijk_count[] = { 1, 4, 10 };
  int num_mom = 0;
  
  switch (mom_type) {
    case GKYL_F_MOMENT_M0:
    case GKYL_F_MOMENT_ENERGY:
    case GKYL_F_MOMENT_M2:
    case GKYL_F_MOMENT_M0_UPPER:
    case GKYL_F_MOMENT_M0_LOWER:
      num_mom = 1;
      break;

    case GKYL_F_MOMENT_M1:
    case GKYL_F_MOMENT_M1_FROM_H:
    case GKYL_F_MOMENT_M3:
    case GKYL_F_MOMENT_ENERGY_FLUX:
      num_mom = vdim;
      break;

    case GKYL_F_MOMENT_M2IJ:
      num_mom = vdim*(vdim+1)/2;
      break;

    case GKYL_F_MOMENT_M3IJK:
      num_mom = m3ijk_count[vdim-1];
      break;

    case GKYL_F_MOMENT_M0M1M2:
      num_mom = vdim+2;
      break;      
      
    default: // Can't happen.
      fprintf(stderr,"Moment option %d not available.\n",mom_type);
      assert(false);
      break;
  }

  return num_mom;
}

__global__
static void
set_cu_ptrs(struct mom_type_vlasov* mom_vlasov, enum gkyl_distribution_moments mom_type,
  enum gkyl_basis_type b_type, int cdim, int vdim, int poly_order, 
  enum gkyl_model_id model_id, enum gkyl_hamil_id hamil_id, const struct gkyl_array *hamil)
{
  int m3ijk_count[] = { 1, 4, 10 };

  // Sparse (separable) vs. dense velocity-space Hamiltonian kernel selection.
  bool hamil_sparse = (hamil_id == GKYL_HAMIL_VEL_SPARSE);

  // choose kernel tables based on basis-function type
  const gkyl_vlasov_mom_kern_list *m0_kernels, *m1i_hamil_vel_kernels, *m1i_hamil_phase_kernels,
    *m2_hamil_vel_kernels, *m2_hamil_phase_kernels, *m3i_hamil_vel_kernels, 
    *m2ij_kernels, *m3ijk_kernels, *five_moments_hamil_vel_kernels, *five_moments_hamil_phase_kernels;

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      m0_kernels = ser_m0_kernels;
      m2ij_kernels = ser_m2ij_kernels;
      m3ijk_kernels = ser_m3ijk_kernels;
      m1i_hamil_vel_kernels = hamil_sparse ? ser_hamil_vel_sparse_m1i_kernels : ser_hamil_vel_dense_m1i_kernels;
      m2_hamil_vel_kernels = hamil_sparse ? ser_hamil_vel_sparse_m2_kernels : ser_hamil_vel_dense_m2_kernels;
      m3i_hamil_vel_kernels = hamil_sparse ? ser_hamil_vel_sparse_m3i_kernels : ser_hamil_vel_dense_m3i_kernels;
      five_moments_hamil_vel_kernels = hamil_sparse ? ser_hamil_vel_sparse_five_moments_kernels : ser_hamil_vel_dense_five_moments_kernels;
      m1i_hamil_phase_kernels = ser_hamil_phase_m1i_kernels;
      m2_hamil_phase_kernels = ser_hamil_phase_m2_kernels;
      five_moments_hamil_phase_kernels = ser_hamil_phase_five_moments_kernels;
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      m0_kernels = tensor_m0_kernels;
      m2ij_kernels = tensor_m2ij_kernels;
      m3ijk_kernels = tensor_m3ijk_kernels;
      m1i_hamil_vel_kernels = hamil_sparse ? tensor_hamil_vel_sparse_m1i_kernels : tensor_hamil_vel_dense_m1i_kernels;
      m2_hamil_vel_kernels = hamil_sparse ? tensor_hamil_vel_sparse_m2_kernels : tensor_hamil_vel_dense_m2_kernels;
      m3i_hamil_vel_kernels = hamil_sparse ? tensor_hamil_vel_sparse_m3i_kernels : tensor_hamil_vel_dense_m3i_kernels;
      five_moments_hamil_vel_kernels = hamil_sparse ? tensor_hamil_vel_sparse_five_moments_kernels : tensor_hamil_vel_dense_five_moments_kernels;
      // Phase-space Hamiltonian moments: only the p=1 tensor hybrid has a
      // phase-space Hamiltonian representation.
      m1i_hamil_phase_kernels = tensor_hamil_phase_m1i_kernels;
      m2_hamil_phase_kernels = tensor_hamil_phase_m2_kernels;
      five_moments_hamil_phase_kernels = tensor_hamil_phase_five_moments_kernels;
      break;

    default:
      assert(false);
      break;    
  }  

  switch (mom_type) {
    case GKYL_F_MOMENT_M0:
      mom_vlasov->momt.kernel = m0_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
      mom_vlasov->momt.num_mom = 1;
      break;

    case GKYL_F_MOMENT_M1:
    case GKYL_F_MOMENT_M1_FROM_H:
      if (hamil_id != GKYL_HAMIL_PHASE) {
        mom_vlasov->momt.kernel = m1i_hamil_vel_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
      }
      else {
        mom_vlasov->momt.kernel = m1i_hamil_phase_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
      }
      mom_vlasov->momt.num_mom = vdim;
      break;

    case GKYL_F_MOMENT_M2:
    case GKYL_F_MOMENT_ENERGY:
      if (hamil_id != GKYL_HAMIL_PHASE) {
        mom_vlasov->momt.kernel = m2_hamil_vel_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
      }
      else {
        mom_vlasov->momt.kernel = m2_hamil_phase_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
      }
      mom_vlasov->momt.num_mom = 1;
      break;

    case GKYL_F_MOMENT_M3:
    case GKYL_F_MOMENT_ENERGY_FLUX:
      if (hamil_id != GKYL_HAMIL_PHASE && model_id != GKYL_MODEL_TRIAD_GR) {
        mom_vlasov->momt.kernel = m3i_hamil_vel_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
      }
      else {
       assert(false); 
      }
      mom_vlasov->momt.num_mom = vdim;
      break;

    case GKYL_F_MOMENT_M2IJ:
      mom_vlasov->momt.kernel = m2ij_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
      mom_vlasov->momt.num_mom = vdim*(vdim+1)/2;
      break;

    case GKYL_F_MOMENT_M3IJK:
      mom_vlasov->momt.kernel = m3ijk_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
      mom_vlasov->momt.num_mom = m3ijk_count[vdim-1];
      break;

    case GKYL_F_MOMENT_M0M1M2:
      if (hamil_id != GKYL_HAMIL_PHASE) {
        mom_vlasov->momt.kernel = five_moments_hamil_vel_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
      }
      else {
        mom_vlasov->momt.kernel = five_moments_hamil_phase_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
      }
      mom_vlasov->momt.num_mom = vdim+2;
      break;

    case GKYL_F_MOMENT_M0_UPPER:
      mom_vlasov->momt.kernel = (b_type == GKYL_BASIS_MODAL_TENSOR) ?
        tensor_m0_upper_kernels[cdim-1].kernels[poly_order] :
        ser_m0_upper_kernels[cdim-1].kernels[poly_order];
      mom_vlasov->momt.num_mom = 1;
      break;

    case GKYL_F_MOMENT_M0_LOWER:
      mom_vlasov->momt.kernel = (b_type == GKYL_BASIS_MODAL_TENSOR) ?
        tensor_m0_lower_kernels[cdim-1].kernels[poly_order] :
        ser_m0_lower_kernels[cdim-1].kernels[poly_order];
      mom_vlasov->momt.num_mom = 1;
      break;

    default: // can't happen
      break;
  }
}

struct gkyl_mom_type*
gkyl_mom_vlasov_cu_dev_inew(const struct gkyl_mom_vlasov_inp *inp)
{
  assert(inp->conf_basis->poly_order == inp->phase_basis->poly_order);

  struct mom_type_vlasov *mom_vlasov = (struct mom_type_vlasov*) gkyl_malloc(sizeof(*mom_vlasov));
  int cdim = inp->conf_basis->ndim, pdim = inp->phase_basis->ndim, vdim = pdim-cdim;
  int poly_order = inp->conf_basis->poly_order;

  mom_vlasov->momt.cdim = cdim;
  mom_vlasov->momt.pdim = pdim;
  mom_vlasov->momt.poly_order = poly_order;
  mom_vlasov->momt.num_config = inp->conf_basis->num_basis;
  mom_vlasov->momt.num_phase = inp->phase_basis->num_basis;

  // Determine Hamiltonian dimensionality and index offset for indexing Hamiltonian
  // from an input phase space index. 
  if (inp->hamil_id == GKYL_HAMIL_PHASE) {
    mom_vlasov->hamil_dim = pdim;
    mom_vlasov->hamil_offset = 0;
  }
  else {
    mom_vlasov->hamil_dim = vdim;
    mom_vlasov->hamil_offset = cdim;
  }
  mom_vlasov->hamil_range = *inp->hamil_range;
  struct gkyl_array *hamil_ho = gkyl_array_acquire(inp->hamil); 
  mom_vlasov->hamil = hamil_ho->on_dev; // store pointer to on_dev for copying over to device. 

  mom_vlasov->vel_range = *inp->vel_range;
  mom_vlasov->vel_map = 0;
  mom_vlasov->vmap = 0;
  mom_vlasov->jacob_vel = 0;
  mom_vlasov->use_vmap = false;
  if (inp->vel_map) {
    mom_vlasov->use_vmap = inp->vel_map->is_mapped;
    if (mom_vlasov->use_vmap) {
      // Unpack raw device pointers for use inside kernels; lifetime is
      // guaranteed by acquiring the vel_map object on the host side below.
      mom_vlasov->vmap = inp->vel_map->vmap->on_dev;
      mom_vlasov->jacob_vel = inp->vel_map->jacob_vel->on_dev;
    }
  }

  // Threshold velocity for integration of moments over a subset of the domain. 
  // Also set threshold for whether we accumulate moment over subset of the domain. 
  mom_vlasov->v_thresh = inp->v_thresh > 0.0 ? inp->v_thresh : 0.0; 
  mom_vlasov->f_thresh = inp->f_thresh > 0.0 ? inp->f_thresh : 0.0; 

  mom_vlasov->momt.num_mom = v_num_mom(vdim, inp->mom_type); // Number of moments.

  mom_vlasov->momt.flags = 0;
  GKYL_SET_CU_ALLOC(mom_vlasov->momt.flags);
  mom_vlasov->momt.ref_count = gkyl_ref_count_init(gkyl_mom_vlasov_free);
  
  // copy struct to device
  struct mom_type_vlasov *mom_vlasov_cu = (struct mom_type_vlasov*) gkyl_cu_malloc(sizeof(*mom_vlasov_cu));
  gkyl_cu_memcpy(mom_vlasov_cu, mom_vlasov, sizeof(struct mom_type_vlasov), GKYL_CU_MEMCPY_H2D);

  set_cu_ptrs<<<1,1>>>(mom_vlasov_cu, inp->mom_type, inp->conf_basis->b_type, 
    cdim, vdim, poly_order, inp->model_id, inp->hamil_id, inp->hamil->on_dev);

  mom_vlasov->momt.on_dev = &mom_vlasov_cu->momt;

  // Host-side moment type object should store host pointers.
  mom_vlasov->vel_map = inp->vel_map ? gkyl_vlasov_velocity_map_acquire(inp->vel_map) : 0;
  mom_vlasov->vmap = (inp->vel_map && mom_vlasov->use_vmap) ? inp->vel_map->vmap : 0;
  mom_vlasov->jacob_vel = (inp->vel_map && mom_vlasov->use_vmap) ? inp->vel_map->jacob_vel : 0;
  mom_vlasov->hamil = hamil_ho; 
  
  return &mom_vlasov->momt;
}

__global__
static void
set_int_cu_ptrs(struct mom_type_vlasov* mom_vlasov, enum gkyl_distribution_moments mom_type,
  enum gkyl_basis_type b_type, int cdim, int vdim, int poly_order, 
  enum gkyl_model_id model_id, enum gkyl_hamil_id hamil_id, const struct gkyl_array *hamil)
{
  // Sparse (separable) vs. dense velocity-space Hamiltonian kernel selection.
  bool hamil_sparse = (hamil_id == GKYL_HAMIL_VEL_SPARSE);

  // Choose kernel tables based on basis-function type.
  const gkyl_vlasov_mom_kern_list *int_five_moments_hamil_vel_kernels, *int_five_moments_hamil_phase_kernels;

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      int_five_moments_hamil_vel_kernels = hamil_sparse ? ser_hamil_vel_sparse_int_five_moments_kernels : ser_hamil_vel_dense_int_five_moments_kernels;
      int_five_moments_hamil_phase_kernels = ser_hamil_phase_int_five_moments_kernels;
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      int_five_moments_hamil_vel_kernels = hamil_sparse ? tensor_hamil_vel_sparse_int_five_moments_kernels : tensor_hamil_vel_dense_int_five_moments_kernels;
      // Phase-space Hamiltonian integrated moments: p=1 tensor hybrid only.
      int_five_moments_hamil_phase_kernels = tensor_hamil_phase_int_five_moments_kernels;
      break;

    default:
      assert(false);
      break;    
  }   

  switch (mom_type) {
    case GKYL_F_MOMENT_M0M1M2:
      if (hamil_id != GKYL_HAMIL_PHASE) {
        mom_vlasov->momt.kernel = int_five_moments_hamil_vel_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
      }
      else {
        mom_vlasov->momt.kernel = int_five_moments_hamil_phase_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
      } 
      mom_vlasov->momt.num_mom = 2+vdim;
      break;

    default:
      assert(false);
      break;
  }
}

struct gkyl_mom_type*
gkyl_int_mom_vlasov_cu_dev_inew(const struct gkyl_mom_vlasov_inp *inp)
{
  assert(inp->conf_basis->poly_order == inp->phase_basis->poly_order);

  struct mom_type_vlasov *mom_vlasov = (struct mom_type_vlasov*) gkyl_malloc(sizeof(*mom_vlasov));
  int cdim = inp->conf_basis->ndim, pdim = inp->phase_basis->ndim, vdim = pdim-cdim;
  int poly_order = inp->conf_basis->poly_order;

  mom_vlasov->momt.cdim = cdim;
  mom_vlasov->momt.pdim = pdim;
  mom_vlasov->momt.poly_order = poly_order;
  mom_vlasov->momt.num_config = inp->conf_basis->num_basis;
  mom_vlasov->momt.num_phase = inp->phase_basis->num_basis;
  // Determine Hamiltonian dimensionality and index offset for indexing Hamiltonian
  // from an input phase space index. 
  if (inp->hamil_id == GKYL_HAMIL_PHASE) {
    mom_vlasov->hamil_dim = pdim;
    mom_vlasov->hamil_offset = 0;
  }
  else {
    mom_vlasov->hamil_dim = vdim;
    mom_vlasov->hamil_offset = cdim;
  }
  mom_vlasov->hamil_range = *inp->hamil_range;
  struct gkyl_array *hamil_ho = gkyl_array_acquire(inp->hamil); 
  mom_vlasov->hamil = hamil_ho->on_dev; // store pointer to on_dev for copying over to device. 

  mom_vlasov->vel_range = *inp->vel_range;
  mom_vlasov->vel_map = 0;
  mom_vlasov->vmap = 0;
  mom_vlasov->jacob_vel = 0;
  mom_vlasov->use_vmap = false;
  if (inp->vel_map) {
    mom_vlasov->use_vmap = inp->vel_map->is_mapped;
    if (mom_vlasov->use_vmap) {
      // Unpack raw device pointers for use inside kernels; lifetime is
      // guaranteed by acquiring the vel_map object on the host side below.
      mom_vlasov->vmap = inp->vel_map->vmap->on_dev;
      mom_vlasov->jacob_vel = inp->vel_map->jacob_vel->on_dev;
    }
  }

  mom_vlasov->momt.num_mom = v_num_mom(vdim, inp->mom_type); // Number of moments.

  mom_vlasov->momt.flags = 0;
  GKYL_SET_CU_ALLOC(mom_vlasov->momt.flags);
  mom_vlasov->momt.ref_count = gkyl_ref_count_init(gkyl_mom_vlasov_free);
  
  // copy struct to device
  struct mom_type_vlasov *mom_vlasov_cu = (struct mom_type_vlasov*) gkyl_cu_malloc(sizeof(*mom_vlasov_cu));
  gkyl_cu_memcpy(mom_vlasov_cu, mom_vlasov, sizeof(struct mom_type_vlasov), GKYL_CU_MEMCPY_H2D);

  set_int_cu_ptrs<<<1,1>>>(mom_vlasov_cu, inp->mom_type, inp->conf_basis->b_type, 
    cdim, vdim, poly_order, inp->model_id, inp->hamil_id, inp->hamil->on_dev);

  mom_vlasov->momt.on_dev = &mom_vlasov_cu->momt;

  // Host-side moment type object should store host pointers.
  mom_vlasov->vel_map = inp->vel_map ? gkyl_vlasov_velocity_map_acquire(inp->vel_map) : 0;
  mom_vlasov->vmap = (inp->vel_map && mom_vlasov->use_vmap) ? inp->vel_map->vmap : 0;
  mom_vlasov->jacob_vel = (inp->vel_map && mom_vlasov->use_vmap) ? inp->vel_map->jacob_vel : 0;
  mom_vlasov->hamil = hamil_ho;  

  return &mom_vlasov->momt;
}
