/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_gk_anomalous_diffusion.h>    
#include <gkyl_gk_anomalous_diffusion_priv.h>
}

#include <cassert>

// CUDA kernel to set pointer to auxiliary fields.
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_gk_anomalous_diffusion_set_auxfields_cu_kernel(const struct gkyl_dg_eqn* eqn,
  const struct gkyl_array* nu, const struct gkyl_array* jacobgeo_inv)
{
  struct gk_anomalous_diffusion* diffusion = container_of(eqn, struct gk_anomalous_diffusion, eqn);
  diffusion->auxfields.nu = nu;
  diffusion->auxfields.jacobgeo_inv = jacobgeo_inv;
}

// Host-side wrapper for set_auxfields_cu_kernel
void
gkyl_gk_anomalous_diffusion_set_auxfields_cu(const struct gkyl_dg_eqn* eqn,
  struct gkyl_gk_anomalous_diffusion_auxfields auxin)
{
  gkyl_gk_anomalous_diffusion_set_auxfields_cu_kernel<<<1,1>>>(eqn, auxin.nu->on_dev, auxin.jacobgeo_inv->on_dev);
}

__global__ void static
gk_anomalous_diffusion_set_cu_dev_ptrs(struct gk_anomalous_diffusion *diffusion, enum gkyl_basis_type b_type,
  int cdim, int vdim, int poly_order, enum gkyl_gyrokinetic_bc_type bc_x_lower, enum gkyl_gyrokinetic_bc_type bc_x_upper)
{
  int pdim = cdim + vdim;

  diffusion->auxfields.nu = 0; 
  diffusion->auxfields.jacobgeo_inv = 0; 

  const gkyl_gk_anomalous_diffusion_vol_kern_list *vol_kernels;
  const gkyl_gk_anomalous_diffusion_surf_kern_list *surfx_kernels;
  const gkyl_gk_anomalous_diffusion_boundary_surf_kern_list *boundary_surfx_lower_kernels, *boundary_surfx_upper_kernels;
  const gkyl_gk_anomalous_diffusion_boundary_surf_kern_list *boundary_diagx_lower_kernels, *boundary_diagx_upper_kernels;

  // Choice of boundary_surf and boundary_diag kernels:
  //   boundary_surf: zero_flux or local
  //   boundary_diag: local or recovery
  // MF 2025/09/10: as of now these options are meant for (here
  // N/A means not applicable):
  //            bound_surf  bound_diag  hyper_dg-zero_flux
  // SKIP:      N/A         recovery    no
  // PERIODIC:  N/A         N/A         no
  // ZERO_FLUX: zero_flux   N/A         yes
  // ELSE:      local       local       yes

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      vol_kernels   = ser_vol_kernels;
      surfx_kernels = ser_gyrokinetic_surfx_kernels;
      if ((bc_x_lower == GKYL_BC_GK_SKIP) ||
          (bc_x_lower == GKYL_BC_GK_SPECIES_PERIODIC)) {
        // Boundary surf kernel not used.
        boundary_surfx_lower_kernels = ser_gyrokinetic_boundary_surfx_lower_zeroflux_kernels;
        boundary_diagx_lower_kernels = ser_gyrokinetic_boundary_diagx_lower_boundrecovery_kernels;
      }
      else if (bc_x_lower == GKYL_BC_GK_SPECIES_ZERO_FLUX) {
        boundary_surfx_lower_kernels = ser_gyrokinetic_boundary_surfx_lower_zeroflux_kernels;
        // Boundary diag kernel not used.
        boundary_diagx_lower_kernels = ser_gyrokinetic_boundary_diagx_lower_boundrecovery_kernels;
      }
      else if ((bc_x_lower == GKYL_BC_GK_SPECIES_ABSORB) ||
               (bc_x_lower == GKYL_BC_GK_SPECIES_FIXED_FUNC)) {
        // Boundary surf kernel not used.
        boundary_surfx_lower_kernels = ser_gyrokinetic_boundary_surfx_lower_zeroflux_kernels;
        boundary_diagx_lower_kernels = ser_gyrokinetic_boundary_diagx_lower_boundrecovery_kernels;
      }
      else {
        boundary_surfx_lower_kernels = ser_gyrokinetic_boundary_surfx_lower_boundlocal_kernels;
        boundary_diagx_lower_kernels = ser_gyrokinetic_boundary_diagx_lower_boundlocal_kernels;
      }

      if ((bc_x_upper == GKYL_BC_GK_SKIP) ||
          (bc_x_upper == GKYL_BC_GK_SPECIES_PERIODIC)) {
        // Boundary surf and diag kernels not used.
        boundary_surfx_upper_kernels = ser_gyrokinetic_boundary_surfx_upper_zeroflux_kernels;
        boundary_diagx_upper_kernels = ser_gyrokinetic_boundary_diagx_upper_boundrecovery_kernels;
      }
      else if (bc_x_upper == GKYL_BC_GK_SPECIES_ZERO_FLUX) {
        boundary_surfx_upper_kernels = ser_gyrokinetic_boundary_surfx_upper_zeroflux_kernels;
        // Boundary diag kernel not used.
        boundary_diagx_upper_kernels = ser_gyrokinetic_boundary_diagx_upper_boundrecovery_kernels;
      }
      else if ((bc_x_upper == GKYL_BC_GK_SPECIES_ABSORB) ||
               (bc_x_upper == GKYL_BC_GK_SPECIES_FIXED_FUNC)) {
        // Boundary surf kernel not used.
        boundary_surfx_upper_kernels = ser_gyrokinetic_boundary_surfx_upper_zeroflux_kernels;
        boundary_diagx_upper_kernels = ser_gyrokinetic_boundary_diagx_upper_boundrecovery_kernels;
      }
      else {
        boundary_surfx_upper_kernels = ser_gyrokinetic_boundary_surfx_upper_boundlocal_kernels;
        boundary_diagx_upper_kernels = ser_gyrokinetic_boundary_diagx_upper_boundlocal_kernels;
      }
      break;

    default:
      assert(false);
      break;
  }

  diffusion->eqn.num_equations = 1;
  diffusion->eqn.surf_term = surf;
  diffusion->eqn.boundary_surf_term = boundary_surf;
  diffusion->eqn.boundary_diag_term = boundary_diag;

  diffusion->eqn.vol_term = CKVOL(vol_kernels, pdim, poly_order);
  diffusion->surf = CKSURF(surfx_kernels, pdim, poly_order);
  diffusion->boundary_surf[0] = CKSURF(boundary_surfx_lower_kernels, pdim,  poly_order);
  diffusion->boundary_surf[1] = CKSURF(boundary_surfx_upper_kernels, pdim,  poly_order);
  diffusion->boundary_diag[0] = CKSURF(boundary_diagx_lower_kernels, pdim,  poly_order);
  diffusion->boundary_diag[1] = CKSURF(boundary_diagx_upper_kernels, pdim,  poly_order);

  // Ensure non-NULL pointers.
  assert(diffusion->eqn.vol_term);
  assert(diffusion->surf);
  for (int i=0; i<2; i++) {
    assert(diffusion->boundary_surf[i]);
    assert(diffusion->boundary_diag[i]);
  }
}

struct gkyl_dg_eqn*
gkyl_gk_anomalous_diffusion_cu_dev_new(const struct gkyl_basis *basis, const struct gkyl_basis *cbasis,
  const struct gkyl_range *conf_range, enum gkyl_gyrokinetic_bc_type bc_x_lower, enum gkyl_gyrokinetic_bc_type bc_x_upper)
{
  struct gk_anomalous_diffusion* diffusion = (struct gk_anomalous_diffusion*) gkyl_malloc(sizeof(struct gk_anomalous_diffusion));

  int cdim = cbasis->ndim;
  int vdim = basis->ndim - cdim;
  int poly_order = cbasis->poly_order;

  diffusion->conf_range = *conf_range;

  diffusion->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(diffusion->eqn.flags);
  diffusion->eqn.ref_count = gkyl_ref_count_init(gkyl_gk_anomalous_diffusion_free);

  // Copy the host struct to device struct.
  struct gk_anomalous_diffusion* diffusion_cu = (struct gk_anomalous_diffusion*) gkyl_cu_malloc(sizeof(struct gk_anomalous_diffusion));
  gkyl_cu_memcpy(diffusion_cu, diffusion, sizeof(struct gk_anomalous_diffusion), GKYL_CU_MEMCPY_H2D);

  gk_anomalous_diffusion_set_cu_dev_ptrs<<<1,1>>>(diffusion_cu, cbasis->b_type, cdim, vdim, poly_order, bc_x_lower, bc_x_upper);

  // Set parent on_dev pointer.
  diffusion->eqn.on_dev = &diffusion_cu->eqn;

  return &diffusion->eqn;
}
