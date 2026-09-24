#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_gk_anomalous_diffusion.h>
#include <gkyl_gk_anomalous_diffusion_priv.h>
#include <gkyl_util.h>

void
gkyl_gk_anomalous_diffusion_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_eqn *base = container_of(ref, struct gkyl_dg_eqn, ref_count);

  if (gkyl_dg_eqn_is_cu_dev(base)) {
    // free inner on_dev object
    struct gk_anomalous_diffusion *diffusion = container_of(base->on_dev, struct gk_anomalous_diffusion, eqn);
    gkyl_cu_free(diffusion);
  }
  struct gk_anomalous_diffusion *diffusion = container_of(base, struct gk_anomalous_diffusion, eqn);
  gkyl_free(diffusion);
}

void
gkyl_gk_anomalous_diffusion_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_gk_anomalous_diffusion_auxfields auxin)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(auxin.nu) && gkyl_array_is_cu_dev(auxin.jacobgeo_inv)) {
    gkyl_gk_anomalous_diffusion_set_auxfields_cu(eqn->on_dev, auxin);
    return;
  }
#endif
  
  struct gk_anomalous_diffusion *diffusion = container_of(eqn, struct gk_anomalous_diffusion, eqn);
  diffusion->auxfields.nu = auxin.nu;
  diffusion->auxfields.jacobgeo_inv = auxin.jacobgeo_inv;
}

struct gkyl_dg_eqn*
gkyl_gk_anomalous_diffusion_new(const struct gkyl_basis *basis, const struct gkyl_basis *cbasis,
  const struct gkyl_range *conf_range, enum gkyl_gyrokinetic_bc_type bc_x_lower, enum gkyl_gyrokinetic_bc_type bc_x_upper, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    return gkyl_gk_anomalous_diffusion_cu_dev_new(basis, cbasis, conf_range, bc_x_lower, bc_x_upper);
#endif
  
  struct gk_anomalous_diffusion *diffusion = gkyl_malloc(sizeof(struct gk_anomalous_diffusion));

  int cdim = cbasis->ndim;
  int vdim = basis->ndim - cdim;
  int pdim = cdim + vdim;
  int poly_order = cbasis->poly_order;

  const gkyl_gk_anomalous_diffusion_vol_kern_list *vol_kernels;
  const gkyl_gk_anomalous_diffusion_surf_kern_list *surfx_kernels;
  const gkyl_gk_anomalous_diffusion_boundary_surf_kern_list *boundary_surfx_lower_kernels, *boundary_surfx_upper_kernels;
  const gkyl_gk_anomalous_diffusion_boundary_surf_kern_list *boundary_diagx_lower_kernels, *boundary_diagx_upper_kernels;

  // Choice of boundary_surf and boundary_diag kernels:
  //   boundary_surf: zero_flux or local
  //   boundary_diag: local or recovery
  // MF 2025/09/10: as of now these options are meant for (here
  // N/A means not applicable):
  //             bound_surf  bound_diag  hyper_dg-zero_flux
  // SKIP:       N/A         recovery    no
  // ABSORB:     N/A         N/A         no
  // PERIODIC:   N/A         N/A         no
  // FIXED_FUNC: N/A         N/A         no
  // ZERO_FLUX:  zero_flux   N/A         yes
  // ELSE:       local       local       yes

  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      vol_kernels            = ser_vol_kernels;
      surfx_kernels          = ser_gyrokinetic_surfx_kernels;
      if ((bc_x_lower == GKYL_BC_GK_SKIP) ||
          (bc_x_lower == GKYL_BC_GK_SPECIES_PERIODIC)) {
        // Boundary surf and diag kernels not used.
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

  diffusion->eqn.vol_term = CKVOL(vol_kernels, cdim+vdim,  poly_order);
  diffusion->surf = CKSURF(surfx_kernels, cdim+vdim,  poly_order);
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

  diffusion->auxfields.nu = 0;
  diffusion->auxfields.jacobgeo_inv = 0;
  diffusion->conf_range = *conf_range;

  diffusion->eqn.flags = 0;
  diffusion->eqn.ref_count = gkyl_ref_count_init(gkyl_gk_anomalous_diffusion_free);
  diffusion->eqn.on_dev = &diffusion->eqn;
  
  return &diffusion->eqn;
}
