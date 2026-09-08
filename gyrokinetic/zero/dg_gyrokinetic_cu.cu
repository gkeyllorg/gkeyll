/* -*- c++ -*- */
#include <float.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_gyrokinetic.h>    
#include <gkyl_dg_gyrokinetic_priv.h>
}

#include <cassert>

// CUDA kernel to set pointer to auxiliary fields.
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_gyrokinetic_set_auxfields_cu_kernel(const struct gkyl_dg_eqn *eqn, 
  const struct gkyl_array *flux_surf, 
  const struct gkyl_array *phi, const struct gkyl_array *apar, const struct gkyl_array *apardot)
{
  struct dg_gyrokinetic *gyrokinetic = container_of(eqn, struct dg_gyrokinetic, eqn);
  gyrokinetic->auxfields.flux_surf = flux_surf;
  gyrokinetic->auxfields.phi = phi;
  gyrokinetic->auxfields.apar = apar;
  gyrokinetic->auxfields.apardot = apardot;
}

// Host-side wrapper for set_auxfields_cu_kernel
void
gkyl_gyrokinetic_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_gyrokinetic_auxfields auxin)
{
  gkyl_gyrokinetic_set_auxfields_cu_kernel<<<1,1>>>(eqn, 
    auxin.flux_surf->on_dev, 
    auxin.phi->on_dev, auxin.apar->on_dev, auxin.apardot->on_dev);
}

// CUDA kernel to set device pointers to range object and gyrokinetic kernel function
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
dg_gyrokinetic_set_cu_dev_ptrs(struct dg_gyrokinetic *gyrokinetic, enum gkyl_basis_type b_type,
  int cv_index, int cdim, int vdim, int poly_order, enum gkyl_gk_collisionless_type collless_type, 
  bool no_by, bool complete_em)
{
  gyrokinetic->auxfields.flux_surf = 0; 
  gyrokinetic->auxfields.phi = 0; 
  gyrokinetic->auxfields.apar = 0; 
  gyrokinetic->auxfields.apardot= 0; 

  gyrokinetic->eqn.surf_term = surf;
  gyrokinetic->eqn.boundary_surf_term = boundary_surf;
  gyrokinetic->eqn.boundary_diag_term = boundary_diag;

  const gkyl_dg_gyrokinetic_vol_es_kern_list *vol_kernels, *vol_no_by_kernels;
  const gkyl_dg_gyrokinetic_vol_add_apar_kern_list *vol_add_apar_kernels;
  const gkyl_dg_gyrokinetic_vol_add_apardot_kern_list *vol_add_apardot_kernels;
  const gkyl_dg_gyrokinetic_surf_kern_list *surf_x_kernels; 
  const gkyl_dg_gyrokinetic_surf_kern_list *surf_y_kernels; 
  const gkyl_dg_gyrokinetic_surf_kern_list *surf_z_kernels; 
  const gkyl_dg_gyrokinetic_surf_kern_list *surf_vpar_kernels; 
  const gkyl_dg_gyrokinetic_boundary_surf_kern_list *boundary_surf_x_kernels; 
  const gkyl_dg_gyrokinetic_boundary_surf_kern_list *boundary_surf_y_kernels; 
  const gkyl_dg_gyrokinetic_boundary_surf_kern_list *boundary_surf_z_kernels; 
  const gkyl_dg_gyrokinetic_boundary_surf_kern_list *boundary_surf_vpar_kernels; 
  
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      vol_kernels = ser_vol_es_kernels;
      vol_add_apar_kernels = ser_add_apar_vol_kernels;
      vol_add_apardot_kernels = ser_add_apardot_vol_kernels;
      surf_x_kernels = ser_surf_x_kernels;
      surf_y_kernels = ser_surf_y_kernels;
      surf_z_kernels = ser_surf_z_kernels;
      surf_vpar_kernels = ser_surf_vpar_kernels;
      boundary_surf_x_kernels = ser_boundary_surf_x_kernels;
      boundary_surf_y_kernels = ser_boundary_surf_y_kernels;
      boundary_surf_z_kernels = ser_boundary_surf_z_kernels;
      boundary_surf_vpar_kernels = ser_boundary_surf_vpar_kernels;

      vol_no_by_kernels = ser_no_by_vol_es_kernels;
      
      break;

    default:
      assert(false);
      break;    
  }  
  
  gyrokinetic->eqn.vol_term = kernel_dg_gyrokinetic_vol;

  if (no_by) {
    gyrokinetic->vol_es_kernel = vol_no_by_kernels[cv_index].kernels[poly_order];
  }
  else {
    gyrokinetic->vol_es_kernel = vol_kernels[cv_index].kernels[poly_order];
  }

  // Setup electromagnetic terms if needed.
  bool is_em = (collless_type == GKYL_GK_COLLISIONLESS_EM);
  gyrokinetic->vol_add_apar_kernel = dg_gyrokinetic_add_apar_vol_none;
  gyrokinetic->vol_add_apardot_kernel = dg_gyrokinetic_add_apardot_vol_none;
  if (is_em) {
    if (complete_em) {
      // We complete the Ohm's law RHS with Apardot contribution to get the full GK RHS.
      gyrokinetic->vol_es_kernel = dg_gyrokinetic_vol_none;
      gyrokinetic->vol_add_apar_kernel = dg_gyrokinetic_add_apar_vol_none;
      gyrokinetic->vol_add_apardot_kernel = vol_add_apardot_kernels[cv_index].kernels[poly_order];
    } else {
      // We build Ohm's law RHS (no Apardot contribution).
      gyrokinetic->vol_es_kernel = gyrokinetic->vol_es_kernel; // no change to ES kernel
      gyrokinetic->vol_add_apar_kernel = vol_add_apar_kernels[cv_index].kernels[poly_order];
      gyrokinetic->vol_add_apardot_kernel = dg_gyrokinetic_add_apardot_vol_none;
    }
  }

  gyrokinetic->surf[0] = surf_x_kernels[cv_index].kernels[poly_order];
  if (cdim>1)
    gyrokinetic->surf[1] = surf_y_kernels[cv_index].kernels[poly_order];
  if (cdim>2)
    gyrokinetic->surf[2] = surf_z_kernels[cv_index].kernels[poly_order];
  gyrokinetic->surf[cdim] = surf_vpar_kernels[cv_index].kernels[poly_order];

  gyrokinetic->boundary_surf[0] = boundary_surf_x_kernels[cv_index].kernels[poly_order];
  if (cdim>1)
    gyrokinetic->boundary_surf[1] = boundary_surf_y_kernels[cv_index].kernels[poly_order];
  if (cdim>2)
    gyrokinetic->boundary_surf[2] = boundary_surf_z_kernels[cv_index].kernels[poly_order];
  gyrokinetic->boundary_surf[cdim] = boundary_surf_vpar_kernels[cv_index].kernels[poly_order];
  
}

struct gkyl_dg_eqn*
gkyl_dg_gyrokinetic_cu_dev_new(const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis,
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, 
  const double charge, const double mass, enum gkyl_gk_collisionless_type collless_type,
  const bool no_by, const bool complete_em, const struct gk_geometry *gk_geom, const struct gkyl_velocity_map *vel_map)
{
  struct dg_gyrokinetic *gyrokinetic = (struct dg_gyrokinetic*) gkyl_malloc(sizeof(*gyrokinetic));

  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  gyrokinetic->cdim = cdim;
  gyrokinetic->pdim = pdim;

  gyrokinetic->charge = charge;
  gyrokinetic->mass = mass;

  gyrokinetic->eqn.num_equations = 1;

  // Acquire pointers to on_dev objects so memcpy below copies those too.
  struct gk_geometry *geom_ho = gkyl_gk_geometry_acquire(gk_geom);
  struct gkyl_velocity_map *vel_map_ho = gkyl_velocity_map_acquire(vel_map);
  
  gyrokinetic->gk_geom = geom_ho->on_dev;
  gyrokinetic->vel_map = vel_map_ho->on_dev;

  gyrokinetic->conf_range = *conf_range;
  gyrokinetic->phase_range = *phase_range;

  gyrokinetic->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(gyrokinetic->eqn.flags);
  gyrokinetic->eqn.ref_count = gkyl_ref_count_init(gkyl_gyrokinetic_free);

  // copy the host struct to device struct
  struct dg_gyrokinetic *gyrokinetic_cu = (struct dg_gyrokinetic*) gkyl_cu_malloc(sizeof(struct dg_gyrokinetic));
  gkyl_cu_memcpy(gyrokinetic_cu, gyrokinetic, sizeof(struct dg_gyrokinetic), GKYL_CU_MEMCPY_H2D);

  dg_gyrokinetic_set_cu_dev_ptrs<<<1,1>>>(gyrokinetic_cu, cbasis->b_type, cv_index[cdim].vdim[vdim],
    cdim, vdim, poly_order, collless_type, no_by, complete_em);

  // set parent on_dev pointer
  gyrokinetic->eqn.on_dev = &gyrokinetic_cu->eqn;
  
  // Updater should store host pointers.
  gyrokinetic->gk_geom = geom_ho; 
  gyrokinetic->vel_map = vel_map_ho;

  return &gyrokinetic->eqn;
}