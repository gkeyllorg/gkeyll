/* -*- c++ -*- */
#include <float.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_lbo_gyrokinetic_drag.h>    
#include <gkyl_dg_lbo_gyrokinetic_drag_priv.h>
}

#include <cassert>

// CUDA kernel to set pointer to nuSum, sum of collisionalities
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_lbo_gyrokinetic_drag_set_auxfields_cu_kernel(const struct gkyl_dg_eqn *eqn, 
  const struct gkyl_array *nuSum, const struct gkyl_array *nuPrimMomsSum, const struct gkyl_array *m2self)
{
  struct dg_lbo_gyrokinetic_drag *lbo = container_of(eqn, struct dg_lbo_gyrokinetic_drag, eqn);
  lbo->auxfields.nuSum = nuSum;
  lbo->auxfields.nuPrimMomsSum = nuPrimMomsSum;
  lbo->auxfields.m2self = m2self;
}

//// Host-side wrapper for device kernels setting nuSum, nuUSum and nuVtSqSum.
void
gkyl_lbo_gyrokinetic_drag_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_lbo_gyrokinetic_drag_auxfields auxin)
{
  gkyl_lbo_gyrokinetic_drag_set_auxfields_cu_kernel<<<1,1>>>(eqn, 
    auxin.nuSum->on_dev, auxin.nuPrimMomsSum->on_dev, auxin.m2self->on_dev);
}

// CUDA kernel to set device pointers to range object and gyrokinetic LBO kernel function
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void
dg_lbo_gyrokinetic_drag_set_cu_dev_ptrs(struct dg_lbo_gyrokinetic_drag *lbo, enum gkyl_basis_type b_type,
  int cv_index, int cdim, int vdim, int poly_order)
{
  lbo->auxfields.nuSum = 0; 
  lbo->auxfields.nuPrimMomsSum = 0; 
  lbo->auxfields.m2self = 0; 

  lbo->eqn.surf_term = surf;
  lbo->eqn.boundary_surf_term = boundary_surf;

  const gkyl_dg_lbo_gyrokinetic_drag_vol_kern_list *vol_kernels;
  const gkyl_dg_lbo_gyrokinetic_drag_surf_kern_list *surf_vpar_kernels, *surf_mu_kernels;
  const gkyl_dg_lbo_gyrokinetic_drag_boundary_surf_kern_list *boundary_surf_vpar_kernels, *boundary_surf_mu_kernels;
  
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      vol_kernels = ser_vol_kernels;
      surf_vpar_kernels = ser_surf_vpar_kernels;
      surf_mu_kernels = ser_surf_mu_kernels;
      boundary_surf_vpar_kernels = ser_boundary_surf_vpar_kernels;
      boundary_surf_mu_kernels = ser_boundary_surf_mu_kernels;      
      break;

    default:
      assert(false);
      break;    
  }  

  lbo->eqn.vol_term = vol_kernels[cv_index].kernels[poly_order];

  lbo->surf[0] = surf_vpar_kernels[cv_index].kernels[poly_order];
  if (vdim>1)
    lbo->surf[1] = surf_mu_kernels[cv_index].kernels[poly_order];

  lbo->boundary_surf[0] = boundary_surf_vpar_kernels[cv_index].kernels[poly_order];
  if (vdim>1)
    lbo->boundary_surf[1] = boundary_surf_mu_kernels[cv_index].kernels[poly_order];

}

struct gkyl_dg_eqn*
gkyl_dg_lbo_gyrokinetic_drag_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const struct gkyl_range* conf_range, const struct gkyl_rect_grid *pgrid,
  double mass, const struct gk_geometry *gk_geom, const struct gkyl_velocity_map *vel_map)
{
  struct dg_lbo_gyrokinetic_drag *lbo =
    (struct dg_lbo_gyrokinetic_drag*) gkyl_malloc(sizeof(*lbo));

  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  lbo->cdim = cdim;
  lbo->pdim = pdim;

  lbo->eqn.num_equations = 1;

  lbo->mass = mass;
  lbo->conf_range = *conf_range;

  // Acquire pointers to on_dev objects so memcpy below copies those too.
  struct gk_geometry *geom_ho = gkyl_gk_geometry_acquire(gk_geom);
  struct gkyl_velocity_map *vel_map_ho = gkyl_velocity_map_acquire(vel_map);
  
  lbo->gk_geom = geom_ho->on_dev;
  lbo->vel_map = vel_map_ho->on_dev;

  lbo->vparMax = GKYL_MAX2(fabs(vel_map->vbounds[0]),vel_map->vbounds[vdim]);
  lbo->vparMaxSq = pow(lbo->vparMax,2);
  lbo->num_cbasis = cbasis->num_basis;

  lbo->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(lbo->eqn.flags);
  lbo->eqn.ref_count = gkyl_ref_count_init(gkyl_lbo_gyrokinetic_drag_free);

  // copy the host struct to device struct
  struct dg_lbo_gyrokinetic_drag *lbo_cu =
    (struct dg_lbo_gyrokinetic_drag*) gkyl_cu_malloc(sizeof(struct dg_lbo_gyrokinetic_drag));

  gkyl_cu_memcpy(lbo_cu, lbo,
    sizeof(struct dg_lbo_gyrokinetic_drag), GKYL_CU_MEMCPY_H2D);

  dg_lbo_gyrokinetic_drag_set_cu_dev_ptrs<<<1,1>>>(lbo_cu,
    cbasis->b_type, cv_index[cdim].vdim[vdim], cdim, vdim, poly_order);

  lbo->eqn.on_dev = &lbo_cu->eqn;  

  // Updater should store host pointers.
  lbo->gk_geom = geom_ho;
  lbo->vel_map = vel_map_ho;
  
  return &lbo->eqn;
}
