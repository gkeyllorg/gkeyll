/* -*- c++ -*- */
#include <float.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_gyrokinetic_passive.h>
#include <gkyl_dg_gyrokinetic_passive_priv.h>
}

#include <cassert>

// CUDA kernel to set pointer to auxiliary fields.
// Required because the equation object lives on the device and its members
// cannot be modified without a __global__ kernel.
__global__ static void
gkyl_gyrokinetic_passive_set_auxfields_cu_kernel(const struct gkyl_dg_eqn *eqn,
  const struct gkyl_array *flux_surf, const struct gkyl_array *speeds)
{
  struct dg_gyrokinetic_passive *gkp = container_of(eqn, struct dg_gyrokinetic_passive, eqn);
  gkp->auxfields.flux_surf = flux_surf;
  gkp->auxfields.speeds    = speeds;
}

// Host-side wrapper for gkyl_gyrokinetic_passive_set_auxfields_cu_kernel.
void
gkyl_gyrokinetic_passive_set_auxfields_cu(const struct gkyl_dg_eqn *eqn,
  struct gkyl_dg_gyrokinetic_passive_auxfields auxin)
{
  gkyl_gyrokinetic_passive_set_auxfields_cu_kernel<<<1,1>>>(eqn,
    auxin.flux_surf->on_dev, auxin.speeds->on_dev);
}

// CUDA kernel to set device function pointers and zero auxfields.
// Doing function-pointer work here avoids troublesome cudaMemcpyFromSymbol.
__global__ static void
dg_gyrokinetic_passive_set_cu_dev_ptrs(struct dg_gyrokinetic_passive *gkp,
  enum gkyl_basis_type b_type, int cv_index, int cdim, int vdim, int poly_order)
{
  gkp->auxfields.flux_surf = 0;
  gkp->auxfields.speeds    = 0;

  gkp->eqn.surf_term          = surf_passive;
  gkp->eqn.boundary_surf_term = boundary_surf_passive;
  gkp->eqn.boundary_diag_term = boundary_diag_passive;

  const gkyl_dg_gyrokinetic_passive_vol_kern_list           *vol_kernels;
  const gkyl_dg_gyrokinetic_passive_surf_kern_list          *surf_x_kernels;
  const gkyl_dg_gyrokinetic_passive_surf_kern_list          *surf_y_kernels;
  const gkyl_dg_gyrokinetic_passive_surf_kern_list          *surf_z_kernels;
  const gkyl_dg_gyrokinetic_passive_boundary_surf_kern_list *boundary_surf_x_kernels;
  const gkyl_dg_gyrokinetic_passive_boundary_surf_kern_list *boundary_surf_y_kernels;
  const gkyl_dg_gyrokinetic_passive_boundary_surf_kern_list *boundary_surf_z_kernels;

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      vol_kernels             = ser_passive_vol_kernels;
      surf_x_kernels          = ser_passive_surf_x_kernels;
      surf_y_kernels          = ser_passive_surf_y_kernels;
      surf_z_kernels          = ser_passive_surf_z_kernels;
      boundary_surf_x_kernels = ser_passive_boundary_surf_x_kernels;
      boundary_surf_y_kernels = ser_passive_boundary_surf_y_kernels;
      boundary_surf_z_kernels = ser_passive_boundary_surf_z_kernels;
      break;

    default:
      assert(false);
      break;
  }

  gkp->eqn.vol_term = vol_kernels[cv_index].kernels[poly_order-1];

  gkp->surf[0]            = surf_x_kernels[cv_index].kernels[poly_order-1];
  gkp->boundary_surf[0]   = boundary_surf_x_kernels[cv_index].kernels[poly_order-1];
  if (cdim > 1) {
    gkp->surf[1]          = surf_y_kernels[cv_index].kernels[poly_order-1];
    gkp->boundary_surf[1] = boundary_surf_y_kernels[cv_index].kernels[poly_order-1];
  }
  if (cdim > 2) {
    gkp->surf[2]          = surf_z_kernels[cv_index].kernels[poly_order-1];
    gkp->boundary_surf[2] = boundary_surf_z_kernels[cv_index].kernels[poly_order-1];
  }
}

struct gkyl_dg_eqn*
gkyl_dg_gyrokinetic_passive_cu_dev_new(const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis,
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const double charge, const double mass,
  const struct gk_geometry *gk_geom, const struct gkyl_velocity_map *vel_map)
{
  struct dg_gyrokinetic_passive *gkp = (struct dg_gyrokinetic_passive*) gkyl_malloc(sizeof(*gkp));

  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim - cdim;
  int poly_order = cbasis->poly_order;

  gkp->cdim = cdim;
  gkp->pdim = pdim;

  gkp->charge = charge;
  gkp->mass   = mass;

  gkp->eqn.num_equations = 1;

  // Acquire on_dev pointers so the memcpy below copies device references.
  struct gk_geometry       *geom_ho    = gkyl_gk_geometry_acquire(gk_geom);
  struct gkyl_velocity_map *vel_map_ho = gkyl_velocity_map_acquire(vel_map);

  gkp->gk_geom  = geom_ho->on_dev;
  gkp->vel_map  = vel_map_ho->on_dev;

  gkp->conf_range  = *conf_range;
  gkp->phase_range = *phase_range;

  gkp->eqn.flags     = 0;
  GKYL_SET_CU_ALLOC(gkp->eqn.flags);
  gkp->eqn.ref_count = gkyl_ref_count_init(gkyl_gyrokinetic_passive_free);

  // Copy host struct to device.
  struct dg_gyrokinetic_passive *gkp_cu =
    (struct dg_gyrokinetic_passive*) gkyl_cu_malloc(sizeof(struct dg_gyrokinetic_passive));
  gkyl_cu_memcpy(gkp_cu, gkp, sizeof(struct dg_gyrokinetic_passive), GKYL_CU_MEMCPY_H2D);

  dg_gyrokinetic_passive_set_cu_dev_ptrs<<<1,1>>>(gkp_cu,
    cbasis->b_type, cv_index_gkp[cdim].vdim[vdim], cdim, vdim, poly_order);

  // Set on_dev pointer to the device struct.
  gkp->eqn.on_dev = &gkp_cu->eqn;

  // Restore host-side geometry pointers for the updater.
  gkp->gk_geom = geom_ho;
  gkp->vel_map = vel_map_ho;

  return &gkp->eqn;
}
