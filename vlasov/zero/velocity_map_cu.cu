/* -*- c++ -*- */
extern "C" {
#include <gkyl_velocity_map.h>
#include <gkyl_velocity_map_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_alloc_flags_priv.h>
}

struct gkyl_velocity_map* gkyl_velocity_map_new_cu_dev(struct gkyl_velocity_map *gvm_ho)
{
  struct gkyl_velocity_map *gvm = (struct gkyl_velocity_map *) gkyl_malloc(sizeof(*gvm));

  gvm->is_identity   = gvm_ho->is_identity;
  gvm->grid          = gvm_ho->grid;
  gvm->grid_vel      = gvm_ho->grid_vel;
  gvm->local         = gvm_ho->local;
  gvm->local_ext     = gvm_ho->local_ext;
  gvm->local_vel     = gvm_ho->local_vel;
  gvm->local_ext_vel = gvm_ho->local_ext_vel;
  gvm->vmap_basis_ho = gvm_ho->vmap_basis_ho;
  memcpy(gvm->vbounds, gvm_ho->vbounds, sizeof(double[2*GKYL_MAX_VDIM]));

  // Copy the host-side initialized object to the device.
  struct gkyl_array *vmap       = mkarr(true, gvm_ho->vmap      ->ncomp, gvm_ho->vmap      ->size);
  struct gkyl_array *vmap_sq    = mkarr(true, gvm_ho->vmap_sq   ->ncomp, gvm_ho->vmap_sq   ->size);
  struct gkyl_array *vmap_prime = mkarr(true, gvm_ho->vmap_prime->ncomp, gvm_ho->vmap_prime->size);
  struct gkyl_array *jacobvel   = mkarr(true, gvm_ho->jacobvel  ->ncomp, gvm_ho->jacobvel  ->size);
  struct gkyl_basis *vmap_basis = gkyl_cart_modal_serendip_cu_dev_new(gvm_ho->vmap_basis->ndim,
    gvm_ho->vmap_basis->poly_order);
  // Need a host copy of vmap for some IC projection options.
  struct gkyl_array *vmap_ho    = mkarr(false, gvm_ho->vmap      ->ncomp, gvm_ho->vmap      ->size);

  gkyl_array_copy(vmap      , gvm_ho->vmap      );
  gkyl_array_copy(vmap_sq   , gvm_ho->vmap_sq   );
  gkyl_array_copy(vmap_prime, gvm_ho->vmap_prime);
  gkyl_array_copy(jacobvel  , gvm_ho->jacobvel  );
  gkyl_array_copy(vmap_ho   , gvm_ho->vmap_ho   );

  gvm->vmap       = vmap      ->on_dev;
  gvm->vmap_sq    = vmap_sq   ->on_dev;
  gvm->vmap_prime = vmap_prime->on_dev;
  gvm->jacobvel   = jacobvel  ->on_dev;
  gvm->vmap_ho    = vmap      ->on_dev; // MF 2024/05/06: I think this is safer.
  gvm->vmap_basis = vmap_basis;

  gvm->flags = 0;
  GKYL_SET_CU_ALLOC(gvm->flags);
  gvm->ref_count = gkyl_ref_count_init(gkyl_velocity_map_free);

  // Initialize the device object.
  struct gkyl_velocity_map *gvm_cu = (struct gkyl_velocity_map*) gkyl_cu_malloc(sizeof(*gvm_cu));
  gkyl_cu_memcpy(gvm_cu, gvm, sizeof(struct gkyl_velocity_map), GKYL_CU_MEMCPY_H2D);
  gvm->on_dev = gvm_cu;

  // The returned object should store host pointers to gkyl_arrays.
  gvm->vmap       = vmap      ;
  gvm->vmap_sq    = vmap_sq   ;
  gvm->vmap_prime = vmap_prime;
  gvm->jacobvel   = jacobvel  ;
  gvm->vmap_ho    = vmap_ho   ;

  return gvm;
}
