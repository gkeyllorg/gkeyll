#include "gkyl_dg_eqn.h"
#include <assert.h>
#include <stdio.h>
#include <float.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_dg_gyrokinetic.h>
#include <gkyl_dg_gyrokinetic_priv.h>
#include <gkyl_util.h>

void
gkyl_gyrokinetic_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_eqn *base = container_of(ref, struct gkyl_dg_eqn, ref_count);
  struct dg_gyrokinetic *gyrokinetic = container_of(base, struct dg_gyrokinetic, eqn);
  gkyl_gk_geometry_release(gyrokinetic->gk_geom);
  gkyl_velocity_map_release(gyrokinetic->vel_map);

  if (gkyl_dg_eqn_is_cu_dev(base)) {
    // free inner on_dev object
    struct dg_gyrokinetic *gyrokinetic_cu = container_of(base->on_dev, struct dg_gyrokinetic, eqn);
    gkyl_cu_free(gyrokinetic_cu);
  }

  gkyl_free(gyrokinetic);
}

void
gkyl_gyrokinetic_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_gyrokinetic_auxfields auxin)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_dg_eqn_is_cu_dev(eqn)) {
    gkyl_gyrokinetic_set_auxfields_cu(eqn->on_dev, auxin);
    return;
  }
#endif

  struct dg_gyrokinetic *gyrokinetic = container_of(eqn, struct dg_gyrokinetic, eqn);
  gyrokinetic->auxfields.flux_surf = auxin.flux_surf;
  gyrokinetic->auxfields.phi = auxin.phi;
  gyrokinetic->auxfields.apar = auxin.apar;
  gyrokinetic->auxfields.apardot = auxin.apardot;
}

struct gkyl_dg_eqn*
gkyl_dg_gyrokinetic_new(const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis,
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, 
  const double charge, const double mass, 
  enum gkyl_gk_collisionless_type collless_type, const bool no_by, const bool complete_em,
  const struct gk_geometry *gk_geom, const struct gkyl_velocity_map *vel_map, bool use_gpu)
{

#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    return gkyl_dg_gyrokinetic_cu_dev_new(cbasis, pbasis, conf_range, phase_range,
      charge, mass, collless_type, no_by, complete_em, gk_geom, vel_map);
#endif

  struct dg_gyrokinetic *gyrokinetic = gkyl_malloc(sizeof(struct dg_gyrokinetic));

  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  gyrokinetic->cdim = cdim;
  gyrokinetic->pdim = pdim;

  gyrokinetic->charge = charge;
  gyrokinetic->mass = mass;

  gyrokinetic->eqn.num_equations = 1;
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

  switch (cbasis->b_type) {
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
    gyrokinetic->vol_es_kernel = CK(vol_no_by_kernels,cdim,vdim,poly_order);
  }
  else {
    gyrokinetic->vol_es_kernel = CK(vol_kernels,cdim,vdim,poly_order);
  }

  // Setup electromagnetic terms if needed.
  bool is_em = (collless_type == GKYL_GK_COLLISIONLESS_EM);
  gyrokinetic->vol_add_apar_kernel = dg_gyrokinetic_add_apar_vol_none;
  gyrokinetic->vol_add_apardot_kernel = dg_gyrokinetic_add_apardot_vol_none;
  if (is_em) {
    if (complete_em) {
      // We complete with Apardot contribution to get the full GK RHS.
      gyrokinetic->vol_es_kernel = dg_gyrokinetic_vol_none;
      gyrokinetic->vol_add_apar_kernel = dg_gyrokinetic_add_apar_vol_none;
      gyrokinetic->vol_add_apardot_kernel = CK(vol_add_apardot_kernels,cdim,vdim,poly_order);
    } else {
      // We build for Ohm's law RHS (no Apardot contribution).
      gyrokinetic->vol_es_kernel = gyrokinetic->vol_es_kernel; // no change to ES kernel
      gyrokinetic->vol_add_apar_kernel = CK(vol_add_apar_kernels,cdim,vdim,poly_order);
      gyrokinetic->vol_add_apardot_kernel = dg_gyrokinetic_add_apardot_vol_none;
    }
  }

  gyrokinetic->surf[0] = CK(surf_x_kernels,cdim,vdim,poly_order);
  if (cdim>1)
    gyrokinetic->surf[1] = CK(surf_y_kernels,cdim,vdim,poly_order);
  if (cdim>2)
    gyrokinetic->surf[2] = CK(surf_z_kernels,cdim,vdim,poly_order);
  gyrokinetic->surf[cdim] = CK(surf_vpar_kernels,cdim,vdim,poly_order);

  gyrokinetic->boundary_surf[0] = CK(boundary_surf_x_kernels,cdim,vdim,poly_order);
  if (cdim>1)
    gyrokinetic->boundary_surf[1] = CK(boundary_surf_y_kernels,cdim,vdim,poly_order);
  if (cdim>2)
    gyrokinetic->boundary_surf[2] = CK(boundary_surf_z_kernels,cdim,vdim,poly_order);
  gyrokinetic->boundary_surf[cdim] = CK(boundary_surf_vpar_kernels,cdim,vdim,poly_order);

  // Ensure non-NULL pointers.
  for (int i=0; i<cdim; ++i) assert(gyrokinetic->surf[i]);
  assert(gyrokinetic->surf[cdim]);
  for (int i=0; i<cdim+1; ++i) assert(gyrokinetic->boundary_surf[i]);

  gyrokinetic->conf_range = *conf_range;
  gyrokinetic->phase_range = *phase_range;
  gyrokinetic->gk_geom = gkyl_gk_geometry_acquire(gk_geom);
  gyrokinetic->vel_map = gkyl_velocity_map_acquire(vel_map);
  gyrokinetic->auxfields.flux_surf = 0;
  gyrokinetic->auxfields.phi = 0;
  gyrokinetic->auxfields.apar = 0;
  gyrokinetic->auxfields.apardot = 0;

  gyrokinetic->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(gyrokinetic->eqn.flags);

  gyrokinetic->eqn.ref_count = gkyl_ref_count_init(gkyl_gyrokinetic_free);
  gyrokinetic->eqn.on_dev = &gyrokinetic->eqn; // CPU eqn obj points to itself

  return &gyrokinetic->eqn;
}
