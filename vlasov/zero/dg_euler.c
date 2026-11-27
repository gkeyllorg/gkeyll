#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_euler.h>
#include <gkyl_dg_euler_priv.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_wv_euler.h>
#include <gkyl_wv_iso_euler.h>
#include <gkyl_util.h>

// "Choose Kernel" based on cdim and polyorder
#define CK(lst, cdim, poly_order) lst[cdim - 1].kernels[poly_order]

void
gkyl_dg_euler_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_eqn *base = container_of(ref, struct gkyl_dg_eqn, ref_count);
  struct dg_euler *euler = container_of(base, struct dg_euler, eqn);
  gkyl_wv_eqn_release(euler->wv_eqn);
  gkyl_wave_geom_release(euler->geom);
  gkyl_array_release(euler->u);
  gkyl_array_release(euler->p);
  gkyl_vlasov_position_map_release(euler->pos_map);

  if (gkyl_dg_eqn_is_cu_dev(base)) {
    // free inner on_dev object
    struct dg_euler *euler_cu = container_of(base->on_dev, struct dg_euler, eqn);
    gkyl_cu_free(euler_cu);
  }
  gkyl_free(euler);
}

struct gkyl_dg_eqn *
gkyl_dg_euler_inew(const struct gkyl_dg_euler_inp *inp)
{
#ifdef GKYL_HAVE_CUDA
  if (inp->use_gpu) {
    return gkyl_dg_euler_cu_dev_inew(inp);
  }
#endif
  struct dg_euler *euler = gkyl_malloc(sizeof(struct dg_euler));

  const struct gkyl_basis *cbasis = inp->cbasis;
  const struct gkyl_wv_eqn *wv_eqn = inp->wv_eqn;
  int cdim = cbasis->ndim;
  int poly_order = cbasis->poly_order;

  const gkyl_dg_euler_vol_kern_list *vol_kernels;
  const gkyl_dg_euler_surf_kern_list *surf_kernels[3];
  const gkyl_dg_euler_wall_flux_kern_list *wall_flux_kernels[3];
  const gkyl_dg_euler_boundary_surf_kern_list *boundary_surf_kernels[3];

  bool is_iso = wv_eqn->type == GKYL_EQN_ISO_EULER;
  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      vol_kernels = is_iso ? ser_iso_vol_kernels : ser_vol_kernels;
      surf_kernels[0] = is_iso ? ser_iso_surf_x_kernels : ser_surf_x_kernels;
      surf_kernels[1] = is_iso ? ser_iso_surf_y_kernels : ser_surf_y_kernels;
      surf_kernels[2] = is_iso ? ser_iso_surf_z_kernels : ser_surf_z_kernels;
      wall_flux_kernels[0] = is_iso ? ser_iso_wall_flux_x_kernels : ser_wall_flux_x_kernels;
      wall_flux_kernels[1] = is_iso ? ser_iso_wall_flux_y_kernels : ser_wall_flux_y_kernels;
      wall_flux_kernels[2] = is_iso ? ser_iso_wall_flux_z_kernels : ser_wall_flux_z_kernels;
      boundary_surf_kernels[0] = is_iso ? ser_iso_boundary_surf_x_kernels :
                                          ser_boundary_surf_x_kernels;
      boundary_surf_kernels[1] = is_iso ? ser_iso_boundary_surf_y_kernels :
                                          ser_boundary_surf_y_kernels;
      boundary_surf_kernels[2] = is_iso ? ser_iso_boundary_surf_z_kernels :
                                          ser_boundary_surf_z_kernels;
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      vol_kernels = is_iso ? ten_iso_vol_kernels : ten_vol_kernels;
      surf_kernels[0] = is_iso ? ten_iso_surf_x_kernels : ten_surf_x_kernels;
      surf_kernels[1] = is_iso ? ten_iso_surf_y_kernels : ten_surf_y_kernels;
      surf_kernels[2] = is_iso ? ten_iso_surf_z_kernels : ten_surf_z_kernels;
      wall_flux_kernels[0] = is_iso ? ten_iso_wall_flux_x_kernels : ten_wall_flux_x_kernels;
      wall_flux_kernels[1] = is_iso ? ten_iso_wall_flux_y_kernels : ten_wall_flux_y_kernels;
      wall_flux_kernels[2] = is_iso ? ten_iso_wall_flux_z_kernels : ten_wall_flux_z_kernels;
      boundary_surf_kernels[0] = is_iso ? ten_iso_boundary_surf_x_kernels :
                                          ten_boundary_surf_x_kernels;
      boundary_surf_kernels[1] = is_iso ? ten_iso_boundary_surf_y_kernels :
                                          ten_boundary_surf_y_kernels;
      boundary_surf_kernels[2] = is_iso ? ten_iso_boundary_surf_z_kernels :
                                          ten_boundary_surf_z_kernels;
      break;

    default:
      assert(false);
      break;
  }

  euler->eqn_type = wv_eqn->type;
  euler->eqn.num_equations = wv_eqn->num_equations;
  euler->wv_eqn = gkyl_wv_eqn_acquire(wv_eqn);
  euler->geom = gkyl_wave_geom_acquire(inp->geom);
  euler->gas_gamma = 0.0;
  euler->vt = 0.0;
  if (is_iso) {
    euler->vt = gkyl_wv_iso_euler_vt(euler->wv_eqn);
    euler->vol_param = euler->vt;
  } else {
    euler->gas_gamma = gkyl_wv_euler_gas_gamma(euler->wv_eqn);
    euler->vol_param = euler->gas_gamma;
  }
  euler->cdim = cdim;
  euler->num_surf = cdim == 1 ? 1 : cbasis->num_basis / (poly_order + 1);

  // The position map: on a non-identity map each cell has its physical size.
  assert(inp->pos_map);
  euler->pos_map = gkyl_vlasov_position_map_acquire(inp->pos_map);
  euler->jacob_pos = inp->pos_map->jacob_pos; // Borrowed; kept alive by the acquired pos_map.
  euler->jacob_stride = poly_order + 1;
  bool is_mapped = !inp->pos_map->is_identity;

  euler->eqn.vol_term = is_mapped ? vol_mapped : vol;
  euler->eqn.surf_term = is_mapped ? surf_mapped : surf;
  euler->eqn.boundary_surf_term = is_mapped ? boundary_surf_mapped : boundary_surf;

  euler->vol = CK(vol_kernels, cdim, poly_order);
  assert(euler->vol);
  for (int d = 0; d < cdim; ++d) {
    euler->surf[d] = CK(surf_kernels[d], cdim, poly_order);
    euler->wall_flux[d] = CK(wall_flux_kernels[d], cdim, poly_order);
    euler->boundary_surf[d] = CK(boundary_surf_kernels[d], cdim, poly_order);
    assert(euler->surf[d] && euler->wall_flux[d] && euler->boundary_surf[d]);
  }

  assert(inp->u && inp->p);
  euler->u = gkyl_array_acquire(inp->u);
  euler->p = gkyl_array_acquire(inp->p);
  euler->conf_range = *inp->conf_range;

  euler->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(euler->eqn.flags);
  euler->eqn.ref_count = gkyl_ref_count_init(gkyl_dg_euler_free);
  euler->eqn.on_dev = &euler->eqn; // CPU eqn obj points to itself

  return &euler->eqn;
}

void
gkyl_dg_euler_wall_flux(
  const struct gkyl_dg_eqn *eqn, const struct gkyl_rect_grid *grid, int dir, int edge,
  const struct gkyl_range *skin_range, const struct gkyl_array *fluid, struct gkyl_array *flux_int
)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(flux_int)) {
    gkyl_dg_euler_wall_flux_cu(eqn, grid, dir, edge, skin_range, fluid, flux_int);
    return;
  }
#endif
  const struct dg_euler *euler = container_of(eqn, struct dg_euler, eqn);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, skin_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(skin_range, iter.idx);
    dg_euler_wall_flux_cell(
      euler, dir, edge, iter.idx, grid->dx, (const double *)gkyl_array_cfetch(fluid, loc),
      (double *)gkyl_array_fetch(flux_int, loc)
    );
  }
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_eqn *
gkyl_dg_euler_cu_dev_inew(const struct gkyl_dg_euler_inp *inp)
{
  assert(false);
  return 0;
}

void
gkyl_dg_euler_wall_flux_cu(
  const struct gkyl_dg_eqn *eqn, const struct gkyl_rect_grid *grid, int dir, int edge,
  const struct gkyl_range *skin_range, const struct gkyl_array *fluid, struct gkyl_array *flux_int
)
{
  assert(false);
}

#endif
