/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_euler.h>
#include <gkyl_dg_euler_priv.h>
#include <gkyl_util.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_wv_euler.h>
#include <gkyl_wv_iso_euler.h>
}

#include <cassert>

#define CK(lst, cdim, poly_order) lst[cdim - 1].kernels[poly_order]

__global__ void static dg_euler_set_cu_dev_ptrs(
  struct dg_euler *euler, enum gkyl_basis_type b_type, int cdim, int poly_order
)
{
  const gkyl_dg_euler_vol_kern_list *vol_kernels;
  const gkyl_dg_euler_surf_kern_list *surf_kernels[3];
  const gkyl_dg_euler_wall_flux_kern_list *wall_flux_kernels[3];
  const gkyl_dg_euler_boundary_surf_kern_list *boundary_surf_kernels[3];

  bool is_iso = euler->eqn_type == GKYL_EQN_ISO_EULER;
  switch (b_type) {
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

  euler->eqn.vol_term = vol;
  euler->eqn.surf_term = surf;
  euler->eqn.boundary_surf_term = boundary_surf;

  euler->vol = CK(vol_kernels, cdim, poly_order);
  for (int d = 0; d < cdim; ++d) {
    euler->surf[d] = CK(surf_kernels[d], cdim, poly_order);
    euler->wall_flux[d] = CK(wall_flux_kernels[d], cdim, poly_order);
    euler->boundary_surf[d] = CK(boundary_surf_kernels[d], cdim, poly_order);
  }
}

struct gkyl_dg_eqn *
gkyl_dg_euler_cu_dev_inew(const struct gkyl_dg_euler_inp *inp)
{
  struct dg_euler *euler = (struct dg_euler *)gkyl_malloc(sizeof(struct dg_euler));

  const struct gkyl_wv_eqn *wv_eqn = inp->wv_eqn;
  int cdim = inp->cbasis->ndim, poly_order = inp->cbasis->poly_order;

  euler->eqn_type = wv_eqn->type;
  euler->eqn.num_equations = wv_eqn->num_equations;
  euler->gas_gamma = 0.0;
  euler->vt = 0.0;
  if (euler->eqn_type == GKYL_EQN_ISO_EULER) {
    euler->vt = gkyl_wv_iso_euler_vt(wv_eqn);
    euler->vol_param = euler->vt;
  } else {
    euler->gas_gamma = gkyl_wv_euler_gas_gamma(wv_eqn);
    euler->vol_param = euler->gas_gamma;
  }
  euler->cdim = cdim;
  euler->num_surf = cdim == 1 ? 1 : inp->cbasis->num_basis / (poly_order + 1);

  // acquire pointer to wave equation object
  struct gkyl_wv_eqn *eqn = gkyl_wv_eqn_acquire(wv_eqn);
  euler->wv_eqn = eqn->on_dev; // this is so the memcpy below has wv_eqn on_dev

  // acquire pointer to wave geometry object
  struct gkyl_wave_geom *geom = gkyl_wave_geom_acquire(inp->geom);
  euler->geom = geom->on_dev; // this is so the memcpy below has geom on_dev

  // Arrays (acquired; device pointers in the device struct).
  assert(inp->u && inp->p && inp->pos_map);
  struct gkyl_array *u_ho = gkyl_array_acquire(inp->u);
  struct gkyl_array *p_ho = gkyl_array_acquire(inp->p);
  euler->u = u_ho->on_dev;
  euler->p = p_ho->on_dev;
  euler->pos_map = 0; // Host pointer set after the device memcpy below.
  euler->jacob_pos = inp->pos_map->jacob_pos->on_dev;
  euler->jacob_stride = poly_order + 1;

  euler->conf_range = *inp->conf_range;

  euler->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(euler->eqn.flags);
  euler->eqn.ref_count = gkyl_ref_count_init(gkyl_dg_euler_free);

  // copy the host struct to device struct
  struct dg_euler *euler_cu = (struct dg_euler *)gkyl_cu_malloc(sizeof(struct dg_euler));
  gkyl_cu_memcpy(euler_cu, euler, sizeof(struct dg_euler), GKYL_CU_MEMCPY_H2D);
  dg_euler_set_cu_dev_ptrs<<<1, 1>>>(euler_cu, inp->cbasis->b_type, cdim, poly_order);

  // set parent on_dev pointer
  euler->eqn.on_dev = &euler_cu->eqn;

  // updater should store host pointers
  euler->wv_eqn = eqn;
  euler->geom = geom;
  euler->u = u_ho;
  euler->p = p_ho;
  euler->pos_map = gkyl_vlasov_position_map_acquire(inp->pos_map);
  euler->jacob_pos = inp->pos_map->jacob_pos;

  return &euler->eqn;
}

__global__ static void
gkyl_dg_euler_wall_flux_cu_kernel(
  const struct gkyl_dg_eqn *eqn, struct gkyl_rect_grid grid, int dir, int edge,
  struct gkyl_range skin_range, const struct gkyl_array *fluid, struct gkyl_array *flux_int
)
{
  const struct dg_euler *euler = container_of(eqn, struct dg_euler, eqn);
  int idx[GKYL_MAX_DIM];
  for (unsigned long linc = threadIdx.x + blockIdx.x * blockDim.x; linc < skin_range.volume;
       linc += blockDim.x * gridDim.x) {
    gkyl_sub_range_inv_idx(&skin_range, linc, idx);
    long loc = gkyl_range_idx(&skin_range, idx);
    dg_euler_wall_flux_cell(
      euler, dir, edge, idx, grid.dx, (const double *)gkyl_array_cfetch(fluid, loc),
      (double *)gkyl_array_fetch(flux_int, loc)
    );
  }
}

void
gkyl_dg_euler_wall_flux_cu(
  const struct gkyl_dg_eqn *eqn, const struct gkyl_rect_grid *grid, int dir, int edge,
  const struct gkyl_range *skin_range, const struct gkyl_array *fluid, struct gkyl_array *flux_int
)
{
  int nblocks = skin_range->nblocks, nthreads = skin_range->nthreads;
  gkyl_dg_euler_wall_flux_cu_kernel<<<nblocks, nthreads>>>(
    eqn->on_dev, *grid, dir, edge, *skin_range, fluid->on_dev, flux_int->on_dev
  );
}
