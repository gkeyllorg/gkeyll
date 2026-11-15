/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_pkpm_kinetic_conf.h>
#include <gkyl_dg_pkpm_kinetic_conf_priv.h>
#include <gkyl_util.h>
}

// CUDA kernel to set device function pointers.
__global__ static void
dg_pkpm_kinetic_conf_set_cu_dev_ptrs(
  struct dg_pkpm_kinetic_conf *pkpm, enum gkyl_basis_type b_type, int cdim, int poly_order
)
{
  pkpm->eqn.vol_term = vol;
  pkpm->eqn.surf_term = surf;
  pkpm->eqn.boundary_surf_term = boundary_surf;
  pkpm->vol_kern = choose_pkpm_kinetic_conf_vol_kern(b_type, cdim, poly_order);
  for (int d = 0; d < cdim; ++d) {
    pkpm->stream_surf[d] = choose_pkpm_stream_surf_kern(b_type, cdim, poly_order, d);
  }
}

struct gkyl_dg_eqn *
gkyl_dg_pkpm_kinetic_conf_cu_dev_inew(const struct gkyl_dg_pkpm_kinetic_conf_inp *inp)
{
  struct dg_pkpm_kinetic_conf *pkpm =
    (struct dg_pkpm_kinetic_conf *)gkyl_malloc(sizeof(struct dg_pkpm_kinetic_conf));

  int cdim = inp->conf_basis->ndim, pdim = inp->phase_basis->ndim;
  pkpm->cdim = cdim;
  pkpm->pdim = pdim;
  pkpm->num_surf_basis = inp->phase_basis->num_basis / 2;
  pkpm->eqn.num_equations = 8;
  pkpm->conf_range = *inp->conf_range;
  pkpm->phase_range = *inp->phase_range;
  // Acquire the host arrays and store their on_dev pointers for copying over to device.
  struct gkyl_array *bvar_ho = gkyl_array_acquire(inp->bvar);
  struct gkyl_array *u_ho = gkyl_array_acquire(inp->u);
  struct gkyl_array *xi_ho = gkyl_array_acquire(inp->xi);
  struct gkyl_array *conf_flux_surf_ho = gkyl_array_acquire(inp->conf_flux_surf);
  pkpm->bvar = bvar_ho->on_dev;
  pkpm->u = u_ho->on_dev;
  pkpm->xi = xi_ho->on_dev;
  pkpm->conf_flux_surf = conf_flux_surf_ho->on_dev;
  assert(inp->pos_map);
  pkpm->pos_map = 0; // host pointer set after the device memcpy below
  pkpm->jacob_pos = inp->pos_map->jacob_pos->on_dev;

  pkpm->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(pkpm->eqn.flags);
  pkpm->eqn.ref_count = gkyl_ref_count_init(gkyl_dg_pkpm_kinetic_conf_free);

  struct dg_pkpm_kinetic_conf *pkpm_cu =
    (struct dg_pkpm_kinetic_conf *)gkyl_cu_malloc(sizeof(struct dg_pkpm_kinetic_conf));
  gkyl_cu_memcpy(pkpm_cu, pkpm, sizeof(struct dg_pkpm_kinetic_conf), GKYL_CU_MEMCPY_H2D);

  dg_pkpm_kinetic_conf_set_cu_dev_ptrs<<<1, 1>>>(
    pkpm_cu, inp->conf_basis->b_type, cdim, inp->conf_basis->poly_order
  );

  pkpm->eqn.on_dev = &pkpm_cu->eqn;

  // Host-side equation object should store host pointers.
  pkpm->bvar = bvar_ho;
  pkpm->u = u_ho;
  pkpm->xi = xi_ho;
  pkpm->conf_flux_surf = conf_flux_surf_ho;
  pkpm->pos_map = gkyl_vlasov_position_map_acquire(inp->pos_map);
  pkpm->jacob_pos = inp->pos_map->jacob_pos;
  return &pkpm->eqn;
}
