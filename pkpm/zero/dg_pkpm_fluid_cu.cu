/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_pkpm_fluid.h>
#include <gkyl_dg_pkpm_fluid_priv.h>
#include <gkyl_util.h>
}

// CUDA kernel to set device function pointers.
__global__ static void
dg_pkpm_fluid_set_cu_dev_ptrs(
  struct dg_pkpm_fluid *pkpm, enum gkyl_basis_type b_type, int cdim, int poly_order
)
{
  pkpm->eqn.vol_term = vol;
  pkpm->eqn.surf_term = surf;
  pkpm->eqn.boundary_surf_term = boundary_surf;
  choose_pkpm_fluid_kerns(pkpm, b_type, cdim, poly_order);
}

struct gkyl_dg_eqn *
gkyl_dg_pkpm_fluid_cu_dev_inew(const struct gkyl_dg_pkpm_fluid_inp *inp)
{
  struct dg_pkpm_fluid *pkpm = (struct dg_pkpm_fluid *)gkyl_malloc(sizeof(struct dg_pkpm_fluid));

  int cdim = inp->conf_basis->ndim;
  pkpm->cdim = cdim;
  pkpm->num_surf_basis = inp->conf_basis->num_basis / 2;
  pkpm->eqn.num_equations = 5;
  pkpm->conf_range = *inp->conf_range;
  // Acquire the host arrays and store their on_dev pointers for copying over to device.
  struct gkyl_array *vlasov_pkpm_moms_ho = gkyl_array_acquire(inp->vlasov_pkpm_moms);
  struct gkyl_array *bvar_ho = gkyl_array_acquire(inp->bvar);
  struct gkyl_array *u_ho = gkyl_array_acquire(inp->u);
  struct gkyl_array *fluid_flux_surf_ho = gkyl_array_acquire(inp->fluid_flux_surf);
  pkpm->vlasov_pkpm_moms = vlasov_pkpm_moms_ho->on_dev;
  pkpm->bvar = bvar_ho->on_dev;
  pkpm->u = u_ho->on_dev;
  pkpm->fluid_flux_surf = fluid_flux_surf_ho->on_dev;
  assert(inp->pos_map);
  pkpm->pos_map = 0; // host pointer set after the device memcpy below
  pkpm->jacob_pos = inp->pos_map->jacob_pos->on_dev;

  pkpm->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(pkpm->eqn.flags);
  pkpm->eqn.ref_count = gkyl_ref_count_init(gkyl_dg_pkpm_fluid_free);

  struct dg_pkpm_fluid *pkpm_cu = (struct dg_pkpm_fluid *)gkyl_cu_malloc(sizeof(struct dg_pkpm_fluid));
  gkyl_cu_memcpy(pkpm_cu, pkpm, sizeof(struct dg_pkpm_fluid), GKYL_CU_MEMCPY_H2D);

  dg_pkpm_fluid_set_cu_dev_ptrs<<<1, 1>>>(
    pkpm_cu, inp->conf_basis->b_type, cdim, inp->conf_basis->poly_order
  );

  pkpm->eqn.on_dev = &pkpm_cu->eqn;

  // Host-side equation object should store host pointers.
  pkpm->vlasov_pkpm_moms = vlasov_pkpm_moms_ho;
  pkpm->bvar = bvar_ho;
  pkpm->u = u_ho;
  pkpm->fluid_flux_surf = fluid_flux_surf_ho;
  pkpm->pos_map = gkyl_vlasov_position_map_acquire(inp->pos_map);
  pkpm->jacob_pos = inp->pos_map->jacob_pos;
  return &pkpm->eqn;
}
