/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_pkpm_kinetic_vel.h>
#include <gkyl_dg_pkpm_kinetic_vel_priv.h>
#include <gkyl_util.h>
}

// CUDA kernel to set device function pointers.
__global__ static void
dg_pkpm_kinetic_vel_set_cu_dev_ptrs(
  struct dg_pkpm_kinetic_vel *pkpm, enum gkyl_basis_type b_type, int cdim, int poly_order
)
{
  pkpm->eqn.vol_term = vol;
  pkpm->eqn.surf_term = surf;
  pkpm->eqn.boundary_surf_term = boundary_surf;
  choose_pkpm_kinetic_vel_kerns(pkpm, b_type, cdim, poly_order);
}

struct gkyl_dg_eqn *
gkyl_dg_pkpm_kinetic_vel_cu_dev_inew(const struct gkyl_dg_pkpm_kinetic_vel_inp *inp)
{
  struct dg_pkpm_kinetic_vel *pkpm =
    (struct dg_pkpm_kinetic_vel *)gkyl_malloc(sizeof(struct dg_pkpm_kinetic_vel));

  int cdim = inp->conf_basis->ndim, pdim = inp->phase_basis->ndim;
  pkpm->cdim = cdim;
  pkpm->pdim = pdim;
  pkpm->eqn.num_equations = 8;
  pkpm->conf_range = *inp->conf_range;
  pkpm->phase_range = *inp->phase_range;
  // Acquire the host arrays and store their on_dev pointers for copying over to device.
  struct gkyl_array *c_0_ho = gkyl_array_acquire(inp->c_0);
  struct gkyl_array *alpha_0_ho = gkyl_array_acquire(inp->alpha_0);
  struct gkyl_array *c_h_ho = gkyl_array_acquire(inp->c_h);
  struct gkyl_array *a_perp_ho = gkyl_array_acquire(inp->a_perp);
  struct gkyl_array *vel_flux_surf_ho = gkyl_array_acquire(inp->vel_flux_surf);
  pkpm->c_0 = c_0_ho->on_dev;
  pkpm->alpha_0 = alpha_0_ho->on_dev;
  pkpm->c_h = c_h_ho->on_dev;
  pkpm->a_perp = a_perp_ho->on_dev;
  pkpm->vel_flux_surf = vel_flux_surf_ho->on_dev;

  pkpm->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(pkpm->eqn.flags);
  pkpm->eqn.ref_count = gkyl_ref_count_init(gkyl_dg_pkpm_kinetic_vel_free);

  struct dg_pkpm_kinetic_vel *pkpm_cu =
    (struct dg_pkpm_kinetic_vel *)gkyl_cu_malloc(sizeof(struct dg_pkpm_kinetic_vel));
  gkyl_cu_memcpy(pkpm_cu, pkpm, sizeof(struct dg_pkpm_kinetic_vel), GKYL_CU_MEMCPY_H2D);

  dg_pkpm_kinetic_vel_set_cu_dev_ptrs<<<1, 1>>>(
    pkpm_cu, inp->conf_basis->b_type, cdim, inp->conf_basis->poly_order
  );

  pkpm->eqn.on_dev = &pkpm_cu->eqn;

  // Host-side equation object should store host pointers.
  pkpm->c_0 = c_0_ho;
  pkpm->alpha_0 = alpha_0_ho;
  pkpm->c_h = c_h_ho;
  pkpm->a_perp = a_perp_ho;
  pkpm->vel_flux_surf = vel_flux_surf_ho;
  return &pkpm->eqn;
}
