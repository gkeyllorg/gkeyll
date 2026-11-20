// Private header: not for direct use
#pragma once

#include <assert.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_pkpm_harm_slots.h>
#include <gkyl_pkpm_kinetic_kernels.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

// Largest number of w-face nodes (5^cdim).
#define PKPM_VF_MAXNODE 125

typedef void (*pkpm_vel_flux_nodes_t)(
  const double *wn, const double *theta, const double *bvar, const double *alpha_0, const double *c_h,
  double *GKYL_RESTRICT nv
);
typedef void (*pkpm_vel_flux_t)(
  const double *w, const double *dxv, const double *nv, const double *f_l, const double *f_r,
  double *GKYL_RESTRICT flux
);

struct gkyl_dg_pkpm_vel_flux_surf {
  struct gkyl_rect_grid phase_grid; // Phase-space grid.
  int cdim; // Configuration-space dimensions.
  int pdim; // Phase-space dimensions.
  int num_surf_basis; // Surface modal basis size per component.
  bool use_gpu; // bool to determine if on GPU.
  pkpm_vel_flux_nodes_t vel_nodes; // Node values of the w flux of a configuration cell.
  pkpm_vel_flux_t vel_flux; // Flux kernel on the v_par face.

  uint32_t flags;
  struct gkyl_dg_pkpm_vel_flux_surf *on_dev; // pointer to itself or device data
};

// Sizes of the updater (host function shared by the host and device constructors).
void gkyl_dg_pkpm_vel_flux_surf_setup(
  struct gkyl_dg_pkpm_vel_flux_surf *up, const struct gkyl_dg_pkpm_vel_flux_surf_inp *inp
);

// Node-value and v_par flux kernels (p=1 serendipity only; indexed by cdim-1)
GKYL_CU_D static const pkpm_vel_flux_nodes_t ser_pkpm_vel_flux_nodes_kernels[] = {
  pkpm_vel_flux_nodes_1x_ser_p1,
  pkpm_vel_flux_nodes_2x_ser_p1,
  pkpm_vel_flux_nodes_3x_ser_p1,
};
GKYL_CU_D static const pkpm_vel_flux_t ser_pkpm_vel_flux_kernels[] = {
  pkpm_vel_flux_vpar_1x1v_ser_p1,
  pkpm_vel_flux_vpar_2x1v_ser_p1,
  pkpm_vel_flux_vpar_3x1v_ser_p1,
};

GKYL_CU_D
static pkpm_vel_flux_t
choose_pkpm_vel_flux_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  assert(poly_order == 1);
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_pkpm_vel_flux_kernels[cdim - 1];
      break;
    default:
      assert(false);
      break;
  }
  return 0;
}
