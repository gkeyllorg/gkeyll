// Private header: not for direct use
#pragma once

#include <assert.h>

#include <gkyl_array.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_pkpm_fluid.h>
#include <gkyl_pkpm_fluid_kernels.h>
#include <gkyl_pkpm_layout.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

typedef double (*pkpm_fluid_vol_t)(
  const double *w, const double *dxv, const double *jacob_pos, const double *vlasov_pkpm_moms,
  const double *bvar, const double *u, const double *fluid, double *GKYL_RESTRICT out
);

typedef double (*pkpm_fluid_surf_t)(
  const double *w, const double *dxv, const double *flux_l, const double *flux_r,
  double *GKYL_RESTRICT out
);

// for use in kernel tables
typedef struct {
  pkpm_fluid_surf_t kernels[3];
} gkyl_dg_pkpm_fluid_surf_kern_list;

struct dg_pkpm_fluid {
  struct gkyl_dg_eqn eqn; // Base object
  int cdim; // Configuration-space dimensions
  int num_surf_basis; // Surface modal basis size per component and direction
  pkpm_fluid_vol_t vol_kern; // Volume kernel
  pkpm_fluid_surf_t fluid_surf[3]; // Surface kernels scattering stored fluxes
  struct gkyl_range conf_range; // Configuration-space range
  const struct gkyl_vlasov_position_map
    *pos_map; // Position map (acquired; host pointer, 0 on the device copy)
  const struct gkyl_array
    *jacob_pos; // Per-direction Jacobians (borrowed from pos_map; device pointer on device)
  struct gkyl_dg_pkpm_fluid_auxfields auxfields; // Auxiliary fields.
};

// Kernel tables (p=1 serendipity only)
GKYL_CU_D static const pkpm_fluid_vol_t ser_pkpm_fluid_vol_kernels[] = {
  pkpm_fluid_vol_1x_ser_p1, pkpm_fluid_vol_2x_ser_p1, pkpm_fluid_vol_3x_ser_p1
};

GKYL_CU_D static const gkyl_dg_pkpm_fluid_surf_kern_list ser_pkpm_fluid_surf_kernels[] = {
  {pkpm_fluid_surfx_1x_ser_p1, NULL, NULL}, // 0
  {pkpm_fluid_surfx_2x_ser_p1, pkpm_fluid_surfy_2x_ser_p1, NULL}, // 1
  {pkpm_fluid_surfx_3x_ser_p1, pkpm_fluid_surfy_3x_ser_p1, pkpm_fluid_surfz_3x_ser_p1} // 2
};

GKYL_CU_D static void
choose_pkpm_fluid_kerns(
  struct dg_pkpm_fluid *pkpm, enum gkyl_basis_type b_type, int cdim, int poly_order
)
{
  assert(poly_order == 1);
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      pkpm->vol_kern = ser_pkpm_fluid_vol_kernels[cdim - 1];
      for (int d = 0; d < cdim; ++d) {
        pkpm->fluid_surf[d] = ser_pkpm_fluid_surf_kernels[cdim - 1].kernels[d];
      }
      break;
    default:
      assert(false);
      break;
  }
}

GKYL_CU_D static double
vol(
  const struct gkyl_dg_eqn *eqn, const double *xc, const double *dx, const int *idx,
  const double *qIn, double *GKYL_RESTRICT qRhsOut
)
{
  struct dg_pkpm_fluid *pkpm = container_of(eqn, struct dg_pkpm_fluid, eqn);
  long cidx = gkyl_range_idx(&pkpm->conf_range, idx);
  return pkpm->vol_kern(
    xc, dx, (const double *)gkyl_array_cfetch(pkpm->jacob_pos, cidx),
    (const double *)gkyl_array_cfetch(pkpm->auxfields.vlasov_pkpm_moms, cidx),
    (const double *)gkyl_array_cfetch(pkpm->auxfields.bvar, cidx),
    (const double *)gkyl_array_cfetch(pkpm->auxfields.u, cidx), qIn, qRhsOut
  );
}

// Each cell owns its lower face; the upper face is the lower face of the right
// neighbor (an interior cell or the upper ghost cell).
GKYL_CU_D static double
surf(
  const struct gkyl_dg_eqn *eqn, int dir, const double *xcL, const double *xcC, const double *xcR,
  const double *dxL, const double *dxC, const double *dxR, const int *idxL, const int *idxC,
  const int *idxR, const double *qInL, const double *qInC, const double *qInR,
  double *GKYL_RESTRICT qRhsOut
)
{
  struct dg_pkpm_fluid *pkpm = container_of(eqn, struct dg_pkpm_fluid, eqn);
  int off = dir * PKPM_FLUID_NUM_COMP * pkpm->num_surf_basis;
  long cidxC = gkyl_range_idx(&pkpm->conf_range, idxC);
  long cidxR = gkyl_range_idx(&pkpm->conf_range, idxR);
  const double *flux_l = (const double *)gkyl_array_cfetch(pkpm->auxfields.fluid_flux_surf, cidxC);
  const double *flux_r = (const double *)gkyl_array_cfetch(pkpm->auxfields.fluid_flux_surf, cidxR);
  return pkpm->fluid_surf[dir](xcC, dxC, &flux_l[off], &flux_r[off], qRhsOut);
}

// The momentum equation has ghost cells in every direction; no zero-flux faces.
GKYL_CU_D static double
boundary_surf(
  const struct gkyl_dg_eqn *eqn, int dir, const double *xcEdge, const double *xcSkin,
  const double *dxEdge, const double *dxSkin, const int *idxEdge, const int *idxSkin,
  const int edge, const double *qInEdge, const double *qInSkin, double *GKYL_RESTRICT qRhsOut
)
{
  return 0.0;
}
