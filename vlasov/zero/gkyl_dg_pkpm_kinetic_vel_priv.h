// Private header: not for direct use
#pragma once

#include <assert.h>

#include <gkyl_array.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_pkpm_kinetic_vel.h>
#include <gkyl_pkpm_kinetic_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

typedef double (*pkpm_kinetic_vel_vol_t)(
  const double *w, const double *dxv, const double *c_0, const double *alpha_0, const double *c_h,
  const double *a_perp, const double *f, double *GKYL_RESTRICT out
);

typedef double (*pkpm_accel_surf_t)(
  const double *w, const double *dxv, const double *flux_l, const double *flux_r,
  double *GKYL_RESTRICT out
);

struct dg_pkpm_kinetic_vel {
  struct gkyl_dg_eqn eqn; // Base object
  int cdim; // Configuration-space dimensions
  int pdim; // Phase-space dimensions
  pkpm_kinetic_vel_vol_t vol_kern; // Volume kernel (closure terms)
  pkpm_accel_surf_t accel_surf; // Surface kernel scattering stored v_par fluxes
  struct gkyl_range conf_range; // Configuration-space range
  struct gkyl_range phase_range; // Phase-space range
  // Arrays the terms read (acquired; device pointers on the device copy).
  const struct gkyl_array *c_0;
  const struct gkyl_array *alpha_0;
  const struct gkyl_array *c_h;
  const struct gkyl_array *a_perp;
  const struct gkyl_array *vel_flux_surf;
};

/**
 * Free the equation object of the PKPM velocity-pass kinetic equations.
 *
 * @param ref Reference counter of the equation object
 */
void gkyl_dg_pkpm_kinetic_vel_free(const struct gkyl_ref_count *ref);

// Kernel tables (p=1 serendipity only; indexed by cdim-1)
GKYL_CU_D static const pkpm_kinetic_vel_vol_t ser_pkpm_kinetic_vel_vol_kernels[] = {
  pkpm_kinetic_vel_vol_1x1v_ser_p1, pkpm_kinetic_vel_vol_2x1v_ser_p1,
  pkpm_kinetic_vel_vol_3x1v_ser_p1
};

GKYL_CU_D static const pkpm_accel_surf_t ser_pkpm_accel_surf_kernels[] = {
  pkpm_accel_surfvpar_1x1v_ser_p1, pkpm_accel_surfvpar_2x1v_ser_p1, pkpm_accel_surfvpar_3x1v_ser_p1
};

GKYL_CU_D static void
choose_pkpm_kinetic_vel_kerns(
  struct dg_pkpm_kinetic_vel *pkpm, enum gkyl_basis_type b_type, int cdim, int poly_order
)
{
  assert(poly_order == 1);
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
    case GKYL_BASIS_MODAL_TENSOR: // the same basis at polynomial order 1
      pkpm->vol_kern = ser_pkpm_kinetic_vel_vol_kernels[cdim - 1];
      pkpm->accel_surf = ser_pkpm_accel_surf_kernels[cdim - 1];
      break;
    default:
      assert(false);
      break;
  }
}

// Volume term: the closure terms.
GKYL_CU_D static double
vol(
  const struct gkyl_dg_eqn *eqn, const double *xc, const double *dx, const int *idx,
  const double *qIn, double *GKYL_RESTRICT qRhsOut
)
{
  struct dg_pkpm_kinetic_vel *pkpm = container_of(eqn, struct dg_pkpm_kinetic_vel, eqn);
  long cidx = gkyl_range_idx(&pkpm->conf_range, idx);
  return pkpm->vol_kern(
    xc, dx, (const double *)gkyl_array_cfetch(pkpm->c_0, cidx),
    (const double *)gkyl_array_cfetch(pkpm->alpha_0, cidx),
    (const double *)gkyl_array_cfetch(pkpm->c_h, cidx),
    (const double *)gkyl_array_cfetch(pkpm->a_perp, cidx), qIn, qRhsOut
  );
}

// Surface term in v_par: each cell owns its lower face; the upper face is the lower
// face of the cell above.
GKYL_CU_D static double
surf(
  const struct gkyl_dg_eqn *eqn, int dir, const double *xcL, const double *xcC, const double *xcR,
  const double *dxL, const double *dxC, const double *dxR, const int *idxL, const int *idxC,
  const int *idxR, const double *qInL, const double *qInC, const double *qInR,
  double *GKYL_RESTRICT qRhsOut
)
{
  struct dg_pkpm_kinetic_vel *pkpm = container_of(eqn, struct dg_pkpm_kinetic_vel, eqn);
  if (dir == pkpm->cdim) {
    long pidxC = gkyl_range_idx(&pkpm->phase_range, idxC);
    long pidxR = gkyl_range_idx(&pkpm->phase_range, idxR);
    const double *flux_l = (const double *)gkyl_array_cfetch(pkpm->vel_flux_surf, pidxC);
    const double *flux_r = (const double *)gkyl_array_cfetch(pkpm->vel_flux_surf, pidxR);
    return pkpm->accel_surf(xcC, dxC, flux_l, flux_r, qRhsOut);
  }
  return 0.0;
}

// Lower velocity domain edge (no ghost cell below w = -w_max): the lowest cell's own slot
// holds the characteristic edge flux written by the closure, so the cell is scattered by
// the regular surface kernel with its slot and the slot of the cell above (the edge
// cell). The upper edge is a regular face against the velocity ghost cell, whose slot the
// closure fills too; it never comes through here.
GKYL_CU_D static double
boundary_surf(
  const struct gkyl_dg_eqn *eqn, int dir, const double *xcEdge, const double *xcSkin,
  const double *dxEdge, const double *dxSkin, const int *idxEdge, const int *idxSkin,
  const int edge, const double *qInEdge, const double *qInSkin, double *GKYL_RESTRICT qRhsOut
)
{
  struct dg_pkpm_kinetic_vel *pkpm = container_of(eqn, struct dg_pkpm_kinetic_vel, eqn);
  if (dir == pkpm->cdim && edge == -1) {
    long pidxSkin = gkyl_range_idx(&pkpm->phase_range, idxSkin);
    const double *flux_l = (const double *)gkyl_array_cfetch(pkpm->vel_flux_surf, pidxSkin);
    const double *flux_r = (const double *)gkyl_array_cfetch(
      pkpm->vel_flux_surf, gkyl_range_idx(&pkpm->phase_range, idxEdge)
    );
    return pkpm->accel_surf(xcSkin, dxSkin, flux_l, flux_r, qRhsOut);
  }
  return 0.0;
}
