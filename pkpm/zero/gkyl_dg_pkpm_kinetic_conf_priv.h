// Private header: not for direct use
#pragma once

#include <assert.h>

#include <gkyl_array.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_pkpm_kinetic_conf.h>
#include <gkyl_pkpm_kinetic_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

typedef double (*pkpm_kinetic_conf_vol_t)(
  const double *w, const double *dxv, const double *jacob_pos, const double *bvar, const double *u,
  const double *xi, const double *f, double *GKYL_RESTRICT out
);

typedef double (*pkpm_stream_surf_t)(
  const double *w, const double *dxv, const double *flux_l, const double *flux_r,
  double *GKYL_RESTRICT out
);

// for use in kernel tables
typedef struct {
  pkpm_stream_surf_t kernels[3];
} gkyl_dg_pkpm_stream_surf_kern_list;

struct dg_pkpm_kinetic_conf {
  struct gkyl_dg_eqn eqn; // Base object
  int cdim; // Configuration-space dimensions
  int pdim; // Phase-space dimensions
  int num_surf_basis; // Surface modal basis size per component and direction
  pkpm_kinetic_conf_vol_t vol_kern; // Volume kernel of the configuration pass
  pkpm_stream_surf_t stream_surf[3]; // Surface kernels scattering stored fluxes
  struct gkyl_range conf_range; // Configuration-space range
  struct gkyl_range phase_range; // Phase-space range
  const struct gkyl_vlasov_position_map *pos_map; // Position map (acquired; host pointer, 0 on the device copy)
  const struct gkyl_array *jacob_pos; // Per-direction Jacobians (borrowed from pos_map; device pointer on device)
  struct gkyl_dg_pkpm_kinetic_conf_auxfields auxfields; // Auxiliary fields.
};

// Volume kernels of the configuration pass (p=1 serendipity only; indexed by cdim-1)
GKYL_CU_D static const pkpm_kinetic_conf_vol_t ser_pkpm_kinetic_conf_vol_kernels[] = {
  pkpm_kinetic_conf_vol_1x1v_ser_p1,
  pkpm_kinetic_conf_vol_2x1v_ser_p1,
  pkpm_kinetic_conf_vol_3x1v_ser_p1,
};

GKYL_CU_D
static pkpm_kinetic_conf_vol_t
choose_pkpm_kinetic_conf_vol_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  assert(poly_order == 1);
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_pkpm_kinetic_conf_vol_kernels[cdim - 1];
      break;
    default:
      assert(false);
      break;
  }
  return 0;
}

// Stored-flux scatter kernels (p=1 serendipity only; indexed [cdim-1].kernels[dir])
GKYL_CU_D static const gkyl_dg_pkpm_stream_surf_kern_list ser_pkpm_stream_surf_kernels[] = {
  {pkpm_stream_surfx_1x1v_ser_p1, NULL, NULL}, // 0
  {pkpm_stream_surfx_2x1v_ser_p1, pkpm_stream_surfy_2x1v_ser_p1, NULL}, // 1
  {pkpm_stream_surfx_3x1v_ser_p1, pkpm_stream_surfy_3x1v_ser_p1, pkpm_stream_surfz_3x1v_ser_p1}, // 2
};

GKYL_CU_D
static pkpm_stream_surf_t
choose_pkpm_stream_surf_kern(enum gkyl_basis_type b_type, int cdim, int poly_order, int dir)
{
  assert(poly_order == 1);
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_pkpm_stream_surf_kernels[cdim - 1].kernels[dir];
      break;
    default:
      assert(false);
      break;
  }
  return 0;
}

// Volume term: streaming of the eight blocks.
GKYL_CU_D
static double
vol(
  const struct gkyl_dg_eqn *eqn, const double *xc, const double *dx, const int *idx,
  const double *qIn, double *GKYL_RESTRICT qRhsOut
)
{
  struct dg_pkpm_kinetic_conf *pkpm = container_of(eqn, struct dg_pkpm_kinetic_conf, eqn);
  long cidx = gkyl_range_idx(&pkpm->conf_range, idx);
  return pkpm->vol_kern(
    xc, dx, (const double *)gkyl_array_cfetch(pkpm->jacob_pos, cidx),
    (const double *)gkyl_array_cfetch(pkpm->auxfields.bvar, cidx),
    (const double *)gkyl_array_cfetch(pkpm->auxfields.u, cidx),
    (const double *)gkyl_array_cfetch(pkpm->auxfields.xi, cidx), qIn, qRhsOut
  );
}

// Surface term: each cell owns its lower face; the upper face is the lower face
// of the right neighbor (an interior cell or the upper ghost cell).
GKYL_CU_D
static double
surf(
  const struct gkyl_dg_eqn *eqn, int dir, const double *xcL, const double *xcC, const double *xcR,
  const double *dxL, const double *dxC, const double *dxR, const int *idxL, const int *idxC,
  const int *idxR, const double *qInL, const double *qInC, const double *qInR,
  double *GKYL_RESTRICT qRhsOut
)
{
  struct dg_pkpm_kinetic_conf *pkpm = container_of(eqn, struct dg_pkpm_kinetic_conf, eqn);
  if (dir < pkpm->cdim) {
    // Face-owned storage, [dir][side][block][mode]: this cell receives side 1 of its lower face and
    // side 0 of its upper face, which the upper neighbor owns.
    int off = 2 * dir * 8 * pkpm->num_surf_basis, side1 = 8 * pkpm->num_surf_basis;
    long pidxC = gkyl_range_idx(&pkpm->phase_range, idxC);
    long pidxR = gkyl_range_idx(&pkpm->phase_range, idxR);
    const double *flux_l = (const double *)gkyl_array_cfetch(pkpm->auxfields.conf_flux_surf, pidxC);
    const double *flux_r = (const double *)gkyl_array_cfetch(pkpm->auxfields.conf_flux_surf, pidxR);
    return pkpm->stream_surf[dir](xcC, dxC, &flux_l[off + side1], &flux_r[off], qRhsOut);
  }
  return 0.0;
}

// The configuration pass has ghost cells in every direction it updates, so no
// zero-flux boundary surfaces are ever requested of it.
GKYL_CU_D
static double
boundary_surf(
  const struct gkyl_dg_eqn *eqn, int dir, const double *xcEdge, const double *xcSkin,
  const double *dxEdge, const double *dxSkin, const int *idxEdge, const int *idxSkin,
  const int edge, const double *qInEdge, const double *qInSkin, double *GKYL_RESTRICT qRhsOut
)
{
  return 0.0;
}
