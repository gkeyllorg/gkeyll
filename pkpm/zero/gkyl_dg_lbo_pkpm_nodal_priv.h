#pragma once

// Private header, not for direct use in user code

#include <math.h>

#include <gkyl_lbo_pkpm_nodal_kernels.h>
#include <gkyl_pkpm_layout.h>

// Types for the kernels in the velocity alone at one configuration-space node
typedef void (*lbo_pkpm_nodal_vol_t)(
  const double *w, const double *dxv, double nu, double nu_u, double nu_vtsq, int stride,
  const double *f, double *GKYL_RESTRICT out
);

typedef void (*lbo_pkpm_nodal_surf_t)(
  const double *w, const double *dxv, double nu, double nu_u, double nu_vtsq, int stride,
  const double *fl, const double *fc, const double *fr, double *GKYL_RESTRICT out
);

typedef void (*lbo_pkpm_nodal_boundary_surf_t)(
  const double *w, const double *dxv, double nu, double nu_u, double nu_vtsq, int stride, int edge,
  const double *fedge, const double *fskin, double *GKYL_RESTRICT out
);

// for use in kernel tables, indexed by vdim
typedef struct {
  lbo_pkpm_nodal_vol_t kernels[3];
} gkyl_dg_lbo_pkpm_nodal_vol_kern_list;

typedef struct {
  lbo_pkpm_nodal_surf_t kernels[3];
} gkyl_dg_lbo_pkpm_nodal_surf_kern_list;

typedef struct {
  lbo_pkpm_nodal_boundary_surf_t kernels[3];
} gkyl_dg_lbo_pkpm_nodal_boundary_surf_kern_list;

// Volume kernel list (p=2 in the velocity): 1v, 2v, 3v
GKYL_CU_D static const gkyl_dg_lbo_pkpm_nodal_vol_kern_list lbo_pkpm_nodal_vol_kernels = {
  lbo_pkpm_nodal_vol_1v_p2, NULL, NULL
};

// Surface kernel list, vpar direction
GKYL_CU_D static const gkyl_dg_lbo_pkpm_nodal_surf_kern_list lbo_pkpm_nodal_surf_vpar_kernels = {
  lbo_pkpm_nodal_surfvpar_1v_p2, NULL, NULL
};

// Boundary surface kernel (zero-flux BCs) list, vpar direction
GKYL_CU_D static const gkyl_dg_lbo_pkpm_nodal_boundary_surf_kern_list
  lbo_pkpm_nodal_boundary_surf_vpar_kernels = {lbo_pkpm_nodal_boundary_surfvpar_1v_p2, NULL, NULL};

struct dg_lbo_pkpm_nodal {
  struct gkyl_dg_eqn eqn; // Base object
  int cdim; // Config-space dimensions
  int num_nodes; // nodes per configuration-space cell
  int num_vel_basis; // velocity modes at each node
  int stride; // offset between the blocks of a cell (number of phase-space basis functions)
  double vMaxSq; // square of the largest velocity of the grid
  lbo_pkpm_nodal_vol_t vol; // Volume term at a node
  lbo_pkpm_nodal_surf_t surf; // Surface term at a node
  lbo_pkpm_nodal_boundary_surf_t boundary_surf; // Boundary surface term at a node
  struct gkyl_range conf_range; // Configuration space range.
  struct gkyl_dg_lbo_pkpm_nodal_auxfields auxfields; // Auxiliary fields.
};

void gkyl_lbo_pkpm_nodal_free(const struct gkyl_ref_count *ref);

// The operator acts at a node only if its thermal speed squared is positive and inside the grid.
GKYL_CU_DH static inline bool
lbo_pkpm_nodal_check_vtsq(const struct dg_lbo_pkpm_nodal *lbo, double nu, double nu_vtsq)
{
  return nu_vtsq > 0.0 && nu_vtsq < nu * lbo->vMaxSq;
}

GKYL_CU_D static double
vol(
  const struct gkyl_dg_eqn *eqn, const double *xc, const double *dx, const int *idx,
  const double *qIn, double *GKYL_RESTRICT qRhsOut
)
{
  struct dg_lbo_pkpm_nodal *lbo = container_of(eqn, struct dg_lbo_pkpm_nodal, eqn);
  long cidx = gkyl_range_idx(&lbo->conf_range, idx);
  const double *nu_p = (const double *)gkyl_array_cfetch(lbo->auxfields.nu, cidx);
  const double *prim_moms_p = (const double *)gkyl_array_cfetch(lbo->auxfields.prim_moms, cidx);
  int cdim = lbo->cdim, num_nodes = lbo->num_nodes, num_vel_basis = lbo->num_vel_basis;

  // The rate of the cell is the largest of the rates of its nodes.
  double rdv = 2.0 / dx[cdim], rate_max = 0.0;
  for (int q = 0; q < num_nodes; ++q) {
    double nu = nu_p[q], nu_u = nu * prim_moms_p[q], nu_vtsq = nu * prim_moms_p[num_nodes + q];
    if (lbo_pkpm_nodal_check_vtsq(lbo, nu, nu_vtsq)) {
      lbo->vol(
        &xc[cdim], &dx[cdim], nu, nu_u, nu_vtsq, lbo->stride, &qIn[q * num_vel_basis],
        &qRhsOut[q * num_vel_basis]
      );
      double rate = 2.5 * rdv * fabs(nu_u - nu * xc[cdim]) + 3.0 * nu + 9.0 * rdv * rdv * nu_vtsq;
      rate_max = fmax(rate_max, rate);
    }
  }
  return rate_max;
}

GKYL_CU_D static double
surf(
  const struct gkyl_dg_eqn *eqn, int dir, const double *xcL, const double *xcC, const double *xcR,
  const double *dxL, const double *dxC, const double *dxR, const int *idxL, const int *idxC,
  const int *idxR, const double *qInL, const double *qInC, const double *qInR,
  double *GKYL_RESTRICT qRhsOut
)
{
  struct dg_lbo_pkpm_nodal *lbo = container_of(eqn, struct dg_lbo_pkpm_nodal, eqn);
  int cdim = lbo->cdim, num_nodes = lbo->num_nodes, num_vel_basis = lbo->num_vel_basis;
  if (dir < cdim) {
    return 0.;
  }
  long cidx = gkyl_range_idx(&lbo->conf_range, idxC);
  const double *nu_p = (const double *)gkyl_array_cfetch(lbo->auxfields.nu, cidx);
  const double *prim_moms_p = (const double *)gkyl_array_cfetch(lbo->auxfields.prim_moms, cidx);

  for (int q = 0; q < num_nodes; ++q) {
    double nu = nu_p[q], nu_u = nu * prim_moms_p[q], nu_vtsq = nu * prim_moms_p[num_nodes + q];
    if (lbo_pkpm_nodal_check_vtsq(lbo, nu, nu_vtsq)) {
      lbo->surf(
        &xcC[cdim], &dxC[cdim], nu, nu_u, nu_vtsq, lbo->stride, &qInL[q * num_vel_basis],
        &qInC[q * num_vel_basis], &qInR[q * num_vel_basis], &qRhsOut[q * num_vel_basis]
      );
    }
  }
  return 0.;
}

GKYL_CU_D static double
boundary_surf(
  const struct gkyl_dg_eqn *eqn, int dir, const double *xcEdge, const double *xcSkin,
  const double *dxEdge, const double *dxSkin, const int *idxEdge, const int *idxSkin,
  const int edge, const double *qInEdge, const double *qInSkin, double *GKYL_RESTRICT qRhsOut
)
{
  struct dg_lbo_pkpm_nodal *lbo = container_of(eqn, struct dg_lbo_pkpm_nodal, eqn);
  int cdim = lbo->cdim, num_nodes = lbo->num_nodes, num_vel_basis = lbo->num_vel_basis;
  if (dir < cdim) {
    return 0.;
  }
  long cidx = gkyl_range_idx(&lbo->conf_range, idxSkin);
  const double *nu_p = (const double *)gkyl_array_cfetch(lbo->auxfields.nu, cidx);
  const double *prim_moms_p = (const double *)gkyl_array_cfetch(lbo->auxfields.prim_moms, cidx);

  for (int q = 0; q < num_nodes; ++q) {
    double nu = nu_p[q], nu_u = nu * prim_moms_p[q], nu_vtsq = nu * prim_moms_p[num_nodes + q];
    if (lbo_pkpm_nodal_check_vtsq(lbo, nu, nu_vtsq)) {
      lbo->boundary_surf(
        &xcSkin[cdim], &dxSkin[cdim], nu, nu_u, nu_vtsq, lbo->stride, edge,
        &qInEdge[q * num_vel_basis], &qInSkin[q * num_vel_basis], &qRhsOut[q * num_vel_basis]
      );
    }
  }
  return 0.;
}
