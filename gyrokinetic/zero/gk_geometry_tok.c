#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_comm.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_util.h>
#include <gkyl_gk_geometry_priv.h>
#include <gkyl_gk_geometry_tok.h>
#include <gkyl_math.h>
#include <gkyl_nodal_ops.h>

#include <gkyl_tok_geo.h>
#include <gkyl_rz_calc_derived_geo.h>
#include <gkyl_calc_metric.h>
#include <gkyl_calc_bmag.h>
#include <assert.h>
#include <float.h>
#include <stdint.h>

static bool
tok_xpt_trial_active(const struct gkyl_tok_geo_grid_inp *ginp)
{
  return ginp->relaxed_xpt_seam_optimizer_trial &&
    ginp->relaxed_xpt_seam_trial_status;
}

static void
tok_xpt_trial_reject(struct gkyl_tok_geo_xpt_seam_trial_status *status,
  enum gkyl_tok_geo_xpt_seam_trial_failure reason)
{
  if (status->first_failure_reason == GKYL_XPT_SEAM_TRIAL_OK)
    status->first_failure_reason = reason;
  switch (reason) {
    case GKYL_XPT_SEAM_TRIAL_NONFINITE_MAP:
      status->finite_map_valid = false;
      break;
    case GKYL_XPT_SEAM_TRIAL_CELL_JACOBIAN:
      status->cell_jacobian_valid = false;
      break;
    case GKYL_XPT_SEAM_TRIAL_METRIC_JACOBIAN:
      status->jacobian_valid = false;
      break;
    case GKYL_XPT_SEAM_TRIAL_REMOTE_FAILURE:
      break;
    default:
      break;
  }
}

static bool
tok_xpt_trial_mapping_valid(
  const struct gkyl_tok_geo_xpt_seam_trial_status *status)
{
  return status->first_failure_reason == GKYL_XPT_SEAM_TRIAL_OK &&
    status->contour_valid && status->branch_valid &&
    status->trace_ordering_valid && status->finite_map_valid &&
    status->cell_jacobian_valid && status->jacobian_valid;
}

static void
tok_xpt_trial_sync_status(struct gkyl_comm *comm,
  struct gkyl_tok_geo_xpt_seam_trial_status *status)
{
  int64_t local_invalid[7] = {
    !status->contour_valid, !status->branch_valid,
    !status->trace_ordering_valid, !status->finite_map_valid,
    !status->cell_jacobian_valid, !status->jacobian_valid,
    status->first_failure_reason,
  };
  int64_t global_invalid[7] = { 0 };
  gkyl_comm_allreduce_host(comm, GKYL_INT_64, GKYL_MAX, 7,
    local_invalid, global_invalid);
  status->contour_valid = !global_invalid[0];
  status->branch_valid = !global_invalid[1];
  status->trace_ordering_valid = !global_invalid[2];
  status->finite_map_valid = !global_invalid[3];
  status->cell_jacobian_valid = !global_invalid[4];
  status->jacobian_valid = !global_invalid[5];
  if (status->first_failure_reason == GKYL_XPT_SEAM_TRIAL_OK &&
      global_invalid[6] != GKYL_XPT_SEAM_TRIAL_OK)
    status->first_failure_reason = global_invalid[6];

  int64_t local_sign = status->jacobian_sign;
  int64_t min_sign = 0, max_sign = 0;
  gkyl_comm_allreduce_host(comm, GKYL_INT_64, GKYL_MIN, 1,
    &local_sign, &min_sign);
  gkyl_comm_allreduce_host(comm, GKYL_INT_64, GKYL_MAX, 1,
    &local_sign, &max_sign);
  if (tok_xpt_trial_mapping_valid(status)) {
    if (min_sign == 0 || min_sign != max_sign)
      tok_xpt_trial_reject(status, GKYL_XPT_SEAM_TRIAL_REMOTE_FAILURE);
    else
      status->jacobian_sign = min_sign;
  }

  double global_margin = 0.0, global_displacement = 0.0;
  gkyl_comm_allreduce_host(comm, GKYL_DOUBLE, GKYL_MIN, 1,
    &status->min_cell_jacobian_margin, &global_margin);
  gkyl_comm_allreduce_host(comm, GKYL_DOUBLE, GKYL_MAX, 1,
    &status->max_realized_displacement, &global_displacement);
  status->min_cell_jacobian_margin = global_margin;
  status->max_realized_displacement = global_displacement;
}

// Sign consistency is checked per representation (corner/interior/surface),
// never across them: the real, always-used metric pipeline in calc_metric.c
// (gkyl_calc_metric_advance_rz for corners, _advance_rz_interior, and
// _advance_rz_surface) each declare and check their own independent
// signed_jacobian_guard_state and never compare against each other -- a
// corner-grid bilinear/affine Jacobian and an interior-grid finite-difference
// Jacobian are different approximations of the same underlying map and are
// not expected to agree in sign at a highly compressed cell (e.g. next to an
// X-point-adjacent seam, where GKYL_PMAP_XPT_COMPRESSION makes cells
// extremely anisotropic and the affine approximation is poor). Confirmed
// directly against calc_metric.c: gkyl_calc_metric_advance_rz's own signed_J
// guard (jac_guard local to that function), _advance_rz_interior's (a
// separate local jac_guard), and _advance_rz_surface's (a third, separate
// local jac_guard) are three independent instances, never cross-checked.
// sign_state is therefore caller-owned (a field in *status for interior,
// which is what later gets compared against the real metric calculation's
// own interior sign -- see gk_geometry_tok_init -- or a call-local variable
// for corner/surface, which have no such downstream comparison).
static bool
tok_xpt_trial_set_jacobian_sign(
  struct gkyl_tok_geo_xpt_seam_trial_status *status, int *sign_state,
  double jacobian, double scale)
{
  double relative = fabs(jacobian)/fmax(scale, DBL_MIN);
  if (!isfinite(jacobian) || !isfinite(relative) ||
      relative <= 256.0*DBL_EPSILON) {
    tok_xpt_trial_reject(status, GKYL_XPT_SEAM_TRIAL_CELL_JACOBIAN);
    return false;
  }
  int sign = jacobian < 0.0 ? -1 : 1;
  if (*sign_state == 0)
    *sign_state = sign;
  else if (sign != *sign_state) {
    tok_xpt_trial_reject(status, GKYL_XPT_SEAM_TRIAL_CELL_JACOBIAN);
    return false;
  }
  status->min_cell_jacobian_margin =
    fmin(status->min_cell_jacobian_margin, relative);
  return true;
}

// Corner_map's bilinear/affine Jacobian (see tok_xpt_trial_preflight_corner_map)
// is a coarse, 4-Gauss-point approximation with no analog anywhere in the
// real metric pipeline (calc_metric.c's own corner-grid Jacobian is a
// finite-difference formula, not this affine reconstruction), so unlike
// interior/surface it has no real-pipeline sign convention to preserve.
// Checks only finiteness and non-degenerate magnitude, no sign tracking --
// see the comment on tok_xpt_trial_set_jacobian_sign for why cross-cell
// (or, as confirmed here, even cross-reference-corner) sign agreement isn't
// a real invariant. Root-caused empirically (resolution-convergence check,
// same method as the "radial-line kink"/Phase 3 investigations referenced
// in the ixf-vs-main-comp memory): doubling theta resolution on both a
// genuine two-block pair (STEP's DN_SOL_OUT_LO<->MID) and a self-periodic
// one (ASDEX's CORE) made the sign disagreement disappear entirely at the
// exact same physical location, the classic signature of a coarse-mesh
// approximation artifact (shrinks under refinement) rather than a real
// fold (which would stay anchored or grow). Confirmed at the finer
// resolution the optimizer went on to find and apply a genuine 8.97%
// improvement on the STEP pair that the coarse-resolution false rejection
// had been hiding entirely.
static bool
tok_xpt_trial_check_cell_jacobian_magnitude(
  struct gkyl_tok_geo_xpt_seam_trial_status *status, double jacobian,
  double scale)
{
  double relative = fabs(jacobian)/fmax(scale, DBL_MIN);
  if (!isfinite(jacobian) || !isfinite(relative) ||
      relative <= 256.0*DBL_EPSILON) {
    tok_xpt_trial_reject(status, GKYL_XPT_SEAM_TRIAL_CELL_JACOBIAN);
    return false;
  }
  status->min_cell_jacobian_margin =
    fmin(status->min_cell_jacobian_margin, relative);
  return true;
}

// Check the p1 DG R-Z map before metric construction.  Its four Gauss values
// determine the bilinear cell map; checking all four reference-cell corners
// catches a nonfinite or degenerate (near-zero) map before metric
// construction. Does not require the four corners' Jacobian signs to agree
// with each other -- see tok_xpt_trial_check_cell_jacobian_magnitude.
static bool
tok_xpt_trial_preflight_corner_map(struct gk_geometry *up,
  struct gkyl_tok_geo_xpt_seam_trial_status *status)
{
  if (up->basis.poly_order != 1 || up->grid.ndim != 3) {
    tok_xpt_trial_reject(status, GKYL_XPT_SEAM_TRIAL_CELL_JACOBIAN);
    return false;
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->nrange_corn);
  while (gkyl_range_iter_next(&iter)) {
    const double *p = gkyl_array_cfetch(up->geo_corn.mc2p_nodal,
      gkyl_range_idx(&up->nrange_corn, iter.idx));
    if (!(p[0] > 0.0) || !isfinite(p[0]) || !isfinite(p[1]) ||
        !isfinite(p[2])) {
      tok_xpt_trial_reject(status, GKYL_XPT_SEAM_TRIAL_NONFINITE_MAP);
      return false;
    }
  }

  double gauss = 1.0/sqrt(3.0);
  int idx[4][3];
  static const double reference_corner[4][2] = {
    { -1.0, -1.0 }, { 1.0, -1.0 },
    { -1.0, 1.0 }, { 1.0, 1.0 },
  };
  int ia = up->nrange_int.lower[1];
  for (int ip=up->nrange_int.lower[0]; ip<=up->nrange_int.upper[0]; ip+=2) {
    for (int it=up->nrange_int.lower[2]; it<=up->nrange_int.upper[2]; it+=2) {
      idx[0][0] = ip;   idx[0][1] = ia; idx[0][2] = it;
      idx[1][0] = ip+1; idx[1][1] = ia; idx[1][2] = it;
      idx[2][0] = ip;   idx[2][1] = ia; idx[2][2] = it+1;
      idx[3][0] = ip+1; idx[3][1] = ia; idx[3][2] = it+1;
      const double *point[4];
      for (int n=0; n<4; ++n) {
        point[n] = gkyl_array_cfetch(up->geo_int.mc2p_nodal,
          gkyl_range_idx(&up->nrange_int, idx[n]));
        if (!(point[n][0] > 0.0) || !isfinite(point[n][0]) ||
            !isfinite(point[n][1]) || !isfinite(point[n][2])) {
          tok_xpt_trial_reject(status, GKYL_XPT_SEAM_TRIAL_NONFINITE_MAP);
          return false;
        }
      }
      double coeff[2][4];
      for (int c=0; c<2; ++c) {
        double xmm = point[0][c], xpm = point[1][c];
        double xmp = point[2][c], xpp = point[3][c];
        coeff[c][0] = 0.25*(xmm+xpm+xmp+xpp);
        coeff[c][1] = (-xmm+xpm-xmp+xpp)/(4.0*gauss);
        coeff[c][2] = (-xmm-xpm+xmp+xpp)/(4.0*gauss);
        coeff[c][3] = (xmm-xpm-xmp+xpp)/(4.0*gauss*gauss);
      }
      for (int n=0; n<4; ++n) {
        double xi = reference_corner[n][0];
        double eta = reference_corner[n][1];
        double dR_dxi = coeff[0][1]+coeff[0][3]*eta;
        double dZ_dxi = coeff[1][1]+coeff[1][3]*eta;
        double dR_deta = coeff[0][2]+coeff[0][3]*xi;
        double dZ_deta = coeff[1][2]+coeff[1][3]*xi;
        double jacobian = dR_dxi*dZ_deta-dR_deta*dZ_dxi;
        double scale = hypot(dR_dxi, dZ_dxi)*hypot(dR_deta, dZ_deta);
        if (!tok_xpt_trial_check_cell_jacobian_magnitude(status, jacobian,
            scale))
          return false;
      }
    }
  }
  return true;
}

// Preflight the exact R-Z Jacobian algebra used by the interior metric loop.
// This catches a nonfinite/degenerate finite-difference map before any metric
// divisions or square roots are evaluated.
static bool
tok_xpt_trial_preflight_interior_map(struct gk_geometry *up,
  struct gkyl_tok_geo_xpt_seam_trial_status *status)
{
  int idx[3];
  for (int ia=up->nrange_int.lower[1]; ia<=up->nrange_int.upper[1]; ++ia) {
    for (int ip=up->nrange_int.lower[0]; ip<=up->nrange_int.upper[0]; ++ip) {
      for (int it=up->nrange_int.lower[2]; it<=up->nrange_int.upper[2]; ++it) {
        idx[0] = ip; idx[1] = ia; idx[2] = it;
        long loc = gkyl_range_idx(&up->nrange_int, idx);
        const double *map = gkyl_array_cfetch(up->geo_int.mc2p_nodal_fd, loc);
        const double *dtheta = gkyl_array_cfetch(up->geo_int.ddtheta_nodal, loc);
        const double *dpsi = gkyl_array_cfetch(up->geo_int.ddpsi_nodal, loc);
        const double *bmag = gkyl_array_cfetch(up->geo_int.bmag_nodal, loc);
        bool finite = isfinite(dtheta[0]) && isfinite(dtheta[1]) &&
          isfinite(dtheta[2]) && isfinite(dpsi[0]) && dpsi[0] != 0.0 &&
          isfinite(bmag[0]) && bmag[0] > 0.0;
        for (int n=0; n<39; ++n)
          finite = finite && isfinite(map[n]);
        double R = map[0];
        finite = finite && R > 0.0;
        if (!finite) {
          tok_xpt_trial_reject(status, GKYL_XPT_SEAM_TRIAL_NONFINITE_MAP);
          return false;
        }
        double dR_dpsi = -(map[3]-map[6])/(2.0*up->dzc[0]);
        double dZ_dpsi = -(map[4]-map[7])/(2.0*up->dzc[0]);
        double jacobian = R*(dR_dpsi*dtheta[1]-dtheta[0]*dZ_dpsi);
        double scale = R*hypot(dR_dpsi, dZ_dpsi)
          *hypot(dtheta[0], dtheta[1]);
        if (!tok_xpt_trial_set_jacobian_sign(status, &status->jacobian_sign,
            jacobian, scale))
          return false;
        double radicand = (jacobian*jacobian*bmag[0]*bmag[0]
            /(dpsi[0]*dpsi[0])-dtheta[0]*dtheta[0]
            -dtheta[1]*dtheta[1])/(R*R);
        if (!isfinite(radicand) || radicand < 0.0) {
          tok_xpt_trial_reject(status,
            GKYL_XPT_SEAM_TRIAL_METRIC_JACOBIAN);
          return false;
        }
      }
    }
  }
  return true;
}

// sign_state is call-local: the real gkyl_calc_metric_advance_rz_surface is
// invoked once per dir with its own fresh, local signed_jacobian_guard_state
// (never shared across dir, let alone with corner/interior), so this checks
// self-consistency within one surface direction only -- see the comment on
// tok_xpt_trial_set_jacobian_sign for the corner/interior rationale, which
// applies identically here.
static bool
tok_xpt_trial_check_metric_jacobian(
  struct gkyl_tok_geo_xpt_seam_trial_status *status, int *sign_state,
  double jacobian, double scale)
{
  double relative = fabs(jacobian)/fmax(scale, DBL_MIN);
  int sign = jacobian < 0.0 ? -1 : 1;
  if (!isfinite(jacobian) || !isfinite(relative) ||
      relative <= 256.0*DBL_EPSILON ||
      (*sign_state != 0 && sign != *sign_state)) {
    tok_xpt_trial_reject(status, GKYL_XPT_SEAM_TRIAL_METRIC_JACOBIAN);
    return false;
  }
  if (*sign_state == 0)
    *sign_state = sign;
  return true;
}

// Check the surface prerequisites used by calc_metric before it performs
// one-sided reconstruction, divisions, or square roots.
static bool
tok_xpt_trial_preflight_surface_maps(struct gk_geometry *up,
  struct gkyl_tok_geo_xpt_seam_trial_status *status)
{
  int idx[3];
  for (int dir=0; dir<up->grid.ndim; ++dir) {
    const struct gkyl_range *range = &up->nrange_surf[dir];
    // Fresh per dir: matches gkyl_calc_metric_advance_rz_surface, which is
    // called once per dir with its own fresh local guard state.
    int surface_jacobian_sign = 0;
    for (int ia=range->lower[1]; ia<=range->upper[1]; ++ia) {
      for (int ip=range->lower[0]; ip<=range->upper[0]; ++ip) {
        for (int it=range->lower[2]; it<=range->upper[2]; ++it) {
          idx[0] = ip; idx[1] = ia; idx[2] = it;
          long loc = gkyl_range_idx(range, idx);
          const double *map =
            gkyl_array_cfetch(up->geo_surf[dir].mc2p_nodal_fd, loc);
          const double *dtheta =
            gkyl_array_cfetch(up->geo_surf[dir].ddtheta_nodal, loc);
          const double *dpsi =
            gkyl_array_cfetch(up->geo_surf[dir].ddpsi_nodal, loc);
          const double *bmag =
            gkyl_array_cfetch(up->geo_surf[dir].bmag_nodal, loc);
          bool finite = isfinite(dtheta[0]) && isfinite(dtheta[1]) &&
            isfinite(dtheta[2]) && isfinite(dpsi[0]) && dpsi[0] != 0.0 &&
            isfinite(bmag[0]) && bmag[0] > 0.0;
          for (int n=0; n<39; ++n)
            finite = finite && isfinite(map[n]);
          double R = map[0];
          finite = finite && R > 0.0;
          if (!finite) {
            tok_xpt_trial_reject(status, GKYL_XPT_SEAM_TRIAL_NONFINITE_MAP);
            return false;
          }

          double dR_dpsi, dZ_dpsi;
          if (dir == 0 && ip == range->lower[0] &&
              up->local.lower[0] == up->global.lower[0]) {
            dR_dpsi = (-3.0*map[0]+4.0*map[6]-map[12])
              /(2.0*up->dzc[0]);
            dZ_dpsi = (-3.0*map[1]+4.0*map[7]-map[13])
              /(2.0*up->dzc[0]);
          }
          else if (dir == 0 && ip == range->upper[0] &&
              up->local.upper[0] == up->global.upper[0]) {
            dR_dpsi = (3.0*map[0]-4.0*map[3]+map[9])
              /(2.0*up->dzc[0]);
            dZ_dpsi = (3.0*map[1]-4.0*map[4]+map[10])
              /(2.0*up->dzc[0]);
          }
          else {
            dR_dpsi = -(map[3]-map[6])/(2.0*up->dzc[0]);
            dZ_dpsi = -(map[4]-map[7])/(2.0*up->dzc[0]);
          }

          bool psi_boundary = dir == 0 &&
            ((ip == range->lower[0] &&
                up->local.lower[0] == up->global.lower[0]) ||
              (ip == range->upper[0] &&
                up->local.upper[0] == up->global.upper[0]));
          if (psi_boundary) {
            double tangent_sq = dtheta[0]*dtheta[0]+dtheta[1]*dtheta[1];
            double g33_exact = tangent_sq+R*R*dtheta[2]*dtheta[2];
            if (!(tangent_sq > 0.0) || !(g33_exact > 0.0) ||
                !isfinite(tangent_sq) || !isfinite(g33_exact)) {
              tok_xpt_trial_reject(status,
                GKYL_XPT_SEAM_TRIAL_METRIC_JACOBIAN);
              return false;
            }
            double cross_exact = status->jacobian_sign
              *fabs(dpsi[0])*sqrt(g33_exact)/(R*bmag[0]);
            double tangent_projection =
              (dR_dpsi*dtheta[0]+dZ_dpsi*dtheta[1])/tangent_sq;
            dR_dpsi = tangent_projection*dtheta[0]
              +cross_exact*dtheta[1]/tangent_sq;
            dZ_dpsi = tangent_projection*dtheta[1]
              -cross_exact*dtheta[0]/tangent_sq;
          }

          double jacobian =
            R*(dR_dpsi*dtheta[1]-dtheta[0]*dZ_dpsi);
          double scale = R*hypot(dR_dpsi, dZ_dpsi)
            *hypot(dtheta[0], dtheta[1]);
          if (!tok_xpt_trial_check_metric_jacobian(status,
              &surface_jacobian_sign, jacobian, scale))
            return false;
          if (!psi_boundary) {
            double radicand = (jacobian*jacobian*bmag[0]*bmag[0]
                /(dpsi[0]*dpsi[0])-dtheta[0]*dtheta[0]
                -dtheta[1]*dtheta[1])/(R*R);
            if (!isfinite(radicand) || radicand < 0.0) {
              tok_xpt_trial_reject(status,
                GKYL_XPT_SEAM_TRIAL_METRIC_JACOBIAN);
              return false;
            }
          }
        }
      }
    }
  }
  return true;
}

static void
tok_geometry_build_maps(struct gk_geometry *up, struct gkyl_tok_geo *geo,
  struct gkyl_tok_geo_grid_inp *ginp,
  struct gkyl_position_map *position_map)
{
  gkyl_tok_geo_calc(up, &up->nrange_corn, geo, ginp, position_map);
  gkyl_tok_geo_calc_interior(up, &up->nrange_int, up->dzc, geo, ginp,
    position_map);
  for (int dir=0; dir<up->grid.ndim; ++dir)
    gkyl_tok_geo_calc_surface(up, dir, &up->nrange_surf[dir], up->dzc, geo,
      ginp, position_map);
}

static bool
tok_xpt_trial_array_finite(const struct gkyl_array *array)
{
  for (long i=0; i<array->size; ++i) {
    const double *values = gkyl_array_cfetch(array, i);
    for (int c=0; c<array->ncomp; ++c)
      if (!isfinite(values[c]))
        return false;
  }
  return true;
}

static bool
tok_xpt_trial_interior_metric_finite(const struct gk_geometry *up)
{
  const struct gkyl_array *arrays[] = {
    up->geo_int.jacobgeo_nodal, up->geo_int.g_ij_nodal,
    up->geo_int.dxdz_nodal, up->geo_int.dzdx_nodal,
    up->geo_int.dualmag_nodal, up->geo_int.b_i_nodal,
    up->geo_int.bcart_nodal, up->geo_int.B3_nodal,
  };
  for (int i=0; i<(int) (sizeof arrays/sizeof arrays[0]); ++i)
    if (!tok_xpt_trial_array_finite(arrays[i]))
      return false;
  return true;
}

static bool
tok_xpt_trial_surface_metric_finite(const struct gk_geometry *up)
{
  for (int dir=0; dir<up->grid.ndim; ++dir) {
    const struct gk_geom_surf *surf = &up->geo_surf[dir];
    const struct gkyl_array *arrays[] = {
      surf->jacobgeo_nodal, surf->jacobgeo_signed_nodal,
      surf->g_ij_nodal, surf->dxdz_nodal, surf->dzdx_nodal,
      surf->dualmag_nodal, surf->b_i_nodal, surf->bcart_nodal,
      surf->B3_nodal, surf->jacobtot_inv_nodal,
      surf->bimpactangle_nodal,
    };
    for (int i=0; i<(int) (sizeof arrays/sizeof arrays[0]); ++i)
      if (!tok_xpt_trial_array_finite(arrays[i]))
        return false;
  }
  return true;
}

struct gk_geometry*
gk_geometry_tok_init(struct gkyl_gk_geometry_inp *geometry_inp)
{

  struct gk_geometry *up = gkyl_calloc(1,sizeof(struct gk_geometry));
  up->geometry_id = geometry_inp->geometry_id;
  up->basis = geometry_inp->geo_basis;

  // Construction forms bmag_inv = 1/|B| with gkyl_dg_inv_op_range, whose
  // kernel table is sparse. Ask before building rather than discovering it as a
  // NULL function pointer several thousand nodes in: at poly_order 2 the DG
  // inverse kernel does not exist in 3x, so a p2 tokamak geometry cannot be
  // completed no matter what the gridding does. Refuse here, while nothing has
  // been allocated and no geometry file has been written.
  //
  // This asks whether the kernel exists rather than testing the order, so the
  // refusal lifts by itself once the kernel is generated.
  if (!gkyl_dg_inv_op_supported(&up->basis)) {
    fprintf(stderr,
      "GKYL_GEOMETRY_UNSUPPORTED_POLY_ORDER poly_order=%d ndim=%d "
      "reason=no_dg_inverse_kernel context=bmag_inv\n",
      up->basis.poly_order, up->basis.ndim);
    gkyl_exit("gk_geometry_tok: no DG inverse kernel for this basis, so "
      "bmag_inv cannot be formed; polynomial order 1 is the implemented "
      "contract");
  }
  up->local = geometry_inp->geo_local;
  up->local_ext = geometry_inp->geo_local_ext;
  up->global = geometry_inp->geo_global;
  up->global_ext = geometry_inp->geo_global_ext;
  up->grid = geometry_inp->geo_grid;
  up->has_LCFS = geometry_inp->has_LCFS;
  if (up->has_LCFS) {
    up->x_LCFS = geometry_inp->x_LCFS;
    // Check that the split happens within the domain.
    assert((up->grid.lower[0] <= up->x_LCFS) && (up->x_LCFS <= up->grid.upper[0]));
    // Check that the split happens at a cell boundary;
    double needint = (up->x_LCFS - up->grid.lower[0])/up->grid.dx[0];
    double rem_floor = fabs(needint-floor(needint));
    double rem_ceil = fabs(needint-ceil(needint));
    if (rem_floor < 1.0e-12) {
      up->idx_LCFS_lo = (int) floor(needint);
    }
    else if (rem_ceil < 1.0e-12) {
      up->idx_LCFS_lo = (int) ceil(needint);
    }
    else {
      fprintf(stderr, "x_LCFS = %.9e must be at a cell boundary.\n", up->x_LCFS);
      assert(false);
    }
  }

  gk_geometry_set_nodal_ranges(up) ;

  // Initialize surface basis
  gkyl_cart_modal_serendip(&up->surf_basis, up->grid.ndim-1, up->basis.poly_order);
  up->num_surf_basis = up->surf_basis.num_basis;

  // Initialize tokamak geometry object from EFIT
  const struct gkyl_efit_inp inp = geometry_inp->efit_info;
  struct gkyl_tok_geo_grid_inp ginp = geometry_inp->tok_grid_info;
  ginp.cgrid = up->grid;
  ginp.cbasis = up->basis;
  struct gkyl_tok_geo *geo = gkyl_tok_geo_new(&inp, &ginp);
  up->geqdsk_sign_convention = geo->efit->sibry > geo->efit->simag ? 0 : 1;
  up->half_domain = ginp.half_domain ? 1 : 0;

  // Allocate nodal and modal arrays for corner, interior, and surface geo
  gk_geometry_corn_alloc_nodal(up);
  gk_geometry_corn_alloc_expansions(up);
  gk_geometry_int_alloc_nodal(up);
  gk_geometry_int_alloc_expansions(up);
  for (int dir=0; dir<up->grid.ndim; ++dir) {
    gk_geometry_surf_alloc_nodal(up, dir);
    gk_geometry_surf_alloc_expansions(up, dir);
  }

  gkyl_position_map_optimize(geometry_inp->position_map, up->grid, up->global);

  bool optimizer_trial = tok_xpt_trial_active(&ginp);
  struct gkyl_tok_geo_xpt_seam_trial_status *trial_status =
    ginp.relaxed_xpt_seam_trial_status;
  bool straight_placeholder = false, geometry_built = false;
  for (int pass=0; pass<2 && !geometry_built; ++pass) {
    struct gkyl_tok_geo_grid_inp pass_ginp = ginp;
    bool candidate_pass = optimizer_trial && !straight_placeholder;
    if (optimizer_trial && straight_placeholder) {
      pass_ginp.relaxed_xpt_seam = false;
      pass_ginp.relaxed_xpt_seam_delta_s_coeff = 0.0;
      pass_ginp.relaxed_xpt_seam_sweep = false;
      pass_ginp.relaxed_xpt_seam_optimize = false;
      pass_ginp.relaxed_xpt_seam_optimizer_trial = false;
      pass_ginp.relaxed_xpt_seam_trial_status = 0;
    }

    // Calculate bmag and mapc2p in cylindrical coordinates at corner,
    // interior, and surface nodes.  Every retry overwrites all three array
    // families so a rejected candidate can never leave a mixed map behind.
    tok_geometry_build_maps(up, geo, &pass_ginp,
      geometry_inp->position_map);
    if (candidate_pass && tok_xpt_trial_mapping_valid(trial_status)) {
      tok_xpt_trial_preflight_corner_map(up, trial_status);
      if (tok_xpt_trial_mapping_valid(trial_status))
        tok_xpt_trial_preflight_interior_map(up, trial_status);
      if (tok_xpt_trial_mapping_valid(trial_status))
        tok_xpt_trial_preflight_surface_maps(up, trial_status);
    }
    if (candidate_pass) {
      tok_xpt_trial_sync_status(geometry_inp->comm, trial_status);
      if (!tok_xpt_trial_mapping_valid(trial_status)) {
        straight_placeholder = true;
        continue;
      }
    }

    struct gkyl_calc_metric *mcalc = gkyl_calc_metric_new(&up->basis,
      &up->grid, &up->global, &up->global_ext, &up->local, &up->local_ext,
      true, false);
    if (candidate_pass)
      gkyl_calc_metric_set_signed_jacobian_guard_nonfatal(mcalc, true);
    gkyl_calc_metric_advance_rz_interior(mcalc, up);
    if (candidate_pass) {
      int metric_sign = gkyl_calc_metric_signed_jacobian_sign(mcalc);
      if (!gkyl_calc_metric_signed_jacobian_valid(mcalc) ||
          metric_sign == 0 ||
          (trial_status->jacobian_sign != 0 &&
            metric_sign != trial_status->jacobian_sign) ||
          !tok_xpt_trial_interior_metric_finite(up))
        tok_xpt_trial_reject(trial_status,
          GKYL_XPT_SEAM_TRIAL_METRIC_JACOBIAN);
      else
        trial_status->jacobian_sign = metric_sign;
      tok_xpt_trial_sync_status(geometry_inp->comm, trial_status);
      if (!tok_xpt_trial_mapping_valid(trial_status)) {
        gkyl_calc_metric_release(mcalc);
        straight_placeholder = true;
        continue;
      }
    }

    gkyl_array_copy(up->geo_int.jacobgeo_ghost, up->geo_int.jacobgeo);
    gkyl_calc_metric_advance_rz_neut_interior(mcalc, up);
    gkyl_rz_calc_derived_geo *jcalculator =
      gkyl_rz_calc_derived_geo_new(&up->basis, &up->grid, 1, false);
    gkyl_rz_calc_derived_geo_advance(jcalculator, &up->local,
      up->geo_int.g_ij, up->geo_int.bmag, up->geo_int.jacobgeo,
      up->geo_int.jacobgeo_inv, up->geo_int.gij, up->geo_int.b_i,
      up->geo_int.cmag, up->geo_int.jacobtot,
      up->geo_int.jacobtot_inv, up->geo_int.gxxj, up->geo_int.gxyj,
      up->geo_int.gyyj, up->geo_int.gxzj, up->geo_int.eps2);
    gkyl_rz_calc_derived_geo_release(jcalculator);
    for (int dir=0; dir<up->grid.ndim; ++dir)
      gkyl_calc_metric_advance_rz_surface(mcalc, dir, up);
    if (candidate_pass) {
      if (!gkyl_calc_metric_signed_jacobian_valid(mcalc) ||
          !tok_xpt_trial_surface_metric_finite(up))
        tok_xpt_trial_reject(trial_status,
          GKYL_XPT_SEAM_TRIAL_METRIC_JACOBIAN);
      tok_xpt_trial_sync_status(geometry_inp->comm, trial_status);
    }
    gkyl_calc_metric_release(mcalc);
    if (candidate_pass && !tok_xpt_trial_mapping_valid(trial_status)) {
      straight_placeholder = true;
      continue;
    }
    geometry_built = true;
  }
  assert(geometry_built);

  // Calculate surface expansions.
  for (int dir = 0; dir <up->grid.ndim; dir++)
    gk_geometry_surf_calc_expansions(up, dir, up->nrange_surf[dir]);

  // Store metadata for I/O.
  char geqdsk_file_name[128];
  get_filename_from_path(geo->efit->filepath, geqdsk_file_name, sizeof(geqdsk_file_name));
  struct gkyl_msgpack_map_elem io_meta_basic[] = {
    { .key = "geometry_type", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = up->geometry_id },
    { .key = "geqdsk_sign_convention", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = up->geqdsk_sign_convention },
    { .key = "geqdsk_file", .elem_type = GKYL_MP_STRING, .cval = geqdsk_file_name},
    { .key = "half_domain", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = up->half_domain },
  };
  up->io_meta_basic_len = sizeof(io_meta_basic)/sizeof(io_meta_basic[0]);
  up->io_meta_basic = gkyl_msgpack_map_elem_clone(up->io_meta_basic_len, io_meta_basic);

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->ref_count = gkyl_ref_count_init(gkyl_gk_geometry_free);
  up->on_dev = up; // CPU eqn obj points to itself.

  gkyl_tok_geo_release(geo);

  return up;
}

struct gk_geometry*
gkyl_gk_geometry_tok_new(struct gkyl_gk_geometry_inp *geometry_inp)
{
  struct gk_geometry* gk_geom_3d;
  struct gk_geometry* gk_geom;

  if (geometry_inp->position_map->id == GKYL_PMAP_XPT_COMPRESSION) {
    double zcenter, zcut, len;
    switch(geometry_inp->tok_grid_info.ftype)
    {
      case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID:
      case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID:
        len = geometry_inp->tok_grid_info.half_domain ? 2.0*(geometry_inp->geo_grid.upper[2] - geometry_inp->geo_grid.lower[2])
                                                      : geometry_inp->geo_grid.upper[2] - geometry_inp->geo_grid.lower[2];
        zcut = len/2.0;
        zcenter = 0.0;
        break;
      case GKYL_GEOMETRY_TOKAMAK_CORE_R:
      case GKYL_GEOMETRY_TOKAMAK_CORE:
      case GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID:
        len = geometry_inp->tok_grid_info.half_domain ? 2.0*(geometry_inp->geo_grid.upper[2] - geometry_inp->geo_grid.lower[2])
                                                      : geometry_inp->geo_grid.upper[2] - geometry_inp->geo_grid.lower[2];
        zcenter = geometry_inp->geo_grid.lower[2] + len/2.0;
        zcut = len/2.0;
        break;
      case GKYL_GEOMETRY_TOKAMAK_CORE_L:
        len = geometry_inp->tok_grid_info.half_domain ? 2.0*(geometry_inp->geo_grid.upper[2] - geometry_inp->geo_grid.lower[2])
                                                      : geometry_inp->geo_grid.upper[2] - geometry_inp->geo_grid.lower[2];
        zcenter = geometry_inp->geo_grid.upper[2] - len/2.0;
        zcut = len/2.0;
        break;
      case GKYL_GEOMETRY_TOKAMAK_PF_LO_R:
      case GKYL_GEOMETRY_TOKAMAK_PF_UP_L:
      case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO:
      case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_UP:
      case GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO:
        len = geometry_inp->geo_grid.upper[2] - geometry_inp->geo_grid.lower[2];
        zcenter = geometry_inp->position_map->xpt_ctx->compress_divertor ? geometry_inp->geo_grid.lower[2] + len/2.0 : geometry_inp->geo_grid.lower[2];
        zcut = geometry_inp->position_map->xpt_ctx->compress_divertor ? len/2.0 : len;
        break;
      case GKYL_GEOMETRY_TOKAMAK_PF_LO_L:
      case GKYL_GEOMETRY_TOKAMAK_PF_UP_R:
      case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP:
      case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO:
      case GKYL_GEOMETRY_TOKAMAK_LSN_SOL_UP:
        len = geometry_inp->geo_grid.upper[2] - geometry_inp->geo_grid.lower[2];
        zcenter = geometry_inp->position_map->xpt_ctx->compress_divertor ? geometry_inp->geo_grid.upper[2] - len/2.0 : geometry_inp->geo_grid.upper[2];
        zcut = geometry_inp->position_map->xpt_ctx->compress_divertor ? len/2.0 : len;
        break;
      default:
        break;
    }

    double w = geometry_inp->geo_grid.upper[0] - geometry_inp->geo_grid.lower[0];
    struct gkyl_efit *efit = gkyl_efit_new(&geometry_inp->efit_info);
    double psisep = efit->psisep;
    gkyl_efit_release(efit);
    gkyl_position_map_set_compression(geometry_inp->position_map, zcut, zcenter, w, psisep);
  }
  else if (geometry_inp->position_map->id == GKYL_PMAP_CONSTANT_DB_POLYNOMIAL || \
           geometry_inp->position_map->id == GKYL_PMAP_CONSTANT_DB_NUMERIC) {
    // First construct the uniform 3d geometry
    struct gkyl_gk_geometry_inp uniform_inp = *geometry_inp;
    uniform_inp.tok_grid_info.relaxed_xpt_seam_optimizer_trial = false;
    uniform_inp.tok_grid_info.relaxed_xpt_seam_trial_status = 0;
    gk_geom_3d = gk_geometry_tok_init(&uniform_inp);
    // The array mc2nu is computed using the uniform geometry, so we need to deflate it
    // Must deflate the 3D uniform geometry in order for the allgather to work
    if(geometry_inp->grid.ndim < 3)
      gk_geom = gkyl_gk_geometry_deflate(gk_geom_3d, geometry_inp);
    else
      gk_geom = gkyl_gk_geometry_acquire(gk_geom_3d);

    gkyl_position_map_set_bmag(geometry_inp->position_map, geometry_inp->comm, \
      gk_geom->geo_int.bmag);

    gkyl_gk_geometry_release(gk_geom_3d); // release temporary 3d geometry
    gkyl_gk_geometry_release(gk_geom); // release 3d geometry
  }
  // Construct the non-uniform grid
  gk_geom_3d = gk_geometry_tok_init(geometry_inp);
  return gk_geom_3d;
}


void
gkyl_gk_geometry_tok_set_grid_extents(struct gkyl_efit_inp efit_info, struct gkyl_tok_geo_grid_inp grid_info, double *theta_lo, double *theta_up) {
  struct gkyl_tok_geo *geo = gkyl_tok_geo_new(&efit_info, &grid_info);
  gkyl_tok_geo_set_extent(&grid_info, geo, theta_lo, theta_up);
  gkyl_tok_geo_release(geo);
}
