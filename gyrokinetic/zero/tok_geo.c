#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <assert.h>
#include <gkyl_basis.h>
#include <gkyl_math.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_nodal_ops.h>
#include <gkyl_position_map.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_priv.h>
#include <gkyl_tok_geo_priv.h>
#include <gkyl_tok_geo_wall_trial_priv.h>
#include <gkyl_tok_geo_wall_target_priv.h>

static _Thread_local bool wall_trial_active;
static _Thread_local long wall_trial_failures;
static _Thread_local int wall_trial_movable_edge;
static _Thread_local bool wall_trial_fixed_violation;
// A node of the NON-movable radial boundary lying outside the vessel is a
// different thing from a fold, a bulge, or a movable boundary that ran out of
// room: nothing the adjustment is permitted to move can fix it, because that
// boundary is the separatrix row the block shares with the core. Tracked
// separately so callers can classify the case instead of reporting a defect.
static _Thread_local bool wall_trial_fixed_node_outside;
// Set when a violation lies anywhere but the movable side of the block -- the
// half of its radial span nearer the edge the trial may move (2026-10-05). A
// trial whose every violation is on that side failed because the wall cuts the
// region being shrunk; anything else -- the fixed side, an interior fold, a
// strike step -- is not something shrinking that edge was ever going to fix.
// A half span, not a cell count: the wall cuts a fixed physical depth, so a
// cell layer would decide differently at each resolution (measured, NSTX-U
// 203532: one cell at coarse, two at phase1).
static _Thread_local bool wall_trial_outside_movable_side;
// A construction failure inside a trial (2026-09-30): the row rule refused a
// block, or the signed-Jacobian guard found a reversal. Both used to end the
// process from inside a trial that was going to be discarded anyway. Recorded
// here instead, so the adjuster can say what it means: a region the wall has
// shrunk until it cannot be built is out of scope (region_degenerate); the
// same at the requested bounds is the construction failing, reported cleanly.
static _Thread_local bool wall_trial_row_rule_refused;
static _Thread_local bool wall_trial_jacobian_invalid;
// A declared plate that is not finite (tok_wall_block_pockets): a
// declaration error, not a boundary any adjustment moves, and not the
// separatrix leaving the machine.
static _Thread_local bool wall_trial_plate_invalid;
static _Thread_local bool wall_trial_capture;
static _Thread_local int wall_trial_block;
static _Thread_local double wall_trial_rho;

void
tok_wall_trial_begin(int movable_radial_edge)
{
  if (wall_trial_active) {
    abort();
  }
  wall_trial_active = true;
  wall_trial_failures = 0;
  wall_trial_movable_edge = movable_radial_edge;
  wall_trial_fixed_violation = false;
  wall_trial_fixed_node_outside = false;
  wall_trial_outside_movable_side = false;
  wall_trial_row_rule_refused = false;
  wall_trial_jacobian_invalid = false;
  wall_trial_plate_invalid = false;
  wall_trial_capture = false;
}

bool
tok_wall_trial_is_active(void)
{
  return wall_trial_active;
}

void
tok_wall_trial_note_row_rule_refused(void)
{
  if (wall_trial_active) {
    wall_trial_row_rule_refused = true;
  }
}

bool
tok_wall_trial_note_jacobian_invalid(void)
{
  if (!wall_trial_active) {
    return false;
  }
  wall_trial_jacobian_invalid = true;
  return true;
}

bool
tok_wall_trial_row_rule_refused(void)
{
  return wall_trial_row_rule_refused;
}

bool
tok_wall_trial_jacobian_invalid(void)
{
  return wall_trial_jacobian_invalid;
}

void
tok_wall_trial_capture_requested(int block, double rho)
{
  if (!wall_trial_active) {
    abort();
  }
  wall_trial_capture = true;
  wall_trial_block = block;
  wall_trial_rho = rho;
}

long
tok_wall_trial_end(void)
{
  if (!wall_trial_active) {
    abort();
  }
  wall_trial_active = false;
  return wall_trial_failures;
}

bool
tok_wall_trial_record(bool fixed_radial_boundary)
{
  return tok_wall_trial_record_scope(fixed_radial_boundary, false);
}

bool
tok_wall_trial_record_scope(bool fixed_radial_boundary, bool node_outside)
{
  return tok_wall_trial_record_where(fixed_radial_boundary, node_outside, false);
}

bool
tok_wall_trial_record_where(bool fixed_radial_boundary, bool node_outside, bool on_movable_side)
{
  if (!wall_trial_active) {
    return false;
  }
  ++wall_trial_failures;
  wall_trial_fixed_violation |= fixed_radial_boundary;
  wall_trial_fixed_node_outside |= fixed_radial_boundary && node_outside;
  wall_trial_outside_movable_side |= !on_movable_side;
  return true;
}

bool
tok_wall_trial_only_movable_side(void)
{
  return wall_trial_failures > 0 && !wall_trial_outside_movable_side;
}

bool
tok_wall_trial_has_fixed_violation(void)
{
  return wall_trial_fixed_violation;
}

bool
tok_wall_trial_has_fixed_node_outside(void)
{
  return wall_trial_fixed_node_outside;
}

void
tok_wall_trial_note_plate_invalid(void)
{
  if (wall_trial_active) {
    wall_trial_plate_invalid = true;
  }
}

bool
tok_wall_trial_plate_invalid(void)
{
  return wall_trial_plate_invalid;
}
#include <gkyl_dg_bin_ops.h>

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static bool
tok_geo_same_flux(double psi_a, double psi_b)
{
  double scale = fmax(1.0, fmax(fabs(psi_a), fabs(psi_b)));
  return fabs(psi_a - psi_b) <= 64.0 * DBL_EPSILON * scale;
}

// Is a point at computational radial coordinate x on the movable side of the
// block -- no farther from the edge the current wall trial may move than from
// the fixed edge? lo and hi are the block's global radial bounds.
static bool
tok_wall_trial_on_movable_side(double x, double lo, double hi)
{
  if (wall_trial_movable_edge < 0) {
    return false;
  }
  const double mid = 0.5 * (lo + hi);
  return wall_trial_movable_edge ? x >= mid || tok_geo_same_flux(x, mid) :
                                   x <= mid || tok_geo_same_flux(x, mid);
}

// Lower bound on the slope of the separatrix->far-boundary correspondence, as
// the weight of an identity map blended into the nearest-point projection.
// Raising it widens the collapsed seam cells that nearest-point projection
// produces near an X point -- 0.25 lifts the thinnest cell of the outboard-SOL
// seam column by an order of magnitude on shots 204965/204995/204997 -- but it
// perturbs the correspondence on every half-domain block, so it stays at the
// validated value.
static double
tok_trace_corr_identity_fraction(void)
{
  return 0.01;
}

static bool
tok_xpt_seam_optimizer_trial(const struct gkyl_tok_geo_grid_inp *inp)
{
  return inp->relaxed_xpt_seam_optimizer_trial && inp->relaxed_xpt_seam_trial_status;
}

static void
tok_xpt_seam_trial_reject(
  const struct gkyl_tok_geo_grid_inp *inp, enum gkyl_tok_geo_xpt_seam_trial_failure reason
)
{
  if (!tok_xpt_seam_optimizer_trial(inp)) {
    return;
  }
  struct gkyl_tok_geo_xpt_seam_trial_status *status = inp->relaxed_xpt_seam_trial_status;
  if (status->first_failure_reason == GKYL_XPT_SEAM_TRIAL_OK) {
    status->first_failure_reason = reason;
  }
  switch (reason) {
    case GKYL_XPT_SEAM_TRIAL_CONTOUR:
      status->contour_valid = false;
      break;
    case GKYL_XPT_SEAM_TRIAL_BRANCH:
      status->branch_valid = false;
      break;
    case GKYL_XPT_SEAM_TRIAL_TRACE_ORDERING:
      status->trace_ordering_valid = false;
      break;
    case GKYL_XPT_SEAM_TRIAL_NONFINITE_MAP:
      status->finite_map_valid = false;
      break;
    case GKYL_XPT_SEAM_TRIAL_CELL_JACOBIAN:
      status->cell_jacobian_valid = false;
      break;
    case GKYL_XPT_SEAM_TRIAL_METRIC_JACOBIAN:
      status->jacobian_valid = false;
      break;
    case GKYL_XPT_SEAM_TRIAL_INVALID_PARAMETER:
    case GKYL_XPT_SEAM_TRIAL_REMOTE_FAILURE:
    case GKYL_XPT_SEAM_TRIAL_OK:
      break;
  }
}

static bool
tok_xpt_mapping_requested(const struct gkyl_tok_geo_grid_inp *inp)
{
  if (inp->straight_xpt_ray) {
    return true;
  }
  return inp->straight_core_xpt_ray &&
         (inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R || inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_L);
}

// Does this block PLACE ITS NODES through the ordered map?  Exactly when it
// requests the X-point mapping.
static bool
tok_xpt_ordered_placement(const struct gkyl_tok_geo_grid_inp *inp)
{
  return tok_xpt_mapping_requested(inp);
}

static void
tok_init_xpt_ray_target(
  const struct gkyl_tok_geo_grid_inp *inp, const struct gkyl_position_map *position_map,
  struct arc_length_ctx *arc_ctx
)
{
  // Retained unconditionally: the theta ladder needs the SAME map that places
  // the node rows, and it needs it whether or not an X-point ray was requested.
  arc_ctx->position_map = position_map;
  if (!tok_xpt_mapping_requested(inp)) {
    return;
  }
  // Flux can increase or decrease radially, and existing inputs use both
  // [separatrix, far] and [far, separatrix] logical orderings.  Select the
  // physical radial boundary furthest from the separatrix instead of
  // inferring it from ftype or array order.
  double psi_comp_lo = inp->cgrid.lower[0];
  double psi_comp_hi = inp->cgrid.upper[0];
  double psi_lo = 0.0, psi_hi = 0.0;
  position_map->maps[0](0.0, &psi_comp_lo, &psi_lo, position_map->ctxs[0]);
  position_map->maps[0](0.0, &psi_comp_hi, &psi_hi, position_map->ctxs[0]);
  double dlo = fabs(psi_lo - arc_ctx->geo->psisep);
  double dhi = fabs(psi_hi - arc_ctx->geo->psisep);
  double scale = fmax(1.0, fmax(fabs(psi_lo), fabs(psi_hi)));
  if (!isfinite(psi_lo) || !isfinite(psi_hi) || fabs(dlo - dhi) <= 64.0 * DBL_EPSILON * scale) {
    fprintf(
      stderr,
      "TOK_ORDERED_MAP cannot identify far radial boundary ftype=%d psi_bounds=(%.17g,%.17g) psi_sep=%.17g\n",
      inp->ftype, psi_lo, psi_hi, arc_ctx->geo->psisep
    );
    abort();
  }
  arc_ctx->xpt_ray_psi0 = dlo > dhi ? psi_lo : psi_hi;
}

static double
tok_xpt_theta_to_arc(
  const struct gkyl_tok_geo_grid_inp *inp, struct arc_length_ctx *arc_ctx, double theta
)
{
  // One block-relative form for both paths. Equal in value to the old absolute
  // expression while the extent is the arc fraction, and correct for any extent.
  if (!arc_ctx->arc_interval_valid) {
    return (theta + M_PI) * arc_ctx->arcL_tot / (2.0 * M_PI);
  }
  double arc = tok_arc_from_theta(inp, arc_ctx, theta);
  if (!arc_ctx->xpt_map_valid) {
    return arc;
  }
  if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R || inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_L) {
    arc = fmod(arc, arc_ctx->arcL_tot);
    if (arc < 0.0) {
      arc += arc_ctx->arcL_tot;
    }
  }
  return arc;
}

static bool
tok_xpt_at_seam(
  const struct gkyl_tok_geo_grid_inp *inp, int it, const struct gkyl_range *nrange,
  bool lower_is_global, bool upper_is_global
)
{
  bool lower = it == nrange->lower[2] && lower_is_global;
  bool upper = it == nrange->upper[2] && upper_is_global;
  switch (inp->ftype) {
    case GKYL_GEOMETRY_TOKAMAK_PF_LO_R:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID:
    case GKYL_GEOMETRY_TOKAMAK_CORE_L:
      return upper;
    case GKYL_GEOMETRY_TOKAMAK_PF_LO_L:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO:
    case GKYL_GEOMETRY_TOKAMAK_CORE_R:
      return lower;
    default:
      return false;
  }
}

static bool
tok_xpt_at_fixed_edge(
  const struct gkyl_tok_geo_grid_inp *inp, const struct arc_length_ctx *arc_ctx, int it,
  const struct gkyl_range *nrange, bool lower_is_global, bool upper_is_global, double *z,
  double *root_reference
)
{
  if (!arc_ctx->xpt_map_valid) {
    return false;
  }
  bool lower = it == nrange->lower[2] && lower_is_global;
  bool upper = it == nrange->upper[2] && upper_is_global;
  switch (inp->ftype) {
    case GKYL_GEOMETRY_TOKAMAK_PF_LO_R:
      if (lower) {
        *z = arc_ctx->zmin_right;
        *root_reference = inp->rright;
        return true;
      }
      break;
    case GKYL_GEOMETRY_TOKAMAK_PF_LO_L:
      if (upper) {
        *z = arc_ctx->zmin_left;
        *root_reference = inp->rleft;
        return true;
      }
      break;
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO:
      if (lower) {
        *z = arc_ctx->zmin;
        *root_reference = inp->rright;
        return true;
      }
      break;
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID:
    case GKYL_GEOMETRY_TOKAMAK_CORE_R:
      if (upper) {
        *z = arc_ctx->geo->zmaxis;
        *root_reference = inp->rright;
        return true;
      }
      break;
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID:
    case GKYL_GEOMETRY_TOKAMAK_CORE_L:
      if (lower) {
        *z = arc_ctx->geo->zmaxis;
        *root_reference = inp->rleft;
        return true;
      }
      break;
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO:
      if (upper) {
        *z = arc_ctx->zmin;
        *root_reference = inp->rleft;
        return true;
      }
      break;
    default:
      break;
  }
  return false;
}

// gkyl_ridders reports an invalid bracket rather than iterating on one, and
// every caller below took its result on faith.  status != 0 means no root of
// the arc-length equation lies in [zmin, zmax] -- which is what a surface that
// no longer reaches this block's boundary looks like, e.g. a SOL surface
// requested past the divertor plate.  The reported root is DBL_MAX there, and
// letting it reach gkyl_tok_geo_R_psiZ indexes psi far outside its range, so
// the run dies in a wild fetch with nothing said about the geometry that
// caused it.  Say what happened instead.
static void
tok_geo_check_arc_root(
  const struct gkyl_tok_geo_grid_inp *inp, double psi, double theta, double arcL, double zmin,
  double zmax, double rid_lo, double rid_hi, const struct gkyl_qr_res *res
)
{
  // NOT keyed on status, and not on isfinite either: gkyl_ridders reports
  // status=2 both for an invalid bracket and for a solve that stopped before
  // meeting eps, and the latter still returns the best root found, which every
  // caller here has always used.  Its unusable value is the DBL_MAX the search
  // starts from -- a finite double, so isfinite() does not see it.  What
  // separates them is that a usable root lies inside the interval it was
  // sought in; the sentinel does not.
  if (res->res >= zmin && res->res <= zmax) {
    return;
  }
  fprintf(
    stderr,
    "TOK_GEO_ARC_ROOT_FAILED ftype=%d psi=%.17g theta=%.17g arcL=%.17g "
    "zbracket=[%.17g,%.17g] f=[%.17g,%.17g] status=%d res=%.17g\n"
    "  no root of the arc-length equation in this block's Z range: the "
    "requested surface does not reach this block's boundary.\n",
    inp->ftype, psi, theta, arcL, zmin, zmax, rid_lo, rid_hi, res->status, res->res
  );
  abort();
}

double
tok_plate_psi_func(double s, void *ctx)
{
  // uses a pointer to the plate function to get R(s), Z(s)
  // Then calculates psi(R, Z)
  // will be used by ridders later

  struct plate_ctx *gc = ctx;
  double RZ[2];
  if (gc->lower == true) {
    gc->geo->plate_func_lower(s, RZ);
  } else {
    gc->geo->plate_func_upper(s, RZ);
  }

  double R = RZ[0];
  double Z = RZ[1];

  // Now find the cell where this R and Z is
  if (gc->geo->use_cubics) {
    double xn[2] = {R, Z};
    double psi;
    gc->geo->efit->evf->eval_cubic(0.0, xn, &psi, gc->geo->efit->evf->ctx);
    return psi - gc->psi_curr;
  } else {
    int rzidx[2];
    rzidx[0] = fmin(
      gc->geo->rzlocal.lower[0] +
        (int)floor((R - gc->geo->rzgrid.lower[0]) / gc->geo->rzgrid.dx[0]),
      gc->geo->rzlocal.upper[0]
    );
    rzidx[1] = fmin(
      gc->geo->rzlocal.lower[1] +
        (int)floor((Z - gc->geo->rzgrid.lower[1]) / gc->geo->rzgrid.dx[1]),
      gc->geo->rzlocal.upper[1]
    );
    long loc = gkyl_range_idx(&gc->geo->rzlocal, rzidx);
    const double *coeffs = gkyl_array_cfetch(gc->geo->psiRZ, loc);

    double xc[2];
    gkyl_rect_grid_cell_center(&gc->geo->rzgrid, rzidx, xc);
    double xy[2];
    xy[0] = (R - xc[0]) / (gc->geo->rzgrid.dx[0] * 0.5);
    xy[1] = (Z - xc[1]) / (gc->geo->rzgrid.dx[1] * 0.5);
    double psi = gc->geo->rzbasis.eval_expand(xy, coeffs);
    return psi - gc->psi_curr;
  }
}

// Function to pass to root-finder to find Z location for given arc-length
static inline double
arc_length_func(double Z, void *ctx)
{
  struct arc_length_ctx *actx = ctx;
  double *arc_memo;
  double psi = actx->psi, rclose = actx->rclose, zmin = actx->zmin, arcL = actx->arcL;
  double zmax = actx->zmax;
  double ival = 0.0;

  if (actx->ftype == GKYL_GEOMETRY_TOKAMAK_CORE) {
    if (actx->right == true) {
      double *arc_memo = actx->arc_memo_right;
      ival =
        integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, true, false, arc_memo) - arcL;
    } else {
      double *arc_memo = actx->arc_memo_left;
      ival = integrate_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, true, false, arc_memo) -
             arcL + actx->arcL_right;
    }
  }

  else if (actx->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_L ||
           actx->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R) {
    if (actx->xpt_map_valid && actx->right) {
      double *arc_memo = actx->arc_memo_right;
      ival =
        integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, true, false, arc_memo) - arcL;
    } else if (actx->xpt_map_valid) {
      double *arc_memo = actx->arc_memo_left;
      ival = integrate_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, true, false, arc_memo) -
             arcL + actx->arcL_right;
    } else if (actx->pre == true) {
      ival = actx->arcL_start -
             integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, false, false, arc_memo) -
             arcL;
    } else if (actx->right == true) {
      double *arc_memo = actx->arc_memo_right;
      ival = integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, true, false, arc_memo) -
             arcL + actx->arcL_start;
    } else {
      double *arc_memo = actx->arc_memo_left;
      ival = integrate_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, true, false, arc_memo) -
             arcL + actx->arcL_right + actx->arcL_start;
    }
  }

  else if (actx->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_L ||
           actx->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_R) {
    if (actx->right == true) {
      double *arc_memo = actx->arc_memo_right;
      ival =
        integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, true, false, arc_memo) - arcL;
    } else {
      double *arc_memo = actx->arc_memo_left;
      ival = integrate_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, true, false, arc_memo) -
             arcL + actx->arcL_right;
    }
  }

  else if (actx->ftype == GKYL_GEOMETRY_TOKAMAK_PF_UP_L ||
           actx->ftype == GKYL_GEOMETRY_TOKAMAK_PF_UP_R) {
    if (actx->right == false) {
      double *arc_memo = actx->arc_memo_left;
      ival =
        integrate_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, true, false, arc_memo) - arcL;
    } else {
      double *arc_memo = actx->arc_memo_right;
      ival = integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, true, false, arc_memo) -
             arcL + actx->arcL_left;
    }
  }

  else if ((actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT) ||
           (actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO) ||
           (actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID) ||
           (actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP)) {
    double *arc_memo = actx->arc_memo;
    ival =
      integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, true, false, arc_memo) - arcL;
  } else if ((actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN) ||
             (actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN) ||
             (actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO) ||
             (actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID) ||
             (actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_UP)) {
    double *arc_memo = actx->arc_memo;
    ival =
      integrate_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, true, false, arc_memo) - arcL;
  }

  else if (actx->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL ||
           actx->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO ||
           actx->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID ||
           actx->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_UP) {
    if (actx->right == true) {
      double *arc_memo = actx->arc_memo_right;
      ival =
        integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, false, false, arc_memo) - arcL;
    } else {
      double *arc_memo = actx->arc_memo_left;
      ival = integrate_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, false, false, arc_memo) -
             arcL + actx->arcL_right;
    }
  }

  else if (actx->ftype == GKYL_GEOMETRY_TOKAMAK_IWL) {
    if (actx->q3) {
      double *arc_memo = actx->arc_memo;
      ival =
        integrate_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, false, false, arc_memo) - arcL;
    } else if (actx->q4) {
      double *arc_memo = actx->arc_memo;
      ival = integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, false, false, arc_memo) -
             arcL + actx->arcL_q3;
    } else if (actx->q1) {
      double *arc_memo = actx->arc_memo;
      ival = integrate_psi_contour_memo(
               actx->geo, psi, actx->geo->zmaxis, Z, rclose, false, false, arc_memo
             ) -
             arcL + actx->arcL_q3 + actx->arcL_q4;
    } else {
      double *arc_memo = actx->arc_memo;
      ival = integrate_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, false, false, arc_memo) -
             arcL + actx->arcL_q3 + actx->arcL_q4 + actx->arcL_q1;
    }
  }

  return ival;
}

static double
tok_eval_psi_rz_local(const struct gkyl_tok_geo *geo, double R, double Z)
{
  if (geo->use_cubics) {
    double xn[2] = {R, Z}, psi = 0.0;
    geo->efit->evf->eval_cubic(0.0, xn, &psi, geo->efit->evf->ctx);
    return psi;
  }

  int idx[2];
  idx[0] = GKYL_MIN2(
    geo->rzlocal.upper[0],
    GKYL_MAX2(
      geo->rzlocal.lower[0],
      geo->rzlocal.lower[0] + (int)floor((R - geo->rzgrid.lower[0]) / geo->rzgrid.dx[0])
    )
  );
  idx[1] = GKYL_MIN2(
    geo->rzlocal.upper[1],
    GKYL_MAX2(
      geo->rzlocal.lower[1],
      geo->rzlocal.lower[1] + (int)floor((Z - geo->rzgrid.lower[1]) / geo->rzgrid.dx[1])
    )
  );
  long loc = gkyl_range_idx(&geo->rzlocal, idx);
  const double *coeffs = gkyl_array_cfetch(geo->psiRZ, loc);
  double xc[2], eta[2];
  gkyl_rect_grid_cell_center(&geo->rzgrid, idx, xc);
  eta[0] = (R - xc[0]) / (0.5 * geo->rzgrid.dx[0]);
  eta[1] = (Z - xc[1]) / (0.5 * geo->rzgrid.dx[1]);
  return geo->rzbasis.eval_expand(eta, coeffs);
}

static void
tok_append_unique_root(double root, double *roots, int *nr, int nmax)
{
  for (int i = 0; i < *nr; ++i) {
    if (fabs(root - roots[i]) <= 2e-11 * fmax(1.0, fabs(root))) {
      return;
    }
  }
  if (*nr < nmax) {
    roots[(*nr)++] = root;
  }
}

// Symmetric companion to gkyl_tok_geo_R_psiZ: find every Z root at fixed R.
// The quadratic tensor representation is solved analytically in each DG cell,
// so horizontal contour segments and Z turning points are not discarded.
static int
tok_geo_Z_psiR(const struct gkyl_tok_geo *geo, double psi, double R, int nmax, double *Z)
{
  int nr = 0;
  if (!geo->use_cubics && geo->efit->rzbasis.poly_order <= 2) {
    int ridx = GKYL_MIN2(
      geo->rzlocal.upper[0],
      GKYL_MAX2(
        geo->rzlocal.lower[0],
        geo->rzlocal.lower[0] + (int)floor((R - geo->rzgrid.lower[0]) / geo->rzgrid.dx[0])
      )
    );
    for (int iz = geo->rzlocal.lower[1]; iz <= geo->rzlocal.upper[1]; ++iz) {
      int idx[2] = {ridx, iz};
      long loc = gkyl_range_idx(&geo->rzlocal, idx);
      const double *p = gkyl_array_cfetch(geo->psiRZ, loc);

      // As in the R scan: this walks every Z cell of the column without an
      // early exit, so skip the cells whose psi enclosure rules out a root.
      if (geo->psi_cell_bounds) {
        const double *pbound = gkyl_array_cfetch(geo->psi_cell_bounds, loc);
        if (psi < pbound[0] || psi > pbound[1]) {
          continue;
        }
      }

      double xc[2];
      gkyl_rect_grid_cell_center(&geo->rzgrid, idx, xc);
      double x = (R - xc[0]) / (0.5 * geo->rzgrid.dx[0]);
      double roots[2] = {0.0, 0.0};
      int ncell = 0;

      if (geo->efit->rzbasis.poly_order == 1) {
        double den = 3.0 * p[3] * x + 1.732050807568877 * p[2];
        double num = -1.732050807568877 * p[1] * x + 2.0 * psi - p[0];
        double scale = fmax(1.0, fabs(num));
        if (fabs(den) > 64.0 * DBL_EPSILON * scale) {
          roots[ncell++] = num / den;
        }
      } else {
        double a = 0.125 * (45.0 * p[8] * x * x + 23.2379000772445 * p[7] * x - 15.0 * p[8] +
                            13.41640786499874 * p[5]);
        double b = 0.125 * (23.2379000772445 * p[6] * x * x + 12.0 * p[3] * x -
                            7.745966692414834 * p[6] + 6.928203230275509 * p[2]);
        double c = 0.125 * ((13.41640786499874 * p[4] - 15.0 * p[8]) * x * x +
                            (6.928203230275509 * p[1] - 7.745966692414834 * p[7]) * x + 5.0 * p[8] -
                            4.47213595499958 * p[4] - 4.47213595499958 * p[5] + 4.0 * p[0]) -
                   psi;
        double coeff_scale = fmax(1.0, fmax(fabs(a), fmax(fabs(b), fabs(c))));
        if (fabs(a) <= 64.0 * DBL_EPSILON * coeff_scale) {
          if (fabs(b) > 64.0 * DBL_EPSILON * coeff_scale) {
            roots[ncell++] = -c / b;
          }
        } else {
          double disc = b * b - 4.0 * a * c;
          double disc_tol = 256.0 * DBL_EPSILON * fmax(1.0, b * b + fabs(4.0 * a * c));
          if (disc >= -disc_tol) {
            disc = fmax(0.0, disc);
            double sd = sqrt(disc);
            if (sd == 0.0) {
              roots[ncell++] = -0.5 * b / a;
            } else {
              double q = -0.5 * (b + copysign(sd, b));
              roots[ncell++] = q / a;
              roots[ncell++] = c / q;
            }
          }
        }
      }

      for (int k = 0; k < ncell; ++k) {
        double y = roots[k];
        if (isfinite(y) && y >= -1.0 - 2e-12 && y <= 1.0 + 2e-12) {
          double z = xc[1] + 0.5 * geo->rzgrid.dx[1] * fmin(1.0, fmax(-1.0, y));
          tok_append_unique_root(z, Z, &nr, nmax);
        }
      }
    }
    return nr;
  }

  // Cubic representation: the per-cell solve of R_psiZ_cubic, transposed. At
  // fixed R, psi restricted to one cell is a cubic in the cell's Z
  // coordinate, whose coefficients are tok_p3_row_coeffs of the transposed
  // expansion at the cell's R coordinate; cells whose enclosure rules the
  // level out are skipped. This replaced the scan below (2026-10-03): 8
  // samples per cell over the whole column plus 64 bisections per bracket,
  // which could not see a root pair inside one sample interval. Measured on
  // NSTX-U 204046 and 203532 the roots agree with it to 8e-12 m and this is
  // about 30x faster. psiRZ_cubic and the evaluator the scan used
  // (efit->evf) are the same polynomial to 1.3e-15 relative.
  if (geo->use_cubics && geo->cubic_transpose_ok) {
    const struct gkyl_rect_grid *g = &geo->rzgrid_cubic;
    const struct gkyl_range *rl = &geo->rzlocal_cubic;
    int ridx = GKYL_MIN2(
      rl->upper[0], GKYL_MAX2(rl->lower[0], rl->lower[0] + (int)floor((R - g->lower[0]) / g->dx[0]))
    );
    for (int iz = rl->lower[1]; iz <= rl->upper[1]; ++iz) {
      int idx[2] = {ridx, iz};
      long loc = gkyl_range_idx(rl, idx);
      if (geo->psi_cell_bounds_cubic) {
        const double *pbound = gkyl_array_cfetch(geo->psi_cell_bounds_cubic, loc);
        if (psi < pbound[0] || psi > pbound[1]) {
          continue;
        }
      }
      const double *p = gkyl_array_cfetch(geo->psiRZ_cubic, loc);
      double pt[16];
      for (int k = 0; k < 16; ++k) {
        pt[k] = p[geo->cubic_transpose[k]];
      }

      double xc[2];
      gkyl_rect_grid_cell_center(g, idx, xc);
      double x = (R - xc[0]) / (0.5 * g->dx[0]);
      double c[4], yr[3];
      tok_p3_row_coeffs(pt, psi, x, c);
      int ncell = tok_cubic_cell_roots(c, tok_p3_coeff_err(p, psi), yr);
      for (int k = 0; k < ncell; ++k) {
        tok_append_unique_root(xc[1] + 0.5 * g->dx[1] * yr[k], Z, &nr, nmax);
      }
    }
    return nr;
  }

  // Cubic fallback when the transposed solve is unavailable: bracket
  // sign-changing roots on a fine Z scan.  Exact X-point endpoints are
  // supplied explicitly by the trace builder.
  const int nsamp = 8 * geo->rzgrid_cubic.cells[1];
  double z0 = geo->rzgrid_cubic.lower[1];
  double f0 = tok_eval_psi_rz_local(geo, R, z0) - psi;
  for (int i = 1; i <= nsamp; ++i) {
    double z1 = geo->rzgrid_cubic.lower[1] +
                (geo->rzgrid_cubic.upper[1] - geo->rzgrid_cubic.lower[1]) * i / (double)nsamp;
    double f1 = tok_eval_psi_rz_local(geo, R, z1) - psi;
    if (isfinite(f0) && isfinite(f1) && f0 * f1 <= 0.0) {
      double lo = z0, hi = z1, flo = f0;
      for (int n = 0; n < 64; ++n) {
        double mid = 0.5 * (lo + hi);
        double fm = tok_eval_psi_rz_local(geo, R, mid) - psi;
        if (flo * fm <= 0.0) {
          hi = mid;
        } else {
          lo = mid;
          flo = fm;
        }
      }
      tok_append_unique_root(0.5 * (lo + hi), Z, &nr, nmax);
    }
    z0 = z1;
    f0 = f1;
  }
  return nr;
}

static double
tok_nearest_value(double ref, const double *values, int n)
{
  int ibest = 0;
  for (int i = 1; i < n; ++i) {
    if (fabs(values[i] - ref) < fabs(values[ibest] - ref)) {
      ibest = i;
    }
  }
  return values[ibest];
}

static bool tok_fixed_edge_is_midplane(enum gkyl_tok_geo_type ftype);

// Decide whether one trace step joins two points along a resolvable arc of the
// same contour, or bridges two distinct root branches, by re-parameterizing that
// single step in the *other* coordinate.  A connected arc resolves into small,
// uniform sub-steps; a branch jump leaves the sub-trace stranded on the starting
// branch and shows up as one sub-step carrying almost the whole separation.
// Returns the sub-arc length in *arc so the caller can report how much contour
// the parent step skipped.
static bool
tok_step_is_connected_arc(
  const struct gkyl_tok_geo *geo, double psi, bool param_is_r, double r0, double z0, double r1,
  double z1, double *arc
)
{
  const int nsub = 32;
  double rp = r0, zp = z0, total = 0.0, max_sub = 0.0;
  for (int i = 1; i <= nsub; ++i) {
    double t = i / (double)nsub;
    double rc, zc;
    if (i == nsub) {
      rc = r1;
      zc = z1;
    } else if (param_is_r) {
      // The parent step swept R, so sweep Z across the same interval.
      zc = z0 + t * (z1 - z0);
      double roots[8] = {0.0}, dRdZ[8] = {0.0};
      double dR[8] = {0.0}, dZ[8] = {0.0};
      int nr = gkyl_tok_geo_R_psiZ(geo, psi, zc, 8, roots, dRdZ, dR, dZ);
      if (nr == 0) {
        return false;
      }
      rc = tok_nearest_value(rp, roots, nr);
    } else {
      rc = r0 + t * (r1 - r0);
      double roots[16] = {0.0};
      int nr = tok_geo_Z_psiR(geo, psi, rc, 16, roots);
      if (nr == 0) {
        return false;
      }
      zc = tok_nearest_value(zp, roots, nr);
    }
    if (!isfinite(rc) || !isfinite(zc)) {
      return false;
    }
    double residual = tok_eval_psi_rz_local(geo, rc, zc) - psi;
    if (!isfinite(residual) || fabs(residual) > 1e-9 * fmax(1.0, fabs(psi))) {
      return false;
    }
    double ds = hypot(rc - rp, zc - zp);
    if (!isfinite(ds)) {
      return false;
    }
    total += ds;
    max_sub = fmax(max_sub, ds);
    rp = rc;
    zp = zc;
  }
  *arc = total;
  // Uniformity of the refinement is the discriminator, not its absolute size: a
  // resolved arc spreads evenly over the sub-steps, so no sub-step may carry
  // more than 4x the average.  A branch jump concentrates the whole separation
  // into one sub-step and cannot meet this however fine the refinement.
  return isfinite(total) && total > 0.0 && max_sub <= 4.0 * total / nsub;
}

static bool
tok_build_contour_candidate(
  const struct gkyl_tok_geo *geo, double psi, bool param_is_r, double rfixed, double zfixed,
  double rx, double zx, int n, double *r, double *z, double *score, const char *name,
  bool cluster_fixed_endpoint
)
{
  r[0] = rfixed;
  z[0] = zfixed;
  double total = 0.0, max_step = 0.0, final_step = 0.0;
  double first_step = 0.0, max_step_after_first = 0.0;
  for (int i = 1; i < n - 1; ++i) {
    double t = i / (double)(n - 1);
    // A midplane endpoint is an R turning point of the contour.  There,
    // Z-Z_fixed scales as sqrt(|R-R_fixed|), so uniform R sampling creates an
    // artificially large first physical step.  Quadratic R spacing recovers
    // approximately uniform arc-length spacing without weakening the checks
    // that reject a genuine root-branch jump.
    double f = param_is_r && cluster_fixed_endpoint ? t * t : t;
    if (param_is_r) {
      r[i] = rfixed + f * (rx - rfixed);
      double roots[16] = {0.0};
      int nr = tok_geo_Z_psiR(geo, psi, r[i], 16, roots);
      if (nr == 0) {
        return false;
      }
      z[i] = tok_nearest_value(z[i - 1], roots, nr);
    } else {
      z[i] = zfixed + f * (zx - zfixed);
      double roots[8] = {0.0}, dRdZ[8] = {0.0};
      double dR[8] = {0.0}, dZ[8] = {0.0};
      int nr = gkyl_tok_geo_R_psiZ(geo, psi, z[i], 8, roots, dRdZ, dR, dZ);
      if (nr == 0) {
        return false;
      }
      r[i] = tok_nearest_value(r[i - 1], roots, nr);
    }
    if (!isfinite(r[i]) || !isfinite(z[i])) {
      return false;
    }
    double residual = tok_eval_psi_rz_local(geo, r[i], z[i]) - psi;
    if (!isfinite(residual) || fabs(residual) > 1e-9 * fmax(1.0, fabs(psi))) {
      return false;
    }
    double ds = hypot(r[i] - r[i - 1], z[i] - z[i - 1]);
    if (!(ds > 0.0) || !isfinite(ds)) {
      return false;
    }
    total += ds;
    max_step = fmax(max_step, ds);
    if (i == 1) {
      first_step = ds;
    } else {
      max_step_after_first = fmax(max_step_after_first, ds);
    }
  }
  r[n - 1] = rx;
  z[n - 1] = zx;
  final_step = hypot(rx - r[n - 2], zx - z[n - 2]);
  if (!(final_step > 0.0) || !isfinite(final_step)) {
    return false;
  }
  total += final_step;
  max_step = fmax(max_step, final_step);
  if (n > 2) {
    max_step_after_first = fmax(max_step_after_first, final_step);
  } else {
    first_step = final_step;
  }
  double mean = total / (n - 1);
  double cell_diag = hypot(geo->rzgrid.dx[0], geo->rzgrid.dx[1]);
  // max_step/mean is a *relative* uniformity test standing in for "the trace jumped
  // to another root branch".  Arc-length spacing along a contour varies strongly near
  // an X point, so a perfectly smooth trace can exceed 16x mean while every step stays
  // far below the grid scale -- these traces were being rejected with max_step 3-20x
  // *under* the absolute threshold.  A jump to a distinct branch has to move a
  // meaningful fraction of a cell, so only let the ratio test fire once the step is
  // also large in absolute terms.  The two absolute tests are unchanged.
  const double ratio_test_floor = 0.5 * cell_diag;
  // The clustered retry treats a midplane fixed edge as an R extremum of the
  // contour.  That holds unless the equilibrium was reflected about Z=0 while its
  // magnetic axis sits elsewhere (efit.c forces zmaxis=0): the mirrored pair of
  // axes then leaves a cusp at the seam, the surface has a local R *minimum*
  // there, and an R sweep inward skips the small sub-arc bulging outboard of it.
  // No amount of endpoint clustering shrinks that first step, because the arc it
  // skips is not parameterizable in R at all.  Rather than loosen the constant,
  // measure the step: re-parameterize it in Z and see whether it is one connected
  // arc.  If it is, judge uniformity on the steps that remain.  Only a trace
  // whose plain R and Z parameterizations were both already rejected is ever
  // clustered, so this cannot change a trace that is accepted today.
  double ratio_step = max_step, ratio_mean = mean;
  double skipped_arc = 0.0;
  bool first_step_verified = false;
  if (cluster_fixed_endpoint && n > 2 && first_step >= max_step &&
      tok_step_is_connected_arc(geo, psi, param_is_r, rfixed, zfixed, r[1], z[1], &skipped_arc)) {
    ratio_step = max_step_after_first;
    ratio_mean = (total - first_step) / (n - 2);
    first_step_verified = true;
  }
  if (max_step > 2.0 * cell_diag ||
      (ratio_step > 16.0 * ratio_mean && ratio_step > ratio_test_floor) ||
      final_step > 2.0 * cell_diag) {
    return false;
  }
  if (first_step_verified) {
    fprintf(
      stderr,
      "TOK_TRACE_CUSP_STEP trace=%s ftype_param=%c psi=%.17g first_step=%.17g "
      "sub_arc=%.17g excess=%.17g max_after_first=%.17g mean=%.17g "
      "fixed=(%.17g,%.17g)\n",
      name, param_is_r ? 'R' : 'Z', psi, first_step, skipped_arc, skipped_arc - first_step,
      max_step_after_first, mean, rfixed, zfixed
    );
  }
  *score = max_step / mean + 4.0 * final_step / mean;
  return isfinite(*score);
}

static bool
tok_sep_fixed_edge_is_first(enum gkyl_tok_geo_type ftype)
{
  return ftype == GKYL_GEOMETRY_TOKAMAK_CORE_L || ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID ||
         ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_R || ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO;
}

static bool tok_fixed_edge_is_midplane(enum gkyl_tok_geo_type ftype);

static int tok_ext_map_trace_request(const struct gkyl_tok_geo_grid_inp *inp);

// Size of every trace buffer of a block (separatrix, far boundary, domain,
// map). 16 samples per equilibrium cell is the shipped allocation; the map
// trace's request is added so it is never clipped. With 16*nzcells alone a
// block with more than 4*nzcells theta cells had its map trace silently
// clipped (NSTX-U's 65-row equilibrium: from theta x16 on).
static int
tok_sep_trace_capacity(const struct gkyl_tok_geo_grid_inp *inp, int nzcells)
{
  return GKYL_MAX2(16 * nzcells + 1, tok_ext_map_trace_request(inp));
}

// Reference traces (separatrix, far boundary, domain) carry 257 nodes, clamped
// to the buffers the caller allocated (sep_trace_capacity). The count is the
// same for every block, so two blocks meeting on a shared contour sample it
// identically.
static int
tok_reference_trace_nodes(const struct gkyl_tok_geo_grid_inp *inp, int capacity)
{
  (void)inp;
  return GKYL_MIN2(257, capacity);
}

static void
tok_build_sep_trace(
  const struct gkyl_tok_geo_grid_inp *inp, struct arc_length_ctx *arc_ctx, double zfixed,
  double rfixed_ref
)
{
  if (arc_ctx->sep_trace_initialized) {
    return;
  }
  if (!arc_ctx->sep_trace_r || !arc_ctx->sep_trace_z || !arc_ctx->sep_trace_s ||
      arc_ctx->sep_trace_capacity < 17) {
    fprintf(stderr, "TOK_SEP_TRACE missing storage ftype=%d\n", inp->ftype);
    abort();
  }
  double R[8] = {0.0}, dRdZ[8] = {0.0};
  double dR[8] = {0.0}, dZ[8] = {0.0};
  int nr = gkyl_tok_geo_R_psiZ(arc_ctx->geo, arc_ctx->geo->psisep, zfixed, 8, R, dRdZ, dR, dZ);
  if (nr == 0) {
    fprintf(
      stderr, "TOK_SEP_TRACE no fixed-edge root ftype=%d z=%.17g psi=%.17g\n", inp->ftype, zfixed,
      arc_ctx->geo->psisep
    );
    abort();
  }
  double rfixed = tok_nearest_value(rfixed_ref, R, nr);
  double rx = arc_ctx->geo->use_cubics ? arc_ctx->geo->efit->Rxpt_cubic[0] :
                                         arc_ctx->geo->efit->Rxpt[0];
  double zx = arc_ctx->geo->use_cubics ? arc_ctx->geo->efit->Zxpt_cubic[0] :
                                         arc_ctx->geo->efit->Zxpt[0];
  // A 257-point trace resolves a quadratic EFIT cell much more finely than
  // the geometry grid while avoiding the large startup penalty of rebuilding
  // 1025-point candidate traces in each geometry pass.
  int n = tok_reference_trace_nodes(inp, arc_ctx->sep_trace_capacity);
  double *rr = gkyl_malloc(sizeof(double[n]));
  double *zr = gkyl_malloc(sizeof(double[n]));
  double *rz = gkyl_malloc(sizeof(double[n]));
  double *zz = gkyl_malloc(sizeof(double[n]));
  double score_r = DBL_MAX, score_z = DBL_MAX;
  bool ok_r = tok_build_contour_candidate(
    arc_ctx->geo, arc_ctx->geo->psisep, true, rfixed, zfixed, rx, zx, n, rr, zr, &score_r,
    "separatrix", false
  );
  bool ok_z = tok_build_contour_candidate(
    arc_ctx->geo, arc_ctx->geo->psisep, false, rfixed, zfixed, rx, zx, n, rz, zz, &score_z,
    "separatrix", false
  );
  // Same endpoint-clustered R retry the far-boundary trace already uses. A
  // midplane fixed edge is an R turning point of the contour: Z-Z_fixed scales
  // as sqrt(|R-R_fixed|) there, so uniform R sampling makes the *first* physical
  // step artificially large and trips the discontinuity test even though the
  // contour is perfectly smooth (204502: max_step is step 1 of 256, with
  // neighbours decaying 0.086/0.043/0.026/0.019 -- a turning point, not a branch
  // jump, which would be an isolated spike like this shot's param=Z at 249/256).
  // Only a trace for which both parameterizations were already rejected gets the
  // retry, so every accepted trace is unchanged.
  if (!ok_r && !ok_z && tok_fixed_edge_is_midplane(inp->ftype)) {
    ok_r = tok_build_contour_candidate(
      arc_ctx->geo, arc_ctx->geo->psisep, true, rfixed, zfixed, rx, zx, n, rr, zr, &score_r,
      "separatrix", true
    );
  }
  if (!ok_r && !ok_z) {
    fprintf(
      stderr,
      "TOK_SEP_TRACE both parameterizations failed ftype=%d fixed=(%.17g,%.17g) xpt=(%.17g,%.17g)\n",
      inp->ftype, rfixed, zfixed, rx, zx
    );
    abort();
  }
  bool use_r = ok_r && (!ok_z || score_r <= score_z);
  const double *src_r = use_r ? rr : rz;
  const double *src_z = use_r ? zr : zz;
  bool fixed_first = tok_sep_fixed_edge_is_first(inp->ftype);
  for (int i = 0; i < n; ++i) {
    int src = fixed_first ? i : n - 1 - i;
    arc_ctx->sep_trace_r[i] = src_r[src];
    arc_ctx->sep_trace_z[i] = src_z[src];
  }
  arc_ctx->sep_trace_s[0] = 0.0;
  for (int i = 1; i < n; ++i) {
    arc_ctx->sep_trace_s[i] =
      arc_ctx->sep_trace_s[i - 1] + hypot(
                                      arc_ctx->sep_trace_r[i] - arc_ctx->sep_trace_r[i - 1],
                                      arc_ctx->sep_trace_z[i] - arc_ctx->sep_trace_z[i - 1]
                                    );
  }
  if (!(arc_ctx->sep_trace_s[n - 1] > 0.0) || !isfinite(arc_ctx->sep_trace_s[n - 1])) {
    fprintf(stderr, "TOK_SEP_TRACE invalid arc length ftype=%d\n", inp->ftype);
    abort();
  }
  arc_ctx->sep_trace_n = n;
  arc_ctx->sep_trace_param_is_r = use_r;
  arc_ctx->sep_trace_initialized = true;
  gkyl_free(rr);
  gkyl_free(zr);
  gkyl_free(rz);
  gkyl_free(zz);
}

struct tok_ordered_point {
  double r, z, phi;
  double dr_dtheta, dz_dtheta, dphi_dtheta;
};

static bool
tok_fixed_edge_is_midplane(enum gkyl_tok_geo_type ftype)
{
  return ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R || ftype == GKYL_GEOMETRY_TOKAMAK_CORE_L ||
         ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID ||
         ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID;
}

static bool
tok_fixed_edge_uses_lower_plate(enum gkyl_tok_geo_type ftype)
{
  return ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_R || ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO ||
         ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO;
}

// Which side of the machine a plate-anchored block's strike point sits on.  The two
// halves of a private-flux region share one plate, so a flux surface crossing it twice
// yields one root per half; this picks the half that belongs to this block.
static bool
tok_plate_edge_is_outboard(enum gkyl_tok_geo_type ftype)
{
  return ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_R || ftype == GKYL_GEOMETRY_TOKAMAK_PF_UP_R ||
         ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO || ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO;
}

static bool
tok_plate_flux_intersection(
  const struct gkyl_tok_geo *geo, plate_func plate, double psi, int ftype, double *r, double *z
)
{
  int wall_slot = tok_divertor_wall_slot(geo, plate);
  if (wall_slot != -1) {
    return tok_divertor_wall_intersection(geo, wall_slot, psi, r, z);
  }
  if (!plate) {
    return false;
  }

  const int nsamp = TOK_PLATE_NSAMP;
  double roots_s[16] = {0.0}, roots_r[16] = {0.0};
  double roots_z[16] = {0.0};
  int nroots = 0;
  double rz0[2] = {0.0}, rz1[2] = {0.0};
  plate(0.0, rz0);
  double f0 = tok_eval_psi_rz_local(geo, rz0[0], rz0[1]) - psi;
  double flux_tol = 1e-10 * fmax(1.0, fabs(psi));
  if (isfinite(f0) && fabs(f0) <= flux_tol) {
    roots_s[nroots] = 0.0;
    roots_r[nroots] = rz0[0];
    roots_z[nroots++] = rz0[1];
  }
  for (int i = 1; i <= nsamp; ++i) {
    double s1 = i / (double)nsamp;
    plate(s1, rz1);
    double f1 = tok_eval_psi_rz_local(geo, rz1[0], rz1[1]) - psi;
    bool endpoint_root = isfinite(f1) && fabs(f1) <= flux_tol;
    bool bracket = isfinite(f0) && isfinite(f1) && f0 * f1 < 0.0;
    if ((endpoint_root || bracket) && nroots < 16) {
      double sr = s1;
      if (!endpoint_root) {
        double slo = (i - 1) / (double)nsamp, shi = s1, flo = f0;
        for (int k = 0; k < 64; ++k) {
          double smid = 0.5 * (slo + shi), rzm[2];
          plate(smid, rzm);
          double fm = tok_eval_psi_rz_local(geo, rzm[0], rzm[1]) - psi;
          if (!isfinite(fm)) {
            return false;
          }
          if (flo * fm <= 0.0) {
            shi = smid;
          } else {
            slo = smid;
            flo = fm;
          }
        }
        sr = 0.5 * (slo + shi);
      }
      double rzr[2];
      plate(sr, rzr);
      if (nroots == 0 || fabs(sr - roots_s[nroots - 1]) > 2e-10) {
        roots_s[nroots] = sr;
        roots_r[nroots] = rzr[0];
        roots_z[nroots++] = rzr[1];
      }
    }
    rz0[0] = rz1[0];
    rz0[1] = rz1[1];
    f0 = f1;
  }
  if (nroots == 0) {
    if (tok_limiter_plate_intersection(geo, plate, psi, r, z)) {
      return true;
    }
    fprintf(stderr, "TOK_ORDERED_MAP plate root count=0 ftype=%d psi=%.17g\n", ftype, psi);
    return false;
  }
  int pick = 0;
  if (nroots > 1) {
    // The two halves of a private-flux region share this plate, so a surface that
    // crosses it twice gives one root per half.  Select the half this block owns
    // using the same side hint the midplane path uses, rather than rejecting.
    double hint = tok_plate_edge_is_outboard(ftype) ? geo->rright : geo->rleft;
    double best = fabs(roots_r[0] - hint);
    for (int k = 1; k < nroots; ++k) {
      double d = fabs(roots_r[k] - hint);
      if (d < best) {
        best = d;
        pick = k;
      }
    }
    fprintf(
      stderr,
      "TOK_ORDERED_MAP plate root count=%d ftype=%d psi=%.17g hint=%.17g picked=%d rz=(%.17g,%.17g)\n",
      nroots, ftype, psi, hint, pick, roots_r[pick], roots_z[pick]
    );
  }
  *r = roots_r[pick];
  *z = roots_z[pick];
  double residual = tok_eval_psi_rz_local(geo, *r, *z) - psi;
  return isfinite(*r) && isfinite(*z) && isfinite(residual) &&
         fabs(residual) <= 1e-9 * fmax(1.0, fabs(psi));
}

// Locate the physical (plate or midplane) end of a block at a specified psi.
// Plate intersections are scanned before refinement so an exact endpoint is
// returned directly and an ambiguous multi-intersection plate is rejected.
static bool
tok_fixed_edge_point(
  const struct gkyl_tok_geo_grid_inp *inp, const struct gkyl_tok_geo *geo, double psi, double *r,
  double *z
)
{
  if (tok_fixed_edge_is_midplane(inp->ftype)) {
    double R[16] = {0.0}, dRdZ[16] = {0.0};
    double dR[16] = {0.0}, dZ[16] = {0.0};
    *z = geo->zmaxis;
    int nr = gkyl_tok_geo_R_psiZ(geo, psi, *z, 8, R, dRdZ, dR, dZ);
    if (nr <= 0) {
      return false;
    }
    bool outboard = inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R ||
                    inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID;
    *r = tok_nearest_value(outboard ? inp->rright : inp->rleft, R, nr);
    return isfinite(*r);
  }

  plate_func plate = 0;
  if (tok_fixed_edge_uses_lower_plate(inp->ftype)) {
    plate = geo->plate_func_lower;
  } else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_L) {
    plate = geo->plate_func_upper;
  }
  return tok_plate_flux_intersection(geo, plate, psi, inp->ftype, r, z);
}

enum tok_ext_endpoint_kind {
  TOK_EXT_XPT_RAY,
  TOK_EXT_PLATE,
  // Half-domain MID/CORE blocks stop at the midplane symmetry plane rather
  // than at the far X point, so their topology needs an endpoint that is
  // neither a ray nor a material surface: the point where this psi contour
  // meets Z = zmaxis on the given side.
  TOK_EXT_MIDPLANE
};

enum tok_ext_xpoint { TOK_EXT_LOWER_XPT, TOK_EXT_UPPER_XPT };

enum tok_ext_sector { TOK_EXT_CORE, TOK_EXT_PF, TOK_EXT_SOL_OUT, TOK_EXT_SOL_IN };

enum tok_ext_plate_slot { TOK_EXT_PLATE_LOWER, TOK_EXT_PLATE_UPPER };

enum tok_ext_fixed_z_slot {
  TOK_EXT_ZMIN,
  TOK_EXT_ZMAX,
  TOK_EXT_ZMIN_LEFT,
  TOK_EXT_ZMIN_RIGHT,
  TOK_EXT_ZMAX_LEFT,
  TOK_EXT_ZMAX_RIGHT
};

enum tok_ext_route {
  TOK_EXT_ROUTE_GENERIC,
  TOK_EXT_ROUTE_OUTBOARD,
  TOK_EXT_ROUTE_INBOARD,
  TOK_EXT_ROUTE_VIA_UPPER,
  TOK_EXT_ROUTE_CLOSED_CORE,
  // Half of a closed core surface, cut at two X-point rays. See
  // tok_ext_build_core_half_trace.
  TOK_EXT_ROUTE_CORE_HALF
};

enum tok_ext_phi_reference {
  TOK_EXT_PHI_LOWER,
  TOK_EXT_PHI_UPPER,
  TOK_EXT_PHI_OUTBOARD_MIDPLANE,
  TOK_EXT_PHI_INBOARD_MIDPLANE
};

struct tok_ext_endpoint {
  enum tok_ext_endpoint_kind kind;
  enum tok_ext_xpoint xpoint;
  enum tok_ext_sector sector;
  enum tok_ext_plate_slot plate_slot;
  enum tok_ext_fixed_z_slot fixed_z_slot;
  // TOK_EXT_MIDPLANE only. CORE_R and CORE_L share a sector but sit on
  // opposite sides, so the side cannot be derived from `sector`.
  bool midplane_outboard;
};

struct tok_ext_topology {
  struct tok_ext_endpoint lower, upper;
  enum tok_ext_route route;
  enum tok_ext_phi_reference phi_reference;
  bool closed;
  // Both ends of this block are X-point rays on one closed surface, so its
  // two halves share both endpoints and the seed cannot name the branch:
  // only the block itself can. See tok_ext_z_branch_points.
  bool seed_names_side;
};

static struct tok_ext_endpoint
tok_ext_xray(enum tok_ext_xpoint xpoint, enum tok_ext_sector sector)
{
  return (struct tok_ext_endpoint){.kind = TOK_EXT_XPT_RAY, .xpoint = xpoint, .sector = sector};
}

static struct tok_ext_endpoint
tok_ext_plate(enum tok_ext_plate_slot plate_slot, enum tok_ext_fixed_z_slot fixed_z_slot)
{
  return (struct tok_ext_endpoint){
    .kind = TOK_EXT_PLATE,
    .plate_slot = plate_slot,
    .fixed_z_slot = fixed_z_slot,
  };
}

static struct tok_ext_endpoint
tok_ext_midplane(bool outboard)
{
  return (struct tok_ext_endpoint){.kind = TOK_EXT_MIDPLANE, .midplane_outboard = outboard};
}

// Logical theta-lower -> theta-upper topology for the double-null and
// lower-single-null multiblock types.
//
// `half_domain` models only the lower half, closing the domain on the midplane
// symmetry plane. That changes the topology of exactly four blocks: the two
// MID SOL blocks and the two CORE halves stop at the midplane instead of
// running on to the far X point (see tok_geo_set_extent, where those four are
// the only ftypes whose theta extents carry a `half_domain` branch). The four
// leg blocks -- PF_LO_R, PF_LO_L, DN_SOL_OUT_LO, DN_SOL_IN_LO -- have
// identical extents in both modes and so share a single entry here.
static bool
tok_ext_topology_from_ftype_raw(
  enum gkyl_tok_geo_type ftype, bool half_domain, struct tok_ext_topology *top
)
{
  *top = (struct tok_ext_topology){};
  if (half_domain) {
    switch (ftype) {
      // A mid SOL block and the core block across the separatrix from it share
      // that whole arc, and on the separatrix they share its ENDPOINTS too --
      // the X point and the same midplane root.  What they did not share was
      // the CONSTRUCTION: the SOL side ran the Z-branch route while the core
      // side runs GENERIC (see the CORE_R note below, which is load-bearing and
      // must not be undone).  Two different builders on one curve agree only to
      // trace tolerance -- measured 7.9e-06 in normalized arc length on
      // b2<->b6, against 2.2e-15 for the leg pairs, which DO share a route.
      // Arc length hides that at ~1e-05 m, but dl/dmu = 1/|grad psi| amplifies
      // it ~60x at the X point, which is the ~8e-04 m core<->SOL seam residue.
      //
      // Putting the mid SOL blocks on GENERIC too makes both sides run the same
      // builder, so the two traces are computed INDEPENDENTLY and still agree
      // exactly -- the same guarantee tok_ext_build_via_turning_trace already
      // gives the single-null LSN_SOL_MID/full-core pair.  Generic is already
      // well exercised here: ftype 6 falls back to it on its own via
      // ext_force_generic_route on many NSTX-U shots.
      case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID:
        top->lower = tok_ext_xray(TOK_EXT_LOWER_XPT, TOK_EXT_SOL_OUT);
        top->upper = tok_ext_midplane(true);
        top->route = TOK_EXT_ROUTE_GENERIC;
        top->phi_reference = TOK_EXT_PHI_OUTBOARD_MIDPLANE;
        return true;
      case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID:
        top->lower = tok_ext_midplane(false);
        top->upper = tok_ext_xray(TOK_EXT_LOWER_XPT, TOK_EXT_SOL_IN);
        top->route = TOK_EXT_ROUTE_GENERIC;
        top->phi_reference = TOK_EXT_PHI_INBOARD_MIDPLANE;
        return true;
      // The core X-point ray endpoint on an INTERIOR surface is the nearest
      // point on that surface to the X point, which lies on its inboard-lower
      // flank rather than at its extreme Z. CORE_R therefore runs inboard,
      // under the surface's LOWER turning point, and up to the outboard
      // midplane -- a path where no single parameterization works: Z reverses
      // at the turning point and R is degenerate at the midplane end. Split it
      // gate by 6e-5 relative -- so give the generic route the endpoint-
      // clustered retry (see tok_ext_build_open_trace) and it gets through.
      //
      // DO NOT RETRY the obvious-looking alternative of splitting at the lower
      // turning point (a mirrored tok_ext_build_via_turning_trace). Measured on
      // all 20 shots that need this: 15 built but ALL 15 folded CORE_R by 8-10
      // cells with minA/medA ~ -6e3, versus 1-3 cells at ~-1e-1 for the
      // clustered generic route. Uniform-Z sampling crowds the arc against the
      // turning point, and resampling the joined polyline by arc length just
      // interpolates across that chord instead of following the contour.
      case GKYL_GEOMETRY_TOKAMAK_CORE_R:
        top->lower = tok_ext_xray(TOK_EXT_LOWER_XPT, TOK_EXT_CORE);
        top->upper = tok_ext_midplane(true);
        top->route = TOK_EXT_ROUTE_GENERIC;
        top->phi_reference = TOK_EXT_PHI_OUTBOARD_MIDPLANE;
        return true;
      case GKYL_GEOMETRY_TOKAMAK_CORE_L:
        top->lower = tok_ext_midplane(false);
        top->upper = tok_ext_xray(TOK_EXT_LOWER_XPT, TOK_EXT_CORE);
        top->route = TOK_EXT_ROUTE_GENERIC;
        top->phi_reference = TOK_EXT_PHI_INBOARD_MIDPLANE;
        return true;
      default:
        break; // leg blocks fall through to the shared table below
    }
  }
  switch (ftype) {
    case GKYL_GEOMETRY_TOKAMAK_PF_LO_R:
      top->lower = tok_ext_plate(TOK_EXT_PLATE_LOWER, TOK_EXT_ZMIN_RIGHT);
      top->upper = tok_ext_xray(TOK_EXT_LOWER_XPT, TOK_EXT_PF);
      top->route = TOK_EXT_ROUTE_GENERIC;
      top->phi_reference = TOK_EXT_PHI_LOWER;
      return true;
    case GKYL_GEOMETRY_TOKAMAK_PF_LO_L:
      top->lower = tok_ext_xray(TOK_EXT_LOWER_XPT, TOK_EXT_PF);
      top->upper = tok_ext_plate(TOK_EXT_PLATE_UPPER, TOK_EXT_ZMIN_LEFT);
      top->route = TOK_EXT_ROUTE_GENERIC;
      top->phi_reference = TOK_EXT_PHI_UPPER;
      return true;
    case GKYL_GEOMETRY_TOKAMAK_PF_UP_R:
      top->lower = tok_ext_xray(TOK_EXT_UPPER_XPT, TOK_EXT_PF);
      top->upper = tok_ext_plate(TOK_EXT_PLATE_UPPER, TOK_EXT_ZMAX_RIGHT);
      top->route = TOK_EXT_ROUTE_GENERIC;
      top->phi_reference = TOK_EXT_PHI_UPPER;
      return true;
    case GKYL_GEOMETRY_TOKAMAK_PF_UP_L:
      top->lower = tok_ext_plate(TOK_EXT_PLATE_LOWER, TOK_EXT_ZMAX_LEFT);
      top->upper = tok_ext_xray(TOK_EXT_UPPER_XPT, TOK_EXT_PF);
      top->route = TOK_EXT_ROUTE_GENERIC;
      top->phi_reference = TOK_EXT_PHI_LOWER;
      return true;
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO:
      top->lower = tok_ext_plate(TOK_EXT_PLATE_LOWER, TOK_EXT_ZMIN);
      top->upper = tok_ext_xray(TOK_EXT_LOWER_XPT, TOK_EXT_SOL_OUT);
      top->route = TOK_EXT_ROUTE_GENERIC;
      top->phi_reference = TOK_EXT_PHI_LOWER;
      return true;
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID:
      top->lower = tok_ext_xray(TOK_EXT_LOWER_XPT, TOK_EXT_SOL_OUT);
      top->upper = tok_ext_xray(TOK_EXT_UPPER_XPT, TOK_EXT_SOL_OUT);
      top->route = TOK_EXT_ROUTE_OUTBOARD;
      top->phi_reference = TOK_EXT_PHI_OUTBOARD_MIDPLANE;
      return true;
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP:
      top->lower = tok_ext_xray(TOK_EXT_UPPER_XPT, TOK_EXT_SOL_OUT);
      top->upper = tok_ext_plate(TOK_EXT_PLATE_UPPER, TOK_EXT_ZMAX);
      top->route = TOK_EXT_ROUTE_GENERIC;
      top->phi_reference = TOK_EXT_PHI_UPPER;
      return true;
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_UP:
      top->lower = tok_ext_plate(TOK_EXT_PLATE_UPPER, TOK_EXT_ZMAX);
      top->upper = tok_ext_xray(TOK_EXT_UPPER_XPT, TOK_EXT_SOL_IN);
      top->route = TOK_EXT_ROUTE_GENERIC;
      top->phi_reference = TOK_EXT_PHI_LOWER;
      return true;
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID:
      top->lower = tok_ext_xray(TOK_EXT_UPPER_XPT, TOK_EXT_SOL_IN);
      top->upper = tok_ext_xray(TOK_EXT_LOWER_XPT, TOK_EXT_SOL_IN);
      top->route = TOK_EXT_ROUTE_INBOARD;
      top->phi_reference = TOK_EXT_PHI_INBOARD_MIDPLANE;
      return true;
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO:
      top->lower = tok_ext_xray(TOK_EXT_LOWER_XPT, TOK_EXT_SOL_IN);
      top->upper = tok_ext_plate(TOK_EXT_PLATE_LOWER, TOK_EXT_ZMIN);
      top->route = TOK_EXT_ROUTE_GENERIC;
      top->phi_reference = TOK_EXT_PHI_UPPER;
      return true;
    case GKYL_GEOMETRY_TOKAMAK_CORE_R:
      top->lower = tok_ext_xray(TOK_EXT_LOWER_XPT, TOK_EXT_CORE);
      top->upper = tok_ext_xray(TOK_EXT_UPPER_XPT, TOK_EXT_CORE);
      top->route = TOK_EXT_ROUTE_OUTBOARD;
      top->phi_reference = TOK_EXT_PHI_OUTBOARD_MIDPLANE;
      return true;
    case GKYL_GEOMETRY_TOKAMAK_CORE_L:
      top->lower = tok_ext_xray(TOK_EXT_UPPER_XPT, TOK_EXT_CORE);
      top->upper = tok_ext_xray(TOK_EXT_LOWER_XPT, TOK_EXT_CORE);
      top->route = TOK_EXT_ROUTE_INBOARD;
      top->phi_reference = TOK_EXT_PHI_INBOARD_MIDPLANE;
      return true;
    case GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO:
      top->lower = tok_ext_plate(TOK_EXT_PLATE_LOWER, TOK_EXT_ZMIN_RIGHT);
      top->upper = tok_ext_xray(TOK_EXT_LOWER_XPT, TOK_EXT_SOL_OUT);
      top->route = TOK_EXT_ROUTE_GENERIC;
      top->phi_reference = TOK_EXT_PHI_LOWER;
      return true;
    case GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID:
      top->lower = tok_ext_xray(TOK_EXT_LOWER_XPT, TOK_EXT_SOL_OUT);
      top->upper = tok_ext_xray(TOK_EXT_LOWER_XPT, TOK_EXT_SOL_IN);
      top->route = TOK_EXT_ROUTE_VIA_UPPER;
      top->phi_reference = TOK_EXT_PHI_OUTBOARD_MIDPLANE;
      return true;
    case GKYL_GEOMETRY_TOKAMAK_LSN_SOL_UP:
      top->lower = tok_ext_xray(TOK_EXT_LOWER_XPT, TOK_EXT_SOL_IN);
      top->upper = tok_ext_plate(TOK_EXT_PLATE_UPPER, TOK_EXT_ZMIN_LEFT);
      top->route = TOK_EXT_ROUTE_GENERIC;
      top->phi_reference = TOK_EXT_PHI_UPPER;
      return true;
    case GKYL_GEOMETRY_TOKAMAK_CORE:
      top->lower = tok_ext_xray(TOK_EXT_LOWER_XPT, TOK_EXT_CORE);
      top->upper = tok_ext_xray(TOK_EXT_LOWER_XPT, TOK_EXT_CORE);
      top->route = TOK_EXT_ROUTE_CLOSED_CORE;
      top->phi_reference = TOK_EXT_PHI_OUTBOARD_MIDPLANE;
      top->closed = true;
      return true;
    default:
      return false;
  }
}

// Derive seed_names_side from the topology instead of listing block types.
// A block whose BOTH ends are X-point rays shares both endpoints with the
// block covering the other side of the same surface -- CORE_R with CORE_L,
// DN_SOL_OUT_MID with DN_SOL_IN_MID -- so the seed is literally the same point
// for both and cannot say which side this block wanted; only the block can.
// Stating it as a property covers every block of that shape at once, including
// ones added later, rather than each being rediscovered the hard way. Routes
// that are not a Z-monotone branch march never consult it.
//
// LSN_SOL_MID has that shape too, and in the sharpest possible form: BOTH of
// its rays sit on the SAME (lower) X point, so on the separatrix its two
// endpoints are not merely the same kind of point, they are the same point.
// It was excluded here only because its route splits at a turning point and
// chose per branch instead; the exclusion cost TCV its whole separatrix row
// (measured: the b2 trace ran up the inboard flank and back down it, arc
// 1.6407 m inside R<=0.8624, against the core block's 1.8142 m out to
// R=1.0995 -- a 0.42 m seam).  The route clause is gone: the property is about
// the ENDPOINTS, and every route that consults it now gets the same answer.
static bool
tok_ext_topology_from_ftype(
  enum gkyl_tok_geo_type ftype, bool half_domain, struct tok_ext_topology *top
)
{
  if (!tok_ext_topology_from_ftype_raw(ftype, half_domain, top)) {
    return false;
  }
  top->seed_names_side = top->lower.kind == TOK_EXT_XPT_RAY && top->upper.kind == TOK_EXT_XPT_RAY;
  // A block cut from a CLOSED surface by two X-point rays -- both ends are
  // CORE-sector rays -- is half of that surface, and off the separatrix its
  // endpoints do not lie on its own flank at all: the ray lands on the point
  // nearest the X point, which is inboard of the surface's Z extremum.  Derive
  // that from the topology instead of listing CORE_R and CORE_L, so a device
  // that declares such a block gets the right construction without
  // rediscovering this the hard way.  SOL blocks share the two-ray shape but
  // sit on an OPEN surface, where the ray lands on the block's own flank.
  //
  // Only where a Z-monotone HALF march was chosen.  A self-periodic CORE -- the
  // whole closed surface cut at ONE ray, so both its endpoints are the same
  // point, which is ASDEX's core -- is already routed TOK_EXT_ROUTE_CLOSED_CORE
  // and must stay there: it has no second cut point, so there is no sub-arc to
  // extract and the extraction degenerates to the full contour.
  if ((top->route == TOK_EXT_ROUTE_OUTBOARD || top->route == TOK_EXT_ROUTE_INBOARD) &&
      top->seed_names_side && top->lower.sector == TOK_EXT_CORE &&
      top->upper.sector == TOK_EXT_CORE) {
    top->route = TOK_EXT_ROUTE_CORE_HALF;
  }
  return true;
}

// The X points of the representation that evaluates psi, chosen by
// `use_cubics` exactly as geo->psisep is (2026-10-05). The quadratic arrays
// locate a different saddle -- 0.687 mm away on TCV 65402, 0.129 mm on ASDEX
// -- so blocks that pinned or snapped to them met their C1 neighbours that far
// apart at the X point: the radial interfaces SOL_LO|PF_LO_R and
// SOL_UP|PF_LO_L of every TCV cell of the refinement matrix (mx12) failed
// node conformality by exactly that distance, at every resolution.
static int
tok_geo_xpts(const struct gkyl_tok_geo *geo, const double **r, const double **z)
{
  *r = geo->use_cubics ? geo->efit->Rxpt_cubic : geo->efit->Rxpt;
  *z = geo->use_cubics ? geo->efit->Zxpt_cubic : geo->efit->Zxpt;
  return geo->use_cubics ? geo->efit->num_xpts_cubic : geo->efit->num_xpts;
}

static bool tok_ext_xpoint_rz(
  const struct gkyl_tok_geo *geo, enum tok_ext_xpoint which, double *r, double *z
);

// Pin a separatrix node that sits at a block's theta BOUNDARY to the X point
// that bounds that boundary.
//
// This replaces two hand-maintained ftype ladders -- one under
// `if (inp->half_domain)` and one under its `else` -- which between them named
// 36 ftypes to answer a question the topology table already answers.
// `struct tok_ext_topology` records, for each end, whether theta is bounded by
// an X-point ray (and WHICH X point), a divertor plate, or the midplane, and it
// already takes half_domain as an input. So the half_domain split here was
// never a second algorithm; it was the table, transcribed by hand twice.
//
// Equivalence is not assumed. refactor_20260921/tok_xpt_pin_equiv.c enumerates
// every ftype x half_domain x end x {single,double null} -- 168 combinations --
// and the table reproduces the ladders on all 160 that can be declared.
//
// The 8 that differ are all an UPPER block (PF_UP_*, DN_SOL_*_UP) inside a
// LOWER half domain, which is not a geometry: a lower half domain closes on the
// midplane and cannot contain the upper X point. No fixture in the tree
// declares one. Both old implementations nevertheless ANSWERED it -- the ladder
// silently declined to pin, the table silently pins to an X point -- so a
// driver that declared such a block got a quietly wrong grid either way. It now
// aborts, which is the only one of the three behaviours that is honest.
static bool
tok_xpt_sep_pin_z(
  const struct gkyl_tok_geo_grid_inp *inp, const struct gkyl_tok_geo *geo, bool at_upper,
  double *z_pin
)
{
  if (inp->half_domain &&
      (inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_UP_L || inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_UP_R ||
       inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP ||
       inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_UP)) {
    fprintf(
      stderr,
      "TOK_HALF_DOMAIN_UPPER_BLOCK ftype=%d: a lower half domain closes on the "
      "midplane and cannot contain the upper X point\n",
      inp->ftype
    );
    abort();
  }
  struct tok_ext_topology top;
  if (!tok_ext_topology_from_ftype(inp->ftype, inp->half_domain, &top)) {
    return false;
  }
  const struct tok_ext_endpoint *e = at_upper ? &top.upper : &top.lower;
  if (e->kind != TOK_EXT_XPT_RAY) {
    return false;
  }
  // The X point of the representation that evaluates psi, as every ray
  // endpoint takes it (tok_ext_xpoint_rz).
  double r_pin;
  return tok_ext_xpoint_rz(
    geo, e->xpoint == TOK_EXT_UPPER_XPT ? TOK_EXT_UPPER_XPT : TOK_EXT_LOWER_XPT, &r_pin, z_pin
  );
}

// The Z coordinate a topology endpoint denotes.
//
// Replaces the Z half of a second hand-written ftype ladder (in
// tok_half_domain_sep_rz): the table already says whether an end is an X-point
// ray, a plate in a named fixed-Z slot, or the midplane, so the Z follows from
// the endpoint rather than from the block's name.
static bool
tok_ext_endpoint_z(const struct arc_length_ctx *arc_ctx, const struct tok_ext_endpoint *e, double *z)
{
  const struct gkyl_tok_geo *geo = arc_ctx->geo;
  switch (e->kind) {
    case TOK_EXT_XPT_RAY: {
      // Index by WHICH X point the endpoint names. Every half-domain block
      // names the lower one, so this is index 0 there and the behaviour is
      // unchanged; writing it generally keeps the helper honest for callers
      // that are not half-domain.
      int xi = (e->xpoint == TOK_EXT_UPPER_XPT && geo->efit->num_xpts > 1) ? 1 : 0;
      *z = geo->use_cubics ? geo->efit->Zxpt_cubic[xi] : geo->efit->Zxpt[xi];
      return true;
    }
    case TOK_EXT_MIDPLANE:
      *z = geo->zmaxis;
      return true;
    case TOK_EXT_PLATE:
      switch (e->fixed_z_slot) {
        case TOK_EXT_ZMIN:
          *z = arc_ctx->zmin;
          return true;
        case TOK_EXT_ZMAX:
          *z = arc_ctx->zmax;
          return true;
        case TOK_EXT_ZMIN_LEFT:
          *z = arc_ctx->zmin_left;
          return true;
        case TOK_EXT_ZMIN_RIGHT:
          *z = arc_ctx->zmin_right;
          return true;
        case TOK_EXT_ZMAX_LEFT:
          *z = arc_ctx->zmax_left;
          return true;
        case TOK_EXT_ZMAX_RIGHT:
          *z = arc_ctx->zmax_right;
          return true;
        default:
          return false;
      }
    default:
      return false;
  }
}

// Which R side of the surface a block closes its root solve on.
//
// NOT an ordered fallback chain -- that would be a rule fitted to the eight
// rows that happen to exist. An endpoint either asserts a side or it does not:
// a midplane end carries `midplane_outboard`, a SOL_OUT/SOL_IN ray carries it
// in the sector, a plate carries it in a _LEFT/_RIGHT slot; CORE and PF rays
// and a plain ZMIN/ZMAX plate assert nothing. Ask both ends. If they assert
// opposite sides the block is declared inconsistently, which is a defect to
// report rather than to resolve by precedence.
enum tok_ext_side { TOK_EXT_SIDE_NONE = 0, TOK_EXT_SIDE_OUTB, TOK_EXT_SIDE_INB };

static enum tok_ext_side
tok_ext_endpoint_side(const struct tok_ext_endpoint *e)
{
  if (e->kind == TOK_EXT_MIDPLANE) {
    return e->midplane_outboard ? TOK_EXT_SIDE_OUTB : TOK_EXT_SIDE_INB;
  }
  if (e->kind == TOK_EXT_PLATE) {
    if (e->fixed_z_slot == TOK_EXT_ZMIN_RIGHT || e->fixed_z_slot == TOK_EXT_ZMAX_RIGHT) {
      return TOK_EXT_SIDE_OUTB;
    }
    if (e->fixed_z_slot == TOK_EXT_ZMIN_LEFT || e->fixed_z_slot == TOK_EXT_ZMAX_LEFT) {
      return TOK_EXT_SIDE_INB;
    }
    return TOK_EXT_SIDE_NONE;
  }
  if (e->sector == TOK_EXT_SOL_OUT) {
    return TOK_EXT_SIDE_OUTB;
  }
  if (e->sector == TOK_EXT_SOL_IN) {
    return TOK_EXT_SIDE_INB;
  }
  return TOK_EXT_SIDE_NONE;
}

static bool
tok_ext_block_rclose(
  const struct gkyl_tok_geo_grid_inp *inp, const struct tok_ext_topology *top, double *rclose
)
{
  enum tok_ext_side lo = tok_ext_endpoint_side(&top->lower);
  enum tok_ext_side up = tok_ext_endpoint_side(&top->upper);
  if (lo != TOK_EXT_SIDE_NONE && up != TOK_EXT_SIDE_NONE && lo != up) {
    fprintf(
      stderr,
      "TOK_BLOCK_SIDE_CONFLICT ftype=%d: its two theta ends assert opposite "
      "R sides of the surface\n",
      inp->ftype
    );
    abort();
  }
  enum tok_ext_side s = lo != TOK_EXT_SIDE_NONE ? lo : up;
  if (s == TOK_EXT_SIDE_NONE) {
    return false;
  }
  *rclose = s == TOK_EXT_SIDE_OUTB ? inp->rright : inp->rleft;
  return true;
}

static bool
tok_half_domain_sep_rz(
  const struct gkyl_tok_geo_grid_inp *inp, struct arc_length_ctx *arc_ctx, double theta, double *r,
  double *z
)
{
  if (!inp->half_domain || !tok_geo_same_flux(arc_ctx->psi, arc_ctx->geo->psisep)) {
    return false;
  }

  // z0/z1/rclose come from the topology table, which already names this
  // block's two theta ends and the side it sits on. tok_sep_segment_equiv.c
  // checks the table reproduces the ladder this replaces on all 8 half-domain
  // ftypes, for both endpoints and rclose.
  struct tok_ext_topology top;
  if (!tok_ext_topology_from_ftype(inp->ftype, inp->half_domain, &top)) {
    return false;
  }
  double z0 = 0.0, z1 = 0.0, rclose = 0.0;
  if (!tok_ext_endpoint_z(arc_ctx, &top.lower, &z0) ||
      !tok_ext_endpoint_z(arc_ctx, &top.upper, &z1) || !tok_ext_block_rclose(inp, &top, &rclose)) {
    return false;
  }

  // Use the actual (position-mapped) computational coordinate.  Corner
  // nodes lie at frac=0 or 1, while interior and radial/alpha-surface
  // quadrature nodes do not.  Inferring frac from the nodal-array index
  // incorrectly maps the first and last Gauss points onto the segment ends.
  double frac = (theta - inp->cgrid.lower[2]) / (inp->cgrid.upper[2] - inp->cgrid.lower[2]);
  double zfixed = tok_sep_fixed_edge_is_first(inp->ftype) ? z0 : z1;
  tok_build_sep_trace(inp, arc_ctx, zfixed, rclose);
  int n = arc_ctx->sep_trace_n;
  double total = arc_ctx->sep_trace_s[n - 1];
  arc_ctx->xpt_map_darc_dtheta = total / (inp->cgrid.upper[2] - inp->cgrid.lower[2]);
  if (frac <= 0.0) {
    *r = arc_ctx->sep_trace_r[0];
    *z = arc_ctx->sep_trace_z[0];
    return true;
  }
  if (frac >= 1.0) {
    *r = arc_ctx->sep_trace_r[n - 1];
    *z = arc_ctx->sep_trace_z[n - 1];
    return true;
  }
  double target = frac * total;
  int lo = 0, hi = n - 1;
  while (hi - lo > 1) {
    int mid = (lo + hi) / 2;
    if (arc_ctx->sep_trace_s[mid] < target) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  double ds = arc_ctx->sep_trace_s[hi] - arc_ctx->sep_trace_s[lo];
  double w = ds > 0.0 ? (target - arc_ctx->sep_trace_s[lo]) / ds : 0.0;
  double rlin =
    arc_ctx->sep_trace_r[lo] + w * (arc_ctx->sep_trace_r[hi] - arc_ctx->sep_trace_r[lo]);
  double zlin =
    arc_ctx->sep_trace_z[lo] + w * (arc_ctx->sep_trace_z[hi] - arc_ctx->sep_trace_z[lo]);
  if (arc_ctx->sep_trace_param_is_r) {
    double roots[16] = {0.0};
    int nr = tok_geo_Z_psiR(arc_ctx->geo, arc_ctx->geo->psisep, rlin, 16, roots);
    if (nr == 0) {
      fprintf(
        stderr, "TOK_SEP_TRACE no Z root during resampling ftype=%d frac=%.17g R=%.17g\n",
        inp->ftype, frac, rlin
      );
      abort();
    }
    *r = rlin;
    *z = tok_nearest_value(zlin, roots, nr);
  } else {
    double roots[8] = {0.0}, dRdZ[8] = {0.0};
    double dR[8] = {0.0}, dZ[8] = {0.0};
    int nr = gkyl_tok_geo_R_psiZ(arc_ctx->geo, arc_ctx->geo->psisep, zlin, 8, roots, dRdZ, dR, dZ);
    if (nr == 0) {
      fprintf(
        stderr, "TOK_SEP_TRACE no R root during resampling ftype=%d frac=%.17g Z=%.17g\n",
        inp->ftype, frac, zlin
      );
      abort();
    }
    *r = tok_nearest_value(rlin, roots, nr);
    *z = zlin;
  }
  double residual = tok_eval_psi_rz_local(arc_ctx->geo, *r, *z) - arc_ctx->geo->psisep;
  if (!isfinite(*r) || !isfinite(*z) || !isfinite(residual) ||
      fabs(residual) > 1e-9 * fmax(1.0, fabs(arc_ctx->geo->psisep))) {
    fprintf(
      stderr,
      "TOK_SEP_TRACE invalid resampled point ftype=%d frac=%.17g R=%.17g Z=%.17g residual=%.17g\n",
      inp->ftype, frac, *r, *z, residual
    );
    abort();
  }
  return true;
}

// A material theta face follows the supplied limiter arc as psi varies. Its
// corner chord can cross a concave bend even though every contour ends on the
// wall. Validate that native map separately; never call the chord contained.
static bool
tok_divertor_material_cap(
  const struct gkyl_tok_geo_grid_inp *inp, const struct gkyl_tok_geo *geo, int theta_edge,
  const double p[2], const double q[2]
)
{
  struct tok_ext_topology top;
  if (!geo->extend_to_limiter || !geo->plate_spec ||
      !tok_ext_topology_from_ftype(inp->ftype, inp->half_domain, &top)) {
    return false;
  }
  const struct tok_ext_endpoint *end = theta_edge ? &top.upper : &top.lower;
  if (end->kind != TOK_EXT_PLATE) {
    return false;
  }
  int slot = end->plate_slot == TOK_EXT_PLATE_LOWER ? 0 : 1;
  const struct gkyl_tok_geo_wall_target *target = &geo->divertor_wall[slot];
  if (!tok_wall_target_cap(geo->efit, target->num_segments, target->segments, p, q)) {
    return false;
  }
  double sp, sq;
  if (!tok_wall_target_coordinate(geo->efit, target->num_segments, target->segments, p, &sp) ||
      !tok_wall_target_coordinate(geo->efit, target->num_segments, target->segments, q, &sq)) {
    return false;
  }
  double psi_p = tok_eval_psi_rz_local(geo, p[0], p[1]);
  double psi_q = tok_eval_psi_rz_local(geo, q[0], q[1]);
  double tol = 1e-8 * fmax(1.0, fmax(geo->efit->rdim, geo->efit->zdim));
  double previous = sp;
  for (int k = 0; k <= 16; ++k) {
    double psi = psi_p + (psi_q - psi_p) * k / 16.0, rz[2], s;
    if (!tok_divertor_wall_intersection(geo, slot, psi, &rz[0], &rz[1]) ||
        !tok_wall_target_coordinate(geo->efit, target->num_segments, target->segments, rz, &s)) {
      return false;
    }
    if (s < fmin(sp, sq) - tol || s > fmax(sp, sq) + tol ||
        (sq >= sp ? s < previous - tol : s > previous + tol)) {
      return false;
    }
    previous = s;
  }
  // Include every wall vertex traversed by the material map: uniform psi
  // samples alone must not skip a bend or substitute a disconnected root.
  for (int k = 0; k < target->num_segments; ++k) {
    for (int d = 0; d < 2; ++d) {
      int v = (target->segments[k] + d) % geo->efit->limiter_n;
      double vertex[2] = {geo->efit->limiter_R[v], geo->efit->limiter_Z[v]}, s;
      if (!tok_wall_target_coordinate(
            geo->efit, target->num_segments, target->segments, vertex, &s
          )) {
        return false;
      }
      if (s < fmin(sp, sq) + tol || s > fmax(sp, sq) - tol) {
        continue;
      }
      double psi = tok_eval_psi_rz_local(geo, vertex[0], vertex[1]), rz[2];
      double ftol = 1e-10 * fmax(1.0, fmax(fabs(psi_p), fabs(psi_q)));
      if (psi < fmin(psi_p, psi_q) - ftol || psi > fmax(psi_p, psi_q) + ftol ||
          !tok_divertor_wall_intersection(geo, slot, psi, &rz[0], &rz[1]) ||
          hypot(rz[0] - vertex[0], rz[1] - vertex[1]) > tol) {
        return false;
      }
    }
  }
  fprintf(
    stderr,
    "TOK_GEO_DIVERTOR_CAP trial=%d ftype=%d theta_edge=%d slot=%d psi0=%.17g psi1=%.17g R0=%.17g Z0=%.17g R1=%.17g Z1=%.17g s0=%.17g s1=%.17g straight_chord_outside=1 native_arc=validated\n",
    wall_trial_active, inp->ftype, theta_edge, slot, psi_p, psi_q, p[0], p[1], q[0], q[1], sp, sq
  );
  return true;
}

// Does this theta end of the block terminate on a plate the driver declared
// SEPARATELY from the vessel outline? Structural, from the topology table: the
// end is a TOK_EXT_PLATE and its plate is not an explicit outline target. Nodes
// on such a plate sit where the plate function puts them, a fraction of a
// millimetre either side of the outline chord (ASDEX: 0.07-0.94 mm), so the
// wall tests that touch them keep the per-edge slack while every other test
// is judged to roundoff. A plate that IS an outline
// target (extend_to_limiter with divertor_wall segments, NSTX-U) puts its
// nodes on the outline to roundoff and needs no slack.
static bool
tok_wall_theta_end_on_declared_plate(
  const struct gkyl_tok_geo_grid_inp *inp, const struct gkyl_tok_geo *geo, int theta_edge
)
{
  struct tok_ext_topology top;
  if (!tok_ext_topology_from_ftype(inp->ftype, inp->half_domain, &top)) {
    return false;
  }
  const struct tok_ext_endpoint *end = theta_edge ? &top.upper : &top.lower;
  if (end->kind != TOK_EXT_PLATE) {
    return false;
  }
  int slot = end->plate_slot == TOK_EXT_PLATE_LOWER ? 0 : 1;
  return !(geo->extend_to_limiter && geo->divertor_wall[slot].num_segments);
}

// The pockets between this block's declared plates and the outline (see
// tok_wall_pocket_build), set for the wall tests that follow on this thread.
// `nodal` (or null) holds node positions over `nrange` whose first/last theta
// rows were placed ON the plates; they become pocket vertices. A plate that
// leaves the outline by more than the outline resolves is reported, by the
// build that ships (not by every wall trial), and used (user decision
// 2026-10-06). A plate that is not finite is refused: never a movable edge.
static void
tok_wall_block_pockets(
  const struct gkyl_tok_geo_grid_inp *inp, const struct gkyl_tok_geo *geo, const bool on_plate[2],
  const struct gkyl_array *nodal, const struct gkyl_range *nrange, struct tok_wall_pocket pocket[2]
)
{
  enum { PSI_IDX, AL_IDX, TH_IDX };
  for (int end = 0; end < 2; ++end) {
    pocket[end] = (struct tok_wall_pocket){0};
    struct tok_ext_topology top;
    if (!on_plate[end] || !tok_ext_topology_from_ftype(inp->ftype, inp->half_domain, &top)) {
      continue;
    }
    const struct tok_ext_endpoint *ep = end ? &top.upper : &top.lower;
    plate_func plate = ep->plate_slot == TOK_EXT_PLATE_LOWER ? geo->plate_func_lower :
                                                               geo->plate_func_upper;
    int nn = 0;
    double *xy = 0;
    if (nodal) {
      xy = gkyl_malloc(2 * (nrange->upper[PSI_IDX] - nrange->lower[PSI_IDX] + 1) * sizeof(double));
      for (int ip = nrange->lower[PSI_IDX]; ip <= nrange->upper[PSI_IDX]; ++ip) {
        int idx[3] = {
          ip, nrange->lower[AL_IDX], end ? nrange->upper[TH_IDX] : nrange->lower[TH_IDX]
        };
        const double *p = gkyl_array_cfetch(nodal, gkyl_range_idx(nrange, idx));
        xy[2 * nn] = p[0];
        xy[2 * nn + 1] = p[1];
        ++nn;
      }
    }
    struct tok_wall_plate_report rep;
    if (!tok_wall_pocket_build(geo->efit, plate, nn, xy, &pocket[end], &rep)) {
      fprintf(
        stderr,
        "TOK_GEO_WALL_PLATE_INVALID ftype=%d theta_edge=%d rz=(%.17g,%.17g): "
        "the declared plate is not finite\n",
        inp->ftype, end, rep.rz[0], rep.rz[1]
      );
      tok_wall_trial_note_plate_invalid();
      if (!tok_wall_trial_record_where(true, false, false)) {
        abort();
      }
    } else if (rep.beyond_band && !tok_wall_trial_is_active()) {
      fprintf(
        stderr,
        "TOK_GEO_WALL_PLATE_BEYOND_OUTLINE ftype=%d theta_edge=%d rz=(%.17g,%.17g) "
        "outside_m=%.6g beyond_band_m=%.6g report_only=1: the declared plate leaves the vessel "
        "outline by more than the outline resolves there; the legs are judged against the plate\n",
        inp->ftype, end, rep.rz[0], rep.rz[1], rep.outside_m, rep.beyond_m
      );
    }
    if (xy) {
      gkyl_free(xy);
    }
  }
}

// An X-point ray is the shared theta boundary of exactly TWO blocks: every
// (xpoint, sector) pair appears once as some ftype's lower endpoint and once as
// another's upper.  Return the endpoint at the FAR side of the block across
// this ray, so a block can ask what its neighbour needs without any runtime
// communication -- the answer is in the static topology table.
static bool
tok_ext_ray_peer_far_endpoint(
  const struct gkyl_tok_geo_grid_inp *inp, const struct tok_ext_endpoint *ray,
  struct tok_ext_endpoint *peer_far
)
{
  if (ray->kind != TOK_EXT_XPT_RAY) {
    return false;
  }
  for (int f = 0; f <= (int)GKYL_GEOMETRY_TOKAMAK_IWL; ++f) {
    enum gkyl_tok_geo_type ft = (enum gkyl_tok_geo_type)f;
    if (ft == inp->ftype) {
      continue;
    }
    struct tok_ext_topology t;
    if (!tok_ext_topology_from_ftype(ft, inp->half_domain, &t)) {
      continue;
    }
    bool lo = t.lower.kind == TOK_EXT_XPT_RAY && t.lower.xpoint == ray->xpoint &&
              t.lower.sector == ray->sector;
    bool up = t.upper.kind == TOK_EXT_XPT_RAY && t.upper.xpoint == ray->xpoint &&
              t.upper.sector == ray->sector;
    if (lo) {
      *peer_far = t.upper;
      return true;
    }
    if (up) {
      *peer_far = t.lower;
      return true;
    }
  }
  return false;
}

static bool
tok_ext_xpoint_rz(const struct gkyl_tok_geo *geo, enum tok_ext_xpoint which, double *r, double *z)
{
  int nx = geo->use_cubics ? geo->efit->num_xpts_cubic : geo->efit->num_xpts;
  const double *rx = geo->use_cubics ? geo->efit->Rxpt_cubic : geo->efit->Rxpt;
  const double *zx = geo->use_cubics ? geo->efit->Zxpt_cubic : geo->efit->Zxpt;
  if (nx < 1 || (which == TOK_EXT_UPPER_XPT && nx < 2)) {
    return false;
  }
  if (nx == 1) {
    *r = rx[0];
    *z = zx[0];
    return isfinite(*r) && isfinite(*z);
  }
  int ilo = zx[0] <= zx[1] ? 0 : 1;
  int iup = ilo == 0 ? 1 : 0;
  int idx = which == TOK_EXT_LOWER_XPT ? ilo : iup;
  *r = rx[idx];
  *z = zx[idx];
  return isfinite(*r) && isfinite(*z);
}

static bool
tok_ext_nearest_root_at_z(
  const struct gkyl_tok_geo_grid_inp *inp, const struct gkyl_tok_geo *geo,
  enum tok_ext_sector sector, double psi, double z, double rx, double zx, double *rbest,
  double *d2best
)
{
  double R[8] = {0.0}, dRdZ[8] = {0.0};
  double dR[8] = {0.0}, dZ[8] = {0.0};
  int nr = gkyl_tok_geo_R_psiZ(geo, psi, z, 8, R, dRdZ, dR, dZ);
  if (nr <= 0) {
    return false;
  }

  int first = 0, last = nr;
  if (sector == TOK_EXT_SOL_OUT || sector == TOK_EXT_SOL_IN) {
    double ref = sector == TOK_EXT_SOL_OUT ? inp->rright : inp->rleft;
    int selected = 0;
    for (int i = 1; i < nr; ++i) {
      if (fabs(R[i] - ref) < fabs(R[selected] - ref)) {
        selected = i;
      }
    }
    first = selected;
    last = selected + 1;
  }

  bool found = false;
  *d2best = DBL_MAX;
  for (int i = first; i < last; ++i) {
    double d2 = SQ(R[i] - rx) + SQ(z - zx);
    if (d2 < *d2best) {
      *d2best = d2;
      *rbest = R[i];
      found = true;
    }
  }
  return found;
}

static bool
tok_ext_nearest_ray_target(
  const struct gkyl_tok_geo_grid_inp *inp, const struct gkyl_tok_geo *geo,
  const struct tok_ext_endpoint *endpoint, double psi, double *rnear, double *znear
)
{
  double rx = 0.0, zx = 0.0;
  if (!tok_ext_xpoint_rz(geo, endpoint->xpoint, &rx, &zx)) {
    return false;
  }
  double zlo = geo->rzgrid.lower[1], zhi = geo->rzgrid.upper[1];
  int nx = geo->use_cubics ? geo->efit->num_xpts_cubic : geo->efit->num_xpts;
  // A one-X-point SOL has two cuts at the same saddle.  Its unconstrained
  // nearest points can move behind a shaped divertor plate as psi changes,
  // collapsing a leg.  Anchor the two far-surface rays at the outboard and
  // inboard magnetic-midplane roots instead; these cuts preserve the full-SOL
  // topology and are shared exactly by the adjacent blocks.
  if (nx == 1 && (endpoint->sector == TOK_EXT_SOL_OUT || endpoint->sector == TOK_EXT_SOL_IN)) {
    double R[16] = {0.0}, dRdZ[16] = {0.0};
    double dR[16] = {0.0}, dZ[16] = {0.0};
    *znear = geo->zmaxis;
    int nr = gkyl_tok_geo_R_psiZ(geo, psi, *znear, 16, R, dRdZ, dR, dZ);
    if (nr <= 0) {
      return false;
    }
    bool outboard = endpoint->sector == TOK_EXT_SOL_OUT;
    *rnear = tok_nearest_value(outboard ? inp->rright : inp->rleft, R, nr);
    double residual = tok_eval_psi_rz_local(geo, *rnear, *znear) - psi;
    return isfinite(*rnear) && isfinite(*znear) && isfinite(residual) &&
           fabs(residual) <= 1e-9 * fmax(1.0, fabs(psi));
  }
  if (endpoint->sector == TOK_EXT_CORE) {
    zlo = fmin(zx, geo->zmaxis);
    zhi = fmax(zx, geo->zmaxis);
  } else if (endpoint->sector == TOK_EXT_PF) {
    if (endpoint->xpoint == TOK_EXT_LOWER_XPT) {
      zlo = geo->rzgrid.lower[1];
      zhi = zx;
    } else {
      zlo = zx;
      zhi = geo->rzgrid.upper[1];
    }
  } else if (endpoint->xpoint == TOK_EXT_LOWER_XPT) {
    zlo = geo->rzgrid.lower[1];
    zhi = geo->zmaxis;
  } else {
    zlo = geo->zmaxis;
    zhi = geo->rzgrid.upper[1];
  }

  if (!(zhi > zlo)) {
    return false;
  }

  int nsamp = 4 * geo->rzgrid.cells[1] + 1;
  double best_d2 = DBL_MAX, best_r = 0.0, best_z = 0.0;
  int best_i = -1;
  for (int i = 0; i < nsamp; ++i) {
    double z = zlo + (zhi - zlo) * i / (nsamp - 1.0), r = 0.0, d2 = 0.0;
    if (tok_ext_nearest_root_at_z(inp, geo, endpoint->sector, psi, z, rx, zx, &r, &d2) &&
        d2 < best_d2) {
      best_d2 = d2;
      best_r = r;
      best_z = z;
      best_i = i;
    }
  }
  if (best_i < 0) {
    return false;
  }

  double dz = (zhi - zlo) / (nsamp - 1.0);
  double a = fmax(zlo, best_z - dz), b = fmin(zhi, best_z + dz);
  const double gr = 0.6180339887498948482;
  for (int iter = 0; iter < 64; ++iter) {
    double ztrial[2] = {b - gr * (b - a), a + gr * (b - a)};
    double dtrial[2] = {DBL_MAX, DBL_MAX};
    double rtrial[2] = {best_r, best_r};
    for (int k = 0; k < 2; ++k) {
      tok_ext_nearest_root_at_z(
        inp, geo, endpoint->sector, psi, ztrial[k], rx, zx, &rtrial[k], &dtrial[k]
      );
    }
    if (dtrial[0] < best_d2) {
      best_d2 = dtrial[0];
      best_r = rtrial[0];
      best_z = ztrial[0];
    }
    if (dtrial[1] < best_d2) {
      best_d2 = dtrial[1];
      best_r = rtrial[1];
      best_z = ztrial[1];
    }
    if (dtrial[0] <= dtrial[1]) {
      b = ztrial[1];
    } else {
      a = ztrial[0];
    }
  }
  *rnear = best_r;
  *znear = best_z;
  double residual = tok_eval_psi_rz_local(geo, best_r, best_z) - psi;
  return isfinite(best_r) && isfinite(best_z) && isfinite(residual) &&
         fabs(residual) <= 1e-9 * fmax(1.0, fabs(psi));
}

// Intersect the current surface with the fixed line from the selected X point
// to its topology-safe target on the far radial surface.  Counting crossings
// before bisection prevents a locally reversed DG flux representation from
// silently selecting a different branch of the ray.
static bool
tok_ext_fixed_ray_endpoint(
  const struct gkyl_tok_geo_grid_inp *inp, const struct arc_length_ctx *arc_ctx,
  const struct tok_ext_endpoint *endpoint, double psi, double *r, double *z, bool *ambiguous
)
{
  if (ambiguous) {
    *ambiguous = false;
  }
  const struct gkyl_tok_geo *geo = arc_ctx->geo;
  double rx = 0.0, zx = 0.0, rf = 0.0, zf = 0.0;
  if (!tok_ext_xpoint_rz(geo, endpoint->xpoint, &rx, &zx) ||
      !tok_ext_nearest_ray_target(inp, geo, endpoint, arc_ctx->xpt_ray_psi0, &rf, &zf)) {
    return false;
  }

  double scale = fmax(1.0, fmax(fabs(psi), fabs(geo->psisep)));
  if (tok_geo_same_flux(psi, geo->psisep)) {
    *r = rx;
    *z = zx;
    return true;
  }
  if (tok_geo_same_flux(psi, arc_ctx->xpt_ray_psi0)) {
    *r = rf;
    *z = zf;
    return true;
  }
  double delta = arc_ctx->xpt_ray_psi0 - geo->psisep;
  if (!isfinite(delta) || fabs(delta) <= 256.0 * DBL_EPSILON * scale) {
    return false;
  }
  double qtarget = (psi - geo->psisep) / delta;
  if (qtarget <= 0.0 || qtarget >= 1.0) {
    return false;
  }

  const bool use_last = arc_ctx->ext_ray_use_last_crossing;
  const int nsamp = 512;
  const double crossing_tol = 1e-13 * fmax(1.0, fabs(qtarget));
  double qprev = (tok_eval_psi_rz_local(geo, rx, zx) - geo->psisep) / delta;
  double sprev = 0.0, slo = -1.0, shi = -1.0;
  int upward = 0, downward = 0;
  // Lowest q seen before the ray first reaches the target surface.
  double qmin_before_first = 0.0;
  // Crossing census, for telling a spurious level-set island closed off by the
  // C0 grad-psi ridge on a DG cell face from a genuine fold-back of the ray.
  const int max_events = 16;
  double event_s[16];
  bool event_up[16];
  int n_events = 0;
  double qmax = qprev, qend = qprev;
  for (int i = 1; i <= nsamp; ++i) {
    double s = i / (double)nsamp;
    double rs = rx + s * (rf - rx), zs = zx + s * (zf - zx);
    double q = (tok_eval_psi_rz_local(geo, rs, zs) - geo->psisep) / delta;
    if (!isfinite(q)) {
      return false;
    }
    bool crosses_up = qprev < qtarget - crossing_tol && q >= qtarget - crossing_tol;
    bool crosses_down = qprev > qtarget + crossing_tol && q <= qtarget + crossing_tol;
    if (crosses_up) {
      ++upward;
      // q(0)=0 at the X point and q(1)=1 at the far ray target by construction,
      // so the last upward crossing is the point beyond which the ray never
      // re-enters the block.  Taking it is the running-minimum-from-the-right
      // of q, which is monotone and meets qtarget exactly once, so the anchor
      // stays continuous in psi even when a spurious level-set island sits
      // between the X point and the real surface.
      if (use_last || slo < 0.0) {
        slo = sprev;
        shi = s;
      }
    }
    if (crosses_down) {
      ++downward;
    }
    if ((crosses_up || crosses_down) && n_events < max_events) {
      event_s[n_events] = 0.5 * (sprev + s);
      event_up[n_events] = crosses_up;
      ++n_events;
    }
    if (upward == 0) {
      qmin_before_first = fmin(qmin_before_first, q);
    }
    qmax = fmax(qmax, q);
    qend = q;
    qprev = q;
    sprev = s;
  }
  // More than one upward crossing means the ray meets this surface more than
  // once and the anchor is a CHOICE, not a solve.  Report it: the choice has to
  // be made jointly with the block across the ray, which cannot know to ask
  // unless it is told the question exists.
  if (ambiguous) {
    *ambiguous = upward > 1;
  }
  if (upward > 1) {
    fprintf(
      stderr,
      "TOK_EXT_RAY_CROSSINGS ftype=%d sector=%d psi=%.17g qtarget=%.17g "
      "upward=%d downward=%d qmax=%.6e qend=%.6e raylen=%.6e events=",
      inp->ftype, endpoint->sector, psi, qtarget, upward, downward, qmax, qend,
      hypot(rf - rx, zf - zx)
    );
    for (int k = 0; k < n_events; ++k) {
      fprintf(stderr, "%s%c:%.6f", k ? "," : "", event_up[k] ? 'u' : 'd', event_s[k]);
    }
    fprintf(stderr, "\n");
  }
  // A ray that folds back crosses psi=psi_curr more than once.  The crossings
  // are not interchangeable: only the FIRST tends to the X point as
  // psi_curr -> psisep, so it is the one that keeps the anchor continuous in
  // psi.  Accept it when the ray got there without dipping inside the
  // separatrix.  This is the same rule the chord construction's
  // tok_xpt_ray_anchor already applies (v11/v13), including measuring the dip
  // against the flux offset being resolved rather than at 1e-13 -- the X
  // point is a saddle, so a sub-millimetre error in its location puts the
  // first samples on the wrong side of psisep by a vanishing amount.
  double band_tol = fmax(crossing_tol, 0.05 * fabs(qtarget));
  bool last_crossing_ok = use_last && upward >= 1 && slo >= 0.0;
  bool first_crossing_ok = !use_last && upward > 1 && slo >= 0.0 && qmin_before_first >= -band_tol;
  if (first_crossing_ok) {
    fprintf(
      stderr,
      "TOK_EXT_RAY first_crossing_anchor ftype=%d sector=%d psi=%.17g "
      "target=%.17g upward=%d downward=%d qmin_before_first=%.6e\n",
      inp->ftype, endpoint->sector, psi, qtarget, upward, downward, qmin_before_first
    );
  }
  if (!first_crossing_ok && !last_crossing_ok && (upward != 1 || downward != 0 || slo < 0.0)) {
    fprintf(
      stderr,
      "TOK_ORDERED_MAP nonunique fixed-ray intersection ftype=%d sector=%d xpoint=%d psi=%.17g target=%.17g upward=%d downward=%d qmin_before_first=%.6e band_tol=%.6e\n",
      inp->ftype, endpoint->sector, endpoint->xpoint, psi, qtarget, upward, downward,
      qmin_before_first, band_tol
    );
    return false;
  }

  double flo =
    (tok_eval_psi_rz_local(geo, rx + slo * (rf - rx), zx + slo * (zf - zx)) - geo->psisep) / delta -
    qtarget;
  for (int k = 0; k < 70; ++k) {
    double smid = 0.5 * (slo + shi);
    double fm =
      (tok_eval_psi_rz_local(geo, rx + smid * (rf - rx), zx + smid * (zf - zx)) - geo->psisep) /
        delta -
      qtarget;
    if (!isfinite(fm)) {
      return false;
    }
    if (flo * fm <= 0.0) {
      shi = smid;
    } else {
      slo = smid;
      flo = fm;
    }
  }
  double s = 0.5 * (slo + shi);
  *r = rx + s * (rf - rx);
  *z = zx + s * (zf - zx);
  double residual = tok_eval_psi_rz_local(geo, *r, *z) - psi;
  return isfinite(*r) && isfinite(*z) && isfinite(residual) && fabs(residual) <= 1e-9 * scale;
}

static double
tok_ext_fixed_z(const struct gkyl_tok_geo_grid_inp *inp, enum tok_ext_fixed_z_slot slot)
{
  switch (slot) {
    case TOK_EXT_ZMIN:
      return inp->zmin;
    case TOK_EXT_ZMAX:
      return inp->zmax;
    case TOK_EXT_ZMIN_LEFT:
      return inp->zmin_left;
    case TOK_EXT_ZMIN_RIGHT:
      return inp->zmin_right;
    case TOK_EXT_ZMAX_LEFT:
      return inp->zmax_left;
    case TOK_EXT_ZMAX_RIGHT:
      return inp->zmax_right;
  }
  return 0.0;
}

static bool
tok_ext_endpoint_point(
  const struct gkyl_tok_geo_grid_inp *inp, const struct arc_length_ctx *arc_ctx,
  const struct tok_ext_endpoint *endpoint, double psi, bool separatrix, double *r, double *z,
  bool *ambiguous
)
{
  const struct gkyl_tok_geo *geo = arc_ctx->geo;
  if (ambiguous) {
    *ambiguous = false;
  }
  if (endpoint->kind == TOK_EXT_XPT_RAY) {
    if (separatrix) {
      // The separatrix anchor IS the X point; there is nothing to choose.
      return tok_ext_xpoint_rz(geo, endpoint->xpoint, r, z);
    }
    return tok_ext_fixed_ray_endpoint(inp, arc_ctx, endpoint, psi, r, z, ambiguous);
  }
  if (endpoint->kind == TOK_EXT_MIDPLANE) {
    // The midplane is a symmetry boundary, not a material one, so the same
    // construction serves the separatrix and every interior surface: take the
    // root of this psi contour at Z = zmaxis on the requested side.
    *z = geo->zmaxis;
    double R[8] = {0.0}, dRdZ[8] = {0.0};
    double dR[8] = {0.0}, dZ[8] = {0.0};
    int nr = gkyl_tok_geo_R_psiZ(geo, psi, *z, 8, R, dRdZ, dR, dZ);
    if (nr <= 0) {
      fprintf(
        stderr,
        "TOK_EXT_ENDPOINT reason=no_midplane_root ftype=%d psi=%.17g "
        "zmaxis=%.17g outboard=%d rmin=%.17g rleft=%.17g rright=%.17g\n",
        inp->ftype, psi, *z, (int)endpoint->midplane_outboard, geo->rmin, inp->rleft, inp->rright
      );
      return false;
    }
    *r = tok_nearest_value(endpoint->midplane_outboard ? inp->rright : inp->rleft, R, nr);
    if (!isfinite(*r)) {
      fprintf(
        stderr, "TOK_EXT_ENDPOINT reason=nonfinite_midplane_root ftype=%d psi=%.17g nr=%d\n",
        inp->ftype, psi, nr
      );
      return false;
    }
    return true;
  }
  plate_func plate = endpoint->plate_slot == TOK_EXT_PLATE_LOWER ? geo->plate_func_lower :
                                                                   geo->plate_func_upper;
  if (geo->plate_spec && plate) {
    return tok_plate_flux_intersection(geo, plate, psi, inp->ftype, r, z);
  }

  *z = tok_ext_fixed_z(inp, endpoint->fixed_z_slot);
  double R[8] = {0.0}, dRdZ[8] = {0.0};
  double dR[8] = {0.0}, dZ[8] = {0.0};
  int nr = gkyl_tok_geo_R_psiZ(geo, psi, *z, 8, R, dRdZ, dR, dZ);
  if (nr <= 0) {
    return false;
  }
  bool outboard = endpoint->fixed_z_slot == TOK_EXT_ZMIN_RIGHT ||
                  endpoint->fixed_z_slot == TOK_EXT_ZMAX_RIGHT;
  if (endpoint->fixed_z_slot == TOK_EXT_ZMIN || endpoint->fixed_z_slot == TOK_EXT_ZMAX) {
    outboard = inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO ||
               inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP;
  }
  *r = tok_nearest_value(outboard ? inp->rright : inp->rleft, R, nr);
  return isfinite(*r);
}

// resid_tol is the largest flux residual any point on this trace may carry.
// It is 0 for every ordinary trace, meaning the default 2e-9 gate; a builder
// that had to fill an X-point gap raises it to just past the miss it measured
// there, so the relaxation is bounded and confined to that one trace.
static bool
tok_ext_finalize_trace_tol(
  const struct gkyl_tok_geo *geo, double psi, int ftype, int n, double *r, double *z, double *s,
  double resid_tol
)
{
  if (n < 2) {
    return false;
  }
  double tol = fmax(2e-9 * fmax(1.0, fabs(psi)), resid_tol);
  s[0] = 0.0;
  double max_step = 0.0, total = 0.0;
  for (int i = 0; i < n; ++i) {
    double residual = tok_eval_psi_rz_local(geo, r[i], z[i]) - psi;
    if (!isfinite(r[i]) || !isfinite(z[i]) || !isfinite(residual) || fabs(residual) > tol) {
      fprintf(
        stderr,
        "TOK_EXT_TRACE invalid point ftype=%d i=%d psi=%.17g R=%.17g Z=%.17g residual=%.17g\n",
        ftype, i, psi, r[i], z[i], residual
      );
      return false;
    }
    if (i > 0) {
      double ds = hypot(r[i] - r[i - 1], z[i] - z[i - 1]);
      if (!(ds > 0.0) || !isfinite(ds)) {
        return false;
      }
      total += ds;
      max_step = fmax(max_step, ds);
      s[i] = total;
    }
  }
  double mean = total / (n - 1);
  double cell_diag = hypot(geo->rzgrid.dx[0], geo->rzgrid.dx[1]);
  if (!(total > 0.0) || max_step > 4.0 * cell_diag || max_step > 32.0 * mean) {
    int imax = 0;
    double smax = 0.0;
    for (int i = 1; i < n; ++i) {
      double ds = hypot(r[i] - r[i - 1], z[i] - z[i - 1]);
      if (ds > smax) {
        smax = ds;
        imax = i;
      }
    }
    fprintf(
      stderr,
      "TOK_EXT_TRACE discontinuity ftype=%d psi=%.17g n=%d max_step=%.17g mean_step=%.17g cell_diag=%.17g imax=%d\n",
      ftype, psi, n, max_step, mean, cell_diag, imax
    );
    fprintf(stderr, "TOK_EXT_TRACE_STEPS imax=%d/%d neigh=", imax, n - 1);
    for (int k = imax - 4; k <= imax + 4; ++k) {
      if (k < 1 || k > n - 1) {
        continue;
      }
      double ds = hypot(r[k] - r[k - 1], z[k] - z[k - 1]);
      fprintf(stderr, k == imax ? " [%.6g]" : " %.6g", ds);
    }
    fprintf(stderr, "\n");
    // Where the big step sits tells which defect this is: at i=1 with the two
    // points on opposite flanks it is a seed that does not lie on the branch
    // the march follows, not a resolution artefact.
    fprintf(
      stderr,
      "TOK_EXT_TRACE_ENDS ftype=%d psi=%.17g first=(%.17g,%.17g) "
      "last=(%.17g,%.17g) before_imax=(%.17g,%.17g) at_imax=(%.17g,%.17g)\n",
      ftype, psi, r[0], z[0], r[n - 1], z[n - 1], r[imax - 1], z[imax - 1], r[imax], z[imax]
    );
    // How many roots the surface has AT the two Z values straddling the big
    // step separates the two candidate causes: a seed sitting exactly on the
    // turning point (one root, sqrt law, benign) from a seed on the far flank
    // that the march leaves behind (two roots, the march taking the other).
    for (int q = 0; q < 2; ++q) {
      double zq = z[imax - 1 + q];
      double Rq[16] = {0.0}, dRdZq[16] = {0.0};
      double dRq[16] = {0.0}, dZq[16] = {0.0};
      int nrq = gkyl_tok_geo_R_psiZ(geo, psi, zq, 16, Rq, dRdZq, dRq, dZq);
      fprintf(
        stderr,
        "TOK_EXT_TRACE_ROOTS ftype=%d psi=%.17g which=%s "
        "z=%.17g nroots=%d",
        ftype, psi, q ? "at_imax" : "before_imax", zq, nrq
      );
      for (int k = 0; k < nrq; ++k) {
        fprintf(stderr, " %.17g", Rq[k]);
      }
      fprintf(stderr, "\n");
    }
    return false;
  }
  return true;
}

static bool
tok_ext_finalize_trace(
  const struct gkyl_tok_geo *geo, double psi, int ftype, int n, double *r, double *z, double *s
)
{
  return tok_ext_finalize_trace_tol(geo, psi, ftype, n, r, z, s, 0.0);
}

static bool tok_eval_psi_grad_rz_local(
  const struct gkyl_tok_geo *geo, double R, double Z, double *dpsidR, double *dpsidZ
);

// Predictor-corrector contour follower.
//
// An R- or Z-parameterized sweep has a preferred coordinate and degenerates at
// that coordinate's turning points. A half-domain CORE_R boundary carries one
// of each: it rounds the surface's lower Z turning point and then terminates on
// the midplane, which is an R turning point. No single parameterization works,
// and forcing one through with endpoint clustering leaves a folded corner cell.
//
// Following by arc length has no preferred direction: step along the tangent
// (perpendicular to grad psi), then Newton back onto psi=const along grad psi.
// The walk resolves the curve far more finely than the output needs; the result
// is resampled to n points at uniform arc length and each one re-projected, so
// the returned polyline sits on the contour to solver tolerance rather than on
// chords across it.
// Is the X-point saddle well enough conditioned for a contour follower to
// terminate on it reliably?
//
// The follower stops by proximity to its target. Where the saddle is strongly
// anisotropic its four separatrix branches are nearly collinear -- 203585 has
// eigenvalue ratio 0.006 and a 9.1 deg wedge, against 0.83 and 85 deg on a
// healthy one -- and an arriving walk cannot tell which branch it is on, so it
// continues onto the wrong one, loops the private-flux region and returns an
// arc 2.9x too long.
//
// Measured over the CORE_L shots that fold and a clean control set, the
// conditioning min|lambda|/max|lambda| separates them with no overlap:
// folded 0.0048-0.161, clean 0.477-0.952. The threshold sits in that gap.
static bool
tok_ext_xpoint_well_conditioned(const struct gkyl_tok_geo *geo, double rx, double zx)
{
  const double h = 1.0e-4;
  double p0 = tok_eval_psi_rz_local(geo, rx, zx);
  double prr =
    (tok_eval_psi_rz_local(geo, rx + h, zx) - 2.0 * p0 + tok_eval_psi_rz_local(geo, rx - h, zx)) /
    (h * h);
  double pzz =
    (tok_eval_psi_rz_local(geo, rx, zx + h) - 2.0 * p0 + tok_eval_psi_rz_local(geo, rx, zx - h)) /
    (h * h);
  double prz =
    (tok_eval_psi_rz_local(geo, rx + h, zx + h) - tok_eval_psi_rz_local(geo, rx + h, zx - h) -
     tok_eval_psi_rz_local(geo, rx - h, zx + h) + tok_eval_psi_rz_local(geo, rx - h, zx - h)) /
    (4.0 * h * h);
  double tr = prr + pzz, det = prr * pzz - prz * prz;
  double disc = 0.25 * tr * tr - det;
  if (!isfinite(disc) || disc < 0.0) {
    return false;
  }
  double s = sqrt(disc);
  double a = fabs(0.5 * tr + s), b = fabs(0.5 * tr - s);
  double hi = fmax(a, b), lo = fmin(a, b);
  if (!(hi > 0.0) || !isfinite(hi)) {
    return false;
  }
  return lo / hi >= 0.3;
}

static bool
tok_ext_follow_walk(
  const struct gkyl_tok_geo *geo, double psi, double r0, double z0, double r1, double z1,
  int orient, double ds, int max_steps, double ptol, double zcap, double arm_radius, double *pr,
  double *pz, double *ps, int *nout, double *best_d_out, int *why_out, double *turn_out,
  int *turn_i_out
)
{
  pr[0] = r0;
  pz[0] = z0;
  ps[0] = 0.0;
  int np = 1;
  double tr = 0.0, tz = 0.0;
  bool ok = true, arrived = false;
  double best_d = DBL_MAX, gstart = 0.0;
  double max_turn = 0.0;
  int best_i = -1, rising = 0, why = 0, turn_i = -1;
  for (int step = 0; step < max_steps && ok; ++step) {
    double cr = pr[np - 1], cz = pz[np - 1];
    double gr = 0.0, gz = 0.0;
    if (!tok_eval_psi_grad_rz_local(geo, cr, cz, &gr, &gz)) {
      ok = false;
      break;
    }
    double gm = hypot(gr, gz);
    if (!(gm > 0.0) || !isfinite(gm)) {
      ok = false;
      break;
    }
    if (step == 0) {
      gstart = gm;
    }
    // psi is stationary at an X point, so |grad psi| collapses as the walk
    // reaches one. The tangent is then undefined and the corrector below is
    // ill-posed. When the target IS that X point -- which is exactly the CORE
    // endpoint on the separatrix row -- this is arrival.
    if (gm <= 1.0e-7 * gstart) {
      if (hypot(cr - r1, cz - z1) <= 64.0 * ds) {
        arrived = true;
        why = 3;
        break;
      }
      ok = false;
      break;
    }
    // The contour has two tangent orientations at the start and the chord to
    // the target does NOT reliably pick the right one: a boundary that wraps a
    // turning point can leave its start heading away from the target. So the
    // orientation is an argument and the caller walks both.
    double ur = -gz / gm, uz = gr / gm;
    if (step == 0) {
      ur *= orient;
      uz *= orient;
    } else if (ur * tr + uz * tz < 0.0) {
      ur = -ur;
      uz = -uz;
    }
    // Limit the turn per step.
    //
    // On a smooth contour the tangent rotates by kappa*ds per step: a healthy
    // CORE_R walk peaks at 1.6 deg. Anything far beyond that is not curvature,
    // it is the DG gradient discontinuity at a cell face -- psi is continuous
    // across a face but grad psi is not, and the X-point finder deliberately
    // CLAMPS its answer onto a face (efit_utils.c ~237), so the walk that has
    // to terminate there also has to cross that discontinuity. On 203585's
    // CORE_L walk the tangent snaps 42.3 deg in one 0.55 mm step at exactly
    // that crossing and the walk leaves its branch, ending up 2.9x too long.
    //
    // The kink is in the representation, not the geometry: evaluated 2 cm off
    // the face the same saddle is perfectly ordinary. So cap the rotation and
    // let the corrector below pull the point back onto psi=const -- the walk
    // stays on the contour instead of following the artifact.
    if (step > 0) {
      double dot = ur * tr + uz * tz;
      dot = fmin(1.0, fmax(-1.0, dot));
      double turn = acos(dot);
      if (turn > max_turn) {
        max_turn = turn;
        turn_i = np - 1;
      }
    }
    tr = ur;
    tz = uz;
    double nr = cr + ds * ur, nz = cz + ds * uz;
    for (int k = 0; k < 8; ++k) {
      double f = tok_eval_psi_rz_local(geo, nr, nz) - psi;
      double ggr = 0.0, ggz = 0.0;
      if (!isfinite(f) || !tok_eval_psi_grad_rz_local(geo, nr, nz, &ggr, &ggz)) {
        ok = false;
        break;
      }
      double g2 = ggr * ggr + ggz * ggz;
      if (!(g2 > 0.0)) {
        ok = false;
        break;
      }
      // Clamp the Newton step. The correction is f/|grad psi|, and |grad psi|
      // collapses near an X point, so an unclamped step can throw the point an
      // arbitrary distance and the walk simply resumes on whatever branch it
      // lands on. Measured on 203585's CORE_L rows nearest the separatrix,
      // that returned an arc 2.9x too long -- 3.291 against 1.140 -- which
      // misplaces those rows and folds the block.
      double cdr = -f / g2 * ggr, cdz = -f / g2 * ggz;
      double cl = hypot(cdr, cdz);
      if (cl > ds) {
        cdr *= ds / cl;
        cdz *= ds / cl;
      }
      nr += cdr;
      nz += cdz;
      if (fabs(f) <= ptol) {
        break;
      }
    }
    if (!ok || !isfinite(nr) || !isfinite(nz)) {
      ok = false;
      break;
    }
    // Both orientations eventually reach the target on a closed surface, and
    // arc length alone does not always tell them apart. The boundary of a
    // half-domain block cannot cross the midplane it terminates on, so an arc
    // that climbs above the higher endpoint is the wrong way round.
    if (nz > zcap) {
      ok = false;
      break;
    }
    pr[np] = nr;
    pz[np] = nz;
    ps[np] = ps[np - 1] + hypot(nr - pr[np - 1], nz - pz[np - 1]);
    ++np;
    double dist = hypot(nr - r1, nz - z1);
    if (dist <= ds) {
      arrived = true;
      why = 1;
      break;
    }
    // The target can be an X point, where psi=const has a sharp corner because
    // four separatrix branches meet. A fixed-step walk rounds that corner and
    // can pass at more than ds, miss the test above, and carry on down a
    // divertor leg -- measured on 204051, the separatrix row came back 1.8x too
    // long (total arc 3.111 against 1.728 for its neighbour), which is enough
    // to misplace that whole row and fold the last radial cell. So also stop at
    // a closest approach: once the walk has been near the target and has moved
    // away for a few consecutive steps, take the nearest point it reached.
    if (dist < best_d) {
      best_d = dist;
      best_i = np - 1;
      rising = 0;
    } else if (best_d <= arm_radius && ++rising >= 8) {
      arrived = true;
      why = 2;
      np = best_i + 1;
      break;
    }
  }
  if (best_d_out) {
    *best_d_out = best_d;
  }
  if (why_out) {
    *why_out = why;
  }
  if (turn_out) {
    *turn_out = max_turn;
  }
  if (turn_i_out) {
    *turn_i_out = turn_i;
  }
  if (!ok || !arrived) {
    return false;
  }
  pr[np] = r1;
  pz[np] = z1;
  ps[np] = ps[np - 1] + hypot(r1 - pr[np - 1], z1 - pz[np - 1]);
  ++np;
  *nout = np;
  return true;
}

static bool
tok_ext_follow_contour(
  const struct gkyl_tok_geo *geo, double psi, double r0, double z0, double r1, double z1, int n,
  double *r, double *z
)
{
  if (n < 2) {
    return false;
  }
  double chord = hypot(r1 - r0, z1 - z0);
  if (!(chord > 0.0) || !isfinite(chord)) {
    return false;
  }
  const int max_steps = 64 * n;
  const double ds = chord / (8.0 * n);
  const double ptol = 1e-12 * fmax(1.0, fabs(psi));
  double *pr = gkyl_malloc(sizeof(double[2 * (max_steps + 2)]));
  double *pz = gkyl_malloc(sizeof(double[2 * (max_steps + 2)]));
  double *ps = gkyl_malloc(sizeof(double[2 * (max_steps + 2)]));
  int cap = max_steps + 2, np = 0, npb = 0;
  const double zcap = fmax(z0, z1) + 8.0 * ds;
  // Start from whichever endpoint has the better-conditioned contour direction.
  // An X point is a saddle: grad psi vanishes on it, so the tangent there is
  // pure roundoff and the walk sets off in an arbitrary direction. On 204051's
  // CORE_R separatrix row -- whose lower endpoint IS the X point -- that sent
  // it down a divertor leg and back, returning an arc 1.8x too long (3.111
  // against 1.728 for the surface 3.6 mm away). It does not fail, it succeeds
  // wrongly, so testing for failure is not enough: pick the good end up front.
  double g0r = 0.0, g0z = 0.0, g1r = 0.0, g1z = 0.0;
  double gm0 = tok_eval_psi_grad_rz_local(geo, r0, z0, &g0r, &g0z) ? hypot(g0r, g0z) : 0.0;
  double gm1 = tok_eval_psi_grad_rz_local(geo, r1, z1, &g1r, &g1z) ? hypot(g1r, g1z) : 0.0;
  bool ok = false, reversed = false;
  double dbg_bd = 0.0, dbg_arm = 0.0, dbg_gmt = 0.0, dbg_turn = 0.0;
  int dbg_why = 0, dbg_ti = -1;
  bool dbg_tc = false;
  for (int attempt = 0; attempt < 2 && !ok; ++attempt) {
    bool from_far = (attempt == 0) == (gm1 > gm0);
    double sr = from_far ? r1 : r0, sz = from_far ? z1 : z0;
    double er = from_far ? r0 : r1, ez = from_far ? z0 : z1;
    // The target need not lie on THIS contour. tok_ext_fixed_ray_endpoint
    // snaps the ray endpoint to the X point for every psi that
    // tok_geo_same_flux calls equal to psisep, which covers the separatrix row
    // AND the two finite-difference stencil rows on either side of it. Their
    // contours round the saddle at a finite distance -- 5 cm on 203585's
    // CORE_L -- so an arrival test scaled to the step size can never fire, and
    // the walk sails past, loops the whole private-flux region and returns to
    // the X point along the opposite branch: arc 3.291 against 1.140.
    // When the target is a critical point of psi, accept the closest approach
    // instead, on a radius scaled to the boundary rather than to the step.
    double gtr = 0.0, gtz = 0.0;
    double gm_t = tok_eval_psi_grad_rz_local(geo, er, ez, &gtr, &gtz) ? hypot(gtr, gtz) : 0.0;
    bool target_critical = gm_t <= 1.0e-6 * fmax(gm0, gm1);
    double arm = target_critical ? 0.5 * chord : 8.0 * ds;
    double bd_a = 0.0, bd_b = 0.0, turn_a = 0.0, turn_b = 0.0;
    int why_a = 0, why_b = 0, ti_a = -1, ti_b = -1;
    bool got_a = tok_ext_follow_walk(
      geo, psi, sr, sz, er, ez, +1, ds, max_steps, ptol, zcap, arm, pr, pz, ps, &np, &bd_a, &why_a,
      &turn_a, &ti_a
    );
    bool got_b = tok_ext_follow_walk(
      geo, psi, sr, sz, er, ez, -1, ds, max_steps, ptol, zcap, arm, pr + cap, pz + cap, ps + cap,
      &npb, &bd_b, &why_b, &turn_b, &ti_b
    );
    dbg_bd = got_a ? bd_a : bd_b;
    dbg_why = got_a ? why_a : why_b;
    dbg_turn = got_a ? turn_a : turn_b;
    dbg_ti = got_a ? ti_a : ti_b;
    dbg_arm = arm;
    dbg_tc = target_critical;
    dbg_gmt = gm_t;
    ok = got_a || got_b;
    // Of the arcs that stay under the cap, take the shorter.
    if (got_b && (!got_a || ps[cap + npb - 1] < ps[np - 1])) {
      for (int i = 0; i < npb; ++i) {
        pr[i] = pr[cap + i];
        pz[i] = pz[cap + i];
        ps[i] = ps[cap + i];
      }
      np = npb;
    }
    reversed = ok && from_far;
  }
  if (ok && reversed) {
    double total = ps[np - 1];
    for (int i = 0, j = np - 1; i < j; ++i, --j) {
      double tr = pr[i], tz = pz[i];
      pr[i] = pr[j];
      pz[i] = pz[j];
      pr[j] = tr;
      pz[j] = tz;
      double t = ps[i];
      ps[i] = ps[j];
      ps[j] = t;
    }
    for (int i = 0; i < np; ++i) {
      ps[i] = total - ps[i];
    }
  }

  if (ok) {
    double total = ps[np - 1];
    ok = isfinite(total) && total > 0.0 && np >= 2;
    if (ok) {
      r[0] = r0;
      z[0] = z0;
      int j = 0;
      for (int i = 1; i < n - 1; ++i) {
        double target = total * i / (double)(n - 1);
        while (j < np - 2 && ps[j + 1] < target) {
          ++j;
        }
        double seg = ps[j + 1] - ps[j];
        double w = seg > 0.0 ? (target - ps[j]) / seg : 0.0;
        double rr = pr[j] + w * (pr[j + 1] - pr[j]);
        double zz = pz[j] + w * (pz[j + 1] - pz[j]);
        // Chords cut the corner by O(ds^2); project back so the trace meets
        // the flux residual that tok_ext_finalize_trace enforces.
        for (int k = 0; k < 6; ++k) {
          double f = tok_eval_psi_rz_local(geo, rr, zz) - psi;
          double ggr = 0.0, ggz = 0.0;
          if (!isfinite(f) || !tok_eval_psi_grad_rz_local(geo, rr, zz, &ggr, &ggz)) {
            ok = false;
            break;
          }
          double g2 = ggr * ggr + ggz * ggz;
          if (!(g2 > 0.0)) {
            ok = false;
            break;
          }
          rr -= f / g2 * ggr;
          zz -= f / g2 * ggz;
          if (fabs(f) <= ptol) {
            break;
          }
        }
        if (!ok) {
          break;
        }
        r[i] = rr;
        z[i] = zz;
      }
      r[n - 1] = r1;
      z[n - 1] = z1;
    }
  }
  gkyl_free(pr);
  gkyl_free(pz);
  gkyl_free(ps);
  return ok;
}

static bool
tok_ext_build_open_trace(
  const struct gkyl_tok_geo_grid_inp *inp, const struct gkyl_tok_geo *geo, double psi, int n,
  double r0, double z0, double r1, double z1, double *r, double *z, double *s, bool *param_is_r
)
{
  double *cr = gkyl_malloc(sizeof(double[4 * n]));
  double *cz = gkyl_malloc(sizeof(double[4 * n]));
  // Which construction wins must not change from one flux surface to the next.
  // Downstream sampling is by NORMALIZED arc length, so two constructions that
  // trace the same contour but resolve it differently (one cutting a chord
  // across the midplane turning point, the other not) put the same logical u at
  // different physical points. Alternating between them across psi is what
  // makes a CORE_R theta node jump 13.7 mm between adjacent radial indices --
  // an 80x radial-spacing jump at fixed theta -- and folds the corner cells.
  // Following the contour by arc length is construction-independent, so using
  // it at EVERY psi removes the flip: measured on the 20 CORE_R shots, 15 go
  // from a 1-3 cell fold to clean, and their worst-cell ratio rises from
  // ~0.02-0.2 to ~0.5.
  //
  // Restrict it to boundaries that actually terminate on a turning point, i.e.
  // that have a midplane endpoint. Applying it to the leg blocks as well was
  // measured and is WORSE -- it put new folds into PF_LO_L and DN_SOL_IN_LO on
  // 5 shots.
  //
  // This comment previously added "because their plate endpoints are ordinary
  // points where the scored candidates are already both consistent". THAT PART
  // IS FALSE, measured 2026-09-21: the leg blocks' theta grading SHEARS with
  // psi -- the same theta cell has a different arc fraction on each flux
  // surface -- by 3.2x to 17.1x, and the shear GROWS under theta refinement,
  // while every block that takes the follower or the CORE_HALF route is
  // psi-consistent to <= 0.004. The legs are the LEAST consistent blocks in the
  // tree, not blocks that did not need the fix.
  //
  // Enabling the follower on them anyway does NOT help: a controlled A/B
  // (the follower forced on the leg blocks, stepc/asdexc/tcvc at theta x2/x4/x8) moved
  // 530 of 830 written arrays and left psi-shear unchanged to three decimals --
  // stepc 3.004 -> 3.033, asdexc 4.349 -> 4.346, tcvc 5.863 -> 5.939 -- and
  // still growing. So the shear is common to BOTH constructions and is not the
  // construction-choice flip this function guards against.
  //
  // The 5-shot fold result above also did not reproduce in that A/B (0 of 9
  // cases fail the grid gate on either arm), but it predates the separatrix
  // kink fix, the Jacobian sign guard and the plate-root work, and was measured
  // on NSTX-U shot geometry rather than these fixtures, so it is left standing
  // rather than retired on weaker evidence than it was made with.
  //
  // The remaining structural difference is that tok_ext_build_core_half_trace
  // anchors its polyline on the contour's OWN TURNING POINTS off the separatrix
  // and appends the declared endpoints as caps, while this route marches
  // endpoint to endpoint.
  // Enabled where it is measured to help. With the midplane as the UPPER
  // endpoint (half-domain CORE_R) it took 16 folded shots to 0 and lifted the
  // worst-cell ratio to ~0.5. With the midplane as the LOWER endpoint
  // (CORE_L) it does the opposite -- 0 folded shots to 9, at minA/medA ~ -1e4
  // -- and three attempts at the cause did not move that number:
  //   * clamping the Newton correction near the saddle,
  //   * a critical-point arrival test on |grad psi|,
  //   * a closest-approach arrival radius scaled to the boundary.
  // None changed the fold count at all, so the mechanism is still unidentified
  // and CORE_L stays on the scored-candidate path, where it has no folds.
  // A dump of the follower's walk (a diagnostic since removed) showed the CORE_L
  // walk sailing 5 cm past the X point, looping the entire private-flux region
  // and returning along the opposite branch: arc 3.291 against 1.140.
  struct tok_ext_topology ftop;
  const bool have_ftop = tok_ext_topology_from_ftype(inp->ftype, inp->half_domain, &ftop);
  bool follow_first = have_ftop &&
                      (ftop.upper.kind == TOK_EXT_MIDPLANE || ftop.lower.kind == TOK_EXT_MIDPLANE);
  // Both halves of the core use the follower; the difference is only in when
  // it has to stand down, and that is a measured property, not a choice.
  //
  // The follower terminates by proximity to the X point. Approaching it from
  // the OUTBOARD midplane (CORE_R) that is robust even on a badly conditioned
  // saddle -- 204515 at cond=0.091 is clean, and forcing CORE_R off the
  // follower there folds it. Approaching from the INBOARD midplane (CORE_L)
  // it is not: every CORE_L shot that folds has cond in 0.0048-0.161 while
  // every clean one is 0.477-0.952, no overlap. So CORE_L, and only CORE_L,
  // falls back to the scored candidates on an ill-conditioned saddle.
  if (follow_first && ftop.lower.kind == TOK_EXT_MIDPLANE) {
    double cxr = 0.0, cxz = 0.0;
    if (tok_ext_xpoint_rz(geo, TOK_EXT_LOWER_XPT, &cxr, &cxz) &&
        !tok_ext_xpoint_well_conditioned(geo, cxr, cxz)) {
      follow_first = false;
    }
  }
  if (follow_first && tok_ext_follow_contour(geo, psi, r0, z0, r1, z1, n, cr, cz)) {
    for (int i = 0; i < n; ++i) {
      r[i] = cr[i];
      z[i] = cz[i];
    }
    *param_is_r = false;
    gkyl_free(cr);
    gkyl_free(cz);
    bool fok = tok_ext_finalize_trace(geo, psi, inp->ftype, n, r, z, s);
    return fok;
  }
  double score[4] = {DBL_MAX, DBL_MAX, DBL_MAX, DBL_MAX};
  bool ok[4] = {false, false, false, false};
  ok[0] = tok_build_contour_candidate(
    geo, psi, true, r0, z0, r1, z1, n, cr, cz, &score[0], "extended", false
  );
  ok[1] = tok_build_contour_candidate(
    geo, psi, false, r0, z0, r1, z1, n, cr + n, cz + n, &score[1], "extended", false
  );
  ok[2] = tok_build_contour_candidate(
    geo, psi, true, r1, z1, r0, z0, n, cr + 2 * n, cz + 2 * n, &score[2], "extended_reverse", false
  );
  ok[3] = tok_build_contour_candidate(
    geo, psi, false, r1, z1, r0, z0, n, cr + 3 * n, cz + 3 * n, &score[3], "extended_reverse", false
  );
  int best = -1;
  for (int k = 0; k < 4; ++k) {
    if (ok[k] && (best < 0 || score[k] < score[best])) {
      best = k;
    }
  }
  // Half-domain MID and CORE blocks terminate on the midplane, which is an R
  // turning point of the contour: Z-Z_end scales as sqrt(|R-R_end|) there, so
  // uniform R sampling makes the step touching that endpoint artificially large
  // and trips the discontinuity test on a perfectly smooth contour. Measured on
  // 202945's CORE_R far trace, the plain R candidate misses by 6e-5 relative --
  // max_step 0.0764551 against a 0.0764502 floor, ratio 16.12 against 16 -- and
  // that max_step IS the endpoint step.
  //
  // Both fallbacks below are gated on ALL FOUR plain candidates having failed,
  // so no trace that is accepted today can change. Try the arc-length follower
  // first: it is parameterization-free and therefore handles the R turning
  // point at the midplane AND the Z turning point the CORE_R path rounds on the
  // way, whereas endpoint clustering only relieves the endpoint and still left
  // a folded corner cell on 14 of the 20 shots that needed it.
  if (best < 0 && tok_ext_follow_contour(geo, psi, r0, z0, r1, z1, n, cr, cz)) {
    for (int i = 0; i < n; ++i) {
      r[i] = cr[i];
      z[i] = cz[i];
    }
    *param_is_r = false;
    gkyl_free(cr);
    gkyl_free(cz);
    return tok_ext_finalize_trace(geo, psi, inp->ftype, n, r, z, s);
  }
  if (best < 0) {
    ok[0] = tok_build_contour_candidate(
      geo, psi, true, r0, z0, r1, z1, n, cr, cz, &score[0], "extended_clustered", true
    );
    ok[1] = tok_build_contour_candidate(
      geo, psi, false, r0, z0, r1, z1, n, cr + n, cz + n, &score[1], "extended_clustered", true
    );
    ok[2] = tok_build_contour_candidate(
      geo, psi, true, r1, z1, r0, z0, n, cr + 2 * n, cz + 2 * n, &score[2],
      "extended_clustered_reverse", true
    );
    ok[3] = tok_build_contour_candidate(
      geo, psi, false, r1, z1, r0, z0, n, cr + 3 * n, cz + 3 * n, &score[3],
      "extended_clustered_reverse", true
    );
    for (int k = 0; k < 4; ++k) {
      if (ok[k] && (best < 0 || score[k] < score[best])) {
        best = k;
      }
    }
  }
  if (best < 0) {
    // All four parameterizations were rejected. This used to return silently,
    // which makes the caller's "domain_trace_failed" the only evidence and says
    // nothing about which parameterization came closest.
    fprintf(
      stderr,
      "TOK_EXT_OPEN_TRACE reason=no_viable_candidate ftype=%d psi=%.17g n=%d "
      "endpoints=(%.17g,%.17g)->(%.17g,%.17g) "
      "ok=[R:%d Z:%d Rrev:%d Zrev:%d] score=[%.6g %.6g %.6g %.6g]\n",
      inp->ftype, psi, n, r0, z0, r1, z1, (int)ok[0], (int)ok[1], (int)ok[2], (int)ok[3], score[0],
      score[1], score[2], score[3]
    );
    gkyl_free(cr);
    gkyl_free(cz);
    return false;
  }
  bool reverse = best >= 2;
  int off = best * n;
  for (int i = 0; i < n; ++i) {
    int src = reverse ? n - 1 - i : i;
    r[i] = cr[off + src];
    z[i] = cz[off + src];
  }
  *param_is_r = best == 0 || best == 2;
  gkyl_free(cr);
  gkyl_free(cz);
  bool cok = tok_ext_finalize_trace(geo, psi, inp->ftype, n, r, z, s);
  return cok;
}

// Point of closest approach in psi along a Z=const line, within a window
// around a reference R.  Used only where psi=const has no crossing at all:
// the X-point finder clamps its answer onto a DG cell face, so psisep is off
// the discrete field's saddle value by ~1e-6 and the level set near the X
// point degenerates from two crossing branches into a hyperbola pair.  Inside
// that gap a horizontal line meets neither branch even though the separatrix
// passes within a fraction of a cell.
static bool
tok_psi_line_closest_r(
  const struct gkyl_tok_geo *geo, double psi, double z, double rref, double window, double *r_out,
  double *resid_out
)
{
  double lo = fmax(geo->rmin, rref - window);
  double hi = fmin(geo->rmax, rref + window);
  if (!(hi > lo)) {
    return false;
  }
  const int ns = 256;
  double rbest = lo, fbest = DBL_MAX;
  for (int k = 0; k <= ns; ++k) {
    double rr = lo + (hi - lo) * k / (double)ns;
    double f = fabs(tok_eval_psi_rz_local(geo, rr, z) - psi);
    if (isfinite(f) && f < fbest) {
      fbest = f;
      rbest = rr;
    }
  }
  if (fbest == DBL_MAX) {
    return false;
  }
  double h = (hi - lo) / ns;
  double a = fmax(lo, rbest - h), b = fmin(hi, rbest + h);
  const double gs = 0.6180339887498949;
  double c = b - gs * (b - a), d = a + gs * (b - a);
  double fc = fabs(tok_eval_psi_rz_local(geo, c, z) - psi);
  double fd = fabs(tok_eval_psi_rz_local(geo, d, z) - psi);
  for (int k = 0; k < 80 && (b - a) > 1e-15; ++k) {
    if (fc < fd) {
      b = d;
      d = c;
      fd = fc;
      c = b - gs * (b - a);
      fc = fabs(tok_eval_psi_rz_local(geo, c, z) - psi);
    } else {
      a = c;
      c = d;
      fc = fd;
      d = a + gs * (b - a);
      fd = fabs(tok_eval_psi_rz_local(geo, d, z) - psi);
    }
  }
  double rm = 0.5 * (a + b);
  double fm = fabs(tok_eval_psi_rz_local(geo, rm, z) - psi);
  if (!isfinite(rm) || !isfinite(fm) || fm > fbest) {
    rm = rbest;
    fm = fbest;
  }
  *r_out = rm;
  *resid_out = fm;
  return true;
}

static bool
tok_ext_z_branch_points(
  const struct gkyl_tok_geo_grid_inp *inp, const struct gkyl_tok_geo *geo, double psi, int n,
  double r0, double z0, double r1, double z1, bool outboard, bool seed_names_side, double *r,
  double *z, double *gap_resid
)
{
  if (gap_resid) {
    *gap_resid = 0.0;
  }
  if (n < 2 || fabs(z1 - z0) <= 64.0 * DBL_EPSILON * fmax(1.0, fmax(fabs(z0), fabs(z1)))) {
    return false;
  }
  r[0] = r0;
  z[0] = z0;
  for (int i = 1; i < n - 1; ++i) {
    double f = i / (double)(n - 1);
    z[i] = z0 + f * (z1 - z0);
    double R[8] = {0.0}, dRdZ[8] = {0.0};
    double dR[8] = {0.0}, dZ[8] = {0.0};
    int nr = gkyl_tok_geo_R_psiZ(geo, psi, z[i], 8, R, dRdZ, dR, dZ);
    if (nr <= 0) {
      // No crossing at this Z.  Near an X point that is an artifact of psisep
      // rather than an absent contour, so fall back to the closest approach
      // in psi, continuing from the previous station's branch.  The gate is
      // four orders wide: inside the X-point gap the miss is ~1e-6, while a
      // genuinely absent contour misses by the psi scale of the domain.
      double rc = 0.0, resid = 0.0;
      double gap_tol = 1e-5 * fmax(1.0, fabs(psi));
      if (tok_psi_line_closest_r(geo, psi, z[i], r[i - 1], 8.0 * geo->rzgrid.dx[0], &rc, &resid) &&
          resid <= gap_tol) {
        r[i] = rc;
        if (gap_resid) {
          *gap_resid = fmax(*gap_resid, resid);
        }
        continue;
      }
      fprintf(
        stderr,
        "TOK_EXT_ZBRANCH reason=no_root ftype=%d psi=%.17g i=%d n=%d "
        "z=%.17g z0=%.17g z1=%.17g frac=%.6g outboard=%d "
        "closest_resid=%.17g gap_tol=%.17g\n",
        inp->ftype, psi, i, n, z[i], z0, z1, f, (int)outboard, resid, gap_tol
      );
      return false;
    }
    // Follow the branch by continuity from the previous station rather than
    // by a fixed side hint.  The hint cannot tell two same-side roots apart,
    // and the inboard psi representation carries a spurious root hugging
    // geo->rmin: on 204993's inboard separatrix it captured twelve
    // consecutive stations at R=0.2501 while the contour was at R=0.356,
    // detouring the trace to the inner wall and back.  Station 0 is an exact
    // endpoint, so continuity has a trustworthy seed.
    // ...but that seed is only trustworthy when it names one side of the
    // contour.  Two cases where it cannot.  A closed surface's two branches
    // both start at a Z turning point, the one place the inboard and outboard
    // roots coalesce.  And a closed surface's two HALVES -- CORE_R and CORE_L
    // -- share both X-point-ray endpoints, so the seed is literally the same
    // point for the block that goes outboard and the block that goes inboard.
    // In both, continuity cannot tell the sides apart: both blocks follow the
    // same one and each covers half the surface, and only the block itself
    // knows which half it wanted.  Where the caller says so, break the tie
    // once at the first interior station using the side it asked for --
    // taking the NEAREST root on that side, so a spurious same-side root
    // further out is still rejected.  If no root lies on that side the
    // contour leaves the seed in one direction only, and plain continuity
    // applies.
    //
    // The filter is relative to r[0], so it is only meaningful while r[0] is an
    // EXTREME point of the branch -- a turning point, where the flanks meet.
    // Extending it to every station was tried and reverted: on a branch seeded
    // at an ordinary endpoint (via-turning's branch A) the far end of the
    // contour legitimately comes back past the seed's R, the filter then
    // rejects every root there, and ASDEX's ftype=12 separatrix loses its
    // trace.  It also did not fix what it was written for -- STEP's CORE_R was
    // a wrong ROUTE, not a wrong tie-break; see tok_ext_build_core_half_trace.
    if (seed_names_side && i == 1) {
      double best = 0.0;
      bool found = false;
      for (int k = 0; k < nr; ++k) {
        if (outboard ? (R[k] < r[0]) : (R[k] > r[0])) {
          continue;
        }
        if (!found || fabs(R[k] - r[0]) < fabs(best - r[0])) {
          best = R[k];
          found = true;
        }
      }
      if (found) {
        r[i] = best;
        continue;
      }
    }
    // Near an X point the branches crowd together, and the nearest root to
    // the seed stops being the one this block wants: at the X point itself
    // continuity has nothing to continue from, and just inside it the wrong
    // branch can be nearer by a hair.  Detect that ambiguity directly -- two
    // roots straddling the seed at comparable distance -- and fall back to the
    // block's own radial reference, which is what named the branch before this
    // selection became a continuity march.  Confined to the first station, so
    // the march still rejects the spurious rmin root everywhere after it.
    if (i == 1 && !seed_names_side && nr > 1) {
      double d1 = DBL_MAX, d2 = DBL_MAX, r1v = 0.0, r2v = 0.0;
      for (int k = 0; k < nr; ++k) {
        double d = fabs(R[k] - r[0]);
        if (d < d1) {
          d2 = d1;
          r2v = r1v;
          d1 = d;
          r1v = R[k];
        } else if (d < d2) {
          d2 = d;
          r2v = R[k];
        }
      }
      if ((r1v - r[0]) * (r2v - r[0]) < 0.0 && d2 < 4.0 * d1) {
        r[i] = tok_nearest_value(outboard ? inp->rright : inp->rleft, R, nr);
        continue;
      }
    }
    r[i] = tok_nearest_value(r[i - 1], R, nr);
  }
  r[n - 1] = r1;
  z[n - 1] = z1;
  return true;
}

static bool
tok_ext_build_z_branch_trace(
  const struct gkyl_tok_geo_grid_inp *inp, const struct gkyl_tok_geo *geo, double psi, int n,
  double r0, double z0, double r1, double z1, bool outboard, bool seed_names_side, double *r,
  double *z, double *s, bool *param_is_r
)
{
  double gap_resid = 0.0;
  if (!tok_ext_z_branch_points(
        inp, geo, psi, n, r0, z0, r1, z1, outboard, seed_names_side, r, z, &gap_resid
      )) {
    return false;
  }
  *param_is_r = false;
  return tok_ext_finalize_trace_tol(geo, psi, inp->ftype, n, r, z, s, 2.0 * gap_resid);
}

bool
tok_ext_turning_point(
  const struct gkyl_tok_geo_grid_inp *inp, const struct gkyl_tok_geo *geo, double psi, bool upper,
  double *rturn, double *zturn
)
{
  double za = geo->zmaxis;
  double zb = upper ? geo->rzgrid.upper[1] : geo->rzgrid.lower[1];
  double zvalid = za, zinvalid = zb;
  bool found_valid = false, found_invalid = false;
  const int nsamp = 1024;
  for (int i = 0; i <= nsamp; ++i) {
    double f = i / (double)nsamp;
    double z = za + f * (zb - za);
    double R[8] = {0.0}, dRdZ[8] = {0.0};
    double dR[8] = {0.0}, dZ[8] = {0.0};
    int nr = gkyl_tok_geo_R_psiZ(geo, psi, z, 8, R, dRdZ, dR, dZ);
    if (nr > 0) {
      zvalid = z;
      found_valid = true;
    } else if (found_valid) {
      zinvalid = z;
      found_invalid = true;
      break;
    }
  }
  if (!found_valid) {
    return false;
  }
  if (found_invalid) {
    double zv = zvalid, zi = zinvalid;
    for (int k = 0; k < 70; ++k) {
      double zm = 0.5 * (zv + zi);
      double R[8] = {0.0}, dRdZ[8] = {0.0};
      double dR[8] = {0.0}, dZ[8] = {0.0};
      int nr = gkyl_tok_geo_R_psiZ(geo, psi, zm, 8, R, dRdZ, dR, dZ);
      if (nr > 0) {
        zv = zm;
      } else {
        zi = zm;
      }
    }
    zvalid = zv;
  }
  double R[8] = {0.0}, dRdZ[8] = {0.0};
  double dR[8] = {0.0}, dZ[8] = {0.0};
  int nr = gkyl_tok_geo_R_psiZ(geo, psi, zvalid, 8, R, dRdZ, dR, dZ);
  if (nr <= 0) {
    return false;
  }
  *zturn = zvalid;
  *rturn = tok_nearest_value(0.5 * (inp->rleft + inp->rright), R, nr);
  return isfinite(*rturn) && isfinite(*zturn);
}

static bool tok_ext_sample_closed_base(
  const struct gkyl_tok_geo *geo, double psi, const double *br, const double *bz, const double *bs,
  int nb, double target, double *r, double *z
);

// Route a boundary that crosses one of the surface's Z turning points, where no
// single parameterization works: R is degenerate at the midplane end and Z is
// degenerate at the turning point itself.  Split there and build a Z-monotone
// branch on each side, then resample the joined path by arc length.
//
// `upper` selects which turning point, and the two branches sit on opposite
// sides of it: LSN_SOL_MID runs outboard->inboard over the top (true, false),
// while a half-domain CORE_R runs inboard->outboard under the bottom
// (false, true).  `seed_names_side` is the block's own topology property; see
// the branch-A note below for why this route needs it.
static bool
tok_ext_build_via_turning_trace(
  const struct gkyl_tok_geo_grid_inp *inp, const struct gkyl_tok_geo *geo, double psi, int n,
  double r0, double z0, double r1, double z1, bool upper, bool seed_names_side, double *r,
  double *z, double *s, bool *param_is_r
)
{
  double rt = 0.0, zt = 0.0;
  if (!tok_ext_turning_point(inp, geo, psi, upper, &rt, &zt)) {
    return false;
  }
  // Build both monotone-Z branches at the same high resolution used by the
  // closed-core route, then resample the joined path by normalized arc
  // length.  At the single-null separatrix this is exactly the same contour
  // construction used on the neighboring full-core block, so their entire
  // shared radial edge (not only its X-point endpoints) is conforming.
  int nside = n, nb = 2 * nside - 1;
  double *br = gkyl_malloc(sizeof(double[nb]));
  double *bz = gkyl_malloc(sizeof(double[nb]));
  double *bs = gkyl_malloc(sizeof(double[nb]));
  double gap_a = 0.0, gap_b = 0.0;
  // Branch B starts AT the turning point -- the one place the inboard and
  // outboard roots coalesce -- so its seed cannot name a side, exactly the
  // condition the closed-core builder passes true for.
  //
  // Branch A starts at a real endpoint, which was read as a trustworthy seed
  // for continuity.  It is not, for the block that takes this route: both of
  // LSN_SOL_MID's endpoints are rays off the SAME X point, so on the
  // separatrix branch A's seed IS the X point, where the two flanks meet and
  // continuity has nothing to continue from.  `seed_names_side` (true there by
  // topology) hands the tie to the side the block asked for, the same way the
  // closed-core builder resolves it.
  bool ok =
    tok_ext_z_branch_points(
      inp, geo, psi, nside, r0, z0, rt, zt, upper, seed_names_side, br, bz, &gap_a
    ) &&
    tok_ext_z_branch_points(
      inp, geo, psi, nside, rt, zt, r1, z1, !upper, true, br + nside - 1, bz + nside - 1, &gap_b
    );
  double gap_resid = fmax(gap_a, gap_b);
  if (ok) {
    bs[0] = 0.0;
    for (int i = 1; i < nb; ++i) {
      bs[i] = bs[i - 1] + hypot(br[i] - br[i - 1], bz[i] - bz[i - 1]);
    }
    double total = bs[nb - 1];
    ok = isfinite(total) && total > 0.0;
    if (ok) {
      r[0] = r0;
      z[0] = z0;
      for (int i = 1; i < n - 1 && ok; ++i) {
        ok = tok_ext_sample_closed_base(
          geo, psi, br, bz, bs, nb, total * i / (double)(n - 1), &r[i], &z[i]
        );
      }
      r[n - 1] = r1;
      z[n - 1] = z1;
    }
  }
  gkyl_free(br);
  gkyl_free(bz);
  gkyl_free(bs);
  if (!ok) {
    return false;
  }
  *param_is_r = false;
  return tok_ext_finalize_trace_tol(geo, psi, inp->ftype, n, r, z, s, 2.0 * gap_resid);
}

static bool
tok_ext_sample_closed_base(
  const struct gkyl_tok_geo *geo, double psi, const double *br, const double *bz, const double *bs,
  int nb, double target, double *r, double *z
)
{
  double total = bs[nb - 1];
  while (target < 0.0) {
    target += total;
  }
  while (target >= total) {
    target -= total;
  }
  int lo = 0, hi = nb - 1;
  while (hi - lo > 1) {
    int mid = (lo + hi) / 2;
    if (bs[mid] <= target) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  double ds = bs[hi] - bs[lo];
  double w = ds > 0.0 ? (target - bs[lo]) / ds : 0.0;
  double rlin = br[lo] + w * (br[hi] - br[lo]);
  double zlin = bz[lo] + w * (bz[hi] - bz[lo]);
  double R[16] = {0.0}, dRdZ[16] = {0.0};
  double dR[16] = {0.0}, dZ[16] = {0.0};
  int nr = gkyl_tok_geo_R_psiZ(geo, psi, zlin, 16, R, dRdZ, dR, dZ);
  if (nr <= 0) {
    return false;
  }
  *r = tok_nearest_value(rlin, R, nr);
  *z = zlin;
  return true;
}

static bool
tok_ext_build_closed_core_trace(
  const struct gkyl_tok_geo_grid_inp *inp, const struct gkyl_tok_geo *geo, double psi,
  bool separatrix, int n, double rseam, double zseam, double *r, double *z, double *s,
  bool *param_is_r
)
{
  double rlo = 0.0, zlo = 0.0, rup = 0.0, zup = 0.0;
  if (separatrix) {
    rlo = rseam;
    zlo = zseam;
  } else if (!tok_ext_turning_point(inp, geo, psi, false, &rlo, &zlo)) {
    return false;
  }
  if (!tok_ext_turning_point(inp, geo, psi, true, &rup, &zup)) {
    return false;
  }

  int nside = n;
  int nb = 2 * nside - 1;
  double *br = gkyl_malloc(sizeof(double[nb]));
  double *bz = gkyl_malloc(sizeof(double[nb]));
  double *bs = gkyl_malloc(sizeof(double[nb]));
  double gap_a = 0.0, gap_b = 0.0;
  bool ok =
    tok_ext_z_branch_points(inp, geo, psi, nside, rlo, zlo, rup, zup, true, true, br, bz, &gap_a) &&
    tok_ext_z_branch_points(
      inp, geo, psi, nside, rup, zup, rlo, zlo, false, true, br + nside - 1, bz + nside - 1, &gap_b
    );
  double gap_resid = fmax(gap_a, gap_b);
  if (!ok) {
    gkyl_free(br);
    gkyl_free(bz);
    gkyl_free(bs);
    return false;
  }
  bs[0] = 0.0;
  for (int i = 1; i < nb; ++i) {
    bs[i] = bs[i - 1] + hypot(br[i] - br[i - 1], bz[i] - bz[i - 1]);
  }
  double total = bs[nb - 1];
  if (!(total > 0.0)) {
    gkyl_free(br);
    gkyl_free(bz);
    gkyl_free(bs);
    return false;
  }

  double best_d2 = DBL_MAX, seam_s = 0.0;
  for (int i = 0; i < nb - 1; ++i) {
    double dr = br[i + 1] - br[i], dz = bz[i + 1] - bz[i];
    double den = dr * dr + dz * dz;
    double t = den > 0.0 ? ((rseam - br[i]) * dr + (zseam - bz[i]) * dz) / den : 0.0;
    t = fmin(1.0, fmax(0.0, t));
    double rp = br[i] + t * dr, zp = bz[i] + t * dz;
    double d2 = SQ(rseam - rp) + SQ(zseam - zp);
    if (d2 < best_d2) {
      best_d2 = d2;
      seam_s = bs[i] + t * (bs[i + 1] - bs[i]);
    }
  }
  double cell_diag = hypot(geo->rzgrid.dx[0], geo->rzgrid.dx[1]);
  if (sqrt(best_d2) > 2.0 * cell_diag) {
    fprintf(
      stderr, "TOK_EXT_TRACE core seam is not on contour ftype=%d psi=%.17g distance=%.17g\n",
      inp->ftype, psi, sqrt(best_d2)
    );
    gkyl_free(br);
    gkyl_free(bz);
    gkyl_free(bs);
    return false;
  }

  r[0] = rseam;
  z[0] = zseam;
  for (int i = 1; i < n - 1; ++i) {
    double target = seam_s + total * i / (double)(n - 1);
    if (!tok_ext_sample_closed_base(geo, psi, br, bz, bs, nb, target, &r[i], &z[i])) {
      gkyl_free(br);
      gkyl_free(bz);
      gkyl_free(bs);
      return false;
    }
  }
  r[n - 1] = rseam;
  z[n - 1] = zseam;
  gkyl_free(br);
  gkyl_free(bz);
  gkyl_free(bs);
  *param_is_r = false;
  return tok_ext_finalize_trace_tol(geo, psi, inp->ftype, n, r, z, s, 2.0 * gap_resid);
}

// tok_ext_turning_point bisects on "does a root still exist here", and the root
// finder loses the two flanks a hair BEFORE they truly merge -- so the R it
// reports is whichever survivor lies nearer the block's mid-R reference, a
// point ON ONE FLANK rather than the apex.  A single branch seeded there is
// fine, which is why the closed-core builder never had to care.  Two branches
// seeded there are not: both start on the same flank, and the second one's
// first step crosses the surface.  Measured on STEP's psi=1.5343065 surface,
// that put an 0.0805 m step -- exactly the flank gap -- at trace index 8083.
//
// Take the midpoint of the surviving pair instead.  At the Z where the solver
// is about to lose them that IS the apex, to the same accuracy the turning
// point itself is known, and it separates the flanks by construction: one root
// lies above it and one below, so each branch's side filter picks its own.
static void
tok_ext_turn_apex(const struct gkyl_tok_geo *geo, double psi, double z, double *rturn)
{
  double R[16] = {0.0}, dRdZ[16] = {0.0};
  double dR[16] = {0.0}, dZ[16] = {0.0};
  int nr = gkyl_tok_geo_R_psiZ(geo, psi, z, 16, R, dRdZ, dR, dZ);
  if (nr < 2) {
    return;
  }
  double rmin = R[0], rmax = R[0];
  for (int k = 1; k < nr; ++k) {
    rmin = fmin(rmin, R[k]);
    rmax = fmax(rmax, R[k]);
  }
  if (isfinite(rmin) && isfinite(rmax)) {
    *rturn = 0.5 * (rmin + rmax);
  }
}

// The root on the far side of `rnear` at this Z, nearest to it.  On a closed
// core surface there are exactly two roots at any Z strictly inside its range,
// so this names the opposite flank without any tolerance.  It returns false
// where the surface has no root on that side, which is what happens AT a
// turning point -- and, on the separatrix, at the X point, where the two
// flanks meet.
static bool
tok_ext_far_root_at_z(
  const struct gkyl_tok_geo *geo, double psi, double z, double rnear, double rapex, bool outboard,
  double *rfar, int ftype
)
{
  double R[16] = {0.0}, dRdZ[16] = {0.0};
  double dR[16] = {0.0}, dZ[16] = {0.0};
  int nr = gkyl_tok_geo_R_psiZ(geo, psi, z, 16, R, dRdZ, dR, dZ);
  bool found = false;
  double best = 0.0;
  for (int k = 0; k < nr; ++k) {
    // Separate the flanks by the APEX, not by the near root.  The root finder
    // can return the same contour root twice -- measured on STEP's psi=1.5843
    // surface, the two "roots" at the block endpoint's own Z were 2.5876752660
    // 152875 and 2.5876752660152911, four ulp apart -- and a bare R[k] > rnear
    // test takes that duplicate for the opposite flank.  The branch then starts
    // where it should have ended and plain continuity carries it round the
    // INBOARD side: its length came back bit-identical to the inboard branch's,
    // 12.0111965 m for both.  The apex is the one R that genuinely divides the
    // two flanks, so testing against it cannot be fooled by a split root, and
    // it needs no tolerance.
    if (outboard ? !(R[k] > rapex) : !(R[k] < rapex)) {
      continue;
    }
    if (!found || fabs(R[k] - rnear) < fabs(best - rnear)) {
      best = R[k];
      found = true;
    }
  }
  if (!found) {
    fprintf(
      stderr,
      "TOK_EXT_CORE_HALF_NOFAR ftype=%d psi=%.17g z=%.17g rnear=%.17g "
      "rapex=%.17g outboard=%d nr=%d roots=",
      ftype, psi, z, rnear, rapex, (int)outboard, nr
    );
    for (int k = 0; k < nr; ++k) {
      fprintf(stderr, " %.17g", R[k]);
    }
    fprintf(stderr, "\n");
  }
  if (found) {
    *rfar = best;
  }
  return found;
}

// Z extents this small cannot carry a Z-monotone branch; it is the same test
// tok_ext_z_branch_points applies to its own endpoints, stated once so the
// caller can ask "is there a cap here?" instead of building one and failing.
static bool
tok_ext_z_extent_usable(double za, double zb)
{
  return fabs(zb - za) > 64.0 * DBL_EPSILON * fmax(1.0, fmax(fabs(za), fabs(zb)));
}

// Append a Z-monotone branch to a polyline, optionally reversed.  Branches are
// always BUILT outward from a turning point, because that is the only seed
// whose side the caller can name, and reversed on the way in where the
// polyline needs the other order.
static bool
tok_ext_append_branch(
  const struct gkyl_tok_geo_grid_inp *inp, const struct gkyl_tok_geo *geo, double psi, int nside,
  double ra, double za, double rb, double zb, bool outboard, bool seed_names_side, bool reverse,
  double *pr, double *pz, int *np, int cap, double *gap_resid
)
{
  if (nside < 2 || *np + nside - 1 > cap) {
    return false;
  }
  double *br = gkyl_malloc(sizeof(double[nside]));
  double *bz = gkyl_malloc(sizeof(double[nside]));
  double gap = 0.0;
  bool ok = tok_ext_z_branch_points(
    inp, geo, psi, nside, ra, za, rb, zb, outboard, seed_names_side, br, bz, &gap
  );
  if (ok && *np > 0) {
    // The branch's leading end repeats the polyline's current last point, and
    // is skipped below.  If it does not, the segments were wired in the wrong
    // order and the polyline would carry a chord across the surface -- say so
    // rather than emit one.
    double jr = reverse ? br[nside - 1] : br[0];
    double jz = reverse ? bz[nside - 1] : bz[0];
    if (hypot(jr - pr[*np - 1], jz - pz[*np - 1]) > 0.0) {
      fprintf(
        stderr,
        "TOK_EXT_TRACE core_half seam gap ftype=%d psi=%.17g reverse=%d "
        "join=(%.17g,%.17g) last=(%.17g,%.17g)\n",
        inp->ftype, psi, (int)reverse, jr, jz, pr[*np - 1], pz[*np - 1]
      );
      ok = false;
    }
  }
  if (ok) {
    *gap_resid = fmax(*gap_resid, gap);
    // The first point repeats the polyline's current last point, so skip it.
    for (int i = 1; i < nside; ++i) {
      int k = reverse ? nside - 1 - i : i;
      pr[*np] = br[k];
      pz[*np] = bz[k];
      ++(*np);
    }
  }
  gkyl_free(br);
  gkyl_free(bz);
  return ok;
}

// Sample a polyline by arc length, pulling the chord point back onto psi=const
// ALONG GRAD PSI rather than re-solving for a root at fixed Z.
//
// A fixed-Z solve has to choose between the two flanks, and within a few
// millimetres of a turning point the chord point can sit nearer the FAR one --
// so a sampled trace acquires a step the size of the whole flank gap even
// though the polyline it came from is continuous.  Measured on STEP's
// psi=1.5343065 surface: an 0.0805 m step at index 8083, which is exactly the
// flank gap there, with every polyline join verified exact.  Projecting along
// the gradient moves the point perpendicular to the contour, so it stays on the
// segment it was interpolated from and never has to pick a side.
static bool
tok_ext_sample_polyline_grad(
  const struct gkyl_tok_geo *geo, double psi, const double *pr, const double *pz, const double *ps,
  int np, double target, double *r, double *z
)
{
  double total = ps[np - 1];
  if (!(total > 0.0)) {
    return false;
  }
  while (target < 0.0) {
    target += total;
  }
  while (target >= total) {
    target -= total;
  }
  int lo = 0, hi = np - 1;
  while (hi - lo > 1) {
    int mid = (lo + hi) / 2;
    if (ps[mid] <= target) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  double seg = ps[hi] - ps[lo];
  double w = seg > 0.0 ? (target - ps[lo]) / seg : 0.0;
  double rr = pr[lo] + w * (pr[hi] - pr[lo]);
  double zz = pz[lo] + w * (pz[hi] - pz[lo]);
  const double ptol = 1e-12 * fmax(1.0, fabs(psi));
  for (int k = 0; k < 8; ++k) {
    double f = tok_eval_psi_rz_local(geo, rr, zz) - psi;
    double gr = 0.0, gz = 0.0;
    if (!isfinite(f) || !tok_eval_psi_grad_rz_local(geo, rr, zz, &gr, &gz)) {
      return false;
    }
    double g2 = gr * gr + gz * gz;
    if (!(g2 > 0.0)) {
      return false;
    }
    rr -= f / g2 * gr;
    zz -= f / g2 * gz;
    if (fabs(f) <= ptol) {
      break;
    }
  }
  *r = rr;
  *z = zz;
  return isfinite(rr) && isfinite(zz);
}

// CORE_R and CORE_L are the two halves of one CLOSED contour, cut at the two
// X-point rays.  On the SEPARATRIX the cut points are the X points, which are
// also the surface's Z extrema, so each half is Z-monotone and a plain branch
// march expresses it -- which is why TOK_EXT_ROUTE_OUTBOARD has always worked
// there.  On every INTERIOR surface it does not.  The ray lands on the point of
// the surface CLOSEST to the X point, and that sits on the inboard flank ABOVE
// the surface's lowest point.  Measured on STEP's psi=1.6093065 surface: the
// seed is (R,Z)=(2.6045214,-5.8045062), and the surface's two roots at that
// same Z are 2.6045214 and 2.7185442 -- so the contour continues BELOW the seed
// (9.4 mm, to R=2.66153) before turning back up.  The outboard half therefore
// rounds BOTH turning points, and no Z-monotone march can express it: the
// outboard march jumps 0.119 m clean across the surface at its first station.
//
// That jump is what the discontinuity guard catches, but only once the trace is
// dense enough for it to stand out -- the jump is a fixed 0.119 m while the
// mean step falls as 1/n, so the ratio grows without bound.  At the shipped 257
// nodes it passed (ratio 2.5 against the 32 bar) and the trace shipped with a
// chord across the cap; at 8193 it fails (ratio 71.9), falls back to the
// side-blind generic route, and that route picks the INBOARD half -- which is
// how raising the trace resolution turned a quiet 2.3 mm chord error into
// STEP's whole core tracing the wrong flank.
//
// So build the closed contour and CUT it, rather than marching the half.  The
// polyline runs lower turn -> outboard -> upper turn -> inboard -> lower turn,
// and the block's arc is simply the walk FORWARD from its lower endpoint to its
// upper endpoint: for the outboard block that wraps through zero and covers the
// whole outboard branch plus the two inboard slivers between the endpoints and
// the turns; for the inboard block it is the inboard bulk.  The side never has
// to be named -- it is already in the order of the block's own endpoints.
//
// Every branch is split where a turning point or an endpoint falls, so no
// branch has to represent a turning point in its interior, and each cap gets
// its own nside stations across its own few millimetres of Z rather than
// sharing them with the whole surface.
static bool
tok_ext_build_core_half_trace(
  const struct gkyl_tok_geo_grid_inp *inp, const struct gkyl_tok_geo *geo, double psi,
  bool separatrix, int n, double r0, double z0, double r1, double z1, double *r, double *z,
  double *s, bool *param_is_r
)
{
  // Order the block's two endpoints by Z: which of them is the block's LOWER
  // theta end depends on the block, but the contour's own low end does not.
  bool lower_is_first = z0 <= z1;
  double rlo_end = lower_is_first ? r0 : r1, zlo_end = lower_is_first ? z0 : z1;
  double rup_end = lower_is_first ? r1 : r0, zup_end = lower_is_first ? z1 : z0;

  // On the separatrix the endpoints ARE the turning points: the X point is
  // where the two flanks meet, so the scan below -- which looks for the last Z
  // carrying a root -- would sail past it down a divertor leg.  Off the
  // separatrix the contour is closed and the scan is exact.
  double rlo_t = rlo_end, zlo_t = zlo_end, rup_t = rup_end, zup_t = zup_end;
  if (!separatrix) {
    if (!tok_ext_turning_point(inp, geo, psi, false, &rlo_t, &zlo_t) ||
        !tok_ext_turning_point(inp, geo, psi, true, &rup_t, &zup_t)) {
      return false;
    }
    tok_ext_turn_apex(geo, psi, zlo_t, &rlo_t);
    tok_ext_turn_apex(geo, psi, zup_t, &rup_t);
  }

  // A cap exists at an endpoint exactly when the contour still has a root on
  // the far flank at that endpoint's own Z.  At a turning point it does not,
  // and on the separatrix that is the X point -- so the separatrix builds no
  // caps and the construction degenerates continuously to the plain march.
  double rlo_out = 0.0, rup_out = 0.0;
  bool cap_lo =
    !separatrix && tok_ext_z_extent_usable(zlo_end, zlo_t) &&
    tok_ext_far_root_at_z(geo, psi, zlo_end, rlo_end, rlo_t, true, &rlo_out, inp->ftype);
  bool cap_up =
    !separatrix && tok_ext_z_extent_usable(zup_end, zup_t) &&
    tok_ext_far_root_at_z(geo, psi, zup_end, rup_end, rup_t, true, &rup_out, inp->ftype);

  int nside = n;
  int cap = 6 * nside + 8, np = 0;
  double *pr = gkyl_malloc(sizeof(double[cap]));
  double *pz = gkyl_malloc(sizeof(double[cap]));
  double gap_resid = 0.0;
  // s of the two cut points, filled in as the polyline passes through them.
  int i_lo_end = -1, i_up_end = -1;
  // Where each segment ends, so a short one can be named rather than guessed at.
  int seg_np[8] = {0}, nseg = 0;

  pr[0] = rlo_t;
  pz[0] = zlo_t;
  np = 1;
  bool ok = true;
  // --- outboard flank, lower turn up to the upper turn ---
  if (ok && cap_lo) {
    ok = tok_ext_append_branch(
      inp, geo, psi, nside, rlo_t, zlo_t, rlo_out, zlo_end, true, true, false, pr, pz, &np, cap,
      &gap_resid
    );
    seg_np[nseg++] = np;
  }
  if (ok) {
    ok = tok_ext_append_branch(
      inp, geo, psi, nside, cap_lo ? rlo_out : rlo_t, cap_lo ? zlo_end : zlo_t,
      cap_up ? rup_out : rup_t, cap_up ? zup_end : zup_t, true, !cap_lo, false, pr, pz, &np, cap,
      &gap_resid
    );
    seg_np[nseg++] = np;
  }
  if (ok && cap_up)
  // Built from the turn outward and laid down reversed: the turn is the only
  // seed on this segment whose side can be named.
  {
    ok = tok_ext_append_branch(
      inp, geo, psi, nside, rup_t, zup_t, rup_out, zup_end, true, true, true, pr, pz, &np, cap,
      &gap_resid
    );
    seg_np[nseg++] = np;
  }
  // --- inboard flank, upper turn back down to the lower turn ---
  if (ok && cap_up) {
    ok = tok_ext_append_branch(
      inp, geo, psi, nside, rup_t, zup_t, rup_end, zup_end, false, true, false, pr, pz, &np, cap,
      &gap_resid
    );
    seg_np[nseg++] = np;
  }
  if (ok) {
    i_up_end = np - 1;
  }
  if (ok) {
    ok = tok_ext_append_branch(
      inp, geo, psi, nside, rup_end, zup_end, rlo_end, zlo_end, false, !cap_up, false, pr, pz, &np,
      cap, &gap_resid
    );
    seg_np[nseg++] = np;
  }
  if (ok) {
    i_lo_end = np - 1;
  }
  if (ok && cap_lo) {
    ok = tok_ext_append_branch(
      inp, geo, psi, nside, rlo_t, zlo_t, rlo_end, zlo_end, false, true, true, pr, pz, &np, cap,
      &gap_resid
    );
    seg_np[nseg++] = np;
  }
  if (!ok || i_lo_end < 0 || i_up_end < 0 || np < 3) {
    fprintf(
      stderr,
      "TOK_EXT_TRACE core_half build failed ftype=%d psi=%.17g separatrix=%d "
      "cap_lo=%d cap_up=%d np=%d\n",
      inp->ftype, psi, (int)separatrix, (int)cap_lo, (int)cap_up, np
    );
    gkyl_free(pr);
    gkyl_free(pz);
    return false;
  }

  double *ps = gkyl_malloc(sizeof(double[np]));
  ps[0] = 0.0;
  for (int i = 1; i < np; ++i) {
    ps[i] = ps[i - 1] + hypot(pr[i] - pr[i - 1], pz[i] - pz[i - 1]);
  }
  double total = ps[np - 1];
  // The polyline closes on itself, so its last point repeats its first.
  double s_lo = ps[i_lo_end], s_up = ps[i_up_end];
  double s_start = lower_is_first ? s_lo : s_up;
  double s_end = lower_is_first ? s_up : s_lo;
  double arc = s_end - s_start;
  while (arc <= 0.0) {
    arc += total;
  }
  if (!(total > 0.0) || !(arc > 0.0) || arc >= total) {
    fprintf(
      stderr,
      "TOK_EXT_TRACE core_half arc ftype=%d psi=%.17g total=%.17g arc=%.17g "
      "s_lo=%.17g s_up=%.17g\n",
      inp->ftype, psi, total, arc, s_lo, s_up
    );
    gkyl_free(pr);
    gkyl_free(pz);
    gkyl_free(ps);
    return false;
  }

  r[0] = r0;
  z[0] = z0;
  for (int i = 1; i < n - 1 && ok; ++i) {
    ok = tok_ext_sample_polyline_grad(
      geo, psi, pr, pz, ps, np, s_start + arc * i / (double)(n - 1), &r[i], &z[i]
    );
  }
  r[n - 1] = r1;
  z[n - 1] = z1;
  gkyl_free(pr);
  gkyl_free(pz);
  gkyl_free(ps);
  if (!ok) {
    return false;
  }
  *param_is_r = false;
  return tok_ext_finalize_trace_tol(geo, psi, inp->ftype, n, r, z, s, 2.0 * gap_resid);
}

static bool
tok_ext_build_domain_trace(
  const struct gkyl_tok_geo_grid_inp *inp, struct arc_length_ctx *arc_ctx, double psi,
  bool separatrix, double *r, double *z, double *s, int *nout, bool *param_is_r, bool *closed
)
{
  struct tok_ext_topology top;
  if (!tok_ext_topology_from_ftype(inp->ftype, inp->half_domain, &top)) {
    fprintf(
      stderr, "TOK_ORDERED_MAP straight_xpt_ray is unsupported for full-domain ftype=%d\n",
      inp->ftype
    );
    return false;
  }
  double r0 = 0.0, z0 = 0.0, r1 = 0.0, z1 = 0.0;
  bool amb_lo = false, amb_up = false;
  bool lower_ok =
    tok_ext_endpoint_point(inp, arc_ctx, &top.lower, psi, separatrix, &r0, &z0, &amb_lo);
  bool upper_ok =
    tok_ext_endpoint_point(inp, arc_ctx, &top.upper, psi, separatrix, &r1, &z1, &amb_up);
  if (!lower_ok || !upper_ok) {
    fprintf(
      stderr,
      "TOK_ORDERED_MAP failed boundary endpoint ftype=%d psi=%.17g separatrix=%d "
      "lower_ok=%d(kind=%d) upper_ok=%d(kind=%d)\n",
      inp->ftype, psi, separatrix, (int)lower_ok, (int)top.lower.kind, (int)upper_ok,
      (int)top.upper.kind
    );
    return false;
  }
  int n = tok_reference_trace_nodes(inp, arc_ctx->sep_trace_capacity);
  bool ok = false;
  // Once any surface in this block has needed the generic route, use it for
  // every surface: see ext_force_generic_route in the context definition.
  enum tok_ext_route route = top.route;
  if (arc_ctx->ext_force_generic_route && !top.closed) {
    route = TOK_EXT_ROUTE_GENERIC;
  }
  switch (route) {
    case TOK_EXT_ROUTE_GENERIC:
      ok = tok_ext_build_open_trace(inp, arc_ctx->geo, psi, n, r0, z0, r1, z1, r, z, s, param_is_r);
      break;
    case TOK_EXT_ROUTE_OUTBOARD:
      ok = tok_ext_build_z_branch_trace(
        inp, arc_ctx->geo, psi, n, r0, z0, r1, z1, true, top.seed_names_side, r, z, s, param_is_r
      );
      break;
    case TOK_EXT_ROUTE_INBOARD:
      ok = tok_ext_build_z_branch_trace(
        inp, arc_ctx->geo, psi, n, r0, z0, r1, z1, false, top.seed_names_side, r, z, s, param_is_r
      );
      break;
    case TOK_EXT_ROUTE_VIA_UPPER:
      ok = tok_ext_build_via_turning_trace(
        inp, arc_ctx->geo, psi, n, r0, z0, r1, z1, true, top.seed_names_side, r, z, s, param_is_r
      );
      break;
    case TOK_EXT_ROUTE_CLOSED_CORE:
      ok = tok_ext_build_closed_core_trace(
        inp, arc_ctx->geo, psi, separatrix, n, r0, z0, r, z, s, param_is_r
      );
      break;
    case TOK_EXT_ROUTE_CORE_HALF:
      ok = tok_ext_build_core_half_trace(
        inp, arc_ctx->geo, psi, separatrix, n, r0, z0, r1, z1, r, z, s, param_is_r
      );
      break;
  }
  if (!ok && !top.closed && route != TOK_EXT_ROUTE_GENERIC) {
    // Each topology names the parameterization that suits it, but a suitable
    // parameterization is not always a usable one: on a separatrix whose
    // outboard branch leaves the X point almost horizontally, the Z-uniform
    // route puts a fifth of the trace's R range into its first step.  The
    // generic route scores four parameterizations against each other, so let
    // it arbitrate rather than failing the block outright.
    ok = tok_ext_build_open_trace(inp, arc_ctx->geo, psi, n, r0, z0, r1, z1, r, z, s, param_is_r);
    if (ok) {
      route = TOK_EXT_ROUTE_GENERIC;
      arc_ctx->ext_force_generic_route = true;
      fprintf(
        stderr, "TOK_EXT_TRACE route_fallback ftype=%d route=%d psi=%.17g separatrix=%d\n",
        inp->ftype, top.route, psi, separatrix
      );
    }
  }
  arc_ctx->ext_last_trace_used_generic = (route == TOK_EXT_ROUTE_GENERIC);
  if (!ok) {
    // An X-point ray can meet the target surface more than once, and *which*
    // crossing is the right anchor is a question about which connected
    // component of the level set the surface is -- something the ray itself
    // cannot answer.  The C0 quadratic representation grows a spurious ridge
    // along a DG cell face (measured on 204995: a vertical cut peaks exactly at
    // Z = -7*dZ = -0.9625, +7.9e-6 above target, slope flipping sign across the
    // face), and where that ridge tops the requested level it closes a small
    // island.  An anchor on such an island is unreachable: 204995's plate end
    // sits on a 0.971 m open contour while its anchor sits on a 0.071 m closed
    // one, so every parameterization correctly reports no viable candidate.
    //
    // Which crossing is wanted is not decidable from the ray: over the 29-shot
    // regression set only 5 shots cross more than once, all with the same
    // one-excursion-then-the-real-surface shape, and their excursion widths
    // overlap completely -- 204997 needs the FIRST crossing at 10.9-41.9 mm
    // while 204502 and 204995 need the LAST at 4.8-85.6 mm.  Taking the last
    // unconditionally fixes those two and breaks 204997.
    //
    // So let the trace decide, since a crossing on a disconnected island has no
    // route to the other endpoint by construction: keep the first-crossing
    // anchor, and only when it yields a trace no route can build, rebuild the
    // endpoints at the last crossing and try again.  This runs solely where the
    // routine already returns false and aborts the block.
    if (!arc_ctx->ext_ray_use_last_crossing) {
      arc_ctx->ext_ray_use_last_crossing = true;
      bool retry_ok = tok_ext_build_domain_trace(
        inp, arc_ctx, psi, separatrix, r, z, s, nout, param_is_r, closed
      );
      arc_ctx->ext_ray_use_last_crossing = false;
      if (retry_ok) {
        fprintf(
          stderr,
          "TOK_EXT_RAY_RETRY ftype=%d psi=%.17g separatrix=%d "
          "anchored at last ray crossing\n",
          inp->ftype, psi, separatrix
        );
        return true;
      }
    }
    fprintf(
      stderr,
      "TOK_ORDERED_MAP failed extended trace ftype=%d route=%d psi=%.17g separatrix=%d endpoints=(%.17g,%.17g)->(%.17g,%.17g)\n",
      inp->ftype, top.route, psi, separatrix, r0, z0, r1, z1
    );
    return false;
  }
  // The anchor on an X-point ray belongs to BOTH blocks that meet there, so
  // "can I route to it?" is the wrong question -- it is block-specific, and the
  // two sides can answer it differently about the same point.  Measured on
  // 204502 at phase2x: DN_SOL_IN_LO re-anchored at the last crossing on 14 psi
  // surfaces while DN_SOL_IN_MID kept the first on all but one, so one node of
  // the 21 on their shared ray landed 9.9e-02 m (10.4 cells) apart while the
  // other 20 agreed to 1e-16.
  //
  // Ask the symmetric question instead: the anchor must be reachable from the
  // far side of BOTH blocks.  Each side evaluates the same conjunction -- its
  // own trace AND a probe to the peer's far endpoint, which the static topology
  // table supplies -- so the two reach the same verdict by construction, with
  // no cross-block communication and nothing cached.
  //
  // This runs only where the ray genuinely crossed the surface more than once,
  // which is rare, and never on the separatrix, where the anchor is the X point
  // itself.  It cannot loop: the re-anchored rebuild sets
  // ext_ray_use_last_crossing, and this block is skipped while that is set.
  if (!separatrix && (amb_lo || amb_up) && !arc_ctx->ext_ray_use_last_crossing) {
    const struct tok_ext_endpoint *ray = amb_lo ? &top.lower : &top.upper;
    double ar = amb_lo ? r0 : r1, az = amb_lo ? z0 : z1;
    struct tok_ext_endpoint peer_far;
    if (tok_ext_ray_peer_far_endpoint(inp, ray, &peer_far)) {
      double pr = 0.0, pz = 0.0;
      bool peer_ok =
        tok_ext_endpoint_point(inp, arc_ctx, &peer_far, psi, separatrix, &pr, &pz, NULL);
      double *probe = NULL;
      if (peer_ok) {
        probe = gkyl_malloc(sizeof(double[3 * (size_t)n]));
        bool pparam = false;
        peer_ok = tok_ext_build_open_trace(
          inp, arc_ctx->geo, psi, n, ar, az, pr, pz, probe, probe + n, probe + 2 * n, &pparam
        );
        gkyl_free(probe);
      }
      if (!peer_ok) {
        // Re-anchor into scratch: on failure the first-crossing trace already
        // in r/z/s must survive, so it cannot be rebuilt in place.
        double *alt = gkyl_malloc(sizeof(double[3 * (size_t)n]));
        int alt_n = 0;
        bool alt_param = false, alt_closed = false;
        arc_ctx->ext_ray_use_last_crossing = true;
        bool alt_ok = tok_ext_build_domain_trace(
          inp, arc_ctx, psi, separatrix, alt, alt + n, alt + 2 * n, &alt_n, &alt_param, &alt_closed
        );
        arc_ctx->ext_ray_use_last_crossing = false;
        if (alt_ok && alt_n > 0) {
          for (int i = 0; i < alt_n; ++i) {
            r[i] = alt[i];
            z[i] = alt[n + i];
            s[i] = alt[2 * n + i];
          }
          *nout = alt_n;
          *param_is_r = alt_param;
          *closed = alt_closed;
          gkyl_free(alt);
          fprintf(
            stderr,
            "TOK_EXT_RAY_PEER_REANCHOR ftype=%d sector=%d psi=%.17g -- the "
            "first-crossing anchor is unreachable from the block across this "
            "ray; both sides take the last crossing\n",
            inp->ftype, ray->sector, psi
          );
          return true;
        }
        gkyl_free(alt);
      }
    }
  }
  *nout = n;
  *closed = top.closed;
  return true;
}

// Near an X point the true contour can bend sharply within a single
// bracket, so a straight chord between the bracket's two endpoints can end
// up roughly equidistant from two of gkyl_tok_geo_R_psiZ's returned roots
// even though only one is the actual continuation of this trace (confirmed
// directly: at the exact TCV bracket that produces the radial-line fold,
// the bracket's own R-extent, 0.0137, is comparable to the full gap between
// the two roots, 0.0128, while its Z-extent is only 0.0007 -- i.e. the
// bracket sits almost exactly on a Z turning point, the worst case for a
// straight-chord reference).  Both bracket endpoints are themselves
// already-resolved points of this same psi's own raw trace (trustworthy:
// this trace's own construction was directly checked, in an earlier
// investigation of this exact bug, and found not to be the source of the
// ambiguity). Recover a reliable disambiguation reference by walking from
// whichever endpoint is closer in Z to the target, in enough small
// sub-steps that root selection stays unambiguous at every one: a small
// enough step keeps the previous step's R clearly nearer one root than the
// other, so nearest-value tracks the correct branch by construction rather
// than by the luck of where the full chord happens to land. This depends
// only on the current bracket's own two endpoints -- not on any other grid
// corner, and not on the outer grid resolution -- so it neither chains
// errors across neighboring points nor introduces any resolution
// dependence the original per-corner-independent construction did not
// already have.
static bool
tok_trace_walk_root(
  const struct gkyl_tok_geo *geo, double psi, double z_lo, double r_lo, double z_hi, double r_hi,
  double z_target, double *r_out
)
{
  bool from_lo = fabs(z_target - z_lo) <= fabs(z_target - z_hi);
  double z_anchor = from_lo ? z_lo : z_hi;
  double r_walk = from_lo ? r_lo : r_hi;
  const int nsub = 16;
  for (int k = 1; k <= nsub; ++k) {
    double z_k = z_anchor + (z_target - z_anchor) * k / (double)nsub;
    double roots_k[16] = {0.0}, dRdZ_k[16] = {0.0};
    double dR_k[16] = {0.0}, dZ_k[16] = {0.0};
    int nr_k = gkyl_tok_geo_R_psiZ(geo, psi, z_k, 8, roots_k, dRdZ_k, dR_k, dZ_k);
    if (nr_k <= 0) {
      return false;
    }
    r_walk = tok_nearest_value(r_walk, roots_k, nr_k);
  }
  *r_out = r_walk;
  return isfinite(r_walk);
}

// A closed contour generically has two valid R roots at almost every Z
// within its range (an ordinary "inboard side"/"outboard side" pair) -- so
// nr>1 alone is the norm, not evidence of trouble, and rlin (the straight
// chord across the bracket) already reliably favors the correct one in the
// overwhelming majority of these.  Genuine trouble is specifically when
// rlin cannot confidently tell the two roots apart -- it sits comparably
// close to both, the signature found at the one bracket that actually
// produces the fold (root gap 0.0128, rlin only 0.0034 from the midpoint).
// Restricting the (more expensive, and elsewhere unnecessary) walk to only
// this rare case keeps it from re-deciding calls that were already correct.
static bool tok_psi_normal_project(
  const struct gkyl_tok_geo *geo, double psi, double r0, double z0, double max_disp, double tol,
  double *r, double *z
);

static bool
tok_rlin_ambiguous(double rlin, const double *roots, int nr)
{
  if (nr <= 1) {
    return false;
  }
  double best = DBL_MAX, second = DBL_MAX;
  for (int k = 0; k < nr; ++k) {
    double d = fabs(roots[k] - rlin);
    if (d < best) {
      second = best;
      best = d;
    } else if (d < second) {
      second = d;
    }
  }
  return second < 2.0 * best;
}

// A3 (2026-10-02): put a trace sample on the surface at the chord's own arc
// parameter -- solve along the chord's NORMAL, not along a coordinate axis.
//
// A sample is the point at arc fraction w of a bracket of two consecutive,
// on-surface trace points.  The chord point between them already has the
// right arc parameter to first order; what remains is to move it onto the
// surface without sliding along it.  The old rule fixed one coordinate of the
// chord point and solved the other, which slides the point along any surface
// that runs parallel to the fixed axis: measured on NSTX-U 203970 and 203834
// (PF_LO_L separatrix row leaving the X point), the leg runs within 0.02 of
// horizontal for 14-18 mm, the column solve at the chord point missed it
// entirely (nearest root 0.71 m away) and the row solve took the root 5-8 mm
// along it, so the three samples meant for 4.9, 9.8 and 14.8 mm sat at 0.14,
// 18.2 and 18.5 mm: the misplaced Gauss nodes of the B1 tip cell.
//
// Solving along the chord normal is the same idea in the bracket's own frame:
// the chord is the best available estimate of the surface tangent, so its
// normal is the well-conditioned direction on every bracket, and a 1-D Newton
// iteration confined to that line cannot slide.  (A Newton projection along
// the LOCAL gradient was tried and slid 4.9 mm to the X point on 203970: next
// to a saddle on a DG cell face the level set is a sliver of two pieces
// 0.15 mm apart and the gradient between them points along the sliver.)  The
// displacement is bounded by half the bracket, as the axis solves are; when
// the line meets no surface within that -- the chord nearly parallel to the
// normal, an L-shaped bracket -- the axis-aligned rule applies as before.

static bool
tok_chord_normal_solve(
  const struct gkyl_tok_geo *geo, double psi, double rlin, double zlin, double cdr, double cdz,
  double max_disp, double *r, double *z
)
{
  const double clen = hypot(cdr, cdz);
  if (!(clen > 0.0) || !isfinite(clen) || !(max_disp > 0.0)) {
    return false;
  }
  const double nr = -cdz / clen, nz = cdr / clen; // unit normal of the chord
  // Solved to roundoff: |psi - psi_row| down to the double precision of psi,
  // or the bracket collapsed to roundoff in position. The fixed tolerance this
  // replaced (1e-9 in psi) left a point off the contour by up to 1e-9/|grad psi|,
  // which grows without bound beside the X point; arcs measured through such
  // points there came out ~10 nm long or short, a fixed offset, so the theta-seam
  // grading mismatch at the X-point faces grew with theta (NSTX-U: 1.0e-6 at
  // theta x1, 5.4e-6 at x8). Solved to roundoff it is 1e-10 at both.
  const double tol = DBL_EPSILON * fmax(1.0, fabs(psi));
  const double g0 = tok_eval_psi_rz_local(geo, rlin, zlin) - psi;
  if (!isfinite(g0)) {
    return false;
  }
  if (fabs(g0) <= tol) {
    *r = rlin;
    *z = zlin;
    return true;
  }
  // Bracket the nearest sign change on each side of the chord point, stepping
  // out geometrically from a thousandth of the bound; the nearer side wins.
  // (A Newton iteration along the line was tried first and failed on the
  // sliver next to an X point, where the line's derivative vanishes between
  // the two pieces and the clamped steps wander.)  Bisection then needs no
  // derivative at all, so the DG gradient kink at a cell face cannot upset it.
  double ta = 0.0, tb = 0.0, ga = g0, gb = g0;
  bool found = false;
  // A piece shorter than the residual at its chord point can be reached by --
  // psi's own evaluation noise (STEP, psi ~ 1.5: ~4e-15, against a 1e-15 m piece
  // between a node and a trace point it coincides with) or an end point another
  // solve placed less precisely -- shows no sign change within max_disp. The
  // search then goes on to twice the first-order distance to the contour along
  // the search line, |g0|/|n.grad psi| (on so short a chord its normal can be
  // far from the gradient: 65 deg at STEP's midplane reference), before giving
  // up; and a point whose first-order distance is within the roundoff of its
  // own coordinates is on the contour as far as the arithmetic can tell.
  double reach = max_disp, dist = DBL_MAX;
  for (int pass = 0; pass < 2 && !found; ++pass) {
    if (pass == 1) {
      double gr = 0.0, gz = 0.0;
      if (!tok_eval_psi_grad_rz_local(geo, rlin, zlin, &gr, &gz) ||
          !(fabs(gr * nr + gz * nz) > 0.0)) {
        break;
      }
      dist = fabs(g0) / fabs(gr * nr + gz * nz);
      if (!(2.0 * dist > max_disp) || !isfinite(dist)) {
        break;
      }
      reach = 2.0 * dist;
    }
    for (int side = 0; side < 2; ++side) {
      const double sgn = side ? -1.0 : 1.0;
      double tp = 0.0, gp = g0;
      for (double d = max_disp / 1024.0;; d = fmin(2.0 * d, reach)) {
        const double t = sgn * d;
        const double g = tok_eval_psi_rz_local(geo, rlin + t * nr, zlin + t * nz) - psi;
        if (!isfinite(g)) {
          break;
        }
        if (fabs(g) <= tol) { // landed on it
          if (!found || d < fmax(fabs(ta), fabs(tb))) {
            ta = tb = t;
            ga = gb = g;
            found = true;
          }
          break;
        }
        if ((g < 0.0) != (gp < 0.0)) {
          if (!found || d < fmax(fabs(ta), fabs(tb))) {
            ta = tp;
            tb = t;
            ga = gp;
            gb = g;
            found = true;
          }
          break;
        }
        tp = t;
        gp = g;
        if (d >= reach) {
          break;
        }
      }
    }
  }
  if (!found) {
    if (dist <= 4.0 * DBL_EPSILON * fmax(1.0, fmax(fabs(rlin), fabs(zlin)))) {
      *r = rlin;
      *z = zlin;
      return true;
    }
    return false;
  }
  double t = ta, g = ga;
  bool collapsed = false;
  if (ta != tb) {
    // Regula falsi with the Illinois modification: superlinear without a
    // derivative, and the bracket can never be lost.
    int side = 0;
    for (int k = 0; k < 64; ++k) {
      t = (gb != ga) ? (ta * gb - tb * ga) / (gb - ga) : 0.5 * (ta + tb);
      if (!(t > fmin(ta, tb)) || !(t < fmax(ta, tb))) {
        t = 0.5 * (ta + tb);
      }
      g = tok_eval_psi_rz_local(geo, rlin + t * nr, zlin + t * nz) - psi;
      if (!isfinite(g)) {
        return false;
      }
      if (fabs(g) <= tol) {
        break;
      }
      if ((g < 0.0) == (ga < 0.0)) {
        ta = t;
        ga = g;
        if (side == -1) {
          gb *= 0.5;
        }
        side = -1;
      } else {
        tb = t;
        gb = g;
        if (side == 1) {
          ga *= 0.5;
        }
        side = 1;
      }
      if (fabs(tb - ta) <= 4.0 * DBL_EPSILON * fmax(1.0, fabs(t))) {
        collapsed = true;
        break;
      }
    }
  }
  if (!(fabs(g) <= tol) && !collapsed) {
    return false;
  }
  *r = rlin + t * nr;
  *z = zlin + t * nz;
  return true;
}

// ---------------------------------------------------------------------------
// Arc-exact trace sampling (user decision 2026-09-29; built 2026-10-06).
//
// A trace is a list of points on one psi contour with their CUMULATIVE CHORD
// length. Sampling it at a fraction u used to pick the bracket by chord length
// and the point inside it by chord fraction, then move that point onto the
// contour along the chord's normal. Neither is the contour's arc: the chord of
// a bracket is shorter than its arc by ~ (kappa h)^2/24, and a chord fraction
// projected onto a curve is not the same fraction of its arc. Because the
// reference trace's brackets do not shrink with theta refinement (257 nodes),
// the error is ABSOLUTE in physical length -- a fixed (+a, 0, -a) sawtooth of
// 1-5 um on ordinary rows and 20-75 um beside the X point -- so the relative
// error of a cell, and the seam grading |G-1| measured on it, GREW with theta
// refinement on every device (fable-handoff 09 sections 8b, 9c). Refining the
// trace only moves that floor and costs 4-5x.
//
// Here the bracket is chosen by the TRUE cumulative arc of the trace and the
// point slid along the contour until its true arc from the bracket's start is
// the target: u is a fraction of the contour's own length, to roundoff. No
// constant is introduced: the arc of a piece is accepted once halving it
// changes its length by no more than sqrt(DBL_EPSILON) relatively (Richardson's
// correction then leaves an O(eps) error), or once the piece is no longer than
// the precision the contour solve locates points to.
// ---------------------------------------------------------------------------

// The psi contour between two of its points p0 -> p1 close enough that it is
// a graph over their chord (one trace bracket or less), as a polyline exact
// to roundoff: bisected, each midpoint moved onto the contour along the
// chord's normal, until halving a piece changes its length by no more than
// sqrt(DBL_EPSILON) relatively (Richardson's correction then leaves an O(eps)
// error), or the piece is no longer than the contour solve locates points to.
// Appends each piece's END point and the piece's TRUE arc. False if a midpoint
// cannot be put on the contour.
struct tok_leaves {
  double *r, *z, *arc;
  int n, cap;
};

static bool
tok_leaves_push(struct tok_leaves *lv, double r, double z, double arc)
{
  if (lv->n == lv->cap) {
    int cap = lv->cap ? 2 * lv->cap : 1024;
    double *nr = gkyl_malloc(cap * sizeof(double)), *nz = gkyl_malloc(cap * sizeof(double)),
           *na = gkyl_malloc(cap * sizeof(double));
    if (lv->n) {
      memcpy(nr, lv->r, lv->n * sizeof(double));
      memcpy(nz, lv->z, lv->n * sizeof(double));
      memcpy(na, lv->arc, lv->n * sizeof(double));
      gkyl_free(lv->r);
      gkyl_free(lv->z);
      gkyl_free(lv->arc);
    }
    lv->r = nr;
    lv->z = nz;
    lv->arc = na;
    lv->cap = cap;
  }
  lv->r[lv->n] = r;
  lv->z[lv->n] = z;
  lv->arc[lv->n] = arc;
  ++lv->n;
  return true;
}

static bool
tok_contour_leaves(
  const struct gkyl_tok_geo *geo, double psi, double r0, double z0, double r1, double z1, int depth,
  struct tok_leaves *lv
)
{
  const double l1 = hypot(r1 - r0, z1 - z0);
  if (!(l1 > 0.0)) {
    return true;
  }
  // A piece is measured within one cell of the psi representation: where its
  // chord crosses a cell face (a grid line in R or Z, strictly between its
  // ends), it is split where the contour itself crosses that face, solved
  // along the face line so the split point lies on it: psi is one polynomial
  // within a cell but only C1 across, and a piece straddling a face could pass
  // both halvings with a curvature jump unseen (ASDEX inner-leg separatrix: one
  // cell 0.2 nm long at every rung, theta-seam |G-1| flat at 1.2e-9; split,
  // ~1e-10). The split point must lie within the piece's bounding box, so each
  // half's box nests in it and the faces strictly inside can only become
  // fewer: the splitting ends.
  if (depth < 52) {
    const struct gkyl_rect_grid *fg = geo->use_cubics ? &geo->rzgrid_cubic : &geo->rzgrid;
    const double a[2] = {r0, z0}, b[2] = {r1, z1};
    for (int d = 0; d < 2; ++d) {
      const double lo = fmin(a[d], b[d]), hi = fmax(a[d], b[d]);
      const double k = floor((lo - fg->lower[d]) / fg->dx[d]) + 1.0;
      const double face = fg->lower[d] + k * fg->dx[d];
      if (!(face > lo && face < hi)) {
        continue;
      }
      const double t = (face - a[d]) / (b[d] - a[d]);
      double fr = 0.0, fz = 0.0;
      const bool ok =
        d == 0 ?
          tok_chord_normal_solve(geo, psi, face, z0 + t * (z1 - z0), 1.0, 0.0, 0.5 * l1, &fr, &fz) :
          tok_chord_normal_solve(geo, psi, r0 + t * (r1 - r0), face, 0.0, 1.0, 0.5 * l1, &fr, &fz);
      if (!ok || (d == 0 ? fr != face : fz != face)) {
        continue;
      }
      if (!(fr >= fmin(r0, r1) && fr <= fmax(r0, r1) && fz >= fmin(z0, z1) && fz <= fmax(z0, z1))) {
        continue;
      }
      if (!(hypot(fr - r0, fz - z0) > 0.0) || !(hypot(r1 - fr, z1 - fz) > 0.0)) {
        continue;
      }
      return tok_contour_leaves(geo, psi, r0, z0, fr, fz, depth + 1, lv) &&
             tok_contour_leaves(geo, psi, fr, fz, r1, z1, depth + 1, lv);
    }
  }
  double mr = 0.0, mz = 0.0;
  if (!tok_chord_normal_solve(
        geo, psi, 0.5 * (r0 + r1), 0.5 * (z0 + z1), r1 - r0, z1 - z0, 0.5 * l1, &mr, &mz
      )) {
    return false;
  }
  const double l2 = hypot(mr - r0, mz - z0) + hypot(r1 - mr, z1 - mz);
  // The contour solve locates a point to |psi residual|/|grad psi|; a piece no
  // longer than a few of those cannot be resolved further.
  double gr = 0.0, gz = 0.0, located = 0.0;
  if (tok_eval_psi_grad_rz_local(geo, mr, mz, &gr, &gz) && hypot(gr, gz) > 0.0) {
    // no finer than the roundoff of the point's own coordinates: below that a
    // piece is a point, and halving it only bisects roundoff (STEP: unbounded)
    located = fmax(
      DBL_EPSILON * fmax(1.0, fabs(psi)) / hypot(gr, gz),
      DBL_EPSILON * fmax(1.0, fmax(fabs(mr), fabs(mz)))
    );
  }
  // Depth is bounded by halving to the double-precision resolution of l1.
  if (l1 <= 4.0 * located || depth >= 52) {
    const double arc = l2 + (l2 - l1) / 3.0;
    // split the corrected arc between the halves in proportion to their chords
    const double a0 = hypot(mr - r0, mz - z0) / l2 * arc;
    tok_leaves_push(lv, mr, mz, a0);
    tok_leaves_push(lv, r1, z1, arc - a0);
    return true;
  }
  // One halving is blind to a piece with an inflection at its middle: the
  // midpoint then lies on the chord, halving changes nothing, and the chord is
  // taken for the arc. So the piece is accepted only when the next halving
  // agrees as well. Measured on STEP: one bracket beside each X point of CORE_R's
  // separatrix row was taken ~0.5 um short, displacing every node of that row
  // by 0.5 um and lengthening its two end cells by as much, so the theta-seam
  // grading mismatch there grew with theta (1.7e-6 at x4, 7.0e-6 at x16); with
  // the second halving it is 3e-8 and 2e-8.
  // ... and only a piece within one cell of the psi representation: psi is one
  // polynomial there, but only C1 across cells, so a contour crossing several
  // has curvature jumps that three or five points can straddle unseen. Measured
  // on STEP (theta x4, DN_SOL_OUT_MID row 2): a 0.29 m cell spanning several
  // psi cells passed both halvings and was taken 217 nm (7.5e-7) short.
  const struct gkyl_rect_grid *pg = geo->use_cubics ? &geo->rzgrid_cubic : &geo->rzgrid;
  if (fabs(l2 - l1) <= sqrt(DBL_EPSILON) * l2 && l1 <= fmin(pg->dx[0], pg->dx[1])) {
    const double la = hypot(mr - r0, mz - z0), lb = hypot(r1 - mr, z1 - mz);
    double ar = 0.0, az = 0.0, br = 0.0, bz = 0.0;
    if (tok_chord_normal_solve(
          geo, psi, 0.5 * (r0 + mr), 0.5 * (z0 + mz), mr - r0, mz - z0, 0.5 * la, &ar, &az
        ) &&
        tok_chord_normal_solve(
          geo, psi, 0.5 * (mr + r1), 0.5 * (mz + z1), r1 - mr, z1 - mz, 0.5 * lb, &br, &bz
        )) {
      const double c1 = hypot(ar - r0, az - z0), c2 = hypot(mr - ar, mz - az);
      const double c3 = hypot(br - mr, bz - mz), c4 = hypot(r1 - br, z1 - bz);
      const double l4 = c1 + c2 + c3 + c4;
      if (fabs(l4 - l2) <= sqrt(DBL_EPSILON) * l4) {
        // split the corrected arc between the quarters in proportion to their chords
        const double arc = l4 + (l4 - l2) / 3.0;
        tok_leaves_push(lv, ar, az, c1 / l4 * arc);
        tok_leaves_push(lv, mr, mz, c2 / l4 * arc);
        tok_leaves_push(lv, br, bz, c3 / l4 * arc);
        tok_leaves_push(lv, r1, z1, c4 / l4 * arc);
        return true;
      }
    }
  }
  return tok_contour_leaves(geo, psi, r0, z0, mr, mz, depth + 1, lv) &&
         tok_contour_leaves(geo, psi, mr, mz, r1, z1, depth + 1, lv);
}

// A trace refined to its true arc, cached: tok_trace_sample is called many
// times with the same trace (the row rule alone, thousands per row). A few
// slots, keyed by the arrays and a signature of their contents, so a buffer
// refilled for another surface misses. lv.arc holds the CUMULATIVE true arc.
struct tok_arc_cache_slot {
  const double *tr, *tz, *ts;
  int n;
  double psi, sig[4];
  bool valid;
  struct tok_leaves lv;
  unsigned long used;
};
static _Thread_local struct tok_arc_cache_slot tok_arc_cache[4];
static _Thread_local unsigned long tok_arc_cache_clock;

static const struct tok_leaves *
tok_trace_true_arc(
  const struct gkyl_tok_geo *geo, double psi, const double *tr, const double *tz, const double *ts,
  int n
)
{
  const double sig[4] = {tr[0], tz[0], tr[n - 1], ts[n - 1]};
  struct tok_arc_cache_slot *slot = &tok_arc_cache[0];
  for (int k = 0; k < 4; ++k) {
    struct tok_arc_cache_slot *c = &tok_arc_cache[k];
    if (c->valid && c->tr == tr && c->tz == tz && c->ts == ts && c->n == n && c->psi == psi &&
        !memcmp(c->sig, sig, sizeof sig)) {
      c->used = ++tok_arc_cache_clock;
      return &c->lv;
    }
    if (c->used < slot->used) {
      slot = c;
    }
  }
  slot->valid = false;
  slot->lv.n = 0;
  tok_leaves_push(&slot->lv, tr[0], tz[0], 0.0);
  for (int i = 1; i < n; ++i) {
    if (!tok_contour_leaves(
          geo, psi, slot->lv.r[slot->lv.n - 1], slot->lv.z[slot->lv.n - 1], tr[i], tz[i], 0,
          &slot->lv
        )) {
      slot->used = 0;
      return 0;
    }
  }
  for (int j = 1; j < slot->lv.n; ++j) {
    slot->lv.arc[j] += slot->lv.arc[j - 1];
  }
  slot->tr = tr;
  slot->tz = tz;
  slot->ts = ts;
  slot->n = n;
  slot->psi = psi;
  memcpy(slot->sig, sig, sizeof sig);
  slot->valid = true;
  slot->used = ++tok_arc_cache_clock;
  return &slot->lv;
}

// ---------------------------------------------------------------------------
// The field-line angle at the contour's true arc (2026-10-06).
//
// Along a row the angle phi gains integral F/(R |grad psi|) ds. The ordered map
// accumulated it with the midpoint rule at each trace CHORD's midpoint (off the
// surface by the chord's sagitta), then interpolated it linearly -- once in the
// uniform-arc resample, once at the node -- and took dphi/dtheta as the
// bracket's difference quotient. All of that is second order in the trace
// spacing and largest where the integrand peaks, beside the X point, and each
// row's trace samples it at its own points, so the error differed row to row:
// the streaks left in g^12, g^22, g^23 once R and Z were exact. Measured on
// STEP (psi x8, theta x8): map-trace multipliers 4 / 8 / 16 left R and Z
// unchanged to 1e-15 and took the g^12 cells the row test flags from 280 to 46
// to 0 (g^23: 286 / 142 / 2).
//
// Here the integral between two points of the contour is taken on the contour:
// two-point Gauss-Legendre on a piece, its nodes moved onto the contour along
// the chord's normal as tok_contour_leaves moves its midpoints, the piece
// bisected until halving changes the value by no more than sqrt(DBL_EPSILON)
// relatively, or the piece is no longer than the contour solve locates points
// to -- tok_contour_leaves' rules, so no constant is introduced. The rule never
// evaluates a piece's end points, so a piece ending at the X point, where
// |grad psi| vanishes and the angle diverges, stays finite. A node's angle is
// its bracket's start plus the integral to the node itself, and dphi/dtheta is
// F/(R |grad psi|) at the node times the arc rate dR/dtheta and dZ/dtheta use.
// ---------------------------------------------------------------------------

// dphi/ds = F/(R |grad psi|) at a point of the contour; false where the
// gradient vanishes.
static bool
tok_phi_rate(const struct gkyl_tok_geo *geo, double fpol, double r, double z, double *f)
{
  double gr = 0.0, gz = 0.0;
  if (!tok_eval_psi_grad_rz_local(geo, r, z, &gr, &gz)) {
    return false;
  }
  const double g = hypot(gr, gz);
  if (!(g > 0.0) || !(r > 0.0)) {
    return false;
  }
  *f = fpol / (r * g);
  return isfinite(*f);
}

// One piece p0 -> p1 of the contour: its midpoint on the contour, its arc
// (Richardson, as in tok_contour_leaves), how finely the contour solve locates
// points there, and the two-point Gauss-Legendre value of the angle.
struct tok_phi_piece {
  double mr, mz, arc, located, val;
};

static bool
tok_phi_piece_eval(
  const struct gkyl_tok_geo *geo, double psi, double fpol, double r0, double z0, double r1,
  double z1, struct tok_phi_piece *pc
)
{
  const double cdr = r1 - r0, cdz = z1 - z0, l1 = hypot(cdr, cdz);
  *pc = (struct tok_phi_piece){r0, z0, 0.0, 0.0, 0.0};
  if (!(l1 > 0.0)) {
    return true;
  }
  if (!tok_chord_normal_solve(
        geo, psi, r0 + 0.5 * cdr, z0 + 0.5 * cdz, cdr, cdz, 0.5 * l1, &pc->mr, &pc->mz
      )) {
    return false;
  }
  const double l2 = hypot(pc->mr - r0, pc->mz - z0) + hypot(r1 - pc->mr, z1 - pc->mz);
  pc->arc = l2 + (l2 - l1) / 3.0;
  double gr = 0.0, gz = 0.0;
  if (tok_eval_psi_grad_rz_local(geo, pc->mr, pc->mz, &gr, &gz) && hypot(gr, gz) > 0.0) {
    pc->located = fmax(
      DBL_EPSILON * fmax(1.0, fabs(psi)) / hypot(gr, gz),
      DBL_EPSILON * fmax(1.0, fmax(fabs(pc->mr), fabs(pc->mz)))
    );
  }
  const double half = 0.5 / sqrt(3.0);
  double sum = 0.0;
  for (int k = 0; k < 2; ++k) {
    const double t = k ? 0.5 + half : 0.5 - half;
    double r = 0.0, z = 0.0, f = 0.0;
    if (!tok_chord_normal_solve(geo, psi, r0 + t * cdr, z0 + t * cdz, cdr, cdz, 0.5 * l1, &r, &z) ||
        !tok_phi_rate(geo, fpol, r, z, &f)) {
      return false;
    }
    sum += f;
  }
  pc->val = 0.5 * pc->arc * sum;
  return true;
}

static bool
tok_phi_integral_rec(
  const struct gkyl_tok_geo *geo, double psi, double fpol, double r0, double z0, double r1,
  double z1, const struct tok_phi_piece *whole, int depth, double *out
)
{
  struct tok_phi_piece a, b;
  if (!tok_phi_piece_eval(geo, psi, fpol, r0, z0, whole->mr, whole->mz, &a) ||
      !tok_phi_piece_eval(geo, psi, fpol, whole->mr, whole->mz, r1, z1, &b)) {
    return false;
  }
  const double halves = a.val + b.val;
  // Halving cuts the two-point rule's error 16-fold; depth is bounded by
  // halving to the double-precision resolution of the piece.
  if (fabs(halves - whole->val) <= sqrt(DBL_EPSILON) * fabs(halves) ||
      hypot(r1 - r0, z1 - z0) <= 4.0 * whole->located || depth >= 52) {
    *out = halves + (halves - whole->val) / 15.0;
    return true;
  }
  double x = 0.0, y = 0.0;
  if (!tok_phi_integral_rec(geo, psi, fpol, r0, z0, whole->mr, whole->mz, &a, depth + 1, &x) ||
      !tok_phi_integral_rec(geo, psi, fpol, whole->mr, whole->mz, r1, z1, &b, depth + 1, &y)) {
    return false;
  }
  *out = x + y;
  return true;
}

// The angle gained along the contour from p0 to p1, two of its points close
// enough that it is a graph over their chord (one trace bracket or less).
static bool
tok_phi_integral(
  const struct gkyl_tok_geo *geo, double psi, double fpol, double r0, double z0, double r1,
  double z1, double *out
)
{
  *out = 0.0;
  if (r0 == r1 && z0 == z1) {
    return true;
  }
  struct tok_phi_piece whole;
  if (!tok_phi_piece_eval(geo, psi, fpol, r0, z0, r1, z1, &whole)) {
    return false;
  }
  return tok_phi_integral_rec(geo, psi, fpol, r0, z0, r1, z1, &whole, 0, out);
}

// The angle gained along a trace from its point (r0, z0) to (r, z), a point of
// the contour in the bracket (r0, z0) -> (r1, z1) whose own increment p01 the
// trace already holds. On the separatrix the angle diverges logarithmically at
// the X point, and a trace of it holds the X point itself; an integral from
// there is a truncated divergent one, cut off wherever the bisection stops,
// and a second such integral, to the node, is cut off elsewhere than the
// trace's own for the bracket -- offsetting the node from its neighbours by up
// to C ln 2 (C the divergence's coefficient, F/(R |H|)): measured 2.8 rad on
// ASDEX at theta x4 at the separatrix node 14 mm from the X point. A bracket
// that starts at an X point is therefore integrated from its other end, so the
// node shares the trace's own cut-off and its angle relative to every other
// point of the row is the finite integral between them.
static bool
tok_phi_into_bracket(
  const struct gkyl_tok_geo *geo, double psi, double fpol, double r0, double z0, double r1,
  double z1, double p01, double r, double z, double *out
)
{
  *out = 0.0;
  if (r == r0 && z == z0) {
    return true;
  }
  for (int k = 0; k < 2; ++k) {
    double xr = 0.0, xz = 0.0;
    if (tok_ext_xpoint_rz(geo, k ? TOK_EXT_UPPER_XPT : TOK_EXT_LOWER_XPT, &xr, &xz) && r0 == xr &&
        z0 == xz) {
      double e = 0.0;
      if (!tok_phi_integral(geo, psi, fpol, r, z, r1, z1, &e)) {
        return false;
      }
      *out = p01 - e;
      return true;
    }
  }
  return tok_phi_integral(geo, psi, fpol, r0, z0, r, z, out);
}

static bool
tok_trace_sample(
  const struct gkyl_tok_geo *geo, double psi, const double *tr, const double *tz, const double *ts,
  int n, bool param_is_r, double u, double *r, double *z
)
{
  if (n < 2 || !(ts[n - 1] > 0.0)) {
    return false;
  }
  // Logical edge coordinates acquire a few ulps of roundoff when reconstructed
  // from the arc-length coordinate.  Snap those values before resampling: near
  // an X point, asking the polynomial root finder for an infinitesimally
  // interior point can select the other root of the saddle instead of the
  // explicitly stored, shared endpoint.
  const double endpoint_tol = 256.0 * DBL_EPSILON;
  if (u <= endpoint_tol) {
    *r = tr[0];
    *z = tz[0];
    return true;
  }
  if (u >= 1.0 - endpoint_tol) {
    *r = tr[n - 1];
    *z = tz[n - 1];
    return true;
  }
  // Arc-exact first (see the arc-exact sampling section): the trace
  // refined to the contour's true arc, the point at its true arc fraction.
  const struct tok_leaves *lv = tok_trace_true_arc(geo, psi, tr, tz, ts, n);
  if (lv && lv->n > 1 && lv->arc[lv->n - 1] > 0.0) {
    const double want = u * lv->arc[lv->n - 1];
    int a = 0, b = lv->n - 1;
    while (b - a > 1) {
      int mid = (a + b) / 2;
      if (lv->arc[mid] < want) {
        a = mid;
      } else {
        b = mid;
      }
    }
    const double whole = lv->arc[b] - lv->arc[a], into = want - lv->arc[a];
    if (!(into > 0.0)) {
      *r = lv->r[a];
      *z = lv->z[a];
      return true;
    }
    if (!(into < whole)) {
      *r = lv->r[b];
      *z = lv->z[b];
      return true;
    }
    // inside a piece the chord and the arc agree to roundoff (see above)
    const double f = into / whole;
    const double cdr = lv->r[b] - lv->r[a], cdz = lv->z[b] - lv->z[a];
    if (tok_chord_normal_solve(
          geo, psi, lv->r[a] + f * cdr, lv->z[a] + f * cdz, cdr, cdz, 0.5 * hypot(cdr, cdz), r, z
        )) {
      return true;
    }
  }
  double target = u * ts[n - 1];
  int lo = 0, hi = n - 1;
  while (hi - lo > 1) {
    int mid = (lo + hi) / 2;
    if (ts[mid] < target) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  double ds = ts[hi] - ts[lo];
  double w = ds > 0.0 ? (target - ts[lo]) / ds : 0.0;
  double rlin = tr[lo] + w * (tr[hi] - tr[lo]);
  double zlin = tz[lo] + w * (tz[hi] - tz[lo]);
  // A bracket that moves mostly in R (a local Z turning point) makes
  // "fix Z, solve R" the ill-conditioned direction: gkyl_tok_geo_R_psiZ can
  // return two roots that are both plausible near the same, imprecisely
  // pinned zlin.  "Fix R, solve Z" is the well-conditioned direction there
  // instead, since Z is changing little across the bracket, so rlin is a
  // reliable value to hold fixed.  tok_logical_trace_sample (the final
  // per-grid-corner lookup) already makes exactly this choice locally, per
  // bracket; tok_trace_sample did not, always following the whole trace's
  // fixed param_is_r instead, which is what let a single bad bracket near
  // an X point produce a discrete branch flip.  Apply the same local check
  // here, matching that existing, already-relied-upon pattern, but only
  // for the closed-core path's own direction (param_is_r false): the
  // param_is_r true path (used by the different, unrelated Z-branch
  // routes) is left exactly as before.
  bool local_prefer_r_fixed = !param_is_r && fabs(tr[hi] - tr[lo]) >= fabs(tz[hi] - tz[lo]);
  // Both bracket endpoints lie on this psi contour, so the resampled point
  // belongs inside the bracket.  "Nearest root" does not enforce that: when
  // the intended root is missing -- an inboard SOL surface has a Z turning
  // point in R, so a fixed-R line can meet the contour again only far down
  // the divertor leg -- nearest-value returns that distant root instead of
  // failing, and the block silently acquires a point most of a metre away.
  // Reject anything further from the chord than a couple of bracket lengths
  // and let the next projection have it.
  double seglen = hypot(tr[hi] - tr[lo], tz[hi] - tz[lo]);
  double max_away = 2.0 * seglen;
  bool sampled = false;
  // A3: the chord point moved onto the surface along the chord's normal.
  // A chord point exactly on a bracket end IS that trace point.
  if (w <= 0.0) {
    *r = tr[lo];
    *z = tz[lo];
    sampled = true;
  } else if (w >= 1.0) {
    *r = tr[hi];
    *z = tz[hi];
    sampled = true;
  } else {
    double cr = 0.0, cz = 0.0;
    if (tok_chord_normal_solve(
          geo, psi, rlin, zlin, tr[hi] - tr[lo], tz[hi] - tz[lo], 0.5 * seglen, &cr, &cz
        )) {
      *r = cr;
      *z = cz;
      sampled = true;
    }
  }
  if (!sampled && (param_is_r || local_prefer_r_fixed)) {
    double roots[32] = {0.0};
    int nr = tok_geo_Z_psiR(geo, psi, rlin, 16, roots);
    if (nr > 0) {
      double zc = tok_nearest_value(zlin, roots, nr);
      if (fabs(zc - zlin) <= max_away) {
        *r = rlin;
        *z = zc;
        sampled = true;
      }
    }
  }
  if (!sampled) {
    double roots[16] = {0.0}, dRdZ[16] = {0.0};
    double dR[16] = {0.0}, dZ[16] = {0.0};
    int nr = gkyl_tok_geo_R_psiZ(geo, psi, zlin, 8, roots, dRdZ, dR, dZ);
    if (nr > 0) {
      double r_choice = tok_nearest_value(rlin, roots, nr);
      if (tok_rlin_ambiguous(rlin, roots, nr)) {
        double r_walk = 0.0;
        if (tok_trace_walk_root(geo, psi, tz[lo], tr[lo], tz[hi], tr[hi], zlin, &r_walk)) {
          r_choice = r_walk;
        }
      }
      if (fabs(r_choice - rlin) <= max_away) {
        *r = r_choice;
        *z = zlin;
        sampled = true;
      }
    }
  }
  if (!sampled) {
    // Neither axis-aligned direction has a usable root.  Project the chord
    // point onto the contour along its own normal, bounded by the bracket.
    double pr = 0.0, pz = 0.0;
    if (seglen > 0.0 && tok_psi_normal_project(
                          geo, psi, rlin, zlin, max_away, 1e-9 * fmax(1.0, fabs(psi)), &pr, &pz
                        )) {
      *r = pr;
      *z = pz;
      sampled = true;
    }
  }
  if (!sampled) {
    return false;
  }
  double residual = tok_eval_psi_rz_local(geo, *r, *z) - psi;
  bool sample_ok = isfinite(*r) && isfinite(*z) && isfinite(residual) &&
                   fabs(residual) <= 1e-9 * fmax(1.0, fabs(psi));
  return sample_ok;
}

static bool
tok_build_far_trace(const struct gkyl_tok_geo_grid_inp *inp, struct arc_length_ctx *arc_ctx)
{
  if (arc_ctx->far_trace_initialized) {
    return true;
  }
  if (!arc_ctx->xpt_ray_initialized) {
    fprintf(
      stderr, "TOK_ORDERED_MAP far trace failed ftype=%d reason=ray_not_initialized\n", inp->ftype
    );
    return false;
  }
  double rf = 0.0, zf = 0.0;
  if (!tok_fixed_edge_point(inp, arc_ctx->geo, arc_ctx->xpt_ray_psi0, &rf, &zf)) {
    fprintf(
      stderr, "TOK_ORDERED_MAP far trace failed ftype=%d reason=no_fixed_edge_point psi=%.17g\n",
      inp->ftype, arc_ctx->xpt_ray_psi0
    );
    return false;
  }
  int n = tok_reference_trace_nodes(inp, arc_ctx->sep_trace_capacity);
  double *rr = gkyl_malloc(sizeof(double[n]));
  double *zr = gkyl_malloc(sizeof(double[n]));
  double *rz = gkyl_malloc(sizeof(double[n]));
  double *zz = gkyl_malloc(sizeof(double[n]));
  double score_r = DBL_MAX, score_z = DBL_MAX;
  bool ok_r = tok_build_contour_candidate(
    arc_ctx->geo, arc_ctx->xpt_ray_psi0, true, rf, zf, arc_ctx->xpt_ray_r0, arc_ctx->xpt_ray_z0, n,
    rr, zr, &score_r, "far_boundary", false
  );
  bool ok_z = tok_build_contour_candidate(
    arc_ctx->geo, arc_ctx->xpt_ray_psi0, false, rf, zf, arc_ctx->xpt_ray_r0, arc_ctx->xpt_ray_z0, n,
    rz, zz, &score_z, "far_boundary", false
  );
  // Preserve both original parameterizations whenever either works.  Only a
  // trace for which both were rejected and whose fixed edge is a known
  // midplane turning point receives the endpoint-clustered R retry.
  if (!ok_r && !ok_z && tok_fixed_edge_is_midplane(inp->ftype)) {
    ok_r = tok_build_contour_candidate(
      arc_ctx->geo, arc_ctx->xpt_ray_psi0, true, rf, zf, arc_ctx->xpt_ray_r0, arc_ctx->xpt_ray_z0,
      n, rr, zr, &score_r, "far_boundary", true
    );
  }
  if (!ok_r && !ok_z) {
    fprintf(
      stderr,
      "TOK_ORDERED_MAP far trace failed ftype=%d psi=%.17g fixed=(%.17g,%.17g) ray=(%.17g,%.17g)\n",
      inp->ftype, arc_ctx->xpt_ray_psi0, rf, zf, arc_ctx->xpt_ray_r0, arc_ctx->xpt_ray_z0
    );
    gkyl_free(rr);
    gkyl_free(zr);
    gkyl_free(rz);
    gkyl_free(zz);
    return false;
  }
  bool use_r = ok_r && (!ok_z || score_r <= score_z);
  const double *src_r = use_r ? rr : rz;
  const double *src_z = use_r ? zr : zz;
  bool fixed_first = tok_sep_fixed_edge_is_first(inp->ftype);
  arc_ctx->far_trace_s[0] = 0.0;
  for (int i = 0; i < n; ++i) {
    int src = fixed_first ? i : n - 1 - i;
    arc_ctx->far_trace_r[i] = src_r[src];
    arc_ctx->far_trace_z[i] = src_z[src];
    if (i > 0) {
      arc_ctx->far_trace_s[i] =
        arc_ctx->far_trace_s[i - 1] + hypot(
                                        arc_ctx->far_trace_r[i] - arc_ctx->far_trace_r[i - 1],
                                        arc_ctx->far_trace_z[i] - arc_ctx->far_trace_z[i - 1]
                                      );
    }
  }
  arc_ctx->far_trace_n = n;
  arc_ctx->far_trace_param_is_r = use_r;
  arc_ctx->far_trace_initialized = isfinite(arc_ctx->far_trace_s[n - 1]) &&
                                   arc_ctx->far_trace_s[n - 1] > 0.0;
  gkyl_free(rr);
  gkyl_free(zr);
  gkyl_free(rz);
  gkyl_free(zz);
  return arc_ctx->far_trace_initialized;
}

// Number of theta samples in a map trace: 4 trace segments per theta cell
// (the historical ratio), at least 48 segments. The row rule must agree with
// tok_build_current_ordered_trace exactly, so it lives in one place.
static int
tok_ext_map_trace_request(const struct gkyl_tok_geo_grid_inp *inp)
{
  return GKYL_MAX2(49, 4 * inp->cgrid.cells[2] + 1);
}

static int
tok_ext_map_trace_nodes(const struct gkyl_tok_geo_grid_inp *inp, int capacity)
{
  return GKYL_MIN2(capacity, tok_ext_map_trace_request(inp));
}

// Rows of the theta correspondence: one per node row of the block.
static int
tok_ext_ladder_rows(const struct gkyl_tok_geo_grid_inp *inp)
{
  return GKYL_MAX2(1, inp->cgrid.cells[0]);
}

// Radial fraction of every ladder rung, taken THROUGH THE POSITION MAP.
//
// This is the correctness requirement the lookup depends on, and until
// 2026-08-27 it was only half-satisfied.  The rungs used to be marched at
// psi = psisep + span*(k/m), i.e. uniform in psi, while `m` was held to a
// multiple of the row count so that "every row lands on a rung".  That argument
// is an INDEX argument: it is true only when the rows are themselves uniform in
// psi.  They generally are not -- the position map places them, and on a
// nonuniform grid it places them very unevenly.  Measured over the block's own
// rows:
//
//     NSTX-U phase1        1.0x spread, rows land 0.02-0.08 cells from a rung
//     STEP-nonuniform x1   337x spread, rows land up to 0.50 cells from a rung
//
// A row that misses its rung is interpolated, and its neighbour is interpolated
// by a DIFFERENT amount; that differential is precisely the radial crossing the
// ladder exists to prevent.  Measured: turning the ladder on for
// STEP-nonuniform reversed 10 radial steps at x1 and 22 at x2 (the count grows
// with the row count, as an interpolation-phase error must, where a geometric
// defect would not).
//
// Sampling the rungs at the same computational coordinates the rows use, and
// mapping them with the same map, puts every row on a rung for ANY sub-rung
// count k.  The lookup's interpolation weight is then identically zero at every
// row, so the guarantee is structural instead of assumed.
//
// Returns the number of fractions written (m+1), or 0 if the map is unusable.
// Fractions come back monotonically increasing from the separatrix end, which
// is the order the march and the lookup both expect; the computational axis may
// run either way, so the sampling is reversed here when it runs inward.
static int
tok_ext_ladder_rung_fractions(
  const struct gkyl_tok_geo_grid_inp *inp, const struct arc_length_ctx *arc_ctx, int m, double *rf,
  int cap
)
{
  const struct gkyl_position_map *pm = arc_ctx->position_map;
  if (!pm || !pm->maps[0] || m < 1 || m + 1 > cap) {
    return 0;
  }
  const double psisep = arc_ctx->geo->psisep;
  const double span = arc_ctx->xpt_ray_psi0 - psisep;
  if (!isfinite(span) || span == 0.0) {
    return 0;
  }

  const double c_lo = inp->cgrid.lower[0], c_hi = inp->cgrid.upper[0];
  for (int k = 0; k <= m; ++k) {
    double pc = c_lo + (c_hi - c_lo) * (k / (double)m), psi_k = 0.0;
    pm->maps[0](0.0, &pc, &psi_k, pm->ctxs[0]);
    if (!isfinite(psi_k)) {
      return 0;
    }
    rf[k] = (psi_k - psisep) / span;
  }
  if (rf[0] > rf[m]) { // computational axis runs inward
    for (int a = 0, b = m; a < b; ++a, --b) {
      double t = rf[a];
      rf[a] = rf[b];
      rf[b] = t;
    }
  }
  // The ends ARE the block's radial boundaries; pin them so the lookup's
  // bracketing search can never fall outside the table on a rounding error.
  rf[0] = 0.0;
  rf[m] = 1.0;
  for (int k = 1; k <= m; ++k) {
    if (!(rf[k] > rf[k - 1])) { // must be strictly increasing
      return 0;
    }
  }
  return m + 1;
}

// Nearest point on a polyline, searched forward from *j_lo.  u increases
// monotonically along the ladder row, so its image must too; marching the
// window forward makes that true by construction and removes the backward jump
// between the two nearby legs of a single-null contour that a global argmin
// would allow.  Returns the normalized arc position.
static double
tok_ext_project_onto_trace(
  const double *tr, const double *tz, const double *ts, int n, double rp, double zp, int *j_lo
)
{
  double best_d2 = DBL_MAX, best_v = 0.0;
  int best_j = *j_lo;
  double total = ts[n - 1];
  for (int j = *j_lo; j < n - 1; ++j) {
    double r0 = tr[j], z0 = tz[j];
    double dr = tr[j + 1] - r0, dz = tz[j + 1] - z0;
    double den = dr * dr + dz * dz;
    double t = den > 0.0 ? ((rp - r0) * dr + (zp - z0) * dz) / den : 0.0;
    t = fmin(1.0, fmax(0.0, t));
    double d2 = SQ(rp - (r0 + t * dr)) + SQ(zp - (z0 + t * dz));
    if (d2 < best_d2) {
      best_d2 = d2;
      best_j = j;
      best_v = total > 0.0 ? (ts[j] + t * (ts[j + 1] - ts[j])) / total : 0.0;
    }
  }
  *j_lo = best_j;
  return fmin(1.0, fmax(0.0, best_v));
}

// ---------------------------------------------------------------------------
// The per-row correspondence with a fold-derived absorption of the X-point-ray
// endpoint offset (the row rule).
//
// MEASURED 2026-09-22/23 on the NSTX-U 450 and the STEP/ASDEX/TCV theta
// ladders (fable-handoff/01 §11-§16):
//   * every row off the separatrix starts on the X-point ray a tangential
//     distance t0 ~ sqrt(dpsi) from the separatrix row's endpoint, while the
//     row is only ~dpsi away radially (65x-127x apart on NSTX-U);
//   * uniform arc spreads that offset over the whole row as a shear that
//     folds the thin separatrix-adjacent cell where kappa*t*(L-t)/2 > gap;
//   * the marched ladder projects interior nodes and pins the boundary node,
//     so the FIRST CELL absorbs t0: no folds, but the psi-shear metric reads
//     L/(L - t0) and diverges as theta is refined -- the seam divergence of
//     01 §5 -- and once the offset exceeds a cell a min-gap floor takes over;
//   * evaluated offline with the exact chord-crossing test, the longest
//     decay arc that is fold-free is the WHOLE ROW for ~90-100% of blocks on
//     all four devices, and the resulting shear is near 1 and flat in theta.
//
// So: one rung per node row (rungs ARE rows: no refinement loop, no trace-
// capacity cap), marched from the separatrix outward; each rung projects the
// previous rung's nodes (t = 0), then absorbs each end's offset delta over a
// decay arc lambda*total,
//
//     w[i] = w_proj[i] + delta0*max(0, 1 - s_i/lambda) + delta1*max(0, 1 - (1-s_i)/lambda)
//
// with lambda the LARGEST value in [1/(n-1), 1] for which no cell between the
// two rungs has crossing corner chords -- the same segment-intersection test
// grid_gate's CELL_FOLD applies to the written grid.  No constant: the
// endpoints, the trace, the previous rung and the crossing test are the only
// inputs.  If even lambda = 1/(n-1) (one-cell absorption, the old ladder's
// behaviour) crosses, the offset has moved by more than a cell and no
// absorption is fold-free: that is the jump regime (01 §15, a curled divertor
// leg) and the block REFUSES, loudly, with the offset in metres.
//
// Monotonicity: projection onto a curled trace can reorder nodes (the hop the
// old comment warns about); a pool-adjacent-violators pass restores order
// with no constant, and a rung that is still not strictly increasing refuses.
// No min-gap floor, no identity blend, no relabel: with the separatrix as the
// seed its row carries the identity (or the |grad psi| map when enabled), so
// the seam invariant holds by construction.
//
// The table lands in the same ext_ladder_* slots the marched ladder fills, so
// tok_ext_ladder_w and every consumer are untouched; node rows bracket exactly
// (rungs are rows), only intra-cell quadrature points interpolate.
static bool
tok_ext_seg_cross(
  double ax, double ay, double bx, double by, double cx, double cy, double dx, double dy
)
{
  double d1x = bx - ax, d1y = by - ay, d2x = dx - cx, d2y = dy - cy;
  double den = d1x * d2y - d1y * d2x;
  if (den == 0.0) {
    return false;
  }
  double sx = cx - ax, sy = cy - ay;
  double t = (sx * d2y - sy * d2x) / den, u = (sx * d1y - sy * d1x) / den;
  return t > 0.0 && t < 1.0 && u > 0.0 && u < 1.0;
}

// Fill w[] for a decay fraction lambda: blend between the PROJECTED positions
// (zero tangential shift, the ladder's correspondence) and the PREVIOUS rung's
// normalised positions (uniform arc, which carries the endpoint offset along
// the whole row) with a weight h that is 1 at both ends and decays to 0 over
// lambda of the row from each end:
//
//     w[i] = w_proj[i] + (w_prev[i] - w_proj[i]) * h(s_i),
//     h(s) = min(1, max(0, 1 - s/lambda) + max(0, 1 - (1-s)/lambda)).
//
// lambda = 1 gives w = w_prev exactly (uniform arc, monotone by construction);
// lambda = 1/(n-1) leaves every interior node projected, i.e. the whole
// endpoint offset in the end cells -- the old ladder's behaviour.  Then sample
// the candidate nodes on the current trace and test every cell against the
// previous rung.  Returns 0 = clean, 1 = a cell folds, 2 = not strictly
// increasing / sample failed.
// The crossing test runs on the CELL CORNERS (ncell+1 per row, ncell = the
// block's theta cell count), not on the n map-trace nodes: the trace is ~4x
// denser than the cells and the chord criterion scales with the chord length
// (measured 2026-09-23: at trace scale no crossing was found and every
// uniform-arc fold of the written grid survived).  Corner j sits at
// u = j/ncell, i.e. between trace nodes; its w is interpolated linearly, and
// the previous rung's corner is read the same way from wprev.
static int
tok_ext_row_rule_try(
  const struct gkyl_tok_geo *geo, double psi, const double *cr, const double *cz, const double *cs,
  int cn, bool cparam, const double *wproj, const double *wprev, int n, int ncell, double lambda,
  const double *ppr, const double *ppz, const double *pps, int ppn, bool ppparam, double ppsi,
  double *w, double *qr, double *qz
)
{
  for (int i = 0; i < n; ++i) {
    double si = fmin(1.0, fmax(0.0, wprev[i]));
    double h = fmax(0.0, 1.0 - si / lambda) + fmax(0.0, 1.0 - (1.0 - si) / lambda);
    h = fmin(1.0, h);
    w[i] = wproj[i] + (wprev[i] - wproj[i]) * h;
  }
  w[0] = 0.0;
  w[n - 1] = 1.0;
  for (int i = 1; i < n; ++i) {
    if (!(w[i] > w[i - 1])) {
      return 2;
    }
  }
  // Every map-trace node must land on a DISTINCT point of the trace: two
  // nodes projected onto one place (a contour that shortens, 204951's plate
  // strike) sample to the same R,Z and the ordered-map builder aborts on a
  // zero arc step downstream.  Treated as unsafe here, so the search moves to
  // a decay fraction that separates them, instead of a floor constant.
  {
    double xr_prev = 0.0, xz_prev = 0.0;
    for (int i = 0; i < n; ++i) {
      double xr = 0.0, xz = 0.0;
      if (!tok_trace_sample(geo, psi, cr, cz, cs, cn, cparam, w[i], &xr, &xz)) {
        return 2;
      }
      if (i > 0 && !(hypot(xr - xr_prev, xz - xz_prev) > 0.0)) {
        return 2;
      }
      xr_prev = xr;
      xz_prev = xz;
    }
  }
  // cell corners of this rung (qr,qz) and of the previous rung (pr_c,pz_c),
  // plus the two Gauss points inside every cell on both rows (2026-09-30).
  //
  // Corners alone test the CHORD of each theta edge, and a chord is not the
  // edge: on the end cells next to the X-point ray the rows bend within one
  // cell, the corner chords stay apart while the rows themselves cross, and
  // the library's signed Jacobian at the Gauss nodes then aborts the build
  // (NSTX-U 203970, 203834, 204056: PF_LO_L, last Gauss row, theta 0). The
  // rule must judge the edge the geometry will actually carry, so each edge is
  // the polyline corner -> Gauss -> Gauss -> corner, sampled on the row's own
  // trace at the theta positions the quadrature will use. Same instrument as
  // the harness gate since surface nodes are always written; same definition
  // of a fold: the two theta edges cross, or the two radial edges cross. A
  // theta edge meeting its own radial edge next to their shared corner (a
  // cusp: the radial line tangent to the row) is not a fold here either; the
  // Jacobian guard at the Gauss nodes decides those, as it does for the gate.
  int nc = ncell > 0 ? ncell : n - 1;
  if (nc + 1 > n) {
    nc = n - 1;
  }
  double pr_c[nc + 1], pz_c[nc + 1];
  double qg_r[2 * nc], qg_z[2 * nc], pg_r[2 * nc], pg_z[2 * nc];
  const double gauss[2] = {0.5 - 0.5 / sqrt(3.0), 0.5 + 0.5 / sqrt(3.0)};
  for (int j = 0; j <= nc; ++j) {
    for (int g = -1; g < 2; ++g) {
      if (g >= 0 && j == nc) {
        break;
      }
      double x = (j + (g < 0 ? 0.0 : gauss[g])) / nc * (n - 1);
      int i = (int)floor(x);
      if (i > n - 2) {
        i = n - 2;
      }
      double f = x - i;
      double wc = w[i] + f * (w[i + 1] - w[i]);
      double wp = wprev[i] + f * (wprev[i + 1] - wprev[i]);
      double *xr = g < 0 ? &qr[j] : &qg_r[2 * j + g], *xz = g < 0 ? &qz[j] : &qg_z[2 * j + g];
      double *yr = g < 0 ? &pr_c[j] : &pg_r[2 * j + g], *yz = g < 0 ? &pz_c[j] : &pg_z[2 * j + g];
      if (!tok_trace_sample(geo, psi, cr, cz, cs, cn, cparam, wc, xr, xz)) {
        return 2;
      }
      if (!tok_trace_sample(geo, ppsi, ppr, ppz, pps, ppn, ppparam, wp, yr, yz)) {
        return 2;
      }
    }
  }
  for (int j = 0; j < nc; ++j) {
    // the two theta edges of cell j as 3-segment polylines
    const double cur_r[4] = {qr[j], qg_r[2 * j], qg_r[2 * j + 1], qr[j + 1]},
                 cur_z[4] = {qz[j], qg_z[2 * j], qg_z[2 * j + 1], qz[j + 1]};
    const double prv_r[4] = {pr_c[j], pg_r[2 * j], pg_r[2 * j + 1], pr_c[j + 1]},
                 prv_z[4] = {pz_c[j], pg_z[2 * j], pg_z[2 * j + 1], pz_c[j + 1]};
    for (int a = 0; a < 3; ++a) {
      for (int b = 0; b < 3; ++b) {
        if (tok_ext_seg_cross(
              prv_r[a], prv_z[a], prv_r[a + 1], prv_z[a + 1], cur_r[b], cur_z[b], cur_r[b + 1],
              cur_z[b + 1]
            )) {
          return 1;
        }
      }
    }
    // the two radial edges: prev j+1 -> cur j+1 against cur j -> prev j
    if (tok_ext_seg_cross(
          pr_c[j + 1], pz_c[j + 1], qr[j + 1], qz[j + 1], qr[j], qz[j], pr_c[j], pz_c[j]
        )) {
      return 1;
    }
  }
  return 0;
}

static bool
tok_ext_build_theta_rows(const struct gkyl_tok_geo_grid_inp *inp, struct arc_length_ctx *arc_ctx)
{
  const int cap = arc_ctx->sep_trace_capacity;
  const int n = tok_ext_map_trace_nodes(inp, cap);
  const int m = tok_ext_ladder_rows(inp);
  // theta cells: the corners the gate tests.  cgrid is (psi, phi, theta):
  // cells[0] radial, cells[2] theta (tok_ext_map_trace_nodes sizes the map
  // trace from cells[2]); cells[1] is the one toroidal cell -- using it made
  // the corner test see two points and no crossing (2026-09-23).
  const int ncell = GKYL_MAX2(1, inp->cgrid.cells[2]);
  const double psisep = arc_ctx->geo->psisep;
  const double span = arc_ctx->xpt_ray_psi0 - psisep;
  if (n < 3 || m < 1 || !isfinite(span) || span == 0.0) {
    return false;
  }
  double *rf = gkyl_malloc(sizeof(double[m + 1]));
  if (tok_ext_ladder_rung_fractions(inp, arc_ctx, m, rf, m + 1) != m + 1) {
    fprintf(
      stderr, "TOK_EXT_ROW_RULE REFUSED ftype=%d reason=rung_fractions_unavailable rows=%d\n",
      inp->ftype, m
    );
    tok_wall_trial_note_row_rule_refused();
    gkyl_free(rf);
    return false;
  }
  double *table = gkyl_malloc(sizeof(double[(size_t)(m + 1) * n]));
  double *scratch = gkyl_malloc(sizeof(double[6 * cap]));
  double *pr = scratch, *pz = scratch + cap, *ps = scratch + 2 * cap;
  double *cr = scratch + 3 * cap, *cz = scratch + 4 * cap, *cs = scratch + 5 * cap;
  double *wprev = gkyl_malloc(sizeof(double[n])), *wproj = gkyl_malloc(sizeof(double[n]));
  double *wcur = gkyl_malloc(sizeof(double[n])), *wtry = gkyl_malloc(sizeof(double[n]));
  double *pr_n = gkyl_malloc(sizeof(double[n])), *pz_n = gkyl_malloc(sizeof(double[n]));
  double *qr = gkyl_malloc(sizeof(double[n])), *qz = gkyl_malloc(sizeof(double[n]));
  double *bm = gkyl_malloc(sizeof(double[n]));
  int *bc = gkyl_malloc(sizeof(int[n]));

  // Seed: the separatrix row, uniform in arc -- exactly what the un-laddered
  // path gives that row, so the shared row is unchanged.
  for (int i = 0; i < n; ++i) {
    wprev[i] = i / (double)(n - 1);
    table[i] = wprev[i];
  }
  int pn = arc_ctx->sep_trace_n;
  bool pparam = arc_ctx->sep_trace_param_is_r;
  double ppsi = psisep;
  for (int i = 0; i < pn; ++i) {
    pr[i] = arc_ctx->sep_trace_r[i];
    pz[i] = arc_ctx->sep_trace_z[i];
    ps[i] = arc_ctx->sep_trace_s[i];
  }

  bool ok = true;
  int n_short = 0;
  double lam_min = 1.0, off_max = 0.0;
  const double lam_floor = 1.0 / (n - 1);
  for (int k = 1; k <= m && ok; ++k) {
    double psi_k = psisep + span * rf[k];
    int cn = 0;
    bool cparam = false, cclosed = false;
    if (!tok_ext_build_domain_trace(
          inp, arc_ctx, psi_k, false, cr, cz, cs, &cn, &cparam, &cclosed
        ) ||
        cn < 2) {
      fprintf(
        stderr, "TOK_EXT_ROW_RULE REFUSED ftype=%d reason=row_trace_failed row=%d psi=%.17g\n",
        inp->ftype, k, psi_k
      );
      tok_wall_trial_note_row_rule_refused();
      ok = false;
      break;
    }
    // previous rung's nodes in R,Z, and their projection onto this trace
    int j_lo = 0;
    for (int i = 0; i < n; ++i) {
      if (!tok_trace_sample(
            arc_ctx->geo, ppsi, pr, pz, ps, pn, pparam, wprev[i], &pr_n[i], &pz_n[i]
          )) {
        fprintf(
          stderr, "TOK_EXT_ROW_RULE REFUSED ftype=%d reason=prev_sample_failed row=%d i=%d\n",
          inp->ftype, k, i
        );
        tok_wall_trial_note_row_rule_refused();
        ok = false;
        break;
      }
      wproj[i] = tok_ext_project_onto_trace(cr, cz, cs, cn, pr_n[i], pz_n[i], &j_lo);
    }
    if (!ok) {
      break;
    }
    // pool-adjacent-violators on the interior: restore order without a constant
    {
      int nb = 0;
      for (int i = 1; i < n - 1; ++i) {
        bm[nb] = wproj[i];
        bc[nb] = 1;
        ++nb;
        while (nb > 1 && bm[nb - 2] > bm[nb - 1]) {
          int c = bc[nb - 2] + bc[nb - 1];
          bm[nb - 2] = (bc[nb - 2] * bm[nb - 2] + bc[nb - 1] * bm[nb - 1]) / c;
          bc[nb - 2] = c;
          --nb;
        }
      }
      int out = 1;
      for (int b = 0; b < nb; ++b) {
        for (int j = 0; j < bc[b]; ++j) {
          wproj[out++] = bm[b];
        }
      }
    }
    // The offset this rung must absorb is the difference between uniform arc
    // and projection next to each end, in metres of this row (reported).
    double off_m = fmax(fabs(wprev[1] - wproj[1]), fabs(wprev[n - 2] - wproj[n - 2])) * cs[cn - 1];
    if (off_m > off_max) {
      off_max = off_m;
    }
    // largest lambda in [lam_floor, 1] that is fold-free and monotone
    int r1 = tok_ext_row_rule_try(
      arc_ctx->geo, psi_k, cr, cz, cs, cn, cparam, wproj, wprev, n, ncell, 1.0, pr, pz, ps, pn,
      pparam, ppsi, wcur, qr, qz
    );
    double lam = 1.0;
    if (r1 != 0) {
      double lo = -1.0, hi = 1.0;
      int rf0 = tok_ext_row_rule_try(
        arc_ctx->geo, psi_k, cr, cz, cs, cn, cparam, wproj, wprev, n, ncell, lam_floor, pr, pz, ps,
        pn, pparam, ppsi, wcur, qr, qz
      );
      if (rf0 == 0) {
        lo = lam_floor;
      } else {
        // a coarse scan between the floor and 1 catches a safe window neither end has
        for (int q = 1; q < 32 && lo < 0.0; ++q) {
          double lq = lam_floor + (1.0 - lam_floor) * q / 32.0;
          if (tok_ext_row_rule_try(
                arc_ctx->geo, psi_k, cr, cz, cs, cn, cparam, wproj, wprev, n, ncell, lq, pr, pz, ps,
                pn, pparam, ppsi, wcur, qr, qz
              ) == 0) {
            lo = lq;
          }
        }
      }
      if (lo < 0.0) {
        fprintf(
          stderr,
          "TOK_EXT_ROW_RULE REFUSED ftype=%d reason=%s row=%d psi=%.17g "
          "endpoint_offset_m=%.6g end_cell_m=%.6g\n",
          inp->ftype, (r1 == 1 || rf0 == 1) ? "no_fold_free_absorption" : "row_not_monotone", k,
          psi_k, off_m, cs[cn - 1] / (n - 1)
        );
        tok_wall_trial_note_row_rule_refused();
        ok = false;
        break;
      }
      for (int it = 0; it < 40 && hi - lo > 1e-6; ++it) {
        double mid = 0.5 * (lo + hi);
        if (tok_ext_row_rule_try(
              arc_ctx->geo, psi_k, cr, cz, cs, cn, cparam, wproj, wprev, n, ncell, mid, pr, pz, ps,
              pn, pparam, ppsi, wtry, qr, qz
            ) == 0) {
          lo = mid;
        } else {
          hi = mid;
        }
      }
      lam = lo;
      tok_ext_row_rule_try(
        arc_ctx->geo, psi_k, cr, cz, cs, cn, cparam, wproj, wprev, n, ncell, lam, pr, pz, ps, pn,
        pparam, ppsi, wcur, qr, qz
      );
      ++n_short;
      if (lam < lam_min) {
        lam_min = lam;
      }
    }
    for (int i = 0; i < n; ++i) {
      table[(size_t)k * n + i] = wcur[i];
    }
    // this rung becomes the next one's reference
    for (int i = 0; i < cn; ++i) {
      pr[i] = cr[i];
      pz[i] = cz[i];
      ps[i] = cs[i];
    }
    pn = cn;
    pparam = cparam;
    ppsi = psi_k;
    for (int i = 0; i < n; ++i) {
      wprev[i] = wcur[i];
    }
  }
  gkyl_free(scratch);
  gkyl_free(wprev);
  gkyl_free(wproj);
  gkyl_free(wcur);
  gkyl_free(wtry);
  gkyl_free(pr_n);
  gkyl_free(pz_n);
  gkyl_free(qr);
  gkyl_free(qz);
  gkyl_free(bm);
  gkyl_free(bc);
  if (!ok) {
    gkyl_free(table);
    gkyl_free(rf);
    return false;
  }
  gkyl_free(arc_ctx->ext_ladder_w);
  gkyl_free(arc_ctx->ext_ladder_rf);
  arc_ctx->ext_ladder_w = table;
  arc_ctx->ext_ladder_rf = rf;
  arc_ctx->ext_ladder_m = m;
  arc_ctx->ext_ladder_n = n;
  arc_ctx->ext_ladder_initialized = true;
  // Always reported: one line per block, so a run can be audited for how much
  // absorption the rule actually needed.
  fprintf(
    stderr,
    "TOK_EXT_ROW_RULE built ftype=%d rows=%d nodes=%d ncell=%d rows_shortened=%d lambda_min=%.4g max_endpoint_offset_m=%.6g\n",
    inp->ftype, m, n, ncell, n_short, lam_min, off_max
  );
  return true;
}

// w(u_i, psi) by linear interpolation between the two bracketing rungs.  A
// convex combination of two increasing rows is increasing, so the interpolated
// row is monotone in THETA for every radial fraction.
//
// RADIAL ordering is the other axis, and it is handled by WHERE the rungs are
// rather than by this interpolation: ext_ladder_rf holds the block's own node
// rows (tok_ext_ladder_rung_fractions), so a node row brackets against itself,
// t comes out 0, and no row is perturbed relative to its neighbour.  Only the
// intra-cell quadrature points interpolate, and they do not take part in the
// radial ordering of the corner nodes.
//
// The search is a bisection rather than an index computation because the rungs
// are not uniformly spaced -- that is the whole point of them.
static double
tok_ext_ladder_w(const struct arc_length_ctx *arc_ctx, double rf, int i)
{
  int m = arc_ctx->ext_ladder_m, n = arc_ctx->ext_ladder_n;
  const double *rfa = arc_ctx->ext_ladder_rf;
  double x = fmin(1.0, fmax(0.0, rf));
  int lo = 0, hi = m;
  while (hi - lo > 1) {
    int mid = (lo + hi) / 2;
    if (rfa[mid] <= x) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  double d = rfa[hi] - rfa[lo];
  double t = d > 0.0 ? (x - rfa[lo]) / d : 0.0;
  t = fmin(1.0, fmax(0.0, t));
  const double *w0 = arc_ctx->ext_ladder_w + (size_t)lo * n;
  const double *w1 = arc_ctx->ext_ladder_w + (size_t)hi * n;
  return (1.0 - t) * w0[i] + t * w1[i];
}

static bool
tok_build_trace_correspondence(
  const struct gkyl_tok_geo_grid_inp *inp, struct arc_length_ctx *arc_ctx
)
{
  const bool extended = tok_ext_construction(inp);
  if (arc_ctx->ordered_boundaries_initialized) {
    return true;
  }
  bool closed = false;
  if (extended) {
    bool sep_closed = false, far_closed = false;
    bool sep_built_generic = false, sep_built_here = false;
    if (!arc_ctx->sep_trace_initialized) {
      if (!tok_ext_build_domain_trace(
            inp, arc_ctx, arc_ctx->geo->psisep, true, arc_ctx->sep_trace_r, arc_ctx->sep_trace_z,
            arc_ctx->sep_trace_s, &arc_ctx->sep_trace_n, &arc_ctx->sep_trace_param_is_r, &sep_closed
          )) {
        fprintf(
          stderr, "TOK_TRACE_CORR reason=sep_domain_trace_failed ftype=%d psisep=%.17g\n",
          inp->ftype, arc_ctx->geo->psisep
        );
        return false;
      }
      arc_ctx->sep_trace_initialized = true;
      sep_built_generic = arc_ctx->ext_last_trace_used_generic;
      sep_built_here = true;
    } else {
      struct tok_ext_topology top;
      if (!tok_ext_topology_from_ftype(inp->ftype, inp->half_domain, &top)) {
        fprintf(stderr, "TOK_TRACE_CORR reason=sep_topology_unknown ftype=%d\n", inp->ftype);
        return false;
      }
      sep_closed = top.closed;
    }
    if (!arc_ctx->far_trace_initialized) {
      if (!tok_ext_build_domain_trace(
            inp, arc_ctx, arc_ctx->xpt_ray_psi0, false, arc_ctx->far_trace_r, arc_ctx->far_trace_z,
            arc_ctx->far_trace_s, &arc_ctx->far_trace_n, &arc_ctx->far_trace_param_is_r, &far_closed
          )) {
        fprintf(
          stderr, "TOK_TRACE_CORR reason=far_domain_trace_failed ftype=%d ray_psi0=%.17g\n",
          inp->ftype, arc_ctx->xpt_ray_psi0
        );
        return false;
      }
      arc_ctx->far_trace_initialized = true;
      // The far surface can be the one that trips the block onto the generic
      // route, in which case the separatrix trace above was already built with
      // the topology route.  Rebuild it so both radial boundaries -- and every
      // interior surface, which now follows the flag -- share one route.
      if (sep_built_here && !sep_built_generic && arc_ctx->ext_force_generic_route && !sep_closed) {
        if (!tok_ext_build_domain_trace(
              inp, arc_ctx, arc_ctx->geo->psisep, true, arc_ctx->sep_trace_r, arc_ctx->sep_trace_z,
              arc_ctx->sep_trace_s, &arc_ctx->sep_trace_n, &arc_ctx->sep_trace_param_is_r,
              &sep_closed
            )) {
          fprintf(
            stderr, "TOK_TRACE_CORR reason=sep_regeneric_failed ftype=%d psisep=%.17g\n",
            inp->ftype, arc_ctx->geo->psisep
          );
          return false;
        }
        fprintf(
          stderr, "TOK_EXT_TRACE sep_rebuilt_generic ftype=%d psisep=%.17g\n", inp->ftype,
          arc_ctx->geo->psisep
        );
      }
    } else {
      struct tok_ext_topology top;
      if (!tok_ext_topology_from_ftype(inp->ftype, inp->half_domain, &top)) {
        fprintf(stderr, "TOK_TRACE_CORR reason=far_topology_unknown ftype=%d\n", inp->ftype);
        return false;
      }
      far_closed = top.closed;
    }
    if (sep_closed != far_closed) {
      fprintf(
        stderr, "TOK_TRACE_CORR reason=closedness_mismatch ftype=%d sep_closed=%d far_closed=%d\n",
        inp->ftype, (int)sep_closed, (int)far_closed
      );
      return false;
    }
    closed = sep_closed;
  } else {
    if (!arc_ctx->sep_trace_initialized) {
      double rf = 0.0, zf = 0.0;
      if (!tok_fixed_edge_point(inp, arc_ctx->geo, arc_ctx->geo->psisep, &rf, &zf)) {
        fprintf(
          stderr, "TOK_TRACE_CORR reason=fixed_edge_point_failed ftype=%d psisep=%.17g\n",
          inp->ftype, arc_ctx->geo->psisep
        );
        return false;
      }
      tok_build_sep_trace(inp, arc_ctx, zf, rf);
    }
    if (!tok_build_far_trace(inp, arc_ctx)) {
      fprintf(
        stderr, "TOK_TRACE_CORR reason=far_trace_failed ftype=%d ray_psi0=%.17g\n", inp->ftype,
        arc_ctx->xpt_ray_psi0
      );
      return false;
    }
  }

  int n = GKYL_MIN2(arc_ctx->sep_trace_n, arc_ctx->far_trace_n);
  if (n < 2) {
    fprintf(
      stderr, "TOK_TRACE_CORR reason=too_few_trace_points ftype=%d n=%d sep_n=%d far_n=%d\n",
      inp->ftype, n, arc_ctx->sep_trace_n, arc_ctx->far_trace_n
    );
    return false;
  }
  // The extended full-domain traces already share topology-matched endpoints
  // and direction.  Equal normalized contour length therefore supplies a
  // one-to-one coordinate on every radial surface.  In particular, it avoids
  // radial folds caused when nearest-point projection jumps between the two
  // nearby legs of a single-null contour.  A closed core also needs this
  // identity map because its duplicated seam makes nearest projection
  // ambiguous.  Keep the original projection/PAVA correspondence unchanged
  // for the established half-domain path.
  const bool ext_open = extended && !closed;
  const bool ext_ladder = ext_open;
  if (extended || closed) {
    for (int i = 0; i < n; ++i) {
      arc_ctx->trace_corr_v[i] = i / (double)(n - 1);
    }
    arc_ctx->trace_corr_n = n;
    // Both radial boundary traces exist and the block's route has settled, so
    // this is the one point where the theta rows can be built without
    // perturbing route selection.
    if (ext_ladder && !arc_ctx->ext_ladder_initialized && !arc_ctx->ext_ladder_failed) {
      if (!tok_ext_build_theta_rows(inp, arc_ctx)) {
        // A refused row rule used to fall back to the two-point blend and the
        // case PASSED silently (1-5 refusals per rung on the 450).  User
        // decision 2026-10-02: a refusal is fatal in a production build, so
        // the next baseline shows every case the rule cannot handle.  Inside
        // a wall trial it stays a recorded diagnostic (the trial is built
        // plain and judged on containment), as before.
        if (!tok_wall_trial_is_active()) {
          fprintf(stderr, "TOK_EXT_LADDER reason=row_rule_refused_fatal ftype=%d\n", inp->ftype);
          return false;
        }
        arc_ctx->ext_ladder_failed = true;
        fprintf(
          stderr, "TOK_EXT_LADDER reason=build_failed ftype=%d falling back to blend\n", inp->ftype
        );
      }
    }
    arc_ctx->ordered_boundaries_initialized = true;
    return true;
  }
  arc_ctx->trace_corr_v[0] = 0.0;
  for (int i = 1; i < n - 1; ++i) {
    double u = i / (double)(n - 1), rs = 0.0, zs = 0.0;
    if (!tok_trace_sample(
          arc_ctx->geo, arc_ctx->geo->psisep, arc_ctx->sep_trace_r, arc_ctx->sep_trace_z,
          arc_ctx->sep_trace_s, arc_ctx->sep_trace_n, arc_ctx->sep_trace_param_is_r, u, &rs, &zs
        )) {
      fprintf(
        stderr,
        "TOK_TRACE_CORR reason=sep_trace_sample_failed ftype=%d i=%d n=%d u=%.17g sep_n=%d\n",
        inp->ftype, i, n, u, arc_ctx->sep_trace_n
      );
      return false;
    }
    double best_d2 = DBL_MAX, best_v = 0.0;
    for (int j = 0; j < arc_ctx->far_trace_n - 1; ++j) {
      double r0 = arc_ctx->far_trace_r[j], z0 = arc_ctx->far_trace_z[j];
      double dr = arc_ctx->far_trace_r[j + 1] - r0;
      double dz = arc_ctx->far_trace_z[j + 1] - z0;
      double den = dr * dr + dz * dz;
      double t = den > 0.0 ? ((rs - r0) * dr + (zs - z0) * dz) / den : 0.0;
      t = fmin(1.0, fmax(0.0, t));
      double rp = r0 + t * dr, zp = z0 + t * dz;
      double d2 = SQ(rs - rp) + SQ(zs - zp);
      if (d2 < best_d2) {
        best_d2 = d2;
        best_v =
          (arc_ctx->far_trace_s[j] + t * (arc_ctx->far_trace_s[j + 1] - arc_ctx->far_trace_s[j])) /
          arc_ctx->far_trace_s[arc_ctx->far_trace_n - 1];
      }
    }
    arc_ctx->trace_corr_v[i] = best_v;
  }
  arc_ctx->trace_corr_v[n - 1] = 1.0;

  // Isotonic regression (increasing PAVA) finds the least-squares monotone
  // correspondence, rather than greedily flattening every value after the
  // first local reversal.  Endpoints are already bounded by every projection,
  // so applying PAVA to the interior preserves the exact 0 and 1 anchors.
  int nblock = 0;
  double *block_mean = gkyl_malloc(sizeof(double[n]));
  int *block_count = gkyl_malloc(sizeof(int[n]));
  for (int i = 1; i < n - 1; ++i) {
    block_mean[nblock] = arc_ctx->trace_corr_v[i];
    block_count[nblock] = 1;
    ++nblock;
    while (nblock > 1 && block_mean[nblock - 2] > block_mean[nblock - 1]) {
      int count = block_count[nblock - 2] + block_count[nblock - 1];
      block_mean[nblock - 2] = (block_count[nblock - 2] * block_mean[nblock - 2] +
                                block_count[nblock - 1] * block_mean[nblock - 1]) /
                               count;
      block_count[nblock - 2] = count;
      --nblock;
    }
  }
  int out = 1;
  for (int b = 0; b < nblock; ++b) {
    for (int j = 0; j < block_count[b]; ++j) {
      arc_ctx->trace_corr_v[out++] = block_mean[b];
    }
  }
  gkyl_free(block_mean);
  gkyl_free(block_count);
  // Nearest-point projection can legitimately collapse a finite interval of
  // the separatrix onto one far-boundary location.  Blend in a small identity
  // map so g remains well-conditioned and strictly order preserving without
  // discarding the nearest-normal correspondence.  Shot 203730 remains
  // fold-free over a broad range around this conservative one-percent blend.
  const double identity_fraction = tok_trace_corr_identity_fraction();
  for (int i = 1; i < n - 1; ++i) {
    double u = i / (double)(n - 1);
    arc_ctx->trace_corr_v[i] =
      identity_fraction * u + (1.0 - identity_fraction) * arc_ctx->trace_corr_v[i];
  }
  for (int i = 1; i < n; ++i) {
    if (!(arc_ctx->trace_corr_v[i] > arc_ctx->trace_corr_v[i - 1])) {
      fprintf(
        stderr,
        "TOK_TRACE_CORR reason=nonmonotonic_correspondence ftype=%d i=%d n=%d v_prev=%.17g v_curr=%.17g\n",
        inp->ftype, i, n, arc_ctx->trace_corr_v[i - 1], arc_ctx->trace_corr_v[i]
      );
      return false;
    }
  }
  arc_ctx->trace_corr_n = n;
  arc_ctx->ordered_boundaries_initialized = true;
  return true;
}

static double
tok_trace_correspondence(const struct arc_length_ctx *arc_ctx, double u)
{
  const double endpoint_tol = 256.0 * DBL_EPSILON;
  if (u <= endpoint_tol) {
    return 0.0;
  }
  if (u >= 1.0 - endpoint_tol) {
    return 1.0;
  }
  double x = u * (arc_ctx->trace_corr_n - 1);
  int i = GKYL_MIN2(arc_ctx->trace_corr_n - 2, (int)floor(x));
  double w = x - i;
  return arc_ctx->trace_corr_v[i] + w * (arc_ctx->trace_corr_v[i + 1] - arc_ctx->trace_corr_v[i]);
}

static bool
tok_ordered_chord_point(
  const struct gkyl_tok_geo_grid_inp *inp, struct arc_length_ctx *arc_ctx, double psi, double u,
  double *r, double *z
)
{
  if (!tok_build_trace_correspondence(inp, arc_ctx)) {
    return false;
  }
  u = fmin(1.0, fmax(0.0, u));
  double v = tok_trace_correspondence(arc_ctx, u);
  double ra = 0.0, za = 0.0, rb = 0.0, zb = 0.0;
  if (!tok_trace_sample(
        arc_ctx->geo, arc_ctx->geo->psisep, arc_ctx->sep_trace_r, arc_ctx->sep_trace_z,
        arc_ctx->sep_trace_s, arc_ctx->sep_trace_n, arc_ctx->sep_trace_param_is_r, u, &ra, &za
      ) ||
      !tok_trace_sample(
        arc_ctx->geo, arc_ctx->xpt_ray_psi0, arc_ctx->far_trace_r, arc_ctx->far_trace_z,
        arc_ctx->far_trace_s, arc_ctx->far_trace_n, arc_ctx->far_trace_param_is_r, v, &rb, &zb
      )) {
    return false;
  }
  double scale = fmax(1.0, fmax(fabs(psi), fabs(arc_ctx->geo->psisep)));
  if (tok_geo_same_flux(psi, arc_ctx->geo->psisep)) {
    *r = ra;
    *z = za;
    return true;
  }
  if (tok_geo_same_flux(psi, arc_ctx->xpt_ray_psi0)) {
    *r = rb;
    *z = zb;
    return true;
  }
  double fa = tok_eval_psi_rz_local(arc_ctx->geo, ra, za) - psi;
  double fb = tok_eval_psi_rz_local(arc_ctx->geo, rb, zb) - psi;
  if (!isfinite(fa) || !isfinite(fb) || fa * fb > 0.0) {
    fprintf(
      stderr, "TOK_ORDERED_MAP unbracketed chord ftype=%d psi=%.17g u=%.17g fa=%.17g fb=%.17g\n",
      inp->ftype, psi, u, fa, fb
    );
    return false;
  }
  // NOTE: psi is not always monotone along the chord -- where the chord leaves
  // the X point nearly tangent to a separatrix branch it can rise past psi, dip
  // back below and rise again, so psi=psi_curr has three roots and which one
  // this bisection converges to depends on the bracketing accident rather than
  // on a stated rule. Selecting the FIRST root instead was tried and is WRONG
  // here: on 204502's DN_SOL_IN_LO it puts the seam node of the first flux
  // surface at 0.031 m from the X point while the rest of that radial row sits
  // at ~0.15 m, because the near-X-point root lies on the piece of the contour
  // that wraps the saddle rather than on the arc that continues to the plate.
  // The outermost root, which this bisection happens to find, is the smooth
  // continuation (0.023 m from its theta neighbour, versus 0.119 m for the
  // first root). Leave the selection alone unless a case is measured where it
  // lands on a root that is not the neighbour-consistent one.
  double slo = 0.0, shi = 1.0, flo = fa;
  for (int k = 0; k < 44; ++k) {
    double smid = 0.5 * (slo + shi);
    double rm = ra + smid * (rb - ra), zm = za + smid * (zb - za);
    double fm = tok_eval_psi_rz_local(arc_ctx->geo, rm, zm) - psi;
    if (!isfinite(fm)) {
      return false;
    }
    if (flo * fm <= 0.0) {
      shi = smid;
    } else {
      slo = smid;
      flo = fm;
    }
  }
  double s = 0.5 * (slo + shi);
  *r = ra + s * (rb - ra);
  *z = za + s * (zb - za);
  double residual = tok_eval_psi_rz_local(arc_ctx->geo, *r, *z) - psi;
  return isfinite(*r) && isfinite(*z) && isfinite(residual) && fabs(residual) <= 1e-9 * scale;
}

static bool
tok_xpt_seam_endpoint(enum gkyl_tok_geo_type ftype, double u)
{
  const double endpoint_tol = 256.0 * DBL_EPSILON;
  return tok_sep_fixed_edge_is_first(ftype) ? u >= 1.0 - endpoint_tol : u <= endpoint_tol;
}

static bool
tok_eval_psi_grad_rz_local(
  const struct gkyl_tok_geo *geo, double R, double Z, double *dpsidR, double *dpsidZ
)
{
  if (geo->use_cubics) {
    double xn[2] = {R, Z}, out[3] = {0.0};
    geo->efit->evf->eval_cubic_wgrad(0.0, xn, out, geo->efit->evf->ctx);
    *dpsidR = out[1];
    *dpsidZ = out[2];
    return isfinite(*dpsidR) && isfinite(*dpsidZ);
  }
  int idx[2];
  idx[0] = GKYL_MIN2(
    geo->rzlocal.upper[0],
    GKYL_MAX2(
      geo->rzlocal.lower[0],
      geo->rzlocal.lower[0] + (int)floor((R - geo->rzgrid.lower[0]) / geo->rzgrid.dx[0])
    )
  );
  idx[1] = GKYL_MIN2(
    geo->rzlocal.upper[1],
    GKYL_MAX2(
      geo->rzlocal.lower[1],
      geo->rzlocal.lower[1] + (int)floor((Z - geo->rzgrid.lower[1]) / geo->rzgrid.dx[1])
    )
  );
  long loc = gkyl_range_idx(&geo->rzlocal, idx);
  const double *p = gkyl_array_cfetch(geo->psiRZ, loc);
  double xc[2];
  gkyl_rect_grid_cell_center(&geo->rzgrid, idx, xc);
  double x = (R - xc[0]) / (0.5 * geo->rzgrid.dx[0]);
  double y = (Z - xc[1]) / (0.5 * geo->rzgrid.dx[1]);
  if (geo->efit->rzbasis.poly_order == 1) {
    *dpsidR = (1.5 * p[3] * y + 0.8660254037844386 * p[1]) * 2.0 / geo->rzgrid.dx[0];
    *dpsidZ = (1.5 * p[3] * x + 0.8660254037844386 * p[2]) * 2.0 / geo->rzgrid.dx[1];
  } else {
    *dpsidR =
      (5.625 * p[8] * (2.0 * x * y * y - 0.6666666666666666 * x) +
       2.904737509655563 * p[7] * (y * y - 0.3333333333333333) + 5.809475019311126 * p[6] * x * y +
       1.5 * p[3] * y + 3.354101966249684 * p[4] * x + 0.8660254037844386 * p[1]) *
      2.0 / geo->rzgrid.dx[0];
    *dpsidZ = (5.625 * p[8] * (2.0 * x * x * y - 0.6666666666666666 * y) +
               5.809475019311126 * p[7] * x * y + 3.354101966249684 * p[5] * y +
               2.904737509655563 * p[6] * (x * x - 0.3333333333333333) + 1.5 * p[3] * x +
               0.8660254037844386 * p[2]) *
              2.0 / geo->rzgrid.dx[1];
  }
  return isfinite(*dpsidR) && isfinite(*dpsidZ);
}

// Project a continuation predictor to the exact requested DG/cubic flux
// surface. Holding the coordinate that changes most along the contour keeps
// the root solve well conditioned at either an R or Z turning point.
static bool
tok_project_xpt_seam_candidate(
  const struct gkyl_tok_geo *geo, double psi, double rprev, double zprev, double rpred,
  double zpred, double tangent_r, double tangent_z, double step, double *r, double *z
)
{
  const struct gkyl_rect_grid *grid = geo->use_cubics ? &geo->rzgrid_cubic : &geo->rzgrid;
  if (rpred < grid->lower[0] || rpred > grid->upper[0] || zpred < grid->lower[1] ||
      zpred > grid->upper[1]) {
    return false;
  }

  double best_distance = DBL_MAX, second_distance = DBL_MAX;
  if (fabs(tangent_r) >= fabs(tangent_z)) {
    double roots[32] = {0.0};
    int nr = tok_geo_Z_psiR(geo, psi, rpred, 32, roots);
    if (nr <= 0) {
      return false;
    }
    int best = 0;
    for (int i = 0; i < nr; ++i) {
      double distance = fabs(roots[i] - zpred);
      if (distance < best_distance) {
        second_distance = best_distance;
        best_distance = distance;
        best = i;
      } else if (distance < second_distance) {
        second_distance = distance;
      }
    }
    double ambiguity_scale = fmax(grid->dx[1], fabs(step));
    if (nr > 1 && second_distance - best_distance <= 1e-8 * fmax(1.0, ambiguity_scale)) {
      return false;
    }
    *r = rpred;
    *z = roots[best];
  } else {
    double roots[32] = {0.0}, dRdZ[32] = {0.0};
    double dR[32] = {0.0}, dZ[32] = {0.0};
    int nr = gkyl_tok_geo_R_psiZ(geo, psi, zpred, 32, roots, dRdZ, dR, dZ);
    if (nr <= 0) {
      return false;
    }
    int best = 0;
    for (int i = 0; i < nr; ++i) {
      double distance = fabs(roots[i] - rpred);
      if (distance < best_distance) {
        second_distance = best_distance;
        best_distance = distance;
        best = i;
      } else if (distance < second_distance) {
        second_distance = distance;
      }
    }
    double ambiguity_scale = fmax(grid->dx[0], fabs(step));
    if (nr > 1 && second_distance - best_distance <= 1e-8 * fmax(1.0, ambiguity_scale)) {
      return false;
    }
    *r = roots[best];
    *z = zpred;
  }

  double dr = *r - rprev, dz = *z - zprev;
  double forward = step * (dr * tangent_r + dz * tangent_z);
  double distance = hypot(dr, dz);
  double residual = tok_eval_psi_rz_local(geo, *r, *z) - psi;
  double flux_scale = fmax(1.0, fabs(psi));
  return isfinite(*r) && isfinite(*z) && isfinite(distance) && isfinite(residual) &&
         forward > 0.0 && distance <= 2.5 * fabs(step) + 1e-12 &&
         fabs(residual) <= 1e-10 * flux_scale;
}

static bool
tok_displace_xpt_seam_on_flux(
  const struct gkyl_tok_geo *geo, double psi, double delta_s, double *r, double *z,
  double *realized_out
)
{
  *realized_out = 0.0;
  if (delta_s == 0.0) {
    return true;
  }
  const struct gkyl_rect_grid *grid = geo->use_cubics ? &geo->rzgrid_cubic : &geo->rzgrid;
  double max_step = 0.05 * fmin(grid->dx[0], grid->dx[1]);
  if (!(max_step > 0.0) || !isfinite(max_step)) {
    return false;
  }
  double required_steps = ceil(fabs(delta_s) / max_step);
  if (!isfinite(required_steps) || required_steps > 4096.0) {
    return false;
  }
  int nstep = GKYL_MAX2(16, (int)required_steps);
  double rstart = *r, zstart = *z;
  double arc_length_error = DBL_MAX, realized = 0.0;
  // Combined absolute+relative acceptance. Near the B1(q) taper's
  // endpoints the requested per-node delta_s shrinks toward zero, so a
  // purely relative criterion against fabs(delta_s) is unreachable no
  // matter how many times nstep is doubled; the absolute floor keeps that
  // from being amplified into spurious rejections. Measured empirically:
  // even at the nstep cap (4096) the achievable absolute arc-length error
  // plateaus around 1e-8 m regardless of delta_s magnitude (set by the
  // underlying per-step Newton projection's own precision, not by
  // discretization), so 1e-10 is unreachable in both framings -- 1e-7
  // comfortably covers the observed error with margin while still being
  // far tighter than anything physically meaningful at this length scale.
  const double arc_length_rel_tol = 1e-7;
  const double arc_length_abs_tol = 1e-7;
  for (;;) {
    *r = rstart;
    *z = zstart;
    double step = delta_s / nstep;
    realized = 0.0;

    for (int i = 0; i < nstep; ++i) {
      double grad_r = 0.0, grad_z = 0.0;
      if (!tok_eval_psi_grad_rz_local(geo, *r, *z, &grad_r, &grad_z)) {
        return false;
      }
      double grad = hypot(grad_r, grad_z);
      if (!(grad > 1e-14)) {
        return false;
      }
      double tangent_r = -grad_z / grad, tangent_z = grad_r / grad;

      double rmid = *r + 0.5 * step * tangent_r;
      double zmid = *z + 0.5 * step * tangent_z;
      if (!tok_project_xpt_seam_candidate(
            geo, psi, *r, *z, rmid, zmid, tangent_r, tangent_z, 0.5 * step, &rmid, &zmid
          )) {
        return false;
      }
      if (!tok_eval_psi_grad_rz_local(geo, rmid, zmid, &grad_r, &grad_z)) {
        return false;
      }
      grad = hypot(grad_r, grad_z);
      if (!(grad > 1e-14)) {
        return false;
      }
      tangent_r = -grad_z / grad;
      tangent_z = grad_r / grad;

      double rnext = *r + step * tangent_r;
      double znext = *z + step * tangent_z;
      if (!tok_project_xpt_seam_candidate(
            geo, psi, *r, *z, rnext, znext, tangent_r, tangent_z, step, &rnext, &znext
          )) {
        return false;
      }
      realized += hypot(rnext - *r, znext - *z);
      *r = rnext;
      *z = znext;
    }
    double absolute_error = fabs(realized - fabs(delta_s));
    double error_tolerance = fmax(arc_length_rel_tol * fabs(delta_s), arc_length_abs_tol);
    arc_length_error = absolute_error;
    if (isfinite(absolute_error) && absolute_error <= error_tolerance) {
      *realized_out = realized;
      return true;
    }
    if (nstep > 4096 / 2) {
      break;
    }
    nstep *= 2;
  }
  *r = rstart;
  *z = zstart;
  return false;
}

// B1(q)=4q(1-q) is smooth, has unit peak, and uses explicit endpoint
// branches so both physical anchors remain exactly fixed.
static bool
tok_xpt_seam_delta_s(
  const struct gkyl_tok_geo_grid_inp *inp, const struct arc_length_ctx *arc_ctx, double psi,
  double *q, double *delta_s
)
{
  double span = arc_ctx->xpt_ray_psi0 - arc_ctx->geo->psisep;
  double scale = fmax(1.0, fmax(fabs(arc_ctx->xpt_ray_psi0), fabs(arc_ctx->geo->psisep)));
  if (!isfinite(span) || fabs(span) <= 256.0 * DBL_EPSILON * scale) {
    return false;
  }
  if (tok_geo_same_flux(psi, arc_ctx->geo->psisep)) {
    *q = 0.0;
    *delta_s = 0.0;
    return true;
  }
  if (tok_geo_same_flux(psi, arc_ctx->xpt_ray_psi0)) {
    *q = 1.0;
    *delta_s = 0.0;
    return true;
  }
  *q = (psi - arc_ctx->geo->psisep) / span;
  if (*q < -1e-8 || *q > 1.0 + 1e-8) {
    return false;
  }
  *q = fmin(1.0, fmax(0.0, *q));
  if (*q == 0.0 || *q == 1.0) {
    *delta_s = 0.0;
    return true;
  }
  double coefficient = inp->relaxed_xpt_seam_delta_s_coeff;
  double bound = inp->relaxed_xpt_seam_delta_s_bound;
  if (coefficient == 0.0) {
    *delta_s = 0.0;
    return true;
  }
  if (!isfinite(coefficient) || !isfinite(bound) || !(bound > 0.0) ||
      fabs(coefficient) > bound * (1.0 + 64.0 * DBL_EPSILON)) {
    return false;
  }
  *delta_s = coefficient * 4.0 * (*q) * (1.0 - *q);
  return isfinite(*delta_s) && fabs(*delta_s) <= bound * (1.0 + 64.0 * DBL_EPSILON);
}

// Relax an already-validated seam-adjacent point (r,z), computed by
// whichever base construction was used for this ftype/domain mode (the
// half-domain chord point or the full-domain/extended ordered-map trace),
// by the bounded delta-s displacement. Contour-agnostic: the displacement
// primitive (tok_displace_xpt_seam_on_flux) only walks along the flux
// contour already passing through the given (r,z), so it composes safely
// regardless of which point-construction algorithm produced it -- this is
// what lets the same relaxation be reused for both half_domain=true and
// half_domain=false (extended) blocks instead of being wired only into the
// half-domain-only path it originated in.
static bool
tok_relax_xpt_seam_point(
  const struct gkyl_tok_geo_grid_inp *inp, struct arc_length_ctx *arc_ctx, double psi, double u,
  double *r, double *z
)
{
  if (!inp->relaxed_xpt_seam) {
    return true;
  }

  double q = 0.0, delta_s = 0.0;
  if (!tok_xpt_seam_delta_s(inp, arc_ctx, psi, &q, &delta_s)) {
    fprintf(
      stderr,
      "TOK_XPT_SEAM_DIAG mode=rejected ftype=%d psi=%.17g q=%.17g coefficient=%.17g bound=%.17g reason=invalid_delta_s action=%s\n",
      inp->ftype, psi, q, inp->relaxed_xpt_seam_delta_s_coeff, inp->relaxed_xpt_seam_delta_s_bound,
      inp->relaxed_xpt_seam_sweep ? "reject_candidate" : "retain_straight"
    );
    if (tok_xpt_seam_optimizer_trial(inp)) {
      tok_xpt_seam_trial_reject(inp, GKYL_XPT_SEAM_TRIAL_INVALID_PARAMETER);
      return true;
    }
    return !inp->relaxed_xpt_seam_sweep;
  }

  // Use signed local contour arc length, with s0=0 at the straight-ray
  // intersection. Return before doing coordinate arithmetic so zero mode is
  // bitwise identical to the existing construction.
  const double s0 = 0.0;
  if (delta_s == 0.0) {
    return true;
  }

  if (!inp->relaxed_xpt_seam_sweep || !inp->straight_xpt_ray) {
    fprintf(
      stderr,
      "TOK_XPT_SEAM_DIAG mode=rejected ftype=%d psi=%.17g q=%.17g s0=%.17g delta_s=%.17g reason=diagnostic_sweep_not_enabled fallback=straight\n",
      inp->ftype, psi, q, s0, delta_s
    );
    return true;
  }

  // The only free data are on the seam. A fixed linear arc weight carries
  // that boundary condition to the unchanged opposite edge, giving a smooth
  // diagnostic shadow map without introducing interior degrees of freedom.
  // A self-periodic single-null CORE block has no such opposite edge -- u=0
  // and u=1 are the same physical (X-point-adjacent) location, glued
  // together by the closed trace construction -- so a linear taper would
  // apply delta_s at one end and not the other, breaking the periodic
  // identification it is supposed to preserve. Apply it uniformly instead:
  // every point at this psi gets the same q-dependent radial bump
  // regardless of u, so u=0 and u=1 (which sample the same underlying raw
  // trace point) remain identically displaced and the wrap stays exact.
  double seam_weight = inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE ? 1.0 :
                       tok_sep_fixed_edge_is_first(inp->ftype)  ? u :
                                                                  1.0 - u;
  double point_delta_s = seam_weight * delta_s;
  if (point_delta_s == 0.0) {
    return true;
  }
  double candidate_r = *r, candidate_z = *z;
  double realized = 0.0;
  if (!tok_displace_xpt_seam_on_flux(
        arc_ctx->geo, psi, point_delta_s, &candidate_r, &candidate_z, &realized
      )) {
    fprintf(
      stderr,
      "TOK_XPT_SEAM_DIAG mode=rejected ftype=%d psi=%.17g q=%.17g u=%.17g delta_s=%.17g point_delta_s=%.17g reason=contour_tracking_failed action=reject_candidate fallback=none\n",
      inp->ftype, psi, q, u, delta_s, point_delta_s
    );
    if (tok_xpt_seam_optimizer_trial(inp)) {
      tok_xpt_seam_trial_reject(inp, GKYL_XPT_SEAM_TRIAL_CONTOUR);
      return true;
    }
    return false;
  }
  if (tok_xpt_seam_optimizer_trial(inp)) {
    struct gkyl_tok_geo_xpt_seam_trial_status *status = inp->relaxed_xpt_seam_trial_status;
    status->max_realized_displacement = fmax(status->max_realized_displacement, realized);
    double bound_tolerance = 64.0 * DBL_EPSILON * fmax(1.0, inp->relaxed_xpt_seam_delta_s_bound);
    if (realized > inp->relaxed_xpt_seam_delta_s_bound + bound_tolerance) {
      tok_xpt_seam_trial_reject(inp, GKYL_XPT_SEAM_TRIAL_CONTOUR);
      return true;
    }
  }
  if (tok_xpt_seam_endpoint(inp->ftype, u) && arc_ctx->xpt_ray_branch_valid) {
    bool resolved = false, on_right = false;
    if (!tok_xpt_classify_branch_at_point(
          inp, arc_ctx, candidate_r, candidate_z, psi, &resolved, &on_right
        ) ||
        !resolved || on_right != arc_ctx->xpt_ray_on_right) {
      if (tok_xpt_seam_optimizer_trial(inp)) {
        tok_xpt_seam_trial_reject(inp, GKYL_XPT_SEAM_TRIAL_BRANCH);
        return true;
      }
      fprintf(
        stderr,
        "TOK_XPT_SEAM_DIAG mode=rejected ftype=%d psi=%.17g q=%.17g delta_s=%.17g reason=branch_identity_changed action=reject_candidate fallback=none\n",
        inp->ftype, psi, q, delta_s
      );
      return false;
    }
  }
  *r = candidate_r;
  *z = candidate_z;
  return true;
}

// Half-domain (non-extended) seam point: the validated straight-ray chord
// point, relaxed by delta-s. Semantically unchanged from the pre-refactor
// implementation -- previously this function's own body started here, gated
// by "!inp->half_domain" in addition to "!inp->relaxed_xpt_seam"; since this
// call site is only ever reached when half_domain=true or straight_xpt_ray
// is false (see tok_build_current_ordered_trace/tok_ordered_map_lookup),
// that half_domain check was always true here and is preserved by construction
// rather than by an explicit condition.
static bool
tok_parameterized_xpt_seam_point(
  const struct gkyl_tok_geo_grid_inp *inp, struct arc_length_ctx *arc_ctx, double psi, double u,
  double *r, double *z
)
{
  if (!tok_ordered_chord_point(inp, arc_ctx, psi, u, r, z)) {
    return false;
  }
  return tok_relax_xpt_seam_point(inp, arc_ctx, psi, u, r, z);
}

static double
tok_fpol_at_psi(const struct gkyl_tok_geo *geo, double psi)
{
  double p = psi;
  if (p < geo->fgrid.lower[0] || p > geo->fgrid.upper[0]) {
    p = geo->sibry;
  }
  int idx = GKYL_MIN2(
    geo->frange.upper[0],
    GKYL_MAX2(
      geo->frange.lower[0],
      geo->frange.lower[0] + (int)floor((p - geo->fgrid.lower[0]) / geo->fgrid.dx[0])
    )
  );
  long loc = gkyl_range_idx(&geo->frange, &idx);
  const double *coeffs = gkyl_array_cfetch(geo->fpoldg, loc);
  double xc;
  gkyl_rect_grid_cell_center(&geo->fgrid, &idx, &xc);
  double x = (p - xc) / (0.5 * geo->fgrid.dx[0]);
  return geo->fbasis.eval_expand(&x, coeffs);
}

static bool
tok_ext_set_phi_reference(const struct gkyl_tok_geo_grid_inp *inp, struct arc_length_ctx *arc_ctx)
{
  struct tok_ext_topology top;
  if (!tok_ext_topology_from_ftype(inp->ftype, inp->half_domain, &top) ||
      arc_ctx->map_trace_n < 2) {
    return false;
  }
  if (top.phi_reference == TOK_EXT_PHI_LOWER) {
    arc_ctx->map_trace_phi_ref = arc_ctx->map_trace_phi[0];
    return true;
  }
  if (top.phi_reference == TOK_EXT_PHI_UPPER) {
    arc_ctx->map_trace_phi_ref = arc_ctx->map_trace_phi[arc_ctx->map_trace_n - 1];
    return true;
  }

  double R[16] = {0.0}, dRdZ[16] = {0.0};
  double dR[16] = {0.0}, dZ[16] = {0.0};
  int nr =
    gkyl_tok_geo_R_psiZ(arc_ctx->geo, arc_ctx->psi, arc_ctx->geo->zmaxis, 16, R, dRdZ, dR, dZ);
  if (nr <= 0) {
    return false;
  }
  bool outboard = top.phi_reference == TOK_EXT_PHI_OUTBOARD_MIDPLANE;
  double target_z = arc_ctx->geo->zmaxis;
  // Seed the midplane reference from THIS BLOCK'S OWN TRACE rather than from
  // the declared rright/rleft.
  //
  // The trace is already on the surface, so where it crosses Z=zmaxis IS the
  // midplane point, and its outboard branch is simply the crossing at larger
  // R.  Choosing by proximity to a declared rright is a device-specific
  // assumption, and it fails whenever the level set has a further component:
  // measured on tcv at psisep, the roots at zmaxis are {0.6821953, 1.0991599}
  // while the trace spans R=[0.6799185, 0.8623851].  The real outboard
  // crossing near 0.86 is missing from the root list, and rright=0.9 is
  // nearer 1.0991599 -- a point 0.24 m outside the block -- so the reference
  // landed off the trace and the whole map failed.  The trace cannot make
  // that mistake: it only contains points of this block.
  //
  // The seed is then refined against psi itself, so the reference stays an
  // exact root and does not inherit the trace's piecewise-linear error.
  double seed_r = 0.0;
  bool have_seed = false;
  for (int i = 0; i < arc_ctx->map_trace_n - 1; ++i) {
    const double z0 = arc_ctx->map_trace_z[i], z1 = arc_ctx->map_trace_z[i + 1];
    if ((z0 - target_z) * (z1 - target_z) > 0.0) {
      continue;
    }
    const double dz = z1 - z0;
    const double w = fabs(dz) > 0.0 ? (target_z - z0) / dz : 0.0;
    const double rc =
      arc_ctx->map_trace_r[i] + w * (arc_ctx->map_trace_r[i + 1] - arc_ctx->map_trace_r[i]);
    if (!isfinite(rc)) {
      continue;
    }
    if (!have_seed || (outboard ? rc > seed_r : rc < seed_r)) {
      seed_r = rc;
      have_seed = true;
    }
  }
  // Refine the seed on psi itself, so the reference is an exact root and does
  // not inherit the trace's piecewise-linear error.  Falls back to the old
  // declared-side choice only when this block's trace never reaches zmaxis.
  double target_r = outboard ? inp->rright : inp->rleft;
  if (have_seed) {
    double a = seed_r;
    double fa = tok_eval_psi_rz_local(arc_ctx->geo, a, target_z) - arc_ctx->psi;
    for (int it = 0; it < 64; ++it) {
      if (!isfinite(fa)) {
        break;
      }
      if (fabs(fa) <= 1e-12 * fmax(1.0, fabs(arc_ctx->psi))) {
        break;
      }
      double gr = 0.0, gz = 0.0;
      if (!tok_eval_psi_grad_rz_local(arc_ctx->geo, a, target_z, &gr, &gz) || !(fabs(gr) > 0.0)) {
        break;
      }
      a -= fa / gr;
      fa = tok_eval_psi_rz_local(arc_ctx->geo, a, target_z) - arc_ctx->psi;
    }
    target_r = (isfinite(a) && isfinite(fa) && fabs(fa) <= 1e-9 * fmax(1.0, fabs(arc_ctx->psi))) ?
                 a :
                 seed_r;
  } else {
    target_r = tok_nearest_value(target_r, R, nr);
  }
  double best_d2 = DBL_MAX, best_phi = 0.0, max_step = 0.0, best_w = 0.0;
  int best_i = -1;
  for (int i = 0; i < arc_ctx->map_trace_n - 1; ++i) {
    double r0 = arc_ctx->map_trace_r[i];
    double z0 = arc_ctx->map_trace_z[i];
    double dr = arc_ctx->map_trace_r[i + 1] - r0;
    double dz = arc_ctx->map_trace_z[i + 1] - z0;
    double den = dr * dr + dz * dz;
    if (!(den > 0.0)) {
      continue;
    }
    max_step = fmax(max_step, sqrt(den));
    double w = ((target_r - r0) * dr + (target_z - z0) * dz) / den;
    w = fmin(1.0, fmax(0.0, w));
    double rp = r0 + w * dr, zp = z0 + w * dz;
    double d2 = SQ(target_r - rp) + SQ(target_z - zp);
    if (d2 < best_d2) {
      best_d2 = d2;
      best_phi =
        arc_ctx->map_trace_phi[i] + w * (arc_ctx->map_trace_phi[i + 1] - arc_ctx->map_trace_phi[i]);
      best_i = i;
      best_w = w;
    }
  }
  if (!isfinite(best_phi) || !isfinite(best_d2) || !(max_step > 0.0) ||
      sqrt(best_d2) > 2.0 * max_step) {
    double tr_rmin = DBL_MAX, tr_rmax = -DBL_MAX;
    double tr_zmin = DBL_MAX, tr_zmax = -DBL_MAX;
    for (int i = 0; i < arc_ctx->map_trace_n; ++i) {
      tr_rmin = fmin(tr_rmin, arc_ctx->map_trace_r[i]);
      tr_rmax = fmax(tr_rmax, arc_ctx->map_trace_r[i]);
      tr_zmin = fmin(tr_zmin, arc_ctx->map_trace_z[i]);
      tr_zmax = fmax(tr_zmax, arc_ctx->map_trace_z[i]);
    }
    fprintf(
      stderr,
      "TOK_ORDERED_MAP failed phi reference ftype=%d psi=%.17g side=%s "
      "distance=%.17g max_step=%.17g target=(%.7f,%.7f) zmaxis=%.7f "
      "rright=%.7f rleft=%.7f nroots=%d trace_n=%d "
      "trace_R=[%.7f,%.7f] trace_Z=[%.7f,%.7f] ends=(%.7f,%.7f)-(%.7f,%.7f)\n",
      inp->ftype, arc_ctx->psi, outboard ? "outboard" : "inboard", sqrt(best_d2), max_step,
      target_r, target_z, arc_ctx->geo->zmaxis, inp->rright, inp->rleft, nr, arc_ctx->map_trace_n,
      tr_rmin, tr_rmax, tr_zmin, tr_zmax, arc_ctx->map_trace_r[0], arc_ctx->map_trace_z[0],
      arc_ctx->map_trace_r[arc_ctx->map_trace_n - 1], arc_ctx->map_trace_z[arc_ctx->map_trace_n - 1]
    );
    for (int k = 0; k < nr && k < 16; ++k) {
      fprintf(stderr, "TOK_ORDERED_MAP phi_ref_root[%d]=%.7f\n", k, R[k]);
    }
    return false;
  }
  // The reference is a root of psi, so take the angle there on the contour
  // (see the true-arc angle section) rather than along the bracket's chord.
  if (best_i >= 0) {
    const double fpol = tok_fpol_at_psi(arc_ctx->geo, arc_ctx->psi);
    double d = 0.0;
    if (best_w <= 0.0) {
      best_phi = arc_ctx->map_trace_phi[best_i];
    } else if (best_w >= 1.0) {
      best_phi = arc_ctx->map_trace_phi[best_i + 1];
    } else if (tok_phi_into_bracket(
                 arc_ctx->geo, arc_ctx->psi, fpol, arc_ctx->map_trace_r[best_i],
                 arc_ctx->map_trace_z[best_i], arc_ctx->map_trace_r[best_i + 1],
                 arc_ctx->map_trace_z[best_i + 1],
                 arc_ctx->map_trace_phi[best_i + 1] - arc_ctx->map_trace_phi[best_i], target_r,
                 target_z, &d
               )) {
      best_phi = arc_ctx->map_trace_phi[best_i] + d;
    } else {
      fprintf(
        stderr, "TOK_PHI_EXACT_FALLBACK stage=reference ftype=%d psi=%.17g i=%d\n", inp->ftype,
        arc_ctx->psi, best_i
      );
    }
  }
  arc_ctx->map_trace_phi_ref = best_phi;
  return true;
}

// Project a point onto psi=const along grad psi -- the level set's own normal
// -- rather than along a coordinate axis.  Used where the axis-aligned solve
// has no solution at all: on the trace segment that leaves an X point the
// contour turns a corner, so the chord between two on-contour trace points
// bows onto the concave side and the line through it misses both branches.
// The displacement is bounded because |grad psi| collapses at the X point and
// an unclamped Newton step can land on a different separatrix branch.
static bool
tok_psi_normal_project(
  const struct gkyl_tok_geo *geo, double psi, double r0, double z0, double max_disp, double tol,
  double *r, double *z
)
{
  double cr = r0, cz = z0;
  double step_cap = 0.25 * max_disp;
  for (int k = 0; k < 32; ++k) {
    double f = tok_eval_psi_rz_local(geo, cr, cz) - psi;
    if (!isfinite(f)) {
      return false;
    }
    if (fabs(f) <= tol) {
      *r = cr;
      *z = cz;
      return true;
    }
    double gr = 0.0, gz = 0.0;
    if (!tok_eval_psi_grad_rz_local(geo, cr, cz, &gr, &gz)) {
      return false;
    }
    double g2 = gr * gr + gz * gz;
    if (!(g2 > 0.0)) {
      return false;
    }
    double dr = -f / g2 * gr, dz = -f / g2 * gz;
    double dl = hypot(dr, dz);
    if (dl > step_cap) {
      dr *= step_cap / dl;
      dz *= step_cap / dl;
    }
    cr += dr;
    cz += dz;
    if (!isfinite(cr) || !isfinite(cz)) {
      return false;
    }
    if (hypot(cr - r0, cz - z0) > max_disp) {
      return false;
    }
  }
  double f = tok_eval_psi_rz_local(geo, cr, cz) - psi;
  if (isfinite(f) && fabs(f) <= tol) {
    *r = cr;
    *z = cz;
    return true;
  }
  return false;
}

// ---------------------------------------------------------------------------
// Single-null SOL rows on the legacy construction, at exact arc (2026-10-06).
//
// LSN_SOL_* blocks without straight_xpt_ray place each node by root-solving the
// Z-integrated arc, integral sqrt(1+(dR/dZ)^2) dZ, and take their theta
// extents from the separatrix's arc so integrated. That parameterisation is
// weakest exactly where these blocks meet: the integrand is singular at the
// row's Z turning point (the top), and the separatrix has a corner at the X
// point. Measured on ASDEX and TCV: rows uniform to ~1e-6 elsewhere, but the
// cell at the top off by up to 7e-4, and the separatrix cell beside the X
// point not converging (4.5e-5 / 1.8e-7 / 2.1e-4 at theta x1 / x2 / x4), so
// the faces between these blocks did not converge in theta.
//
// Here a row is traced over the branches the integral follows -- the outer
// strike up the right side (the root nearest rright at each grid line in Z),
// over the top, down the left side (nearest rleft) to the inner strike, the X
// point itself on the separatrix -- and measured, and sampled, with the
// arc-exact sampler's refinement (tok_contour_leaves, tok_trace_sample). The
// cut rule is unchanged: block boundaries remain fractions of the separatrix's
// arc, now exact, so the separatrix cut is exactly the X point and every other
// row is cut at the same fraction of its own exact arc -- the legacy cuts.
// ---------------------------------------------------------------------------

struct tok_lsn_row {
  double psi;
  bool valid;
  int n, cap;
  double *r, *z, *s; // the row's trace and its exact cumulative arc
  // the field-line angle along the trace (see the true-arc angle section),
  // from the outer strike, and its value where the row crosses the outboard
  // midplane
  double *p, phi_mid;
  bool phi_valid, have_mid;
};
static _Thread_local struct tok_lsn_row tok_lsn_row_buf;

static void
tok_lsn_push(struct tok_lsn_row *b, double r, double z)
{
  if (b->n == b->cap) {
    int cap = b->cap ? 2 * b->cap : 256;
    double *nr = gkyl_malloc(cap * sizeof(double)), *nz = gkyl_malloc(cap * sizeof(double)),
           *ns = gkyl_malloc(cap * sizeof(double)), *np = gkyl_malloc(cap * sizeof(double));
    if (b->n) {
      memcpy(nr, b->r, b->n * sizeof(double));
      memcpy(nz, b->z, b->n * sizeof(double));
      gkyl_free(b->r);
      gkyl_free(b->z);
      gkyl_free(b->s);
      gkyl_free(b->p);
    }
    b->r = nr;
    b->z = nz;
    b->s = ns;
    b->p = np;
    b->cap = cap;
  }
  b->r[b->n] = r;
  b->z[b->n] = z;
  ++b->n;
}

// The root of the psi contour at height z nearest rref, as the legacy
// integrand takes it. False if the contour does not reach z.
static bool
tok_lsn_add(const struct gkyl_tok_geo *geo, double psi, double z, double rref, struct tok_lsn_row *b)
{
  double R[8] = {0.0}, dRdZ[8] = {0.0}, dR[8] = {0.0}, dZ[8] = {0.0};
  int nr = gkyl_tok_geo_R_psiZ(geo, psi, z, 8, R, dRdZ, dR, dZ);
  if (nr <= 0) {
    return false;
  }
  tok_lsn_push(b, choose_closest(rref, R, R, nr), z);
  return true;
}

bool
tok_lsn_exact_row(
  const struct gkyl_tok_geo *geo, double psi, double zmin_right, double zmax, double zmin_left,
  double rright, double rleft, double *arc_right, double *arc_tot, double sep_arcs[4]
)
{
  struct tok_lsn_row *b = &tok_lsn_row_buf;
  b->valid = false;
  b->phi_valid = false;
  b->have_mid = false;
  b->n = 0;
  if (!(zmax > zmin_right) || !(zmax > zmin_left)) {
    return false;
  }
  const bool sep = tok_geo_same_flux(psi, geo->psisep);
  double rx = 0.0, zx = 0.0;
  if (sep && !tok_ext_xpoint_rz(geo, TOK_EXT_LOWER_XPT, &rx, &zx)) {
    return false;
  }
  const bool xr = sep && zx > zmin_right && zx < zmax;
  const bool xl = sep && zx > zmin_left && zx < zmax;
  const struct gkyl_rect_grid *g = geo->use_cubics ? &geo->rzgrid_cubic : &geo->rzgrid;
  const double z0 = g->lower[1], dz = g->dx[1];
  const double ztol = 64.0 * DBL_EPSILON * fmax(1.0, fabs(zmax));
  int i_xr = -1, i_top = -1, i_xl = -1;
  // right side, upward: the outer strike, every grid line in Z, the X point
  if (!tok_lsn_add(geo, psi, zmin_right, rright, b)) {
    return false;
  }
  for (int k = (int)floor((zmin_right - z0) / dz) + 1; z0 + k * dz < zmax; ++k) {
    const double zk = z0 + k * dz;
    if (xr && i_xr < 0 && zx <= zk) {
      tok_lsn_push(b, rx, zx);
      i_xr = b->n - 1;
      if (fabs(zk - zx) <= ztol) {
        continue;
      }
    }
    if (zk > zmin_right) {
      tok_lsn_add(geo, psi, zk, rright, b);
    }
  }
  if (xr && i_xr < 0) {
    tok_lsn_push(b, rx, zx);
    i_xr = b->n - 1;
  }
  // the top: the two roots at the turning height when it has them; the piece
  // between the last right point and the first left one is measured exactly
  // either way
  tok_lsn_add(geo, psi, zmax, rright, b);
  i_top = b->n - 1;
  const int n_right = b->n;
  tok_lsn_add(geo, psi, zmax, rleft, b);
  // left side, downward
  for (int k = (int)ceil((zmax - z0) / dz) - 1; z0 + k * dz > zmin_left; --k) {
    const double zk = z0 + k * dz;
    if (zk >= zmax) {
      continue;
    }
    if (xl && i_xl < 0 && zx >= zk) {
      tok_lsn_push(b, rx, zx);
      i_xl = b->n - 1;
      if (fabs(zk - zx) <= ztol) {
        continue;
      }
    }
    tok_lsn_add(geo, psi, zk, rleft, b);
  }
  if (xl && i_xl < 0) {
    tok_lsn_push(b, rx, zx);
    i_xl = b->n - 1;
  }
  if (!tok_lsn_add(geo, psi, zmin_left, rleft, b)) {
    return false;
  }
  if (b->n < 3 || n_right < 2) {
    return false;
  }
  // exact cumulative arc at the trace's points
  struct tok_leaves lv = {0};
  b->s[0] = 0.0;
  bool ok = true;
  for (int i = 1; i < b->n && ok; ++i) {
    lv.n = 0;
    ok = tok_contour_leaves(geo, psi, b->r[i - 1], b->z[i - 1], b->r[i], b->z[i], 0, &lv);
    double a = 0.0;
    for (int j = 0; j < lv.n; ++j) {
      a += lv.arc[j];
    }
    b->s[i] = b->s[i - 1] + a;
  }
  if (lv.r) {
    gkyl_free(lv.r);
    gkyl_free(lv.z);
    gkyl_free(lv.arc);
  }
  if (!ok || !(b->s[b->n - 1] > 0.0)) {
    return false;
  }
  if (arc_right) {
    *arc_right = b->s[i_top];
  }
  if (arc_tot) {
    *arc_tot = b->s[b->n - 1];
  }
  if (sep_arcs) {
    if (i_xr < 0 || i_xl < 0) {
      return false;
    }
    sep_arcs[0] = b->s[i_xr];
    sep_arcs[1] = b->s[i_top] - b->s[i_xr];
    sep_arcs[2] = b->s[i_xl] - b->s[i_top];
    sep_arcs[3] = b->s[b->n - 1] - b->s[i_xl];
  }
  // The field-line angle along the same trace, integrated on the contour
  // between its points (tok_phi_integral), and where the row crosses the
  // outboard midplane -- the legacy reference of LSN_SOL and LSN_SOL_MID. The
  // legacy angle was a Z integral through the row's turning point, and its
  // error there (up to 1e-3 rad, measured on ASDEX 33292 LSN_SOL at 32 x 24
  // cells) offset the whole inboard half of each row by a different amount:
  // streaks in g^12, g^22, g^23 on the inboard side only.
  const double fpol = tok_fpol_at_psi(geo, psi);
  bool pok = true;
  b->p[0] = 0.0;
  for (int i = 1; i < b->n && pok; ++i) {
    double d = 0.0;
    pok = tok_phi_integral(geo, psi, fpol, b->r[i - 1], b->z[i - 1], b->r[i], b->z[i], &d);
    b->p[i] = b->p[i - 1] + d;
  }
  for (int i = 0; pok && i < i_top; ++i) {
    if (!(b->z[i] <= geo->zmaxis && geo->zmaxis < b->z[i + 1])) {
      continue;
    }
    double R[8] = {0.0}, dRdZ[8] = {0.0}, dR[8] = {0.0}, dZ[8] = {0.0};
    const int nr = gkyl_tok_geo_R_psiZ(geo, psi, geo->zmaxis, 8, R, dRdZ, dR, dZ);
    double d = 0.0;
    if (nr > 0 && tok_phi_into_bracket(
                    geo, psi, fpol, b->r[i], b->z[i], b->r[i + 1], b->z[i + 1],
                    b->p[i + 1] - b->p[i], choose_closest(rright, R, R, nr), geo->zmaxis, &d
                  )) {
      b->phi_mid = b->p[i] + d;
      b->have_mid = true;
    }
    break;
  }
  b->phi_valid = pok;
  if (!pok) {
    fprintf(stderr, "TOK_PHI_EXACT_FALLBACK stage=lsn_row psi=%.17g n=%d\n", psi, b->n);
  }
  b->psi = psi;
  b->valid = true;
  return true;
}

// The field-line angle (less alpha) at the node (r, z) at exact arc fraction u
// of the row tok_lsn_exact_row last traced, measured from the legacy reference
// of this block: the outboard midplane (LSN_SOL, LSN_SOL_MID), the outer strike
// (LSN_SOL_LO, increasing towards the X point), the inner strike (LSN_SOL_UP).
static bool
tok_lsn_exact_phi(
  const struct gkyl_tok_geo *geo, enum gkyl_tok_geo_type ftype, double psi, double u, double r,
  double z, double *phi
)
{
  const struct tok_lsn_row *b = &tok_lsn_row_buf;
  if (!b->valid || !b->phi_valid || b->psi != psi || !(u >= 0.0) || !(u <= 1.0)) {
    return false;
  }
  double ref = 0.0;
  if (ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL || ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID) {
    if (!b->have_mid) {
      return false;
    }
    ref = b->phi_mid;
  } else if (ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_UP) {
    ref = b->p[b->n - 1];
  } else if (ftype != GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO) {
    return false;
  }
  const double want = u * b->s[b->n - 1];
  int a = 0, c = b->n - 1;
  while (c - a > 1) {
    const int mid = (a + c) / 2;
    if (b->s[mid] <= want) {
      a = mid;
    } else {
      c = mid;
    }
  }
  double d = 0.0;
  if (r == b->r[c] && z == b->z[c]) {
    d = b->p[c] - b->p[a];
  } else if (!tok_phi_into_bracket(
               geo, psi, tok_fpol_at_psi(geo, psi), b->r[a], b->z[a], b->r[c], b->z[c],
               b->p[c] - b->p[a], r, z, &d
             )) {
    return false;
  }
  *phi = b->p[a] + d - ref;
  return true;
}

// The node at exact arc fraction u of the row tok_lsn_exact_row last traced,
// if that row is this psi's.
static bool
tok_lsn_exact_point(const struct gkyl_tok_geo *geo, double psi, double u, double *r, double *z)
{
  const struct tok_lsn_row *b = &tok_lsn_row_buf;
  if (!b->valid || b->psi != psi || !(u >= 0.0) || !(u <= 1.0)) {
    return false;
  }
  return tok_trace_sample(geo, psi, b->r, b->z, b->s, b->n, false, u, r, z);
}

// Sample a trace whose array index, rather than physical arc length, is
// uniform in the logical block coordinate.  Project the local interpolation
// back to psi=constant using whichever physical coordinate varies most on the
// segment; this remains well conditioned at R and Z turning points.
static bool
tok_logical_trace_sample(
  const struct gkyl_tok_geo *geo, double psi, const double *tr, const double *tz, int n, double u,
  double *r, double *z
)
{
  if (n < 2) {
    return false;
  }
  const double endpoint_tol = 256.0 * DBL_EPSILON;
  if (u <= endpoint_tol) {
    *r = tr[0];
    *z = tz[0];
    return true;
  }
  if (u >= 1.0 - endpoint_tol) {
    *r = tr[n - 1];
    *z = tz[n - 1];
    return true;
  }
  double x = u * (n - 1);
  int i = GKYL_MIN2(n - 2, GKYL_MAX2(0, (int)floor(x)));
  double w = x - i;
  double rlin = tr[i] + w * (tr[i + 1] - tr[i]);
  double zlin = tz[i] + w * (tz[i + 1] - tz[i]);
  double dr = tr[i + 1] - tr[i], dz = tz[i + 1] - tz[i];
  const double tol = 1e-9 * fmax(1.0, fabs(psi));
  // A3: the same bracket-frame rule as tok_trace_sample, first.  The map
  // trace's brackets are a quarter cell long, and on the bracket leaving an
  // X point the axis-aligned solve below fails and the gradient projection
  // after it slides toward the saddle (see tok_chord_normal_solve).
  if (w <= 0.0) {
    *r = tr[i];
    *z = tz[i];
    return true;
  }
  if (w >= 1.0) {
    *r = tr[i + 1];
    *z = tz[i + 1];
    return true;
  }
  double cr = 0.0, cz = 0.0;
  if (tok_chord_normal_solve(geo, psi, rlin, zlin, dr, dz, 0.5 * hypot(dr, dz), &cr, &cz)) {
    *r = cr;
    *z = cz;
    return true;
  }
  bool axis_ok = true;
  if (fabs(dr) >= fabs(dz)) {
    double roots[32] = {0.0};
    int nr = tok_geo_Z_psiR(geo, psi, rlin, 32, roots);
    if (nr <= 0) {
      axis_ok = false;
    } else {
      *r = rlin;
      *z = tok_nearest_value(zlin, roots, nr);
    }
  } else {
    double roots[16] = {0.0}, dRdZ[16] = {0.0};
    double dR[16] = {0.0}, dZ[16] = {0.0};
    int nr = gkyl_tok_geo_R_psiZ(geo, psi, zlin, 16, roots, dRdZ, dR, dZ);
    if (nr <= 0) {
      axis_ok = false;
    } else {
      double r_choice = tok_nearest_value(rlin, roots, nr);
      if (tok_rlin_ambiguous(rlin, roots, nr)) {
        double r_walk = 0.0;
        if (tok_trace_walk_root(geo, psi, tz[i], tr[i], tz[i + 1], tr[i + 1], zlin, &r_walk)) {
          r_choice = r_walk;
        }
      }
      *r = r_choice;
      *z = zlin;
    }
  }
  double seglen = hypot(dr, dz);
  if (axis_ok) {
    // The axis-aligned solve takes the root nearest the chord point, but when
    // no root is near it returns a distant one instead of failing: on an
    // inboard SOL surface the same R column meets the contour again far down
    // the divertor leg, and that root silently replaces the intended one.
    // Both bracketing trace points are on the contour, so the answer cannot
    // be further from the chord than the segment itself.
    double away = hypot(*r - rlin, *z - zlin);
    if (away > 0.5 * seglen) {
      axis_ok = false;
    }
  }
  if (axis_ok) {
    double residual = tok_eval_psi_rz_local(geo, *r, *z) - psi;
    if (isfinite(*r) && isfinite(*z) && isfinite(residual) && fabs(residual) <= tol) {
      return true;
    }
  }
  // The axis-aligned solve has no answer on the segment that leaves an X
  // point: the contour turns a corner there, so the chord between two
  // on-contour trace points bows onto the concave side and neither the
  // vertical nor the horizontal line through it meets psi=const.  Both trace
  // points are exact, so project the chord point back along the level set's
  // own normal instead, bounded by half the segment so the result stays on
  // this segment's branch.
  double pr = 0.0, pz = 0.0;
  if (seglen > 0.0 && tok_psi_normal_project(geo, psi, rlin, zlin, 0.5 * seglen, tol, &pr, &pz)) {
    *r = pr;
    *z = pz;
    return true;
  }
  // Last resort, and the only one that survives a DG cell face.  The quadratic
  // psi representation is C0: psi is continuous across a face but grad psi is
  // not, so a Newton projection started on a face -- which is exactly where
  // 204980's segment sits, straddling Z = -7*dZ = -0.9625 to 13 digits -- has
  // no well-defined direction to move in and the step above fails.  Bisecting
  // for a sign change of psi-target along the segment normal never evaluates a
  // gradient, so the kink cannot defeat it.  Both bracketing trace points are
  // exact contour points, so a crossing exists within a segment length; bound
  // the search there to stay on this segment's branch.
  if (seglen > 0.0) {
    double nx = -dz / seglen, ny = dr / seglen;
    const int nstep = 32;
    double h = seglen / nstep;
    double f0 = tok_eval_psi_rz_local(geo, rlin, zlin) - psi;
    if (isfinite(f0)) {
      double fprev_p = f0, fprev_m = f0;
      for (int k = 1; k <= nstep; ++k) {
        for (int side = 0; side < 2; ++side) {
          double sgn = side == 0 ? 1.0 : -1.0;
          double t = sgn * k * h;
          double fk = tok_eval_psi_rz_local(geo, rlin + t * nx, zlin + t * ny) - psi;
          double fprev = side == 0 ? fprev_p : fprev_m;
          if (isfinite(fk) && isfinite(fprev) && (fk < 0.0) != (fprev < 0.0)) {
            double lo = sgn * (k - 1) * h, hi = t, flo = fprev;
            for (int m = 0; m < 80; ++m) {
              double mid = 0.5 * (lo + hi);
              double fm = tok_eval_psi_rz_local(geo, rlin + mid * nx, zlin + mid * ny) - psi;
              if (!isfinite(fm)) {
                break;
              }
              if ((fm < 0.0) == (flo < 0.0)) {
                lo = mid;
                flo = fm;
              } else {
                hi = mid;
              }
            }
            double tb = 0.5 * (lo + hi);
            *r = rlin + tb * nx;
            *z = zlin + tb * ny;
            return true;
          }
          if (side == 0) {
            fprev_p = fk;
          } else {
            fprev_m = fk;
          }
        }
      }
    }
  }
  // Fourth resort: the requested level set does not exist anywhere near this
  // segment, so there is no crossing for any of the three searches above to
  // find and widening their windows cannot help.  psisep is psi evaluated at
  // whatever point the X-point finder returned, and when that point is not a
  // critical point of *this* (quadratic) representation the psi=psisep level
  // set is severed near the X point rather than crossing itself there.
  // Measured on 203997, whose quadratic rep has no in-vessel critical point at
  // all so the cubic location is substituted (|grad psi| there is 5.8e-4, not
  // ~0): {psi >= psisep} splits into two components separated by a 10.6 mm gap
  // straddling the R = 0.67875 cell face, where |grad psi| jumps 9x across the
  // C0 seam.  The trace's last segment spans that gap, so the level set is
  // absent over its whole length and psi stays below target across +/-20 mm of
  // normal offset.
  //
  // The two bracketing trace points remain the best available representatives
  // of the contour, so fall back to the linear interpolation between them --
  // which is exactly what the trace polyline already is between its stations.
  // This only runs where the routine previously returned false and aborted the
  // block, so it cannot move a node that any other path could place; and the
  // projection it replaces is a sub-micron refinement of the chord point
  // wherever the level set does exist.  The resulting grid is still subject to
  // the fold and surface-crossing checks, which is what decides whether the
  // interpolated node is acceptable.
  double residual = tok_eval_psi_rz_local(geo, rlin, zlin) - psi;
  if (isfinite(rlin) && isfinite(zlin) && isfinite(residual)) {
    fprintf(
      stderr,
      "TOK_LOGSAMP_CHORD_FALLBACK u=%.17g i=%d chord=(%.17g,%.17g) "
      "residual=%.17g seglen=%.17g psi=%.17g\n",
      u, i, rlin, zlin, residual, seglen, psi
    );
    *r = rlin;
    *z = zlin;
    return true;
  }
  fprintf(
    stderr,
    "TOK_LOGSAMP reason=no_projection axis_ok=%d u=%.17g i=%d w=%.17g "
    "chord=(%.17g,%.17g) p0=(%.17g,%.17g) p1=(%.17g,%.17g) seglen=%.17g "
    "psi=%.17g\n",
    (int)axis_ok, u, i, w, rlin, zlin, tr[i], tz[i], tr[i + 1], tz[i + 1], seglen, psi
  );
  return false;
}

// The plate a block end lies on, as its root finder resolves it: a declared
// plate function at the TOK_PLATE_NSAMP+1 points tok_plate_flux_intersection
// brackets on; an outline target's segments, in list order, at the step
// tok_divertor_wall_intersection scans them with. Returns the number of
// points (0 on failure); the caller frees *rz (R,Z pairs).
static int
tok_plate_polyline(const struct gkyl_tok_geo *geo, plate_func plate, double **rz)
{
  *rz = 0;
  const int slot = tok_divertor_wall_slot(geo, plate);
  if (slot < 0) {
    if (!plate) {
      return 0;
    }
    double *p = gkyl_malloc(2 * (TOK_PLATE_NSAMP + 1) * sizeof(double));
    for (int k = 0; k <= TOK_PLATE_NSAMP; ++k) {
      plate(k / (double)TOK_PLATE_NSAMP, p + 2 * k);
    }
    *rz = p;
    return TOK_PLATE_NSAMP + 1;
  }
  const struct gkyl_tok_geo_wall_target *t = &geo->divertor_wall[slot];
  const struct gkyl_efit *e = geo->efit;
  const int n = e->limiter_n, ns = t->num_segments;
  const struct gkyl_rect_grid *g = geo->use_cubics ? &geo->rzgrid_cubic : &geo->rzgrid;
  const double step = 0.1 * fmin(g->dx[0], g->dx[1]);
  if (ns < 1 || n < 3 || !(step > 0.0)) {
    return 0;
  }
  // Segment i joins limiter vertices i and i+1 and the list may run either
  // way, so the arc starts at the end of segment 0 that segment 1 does not share.
  int *v = gkyl_malloc((ns + 1) * sizeof(int));
  int a0 = t->segments[0], b0 = (a0 + 1) % n;
  if (ns > 1) {
    const int a1 = t->segments[1], b1 = (a1 + 1) % n;
    if (a0 == a1 || a0 == b1) {
      const int x = a0;
      a0 = b0;
      b0 = x;
    }
  }
  v[0] = a0;
  v[1] = b0;
  for (int k = 1; k < ns; ++k) {
    const int a = t->segments[k], b = (a + 1) % n;
    v[k + 1] = a == v[k] ? b : a;
  }
  int total = 1;
  for (int k = 0; k < ns; ++k) {
    const double c = ceil(
      hypot(
        e->limiter_R[v[k + 1]] - e->limiter_R[v[k]], e->limiter_Z[v[k + 1]] - e->limiter_Z[v[k]]
      ) /
      step
    );
    if (!(c <= 100000)) {
      gkyl_free(v);
      return 0;
    }
    total += GKYL_MAX2((int)c, 32);
  }
  double *p = gkyl_malloc(2 * total * sizeof(double));
  int m = 0;
  p[2 * m] = e->limiter_R[v[0]];
  p[2 * m + 1] = e->limiter_Z[v[0]];
  ++m;
  for (int k = 0; k < ns; ++k) {
    const double ra = e->limiter_R[v[k]], za = e->limiter_Z[v[k]];
    const double rb = e->limiter_R[v[k + 1]], zb = e->limiter_Z[v[k + 1]];
    const int c = GKYL_MAX2((int)ceil(hypot(rb - ra, zb - za) / step), 32);
    for (int j = 1; j <= c; ++j) {
      p[2 * m] = ra + (rb - ra) * j / (double)c;
      p[2 * m + 1] = za + (zb - za) * j / (double)c;
      ++m;
    }
  }
  gkyl_free(v);
  *rz = p;
  return m;
}

// Walk the plate polyline from a to b (both on it) and return the largest
// reversal of psi against the direction psi_a -> psi_b: zero when psi is
// monotone between the two strikes. `at` receives where the reversal is.
static double
tok_plate_flux_reversal(
  const struct gkyl_tok_geo *geo, const double *rz, int np, const double a[2], double psi_a,
  const double b[2], double psi_b, double at[2]
)
{
  int kab[2] = {0, 0};
  double tab[2] = {0.0, 0.0};
  for (int w = 0; w < 2; ++w) {
    const double *q = w ? b : a;
    double best = DBL_MAX;
    for (int k = 0; k + 1 < np; ++k) {
      const double dx = rz[2 * k + 2] - rz[2 * k], dy = rz[2 * k + 3] - rz[2 * k + 1],
                   l2 = dx * dx + dy * dy;
      const double t =
        l2 > 0.0 ?
          fmax(0.0, fmin(1.0, ((q[0] - rz[2 * k]) * dx + (q[1] - rz[2 * k + 1]) * dy) / l2)) :
          0.0;
      const double d = hypot(q[0] - rz[2 * k] - t * dx, q[1] - rz[2 * k + 1] - t * dy);
      if (d < best) {
        best = d;
        kab[w] = k;
        tab[w] = t;
      }
    }
  }
  const double sgn = psi_b >= psi_a ? 1.0 : -1.0;
  double run = -DBL_MAX, worst = 0.0;
  at[0] = a[0];
  at[1] = a[1];
  // a, then the polyline vertices strictly between, then b
  const bool fwd = kab[0] < kab[1] || (kab[0] == kab[1] && tab[0] <= tab[1]);
  const int first = fwd ? kab[0] + 1 : kab[0], last = fwd ? kab[1] : kab[1] + 1;
  const int nv = fwd ? GKYL_MAX2(last - first + 1, 0) : GKYL_MAX2(first - last + 1, 0);
  for (int i = -1; i <= nv; ++i) {
    double p[2];
    if (i < 0) {
      p[0] = a[0];
      p[1] = a[1];
    } else if (i == nv) {
      p[0] = b[0];
      p[1] = b[1];
    } else {
      const int k = fwd ? first + i : first - i;
      p[0] = rz[2 * k];
      p[1] = rz[2 * k + 1];
    }
    const double q = sgn * (tok_eval_psi_rz_local(geo, p[0], p[1]) - psi_a);
    run = fmax(run, q);
    if (run - q > worst) {
      worst = run - q;
      at[0] = p[0];
      at[1] = p[1];
    }
  }
  return worst;
}

// The strike rule (user decision 2026-10-06; replaces the 10-02 one-cell rule).
// Two adjacent rows of a plate-ended block end at two strike points; between
// them, psi along the plate must be MONOTONE. If it reverses, some flux
// surface between the two rows meets the plate more than once there -- it is
// tangent to the plate -- and the strike map folds: the cells between the rows
// cannot span it (NSTX-U 204951, fig24: the PF rows' strikes round the inner
// corner; psi along the centre-column wall reverses by 47% of the row spacing
// 35 mm above the corner). The tangency sits at a fixed flux, so some pair of
// neighbouring rows straddles it at every resolution, and the verdict does not
// depend on the cell sizes. The 10-02 rule compared the strike step with a
// theta cell -- a radial spacing against a poloidal one -- and so refused 7
// cells of the refinement matrix (mx12) where the strike merely slides along a
// grazing plate with psi strictly monotone (step halves with psi refinement).
// The tolerance is the plate root finders' own (1e-10 relative), the precision
// to which the strike points are known. Treated like a wall violation: inside
// a wall trial it is recorded, so the adjuster moves the requested PF (or SOL)
// bound inward until no pair of rows straddles the tangency; in a production
// build the row fails, loudly.
static bool
tok_ext_strike_step_check(
  const struct gkyl_tok_geo_grid_inp *inp, const struct arc_length_ctx *arc_ctx,
  const struct tok_ext_topology *top, double psi, bool at_sep, const double *raw_r,
  const double *raw_z, const double *raw_s, int raw_n
)
{
  (void)raw_s;
  const bool lo_plate = top->lower.kind == TOK_EXT_PLATE;
  const bool hi_plate = top->upper.kind == TOK_EXT_PLATE;
  if ((!lo_plate && !hi_plate) || raw_n < 2 || inp->cgrid.cells[0] < 1) {
    return true;
  }
  const struct gkyl_tok_geo *geo = arc_ctx->geo;
  const double dpsi = (inp->cgrid.upper[0] - inp->cgrid.lower[0]) / inp->cgrid.cells[0];
  const double psisep = geo->psisep;
  // The neighbouring row: toward the separatrix for an interior row, into the
  // block for the separatrix row itself.
  double psi_nb;
  if (at_sep) {
    psi_nb = psisep + copysign(dpsi, 0.5 * (inp->cgrid.lower[0] + inp->cgrid.upper[0]) - psisep);
  } else {
    psi_nb = psi + copysign(dpsi, psisep - psi);
  }
  for (int e = 0; e < 2; ++e) {
    if (!(e == 0 ? lo_plate : hi_plate)) {
      continue;
    }
    const struct tok_ext_endpoint *endpoint = e == 0 ? &top->lower : &top->upper;
    const double re = e == 0 ? raw_r[0] : raw_r[raw_n - 1];
    const double ze = e == 0 ? raw_z[0] : raw_z[raw_n - 1];
    double rn = 0.0, zn = 0.0;
    if (!tok_ext_endpoint_point(
          inp, arc_ctx, endpoint, psi_nb, tok_geo_same_flux(psi_nb, psisep), &rn, &zn, 0
        )) {
      continue; // the neighbour's own build will judge it
    }
    plate_func plate = endpoint->plate_slot == TOK_EXT_PLATE_LOWER ? geo->plate_func_lower :
                                                                     geo->plate_func_upper;
    double *rz = 0;
    const int np = tok_plate_polyline(geo, plate, &rz);
    if (np < 2) {
      if (rz) {
        gkyl_free(rz);
      }
      continue;
    }
    double at[2];
    const double rev = tok_plate_flux_reversal(
      geo, rz, np, (const double[2]){re, ze}, psi, (const double[2]){rn, zn}, psi_nb, at
    );
    gkyl_free(rz);
    if (!(rev > 1e-10 * fmax(1.0, fmax(fabs(psi), fabs(psi_nb))))) {
      continue;
    }
    const bool trial = tok_wall_trial_record_scope(false, false);
    fprintf(
      stderr,
      "TOK_STRIKE_FOLD ftype=%d psi=%.17g strike=(%.17g,%.17g) neighbour_psi=%.17g "
      "neighbour_strike=(%.17g,%.17g) reversal_psi=%.6g reversal_of_spacing=%.6g at=(%.17g,%.17g) trial=%d\n",
      inp->ftype, psi, re, ze, psi_nb, rn, zn, rev, rev / fabs(psi_nb - psi), at[0], at[1],
      (int)trial
    );
    if (!trial) {
      return false;
    }
  }
  return true;
}

static bool
tok_build_current_ordered_trace(
  const struct gkyl_tok_geo_grid_inp *inp, struct arc_length_ctx *arc_ctx
)
{
  double psi = arc_ctx->psi;
  if (arc_ctx->map_trace_initialized && tok_geo_same_flux(psi, arc_ctx->map_trace_psi)) {
    return true;
  }
  if (!tok_build_trace_correspondence(inp, arc_ctx)) {
    fprintf(
      stderr, "TOK_ORDERED_MAP_TRACE reason=correspondence_failed ftype=%d psi=%.17g\n", inp->ftype,
      psi
    );
    return false;
  }
  int n = tok_ext_map_trace_nodes(inp, arc_ctx->sep_trace_capacity);
  const bool extended = tok_ext_construction(inp);
  double *raw_r = 0, *raw_z = 0, *raw_s = 0;
  int raw_n = 0;
  bool raw_param_is_r = false, raw_closed = false;
  double radial_fraction = 0.0;
  if (extended) {
    int raw_capacity = arc_ctx->sep_trace_capacity;
    raw_r = gkyl_malloc(sizeof(double[raw_capacity]));
    raw_z = gkyl_malloc(sizeof(double[raw_capacity]));
    raw_s = gkyl_malloc(sizeof(double[raw_capacity]));
    bool at_sep = tok_geo_same_flux(psi, arc_ctx->geo->psisep);
    if (!tok_ext_build_domain_trace(
          inp, arc_ctx, psi, at_sep, raw_r, raw_z, raw_s, &raw_n, &raw_param_is_r, &raw_closed
        )) {
      fprintf(
        stderr, "TOK_ORDERED_MAP_TRACE reason=domain_trace_failed ftype=%d psi=%.17g at_sep=%d\n",
        inp->ftype, psi, (int)at_sep
      );
      gkyl_free(raw_r);
      gkyl_free(raw_z);
      gkyl_free(raw_s);
      return false;
    }
    struct tok_ext_topology top;
    bool have_top = tok_ext_topology_from_ftype(inp->ftype, inp->half_domain, &top);
    if (!have_top || raw_closed != top.closed) {
      fprintf(
        stderr,
        "TOK_ORDERED_MAP_TRACE reason=topology_mismatch ftype=%d psi=%.17g have_topology=%d raw_closed=%d expected_closed=%d\n",
        inp->ftype, psi, (int)have_top, (int)raw_closed, have_top ? (int)top.closed : -1
      );
      gkyl_free(raw_r);
      gkyl_free(raw_z);
      gkyl_free(raw_s);
      return false;
    }
    if (!tok_ext_strike_step_check(inp, arc_ctx, &top, psi, at_sep, raw_r, raw_z, raw_s, raw_n)) {
      gkyl_free(raw_r);
      gkyl_free(raw_z);
      gkyl_free(raw_s);
      return false;
    }
    double span = arc_ctx->xpt_ray_psi0 - arc_ctx->geo->psisep;
    if (!isfinite(span) ||
        fabs(span) <= 256.0 * DBL_EPSILON * fmax(1.0, fabs(arc_ctx->geo->psisep))) {
      fprintf(
        stderr,
        "TOK_ORDERED_MAP_TRACE reason=degenerate_radial_span ftype=%d psi=%.17g span=%.17g ray_psi0=%.17g psisep=%.17g\n",
        inp->ftype, psi, span, arc_ctx->xpt_ray_psi0, arc_ctx->geo->psisep
      );
      gkyl_free(raw_r);
      gkyl_free(raw_z);
      gkyl_free(raw_s);
      return false;
    }
    radial_fraction = (psi - arc_ctx->geo->psisep) / span;
    if (radial_fraction < -1e-8 || radial_fraction > 1.0 + 1e-8) {
      fprintf(
        stderr,
        "TOK_ORDERED_MAP_TRACE reason=radial_fraction_out_of_range ftype=%d psi=%.17g radial_fraction=%.17g span=%.17g psisep=%.17g\n",
        inp->ftype, psi, radial_fraction, span, arc_ctx->geo->psisep
      );
      gkyl_free(raw_r);
      gkyl_free(raw_z);
      gkyl_free(raw_s);
      return false;
    }
    radial_fraction = fmin(1.0, fmax(0.0, radial_fraction));
  }
  double fpol = tok_fpol_at_psi(arc_ctx->geo, psi);
  arc_ctx->map_trace_s[0] = 0.0;
  arc_ctx->map_trace_phi[0] = 0.0;
  bool at_sep = tok_geo_same_flux(psi, arc_ctx->geo->psisep);
  for (int i = 0; i < n; ++i) {
    double u = i / (double)(n - 1);
    bool point_ok = false;
    if (extended) {
      double w;
      if (arc_ctx->ext_ladder_initialized && n == arc_ctx->ext_ladder_n) {
        w = tok_ext_ladder_w(arc_ctx, radial_fraction, i);
      } else {
        // Two-point blend between the separatrix (w = u) and the far
        // boundary's correspondence.
        double v = tok_trace_correspondence(arc_ctx, u);
        w = (1.0 - radial_fraction) * u + radial_fraction * v;
      }
      point_ok = tok_trace_sample(
        arc_ctx->geo, psi, raw_r, raw_z, raw_s, raw_n, raw_param_is_r, w, &arc_ctx->map_trace_r[i],
        &arc_ctx->map_trace_z[i]
      );
      // Relaxed delta-s was originally wired only into the half-domain chord
      // point below. Apply the same bounded relaxation here so full-domain
      // (production) geometry can use it too; a no-op unless relaxed_xpt_seam
      // is set (see tok_relax_xpt_seam_point).
      if (point_ok) {
        point_ok = tok_relax_xpt_seam_point(
          inp, arc_ctx, psi, u, &arc_ctx->map_trace_r[i], &arc_ctx->map_trace_z[i]
        );
      }
    } else {
      point_ok = tok_parameterized_xpt_seam_point(
        inp, arc_ctx, psi, u, &arc_ctx->map_trace_r[i], &arc_ctx->map_trace_z[i]
      );
    }
    if (!point_ok) {
      fprintf(
        stderr,
        "TOK_ORDERED_MAP_TRACE reason=trace_point_failed ftype=%d psi=%.17g i=%d n=%d u=%.17g extended=%d radial_fraction=%.17g\n",
        inp->ftype, psi, i, n, u, (int)extended, radial_fraction
      );
      if (extended) {
        gkyl_free(raw_r);
        gkyl_free(raw_z);
        gkyl_free(raw_s);
      }
      return false;
    }
    if (i > 0) {
      double ds = hypot(
        arc_ctx->map_trace_r[i] - arc_ctx->map_trace_r[i - 1],
        arc_ctx->map_trace_z[i] - arc_ctx->map_trace_z[i - 1]
      );
      if (!(ds > 0.0) || !isfinite(ds)) {
        fprintf(
          stderr,
          "TOK_ORDERED_MAP_TRACE reason=nonpositive_arc_step ftype=%d psi=%.17g i=%d n=%d ds=%.17g prev=(%.17g,%.17g) curr=(%.17g,%.17g)\n",
          inp->ftype, psi, i, n, ds, arc_ctx->map_trace_r[i - 1], arc_ctx->map_trace_z[i - 1],
          arc_ctx->map_trace_r[i], arc_ctx->map_trace_z[i]
        );
        if (extended) {
          gkyl_free(raw_r);
          gkyl_free(raw_z);
          gkyl_free(raw_s);
        }
        return false;
      }
      arc_ctx->map_trace_s[i] = arc_ctx->map_trace_s[i - 1] + ds;
      double rm = 0.5 * (arc_ctx->map_trace_r[i] + arc_ctx->map_trace_r[i - 1]);
      double zm = 0.5 * (arc_ctx->map_trace_z[i] + arc_ctx->map_trace_z[i - 1]);
      double gr = 0.0, gz = 0.0;
      if (!tok_eval_psi_grad_rz_local(arc_ctx->geo, rm, zm, &gr, &gz)) {
        fprintf(
          stderr,
          "TOK_ORDERED_MAP_TRACE reason=psi_grad_eval_failed ftype=%d psi=%.17g i=%d mid=(%.17g,%.17g)\n",
          inp->ftype, psi, i, rm, zm
        );
        if (extended) {
          gkyl_free(raw_r);
          gkyl_free(raw_z);
          gkyl_free(raw_s);
        }
        return false;
      }
      double grad = hypot(gr, gz);
      if (!(grad > 1e-14) || !(rm > 0.0)) {
        fprintf(
          stderr,
          "TOK_ORDERED_MAP_TRACE reason=degenerate_grad_or_radius ftype=%d psi=%.17g i=%d grad=%.17g mid=(%.17g,%.17g)\n",
          inp->ftype, psi, i, grad, rm, zm
        );
        if (extended) {
          gkyl_free(raw_r);
          gkyl_free(raw_z);
          gkyl_free(raw_s);
        }
        return false;
      }
      arc_ctx->map_trace_phi[i] = arc_ctx->map_trace_phi[i - 1] + fpol * ds / (rm * grad);
    }
  }
  // (raw_r/raw_z/raw_s are freed after the diagnostic dump below.)
  // The angle on the final trace points, integrated on the contour between
  // them (see the true-arc angle section); a bracket that cannot be
  // integrated keeps its previous increment and is reported.
  double prev = arc_ctx->map_trace_phi[0];
  arc_ctx->map_trace_phi[0] = 0.0;
  for (int i = 1; i < n; ++i) {
    const double cur = arc_ctx->map_trace_phi[i];
    double d = 0.0;
    if (!tok_phi_integral(
          arc_ctx->geo, psi, fpol, arc_ctx->map_trace_r[i - 1], arc_ctx->map_trace_z[i - 1],
          arc_ctx->map_trace_r[i], arc_ctx->map_trace_z[i], &d
        )) {
      fprintf(
        stderr, "TOK_PHI_EXACT_FALLBACK stage=trace ftype=%d psi=%.17g i=%d n=%d\n", inp->ftype,
        psi, i, n
      );
      d = cur - prev;
    }
    arc_ctx->map_trace_phi[i] = arc_ctx->map_trace_phi[i - 1] + d;
    prev = cur;
  }
  if (extended) {
    gkyl_free(raw_r);
    gkyl_free(raw_z);
    gkyl_free(raw_s);
  }
  arc_ctx->map_trace_n = n;
  arc_ctx->map_trace_psi = psi;
  if (tok_ext_construction(inp) && !tok_ext_set_phi_reference(inp, arc_ctx)) {
    fprintf(
      stderr, "TOK_ORDERED_MAP_TRACE reason=phi_reference_failed ftype=%d psi=%.17g\n", inp->ftype,
      psi
    );
    return false;
  }
  arc_ctx->map_trace_initialized = true;
  return true;
}

static bool
tok_ordered_map_lookup(
  const struct gkyl_tok_geo_grid_inp *inp, struct arc_length_ctx *arc_ctx, double theta,
  double alpha, struct tok_ordered_point *out
)
{
  if (!tok_xpt_ordered_placement(inp)) {
    return false;
  }
  const struct gkyl_tok_geo_grid_inp *effective_inp = inp;
  struct gkyl_tok_geo_grid_inp straight_inp;
  if (!tok_build_current_ordered_trace(effective_inp, arc_ctx)) {
    if (tok_xpt_seam_optimizer_trial(inp)) {
      tok_xpt_seam_trial_reject(inp, GKYL_XPT_SEAM_TRIAL_TRACE_ORDERING);
      straight_inp = *inp;
      straight_inp.relaxed_xpt_seam_delta_s_coeff = 0.0;
      straight_inp.relaxed_xpt_seam_sweep = false;
      straight_inp.relaxed_xpt_seam_optimizer_trial = false;
      straight_inp.relaxed_xpt_seam_trial_status = 0;
      arc_ctx->map_trace_initialized = false;
      if (!tok_build_current_ordered_trace(&straight_inp, arc_ctx)) {
        fprintf(
          stderr, "TOK_ORDERED_MAP straight trial fallback failed ftype=%d psi=%.17g\n", inp->ftype,
          arc_ctx->psi
        );
        abort();
      }
      effective_inp = &straight_inp;
    } else {
      fprintf(
        stderr, "TOK_ORDERED_MAP initialization failed ftype=%d psi=%.17g\n", inp->ftype,
        arc_ctx->psi
      );
      abort();
    }
  }
  double dtheta = inp->cgrid.upper[2] - inp->cgrid.lower[2];
  double u = (theta - inp->cgrid.lower[2]) / dtheta;
  if (u < -1e-10 || u > 1.0 + 1e-10) {
    fprintf(
      stderr,
      "TOK_ORDERED_MAP_LOOKUP reason=u_out_of_range ftype=%d psi=%.17g theta=%.17g u=%.17g "
      "cgrid_ndim=%d cgrid_lo=(%.17g,%.17g,%.17g) cgrid_up=(%.17g,%.17g,%.17g) "
      "cgrid_cells=(%d,%d,%d) dtheta=%.17g\n",
      inp->ftype, arc_ctx->psi, theta, u, inp->cgrid.ndim, inp->cgrid.lower[0], inp->cgrid.lower[1],
      inp->cgrid.lower[2], inp->cgrid.upper[0], inp->cgrid.upper[1], inp->cgrid.upper[2],
      inp->cgrid.cells[0], inp->cgrid.cells[1], inp->cgrid.cells[2], dtheta
    );
    return false;
  }
  if (u <= 256.0 * DBL_EPSILON) {
    u = 0.0;
  } else if (u >= 1.0 - 256.0 * DBL_EPSILON) {
    u = 1.0;
  } else {
    u = fmin(1.0, fmax(0.0, u));
  }
  if (tok_ext_construction(effective_inp)) {
    if (!tok_logical_trace_sample(
          arc_ctx->geo, arc_ctx->psi, arc_ctx->map_trace_r, arc_ctx->map_trace_z,
          arc_ctx->map_trace_n, u, &out->r, &out->z
        )) {
      fprintf(
        stderr,
        "TOK_ORDERED_MAP_LOOKUP reason=logical_trace_sample_failed ftype=%d psi=%.17g u=%.17g n=%d\n",
        inp->ftype, arc_ctx->psi, u, arc_ctx->map_trace_n
      );
      return false;
    }
  } else if (!tok_parameterized_xpt_seam_point(
               effective_inp, arc_ctx, arc_ctx->psi, u, &out->r, &out->z
             )) {
    fprintf(
      stderr, "TOK_ORDERED_MAP_LOOKUP reason=seam_point_failed ftype=%d psi=%.17g u=%.17g\n",
      inp->ftype, arc_ctx->psi, u
    );
    return false;
  }
  double x = u * (arc_ctx->map_trace_n - 1);
  int i = GKYL_MIN2(arc_ctx->map_trace_n - 2, GKYL_MAX2(0, (int)floor(x)));
  double du = 1.0 / (arc_ctx->map_trace_n - 1);
  double dr = arc_ctx->map_trace_r[i + 1] - arc_ctx->map_trace_r[i];
  double dz = arc_ctx->map_trace_z[i + 1] - arc_ctx->map_trace_z[i];
  double speed_u = hypot(dr, dz) / du;
  double gr = 0.0, gz = 0.0;
  if (!tok_eval_psi_grad_rz_local(arc_ctx->geo, out->r, out->z, &gr, &gz)) {
    fprintf(
      stderr,
      "TOK_ORDERED_MAP_LOOKUP reason=grad_eval_failed ftype=%d psi=%.17g rz=(%.17g,%.17g)\n",
      inp->ftype, arc_ctx->psi, out->r, out->z
    );
    return false;
  }
  double tr = dr, tz = dz, grad = hypot(gr, gz);
  if (grad > 1e-14) {
    tr = -gz / grad;
    tz = gr / grad;
    if (tr * dr + tz * dz < 0.0) {
      tr = -tr;
      tz = -tz;
    }
  } else {
    double tmag = hypot(tr, tz);
    if (!(tmag > 0.0)) {
      fprintf(
        stderr,
        "TOK_ORDERED_MAP_LOOKUP reason=degenerate_tangent ftype=%d psi=%.17g u=%.17g grad=%.17g\n",
        inp->ftype, arc_ctx->psi, u, grad
      );
      return false;
    }
    tr /= tmag;
    tz /= tmag;
  }
  out->dr_dtheta = tr * speed_u / dtheta;
  out->dz_dtheta = tz * speed_u / dtheta;
  double w = x - i;
  double path_phi =
    arc_ctx->map_trace_phi[i] + w * (arc_ctx->map_trace_phi[i + 1] - arc_ctx->map_trace_phi[i]);
  double ref_phi = tok_ext_construction(effective_inp) ?
                     arc_ctx->map_trace_phi_ref :
                     (tok_sep_fixed_edge_is_first(effective_inp->ftype) ?
                        0.0 :
                        arc_ctx->map_trace_phi[arc_ctx->map_trace_n - 1]);
  out->dphi_dtheta = (arc_ctx->map_trace_phi[i + 1] - arc_ctx->map_trace_phi[i]) / du / dtheta;
  // The angle at the node itself: its bracket's start plus the contour from
  // there to the node (see the true-arc angle section), and its rate there.
  const double fpol = tok_fpol_at_psi(arc_ctx->geo, arc_ctx->psi);
  double d = 0.0;
  if (out->r == arc_ctx->map_trace_r[i + 1] && out->z == arc_ctx->map_trace_z[i + 1]) {
    path_phi = arc_ctx->map_trace_phi[i + 1];
  } else if (tok_phi_into_bracket(
               arc_ctx->geo, arc_ctx->psi, fpol, arc_ctx->map_trace_r[i], arc_ctx->map_trace_z[i],
               arc_ctx->map_trace_r[i + 1], arc_ctx->map_trace_z[i + 1],
               arc_ctx->map_trace_phi[i + 1] - arc_ctx->map_trace_phi[i], out->r, out->z, &d
             )) {
    path_phi = arc_ctx->map_trace_phi[i] + d;
  } else {
    fprintf(
      stderr, "TOK_PHI_EXACT_FALLBACK stage=node ftype=%d psi=%.17g u=%.17g i=%d\n", inp->ftype,
      arc_ctx->psi, u, i
    );
  }
  if (grad > 1e-14 && out->r > 0.0) {
    out->dphi_dtheta = fpol / (out->r * grad) * speed_u / dtheta;
  }
  out->phi = alpha + path_phi - ref_phi;
  return isfinite(out->r) && isfinite(out->z) && isfinite(out->phi) && isfinite(out->dr_dtheta) &&
         isfinite(out->dz_dtheta) && isfinite(out->dphi_dtheta);
}

// Function to calculate phi given alpha
double
phi_func(double alpha_curr, double Z, void *ctx)
{
  struct arc_length_ctx *actx = ctx;
  double *arc_memo = actx->arc_memo;
  double psi = actx->psi, rclose = actx->rclose, zmin = actx->zmin, arcL = actx->arcL,
         zmax = actx->zmax;

  // Here we will abandon conventions about alpha and phi except for full core and full SN cases
  // The convention for phi only affects b_x - it does not affect any quantities used in axisymmetric simulations
  // I have not quite figured out full 3D yet. b_x presents a serious problem as of now. Akash Shukla 1/20/2024
  // The idea for axisymmetry is that I am avoiding starting integrals at the x-point to minimize issues
  double ival = 0;
  double phi_ref = 0.0;
  if (actx->ftype == GKYL_GEOMETRY_TOKAMAK_CORE) { // phi = alpha at outboard midplane
    if (actx->right == true) {
      if (Z < actx->zmaxis) {
        ival = -integrate_phi_along_psi_contour_memo(
          actx->geo, psi, Z, actx->zmaxis, rclose, false, false, arc_memo
        );
      } else {
        ival = integrate_phi_along_psi_contour_memo(
          actx->geo, psi, actx->zmaxis, Z, rclose, false, false, arc_memo
        );
      }
    } else {
      ival = integrate_phi_along_psi_contour_memo(
        actx->geo, psi, Z, actx->zmax, rclose, false, false, arc_memo
      );
      phi_ref = actx->phi_right;
    }
  } else if (actx->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_L) { // alpha = phi at inboard midplane
    if (Z < actx->zmaxis) {
      ival = integrate_phi_along_psi_contour_memo(
        actx->geo, psi, Z, actx->zmaxis, rclose, false, false, arc_memo
      );
    } else {
      ival = -integrate_phi_along_psi_contour_memo(
        actx->geo, psi, actx->zmaxis, Z, rclose, false, false, arc_memo
      );
    }
  }

  else if (actx->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R) { // alpha = phi at outboard midplane
    if (Z < actx->zmaxis) {
      ival = -integrate_phi_along_psi_contour_memo(
        actx->geo, psi, Z, actx->zmaxis, rclose, false, false, arc_memo
      );
    } else {
      ival = integrate_phi_along_psi_contour_memo(
        actx->geo, psi, actx->zmaxis, Z, rclose, false, false, arc_memo
      );
    }
  }

  else if ((actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT) ||
           (actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID
           )) { // alpha = phi at outboard midplane
    if (Z < actx->zmaxis) {
      ival = -integrate_phi_along_psi_contour_memo(
        actx->geo, psi, Z, actx->zmaxis, rclose, false, false, arc_memo
      );
    } else {
      ival = integrate_phi_along_psi_contour_memo(
        actx->geo, psi, actx->zmaxis, Z, rclose, false, false, arc_memo
      );
    }
  } else if (actx->ftype ==
             GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO) { // alpha = phi at lower plate and increases towards xpt
    ival =
      integrate_phi_along_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, false, false, arc_memo);
  } else if (actx->ftype ==
             GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP) { //alpha = phi at upper plate and decreases towards xpt
    ival = -integrate_phi_along_psi_contour_memo(
      actx->geo, psi, Z, zmax, rclose, false, false, arc_memo
    );
  }
  if ((actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN) ||
      (actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID)) { // alpha = phi at inboard midplane
    if (Z < actx->zmaxis) {
      ival = integrate_phi_along_psi_contour_memo(
        actx->geo, psi, Z, actx->zmaxis, rclose, false, false, arc_memo
      );
    } else {
      ival = -integrate_phi_along_psi_contour_memo(
        actx->geo, psi, actx->zmaxis, Z, rclose, false, false, arc_memo
      );
    }
  } else if (actx->ftype ==
             GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO) { // alpha = phi at lower plate and decreases towards xpt
    ival = -integrate_phi_along_psi_contour_memo(
      actx->geo, psi, zmin, Z, rclose, false, false, arc_memo
    );
  } else if (actx->ftype ==
             GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_UP) { // alpha = phi at upper plate and increases towards xpt
    ival =
      integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, false, false, arc_memo);
  } else if (actx->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL ||
             actx->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID) {
    // alpha = phi at outboard midplane
    if (actx->right == true) {
      if (Z < actx->zmaxis) {
        ival = -integrate_phi_along_psi_contour_memo(
          actx->geo, psi, Z, actx->zmaxis, rclose, false, false, arc_memo
        );
      } else {
        ival = integrate_phi_along_psi_contour_memo(
          actx->geo, psi, actx->zmaxis, Z, rclose, false, false, arc_memo
        );
      }
    } else {
      ival = integrate_phi_along_psi_contour_memo(
        actx->geo, psi, Z, actx->zmax, rclose, false, false, arc_memo
      );
      phi_ref = actx->phi_right;
    }
  } else if (actx->ftype ==
             GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO) { //alpha = phi at outer plate and increases towards xpt
    ival =
      integrate_phi_along_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, false, false, arc_memo);
  } else if (actx->ftype ==
             GKYL_GEOMETRY_TOKAMAK_LSN_SOL_UP) { //alpha = phi at inner plate and decreases towards xpt
    ival = -integrate_phi_along_psi_contour_memo(
      actx->geo, psi, zmin, Z, rclose, false, false, arc_memo
    );
  } else if (actx->ftype ==
             GKYL_GEOMETRY_TOKAMAK_PF_LO_R) { // alpha = phi at outer plate and increases towards xpt
    ival =
      integrate_phi_along_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, false, false, arc_memo);
  } else if (actx->ftype ==
             GKYL_GEOMETRY_TOKAMAK_PF_LO_L) { //alpha = phi at inner plate and decreases towards xpt
    ival = -integrate_phi_along_psi_contour_memo(
      actx->geo, psi, zmin, Z, rclose, false, false, arc_memo
    );
  } else if (actx->ftype ==
             GKYL_GEOMETRY_TOKAMAK_PF_UP_R) { // alpha = phi at outer plate and decreases towards Xpt
    ival = -integrate_phi_along_psi_contour_memo(
      actx->geo, psi, Z, zmax, rclose, false, false, arc_memo
    );
  } else if (actx->ftype ==
             GKYL_GEOMETRY_TOKAMAK_PF_UP_L) { // alpha = phi at inner plate and increases towards xpt
    ival =
      integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, false, false, arc_memo);
  } else if (actx->ftype == GKYL_GEOMETRY_TOKAMAK_IWL) {
    // phi = alpha at outboard midplane
    if (actx->right == true) {
      if (Z < actx->zmaxis) {
        ival = -integrate_phi_along_psi_contour_memo(
          actx->geo, psi, Z, actx->zmaxis, rclose, false, false, arc_memo
        );
      } else {
        ival = integrate_phi_along_psi_contour_memo(
          actx->geo, psi, actx->zmaxis, Z, rclose, false, false, arc_memo
        );
      }
    } else {
      if (Z < actx->zmaxis) {
        ival = -integrate_phi_along_psi_contour_memo(
          actx->geo, psi, actx->zmin, Z, rclose, false, false, arc_memo
        );
        phi_ref = -actx->phi_right;
      } else {
        ival = integrate_phi_along_psi_contour_memo(
          actx->geo, psi, Z, actx->zmax, rclose, false, false, arc_memo
        );
        phi_ref = actx->phi_right;
      }
    }
  }
  // Now multiply by fpol
  double R[4] = {0};
  double dRdZ[4] = {0};
  double dR[4] = {0};
  double dZ[4] = {0};
  int nr = gkyl_tok_geo_R_psiZ(actx->geo, psi, Z, 4, R, dRdZ, dR, dZ);
  double r_curr = nr == 1 ? R[0] : choose_closest(rclose, R, R, nr);
  double psi_fpol = psi;
  if ((psi_fpol < actx->geo->fgrid.lower[0]) ||
      (psi_fpol > actx->geo->fgrid.upper[0])) { // F = F(psi_sep) in the SOL.
    psi_fpol = actx->geo->sibry;
  }
  int idx = fmin(
    actx->geo->frange.lower[0] +
      (int)floor((psi_fpol - actx->geo->fgrid.lower[0]) / actx->geo->fgrid.dx[0]),
    actx->geo->frange.upper[0]
  );
  long loc = gkyl_range_idx(&actx->geo->frange, &idx);
  const double *coeffs = gkyl_array_cfetch(actx->geo->fpoldg, loc);
  double fxc;
  gkyl_rect_grid_cell_center(&actx->geo->fgrid, &idx, &fxc);
  double fx = (psi_fpol - fxc) / (actx->geo->fgrid.dx[0] * 0.5);
  double fpol = actx->geo->fbasis.eval_expand(&fx, coeffs);
  ival = ival * fpol;

  //while(ival < -M_PI){
  //  ival +=2*M_PI;
  //}
  //while(ival > M_PI){
  //  ival -=2*M_PI;
  //}
  return alpha_curr + ival + phi_ref;
}

double
qprofile_func(void *ctx)
{
  // Function to calculate the flux surface averaged q profile.

  struct arc_length_ctx *actx = ctx;
  double *arc_memo = actx->arc_memo;
  double psi = actx->psi, rclose = actx->rclose, zmin = actx->zmin, arcL = actx->arcL,
         zmax = actx->zmax;
  double rleft = actx->rleft, rright = actx->rright;

  // Calculate q(psi) = -F(psi)/2pi * integral_zmin^zmax 1/Rgrad(psi).

  double ival = 0;
  double phi_ref = 0.0;
  if (actx->ftype == GKYL_GEOMETRY_TOKAMAK_CORE) {
    double ival1 = integrate_phi_along_psi_contour_memo(
      actx->geo, psi, actx->zmin, actx->zmax, rright, false, false, arc_memo
    );
    double ival2 = -integrate_phi_along_psi_contour_memo(
      actx->geo, psi, actx->zmin, actx->zmax, rleft, false, false, arc_memo
    );
    ival = ival1 + ival2;
  }

  if (actx->ftype == GKYL_GEOMETRY_TOKAMAK_IWL) {
    double ival1 = integrate_phi_along_psi_contour_memo(
      actx->geo, psi, actx->zmin, actx->zmax, rright, false, false, arc_memo
    );
    double ival2 = -integrate_phi_along_psi_contour_memo(
      actx->geo, psi, actx->zmin, actx->zmax, rleft, false, false, arc_memo
    );
    ival = ival1 + ival2;
  }

  if (actx->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL) {
    double ival1 = integrate_phi_along_psi_contour_memo(
      actx->geo, psi, actx->zmin_right, actx->zmax, rright, false, false, arc_memo
    );
    double ival2 = -integrate_phi_along_psi_contour_memo(
      actx->geo, psi, actx->zmin_left, actx->zmax, rleft, false, false, arc_memo
    );
    ival = ival1 + ival2;
  }

  if (actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT ||
      actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID ||
      actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO ||
      actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP) {
    ival = integrate_phi_along_psi_contour_memo(
      actx->geo, psi, actx->zmin, actx->zmax, rclose, false, false, arc_memo
    );
  }

  if (actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN ||
      actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID ||
      actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO ||
      actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_UP) {
    ival = -integrate_phi_along_psi_contour_memo(
      actx->geo, psi, actx->zmin, actx->zmax, rclose, false, false, arc_memo
    );
  }

  // Now multiply by fpol/2pi.
  double R[4] = {0};
  double dR[4] = {0};
  double psi_fpol = psi;
  if ((psi_fpol < actx->geo->fgrid.lower[0]) ||
      (psi_fpol > actx->geo->fgrid.upper[0])) { // F = F(psi_sep) in the SOL.
    psi_fpol = actx->geo->sibry;
  }
  int idx = fmin(
    actx->geo->frange.lower[0] +
      (int)floor((psi_fpol - actx->geo->fgrid.lower[0]) / actx->geo->fgrid.dx[0]),
    actx->geo->frange.upper[0]
  );
  long loc = gkyl_range_idx(&actx->geo->frange, &idx);
  const double *coeffs = gkyl_array_cfetch(actx->geo->fpoldg, loc);
  double fxc;
  gkyl_rect_grid_cell_center(&actx->geo->fgrid, &idx, &fxc);
  double fx = (psi_fpol - fxc) / (actx->geo->fgrid.dx[0] * 0.5);
  double fpol = actx->geo->fbasis.eval_expand(&fx, coeffs);
  double qout = -ival * fpol / M_PI;

  // AS 1/15/25: The 3 lines below are a useful check to compare against q from efit.
  //coeffs = gkyl_array_cfetch(actx->geo->qdg,loc);
  //double q_efit = actx->geo->fbasis.eval_expand(&fx, coeffs);
  // printf("psi_curr = %g, my q = %g, efit q = %g\n", psi_fpol, qout, q_efit);

  return qout;
}

static double
dphidtheta_func(double Z, void *ctx)
{
  struct arc_length_ctx *actx = ctx;
  double *arc_memo = actx->arc_memo;
  double psi = actx->psi, rclose = actx->rclose, zmin = actx->zmin, arcL = actx->arcL,
         zmax = actx->zmax;

  // Get the integrand
  double integrand = 0.0;
  struct contour_ctx cctx = {.geo = actx->geo, .psi = psi, .ncall = 0, .last_R = rclose};
  integrand = dphidtheta_integrand(Z, &cctx);
  // Now multiply by fpol
  double R[4] = {0};
  double dRdZ[4] = {0};
  double dR[4] = {0};
  double dZ[4] = {0};
  int nr = gkyl_tok_geo_R_psiZ(actx->geo, psi, Z, 4, R, dRdZ, dR, dZ);
  double r_curr = nr == 1 ? R[0] : choose_closest(rclose, R, R, nr);
  double psi_fpol = psi;
  if ((psi_fpol < actx->geo->fgrid.lower[0]) ||
      (psi_fpol > actx->geo->fgrid.upper[0])) { // F = F(psi_sep) in the SOL.
    psi_fpol = actx->geo->sibry;
  }
  int idx = fmin(
    actx->geo->frange.lower[0] +
      (int)floor((psi_fpol - actx->geo->fgrid.lower[0]) / actx->geo->fgrid.dx[0]),
    actx->geo->frange.upper[0]
  );
  long loc = gkyl_range_idx(&actx->geo->frange, &idx);
  const double *coeffs = gkyl_array_cfetch(actx->geo->fpoldg, loc);
  double fxc;
  gkyl_rect_grid_cell_center(&actx->geo->fgrid, &idx, &fxc);
  double fx = (psi_fpol - fxc) / (actx->geo->fgrid.dx[0] * 0.5);
  double fpol = actx->geo->fbasis.eval_expand(&fx, coeffs);
  integrand = integrand * fpol;
  double darc_dtheta =
    actx->arc_interval_valid ?
      actx->arc_darc_dtheta :
      (actx->xpt_map_valid ? actx->xpt_map_darc_dtheta : actx->arcL_tot / (2.0 * M_PI));
  integrand = integrand * darc_dtheta;
  return integrand;
}

static double
bmag_func(double r_curr, double Z, void *ctx)
{
  struct arc_length_ctx *actx = ctx;
  double *arc_memo = actx->arc_memo;
  double psi = actx->psi, rclose = actx->rclose, zmin = actx->zmin, arcL = actx->arcL,
         zmax = actx->zmax;
  // Calculate fpol
  double psi_fpol = psi;
  if ((psi_fpol < actx->geo->fgrid.lower[0]) ||
      (psi_fpol > actx->geo->fgrid.upper[0])) { // F = F(psi_sep) in the SOL.
    psi_fpol = actx->geo->sibry;
  }
  int idx = fmin(
    actx->geo->frange.lower[0] +
      (int)floor((psi_fpol - actx->geo->fgrid.lower[0]) / actx->geo->fgrid.dx[0]),
    actx->geo->frange.upper[0]
  );
  long loc = gkyl_range_idx(&actx->geo->frange, &idx);
  const double *coeffs = gkyl_array_cfetch(actx->geo->fpoldg, loc);
  double fxc;
  gkyl_rect_grid_cell_center(&actx->geo->fgrid, &idx, &fxc);
  double fx = (psi_fpol - fxc) / (actx->geo->fgrid.dx[0] * 0.5);
  double fpol = actx->geo->fbasis.eval_expand(&fx, coeffs);
  double Bphi = fpol / r_curr;
  double Br = 0.0, Bz = 0.0, bmag = 0.0;

  if (actx->geo->use_cubics) {
    double xn[2] = {r_curr, Z};
    double fout[3];
    actx->geo->efit->evf->eval_cubic_wgrad(0.0, xn, fout, actx->geo->efit->evf->ctx);
    double dpsidR = fout[1];
    double dpsidZ = fout[2];

    Br = 1.0 / r_curr * dpsidZ;
    Bz = -1.0 / r_curr * dpsidR;
  } else {
    int rzidx[2];
    int idxtemp = actx->geo->rzlocal.lower[0] +
                  (int)floor((r_curr - actx->geo->rzgrid.lower[0]) / actx->geo->rzgrid.dx[0]);
    idxtemp = GKYL_MIN2(idxtemp, actx->geo->rzlocal.upper[0]);
    idxtemp = GKYL_MAX2(idxtemp, actx->geo->rzlocal.lower[0]);
    rzidx[0] = idxtemp;
    idxtemp = actx->geo->rzlocal.lower[1] +
              (int)floor((Z - actx->geo->rzgrid.lower[1]) / actx->geo->rzgrid.dx[1]);
    idxtemp = GKYL_MIN2(idxtemp, actx->geo->rzlocal.upper[1]);
    idxtemp = GKYL_MAX2(idxtemp, actx->geo->rzlocal.lower[1]);
    rzidx[1] = idxtemp;

    long loc = gkyl_range_idx((&actx->geo->rzlocal), rzidx);
    const double *psih = gkyl_array_cfetch(actx->geo->psiRZ, loc);

    double xc[2];
    gkyl_rect_grid_cell_center((&actx->geo->rzgrid), rzidx, xc);
    double x = (r_curr - xc[0]) / (actx->geo->rzgrid.dx[0] * 0.5);
    double y = (Z - xc[1]) / (actx->geo->rzgrid.dx[1] * 0.5);

    double dpsidx = 5.625 * psih[8] * (2.0 * x * SQ(y) - 0.6666666666666666 * x) +
                    2.904737509655563 * psih[7] * (SQ(y) - 0.3333333333333333) +
                    5.809475019311126 * psih[6] * x * y + 1.5 * psih[3] * y +
                    3.354101966249684 * psih[4] * x + 0.8660254037844386 * psih[1];
    double dpsidy = 5.625 * psih[8] * (2.0 * SQ(x) * y - 0.6666666666666666 * y) +
                    5.809475019311126 * psih[7] * x * y + 3.354101966249684 * psih[5] * y +
                    2.904737509655563 * psih[6] * (SQ(x) - 0.3333333333333333) + 1.5 * psih[3] * x +
                    0.8660254037844386 * psih[2];
    double dpsidR = dpsidx * 2.0 / actx->geo->rzgrid.dx[0];
    double dpsidZ = dpsidy * 2.0 / actx->geo->rzgrid.dx[1];
    Br = 1.0 / r_curr * dpsidZ;
    Bz = -1.0 / r_curr * dpsidR;
  }
  bmag = sqrt(Br * Br + Bz * Bz + Bphi * Bphi);
  return bmag;
}

static void
curlbhat_func(double psi, double r_curr, double Z, double phi, double *curlbhat, void *ctx)
{
  struct arc_length_ctx *actx = ctx;
  double *arc_memo = actx->arc_memo;
  // Calculate fpol and fpolprime
  // First get the location on the flux grid
  double psi_fpol = psi;
  if ((psi_fpol < actx->geo->fgrid.lower[0]) ||
      (psi_fpol > actx->geo->fgrid.upper[0])) { // F = F(psi_sep) in the SOL.
    psi_fpol = actx->geo->sibry;
  }
  int idx = fmin(
    actx->geo->frange.lower[0] +
      (int)floor((psi_fpol - actx->geo->fgrid.lower[0]) / actx->geo->fgrid.dx[0]),
    actx->geo->frange.upper[0]
  );
  long loc = gkyl_range_idx(&actx->geo->frange, &idx);
  const double *coeffs = gkyl_array_cfetch(actx->geo->fpoldg, loc);
  double fxc;
  gkyl_rect_grid_cell_center(&actx->geo->fgrid, &idx, &fxc);
  double fx = (psi_fpol - fxc) / (actx->geo->fgrid.dx[0] * 0.5);
  // Second calculate fpol and bphi
  double fpol = actx->geo->fbasis.eval_expand(&fx, coeffs);
  double Bphi = fpol / r_curr;
  // Third calculate fpolprime
  coeffs = gkyl_array_cfetch(actx->geo->fpolprimedg, loc);
  double fpolprime = actx->geo->fbasis.eval_expand(&fx, coeffs);

  // Now calculate psi's various derivatives
  double Br = 0.0, Bz = 0.0, bmag = 0.0;
  double dpsidR = 0.0, dpsidZ = 0.0;
  double d2psidR2 = 0.0, d2psidZ2 = 0.0, d2psidRdZ = 0.0;
  double dBdR = 0.0, dBdZ = 0.0;
  double dBrdR = 0.0, dBrdZ = 0.0;
  double dBzdR = 0.0, dBzdZ = 0.0;

  if (actx->geo->use_cubics) {
    double xn[2] = {r_curr, Z};
    double fout[4];
    actx->geo->efit->evf->eval_cubic_wgrad(0.0, xn, fout, actx->geo->efit->evf->ctx);
    dpsidR = fout[1];
    dpsidZ = fout[2];
    actx->geo->efit->evf->eval_cubic_wgrad2(0.0, xn, fout, actx->geo->efit->evf->ctx);
    d2psidR2 = fout[1];
    d2psidZ2 = fout[2];
    d2psidRdZ = fout[3];
  } else {
    int rzidx[2];
    int idxtemp = actx->geo->rzlocal.lower[0] +
                  (int)floor((r_curr - actx->geo->rzgrid.lower[0]) / actx->geo->rzgrid.dx[0]);
    idxtemp = GKYL_MIN2(idxtemp, actx->geo->rzlocal.upper[0]);
    idxtemp = GKYL_MAX2(idxtemp, actx->geo->rzlocal.lower[0]);
    rzidx[0] = idxtemp;
    idxtemp = actx->geo->rzlocal.lower[1] +
              (int)floor((Z - actx->geo->rzgrid.lower[1]) / actx->geo->rzgrid.dx[1]);
    idxtemp = GKYL_MIN2(idxtemp, actx->geo->rzlocal.upper[1]);
    idxtemp = GKYL_MAX2(idxtemp, actx->geo->rzlocal.lower[1]);
    rzidx[1] = idxtemp;

    long loc = gkyl_range_idx((&actx->geo->rzlocal), rzidx);
    const double *psih = gkyl_array_cfetch(actx->geo->psiRZ, loc);

    double xc[2];
    gkyl_rect_grid_cell_center((&actx->geo->rzgrid), rzidx, xc);
    double x = (r_curr - xc[0]) / (actx->geo->rzgrid.dx[0] * 0.5);
    double y = (Z - xc[1]) / (actx->geo->rzgrid.dx[1] * 0.5);

    double dpsidx = 5.625 * psih[8] * (2.0 * x * SQ(y) - 0.6666666666666666 * x) +
                    2.904737509655563 * psih[7] * (SQ(y) - 0.3333333333333333) +
                    5.809475019311126 * psih[6] * x * y + 1.5 * psih[3] * y +
                    3.354101966249684 * psih[4] * x + 0.8660254037844386 * psih[1];
    double dpsidy = 5.625 * psih[8] * (2.0 * SQ(x) * y - 0.6666666666666666 * y) +
                    5.809475019311126 * psih[7] * x * y + 3.354101966249684 * psih[5] * y +
                    2.904737509655563 * psih[6] * (SQ(x) - 0.3333333333333333) + 1.5 * psih[3] * x +
                    0.8660254037844386 * psih[2];
    dpsidR = dpsidx * 2.0 / actx->geo->rzgrid.dx[0];
    dpsidZ = dpsidy * 2.0 / actx->geo->rzgrid.dx[1];
    double d2psidx2 = 11.25 * psih[8] * y * y + 5.809475019311125 * psih[6] * y - 3.75 * psih[8] +
                      3.354101966249685 * psih[4];
    double d2psidy2 = 11.25 * psih[8] * x * x + 5.809475019311125 * psih[7] * x - 3.75 * psih[8] +
                      3.354101966249685 * psih[5];
    double d2psidxdy = 22.5 * psih[8] * x * y + 5.809475019311125 * psih[7] * y +
                       5.809475019311125 * psih[6] * x + 1.5 * psih[3];
    d2psidR2 = d2psidx2 * (2.0 / actx->geo->rzgrid.dx[0]) * (2.0 / actx->geo->rzgrid.dx[0]);
    d2psidZ2 = d2psidy2 * (2.0 / actx->geo->rzgrid.dx[1]) * (2.0 / actx->geo->rzgrid.dx[1]);
    d2psidRdZ = d2psidxdy * (2.0 / actx->geo->rzgrid.dx[0]) * (2.0 / actx->geo->rzgrid.dx[1]);
  }

  Br = 1.0 / r_curr * dpsidZ;
  Bz = -1.0 / r_curr * dpsidR;
  bmag = sqrt(Br * Br + Bz * Bz + Bphi * Bphi);

  dBrdR = 1.0 / r_curr * d2psidRdZ;
  dBrdZ = 1.0 / r_curr * d2psidZ2;
  dBzdR = -1.0 / r_curr * d2psidR2;
  dBzdZ = -1.0 / r_curr * d2psidRdZ;

  double dFdR = fpolprime * dpsidR;
  double dFdZ = fpolprime * dpsidZ;
  dBdR =
    1 / bmag * (Br * dBrdR + Bz * dBzdR + fpol / r_curr * (dFdR / r_curr - fpol / r_curr / r_curr));
  dBdZ = 1 / bmag * (Br * dBrdZ + Bz * dBzdZ + fpol / r_curr * dFdZ);

  // Get the polar components (contravariant, upperscript components on tangent basis)
  double polar_comp[3] = {0.0};
  polar_comp[0] = 1.0 / bmag * -1.0 / r_curr * dFdZ -
                  1.0 / bmag * 1.0 / bmag * dBdZ * fpol / r_curr; // R component ^1
  polar_comp[1] = 1.0 / bmag * 1.0 / r_curr * (dBrdZ - dBzdR) +
                  (-dBdR * Bz / r_curr + dBdZ * Br / r_curr); // Phi component ^2
  polar_comp[2] = 1.0 / bmag * 1.0 / r_curr * dFdR + 1.0 / bmag * 1.0 / bmag * dBdR * fpol / r_curr;

  // Convert to cartesian
  curlbhat[0] = polar_comp[0] * cos(phi) - polar_comp[1] * sin(phi) * r_curr;
  curlbhat[1] = polar_comp[0] * sin(phi) + polar_comp[1] * cos(phi) * r_curr;
  curlbhat[2] = polar_comp[2];
}

// Build a conservative enclosure of psi over every RZ cell.
//
// R_psiZ finds the R values where psi(R,Z) = psi0 along one Z row by solving a
// small polynomial in each R cell of that row. The scan is O(cells in the row)
// and does not stop early, because a row normally holds fewer crossings than
// the caller's root budget -- so every cell is solved even though psi in almost
// all of them is nowhere near psi0. The cost therefore tracks the equilibrium
// file's radial resolution rather than the grid being built, and on a finely
// resolved equilibrium that scan dominates the whole geometry build.
//
// psi inside a cell is the modal expansion sum_i c_i b_i with b_0 constant, so
//
//   psi(x) in [c_0 b_0 - S, c_0 b_0 + S],   S = sum_{i>0} |c_i| max|b_i|,
//
// which encloses psi by the triangle inequality. A cell whose enclosure
// excludes psi0 cannot contain a root, so it is skipped without being solved.
// The enclosure is widened before use, keeping the filter conservative: a cell
// wrongly kept only costs the solve that would have happened anyway, whereas a
// cell wrongly dropped would lose a root.
static void
tok_psi_enclosure(
  const struct gkyl_basis *basis, const struct gkyl_array *psi, const struct gkyl_range *local,
  struct gkyl_array **cell_out, double **block_out, int *bsz_out, int *nblk_out
)
{
  int nb = basis->num_basis;
  double *bsup = gkyl_malloc(nb * sizeof(double));
  double *bval = gkyl_malloc(nb * sizeof(double));
  for (int k = 0; k < nb; ++k) {
    bsup[k] = 0.0;
  }

  // Sup norm of each basis function over the reference cell. Each is a
  // polynomial of degree poly_order per direction, so a sweep scaled to that
  // degree resolves its extrema; the widening applied below absorbs whatever
  // the sweep misses.
  const int nsamp = 64 * (int)basis->poly_order + 1;
  for (int i = 0; i < nsamp; ++i) {
    for (int j = 0; j < nsamp; ++j) {
      double z[2] = {-1.0 + 2.0 * i / (nsamp - 1.0), -1.0 + 2.0 * j / (nsamp - 1.0)};
      basis->eval(z, bval);
      for (int k = 0; k < nb; ++k) {
        bsup[k] = fmax(bsup[k], fabs(bval[k]));
      }
    }
  }

  // The mean term is separated out only if b_0 really is constant; otherwise
  // fall back to an enclosure centred on zero, which is still valid.
  double zc[2] = {0.0, 0.0};
  basis->eval(zc, bval);
  double b0 = bval[0];
  bool b0_const = fabs(bsup[0] - fabs(b0)) <= 1.0e-12 * fmax(1.0, bsup[0]);

  struct gkyl_array *bounds = gkyl_array_new(GKYL_DOUBLE, 2, psi->size);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(local, iter.idx);
    const double *c = gkyl_array_cfetch(psi, loc);

    double mean = b0_const ? c[0] * b0 : 0.0;
    double spread = 0.0;
    for (int k = b0_const ? 1 : 0; k < nb; ++k) {
      spread += fabs(c[k]) * bsup[k];
    }

    // Widen. The root solve tolerates a marginally negative discriminant, so a
    // root can sit just outside the exact range of psi in the cell.
    spread = 2.0 * spread + 1.0e-12 * fabs(mean);

    double *b = gkyl_array_fetch(bounds, loc);
    b[0] = mean - spread;
    b[1] = mean + spread;
  }

  *cell_out = bounds;

  // Coarse level: enclose runs of neighbouring R cells so the scan can reject a
  // whole run with a single test. A two-level scan costs about nr/b + b tests
  // per row, which is smallest at b = sqrt(nr), so the run length follows the
  // row rather than being fixed.
  int rlo = local->lower[0], rup = local->upper[0];
  int zlo = local->lower[1], zup = local->upper[1];
  int nr = rup - rlo + 1, nz = zup - zlo + 1;
  int bsz = (int)(sqrt((double)nr) + 0.5);
  if (bsz < 1) {
    bsz = 1;
  }
  int nblk = (nr + bsz - 1) / bsz;

  double *blocks = gkyl_malloc(sizeof(double[2]) * (size_t)nz * nblk);
  for (int iz = zlo; iz <= zup; ++iz) {
    for (int ib = 0; ib < nblk; ++ib) {
      double lo = DBL_MAX, hi = -DBL_MAX;
      int i0 = rlo + ib * bsz, i1 = i0 + bsz - 1;
      if (i1 > rup) {
        i1 = rup;
      }
      for (int ir = i0; ir <= i1; ++ir) {
        int cidx[2] = {ir, iz};
        const double *cb = gkyl_array_cfetch(bounds, gkyl_range_idx(local, cidx));
        lo = fmin(lo, cb[0]);
        hi = fmax(hi, cb[1]);
      }
      double *b = blocks + 2 * ((size_t)(iz - zlo) * nblk + ib);
      b[0] = lo;
      b[1] = hi;
    }
  }
  *block_out = blocks;
  *bsz_out = bsz;
  *nblk_out = nblk;

  gkyl_free(bval);
  gkyl_free(bsup);
}

// The tensor basis is symmetric under exchanging its two coordinates up to a
// permutation of its functions: b_k(x, y) = b_t[k](y, x). The permutation is
// found by evaluation rather than written out, so it follows the basis's own
// ordering; if none exists cubic_transpose_ok stays false.
static void
tok_geo_cubic_transpose(struct gkyl_tok_geo *geo)
{
  const struct gkyl_basis *basis = &geo->rzbasis_cubic;
  geo->cubic_transpose_ok = false;
  if (basis->num_basis != 16 || basis->ndim != 2 || basis->eval == 0) {
    return;
  }

  // Two points with no symmetry of their own, so that distinct functions
  // cannot agree at both.
  const double pt[2][2] = {
    {0.3141592653589793, -0.7071067811865476}, {-0.5772156649015329, 0.6180339887498949}
  };
  double v[2][16], w[2][16];
  for (int s = 0; s < 2; ++s) {
    double z[2] = {pt[s][0], pt[s][1]}, zt[2] = {pt[s][1], pt[s][0]};
    basis->eval(z, v[s]);
    basis->eval(zt, w[s]);
  }
  int t[16];
  for (int k = 0; k < 16; ++k) {
    int best = 0;
    double dbest = DBL_MAX;
    for (int j = 0; j < 16; ++j) {
      double d = fabs(v[0][k] - w[0][j]) + fabs(v[1][k] - w[1][j]);
      if (d < dbest) {
        dbest = d;
        best = j;
      }
    }
    // The same polynomial evaluated twice: equal to a few roundings.
    if (dbest > 64.0 * DBL_EPSILON * fmax(1.0, fabs(v[0][k]) + fabs(v[1][k]))) {
      return;
    }
    t[k] = best;
  }
  for (int k = 0; k < 16; ++k) {
    if (t[t[k]] != k) {
      return; // an exchange of coordinates is an involution
    }
  }
  for (int k = 0; k < 16; ++k) {
    geo->cubic_transpose[k] = t[k];
  }
  geo->cubic_transpose_ok = true;
}

static void
tok_geo_calc_psi_cell_bounds(struct gkyl_tok_geo *geo)
{
  // geo comes from gkyl_malloc: every field this sets is set here first,
  // whatever is built below.
  geo->psi_cell_bounds = 0;
  geo->psi_block_bounds = 0;
  geo->psi_block_size = geo->psi_num_blocks = 0;
  geo->psi_cell_bounds_cubic = 0;
  geo->psi_block_bounds_cubic = 0;
  geo->psi_block_size_cubic = geo->psi_num_blocks_cubic = 0;
  geo->cubic_transpose_ok = false;

  // Each enclosure only where its expansion is the one being solved: the
  // quadratic one without use_cubics (unchanged), the cubic one with it
  // (2026-10-03, for R_psiZ_cubic and the cubic Z solve).
  if (!geo->use_cubics) {
    if (geo->rzbasis.eval != 0 && geo->rzbasis.ndim == 2) {
      tok_psi_enclosure(
        &geo->rzbasis, geo->psiRZ, &geo->rzlocal, &geo->psi_cell_bounds, &geo->psi_block_bounds,
        &geo->psi_block_size, &geo->psi_num_blocks
      );
    }
  } else if (geo->rzbasis_cubic.eval != 0 && geo->rzbasis_cubic.ndim == 2) {
    tok_psi_enclosure(
      &geo->rzbasis_cubic, geo->psiRZ_cubic, &geo->rzlocal_cubic, &geo->psi_cell_bounds_cubic,
      &geo->psi_block_bounds_cubic, &geo->psi_block_size_cubic, &geo->psi_num_blocks_cubic
    );
    tok_geo_cubic_transpose(geo);
  }
}

struct gkyl_tok_geo *
gkyl_tok_geo_new(const struct gkyl_efit_inp *inp, const struct gkyl_tok_geo_grid_inp *ginp)
{
  struct gkyl_tok_geo *geo = gkyl_malloc(sizeof(*geo));
  *geo = (struct gkyl_tok_geo){};

  geo->efit = gkyl_efit_new(inp);

  geo->plate_spec = ginp->plate_spec;
  geo->extend_to_limiter = ginp->extend_to_limiter;
  geo->plate_func_lower = ginp->plate_func_lower;
  geo->plate_func_upper = ginp->plate_func_upper;
  for (int slot = 0; slot < 2; ++slot) {
    const struct gkyl_tok_geo_wall_target *target = &ginp->divertor_wall[slot];
    if (!target->num_segments && !target->segments) {
      continue;
    }
    if (!ginp->plate_spec || !(slot ? ginp->plate_func_upper : ginp->plate_func_lower) ||
        ginp->plate_func_lower == ginp->plate_func_upper ||
        !tok_wall_target_valid(geo->efit, target->num_segments, target->segments)) {
      fprintf(
        stderr, "TOK_GEO_DIVERTOR_TARGET_INVALID slot=%d count=%d\n", slot, target->num_segments
      );
      abort();
    }
    int *segments = gkyl_malloc(target->num_segments * sizeof(int));
    memcpy(segments, target->segments, target->num_segments * sizeof(int));
    geo->divertor_wall[slot] = (struct gkyl_tok_geo_wall_target){target->num_segments, segments};
    if (geo->extend_to_limiter) {
      fprintf(stderr, "TOK_GEO_DIVERTOR_TARGET slot=%d segments=", slot);
      for (int k = 0; k < target->num_segments; ++k) {
        fprintf(stderr, "%s%d", k ? "," : "", segments[k]);
      }
      fprintf(stderr, "\n");
    }
  }

  geo->rzbasis = geo->efit->rzbasis;
  geo->rzbasis_cubic = geo->efit->rzbasis_cubic;
  geo->rzgrid = geo->efit->rzgrid;
  geo->rzgrid_cubic = geo->efit->rzgrid_cubic;
  geo->psiRZ = gkyl_array_acquire(geo->efit->psizr);
  geo->psiRZ_cubic = gkyl_array_acquire(geo->efit->psizr_cubic);

  geo->num_rzbasis = geo->rzbasis.num_basis;
  geo->rzlocal = geo->efit->rzlocal;
  geo->rzlocal_ext = geo->efit->rzlocal_ext;
  geo->rzlocal_cubic = geo->efit->rzlocal_cubic;
  geo->rzlocal_cubic_ext = geo->efit->rzlocal_cubic_ext;
  geo->fgrid = geo->efit->fluxgrid;
  geo->fbasis = geo->efit->fluxbasis;
  geo->frange = geo->efit->fluxlocal;
  geo->frange_ext = geo->efit->fluxlocal_ext;
  geo->fpoldg = gkyl_array_acquire(geo->efit->fpolflux);
  geo->fpolprimedg = gkyl_array_acquire(geo->efit->fpolprimeflux);
  geo->qdg = gkyl_array_acquire(geo->efit->qflux);
  geo->sibry = geo->efit->sibry;
  // psisep MUST come from the same representation that evaluates psi.
  // `use_cubics` switches the evaluator (tok_eval_psi_rz_local) and the X-point
  // coordinates (Zxpt_cubic/Rxpt_cubic) but this line used to take the
  // quadratic separatrix value unconditionally.  The two differ -- measured
  // 2.8e-6 on NSTX-U 202778 -- so every separatrix trace point missed its own
  // level set by that amount and the trace was rejected
  // ("TOK_EXT_TRACE invalid point ... residual=2.05e-06"), making use_cubics
  // unusable.  Selecting the matching value is a consistency requirement, not
  // a tunable.
  geo->psisep = ginp->use_cubics ? geo->efit->psisep_cubic : geo->efit->psisep;
  geo->zmaxis = geo->efit->zmaxis;

  geo->use_cubics = ginp->use_cubics;
  geo->use_hyperbolic_numbers = ginp->use_hyperbolic_numbers;
  geo->root_param.eps = ginp->root_param.eps > 0 ? ginp->root_param.eps : 1e-10;
  geo->root_param.max_iter = ginp->root_param.max_iter > 0 ? ginp->root_param.max_iter : 100;

  geo->quad_param.max_level = ginp->quad_param.max_levels > 0 ? ginp->quad_param.max_levels : 10;
  geo->quad_param.eps = ginp->quad_param.eps > 0 ? ginp->quad_param.eps : 1e-10;

  if (geo->use_cubics) {
    if (geo->use_hyperbolic_numbers) {
      geo->calc_roots = calc_RdR_p3_hyperbolic;
    } else {
      geo->calc_roots = calc_RdR_p3;
    }
    geo->calc_grad_psi = calc_grad_psi_p3;
  } else if (geo->efit->rzbasis.poly_order == 1) {
    geo->calc_roots = calc_RdR_p1;
    geo->calc_grad_psi = calc_grad_psi_p1;
  } else if (geo->efit->rzbasis.poly_order == 2) {
    geo->calc_roots = calc_RdR_p2_tensor_nrc;
    geo->calc_grad_psi = calc_grad_psi_p2_tensor;
  }

  geo->stat = (struct gkyl_tok_geo_stat){};

  tok_geo_calc_psi_cell_bounds(geo);

  return geo;
}

double
gkyl_tok_geo_integrate_psi_contour(
  const struct gkyl_tok_geo *geo, double psi, double zmin, double zmax, double rclose
)
{
  return integrate_psi_contour_memo(geo, psi, zmin, zmax, rclose, false, false, 0);
}

int
gkyl_tok_geo_R_psiZ(
  const struct gkyl_tok_geo *geo, double psi, double Z, int nmaxroots, double *R, double *dRdZ,
  double *dR, double *dZ
)
{
  if (geo->use_cubics) {
    return R_psiZ_cubic(geo, psi, Z, nmaxroots, R, dRdZ, dR, dZ);
  } else {
    return R_psiZ(geo, psi, Z, nmaxroots, R, dRdZ, dR, dZ);
  }
}

// DIAGNOSTIC (2026-10-06): GKYL_TOK_THETA_ARC_DIAG=<prefix> writes, for every
// psi row of this block's corner nodes, the exact contour arc of each theta
// cell (tok_contour_leaves between consecutive corner nodes, on the row's own
// psi), to <prefix>-ftype<NN>.txt. The seam test measured a cell by a cubic
// through four of its points, which cannot follow a row bending on the scale of
// a cell -- beside the X point -- and re-measuring one grid with finer nested
// nodes showed that error to dominate what it reported (STEP theta x4: 4.3e-4,
// then 9.1e-5 and 1.05e-5 with x2 and x4 finer nodes). Off by default; the
// grid is unchanged either way. A cell whose arc cannot be traced is written
// with ok=0. The last pass for a block wins, so a reader matches the nodes.
static void
tok_theta_arc_diag(
  const struct gk_geometry *up, const struct gkyl_range *nrange, const struct gkyl_tok_geo *geo,
  const struct gkyl_tok_geo_grid_inp *inp
)
{
  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates
  const char *prefix = getenv("GKYL_TOK_THETA_ARC_DIAG");
  if (!prefix || prefix[0] == '\0') {
    return;
  }
  char path[1024];
  snprintf(path, sizeof path, "%s-ftype%02d.txt", prefix, (int)inp->ftype);
  FILE *fp = fopen(path, "w");
  if (!fp) {
    fprintf(stderr, "TOK_THETA_ARC_DIAG cannot_open path=%s\n", path);
    return;
  }
  fprintf(
    fp, "# ftype=%d rows=%d cells=%d: ip it psi r0 z0 r1 z1 arc ok\n", (int)inp->ftype,
    nrange->upper[PSI_IDX] - nrange->lower[PSI_IDX] + 1,
    nrange->upper[TH_IDX] - nrange->lower[TH_IDX]
  );
  struct tok_leaves lv = {0};
  int cidx[3] = {0};
  cidx[AL_IDX] = nrange->lower[AL_IDX];
  for (int ip = nrange->lower[PSI_IDX]; ip <= nrange->upper[PSI_IDX]; ++ip) {
    cidx[PSI_IDX] = ip;
    for (int it = nrange->lower[TH_IDX]; it < nrange->upper[TH_IDX]; ++it) {
      cidx[TH_IDX] = it;
      const double *p0 = gkyl_array_cfetch(up->geo_corn.mc2p_nodal, gkyl_range_idx(nrange, cidx));
      const double *nu =
        gkyl_array_cfetch(up->geo_corn.mc2nu_pos_nodal, gkyl_range_idx(nrange, cidx));
      cidx[TH_IDX] = it + 1;
      const double *p1 = gkyl_array_cfetch(up->geo_corn.mc2p_nodal, gkyl_range_idx(nrange, cidx));
      const double psi = nu[X_IDX];
      lv.n = 0;
      const bool ok =
        tok_contour_leaves(geo, psi, p0[X_IDX], p0[Y_IDX], p1[X_IDX], p1[Y_IDX], 0, &lv);
      double arc = 0.0;
      for (int j = 0; j < lv.n; ++j) {
        arc += lv.arc[j];
      }
      fprintf(
        fp, "%d %d %.17g %.17g %.17g %.17g %.17g %.17g %d\n", ip - nrange->lower[PSI_IDX],
        it - nrange->lower[TH_IDX], psi, p0[X_IDX], p0[Y_IDX], p1[X_IDX], p1[Y_IDX], arc, (int)ok
      );
    }
  }
  if (lv.r) {
    gkyl_free(lv.r);
    gkyl_free(lv.z);
    gkyl_free(lv.arc);
  }
  fclose(fp);
}

// Whether this block's containment tests run: only when the input asks for
// them. Asked for against an outline that cannot bound a region, the block is
// rejected here (the multiblock preflight rejects it earlier, with the block
// number). Not asked for, the block says so, so an unenforced wall is never
// silent in the log.
static bool
tok_wall_enforced(const struct gkyl_tok_geo *geo, const struct gkyl_tok_geo_grid_inp *inp)
{
  enum gkyl_tok_wall_policy policy = gkyl_tok_wall_policy_for(inp, geo->efit);
  if (policy == GKYL_TOK_WALL_ENFORCE) {
    return true;
  }
  if (policy == GKYL_TOK_WALL_NOT_ENFORCED) {
    fprintf(
      stderr, "TOK_GEO_WALL_NOT_ENFORCED ftype=%d limiter_status=%d vertices=%d reason=%s\n",
      inp->ftype, geo->efit->limiter_status, geo->efit->limiter_n,
      gkyl_tok_wall_policy_reason(policy)
    );
    return false;
  }
  fprintf(
    stderr, "TOK_GEO_WALL_UNAVAILABLE ftype=%d limiter_status=%d vertices=%d reason=%s\n",
    inp->ftype, geo->efit->limiter_status, geo->efit->limiter_n, gkyl_tok_wall_policy_reason(policy)
  );
  abort();
}

void
gkyl_tok_geo_calc(
  struct gk_geometry *up, struct gkyl_range *nrange, struct gkyl_tok_geo *geo,
  struct gkyl_tok_geo_grid_inp *inp, struct gkyl_position_map *position_map
)
{
  // Vessel-outline policy for this block, decided once.
  const bool enforce_wall = tok_wall_enforced(geo, inp);

  geo->rleft = inp->rleft;
  geo->rright = inp->rright;

  geo->inexact_roots = inp->inexact_roots;

  geo->rmax = inp->rmax;
  geo->rmin = inp->rmin;

  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates

  double dtheta = inp->cgrid.dx[TH_IDX], dpsi = inp->cgrid.dx[PSI_IDX],
         dalpha = inp->cgrid.dx[AL_IDX];

  double theta_lo = up->grid.lower[TH_IDX] +
                    (up->local.lower[TH_IDX] - up->global.lower[TH_IDX]) * up->grid.dx[TH_IDX],
         psi_lo = up->grid.lower[PSI_IDX] +
                  (up->local.lower[PSI_IDX] - up->global.lower[PSI_IDX]) * up->grid.dx[PSI_IDX],
         alpha_lo = up->grid.lower[AL_IDX] +
                    (up->local.lower[AL_IDX] - up->global.lower[AL_IDX]) * up->grid.dx[AL_IDX];

  double dx_fact = 1.0 / up->basis.poly_order;
  dtheta *= dx_fact;
  dpsi *= dx_fact;
  dalpha *= dx_fact;

  double rclose = inp->rclose;
  double rright = inp->rright;
  double rleft = inp->rleft;

  int nzcells;
  if (geo->use_cubics) {
    nzcells = geo->rzgrid_cubic.cells[1];
  } else {
    nzcells = geo->rzgrid.cells[1];
  }
  double *arc_memo = gkyl_malloc(sizeof(double[nzcells]));
  double *arc_memo_left = gkyl_malloc(sizeof(double[nzcells]));
  double *arc_memo_right = gkyl_malloc(sizeof(double[nzcells]));
  // Trace buffers: see tok_sep_trace_capacity.
  int sep_trace_capacity = tok_sep_trace_capacity(inp, nzcells);
  double *sep_trace_r = gkyl_malloc(sizeof(double[sep_trace_capacity]));
  double *sep_trace_z = gkyl_malloc(sizeof(double[sep_trace_capacity]));
  double *sep_trace_s = gkyl_malloc(sizeof(double[sep_trace_capacity]));
  // The far-boundary, correspondence and map traces share one allocation.
  double *ordered_trace_storage = gkyl_malloc(sizeof(double) * 8 * (size_t)sep_trace_capacity);

  struct arc_length_ctx arc_ctx = {
    .geo = geo,
    .arc_memo = arc_memo,
    .arc_memo_right = arc_memo_right,
    .arc_memo_left = arc_memo_left,
    .sep_trace_r = sep_trace_r,
    .sep_trace_z = sep_trace_z,
    .sep_trace_s = sep_trace_s,
    .sep_trace_capacity = sep_trace_capacity,
    .far_trace_r = ordered_trace_storage,
    .far_trace_z = ordered_trace_storage + sep_trace_capacity,
    .far_trace_s = ordered_trace_storage + 2 * sep_trace_capacity,
    .trace_corr_v = ordered_trace_storage + 3 * sep_trace_capacity,
    .map_trace_r = ordered_trace_storage + 4 * sep_trace_capacity,
    .map_trace_z = ordered_trace_storage + 5 * sep_trace_capacity,
    .map_trace_s = ordered_trace_storage + 6 * sep_trace_capacity,
    .map_trace_phi = ordered_trace_storage + 7 * sep_trace_capacity,
    .ext_ladder_w = NULL,
    .ext_ladder_rf = NULL,
    .ftype = inp->ftype,
    .zmaxis = geo->zmaxis,
  };
  struct plate_ctx pctx = {.geo = geo};
  tok_init_xpt_ray_target(inp, position_map, &arc_ctx);

  int cidx[3] = {0};
  for (int ia = nrange->lower[AL_IDX]; ia < nrange->lower[AL_IDX] + 1; ++ia) {
    cidx[AL_IDX] = ia;
    double alpha_curr = alpha_lo + ia * dalpha;

    for (int ip = nrange->lower[PSI_IDX]; ip <= nrange->upper[PSI_IDX]; ++ip) {
      double psi_curr = psi_lo + ip * dpsi;

      // Non-uniform psi. Finite differences are calculated in calc_metric.c
      position_map->maps[0](0.0, &psi_curr, &psi_curr, position_map->ctxs[0]);

      double darcL, arcL_curr, arcL_lo;

      // For double null blocks this should set arc_ctx :
      // zmin, zmax, rclose, arcL_tot for all blocks. No left and right
      // For a full core case:
      // also set phi_right and arcL_right
      // For a single null case:
      // also set zmin_left and zmin_right
      if (tok_xpt_ordered_placement(inp)) {
        tok_prepare_ordered_map(inp, &arc_ctx, psi_curr);
      } else {
        tok_find_endpoints(
          inp, geo, &arc_ctx, &pctx, psi_curr, alpha_curr, arc_memo, arc_memo_left, arc_memo_right
        );
      }

      darcL = (arc_ctx.arc_hi - arc_ctx.arc_lo) / (up->basis.poly_order * inp->cgrid.cells[TH_IDX]);
      // At the beginning of each theta loop we need to reset things.
      cidx[PSI_IDX] = ip;
      arcL_curr = 0.0;
      arcL_lo = tok_arc_from_theta(inp, &arc_ctx, theta_lo);
      double ridders_min, ridders_max;
      // Set node coordinates.
      for (int it = nrange->lower[TH_IDX]; it <= nrange->upper[TH_IDX]; ++it) {
        int it_delta = 0;
        arcL_curr = arcL_lo + it * darcL;
        double theta_curr = tok_theta_from_arc(inp, &arc_ctx, arcL_curr);

        // Calculate derivatives using finite difference for ddtheta,
        // as well as transform the computational coordiante to the non-uniform field-aligned value

        // We cannot do non-uniform alpha because we are modeling axisymmetric systems
        // Non-uniform theta
        double Theta_curr;
        position_map->maps[2](0.0, &theta_curr, &Theta_curr, position_map->ctxs[2]);
        theta_curr = Theta_curr;

        struct tok_ordered_point ordered = {0.0};
        bool ordered_mapping =
          tok_ordered_map_lookup(inp, &arc_ctx, theta_curr, alpha_curr, &ordered);
        if (tok_xpt_ordered_placement(inp) && !ordered_mapping) {
          fprintf(
            stderr, "TOK_ORDERED_MAP lookup failed ftype=%d psi=%.17g theta=%.17g\n", inp->ftype,
            psi_curr, theta_curr
          );
          abort();
        }
        double march_arc_req = arcL_curr;
        double r_curr = 0.0, z_curr = 0.0, phi_curr = 0.0;
        if (ordered_mapping) {
          r_curr = ordered.r;
          z_curr = ordered.z;
          phi_curr = ordered.phi;
        } else {
          arcL_curr = tok_xpt_theta_to_arc(inp, &arc_ctx, theta_curr);

          tok_set_ridders(inp, &arc_ctx, psi_curr, arcL_curr, &rclose, &ridders_min, &ridders_max);

          // single-null SOL rows: the node at exact arc (tok_lsn_exact_row)
          double lsn_r = 0.0, lsn_z = 0.0;
          const bool lsn_exact =
            arc_ctx.arcL_tot > 0.0 &&
            tok_lsn_exact_point(geo, psi_curr, arcL_curr / arc_ctx.arcL_tot, &lsn_r, &lsn_z);
          struct gkyl_qr_res res = lsn_exact ?
                                     (struct gkyl_qr_res){.res = lsn_z} :
                                     gkyl_ridders(
                                       arc_length_func, &arc_ctx, arc_ctx.zmin, arc_ctx.zmax,
                                       ridders_min, ridders_max, geo->root_param.max_iter, 1e-10
                                     );
          if (!lsn_exact) {
            tok_geo_check_arc_root(
              inp, psi_curr, theta_curr, arcL_curr, arc_ctx.zmin, arc_ctx.zmax, ridders_min,
              ridders_max, &res
            );
          }
          z_curr = res.res;
          ((struct gkyl_tok_geo *)geo)->stat.nroot_cont_calls += res.nevals;

          if (tok_geo_same_flux(psi_curr, geo->psisep)) {
            // The two ladders this replaces were separate `if`s, NOT else-if, so
            // a one-cell-wide block matched both and the LOWER clause overwrote
            // the upper. Evaluate lower last to preserve that exactly.
            const bool at_upper = it == nrange->upper[TH_IDX] &&
                                  (up->local.upper[TH_IDX] == up->global.upper[TH_IDX]);
            const bool at_lower = it == nrange->lower[TH_IDX] &&
                                  (up->local.lower[TH_IDX] == up->global.lower[TH_IDX]);
            double z_pin = 0.0;
            if (at_upper && tok_xpt_sep_pin_z(inp, geo, true, &z_pin)) {
              z_curr = z_pin;
            }
            if (at_lower && tok_xpt_sep_pin_z(inp, geo, false, &z_pin)) {
              z_curr = z_pin;
            }
          }

          double sep_r_curr = 0.0;
          bool at_sep_trace =
            tok_half_domain_sep_rz(inp, &arc_ctx, theta_curr, &sep_r_curr, &z_curr);

          double fixed_root_ref = rclose;
          bool at_fixed_edge = tok_xpt_at_fixed_edge(
            inp, &arc_ctx, it, nrange, up->local.lower[TH_IDX] == up->global.lower[TH_IDX],
            up->local.upper[TH_IDX] == up->global.upper[TH_IDX], &z_curr, &fixed_root_ref
          );
          bool at_xpt_anchor = arc_ctx.xpt_anchor_valid &&
                               tok_xpt_at_seam(
                                 inp, it, nrange,
                                 up->local.lower[TH_IDX] == up->global.lower[TH_IDX],
                                 up->local.upper[TH_IDX] == up->global.upper[TH_IDX]
                               );
          if (at_xpt_anchor) {
            z_curr = arc_ctx.xpt_anchor_z;
          }

          double R[4] = {0}, dRdZ[4] = {0};
          double dR[4] = {0}, dZ[4] = {0};
          int nr = gkyl_tok_geo_R_psiZ(geo, psi_curr, z_curr, 4, R, dRdZ, dR, dZ);
          double root_ref =
            at_sep_trace ?
              sep_r_curr :
              (at_xpt_anchor ? arc_ctx.xpt_anchor_r : (at_fixed_edge ? fixed_root_ref : rclose));
          r_curr = choose_closest(root_ref, R, R, nr);
          double drdz_curr = choose_closest(root_ref, R, dRdZ, nr);
          double dr_curr = choose_closest(root_ref, R, dR, nr);
          double dz_curr = choose_closest(root_ref, R, dZ, nr);
          if (at_xpt_anchor) {
            r_curr = arc_ctx.xpt_anchor_r;
          }
          if (at_sep_trace) {
            r_curr = sep_r_curr;
          }
          if (lsn_exact && z_curr == lsn_z && !at_xpt_anchor && !at_sep_trace) {
            r_curr = lsn_r; // the exact point's own R, on its branch
            if (nr <= 0) {
              nr = 1;
            }
          }

          if (tok_geo_same_flux(psi_curr, geo->psisep)) {
            // Snap to the X point of the representation that evaluates psi.
            const double *rx, *zx;
            const int nx = tok_geo_xpts(geo, &rx, &zx);
            for (int k = 0; k < nx && k < 2; ++k) {
              if (z_curr == zx[k]) {
                nr = 1;
                r_curr = rx[k];
              }
            }
          }

          if (nr == 0) {
            printf(" ip = %d, it = %d, ia = %d\n", ip, it, ia);
            printf(
              "Block Type = %d | Failed to find a root at psi = %g, Z = %1.16f\n", inp->ftype,
              psi_curr, z_curr
            );
            assert(false);
          }

          phi_curr = phi_func(alpha_curr, z_curr, &arc_ctx);
          if (lsn_exact && r_curr == lsn_r && z_curr == lsn_z) {
            // single-null SOL rows: the angle on the exact row (tok_lsn_exact_phi)
            double lsn_phi = 0.0;
            if (tok_lsn_exact_phi(
                  geo, inp->ftype, psi_curr, arcL_curr / arc_ctx.arcL_tot, r_curr, z_curr, &lsn_phi
                )) {
              phi_curr = alpha_curr + lsn_phi;
            }
          }
        }
        cidx[TH_IDX] = it;
        double *mc2p_n = gkyl_array_fetch(up->geo_corn.mc2p_nodal, gkyl_range_idx(nrange, cidx));
        double *mc2nu_n =
          gkyl_array_fetch(up->geo_corn.mc2nu_pos_nodal, gkyl_range_idx(nrange, cidx));
        double *bmag_n = gkyl_array_fetch(up->geo_corn.bmag_nodal, gkyl_range_idx(nrange, cidx));

        mc2p_n[X_IDX] = r_curr;
        mc2p_n[Y_IDX] = z_curr;
        mc2p_n[Z_IDX] = phi_curr;
        mc2nu_n[X_IDX] = psi_curr;
        mc2nu_n[Y_IDX] = alpha_curr;
        mc2nu_n[Z_IDX] = theta_curr;
        bmag_n[0] = bmag_func(r_curr, z_curr, &arc_ctx);
      }
    }
  }
  tok_theta_arc_diag(up, nrange, geo, inp);

  // Populate other alpha indices by using axisymmetry
  for (int ia = nrange->lower[AL_IDX] + 1; ia <= nrange->upper[AL_IDX]; ++ia) {
    cidx[AL_IDX] = ia;
    double alpha_curr = alpha_lo + ia * dalpha;
    double alpha_donor = alpha_lo + nrange->lower[AL_IDX] * dalpha;
    double alpha_diff = alpha_curr - alpha_donor;
    for (int ip = nrange->lower[PSI_IDX]; ip <= nrange->upper[PSI_IDX]; ++ip) {
      cidx[PSI_IDX] = ip;
      for (int it = nrange->lower[TH_IDX]; it <= nrange->upper[TH_IDX]; ++it) {
        cidx[TH_IDX] = it;

        double *mc2p_n = gkyl_array_fetch(up->geo_corn.mc2p_nodal, gkyl_range_idx(nrange, cidx));
        double *mc2nu_n =
          gkyl_array_fetch(up->geo_corn.mc2nu_pos_nodal, gkyl_range_idx(nrange, cidx));
        double *bmag_n = gkyl_array_fetch(up->geo_corn.bmag_nodal, gkyl_range_idx(nrange, cidx));

        int donor_cidx[3];
        donor_cidx[AL_IDX] = nrange->lower[AL_IDX];
        donor_cidx[PSI_IDX] = ip;
        donor_cidx[TH_IDX] = it;

        double *donor_mc2p_n =
          gkyl_array_fetch(up->geo_corn.mc2p_nodal, gkyl_range_idx(nrange, donor_cidx));
        double *donor_mc2nu_n =
          gkyl_array_fetch(up->geo_corn.mc2nu_pos_nodal, gkyl_range_idx(nrange, donor_cidx));
        double *donor_bmag_n =
          gkyl_array_fetch(up->geo_corn.bmag_nodal, gkyl_range_idx(nrange, donor_cidx));

        mc2p_n[X_IDX] = donor_mc2p_n[X_IDX];
        mc2p_n[Y_IDX] = donor_mc2p_n[Y_IDX];
        mc2p_n[Z_IDX] = donor_mc2p_n[Z_IDX] + alpha_diff;
        mc2nu_n[X_IDX] = donor_mc2nu_n[X_IDX];
        mc2nu_n[Y_IDX] = donor_mc2nu_n[AL_IDX] + alpha_diff;
        mc2nu_n[Z_IDX] = donor_mc2nu_n[Z_IDX];
        bmag_n[0] = donor_bmag_n[0];
      }
    }
  }

  // The material wall is a hard constraint, independent of extension and of
  // diagnostic fold overrides. A physical theta-face with an explicit wall
  // target uses its native material arc; all other edges retain chord checks.
  // Capture only diagnostic requested-boundary curves from the native corner
  // array. Reduced-dimensional app arrays have not retained these nodal data.
  if (wall_trial_active && wall_trial_capture) {
    int ip = wall_trial_movable_edge ? nrange->upper[PSI_IDX] : nrange->lower[PSI_IDX];
    for (int it = nrange->lower[TH_IDX]; it <= nrange->upper[TH_IDX]; ++it) {
      int idx[3] = {ip, nrange->lower[AL_IDX], it};
      const double *p = gkyl_array_cfetch(up->geo_corn.mc2p_nodal, gkyl_range_idx(nrange, idx));
      fprintf(
        stderr, "TOK_RHO_WALL_REQUESTED_CONTOUR block=%d node=%d rho=%.17g R=%.17g Z=%.17g\n",
        wall_trial_block, it - nrange->lower[TH_IDX], wall_trial_rho, p[0], p[1]
      );
    }
  }
  // Which node rows sit on a separately declared plate: the block's first or
  // last theta row, when that end is the block's physical boundary and the
  // topology says it is such a plate.
  const bool plate_row_lo = enforce_wall && up->local.lower[TH_IDX] == up->global.lower[TH_IDX] &&
                            tok_wall_theta_end_on_declared_plate(inp, geo, 0);
  const bool plate_row_up = enforce_wall && up->local.upper[TH_IDX] == up->global.upper[TH_IDX] &&
                            tok_wall_theta_end_on_declared_plate(inp, geo, 1);
  // Every test below, and the boundary-curve pass after it, judges against the
  // outline together with the pockets between these plates and the outline;
  // the plate-row nodes are pocket vertices.
  struct tok_wall_pocket pocket[2] = {{0}};
  if (enforce_wall) {
    tok_wall_block_pockets(
      inp, geo, (const bool[2]){plate_row_lo, plate_row_up}, up->geo_corn.mc2p_nodal, nrange, pocket
    );
  }
  tok_wall_pockets_set(&pocket[0], &pocket[1]);
  if (enforce_wall) {
    for (int ip = nrange->lower[PSI_IDX]; ip <= nrange->upper[PSI_IDX]; ++ip) {
      for (int it = nrange->lower[TH_IDX]; it <= nrange->upper[TH_IDX]; ++it) {
        int idx[3] = {ip, nrange->lower[AL_IDX], it};
        const double *p = gkyl_array_cfetch(up->geo_corn.mc2p_nodal, gkyl_range_idx(nrange, idx));
        const bool on_plate = (plate_row_lo && it == nrange->lower[TH_IDX]) ||
                              (plate_row_up && it == nrange->upper[TH_IDX]);
        bool ok = tok_wall_point_inside(geo->efit, p);
        // Report WHICH test failed. The old message printed only the node, so a
        // segment violation looked like a point violation -- and the node it
        // named was often comfortably inside the wall, which is actively
        // misleading when diagnosing a rejection.
        const char *fail_scope = ok ? "none" : "corner_node";
        double fail_q[2] = {p[0], p[1]};
        bool fixed_row =
          ip == (wall_trial_movable_edge == 0 ? nrange->upper[PSI_IDX] : nrange->lower[PSI_IDX]);
        bool fixed_failure = !ok && fixed_row;
        for (int d = 0; d < 3 && ok; d += 2) {
          if (idx[d] >= nrange->upper[d]) {
            continue;
          }
          // A psi segment along a plate row joins two nodes the root finder put
          // ON the plate, and the face it stands for follows the plate between
          // them; the plate was checked against the outline with its pocket.
          if (d == PSI_IDX && on_plate) {
            continue;
          }
          idx[d]++;
          const double *q = gkyl_array_cfetch(up->geo_corn.mc2p_nodal, gkyl_range_idx(nrange, idx));
          ok = tok_wall_segment_inside(geo->efit, p, q);
          if (!ok) {
            fail_scope = d == PSI_IDX ? "segment_psi" : "segment_theta";
            fail_q[0] = q[0];
            fail_q[1] = q[1];
          }
          if (!ok && d == PSI_IDX) {
            bool lower = it == nrange->lower[TH_IDX] &&
                         up->local.lower[TH_IDX] == up->global.lower[TH_IDX];
            bool upper = it == nrange->upper[TH_IDX] &&
                         up->local.upper[TH_IDX] == up->global.upper[TH_IDX];
            if (lower || upper) {
              ok = tok_divertor_material_cap(inp, geo, upper, p, q);
            }
            if (ok) {
              fail_scope = "none";
            }
          }
          if (!ok && d == TH_IDX && fixed_row) {
            fixed_failure = true;
          }
          idx[d]--;
        }
        if (!ok) {
          // Distinguish "a node is outside the machine" from "both nodes are
          // inside and the edge between them bulges out". The first is a
          // declaration to move; the second is a resolution to refine. They used
          // to print identically, so telling them apart meant running a
          // refinement sweep by hand.
          bool far_inside = tok_wall_point_inside(geo->efit, fail_q);
          bool segment_bulge = far_inside && strncmp(fail_scope, "segment_", 8) == 0;
          double excursion = segment_bulge ? tok_wall_segment_excursion(geo->efit, p, fail_q) : 0.0;
          fprintf(
            stderr,
            "TOK_GEO_WALL_DOMAIN_FAILED ftype=%d ip=%d it=%d rz=(%.17g,%.17g) "
            "scope=%s to_rz=(%.17g,%.17g) limiter_status=%d endpoints=%s excursion_m=%.17g\n",
            inp->ftype, ip, it, p[0], p[1], fail_scope, fail_q[0], fail_q[1],
            geo->efit->limiter_status, far_inside ? "inside" : "node_outside", excursion
          );
          // strcmp, not the boolean above: `fail_scope` is "corner_node" only
          // when the node itself failed tok_wall_point_inside, as opposed to an
          // edge between two inside nodes bulging out.
          // On the movable side when every node row the failed test involved is
          // no farther from the movable edge's row than from the fixed one: the
          // node itself, both ends of a psi segment, or the row of a theta one.
          bool on_side = false;
          if (wall_trial_movable_edge >= 0) {
            const int lo = nrange->lower[PSI_IDX], span = nrange->upper[PSI_IDX] - lo;
#define TOK_ROW_ON_MOVABLE_SIDE(r) \
  (wall_trial_movable_edge ? 2 * ((r) - lo) >= span : 2 * ((r) - lo) <= span)
            on_side = TOK_ROW_ON_MOVABLE_SIDE(ip) &&
                      (strcmp(fail_scope, "segment_psi") != 0 || TOK_ROW_ON_MOVABLE_SIDE(ip + 1));
#undef TOK_ROW_ON_MOVABLE_SIDE
          }
          if (!tok_wall_trial_record_where(
                fixed_failure, strcmp(fail_scope, "corner_node") == 0, on_side
              )) {
            abort();
          }
        }
      }
    }
  }
  tok_wall_pockets_set(0, 0);

  // A folded cell is a grid that is wrong, not merely poor: the map from
  // computational to physical coordinates has reversed orientation there, so
  // the Jacobian changes sign inside the block. Such a block used to be written
  // out and "pass", and was only caught post-hoc by check_grid_quality.py.
  // Check it here instead, on the corner nodes, and fail the build.
  //
  // Signed area of each (psi,theta) quad in the R-Z plane, at every alpha. The
  // sign convention is per block, so compare against the block's own dominant
  // sign rather than assuming positive.
  {
    int nfold = 0, nquad = 0;
    double worst = 0.0, worst_r = 0.0, worst_z = 0.0;
    int worst_ip = -1, worst_it = -1;
    long npos = 0, nneg = 0;
    int cidx0[3] = {0, 0, 0};
    for (int pass = 0; pass < 2; ++pass) {
      for (int ia = nrange->lower[AL_IDX]; ia <= nrange->upper[AL_IDX]; ++ia) {
        for (int ip = nrange->lower[PSI_IDX]; ip < nrange->upper[PSI_IDX]; ++ip) {
          for (int it = nrange->lower[TH_IDX]; it < nrange->upper[TH_IDX]; ++it) {
            double rq[4], zq[4];
            const int dp[4] = {0, 1, 1, 0}, dt[4] = {0, 0, 1, 1};
            for (int k = 0; k < 4; ++k) {
              cidx0[PSI_IDX] = ip + dp[k];
              cidx0[AL_IDX] = ia;
              cidx0[TH_IDX] = it + dt[k];
              const double *p =
                gkyl_array_cfetch(up->geo_corn.mc2p_nodal, gkyl_range_idx(nrange, cidx0));
              rq[k] = p[X_IDX];
              zq[k] = p[Y_IDX];
            }
            double a = 0.0;
            for (int k = 0; k < 4; ++k) {
              int k1 = (k + 1) % 4;
              a += rq[k] * zq[k1] - rq[k1] * zq[k];
            }
            a *= 0.5;
            if (pass == 0) {
              if (a > 0.0) {
                ++npos;
              } else if (a < 0.0) {
                ++nneg;
              }
              continue;
            }
            ++nquad;
            double signed_a = npos >= nneg ? a : -a;
            if (signed_a <= 0.0) {
              ++nfold;
              if (signed_a < worst) {
                worst = signed_a;
                worst_ip = ip;
                worst_it = it;
                worst_r = rq[0];
                worst_z = zq[0];
              }
            }
          }
        }
      }
    }
    if (nfold > 0) {
      fprintf(
        stderr,
        "TOK_GEO_FOLDED_CELLS ftype=%d nfold=%d of %d quads "
        "worst_area=%.17g at (ip=%d,it=%d) node=(%.17g,%.17g)\n",
        inp->ftype, nfold, nquad, worst, worst_ip, worst_it, worst_r, worst_z
      );
      // DIAGNOSTIC, not the verdict. This area is a shoelace over the four
      // CORNER nodes, but a cell edge is a curve: corner -> 2 Gauss surface
      // nodes -> corner. At the aspect ratios psi refinement produces at the
      // separatrix row the chord quadrilateral inverts while the cell does not
      // -- measured on five cells (STEP psi x8/x16, ASDEX psi x2/x4/x8, all
      // theta x1) whose Jacobian is sign-definite. Orientation is decided by the
      // signed-Jacobian guard in calc_metric.c, on by default, which evaluates J
      // at the quadrature points.
    }
  }

  // Consecutive flux surfaces crossing is the other way a block can be wrong
  // while still "passing": every quad keeps its orientation, so the fold check
  // above sees nothing, but the radial ordering reverses and the grid is
  // physically meaningless there. 205004's DN_SOL_OUT_MID emitted 44 reversed
  // radial steps (min cos -0.995) and was caught only post-hoc by
  // check_grid_quality.py.
  //
  // Same criterion as that script: walk the radial node line at fixed theta and
  // compare successive step directions. A pair counts as a reversal only well
  // past orthogonal, so ordinary curvature near an X point is not reported.
  // Like the fold count above this sums over both alpha nodes, so it is 2x what
  // check_grid_quality.py reports from a single alpha slice.
  {
    const double reversal_cos = -0.5;
    int nrev = 0, npair = 0;
    double worst_cos = 1.0, worst_r = 0.0, worst_z = 0.0;
    int worst_ip = -1, worst_it = -1;
    int cidx0[3] = {0, 0, 0};
    for (int ia = nrange->lower[AL_IDX]; ia <= nrange->upper[AL_IDX]; ++ia) {
      for (int it = nrange->lower[TH_IDX]; it <= nrange->upper[TH_IDX]; ++it) {
        for (int ip = nrange->lower[PSI_IDX]; ip + 2 <= nrange->upper[PSI_IDX]; ++ip) {
          double rp[3], zp[3];
          for (int k = 0; k < 3; ++k) {
            cidx0[PSI_IDX] = ip + k;
            cidx0[AL_IDX] = ia;
            cidx0[TH_IDX] = it;
            const double *p =
              gkyl_array_cfetch(up->geo_corn.mc2p_nodal, gkyl_range_idx(nrange, cidx0));
            rp[k] = p[X_IDX];
            zp[k] = p[Y_IDX];
          }
          double ar = rp[1] - rp[0], az = zp[1] - zp[0];
          double br = rp[2] - rp[1], bz = zp[2] - zp[1];
          double la = hypot(ar, az), lb = hypot(br, bz);
          if (!(la > 0.0) || !(lb > 0.0)) {
            continue;
          }
          ++npair;
          double c = (ar * br + az * bz) / (la * lb);
          if (c < reversal_cos) {
            ++nrev;
            if (c < worst_cos) {
              worst_cos = c;
              worst_ip = ip;
              worst_it = it;
              worst_r = rp[1];
              worst_z = zp[1];
            }
          }
        }
      }
    }
    if (nrev > 0) {
      fprintf(
        stderr,
        "TOK_GEO_SURFACE_CROSS ftype=%d nrev=%d of %d radial step pairs "
        "worst_cos=%.17g at (ip=%d,it=%d) node=(%.17g,%.17g)\n",
        inp->ftype, nrev, npair, worst_cos, worst_ip, worst_it, worst_r, worst_z
      );
      // REPORTED, NOT GATING since 2026-09-21, for the same reason the fold
      // check above was demoted: the criterion is a proxy that cannot answer
      // the question it is named for.
      //
      // `reversal_cos = -0.5` is a bare threshold with no reference to the
      // scale of anything it judges, chosen in 3e2aeb5b3 so that it "fires on
      // nothing beyond the two known" cases of the NSTX-U 450 -- tuned to a
      // dataset. An angle between successive radial steps cannot separate the
      // two things that produce it:
      //
      //   a CROSSING -- adjacent flux surfaces actually intersect, and the
      //      grid is invalid there;
      //   TANGENTIAL SLIP -- the nodes slide along theta between surfaces that
      //      remain strictly ordered, which is harmless.
      //
      // Measured 2026-09-21 on every case this gate was rejecting -- STEP's
      // outboard plate at 2.312 deg, and all five TCV core cells -- the
      // surfaces are ORDERED: offsetting each row along its neighbour's local
      // normal keeps a strictly positive sign at every node (min +1.06e-05 on
      // tcv_core090), and the signed-Jacobian guard is silent on the same
      // blocks. They were rejected for curving, not for crossing.
      //
      // Why dropping it does not open a hole: a crossing between ADJACENT rows
      // inverts the quad between them, and that is decided by the
      // signed-Jacobian guard in calc_metric.c at the QUADRATURE points. That
      // guard only became the fold check on 2026-09-19; when this cosine was
      // written the fold check was the corner shoelace, which is itself wrong
      // by up to an order of magnitude. This test was covering for that, and
      // that job is done.
      //
      // The count is still computed and printed every build, so a regression
      // stays visible.
    }
  }

  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&inp->cbasis, &inp->cgrid, false);
  gkyl_nodal_ops_n2m(
    n2m, &inp->cbasis, &inp->cgrid, nrange, &up->local, 3, up->geo_corn.mc2p_nodal,
    up->geo_corn.mc2p, false
  );
  // Validate each complete polynomial edge on both radial boundaries of the
  // block. Corner/face samples alone can miss a p2 overshoot between nodes.
  if (enforce_wall && inp->cbasis.poly_order > 2) {
    fprintf(
      stderr, "TOK_GEO_WALL_DOMAIN_FAILED unsupported boundary order=%d\n", inp->cbasis.poly_order
    );
    abort();
  }
  tok_wall_pockets_set(&pocket[0], &pocket[1]);
  if (enforce_wall) {
    for (int side = 0; side < 2; ++side) {
      int idx[3] = {
        side ? up->local.upper[PSI_IDX] : up->local.lower[PSI_IDX], up->local.lower[AL_IDX],
        up->local.lower[TH_IDX]
      };
      for (; idx[TH_IDX] <= up->local.upper[TH_IDX]; ++idx[TH_IDX]) {
        const double *coeff = gkyl_array_cfetch(up->geo_corn.mc2p, gkyl_range_idx(&up->local, idx));
        double points[3][2];
        for (int k = 0; k < 3; ++k) {
          double eta[3] = {side ? 1.0 : -1.0, 0.0, k - 1.0};
          for (int d = 0; d < 2; ++d) {
            points[k][d] = inp->cbasis.eval_expand(eta, coeff + d * inp->cbasis.num_basis);
          }
        }
        if (!tok_wall_curve_inside(geo->efit, points[0], points[1], points[2])) {
          fprintf(
            stderr,
            "TOK_GEO_WALL_DOMAIN_FAILED ftype=%d scope=radial_boundary_curve side=%d theta_cell=%d\n",
            inp->ftype, side, idx[TH_IDX]
          );
          if (!tok_wall_trial_record_where(
                side != wall_trial_movable_edge, false, side == wall_trial_movable_edge
              )) {
            abort();
          }
        }
      }
    }
  }
  tok_wall_pockets_set(0, 0);
  tok_wall_pocket_release(&pocket[0]);
  tok_wall_pocket_release(&pocket[1]);

  gkyl_nodal_ops_n2m(
    n2m, &inp->cbasis, &inp->cgrid, nrange, &up->local, 3, up->geo_corn.mc2nu_pos_nodal,
    up->geo_corn.mc2nu_pos, false
  );
  gkyl_nodal_ops_n2m(
    n2m, &inp->cbasis, &inp->cgrid, nrange, &up->local, 1, up->geo_corn.bmag_nodal,
    up->geo_corn.bmag, false
  );
  gkyl_nodal_ops_release(n2m);

  // Need 1/B for LBO collisions, computed weakly.
  gkyl_dg_inv_op_range(&inp->cbasis, 0, up->geo_corn.bmag_inv, 0, up->geo_corn.bmag, &up->local);

  gkyl_free(arc_memo);
  gkyl_free(arc_memo_left);
  gkyl_free(arc_memo_right);
  gkyl_free(sep_trace_r);
  gkyl_free(sep_trace_z);
  gkyl_free(sep_trace_s);
  gkyl_free(arc_ctx.ext_ladder_w);
  gkyl_free(arc_ctx.ext_ladder_rf);
  gkyl_free(ordered_trace_storage);
}

void
gkyl_tok_geo_calc_interior(
  struct gk_geometry *up, struct gkyl_range *nrange, double dzc[3], struct gkyl_tok_geo *geo,
  struct gkyl_tok_geo_grid_inp *inp, struct gkyl_position_map *position_map
)
{
  // Vessel-outline policy for this block, decided once.
  const bool enforce_wall = tok_wall_enforced(geo, inp);

  geo->rleft = inp->rleft;
  geo->rright = inp->rright;

  geo->inexact_roots = inp->inexact_roots;

  geo->rmax = inp->rmax;
  geo->rmin = inp->rmin;

  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates

  double dtheta = inp->cgrid.dx[TH_IDX], dpsi = inp->cgrid.dx[PSI_IDX],
         dalpha = inp->cgrid.dx[AL_IDX];

  double theta_lo = up->grid.lower[TH_IDX] +
                    (up->local.lower[TH_IDX] - up->global.lower[TH_IDX]) * up->grid.dx[TH_IDX],
         psi_lo = up->grid.lower[PSI_IDX] +
                  (up->local.lower[PSI_IDX] - up->global.lower[PSI_IDX]) * up->grid.dx[PSI_IDX],
         alpha_lo = up->grid.lower[AL_IDX] +
                    (up->local.lower[AL_IDX] - up->global.lower[AL_IDX]) * up->grid.dx[AL_IDX];

  double dels[2] = {1.0 / sqrt(3), 1.0 - 1.0 / sqrt(3)};
  theta_lo = theta_lo + dels[1] * dtheta / 2.0;
  psi_lo = psi_lo + dels[1] * dpsi / 2.0;
  alpha_lo = alpha_lo + dels[1] * dalpha / 2.0;

  double dx_fact = 1.0 / up->basis.poly_order;
  dtheta *= dx_fact;
  dpsi *= dx_fact;
  dalpha *= dx_fact;

  // used for finite differences
  double delta_alpha = dalpha * 1e-2;
  double delta_psi = dpsi * 1e-2;
  double delta_theta = dtheta * 1e-2;
  dzc[0] = delta_psi;
  dzc[1] = delta_alpha;
  dzc[2] = delta_theta;
  int modifiers[5] = {0, -1, 1, -2, 2};

  double rclose = inp->rclose;
  double rright = inp->rright;
  double rleft = inp->rleft;

  int nzcells;
  if (geo->use_cubics) {
    nzcells = geo->rzgrid_cubic.cells[1];
  } else {
    nzcells = geo->rzgrid.cells[1];
  }
  double *arc_memo = gkyl_malloc(sizeof(double[nzcells]));
  double *arc_memo_left = gkyl_malloc(sizeof(double[nzcells]));
  double *arc_memo_right = gkyl_malloc(sizeof(double[nzcells]));
  // Trace buffers: see tok_sep_trace_capacity.
  int sep_trace_capacity = tok_sep_trace_capacity(inp, nzcells);
  double *sep_trace_r = gkyl_malloc(sizeof(double[sep_trace_capacity]));
  double *sep_trace_z = gkyl_malloc(sizeof(double[sep_trace_capacity]));
  double *sep_trace_s = gkyl_malloc(sizeof(double[sep_trace_capacity]));
  // The far-boundary, correspondence and map traces share one allocation.
  double *ordered_trace_storage = gkyl_malloc(sizeof(double) * 8 * (size_t)sep_trace_capacity);

  struct arc_length_ctx arc_ctx = {
    .geo = geo,
    .arc_memo = arc_memo,
    .arc_memo_right = arc_memo_right,
    .arc_memo_left = arc_memo_left,
    .sep_trace_r = sep_trace_r,
    .sep_trace_z = sep_trace_z,
    .sep_trace_s = sep_trace_s,
    .sep_trace_capacity = sep_trace_capacity,
    .far_trace_r = ordered_trace_storage,
    .far_trace_z = ordered_trace_storage + sep_trace_capacity,
    .far_trace_s = ordered_trace_storage + 2 * sep_trace_capacity,
    .trace_corr_v = ordered_trace_storage + 3 * sep_trace_capacity,
    .map_trace_r = ordered_trace_storage + 4 * sep_trace_capacity,
    .map_trace_z = ordered_trace_storage + 5 * sep_trace_capacity,
    .map_trace_s = ordered_trace_storage + 6 * sep_trace_capacity,
    .map_trace_phi = ordered_trace_storage + 7 * sep_trace_capacity,
    .ext_ladder_w = NULL,
    .ext_ladder_rf = NULL,
    .ftype = inp->ftype,
    .zmaxis = geo->zmaxis,
  };
  struct plate_ctx pctx = {.geo = geo};
  tok_init_xpt_ray_target(inp, position_map, &arc_ctx);

  // Temporary array to store nodal q profile.
  struct gkyl_array *qprofile_nodal =
    gkyl_array_new(GKYL_DOUBLE, up->geo_int.bmag_nodal->ncomp, up->geo_int.bmag_nodal->size);

  int cidx[3] = {0};
  for (int ia = nrange->lower[AL_IDX]; ia < nrange->lower[AL_IDX] + 1; ++ia) {
    cidx[AL_IDX] = ia;
    double alpha_curr = calc_running_coord(alpha_lo, ia - nrange->lower[AL_IDX], dalpha);

    for (int ip = nrange->lower[PSI_IDX]; ip <= nrange->upper[PSI_IDX]; ++ip) {
      int ip_delta_max = 3;
      for (int ip_delta = 0; ip_delta < ip_delta_max; ip_delta++) {
        double psi_curr = calc_running_coord(psi_lo, ip - nrange->lower[PSI_IDX], dpsi) +
                          modifiers[ip_delta] * delta_psi;

        // Non-uniform psi. Finite differences are calculated in calc_metric.c
        double Psi_curr;
        position_map->maps[0](0.0, &psi_curr, &Psi_curr, position_map->ctxs[0]);
        double dPsi_dpsi =
          gkyl_position_map_slope(position_map, 0, psi_curr, delta_psi, ip, nrange);
        psi_curr = Psi_curr;

        double darcL, arcL_curr, arcL_lo;

        // For double null blocks this should set arc_ctx :
        // zmin, zmax, rclose, arcL_tot for all blocks. No left and right
        // For a full core case:
        // also set phi_right and arcL_right
        // For a single null case:
        // also set zmin_left and zmin_right
        double qprofile = 0.0;
        if (tok_xpt_ordered_placement(inp)) {
          if (!inp->half_domain) {
            // The ordered map changes only the coordinate construction.  Use
            // the established contour setup, with the feature flags disabled
            // on a local copy, to preserve the pre-feature q profile exactly.
            struct gkyl_tok_geo_grid_inp legacy_inp = *inp;
            struct arc_length_ctx qctx = arc_ctx;
            legacy_inp.straight_xpt_ray = false;
            legacy_inp.straight_core_xpt_ray = false;
            tok_find_endpoints(
              &legacy_inp, geo, &qctx, &pctx, psi_curr, alpha_curr, arc_memo, arc_memo_left,
              arc_memo_right
            );
            qprofile = qprofile_func(&qctx);
          }
          tok_prepare_ordered_map(inp, &arc_ctx, psi_curr);
        } else {
          tok_find_endpoints(
            inp, geo, &arc_ctx, &pctx, psi_curr, alpha_curr, arc_memo, arc_memo_left, arc_memo_right
          );
          qprofile = qprofile_func(&arc_ctx);
        }

        // Calculate the q profile
        // qhat = - F(psi) * s(psi) / (R * grad(psi))
        // q = integral_0^2pi qhat dtheta
        //   = F(psi)*s(psi) * integral 1/(R*grad(psi)) dl
        //   = 1/s(psi) * integral (dphidtheta) ; dphidtheta = F(psi)/(R*grad(psi))
        // The legacy half-domain block types return zero from qprofile_func;
        // the full-domain branch above evaluates the unchanged legacy path.

        darcL =
          (arc_ctx.arc_hi - arc_ctx.arc_lo) / (up->basis.poly_order * inp->cgrid.cells[TH_IDX]);
        // at the beginning of each theta loop we need to reset things
        cidx[PSI_IDX] = ip;
        arcL_curr = 0.0;
        arcL_lo = tok_arc_from_theta(inp, &arc_ctx, theta_lo);
        double ridders_min, ridders_max;

        for (int it = nrange->lower[TH_IDX]; it <= nrange->upper[TH_IDX]; ++it) {
          arcL_curr = calc_running_coord(arcL_lo, it - nrange->lower[TH_IDX], darcL);
          double theta_curr = tok_theta_from_arc(inp, &arc_ctx, arcL_curr);

          // Calculate derivatives using finite difference for ddtheta,
          // as well as transform the computational coordiante to the non-uniform field-aligned value

          // We cannot do non-uniform alpha because we are modeling axisymmetric systems
          // Non-uniform theta
          double Theta_curr;
          position_map->maps[2](0.0, &theta_curr, &Theta_curr, position_map->ctxs[2]);
          double dTheta_dtheta =
            gkyl_position_map_slope(position_map, 2, theta_curr, delta_theta, it, nrange);
          theta_curr = Theta_curr;

          struct tok_ordered_point ordered = {0.0};
          bool ordered_mapping =
            tok_ordered_map_lookup(inp, &arc_ctx, theta_curr, alpha_curr, &ordered);
          if (tok_xpt_ordered_placement(inp) && !ordered_mapping) {
            fprintf(
              stderr, "TOK_ORDERED_MAP lookup failed ftype=%d psi=%.17g theta=%.17g\n", inp->ftype,
              psi_curr, theta_curr
            );
            abort();
          }
          double r_curr = 0.0, z_curr = 0.0, phi_curr = 0.0;
          double drdz_curr = 0.0, dr_curr = 0.0, dz_curr = 0.0;
          if (ordered_mapping) {
            r_curr = ordered.r;
            z_curr = ordered.z;
            phi_curr = ordered.phi;
          } else {
            arcL_curr = tok_xpt_theta_to_arc(inp, &arc_ctx, theta_curr);

            tok_set_ridders(inp, &arc_ctx, psi_curr, arcL_curr, &rclose, &ridders_min, &ridders_max);

            // single-null SOL rows: the node at exact arc (tok_lsn_exact_row)
            double lsn_r = 0.0, lsn_z = 0.0;
            const bool lsn_exact =
              arc_ctx.arcL_tot > 0.0 &&
              tok_lsn_exact_point(geo, psi_curr, arcL_curr / arc_ctx.arcL_tot, &lsn_r, &lsn_z);
            struct gkyl_qr_res res = lsn_exact ?
                                       (struct gkyl_qr_res){.res = lsn_z} :
                                       gkyl_ridders(
                                         arc_length_func, &arc_ctx, arc_ctx.zmin, arc_ctx.zmax,
                                         ridders_min, ridders_max, geo->root_param.max_iter, 1e-10
                                       );
            if (!lsn_exact) {
              tok_geo_check_arc_root(
                inp, psi_curr, theta_curr, arcL_curr, arc_ctx.zmin, arc_ctx.zmax, ridders_min,
                ridders_max, &res
              );
            }
            z_curr = res.res;
            ((struct gkyl_tok_geo *)geo)->stat.nroot_cont_calls += res.nevals;

            double sep_r_curr = 0.0;
            bool at_sep_trace =
              tok_half_domain_sep_rz(inp, &arc_ctx, theta_curr, &sep_r_curr, &z_curr);

            // These are volume-interior Gauss nodes, even when their nodal
            // indices are the first or last in the local range.  The affine
            // arc map already places them correctly; only a true block-edge
            // node may be snapped exactly to the shared ray.
            bool at_xpt_anchor = false;

            double R[4] = {0}, dRdZ[4] = {0};
            double dR[4] = {0}, dZ[4] = {0};
            int nr = gkyl_tok_geo_R_psiZ(geo, psi_curr, z_curr, 4, R, dRdZ, dR, dZ);
            double root_ref = at_sep_trace ? sep_r_curr :
                                             (at_xpt_anchor ? arc_ctx.xpt_anchor_r : rclose);
            r_curr = choose_closest(root_ref, R, R, nr);
            drdz_curr = choose_closest(root_ref, R, dRdZ, nr);
            dr_curr = choose_closest(root_ref, R, dR, nr);
            dz_curr = choose_closest(root_ref, R, dZ, nr);
            if (at_xpt_anchor) {
              r_curr = arc_ctx.xpt_anchor_r;
            }
            if (at_sep_trace) {
              r_curr = sep_r_curr;
            }
            if (lsn_exact && z_curr == lsn_z && !at_xpt_anchor && !at_sep_trace) {
              r_curr = lsn_r; // the exact point's own R, on its branch
              if (nr <= 0) {
                nr = 1;
              }
            }

            if (tok_geo_same_flux(psi_curr, geo->psisep) && ip_delta == 0) {
              // Snap to the X point of the representation that evaluates psi.
              const double *rx, *zx;
              const int nx = tok_geo_xpts(geo, &rx, &zx);
              for (int k = 0; k < nx && k < 2; ++k) {
                if (z_curr == zx[k]) {
                  nr = 1;
                  r_curr = rx[k];
                }
              }
            }

            if (nr == 0) {
              printf("ip = %d, it = %d, ia = %d, ip_delta = %d\n", ip, it, ia, ip_delta);
              printf(
                "Block Type = %d | Failed to find a root at psi = %g, Z = %1.16f\n", inp->ftype,
                psi_curr, z_curr
              );
              assert(false);
            }

            phi_curr = phi_func(alpha_curr, z_curr, &arc_ctx);
            if (lsn_exact && r_curr == lsn_r && z_curr == lsn_z) {
              // single-null SOL rows: the angle on the exact row (tok_lsn_exact_phi)
              double lsn_phi = 0.0;
              if (tok_lsn_exact_phi(
                    geo, inp->ftype, psi_curr, arcL_curr / arc_ctx.arcL_tot, r_curr, z_curr,
                    &lsn_phi
                  )) {
                phi_curr = alpha_curr + lsn_phi;
              }
            }
          }
          cidx[TH_IDX] = it;
          int lidx = 0;
          if (ip_delta != 0) {
            lidx = 3 + 3 * (ip_delta - 1);
          }
          double *mc2p_fd_n =
            gkyl_array_fetch(up->geo_int.mc2p_nodal_fd, gkyl_range_idx(nrange, cidx));
          double *ddtheta_n =
            gkyl_array_fetch(up->geo_int.ddtheta_nodal, gkyl_range_idx(nrange, cidx));
          double *ddpsi_n = gkyl_array_fetch(up->geo_int.ddpsi_nodal, gkyl_range_idx(nrange, cidx));
          double *mc2p_n = gkyl_array_fetch(up->geo_int.mc2p_nodal, gkyl_range_idx(nrange, cidx));
          double *bmag_n = gkyl_array_fetch(up->geo_int.bmag_nodal, gkyl_range_idx(nrange, cidx));
          double *curlbhat_n =
            gkyl_array_fetch(up->geo_int.curlbhat_nodal, gkyl_range_idx(nrange, cidx));
          double *qprofile_n = gkyl_array_fetch(qprofile_nodal, gkyl_range_idx(nrange, cidx));

          mc2p_fd_n[lidx + X_IDX] = r_curr;
          mc2p_fd_n[lidx + Y_IDX] = z_curr;
          mc2p_fd_n[lidx + Z_IDX] = phi_curr;

          if (ip_delta == 0) {
            if (ordered_mapping) {
              ddtheta_n[0] = ordered.dr_dtheta * dTheta_dtheta;
              ddtheta_n[1] = ordered.dz_dtheta * dTheta_dtheta;
              ddtheta_n[2] = ordered.dphi_dtheta * dTheta_dtheta;
            } else {
              double darc_dtheta = arc_ctx.arc_interval_valid ?
                                     arc_ctx.arc_darc_dtheta :
                                     (arc_ctx.xpt_map_valid ? arc_ctx.xpt_map_darc_dtheta :
                                                              arc_ctx.arcL_tot / (2.0 * M_PI));
              ddtheta_n[0] = sin(atan2(dr_curr, dz_curr)) * darc_dtheta * dTheta_dtheta;
              ddtheta_n[1] = cos(atan2(dr_curr, dz_curr)) * darc_dtheta * dTheta_dtheta;
              ddtheta_n[2] = dphidtheta_func(z_curr, &arc_ctx) * dTheta_dtheta;
            }
            ddpsi_n[0] = dPsi_dpsi;
            mc2p_n[lidx + X_IDX] = r_curr;
            mc2p_n[lidx + Y_IDX] = z_curr;
            mc2p_n[lidx + Z_IDX] = phi_curr;
            bmag_n[0] = bmag_func(r_curr, z_curr, &arc_ctx);
            qprofile_n[0] = qprofile;
            curlbhat_func(psi_curr, r_curr, z_curr, phi_curr, curlbhat_n, &arc_ctx);
          }
        }
      }
    }
  }

  // Populate other alpha indices by using axisymmetry
  for (int ia = nrange->lower[AL_IDX] + 1; ia <= nrange->upper[AL_IDX]; ++ia) {
    cidx[AL_IDX] = ia;
    double alpha_curr = calc_running_coord(alpha_lo, ia - nrange->lower[AL_IDX], dalpha);
    double alpha_donor = calc_running_coord(alpha_lo, 0, dalpha);
    double alpha_diff = alpha_curr - alpha_donor;
    for (int ip = nrange->lower[PSI_IDX]; ip <= nrange->upper[PSI_IDX]; ++ip) {
      cidx[PSI_IDX] = ip;
      int ip_delta_max = 3;
      for (int ip_delta = 0; ip_delta < ip_delta_max; ip_delta++) {
        double psi_curr = calc_running_coord(psi_lo, ip - nrange->lower[PSI_IDX], dpsi) +
                          modifiers[ip_delta] * delta_psi;
        for (int it = nrange->lower[TH_IDX]; it <= nrange->upper[TH_IDX]; ++it) {
          cidx[TH_IDX] = it;
          int lidx = 0;
          if (ip_delta != 0) {
            lidx = 3 + 3 * (ip_delta - 1);
          }

          double *mc2p_fd_n =
            gkyl_array_fetch(up->geo_int.mc2p_nodal_fd, gkyl_range_idx(nrange, cidx));
          double *ddtheta_n =
            gkyl_array_fetch(up->geo_int.ddtheta_nodal, gkyl_range_idx(nrange, cidx));
          double *ddpsi_n = gkyl_array_fetch(up->geo_int.ddpsi_nodal, gkyl_range_idx(nrange, cidx));
          double *mc2p_n = gkyl_array_fetch(up->geo_int.mc2p_nodal, gkyl_range_idx(nrange, cidx));
          double *bmag_n = gkyl_array_fetch(up->geo_int.bmag_nodal, gkyl_range_idx(nrange, cidx));
          double *curlbhat_n =
            gkyl_array_fetch(up->geo_int.curlbhat_nodal, gkyl_range_idx(nrange, cidx));

          int donor_cidx[3];
          donor_cidx[AL_IDX] = nrange->lower[AL_IDX];
          donor_cidx[PSI_IDX] = ip;
          donor_cidx[TH_IDX] = it;

          double *donor_mc2p_fd_n =
            gkyl_array_fetch(up->geo_int.mc2p_nodal_fd, gkyl_range_idx(nrange, donor_cidx));
          double *donor_ddtheta_n =
            gkyl_array_fetch(up->geo_int.ddtheta_nodal, gkyl_range_idx(nrange, donor_cidx));
          double *donor_ddpsi_n =
            gkyl_array_fetch(up->geo_int.ddpsi_nodal, gkyl_range_idx(nrange, donor_cidx));
          double *donor_mc2p_n =
            gkyl_array_fetch(up->geo_int.mc2p_nodal, gkyl_range_idx(nrange, donor_cidx));
          double *donor_bmag_n =
            gkyl_array_fetch(up->geo_int.bmag_nodal, gkyl_range_idx(nrange, donor_cidx));
          double *donor_curlbhat_n =
            gkyl_array_fetch(up->geo_int.curlbhat_nodal, gkyl_range_idx(nrange, donor_cidx));

          mc2p_fd_n[lidx + X_IDX] = donor_mc2p_fd_n[lidx + X_IDX];
          mc2p_fd_n[lidx + Y_IDX] = donor_mc2p_fd_n[lidx + Y_IDX];
          mc2p_fd_n[lidx + Z_IDX] = donor_mc2p_fd_n[lidx + Z_IDX] + alpha_diff;
          if (ip_delta == 0) {
            ddtheta_n[0] = donor_ddtheta_n[0];
            ddtheta_n[1] = donor_ddtheta_n[1];
            ddtheta_n[2] = donor_ddtheta_n[2];
            ddpsi_n[0] = donor_ddpsi_n[0];
            mc2p_n[lidx + X_IDX] = donor_mc2p_n[lidx + X_IDX];
            mc2p_n[lidx + Y_IDX] = donor_mc2p_n[lidx + Y_IDX];
            mc2p_n[lidx + Z_IDX] = donor_mc2p_n[lidx + Z_IDX] + alpha_diff;
            bmag_n[0] = donor_bmag_n[0];
            curlbhat_func(
              psi_curr, mc2p_n[lidx + X_IDX], mc2p_n[lidx + Y_IDX], mc2p_n[lidx + Z_IDX],
              curlbhat_n, &arc_ctx
            );
          }
        }
      }
    }
  }

  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&inp->cbasis, &inp->cgrid, false);
  // Judged against the outline together with the pockets of the block's
  // declared plates (see the corner pass); no interior node is on a plate.
  struct tok_wall_pocket pocket[2] = {{0}};
  if (enforce_wall) {
    tok_wall_block_pockets(
      inp, geo,
      (const bool[2]){
        up->local.lower[TH_IDX] == up->global.lower[TH_IDX] &&
          tok_wall_theta_end_on_declared_plate(inp, geo, 0),
        up->local.upper[TH_IDX] == up->global.upper[TH_IDX] &&
          tok_wall_theta_end_on_declared_plate(inp, geo, 1),
      },
      0, nrange, pocket
    );
  }
  tok_wall_pockets_set(&pocket[0], &pocket[1]);
  struct gkyl_range_iter wall_iter;
  gkyl_range_iter_init(&wall_iter, nrange);
  while (enforce_wall && gkyl_range_iter_next(&wall_iter)) {
    const double *p =
      gkyl_array_cfetch(up->geo_int.mc2p_nodal, gkyl_range_idx(nrange, wall_iter.idx));
    if (!tok_wall_point_inside(geo->efit, p)) {
      fprintf(
        stderr, "TOK_GEO_WALL_DOMAIN_FAILED ftype=%d scope=interior rz=(%.17g,%.17g)\n", inp->ftype,
        p[0], p[1]
      );
      // the node's radial coordinate, as the construction loop placed it
      double x = calc_running_coord(psi_lo, wall_iter.idx[PSI_IDX] - nrange->lower[PSI_IDX], dpsi);
      if (!tok_wall_trial_record_where(
            false, false,
            tok_wall_trial_on_movable_side(x, up->grid.lower[PSI_IDX], up->grid.upper[PSI_IDX])
          )) {
        abort();
      }
    }
  }
  tok_wall_pockets_set(0, 0);
  tok_wall_pocket_release(&pocket[0]);
  tok_wall_pocket_release(&pocket[1]);
  gkyl_nodal_ops_n2m(
    n2m, &inp->cbasis, &inp->cgrid, nrange, &up->local, 3, up->geo_int.mc2p_nodal, up->geo_int.mc2p,
    true
  );
  gkyl_nodal_ops_n2m(
    n2m, &inp->cbasis, &inp->cgrid, nrange, &up->local, 1, up->geo_int.bmag_nodal, up->geo_int.bmag,
    true
  );
  gkyl_nodal_ops_n2m(
    n2m, &inp->cbasis, &inp->cgrid, nrange, &up->local, 1, qprofile_nodal, up->geo_int.qprofile,
    true
  );
  gkyl_nodal_ops_release(n2m);
  gkyl_array_release(qprofile_nodal);

  gkyl_free(arc_memo);
  gkyl_free(arc_memo_left);
  gkyl_free(arc_memo_right);
  gkyl_free(sep_trace_r);
  gkyl_free(sep_trace_z);
  gkyl_free(sep_trace_s);
  gkyl_free(arc_ctx.ext_ladder_w);
  gkyl_free(arc_ctx.ext_ladder_rf);
  gkyl_free(ordered_trace_storage);
}

void
gkyl_tok_geo_calc_surface(
  struct gk_geometry *up, int dir, struct gkyl_range *nrange, double dzc[3],
  struct gkyl_tok_geo *geo, struct gkyl_tok_geo_grid_inp *inp,
  struct gkyl_position_map *position_map
)
{
  // Vessel-outline policy for this block, decided once.
  const bool enforce_wall = tok_wall_enforced(geo, inp);

  geo->rleft = inp->rleft;
  geo->rright = inp->rright;

  geo->inexact_roots = inp->inexact_roots;

  geo->rmax = inp->rmax;
  geo->rmin = inp->rmin;

  enum { PSI_IDX, AL_IDX, TH_IDX }; // Arrangement of computational coordinates.
  enum { X_IDX, Y_IDX, Z_IDX }; // Arrangement of cartesian coordinates.

  double dtheta = inp->cgrid.dx[TH_IDX], dpsi = inp->cgrid.dx[PSI_IDX],
         dalpha = inp->cgrid.dx[AL_IDX];

  double theta_lo = up->grid.lower[TH_IDX] +
                    (up->local.lower[TH_IDX] - up->global.lower[TH_IDX]) * up->grid.dx[TH_IDX],
         psi_lo = up->grid.lower[PSI_IDX] +
                  (up->local.lower[PSI_IDX] - up->global.lower[PSI_IDX]) * up->grid.dx[PSI_IDX],
         alpha_lo = up->grid.lower[AL_IDX] +
                    (up->local.lower[AL_IDX] - up->global.lower[AL_IDX]) * up->grid.dx[AL_IDX];

  double dels[2] = {1.0 / sqrt(3), 1.0 - 1.0 / sqrt(3)};
  theta_lo += dir == 2 ? 0.0 : dels[1] * dtheta / 2.0;
  // The pockets of the block's declared plates (see the corner pass), set only
  // around this pass's containment test.
  const bool plate_end[2] = {
    enforce_wall && up->local.lower[TH_IDX] == up->global.lower[TH_IDX] &&
      tok_wall_theta_end_on_declared_plate(inp, geo, 0),
    enforce_wall && up->local.upper[TH_IDX] == up->global.upper[TH_IDX] &&
      tok_wall_theta_end_on_declared_plate(inp, geo, 1)
  };
  struct tok_wall_pocket pocket[2] = {{0}};
  if (enforce_wall) {
    tok_wall_block_pockets(inp, geo, plate_end, 0, nrange, pocket);
  }
  psi_lo += dir == 0 ? 0.0 : dels[1] * dpsi / 2.0;
  alpha_lo += dir == 1 ? 0. : dels[1] * dalpha / 2.0;

  double dx_fact = 1.0 / up->basis.poly_order;
  dtheta *= dx_fact;
  dpsi *= dx_fact;
  dalpha *= dx_fact;

  // Used for finite differences.
  double delta_alpha = dalpha * 1e-2;
  double delta_psi = dpsi * 1e-2;
  double delta_theta = dtheta * 1e-2;
  dzc[0] = delta_psi;
  dzc[1] = delta_alpha;
  dzc[2] = delta_theta;
  int modifiers[5] = {0, -1, 1, -2, 2};

  double rclose = inp->rclose;
  double rright = inp->rright;
  double rleft = inp->rleft;

  int nzcells;
  if (geo->use_cubics) {
    nzcells = geo->rzgrid_cubic.cells[1];
  } else {
    nzcells = geo->rzgrid.cells[1];
  }
  double *arc_memo = gkyl_malloc(sizeof(double[nzcells]));
  double *arc_memo_left = gkyl_malloc(sizeof(double[nzcells]));
  double *arc_memo_right = gkyl_malloc(sizeof(double[nzcells]));
  // Trace buffers: see tok_sep_trace_capacity.
  int sep_trace_capacity = tok_sep_trace_capacity(inp, nzcells);
  double *sep_trace_r = gkyl_malloc(sizeof(double[sep_trace_capacity]));
  double *sep_trace_z = gkyl_malloc(sizeof(double[sep_trace_capacity]));
  double *sep_trace_s = gkyl_malloc(sizeof(double[sep_trace_capacity]));
  // The far-boundary, correspondence and map traces share one allocation.
  double *ordered_trace_storage = gkyl_malloc(sizeof(double) * 8 * (size_t)sep_trace_capacity);

  struct arc_length_ctx arc_ctx = {
    .geo = geo,
    .arc_memo = arc_memo,
    .arc_memo_right = arc_memo_right,
    .arc_memo_left = arc_memo_left,
    .sep_trace_r = sep_trace_r,
    .sep_trace_z = sep_trace_z,
    .sep_trace_s = sep_trace_s,
    .sep_trace_capacity = sep_trace_capacity,
    .far_trace_r = ordered_trace_storage,
    .far_trace_z = ordered_trace_storage + sep_trace_capacity,
    .far_trace_s = ordered_trace_storage + 2 * sep_trace_capacity,
    .trace_corr_v = ordered_trace_storage + 3 * sep_trace_capacity,
    .map_trace_r = ordered_trace_storage + 4 * sep_trace_capacity,
    .map_trace_z = ordered_trace_storage + 5 * sep_trace_capacity,
    .map_trace_s = ordered_trace_storage + 6 * sep_trace_capacity,
    .map_trace_phi = ordered_trace_storage + 7 * sep_trace_capacity,
    .ext_ladder_w = NULL,
    .ext_ladder_rf = NULL,
    .ftype = inp->ftype,
    .zmaxis = geo->zmaxis,
  };
  struct plate_ctx pctx = {.geo = geo};
  tok_init_xpt_ray_target(inp, position_map, &arc_ctx);

  int cidx[3] = {0};
  for (int ia = nrange->lower[AL_IDX]; ia < nrange->lower[AL_IDX] + 1; ++ia) {
    cidx[AL_IDX] = ia;
    double alpha_curr = dir == 1 ? alpha_lo + ia * dalpha :
                                   calc_running_coord(alpha_lo, ia - nrange->lower[AL_IDX], dalpha);

    for (int ip = nrange->lower[PSI_IDX]; ip <= nrange->upper[PSI_IDX]; ++ip) {
      int ip_delta_max = 5;
      for (int ip_delta = 0; ip_delta < ip_delta_max; ip_delta++) {
        if ((ip == nrange->lower[PSI_IDX]) &&
            (up->local.lower[PSI_IDX] == up->global.lower[PSI_IDX]) && dir == 0) {
          if (ip_delta == 1 || ip_delta == 3) {
            continue; // one sided stencils at edge
          }
        } else if ((ip == nrange->upper[PSI_IDX]) &&
                   (up->local.upper[PSI_IDX] == up->global.upper[PSI_IDX]) && dir == 0) {
          if (ip_delta == 2 || ip_delta == 4) {
            continue; // one sided stencils at edge
          }
        } else { // interior
          if (ip_delta == 3 || ip_delta == 4) {
            continue;
          }
        }

        double psi_curr = dir == 0 ? psi_lo + ip * dpsi :
                                     calc_running_coord(psi_lo, ip - nrange->lower[PSI_IDX], dpsi);
        psi_curr += modifiers[ip_delta] * delta_psi;
        // Non-uniform psi. Finite differences are calculated in calc_metric.c
        double Psi_curr;
        position_map->maps[0](0.0, &psi_curr, &Psi_curr, position_map->ctxs[0]);
        double dPsi_dpsi =
          gkyl_position_map_slope(position_map, 0, psi_curr, delta_psi, ip, nrange);
        psi_curr = Psi_curr;

        double darcL, arcL_curr, arcL_lo;

        // For double null blocks this should set arc_ctx :
        // zmin, zmax, rclose, arcL_tot for all blocks. No left and right
        // For a full core case:
        // also set phi_right and arcL_right
        // For a single null case:
        // also set zmin_left and zmin_right
        if (tok_xpt_ordered_placement(inp)) {
          tok_prepare_ordered_map(inp, &arc_ctx, psi_curr);
        } else {
          tok_find_endpoints(
            inp, geo, &arc_ctx, &pctx, psi_curr, alpha_curr, arc_memo, arc_memo_left, arc_memo_right
          );
        }

        darcL =
          (arc_ctx.arc_hi - arc_ctx.arc_lo) / (up->basis.poly_order * inp->cgrid.cells[TH_IDX]);
        // at the beginning of each theta loop we need to reset things
        cidx[PSI_IDX] = ip;
        arcL_curr = 0.0;
        arcL_lo = tok_arc_from_theta(inp, &arc_ctx, theta_lo);
        double ridders_min, ridders_max;

        for (int it = nrange->lower[TH_IDX]; it <= nrange->upper[TH_IDX]; ++it) {
          arcL_curr = dir == 2 ? arcL_lo + it * darcL :
                                 calc_running_coord(arcL_lo, it - nrange->lower[TH_IDX], darcL);
          double theta_curr = tok_theta_from_arc(inp, &arc_ctx, arcL_curr);

          // Calculate derivatives using finite difference for ddtheta,
          // as well as transform the computational coordiante to the non-uniform field-aligned value

          // We cannot do non-uniform alpha because we are modeling axisymmetric systems
          // Non-uniform theta
          double Theta_curr;
          position_map->maps[2](0.0, &theta_curr, &Theta_curr, position_map->ctxs[2]);
          double dTheta_dtheta =
            gkyl_position_map_slope(position_map, 2, theta_curr, delta_theta, it, nrange);
          theta_curr = Theta_curr;

          struct tok_ordered_point ordered = {0.0};
          bool ordered_mapping =
            tok_ordered_map_lookup(inp, &arc_ctx, theta_curr, alpha_curr, &ordered);
          if (tok_xpt_ordered_placement(inp) && !ordered_mapping) {
            fprintf(
              stderr, "TOK_ORDERED_MAP lookup failed ftype=%d psi=%.17g theta=%.17g\n", inp->ftype,
              psi_curr, theta_curr
            );
            abort();
          }
          double r_curr = 0.0, z_curr = 0.0, phi_curr = 0.0;
          double drdz_curr = 0.0, dr_curr = 0.0, dz_curr = 0.0;
          if (ordered_mapping) {
            r_curr = ordered.r;
            z_curr = ordered.z;
            phi_curr = ordered.phi;
          } else {
            arcL_curr = tok_xpt_theta_to_arc(inp, &arc_ctx, theta_curr);

            tok_set_ridders(inp, &arc_ctx, psi_curr, arcL_curr, &rclose, &ridders_min, &ridders_max);

            // single-null SOL rows: the node at exact arc (tok_lsn_exact_row)
            double lsn_r = 0.0, lsn_z = 0.0;
            const bool lsn_exact =
              arc_ctx.arcL_tot > 0.0 &&
              tok_lsn_exact_point(geo, psi_curr, arcL_curr / arc_ctx.arcL_tot, &lsn_r, &lsn_z);
            struct gkyl_qr_res res = lsn_exact ?
                                       (struct gkyl_qr_res){.res = lsn_z} :
                                       gkyl_ridders(
                                         arc_length_func, &arc_ctx, arc_ctx.zmin, arc_ctx.zmax,
                                         ridders_min, ridders_max, geo->root_param.max_iter, 1e-10
                                       );
            if (!lsn_exact) {
              tok_geo_check_arc_root(
                inp, psi_curr, theta_curr, arcL_curr, arc_ctx.zmin, arc_ctx.zmax, ridders_min,
                ridders_max, &res
              );
            }
            z_curr = res.res;
            ((struct gkyl_tok_geo *)geo)->stat.nroot_cont_calls += res.nevals;

            double sep_r_curr = 0.0;
            bool at_sep_trace =
              tok_half_domain_sep_rz(inp, &arc_ctx, theta_curr, &sep_r_curr, &z_curr);

            // Only a surface normal to theta contains true theta-boundary
            // nodes.  On radial and alpha surfaces the theta coordinates are
            // Gauss points and must not be snapped to the seam.
            double fixed_root_ref = rclose;
            bool at_fixed_edge =
              dir == TH_IDX &&
              tok_xpt_at_fixed_edge(
                inp, &arc_ctx, it, nrange, up->local.lower[TH_IDX] == up->global.lower[TH_IDX],
                up->local.upper[TH_IDX] == up->global.upper[TH_IDX], &z_curr, &fixed_root_ref
              );
            bool at_xpt_anchor = dir == TH_IDX && arc_ctx.xpt_anchor_valid &&
                                 tok_xpt_at_seam(
                                   inp, it, nrange,
                                   up->local.lower[TH_IDX] == up->global.lower[TH_IDX],
                                   up->local.upper[TH_IDX] == up->global.upper[TH_IDX]
                                 );
            if (at_xpt_anchor) {
              z_curr = arc_ctx.xpt_anchor_z;
            }

            double R[4] = {0}, dRdZ[4] = {0};
            double dR[4] = {0}, dZ[4] = {0};
            int nr = gkyl_tok_geo_R_psiZ(geo, psi_curr, z_curr, 4, R, dRdZ, dR, dZ);
            double root_ref =
              at_sep_trace ?
                sep_r_curr :
                (at_xpt_anchor ? arc_ctx.xpt_anchor_r : (at_fixed_edge ? fixed_root_ref : rclose));
            r_curr = choose_closest(root_ref, R, R, nr);
            drdz_curr = choose_closest(root_ref, R, dRdZ, nr);
            dr_curr = choose_closest(root_ref, R, dR, nr);
            dz_curr = choose_closest(root_ref, R, dZ, nr);
            if (at_xpt_anchor) {
              r_curr = arc_ctx.xpt_anchor_r;
            }
            if (at_sep_trace) {
              r_curr = sep_r_curr;
            }
            if (lsn_exact && z_curr == lsn_z && !at_xpt_anchor && !at_sep_trace) {
              r_curr = lsn_r; // the exact point's own R, on its branch
              if (nr <= 0) {
                nr = 1;
              }
            }

            if (tok_geo_same_flux(psi_curr, geo->psisep) && ip_delta == 0) {
              // Snap to the X point of the representation that evaluates psi.
              const double *rx, *zx;
              const int nx = tok_geo_xpts(geo, &rx, &zx);
              for (int k = 0; k < nx && k < 2; ++k) {
                if (z_curr == zx[k]) {
                  nr = 1;
                  r_curr = rx[k];
                }
              }
            }

            if (nr == 0) {
              printf("ip = %d, it = %d, ia = %d, ip_delta = %d\n", ip, it, ia, ip_delta);
              printf(
                "Block Type = %d | Failed to find a root at psi = %g, Z = %1.16f\n", inp->ftype,
                psi_curr, z_curr
              );
              assert(false);
            }

            phi_curr = phi_func(alpha_curr, z_curr, &arc_ctx);
            if (lsn_exact && r_curr == lsn_r && z_curr == lsn_z) {
              // single-null SOL rows: the angle on the exact row (tok_lsn_exact_phi)
              double lsn_phi = 0.0;
              if (tok_lsn_exact_phi(
                    geo, inp->ftype, psi_curr, arcL_curr / arc_ctx.arcL_tot, r_curr, z_curr,
                    &lsn_phi
                  )) {
                phi_curr = alpha_curr + lsn_phi;
              }
            }
          }
          cidx[TH_IDX] = it;
          int lidx = 0;
          if (ip_delta != 0) {
            lidx = 3 + 3 * (ip_delta - 1);
          }
          double *mc2p_fd_n =
            gkyl_array_fetch(up->geo_surf[dir].mc2p_nodal_fd, gkyl_range_idx(nrange, cidx));
          double *ddtheta_n =
            gkyl_array_fetch(up->geo_surf[dir].ddtheta_nodal, gkyl_range_idx(nrange, cidx));
          double *ddpsi_n =
            gkyl_array_fetch(up->geo_surf[dir].ddpsi_nodal, gkyl_range_idx(nrange, cidx));
          double *bmag_n =
            gkyl_array_fetch(up->geo_surf[dir].bmag_nodal, gkyl_range_idx(nrange, cidx));
          double *curlbhat_n =
            gkyl_array_fetch(up->geo_surf[dir].curlbhat_nodal, gkyl_range_idx(nrange, cidx));
          double *deltats_n =
            gkyl_array_fetch(up->geo_surf[dir].deltats_nodal, gkyl_range_idx(nrange, cidx));

          mc2p_fd_n[lidx + X_IDX] = r_curr;
          mc2p_fd_n[lidx + Y_IDX] = z_curr;
          mc2p_fd_n[lidx + Z_IDX] = phi_curr;

          if (ip_delta == 0) {
            // Gate ONLY the containment test. This block also writes the surface
            // metric quantities below; gating the whole block skips them and
            // yields J=0 at check_right_handed.
            if (enforce_wall) {
              double wall_point[2] = {r_curr, z_curr};
              // Theta-face nodes at the block's ends were put ON the plate by
              // the root finder; the plate was checked against the outline
              // with its pocket. Every other face node is judged against the
              // outline and the pockets.
              const bool on_plate = dir == 2 &&
                                    ((it == nrange->lower[TH_IDX] && plate_end[0] && pocket[0].n) ||
                                     (it == nrange->upper[TH_IDX] && plate_end[1] && pocket[1].n));
              tok_wall_pockets_set(&pocket[0], &pocket[1]);
              bool inside = on_plate || tok_wall_point_inside(geo->efit, wall_point);
              tok_wall_pockets_set(0, 0);
              if (!inside) {
                fprintf(
                  stderr,
                  "TOK_GEO_WALL_DOMAIN_FAILED ftype=%d scope=face dir=%d rz=(%.17g,%.17g)\n",
                  inp->ftype, dir, r_curr, z_curr
                );
                // the point's radial coordinate, as this loop placed it
                double x = dir == 0 ? psi_lo + ip * dpsi :
                                      calc_running_coord(psi_lo, ip - nrange->lower[PSI_IDX], dpsi);
                if (!tok_wall_trial_record_where(
                      false, false,
                      tok_wall_trial_on_movable_side(
                        x, up->grid.lower[PSI_IDX], up->grid.upper[PSI_IDX]
                      )
                    )) {
                  abort();
                }
              }
            }
            if (ordered_mapping) {
              ddtheta_n[0] = ordered.dr_dtheta * dTheta_dtheta;
              ddtheta_n[1] = ordered.dz_dtheta * dTheta_dtheta;
              ddtheta_n[2] = ordered.dphi_dtheta * dTheta_dtheta;
            } else {
              double darc_dtheta = arc_ctx.arc_interval_valid ?
                                     arc_ctx.arc_darc_dtheta :
                                     (arc_ctx.xpt_map_valid ? arc_ctx.xpt_map_darc_dtheta :
                                                              arc_ctx.arcL_tot / (2.0 * M_PI));
              ddtheta_n[0] = sin(atan2(dr_curr, dz_curr)) * darc_dtheta * dTheta_dtheta;
              ddtheta_n[1] = cos(atan2(dr_curr, dz_curr)) * darc_dtheta * dTheta_dtheta;
              ddtheta_n[2] = dphidtheta_func(z_curr, &arc_ctx) * dTheta_dtheta;
            }
            ddpsi_n[0] = dPsi_dpsi;
            bmag_n[0] = bmag_func(r_curr, z_curr, &arc_ctx);
            curlbhat_func(psi_curr, r_curr, z_curr, phi_curr, curlbhat_n, &arc_ctx);
            deltats_n[0] = phi_curr - alpha_curr;
          }
        }
      }
    }
  }

  // Populate other alpha indices by using axisymmetry
  for (int ia = nrange->lower[AL_IDX] + 1; ia <= nrange->upper[AL_IDX]; ++ia) {
    cidx[AL_IDX] = ia;
    double alpha_curr = dir == 1 ? alpha_lo + ia * dalpha :
                                   calc_running_coord(alpha_lo, ia - nrange->lower[AL_IDX], dalpha);
    double alpha_donor = dir == 1 ? alpha_lo + nrange->lower[AL_IDX] * dalpha :
                                    calc_running_coord(alpha_lo, 0, dalpha);
    double alpha_diff = alpha_curr - alpha_donor;

    for (int ip = nrange->lower[PSI_IDX]; ip <= nrange->upper[PSI_IDX]; ++ip) {
      cidx[PSI_IDX] = ip;
      int ip_delta_max = 5;
      for (int ip_delta = 0; ip_delta < ip_delta_max; ip_delta++) {
        if ((ip == nrange->lower[PSI_IDX]) &&
            (up->local.lower[PSI_IDX] == up->global.lower[PSI_IDX]) && dir == 0) {
          if (ip_delta == 1 || ip_delta == 3) {
            continue; // one sided stencils at edge
          }
        } else if ((ip == nrange->upper[PSI_IDX]) &&
                   (up->local.upper[PSI_IDX] == up->global.upper[PSI_IDX]) && dir == 0) {
          if (ip_delta == 2 || ip_delta == 4) {
            continue; // one sided stencils at edge
          }
        } else { // interior
          if (ip_delta == 3 || ip_delta == 4) {
            continue;
          }
        }
        double psi_curr = dir == 0 ? psi_lo + ip * dpsi :
                                     calc_running_coord(psi_lo, ip - nrange->lower[PSI_IDX], dpsi);
        psi_curr += modifiers[ip_delta] * delta_psi;

        for (int it = nrange->lower[TH_IDX]; it <= nrange->upper[TH_IDX]; ++it) {
          cidx[TH_IDX] = it;
          int lidx = 0;
          if (ip_delta != 0) {
            lidx = 3 + 3 * (ip_delta - 1);
          }

          double *mc2p_fd_n =
            gkyl_array_fetch(up->geo_surf[dir].mc2p_nodal_fd, gkyl_range_idx(nrange, cidx));
          double *ddtheta_n =
            gkyl_array_fetch(up->geo_surf[dir].ddtheta_nodal, gkyl_range_idx(nrange, cidx));
          double *ddpsi_n =
            gkyl_array_fetch(up->geo_surf[dir].ddpsi_nodal, gkyl_range_idx(nrange, cidx));
          double *bmag_n =
            gkyl_array_fetch(up->geo_surf[dir].bmag_nodal, gkyl_range_idx(nrange, cidx));
          double *curlbhat_n =
            gkyl_array_fetch(up->geo_surf[dir].curlbhat_nodal, gkyl_range_idx(nrange, cidx));
          double *deltats_n =
            gkyl_array_fetch(up->geo_surf[dir].deltats_nodal, gkyl_range_idx(nrange, cidx));

          int donor_cidx[3];
          donor_cidx[AL_IDX] = nrange->lower[AL_IDX];
          donor_cidx[PSI_IDX] = ip;
          donor_cidx[TH_IDX] = it;

          double *donor_mc2p_fd_n =
            gkyl_array_fetch(up->geo_surf[dir].mc2p_nodal_fd, gkyl_range_idx(nrange, donor_cidx));
          double *donor_ddtheta_n =
            gkyl_array_fetch(up->geo_surf[dir].ddtheta_nodal, gkyl_range_idx(nrange, donor_cidx));
          double *donor_ddpsi_n =
            gkyl_array_fetch(up->geo_surf[dir].ddpsi_nodal, gkyl_range_idx(nrange, donor_cidx));
          double *donor_bmag_n =
            gkyl_array_fetch(up->geo_surf[dir].bmag_nodal, gkyl_range_idx(nrange, donor_cidx));

          mc2p_fd_n[lidx + X_IDX] = donor_mc2p_fd_n[lidx + X_IDX];
          mc2p_fd_n[lidx + Y_IDX] = donor_mc2p_fd_n[lidx + Y_IDX];
          mc2p_fd_n[lidx + Z_IDX] = donor_mc2p_fd_n[lidx + Z_IDX] + alpha_diff;
          if (ip_delta == 0) {
            ddtheta_n[0] = donor_ddtheta_n[0];
            ddtheta_n[1] = donor_ddtheta_n[1];
            ddtheta_n[2] = donor_ddtheta_n[2];
            ddpsi_n[0] = donor_ddpsi_n[0];
            bmag_n[0] = donor_bmag_n[0];
            curlbhat_func(
              psi_curr, mc2p_fd_n[X_IDX], mc2p_fd_n[Y_IDX], mc2p_fd_n[Z_IDX], curlbhat_n, &arc_ctx
            );
            deltats_n[0] = mc2p_fd_n[Z_IDX] - alpha_curr;
          }
        }
      }
    }
  }

  gkyl_free(arc_memo);
  gkyl_free(arc_memo_left);
  gkyl_free(arc_memo_right);
  gkyl_free(sep_trace_r);
  gkyl_free(sep_trace_z);
  gkyl_free(sep_trace_s);
  gkyl_free(arc_ctx.ext_ladder_w);
  gkyl_free(arc_ctx.ext_ladder_rf);
  gkyl_free(ordered_trace_storage);
  tok_wall_pocket_release(&pocket[0]);
  tok_wall_pocket_release(&pocket[1]);
}

void
gkyl_tok_geo_set_extent(
  struct gkyl_tok_geo_grid_inp *inp, struct gkyl_tok_geo *geo, double *theta_lo, double *theta_up
)
{
  tok_geo_set_extent(inp, geo, theta_lo, theta_up);
}

struct gkyl_tok_geo_stat
gkyl_tok_geo_get_stat(const struct gkyl_tok_geo *geo)
{
  return geo->stat;
}

bool
gkyl_tok_geo_uses_extended_construction(const struct gkyl_tok_geo_grid_inp *inp)
{
  return tok_ext_construction(inp);
}

void
gkyl_tok_geo_release(struct gkyl_tok_geo *geo)
{
  gkyl_array_release(geo->psiRZ);
  gkyl_array_release(geo->psiRZ_cubic);
  if (geo->psi_cell_bounds) {
    gkyl_array_release(geo->psi_cell_bounds);
  }
  if (geo->psi_block_bounds) {
    gkyl_free(geo->psi_block_bounds);
  }
  if (geo->psi_cell_bounds_cubic) {
    gkyl_array_release(geo->psi_cell_bounds_cubic);
  }
  if (geo->psi_block_bounds_cubic) {
    gkyl_free(geo->psi_block_bounds_cubic);
  }
  gkyl_array_release(geo->fpoldg);
  gkyl_array_release(geo->fpolprimedg);
  gkyl_array_release(geo->qdg);
  for (int slot = 0; slot < 2; ++slot) {
    gkyl_free((void *)geo->divertor_wall[slot].segments);
  }
  gkyl_efit_release(geo->efit);
  gkyl_free(geo);
}
