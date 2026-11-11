#include <gkyl_cocos.h>
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
#include <gkyl_dg_bin_ops.h>

#include <math.h>
#include <string.h>

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
    if (actx->pre == true) {
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

// Once contour's end points and branch (i.e. left or right branch of psi) are
// chosen, this function verifies traversal against the local flux gradient on a regular branch.
// This includes private-flux contours: a global outward-flux flag alone is not a
// sufficient check of their selected branch. The probe does not alter the map context.
static void
tok_verify_contour_orientation(struct gkyl_tok_geo_grid_inp *inp, struct arc_length_ctx *ctx)
{
  struct arc_length_ctx probe = *ctx;
  double close, flo, fhi;
  tok_set_ridders(inp, &probe, probe.psi, 0.371 * probe.arcL_tot, &close, &flo, &fhi);
  struct gkyl_qr_res root = gkyl_ridders(
    arc_length_func, &probe, probe.zmin, probe.zmax, flo, fhi, probe.geo->root_param.max_iter, 1e-10
  );
  double radius[4], slope[4], tangent_R[4], tangent_Z[4];
  int count =
    gkyl_tok_geo_R_psiZ(probe.geo, probe.psi, root.res, 4, radius, slope, tangent_R, tangent_Z);
  if (!count) {
    fprintf(stderr, "TOKAMAK: cannot determine contour orientation.\n");
    abort();
  }
  double field_Z = choose_closest(close, radius, tangent_Z, count);
  int reference_Z_sign = probe.right ? 1 : -1;
  if (reference_Z_sign * probe.geo->theta_orientation * field_Z <= 0) {
    fprintf(
      stderr, "TOKAMAK: contour branch orientation disagrees with flux (block %d, psi=%g).\n",
      inp->ftype, probe.psi
    );
    abort();
  }
}

// Function to calculate phi given alpha
double
phi_func(double alpha_curr, double Z, void *ctx)
{
  struct arc_length_ctx *actx = ctx;
  double *arc_memo = actx->arc_memo;
  double psi = actx->psi, rclose = actx->rclose, zmin = actx->zmin, arcL = actx->arcL,
         zmax = actx->zmax;

  // alpha labels a field line. Reference locations may differ between blocks;
  // the oriented integral must nevertheless give B dot grad(alpha) = 0.
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

  return alpha_curr + actx->geo->theta_orientation * ival + phi_ref;
}

// Signed q on closed contours; full-connection toroidal winding on open contours.
double
qprofile_func(void *ctx)
{
  struct arc_length_ctx *actx = ctx;
  const struct gkyl_tok_geo *geo = actx->geo;
  if (geo->efit->num_xpts > 0 && actx->psi == geo->psisep) {
    return NAN;
  }

  // The endpoint context spans the complete branch, even for split blocks.
  // Integrals are positive arc integrals; COCOS supplies the sign below.
  double integral = 0;
  switch (actx->ftype) {
    case GKYL_GEOMETRY_TOKAMAK_CORE:
    case GKYL_GEOMETRY_TOKAMAK_CORE_L:
    case GKYL_GEOMETRY_TOKAMAK_CORE_R:
      integral = integrate_phi_along_psi_contour_memo(
                   geo, actx->psi, actx->zmin, actx->zmax, actx->rright, false, false, NULL
                 ) +
                 integrate_phi_along_psi_contour_memo(
                   geo, actx->psi, actx->zmin, actx->zmax, actx->rleft, false, false, NULL
                 );
      break;
    case GKYL_GEOMETRY_TOKAMAK_IWL:
      integral = integrate_phi_along_psi_contour_memo(
        geo, actx->psi, actx->zmin, actx->zmax, actx->rright, false, false, NULL
      );
      if (geo->theta_orientation * (actx->psi - geo->sibry) <= 0) {
        integral += integrate_phi_along_psi_contour_memo(
          geo, actx->psi, actx->zmin, actx->zmax, actx->rleft, false, false, NULL
        );
      } else {
        integral +=
          integrate_phi_along_psi_contour_memo(
            geo, actx->psi, actx->zmin_iwl_plate, actx->zmax, actx->rleft, false, false, NULL
          ) +
          integrate_phi_along_psi_contour_memo(
            geo, actx->psi, actx->zmin, actx->zmax_iwl_plate, actx->rleft, false, false, NULL
          );
      }
      break;
    case GKYL_GEOMETRY_TOKAMAK_LSN_SOL:
    case GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO:
    case GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID:
    case GKYL_GEOMETRY_TOKAMAK_LSN_SOL_UP:
    case GKYL_GEOMETRY_TOKAMAK_PF_LO_L:
    case GKYL_GEOMETRY_TOKAMAK_PF_LO_R:
      integral = integrate_phi_along_psi_contour_memo(
                   geo, actx->psi, actx->zmin_right, actx->zmax, actx->rright, false, false, NULL
                 ) +
                 integrate_phi_along_psi_contour_memo(
                   geo, actx->psi, actx->zmin_left, actx->zmax, actx->rleft, false, false, NULL
                 );
      break;
    case GKYL_GEOMETRY_TOKAMAK_PF_UP_L:
    case GKYL_GEOMETRY_TOKAMAK_PF_UP_R:
      integral = integrate_phi_along_psi_contour_memo(
                   geo, actx->psi, actx->zmin, actx->zmax_right, actx->rright, false, false, NULL
                 ) +
                 integrate_phi_along_psi_contour_memo(
                   geo, actx->psi, actx->zmin, actx->zmax_left, actx->rleft, false, false, NULL
                 );
      break;
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_UP:
      integral = integrate_phi_along_psi_contour_memo(
        geo, actx->psi, actx->zmin, actx->zmax, actx->rclose, false, false, NULL
      );
      break;
    default:
      return NAN;
  }
  double psi_fpol = actx->psi;
  if (psi_fpol < geo->fgrid.lower[0] || psi_fpol > geo->fgrid.upper[0]) {
    psi_fpol = geo->sibry;
  }
  int idx = fmin(
    geo->frange.lower[0] + (int)floor((psi_fpol - geo->fgrid.lower[0]) / geo->fgrid.dx[0]),
    geo->frange.upper[0]
  );
  const double *coeffs = gkyl_array_cfetch(geo->fpoldg, gkyl_range_idx(&geo->frange, &idx));
  double center;
  gkyl_rect_grid_cell_center(&geo->fgrid, &idx, &center);
  double eta = 2 * (psi_fpol - center) / geo->fgrid.dx[0];
  double fpol = geo->fbasis.eval_expand(&eta, coeffs);
  return geo->theta_orientation * fpol * integral / (2 * M_PI);
}

void
gkyl_tok_geo_calc_qprofile(
  struct gkyl_tok_geo *geo, const struct gkyl_tok_geo_grid_inp *inp,
  const struct gkyl_rect_grid *grid, const struct gkyl_basis *basis, struct gkyl_array *qprofile
)
{
  assert(grid->ndim == 1 && basis->ndim == 1);

  struct gkyl_tok_geo_grid_inp ginp = *inp;
  geo->rleft = inp->rleft;
  geo->rright = inp->rright;
  geo->rmin = inp->rmin;
  geo->rmax = inp->rmax;
  geo->inexact_roots = inp->inexact_roots;
  int nz = geo->use_cubics ? geo->rzgrid_cubic.cells[1] : geo->rzgrid.cells[1];

  double *memo = gkyl_malloc(3 * nz * sizeof(double));

  struct arc_length_ctx actx = {
    .geo = geo,
    .ftype = inp->ftype,
    .zmaxis = geo->zmaxis,
    .arc_memo = memo,
    .arc_memo_left = memo + nz,
    .arc_memo_right = memo + 2 * nz,
  };
  struct plate_ctx pctx = {.geo = geo};
  int nquad = basis->poly_order + 1;
  gkyl_array_clear(qprofile, 0.0);

  for (int cell = 0; cell < grid->cells[0]; ++cell) {
    double *coeff = gkyl_array_fetch(qprofile, cell);
    for (int node = 0; node < nquad; ++node) {
      double eta = gkyl_gauss_ordinates[nquad][node];
      double psi = grid->lower[0] + (cell + 0.5 * (1 + eta)) * grid->dx[0];

      tok_find_endpoints(
        &ginp, geo, &actx, &pctx, psi, 0.0, actx.arc_memo, actx.arc_memo_left, actx.arc_memo_right
      );

      double qval = qprofile_func(&actx), values[basis->num_basis];
      basis->eval(&eta, values);

      for (int k = 0; k < basis->num_basis; ++k) {
        coeff[k] += gkyl_gauss_weights[nquad][node] * qval * values[k];
      }
    }
  }
  gkyl_free(memo);
}

double
gkyl_tok_geo_calc_Cy(
  struct gkyl_tok_geo *geo, const struct gkyl_tok_geo_grid_inp *inp,
  const struct gkyl_rect_grid *grid, int poly_order, double psi_ref
)
{
  assert(psi_ref >= grid->lower[0] && psi_ref <= grid->upper[0]);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 1, poly_order);

  struct gkyl_array *profile = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, grid->cells[0]);
  gkyl_tok_geo_calc_qprofile(geo, inp, grid, &basis, profile);

  int cell = fmin((int)((psi_ref - grid->lower[0]) / grid->dx[0]), grid->cells[0] - 1);
  double eta = 2 * (psi_ref - grid->lower[0]) / grid->dx[0] - 2 * cell - 1;
  double qref = basis.eval_expand(&eta, gkyl_array_cfetch(profile, cell));
  gkyl_array_release(profile);

  double roots[16], slope[16], tangent_R[16], tangent_Z[16];
  int nroots =
    gkyl_tok_geo_R_psiZ(geo, psi_ref, geo->zmaxis, 16, roots, slope, tangent_R, tangent_Z);

  double Raxis = geo->efit->rmaxis, Romp = HUGE_VAL;
  for (int i = 0; i < nroots; ++i) {
    if (roots[i] > Raxis && roots[i] < Romp) {
      Romp = roots[i];
    }
  }
  return (Romp - Raxis) / qref;
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
  integrand = integrand * actx->arcL_tot / 2 / M_PI;
  return integrand;
}

// Supply derivatives of the field-aligned contour map, including the position-map chain rule.
static void
tok_contour_derivatives(
  double Zvert, double tangent_R, double tangent_Z, double dTheta_dtheta,
  struct arc_length_ctx *ctx, double deriv[3]
)
{
  double norm = hypot(tangent_R, tangent_Z);
  assert(norm > 0 && dTheta_dtheta > 0);
  double scale = ctx->arcL_tot / (2 * M_PI) * dTheta_dtheta;
  deriv[0] = tangent_R / norm * scale;
  deriv[1] = tangent_Z / norm * scale;
  deriv[2] = dphidtheta_func(Zvert, ctx) * dTheta_dtheta;
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

    Br = actx->geo->cocos_psi_factor * dpsidZ / r_curr;
    Bz = -actx->geo->cocos_psi_factor * dpsidR / r_curr;
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
    Br = actx->geo->cocos_psi_factor * dpsidZ / r_curr;
    Bz = -actx->geo->cocos_psi_factor * dpsidR / r_curr;
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
  double fpolprime = psi < actx->geo->fgrid.lower[0] || psi > actx->geo->fgrid.upper[0] ?
                       0.0 :
                       actx->geo->fbasis.eval_expand(&fx, coeffs);

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

  Br = actx->geo->cocos_psi_factor * dpsidZ / r_curr;
  Bz = -actx->geo->cocos_psi_factor * dpsidR / r_curr;
  bmag = sqrt(Br * Br + Bz * Bz + Bphi * Bphi);

  dBrdR = actx->geo->cocos_psi_factor * d2psidRdZ / r_curr - Br / r_curr;
  dBrdZ = actx->geo->cocos_psi_factor * d2psidZ2 / r_curr;
  dBzdR = -actx->geo->cocos_psi_factor * d2psidR2 / r_curr - Bz / r_curr;
  dBzdZ = -actx->geo->cocos_psi_factor * d2psidRdZ / r_curr;

  double dFdR = fpolprime * dpsidR;
  double dFdZ = fpolprime * dpsidZ;
  dBdR =
    1 / bmag * (Br * dBrdR + Bz * dBzdR + fpol / r_curr * (dFdR / r_curr - fpol / r_curr / r_curr));
  dBdZ = 1 / bmag * (Br * dBrdZ + Bz * dBzdZ + fpol / r_curr * dFdZ / r_curr);

  // Get contravariant R-Z-phi components of curlbhat (upperscript components on tangent basis).
  double polar_comp[3] = {0.0};
  polar_comp[0] = -dFdZ / (bmag * r_curr) + dBdZ * fpol / (bmag * bmag * r_curr); // R component ^1
  polar_comp[1] = (dBrdZ - dBzdR) / (bmag * r_curr) +
                  (dBdR * Bz - dBdZ * Br) / (bmag * bmag * r_curr); // Phi component ^2
  polar_comp[2] = dFdR / (bmag * r_curr) - dBdR * fpol / (bmag * bmag * r_curr);

  // Convert to cartesian
  curlbhat[0] = polar_comp[0] * cos(phi) - polar_comp[1] * sin(phi) * r_curr;
  curlbhat[1] = polar_comp[0] * sin(phi) + polar_comp[1] * cos(phi) * r_curr;
  curlbhat[2] = polar_comp[2];
}

struct gkyl_tok_geo *
gkyl_tok_geo_new(const struct gkyl_efit_inp *inp, const struct gkyl_tok_geo_grid_inp *ginp)
{
  struct gkyl_tok_geo *geo = gkyl_malloc(sizeof(*geo));
  *geo = (struct gkyl_tok_geo){};

  geo->efit = gkyl_efit_new(inp);
  const struct gkyl_cocos *cocos = &gkyl_cocos_tokamak;
  geo->cocos_psi_factor =
    cocos->sigma_RvarphiZ_cocos * cocos->sigma_Bp_cocos / pow(2 * M_PI, cocos->e_Bp_cocos);
  geo->theta_orientation = geo->efit->sibry > geo->efit->simag ? 1 : -1;

  geo->plate_spec = ginp->plate_spec;
  geo->plate_func_lower = ginp->plate_func_lower;
  geo->plate_func_upper = ginp->plate_func_upper;

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
  geo->psisep = geo->efit->psisep;
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

void
gkyl_tok_geo_calc(
  struct gk_geometry *up, struct gkyl_range *nrange, struct gkyl_tok_geo *geo,
  struct gkyl_tok_geo_grid_inp *inp, struct gkyl_position_map *position_map
)
{
  geo->rleft = inp->rleft;
  geo->rright = inp->rright;

  geo->inexact_roots = inp->inexact_roots;

  geo->rmax = inp->rmax;
  geo->rmin = inp->rmin;

  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates

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

  struct arc_length_ctx arc_ctx = {
    .geo = geo,
    .arc_memo = arc_memo,
    .arc_memo_right = arc_memo_right,
    .arc_memo_left = arc_memo_left,
    .ftype = inp->ftype,
    .zmaxis = geo->zmaxis,
  };
  struct plate_ctx pctx = {.geo = geo};

  int cidx[3] = {0};
  for (int ia = nrange->lower[AL_IDX]; ia < nrange->lower[AL_IDX] + 1; ++ia) {
    cidx[AL_IDX] = ia;
    double alpha_curr = gk_geometry_node_coord(
      &up->grid, &up->local, &up->global, nrange, 1, ia, up->basis.poly_order, false
    );

    for (int ip = nrange->lower[PSI_IDX]; ip <= nrange->upper[PSI_IDX]; ++ip) {
      double psi_curr = gk_geometry_node_coord(
        &up->grid, &up->local, &up->global, nrange, 0, ip, up->basis.poly_order, false
      );

      // Non-uniform psi. Finite differences are calculated in calc_metric.c
      double psi_curr_mapped;
      position_map->maps[0](0.0, &psi_curr, &psi_curr_mapped, position_map->ctxs[0]);
      psi_curr = psi_curr_mapped;

      double arcL_curr;

      // For double null blocks this should set arc_ctx :
      // zmin, zmax, rclose, arcL_tot for all blocks. No left and right
      // For a full core case:
      // also set phi_right and arcL_right
      // For a single null case:
      // also set zmin_left and zmin_right
      tok_find_endpoints(
        inp, geo, &arc_ctx, &pctx, psi_curr, alpha_curr, arc_memo, arc_memo_left, arc_memo_right
      );
      tok_verify_contour_orientation(inp, &arc_ctx);

      // At the beginning of each theta loop we need to reset things.
      cidx[PSI_IDX] = ip;
      arcL_curr = 0.0;

      double ridders_min, ridders_max;
      // Set node coordinates.
      for (int it = nrange->lower[TH_IDX]; it <= nrange->upper[TH_IDX]; ++it) {
        int it_delta = 0;

        double theta_curr = gk_geometry_node_coord(
          &up->grid, &up->local, &up->global, nrange, 2, it, up->basis.poly_order, false
        );

        // Calculate derivatives using finite difference for ddtheta,
        // as well as transform the computational coordiante to the non-uniform field-aligned value

        // We cannot do non-uniform alpha because we are modeling axisymmetric systems
        // Non-uniform theta
        double Theta_curr;
        position_map->maps[2](0.0, &theta_curr, &Theta_curr, position_map->ctxs[2]);
        theta_curr = Theta_curr;
        arcL_curr = (geo->theta_orientation * theta_curr + M_PI) / (2 * M_PI / arc_ctx.arcL_tot);

        tok_set_ridders(inp, &arc_ctx, psi_curr, arcL_curr, &rclose, &ridders_min, &ridders_max);

        struct gkyl_qr_res res = gkyl_ridders(
          arc_length_func, &arc_ctx, arc_ctx.zmin, arc_ctx.zmax, ridders_min, ridders_max,
          geo->root_param.max_iter, 1e-10
        );
        double z_curr = res.res;
        ((struct gkyl_tok_geo *)geo)->stat.nroot_cont_calls += res.nevals;

        bool at_lower = it == nrange->lower[TH_IDX] &&
                        up->local.lower[TH_IDX] == up->global.lower[TH_IDX];
        bool at_upper = it == nrange->upper[TH_IDX] &&
                        up->local.upper[TH_IDX] == up->global.upper[TH_IDX];
        bool at_ref_lower = geo->theta_orientation > 0 ? at_lower : at_upper;
        bool at_ref_upper = geo->theta_orientation > 0 ? at_upper : at_lower;

        if (inp->half_domain) { //Alternative for half domain
          if (psi_curr == geo->psisep) {
            if (at_ref_upper) {
              if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_R ||
                  inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_L ||
                  inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO ||
                  inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID) {
                z_curr = geo->efit->Zxpt[0];
              } else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO ||
                         inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID ||
                         inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE) {
                z_curr = geo->efit->Zxpt[0];
              }
            }
            if (at_ref_lower) {
              if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_L ||
                  inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R ||
                  inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID ||
                  inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO) {
                z_curr = geo->efit->Zxpt[0];
              } else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_UP ||
                         inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID ||
                         inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE) {
                z_curr = geo->efit->Zxpt[0];
              }
            }
          }
        } else { // For full domain
          if (psi_curr == geo->psisep) {
            if (at_ref_upper) {
              if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_UP_L ||
                  inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R ||
                  inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID ||
                  inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_UP) {
                z_curr = geo->efit->Zxpt[1];
              } else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_R ||
                         inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_L ||
                         inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO ||
                         inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID) {
                z_curr = geo->efit->Zxpt[0];
              } else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO ||
                         inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID ||
                         inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE) {
                z_curr = geo->efit->Zxpt[0];
              }
            }
            if (at_ref_lower) {
              if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_UP_R ||
                  inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_L ||
                  inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP ||
                  inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID) {
                z_curr = geo->efit->Zxpt[1];
              } else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_L ||
                         inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R ||
                         inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID ||
                         inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO) {
                z_curr = geo->efit->Zxpt[0];
              } else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_UP ||
                         inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID ||
                         inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE) {
                z_curr = geo->efit->Zxpt[0];
              }
            }
          }
        }

        // Ensure that node at the lower Xpt (lower end of right core/upper end of left core) is at the same location
        if (psi_curr != geo->psisep) {
          if (at_ref_upper && inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_L) {
            z_curr = arc_ctx.zmin;
          }
          if (at_ref_lower && inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R) {
            z_curr = arc_ctx.zmin;
          }
        }

        double R[4] = {0}, dRdZ[4] = {0};
        double dR[4] = {0}, dZ[4] = {0};
        int nr = gkyl_tok_geo_R_psiZ(geo, psi_curr, z_curr, 4, R, dRdZ, dR, dZ);
        double r_curr = choose_closest(rclose, R, R, nr);
        double drdz_curr = choose_closest(rclose, R, dRdZ, nr);
        double dr_curr = choose_closest(rclose, R, dR, nr);
        double dz_curr = choose_closest(rclose, R, dZ, nr);

        if (psi_curr == geo->psisep) {
          if (z_curr == geo->efit->Zxpt[0]) {
            nr = 1;
            r_curr = geo->efit->Rxpt[0];
          }
          if (z_curr == geo->efit->Zxpt[1]) {
            nr = 1;
            r_curr = geo->efit->Rxpt[1];
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

        cidx[TH_IDX] = it;
        int lidx = 0;

        double phi_curr = phi_func(alpha_curr, z_curr, &arc_ctx);
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

  // Populate other alpha indices by using axisymmetry
  for (int ia = nrange->lower[AL_IDX] + 1; ia <= nrange->upper[AL_IDX]; ++ia) {
    cidx[AL_IDX] = ia;
    double alpha_curr = gk_geometry_node_coord(
      &up->grid, &up->local, &up->global, nrange, 1, ia, up->basis.poly_order, false
    );
    double alpha_donor = gk_geometry_node_coord(
      &up->grid, &up->local, &up->global, nrange, 1, nrange->lower[AL_IDX], up->basis.poly_order,
      false
    );
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

  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&inp->cbasis, &inp->cgrid, false);
  gkyl_nodal_ops_n2m(
    n2m, &inp->cbasis, &inp->cgrid, nrange, &up->local, 3, up->geo_corn.mc2p_nodal,
    up->geo_corn.mc2p, false
  );
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
}

void
gkyl_tok_geo_calc_interior(
  struct gk_geometry *up, struct gkyl_range *nrange, double dzc[3], struct gkyl_tok_geo *geo,
  struct gkyl_tok_geo_grid_inp *inp, struct gkyl_position_map *position_map
)
{
  geo->rleft = inp->rleft;
  geo->rright = inp->rright;

  geo->inexact_roots = inp->inexact_roots;

  geo->rmax = inp->rmax;
  geo->rmin = inp->rmin;

  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates

  double dtheta = inp->cgrid.dx[TH_IDX], dpsi = inp->cgrid.dx[PSI_IDX],
         dalpha = inp->cgrid.dx[AL_IDX];

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

  struct arc_length_ctx arc_ctx = {
    .geo = geo,
    .arc_memo = arc_memo,
    .arc_memo_right = arc_memo_right,
    .arc_memo_left = arc_memo_left,
    .ftype = inp->ftype,
    .zmaxis = geo->zmaxis,
  };
  struct plate_ctx pctx = {.geo = geo};

  // Temporary array to store nodal q profile.
  struct gkyl_array *qprofile_nodal =
    gkyl_array_new(GKYL_DOUBLE, up->geo_int.bmag_nodal->ncomp, up->geo_int.bmag_nodal->size);

  int cidx[3] = {0};
  for (int ia = nrange->lower[AL_IDX]; ia < nrange->lower[AL_IDX] + 1; ++ia) {
    cidx[AL_IDX] = ia;
    double alpha_curr = gk_geometry_node_coord(
      &up->grid, &up->local, &up->global, nrange, 1, ia, up->basis.poly_order + 1, true
    );

    for (int ip = nrange->lower[PSI_IDX]; ip <= nrange->upper[PSI_IDX]; ++ip) {
      int ip_delta_max = 3;
      for (int ip_delta = 0; ip_delta < ip_delta_max; ip_delta++) {
        double psi_curr =
          gk_geometry_node_coord(
            &up->grid, &up->local, &up->global, nrange, 0, ip, up->basis.poly_order + 1, true
          ) +
          modifiers[ip_delta] * delta_psi;

        // Non-uniform psi. Finite differences are calculated in calc_metric.c
        double Psi_curr;
        position_map->maps[0](0.0, &psi_curr, &Psi_curr, position_map->ctxs[0]);
        double dPsi_dpsi = gkyl_position_map_slope(
          position_map, 0, psi_curr, delta_psi, up->grid.lower[0], up->grid.upper[0]
        );
        psi_curr = Psi_curr;

        double arcL_curr;

        // For double null blocks this should set arc_ctx :
        // zmin, zmax, rclose, arcL_tot for all blocks. No left and right
        // For a full core case:
        // also set phi_right and arcL_right
        // For a single null case:
        // also set zmin_left and zmin_right
        tok_find_endpoints(
          inp, geo, &arc_ctx, &pctx, psi_curr, alpha_curr, arc_memo, arc_memo_left, arc_memo_right
        );
        tok_verify_contour_orientation(inp, &arc_ctx);

        // Diagnostic signed safety factor; independent of the theta subdomain.
        double qprofile = ip_delta == 0 ? qprofile_func(&arc_ctx) : 0.0;

        // at the beginning of each theta loop we need to reset things
        cidx[PSI_IDX] = ip;
        arcL_curr = 0.0;

        double ridders_min, ridders_max;

        for (int it = nrange->lower[TH_IDX]; it <= nrange->upper[TH_IDX]; ++it) {
          double theta_curr = gk_geometry_node_coord(
            &up->grid, &up->local, &up->global, nrange, 2, it, up->basis.poly_order + 1, true
          );

          // Calculate derivatives using finite difference for ddtheta,
          // as well as transform the computational coordiante to the non-uniform field-aligned value

          // We cannot do non-uniform alpha because we are modeling axisymmetric systems
          // Non-uniform theta
          double Theta_curr;
          position_map->maps[2](0.0, &theta_curr, &Theta_curr, position_map->ctxs[2]);
          double dTheta_dtheta = gkyl_position_map_slope(
            position_map, 2, theta_curr, delta_theta, up->grid.lower[2], up->grid.upper[2]
          );
          theta_curr = Theta_curr;
          arcL_curr = (geo->theta_orientation * theta_curr + M_PI) / (2 * M_PI / arc_ctx.arcL_tot);

          tok_set_ridders(inp, &arc_ctx, psi_curr, arcL_curr, &rclose, &ridders_min, &ridders_max);

          struct gkyl_qr_res res = gkyl_ridders(
            arc_length_func, &arc_ctx, arc_ctx.zmin, arc_ctx.zmax, ridders_min, ridders_max,
            geo->root_param.max_iter, 1e-10
          );
          double z_curr = res.res;
          ((struct gkyl_tok_geo *)geo)->stat.nroot_cont_calls += res.nevals;

          double R[4] = {0}, dRdZ[4] = {0};
          double dR[4] = {0}, dZ[4] = {0};
          int nr = gkyl_tok_geo_R_psiZ(geo, psi_curr, z_curr, 4, R, dRdZ, dR, dZ);
          double r_curr = choose_closest(rclose, R, R, nr);
          double drdz_curr = choose_closest(rclose, R, dRdZ, nr);
          double dr_curr = choose_closest(rclose, R, dR, nr);
          double dz_curr = choose_closest(rclose, R, dZ, nr);

          if (psi_curr == geo->psisep && ip_delta == 0) {
            if (z_curr == geo->efit->Zxpt[0]) {
              nr = 1;
              r_curr = geo->efit->Rxpt[0];
            }
            if (z_curr == geo->efit->Zxpt[1]) {
              nr = 1;
              r_curr = geo->efit->Rxpt[1];
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

          cidx[TH_IDX] = it;
          int lidx = 0;
          if (ip_delta != 0) {
            lidx = 3 + 3 * (ip_delta - 1);
          }

          double phi_curr = phi_func(alpha_curr, z_curr, &arc_ctx);
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
            tok_contour_derivatives(z_curr, dr_curr, dz_curr, dTheta_dtheta, &arc_ctx, ddtheta_n);
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
    double alpha_curr = gk_geometry_node_coord(
      &up->grid, &up->local, &up->global, nrange, 1, ia, up->basis.poly_order + 1, true
    );
    double alpha_donor = gk_geometry_node_coord(
      &up->grid, &up->local, &up->global, nrange, 1, nrange->lower[AL_IDX],
      up->basis.poly_order + 1, true
    );
    double alpha_diff = alpha_curr - alpha_donor;
    for (int ip = nrange->lower[PSI_IDX]; ip <= nrange->upper[PSI_IDX]; ++ip) {
      cidx[PSI_IDX] = ip;
      int ip_delta_max = 3;
      for (int ip_delta = 0; ip_delta < ip_delta_max; ip_delta++) {
        double psi_comp =
          gk_geometry_node_coord(
            &up->grid, &up->local, &up->global, nrange, 0, ip, up->basis.poly_order + 1, true
          ) +
          modifiers[ip_delta] * delta_psi;
        double psi_curr;
        position_map->maps[0](0.0, &psi_comp, &psi_curr, position_map->ctxs[0]);
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
            double *qprofile_n = gkyl_array_fetch(qprofile_nodal, gkyl_range_idx(nrange, cidx));
            const double *donor_qprofile_n =
              gkyl_array_cfetch(qprofile_nodal, gkyl_range_idx(nrange, donor_cidx));
            qprofile_n[0] = donor_qprofile_n[0];
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
}

void
gkyl_tok_geo_calc_surface(
  struct gk_geometry *up, int dir, struct gkyl_range *nrange, double dzc[3],
  struct gkyl_tok_geo *geo, struct gkyl_tok_geo_grid_inp *inp,
  struct gkyl_position_map *position_map
)
{
  geo->rleft = inp->rleft;
  geo->rright = inp->rright;

  geo->inexact_roots = inp->inexact_roots;

  geo->rmax = inp->rmax;
  geo->rmin = inp->rmin;

  enum { PSI_IDX, AL_IDX, TH_IDX }; // Arrangement of computational coordinates.
  enum { X_IDX, Y_IDX, Z_IDX }; // Arrangement of cartesian coordinates.

  double dtheta = inp->cgrid.dx[TH_IDX], dpsi = inp->cgrid.dx[PSI_IDX],
         dalpha = inp->cgrid.dx[AL_IDX];

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

  struct arc_length_ctx arc_ctx = {
    .geo = geo,
    .arc_memo = arc_memo,
    .arc_memo_right = arc_memo_right,
    .arc_memo_left = arc_memo_left,
    .ftype = inp->ftype,
    .zmaxis = geo->zmaxis,
  };
  struct plate_ctx pctx = {.geo = geo};

  int cidx[3] = {0};
  for (int ia = nrange->lower[AL_IDX]; ia < nrange->lower[AL_IDX] + 1; ++ia) {
    cidx[AL_IDX] = ia;
    double alpha_curr = gk_geometry_node_coord(
      &up->grid, &up->local, &up->global, nrange, 1, ia, dir == 1 ? 1 : up->basis.poly_order + 1,
      dir != 1
    );

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

        double psi_curr = gk_geometry_node_coord(
          &up->grid, &up->local, &up->global, nrange, 0, ip,
          dir == 0 ? 1 : up->basis.poly_order + 1, dir != 0
        );
        psi_curr += modifiers[ip_delta] * delta_psi;
        // Non-uniform psi. Finite differences are calculated in calc_metric.c
        double Psi_curr;
        position_map->maps[0](0.0, &psi_curr, &Psi_curr, position_map->ctxs[0]);
        double dPsi_dpsi = gkyl_position_map_slope(
          position_map, 0, psi_curr, delta_psi, up->grid.lower[0], up->grid.upper[0]
        );
        psi_curr = Psi_curr;

        double arcL_curr;

        // For double null blocks this should set arc_ctx :
        // zmin, zmax, rclose, arcL_tot for all blocks. No left and right
        // For a full core case:
        // also set phi_right and arcL_right
        // For a single null case:
        // also set zmin_left and zmin_right
        tok_find_endpoints(
          inp, geo, &arc_ctx, &pctx, psi_curr, alpha_curr, arc_memo, arc_memo_left, arc_memo_right
        );
        tok_verify_contour_orientation(inp, &arc_ctx);

        // at the beginning of each theta loop we need to reset things
        cidx[PSI_IDX] = ip;
        arcL_curr = 0.0;

        double ridders_min, ridders_max;

        for (int it = nrange->lower[TH_IDX]; it <= nrange->upper[TH_IDX]; ++it) {
          double theta_curr = gk_geometry_node_coord(
            &up->grid, &up->local, &up->global, nrange, 2, it,
            dir == 2 ? 1 : up->basis.poly_order + 1, dir != 2
          );

          // Calculate derivatives using finite difference for ddtheta,
          // as well as transform the computational coordiante to the non-uniform field-aligned value

          // We cannot do non-uniform alpha because we are modeling axisymmetric systems
          // Non-uniform theta
          double Theta_curr;
          position_map->maps[2](0.0, &theta_curr, &Theta_curr, position_map->ctxs[2]);
          double dTheta_dtheta = gkyl_position_map_slope(
            position_map, 2, theta_curr, delta_theta, up->grid.lower[2], up->grid.upper[2]
          );
          theta_curr = Theta_curr;
          arcL_curr = (geo->theta_orientation * theta_curr + M_PI) / (2 * M_PI / arc_ctx.arcL_tot);

          tok_set_ridders(inp, &arc_ctx, psi_curr, arcL_curr, &rclose, &ridders_min, &ridders_max);

          struct gkyl_qr_res res = gkyl_ridders(
            arc_length_func, &arc_ctx, arc_ctx.zmin, arc_ctx.zmax, ridders_min, ridders_max,
            geo->root_param.max_iter, 1e-10
          );
          double z_curr = res.res;
          ((struct gkyl_tok_geo *)geo)->stat.nroot_cont_calls += res.nevals;

          // Ensure that node at the lower Xpt (lower end of right core/upper end of left core) is at the same location
          if (psi_curr != geo->psisep) {
            if (it == nrange->upper[TH_IDX] && inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_L &&
                up->local.upper[TH_IDX] == up->global.upper[TH_IDX]) {
              z_curr = arc_ctx.zmin;
            }
            if (it == nrange->lower[TH_IDX] && inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R &&
                up->local.lower[TH_IDX] == up->global.lower[TH_IDX]) {
              z_curr = arc_ctx.zmin;
            }
          }

          double R[4] = {0}, dRdZ[4] = {0};
          double dR[4] = {0}, dZ[4] = {0};
          int nr = gkyl_tok_geo_R_psiZ(geo, psi_curr, z_curr, 4, R, dRdZ, dR, dZ);
          double r_curr = choose_closest(rclose, R, R, nr);
          double drdz_curr = choose_closest(rclose, R, dRdZ, nr);
          double dr_curr = choose_closest(rclose, R, dR, nr);
          double dz_curr = choose_closest(rclose, R, dZ, nr);

          if (psi_curr == geo->psisep && ip_delta == 0) {
            if (z_curr == geo->efit->Zxpt[0]) {
              nr = 1;
              r_curr = geo->efit->Rxpt[0];
            }
            if (z_curr == geo->efit->Zxpt[1]) {
              nr = 1;
              r_curr = geo->efit->Rxpt[1];
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

          cidx[TH_IDX] = it;
          int lidx = 0;
          if (ip_delta != 0) {
            lidx = 3 + 3 * (ip_delta - 1);
          }

          double phi_curr = phi_func(alpha_curr, z_curr, &arc_ctx);
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
            tok_contour_derivatives(z_curr, dr_curr, dz_curr, dTheta_dtheta, &arc_ctx, ddtheta_n);
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
    double alpha_curr = gk_geometry_node_coord(
      &up->grid, &up->local, &up->global, nrange, 1, ia, dir == 1 ? 1 : up->basis.poly_order + 1,
      dir != 1
    );
    double alpha_donor = gk_geometry_node_coord(
      &up->grid, &up->local, &up->global, nrange, 1, nrange->lower[AL_IDX],
      dir == 1 ? 1 : up->basis.poly_order + 1, dir != 1
    );
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
        double psi_comp = gk_geometry_node_coord(
          &up->grid, &up->local, &up->global, nrange, 0, ip,
          dir == 0 ? 1 : up->basis.poly_order + 1, dir != 0
        );
        psi_comp += modifiers[ip_delta] * delta_psi;

        double psi_curr;
        position_map->maps[0](0.0, &psi_comp, &psi_curr, position_map->ctxs[0]);
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

void
gkyl_tok_geo_release(struct gkyl_tok_geo *geo)
{
  gkyl_array_release(geo->psiRZ);
  gkyl_array_release(geo->psiRZ_cubic);
  gkyl_array_release(geo->fpoldg);
  gkyl_array_release(geo->fpolprimedg);
  gkyl_array_release(geo->qdg);
  gkyl_efit_release(geo->efit);
  gkyl_free(geo);
}
