#include <gkyl_efit_priv.h>
#include <float.h>
#include <math.h>
#include <string.h>

static bool newton_raphson_seed(struct gkyl_efit *up, const double *coeffs,
  const double *xinit, double *xsol, bool cubics);

static inline double
efit_sq(double x)
{
  return x*x;
}

static inline double
efit_dist2(double R0, double Z0, double R1, double Z1)
{
  return efit_sq(R0-R1) + efit_sq(Z0-Z1);
}

static inline double
efit_clamp(double x, double lo, double hi)
{
  return x < lo ? lo : (x > hi ? hi : x);
}


// Is an X-point candidate inside the caller-supplied vessel outline?  Candidates are
// tested at -|Z| against the lower-half polygon, since a reflected equilibrium yields
// X-points in +/- pairs.  Always true when no bounds were supplied, so the default
// behaviour is unchanged.
static bool
efit_xpt_in_bounds(const struct gkyl_efit *up, double R, double Z)
{
  if (up->xpt_bound_n <= 0 || !up->xpt_bound_R || !up->xpt_bound_Z)
    return true;

  double z = -fabs(Z);
  bool inside = false;
  for (int i = 0, j = up->xpt_bound_n-1; i < up->xpt_bound_n; j = i++) {
    double Ri = up->xpt_bound_R[i], Zi = up->xpt_bound_Z[i];
    double Rj = up->xpt_bound_R[j], Zj = up->xpt_bound_Z[j];
    // Standard crossing-number test: count edges the rightward ray from (R,z) crosses.
    if ((Zi > z) != (Zj > z)) {
      double Rcross = Ri + (z-Zi)/(Zj-Zi)*(Rj-Ri);
      if (R < Rcross)
        inside = !inside;
    }
  }
  return inside;
}

// Evaluate psi from the quadratic DG representation at (R,Z).  This mirrors
// tok_eval_psi_rz()'s use_cubics==false branch exactly -- same cell lookup, same
// clamping, same logical coordinates -- so a psisep produced here is bit-identical
// to what the geometry code will later evaluate at the same point.
static double
efit_eval_psi_quad(const struct gkyl_efit *up, double R, double Z)
{
  int rzidx[2];
  rzidx[0] = fmin(up->rzlocal.lower[0]
      + (int) floor((R-up->rzgrid.lower[0])/up->rzgrid.dx[0]),
    up->rzlocal.upper[0]);
  rzidx[1] = fmin(up->rzlocal.lower[1]
      + (int) floor((Z-up->rzgrid.lower[1])/up->rzgrid.dx[1]),
    up->rzlocal.upper[1]);
  rzidx[0] = fmax(rzidx[0], up->rzlocal.lower[0]);
  rzidx[1] = fmax(rzidx[1], up->rzlocal.lower[1]);

  long loc = gkyl_range_idx(&up->rzlocal, rzidx);
  const double *coeffs = gkyl_array_cfetch(up->psizr, loc);
  double xc[2], xy[2];
  gkyl_rect_grid_cell_center(&up->rzgrid, rzidx, xc);
  xy[0] = (R-xc[0])/(0.5*up->rzgrid.dx[0]);
  xy[1] = (Z-xc[1])/(0.5*up->rzgrid.dx[1]);
  return up->rzbasis.eval_expand(xy, coeffs);
}

static double
eval_laplacian_expand_2d_tensor_p2(int dir, const double *z, const double *f )
{
  const double z0 = z[0];
  const double z1 = z[1];
  if (dir == 0)
    return 11.25*f[8]*z1*z1+5.809475019311125*f[6]*z1-3.75*f[8]+3.354101966249685*f[4];
  if (dir == 1)
    return 11.25*f[8]*z0*z0+5.809475019311125*f[7]*z0-3.75*f[8]+3.354101966249685*f[5];

  return 0.0; // can't happen, suppresses warning

}

static double
eval_mixedpartial_expand_2d_tensor_p2(const double *z, const double *f )
{
  const double z0 = z[0];
  const double z1 = z[1];
  return 22.5*f[8]*z0*z1+5.809475019311125*f[7]*z1+5.809475019311125*f[6]*z0+1.5*f[3]; 

}

static void print_result(int n, double x[], double dx[], double errx, double errf, int niter)
{
  double x0 = x[0];
  double y0 = x[1];
  if (x0 >= -1 && x0 <= 1 && y0 >= -1 && y0 <= 1 ) {
    printf("x = ");
    for(int i=0; i<n; i++) printf(" %g ",  x[i]);
    printf("\n");
    printf("dx = ");
    for(int i=0; i<n; i++) printf("= %g ",  dx[i]);
    printf("\n");
    printf("errx = %g\n", errx);
    printf("errf = %g\n", errf);
  }

}

bool 
newton_raphson(struct gkyl_efit *up, const double *coeffs, double *xsol, bool cubics)
{
  return newton_raphson_seed(up, coeffs, (double[2]) { 0.0, 0.0 }, xsol, cubics);
}

static bool
newton_raphson_seed(struct gkyl_efit *up, const double *coeffs,
  const double *xinit, double *xsol, bool cubics)
{
  int n = 2;
  double x[2] = {xinit[0],xinit[1]};
  double dx[2] = {0.0,0.0};
  double fjac[2][2];
  double fjac_inv[2][2];
  double fvec[2];
  double p[2];
  int ntrial = 100;
  double errx = 0.0;
  double errf = 0.0;
  for (int i=0;i<n;i++) xsol[i] = 0.0;

  for (int niter = 0; niter < ntrial; niter++) {
    if (cubics) {
      for(int i=0; i<n; i++) fvec[i] = up->rzbasis_cubic.eval_grad_expand(i,x,coeffs);
      fjac[0][0] = up->evf->eval_cubic_laplacian(0,x,coeffs);
      fjac[0][1] = up->evf->eval_cubic_mixedpartial(x,coeffs);
      fjac[1][0] = up->evf->eval_cubic_mixedpartial(x,coeffs);
      fjac[1][1] = up->evf->eval_cubic_laplacian(1,x,coeffs);
    }
    else {
      for(int i=0; i<n; i++) fvec[i] = up->rzbasis.eval_grad_expand(i,x,coeffs);
      fjac[0][0] = eval_laplacian_expand_2d_tensor_p2(0,x,coeffs);
      fjac[0][1] = eval_mixedpartial_expand_2d_tensor_p2(x,coeffs);
      fjac[1][0] = eval_mixedpartial_expand_2d_tensor_p2(x,coeffs);
      fjac[1][1] = eval_laplacian_expand_2d_tensor_p2(1,x,coeffs);
    }
    errf = 0.0;
    for (int i=0;i<n;i++) errf += fvec[i]*fvec[i];
    errf = sqrt(errf);
    if (errf <= 1e-18) {
      for (int i=0;i<n;i++) xsol[i] = x[i];
      return true;
    }
    for (int i=0;i<n;i++) p[i] = -fvec[i];

    double det = fjac[0][0]*fjac[1][1] - fjac[0][1]*fjac[1][0];
    fjac_inv[0][0] = fjac[1][1]/det;
    fjac_inv[1][1] = fjac[0][0]/det;
    fjac_inv[0][1] = -fjac[0][1]/det;
    fjac_inv[1][0] = -fjac[1][0]/det;

    dx[0] = fjac_inv[0][0]*p[0] + fjac_inv[0][1]*p[1];
    dx[1] = fjac_inv[1][0]*p[0] + fjac_inv[1][1]*p[1];

    errx = 0.0;
    for (int i=0;i<n;i++) {
      errx += fabs(dx[i]);
      x[i] += dx[i];
    }

    if (errx<=1e-18) {
      for (int i=0;i<n;i++) xsol[i] = x[i];
      return true;
    }

  }
  return false;

}

static bool
find_quadratic_xpt_near_target(gkyl_efit *up, double Rtar, double Ztar,
  double *Rout, double *Zout, double *psiout, double *metric_out)
{
  // A true X-point often sits within rounding of a cell interface, in which case the
  // piecewise-quadratic patch on *each* side places the critical point just outside
  // itself (e.g. z=+1.026 from one cell, z=-1.101 from its neighbour).  With a
  // rejection tolerance of 1e-10 every candidate is thrown away and the caller falls
  // through to a global scan that picks a spurious saddle several cells away.  Accept
  // solutions slightly outside the reference cell and let the distance-to-target
  // metric below choose; the clamp keeps the evaluation point inside the cell.
  const double eps_ref = 0.25;

  double point[2] = {Rtar, Ztar};
  bool pick_lower[2] = {true, true};
  int base_idx[2] = {0, 0};
  gkyl_rect_grid_find_cell(&up->rzgrid, point, pick_lower, (int[2]) {-1, -1}, base_idx);

  // Two selections are tracked at once: bd_* is the historical nearest-to-target
  // rule and bu_* additionally outranks clamped candidates. bu_* is returned;
  // bd_* is kept to report every equilibrium where the two disagree.
  bool found = false;
  bool bd_clamped = false, bu_clamped = false;
  double bd_m = DBL_MAX, bu_m = DBL_MAX;
  double bd_R = 0.0, bd_Z = 0.0, bd_psi = DBL_MAX;
  double bu_R = 0.0, bu_Z = 0.0, bu_psi = DBL_MAX;

  for (int dR = -2; dR <= 2; dR++) {
    for (int dZ = -2; dZ <= 2; dZ++) {
      int idx[2] = {base_idx[0] + dR, base_idx[1] + dZ};
      if (!gkyl_range_contains_idx(&up->rzlocal, idx))
        continue;

      const double* psi = gkyl_array_cfetch(up->psizr, gkyl_range_idx(&up->rzlocal, idx));

      double xc[2];
      gkyl_rect_grid_cell_center(&up->rzgrid, idx, xc);

      double xinit[2] = {
        2.0*(Rtar - xc[0])/up->rzgrid.dx[0],
        2.0*(Ztar - xc[1])/up->rzgrid.dx[1]
      };
      xinit[0] = efit_clamp(xinit[0], -0.9, 0.9);
      xinit[1] = efit_clamp(xinit[1], -0.9, 0.9);

      double xsol[2];
      bool status = newton_raphson_seed(up, psi, xinit, xsol, false);
      double xr = xsol[0];
      double zr = xsol[1];
      if (!(status && xr >= -1.0-eps_ref && xr <= 1.0+eps_ref &&
            zr >= -1.0-eps_ref && zr <= 1.0+eps_ref))
        continue;

      // Whether this candidate had to be pulled onto its cell boundary.  A
      // clamped candidate is a critical point of no polynomial at all: it is
      // the point on the face nearest one patch's critical point, and grad psi
      // there is not small.  Measured on 204502, whose selected X point sits
      // exactly on Z = -7*dZ = -0.9625 with |grad psi| = 4.4e-3, while a
      // genuine interior saddle of the neighbouring patch (|grad psi| = 2.6e-5,
      // Hessian eigenvalues -0.46/+0.18) sits 16.3 mm away at (0.703675,-0.947)
      // and loses only because it is further from the cubic target.
      bool clamped = (xr < -1.0 || xr > 1.0 || zr < -1.0 || zr > 1.0);

      xr = efit_clamp(xr, -1.0, 1.0);
      zr = efit_clamp(zr, -1.0, 1.0);
      double R0 = up->rzgrid.dx[0]*xr/2.0 + xc[0];
      double Z0 = up->rzgrid.dx[1]*zr/2.0 + xc[1];
      if (!efit_xpt_in_bounds(up, R0, Z0))
        continue;
      // Evaluate through the same lookup the geometry uses, so psi(Rxpt,Zxpt)==psisep
      // holds exactly even when the clamp puts (R0,Z0) on a cell interface and the
      // floor() lookup resolves to the neighbouring cell.
      double psi0 = efit_eval_psi_quad(up, R0, Z0);
      double m = efit_dist2(R0, Z0, Rtar, Ztar);

      if (!found || m < bd_m) {
        bd_m = m; bd_R = R0; bd_Z = Z0; bd_psi = psi0; bd_clamped = clamped;
      }
      // Rank unclamped candidates ahead of clamped ones, and only fall back to
      // distance within the same class.  Where the X point genuinely straddles
      // a cell interface -- the case eps_ref above was widened for -- every
      // candidate is clamped and this is exactly the old distance rule.
      if (!found || (clamped != bu_clamped ? !clamped : m < bu_m)) {
        bu_m = m; bu_R = R0; bu_Z = Z0; bu_psi = psi0; bu_clamped = clamped;
      }
      found = true;
    }
  }

  if (found) {
    // A genuine critical point interior to its own cell is preferred over one
    // that had to be clamped onto a cell face. The distance-only pick is
    // reported where it differs.
    if (bd_R != bu_R || bd_Z != bu_Z)
      fprintf(stderr,
        "EFIT_XPT_CLAMP_PREF target=(%.17g,%.17g) "
        "dist_pick=(%.17g,%.17g) dist_clamped=%d dist_psi=%.17g "
        "unclamped_pick=(%.17g,%.17g) unclamped_psi=%.17g "
        "separation=%.17g dpsi=%.17g active=unclamped\n",
        Rtar, Ztar, bd_R, bd_Z, (int) bd_clamped, bd_psi,
        bu_R, bu_Z, bu_psi,
        sqrt(efit_dist2(bd_R, bd_Z, bu_R, bu_Z)), bu_psi - bd_psi);

    *Rout = bu_R;
    *Zout = bu_Z;
    *psiout = bu_psi;
    *metric_out = bu_m;
  }
  return found;
}

int
find_xpts(gkyl_efit* up, int num_cubic_xpts, const double *Rxpt_cubic,
  const double *Zxpt_cubic, double *Rxpt, double *Zxpt)
{
  up->xpt_diag_quadratic_found = false;
  up->xpt_diag_fallback_to_cubic = false;
  up->xpt_diag_quad_R = 0.0;
  up->xpt_diag_quad_Z = 0.0;
  up->xpt_diag_quad_psi = DBL_MAX;
  up->xpt_diag_quad_dist_cell = DBL_MAX;

  bool have_cubic = num_cubic_xpts > 0;
  double Rtarget_top = 0.0, Ztarget_top = 0.0;
  double Rtarget_bot = 0.0, Ztarget_bot = 0.0;

  if (have_cubic) {
    Rtarget_top = Rxpt_cubic[0];
    Ztarget_top = up->reflect ? fabs(Zxpt_cubic[0]) : Zxpt_cubic[0];
    Rtarget_bot = Rtarget_top;
    Ztarget_bot = -Ztarget_top;
  }

  bool found_any = false;
  bool found_top = false, found_bot = false;
  double Rbest = 0.0, Zbest = 0.0, psibest = DBL_MAX, metric_best = DBL_MAX;
  double Rtop = 0.0, Ztop = 0.0, psitop = DBL_MAX, metric_top = DBL_MAX;
  double Rbot = 0.0, Zbot = 0.0, psibot = DBL_MAX, metric_bot = DBL_MAX;

  if (have_cubic) {
    if (up->reflect) {
      found_top = find_quadratic_xpt_near_target(up, Rtarget_top, Ztarget_top,
        &Rtop, &Ztop, &psitop, &metric_top);
      found_bot = find_quadratic_xpt_near_target(up, Rtarget_bot, Ztarget_bot,
        &Rbot, &Zbot, &psibot, &metric_bot);
      found_any = found_top || found_bot;
    }
    else {
      found_any = find_quadratic_xpt_near_target(up, Rtarget_top, Ztarget_top,
        &Rbest, &Zbest, &psibest, &metric_best);
    }
  }

  if (!found_any) {
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &up->rzlocal);
    while (gkyl_range_iter_next(&iter)) {
      const double* psi = gkyl_array_cfetch(up->psizr, gkyl_range_idx(&up->rzlocal, iter.idx));
      double xsol[2];
      bool status = newton_raphson(up, psi, xsol, false);
      double x0 = xsol[0];
      double y0 = xsol[1];
      if (!(status && x0 >= -1 && x0 <= 1 && y0 >= -1 && y0 <= 1))
        continue;

      double psi0 = up->rzbasis.eval_expand(xsol, psi);
      double xc[2];
      gkyl_rect_grid_cell_center(&up->rzgrid, iter.idx, xc);
      double R0 = up->rzgrid.dx[0]*x0/2.0 + xc[0];
      double Z0 = up->rzgrid.dx[1]*y0/2.0 + xc[1];
      if (!efit_xpt_in_bounds(up, R0, Z0))
        continue;

      double metric = fabs(psi0 - up->sibry);
      if (have_cubic) {
        if (up->reflect)
          metric = Z0 >= 0.0 ? efit_dist2(R0, Z0, Rtarget_top, Ztarget_top) :
            efit_dist2(R0, Z0, Rtarget_bot, Ztarget_bot);
        else
          metric = efit_dist2(R0, Z0, Rtarget_top, Ztarget_top);
      }

      found_any = true;
      if (!up->reflect) {
        if (metric < metric_best) {
          Rbest = R0;
          Zbest = Z0;
          psibest = psi0;
          metric_best = metric;
        }
      }
      else if (Z0 >= 0.0) {
        if (metric < metric_top) {
          Rtop = R0;
          Ztop = Z0;
          psitop = psi0;
          metric_top = metric;
          found_top = true;
        }
      }
      else {
        if (metric < metric_bot) {
          Rbot = R0;
          Zbot = Z0;
          psibot = psi0;
          metric_bot = metric;
          found_bot = true;
        }
      }
    }
  }

  // The quadratic representation can have no critical point inside the vessel even
  // when the cubic one does (the two reps resolve a shallow saddle differently).
  // Rather than report no X-point at all -- which leaves Rxpt/Zxpt unwritten for
  // callers to read -- fall back to the in-bounds cubic location, taking psi through
  // the same lookup the geometry uses so psi(Rxpt,Zxpt)==psisep still holds.
  if (!found_any && have_cubic &&
      efit_xpt_in_bounds(up, Rtarget_top, Ztarget_top)) {
    found_any = true;
    if (up->reflect) {
      found_top = true;
      found_bot = false;
      Rtop = Rtarget_top;
      Ztop = Ztarget_top;
      psitop = efit_eval_psi_quad(up, Rtarget_top, -Ztarget_top);
      metric_top = 0.0;
    }
    else {
      Rbest = Rtarget_top;
      Zbest = Ztarget_top;
      psibest = efit_eval_psi_quad(up, Rtarget_top, Ztarget_top);
      metric_best = 0.0;
    }
  }

  int num_xpts = 0;
  if (found_any) {
    if (up->reflect) {
      if (found_top || found_bot) {
        double Rsel = 0.0, Zsel = 0.0, psisel = DBL_MAX, msel = DBL_MAX;
        if (found_top) {
          Rsel = Rtop;
          Zsel = Ztop;
          psisel = psitop;
          msel = metric_top;
        }
        if (found_bot && (!found_top || metric_bot < msel)) {
          Rsel = Rbot;
          Zsel = Zbot;
          psisel = psibot;
          msel = metric_bot;
        }

        double absZ = fabs(Zsel);
        if (have_cubic) {
          up->xpt_diag_quadratic_found = true;
          up->xpt_diag_quad_R = Rsel;
          up->xpt_diag_quad_Z = Zsel;
          up->xpt_diag_quad_psi = psisel;
          up->xpt_diag_quad_dist_cell = sqrt(
            efit_sq((Rsel - Rtarget_top)/up->rzgrid.dx[0]) +
            efit_sq((absZ - Ztarget_top)/up->rzgrid.dx[1]));

          double tol2 = efit_sq(0.25*up->rzgrid.dx[0]) + efit_sq(0.25*up->rzgrid.dx[1]);
          double d2 = efit_dist2(Rsel, absZ, Rtarget_top, Ztarget_top);
          if (d2 > tol2) {
            Rsel = Rtarget_top;
            absZ = Ztarget_top;
            // psisep must be the separatrix flux *of this (quadratic) DG
            // representation*, so evaluate it at the substituted cubic X-point
            // rather than copying psisep_cubic.  Copying the cubic value leaves
            // psi(Rxpt,Zxpt) != psisep by ~1e-5, which is a larger fraction of a
            // flux-thin block's span than the first radial surface's offset and
            // makes the X-point ray anchor search start past its own target.
            psisel = efit_eval_psi_quad(up, Rtarget_top, -Ztarget_top);
            up->xpt_diag_fallback_to_cubic = true;
          }
        }
        else {
          up->xpt_diag_quadratic_found = true;
          up->xpt_diag_quad_R = Rsel;
          up->xpt_diag_quad_Z = Zsel;
          up->xpt_diag_quad_psi = psisel;
          up->xpt_diag_quad_dist_cell = -1.0;
        }

        num_xpts = 2;
        Rxpt[0] = Rsel;
        Rxpt[1] = Rsel;
        Zxpt[0] = -absZ;
        Zxpt[1] = absZ;
        up->psisep = psisel;
      }
    }
    else {
      if (have_cubic) {
        up->xpt_diag_quadratic_found = true;
        up->xpt_diag_quad_R = Rbest;
        up->xpt_diag_quad_Z = Zbest;
        up->xpt_diag_quad_psi = psibest;
        up->xpt_diag_quad_dist_cell = sqrt(
          efit_sq((Rbest - Rtarget_top)/up->rzgrid.dx[0]) +
          efit_sq((Zbest - Ztarget_top)/up->rzgrid.dx[1]));

        double tol2 = efit_sq(0.25*up->rzgrid.dx[0]) + efit_sq(0.25*up->rzgrid.dx[1]);
        double d2 = efit_dist2(Rbest, Zbest, Rtarget_top, Ztarget_top);
        if (d2 > tol2) {
          Rbest = Rtarget_top;
          Zbest = Ztarget_top;
          // See the reflect branch above: keep psisep consistent with the
          // quadratic representation the geometry actually evaluates.
          psibest = efit_eval_psi_quad(up, Rtarget_top, Ztarget_top);
          up->xpt_diag_fallback_to_cubic = true;
        }
      }
      else {
        up->xpt_diag_quadratic_found = true;
        up->xpt_diag_quad_R = Rbest;
        up->xpt_diag_quad_Z = Zbest;
        up->xpt_diag_quad_psi = psibest;
        up->xpt_diag_quad_dist_cell = -1.0;
      }

      num_xpts = 1;
      Rxpt[0] = Rbest;
      Zxpt[0] = Zbest;
      up->psisep = psibest;
    }
  }
  return num_xpts;
}

int
find_xpts_cubic(gkyl_efit* up, double *Rxpt, double *Zxpt)
{
  bool found_xpt = false;
  double Rsep, Zsep;
  double psisep = DBL_MAX;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->rzlocal_cubic);
  while (gkyl_range_iter_next(&iter)) {
    if ((iter.idx[1] < gkyl_range_shape(&up->rzlocal_cubic,1)/2 + 1) || (!up->reflect)) {
      const double* psi = gkyl_array_cfetch(up->psizr_cubic, gkyl_range_idx(&up->rzlocal_cubic, iter.idx));
      double xsol[2];
      bool status = newton_raphson(up, psi, xsol, true);
      double x0 = xsol[0];
      double y0 = xsol[1];
      double psi0 = up->rzbasis_cubic.eval_expand(xsol, psi);
      if (x0 >= -1 && x0 <= 1 && y0 >= -1 && y0 <= 1 && status) {
        double xc[2];
        gkyl_rect_grid_cell_center(&up->rzgrid_cubic, iter.idx, xc);
        double R0 = up->rzgrid_cubic.dx[0]*x0/2.0 + xc[0];
        double Z0 = up->rzgrid_cubic.dx[1]*y0/2.0 + xc[1];
        if (!efit_xpt_in_bounds(up, R0, Z0))
          continue;
        // Only count a candidate as found once it is actually selectable, so an
        // all-rejected sweep cannot report success with unset Rsep/Zsep.
        if (fabs(psi0 - up->sibry) <= fabs(psisep - up->sibry)) {
          found_xpt = true;
          Rsep = R0;
          Zsep = Z0;
          psisep = psi0;
        }
      }
    }
  }

  int num_xpts = 0;
  if (found_xpt) {
    if (up->reflect) {
      num_xpts = 2;
      Rxpt[0] = Rsep;
      Rxpt[1] = Rsep;
      Zxpt[0] = Zsep;
      Zxpt[1] = -Zsep;
      up->psisep_cubic = psisep;
    }
    else {
      num_xpts = 1;
      Rxpt[0] = Rsep;
      Zxpt[0] = Zsep;
      up->psisep_cubic = psisep;
    }
  }
  return num_xpts;
}


// Every critical point of the flux (grad psi = 0), from the cubic
// representation, classified by the sign of the Hessian determinant: negative
// is a saddle (an X point), positive an extremum (an O point).
//
// The cubic array is built from the nodal data as supplied and is never
// reflected (only the quadratic copy is). With reflect the geometry uses the
// lower half mirrored, so only that half is searched, exactly as
// find_xpts_cubic does: a critical point in the other half is not in any
// domain the grid can cover (204554: an extremum near the upper X point
// refused a shot whose grid never touches that half). The cubic is C1, so a
// critical point on a cell interface is a genuine zero of the patch on either
// side and is found from both; duplicates are dropped.
// Twenty-five Newton seeds per cell rather than one, because an O point and an
// X point 25-50 mm apart share a 33 mm cell on NSTX-U (204502) and a single
// seed at the cell centre converges to only one of them.
//
// Returns the number found; fills the first nmax.
int
gkyl_efit_critical_points(gkyl_efit *up, int nmax, double *Rc, double *Zc,
  double *psic, int *is_opoint)
{
  const double seeds[5] = { -0.8, -0.4, 0.0, 0.4, 0.8 };
  const double dmin = 1e-6*fmin(up->rzgrid_cubic.dx[0], up->rzgrid_cubic.dx[1]);
  int n = 0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->rzlocal_cubic);
  while (gkyl_range_iter_next(&iter)) {
    if (up->reflect && iter.idx[1] >= gkyl_range_shape(&up->rzlocal_cubic, 1)/2 + 1) continue;
    const double *psi = gkyl_array_cfetch(up->psizr_cubic, gkyl_range_idx(&up->rzlocal_cubic, iter.idx));
    double xc[2];
    gkyl_rect_grid_cell_center(&up->rzgrid_cubic, iter.idx, xc);
    for (int a = 0; a < 5; ++a) for (int b = 0; b < 5; ++b) {
      double xinit[2] = { seeds[a], seeds[b] }, xsol[2];
      if (!newton_raphson_seed(up, psi, xinit, xsol, true)) continue;
      if (xsol[0] < -1.0 || xsol[0] > 1.0 || xsol[1] < -1.0 || xsol[1] > 1.0) continue;
      double R0 = xc[0] + 0.5*up->rzgrid_cubic.dx[0]*xsol[0];
      double Z0 = xc[1] + 0.5*up->rzgrid_cubic.dx[1]*xsol[1];
      bool dup = false;
      for (int k = 0; k < n && k < nmax && !dup; ++k)
        dup = efit_dist2(R0, Z0, Rc[k], Zc[k]) <= dmin*dmin;
      if (dup) continue;
      if (n < nmax) {
        double hrr = up->evf->eval_cubic_laplacian(0, xsol, psi);
        double hzz = up->evf->eval_cubic_laplacian(1, xsol, psi);
        double hrz = up->evf->eval_cubic_mixedpartial(xsol, psi);
        Rc[n] = R0; Zc[n] = Z0;
        psic[n] = up->rzbasis_cubic.eval_expand(xsol, psi);
        is_opoint[n] = hrr*hzz - hrz*hrz > 0.0;
      }
      ++n;
    }
  }
  return n;
}

void
get_stripped_filename(const char *filepath, char *out_buffer) {
  const char *last_slash = strrchr(filepath, '/');
  const char *filename_start = (last_slash) ? last_slash + 1 : filepath;

  const char *last_dot = strrchr(filename_start, '.');

  if (!last_dot || last_dot == filename_start) {
    strcpy(out_buffer, filename_start);
  }
  else {
    size_t len = last_dot - filename_start;
    strncpy(out_buffer, filename_start, len);
    out_buffer[len] = '\0';
  }
}
