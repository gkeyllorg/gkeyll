#include <gkyl_tok_geo_priv.h>
#include <assert.h>

// Helper functions for finding turning points when necessary


// This function will set zmax to be the upper turning point location
void find_upper_turning_point(struct gkyl_tok_geo *geo, double psi_curr, double zlo, double *zmax, double tolerance)
{
    double tol = tolerance ? tolerance : 1e-12 ;
    //Find the turning points
    double zlo_last;
    double zup=*zmax;
    zlo_last = zlo;
    double R[4], dRdZ[4];
    double dR[4], dZ[4];
    double Rup[4], dRdZup[4];
    double dRup[4], dZup[4];
    while(true){
      int nlo = R_psiZ(geo, psi_curr, zlo, 4, R, dRdZ, dR, dZ);
      int nup = R_psiZ(geo, psi_curr, zup, 4, Rup, dRdZup, dRup, dZup);
      //printf("nlo, nup = %d %d; zlo, zup = %g %g\n", nlo, nup, zlo, zup);
      if (nup > 0) { // This is for the PF_LO regions. Does not seem to break core_L or core_R
                  // However I need to think thos through more. I think it is ok only when xpt
                  // is known quite precisely.
        *zmax = zup;
        break;
      }
      if (nlo>=1){
        if(fabs(zlo-zup)<tol){
          *zmax = zlo;
          break;
        }
        zlo_last = zlo;
        zlo = (zlo+zup)/2;
      }
      if(nlo==0){
        zup = zlo;
        zlo = zlo_last;
      }
    }
}

// This function will set zmin to be the upper turning point location
void find_lower_turning_point(struct gkyl_tok_geo *geo, double psi_curr, double zup, double *zmin, double tolerance)
{
    double tol = tolerance ? tolerance : 1e-12 ;
    int nup = 0;
    double zlo=*zmin;
    double zup_last = zup;
    double R[4], dRdZ[4];
    double dR[4], dZ[4];
    double Rlo[4], dRdZlo[4];
    double dRlo[4], dZlo[4];
    while(true){
      int nup = R_psiZ(geo, psi_curr, zup, 4, R, dRdZ, dR, dZ);
      int nlo = R_psiZ(geo, psi_curr, zlo, 4, Rlo, dRdZlo, dRlo, dZlo);
      //printf("psi = %g;lo, nup = %d %d; zlo, zup = %g %g\n", psi_curr, nlo, nup, zlo, zup);
      if (nlo > 0) {
        *zmin = zlo;
        break;
      }
      if(nup>=1){
        if(fabs(zlo-zup)<tol){
          *zmin = zup;
          break;
        }
        zup_last = zup;
        zup = (zlo+zup)/2;
      }
      if(nup==0){
        zlo = zup;
        zup = zup_last;
      }
    }

}

// This function will set zmin to be the upper lower point location
void find_lower_turning_point_pf_up(struct gkyl_tok_geo *geo, double psi_curr, double zup, double *zmin)
{
    int nup = 0;
    double zlo=*zmin;
    double zup_last = zup;
    double R[4], dRdZ[4];
    double dR[4], dZ[4];
    double Rlo[4], dRdZlo[4];
    double dRlo[4], dZlo[4];
    while(true){
      int nup = R_psiZ(geo, psi_curr, zup, 4, R, dRdZ, dR, dZ);
      int nlo = R_psiZ(geo, psi_curr, zlo, 4, Rlo, dRlo, dRlo, dZlo);
      //if(nlo==1){
      //  if (Rlo[0] < geo->rleft)
      //    nlo=0;
      //}
      if (nlo > 0) {
        *zmin = zlo;
        break;
      }
      if(nup>=2){
        if(fabs(zlo-zup)<1e-12){
          *zmin = zup;
          break;
        }
        zup_last = zup;
        zup = (zlo+zup)/2;
      }
      if(nup==1){
        zlo = zup;
        zup = zup_last;
      }
      if(nup==0){
        zlo = zup;
        zup = zup_last;
      }
    }
}

// This function will set zmax to be the upper turning point location
void find_upper_turning_point_pf_lo(struct gkyl_tok_geo *geo, double psi_curr, double zlo, double *zmax)
{
    //Find the turning points
    double zlo_last;
    double zup=*zmax;
    zlo_last = zlo;
    double R[4], dRdZ[4];
    double dR[4], dZ[4];
    double Rup[4], dRdZup[4];
    double dRup[4], dZup[4];
    while(true){
      int nlo = R_psiZ(geo, psi_curr, zlo, 4, R, dRdZ, dR, dZ);
      int nup = R_psiZ(geo, psi_curr, zup, 4, Rup, dRdZup, dRup, dZup);
      if(nup==1){
        if (Rup[0] < geo->rleft)
          nup=0;
      }
      if (nup > 0) { // This is for the PF_LO regions. Does not seem to break core_L or core_R
                  // However I need to think thos through more. I think it is ok only when xpt
                  // is known quite precisely.
        *zmax = zup;
        break;
      }
      if (nlo>=2){
        if(fabs(zlo-zup)<1e-12){
          *zmax = zlo;
          break;
        }
        zlo_last = zlo;
        zlo = (zlo+zup)/2;
      }
      if(nlo==1){
        zup = zlo;
        zlo = zlo_last;
      }
      if(nlo==0){
        zup = zlo;
        zlo = zlo_last;
      }
    }
}




// Sets zmax if plate is specified
void set_upper_plate(struct gkyl_tok_geo *geo, struct arc_length_ctx* arc_ctx, struct plate_ctx* pctx, double psi_curr)
{
      double rzplate[2];
      pctx->psi_curr = psi_curr;
      pctx->lower=false;
      double a = 0;
      double b = 1;
      double fa = tok_plate_psi_func(a, pctx);
      double fb = tok_plate_psi_func(b, pctx);
      struct gkyl_qr_res res = gkyl_ridders(tok_plate_psi_func, pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smax = res.res;
      geo->plate_func_upper(smax, rzplate);
      arc_ctx->zmax = rzplate[1];
}

// Sets zmin if plate is specified
void set_lower_plate(struct gkyl_tok_geo *geo, struct arc_length_ctx* arc_ctx, struct plate_ctx* pctx, double psi_curr)
{
      double rzplate[2];
      pctx->psi_curr = psi_curr;
      pctx->lower=true;
      double a = 0;
      double b = 1;
      double fa = tok_plate_psi_func(a, pctx);
      double fb = tok_plate_psi_func(b, pctx);
      struct gkyl_qr_res res = gkyl_ridders(tok_plate_psi_func, pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smin = res.res;
      geo->plate_func_lower(smin, rzplate);
      arc_ctx->zmin = rzplate[1];
}

// Sets zmax if plate is specified
void set_upper_iwl_plate(struct gkyl_tok_geo *geo, struct arc_length_ctx* arc_ctx, struct plate_ctx* pctx, double psi_curr)
{
      double rzplate[2];
      pctx->psi_curr = psi_curr;
      pctx->lower=false;
      double a = 0;
      double b = 1;
      double fa = tok_plate_psi_func(a, pctx);
      double fb = tok_plate_psi_func(b, pctx);
      struct gkyl_qr_res res = gkyl_ridders(tok_plate_psi_func, pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smax = res.res;
      geo->plate_func_upper(smax, rzplate);
      arc_ctx->zmax_iwl_plate = rzplate[1];
      geo->rmin = rzplate[0];
}

// Sets zmin if plate is specified
void set_lower_iwl_plate(struct gkyl_tok_geo *geo, struct arc_length_ctx* arc_ctx, struct plate_ctx* pctx, double psi_curr)
{
      double rzplate[2];
      pctx->psi_curr = psi_curr;
      pctx->lower=true;
      double a = 0;
      double b = 1;
      double fa = tok_plate_psi_func(a, pctx);
      double fb = tok_plate_psi_func(b, pctx);
      struct gkyl_qr_res res = gkyl_ridders(tok_plate_psi_func, pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smin = res.res;
      geo->plate_func_lower(smin, rzplate);
      arc_ctx->zmin_iwl_plate = rzplate[1];
      geo->rmin = rzplate[0];
}

void 
tok_geo_set_extent(struct gkyl_tok_geo_grid_inp* inp, struct gkyl_tok_geo *geo, double *theta_lo, double *theta_up)
{
  geo->rleft = inp->rleft;
  geo->rright = inp->rright;

  geo->inexact_roots = inp->inexact_roots;

  geo->rmax = inp->rmax;
  geo->rmin = inp->rmin;
  int nzcells;
  if(geo->use_cubics)
    nzcells = geo->rzgrid_cubic.cells[1];
  else
    nzcells = geo->rzgrid.cells[1];
  double *arc_memo = gkyl_malloc(sizeof(double[nzcells]));
  double *arc_memo_left = gkyl_malloc(sizeof(double[nzcells]));
  double *arc_memo_right = gkyl_malloc(sizeof(double[nzcells]));

  struct arc_length_ctx arc_ctx = {
    .geo = geo,
    .arc_memo = arc_memo,
    .arc_memo_right = arc_memo_right,
    .arc_memo_left = arc_memo_left,
    .ftype = inp->ftype,
    .zmaxis = geo->zmaxis
  };
  struct plate_ctx pctx = {
    .geo = geo
  };

  double del = 1.0e-14;

  if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT || inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO || inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID || inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP) {
    // Immediately set rclose
    arc_ctx.rclose = inp->rright;
    // Set zmin and zmax either fixed or with plate
    if (geo->plate_spec){
      set_upper_plate(geo, &arc_ctx, &pctx, geo->psisep);
      set_lower_plate(geo, &arc_ctx, &pctx, geo->psisep);
    }
    else{
      arc_ctx.zmin = inp->zmin;
      arc_ctx.zmax = inp->zmax;
    }
    double zxpt_up = geo->use_cubics ? geo->efit->Zxpt_cubic[1] : geo->efit->Zxpt[1];
    double zxpt_lo = geo->use_cubics ? geo->efit->Zxpt_cubic[0] : geo->efit->Zxpt[0];
    // Set the arc length
    double arcL_tot = integrate_psi_contour_memo(geo, geo->psisep, arc_ctx.zmin, arc_ctx.zmax, arc_ctx.rclose, false, false, arc_memo);
    double arcL_lo = integrate_psi_contour_memo(geo, geo->psisep, arc_ctx.zmin, zxpt_lo, arc_ctx.rclose, false, false, arc_memo);
    double arcL_mid = integrate_psi_contour_memo(geo, geo->psisep, zxpt_lo, zxpt_up, arc_ctx.rclose, false, false, arc_memo);
    double arcL_up = integrate_psi_contour_memo(geo, geo->psisep, zxpt_up, arc_ctx.zmax, arc_ctx.rclose, false, false, arc_memo);
    if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT) {
      *theta_lo = -M_PI+del;
      *theta_up = M_PI-del;
    }
    else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO) {
      *theta_lo = -M_PI+del;
      *theta_up = -M_PI+del + arcL_lo/arcL_tot*2.0*M_PI;
    }
    else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID) {
      *theta_lo = -M_PI+del + arcL_lo/arcL_tot*2.0*M_PI;
      *theta_up = inp->half_domain ? 0.0 : M_PI-del - arcL_up/arcL_tot*2.0*M_PI;
    }
    else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP) {
      *theta_lo = M_PI-del - arcL_up/arcL_tot*2.0*M_PI;
      *theta_up = M_PI-del;
    }
  }

  else if(inp->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN || inp->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO || inp->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID || inp->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_UP){
    // Immediately set rclose
    arc_ctx.rclose = inp->rleft;
    // Set zmin and zmax either fixed or with plate
    if (geo->plate_spec){
      set_upper_plate(geo, &arc_ctx, &pctx, geo->psisep);
      set_lower_plate(geo, &arc_ctx, &pctx, geo->psisep);
    }
    else{
      arc_ctx.zmin = inp->zmin;
      arc_ctx.zmax = inp->zmax;
    }
    double zxpt_up = geo->use_cubics ? geo->efit->Zxpt_cubic[1] : geo->efit->Zxpt[1];
    double zxpt_lo = geo->use_cubics ? geo->efit->Zxpt_cubic[0] : geo->efit->Zxpt[0];
    // Set the arc Length
    double arcL_tot = integrate_psi_contour_memo(geo, geo->psisep, arc_ctx.zmin, arc_ctx.zmax, arc_ctx.rclose, false, false, arc_memo);
    double arcL_lo = integrate_psi_contour_memo(geo, geo->psisep, arc_ctx.zmin, zxpt_lo, arc_ctx.rclose, false, false, arc_memo);
    double arcL_mid = integrate_psi_contour_memo(geo, geo->psisep, zxpt_lo, zxpt_up, arc_ctx.rclose, false, false, arc_memo);
    double arcL_up = integrate_psi_contour_memo(geo, geo->psisep, zxpt_up, arc_ctx.zmax, arc_ctx.rclose, false, false, arc_memo);
    if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN) {
      *theta_lo = -M_PI+del;
      *theta_up = M_PI-del;
    }
    else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_UP) {
      *theta_lo = -M_PI+del;
      *theta_up = -M_PI+del + arcL_lo/arcL_tot*2.0*M_PI;
    }
    else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID) {
      *theta_lo = inp->half_domain ? 0.0 : -M_PI+del + arcL_lo/arcL_tot*2.0*M_PI;
      *theta_up = M_PI-del - arcL_up/arcL_tot*2.0*M_PI;
    }
    else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO) {
      *theta_lo = M_PI-del - arcL_up/arcL_tot*2.0*M_PI;
      *theta_up = M_PI-del;
    }
  }
  else if(inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE || inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R || inp->ftype ==  GKYL_GEOMETRY_TOKAMAK_CORE_L){
    // Immediately set rleft and rright. Will need both
    arc_ctx.rright = inp->rright;
    arc_ctx.rleft = inp->rleft;

    double zxpt_up = geo->use_cubics ? geo->efit->Zxpt_cubic[1] : geo->efit->Zxpt[1];
    double zxpt_lo = geo->use_cubics ? geo->efit->Zxpt_cubic[0] : geo->efit->Zxpt[0];
    arc_ctx.zmax = inp->zmax ? inp->zmax : zxpt_up; // Initial guess.
                                                  // zmax is specified for single null full core
    double zlo = geo->zmaxis;
    find_upper_turning_point(geo, geo->psisep, zlo, &arc_ctx.zmax, 0);
    arc_ctx.zmin = zxpt_lo; // Initial guess
    double zup = geo->zmaxis;
    find_lower_turning_point(geo, geo->psisep, zup, &arc_ctx.zmin, 0);
    // Done finding turning points
    arc_ctx.right = true;
    double arcL_r = integrate_psi_contour_memo(geo, geo->psisep, arc_ctx.zmin, arc_ctx.zmax, arc_ctx.rright,
      false, false, arc_memo_right);
    arc_ctx.right = false;
    double arcL_l = integrate_psi_contour_memo(geo, geo->psisep, arc_ctx.zmin, arc_ctx.zmax, arc_ctx.rleft,
      false, false, arc_memo_left);
    double arcL_tot = arcL_l + arcL_r;

    if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE) {
      *theta_lo = -M_PI+del;
      *theta_up = M_PI-del;
    }
    else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R) {
      *theta_lo = -M_PI+del;
      *theta_up = inp->half_domain? -M_PI+del + arcL_r/arcL_tot*2.0*M_PI/2.0 : -M_PI+del + arcL_r/arcL_tot*2.0*M_PI;
    }
    else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_L) {
      *theta_lo = inp->half_domain ? M_PI-del - arcL_l/arcL_tot*2.0*M_PI/2.0 : M_PI-del - arcL_l/arcL_tot*2.0*M_PI;
      *theta_up = M_PI-del;
    }
  }
  else if(inp->ftype==GKYL_GEOMETRY_TOKAMAK_LSN_SOL || inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO || inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID || inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_UP){
    // Immediately set rleft and rright. Will need both
    arc_ctx.rright = inp->rright;
    arc_ctx.rleft = inp->rleft;

    //Find the  upper turning point
    arc_ctx.zmax = inp->zmax; // Initial guess
    double zlo = fmax(inp->zmin_left, inp->zmin_right);
    find_upper_turning_point(geo, geo->psisep, zlo, &arc_ctx.zmax, 0);

    // Set zmin left and zmin right wither with plate or fixed
    // This one can't be used with the general func for setting upper and lower plates because it uses zmin left and zmin right
    if (geo->plate_spec){
      double rzplate[2];
      pctx.psi_curr = geo->psisep;
      pctx.lower=false;
      double a = 0;
      double b = 1;
      double fa = tok_plate_psi_func(a, &pctx);
      double fb = tok_plate_psi_func(b, &pctx);
      struct gkyl_qr_res res = gkyl_ridders(tok_plate_psi_func, &pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smax = res.res;
      geo->plate_func_upper(smax, rzplate);
      arc_ctx.zmin_left = rzplate[1];

      pctx.lower=true;
      a = 0;
      b = 1;
      fa = tok_plate_psi_func(a, &pctx);
      fb = tok_plate_psi_func(b, &pctx);
      res = gkyl_ridders(tok_plate_psi_func, &pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smin = res.res;
      geo->plate_func_lower(smin, rzplate);
      arc_ctx.zmin_right = rzplate[1];
    }
    else{
      arc_ctx.zmin_left = inp->zmin_left;
      arc_ctx.zmin_right = inp->zmin_right;
    }

    // Done finding turning points
    double zxpt = geo->use_cubics ? geo->efit->Zxpt_cubic[0] : geo->efit->Zxpt[0];

    arc_ctx.right = true;
    double arcL_mid_r = integrate_psi_contour_memo(geo, geo->psisep, zxpt, arc_ctx.zmax, arc_ctx.rright,
      false, false, arc_memo_right);
    double arcL_lo = integrate_psi_contour_memo(geo, geo->psisep, arc_ctx.zmin_right, zxpt, arc_ctx.rright,
      false, false, arc_memo_right);
    arc_ctx.right = false;
    double arcL_mid_l = integrate_psi_contour_memo(geo, geo->psisep, zxpt, arc_ctx.zmax, arc_ctx.rleft,
      false, false, arc_memo_right);
    double arcL_up = integrate_psi_contour_memo(geo, geo->psisep, arc_ctx.zmin_left, zxpt, arc_ctx.rleft,
      false, false, arc_memo_left);
    double arcL_tot = arcL_lo + arcL_mid_l + arcL_mid_r + arcL_up;

    if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL) {
      *theta_lo = -M_PI+del;
      *theta_up = M_PI-del;
    }
    else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO) {
      *theta_lo = -M_PI+del;
      *theta_up = -M_PI+del + arcL_lo/arcL_tot*2.0*M_PI;
    }
    else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID) {
      *theta_lo = -M_PI+del + arcL_lo/arcL_tot*2.0*M_PI;
      *theta_up = M_PI+del - arcL_up/arcL_tot*2.0*M_PI;
    }
    else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_UP) {
      *theta_lo = M_PI+del - arcL_up/arcL_tot*2.0*M_PI;
      *theta_up = M_PI-del;
    }
  }

  else if(inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_R || inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_L){
    arc_ctx.rright = inp->rright;
    arc_ctx.rleft = inp->rleft;

    //Find the  upper turning point to set zmax
    double zxpt_lo = geo->use_cubics ? geo->efit->Zxpt_cubic[0] : geo->efit->Zxpt[0];
    arc_ctx.zmax = zxpt_lo; // Initial guess
    double zlo = fmin(inp->zmin_left, inp->zmin_right);
    find_upper_turning_point(geo, geo->psisep, zlo, &arc_ctx.zmax, 1e-15);

    // Set zmin left and zmin right wither with plate or fixed
    // This one can't be used with the general func for setting upper and lower plates because it uses zmin left and zmin right
    if (geo->plate_spec){
      double rzplate[2];
      pctx.psi_curr = geo->psisep;
      pctx.lower=false;
      double a = 0;
      double b = 1;
      double fa = tok_plate_psi_func(a, &pctx);
      double fb = tok_plate_psi_func(b, &pctx);
      struct gkyl_qr_res res = gkyl_ridders(tok_plate_psi_func, &pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smax = res.res;
      geo->plate_func_upper(smax, rzplate);
      arc_ctx.zmin_left = rzplate[1];


      pctx.lower=true;
      a = 0;
      b = 1;
      fa = tok_plate_psi_func(a, &pctx);
      fb = tok_plate_psi_func(b, &pctx);
      res = gkyl_ridders(tok_plate_psi_func, &pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smin = res.res;
      geo->plate_func_lower(smin, rzplate);
      arc_ctx.zmin_right = rzplate[1];
    }
    else{
      arc_ctx.zmin_left = inp->zmin_left;
      arc_ctx.zmin_right = inp->zmin_right;
    }

    // Set arc length
    arc_ctx.rclose = inp->rright;
    arc_ctx.right = true;
    double arcL_r = integrate_psi_contour_memo(geo, geo->psisep, arc_ctx.zmin_right, arc_ctx.zmax, arc_ctx.rright,
      false, false, arc_memo_right);

    // Immediately set rclose
    arc_ctx.rclose = inp->rleft;
    arc_ctx.right = false;
    double arcL_l = integrate_psi_contour_memo(geo, geo->psisep, arc_ctx.zmin_left, arc_ctx.zmax, arc_ctx.rleft,
      false, false, arc_memo_left);
    double arcL_tot = arcL_l + arcL_r;

    if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_R) {
      *theta_lo = -M_PI+del;
      *theta_up = -M_PI+del + arcL_r/arcL_tot*2.0*M_PI;
    }
    else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_L) {
      *theta_lo = M_PI-del - arcL_l/arcL_tot*2.0*M_PI;
      *theta_up = M_PI-del;
    }
  }

  else if(inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_UP_L || inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_UP_R){
    arc_ctx.rright = inp->rright;
    arc_ctx.rleft = inp->rleft;
    //Find the lower turning point to set zmin
    double zxpt_up = geo->use_cubics ? geo->efit->Zxpt_cubic[1] : geo->efit->Zxpt[1];
    arc_ctx.zmin = zxpt_up; // Initial guess
    double zup = fmax(inp->zmax_left,  inp->zmax_right);
    find_lower_turning_point(geo, geo->psisep, zup, &arc_ctx.zmin, 1e-15);
    // Done finding turning point

    // Set zmax left and zmax right wither with plate or fixed
    // This one can't be used with the general func for setting upper and lower plates because it uses zmax left and zmax right
    if (geo->plate_spec){
      double rzplate[2];
      pctx.psi_curr = geo->psisep;
      pctx.lower=false;
      double a = 0;
      double b = 1;
      double fa = tok_plate_psi_func(a, &pctx);
      double fb = tok_plate_psi_func(b, &pctx);
      struct gkyl_qr_res res = gkyl_ridders(tok_plate_psi_func, &pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smax = res.res;
      geo->plate_func_upper(smax, rzplate);
      arc_ctx.zmax_right= rzplate[1];

      pctx.lower=true;
      a = 0;
      b = 1;
      fa = tok_plate_psi_func(a, &pctx);
      fb = tok_plate_psi_func(b, &pctx);
      res = gkyl_ridders(tok_plate_psi_func, &pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smin = res.res;
      geo->plate_func_lower(smin, rzplate);
      arc_ctx.zmax_left= rzplate[1];
    }
    else{
      arc_ctx.zmax_left = inp->zmax_left;
      arc_ctx.zmax_right = inp->zmax_right;
    }

    // Immediately set rclose
    arc_ctx.rclose = inp->rleft;
    arc_ctx.right = false;
    double arcL_l = integrate_psi_contour_memo(geo, geo->psisep, arc_ctx.zmin, arc_ctx.zmax_left, arc_ctx.rleft,
      false, false, arc_memo_left);

    // Immediately set rclose
    arc_ctx.rclose = inp->rright;
    arc_ctx.right = true;
    double arcL_r = integrate_psi_contour_memo(geo, geo->psisep, arc_ctx.zmin, arc_ctx.zmax_right, arc_ctx.rright,
      false, false, arc_memo_right);
    double arcL_tot = arcL_r + arcL_l;

    if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_UP_L) {
      *theta_lo = -M_PI+del;
      *theta_up = -M_PI+del + arcL_l/arcL_tot*2.0*M_PI;
    }
    else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_UP_R) {
      *theta_lo = M_PI-del - arcL_r/arcL_tot*2.0*M_PI;
      *theta_up = M_PI-del;
    }

  }

  gkyl_free(arc_memo);
  gkyl_free(arc_memo_left);
  gkyl_free(arc_memo_right);

}

void
tok_find_endpoints(struct gkyl_tok_geo_grid_inp* inp, struct gkyl_tok_geo *geo, struct arc_length_ctx* arc_ctx, struct plate_ctx* pctx, double psi_curr, double alpha_curr, double* arc_memo, double* arc_memo_left, double* arc_memo_right){
  enum { PH_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates


  // Set psicurr no matter what
  arc_ctx->psi = psi_curr;

  if(inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE || inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R || inp->ftype ==  GKYL_GEOMETRY_TOKAMAK_CORE_L){
    // Immediately set rleft and rright. Will need both
    arc_ctx->rright = inp->rright;
    arc_ctx->rleft = inp->rleft;

    double zxpt_up = geo->use_cubics ? geo->efit->Zxpt_cubic[1] : geo->efit->Zxpt[1];
    double zxpt_lo = geo->use_cubics ? geo->efit->Zxpt_cubic[0] : geo->efit->Zxpt[0];
    arc_ctx->zmax = inp->zmax ? inp->zmax : zxpt_up; // Initial guess.
                                                  // zmax is specified for single null full core
    double zlo = geo->zmaxis;
    find_upper_turning_point(geo, psi_curr, zlo, &arc_ctx->zmax, 0);
    arc_ctx->zmin = zxpt_lo; // Initial guess
    double zup = geo->zmaxis;
    find_lower_turning_point(geo, psi_curr, zup, &arc_ctx->zmin, 0);
    // Done finding turning points
    arc_ctx->arcL_right = integrate_psi_contour_memo(geo, psi_curr, arc_ctx->zmin, arc_ctx->zmax, arc_ctx->rright,
      true, true, arc_memo_right);
    arc_ctx->right = false;
    arc_ctx->arcL_left = integrate_psi_contour_memo(geo, psi_curr, arc_ctx->zmin, arc_ctx->zmax, arc_ctx->rleft,
      true, true, arc_memo_left);
    arc_ctx->arcL_tot = arc_ctx->arcL_left + arc_ctx->arcL_right;


    // Adjust the starting point (theta=0) so that the core blocks
    // are up-down symmetric
    if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R) {
      double theta_extent = inp->half_domain ?  2.0*(inp->cgrid.upper[TH_IDX] - inp->cgrid.lower[TH_IDX]) : inp->cgrid.upper[TH_IDX] - inp->cgrid.lower[TH_IDX];
      double arcL_extent = theta_extent/(2.0*M_PI)*arc_ctx->arcL_tot;
      double extra_arcL = arcL_extent - arc_ctx->arcL_right;
      arc_ctx->arcL_start = extra_arcL/2.0;
    }
    else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_L){
      double theta_extent = inp->half_domain ? 2.0*M_PI - 2.0*(inp->cgrid.upper[TH_IDX] - inp->cgrid.lower[TH_IDX]) : 2.0*M_PI - (inp->cgrid.upper[TH_IDX] - inp->cgrid.lower[TH_IDX]);
      double arcL_extent = theta_extent/(2.0*M_PI)*arc_ctx->arcL_tot;
      double extra_arcL = arcL_extent - arc_ctx->arcL_right;
      arc_ctx->arcL_start = extra_arcL/2.0;
    }



    arc_ctx->right = true;
    arc_ctx->phi_right = 0.0;
    arc_ctx->rclose = arc_ctx->rright;
    arc_ctx->phi_right = phi_func(alpha_curr, arc_ctx->zmax, arc_ctx) - alpha_curr;

    if(inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_L) {
      arc_ctx->right = false;
      arc_ctx->rclose = inp->rleft;
    }

  }

  else if(inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_R || inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_L){
    arc_ctx->rright = inp->rright;
    arc_ctx->rleft = inp->rleft;

    //Find the  upper turning point to set zmax
    double zxpt_lo = geo->use_cubics ? geo->efit->Zxpt_cubic[0] : geo->efit->Zxpt[0];
    arc_ctx->zmax = zxpt_lo; // Initial guess
    double zlo = fmin(inp->zmin_left, inp->zmin_right);
    find_upper_turning_point(geo, psi_curr, zlo, &arc_ctx->zmax, 1e-15);

    // Set zmin left and zmin right wither with plate or fixed
    // This one can't be used with the general func for setting upper and lower plates because it uses zmin left and zmin right
    if (geo->plate_spec){
      double rzplate[2];
      pctx->psi_curr = psi_curr;
      pctx->lower=false;
      double a = 0;
      double b = 1;
      double fa = tok_plate_psi_func(a, pctx);
      double fb = tok_plate_psi_func(b, pctx);
      struct gkyl_qr_res res = gkyl_ridders(tok_plate_psi_func, pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smax = res.res;
      geo->plate_func_upper(smax, rzplate);
      arc_ctx->zmin_left = rzplate[1];

      pctx->lower=true;
      a = 0;
      b = 1;
      fa = tok_plate_psi_func(a, pctx);
      fb = tok_plate_psi_func(b, pctx);
      res = gkyl_ridders(tok_plate_psi_func, pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smin = res.res;
      geo->plate_func_lower(smin, rzplate);
      arc_ctx->zmin_right = rzplate[1];
    }
    else{
      arc_ctx->zmin_left = inp->zmin_left;
      arc_ctx->zmin_right = inp->zmin_right;
    }

    // Immediately set rclose
    arc_ctx->rclose = inp->rright;
    arc_ctx->right = true;
    // Set arc length
    arc_ctx->arcL_right = integrate_psi_contour_memo(geo, psi_curr, arc_ctx->zmin_right, arc_ctx->zmax, arc_ctx->rright,
      true, true, arc_memo_right);
    // Immediately set rclose
    arc_ctx->right = false;
    arc_ctx->rclose = inp->rleft;
    double arcL_l = integrate_psi_contour_memo(geo, psi_curr, arc_ctx->zmin_left, arc_ctx->zmax, arc_ctx->rleft,
      true, true, arc_memo_left);
    arc_ctx->arcL_tot = arcL_l + arc_ctx->arcL_right;

    if(inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_R) {
      arc_ctx->right = true;
      arc_ctx->rclose = inp->rright;
    }
    else if(inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_L) {
      arc_ctx->right = false;
      arc_ctx->rclose = inp->rleft;
    }

  }

  else if(inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_UP_L || inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_UP_R){
    arc_ctx->rright = inp->rright;
    arc_ctx->rleft = inp->rleft;
    //Find the lower turning point to set zmin
    double zxpt_up = geo->use_cubics ? geo->efit->Zxpt_cubic[1] : geo->efit->Zxpt[1];
    arc_ctx->zmin = zxpt_up; // Initial guess
    double zup = fmax(inp->zmax_left,  inp->zmax_right);
    find_lower_turning_point(geo, psi_curr, zup, &arc_ctx->zmin, 1e-15);
    // Done finding turning point

    // Set zmax left and zmax right wither with plate or fixed
    // This one can't be used with the general func for setting upper and lower plates because it uses zmax left and zmax right
    if (geo->plate_spec){
      double rzplate[2];
      pctx->psi_curr = psi_curr;
      pctx->lower=false;
      double a = 0;
      double b = 1;
      double fa = tok_plate_psi_func(a, pctx);
      double fb = tok_plate_psi_func(b, pctx);
      struct gkyl_qr_res res = gkyl_ridders(tok_plate_psi_func, pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smax = res.res;
      geo->plate_func_upper(smax, rzplate);
      arc_ctx->zmax_right= rzplate[1];

      pctx->lower=true;
      a = 0;
      b = 1;
      fa = tok_plate_psi_func(a, pctx);
      fb = tok_plate_psi_func(b, pctx);
      res = gkyl_ridders(tok_plate_psi_func, pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smin = res.res;
      geo->plate_func_lower(smin, rzplate);
      arc_ctx->zmax_left= rzplate[1];
    }
    else{
      arc_ctx->zmax_left = inp->zmax_left;
      arc_ctx->zmax_right = inp->zmax_right;
    }


    // Immediately set rclose
    arc_ctx->rclose = inp->rleft;
    arc_ctx->right = false;
    arc_ctx->arcL_left = integrate_psi_contour_memo(geo, psi_curr, arc_ctx->zmin, arc_ctx->zmax_left, arc_ctx->rleft,
      true, true, arc_memo_left);

    // Immediately set rclose
    arc_ctx->rclose = inp->rright;
    arc_ctx->right = true;
    double arcL_r = integrate_psi_contour_memo(geo, psi_curr, arc_ctx->zmin, arc_ctx->zmax_right, arc_ctx->rright,
      true, true, arc_memo_right);
    arc_ctx->arcL_tot = arcL_r + arc_ctx->arcL_left;

    if(inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_UP_R) {
      arc_ctx->right = true;
      arc_ctx->rclose = inp->rright;
    }
    else if(inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_UP_L) {
      arc_ctx->right = false;
      arc_ctx->rclose = inp->rleft;
    }
  }

  else if(inp->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT || inp->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO || inp->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID || inp->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP){
    // Immediately set rclose
    arc_ctx->rclose = inp->rright;
    // Set zmin and zmax either fixed or with plate
    if (geo->plate_spec){
      set_upper_plate(geo, arc_ctx, pctx, arc_ctx->psi);
      set_lower_plate(geo, arc_ctx, pctx, arc_ctx->psi);
    }
    else{
      arc_ctx->zmin = inp->zmin;
      arc_ctx->zmax = inp->zmax;
    }
    // Set the arc length
    arc_ctx->arcL_tot = integrate_psi_contour_memo(geo, psi_curr, arc_ctx->zmin, arc_ctx->zmax, arc_ctx->rclose, true, true, arc_memo);
  }

  else if(inp->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN || inp->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO || inp->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID || inp->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_UP){
    // Immediately set rclose
    arc_ctx->rclose = inp->rleft;
    // Set zmin and zmax either fixed or with plate
    if (geo->plate_spec){
      set_upper_plate(geo, arc_ctx, pctx, arc_ctx->psi);
      set_lower_plate(geo, arc_ctx, pctx, arc_ctx->psi);
    }
    else{
      arc_ctx->zmin = inp->zmin;
      arc_ctx->zmax = inp->zmax;
    }
    // Set the arc Length
    arc_ctx->arcL_tot = integrate_psi_contour_memo(geo, psi_curr, arc_ctx->zmin, arc_ctx->zmax, arc_ctx->rclose, true, true, arc_memo);
  }

  else if(inp->ftype==GKYL_GEOMETRY_TOKAMAK_LSN_SOL || inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO || inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID || inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_UP){
    // Immediately set rleft and rright. Will need both
    arc_ctx->rright = inp->rright;
    arc_ctx->rleft = inp->rleft;
    //Find the  upper turning point
    arc_ctx->zmax = inp->zmax; // Initial guess
    double zlo = fmax(inp->zmin_left, inp->zmin_right);
    find_upper_turning_point(geo, psi_curr, zlo, &arc_ctx->zmax, 0);

    // Set zmin left and zmin right wither with plate or fixed
    // This one can't be used with the general func for setting upper and lower plates because it uses zmin left and zmin right
    if (geo->plate_spec){
      double rzplate[2];
      pctx->psi_curr = psi_curr;
      pctx->lower=false;
      double a = 0;
      double b = 1;
      double fa = tok_plate_psi_func(a, pctx);
      double fb = tok_plate_psi_func(b, pctx);
      struct gkyl_qr_res res = gkyl_ridders(tok_plate_psi_func, pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smax = res.res;
      geo->plate_func_upper(smax, rzplate);
      arc_ctx->zmin_left = rzplate[1];

      pctx->lower=true;
      a = 0;
      b = 1;
      fa = tok_plate_psi_func(a, pctx);
      fb = tok_plate_psi_func(b, pctx);
      res = gkyl_ridders(tok_plate_psi_func, pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smin = res.res;
      geo->plate_func_lower(smin, rzplate);
      arc_ctx->zmin_right = rzplate[1];
    }
    else{
      arc_ctx->zmin_left = inp->zmin_left;
      arc_ctx->zmin_right = inp->zmin_right;
    }

    // Done finding turning point
    arc_ctx->arcL_right = integrate_psi_contour_memo(geo, psi_curr, arc_ctx->zmin_right, arc_ctx->zmax, arc_ctx->rright,
      true, true, arc_memo_right);
    arc_ctx->right = false;
    double arcL_l = integrate_psi_contour_memo(geo, psi_curr, arc_ctx->zmin_left, arc_ctx->zmax, arc_ctx->rleft,
      true, true, arc_memo_left);
    arc_ctx->arcL_tot = arcL_l + arc_ctx->arcL_right;

    arc_ctx->right = true;
    arc_ctx->phi_right = 0.0;
    arc_ctx->rclose = arc_ctx->rright;
    arc_ctx->psi = psi_curr;
    arc_ctx->zmin = arc_ctx->zmin_right;
    arc_ctx->phi_right = phi_func(alpha_curr, arc_ctx->zmax, arc_ctx) - alpha_curr;
  }

  else if(inp->ftype == GKYL_GEOMETRY_TOKAMAK_IWL){
    // Immediately set rleft and rright. Will need both
    arc_ctx->rright = inp->rright;
    arc_ctx->rleft = inp->rleft;

    double R_lcfs[4], dRdZ_lcfs[4];
    double dR_lcfs[4], dZ_lcfs[4];
    int nr_lcfs = gkyl_tok_geo_R_psiZ(geo, geo->efit->sibry, geo->efit->zmaxis, 4, R_lcfs, dRdZ_lcfs, dR_lcfs, dZ_lcfs);
    double r_lcfs = nr_lcfs == 1 ? R_lcfs[0] : choose_closest(arc_ctx->rleft, R_lcfs, R_lcfs, nr_lcfs);
    double rz_lcfs[2];

    double rmin_old = geo->rmin;
    if (geo->plate_spec) {
      geo->plate_func_upper(0.0, rz_lcfs);
      if(fabs(rz_lcfs[0] - r_lcfs) > 1e-6) {
        fprintf(stderr, "The upper plate function has an error. It must return (R(s=0),Z(s=0)) = (%1.16f, %1.16f). \n", R_lcfs[0], geo->efit->zmaxis);
        assert(false);
      }
      geo->plate_func_lower(0.0, rz_lcfs);
      if(fabs(rz_lcfs[0] - r_lcfs) > 1e-6) {
        fprintf(stderr, "The lower plate function has an error. It must return (R(s=0),Z(s=0)) = (%1.16f, %1.16f). \n", R_lcfs[0], geo->efit->zmaxis);
        assert(false);
      }
  
      if ( (arc_ctx->psi >= geo->efit->sibry && geo->efit->sibry >= geo->efit->simag) ||
           (arc_ctx->psi <= geo->efit->sibry && geo->efit->sibry <= geo->efit->simag) ) {
        set_upper_iwl_plate(geo, arc_ctx, pctx, arc_ctx->psi);
        set_lower_iwl_plate(geo, arc_ctx, pctx, arc_ctx->psi);
      }
    }
    else {
      arc_ctx->zmin_iwl_plate = geo->zmaxis;
      arc_ctx->zmax_iwl_plate = geo->zmaxis;
    }

    arc_ctx->zmax = inp->zmax; // Initial guess.

    double zlo = geo->zmaxis;
    find_upper_turning_point(geo, psi_curr, zlo, &arc_ctx->zmax, 0);
    arc_ctx->zmin = inp->zmin; // Initial guess
    double zup = geo->zmaxis;
    find_lower_turning_point(geo, psi_curr, zup, &arc_ctx->zmin, 0);
    // Done finding turning points
    arc_ctx->zmin_iwl = arc_ctx->zmin;
    arc_ctx->zmax_iwl = arc_ctx->zmax;
    geo->rmin = rmin_old;

    arc_ctx->right = true;
    arc_ctx->arcL_q1 = integrate_psi_contour_memo(geo, psi_curr, geo->zmaxis, arc_ctx->zmax, arc_ctx->rright, false, false, arc_memo);
    arc_ctx->right = false;
    arc_ctx->arcL_q2 = integrate_psi_contour_memo(geo, psi_curr, arc_ctx->zmin_iwl_plate, arc_ctx->zmax, arc_ctx->rleft, false, false, arc_memo);
    arc_ctx->right = false;
    arc_ctx->arcL_q3 = integrate_psi_contour_memo(geo, psi_curr, arc_ctx->zmin, arc_ctx->zmax_iwl_plate, arc_ctx->rleft, false, false, arc_memo);
    arc_ctx->right = true;
    arc_ctx->arcL_q4 = integrate_psi_contour_memo(geo, psi_curr, arc_ctx->zmin, geo->zmaxis, arc_ctx->rright, false, false, arc_memo);

    arc_ctx->arcL_tot = arc_ctx->arcL_q1 + arc_ctx->arcL_q2 + arc_ctx->arcL_q3 + arc_ctx->arcL_q4;

    arc_ctx->right = true;
    arc_ctx->rclose = arc_ctx->rright;
    arc_ctx->phi_right = phi_func(alpha_curr, arc_ctx->zmax, arc_ctx) - alpha_curr;
  }

}



void
tok_set_ridders(struct gkyl_tok_geo_grid_inp* inp, struct arc_length_ctx* arc_ctx, double psi_curr, double arcL_curr,double* rclose, double *ridders_min, double* ridders_max){


  if(inp->ftype==GKYL_GEOMETRY_TOKAMAK_CORE){
    if(arcL_curr <= arc_ctx->arcL_right){
      *rclose = arc_ctx->rright;
      arc_ctx->right = true;
      *ridders_min = -arcL_curr;
      *ridders_max = arc_ctx->arcL_tot-arcL_curr;
    }
    else{
      *rclose = arc_ctx->rleft;
      arc_ctx->right = false;
      *ridders_min = arc_ctx->arcL_tot - arcL_curr;
      *ridders_max = -arcL_curr + arc_ctx->arcL_right;
    }
  }

  if(inp->ftype==GKYL_GEOMETRY_TOKAMAK_CORE_R || inp->ftype==GKYL_GEOMETRY_TOKAMAK_CORE_L){
    if(arcL_curr <= arc_ctx->arcL_start){
      *rclose = arc_ctx->rleft;
      arc_ctx->right = false;
      arc_ctx->pre= true;
      *ridders_min = arc_ctx->arcL_start-arcL_curr;
      *ridders_max = arc_ctx->arcL_start - arc_ctx->arcL_left - arcL_curr;
    }
    else if(arcL_curr <= arc_ctx->arcL_right + arc_ctx->arcL_start){
      *rclose = arc_ctx->rright;
      arc_ctx->right = true;
      arc_ctx->pre= false;
      *ridders_min = arc_ctx->arcL_start -arcL_curr;
      *ridders_max = arc_ctx->arcL_right + arc_ctx->arcL_start - arcL_curr;
    }
    else{
      *rclose = arc_ctx->rleft;
      arc_ctx->right = false;
      arc_ctx->pre= false;
      *ridders_min = arc_ctx->arcL_tot + arc_ctx->arcL_start - arcL_curr;
      *ridders_max = -arcL_curr + arc_ctx->arcL_right+arc_ctx->arcL_start;
    }
  }


  else if(inp->ftype==GKYL_GEOMETRY_TOKAMAK_PF_LO_R || inp->ftype==GKYL_GEOMETRY_TOKAMAK_PF_LO_L){
    if(arcL_curr <= arc_ctx->arcL_right){
      *rclose = arc_ctx->rright;
      arc_ctx->right = true;
      *ridders_min = -arcL_curr;
      *ridders_max = arc_ctx->arcL_tot-arcL_curr;
      arc_ctx->zmin = arc_ctx->zmin_right;
    }
    else{
      *rclose = arc_ctx->rleft;
      arc_ctx->right = false;
      *ridders_min = arc_ctx->arcL_tot - arcL_curr;
      *ridders_max = -arcL_curr + arc_ctx->arcL_right;
      arc_ctx->zmin = arc_ctx->zmin_left;
    }
  }

  else if(inp->ftype==GKYL_GEOMETRY_TOKAMAK_PF_UP_R || inp->ftype==GKYL_GEOMETRY_TOKAMAK_PF_UP_L){
    if(arcL_curr <= arc_ctx->arcL_left){
      *rclose = arc_ctx->rleft;
      arc_ctx->right = false;
      *ridders_min = -arcL_curr + arc_ctx->arcL_left;
      *ridders_max = -arcL_curr;
      arc_ctx->zmax = arc_ctx->zmax_left;
    }
    else{
      *rclose = arc_ctx->rright;
      arc_ctx->right = true;
      *ridders_min = arc_ctx->arcL_left - arcL_curr;
      *ridders_max = -arcL_curr + arc_ctx->arcL_tot;
      arc_ctx->zmax = arc_ctx->zmax_right;
    }
  }

  else if( (arc_ctx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT) || (arc_ctx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT) || (arc_ctx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO) || (arc_ctx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID) || (arc_ctx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP) ){
    *ridders_min = -arcL_curr;
    *ridders_max = arc_ctx->arcL_tot-arcL_curr;
    *rclose = arc_ctx->rclose;
  }
  else if( (arc_ctx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN) || (arc_ctx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN) || (arc_ctx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO) || (arc_ctx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID) || (arc_ctx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_UP) ){
    *ridders_min = arc_ctx->arcL_tot-arcL_curr;
    *ridders_max = -arcL_curr;
    *rclose = arc_ctx->rclose;
  }
  else if(arc_ctx->ftype==GKYL_GEOMETRY_TOKAMAK_LSN_SOL || arc_ctx->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO || arc_ctx->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID || arc_ctx->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_UP){
    if(arcL_curr <= arc_ctx->arcL_right){
      *rclose = arc_ctx->rright;
      arc_ctx->right = true;
      *ridders_min = -arcL_curr;
      *ridders_max = arc_ctx->arcL_tot-arcL_curr;
      arc_ctx->zmin = arc_ctx->zmin_right;
    }
    else{
      *rclose = arc_ctx->rleft;
      arc_ctx->right = false;
      *ridders_min = arc_ctx->arcL_tot - arcL_curr;
      *ridders_max = -arcL_curr + arc_ctx->arcL_right;
      arc_ctx->zmin = arc_ctx->zmin_left;
    }
  }

  else if (arc_ctx->ftype==GKYL_GEOMETRY_TOKAMAK_IWL) {
    arc_ctx->q1 = false;
    arc_ctx->q2 = false;
    arc_ctx->q3 = false;
    arc_ctx->q4 = false;
    if(arcL_curr <= arc_ctx->arcL_q3) {
      arc_ctx->q3 = true;
      *rclose = arc_ctx->rleft;
      arc_ctx->right = false;
      *ridders_min = arc_ctx->arcL_q3 - arcL_curr;
      *ridders_max = -arcL_curr;
      arc_ctx->zmin = arc_ctx->zmin_iwl;
      arc_ctx->zmax = arc_ctx->zmax_iwl_plate;
    }

    else if(arcL_curr <= ( arc_ctx->arcL_q3 + arc_ctx->arcL_q4) ) {
      arc_ctx->q4 = true;
      *rclose = arc_ctx->rright;
      arc_ctx->right = true;
      *ridders_min = arc_ctx->arcL_q3 - arcL_curr;
      *ridders_max = arc_ctx->arcL_q3 + arc_ctx->arcL_q4 - arcL_curr;
      arc_ctx->zmin = arc_ctx->zmin_iwl;
      arc_ctx->zmax = arc_ctx->geo->zmaxis;
    }

    else if(arcL_curr <= (arc_ctx->arcL_q3 + arc_ctx->arcL_q4 + arc_ctx->arcL_q1) ) {
      arc_ctx->q1 = true;
      *rclose = arc_ctx->rright;
      arc_ctx->right = true;
      *ridders_min = arc_ctx->arcL_q3 + arc_ctx->arcL_q4 - arcL_curr;
      *ridders_max = arc_ctx->arcL_q3 + arc_ctx->arcL_q4 + arc_ctx->arcL_q1 - arcL_curr;
      arc_ctx->zmin = arc_ctx->geo->zmaxis;
      arc_ctx->zmax = arc_ctx->zmax_iwl;
    }

    else {
      arc_ctx->q2 = true;
      *rclose = arc_ctx->rleft;
      arc_ctx->right = false;
      *ridders_min = arc_ctx->arcL_q3 + arc_ctx->arcL_q4 + arc_ctx->arcL_q1 + arc_ctx->arcL_q2 - arcL_curr;
      *ridders_max = arc_ctx->arcL_q3 + arc_ctx->arcL_q4 + arc_ctx->arcL_q1 - arcL_curr;
      arc_ctx->zmin = arc_ctx->zmin_iwl_plate;
      arc_ctx->zmax = arc_ctx->zmax_iwl;
    }
  }

  arc_ctx->arcL = arcL_curr;
  arc_ctx->rclose = *rclose; // This would be unnecessary for all double null block cases. Only needed for SN and full core
}
