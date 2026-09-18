#pragma once

// Experimental C1 arc map. The topology plan supplies physical endpoint
// speeds. Match both these derivatives and the actual endpoint CELL chords.
// All interpolation coefficients are derived from row lengths and cell count.

// Legacy blocks do not have a stored trace of the row: their nodes come from
// root-finding Z at a requested arc length and then taking the psi-contour root
// nearest rclose.  Ask the SAME mapper for the chord samples.  Building a
// separate domain trace instead is both a different parameterization of the row
// and, for rows whose topology endpoints the extended construction cannot place
// (measured on the ASDEX LSN_SOL_LO block off the separatrix), not buildable at
// all.  tok_set_ridders rewrites the bracket fields it is handed, so the probe
// runs on a copy and leaves the live context to the caller that follows it.
static bool
tok_row_arc_native_point(const struct gkyl_tok_geo_grid_inp *inp,
  const struct arc_length_ctx *ctx, double fraction, double *r, double *z)
{
  struct gkyl_tok_geo_grid_inp probe_inp = *inp;
  struct arc_length_ctx probe = *ctx;
  const double lo = inp->cgrid.lower[2], w = inp->cgrid.upper[2]-lo;
  if (!(w > 0.0) || !isfinite(fraction)) return false;
  double arcL = tok_xpt_theta_to_arc(inp, &probe, lo+fraction*w);
  double rclose = probe.rclose, rid_lo = 0.0, rid_hi = 0.0;
  tok_set_ridders(&probe_inp, &probe, probe.psi, arcL, &rclose, &rid_lo, &rid_hi);
  struct gkyl_qr_res res = gkyl_ridders(arc_length_func, &probe,
    probe.zmin, probe.zmax, rid_lo, rid_hi, probe.geo->root_param.max_iter, 1e-10);
  // Same test tok_geo_check_arc_root applies to the construction itself: a
  // usable root lies inside the Z interval it was sought in, and the sentinel
  // the search starts from does not.  Reported, not aborted, so the caller can
  // name the map it failed in.
  if (!(res.res >= probe.zmin && res.res <= probe.zmax)) return false;
  double R[4] = { 0 }, dRdZ[4] = { 0 }, dR[4] = { 0 }, dZ[4] = { 0 };
  int nr = gkyl_tok_geo_R_psiZ(probe.geo, probe.psi, res.res, 4, R, dRdZ, dR, dZ);
  if (nr == 0) return false;
  *r = choose_closest(rclose, R, R, nr);
  *z = res.res;
  return isfinite(*r) && isfinite(*z);
}

// Arclength fraction along the row -> the parameter the sampler wants. Split
// out of tok_row_arc_point so a diagnostic can see the two halves separately:
// a jump in the sampled point is either this conversion or the consumer, and
// they need different fixes.
//
// NOTE this is a SECOND copy of the same branch that lives in
// tok_ordered_map_lookup. The separatrix-ruler fix patched that one; this one
// still skips the conversion on the separatrix row.
static bool
tok_row_arc_fraction_to_u(const struct gkyl_tok_geo_grid_inp *inp,
  const struct arc_length_ctx *ctx, double fraction, double *u)
{
  *u=fraction;
  if (tok_geo_same_flux(ctx->psi,ctx->geo->psisep) && tok_row_arc_sep_skip())
    return true;
  int lo=0, hi=ctx->map_trace_n-1;
  double target=fraction*ctx->map_trace_s[hi];
  while (hi-lo > 1) {
    int mid=lo+(hi-lo)/2;
    if (ctx->map_trace_s[mid] <= target) lo=mid; else hi=mid;
  }
  double ds=ctx->map_trace_s[lo+1]-ctx->map_trace_s[lo];
  if (!(ds > 0.0)) return false;
  *u=(lo+(target-ctx->map_trace_s[lo])/ds)/(ctx->map_trace_n-1);
  return true;
}

static bool
tok_row_arc_point(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *ctx, double fraction, double *r, double *z)
{
  if (!tok_xpt_ordered_placement(inp))
    return tok_row_arc_native_point(inp,ctx,fraction,r,z);
  double u;
  if (!tok_row_arc_fraction_to_u(inp,ctx,fraction,&u)) return false;
  if (tok_ext_construction(inp))
    return tok_logical_trace_sample(ctx->geo,ctx->psi,ctx->map_trace_r,
      ctx->map_trace_z,ctx->map_trace_n,u,r,z);
  return tok_parameterized_xpt_seam_point(inp,ctx,ctx->psi,u,r,z);
}

// Resolution of the coarse scan that brackets the end-cell chord.  The ordered
// path walks its stored trace; the legacy path has none, so it scans the row at
// the resolution those traces are built with.
static int
tok_row_arc_scan_nodes(const struct gkyl_tok_geo_grid_inp *inp,
  const struct arc_length_ctx *ctx)
{
  if (tok_xpt_ordered_placement(inp)) return ctx->map_trace_n;
  return GKYL_MIN2(257, ctx->sep_trace_capacity);
}

static bool tok_row_arc_plan_dump(void);

// Every refusal names itself. These used to return false mute, and a caller
// three frames up aborted with no reason at all.
static bool
tok_row_arc_end_distance(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *ctx, int end, double chord, double length, double *distance)
{
  double r0,z0;
  if (!(chord > 0.0)) {
    fprintf(stderr,"TOK_ROW_ARC_FAILED reason=end_chord_nonpositive ftype=%d psi=%.17g end=%d chord=%.17g\n",
      inp->ftype,ctx->psi,end,chord);
    return false;
  }
  if (!tok_row_arc_point(inp,ctx,end,&r0,&z0)) {
    fprintf(stderr,"TOK_ROW_ARC_FAILED reason=end_point_unmappable ftype=%d psi=%.17g end=%d\n",
      inp->ftype,ctx->psi,end);
    return false;
  }
  const int n=tok_row_arc_scan_nodes(inp,ctx);
  if (n < 2) {
    fprintf(stderr,"TOK_ROW_ARC_FAILED reason=no_scan_resolution ftype=%d psi=%.17g nodes=%d\n",
      inp->ftype,ctx->psi,n);
    return false;
  }
  // Accumulate the REAL arclength between successive samples. The parameter is
  // not proportional to arclength on every path -- the chord construction
  // grades its own logical u, and the trace sampler walks index space -- so
  // scaling the parameter by the row length measured one row at 4.74x its own
  // end-cell arc. Summing chords is correct for all three parameterizations,
  // which also makes the answer independent of which one sampled the points.
  double lo=0.0, hi=0.0;
  double arc_lo=0.0, arc=0.0;
  double pr=r0, pz=z0;
  bool bracket=false;
  for (int i=1; i<n; ++i) {
    hi=i/(double)(n-1);
    double r,z;
    if (!tok_row_arc_point(inp,ctx,end ? 1.0-hi:hi,&r,&z)) {
      fprintf(stderr,"TOK_ROW_ARC_FAILED reason=scan_point_unmappable ftype=%d psi=%.17g end=%d fraction=%.17g\n",
        inp->ftype,ctx->psi,end,end ? 1.0-hi:hi);
      return false;
    }
    arc += hypot(r-pr,z-pz);
    pr=r; pz=z;
    if (hypot(r-r0,z-z0) >= chord) { bracket=true; break; }
    lo=hi; arc_lo=arc;
  }
  if (!bracket) {
    // The whole row is shorter than the cell chord the plan asked this end to
    // span, so no interior point is far enough away to bracket it.
    fprintf(stderr,"TOK_ROW_ARC_FAILED reason=chord_exceeds_row ftype=%d psi=%.17g end=%d chord=%.17g length=%.17g\n",
      inp->ftype,ctx->psi,end,chord,length);
    return false;
  }
  // The bisection moves `lo`, but `arc_lo` was accumulated at the scan node that
  // opened the bracket. Keep that node: the refinement below has to integrate
  // from it, not from wherever the bisection leaves `lo`.
  const double bracket_lo=lo;
  for (int k=0; k<=DBL_MANT_DIG; ++k) {
    double mid=0.5*(lo+hi), r,z;
    if (mid == lo || mid == hi) break;
    if (!tok_row_arc_point(inp,ctx,end ? 1.0-mid:mid,&r,&z)) {
      fprintf(stderr,"TOK_ROW_ARC_FAILED reason=bisect_point_unmappable ftype=%d psi=%.17g end=%d fraction=%.17g\n",
        inp->ftype,ctx->psi,end,end ? 1.0-mid:mid);
      return false;
    }
    if (hypot(r-r0,z-z0) < chord) lo=mid; else hi=mid;
  }
  // Integrate the last sub-interval at a step finer than the scan's, and add it
  // to the arclength already accumulated at the bracketing scan node. Anchoring
  // on `lo` here dropped the whole sub-interval -- the bisection had narrowed it
  // to one ulp -- so the answer was the arc at the previous SCAN node, short by
  // up to one scan step. That is why a well-sampled row came back a few 1e-5
  // under its own chord and a row whose trace is 10x non-uniform at its end came
  // back 16% under.
  const double target=0.5*(lo+hi);
  const int nsub=8;
  double pr2=0.0, pz2=0.0;
  if (!tok_row_arc_point(inp,ctx,end ? 1.0-bracket_lo:bracket_lo,&pr2,&pz2)) {
    fprintf(stderr,"TOK_ROW_ARC_FAILED reason=refine_anchor_unmappable ftype=%d psi=%.17g end=%d\n",
      inp->ftype,ctx->psi,end);
    return false;
  }
  *distance=arc_lo;
  for (int k=1; k<=nsub; ++k) {
    double f=bracket_lo+(target-bracket_lo)*k/(double) nsub, r,z;
    if (!tok_row_arc_point(inp,ctx,end ? 1.0-f:f,&r,&z)) {
      fprintf(stderr,"TOK_ROW_ARC_FAILED reason=refine_point_unmappable ftype=%d psi=%.17g end=%d fraction=%.17g\n",
        inp->ftype,ctx->psi,end,end ? 1.0-f:f);
      return false;
    }
    *distance += hypot(r-pr2,z-pz2);
    pr2=r; pz2=z;
  }
  // An arc can never be shorter than the chord it subtends. If it is, the row
  // was sampled inconsistently and the map must not be built on it.
  //
  // Subtended by the point the arc was actually measured to, not by the chord
  // that was requested: the bisection stops when the parameter interval closes,
  // so the two differ by its residual, and a guard with no tolerance turns that
  // residual into a refusal. This comparison needs no tolerance -- both sides
  // come from the same sample, so the triangle inequality holds exactly.
  const double chord_at_target=hypot(pr2-r0,pz2-z0);
  if (*distance < chord_at_target) {
    fprintf(stderr,"TOK_ROW_ARC_FAILED reason=arc_below_chord ftype=%d psi=%.17g end=%d arc=%.17g chord=%.17g requested=%.17g\n",
      inp->ftype,ctx->psi,end,*distance,chord_at_target,chord);
    tok_chord_stage_dump=false;
    return false;
  }
  if (*distance >= length) {
    fprintf(stderr,"TOK_ROW_ARC_FAILED reason=end_cell_spans_row ftype=%d psi=%.17g end=%d arc=%.17g length=%.17g\n",
      inp->ftype,ctx->psi,end,*distance,length);
    return false;
  }
  if (!isfinite(*distance) || !(*distance > 0.0)) {
    fprintf(stderr,"TOK_ROW_ARC_FAILED reason=end_distance_nonpositive ftype=%d psi=%.17g end=%d distance=%.17g "
      "chord=%.17g bracket_lo=%.17g target=%.17g arc_lo=%.17g scan_n=%d\n",
      inp->ftype,ctx->psi,end,*distance,chord,bracket_lo,target,arc_lo,n);
    // Walk out from the end so the sampler's behaviour near u=end is visible:
    // a distance of zero means the bisection could not find ANY point closer
    // than `chord`, which is a statement about the sampler, not the geometry.
    tok_chord_stage_dump=true;
    for (int k=0; k<=8; ++k) {
      const double f=k/(8.0*(n-1)), r=0.0, z=0.0; (void) r; (void) z;
      double rr,zz,uu=0.0;
      const double frac=end ? 1.0-f:f;
      tok_row_arc_fraction_to_u(inp,ctx,frac,&uu);
      if (!tok_row_arc_point(inp,ctx,frac,&rr,&zz)) break;
      fprintf(stderr,"  TOK_ROW_ARC_ENDWALK k=%d f=%.17g fraction=%.17g u=%.17g r=%.17g z=%.17g chord_from_end=%.17g\n",
        k,f,frac,uu,rr,zz,hypot(rr-r0,zz-z0));
    }
    return false;
  }
  // Did the bisection actually hit the chord it was asked for? The seam is
  // element = chord(first cell)/dtheta, so `chord_at_target` IS the quantity the
  // seam measures; if it differs from `chord`, the end cell is not the size the
  // plan asked for and every block misses by its own amount.
  if (tok_row_arc_plan_dump())
    fprintf(stderr,
      "TOK_ROW_ARC_END ftype=%d psi=%.17g end=%d requested_chord=%.17g "
      "realized_chord=%.17g arc=%.17g realized/requested=%.17g arc/chord=%.17g\n",
      inp->ftype, ctx->psi, end, chord, chord_at_target, *distance,
      chord_at_target/chord, *distance/chord_at_target);
  return true;
}

// One bounded dump at the point of refusal: walk the end scan again and show
// where the chord stops agreeing with the arc the trace says was travelled.
// chord/arc > 1 is impossible for a true arc, so whichever step first breaks it
// names the cause -- a smooth drift means map_trace_s (a chord sum between trace
// NODES) is under-measuring a curved row; a jump means the sampler's root
// projection changed contour branch.
static void
tok_row_arc_dump_end_scan(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *ctx, int end, double chord, double length)
{
  double r0,z0;
  if (!tok_row_arc_point(inp,ctx,end,&r0,&z0)) {
    fprintf(stderr,"TOK_ROW_ARC_SCAN ftype=%d end=%d endpoint_unmappable\n",inp->ftype,end);
    return;
  }
  const int n=tok_row_arc_scan_nodes(inp,ctx);
  fprintf(stderr,"TOK_ROW_ARC_SCAN ftype=%d psi=%.17g end=%d chord_target=%.17g length=%.17g scan_nodes=%d r0=%.17g z0=%.17g\n",
    inp->ftype,ctx->psi,end,chord,length,n,r0,z0);
  for (int i=1; i<n && i<=32; ++i) {
    double f=i/(double)(n-1), r,z;
    if (!tok_row_arc_point(inp,ctx,end ? 1.0-f:f,&r,&z)) {
      fprintf(stderr,"  i=%d fraction=%.17g UNMAPPABLE\n",i,end ? 1.0-f:f);
      return;
    }
    double c=hypot(r-r0,z-z0), a=f*length;
    fprintf(stderr,"  i=%2d arc=%.9e chord=%.9e chord/arc=%.6f r=%.9f z=%.9f\n",
      i,a,c,a > 0.0 ? c/a:0.0,r,z);
    if (c >= chord) { fprintf(stderr,"  bracket fires at i=%d\n",i); break; }
  }
  const int m=ctx->map_trace_n;
  fprintf(stderr,"  trace nodes inward from end=%d (map_trace_n=%d):\n",end,m);
  for (int k=0; k<8 && k<m; ++k) {
    int i=end ? m-1-k:k;
    fprintf(stderr,"    idx=%4d s=%.9e r=%.9f z=%.9f\n",
      i,ctx->map_trace_s[i],ctx->map_trace_r[i],ctx->map_trace_z[i]);
  }
}

// Width of the two OUTER segments of the end-cell map, as a fraction of the row
// in u.
//
// Historically exactly 1/nc -- one cell. That is why the seam stops converging:
// the metric measures the first CELL's average slope while the plan constrains
// the map's ENDPOINT slope, and with a one-cell segment the two are the same
// thing only if the cubic does not bend. It does bend, by the fixed fractional
// element jump between a block's own grading and its neighbour's (measured
// ~8% on NSTX-U b3, constant while the cell shrank 7.4x). Widen the segment and
// the first cell becomes a small part of it, so the cell average tends to the
// endpoint slope as theta refines -- which is exactly the quantity the seam
// constraint equalises.
//
// Returns 1/nc unless GKYL_TOK_ROW_ARC_BLEND_FRACTION is set, so the default is
// bit-identical to before.
static double
tok_row_arc_blend_fraction(int nc)
{
  const double cell=1.0/nc;
  const char *e=getenv("GKYL_TOK_ROW_ARC_BLEND_FRACTION");
  // NOT the on/off idiom: "0.125" starts with '0', and testing e[0]=='0'
  // silently rejected every fraction below one. Parse, then judge the value.
  if (!e || !e[0]) return cell;
  double f=atof(e);
  if (!(f > 0.0)) return cell;
  // Never narrower than a cell, and the middle segment must survive.
  return fmin(0.45, fmax(cell, f));
}

// Set GKYL_TOK_ROW_ARC_PLAN_DUMP=1 to emit, once per fitted row, the inputs the
// C1 cubic is built from. The seam is element = chord(first cell)/dtheta, and
// with blend = 1/nc the first knot IS the first node, so the first cell's chord
// is whatever end_distance was asked for -- sigma*dtheta. If the two sides of an
// interface see the same sigma, the seam must be 1. This says whether they do.
static bool
tok_row_arc_plan_dump(void)
{
  const char *e = getenv("GKYL_TOK_ROW_ARC_PLAN_DUMP");
  return e && e[0] != '\0' && e[0] != '0';
}

static bool
tok_row_arc_c1_fraction(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *ctx, double length, const double sigma[2],
  double u, double *v, double *derivative)
{
  int nc=inp->cgrid.cells[2];
  if (nc < 3) {
    fprintf(stderr,"TOK_ROW_ARC_FAILED reason=c1_prototype_needs_three_cells cells=%d\n",nc);
    return false;
  }
  double width=inp->cgrid.upper[2]-inp->cgrid.lower[2];
  if (!ctx->row_arc_fit_ready || ctx->row_arc_fit_psi != ctx->psi ||
      ctx->row_arc_fit_length != length || ctx->row_arc_fit_sigma[0] != sigma[0] ||
      ctx->row_arc_fit_sigma[1] != sigma[1]) {
    const double blend=tok_row_arc_blend_fraction(nc);
    double s0,s1;
    // The chord asked of each end is the endpoint slope carried across the
    // segment, so it scales with the segment width, not with one cell.
    if (!tok_row_arc_end_distance(inp,ctx,0,sigma[0]*width*blend,length,&s0) ||
        !tok_row_arc_end_distance(inp,ctx,1,sigma[1]*width*blend,length,&s1) ||
        !(s0+s1 < length)) {
      fprintf(stderr,"TOK_ROW_ARC_FAILED reason=end_cell_infeasible ftype=%d psi=%.17g length=%.17g sigma=(%.17g,%.17g)\n",
        inp->ftype,ctx->psi,length,sigma[0],sigma[1]);
      return false;
    }
    double *y=ctx->row_arc_fit_y, *m=ctx->row_arc_fit_m;
    y[0]=0.0; y[1]=s0; y[2]=length-s1; y[3]=length;
    const double h[3]={blend,1.0-2.0*blend,blend};
    double d[3];
    for (int i=0; i<3; ++i) d[i]=(y[i+1]-y[i])/h[i];
    m[0]=sigma[0]*width; m[3]=sigma[1]*width;
    for (int i=1; i<3; ++i) {
      double w1=2.0*h[i]+h[i-1], w2=h[i]+2.0*h[i-1];
      m[i]=(w1+w2)/(w1/d[i-1]+w2/d[i]);
    }
    for (int i=0; i<3; ++i) {
      double b=2.0*(3.0*d[i]-2.0*m[i]-m[i+1]);
      double a=3.0*(m[i]+m[i+1]-2.0*d[i]);
      double minimum=fmin(m[i],m[i+1]);
      if (a > 0.0) {
        double t=-b/(2.0*a);
        if (t > 0.0 && t < 1.0) minimum=fmin(minimum,m[i]+b*t+a*t*t);
      }
      if (!(minimum > 0.0) || !isfinite(minimum)) {
        // The piecewise cubic through the four knots is not monotone on this
        // segment, so the requested endpoint rates cannot be met by a positive
        // map through the measured end-cell arcs.
        // length is the row the ORDERED MAP samples; sigma*width is the rate
        // the PLAN measured on its own probe of that row. If those disagree the
        // cubic is being asked to span one row at another row's speed.
        fprintf(stderr,"TOK_ROW_ARC_FAILED reason=cubic_not_monotone ftype=%d psi=%.17g segment=%d minimum=%.17g m=(%.17g,%.17g) secant=%.17g length=%.17g cells=%d s0=%.17g s1=%.17g plan_rate_lo=%.17g plan_rate_hi=%.17g\n",
          inp->ftype,ctx->psi,i,minimum,m[i],m[i+1],d[i],length,nc,y[1],length-y[2],
          sigma[0]*width,sigma[1]*width);
        tok_row_arc_dump_end_scan(inp,ctx,0,sigma[0]*width*blend,length);
        tok_row_arc_dump_end_scan(inp,ctx,1,sigma[1]*width*blend,length);
        return false;
      }
    }
    ctx->row_arc_fit_psi=ctx->psi; ctx->row_arc_fit_length=length;
    ctx->row_arc_fit_sigma[0]=sigma[0]; ctx->row_arc_fit_sigma[1]=sigma[1];
    ctx->row_arc_fit_ready=true;
    if (tok_row_arc_plan_dump())
      fprintf(stderr,
        "TOK_ROW_ARC_PLAN ftype=%d psi=%.17g nc=%d width=%.17g blend=%.17g "
        "sigma_lo=%.17g sigma_hi=%.17g length=%.17g s0=%.17g s1=%.17g "
        "chord_target_lo=%.17g chord_target_hi=%.17g dtheta=%.17g "
        "ordered=%d ext=%d ctx_arc=%.17g v_first=%.17g\n",
        inp->ftype, ctx->psi, nc, width, blend, sigma[0], sigma[1], length,
        y[1], length-y[2], sigma[0]*width*blend, sigma[1]*width*blend,
        width/nc,
        (int) tok_xpt_ordered_placement(inp), (int) tok_ext_construction(inp),
        ctx->arc_hi-ctx->arc_lo, y[1]/length);
    // Does the trace the cubic is scaled by actually span the BLOCK? Print its
    // two endpoints; compare against the grid row's first and last node. If they
    // differ, the trace is longer than the row and `length` is the wrong total;
    // if they match, the row/length gap is only the grid's chord-sum shortfall.
    if (tok_row_arc_plan_dump() && ctx->map_trace_n > 1)
      fprintf(stderr,
        "TOK_ROW_ARC_TRACE ftype=%d psi=%.17g n=%d first=(%.17g,%.17g) "
        "last=(%.17g,%.17g) s_total=%.17g\n",
        inp->ftype, ctx->psi, ctx->map_trace_n,
        ctx->map_trace_r[0], ctx->map_trace_z[0],
        ctx->map_trace_r[ctx->map_trace_n-1], ctx->map_trace_z[ctx->map_trace_n-1],
        ctx->map_trace_s[ctx->map_trace_n-1]);
  }
  if (u <= 0.0) { *v=0.0; *derivative=sigma[0]*width/length; return true; }
  if (u >= 1.0) { *v=1.0; *derivative=sigma[1]*width/length; return true; }
  const double bl=tok_row_arc_blend_fraction(nc);
  const double x[4]={0.0,bl,1.0-bl,1.0};
  int i=u < x[1] ? 0 : (u < x[2] ? 1:2);
  double h=x[i+1]-x[i], t=(u-x[i])/h;
  const double *y=ctx->row_arc_fit_y, *m=ctx->row_arc_fit_m;
  double d=(y[i+1]-y[i])/h, c=3.0*d-2.0*m[i]-m[i+1], b=m[i]+m[i+1]-2.0*d;
  *v=(y[i]+h*t*(m[i]+t*(c+t*b)))/length;
  *derivative=(m[i]+t*(2.0*c+3.0*t*b))/length;
  if (!isfinite(*v) || !isfinite(*derivative) || !(*derivative > 0.0)) {
    fprintf(stderr,"TOK_ROW_ARC_FAILED reason=nonpositive_derivative ftype=%d psi=%.17g u=%.17g segment=%d v=%.17g derivative=%.17g\n",
      inp->ftype,ctx->psi,u,i,*v,*derivative);
    return false;
  }
  return true;
}
