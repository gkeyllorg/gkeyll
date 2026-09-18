#pragma once

// Included in gyrokinetic_multib.c after singleb_app_new_geom_from_block.
// One complete no-write probe per block; no optimizer iteration and no ftype
// singleton. Final geometry consumes shared physical ds/dtheta target curves.
struct row_arc_sample { double psi, length; };
struct row_arc_curve { int n; double *psi, *sigma; };
struct row_arc_block {
  int bid, n, capacity, cells;
  double width;
  double physical_psi[2]; // Actual radial boundaries, distinct from FD halos.
  // How close a psi has to be to one of those boundaries to BE it. Grid
  // arithmetic reaches a boundary as lower + i*(upper-lower)/cells, which lands
  // on the declared double only when the two ends are of comparable magnitude.
  // On TCV psi_sep is 1.4e-5 against a 9.1e-3 span, and the extreme node misses
  // it by 511 ulps OF THE BOUNDARY -- still under one eps of the SPAN the
  // arithmetic was done in. Sized from that span, this is ~1e-14 of one radial
  // cell, so it can only ever absorb rounding.
  double tol;
  bool active, invalid;
  struct row_arc_sample *sample;
  struct row_arc_curve *edge[2];
  // Log corrections at lower/upper radial boundaries, per theta endpoint.
  // Radial partners share normalized H there, accepting the explicitly
  // diagnosed physical-rate residual due to their differing L/width.
  double log_correction[2][2];
  // Exact per-row corner target. The whole boundary solve below runs at the two
  // radial boundaries only, and row_arc_rates blends the two log corrections
  // linearly in t. Measured on NSTX-U, that blend is what the theta seam is: the
  // element ratio is ~1 at both radial ends and bulges to 1.0437 between them,
  // identically at 151 and 601 rows. The constraint graph's potentials are
  // psi-INDEPENDENT -- its offsets are pure dtheta ratios -- so only one scalar
  // per connected component varies with psi, and the solve can be evaluated at
  // ANY row for the cost of one min over the component.
  double theta_potential[2];
  int theta_component[2];
  struct row_arc_plan *plan;
};
struct row_arc_plan {
  int nblocks, ncurves;
  struct row_arc_block *block;
  struct row_arc_curve *curve;
  // Theta-adjacency components, as (block,theta edge) pairs encoded 2*b+e.
  int ncomponents;
  int *comp_member, *comp_begin;
};
struct row_arc_constraint { int from, to; double offset; };

static bool
row_arc_plan_enabled(void)
{
  const char *e = getenv("GKYL_TOK_ROW_ARC");
  return e && e[0] && e[0] != '0';
}

// Evaluate the corner solve at every row instead of at the two radial
// boundaries and blending. Default off while it is under measurement.
static bool
row_arc_exact_target(void)
{
  const char *e = getenv("GKYL_TOK_ROW_ARC_EXACT_TARGET");
  return e && e[0] && e[0] != '0';
}

static void
row_arc_plan_release(void *vplan)
{
  struct row_arc_plan *p = vplan;
  for (int b=0; b<p->nblocks; ++b) gkyl_free(p->block[b].sample);
  for (int c=0; c<p->ncurves; ++c) {
    gkyl_free(p->curve[c].psi); gkyl_free(p->curve[c].sigma);
  }
  gkyl_free(p->comp_member); gkyl_free(p->comp_begin);
  gkyl_free(p->curve); gkyl_free(p->block); gkyl_free(p);
}

// Name a boundary row by its declared psi wherever one is matched, so the
// declared value and the value the grid actually reaches are one knot and not
// two 511 ulps apart. Two knots that close put a step in the rate curve right
// where the graph writes its endpoint target, and the construction then reads
// the uncorrected side of it.
static double
row_arc_canonical_psi(const struct row_arc_block *b, double psi)
{
  for (int r=0; r<2; ++r)
    if (fabs(psi-b->physical_psi[r]) <= b->tol) return b->physical_psi[r];
  return psi;
}

static void
row_arc_capture(void *ctx, double psi, double length, double width, int cells)
{
  struct row_arc_block *b = ctx;
  psi = row_arc_canonical_psi(b, psi);
  if (!isfinite(psi) || !isfinite(length) || !isfinite(width) ||
      !(length > 0.0) || !(width > 0.0) || cells < 1) {
    b->invalid = true;
    return;
  }
  if (b->n && (width != b->width || cells != b->cells)) {
    b->invalid = true;
    return;
  }
  b->width = width; b->cells = cells;
  // The legacy hook is called at every theta point. Store a row once when
  // successive calls report the same psi and length.
  if (b->n && b->sample[b->n-1].psi == psi &&
      b->sample[b->n-1].length == length) return;
  if (b->n == b->capacity) {
    int capacity = b->capacity ? 2*b->capacity : 1;
    struct row_arc_sample *s = gkyl_malloc(capacity*sizeof(*s));
    if (b->n) memcpy(s, b->sample, b->n*sizeof(*s));
    gkyl_free(b->sample); b->sample = s; b->capacity = capacity;
  }
  b->sample[b->n++] = (struct row_arc_sample) { psi, length };
}

static int
row_arc_sample_compare(const void *pa, const void *pb)
{
  const struct row_arc_sample *a=pa, *b=pb;
  if (a->psi != b->psi) return (a->psi > b->psi) - (a->psi < b->psi);
  return (a->length > b->length) - (a->length < b->length);
}

static int
row_arc_double_compare(const void *pa, const void *pb)
{
  double a=*(const double *)pa, b=*(const double *)pb;
  return (a>b)-(a<b);
}

static bool
row_arc_same_psi(double a, double b)
{
  return a == b || nextafter(a,b) == b;
}

// Positive interpolation, with exact values at captured rows. No extrapolation:
// all physical, Gauss, and finite-difference psi requests came through the
// probe. Only one-ulp endpoint differences in coordinate arithmetic are snapped.
static bool
row_arc_curve_eval(const struct row_arc_curve *c, double psi, double *sigma)
{
  if (!isfinite(psi) || c->n < 2) return false;
  int lo=0, hi=c->n-1;
  if (row_arc_same_psi(psi,c->psi[lo])) { *sigma=c->sigma[lo]; return true; }
  if (row_arc_same_psi(psi,c->psi[hi])) { *sigma=c->sigma[hi]; return true; }
  if (psi < c->psi[lo] || psi > c->psi[hi]) return false;
  while (hi-lo > 1) {
    int mid=lo+(hi-lo)/2;
    if (c->psi[mid] <= psi) lo=mid; else hi=mid;
  }
  if (row_arc_same_psi(psi,c->psi[lo])) { *sigma=c->sigma[lo]; return true; }
  if (row_arc_same_psi(psi,c->psi[hi])) { *sigma=c->sigma[hi]; return true; }
  double t=(psi-c->psi[lo])/(c->psi[hi]-c->psi[lo]);
  *sigma=exp((1.0-t)*log(c->sigma[lo])+t*log(c->sigma[hi]));
  return isfinite(*sigma) && *sigma > 0.0;
}

static bool
row_arc_base_rate(const struct row_arc_block *b, double psi, double *rate)
{
  if (b->n < 2 || !isfinite(psi)) return false;
  int lo=0, hi=b->n-1;
  if (row_arc_same_psi(psi,b->sample[lo].psi)) {
    *rate=b->sample[lo].length/b->width; return true;
  }
  if (row_arc_same_psi(psi,b->sample[hi].psi)) {
    *rate=b->sample[hi].length/b->width; return true;
  }
  // The common theta target uses the union of both peers' captured FD rows.
  // If one peer's halo is shorter, extend its positive local rate model with
  // the adjacent segment only during plan construction. Final callbacks never
  // extrapolate beyond the explicitly stored union of actual probe queries.
  if (psi < b->sample[lo].psi) hi=1;
  else if (psi > b->sample[hi].psi) lo=hi-1;
  while (hi-lo > 1) {
    int mid=lo+(hi-lo)/2;
    if (b->sample[mid].psi <= psi) lo=mid; else hi=mid;
  }
  if (row_arc_same_psi(psi,b->sample[lo].psi)) {
    *rate=b->sample[lo].length/b->width; return true;
  }
  if (row_arc_same_psi(psi,b->sample[hi].psi)) {
    *rate=b->sample[hi].length/b->width; return true;
  }
  double t=(psi-b->sample[lo].psi)/(b->sample[hi].psi-b->sample[lo].psi);
  *rate=exp((1.0-t)*log(b->sample[lo].length)+t*log(b->sample[hi].length))/b->width;
  return isfinite(*rate) && *rate > 0.0;
}

// T_{b,e}(psi) = exp(scale(psi) + potential_{b,e}), the largest common end-cell
// arc that satisfies every theta-adjacency constraint arc_peer/arc_self =
// dtheta_peer/dtheta_self and stays inside each member's own uniform cell.
// Exactly the quantity the boundary solve computes, evaluated at psi instead of
// at a radial boundary. sigma = cells*T/width, so the block's own row length
// cancels out of everything but `scale`.
static bool
row_arc_exact_sigma(const struct row_arc_block *b, double psi, int e, double *sigma)
{
  const struct row_arc_plan *p=b->plan;
  if (!p || b->theta_component[e] < 0) return false;
  const int c=b->theta_component[e];
  double scale=DBL_MAX;
  int members=0;
  for (int k=p->comp_begin[c]; k<p->comp_begin[c+1]; ++k) {
    const struct row_arc_block *m=&p->block[p->comp_member[k]/2];
    const int me=p->comp_member[k]%2;
    double rate;
    // A member whose captured support does not reach this row cannot bound it.
    // Theta mates share a radial extent by construction (the plan refuses
    // otherwise, reason=theta_peer_flux_span), so in practice every member is
    // present; this only guards the halo edges.
    if (!row_arc_base_rate(m,psi,&rate)) continue;
    const double cap=rate*m->width/m->cells;   // that member's uniform cell arc
    if (!(cap > 0.0) || !isfinite(cap)) continue;
    const double s=log(cap)-m->theta_potential[me];
    if (s < scale) scale=s;
    ++members;
  }
  if (!members || scale == DBL_MAX) return false;
  *sigma=(b->cells/b->width)*exp(scale+b->theta_potential[e]);
  return isfinite(*sigma) && *sigma > 0.0;
}

static bool
row_arc_rates(void *ctx, double psi, double *sigma_lo, double *sigma_hi)
{
  const struct row_arc_block *b=ctx;
  psi=row_arc_canonical_psi(b,psi);
  if (row_arc_exact_target()) {
    if (!row_arc_exact_sigma(b,psi,0,sigma_lo) ||
        !row_arc_exact_sigma(b,psi,1,sigma_hi)) {
      fprintf(stderr,"TOK_ROW_ARC_RATES_FAILED block=%d psi=%.17g reason=exact_target_unavailable\n",
        b->bid,psi);
      return false;
    }
    if (getenv("GKYL_TOK_ROW_ARC_TRACE_SIGMA")) {
      double own=0.0; bool have=row_arc_base_rate(b,psi,&own);
      fprintf(stderr,"TOK_ROW_ARC_SIGMA block=%d psi=%.17g t=%.17g sigma_lo=%.17g sigma_hi=%.17g own=%.17g\n",
        b->bid,psi,(psi-b->physical_psi[0])/(b->physical_psi[1]-b->physical_psi[0]),
        *sigma_lo,*sigma_hi,have?own:-1.0);
    }
    return true;
  }
  if (!b->edge[0] || !b->edge[1] ||
      !row_arc_curve_eval(b->edge[0],psi,sigma_lo) ||
      !row_arc_curve_eval(b->edge[1],psi,sigma_hi)) {
    fprintf(stderr,"TOK_ROW_ARC_RATES_FAILED block=%d psi=%.17g reason=missing_row_or_extrapolation\n",
      b->bid,psi);
    return false;
  }
  double t=(psi-b->physical_psi[0])/(b->physical_psi[1]-b->physical_psi[0]);
  // Continue the same smooth correction through captured FD halos. Clamping
  // at a physical boundary would change its radial derivative.
  *sigma_lo *= exp((1.0-t)*b->log_correction[0][0]+t*b->log_correction[1][0]);
  *sigma_hi *= exp((1.0-t)*b->log_correction[0][1]+t*b->log_correction[1][1]);
  // Opt-in trace of what this block is actually told to do, against the rate it
  // would have chosen alone. The difference is pinned to zero where theta mates
  // share a row and jumps off it; that derivative is the thing under test, and
  // it cannot be read from the written grid because the grid only shows the
  // imposed result.
  if (getenv("GKYL_TOK_ROW_ARC_TRACE_SIGMA")) {
    double own=0.0;
    bool have=row_arc_base_rate(b,psi,&own);
    fprintf(stderr,"TOK_ROW_ARC_SIGMA block=%d psi=%.17g t=%.17g sigma_lo=%.17g sigma_hi=%.17g own=%.17g\n",
      b->bid,psi,t,*sigma_lo,*sigma_hi,have?own:-1.0);
  }
  return isfinite(*sigma_lo) && isfinite(*sigma_hi) && *sigma_lo > 0.0 && *sigma_hi > 0.0;
}

static int
row_arc_root(int *parent, int x)
{
  while (parent[x] != x) { parent[x]=parent[parent[x]]; x=parent[x]; }
  return x;
}

static void
row_arc_join(int *parent, int a, int b)
{
  a=row_arc_root(parent,a); b=row_arc_root(parent,b);
  if (a != b) parent[a>b ? a:b]=a>b ? b:a;
}

static int
row_arc_target_side(int edge)
{
  switch (edge) {
    case GKYL_LOWER_POSITIVE: case GKYL_LOWER_NEGATIVE: return 0;
    case GKYL_UPPER_POSITIVE: case GKYL_UPPER_NEGATIVE: return 1;
    default: return -1;
  }
}

static bool
row_arc_make_curves(struct row_arc_plan *p, struct gkyl_gk_block_geom *bg)
{
  const int n=p->nblocks;
  int *ends=gkyl_malloc(2*n*sizeof(*ends));
  int *corners=gkyl_malloc(4*n*sizeof(*corners));
  int *radial_corners=gkyl_malloc(4*n*sizeof(*radial_corners));
  int *radial_rows=gkyl_malloc(2*n*sizeof(*radial_rows));
  struct row_arc_constraint *constraint=gkyl_malloc(4*n*sizeof(*constraint));
  int nc=0;
  double *canonical_length=gkyl_calloc(2*n,sizeof(*canonical_length));
  double *potential=gkyl_calloc(4*n,sizeof(*potential));
  double *cap=gkyl_malloc(4*n*sizeof(*cap));
  double *target=gkyl_calloc(4*n,sizeof(*target));
  double *desired_sigma=gkyl_calloc(4*n,sizeof(*desired_sigma));
  bool *seen=gkyl_calloc(4*n,sizeof(*seen));
  int *queue=gkyl_malloc(4*n*sizeof(*queue));
  bool ok=false;
  for (int i=0; i<2*n; ++i) { ends[i]=i; radial_rows[i]=i; }
  for (int i=0; i<4*n; ++i) { corners[i]=i; radial_corners[i]=i; cap[i]=DBL_MAX; }
  for (int b=0; b<n; ++b) {
    if (!p->block[b].active) continue;
    const struct gkyl_gk_block_geom_info *bi=gkyl_gk_block_geom_get_block(bg,b);
    for (int d=0; d<2; ++d) for (int e=0; e<2; ++e) {
      const struct gkyl_target_edge *te=&bi->connections[d][e];
      if (te->edge == GKYL_PHYSICAL) continue;
      int c=te->bid, ce=row_arc_target_side(te->edge);
      if (c < 0 || c >= n || te->dir != d || ce < 0 || !p->block[c].active) {
        fprintf(stderr,"TOK_ROW_ARC_PLAN_FAILED reason=unsupported_connection block=%d dir=%d edge=%d\n",b,d,e);
        goto done;
      }
      bool reverse=te->edge == GKYL_LOWER_NEGATIVE || te->edge == GKYL_UPPER_NEGATIVE;
      if (d == 1) row_arc_join(ends,2*b+e,2*c+ce);
      else {
        if (p->block[b].cells != p->block[c].cells) {
          fprintf(stderr,"TOK_ROW_ARC_PLAN_FAILED reason=radial_theta_cell_count block=%d peer=%d\n",b,c);
          goto done;
        }
        row_arc_join(radial_rows,2*b+e,2*c+ce);
      }
      for (int t=0; t<2; ++t) {
        int ct=reverse ? 1-t:t;
        int br=d == 0 ? e:t, bt=d == 0 ? t:e;
        int cr=d == 0 ? ce:ct, tt=d == 0 ? ct:ce;
        row_arc_join(corners,4*b+2*br+bt,4*c+2*cr+tt);
        if (d == 0) row_arc_join(radial_corners,4*b+2*br+bt,4*c+2*cr+tt);
      }
    }
  }
  for (int root=0; root<2*n; ++root) {
    if (!p->block[root/2].active || row_arc_root(ends,root) != root) continue;
    int total=0, members=0;
    double psi_lo=0.0, psi_hi=0.0;
    for (int e=0; e<2*n; ++e) if (row_arc_root(ends,e) == root) {
      const struct row_arc_block *b=&p->block[e/2];
      if (!members) { psi_lo=b->physical_psi[0]; psi_hi=b->physical_psi[1]; }
      else if (!row_arc_same_psi(psi_lo,b->physical_psi[0]) ||
               !row_arc_same_psi(psi_hi,b->physical_psi[1])) {
        fprintf(stderr,"TOK_ROW_ARC_PLAN_FAILED reason=theta_peer_flux_span block=%d\n",b->bid);
        goto done;
      }
      total+=b->n+2; ++members;
    }
    struct row_arc_curve *curve=&p->curve[p->ncurves++];
    curve->psi=gkyl_malloc(total*sizeof(*curve->psi));
    curve->sigma=gkyl_malloc(total*sizeof(*curve->sigma));
    for (int e=0; e<2*n; ++e) if (row_arc_root(ends,e) == root) {
      struct row_arc_block *b=&p->block[e/2];
      b->edge[e%2]=curve;
      for (int j=0; j<b->n; ++j) curve->psi[curve->n++]=b->sample[j].psi;
      curve->psi[curve->n++]=b->physical_psi[0];
      curve->psi[curve->n++]=b->physical_psi[1];
    }
    qsort(curve->psi,curve->n,sizeof(*curve->psi),row_arc_double_compare);
    int count=0;
    for (int j=0; j<curve->n; ++j)
      if (!count || curve->psi[j] != curve->psi[count-1]) curve->psi[count++]=curve->psi[j];
    curve->n=count;
    for (int j=0; j<curve->n; ++j) {
      double smallest=DBL_MAX;
      for (int e=0; e<2*n; ++e) if (row_arc_root(ends,e) == root) {
        double rate;
        const struct row_arc_block *b=&p->block[e/2];
        if (curve->psi[j] < b->sample[0].psi || curve->psi[j] > b->sample[b->n-1].psi)
          fprintf(stderr,"TOK_ROW_ARC_MODEL_EXTEND block=%d psi=%.17g captured_lo=%.17g captured_hi=%.17g\n",
            b->bid,curve->psi[j],b->sample[0].psi,b->sample[b->n-1].psi);
        if (!row_arc_base_rate(&p->block[e/2],curve->psi[j],&rate)) {
          fprintf(stderr,"TOK_ROW_ARC_PLAN_FAILED reason=no_base_rate stage=theta_curve block=%d psi=%.17g\n",
            b->bid,curve->psi[j]);
          goto done;
        }
        smallest=fmin(smallest,rate);
      }
      // Largest common physical rate no greater than either uniform row rate.
      curve->sigma[j]=smallest;
      if (getenv("GKYL_TOK_ROW_ARC_TRACE_SIGMA"))
        fprintf(stderr,"TOK_ROW_ARC_KNOT curve=%d j=%d psi=%.17g sigma=%.17g\n",
          p->ncurves,j,curve->psi[j],curve->sigma[j]);
      if (!(curve->sigma[j] > 0.0) || !isfinite(curve->sigma[j])) {
        fprintf(stderr,"TOK_ROW_ARC_PLAN_FAILED reason=nonpositive_common_rate psi=%.17g sigma=%.17g\n",
          curve->psi[j],curve->sigma[j]);
        goto done;
      }
    }
  }
  // Check physical corner identity, separately from its rate constraints.
  // A radial join identifies first/last cell ARCS, not physical ds/dtheta.
  for (int root=0; root<4*n; ++root) {
    if (!p->block[root/4].active || row_arc_root(corners,root) != root) continue;
    int members=0;
    double psi=0.0;
    for (int node=0; node<4*n; ++node) if (row_arc_root(corners,node) == root) {
      const struct row_arc_block *b=&p->block[node/4];
      int radial=(node%4)/2;
      double at=b->physical_psi[radial];
      if (members && !row_arc_same_psi(psi,at)) {
        fprintf(stderr,"TOK_ROW_ARC_PLAN_FAILED reason=corner_flux_mismatch block=%d psi=%.17g expected=%.17g\n",b->bid,at,psi);
        goto done;
      }
      psi=at; ++members;
    }
  }
  // Use the lowest-index representative's row length for radial equivalence.
  // The actual curves are not overwritten here: any pre-existing R,Z mismatch
  // remains visible to iface_gate. Both partners receive the same normalized H.
  for (int node=0; node<2*n; ++node) if (p->block[node/2].active) {
    int rep=row_arc_root(radial_rows,node);
    const struct row_arc_block *b=&p->block[rep/2];
    double rate;
    if (!row_arc_base_rate(b,b->physical_psi[rep%2],&rate)) {
      fprintf(stderr,"TOK_ROW_ARC_PLAN_FAILED reason=no_base_rate stage=canonical_length block=%d psi=%.17g\n",
        b->bid,b->physical_psi[rep%2]);
      goto done;
    }
    canonical_length[node]=rate*b->width;
  }
  for (int node=0; node<4*n; ++node) if (p->block[node/4].active) {
    int rep=row_arc_root(radial_corners,node);
    cap[rep]=fmin(cap[rep],canonical_length[node/2]/p->block[node/4].cells);
  }
  // Reduced graph: radial endpoint identities were eliminated exactly above.
  // Theta adjacency requires arc_peer/arc_self=dtheta_peer/dtheta_self.
  for (int b=0; b<n; ++b) if (p->block[b].active) {
    const struct gkyl_gk_block_geom_info *bi=gkyl_gk_block_geom_get_block(bg,b);
    for (int e=0; e<2; ++e) {
      const struct gkyl_target_edge *te=&bi->connections[1][e];
      if (te->edge == GKYL_PHYSICAL) continue;
      int c=te->bid, ce=row_arc_target_side(te->edge);
      bool reverse=te->edge == GKYL_LOWER_NEGATIVE || te->edge == GKYL_UPPER_NEGATIVE;
      double dt=p->block[b].width/p->block[b].cells;
      double peer_dt=p->block[c].width/p->block[c].cells;
      for (int r=0; r<2; ++r) {
        int cr=reverse ? 1-r:r;
        constraint[nc++]=(struct row_arc_constraint) {
          row_arc_root(radial_corners,4*b+2*r+e),
          row_arc_root(radial_corners,4*c+2*cr+ce),log(peer_dt/dt)
        };
      }
    }
  }
  double cycle=0.0;
  for (int root=0; root<4*n; ++root) {
    if (cap[root] == DBL_MAX || seen[root]) continue;
    int head=0, tail=1;
    queue[0]=root; seen[root]=true; potential[root]=0.0;
    while (head < tail) {
      int current=queue[head++];
      for (int k=0; k<nc; ++k) {
        const struct row_arc_constraint *edge=&constraint[k];
        int peer=-1;
        double expected=0.0;
        if (edge->from == current) { peer=edge->to; expected=potential[current]+edge->offset; }
        else if (edge->to == current) { peer=edge->from; expected=potential[current]-edge->offset; }
        if (peer < 0) continue;
        if (seen[peer]) cycle=fmax(cycle,fabs(potential[peer]-expected));
        else { seen[peer]=true; potential[peer]=expected; queue[tail++]=peer; }
      }
    }
    double scale=DBL_MAX;
    for (int j=0; j<tail; ++j) scale=fmin(scale,log(cap[queue[j]])-potential[queue[j]]);
    for (int j=0; j<tail; ++j) target[queue[j]]=exp(scale+potential[queue[j]]);
  }
  fprintf(stderr,"TOK_ROW_ARC_GRAPH max_log_cycle_residual=%.17g radial_endpoint_identity=exact\n",cycle);
  for (int node=0; node<4*n; ++node) if (p->block[node/4].active) {
    const struct row_arc_block *b=&p->block[node/4];
    int r=(node%4)/2, rep=row_arc_root(radial_corners,node);
    double rate;
    if (!row_arc_base_rate(b,b->physical_psi[r],&rate)) {
      fprintf(stderr,"TOK_ROW_ARC_PLAN_FAILED reason=no_base_rate stage=boundary_sigma block=%d psi=%.17g\n",
        b->bid,b->physical_psi[r]);
      goto done;
    }
    double normalized=b->cells*target[rep]/canonical_length[node/2];
    desired_sigma[node]=rate*normalized;
    if (!isfinite(desired_sigma[node]) || !(desired_sigma[node] > 0.0)) {
      fprintf(stderr,"TOK_ROW_ARC_PLAN_FAILED reason=nonpositive_boundary_sigma block=%d radial_edge=%d theta_edge=%d sigma=%.17g\n",
        b->bid,r,node%2,desired_sigma[node]);
      goto done;
    }
    fprintf(stderr,"TOK_ROW_ARC_BOUNDARY block=%d radial_edge=%d theta_edge=%d psi=%.17g canonical_length=%.17g endpoint_arc=%.17g normalized_slope=%.17g sigma=%.17g\n",
      b->bid,r,node%2,b->physical_psi[r],canonical_length[node/2],target[rep],normalized,desired_sigma[node]);
  }
  // At each theta boundary, use the graph's common physical target. Only actual
  // theta mates are averaged; radial edges never impose equal physical sigma.
  // Any tiny cycle/private-trace discrepancy is carried in a per-block factor
  // so the normalized radial maps remain identical.
  for (int node=0; node<4*n; ++node) if (p->block[node/4].active) {
    struct row_arc_block *b=&p->block[node/4];
    int r=(node%4)/2, e=node%2;
    double psi=b->physical_psi[r], log_sum=0.0;
    struct row_arc_curve *c=b->edge[e];
    int members=0, index=-1;
    for (int other=0; other<4*n; ++other) if (p->block[other/4].active) {
      const struct row_arc_block *peer=&p->block[other/4];
      if (peer->edge[other%2] == c && row_arc_same_psi(peer->physical_psi[(other%4)/2],psi)) {
        log_sum+=log(desired_sigma[other]); ++members;
      }
    }
    for (int j=0; j<c->n; ++j) if (row_arc_same_psi(psi,c->psi[j])) { index=j; break; }
    if (index < 0 || members < 1) {
      fprintf(stderr,"TOK_ROW_ARC_PLAN_FAILED reason=boundary_row_not_on_curve block=%d radial_edge=%d theta_edge=%d psi=%.17g members=%d\n",
        b->bid,r,e,psi,members);
      goto done;
    }
    c->sigma[index]=exp(log_sum/members);
  }
  for (int node=0; node<4*n; ++node) if (p->block[node/4].active) {
    struct row_arc_block *b=&p->block[node/4];
    int r=(node%4)/2, e=node%2;
    double sigma;
    if (!row_arc_curve_eval(b->edge[e],b->physical_psi[r],&sigma)) {
      fprintf(stderr,"TOK_ROW_ARC_PLAN_FAILED reason=curve_eval stage=log_correction block=%d radial_edge=%d theta_edge=%d psi=%.17g\n",
        b->bid,r,e,b->physical_psi[r]);
      goto done;
    }
    b->log_correction[r][e]=log(desired_sigma[node]/sigma);
  }
  // Theta-adjacency components and their potentials, for the exact per-row
  // target. Same constraint the boundary graph uses -- arc_peer/arc_self =
  // dtheta_peer/dtheta_self -- but over (block,theta edge) nodes only: the
  // radial identifications in `radial_corners` exist at radial boundaries and
  // have no counterpart at an interior row. The offsets are pure dtheta ratios,
  // so the potentials do not depend on psi and are solved once here.
  {
    int *tparent=gkyl_malloc(2*n*sizeof(*tparent));
    bool *tseen=gkyl_calloc(2*n,sizeof(*tseen));
    for (int i=0; i<2*n; ++i) { tparent[i]=i; p->block[i/2].theta_component[i%2]=-1; }
    for (int b=0; b<n; ++b) if (p->block[b].active) {
      p->block[b].plan=p;
      const struct gkyl_gk_block_geom_info *bi=gkyl_gk_block_geom_get_block(bg,b);
      for (int e=0; e<2; ++e) {
        const struct gkyl_target_edge *te=&bi->connections[1][e];
        if (te->edge == GKYL_PHYSICAL) continue;
        int c=te->bid, ce=row_arc_target_side(te->edge);
        if (c < 0 || c >= n || ce < 0 || !p->block[c].active) continue;
        row_arc_join(tparent,2*b+e,2*c+ce);
      }
    }
    p->comp_member=gkyl_malloc(2*n*sizeof(*p->comp_member));
    p->comp_begin=gkyl_malloc((2*n+1)*sizeof(*p->comp_begin));
    p->ncomponents=0;
    int filled=0;
    for (int root=0; root<2*n; ++root) {
      if (!p->block[root/2].active || row_arc_root(tparent,root) != root) continue;
      const int cid=p->ncomponents++;
      p->comp_begin[cid]=filled;
      // BFS the component, carrying log(dtheta) offsets.
      int head=filled;
      p->comp_member[filled++]=root;
      tseen[root]=true;
      p->block[root/2].theta_potential[root%2]=0.0;
      p->block[root/2].theta_component[root%2]=cid;
      while (head < filled) {
        const int cur=p->comp_member[head++];
        const int cb=cur/2, cel=cur%2;
        const struct gkyl_gk_block_geom_info *bi=gkyl_gk_block_geom_get_block(bg,cb);
        const struct gkyl_target_edge *te=&bi->connections[1][cel];
        if (te->edge == GKYL_PHYSICAL) continue;
        const int c=te->bid, ce=row_arc_target_side(te->edge);
        if (c < 0 || c >= n || ce < 0 || !p->block[c].active) continue;
        const int peer=2*c+ce;
        if (tseen[peer]) continue;
        const double dt=p->block[cb].width/p->block[cb].cells;
        const double peer_dt=p->block[c].width/p->block[c].cells;
        tseen[peer]=true;
        p->block[c].theta_potential[ce]=
          p->block[cb].theta_potential[cel]+log(peer_dt/dt);
        p->block[c].theta_component[ce]=cid;
        p->comp_member[filled++]=peer;
      }
    }
    p->comp_begin[p->ncomponents]=filled;
    for (int b=0; b<n; ++b) if (p->block[b].active) for (int e=0; e<2; ++e)
      if (p->block[b].theta_component[e] >= 0)
        fprintf(stderr,"TOK_ROW_ARC_THETA_COMPONENT block=%d edge=%d component=%d potential=%.17g dtheta=%.17g\n",
          b,e,p->block[b].theta_component[e],p->block[b].theta_potential[e],
          p->block[b].width/p->block[b].cells);
    gkyl_free(tparent); gkyl_free(tseen);
  }
  for (int root=0; root<2*n; ++root) {
    if (!p->block[root/2].active || row_arc_root(ends,root) != root) continue;
    const struct row_arc_curve *curve=p->block[root/2].edge[root%2];
    double max_ratio=1.0, at_psi=curve->psi[0];
    int members=0;
    for (int e=0; e<2*n; ++e) if (row_arc_root(ends,e) == root) ++members;
    if (members < 2) continue;
    for (int j=0; j<curve->n; ++j) {
      double lo=DBL_MAX, hi=0.0;
      for (int e=0; e<2*n; ++e) if (row_arc_root(ends,e) == root) {
        double lower,upper;
        if (!row_arc_rates(&p->block[e/2],curve->psi[j],&lower,&upper)) goto done;
        double rate=e%2 ? upper:lower;
        lo=fmin(lo,rate); hi=fmax(hi,rate);
      }
      if (hi/lo > max_ratio) { max_ratio=hi/lo; at_psi=curve->psi[j]; }
    }
    fprintf(stderr,"TOK_ROW_ARC_THETA_RESIDUAL root=%d members=%d max_ratio=%.17g psi=%.17g radial_map_preserved=1\n",
      root,members,max_ratio,at_psi);
  }
  ok=true;
done:
  gkyl_free(queue); gkyl_free(seen); gkyl_free(desired_sigma);
  gkyl_free(target); gkyl_free(cap); gkyl_free(potential);
  gkyl_free(canonical_length); gkyl_free(constraint); gkyl_free(radial_rows);
  gkyl_free(radial_corners); gkyl_free(corners); gkyl_free(ends);
  return ok;
}

static bool
gyrokinetic_multib_prepare_row_arc(const struct gkyl_gyrokinetic_multib *inp,
  struct gkyl_gk_block_geom *bg)
{
  int n=gkyl_gk_block_geom_num_blocks(bg), ndim=gkyl_gk_block_geom_ndim(bg), ranks=0;
  gkyl_comm_get_size(inp->comm,&ranks);
  if (n < 1 || ndim != 2 || inp->poly_order != 1 || ranks != 1) {
    fprintf(stderr,"TOK_ROW_ARC_PLAN_FAILED reason=prototype_scope blocks=%d ndim=%d poly_order=%d ranks=%d\n",
      n,ndim,inp->poly_order,ranks);
    return false;
  }
  struct row_arc_plan *p=gkyl_calloc(1,sizeof(*p));
  p->nblocks=n;
  p->block=gkyl_calloc(n,sizeof(*p->block));
  p->curve=gkyl_calloc(2*n,sizeof(*p->curve));
  struct gkyl_gyrokinetic_multib *probe_inp=gkyl_malloc(sizeof(*probe_inp));
  *probe_inp=*inp; probe_inp->gk_block_geom=bg; probe_inp->use_gpu=false;
  struct gkyl_gyrokinetic_multib_app *probe=gkyl_calloc(1,sizeof(*probe));
  probe->gk_block_geom=bg;
  probe->block_comms=gkyl_calloc(n,sizeof(*probe->block_comms));
  struct gkyl_comm *serial=gkyl_null_comm_inew(&(struct gkyl_null_comm_inp) { .use_gpu=false });
  bool ok=false;
  for (int b=0; b<n; ++b) probe->block_comms[b]=serial;
  for (int b=0; b<n; ++b) {
    struct gkyl_gk_block_geom_info bi=*gkyl_gk_block_geom_get_block(bg,b);
    struct row_arc_block *capture=&p->block[b];
    capture->bid=b;
    if (bi.geometry.geometry_id != GKYL_GEOMETRY_TOKAMAK) continue;
    const struct gkyl_position_map_inp *pm=&bi.geometry.position_map_info;
    if ((pm->id != GKYL_PMAP_USER_INPUT && pm->id != GKYL_PMAP_USER_INPUT_W_DERIVATIVE) || pm->maps[0]) {
      fprintf(stderr,"TOK_ROW_ARC_PLAN_FAILED reason=prototype_requires_identity_radial_map block=%d map_id=%d\n",b,pm->id);
      goto done;
    }
    capture->physical_psi[0]=bi.lower[0]; capture->physical_psi[1]=bi.upper[0];
    capture->tol=(bi.cells[0]+3)*DBL_EPSILON*
      fmax(fabs(bi.lower[0]),fabs(bi.upper[0]));
    // Which radial END carries the larger psi is the device's, not ours: psi
    // falls outward on TCV, so its CORE block is declared from psi_max_core
    // down to psi_sep.  Everything downstream indexes by radial edge and looks
    // rows up by psi, so only a degenerate extent is a defect.  Requiring an
    // increasing declaration here rejected that block with no message at all.
    if (!isfinite(bi.lower[0]) || !isfinite(bi.upper[0]) ||
        bi.upper[0] == bi.lower[0]) {
      fprintf(stderr,"TOK_ROW_ARC_PLAN_FAILED reason=degenerate_radial_extent block=%d lower=%.17g upper=%.17g\n",
        b,bi.lower[0],bi.upper[0]);
      goto done;
    }
    capture->active=true;
    for (int d=0; d<ndim; ++d) bi.cuts[d]=1;
    bi.geometry.tok_grid_info.row_arc_ctx=capture;
    bi.geometry.tok_grid_info.row_arc_capture=row_arc_capture;
    bi.geometry.tok_grid_info.row_arc_rates=0;
    fprintf(stderr,"TOK_ROW_ARC_PROBE block=%d event=begin\n",b);
    struct gkyl_gyrokinetic_app *app=singleb_app_new_geom_from_block(probe_inp,b,probe,&bi,false);
    if (!app) {
      fprintf(stderr,"TOK_ROW_ARC_PLAN_FAILED reason=probe_geometry_build block=%d\n",b);
      goto done;
    }
    gkyl_gyrokinetic_app_release_geom(app);
    if (capture->invalid || capture->n < 2) {
      fprintf(stderr,"TOK_ROW_ARC_PLAN_FAILED reason=missing_or_invalid_capture block=%d rows=%d\n",b,capture->n);
      goto done;
    }
    qsort(capture->sample,capture->n,sizeof(*capture->sample),row_arc_sample_compare);
    int count=0;
    double duplicate_ratio=1.0, duplicate_psi=capture->sample[0].psi;
    for (int j=0; j<capture->n;) {
      int k=j+1;
      long double total=capture->sample[j].length;
      while (k<capture->n && capture->sample[k].psi == capture->sample[j].psi)
        total+=capture->sample[k++].length;
      double ratio=capture->sample[k-1].length/capture->sample[j].length;
      if (ratio > duplicate_ratio) { duplicate_ratio=ratio; duplicate_psi=capture->sample[j].psi; }
      capture->sample[count++]=(struct row_arc_sample) {capture->sample[j].psi,(double)(total/(k-j))};
      j=k;
    }
    capture->n=count;
    if (count < 2) {
      fprintf(stderr,"TOK_ROW_ARC_PLAN_FAILED reason=single_distinct_row block=%d rows=%d\n",b,count);
      goto done;
    }
    // The samples are sorted ascending in psi; the declaration may run either
    // way, so the census is against the extent's min and max.
    double psi_min=fmin(capture->physical_psi[0],capture->physical_psi[1]);
    double psi_max=fmax(capture->physical_psi[0],capture->physical_psi[1]);
    if ((capture->sample[0].psi > psi_min && !row_arc_same_psi(capture->sample[0].psi,psi_min)) ||
        (capture->sample[count-1].psi < psi_max && !row_arc_same_psi(capture->sample[count-1].psi,psi_max))) {
      fprintf(stderr,"TOK_ROW_ARC_PLAN_FAILED reason=physical_boundary_not_captured block=%d captured=[%.17g,%.17g] physical=[%.17g,%.17g]\n",
        b,capture->sample[0].psi,capture->sample[count-1].psi,psi_min,psi_max);
      goto done;
    }
    fprintf(stderr,"TOK_ROW_ARC_PROBE block=%d event=end rows=%d support_lo=%.17g support_hi=%.17g physical_lo=%.17g physical_hi=%.17g width=%.17g cells=%d duplicate_length_ratio=%.17g duplicate_psi=%.17g\n",
      b,count,capture->sample[0].psi,capture->sample[count-1].psi,capture->physical_psi[0],capture->physical_psi[1],
      capture->width,capture->cells,duplicate_ratio,duplicate_psi);
  }
  if (!row_arc_make_curves(p,bg)) goto done;
  for (int b=0; b<n; ++b) if (p->block[b].active) {
    struct gkyl_gk_block_geom_info bi=*gkyl_gk_block_geom_get_block(bg,b);
    bi.geometry.tok_grid_info.row_arc_ctx=&p->block[b];
    bi.geometry.tok_grid_info.row_arc_capture=0;
    bi.geometry.tok_grid_info.row_arc_rates=row_arc_rates;
    gkyl_gk_block_geom_set_block(bg,b,&bi);
  }
  gkyl_gk_block_geom_set_row_arc_owner(bg,p,row_arc_plan_release);
  fprintf(stderr,"TOK_ROW_ARC_PLAN_READY blocks=%d curves=%d probes_per_block=1\n",n,p->ncurves);
  ok=true;
done:
  gkyl_comm_release(serial);
  gkyl_free(probe->block_comms); gkyl_free(probe); gkyl_free(probe_inp);
  if (!ok) row_arc_plan_release(p);
  return ok;
}
