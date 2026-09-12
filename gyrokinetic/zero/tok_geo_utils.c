#include <gkyl_tok_geo_priv.h>
#include <gkyl_tok_geo_wall_target_priv.h>
#include <assert.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

// Helper functions for finding turning points when necessary

void find_lower_turning_point(struct gkyl_tok_geo *geo, double psi_curr,
  double zup, double *zmin, double tolerance);

static double
tok_eval_psi_rz(const struct gkyl_tok_geo *geo, double R, double Z)
{
  if (geo->use_cubics) {
    double xn[2] = { R, Z }, psi;
    geo->efit->evf->eval_cubic(0.0, xn, &psi, geo->efit->evf->ctx);
    return psi;
  }

  int rzidx[2];
  rzidx[0] = fmin(geo->rzlocal.lower[0]
      + (int) floor((R-geo->rzgrid.lower[0])/geo->rzgrid.dx[0]),
    geo->rzlocal.upper[0]);
  rzidx[1] = fmin(geo->rzlocal.lower[1]
      + (int) floor((Z-geo->rzgrid.lower[1])/geo->rzgrid.dx[1]),
    geo->rzlocal.upper[1]);
  rzidx[0] = fmax(rzidx[0], geo->rzlocal.lower[0]);
  rzidx[1] = fmax(rzidx[1], geo->rzlocal.lower[1]);

  long loc = gkyl_range_idx(&geo->rzlocal, rzidx);
  const double *coeffs = gkyl_array_cfetch(geo->psiRZ, loc);
  double xc[2], xy[2];
  gkyl_rect_grid_cell_center(&geo->rzgrid, rzidx, xc);
  xy[0] = (R-xc[0])/(0.5*geo->rzgrid.dx[0]);
  xy[1] = (Z-xc[1])/(0.5*geo->rzgrid.dx[1]);
  return geo->rzbasis.eval_expand(xy, coeffs);
}

static double
tok_wall_tolerance(const struct gkyl_efit *e)
{
  return 1e-8*fmax(1.0,fmax(e->rdim,e->zdim));
}

bool
tok_wall_point_inside(const struct gkyl_efit *e, const double p[2])
{
  // "No usable outline" answers false here because the question -- is this point
  // inside the wall? -- has no answer without one. Callers must therefore not
  // reach this when the outline is acknowledged absent; they gate on
  // gkyl_tok_wall_policy_for first.
  if (!gkyl_tok_wall_usable(e) || !tok_geo_finite(p[0]) || !tok_geo_finite(p[1])) return false;
  for (int i=0; i<e->limiter_n; ++i)
    if (!tok_geo_finite(e->limiter_R[i]) || !tok_geo_finite(e->limiter_Z[i])) return false;
  bool inside=false;
  double tol=tok_wall_tolerance(e);
  // EQDSK supplies an ordered vessel outline. Its last-to-first edge closes
  // the polygon; a repeated first vertex is optional.
  for (int i=0,j=e->limiter_n-1; i<e->limiter_n; j=i++) {
    double ax=e->limiter_R[j],ay=e->limiter_Z[j];
    double bx=e->limiter_R[i],by=e->limiter_Z[i];
    if (!tok_geo_finite(ax) || !tok_geo_finite(ay) || !tok_geo_finite(bx) || !tok_geo_finite(by)) return false;
    double dx=bx-ax,dy=by-ay,l2=dx*dx+dy*dy;
    if (l2>0.0) {
      double t=fmax(0.0,fmin(1.0,((p[0]-ax)*dx+(p[1]-ay)*dy)/l2));
      if (hypot(p[0]-ax-t*dx,p[1]-ay-t*dy)<=tol) return true;
    }
    if ((ay>p[1])!=(by>p[1]) && p[0]<ax+(p[1]-ay)*dx/(by-ay)) inside=!inside;
  }
  return inside;
}

static int
tok_wall_compare_double(const void *a, const void *b)
{
  double x=*(const double *)a,y=*(const double *)b;
  return (x>y)-(x<y);
}

bool
tok_wall_segment_inside(const struct gkyl_efit *e,
  const double a[2], const double b[2])
{
  if (!tok_wall_point_inside(e,a) || !tok_wall_point_inside(e,b)) return false;
  double dx=b[0]-a[0],dy=b[1]-a[1],l2=dx*dx+dy*dy;
  double tol=tok_wall_tolerance(e);
  if (l2<=tol*tol) return true;
  // Partition at EVERY intersection (including collinear wall vertices), then
  // test each open interval. Endpoints alone miss excursions across a concavity.
  double *cuts=gkyl_malloc((2*e->limiter_n+2)*sizeof(double));
  int n=0;cuts[n++]=0.0;cuts[n++]=1.0;
  for (int i=0,j=e->limiter_n-1; i<e->limiter_n; j=i++) {
    double px=e->limiter_R[j]-a[0],py=e->limiter_Z[j]-a[1];
    double ex=e->limiter_R[i]-e->limiter_R[j],ey=e->limiter_Z[i]-e->limiter_Z[j];
    double den=dx*ey-dy*ex;
    double roundoff=64.0*DBL_EPSILON*fmax(DBL_MIN,sqrt(l2)*hypot(ex,ey));
    if (fabs(den)>roundoff) {
      double t=(px*ey-py*ex)/den,u=(px*dy-py*dx)/den;
      if (t>0.0 && t<1.0 && u>=0.0 && u<=1.0) cuts[n++]=t;
    }
    else if (fabs(px*dy-py*dx)<=tol*sqrt(l2)) {
      double t=(px*dx+py*dy)/l2;
      if (t>0.0 && t<1.0) cuts[n++]=t;
      t=((px+ex)*dx+(py+ey)*dy)/l2;
      if (t>0.0 && t<1.0) cuts[n++]=t;
    }
  }
  qsort(cuts,n,sizeof(double),tok_wall_compare_double);
  bool ok=true;
  for (int i=1; i<n && ok; ++i) {
    double t=0.5*(cuts[i-1]+cuts[i]);
    double p[2]={a[0]+t*dx,a[1]+t*dy};
    ok=tok_wall_point_inside(e,p);
  }
  gkyl_free(cuts);
  return ok;
}

// Distance from a point to the wall outline, as a positive number. Only
// meaningful for a point already known to be outside.
static double
tok_wall_distance_to_outline(const struct gkyl_efit *e, const double p[2])
{
  double best = DBL_MAX;
  for (int i=0,j=e->limiter_n-1; i<e->limiter_n; j=i++) {
    double ax=e->limiter_R[j],ay=e->limiter_Z[j];
    double bx=e->limiter_R[i],by=e->limiter_Z[i];
    double dx=bx-ax,dy=by-ay,l2=dx*dx+dy*dy;
    double t = l2>0.0 ? fmax(0.0,fmin(1.0,((p[0]-ax)*dx+(p[1]-ay)*dy)/l2)) : 0.0;
    double d = hypot(p[0]-ax-t*dx, p[1]-ay-t*dy);
    if (d < best) best = d;
  }
  return best==DBL_MAX ? 0.0 : best;
}

// REPORTING ONLY: how far outside the wall this segment gets, in metres.
//
// This does not decide anything. tok_wall_segment_inside above is the verdict,
// and it is exact; this exists so a rejection can say whether it is a
// millimetre of discretisation overshoot or a quarter metre of misdeclared
// domain. It reuses the same partition so the outside intervals are found
// exactly, then samples within them -- so the figure is a lower bound on the
// true maximum, tight enough to tell those two cases apart.
//
// Runs only on the rejection path, so an O(samples x vertices) scan is cheap.
double
tok_wall_segment_excursion(const struct gkyl_efit *e,
  const double a[2], const double b[2])
{
  if (!gkyl_tok_wall_usable(e)) return 0.0;
  double dx=b[0]-a[0],dy=b[1]-a[1],l2=dx*dx+dy*dy;
  double tol=tok_wall_tolerance(e);
  if (l2<=tol*tol) return 0.0;

  double *cuts=gkyl_malloc((2*e->limiter_n+2)*sizeof(double));
  int n=0;cuts[n++]=0.0;cuts[n++]=1.0;
  for (int i=0,j=e->limiter_n-1; i<e->limiter_n; j=i++) {
    double px=e->limiter_R[j]-a[0],py=e->limiter_Z[j]-a[1];
    double ex=e->limiter_R[i]-e->limiter_R[j],ey=e->limiter_Z[i]-e->limiter_Z[j];
    double den=dx*ey-dy*ex;
    double roundoff=64.0*DBL_EPSILON*fmax(DBL_MIN,sqrt(l2)*hypot(ex,ey));
    if (fabs(den)>roundoff) {
      double t=(px*ey-py*ex)/den,u=(px*dy-py*dx)/den;
      if (t>0.0 && t<1.0 && u>=0.0 && u<=1.0) cuts[n++]=t;
    }
  }
  qsort(cuts,n,sizeof(double),tok_wall_compare_double);

  double worst=0.0;
  const int nsample=64;
  for (int i=1; i<n; ++i) {
    double t0=cuts[i-1],t1=cuts[i];
    double tm=0.5*(t0+t1);
    double pm[2]={a[0]+tm*dx,a[1]+tm*dy};
    if (tok_wall_point_inside(e,pm)) continue;   // this interval is inside
    for (int k=0; k<=nsample; ++k) {
      double t=t0+(t1-t0)*k/(double)nsample;
      double p[2]={a[0]+t*dx,a[1]+t*dy};
      if (tok_wall_point_inside(e,p)) continue;
      double d=tok_wall_distance_to_outline(e,p);
      if (d>worst) worst=d;
    }
  }
  gkyl_free(cuts);
  return worst;
}

// Roots of a*t^2+b*t+c in [0,1], with a stable quadratic formula.
static int
tok_wall_quadratic_roots(double a, double b, double c, double roots[2])
{
  double scale=fmax(fabs(a),fmax(fabs(b),fabs(c)));
  if (scale==0.0) return 0;
  a/=scale;b/=scale;c/=scale;
  int n=0;
  if (fabs(a)<=64*DBL_EPSILON) {
    if (fabs(b)>64*DBL_EPSILON) {
      double t=-c/b;if(t>0.0 && t<1.0)roots[n++]=t;
    }
  }
  else {
    double disc=b*b-4*a*c;
    if (disc>=0.0) {
      double q=-0.5*(b+copysign(sqrt(disc),b));
      double t=q/a;
      if(t>0.0 && t<1.0)roots[n++]=t;
      if(q!=0.0) {
        t=c/q;if(t>0.0 && t<1.0)roots[n++]=t;
      }
    }
  }
  return n;
}

// Exact containment test for the represented p1/p2 boundary curve through
// t=0,1/2,1. This catches an excursion between all three interpolation nodes.
bool
tok_wall_curve_inside(const struct gkyl_efit *e,
  const double p0[2], const double pm[2], const double p1[2])
{
  if (!tok_wall_point_inside(e,p0) || !tok_wall_point_inside(e,pm) ||
      !tok_wall_point_inside(e,p1)) return false;
  double a[2],b[2];
  for(int d=0;d<2;++d) {
    a[d]=2*(p1[d]+p0[d]-2*pm[d]);
    b[d]=4*pm[d]-3*p0[d]-p1[d];
  }
  double *cuts=gkyl_malloc((4*e->limiter_n+2)*sizeof(double));
  int nc=0;cuts[nc++]=0.0;cuts[nc++]=1.0;
  double tol=tok_wall_tolerance(e);
  for(int i=0,j=e->limiter_n-1;i<e->limiter_n;j=i++) {
    double ex=e->limiter_R[i]-e->limiter_R[j],ey=e->limiter_Z[i]-e->limiter_Z[j];
    double l2=ex*ex+ey*ey;
    if(l2==0.0)continue;
    double cx=p0[0]-e->limiter_R[j],cy=p0[1]-e->limiter_Z[j];
    double qa=a[0]*ey-a[1]*ex,qb=b[0]*ey-b[1]*ex,qc=cx*ey-cy*ex;
    double roots[2];
    double threshold=64*DBL_EPSILON*sqrt(l2)*fmax(1.0,hypot(a[0],a[1])+hypot(b[0],b[1])+hypot(cx,cy));
    if(fmax(fabs(qa),fmax(fabs(qb),fabs(qc)))<=threshold) {
      // Curve lies on this wall line. Its passage past either end of the
      // finite wall segment still partitions inside from outside intervals.
      int d=fabs(ex)>=fabs(ey) ? 0 : 1;
      for(int end=0;end<2;++end) {
        int vertex=end ? i : j;
        double value=d==0 ? e->limiter_R[vertex] : e->limiter_Z[vertex];
        int nr=tok_wall_quadratic_roots(a[d],b[d],p0[d]-value,roots);
        for(int k=0;k<nr;++k)cuts[nc++]=roots[k];
      }
    }
    else {
      int nr=tok_wall_quadratic_roots(qa,qb,qc,roots);
      for(int k=0;k<nr;++k) {
        double t=roots[k],r=p0[0]+t*(b[0]+t*a[0]),z=p0[1]+t*(b[1]+t*a[1]);
        double u=((r-e->limiter_R[j])*ex+(z-e->limiter_Z[j])*ey)/l2;
        if(u>=-tol/sqrt(l2) && u<=1.0+tol/sqrt(l2))cuts[nc++]=t;
      }
    }
  }
  qsort(cuts,nc,sizeof(double),tok_wall_compare_double);
  bool ok=true;
  for(int i=1;i<nc && ok;++i) {
    double t=0.5*(cuts[i-1]+cuts[i]);
    double p[2]={p0[0]+t*(b[0]+t*a[0]),p0[1]+t*(b[1]+t*a[1])};
    ok=tok_wall_point_inside(e,p);
  }
  gkyl_free(cuts);return ok;
}

// A join is geometric, not a device-specific distance allowance. Require the
// callback endpoint (and a nearby interior point) to lie on the stored wall.
static bool
tok_wall_contains(const struct gkyl_efit *efit, const double p[2], double tol)
{
  for (int j=0; j<efit->limiter_n; ++j) {
    int next=(j+1)%efit->limiter_n;
    double dr=efit->limiter_R[next]-efit->limiter_R[j];
    double dz=efit->limiter_Z[next]-efit->limiter_Z[j];
    double l2=dr*dr+dz*dz;
    if (l2 == 0.0) continue;
    double t=((p[0]-efit->limiter_R[j])*dr+(p[1]-efit->limiter_Z[j])*dz)/l2;
    t=fmax(0.0, fmin(1.0, t));
    if (hypot(p[0]-efit->limiter_R[j]-t*dr,
        p[1]-efit->limiter_Z[j]-t*dz) <= tol) return true;
  }
  return false;
}

static bool
tok_wall_in_equilibrium(const struct gkyl_tok_geo *geo, const double p[2])
{
  const struct gkyl_rect_grid *g=geo->use_cubics ? &geo->rzgrid_cubic : &geo->rzgrid;
  return p[0]>=g->lower[0] && p[0]<=g->upper[0] &&
    p[1]>=g->lower[1] && p[1]<=g->upper[1];
}

int
tok_divertor_wall_slot(const struct gkyl_tok_geo *geo, plate_func plate)
{
  if (!geo->extend_to_limiter || !plate) return -1;
  int slot=-1;
  if (plate==geo->plate_func_lower && geo->divertor_wall[0].num_segments) slot=0;
  if (plate==geo->plate_func_upper && geo->divertor_wall[1].num_segments)
    slot=slot>=0 ? -2 : 1;
  return slot;
}

// Solve only on the explicitly named connected material target. Requiring a
// unique native root prevents a remote leg/root from replacing an earlier
// contact. The retained contour and all nonmaterial edges are still checked
// against the entire wall after construction.
bool
tok_divertor_wall_intersection(const struct gkyl_tok_geo *geo, int slot,
  double psi, double *r, double *z)
{
  if (slot<0 || slot>1 || !tok_geo_finite(psi)) return false;
  const struct gkyl_tok_geo_wall_target *target=&geo->divertor_wall[slot];
  const struct gkyl_efit *e=geo->efit;
  const struct gkyl_rect_grid *grid=geo->use_cubics ? &geo->rzgrid_cubic : &geo->rzgrid;
  if (!tok_wall_target_valid(e,target->num_segments,target->segments)) return false;
  double step=0.1*fmin(grid->dx[0],grid->dx[1]);
  if (!(step>0.0) || !tok_geo_finite(step)) return false;
  double flux_tol=1e-10*fmax(1.0,fabs(psi)), root_r=0.0,root_z=0.0;
  int roots=0;
  for (int k=0;k<target->num_segments;++k) {
    int i=target->segments[k],j=(i+1)%e->limiter_n;
    double a[2]={e->limiter_R[i],e->limiter_Z[i]};
    double b[2]={e->limiter_R[j],e->limiter_Z[j]};
    if (!tok_wall_in_equilibrium(geo,a) || !tok_wall_in_equilibrium(geo,b)) return false;
    double count=ceil(hypot(b[0]-a[0],b[1]-a[1])/step);
    if (!tok_geo_finite(count) || count>100000) return false;
    int ns=(int)count;
    ns=GKYL_MAX2(ns,32);
    if (ns>100000) return false;
    double f0=tok_eval_psi_rz(geo,a[0],a[1])-psi;
    if (!tok_geo_finite(f0)) return false;
    for (int t=0;t<=ns;++t) {
      double u=t/(double)ns;
      double p[2]={a[0]+u*(b[0]-a[0]),a[1]+u*(b[1]-a[1])};
      double f=tok_eval_psi_rz(geo,p[0],p[1])-psi;
      if (!tok_geo_finite(f)) return false;
      bool exact=fabs(f)<=flux_tol;
      if (exact || (t && f*f0<0.0)) {
        if (!exact) {
          double lo=(t-1)/(double)ns,hi=u,flo=f0;
          for (int iter=0;iter<64;++iter) {
            double mid=0.5*(lo+hi);
            double fm=tok_eval_psi_rz(geo,a[0]+mid*(b[0]-a[0]),a[1]+mid*(b[1]-a[1]))-psi;
            if (!tok_geo_finite(fm)) return false;
            if (flo*fm<=0.0) hi=mid;
            else {lo=mid;flo=fm;}
          }
          u=0.5*(lo+hi);
          p[0]=a[0]+u*(b[0]-a[0]);p[1]=a[1]+u*(b[1]-a[1]);
        }
        if (fabs(tok_eval_psi_rz(geo,p[0],p[1])-psi)>flux_tol) return false;
        if (roots && hypot(p[0]-root_r,p[1]-root_z)>tok_wall_tolerance(e)) {
          fprintf(stderr,"TOK_GEO_DIVERTOR_ROOT_FAILED reason=ambiguous slot=%d psi=%.17g\n",slot,psi);
          return false;
        }
        root_r=p[0];root_z=p[1];++roots;
      }
      f0=f;
    }
  }
  if (!roots) return false;
  *r=root_r;*z=root_z;
  return true;
}

// Follow one outward endpoint continuation in file order or its reverse.
// Stop if flux turns away before reaching the target, or if wall data runs
// out. This deliberately refuses to jump to another leg elsewhere on the wall.
static bool
tok_wall_walk(const struct gkyl_tok_geo *geo, int seg, int dir,
  const double start[2], double psi, double *r, double *z, double *distance)
{
  const struct gkyl_efit *e=geo->efit;
  double p[2]={start[0],start[1]}, travelled=0.0;
  if (!tok_wall_in_equilibrium(geo,p)) return false;
  double f=tok_eval_psi_rz(geo,p[0],p[1])-psi;
  double tol=1e-10*fmax(1.0,fabs(psi));
  const struct gkyl_rect_grid *g=geo->use_cubics ? &geo->rzgrid_cubic : &geo->rzgrid;
  double step=0.05*fmin(g->dx[0],g->dx[1]);
  for (int visited=0; visited<e->limiter_n; ++visited) {
    int v=dir>0 ? (seg+1)%e->limiter_n : seg;
    double end[2]={e->limiter_R[v],e->limiter_Z[v]};
    double base[2]={p[0],p[1]}, length=hypot(end[0]-p[0],end[1]-p[1]);
    int ns=(int)ceil(length/step);
    ns=GKYL_MAX2(ns,32);
    if (ns>100000) return false;
    for (int k=1; k<=ns; ++k) {
      double q[2]={base[0]+(end[0]-base[0])*k/ns,
        base[1]+(end[1]-base[1])*k/ns};
      if (!tok_wall_in_equilibrium(geo,q)) return false;
      double fq=tok_eval_psi_rz(geo,q[0],q[1])-psi;
      if (!tok_geo_finite(fq) || !tok_geo_finite(f)) return false;
      if (f*fq<=0.0 || fabs(fq)<=tol) {
        double lo=0.0, hi=1.0, flo=f;
        if (fabs(fq)>tol) {
          for (int it=0; it<64; ++it) {
            double t=0.5*(lo+hi);
            double fm=tok_eval_psi_rz(geo,p[0]+t*(q[0]-p[0]),p[1]+t*(q[1]-p[1]))-psi;
            if (!tok_geo_finite(fm)) return false;
            if (flo*fm<=0.0) hi=t;
            else {lo=t;flo=fm;}
          }
        }
        double t=fabs(fq)<=tol ? 1.0 : 0.5*(lo+hi);
        *r=p[0]+t*(q[0]-p[0]); *z=p[1]+t*(q[1]-p[1]);
        *distance=travelled+hypot(*r-p[0],*z-p[1]);
        return fabs(tok_eval_psi_rz(geo,*r,*z)-psi)<=tol;
      }
      if (fabs(fq)>fabs(f)+tol) return false;
      travelled+=hypot(q[0]-p[0],q[1]-p[1]);
      p[0]=q[0];p[1]=q[1];f=fq;
    }
    seg=(seg+dir+e->limiter_n)%e->limiter_n;
  }
  return false;
}

bool
tok_limiter_plate_intersection(const struct gkyl_tok_geo *geo,
  plate_func plate, double psi, double *r, double *z)
{
  const struct gkyl_efit *e=geo->efit;
  if (!geo->extend_to_limiter || !plate || e->limiter_status!=1) return false;
  double tol=1e-8*fmax(1.0,fmax(e->rdim,e->zdim));
  // Extend the endpoint through which this flux family leaves the original
  // plate's range. The remote endpoint can otherwise reach the same psi after
  // walking around the vessel onto the other divertor leg.
  double closest_flux=DBL_MAX;
  for (int i=0; i<=512; ++i) {
    double p[2];plate(i/512.0,p);
    if (!tok_wall_in_equilibrium(geo,p)) return false;
    double f=fabs(tok_eval_psi_rz(geo,p[0],p[1])-psi);
    if (!tok_geo_finite(f)) return false;
    closest_flux=fmin(closest_flux,f);
  }
  double flux_tol=1e-10*fmax(1.0,fabs(psi));
  int found=0, picked_end=-1;
  double best_r=0.0,best_z=0.0,best_dist=0.0;
  for (int edge=0; edge<2; ++edge) {
    double p[2],inside[2];
    plate(edge,p);plate(edge ? 1.0-1e-4 : 1e-4,inside);
    if (fabs(tok_eval_psi_rz(geo,p[0],p[1])-psi)>closest_flux+flux_tol) continue;
    if (!tok_wall_contains(e,p,tol) || !tok_wall_contains(e,inside,tol)) continue;
    double tr=p[0]-inside[0],tz=p[1]-inside[1],tn=hypot(tr,tz);
    if (tn<=tol) continue;
    int start_seg=-1,start_dir=0,joins=0;
    for (int j=0; j<e->limiter_n; ++j) {
      int next=(j+1)%e->limiter_n;
      double dr=e->limiter_R[next]-e->limiter_R[j];
      double dz=e->limiter_Z[next]-e->limiter_Z[j],ln=hypot(dr,dz);
      if (ln<=tol) continue;
      double t=((p[0]-e->limiter_R[j])*dr+(p[1]-e->limiter_Z[j])*dz)/(ln*ln);
      if (t < -tol/ln || t > 1.0+tol/ln) continue;
      if (hypot(p[0]-e->limiter_R[j]-t*dr,p[1]-e->limiter_Z[j]-t*dz)>tol) continue;
      double align=(tr*dr+tz*dz)/(tn*ln);
      if (fabs(align)<0.5) continue;
      int dir=align>0.0 ? 1 : -1;
      if ((dir>0 && (1.0-t)*ln<=tol) || (dir<0 && t*ln<=tol)) continue;
      start_seg=j;start_dir=dir;++joins;
    }
    if (joins!=1) continue;
    double rr,zz,dist;
    if (tok_wall_walk(geo,start_seg,start_dir,p,psi,&rr,&zz,&dist)) {
      if (found && hypot(rr-best_r,zz-best_z)>tol) return false;
      best_r=rr;best_z=zz;best_dist=dist;picked_end=edge;++found;
    }
  }
  if (!found) return false;
  *r=best_r;*z=best_z;
  static int nreport=0;
  if (nreport++<8)
    fprintf(stderr,"TOK_GEO_PLATE_EXTENDED psi=%.17g endpoint=%d rz=(%.17g,%.17g) wall_distance=%.17g\n",
      psi,picked_end,*r,*z,best_dist);
  return true;
}

// Find the unique intersection of psi=psi0 with a shaped plate.  The search
// uses only the fixed target surface, so the cached X-point ray is identical
// in corner, interior, and surface geometry passes.
static bool
tok_plate_intersection(const struct gkyl_tok_geo *geo, plate_func plate,
  double psi0, double *rplate, double *zplate)
{
  int slot=tok_divertor_wall_slot(geo,plate);
  if (slot!=-1) return tok_divertor_wall_intersection(geo,slot,psi0,rplate,zplate);
  if (!plate)
    return false;
  const int nsamp = 512;
  const double flux_tol = 1e-10*fmax(1.0, fabs(psi0));
  int nroot = 0;
  double last_root = -DBL_MAX, root_s = 0.0;
  double roots[8];
  double rz0[2], rz1[2];
  plate(0.0, rz0);
  double f0 = tok_eval_psi_rz(geo, rz0[0], rz0[1])-psi0;
  if (isfinite(f0) && fabs(f0) <= flux_tol) {
    root_s = 0.0; last_root = 0.0; roots[0] = 0.0; nroot = 1;
  }
  for (int i=1; i<=nsamp; ++i) {
    double s1 = i/(double) nsamp;
    plate(s1, rz1);
    double f1 = tok_eval_psi_rz(geo, rz1[0], rz1[1])-psi0;
    bool endpoint = isfinite(f1) && fabs(f1) <= flux_tol;
    bool bracket = isfinite(f0) && isfinite(f1) && f0*f1 < 0.0;
    if (endpoint || bracket) {
      double sr = s1;
      if (!endpoint) {
        double lo = (i-1)/(double) nsamp, hi = s1, flo = f0;
        for (int k=0; k<64; ++k) {
          double mid = 0.5*(lo+hi), rzm[2];
          plate(mid, rzm);
          double fm = tok_eval_psi_rz(geo, rzm[0], rzm[1])-psi0;
          if (!isfinite(fm))
            return false;
          if (flo*fm <= 0.0)
            hi = mid;
          else { lo = mid; flo = fm; }
        }
        sr = 0.5*(lo+hi);
      }
      if (fabs(sr-last_root) > 2e-10) {
        root_s = sr; last_root = sr;
        if (nroot < (int)(sizeof(roots)/sizeof(roots[0])))
          roots[nroot] = sr;
        ++nroot;
      }
    }
    rz0[0] = rz1[0]; rz0[1] = rz1[1]; f0 = f1;
  }
  if (nroot > 1) {
    // A plate long enough to reach past the X-point (an inner plate continued
    // onto the floor) is crossed twice by one surface: once by the main-SOL
    // boundary and once by the divertor leg. Rejecting outright would discard
    // the perfectly good leg strike, so select it the same way
    // tok_plate_select_leg_root does -- the crossing beyond the X-point in |Z|
    // that is nearest to it. Falls back to the old rejection when no crossing
    // lies beyond the X-point, and the nroot==1 path below is untouched.
    int nbeyond = 0;
    double best_s = 0.0, best_absz = 0.0;
    int navail = nroot < (int)(sizeof(roots)/sizeof(roots[0]))
      ? nroot : (int)(sizeof(roots)/sizeof(roots[0]));
    for (int k=0; k<navail; ++k) {
      double rzk[2];
      plate(roots[k], rzk);
      double zxpt = rzk[1] < 0.0
        ? (geo->use_cubics ? geo->efit->Zxpt_cubic[0] : geo->efit->Zxpt[0])
        : ((geo->efit->num_xpts > 1 ? (geo->use_cubics ? geo->efit->Zxpt_cubic[1] : geo->efit->Zxpt[1]) : (geo->use_cubics ? geo->efit->Zxpt_cubic[0] : geo->efit->Zxpt[0])));
      if (fabs(rzk[1]) > fabs(zxpt)) {
        if (nbeyond == 0 || fabs(rzk[1]) < best_absz) {
          best_absz = fabs(rzk[1]);
          best_s = roots[k];
        }
        ++nbeyond;
      }
    }
    if (nbeyond == 0)
      return false;
    root_s = best_s;
  }
  else if (nroot != 1)
    return tok_limiter_plate_intersection(geo,plate,psi0,rplate,zplate);

  // The sign-change test cannot tell a crossing from a DISCONTINUITY in the
  // plate. plate_func is not required to be continuous -- two separate divertor
  // tiles are an ordinary thing to describe -- and where two disjoint pieces
  // sit on opposite sides of the surface, f changes sign with no root between
  // them. Bisection then converges onto the gap, and every value it evaluates
  // stays on one piece, so it never notices.
  //
  // Measured before this check existed: a plate split by a gap spanning 40% of
  // its length reported a strike point about 55 mm from the requested surface
  // AND reported success -- a false acceptance, where every other failure mode
  // in this routine refuses. So confirm the point actually lies on the surface
  // before accepting it. A genuine bisected root converges to ~1e-23 here and a
  // root found by the endpoint test satisfies this by construction, so the
  // existing tolerance is the right bar rather than a new tuned one.
  double rz[2];
  plate(root_s, rz);
  double resid = tok_eval_psi_rz(geo, rz[0], rz[1]) - psi0;
  if (!isfinite(resid) || fabs(resid) > flux_tol) {
    fprintf(stderr,
      "TOK_GEO_PLATE_ROOT_REJECTED reason=residual s=%.17g psi_target=%.17g "
      "residual=%.17g tol=%.17g\n", root_s, psi0, resid, flux_tol);
    return false;
  }
  *rplate = rz[0]; *zplate = rz[1];
  return isfinite(*rplate) && isfinite(*zplate);
}

// Only material endpoints actually owned by this block are mandatory. In
// particular, half-domain PF setup also evaluates the other, unused plate.
static bool
tok_plate_slot_required(const struct gkyl_tok_geo *geo,
  enum gkyl_tok_geo_type ftype, double psi, int side)
{
  bool lower=false,upper=false;
  switch (ftype) {
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN:
    case GKYL_GEOMETRY_TOKAMAK_LSN_SOL: lower=upper=true;break;
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO:
    case GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO:
    case GKYL_GEOMETRY_TOKAMAK_PF_LO_R:
    case GKYL_GEOMETRY_TOKAMAK_PF_UP_L: lower=true;break;
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_UP:
    case GKYL_GEOMETRY_TOKAMAK_LSN_SOL_UP:
    case GKYL_GEOMETRY_TOKAMAK_PF_LO_L:
    case GKYL_GEOMETRY_TOKAMAK_PF_UP_R: upper=true;break;
    case GKYL_GEOMETRY_TOKAMAK_IWL:
      lower=upper=(psi-geo->sibry)*(geo->sibry-geo->efit->simag)>=0.0;break;
    default: break;
  }
  return side ? upper : lower;
}

int
tok_plate_coverage_status(const struct gkyl_tok_geo *geo,
  const struct gkyl_tok_geo_grid_inp *inp, double psi)
{
  if (!geo->plate_spec) return 0;
  int status=0;
  for (int side=0; side<2; ++side) {
    if (!tok_plate_slot_required(geo,inp->ftype,psi,side)) continue;
    double r,z;
    if (!tok_plate_intersection(geo,side ? geo->plate_func_upper : geo->plate_func_lower,psi,&r,&z)) {
      fprintf(stderr,"TOK_GEO_ACTIVE_PLATE_FAILED ftype=%d plate=%s psi=%.17g EXTEND_TO_LIMITER=%d limiter_status=%d\n",
        inp->ftype,side ? "upper" : "lower",psi,geo->extend_to_limiter,geo->efit->limiter_status);
      int failure=geo->extend_to_limiter && geo->divertor_wall[side].num_segments ? 1 : 2;
      status=GKYL_MAX2(status,failure);
    }
  }
  return status;
}

bool
gkyl_tok_geo_check_plate_coverage(const struct gkyl_tok_geo *geo,
  const struct gkyl_tok_geo_grid_inp *inp, double psi)
{
  return tok_plate_coverage_status(geo,inp,psi)==0;
}

enum tok_xpt_sector {
  TOK_XPT_CORE,
  TOK_XPT_PF,
  TOK_XPT_SOL_OUT,
  TOK_XPT_SOL_IN,
  TOK_XPT_UNSUPPORTED,
};

static enum tok_xpt_sector
tok_xpt_sector_from_ftype(enum gkyl_tok_geo_type ftype)
{
  switch (ftype) {
    case GKYL_GEOMETRY_TOKAMAK_CORE_R:
    case GKYL_GEOMETRY_TOKAMAK_CORE_L:
      return TOK_XPT_CORE;
    case GKYL_GEOMETRY_TOKAMAK_PF_LO_R:
    case GKYL_GEOMETRY_TOKAMAK_PF_LO_L:
      return TOK_XPT_PF;
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID:
      return TOK_XPT_SOL_OUT;
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO:
      return TOK_XPT_SOL_IN;
    default:
      return TOK_XPT_UNSUPPORTED;
  }
}

static bool
tok_xpt_ray_enabled(const struct gkyl_tok_geo_grid_inp *inp)
{
  enum tok_xpt_sector sector = tok_xpt_sector_from_ftype(inp->ftype);
  if (sector == TOK_XPT_UNSUPPORTED)
    return false;
  if (inp->straight_xpt_ray)
    return true;
  return inp->straight_core_xpt_ray && sector == TOK_XPT_CORE;
}

static bool
tok_nearest_root_at_z(const struct gkyl_tok_geo_grid_inp *inp,
  const struct gkyl_tok_geo *geo,
  enum tok_xpt_sector sector, double psi0, double z, double rx, double zx,
  double *rbest, double *d2best)
{
  double R[4] = { 0 }, dRdZ[4] = { 0 }, dR[4] = { 0 }, dZ[4] = { 0 };
  int nr = gkyl_tok_geo_R_psiZ(geo, psi0, z, 4, R, dRdZ, dR, dZ);
  bool found = false;
  *d2best = DBL_MAX;

  // An R threshold around the X point lets the inboard and outboard tests
  // overlap whenever the two roots approach the saddle.  Select the root
  // belonging to the requested SOL branch by its established branch seed
  // instead.  Core and PF searches remain global because their nearest point
  // can lie on either side of a turning point.
  int nlo = 0, nhi = nr;
  int selected = -1;
  if (sector == TOK_XPT_SOL_OUT || sector == TOK_XPT_SOL_IN) {
    double reference = sector == TOK_XPT_SOL_OUT ? inp->rright : inp->rleft;
    double best_reference_distance = DBL_MAX;
    for (int n=0; n<nr; ++n) {
      double distance = fabs(R[n]-reference);
      if (distance < best_reference_distance) {
        best_reference_distance = distance;
        selected = n;
      }
    }
    nlo = selected;
    nhi = selected+1;
  }
  for (int n=0; n<nr; ++n) {
    if (n < nlo || n >= nhi)
      continue;
    double d2 = SQ(R[n]-rx)+SQ(z-zx);
    if (d2 < *d2best) {
      *d2best = d2;
      *rbest = R[n];
      found = true;
    }
  }
  return found;
}

static bool
tok_nearest_point_on_surface(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *arc_ctx, enum tok_xpt_sector sector, double psi0,
  double rx, double zx, double *rnear, double *znear)
{
  const struct gkyl_tok_geo *geo = arc_ctx->geo;
  double zlo, zhi;
  if (sector == TOK_XPT_PF) {
    // These input bounds are fixed for the block.  Do not use the current
    // surface's plate intersections here: corner, interior, and surface
    // calculations begin at different radial quadrature points, but must all
    // initialize exactly the same ray on xpt_ray_psi0.
    zlo = fmax(geo->rzgrid.lower[1],
      fmin(inp->zmin_left, inp->zmin_right));
    zhi = zx;
  }
  else if (sector == TOK_XPT_SOL_OUT || sector == TOK_XPT_SOL_IN) {
    // A lower-SOL half surface runs from its physical divertor plate to the
    // magnetic-axis midplane.  The nearest point to the X point is allowed to
    // lie below the X point; restricting this search to X->midplane selected a
    // visibly non-nearest ray in shot 203730.
    double rplate = 0.0, zplate = 0.0;
    if (!tok_plate_intersection(geo, geo->plate_func_lower, psi0,
        &rplate, &zplate))
      return false;
    zlo = fmin(zplate, geo->zmaxis);
    zhi = fmax(zplate, geo->zmaxis);
  }
  else {
    zlo = fmin(zx, geo->zmaxis);
    zhi = fmax(zx, geo->zmaxis);
  }
  if (!(zhi > zlo))
    return false;

  // This search is intentionally global in Z. It is run only once per block,
  // and avoids making the selected ray depend on a branch-local Newton guess.
  const int nsamp = 4*geo->rzgrid.cells[1]+1;
  double best_d2 = DBL_MAX, best_r = 0.0, best_z = 0.0;
  int best_i = -1;
  for (int i=0; i<nsamp; ++i) {
    double z = zlo+(zhi-zlo)*i/(nsamp-1.0), r, d2;
    if (tok_nearest_root_at_z(inp, geo, sector, psi0, z, rx, zx, &r, &d2)
        && d2 < best_d2) {
      best_d2 = d2; best_r = r; best_z = z; best_i = i;
    }
  }
  if (best_i < 0)
    return false;

  // Refine the sampled minimum without differentiating through DG-cell
  // interfaces or assuming that only one R root exists at a given Z.
  double dzs = (zhi-zlo)/(nsamp-1.0);
  double a = fmax(zlo, best_z-dzs), b = fmin(zhi, best_z+dzs);
  const double gr = 0.6180339887498948482;
  for (int iter=0; iter<64; ++iter) {
    double ztrial[2] = { b-gr*(b-a), a+gr*(b-a) };
    double dtrial[2] = { DBL_MAX, DBL_MAX }, rtrial[2] = { best_r, best_r };
    for (int k=0; k<2; ++k)
      tok_nearest_root_at_z(inp, geo, sector, psi0, ztrial[k], rx, zx,
        &rtrial[k], &dtrial[k]);
    if (dtrial[0] < best_d2) {
      best_d2 = dtrial[0]; best_r = rtrial[0]; best_z = ztrial[0];
    }
    if (dtrial[1] < best_d2) {
      best_d2 = dtrial[1]; best_r = rtrial[1]; best_z = ztrial[1];
    }
    if (dtrial[0] <= dtrial[1])
      b = ztrial[1];
    else
      a = ztrial[0];
  }

  *rnear = best_r;
  *znear = best_z;
  return isfinite(best_r) && isfinite(best_z);
}

static bool
tok_xpt_branch_diag_enabled(void)
{
  const char *d = getenv("GKYL_TOK_XPT_BRANCH_DIAG");
  return d && d[0] != '\0' && d[0] != '0';
}

bool
tok_xpt_classify_branch_at_point(const struct gkyl_tok_geo_grid_inp *inp,
  const struct arc_length_ctx *arc_ctx, double Rpoint, double Zpoint,
  double psi, bool *resolved, bool *on_right)
{
  double R[4] = { 0 }, dRdZ[4] = { 0 }, dR[4] = { 0 }, dZ[4] = { 0 };
  int nr = gkyl_tok_geo_R_psiZ(arc_ctx->geo, psi, Zpoint,
    4, R, dRdZ, dR, dZ);
  *resolved = false;
  if (nr < 2) {
    if (tok_xpt_branch_diag_enabled())
      fprintf(stderr,
        "TOK_XPT_BRANCH_DIAG outcome=unresolved_too_few_roots ftype=%d nr=%d "
        "point=(%.17g,%.17g) psi=%.17g\n",
        inp->ftype, nr, Rpoint, Zpoint, psi);
    return true;
  }

  double rr = choose_closest(inp->rright, R, R, nr);
  double rl = choose_closest(inp->rleft, R, R, nr);
  double scale = fmax(1.0, fmax(fabs(rr), fabs(rl)));
  if (fabs(rr-rl) <= 1e-8*scale) {
    if (tok_xpt_branch_diag_enabled())
      fprintf(stderr,
        "TOK_XPT_BRANCH_DIAG outcome=unresolved_coalesced ftype=%d nr=%d "
        "point=(%.17g,%.17g) psi=%.17g rr=%.17g rl=%.17g\n",
        inp->ftype, nr, Rpoint, Zpoint, psi, rr, rl);
    return true;
  }

  double error_right = fabs(Rpoint-rr), error_left = fabs(Rpoint-rl);
  *on_right = error_right <= error_left;
  *resolved = true;
  if (tok_xpt_branch_diag_enabled()) {
    // margin is the fraction by which the winning side wins; a value near 0 means
    // the probe sits nearly equidistant from both candidate roots, so the label is
    // a coin flip rather than a real branch identity.
    double denom = fmax(error_right+error_left, 1e-300);
    double margin = fabs(error_right-error_left)/denom;
    fprintf(stderr,
      "TOK_XPT_BRANCH_DIAG outcome=resolved ftype=%d nr=%d point=(%.17g,%.17g) "
      "psi=%.17g roots=[%.17g,%.17g,%.17g,%.17g] rr=%.17g rl=%.17g "
      "err_r=%.17g err_l=%.17g margin=%.6e on_right=%d\n",
      inp->ftype, nr, Rpoint, Zpoint, psi,
      R[0], R[1], nr > 2 ? R[2] : 0.0/0.0, nr > 3 ? R[3] : 0.0/0.0,
      rr, rl, error_right, error_left, margin, (int) *on_right);
  }
  return true;
}

static bool
tok_xpt_initialize_branch(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *arc_ctx, enum tok_xpt_sector sector,
  double rx, double zx)
{
  if (sector == TOK_XPT_SOL_OUT || sector == TOK_XPT_SOL_IN) {
    arc_ctx->xpt_ray_on_right = sector == TOK_XPT_SOL_OUT;
    arc_ctx->xpt_ray_branch_valid = true;
    return true;
  }

  // A nearest core/PF endpoint is commonly a turning point, where the two R
  // roots coalesce and cannot label the branch.  Probe the already-fixed ray
  // away from both endpoints, using psi evaluated at the probe itself.  This
  // makes branch identity independent of which corner/Gauss surface happened
  // to initialize this arc-length context first.
  // Each probe evaluates psi at its own location, so the probes sample different
  // flux surfaces, and the classifier labels a point by whether it is the larger or
  // smaller of the two R roots at that Z.  That rank is not invariant along a
  // straight ray: the contour legs are curved, so their ordering relative to the ray
  // can swap with no change in the ray's actual branch identity.  Treating a single
  // disagreement as fatal therefore rejects otherwise sound geometry.  Take the
  // majority over all probes instead; on a tie, trust the probe farthest from the
  // X point, where the two roots are best separated.  Shots whose probes already
  // agree unanimously are unaffected.
  const double probes[] = { 0.5, 0.25, 0.75, 0.125, 0.875 };
  int n_right = 0, n_left = 0;
  double smax_right = -1.0, smax_left = -1.0;
  for (unsigned i=0; i<sizeof(probes)/sizeof(probes[0]); ++i) {
    double s = probes[i];
    double R = rx+s*(arc_ctx->xpt_ray_r0-rx);
    double Z = zx+s*(arc_ctx->xpt_ray_z0-zx);
    double psi = tok_eval_psi_rz(arc_ctx->geo, R, Z);
    bool resolved = false, on_right = false;
    if (!tok_xpt_classify_branch_at_point(inp, arc_ctx, R, Z, psi,
        &resolved, &on_right))
      return false;
    if (!resolved)
      continue;
    if (on_right) { ++n_right; smax_right = fmax(smax_right, s); }
    else { ++n_left; smax_left = fmax(smax_left, s); }
  }
  if (n_right == 0 && n_left == 0) {
    fprintf(stderr,
      "TOK_XPT_RAY_DIAG reason=unresolved_ray_branch ftype=%d endpoint=(%.17g,%.17g)\n",
      inp->ftype, arc_ctx->xpt_ray_r0, arc_ctx->xpt_ray_z0);
    return false;
  }
  bool selected_right = n_right != n_left ? n_right > n_left
    : smax_right > smax_left;
  if (n_right > 0 && n_left > 0)
    fprintf(stderr,
      "TOK_XPT_RAY_DIAG reason=ray_branch_probe_disagreement ftype=%d "
      "n_right=%d n_left=%d smax_right=%.17g smax_left=%.17g selected_right=%d\n",
      inp->ftype, n_right, n_left, smax_right, smax_left, (int) selected_right);
  arc_ctx->xpt_ray_on_right = selected_right;
  arc_ctx->xpt_ray_branch_valid = true;
  return true;
}

// Walk a straight ray leaving the X point at a fixed angle and report where it
// first reaches the far flux surface, together with how far outside the flux
// band [psisep, psi0] it strayed on the way.  q is the normalized flux
// coordinate, 0 on the separatrix and 1 on the far surface, so it is sign
// convention independent: some equilibria have psi increasing outward and
// others decreasing.
static bool
tok_xpt_ray_walk(const struct gkyl_tok_geo *geo, double rx, double zx,
  double angle, double smax, double psisep, double delta,
  double *s_cross, double *qmin)
{
  const int nsamp = 2048;
  double ca = cos(angle), sa = sin(angle);
  double q_prev = 0.0, s_prev = 0.0;
  *qmin = 0.0; *s_cross = 0.0;
  for (int i=1; i<=nsamp; ++i) {
    double s = smax*i/nsamp;
    double q = (tok_eval_psi_rz(geo, rx+s*ca, zx+s*sa)-psisep)/delta;
    if (!isfinite(q))
      return false;
    *qmin = fmin(*qmin, q);
    if (q >= 1.0) {
      // Bisect the bracketing interval for the crossing itself.
      double lo = s_prev, hi = s;
      for (int k=0; k<60; ++k) {
        double mid = 0.5*(lo+hi);
        double qm = (tok_eval_psi_rz(geo, rx+mid*ca, zx+mid*sa)-psisep)/delta;
        if (!isfinite(qm))
          return false;
        if (qm >= 1.0) hi = mid; else lo = mid;
      }
      *s_cross = 0.5*(lo+hi);
      return *s_cross > 0.0;
    }
    q_prev = q; s_prev = s;
  }
  (void) q_prev;
  return false;
}

// Replace an X-point ray that leaves its own flux band.  The anchor is normally
// the nearest point on the far surface, which points the ray along the saddle's
// principal axis -- correct where the SOL channel is straight.  Where the
// channel bends sharply away from that axis the straight ray crosses a
// separatrix branch, runs outside the block's own flux region and re-enters,
// and the seam it defines is no longer a boundary of the block.  Sweep the ray
// angle and take the nearest anchor whose entire ray stays in the band.
static bool
tok_xpt_ray_reanchor_in_band(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *arc_ctx, double rx, double zx, double delta,
  double band_tol)
{
  const struct gkyl_tok_geo *geo = arc_ctx->geo;
  double angle0 = atan2(arc_ctx->xpt_ray_z0-zx, arc_ctx->xpt_ray_r0-rx);
  double s0 = hypot(arc_ctx->xpt_ray_r0-rx, arc_ctx->xpt_ray_z0-zx);
  if (!(s0 > 0.0) || !isfinite(angle0))
    return false;
  const double sweep = 60.0*M_PI/180.0;
  const int nangle = 961;
  // Accept only rays that stay in the band to within DG round-off, not merely
  // to within the trigger threshold: an anchor taken from the edge of the valid
  // angular range grazes the separatrix branch and folds the corner cell just
  // as the original did.  The trigger stays loose so that rays with a harmless
  // wiggle are left exactly as they are.
  const double accept_tol = 1.0e-6*fmax(1.0, band_tol/1.0e-3);
  double best_s = DBL_MAX, best_angle = 0.0, best_qmin = 0.0;
  for (int i=0; i<nangle; ++i) {
    double angle = angle0+sweep*(2.0*i/(nangle-1.0)-1.0);
    double s_cross = 0.0, qmin = 0.0;
    if (!tok_xpt_ray_walk(geo, rx, zx, angle, 4.0*s0, geo->psisep, delta,
        &s_cross, &qmin))
      continue;
    if (qmin < -accept_tol || s_cross >= best_s)
      continue;
    best_s = s_cross; best_angle = angle; best_qmin = qmin;
  }
  if (best_s == DBL_MAX)
    return false;
  double r_new = rx+best_s*cos(best_angle);
  double z_new = zx+best_s*sin(best_angle);
  fprintf(stderr,
    "TOK_XPT_RAY_DIAG reason=reanchored_in_band ftype=%d old_angle=%.6f "
    "new_angle=%.6f old_s=%.17g new_s=%.17g new_qmin=%.6e "
    "old_endpoint=(%.17g,%.17g) new_endpoint=(%.17g,%.17g)\n",
    inp->ftype, angle0*180.0/M_PI, best_angle*180.0/M_PI, s0, best_s,
    best_qmin, arc_ctx->xpt_ray_r0, arc_ctx->xpt_ray_z0, r_new, z_new);
  arc_ctx->xpt_ray_r0 = r_new;
  arc_ctx->xpt_ray_z0 = z_new;
  return true;
}

static bool
tok_xpt_ray_anchor(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *arc_ctx, double psi_curr)
{
  const struct gkyl_tok_geo *geo = arc_ctx->geo;
  enum tok_xpt_sector sector = tok_xpt_sector_from_ftype(inp->ftype);
  arc_ctx->xpt_anchor_valid = false;
  double rx = geo->use_cubics ? geo->efit->Rxpt_cubic[0] : geo->efit->Rxpt[0];
  double zx = geo->use_cubics ? geo->efit->Zxpt_cubic[0] : geo->efit->Zxpt[0];

  if (!arc_ctx->xpt_ray_initialized) {
    if (!tok_nearest_point_on_surface(inp, arc_ctx, sector,
        arc_ctx->xpt_ray_psi0, rx, zx,
        &arc_ctx->xpt_ray_r0, &arc_ctx->xpt_ray_z0)) {
      fprintf(stderr, "TOK_XPT_RAY_DIAG reason=no_far_surface_point ftype=%d target_psi=%.17g\n",
        inp->ftype, arc_ctx->xpt_ray_psi0);
      return false;
    }

    double psi_endpoint = tok_eval_psi_rz(geo, arc_ctx->xpt_ray_r0,
      arc_ctx->xpt_ray_z0);
    double flux_scale = fmax(1.0,
      fmax(fabs(arc_ctx->xpt_ray_psi0), fabs(geo->psisep)));
    if (!isfinite(psi_endpoint) ||
        fabs(psi_endpoint-arc_ctx->xpt_ray_psi0) > 1e-9*flux_scale) {
      fprintf(stderr,
        "TOK_XPT_RAY_DIAG reason=far_surface_flux_residual ftype=%d residual=%.17g requested_psi=%.17g evaluated_psi=%.17g endpoint=(%.17g,%.17g)\n",
        inp->ftype, psi_endpoint-arc_ctx->xpt_ray_psi0,
        arc_ctx->xpt_ray_psi0, psi_endpoint, arc_ctx->xpt_ray_r0,
        arc_ctx->xpt_ray_z0);
      return false;
    }
    double delta = arc_ctx->xpt_ray_psi0-geo->psisep;
    if (!isfinite(delta) || fabs(delta) <= 256.0*DBL_EPSILON*flux_scale) {
      fprintf(stderr, "TOK_XPT_RAY_DIAG reason=degenerate_flux_span ftype=%d psi_x=%.17g psi_far=%.17g\n",
        inp->ftype, geo->psisep, arc_ctx->xpt_ray_psi0);
      return false;
    }

    // The quadratic DG representation can have a small local reversal at a
    // cell boundary even when the underlying equilibrium is smooth.  A global
    // monotonicity test therefore rejects otherwise unambiguous rays.  Check
    // finiteness here; uniqueness is tested below for every actual geometry
    // and finite-difference target before an anchor is accepted.
    double ray_qmin = 0.0;
    for (int i=1; i<=256; ++i) {
      double s = i/256.0;
      double R = rx+s*(arc_ctx->xpt_ray_r0-rx);
      double Z = zx+s*(arc_ctx->xpt_ray_z0-zx);
      double q = (tok_eval_psi_rz(geo, R, Z)-geo->psisep)/delta;
      if (!isfinite(q)) {
        fprintf(stderr,
          "TOK_XPT_RAY_DIAG reason=nonfinite_initial_ray ftype=%d sample=%d endpoint=(%.17g,%.17g)\n",
          inp->ftype, i, arc_ctx->xpt_ray_r0, arc_ctx->xpt_ray_z0);
        return false;
      }
      ray_qmin = fmin(ray_qmin, q);
    }
    // A ray that dips out of its own flux band has crossed a separatrix branch
    // and no longer bounds this block: the seam it defines cuts across the
    // neighbouring region, and the last interior theta column overtakes it,
    // folding the corner cells.  Measured over 450 NSTX-U equilibria, the
    // outboard-SOL blocks that fold are exactly the three with
    // ray_qmin <= -6.3e-3, while every non-folding one stays above -7.9e-4 --
    // an order of magnitude of separation either side of this threshold.
    const char *tol_env = getenv("GKYL_TOK_XPT_RAY_BAND_TOL");
    const double ray_band_tol =
      tol_env && tol_env[0] != '\0' ? atof(tol_env) : 3.0e-3;
    if (ray_qmin < -ray_band_tol &&
        !tok_xpt_ray_reanchor_in_band(inp, arc_ctx, rx, zx, delta,
          ray_band_tol))
      fprintf(stderr,
        "TOK_XPT_RAY_DIAG reason=no_in_band_reanchor ftype=%d ray_qmin=%.6e "
        "endpoint=(%.17g,%.17g)\n",
        inp->ftype, ray_qmin, arc_ctx->xpt_ray_r0, arc_ctx->xpt_ray_z0);
    if (!tok_xpt_initialize_branch(inp, arc_ctx, sector, rx, zx))
      return false;
    arc_ctx->xpt_ray_initialized = true;
  }

  double psi1 = arc_ctx->xpt_ray_psi0;
  double delta = psi1-geo->psisep;
  double scale = fmax(1.0, fmax(fabs(psi_curr), fabs(geo->psisep)));
  double tol = 128.0*DBL_EPSILON*scale;
  if (fabs(psi_curr-geo->psisep) <= tol) {
    arc_ctx->xpt_anchor_r = rx;
    arc_ctx->xpt_anchor_z = zx;
    arc_ctx->xpt_anchor_valid = true;
    return true;
  }
  if (fabs(psi_curr-psi1) <= tol) {
    arc_ctx->xpt_anchor_r = arc_ctx->xpt_ray_r0;
    arc_ctx->xpt_anchor_z = arc_ctx->xpt_ray_z0;
    arc_ctx->xpt_anchor_valid = true;
    return true;
  }

  double dr = arc_ctx->xpt_ray_r0-rx, dz = arc_ctx->xpt_ray_z0-zx;
  double smax = DBL_MAX;
  if (dr > 0.0)
    smax = fmin(smax, (geo->rzgrid.upper[0]-rx)/dr);
  else if (dr < 0.0)
    smax = fmin(smax, (geo->rzgrid.lower[0]-rx)/dr);
  if (dz > 0.0)
    smax = fmin(smax, (geo->rzgrid.upper[1]-zx)/dz);
  else if (dz < 0.0)
    smax = fmin(smax, (geo->rzgrid.lower[1]-zx)/dz);
  smax *= 1.0-64.0*DBL_EPSILON;
  if (!isfinite(smax) || smax <= 0.0) {
    fprintf(stderr, "TOK_XPT_RAY_DIAG reason=invalid_ray_extent ftype=%d smax=%.17g\n",
      inp->ftype, smax);
    return false;
  }

  double qtarget = (psi_curr-geo->psisep)/delta;
  if (qtarget < -1e-8 || qtarget > 1.0+1e-8) {
    fprintf(stderr,
      "TOK_XPT_RAY_DIAG reason=target_outside_ray_span ftype=%d qtarget=%.17g psi=%.17g ray_psi=(%.17g,%.17g)\n",
      inp->ftype, qtarget, psi_curr, geo->psisep, psi1);
    return false;
  }
  double search_max = fmin(1.0, smax);
  double qprev = (tok_eval_psi_rz(geo, rx, zx)-geo->psisep)/delta;
  double qpeak = qprev;
  double sprev = 0.0, slo = -1.0, sup = -1.0;
  // Track the last upward crossing and the deepest dip after the first one, so a
  // spurious local reversal in the DG representation (all crossings clustered in a
  // tiny s window) can be told apart from a genuinely folded ray.
  double slast_lo = -1.0, slast_up = -1.0, qmin_after_first = DBL_MAX;
  // Lowest q seen before the ray first reaches the target surface. If the ray
  // never dips inside the separatrix on that stretch, the first crossing is
  // joined to the X point through the block's own flux band.
  double qmin_before_first = 0.0;
  int upward_crossings = 0, downward_crossings = 0;
  double crossing_tol = 1e-13*fmax(1.0, fabs(qtarget));
  const int nsamp = 512;
  for (int i=1; i<=nsamp; ++i) {
    double s = search_max*i/(double) nsamp;
    double R = rx+s*dr, Z = zx+s*dz;
    double q = (tok_eval_psi_rz(geo, R, Z)-geo->psisep)/delta;
    if (!isfinite(q)) {
      fprintf(stderr,
        "TOK_XPT_RAY_DIAG reason=nonfinite_anchor_search ftype=%d sample=%d qtarget=%.17g\n",
        inp->ftype, i, qtarget);
      return false;
    }
    bool upward = qprev < qtarget-crossing_tol &&
      q >= qtarget-crossing_tol;
    bool downward = qprev > qtarget+crossing_tol &&
      q <= qtarget+crossing_tol;
    if (upward) {
      upward_crossings++;
      if (slo < 0.0) {
        slo = sprev;
        sup = s;
      }
      slast_lo = sprev;
      slast_up = s;
    }
    if (downward)
      downward_crossings++;
    if (upward_crossings == 0)
      qmin_before_first = fmin(qmin_before_first, q);
    if (upward_crossings > 0)
      qmin_after_first = fmin(qmin_after_first, q);
    qprev = q;
    qpeak = fmax(qpeak, q);
    sprev = s;
  }
  if (slo < 0.0) {
    fprintf(stderr,
      "TOK_XPT_RAY_DIAG reason=no_anchor_bracket ftype=%d qtarget=%.17g qpeak=%.17g qend=%.17g smax=%.17g\n",
      inp->ftype, qtarget, qpeak, qprev, search_max);
    return false;
  }
  // A ray that folds back can cross psi=psi_curr more than once, and the
  // residual check accepts either crossing, so this used to bail out. But the
  // crossings are not interchangeable: the anchor has to vary continuously with
  // psi_curr, because psi=psisep anchors exactly at the X point above and the
  // near-separatrix samples here are the finite-difference stencil for the
  // radial derivative there. Only the FIRST crossing tends to the X point as
  // psi_curr -> psisep; a later one stays a finite distance away and would put a
  // step in the derivative.
  //
  // Accept the first crossing when the ray reaches it without ever dipping
  // inside the separatrix, i.e. it is joined to the X point through the block's
  // own flux band. A ray that wandered out of that band before arriving is still
  // rejected. upward_crossings==1 never enters here, so shots that already
  // anchor unambiguously are untouched.
  // The X point is a saddle, so psi is stationary there and a sub-millimetre
  // error in its location -- the cubic and quadratic reps routinely disagree by
  // ~1e-3 m -- puts the first samples on the wrong side of psisep by a
  // vanishing amount. Judging "did the ray leave the flux band" at 1e-13 would
  // call that noise a genuine excursion, so measure the dip against the flux
  // offset actually being resolved. Measured separation is wide: 204502 dips to
  // 0.22% of qtarget (4.8e-8 in psi, 6e-7 of psisep), whereas a ray that really
  // crosses into the confined region dips by 1800% of qtarget.
  double band_tol = fmax(crossing_tol, 0.05*fabs(qtarget));
  if (upward_crossings > 1 && qmin_before_first >= -band_tol) {
    fprintf(stderr,
      "TOK_XPT_RAY_DIAG reason=first_crossing_anchor ftype=%d qtarget=%.17g "
      "upward_crossings=%d bracket=(%.17g,%.17g) qmin_before_first=%.6e\n",
      inp->ftype, qtarget, upward_crossings, slo, sup, qmin_before_first);
  }
  else if (upward_crossings != 1) {
    double crossing_span = slast_up-slo;
    double dip_below = qtarget-qmin_after_first;
    // How high q rises on the near-X-point lobe before it dips back decides how
    // many radial nodes are ambiguous: every node with qtarget < qlobe_max sees
    // two crossings, every node above it sees only the far one.
    {
      double qlobe_max = -DBL_MAX, s_lobe = -1.0, qp = 0.0;
      for (int i=1; i<=nsamp; ++i) {
        double s = search_max*i/(double) nsamp;
        double q = (tok_eval_psi_rz(geo, rx+s*dr, zx+s*dz)-geo->psisep)/delta;
        if (i > 1 && q < qp) break;   // first turnover ends the lobe
        if (q > qlobe_max) { qlobe_max = q; s_lobe = s; }
        qp = q;
      }
      fprintf(stderr,
        "TOK_XPT_RAY_LOBE ftype=%d qtarget=%.17g qlobe_max=%.17g s_lobe=%.17g "
        "psi_curr=%.17g psisep=%.17g psi_far=%.17g xpt=(%.17g,%.17g) far=(%.17g,%.17g) "
        "qmin_before_first=%.6e (dpsi=%.6e) qmin_after_first=%.6e\n",
        inp->ftype, qtarget, qlobe_max, s_lobe, psi_curr, geo->psisep, psi1,
        rx, zx, arc_ctx->xpt_ray_r0, arc_ctx->xpt_ray_z0,
        qmin_before_first, qmin_before_first*delta, qmin_after_first);
    }
    fprintf(stderr,
      "TOK_XPT_RAY_DIAG reason=ambiguous_anchor_crossing ftype=%d qtarget=%.17g upward_crossings=%d downward_crossings=%d first_bracket=(%.17g,%.17g) last_bracket=(%.17g,%.17g) crossing_span=%.17g dip_below_target=%.6e search_max=%.17g\n",
      inp->ftype, qtarget, upward_crossings, downward_crossings, slo, sup,
      slast_lo, slast_up, crossing_span, dip_below, search_max);
    return false;
  }

  for (int n=0; n<80; ++n) {
    double smid = 0.5*(slo+sup);
    double R = rx+smid*dr, Z = zx+smid*dz;
    double qmid = (tok_eval_psi_rz(geo, R, Z)-geo->psisep)/delta;
    if (!isfinite(qmid)) {
      fprintf(stderr,
        "TOK_XPT_RAY_DIAG reason=nonfinite_anchor_refinement ftype=%d qtarget=%.17g s=%.17g\n",
        inp->ftype, qtarget, smid);
      return false;
    }
    if (qmid-qtarget >= 0.0)
      sup = smid;
    else
      slo = smid;
  }
  double s = 0.5*(slo+sup);
  arc_ctx->xpt_anchor_r = rx+s*dr;
  arc_ctx->xpt_anchor_z = zx+s*dz;
  double residual = tok_eval_psi_rz(geo, arc_ctx->xpt_anchor_r,
    arc_ctx->xpt_anchor_z)-psi_curr;
  if (!isfinite(residual) || fabs(residual) > 1e-9*scale) {
    fprintf(stderr, "TOK_XPT_RAY_DIAG reason=anchor_flux_residual ftype=%d residual=%.17g anchor=(%.17g,%.17g)\n",
      inp->ftype, residual, arc_ctx->xpt_anchor_r,
      arc_ctx->xpt_anchor_z);
    return false;
  }
  arc_ctx->xpt_anchor_valid = true;
  return true;
}

static bool
tok_xpt_anchor_on_right(const struct gkyl_tok_geo_grid_inp *inp,
  const struct arc_length_ctx *arc_ctx)
{
  if (!arc_ctx->xpt_ray_branch_valid) {
    fprintf(stderr,
      "TOK_XPT_RAY_DIAG reason=missing_cached_branch ftype=%d\n",
      inp->ftype);
    abort();
  }
  bool resolved = false, on_right = false;
  if (!tok_xpt_classify_branch_at_point(inp, arc_ctx,
      arc_ctx->xpt_anchor_r, arc_ctx->xpt_anchor_z, arc_ctx->psi,
      &resolved, &on_right))
    abort();
  // The cached branch is a majority decision over probes spanning the ray (see
  // tok_xpt_initialize_branch).  A single anchor sample disagreeing with it carries
  // the same rank-swap ambiguity as any individual probe, so it is not grounds to
  // abort the whole build; the cached value is the better-supported answer and is
  // what this function returns either way.
  if (resolved && on_right != arc_ctx->xpt_ray_on_right)
    fprintf(stderr,
      "TOK_XPT_RAY_DIAG reason=anchor_branch_disagrees_with_cached ftype=%d psi=%.17g anchor=(%.17g,%.17g) cached_right=%d resolved_right=%d action=keep_cached\n",
      inp->ftype, arc_ctx->psi, arc_ctx->xpt_anchor_r,
      arc_ctx->xpt_anchor_z, arc_ctx->xpt_ray_on_right, on_right);
  return arc_ctx->xpt_ray_on_right;
}

static void
tok_xpt_map_interval(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *arc_ctx, double lo, double hi)
{
  double dtheta = inp->cgrid.upper[2]-inp->cgrid.lower[2];
  if (!isfinite(lo) || !isfinite(hi) || !(hi > lo) || !(dtheta > 0.0)) {
    fprintf(stderr,
      "TOK_XPT_RAY invalid arc interval ftype=%d psi=%.17g lo=%.17g hi=%.17g dtheta=%.17g\n",
      inp->ftype, arc_ctx->psi, lo, hi, dtheta);
    abort();
  }
  arc_ctx->xpt_map_arc_lo = lo;
  arc_ctx->xpt_map_arc_hi = hi;
  arc_ctx->xpt_map_darc_dtheta = (hi-lo)/dtheta;
  arc_ctx->xpt_map_valid = true;
}

static void
tok_configure_xpt_map(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *arc_ctx)
{
  arc_ctx->xpt_map_valid = false;
  arc_ctx->xpt_anchor_valid = false;
  if (!tok_xpt_ray_enabled(inp))
    return;
  // Blocks on the extended construction get their poloidal origin from the
  // topology-aware boundary traces built in tok_geo.c, which can carry two
  // distinct X-point rays or a closed core seam. The scalar lower-X-point
  // anchor and arc-interval map below are specific to the chord construction,
  // so they are simply not built for those blocks (matching the early return
  // in tok_prepare_ordered_map).
  if (tok_ext_construction(inp))
    return;
  if (!inp->half_domain) {
    fprintf(stderr,
      "TOK_XPT_RAY ftype=%d currently requires half_domain=true\n", inp->ftype);
    abort();
  }
  if (!tok_xpt_ray_anchor(inp, arc_ctx, arc_ctx->psi)) {
    fprintf(stderr,
      "TOK_XPT_RAY failed ftype=%d psi=%.17g target_psi=%.17g\n",
      inp->ftype, arc_ctx->psi, arc_ctx->xpt_ray_psi0);
    abort();
  }

  enum tok_xpt_sector sector = tok_xpt_sector_from_ftype(inp->ftype);
  double cut = 0.0, lo = 0.0, hi = 0.0;
  bool at_sep = fabs(arc_ctx->psi-arc_ctx->geo->psisep)
    <= 128.0*DBL_EPSILON*fmax(1.0, fabs(arc_ctx->geo->psisep));

  if (sector == TOK_XPT_CORE) {
    if (at_sep)
      cut = 0.0;
    else if (tok_xpt_anchor_on_right(inp, arc_ctx))
      cut = integrate_psi_contour_memo(arc_ctx->geo, arc_ctx->psi,
        arc_ctx->zmin, arc_ctx->xpt_anchor_z, inp->rright,
        true, false, arc_ctx->arc_memo_right);
    else
      cut = arc_ctx->arcL_right
        + integrate_psi_contour_memo(arc_ctx->geo, arc_ctx->psi,
          arc_ctx->xpt_anchor_z, arc_ctx->zmax, inp->rleft,
          true, false, arc_ctx->arc_memo_left);
    double out_mid = integrate_psi_contour_memo(arc_ctx->geo, arc_ctx->psi,
      arc_ctx->zmin, arc_ctx->geo->zmaxis, inp->rright,
      true, false, arc_ctx->arc_memo_right);
    double in_mid = arc_ctx->arcL_right
      + integrate_psi_contour_memo(arc_ctx->geo, arc_ctx->psi,
        arc_ctx->geo->zmaxis, arc_ctx->zmax, inp->rleft,
        true, false, arc_ctx->arc_memo_left);
    if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R) {
      lo = cut; hi = out_mid;
    }
    else {
      lo = in_mid; hi = cut;
    }
    while (hi <= lo)
      hi += arc_ctx->arcL_tot;
  }
  else if (sector == TOK_XPT_PF) {
    if (at_sep)
      cut = arc_ctx->arcL_right;
    else if (tok_xpt_anchor_on_right(inp, arc_ctx))
      cut = integrate_psi_contour_memo(arc_ctx->geo, arc_ctx->psi,
        arc_ctx->zmin_right, arc_ctx->xpt_anchor_z, inp->rright,
        true, false, arc_ctx->arc_memo_right);
    else
      cut = arc_ctx->arcL_right
        + integrate_psi_contour_memo(arc_ctx->geo, arc_ctx->psi,
          arc_ctx->xpt_anchor_z, arc_ctx->zmax, inp->rleft,
          true, false, arc_ctx->arc_memo_left);
    if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_R) {
      lo = 0.0; hi = cut;
    }
    else {
      lo = cut; hi = arc_ctx->arcL_tot;
    }
  }
  else if (sector == TOK_XPT_SOL_OUT) {
    cut = at_sep
      ? integrate_psi_contour_memo(arc_ctx->geo, arc_ctx->psi,
          arc_ctx->zmin, arc_ctx->xpt_anchor_z, inp->rright,
          true, false, arc_ctx->arc_memo)
      : integrate_psi_contour_memo(arc_ctx->geo, arc_ctx->psi,
          arc_ctx->zmin, arc_ctx->xpt_anchor_z, inp->rright,
          true, false, arc_ctx->arc_memo);
    double mid = integrate_psi_contour_memo(arc_ctx->geo, arc_ctx->psi,
      arc_ctx->zmin, arc_ctx->geo->zmaxis, inp->rright,
      true, false, arc_ctx->arc_memo);
    if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO) {
      lo = 0.0; hi = cut;
    }
    else {
      lo = cut; hi = mid;
    }
  }
  else if (sector == TOK_XPT_SOL_IN) {
    cut = integrate_psi_contour_memo(arc_ctx->geo, arc_ctx->psi,
      arc_ctx->xpt_anchor_z, arc_ctx->zmax, inp->rleft,
      true, false, arc_ctx->arc_memo);
    double mid = integrate_psi_contour_memo(arc_ctx->geo, arc_ctx->psi,
      arc_ctx->geo->zmaxis, arc_ctx->zmax, inp->rleft,
      true, false, arc_ctx->arc_memo);
    if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID) {
      lo = mid; hi = cut;
    }
    else {
      lo = cut; hi = arc_ctx->arcL_tot;
    }
  }
  tok_xpt_map_interval(inp, arc_ctx, lo, hi);
}

void
tok_prepare_ordered_map(struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *arc_ctx, double psi_curr)
{
  arc_ctx->psi = psi_curr;
  arc_ctx->xpt_map_valid = false;
  arc_ctx->xpt_anchor_valid = false;
  // Preserve the old coordinate algebra in callers while the ordered map
  // consumes theta directly; no physical contour length is needed here.
  arc_ctx->arcL_tot = 2.0*M_PI;
  // Full-domain single- and double-null blocks use topology-aware boundary
  // traces built lazily in tok_geo.c.  They can have two distinct X-point
  // rays (or a closed core seam), so the scalar lower-X-point anchor used by
  // the original half-domain path is intentionally bypassed here.
  if (tok_ext_construction(inp))
    return;
  if (!tok_xpt_ray_enabled(inp) ||
      !tok_xpt_ray_anchor(inp, arc_ctx, psi_curr)) {
    fprintf(stderr,
      "TOK_ORDERED_MAP failed to initialize ray ftype=%d psi=%.17g target_psi=%.17g\n",
      inp->ftype, psi_curr, arc_ctx->xpt_ray_psi0);
    abort();
  }
}


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




static bool
tok_extent_diag_enabled(void)
{
  static int cached = -1;
  if (cached < 0) {
    const char *env = getenv("GKYL_TOK_EXTENT_DIAG");
    cached = (env && env[0] != '\0' && env[0] != '0') ? 1 : 0;
  }
  return cached == 1;
}

// The DN block decomposition assumes the plate-derived span brackets both
// X-points: zmin < zxpt_lo < zxpt_up < zmax. When it does not, the divertor-leg
// arc integrals below collapse to zero and the corresponding _LO/_UP block is
// handed a zero-width theta range, which only surfaces much later as a
// divide-by-zero in tok_ordered_map_lookup. Report it here, where the cause is
// still visible. Diagnostic only -- no behaviour change when the span is sane.
static void
tok_check_xpt_span(enum gkyl_tok_geo_type ftype, double zmin, double zmax,
  double zxpt_lo, double zxpt_up, double arcL_tot, double arcL_lo,
  double arcL_up)
{
  if (zmin < zxpt_lo && zxpt_up < zmax)
    return;
  fprintf(stderr,
    "TOK_EXTENT reason=xpt_outside_plate_span ftype=%d zmin=%.17g zxpt_lo=%.17g "
    "zxpt_up=%.17g zmax=%.17g margin_lo=%.17g margin_up=%.17g "
    "arcL_tot=%.17g arcL_lo=%.17g arcL_up=%.17g\n",
    ftype, zmin, zxpt_lo, zxpt_up, zmax, zxpt_lo-zmin, zmax-zxpt_up,
    arcL_tot, arcL_lo, arcL_up);
}

// Diagnostic only: enumerate every sign change of psi(plate(s))-psi_curr over
// s in [0,1], so a plate that the flux surface crosses more than once (where
// gkyl_ridders returns just one of the roots) is visible.
static void
tok_plate_scan_diag(struct gkyl_tok_geo *geo, struct plate_ctx *pctx,
  const char *side, int ftype, double psi_curr)
{
  enum { NSAMP = 4001 };
  double rzplate[2];
  double s_prev = 0.0, f_prev = tok_plate_psi_func(0.0, pctx);
  int ncross = 0;
  for (int i=1; i<NSAMP; ++i) {
    double s = (double) i/(NSAMP-1);
    double f = tok_plate_psi_func(s, pctx);
    if (f_prev == 0.0 || (f_prev < 0.0) != (f < 0.0)) {
      // Bisect the bracket to locate the crossing well enough to report R,Z.
      double a = s_prev, b = s, fa = f_prev;
      for (int it=0; it<60; ++it) {
        double m = 0.5*(a+b), fm = tok_plate_psi_func(m, pctx);
        if ((fa < 0.0) != (fm < 0.0)) { b = m; } else { a = m; fa = fm; }
      }
      double sroot = 0.5*(a+b);
      if (pctx->lower) geo->plate_func_lower(sroot, rzplate);
      else geo->plate_func_upper(sroot, rzplate);
      fprintf(stderr,
        "TOK_PLATE_SCAN side=%s ftype=%d idx=%d s=%.17g rz=(%.17g,%.17g)\n",
        side, ftype, ncross, sroot, rzplate[0], rzplate[1]);
      ncross++;
    }
    s_prev = s; f_prev = f;
  }
  fprintf(stderr, "TOK_PLATE_SCAN side=%s ftype=%d ncross=%d psi=%.17g\n",
    side, ftype, ncross, psi_curr);
}

// A divertor plate that reaches past the X-point in |Z| (e.g. an inner plate
// continued onto the floor) can be crossed by one flux surface more than once:
// once by the main-SOL boundary and once by the divertor leg itself. gkyl_ridders
// is only given the whole [0,1] bracket and returns an arbitrary one of them.
//
// Select deterministically instead: among crossings that lie *beyond* the
// X-point in |Z| -- the only ones that can be a leg strike -- take the one
// nearest the X-point, i.e. the first material surface the leg reaches.
//
// This is a no-op unless the plate really is crossed more than once, and even
// then it can only ever prefer a crossing at smaller |Z|, so adding plate
// length further out (larger |Z|) cannot change a shot that already had a valid
// crossing. Returns true and writes *s_out when it selects a root; false leaves
// the caller's existing gkyl_ridders answer in place.
// Note the X-point to compare against is chosen from each candidate root's own
// Z sign, not from the caller's upper/lower naming: PF_LO_L legitimately passes
// the *lower* inner plate as its plate_func_upper.
static bool
tok_plate_select_leg_root(struct gkyl_tok_geo *geo, struct plate_ctx *pctx,
  double s_ridders, double *s_out)
{
  double zxpt_lo = geo->use_cubics ? geo->efit->Zxpt_cubic[0] : geo->efit->Zxpt[0];
  double zxpt_up = (geo->efit->num_xpts > 1 ? (geo->use_cubics ? geo->efit->Zxpt_cubic[1] : geo->efit->Zxpt[1]) : (geo->use_cubics ? geo->efit->Zxpt_cubic[0] : geo->efit->Zxpt[0]));
  enum { NSAMP = 801 };
  double rzplate[2];
  double s_prev = 0.0, f_prev = tok_plate_psi_func(0.0, pctx);
  double best_s = 0.0, best_absz = 0.0;
  int nroots = 0, nbeyond = 0;
  for (int i=1; i<NSAMP; ++i) {
    double s = (double) i/(NSAMP-1);
    double f = tok_plate_psi_func(s, pctx);
    if (f_prev == 0.0 || (f_prev < 0.0) != (f < 0.0)) {
      double a = s_prev, b = s, fa = f_prev;
      for (int it=0; it<60; ++it) {
        double m = 0.5*(a+b), fm = tok_plate_psi_func(m, pctx);
        if ((fa < 0.0) != (fm < 0.0)) { b = m; } else { a = m; fa = fm; }
      }
      double sroot = 0.5*(a+b);
      nroots++;
      if (pctx->lower) geo->plate_func_lower(sroot, rzplate);
      else geo->plate_func_upper(sroot, rzplate);
      // "Beyond the X-point" means further from the midplane on this root's own
      // side of the machine.
      double zxpt = rzplate[1] < 0.0 ? zxpt_lo : zxpt_up;
      if (fabs(rzplate[1]) > fabs(zxpt)) {
        if (nbeyond == 0 || fabs(rzplate[1]) < best_absz) {
          best_absz = fabs(rzplate[1]);
          best_s = sroot;
        }
        nbeyond++;
      }
    }
    s_prev = s; f_prev = f;
  }
  if (nroots <= 1 || nbeyond == 0)
    return false;
  // Only override when we actually land on a *different* crossing than
  // gkyl_ridders found. Where the two agree, returning false keeps the original
  // value bit-for-bit, so this cannot perturb zmin/zmax by root-finder rounding
  // on the shots that were already correct.
  if (pctx->lower) geo->plate_func_lower(s_ridders, rzplate);
  else geo->plate_func_upper(s_ridders, rzplate);
  if (fabs(fabs(rzplate[1]) - best_absz) < 1e-8)
    return false;
  *s_out = best_s;
  return true;
}

// Sets zmax if plate is specified
// gkyl_ridders reports an invalid bracket rather than iterating on one, and
// every plate solve below took its result on faith.  status != 0 means
// psi(plate(s)) never reaches the requested surface's psi for any s in [0,1]
// -- the divertor plate does not span that surface, which is exactly what
// asking the SOL to reach past the plate looks like.  The reported root is
// DBL_MAX there, and it was fed straight into the plate parameterization: on
// a straight plate Z_lo + (Z_up-Z_lo)*s that returns DBL_MAX scaled by the
// plate's extent, poisoning zmin/zmax.  The run then died far downstream in a
// wild psi fetch, with nothing said about the geometry that caused it.
static double
tok_plate_root_s(const struct gkyl_qr_res *res, bool lower, double psi)
{
  // NOT keyed on status, and not on isfinite either: gkyl_ridders reports
  // status=2 both for an invalid bracket and for a solve that stopped before
  // meeting eps, and the latter still returns a usable root these callers have
  // always relied on.  Its unusable value is the DBL_MAX the search starts
  // from -- a finite double, so isfinite() does not see it.  Every plate solve
  // here brackets s on [0,1], so a usable root lies in it and the sentinel
  // does not.
  if (res->res >= 0.0 && res->res <= 1.0)
    return res->res;
  // The plate does not span this surface.  CLAMP to the nearest plate end
  // rather than abort: the meaningful answer is the furthest the plate
  // actually reaches, and the caller consumes that endpoint's R or Z.
  //
  // DO NOT go back to aborting here.  Measured 2026-08-26 on the 450-shot
  // NSTX-U suite: aborting failed 12 shots at every resolution
  // (203963_ms457, 204046_ms433, 204051_ms450, 204076_ms193, 204077_ms187,
  // 204080_ms193, 204112_ms907, 204118_ms553, 204170_ms580, 205017_ms205,
  // 205050_ms219, 205088_ms195), all with plate=upper, status=2,
  // res=DBL_MAX at psisep from the PF_LO extent setup.  Every one of those
  // shots builds a CLEAN 8-block grid when the root is clamped, and did so
  // historically too.  What the abort caught was a value the half-domain
  // path does not depend on.  The clamp still removes the DBL_MAX that
  // poisoned TCV: it never reaches the plate parameterization.
  //
  // GKYL_TOK_PLATE_ROOT_STRICT=1 restores the hard abort for diagnosis.
  double clamped = res->res > 1.0 || !(res->res == res->res) ? 1.0
    : (res->res < 0.0 ? 0.0 : res->res);
  static int nwarn = 0;
  const char *strict = getenv("GKYL_TOK_PLATE_ROOT_STRICT");
  if (strict && strict[0] != '\0' && strict[0] != '0') {
    fprintf(stderr,
      "TOK_GEO_PLATE_ROOT_FAILED plate=%s psi=%.17g status=%d res=%.17g\n"
      "  psi(plate(s)) does not reach this surface for any s in [0,1].\n",
      lower ? "lower" : "upper", psi, res->status, res->res);
    abort();
  }
  if (nwarn < 8) {
    ++nwarn;
    fprintf(stderr,
      "TOK_GEO_PLATE_ROOT_CLAMPED plate=%s psi=%.17g status=%d res=%.17g -> s=%g\n"
      "  the divertor plate does not span this flux surface; using the\n"
      "  nearest plate end.%s\n",
      lower ? "lower" : "upper", psi, res->status, res->res, clamped,
      nwarn == 8 ? "  (further occurrences suppressed)" : "");
  }
  return clamped;
}

static void
tok_legacy_plate_point(struct gkyl_tok_geo *geo, bool lower, enum gkyl_tok_geo_type ftype, double psi,
  const struct gkyl_qr_res *res, double *s, double rz[2])
{
  plate_func plate=lower ? geo->plate_func_lower : geo->plate_func_upper;
  int slot=lower ? 0 : 1;
  if (geo->extend_to_limiter && geo->divertor_wall[slot].num_segments) {
    if (tok_divertor_wall_intersection(geo,slot,psi,&rz[0],&rz[1])) {
      *s=-1.0;
      return;
    }
    if (tok_plate_slot_required(geo,ftype,psi,slot)) {
      fprintf(stderr,"TOK_GEO_DIVERTOR_ROOT_FAILED slot=%d psi=%.17g\n",slot,psi);
      abort();
    }
    // Inactive half-domain endpoints retain the existing auxiliary fallback;
    // this point is never granted ownership of a material face.

  }
  if (res->res>=0.0 && res->res<=1.0) {
    *s=res->res;plate(*s,rz);return;
  }
  if (tok_plate_intersection(geo,plate,psi,&rz[0],&rz[1])) {
    *s=-1.0; // resolved R,Z may be on the wall, outside the callback's domain
    return;
  }
  // This can only serve an auxiliary endpoint: active endpoints are checked
  // before contour construction. Keep the historical half-domain clamp here.
  *s=tok_plate_root_s(res,lower,psi);plate(*s,rz);
}

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
      double smax;
      tok_legacy_plate_point(geo, false, arc_ctx->ftype, psi_curr, &res, &smax, rzplate);
      double s_leg;
      if (smax>=0.0 && smax<=1.0 && tok_plate_select_leg_root(geo, pctx, smax, &s_leg)) {
        smax = s_leg;
        geo->plate_func_upper(smax, rzplate);
      }
      arc_ctx->zmax = rzplate[1];
      if (tok_extent_diag_enabled()) {
        fprintf(stderr,
          "TOK_EXTENT_PLATE side=upper ftype=%d psi=%.17g fa=%.17g fb=%.17g "
          "bracketed=%d s=%.17g status=%d rz=(%.17g,%.17g)\n",
          arc_ctx->ftype, psi_curr, fa, fb, (fa*fb <= 0.0), smax, res.status,
          rzplate[0], rzplate[1]);
      }
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
      double smin;
      tok_legacy_plate_point(geo, true, arc_ctx->ftype, psi_curr, &res, &smin, rzplate);
      double s_leg;
      if (smin>=0.0 && smin<=1.0 && tok_plate_select_leg_root(geo, pctx, smin, &s_leg)) {
        smin = s_leg;
        geo->plate_func_lower(smin, rzplate);
      }
      arc_ctx->zmin = rzplate[1];
      if (tok_extent_diag_enabled()) {
        fprintf(stderr,
          "TOK_EXTENT_PLATE side=lower ftype=%d psi=%.17g fa=%.17g fb=%.17g "
          "bracketed=%d s=%.17g status=%d rz=(%.17g,%.17g)\n",
          arc_ctx->ftype, psi_curr, fa, fb, (fa*fb <= 0.0), smin, res.status,
          rzplate[0], rzplate[1]);
        tok_plate_scan_diag(geo, pctx, "lower", arc_ctx->ftype, psi_curr);
      }
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
      double smax;
      tok_legacy_plate_point(geo, false, arc_ctx->ftype, psi_curr, &res, &smax, rzplate);
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
      double smin;
      tok_legacy_plate_point(geo, true, arc_ctx->ftype, psi_curr, &res, &smin, rzplate);
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
    double zxpt_up = (geo->efit->num_xpts > 1 ? (geo->use_cubics ? geo->efit->Zxpt_cubic[1] : geo->efit->Zxpt[1]) : (geo->use_cubics ? geo->efit->Zxpt_cubic[0] : geo->efit->Zxpt[0]));
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
    tok_check_xpt_span(inp->ftype, arc_ctx.zmin, arc_ctx.zmax, zxpt_lo, zxpt_up,
      arcL_tot, arcL_lo, arcL_up);
    if (tok_extent_diag_enabled()) {
      fprintf(stderr,
        "TOK_EXTENT_ARCS ftype=%d rclose=%.17g zmin=%.17g zmax=%.17g "
        "zxpt_lo=%.17g zxpt_up=%.17g arcL_tot=%.17g arcL_lo=%.17g "
        "arcL_mid=%.17g arcL_up=%.17g theta_lo=%.17g theta_up=%.17g\n",
        inp->ftype, arc_ctx.rclose, arc_ctx.zmin, arc_ctx.zmax,
        zxpt_lo, zxpt_up, arcL_tot, arcL_lo, arcL_mid, arcL_up,
        *theta_lo, *theta_up);
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
    double zxpt_up = (geo->efit->num_xpts > 1 ? (geo->use_cubics ? geo->efit->Zxpt_cubic[1] : geo->efit->Zxpt[1]) : (geo->use_cubics ? geo->efit->Zxpt_cubic[0] : geo->efit->Zxpt[0]));
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
    tok_check_xpt_span(inp->ftype, arc_ctx.zmin, arc_ctx.zmax, zxpt_lo, zxpt_up,
      arcL_tot, arcL_lo, arcL_up);
    if (tok_extent_diag_enabled()) {
      fprintf(stderr,
        "TOK_EXTENT_ARCS ftype=%d rclose=%.17g zmin=%.17g zmax=%.17g "
        "zxpt_lo=%.17g zxpt_up=%.17g arcL_tot=%.17g arcL_lo=%.17g "
        "arcL_mid=%.17g arcL_up=%.17g theta_lo=%.17g theta_up=%.17g\n",
        inp->ftype, arc_ctx.rclose, arc_ctx.zmin, arc_ctx.zmax,
        zxpt_lo, zxpt_up, arcL_tot, arcL_lo, arcL_mid, arcL_up,
        *theta_lo, *theta_up);
    }
  }
  else if(inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE || inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R || inp->ftype ==  GKYL_GEOMETRY_TOKAMAK_CORE_L){
    // Immediately set rleft and rright. Will need both
    arc_ctx.rright = inp->rright;
    arc_ctx.rleft = inp->rleft;

    double zxpt_up = (geo->efit->num_xpts > 1 ? (geo->use_cubics ? geo->efit->Zxpt_cubic[1] : geo->efit->Zxpt[1]) : (geo->use_cubics ? geo->efit->Zxpt_cubic[0] : geo->efit->Zxpt[0]));
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
      double smax;
      if (geo->extend_to_limiter && geo->divertor_wall[1].num_segments)
        tok_legacy_plate_point(geo,false,inp->ftype,geo->psisep,&res,&smax,rzplate);
      else {
        smax=tok_plate_root_s(&res,false,geo->psisep);
        geo->plate_func_upper(smax,rzplate);
      }
      arc_ctx.zmin_left = rzplate[1];

      pctx.lower=true;
      a = 0;
      b = 1;
      fa = tok_plate_psi_func(a, &pctx);
      fb = tok_plate_psi_func(b, &pctx);
      res = gkyl_ridders(tok_plate_psi_func, &pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smin;
      if (geo->extend_to_limiter && geo->divertor_wall[0].num_segments)
        tok_legacy_plate_point(geo,true,inp->ftype,geo->psisep,&res,&smin,rzplate);
      else {
        smin=tok_plate_root_s(&res,true,geo->psisep);
        geo->plate_func_lower(smin,rzplate);
      }
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
      double theta_pi = inp->straight_xpt_ray && !inp->half_domain
        ? M_PI-del : M_PI+del;
      *theta_up = theta_pi - arcL_up/arcL_tot*2.0*M_PI;
    }
    else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_UP) {
      double theta_pi = inp->straight_xpt_ray && !inp->half_domain
        ? M_PI-del : M_PI+del;
      *theta_lo = theta_pi - arcL_up/arcL_tot*2.0*M_PI;
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
      double smax;
      if (geo->extend_to_limiter && geo->divertor_wall[1].num_segments)
        tok_legacy_plate_point(geo,false,inp->ftype,geo->psisep,&res,&smax,rzplate);
      else {
        smax=tok_plate_root_s(&res,false,geo->psisep);
        geo->plate_func_upper(smax,rzplate);
      }
      arc_ctx.zmin_left = rzplate[1];


      pctx.lower=true;
      a = 0;
      b = 1;
      fa = tok_plate_psi_func(a, &pctx);
      fb = tok_plate_psi_func(b, &pctx);
      res = gkyl_ridders(tok_plate_psi_func, &pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smin;
      if (geo->extend_to_limiter && geo->divertor_wall[0].num_segments)
        tok_legacy_plate_point(geo,true,inp->ftype,geo->psisep,&res,&smin,rzplate);
      else {
        smin=tok_plate_root_s(&res,true,geo->psisep);
        geo->plate_func_lower(smin,rzplate);
      }
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
    double zxpt_up = (geo->efit->num_xpts > 1 ? (geo->use_cubics ? geo->efit->Zxpt_cubic[1] : geo->efit->Zxpt[1]) : (geo->use_cubics ? geo->efit->Zxpt_cubic[0] : geo->efit->Zxpt[0]));
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
      double smax;
      if (geo->extend_to_limiter && geo->divertor_wall[1].num_segments)
        tok_legacy_plate_point(geo,false,inp->ftype,geo->psisep,&res,&smax,rzplate);
      else {
        smax=tok_plate_root_s(&res,false,geo->psisep);
        geo->plate_func_upper(smax,rzplate);
      }
      arc_ctx.zmax_right= rzplate[1];

      pctx.lower=true;
      a = 0;
      b = 1;
      fa = tok_plate_psi_func(a, &pctx);
      fb = tok_plate_psi_func(b, &pctx);
      res = gkyl_ridders(tok_plate_psi_func, &pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smin;
      if (geo->extend_to_limiter && geo->divertor_wall[0].num_segments)
        tok_legacy_plate_point(geo,true,inp->ftype,geo->psisep,&res,&smin,rzplate);
      else {
        smin=tok_plate_root_s(&res,true,geo->psisep);
        geo->plate_func_lower(smin,rzplate);
      }
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


  if (!gkyl_tok_geo_check_plate_coverage(geo,inp,psi_curr))
    abort();

  // Set psicurr no matter what
  arc_ctx->psi = psi_curr;
  arc_ctx->xpt_map_valid = false;
  arc_ctx->xpt_anchor_valid = false;

  if(inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE || inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R || inp->ftype ==  GKYL_GEOMETRY_TOKAMAK_CORE_L){
    // Immediately set rleft and rright. Will need both
    arc_ctx->rright = inp->rright;
    arc_ctx->rleft = inp->rleft;

    double zxpt_up = (geo->efit->num_xpts > 1 ? (geo->use_cubics ? geo->efit->Zxpt_cubic[1] : geo->efit->Zxpt[1]) : (geo->use_cubics ? geo->efit->Zxpt_cubic[0] : geo->efit->Zxpt[0]));
    double zxpt_lo = geo->use_cubics ? geo->efit->Zxpt_cubic[0] : geo->efit->Zxpt[0];
    arc_ctx->zmax = inp->zmax ? inp->zmax : zxpt_up; // Initial guess.
                                                  // zmax is specified for single null full core
    double zlo = geo->zmaxis;
    find_upper_turning_point(geo, psi_curr, zlo, &arc_ctx->zmax, 0);
    arc_ctx->zmin = zxpt_lo; // Initial guess
    double zup = geo->zmaxis;
    find_lower_turning_point(geo, psi_curr, zup, &arc_ctx->zmin, 0);
    if (fabs(psi_curr-geo->psisep) <=
        128.0*DBL_EPSILON*fmax(1.0, fabs(geo->psisep)))
      arc_ctx->zmin = zxpt_lo;
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
    if (fabs(psi_curr-geo->psisep) <=
        128.0*DBL_EPSILON*fmax(1.0, fabs(geo->psisep)))
      arc_ctx->zmax = zxpt_lo;

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
      double smax;
      tok_legacy_plate_point(geo, false, arc_ctx->ftype, psi_curr, &res, &smax, rzplate);
      arc_ctx->zmin_left = rzplate[1];

      pctx->lower=true;
      a = 0;
      b = 1;
      fa = tok_plate_psi_func(a, pctx);
      fb = tok_plate_psi_func(b, pctx);
      res = gkyl_ridders(tok_plate_psi_func, pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smin;
      tok_legacy_plate_point(geo, true, arc_ctx->ftype, psi_curr, &res, &smin, rzplate);
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
    double zxpt_up = (geo->efit->num_xpts > 1 ? (geo->use_cubics ? geo->efit->Zxpt_cubic[1] : geo->efit->Zxpt[1]) : (geo->use_cubics ? geo->efit->Zxpt_cubic[0] : geo->efit->Zxpt[0]));
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
      double smax;
      tok_legacy_plate_point(geo, false, arc_ctx->ftype, psi_curr, &res, &smax, rzplate);
      arc_ctx->zmax_right= rzplate[1];

      pctx->lower=true;
      a = 0;
      b = 1;
      fa = tok_plate_psi_func(a, pctx);
      fb = tok_plate_psi_func(b, pctx);
      res = gkyl_ridders(tok_plate_psi_func, pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smin;
      tok_legacy_plate_point(geo, true, arc_ctx->ftype, psi_curr, &res, &smin, rzplate);
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
    const char *bounded_turn = getenv("GKYL_TOK_LEGACY_BOUNDED_TURN");
    if (bounded_turn && bounded_turn[0] && bounded_turn[0] != '0') {
      double rturn;
      if (!tok_ext_turning_point(inp, geo, psi_curr, true, &rturn, &arc_ctx->zmax)) {
        fprintf(stderr, "TOK_LEGACY_BOUNDED_TURN failed ftype=%d psi=%.17g\n", inp->ftype, psi_curr);
        abort();
      }
      if (getenv("GKYL_TOK_LEGACY_TURN_DIAG"))
        fprintf(stderr, "TOK_LEGACY_BOUNDED_TURN ftype=%d psi=%.17g guess=%.17g turn=(%.17g,%.17g)\n",
          inp->ftype, psi_curr, inp->zmax, rturn, arc_ctx->zmax);
    }
    else
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
      double smax;
      tok_legacy_plate_point(geo, false, arc_ctx->ftype, psi_curr, &res, &smax, rzplate);
      arc_ctx->zmin_left = rzplate[1];

      pctx->lower=true;
      a = 0;
      b = 1;
      fa = tok_plate_psi_func(a, pctx);
      fb = tok_plate_psi_func(b, pctx);
      res = gkyl_ridders(tok_plate_psi_func, pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smin;
      tok_legacy_plate_point(geo, true, arc_ctx->ftype, psi_curr, &res, &smin, rzplate);
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
    if (nr_lcfs <= 0) {
      fprintf(stderr,
        "TOK_GEOMETRY IWL failed to find LCFS root at psi=%.17g Z=%.17g\n",
        geo->efit->sibry, geo->efit->zmaxis);
      abort();
    }
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

  tok_configure_xpt_map(inp, arc_ctx);
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
    if (arc_ctx->xpt_map_valid) {
      while (arcL_curr < 0.0)
        arcL_curr += arc_ctx->arcL_tot;
      while (arcL_curr >= arc_ctx->arcL_tot)
        arcL_curr -= arc_ctx->arcL_tot;
      arc_ctx->pre = false;
      if (arcL_curr <= arc_ctx->arcL_right) {
        *rclose = arc_ctx->rright;
        arc_ctx->right = true;
        *ridders_min = -arcL_curr;
        *ridders_max = arc_ctx->arcL_tot-arcL_curr;
      }
      else {
        *rclose = arc_ctx->rleft;
        arc_ctx->right = false;
        *ridders_min = arc_ctx->arcL_tot-arcL_curr;
        *ridders_max = arc_ctx->arcL_right-arcL_curr;
      }
    }
    else if(arcL_curr <= arc_ctx->arcL_start){
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
