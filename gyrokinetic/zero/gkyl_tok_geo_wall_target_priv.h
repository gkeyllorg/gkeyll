#pragma once

#include <gkyl_efit.h>
#include <gkyl_tok_geo.h>   // gkyl_tok_wall_usable: one shared usability predicate

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

// Private predicates for an explicitly supplied material target. Segment i
// always denotes limiter vertex i -> (i+1)%limiter_n. The list may traverse
// those segments in either direction; it must describe one simple open arc.
// No predicate here grants an exemption to an interior point or flux contour.

static inline bool
tok_wall_target_finite(double x)
{
  // This project uses -ffast-math, which may remove isfinite() checks.
  uint64_t bits;
  memcpy(&bits, &x, sizeof bits);
  return (bits & UINT64_C(0x7ff0000000000000)) !=
    UINT64_C(0x7ff0000000000000);
}

static inline double
tok_wall_target_tolerance(const struct gkyl_efit *e)
{
  return 1e-8*fmax(1.0, fmax(e->rdim, e->zdim));
}

static inline bool
tok_wall_target_point_on_segment(const double p[2], const double a[2],
  const double b[2], double tol)
{
  if (!p || !tok_wall_target_finite(p[0]) || !tok_wall_target_finite(p[1]))
    return false;
  double dx=b[0]-a[0], dy=b[1]-a[1], length=hypot(dx,dy);
  double px=p[0]-a[0], py=p[1]-a[1];
  if (!tok_wall_target_finite(length) || !(length>tol) ||
      !tok_wall_target_finite(px) || !tok_wall_target_finite(py)) return false;
  double ux=dx/length, uy=dy/length, along=px*ux+py*uy;
  if (!tok_wall_target_finite(along)) return false;
  along=fmax(0.0,fmin(length,along));
  double distance=hypot(px-along*ux,py-along*uy);
  return tok_wall_target_finite(distance) && distance<=tol;
}

// Touching non-neighbor segments, including coincident endpoints, prevents a
// unique simple subarc. Numerical overflow is also an invalid geometry here.
static inline bool
tok_wall_target_segments_touch(const double a[2], const double b[2],
  const double c[2], const double d[2], double tol)
{
  if (tok_wall_target_point_on_segment(a,c,d,tol) ||
      tok_wall_target_point_on_segment(b,c,d,tol) ||
      tok_wall_target_point_on_segment(c,a,b,tol) ||
      tok_wall_target_point_on_segment(d,a,b,tol)) return true;
  double la=hypot(b[0]-a[0],b[1]-a[1]);
  double lb=hypot(d[0]-c[0],d[1]-c[1]);
  double ux=(b[0]-a[0])/la, uy=(b[1]-a[1])/la;
  double vx=(d[0]-c[0])/lb, vy=(d[1]-c[1])/lb;
  double ac=ux*(c[1]-a[1])-uy*(c[0]-a[0]);
  double ad=ux*(d[1]-a[1])-uy*(d[0]-a[0]);
  double ca=vx*(a[1]-c[1])-vy*(a[0]-c[0]);
  double cb=vx*(b[1]-c[1])-vy*(b[0]-c[0]);
  if (!tok_wall_target_finite(ac) || !tok_wall_target_finite(ad) ||
      !tok_wall_target_finite(ca) || !tok_wall_target_finite(cb)) return true;
  return ((ac<0.0 && ad>0.0) || (ac>0.0 && ad<0.0)) &&
    ((ca<0.0 && cb>0.0) || (ca>0.0 && cb<0.0));
}

static inline bool
tok_wall_target_valid(const struct gkyl_efit *e, int n, const int *segments)
{
  if (!gkyl_tok_wall_usable(e) ||
      !e->limiter_R || !e->limiter_Z || !segments || n<1 || n>=e->limiter_n ||
      !tok_wall_target_finite(e->rdim) || !tok_wall_target_finite(e->zdim) ||
      e->rdim<0.0 || e->zdim<0.0) return false;
  const int m=e->limiter_n;
  const double tol=tok_wall_target_tolerance(e);
  for (int i=0; i<m; ++i)
    if (!tok_wall_target_finite(e->limiter_R[i]) ||
        !tok_wall_target_finite(e->limiter_Z[i])) return false;

  int direction=0;
  double total_length=0.0;
  for (int i=0; i<n; ++i) {
    int s=segments[i];
    if (s<0 || s>=m) return false;
    for (int j=0; j<i; ++j) if (segments[j]==s) return false;
    if (i>0) {
      int prev=segments[i-1];
      int step=s==(prev+1)%m ? 1 : s==(prev==0 ? m-1 : prev-1) ? -1 : 0;
      if (!step || (direction && direction!=step)) return false;
      direction=step;
    }
    int t=(s+1)%m;
    double length=hypot(e->limiter_R[t]-e->limiter_R[s],
      e->limiter_Z[t]-e->limiter_Z[s]);
    if (!tok_wall_target_finite(length) || !(length>tol)) return false;
    total_length+=length;
    if (!tok_wall_target_finite(total_length)) return false;
  }

  for (int i=0; i<n; ++i) for (int j=i+1; j<n; ++j) {
    int ia=segments[i], ib=(ia+1)%m, ic=segments[j], id=(ic+1)%m;
    double a[2]={e->limiter_R[ia],e->limiter_Z[ia]};
    double b[2]={e->limiter_R[ib],e->limiter_Z[ib]};
    double c[2]={e->limiter_R[ic],e->limiter_Z[ic]};
    double d[2]={e->limiter_R[id],e->limiter_Z[id]};
    if (j==i+1) {
      // Adjacent edges may share their declared vertex, but may not double
      // back and overlap. Ordinary collinear forward continuation is valid.
      const double *outer_a=ib==ic ? a : b;
      const double *outer_b=ib==ic ? d : c;
      if (tok_wall_target_point_on_segment(outer_a,c,d,tol) ||
          tok_wall_target_point_on_segment(outer_b,a,b,tol)) return false;
    }
    else if (tok_wall_target_segments_touch(a,b,c,d,tol)) return false;
  }
  return true;
}

static inline bool
tok_wall_target_coordinate_validated(const struct gkyl_efit *e, int n,
  const int *segments, const double p[2], double *coordinate)
{
  if (!p || !coordinate || !tok_wall_target_finite(p[0]) ||
      !tok_wall_target_finite(p[1]))
    return false;
  double tol=tok_wall_target_tolerance(e);
  double before=0.0;
  bool reverse=n>1 && segments[1]!=(segments[0]+1)%e->limiter_n;
  for (int i=0; i<n; ++i) {
    int s=segments[i], t=(s+1)%e->limiter_n;
    if (reverse) { int tmp=s; s=t; t=tmp; }
    double a[2]={e->limiter_R[s],e->limiter_Z[s]};
    double b[2]={e->limiter_R[t],e->limiter_Z[t]};
    double dx=b[0]-a[0], dy=b[1]-a[1], length=hypot(dx,dy);
    if (tok_wall_target_point_on_segment(p,a,b,tol)) {
      double along=(p[0]-a[0])*(dx/length)+(p[1]-a[1])*(dy/length);
      if (!tok_wall_target_finite(along)) return false;
      double result=before+fmax(0.0,fmin(length,along));
      if (!tok_wall_target_finite(result)) return false;
      *coordinate=result;
      return true;
    }
    before+=length;
  }
  return false;
}

// Arclength in metres from the start of the ordered target list. For a
// one-segment list, use that segment's native i -> i+1 direction. At a join,
// prefer the first containing segment within the shared wall tolerance.
// Failure leaves *coordinate untouched; no nearest-wall projection is granted
// to a point farther than the membership tolerance from this supplied arc.
static inline bool
tok_wall_target_coordinate(const struct gkyl_efit *e, int n,
  const int *segments, const double p[2], double *coordinate)
{
  return tok_wall_target_valid(e,n,segments) &&
    tok_wall_target_coordinate_validated(e,n,segments,p,coordinate);
}

static inline bool
tok_wall_target_contains(const struct gkyl_efit *e, int n,
  const int *segments, const double p[2])
{
  double coordinate;
  return tok_wall_target_coordinate(e,n,segments,p,&coordinate);
}

// Both endpoints belong to the same explicitly declared material subarc.
// This does NOT certify the straight chord a--b, a polynomial cap, native
// flux-root continuity, or the intervening volume map. Callers must validate
// those separately and use this only for a physical material theta face.
static inline bool
tok_wall_target_cap(const struct gkyl_efit *e, int n, const int *segments,
  const double a[2], const double b[2])
{
  double sa,sb;
  return tok_wall_target_valid(e,n,segments) &&
    tok_wall_target_coordinate_validated(e,n,segments,a,&sa) &&
    tok_wall_target_coordinate_validated(e,n,segments,b,&sb);
}
