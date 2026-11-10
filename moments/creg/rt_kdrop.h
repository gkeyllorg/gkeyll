#pragma once
#include <math.h>
#include <stdbool.h>

// Fixed-velocity gamma-law remap: E = W^2 rho + [G W^2/(G-1)-1] K rho^G.
// W^2=1 gives the fluid-frame energy-preserving variant. The bracket is
// monotone for G>1; do not hide a failed solve with a density/pressure floor.
static inline bool
rt_kdrop_density(double G, double K, double W2, double E, double *rho)
{
  if (!(isfinite(G) && G>1 && isfinite(K) && K>0
      && isfinite(W2) && W2>=1 && isfinite(E) && E>0)) return false;
  double lo=0, hi=E/W2, c=G*W2/(G-1)-1;
  for (int it=0; it<100; ++it) {
    double mid=.5*(lo+hi), value=W2*mid+c*K*pow(mid,G);
    if (!isfinite(value)) return false;
    if (value>E) hi=mid; else lo=mid;
  }
  *rho=.5*(lo+hi);
  return *rho>0 && fabs(W2*(*rho)+c*K*pow(*rho,G)-E)<=1e-13*E;
}

// Preserve both Eulerian constraint sources at fixed metric, allowing velocity
// to change: S_i=(E+p)v_i and rho=E-p/(G-1)-S^2/(E+p). For 1<G<=2 this
// rho is decreasing in p, so p-K*rho(p)^G has one root in [0,(G-1)E].
static inline bool
rt_kdrop_es_primitives(double G, double K, double E, double S2,
  double *rho, double *pressure, double *W)
{
  if (!(isfinite(G) && G>1 && G<=2 && isfinite(K) && K>0
      && isfinite(E) && E>0 && isfinite(S2) && S2>=0 && S2<E*E)) return false;
  double lo=0, hi=(G-1)*E;
  for (int it=0; it<100; ++it) {
    double p=.5*(lo+hi), r=E-p/(G-1)-S2/(E+p);
    double f=p-K*pow(fmax(r,0),G);
    if (!isfinite(f)) return false;
    if (f>0) hi=p; else lo=p;
  }
  double p=.5*(lo+hi), r=E-p/(G-1)-S2/(E+p);
  double v2=S2/((E+p)*(E+p));
  if (!(r>0 && p>0 && v2<1 && fabs(p-K*pow(r,G))<=1e-12*fmax(p,K*pow(r,G)))) return false;
  *rho=r; *pressure=p; *W=1/sqrt(1-v2);
  return true;
}

// q[1..3] store densitized COVARIANT momenta. Write only D and tau, and
// only after a valid solve; geometry, momentum and WB reference slots stay put.
static inline bool
rt_kdrop_es_conserved(double G, double K, double sqrt_det, const double inv[3][3],
  double rho_floor, double p_floor, double *q, double *pressure)
{
  if (!(isfinite(sqrt_det) && sqrt_det>0)) return false;
  double energy=q[0]+q[4], E=energy/sqrt_det, S2=0, rho, p, W;
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j) S2+=inv[i][j]*(q[1+i]/sqrt_det)*(q[1+j]/sqrt_det);
  if (!rt_kdrop_es_primitives(G,K,E,S2,&rho,&p,&W) || rho<rho_floor || p<p_floor) return false;
  double D=sqrt_det*rho*W;
  if (!(isfinite(D) && D>0 && isfinite(energy-D))) return false;
  q[0]=D; q[4]=energy-D; *pressure=p;
  return true;
}
