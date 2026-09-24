#pragma once

#include <math.h>

// Three-component vector
struct  gkyl_vec3 { double x[3]; };

// new vector with zeros 
static inline struct gkyl_vec3
gkyl_vec3_zeros()
{
  return (struct gkyl_vec3) { .x = { 0.0, 0.0, 0.0} };
}

// new vector with x,y,z components
static inline struct gkyl_vec3
gkyl_vec3_new(double x, double y, double z)
{
  return (struct gkyl_vec3) { .x = { x, y, z} };
}

// scalar*a
static inline struct gkyl_vec3
gkyl_vec3_scale(double scalar, struct gkyl_vec3 a)
{
  return (struct gkyl_vec3) { .x = { scalar*a.x[0], scalar*a.x[1], scalar*a.x[2] } };
}

// a+b
static inline struct gkyl_vec3
gkyl_vec3_add(struct gkyl_vec3 a, struct gkyl_vec3 b)
{
  return (struct gkyl_vec3) { .x = { a.x[0]+b.x[0], a.x[1]+b.x[1], a.x[2]+b.x[2] } };
}

// a-b
static inline struct gkyl_vec3
gkyl_vec3_sub(struct gkyl_vec3 a, struct gkyl_vec3 b)
{
  return (struct gkyl_vec3) { .x = { a.x[0]-b.x[0], a.x[1]-b.x[1], a.x[2]-b.x[2] } };
}

// |a|
static inline double
gkyl_vec3_len(struct gkyl_vec3 a)
{
  return sqrt( a.x[0]*a.x[0] + a.x[1]*a.x[1] + a.x[2]*a.x[2] );
}

// normalize a
static inline struct gkyl_vec3
gkyl_vec3_norm(struct gkyl_vec3 a)
{
  double len = gkyl_vec3_len(a);
  return (struct gkyl_vec3) { .x = { a.x[0]/len, a.x[1]/len, a.x[2]/len } };
}

// a \dot b
static inline double
gkyl_vec3_dot(struct gkyl_vec3 a, struct gkyl_vec3 b)
{
  return a.x[0]*b.x[0] + a.x[1]*b.x[1] + a.x[2]*b.x[2];
}

// a \times b
static inline struct gkyl_vec3
gkyl_vec3_cross(struct gkyl_vec3 a, struct gkyl_vec3 b)
{
  return (struct gkyl_vec3) { .x = {
      a.x[1]*b.x[2]-a.x[2]*b.x[1],
      a.x[2]*b.x[0]-a.x[0]*b.x[2],
      a.x[0]*b.x[1]-a.x[1]*b.x[0]
    }
  };
}

// a \dot (b \times c)
static inline double
gkyl_vec3_triple(struct gkyl_vec3 a, struct gkyl_vec3 b, struct gkyl_vec3 c)
{
  return gkyl_vec3_dot(a, gkyl_vec3_cross(b, c));
}

/**
 * Convert contravariant components of a vector in polar coordinate
 * basis to Cartesian basis components. Note the polar basis are not
 * normalized.
 *
 * @param r Radial coordinate
 * @param phi Polar angle
 * @param pin Input contravariant vector components
 * @return Cartesian component for @a pin
 */
static inline struct gkyl_vec3
gkyl_vec3_polar_con_to_cart(double r, double phi, struct gkyl_vec3 pin)
{
  double c0 = pin.x[0] * cos(phi) - pin.x[1] * sin(phi) * r;
  double c1 = pin.x[0] * sin(phi) + pin.x[1] * cos(phi) * r;
  double c2 = pin.x[2];
  return gkyl_vec3_new(c0, c1, c2);
}

/**
 * Convert contravariant components of a vector in polar coordinate
 * basis to Cartesian basis components. Note the polar basis are not
 * normalized.
 *
 * @param r Radial coordinate
 * @param phi Polar angle
 * @param pin Input contravariant vector components
 * @return Cartesian component for @a pin
 */
static inline struct gkyl_vec3
gkyl_vec3_polar_cov_to_cart(double r, double phi, struct gkyl_vec3 pin)
{
  double c0 = pin.x[0]*cos(phi)-(pin.x[1]*sin(phi))/r;
  double c1 = pin.x[0]*sin(phi)+(pin.x[1]*cos(phi))/r;
  double c2 = pin.x[2];
  return gkyl_vec3_new(c0, c1, c2);
}    

/**
 * Convert contravariant components of a vector in polar coordinate
 * basis to Cartesian basis components. Note the polar basis are not
 * normalized.
 *
 * @param r Radial coordinate
 * @param phi Polar angle
 * @param pin Input contravariant vector components
 * @return Cartesian component for @a pin
 */
static inline struct gkyl_vec3
gkyl_vec3_polar_con_to_cov(double r, struct gkyl_vec3 pin)
{
  double c0 = pin.x[0];
  double c1 = pin.x[1] * r * r;
  double c2 = pin.x[2];
  return gkyl_vec3_new(c0, c1, c2);
}

// minmod(a,b,c, ...) returns min of all parameters are positive, max
// of all parameters if all are negative, and zero otherwise.


// minmod(x,y)
static inline double
gkyl_minmod_2(double x, double y)
{
  if (x>0 && y>0)
    return fmin(x,y);
  if (x<0 && y<0)
    return fmax(x,y);
  return 0.0;
}

// minmod(x,y)
static inline double
gkyl_minmod_3(double x, double y, double z)
{
  if (x>0 && y>0 && z>0)
    return fmin(fmin(x,y),z);
  if (x<0 && y<0 && z<0)
    return fmax(fmax(x,y),z);
  return 0.0;
}

// minmod(x,y,z,w)
static inline double
gkyl_minmod_4(double x, double y, double z, double w)
{
  if (x>0 && y>0 && z>0 && w>0)
    return fmin(fmin(x,y),fmin(z,w));
  if (x<0 && y<0 && z<0 && w<0)
    return fmax(fmax(x,y),fmax(z,w));
  return 0.0;
}

// median(x,y,z) returns the value between two others
static inline double
gkyl_median(double x, double y, double z)
{
  return x + gkyl_minmod_2(y-x,z-x);
}

// min(x,y,z)
static inline double
gkyl_min_3(double x, double y, double z)
{
  return fmin(x,fmin(y,z));
}

// max(x,y,z)
static inline double
gkyl_max_3(double x, double y, double z)
{
  return fmax(x,fmax(y,z));
}

// Result from numerical quadrature
struct gkyl_qr_res {
  double res; // result of quadrature
  double error; // error estimate
  int status; // 0 - success.
  int nevals; // number of function evaluations
  int nlevels; // number of levels
};

// polynomial order for low-order polynomial
enum gkyl_lo_poly_order {
  GKYL_LO_POLY_2 = 2,
  GKYL_LO_POLY_3,
  GKYL_LO_POLY_4,  
};  

// Polynomial roots for low-order polynomials (upto 4th order)
struct gkyl_lo_poly_roots {
  double rpart[4], impart[4]; // real and imaginary part of roots
  double err[4]; // estimated error for each root
  int niter; // number of iterations
};

// Polynomial roots for arbitrary order polynomials: allocate an free
// using new and release methods as usual
struct gkyl_poly_roots {
  int poly_order; // polynomial order
  double *rpart, *impart; // real and imaginary part of roots
  double *err; // estimated error for each root
  int niter;   // number of iterations
  void *work; // some memory needed internally. Do not muck or access!
};  

// Quartic polynomial to hand off to ridders
struct gkyl_quartic_polynomial {
  double p[4];
};

// Result polynomials from Sturn Chain
struct gkyl_sturm_polynomials {
  double p0[4];
  double p1[4];
  double p2[4];
  double p3[4];
  double p4[4];
};

// Result from binary search/ strun chain
struct gkyl_root_intervals {
  double root_bound_lower[4]; // lower root bound
  double root_bound_upper[4]; // upper root bound
  int status; // 0 - success.
  int status_refinement[4]; // 0 - success.
  int niter; // number of iterations
  int niter_refinement[4]; // number of iterations
  int nroots; // number of distinct-real-roots
  struct gkyl_sturm_polynomials sturn_chain;
  double real_roots_ridders[4]; // Output of the root finding algorithm
  int status_ridders[4]; // Status of the ridders' ability to find the roots 
};

/**
 * Integrate func(x,ctx) in interval [a,b], using maximum levels
 * "n". Typical value of levels is 7. The status flag of the return
 * struct is set to 1 if quadrature failed due to number of levels
 * exceeded.
 *
 * @param func Function to integrate
 * @param ctx Context to pass to function
 * @param a Lower limit of interval
 * @param b Upper limit of interval
 * @param n Number of levels
 * @param eps Error to compute integral
 * @return Result of quadrature
 */
struct gkyl_qr_res gkyl_dbl_exp(double (*func)(double, void *), void *ctx,
  double a, double b, int n, double eps);

/**
 * Compute single real root of a function using Ridders' method. See
 * IEEE Tran. Circuit and Systems, vol CAS-26 No 11, Pg 976 1976. The
 * root must be inside the interval specified, and there must only be
 * one such root in the interval.
 *
 * @param func Function to integrate
 * @param ctx Context to pass to function
 * @param x1 Lower limit of interval
 * @param x2 Upper limit of interval
 * @param f1 Function at x1
 * @param f2 Function at x2
 * @param max_iter Maximum number of iterations
 * @param eps Desired error bound
 * @return Result of root finder
 */
struct gkyl_qr_res gkyl_ridders(double (*func)(double,void*), void *ctx,
  double x1, double x2, double f1, double f2, int max_iter, double eps);

/**
 * Compute all simple roots of the low-order polynomial with monomial
 * coefficients given by @a coeff. The leading coefficient is always *
 * assumed be 1.0 and so the coeff[i] give the coefficient for the *
 * monomial x^i. For example for the quartic:
 *
 * p(x) = x^4 + 2 x^3 + 4
 *
 * coeff[4] = { 4.0, 0.0, 0.0, 2.0 };
 *
 * @param order Order of polynomial
 * @param coeff Coefficients of the polynomial
 * @return Roots of the polynomial
 */
struct gkyl_lo_poly_roots gkyl_calc_lo_poly_roots(enum gkyl_lo_poly_order order,
  double coeff[4]);

/**
 * Allocate memory to store polynomial roots.
 *
 * @param poly_order Polynomial order
 * @return newly allocated memory. Free using release method.
 */
struct gkyl_poly_roots* gkyl_poly_roots_new(int poly_order);

/**
 * Compute all simple roots of polymomial with real coefficients. The
 * leading coefficient must be 1.0 and the coeff[i] is the coefficient
 * of x^i.
 *
 * @param pr On output, contains all the roots. Must be preallocated
 * @param coeff Coefficients of monomials
 */
void gkyl_calc_poly_roots(struct gkyl_poly_roots *pr, const double *coeff);

/**
 * Compute all *real-district-root* intervals to a *quartic* equation via Sturn's sequence
 * and isolate individual intervals via bisection search. Once every real-distict root 
 * interval is isolated, then the root intervals are returned. (There is no further 
 * iterative isolation)
 * 
 * The leading coefficient 
 * is always assumed be 1.0 and so the coeff[i] give the coefficient for the
 * monomial x^i. For example:
 *
 * p(x) = x^4 + 2 x^3 + 4
 *
 * coeff[4] = { 4.0, 0.0, 0.0, 2.0 };
 *
 * @param coeff Coefficients of the polynomial
 * @param domain Interval to find the real-district roots within
 * @param tol Tolerance of the quartic solve
 * @return Roots of the polynomial
 */
struct gkyl_root_intervals gkyl_calc_quartic_root_intervals(
  double coeff[4], double domain[2], double tol);

/**
 * Refine the result of gkyl_calc_quartic_root_intervals() using bisection search
 * to narrow down the domain containing the root. 
 * 
 * The leading coefficient 
 * is always assumed be 1.0 and so the coeff[i] give the coefficient for the
 * monomial x^i. For example:
 *
 * p(x) = x^4 + 2 x^3 + 4
 *
 * coeff[4] = { 4.0, 0.0, 0.0, 2.0 };
 *
 * @param root_intervals Object containing the number of roots and root intervals
 * @param tol Tolerance of the interval isolation
 * @return Roots of the polynomial
 */
void gkyl_refine_root_intervals_bisection(struct gkyl_root_intervals *root_intervals,
  double tol);


/**
 * Compute the roots from the intervals given by gkyl_calc_quartic_root_intervals() 
 * using ridders' algorithm. This can handle either refined or unrefined intervals.
 * 
 * The leading coefficient 
 * is always assumed be 1.0 and so the coeff[i] give the coefficient for the
 * monomial x^i. For example:
 *
 * p(x) = x^4 + 2 x^3 + 4
 *
 * coeff[4] = { 4.0, 0.0, 0.0, 2.0 };
 *
 * @param root_intervals Object containing the number of roots and root intervals
 * @param tol Tolerance of the interval isolation
 * @return Roots of the polynomial
 */
void gkyl_root_isolation_from_intervals_via_ridders(struct gkyl_root_intervals *root_intervals,
  double tol);

/**
 * Release memory for use in polynomial root finder.
 *
 * @param pr Memory to release.
 */
void gkyl_poly_roots_release(struct gkyl_poly_roots *pr);

