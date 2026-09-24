#include <gkyl_alloc.h>
#include <gkyl_math.h>

#include <float.h>
#include <stdbool.h>
#include <complex.h>

// Temporary for diagnostics
#include <stdio.h>
static double ROOT_EPS = 1e-14;

// This implementation is taken from the Appendix of "Improving the
// Double Exponential Quadrature Tanh-Sinh, Sinh-Sinh and Exp-Sinh
// Formulas", Robert A. van Engelen. See
// https://www.genivia.com/qthsh.html
//
// PS: I do not think the improvements claimed in the above note are
// actually correct. Hence, at present I am using the wp34s
// implementation in Appendix B of the note.

struct gkyl_qr_res
gkyl_dbl_exp(double (*func)(double, void *), void *ctx,
  double a, double b, int n, double eps)
{
  int nev = 0;  
  double thr = 10*sqrt(eps); // too generous for larger eps, e.g. eps=1e-9
  //double thr = eps; // too generous for larger eps, e.g. eps=1e-9
  double c = (a+b)/2; // center (mean)
  double d = (b-a)/2; // half distance
  double s = func(c, ctx); nev += 1;
  double fp = 0, fm = 0;
  double p, e, v, h = 2;
  double tmax = log(2/M_PI * log((d < 1 ? 2*d : 2) / eps));
  int k = 0; // level
  do {
    double q, t;
    int j = 1;
    v = s*d*M_PI/2*h; // last sum
    p = 0;
    h /= 2;
    t = h;
    do {
      double ch = cosh(t);
      double ecs = cosh(M_PI/2 * sqrt(ch*ch - 1)); // = cosh(pi/2*sinh(t))
      double w = 1/(ecs*ecs);
      double r = sqrt(ecs*ecs - 1)/ecs;
      double x = d*r;
      if (c+x > a) {
        double y = func(c+x, ctx); nev += 1;
        if (isfinite(y))
          fp = y;
      }
      if (c-x < b) {
        double y = func(c-x, ctx); nev += 1;
        if (isfinite(y)) 
          fm = y;
      }
      q = ch*w*(fp+fm);
      p += q;
      j += 1+(k>0);
      t = j*h;
    } while (t <= tmax && fabs(q) > eps*fabs(p));
    s += p;
    ++k;
  } while (s && fabs(2*fabs(p) - fabs(s)) >= fabs(thr*s) && k <= n);
  s *= d*M_PI/2*h;
  e = fabs(v-s);
  if (10*e >= fabs(s)) {
    e += fabs(s);
    s = 0;
  }
  
  return (struct gkyl_qr_res) {
    .error = e,
    .res = s,
    .nevals = nev,
    .status = k>n ? 1 : 0,
    .nlevels = k
  };
}

// Helper functions for ridders
static inline double dsign(double x) { return x >= 0 ? 1 : -1; }

// See IEEE Tran. Circuit and Systems, vol CAS-26 No 11, Pg 976
// 1976. The following is almost direct implementation from the
// original paper
struct gkyl_qr_res
gkyl_ridders(double (*func)(double,void*), void *ctx,
  double xl, double xr, double fl, double fr, int max_iter, double eps)
{
  double x0 = xl, x2 = xr, f0 = fl, f2 = fr;
  double res = DBL_MAX, err = DBL_MAX;

  int nev = 0, nitr = 0, iterating = 1;
  while (iterating && nitr <= max_iter) {
    double x1 = 0.5*(x0+x2);
    double f1 = func(x1, ctx); nev += 1;
    double W = f1*f1 - f0*f2;

    double d = x2-x1;
    double x3 = x1 + dsign(f0)*f1*d/sqrt(W);
    double f3 = func(x3, ctx); nev += 1;

    if (fabs(res-x3) < eps) {
      err = fabs(res-x3);
      iterating = 0;
    }
    res = x3;

    if (f3*f0 < 0) {
      x2 = x3;
      f2 = f3;
    }
    else if (f3*f1 < 0) {
      x0 = x1 < x3 ? x1 : x3;
      f0 = x1 < x3 ? f1 : f3;

      x2 = x1 < x3 ? x3 : x1;
      f2 = x1 < x3 ? f3 : f1;
    }
    else if (f3*f2 < 0 ) {
      x0 = x3;
      f0 = f3;
    }
    nitr++;
  }
  
  return (struct gkyl_qr_res) {
    .error = err,
    .res = res,
    .nevals = nev,
    .status = nitr>max_iter ? 1 : 0,
  };
}

/////// roots of a quadratic polynomial
static struct gkyl_lo_poly_roots
quad_poly_roots(double coeff[4])
{
  double c = coeff[0], b = coeff[1], a = 1.0;

  double complex x1 = 0.0, x2 = 0.0;
  if (b>=0.0) {
    x1 = (-b-csqrt(b*b-4*a*c))/(2*a);
    x2 = 2*c/(-b-csqrt(b*b-4*a*c));
  }
  else {
    x1 = 2*c/(-b+csqrt(b*b-4*a*c));
    x2 = (-b+csqrt(b*b-4*a*c))/(2*a);
  }
  
  return (struct gkyl_lo_poly_roots) {
    .niter = 0,
    .err = { 0.0, 0.0 },
    .rpart = { creal(x1), creal(x2) },
    .impart = { cimag(x1), cimag(x2) }
  };
}

/////// roots of a cubic
static inline double complex
eval_poly3(const double coeff[4], double complex x)
{
  double complex x2 = x*x;
  double complex x3 = x2*x;
  return coeff[0] + coeff[1]*x + coeff[2]*x2 + x3;
}

static inline bool
check_converged3(double complex c1, double complex c2, double complex c3, double err[3])
{
  double eps = ROOT_EPS;
  err[0] = cabs(c1);
  err[1] = cabs(c2);
  err[2] = cabs(c3);
  
  return (err[0]+err[1]+err[2])/3.0 < eps;
}

// roots of a cubic polynomial (Durand-Kerner method)
static struct gkyl_lo_poly_roots
cubic_poly_roots(double coeff[4])
{
  double complex r1 = 0.4+0.9*I; // arbitrary complex number, not a root of unity
  double complex pn1 = r1;
  double complex qn1 = pn1*r1;
  double complex rn1 = qn1*r1;
  double complex pn = 0.0, qn = 0.0, rn = 0.0;

  double err[3] = { 0.0 };
  int max_iter = 100;
  int niter = 0;
  do {
    pn = pn1; qn = qn1; rn = rn1;
    
    pn1 = pn1 - eval_poly3(coeff, pn1)/( (pn1-qn1)*(pn1-rn1) );
    qn1 = qn1 - eval_poly3(coeff, qn1)/( (qn1-pn1)*(qn1-rn1) );
    rn1 = rn1 - eval_poly3(coeff, rn1)/( (rn1-pn1)*(rn1-qn1) );

    niter += 1;
    
  } while( !check_converged3(pn-pn1, qn-qn1, rn-rn1, err) && niter < max_iter );

  return (struct gkyl_lo_poly_roots) {
    .niter = niter,
    .err = { err[0], err[1], err[2] },
    .rpart = { creal(pn1), creal(qn1), creal(rn1) },
    .impart = { cimag(pn1), cimag(qn1), cimag(rn1) }
  };
}

/////// roots of a quartic
static inline double complex
eval_poly4(const double coeff[4], double complex x)
{
  double complex x2 = x*x;
  double complex x3 = x2*x;
  double complex x4 = x2*x2;
  return coeff[0] + coeff[1]*x + coeff[2]*x2 + coeff[3]*x3 + x4;
}

static inline bool
check_converged4(double complex c1, double complex c2, double complex c3, double complex c4,
  double err[4])
{
  double eps = ROOT_EPS;
  err[0] = cabs(c1);
  err[1] = cabs(c2);
  err[2] = cabs(c3);
  err[3] = cabs(c4);  
  
  return (err[0]+err[1]+err[2]+err[3])/4.0 < eps;
}

// roots of a quartic polynomial (Durand-Kerner method)
static struct gkyl_lo_poly_roots
quart_poly_roots(double coeff[4])
{
  double complex r1 = 0.4+0.9*I; // arbitrary complex number, not a root of unity
  double complex pn1 = r1;
  double complex qn1 = pn1*r1;
  double complex rn1 = qn1*r1;
  double complex sn1 = rn1*r1;
  double complex pn = 0.0, qn = 0.0, rn = 0.0, sn = 0.0;

  double err[4] = { 0.0 };
  int max_iter = 100;
  int niter = 0;
  do {
    pn = pn1; qn = qn1; rn = rn1; sn = sn1;
    
    pn1 = pn1 - eval_poly4(coeff, pn1)/( (pn1-qn1)*(pn1-rn1)*(pn1-sn1) );
    qn1 = qn1 - eval_poly4(coeff, qn1)/( (qn1-pn1)*(qn1-rn1)*(qn1-sn1) );
    rn1 = rn1 - eval_poly4(coeff, rn1)/( (rn1-pn1)*(rn1-qn1)*(rn1-sn1) );
    sn1 = sn1 - eval_poly4(coeff, sn1)/( (sn1-pn1)*(sn1-qn1)*(sn1-rn1) );

    niter += 1;
    
  } while( !check_converged4(pn-pn1, qn-qn1, rn-rn1, sn-sn1, err) && niter < max_iter );

  return (struct gkyl_lo_poly_roots) {
    .niter = niter,
    .err = { err[0], err[1], err[2], err[3] },
    .rpart = { creal(pn1), creal(qn1), creal(rn1), creal(sn1) },
    .impart = { cimag(pn1), cimag(qn1), cimag(rn1), cimag(sn1) }
  };
}

struct gkyl_lo_poly_roots
gkyl_calc_lo_poly_roots(enum gkyl_lo_poly_order order, double coeff[4])
{
  struct gkyl_lo_poly_roots proots = { };
  switch(order) {
    case GKYL_LO_POLY_2:
      proots = quad_poly_roots(coeff);
      break;
    case GKYL_LO_POLY_3:
      proots = cubic_poly_roots(coeff);
      break;
    case GKYL_LO_POLY_4:
      proots = quart_poly_roots(coeff);
      break;
  }
  return proots;
}

static inline double 
eval_poly_lo(double *p, double x)
{

  // eval_poly_lo() assumes a fourth order quartic form with 1.0*x^4 leading term is 0

  //compute the polynomial evaluated at x
  double res = 0.0;

  for (int i=0; i<=3; ++i) {
    res = res + p[i]*pow(x,i);
  }

  // return the resulted evalued polynomial
  return res;
}


static inline double 
eval_poly(double x, void *ctx)
{

  // eval_poly() assumes a fourth order quartic form with 1.0*x^4 leading term,
  // bounds must be finite
  double *p = (double*)ctx; 

  //compute the polynomial evaluated at x
  double res = 0.0;

  for (int i=0; i<=3; ++i) {
    res = res + p[i]*pow(x,i);
  }
  res = res + pow(x,4);

  // return the resulted evalued polynomial
  return res;
}


static inline int
sign_changes(int *signs)
{
  int res = 0;
  for (int i=0; i<=3; ++i) {
    if (signs[i]*signs[i+1] < 0.0){
      res = res + 1;
    }
  }
  return res;
}


static void
signs_strum_chain(double *eval_sturn_chain, int *signs)
{
  int iter = 0;
  for (int i=0; i<=4; ++i) {
    signs[i] = 0;
    if (eval_sturn_chain[i] > 0.0){
      signs[iter] = 1;
      iter = iter + 1;
    } else if (eval_sturn_chain[i] < 0.0){
      signs[iter] = -1;
      iter = iter + 1;
    }
  }
}


static void
eval_strum_chain(struct gkyl_sturm_polynomials *sturn_chain, double x, double *eval)
{
  eval[0] = eval_poly(x,sturn_chain->p0);
  eval[1] = eval_poly_lo(sturn_chain->p1,x);
  eval[2] = eval_poly_lo(sturn_chain->p2,x);
  eval[3] = eval_poly_lo(sturn_chain->p3,x);
  eval[4] = eval_poly_lo(sturn_chain->p4,x);
}


static double
eval_num_roots(struct gkyl_sturm_polynomials *sturn_chain, double domain[2])
{
  // Evaluate the sturn chain
  double eval_sturn_chain_l[5];
  double eval_sturn_chain_r[5];
  eval_strum_chain(sturn_chain,domain[0],eval_sturn_chain_l);
  eval_strum_chain(sturn_chain,domain[1],eval_sturn_chain_r);

  // Compute the signs 
  int signs_l[5];
  int signs_r[5];
  signs_strum_chain(eval_sturn_chain_l,signs_l);
  signs_strum_chain(eval_sturn_chain_r,signs_r);

  // Compute the number of roots
  int num_roots = sign_changes(signs_l) - sign_changes(signs_r);

  // Return the number of roots
  return num_roots;
}


static void
check_poly_full_domain(double *p, double *domain, double tol)
{
  double left_bound = domain[0];
  double right_bound = domain[1];
  double eval_left_bound = eval_poly(left_bound,p);
  double eval_right_bound = eval_poly(right_bound,p);

  // If the edges are zero, then shift the boundary
  int iter = 0;
  while(eval_left_bound == 0.0 && iter <= 3){
    left_bound = left_bound - tol;
    eval_left_bound = eval_poly(left_bound,p);
    iter = iter + 1;
  }
  iter = 0;
  while(eval_right_bound == 0.0 && iter <= 3){
    right_bound = right_bound + tol;
    eval_right_bound = eval_poly(right_bound,p);
    iter = iter + 1;
  }

  // Update the boundaries
  domain[0] = left_bound;
  domain[1] = right_bound;;
}


static double 
check_poly_bounded(double *p, double left_bound, double middle_bound, double right_bound)
{

  // Compute the updated middle bound in the polynomial evaluation falls directly on 
  // p(x) = 0
  double updated_middle_bound = middle_bound;
  double eval_middle_bound = eval_poly(middle_bound,p);

  int iter = 0;
  while(eval_middle_bound == 0.0 && iter <= 3){
      updated_middle_bound = (right_bound+left_bound)/2.0 + (right_bound-left_bound)*(0.1*(iter+1));
      eval_middle_bound = eval_poly(updated_middle_bound,p);
      iter = iter + 1;
  }

  // Check the new V
  bool bounded = (left_bound < updated_middle_bound && updated_middle_bound < right_bound);
  if (bounded != 1){
    // return original boundary 
    updated_middle_bound = middle_bound;
  }

  return updated_middle_bound;
}


static struct gkyl_root_intervals
bisection_root_search(struct gkyl_sturm_polynomials *sturn_chain, double domain[2], double tol)
{

  // Initialize 
  int status = 0;
  int niter = 0;
  double lower_bound[4] = {0.0, 0.0, 0.0, 0.0};
  double upper_bound[4] = {0.0, 0.0, 0.0, 0.0};

  // Check the polynomial evaluatation at the bounds is not zero
  check_poly_full_domain(sturn_chain->p0, domain, tol);

  // Compute number of real-distinct-roots
  int total_roots = eval_num_roots(sturn_chain,domain);

  // Check cases where there are no roots/only one root
  if (total_roots == 1){

    // Complete
    lower_bound[0] = domain[0];
    upper_bound[0] = domain[1];

  } else if (total_roots == 0){

    // Complete

  // Bisection search begins
  } else if (total_roots > 1){

    // Initialize additional variables needed in the search
    int roots_isolated = 0;
    double return_to_this_domain[2] = {0.0,0.0};

    // Break the domain into two parts
    double left_bound = domain[0];
    double right_bound = domain[1];
    double middle_bound = (left_bound+right_bound)/2.0;
    middle_bound = check_poly_bounded(sturn_chain->p0, left_bound, middle_bound, right_bound);
    int not_all_roots_are_isolated = 1;
    int isolated_num_roots = 0;
    double domain_left[2];
    double domain_right[2];
    int left_num_roots, right_num_roots;

    // Iterate while not all roots are isolated
    while(not_all_roots_are_isolated && niter < 1000){

      // update niter
      niter = niter + 1;

      // Compute the roots (L)
      domain_left[0] = left_bound;
      domain_left[1] = middle_bound;
      left_num_roots = eval_num_roots(sturn_chain, domain_left);

      // Compute the roots (R)
      domain_right[0] = middle_bound;
      domain_right[1] = right_bound;
      right_num_roots = eval_num_roots(sturn_chain, domain_right);


      //Cases of different number of roots in each side
      if (left_num_roots == total_roots - isolated_num_roots){

        // update to the left domain
        right_bound = middle_bound;
        middle_bound = (left_bound+right_bound)/2.0;
        middle_bound = check_poly_bounded(sturn_chain->p0, left_bound, middle_bound, right_bound);

      } else if (right_num_roots == total_roots - isolated_num_roots){

        // update to the right domain
        left_bound = middle_bound;
        middle_bound = (left_bound+right_bound)/2.0;
        middle_bound = check_poly_bounded(sturn_chain->p0, left_bound, middle_bound, right_bound);

      // We succussfully split the domain to some degree
      } else { 

        //if the left or right is 1 then we can add a domain to our result
        if (left_num_roots == 1 || right_num_roots == 1){

          if (left_num_roots == 1){
            lower_bound[roots_isolated] = left_bound;
            upper_bound[roots_isolated] = middle_bound;
            roots_isolated = roots_isolated + 1;
            isolated_num_roots = roots_isolated;
          }

          if (right_num_roots == 1){
            lower_bound[roots_isolated] = middle_bound;
            upper_bound[roots_isolated] = right_bound;
            roots_isolated = roots_isolated + 1;
            isolated_num_roots = roots_isolated;
          }

          //update to the right-side of the domain
          if (left_num_roots == 1){
            left_bound = middle_bound;
            middle_bound = (left_bound+right_bound)/2.0;
            middle_bound = check_poly_bounded(sturn_chain->p0, left_bound, middle_bound, right_bound);
          }

          // update to the left-side of the domain
          if (right_num_roots == 1){
            right_bound = middle_bound;
            middle_bound = (left_bound+right_bound)/2.0;
            middle_bound = check_poly_bounded(sturn_chain->p0, left_bound, middle_bound, right_bound);
          }

          if (right_num_roots == 1 && left_num_roots == 1){

            // Complete if all roots are isolated
            if (total_roots == roots_isolated){
              not_all_roots_are_isolated = 0;
            } else {

              // Otherwise, load the second domain
              left_bound = return_to_this_domain[0];
              right_bound = return_to_this_domain[1];
              middle_bound = (left_bound+right_bound)/2.0;
              middle_bound = check_poly_bounded(sturn_chain->p0, left_bound, middle_bound, right_bound);
              isolated_num_roots = roots_isolated;

            } // end returning to isolated domain
          }

        } else { // we've exactly split 2v2 in the domains

          if (right_num_roots == 2 && left_num_roots == 2){

            // Save a domain to return to then
            return_to_this_domain[0] = middle_bound;
            return_to_this_domain[1] = right_bound;
            isolated_num_roots = right_num_roots;

            // update to the left domain
            right_bound = middle_bound;
            middle_bound = (left_bound+right_bound)/2.0;
            middle_bound = check_poly_bounded(sturn_chain->p0, left_bound, middle_bound, right_bound);

          } else {
            status = 1;
          }

        } // end if for number of root selected

      } // end else due to successful domain divide

    } // end while

  } //end roots conditions

  // Throw an error if the while loop stopped due to unending iterations
  if (niter == 1000){
    printf("Couldn't isolate the root intervals: Total_roots: %d!\n",total_roots);
  }

  return (struct gkyl_root_intervals) {
    .status = status,
    .niter = niter,
    .nroots = total_roots,
    .root_bound_lower = { lower_bound[0], lower_bound[1], lower_bound[2], lower_bound[3] },
    .root_bound_upper = { upper_bound[0], upper_bound[1], upper_bound[2], upper_bound[3] },
    .sturn_chain = *sturn_chain
  };
}

static inline int
deg_modified(double *p)
{
  // If everything is zero, the degree is zero
  int res = 0;

  // Compute the degree
  for (int i=0; i<=3; ++i){
    if (p[i] != 0.0){
      res = i; 
    }
  }
  
  // return the result
  return res;
}


static void
minus_euclidean_division_rem(double *p0, double *p1, int p0_deg, double *res)
{

  // Initialize the problem 
  double r[4];
  for (int i=0; i<=3; ++i) r[i] = p0[i];
  int r_deg = p0_deg;
  double rmax;
  if (p0_deg < 4) {
    r_deg = deg_modified(r);
    rmax = r[r_deg];
  } else if (p0_deg == 4) {
    r_deg = p0_deg;
    rmax = 1.0;
  }
  int p1_deg = deg_modified(p1);
  int iter_max = r_deg-p1_deg; // TODO: check index
  double q[iter_max];
  double p1g[r_deg+1];
  for (int i=0; i<=r_deg; ++i) p1g[i] = 0.0;
  int iter = 0;

  // Iterate while the degree is still higher the the divisor
  while(r_deg >= p1_deg) {

    // compute the quotient
    q[iter_max - iter] = rmax/p1[p1_deg];
 
    // Multiply p1 by q
    for (int i=0; i<=p1_deg; ++i) {
        p1g[i+(iter_max - iter)] = p1[i]*q[iter_max - iter];
    }

    // Subtract p1g from r
    for (int i=0; i<=fmin(r_deg,3); ++i) {
      if (i != r_deg) {
        r[i] = r[i] - p1g[i];
      }
      else {
        r[i] = 0.0;
      }
    } 

    // update the iter, degree of the remaining poly
    iter = iter + 1;
    r_deg = r_deg - 1;
    rmax = r[r_deg];

    // set p1g to zero
    for (int i=0; i<=r_deg+1; ++i) p1g[i] = 0.0*p1g[i];
  }

  // We require the minus of the returned remainder
   for (int i=0; i<=3; ++i)  res[i] = -r[i];
}


static struct gkyl_sturm_polynomials 
compute_sturn_chain(double p0[4])
{

  // Compute the first derivative of p0 to get p1
  double p1[4] = { p0[1], 2.0*p0[2], 3.0*p0[3], 4.0*1.0 };

  // Compute the Euclidean division of p0 by p1, return the -remainder
  double p2[4], p3[4], p4[4];
  minus_euclidean_division_rem(p0,p1,4,p2);

  // Compute the next division p1 by p2
  minus_euclidean_division_rem(p1,p2,0,p3);

  // Compute the next division p2 by p3
  minus_euclidean_division_rem(p2,p3,0,p4);

  return (struct gkyl_sturm_polynomials) {
    .p0 = { p0[0], p0[1], p0[2], p0[3] },
    .p1 = { p1[0], p1[1], p1[2], p1[3] },
    .p2 = { p2[0], p2[1], p2[2], p2[3] },
    .p3 = { p3[0], p3[1], p3[2], p3[3] },
    .p4 = { p4[0], p4[1], p4[2], p4[3] }
  };
}

void 
gkyl_refine_root_intervals_bisection(struct gkyl_root_intervals *root_intervals, double tol)
{
  // Compute the number of domains to refine
  int nroots = root_intervals->nroots; 

  // For all domains, refine the result
  for (int i=0; i<nroots; ++i) {

    // Grab the domain
    double left_bound = root_intervals->root_bound_lower[i];
    double right_bound = root_intervals->root_bound_upper[i];
    double error = 2.0*tol;
    int iter = 0;
    double middle_bound;
    double domain_left[2];
    double domain_right[2];
    int left_num_roots, right_num_roots;
    int status = 0;

    // Iterate on the domain for some tolerance
    while(error > tol && iter < 100){

      // Grab middle location in the domain
      middle_bound = (left_bound+right_bound)/2.0;
      middle_bound = check_poly_bounded(root_intervals->sturn_chain.p0, left_bound, middle_bound, right_bound);

      // Compute the roots (L)
      domain_left[0] = left_bound;
      domain_left[1] = middle_bound;
      left_num_roots = eval_num_roots(&root_intervals->sturn_chain, domain_left);

      // Compute the roots (R)
      domain_right[0] = middle_bound;
      domain_right[1] = right_bound;
      right_num_roots = eval_num_roots(&root_intervals->sturn_chain, domain_right);


      // Cases of different number of roots in each side
      if (left_num_roots == 1){

          // update to the left domain
          right_bound = middle_bound;

      } else if (right_num_roots == 1){

          // update to the right domain
          left_bound = middle_bound;

      } else {

        status = 1;

      }

      // update the error
      error = fabs(left_bound - right_bound);
      iter = iter + 1;
    }
    
    // Save the restricted domain
    root_intervals->root_bound_lower[i] = left_bound;
    root_intervals->root_bound_upper[i] = right_bound;
    root_intervals->status_refinement[i] = status;
    root_intervals->niter_refinement[i] = iter;
  }
} 

void 
gkyl_root_isolation_from_intervals_via_ridders(struct gkyl_root_intervals *root_intervals, double tol)
{

  // Compute the number of domains to compute roots on
  int nroots = root_intervals->nroots; 
  double root;
  int status;

  // For all domains, refine the result
  for (int i=0; i<nroots; ++i) {

    // Compute the roots via ridders
    double x1 = root_intervals->root_bound_lower[i];
    double x2 = root_intervals->root_bound_upper[i];
    double f1 = eval_poly(x1,root_intervals->sturn_chain.p0);
    double f2 = eval_poly(x2,root_intervals->sturn_chain.p0);
    struct gkyl_qr_res res = gkyl_ridders(eval_poly, root_intervals->sturn_chain.p0, x1, x2, f1, f2, 100, tol);

    // Save the roots
    root_intervals->real_roots_ridders[i] = res.res;
    root_intervals->status_ridders[i] = res.status;
  }
}


struct gkyl_root_intervals 
gkyl_calc_quartic_root_intervals(double coeff[4], double domain[2], double tol)
{

  // Compute the sturn chain
  struct gkyl_sturm_polynomials sturn_chain = compute_sturn_chain(coeff);

  // Isolate the intervals
  struct gkyl_root_intervals root_intervals = bisection_root_search(&sturn_chain,domain, tol);

  // Retrun the structure
  return root_intervals;

} 
struct gkyl_poly_roots*
gkyl_poly_roots_new(int poly_order)
{
  struct gkyl_poly_roots *pr = gkyl_malloc(sizeof *pr);
  pr->impart = gkyl_malloc(sizeof(double[poly_order]));
  pr->rpart = gkyl_malloc(sizeof(double[poly_order]));
  pr->err = gkyl_malloc(sizeof(double[poly_order]));

  pr->work = gkyl_malloc(sizeof(double complex[2*poly_order]));

  pr->poly_order = poly_order;

  return pr;
}

static inline double complex
eval_poly2(int poly_order, const double *coeff, double complex x)
{
  double complex xn = x;
  double complex res = coeff[0];
  for (int i=1; i<poly_order; ++i) {
    res += coeff[i]*xn;
    xn = xn*x;
  }
  return res + xn;
}

static inline bool
check_converged(int n, double complex *p1, double complex *p2, double *err)
{
  double eps = ROOT_EPS;
  double tot_err = 0.0;
  for (int i=0; i<n; ++i) {
    err[i] = cabs(p1[i]-p2[i]);
    tot_err += err[i];
  }
  return tot_err/n < eps;
}

void
gkyl_calc_poly_roots(struct gkyl_poly_roots *pr, const double *coeff)
{
  int poly_order = pr->poly_order;
  double complex *pn1 = pr->work;
  double complex *pn = pn1 + poly_order;

  double complex r1 = 0.4+0.9*I; // arbitrary complex number, not a root of unity  
  pn1[0] = r1;
  for (int i=1; i<poly_order; ++i)
    pn1[i] = pn1[i-1]*r1;

  int max_iter = 100;
  int niter = 0;
  do {
    for (int i=0; i<poly_order; ++i) {
      pn[i] = pn1[i];
      
      double complex denom = 1.0;
      for (int j=0; j<poly_order; ++j)
        if (i != j)
          denom = denom*(pn1[i]-pn1[j]);
      pn1[i] = pn1[i] - eval_poly2(poly_order, coeff, pn1[i])/denom;
    }
    niter += 1;
  } while (!check_converged(poly_order, pn, pn1, pr->err) && niter < max_iter);

  pr->niter = niter;
  for (int i=0; i<poly_order; ++i) {
    pr->rpart[i] = creal(pn1[i]);
    pr->impart[i] = cimag(pn1[i]);
  }
}

void
gkyl_poly_roots_release(struct gkyl_poly_roots *pr)
{
  gkyl_free(pr->impart);
  gkyl_free(pr->rpart);
  gkyl_free(pr->err);
  gkyl_free(pr->work);
  gkyl_free(pr);
}
