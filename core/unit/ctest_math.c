#include <acutest.h>
#include <gkyl_math.h>

#include <complex.h>
#include <math.h>

// Table on Page 5 of van Engelen paper. See math.c for full reference
double func_1(double x, void *ctx) { return x*x*x-2*x*x+x; }

double func_2(double x, void *ctx) { return 1/(1+x); }

double func_12(double x, void *ctx) { return 1/sqrt(sin(M_PI*x)); }

double func_13(double x, void *ctx) { return pow(sin(M_PI*x), -0.8); }

double
func_cir(double x, void *ctx)
{
  return sqrt(1+x*x/(1-x*x));
}


static void
show_qr_res(struct gkyl_qr_res res, const char *msg)
{
  fprintf(stdout, "%s\n", msg);
  fprintf(stdout, ">> Status = %d. Res = %.15lg, Error = %g, Neval = %d. Nlevel = %d\n",
    res.status, res.res, res.error, res.nevals, res.nlevels);
}

void
test_dbl_exp(void)
{
  do {
    struct gkyl_qr_res res = gkyl_dbl_exp(func_1, 0, 0.0, 1.0, 10, 1e-11);
    //show_quad_res(res, "func_1");
    TEST_CHECK( gkyl_compare(res.res, 1.0/12.0, 1e-10) );
  } while(0);

  do {
    struct gkyl_qr_res res = gkyl_dbl_exp(func_2, 0, 0.0, 1.0, 10, 1e-11);
    //show_quad_res(res, "func_2");
    TEST_CHECK( gkyl_compare(res.res, log(2.0), 1e-10) );
  } while(0);

  do {
    // seems the following only gets to 1e-7 accuracy
    struct gkyl_qr_res res = gkyl_dbl_exp(func_12, 0, 0.0, 1.0, 10, 1e-15);
    //show_quad_res(res, "func_12");
    TEST_CHECK( gkyl_compare(res.res, 1.669253683348149, 1e-7) );
  } while(0);

  do {
    // The following integral 13 does not converge due to singularities
    // at the end point. van Engelen claims it does. I can't make it work.
    
    /* struct gkyl_quad_res res = gkyl_dbl_exp(func_13, 0, 0.0, 1.0, 10, 1e-11); */
    /* show_quad_res(res, "func_13"); */
    /* TEST_CHECK( gkyl_compare(res.quad_res, 3.604250526330095, 1e-10) ); */
  } while(0);

  do {
    struct gkyl_qr_res res = gkyl_dbl_exp(func_cir, 0, -1.0, 1.0, 10, 1e-16);
    //show_qr_res(res, "func_2");
    TEST_CHECK( gkyl_compare(res.res, M_PI, 1e-7) );
  } while(0);
}

double rfunc_1(double x, void *ctx)
{
  return x*x-1;
}

double rfunc_2(double x, void *ctx)
{
  // flat around the root x = 0
  return x*sin(x*x);
}

double rfunc_3(double x, void *ctx)
{
  return x*exp(x)-10;
}

void
test_ridders(void)
{
  do {
    double x1 = 0.5, x2 = 2.0;
    double f1 = rfunc_1(x1, 0), f2 = rfunc_1(x2, 0);
    struct gkyl_qr_res res = gkyl_ridders(rfunc_1, 0, x1, x2, f1, f2,
      100, 1e-12);
//    show_qr_res(res, "rfunc_1");
    TEST_CHECK( gkyl_compare(res.res, 1.0, 1e-10) );
  } while(0);

  do {
    double x1 = -2.0, x2 = 0.0;
    double f1 = rfunc_1(x1, 0), f2 = rfunc_1(x2, 0);
    struct gkyl_qr_res res = gkyl_ridders(rfunc_1, 0, x1, x2, f1, f2,
      100, 1e-12);
//    show_qr_res(res, "rfunc_1");
    TEST_CHECK( gkyl_compare(res.res, -1.0, 1e-10) );
  } while(0);  

  do {
    double x1 = -0.5, x2 = 1.0;
    double f1 = rfunc_2(x1, 0), f2 = rfunc_2(x2, 0);
    struct gkyl_qr_res res = gkyl_ridders(rfunc_2, 0, x1, x2, f1, f2,
      100, 1e-12);
//    show_qr_res(res, "rfunc_2");
    TEST_CHECK( gkyl_compare(res.res, 0.0, 1e-10) );
  } while(0);

  do {
    double x1 = -1.0, x2 = 5.0;
    double f1 = rfunc_3(x1, 0), f2 = rfunc_3(x2, 0);
    struct gkyl_qr_res res = gkyl_ridders(rfunc_3, 0, x1, x2, f1, f2,
      100, 1e-12);
//    show_qr_res(res, "rfunc_3");
    TEST_CHECK( gkyl_compare(res.res, 1.745528002740699, 1e-10) );
  } while(0);  
}

struct idx_status {
  int idx;
  bool status;
};

static struct idx_status
check_in_list(int nvals, const double complex *vals, double complex tocheck, double eps)
{
  for (int i=0; i<nvals; ++i)
    if (gkyl_compare_double(creal(vals[i]), creal(tocheck), eps) &&
      gkyl_compare_double(cimag(vals[i]), cimag(tocheck), eps) )
      return (struct idx_status) { .idx = i, .status = true };
  return (struct idx_status) { .idx = nvals, .status = false };
}

void
test_poly2_roots(void)
{
  struct gkyl_lo_poly_roots rts;

  do {
    static struct idx_status istat;
    
    double c1[4] = { 10.0, -7.0 };
    rts = gkyl_calc_lo_poly_roots(GKYL_LO_POLY_2, c1);
    // checks
    istat = check_in_list(2, (double complex[]) { 2.0, 5.0 },
      rts.rpart[0]+I*rts.impart[0],
      1e-15);
    TEST_CHECK( istat.status );
    
    istat = check_in_list(2, (double complex[]) { 2.0, 5.0 },
      rts.rpart[1]+I*rts.impart[1],
      1e-15);
    TEST_CHECK( istat.status );
  } while (0);

  do {
    struct idx_status istat;
    double c1[4] = { 13.0, -4.0 };
    rts = gkyl_calc_lo_poly_roots(GKYL_LO_POLY_2, c1);
    
    istat = check_in_list(2, (double complex[]) { 2.0-3.0*I, 2.0+3.0*I },
      rts.rpart[0] + I*rts.impart[0],
      1e-15);
    TEST_CHECK( istat.status );

    istat = check_in_list(2, (double complex[]) { 2.0-3.0*I, 2.0+3.0*I },
      rts.rpart[1] + I*rts.impart[1],
      1e-15);
    TEST_CHECK( istat.status );
  } while (0);
}

void
test_poly3_roots(void)
{
  struct gkyl_lo_poly_roots rts;

  do {
    static struct idx_status istat;

    double c1[4] = { 3.0, -5.5, -1.5 };
    rts = gkyl_calc_lo_poly_roots(GKYL_LO_POLY_3, c1);

    double complex res[4] = { 3.0, 0.5, -2.0 };

    for (int i=0; i<3; ++i)
      TEST_CHECK(
        check_in_list(3, res, rts.rpart[i]+I*rts.impart[i], 1e-14).status
      );
  } while (0);

  do {
    static struct idx_status istat;

    double c1[4] = { -7.5, 11.0, -5.5 };
    rts = gkyl_calc_lo_poly_roots(GKYL_LO_POLY_3, c1);

    double complex res[4] = { 1.5, 2+I, 2-I };

    for (int i=0; i<3; ++i)
      TEST_CHECK(
        check_in_list(3, res, rts.rpart[i]+I*rts.impart[i], 1e-14).status
      );
  } while (0);  

}

void
test_poly4_roots(void)
{
  struct gkyl_lo_poly_roots rts;

  do {
    double c1[4] = { 120.0, -26.0, -25.0, 2.0 };
    rts = gkyl_calc_lo_poly_roots(GKYL_LO_POLY_4, c1);

    double complex res[4] = { 4.0, 2.0, -3.0, -5.0 };

    for (int i=0; i<4; ++i)
      TEST_CHECK(
        check_in_list(4, res, rts.rpart[i]+I*rts.impart[i], 1e-14).status
      );
    
  } while (0);

  do {
    double c1[4] = { 520.0, 22.0, -3.0, 10.0 };
    rts = gkyl_calc_lo_poly_roots(GKYL_LO_POLY_4, c1);

    double complex res[4] = {
      -4.0, -10.0,
      2.0 + 3.0*I,
      2.0 - 3.0*I
    };

    for (int i=0; i<4; ++i)
      TEST_CHECK(
        check_in_list(4, res, rts.rpart[i]+I*rts.impart[i], 1e-14).status
      );
    
  } while (0);
}

void
test_polyn_roots(void)
{

  do {
    struct gkyl_poly_roots *rts = gkyl_poly_roots_new(4);
    
    double c1[] = { 120.0, -26.0, -25.0, 2.0 };
    gkyl_calc_poly_roots(rts, c1);

    double complex res[4] = { 4.0, 2.0, -3.0, -5.0 };

    for (int i=0; i<4; ++i)
      TEST_CHECK(
        check_in_list(4, res, rts->rpart[i]+I*rts->impart[i], 1e-14).status
      );

    gkyl_poly_roots_release(rts);
  } while (0);

  do {
    struct gkyl_poly_roots *rts = gkyl_poly_roots_new(5);
    
    double c1[] = { 14400.0, -11400.0, 1174.0, 711.0, -126.0 };
    gkyl_calc_poly_roots(rts, c1);

    double complex res[5] = { 120.0, 5.0, 3.0, 2.0, -4.0 };

    for (int i=0; i<5; ++i)
      TEST_CHECK(
        check_in_list(5, res, rts->rpart[i]+I*rts->impart[i], 1e-14).status
      );

    gkyl_poly_roots_release(rts);
  } while (0);  

}


void
test_sturn_root_intervals(void)
{

  // Test from wiki example: 2 real, distinct roots
  // https://en.wikipedia.org/wiki/Sturm%27s_theorem
  do {
    struct gkyl_root_intervals root_intervals; 

    // Setup the specific test
    double coeff[4] = {-1.0, -1.0, 0.0, 1.0};
    double domain[2] = {-3.0, 3.0};
    double tol = 1e-13;

    // compute root inverals
    root_intervals = gkyl_calc_quartic_root_intervals( coeff, domain, tol);

    // Check the outputs
    double lower[4] = {-3.0,0.0,0.0,0.0};
    double upper[4] = {0.0,3.0,0.0,0.0};
    TEST_CHECK( root_intervals.status == 0 );
    TEST_CHECK( root_intervals.nroots == 2 );
    TEST_CHECK( root_intervals.niter > 0 );
    for (int i=0; i<4; ++i){
      TEST_CHECK(gkyl_compare_double(upper[i], root_intervals.root_bound_upper[i], 1e-12));
      TEST_CHECK(gkyl_compare_double(lower[i], root_intervals.root_bound_lower[i], 1e-12));
    }

    // Compute the roots via ridders
    gkyl_root_isolation_from_intervals_via_ridders(&root_intervals, tol);

    // Check the outputs of the roots via ridders
    double roots[4] = {-1.0,1.0,0.0,0.0};
    for (int i=0; i<root_intervals.nroots; ++i){
      TEST_CHECK(gkyl_compare_double(roots[i], root_intervals.real_roots_ridders[i], 1e-12));
      TEST_CHECK( root_intervals.status_ridders[i] == 0 );
    }

    // test refined root intervals
    gkyl_refine_root_intervals_bisection(&root_intervals, tol);

    // Roots (maxima): [x=1.0,x=-1.0,x=-0.5*(1.732050807568877*%i+1.0),x=0.5*(1.732050807568877*%i-1.0)]

    // Check the outputs of the refinement pass
    double lower_refined[4] = {-1.0000000000000284e+00,9.9999999999994316e-01,0.0,0.0};
    double upper_refined[4] = {-9.9999999999994316e-01,1.0000000000000284e+00,0.0,0.0};
    for (int i=0; i<root_intervals.nroots; ++i){
      TEST_CHECK(gkyl_compare_double(upper_refined[i], root_intervals.root_bound_upper[i], 1e-12));
      TEST_CHECK(gkyl_compare_double(lower_refined[i], root_intervals.root_bound_lower[i], 1e-12));
      TEST_CHECK( root_intervals.status_refinement[i] == 0 );
      TEST_CHECK( root_intervals.niter_refinement[i] > 0 );
    }
  } while (0);  

  // Test: 4 real roots
  do {
    struct gkyl_root_intervals root_intervals; 

    // Setup the specific test
    double coeff[4] = {0.1000, 0.0, -1.0000, 0.0};
    double domain[2] = {-3.0, 3.0};
    double tol = 1e-13;

    // compute root inverals
    root_intervals = gkyl_calc_quartic_root_intervals( coeff, domain, tol);

    // Check the outputs
    double lower[4] = {-1.5,-0.75,0.0,0.75};
    double upper[4] = {-0.75,0.0,0.75,1.5};
    TEST_CHECK( root_intervals.status == 0 );
    TEST_CHECK( root_intervals.nroots == 4 );
    TEST_CHECK( root_intervals.niter > 0 );
    for (int i=0; i<4; ++i){
      TEST_CHECK(gkyl_compare_double(upper[i], root_intervals.root_bound_upper[i], 1e-12));
      TEST_CHECK(gkyl_compare_double(lower[i], root_intervals.root_bound_lower[i], 1e-12));
    }

    // Compute the roots via ridders
    gkyl_root_isolation_from_intervals_via_ridders(&root_intervals, tol);

    // Check the outputs of the roots via ridders
    double roots[4] = {-0.9419651451198933,-0.3357106870197288,0.3357106870197288,0.9419651451198933};
    for (int i=0; i<root_intervals.nroots; ++i){
      TEST_CHECK(gkyl_compare_double(roots[i], root_intervals.real_roots_ridders[i], 1e-12));
      TEST_CHECK( root_intervals.status_ridders[i] == 0 );
    }

    // test refined root intervals
    gkyl_refine_root_intervals_bisection(&root_intervals, tol);

    // Maxima roots: 	[x=-0.9419651451198933,x=0.9419651451198933,x=-0.3357106870197288,x=0.3357106870197288]

    // Check the outputs of the refinement pass
    double lower_refined[4] = {-9.4196514511995133e-01,-3.3571068701979812e-01,3.3571068701971285e-01,9.4196514511986607e-01};
    double upper_refined[4] = {-9.4196514511986607e-01,-3.3571068701971285e-01,3.3571068701979812e-01,9.4196514511995133e-01};
    for (int i=0; i<root_intervals.nroots; ++i){
      TEST_CHECK(gkyl_compare_double(upper_refined[i], root_intervals.root_bound_upper[i], 1e-12));
      TEST_CHECK(gkyl_compare_double(lower_refined[i], root_intervals.root_bound_lower[i], 1e-12));
      TEST_CHECK( root_intervals.status_refinement[i] == 0 );
      TEST_CHECK( root_intervals.niter_refinement[i] > 0 );
    }
  } while (0);  


  // Test: 4 real roots, more complex polynomial
  do {
    struct gkyl_root_intervals root_intervals; 

    // Setup the specific test
    double coeff[4] = {-0.5170, 1.2377, 0.0354, -1.7561};
    double domain[2] = {-3.0, 3.0};
    double tol = 1e-13;

    // compute root inverals
    root_intervals = gkyl_calc_quartic_root_intervals( coeff, domain, tol);

    // Check the outputs
    double lower[4] = {-3.0,0.0,0.75,0.9375};
    double upper[4] = {0,0.75,0.9375,1.1250};
    TEST_CHECK( root_intervals.status == 0 );
    TEST_CHECK( root_intervals.nroots == 4 );
    TEST_CHECK( root_intervals.niter > 0 );
    for (int i=0; i<4; ++i){
      TEST_CHECK(gkyl_compare_double(upper[i], root_intervals.root_bound_upper[i], 1e-12));
      TEST_CHECK(gkyl_compare_double(lower[i], root_intervals.root_bound_lower[i], 1e-12));
    }

    // Compute the roots via ridders
    gkyl_root_isolation_from_intervals_via_ridders(&root_intervals, tol);

    // Check the outputs of the roots via ridders
    double roots[4] = {-8.3856473883897131e-01,6.5873479706025262e-01,9.3592994177871747e-01,1.0000000000000040e+00};
    for (int i=0; i<root_intervals.nroots; ++i){
      TEST_CHECK(gkyl_compare_double(roots[i], root_intervals.real_roots_ridders[i], 1e-12));
      TEST_CHECK( root_intervals.status_ridders[i] == 0 );
    }

    // test refined root intervals
    gkyl_refine_root_intervals_bisection(&root_intervals, tol);

    // Maxima roots: 	... fails to simplify expr, verified by plotting / Matlab testcode

    // Check the outputs of the refinement pass
    double lower_refined[4] = {-8.3856473883901117e-01,6.5873479706024796e-01,9.3592994177868438e-01,9.9999999999994316e-01};
    double upper_refined[4] = {-8.3856473883892590e-01,6.5873479706033322e-01,9.3592994177874400e-01,1.0000000000000284e+00};
    for (int i=0; i<root_intervals.nroots; ++i){
      TEST_CHECK(gkyl_compare_double(upper_refined[i], root_intervals.root_bound_upper[i], 1e-12));
      TEST_CHECK(gkyl_compare_double(lower_refined[i], root_intervals.root_bound_lower[i], 1e-12));
      TEST_CHECK( root_intervals.status_refinement[i] == 0 );
      TEST_CHECK( root_intervals.niter_refinement[i] > 0 );
    }

    // Compute the roots via ridders using the refined domains
    gkyl_root_isolation_from_intervals_via_ridders(&root_intervals, tol);

    // Check the outputs of the roots via ridders
    for (int i=0; i<root_intervals.nroots; ++i){
      TEST_CHECK(gkyl_compare_double(roots[i], root_intervals.real_roots_ridders[i], 1e-12));
      TEST_CHECK( root_intervals.status_ridders[i] == 0 );
    }
  } while (0); 


  // Test: 3-distinct roots, search interval falls on x = 0 root
  do {
    struct gkyl_root_intervals root_intervals; 

    // Setup the specific test
    double coeff[4] = {0.0, 0.0, -1.0, 0.0};
    double domain[2] = {-3.0, 3.0};
    double tol = 1e-13;

    // compute root inverals
    root_intervals = gkyl_calc_quartic_root_intervals( coeff, domain, tol);

    // Check the outputs
    double lower[4] = {0.6,-1.2,-0.3,0.0};
    double upper[4] = {3.0,-0.3,0.6,0.0};
    TEST_CHECK( root_intervals.status == 0 );
    TEST_CHECK( root_intervals.nroots == 3 );
    TEST_CHECK( root_intervals.niter > 0 );
    for (int i=0; i<4; ++i){
      TEST_CHECK(gkyl_compare_double(upper[i], root_intervals.root_bound_upper[i], 1e-12));
      TEST_CHECK(gkyl_compare_double(lower[i], root_intervals.root_bound_lower[i], 1e-12));
    }

    // Compute the roots via ridders using the refined domains
    gkyl_root_isolation_from_intervals_via_ridders(&root_intervals, tol);

    // Ridders fails to find the repeated root! (root 3), returns status_ridders == 1 

    // Check the outputs of the roots via ridders
    double roots[4] = {1.0,-1.0,0.0,0.0};
    for (int i=0; i<2; ++i){
      TEST_CHECK(gkyl_compare_double(roots[i], root_intervals.real_roots_ridders[i], 1e-12));
      TEST_CHECK( root_intervals.status_ridders[i] == 0 );
    }

   // test refined root intervals
    gkyl_refine_root_intervals_bisection(&root_intervals, tol);

    // Roots: [0,0,+1,-1]

    // Check the outputs of the refinement pass
    double lower_refined[4] = {9.9999999999997746e-01,-1.0000000000000455e+00,-1.6979010789934061e-14,0.0};
    double upper_refined[4] = {1.0000000000000457e+00,-9.9999999999999434e-01,3.4180066184793154e-14,0.0};
    for (int i=0; i<root_intervals.nroots; ++i){
      TEST_CHECK(gkyl_compare_double(upper_refined[i], root_intervals.root_bound_upper[i], 1e-12));
      TEST_CHECK(gkyl_compare_double(lower_refined[i], root_intervals.root_bound_lower[i], 1e-12));
      TEST_CHECK( root_intervals.status_refinement[i] == 0 );
      TEST_CHECK( root_intervals.niter_refinement[i] > 0 );
    }
  } while (0); 


 // Test: 3-distinct roots, large ordering
 // Maxima res: [x=-1.0*10^-10,x=1.0*10^-10,x=0.0] 
  do {
    struct gkyl_root_intervals root_intervals; 

    // Setup the specific test
    double coeff[4] = {0.0, 0.0, -1.0e-20, 0.0};
    double domain[2] = {-3.0, 3.0};
    double tol = 1e-13;

    // compute root inverals
    root_intervals = gkyl_calc_quartic_root_intervals( coeff, domain, tol);

    // Check the outputs
    double lower[4] = {-2.793966983697753e-10,-6.984911908129258e-11,3.492467056294875e-11,0.0};
    double upper[4] = {-6.984911908129258e-11,3.492467056294875e-11,1.396984602071901e-10,0.0};
    TEST_CHECK( root_intervals.status == 0 );
    TEST_CHECK( root_intervals.nroots == 3 );
    TEST_CHECK( root_intervals.niter > 0 );
    for (int i=0; i<4; ++i){
      TEST_CHECK(gkyl_compare_double(upper[i], root_intervals.root_bound_upper[i], 1e-12));
      TEST_CHECK(gkyl_compare_double(lower[i], root_intervals.root_bound_lower[i], 1e-12));
    }

       // Compute the roots via ridders using the refined domains
    gkyl_root_isolation_from_intervals_via_ridders(&root_intervals, tol);

    // Ridders fails to find the repeated root! (root 2), returns status_ridders == 1 

    // Check the outputs of the roots via ridders
    double roots[4] = {-1.0*10e-10,0.0,1.0*10e-10,0.0};
    for (int i=0; i<root_intervals.nroots; ++i){
      if (i != 1){
        TEST_CHECK(gkyl_compare_double(roots[i], root_intervals.real_roots_ridders[i], 1e-2));
        TEST_CHECK( root_intervals.status_ridders[i] == 0 );
      }
    }

    // test refined root intervals
    gkyl_refine_root_intervals_bisection(&root_intervals, tol);

    // Maxima res: [x=-1.0*10^-10,x=1.0*10^-10,x=0.0] 

    // Check the outputs of the refinement pass
    double lower_refined[4] = {-1.0003297449638165e-10,-1.6979010789934061e-14,9.9999016474801779e-11,0.0};
    double upper_refined[4] = {-9.9981815419406931e-11,3.4180066184793154e-14,1.0005017555177650e-10,0.0};
    for (int i=0; i<root_intervals.nroots; ++i){
      TEST_CHECK(gkyl_compare_double(upper_refined[i], root_intervals.root_bound_upper[i], 1e-12));
      TEST_CHECK(gkyl_compare_double(lower_refined[i], root_intervals.root_bound_lower[i], 1e-12));
      TEST_CHECK( root_intervals.status_refinement[i] == 0 );
      TEST_CHECK( root_intervals.niter_refinement[i] > 0 );
    }
  } while (0); 


  // Test: 3-distinct roots, shifted (x+1.1)^4 - (x+1.1)^2 = 0
  // Maxima res: [x=-2.1,x=-0.1,x=-1.1]
  // FAILS TO FIND REPEATED ROOT AT x = -0.1!
  do {
    struct gkyl_root_intervals root_intervals; 

    // Setup the specific test
    double coeff[4] = {0.2541, 3.1240, 6.2600, 4.4000};
    double domain[2] = {-3.0, 3.0};
    double tol = 1e-13;

    // compute root inverals
    root_intervals = gkyl_calc_quartic_root_intervals( coeff, domain, tol);

    // Check we have the right number of roots etc
    //printf("\nnum-roots: %d\n",root_intervals.nroots);
    //printf("num-iterations: %d\n",root_intervals.niter );
    //printf("Status: %d\n",root_intervals.status );
    //for (int i=0; i<4; ++i){
    //  printf("Root bounds: [L,R]: [%1.16e,%1.16e]\n",root_intervals.root_bound_lower[i],
    //  root_intervals.root_bound_upper[i]);
    //}

    // Check the outputs
    double lower[4] = {-3.0,-1.5,0.0,0.0};
    double upper[4] = {-1.5,0.0,0.0,0.0};
    TEST_CHECK( root_intervals.status == 0 );
    TEST_CHECK( root_intervals.nroots == 2 );
    TEST_CHECK( root_intervals.niter > 0 );
    for (int i=0; i<4; ++i){
      TEST_CHECK(gkyl_compare_double(upper[i], root_intervals.root_bound_upper[i], 1e-12));
      TEST_CHECK(gkyl_compare_double(lower[i], root_intervals.root_bound_lower[i], 1e-12));
    }

  // Compute the roots via ridders using the refined domains
    gkyl_root_isolation_from_intervals_via_ridders(&root_intervals, tol);

    // Check we have the right number of roots etc
    //printf("\nnum-roots: %d\n",root_intervals.nroots);
    //for (int i=0; i<root_intervals.nroots; ++i){
    //  printf("Root bounds: [L,R]: [%1.16e,%1.16e]\n",root_intervals.root_bound_lower[i],
    //  root_intervals.root_bound_upper[i]);
    //  printf("Roots: [%1.16e]\n",root_intervals.real_roots_ridders[i] );
    //  printf("num-iterations: %d\n",root_intervals.niter_refinement[i] );
    //  printf("Status Ridders: %d\n",root_intervals.status_ridders[i] );
    //}

    // Check the outputs of the roots via ridders
    double roots[4] = {-2.1,-0.1,0.0,0.0};
    for (int i=0; i<root_intervals.nroots; ++i){
      TEST_CHECK(gkyl_compare_double(roots[i], root_intervals.real_roots_ridders[i], 1e-12));
      TEST_CHECK( root_intervals.status_ridders[i] == 0 );
    }

    // test refined root intervals
    gkyl_refine_root_intervals_bisection(&root_intervals, tol);

    // Check we have the right number of roots etc
    //printf("\nnum-roots: %d\n",root_intervals.nroots);
    //for (int i=0; i<root_intervals.nroots; ++i){
    //  printf("Root bounds: [L,R]: [%1.16e,%1.16e]\n",root_intervals.root_bound_lower[i],
    //  root_intervals.root_bound_upper[i]);
    //  printf("num-iterations: %d\n",root_intervals.niter_refinement[i] );
    //  printf("Status: %d\n",root_intervals.status_refinement[i] );
    //}

    // Maxima res: [x=-2.1,x=-0.1,x=-1.1 (x2)]

    // Check the outputs of the refinement pass
    double lower_refined[4] = {-2.1000000000000512e+00,-1.0000000000007958e-01,0.0,0.0};
    double upper_refined[4] = {-2.0999999999999659e+00,-9.9999999999994316e-02,0.0,0.0};
    for (int i=0; i<root_intervals.nroots; ++i){
      TEST_CHECK(gkyl_compare_double(upper_refined[i], root_intervals.root_bound_upper[i], 1e-12));
      TEST_CHECK(gkyl_compare_double(lower_refined[i], root_intervals.root_bound_lower[i], 1e-12));
      TEST_CHECK( root_intervals.status_refinement[i] == 0 );
      TEST_CHECK( root_intervals.niter_refinement[i] > 0 );
    }
  } while (0);

}

TEST_LIST = {
  { "dbl_exp", test_dbl_exp },
  { "ridders", test_ridders },
  { "poly2_roots", test_poly2_roots },
  { "poly3_roots", test_poly3_roots },
  { "poly4_roots", test_poly4_roots },
  { "strun_root_intervals", test_sturn_root_intervals },
  { "polyn_roots", test_polyn_roots },  
  { NULL, NULL },
};
