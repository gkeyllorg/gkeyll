#include <acutest.h>

#include <stdio.h>

#include <gkyl_dual_num.h>
#include <gkyl_math.h>

static inline bool
cmp_dn(struct gkyl_dn x, struct gkyl_dn xcmp)
{
  return (x.x[0] == xcmp.x[0]) && (x.x[1] == xcmp.x[1]);
}

static inline bool
cmp_dn2(struct gkyl_dn2 x, struct gkyl_dn2 xcmp)
{
  return (x.x[0] == xcmp.x[0]) && (x.x[1] == xcmp.x[1]) && (x.x[2] == xcmp.x[2]);
}

static void
show_dn(FILE *fp, const char *msg, struct gkyl_dn d1)
{
  fprintf(fp, "%s: [%lg, %lg]\n", msg, d1.x[0], d1.x[1]);
}

static inline
struct gkyl_dn
func_1(struct gkyl_dn x)
{
  // -2*x + x^2/(3+x^3) + 1/x
  return gdn_add(
    gdn_smul(-2.0, x),
    gdn_add(
      gdn_div(gdn_sq(x), gdn_sadd(3, gdn_cube(x))),
      gdn_inv(x)
    )
  );
}

double func_1_0(double x)
{
  return -2*x+x*x/(3+x*x*x) + 1/x;
}

double func_1_1(double x)
{
  return 2*x/(3+x*x*x) - 3*x*x*x*x/((3+x*x*x)*(3+x*x*x)) - 1/(x*x) - 2;
}

void test_basic(void)
{
  double x10 = 2.5;
  struct gkyl_dn x1 = gdn_new1(x10);

  struct gkyl_dn res = { };

  res = gdn_sq(x1);
  TEST_CHECK( cmp_dn(res, gdn_new(x10*x10, 2*x10)) );

  res = gdn_cube(x1);
  TEST_CHECK( cmp_dn(res, gdn_new(x10*x10*x10, 3*x10*x10)) );

  res = gdn_npow(x1, 5);
  TEST_CHECK( cmp_dn(res, gdn_new(pow(x10,5), 5*pow(x10,4))) );

  res = gdn_sqrt(x1);
  TEST_CHECK( cmp_dn(res, gdn_new(sqrt(x10), 0.5/sqrt(x10))) );

  res = gdn_cos(x1);
  TEST_CHECK( cmp_dn(res, gdn_new(cos(x10), -sin(x10))) );

  res = gdn_sin(x1);
  TEST_CHECK( cmp_dn(res, gdn_new(sin(x10), cos(x10))) );

  res = gdn_tan(x1);
  TEST_CHECK( res.x[0] == tan(x10) );
  TEST_CHECK( gkyl_compare_double(res.x[1], 1/(cos(x10)*cos(x10)), 1e-15) );

  res = gdn_log(x1);
  TEST_CHECK( res.x[0] == log(x10) );
  TEST_CHECK( gkyl_compare_double(res.x[1], 1/x10, 1e-15) );

  res = func_1(x1);
  TEST_CHECK( cmp_dn(res, gdn_new(func_1_0(x10), func_1_1(x10))) );

}

void test_basic2(void)
{
  double x10 = 2.5;
  struct gkyl_dn2 x1 = gdn2_new(x10, 1.0, 2.0);

  struct gkyl_dn2 res = { };

  res = gdn2_sq(x1);
  TEST_CHECK( cmp_dn2(res, gdn2_new(x10*x10, 2*x10, 2.0*2*x10)) );

  res = gdn2_cube(x1);
  TEST_CHECK( cmp_dn2(res, gdn2_new(x10*x10*x10, 3*x10*x10, 2.0*3*x10*x10)) );

  res = gdn2_npow(x1, 5);
  TEST_CHECK( cmp_dn2(res, gdn2_new(pow(x10,5), 5*pow(x10,4), 2.0*5*pow(x10,4))) );

  res = gdn2_sqrt(x1);
  TEST_CHECK( cmp_dn2(res, gdn2_new(sqrt(x10), 0.5/sqrt(x10), 2.0*0.5/sqrt(x10))) );

  res = gdn2_cos(x1);
  TEST_CHECK( cmp_dn2(res, gdn2_new(cos(x10), -sin(x10), -2.0*sin(x10))) );

  res = gdn2_sin(x1);
  TEST_CHECK( cmp_dn2(res, gdn2_new(sin(x10), cos(x10), 2.0*cos(x10))) );

  res = gdn2_log(x1);
  TEST_CHECK( cmp_dn2(res, gdn2_new(log(x10), 1/x10, 2.0/x10)) );  

  res = gdn2_tan(x1);
  TEST_CHECK( res.x[0] == tan(x10) );
  TEST_CHECK( gkyl_compare_double(res.x[1], 1/(cos(x10)*cos(x10)), 1e-15) );
  TEST_CHECK( gkyl_compare_double(res.x[2], 2.0/(cos(x10)*cos(x10)), 1e-15)  );
}

static struct gkyl_dn2
fxy( struct gkyl_dn2 x, struct gkyl_dn2 y)
{
  // cos(x/y)
  return gdn2_cos( gdn2_div(x,y) );
}

static struct gkyl_dn2
psixy( struct gkyl_dn2 x, struct gkyl_dn2 y)
{
  // sqrt(x^2 + y^2)
  return gdn2_sqrt( gdn2_add(gdn2_sq(x), gdn2_sq(y)) );
}

static struct gkyl_dn2
custom_xy( struct gkyl_dn2 x, struct gkyl_dn2 y)
{
  // example of custom function of x,y with hand-computed gradient
  
  double g0 = (x.x[0]*x.x[0])/(y.x[0]*y.x[0]);
  double g0x = 2.0*x.x[0]/(y.x[0]*y.x[0]);
  double g0y = -2.0*(x.x[0]*x.x[0])/(y.x[0]*y.x[0]*y.x[0]);
  
  return (struct gkyl_dn2) {
    g0,
    x.x[1]*g0x + y.x[1]*g0y,
    x.x[2]*g0x + y.x[2]*g0y
  };
}

void test_xy(void)
{
  struct gkyl_dn2 x = gdn2_new10(2.5), y = gdn2_new01(1.5);
  struct gkyl_dn2 res = { };

  res = fxy(x, y);
  TEST_CHECK( gkyl_compare_double(res.x[0], cos(x.x[0]/y.x[0]), 1e-14) );
  TEST_CHECK( gkyl_compare_double(res.x[1], -sin(x.x[0]/y.x[0])/y.x[0], 1e-14) );
  TEST_CHECK( gkyl_compare_double(res.x[2], x.x[0]*sin(x.x[0]/y.x[0])/(y.x[0]*y.x[0]), 1e-14) );

  double psi0 = sqrt(x.x[0]*x.x[0] + y.x[0]*y.x[0]);
  res = psixy(x, y);
  TEST_CHECK( gkyl_compare_double(res.x[0], psi0, 1e-14) );
  TEST_CHECK( gkyl_compare_double(res.x[1], x.x[0]/psi0, 1e-14) );
  TEST_CHECK( gkyl_compare_double(res.x[2], y.x[0]/psi0, 1e-14) );

  res = gdn2_cos( custom_xy(x, y) ); // cos of custom defined function
  
  double cus0 = cos(x.x[0]*x.x[0]/(y.x[0]*y.x[0]));
  double cusx = -2.0*x.x[0]*sin(x.x[0]*x.x[0]/(y.x[0]*y.x[0]))/(y.x[0]*y.x[0]);
  double cusy = 2.0*x.x[0]*x.x[0]*sin(x.x[0]*x.x[0]/(y.x[0]*y.x[0]))/(y.x[0]*y.x[0]*y.x[0]);
  TEST_CHECK( gkyl_compare_double(res.x[0], cus0, 1e-14) );
  TEST_CHECK( gkyl_compare_double(res.x[1], cusx, 1e-14) );
  TEST_CHECK( gkyl_compare_double(res.x[2], cusy, 1e-14) );
}

// for testing 1D mapc2p determine by its inverse mapping
struct inv_map_ctx { double s0; };

static inline struct gkyl_dn
gdn_inv_mapc2p(struct gkyl_dn x, void *ctx)
{
  struct inv_map_ctx *imctx = ctx;
  double s0 = imctx->s0;
  // sqrt(x+4) - 2 - s0
  return gdn_sadd(-2.0-s0, gdn_sqrt(gdn_sadd(4.0, x)));
}

static inline double
inv_mapc2p(double x, void *ctx)
{
  return gdn_inv_mapc2p(gdn_new0(x),ctx).x[0];
}

void test_inv_mapc2p(void)
{
  struct inv_map_ctx imctx = { 0.75 };
  double xl = 0.0, xr = 5.0;
  double fl = inv_mapc2p(xl, &imctx);
  double fr = inv_mapc2p(xr, &imctx);

  struct gkyl_qr_res root =
    gkyl_ridders(inv_mapc2p, &imctx, xl, xr, fl, fr, 100, 1e-10);

  // compute derivative at the root
  struct gkyl_dn res = gdn_inv_mapc2p(gdn_new1(root.res), &imctx);
  // at this point 1/res.x[1] is the derivative of the mapping wrt s
  // at the root
  
  TEST_CHECK( gkyl_compare_double(root.res, 3.5625, 1e-14) );
  TEST_CHECK( gkyl_compare_double(1/res.x[1], 5.5, 1e-14) );
}

static void
mapc2p(const struct gkyl_dn xc[2], struct gkyl_dn xp[2])
{
  // mapping from r,theta -> x,y
  xp[0] = gdn_mul(xc[0], gdn_cos(xc[1]));
  xp[1] = gdn_mul(xc[0], gdn_sin(xc[1]));
}

void test_mapc2p(void)
{
  double r = 1.5, theta = M_PI/3;
  struct gkyl_dn xc[2], xp[2];

  // compute gradient wrt r (notice new0 used for theta)
  xc[0] = gdn_new1(r); xc[1] = gdn_new0(theta);
  mapc2p(xc, xp);

  TEST_CHECK( xp[0].x[0] == r*cos(theta) );
  TEST_CHECK( xp[1].x[0] == r*sin(theta) );

  TEST_CHECK( xp[0].x[1] == cos(theta) );
  TEST_CHECK( xp[1].x[1] == sin(theta) );

  // compute gradient wrt theta (notice new0 used for r)
  xc[0] = gdn_new0(r); xc[1] = gdn_new1(theta);
  mapc2p(xc, xp);

  TEST_CHECK( xp[0].x[0] == r*cos(theta) );
  TEST_CHECK( xp[1].x[0] == r*sin(theta) );

  TEST_CHECK( xp[0].x[1] == -r*sin(theta) );
  TEST_CHECK( xp[1].x[1] == r*cos(theta) );
}

static void
mapc2p_2(const struct gkyl_dn2 xc[2], struct gkyl_dn2 xp[2])
{
  // mapping from r,theta -> x,y
  xp[0] = gdn2_mul(xc[0], gdn2_cos(xc[1]));
  xp[1] = gdn2_mul(xc[0], gdn2_sin(xc[1]));
}

void test_mapc2p_2(void)
{
  double r = 1.5, theta = M_PI/3;
  struct gkyl_dn2 xc[2], xp[2];

  // compute function and gradients
  xc[0] = gdn2_new10(r); xc[1] = gdn2_new01(theta);
  mapc2p_2(xc, xp);

  TEST_CHECK( xp[0].x[0] == r*cos(theta) );
  TEST_CHECK( xp[0].x[1] == cos(theta) );
  TEST_CHECK( xp[0].x[2] == -r*sin(theta) );

  TEST_CHECK( xp[1].x[0] == r*sin(theta) );
  TEST_CHECK( xp[1].x[1] == sin(theta) );
  TEST_CHECK( xp[1].x[2] == r*cos(theta) );
}


static inline double sq(double x) { return x * x; }

static inline struct gkyl_dn2
RpsiZ_ellip(const struct gkyl_dn2 psiZ[2])
{
  // psi*sin(Z)
  return gdn2_mul(psiZ[0], gdn2_sin(psiZ[1]));
}

/* static inline struct gkyl_dn2 */
/* dRdZ_ellip(const struct gkyl_dn2 psiZ[2]) */
/* { */
/*   struct gkyl_dn2 RpsiZ = RpsiZ_ellip(psiZ); */
/*   double dRdZ = RpsiZ.x[2]; */
/* } */

void
test_psi_mapping(void)
{
  double pz[2] = { 1.0, 2.0 };
  struct gkyl_dn2 psiZ[2] = { gdn2_new10(pz[0]), gdn2_new01(pz[1]) };
  struct gkyl_dn2 RpsiZ = RpsiZ_ellip(psiZ);

  TEST_CHECK( RpsiZ.x[0] == pz[0]*sin(pz[1]) );
  TEST_CHECK( RpsiZ.x[1] == sin(pz[1]) ); // dR/dpsi
  TEST_CHECK( RpsiZ.x[2] == pz[0]*cos(pz[1]) ); // dR/dZ
}

TEST_LIST = {
  { "test_basic", test_basic },
  { "test_basic2", test_basic2 },
  { "test_xy", test_xy },
  { "test_inv_mapc2p", test_inv_mapc2p },
  { "test_mapc2p", test_mapc2p },
  { "test_mapc2p_2", test_mapc2p_2 },
  { "test_psi_mapping", test_psi_mapping  },
  { NULL, NULL },  
};
