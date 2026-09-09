// Unit tests for the inline vector and limiter helpers in gkyl_math.h
#include <acutest.h>
#include <gkyl_math.h>
#include <gkyl_util.h>
#include <math.h>

void
test_vec3_basic()
{
  struct gkyl_vec3 z = gkyl_vec3_zeros();
  TEST_CHECK( z.x[0] == 0.0 );
  TEST_CHECK( z.x[1] == 0.0 );
  TEST_CHECK( z.x[2] == 0.0 );

  struct gkyl_vec3 a = gkyl_vec3_new(1.0, 2.0, 3.0);
  TEST_CHECK( a.x[0] == 1.0 );
  TEST_CHECK( a.x[1] == 2.0 );
  TEST_CHECK( a.x[2] == 3.0 );
}

void
test_vec3_scale()
{
  struct gkyl_vec3 a = gkyl_vec3_new(1.0, -2.0, 3.5);
  struct gkyl_vec3 b = gkyl_vec3_scale(2.0, a);
  TEST_CHECK( b.x[0] == 2.0 );
  TEST_CHECK( b.x[1] == -4.0 );
  TEST_CHECK( b.x[2] == 7.0 );

  struct gkyl_vec3 c = gkyl_vec3_scale(0.0, a);
  TEST_CHECK( c.x[0] == 0.0 && c.x[1] == 0.0 && c.x[2] == 0.0 );
}

void
test_vec3_add_sub()
{
  struct gkyl_vec3 a = gkyl_vec3_new(1.0, 2.0, 3.0);
  struct gkyl_vec3 b = gkyl_vec3_new(4.0, -1.0, 0.5);

  struct gkyl_vec3 s = gkyl_vec3_add(a, b);
  TEST_CHECK( s.x[0] == 5.0 );
  TEST_CHECK( s.x[1] == 1.0 );
  TEST_CHECK( s.x[2] == 3.5 );

  struct gkyl_vec3 d = gkyl_vec3_sub(a, b);
  TEST_CHECK( d.x[0] == -3.0 );
  TEST_CHECK( d.x[1] == 3.0 );
  TEST_CHECK( d.x[2] == 2.5 );

  // a - a = 0
  struct gkyl_vec3 z = gkyl_vec3_sub(a, a);
  TEST_CHECK( z.x[0] == 0.0 && z.x[1] == 0.0 && z.x[2] == 0.0 );
}

void
test_vec3_len_norm()
{
  struct gkyl_vec3 a = gkyl_vec3_new(3.0, 4.0, 0.0);
  TEST_CHECK( gkyl_compare_double(gkyl_vec3_len(a), 5.0, 1e-15) );

  struct gkyl_vec3 n = gkyl_vec3_norm(a);
  TEST_CHECK( gkyl_compare_double(gkyl_vec3_len(n), 1.0, 1e-15) );
  TEST_CHECK( gkyl_compare_double(n.x[0], 0.6, 1e-15) );
  TEST_CHECK( gkyl_compare_double(n.x[1], 0.8, 1e-15) );

  struct gkyl_vec3 b = gkyl_vec3_new(1.0, 2.0, 2.0);
  TEST_CHECK( gkyl_compare_double(gkyl_vec3_len(b), 3.0, 1e-15) );
}

void
test_vec3_dot()
{
  struct gkyl_vec3 a = gkyl_vec3_new(1.0, 2.0, 3.0);
  struct gkyl_vec3 b = gkyl_vec3_new(4.0, 5.0, 6.0);
  TEST_CHECK( gkyl_compare_double(gkyl_vec3_dot(a, b), 32.0, 1e-15) );

  // orthogonal vectors
  struct gkyl_vec3 e1 = gkyl_vec3_new(1.0, 0.0, 0.0);
  struct gkyl_vec3 e2 = gkyl_vec3_new(0.0, 1.0, 0.0);
  TEST_CHECK( gkyl_vec3_dot(e1, e2) == 0.0 );

  // dot with self == len^2
  double l = gkyl_vec3_len(a);
  TEST_CHECK( gkyl_compare_double(gkyl_vec3_dot(a, a), l*l, 1e-13) );
}

void
test_vec3_cross()
{
  struct gkyl_vec3 e1 = gkyl_vec3_new(1.0, 0.0, 0.0);
  struct gkyl_vec3 e2 = gkyl_vec3_new(0.0, 1.0, 0.0);
  struct gkyl_vec3 e3 = gkyl_vec3_cross(e1, e2);
  TEST_CHECK( gkyl_compare_double(e3.x[0], 0.0, 1e-15) );
  TEST_CHECK( gkyl_compare_double(e3.x[1], 0.0, 1e-15) );
  TEST_CHECK( gkyl_compare_double(e3.x[2], 1.0, 1e-15) );

  // anti-commutativity: a x b = -(b x a)
  struct gkyl_vec3 a = gkyl_vec3_new(2.0, -1.0, 3.0);
  struct gkyl_vec3 b = gkyl_vec3_new(0.5, 4.0, -2.0);
  struct gkyl_vec3 ab = gkyl_vec3_cross(a, b);
  struct gkyl_vec3 ba = gkyl_vec3_cross(b, a);
  TEST_CHECK( gkyl_compare_double(ab.x[0], -ba.x[0], 1e-14) );
  TEST_CHECK( gkyl_compare_double(ab.x[1], -ba.x[1], 1e-14) );
  TEST_CHECK( gkyl_compare_double(ab.x[2], -ba.x[2], 1e-14) );

  // cross product is orthogonal to both inputs
  TEST_CHECK( gkyl_compare_double(gkyl_vec3_dot(ab, a), 0.0, 1e-13) );
  TEST_CHECK( gkyl_compare_double(gkyl_vec3_dot(ab, b), 0.0, 1e-13) );

  // a x a = 0
  struct gkyl_vec3 z = gkyl_vec3_cross(a, a);
  TEST_CHECK( gkyl_compare_double(gkyl_vec3_len(z), 0.0, 1e-14) );
}

void
test_vec3_triple()
{
  struct gkyl_vec3 e1 = gkyl_vec3_new(1.0, 0.0, 0.0);
  struct gkyl_vec3 e2 = gkyl_vec3_new(0.0, 1.0, 0.0);
  struct gkyl_vec3 e3 = gkyl_vec3_new(0.0, 0.0, 1.0);
  // triple product of unit basis is the determinant == 1
  TEST_CHECK( gkyl_compare_double(gkyl_vec3_triple(e1, e2, e3), 1.0, 1e-15) );

  // general: a.(b x c) equals scalar triple via determinant
  struct gkyl_vec3 a = gkyl_vec3_new(1.0, 2.0, 3.0);
  struct gkyl_vec3 b = gkyl_vec3_new(4.0, 5.0, 6.0);
  struct gkyl_vec3 c = gkyl_vec3_new(7.0, 8.0, 10.0);
  double det = 1.0*(5.0*10.0-6.0*8.0) - 2.0*(4.0*10.0-6.0*7.0) + 3.0*(4.0*8.0-5.0*7.0);
  TEST_CHECK( gkyl_compare_double(gkyl_vec3_triple(a, b, c), det, 1e-12) );

  // coplanar vectors give zero triple product
  TEST_CHECK( gkyl_compare_double(gkyl_vec3_triple(e1, e2, gkyl_vec3_add(e1,e2)), 0.0, 1e-14) );
}

void
test_vec3_polar()
{
  // At phi=0 the contravariant->cartesian transform should be near-identity in x.
  struct gkyl_vec3 pin = gkyl_vec3_new(2.0, 0.0, 5.0);
  struct gkyl_vec3 out = gkyl_vec3_polar_con_to_cart(1.0, 0.0, pin);
  TEST_CHECK( gkyl_compare_double(out.x[0], 2.0, 1e-14) );
  TEST_CHECK( gkyl_compare_double(out.x[1], 0.0, 1e-14) );
  TEST_CHECK( gkyl_compare_double(out.x[2], 5.0, 1e-14) );

  // con_to_cov scales the angular component by r^2 and leaves others.
  struct gkyl_vec3 cov = gkyl_vec3_polar_con_to_cov(2.0, gkyl_vec3_new(1.0, 3.0, 7.0));
  TEST_CHECK( gkyl_compare_double(cov.x[0], 1.0, 1e-14) );
  TEST_CHECK( gkyl_compare_double(cov.x[1], 3.0*4.0, 1e-14) );
  TEST_CHECK( gkyl_compare_double(cov.x[2], 7.0, 1e-14) );
}

void
test_minmod_2()
{
  TEST_CHECK( gkyl_minmod_2(2.0, 3.0) == 2.0 );
  TEST_CHECK( gkyl_minmod_2(3.0, 2.0) == 2.0 );
  TEST_CHECK( gkyl_minmod_2(-2.0, -3.0) == -2.0 );
  TEST_CHECK( gkyl_minmod_2(-3.0, -2.0) == -2.0 );
  // opposite signs -> 0
  TEST_CHECK( gkyl_minmod_2(2.0, -3.0) == 0.0 );
  TEST_CHECK( gkyl_minmod_2(-2.0, 3.0) == 0.0 );
  TEST_CHECK( gkyl_minmod_2(0.0, 5.0) == 0.0 );
}

void
test_minmod_3()
{
  TEST_CHECK( gkyl_minmod_3(2.0, 3.0, 4.0) == 2.0 );
  TEST_CHECK( gkyl_minmod_3(-2.0, -3.0, -4.0) == -2.0 );
  TEST_CHECK( gkyl_minmod_3(2.0, -3.0, 4.0) == 0.0 );
  TEST_CHECK( gkyl_minmod_3(1.0, 2.0, 0.0) == 0.0 );
}

void
test_minmod_4()
{
  TEST_CHECK( gkyl_minmod_4(2.0, 3.0, 4.0, 5.0) == 2.0 );
  TEST_CHECK( gkyl_minmod_4(-2.0, -3.0, -4.0, -5.0) == -2.0 );
  TEST_CHECK( gkyl_minmod_4(2.0, 3.0, -4.0, 5.0) == 0.0 );
}

void
test_median()
{
  TEST_CHECK( gkyl_compare_double(gkyl_median(1.0, 2.0, 3.0), 2.0, 1e-15) );
  TEST_CHECK( gkyl_compare_double(gkyl_median(3.0, 1.0, 2.0), 2.0, 1e-15) );
  TEST_CHECK( gkyl_compare_double(gkyl_median(2.0, 3.0, 1.0), 2.0, 1e-15) );
}

void
test_min_max_3()
{
  TEST_CHECK( gkyl_min_3(1.0, 2.0, 3.0) == 1.0 );
  TEST_CHECK( gkyl_min_3(3.0, -1.0, 2.0) == -1.0 );
  TEST_CHECK( gkyl_max_3(1.0, 2.0, 3.0) == 3.0 );
  TEST_CHECK( gkyl_max_3(3.0, -1.0, 2.0) == 3.0 );
  TEST_CHECK( gkyl_max_3(-5.0, -2.0, -9.0) == -2.0 );
}

TEST_LIST = {
  { "vec3_basic", test_vec3_basic },
  { "vec3_scale", test_vec3_scale },
  { "vec3_add_sub", test_vec3_add_sub },
  { "vec3_len_norm", test_vec3_len_norm },
  { "vec3_dot", test_vec3_dot },
  { "vec3_cross", test_vec3_cross },
  { "vec3_triple", test_vec3_triple },
  { "vec3_polar", test_vec3_polar },
  { "minmod_2", test_minmod_2 },
  { "minmod_3", test_minmod_3 },
  { "minmod_4", test_minmod_4 },
  { "median", test_median },
  { "min_max_3", test_min_max_3 },
  { NULL, NULL },
};
