#include <math.h>
#include <acutest.h>
#include <../zero/wave_geom.c>

static void
test_tri()
{
  struct gkyl_vec3 p1 = gkyl_vec3_new(1, 1, 2);
  struct gkyl_vec3 p2 = gkyl_vec3_new(-2, 4, 3);
  struct gkyl_vec3 p3 = gkyl_vec3_new(-2, -2, 4);

  double area = triangle_area(p1, p2, p3);

  TEST_CHECK( gkyl_compare_double( area, 10.173494974687902, 1e-15) );
}

static void
test_planar_quad_1()
{
  struct gkyl_vec3 p1 = gkyl_vec3_new(0, 0, 0);
  struct gkyl_vec3 p2 = gkyl_vec3_new(1, 2, 3);
  struct gkyl_vec3 p3 = gkyl_vec3_new(5, 7, 9);
  struct gkyl_vec3 p4 = gkyl_vec3_new(4, 5, 6);

  struct gkyl_vec3 norm;

  double area = planar_quad_area_norm(
      p1, p2, p3, p4, &norm);

  TEST_CHECK( gkyl_compare_double( area, 3*sqrt(6), 1e-15) );

  TEST_CHECK( gkyl_compare_double( gkyl_vec3_len(norm), 1, 1e-15) );

  struct gkyl_vec3 v12 = gkyl_vec3_sub(p2, p1);
  struct gkyl_vec3 v23 = gkyl_vec3_sub(p3, p2);
  struct gkyl_vec3 v34 = gkyl_vec3_sub(p4, p3);
  struct gkyl_vec3 v41 = gkyl_vec3_sub(p1, p2);
  TEST_CHECK( gkyl_compare_double( gkyl_vec3_dot(norm, v12), 0, 1e-15) );
  TEST_CHECK( gkyl_compare_double( gkyl_vec3_dot(norm, v23), 0, 1e-15) );
  TEST_CHECK( gkyl_compare_double( gkyl_vec3_dot(norm, v34), 0, 1e-15) );
  TEST_CHECK( gkyl_compare_double( gkyl_vec3_dot(norm, v41), 0, 1e-15) );
}

static void
test_quad_1()
{
  struct gkyl_vec3 p1 = gkyl_vec3_new(0, 0, 0);
  struct gkyl_vec3 p2 = gkyl_vec3_new(1, 2, 3);
  struct gkyl_vec3 p3 = gkyl_vec3_new(5, 7, 9);
  struct gkyl_vec3 p4 = gkyl_vec3_new(4, 5, 6);

  struct gkyl_vec3 norm, tau1, tau2;

  double area = quad_area_norm_tang(p1, p2, p3, p4, &norm, &tau1, &tau2);

  TEST_CHECK( gkyl_compare_double( area, 3*sqrt(6), 1e-15) );

  TEST_CHECK( gkyl_compare_double( gkyl_vec3_len(norm), 1, 1e-15) );
  TEST_CHECK( gkyl_compare_double( gkyl_vec3_len(tau1), 1, 1e-15) );
  TEST_CHECK( gkyl_compare_double( gkyl_vec3_len(tau2), 1, 1e-15) );

  TEST_CHECK( gkyl_compare_double( gkyl_vec3_dot(norm, tau1), 0, 1e-15) );
  TEST_CHECK( gkyl_compare_double( gkyl_vec3_dot(norm, tau2), 0, 1e-15) );
  TEST_CHECK( gkyl_compare_double( gkyl_vec3_dot(tau1, tau2), 0, 1e-15) );

  struct gkyl_vec3 cross = gkyl_vec3_cross(tau1, tau2);
  TEST_CHECK( gkyl_compare_double( cross.x[0], norm.x[0], 1e-15) );
  TEST_CHECK( gkyl_compare_double( cross.x[1], norm.x[1], 1e-15) );
  TEST_CHECK( gkyl_compare_double( cross.x[2], norm.x[2], 1e-15) );

  struct gkyl_vec3 v12 = gkyl_vec3_sub(p2, p1);
  struct gkyl_vec3 v23 = gkyl_vec3_sub(p3, p2);
  struct gkyl_vec3 v34 = gkyl_vec3_sub(p4, p3);
  struct gkyl_vec3 v41 = gkyl_vec3_sub(p1, p2);
  TEST_CHECK( gkyl_compare_double( gkyl_vec3_dot(norm, v12), 0, 1e-15) );
  TEST_CHECK( gkyl_compare_double( gkyl_vec3_dot(norm, v23), 0, 1e-15) );
  TEST_CHECK( gkyl_compare_double( gkyl_vec3_dot(norm, v34), 0, 1e-15) );
  TEST_CHECK( gkyl_compare_double( gkyl_vec3_dot(norm, v41), 0, 1e-15) );
}

static void
test_quad_2()
{
  struct gkyl_vec3 p1 = gkyl_vec3_new(9., 3., 4.);
  struct gkyl_vec3 p2 = gkyl_vec3_new(9., 8., 7.);
  struct gkyl_vec3 p3 = gkyl_vec3_new(5., 2., 4.);
  struct gkyl_vec3 p4 = gkyl_vec3_new(7., 9., 9.);

  struct gkyl_vec3 norm, tau1, tau2;

  double area = quad_area_norm_tang(p1, p2, p3, p4, &norm, &tau1, &tau2);

  TEST_CHECK( gkyl_compare_double( area, 5.099019513592786, 1e-15) );

  TEST_CHECK( gkyl_compare_double( gkyl_vec3_len(norm), 1, 1e-15) );
  TEST_CHECK( gkyl_compare_double( gkyl_vec3_len(tau1), 1, 1e-15) );
  TEST_CHECK( gkyl_compare_double( gkyl_vec3_len(tau2), 1, 1e-15) );

  TEST_CHECK( gkyl_compare_double( norm.x[0], -0.196116, 1e-5) );
  TEST_CHECK( gkyl_compare_double( norm.x[1],  0.784465, 1e-5) );
  TEST_CHECK( gkyl_compare_double( norm.x[2], -0.588348, 1e-5) );

  TEST_CHECK( gkyl_compare_double( gkyl_vec3_dot(norm, tau1), 0, 1e-15) );
  TEST_CHECK( gkyl_compare_double( gkyl_vec3_dot(norm, tau2), 0, 1e-15) );
  TEST_CHECK( gkyl_compare_double( gkyl_vec3_dot(tau1, tau2), 0, 1e-15) );

  struct gkyl_vec3 cross = gkyl_vec3_cross(tau1, tau2);
  TEST_CHECK( gkyl_compare_double( cross.x[0], norm.x[0], 1e-15) );
  TEST_CHECK( gkyl_compare_double( cross.x[1], norm.x[1], 1e-15) );
  TEST_CHECK( gkyl_compare_double( cross.x[2], norm.x[2], 1e-15) );
}

static void
test_vol_tetra_1()
{
  struct gkyl_vec3 p1 = gkyl_vec3_new(0, 0, 0);
  struct gkyl_vec3 p2 = gkyl_vec3_new(1, 0, 0);
  struct gkyl_vec3 p3 = gkyl_vec3_new(1, 1, 0);
  struct gkyl_vec3 p4 = gkyl_vec3_new(1, 1, 1);
  double vol = vol_tetra(p1, p2, p3, p4);
  TEST_CHECK( gkyl_compare_double( vol, 1./6., 1e-15) );
}

static void
test_vol_tetra_2()
{
  struct gkyl_vec3 p1 = gkyl_vec3_new(0.37, 0.07, 0.29);
  struct gkyl_vec3 p2 = gkyl_vec3_new(1.03, 0.08, 0.05);
  struct gkyl_vec3 p3 = gkyl_vec3_new(1.1 , 1.2 , 0.32);
  struct gkyl_vec3 p4 = gkyl_vec3_new(1.17, 1.18, 1.31);
  double vol = vol_tetra(p1, p2, p3, p4);
  TEST_CHECK( gkyl_compare_double( vol, 0.12567, 1e-15) );
}

static void
test_vol_hexa_1()
{
  struct gkyl_vec3 verts[8] = {
    {0, 0, 0},
    {1, 0, 0},
    {1, 1, 0},
    {0, 1, 0},
    {0, 0, 1},
    {1, 0, 1},
    {1, 1, 1},
    {0, 1, 1},
  };
  
  double vol = vol_hexa(verts);
  TEST_CHECK( gkyl_compare_double( vol, 1, 1e-15) );
}

static void
test_vol_hexa_2()
{
  struct gkyl_vec3 verts[8] = {
    {0.37, 0.07, 0.21},
    {1.16, 0.2 , 0.29},
    {1.33, 1.18, 0.23},
    {0.01, 1.21, 0.28},
    {0.38, 0.36, 1.02},
    {1.19, 0.11, 1.01},
    {1.04, 1.23, 1.14},
    {0.22, 1.2 , 1.18},
  };
  
  double vol = vol_hexa(verts);
  TEST_CHECK( gkyl_compare_double( vol, 0.746420666666667, 1e-15) );
}

TEST_LIST = {
  { "triangle", test_tri },
  { "parallelogram_as_planar_quad", test_planar_quad_1 },
  { "parallelogram_as_quad", test_quad_1 },
  { "quad", test_quad_2 },
  { "vol_tetra_1", test_vol_tetra_1 },
  { "vol_tetra_2", test_vol_tetra_2 },
  { "vol_hexa_1", test_vol_hexa_1 },
  { "vol_hexa_2", test_vol_hexa_2 },
  { NULL, NULL },
};
