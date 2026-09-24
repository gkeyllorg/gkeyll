#include <acutest.h>
#include <gkyl_mat_triples.h>
#include <gkyl_util.h>
#include <gkyl_range.h>

void test_tri_1()
{
  gkyl_mat_triples *tri = gkyl_mat_triples_new(5, 5);

  /*  A : matrix([s,0,u,u,0],[l,u,0,0,0],[0,l,p,0,0],[0,0,0,e,u],[l,l,0,0,r]); */

  // row 0
  gkyl_mat_triples_insert(tri, 0, 0, 19.0);
  gkyl_mat_triples_insert(tri, 0, 2, 21.0);
  gkyl_mat_triples_insert(tri, 0, 3, 21.0);

  // row 1
  gkyl_mat_triples_insert(tri, 1, 0, 21.0);
  gkyl_mat_triples_insert(tri, 1, 1, 21.0);

  // row 2
  gkyl_mat_triples_insert(tri, 2, 1, 12.0);
  gkyl_mat_triples_insert(tri, 2, 2, 16.0);

  // row 3
  gkyl_mat_triples_insert(tri, 3, 3, 5.0);
  gkyl_mat_triples_insert(tri, 3, 4, 21.0);

  // row 4
  gkyl_mat_triples_insert(tri, 4, 0, 12.0);
  gkyl_mat_triples_insert(tri, 4, 1, 12.0);
  gkyl_mat_triples_insert(tri, 4, 4, 18.0);

  TEST_CHECK( gkyl_mat_triples_get(tri, 4, 4) == 18.0 );
  TEST_CHECK( gkyl_mat_triples_get(tri, 4, 0) == 12.0 );

  TEST_CHECK( gkyl_mat_triples_get(tri, 0, 1) == 0.0 ); // zero element

  TEST_CHECK( gkyl_mat_triples_size(tri) == 12 );

  gkyl_mat_triples_release(tri);
}

void test_tri_2()
{
  gkyl_mat_triples *tri = gkyl_mat_triples_new(3, 3);

  // row 1
  gkyl_mat_triples_insert(tri, 1, 1, 2.0);
  gkyl_mat_triples_insert(tri, 1, 0, 2.1);

  // row 0
  gkyl_mat_triples_insert(tri, 0, 0, 1.0);
  gkyl_mat_triples_insert(tri, 0, 2, 1.1);

  // row 2
  gkyl_mat_triples_insert(tri, 2, 2, 3.0);

  TEST_CHECK( gkyl_mat_triples_size(tri) == 5 );

  // order in which col-maj sorting should return indices
  size_t cm_idx[][2] = {
    {0,0}, {1,0}, {1,1}, {0,2}, {2,2}
  };

  double vals[] = { 1.0, 2.1, 2.0, 1.1, 3.0 };

  int i = 0;
  gkyl_mat_triples_iter *iter = gkyl_mat_triples_iter_new(tri);
  while (gkyl_mat_triples_iter_next(iter)) {
    struct gkyl_mtriple mt = gkyl_mat_triples_iter_at(iter);
    
    TEST_CHECK( (mt.row == cm_idx[i][0]) && (mt.col == cm_idx[i][1]) );
    TEST_CHECK( mt.val == vals[i] );
    
    i += 1;
  }

  gkyl_mat_triples_iter_release(iter);
  gkyl_mat_triples_release(tri);
}

void test_tri_3()
{
  gkyl_mat_triples *tri = gkyl_mat_triples_new(3, 3);

  double vals[] = { 1.0, 2.1, 2.0, 1.1, 3.0 };

  // row 1
  gkyl_mat_triples_insert(tri, 1, 1, 2.0);
  gkyl_mat_triples_insert(tri, 1, 0, 2.1);
  // row 0
  gkyl_mat_triples_insert(tri, 0, 0, 1.0);
  gkyl_mat_triples_insert(tri, 0, 2, 1.1);
  // row 2
  gkyl_mat_triples_insert(tri, 2, 2, 3.0);

  int i = 0;
  gkyl_mat_triples_iter *iter = gkyl_mat_triples_iter_new(tri);
  while (gkyl_mat_triples_iter_next(iter)) {
    struct gkyl_mtriple mt = gkyl_mat_triples_iter_at(iter);
    TEST_CHECK( mt.val == vals[i] );
    i += 1;
  }

  // Zero out the values and check.
  gkyl_mat_triples_clear(tri, 0.);
  gkyl_mat_triples_iter_init(iter, tri);
  while (gkyl_mat_triples_iter_next(iter)) {
    struct gkyl_mtriple mt = gkyl_mat_triples_iter_at(iter);
    TEST_CHECK( mt.val == 0. );
    TEST_MSG("Expected: %.13e in (%zu,%zu)", 0., mt.row, mt.col);
    TEST_MSG("Produced: %.13e", mt.val);
  }

  // Test that assigning/accumulating again works.
  gkyl_mat_triples_accum(tri, 1, 1, 2.0);
  gkyl_mat_triples_accum(tri, 1, 0, 2.1);
  gkyl_mat_triples_accum(tri, 0, 0, 1.0);
  gkyl_mat_triples_accum(tri, 0, 2, 1.1);
  gkyl_mat_triples_accum(tri, 2, 2, 3.0);
  i = 0;
  gkyl_mat_triples_iter_init(iter, tri);
  while (gkyl_mat_triples_iter_next(iter)) {
    struct gkyl_mtriple mt = gkyl_mat_triples_iter_at(iter);
    TEST_CHECK( mt.val == vals[i] );
    i += 1;
  }

  gkyl_mat_triples_iter_release(iter);
  gkyl_mat_triples_release(tri);
}

TEST_LIST = {
  { "tri_1", test_tri_1 },
  { "tri_2", test_tri_2 },
  { "tri_3", test_tri_3 },
  { NULL, NULL }
};
