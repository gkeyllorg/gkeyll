// Additional unit tests for gkyl_mat (dense column-major matrices)
#include <acutest.h>
#include <gkyl_mat.h>
#include <gkyl_util.h>

void
test_mat_new_set_get()
{
  struct gkyl_mat *m = gkyl_mat_new(3, 2, 7.0);
  TEST_CHECK( m->nr == 3 );
  TEST_CHECK( m->nc == 2 );
  // initial value
  for (size_t r=0; r<3; ++r)
    for (size_t c=0; c<2; ++c)
      TEST_CHECK( gkyl_mat_get(m, r, c) == 7.0 );

  gkyl_mat_set(m, 0, 0, 1.0);
  gkyl_mat_set(m, 2, 1, -3.5);
  TEST_CHECK( gkyl_mat_get(m, 0, 0) == 1.0 );
  TEST_CHECK( gkyl_mat_get(m, 2, 1) == -3.5 );
  // unchanged entry
  TEST_CHECK( gkyl_mat_get(m, 1, 0) == 7.0 );

  gkyl_mat_release(m);
}

void
test_mat_clear()
{
  struct gkyl_mat *m = gkyl_mat_new(4, 4, 1.0);
  gkyl_mat_clear(m, 0.0);
  for (size_t r=0; r<4; ++r)
    for (size_t c=0; c<4; ++c)
      TEST_CHECK( gkyl_mat_get(m, r, c) == 0.0 );
  gkyl_mat_release(m);
}

void
test_mat_diag()
{
  struct gkyl_mat *m = gkyl_mat_new(3, 3, 9.0);
  gkyl_mat_diag(m, 2.0);
  for (size_t r=0; r<3; ++r)
    for (size_t c=0; c<3; ++c)
      TEST_CHECK( gkyl_mat_get(m, r, c) == (r==c ? 2.0 : 0.0) );
  gkyl_mat_release(m);
}

void
test_mat_mm_identity()
{
  // A * I == A
  struct gkyl_mat *A = gkyl_mat_new(2, 2, 0.0);
  gkyl_mat_set(A, 0, 0, 1.0); gkyl_mat_set(A, 0, 1, 2.0);
  gkyl_mat_set(A, 1, 0, 3.0); gkyl_mat_set(A, 1, 1, 4.0);

  struct gkyl_mat *I = gkyl_mat_new(2, 2, 0.0);
  gkyl_mat_diag(I, 1.0);

  struct gkyl_mat *C = gkyl_mat_new(2, 2, 0.0);
  gkyl_mat_mm(1.0, 0.0, GKYL_NO_TRANS, A, GKYL_NO_TRANS, I, C, false);

  for (size_t r=0; r<2; ++r)
    for (size_t c=0; c<2; ++c)
      TEST_CHECK( gkyl_compare_double(gkyl_mat_get(C, r, c), gkyl_mat_get(A, r, c), 1e-14) );

  gkyl_mat_release(A); gkyl_mat_release(I); gkyl_mat_release(C);
}

void
test_mat_mm_known()
{
  // [1 2; 3 4] * [5 6; 7 8] = [19 22; 43 50]
  struct gkyl_mat *A = gkyl_mat_new(2, 2, 0.0);
  gkyl_mat_set(A, 0, 0, 1.0); gkyl_mat_set(A, 0, 1, 2.0);
  gkyl_mat_set(A, 1, 0, 3.0); gkyl_mat_set(A, 1, 1, 4.0);

  struct gkyl_mat *B = gkyl_mat_new(2, 2, 0.0);
  gkyl_mat_set(B, 0, 0, 5.0); gkyl_mat_set(B, 0, 1, 6.0);
  gkyl_mat_set(B, 1, 0, 7.0); gkyl_mat_set(B, 1, 1, 8.0);

  struct gkyl_mat *C = gkyl_mat_new(2, 2, 0.0);
  gkyl_mat_mm(1.0, 0.0, GKYL_NO_TRANS, A, GKYL_NO_TRANS, B, C, false);

  TEST_CHECK( gkyl_compare_double(gkyl_mat_get(C, 0, 0), 19.0, 1e-13) );
  TEST_CHECK( gkyl_compare_double(gkyl_mat_get(C, 0, 1), 22.0, 1e-13) );
  TEST_CHECK( gkyl_compare_double(gkyl_mat_get(C, 1, 0), 43.0, 1e-13) );
  TEST_CHECK( gkyl_compare_double(gkyl_mat_get(C, 1, 1), 50.0, 1e-13) );

  // alpha scaling: 2*(A*B)
  struct gkyl_mat *C2 = gkyl_mat_new(2, 2, 0.0);
  gkyl_mat_mm(2.0, 0.0, GKYL_NO_TRANS, A, GKYL_NO_TRANS, B, C2, false);
  TEST_CHECK( gkyl_compare_double(gkyl_mat_get(C2, 0, 0), 38.0, 1e-13) );
  TEST_CHECK( gkyl_compare_double(gkyl_mat_get(C2, 1, 1), 100.0, 1e-13) );

  gkyl_mat_release(A); gkyl_mat_release(B);
  gkyl_mat_release(C); gkyl_mat_release(C2);
}

void
test_mat_mm_transpose()
{
  // A^T with A = [1 2; 3 4] is [1 3; 2 4]; (A^T)*I checks transpose handling
  struct gkyl_mat *A = gkyl_mat_new(2, 2, 0.0);
  gkyl_mat_set(A, 0, 0, 1.0); gkyl_mat_set(A, 0, 1, 2.0);
  gkyl_mat_set(A, 1, 0, 3.0); gkyl_mat_set(A, 1, 1, 4.0);

  struct gkyl_mat *I = gkyl_mat_new(2, 2, 0.0);
  gkyl_mat_diag(I, 1.0);

  struct gkyl_mat *C = gkyl_mat_new(2, 2, 0.0);
  gkyl_mat_mm(1.0, 0.0, GKYL_TRANS, A, GKYL_NO_TRANS, I, C, false);

  TEST_CHECK( gkyl_compare_double(gkyl_mat_get(C, 0, 0), 1.0, 1e-14) );
  TEST_CHECK( gkyl_compare_double(gkyl_mat_get(C, 0, 1), 3.0, 1e-14) );
  TEST_CHECK( gkyl_compare_double(gkyl_mat_get(C, 1, 0), 2.0, 1e-14) );
  TEST_CHECK( gkyl_compare_double(gkyl_mat_get(C, 1, 1), 4.0, 1e-14) );

  gkyl_mat_release(A); gkyl_mat_release(I); gkyl_mat_release(C);
}

void
test_nmat_new()
{
  // Batched matrices: 3 matrices each 2x2.
  struct gkyl_nmat *nm = gkyl_nmat_new(3, 2, 2);
  TEST_CHECK( nm->num == 3 );
  TEST_CHECK( nm->nr == 2 );
  TEST_CHECK( nm->nc == 2 );
  TEST_CHECK( !gkyl_nmat_is_cu_dev(nm) );

  // Fill each matrix differently and read back.
  for (size_t k=0; k<3; ++k) {
    struct gkyl_mat mk = gkyl_nmat_get(nm, k);
    gkyl_mat_clear(&mk, 0.0);
    gkyl_mat_set(&mk, 0, 0, (double) k+1);
  }
  for (size_t k=0; k<3; ++k) {
    struct gkyl_mat mk = gkyl_nmat_get(nm, k);
    TEST_CHECK( gkyl_mat_get(&mk, 0, 0) == (double) k+1 );
  }

  gkyl_nmat_release(nm);
}

TEST_LIST = {
  { "mat_new_set_get", test_mat_new_set_get },
  { "mat_clear", test_mat_clear },
  { "mat_diag", test_mat_diag },
  { "mat_mm_identity", test_mat_mm_identity },
  { "mat_mm_known", test_mat_mm_known },
  { "mat_mm_transpose", test_mat_mm_transpose },
  { "nmat_new", test_nmat_new },
  { NULL, NULL },
};
