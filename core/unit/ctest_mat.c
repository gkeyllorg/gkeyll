#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_mat.h>
#include <gkyl_mat_priv.h>

void
test_mat_base()
{
  struct gkyl_mat *m = gkyl_mat_new(10, 20, 0.25);

  TEST_CHECK( 10 == m->nr );
  TEST_CHECK( 20 == m->nc );

  for (size_t j=0; j<m->nc; ++j)
    for (size_t i=0; i<m->nr; ++i)
      TEST_CHECK ( 0.25 == gkyl_mat_get(m, i, j) );

  gkyl_mat_clear(m, 0.1);

  for (size_t j=0; j<m->nc; ++j)
    for (size_t i=0; i<m->nr; ++i)
      TEST_CHECK ( 0.1 == gkyl_mat_get(m, i, j) );

  size_t count = 0;
  for (size_t j=0; j<m->nc; ++j)
    for (size_t i=0; i<m->nr; ++i)
      gkyl_mat_set(m, i, j, count++);

  count = 0;
  for (size_t j=0; j<m->nc; ++j) {
    const double* col = gkyl_mat_get_ccol(m, j);
    for (size_t i=0; i<m->nr; ++i)
      TEST_CHECK( col[i] == count++ );
  }

  count = 0;
  for (size_t i=0; i<m->nr*m->nc; ++i)
    TEST_CHECK( m->data[i] == i );

  gkyl_mat_diag(m, 1.0);
  
  for (size_t j=0; j<m->nc; ++j)
    for (size_t i=0; i<m->nr; ++i)
      TEST_CHECK( gkyl_mat_get(m, i, j) == ( i==j ? 1.0 : 0.0 ) );

  struct gkyl_mat *m2 = gkyl_mat_clone(m);

  TEST_CHECK( m2->nr == m->nr );
  TEST_CHECK( m2->nc == m->nc );

  for (size_t j=0; j<m->nc; ++j)
    for (size_t i=0; i<m->nr; ++i)
      TEST_CHECK( gkyl_mat_get(m, i, j) == gkyl_mat_get(m2, i, j) );

  double old = gkyl_mat_get(m, 3, 4);
  double inc = -123.0;
  gkyl_mat_inc(m, 3, 4, inc);
  TEST_CHECK( gkyl_mat_get(m, 3, 4) == old + inc );

  gkyl_mat_release(m);
  gkyl_mat_release(m2);
}

void
test_mat_mm_op()
{
  struct gkyl_mat *A = gkyl_mat_new(2, 3, 0.0);
  struct gkyl_mat *B = gkyl_mat_new(3, 2, 0.0);
  struct gkyl_mat *C = gkyl_mat_new(2, 2, 0.0);
  struct gkyl_mat *D = gkyl_mat_new(3, 3, 0.0);

  double val = 1.0;

  // A : matrix( [1,2,3], [4,5,6] );
  for (int i=0; i<A->nr; ++i)  
    for (int j=0; j<A->nc; ++j) {
      gkyl_mat_set(A,i,j,val);
      val += 1.0;
    }

  // B : matrix( [7,8], [9,10], [11,12] );
  for (int i=0; i<B->nr; ++i)  
    for (int j=0; j<B->nc; ++j) {
      gkyl_mat_set(B,i,j,val);
      val += 1.0;
    }

  // C = 0.5*A*B + 0.0*C
  gkyl_mat_mm(0.5, 0.0, GKYL_NO_TRANS, A, GKYL_NO_TRANS, B, C, false);

  // C : matrix( [29.0, 32.0], [69.5, 77.0] )
  TEST_CHECK( gkyl_mat_get(C, 0, 0) == 29.0 );
  TEST_CHECK( gkyl_mat_get(C, 0, 1) == 32.0 );
  TEST_CHECK( gkyl_mat_get(C, 1, 0) == 69.5 );
  TEST_CHECK( gkyl_mat_get(C, 1, 1) == 77.0 );

  // D = 0.5*A'*B'
  gkyl_mat_mm(0.5, 0.0, GKYL_TRANS, A, GKYL_TRANS, B, D, false);

  // D : matrix( [ 19.5  24.5  29.5 ], [ 27.0  34.0  41.0 ], [ 34.5  43.5  52.5 ] )
  TEST_CHECK( gkyl_mat_get(D, 0, 0) == 19.5 );
  TEST_CHECK( gkyl_mat_get(D, 0, 1) == 24.5 );
  TEST_CHECK( gkyl_mat_get(D, 0, 2) == 29.5 );

  TEST_CHECK( gkyl_mat_get(D, 1, 0) == 27.0 );
  TEST_CHECK( gkyl_mat_get(D, 1, 1) == 34.0 );
  TEST_CHECK( gkyl_mat_get(D, 1, 2) == 41.0 );

  TEST_CHECK( gkyl_mat_get(D, 2, 0) == 34.5 );
  TEST_CHECK( gkyl_mat_get(D, 2, 1) == 43.5 );
  TEST_CHECK( gkyl_mat_get(D, 2, 2) == 52.5 );
  
  gkyl_mat_release(A);
  gkyl_mat_release(B);
  gkyl_mat_release(C);
  gkyl_mat_release(D);
}

void
test_mat_linsolve()
{
  struct gkyl_mat *A = gkyl_mat_new(3, 3, 0.0);
  struct gkyl_mat *x = gkyl_mat_new(3, 1, 0.0);

  double val = 1.0;

  // A : matrix( [1,2,3], [4,5,6], [7,8,10] );
  for (int i=0; i<A->nr; ++i)  
    for (int j=0; j<A->nc; ++j) {
      gkyl_mat_set(A,i,j,val);
      val += 1.0;
    }
  gkyl_mat_set(A,2,2,10.0); // ensures determinant is not zero
  struct gkyl_mat *AA = gkyl_mat_clone(A);

  //gkyl_mat_show("A", stdout, A);

  // x = matrix( [1, 1, 1] );
  gkyl_mat_clear(x, 1.0);

  // memory for pivot vector
  gkyl_mem_buff ipiv = gkyl_mem_buff_new(sizeof(long[3]));

  // solve linear system: sol : matrix( [-1, 1, 0] )
  bool status = gkyl_mat_linsolve_lu(A, x, gkyl_mem_buff_data(ipiv));

  TEST_CHECK( status );

  //gkyl_mat_show("A", stdout, A);
  //gkyl_mat_show("x", stdout, x);

  TEST_CHECK( gkyl_compare(gkyl_mat_get(x,0,0), -1.0, 1e-15) );
  TEST_CHECK( gkyl_compare(gkyl_mat_get(x,1,0), 1.0, 1e-15) );
  TEST_CHECK( gkyl_compare(gkyl_mat_get(x,2,0), 0.0, 1e-15) );
 
  // trivial extension of the test above; rhs is two column vectors
  struct gkyl_mat *xx = gkyl_mat_new(3, 2, 1.0);
  gkyl_mat_set(xx, 0, 1, 2.0);
  gkyl_mat_set(xx, 1, 1, 2.0);
  gkyl_mat_set(xx, 2, 1, 2.0);
  status = gkyl_mat_linsolve_lu(AA, xx, gkyl_mem_buff_data(ipiv));

  TEST_CHECK( gkyl_compare(gkyl_mat_get(xx,0,0), -1.0, 1e-15) );
  TEST_CHECK( gkyl_compare(gkyl_mat_get(xx,1,0), 1.0, 1e-15) );
  TEST_CHECK( gkyl_compare(gkyl_mat_get(xx,2,0), 0.0, 1e-15) );
  TEST_CHECK( gkyl_compare(gkyl_mat_get(xx,0,1), -2.0, 1e-15) );
  TEST_CHECK( gkyl_compare(gkyl_mat_get(xx,1,1), 2.0, 1e-15) );
  TEST_CHECK( gkyl_compare(gkyl_mat_get(xx,2,1), 0.0, 1e-15) );

  gkyl_mat_release(AA);
  gkyl_mat_release(xx);
  gkyl_mat_release(A);
  gkyl_mat_release(x);
  gkyl_mem_buff_release(ipiv);
}

void
test_nmat_base()
{
  // 5 matrices with shape 10x20
  struct gkyl_nmat *nmat = gkyl_nmat_new(5, 10, 20);

  TEST_CHECK( false == gkyl_nmat_is_cu_dev(nmat) );

  TEST_CHECK( 5 == nmat->num );
  TEST_CHECK( 10 == nmat->nr );
  TEST_CHECK( 20 == nmat->nc );

  struct gkyl_mat m = gkyl_nmat_get(nmat, 0);
  TEST_CHECK( 10 == m.nr );
  TEST_CHECK( 20 == m.nc );

  for (size_t n=0; n<nmat->num; ++n) {
    struct gkyl_mat m = gkyl_nmat_get(nmat, n);
    
    for (size_t j=0; j<nmat->nc; ++j)
      for (size_t i=0; i<nmat->nr; ++i)
        gkyl_mat_set(&m, i, j, n*0.5);
  }

  for (size_t n=0; n<nmat->num; ++n) {
    struct gkyl_mat m = gkyl_nmat_get(nmat, n);
    
    for (size_t j=0; j<nmat->nc; ++j)
      for (size_t i=0; i<nmat->nr; ++i)
        TEST_CHECK ( n*0.5 == gkyl_mat_get(&m, i, j) );
  }  

  for (size_t n=0; n<nmat->num; ++n)
    for (size_t i=0; i<nmat->nr*nmat->nc; ++i)
      TEST_CHECK( nmat->mptr[n][i] == n*0.5 );

  // copy matrix
  struct gkyl_nmat *ncpy = gkyl_nmat_new(nmat->num, nmat->nr, nmat->nc);
  gkyl_nmat_copy(ncpy, nmat);

  for (size_t n=0; n<ncpy->num; ++n) {
    struct gkyl_mat m = gkyl_nmat_get(ncpy, n);
    
    for (size_t j=0; j<ncpy->nc; ++j)
      for (size_t i=0; i<ncpy->nr; ++i)
        TEST_CHECK ( n*0.5 == gkyl_mat_get(&m, i, j) );
  }  

  for (size_t n=0; n<ncpy->num; ++n)
    for (size_t i=0; i<ncpy->nr*ncpy->nc; ++i)
      TEST_CHECK( ncpy->mptr[n][i] == n*0.5 );  

  gkyl_nmat_release(nmat);
  gkyl_nmat_release(ncpy);
}

void
test_nmat_linsolve_(bool pre_alloc)
{
  struct gkyl_nmat *As = gkyl_nmat_new(5, 3, 3);
  struct gkyl_nmat *xs = gkyl_nmat_new(5, 3, 1);

  for (int n=0; n<As->num; ++n) {
    struct gkyl_mat A = gkyl_nmat_get(As, n);
    
    double val = 1.0;
    // A : matrix( [1,2,3], [4,5,6], [7,8,10] );
    for (int i=0; i<A.nr; ++i)  
      for (int j=0; j<A.nc; ++j) {
        gkyl_mat_set(&A,i,j,val);
        val += 1.0;
      }
    gkyl_mat_set(&A,2,2,10.0); // ensures determinant is not zero
  }

  // x = matrix( [1, 1, 1] );
  for (size_t n=0; n<As->num; ++n) {
    struct gkyl_mat x = gkyl_nmat_get(xs, n);
    gkyl_mat_clear(&x, 1.0);
  }

  // solve all linear systems: sol : matrix( [-1, 1, 0] )
  bool status = false;

  if (pre_alloc) {
    gkyl_nmat_mem *mem = gkyl_nmat_linsolve_lu_new(As->num, As->nr);
    status = gkyl_nmat_linsolve_lu_pa(mem, As, xs);
    gkyl_nmat_linsolve_lu_release(mem);
  }
  else {
    status = gkyl_nmat_linsolve_lu(As, xs);
  }

  TEST_CHECK( status );

  //gkyl_mat_show("A", stdout, A);
  //gkyl_mat_show("x", stdout, x);

  for (int n=0; n<xs->num; ++n) {
    struct gkyl_mat x = gkyl_nmat_get(xs, n);

    TEST_CHECK( gkyl_compare(gkyl_mat_get(&x,0,0), -1.0, 1e-15) );
    TEST_CHECK( gkyl_compare(gkyl_mat_get(&x,1,0), 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare(gkyl_mat_get(&x,2,0), 0.0, 1e-15) );
  }
  
  gkyl_nmat_release(As);
  gkyl_nmat_release(xs);
}

void test_nmat_linsolve() { test_nmat_linsolve_(false); }
void test_nmat_linsolve_pa() { test_nmat_linsolve_(true); }

#ifdef GKYL_HAVE_CUDA

void
test_cu_nmat_base()
{
  // 5 matrices with shape 10x20
  struct gkyl_nmat *nmat = gkyl_nmat_cu_dev_new(5, 10, 20);

  TEST_CHECK( 5 == nmat->num );
  TEST_CHECK( 10 == nmat->nr );
  TEST_CHECK( 20 == nmat->nc );

  TEST_CHECK( gkyl_nmat_is_cu_dev(nmat) == true );

  // create host-side matrix
  struct gkyl_nmat *h1 = gkyl_nmat_new(5, 10, 20);
  for (size_t n=0; n<h1->num; ++n) {
    struct gkyl_mat m = gkyl_nmat_get(h1, n);
    
    for (size_t j=0; j<h1->nc; ++j)
      for (size_t i=0; i<h1->nr; ++i)
        gkyl_mat_set(&m, i, j, n*0.5);
  }

  // copy to device
  gkyl_nmat_copy(nmat, h1);

  // copy it back to host
  struct gkyl_nmat *h2 = gkyl_nmat_new(5, 10, 20);
  gkyl_nmat_copy(h2, nmat);

  // check
  for (size_t n=0; n<h1->num; ++n) {
    struct gkyl_mat m1 = gkyl_nmat_get(h1, n);
    struct gkyl_mat m2 = gkyl_nmat_get(h2, n);
    
    for (size_t j=0; j<h1->nc; ++j)
      for (size_t i=0; i<h1->nr; ++i)
        TEST_CHECK( gkyl_mat_get(&m1, i, j) == gkyl_mat_get(&m2, i, j) );
  }

  gkyl_nmat_release(nmat);
  gkyl_nmat_release(h1);
  gkyl_nmat_release(h2);
}

void
test_cu_nmat_linsolve_(bool pre_alloc)
{
  struct gkyl_nmat *As = gkyl_nmat_new(5, 3, 3);
  struct gkyl_nmat *xs = gkyl_nmat_new(5, 3, 1);

  for (int n=0; n<As->num; ++n) {
    struct gkyl_mat A = gkyl_nmat_get(As, n);
    
    double val = 1.0;
    // A : matrix( [1,2,3], [4,5,6], [7,8,10] );
    for (int i=0; i<A.nr; ++i)  
      for (int j=0; j<A.nc; ++j) {
        gkyl_mat_set(&A,i,j,val);
        val += 1.0;
      }
    gkyl_mat_set(&A,2,2,10.0); // ensures determinant is not zero
  }

  // x = matrix( [1, 1, 1] );
  for (size_t n=0; n<As->num; ++n) {
    struct gkyl_mat x = gkyl_nmat_get(xs, n);
    gkyl_mat_clear(&x, 1.0);
  }

  // copy data to GPU matrices
  struct gkyl_nmat *As_d = gkyl_nmat_cu_dev_new(5, 3, 3);
  struct gkyl_nmat *xs_d = gkyl_nmat_cu_dev_new(5, 3, 1);

  gkyl_nmat_copy(As_d, As);
  gkyl_nmat_copy(xs_d, xs);

  // solve all linear systems: sol : matrix( [-1, 1, 0] )
  bool status = false;

  if (pre_alloc) {
    gkyl_nmat_mem *mem = gkyl_nmat_linsolve_lu_cu_dev_new(As->num, As->nr);
    status = gkyl_nmat_linsolve_lu_pa(mem, As_d, xs_d);
    gkyl_nmat_linsolve_lu_release(mem);
  }
  else {
    status = gkyl_nmat_linsolve_lu(As_d, xs_d);
  }

  TEST_CHECK( status );

  // copy solution back
  gkyl_nmat_copy(As, As_d);
  gkyl_nmat_copy(xs, xs_d);

  for (int n=0; n<xs->num; ++n) {
    struct gkyl_mat x = gkyl_nmat_get(xs, n);

    TEST_CHECK( gkyl_compare(gkyl_mat_get(&x,0,0), -1.0, 1e-15) );
    TEST_CHECK( gkyl_compare(gkyl_mat_get(&x,1,0), 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare(gkyl_mat_get(&x,2,0), 0.0, 1e-15) );
  }
  
  gkyl_nmat_release(As);
  gkyl_nmat_release(xs);

  gkyl_nmat_release(As_d);
  gkyl_nmat_release(xs_d);
}

void test_cu_nmat_linsolve() { test_cu_nmat_linsolve_(false); }
void test_cu_nmat_linsolve_pa() { test_cu_nmat_linsolve_(true); }

#endif

void test_mat_mv()
{

  struct gkyl_mat *A = gkyl_mat_new(4, 4, 1);
  struct gkyl_mat *x = gkyl_mat_new(4, 1, 2);
  struct gkyl_mat *y = gkyl_mat_new(4, 1, 1);
  // y = 1.0*A*x + 0.0*C
  gkyl_mat_mv(1.0, 0.0, GKYL_NO_TRANS, A, x, y);
  for (size_t j=0; j<y->nc; ++j)
    for (size_t i=0; i<y->nr; ++i)
      TEST_CHECK ( 8 == gkyl_mat_get(y, i, j) );

  gkyl_mat_release(A);
  gkyl_mat_release(x);
  gkyl_mat_release(y);
}

void test_nmat_mv()
{

  // n_do matrices with shape 4x4
  int n_do = 3;
  struct gkyl_nmat *nmat_A = gkyl_nmat_new(n_do, 4, 4);
  struct gkyl_nmat *nmat_x = gkyl_nmat_new(n_do, 4, 1);
  struct gkyl_nmat *nmat_y = gkyl_nmat_new(n_do, 4, 1);

  // fill each A matrix with its number
  for (size_t n=0; n<nmat_A->num; ++n) {
    struct gkyl_mat m = gkyl_nmat_get(nmat_A, n);
    for (size_t j=0; j<nmat_A->nc; ++j)
      for (size_t i=0; i<nmat_A->nr; ++i)
        gkyl_mat_set(&m, i, j, n);
  }

  // fill each x vecotr with ones
  for (size_t n=0; n<nmat_x->num; ++n) {
    struct gkyl_mat m = gkyl_nmat_get(nmat_x, n);
    for (size_t j=0; j<nmat_x->nc; ++j)
      for (size_t i=0; i<nmat_x->nr; ++i)
        gkyl_mat_set(&m, i, j, 1);
  }

  // fill each y matrix with 0
  for (size_t n=0; n<nmat_y->num; ++n) {
    struct gkyl_mat m = gkyl_nmat_get(nmat_y, n);
    for (size_t j=0; j<nmat_y->nc; ++j)
      for (size_t i=0; i<nmat_y->nr; ++i)
        gkyl_mat_set(&m, i, j, 0);
  }

  enum gkyl_mat_trans transa = GKYL_NO_TRANS;
  double alpha = 1.0;
  double beta = 0.0;

  gkyl_nmat_mv(alpha, beta, transa, nmat_A, nmat_x, nmat_y);


  // check the expected result
  for (size_t n=0; n<nmat_y->num; ++n) {
    struct gkyl_mat y = gkyl_nmat_get(nmat_y,n);
    double expected = 4*n;
    for (size_t j=0; j<y.nc; ++j)
      for (size_t i=0; i<y.nr; ++i)
        TEST_CHECK ( expected == gkyl_mat_get(&y, i, j) );
  }

  gkyl_nmat_release(nmat_A);
  gkyl_nmat_release(nmat_x);
  gkyl_nmat_release(nmat_y);
}

void test_nmat_mm()
{

  int n_do = 3;
  struct gkyl_nmat *nmat_A = gkyl_nmat_new(n_do, 4, 3);
  struct gkyl_nmat *nmat_x = gkyl_nmat_new(n_do, 3, 2);
  struct gkyl_nmat *nmat_y = gkyl_nmat_new(n_do, 4, 2);

  // fill each A matrix with its number
  for (size_t n=0; n<nmat_A->num; ++n) {
    struct gkyl_mat m = gkyl_nmat_get(nmat_A, n);
    for (size_t j=0; j<nmat_A->nc; ++j)
      for (size_t i=0; i<nmat_A->nr; ++i)
        gkyl_mat_set(&m, i, j, n);
  }

  // fill each x vector with ones
  for (size_t n=0; n<nmat_x->num; ++n) {
    struct gkyl_mat m = gkyl_nmat_get(nmat_x, n);
    for (size_t j=0; j<nmat_x->nc; ++j)
      for (size_t i=0; i<nmat_x->nr; ++i)
        gkyl_mat_set(&m, i, j, 1);
  }

  // fill each y matrix with 0
  for (size_t n=0; n<nmat_y->num; ++n) {
    struct gkyl_mat m = gkyl_nmat_get(nmat_y, n);
    for (size_t j=0; j<nmat_y->nc; ++j)
      for (size_t i=0; i<nmat_y->nr; ++i)
        gkyl_mat_set(&m, i, j, 0);
  }

  enum gkyl_mat_trans transa = GKYL_NO_TRANS;
  enum gkyl_mat_trans transb = GKYL_NO_TRANS;
  double alpha = 1.0;
  double beta = 0.0;

  gkyl_nmat_mm(alpha, beta, transa, nmat_A, transb, nmat_x, nmat_y);


  // check the expected result
  for (size_t n=0; n<nmat_y->num; ++n) {
    struct gkyl_mat y = gkyl_nmat_get(nmat_y,n);
    double expected = 3*n;
    for (size_t j=0; j<y.nc; ++j)
      for (size_t i=0; i<y.nr; ++i)
        TEST_CHECK ( expected == gkyl_mat_get(&y, i, j) );
  }

  gkyl_nmat_release(nmat_A);
  gkyl_nmat_release(nmat_x);
  gkyl_nmat_release(nmat_y);  
}


void test_mat_mm_arrays()
{
  struct gkyl_mat_mm_array_mem *ctest_prob_mem; 
  ctest_prob_mem = gkyl_mat_mm_array_mem_new(4, 3, 1.0, 0.0, 
    GKYL_NO_TRANS, GKYL_NO_TRANS, false);

  struct gkyl_mat *mat_A = ctest_prob_mem->A;
  struct gkyl_array *array_x = gkyl_array_new(GKYL_DOUBLE, 3, 2);
  struct gkyl_array *array_y = gkyl_array_new(GKYL_DOUBLE, 4, 2);

  // fill each A matrix with its number
  for (size_t j=0; j<mat_A->nc; ++j){
    for (size_t i=0; i<mat_A->nr; ++i){
      double a_val = i*3 + j;
      gkyl_mat_set(mat_A, i, j, a_val);
      double actual = gkyl_mat_get(mat_A, i, j);
      //printf("A(m=%d,n=%d): %1.2e,\n", i, j, actual);
    }
  }

  // fill each x vector with ones
  for (size_t j=0; j<array_x->size; ++j){
    double *x = gkyl_array_fetch(array_x,j);
    for (size_t i=0; i<array_x->ncomp; ++i){
      x[i] = i*2 + j;
      //printf("B(m=%d,n=%d): %1.2e,\n", i, j, x[i]);
    }
  } 

  for (size_t j=0; j<array_y->size; ++j){
    double *y = gkyl_array_fetch(array_y,j);
    for (size_t i=0; i<array_y->ncomp; ++i){
      y[i] = 0.0;
      //printf("C(m=%d,n=%d): %1.2e,\n", i, j, 0.0);
    }
  } 

  // Preform the matrix multiply
  gkyl_mat_mm_array(ctest_prob_mem, array_x, array_y);

  // check the expected result
  double expected_array[8] = {10.0, 13.0, 28.0, 40.0, 46.0, 67.0, 64.0, 94.0};
  for (size_t j=0; j<array_y->size; ++j){
    double *y = gkyl_array_fetch(array_y,j);
    for (size_t i=0; i<array_y->ncomp; ++i){
      double actual = y[i];
      double expected = expected_array[i*2 + j];
      //printf("expected: %1.2e, actual: %1.2e\n", expected, actual);
      TEST_CHECK ( expected == actual );
    }
  }

  gkyl_array_release(array_x);
  gkyl_array_release(array_y);
  gkyl_mat_mm_array_mem_release(ctest_prob_mem);

}


#ifdef GKYL_HAVE_CUDA
void test_cu_nmat_mv()
{

  // n_do matrices with shape 4x4
  int n_do = 3;
  struct gkyl_nmat *nmat_A = gkyl_nmat_new(n_do, 4, 4);
  struct gkyl_nmat *nmat_x = gkyl_nmat_new(n_do, 4, 1);
  struct gkyl_nmat *nmat_y = gkyl_nmat_new(n_do, 4, 1);

  struct gkyl_nmat *nmat_Acu = gkyl_nmat_cu_dev_new(n_do, 4, 4);
  struct gkyl_nmat *nmat_xcu = gkyl_nmat_cu_dev_new(n_do, 4, 1);
  struct gkyl_nmat *nmat_ycu = gkyl_nmat_cu_dev_new(n_do, 4, 1);

  // fill each A matrix with its number
  for (size_t n=0; n<nmat_A->num; ++n) {
    struct gkyl_mat m = gkyl_nmat_get(nmat_A, n);
    for (size_t j=0; j<nmat_A->nc; ++j)
      for (size_t i=0; i<nmat_A->nr; ++i)
        gkyl_mat_set(&m, i, j, n);
  }

  // fill each x vecotr with ones
  for (size_t n=0; n<nmat_x->num; ++n) {
    struct gkyl_mat m = gkyl_nmat_get(nmat_x, n);
    for (size_t j=0; j<nmat_x->nc; ++j)
      for (size_t i=0; i<nmat_x->nr; ++i)
        gkyl_mat_set(&m, i, j, 1);
  }

  // fill each y matrix with 0
  for (size_t n=0; n<nmat_y->num; ++n) {
    struct gkyl_mat m = gkyl_nmat_get(nmat_y, n);
    for (size_t j=0; j<nmat_y->nc; ++j)
      for (size_t i=0; i<nmat_y->nr; ++i)
        gkyl_mat_set(&m, i, j, 0);
  }

  enum gkyl_mat_trans transa = GKYL_NO_TRANS;
  enum gkyl_mat_trans transb = GKYL_NO_TRANS;
  double alpha = 1.0;
  double beta = 0.0;


  // copy to device
  gkyl_nmat_copy(nmat_Acu, nmat_A);
  gkyl_nmat_copy(nmat_xcu, nmat_x);
  gkyl_nmat_copy(nmat_ycu, nmat_y);

  gkyl_nmat_mm(alpha, beta, transa, nmat_Acu, transb, nmat_xcu, nmat_ycu);

  // copy to host
  gkyl_nmat_copy(nmat_A, nmat_Acu);
  gkyl_nmat_copy(nmat_x, nmat_xcu);
  gkyl_nmat_copy(nmat_y, nmat_ycu);


  // check the expected result
  for (size_t n=0; n<nmat_y->num; ++n) {
    struct gkyl_mat y = gkyl_nmat_get(nmat_y,n);
    double expected = 4*n;
    for (size_t j=0; j<y.nc; ++j){
      for (size_t i=0; i<y.nr; ++i){
        double actual = gkyl_mat_get(&y, i, j);
        TEST_CHECK ( expected == actual );
      }
    }
  }
  gkyl_nmat_release(nmat_A);
  gkyl_nmat_release(nmat_x);
  gkyl_nmat_release(nmat_y);
  gkyl_nmat_release(nmat_Acu);
  gkyl_nmat_release(nmat_xcu);
  gkyl_nmat_release(nmat_ycu);


}


void test_cu_mat_mm()
{

  struct gkyl_mat *mat_A = gkyl_mat_new(4, 3, 0);
  struct gkyl_mat *mat_x = gkyl_mat_new(3, 2, 0);
  struct gkyl_mat *mat_y = gkyl_mat_new(4, 2, 0);

  struct gkyl_mat *mat_Acu = gkyl_mat_cu_dev_new(4, 3);
  struct gkyl_mat *mat_xcu = gkyl_mat_cu_dev_new(3, 2);
  struct gkyl_mat *mat_ycu = gkyl_mat_cu_dev_new(4, 2);

  // fill each A matrix with its number
  for (size_t j=0; j<mat_A->nc; ++j)
    for (size_t i=0; i<mat_A->nr; ++i)
      gkyl_mat_set(mat_A, i, j, 1);

  // fill each x vecotr with ones
  for (size_t j=0; j<mat_x->nc; ++j)
    for (size_t i=0; i<mat_x->nr; ++i)
      gkyl_mat_set(mat_x, i, j, 1);


  // fill each y matrix with 0
  for (size_t j=0; j<mat_y->nc; ++j)
    for (size_t i=0; i<mat_y->nr; ++i)
      gkyl_mat_set(mat_y, i, j, 0);

  enum gkyl_mat_trans transa = GKYL_NO_TRANS;
  enum gkyl_mat_trans transb = GKYL_NO_TRANS;
  double alpha = 1.0;
  double beta = 0.0;


  // copy to device
  gkyl_mat_copy(mat_Acu, mat_A);
  gkyl_mat_copy(mat_xcu, mat_x);
  gkyl_mat_copy(mat_ycu, mat_y);

  gkyl_mat_mm(alpha, beta, transa, mat_Acu, transb, mat_xcu, mat_ycu, true);

  // copy to host
  gkyl_mat_copy(mat_A, mat_Acu);
  gkyl_mat_copy(mat_x, mat_xcu);
  gkyl_mat_copy(mat_y, mat_ycu);


  // check the expected result
  double expected = 3;
  for (size_t j=0; j<mat_y->nc; ++j){
    for (size_t i=0; i<mat_y->nr; ++i){
      double actual = gkyl_mat_get(mat_y, i, j);
      TEST_CHECK ( expected == actual );
    }
  }
  gkyl_mat_release(mat_A);
  gkyl_mat_release(mat_x);
  gkyl_mat_release(mat_y);
  gkyl_mat_release(mat_Acu);
  gkyl_mat_release(mat_xcu);
  gkyl_mat_release(mat_ycu);


}

void test_cu_nmat_mm()
{

  int n_do = 3;
  struct gkyl_nmat *nmat_A = gkyl_nmat_new(n_do, 4, 3);
  struct gkyl_nmat *nmat_x = gkyl_nmat_new(n_do, 3, 2);
  struct gkyl_nmat *nmat_y = gkyl_nmat_new(n_do, 4, 2);

  struct gkyl_nmat *nmat_Acu = gkyl_nmat_cu_dev_new(n_do, 4, 3);
  struct gkyl_nmat *nmat_xcu = gkyl_nmat_cu_dev_new(n_do, 3, 2);
  struct gkyl_nmat *nmat_ycu = gkyl_nmat_cu_dev_new(n_do, 4, 2);

  // fill each A matrix with its number
  for (size_t n=0; n<nmat_A->num; ++n) {
    struct gkyl_mat m = gkyl_nmat_get(nmat_A, n);
    for (size_t j=0; j<nmat_A->nc; ++j)
      for (size_t i=0; i<nmat_A->nr; ++i)
        gkyl_mat_set(&m, i, j, n);
  }

  // fill each x vecotr with ones
  for (size_t n=0; n<nmat_x->num; ++n) {
    struct gkyl_mat m = gkyl_nmat_get(nmat_x, n);
    for (size_t j=0; j<nmat_x->nc; ++j)
      for (size_t i=0; i<nmat_x->nr; ++i)
        gkyl_mat_set(&m, i, j, 1);
  }

  // fill each y matrix with 0
  for (size_t n=0; n<nmat_y->num; ++n) {
    struct gkyl_mat m = gkyl_nmat_get(nmat_y, n);
    for (size_t j=0; j<nmat_y->nc; ++j)
      for (size_t i=0; i<nmat_y->nr; ++i)
        gkyl_mat_set(&m, i, j, 0);
  }

  enum gkyl_mat_trans transa = GKYL_NO_TRANS;
  double alpha = 1.0;
  double beta = 0.0;


  // copy to device
  gkyl_nmat_copy(nmat_Acu, nmat_A);
  gkyl_nmat_copy(nmat_xcu, nmat_x);
  gkyl_nmat_copy(nmat_ycu, nmat_y);

  gkyl_nmat_mv(alpha, beta, transa, nmat_Acu, nmat_xcu, nmat_ycu);

  // copy to host
  gkyl_nmat_copy(nmat_A, nmat_Acu);
  gkyl_nmat_copy(nmat_x, nmat_xcu);
  gkyl_nmat_copy(nmat_y, nmat_ycu);


  // check the expected result
  for (size_t n=0; n<nmat_y->num; ++n) {
    struct gkyl_mat y = gkyl_nmat_get(nmat_y,n);
    double expected = 3*n;
    for (size_t j=0; j<y.nc; ++j){
      for (size_t i=0; i<y.nr; ++i){
        double actual = gkyl_mat_get(&y, i, j);
        TEST_CHECK ( expected == actual );
      }
    }
  }
  gkyl_nmat_release(nmat_A);
  gkyl_nmat_release(nmat_x);
  gkyl_nmat_release(nmat_y);
  gkyl_nmat_release(nmat_Acu);
  gkyl_nmat_release(nmat_xcu);
  gkyl_nmat_release(nmat_ycu);


}


void test_cu_mat_mm_arrays()
{
  struct gkyl_mat_mm_array_mem *ctest_prob_mem_ho, *ctest_prob_mem_cu; 
  ctest_prob_mem_ho = gkyl_mat_mm_array_mem_new(4, 3, 1.0, 0.0, 
    GKYL_NO_TRANS, GKYL_NO_TRANS, false);
  ctest_prob_mem_cu = gkyl_mat_mm_array_mem_new(4, 3, 1.0, 0.0, 
    GKYL_NO_TRANS, GKYL_NO_TRANS, true);

  struct gkyl_mat *mat_A = ctest_prob_mem_ho->A;
  struct gkyl_array *array_x = gkyl_array_new(GKYL_DOUBLE, 3, 2);
  struct gkyl_array *array_y = gkyl_array_new(GKYL_DOUBLE, 4, 2);

  struct gkyl_mat *mat_Acu = ctest_prob_mem_cu->A;
  struct gkyl_array *array_xcu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3, 2);
  struct gkyl_array *array_ycu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 4, 2);

  // fill each A matrix with its number
  for (size_t j=0; j<mat_A->nc; ++j){
    for (size_t i=0; i<mat_A->nr; ++i){
      double a_val = i*3 + j;
      gkyl_mat_set(mat_A, i, j, a_val);
      double actual = gkyl_mat_get(mat_A, i, j);
      //printf("A(m=%d,n=%d): %1.2e,\n", i, j, actual);
    }
  }

  // fill each x vector with ones
  for (size_t j=0; j<array_x->size; ++j){
    double *x = gkyl_array_fetch(array_x,j);
    for (size_t i=0; i<array_x->ncomp; ++i){
      x[i] = i*2 + j;
      //printf("B(m=%d,n=%d): %1.2e,\n", i, j, x[i]);
    }
  } 

  for (size_t j=0; j<array_y->size; ++j){
    double *y = gkyl_array_fetch(array_y,j);
    for (size_t i=0; i<array_y->ncomp; ++i){
      y[i] = 0.0;
      //printf("C(m=%d,n=%d): %1.2e,\n", i, j, 0.0);
    }
  } 

  // copy to device
  gkyl_mat_copy(mat_Acu, mat_A);
  gkyl_array_copy(array_xcu, array_x);
  gkyl_array_copy(array_ycu, array_y);

  gkyl_mat_mm_array(ctest_prob_mem_cu, array_xcu, array_ycu);

  // copy to host
  gkyl_mat_copy(mat_A, mat_Acu);
  gkyl_array_copy(array_x, array_xcu);
  gkyl_array_copy(array_y, array_ycu);


  // check the expected result
  double expected_array[8] = {10.0, 13.0, 28.0, 40.0, 46.0, 67.0, 64.0, 94.0};
  for (size_t j=0; j<array_y->size; ++j){
    double *y = gkyl_array_fetch(array_y,j);
    for (size_t i=0; i<array_y->ncomp; ++i){
      double actual = y[i];
      double expected = expected_array[i*2 + j];
      //printf("expected: %1.2e, actual: %1.2e\n", expected, actual);
      TEST_CHECK ( expected == actual );
    }
  }
  gkyl_array_release(array_x);
  gkyl_array_release(array_y);
  gkyl_array_release(array_xcu);
  gkyl_array_release(array_ycu);
  gkyl_mat_mm_array_mem_release(ctest_prob_mem_ho);
  gkyl_mat_mm_array_mem_release(ctest_prob_mem_cu);

}

#endif




TEST_LIST = {
  { "mat_base", test_mat_base },
  { "mat_mm_op", test_mat_mm_op },
  { "mat_linsolve", test_mat_linsolve },
  { "nmat_base", test_nmat_base },
  { "nmat_linsolve", test_nmat_linsolve },
  { "nmat_linsolve_pa", test_nmat_linsolve_pa },
  { "mv", test_mat_mv},
  { "nmat_mv", test_nmat_mv},
  { "nmat_mm", test_nmat_mm},
  { "mat_mm_arrays", test_mat_mm_arrays},
#ifdef GKYL_HAVE_CUDA
  { "cu_nmat_base", test_cu_nmat_base },
  { "cu_nmat_linsolve", test_cu_nmat_linsolve },
  { "cu_nmat_linsolve_pa", test_cu_nmat_linsolve_pa },
  { "cu_nmat_mv", test_cu_nmat_mv},
  { "cu_mat_mm", test_cu_mat_mm},
  { "cu_nmat_mm", test_cu_nmat_mm},
  { "cu_mat_mm_arrays", test_cu_mat_mm_arrays},
#endif
  { NULL, NULL },
};

