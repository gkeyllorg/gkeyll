#ifndef GKYL_HAVE_CUDSS

#include <cusparse.h>
#include <cusolverSp.h>
#include <cusolverRf.h>
#include <cusolverSp_LOWLEVEL_PREVIEW.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_mat_triples.h>
#include <gkyl_array_ops.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_culinsolver_ops.h>
}

extern "C" {
void test_cusolver_qr();
void test_cusolver_rf();
void test_cusolver_ops();
void test_cusolver_ops_multiple_rhs();
void test_cusolver_ops_multiple_prob();
}

void test_cusolver_qr()
{
/*  
 * This is the small 5x5 example used in the Sections 2 and 3 of the 
 * SuperLU Users' Guide to illustrate how to solve a linear problem
 * with cuSolver.
 *
 */
  double   s, u, p, e, r, l;
  int      nrhs, m, n, nnz;
  cusolverSpHandle_t cusolverH = NULL;
  cusparseHandle_t cusparseH = NULL;
  csrqrInfo_t info = NULL;
  cusparseMatDescr_t A = NULL;
  cudaStream_t stream = NULL;
  size_t size_qr = 0;
  size_t size_internal = 0;
  void *buffer_qr = nullptr; // working space for numerical factorization
  

  /* Initialize matrix A. */
  /*  A : matrix([s,0,u,u,0],[l,u,0,0,0],[0,l,p,0,0],[0,0,0,e,u],[l,l,0,0,r]); */
  m = n = 5;
  nnz = 12;
  nrhs = 1;

  s = 19.0; u = 21.0; p = 16.0; e = 5.0; r = 18.0; l = 12.0;
  int cooRowInd[] = {0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4};
  int colInd[] = {0, 2, 3, 0, 1, 1, 2, 3, 4, 0, 1, 4};
  double Aval[] = {s, u, u, l, u, l, p, e, u, l, l, r};

  cusolverSpCreate(&cusolverH);
  cusparseCreate(&cusparseH);

  cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
  cusolverSpSetStream(cusolverH, stream);

  // allocate vectors on device, and copy to device
  int *cooRowInd_cu, *colInd_cu;
  double *Aval_cu;
  cooRowInd_cu = (int*) gkyl_cu_malloc(sizeof(int)*nnz);
  colInd_cu = (int*) gkyl_cu_malloc(sizeof(int)*nnz);
  Aval_cu = (double*) gkyl_cu_malloc(sizeof(double)*nnz);

  gkyl_cu_memcpy(cooRowInd_cu, cooRowInd, sizeof(int)*nnz, GKYL_CU_MEMCPY_H2D);
  gkyl_cu_memcpy(colInd_cu, colInd, sizeof(int)*nnz, GKYL_CU_MEMCPY_H2D);
  gkyl_cu_memcpy(Aval_cu, Aval, sizeof(double)*nnz, GKYL_CU_MEMCPY_H2D);

  // convert to csr format
  int *csrRowPtr_cu;
  csrRowPtr_cu = (int*) gkyl_cu_malloc(sizeof(int)*(m+1));
  cusparseXcoo2csr(cusparseH, cooRowInd_cu, nnz, m, csrRowPtr_cu, CUSPARSE_INDEX_BASE_ZERO);

  int *csrRowPtr = (int*) gkyl_malloc(sizeof(int)*(m+1));
  gkyl_cu_memcpy(csrRowPtr, csrRowPtr_cu, sizeof(int)*(m+1), GKYL_CU_MEMCPY_D2H);

  cusparseCreateMatDescr(&A);

  cusparseSetMatType(A, CUSPARSE_MATRIX_TYPE_GENERAL);
  cusparseSetMatIndexBase(A, CUSPARSE_INDEX_BASE_ZERO); 

  cusolverSpCreateCsrqrInfo(&info);
  
  // allocate rhs vector
  double b[] = {1.0, 1.0, 1.0, 1.0, 1.0};
  double *b_cu;
  b_cu = (double*) gkyl_cu_malloc(sizeof(double)*m);
  gkyl_cu_memcpy(b_cu, b, sizeof(double)*m, GKYL_CU_MEMCPY_H2D);

  // allocate solution vector
  double *x = (double*) gkyl_malloc(sizeof(double)*n);
  double *x_cu = (double*) gkyl_cu_malloc(sizeof(double)*n);

  // symbolic analysis
  cusolverSpXcsrqrAnalysisBatched(cusolverH, m, n, nnz, A, csrRowPtr_cu, colInd_cu, info);

  // prepare working space
  cusolverSpDcsrqrBufferInfoBatched(cusolverH, m, n, nnz, A, Aval_cu,
                                    csrRowPtr_cu, colInd_cu, nrhs, info,
                                    &size_internal, &size_qr);
  
  // numerical factorization
  buffer_qr = (void*) gkyl_cu_malloc(size_qr);
  cusolverSpDcsrqrsvBatched(cusolverH, m, n, nnz, A, Aval_cu, csrRowPtr_cu,
                            colInd_cu, b_cu, x_cu, nrhs, info, buffer_qr);

  cudaStreamSynchronize(stream);

  gkyl_cu_memcpy(x, x_cu, sizeof(double)*n, GKYL_CU_MEMCPY_D2H);

  int nfail = 0;
  /* Solution is: [-1/32, 11/168, 3/224, 1/16, 11/336] */
  GKYL_CU_CHECK( gkyl_compare_double(-1.0/32.0, x[0], 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 11.0/168.0, x[1], 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 3.0/224.0, x[2], 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 1.0/16.0, x[3], 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 11.0/336.0, x[4], 1e-14), &nfail );

  gkyl_cu_free(x_cu);
  gkyl_cu_free(b_cu);
  gkyl_cu_free(buffer_qr);
  gkyl_cu_free(cooRowInd_cu);
  gkyl_cu_free(colInd_cu);
  gkyl_cu_free(Aval_cu);
  gkyl_cu_free(csrRowPtr_cu);

  cusolverSpDestroy(cusolverH);
  cudaStreamDestroy(stream);
}

void test_cusolver_ops()
{
  double s, u, p, e, r, l;
  int    nrhs, m, n;

  /* Initialize matrix A. */
  /*  A : matrix([s,0,u,u,0],[l,u,0,0,0],[0,l,p,0,0],[0,0,0,e,u],[l,l,0,0,r]); */
  m = n = 5;
  nrhs = 1;
  
  s = 19.0; u = 21.0; p = 16.0; e = 5.0; r = 18.0; l = 12.0;
  /*  A : matrix([s,0,u,u,0],[l,u,0,0,0],[0,l,p,0,0],[0,0,0,e,u],[l,l,0,0,r]); */
  struct gkyl_mat_triples **tri_arr = (struct gkyl_mat_triples **) gkyl_malloc(sizeof(struct gkyl_mat_triples *));
  tri_arr[0] = gkyl_mat_triples_new(m, n);
  struct gkyl_mat_triples *tri = tri_arr[0];
  gkyl_mat_triples_set_rowmaj_order(tri);
  // row 0
  gkyl_mat_triples_insert(tri, 0, 0, s);
  gkyl_mat_triples_insert(tri, 0, 2, u);
  gkyl_mat_triples_insert(tri, 0, 3, u);
  // row 1
  gkyl_mat_triples_insert(tri, 1, 0, l);
  gkyl_mat_triples_insert(tri, 1, 1, u);
  // row 2
  gkyl_mat_triples_insert(tri, 2, 1, l);
  gkyl_mat_triples_insert(tri, 2, 2, p);
  // row 3
  gkyl_mat_triples_insert(tri, 3, 3, e);
  gkyl_mat_triples_insert(tri, 3, 4, u);
  // row 4
  gkyl_mat_triples_insert(tri, 4, 0, l);
  gkyl_mat_triples_insert(tri, 4, 1, l);
  gkyl_mat_triples_insert(tri, 4, 4, r);

  // Create the cuSolver linear problem setup.
  gkyl_culinsolver_prob *prob = gkyl_culinsolver_prob_new(1, m, n, nrhs);

  // Allocate the A matrix from triples.
  gkyl_culinsolver_amat_from_triples(prob, tri_arr);
  gkyl_mat_triples_release(tri);
  gkyl_free(tri_arr);

  // Create right-hand side matrix B = transpose([1,1,1,1,1]).
  gkyl_mat_triples *triRHS = gkyl_mat_triples_new(m, nrhs);
  gkyl_mat_triples_insert(triRHS, 0, 0, 1.0);
  gkyl_mat_triples_insert(triRHS, 1, 0, 1.0);
  gkyl_mat_triples_insert(triRHS, 2, 0, 1.0);
  gkyl_mat_triples_insert(triRHS, 3, 0, 1.0);
  gkyl_mat_triples_insert(triRHS, 4, 0, 1.0);
  gkyl_culinsolver_brhs_from_triples(prob, triRHS);
  gkyl_mat_triples_release(triRHS);

  gkyl_culinsolver_solve(prob);
  gkyl_culinsolver_finish_host(prob);

  int nfail = 0;  
  // Solution is: [-1/32, 11/168, 3/224, 1/16, 11/336].
  GKYL_CU_CHECK( gkyl_compare_double(-1.0/32.0,   gkyl_culinsolver_get_sol_lin(prob,0), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 11.0/168.0, gkyl_culinsolver_get_sol_lin(prob,1), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 3.0/224.0,  gkyl_culinsolver_get_sol_lin(prob,2), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 1.0/16.0,   gkyl_culinsolver_get_sol_lin(prob,3), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 11.0/336.0, gkyl_culinsolver_get_sol_lin(prob,4), 1e-14), &nfail );

  gkyl_culinsolver_prob_release(prob);
}

void test_cusolver_ops_multiple_rhs()
{
  double s, u, p, e, r, l;
  int    nrhs, m, n;

  /* Initialize matrix A. */
  /*  A : matrix([s,0,u,u,0],[l,u,0,0,0],[0,l,p,0,0],[0,0,0,e,u],[l,l,0,0,r]); */
  m = n = 5;
  nrhs = 3;
  
  s = 19.0; u = 21.0; p = 16.0; e = 5.0; r = 18.0; l = 12.0;
  /*  A : matrix([s,0,u,u,0],[l,u,0,0,0],[0,l,p,0,0],[0,0,0,e,u],[l,l,0,0,r]); */
  struct gkyl_mat_triples **tri_arr = (struct gkyl_mat_triples **) gkyl_malloc(sizeof(struct gkyl_mat_triples *));
  tri_arr[0] = gkyl_mat_triples_new(m, n);
  struct gkyl_mat_triples *tri = tri_arr[0];
  gkyl_mat_triples_set_rowmaj_order(tri);
  // row 0
  gkyl_mat_triples_insert(tri, 0, 0, s);
  gkyl_mat_triples_insert(tri, 0, 2, u);
  gkyl_mat_triples_insert(tri, 0, 3, u);
  // row 1
  gkyl_mat_triples_insert(tri, 1, 0, l);
  gkyl_mat_triples_insert(tri, 1, 1, u);
  // row 2
  gkyl_mat_triples_insert(tri, 2, 1, l);
  gkyl_mat_triples_insert(tri, 2, 2, p);
  // row 3
  gkyl_mat_triples_insert(tri, 3, 3, e);
  gkyl_mat_triples_insert(tri, 3, 4, u);
  // row 4
  gkyl_mat_triples_insert(tri, 4, 0, l);
  gkyl_mat_triples_insert(tri, 4, 1, l);
  gkyl_mat_triples_insert(tri, 4, 4, r);

  // Create the cuSolver linear problem setup.
  gkyl_culinsolver_prob *prob = gkyl_culinsolver_prob_new(1, m, n, nrhs);

  // Allocate the A matrix from triples.
  gkyl_culinsolver_amat_from_triples(prob, tri_arr);
  gkyl_mat_triples_release(tri);
  gkyl_free(tri_arr);

  // Create right-hand side matrix B = transpose([1,1,1,1,1]).
  gkyl_mat_triples *triRHS = gkyl_mat_triples_new(m, nrhs);
  for (int k=0; k<nrhs; k++) {
    gkyl_mat_triples_insert(triRHS, 0, k, 1.0);
    gkyl_mat_triples_insert(triRHS, 1, k, 1.0);
    gkyl_mat_triples_insert(triRHS, 2, k, 1.0);
    gkyl_mat_triples_insert(triRHS, 3, k, 1.0);
    gkyl_mat_triples_insert(triRHS, 4, k, 1.0);
  }
  gkyl_culinsolver_brhs_from_triples(prob, triRHS);
  gkyl_mat_triples_release(triRHS);

  gkyl_culinsolver_solve(prob);
  gkyl_culinsolver_finish_host(prob);

  int nfail = 0;
  // Solution is: [-1/32, 11/168, 3/224, 1/16, 11/336].
  // 1st problem
  GKYL_CU_CHECK( gkyl_compare_double(-1.0/32.0,   gkyl_culinsolver_get_sol_lin(prob,0), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 11.0/168.0, gkyl_culinsolver_get_sol_lin(prob,1), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 3.0/224.0,  gkyl_culinsolver_get_sol_lin(prob,2), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 1.0/16.0,   gkyl_culinsolver_get_sol_lin(prob,3), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 11.0/336.0, gkyl_culinsolver_get_sol_lin(prob,4), 1e-14), &nfail );
  // 2nd problem
  GKYL_CU_CHECK( gkyl_compare_double(-1.0/32.0,   gkyl_culinsolver_get_sol_lin(prob,5), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 11.0/168.0, gkyl_culinsolver_get_sol_lin(prob,6), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 3.0/224.0,  gkyl_culinsolver_get_sol_lin(prob,7), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 1.0/16.0,   gkyl_culinsolver_get_sol_lin(prob,8), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 11.0/336.0, gkyl_culinsolver_get_sol_lin(prob,9), 1e-14), &nfail );
  // 3rd problem
  GKYL_CU_CHECK( gkyl_compare_double(-1.0/32.0,   gkyl_culinsolver_get_sol_lin(prob,10), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 11.0/168.0, gkyl_culinsolver_get_sol_lin(prob,11), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 3.0/224.0,  gkyl_culinsolver_get_sol_lin(prob,12), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 1.0/16.0,   gkyl_culinsolver_get_sol_lin(prob,13), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 11.0/336.0, gkyl_culinsolver_get_sol_lin(prob,14), 1e-14), &nfail );

  gkyl_culinsolver_prob_release(prob);
}

double superlu_test_answer(double s, double u, double p, double e, double r, double l, int idx) {
  // Solution is: [-1/32, 11/168, 3/224, 1/16, 11/336].
  // for a unit RHS vector and 
  //  s = 19.0; u = 21.0; p = 16.0; e = 5.0; r = 18.0; l = 12.0;
  double sol;
  switch (idx) {
    case 0:
      sol = (e*l*r + e*p*r - l*p*u - e*r*u - p*r*u + p*pow(u,2))/
        (e*pow(l,2)*r + e*p*r*s - pow(l,2)*p*u + l*p*pow(u,2));
      break;
    case 1:
      sol = (r*(-(e*l*p) + e*p*s + e*l*u + l*p*u))/
        (u*(e*pow(l,2)*r + e*p*r*s - pow(l,2)*p*u + l*p*pow(u,2)));
      break;
    case 2:
      sol = -((-(e*pow(l,2)*r) + e*l*r*s + pow(l,2)*r*u - e*r*s*u + pow(l,2)*pow(u,2) -
        l*pow(u,3))/(u*(e*pow(l,2)*r + e*p*r*s - pow(l,2)*p*u + l*p*pow(u,2))));
      break;
    case 3:
      sol = (-(pow(l,2)*p) + pow(l,2)*r + l*p*s + p*r*s + pow(l,2)*u + l*p*u - p*s*u -
      l*pow(u,2))/(e*pow(l,2)*r + e*p*r*s - pow(l,2)*p*u + l*p*pow(u,2));
      break;
    case 4:
      sol = (e*pow(l,2)*p - e*l*p*s - e*pow(l,2)*u - e*l*p*u - pow(l,2)*p*u + e*p*s*u +
        e*l*pow(u,2) + l*p*pow(u,2))/
        (u*(e*pow(l,2)*r + e*p*r*s - pow(l,2)*p*u + l*p*pow(u,2)));
      break;
  }
  return sol;
};

void test_cusolver_ops_multiple_prob()
{
  double s, u, p, e, r, l;
  int    nprob, m, n;

  /* Initialize matrix A. */
  /*  A : matrix([s,0,u,u,0],[l,u,0,0,0],[0,l,p,0,0],[0,0,0,e,u],[l,l,0,0,r]); */
  m = n = 5;
  nprob = 7;

  /*  A : matrix([s,0,u,u,0],[l,u,0,0,0],[0,l,p,0,0],[0,0,0,e,u],[l,l,0,0,r]); */
  struct gkyl_mat_triples **tri_arr = (struct gkyl_mat_triples **) gkyl_malloc(nprob*sizeof(struct gkyl_mat_triples *));
  for (size_t k=0; k<nprob; k++) {
    tri_arr[k] = gkyl_mat_triples_new(m, n);
    struct gkyl_mat_triples *tri = tri_arr[k];
    gkyl_mat_triples_set_rowmaj_order(tri);

    s = 19.0*(k+1)/nprob; u = 21.0*(k+1)/nprob; p = 16.0*(k+1)/nprob; e = 5.0*(k+1)/nprob; r = 18.0*(k+1)/nprob; l = 12.0*(k+1)/nprob;

    // row 0
    gkyl_mat_triples_insert(tri, 0, 0, s);
    gkyl_mat_triples_insert(tri, 0, 2, u);
    gkyl_mat_triples_insert(tri, 0, 3, u);
    // row 1
    gkyl_mat_triples_insert(tri, 1, 0, l);
    gkyl_mat_triples_insert(tri, 1, 1, u);
    // row 2
    gkyl_mat_triples_insert(tri, 2, 1, l);
    gkyl_mat_triples_insert(tri, 2, 2, p);
    // row 3
    gkyl_mat_triples_insert(tri, 3, 3, e);
    gkyl_mat_triples_insert(tri, 3, 4, u);
    // row 4
    gkyl_mat_triples_insert(tri, 4, 0, l);
    gkyl_mat_triples_insert(tri, 4, 1, l);
    gkyl_mat_triples_insert(tri, 4, 4, r);
  }

  // Create the cuSolver linear problem setup.
  gkyl_culinsolver_prob *prob = gkyl_culinsolver_prob_new(nprob, m, n, 1);

  // Allocate the A matrix from triples.
  gkyl_culinsolver_amat_from_triples(prob, tri_arr);
  for (size_t k=0; k<nprob; k++)
    gkyl_mat_triples_release(tri_arr[k]);
  gkyl_free(tri_arr);

  // Create right-hand side matrix B = transpose([1,1,1,1,1]).
  gkyl_mat_triples *triRHS = gkyl_mat_triples_new(m, nprob);
  for (int k=0; k<nprob; k++) {
    gkyl_mat_triples_insert(triRHS, 0, k, 1.0);
    gkyl_mat_triples_insert(triRHS, 1, k, 1.0);
    gkyl_mat_triples_insert(triRHS, 2, k, 1.0);
    gkyl_mat_triples_insert(triRHS, 3, k, 1.0);
    gkyl_mat_triples_insert(triRHS, 4, k, 1.0);
  }
  gkyl_culinsolver_brhs_from_triples(prob, triRHS);
  gkyl_mat_triples_release(triRHS);

  gkyl_culinsolver_solve(prob);
  gkyl_culinsolver_finish_host(prob);

  int nfail = 0;  
  for (int k=0; k<nprob; k++) {
    s = 19.0*(k+1)/nprob; u = 21.0*(k+1)/nprob; p = 16.0*(k+1)/nprob; e = 5.0*(k+1)/nprob; r = 18.0*(k+1)/nprob; l = 12.0*(k+1)/nprob;
    for (int i=0; i<m; i++)
      GKYL_CU_CHECK( gkyl_compare_double( superlu_test_answer(s,u,p,e,r,l,i), gkyl_culinsolver_get_sol_lin(prob,k*5+i), 1e-10), &nfail );
  }

  gkyl_culinsolver_prob_release(prob);
}

void test_cusolver_rf()
{
/*  
 * This is the small 5x5 example used in the Sections 2 and 3 of the 
 * SuperLU Users' Guide to illustrate how to solve a linear problem
 * with cuSolverRF, which allows us to update the coefficients of the
 * A matrix and reuse some of the information from the first LU
 * decomposition.
 * See cuda-samples/7_CUDALibraries/cuSolverRf/cuSolverRf.cpp
 *
 */
  double   s, u, p, e, r, l;

  /* Initialize matrix A. */
  /*  A : matrix([s,0,u,u,0],[l,u,0,0,0],[0,l,p,0,0],[0,0,0,e,u],[l,l,0,0,r]); */
  s = 19.0; u = 21.0; p = 16.0; e = 5.0; r = 18.0; l = 12.0;
  double h_csrValA[] = {s, u, u, l, u, l, p, e, u, l, l, r};
  int h_csrRowIndA[] = {0, 3, 5, 7, 9, 12};
  int h_csrColIndA[] = {0, 2, 3, 0, 1, 1, 2, 3, 4, 0, 1, 4};

  int rowsA = 5; // number of rows of A
  int colsA = rowsA; // number of columns of A
  int nnzA  = sizeof(h_csrValA)/sizeof(h_csrValA[0]); // number of nonzeros of A


  // Constants used in cusolverRf:
  //   nzero is the value below which zero pivot is flagged.
  //   nboost is the value which is substitured for zero pivot.
  double nzero = 0.0;
  double nboost= 0.0;
  // Constants used in cusolverSp:
  //   singularity is -1 if A is invertible under tol
  //   tol determines the condition of singularity
  //   pivot_threshold decides pivoting strategy
  int singularity = 0;
  const double tol = 1.e-14;
  const double pivot_threshold = 1.0;
  // Constants used in cusolverRf:
  const cusolverRfFactorization_t fact_alg = CUSOLVERRF_FACTORIZATION_ALG0; // default
  const cusolverRfTriangularSolve_t solve_alg = CUSOLVERRF_TRIANGULAR_SOLVE_ALG1; // default

  cusolverSpHandle_t cusolverSpH = NULL; // reordering, permutation and 1st LU factorization
  cusparseHandle_t   cusparseH = NULL;   // residual evaluation
  cudaStream_t stream = NULL;
  cusolverSpCreate(&cusolverSpH);
  cusparseCreate(&cusparseH);
  cudaStreamCreate(&stream);
  checkCuda(cudaGetLastError());

  cusolverSpSetStream(cusolverSpH, stream);
  cusparseSetStream(cusparseH, stream);
  checkCuda(cudaGetLastError());

  cusparseMatDescr_t descrA = NULL; // A is a base-0 general matrix
  cusparseCreateMatDescr(&descrA);
  cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL);
  cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ZERO);
  checkCuda(cudaGetLastError());

  // Allocate rhs vector (host only for now).
  double h_b[] = {1.0, 1.0, 1.0, 1.0, 1.0};

  // reorder to reduce zero fill-in
  // Qreorder = symrcm(A) or Qreroder = symamd(A)
  int *h_Qreorder = (int*)gkyl_malloc(sizeof(int)*colsA);
  cusolverSpXcsrsymrcmHost(cusolverSpH, rowsA, nnzA,
    descrA, h_csrRowIndA, h_csrColIndA, h_Qreorder);
//  cusolverSpXcsrsymamdHost(cusolverSpH, rowsA, nnzA,
//    descrA, h_csrRowIndA, h_csrColIndA, h_Qreorder);

  // ............... Compute B = Q*A*Q^T ................... //
  int *h_csrRowIndB = (int*) gkyl_malloc(sizeof(int)*(rowsA+1));
  int *h_csrColIndB = (int*) gkyl_malloc(sizeof(int)*nnzA);
  memcpy(h_csrRowIndB, h_csrRowIndA, sizeof(int)*(rowsA+1));
  memcpy(h_csrColIndB, h_csrColIndA, sizeof(int)*nnzA);

  size_t size_perm = 0;
  cusolverSpXcsrperm_bufferSizeHost(cusolverSpH, rowsA, colsA, nnzA,
    descrA, h_csrRowIndB, h_csrColIndB, h_Qreorder, h_Qreorder, &size_perm);

  void *buffer_cpu = NULL; // working space for permutation (B = Q*A*Q^T) and LU w/ partial pivoting in cusolverSp.
  buffer_cpu = (void*) gkyl_malloc(sizeof(char)*size_perm);

  // h_mapBfromA = Identity
  int *h_mapBfromA  = (int*) gkyl_malloc(sizeof(int)*nnzA);
  for(int j = 0 ; j < nnzA ; j++) h_mapBfromA[j] = j;
  cusolverSpXcsrpermHost(cusolverSpH, rowsA, colsA, nnzA, descrA,
    h_csrRowIndB, h_csrColIndB, h_Qreorder, h_Qreorder, h_mapBfromA, buffer_cpu);

  // B = A( mapBfromA )
  double *h_csrValB = (double*) gkyl_malloc(sizeof(double)*nnzA);
  for(int j = 0 ; j < nnzA ; j++) h_csrValB[j] = h_csrValA[ h_mapBfromA[j] ];

  // ................ Solve A*x = b by LU(B) in cusolverSp ................ //

  // Create opaque info structure.
  csrluInfoHost_t info = NULL; // opaque info structure for LU with parital pivoting
  cusolverSpCreateCsrluInfoHost(&info);

  // Analyze LU(B) to know structure of Q and R, and upper bound for nnz(L+U).
  cusolverSpXcsrluAnalysisHost(cusolverSpH, rowsA, nnzA,
    descrA, h_csrRowIndB, h_csrColIndB, info);

  // Workspace for LU(B).
  size_t size_internal = 0;
  size_t size_lu = 0; // Size of working space for csrlu.
  cusolverSpDcsrluBufferInfoHost(cusolverSpH, rowsA, nnzA, descrA,
    h_csrValB, h_csrRowIndB, h_csrColIndB, info, &size_internal, &size_lu);

  if (buffer_cpu) free(buffer_cpu);
  buffer_cpu = (void*)gkyl_malloc(sizeof(char)*size_lu);

  // Compute Ppivot*B = L*U.
  cusolverSpDcsrluFactorHost(cusolverSpH, rowsA, nnzA, descrA,
    h_csrValB, h_csrRowIndB, h_csrColIndB, info, pivot_threshold, buffer_cpu);

  // Check if the matrix is singular \n");
  cusolverSpDcsrluZeroPivotHost(cusolverSpH, info, tol, &singularity);
  if ( 0 <= singularity){
    fprintf(stderr, "Error: A is not invertible, singularity=%d\n", singularity);
    assert(false);
  }

  // Solve A*x = b, i.e. solve B*(Qx) = Q*b.
  double *h_bhat = (double*) gkyl_malloc(sizeof(double)*rowsA); // b_hat = Q*b.
  double *h_xhat = (double*) gkyl_malloc(sizeof(double)*colsA); // Q*x_hat = x.
  for(int j = 0 ; j < rowsA ; j++) h_bhat[j] = h_b[h_Qreorder[j]]; // b_hat = Q*b
  // B*x_hat = b_hat.
  cusolverSpDcsrluSolveHost(cusolverSpH, rowsA, h_bhat, h_xhat, info, buffer_cpu);

  // x = Q^T * x_hat
  double *h_x = (double*) gkyl_malloc(sizeof(double)*colsA); // x = A \ b
  for (int j = 0 ; j < rowsA ; j++) h_x[h_Qreorder[j]] = h_xhat[j];

  // .............. Extract P, Q, L and U from P*B*Q^T = L*U .............. //

  // L has implicit unit diagonal.
  int nnzL = 0, nnzU = 0;
  cusolverSpXcsrluNnzHost(cusolverSpH, &nnzL, &nnzU, info);

  int *h_Plu = (int*) gkyl_malloc(sizeof(int)*rowsA);
  int *h_Qlu = (int*) gkyl_malloc(sizeof(int)*colsA);
  double *h_csrValL = (double*)gkyl_malloc(sizeof(double)*nnzL);
  int *h_csrRowIndL = (int*)gkyl_malloc(sizeof(int)*(rowsA+1)); 
  int *h_csrColIndL = (int*)gkyl_malloc(sizeof(int)*nnzL);
  double *h_csrValU = (double*)gkyl_malloc(sizeof(double)*nnzU);
  int *h_csrRowIndU = (int*)gkyl_malloc(sizeof(int)*(rowsA+1)); 
  int *h_csrColIndU = (int*)gkyl_malloc(sizeof(int)*nnzU);

  cusolverSpDcsrluExtractHost(cusolverSpH, h_Plu, h_Qlu, descrA,
    h_csrValL, h_csrRowIndL,h_csrColIndL, descrA,
    h_csrValU, h_csrRowIndU,h_csrColIndU, info, buffer_cpu);

  /*  B = Qreorder*A*Qreorder^T
   *  Plu*B*Qlu^T = L*U
   *
   *  (Plu*Qreorder)*A*(Qlu*Qreorder)^T = L*U
   *  
   *  Let P = Plu*Qreroder, Q = Qlu*Qreorder, 
   *  then we have
   *      P*A*Q^T = L*U
   *  which is the fundamental relation in cusolverRf.
   */
  // ............ Form P*A*Q^T = L*U .................. //
  int *h_P = (int*)gkyl_malloc(sizeof(int)*rowsA);
  int *h_Q = (int*)gkyl_malloc(sizeof(int)*colsA);

  // P = Plu*Qreroder.
  // Gather operation, P = Qreorder(Plu).
  for(int j = 0 ; j < rowsA ; j++) h_P[j] = h_Qreorder[h_Plu[j]];

  // Q = Qlu*Qreorder.
  // Gather operation, Q = Qreorder(Qlu).
  for(int j = 0 ; j < colsA ; j++) h_Q[j] = h_Qreorder[h_Qlu[j]];

  // ............... Create cusolverRf handle ................ //
  cusolverRfHandle_t cusolverRfH = NULL; // Refactorization object.
  cusolverRfCreate(&cusolverRfH);

  // ............... Set parameters for cusolverRf ................ //
  // Numerical values for checking "zeros" and for boosting.
  cusolverRfSetNumericProperties(cusolverRfH, nzero, nboost);

  // Choose algorithm for refactorization and solve
  cusolverRfSetAlgs(cusolverRfH, fact_alg, solve_alg);

  // Matrix mode: L and U are CSR format, and L has implicit unit diagonal
  cusolverRfSetMatrixFormat(cusolverRfH, CUSOLVERRF_MATRIX_FORMAT_CSR, CUSOLVERRF_UNIT_DIAGONAL_ASSUMED_L);

  // Fast mode for matrix assembling
  cusolverRfSetResetValuesFastMode(cusolverRfH, CUSOLVERRF_RESET_VALUES_FAST_MODE_ON);

  // ............... Assemble P*A*Q = L*U .................. //
  cusolverRfSetupHost(rowsA, nnzA, h_csrRowIndA, h_csrColIndA, h_csrValA,
    nnzL, h_csrRowIndL, h_csrColIndL, h_csrValL,
    nnzU, h_csrRowIndU, h_csrColIndU, h_csrValU, h_P, h_Q, cusolverRfH);

  cudaDeviceSynchronize();

  // ................ Analyze to extract parallelism ............ //
  cusolverRfAnalyze(cusolverRfH);

  // ................ Import A to cusolverRf ................... //
  int *d_csrRowIndA = (int*) gkyl_cu_malloc(sizeof(int)*(rowsA+1));
  int *d_csrColIndA = (int*) gkyl_cu_malloc(sizeof(int)*nnzA);
  double *d_csrValA = (double*) gkyl_cu_malloc(sizeof(double)*nnzA);
  int *d_P = (int*) gkyl_cu_malloc(sizeof(int)*rowsA); // P*A*Q^T = L*U
  int *d_Q = (int*) gkyl_cu_malloc(sizeof(int)*colsA);
  gkyl_cu_memcpy(d_csrRowIndA, h_csrRowIndA, sizeof(int)*(rowsA+1), GKYL_CU_MEMCPY_H2D);
  gkyl_cu_memcpy(d_csrColIndA, h_csrColIndA, sizeof(int)*nnzA     , GKYL_CU_MEMCPY_H2D);
  gkyl_cu_memcpy(d_csrValA   , h_csrValA   , sizeof(double)*nnzA  , GKYL_CU_MEMCPY_H2D);
  gkyl_cu_memcpy(d_P, h_P, sizeof(int)*rowsA, GKYL_CU_MEMCPY_H2D);
  gkyl_cu_memcpy(d_Q, h_Q, sizeof(int)*colsA, GKYL_CU_MEMCPY_H2D);

  cusolverRfResetValues(rowsA, nnzA, d_csrRowIndA, d_csrColIndA, d_csrValA,
    d_P, d_Q, cusolverRfH);

  cudaDeviceSynchronize();

  //................... Refactorization .................... //

  cusolverRfRefactor(cusolverRfH);

  cudaDeviceSynchronize();

  // .................. Solve A*x = b .................... ..
  double *d_x = (double*) gkyl_cu_malloc(sizeof(double)*colsA);
  gkyl_cu_memcpy(d_x, h_b, sizeof(double)*rowsA, GKYL_CU_MEMCPY_H2D);

  double *d_T = (double*) gkyl_cu_malloc(sizeof(double)*rowsA*1); // Working space in cusolverRfSolve, |d_T| = n * nrhs.
  cusolverRfSolve(cusolverRfH, d_P, d_Q, 1, d_T, rowsA, d_x, rowsA);

  cudaDeviceSynchronize();

  gkyl_cu_memcpy(h_x, d_x, sizeof(double)*colsA, GKYL_CU_MEMCPY_D2H);

  int nfail = 0;
  /* Solution is: [-1/32, 11/168, 3/224, 1/16, 11/336] */
  GKYL_CU_CHECK( gkyl_compare_double(-1.0/32.0  , h_x[0], 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 11.0/168.0, h_x[1], 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 3.0/224.0 , h_x[2], 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 1.0/16.0  , h_x[3], 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 11.0/336.0, h_x[4], 1e-14), &nfail );

  cusolverRfDestroy(cusolverRfH);
  cusolverSpDestroy(cusolverSpH);
  cusparseDestroy(cusparseH);
  cudaStreamDestroy(stream);
  cusparseDestroyMatDescr(descrA);
  cusolverSpDestroyCsrluInfoHost(info);

  gkyl_free(h_Qreorder);

  gkyl_free(h_csrRowIndB);
  gkyl_free(h_csrColIndB);
  gkyl_free(h_csrValB   );
  gkyl_free(h_mapBfromA );

  gkyl_free(h_x);
  gkyl_free(h_xhat);
  gkyl_free(h_bhat);

  gkyl_free(buffer_cpu);

  gkyl_free(h_Plu);
  gkyl_free(h_Qlu);
  gkyl_free(h_csrRowIndL);
  gkyl_free(h_csrColIndL);
  gkyl_free(h_csrValL   );
  gkyl_free(h_csrRowIndU);
  gkyl_free(h_csrColIndU);
  gkyl_free(h_csrValU   );

  gkyl_free(h_P);
  gkyl_free(h_Q);

  gkyl_cu_free(d_csrValA);
  gkyl_cu_free(d_csrRowIndA);
  gkyl_cu_free(d_csrColIndA);
  gkyl_cu_free(d_x);
  gkyl_cu_free(d_P);
  gkyl_cu_free(d_Q);
  gkyl_cu_free(d_T);

}

// End ifndef GKYL_HAVE_CUDSS statement.
#endif
