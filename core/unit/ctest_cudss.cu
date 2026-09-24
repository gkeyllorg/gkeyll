#ifdef GKYL_HAVE_CUDSS

#include <cudss.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_culinsolver_ops.h>
#include <assert.h>
}

extern "C" {
void test_cudss_simple();
void test_cudss_ops();
void test_cudss_ops_update_amat();
void test_cudss_ops_multiple_rhs();
}

#define checkCUDSS(call, status, msg) \
do { \
    status = call; \
    if (status != CUDSS_STATUS_SUCCESS) { \
        printf("Example FAILED: CUDSS call ended unsuccessfully with status = %d, details: " #msg "\n", status); \
        exit(EXIT_FAILURE); \
    } \
} while(0);

void test_cudss_simple()
{
  // This is meant to replicate the "simple" example in the cuDSS folder of the CUDA samples repo:
  //   https://github.com/NVIDIA/CUDALibrarySamples/blob/master/cuDSS/simple/simple.cpp

  // ---------------------------------------------------------
  // cuDSS example: solving a real linear 5x5 system
  // with a symmetric positive-definite matrix
  // ---------------------------------------------------------
  cudssStatus_t status = CUDSS_STATUS_SUCCESS;

  int n = 5;
  int nnz = 8;
  int nrhs = 1;

  int *csr_offsets_h = NULL;
  int *csr_columns_h = NULL;
  double *csr_values_h = NULL;
  double *x_values_h = NULL, *b_values_h = NULL;

  int *csr_offsets_d = NULL;
  int *csr_columns_d = NULL;
  double *csr_values_d = NULL;
  double *x_values_d = NULL, *b_values_d = NULL;

  csr_offsets_h = (int*)malloc((n + 1) * sizeof(int));
  csr_columns_h = (int*)malloc(nnz * sizeof(int));
  csr_values_h = (double*)malloc(nnz * sizeof(double));
  x_values_h = (double*)malloc(nrhs * n * sizeof(double));
  b_values_h = (double*)malloc(nrhs * n * sizeof(double));

  if (!csr_offsets_h || ! csr_columns_h || !csr_values_h ||
      !x_values_h || !b_values_h) {
      printf("Error: host memory allocation failed\n");
      assert(false);
  }

  /* Initialize host memory for A and b */
  int i = 0;
  csr_offsets_h[i++] = 0;
  csr_offsets_h[i++] = 2;
  csr_offsets_h[i++] = 4;
  csr_offsets_h[i++] = 6;
  csr_offsets_h[i++] = 7;
  csr_offsets_h[i++] = 8;

  i = 0;
  csr_columns_h[i++] = 0; csr_columns_h[i++] = 2;
  csr_columns_h[i++] = 1; csr_columns_h[i++] = 2;
  csr_columns_h[i++] = 2; csr_columns_h[i++] = 4;
  csr_columns_h[i++] = 3;
  csr_columns_h[i++] = 4;

  i = 0;
  csr_values_h[i++] = 4.0; csr_values_h[i++] = 1.0;
  csr_values_h[i++] = 3.0; csr_values_h[i++] = 2.0;
  csr_values_h[i++] = 5.0; csr_values_h[i++] = 1.0;
  csr_values_h[i++] = 1.0;
  csr_values_h[i++] = 2.0;

  /* Note: Right-hand side b is initialized with values which correspond
     to the exact solution vector {1, 2, 3, 4, 5} */
  i = 0;
  b_values_h[i++] = 7.0;
  b_values_h[i++] = 12.0;
  b_values_h[i++] = 25.0;
  b_values_h[i++] = 4.0;
  b_values_h[i++] = 13.0;

  /* Allocate device memory for A, x and b */
  checkCuda(cudaMalloc(&csr_offsets_d, (n + 1) * sizeof(int)));
  checkCuda(cudaMalloc(&csr_columns_d, nnz * sizeof(int)));
  checkCuda(cudaMalloc(&csr_values_d, nnz * sizeof(double)));
  checkCuda(cudaMalloc(&b_values_d, nrhs * n * sizeof(double)));
  checkCuda(cudaMalloc(&x_values_d, nrhs * n * sizeof(double)));

  /* Copy host memory to device for A and b */
  checkCuda(cudaMemcpy(csr_offsets_d, csr_offsets_h, (n + 1) * sizeof(int), cudaMemcpyHostToDevice));
  checkCuda(cudaMemcpy(csr_columns_d, csr_columns_h, nnz * sizeof(int), cudaMemcpyHostToDevice));
  checkCuda(cudaMemcpy(csr_values_d, csr_values_h, nnz * sizeof(double), cudaMemcpyHostToDevice));
  checkCuda(cudaMemcpy(b_values_d, b_values_h, nrhs * n * sizeof(double), cudaMemcpyHostToDevice));

  /* Create a CUDA stream */
  cudaStream_t stream = NULL;
  checkCuda(cudaStreamCreate(&stream));

  /* Creating the cuDSS library handle */
  cudssHandle_t handle;

  checkCUDSS(cudssCreate(&handle), status, "cudssCreate");

  /* (optional) Setting the custom stream for the library handle */
  checkCUDSS(cudssSetStream(handle, stream), status, "cudssSetStream");

  /* Creating cuDSS solver configuration and data objects */
  cudssConfig_t solverConfig;
  cudssData_t solverData;

  checkCUDSS(cudssConfigCreate(&solverConfig), status, "cudssConfigCreate");
  checkCUDSS(cudssDataCreate(handle, &solverData), status, "cudssDataCreate");

  /* Create matrix objects for the right-hand side b and solution x (as dense matrices). */
  cudssMatrix_t x, b;

  int64_t nrows = n, ncols = n;
  int ldb = ncols, ldx = nrows;
  checkCUDSS(cudssMatrixCreateDn(&b, ncols, nrhs, ldb, b_values_d, CUDA_R_64F, CUDSS_LAYOUT_COL_MAJOR),
    status, "cudssMatrixCreateDn for b");
  checkCUDSS(cudssMatrixCreateDn(&x, nrows, nrhs, ldx, x_values_d, CUDA_R_64F, CUDSS_LAYOUT_COL_MAJOR),
    status, "cudssMatrixCreateDn for x");

  /* Create a matrix object for the sparse input matrix. */
  cudssMatrix_t A;
  cudssMatrixType_t mtype     = CUDSS_MTYPE_SPD;
  cudssMatrixViewType_t mview = CUDSS_MVIEW_UPPER;
  cudssIndexBase_t base       = CUDSS_BASE_ZERO;
  checkCUDSS(cudssMatrixCreateCsr(&A, nrows, ncols, nnz, csr_offsets_d, NULL,
    csr_columns_d, csr_values_d, CUDA_R_32I, CUDA_R_64F, mtype, mview,
    base), status, "cudssMatrixCreateCsr");

  /* Symbolic factorization */
  checkCUDSS(cudssExecute(handle, CUDSS_PHASE_ANALYSIS, solverConfig, solverData,
    A, x, b), status, "cudssExecute for analysis");

  /* Factorization */
  checkCUDSS(cudssExecute(handle, CUDSS_PHASE_FACTORIZATION, solverConfig,
    solverData, A, x, b), status, "cudssExecute for factor");

  /* Solving */
  checkCUDSS(cudssExecute(handle, CUDSS_PHASE_SOLVE, solverConfig, solverData,
    A, x, b), status, "cudssExecute for solve");

  /* Destroying opaque objects, matrix wrappers and the cuDSS library handle */
  checkCUDSS(cudssMatrixDestroy(A), status, "cudssMatrixDestroy for A");
  checkCUDSS(cudssMatrixDestroy(b), status, "cudssMatrixDestroy for b");
  checkCUDSS(cudssMatrixDestroy(x), status, "cudssMatrixDestroy for x");
  checkCUDSS(cudssDataDestroy(handle, solverData), status, "cudssDataDestroy");
  checkCUDSS(cudssConfigDestroy(solverConfig), status, "cudssConfigDestroy");
  checkCUDSS(cudssDestroy(handle), status, "cudssHandleDestroy");

  checkCuda(cudaStreamSynchronize(stream));

  /* Print the solution and compare against the exact solution */
  checkCuda(cudaMemcpy(x_values_h, x_values_d, nrhs * n * sizeof(double), cudaMemcpyDeviceToHost));

  // int passed = 1;
  // for (int i = 0; i < n; i++) {
  //   // printf("x[%d] = %1.4f expected %1.4f\n", i, x_values_h[i], double(i+1));
  //   if (fabs(x_values_h[i] - (i + 1)) > 2.e-15)
  //     passed = 0;
  // }

  /* Release the data allocated on the user side */
  free(csr_offsets_h);
  free(csr_columns_h);
  free(csr_values_h);
  free(x_values_h);
  free(b_values_h);
  cudaFree(csr_offsets_d);
  cudaFree(csr_columns_d);
  cudaFree(csr_values_d);
  cudaFree(x_values_d);
  cudaFree(b_values_d);

  // if (status == CUDSS_STATUS_SUCCESS && passed)
  //   printf("Example PASSED\n");
  // else
  //   printf("Example FAILED\n");
}

void test_cudss_ops()
{
  int nfail = 0;  
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

  // Solution is: [-1/32, 11/168, 3/224, 1/16, 11/336].
  GKYL_CU_CHECK( gkyl_compare_double(-1.0/32.0,   gkyl_culinsolver_get_sol_lin(prob,0), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 11.0/168.0, gkyl_culinsolver_get_sol_lin(prob,1), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 3.0/224.0,  gkyl_culinsolver_get_sol_lin(prob,2), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 1.0/16.0,   gkyl_culinsolver_get_sol_lin(prob,3), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 11.0/336.0, gkyl_culinsolver_get_sol_lin(prob,4), 1e-14), &nfail );

  gkyl_culinsolver_prob_release(prob);
}

void test_cudss_ops_update_amat()
{
  int nfail = 0;  
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

  // Create right-hand side matrix B = transpose([1,1,1,1,1]).
  gkyl_mat_triples *triRHS = gkyl_mat_triples_new(m, nrhs);
  gkyl_mat_triples_insert(triRHS, 0, 0, 1.0);
  gkyl_mat_triples_insert(triRHS, 1, 0, 1.0);
  gkyl_mat_triples_insert(triRHS, 2, 0, 1.0);
  gkyl_mat_triples_insert(triRHS, 3, 0, 1.0);
  gkyl_mat_triples_insert(triRHS, 4, 0, 1.0);
  gkyl_culinsolver_brhs_from_triples(prob, triRHS);

  gkyl_culinsolver_solve(prob);
  gkyl_culinsolver_finish_host(prob);

  // Solution is: [-1/32, 11/168, 3/224, 1/16, 11/336].
  GKYL_CU_CHECK( gkyl_compare_double(-1.0/32.0  , gkyl_culinsolver_get_sol_lin(prob,0), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 11.0/168.0, gkyl_culinsolver_get_sol_lin(prob,1), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 3.0/224.0 , gkyl_culinsolver_get_sol_lin(prob,2), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 1.0/16.0  , gkyl_culinsolver_get_sol_lin(prob,3), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 11.0/336.0, gkyl_culinsolver_get_sol_lin(prob,4), 1e-14), &nfail );

  // Now update the LHS matrix. Multiply it by a constant so the solution should be the same as before but divided by that constant.
  double prob_fac = 1.3;

  gkyl_mat_triples_set_rowmaj_order(tri);
  // row 0
  gkyl_mat_triples_insert(tri, 0, 0, s*prob_fac);
  gkyl_mat_triples_insert(tri, 0, 2, u*prob_fac);
  gkyl_mat_triples_insert(tri, 0, 3, u*prob_fac);
  // row 1
  gkyl_mat_triples_insert(tri, 1, 0, l*prob_fac);
  gkyl_mat_triples_insert(tri, 1, 1, u*prob_fac);
  // row 2
  gkyl_mat_triples_insert(tri, 2, 1, l*prob_fac);
  gkyl_mat_triples_insert(tri, 2, 2, p*prob_fac);
  // row 3
  gkyl_mat_triples_insert(tri, 3, 3, e*prob_fac);
  gkyl_mat_triples_insert(tri, 3, 4, u*prob_fac);
  // row 4
  gkyl_mat_triples_insert(tri, 4, 0, l*prob_fac);
  gkyl_mat_triples_insert(tri, 4, 1, l*prob_fac);
  gkyl_mat_triples_insert(tri, 4, 4, r*prob_fac);

  gkyl_culinsolver_amat_update_from_triples(prob, tri_arr);

  // Reset RHS for good measure.
  gkyl_culinsolver_brhs_from_triples(prob, triRHS);

  gkyl_culinsolver_solve(prob);
  gkyl_culinsolver_finish_host(prob);

  // Solution is: (1/prob_fac)*[-1/32, 11/168, 3/224, 1/16, 11/336].
  GKYL_CU_CHECK( gkyl_compare_double((1.0/prob_fac)*(-1.0/32.0  ), gkyl_culinsolver_get_sol_lin(prob,0), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double((1.0/prob_fac)*( 11.0/168.0), gkyl_culinsolver_get_sol_lin(prob,1), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double((1.0/prob_fac)*( 3.0/224.0 ), gkyl_culinsolver_get_sol_lin(prob,2), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double((1.0/prob_fac)*( 1.0/16.0  ), gkyl_culinsolver_get_sol_lin(prob,3), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double((1.0/prob_fac)*( 11.0/336.0), gkyl_culinsolver_get_sol_lin(prob,4), 1e-14), &nfail );

  gkyl_mat_triples_release(tri);
  gkyl_free(tri_arr);
  gkyl_mat_triples_release(triRHS);
  gkyl_culinsolver_prob_release(prob);
}

void test_cudss_ops_multiple_rhs()
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
  GKYL_CU_CHECK( gkyl_compare_double(-1.0/32.0  , gkyl_culinsolver_get_sol_lin(prob,0), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 11.0/168.0, gkyl_culinsolver_get_sol_lin(prob,1), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 3.0/224.0 , gkyl_culinsolver_get_sol_lin(prob,2), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 1.0/16.0  , gkyl_culinsolver_get_sol_lin(prob,3), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 11.0/336.0, gkyl_culinsolver_get_sol_lin(prob,4), 1e-14), &nfail );
  // 2nd problem
  GKYL_CU_CHECK( gkyl_compare_double(-1.0/32.0  , gkyl_culinsolver_get_sol_lin(prob,5), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 11.0/168.0, gkyl_culinsolver_get_sol_lin(prob,6), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 3.0/224.0 , gkyl_culinsolver_get_sol_lin(prob,7), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 1.0/16.0  , gkyl_culinsolver_get_sol_lin(prob,8), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 11.0/336.0, gkyl_culinsolver_get_sol_lin(prob,9), 1e-14), &nfail );
  // 3rd problem
  GKYL_CU_CHECK( gkyl_compare_double(-1.0/32.0  , gkyl_culinsolver_get_sol_lin(prob,10), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 11.0/168.0, gkyl_culinsolver_get_sol_lin(prob,11), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 3.0/224.0 , gkyl_culinsolver_get_sol_lin(prob,12), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 1.0/16.0  , gkyl_culinsolver_get_sol_lin(prob,13), 1e-14), &nfail );
  GKYL_CU_CHECK( gkyl_compare_double( 11.0/336.0, gkyl_culinsolver_get_sol_lin(prob,14), 1e-14), &nfail );

  gkyl_culinsolver_prob_release(prob);
}

void test_cudss_ops_update_amat_multiple_rhs()
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

  gkyl_culinsolver_solve(prob);
  gkyl_culinsolver_finish_host(prob);

  int nfail = 0;
  
  // Solution is: [-1/32, 11/168, 3/224, 1/16, 11/336].
  for (int k=0; k<nrhs; k++) {
    int off = m*k;
    GKYL_CU_CHECK( gkyl_compare_double(-1.0/32.0  , gkyl_culinsolver_get_sol_lin(prob,off+0), 1e-14), &nfail );
    GKYL_CU_CHECK( gkyl_compare_double( 11.0/168.0, gkyl_culinsolver_get_sol_lin(prob,off+1), 1e-14), &nfail );
    GKYL_CU_CHECK( gkyl_compare_double( 3.0/224.0 , gkyl_culinsolver_get_sol_lin(prob,off+2), 1e-14), &nfail );
    GKYL_CU_CHECK( gkyl_compare_double( 1.0/16.0  , gkyl_culinsolver_get_sol_lin(prob,off+3), 1e-14), &nfail );
    GKYL_CU_CHECK( gkyl_compare_double( 11.0/336.0, gkyl_culinsolver_get_sol_lin(prob,off+4), 1e-14), &nfail );
  }

  // Now update the LHS matrix. Multiply it by a constant, and multiply the RHS
  // by 2X that constant, so the solution should be the same but divided by 2.
  double prob_fac = 1.3;
  double *prob_fac_per_rhs = (double *) gkyl_malloc(nrhs*sizeof(double));
  for (int k=0; k<nrhs; k++)
    prob_fac_per_rhs[k] = (k+1)*prob_fac;

  gkyl_mat_triples_set_rowmaj_order(tri);
  // row 0
  gkyl_mat_triples_insert(tri, 0, 0, s*prob_fac);
  gkyl_mat_triples_insert(tri, 0, 2, u*prob_fac);
  gkyl_mat_triples_insert(tri, 0, 3, u*prob_fac);
  // row 1
  gkyl_mat_triples_insert(tri, 1, 0, l*prob_fac);
  gkyl_mat_triples_insert(tri, 1, 1, u*prob_fac);
  // row 2
  gkyl_mat_triples_insert(tri, 2, 1, l*prob_fac);
  gkyl_mat_triples_insert(tri, 2, 2, p*prob_fac);
  // row 3
  gkyl_mat_triples_insert(tri, 3, 3, e*prob_fac);
  gkyl_mat_triples_insert(tri, 3, 4, u*prob_fac);
  // row 4
  gkyl_mat_triples_insert(tri, 4, 0, l*prob_fac);
  gkyl_mat_triples_insert(tri, 4, 1, l*prob_fac);
  gkyl_mat_triples_insert(tri, 4, 4, r*prob_fac);

  // Create the cuSolver linear problem setup.
  gkyl_culinsolver_amat_update_from_triples(prob, tri_arr);

  // Reset RHS.
  for (int k=0; k<nrhs; k++) {
    gkyl_mat_triples_insert(triRHS, 0, k, 1.0*prob_fac_per_rhs[k]);
    gkyl_mat_triples_insert(triRHS, 1, k, 1.0*prob_fac_per_rhs[k]);
    gkyl_mat_triples_insert(triRHS, 2, k, 1.0*prob_fac_per_rhs[k]);
    gkyl_mat_triples_insert(triRHS, 3, k, 1.0*prob_fac_per_rhs[k]);
    gkyl_mat_triples_insert(triRHS, 4, k, 1.0*prob_fac_per_rhs[k]);
  }
  gkyl_culinsolver_brhs_from_triples(prob, triRHS);

  gkyl_culinsolver_solve(prob);
  gkyl_culinsolver_finish_host(prob);
  
  // Solution is: (k+1)*[-1/32, 11/168, 3/224, 1/16, 11/336].
  for (int k=0; k<nrhs; k++) {
    int off = m*k;
    double fac = prob_fac_per_rhs[k]/prob_fac;
    GKYL_CU_CHECK( gkyl_compare_double(fac*(-1.0/32.0  ), gkyl_culinsolver_get_sol_lin(prob,off+0), 1e-14), &nfail );
    GKYL_CU_CHECK( gkyl_compare_double(fac*( 11.0/168.0), gkyl_culinsolver_get_sol_lin(prob,off+1), 1e-14), &nfail );
    GKYL_CU_CHECK( gkyl_compare_double(fac*( 3.0/224.0 ), gkyl_culinsolver_get_sol_lin(prob,off+2), 1e-14), &nfail );
    GKYL_CU_CHECK( gkyl_compare_double(fac*( 1.0/16.0  ), gkyl_culinsolver_get_sol_lin(prob,off+3), 1e-14), &nfail );
    GKYL_CU_CHECK( gkyl_compare_double(fac*( 11.0/336.0), gkyl_culinsolver_get_sol_lin(prob,off+4), 1e-14), &nfail );
  }

  gkyl_free(prob_fac_per_rhs);
  gkyl_mat_triples_release(tri);
  gkyl_free(tri_arr);
  gkyl_mat_triples_release(triRHS);
  gkyl_culinsolver_prob_release(prob);
}

// End ifdef GKYL_HAVE_CUDSS statement.
#endif
