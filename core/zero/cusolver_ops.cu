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

// ..................................
// MF 2022/06/24: unfortunately cusolverRf doesn't work for multiple RHS columns.
//                so for now we loop over the columns.
//                This will hopefully be fixed/optimized in the future so we get
//                better parallelism.
//

struct gkyl_culinsolver_prob {
  double *rhs, *rhs_cu; // right-hand side vector (reused to store the answer x). 
  double *x;
  double *csrvalA_cu;
  int *csrrowptrA_cu, *csrcolindA_cu;
  int nprob; // number of problems to solve.
  int mrow, ncol; // A is a mrow x ncol matrix.
  int nnz; // number of non-zero entries in A.
  int nrhs; // number of columns in B (B is an mrow x nrhs matrix).
  cusolverSpHandle_t cusolverSpH;
  cusparseMatDescr_t A;
  cudaStream_t stream;
  size_t size_internal;

  // cusolverrf objects
  cusparseHandle_t cusparseH;
  csrluInfoHost_t infolu;
  cusolverRfHandle_t cusolverRfH; // Refactorization object.
  // cusolverrf parameters.
  double nzero;
  double nboost;
  double tol;
  double pivot_threshold;
  cusolverRfFactorization_t fact_alg;
  cusolverRfTriangularSolve_t solve_alg;
  double *d_T; // Working space in cusolverRfSolve.
  int *d_P; // P*A*Q^T = L*U
  int *d_Q;

  // for cusolverRfBatch:
  double **rhspointers_cu; // array of pointers to rhs vectors.
  double **csrvalApointers_cu; // array of pointers to LHS A matrices.
};

gkyl_culinsolver_prob*
gkyl_culinsolver_prob_new(int nprob, int mrow, int ncol, int nrhs)
{
  assert((nprob==1) || (nrhs==1));

  struct gkyl_culinsolver_prob *prob = (struct gkyl_culinsolver_prob*) gkyl_malloc(sizeof(*prob));

  prob->nprob = nprob;
  prob->mrow = mrow;
  prob->ncol = ncol;
  prob->nrhs = GKYL_MAX2(nprob,nrhs);

  prob->rhs = (double*) gkyl_malloc(mrow*prob->nrhs*sizeof(double));
  prob->rhs_cu = (double*) gkyl_cu_malloc(mrow*prob->nrhs*sizeof(double));

  if (prob->nrhs > 1) {
    double **rhspointers = (double**) gkyl_malloc(prob->nrhs*sizeof(double*));
    prob->rhspointers_cu = (double**) gkyl_cu_malloc(prob->nrhs*sizeof(double*));
    for (size_t k=0; k<prob->nrhs; k++)
      rhspointers[k] = &prob->rhs_cu[k*mrow];
    gkyl_cu_memcpy(prob->rhspointers_cu, rhspointers, prob->nrhs*sizeof(double*), GKYL_CU_MEMCPY_H2D);
    gkyl_free(rhspointers);
  }

  // create solver
  cusolverSpCreate(&prob->cusolverSpH);

  cudaStreamCreate(&prob->stream); // this stream will block w.r.t. default stream
  //cudaStreamCreateWithFlags(&prob->stream, cudaStreamNonBlocking); // non-blocking w.r.t. default stream
  cusolverSpSetStream(prob->cusolverSpH, prob->stream);

  // create matrix A
  cusparseCreateMatDescr(&prob->A);

  cusparseSetMatType(prob->A, CUSPARSE_MATRIX_TYPE_GENERAL);
  cusparseSetMatIndexBase(prob->A, CUSPARSE_INDEX_BASE_ZERO); 

  cusparseCreate(&prob->cusparseH);
  cusparseSetStream(prob->cusparseH, prob->stream);

  // Create opaque info structure.
  cusolverSpCreateCsrluInfoHost(&prob->infolu);

  // Constants used in cusolverRf:
  //   nzero is the value below which zero pivot is flagged.
  //   nboost is the value which is substitured for zero pivot.
  prob->nzero = 0.0;
  prob->nboost= 0.0;
  // Constants used in cusolverSp:
  //   singularity is -1 if A is invertible under tol
  //   tol determines the condition of singularity
  //   pivot_threshold decides pivoting strategy
  prob->tol = 1.e-16;
  prob->pivot_threshold = 1.0;
  // Constants used in cusolverRf:
  prob->fact_alg = CUSOLVERRF_FACTORIZATION_ALG0; // default
  prob->solve_alg = CUSOLVERRF_TRIANGULAR_SOLVE_ALG1; // default

  return prob;
}

void
gkyl_culinsolver_amat_from_triples(struct gkyl_culinsolver_prob *prob, struct gkyl_mat_triples **tri)
{
  prob->nnz = gkyl_mat_triples_size(tri[0]);
  for (size_t k=0; k<prob->nprob; k++) {
    assert(gkyl_mat_triples_size(tri[k]) == prob->nnz);  // No. of nonzeros must be the same for every problem.
    assert(gkyl_mat_triples_is_rowmaj(tri[k]));  // Triples must be in rowmaj order for cusolver.
  }

  // Convert triples to CSR arrays on device.
  // Use CSR format
  double *csrvalA = (double*) gkyl_malloc(prob->nprob*prob->nnz*sizeof(double)); // non-zero matrix elements.
  int *csrcolindA = (int*) gkyl_malloc(sizeof(int)*prob->nnz); // col index of entries in csrvalA.
  int *csrrowptrA = (int*) gkyl_malloc(sizeof(int)*(prob->mrow+1)); // 1st entry of each row as index in csrvalA.

  bool *csrrowptrA_assigned = (bool*) gkyl_malloc(sizeof(bool)*prob->mrow);
  for (size_t i=0; i<prob->mrow; i++) csrrowptrA_assigned[i] = false;

  // Sorted (row-major order) keys (linear indices to flattened matrix).
  for (size_t k=0; k<prob->nprob; k++) {
    gkyl_mat_triples_iter *iter = gkyl_mat_triples_iter_new(tri[k]);
    for (size_t i=0; i<prob->nnz; ++i) {
      gkyl_mat_triples_iter_next(iter); // bump iterator.
      struct gkyl_mtriple mt = gkyl_mat_triples_iter_at(iter);
      size_t idx[2] = { mt.row, mt.col };
      
      csrvalA[k*prob->nnz+i] = mt.val;
      if (k==0) {
        csrcolindA[i] = idx[1];
        if (!csrrowptrA_assigned[idx[0]]) {
          csrrowptrA[idx[0]] = i;
          csrrowptrA_assigned[idx[0]] = true;
        }
      }
    }
    gkyl_mat_triples_iter_release(iter);
  }
  csrrowptrA[prob->mrow] = prob->nnz;
  gkyl_free(csrrowptrA_assigned);

  // copy arrays to device
  prob->csrvalA_cu = (double*) gkyl_cu_malloc(prob->nprob*prob->nnz*sizeof(double)); // non-zero matrix elements.
  prob->csrcolindA_cu = (int*) gkyl_cu_malloc(sizeof(int)*prob->nnz); // col index of entries in csrvalA.
  prob->csrrowptrA_cu = (int*) gkyl_cu_malloc(sizeof(int)*(prob->mrow+1)); // 1st entry of each row as index in csrvalA.
  gkyl_cu_memcpy(prob->csrvalA_cu, csrvalA, prob->nprob*prob->nnz*sizeof(double), GKYL_CU_MEMCPY_H2D);
  gkyl_cu_memcpy(prob->csrcolindA_cu, csrcolindA, sizeof(int)*prob->nnz, GKYL_CU_MEMCPY_H2D);
  gkyl_cu_memcpy(prob->csrrowptrA_cu, csrrowptrA, sizeof(int)*(prob->mrow+1), GKYL_CU_MEMCPY_H2D);

  double **csrvalApointers;
  if (prob->nrhs > 1) {
    // cusolverRfBatch also needs an array of pointers to
    // the various A matrices (all the same if nprob=1).
    csrvalApointers = (double**) gkyl_malloc(prob->nrhs*sizeof(double*));
    prob->csrvalApointers_cu = (double**) gkyl_cu_malloc(prob->nrhs*sizeof(double*));
    for (size_t k=0; k<prob->nrhs; k++)
      csrvalApointers[k] = prob->nprob == 1? &prob->csrvalA_cu[0] : &prob->csrvalA_cu[k*prob->nnz];
    gkyl_cu_memcpy(prob->csrvalApointers_cu, csrvalApointers, sizeof(double*)*prob->nrhs, GKYL_CU_MEMCPY_H2D);
  }

  // Use CusolverRf

  // reorder to reduce zero fill-in
  // Qreorder = symrcm(A) or Qreroder = symamd(A)
  int *h_Qreorder = (int*)gkyl_malloc(sizeof(int)*prob->ncol);
  // RCM reordering -- seems much slower than others!
//  cusolverSpXcsrsymrcmHost(prob->cusolverSpH, prob->mrow, prob->nnz,
//    prob->A, csrrowptrA, csrcolindA, h_Qreorder);
  // AMD reordering
  cusolverSpXcsrsymamdHost(prob->cusolverSpH, prob->mrow, prob->nnz,
    prob->A, csrrowptrA, csrcolindA, h_Qreorder);
  // MDQ reordering
//  cusolverSpXcsrsymmdqHost(prob->cusolverSpH, prob->mrow, prob->nnz,
//    prob->A, csrrowptrA, csrcolindA, h_Qreorder);
  // METIS reordering
//  cusolverSpXcsrmetisndHost(prob->cusolverSpH, prob->mrow, prob->nnz,
//    prob->A, csrrowptrA, csrcolindA, NULL, h_Qreorder);

  // ............... Compute B = Q*A*Q^T ................... //
  int *h_csrRowIndB = (int*) gkyl_malloc(sizeof(int)*(prob->mrow+1));
  int *h_csrColIndB = (int*) gkyl_malloc(sizeof(int)*prob->nnz);
  memcpy(h_csrRowIndB, csrrowptrA, sizeof(int)*(prob->mrow+1));
  memcpy(h_csrColIndB, csrcolindA, sizeof(int)*prob->nnz);

  size_t size_perm = 0;
  cusolverSpXcsrperm_bufferSizeHost(prob->cusolverSpH, prob->mrow, prob->ncol, prob->nnz,
    prob->A, h_csrRowIndB, h_csrColIndB, h_Qreorder, h_Qreorder, &size_perm);

  void *buffer_cpu = NULL; // working space for permutation (B = Q*A*Q^T) and LU w/ partial pivoting in cusolverSp.
  buffer_cpu = (void*) gkyl_malloc(sizeof(char)*size_perm);

  // h_mapBfromA = Identity
  int *h_mapBfromA  = (int*) gkyl_malloc(sizeof(int)*prob->nnz);
  for(int j = 0 ; j < prob->nnz ; j++) h_mapBfromA[j] = j;
  cusolverSpXcsrpermHost(prob->cusolverSpH, prob->mrow, prob->ncol, prob->nnz, prob->A,
    h_csrRowIndB, h_csrColIndB, h_Qreorder, h_Qreorder, h_mapBfromA, buffer_cpu);

  // B = A( mapBfromA )
  double *h_csrValB = (double*) gkyl_malloc(sizeof(double)*prob->nnz);
  for(int j = 0 ; j < prob->nnz ; j++) h_csrValB[j] = csrvalA[ h_mapBfromA[j] ];

  // ................ Solve A*x = b by LU(B) in cusolverSp ................ //

  // Analyze LU(B) to know structure of Q and R, and upper bound for nnz(L+U).
  cusolverSpXcsrluAnalysisHost(prob->cusolverSpH, prob->mrow, prob->nnz,
    prob->A, h_csrRowIndB, h_csrColIndB, prob->infolu);

  // Workspace for LU(B).
  size_t size_lu = 0; // Size of working space for csrlu.
  cusolverSpDcsrluBufferInfoHost(prob->cusolverSpH, prob->mrow, prob->nnz, prob->A,
    h_csrValB, h_csrRowIndB, h_csrColIndB, prob->infolu, &prob->size_internal, &size_lu);

  if (buffer_cpu) free(buffer_cpu);
  buffer_cpu = (void*)gkyl_malloc(sizeof(char)*size_lu);

  // Compute Ppivot*B = L*U.
  cusolverSpDcsrluFactorHost(prob->cusolverSpH, prob->mrow, prob->nnz, prob->A,
    h_csrValB, h_csrRowIndB, h_csrColIndB, prob->infolu, prob->pivot_threshold, buffer_cpu);

  // Check if the matrix is singular \n");
  int singularity = 0;
  cusolverSpDcsrluZeroPivotHost(prob->cusolverSpH, prob->infolu, prob->tol, &singularity);
  if ( 0 <= singularity){
    fprintf(stderr, "Error: A is not invertible, singularity=%d\n", singularity);
    assert(false);
  }

  // Solve A*x = b, i.e. solve B*(Qx) = Q*b.
  double *h_b = (double*) gkyl_malloc(sizeof(double)*prob->mrow);
  for (size_t i=0; i<prob->mrow; i++) h_b[i] = 1.0;  // Use arbitrary RHS for now.
  double *h_bhat = (double*) gkyl_malloc(sizeof(double)*prob->mrow); // b_hat = Q*b.
  double *h_xhat = (double*) gkyl_malloc(sizeof(double)*prob->ncol); // Q*x_hat = x.
  for(int j = 0 ; j < prob->mrow ; j++) h_bhat[j] = h_b[h_Qreorder[j]]; // b_hat = Q*b
  // B*x_hat = b_hat.
  cusolverSpDcsrluSolveHost(prob->cusolverSpH, prob->mrow, h_bhat, h_xhat, prob->infolu, buffer_cpu);

  // x = Q^T * x_hat
  double *h_x = (double*) gkyl_malloc(sizeof(double)*prob->ncol); // x = A \ b
  for (int j = 0 ; j < prob->mrow ; j++) h_x[h_Qreorder[j]] = h_xhat[j];

  // .............. Extract P, Q, L and U from P*B*Q^T = L*U .............. //

  // L has implicit unit diagonal.
  int nnzL = 0, nnzU = 0;
  cusolverSpXcsrluNnzHost(prob->cusolverSpH, &nnzL, &nnzU, prob->infolu);

  int *h_Plu = (int*) gkyl_malloc(sizeof(int)*prob->mrow);
  int *h_Qlu = (int*) gkyl_malloc(sizeof(int)*prob->ncol);
  double *h_csrValL = (double*)gkyl_malloc(sizeof(double)*nnzL);
  int *h_csrRowIndL = (int*)gkyl_malloc(sizeof(int)*(prob->mrow+1));
  int *h_csrColIndL = (int*)gkyl_malloc(sizeof(int)*nnzL);
  double *h_csrValU = (double*)gkyl_malloc(sizeof(double)*nnzU);
  int *h_csrRowIndU = (int*)gkyl_malloc(sizeof(int)*(prob->mrow+1));
  int *h_csrColIndU = (int*)gkyl_malloc(sizeof(int)*nnzU);

  cusolverSpDcsrluExtractHost(prob->cusolverSpH, h_Plu, h_Qlu, prob->A,
    h_csrValL, h_csrRowIndL, h_csrColIndL, prob->A,
    h_csrValU, h_csrRowIndU, h_csrColIndU, prob->infolu, buffer_cpu);

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
  int *h_P = (int*)gkyl_malloc(sizeof(int)*prob->mrow);
  int *h_Q = (int*)gkyl_malloc(sizeof(int)*prob->ncol);

  // P = Plu*Qreroder.
  // Gather operation, P = Qreorder(Plu).
  for(int j = 0 ; j < prob->mrow ; j++) h_P[j] = h_Qreorder[h_Plu[j]];

  // Q = Qlu*Qreorder.
  // Gather operation, Q = Qreorder(Qlu).
  for(int j = 0 ; j < prob->ncol ; j++) h_Q[j] = h_Qreorder[h_Qlu[j]];

  // ............... Create cusolverRf handle ................ //
  cusolverRfCreate(&prob->cusolverRfH);

  // ............... Set parameters for cusolverRf ................ //
  // Numerical values for checking "zeros" and for boosting.
  cusolverRfSetNumericProperties(prob->cusolverRfH, prob->nzero, prob->nboost);

  // Choose algorithm for refactorization and solve
  cusolverRfSetAlgs(prob->cusolverRfH, prob->fact_alg, prob->solve_alg);

  // Matrix mode: L and U are CSR format, and L has implicit unit diagonal
  cusolverRfSetMatrixFormat(prob->cusolverRfH, CUSOLVERRF_MATRIX_FORMAT_CSR, CUSOLVERRF_UNIT_DIAGONAL_ASSUMED_L);

  // Fast mode for matrix assembling
  cusolverRfSetResetValuesFastMode(prob->cusolverRfH, CUSOLVERRF_RESET_VALUES_FAST_MODE_ON);

  // ............... Assemble P*A*Q = L*U .................. //
  if ((prob->nprob == 1) && (prob->nrhs == 1)) {
    cusolverRfSetupHost(prob->mrow, prob->nnz, csrrowptrA, csrcolindA, csrvalA,
      nnzL, h_csrRowIndL, h_csrColIndL, h_csrValL,
      nnzU, h_csrRowIndU, h_csrColIndU, h_csrValU, h_P, h_Q, prob->cusolverRfH);
  } else {
    for (size_t k=0; k<prob->nrhs; k++)
      csrvalApointers[k] = prob->nprob == 1? &csrvalA[0] : &csrvalA[k*prob->nnz];
    cusolverRfBatchSetupHost(prob->nrhs, prob->mrow, prob->nnz, csrrowptrA, csrcolindA, csrvalApointers,
      nnzL, h_csrRowIndL, h_csrColIndL, h_csrValL,
      nnzU, h_csrRowIndU, h_csrColIndU, h_csrValU, h_P, h_Q, prob->cusolverRfH);
  }

  cudaDeviceSynchronize();

  // ................ Analyze to extract parallelism ............ //
  if (prob->nrhs == 1)
    cusolverRfAnalyze(prob->cusolverRfH);
  else
    cusolverRfBatchAnalyze(prob->cusolverRfH);

  // ................ Import A to cusolverRf ................... //
  prob->d_P = (int*) gkyl_cu_malloc(sizeof(int)*prob->mrow); // P*A*Q^T = L*U
  prob->d_Q = (int*) gkyl_cu_malloc(sizeof(int)*prob->ncol);
  gkyl_cu_memcpy(prob->d_P, h_P, sizeof(int)*prob->mrow, GKYL_CU_MEMCPY_H2D);
  gkyl_cu_memcpy(prob->d_Q, h_Q, sizeof(int)*prob->ncol, GKYL_CU_MEMCPY_H2D);
  
  if (prob->nrhs == 1)
    cusolverRfResetValues(prob->mrow, prob->nnz, prob->csrrowptrA_cu, prob->csrcolindA_cu,
      prob->csrvalA_cu, prob->d_P, prob->d_Q, prob->cusolverRfH);
  else
    cusolverRfBatchResetValues(prob->nrhs, prob->mrow, prob->nnz, prob->csrrowptrA_cu, prob->csrcolindA_cu,
      prob->csrvalApointers_cu, prob->d_P, prob->d_Q, prob->cusolverRfH);
  
  cudaDeviceSynchronize();
  
  //................... Refactorization .................... //
  
  if (prob->nrhs == 1)
    cusolverRfRefactor(prob->cusolverRfH);
  else
    cusolverRfBatchRefactor(prob->cusolverRfH);
  
  cudaDeviceSynchronize();

  if (prob->nrhs == 1)
    prob->d_T = (double*) gkyl_cu_malloc(sizeof(double)*prob->mrow); // Working space in cusolverRfSolve, |d_T| = n * nrhs.
  else
    prob->d_T = (double*) gkyl_cu_malloc(sizeof(double)*prob->mrow*prob->nrhs*2); // Working space in cusolverRfSolve, |d_T| = 2*n*nrhs*batchSize.

  gkyl_free(h_Qreorder);

  gkyl_free(h_csrRowIndB);
  gkyl_free(h_csrColIndB);
  gkyl_free(h_csrValB   );
  gkyl_free(h_mapBfromA );
  
  gkyl_free(h_b);
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

  gkyl_free(csrcolindA);
  gkyl_free(csrrowptrA);
  gkyl_free(csrvalA);
  if (prob->nrhs > 1) gkyl_free(csrvalApointers);

}

void
gkyl_culinsolver_brhs_from_triples(struct gkyl_culinsolver_prob *prob, gkyl_mat_triples *tri)
{
  long nnz_rhs = gkyl_mat_triples_size(tri);  // number of non-zero entries in RHS matrix B
  
  // sorted (column-major order) keys (linear indices to flattened matrix)
  gkyl_mat_triples_iter *iter = gkyl_mat_triples_iter_new(tri);
  for (size_t i=0; i<nnz_rhs; i++) {
    gkyl_mat_triples_iter_next(iter); // bump iterator
    struct gkyl_mtriple mt = gkyl_mat_triples_iter_at(iter);    
    prob->rhs[i] = mt.val;
  }
  gkyl_mat_triples_iter_release(iter);
  
  gkyl_cu_memcpy(prob->rhs_cu, prob->rhs, sizeof(double)*prob->mrow*prob->nrhs, GKYL_CU_MEMCPY_H2D);
}

void
gkyl_culinsolver_solve(struct gkyl_culinsolver_prob *prob)
{
  // MF 2023/05/25: the 1 below is nrhs, and cuSolver docs say only 1 is supported. To me it is not
  // clear whether this means one can only solve 1 system, or whether we can solve multiple systems
  // but each system can only have nrhs=1. I think it's the latter.
  if (prob->nrhs==1)
    cusolverRfSolve(prob->cusolverRfH, prob->d_P, prob->d_Q, 1, prob->d_T, prob->mrow, prob->rhs_cu, prob->mrow);
  else
    cusolverRfBatchSolve(prob->cusolverRfH, prob->d_P, prob->d_Q, 1, prob->d_T, prob->mrow, prob->rhspointers_cu, prob->mrow);
}

void
gkyl_culinsolver_finish_host(struct gkyl_culinsolver_prob *prob)
{
  //cudaStreamSynchronize(prob->stream); // not needed when using blocking stream
  gkyl_cu_memcpy(prob->rhs, prob->rhs_cu, sizeof(double)*prob->mrow*prob->nrhs, GKYL_CU_MEMCPY_D2H);
}

void
gkyl_culinsolver_clear_rhs(struct gkyl_culinsolver_prob *prob, double val)
{
  gkyl_cu_memset(prob->rhs_cu, val, prob->mrow*prob->nrhs*sizeof(double));
}

double*
gkyl_culinsolver_get_rhs_ptr(struct gkyl_culinsolver_prob *prob, long loc)
{
  return prob->rhs_cu+loc;
}

double*
gkyl_culinsolver_get_sol_ptr(struct gkyl_culinsolver_prob *prob, long loc)
{
  return prob->rhs_cu+loc;
}

double
gkyl_culinsolver_get_sol_ij(struct gkyl_culinsolver_prob *prob, long ielement, long jprob)
{
  return prob->rhs[jprob*prob->mrow+ielement];
}

double
gkyl_culinsolver_get_sol_lin(struct gkyl_culinsolver_prob *prob, long loc)
{
  return prob->rhs[loc];
}

void
gkyl_culinsolver_prob_release(struct gkyl_culinsolver_prob *prob)
{
  gkyl_cu_free(prob->rhs_cu);
  gkyl_cu_free(prob->csrcolindA_cu);
  gkyl_cu_free(prob->csrrowptrA_cu);
  gkyl_cu_free(prob->csrvalA_cu);
  if (prob->nrhs > 1) {
    gkyl_cu_free(prob->rhspointers_cu);
    gkyl_cu_free(prob->csrvalApointers_cu);
  }
  gkyl_free(prob->rhs);

  gkyl_cu_free(prob->d_P);
  gkyl_cu_free(prob->d_Q);
  gkyl_cu_free(prob->d_T);
  cusolverRfDestroy(prob->cusolverRfH);
  cusparseDestroy(prob->cusparseH);
  cusolverSpDestroyCsrluInfoHost(prob->infolu);

  cusparseDestroyMatDescr(prob->A);
  cusolverSpDestroy(prob->cusolverSpH);
  cudaStreamDestroy(prob->stream);
  gkyl_free(prob);
}

// End ifndef GKYL_HAVE_CUDSS statement.
#endif
