#ifdef GKYL_HAVE_CUDSS

#include <cudss.h>
#include <cusparse.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_mat_triples.h>
#include <gkyl_array_ops.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_culinsolver_ops.h>
}

#define checkCUDSS(call, status, msg) \
do { \
    status = call; \
    if (status != CUDSS_STATUS_SUCCESS) { \
        fprintf(stderr, "Example FAILED: CUDSS call ended unsuccessfully with status = %d, details: " #msg "\n", status); \
        exit(EXIT_FAILURE); \
    } \
} while(0);

struct gkyl_culinsolver_prob {
  double *rhs_ho, *rhs_cu; // right-hand side vector.
  double *x_ho, *x_cu; // solution vector.
  int nprob; // number of problems to solve.
  int mrow, ncol; // A is a mrow x ncol matrix.
  int nnz; // number of non-zero entries in A.
  int nrhs; // number of columns in B (B is an mrow x nrhs matrix).

  cudaStream_t stream; // CUDA stream cuDSS runs on.
  cudssHandle_t handle; // cuDSS handle.
  cudssConfig_t *solverConfig;
  cudssData_t *solverData;

  cudssMatrix_t *A; // cuDSS object holding the LHS matrix.
  cudssMatrix_t *x, *b; // cuDSS objects holding the unknowns vector and RHS vector.

  // Arrays used to populate the LHS matrix in CSR format.
  double *csr_val_ho;
  double *csr_val_cu;
  int *csr_rowptr_cu, *csr_colind_cu;
};

gkyl_culinsolver_prob*
gkyl_culinsolver_prob_new(int nprob, int mrow, int ncol, int nrhs)
{
  struct gkyl_culinsolver_prob *prob = (struct gkyl_culinsolver_prob*) gkyl_malloc(sizeof(*prob));

  prob->nprob = nprob;
  prob->mrow = mrow;
  prob->ncol = ncol;
  prob->nrhs = nrhs;

  prob->rhs_ho = (double*) gkyl_malloc(nprob * nrhs * mrow * sizeof(double));
  prob->rhs_cu = (double*) gkyl_cu_malloc(nprob * nrhs * mrow * sizeof(double));
  prob->x_ho = (double*) gkyl_malloc(nprob * nrhs * mrow * sizeof(double));
  prob->x_cu = (double*) gkyl_cu_malloc(nprob * nrhs * mrow * sizeof(double));

  cudssStatus_t status = CUDSS_STATUS_SUCCESS;

  /* Create a CUDA stream */
  prob->stream = NULL;
  checkCuda(cudaStreamCreate(&prob->stream));

  /* Creating the cuDSS library handle */
  checkCUDSS(cudssCreate(&prob->handle), status, "cudssCreate");

  /* (optional) Setting the custom stream for the library handle */
  checkCUDSS(cudssSetStream(prob->handle, prob->stream), status, "cudssSetStream");

  /* Creating cuDSS solver configuration and data objects */
  prob->solverConfig = (cudssConfig_t *) gkyl_malloc(nprob * sizeof(cudssConfig_t));
  prob->solverData = (cudssData_t *) gkyl_malloc(nprob * sizeof(cudssData_t));
  for (int i=0; i<nprob; i++) {
    checkCUDSS(cudssConfigCreate(&prob->solverConfig[i]), status, "cudssConfigCreate");
    checkCUDSS(cudssDataCreate(prob->handle, &prob->solverData[i]), status, "cudssDataCreate");
  }

  /* Create matrix objects for the right-hand side b and solution x (as dense matrices). */
  int64_t mrow_64 = mrow, ncol_64 = ncol;
  int ldb = ncol_64, ldx = mrow_64;
  for (int i=0; i<nprob * nrhs * mrow; i++)
    prob->rhs_ho[i] = 1.0;
  gkyl_cu_memcpy(prob->rhs_cu, prob->rhs_ho, nprob * nrhs * mrow * sizeof(double), GKYL_CU_MEMCPY_H2D);
  gkyl_cu_memcpy(prob->x_cu, prob->rhs_ho, nprob * nrhs * mrow * sizeof(double), GKYL_CU_MEMCPY_H2D);

  prob->b = (cudssMatrix_t *) gkyl_malloc(nprob * sizeof(cudssMatrix_t));
  prob->x = (cudssMatrix_t *) gkyl_malloc(nprob * sizeof(cudssMatrix_t));

  for (int i=0; i<nprob; i++) {
    long off = i * nrhs * mrow;
    checkCUDSS(cudssMatrixCreateDn(&prob->b[i], ncol_64, nrhs, ldb, prob->rhs_cu+off, CUDA_R_64F, CUDSS_LAYOUT_COL_MAJOR),
      status, "cudssMatrixCreateDn for b");
    checkCUDSS(cudssMatrixCreateDn(&prob->x[i], mrow_64, nrhs, ldx, prob->x_cu+off, CUDA_R_64F, CUDSS_LAYOUT_COL_MAJOR),
      status, "cudssMatrixCreateDn for x");
  }

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
  prob->csr_val_ho = (double*) gkyl_malloc(prob->nprob*prob->nnz*sizeof(double)); // non-zero matrix elements.
  int *csr_colind = (int*) gkyl_malloc(sizeof(int)*prob->nnz); // col index of entries in csrvalA.
  int *csr_rowptr = (int*) gkyl_malloc(sizeof(int)*(prob->mrow+1)); // 1st entry of each row as index in csrvalA.

  bool *csr_rowptr_assigned = (bool*) gkyl_malloc(sizeof(bool)*prob->mrow);
  for (size_t i=0; i<prob->mrow; i++) csr_rowptr_assigned[i] = false;

  // Sorted (row-major order) keys (linear indices to flattened matrix).
  for (size_t k=0; k<prob->nprob; k++) {
    gkyl_mat_triples_iter *iter = gkyl_mat_triples_iter_new(tri[k]);
    for (size_t i=0; i<prob->nnz; ++i) {
      gkyl_mat_triples_iter_next(iter); // bump iterator.
      struct gkyl_mtriple mt = gkyl_mat_triples_iter_at(iter);
      size_t idx[2] = { mt.row, mt.col };

      prob->csr_val_ho[k*prob->nnz+i] = mt.val;
      if (k==0) {
        csr_colind[i] = idx[1];
        if (!csr_rowptr_assigned[idx[0]]) {
          csr_rowptr[idx[0]] = i;
          csr_rowptr_assigned[idx[0]] = true;
        }
      }
    }
    gkyl_mat_triples_iter_release(iter);
  }
  csr_rowptr[prob->mrow] = prob->nnz;
  gkyl_free(csr_rowptr_assigned);

  // Copy arrays to device.
  prob->csr_val_cu = (double*) gkyl_cu_malloc(prob->nprob*prob->nnz*sizeof(double)); // Non-zero matrix elements.
  prob->csr_colind_cu = (int*) gkyl_cu_malloc(sizeof(int)*prob->nnz); // Col index of entries in csrvalA.
  prob->csr_rowptr_cu = (int*) gkyl_cu_malloc(sizeof(int)*(prob->mrow+1)); // 1st entry of each row as index in csrvalA.
  gkyl_cu_memcpy(prob->csr_val_cu, prob->csr_val_ho, prob->nprob*prob->nnz*sizeof(double), GKYL_CU_MEMCPY_H2D);
  gkyl_cu_memcpy(prob->csr_colind_cu, csr_colind, sizeof(int)*prob->nnz, GKYL_CU_MEMCPY_H2D);
  gkyl_cu_memcpy(prob->csr_rowptr_cu, csr_rowptr, sizeof(int)*(prob->mrow+1), GKYL_CU_MEMCPY_H2D);

  // Create a matrix object for the sparse input matrix.
  cudssStatus_t status = CUDSS_STATUS_SUCCESS;
//  cudssMatrixType_t mtype     = CUDSS_MTYPE_SPD;
  cudssMatrixType_t mtype     = CUDSS_MTYPE_GENERAL;
  cudssMatrixViewType_t mview = CUDSS_MVIEW_UPPER;
  cudssIndexBase_t base       = CUDSS_BASE_ZERO;
  prob->A = (cudssMatrix_t *) gkyl_malloc(prob->nprob * sizeof(cudssMatrix_t));
  for (int i=0; i<prob->nprob; i++) {
    long off = i * prob->nnz;

    checkCUDSS(cudssMatrixCreateCsr(&prob->A[i], prob->mrow, prob->ncol, prob->nnz, prob->csr_rowptr_cu, NULL,
      prob->csr_colind_cu, prob->csr_val_cu+off, CUDA_R_32I, CUDA_R_64F, mtype, mview,
      base), status, "cudssMatrixCreateCsr");
  
    // Symbolic factorization.
    checkCUDSS(cudssExecute(prob->handle, CUDSS_PHASE_ANALYSIS, prob->solverConfig[i], prob->solverData[i],
      prob->A[i], prob->x[i], prob->b[i]), status, "cudssExecute for analysis");
  
    // Factorization.
    checkCUDSS(cudssExecute(prob->handle, CUDSS_PHASE_FACTORIZATION, prob->solverConfig[i],
      prob->solverData[i], prob->A[i], prob->x[i], prob->b[i]), status, "cudssExecute for factor");
  }

  gkyl_free(csr_colind);
  gkyl_free(csr_rowptr);
}

void
gkyl_culinsolver_amat_update_from_triples(struct gkyl_culinsolver_prob *prob, struct gkyl_mat_triples **tri)
{
  // Convert triples to CSR arrays on device.
  // Sorted (row-major order) keys (linear indices to flattened matrix).
  for (size_t k=0; k<prob->nprob; k++) {
    gkyl_mat_triples_iter *iter = gkyl_mat_triples_iter_new(tri[k]);
    for (size_t i=0; i<prob->nnz; ++i) {
      gkyl_mat_triples_iter_next(iter); // bump iterator.
      struct gkyl_mtriple mt = gkyl_mat_triples_iter_at(iter);
      prob->csr_val_ho[k*prob->nnz+i] = mt.val;
    }
    gkyl_mat_triples_iter_release(iter);
  }

  // Copy arrays to device.
  gkyl_cu_memcpy(prob->csr_val_cu, prob->csr_val_ho, prob->nprob*prob->nnz*sizeof(double), GKYL_CU_MEMCPY_H2D);

  for (int i=0; i<prob->nprob; i++) {
    long off = i * prob->nnz;

    cudssStatus_t status = CUDSS_STATUS_SUCCESS;
    // Set matrix values.
    checkCUDSS(cudssMatrixSetValues(prob->A[i], prob->csr_val_cu+off), status, "cudssMatrixSetValues for resetting A.");

    // Factorize.
    checkCUDSS(cudssExecute(prob->handle, CUDSS_PHASE_FACTORIZATION, prob->solverConfig[i],
      prob->solverData[i], prob->A[i], prob->x[i], prob->b[i]), status, "cudssExecute for factor");
  }
}

void
gkyl_culinsolver_amat_update(struct gkyl_culinsolver_prob *prob, double *csr_values)
{
  for (int i=0; i<prob->nprob; i++) {
    long off = i * prob->nnz;

    cudssStatus_t status = CUDSS_STATUS_SUCCESS;
    // Set matrix values.
    checkCUDSS(cudssMatrixSetValues(prob->A[i], csr_values+off), status, "cudssMatrixSetValues for resetting A.");

    // Factorize.
    checkCUDSS(cudssExecute(prob->handle, CUDSS_PHASE_FACTORIZATION, prob->solverConfig[i],
      prob->solverData[i], prob->A[i], prob->x[i], prob->b[i]), status, "cudssExecute for factor");
  }
}

void
gkyl_culinsolver_brhs_from_triples(struct gkyl_culinsolver_prob *prob, gkyl_mat_triples *tri)
{
  long nnz_rhs = gkyl_mat_triples_size(tri);  // Number of non-zero entries in RHS matrix B.

  // Sorted (column-major order) keys (linear indices to flattened matrix).
  gkyl_mat_triples_iter *iter = gkyl_mat_triples_iter_new(tri);
  for (size_t i=0; i<nnz_rhs; i++) {
    gkyl_mat_triples_iter_next(iter); // bump iterator
    struct gkyl_mtriple mt = gkyl_mat_triples_iter_at(iter);
    prob->rhs_ho[i] = mt.val;
  }
  gkyl_mat_triples_iter_release(iter);

  gkyl_cu_memcpy(prob->rhs_cu, prob->rhs_ho, prob->nprob*prob->mrow*prob->nrhs*sizeof(double), GKYL_CU_MEMCPY_H2D);

  cudssStatus_t status = CUDSS_STATUS_SUCCESS;
  for (size_t i=0; i<prob->nprob; i++) {
    long off = i * prob->mrow*prob->nrhs;
    checkCUDSS(cudssMatrixSetValues(prob->b[i], prob->rhs_cu+off),
      status, "cudssMatrixSetValues for setting brhs_from_triples");
  }
}

void
gkyl_culinsolver_solve(struct gkyl_culinsolver_prob *prob)
{
  cudssStatus_t status = CUDSS_STATUS_SUCCESS;

  for (size_t i=0; i<prob->nprob; i++) {
    checkCUDSS(cudssExecute(prob->handle, CUDSS_PHASE_SOLVE, prob->solverConfig[i], prob->solverData[i],
      prob->A[i], prob->x[i], prob->b[i]), status, "cudssExecute for solve");
  }
}

void
gkyl_culinsolver_sync(struct gkyl_culinsolver_prob *prob)
{
  cudaStreamSynchronize(prob->stream);
}

void
gkyl_culinsolver_mat_vec(struct gkyl_culinsolver_prob *prob, const double *x, double *y)
{
  cusparseHandle_t handle;
  cusparseCreate(&handle);
  cusparseSetStream(handle, prob->stream);

  double alpha = 1.0, beta = 0.0;
  void *buf = NULL;
  size_t buf_sz = 0;
  for (int k=0; k<prob->nprob; k++) {
    cusparseSpMatDescr_t matA;
    cusparseCreateCsr(&matA, prob->mrow, prob->ncol, prob->nnz,
      prob->csr_rowptr_cu, prob->csr_colind_cu, prob->csr_val_cu+k*prob->nnz,
      CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F);

    cusparseDnVecDescr_t vecX, vecY;
    cusparseCreateDnVec(&vecX, prob->ncol, (void*) &x[k*prob->ncol], CUDA_R_64F);
    cusparseCreateDnVec(&vecY, prob->mrow, &y[k*prob->mrow], CUDA_R_64F);

    if (k == 0) {
      // All problems have the same sparsity pattern, so one buffer suffices.
      cusparseSpMV_bufferSize(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecX, &beta, vecY,
        CUDA_R_64F, CUSPARSE_SPMV_ALG_DEFAULT, &buf_sz);
      buf = buf_sz > 0? gkyl_cu_malloc(buf_sz) : NULL;
    }

    cusparseSpMV(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecX, &beta, vecY,
      CUDA_R_64F, CUSPARSE_SPMV_ALG_DEFAULT, buf);

    cusparseDestroyDnVec(vecX);
    cusparseDestroyDnVec(vecY);
    cusparseDestroySpMat(matA);
  }
  cudaStreamSynchronize(prob->stream);

  if (buf) gkyl_cu_free(buf);
  cusparseDestroy(handle);
}

void
gkyl_culinsolver_finish_host(struct gkyl_culinsolver_prob *prob)
{
  //cudaStreamSynchronize(prob->stream); // not needed when using blocking stream
  gkyl_cu_memcpy(prob->x_ho, prob->x_cu, prob->nprob*prob->mrow*prob->nrhs*sizeof(double), GKYL_CU_MEMCPY_D2H);
}

void
gkyl_culinsolver_clear_rhs(struct gkyl_culinsolver_prob *prob, double val)
{
  gkyl_cu_memset(prob->rhs_cu, val, prob->nprob*prob->mrow*prob->nrhs*sizeof(double));
}

void
gkyl_culinsolver_clear_csr_values(struct gkyl_culinsolver_prob *prob, double val)
{
  gkyl_cu_memset(prob->csr_val_cu, val, prob->nprob*prob->nnz*sizeof(double));
}

double*
gkyl_culinsolver_get_rhs_ptr(struct gkyl_culinsolver_prob *prob, long loc)
{
  return prob->rhs_cu+loc;
}

double*
gkyl_culinsolver_get_sol_ptr(struct gkyl_culinsolver_prob *prob, long loc)
{
  return prob->x_cu+loc;
}

double*
gkyl_culinsolver_get_csr_values_ptr(struct gkyl_culinsolver_prob *prob, long loc)
{
  return prob->csr_val_cu+loc;
}

double
gkyl_culinsolver_get_sol_lin(struct gkyl_culinsolver_prob *prob, long loc)
{
  return prob->x_ho[loc];
}

void
gkyl_culinsolver_prob_release(struct gkyl_culinsolver_prob *prob)
{
  cudssStatus_t status = CUDSS_STATUS_SUCCESS;

  for (size_t i=0; i<prob->nprob; i++) {
    checkCUDSS(cudssMatrixDestroy(prob->A[i]), status, "cudssMatrixDestroy for A");
    checkCUDSS(cudssMatrixDestroy(prob->b[i]), status, "cudssMatrixDestroy for b");
    checkCUDSS(cudssMatrixDestroy(prob->x[i]), status, "cudssMatrixDestroy for x");
    checkCUDSS(cudssDataDestroy(prob->handle, prob->solverData[i]), status, "cudssDataDestroy");
    checkCUDSS(cudssConfigDestroy(prob->solverConfig[i]), status, "cudssConfigDestroy");
  }
  checkCUDSS(cudssDestroy(prob->handle), status, "cudssHandleDestroy");
  gkyl_free(prob->A);
  gkyl_free(prob->b);
  gkyl_free(prob->x);
  gkyl_free(prob->solverData);
  gkyl_free(prob->solverConfig);

  gkyl_cu_free(prob->csr_colind_cu);
  gkyl_cu_free(prob->csr_rowptr_cu);
  gkyl_cu_free(prob->csr_val_cu);
  gkyl_free(prob->csr_val_ho);

  checkCuda(cudaStreamSynchronize(prob->stream));
  cudaStreamDestroy(prob->stream);

  gkyl_free(prob->rhs_ho);
  gkyl_free(prob->x_ho);
  gkyl_cu_free(prob->rhs_cu);
  gkyl_cu_free(prob->x_cu);

  gkyl_free(prob);
}

// End ifdef GKYL_HAVE_CUDSS statement.
#endif
