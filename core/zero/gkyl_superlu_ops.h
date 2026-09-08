#pragma once

#include <gkyl_mat_triples.h>
#include <gkyl_superlu.h>

// Object type
typedef struct gkyl_superlu_prob gkyl_superlu_prob;

/**
 * Create a new SuperLU Ax=B problem.
 * Create a new SuperLU A_i x_i = B_i where i \in {0,1,...,nprob-1},
 * for solving nprob linear problems. This module can be used in two ways:
 *   a) Using nprob=1 and nrhs>=1, so that one can solve nrhs linear problems
 *      with the same left-side matrix A_i for each problem, i.e. B is a matrix
 *      of nprob columns.
 *   b) Using nprob=>1 and nrhs=1, so we are solving nprob separate linear
 *      problems (each with a different left-side Matrix A_i, but with the same
 *      dimensions and sparsity pattern) where each problem only has a
 *      right-side vector with a single column.
 */
struct gkyl_superlu_prob *gkyl_superlu_prob_new(int nprob, int mrow, int ncol, int nrhs);

/**
 * Initialize SuperLU matrix A in Ax=B problem from a list of triples.
 *
 * @param prob SuperLu struct holding arrays used in problem.
 * @param tri coordinates & values of non-zero entries in A matrix (triplets).
 */
void gkyl_superlu_amat_from_triples(struct gkyl_superlu_prob *prob, struct gkyl_mat_triples **tri);

/**
 * Method to print the matrix A to screen.
 */
void gkyl_superlu_print_amat(struct gkyl_superlu_prob *prob);

/**
 * Perform the LU decomposition of the A matrix.
 * The _solve method will use these if they are pre-computed.
 */
void gkyl_superlu_ludecomp(struct gkyl_superlu_prob *prob);

/**
 * Initialize right-hand-side SuperLU matrix B in Ax=B problem from a list of
 * triples.
 *
 * @param prob SuperLu struct holding arrays used in problem.
 * @param tri coordinates & values of non-zero entries in B matrix (triplets).
 */
void gkyl_superlu_brhs_from_triples(struct gkyl_superlu_prob *prob, struct gkyl_mat_triples *tri);

/**
 * Initialize right-hand-side SuperLU matrix B in Ax=B problem from a simple
 * array. This is enough especially when solving a single Ax=B problem (with B
 * a vector) and makes the GPU interface implementation easier.
 *
 * @param prob SuperLu struct holding arrays used in problem.
 * @param bin column array for the right side source.
 */
void gkyl_superlu_brhs_from_array(struct gkyl_superlu_prob *prob, const double *bin);

/**
 * Solve Ax=B problem.
 *
 * @param prob SuperLu struct holding arrays used in problem.
 */
void gkyl_superlu_solve(struct gkyl_superlu_prob *prob);

/**
 * Update the values of the left-hand-side SuperLU matrix A in Ax=B problem from a list of
 * triples. Assumes the sparsity pattern didn't change.
 *
 * @param prob SuperLu struct holding arrays used in problem.
 * @param tri coordinates & values of non-zero entries in A matrix (triplets).
 */
void gkyl_superlu_amat_update_from_triples(
  struct gkyl_superlu_prob *prob, struct gkyl_mat_triples **tri);

/**
 * Obtain the RHS ielement-th value of the jprob-th linear problem.
 * Recall the RHS is a mrow x nrhs matrix.
 *
 * @param linear index into the RHS flattened array for the desired value.
 * @return RHS value.
 */
double gkyl_superlu_get_rhs_ij(struct gkyl_superlu_prob *prob, long ielement, long jprob);

/**
 * Obtain the RHS value at location loc (a linear index into the RHS matrix).
 *
 * @param linear index into the RHS flattened array for the desired value.
 * @return RHS value.
 */
double gkyl_superlu_get_rhs_lin(struct gkyl_superlu_prob *prob, long loc);

double *gkyl_superlu_get_rhs_ptr(struct gkyl_superlu_prob *prob, long loc);

/**
 * Release SuperLU problem
 *
 * @param prob Pointer to SuperLU problem to release.
 */
void gkyl_superlu_prob_release(struct gkyl_superlu_prob *prob);
