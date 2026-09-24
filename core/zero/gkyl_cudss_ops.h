#pragma once
#ifdef GKYL_HAVE_CUDSS

#include <assert.h>
#include <stdbool.h>

#include <gkyl_alloc.h>
#include <gkyl_mat_triples.h>
#include <gkyl_array_ops.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

// Object type
typedef struct gkyl_culinsolver_prob gkyl_culinsolver_prob;

/**
 * Create a new cuDSS object to solve the linear problem
 *   A_i x_i = B_i
 * where i \in {0,1,...,nprob-1}. Each B_i can be a matrix (nrhs>1),
 * meaning that the result x_i is a matrix with the same number of columns.
 * This solver assumes that nrhs is the same for all nprob problems,
 * and that the sparsity pattern of all A_i's is the same.
 */
struct gkyl_culinsolver_prob* gkyl_culinsolver_prob_new(int nprob, int mrow, int ncol, int nrhs);

/**
 * Initialize cuDSS matrix A in Ax=B problem from a list of triples.
 *
 * @param prob cuDSS struct holding arrays used in problem.
 * @param tri (array of) coordinates & values of non-zero entries in A matrix (triplets).
 */
void gkyl_culinsolver_amat_from_triples(struct gkyl_culinsolver_prob *prob, struct gkyl_mat_triples **tri);

/**
 * Update the cuDSS matrix A in Ax=B problem using an array of values (on the device).
 * NOTE: it assumes the sparsity pattern hasn't changed.
 *
 * @param prob cuDSS struct holding arrays used in problem.
 * @param csr_values Array of new values in CSR format, on the device.
 */
void gkyl_culinsolver_amat_update(struct gkyl_culinsolver_prob *prob, double *csr_values);

/**
 * Update the cuDSS matrix A in Ax=B problem using values from a list of triples.
 * NOTE: it assumes the sparsity pattern hasn't changed.
 *
 * @param prob cuDSS struct holding arrays used in problem.
 * @param tri (array of) coordinates & values of non-zero entries in A matrix (triplets).
 */
void gkyl_culinsolver_amat_update_from_triples(struct gkyl_culinsolver_prob *prob, struct gkyl_mat_triples **tri);

/**
 * Initialize right-hand-side cuDSS matrix B in Ax=B problem from a list of
 * triples.
 *
 * @param prob cuDSS struct holding arrays used in problem.
 * @param tri coordinates & values of non-zero entries in B matrix (triplets).
 */
void gkyl_culinsolver_brhs_from_triples(struct gkyl_culinsolver_prob *prob, gkyl_mat_triples *tri);

/**
 * Solve Ax=B problem.
 *
 * @param prob cuDSS struct holding arrays used in problem.
 */
void gkyl_culinsolver_solve(struct gkyl_culinsolver_prob *prob);

/**
 * Copy solution back to host
 *
 * @param prob cuDSS struct holding arrays used in problem.
 */
void gkyl_culinsolver_finish_host(struct gkyl_culinsolver_prob *prob);

/**
 * Synchronize the stream cuDSS ran on.
 *
 * @param prob cuDSS struct holding arrays used in problem.
 */
void gkyl_culinsolver_sync(struct gkyl_culinsolver_prob *prob);

/**
 * Clear the RHS vector by setting all its elements to a value (e.g. 0.).
 *
 * @param prob cuDSS struct holding arrays used in problem.
 * @param val value to set entries of RHS vector to.
 */
void gkyl_culinsolver_clear_rhs(struct gkyl_culinsolver_prob *prob, double val);

/**
 * Clear the array holding (in CSR format) the nonzero values for the LHS matrix.
 *
 * @param prob cuDSS struct holding arrays used in problem.
 * @param val value to set entries of the CSR values array to.
 */
void gkyl_culinsolver_clear_csr_values(struct gkyl_culinsolver_prob *prob, double val);

/**
 * Get a pointer to the element of the RHS vector at a given location.
 *
 * @param prob cuDSS struct holding arrays used in problem.
 * @param loc element we wish to return a pointer to.
 * @return pointer to loc-th element in RHS vector.
 */
double* gkyl_culinsolver_get_rhs_ptr(struct gkyl_culinsolver_prob *prob, long loc);

/**
 * Get a pointer to the element of the solution vector at a given location.
 * Note that the solution vector is the RHS vector after the problem was solved,
 * so this function is equivalent to gkyl_culinsolver_get_rhs_ptr.
 *
 * @param prob cuDSS struct holding arrays used in problem.
 * @param loc element we wish to return a pointer to.
 * @return pointer to loc-th element in solution vector.
 */
double* gkyl_culinsolver_get_sol_ptr(struct gkyl_culinsolver_prob *prob, long loc);

/**
 * Obtain the RHS value at location loc (a linear index into the RHS matrix).
 *
 * @param prob cuDSS struct holding arrays used in problem.
 * @param loc linear index into the RHS flattened array for the desired value.
 * @return RHS value.
 */
double gkyl_culinsolver_get_sol_lin(struct gkyl_culinsolver_prob *prob, long loc);

/**
 * Obtain the a pointer to the array where LHS matrix values at packed in CSR format, 
 * with the pointer offset by @ loc.
 *
 * @param prob cuDSS struct holding arrays used in problem.
 * @param loc linear index into the csr_val array.
 * @return RHS value.
 */
double* gkyl_culinsolver_get_csr_values_ptr(struct gkyl_culinsolver_prob *prob, long loc);

/**
 * Release cuDSS problem
 *
 * @param prob Pointer to cuDSS problem to release.
 */
void gkyl_culinsolver_prob_release(struct gkyl_culinsolver_prob *prob);

// End ifdef GKYL_HAVE_CUDSS statement.
#endif
