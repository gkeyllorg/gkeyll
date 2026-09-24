#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_evalf_def.h>
#include <gkyl_range.h>
#include <gkyl_ref_count.h>
#include <gkyl_rect_grid.h>

// Type for storing preallocating memory needed in various operations
typedef struct gkyl_dg_basis_op_mem gkyl_dg_basis_op_mem;

/**
 * Evaluate the an array containing a DG field
 * at a specific coordinate in the grid. 
 *
 * @param arr Array with DG expansion.
 * @param coord Coordinate to evaluate at.
 * @param basis DG basis object (must be on the GPU if arr is on GPU).
 * @param grid Grid object.
 * @param rng Range object (should contain coord).
 * @param out Evaluation of the array at coord (on GPU if arr is on GPU).
 **/
void gkyl_dg_basis_ops_eval_array_at_coord_comp(const struct gkyl_array *arr, const double *coord,
  const struct gkyl_basis *basis, const struct gkyl_rect_grid *grid, const struct gkyl_range *rng,
  double *out);

/**
 * Evaluate laplacian, given expansion at point in the logical cell
 * (hypercube)
 *
 * @param dir Direction to compute gradient
 * @param z Location to evaluate exansion. z \in [-1,1]^n
 * @param f Expansion coefficients
 * @return Expansion evaluated at z
 */
typedef double (*eval_laplacian_expand)(int dir, const double *z, const double *f);

/**
 * Evaluate mixed partial derivative, given expansion at point in the logical cell
 * (hypercube)
 *
 * @param z Location to evaluate exansion. z \in [-1,1]^n
 * @param f Expansion coefficients
 * @return Expansion evaluated at z
 */
typedef double (*eval_mixedpartial_expand)(const double *z, const double *f);

// Function pointer and context for use in updater that expect evalf_t
// function pointers. This struct essentially wraps the cubic
// interpolation in 1D and 2D.
struct gkyl_basis_ops_evalf {
  void *ctx; // function context
  evalf_t eval_cubic; // function pointer to evaluate the cubic
  evalf_t eval_cubic_wgrad; // function pointer to evaluate the cubic & its gradient with signature
                            // void (*evalf_t)(double t, const double *xn, double *fout, void *ctx);
                            // On return, fout[0], fout[1], and fout[2] are the value, gradient in direction 0,
                            // and gradient in direction 1 of the cubic evaluated at grid coordinates xn.
  evalf_t eval_cubic_wgrad2; // function pointer to evaluate the cubic & its 2nd derivatives
                             // On return, fout[0], fout[1], fout[2], and fout[3] are the value, second derivative
                             // in direction 0, second derivative in direction 1, and mixed partial derivative
                             // of the cubic evaluated at grid coordinates xn.
  eval_laplacian_expand eval_cubic_laplacian; // function pointer to evaluate the laplacian
  eval_mixedpartial_expand eval_cubic_mixedpartial; // function pointer to evaluate the mixed partial
  struct gkyl_ref_count ref_count;   
};  

/**
 * Given values and gradients at the corner of a 1D cell, compute the
 * expansion cofficients of modal cubic (p=3) basis functions. The
 * output @a coeff must be preallocated and have at least 4 elements.
 *
 * IMPORTANT: The input gradients are computed on the unit cell
 * [-1,1]. Hence you MUST properly transform your gradients from
 * gradients on your cell to the unit cell.
 *
 * @param val val[0] and val[1] are the values at left/right edges
 * @param grad grad[0] and grad[1] are the gradients at left/right edges
 * @param coeff On output, the DG expansion coefficients for p=3 basis.
 */
void gkyl_dg_calc_cubic_1d(const double val[2], const double grad[2], double *coeff);

/**
 * Given values and gradients at the corners of a 2D cell, compute the
 * expansion cofficients of modal cubic (p=3) tensor basis functions. The
 * output @a coeff must be preallocated and have at least 16 elements. The nodes
 * are numbered as follows.
 *
 *  2       4
 *  *-------*
 *  |       |
 *  |       |
 *  |       |
 *  *-------*
 *  1       3
 *
 * IMPORTANT: The input gradients are computed on the unit cell [-1,1]
 * x [-1,1]. Hence you MUST properly transform your gradients from
 * gradients on your cell to the unit cell.
 *
 * @param val val[i] is the value at node i
 * @param gradx gradx[i] is the x-derivative at node i
 * @param grady grady[i] is the y-derivative at node i
 * @param gradxy gradxy[i] is the xy- (cross) derivative at node i
 * @param coeff On output, the DG expansion coefficients for p=3 tensor basis.
 */
void gkyl_dg_calc_cubic_2d(const double val[4],
  const double gradx[4], const double grady[4], const double gradxy[4],
  double *coeff);

/**
 * Allocate memory for use in the computing 1D cubic reconstruction
 * from nodal data. Free memory using gkyl_dg_basis_op_mem_release method.
 *
 * @param cells Number of cells in grid
 * @return Newly alloicated memory.
 */
gkyl_dg_basis_op_mem *gkyl_dg_alloc_cubic_1d(int cells);

/**
 * Allocate memory for use in the computing 2D cubic reconstruction
 * from nodal data. Free memory using gkyl_dg_basis_op_mem_release method.
 *
 * @param cells Number of cells in each direction in 2D grid
 * @return Newly allocated memory.
 */
gkyl_dg_basis_op_mem *gkyl_dg_alloc_cubic_2d(int cells[2]);

/**
 * Release memory allocated for use in DG basis ops.
 *
 * @param mem Memory to release
 */
void gkyl_dg_basis_op_mem_release(gkyl_dg_basis_op_mem *mem);

/**
 * Compute cubic expansion from 1D nodal values. Note that the
 * nodal_vals and cubic arrays must not have any ghost cells and must
 * be properly allocated by the caller. Further, mem must also be
 * constructed from a call to gkyl_dg_alloc_cubic_1d.
 *
 * @param mem Memory space to use
 * @param cells Number of cells in grid
 * @param dx Cell spacing
 * @param nodal_vals Array holding nodal values
 * @param cubic On output, DG expansions of cubic
 */
void gkyl_dg_calc_cubic_1d_from_nodal_vals(gkyl_dg_basis_op_mem *mem, int cells, double dx,
  const struct gkyl_array *nodal_vals, struct gkyl_array *cubic);

/**
 * Compute cubic expansion from 2D nodal values. Note that the
 * nodal_vals and cubic arrays must not have any ghost cells and must
 * be properly allocated by the caller. Further, mem must also be
 * constructed from a call to gkyl_dg_alloc_cubic_2d.
 *
 * @param mem Memory space to use
 * @param cells Number of cells along each direction
 * @param dx Cell spacing in each direction
 * @param nodal_vals Array holding nodal values
 * @param cubic On output, DG expansions of cubic
 */
void gkyl_dg_calc_cubic_2d_from_nodal_vals(gkyl_dg_basis_op_mem *mem, int cells[2], double dx[2],
  const struct gkyl_array *nodal_vals, struct gkyl_array *cubic);

/**
 * Create a wrapper to the cubic interpolation from a given set of
 * nodal data in 1D or 2D. The returned context and function pointer
 * can be passed to updaters that expect evalf_t pointers.
 *
 * @param Grid on which interpolation is needed
 * @param nodal_vals Array holding nodal values
 * @return Struct with context and function pointer
 */
struct gkyl_basis_ops_evalf* gkyl_dg_basis_ops_evalf_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_array *nodal_vals);

/**
 * Write the cubic DG representation to specified file.
 *
 * @param evf Interpolator object
 * @param fname Output file name
 * @return If file was written to disk
 */
bool gkyl_dg_basis_ops_evalf_write_cubic(const struct gkyl_basis_ops_evalf *evf, const char *fname);

/**
 * Acquire pointer to memory allocated for cubic interpolation.
 *
 * @param evf Memory to acquire
 * @retrun Pointer to newly acquired memory
 */
struct gkyl_basis_ops_evalf* gkyl_dg_basis_ops_evalf_acquire(const struct gkyl_basis_ops_evalf *evf);

/**
 * Release memory allocated for cubic interpolation.
 *
 * @param evf Memory to release
 */
void gkyl_dg_basis_ops_evalf_release(struct gkyl_basis_ops_evalf *evf);
