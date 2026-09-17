#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_array_average gkyl_array_average;

// Input of the new routine is packaged as a struct
typedef struct gkyl_array_average_inp gkyl_array_average_inp;
/**
 * The input structure contains the parameters required to perform a weighted array average:
 * 
 * @param grid Pointer to the computational grid, used to compute the surface element.
 * @param basis Total basis, describes the full dimensionality (ndim), polynomial type, and order.
 * @param basis_avg Subset basis, describes the reduced dimensionality, polynomial type, and order for the output.
 * @param local Full range of the input array, covering all dimensions of the total basis.
 * @param local_avg Reduced range of the output array, covering only the non-averaged dimensions.
 * @param local_avg_ext Extended reduced range of the output array, used only to define the integrated weight.
 * @param weight Pointer to the array containing weight for the averaging process. (set it to NULL for integral)
 * @param avg_dim Flag array to set which dimension is averaged
 * @param use_gpu Boolean flag indicating whether the computation should be performed on a GPU.
 */
struct gkyl_array_average_inp {
  const struct gkyl_rect_grid *grid;
  const struct gkyl_basis basis;
  const struct gkyl_basis basis_avg;
  const struct gkyl_range *local;
  const struct gkyl_range *local_avg;
  const struct gkyl_range *local_avg_ext;
  const struct gkyl_array *weight;
  const int *avg_dim;
  bool use_gpu;
};

/**
 * Create a new updater that computes the weighted average of a gkyl_array.
 * The average can be performed in a subset of the dimensions of the input array and is written generally as:
 * ```math
 *  \int f(x_k) w(x_i) dx^i / \int w(x_k) dx^i
 * ```
 * where $f(x_k)$ is the input array, $w(x_i)$ is the weight.
 * The input are arrays in the configuration space coordinates $\{x_k\}$ and the integral is performed over the dimensions $i$.
 * 
 * If no weighted average is required (weight = NULL), the method will compute the average of the input array, 
 * divided by the volume of the averaging space, i.e.:
 * ```math
 * \int f(x_i) dx^i / \int dx^i
 * ```
 * 
 * If a weight is provided the method will compute the weighted average in two steps:
 *  1. Compute the non-weighted average of the weight by a recursive call with weight = NULL, i.e.:
 *    ```math 
 *    A_w = \int w(x_i) dx^i / \int dx^i
 *    ```
 *  2. Compute the non-weighted average of the array multiplied by the weight, i.e.:
 *    ```math
 *    A_f = \int f(x_i) w(x_i) dx^i / \int dx^i
 *    ```
 *  3. Return the weighted average as:
 *    ```math
 *    A_f / A_w = \int f(x_i) w(x_i) dx^i / \int w(x_i) dx^i
 *    ```
 * 
 * This updater only operates on the local range it is given, it does not communicate.
 * If one of the averaged dimensions is decomposed, the integrals must be summed across
 * the ranks that share it, which the caller does (see gkyl_array_average_advance_range
 * and gkyl_array_average_normalize).
 * 
 * @param inp see gkyl_array_average_inp structure
 */
struct gkyl_array_average*
gkyl_array_average_inew(const struct gkyl_array_average_inp *inp);

/**
 * Create a new updater that computes the weighted average of a gkyl_array.
 * 
 * @param grid Pointer to the computational grid, used to compute the surface element.
 * @param basis Total basis, describes the full dimensionality (ndim), polynomial type, and order.
 * @param basis_avg Subset basis, describes the reduced dimensionality, polynomial type, and order for the output.
 * @param local Full range of the input array, covering all dimensions of the total basis.
 * @param local_avg Reduced range of the output array, covering only the non-averaged dimensions.
 * @param local_avg_ext Extended reduced range of the output array, used only to define the integrated weight.
 * @param weight Pointer to the array containing weight for the averaging process. (set it to NULL for integral)
 * @param avg_dim Flag array to set which dimension is averaged
 * @param use_gpu Boolean flag indicating whether the computation should be performed on a GPU.
 */
struct gkyl_array_average*
gkyl_array_average_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *basis,
  const struct gkyl_basis *basis_avg, const struct gkyl_range *local, const struct gkyl_range *local_avg,
  const struct gkyl_range *local_avg_ext, const struct gkyl_array *weight, const int *avg_dim, bool use_gpu);

/**
 * Compute the array average over the local range. Note: the weight is linked to the updater.
 * This is only the complete average if no averaged dimension is decomposed, since the
 * division by the integrated weight cannot be undone by a later reduction. On a decomposed
 * domain use gkyl_array_average_advance_range and gkyl_array_average_normalize instead.
 *
 * @param up array_average updater.
 * @param fin input gkyl_array
 * @param avgout Output gkyl_array
 */
void gkyl_array_average_advance(const struct gkyl_array_average *up, 
  const struct gkyl_array *fin, struct gkyl_array *avgout);

/**
 * Accumulate the integral of the input (times the weight, if any) over the averaged
 * dimensions of the local range, divided by the volume of the averaging space, i.e.:
 * ```math
 *  \int f(x_k) w(x_i) dx^i / \int dx^i
 * ```
 * The division by the integrated weight is *not* performed, so the result is linear in
 * the input and can be summed across ranks.
 * 
 * The result is written at the indices of out_range, which may be larger than the range
 * of this rank (e.g. the output range of the whole domain), so that the contributions of
 * all ranks can be summed with a single gkyl_comm_allreduce. Cells outside of the range of
 * this rank are left untouched, hence out must be zeroed by the caller beforehand.
 *
 * On a decomposed domain, a weighted average is obtained by
 *  1. calling this routine with the weighted updater to get the numerator,
 *  2. calling this routine with an unweighted updater applied to the weight to get the
 *     denominator (both carry the same volume normalization, which cancels),
 *  3. summing both across the ranks that share the averaged dimensions,
 *  4. calling gkyl_array_average_normalize to perform the division.
 *
 * @param up array_average updater.
 * @param fin input gkyl_array
 * @param out_range Range used to index out (local_avg or a range containing it).
 * @param out Output gkyl_array, accumulated into.
 */
void gkyl_array_average_advance_range(const struct gkyl_array_average *up,
  const struct gkyl_array *fin, const struct gkyl_range *out_range, struct gkyl_array *out);

/**
 * Divide a (possibly reduced) integral by a (possibly reduced) integrated weight, on the
 * local_avg range of the updater. The updater must have been created with a weight.
 *
 * @param up array_average updater.
 * @param wint Integrated weight, as produced by an unweighted updater applied to the weight.
 * @param avgout Array holding the integral on entry and the average on exit.
 */
void gkyl_array_average_normalize(const struct gkyl_array_average *up,
  const struct gkyl_array *wint, struct gkyl_array *avgout);

/**
 * Release memory associated with this updater.
 * @param up array_average updater.
 */
void gkyl_array_average_release(struct gkyl_array_average *up);