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
 * Compute the array average. Note: the weight is linked to the updater.
 *
 * @param up array_average updater.
 * @param fin input gkyl_array
 * @param avgout Output gkyl_array
 */
void gkyl_array_average_advance(const struct gkyl_array_average *up, 
  const struct gkyl_array *fin, struct gkyl_array *avgout);

/**
 * Release memory associated with this updater.
 * @param up array_average updater.
 */
void gkyl_array_average_release(struct gkyl_array_average *up);