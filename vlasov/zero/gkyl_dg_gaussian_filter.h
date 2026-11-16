#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_dg_gaussian_filter gkyl_dg_gaussian_filter;

// input packaged as a struct
struct gkyl_dg_gaussian_filter_inp {
  const struct gkyl_rect_grid *conf_grid; // Configuration-space grid
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions
  const struct gkyl_range *conf_range; // Configuration-space range
  const struct gkyl_range
    *conf_range_ext; // Extended configuration-space range (for internal memory allocations)
  bool extend_filter; // bool for whether to extend filter to 5-cell filter
  bool use_gpu; // bool for gpu useage
};

/**
 * Create new updater to apply a Gaussian filter to an input field. 
 * By default, performs a three-cell filter where the value of output 
 * at the ith configuration space quadrature point is:
 * conf_arr_filter[i] = sum_j [conf_arr_l[j]*exp(-(xmu_i - xmu_l_j)/2 dx^2) +
 * conf_arr_c[j]*exp(-(xmu_i - xmu_c_j)/2 dx^2) + conf_arr_r[j]*exp(-(xmu_i - xmu_r_j)/2 dx^2) ]
 * where l,c,r are the lower, central, and upper cells we are fetching quadrature point information 
 * from for performing the Gaussian filter, and xmu are the quadrature point abscissas.
 * 
 * Optionally extends the filter to be a five-cell filter, accumulating the contribution of the
 * filter from the lower-lower (i-2 index) and upper-upper (i+2) index for every cell in the
 * interior that does not abut the boundary (i.e., not the skin cells). 
 *
 * @param inp Input parameters defined in gkyl_dg_gaussian_filter_inp struct.
 * @return New updater pointer.
 */
struct gkyl_dg_gaussian_filter *gkyl_dg_gaussian_filter_inew(
  const struct gkyl_dg_gaussian_filter_inp *inp
);

/**
 * Compute the Gaussian filter on input conf_arry array. 
 *
 * @param up Project on basis updater to run
 * @param conf_range Configuration-space range
 * @param conf_arr Input array to compute Gaussian filter of. Result of Gaussian filter
 *                 stored in this array as output. 
 */
void gkyl_dg_gaussian_filter_advance(
  gkyl_dg_gaussian_filter *up, const struct gkyl_range *conf_range, struct gkyl_array *conf_arr
);

/**
 * Host-side wrapper for computing Gaussian filter on device.
 */
void gkyl_dg_gaussian_filter_advance_cu(
  gkyl_dg_gaussian_filter *up, const struct gkyl_range *conf_range, struct gkyl_array *conf_arr
);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_gaussian_filter_release(gkyl_dg_gaussian_filter *up);
