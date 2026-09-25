#pragma once

#include "gkyl_mom_type.h"
#include <gkyl_ref_count.h>
#include <stdint.h>


// Define options for projecting a distribution function.

enum gkyl_distribution_projection {
    GKYL_DIST_TYPE_LTE = 0 ,  // Local thermodynamic equilibirum
    GKYL_DIST_TYPE_BIMAX,  // Bi-Maxwellain
    GKYL_DIST_TYPE_KAPPA,  // Kappa 
};

// Forward declare for use in function pointers
struct gkyl_dist_proj_type;

/**
 * Function pointer type to compute the projected distribution.
 */
typedef void (*dist_proj_t)(const struct gkyl_dist_proj_type *dist_proj,
  const double *xc, const double *dx,   
  const int *idx, const double *moms,
  double* out);

struct gkyl_dist_proj_type {
  enum gkyl_distribution_projection dist_id; // selected distribution

  int cdim; // config-space dim
  int pdim; // phase-space dim
  int poly_order; // polynomal order
  int num_config; // number of basis functions in config-space
  int num_phase; // number of basis functions in phase-space
  int num_mom; // number of moments needed for projection

  dist_proj_t proj_kernel; // distribution projection kernel
  struct gkyl_ref_count ref_count; // reference count

  uint32_t flags;
  struct gkyl_dist_proj_type *on_dev; // pointer to itself or device data
};


/**
 * Check if distribution projection type is on device.
 *
 * @param dist_proj distribution projection type to check
 * @return true if dist_proj on device, false otherwise
 */
bool gkyl_dist_proj_type_is_cu_dev(const struct gkyl_dist_proj_type *dist_proj);


/**
 * Acquire pointer to distribution projection object. Delete using the release() method
 *
 * @param dist_proj Distribution projection object to get pointer from.
 * @return acquired object
 */
struct gkyl_dist_proj_type* gkyl_dist_proj_type_acquire(const struct gkyl_dist_proj_type* dist_proj);



/**
 * Delete distribution projection object
 *
 * @param dist_proj Moment object to delete.
 */
void gkyl_dist_proj_type_release(const struct gkyl_dist_proj_type* dist_proj);


/**
 * Project a distribution function using the specified distribution
 * projection object and its required moments.
 *
 * @param dist_proj Distribution projection object
 * @param xc Cell center coordinates
 * @param dx Cell size in each direction
 * @param idx Index into phase-space cell
 * @param moms Input moments required to construct the distribution
 * @param out On output, projected distribution function
 */

void gkyl_dist_proj_type_calc(
  const struct gkyl_dist_proj_type *dist_proj,
  const double *xc, const double *dx,
  const int *idx, const double *moms,
  double *GKYL_RESTRICT out);

/**
 * Get the number of moments required by a distribution projection object.
 *
 * @param dist_proj Distribution projection object
 * @return Number of moments required to construct the distribution
 */
int gkyl_dist_proj_type_num_mom(
  const struct gkyl_dist_proj_type *dist_proj);
