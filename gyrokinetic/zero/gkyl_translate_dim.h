#pragma once

#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_array.h>

// Object type.
typedef struct gkyl_translate_dim gkyl_translate_dim;

/**
 * Create a new updater that translates the DG coefficients of a
 * donor field to those of a target field with a different
 * dimensionality by projecting onto the basis of the target field.
 *
 * This is meant for projecting:
 *   - 1x:
 *       x -> x,y
 *   - 2x:
 *       x,y -> x
 *       x,y -> y
 *       x,y -> x,y,z
 *   - 3x:
 *       x,y,z -> y,z
 *       x,y,z -> x,z
 *       x,y,z -> x,y
 *   - 1x2v:
 *       z,vpar,mu -> 2x2v x,z,vpar,mu
 *   - 1x2v:
 *       z,vpar,mu -> 3x2v x,y,z,vpar,mu
 *   - 2x2v:
 *       x,z,vpar,mu -> 3x2v x,y,z,vpar,mu
 * In downprojecting cases one must choose which dimension to remove
 * via the 'dir' argument, in which case the (logical) variable in
 * that direction is evaluated at -1,0 or 1 within a cell depending
 * on the value of the 'edge' argument.
 *
 * @param cdim_do Configuration space dimension of the donor field.
 * @param basis_do Basis of the donor field.
 * @param cdim_tar Configuration space dimension of the target field.
 * @param basis_tar Basis of the target field.
 * @param dir Direction to remove if cdim_tar < cdim_do.
 * @param edge Where to evaluate the dimension to be removed within a cell if 
 *             cdim_tar < cdim_do (lower boundary, center, or upper boundary).
 * @param use_gpu Whether to run it on the GPU or not.
 */
struct gkyl_translate_dim *gkyl_translate_dim_new(int cdim_do, struct gkyl_basis basis_do,
  int cdim_tar, struct gkyl_basis basis_tar, int dir, enum gkyl_edge_loc edge, bool use_gpu);

/**
 * Run the updater that translates the DG coefficients of a donor
 * field to those of a target field of a different dimensionality.
 *
 * @param up Updater object. 
 * @param rng_do Range of the donor field.
 * @param rng_tar Range of the target field.
 * @param fdo Donor field.
 * @param ncomp Number of scalar fields in fdo/ftar.
 * @param ftar target field.
 */
void gkyl_translate_dim_advance(gkyl_translate_dim *up, const struct gkyl_range *rng_do,
  const struct gkyl_range *rng_tar, const struct gkyl_array *GKYL_RESTRICT fdo, int ncomp,
  struct gkyl_array *GKYL_RESTRICT ftar);

/**
 * Release the memory associated with the translate_dim updater. 
 *
 * @param up translate_dim updater.
 */
void gkyl_translate_dim_release(gkyl_translate_dim *up);
