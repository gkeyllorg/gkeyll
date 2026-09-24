#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_calc_derived_geo gkyl_calc_derived_geo;

/**
 * Create new updater to compute the derived_geo coefficients
 *
 * @param cbasis Configuration-space basis object.
 * @param grid Configuration-space grid
 * @param node_type 0 for serendipity, 1 for gauss-legendre quadrature nodes
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_calc_derived_geo* gkyl_calc_derived_geo_new(const struct gkyl_basis *cbasis,
  const struct gkyl_rect_grid *grid, int node_type, bool use_gpu);

/**
 * Advance calc_derived_geo (compute the derived_geo coefficients) under the constraint
 * that the representations of these quantities are continuous.
 *
 * @param up calc_derived_geo updater object.
 * @param crange Configuration-space range.
 * @param gFld DG representation of the metric coefficients
 * @param bmagFld DG representation of the magnetic field magnitude
 * @param jFld DG representation of the Jacobian from the configuration space transformation (jacobgeo)
 * @param jinvFld DG representation of the inverse of the Jacobian from the configuration space transformation
 * @param grFld DG representation of duals of the metric g^ij
 * @param biFld DG representation of covariant components of the magnetic field unit vector, b_i
 * @param cmagFld DG representation of C factor in field-line following coordinate formulation
 * @param jtotFld DG representation of total Jacobian (jacobgeo*bmag)
 * @param jtotinvFld DG represenation of inverse of the total Jacobian (1/(jacobgeo*B))
 * @param gxxJFld DG representation of gxx*J (appears in Poisson and diffusion equation)
 * @param gxyJFld DG representation of gxy*J (appears in Poisson equation)
 * @param gyyJFld DG representation of gyy*J (appears in Poisson equation)
 * @param gxzJFld DG representation of gxz*J (can appear in Poisson and diffusion equation depending on formulation)
 * @param eps2Fld DG representation of eps2 = Jg^33 - J/g_33
 */

void gkyl_calc_derived_geo_advance(const gkyl_calc_derived_geo *up, const struct gkyl_range *crange,
  struct gkyl_array *gFld, struct gkyl_array *bmagFld, struct gkyl_array *jFld, struct gkyl_array *jinvFld,
  struct gkyl_array *grFld, struct gkyl_array *biFld, struct gkyl_array *cmagFld, struct gkyl_array *jtotFld, 
  struct gkyl_array *jtotinvFld, struct gkyl_array *gxxJFld,  struct gkyl_array *gxyJFld,
  struct gkyl_array *gyyJFld, struct gkyl_array *gxzJFld, struct gkyl_array *eps2Fld);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_calc_derived_geo_release(gkyl_calc_derived_geo* up);
