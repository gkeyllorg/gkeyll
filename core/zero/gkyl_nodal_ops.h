#pragma once

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_basis.h>

// Struct definition. Used to store nodal values 
struct gkyl_nodal_ops {
  struct gkyl_array *nodes;
  int poly_order;
};

/**
 * Create new nodal_ops struct for transforming between nodal and modal representations.
 *
 * @param cbasis Configuration-space basis
 * @param grid Configuration-space grid
 * @param use_gpu Boolean for whether nodes are stored on device
 * Returns pointer to gkyl_nodal_ops struct.
 */
struct gkyl_nodal_ops* gkyl_nodal_ops_new(const struct gkyl_basis *cbasis, const struct gkyl_rect_grid *grid, bool use_gpu);

/**
 * Transform nodal representation to modal representation
 *
 * @param nodal_ops Nodal operations struct
 * @param cbasis Configuration-space basis
 * @param grid Configuration-space grid
 * @param nrange Nodal range
 * @param update_range Configuration-space range on which we are operating
 * @param num_comp Number of components
 * @param nodal_fld Input nodal representation
 * @param modal_fld Output modal representation
 * @param use_quad flag to indicate whether to use interior quadrature points
 * if false, serendipity nodes are used
 */
void gkyl_nodal_ops_n2m(const struct gkyl_nodal_ops *nodal_ops, 
  const struct gkyl_basis *cbasis, const struct gkyl_rect_grid *grid, 
  const struct gkyl_range *nrange, const struct gkyl_range *update_range, int num_comp, 
  const struct gkyl_array *nodal_fld, struct gkyl_array *modal_fld, bool use_quad);

/**
 * Transform nodal representation to modal representation
 * at surface quadraure points
 *
 * @param nodal_ops Nodal operations struct
 * @param cbasis Configuration-space basis
 * @param grid Configuration-space grid
 * @param nrange Nodal range
 * @param update_range Configuration-space range on which we are operating
 * @param num_comp Number of components
 * @param nodal_fld Input nodal representation
 * @param modal_fld Output modal representation
 * @param dir direction of surfaces
 * if false, serendipity nodes are used
 */
void 
gkyl_nodal_ops_n2m_surface(const struct gkyl_nodal_ops *nodal_ops, 
  const struct gkyl_basis *cbasis, const struct gkyl_rect_grid *grid, 
  const struct gkyl_range *nrange, const struct gkyl_range *update_range, int num_comp, 
  const struct gkyl_array *nodal_fld, struct gkyl_array *modal_fld, int dir) ;

/**
 * Transform modal representation to nodal representation
 * at surface quadraure points
 *
 * @param nodal_ops Nodal operations struct
 * @param cbasis Configuration-space basis
 * @param grid Configuration-space grid
 * @param nrange Nodal range
 * @param update_range Configuration-space range on which we are operating
 * @param num_comp Number of components
 * @param nodal_fld output nodal representation
 * @param modal_fld input modal representation
 * @param dir direction of surfaces
 * if false, serendipity nodes are used
 */
void 
gkyl_nodal_ops_m2n_surface(const struct gkyl_nodal_ops *nodal_ops, 
  const struct gkyl_basis *cbasis, const struct gkyl_rect_grid *grid, 
  const struct gkyl_range *nrange, const struct gkyl_range *update_range, int num_comp, 
  struct gkyl_array *nodal_fld, const struct gkyl_array *modal_fld, int dir) ;

/**
 * Transform modal representation to nodal representation
 *
 * @param nodal_ops Nodal operations struct
 * @param cbasis Configuration-space basis
 * @param grid Configuration-space grid
 * @param nrange Nodal range
 * @param update_range Configuration-space range on which we are operating
 * @param num_comp Number of components
 * @param nodal_fld Output nodal representation
 * @param modal_fld Input modal representation
 * @param use_quad flag to indicate whether to use interior quadrature points
 * if false, serendipity nodes are used
 */
void gkyl_nodal_ops_m2n(const struct gkyl_nodal_ops *nodal_ops, 
  const struct gkyl_basis *cbasis, const struct gkyl_rect_grid *grid, 
  const struct gkyl_range *nrange, const struct gkyl_range *update_range, int num_comp, 
  struct gkyl_array *nodal_fld, const struct gkyl_array *modal_fld, bool use_quad);

/**
 * Transform modal representation of dim = d-1 of to nodal representation
 * along one slice (last dim removed) of the nodal array of dim = d
 *
 * @param nodal_ops Nodal operations struct should have dim = d-1
 * @param deflated_cbasis Configuration-space basis with dim = d-1
 * @param deflated_grid Configuration-space grid
 * @param nrange Nodal range with dim = d
 * @param deflated nrange Nodal range with dim = d-1. Should not be a subrange of nrange
 * @param deflated_update_range Configuration-space range on which we are operating with dim = d-1
 * @param num_comp Number of components
 * @param nodal_fld Output nodal representation with dim = d
 * @param deflated_modal_fld Input modal representation with dim = d-1
 */

void gkyl_nodal_ops_m2n_deflated(const struct gkyl_nodal_ops *nodal_ops, 
  const struct gkyl_basis *deflated_cbasis, const struct gkyl_rect_grid *deflated_grid, 
  const struct gkyl_range *nrange, const struct gkyl_range *deflated_nrange, const struct gkyl_range *deflated_update_range, int num_comp, 
  struct gkyl_array *nodal_fld, const struct gkyl_array *deflated_modal_fld, int extra_idx) ;

/**
 * Delete pointer to gkyl_nodal_ops struct.
 *
 * @param up Struct to delete.
 */
void gkyl_nodal_ops_release(struct gkyl_nodal_ops *up);

/**
 * Host-side wrappers for nodal-modal transformation operations on device
 */

void gkyl_nodal_ops_n2m_cu(const struct gkyl_nodal_ops *nodal_ops, 
  const struct gkyl_basis *cbasis, const struct gkyl_rect_grid *grid, 
  const struct gkyl_range *nrange, const struct gkyl_range *update_range, int num_comp, 
  const struct gkyl_array *nodal_fld, struct gkyl_array *modal_fld);

void gkyl_nodal_ops_m2n_cu(const struct gkyl_nodal_ops *nodal_ops, 
  const struct gkyl_basis *cbasis, const struct gkyl_rect_grid *grid, 
  const struct gkyl_range *nrange, const struct gkyl_range *update_range, int num_comp, 
  struct gkyl_array *nodal_fld, const struct gkyl_array *modal_fld);

void gkyl_nodal_ops_m2n_deflated_cu(const struct gkyl_nodal_ops *nodal_ops, 
  const struct gkyl_basis *deflated_cbasis, const struct gkyl_rect_grid *deflated_grid, 
  const struct gkyl_range *nrange, const struct gkyl_range *deflated_nrange, const struct gkyl_range *deflated_update_range, int num_comp, 
  struct gkyl_array *nodal_fld, const struct gkyl_array *deflated_modal_fld, int extra_idx) ;



