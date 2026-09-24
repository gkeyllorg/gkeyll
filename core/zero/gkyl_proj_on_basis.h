#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_evalf_def.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_proj_on_basis gkyl_proj_on_basis;

// Types for various kernels.
typedef void (*proj_on_basis_c2p_t)(const double *xcomp, double *xphys, void *ctx);

// input packaged as a struct
struct gkyl_proj_on_basis_inp {
  const struct gkyl_rect_grid *grid; // grid on which to project
  const struct gkyl_basis *basis; // basis functions

  enum gkyl_quad_type qtype; // quadrature to use
  
  int num_quad; // number of quadrature points
  int num_ret_vals; // number of return values in eval function
  evalf_t eval; // function to project
  void *ctx; // function context

  proj_on_basis_c2p_t c2p_func; // Function that transforms a set of ndim
                                // computational coordinates to physical ones.
  void *c2p_func_ctx; // Context for c2p_func.
};

/**
 * Create new updater to project function on basis functions on a
 * grid. Free using gkyl_proj_on_basis_release method.
 *
 * @param inp Input parameters
 * @return New updater pointer.
 */
gkyl_proj_on_basis* gkyl_proj_on_basis_inew(const struct gkyl_proj_on_basis_inp *inp);

/**
 * Create new updater to project function on basis functions on a
 * grid, passing inputs separately. Free using gkyl_proj_on_basis_release method.
 *
 * @param grid Grid object
 * @param basis Basis functions to project on
 * @param num_quad Number of quadrature nodes
 * @param num_ret_vals Number of values 'eval' sets
 * @param eval Function to project.
 * @param ctx Context for function evaluation. Can be NULL.
 * @return New updater pointer.
 */
gkyl_proj_on_basis *gkyl_proj_on_basis_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *basis, int num_quad, int num_ret_vals, evalf_t eval, void *ctx);

/**
 * Compute projection on basis. The update_rng MUST be a sub-range of
 * the range on which the array is defined. That is, it must be either
 * the same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param pob Project on basis updater to run
 * @param tm Time at which projection must be computed
 * @param update_rng Range on which to run projection.
 * @param out Output array
 */
void gkyl_proj_on_basis_advance(const gkyl_proj_on_basis *pob,
  double tm, const struct gkyl_range *update_rng, struct gkyl_array *out);

/**
 * Perform the quadrature in the proj_on_basis procedure.
 * Intended for systems that can't perform the whole procedure in _advance.
 *
 * @param up Project on basis updater.
 * @param fun_at_ords Function evaluated at ordinates in one cell.
 * @param f Output projected function in one cell.
 */
void gkyl_proj_on_basis_quad(const gkyl_proj_on_basis *up, const struct gkyl_array *fun_at_ords, double* f);

/**
 * Return the total number of quadrature points/ordinates.
 *
 * @param up Project on basis updater.
 * @return Number of ordinates.
 */
int gkyl_proj_on_basis_get_tot_quad(const gkyl_proj_on_basis *up);

/**
 * Get the coordinates of a given ordinate.
 *
 * @param up Project on basis updater.
 * @param node Index indicate the desired node.
 * @return Node coordinates.
 */
double* gkyl_proj_on_basis_fetch_ordinate(const gkyl_proj_on_basis *up, long node);

/**
 * Delete updater.
 *
 * @param pob Updater to delete.
 */
void gkyl_proj_on_basis_release(gkyl_proj_on_basis* pob);
