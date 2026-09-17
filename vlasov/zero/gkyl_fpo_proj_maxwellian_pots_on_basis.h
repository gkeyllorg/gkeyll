#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_proj_maxwellian_pots_on_basis gkyl_proj_maxwellian_pots_on_basis;

/**
 * Create a new updater to project H and G potentials of a Mxwellian
 * on basis functions. Free after use with gkyl_proj_maxwellian_pots_on_basis_release.
 *
 * @param grid Grid object
 * @param conf_basis Configuration space basis functions
 * @param phase_basis Phase space basis functions
 * @param surf_basis Phase space surface basis functions
 * @param num_quad Number of quadrature nodes
 * @return New updater pointer
*/
struct gkyl_proj_maxwellian_pots_on_basis* 
gkyl_proj_maxwellian_pots_on_basis_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis,
  int num_quad, bool use_gpu);

/**
 * Compute the Rosenbluth potentials from the analytic solution for the equivalent
 * Maxwellian distribution function. Computes both h and g in the volume, as well as
 * a number of needed surface expansions at the edge of velocity space for boundary 
 * conditions of the drag and diffusion coefficients.
 * 
 * @param up Project on basis updater to run
 * @param phase_range Phase space range
 * @param conf_range Configuration space range
 * @param prim_moms Primitive moments (n, u_i, vtsq)
 * @param fpo_h Potential H output array
 * @param fpo_g Potential G output array
 * @param fpo_h_surf Potential H surface expansion output array
 * @param fpo_g_surf Potential G surface expansion output array
 * @param fpo_dhdv_surf First derivative of potential H surface expansion output array
 * @param fpo_dgdv_surf First derivative of potential G surface expansion output array
 * @param fpo_d2gdv2_surf Second derivative of potential G surface expansion output array
*/
void 
gkyl_proj_maxwellian_pots_on_basis_advance(const gkyl_proj_maxwellian_pots_on_basis *up,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array* prim_moms,
  struct gkyl_array *fpo_h, struct gkyl_array *fpo_g,
  struct gkyl_array *fpo_h_surf, struct gkyl_array *fpo_g_surf,
  struct gkyl_array *fpo_dhdv_surf, struct gkyl_array *fpo_dgdv_surf,
  struct gkyl_array *fpo_d2gdv2_surf);

/**
 * Compute the drag coefficient a_i = dh/dv_i and diffusion tensor D_ij = d^2g/dv_i dv_j 
 * from the the analytic solution for the equivalent Maxwellian distribution function. 
 * 
 * @param up Project on basis updater to run
 * @param phase_range Phase space range
 * @param conf_range Configuration space range
 * @param prim_moms Primitive moments (n, u_i, vtsq)
 * @param fpo_dhdv First derivatives of potential H output array
 * @param fpo_d2gdv2 Second derivatives of potential G output array
*/
void 
gkyl_proj_maxwellian_pots_deriv_on_basis_advance(const gkyl_proj_maxwellian_pots_on_basis *up,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *prim_moms,
  struct gkyl_array *fpo_dhdv, struct gkyl_array *fpo_d2gdv2);

/**
 * Delete updater.
 *
 *@param up Updater to delete.
*/
void gkyl_proj_maxwellian_pots_on_basis_release(gkyl_proj_maxwellian_pots_on_basis *up);

struct gkyl_proj_maxwellian_pots_on_basis*
gkyl_proj_maxwellian_pots_on_basis_cu_ho_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis,
  int num_quad);
