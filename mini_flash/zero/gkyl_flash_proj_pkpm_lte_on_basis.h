#pragma once

#include <gkyl_basis.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

typedef struct gkyl_flash_proj_pkpm_lte_on_basis gkyl_flash_proj_pkpm_lte_on_basis;

struct gkyl_flash_proj_pkpm_lte_on_basis_inp {
  const struct gkyl_rect_grid *phase_grid; // Phase-space grid on which to compute moments
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions
  const struct gkyl_basis *phase_basis; // Phase-space basis functions
  bool use_gpu; // bool for gpu useage (not supported)

  enum gkyl_quad_type quad_type;
};

/**
 * Initialize an updater to project a Maxwellian in the moving reference frame
 * on basis with Gaussian quadrature.
 * 
 * @param inp 
 * @return struct gkyl_flash_proj_pkpm_lte_on_basis* 
 */
struct gkyl_flash_proj_pkpm_lte_on_basis*
gkyl_flash_proj_pkpm_lte_on_basis_new(struct gkyl_flash_proj_pkpm_lte_on_basis_inp *inp);

/**
 * Project a Maxwellian in the moving reference frame from input moments.
 * The input `moms` array stores [n, T/m].
 * 
 * @param up 
 * @param phase_range 
 * @param conf_range 
 * @param moms 
 * @param f_lte 
 */
void
gkyl_flash_proj_pkpm_lte_on_basis_advance(gkyl_flash_proj_pkpm_lte_on_basis *up,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *moms, struct gkyl_array *f_lte);

/**
 * Release memory allocated for updater. 
 * 
 * @param up 
*/
void
gkyl_flash_proj_pkpm_lte_on_basis_release(gkyl_flash_proj_pkpm_lte_on_basis *up);