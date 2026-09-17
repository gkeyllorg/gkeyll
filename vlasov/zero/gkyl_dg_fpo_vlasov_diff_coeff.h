#pragma once

#include <gkyl_array.h>
#include <gkyl_fpo_vlasov_coeff_recovery.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_basis.h>

/**
 * Compute the diffusion tensor D_ij for the diffusion term of the FPO
 * from an input potential G.
 * 
 * @param grid Grid (for getting cell spacing)
 * @param cbasis Basis functions in configuration space
 * @param pbasis Basis functions in phase space
 * @param range Range to calculate gradient
 * @param fpo_g Input potential
 * @param fpo_diff_coeff Output array of diffusion coefficient
 * @param fpo_diff_coeff_surf Output array of surface expansion of recovered diffusion coefficient at lower cell boundary
*/
void gkyl_calc_fpo_diff_coeff_recovery(const struct gkyl_fpo_vlasov_coeff_recovery* coeff_recovery,
    const struct gkyl_rect_grid *grid, 
    struct gkyl_basis pbasis, const struct gkyl_range *range, const struct gkyl_range *conf_range,
    const struct gkyl_array *gamma, 
    const struct gkyl_array *fpo_g, const struct gkyl_array *fpo_g_surf,
    const struct gkyl_array *fpo_dgdv_surf, const struct gkyl_array *fpo_d2gdv2_surf,
    struct gkyl_array *fpo_diff_coeff, struct gkyl_array *fpo_diff_coeff_surf, bool use_gpu);

void gkyl_calc_fpo_diff_coeff_recovery_cu(const struct gkyl_fpo_vlasov_coeff_recovery* coeff_recovery,
    const struct gkyl_rect_grid *grid, 
    struct gkyl_basis pbasis, const struct gkyl_range *range, const struct gkyl_range *conf_range,
    const struct gkyl_array *gamma, 
    const struct gkyl_array *fpo_g, const struct gkyl_array *fpo_g_surf,
    const struct gkyl_array *fpo_dgdv_surf, const struct gkyl_array *fpo_d2gdv2_surf,
    struct gkyl_array *fpo_diff_coeff, struct gkyl_array *fpo_diff_coeff_surf);
