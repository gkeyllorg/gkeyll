#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_pkpm_lte gkyl_pkpm_lte;

// Input of the updater, packaged as a struct.
struct gkyl_pkpm_lte_inp {
  const struct gkyl_rect_grid *vel_grid; // 1V grid of the normalized velocity w
  const struct gkyl_range *vel_range; // velocity range of the output arrays
  int num_quad; // quadrature points per velocity cell
  double eps; // tolerance of the correction (relative error of the moments)
  int max_iter; // maximum number of correction iterations
};

// Status of a correction.
struct gkyl_pkpm_lte_status {
  bool converged; // whether the moments reached the tolerance
  int num_iter; // number of iterations taken
  double error; // largest relative error of the moments
};

/**
 * Create an updater for a Maxwellian in w at one configuration-space node: its projection on the
 * p=2 velocity basis, its moments and the correction of those moments.
 *
 * @param inp Input parameters.
 * @return New updater pointer.
 */
struct gkyl_pkpm_lte *gkyl_pkpm_lte_inew(const struct gkyl_pkpm_lte_inp *inp);

/**
 * Project the Maxwellian with primitive moments [n, u, vtsq] on the p=2 velocity basis.
 *
 * @param up Updater.
 * @param moms Primitive moments [n, u, vtsq].
 * @param f_lte Output: 3 coefficients per velocity cell.
 */
void gkyl_pkpm_lte_proj(
  const struct gkyl_pkpm_lte *up, const double *moms, struct gkyl_array *f_lte
);

/**
 * Moments [M0, M1, M2] of f over the velocity domain, integrating the p=2 expansion exactly.
 *
 * @param up Updater.
 * @param f_lte Input: 3 coefficients per velocity cell.
 * @param moms Output moments [M0, M1, M2].
 */
void gkyl_pkpm_lte_moments(
  const struct gkyl_pkpm_lte *up, const struct gkyl_array *f_lte, double *moms
);

/**
 * Project a Maxwellian whose discrete moments match the primitive moments [n, u, vtsq], by a
 * fixed-point iteration on the moments of the projected Maxwellian. The result stays an exact
 * Maxwellian at the quadrature points.
 *
 * @param up Updater.
 * @param moms Target primitive moments [n, u, vtsq].
 * @param f_lte Output: 3 coefficients per velocity cell.
 * @return Status of the correction.
 */
struct gkyl_pkpm_lte_status gkyl_pkpm_lte_correct(
  const struct gkyl_pkpm_lte *up, const double *moms, struct gkyl_array *f_lte
);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_pkpm_lte_release(struct gkyl_pkpm_lte *up);
