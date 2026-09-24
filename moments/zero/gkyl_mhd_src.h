#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_wv_mhd.h>

struct gkyl_mhd_src_inp {
  const struct gkyl_rect_grid *grid; // grid on which to solve equations
  enum gkyl_wv_mhd_div_constraint divergence_constraint;
  double glm_ch;
  double glm_alpha; // Mignone & Tzeferacos, JCP (2010) 229, 2117, Eq (27).
  double dxyz_min;
};

// Object type
struct gkyl_mhd_src {
  struct gkyl_rect_grid grid;
  int ndim;
  enum gkyl_wv_mhd_div_constraint divergence_constraint;
  double glm_ch;
  double glm_alpha; // Mignone & Tzeferacos, JCP (2010) 229, 2117, Eq (27).
  double dxyz_min;

  struct gkyl_array *divB_array;
  struct gkyl_array *B_dot_gradPsi_array;
};
typedef struct gkyl_mhd_src gkyl_mhd_src;

/**
 * Create new updater.
 *
 * @param inp Input parameters to updater
 */
gkyl_mhd_src *gkyl_mhd_src_new(struct gkyl_mhd_src_inp inp,
  const struct gkyl_range *local_ext);

/**
 * @param mes
 * @param dt
 */
void gkyl_mhd_src_advance(const gkyl_mhd_src *mes, double dt,
  const struct gkyl_range *update_rng,
  struct gkyl_array *q,
  const struct gkyl_array *app_accel);

/**
 * Delete updater.
 *
 * @param mes Updater to delete.
 */
void gkyl_mhd_src_release(gkyl_mhd_src *mes);

/**
 * Set ch parameter used by the GLM divergence constraint.
 *
 * @param up MHD source updater
 * @param ch number
 */
void gkyl_mhd_src_set_glm_ch(struct gkyl_mhd_src* up, double glm_ch);

/**
 * Get average |div(B)| for diagnostic usage.
 *
 * @param up MHD source updater
 * @param update_range
 * @param q_array input field
 */
double gkyl_mhd_src_calc_divB(const gkyl_mhd_src *up,
  const struct gkyl_range *update_range,
  struct gkyl_array *q_array);
