#pragma once

#include <gkyl_wv_eqn.h>

// Type of Rieman problem solver to use
enum gkyl_wv_mhd_rp {
  WV_MHD_RP_ROE = 0, // default
  WV_MHD_RP_HLLD,
  WV_MHD_RP_LAX
};

// flags to indicate which divergence constraint scheme to use
enum gkyl_wv_mhd_div_constraint {
  GKYL_MHD_DIVB_NONE,
  GKYL_MHD_DIVB_EIGHT_WAVES,
  GKYL_MHD_DIVB_GLM
};

struct gkyl_wv_mhd_inp {
  double gas_gamma; // gas adiabatic constant
  enum gkyl_wv_mhd_rp rp_type; // Riemann problem solver
  enum gkyl_wv_mhd_div_constraint divergence_constraint; // divB correction
  double glm_ch; // factor to use in GLM scheme
  double glm_alpha; // Mignone & Tzeferacos, JCP (2010) 229, 2117, Eq (27).
};

/**
 * Create a new ideal MHD equation object.
 *
 * @param gas_gamma Gas adiabatic constant
 * @param divb Divergence constraint method
 * @return Pointer to mhd equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_mhd_new(const struct gkyl_wv_mhd_inp *inp);

/**
 * Get gas adiabatic constant.
 *
 * @param wv mhd equation object
 * @return Get gas adiabatic constant
 */
double gkyl_wv_mhd_gas_gamma(const struct gkyl_wv_eqn* wv);

/**
 * Get ch parameter used by the GLM divergence constraint.
 *
 * @param wv mhd equation object
 * @return divergence constraint type
 */
double gkyl_wv_mhd_divergence_constraint(const struct gkyl_wv_eqn* wv);

/**
 * Get ch parameter used by the GLM divergence constraint.
 *
 * @param wv mhd equation object
 * @return ch number
 */
double gkyl_wv_mhd_glm_ch(const struct gkyl_wv_eqn* wv);

/**
 * Get alpha parameter used by the GLM divergence constraint.
 *
 * @param wv mhd equation object
 * @return alpha number
 */
double gkyl_wv_mhd_glm_alpha(const struct gkyl_wv_eqn* wv);

/**
 * Set ch parameter used by the GLM divergence constraint.
 *
 * @param wv mhd equation object
 * @param ch number
 */
void gkyl_wv_mhd_set_glm_ch(struct gkyl_wv_eqn* wv, double glm_ch);

