#pragma once

// Private header for bc_emission_spectrum updater, not for direct use in user code.

#include <gkyl_bc_emission_spectrum.h>
#include <gkyl_array_ops.h>
#include <math.h>
#include <assert.h>
#include <gkyl_alloc.h>

// Primary struct for the updater
struct gkyl_bc_emission_spectrum {
  int dir, cdim, vdim;
  enum gkyl_edge_loc edge;
  double charge;
  double mass;
  struct gkyl_rect_grid *grid;
  struct gkyl_emission_spectrum_model *spectrum_model;
  struct gkyl_emission_yield_model *yield_model;
  bool use_gpu;
};

// Function to calculate the weighted mean of the SE yield
GKYL_CU_D
static void
bc_weighted_delta(const double *inp, int cdim, int dir, enum gkyl_edge_loc edge, double xc[GKYL_MAX_DIM], const double *gain, double *weight)
{
  if ((edge == GKYL_LOWER_EDGE && xc[cdim+dir] < 0) || (edge == GKYL_UPPER_EDGE && xc[cdim+dir] > 0)) {
    weight[0] += inp[0]*gain[0];
    weight[1] += inp[0];
  }
}

#ifdef GKYL_HAVE_CUDA

void
gkyl_bc_emission_spectrum_set_extern_params_cu(const struct gkyl_bc_emission_spectrum *up,
  int cdim, int vdim, double mass_in, double mass_out);

/**
 * CUDA device function to set up function to apply boundary conditions.
 *
 * @param up BC updater
 * @param f_skin Skin cell distribution
 * @param f_proj Projected spectrum distribution
 * @param f_buff Distribution buffer array
 * @param weight Weighting coefficients
 * @param k Normalization factor
 * @param flux Flux into boundary
 * @param grid Domain grid
 * @param gamma SE yield values on incoming ghost space
 * @param skin_r Incoming skin space range
 * @param ghost_r Incoming ghost space range
 * @param conf_r Configuration space range
 * @param buff_r Buffer array range
 */
void
gkyl_bc_emission_spectrum_advance_cu(const struct gkyl_bc_emission_spectrum *up,
  struct gkyl_range *impact_buff_r, struct gkyl_range *impact_cbuff_r,
  struct gkyl_range *emit_buff_r, struct gkyl_array *bflux, struct gkyl_array *f_emit,
  struct gkyl_array *yield, struct gkyl_array *spectrum, struct gkyl_array *weight,
  struct gkyl_array *flux, struct gkyl_array *k);

/**
 * CUDA device function to set up function to calculate SEY
 *
 * @param up BC updater
 * @param grid Domain grid
 * @param gamma SE yield values on incoming ghost space
 * @param ghost_r Incoming ghost space range
 */
void
gkyl_bc_emission_spectrum_sey_calc_cu(const struct gkyl_bc_emission_spectrum *up,
  struct gkyl_array *yield, struct gkyl_rect_grid *grid, const struct gkyl_range *gamma_r);

#endif
