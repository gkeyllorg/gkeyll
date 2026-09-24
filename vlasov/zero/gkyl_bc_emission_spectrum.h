#pragma once

#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_array.h>
#include <gkyl_rect_grid.h>
#include <gkyl_emission_spectrum_model.h>
#include <gkyl_emission_yield_model.h>
#include <gkyl_emission_elastic_model.h>

// Object type
typedef struct gkyl_bc_emission_spectrum gkyl_bc_emission_spectrum;

/**
 * Create a new updater to apply emitting wall spectrum boundary conditions.
 *
 * @param spectrum_model Spectrum model type
 * @param yield_model Yield model type
 * @param yield Array of calculated yield values at cell centers
 * @param spectrum Emission spectrum projected onto basis
 * @param dir Direction in which to apply BC
 * @param edge Lower or upper edge at which to apply BC (emission_spectrum gkyl_edge_loc)
 * @param cdim Configuration space dimensions
 * @param vdim Velocity space dimensions
 * @param mass_in Impacting species mass
 * @param mass_out Emitted species mass
 * @param impact_buff_r Range over the impacting species buffer array
 * @param emit_buff_r Range over the emitting species buffer array
 * @param impact_grid Impacting species boundary grid
 * @param emit_grid Emitted species boundary grid
 * @param poly_order Polynomial order of basis functions.
 * @param basis Basis functions
 * @param proj_buffer Host array to temporarily store projection of emission spectrum
 * @param use_gpu Boolean to indicate whether to use the GPU
 * @return New updater pointer
 */
struct gkyl_bc_emission_spectrum*
gkyl_bc_emission_spectrum_new(struct gkyl_emission_spectrum_model *spectrum_model,
  struct gkyl_emission_yield_model *yield_model, struct gkyl_array *yield,
  struct gkyl_array *spectrum, int dir, enum gkyl_edge_loc edge, int cdim, int vdim,
  double mass_in, double mass_out, struct gkyl_range *impact_buff_r, struct gkyl_range *emit_buff_r,
  struct gkyl_rect_grid *impact_grid, struct gkyl_rect_grid *emit_grid, int poly_order,
  struct gkyl_basis *basis, struct gkyl_array *proj_buffer, bool use_gpu);

/**
 * @param up BC updater
 * @param impact_buff_r Range over the impacting species buffer array
 * @param impact_cbuff_r Configuration space range over the impacting species buffer array
 * @param emit_buff_r Range over the emitting species buffer array
 * @param bflux Boundary flux df/dt in ghost cells
 * @param f_emit Emitted distribution
 * @param yield Array of calculated yield values at cell centers
 * @param spectrum Emission spectrum projected onto basis
 * @param weight Weighting coefficients
 * @param flux Flux into boundary
 * @param k Normalization factor
 */
void
gkyl_bc_emission_spectrum_advance(const struct gkyl_bc_emission_spectrum *up,
  struct gkyl_range *impact_buff_r, struct gkyl_range *impact_cbuff_r,
  struct gkyl_range *emit_buff_r, struct gkyl_array *bflux, struct gkyl_array *f_emit,
  struct gkyl_array *yield, struct gkyl_array *spectrum, struct gkyl_array *weight,
  struct gkyl_array *flux, struct gkyl_array *k);

/**
 * Loop over impacting species velocity space and calculate SEY at cell centers
 *
 * @param up BC updater
 * @param yield Array of calculated yield values at cell centers
 * @param grid Impacting species boundary grid
 * @param impact_buff_r Range over the impacting species buffer array
 */
void
gkyl_bc_emission_spectrum_sey_calc(const struct gkyl_bc_emission_spectrum *up,
  struct gkyl_array *yield, struct gkyl_rect_grid *grid, const struct gkyl_range *impact_buffer_r);

/**
 * Create range over velocity space into wall.
 *
 * @param flux_r Output range over impacting velocity space.
 * @param dir Direction in which to apply BC.
 * @param parent Input range over all of velocity space.
 * @param nghost Number of ghost cells.
 * @param edge Lower or upper edge at which to apply BC (emission_spectrum gkyl_edge_loc).
 */
void
gkyl_bc_emission_flux_ranges(struct gkyl_range *impact_buff_r, int dir,
  const struct gkyl_range *parent, const int *nghost, enum gkyl_edge_loc edge);

/**
 * Free memory associated with bc_emission_spectrum updater.
 *
 * @param up BC updater.
 */
void
gkyl_bc_emission_spectrum_release(struct gkyl_bc_emission_spectrum *up);
