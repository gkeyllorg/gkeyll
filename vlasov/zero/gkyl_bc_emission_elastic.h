#pragma once

#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_array.h>
#include <gkyl_rect_grid.h>
#include <gkyl_emission_elastic_model.h>

// Object type
typedef struct gkyl_bc_emission_elastic gkyl_bc_emission_elastic;

/**
 * Create a new updater to apply emitting wall spectrum boundary conditions.
 *
 * @param elastic_model Elastic model type
 * @param elastic_yield Projection of elastic yield model onto basis
 * @param dir Direction in which to apply BC
 * @param edge Lower or upper edge at which to apply BC (emission_spectrum gkyl_edge_loc)
 * @param cdim Configuration space dimensions
 * @param vdim Velocity space dimensions
 * @param mass Mass of species
 * @param ncomp Number of components
 * @param grid Impacting species boundary grid
 * @param emit_buff_r Range over the emitting species buffer array
 * @param poly_order Polynomial order of basis functions
 * @param dev_basis Pointer to basis functions on device
 * @param basis Pointer to basis functions on host
 * @param proj_buffer Host array to temporarily store projection of emission spectrum
 * @param use_gpu Boolean to indicate whether to use the GPU
 * @return New updater pointer
 */
struct gkyl_bc_emission_elastic*
gkyl_bc_emission_elastic_new(struct gkyl_emission_elastic_model *elastic_model,
  struct gkyl_array *elastic_yield, int dir, enum gkyl_edge_loc edge, int cdim,
  int vdim, double mass, int ncomp, struct gkyl_rect_grid *grid, struct gkyl_range *emit_buff_r,
  int poly_order, const struct gkyl_basis *dev_basis, struct gkyl_basis *basis,
  struct gkyl_array *proj_buffer, bool use_gpu);

/**
 * @param up BC updater
 * @param emit_skin_r Range over the species skin cells
 * @param buff_arr BC buffer array
 * @param f_skin Skin cell distribution
 * @param f_emit Emitted distribution
 * @param elastic_yield Projection of elastic yield model onto basis
 * @param basis Pointer to basis functions on host
 */
void
gkyl_bc_emission_elastic_advance(const struct gkyl_bc_emission_elastic *up,
  struct gkyl_range *emit_skin_r, struct gkyl_array *buff_arr, struct gkyl_array *f_skin,
  struct gkyl_array *f_emit, struct gkyl_array *elastic_yield, struct gkyl_basis *basis);

/**
 * @param dir Direction in which to apply BC
 * @param cdim Configuration space dimensions
 * @param basis Pointer to basis functions on device
 * @param ncomp Number of components
 * @param use_gpu Boolean to indicate whether to use the GPU
 */
struct gkyl_array_copy_func*
gkyl_bc_emission_elastic_create_arr_copy_func(int dir, int cdim,
  const struct gkyl_basis *basis, int ncomp, bool use_gpu);

/**
 * Free memory associated with bc_emission_elastic updater.
 *
 * @param up BC updater.
 */
void
gkyl_bc_emission_elastic_release(struct gkyl_bc_emission_elastic *up);
