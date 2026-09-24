#pragma once

#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_array.h>
#include <gkyl_gk_bc_type.h>

// Object type
typedef struct gkyl_bc_basic_gyrokinetic gkyl_bc_basic_gyrokinetic;

/**
 * Create new updater to apply basic BCs to a field
 * in a gkyl_array. Basic BCs are those in which the
 * ghost cell depends solely on the skin cell next to it
 * via a function of type array_copy_func_t (e.g. absorb, reflect).
 *
 * @param dir Direction in which to apply BC.
 * @param edge Lower or upper edge at which to apply BC (see gkyl_edge_loc).
 * @param bctype BC type (see gkyl_bc_basic_gyrokinetic_type).
 * @param basis Basis on which coefficients in array are expanded.
 * @param skin_r Skin range.
 * @param ghost_r Ghost range.
 * @param cdim Configuration space dimensions.
 * @param num_comp Number of components (DOFs) within a cell.
 * @param use_gpu Boolean to indicate whether to use the GPU.
 * @return New updater pointer.
 */
struct gkyl_bc_basic_gyrokinetic* gkyl_bc_basic_gyrokinetic_new(int dir, enum gkyl_edge_loc edge,
  enum gkyl_gyrokinetic_bc_type bctype, const struct gkyl_basis *basis, const struct gkyl_range *skin_r,
  const struct gkyl_range *ghost_r, int num_comp, int cdim, bool use_gpu);

/**
 * Advance boundary conditions *in special case where buffer is fixed in time*. 
 * *Only* fills the buffer based on the input function and array f_arr. 
 * Should be called only when bctype = GKYL_BC_FIXED_FUNC
 *
 * @param up BC updater.
 * @param buff_arr Buffer array, big enough for ghost cells at this boundary.
 * @param f_arr Field array to apply BC to.
 */
void gkyl_bc_basic_gyrokinetic_buffer_fixed_func(const struct gkyl_bc_basic_gyrokinetic *up,
  struct gkyl_array *buff_arr, struct gkyl_array *f_arr);

/**
 * Advance boundary conditions. Fill buffer array based on boundary conditions and copy
 * contents to ghost cells of input f_arr
 *
 * @param up BC updater.
 * @param buff_arr Buffer array, big enough for ghost cells at this boundary.
 * @param f_arr Field array to apply BC to.
 */
void gkyl_bc_basic_gyrokinetic_advance(const struct gkyl_bc_basic_gyrokinetic *up,
  struct gkyl_array *buff_arr, struct gkyl_array *f_arr);

/**
 * Free memory associated with bc_basic_gyrokinetic updater.
 *
 * @param up BC updater.
 */
void gkyl_bc_basic_gyrokinetic_release(struct gkyl_bc_basic_gyrokinetic *up);
