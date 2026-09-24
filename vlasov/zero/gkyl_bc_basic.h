#pragma once

#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_array.h>

// BC types in this updater.
enum gkyl_bc_basic_type { 
  GKYL_BC_COPY = 0, 
  GKYL_BC_ABSORB, 
  GKYL_BC_REFLECT, 
  GKYL_BC_DISTF_REFLECT, 
  GKYL_BC_MAXWELL_PEC, 
  GKYL_BC_MAXWELL_SYM, 
  GKYL_BC_MAXWELL_RESERVOIR, 
  GKYL_BC_FIXED_FUNC,
  GKYL_BC_PKPM_SPECIES_REFLECT,
  GKYL_BC_PKPM_MOM_REFLECT, 
  GKYL_BC_PKPM_MOM_NO_SLIP,
  GKYL_BC_EULER_REFLECT, 
  GKYL_BC_EULER_NO_SLIP,
  GKYL_BC_CONF_BOUNDARY_VALUE, 
};

// Object type
typedef struct gkyl_bc_basic gkyl_bc_basic;

/**
 * Create new updater to apply basic BCs to a field
 * in a gkyl_array. Basic BCs are those in which the
 * ghost cell depends solely on the skin cell next to it
 * via a function of type array_copy_func_t (e.g. absorb, reflect).
 *
 * @param dir Direction in which to apply BC.
 * @param edge Lower or upper edge at which to apply BC (see gkyl_edge_loc).
 * @param bctype BC type (see gkyl_bc_basic_type).
 * @param basis Basis on which coefficients in array are expanded.
 * @param skin_r Skin range.
 * @param ghost_r Ghost range.
 * @param cdim Configuration space dimensions.
 * @param num_comp Number of components (DOFs) within a cell.
 * @param use_gpu Boolean to indicate whether to use the GPU.
 * @return New updater pointer.
 */
struct gkyl_bc_basic* gkyl_bc_basic_new(int dir, enum gkyl_edge_loc edge, enum gkyl_bc_basic_type bctype,
  const struct gkyl_basis *basis, const struct gkyl_range *skin_r,
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
void gkyl_bc_basic_buffer_fixed_func(const struct gkyl_bc_basic *up, struct gkyl_array *buff_arr, struct gkyl_array *f_arr);

/**
 * Advance boundary conditions. Fill buffer array based on boundary conditions and copy
 * contents to ghost cells of input f_arr
 *
 * @param up BC updater.
 * @param buff_arr Buffer array, big enough for ghost cells at this boundary.
 * @param f_arr Field array to apply BC to.
 */
void gkyl_bc_basic_advance(const struct gkyl_bc_basic *up, struct gkyl_array *buff_arr, struct gkyl_array *f_arr);

/**
 * Free memory associated with bc_basic updater.
 *
 * @param up BC updater.
 */
void gkyl_bc_basic_release(struct gkyl_bc_basic *up);
