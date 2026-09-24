#pragma once

#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_array.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_skin_surf_from_ghost gkyl_skin_surf_from_ghost;

/**
 * Create an updater that sets the values of the skin cell at the skin-ghost boundary to equal 
 * the ghost cell evaluated on that boundary, while maintaining the other skin nodal values.
 * The updater alters the DG coefficient of the skin cell such that
 *  ```math
 *  f_{new}(x_s) = f_g(x_s) if x_s \in \partial\Omega_{sg}
 *  ```
 * or
 *  ```math
 *  f_{new}(x_s) = f(x_s) if x_s \in \partial \Omega_s  \setminus \partial\Omega_{sg}
 *  ```
 * where $f$ is the manipulated field, $f_g$ is the ghost field, $\partial\Omega_{sg}$ is the 
 * skin cell boundary, and $\partial\Omega_{sg}$ is the skin-ghost boundary.
 * Note that the cell average may be altered through this operation.
 *
 * @param dir Direction along which to apply the skin-ghost boundary copy.
 * @param edge Lower or upper edge where the boundary condition is applied (see gkyl_edge_loc).
 * @param basis Basis on which the array coefficients are expanded.
 * @param skin_r Range representing the skin (boundary) region.
 * @param ghost_r Range representing the ghost (outer boundary) region.
 * @param use_gpu Boolean flag to indicate whether GPU computation should be used.
 * @return Pointer to the newly created updater.
 */
struct gkyl_skin_surf_from_ghost* gkyl_skin_surf_from_ghost_new(int dir, enum gkyl_edge_loc edge,
  const struct gkyl_basis basis, const struct gkyl_range *skin_r, const struct gkyl_range *ghost_r,
  bool use_gpu);

/**
 * Enforce that the value of the skin cell at the node facing the ghost cell is equal to the ghost value
 * at that same node. The other nodal value of the skin cell remains unchanged (see new routine description above).
 *
 * @param up Pointer to the boundary condition updater.
 * @param field Array representing the field values to update (currently works only in configuration space).
 */
void gkyl_skin_surf_from_ghost_advance(const struct gkyl_skin_surf_from_ghost *up, struct gkyl_array *field);

/**
 * Free the memory associated with the skin_surf_from_ghost updater.
 *
 * @param up Pointer to the boundary condition updater to be released.
 */
void gkyl_skin_surf_from_ghost_release(struct gkyl_skin_surf_from_ghost *up);