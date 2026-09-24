#include <gkyl_array.h>
#include <gkyl_rect_grid.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>

typedef struct bc_block_tensor  bc_block_tensor;


/**
 * Create a new updater to compute the transformation required to pass fluxes from another block with different
 * geometry to this block
 *
 * @param grid of this block
 * @param range, range_ext : range and extended range of this block
 * @param basis configuration space basis
 * @param use_gpu whether or not to use a gpu
 */
struct bc_block_tensor*
gkyl_bc_block_tensor_new(const struct gkyl_rect_grid* grid, const struct gkyl_range *range, const struct gkyl_range* range_ext, 
  const struct gkyl_basis* basis, bool use_gpu);


/**
 * Take in modal expansions of duals of one block and tangents of the other (cartesian components)
 * in a single cell
 * 1 denotes the block which fluxes are leaving
 * 2 denotes the block which fluxes are entering
 * and calculate T^j'_i = e^j' \dot e_i at the quadrature nodes of the interface
 * @param up bc_block_tensor object. Stored with i changing fastest, then quad node location, then j'
 * @param edge1 edge of block which fluxes leave (0 is lower, 1 is upper)
 * @param edge2 edge of block which fluxes enter(0 is lower, 1 is upper)
 * @param ej duals of block which fluxes enter
 * @param e_i tangent vectors of block which fluxes leave
 */
void calc_tensor(struct bc_block_tensor *up, int dir, int edge1, int edge2, const double *ej, const double *e_i, double *tj_i);


/**
 * Take in modal expansions of duals of one block and tangents of the other (cartesian components)
 * and calculate T^j'_i = e^j' \dot e_i at the quadrature nodes of the interface
 * @param up bc_block_tensor object. Tensor Stored with i changing fastest, then quad node location, then j'
 * @param edge1 edge of block which fluxes leave (0 is lower, 1 is upper)
 * @param edge2 edge of block which fluxes enter(0 is lower, 1 is upper)
 * @param ej duals of block which fluxes enter
 * @param e_i tangent vectors of block which fluxes leave
 */
void gkyl_bc_block_tensor_advance(struct bc_block_tensor* up, int dir, int edge1, int edge2,
    struct gkyl_array* dxdz1, struct gkyl_array* dzdx2, struct gkyl_range *range1, struct gkyl_range *range2);



/**
 * Free the bc_block_tensor updater
 * @param up updater to be freed
 */
void gkyl_bc_block_tensor_release(struct bc_block_tensor* up);
