#include <gkyl_dg_basis_ops.h>

#ifdef GKYL_HAVE_CUDA
/**
 * Evaluate the an array containing a DG field
 * at a specific coordinate in the grid on the NVIDIA gpu. 
 *
 * @param arr Array with DG expansion.
 * @param coord Coordinate to evaluate at.
 * @param basis DG basis object (must be on the GPU if arr is on GPU).
 * @param grid Grid object.
 * @param rng Range object (should contain coord).
 * @param out Evaluation of the array at coord (on GPU if arr is on GPU).
 **/
void
gkyl_dg_basis_ops_eval_array_at_coord_comp_cu(const struct gkyl_array *arr, const double *coord,
  const struct gkyl_basis *basis, const struct gkyl_rect_grid *grid, const struct gkyl_range *rng,
  double *out);
#endif
