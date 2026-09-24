#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_basis.h>
#include <gkyl_deflate_zsurf.h>
#include <gkyl_nodal_ops.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_fem_poisson.h>


// Object type
typedef struct gkyl_deflated_fem_poisson gkyl_deflated_fem_poisson;

/**
 * Create new updater to solve a poisson equation
 *   - nabla . (epsilon * nabla phi) - kSq * phi = rho
 * on lines or planes using FEM to ensure phi is continuous.
 * The inputs are a DG field rho, epsilon and kSq, and the
 * output is a DG field phi.
 *
 * rho, epsilon and kSq are evaluated at surfaces that are
 * constant in the last coordinate and then the FEM poisson
 * solver is called to solve the poisson equation at each
 * surface. The surfaces are lines (1d) if the initial
 * fields are 2d and planes (2d) if the initial fields are
 * 3d). Finally, the surfaces are assembled and converted
 * into a solution with the initial dimenionality.
 * Free using gkyl_deflated_fem_poisson_release method.
 *
 * @param grid Grid object
 * @param basis_on_dev Basis functions of the DG field stored on the gpu.
 * @param basis Basis functions of the DG field.
 * @param local Range on which the poisson problem should be solved.
 * @param global_sub_range Local range as a sub-range of the global range.
 * @param bcs Boundary conditions.
 * @param epsilon Permittivity tensor. Defined over the extended range.
 * @param kSq Field multiplying phi in Helmholtz equation. Defined over the extended range.
 * @param use_gpu Boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
struct gkyl_deflated_fem_poisson* gkyl_deflated_fem_poisson_new(struct gkyl_rect_grid grid, 
  struct gkyl_basis *basis_on_dev, struct gkyl_basis basis, struct gkyl_range local, 
  struct gkyl_range global_sub_range, struct gkyl_array *epsilon, struct gkyl_array *kSq,
  struct gkyl_poisson_bc poisson_bc, struct gkyl_poisson_bias_plane_list *bias_plane_list, bool use_gpu);

/**
 * Solve the poisson equation for the given charge density
 *
 * @param up Deflated FEM poisson updater to run.
 * @param rhs DG field to set as RHS source (charge density rho).
 * @param phibc Spatially varying BC as a DG (volume) field, defined in the whole domain.
 * @param phi DG field solution to poison problem (phi).
 */
void gkyl_deflated_fem_poisson_advance(struct gkyl_deflated_fem_poisson* up, struct gkyl_array *rhs,
  struct gkyl_array *phibc, struct gkyl_array* phi);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_deflated_fem_poisson_release(struct gkyl_deflated_fem_poisson* up);
