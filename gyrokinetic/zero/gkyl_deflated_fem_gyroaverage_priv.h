#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_deflate_surf_kernels.h>
#include <gkyl_deflate_zsurf.h>
#include <gkyl_deflated_fem_gyroaverage.h>
#include <gkyl_fem_gyroaverage.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

// Maximum number of basis functions in configuration space and in a plane.
#define DEFL_GYROAVG_MAX_BASIS 8
#define DEFL_GYROAVG_MAX_DEFL_BASIS 4

// Kernel that evaluates a field at the lower or upper z surface of a cell.
typedef void (*defl_gyroavg_deflate_kern_t)(const double *fld, double *fld_deflated);

// Updater type.
struct gkyl_deflated_fem_gyroaverage {
  int cdim; // Number of configuration space dimensions.
  int num_basis; // Number of configuration space basis functions.
  int defl_num_basis; // Number of basis functions in a plane.
  int num_mu; // Number of magnetic moments.
  int num_planes; // Number of planes of constant z (z nodes) in the local range.
  int num_cells_z; // Number of cells along z in the local range.

  struct gkyl_rect_grid deflated_grid; // Grid of a plane.
  struct gkyl_basis basis; // Configuration space basis.
  struct gkyl_basis deflated_basis; // Basis of a plane.
  struct gkyl_range local; // Local configuration space range.
  struct gkyl_range deflated_local, deflated_local_ext; // Ranges of a plane.

  // Evaluate a field at the lower and upper z surfaces of a cell (on the host).
  struct gkyl_deflate_zsurf *deflator_lo, *deflator_up;
  // Cell (z index) and deflator used to evaluate a field on each plane.
  int *plane_zidx;
  struct gkyl_deflate_zsurf **plane_deflator;

  struct gkyl_fem_gyroaverage **gyroavg; // Gyroaverage updater of each plane (on the host).

  // Configuration space coefficients (second index) of the product of each
  // basis function of a plane (third index) and the linear function of z that
  // is 1 at the lower/upper (first index 0/1) z surface of a cell and 0 at the other.
  double inflate[2][DEFL_GYROAVG_MAX_BASIS][DEFL_GYROAVG_MAX_DEFL_BASIS];

  // Buffers used on the host.
  struct gkyl_array *plane_fld; // Field on a plane.
  struct gkyl_array **plane_fld_gyroavg; // Gyroaveraged field on each plane.
  struct gkyl_array *plane_dens_lo, *plane_dens_up, *plane_dens_buff; // Densities on a plane.
  struct gkyl_array *plane_adj_lo, *plane_adj_up; // Adjoint gyroaverages on a plane.

  bool use_gpu; // Whether to run on the GPU.
  // On the GPU all the planes are processed at once. The fields on the planes are stored one
  // after the other in an array (see gkyl_fem_gyroaverage_batch_advance).
  // Kernels evaluating a field at the lower (0) and upper (1) z surfaces of a cell.
  defl_gyroavg_deflate_kern_t deflate[2];
  struct gkyl_fem_gyroaverage_batch *batch; // Gyroaverages of all the planes.
  // Gyroaverages at the lower and upper z surfaces of all the cells along z (for the adjoint):
  // the fields 2*c and 2*c+1 are at the lower and upper z surfaces of the c-th cell along z.
  struct gkyl_fem_gyroaverage_batch *batch_adj;
  struct gkyl_array *planes_fld, *planes_fld_gyroavg; // Field and its gyroaverage on the planes.
  // Densities and their adjoint gyroaverages at the z surfaces of the cells.
  struct gkyl_array *planes_dens, *planes_adj;
  struct gkyl_deflated_fem_gyroaverage *on_dev; // Pointer to itself or to a copy on the device.
};

// Kernels evaluating a field at the lower and upper z surfaces of a cell, for each number of
// configuration space dimensions.
GKYL_CU_D static const defl_gyroavg_deflate_kern_t defl_gyroavg_deflate_kernels[][2] = {
  {NULL, NULL},
  {NULL, NULL},
  {deflate_surfy_lower_2x_ser_p1, deflate_surfy_upper_2x_ser_p1},
  {deflate_surfz_lower_3x_ser_p1, deflate_surfz_upper_3x_ser_p1}
};

GKYL_CU_D static defl_gyroavg_deflate_kern_t
defl_gyroavg_choose_deflate_kernel(int cdim, int edge)
{
  return defl_gyroavg_deflate_kernels[cdim][edge];
}

GKYL_CU_DH static inline void
defl_gyroavg_inflate_cell(
  const struct gkyl_deflated_fem_gyroaverage *up, double wgt_same, double wgt_other,
  const double *plane_lo, const double *plane_up, double *out
)
{
  // Configuration space coefficients in a cell of the field that is linear in z and equals
  // wgt_same*plane_lo + wgt_other*plane_up at the lower z surface of the cell, and
  // wgt_same*plane_up + wgt_other*plane_lo at the upper one.
  for (int k = 0; k < up->num_basis; k++) {
    double coeff = 0.0;
    for (int kd = 0; kd < up->defl_num_basis; kd++) {
      coeff += up->inflate[0][k][kd] * (wgt_same * plane_lo[kd] + wgt_other * plane_up[kd]) +
               up->inflate[1][k][kd] * (wgt_same * plane_up[kd] + wgt_other * plane_lo[kd]);
    }
    out[k] = coeff;
  }
}

#ifdef GKYL_HAVE_CUDA
/**
 * Create a copy of the updater on the device, and store a pointer to it in the updater.
 *
 * @param up Updater (on the host).
 */
void defl_gyroavg_dev_new(struct gkyl_deflated_fem_gyroaverage *up);

/**
 * Gyroaverage a field at each magnetic moment on the GPU. See
 * gkyl_deflated_fem_gyroaverage_advance for documentation.
 */
void gkyl_deflated_fem_gyroaverage_advance_cu(
  struct gkyl_deflated_fem_gyroaverage *up, const struct gkyl_array *fin, struct gkyl_array *fout
);

/**
 * Gyroaverage a field that may be discontinuous along z, cell by cell along z, on the GPU. See
 * gkyl_deflated_fem_gyroaverage_advance_dg for documentation.
 */
void gkyl_deflated_fem_gyroaverage_advance_dg_cu(
  struct gkyl_deflated_fem_gyroaverage *up, const struct gkyl_array *fin, struct gkyl_array *fout
);

/**
 * Apply the adjoint of the gyroaverage on the GPU. See
 * gkyl_deflated_fem_gyroaverage_advance_adjoint for documentation.
 */
void gkyl_deflated_fem_gyroaverage_advance_adjoint_cu(
  struct gkyl_deflated_fem_gyroaverage *up, const struct gkyl_array *fin, struct gkyl_array *fout
);
#endif
