#pragma once

#include <gkyl_gk_dens_mu_cells.h>
#include <gkyl_mom_gyrokinetic_kernels.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

// Maximum number of configuration space basis functions, and of basis
// functions of the (configuration space, mu) expansion.
#define GK_DENS_MU_MAX_BASIS 8
#define GK_DENS_MU_MAX_CMU_BASIS 16

// Kernel that integrates the distribution function over vpar.
typedef void (*gk_dens_mu_cells_kern_t)(
  const double *dxv, const double *vmap, double m_, const double *bmag, const double *f,
  double *GKYL_RESTRICT out
);

// Updater type.
struct gkyl_gk_dens_mu_cells {
  int cdim, pdim; // Configuration and phase space dimensions.
  int num_basis; // Number of configuration space basis functions.
  int num_cmu_basis; // Number of basis functions of the (configuration space, mu) expansion.
  struct gkyl_rect_grid phase_grid; // Phase space grid.
  double mass; // Species mass.
  gk_dens_mu_cells_kern_t kernel; // Kernel integrating over vpar.
  // Configuration space coefficients (first index) of the integral over a mu
  // cell of each (configuration space, mu) basis function (second index).
  double mu_weight[GK_DENS_MU_MAX_BASIS][GK_DENS_MU_MAX_CMU_BASIS];
  bool use_gpu; // Whether to run on the GPU.
  struct gkyl_gk_dens_mu_cells *on_dev; // Pointer to itself or to a copy on the device.
};

// Kernels integrating over vpar, for each number of configuration space dimensions.
GKYL_CU_D static const gk_dens_mu_cells_kern_t gk_dens_mu_cells_kernels[] = {
  NULL, gyrokinetic_M0_step1_1x2v_ser_p1, gyrokinetic_M0_step1_2x2v_ser_p1,
  gyrokinetic_M0_step1_3x2v_ser_p1
};

GKYL_CU_D static gk_dens_mu_cells_kern_t
gk_dens_mu_cells_choose_kernel(int cdim)
{
  return gk_dens_mu_cells_kernels[cdim];
}

#ifdef GKYL_HAVE_CUDA
/**
 * Create a copy of the updater on the device, and store a pointer to it in the updater.
 *
 * @param up Updater (on the host).
 */
void gk_dens_mu_cells_dev_new(struct gkyl_gk_dens_mu_cells *up);

/**
 * Compute the mu-resolved density on the GPU. See gkyl_gk_dens_mu_cells_advance for
 * documentation.
 */
void gkyl_gk_dens_mu_cells_advance_cu(
  const struct gkyl_gk_dens_mu_cells *up, const struct gkyl_range *phase_range,
  const struct gkyl_range *conf_range, const struct gkyl_array *fin, struct gkyl_array *out
);
#endif
