#pragma once

// Private header for gyrokinetic_pol_density updater, not for direct use in user code.

#include <gkyl_gyrokinetic_pol_density.h>
#include <gkyl_gyrokinetic_pol_density_kernels.h>
#include <gkyl_util.h>
#include <assert.h>

// Function pointer type for sheath reflection kernels.
typedef void (*gkpolden_t)(const double *dx,
  const double *epsilon, const double *phi, double *out);

typedef struct { gkpolden_t kernels[3]; } gk_pol_den_kern_list;  // For use in kernel tables.

// Serendipity  kernels.
GKYL_CU_D
static const gk_pol_den_kern_list gk_pol_density_kern_list_ser[] = {
  { gkyl_gyrokinetic_pol_density_1x_ser_p1, NULL, NULL },
  { gkyl_gyrokinetic_pol_density_2x_ser_p1, NULL, NULL },
  { gkyl_gyrokinetic_pol_density_3x_ser_p1, NULL, NULL },
};

struct gkyl_gyrokinetic_pol_density_kernels {
  gkpolden_t pol_den;  // Kernel that computes the polarization density.
};

// Primary struct in this updater.
struct gkyl_gyrokinetic_pol_density {
  struct gkyl_rect_grid grid;  // Phase-space grid.
  bool use_gpu;
  struct gkyl_gyrokinetic_pol_density_kernels *kernels;
};

#ifdef GKYL_HAVE_CUDA
// Declaration of cuda device functions.

void gk_pol_den_choose_kernel_cu(struct gkyl_gyrokinetic_pol_density_kernels *kernels,
  struct gkyl_basis cbasis);

void gkyl_gyrokinetic_pol_density_advance_cu(gkyl_gyrokinetic_pol_density* up,
  const struct gkyl_range *conf_rng, const struct gkyl_array *GKYL_RESTRICT pol_weight,
  const struct gkyl_array *GKYL_RESTRICT phi, struct gkyl_array *GKYL_RESTRICT npol);
#endif

GKYL_CU_D
static void gk_pol_den_choose_kernel(struct gkyl_gyrokinetic_pol_density_kernels *kernels,
  struct gkyl_basis cbasis, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    gk_pol_den_choose_kernel_cu(kernels, cbasis);
    return;
  }
#endif

  enum gkyl_basis_type basis_type = cbasis.b_type;
  int cdim = cbasis.ndim;
  int poly_order = cbasis.poly_order;

  switch (basis_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kernels->pol_den = gk_pol_density_kern_list_ser[cdim-1].kernels[poly_order-1];
      break;
    default:
      assert(false);
      break;
  }
}
