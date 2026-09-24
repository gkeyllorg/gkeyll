// Private header: not for direct use
#pragma once

#include <math.h>
#include <gkyl_deflate_surf_kernels.h>
#include <assert.h>

typedef void (*deflate_zsurf_kernel)(const double *fld, double* deflated_fld);

typedef struct { deflate_zsurf_kernel kernels[3]; } deflate_zsurf_kernel_list;
typedef struct { deflate_zsurf_kernel_list list[2]; } deflate_zsurf_kernel_dim_list;

GKYL_CU_D
static const deflate_zsurf_kernel_dim_list ser_deflate_zsurf_kernel_dim_list[] = {
  { .list = {
      {NULL, NULL, NULL},
      {NULL, NULL, NULL},
    }
  },
  { .list = {
      {NULL, NULL, NULL},
      {NULL, NULL, NULL},
    }
  },
  { .list = {
      {NULL, deflate_surfy_lower_2x_ser_p1, deflate_surfy_lower_2x_ser_p2},
      {NULL, deflate_surfy_upper_2x_ser_p1, deflate_surfy_upper_2x_ser_p2},
    }
  },
  { .list = {
      {NULL, deflate_surfz_lower_3x_ser_p1, deflate_surfz_lower_3x_ser_p2},
      {NULL, deflate_surfz_upper_3x_ser_p1, deflate_surfz_upper_3x_ser_p2},
    }
  },
};


struct gkyl_deflate_zsurf {
  int num_basis; // Number of basis functions in full basis
  int num_deflated_basis; // Number of basis functions in deflated basis
  int cdim; // Dimension of un-deflated grid
  deflate_zsurf_kernel kernel;
  uint32_t flags;
  struct gkyl_deflate_zsurf *on_dev; // pointer to itself or device data
  bool use_gpu;
};

GKYL_CU_D
static deflate_zsurf_kernel
deflate_zsurf_choose_kernel(enum gkyl_basis_type basis_type, int dim, int edge, int poly_order)
{
  switch (basis_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_deflate_zsurf_kernel_dim_list[dim].list[edge].kernels[poly_order];

    default:
      assert(false);
      break;
  }

  return 0;
}

/**
 * Create new updater deflate a 2d (x,z)  or 3d (x,y,z) modal expansion to a 1d (x) modal expansion or 2d (x,y)
 * NV-GPU. See new() method for documentation.
 */
struct gkyl_deflate_zsurf* 
gkyl_deflate_zsurf_cu_dev_new(const struct gkyl_basis *cbasis, 
  const struct gkyl_basis *deflated_cbasis, int edge);
