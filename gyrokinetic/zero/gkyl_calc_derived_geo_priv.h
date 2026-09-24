#include <gkyl_calc_derived_geo.h>
#include <gkyl_rect_grid.h>
#include <gkyl_calc_derived_geo_kernels.h>
#include <assert.h>

typedef void (*derived_geo_kernel)(const double *gij, const double *bmag, double *J, double *Jinv, double *grij, double *bi, double *cmag, double *Jtot, double *Jtotinv, double *gxxJ, double *gxyJ, double *gyyJ, double *gxzJ, double *eps2);

typedef struct { derived_geo_kernel kernels[3]; } derived_geo_kernel_list;  // For use in kernel tables.
typedef struct { derived_geo_kernel_list list[4]; } derived_geo_node_list;  // For use in kernel tables.

GKYL_CU_DH
static const derived_geo_node_list ser_derived_geo_kernel_list[] = {
  { .list =  {
    { NULL, NULL, NULL }, // 0x No 0D basis functions
    { NULL, NULL, NULL}, // 1x Not tested yet
    { NULL, NULL, NULL}, // 2x Not tested yet
    { NULL, derived_geo_3x_Ser_p1, derived_geo_3x_Ser_p2}
  }
  },
  { .list =  {
    { NULL, NULL, NULL }, // 0x No 0D basis functions
    { NULL, NULL, NULL}, // 1x Not tested yet
    { NULL, NULL, NULL}, // 2x Not tested yet
    { NULL, derived_geo_quad_3x_Ser_p1, NULL}
  }
  },

};

struct gkyl_calc_derived_geo{
  unsigned cdim; // Configuration-space dimension.
  unsigned cnum_basis; // Number of conf-space basis functions.
  unsigned poly_order; // Polynomial order of the basis.
  struct gkyl_basis cbasis; // configuration space basis
  const struct gkyl_rect_grid* grid;
  bool use_gpu;
  derived_geo_kernel kernel;
};

GKYL_CU_DH
static derived_geo_kernel
derived_geo_choose_kernel(int dim, int basis_type, int node_type, int poly_order)
{
  switch (basis_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_derived_geo_kernel_list[node_type].list[dim].kernels[poly_order];
    default:
      assert(false);
      break;
  }
  return 0;
}




