#include <gkyl_deflate_geo.h>
#include <gkyl_rect_grid.h>
#include <assert.h>

#include <gkyl_deflate_geo_kernels.h>
#include <gkyl_deflate_geo_surf_kernels.h>

typedef void (*deflate_geo_kernel)(const double *fld, double* deflated_fld);


typedef struct { deflate_geo_kernel kernels[3]; } deflate_geo_kernel_list;
typedef struct { deflate_geo_kernel_list list[4]; } deflate_geo_kernel_remy_list;
typedef struct { deflate_geo_kernel_remy_list list[2]; } deflate_geo_kernel_remx_list;

typedef struct { deflate_geo_kernel_list list[3]; } deflate_geo_surf_kernel_dim_list;


static const deflate_geo_surf_kernel_dim_list ser_deflate_surf_geo_kernel_list[] = {
  { .list = {
      {NULL, NULL, NULL},
      {NULL, NULL, NULL },
      {NULL, NULL, NULL},
    }
  },
  { .list = {
      {NULL, deflate_geo_surfx_1x_ser_p1, NULL},
      {NULL, NULL, NULL },
      {NULL, NULL, NULL},
    }
  },
  { .list = {
      {NULL, deflate_geo_surfx_2x_ser_p1, NULL},
      {NULL, deflate_geo_surfy_2x_ser_p1, NULL },
      {NULL, NULL, NULL},
    }
  },
};

GKYL_CU_DH
static const deflate_geo_kernel_remx_list ser_deflate_geo_kernel_rem_list[] = {
// don't remove x
{ .list = {
    // don't remove y 
    { .list = {
        {NULL, NULL, NULL},
        {NULL, NULL, NULL},
        {NULL, NULL, NULL},
        {NULL, NULL, NULL}
      }
    },

    // remove y 
    { .list = {
      {NULL, NULL, NULL},
      {NULL, NULL, NULL},
      {NULL, deflate_geo_2x_ser_p1_remy, deflate_geo_2x_ser_p2_remy},
      {NULL, NULL, NULL}
      }
    },


    }
  },
// remove x
{ .list = {
    // don't remove y 
    { .list = {
        {NULL, NULL, NULL},
        {NULL, NULL, NULL},
        {NULL, deflate_geo_2x_ser_p1_remx, deflate_geo_2x_ser_p2_remx},
        {NULL, NULL, NULL}
      }
    },

    // remove y 
    { .list = {
      {NULL, NULL, NULL},
      {NULL, deflate_geo_1x_ser_p1_remxy, deflate_geo_1x_ser_p2_remxy},
      {NULL, NULL, NULL},
      {NULL, NULL, NULL}
      }
    },


    }
  },


};


struct gkyl_deflate_geo{
  const struct gkyl_basis *basis;
  const struct gkyl_basis *deflated_basis;
  const struct gkyl_rect_grid* grid;
  const struct gkyl_rect_grid* deflated_grid;
  deflate_geo_kernel kernel;
  int *rem_dirs;
  bool use_gpu;
};

struct gkyl_deflate_geo_surf{
  const struct gkyl_basis *basis;
  int deflated_num_basis;
  const struct gkyl_rect_grid* grid;
  const struct gkyl_rect_grid* deflated_grid;
  deflate_geo_kernel kernel;
  int *rem_dirs;
  int dir;
  bool use_gpu;
};

GKYL_CU_DH
static deflate_geo_kernel
deflate_geo_choose_kernel(const int *rem_dirs, const int cdim, const int basis_type, const int poly_order)
{
  switch (basis_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_deflate_geo_kernel_rem_list[rem_dirs[0]].list[rem_dirs[1]].list[cdim].kernels[poly_order];

    default:
      assert(false);
      break;
  }
}

GKYL_CU_DH
static deflate_geo_kernel
deflate_geo_surf_choose_kernel(const int dir, const int cdim, const int basis_type, const int poly_order)
{
  switch (basis_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_deflate_surf_geo_kernel_list[cdim].list[dir].kernels[poly_order];

    default:
      assert(false);
      break;
  }
}






