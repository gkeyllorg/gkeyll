#pragma once

// Private header for skin_surf_from_ghost updater, not for direct use in user code.

#include <gkyl_skin_surf_from_ghost.h>
#include <gkyl_skin_surf_from_ghost_kernels.h>
#include <assert.h>

// Function pointer type for the kernel setting the skin value equal to the ghost value at surface
typedef void (*skin_surf_from_ghost_t)(const double *fghost, double *fskin);

typedef struct { skin_surf_from_ghost_t kernels[2]; } skin_surf_from_ghost_kern_list;  // For use in kernel tables.
typedef struct { skin_surf_from_ghost_kern_list dirlist[3]; } dir_skin_surf_from_ghost_kern_list;
typedef struct { dir_skin_surf_from_ghost_kern_list edgedlist[3]; } edged_skin_surf_from_ghost_kern_list;

// Serendipity  kernels.
GKYL_CU_D
static const edged_skin_surf_from_ghost_kern_list ser_skin_surf_from_ghost_list[] = {
  { 
    .edgedlist=
      { 
        { 
          .dirlist=
            {
              { skin_surf_from_ghost_lowerx_1x_ser_p1, NULL },
              { NULL, NULL },
              { NULL, NULL },
            }
        },
        {
          .dirlist=
            {
              { skin_surf_from_ghost_lowerx_2x_ser_p1, NULL },
              { skin_surf_from_ghost_lowery_2x_ser_p1, NULL },
              { NULL, NULL },
            }
        },
        {
          .dirlist=
            {
              { skin_surf_from_ghost_lowerx_3x_ser_p1, NULL },
              { skin_surf_from_ghost_lowery_3x_ser_p1, NULL },
              { skin_surf_from_ghost_lowerz_3x_ser_p1, NULL },
            }
        }
      }
  },
  { 
    .edgedlist=
      { 
        { 
          .dirlist=
            {
              { skin_surf_from_ghost_upperx_1x_ser_p1, NULL },
              { NULL, NULL },
              { NULL, NULL },
            }
        },
        {
          .dirlist=
            {
              { skin_surf_from_ghost_upperx_2x_ser_p1, NULL },
              { skin_surf_from_ghost_uppery_2x_ser_p1, NULL },
              { NULL, NULL },
            }
        },
        {
          .dirlist=
            {
              { skin_surf_from_ghost_upperx_3x_ser_p1, NULL },
              { skin_surf_from_ghost_uppery_3x_ser_p1, NULL },
              { skin_surf_from_ghost_upperz_3x_ser_p1, NULL },
            }
        }
      }
  }
};
 
struct gkyl_skin_surf_from_ghost_kernels {
  skin_surf_from_ghost_t ghost_to_skin;
};

// Primary struct in this updater.
struct gkyl_skin_surf_from_ghost {
  int dir, cdim;
  enum gkyl_edge_loc edge;
  bool use_gpu;
  struct gkyl_skin_surf_from_ghost_kernels *kernels;  // skin surface from ghost kernel.
  struct gkyl_skin_surf_from_ghost_kernels *kernels_cu;  // device copy.
  const struct gkyl_range *skin_r, *ghost_r;
};

#ifdef GKYL_HAVE_CUDA

void skin_surf_from_ghost_choose_kernel_cu(const struct gkyl_basis basis, enum gkyl_edge_loc edge, 
  int dir, struct gkyl_skin_surf_from_ghost_kernels *kers);

void skin_surf_from_ghost_advance_cu(const struct gkyl_skin_surf_from_ghost *up, struct gkyl_array *field);
#endif

GKYL_CU_D
static void skin_surf_from_ghost_choose_kernel(const struct gkyl_basis basis, enum gkyl_edge_loc edge, 
  int dir, bool use_gpu, struct gkyl_skin_surf_from_ghost_kernels *kernels)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    skin_surf_from_ghost_choose_kernel_cu(basis, edge, dir, kernels);
    return;
  }
#endif

  int dim = basis.ndim;
  enum gkyl_basis_type basis_type = basis.b_type;
  int poly_order = basis.poly_order;
  switch (basis_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kernels->ghost_to_skin = ser_skin_surf_from_ghost_list[edge].edgedlist[dim-1].dirlist[dir].kernels[poly_order-1];
      break;
    default:
      assert(false);
      break;
  }
}
