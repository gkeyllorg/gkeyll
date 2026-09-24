#include <gkyl_bc_block_tensor.h>
#include <gkyl_rect_grid.h>
#include <assert.h>
#include <gkyl_basis_ser_2x_p1_surfx1_eval_quad.h>
#include <gkyl_basis_ser_2x_p1_surfx2_eval_quad.h>
#include <gkyl_basis_ser_3x_p1_surfx1_eval_quad.h>
#include <gkyl_basis_ser_3x_p1_surfx2_eval_quad.h>
#include <gkyl_basis_ser_3x_p1_surfx3_eval_quad.h>

typedef double (*modal_to_quad_kernel)(const double *f);


typedef struct { modal_to_quad_kernel kernels[4]; } modal_to_quad_kernel_list;
typedef struct { modal_to_quad_kernel_list list[6]; } modal_to_quad_kernel_dim_list;



GKYL_CU_D
static const modal_to_quad_kernel_dim_list ser_modal_to_quad_kernel_dim_list[] = {
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
  // 2X
  { .list = {
      {ser_2x_p1_surfx1_eval_quad_node_0_l, ser_2x_p1_surfx1_eval_quad_node_1_l, NULL, NULL},
      {ser_2x_p1_surfx1_eval_quad_node_0_r, ser_2x_p1_surfx1_eval_quad_node_1_r, NULL, NULL},
      {ser_2x_p1_surfx2_eval_quad_node_0_l, ser_2x_p1_surfx2_eval_quad_node_1_l, NULL, NULL},
      {ser_2x_p1_surfx2_eval_quad_node_0_r, ser_2x_p1_surfx2_eval_quad_node_1_r, NULL, NULL},
      {NULL, NULL, NULL , NULL},
      {NULL, NULL, NULL , NULL},
    }
  },
  // 3X
  { .list = {
      {ser_3x_p1_surfx1_eval_quad_node_0_l, ser_3x_p1_surfx1_eval_quad_node_1_l, ser_3x_p1_surfx1_eval_quad_node_2_l, ser_3x_p1_surfx1_eval_quad_node_3_l},
      {ser_3x_p1_surfx1_eval_quad_node_0_r, ser_3x_p1_surfx1_eval_quad_node_1_r, ser_3x_p1_surfx1_eval_quad_node_2_r, ser_3x_p1_surfx1_eval_quad_node_3_r},
      {ser_3x_p1_surfx2_eval_quad_node_0_l, ser_3x_p1_surfx2_eval_quad_node_1_l, ser_3x_p1_surfx2_eval_quad_node_2_l, ser_3x_p1_surfx2_eval_quad_node_3_l},
      {ser_3x_p1_surfx2_eval_quad_node_0_r, ser_3x_p1_surfx2_eval_quad_node_1_r, ser_3x_p1_surfx2_eval_quad_node_2_r, ser_3x_p1_surfx2_eval_quad_node_3_r},
      {ser_3x_p1_surfx3_eval_quad_node_0_l, ser_3x_p1_surfx3_eval_quad_node_1_l, ser_3x_p1_surfx3_eval_quad_node_2_l, ser_3x_p1_surfx3_eval_quad_node_3_l},
      {ser_3x_p1_surfx3_eval_quad_node_0_r, ser_3x_p1_surfx3_eval_quad_node_1_r, ser_3x_p1_surfx3_eval_quad_node_2_r, ser_3x_p1_surfx3_eval_quad_node_3_r},
    }
  },
};


struct bc_block_tensor {
  struct gkyl_range range;
  struct gkyl_range range_ext;
  struct gkyl_basis basis;
  struct gkyl_rect_grid grid;
  int cdim;
  int poly_order;
  int num_surf_nodes;
  struct gkyl_array* tensor;
  modal_to_quad_kernel_list kernels_lo;
  modal_to_quad_kernel_list kernels_up;
};


GKYL_CU_D
static modal_to_quad_kernel
bc_block_tensor_choose_kernel(int cdim, int edge, int dir, int node_num)
{
  return ser_modal_to_quad_kernel_dim_list[cdim].list[2*dir + edge].kernels[node_num];
}




