// Private header for fem_parproj updater.
#pragma once
#include <gkyl_fem_parproj_kernels.h>
#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_mat_triples.h>
#include <gkyl_superlu_ops.h>
#ifdef GKYL_HAVE_CUDA
#include <gkyl_culinsolver_ops.h>
#endif

static long
gkyl_fem_parproj_global_num_nodes(const struct gkyl_basis *basis, bool isperiodic, int parnum_cells)
{
  int dim = basis->ndim;
  int poly_order = basis->poly_order;
  enum gkyl_basis_type basis_type = basis->b_type;

  if (dim==1) {
    if (poly_order == 1) {
      return isperiodic? fem_parproj_num_nodes_global_1x_ser_p1_periodicx(parnum_cells)
                       : fem_parproj_num_nodes_global_1x_ser_p1_nonperiodicx(parnum_cells);
    } else if (poly_order == 2) {
      return isperiodic? fem_parproj_num_nodes_global_1x_ser_p2_periodicx(parnum_cells)
                       : fem_parproj_num_nodes_global_1x_ser_p2_nonperiodicx(parnum_cells);
    }
  } else if (dim==2) {
    if (basis_type == GKYL_BASIS_MODAL_SERENDIPITY) {
      if (poly_order == 1) {
        return isperiodic? fem_parproj_num_nodes_global_2x_ser_p1_periodicy(parnum_cells)
                         : fem_parproj_num_nodes_global_2x_ser_p1_nonperiodicy(parnum_cells);
      } else if (poly_order == 2) {
        return isperiodic? fem_parproj_num_nodes_global_2x_ser_p2_periodicy(parnum_cells)
                         : fem_parproj_num_nodes_global_2x_ser_p2_nonperiodicy(parnum_cells);
      }
    }
  } else if (dim==3) {
    if (basis_type == GKYL_BASIS_MODAL_SERENDIPITY) {
      if (poly_order == 1) {
        return isperiodic? fem_parproj_num_nodes_global_3x_ser_p1_periodicz(parnum_cells)
                         : fem_parproj_num_nodes_global_3x_ser_p1_nonperiodicz(parnum_cells);
      } else if (poly_order == 2) {
        return isperiodic? fem_parproj_num_nodes_global_3x_ser_p2_periodicz(parnum_cells)
                         : fem_parproj_num_nodes_global_3x_ser_p2_nonperiodicz(parnum_cells);
      }
    }
  }
  assert(false);  // Other dimensionalities or basis not supported.
  return -1;
}

// Function pointer type for local-to-global mapping.
typedef void (*local2global_t)(int numCellsPar, int parIdx, long *globalIdxs);

// For use in kernel tables.
typedef struct { local2global_t kernels[2]; } local2global_kern_loc_list;
typedef struct { local2global_kern_loc_list list[2]; } local2global_kern_bc_list;
typedef struct { local2global_kern_bc_list list[2]; } local2global_kern_list;

// Serendipity local-to-global kernels.
GKYL_CU_D
static const local2global_kern_list ser_loc2glob_list[] = {
  // 1x
  {.list={ 
          // periodicx
          {.list={
                  {fem_parproj_local_to_global_1x_ser_p1_inx_periodicx, fem_parproj_local_to_global_1x_ser_p1_upx_periodicx,}, 
                  {fem_parproj_local_to_global_1x_ser_p2_inx_periodicx, fem_parproj_local_to_global_1x_ser_p2_upx_periodicx,}, 
                 },
          },
          // nonperiodicx
          {.list={ 
                  {fem_parproj_local_to_global_1x_ser_p1_inx_nonperiodicx, fem_parproj_local_to_global_1x_ser_p1_upx_nonperiodicx,}, 
                  {fem_parproj_local_to_global_1x_ser_p2_inx_nonperiodicx, fem_parproj_local_to_global_1x_ser_p2_upx_nonperiodicx,}, 
                 },
          },
         }
  },
  // 2x
  {.list={ 
          // periodicy
          {.list={
                  {fem_parproj_local_to_global_2x_ser_p1_iny_periodicy, fem_parproj_local_to_global_2x_ser_p1_upy_periodicy,}, 
                  {fem_parproj_local_to_global_2x_ser_p2_iny_periodicy, fem_parproj_local_to_global_2x_ser_p2_upy_periodicy,}, 
                 },
          },
          // nonperiodicy
          {.list={ 
                  {fem_parproj_local_to_global_2x_ser_p1_iny_nonperiodicy, fem_parproj_local_to_global_2x_ser_p1_upy_nonperiodicy,}, 
                  {fem_parproj_local_to_global_2x_ser_p2_iny_nonperiodicy, fem_parproj_local_to_global_2x_ser_p2_upy_nonperiodicy,}, 
                 },
          },
         }
  },
  // 3x
  {.list={ 
          // periodicz
          {.list={
                  {fem_parproj_local_to_global_3x_ser_p1_inz_periodicz, fem_parproj_local_to_global_3x_ser_p1_upz_periodicz,}, 
                  {fem_parproj_local_to_global_3x_ser_p2_inz_periodicz, fem_parproj_local_to_global_3x_ser_p2_upz_periodicz,}, 
                 },
          },
          // nonperiodicz
          {.list={ 
                  {fem_parproj_local_to_global_3x_ser_p1_inz_nonperiodicz, fem_parproj_local_to_global_3x_ser_p1_upz_nonperiodicz,}, 
                  {fem_parproj_local_to_global_3x_ser_p2_inz_nonperiodicz, fem_parproj_local_to_global_3x_ser_p2_upz_nonperiodicz,}, 
                 },
          },
         }
  },
};

// Function pointer type for lhs kernels.
typedef void (*lhsstencil_t)(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);

// For use in kernel tables.
typedef struct { lhsstencil_t kernels[3]; } lhsstencil_kern_loc_list;
typedef struct { lhsstencil_kern_loc_list list[2]; } lhsstencil_kern_bc_list;
typedef struct { lhsstencil_kern_bc_list list[2]; } lhsstencil_kern_list;

// Serendipity unweighted lhs kernels.
static const lhsstencil_kern_list ser_lhsstencil_list_noweight[] = {
  // 1x
  {.list={ 
          // nondirichletx
          {.list={
                  {fem_parproj_lhs_stencil_noweight_1x_ser_p1_inx_nondirichletx, fem_parproj_lhs_stencil_noweight_1x_ser_p1_lox_nondirichletx, fem_parproj_lhs_stencil_noweight_1x_ser_p1_upx_nondirichletx,}, 
                  {fem_parproj_lhs_stencil_noweight_1x_ser_p2_inx_nondirichletx, fem_parproj_lhs_stencil_noweight_1x_ser_p2_lox_nondirichletx, fem_parproj_lhs_stencil_noweight_1x_ser_p2_upx_nondirichletx,}, 
                 },
          },
          // dirichletx
          {.list={ 
                  {fem_parproj_lhs_stencil_noweight_1x_ser_p1_inx_nondirichletx, fem_parproj_lhs_stencil_noweight_1x_ser_p1_lox_dirichletx, fem_parproj_lhs_stencil_noweight_1x_ser_p1_upx_dirichletx,}, 
                  {fem_parproj_lhs_stencil_noweight_1x_ser_p2_inx_nondirichletx, fem_parproj_lhs_stencil_noweight_1x_ser_p2_lox_dirichletx, fem_parproj_lhs_stencil_noweight_1x_ser_p2_upx_dirichletx,}, 
                 },
          },
         },
  },
  // 2x
  {.list={ 
          // nondirichlety
          {.list={
                  {fem_parproj_lhs_stencil_noweight_2x_ser_p1_iny_nondirichlety, fem_parproj_lhs_stencil_noweight_2x_ser_p1_loy_nondirichlety, fem_parproj_lhs_stencil_noweight_2x_ser_p1_upy_nondirichlety,}, 
                  {fem_parproj_lhs_stencil_noweight_2x_ser_p2_iny_nondirichlety, fem_parproj_lhs_stencil_noweight_2x_ser_p2_loy_nondirichlety, fem_parproj_lhs_stencil_noweight_2x_ser_p2_upy_nondirichlety,}, 
                 },
          },
          // dirichlety
          {.list={ 
                  {fem_parproj_lhs_stencil_noweight_2x_ser_p1_iny_nondirichlety, fem_parproj_lhs_stencil_noweight_2x_ser_p1_loy_dirichlety, fem_parproj_lhs_stencil_noweight_2x_ser_p1_upy_dirichlety,}, 
                  {fem_parproj_lhs_stencil_noweight_2x_ser_p2_iny_nondirichlety, fem_parproj_lhs_stencil_noweight_2x_ser_p2_loy_dirichlety, fem_parproj_lhs_stencil_noweight_2x_ser_p2_upy_dirichlety,}, 
                 },
          },
         }
  },
  // 3x
  {.list={ 
          // nondirichletz
          {.list={
                  {fem_parproj_lhs_stencil_noweight_3x_ser_p1_inz_nondirichletz, fem_parproj_lhs_stencil_noweight_3x_ser_p1_loz_nondirichletz, fem_parproj_lhs_stencil_noweight_3x_ser_p1_upz_nondirichletz,}, 
                  {fem_parproj_lhs_stencil_noweight_3x_ser_p2_inz_nondirichletz, fem_parproj_lhs_stencil_noweight_3x_ser_p2_loz_nondirichletz, fem_parproj_lhs_stencil_noweight_3x_ser_p2_upz_nondirichletz,}, 
                 },
          },
          // dirichletz
          {.list={ 
                  {fem_parproj_lhs_stencil_noweight_3x_ser_p1_inz_nondirichletz, fem_parproj_lhs_stencil_noweight_3x_ser_p1_loz_dirichletz, fem_parproj_lhs_stencil_noweight_3x_ser_p1_upz_dirichletz,}, 
                  {fem_parproj_lhs_stencil_noweight_3x_ser_p2_inz_nondirichletz, fem_parproj_lhs_stencil_noweight_3x_ser_p2_loz_dirichletz, fem_parproj_lhs_stencil_noweight_3x_ser_p2_upz_dirichletz,}, 
                 },
          },
         }
  },
};

// Serendipity weighted lhs kernels.
static const lhsstencil_kern_list ser_lhsstencil_list_weighted[] = {
  // 1x
  {.list={ 
          // nondirichletx
          {.list={
                  {fem_parproj_lhs_stencil_weighted_1x_ser_p1_inx_nondirichletx, fem_parproj_lhs_stencil_weighted_1x_ser_p1_lox_nondirichletx, fem_parproj_lhs_stencil_weighted_1x_ser_p1_upx_nondirichletx,}, 
                  {fem_parproj_lhs_stencil_weighted_1x_ser_p2_inx_nondirichletx, fem_parproj_lhs_stencil_weighted_1x_ser_p2_lox_nondirichletx, fem_parproj_lhs_stencil_weighted_1x_ser_p2_upx_nondirichletx,}, 
                 },
          },
          // dirichletx
          {.list={ 
                  {fem_parproj_lhs_stencil_weighted_1x_ser_p1_inx_nondirichletx, fem_parproj_lhs_stencil_weighted_1x_ser_p1_lox_dirichletx, fem_parproj_lhs_stencil_weighted_1x_ser_p1_upx_dirichletx,}, 
                  {fem_parproj_lhs_stencil_weighted_1x_ser_p2_inx_nondirichletx, fem_parproj_lhs_stencil_weighted_1x_ser_p2_lox_dirichletx, fem_parproj_lhs_stencil_weighted_1x_ser_p2_upx_dirichletx,}, 
                 },
          },
         }
  },
  // 2x
  {.list={ 
          // nondirichlety
          {.list={
                  {fem_parproj_lhs_stencil_weighted_2x_ser_p1_iny_nondirichlety, fem_parproj_lhs_stencil_weighted_2x_ser_p1_loy_nondirichlety, fem_parproj_lhs_stencil_weighted_2x_ser_p1_upy_nondirichlety,}, 
                  {fem_parproj_lhs_stencil_weighted_2x_ser_p2_iny_nondirichlety, fem_parproj_lhs_stencil_weighted_2x_ser_p2_loy_nondirichlety, fem_parproj_lhs_stencil_weighted_2x_ser_p2_upy_nondirichlety,}, 
                 },
          },
          // dirichlety
          {.list={ 
                  {fem_parproj_lhs_stencil_weighted_2x_ser_p1_iny_nondirichlety, fem_parproj_lhs_stencil_weighted_2x_ser_p1_loy_dirichlety, fem_parproj_lhs_stencil_weighted_2x_ser_p1_upy_dirichlety,}, 
                  {fem_parproj_lhs_stencil_weighted_2x_ser_p2_iny_nondirichlety, fem_parproj_lhs_stencil_weighted_2x_ser_p2_loy_dirichlety, fem_parproj_lhs_stencil_weighted_2x_ser_p2_upy_dirichlety,}, 
                 },
          },
         }
  },
  // 3x
  {.list={ 
          // nondirichletz
          {.list={
                  {fem_parproj_lhs_stencil_weighted_3x_ser_p1_inz_nondirichletz, fem_parproj_lhs_stencil_weighted_3x_ser_p1_loz_nondirichletz, fem_parproj_lhs_stencil_weighted_3x_ser_p1_upz_nondirichletz,}, 
                  {fem_parproj_lhs_stencil_weighted_3x_ser_p2_inz_nondirichletz, fem_parproj_lhs_stencil_weighted_3x_ser_p2_loz_nondirichletz, fem_parproj_lhs_stencil_weighted_3x_ser_p2_upz_nondirichletz,}, 
                 },
          },
          // dirichletz
          {.list={ 
                  {fem_parproj_lhs_stencil_weighted_3x_ser_p1_inz_nondirichletz, fem_parproj_lhs_stencil_weighted_3x_ser_p1_loz_dirichletz, fem_parproj_lhs_stencil_weighted_3x_ser_p1_upz_dirichletz,}, 
                  {fem_parproj_lhs_stencil_weighted_3x_ser_p2_inz_nondirichletz, fem_parproj_lhs_stencil_weighted_3x_ser_p2_loz_dirichletz, fem_parproj_lhs_stencil_weighted_3x_ser_p2_upz_dirichletz,}, 
                 },
          },
         }
  },
};

// Function pointer type for rhs source kernels.
typedef void (*srcstencil_t)(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs,
  double *bsrc);

typedef struct { srcstencil_t kernels[3]; } srcstencil_kern_loc_list;  // For use in kernel tables.
typedef struct { srcstencil_kern_loc_list list[2]; } srcstencil_kern_bc_list;  // For use in kernel tables.
typedef struct { srcstencil_kern_bc_list list[3]; } srcstencil_kern_list;  // For use in kernel tables.

// Serendipity src kernels.
GKYL_CU_D
static const srcstencil_kern_list ser_srcstencil_list_noweight[] = {
  // 1x
  {.list={
          // nondirichletx
          {.list={
                  {fem_parproj_src_stencil_noweight_1x_ser_p1_inx_nondirichletx, fem_parproj_src_stencil_noweight_1x_ser_p1_lox_nondirichletx, fem_parproj_src_stencil_noweight_1x_ser_p1_upx_nondirichletx,},
                  {fem_parproj_src_stencil_noweight_1x_ser_p2_inx_nondirichletx, fem_parproj_src_stencil_noweight_1x_ser_p2_lox_nondirichletx, fem_parproj_src_stencil_noweight_1x_ser_p2_upx_nondirichletx,},
                 },
          },
          // dirichletx ghost
          {.list={
                  {fem_parproj_src_stencil_noweight_1x_ser_p1_inx_nondirichletx, fem_parproj_src_stencil_noweight_1x_ser_p1_lox_dirichlet_ghostx, fem_parproj_src_stencil_noweight_1x_ser_p1_upx_dirichlet_ghostx,},
                  {fem_parproj_src_stencil_noweight_1x_ser_p2_inx_nondirichletx, fem_parproj_src_stencil_noweight_1x_ser_p2_lox_dirichlet_ghostx, fem_parproj_src_stencil_noweight_1x_ser_p2_upx_dirichlet_ghostx,},
                 },
          },
          // dirichletx skin
          {.list={
                  {fem_parproj_src_stencil_noweight_1x_ser_p1_inx_nondirichletx, fem_parproj_src_stencil_noweight_1x_ser_p1_lox_dirichlet_skinx, fem_parproj_src_stencil_noweight_1x_ser_p1_upx_dirichlet_skinx,},
                  {fem_parproj_src_stencil_noweight_1x_ser_p2_inx_nondirichletx, fem_parproj_src_stencil_noweight_1x_ser_p2_lox_dirichlet_skinx, fem_parproj_src_stencil_noweight_1x_ser_p2_upx_dirichlet_skinx,},
                 },
          },
         }
  },
  // 2x
  {.list={
          // nondirichlety
          {.list={
                  {fem_parproj_src_stencil_noweight_2x_ser_p1_iny_nondirichlety, fem_parproj_src_stencil_noweight_2x_ser_p1_loy_nondirichlety, fem_parproj_src_stencil_noweight_2x_ser_p1_upy_nondirichlety,},
                  {fem_parproj_src_stencil_noweight_2x_ser_p2_iny_nondirichlety, fem_parproj_src_stencil_noweight_2x_ser_p2_loy_nondirichlety, fem_parproj_src_stencil_noweight_2x_ser_p2_upy_nondirichlety,},
                 },
          },
          // dirichlety ghost
          {.list={
                  {fem_parproj_src_stencil_noweight_2x_ser_p1_iny_nondirichlety, fem_parproj_src_stencil_noweight_2x_ser_p1_loy_dirichlet_ghosty, fem_parproj_src_stencil_noweight_2x_ser_p1_upy_dirichlet_ghosty,},
                  {fem_parproj_src_stencil_noweight_2x_ser_p2_iny_nondirichlety, fem_parproj_src_stencil_noweight_2x_ser_p2_loy_dirichlet_ghosty, fem_parproj_src_stencil_noweight_2x_ser_p2_upy_dirichlet_ghosty,},
                 },
          },
          // dirichlety skin
          {.list={
                  {fem_parproj_src_stencil_noweight_2x_ser_p1_iny_nondirichlety, fem_parproj_src_stencil_noweight_2x_ser_p1_loy_dirichlet_skiny, fem_parproj_src_stencil_noweight_2x_ser_p1_upy_dirichlet_skiny,},
                  {fem_parproj_src_stencil_noweight_2x_ser_p2_iny_nondirichlety, fem_parproj_src_stencil_noweight_2x_ser_p2_loy_dirichlet_skiny, fem_parproj_src_stencil_noweight_2x_ser_p2_upy_dirichlet_skiny,},
                 },
          },
         }
  },
  // 3x
  {.list={
          // nondirichletz
          {.list={
                  {fem_parproj_src_stencil_noweight_3x_ser_p1_inz_nondirichletz, fem_parproj_src_stencil_noweight_3x_ser_p1_loz_nondirichletz, fem_parproj_src_stencil_noweight_3x_ser_p1_upz_nondirichletz,},
                  {fem_parproj_src_stencil_noweight_3x_ser_p2_inz_nondirichletz, fem_parproj_src_stencil_noweight_3x_ser_p2_loz_nondirichletz, fem_parproj_src_stencil_noweight_3x_ser_p2_upz_nondirichletz,},
                 },
          },
          // dirichletz ghost
          {.list={
                  {fem_parproj_src_stencil_noweight_3x_ser_p1_inz_nondirichletz, fem_parproj_src_stencil_noweight_3x_ser_p1_loz_dirichlet_ghostz, fem_parproj_src_stencil_noweight_3x_ser_p1_upz_dirichlet_ghostz,},
                  {fem_parproj_src_stencil_noweight_3x_ser_p2_inz_nondirichletz, fem_parproj_src_stencil_noweight_3x_ser_p2_loz_dirichlet_ghostz, fem_parproj_src_stencil_noweight_3x_ser_p2_upz_dirichlet_ghostz,},
                 },
          },
          // dirichletz skin
          {.list={
                  {fem_parproj_src_stencil_noweight_3x_ser_p1_inz_nondirichletz, fem_parproj_src_stencil_noweight_3x_ser_p1_loz_dirichlet_skinz, fem_parproj_src_stencil_noweight_3x_ser_p1_upz_dirichlet_skinz,},
                  {fem_parproj_src_stencil_noweight_3x_ser_p2_inz_nondirichletz, fem_parproj_src_stencil_noweight_3x_ser_p2_loz_dirichlet_skinz, fem_parproj_src_stencil_noweight_3x_ser_p2_upz_dirichlet_skinz,},
                 },
          },
         }
  },
};

GKYL_CU_D
static const srcstencil_kern_list ser_srcstencil_list_weighted[] = {
  // 1x
  {.list={
          // nondirichletx
          {.list={
                  {fem_parproj_src_stencil_weighted_1x_ser_p1_inx_nondirichletx, fem_parproj_src_stencil_weighted_1x_ser_p1_lox_nondirichletx, fem_parproj_src_stencil_weighted_1x_ser_p1_upx_nondirichletx,},
                  {fem_parproj_src_stencil_weighted_1x_ser_p2_inx_nondirichletx, fem_parproj_src_stencil_weighted_1x_ser_p2_lox_nondirichletx, fem_parproj_src_stencil_weighted_1x_ser_p2_upx_nondirichletx,},
                 },
          },
          // dirichletx ghost
          {.list={
                  {fem_parproj_src_stencil_weighted_1x_ser_p1_inx_nondirichletx, fem_parproj_src_stencil_weighted_1x_ser_p1_lox_dirichlet_ghostx, fem_parproj_src_stencil_weighted_1x_ser_p1_upx_dirichlet_ghostx,},
                  {fem_parproj_src_stencil_weighted_1x_ser_p2_inx_nondirichletx, fem_parproj_src_stencil_weighted_1x_ser_p2_lox_dirichlet_ghostx, fem_parproj_src_stencil_weighted_1x_ser_p2_upx_dirichlet_ghostx,},
                 },
          },
          // dirichletx skin
          {.list={
                  {fem_parproj_src_stencil_weighted_1x_ser_p1_inx_nondirichletx, fem_parproj_src_stencil_weighted_1x_ser_p1_lox_dirichlet_skinx, fem_parproj_src_stencil_weighted_1x_ser_p1_upx_dirichlet_skinx,},
                  {fem_parproj_src_stencil_weighted_1x_ser_p2_inx_nondirichletx, fem_parproj_src_stencil_weighted_1x_ser_p2_lox_dirichlet_skinx, fem_parproj_src_stencil_weighted_1x_ser_p2_upx_dirichlet_skinx,},
                 },
          },
         }
  },
  // 2x
  {.list={
          // nondirichlety
          {.list={
                  {fem_parproj_src_stencil_weighted_2x_ser_p1_iny_nondirichlety, fem_parproj_src_stencil_weighted_2x_ser_p1_loy_nondirichlety, fem_parproj_src_stencil_weighted_2x_ser_p1_upy_nondirichlety,},
                  {fem_parproj_src_stencil_weighted_2x_ser_p2_iny_nondirichlety, fem_parproj_src_stencil_weighted_2x_ser_p2_loy_nondirichlety, fem_parproj_src_stencil_weighted_2x_ser_p2_upy_nondirichlety,},
                 },
          },
          // dirichlety ghost
          {.list={
                  {fem_parproj_src_stencil_weighted_2x_ser_p1_iny_nondirichlety, fem_parproj_src_stencil_weighted_2x_ser_p1_loy_dirichlet_ghosty, fem_parproj_src_stencil_weighted_2x_ser_p1_upy_dirichlet_ghosty,},
                  {fem_parproj_src_stencil_weighted_2x_ser_p2_iny_nondirichlety, fem_parproj_src_stencil_weighted_2x_ser_p2_loy_dirichlet_ghosty, fem_parproj_src_stencil_weighted_2x_ser_p2_upy_dirichlet_ghosty,},
                 },
          },
          // dirichlety skin
          {.list={
                  {fem_parproj_src_stencil_weighted_2x_ser_p1_iny_nondirichlety, fem_parproj_src_stencil_weighted_2x_ser_p1_loy_dirichlet_skiny, fem_parproj_src_stencil_weighted_2x_ser_p1_upy_dirichlet_skiny,},
                  {fem_parproj_src_stencil_weighted_2x_ser_p2_iny_nondirichlety, fem_parproj_src_stencil_weighted_2x_ser_p2_loy_dirichlet_skiny, fem_parproj_src_stencil_weighted_2x_ser_p2_upy_dirichlet_skiny,},
                 },
          },
         }
  },
  // 3x
  {.list={
          // nondirichletz
          {.list={
                  {fem_parproj_src_stencil_weighted_3x_ser_p1_inz_nondirichletz, fem_parproj_src_stencil_weighted_3x_ser_p1_loz_nondirichletz, fem_parproj_src_stencil_weighted_3x_ser_p1_upz_nondirichletz,},
                  {fem_parproj_src_stencil_weighted_3x_ser_p2_inz_nondirichletz, fem_parproj_src_stencil_weighted_3x_ser_p2_loz_nondirichletz, fem_parproj_src_stencil_weighted_3x_ser_p2_upz_nondirichletz,},
                 },
          },
          // dirichletz ghost
          {.list={
                  {fem_parproj_src_stencil_weighted_3x_ser_p1_inz_nondirichletz, fem_parproj_src_stencil_weighted_3x_ser_p1_loz_dirichlet_ghostz, fem_parproj_src_stencil_weighted_3x_ser_p1_upz_dirichlet_ghostz,},
                  {fem_parproj_src_stencil_weighted_3x_ser_p2_inz_nondirichletz, fem_parproj_src_stencil_weighted_3x_ser_p2_loz_dirichlet_ghostz, fem_parproj_src_stencil_weighted_3x_ser_p2_upz_dirichlet_ghostz,},
                 },
          },
          // dirichletz skin
          {.list={
                  {fem_parproj_src_stencil_weighted_3x_ser_p1_inz_nondirichletz, fem_parproj_src_stencil_weighted_3x_ser_p1_loz_dirichlet_skinz, fem_parproj_src_stencil_weighted_3x_ser_p1_upz_dirichlet_skinz,},
                  {fem_parproj_src_stencil_weighted_3x_ser_p2_inz_nondirichletz, fem_parproj_src_stencil_weighted_3x_ser_p2_loz_dirichlet_skinz, fem_parproj_src_stencil_weighted_3x_ser_p2_upz_dirichlet_skinz,},
                 },
          },
         }
  },
};

// Function pointer type for kernels that convert the solution from nodal to
// modal.
typedef void (*solstencil_t)(const double *sol_nodal_global, long nodeOff,
  const long *globalIdxs, double *sol_modal_local);

typedef struct { solstencil_t kernels[3]; } solstencil_kern_list;  // For use in kernel tables.

// Serendipity sol kernels.
GKYL_CU_D
static const solstencil_kern_list ser_solstencil_list[] = {
  { fem_parproj_sol_stencil_1x_ser_p1, fem_parproj_sol_stencil_1x_ser_p2 },
  { fem_parproj_sol_stencil_2x_ser_p1, fem_parproj_sol_stencil_2x_ser_p2 },
  { fem_parproj_sol_stencil_3x_ser_p1, fem_parproj_sol_stencil_3x_ser_p2 }
};

// Function pointer type for kernels that enforce biasing in LHS matrix.
typedef void (*bias_lhs_t)(const int *edge, const int *perp_dirs, const long *globalIdxs, gkyl_mat_triples *tri);

// For use in kernel tables.
typedef struct { bias_lhs_t kernels[2]; } bias_lhs_kern_loc_list;
typedef struct { bias_lhs_kern_loc_list list[2]; } bias_lhs_kern_bc_list;
typedef struct { bias_lhs_kern_bc_list list[2]; } bias_lhs_kern_dim_list;

// Serendipity bias_lhs kernels.
static const bias_lhs_kern_dim_list ser_bias_lhs_list[] = {
  // 1x
  {.list = 
    { 
      // periodicy
      { .list = {
          {NULL, NULL},
          {NULL, NULL},
        },
      },
      // nonperiodicy
      { .list = {
          {NULL, NULL},
          {NULL, NULL},
        },
      },
    },
  },
  // 2x
  {.list = 
    { 
      // periodicy
      { .list = {
          {fem_parproj_bias_line_lhs_2x_ser_p1_iny_periodicy, fem_parproj_bias_line_lhs_2x_ser_p1_upy_periodicy},
          {NULL, NULL},
        },
      },
      // nonperiodicy
      { .list = {
          {fem_parproj_bias_line_lhs_2x_ser_p1_iny_periodicy, fem_parproj_bias_line_lhs_2x_ser_p1_upy_nonperiodicy},
          {NULL, NULL},
        },
      },
    },
  },
  // 3x
  {.list = 
    { 
      // periodicz
      { .list = {
          {fem_parproj_bias_line_lhs_3x_ser_p1_inz_periodicz, fem_parproj_bias_line_lhs_3x_ser_p1_upz_periodicz},
          {NULL, NULL},
        },
      },
      // nonperiodicz
      { .list = {
          {fem_parproj_bias_line_lhs_3x_ser_p1_inz_periodicz, fem_parproj_bias_line_lhs_3x_ser_p1_upz_nonperiodicz},
          {NULL, NULL},
        },
      },
    },
  },
            
};

// Function pointer type for kernels that enforce biasing in RHS source.
typedef void (*bias_src_t)(const int *edge, const int *perp_dirs, double val, long nodeOff, const long *globalIdxs, double *bsrc);

// For use in kernel tables.
typedef struct { bias_src_t kernels[2]; } bias_src_kern_loc_list;
typedef struct { bias_src_kern_loc_list list[2]; } bias_src_kern_bc_list;
typedef struct { bias_src_kern_bc_list list[2]; } bias_src_kern_dim_list;

// Serendipity bias_src kernels.
GKYL_CU_D
static const bias_src_kern_dim_list ser_bias_src_list[] = {
  // 1x
  {.list = 
    { 
      // periodicy
      { .list =
        {
          {NULL, NULL},
          {NULL, NULL},
        },
      },
      // nonperiodicy
      { .list =
        {
          {NULL, NULL},
          {NULL, NULL},
        },
      },
    },
  },
  // 2x
  {.list = 
    { 
      // periodicy
      { .list =
        {
          {fem_parproj_bias_line_src_2x_ser_p1_iny_periodicy, fem_parproj_bias_line_src_2x_ser_p1_upy_periodicy},
          {NULL, NULL},
        },
      },
      // nonperiodicy
      { .list =
        {
          {fem_parproj_bias_line_src_2x_ser_p1_iny_periodicy, fem_parproj_bias_line_src_2x_ser_p1_upy_nonperiodicy},
          {NULL, NULL},
        },
      },
    },
  },
  // 3x
  {.list = 
    { 
      // periodicz
      { .list =
        {
          {fem_parproj_bias_line_src_3x_ser_p1_inz_periodicz, fem_parproj_bias_line_src_3x_ser_p1_upz_periodicz},
          {NULL, NULL},
        },
      },
      // nonperiodicz
      { .list =
        {
          {fem_parproj_bias_line_src_3x_ser_p1_inz_periodicz, fem_parproj_bias_line_src_3x_ser_p1_upz_nonperiodicz},
          {NULL, NULL},
        },
      }
    }
  },
};

// Functions that return the value to impose as Dirichlet BC.
typedef const double *(*get_diri_val_t)(int par_dir, int par_num_cells,
  const int *idx, const struct gkyl_range *solve_range, const struct gkyl_array *phibc);

// No Dirichlet BC.
GKYL_CU_D
static const double *get_dirichlet_value_disabled(int par_dir, int par_num_cells,
  const int *idx, const struct gkyl_range *solve_range, const struct gkyl_array *phibc)
{
  return 0;
}

// Dirichlet BC using the ghost value.
GKYL_CU_D
static const double *get_dirichlet_value_enabled_ghost(int par_dir, int par_num_cells,
  const int *idx, const struct gkyl_range *solve_range, const struct gkyl_array *phibc)
{
  int dirichlet_idx[GKYL_MAX_CDIM];
  for (size_t d=0; d<par_dir+1; d++)
    dirichlet_idx[d] = idx[d];

  dirichlet_idx[par_dir] = dirichlet_idx[par_dir] == par_num_cells? dirichlet_idx[par_dir]+1 : dirichlet_idx[par_dir]-1;
  return (const double *) gkyl_array_cfetch(phibc, gkyl_range_idx(solve_range, dirichlet_idx));
}

// Dirichlet BC using the skin value.
GKYL_CU_D
static const double *get_dirichlet_value_enabled_skin(int par_dir, int par_num_cells,
  const int *idx, const struct gkyl_range *solve_range, const struct gkyl_array *phibc)
{
  int dirichlet_idx[GKYL_MAX_CDIM];
  for (size_t d=0; d<par_dir+1; d++)
    dirichlet_idx[d] = idx[d];

  return (const double *) gkyl_array_cfetch(phibc, gkyl_range_idx(solve_range, dirichlet_idx));
}

// Struct containing pointers to the various kernels. Needed to create a similar struct on the GPU.
struct gkyl_fem_parproj_kernels {
  local2global_t l2g[2];  // Pointer to local-to-global kernel.

  lhsstencil_t lhsker[3];  // Weighted LHS kernel.

  srcstencil_t srcker[3];  // RHS source kernel.

  solstencil_t solker;  // Kernel that takes the solution and converts it to modal.

  get_diri_val_t get_dirichlet_value; // Gets value to impose as Dirichlet BC.

  // Pointer to kernels that enforce biasing, 2 (interior and upper).
  bias_lhs_t bias_lhs_ker[2];
  bias_src_t bias_src_ker[2];
};

// Type of function used to enforce biasing in the RHS src.
typedef void (*bias_src_func_t)(gkyl_fem_parproj* up, const struct gkyl_array *rhsin);

struct gkyl_fem_parproj {
  int ndim; // Grid's number of dimensions.
  int num_basis; // Number of basis functions.
  enum gkyl_basis_type basis_type; // Type of DG basis.
  int poly_order; // Polynomial object of DG basis.
  struct gkyl_rect_grid grid; // Grid object.
  int pardir; // Parallel (z) direction.
  int parnum_cells; // Number of cells in parallel (z) direction.
  bool isperiodic; // =true if parallel direction is periodic.
  bool isdirichlet; // =true if parallel direction has periodic BCs.
  bool has_weight_rhs; // Whether there's a weight on the RHS.
  struct gkyl_array *weight_rhs; // The RHS weight.

  const struct gkyl_range *solve_range;
  struct gkyl_range perp_range2d; // 2D range of perpendicular cells.
  struct gkyl_range par_range1d; // 1D range of parallel cells.
  struct gkyl_range_iter perp_iter2d;
  struct gkyl_range_iter par_iter1d;

  int numnodes_local;
  long numnodes_global;

  struct gkyl_array *brhs;

  struct gkyl_superlu_prob* prob;
#ifdef GKYL_HAVE_CUDA
  struct gkyl_culinsolver_prob* prob_cu;
  struct gkyl_array *brhs_cu;
#endif

  long *globalidx;

  struct gkyl_fem_parproj_kernels *kernels;
  bool use_gpu;  

  int bl_ndim_perp; // Number of perpendicular directions of bias lines.
  int num_bias_line; // Number of biased lines.
  struct gkyl_poisson_bias_line *bias_lines; // Biased lines.
  bias_src_func_t bias_line_src; // Function to enforce biasing in RHS source.
};

// "Choose Kernel" based on polyorder, stencil location and BCs.
#define CK(lst,dim,bc,poly_order,loc) lst[dim-1].list[bc].list[poly_order-1].kernels[loc]

#ifdef GKYL_HAVE_CUDA
void fem_parproj_choose_kernels_cu(const struct gkyl_basis* basis, bool has_weight_lhs, bool has_weight_rhs,
  enum gkyl_fem_parproj_bc_type bctype, struct gkyl_fem_parproj_kernels *kers);

/**
 * Assign the right-side vector with the discontinuous (DG) source field
 * on the NVIDIA GPU.
 *
 * @param up FEM project updater to run.
 * @param rhsin DG field to set as RHS source.
 * @param phibc Potential to use for Dirichlet BCs (only use ghost cells).
 */
void gkyl_fem_parproj_set_rhs_cu(struct gkyl_fem_parproj *up, const struct gkyl_array *rhsin, const struct gkyl_array *phibc);

/**
 * Replace the entries in the RHS src vector with the biased potential values.
 *
 * @param up FEM parallel projection updater to run.
 * @param rhsin DG field to set as RHS source.
 */
void
gkyl_fem_parproj_bias_src_enabled_cu(gkyl_fem_parproj *up, const struct gkyl_array *rhsin);

/**
 * Solve the linear problem
 * on the NVIDIA GPU.
 *
 * @param up FEM project updater to run.
 */
void gkyl_fem_parproj_solve_cu(struct gkyl_fem_parproj* up, struct gkyl_array *phiout);
#endif

GKYL_CU_D
static void
fem_parproj_choose_local2global_kernel(const struct gkyl_basis *basis,
  enum gkyl_fem_parproj_bc_type bctype, local2global_t *l2gout)
{
  int bckey[1] = {-1};
  bckey[0] = bctype == GKYL_FEM_PARPROJ_PERIODIC? 0 : 1;

  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      for (int k=0; k<2; k++)
        l2gout[k] = CK(ser_loc2glob_list, basis->ndim, bckey[0], basis->poly_order, k);
      break;
    default:
      assert(false);
      break;
  }
}

GKYL_CU_D
static void 
fem_parproj_choose_lhs_kernel(const struct gkyl_basis *basis,
  enum gkyl_fem_parproj_bc_type bctype, bool isweighted, lhsstencil_t *lhsout)
{
  int bckey[1] = {-1};
  if (bctype == GKYL_FEM_PARPROJ_DIRICHLET_GHOST || bctype == GKYL_FEM_PARPROJ_DIRICHLET_SKIN)
    bckey[0] = 1;
  else
    bckey[0] = 0;

  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      for (int k=0; k<3; k++)
        lhsout[k] = isweighted? CK(ser_lhsstencil_list_weighted, basis->ndim, bckey[0], basis->poly_order, k)
                              : CK(ser_lhsstencil_list_noweight, basis->ndim, bckey[0], basis->poly_order, k); 
      break;
    default:
      assert(false);
      break;
  }
}

GKYL_CU_D
static void
fem_parproj_choose_srcstencil_kernel(const struct gkyl_basis *basis,
  enum gkyl_fem_parproj_bc_type bctype, bool isweighted, srcstencil_t *srcout)
{
  int bckey[1] = {-1};
  if (bctype == GKYL_FEM_PARPROJ_DIRICHLET_GHOST)
    bckey[0] = 1;
  else if (bctype == GKYL_FEM_PARPROJ_DIRICHLET_SKIN)
    bckey[0] = 2;
  else
    bckey[0] = 0;

  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      for (int k=0; k<3; k++)
        srcout[k] = isweighted? CK(ser_srcstencil_list_weighted, basis->ndim, bckey[0], basis->poly_order, k)
                              : CK(ser_srcstencil_list_noweight, basis->ndim, bckey[0], basis->poly_order, k);
      break;
    default:
      assert(false);
      break;
  }
}

GKYL_CU_D
static solstencil_t
fem_parproj_choose_solstencil_kernel(const struct gkyl_basis *basis)
{
  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_solstencil_list[basis->ndim-1].kernels[basis->poly_order-1];
    default:
      assert(false);
      break;
  }
  return 0;
}

GKYL_CU_D
static void
fem_parproj_choose_bias_lhs_kernels(const struct gkyl_basis* basis,
  enum gkyl_fem_parproj_bc_type bctype, bias_lhs_t *blhs_out)
{
  int poly_order = basis->poly_order;
  int ndim = basis->ndim;

  int bckey[1];
  bckey[0] = bctype == GKYL_FEM_PARPROJ_PERIODIC? 0 : 1;

  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      for (int k=0; k<2; k++)
        blhs_out[k] = CK(ser_bias_lhs_list, ndim, bckey[0], poly_order, k);

      break;
//    case GKYL_BASIS_MODAL_TENSOR:
//      break;
    default:
      assert(false);
      break;
  }
}

GKYL_CU_D
static void
fem_parproj_choose_bias_src_kernels(const struct gkyl_basis* basis,
  enum gkyl_fem_parproj_bc_type bctype, bias_src_t *bsrc_out)
{
  int poly_order = basis->poly_order;
  int ndim = basis->ndim;

  int bckey[1];
  bckey[0] = bctype == GKYL_FEM_PARPROJ_PERIODIC? 0 : 1;

  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      for (int k=0; k<2; k++)
        bsrc_out[k] = CK(ser_bias_src_list, ndim, bckey[0], poly_order, k); 

      break;
//    case GKYL_BASIS_MODAL_TENSOR:
//      break;
    default:
      assert(false);
      break;
  }
}

GKYL_CU_D
static void
fem_parproj_choose_kernels(const struct gkyl_basis* basis, bool has_weight_lhs, bool has_weight_rhs,
  enum gkyl_fem_parproj_bc_type bctype, bool use_gpu, struct gkyl_fem_parproj_kernels *kers)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    fem_parproj_choose_kernels_cu(basis, has_weight_lhs, has_weight_rhs, bctype, kers);
    return;
  }
#endif

  // Select local-to-global mapping kernel:
  fem_parproj_choose_local2global_kernel(basis, bctype, kers->l2g);

  // Select weighted LHS kernel (not always used):
  fem_parproj_choose_lhs_kernel(basis, bctype, has_weight_lhs, kers->lhsker);

  // Select RHS source kernel:
  fem_parproj_choose_srcstencil_kernel(basis, bctype, has_weight_rhs, kers->srcker);

  // Select kernel that fetches the solution:
  kers->solker = fem_parproj_choose_solstencil_kernel(basis);

  // Select function that obtains the value to impose as Dirichlet BC.
  if (bctype == GKYL_FEM_PARPROJ_DIRICHLET_GHOST)
    kers->get_dirichlet_value = get_dirichlet_value_enabled_ghost;
  else if (bctype == GKYL_FEM_PARPROJ_DIRICHLET_SKIN)
    kers->get_dirichlet_value = get_dirichlet_value_enabled_skin;
  else
    kers->get_dirichlet_value = get_dirichlet_value_disabled;

  // Select biasing kernels:
  fem_parproj_choose_bias_lhs_kernels(basis, bctype, kers->bias_lhs_ker);
  fem_parproj_choose_bias_src_kernels(basis, bctype, kers->bias_src_ker);
}

GKYL_CU_DH
static inline int idx_to_inloup_ker(int num_cells, int idx) {
  // Return the index of the kernel (in the array of kernels) needed given the grid index.
  // This function is for kernels that differentiate between lower, interior
  // and upper cells.
  int iout = 0;
  if (idx == 1)
    iout = 1;
  else if (idx == num_cells)
    iout = 2;
  return iout;
}
