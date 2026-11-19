#pragma once

#include <assert.h>
#include <gkyl_fem_proj.h>
#include <gkyl_fem_proj_kernels.h>
#include <gkyl_superlu_ops.h>
#ifdef GKYL_HAVE_CUDA
#include <gkyl_culinsolver_ops.h>
#endif

// Largest serendipity basis planned here: 3D p2 (20 nodes).
#define FEM_PROJ_MAX_NODES 20

typedef long (*fem_proj_num_nodes_t)(const int *, const bool *);
typedef void (*fem_proj_local_to_global_t)(const int *, const int *, long *);
typedef void (*fem_proj_mass_t)(const double *, const long *, gkyl_mat_triples *);
typedef void (*fem_proj_src_t)(const double *, const double *, const double *, const double *, const long *, double *);
typedef void (*fem_proj_bias_lhs_t)(const long *, gkyl_mat_triples *);
typedef void (*fem_proj_bias_src_t)(double, const long *, double *);
typedef void (*fem_proj_transform_t)(const double *, double *);

struct fem_proj_kernels {
  fem_proj_num_nodes_t num_nodes; // Global node count; host initialization only.
  fem_proj_local_to_global_t
    local_to_global[1 << GKYL_MAX_CDIM]; // Compact indices by upper-boundary flags.
  fem_proj_mass_t mass[1 << (2 * GKYL_MAX_CDIM)]; // Host matrix assembly by cell location.
  fem_proj_src_t src[1 << (2 * GKYL_MAX_CDIM)]; // RHS assembly by cell location.
  fem_proj_bias_lhs_t bias_lhs[2 * GKYL_MAX_CDIM]; // Host constraints by cell face.
  fem_proj_bias_src_t bias_src[2 * GKYL_MAX_CDIM]; // RHS constraints by cell face.
  fem_proj_transform_t sol; // Nodal-to-modal reconstruction.
};

// Local-to-global kernels indexed by periodicity mask, polynomial order, and upper-boundary mask.
GKYL_CU_D static const fem_proj_local_to_global_t fem_proj_l2g_1x[2][3][2] = {
  // nonperiodic x.
  {{NULL}, // p0 is unsupported.
   {// p1.
    fem_proj_1x_ser_p1_local_to_global_inx_nonperiodicx,
    fem_proj_1x_ser_p1_local_to_global_upx_nonperiodicx
   },
   {// p2.
    fem_proj_1x_ser_p2_local_to_global_inx_nonperiodicx,
    fem_proj_1x_ser_p2_local_to_global_upx_nonperiodicx
   }},
  // periodic x.
  {{NULL}, // p0 is unsupported.
   {// p1.
    fem_proj_1x_ser_p1_local_to_global_inx_periodicx,
    fem_proj_1x_ser_p1_local_to_global_upx_periodicx
   },
   {// p2.
    fem_proj_1x_ser_p2_local_to_global_inx_periodicx,
    fem_proj_1x_ser_p2_local_to_global_upx_periodicx
   }}
};

// Local-to-global kernels indexed by periodicity mask, polynomial order, and upper-boundary mask.
GKYL_CU_D static const fem_proj_local_to_global_t fem_proj_l2g_2x[4][3][4] = {
  // nonperiodic x, nonperiodic y.
  {{NULL}, // p0 is unsupported.
   {// p1.
    fem_proj_2x_ser_p1_local_to_global_inx_nonperiodicx_iny_nonperiodicy,
    fem_proj_2x_ser_p1_local_to_global_upx_nonperiodicx_iny_nonperiodicy,
    fem_proj_2x_ser_p1_local_to_global_inx_nonperiodicx_upy_nonperiodicy,
    fem_proj_2x_ser_p1_local_to_global_upx_nonperiodicx_upy_nonperiodicy
   },
   {// p2.
    fem_proj_2x_ser_p2_local_to_global_inx_nonperiodicx_iny_nonperiodicy,
    fem_proj_2x_ser_p2_local_to_global_upx_nonperiodicx_iny_nonperiodicy,
    fem_proj_2x_ser_p2_local_to_global_inx_nonperiodicx_upy_nonperiodicy,
    fem_proj_2x_ser_p2_local_to_global_upx_nonperiodicx_upy_nonperiodicy
   }},
  // periodic x, nonperiodic y.
  {{NULL}, // p0 is unsupported.
   {// p1.
    fem_proj_2x_ser_p1_local_to_global_inx_periodicx_iny_nonperiodicy,
    fem_proj_2x_ser_p1_local_to_global_upx_periodicx_iny_nonperiodicy,
    fem_proj_2x_ser_p1_local_to_global_inx_periodicx_upy_nonperiodicy,
    fem_proj_2x_ser_p1_local_to_global_upx_periodicx_upy_nonperiodicy
   },
   {// p2.
    fem_proj_2x_ser_p2_local_to_global_inx_periodicx_iny_nonperiodicy,
    fem_proj_2x_ser_p2_local_to_global_upx_periodicx_iny_nonperiodicy,
    fem_proj_2x_ser_p2_local_to_global_inx_periodicx_upy_nonperiodicy,
    fem_proj_2x_ser_p2_local_to_global_upx_periodicx_upy_nonperiodicy
   }},
  // nonperiodic x, periodic y.
  {{NULL}, // p0 is unsupported.
   {// p1.
    fem_proj_2x_ser_p1_local_to_global_inx_nonperiodicx_iny_periodicy,
    fem_proj_2x_ser_p1_local_to_global_upx_nonperiodicx_iny_periodicy,
    fem_proj_2x_ser_p1_local_to_global_inx_nonperiodicx_upy_periodicy,
    fem_proj_2x_ser_p1_local_to_global_upx_nonperiodicx_upy_periodicy
   },
   {// p2.
    fem_proj_2x_ser_p2_local_to_global_inx_nonperiodicx_iny_periodicy,
    fem_proj_2x_ser_p2_local_to_global_upx_nonperiodicx_iny_periodicy,
    fem_proj_2x_ser_p2_local_to_global_inx_nonperiodicx_upy_periodicy,
    fem_proj_2x_ser_p2_local_to_global_upx_nonperiodicx_upy_periodicy
   }},
  // periodic x, periodic y.
  {{NULL}, // p0 is unsupported.
   {// p1.
    fem_proj_2x_ser_p1_local_to_global_inx_periodicx_iny_periodicy,
    fem_proj_2x_ser_p1_local_to_global_upx_periodicx_iny_periodicy,
    fem_proj_2x_ser_p1_local_to_global_inx_periodicx_upy_periodicy,
    fem_proj_2x_ser_p1_local_to_global_upx_periodicx_upy_periodicy
   },
   {// p2.
    fem_proj_2x_ser_p2_local_to_global_inx_periodicx_iny_periodicy,
    fem_proj_2x_ser_p2_local_to_global_upx_periodicx_iny_periodicy,
    fem_proj_2x_ser_p2_local_to_global_inx_periodicx_upy_periodicy,
    fem_proj_2x_ser_p2_local_to_global_upx_periodicx_upy_periodicy
   }}
};

// Matrix assembly is host-only, including for GPU solves.
#ifndef __CUDA_ARCH__
// 1D mass kernels: polynomial order, weight flag, then active face code.
// Face code uses radix 2, ordered lower x, upper x, lower y, upper y.
static const fem_proj_mass_t fem_proj_mass_1x[3][2][4] = {
  {{NULL}}, // p0 is unsupported.
  {// p1.
   {// noweight.
    fem_proj_1x_ser_p1_mass_noweight_inx, fem_proj_1x_ser_p1_mass_noweight_lox_dirichletx,
    fem_proj_1x_ser_p1_mass_noweight_upx_dirichletx,
    fem_proj_1x_ser_p1_mass_noweight_lox_dirichletx_upx_dirichletx
   },
   {// weighted.
    fem_proj_1x_ser_p1_mass_weighted_inx, fem_proj_1x_ser_p1_mass_weighted_lox_dirichletx,
    fem_proj_1x_ser_p1_mass_weighted_upx_dirichletx,
    fem_proj_1x_ser_p1_mass_weighted_lox_dirichletx_upx_dirichletx
   }
  },
  {// p2.
   {// noweight.
    fem_proj_1x_ser_p2_mass_noweight_inx, fem_proj_1x_ser_p2_mass_noweight_lox_dirichletx,
    fem_proj_1x_ser_p2_mass_noweight_upx_dirichletx,
    fem_proj_1x_ser_p2_mass_noweight_lox_dirichletx_upx_dirichletx
   },
   {// weighted.
    fem_proj_1x_ser_p2_mass_weighted_inx, fem_proj_1x_ser_p2_mass_weighted_lox_dirichletx,
    fem_proj_1x_ser_p2_mass_weighted_upx_dirichletx,
    fem_proj_1x_ser_p2_mass_weighted_lox_dirichletx_upx_dirichletx
   }
  }
};

// 2D mass kernels: polynomial order, weight flag, then active face code.
// Face code uses radix 2, ordered lower x, upper x, lower y, upper y.
static const fem_proj_mass_t fem_proj_mass_2x[3][2][16] = {
  {{NULL}}, // p0 is unsupported.
  {// p1.
   {// noweight.
    fem_proj_2x_ser_p1_mass_noweight_inx_iny, fem_proj_2x_ser_p1_mass_noweight_lox_dirichletx_iny,
    fem_proj_2x_ser_p1_mass_noweight_upx_dirichletx_iny,
    fem_proj_2x_ser_p1_mass_noweight_lox_dirichletx_upx_dirichletx_iny,
    fem_proj_2x_ser_p1_mass_noweight_inx_loy_dirichlety,
    fem_proj_2x_ser_p1_mass_noweight_lox_dirichletx_loy_dirichlety,
    fem_proj_2x_ser_p1_mass_noweight_upx_dirichletx_loy_dirichlety,
    fem_proj_2x_ser_p1_mass_noweight_lox_dirichletx_upx_dirichletx_loy_dirichlety,
    fem_proj_2x_ser_p1_mass_noweight_inx_upy_dirichlety,
    fem_proj_2x_ser_p1_mass_noweight_lox_dirichletx_upy_dirichlety,
    fem_proj_2x_ser_p1_mass_noweight_upx_dirichletx_upy_dirichlety,
    fem_proj_2x_ser_p1_mass_noweight_lox_dirichletx_upx_dirichletx_upy_dirichlety,
    fem_proj_2x_ser_p1_mass_noweight_inx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p1_mass_noweight_lox_dirichletx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p1_mass_noweight_upx_dirichletx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p1_mass_noweight_lox_dirichletx_upx_dirichletx_loy_dirichlety_upy_dirichlety
   },
   {// weighted.
    fem_proj_2x_ser_p1_mass_weighted_inx_iny, fem_proj_2x_ser_p1_mass_weighted_lox_dirichletx_iny,
    fem_proj_2x_ser_p1_mass_weighted_upx_dirichletx_iny,
    fem_proj_2x_ser_p1_mass_weighted_lox_dirichletx_upx_dirichletx_iny,
    fem_proj_2x_ser_p1_mass_weighted_inx_loy_dirichlety,
    fem_proj_2x_ser_p1_mass_weighted_lox_dirichletx_loy_dirichlety,
    fem_proj_2x_ser_p1_mass_weighted_upx_dirichletx_loy_dirichlety,
    fem_proj_2x_ser_p1_mass_weighted_lox_dirichletx_upx_dirichletx_loy_dirichlety,
    fem_proj_2x_ser_p1_mass_weighted_inx_upy_dirichlety,
    fem_proj_2x_ser_p1_mass_weighted_lox_dirichletx_upy_dirichlety,
    fem_proj_2x_ser_p1_mass_weighted_upx_dirichletx_upy_dirichlety,
    fem_proj_2x_ser_p1_mass_weighted_lox_dirichletx_upx_dirichletx_upy_dirichlety,
    fem_proj_2x_ser_p1_mass_weighted_inx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p1_mass_weighted_lox_dirichletx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p1_mass_weighted_upx_dirichletx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p1_mass_weighted_lox_dirichletx_upx_dirichletx_loy_dirichlety_upy_dirichlety
   }
  },
  {// p2.
   {// noweight.
    fem_proj_2x_ser_p2_mass_noweight_inx_iny, fem_proj_2x_ser_p2_mass_noweight_lox_dirichletx_iny,
    fem_proj_2x_ser_p2_mass_noweight_upx_dirichletx_iny,
    fem_proj_2x_ser_p2_mass_noweight_lox_dirichletx_upx_dirichletx_iny,
    fem_proj_2x_ser_p2_mass_noweight_inx_loy_dirichlety,
    fem_proj_2x_ser_p2_mass_noweight_lox_dirichletx_loy_dirichlety,
    fem_proj_2x_ser_p2_mass_noweight_upx_dirichletx_loy_dirichlety,
    fem_proj_2x_ser_p2_mass_noweight_lox_dirichletx_upx_dirichletx_loy_dirichlety,
    fem_proj_2x_ser_p2_mass_noweight_inx_upy_dirichlety,
    fem_proj_2x_ser_p2_mass_noweight_lox_dirichletx_upy_dirichlety,
    fem_proj_2x_ser_p2_mass_noweight_upx_dirichletx_upy_dirichlety,
    fem_proj_2x_ser_p2_mass_noweight_lox_dirichletx_upx_dirichletx_upy_dirichlety,
    fem_proj_2x_ser_p2_mass_noweight_inx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p2_mass_noweight_lox_dirichletx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p2_mass_noweight_upx_dirichletx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p2_mass_noweight_lox_dirichletx_upx_dirichletx_loy_dirichlety_upy_dirichlety
   },
   {// weighted.
    fem_proj_2x_ser_p2_mass_weighted_inx_iny, fem_proj_2x_ser_p2_mass_weighted_lox_dirichletx_iny,
    fem_proj_2x_ser_p2_mass_weighted_upx_dirichletx_iny,
    fem_proj_2x_ser_p2_mass_weighted_lox_dirichletx_upx_dirichletx_iny,
    fem_proj_2x_ser_p2_mass_weighted_inx_loy_dirichlety,
    fem_proj_2x_ser_p2_mass_weighted_lox_dirichletx_loy_dirichlety,
    fem_proj_2x_ser_p2_mass_weighted_upx_dirichletx_loy_dirichlety,
    fem_proj_2x_ser_p2_mass_weighted_lox_dirichletx_upx_dirichletx_loy_dirichlety,
    fem_proj_2x_ser_p2_mass_weighted_inx_upy_dirichlety,
    fem_proj_2x_ser_p2_mass_weighted_lox_dirichletx_upy_dirichlety,
    fem_proj_2x_ser_p2_mass_weighted_upx_dirichletx_upy_dirichlety,
    fem_proj_2x_ser_p2_mass_weighted_lox_dirichletx_upx_dirichletx_upy_dirichlety,
    fem_proj_2x_ser_p2_mass_weighted_inx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p2_mass_weighted_lox_dirichletx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p2_mass_weighted_upx_dirichletx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p2_mass_weighted_lox_dirichletx_upx_dirichletx_loy_dirichlety_upy_dirichlety
   }
  }
};

#endif

// 1D src kernels: polynomial order, weight flag, then active face code.
// Face code uses radix 3, ordered lower x, upper x, lower y, upper y.
GKYL_CU_D static const fem_proj_src_t fem_proj_src_1x[3][2][9] = {
  {{NULL}}, // p0 is unsupported.
  {// p1.
   {// noweight.
    fem_proj_1x_ser_p1_src_noweight_inx, fem_proj_1x_ser_p1_src_noweight_lox_dirichletx,
    fem_proj_1x_ser_p1_src_noweight_lox_dirichletvarx,
    fem_proj_1x_ser_p1_src_noweight_upx_dirichletx,
    fem_proj_1x_ser_p1_src_noweight_lox_dirichletx_upx_dirichletx,
    fem_proj_1x_ser_p1_src_noweight_lox_dirichletvarx_upx_dirichletx,
    fem_proj_1x_ser_p1_src_noweight_upx_dirichletvarx,
    fem_proj_1x_ser_p1_src_noweight_lox_dirichletx_upx_dirichletvarx,
    fem_proj_1x_ser_p1_src_noweight_lox_dirichletvarx_upx_dirichletvarx
   },
   {// weighted.
    fem_proj_1x_ser_p1_src_weighted_inx, fem_proj_1x_ser_p1_src_weighted_lox_dirichletx,
    fem_proj_1x_ser_p1_src_weighted_lox_dirichletvarx,
    fem_proj_1x_ser_p1_src_weighted_upx_dirichletx,
    fem_proj_1x_ser_p1_src_weighted_lox_dirichletx_upx_dirichletx,
    fem_proj_1x_ser_p1_src_weighted_lox_dirichletvarx_upx_dirichletx,
    fem_proj_1x_ser_p1_src_weighted_upx_dirichletvarx,
    fem_proj_1x_ser_p1_src_weighted_lox_dirichletx_upx_dirichletvarx,
    fem_proj_1x_ser_p1_src_weighted_lox_dirichletvarx_upx_dirichletvarx
   }
  },
  {// p2.
   {// noweight.
    fem_proj_1x_ser_p2_src_noweight_inx, fem_proj_1x_ser_p2_src_noweight_lox_dirichletx,
    fem_proj_1x_ser_p2_src_noweight_lox_dirichletvarx,
    fem_proj_1x_ser_p2_src_noweight_upx_dirichletx,
    fem_proj_1x_ser_p2_src_noweight_lox_dirichletx_upx_dirichletx,
    fem_proj_1x_ser_p2_src_noweight_lox_dirichletvarx_upx_dirichletx,
    fem_proj_1x_ser_p2_src_noweight_upx_dirichletvarx,
    fem_proj_1x_ser_p2_src_noweight_lox_dirichletx_upx_dirichletvarx,
    fem_proj_1x_ser_p2_src_noweight_lox_dirichletvarx_upx_dirichletvarx
   },
   {// weighted.
    fem_proj_1x_ser_p2_src_weighted_inx, fem_proj_1x_ser_p2_src_weighted_lox_dirichletx,
    fem_proj_1x_ser_p2_src_weighted_lox_dirichletvarx,
    fem_proj_1x_ser_p2_src_weighted_upx_dirichletx,
    fem_proj_1x_ser_p2_src_weighted_lox_dirichletx_upx_dirichletx,
    fem_proj_1x_ser_p2_src_weighted_lox_dirichletvarx_upx_dirichletx,
    fem_proj_1x_ser_p2_src_weighted_upx_dirichletvarx,
    fem_proj_1x_ser_p2_src_weighted_lox_dirichletx_upx_dirichletvarx,
    fem_proj_1x_ser_p2_src_weighted_lox_dirichletvarx_upx_dirichletvarx
   }
  }
};

// 2D src kernels: polynomial order, weight flag, then active face code.
// Face code uses radix 3, ordered lower x, upper x, lower y, upper y.
GKYL_CU_D static const fem_proj_src_t fem_proj_src_2x[3][2][81] = {
  {{NULL}}, // p0 is unsupported.
  {// p1.
   {// noweight.
    fem_proj_2x_ser_p1_src_noweight_inx_iny,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletx_iny,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletvarx_iny,
    fem_proj_2x_ser_p1_src_noweight_upx_dirichletx_iny,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletx_upx_dirichletx_iny,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletvarx_upx_dirichletx_iny,
    fem_proj_2x_ser_p1_src_noweight_upx_dirichletvarx_iny,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletx_upx_dirichletvarx_iny,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletvarx_upx_dirichletvarx_iny,
    fem_proj_2x_ser_p1_src_noweight_inx_loy_dirichlety,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletx_loy_dirichlety,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletvarx_loy_dirichlety,
    fem_proj_2x_ser_p1_src_noweight_upx_dirichletx_loy_dirichlety,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletx_upx_dirichletx_loy_dirichlety,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletvarx_upx_dirichletx_loy_dirichlety,
    fem_proj_2x_ser_p1_src_noweight_upx_dirichletvarx_loy_dirichlety,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletx_upx_dirichletvarx_loy_dirichlety,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletvarx_upx_dirichletvarx_loy_dirichlety,
    fem_proj_2x_ser_p1_src_noweight_inx_loy_dirichletvary,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletx_loy_dirichletvary,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletvarx_loy_dirichletvary,
    fem_proj_2x_ser_p1_src_noweight_upx_dirichletx_loy_dirichletvary,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletx_upx_dirichletx_loy_dirichletvary,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletvarx_upx_dirichletx_loy_dirichletvary,
    fem_proj_2x_ser_p1_src_noweight_upx_dirichletvarx_loy_dirichletvary,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletx_upx_dirichletvarx_loy_dirichletvary,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletvarx_upx_dirichletvarx_loy_dirichletvary,
    fem_proj_2x_ser_p1_src_noweight_inx_upy_dirichlety,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletx_upy_dirichlety,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletvarx_upy_dirichlety,
    fem_proj_2x_ser_p1_src_noweight_upx_dirichletx_upy_dirichlety,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletx_upx_dirichletx_upy_dirichlety,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletvarx_upx_dirichletx_upy_dirichlety,
    fem_proj_2x_ser_p1_src_noweight_upx_dirichletvarx_upy_dirichlety,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletx_upx_dirichletvarx_upy_dirichlety,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletvarx_upx_dirichletvarx_upy_dirichlety,
    fem_proj_2x_ser_p1_src_noweight_inx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletvarx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p1_src_noweight_upx_dirichletx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletx_upx_dirichletx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletvarx_upx_dirichletx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p1_src_noweight_upx_dirichletvarx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletx_upx_dirichletvarx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletvarx_upx_dirichletvarx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p1_src_noweight_inx_loy_dirichletvary_upy_dirichlety,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletx_loy_dirichletvary_upy_dirichlety,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletvarx_loy_dirichletvary_upy_dirichlety,
    fem_proj_2x_ser_p1_src_noweight_upx_dirichletx_loy_dirichletvary_upy_dirichlety,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletx_upx_dirichletx_loy_dirichletvary_upy_dirichlety,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletvarx_upx_dirichletx_loy_dirichletvary_upy_dirichlety,
    fem_proj_2x_ser_p1_src_noweight_upx_dirichletvarx_loy_dirichletvary_upy_dirichlety,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletx_upx_dirichletvarx_loy_dirichletvary_upy_dirichlety,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletvarx_upx_dirichletvarx_loy_dirichletvary_upy_dirichlety,
    fem_proj_2x_ser_p1_src_noweight_inx_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletx_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletvarx_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_noweight_upx_dirichletx_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletx_upx_dirichletx_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletvarx_upx_dirichletx_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_noweight_upx_dirichletvarx_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletx_upx_dirichletvarx_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletvarx_upx_dirichletvarx_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_noweight_inx_loy_dirichlety_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletx_loy_dirichlety_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletvarx_loy_dirichlety_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_noweight_upx_dirichletx_loy_dirichlety_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletx_upx_dirichletx_loy_dirichlety_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletvarx_upx_dirichletx_loy_dirichlety_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_noweight_upx_dirichletvarx_loy_dirichlety_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletx_upx_dirichletvarx_loy_dirichlety_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletvarx_upx_dirichletvarx_loy_dirichlety_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_noweight_inx_loy_dirichletvary_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletx_loy_dirichletvary_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletvarx_loy_dirichletvary_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_noweight_upx_dirichletx_loy_dirichletvary_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletx_upx_dirichletx_loy_dirichletvary_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletvarx_upx_dirichletx_loy_dirichletvary_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_noweight_upx_dirichletvarx_loy_dirichletvary_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletx_upx_dirichletvarx_loy_dirichletvary_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_noweight_lox_dirichletvarx_upx_dirichletvarx_loy_dirichletvary_upy_dirichletvary
   },
   {// weighted.
    fem_proj_2x_ser_p1_src_weighted_inx_iny,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletx_iny,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletvarx_iny,
    fem_proj_2x_ser_p1_src_weighted_upx_dirichletx_iny,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletx_upx_dirichletx_iny,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletvarx_upx_dirichletx_iny,
    fem_proj_2x_ser_p1_src_weighted_upx_dirichletvarx_iny,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletx_upx_dirichletvarx_iny,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletvarx_upx_dirichletvarx_iny,
    fem_proj_2x_ser_p1_src_weighted_inx_loy_dirichlety,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletx_loy_dirichlety,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletvarx_loy_dirichlety,
    fem_proj_2x_ser_p1_src_weighted_upx_dirichletx_loy_dirichlety,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletx_upx_dirichletx_loy_dirichlety,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletvarx_upx_dirichletx_loy_dirichlety,
    fem_proj_2x_ser_p1_src_weighted_upx_dirichletvarx_loy_dirichlety,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletx_upx_dirichletvarx_loy_dirichlety,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletvarx_upx_dirichletvarx_loy_dirichlety,
    fem_proj_2x_ser_p1_src_weighted_inx_loy_dirichletvary,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletx_loy_dirichletvary,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletvarx_loy_dirichletvary,
    fem_proj_2x_ser_p1_src_weighted_upx_dirichletx_loy_dirichletvary,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletx_upx_dirichletx_loy_dirichletvary,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletvarx_upx_dirichletx_loy_dirichletvary,
    fem_proj_2x_ser_p1_src_weighted_upx_dirichletvarx_loy_dirichletvary,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletx_upx_dirichletvarx_loy_dirichletvary,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletvarx_upx_dirichletvarx_loy_dirichletvary,
    fem_proj_2x_ser_p1_src_weighted_inx_upy_dirichlety,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletx_upy_dirichlety,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletvarx_upy_dirichlety,
    fem_proj_2x_ser_p1_src_weighted_upx_dirichletx_upy_dirichlety,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletx_upx_dirichletx_upy_dirichlety,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletvarx_upx_dirichletx_upy_dirichlety,
    fem_proj_2x_ser_p1_src_weighted_upx_dirichletvarx_upy_dirichlety,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletx_upx_dirichletvarx_upy_dirichlety,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletvarx_upx_dirichletvarx_upy_dirichlety,
    fem_proj_2x_ser_p1_src_weighted_inx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletvarx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p1_src_weighted_upx_dirichletx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletx_upx_dirichletx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletvarx_upx_dirichletx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p1_src_weighted_upx_dirichletvarx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletx_upx_dirichletvarx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletvarx_upx_dirichletvarx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p1_src_weighted_inx_loy_dirichletvary_upy_dirichlety,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletx_loy_dirichletvary_upy_dirichlety,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletvarx_loy_dirichletvary_upy_dirichlety,
    fem_proj_2x_ser_p1_src_weighted_upx_dirichletx_loy_dirichletvary_upy_dirichlety,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletx_upx_dirichletx_loy_dirichletvary_upy_dirichlety,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletvarx_upx_dirichletx_loy_dirichletvary_upy_dirichlety,
    fem_proj_2x_ser_p1_src_weighted_upx_dirichletvarx_loy_dirichletvary_upy_dirichlety,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletx_upx_dirichletvarx_loy_dirichletvary_upy_dirichlety,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletvarx_upx_dirichletvarx_loy_dirichletvary_upy_dirichlety,
    fem_proj_2x_ser_p1_src_weighted_inx_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletx_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletvarx_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_weighted_upx_dirichletx_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletx_upx_dirichletx_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletvarx_upx_dirichletx_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_weighted_upx_dirichletvarx_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletx_upx_dirichletvarx_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletvarx_upx_dirichletvarx_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_weighted_inx_loy_dirichlety_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletx_loy_dirichlety_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletvarx_loy_dirichlety_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_weighted_upx_dirichletx_loy_dirichlety_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletx_upx_dirichletx_loy_dirichlety_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletvarx_upx_dirichletx_loy_dirichlety_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_weighted_upx_dirichletvarx_loy_dirichlety_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletx_upx_dirichletvarx_loy_dirichlety_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletvarx_upx_dirichletvarx_loy_dirichlety_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_weighted_inx_loy_dirichletvary_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletx_loy_dirichletvary_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletvarx_loy_dirichletvary_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_weighted_upx_dirichletx_loy_dirichletvary_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletx_upx_dirichletx_loy_dirichletvary_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletvarx_upx_dirichletx_loy_dirichletvary_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_weighted_upx_dirichletvarx_loy_dirichletvary_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletx_upx_dirichletvarx_loy_dirichletvary_upy_dirichletvary,
    fem_proj_2x_ser_p1_src_weighted_lox_dirichletvarx_upx_dirichletvarx_loy_dirichletvary_upy_dirichletvary
   }
  },
  {// p2.
   {// noweight.
    fem_proj_2x_ser_p2_src_noweight_inx_iny,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletx_iny,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletvarx_iny,
    fem_proj_2x_ser_p2_src_noweight_upx_dirichletx_iny,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletx_upx_dirichletx_iny,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletvarx_upx_dirichletx_iny,
    fem_proj_2x_ser_p2_src_noweight_upx_dirichletvarx_iny,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletx_upx_dirichletvarx_iny,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletvarx_upx_dirichletvarx_iny,
    fem_proj_2x_ser_p2_src_noweight_inx_loy_dirichlety,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletx_loy_dirichlety,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletvarx_loy_dirichlety,
    fem_proj_2x_ser_p2_src_noweight_upx_dirichletx_loy_dirichlety,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletx_upx_dirichletx_loy_dirichlety,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletvarx_upx_dirichletx_loy_dirichlety,
    fem_proj_2x_ser_p2_src_noweight_upx_dirichletvarx_loy_dirichlety,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletx_upx_dirichletvarx_loy_dirichlety,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletvarx_upx_dirichletvarx_loy_dirichlety,
    fem_proj_2x_ser_p2_src_noweight_inx_loy_dirichletvary,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletx_loy_dirichletvary,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletvarx_loy_dirichletvary,
    fem_proj_2x_ser_p2_src_noweight_upx_dirichletx_loy_dirichletvary,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletx_upx_dirichletx_loy_dirichletvary,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletvarx_upx_dirichletx_loy_dirichletvary,
    fem_proj_2x_ser_p2_src_noweight_upx_dirichletvarx_loy_dirichletvary,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletx_upx_dirichletvarx_loy_dirichletvary,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletvarx_upx_dirichletvarx_loy_dirichletvary,
    fem_proj_2x_ser_p2_src_noweight_inx_upy_dirichlety,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletx_upy_dirichlety,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletvarx_upy_dirichlety,
    fem_proj_2x_ser_p2_src_noweight_upx_dirichletx_upy_dirichlety,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletx_upx_dirichletx_upy_dirichlety,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletvarx_upx_dirichletx_upy_dirichlety,
    fem_proj_2x_ser_p2_src_noweight_upx_dirichletvarx_upy_dirichlety,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletx_upx_dirichletvarx_upy_dirichlety,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletvarx_upx_dirichletvarx_upy_dirichlety,
    fem_proj_2x_ser_p2_src_noweight_inx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletvarx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p2_src_noweight_upx_dirichletx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletx_upx_dirichletx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletvarx_upx_dirichletx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p2_src_noweight_upx_dirichletvarx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletx_upx_dirichletvarx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletvarx_upx_dirichletvarx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p2_src_noweight_inx_loy_dirichletvary_upy_dirichlety,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletx_loy_dirichletvary_upy_dirichlety,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletvarx_loy_dirichletvary_upy_dirichlety,
    fem_proj_2x_ser_p2_src_noweight_upx_dirichletx_loy_dirichletvary_upy_dirichlety,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletx_upx_dirichletx_loy_dirichletvary_upy_dirichlety,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletvarx_upx_dirichletx_loy_dirichletvary_upy_dirichlety,
    fem_proj_2x_ser_p2_src_noweight_upx_dirichletvarx_loy_dirichletvary_upy_dirichlety,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletx_upx_dirichletvarx_loy_dirichletvary_upy_dirichlety,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletvarx_upx_dirichletvarx_loy_dirichletvary_upy_dirichlety,
    fem_proj_2x_ser_p2_src_noweight_inx_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletx_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletvarx_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_noweight_upx_dirichletx_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletx_upx_dirichletx_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletvarx_upx_dirichletx_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_noweight_upx_dirichletvarx_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletx_upx_dirichletvarx_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletvarx_upx_dirichletvarx_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_noweight_inx_loy_dirichlety_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletx_loy_dirichlety_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletvarx_loy_dirichlety_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_noweight_upx_dirichletx_loy_dirichlety_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletx_upx_dirichletx_loy_dirichlety_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletvarx_upx_dirichletx_loy_dirichlety_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_noweight_upx_dirichletvarx_loy_dirichlety_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletx_upx_dirichletvarx_loy_dirichlety_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletvarx_upx_dirichletvarx_loy_dirichlety_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_noweight_inx_loy_dirichletvary_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletx_loy_dirichletvary_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletvarx_loy_dirichletvary_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_noweight_upx_dirichletx_loy_dirichletvary_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletx_upx_dirichletx_loy_dirichletvary_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletvarx_upx_dirichletx_loy_dirichletvary_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_noweight_upx_dirichletvarx_loy_dirichletvary_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletx_upx_dirichletvarx_loy_dirichletvary_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_noweight_lox_dirichletvarx_upx_dirichletvarx_loy_dirichletvary_upy_dirichletvary
   },
   {// weighted.
    fem_proj_2x_ser_p2_src_weighted_inx_iny,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletx_iny,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletvarx_iny,
    fem_proj_2x_ser_p2_src_weighted_upx_dirichletx_iny,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletx_upx_dirichletx_iny,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletvarx_upx_dirichletx_iny,
    fem_proj_2x_ser_p2_src_weighted_upx_dirichletvarx_iny,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletx_upx_dirichletvarx_iny,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletvarx_upx_dirichletvarx_iny,
    fem_proj_2x_ser_p2_src_weighted_inx_loy_dirichlety,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletx_loy_dirichlety,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletvarx_loy_dirichlety,
    fem_proj_2x_ser_p2_src_weighted_upx_dirichletx_loy_dirichlety,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletx_upx_dirichletx_loy_dirichlety,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletvarx_upx_dirichletx_loy_dirichlety,
    fem_proj_2x_ser_p2_src_weighted_upx_dirichletvarx_loy_dirichlety,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletx_upx_dirichletvarx_loy_dirichlety,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletvarx_upx_dirichletvarx_loy_dirichlety,
    fem_proj_2x_ser_p2_src_weighted_inx_loy_dirichletvary,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletx_loy_dirichletvary,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletvarx_loy_dirichletvary,
    fem_proj_2x_ser_p2_src_weighted_upx_dirichletx_loy_dirichletvary,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletx_upx_dirichletx_loy_dirichletvary,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletvarx_upx_dirichletx_loy_dirichletvary,
    fem_proj_2x_ser_p2_src_weighted_upx_dirichletvarx_loy_dirichletvary,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletx_upx_dirichletvarx_loy_dirichletvary,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletvarx_upx_dirichletvarx_loy_dirichletvary,
    fem_proj_2x_ser_p2_src_weighted_inx_upy_dirichlety,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletx_upy_dirichlety,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletvarx_upy_dirichlety,
    fem_proj_2x_ser_p2_src_weighted_upx_dirichletx_upy_dirichlety,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletx_upx_dirichletx_upy_dirichlety,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletvarx_upx_dirichletx_upy_dirichlety,
    fem_proj_2x_ser_p2_src_weighted_upx_dirichletvarx_upy_dirichlety,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletx_upx_dirichletvarx_upy_dirichlety,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletvarx_upx_dirichletvarx_upy_dirichlety,
    fem_proj_2x_ser_p2_src_weighted_inx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletvarx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p2_src_weighted_upx_dirichletx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletx_upx_dirichletx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletvarx_upx_dirichletx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p2_src_weighted_upx_dirichletvarx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletx_upx_dirichletvarx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletvarx_upx_dirichletvarx_loy_dirichlety_upy_dirichlety,
    fem_proj_2x_ser_p2_src_weighted_inx_loy_dirichletvary_upy_dirichlety,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletx_loy_dirichletvary_upy_dirichlety,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletvarx_loy_dirichletvary_upy_dirichlety,
    fem_proj_2x_ser_p2_src_weighted_upx_dirichletx_loy_dirichletvary_upy_dirichlety,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletx_upx_dirichletx_loy_dirichletvary_upy_dirichlety,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletvarx_upx_dirichletx_loy_dirichletvary_upy_dirichlety,
    fem_proj_2x_ser_p2_src_weighted_upx_dirichletvarx_loy_dirichletvary_upy_dirichlety,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletx_upx_dirichletvarx_loy_dirichletvary_upy_dirichlety,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletvarx_upx_dirichletvarx_loy_dirichletvary_upy_dirichlety,
    fem_proj_2x_ser_p2_src_weighted_inx_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletx_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletvarx_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_weighted_upx_dirichletx_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletx_upx_dirichletx_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletvarx_upx_dirichletx_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_weighted_upx_dirichletvarx_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletx_upx_dirichletvarx_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletvarx_upx_dirichletvarx_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_weighted_inx_loy_dirichlety_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletx_loy_dirichlety_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletvarx_loy_dirichlety_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_weighted_upx_dirichletx_loy_dirichlety_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletx_upx_dirichletx_loy_dirichlety_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletvarx_upx_dirichletx_loy_dirichlety_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_weighted_upx_dirichletvarx_loy_dirichlety_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletx_upx_dirichletvarx_loy_dirichlety_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletvarx_upx_dirichletvarx_loy_dirichlety_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_weighted_inx_loy_dirichletvary_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletx_loy_dirichletvary_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletvarx_loy_dirichletvary_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_weighted_upx_dirichletx_loy_dirichletvary_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletx_upx_dirichletx_loy_dirichletvary_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletvarx_upx_dirichletx_loy_dirichletvary_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_weighted_upx_dirichletvarx_loy_dirichletvary_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletx_upx_dirichletvarx_loy_dirichletvary_upy_dirichletvary,
    fem_proj_2x_ser_p2_src_weighted_lox_dirichletvarx_upx_dirichletvarx_loy_dirichletvary_upy_dirichletvary
   }
  }
};

// Reconstruction kernels indexed by dimension minus one and polynomial order.
GKYL_CU_D static const fem_proj_transform_t fem_proj_sol_list[2][3] = {
  {NULL, fem_proj_1x_ser_p1_sol, fem_proj_1x_ser_p2_sol},
  {NULL, fem_proj_2x_ser_p1_sol, fem_proj_2x_ser_p2_sol}
};
// Node counts are used only during host initialization.
#ifndef __CUDA_ARCH__
static const fem_proj_num_nodes_t fem_proj_num_nodes_list[2][3] = {
  {NULL, fem_proj_1x_ser_p1_num_nodes, fem_proj_1x_ser_p2_num_nodes},
  {NULL, fem_proj_2x_ser_p1_num_nodes, fem_proj_2x_ser_p2_num_nodes}
};
#endif

// Bias kernels indexed by dimension minus one, polynomial order, and cell face.
#ifndef __CUDA_ARCH__
static const fem_proj_bias_lhs_t fem_proj_bias_lhs_list[2][3][4] = {
  {{NULL},
   {fem_proj_1x_ser_p1_bias_lhs_lox_dirichletx, fem_proj_1x_ser_p1_bias_lhs_upx_dirichletx},
   {fem_proj_1x_ser_p2_bias_lhs_lox_dirichletx, fem_proj_1x_ser_p2_bias_lhs_upx_dirichletx}},
  {{NULL},
   {fem_proj_2x_ser_p1_bias_lhs_lox_dirichletx_iny, fem_proj_2x_ser_p1_bias_lhs_upx_dirichletx_iny,
    fem_proj_2x_ser_p1_bias_lhs_inx_loy_dirichlety, fem_proj_2x_ser_p1_bias_lhs_inx_upy_dirichlety},
   {fem_proj_2x_ser_p2_bias_lhs_lox_dirichletx_iny, fem_proj_2x_ser_p2_bias_lhs_upx_dirichletx_iny,
    fem_proj_2x_ser_p2_bias_lhs_inx_loy_dirichlety, fem_proj_2x_ser_p2_bias_lhs_inx_upy_dirichlety}}
};
#endif

GKYL_CU_D static const fem_proj_bias_src_t fem_proj_bias_src_list[2][3][4] = {
  {{NULL},
   {fem_proj_1x_ser_p1_bias_src_lox_dirichletx, fem_proj_1x_ser_p1_bias_src_upx_dirichletx},
   {fem_proj_1x_ser_p2_bias_src_lox_dirichletx, fem_proj_1x_ser_p2_bias_src_upx_dirichletx}},
  {{NULL},
   {fem_proj_2x_ser_p1_bias_src_lox_dirichletx_iny, fem_proj_2x_ser_p1_bias_src_upx_dirichletx_iny,
    fem_proj_2x_ser_p1_bias_src_inx_loy_dirichlety, fem_proj_2x_ser_p1_bias_src_inx_upy_dirichlety},
   {fem_proj_2x_ser_p2_bias_src_lox_dirichletx_iny, fem_proj_2x_ser_p2_bias_src_upx_dirichletx_iny,
    fem_proj_2x_ser_p2_bias_src_inx_loy_dirichlety, fem_proj_2x_ser_p2_bias_src_inx_upy_dirichlety}}
};

/**
 * Choose kernels once at construction in the CPU or CUDA execution space.
 *
 * @param dim Configuration-space dimension (1 or 2).
 * @param porder Polynomial order (1 or 2).
 * @param weighted_left Whether the mass matrix has a DG weight.
 * @param weighted_right Whether the source has a DG weight.
 * @param periodic_mask Bit d is set when direction d is periodic.
 * @param bcs Lower/upper boundary conditions for each direction.
 * @return Selected kernels; num_nodes, mass and bias_lhs are host-only.
 */
GKYL_CU_D static inline struct fem_proj_kernels
fem_proj_choose_kernels(
  int dim, int porder, bool weighted_left, bool weighted_right, int periodic_mask,
  const struct gkyl_fem_proj_bc *bcs
)
{
  struct fem_proj_kernels kers = {0};

  // Consistency checks.
  assert(dim >= 1 && dim <= 2);
  assert(porder >= 1 && porder <= 2);
  assert(periodic_mask >= 0 && periodic_mask < (1 << dim));

  // Local to global index translation kernel.
  for (int upper = 0; upper < (1 << dim); upper++) {
    if (dim == 1) {
      kers.local_to_global[upper] = fem_proj_l2g_1x[periodic_mask][porder][upper];
    } else if (dim == 2) {
      kers.local_to_global[upper] = fem_proj_l2g_2x[periodic_mask][porder][upper];
    }
  }

#ifndef __CUDA_ARCH__
  // Number of nodes is only needed on the CPU.
  kers.num_nodes = fem_proj_num_nodes_list[dim - 1][porder];
#endif

  // Each direction has interior, lower, upper and both-boundaries locations.
  // Map the active faces to a specialized kernel once, during construction.
  for (int location = 0; location < (1 << (2 * dim)); location++) {
    int mass_code = 0, src_code = 0, src_stride = 1;
    for (int d = 0; d < dim; d++) {
      int position = (location >> (2 * d)) & 3;
      for (int side = 0; side < 2; side++) {
        int type = side ? bcs->up_type[d] : bcs->lo_type[d];
        bool active = position & (1 << side);
        int face_type = active && type == GKYL_FEM_PROJ_DIRICHLET         ? 1 :
                        active && type == GKYL_FEM_PROJ_DIRICHLET_VARYING ? 2 :
                                                                            0;
        mass_code |= (face_type != 0) << (2 * d + side);
        src_code += face_type * src_stride;
        src_stride *= 3;
      }
    }

#ifndef __CUDA_ARCH__
    // Mass kernel only needed on CPU.
    kers.mass[location] = dim == 1 ? fem_proj_mass_1x[porder][weighted_left][mass_code] :
                                     fem_proj_mass_2x[porder][weighted_left][mass_code];
#endif

    kers.src[location] = dim == 1 ? fem_proj_src_1x[porder][weighted_right][src_code] :
                                    fem_proj_src_2x[porder][weighted_right][src_code];
  }
  for (int face = 0; face < 2 * dim; face++) {
#ifndef __CUDA_ARCH__
    kers.bias_lhs[face] = fem_proj_bias_lhs_list[dim - 1][porder][face];
#endif
    kers.bias_src[face] = fem_proj_bias_src_list[dim - 1][porder][face];
  }
  kers.sol = fem_proj_sol_list[dim - 1][porder];
  return kers;
}

// Compute compact global node indices using a cell's parent-range coordinates.
// A one-cell direction is an upper-boundary direction as well as a lower one.
GKYL_CU_DH static inline void
fem_proj_global_indices(
  const struct fem_proj_kernels *kernels, const struct gkyl_range *range, const int *num_cells,
  const int *idx, long *global_idx
)
{
  int relative[GKYL_MAX_CDIM], upper_mask = 0;
  for (int d = 0; d < range->ndim; d++) {
    relative[d] = idx[d] - range->lower[d];
    upper_mask |= (relative[d] == num_cells[d] - 1) << d;
  }
  kernels->local_to_global[upper_mask](num_cells, relative, global_idx);
}

/**
 * Classify a cell as interior (0), lower (1), upper (2), or both (3) per direction.
 * @param range Solve range, including its possibly non-unit lower indices.
 * @param idx Cell indices in the parent range.
 * @return Location index, with x occupying the lowest two bits.
 */
GKYL_CU_DH static inline int
fem_proj_cell_location(const struct gkyl_range *range, const int *idx)
{
  int location = 0;
  for (int d = 0; d < range->ndim; d++) {
    int position = (idx[d] == range->lower[d]) | ((idx[d] == range->upper[d]) << 1);
    location |= position << (2 * d);
  }
  return location;
}

// One side of a bias plane; internal planes have two incident cell ranges.
struct fem_proj_bias_face {
  struct gkyl_range range; // Incident cells, retaining solve-array indexing.
  int face; // Lower/upper face index, 2*dir + side.
  double value; // Constant prescribed potential.
};

struct gkyl_fem_proj {
  struct gkyl_rect_grid grid; // Projection grid.
  struct gkyl_range range; // Solve range with parent-array indexing.
  struct gkyl_basis basis; // Modal DG basis.
  struct gkyl_fem_proj_bc bcs; // Per-direction boundary conditions.
  bool use_gpu; // Execute on the GPU.
  bool varying; // Any boundary requires a Dirichlet DG field.
  bool rhs_ready; // RHS assembled and available for the next solve.
  long num_nodes; // Number of global FEM nodes.
  size_t expect_arr_sz; // Minimum array size required by the solve range.
  int num_cells[GKYL_MAX_CDIM]; // Solve-range cell counts.
  int periodic_mask; // Bit d identifies a periodic direction.
  struct gkyl_array *weight_left; // Fixed mass-matrix weight; NULL means unity.
  struct gkyl_array *weight_right; // Fixed RHS weight; NULL means unity.
  double bc_values[2 * GKYL_MAX_CDIM]; // Constant values ordered lower/upper in each direction.
  struct fem_proj_kernels kernels; // Host kernel function pointers.
  int num_bias_face; // Number of incident cell ranges for fixed bias planes.
  struct fem_proj_bias_face *bias_faces; // Host bias metadata.
  void (*bias_src)(struct gkyl_fem_proj *, double *); // Assign biases after source assembly.
  double *rhs; // Assembled global RHS buffer for CPU solves.
  struct gkyl_superlu_prob *prob; // CPU sparse linear solver.
#ifdef GKYL_HAVE_CUDA
  struct fem_proj_kernels *kernels_cu; // Device kernel function pointers.
  int num_bias_blocks; // Fixed launch size for the largest incident face range.
  struct fem_proj_bias_face *bias_faces_cu; // Device copy of bias metadata.
  double *bc_values_cu; // Device copy of constant boundary values.
  struct gkyl_culinsolver_prob *prob_cu; // GPU sparse linear solver.
#endif
};

#ifdef GKYL_HAVE_CUDA
// Allocate device kernel pointers during construction.
void fem_proj_init_cu(struct gkyl_fem_proj *up);

// Assemble the RHS and enforce boundary values on device.
void fem_proj_set_rhs_cu(
  struct gkyl_fem_proj *up, const struct gkyl_array *rhs, const struct gkyl_array *boundary_values
);

// Assign fixed bias values on device after the volume source kernel.
void fem_proj_bias_src_cu(struct gkyl_fem_proj *up, double *rhs);

// Solve and reconstruct without host transfers.
void fem_proj_solve_cu(struct gkyl_fem_proj *up, struct gkyl_array *out);
#endif
