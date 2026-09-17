// Private header for dg_eval_at_coord_proj updater. Not for direct use.
#pragma once

#include <assert.h>
#include <stdbool.h>

#include <gkyl_alloc.h>
#include <gkyl_dg_eval_at_coord_proj.h>
#include <gkyl_dg_eval_at_coord_proj_kernels.h>
#include <gkyl_util.h>

// Function pointer type for evaluation-projection kernels.
typedef void (*eval_at_coord_t)(const double *coords, const double *fdo, double *ftar);

// Function pointer type for kernels returning target basis.
typedef void (*basis_tar_t)(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis);

// For use in kernel tables.
typedef struct { eval_at_coord_t kernels[3]; } eval_at_coord_kern_list;
typedef struct { basis_tar_t kernels[3]; } basis_tar_kern_list;

// Struct with double array for passing coordinates to CUDA kernels.
typedef struct { double c[GKYL_MAX_DIM]; } dg_evproj_struct_double_t;

// Struct with int array for passing cell index to CUDA kernels.
typedef struct { int c[GKYL_MAX_DIM]; } dg_evproj_struct_int_t;

// Struct with bool array for passing is_dir_eval to CUDA kernels.
typedef struct { bool c[GKYL_MAX_DIM]; } dg_evproj_struct_bool_t;

// Serendipity kernels.
// Indexed as [ndim_do-1][dir_bitmask-1].kernels[poly_order-1]
//   ndim_do = 1 -> index 0   (1  combos)
//   ndim_do = 2 -> index 1   (3  combos)
//   ndim_do = 3 -> index 2   (7  combos)
//   ndim_do = 4 -> index 3   (15 combos)
//   ndim_do = 5 -> index 4   (31 combos)
//   ndim_do = 6 -> index 5   (63 combos)
GKYL_CU_D
static const eval_at_coord_kern_list ser_eval_at_coord_list[6][63] = {
  // ndim_do = 1.
  {
    { gkyl_dg_eval_at_coord_proj_1x_ser_p1_eval_dirs_0, gkyl_dg_eval_at_coord_proj_1x_ser_p2_eval_dirs_0, NULL },
  },
  // ndim_do = 2.
  {
    { gkyl_dg_eval_at_coord_proj_2x_ser_p1_eval_dirs_0 , gkyl_dg_eval_at_coord_proj_2x_ser_p2_eval_dirs_0 , NULL },
    { gkyl_dg_eval_at_coord_proj_2x_ser_p1_eval_dirs_1 , gkyl_dg_eval_at_coord_proj_2x_ser_p2_eval_dirs_1 , NULL },
    { gkyl_dg_eval_at_coord_proj_2x_ser_p1_eval_dirs_01, gkyl_dg_eval_at_coord_proj_2x_ser_p2_eval_dirs_01, NULL },
  },
  // ndim_do = 3.
  {
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_0  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_1  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_01 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_2  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_02 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_12 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_012, NULL, NULL },
  },
  // ndim_do = 4.
  {
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_0   , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_0   , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_1   , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_1   , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_01  , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_01  , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_2   , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_2   , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_02  , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_02  , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_12  , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_12  , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_012 , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_012 , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_3   , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_3   , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_03  , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_03  , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_13  , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_13  , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_013 , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_013 , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_23  , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_23  , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_023 , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_023 , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_123 , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_123 , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_0123, gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_0123, NULL },
  },
  // ndim_do = 5.
  {
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_0    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_1    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_01   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_2    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_02   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_12   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_012  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_3    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_03   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_13   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_013  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_23   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_023  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_123  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_0123 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_4    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_04   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_14   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_014  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_24   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_024  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_124  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_0124 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_34   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_034  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_134  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_0134 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_234  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_0234 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_1234 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_01234, NULL, NULL },
  },
  // ndim_do = 6.
  {
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0     , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_1     , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_01    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_2     , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_02    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_12    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_012   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_3     , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_03    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_13    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_013   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_23    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_023   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_123   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0123  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_4     , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_04    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_14    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_014   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_24    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_024   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_124   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0124  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_34    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_034   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_134   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0134  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_234   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0234  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_1234  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_01234 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_5     , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_05    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_15    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_015   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_25    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_025   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_125   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0125  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_35    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_035   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_135   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0135  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_235   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0235  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_1235  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_01235 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_45    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_045   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_145   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0145  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_245   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0245  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_1245  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_01245 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_345   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0345  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_1345  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_01345 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_2345  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_02345 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_12345 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_012345, NULL, NULL },
  },
};

GKYL_CU_D
static const basis_tar_kern_list ser_basis_tar_list[6][63] = {
  // ndim_do = 1.
  {
    { gkyl_dg_eval_at_coord_proj_1x_ser_p1_eval_dirs_0_target_basis, gkyl_dg_eval_at_coord_proj_1x_ser_p2_eval_dirs_0_target_basis, NULL },
  },
  // ndim_do = 2.
  {
    { gkyl_dg_eval_at_coord_proj_2x_ser_p1_eval_dirs_0_target_basis , gkyl_dg_eval_at_coord_proj_2x_ser_p2_eval_dirs_0_target_basis , NULL },
    { gkyl_dg_eval_at_coord_proj_2x_ser_p1_eval_dirs_1_target_basis , gkyl_dg_eval_at_coord_proj_2x_ser_p2_eval_dirs_1_target_basis , NULL },
    { gkyl_dg_eval_at_coord_proj_2x_ser_p1_eval_dirs_01_target_basis, gkyl_dg_eval_at_coord_proj_2x_ser_p2_eval_dirs_01_target_basis, NULL },
  },
  // ndim_do = 3.
  {
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_0_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_1_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_01_target_basis , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_2_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_02_target_basis , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_12_target_basis , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_012_target_basis, NULL, NULL },
  },
  // ndim_do = 4.
  {
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_0_target_basis   , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_0_target_basis   , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_1_target_basis   , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_1_target_basis   , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_01_target_basis  , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_01_target_basis  , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_2_target_basis   , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_2_target_basis   , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_02_target_basis  , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_02_target_basis  , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_12_target_basis  , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_12_target_basis  , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_012_target_basis , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_012_target_basis , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_3_target_basis   , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_3_target_basis   , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_03_target_basis  , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_03_target_basis  , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_13_target_basis  , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_13_target_basis  , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_013_target_basis , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_013_target_basis , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_23_target_basis  , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_23_target_basis  , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_023_target_basis , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_023_target_basis , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_123_target_basis , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_123_target_basis , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_0123_target_basis, gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_0123_target_basis, NULL },
  },
  // ndim_do = 5.
  {
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_0_target_basis    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_1_target_basis    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_01_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_2_target_basis    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_02_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_12_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_012_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_3_target_basis    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_03_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_13_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_013_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_23_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_023_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_123_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_0123_target_basis , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_4_target_basis    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_04_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_14_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_014_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_24_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_024_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_124_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_0124_target_basis , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_34_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_034_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_134_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_0134_target_basis , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_234_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_0234_target_basis , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_1234_target_basis , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_01234_target_basis, NULL, NULL },
  },
  // ndim_do = 6.
  {
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0_target_basis     , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_1_target_basis     , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_01_target_basis    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_2_target_basis     , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_02_target_basis    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_12_target_basis    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_012_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_3_target_basis     , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_03_target_basis    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_13_target_basis    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_013_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_23_target_basis    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_023_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_123_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0123_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_4_target_basis     , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_04_target_basis    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_14_target_basis    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_014_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_24_target_basis    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_024_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_124_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0124_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_34_target_basis    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_034_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_134_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0134_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_234_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0234_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_1234_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_01234_target_basis , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_5_target_basis     , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_05_target_basis    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_15_target_basis    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_015_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_25_target_basis    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_025_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_125_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0125_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_35_target_basis    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_035_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_135_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0135_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_235_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0235_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_1235_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_01235_target_basis , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_45_target_basis    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_045_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_145_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0145_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_245_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0245_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_1245_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_01245_target_basis , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_345_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0345_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_1345_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_01345_target_basis , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_2345_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_02345_target_basis , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_12345_target_basis , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_012345_target_basis, NULL, NULL },
  },
};

// Tensor kernels.
GKYL_CU_D
static const eval_at_coord_kern_list ten_eval_at_coord_list[3][7] = {
  // ndim_do = 1.
  {
    { gkyl_dg_eval_at_coord_proj_1x_ser_p1_eval_dirs_0, NULL, NULL },
    { NULL, NULL, NULL },
    { NULL, NULL, NULL },
    { NULL, NULL, NULL },
    { NULL, NULL, NULL },
    { NULL, NULL, NULL },
    { NULL, NULL, NULL },
  },
  // ndim_do = 2.
  {
    { gkyl_dg_eval_at_coord_proj_2x_ser_p1_eval_dirs_0 , gkyl_dg_eval_at_coord_proj_2x_tensor_p2_eval_dirs_0 , NULL },
    { gkyl_dg_eval_at_coord_proj_2x_ser_p1_eval_dirs_1 , gkyl_dg_eval_at_coord_proj_2x_tensor_p2_eval_dirs_1 , NULL },
    { gkyl_dg_eval_at_coord_proj_2x_ser_p1_eval_dirs_01, gkyl_dg_eval_at_coord_proj_2x_tensor_p2_eval_dirs_01, NULL },
    { NULL, NULL, NULL },
    { NULL, NULL, NULL },
    { NULL, NULL, NULL },
    { NULL, NULL, NULL },
  },
  // ndim_do = 3.
  {
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_0  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_1  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_01 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_2  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_02 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_12 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_012, NULL, NULL },
  },
};

GKYL_CU_D
static const basis_tar_kern_list ten_basis_tar_list[3][7] = {
  // ndim_do = 1.
  {
    { gkyl_dg_eval_at_coord_proj_1x_ser_p1_eval_dirs_0_target_basis, NULL, NULL },
    { NULL, NULL, NULL },
    { NULL, NULL, NULL },
    { NULL, NULL, NULL },
    { NULL, NULL, NULL },
    { NULL, NULL, NULL },
    { NULL, NULL, NULL },
  },
  // ndim_do = 2.
  {
    { gkyl_dg_eval_at_coord_proj_2x_ser_p1_eval_dirs_0_target_basis , gkyl_dg_eval_at_coord_proj_2x_tensor_p2_eval_dirs_0_target_basis , NULL },
    { gkyl_dg_eval_at_coord_proj_2x_ser_p1_eval_dirs_1_target_basis , gkyl_dg_eval_at_coord_proj_2x_tensor_p2_eval_dirs_1_target_basis , NULL },
    { gkyl_dg_eval_at_coord_proj_2x_ser_p1_eval_dirs_01_target_basis, gkyl_dg_eval_at_coord_proj_2x_tensor_p2_eval_dirs_01_target_basis, NULL },
    { NULL, NULL, NULL },
    { NULL, NULL, NULL },
    { NULL, NULL, NULL },
    { NULL, NULL, NULL },
  },
  // ndim_do = 3.
  {
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_0_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_1_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_01_target_basis , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_2_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_02_target_basis , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_12_target_basis , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_012_target_basis, NULL, NULL },
  },
};

// GK-hybrid kernels.
GKYL_CU_D
static const eval_at_coord_kern_list gkhyb_eval_at_coord_list[4][31] = {
  // 1x1v.
  {
    { gkyl_dg_eval_at_coord_proj_1x1v_gkhyb_p1_eval_dirs_0 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_1x1v_gkhyb_p1_eval_dirs_1 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_1x1v_gkhyb_p1_eval_dirs_01, NULL, NULL },
  },
  // 1x2v.
  {
    { gkyl_dg_eval_at_coord_proj_1x2v_gkhyb_p1_eval_dirs_0  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_1x2v_gkhyb_p1_eval_dirs_1  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_1x2v_gkhyb_p1_eval_dirs_01 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_1x2v_gkhyb_p1_eval_dirs_2  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_1x2v_gkhyb_p1_eval_dirs_02 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_1x2v_gkhyb_p1_eval_dirs_12 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_1x2v_gkhyb_p1_eval_dirs_012, NULL, NULL },
  },
  // 2x2v.
  {
    { gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_0   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_1   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_01  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_2   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_02  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_12  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_012 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_3   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_03  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_13  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_013 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_23  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_023 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_123 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_0123, NULL, NULL },
  },
  // 3x2v.
  {
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_0    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_1    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_01   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_2    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_02   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_12   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_012  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_3    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_03   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_13   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_013  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_23   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_023  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_123  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_0123 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_4    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_04   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_14   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_014  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_24   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_024  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_124  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_0124 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_34   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_034  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_134  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_0134 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_234  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_0234 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_1234 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_01234, NULL, NULL },
  },
};

GKYL_CU_D
static const basis_tar_kern_list gkhyb_basis_tar_list[4][31] = {
  // 1x1v.
  {
    { gkyl_dg_eval_at_coord_proj_1x1v_gkhyb_p1_eval_dirs_0_target_basis , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_1x1v_gkhyb_p1_eval_dirs_1_target_basis , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_1x1v_gkhyb_p1_eval_dirs_01_target_basis, NULL, NULL },
  },
  // 1x2v.
  {
    { gkyl_dg_eval_at_coord_proj_1x2v_gkhyb_p1_eval_dirs_0_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_1x2v_gkhyb_p1_eval_dirs_1_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_1x2v_gkhyb_p1_eval_dirs_01_target_basis , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_1x2v_gkhyb_p1_eval_dirs_2_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_1x2v_gkhyb_p1_eval_dirs_02_target_basis , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_1x2v_gkhyb_p1_eval_dirs_12_target_basis , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_1x2v_gkhyb_p1_eval_dirs_012_target_basis, NULL, NULL },
  },
  // 2x2v.
  {
    { gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_0_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_1_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_01_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_2_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_02_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_12_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_012_target_basis , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_3_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_03_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_13_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_013_target_basis , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_23_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_023_target_basis , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_123_target_basis , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_0123_target_basis, NULL, NULL },
  },
  // 3x2v.
  {
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_0_target_basis    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_1_target_basis    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_01_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_2_target_basis    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_02_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_12_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_012_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_3_target_basis    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_03_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_13_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_013_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_23_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_023_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_123_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_0123_target_basis , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_4_target_basis    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_04_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_14_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_014_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_24_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_024_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_124_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_0124_target_basis , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_34_target_basis   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_034_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_134_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_0134_target_basis , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_234_target_basis  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_0234_target_basis , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_1234_target_basis , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_01234_target_basis, NULL, NULL },
  },
};

GKYL_CU_DH static void
eval_at_coord_get_idx_do(const bool *is_eval, int ndim_do,
  const int *idx_tar, const int *cell_idx, int *idx_do)
{
  // Translate a target index to the corresponding donor index.
  int c = 0;
  for (int d=0; d<ndim_do; d++) {
    if (!is_eval[d])
      idx_do[d] = idx_tar[c++];
    else
      idx_do[d] = cell_idx[d];
  }
}

struct dg_ev_proj_kernels {
  eval_at_coord_t ev_ker; // Projection kernel.
  basis_tar_t basis_ker; // Target basis kernel.
};

// Primary struct for this updater.
struct gkyl_dg_eval_at_coord_proj {
  int ndim_do;                 // Number of dimensions in donor field.
  int ndim_tar;                // Number of dimensions in target field.
  int num_basis_do;            // Number of monomials in donor basis.
  int num_basis_tar;           // Number of monomials in target basis.
  int num_eval_dirs;           // Number of directions being evaluated.
  int eval_dirs[GKYL_MAX_DIM]; // Which directions are evaluated.
  bool is_eval[GKYL_MAX_DIM];  // Is direction evaluated.
  bool use_gpu;                // Whether to run on GPU.
  struct dg_ev_proj_kernels *kers; // Projection kernel.
};

#ifdef GKYL_HAVE_CUDA
// Host-side wrapper for the GPU implementation of gkyl_dg_eval_at_coord_proj_advance.
void gkyl_dg_eval_at_coord_proj_advance_cu(struct gkyl_dg_eval_at_coord_proj *up, const double *eval_coords,
  const struct gkyl_rect_grid *grid, const bool *pick_lower, const int *known_index,
  const struct gkyl_range *rng_do, const struct gkyl_range *rng_tar,
  const struct gkyl_array *fdo, struct gkyl_array *ftar);

struct dg_ev_proj_kernels* dg_eval_at_coord_choose_ker_cu(int cdim, int ndim,
  const struct gkyl_basis *basis, int num_eval_dirs, const int *eval_dirs);
#endif

GKYL_CU_DH static int
eval_dirs_to_mask(int num_eval_dirs, const int *eval_dirs)
{
  // Encode eval_dirs as a bitmask: bit d is set iff direction d is evaluated.
  // Used to index the dispatch table as [ndim_do-1][mask-1].
  int mask = 0;
  for (int i=0; i<num_eval_dirs; i++)
    mask |= (1 << eval_dirs[i]);
  return mask;
}

GKYL_CU_D static struct dg_ev_proj_kernels*
dg_eval_at_coord_choose_ker(bool use_gpu, int cdim, int ndim, const struct gkyl_basis *basis, int num_eval_dirs,
  const int *eval_dirs)
{
  // Choose the projection kernel.

#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    return dg_eval_at_coord_choose_ker_cu(cdim, ndim, basis, num_eval_dirs, eval_dirs);
#endif

  struct dg_ev_proj_kernels *kers = (struct dg_ev_proj_kernels *) gkyl_calloc(1, sizeof(struct dg_ev_proj_kernels));

  int dir_mask = eval_dirs_to_mask(num_eval_dirs, eval_dirs);

  int poly_order = basis->poly_order;

  assert(ndim >= 1 && ndim <= 6);
  assert(dir_mask >= 1 && dir_mask < (1 << ndim));
  assert(poly_order >= 1 && poly_order <= 3);

  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kers->ev_ker = ser_eval_at_coord_list[ndim-1][dir_mask-1].kernels[poly_order-1];
      kers->basis_ker = ser_basis_tar_list[ndim-1][dir_mask-1].kernels[poly_order-1];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      kers->ev_ker = ten_eval_at_coord_list[ndim-1][dir_mask-1].kernels[poly_order-1];
      kers->basis_ker = ten_basis_tar_list[ndim-1][dir_mask-1].kernels[poly_order-1];
      break;
    case GKYL_BASIS_MODAL_GKHYBRID:
      kers->ev_ker = gkhyb_eval_at_coord_list[ndim-2][dir_mask-1].kernels[poly_order-1];
      kers->basis_ker = gkhyb_basis_tar_list[ndim-2][dir_mask-1].kernels[poly_order-1];
      break;
    default:
      assert(false);
      break;
  }
  assert(kers->ev_ker);

  return kers;
}
