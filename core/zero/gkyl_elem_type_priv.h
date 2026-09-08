#pragma once

#include <gkyl_elem_type.h>

// code for array datatype for use in IO
static const uint64_t gkyl_array_data_type[] = { [GKYL_INT] = 0,
  [GKYL_FLOAT] = 1,
  [GKYL_DOUBLE] = 2,
  [GKYL_INT_64] = 3,
  [GKYL_LONG] = 4,
  [GKYL_USER] = 32 };

// mapping of code to datatype: MUST be consistent with the
// gkyl_array_data_type array above
static const int gkyl_array_code_to_data_type[] = { [0] = GKYL_INT,
  [1] = GKYL_FLOAT,
  [2] = GKYL_DOUBLE,
  [3] = GKYL_INT_64,
  [4] = GKYL_LONG,
  [32] = GKYL_USER };

// size in bytes for various data-types
static const size_t gkyl_elem_type_size[] = { [GKYL_INT] = sizeof(int),
  [GKYL_FLOAT] = sizeof(float),
  [GKYL_DOUBLE] = sizeof(double),
  [GKYL_INT_64] = sizeof(int64_t),
  [GKYL_LONG] = sizeof(long),
  [GKYL_USER] = 1 };

static const uint64_t gkyl_file_type_int[] = { [GKYL_FIELD_DATA_FILE] = 1,
  [GKYL_DYNVEC_DATA_FILE] = 2,
  [GKYL_MULTI_RANGE_DATA_FILE] = 3,
  [GKYL_BLOCK_TOPO_DATA_FILE] = 4,
  [GKYL_MULTI_BLOCK_DATA_FILE] = 5 };
