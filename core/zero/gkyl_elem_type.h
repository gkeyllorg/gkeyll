#pragma once

#include <stdio.h>
#include <stdint.h>

// Type of element stored in array
enum gkyl_elem_type { GKYL_INT, GKYL_FLOAT, GKYL_DOUBLE, GKYL_INT_64, GKYL_LONG, GKYL_USER };

// Array reduce operators
enum gkyl_array_op { GKYL_MIN, GKYL_MAX, GKYL_SUM, GKYL_SUM_ABS };

// file types for raw IO of Gkeyll data
enum gkyl_file_type {
  GKYL_FIELD_DATA_FILE,
  GKYL_DYNVEC_DATA_FILE,
  GKYL_MULTI_RANGE_DATA_FILE,
  GKYL_BLOCK_TOPO_DATA_FILE,
  GKYL_MULTI_BLOCK_DATA_FILE
};
