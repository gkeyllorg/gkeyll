// Private header: not for direct use.
#pragma once

#include <stdbool.h>

#include <gkyl_dg_gr_maxwell_kernels.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

typedef void (*gr_maxwell_geom_source_t)(const gkyl_dg_gr_maxwell_inp *meq,
  const double *w, const double *dx,
  const double *geom_factor_con_nodal, const double *Jc_fields, double* GKYL_RESTRICT out);

// for use in kernel tables
typedef struct { gr_maxwell_geom_source_t kernels[4]; } gkyl_dg_gr_maxwell_geom_source_kern_list;

struct gkyl_dg_gr_maxwell_geom_source {
  struct gkyl_rect_grid conf_grid; // Configuration-space grid.
  gkyl_dg_gr_maxwell_inp gr_maxwell_data; // Parameters needed by GR-Maxwell kernels.
  int cdim; // Configuration-space dimensions.
  gr_maxwell_geom_source_t geom_source; // Kernel for geometric source contribution.

  uint32_t flags;
  bool use_gpu;
  struct gkyl_dg_gr_maxwell_geom_source *on_dev; // pointer to itself or device data.
};

// Geometric source kernels (Serendipity basis).
GKYL_CU_D
static const gkyl_dg_gr_maxwell_geom_source_kern_list ser_geom_source_kernels[] = {
  // 1x kernels
  { NULL, gr_maxwell_geom_source_1x_ser_p1, gr_maxwell_geom_source_1x_ser_p2, NULL }, // 0
  // 2x kernels
  { NULL, gr_maxwell_geom_source_2x_ser_p1, gr_maxwell_geom_source_2x_ser_p2, NULL }, // 1
  // 3x kernels
  { NULL, gr_maxwell_geom_source_3x_ser_p1, gr_maxwell_geom_source_3x_ser_p2, NULL }, // 2
};

// Geometric source kernels (Tensor basis).
GKYL_CU_D
static const gkyl_dg_gr_maxwell_geom_source_kern_list ten_geom_source_kernels[] = {
  // 1x kernels
  { NULL, gr_maxwell_geom_source_1x_tensor_p1, gr_maxwell_geom_source_1x_tensor_p2, NULL }, // 0
  // 2x kernels
  { NULL, gr_maxwell_geom_source_2x_tensor_p1, gr_maxwell_geom_source_2x_tensor_p2, NULL }, // 1
  // 3x kernels
  { NULL, gr_maxwell_geom_source_3x_tensor_p1, gr_maxwell_geom_source_3x_tensor_p2, NULL }, // 2
};
