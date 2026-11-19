#pragma once

#include <gkyl_mom_type.h>
#include <gkyl_rect_grid.h>

struct gkyl_mom_calc {
  struct gkyl_rect_grid grid;
  const struct gkyl_mom_type *momt;

  uint32_t flags;
  struct gkyl_mom_calc *on_dev; // pointer to itself or device data
  bool use_gpu; // Whether to run on the gpu.
};
