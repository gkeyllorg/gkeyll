#include <gkyl_calc_bmag.h>
#include <gkyl_rect_grid.h>
#include <assert.h>

struct gkyl_calc_bmag {
  const struct gkyl_basis* cbasis; //comp basis
  const struct gkyl_basis* pbasis; //physical RZ basis
  const struct gkyl_rect_grid* cgrid; // computational grid
  const struct gkyl_rect_grid* pgrid; // physical RZ grid
  bool use_gpu;
  struct gkyl_bmag_ctx* bmag_ctx;
};

