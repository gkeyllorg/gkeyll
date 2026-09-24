/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_math.h>
#include <gkyl_util.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wave_geom_priv.h>
}

#include <cassert>

// CPU interface to create and track a GPU object
struct gkyl_wave_geom*
gkyl_wave_geom_cu_dev_new(const struct gkyl_rect_grid *grid, struct gkyl_range *range,
  evalf_t mapc2p, void *ctx)
{
  struct gkyl_wave_geom *wg =(struct gkyl_wave_geom*) gkyl_malloc(sizeof(struct gkyl_wave_geom));

  wg->range = *range;

  // Initialize the geometry object on the host side
  struct gkyl_array *geom = gkyl_array_new(GKYL_USER, sizeof(struct gkyl_wave_cell_geom), range->volume);
  double xc[GKYL_MAX_CDIM];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_rect_grid_cell_center(grid, iter.idx, xc);
    struct gkyl_wave_cell_geom *geo =(struct gkyl_wave_cell_geom*) gkyl_array_fetch(geom, gkyl_range_idx(range, iter.idx));
    switch (grid->ndim) {
      case 1:
        calc_geom_1d(grid->dx, xc, mapc2p ? mapc2p : nomapc2p, ctx, geo);
        break;
      case 2:
        calc_geom_2d(grid->dx, xc, mapc2p ? mapc2p : nomapc2p, ctx, geo);
        break;
      case 3:
        calc_geom_3d(grid->dx, xc, mapc2p ? mapc2p : nomapc2p, ctx, geo);
        break;
    };
  }
  // Copy the host-side initialized geometry object to the device
  struct gkyl_array *geom_dev = gkyl_array_cu_dev_new(GKYL_USER, sizeof(struct gkyl_wave_cell_geom), range->volume);
  gkyl_array_copy(geom_dev, geom);
  gkyl_array_release(geom);

  wg->geom = geom_dev->on_dev; // this is so the memcpy below has wv_eqn on_dev

  wg->flags = 0;
  GKYL_SET_CU_ALLOC(wg->flags);
  wg->ref_count = gkyl_ref_count_init(gkyl_wave_geom_free);

  // Initialize the device geometry object
  struct gkyl_wave_geom *wg_cu = (struct gkyl_wave_geom*) gkyl_cu_malloc(sizeof(struct gkyl_wave_geom));
  gkyl_cu_memcpy(wg_cu, wg, sizeof(struct gkyl_wave_geom), GKYL_CU_MEMCPY_H2D);

  wg->on_dev = wg_cu;

  wg->geom = geom_dev; // geometry object should store host pointer
  
  return wg;
}
