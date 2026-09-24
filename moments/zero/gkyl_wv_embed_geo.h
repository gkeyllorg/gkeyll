#pragma once

#include <gkyl_array.h>
#include <gkyl_evalf_def.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wv_eqn.h>

// BC types in this updater.
enum gkyl_embed_type {  
  GKYL_EMBED_ABSORB = 0, 
  GKYL_EMBED_REFLECT,
  GKYL_EMBED_COPY_B,
  GKYL_EMBED_PEC,
  GKYL_EMBED_FUNC,
};

// Object type
typedef struct gkyl_wv_embed_geo gkyl_wv_embed_geo;

struct gkyl_wv_embed_geo {
  void *ctx;
  void *mask_func;
  enum gkyl_embed_type type;
  wv_embed_func_t embed_func;
};

/**
 * Create new updater to apply set of boundary conditions.
 *
 * @param mask_func
 * @param type
 * @param ctx Context to pass to bcfunc.
 * @return New updater pointer.
 */
gkyl_wv_embed_geo* gkyl_wv_embed_geo_new(enum gkyl_embed_type type, void *mask_func,
  wv_embed_func_t embed_func, void *ctx);

void gkyl_wv_embed_geo_new_mask(struct gkyl_wv_embed_geo *geo,
  struct gkyl_rect_grid *grid, struct gkyl_range *rng, struct gkyl_array *mask);

/**
 * Delete structure.
 *
 * @param embed_geo Structure to delete.
 */
void gkyl_wv_embed_geo_release(gkyl_wv_embed_geo *embed_geo);
