#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_fv_proj.h>
#include <gkyl_wv_embed_geo.h>

struct gkyl_wv_embed_geo*
gkyl_wv_embed_geo_new(enum gkyl_embed_type type, void *mask_func,
  wv_embed_func_t embed_func, void *ctx)
{
  gkyl_wv_embed_geo *geo = gkyl_malloc(sizeof(gkyl_wv_embed_geo));

  geo->type = type;
  geo->ctx = ctx;
  geo->embed_func = embed_func;
  geo->mask_func = mask_func;

  return geo;
}

void
gkyl_wv_embed_geo_new_mask(struct gkyl_wv_embed_geo *geo, struct gkyl_rect_grid *grid,
  struct gkyl_range *rng, struct gkyl_array *mask)
{
  gkyl_fv_proj *proj = gkyl_fv_proj_new(grid, 1, 1, geo->mask_func,
    geo->ctx);
  gkyl_fv_proj_advance(proj, 0.0, rng, mask);
  gkyl_fv_proj_release(proj);
}

void
gkyl_wv_embed_geo_release(gkyl_wv_embed_geo* geo)
{
  gkyl_free(geo);
}
