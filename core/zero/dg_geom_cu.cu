extern "C" {
#include <assert.h>
#include <math.h>
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_dg_geom.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_math.h>
#include <gkyl_util.h>
}

struct gkyl_dg_geom *
gkyl_dg_geom_cu_dev_new_from_host(const struct gkyl_dg_geom_inp *inp, struct gkyl_dg_geom *up_host)
{
  struct gkyl_dg_geom *dgg = (struct gkyl_dg_geom *) gkyl_malloc(sizeof *dgg);

  dgg->range = *inp->range;

  int ndim = dgg->range.ndim;
  int shape[GKYL_MAX_CDIM];
  for (int d=0; d<ndim; ++d) shape[d] = inp->nquad;

  // NOTE: surfaces are ndim-1 objects
  gkyl_range_init_from_shape(&dgg->surf_quad_range, ndim-1, shape);
  gkyl_range_init_from_shape(&dgg->vol_quad_range, ndim, shape);

  // Copy the host-side initialized geometry object to the device
  struct gkyl_array *vol_geom_dev = gkyl_array_cu_dev_new(GKYL_USER,
    sizeof(struct gkyl_dg_vol_geom[dgg->vol_quad_range.volume]), dgg->range.volume);

  struct gkyl_array *surf_geom_dev[ndim];
  for (int dir=0; dir<ndim; ++dir) {
    surf_geom_dev[dir] = gkyl_array_cu_dev_new(GKYL_USER,
      sizeof(struct gkyl_dg_surf_geom[dgg->surf_quad_range.volume]), dgg->range.volume);
  }

  gkyl_array_copy(vol_geom_dev, up_host->vol_geom);

  for (int dir=0; dir<ndim; ++dir)
    gkyl_array_copy(surf_geom_dev[dir], up_host->surf_geom[dir]);

  // this is for the memcpy below
  dgg->vol_geom = vol_geom_dev->on_dev;

  for (int dir=0; dir<ndim; ++dir)
     dgg->surf_geom[dir] = surf_geom_dev[dir]->on_dev;

  dgg->flags = 0;
  GKYL_SET_CU_ALLOC(dgg->flags);
  dgg->ref_count = gkyl_ref_count_init(dg_geom_free);

  // Initialize the device geometry object
  struct gkyl_dg_geom *dgg_cu = (struct gkyl_dg_geom*) gkyl_cu_malloc(sizeof(struct gkyl_dg_geom));
  gkyl_cu_memcpy(dgg_cu, dgg, sizeof(struct gkyl_dg_geom), GKYL_CU_MEMCPY_H2D);
  dgg->on_dev = dgg_cu;

  // geometry object should store host pointer
  dgg->vol_geom = vol_geom_dev;
  for (int dir=0; dir<ndim; ++dir)
    dgg->surf_geom[dir] = surf_geom_dev[dir];

  
  return dgg;
}
