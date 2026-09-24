#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_dg_geom.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_math.h>
#include <gkyl_util.h>

static bool
dg_geom_is_cu_dev(const struct gkyl_dg_geom* dgg)
{
  return GKYL_IS_CU_ALLOC(dgg->flags);
}

void
dg_geom_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_geom *dgg = container_of(ref, struct gkyl_dg_geom, ref_count);

  for (int d=0; d<dgg->range.ndim; ++d)
    gkyl_array_release(dgg->surf_geom[d]);

  gkyl_array_release(dgg->vol_geom);

  if (!dg_geom_is_cu_dev(dgg)) {
    gkyl_free(dgg->vol_weights);
    gkyl_free(dgg->vol_ords);

    gkyl_free(dgg->surf_weights);
    gkyl_free(dgg->surf_ords);
  }
  
  if (dg_geom_is_cu_dev(dgg)) 
    gkyl_cu_free(dgg->on_dev); 

  gkyl_free(dgg);
}


struct gkyl_dg_geom *
gkyl_dg_geom_new(const struct gkyl_dg_geom_inp *inp)
{
  struct gkyl_dg_geom *dgg = gkyl_malloc(sizeof *dgg);

  dgg->range = *inp->range;

  int ndim = dgg->range.ndim;
  int shape[GKYL_MAX_CDIM];
  for (int d=0; d<ndim; ++d) shape[d] = inp->nquad;

  // NOTE: surfaces are ndim-1 objects
  gkyl_range_init_from_shape(&dgg->surf_quad_range, ndim-1, shape);
  gkyl_range_init_from_shape(&dgg->vol_quad_range, ndim, shape);

  for (int d=0; d<ndim; ++d)
    dgg->surf_geom[d] = gkyl_array_new(GKYL_USER,
      sizeof(struct gkyl_dg_surf_geom[dgg->surf_quad_range.volume]), dgg->range.volume);

  dgg->vol_geom = gkyl_array_new(GKYL_USER,
    sizeof(struct gkyl_dg_vol_geom[dgg->vol_quad_range.volume]), dgg->range.volume);
  
  // compute surface and volume quadrature weights & ordinates
  long nsq = dgg->surf_quad_range.volume;
  dgg->surf_weights = gkyl_malloc(sizeof(double)*nsq);
  dgg->surf_ords = gkyl_malloc(sizeof(double)*nsq*(ndim-1));
  gkyl_ndim_ordinates_weights(ndim-1, dgg->surf_ords, dgg->surf_weights, inp->nquad);

  long nvq = dgg->vol_quad_range.volume;
  dgg->vol_weights = gkyl_malloc(sizeof(double)*nvq);
  dgg->vol_ords = gkyl_malloc(sizeof(double)*nvq*ndim);
  gkyl_ndim_ordinates_weights(ndim, dgg->vol_ords, dgg->vol_weights, inp->nquad);
  
  dgg->flags = 0;
  GKYL_CLEAR_CU_ALLOC(dgg->flags);
  dgg->ref_count = gkyl_ref_count_init(dg_geom_free);
  dgg->on_dev = dgg; // CPU eqn obj points to itself
  
  return dgg;
}

struct gkyl_dg_geom *
gkyl_dg_geom_new_from_host(const struct gkyl_dg_geom_inp *inp, struct gkyl_dg_geom *up_host, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    return gkyl_dg_geom_cu_dev_new_from_host(inp, up_host);
  } 
#endif 
  return up_host;
}


struct gkyl_dg_geom*
gkyl_dg_geom_acquire(const struct gkyl_dg_geom* dgg)
{
  gkyl_ref_count_inc(&dgg->ref_count);
  return (struct gkyl_dg_geom*) dgg;
}

void
gkyl_dg_geom_write(const struct gkyl_dg_geom* dgg, const char *fname)
{
  
}

void
gkyl_dg_geom_release(const struct gkyl_dg_geom *dgg)
{
  gkyl_ref_count_dec(&dgg->ref_count);
}
