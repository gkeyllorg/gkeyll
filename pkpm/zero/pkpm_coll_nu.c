#include <assert.h>
#include <math.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_pkpm_coll_nu.h>
#include <gkyl_pkpm_coll_nu_priv.h>

static void
pkpm_coll_nu_advance_ref(
  const struct gkyl_pkpm_coll_nu *up, const struct gkyl_range *range,
  const struct gkyl_array *nu_ref, const struct gkyl_array *prim_self,
  const struct gkyl_array *prim_other, struct gkyl_array *nu
)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);
    memcpy(gkyl_array_fetch(nu, loc), gkyl_array_cfetch(nu_ref, loc), sizeof(double[up->num_nodes]));
  }
}

static void
pkpm_coll_nu_advance_spitzer(
  const struct gkyl_pkpm_coll_nu *up, const struct gkyl_range *range,
  const struct gkyl_array *nu_ref, const struct gkyl_array *prim_self,
  const struct gkyl_array *prim_other, struct gkyl_array *nu
)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);
    const double *ref_d = gkyl_array_cfetch(nu_ref, loc);
    const double *n_d = gkyl_array_cfetch(prim_other, loc);
    const double *vs_d =
      (const double *)gkyl_array_cfetch(prim_self, loc) + PKPM_COLL_NU_VTSQ * up->num_nodes;
    const double *vr_d =
      (const double *)gkyl_array_cfetch(prim_other, loc) + PKPM_COLL_NU_VTSQ * up->num_nodes;
    double *nu_d = gkyl_array_fetch(nu, loc);
    for (int q = 0; q < up->num_nodes; ++q) {
      double vtsq_sum = fmax(vs_d[q], up->vtsq_min) + fmax(vr_d[q], up->vtsq_min);
      nu_d[q] = (n_d[q] > 0.0 && vtsq_sum > 0.0) ? ref_d[q] * n_d[q] / (vtsq_sum * sqrt(vtsq_sum)) :
                                                   0.0;
    }
  }
}

struct gkyl_pkpm_coll_nu *
gkyl_pkpm_coll_nu_new(
  const struct gkyl_basis *conf_basis, bool spitzer, double vtsq_min, bool use_gpu
)
{
  assert(!use_gpu);
  assert(conf_basis->b_type == GKYL_BASIS_MODAL_SERENDIPITY && conf_basis->poly_order == 1);

  struct gkyl_pkpm_coll_nu *up = gkyl_malloc(sizeof(*up));
  up->num_nodes = conf_basis->num_basis;
  up->vtsq_min = vtsq_min;
  up->advance = spitzer ? pkpm_coll_nu_advance_spitzer : pkpm_coll_nu_advance_ref;
  return up;
}

void
gkyl_pkpm_coll_nu_advance(
  const struct gkyl_pkpm_coll_nu *up, const struct gkyl_range *range,
  const struct gkyl_array *nu_ref, const struct gkyl_array *prim_self,
  const struct gkyl_array *prim_other, struct gkyl_array *nu
)
{
  up->advance(up, range, nu_ref, prim_self, prim_other, nu);
}

void
gkyl_pkpm_coll_nu_release(struct gkyl_pkpm_coll_nu *up)
{
  gkyl_free(up);
}
