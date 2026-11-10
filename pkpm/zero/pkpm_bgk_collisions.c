#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_pkpm_bgk_collisions.h>
#include <gkyl_pkpm_bgk_collisions_priv.h>

struct gkyl_pkpm_bgk_collisions *
gkyl_pkpm_bgk_collisions_new(
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, bool use_gpu
)
{
  assert(!use_gpu);
  assert(conf_basis->b_type == GKYL_BASIS_MODAL_SERENDIPITY && conf_basis->poly_order == 1);
  assert(phase_basis->num_basis == PKPM_BGK_NUM_MODES * conf_basis->num_basis);

  struct gkyl_pkpm_bgk_collisions *up = gkyl_malloc(sizeof(*up));
  up->cdim = conf_basis->ndim;
  up->num_nodes = conf_basis->num_basis;
  up->phase_to_nodal = pkpm_bgk_phase_to_nodal_list[up->cdim - 1];
  up->nodal_to_phase = pkpm_bgk_nodal_to_phase_list[up->cdim - 1];
  return up;
}

void
gkyl_pkpm_bgk_collisions_advance(
  const struct gkyl_pkpm_bgk_collisions *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *vel_range, const struct gkyl_range *phase_range,
  const struct gkyl_array *rate, const struct gkyl_array *jacob, const struct gkyl_array *f_lte,
  const struct gkyl_array *fin, struct gkyl_array *rhs, struct gkyl_array *cflrate
)
{
  enum { NUM_NODAL = PKPM_BGK_MAX_NODES * PKPM_BGK_NUM_BLOCKS * PKPM_BGK_NUM_MODES };
  double g[NUM_NODAL], incr[NUM_NODAL];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);
  while (gkyl_range_iter_next(&iter)) {
    long ploc = gkyl_range_idx(phase_range, iter.idx);
    long cloc = gkyl_range_idx(conf_range, iter.idx);
    long vloc = gkyl_range_idx(vel_range, &iter.idx[up->cdim]);

    const double *rate_d = gkyl_array_cfetch(rate, cloc);
    const double *jacob_d = gkyl_array_cfetch(jacob, cloc);
    const double *f_lte_d = gkyl_array_cfetch(f_lte, vloc);

    up->phase_to_nodal(gkyl_array_cfetch(fin, ploc), g);

    double rate_max = 0.0;
    int k = 0;
    for (int q = 0; q < up->num_nodes; ++q) {
      rate_max = fmax(rate_max, rate_d[q]);
      for (int b = 0; b < PKPM_BGK_NUM_BLOCKS; ++b) {
        for (int j = 0; j < PKPM_BGK_NUM_MODES; ++j) {
          // F and G relax to the Maxwellian, the first-harmonic blocks to zero.
          double target = b < 2 ? jacob_d[0] * f_lte_d[j] : 0.0;
          incr[k] = rate_d[q] * (target - g[k]);
          k += 1;
        }
      }
    }

    up->nodal_to_phase(incr, gkyl_array_fetch(rhs, ploc));
    double *cflrate_d = gkyl_array_fetch(cflrate, ploc);
    cflrate_d[0] += rate_max;
  }
}

void
gkyl_pkpm_bgk_collisions_implicit_advance(
  const struct gkyl_pkpm_bgk_collisions *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *vel_range, const struct gkyl_range *phase_range,
  const struct gkyl_array *shape, const struct gkyl_array *jacob, const struct gkyl_array *f_lte,
  struct gkyl_array *f
)
{
  enum { NUM_NODAL = PKPM_BGK_MAX_NODES * PKPM_BGK_NUM_BLOCKS * PKPM_BGK_NUM_MODES };
  double g[NUM_NODAL], incr[NUM_NODAL];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);
  while (gkyl_range_iter_next(&iter)) {
    long ploc = gkyl_range_idx(phase_range, iter.idx);
    long cloc = gkyl_range_idx(conf_range, iter.idx);
    long vloc = gkyl_range_idx(vel_range, &iter.idx[up->cdim]);

    const double *shape_d = gkyl_array_cfetch(shape, cloc);
    const double *jacob_d = gkyl_array_cfetch(jacob, cloc);
    const double *f_lte_d = gkyl_array_cfetch(f_lte, vloc);
    double *f_d = gkyl_array_fetch(f, ploc);

    up->phase_to_nodal(f_d, g);

    int k = 0;
    for (int q = 0; q < up->num_nodes; ++q) {
      double relax = 1.0 - shape_d[q];
      for (int b = 0; b < PKPM_BGK_NUM_BLOCKS; ++b) {
        for (int j = 0; j < PKPM_BGK_NUM_MODES; ++j) {
          double target = b < 2 ? jacob_d[0] * f_lte_d[j] : 0.0;
          incr[k] = relax * (target - g[k]);
          k += 1;
        }
      }
    }

    up->nodal_to_phase(incr, f_d);
  }
}

void
gkyl_pkpm_bgk_collisions_release(struct gkyl_pkpm_bgk_collisions *up)
{
  gkyl_free(up);
}
