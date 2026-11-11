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
  up->num_phase_basis = phase_basis->num_basis;
  return up;
}

// Accumulate factor_q (target - f_q) onto the increment of one phase-space cell: F and G relax
// to the Maxwellian jacob f_lte, the first-harmonic blocks to zero.
static inline void
pkpm_bgk_collisions_relax(
  const struct gkyl_pkpm_bgk_collisions *up, const double *factor, double jacob,
  const double *f_lte, const double *fnodal, double *incr
)
{
  for (int b = 0; b < PKPM_BGK_NUM_BLOCKS; ++b) {
    const double *fnodal_b = &fnodal[b * up->num_phase_basis];
    double *incr_b = &incr[b * up->num_phase_basis];
    double target_fac = b < 2 ? jacob : 0.0;
    for (int q = 0; q < up->num_nodes; ++q) {
      for (int j = 0; j < PKPM_BGK_NUM_MODES; ++j) {
        int k = q * PKPM_BGK_NUM_MODES + j;
        incr_b[k] += factor[q] * (target_fac * f_lte[j] - fnodal_b[k]);
      }
    }
  }
}

void
gkyl_pkpm_bgk_collisions_advance(
  const struct gkyl_pkpm_bgk_collisions *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *vel_range, const struct gkyl_range *phase_range,
  const struct gkyl_array *rate, const struct gkyl_array *jacob, const struct gkyl_array *f_lte,
  const struct gkyl_array *fnodal, struct gkyl_array *rhs_nodal, struct gkyl_array *cflrate
)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);
  while (gkyl_range_iter_next(&iter)) {
    long ploc = gkyl_range_idx(phase_range, iter.idx);
    long cloc = gkyl_range_idx(conf_range, iter.idx);
    long vloc = gkyl_range_idx(vel_range, &iter.idx[up->cdim]);

    const double *rate_d = gkyl_array_cfetch(rate, cloc);
    const double *jacob_d = gkyl_array_cfetch(jacob, cloc);
    pkpm_bgk_collisions_relax(
      up, rate_d, jacob_d[0], gkyl_array_cfetch(f_lte, vloc), gkyl_array_cfetch(fnodal, ploc),
      gkyl_array_fetch(rhs_nodal, ploc)
    );

    double rate_max = 0.0;
    for (int q = 0; q < up->num_nodes; ++q) {
      rate_max = fmax(rate_max, rate_d[q]);
    }
    double *cflrate_d = gkyl_array_fetch(cflrate, ploc);
    cflrate_d[0] += rate_max;
  }
}

void
gkyl_pkpm_bgk_collisions_implicit_advance(
  const struct gkyl_pkpm_bgk_collisions *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *vel_range, const struct gkyl_range *phase_range,
  const struct gkyl_array *shape, const struct gkyl_array *jacob, const struct gkyl_array *f_lte,
  const struct gkyl_array *fnodal, struct gkyl_array *incr_nodal
)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);
  while (gkyl_range_iter_next(&iter)) {
    long ploc = gkyl_range_idx(phase_range, iter.idx);
    long cloc = gkyl_range_idx(conf_range, iter.idx);
    long vloc = gkyl_range_idx(vel_range, &iter.idx[up->cdim]);

    const double *shape_d = gkyl_array_cfetch(shape, cloc);
    const double *jacob_d = gkyl_array_cfetch(jacob, cloc);
    double relax[PKPM_BGK_MAX_NODES];
    for (int q = 0; q < up->num_nodes; ++q) {
      relax[q] = 1.0 - shape_d[q];
    }
    pkpm_bgk_collisions_relax(
      up, relax, jacob_d[0], gkyl_array_cfetch(f_lte, vloc), gkyl_array_cfetch(fnodal, ploc),
      gkyl_array_fetch(incr_nodal, ploc)
    );
  }
}

void
gkyl_pkpm_bgk_collisions_release(struct gkyl_pkpm_bgk_collisions *up)
{
  gkyl_free(up);
}
