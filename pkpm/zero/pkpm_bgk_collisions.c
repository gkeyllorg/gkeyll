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
  assert(phase_basis->num_basis == PKPM_NUM_VEL_MODES * conf_basis->num_basis);

  struct gkyl_pkpm_bgk_collisions *up = gkyl_malloc(sizeof(*up));
  up->cdim = conf_basis->ndim;
  up->num_nodes = conf_basis->num_basis;
  up->num_phase_basis = phase_basis->num_basis;
  up->conf_to_modal = pkpm_bgk_conf_to_modal_list[up->cdim - 1];
  up->phase_to_modal = pkpm_bgk_phase_to_modal_list[up->cdim - 1];
  up->conf_phase_mul = pkpm_bgk_conf_phase_mul_list[up->cdim - 1];
  up->conf_phase_mul_accumulate = pkpm_bgk_conf_phase_mul_accumulate_list[up->cdim - 1];
  return up;
}

// The expansion of factor_q jacob f_lte, the target times a factor given at the nodes, accumulated
// onto one block of a phase-space cell (a p=1 weak product equals the product at the nodes).
static inline void
pkpm_bgk_collisions_target(
  const struct gkyl_pkpm_bgk_collisions *up, const double *factor, double jacob,
  const double *f_lte, double *out
)
{
  double target_nodal[PKPM_BGK_MAX_PHASE_BASIS];
  for (int q = 0; q < up->num_nodes; ++q) {
    for (int j = 0; j < PKPM_NUM_VEL_MODES; ++j) {
      target_nodal[q * PKPM_NUM_VEL_MODES + j] = factor[q] * jacob * f_lte[j];
    }
  }
  up->phase_to_modal(target_nodal, out);
}

void
gkyl_pkpm_bgk_collisions_advance(
  const struct gkyl_pkpm_bgk_collisions *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *vel_range, const struct gkyl_range *phase_range,
  const struct gkyl_array *rate, const struct gkyl_array *jacob, const struct gkyl_array *f_lte,
  const struct gkyl_array *f, struct gkyl_array *rhs, struct gkyl_array *cflrate
)
{
  int num_nodes = up->num_nodes, num_phase_basis = up->num_phase_basis;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);
  while (gkyl_range_iter_next(&iter)) {
    long ploc = gkyl_range_idx(phase_range, iter.idx);
    long cloc = gkyl_range_idx(conf_range, iter.idx);
    long vloc = gkyl_range_idx(vel_range, &iter.idx[up->cdim]);

    const double *rate_d = gkyl_array_cfetch(rate, cloc);
    const double *jacob_d = gkyl_array_cfetch(jacob, cloc);
    const double *f_d = gkyl_array_cfetch(f, ploc);
    double *rhs_d = gkyl_array_fetch(rhs, ploc);

    double rate_modal[PKPM_MAX_CONF_BASIS], target[PKPM_BGK_MAX_PHASE_BASIS];
    up->conf_to_modal(rate_d, rate_modal);
    for (int k = 0; k < num_phase_basis; ++k) {
      target[k] = 0.0;
    }
    pkpm_bgk_collisions_target(up, rate_d, jacob_d[0], gkyl_array_cfetch(f_lte, vloc), target);

    // F and G relax to the Maxwellian, the first-harmonic blocks to zero
    for (int b = 0; b < PKPM_DIST_NUM_BLOCKS; ++b) {
      double *rhs_b = &rhs_d[b * num_phase_basis];
      up->conf_phase_mul_accumulate(-1.0, rate_modal, &f_d[b * num_phase_basis], rhs_b);
      double target_fac = b < PKPM_DIST_M0 ? 1.0 : 0.0;
      for (int k = 0; k < num_phase_basis; ++k) {
        rhs_b[k] += target_fac * target[k];
      }
    }

    double rate_max = 0.0;
    for (int q = 0; q < num_nodes; ++q) {
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
  struct gkyl_array *f
)
{
  int num_nodes = up->num_nodes, num_phase_basis = up->num_phase_basis;

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

    // expansions of S and of (1 - S) J f_t
    double shape_modal[PKPM_MAX_CONF_BASIS], relax[PKPM_MAX_CONF_BASIS];
    double target[PKPM_BGK_MAX_PHASE_BASIS];
    up->conf_to_modal(shape_d, shape_modal);
    for (int q = 0; q < num_nodes; ++q) {
      relax[q] = 1.0 - shape_d[q];
    }
    for (int k = 0; k < num_phase_basis; ++k) {
      target[k] = 0.0;
    }
    pkpm_bgk_collisions_target(up, relax, jacob_d[0], f_lte_d, target);

    // F and G relax to the Maxwellian, the first-harmonic blocks to zero
    for (int b = 0; b < PKPM_DIST_NUM_BLOCKS; ++b) {
      double *f_b = &f_d[b * num_phase_basis];
      double prod[PKPM_BGK_MAX_PHASE_BASIS];
      up->conf_phase_mul(shape_modal, f_b, prod);
      double target_fac = b < PKPM_DIST_M0 ? 1.0 : 0.0;
      for (int k = 0; k < num_phase_basis; ++k) {
        f_b[k] = prod[k] + target_fac * target[k];
      }
    }
  }
}

void
gkyl_pkpm_bgk_collisions_release(struct gkyl_pkpm_bgk_collisions *up)
{
  gkyl_free(up);
}
