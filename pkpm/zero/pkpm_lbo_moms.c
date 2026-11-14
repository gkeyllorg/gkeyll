#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_pkpm_lbo_moms.h>
#include <gkyl_pkpm_lbo_moms_priv.h>

// F at the velocity edge of a cell (sign = -1 lower, +1 upper) at each node, from the nodal
// form of its leading block, fnodal[q*3 + j]: the three velocity modes at w = sign.
static inline void
pkpm_lbo_moms_edge_values(int num_nodes, double sign, const double *fnodal, double *edge)
{
  for (int q = 0; q < num_nodes; ++q) {
    const double *mode = &fnodal[q * PKPM_NUM_VEL_MODES];
    edge[q] = 0.7071067811865475 * mode[0] + sign * 1.224744871391589 * mode[1] +
              1.5811388300841898 * mode[2];
  }
}

struct gkyl_pkpm_lbo_moms *
gkyl_pkpm_lbo_moms_new(
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis,
  const double *vel_bounds, bool use_gpu
)
{
  assert(!use_gpu);
  assert(conf_basis->b_type == GKYL_BASIS_MODAL_SERENDIPITY && conf_basis->poly_order == 1);
  assert(phase_basis->b_type == GKYL_BASIS_MODAL_HYBRID);
  assert(phase_basis->ndim == conf_basis->ndim + 1);
  assert(phase_basis->num_basis == PKPM_NUM_VEL_MODES * conf_basis->num_basis);

  struct gkyl_pkpm_lbo_moms *up = gkyl_malloc(sizeof(*up));
  up->cdim = conf_basis->ndim;
  up->num_nodes = conf_basis->num_basis;
  up->vel_bounds[0] = vel_bounds[0];
  up->vel_bounds[1] = vel_bounds[1];
  up->conf_to_nodal = pkpm_lbo_moms_conf_to_nodal_list[up->cdim - 1];
  return up;
}

void
gkyl_pkpm_lbo_moms_advance(
  const struct gkyl_pkpm_lbo_moms *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *moms,
  const struct gkyl_array *fnodal, struct gkyl_array *prim_moms
)
{
  int cdim = up->cdim, num_nodes = up->num_nodes;
  double w_lower = up->vel_bounds[0], w_upper = up->vel_bounds[1];

  int pidx[GKYL_MAX_DIM];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);
    const double *moms_d = gkyl_array_cfetch(moms, loc);
    double *prim_d = gkyl_array_fetch(prim_moms, loc);

    // [M0, M2, M0_G, M1] at the nodes
    double m0[PKPM_MAX_CONF_BASIS], m2[PKPM_MAX_CONF_BASIS];
    double m0_g[PKPM_MAX_CONF_BASIS], m1[PKPM_MAX_CONF_BASIS];
    up->conf_to_nodal(&moms_d[0], m0);
    up->conf_to_nodal(&moms_d[num_nodes], m2);
    up->conf_to_nodal(&moms_d[2 * num_nodes], m0_g);
    up->conf_to_nodal(&moms_d[3 * num_nodes], m1);

    // F at the two edges of the velocity grid at the nodes
    double f_lower[PKPM_MAX_CONF_BASIS], f_upper[PKPM_MAX_CONF_BASIS];
    for (int d = 0; d < cdim; ++d) {
      pidx[d] = iter.idx[d];
    }
    pidx[cdim] = phase_range->lower[cdim];
    pkpm_lbo_moms_edge_values(
      num_nodes, -1.0, gkyl_array_cfetch(fnodal, gkyl_range_idx(phase_range, pidx)), f_lower
    );
    pidx[cdim] = phase_range->upper[cdim];
    pkpm_lbo_moms_edge_values(
      num_nodes, 1.0, gkyl_array_cfetch(fnodal, gkyl_range_idx(phase_range, pidx)), f_upper
    );

    for (int q = 0; q < num_nodes; ++q) {
      double corr_m = f_upper[q] - f_lower[q];
      double corr_e = w_upper * f_upper[q] - w_lower * f_lower[q];
      double energy = m2[q] + 2.0 * m0_g[q];
      double diag = 3.0 * m0[q] - corr_e;
      double det = m0[q] * diag + corr_m * m1[q];
      prim_d[q] = (m1[q] * diag + corr_m * energy) / det;
      prim_d[num_nodes + q] = (m0[q] * energy - m1[q] * m1[q]) / det;
    }
  }
}

void
gkyl_pkpm_lbo_moms_release(struct gkyl_pkpm_lbo_moms *up)
{
  gkyl_free(up);
}
