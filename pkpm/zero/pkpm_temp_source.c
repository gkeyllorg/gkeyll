#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_pkpm_temp_source.h>
#include <gkyl_pkpm_temp_source_priv.h>

// Thermal state of one cell at the nodes.
struct pkpm_temp_source_nodal {
  double eps[PKPM_MAX_CONF_BASIS]; // thermal energy E - |M|^2/(2 rho)
  double eps_target[PKPM_MAX_CONF_BASIS]; // its target (3/2) rho vtsq_Q
  double ratio[PKPM_MAX_CONF_BASIS]; // eps_target/eps (1 at an inadmissible node)
  bool has_density[PKPM_MAX_CONF_BASIS]; // whether the density is positive
};

static void
pkpm_temp_source_nodal_state(
  const struct gkyl_pkpm_temp_source *up, const double *fluid_d, const double *vtsq_d,
  struct pkpm_temp_source_nodal *state
)
{
  int num_nodes = up->num_nodes;
  double fluid_n[PKPM_FLUID_NUM_COMP][PKPM_MAX_CONF_BASIS];
  for (int c = 0; c < PKPM_FLUID_NUM_COMP; ++c) {
    up->conf_to_nodal(&fluid_d[c * num_nodes], fluid_n[c]);
  }
  for (int q = 0; q < num_nodes; ++q) {
    double rho = fluid_n[PKPM_FLUID_RHO][q];
    const double mom[3] = {
      fluid_n[PKPM_FLUID_MOM][q], fluid_n[PKPM_FLUID_MOM + 1][q], fluid_n[PKPM_FLUID_MOM + 2][q]
    };
    double mom_sq = mom[0] * mom[0] + mom[1] * mom[1] + mom[2] * mom[2];
    state->has_density[q] = rho > 0.0;
    state->eps[q] = state->has_density[q] ? fluid_n[PKPM_FLUID_ENERGY][q] - 0.5 * mom_sq / rho :
                                            0.0;
    state->eps_target[q] = 1.5 * rho * vtsq_d[q];
    state->ratio[q] = state->eps[q] > 0.0 ? state->eps_target[q] / state->eps[q] : 1.0;
  }
}

struct gkyl_pkpm_temp_source *
gkyl_pkpm_temp_source_new(const struct gkyl_basis *conf_basis, bool use_gpu)
{
  assert(!use_gpu);
  assert(conf_basis->b_type == GKYL_BASIS_MODAL_SERENDIPITY && conf_basis->poly_order == 1);

  struct gkyl_pkpm_temp_source *up = gkyl_malloc(sizeof(*up));
  up->num_nodes = conf_basis->num_basis;
  up->conf_to_nodal = pkpm_temp_source_conf_to_nodal_list[conf_basis->ndim - 1];
  up->nodal_to_conf = pkpm_temp_source_nodal_to_conf_list[conf_basis->ndim - 1];
  return up;
}

void
gkyl_pkpm_temp_source_advance(
  const struct gkyl_pkpm_temp_source *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *rate, const struct gkyl_array *vtsq, const struct gkyl_array *fluid,
  struct gkyl_array *rhs, struct gkyl_array *cflrate, struct gkyl_array *shape_rate
)
{
  int num_nodes = up->num_nodes;
  struct pkpm_temp_source_nodal state;
  double incr[PKPM_MAX_CONF_BASIS], incr_modal[PKPM_MAX_CONF_BASIS];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);
    const double *rate_d = gkyl_array_cfetch(rate, loc);
    double *shape_rate_d = gkyl_array_fetch(shape_rate, loc);

    pkpm_temp_source_nodal_state(
      up, gkyl_array_cfetch(fluid, loc), gkyl_array_cfetch(vtsq, loc), &state
    );

    double rate_max = 0.0;
    for (int q = 0; q < num_nodes; ++q) {
      rate_max = fmax(rate_max, rate_d[q]);
      shape_rate_d[q] = rate_d[q] * state.ratio[q];
      incr[q] = state.has_density[q] ? rate_d[q] * (state.eps_target[q] - state.eps[q]) : 0.0;
    }

    up->nodal_to_conf(incr, incr_modal);
    double *rhs_energy = (double *)gkyl_array_fetch(rhs, loc) + PKPM_FLUID_ENERGY * num_nodes;
    for (int k = 0; k < num_nodes; ++k) {
      rhs_energy[k] += incr_modal[k];
    }
    double *cflrate_d = gkyl_array_fetch(cflrate, loc);
    cflrate_d[0] += rate_max;
  }
}

void
gkyl_pkpm_temp_source_implicit_advance(
  const struct gkyl_pkpm_temp_source *up, const struct gkyl_range *conf_range, double dt,
  const struct gkyl_array *rate, const struct gkyl_array *vtsq, struct gkyl_array *fluid,
  struct gkyl_array *shape, struct gkyl_array *temp_ratio
)
{
  int num_nodes = up->num_nodes;
  struct pkpm_temp_source_nodal state;
  double incr[PKPM_MAX_CONF_BASIS], incr_modal[PKPM_MAX_CONF_BASIS];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);
    const double *rate_d = gkyl_array_cfetch(rate, loc);
    double *fluid_d = gkyl_array_fetch(fluid, loc);
    double *shape_d = gkyl_array_fetch(shape, loc);
    double *temp_ratio_d = gkyl_array_fetch(temp_ratio, loc);

    pkpm_temp_source_nodal_state(up, fluid_d, gkyl_array_cfetch(vtsq, loc), &state);

    for (int q = 0; q < num_nodes; ++q) {
      double nu_dt = rate_d[q] * dt;
      shape_d[q] = 1.0 / (1.0 + nu_dt * state.ratio[q]);
      incr[q] =
        state.has_density[q] ? nu_dt * (state.eps_target[q] - state.eps[q]) / (1.0 + nu_dt) : 0.0;
      // energy before the step over that after it (the density does not change)
      temp_ratio_d[q] = state.eps[q] > 0.0 ? state.eps[q] / (state.eps[q] + incr[q]) : 1.0;
    }

    up->nodal_to_conf(incr, incr_modal);
    double *energy = fluid_d + PKPM_FLUID_ENERGY * num_nodes;
    for (int k = 0; k < num_nodes; ++k) {
      energy[k] += incr_modal[k];
    }
  }
}

void
gkyl_pkpm_temp_source_release(struct gkyl_pkpm_temp_source *up)
{
  gkyl_free(up);
}
