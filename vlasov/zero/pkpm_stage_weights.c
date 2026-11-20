#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_pkpm_stage_weights.h>
#include <gkyl_pkpm_stage_weights_priv.h>

struct gkyl_pkpm_stage_weights *
gkyl_pkpm_stage_weights_inew(const struct gkyl_pkpm_stage_weights_inp *inp)
{
  assert(!inp->use_gpu);
  assert(
    inp->conf_basis->b_type == GKYL_BASIS_MODAL_SERENDIPITY && inp->conf_basis->poly_order == 1
  );
  assert(inp->phase_basis->b_type == GKYL_BASIS_MODAL_HYBRID);
  assert(inp->phase_basis->num_basis == PKPM_NUM_VEL_MODES * inp->conf_basis->num_basis);

  struct gkyl_pkpm_stage_weights *up = gkyl_malloc(sizeof(*up));
  int cdim = inp->conf_basis->ndim, vdim = inp->phase_basis->ndim - cdim;
  up->num_nodes = inp->conf_basis->num_basis;
  up->num_phase_basis = inp->phase_basis->num_basis;
  up->conf_to_modal = dg_conf_nodal_conf_to_modal_kernels.kernels[cdim - 1];
  up->conf_phase_mul_accumulate = choose_mul_conf_phase_accumulate_kern(
    inp->phase_basis->b_type, cdim, vdim, inp->conf_basis->poly_order
  );
  assert(up->conf_phase_mul_accumulate);
  return up;
}

void
gkyl_pkpm_stage_weights_advance(
  const struct gkyl_pkpm_stage_weights *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, double coef_fe, const struct gkyl_array *f_old,
  const struct gkyl_array *f_fe, struct gkyl_array *lam, struct gkyl_array *f
)
{
  int num_nodes = up->num_nodes, num_phase_basis = up->num_phase_basis;
  int num_blocks = f->ncomp / num_phase_basis;
  double coef_old = 1.0 - coef_fe;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);
  while (gkyl_range_iter_next(&iter)) {
    long ploc = gkyl_range_idx(phase_range, iter.idx);
    const double *lam_d = gkyl_array_cfetch(lam, gkyl_range_idx(conf_range, iter.idx));
    const double *f_old_d = gkyl_array_cfetch(f_old, ploc);
    const double *f_fe_d = gkyl_array_cfetch(f_fe, ploc);
    double *f_d = gkyl_array_fetch(f, ploc);

    // w_Y - c_Y at the nodes, as an expansion
    double dweight[PKPM_MAX_CONF_BASIS], dweight_modal[PKPM_MAX_CONF_BASIS];
    for (int q = 0; q < num_nodes; ++q) {
      dweight[q] = coef_fe * lam_d[q] / (coef_old + coef_fe * lam_d[q]) - coef_fe;
    }
    up->conf_to_modal(dweight, dweight_modal);

    for (int b = 0; b < num_blocks; ++b) {
      double *f_b = &f_d[b * num_phase_basis];
      up->conf_phase_mul_accumulate(1.0, dweight_modal, &f_fe_d[b * num_phase_basis], f_b);
      up->conf_phase_mul_accumulate(-1.0, dweight_modal, &f_old_d[b * num_phase_basis], f_b);
    }
  }

  // temperature factor of the combined state
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    double *lam_d = gkyl_array_fetch(lam, gkyl_range_idx(conf_range, iter.idx));
    for (int q = 0; q < num_nodes; ++q) {
      lam_d[q] = coef_old + coef_fe * lam_d[q];
    }
  }
}

void
gkyl_pkpm_stage_weights_scale(
  const struct gkyl_pkpm_stage_weights *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *temp_ratio, struct gkyl_array *lam
)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);
    const double *ratio_d = gkyl_array_cfetch(temp_ratio, loc);
    double *lam_d = gkyl_array_fetch(lam, loc);
    for (int q = 0; q < up->num_nodes; ++q) {
      lam_d[q] /= ratio_d[q];
    }
  }
}

void
gkyl_pkpm_stage_weights_release(struct gkyl_pkpm_stage_weights *up)
{
  gkyl_free(up);
}
