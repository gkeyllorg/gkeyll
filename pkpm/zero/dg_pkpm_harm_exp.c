#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_dg_pkpm_harm_exp.h>
#include <gkyl_dg_pkpm_harm_exp_priv.h>
#include <gkyl_range.h>

// Sizes and parameters shared by the host and device constructors.
static void
dg_pkpm_harm_exp_setup(struct gkyl_dg_pkpm_harm_exp *up, const struct gkyl_dg_pkpm_harm_exp_inp *inp)
{
  assert(
    inp->conf_basis->poly_order == 1 && inp->conf_basis->b_type == GKYL_BASIS_MODAL_SERENDIPITY
  );
  up->cdim = inp->conf_basis->ndim;
  up->num_conf_nodes = 1 << up->cdim;
  up->num_phase_nodes = PKPM_NUM_VEL_MODES * up->num_conf_nodes;
  up->qbym = inp->qbym;
  up->mix = inp->mix > 0.0 ? inp->mix : 0.0;
  up->half_dw = 0.5 * inp->phase_grid->dx[up->cdim];
  up->use_gpu = inp->use_gpu;
  up->flags = 0;
}

struct gkyl_dg_pkpm_harm_exp *
gkyl_dg_pkpm_harm_exp_new(const struct gkyl_dg_pkpm_harm_exp_inp *inp)
{
#ifdef GKYL_HAVE_CUDA
  if (inp->use_gpu) {
    return gkyl_dg_pkpm_harm_exp_cu_dev_new(inp);
  }
#endif
  struct gkyl_dg_pkpm_harm_exp *up = gkyl_malloc(sizeof(*up));
  dg_pkpm_harm_exp_setup(up, inp);
  up->use_gpu = false;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host
  return up;
}

int
gkyl_dg_pkpm_harm_exp_num_coef(const struct gkyl_dg_pkpm_harm_exp *up)
{
  return up->num_conf_nodes * PKPM_HX_NC;
}

int
gkyl_dg_pkpm_harm_exp_num_rate(const struct gkyl_dg_pkpm_harm_exp *up)
{
  return up->num_phase_nodes;
}

int
gkyl_dg_pkpm_harm_exp_num_col(const struct gkyl_dg_pkpm_harm_exp *up)
{
  return up->num_conf_nodes * 4;
}

void
gkyl_dg_pkpm_harm_exp_coef(
  const struct gkyl_dg_pkpm_harm_exp *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *em, struct gkyl_array *coef
)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(coef)) {
    gkyl_dg_pkpm_harm_exp_coef_cu(up, conf_range, em, coef);
    return;
  }
#endif
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);
    pkpm_hx_coef_cell(up, gkyl_array_cfetch(em, loc), gkyl_array_fetch(coef, loc));
  }
}

void
gkyl_dg_pkpm_harm_exp_advance(
  const struct gkyl_dg_pkpm_harm_exp *up, int stage, double h, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *coef, const struct gkyl_array *f0,
  const struct gkyl_array *lam, struct gkyl_array *nu, struct gkyl_array *act,
  struct gkyl_array *out
)
{
  if (up->mix <= 0.0) {
    return;
  }
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) {
    // the weighted combination is not available on the device
    assert(!lam);
    gkyl_dg_pkpm_harm_exp_advance_cu(up, stage, h, conf_range, phase_range, coef, f0, nu, act, out);
    return;
  }
#endif
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);
  while (gkyl_range_iter_next(&iter)) {
    long cloc = gkyl_range_idx(conf_range, iter.idx);
    long ploc = gkyl_range_idx(phase_range, iter.idx);
    pkpm_hx_advance_cell(
      up, stage, h, gkyl_array_cfetch(coef, cloc), gkyl_array_cfetch(f0, ploc),
      lam ? gkyl_array_cfetch(lam, cloc) : 0, gkyl_array_fetch(nu, ploc),
      gkyl_array_fetch(act, cloc), gkyl_array_fetch(out, ploc)
    );
  }
}

void
gkyl_dg_pkpm_harm_exp_restore(
  const struct gkyl_dg_pkpm_harm_exp *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *act, struct gkyl_array *col,
  struct gkyl_array *out
)
{
  if (up->mix <= 0.0) {
    return;
  }
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) {
    gkyl_dg_pkpm_harm_exp_restore_cu(up, conf_range, phase_range, act, col, out);
    return;
  }
#endif
  int cdim = up->cdim, nc = 4 * up->num_conf_nodes;
  // the w integrals of the flagged columns
  struct gkyl_range_iter citer;
  gkyl_range_iter_init(&citer, conf_range);
  while (gkyl_range_iter_next(&citer)) {
    long cloc = gkyl_range_idx(conf_range, citer.idx);
    if (*(const double *)gkyl_array_cfetch(act, cloc) == 0.0) {
      continue;
    }
    double *c = gkyl_array_fetch(col, cloc);
    for (int i = 0; i < nc; ++i) {
      c[i] = 0.0;
    }
    int pidx[GKYL_MAX_DIM];
    for (int d = 0; d < cdim; ++d) {
      pidx[d] = citer.idx[d];
    }
    for (int iv = phase_range->lower[cdim]; iv <= phase_range->upper[cdim]; ++iv) {
      pidx[cdim] = iv;
      pkpm_hx_col_add(up, gkyl_array_cfetch(out, gkyl_range_idx(phase_range, pidx)), c);
    }
  }
  // the frame shift of every phase cell of those columns
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);
  while (gkyl_range_iter_next(&iter)) {
    long cloc = gkyl_range_idx(conf_range, iter.idx);
    if (*(const double *)gkyl_array_cfetch(act, cloc) == 0.0) {
      continue;
    }
    pkpm_hx_restore_cell(
      up, gkyl_array_cfetch(col, cloc), gkyl_array_fetch(out, gkyl_range_idx(phase_range, iter.idx))
    );
  }
}

void
gkyl_dg_pkpm_harm_exp_release(struct gkyl_dg_pkpm_harm_exp *up)
{
  if (GKYL_IS_CU_ALLOC(up->flags)) {
    gkyl_cu_free(up->on_dev);
  }
  gkyl_free(up);
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_pkpm_harm_exp *
gkyl_dg_pkpm_harm_exp_cu_dev_new(const struct gkyl_dg_pkpm_harm_exp_inp *inp)
{
  assert(false);
  return 0;
}

#endif
