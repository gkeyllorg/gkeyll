#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_dg_conf_nodal.h>
#include <gkyl_dg_conf_nodal_priv.h>

// Apply a per-expansion kernel to every block of every cell of the range.
static void
dg_conf_nodal_apply(
  dg_conf_nodal_kern_t kernel, int num_basis, const struct gkyl_range *range,
  const struct gkyl_array *fin, struct gkyl_array *fout
)
{
  assert(kernel);
  assert(fin->ncomp == fout->ncomp && fin->ncomp % num_basis == 0);
  int num_blocks = fin->ncomp / num_basis;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);
    const double *fin_d = gkyl_array_cfetch(fin, loc);
    double *fout_d = gkyl_array_fetch(fout, loc);
    for (int b = 0; b < num_blocks; ++b) {
      kernel(&fin_d[b * num_basis], &fout_d[b * num_basis]);
    }
  }
}

struct gkyl_dg_conf_nodal *
gkyl_dg_conf_nodal_new(
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, bool use_gpu
)
{
  assert(!use_gpu);
  assert(conf_basis->b_type == GKYL_BASIS_MODAL_SERENDIPITY && conf_basis->poly_order == 1);

  struct gkyl_dg_conf_nodal *up = gkyl_malloc(sizeof(*up));
  int cdim = conf_basis->ndim;
  up->num_conf_basis = conf_basis->num_basis;
  up->conf_to_nodal = dg_conf_nodal_conf_to_nodal_kernels.kernels[cdim - 1];
  up->conf_to_modal = dg_conf_nodal_conf_to_modal_kernels.kernels[cdim - 1];

  up->num_phase_basis = 0;
  up->phase_to_nodal = 0;
  up->phase_to_modal = 0;
  if (phase_basis) {
    assert(phase_basis->b_type == GKYL_BASIS_MODAL_HYBRID && phase_basis->poly_order == 1);
    int vdim = phase_basis->ndim - cdim;
    up->num_phase_basis = phase_basis->num_basis;
    up->phase_to_nodal = dg_conf_nodal_phase_to_nodal_kernels[cdim - 1].kernels[vdim - 1];
    up->phase_to_modal = dg_conf_nodal_phase_to_modal_kernels[cdim - 1].kernels[vdim - 1];
    assert(up->phase_to_nodal && up->phase_to_modal);
  }
  return up;
}

void
gkyl_dg_conf_nodal_conf_to_nodal(
  const struct gkyl_dg_conf_nodal *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *fmodal, struct gkyl_array *fnodal
)
{
  dg_conf_nodal_apply(up->conf_to_nodal, up->num_conf_basis, conf_range, fmodal, fnodal);
}

void
gkyl_dg_conf_nodal_conf_to_modal(
  const struct gkyl_dg_conf_nodal *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *fnodal, struct gkyl_array *fmodal
)
{
  dg_conf_nodal_apply(up->conf_to_modal, up->num_conf_basis, conf_range, fnodal, fmodal);
}

void
gkyl_dg_conf_nodal_phase_to_nodal(
  const struct gkyl_dg_conf_nodal *up, const struct gkyl_range *phase_range,
  const struct gkyl_array *fmodal, struct gkyl_array *fnodal
)
{
  dg_conf_nodal_apply(up->phase_to_nodal, up->num_phase_basis, phase_range, fmodal, fnodal);
}

void
gkyl_dg_conf_nodal_phase_to_modal_accumulate(
  const struct gkyl_dg_conf_nodal *up, const struct gkyl_range *phase_range,
  const struct gkyl_array *fnodal, struct gkyl_array *fmodal
)
{
  dg_conf_nodal_apply(up->phase_to_modal, up->num_phase_basis, phase_range, fnodal, fmodal);
}

void
gkyl_dg_conf_nodal_release(struct gkyl_dg_conf_nodal *up)
{
  gkyl_free(up);
}
