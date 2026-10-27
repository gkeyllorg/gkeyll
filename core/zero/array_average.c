#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_average.h>
#include <gkyl_array_average_priv.h>
#include <gkyl_array_ops.h>

#include <assert.h>
#include <string.h>

struct gkyl_array_average *
gkyl_array_average_new(
  const struct gkyl_rect_grid *grid, const struct gkyl_basis *basis,
  const struct gkyl_basis *basis_avg, const struct gkyl_range *local,
  const struct gkyl_range *local_avg, const struct gkyl_range *local_avg_ext,
  const struct gkyl_array *weight, const int *avg_dim, bool use_gpu
)
{
  return gkyl_array_average_inew(&(struct gkyl_array_average_inp){
    .grid = grid,
    .basis = *basis,
    .basis_avg = *basis_avg,
    .local = local,
    .local_avg = local_avg,
    .local_avg_ext = local_avg_ext,
    .weight = weight,
    .avg_dim = avg_dim,
    .use_gpu = use_gpu,
  });
}

struct gkyl_array_average *
gkyl_array_average_inew(const struct gkyl_array_average_inp *inp)
{
  struct gkyl_array_average *up = gkyl_malloc(sizeof(struct gkyl_array_average));

  up->use_gpu = inp->use_gpu;
  up->ndim = inp->basis.ndim;
  up->vdim = inp->vdim;
  up->basis = inp->basis;
  up->local = *inp->local;
  up->full_dim_out = inp->full_dim_out;
  if (up->full_dim_out) {
    up->basis_avg = inp->basis;
    up->local_avg = *inp->local;
  } else {
    up->basis_avg = inp->basis_avg;
    up->local_avg = *inp->local_avg;
  }

  up->num_avg_dim = 0;
  for (int d = 0; d < up->ndim; ++d) {
    up->avg_dim[d] = inp->avg_dim[d];
    up->dim_remains[d] = 1 - inp->avg_dim[d];
    up->num_avg_dim += inp->avg_dim[d];
  }
  assert(up->num_avg_dim > 0);
  up->num_dim_remain = up->ndim - up->num_avg_dim;

  // Inverse of the domain volume and cell volume of the averaged dimensions.
  up->vol_avg_inv = 1.0;
  up->subvol = 1.0;
  for (int d = 0; d < up->ndim; d++) {
    if (up->avg_dim[d]) {
      up->vol_avg_inv *= inp->grid->upper[d] - inp->grid->lower[d];
      up->subvol *= 0.5 * inp->grid->dx[d];
    }
  }
  up->vol_avg_inv = 1.0 / up->vol_avg_inv;

  up->isweighted = false;
  up->weight = NULL;
  up->weight_avg = NULL;
  up->div_mem = NULL;
  if (inp->weight) {
    // The weak division at the end of the advance is only available for conf-space p<=2.
    assert(up->vdim == 0 && up->basis.poly_order <= 2);
    up->isweighted = true;
    // Integrate the weight over the averaged dimensions with an unweighted average.
    up->weight_avg =
      up->use_gpu ?
        gkyl_array_cu_dev_new(GKYL_DOUBLE, up->basis_avg.num_basis, inp->local_avg_ext->volume) :
        gkyl_array_new(GKYL_DOUBLE, up->basis_avg.num_basis, inp->local_avg_ext->volume);
    struct gkyl_array_average_inp inp_integral = {
      .grid = inp->grid,
      .basis = inp->basis,
      .basis_avg = inp->basis_avg,
      .vdim = inp->vdim,
      .full_dim_out = inp->full_dim_out,
      .local = inp->local,
      .local_avg = inp->local_avg,
      .weight = NULL,
      .avg_dim = inp->avg_dim,
      .use_gpu = inp->use_gpu,
    };
    struct gkyl_array_average *int_w = gkyl_array_average_inew(&inp_integral);
    gkyl_array_average_advance(int_w, inp->weight, up->weight_avg);
    gkyl_array_average_release(int_w);
    gkyl_array_scale(up->weight_avg, 1.0 / up->vol_avg_inv);

    up->div_mem = up->use_gpu ?
                    gkyl_dg_bin_op_mem_cu_dev_new(up->local_avg.volume, up->basis_avg.num_basis) :
                    gkyl_dg_bin_op_mem_new(up->local_avg.volume, up->basis_avg.num_basis);
    up->weight = gkyl_array_acquire(inp->weight);
  } else {
    // Unit weight, and divide by the domain volume through subvol.
    up->weight = up->use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, up->basis.num_basis, 1) :
                               gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, 1);
    gkyl_array_shiftc(up->weight, pow(sqrt(2.0), up->ndim), 0);
    up->subvol *= up->vol_avg_inv;
  }

  array_average_choose_kernel(up);

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    return gkyl_array_average_cu_dev_new(up);
  }
#endif

  return up;
}

void
gkyl_array_average_advance(
  const struct gkyl_array_average *up, const struct gkyl_array *fin, struct gkyl_array *avgout
)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_array_average_advance_cu(up, fin, avgout);
    return;
  }
#endif

  gkyl_array_clear_range(avgout, 0.0, &up->local_avg);

  // Accumulate the contribution of every input cell into its output cell.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->local);
  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(&up->local, iter.idx);
    long lidx_out = array_average_out_lidx(up, iter.idx);

    const double *fin_c = gkyl_array_cfetch(fin, lidx);
    const double *win_c = gkyl_array_cfetch(up->weight, up->isweighted ? lidx : 0);
    double *out_c = gkyl_array_fetch(avgout, lidx_out);

    up->kernel(up->subvol, win_c, fin_c, out_c);
  }

  if (up->full_dim_out) {
    // Copy the average into every other cell along the averaged dimensions.
    gkyl_range_iter_init(&iter, &up->local_avg);
    while (gkyl_range_iter_next(&iter)) {
      long lidx = gkyl_range_idx(&up->local_avg, iter.idx);
      long lidx_rep = array_average_out_lidx(up, iter.idx);
      if (lidx != lidx_rep) {
        memcpy(gkyl_array_fetch(avgout, lidx), gkyl_array_cfetch(avgout, lidx_rep), avgout->esznc);
      }
    }
  }

  if (up->isweighted) {
    gkyl_dg_div_op_range(
      up->div_mem, &up->basis_avg, 0, avgout, 0, avgout, 0, up->weight_avg, &up->local_avg
    );
  }
}

void
gkyl_array_average_release(struct gkyl_array_average *up)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_cu_free(up->on_dev);
  }
#endif

  gkyl_array_release(up->weight);
  if (up->isweighted) {
    gkyl_array_release(up->weight_avg);
    gkyl_dg_bin_op_mem_release(up->div_mem);
  }

  gkyl_free(up);
}
