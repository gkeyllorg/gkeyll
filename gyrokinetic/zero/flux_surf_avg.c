#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_flux_surf_avg.h>
#include <gkyl_flux_surf_avg_priv.h>

#include <assert.h>

static inline struct gkyl_array *
mkarr(bool use_gpu, long ncomp, long size)
{
  return use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, ncomp, size) :
                   gkyl_array_new(GKYL_DOUBLE, ncomp, size);
}

static void
flux_surf_avg_reduce_comm(
  struct gkyl_flux_surf_avg *up, const struct gkyl_array *avg_jf, struct gkyl_array *avg_jf_red
)
{
  gkyl_comm_allreduce(
    up->comm, GKYL_DOUBLE, GKYL_SUM, avg_jf->ncomp * avg_jf->size, avg_jf->data, avg_jf_red->data
  );
}

static void
flux_surf_avg_reduce_none(
  struct gkyl_flux_surf_avg *up, const struct gkyl_array *avg_jf, struct gkyl_array *avg_jf_red
)
{
  gkyl_array_copy(avg_jf_red, avg_jf);
}

static void
flux_surf_avg_inflate_2x(
  struct gkyl_flux_surf_avg *up, const struct gkyl_array *avg, struct gkyl_array *out
)
{
  gkyl_translate_dim_advance(up->inflate_lo, &up->local_x, &up->local, avg, 1, out);
}

static void
flux_surf_avg_inflate_3x(
  struct gkyl_flux_surf_avg *up, const struct gkyl_array *avg, struct gkyl_array *out
)
{
  gkyl_translate_dim_advance(up->inflate_lo, &up->local_x, &up->local_xy, avg, 1, up->inflate_buff);
  gkyl_translate_dim_advance(up->inflate_up, &up->local_xy, &up->local, up->inflate_buff, 1, out);
}

// Flux-surface average of fin into the 1D field avg: <f> = (sum over ranks of the local (y,z)
// average of J*f) / (global (y,z) average of J).
static void
flux_surf_avg_advance_deflated(
  struct gkyl_flux_surf_avg *up, const struct gkyl_array *fin, struct gkyl_array *avg
)
{
  gkyl_dg_mul_op_range(&up->basis, 0, up->jf, 0, up->jacobgeo, 0, fin, &up->local);
  gkyl_array_average_advance(up->avg, up->jf, up->avg_jf);
  up->reduce_func(up, up->avg_jf, up->avg_jf_red);
  gkyl_dg_div_op_range(
    up->div_mem, &up->basis_x, 0, avg, 0, up->avg_jf_red, 0, up->avg_j, &up->local_x
  );
}

static void
flux_surf_avg_advance_inflated(
  struct gkyl_flux_surf_avg *up, const struct gkyl_array *fin, struct gkyl_array *out
)
{
  flux_surf_avg_advance_deflated(up, fin, up->avg_buff);
  up->inflate_func(up, up->avg_buff, out);
}

struct gkyl_flux_surf_avg *
gkyl_flux_surf_avg_inew(const struct gkyl_flux_surf_avg_inp *inp)
{
  struct gkyl_flux_surf_avg *up = gkyl_malloc(sizeof(*up));

  int cdim = inp->basis.ndim;
  assert(cdim == 2 || cdim == 3);
  assert(inp->basis.poly_order == 1); // array_average and translate_dim kernels are p=1 only.
  assert(inp->basis_x.ndim == 1);

  up->basis = inp->basis;
  up->basis_x = inp->basis_x;
  up->local = *inp->local;
  up->local_x = *inp->local_x;
  up->use_gpu = inp->use_gpu;
  up->jacobgeo = gkyl_array_acquire(inp->jacobgeo);

  up->comm = 0;
  up->reduce_func = flux_surf_avg_reduce_none;
  if (inp->comm) {
    up->comm = gkyl_comm_acquire(inp->comm);
    up->reduce_func = flux_surf_avg_reduce_comm;
  }

  int nb = up->basis.num_basis, nb_x = up->basis_x.num_basis;
  long vol_x = inp->local_x_ext->volume;
  up->jf = mkarr(up->use_gpu, nb, inp->local_ext->volume);
  up->avg_jf = mkarr(up->use_gpu, nb_x, vol_x);
  up->avg_jf_red = mkarr(up->use_gpu, nb_x, vol_x);
  up->avg_j = mkarr(up->use_gpu, nb_x, vol_x);
  up->div_mem = up->use_gpu ? gkyl_dg_bin_op_mem_cu_dev_new(up->local_x.volume, nb_x) :
                              gkyl_dg_bin_op_mem_new(up->local_x.volume, nb_x);

  // Unweighted average over all directions but x. It divides by the global (y,z) extent, so
  // the sum over ranks is the global average and the weight J is divided out afterwards.
  int avg_dim[GKYL_MAX_CDIM] = {0};
  for (int d = 1; d < cdim; d++) {
    avg_dim[d] = 1;
  }
  up->avg = gkyl_array_average_inew(&(struct gkyl_array_average_inp){
    .grid = inp->grid,
    .basis = up->basis,
    .basis_avg = up->basis_x,
    .local = &up->local,
    .local_avg = &up->local_x,
    .local_avg_ext = inp->local_x_ext,
    .weight = NULL,
    .avg_dim = avg_dim,
    .use_gpu = up->use_gpu,
  });

  // Global (y,z) average of J (time independent).
  gkyl_array_average_advance(up->avg, up->jacobgeo, up->avg_jf);
  up->reduce_func(up, up->avg_jf, up->avg_j);

  // Updaters extending a function of x to cdim dimensions.
  up->inflate_up = 0;
  up->inflate_buff = 0;
  if (cdim == 2) {
    up->inflate_lo =
      gkyl_translate_dim_new(1, up->basis_x, 2, up->basis, 0, GKYL_NO_EDGE, up->use_gpu);
    up->inflate_func = flux_surf_avg_inflate_2x;
  } else {
    struct gkyl_basis basis_xy;
    gkyl_cart_modal_serendip(&basis_xy, 2, up->basis.poly_order);
    gkyl_range_init(&up->local_xy_ext, 2, inp->local_ext->lower, inp->local_ext->upper);
    gkyl_sub_range_init(&up->local_xy, &up->local_xy_ext, up->local.lower, up->local.upper);
    up->inflate_buff = mkarr(up->use_gpu, basis_xy.num_basis, up->local_xy_ext.volume);
    up->inflate_lo =
      gkyl_translate_dim_new(1, up->basis_x, 2, basis_xy, 0, GKYL_NO_EDGE, up->use_gpu);
    up->inflate_up =
      gkyl_translate_dim_new(2, basis_xy, 3, up->basis, 0, GKYL_NO_EDGE, up->use_gpu);
    up->inflate_func = flux_surf_avg_inflate_3x;
  }

  up->avg_buff = 0;
  up->advance_func = flux_surf_avg_advance_deflated;
  if (!inp->is_deflated) {
    up->avg_buff = mkarr(up->use_gpu, nb_x, vol_x);
    up->advance_func = flux_surf_avg_advance_inflated;
  }

  return up;
}

void
gkyl_flux_surf_avg_advance(
  struct gkyl_flux_surf_avg *up, const struct gkyl_array *fin, struct gkyl_array *out
)
{
  up->advance_func(up, fin, out);
}

void
gkyl_flux_surf_avg_inflate(
  struct gkyl_flux_surf_avg *up, const struct gkyl_array *avg, struct gkyl_array *out
)
{
  up->inflate_func(up, avg, out);
}

void
gkyl_flux_surf_avg_release(struct gkyl_flux_surf_avg *up)
{
  if (up->comm) {
    gkyl_comm_release(up->comm);
  }
  gkyl_array_release(up->jacobgeo);
  gkyl_array_average_release(up->avg);
  gkyl_array_release(up->jf);
  gkyl_array_release(up->avg_jf);
  gkyl_array_release(up->avg_jf_red);
  gkyl_array_release(up->avg_j);
  gkyl_dg_bin_op_mem_release(up->div_mem);
  gkyl_translate_dim_release(up->inflate_lo);
  if (up->inflate_up) {
    gkyl_translate_dim_release(up->inflate_up);
    gkyl_array_release(up->inflate_buff);
  }
  if (up->avg_buff) {
    gkyl_array_release(up->avg_buff);
  }
  gkyl_free(up);
}
