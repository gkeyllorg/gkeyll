/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_cx.h>
#include <gkyl_dg_cx_priv.h>
#include <gkyl_util.h>
#include <gkyl_const.h>
}

// CUDA kernel to set device pointers to kernels.
__global__ static void
gkyl_dg_cx_set_cu_dev_ptrs(struct gkyl_dg_cx *up, struct gkyl_basis cbasis)
{
  up->react_rate = choose_kern(cbasis);
};

__global__ static void
gkyl_cx_react_rate_cu_ker(struct gkyl_dg_cx *up, const struct gkyl_range conf_rng, 
  const struct gkyl_array *maxwellian_moms_ion, const struct gkyl_array *maxwellian_moms_neut,
  const struct gkyl_array *upar_b_i,  double vt_sq_ion_min, double vt_sq_neut_min,
  struct gkyl_array *coef_cx, double a, double b)
{
  int cidx[GKYL_MAX_CDIM];
  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < conf_rng.volume; tid += blockDim.x*gridDim.x) {

    gkyl_sub_range_inv_idx(&conf_rng, tid, cidx);
    long linidx = gkyl_range_idx(&conf_rng, cidx);

    const double *maxwellian_moms_ion_d = (const double*) gkyl_array_cfetch(maxwellian_moms_ion, linidx);
    const double *maxwellian_moms_neut_d = (const double*) gkyl_array_cfetch(maxwellian_moms_neut, linidx);
    const double *upar_b_i_d = (const double*) gkyl_array_cfetch(upar_b_i, linidx);

    double *coef_cx_d = (double*) gkyl_array_fetch(coef_cx, linidx);

    // call the cx kernel
    double cflr = up->react_rate(a, b, vt_sq_ion_min, vt_sq_neut_min, 
      maxwellian_moms_ion_d, maxwellian_moms_neut_d, upar_b_i_d, coef_cx_d);
  }
}

void gkyl_dg_cx_coll_cu(const struct gkyl_dg_cx *up, 
  struct gkyl_array *maxwellian_moms_ion, struct gkyl_array *maxwellian_moms_neut,
  struct gkyl_array *upar_b_i, struct gkyl_array *coef_cx, struct gkyl_array *cflrate)
{  
  gkyl_cx_react_rate_cu_ker<<<up->conf_rng->nblocks, up->conf_rng->nthreads>>>(up->on_dev, *up->conf_rng,
    maxwellian_moms_ion->on_dev, maxwellian_moms_neut->on_dev, upar_b_i->on_dev, 
    up->vt_sq_ion_min, up->vt_sq_neut_min, coef_cx->on_dev, up->a, up->b);
}

gkyl_dg_cx*
gkyl_dg_cx_cu_dev_new(struct gkyl_dg_cx_inp *inp)
{
  gkyl_dg_cx *up = (struct gkyl_dg_cx*) gkyl_malloc(sizeof(*up));

  up->cbasis = inp->cbasis;
  up->conf_rng = inp->conf_rng;
  up->type_ion = inp->type_ion;
  up->vt_sq_ion_min = inp->vt_sq_ion_min;
  up->vt_sq_neut_min = inp->vt_sq_neut_min;

  fit_param(up->type_ion, &up->a, &up->b);
  
  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_dg_cx *up_cu = (struct gkyl_dg_cx*) gkyl_cu_malloc(sizeof(*up_cu));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_dg_cx), GKYL_CU_MEMCPY_H2D);

  gkyl_dg_cx_set_cu_dev_ptrs<<<1,1>>>(up_cu, *up->cbasis);

  // Set parent on_dev pointer.
  up->on_dev = up_cu;
  
  return up;
}
