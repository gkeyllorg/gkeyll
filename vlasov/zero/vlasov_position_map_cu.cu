/* -*- c++ -*- */

#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_position_map.h>
#include <gkyl_vlasov_position_map_priv.h>
}

__global__ void
gkyl_vlasov_position_map_divide_jacobpos_cu_kernel(
  int cdim, int num_basis, struct gkyl_range conf_range, struct gkyl_range phase_range,
  const struct gkyl_array *jacob_pos_gauss, const struct gkyl_array *Jf, struct gkyl_array *f_no_J
)
{
  int idx[GKYL_MAX_DIM], idx_conf[GKYL_MAX_CDIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < phase_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&phase_range, linc1, idx);

    for (int i = 0; i < cdim; ++i) {
      idx_conf[i] = idx[i];
    }
    long cidx = gkyl_range_idx(&conf_range, idx_conf);
    long pidx = gkyl_range_idx(&phase_range, idx);

    double jacob_inv = 1.0 / ((const double *)gkyl_array_cfetch(jacob_pos_gauss, cidx))[0];
    const double *Jf_d = (const double *)gkyl_array_cfetch(Jf, pidx);
    double *f_no_J_d = (double *)gkyl_array_fetch(f_no_J, pidx);
    for (int k = 0; k < num_basis; ++k) {
      f_no_J_d[k] = jacob_inv * Jf_d[k];
    }
  }
}

// Configuration-space field scaled per cell by 1/J (first num_coeff_divide
// coefficients) or by J (all coefficients), only on conf_range.
__global__ void
gkyl_vlasov_position_map_scale_jacobpos_conf_cu_kernel(
  struct gkyl_range jac_range, struct gkyl_range conf_range, int ncomp, int num_coeff_scale,
  bool divide, const struct gkyl_array *jacob_pos_gauss, const struct gkyl_array *in,
  struct gkyl_array *out
)
{
  int idx[GKYL_MAX_CDIM];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < conf_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);
    long cidx_jac = gkyl_range_idx(&jac_range, idx);
    long cidx = gkyl_range_idx(&conf_range, idx);

    double jacob = ((const double *)gkyl_array_cfetch(jacob_pos_gauss, cidx_jac))[0];
    double fac = divide ? 1.0 / jacob : jacob;
    const double *in_d = (const double *)gkyl_array_cfetch(in, cidx);
    double *out_d = (double *)gkyl_array_fetch(out, cidx);
    for (int k = 0; k < num_coeff_scale; ++k) {
      out_d[k] = fac * in_d[k];
    }
    for (int k = num_coeff_scale; k < ncomp; ++k) {
      out_d[k] = in_d[k];
    }
  }
}

void
gkyl_vlasov_position_map_divide_jacobpos_conf_cu(
  const struct gkyl_range *jac_range, const struct gkyl_range *conf_range, int num_coeff_divide,
  const struct gkyl_array *jacob_pos_gauss, const struct gkyl_array *Jmom,
  struct gkyl_array *mom_no_J
)
{
  gkyl_vlasov_position_map_scale_jacobpos_conf_cu_kernel<<<
    conf_range->nblocks, conf_range->nthreads>>>(
    *jac_range, *conf_range, Jmom->ncomp, num_coeff_divide, true, jacob_pos_gauss->on_dev,
    Jmom->on_dev, mom_no_J->on_dev
  );
}

void
gkyl_vlasov_position_map_rescale_jacobpos_conf_cu(
  const struct gkyl_range *jac_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *jacob_pos_gauss, const struct gkyl_array *a_no_J, struct gkyl_array *Ja
)
{
  gkyl_vlasov_position_map_scale_jacobpos_conf_cu_kernel<<<
    conf_range->nblocks, conf_range->nthreads>>>(
    *jac_range, *conf_range, a_no_J->ncomp, a_no_J->ncomp, false, jacob_pos_gauss->on_dev,
    a_no_J->on_dev, Ja->on_dev
  );
}

// Host-side wrapper for dividing out the configuration-space Jacobian from Jf.
void
gkyl_vlasov_position_map_divide_jacobpos_cu(
  int cdim, const struct gkyl_basis *phase_basis, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *jacob_pos_gauss,
  const struct gkyl_array *Jf, struct gkyl_array *f_no_J
)
{
  int nblocks = phase_range->nblocks;
  int nthreads = phase_range->nthreads;
  gkyl_vlasov_position_map_divide_jacobpos_cu_kernel<<<nblocks, nthreads>>>(
    cdim, phase_basis->num_basis, *conf_range, *phase_range, jacob_pos_gauss->on_dev, Jf->on_dev,
    f_no_J->on_dev
  );
}

__global__ void
gkyl_vlasov_position_map_rescale_jacobpos_cu_kernel(
  int cdim, int num_basis, struct gkyl_range conf_range, struct gkyl_range phase_range,
  const struct gkyl_array *jacob_pos_gauss, const struct gkyl_array *f_no_J, struct gkyl_array *Jf
)
{
  int idx[GKYL_MAX_DIM], idx_conf[GKYL_MAX_CDIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < phase_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&phase_range, linc1, idx);

    for (int i = 0; i < cdim; ++i) {
      idx_conf[i] = idx[i];
    }
    long cidx = gkyl_range_idx(&conf_range, idx_conf);
    long pidx = gkyl_range_idx(&phase_range, idx);

    double jacob = ((const double *)gkyl_array_cfetch(jacob_pos_gauss, cidx))[0];
    const double *f_no_J_d = (const double *)gkyl_array_cfetch(f_no_J, pidx);
    double *Jf_d = (double *)gkyl_array_fetch(Jf, pidx);
    for (int k = 0; k < num_basis; ++k) {
      Jf_d[k] = jacob * f_no_J_d[k];
    }
  }
}

// Host-side wrapper for multiplying the configuration-space Jacobian into f to obtain Jf.
void
gkyl_vlasov_position_map_rescale_jacobpos_cu(
  int cdim, const struct gkyl_basis *phase_basis, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *jacob_pos_gauss,
  const struct gkyl_array *f_no_J, struct gkyl_array *Jf
)
{
  int nblocks = phase_range->nblocks;
  int nthreads = phase_range->nthreads;
  gkyl_vlasov_position_map_rescale_jacobpos_cu_kernel<<<nblocks, nthreads>>>(
    cdim, phase_basis->num_basis, *conf_range, *phase_range, jacob_pos_gauss->on_dev,
    f_no_J->on_dev, Jf->on_dev
  );
}
