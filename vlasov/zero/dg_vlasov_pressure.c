
#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_bin_ops.h>

#include <gkyl_dg_vlasov_pressure.h>
#include <gkyl_dg_vlasov_pressure_priv.h>

// Index of tensor components
static inline int
tensor_idx(int i, int j, int vdim)
{
  return i*vdim - i*(i-1)/2 + (j-i);
}

struct gkyl_dg_vlasov_pressure *gkyl_dg_vlasov_pressure_new(
  const struct gkyl_basis *conf_basis,
  const struct gkyl_range *conf_range,
  const struct gkyl_range *conf_range_ext,
  int vdim, double mass, bool use_gpu
)
{
  assert(vdim > 0);
  assert(mass > 0.0);

  struct gkyl_dg_vlasov_pressure *up = gkyl_malloc(sizeof(struct gkyl_dg_vlasov_pressure));

  up->conf_basis = *conf_basis;
  up->conf_range = *conf_range;
  up->num_conf_basis = conf_basis->num_basis;
  up->vdim = vdim;
  up->mass = mass;
  up->use_gpu = use_gpu;

  int num_conf_basis = up->num_conf_basis;
  int num_tensor_comp = vdim*(vdim+1)/2;


  if (use_gpu) {
    up->V_drift_dot_M1i = gkyl_array_cu_dev_new(GKYL_DOUBLE, num_conf_basis, conf_range_ext->volume);
    up->V_drift_M1i = gkyl_array_cu_dev_new(GKYL_DOUBLE, num_tensor_comp*num_conf_basis, conf_range_ext->volume);
    up->V_drift_M1i_transpose = gkyl_array_cu_dev_new(GKYL_DOUBLE, num_tensor_comp*num_conf_basis, conf_range_ext->volume);
  }

  else {
    up->V_drift_dot_M1i = gkyl_array_new(GKYL_DOUBLE, num_conf_basis, conf_range_ext->volume);
    up->V_drift_M1i = gkyl_array_new(GKYL_DOUBLE, num_tensor_comp*num_conf_basis, conf_range_ext->volume);
    up->V_drift_M1i_transpose = gkyl_array_new(GKYL_DOUBLE, num_tensor_comp*num_conf_basis, conf_range_ext->volume);
  }

  return up;
}

void
gkyl_dg_vlasov_pressure_scalar(
struct gkyl_dg_vlasov_pressure *up,
  const struct gkyl_range *conf_local,
  const struct gkyl_array *M2,
  const struct gkyl_array *M1i,
  const struct gkyl_array *V_drift,
  struct gkyl_array *pressure
)
{
  // Compute V_drift dot M1i.
  gkyl_dg_dot_product_op_range(&up->conf_basis, up->V_drift_dot_M1i, V_drift, M1i, conf_local);

  // pressure = M2 - V_drift dot M1i.
  gkyl_array_set_range(pressure, 1.0, M2, conf_local);
  gkyl_array_accumulate_range(pressure, -1.0, up->V_drift_dot_M1i, conf_local);
  gkyl_array_scale_range(pressure, up->mass/up->vdim, conf_local);
}

void
gkyl_dg_vlasov_pressure_tensor(
  struct gkyl_dg_vlasov_pressure *up,
  const struct gkyl_range *conf_local,
  const struct gkyl_array *M2ij,
  const struct gkyl_array *M1i,
  const struct gkyl_array *V_drift,
  struct gkyl_array *pressure_tensor
)
{
  const int vdim = up->vdim;

  /*
   * Compute the symmetric pressure tensor:
   *
   * P_ij = m * [M2ij - 0.5*(V_i*M1_j + V_j*M1_i)]
   */

  for (int i = 0; i < vdim; ++i) {
    for (int j = i; j < vdim; ++j) {

      int comp = tensor_idx(i, j, vdim);

      // Calculate V_i * M1_j.
      gkyl_dg_mul_op_range(&up->conf_basis, comp, up->V_drift_M1i, i, V_drift, j, M1i, conf_local);
      gkyl_dg_mul_op_range(&up->conf_basis, comp, up->V_drift_M1i_transpose, j, V_drift, i, M1i, conf_local);
    }
  }

  gkyl_array_set_range(pressure_tensor, 1.0, M2ij, conf_local);

  // Subtract V_i*M1_j and V_j*M1_i
  gkyl_array_accumulate_range(pressure_tensor, -0.5, up->V_drift_M1i, conf_local);
  gkyl_array_accumulate_range(pressure_tensor, -0.5, up->V_drift_M1i_transpose, conf_local);
  gkyl_array_scale_range(pressure_tensor, up->mass, conf_local);
}

void
gkyl_dg_vlasov_pressure_release(struct gkyl_dg_vlasov_pressure *up)
{
  gkyl_array_release(up->V_drift_dot_M1i);
  gkyl_array_release(up->V_drift_M1i);
  gkyl_array_release(up->V_drift_M1i_transpose);
  gkyl_free(up);
}
