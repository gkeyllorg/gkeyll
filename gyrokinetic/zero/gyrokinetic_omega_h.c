#include <assert.h>
#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_gyrokinetic_omega_h_priv.h>

struct gkyl_gyrokinetic_omega_h *
gkyl_gyrokinetic_omega_h_new(
  const struct gkyl_basis *basis, const struct gkyl_range *range,
  const struct gkyl_array *geo_factor, bool use_gpu
)
{
  assert(basis->num_quad > 0 && geo_factor->ncomp == basis->num_basis);
  assert(use_gpu == gkyl_array_is_cu_dev(geo_factor));
  struct gkyl_gyrokinetic_omega_h *up = gkyl_malloc(sizeof(*up));
  up->range = *range;
  up->use_gpu = use_gpu;
  up->max_dev = NULL;
  int num_samples = basis->num_quad + 1;
  struct gkyl_array *basis_at_nodes = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, num_samples);
  struct gkyl_array *weight_sq = gkyl_array_new(GKYL_DOUBLE, num_samples, geo_factor->size);

  // Evaluate each basis mode at the basis's quadrature nodes and at the center.
  double *unit_mode = gkyl_calloc(basis->num_basis, sizeof(double));
  double *nodal = gkyl_malloc(basis->num_quad * sizeof(double));
  for (int mode = 0; mode < basis->num_basis; ++mode) {
    unit_mode[mode] = 1.0;
    for (int node = 0; node < basis->num_quad; ++node) {
      basis->modal_to_quad_nodal(unit_mode, nodal, node);
      double *basis_val = gkyl_array_fetch(basis_at_nodes, node);
      basis_val[mode] = nodal[node];
    }
    unit_mode[mode] = 0.0;
  }
  double center[GKYL_MAX_DIM] = {0.0};
  basis->eval(center, gkyl_array_fetch(basis_at_nodes, basis->num_quad));
  gkyl_free(unit_mode);
  gkyl_free(nodal);

  // Initialization alone transfers the modal geometry to the host for caching.
  struct gkyl_array *geo_ho = use_gpu ?
                                gkyl_array_new(GKYL_DOUBLE, geo_factor->ncomp, geo_factor->size) :
                                gkyl_array_acquire(geo_factor);
  if (use_gpu) {
    gkyl_array_copy(geo_ho, geo_factor);
  }
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);
    const double *geo = gkyl_array_cfetch(geo_ho, loc);
    double *weight = gkyl_array_fetch(weight_sq, loc);
    for (int node = 0; node < num_samples; ++node) {
      const double *basis_val = gkyl_array_cfetch(basis_at_nodes, node);
      double value = 0.0;
      for (int mode = 0; mode < basis->num_basis; ++mode) {
        value += geo[mode] * basis_val[mode];
      }
      weight[node] = value * value;
    }
  }
  gkyl_array_release(geo_ho);

  up->basis_at_nodes = basis_at_nodes;
  up->weight_sq = weight_sq;
  if (use_gpu) {
    up->basis_at_nodes =
      gkyl_array_cu_dev_new(GKYL_DOUBLE, basis_at_nodes->ncomp, basis_at_nodes->size);
    up->weight_sq = gkyl_array_cu_dev_new(GKYL_DOUBLE, weight_sq->ncomp, weight_sq->size);
    gkyl_array_copy(up->basis_at_nodes, basis_at_nodes);
    gkyl_array_copy(up->weight_sq, weight_sq);
    gkyl_array_release(basis_at_nodes);
    gkyl_array_release(weight_sq);
    up->max_dev = gkyl_cu_malloc(sizeof(double));
  }
  return up;
}

double
gkyl_gyrokinetic_omega_h_advance(
  struct gkyl_gyrokinetic_omega_h *up, const struct gkyl_array *density
)
{
  assert(density->ncomp == up->basis_at_nodes->ncomp);
  assert(density->size == up->weight_sq->size);
  assert(up->use_gpu == gkyl_array_is_cu_dev(density));
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    return gkyl_gyrokinetic_omega_h_advance_cu(up, density);
  }
#endif
  double max_value = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&up->range, iter.idx);
    max_value = fmax(
      max_value,
      gyrokinetic_omega_h_cell_max(
        up->basis_at_nodes, gkyl_array_cfetch(up->weight_sq, loc), gkyl_array_cfetch(density, loc)
      )
    );
  }
  return max_value;
}

void
gkyl_gyrokinetic_omega_h_release(struct gkyl_gyrokinetic_omega_h *up)
{
  if (!up) {
    return;
  }
  gkyl_array_release(up->basis_at_nodes);
  gkyl_array_release(up->weight_sq);
  if (up->use_gpu) {
    gkyl_cu_free(up->max_dev);
  }
  gkyl_free(up);
}
