#include <assert.h>
#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_gyrokinetic_omega_h_priv.h>

struct gkyl_gyrokinetic_omega_h *
gkyl_gyrokinetic_omega_h_new(
  const struct gkyl_basis *basis, const struct gkyl_range *range,
  const struct gkyl_array *omegaH_geo_fac, bool use_gpu
)
{
  assert(basis->num_quad > 0 && omegaH_geo_fac->ncomp == basis->num_basis);
  assert(use_gpu == gkyl_array_is_cu_dev(omegaH_geo_fac));
  struct gkyl_gyrokinetic_omega_h *up = gkyl_malloc(sizeof(*up));
  up->range = *range;
  up->use_gpu = use_gpu;
  up->max_dev = NULL;
  int num_samples = basis->num_quad + 1;
  struct gkyl_array *basis_at_nodes = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, num_samples);
  struct gkyl_array *omegaH_geo_fac_sq =
    gkyl_array_new(GKYL_DOUBLE, num_samples, omegaH_geo_fac->size);

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
  struct gkyl_array *omegaH_geo_fac_ho =
    use_gpu ? gkyl_array_new(GKYL_DOUBLE, omegaH_geo_fac->ncomp, omegaH_geo_fac->size) :
              gkyl_array_acquire(omegaH_geo_fac);
  if (use_gpu) {
    gkyl_array_copy(omegaH_geo_fac_ho, omegaH_geo_fac);
  }
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);
    const double *omegaH_geo_fac_d = gkyl_array_cfetch(omegaH_geo_fac_ho, loc);
    double *omegaH_geo_fac_sq_d = gkyl_array_fetch(omegaH_geo_fac_sq, loc);
    for (int node = 0; node < num_samples; ++node) {
      const double *basis_val = gkyl_array_cfetch(basis_at_nodes, node);
      double value = 0.0;
      for (int mode = 0; mode < basis->num_basis; ++mode) {
        value += omegaH_geo_fac_d[mode] * basis_val[mode];
      }
      omegaH_geo_fac_sq_d[node] = value * value;
    }
  }
  gkyl_array_release(omegaH_geo_fac_ho);

  up->basis_at_nodes = basis_at_nodes;
  up->omegaH_geo_fac_sq = omegaH_geo_fac_sq;
  if (use_gpu) {
    up->basis_at_nodes =
      gkyl_array_cu_dev_new(GKYL_DOUBLE, basis_at_nodes->ncomp, basis_at_nodes->size);
    up->omegaH_geo_fac_sq =
      gkyl_array_cu_dev_new(GKYL_DOUBLE, omegaH_geo_fac_sq->ncomp, omegaH_geo_fac_sq->size);
    gkyl_array_copy(up->basis_at_nodes, basis_at_nodes);
    gkyl_array_copy(up->omegaH_geo_fac_sq, omegaH_geo_fac_sq);
    gkyl_array_release(basis_at_nodes);
    gkyl_array_release(omegaH_geo_fac_sq);
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
  assert(density->size == up->omegaH_geo_fac_sq->size);
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
      max_value, gyrokinetic_omega_h_cell_max(
                   up->basis_at_nodes, gkyl_array_cfetch(up->omegaH_geo_fac_sq, loc),
                   gkyl_array_cfetch(density, loc)
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
  gkyl_array_release(up->omegaH_geo_fac_sq);
  if (up->use_gpu) {
    gkyl_cu_free(up->max_dev);
  }
  gkyl_free(up);
}
