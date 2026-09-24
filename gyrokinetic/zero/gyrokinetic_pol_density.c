#include <gkyl_gyrokinetic_pol_density.h>
#include <gkyl_gyrokinetic_pol_density_priv.h>
#include <gkyl_alloc.h>

struct gkyl_gyrokinetic_pol_density*
gkyl_gyrokinetic_pol_density_new(struct gkyl_basis cbasis, struct gkyl_rect_grid cgrid,
  bool use_gpu)
{
  // Allocate space for new updater.
  struct gkyl_gyrokinetic_pol_density *up = gkyl_malloc(sizeof(*up));

  up->grid = cgrid;
  up->use_gpu = use_gpu;

  if (!use_gpu) {
    up->kernels = gkyl_malloc(sizeof(struct gkyl_gyrokinetic_pol_density_kernels));
  }

#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    up->kernels = gkyl_cu_malloc(sizeof(struct gkyl_gyrokinetic_pol_density_kernels));
  }
#endif

  // Choose kernels that shift f and compute int moms of Deltaf.
  gk_pol_den_choose_kernel(up->kernels, cbasis, use_gpu);

  return up;
}

void
gkyl_gyrokinetic_pol_density_advance(gkyl_gyrokinetic_pol_density* up,
  const struct gkyl_range *conf_rng, const struct gkyl_array *GKYL_RESTRICT pol_weight,
  const struct gkyl_array *GKYL_RESTRICT phi, struct gkyl_array *GKYL_RESTRICT npol)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_gyrokinetic_pol_density_advance_cu(up, conf_rng, pol_weight, phi, npol);
    return;
  }
#endif

  struct gkyl_range_iter conf_iter;
  gkyl_range_iter_init(&conf_iter, conf_rng);
  while (gkyl_range_iter_next(&conf_iter)) {
    long linidx = gkyl_range_idx(conf_rng, conf_iter.idx);

    const double *pol_weight_d = gkyl_array_cfetch(pol_weight, linidx);
    const double *phi_d = gkyl_array_cfetch(phi, linidx);
    double *npol_d = gkyl_array_fetch(npol, linidx);

    up->kernels->pol_den(up->grid.dx, pol_weight_d, phi_d, npol_d);
  }

}

void
gkyl_gyrokinetic_pol_density_release(gkyl_gyrokinetic_pol_density* up)
{
  // Release memory associated with this updater.
  if (!up->use_gpu) {
    gkyl_free(up->kernels);
  }
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_cu_free(up->kernels);
  }
#endif
  gkyl_free(up);
}
