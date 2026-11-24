#include <gkyl_alloc.h>
#include <gkyl_mom_type.h>
#include <gkyl_mom_gyrokinetic.h>
#include <gkyl_mom_gyrokinetic_priv.h>
#include <gkyl_util.h>

struct gkyl_mom_gyrokinetic *
gkyl_mom_gyrokinetic_new(
  double mass, double charge, const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis,
  const struct gkyl_rect_grid *phase_grid, const struct gkyl_velocity_map *vel_map,
  const struct gk_geometry *gk_geom, enum gkyl_distribution_moments mom_type, bool is_integrated,
  bool use_gpu
)
{
#ifndef GKYL_HAVE_CUDA
  if (use_gpu) {
    gkyl_exit("CUDA gyrokinetic moments requested in a CPU build.");
  }
#endif
  if (pbasis->ndim - cbasis->ndim == 1 && mom_type == GKYL_F_MOMENT_M0M1M2PARM2PERP) {
    mom_type = GKYL_F_MOMENT_M0M1M2;
  }
  struct gkyl_mom_gyrokinetic *up = gkyl_malloc(sizeof(*up));

  up->use_gpu = use_gpu;
  up->phase_grid = *phase_grid;
  up->mass = mass;
  up->charge = charge;
  up->gk_geom = gkyl_gk_geometry_acquire(gk_geom);
  up->vel_map = gkyl_velocity_map_acquire(vel_map);

  // Choose kernels that compute the moment and determine the number of moments.
  if (!use_gpu) {
    up->kernels = gkyl_malloc(sizeof(struct gkyl_mom_gyrokinetic_kernels));
    mom_gk_choose_kernel(up->kernels, &up->num_mom, cbasis, pbasis, mom_type, is_integrated);
  }
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    up->kernels = gkyl_cu_malloc(sizeof(struct gkyl_mom_gyrokinetic_kernels));
    mom_gk_choose_kernel_cu(up->kernels, &up->num_mom, cbasis, pbasis, mom_type, is_integrated);
  }
#endif

  return up;
}

int
gkyl_mom_gyrokinetic_num_mom(const gkyl_mom_gyrokinetic *up)
{
  return up->num_mom;
}

void
gkyl_mom_gyrokinetic_advance(
  struct gkyl_mom_gyrokinetic *up, const struct gkyl_range *phase_rng,
  const struct gkyl_range *conf_rng, const struct gkyl_array *GKYL_RESTRICT phi,
  const struct gkyl_array *GKYL_RESTRICT fin, struct gkyl_array *GKYL_RESTRICT mout
)
{
  if (conf_rng->volume == 0) {
    return;
  }

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_mom_gyrokinetic_advance_cu(up, phase_rng, conf_rng, phi, fin, mout);
    return;
  }
#endif

  struct gkyl_range vel_rng;
  struct gkyl_range_iter conf_iter, vel_iter;

  int rem_dir[GKYL_MAX_DIM] = {0};
  for (int d = 0; d < conf_rng->ndim; ++d) {
    rem_dir[d] = 1;
  }

  gkyl_array_clear_range(mout, 0.0, conf_rng);

  // the outer loop is over configuration space cells; for each
  // config-space cell the inner loop walks over the velocity space
  // computing the contribution to the moment
  gkyl_range_iter_init(&conf_iter, conf_rng);
  while (gkyl_range_iter_next(&conf_iter)) {
    long linidx_conf = gkyl_range_idx(conf_rng, conf_iter.idx);

    const double *bmag_c = gkyl_array_cfetch(up->gk_geom->geo_corn.bmag, linidx_conf);
    const double *phi_c = phi ? gkyl_array_cfetch(phi, linidx_conf) : 0;
    double *mom_c = gkyl_array_fetch(mout, linidx_conf);

    gkyl_range_deflate(&vel_rng, phase_rng, rem_dir, conf_iter.idx);
    gkyl_range_iter_no_split_init(&vel_iter, &vel_rng);

    while (gkyl_range_iter_next(&vel_iter)) {
      long linidx_vel = gkyl_range_idx(&up->vel_map->local_vel, vel_iter.idx);
      long linidx_phase = gkyl_range_idx(&vel_rng, vel_iter.idx);

      const double *vmap_c = gkyl_array_cfetch(up->vel_map->vmap, linidx_vel);
      const double *f_c = gkyl_array_cfetch(fin, linidx_phase);

      up->kernels->wmom(up->phase_grid.dx, vmap_c, up->mass, up->charge, bmag_c, phi_c, f_c, mom_c);
    }
  }
}

void
gkyl_mom_gyrokinetic_release(struct gkyl_mom_gyrokinetic *up)
{
  // Release memory associated with this updater.
  gkyl_gk_geometry_release(up->gk_geom);
  gkyl_velocity_map_release(up->vel_map);
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
