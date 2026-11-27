extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_gk_particle_flux_parallel_priv.h>
}

__global__ static void
particle_flux_parallel_init_cu_ker(struct gkyl_gk_particle_flux_parallel *up)
{
  if (threadIdx.x == 0 && blockIdx.x == 0) {
    up->integral = particle_flux_parallel_integral_kernels[up->cdim - 1][up->vdim - 1];
    up->kernel = up->edge == GKYL_LOWER_EDGE ?
                   choose_gk_collisionless_flux_surf_conf_kern(
                     up->cdim - 1, up->cdim, up->vdim, 1, GKYL_BC_GK_SPECIES_ABSORB
                   ) :
                   choose_gk_collisionless_flux_edge_surf_conf_kern(
                     up->cdim - 1, up->cdim, up->vdim, 1, GKYL_BC_GK_SPECIES_ABSORB
                   );
  }
  for (long i = threadIdx.x + blockIdx.x * blockDim.x; i < up->compact.volume;
       i += blockDim.x * gridDim.x) {
    particle_flux_parallel_init_geometry(up, i);
  }
}

void
particle_flux_parallel_cu_init(struct gkyl_gk_particle_flux_parallel *up)
{
  struct gkyl_gk_particle_flux_parallel device = *up;
  device.gk_geom = up->gk_geom->on_dev;
  device.dg_geom = up->dg_geom->on_dev;
  device.gk_dg_geom = up->gk_dg_geom->on_dev;
  device.vel_map = up->vel_map->on_dev;
  device.surf_geo_fields = up->surf_geo_fields->on_dev;
  device.area = up->area->on_dev;
  up->on_dev = (struct gkyl_gk_particle_flux_parallel *)gkyl_cu_malloc(sizeof(*up));
  device.on_dev = up->on_dev;
  gkyl_cu_memcpy(up->on_dev, &device, sizeof(device), GKYL_CU_MEMCPY_H2D);
  particle_flux_parallel_init_cu_ker<<<up->compact.nblocks, up->compact.nthreads>>>(up->on_dev);
}

__global__ static void
particle_flux_parallel_advance_cu_ker(
  const struct gkyl_gk_particle_flux_parallel *up, const struct gkyl_array *fin,
  struct gkyl_array *flux
)
{
  int pidx[GKYL_MAX_DIM];
  for (long i = threadIdx.x + blockIdx.x * blockDim.x; i < up->phase_skin.volume;
       i += blockDim.x * gridDim.x) {
    gkyl_sub_range_inv_idx(&up->phase_skin, i, pidx);
    double val[4];
    particle_flux_parallel_cell(up, pidx, fin, val);
    double *out = (double *)gkyl_array_fetch(flux, gkyl_range_idx(&up->compact, pidx));
    for (int q = 0; q < up->nq; ++q) {
      atomicAdd(out + q, val[q]);
    }
  }
}

void
particle_flux_parallel_advance_cu(
  struct gkyl_gk_particle_flux_parallel *up, const struct gkyl_array *fin, struct gkyl_array *flux
)
{
  particle_flux_parallel_advance_cu_ker<<<up->phase_skin.nblocks, up->phase_skin.nthreads>>>(
    up->on_dev, fin->on_dev, flux->on_dev
  );
}
