#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_gk_particle_flux_parallel_priv.h>

struct gkyl_gk_particle_flux_parallel *
gkyl_gk_particle_flux_parallel_inew(const struct gkyl_gk_particle_flux_parallel_inp *inp)
{
  int cdim = inp->cbasis->ndim, vdim = inp->pbasis->ndim - cdim;

  if (inp->cbasis->poly_order != 1 || inp->pbasis->b_type != GKYL_BASIS_MODAL_GKHYBRID ||
      !((cdim == 1 && (vdim == 1 || vdim == 2)) || ((cdim == 2 || cdim == 3) && vdim == 2))) {
    gkyl_exit("Parallel flux requires a supported p1 gyrokinetic basis.");
  }

  struct gkyl_gk_particle_flux_parallel *up = gkyl_malloc(sizeof(*up));

  *up = (struct gkyl_gk_particle_flux_parallel){
    .cdim = cdim,
    .vdim = vdim,
    .nq = 1 << (cdim - 1),
    .nv = vdim == 1 ? 3 : 6,
    .edge = inp->edge,
    .use_gpu = inp->use_gpu,
    .charge = inp->charge,
    .mass = inp->mass,
    .grid = *inp->phase_grid,
    .conf_skin = *inp->conf_skin,
    .phase_skin = *inp->phase_skin,
    .gk_geom = gkyl_gk_geometry_acquire(inp->gk_geom),
    .dg_geom = gkyl_dg_geom_acquire(inp->dg_geom),
    .gk_dg_geom = gkyl_gk_dg_geom_acquire(inp->gk_dg_geom),
    .vel_map = gkyl_velocity_map_acquire(inp->vel_map),
  };

  gkyl_range_init(&up->compact, cdim, inp->conf_skin->lower, inp->conf_skin->upper);
  up->factor = vdim == 2 ? 2.0 * M_PI / inp->mass : 1.0;
  for (int d = 0; d < vdim; ++d) {
    up->factor *= inp->phase_grid->dx[cdim + d] / 2.0;
  }

  up->on_dev = up;
  if (up->use_gpu) {
#ifdef GKYL_HAVE_CUDA
    up->surf_geo_fields = gkyl_array_cu_dev_new(
      GKYL_USER, up->nq * sizeof(struct gkyl_gk_dg_surf_geom), up->compact.volume
    );
    up->area = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->nq, up->compact.volume);
    particle_flux_parallel_cu_init(up);
#else
    gkyl_exit("CUDA parallel flux requested in a CPU build.");
#endif
  } else {
    up->surf_geo_fields =
      gkyl_array_new(GKYL_USER, up->nq * sizeof(struct gkyl_gk_dg_surf_geom), up->compact.volume);
    up->area = gkyl_array_new(GKYL_DOUBLE, up->nq, up->compact.volume);
    up->kernel = inp->edge == GKYL_LOWER_EDGE ? choose_gk_collisionless_flux_surf_conf_kern(
                                                  cdim - 1, cdim, vdim, 1, GKYL_BC_GK_SPECIES_ABSORB
                                                ) :
                                                choose_gk_collisionless_flux_edge_surf_conf_kern(
                                                  cdim - 1, cdim, vdim, 1, GKYL_BC_GK_SPECIES_ABSORB
                                                );
    up->integral = particle_flux_parallel_integral_kernels[cdim - 1][vdim - 1];
    for (long i = 0; i < up->compact.volume; ++i) {
      particle_flux_parallel_init_geometry(up, i);
    }
  }

  return up;
}

void
gkyl_gk_particle_flux_parallel_advance(
  struct gkyl_gk_particle_flux_parallel *up, const struct gkyl_array *fin, struct gkyl_array *flux
)
{
  gkyl_array_clear(flux, 0.0);
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    return particle_flux_parallel_advance_cu(up, fin, flux);
  }
#endif
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->phase_skin);
  while (gkyl_range_iter_next(&iter)) {
    double val[4];
    particle_flux_parallel_cell(up, iter.idx, fin, val);

    double *out = gkyl_array_fetch(flux, gkyl_range_idx(&up->compact, iter.idx));
    for (int q = 0; q < up->nq; ++q) {
      out[q] += val[q];
    }
  }
}

const struct gkyl_array *
gkyl_gk_particle_flux_parallel_area(const struct gkyl_gk_particle_flux_parallel *up)
{
  return up->area;
}

void
gkyl_gk_particle_flux_parallel_release(struct gkyl_gk_particle_flux_parallel *up)
{
  gkyl_array_release(up->surf_geo_fields);
  gkyl_array_release(up->area);
  gkyl_gk_geometry_release(up->gk_geom);
  gkyl_dg_geom_release(up->dg_geom);
  gkyl_gk_dg_geom_release(up->gk_dg_geom);
  gkyl_velocity_map_release(up->vel_map);
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_cu_free(up->on_dev);
  }
#endif
  gkyl_free(up);
}
