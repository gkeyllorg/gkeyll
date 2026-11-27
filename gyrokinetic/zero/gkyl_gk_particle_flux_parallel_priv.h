#pragma once

#include <gkyl_gk_particle_flux_parallel.h>
#include <gkyl_gk_collisionless_flux.h>
#include <gkyl_gk_collisionless_flux_priv.h>
#include <gkyl_gk_particle_flux_parallel_kernels.h>

typedef void (*gk_particle_flux_parallel_integral_t)(double factor, const double *flux, double *out);

// Velocity-integral kernels indexed by cdim-1 and vdim-1 (p1 hybrid basis).
GKYL_CU_D static const gk_particle_flux_parallel_integral_t
  particle_flux_parallel_integral_kernels[3][2] = {
    {gk_particle_flux_parallel_integral_1x1v_ser_p1, gk_particle_flux_parallel_integral_1x2v_ser_p1
    },
    {NULL, gk_particle_flux_parallel_integral_2x2v_ser_p1},
    {NULL, gk_particle_flux_parallel_integral_3x2v_ser_p1}
};

struct gkyl_gk_particle_flux_parallel {
  int cdim; // Configuration-space dimension.
  int vdim; // Velocity-space dimension.
  int nq; // Configuration-space quadrature nodes per face.
  int nv; // Velocity quadrature nodes per velocity cell.
  enum gkyl_edge_loc edge; // Lower or upper parallel boundary.
  double charge; // Species charge.
  double mass; // Species mass.
  double factor; // Velocity integration prefactor.
  struct gkyl_rect_grid grid; // Phase-space grid.
  struct gkyl_range conf_skin; // Configuration-space skin range.
  struct gkyl_range phase_skin; // Phase-space skin range.
  struct gkyl_range compact; // Contiguous indexing range for surface arrays.
  const struct gk_geometry *gk_geom; // Modal gyrokinetic geometry.
  const struct gkyl_dg_geom *dg_geom; // Nodal metric geometry.
  const struct gkyl_gk_dg_geom *gk_dg_geom; // Nodal magnetic geometry.
  const struct gkyl_velocity_map *vel_map; // Computational-to-physical velocity map.
  struct gkyl_array *surf_geo_fields; // Surface gkyl_gk_dg_surf_geom records with drift terms zeroed.
  struct gkyl_array *area; // Nodal face area times absolute magnetic incidence.
  gk_collisionless_flux_surf_t kernel; // Absorbing-boundary numerical flux kernel.
  gk_particle_flux_parallel_integral_t
    integral; // Exact velocity integral of the nodal flux interpolant.
  struct gkyl_gk_particle_flux_parallel *on_dev; // Device copy, or self for CPU execution.
  bool use_gpu; // Run on the GPU.
};

GKYL_CU_DH static inline void
particle_flux_parallel_init_geometry(struct gkyl_gk_particle_flux_parallel *up, long linc)
{
  int idx[GKYL_MAX_DIM];
  gkyl_range_inv_idx(&up->compact, linc, idx);
  if (up->edge == GKYL_UPPER_EDGE) {
    idx[up->cdim - 1]++;
  }
  const struct gkyl_dg_surf_geom *dg = gkyl_dg_geom_get_surf(up->dg_geom, up->cdim - 1, idx);
  const struct gkyl_gk_dg_surf_geom *gk =
    gkyl_gk_dg_geom_get_surf(up->gk_dg_geom, up->cdim - 1, idx);
  struct gkyl_gk_dg_surf_geom *copy =
    (struct gkyl_gk_dg_surf_geom *)gkyl_array_fetch(up->surf_geo_fields, linc);
  double *area = (double *)gkyl_array_fetch(up->area, linc);
  for (int q = 0; q < up->nq; ++q) {
    copy[q] = gk[q];
    copy[q].bhat.x[0] = copy[q].bhat.x[1] = 0.0;
    copy[q].normcurlbhat = 0.0;
    area[q] = dg[q].area_elem * fabs(gk[q].B3 / gk[q].bmag);
  }
}

// Integrate one velocity cell's existing numerical surface flux.
GKYL_CU_DH static inline void
particle_flux_parallel_cell(
  const struct gkyl_gk_particle_flux_parallel *up, const int *pidx, const struct gkyl_array *fin,
  double *out
)
{
  int vidx[GKYL_MAX_VDIM], face[GKYL_MAX_DIM];
  gkyl_copy_int_arr(up->cdim, pidx, face);
  if (up->edge == GKYL_UPPER_EDGE) {
    face[up->cdim - 1]++;
  }
  for (int d = 0; d < up->vdim; ++d) {
    vidx[d] = pidx[up->cdim + d];
  }

  long cell = gkyl_range_idx(&up->compact, pidx);
  long vel = gkyl_range_idx(&up->vel_map->local_vel, vidx);

  const double *vmap = (const double *)gkyl_array_cfetch(up->vel_map->vmap, vel);
  const double *vmap_sq = (const double *)gkyl_array_cfetch(up->vel_map->vmap_sq, vel);
  const double *bmag = (const double *)gkyl_array_cfetch(
    up->gk_geom->geo_corn.bmag, gkyl_range_idx(&up->conf_skin, pidx)
  );
  const double *skin =
    (const double *)gkyl_array_cfetch(fin, gkyl_range_idx(&up->phase_skin, pidx));
  const struct gkyl_dg_surf_geom *dg = gkyl_dg_geom_get_surf(up->dg_geom, up->cdim - 1, face);
  const struct gkyl_gk_dg_surf_geom *gk =
    (const struct gkyl_gk_dg_surf_geom *)gkyl_array_cfetch(up->surf_geo_fields, cell);

  double zero[48] = {0.0}, nodal[72] = {0.0}, xc[GKYL_MAX_DIM];
  gkyl_rect_grid_cell_center(&up->grid, pidx, xc);

  bool lower = up->edge == GKYL_LOWER_EDGE;
  up->kernel(
    xc, up->grid.dx, vmap, vmap_sq, up->charge, up->mass, dg, gk, bmag, zero, zero, zero,
    lower ? zero : skin, lower ? skin : zero, nodal
  );

  int offset = (up->cdim - 1) * up->nq * up->nv;
  up->integral((lower ? -1.0 : 1.0) * up->factor, nodal + offset, out);
}

#ifdef GKYL_HAVE_CUDA
void particle_flux_parallel_cu_init(struct gkyl_gk_particle_flux_parallel *up);
void particle_flux_parallel_advance_cu(
  struct gkyl_gk_particle_flux_parallel *up, const struct gkyl_array *fin, struct gkyl_array *flux
);
#endif
