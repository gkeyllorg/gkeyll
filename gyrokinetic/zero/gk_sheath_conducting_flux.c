#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_gk_sheath_conducting_flux.h>
#include <gkyl_gk_sheath_conducting_flux_priv.h>

#include <assert.h>
#include <string.h>

static bool
gk_sheath_conducting_flux_eligible_flux(
  const struct gkyl_gk_sheath_conducting_flux *up, const double *flux, const double *alpha_out,
  const double *alpha_in, double *eligible_flux
)
{
  int face_basis = up->surf_basis->num_basis;
  double edge_sign = up->edge == GKYL_LOWER_EDGE ? -1.0 : 1.0;
  double flux_quad[face_basis], alpha_in_modal[face_basis], alpha_in_reflected[face_basis];
  double alpha_in_quad[face_basis], eligible_modal[face_basis];

  for (int node = 0; node < face_basis; ++node) {
    up->surf_basis->modal_to_quad_nodal(flux, flux_quad, node);
    up->surf_basis->quad_nodal_to_modal(alpha_in, alpha_in_modal, node);
  }
  up->surf_basis->flip_odd_sign(up->cdim - 1, alpha_in_modal, alpha_in_reflected);
  for (int node = 0; node < face_basis; ++node) {
    up->surf_basis->modal_to_quad_nodal(alpha_in_reflected, alpha_in_quad, node);
  }

  bool has_eligible_flux = false;
  for (int node = 0; node < face_basis; ++node) {
    bool outgoing = edge_sign * alpha_out[node] > 0.0;
    bool incoming_partner = edge_sign * alpha_in_quad[node] < 0.0;
    eligible_flux[node] = outgoing && incoming_partner ? flux_quad[node] : 0.0;
    has_eligible_flux = has_eligible_flux || (outgoing && incoming_partner);
  }
  for (int basis_idx = 0; basis_idx < face_basis; ++basis_idx) {
    up->surf_basis->quad_nodal_to_modal(eligible_flux, eligible_modal, basis_idx);
  }
  memcpy(eligible_flux, eligible_modal, sizeof(eligible_modal));

  return has_eligible_flux;
}

void
gkyl_gk_sheath_conducting_flux_advance(
  gkyl_gk_sheath_conducting_flux *up, const struct gkyl_array *phi,
  const struct gkyl_array *gyro_phi, const struct gkyl_array *phi_wall,
  const struct gkyl_array *fin, struct gkyl_array *flux_surf
)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_gk_sheath_conducting_flux_advance_cu(up, phi, gyro_phi, phi_wall, fin, flux_surf);
    return;
  }
#endif
  const struct gkyl_gk_collisionless_flux *flux_op = up->flux_op;
  int pdim = up->phase_basis->ndim;
  int vpar_dir = up->cdim;
  int face_basis = up->surf_basis->num_basis;
  int flux_ncomp = flux_surf->ncomp;
  int vpar_sum = up->phase_ext_range.lower[vpar_dir] + up->phase_ext_range.upper[vpar_dir];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->skin_range);
  while (gkyl_range_iter_next(&iter)) {
    int face_idx[GKYL_MAX_DIM], partner_idx[GKYL_MAX_DIM], vel_idx[2], partner_vel_idx[2];
    gkyl_copy_int_arr(pdim, iter.idx, face_idx);
    if (up->edge == GKYL_UPPER_EDGE) {
      face_idx[up->dir] += 1;
    }
    gkyl_copy_int_arr(pdim, face_idx, partner_idx);
    partner_idx[vpar_dir] = vpar_sum - iter.idx[vpar_dir];

    for (int d = up->cdim; d < pdim; ++d) {
      vel_idx[d - up->cdim] = iter.idx[d];
      partner_vel_idx[d - up->cdim] = partner_idx[d];
    }

    long conf_loc = gkyl_range_idx(&up->conf_range, iter.idx);
    long vel_loc = gkyl_range_idx(&flux_op->vel_map->local_vel, vel_idx);
    long partner_vel_loc = gkyl_range_idx(&flux_op->vel_map->local_vel, partner_vel_idx);
    long skin_loc = gkyl_range_idx(&up->phase_ext_range, iter.idx);
    long partner_loc = gkyl_range_idx(&up->phase_ext_range, partner_idx);

    int geom_idx[GKYL_MAX_DIM];
    gkyl_copy_int_arr(pdim, iter.idx, geom_idx);
    if (up->edge == GKYL_UPPER_EDGE) {
      geom_idx[up->dir] += 1;
    }
    long conf_face_loc = gkyl_range_idx(&up->conf_ext_range, geom_idx);

    const double *f = gkyl_array_cfetch(fin, skin_loc);
    const double *phi_cut = gkyl_array_cfetch(phi, conf_loc);
    const double *phi = gkyl_array_cfetch(gyro_phi, conf_loc);
    const double *wall = gkyl_array_cfetch(phi_wall, conf_loc);
    const double *bmag = gkyl_array_cfetch(flux_op->gk_geom->geo_corn.bmag, conf_loc);
    const double *vmap = gkyl_array_cfetch(flux_op->vel_map->vmap, vel_loc);
    const double *vmap_sq = gkyl_array_cfetch(flux_op->vel_map->vmap_sq, vel_loc);
    const double *partner_vmap = gkyl_array_cfetch(flux_op->vel_map->vmap, partner_vel_loc);
    const double *partner_vmap_sq = gkyl_array_cfetch(flux_op->vel_map->vmap_sq, partner_vel_loc);
    const struct gkyl_dg_surf_geom *dgs =
      gkyl_dg_geom_get_surf(flux_op->dg_geom, up->dir, geom_idx);
    const struct gkyl_gk_dg_surf_geom *gkdgs =
      gkyl_gk_dg_geom_get_surf(flux_op->gk_dg_geom, up->dir, geom_idx);
    const double *jac_l;
    const double *jac_r;
    if (up->edge == GKYL_LOWER_EDGE) {
      int left_idx[GKYL_MAX_DIM];
      gkyl_copy_int_arr(pdim, iter.idx, left_idx);
      left_idx[up->dir] -= 1;
      jac_l = gkyl_array_cfetch(
        flux_op->gk_geom->geo_surf[up->dir].jacobgeo_ratio,
        gkyl_range_idx(&up->conf_ext_range, left_idx)
      );
      jac_r = gkyl_array_cfetch(flux_op->gk_geom->geo_surf[up->dir].jacobgeo_ratio, conf_loc);
    } else {
      jac_l = gkyl_array_cfetch(flux_op->gk_geom->geo_surf[up->dir].jacobgeo_ratio, conf_loc);
      jac_r = gkyl_array_cfetch(flux_op->gk_geom->geo_surf[up->dir].jacobgeo_ratio, conf_face_loc);
    }

    double fhat[up->phase_basis->num_basis];
    double zero[up->phase_basis->num_basis];
    double low_flux[flux_ncomp];
    double alpha_out[face_basis];
    double alpha_in[face_basis];
    double xc[GKYL_MAX_DIM], partner_xc[GKYL_MAX_DIM];
    memset(zero, 0, sizeof(zero));
    memset(low_flux, 0, sizeof(low_flux));
    memset(alpha_out, 0, sizeof(alpha_out));
    memset(alpha_in, 0, sizeof(alpha_in));
    gkyl_rect_grid_cell_center(&up->phase_grid, iter.idx, xc);
    gkyl_rect_grid_cell_center(&up->phase_grid, partner_idx, partner_xc);

    up->reflectedf(vmap, up->q2Dm, phi_cut, wall, f, fhat);
    if (up->edge == GKYL_LOWER_EDGE) {
      flux_op->flux_surf_edge_lo[up->dir](
        xc, up->phase_grid.dx, vmap, vmap_sq, flux_op->charge, flux_op->mass, dgs, gkdgs, bmag,
        jac_l, jac_r, phi, zero, fhat, low_flux
      );
    } else {
      flux_op->flux_surf_edge_up[up->dir](
        xc, up->phase_grid.dx, vmap, vmap_sq, flux_op->charge, flux_op->mass, dgs, gkdgs, bmag,
        jac_l, jac_r, phi, fhat, zero, low_flux
      );
    }
    gk_collisionless_flux_alpha_surf_t alpha_surf = up->edge == GKYL_LOWER_EDGE ?
                                                      flux_op->alpha_surf_lo[up->dir] :
                                                      flux_op->alpha_surf_up[up->dir];
    alpha_surf(
      xc, up->phase_grid.dx, vmap, vmap_sq, flux_op->charge, flux_op->mass, dgs, gkdgs, bmag, phi,
      alpha_out
    );
    alpha_surf(
      partner_xc, up->phase_grid.dx, partner_vmap, partner_vmap_sq, flux_op->charge, flux_op->mass,
      dgs, gkdgs, bmag, phi, alpha_in
    );

    double eligible_flux[face_basis];
    if (gk_sheath_conducting_flux_eligible_flux(
          up, low_flux + up->dir * face_basis, alpha_out, alpha_in, eligible_flux
        )) {
      double reflected_flux[face_basis];
      double *target_flux = gkyl_array_fetch(flux_surf, partner_loc) + up->dir * face_basis;
      up->surf_basis->flip_odd_sign(up->cdim - 1, eligible_flux, reflected_flux);
      for (int basis_idx = 0; basis_idx < face_basis; ++basis_idx) {
        target_flux[basis_idx] -= reflected_flux[basis_idx];
      }
    }
  }
}

gkyl_gk_sheath_conducting_flux *
gkyl_gk_sheath_conducting_flux_inew(const struct gkyl_gk_sheath_conducting_flux_inp *inp)
{
  assert(inp->phase_basis->poly_order == 1);
  assert((int)inp->phase_basis->ndim == inp->cdim + 2);
  assert(inp->dir == inp->cdim - 1);

  gkyl_gk_sheath_conducting_flux *up = gkyl_malloc(sizeof(*up));
  up->dir = inp->dir;
  up->cdim = inp->cdim;
  up->edge = inp->edge;
  up->q2Dm = 2.0 * inp->charge / inp->mass;
  up->use_gpu = inp->use_gpu;
  up->phase_grid = *inp->phase_grid;
  up->phase_basis = inp->phase_basis;
  up->surf_basis = gkyl_cart_modal_gkhybrid_new(inp->cdim - 1, 2);
  up->conf_range = *inp->conf_range;
  up->conf_ext_range = *inp->conf_ext_range;
  up->phase_ext_range = *inp->phase_ext_range;
  up->skin_range = *inp->skin_range;
  up->flux_op = inp->flux_op;
  up->reflectedf = bc_gksheath_choose_reflectedf_kernel(inp->phase_basis, inp->edge);
  assert(up->reflectedf);
  assert(up->flux_op->alpha_surf_lo[up->dir]);
  assert(up->flux_op->alpha_surf_up[up->dir]);
  up->on_dev = up;
#ifdef GKYL_HAVE_CUDA
  if (inp->use_gpu) {
    const struct gkyl_basis *phase_basis_cu = inp->phase_basis_cu;
    assert(phase_basis_cu);
    struct gkyl_basis *surf_basis_cu = gkyl_cart_modal_gkhybrid_cu_dev_new(inp->cdim - 1, 2);
    up->on_dev = gkyl_cu_malloc(sizeof(*up));
    gkyl_cu_memcpy(up->on_dev, up, sizeof(*up), GKYL_CU_MEMCPY_H2D);
    gkyl_cu_memcpy(
      &up->on_dev->phase_basis, &phase_basis_cu, sizeof(phase_basis_cu), GKYL_CU_MEMCPY_H2D
    );
    gkyl_cu_memcpy(
      &up->on_dev->surf_basis, &surf_basis_cu, sizeof(surf_basis_cu), GKYL_CU_MEMCPY_H2D
    );
    gkyl_cu_memcpy(
      &up->on_dev->flux_op, &inp->flux_op->on_dev, sizeof(inp->flux_op->on_dev), GKYL_CU_MEMCPY_H2D
    );
    gkyl_gk_sheath_conducting_flux_set_cu_ptrs(up->on_dev, phase_basis_cu);
  }
#endif
  return up;
}

void
gkyl_gk_sheath_conducting_flux_release(gkyl_gk_sheath_conducting_flux *up)
{
  if (up) {
#ifdef GKYL_HAVE_CUDA
    if (up->use_gpu) {
      struct gkyl_basis *surf_basis_cu;
      gkyl_cu_memcpy(
        &surf_basis_cu, &up->on_dev->surf_basis, sizeof(surf_basis_cu), GKYL_CU_MEMCPY_D2H
      );
      gkyl_cu_free(surf_basis_cu);
      gkyl_cu_free(up->on_dev);
    }
#endif
    gkyl_free((void *)up->surf_basis);
    gkyl_free(up);
  }
}
