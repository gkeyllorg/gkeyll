/* -*- c++ -*- */

extern "C" {
#include <gkyl_gk_sheath_conducting_flux.h>
#include <gkyl_gk_sheath_conducting_flux_priv.h>
}

enum {
  GKYL_GK_SHEATH_CONDUCTING_FLUX_MAX_BASIS = 48,
  GKYL_GK_SHEATH_CONDUCTING_FLUX_MAX_NCOMP = 128
};

GKYL_CU_D static bool
gk_sheath_conducting_flux_eligible_flux_cu(
  const struct gkyl_gk_sheath_conducting_flux *up, const double *flux, const double *alpha_out,
  const double *alpha_in, double *eligible_flux
)
{
  int face_basis = up->surf_basis->num_basis;
  double edge_sign = up->edge == GKYL_LOWER_EDGE ? -1.0 : 1.0;
  double flux_quad[GKYL_GK_SHEATH_CONDUCTING_FLUX_MAX_BASIS];
  double alpha_in_modal[GKYL_GK_SHEATH_CONDUCTING_FLUX_MAX_BASIS];
  double alpha_in_reflected[GKYL_GK_SHEATH_CONDUCTING_FLUX_MAX_BASIS];
  double alpha_in_quad[GKYL_GK_SHEATH_CONDUCTING_FLUX_MAX_BASIS];
  double eligible_modal[GKYL_GK_SHEATH_CONDUCTING_FLUX_MAX_BASIS];

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
  for (int basis_idx = 0; basis_idx < face_basis; ++basis_idx) {
    eligible_flux[basis_idx] = eligible_modal[basis_idx];
  }

  return has_eligible_flux;
}

__global__ static void
gkyl_gk_sheath_conducting_flux_advance_cu_kernel(
  const struct gkyl_gk_sheath_conducting_flux *up, const struct gkyl_array *phi,
  const struct gkyl_array *gyro_phi, const struct gkyl_array *phi_wall,
  const struct gkyl_array *fin, struct gkyl_array *flux_surf
)
{
  const struct gkyl_gk_collisionless_flux *flux_op = up->flux_op;
  int pdim = up->phase_basis->ndim;
  int vpar_dir = up->cdim;
  int face_basis = up->surf_basis->num_basis;
  int vpar_sum = up->phase_ext_range.lower[vpar_dir] + up->phase_ext_range.upper[vpar_dir];

  for (unsigned long linc = threadIdx.x + blockIdx.x * blockDim.x; linc < up->skin_range.volume;
       linc += blockDim.x * gridDim.x) {
    int face_idx[GKYL_MAX_DIM], partner_idx[GKYL_MAX_DIM], vel_idx[2], partner_vel_idx[2];
    int skin_idx[GKYL_MAX_DIM], geom_idx[GKYL_MAX_DIM];
    gkyl_sub_range_inv_idx(&up->skin_range, linc, skin_idx);
    gkyl_copy_int_arr(pdim, skin_idx, face_idx);
    if (up->edge == GKYL_UPPER_EDGE) {
      face_idx[up->dir] += 1;
    }
    gkyl_copy_int_arr(pdim, face_idx, partner_idx);
    partner_idx[vpar_dir] = vpar_sum - skin_idx[vpar_dir];

    for (int d = up->cdim; d < pdim; ++d) {
      vel_idx[d - up->cdim] = skin_idx[d];
      partner_vel_idx[d - up->cdim] = partner_idx[d];
    }

    long conf_loc = gkyl_range_idx(&up->conf_range, skin_idx);
    long vel_loc = gkyl_range_idx(&flux_op->vel_map->local_vel, vel_idx);
    long partner_vel_loc = gkyl_range_idx(&flux_op->vel_map->local_vel, partner_vel_idx);
    long skin_loc = gkyl_range_idx(&up->phase_ext_range, skin_idx);
    long partner_loc = gkyl_range_idx(&up->phase_ext_range, partner_idx);

    gkyl_copy_int_arr(pdim, skin_idx, geom_idx);
    if (up->edge == GKYL_UPPER_EDGE) {
      geom_idx[up->dir] += 1;
    }
    long conf_face_loc = gkyl_range_idx(&up->conf_ext_range, geom_idx);

    const double *f = (const double *)gkyl_array_cfetch(fin, skin_loc);
    const double *phi_cut = (const double *)gkyl_array_cfetch(phi, conf_loc);
    const double *phi = (const double *)gkyl_array_cfetch(gyro_phi, conf_loc);
    const double *wall = (const double *)gkyl_array_cfetch(phi_wall, conf_loc);
    const double *bmag =
      (const double *)gkyl_array_cfetch(flux_op->gk_geom->geo_corn.bmag, conf_loc);
    const double *vmap = (const double *)gkyl_array_cfetch(flux_op->vel_map->vmap, vel_loc);
    const double *vmap_sq = (const double *)gkyl_array_cfetch(flux_op->vel_map->vmap_sq, vel_loc);
    const double *partner_vmap =
      (const double *)gkyl_array_cfetch(flux_op->vel_map->vmap, partner_vel_loc);
    const double *partner_vmap_sq =
      (const double *)gkyl_array_cfetch(flux_op->vel_map->vmap_sq, partner_vel_loc);
    const struct gkyl_dg_surf_geom *dgs =
      gkyl_dg_geom_get_surf(flux_op->dg_geom, up->dir, geom_idx);
    const struct gkyl_gk_dg_surf_geom *gkdgs =
      gkyl_gk_dg_geom_get_surf(flux_op->gk_dg_geom, up->dir, geom_idx);
    const double *jac_l;
    const double *jac_r;
    if (up->edge == GKYL_LOWER_EDGE) {
      int left_idx[GKYL_MAX_DIM];
      gkyl_copy_int_arr(pdim, skin_idx, left_idx);
      left_idx[up->dir] -= 1;
      jac_l = (const double *)gkyl_array_cfetch(
        flux_op->gk_geom->geo_surf[up->dir].jacobgeo_ratio,
        gkyl_range_idx(&up->conf_ext_range, left_idx)
      );
      jac_r = (const double *)gkyl_array_cfetch(
        flux_op->gk_geom->geo_surf[up->dir].jacobgeo_ratio, conf_loc
      );
    } else {
      jac_l = (const double *)gkyl_array_cfetch(
        flux_op->gk_geom->geo_surf[up->dir].jacobgeo_ratio, conf_loc
      );
      jac_r = (const double *)gkyl_array_cfetch(
        flux_op->gk_geom->geo_surf[up->dir].jacobgeo_ratio, conf_face_loc
      );
    }

    double fhat[GKYL_GK_SHEATH_CONDUCTING_FLUX_MAX_BASIS] = {0.0};
    double zero[GKYL_GK_SHEATH_CONDUCTING_FLUX_MAX_BASIS] = {0.0};
    double low_flux[GKYL_GK_SHEATH_CONDUCTING_FLUX_MAX_NCOMP] = {0.0};
    double alpha_out[GKYL_GK_SHEATH_CONDUCTING_FLUX_MAX_BASIS] = {0.0};
    double alpha_in[GKYL_GK_SHEATH_CONDUCTING_FLUX_MAX_BASIS] = {0.0};
    double xc[GKYL_MAX_DIM], partner_xc[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(&up->phase_grid, skin_idx, xc);
    gkyl_rect_grid_cell_center(&up->phase_grid, partner_idx, partner_xc);

    up->reflectedf(vmap, up->q2Dm, phi_cut, wall, f, fhat);
    if (up->edge == GKYL_LOWER_EDGE) {
      up->flux_op->flux_surf_edge_lo[up->dir](
        xc, up->phase_grid.dx, vmap, vmap_sq, flux_op->charge, flux_op->mass, dgs, gkdgs, bmag,
        jac_l, jac_r, phi, zero, fhat, low_flux
      );
    } else {
      up->flux_op->flux_surf_edge_up[up->dir](
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

    double eligible_flux[GKYL_GK_SHEATH_CONDUCTING_FLUX_MAX_BASIS];
    if (gk_sheath_conducting_flux_eligible_flux_cu(
          up, low_flux + up->dir * face_basis, alpha_out, alpha_in, eligible_flux
        )) {
      double reflected_flux[GKYL_GK_SHEATH_CONDUCTING_FLUX_MAX_BASIS];
      double *target_flux =
        (double *)gkyl_array_fetch(flux_surf, partner_loc) + up->dir * face_basis;
      up->surf_basis->flip_odd_sign(up->cdim - 1, eligible_flux, reflected_flux);
      for (int basis_idx = 0; basis_idx < face_basis; ++basis_idx) {
        target_flux[basis_idx] -= reflected_flux[basis_idx];
      }
    }
  }
}

__global__ static void
gkyl_gk_sheath_conducting_flux_set_cu_ptrs_kernel(
  struct gkyl_gk_sheath_conducting_flux *up, const struct gkyl_basis *phase_basis
)
{
  up->phase_basis = phase_basis;
  up->reflectedf = bc_gksheath_choose_reflectedf_kernel(phase_basis, up->edge);
}

void
gkyl_gk_sheath_conducting_flux_set_cu_ptrs(
  struct gkyl_gk_sheath_conducting_flux *up, const struct gkyl_basis *phase_basis
)
{
  gkyl_gk_sheath_conducting_flux_set_cu_ptrs_kernel<<<1, 1>>>(up, phase_basis);
}

void
gkyl_gk_sheath_conducting_flux_advance_cu(
  const struct gkyl_gk_sheath_conducting_flux *up, const struct gkyl_array *phi,
  const struct gkyl_array *gyro_phi, const struct gkyl_array *phi_wall,
  const struct gkyl_array *fin, struct gkyl_array *flux_surf
)
{
  if (up->skin_range.volume > 0) {
    gkyl_gk_sheath_conducting_flux_advance_cu_kernel<<<
      up->skin_range.nblocks, up->skin_range.nthreads>>>(
      up->on_dev, phi->on_dev, gyro_phi->on_dev, phi_wall->on_dev, fin->on_dev, flux_surf->on_dev
    );
  }
}
