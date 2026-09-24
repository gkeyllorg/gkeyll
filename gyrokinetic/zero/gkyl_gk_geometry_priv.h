#include <gkyl_gk_geometry.h>
#include <gkyl_nodal_ops.h>

static inline void
get_filename_from_path(const char *filepath, char *out_buffer, size_t out_size)
{
  // Remove path from 'filepath' and output only the file name (including its extension)
  // to the buffer 'out_buffer' of size 'out_size'.
  const char *last_slash = strrchr(filepath, '/');
  const char *filename_start = (last_slash) ? last_slash + 1 : filepath;

  if (out_size > 0) {
    strncpy(out_buffer, filename_start, out_size - 1);
    out_buffer[out_size - 1] = '\0';
  }
}

static double calc_running_coord(double coord_lo, int i, double dx) {
  double dels[2] = {1.0/sqrt(3), 1.0-1.0/sqrt(3) };
  double coord = coord_lo;
  for(int j = 0; j < i; j++)
    coord+=dels[j%2]*dx;
  return coord;
}

static double calc_running_surf_coord(double coord_lo, int i, double dx) {
  double dels[3] = {(1.0-1.0/sqrt(3))/2.0, 1.0/sqrt(3), (1.0-1.0/sqrt(3))/2.0 };
  double coord = coord_lo;
  for(int j = 0; j < i; j++)
    coord+=dels[j%3]*dx;
  return coord;
}

static void gk_geometry_set_nodal_ranges(struct gk_geometry* up) 
{
  int poly_order = up->basis.poly_order;

  // nodes tensor
  int num_nodes_corners[GKYL_MAX_CDIM];
  if (poly_order == 1) {
    for (int d=0; d<up->grid.ndim; ++d)
      num_nodes_corners[d] = gkyl_range_shape(&up->local, d) + 1;
  }
  if (poly_order == 2) {
    for (int d=0; d<up->grid.ndim; ++d)
      num_nodes_corners[d] = 2*gkyl_range_shape(&up->local, d) + 1;
  }

  int num_quad_points = poly_order+1;

  int num_nodes_quad_interior[GKYL_MAX_CDIM];
  for (int d=0; d<up->grid.ndim; ++d)
    num_nodes_quad_interior[d] = gkyl_range_shape(&up->local, d)*num_quad_points;

  int num_nodes_quad_surf_in_dir[up->grid.ndim][GKYL_MAX_CDIM];
  for (int dir=0; dir<up->grid.ndim; ++dir)
    for (int d=0; d<up->grid.ndim; ++d)
      num_nodes_quad_surf_in_dir[dir][d] = d == dir ? gkyl_range_shape(&up->local, d)+1 : gkyl_range_shape(&up->local, d)*num_quad_points;

  gkyl_range_init_from_shape(&up->nrange_corn, up->grid.ndim, num_nodes_corners);
  gkyl_range_init_from_shape(&up->nrange_int, up->grid.ndim, num_nodes_quad_interior);
  for (int dir=0; dir<up->grid.ndim; ++dir)
    gkyl_range_init_from_shape(&up->nrange_surf[dir], up->grid.ndim, num_nodes_quad_surf_in_dir[dir]);
}

static void
gk_geometry_surf_alloc_nodal(struct gk_geometry* gk_geom, int dir)
{
  // mapc2p for calculations of tangents
  int num_fd_nodes = 13;
  gk_geom->geo_surf[dir].mc2p_nodal_fd = gkyl_array_new(GKYL_DOUBLE, 3*num_fd_nodes, gk_geom->nrange_surf[dir].volume);
  gk_geom->geo_surf[dir].mc2p_nodal = gkyl_array_new(GKYL_DOUBLE, 3, gk_geom->nrange_surf[dir].volume);
  // bmag,metrics and derived geo quantities
  gk_geom->geo_surf[dir].bmag_nodal = gkyl_array_new(GKYL_DOUBLE, 1, gk_geom->nrange_surf[dir].volume);
  gk_geom->geo_surf[dir].ddtheta_nodal = gkyl_array_new(GKYL_DOUBLE, 3, gk_geom->nrange_surf[dir].volume);
  gk_geom->geo_surf[dir].ddpsi_nodal = gkyl_array_new(GKYL_DOUBLE, 1, gk_geom->nrange_surf[dir].volume);
  gk_geom->geo_surf[dir].curlbhat_nodal = gkyl_array_new(GKYL_DOUBLE, 3, gk_geom->nrange_surf[dir].volume);
  gk_geom->geo_surf[dir].normcurlbhat_nodal = gkyl_array_new(GKYL_DOUBLE, 1, gk_geom->nrange_surf[dir].volume);
  gk_geom->geo_surf[dir].jacobgeo_nodal = gkyl_array_new(GKYL_DOUBLE, 1, gk_geom->nrange_surf[dir].volume);
  gk_geom->geo_surf[dir].b_i_nodal = gkyl_array_new(GKYL_DOUBLE, 3, gk_geom->nrange_surf[dir].volume);
  gk_geom->geo_surf[dir].b_i_nodal_fd = gkyl_array_new(GKYL_DOUBLE, 3*num_fd_nodes, gk_geom->nrange_surf[dir].volume);
  gk_geom->geo_surf[dir].cmag_nodal = gkyl_array_new(GKYL_DOUBLE, 1, gk_geom->nrange_surf[dir].volume);
  gk_geom->geo_surf[dir].jacobtot_inv_nodal = gkyl_array_new(GKYL_DOUBLE, 1, gk_geom->nrange_surf[dir].volume);
  gk_geom->geo_surf[dir].g_ij_nodal = gkyl_array_new(GKYL_DOUBLE, 6, gk_geom->nrange_surf[dir].volume);
  gk_geom->geo_surf[dir].dxdz_nodal = gkyl_array_new(GKYL_DOUBLE, 9, gk_geom->nrange_surf[dir].volume);
  gk_geom->geo_surf[dir].dzdx_nodal = gkyl_array_new(GKYL_DOUBLE, 9, gk_geom->nrange_surf[dir].volume);
  gk_geom->geo_surf[dir].normals_nodal = gkyl_array_new(GKYL_DOUBLE, 9, gk_geom->nrange_surf[dir].volume);
  gk_geom->geo_surf[dir].dualmag_nodal = gkyl_array_new(GKYL_DOUBLE, 3, gk_geom->nrange_surf[dir].volume);
  gk_geom->geo_surf[dir].bcart_nodal = gkyl_array_new(GKYL_DOUBLE, 3, gk_geom->nrange_surf[dir].volume);
  gk_geom->geo_surf[dir].B3_nodal = gkyl_array_new(GKYL_DOUBLE, 1, gk_geom->nrange_surf[dir].volume);
  gk_geom->geo_surf[dir].lenr_nodal = gkyl_array_new(GKYL_DOUBLE, 1, gk_geom->nrange_surf[dir].volume);
  gk_geom->geo_surf[dir].bimpactangle_nodal = gkyl_array_new(GKYL_DOUBLE, 1, gk_geom->nrange_surf[dir].volume);
  gk_geom->geo_surf[dir].deltats_nodal = gkyl_array_new(GKYL_DOUBLE, 1, gk_geom->nrange_surf[dir].volume);
}

static void
gk_geometry_surf_alloc_expansions(struct gk_geometry* up, int dir)
{
  up->geo_surf[dir].bmag = gkyl_array_new(GKYL_DOUBLE, 1*up->num_surf_basis, up->local_ext.volume);
  up->geo_surf[dir].jacobgeo = gkyl_array_new(GKYL_DOUBLE, 1*up->num_surf_basis, up->local_ext.volume);
  up->geo_surf[dir].jacobgeo_ratio = gkyl_array_new(GKYL_DOUBLE, 1*up->num_surf_basis, up->local_ext.volume);
  up->geo_surf[dir].b_i = gkyl_array_new(GKYL_DOUBLE, 3*up->num_surf_basis, up->local_ext.volume);
  up->geo_surf[dir].cmag = gkyl_array_new(GKYL_DOUBLE, 1*up->num_surf_basis, up->local_ext.volume);
  up->geo_surf[dir].jacobtot_inv = gkyl_array_new(GKYL_DOUBLE, 1*up->num_surf_basis, up->local_ext.volume);
  up->geo_surf[dir].B3 = gkyl_array_new(GKYL_DOUBLE, 1*up->num_surf_basis, up->local_ext.volume);
  up->geo_surf[dir].normcurlbhat = gkyl_array_new(GKYL_DOUBLE, 1*up->num_surf_basis, up->local_ext.volume);
  up->geo_surf[dir].normals = gkyl_array_new(GKYL_DOUBLE, 9*up->num_surf_basis, up->local_ext.volume);
  up->geo_surf[dir].lenr = gkyl_array_new(GKYL_DOUBLE, 1*up->num_surf_basis, up->local_ext.volume);
  up->geo_surf[dir].bimpactangle = gkyl_array_new(GKYL_DOUBLE, 1*up->num_surf_basis, up->local_ext.volume);
  up->geo_surf[dir].deltats = gkyl_array_new(GKYL_DOUBLE, 1*up->num_surf_basis, up->local_ext.volume);
}

static void
gk_geometry_surf_release_nodal(struct gk_geometry* gk_geom, int dir)
{
  gkyl_array_release(gk_geom->geo_surf[dir].mc2p_nodal_fd);
  gkyl_array_release(gk_geom->geo_surf[dir].mc2p_nodal);
  gkyl_array_release(gk_geom->geo_surf[dir].bmag_nodal);
  gkyl_array_release(gk_geom->geo_surf[dir].curlbhat_nodal);
  gkyl_array_release(gk_geom->geo_surf[dir].normcurlbhat_nodal);
  gkyl_array_release(gk_geom->geo_surf[dir].ddtheta_nodal);
  gkyl_array_release(gk_geom->geo_surf[dir].ddpsi_nodal);
  gkyl_array_release(gk_geom->geo_surf[dir].jacobgeo_nodal);
  gkyl_array_release(gk_geom->geo_surf[dir].b_i_nodal);
  gkyl_array_release(gk_geom->geo_surf[dir].b_i_nodal_fd);
  gkyl_array_release(gk_geom->geo_surf[dir].cmag_nodal);
  gkyl_array_release(gk_geom->geo_surf[dir].jacobtot_inv_nodal);
  gkyl_array_release(gk_geom->geo_surf[dir].g_ij_nodal);
  gkyl_array_release(gk_geom->geo_surf[dir].dxdz_nodal);
  gkyl_array_release(gk_geom->geo_surf[dir].dzdx_nodal);
  gkyl_array_release(gk_geom->geo_surf[dir].normals_nodal);
  gkyl_array_release(gk_geom->geo_surf[dir].dualmag_nodal);
  gkyl_array_release(gk_geom->geo_surf[dir].bcart_nodal);
  gkyl_array_release(gk_geom->geo_surf[dir].B3_nodal);
  gkyl_array_release(gk_geom->geo_surf[dir].lenr_nodal);
  gkyl_array_release(gk_geom->geo_surf[dir].bimpactangle_nodal);
  gkyl_array_release(gk_geom->geo_surf[dir].deltats_nodal);
}

static void
gk_geometry_int_alloc_nodal(struct gk_geometry* gk_geom)
{
  // mapc2p for calculations of tangents
  int num_fd_nodes = 13;
  gk_geom->geo_int.mc2p_nodal_fd = gkyl_array_new(GKYL_DOUBLE, 3*num_fd_nodes, gk_geom->nrange_int.volume);
  gk_geom->geo_int.mc2p_nodal = gkyl_array_new(GKYL_DOUBLE, 3, gk_geom->nrange_int.volume);
  // bmag.metrics and derived geo quantities
  gk_geom->geo_int.bmag_nodal = gkyl_array_new(GKYL_DOUBLE, 1, gk_geom->nrange_int.volume);
  gk_geom->geo_int.ddtheta_nodal = gkyl_array_new(GKYL_DOUBLE, 3, gk_geom->nrange_int.volume);
  gk_geom->geo_int.ddpsi_nodal = gkyl_array_new(GKYL_DOUBLE, 1, gk_geom->nrange_int.volume);
  gk_geom->geo_int.curlbhat_nodal = gkyl_array_new(GKYL_DOUBLE, 3, gk_geom->nrange_int.volume);
  gk_geom->geo_int.dualcurlbhat_nodal = gkyl_array_new(GKYL_DOUBLE, 3, gk_geom->nrange_int.volume);
  gk_geom->geo_int.jacobgeo_nodal = gkyl_array_new(GKYL_DOUBLE, 1, gk_geom->nrange_int.volume);
  gk_geom->geo_int.g_ij_nodal = gkyl_array_new(GKYL_DOUBLE, 6, gk_geom->nrange_int.volume);
  gk_geom->geo_int.g_ij_neut_nodal = gkyl_array_new(GKYL_DOUBLE, 6, gk_geom->nrange_int.volume);
  gk_geom->geo_int.dxdz_nodal = gkyl_array_new(GKYL_DOUBLE, 9, gk_geom->nrange_int.volume);
  gk_geom->geo_int.dzdx_nodal = gkyl_array_new(GKYL_DOUBLE, 9, gk_geom->nrange_int.volume);
  gk_geom->geo_int.dualmag_nodal = gkyl_array_new(GKYL_DOUBLE, 3, gk_geom->nrange_int.volume);
  gk_geom->geo_int.normals_nodal = gkyl_array_new(GKYL_DOUBLE, 9, gk_geom->nrange_int.volume);
  gk_geom->geo_int.gij_neut_nodal = gkyl_array_new(GKYL_DOUBLE, 6, gk_geom->nrange_int.volume);
  gk_geom->geo_int.b_i_nodal = gkyl_array_new(GKYL_DOUBLE, 3, gk_geom->nrange_int.volume);
  gk_geom->geo_int.b_i_nodal_fd = gkyl_array_new(GKYL_DOUBLE, 3*num_fd_nodes, gk_geom->nrange_int.volume);
  gk_geom->geo_int.bcart_nodal = gkyl_array_new(GKYL_DOUBLE, 3, gk_geom->nrange_int.volume);
  gk_geom->geo_int.B3_nodal = gkyl_array_new(GKYL_DOUBLE, 1, gk_geom->nrange_int.volume);
  gk_geom->geo_int.dualcurlbhatoverB_nodal = gkyl_array_new(GKYL_DOUBLE, 3, gk_geom->nrange_int.volume);
  gk_geom->geo_int.rtg33inv_nodal = gkyl_array_new(GKYL_DOUBLE, 1, gk_geom->nrange_int.volume);
  gk_geom->geo_int.bioverJB_nodal = gkyl_array_new(GKYL_DOUBLE, 3, gk_geom->nrange_int.volume);
}

static void
gk_geometry_int_alloc_expansions(struct gk_geometry* up)
{
  // mapc2p
  up->geo_int.mc2p = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
  // bmag, metrics and derived geo quantities
  up->geo_int.bmag = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->geo_int.g_ij = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->local_ext.volume);
  up->geo_int.g_ij_neut = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->local_ext.volume);
  up->geo_int.dxdz = gkyl_array_new(GKYL_DOUBLE, 9*up->basis.num_basis, up->local_ext.volume);
  up->geo_int.dzdx = gkyl_array_new(GKYL_DOUBLE, 9*up->basis.num_basis, up->local_ext.volume);
  up->geo_int.dualmag = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
  up->geo_int.normals = gkyl_array_new(GKYL_DOUBLE, 9*up->basis.num_basis, up->local_ext.volume);
  up->geo_int.jacobgeo = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->geo_int.jacobgeo_ghost = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->geo_int.jacobgeo_inv = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->geo_int.gij = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->local_ext.volume);
  up->geo_int.gij_neut = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->local_ext.volume);
  up->geo_int.b_i = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
  up->geo_int.bcart = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
  up->geo_int.cmag = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->geo_int.jacobtot = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->geo_int.jacobtot_inv = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->geo_int.gxxj= gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->geo_int.gxyj= gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->geo_int.gyyj= gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->geo_int.gxzj= gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->geo_int.eps2= gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->geo_int.dualcurlbhat = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
  up->geo_int.dualcurlbhatoverB = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
  up->geo_int.rtg33inv = gkyl_array_new(GKYL_DOUBLE, 1*up->basis.num_basis, up->local_ext.volume);
  up->geo_int.bioverJB = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
  up->geo_int.B3 = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->geo_int.qprofile = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
}

static void
gk_geometry_int_release_nodal(struct gk_geometry* gk_geom)
{
  gkyl_array_release(gk_geom->geo_int.mc2p_nodal_fd);
  gkyl_array_release(gk_geom->geo_int.mc2p_nodal);
  gkyl_array_release(gk_geom->geo_int.bmag_nodal);
  gkyl_array_release(gk_geom->geo_int.ddtheta_nodal);
  gkyl_array_release(gk_geom->geo_int.ddpsi_nodal);
  gkyl_array_release(gk_geom->geo_int.curlbhat_nodal);
  gkyl_array_release(gk_geom->geo_int.dualcurlbhat_nodal);
  gkyl_array_release(gk_geom->geo_int.jacobgeo_nodal);
  gkyl_array_release(gk_geom->geo_int.g_ij_nodal);
  gkyl_array_release(gk_geom->geo_int.g_ij_neut_nodal);
  gkyl_array_release(gk_geom->geo_int.dxdz_nodal);
  gkyl_array_release(gk_geom->geo_int.dzdx_nodal);
  gkyl_array_release(gk_geom->geo_int.dualmag_nodal);
  gkyl_array_release(gk_geom->geo_int.normals_nodal);
  gkyl_array_release(gk_geom->geo_int.gij_neut_nodal);
  gkyl_array_release(gk_geom->geo_int.b_i_nodal);
  gkyl_array_release(gk_geom->geo_int.b_i_nodal_fd);
  gkyl_array_release(gk_geom->geo_int.bcart_nodal);
  gkyl_array_release(gk_geom->geo_int.B3_nodal);
  gkyl_array_release(gk_geom->geo_int.dualcurlbhatoverB_nodal);
  gkyl_array_release(gk_geom->geo_int.rtg33inv_nodal);
  gkyl_array_release(gk_geom->geo_int.bioverJB_nodal);
}

static void
gk_geometry_corn_alloc_nodal(struct gk_geometry* gk_geom)
{
  // mapc2p
  gk_geom->geo_corn.mc2p_nodal = gkyl_array_new(GKYL_DOUBLE, 3, gk_geom->nrange_corn.volume);
  gk_geom->geo_corn.mc2nu_pos_nodal = gkyl_array_new(GKYL_DOUBLE, 3, gk_geom->nrange_corn.volume);
  gk_geom->geo_corn.bmag_nodal = gkyl_array_new(GKYL_DOUBLE, 1, gk_geom->nrange_corn.volume);
}

static void
gk_geometry_corn_alloc_expansions(struct gk_geometry* up)
{
  // mapc2p
  up->geo_corn.mc2p = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
  up->geo_corn.mc2nu_pos = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
  // bmag
  up->geo_corn.bmag = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->geo_corn.bmag_inv = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  // deflated quantities for plotting
  up->geo_corn.mc2p_deflated = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*up->basis.num_basis, up->local_ext.volume);
  up->geo_corn.mc2nu_pos_deflated = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*up->basis.num_basis, up->local_ext.volume);
}

static void
gk_geometry_corn_release_nodal(struct gk_geometry* gk_geom)
{
  gkyl_array_release(gk_geom->geo_corn.mc2p_nodal);
  gkyl_array_release(gk_geom->geo_corn.mc2nu_pos_nodal);
  gkyl_array_release(gk_geom->geo_corn.bmag_nodal);
}

static void
gk_geometry_surf_calc_expansions(struct gk_geometry* gk_geom, int dir, 
  struct gkyl_range nrange_quad_surf)
{
  struct gk_geom_surf up_surf = gk_geom->geo_surf[dir];
  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&gk_geom->basis, &gk_geom->grid, false);

  struct gkyl_range local_ext_in_dir;
  int lower[3] = {gk_geom->local.lower[0], gk_geom->local.lower[1], gk_geom->local.lower[2]};
  int upper[3] = {gk_geom->local.upper[0], gk_geom->local.upper[1], gk_geom->local.upper[2]};
  upper[dir] += 1;
  gkyl_sub_range_init(&local_ext_in_dir, &gk_geom->local_ext, lower, upper);

  gkyl_nodal_ops_n2m_surface(n2m, &gk_geom->surf_basis, &gk_geom->grid, &nrange_quad_surf, &local_ext_in_dir, 1, up_surf.bmag_nodal, up_surf.bmag, dir);
  gkyl_nodal_ops_n2m_surface(n2m, &gk_geom->surf_basis, &gk_geom->grid, &nrange_quad_surf, &local_ext_in_dir, 1, up_surf.jacobgeo_nodal, up_surf.jacobgeo, dir);
  gkyl_nodal_ops_n2m_surface(n2m, &gk_geom->surf_basis, &gk_geom->grid, &nrange_quad_surf, &local_ext_in_dir, 3, up_surf.b_i_nodal, up_surf.b_i, dir);
  gkyl_nodal_ops_n2m_surface(n2m, &gk_geom->surf_basis, &gk_geom->grid, &nrange_quad_surf, &local_ext_in_dir, 1, up_surf.cmag_nodal, up_surf.cmag, dir);
  gkyl_nodal_ops_n2m_surface(n2m, &gk_geom->surf_basis, &gk_geom->grid, &nrange_quad_surf, &local_ext_in_dir, 1, up_surf.jacobtot_inv_nodal, up_surf.jacobtot_inv, dir);
  gkyl_nodal_ops_n2m_surface(n2m, &gk_geom->surf_basis, &gk_geom->grid, &nrange_quad_surf, &local_ext_in_dir, 1, up_surf.B3_nodal, up_surf.B3, dir);
  gkyl_nodal_ops_n2m_surface(n2m, &gk_geom->surf_basis, &gk_geom->grid, &nrange_quad_surf, &local_ext_in_dir, 1, up_surf.normcurlbhat_nodal, up_surf.normcurlbhat, dir);
  gkyl_nodal_ops_n2m_surface(n2m, &gk_geom->surf_basis, &gk_geom->grid, &nrange_quad_surf, &local_ext_in_dir, 9, up_surf.normals_nodal, up_surf.normals, dir);
  gkyl_nodal_ops_n2m_surface(n2m, &gk_geom->surf_basis, &gk_geom->grid, &nrange_quad_surf, &local_ext_in_dir, 1, up_surf.lenr_nodal, up_surf.lenr, dir);
  gkyl_nodal_ops_n2m_surface(n2m, &gk_geom->surf_basis, &gk_geom->grid, &nrange_quad_surf, &local_ext_in_dir, 1, up_surf.bimpactangle_nodal, up_surf.bimpactangle, dir);
  gkyl_nodal_ops_n2m_surface(n2m, &gk_geom->surf_basis, &gk_geom->grid, &nrange_quad_surf, &local_ext_in_dir, 1, up_surf.deltats_nodal, up_surf.deltats, dir);

  // jacobgeo_ratio is not used in single block.
  int cdim = gk_geom->grid.ndim;
  gkyl_array_clear(up_surf.jacobgeo_ratio, 0.0);
  gkyl_array_shiftc(up_surf.jacobgeo_ratio, pow(sqrt(2.0),cdim), 0);

  gkyl_nodal_ops_release(n2m);
}


