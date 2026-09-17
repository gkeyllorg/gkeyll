#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_comm.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_priv.h>
#include <gkyl_gk_geometry_tok.h>
#include <gkyl_math.h>
#include <gkyl_nodal_ops.h>

#include <gkyl_tok_geo.h>
#include <gkyl_rz_calc_derived_geo.h>
#include <gkyl_calc_metric.h>
#include <gkyl_calc_bmag.h>
#include <assert.h>

struct gk_geometry*
gk_geometry_tok_init(struct gkyl_gk_geometry_inp *geometry_inp)
{

  struct gk_geometry *up = gkyl_calloc(1,sizeof(struct gk_geometry));
  up->geometry_id = geometry_inp->geometry_id;
  up->basis = geometry_inp->geo_basis;
  up->local = geometry_inp->geo_local;
  up->local_ext = geometry_inp->geo_local_ext;
  up->global = geometry_inp->geo_global;
  up->global_ext = geometry_inp->geo_global_ext;
  up->grid = geometry_inp->geo_grid;
  up->has_LCFS = geometry_inp->has_LCFS;
  if (up->has_LCFS) {
    up->x_LCFS = geometry_inp->x_LCFS;
    // Check that the split happens within the domain.
    assert((up->grid.lower[0] <= up->x_LCFS) && (up->x_LCFS <= up->grid.upper[0]));
    // Check that the split happens at a cell boundary;
    double needint = (up->x_LCFS - up->grid.lower[0])/up->grid.dx[0];
    double rem_floor = fabs(needint-floor(needint));
    double rem_ceil = fabs(needint-ceil(needint));
    if (rem_floor < 1.0e-12) {
      up->idx_LCFS_lo = (int) floor(needint);
    }
    else if (rem_ceil < 1.0e-12) {
      up->idx_LCFS_lo = (int) ceil(needint);
    }
    else {
      fprintf(stderr, "x_LCFS = %.9e must be at a cell boundary.\n", up->x_LCFS);
      assert(false);
    }
  }

  gk_geometry_set_nodal_ranges(up) ;

  // Initialize surface basis
  gkyl_cart_modal_serendip(&up->surf_basis, up->grid.ndim-1, up->basis.poly_order);
  up->num_surf_basis = up->surf_basis.num_basis;

  // Initialize tokamak geometry object from EFIT
  const struct gkyl_efit_inp inp = geometry_inp->efit_info;
  struct gkyl_tok_geo_grid_inp ginp = geometry_inp->tok_grid_info;
  ginp.cgrid = up->grid;
  ginp.cbasis = up->basis;
  struct gkyl_tok_geo *geo = gkyl_tok_geo_new(&inp, &ginp);
  up->geqdsk_sign_convention = geo->efit->sibry > geo->efit->simag ? 0 : 1;
  up->half_domain = ginp.half_domain ? 1 : 0;

  // Allocate nodal and modal arrays for corner, interior, and surface geo
  gk_geometry_corn_alloc_nodal(up);
  gk_geometry_corn_alloc_expansions(up);
  gk_geometry_int_alloc_nodal(up);
  gk_geometry_int_alloc_expansions(up);
  for (int dir=0; dir<up->grid.ndim; ++dir) {
    gk_geometry_surf_alloc_nodal(up, dir);
    gk_geometry_surf_alloc_expansions(up, dir);
  }

  gkyl_position_map_optimize(geometry_inp->position_map, up->grid, up->global);

  // Calculate bmag and mapc2p in cylindrical coords at corner nodes for
  // getting cell coordinates (used only for plotting).
  gkyl_tok_geo_calc(up, &up->nrange_corn, geo, &ginp, geometry_inp->position_map);
  // Calculate bmag and mapc2p in cylindrical coords at interior nodes for
  // calculating geo quantity volume expansions.
  gkyl_tok_geo_calc_interior(up, &up->nrange_int, up->dzc, geo, &ginp, geometry_inp->position_map);
  // Calculate bmag and mapc2p in cylindrical coords at surfaces.
  for (int dir = 0; dir <up->grid.ndim; dir++)
    gkyl_tok_geo_calc_surface(up, dir, &up->nrange_surf[dir], up->dzc, geo, &ginp, geometry_inp->position_map);

  // Now calculate the metrics at interior nodes.
  struct gkyl_calc_metric* mcalc = gkyl_calc_metric_new(&up->basis, &up->grid,
    &up->global, &up->global_ext, &up->local, &up->local_ext, true, false);
  gkyl_calc_metric_advance_rz_interior(mcalc, up);
  gkyl_array_copy(up->geo_int.jacobgeo_ghost, up->geo_int.jacobgeo);
  // Calculate neutral metrics at interior nodes.
  gkyl_calc_metric_advance_rz_neut_interior(mcalc, up);
  // Calculate the derived geometric quantities at interior nodes.
  gkyl_rz_calc_derived_geo *jcalculator = gkyl_rz_calc_derived_geo_new(&up->basis, &up->grid, 1, false);
  gkyl_rz_calc_derived_geo_advance(jcalculator, &up->local, up->geo_int.g_ij, up->geo_int.bmag, 
    up->geo_int.jacobgeo, up->geo_int.jacobgeo_inv, up->geo_int.gij, up->geo_int.b_i, up->geo_int.cmag, up->geo_int.jacobtot,
    up->geo_int.jacobtot_inv, up->geo_int.gxxj, up->geo_int.gxyj, up->geo_int.gyyj, up->geo_int.gxzj, up->geo_int.eps2);
  gkyl_rz_calc_derived_geo_release(jcalculator);
  // Calculate metrics/derived geo quantities at surface.
  for (int dir = 0; dir <up->grid.ndim; dir++) {
    gkyl_calc_metric_advance_rz_surface(mcalc, dir,  up);
  }
  gkyl_calc_metric_release(mcalc);
  // Calculate surface expansions.
  for (int dir = 0; dir <up->grid.ndim; dir++)
    gk_geometry_surf_calc_expansions(up, dir, up->nrange_surf[dir]);

  // Store metadata for I/O.
  char geqdsk_file_name[128];
  get_filename_from_path(geo->efit->filepath, geqdsk_file_name, sizeof(geqdsk_file_name));
  struct gkyl_msgpack_map_elem io_meta_basic[] = {
    { .key = "geometry_type", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = up->geometry_id },
    { .key = "geqdsk_sign_convention", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = up->geqdsk_sign_convention },
    { .key = "geqdsk_file", .elem_type = GKYL_MP_STRING, .cval = geqdsk_file_name},
    { .key = "half_domain", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = up->half_domain },
  };
  up->io_meta_basic_len = sizeof(io_meta_basic)/sizeof(io_meta_basic[0]);
  up->io_meta_basic = gkyl_msgpack_map_elem_clone(up->io_meta_basic_len, io_meta_basic);

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->ref_count = gkyl_ref_count_init(gkyl_gk_geometry_free);
  up->on_dev = up; // CPU eqn obj points to itself.

  gkyl_tok_geo_release(geo);

  return up;
}

struct gk_geometry*
gkyl_gk_geometry_tok_new(struct gkyl_gk_geometry_inp *geometry_inp)
{
  struct gk_geometry* gk_geom_3d;
  struct gk_geometry* gk_geom;

  if (geometry_inp->position_map->id == GKYL_PMAP_XPT_COMPRESSION) {
    double zcenter, zcut, len;
    switch(geometry_inp->tok_grid_info.ftype)
    {
      case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID:
      case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID:
        len = geometry_inp->tok_grid_info.half_domain ? 2.0*(geometry_inp->geo_grid.upper[2] - geometry_inp->geo_grid.lower[2])
                                                      : geometry_inp->geo_grid.upper[2] - geometry_inp->geo_grid.lower[2];
        zcut = len/2.0;
        zcenter = 0.0;
        break;
      case GKYL_GEOMETRY_TOKAMAK_CORE_R:
      case GKYL_GEOMETRY_TOKAMAK_CORE:
      case GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID:
        len = geometry_inp->tok_grid_info.half_domain ? 2.0*(geometry_inp->geo_grid.upper[2] - geometry_inp->geo_grid.lower[2])
                                                      : geometry_inp->geo_grid.upper[2] - geometry_inp->geo_grid.lower[2];
        zcenter = geometry_inp->geo_grid.lower[2] + len/2.0;
        zcut = len/2.0;
        break;
      case GKYL_GEOMETRY_TOKAMAK_CORE_L:
        len = geometry_inp->tok_grid_info.half_domain ? 2.0*(geometry_inp->geo_grid.upper[2] - geometry_inp->geo_grid.lower[2])
                                                      : geometry_inp->geo_grid.upper[2] - geometry_inp->geo_grid.lower[2];
        zcenter = geometry_inp->geo_grid.upper[2] - len/2.0;
        zcut = len/2.0;
        break;
      case GKYL_GEOMETRY_TOKAMAK_PF_LO_R:
      case GKYL_GEOMETRY_TOKAMAK_PF_UP_L:
      case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO:
      case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_UP:
      case GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO:
        len = geometry_inp->geo_grid.upper[2] - geometry_inp->geo_grid.lower[2];
        zcenter = geometry_inp->position_map->xpt_ctx->compress_divertor ? geometry_inp->geo_grid.lower[2] + len/2.0 : geometry_inp->geo_grid.lower[2];
        zcut = geometry_inp->position_map->xpt_ctx->compress_divertor ? len/2.0 : len;
        break;
      case GKYL_GEOMETRY_TOKAMAK_PF_LO_L:
      case GKYL_GEOMETRY_TOKAMAK_PF_UP_R:
      case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP:
      case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO:
      case GKYL_GEOMETRY_TOKAMAK_LSN_SOL_UP:
        len = geometry_inp->geo_grid.upper[2] - geometry_inp->geo_grid.lower[2];
        zcenter = geometry_inp->position_map->xpt_ctx->compress_divertor ? geometry_inp->geo_grid.upper[2] - len/2.0 : geometry_inp->geo_grid.upper[2];
        zcut = geometry_inp->position_map->xpt_ctx->compress_divertor ? len/2.0 : len;
        break;
      default:
        break;
    }

    double w = geometry_inp->geo_grid.upper[0] - geometry_inp->geo_grid.lower[0];
    struct gkyl_efit *efit = gkyl_efit_new(&geometry_inp->efit_info);
    double psisep = efit->psisep;
    gkyl_efit_release(efit);
    gkyl_position_map_set_compression(geometry_inp->position_map, zcut, zcenter, w, psisep);
  }
  else if (geometry_inp->position_map->id == GKYL_PMAP_CONSTANT_DB_POLYNOMIAL || \
           geometry_inp->position_map->id == GKYL_PMAP_CONSTANT_DB_NUMERIC) {
    // First construct the uniform 3d geometry
    gk_geom_3d = gk_geometry_tok_init(geometry_inp);
    // The array mc2nu is computed using the uniform geometry, so we need to deflate it
    // Must deflate the 3D uniform geometry in order for the allgather to work
    if(geometry_inp->grid.ndim < 3)
      gk_geom = gkyl_gk_geometry_deflate(gk_geom_3d, geometry_inp);
    else
      gk_geom = gkyl_gk_geometry_acquire(gk_geom_3d);

    gkyl_position_map_set_bmag(geometry_inp->position_map, geometry_inp->comm, \
      gk_geom->geo_int.bmag);

    gkyl_gk_geometry_release(gk_geom_3d); // release temporary 3d geometry
    gkyl_gk_geometry_release(gk_geom); // release 3d geometry
  }
  // Construct the non-uniform grid
  gk_geom_3d = gk_geometry_tok_init(geometry_inp);
  return gk_geom_3d;
}


void
gkyl_gk_geometry_tok_set_grid_extents(struct gkyl_efit_inp efit_info, struct gkyl_tok_geo_grid_inp grid_info, double *theta_lo, double *theta_up) {
  struct gkyl_tok_geo *geo = gkyl_tok_geo_new(&efit_info, &grid_info);
  gkyl_tok_geo_set_extent(&grid_info, geo, theta_lo, theta_up);
  gkyl_tok_geo_release(geo);
}
