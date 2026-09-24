#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <assert.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_priv.h>
#include <gkyl_gk_geometry_mirror.h>
#include <gkyl_position_map.h>
#include <gkyl_mirror_grid_gen.h>
#include <gkyl_rz_calc_derived_geo.h>
#include <gkyl_calc_metric_mirror.h>

struct gk_geometry*
gk_geometry_mirror_init(struct gkyl_gk_geometry_inp *geometry_inp)
{
  struct gk_geometry *up = gkyl_malloc(sizeof(struct gk_geometry));
  up->geometry_id = geometry_inp->geometry_id;
  up->basis = geometry_inp->geo_basis;
  up->local = geometry_inp->geo_local;
  up->local_ext = geometry_inp->geo_local_ext;
  up->global = geometry_inp->geo_global;
  up->global_ext = geometry_inp->geo_global_ext;
  up->grid = geometry_inp->geo_grid;
  up->geqdsk_sign_convention = 0.0; // Hardcoded 0. Means psi increases from axis. Always true for mirror geometry.
  up->half_domain = 0;
  up->has_LCFS = false; // MF 2025/09/19: Hardcodded for now.

  // Initialize nodal ranges for corner, interior, and surface geometry
  gk_geometry_set_nodal_ranges(up) ;
  // Initialize surface basis
  gkyl_cart_modal_serendip(&up->surf_basis, up->grid.ndim-1, up->basis.poly_order);
  up->num_surf_basis = up->surf_basis.num_basis;
  // Allocate nodal and modal arrays for corner, interior, and surface geo
  gk_geometry_corn_alloc_nodal(up);
  gk_geometry_corn_alloc_expansions(up);
  gk_geometry_int_alloc_nodal(up);
  gk_geometry_int_alloc_expansions(up);
  for (int dir=0; dir<up->grid.ndim; ++dir) {
    gk_geometry_surf_alloc_nodal(up, dir);
    gk_geometry_surf_alloc_expansions(up, dir);
  }

  // read psi(R,Z) from file
  struct gkyl_rect_grid psi_grid;
  struct gkyl_array *psi = gkyl_grid_array_new_from_file(&psi_grid, geometry_inp->mirror_grid_info.filename_psi);

  // Check if file exists using gkyl_check_file_exists and handle error only on rank 0
  if (!gkyl_check_file_exists(geometry_inp->mirror_grid_info.filename_psi)) {
    fprintf(stderr, "Failed to open the eqdsk file: %s\n", geometry_inp->mirror_grid_info.filename_psi);
    assert(false);
  }

  gkyl_position_map_optimize(geometry_inp->position_map, up->grid, up->global);

  // Create mirror geometry for corners.
  struct gkyl_mirror_grid_gen *mirror_grid_corn =
    gkyl_mirror_grid_gen_inew(&(struct gkyl_mirror_grid_gen_inp) {
        .comp_grid = &up->grid,
        .nrange = up->nrange_corn,
        .local = up->local,
        .global = up->global,
        .position_map = geometry_inp->position_map,
        
        .R = { psi_grid.lower[0], psi_grid.upper[0] },
        .Z = { psi_grid.lower[1], psi_grid.upper[1] },
        
        // psi(R,Z) grid size
        .nrcells = psi_grid.cells[0]-1, // Cells and not nodes.
        .nzcells = psi_grid.cells[1]-1, // Cells and not nodes.

        .psiRZ = psi,
        .fl_coord = geometry_inp->mirror_grid_info.fl_coord,
        .include_axis = geometry_inp->mirror_grid_info.include_axis,
        .write_psi_cubic = false,
      }
    );

  // Create mirror geometry for interior.
  struct gkyl_mirror_grid_gen *mirror_grid_int =
    gkyl_mirror_grid_gen_int_inew(&(struct gkyl_mirror_grid_gen_inp) {
        .comp_grid = &up->grid,
        .nrange = up->nrange_int,
        .local = up->local,
        .global = up->global,
        .position_map = geometry_inp->position_map,
        
        .R = { psi_grid.lower[0], psi_grid.upper[0] },
        .Z = { psi_grid.lower[1], psi_grid.upper[1] },
        
        // psi(R,Z) grid size.
        .nrcells = psi_grid.cells[0]-1, // Cells and not nodes.
        .nzcells = psi_grid.cells[1]-1, // Cells and not nodes.

        .psiRZ = psi,
        .fl_coord = geometry_inp->mirror_grid_info.fl_coord,
        .include_axis = geometry_inp->mirror_grid_info.include_axis,
        .write_psi_cubic = false,
      }
    );

  // create mirror geometry for surfaces
  struct gkyl_mirror_grid_gen *mirror_grid_surf[3];
  for (int dir = 0; dir <up->grid.ndim; dir++) {
    mirror_grid_surf[dir] =
      gkyl_mirror_grid_gen_surf_inew(&(struct gkyl_mirror_grid_gen_inp) {
          .comp_grid = &up->grid,
          .nrange = up->nrange_surf[dir],
          .local = up->local,
          .global = up->global,
          .dir = dir,
          .position_map = geometry_inp->position_map,
          
          .R = { psi_grid.lower[0], psi_grid.upper[0] },
          .Z = { psi_grid.lower[1], psi_grid.upper[1] },
          
          // psi(R,Z) grid size.
          .nrcells = psi_grid.cells[0]-1, // Cells and not nodes.
          .nzcells = psi_grid.cells[1]-1, // Cells and not nodes.

          .psiRZ = psi,
          .fl_coord = geometry_inp->mirror_grid_info.fl_coord,
          .include_axis = geometry_inp->mirror_grid_info.include_axis,
          .write_psi_cubic = false,
        }
      );
  }

  // Now calculate the derived geometric coefficients at necessary nodes and compute modal expansions where required.
  // Now calculate the metrics at corner and interior nodes.
  struct gkyl_calc_metric_mirror* mcalc = gkyl_calc_metric_mirror_new(&up->basis, &up->grid, &up->local, &up->local_ext, false);
  gkyl_calc_metric_mirror_advance(mcalc, up, mirror_grid_corn);
  gkyl_calc_metric_mirror_advance_interior(mcalc, up, mirror_grid_int);
  // Calculate the derived geometric quantities at interior nodes.
  gkyl_rz_calc_derived_geo *jcalculator = gkyl_rz_calc_derived_geo_new(&up->basis, &up->grid, 1, false);
  gkyl_rz_calc_derived_geo_advance(jcalculator, &up->local, up->geo_int.g_ij, up->geo_int.bmag, 
    up->geo_int.jacobgeo, up->geo_int.jacobgeo_inv, up->geo_int.gij, up->geo_int.b_i, up->geo_int.cmag, up->geo_int.jacobtot,
    up->geo_int.jacobtot_inv, up->geo_int.gxxj, up->geo_int.gxyj, up->geo_int.gyyj, up->geo_int.gxzj, up->geo_int.eps2);
  gkyl_rz_calc_derived_geo_release(jcalculator);
  // Calculate metrics/derived geo quantities at surface.
  for (int dir = 0; dir <up->grid.ndim; dir++) {
    gkyl_calc_metric_mirror_advance_surface(mcalc, dir,  up, mirror_grid_surf[dir]);
  }
  gkyl_calc_metric_mirror_release(mcalc);
  // Calculate surface expansions.
  for (int dir = 0; dir <up->grid.ndim; dir++)
    gk_geometry_surf_calc_expansions(up, dir, up->nrange_surf[dir]);

  // Store metadata for I/O.
  struct gkyl_msgpack_map_elem io_meta_basic[] = {
    { .key = "geometry_type", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = up->geometry_id },
    { .key = "geqdsk_sign_convention", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = up->geqdsk_sign_convention },
  };
  up->io_meta_basic_len = sizeof(io_meta_basic)/sizeof(io_meta_basic[0]);
  up->io_meta_basic = gkyl_msgpack_map_elem_clone(up->io_meta_basic_len, io_meta_basic);

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->ref_count = gkyl_ref_count_init(gkyl_gk_geometry_free);
  up->on_dev = up; // CPU eqn obj points to itself.

  gkyl_mirror_grid_gen_release(mirror_grid_corn);
  gkyl_mirror_grid_gen_release(mirror_grid_int);
  for (int dir = 0; dir <up->grid.ndim; dir++)
    gkyl_mirror_grid_gen_release(mirror_grid_surf[dir]);
  gkyl_array_release(psi);

  return up;
}

struct gk_geometry*
gkyl_gk_geometry_mirror_new(struct gkyl_gk_geometry_inp *geometry_inp)
{
  struct gk_geometry* gk_geom_3d;
  struct gk_geometry* gk_geom;

  // First construct the uniform 3d geometry
  gk_geom_3d = gk_geometry_mirror_init(geometry_inp);
  if (geometry_inp->position_map->id == GKYL_PMAP_CONSTANT_DB_POLYNOMIAL || \
      geometry_inp->position_map->id == GKYL_PMAP_CONSTANT_DB_NUMERIC) {
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

    // Construct the non-uniform grid
    gk_geom_3d = gk_geometry_mirror_init(geometry_inp);
  }
  return gk_geom_3d;
}
