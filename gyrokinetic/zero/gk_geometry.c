#include <gkyl_range.h>
#include <gkyl_alloc.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_math.h>
#include <gkyl_basis.h>
#include <gkyl_deflate_geo.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_priv.h>
#include <gkyl_nodal_ops.h>
#include <gkyl_alloc_flags_priv.h>
#include <assert.h>
#include <float.h>
#include <ctype.h>


struct gk_geometry*
gkyl_gk_geometry_new(struct gk_geometry* geo_host, struct gkyl_gk_geometry_inp *geometry_inp, bool use_gpu)
{

#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    return gkyl_gk_geometry_cu_dev_new(geo_host, geometry_inp);
  } 
#endif 

  struct gk_geometry *up = gkyl_malloc(sizeof(struct gk_geometry));
  up->geometry_id = geometry_inp->geometry_id;
  up->basis = geometry_inp->basis;
  up->local = geometry_inp->local;
  up->local_ext = geometry_inp->local_ext;
  up->global = geometry_inp->global;
  up->global_ext = geometry_inp->global_ext;
  up->grid = geometry_inp->grid;
  gk_geometry_set_nodal_ranges(up) ;

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

  // Function pointers to twistshift function.
  up->parallel_lower_bc_shift_func = geometry_inp->parallel_lower_bc_shift_func;
  up->parallel_upper_bc_shift_func = geometry_inp->parallel_upper_bc_shift_func;
  up->parallel_lower_bc_shift_ctx  = geometry_inp->parallel_lower_bc_shift_ctx;
  up->parallel_upper_bc_shift_ctx  = geometry_inp->parallel_upper_bc_shift_ctx;

  if (up->grid.ndim > 1) {
    gkyl_cart_modal_serendip(&up->surf_basis, up->grid.ndim-1, up->basis.poly_order);
    up->num_surf_basis = up->surf_basis.num_basis;
  }
  else {
    up->num_surf_basis = 1;
  }

  gk_geometry_corn_alloc_expansions(up);
  gk_geometry_corn_alloc_nodal(up);
  gk_geometry_int_alloc_expansions(up);
  gk_geometry_int_alloc_nodal(up);
  for (int dir=0; dir<up->grid.ndim; ++dir) {
    gk_geometry_surf_alloc_expansions(up, dir);
    gk_geometry_surf_alloc_nodal(up, dir);
  }

  // Store metadata for I/O.
  if (up->geometry_id == GKYL_GEOMETRY_TOKAMAK) {
    struct gkyl_msgpack_map_elem io_meta_basic[] = {
      { .key = "geometry_type", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = up->geometry_id },
      { .key = "geqdsk_sign_convention", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = up->geqdsk_sign_convention },
      { .key = "half_domain", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = up->half_domain },
    };
    up->io_meta_basic_len = sizeof(io_meta_basic)/sizeof(io_meta_basic[0]);
    up->io_meta_basic = gkyl_msgpack_map_elem_clone(up->io_meta_basic_len, io_meta_basic);
  }
  else if (up->geometry_id == GKYL_GEOMETRY_MIRROR) {
    struct gkyl_msgpack_map_elem io_meta_basic[] = {
      { .key = "geometry_type", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = up->geometry_id },
      { .key = "geqdsk_sign_convention", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = up->geqdsk_sign_convention },
    };
    up->io_meta_basic_len = sizeof(io_meta_basic)/sizeof(io_meta_basic[0]);
    up->io_meta_basic = gkyl_msgpack_map_elem_clone(up->io_meta_basic_len, io_meta_basic);
  }
  else {
    struct gkyl_msgpack_map_elem io_meta_basic[] = {
      { .key = "geometry_type", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = up->geometry_id },
    };
    up->io_meta_basic_len = sizeof(io_meta_basic)/sizeof(io_meta_basic[0]);
    up->io_meta_basic = gkyl_msgpack_map_elem_clone(up->io_meta_basic_len, io_meta_basic);
  }

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->ref_count = gkyl_ref_count_init(gkyl_gk_geometry_free);
  up->on_dev = up; // CPU eqn obj points to itself
                   
  return up;
}

void gkyl_gk_geometry_reset_io_meta(struct gk_geometry *up)
{
  gkyl_msgpack_map_elem_set_uint(up->io_meta_basic_len, up->io_meta_basic, "geometry_type", up->geometry_id);

  if (up->geometry_id == GKYL_GEOMETRY_TOKAMAK || up->geometry_id == GKYL_GEOMETRY_MIRROR) {
    if (gkyl_msgpack_map_elem_has_key(up->io_meta_basic_len, up->io_meta_basic, "geqdsk_sign_convention")) {
      // Element list has this key. Update its value.
      gkyl_msgpack_map_elem_set_uint(up->io_meta_basic_len, up->io_meta_basic, "geqdsk_sign_convention", up->geqdsk_sign_convention);
    }
    else {
      // Element list doesn't have this key. Create a new list with it.
      struct gkyl_msgpack_map_elem io_meta_basic_new[] = {
        { .key = "geqdsk_sign_convention", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = up->geqdsk_sign_convention },
      };
      int io_meta_basic_new_len = sizeof(io_meta_basic_new)/sizeof(io_meta_basic_new[0]);

      struct gkyl_msgpack_map_elem *io_meta_basic_buffer = gkyl_msgpack_map_elem_clone(up->io_meta_basic_len, up->io_meta_basic);
      int io_meta_basic_buffer_len = up->io_meta_basic_len;
      gkyl_msgpack_map_elem_release(io_meta_basic_buffer_len, up->io_meta_basic); 

      int io_meta_basic_list_len[] = {io_meta_basic_new_len, io_meta_basic_buffer_len};
      const struct gkyl_msgpack_map_elem* io_meta_basic_list[] = {io_meta_basic_new, io_meta_basic_buffer};
      up->io_meta_basic = gkyl_msgpack_map_elem_union(sizeof(io_meta_basic_list_len)/sizeof(int),
        io_meta_basic_list_len, io_meta_basic_list, &up->io_meta_basic_len);

      gkyl_msgpack_map_elem_release(io_meta_basic_buffer_len, io_meta_basic_buffer); 
    }
  }

  if (up->geometry_id == GKYL_GEOMETRY_TOKAMAK) {
    if (gkyl_msgpack_map_elem_has_key(up->io_meta_basic_len, up->io_meta_basic, "half_domain")) {
      // Element list has this key. Update its value.
      gkyl_msgpack_map_elem_set_uint(up->io_meta_basic_len, up->io_meta_basic, "half_domain", up->half_domain);
    }
    else {
      // Element list doesn't have this key. Create a new list with it.
      struct gkyl_msgpack_map_elem io_meta_basic_new[] = {
        { .key = "half_domain", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = up->half_domain },
      };
      int io_meta_basic_new_len = sizeof(io_meta_basic_new)/sizeof(io_meta_basic_new[0]);

      struct gkyl_msgpack_map_elem *io_meta_basic_buffer = gkyl_msgpack_map_elem_clone(up->io_meta_basic_len, up->io_meta_basic);
      int io_meta_basic_buffer_len = up->io_meta_basic_len;
      gkyl_msgpack_map_elem_release(io_meta_basic_buffer_len, up->io_meta_basic); 

      int io_meta_basic_list_len[] = {io_meta_basic_new_len, io_meta_basic_buffer_len};
      const struct gkyl_msgpack_map_elem* io_meta_basic_list[] = {io_meta_basic_new, io_meta_basic_buffer};
      up->io_meta_basic = gkyl_msgpack_map_elem_union(sizeof(io_meta_basic_list_len)/sizeof(int),
        io_meta_basic_list_len, io_meta_basic_list, &up->io_meta_basic_len);

      gkyl_msgpack_map_elem_release(io_meta_basic_buffer_len, io_meta_basic_buffer); 
    }
  }
}

void gkyl_gk_geometry_populate_nodal(struct gk_geometry *gk_geom)
{
  // Populate nodal arrays for vpar surface terms
  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&gk_geom->basis, &gk_geom->grid, false);
  gkyl_nodal_ops_m2n(n2m, &gk_geom->basis, &gk_geom->grid, &gk_geom->nrange_int, &gk_geom->local, 1, gk_geom->geo_int.jacobgeo_nodal, gk_geom->geo_int.jacobgeo, true);
  gkyl_nodal_ops_m2n(n2m, &gk_geom->basis, &gk_geom->grid, &gk_geom->nrange_int, &gk_geom->local, 9, gk_geom->geo_int.dxdz_nodal, gk_geom->geo_int.dxdz, true);
  gkyl_nodal_ops_m2n(n2m, &gk_geom->basis, &gk_geom->grid, &gk_geom->nrange_int, &gk_geom->local, 9, gk_geom->geo_int.dzdx_nodal, gk_geom->geo_int.dzdx, true);
  gkyl_nodal_ops_m2n(n2m, &gk_geom->basis, &gk_geom->grid, &gk_geom->nrange_int, &gk_geom->local, 1, gk_geom->geo_int.bmag_nodal, gk_geom->geo_int.bmag, true);
  gkyl_nodal_ops_m2n(n2m, &gk_geom->basis, &gk_geom->grid, &gk_geom->nrange_int, &gk_geom->local, 1, gk_geom->geo_int.B3_nodal, gk_geom->geo_int.B3, true);
  gkyl_nodal_ops_m2n(n2m, &gk_geom->basis, &gk_geom->grid, &gk_geom->nrange_int, &gk_geom->local, 3, gk_geom->geo_int.dualcurlbhat_nodal, gk_geom->geo_int.dualcurlbhat, true);

  // Populate nodal arrays for config space surface terms
  int lower[gk_geom->grid.ndim];
  int upper[gk_geom->grid.ndim];
  for (int dir =0; dir < gk_geom->grid.ndim; dir++) {
    struct gk_geom_surf up_surf = gk_geom->geo_surf[dir];

    struct gkyl_range local_ext_in_dir;
    for (int i = 0; i <gk_geom->grid.ndim; i++) {
      lower[i] = gk_geom->local.lower[i];
      upper[i] = gk_geom->local.upper[i];
    }
    upper[dir]+=1;
    gkyl_sub_range_init(&local_ext_in_dir, &gk_geom->local_ext, lower, upper);

    gkyl_nodal_ops_m2n_surface(n2m, &gk_geom->surf_basis, &gk_geom->grid, &gk_geom->nrange_surf[dir], &local_ext_in_dir, 1, up_surf.bmag_nodal, up_surf.bmag, dir);
    gkyl_nodal_ops_m2n_surface(n2m, &gk_geom->surf_basis, &gk_geom->grid, &gk_geom->nrange_surf[dir], &local_ext_in_dir, 1, up_surf.jacobgeo_nodal, up_surf.jacobgeo, dir);
    gkyl_nodal_ops_m2n_surface(n2m, &gk_geom->surf_basis, &gk_geom->grid, &gk_geom->nrange_surf[dir], &local_ext_in_dir, 3, up_surf.b_i_nodal, up_surf.b_i, dir);
    gkyl_nodal_ops_m2n_surface(n2m, &gk_geom->surf_basis, &gk_geom->grid, &gk_geom->nrange_surf[dir], &local_ext_in_dir, 1, up_surf.cmag_nodal, up_surf.cmag, dir);
    gkyl_nodal_ops_m2n_surface(n2m, &gk_geom->surf_basis, &gk_geom->grid, &gk_geom->nrange_surf[dir], &local_ext_in_dir, 1, up_surf.jacobtot_inv_nodal, up_surf.jacobtot_inv, dir);
    gkyl_nodal_ops_m2n_surface(n2m, &gk_geom->surf_basis, &gk_geom->grid, &gk_geom->nrange_surf[dir], &local_ext_in_dir, 1, up_surf.B3_nodal, up_surf.B3, dir);
    gkyl_nodal_ops_m2n_surface(n2m, &gk_geom->surf_basis, &gk_geom->grid, &gk_geom->nrange_surf[dir], &local_ext_in_dir, 1, up_surf.normcurlbhat_nodal, up_surf.normcurlbhat, dir);
    gkyl_nodal_ops_m2n_surface(n2m, &gk_geom->surf_basis, &gk_geom->grid, &gk_geom->nrange_surf[dir], &local_ext_in_dir, 9, up_surf.normals_nodal, up_surf.normals, dir);
    gkyl_nodal_ops_m2n_surface(n2m, &gk_geom->surf_basis, &gk_geom->grid, &gk_geom->nrange_surf[dir], &local_ext_in_dir, 1, up_surf.lenr_nodal, up_surf.lenr, dir);
    gkyl_nodal_ops_m2n_surface(n2m, &gk_geom->surf_basis, &gk_geom->grid, &gk_geom->nrange_surf[dir], &local_ext_in_dir, 1, up_surf.deltats_nodal, up_surf.deltats, dir);
    gkyl_nodal_ops_m2n_surface(n2m, &gk_geom->surf_basis, &gk_geom->grid, &gk_geom->nrange_surf[dir], &local_ext_in_dir, 1, up_surf.bimpactangle_nodal, up_surf.bimpactangle, dir);
  }
  gkyl_nodal_ops_release(n2m);

}

bool
gkyl_gk_geometry_is_cu_dev(const struct gk_geometry* up)
{
  return GKYL_IS_CU_ALLOC(up->flags);
}

struct gkyl_rect_grid gkyl_gk_geometry_augment_grid(struct gkyl_rect_grid grid, struct gkyl_gk_geometry_inp geometry)
{
  struct gkyl_rect_grid augmented_grid;
  int cells[3];
  double lower[3];
  double upper[3];

  if (grid.ndim==1) {
    cells[0] = 1;
    cells[1] = 1;
    cells[2] = grid.cells[0];

    lower[0] = geometry.world[0] - (geometry.world[0]>1e-14? fmin(1e-5, geometry.world[0]*0.1) : 1e-5);
    lower[1] = geometry.world[1] - 1e-1;
    lower[2] = grid.lower[0];

    upper[0] = geometry.world[0] + (geometry.world[0]>1e-14? fmin(1e-5, geometry.world[0]*0.1) : 1e-5);
    upper[1] = geometry.world[1] + 1e-1;
    upper[2] = grid.upper[0];
  }
  else if (grid.ndim==2) {
    cells[0] = grid.cells[0];
    cells[1] = 1;
    cells[2] = grid.cells[1];

    lower[0] = grid.lower[0];
    lower[1] = geometry.world[0] - 1e-5;
    lower[2] = grid.lower[1];

    upper[0] = grid.upper[0];
    upper[1] = geometry.world[0] + 1e-5;
    upper[2] = grid.upper[1];
  }

  gkyl_rect_grid_init(&augmented_grid, 3, lower, upper, cells);
  return augmented_grid;
}


void gkyl_gk_geometry_augment_local(const struct gkyl_range *inrange,
  const int *nghost, struct gkyl_range *ext_range, struct gkyl_range *range)
{
  if (inrange->ndim == 2) {
    int lower_ext[GKYL_MAX_DIM], upper_ext[GKYL_MAX_DIM];
    int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
    
    lower_ext[0] = inrange->lower[0]-nghost[0];
    upper_ext[0] = inrange->upper[0]+nghost[0];
    lower[0] = inrange->lower[0];
    upper[0] = inrange->upper[0];

    lower_ext[1] = 1 - 1;
    upper_ext[1] = 1 + 1;
    lower[1] = 1;
    upper[1] = 1;

    lower_ext[2] = inrange->lower[1]-nghost[1];
    upper_ext[2] = inrange->upper[1]+nghost[1];
    lower[2] = inrange->lower[1];
    upper[2] = inrange->upper[1];


    gkyl_range_init(ext_range, inrange->ndim+1, lower_ext, upper_ext);
    gkyl_sub_range_init(range, ext_range, lower, upper);  
  }
  else if (inrange->ndim == 1) {
    int lower_ext[GKYL_MAX_DIM], upper_ext[GKYL_MAX_DIM];
    int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
    
    lower_ext[0] = 1 - 1;
    upper_ext[0] = 1 + 1;
    lower[0] = 1;
    upper[0] = 1;

    lower_ext[1] = 1 - 1;
    upper_ext[1] = 1 + 1;
    lower[1] = 1;
    upper[1] = 1;

    lower_ext[2] = inrange->lower[0]-nghost[0];
    upper_ext[2] = inrange->upper[0]+nghost[0];
    lower[2] = inrange->lower[0];
    upper[2] = inrange->upper[0];



    gkyl_range_init(ext_range, inrange->ndim+2, lower_ext, upper_ext);
    gkyl_sub_range_init(range, ext_range, lower, upper);  
  }
}

double
gkyl_gk_geometry_reduce_bmag(struct gk_geometry* up, enum gkyl_array_op op)
{
  int cdim = up->grid.ndim;
  double b_m;
  if (op == GKYL_MIN)
    b_m = DBL_MAX;
  else if (op == GKYL_MAX)
    b_m = -DBL_MAX;
  else
    assert(false);

  struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, cdim, up->basis.num_basis);
  up->basis.node_list(gkyl_array_fetch(nodes, 0));

  struct gkyl_array *bmag_ho = gkyl_array_new(GKYL_DOUBLE, up->geo_int.bmag->ncomp, up->geo_int.bmag->size);
  gkyl_array_copy(bmag_ho, up->geo_int.bmag);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&up->local, iter.idx);
    double *b_d = gkyl_array_fetch(bmag_ho, linidx);
    double nod_log[cdim];
    for (int n = 0; n < up->basis.num_basis; n++) {
      const double *nod_log = gkyl_array_cfetch(nodes,n);
      double b = up->basis.eval_expand(nod_log, b_d);
      if (op == GKYL_MIN)
        b_m = GKYL_MIN2(b_m, b);
      else if (op == GKYL_MAX)
        b_m = GKYL_MAX2(b_m, b);
    }
  }

  gkyl_array_release(nodes);
  gkyl_array_release(bmag_ho);

  return b_m;
}

static inline void
log_to_comp(int ndim, const double *eta,
  const double * GKYL_RESTRICT dx, const double * GKYL_RESTRICT xc,
  double* GKYL_RESTRICT xout)
{
  // Convert logical to computational coordinates.
  for (int d=0; d<ndim; ++d) xout[d] = 0.5*dx[d]*eta[d]+xc[d];
}

double
gkyl_gk_geometry_reduce_arg_bmag(struct gk_geometry* up, enum gkyl_array_op op, double *coord)
{
  int cdim = up->grid.ndim;
  for (int d=0; d<cdim; d++)
    coord[d] = DBL_MAX; 

  double b_m;
  if (op == GKYL_MIN)
    b_m = DBL_MAX;
  else if (op == GKYL_MAX)
    b_m = -DBL_MAX;
  else
    assert(false);

  struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, cdim, up->basis.num_basis);
  up->basis.node_list(gkyl_array_fetch(nodes, 0));

  struct gkyl_array *bmag_ho = gkyl_array_new(GKYL_DOUBLE, up->geo_int.bmag->ncomp, up->geo_int.bmag->size);
  gkyl_array_copy(bmag_ho, up->geo_int.bmag);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&up->local, iter.idx);
    double *b_d = gkyl_array_fetch(bmag_ho, linidx);

    double xc[cdim];
    gkyl_rect_grid_cell_center(&up->grid, iter.idx, xc);

    double nod_log[cdim], nod_phys[cdim];
    for (int n = 0; n < up->basis.num_basis; n++) {
      const double *nod_log = gkyl_array_cfetch(nodes,n);
      double b = up->basis.eval_expand(nod_log, b_d);

      log_to_comp(cdim, nod_log, up->grid.dx, xc, nod_phys);

      if (op == GKYL_MIN) {
        if (b < b_m) {
          b_m = b;
          for (int d=0; d<cdim; d++)
            coord[d] = nod_phys[d];
        }
      }
      else if (op == GKYL_MAX) {
        if (b_m < b) {
          b_m = b;
          for (int d=0; d<cdim; d++)
            coord[d] = nod_phys[d];
        }
      }
    }
  }

  gkyl_array_release(nodes);
  gkyl_array_release(bmag_ho);

  return b_m;
}

void
gkyl_gk_geometry_init_nodal_range( struct gkyl_range *nrange, struct gkyl_range *range, int poly_order)
{
    int nodes[GKYL_MAX_DIM];
    if (poly_order == 1) {
      for (int d=0; d<range->ndim; ++d)
        nodes[d] = gkyl_range_shape(range, d) + 1;
    }
    if (poly_order == 2) {
      for (int d=0; d<range->ndim; ++d)
        nodes[d] = 2*gkyl_range_shape(range, d) + 1;
    }
    gkyl_range_init_from_shape(nrange, range->ndim, nodes);

}

void
gkyl_gk_geometry_init_nodal_grid(struct gkyl_rect_grid *ngrid, struct gkyl_rect_grid *grid, struct gkyl_range *nrange)
{
    double lower[GKYL_MAX_DIM];
    double upper[GKYL_MAX_DIM];
    int cells[GKYL_MAX_DIM];
    for (int i=0; i<nrange->ndim; ++i) {
      lower[i] = grid->lower[i];
      upper[i] = grid->upper[i];
      cells[i] = gkyl_range_shape(nrange, i);
    }
    gkyl_rect_grid_init(ngrid, nrange->ndim, lower, upper, cells);
}

struct gk_geometry*
gkyl_gk_geometry_deflate(const struct gk_geometry* up_3d, struct gkyl_gk_geometry_inp *geometry_inp)
{
  struct gk_geometry *up = gkyl_malloc(sizeof(struct gk_geometry));
  up->geometry_id = geometry_inp->geometry_id;
  up->basis = geometry_inp->basis;
  up->local = geometry_inp->local;
  up->local_ext = geometry_inp->local_ext;
  up->grid = geometry_inp->grid;
  gk_geometry_set_nodal_ranges(up) ;
  if (up->grid.ndim > 1) {
    gkyl_cart_modal_serendip(&up->surf_basis, up->grid.ndim-1, up->basis.poly_order);
    up->num_surf_basis = up->surf_basis.num_basis;
  }
  else {
    up->num_surf_basis = 1;
  }
  up->geqdsk_sign_convention = up_3d->geqdsk_sign_convention;
  up->half_domain = up_3d->half_domain;
  up->has_LCFS = up_3d->has_LCFS;
  up->x_LCFS = up_3d->x_LCFS;
  up->idx_LCFS_lo = up_3d->idx_LCFS_lo;

  gk_geometry_corn_alloc_expansions(up);
  gk_geometry_corn_alloc_nodal(up);
  gk_geometry_int_alloc_expansions(up);
  gk_geometry_int_alloc_nodal(up);
  for (int dir=0; dir<up->grid.ndim; ++dir) {
    gk_geometry_surf_alloc_expansions(up, dir);
    gk_geometry_surf_alloc_nodal(up, dir);
  }

  // Now fill the arrays by deflation
  int rem_dirs[3] = {0};
  if (up->grid.ndim==1) {
    rem_dirs[0] = 1;
    rem_dirs[1] = 1;
  }
  else if (up->grid.ndim==2) {
    rem_dirs[1] = 1;
  }
  struct gkyl_deflate_geo* deflator = gkyl_deflate_geo_new(&up_3d->basis, &up->basis, &up_3d->grid, &up->grid, rem_dirs, false);

  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->geo_corn.mc2p, up->geo_corn.mc2p, 3);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->geo_corn.mc2nu_pos, up->geo_corn.mc2nu_pos, 3);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->geo_corn.bmag, up->geo_corn.bmag, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->geo_corn.bmag_inv, up->geo_corn.bmag_inv, 1);

  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->geo_int.mc2p, up->geo_int.mc2p, 3);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->geo_int.bmag, up->geo_int.bmag, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->geo_int.g_ij, up->geo_int.g_ij, 6);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->geo_int.g_ij_neut, up->geo_int.g_ij_neut, 6);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->geo_int.dxdz, up->geo_int.dxdz, 9);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->geo_int.dzdx, up->geo_int.dzdx, 9);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->geo_int.dualmag, up->geo_int.dualmag, 3);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->geo_int.normals, up->geo_int.normals, 9);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->geo_int.jacobgeo, up->geo_int.jacobgeo, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->geo_int.jacobgeo_ghost, up->geo_int.jacobgeo_ghost, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->geo_int.jacobgeo_inv, up->geo_int.jacobgeo_inv, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->geo_int.gij, up->geo_int.gij, 6);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->geo_int.gij_neut, up->geo_int.gij_neut, 6);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->geo_int.b_i, up->geo_int.b_i, 3);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->geo_int.bcart, up->geo_int.bcart, 3);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->geo_int.cmag, up->geo_int.cmag, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->geo_int.jacobtot, up->geo_int.jacobtot, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->geo_int.jacobtot_inv, up->geo_int.jacobtot_inv, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->geo_int.gxxj, up->geo_int.gxxj, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->geo_int.gxyj, up->geo_int.gxyj, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->geo_int.gyyj, up->geo_int.gyyj, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->geo_int.gxzj, up->geo_int.gxzj, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->geo_int.eps2, up->geo_int.eps2, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->geo_int.dualcurlbhat, up->geo_int.dualcurlbhat, 3);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->geo_int.dualcurlbhatoverB, up->geo_int.dualcurlbhatoverB, 3);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->geo_int.rtg33inv, up->geo_int.rtg33inv, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->geo_int.bioverJB, up->geo_int.bioverJB, 3);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->geo_int.B3, up->geo_int.B3, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->geo_int.qprofile, up->geo_int.qprofile, 1);
  // Done deflating modal

  // Deflate nodal quantities
  gkyl_deflate_geo_advance_nodal(deflator, &up_3d->nrange_int, &up->nrange_int, up_3d->geo_int.jacobgeo_nodal, up->geo_int.jacobgeo_nodal, 1);
  gkyl_deflate_geo_advance_nodal(deflator, &up_3d->nrange_int, &up->nrange_int, up_3d->geo_int.dxdz_nodal, up->geo_int.dxdz_nodal, 9);
  gkyl_deflate_geo_advance_nodal(deflator, &up_3d->nrange_int, &up->nrange_int, up_3d->geo_int.dzdx_nodal, up->geo_int.dzdx_nodal, 9);
  gkyl_deflate_geo_advance_nodal(deflator, &up_3d->nrange_int, &up->nrange_int, up_3d->geo_int.B3_nodal, up->geo_int.B3_nodal, 1);
  gkyl_deflate_geo_advance_nodal(deflator, &up_3d->nrange_int, &up->nrange_int, up_3d->geo_int.bmag_nodal, up->geo_int.bmag_nodal, 1);
  gkyl_deflate_geo_advance_nodal(deflator, &up_3d->nrange_int, &up->nrange_int, up_3d->geo_int.dualcurlbhat_nodal, up->geo_int.dualcurlbhat_nodal, 3);
  gkyl_deflate_geo_release(deflator);

  if (up->grid.ndim==1) {
    // In 1D geometry, make mapc2p a function of only Z and mc2nu_pos only a function of length along field line
    gkyl_array_set_offset(up->geo_corn.mc2p_deflated, 1.0, up->geo_corn.mc2p, 1 * up->basis.num_basis);
    gkyl_array_set_offset(up->geo_corn.mc2nu_pos_deflated, 1.0, up->geo_corn.mc2nu_pos, 2 * up->basis.num_basis);
  }
  else if (up->grid.ndim==2) {
    // In 2D geometry, make mapc2p a function of only R and Z 
    // and mc2nu_pos only a function of psi and length along field line
    struct gkyl_array *temp = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
    gkyl_array_set_offset(up->geo_corn.mc2p_deflated     , 1.0, up->geo_corn.mc2p     , 0*up->basis.num_basis);
    gkyl_array_set_offset(temp                           , 1.0, up->geo_corn.mc2p     , 1*up->basis.num_basis);
    gkyl_array_set_offset(up->geo_corn.mc2p_deflated     , 1.0, temp                  , 1*up->basis.num_basis);
    gkyl_array_set_offset(up->geo_corn.mc2nu_pos_deflated, 1.0, up->geo_corn.mc2nu_pos, 0*up->basis.num_basis);
    gkyl_array_set_offset(temp                           , 1.0, up->geo_corn.mc2nu_pos, 2*up->basis.num_basis);
    gkyl_array_set_offset(up->geo_corn.mc2nu_pos_deflated, 1.0, temp                  , 1*up->basis.num_basis);
    gkyl_array_release(temp);
  }
  else if (up->grid.ndim==3) {
    // In 3D geoemtry, these are identical
    gkyl_array_copy(up->geo_corn.mc2p_deflated, up->geo_corn.mc2p);
    gkyl_array_copy(up->geo_corn.mc2nu_pos_deflated, up->geo_corn.mc2nu_pos);
  }
  // Deflate surface geo
  int count = 0;
  for (int dir = 0; dir < 3; dir++) {
    if(rem_dirs[dir] == 0) {
      struct gkyl_range local_ext_in_dir_3d;
      int lower_3d[3] = {up_3d->local.lower[0], up_3d->local.lower[1], up_3d->local.lower[2]};
      int upper_3d[3] = {up_3d->local.upper[0], up_3d->local.upper[1], up_3d->local.upper[2]};
      upper_3d[dir]+=1;
      gkyl_sub_range_init(&local_ext_in_dir_3d, &up_3d->local_ext, lower_3d, upper_3d);

      struct gkyl_range local_ext_in_dir;
      int lower[up->grid.ndim];
      int upper[up->grid.ndim];
      for(int j=0; j<up->grid.ndim; j++) {
        lower[j] = up->local.lower[j];
        upper[j] = up->local.upper[j];
      }
      upper[count]+=1;
      gkyl_sub_range_init(&local_ext_in_dir, &up->local_ext, lower, upper);

      struct gkyl_deflate_geo_surf* deflator_surf = gkyl_deflate_geo_surf_new(&up_3d->surf_basis, up->num_surf_basis, &up_3d->grid, &up->grid, rem_dirs, count, false);
      gkyl_deflate_geo_surf_advance(deflator_surf, &local_ext_in_dir_3d, &local_ext_in_dir, up_3d->geo_surf[dir].bmag, up->geo_surf[count].bmag, 1);
      gkyl_deflate_geo_surf_advance(deflator_surf, &local_ext_in_dir_3d, &local_ext_in_dir, up_3d->geo_surf[dir].jacobgeo, up->geo_surf[count].jacobgeo, 1);
      gkyl_deflate_geo_surf_advance(deflator_surf, &local_ext_in_dir_3d, &local_ext_in_dir, up_3d->geo_surf[dir].jacobgeo_ratio, up->geo_surf[count].jacobgeo_ratio, 1);
      gkyl_deflate_geo_surf_advance(deflator_surf, &local_ext_in_dir_3d, &local_ext_in_dir, up_3d->geo_surf[dir].jacobtot_inv, up->geo_surf[count].jacobtot_inv, 1);
      gkyl_deflate_geo_surf_advance(deflator_surf, &local_ext_in_dir_3d, &local_ext_in_dir, up_3d->geo_surf[dir].b_i, up->geo_surf[count].b_i, 3);
      gkyl_deflate_geo_surf_advance(deflator_surf, &local_ext_in_dir_3d, &local_ext_in_dir, up_3d->geo_surf[dir].cmag, up->geo_surf[count].cmag, 1);
      gkyl_deflate_geo_surf_advance(deflator_surf, &local_ext_in_dir_3d, &local_ext_in_dir, up_3d->geo_surf[dir].B3, up->geo_surf[count].B3, 1);
      gkyl_deflate_geo_surf_advance(deflator_surf, &local_ext_in_dir_3d, &local_ext_in_dir, up_3d->geo_surf[dir].normcurlbhat, up->geo_surf[count].normcurlbhat, 1);
      gkyl_deflate_geo_surf_advance(deflator_surf, &local_ext_in_dir_3d, &local_ext_in_dir, up_3d->geo_surf[dir].normals, up->geo_surf[count].normals, 9);
      gkyl_deflate_geo_surf_advance(deflator_surf, &local_ext_in_dir_3d, &local_ext_in_dir, up_3d->geo_surf[dir].lenr, up->geo_surf[count].lenr, 1);
      gkyl_deflate_geo_surf_advance(deflator_surf, &local_ext_in_dir_3d, &local_ext_in_dir, up_3d->geo_surf[dir].deltats, up->geo_surf[count].deltats, 1);
      gkyl_deflate_geo_surf_advance(deflator_surf, &local_ext_in_dir_3d, &local_ext_in_dir, up_3d->geo_surf[dir].bimpactangle, up->geo_surf[count].bimpactangle, 1);
      // deflate nodal quantities 
      gkyl_deflate_geo_surf_advance_nodal(deflator_surf, &up_3d->nrange_surf[dir], &up->nrange_surf[count], up_3d->geo_surf[dir].bmag_nodal, up->geo_surf[count].bmag_nodal, 1);
      gkyl_deflate_geo_surf_advance_nodal(deflator_surf, &up_3d->nrange_surf[dir], &up->nrange_surf[count], up_3d->geo_surf[dir].jacobgeo_nodal, up->geo_surf[count].jacobgeo_nodal, 1);
      gkyl_deflate_geo_surf_advance_nodal(deflator_surf, &up_3d->nrange_surf[dir], &up->nrange_surf[count], up_3d->geo_surf[dir].B3_nodal, up->geo_surf[count].B3_nodal, 1);
      gkyl_deflate_geo_surf_advance_nodal(deflator_surf, &up_3d->nrange_surf[dir], &up->nrange_surf[count], up_3d->geo_surf[dir].normcurlbhat_nodal, up->geo_surf[count].normcurlbhat_nodal, 1);
      gkyl_deflate_geo_surf_advance_nodal(deflator_surf, &up_3d->nrange_surf[dir], &up->nrange_surf[count], up_3d->geo_surf[dir].b_i_nodal, up->geo_surf[count].b_i_nodal, 3);
      gkyl_deflate_geo_surf_advance_nodal(deflator_surf, &up_3d->nrange_surf[dir], &up->nrange_surf[count], up_3d->geo_surf[dir].normals_nodal, up->geo_surf[count].normals_nodal, 9);
      gkyl_deflate_geo_surf_advance_nodal(deflator_surf, &up_3d->nrange_surf[dir], &up->nrange_surf[count], up_3d->geo_surf[dir].lenr_nodal, up->geo_surf[count].lenr_nodal, 1);
      gkyl_deflate_geo_surf_advance_nodal(deflator_surf, &up_3d->nrange_surf[dir], &up->nrange_surf[count], up_3d->geo_surf[dir].deltats_nodal, up->geo_surf[count].deltats_nodal, 1);
      gkyl_deflate_geo_surf_advance_nodal(deflator_surf, &up_3d->nrange_surf[dir], &up->nrange_surf[count], up_3d->geo_surf[dir].bimpactangle_nodal, up->geo_surf[count].bimpactangle_nodal, 1);
      count+=1;
      gkyl_deflate_geo_surf_release(deflator_surf);
    }
  }

  // Copy metadata.
  up->io_meta_basic = gkyl_msgpack_map_elem_clone(up_3d->io_meta_basic_len, up_3d->io_meta_basic);
  up->io_meta_basic_len = up_3d->io_meta_basic_len;
 
  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->ref_count = gkyl_ref_count_init(gkyl_gk_geometry_free);
  up->on_dev = up; // CPU eqn obj points to itself

  return up;
}

void 
gkyl_gk_geometry_write_efit(struct gkyl_gk_geometry_inp *geometry_inp,
  struct gkyl_msgpack_map_elem* io_meta_basic_basic, int io_meta_basic_basic_len)
{
  struct gkyl_efit *efit = gkyl_efit_new(&geometry_inp->efit_info);
  const char *fmt = "%s-psi.gkyl";
  int sz = gkyl_calc_strlen(fmt, efit->name);
  char fileNm[sz+1];
  snprintf(fileNm, sizeof fileNm, fmt, efit->name);

  // Set extra metadata
  char geqdsk_file_name[128];
  get_filename_from_path(efit->filepath, geqdsk_file_name, sizeof(geqdsk_file_name));
  struct gkyl_msgpack_map_elem io_meta_basic_rz[] = {
    { .key = "poly_order", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = efit->rzbasis.poly_order},
    { .key = "basis_type", .elem_type = GKYL_MP_STRING, .cval = efit->rzbasis.id},
    { .key = "geqdsk_file", .elem_type = GKYL_MP_STRING, .cval = geqdsk_file_name},
    { .key = "psisep", .elem_type = GKYL_MP_DOUBLE, .dval = efit->psisep},
    { .key = "sibry", .elem_type = GKYL_MP_DOUBLE, .dval = efit->sibry},
    { .key = "simag", .elem_type = GKYL_MP_DOUBLE, .dval = efit->simag},
    { .key = "bcentr", .elem_type = GKYL_MP_DOUBLE, .dval = efit->bcentr},
    { .key = "current", .elem_type = GKYL_MP_DOUBLE, .dval = efit->current},
    { .key = "rmaxis", .elem_type = GKYL_MP_DOUBLE, .dval = efit->rmaxis},
    { .key = "rcentr", .elem_type = GKYL_MP_DOUBLE, .dval = efit->zmaxis},
    { .key = "rleft", .elem_type = GKYL_MP_DOUBLE, .dval = efit->rleft},
    { .key = "rdim", .elem_type = GKYL_MP_DOUBLE, .dval = efit->rdim},
    { .key = "zmaxis", .elem_type = GKYL_MP_DOUBLE, .dval = efit->zmaxis},
    { .key = "zmid", .elem_type = GKYL_MP_DOUBLE, .dval = efit->zmid},
    { .key = "zdim", .elem_type = GKYL_MP_DOUBLE, .dval = efit->zdim},
  };
  int io_meta_basic_rz_len = sizeof(io_meta_basic_rz)/sizeof(io_meta_basic_rz[0]);
  int io_meta_basic_len[] = {io_meta_basic_basic_len, io_meta_basic_rz_len};
  const struct gkyl_msgpack_map_elem* io_meta_basic[] = {io_meta_basic_basic, io_meta_basic_rz};
  struct gkyl_msgpack_data *mt = gkyl_msgpack_create_union(sizeof(io_meta_basic_len)/sizeof(int), io_meta_basic_len, io_meta_basic);
  gkyl_grid_sub_array_write(&efit->rzgrid, &efit->rzlocal, mt, efit->psizr, fileNm);

  gkyl_msgpack_data_release(mt);
}

void
gkyl_gk_geometry_free(const struct gkyl_ref_count *ref)
{
  struct gk_geometry *up = container_of(ref, struct gk_geometry, ref_count);
  gkyl_array_release(up->geo_corn.mc2p);
  gkyl_array_release(up->geo_corn.mc2nu_pos);
  gkyl_array_release(up->geo_corn.bmag);
  gkyl_array_release(up->geo_corn.bmag_inv);
  gkyl_array_release(up->geo_corn.mc2p_deflated);
  gkyl_array_release(up->geo_corn.mc2nu_pos_deflated);

  gkyl_array_release(up->geo_int.mc2p);
  gkyl_array_release(up->geo_int.bmag);
  gkyl_array_release(up->geo_int.g_ij);
  gkyl_array_release(up->geo_int.g_ij_neut);
  gkyl_array_release(up->geo_int.jacobgeo);
  gkyl_array_release(up->geo_int.jacobgeo_ghost);
  gkyl_array_release(up->geo_int.jacobgeo_inv);
  gkyl_array_release(up->geo_int.dxdz);
  gkyl_array_release(up->geo_int.dzdx);
  gkyl_array_release(up->geo_int.dualmag);
  gkyl_array_release(up->geo_int.normals);
  gkyl_array_release(up->geo_int.gij);
  gkyl_array_release(up->geo_int.gij_neut);
  gkyl_array_release(up->geo_int.b_i);
  gkyl_array_release(up->geo_int.bcart);
  gkyl_array_release(up->geo_int.cmag);
  gkyl_array_release(up->geo_int.jacobtot);
  gkyl_array_release(up->geo_int.jacobtot_inv);
  gkyl_array_release(up->geo_int.gxxj);
  gkyl_array_release(up->geo_int.gxyj);
  gkyl_array_release(up->geo_int.gyyj);
  gkyl_array_release(up->geo_int.gxzj);
  gkyl_array_release(up->geo_int.eps2);
  gkyl_array_release(up->geo_int.dualcurlbhat);
  gkyl_array_release(up->geo_int.dualcurlbhatoverB);
  gkyl_array_release(up->geo_int.rtg33inv);
  gkyl_array_release(up->geo_int.bioverJB);
  gkyl_array_release(up->geo_int.B3);
  gkyl_array_release(up->geo_int.qprofile);

  for (int dir = 0; dir < up->grid.ndim; dir++) {
    gkyl_array_release(up->geo_surf[dir].jacobgeo);
    gkyl_array_release(up->geo_surf[dir].jacobgeo_ratio);
    gkyl_array_release(up->geo_surf[dir].bmag);
    gkyl_array_release(up->geo_surf[dir].b_i);
    gkyl_array_release(up->geo_surf[dir].cmag);
    gkyl_array_release(up->geo_surf[dir].jacobtot_inv);
    gkyl_array_release(up->geo_surf[dir].B3);
    gkyl_array_release(up->geo_surf[dir].normcurlbhat);
    gkyl_array_release(up->geo_surf[dir].normals);
    gkyl_array_release(up->geo_surf[dir].lenr);
    gkyl_array_release(up->geo_surf[dir].deltats);
    gkyl_array_release(up->geo_surf[dir].bimpactangle);
  }

  // Release nodal data
  gk_geometry_corn_release_nodal(up);
  gk_geometry_int_release_nodal(up);
  for (int dir=0; dir<up->grid.ndim; ++dir)
    gk_geometry_surf_release_nodal(up, dir);

  gkyl_msgpack_map_elem_release(up->io_meta_basic_len, up->io_meta_basic); 

  if (gkyl_gk_geometry_is_cu_dev(up)) 
    gkyl_cu_free(up->on_dev); 

  gkyl_free(up);
}

struct gk_geometry*
gkyl_gk_geometry_acquire(const struct gk_geometry* up)
{
  gkyl_ref_count_inc(&up->ref_count);
  return (struct gk_geometry*) up;
}

void
gkyl_gk_geometry_release(const struct gk_geometry *up)
{
  gkyl_ref_count_dec(&up->ref_count);
}



