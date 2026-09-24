#include <assert.h>

#include <gkyl_deflated_fem_poisson.h>
#include <gkyl_deflated_fem_poisson_priv.h>

// allocate array (filled with zeros)
static inline struct gkyl_array*
mkarr(bool use_gpu, long nc, long size)
{
  struct gkyl_array* a = use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
                                : gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

struct gkyl_deflated_fem_poisson* 
gkyl_deflated_fem_poisson_new(struct gkyl_rect_grid grid, struct gkyl_basis *basis_on_dev, struct gkyl_basis basis, 
  struct gkyl_range local, struct gkyl_range global_sub_range, struct gkyl_array *epsilon, struct gkyl_array *kSq,
  struct gkyl_poisson_bc poisson_bc, struct gkyl_poisson_bias_plane_list *bias_plane_list, bool use_gpu)
{
  struct gkyl_deflated_fem_poisson *up = gkyl_malloc(sizeof(*up));
  up->use_gpu = use_gpu;
  up->grid = grid;
  up->basis = basis;
  up->basis_on_dev = basis_on_dev;

  // Assumes domain decomposition *only* in z.
  // local is range of potential; global_sub_range is range for indexing global charge density
  up->local = local;
  up->global_sub_range = global_sub_range;
  assert(up->local.volume == up->global_sub_range.volume);

  up->poisson_bc = poisson_bc;
  up->bias_plane_list = bias_plane_list;
  up->cdim = grid.ndim;
  up->num_solves_z = up->local.upper[up->cdim-1] - up->local.lower[up->cdim-1] + 2;
  up->d_fem_data = gkyl_malloc(sizeof(struct deflated_fem_data[up->num_solves_z]));

  // Check if one of the boundaries needs a spatially varying Dirichlet BC.
  up->isdirichletvar = false;
  for (int d=0; d<up->cdim; d++) up->isdirichletvar = up->isdirichletvar ||
    (poisson_bc.lo_type[d] == GKYL_POISSON_DIRICHLET_VARYING ||
     poisson_bc.up_type[d] == GKYL_POISSON_DIRICHLET_VARYING);

  int poly_order = up->basis.poly_order;

  // Create 2d/3d nodal range nodal array to be populated
  int nodes[GKYL_MAX_DIM];
  if (poly_order == 1) {
    for (int d=0; d<up->cdim; ++d)
      nodes[d] = gkyl_range_shape(&up->local, d) + 1;
  }
  if (poly_order == 2) {
    for (int d=0; d<up->cdim; ++d)
      nodes[d] = 2*gkyl_range_shape(&up->local, d) + 1;
  }
  gkyl_range_init_from_shape(&up->nrange, up->cdim, nodes);

  // Create deflated 1d/2d grid, ranges, basis, and nodal range
  double deflated_lower[GKYL_MAX_DIM] = { 0.0 } ;
  double deflated_upper[GKYL_MAX_DIM] = { 0.0 };
  int deflated_cells[GKYL_MAX_DIM] = { 0 };
  for(int i = 0; i < up->cdim-1; i++){
    deflated_lower[i] = up->grid.lower[i];
    deflated_upper[i] = up->grid.upper[i];
    deflated_cells[i] = up->grid.cells[i];
  }

  gkyl_rect_grid_init(&up->deflated_grid, up->cdim-1, deflated_lower, deflated_upper, deflated_cells);
  int deflated_nghost[GKYL_MAX_CDIM] = { 1 };
  gkyl_create_grid_ranges(&up->deflated_grid, deflated_nghost, &up->deflated_local_ext, &up->deflated_local);
  gkyl_cart_modal_serendip(&up->deflated_basis, up->cdim-1, poly_order);

  if (up->use_gpu) {
    up->deflated_basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
    gkyl_cart_modal_serendip_cu_dev(up->deflated_basis_on_dev, up->cdim-1, poly_order);
  }
  else {
    up->deflated_basis_on_dev = &up->deflated_basis;
  }
  up->nodal_fld = mkarr(up->use_gpu, up->cdim, up->nrange.volume);

  int deflated_nodes[up->deflated_grid.ndim];
  if (poly_order == 1){
    for (int d=0; d<up->deflated_grid.ndim; ++d)
      deflated_nodes[d] = up->deflated_grid.cells[d] + 1;
  }
  if (poly_order == 2){
    for (int d=0; d<up->deflated_grid.ndim; ++d)
      deflated_nodes[d] = 2*(up->deflated_grid.cells[d]) + 1;
  }
  gkyl_range_init_from_shape(&up->deflated_nrange, up->deflated_grid.ndim, deflated_nodes);

  // Allocate a deflation object for the lower and upper edge
  up->deflator_lo = gkyl_deflate_zsurf_new(&up->basis, &up->deflated_basis, 0, use_gpu);
  up->deflator_up = gkyl_deflate_zsurf_new(&up->basis, &up->deflated_basis, 1, use_gpu);

  // Allocate the 1d and 2d nodal_ops updater to be used in the advance method
  up->n2m = gkyl_nodal_ops_new(&up->basis, &up->grid, use_gpu);
  up->n2m_deflated = gkyl_nodal_ops_new(&up->deflated_basis, &up->deflated_grid, use_gpu);

  if (kSq) {
    up->ishelmholtz = true;
  } else {
    up->ishelmholtz = false;
  }

  // Allocate necessary fields and solvers for each z slice
  int ctr = 0;
  for (int zidx = up->local.lower[up->cdim-1]; zidx <= up->local.upper[up->cdim-1]+1; zidx++) {
    int defl_num_basis = up->deflated_basis.num_basis;
    up->d_fem_data[ctr].deflated_rhs     = mkarr(up->use_gpu, defl_num_basis, up->deflated_local_ext.volume);
    up->d_fem_data[ctr].deflated_phi     = mkarr(up->use_gpu, defl_num_basis, up->deflated_local_ext.volume);
    up->d_fem_data[ctr].deflated_epsilon = mkarr(up->use_gpu, (2*up->deflated_grid.ndim-1)*defl_num_basis, up->deflated_local_ext.volume);
    up->d_fem_data[ctr].deflated_nodal_fld = mkarr(up->use_gpu, up->deflated_grid.ndim, up->deflated_nrange.volume);
    up->d_fem_data[ctr].deflated_phibc = up->isdirichletvar?
      mkarr(up->use_gpu, defl_num_basis, up->deflated_local_ext.volume) : 0;
    up->d_fem_data[ctr].deflated_kSq = up->ishelmholtz?
      mkarr(up->use_gpu, (2*up->deflated_grid.ndim-1)*defl_num_basis, up->deflated_local_ext.volume) : 0;

    if (zidx == up->local.upper[up->cdim-1] + 1 ) {
      gkyl_deflate_zsurf_advance(up->deflator_up, zidx-1, &up->local, &up->deflated_local,
        epsilon, up->d_fem_data[ctr].deflated_epsilon,  2*up->deflated_grid.ndim-1);
      if (up->ishelmholtz)
        gkyl_deflate_zsurf_advance(up->deflator_up, zidx-1, &up->local, &up->deflated_local,
          kSq, up->d_fem_data[ctr].deflated_kSq,  2*up->deflated_grid.ndim-1);
    }
    else {
      gkyl_deflate_zsurf_advance(up->deflator_lo, zidx, &up->local, &up->deflated_local,
        epsilon, up->d_fem_data[ctr].deflated_epsilon, 2*up->deflated_grid.ndim-1);
      if (up->ishelmholtz)
        gkyl_deflate_zsurf_advance(up->deflator_lo, zidx, &up->local, &up->deflated_local,
          kSq, up->d_fem_data[ctr].deflated_kSq, 2*up->deflated_grid.ndim-1);
    }

    // we check if we are at an extremal global index of z to apply the bias plane at extremal z values only
    int global_zidx = zidx + up->global_sub_range.lower[up->cdim-1]; 
    // check if the global index is equal to the lower or upper limit of the z grid
    bool is_lower_z_edge = (zidx == up->local.lower[up->cdim-1])   && up->poisson_bc.contains_lower_z_edge;
    bool is_upper_z_edge = (zidx == up->local.upper[up->cdim-1]+1) && up->poisson_bc.contains_upper_z_edge;
    // store it in pisson_bc to pass it to gkyl_fem_poisson_new
    struct gkyl_poisson_bias_plane_list *bias_plane_list;
    if (is_upper_z_edge || is_lower_z_edge) {
      bias_plane_list = up->bias_plane_list;
    } else {
      bias_plane_list = NULL;
    }

    up->d_fem_data[ctr].fem_poisson = gkyl_fem_poisson_new(&up->deflated_local, &up->deflated_grid,
      up->deflated_basis, &up->poisson_bc, bias_plane_list, up->d_fem_data[ctr].deflated_epsilon,
      up->d_fem_data[ctr].deflated_kSq, false, use_gpu);
    ctr += 1;
  }

  return up;
}

void 
gkyl_deflated_fem_poisson_advance(struct gkyl_deflated_fem_poisson *up, struct gkyl_array *rhs,
  struct gkyl_array *phibc, struct gkyl_array* phi)
{
  int ctr = 0;
  int local_range_ctr = up->local.lower[up->cdim-1];
  for (int zidx = up->global_sub_range.lower[up->cdim-1]; zidx <= up->global_sub_range.upper[up->cdim-1]; zidx++) {
    // Deflate rhs indexing global sub-range to fetch correct place in z
    gkyl_deflate_zsurf_advance(up->deflator_lo, zidx, 
      &up->global_sub_range, &up->deflated_local, rhs, up->d_fem_data[ctr].deflated_rhs, 1);
    if (up->isdirichletvar) {
      // Deflate the BC field.
      gkyl_deflate_zsurf_advance(up->deflator_lo, zidx, 
        &up->global_sub_range, &up->deflated_local, phibc, up->d_fem_data[ctr].deflated_phibc, 1);
    }
    // Do the poisson solve 
    gkyl_fem_poisson_set_rhs(up->d_fem_data[ctr].fem_poisson, up->d_fem_data[ctr].deflated_rhs, up->d_fem_data[ctr].deflated_phibc);
    gkyl_fem_poisson_solve(up->d_fem_data[ctr].fem_poisson, up->d_fem_data[ctr].deflated_phi);
    // Modal to Nodal in 1d -> Store the result in the 2d nodal field
    gkyl_nodal_ops_m2n_deflated(up->n2m_deflated, up->deflated_basis_on_dev, 
      &up->deflated_grid, &up->nrange, &up->deflated_nrange, &up->deflated_local, 1, 
      up->nodal_fld, up->d_fem_data[ctr].deflated_phi, ctr);
    ctr += 1;
    local_range_ctr += 1;
    if (zidx == up->global_sub_range.upper[up->cdim-1]) {
      // Deflate rhs indexing global sub-range to fetch correct place in z
      gkyl_deflate_zsurf_advance(up->deflator_up, zidx, 
        &up->global_sub_range, &up->deflated_local, rhs, up->d_fem_data[ctr].deflated_rhs, 1);
      if (up->isdirichletvar) {
        // Deflate the BC field.
        gkyl_deflate_zsurf_advance(up->deflator_up, zidx, 
          &up->global_sub_range, &up->deflated_local, phibc, up->d_fem_data[ctr].deflated_phibc, 1);
      }
      // Do the poisson solve 
      gkyl_fem_poisson_set_rhs(up->d_fem_data[ctr].fem_poisson, up->d_fem_data[ctr].deflated_rhs, up->d_fem_data[ctr].deflated_phibc);
      gkyl_fem_poisson_solve(up->d_fem_data[ctr].fem_poisson, up->d_fem_data[ctr].deflated_phi);
      // Modal to Nodal in 1d -> Store the result in the 2d nodal rhs.
      gkyl_nodal_ops_m2n_deflated(up->n2m_deflated, up->deflated_basis_on_dev, 
        &up->deflated_grid, &up->nrange, &up->deflated_nrange, &up->deflated_local, 1, 
        up->nodal_fld, up->d_fem_data[ctr].deflated_phi, ctr);
    }
  }
  gkyl_nodal_ops_n2m(up->n2m, up->basis_on_dev, &up->grid, &up->nrange, &up->local, 1, up->nodal_fld, phi, false);
}

void gkyl_deflated_fem_poisson_release(struct gkyl_deflated_fem_poisson* up){
  gkyl_array_release(up->nodal_fld);
  gkyl_nodal_ops_release(up->n2m);
  gkyl_nodal_ops_release(up->n2m_deflated);
  gkyl_deflate_zsurf_release(up->deflator_lo);
  gkyl_deflate_zsurf_release(up->deflator_up);
  int ctr = 0;
  for (int zidx = up->local.lower[up->cdim-1]; zidx <= up->local.upper[up->cdim-1] + 1; zidx++) {
    gkyl_array_release(up->d_fem_data[ctr].deflated_rhs);
    gkyl_array_release(up->d_fem_data[ctr].deflated_phi);
    gkyl_array_release(up->d_fem_data[ctr].deflated_epsilon);
    if (up->isdirichletvar)
      gkyl_array_release(up->d_fem_data[ctr].deflated_phibc);
    if (up->ishelmholtz)
      gkyl_array_release(up->d_fem_data[ctr].deflated_kSq);
    gkyl_array_release(up->d_fem_data[ctr].deflated_nodal_fld);
    gkyl_fem_poisson_release(up->d_fem_data[ctr].fem_poisson);
    ctr += 1;
  }
  if (up->use_gpu) {
    gkyl_cu_free(up->deflated_basis_on_dev);
  }

  gkyl_free(up->d_fem_data);
  gkyl_free(up);
}


