#include <assert.h>
#include <math.h>

#include <gkyl_deflated_fem_poisson.h>
#include <gkyl_deflated_fem_poisson_priv.h>

// allocate array (filled with zeros)
static inline struct gkyl_array *
mkarr(bool use_gpu, long nc, long size)
{
  struct gkyl_array *a = use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size) :
                                   gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

struct gkyl_deflated_fem_poisson *
gkyl_deflated_fem_poisson_new(
  struct gkyl_rect_grid grid, struct gkyl_basis *basis_on_dev, struct gkyl_basis basis,
  struct gkyl_range local, struct gkyl_range global_sub_range, struct gkyl_array *epsilon,
  struct gkyl_array *kSq, struct gkyl_poisson_bc poisson_bc,
  struct gkyl_poisson_bias_plane_list *bias_plane_list, bool use_gpu
)
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
  up->num_solves_z = up->local.upper[up->cdim - 1] - up->local.lower[up->cdim - 1] + 2;
  up->d_fem_data = gkyl_malloc(sizeof(struct deflated_fem_data[up->num_solves_z]));

  // Check if one of the boundaries needs a spatially varying Dirichlet BC.
  up->isdirichletvar = false;
  for (int d = 0; d < up->cdim; d++) {
    up->isdirichletvar = up->isdirichletvar ||
                         (poisson_bc.lo_type[d] == GKYL_POISSON_DIRICHLET_VARYING ||
                          poisson_bc.up_type[d] == GKYL_POISSON_DIRICHLET_VARYING);
  }

  int poly_order = up->basis.poly_order;

  // Create 2d/3d nodal range nodal array to be populated
  int nodes[GKYL_MAX_DIM];
  if (poly_order == 1) {
    for (int d = 0; d < up->cdim; ++d) {
      nodes[d] = gkyl_range_shape(&up->local, d) + 1;
    }
  }
  if (poly_order == 2) {
    for (int d = 0; d < up->cdim; ++d) {
      nodes[d] = 2 * gkyl_range_shape(&up->local, d) + 1;
    }
  }
  gkyl_range_init_from_shape(&up->nrange, up->cdim, nodes);

  // Create deflated 1d/2d grid, ranges, basis, and nodal range
  double deflated_lower[GKYL_MAX_DIM] = {0.0};
  double deflated_upper[GKYL_MAX_DIM] = {0.0};
  int deflated_cells[GKYL_MAX_DIM] = {0};
  for (int i = 0; i < up->cdim - 1; i++) {
    deflated_lower[i] = up->grid.lower[i];
    deflated_upper[i] = up->grid.upper[i];
    deflated_cells[i] = up->grid.cells[i];
  }

  gkyl_rect_grid_init(
    &up->deflated_grid, up->cdim - 1, deflated_lower, deflated_upper, deflated_cells
  );
  int deflated_nghost[GKYL_MAX_CDIM] = {1};
  gkyl_create_grid_ranges(
    &up->deflated_grid, deflated_nghost, &up->deflated_local_ext, &up->deflated_local
  );
  gkyl_cart_modal_serendip(&up->deflated_basis, up->cdim - 1, poly_order);

  if (up->use_gpu) {
    up->deflated_basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
    gkyl_cart_modal_serendip_cu_dev(up->deflated_basis_on_dev, up->cdim - 1, poly_order);
  } else {
    up->deflated_basis_on_dev = &up->deflated_basis;
  }
  up->nodal_fld = mkarr(up->use_gpu, up->cdim, up->nrange.volume);
  // Buffers of gkyl_deflated_fem_poisson_advance_dg (one coefficient per cell of the extended
  // range, on which epsilon is defined).
  up->diff_even = mkarr(up->use_gpu, 1, epsilon->size);
  up->diff_odd = mkarr(up->use_gpu, 1, epsilon->size);

  int deflated_nodes[up->deflated_grid.ndim];
  if (poly_order == 1) {
    for (int d = 0; d < up->deflated_grid.ndim; ++d) {
      deflated_nodes[d] = up->deflated_grid.cells[d] + 1;
    }
  }
  if (poly_order == 2) {
    for (int d = 0; d < up->deflated_grid.ndim; ++d) {
      deflated_nodes[d] = 2 * (up->deflated_grid.cells[d]) + 1;
    }
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
  for (int zidx = up->local.lower[up->cdim - 1]; zidx <= up->local.upper[up->cdim - 1] + 1;
       zidx++) {
    int defl_num_basis = up->deflated_basis.num_basis;
    up->d_fem_data[ctr].deflated_rhs =
      mkarr(up->use_gpu, defl_num_basis, up->deflated_local_ext.volume);
    up->d_fem_data[ctr].deflated_phi =
      mkarr(up->use_gpu, defl_num_basis, up->deflated_local_ext.volume);
    up->d_fem_data[ctr].deflated_epsilon = mkarr(
      up->use_gpu, (2 * up->deflated_grid.ndim - 1) * defl_num_basis, up->deflated_local_ext.volume
    );
    up->d_fem_data[ctr].deflated_nodal_fld =
      mkarr(up->use_gpu, up->deflated_grid.ndim, up->deflated_nrange.volume);
    up->d_fem_data[ctr].deflated_phibc =
      up->isdirichletvar ? mkarr(up->use_gpu, defl_num_basis, up->deflated_local_ext.volume) : 0;
    up->d_fem_data[ctr].deflated_kSq =
      up->ishelmholtz ? mkarr(up->use_gpu, defl_num_basis, up->deflated_local_ext.volume) : 0;

    if (zidx == up->local.upper[up->cdim - 1] + 1) {
      gkyl_deflate_zsurf_advance(
        up->deflator_up, zidx - 1, &up->local, &up->deflated_local, epsilon,
        up->d_fem_data[ctr].deflated_epsilon, 2 * up->deflated_grid.ndim - 1
      );
      if (up->ishelmholtz) {
        gkyl_deflate_zsurf_advance(
          up->deflator_up, zidx - 1, &up->local, &up->deflated_local, kSq,
          up->d_fem_data[ctr].deflated_kSq, 1
        );
      }
    } else {
      gkyl_deflate_zsurf_advance(
        up->deflator_lo, zidx, &up->local, &up->deflated_local, epsilon,
        up->d_fem_data[ctr].deflated_epsilon, 2 * up->deflated_grid.ndim - 1
      );
      if (up->ishelmholtz) {
        gkyl_deflate_zsurf_advance(
          up->deflator_lo, zidx, &up->local, &up->deflated_local, kSq,
          up->d_fem_data[ctr].deflated_kSq, 1
        );
      }
    }

    // we check if we are at an extremal global index of z to apply the bias plane at extremal z values only
    int global_zidx = zidx + up->global_sub_range.lower[up->cdim - 1];
    // check if the global index is equal to the lower or upper limit of the z grid
    bool is_lower_z_edge = (zidx == up->local.lower[up->cdim - 1]) &&
                           up->poisson_bc.contains_lower_z_edge;
    bool is_upper_z_edge = (zidx == up->local.upper[up->cdim - 1] + 1) &&
                           up->poisson_bc.contains_upper_z_edge;
    // store it in pisson_bc to pass it to gkyl_fem_poisson_new
    struct gkyl_poisson_bias_plane_list *bias_plane_list;
    if (is_upper_z_edge || is_lower_z_edge) {
      bias_plane_list = up->bias_plane_list;
    } else {
      bias_plane_list = NULL;
    }

    up->d_fem_data[ctr].fem_poisson = gkyl_fem_poisson_new(
      &up->deflated_local, &up->deflated_grid, up->deflated_basis, &up->poisson_bc, bias_plane_list,
      up->d_fem_data[ctr].deflated_epsilon, up->d_fem_data[ctr].deflated_kSq, false, use_gpu
    );
    ctr += 1;
  }

  return up;
}

static void
deflated_fem_poisson_solve_plane(
  struct gkyl_deflated_fem_poisson *up, struct gkyl_deflate_zsurf *deflator, int zidx, int plane,
  struct gkyl_array *rhs, struct gkyl_array *phibc
)
{
  // Solve on the z plane with index plane, with the RHS (and the BC field) evaluated at the
  // lower or upper (depending on the deflator) z surface of the cells with z index zidx, and
  // store the solution in that plane of the nodal field.
  struct deflated_fem_data *plane_data = &up->d_fem_data[plane];
  // Deflate rhs indexing global sub-range to fetch correct place in z
  gkyl_deflate_zsurf_advance(
    deflator, zidx, &up->global_sub_range, &up->deflated_local, rhs, plane_data->deflated_rhs, 1
  );
  if (up->isdirichletvar) {
    // Deflate the BC field.
    gkyl_deflate_zsurf_advance(
      deflator, zidx, &up->global_sub_range, &up->deflated_local, phibc, plane_data->deflated_phibc,
      1
    );
  }
  // Do the poisson solve
  gkyl_fem_poisson_set_rhs(
    plane_data->fem_poisson, plane_data->deflated_rhs, plane_data->deflated_phibc
  );
  gkyl_fem_poisson_solve(plane_data->fem_poisson, plane_data->deflated_phi);
  // Modal to Nodal in 1d -> Store the result in the 2d nodal field
  gkyl_nodal_ops_m2n_deflated(
    up->n2m_deflated, up->deflated_basis_on_dev, &up->deflated_grid, &up->nrange,
    &up->deflated_nrange, &up->deflated_local, 1, up->nodal_fld, plane_data->deflated_phi, plane
  );
}

void
gkyl_deflated_fem_poisson_advance(
  struct gkyl_deflated_fem_poisson *up, struct gkyl_array *rhs, struct gkyl_array *phibc,
  struct gkyl_array *phi
)
{
  int zidx_lo = up->global_sub_range.lower[up->cdim - 1];
  int zidx_up = up->global_sub_range.upper[up->cdim - 1];
  for (int zidx = zidx_lo; zidx <= zidx_up; zidx++) {
    deflated_fem_poisson_solve_plane(up, up->deflator_lo, zidx, zidx - zidx_lo, rhs, phibc);
  }
  deflated_fem_poisson_solve_plane(up, up->deflator_up, zidx_up, zidx_up - zidx_lo + 1, rhs, phibc);
  gkyl_nodal_ops_n2m(
    up->n2m, up->basis_on_dev, &up->grid, &up->nrange, &up->local, 1, up->nodal_fld, phi, false
  );
}

void
gkyl_deflated_fem_poisson_advance_dg(
  struct gkyl_deflated_fem_poisson *up, struct gkyl_array *rhs, struct gkyl_array *phibc,
  struct gkyl_array *phi, struct gkyl_array *buffer
)
{
  // In each cell the RHS is r = r_lo*h_lo(z) + r_up*h_up(z), with h_lo/h_up the linear functions
  // of z that are 1 at the lower/upper z surface of the cell and 0 at the other, and the solution
  // is S_lo(r_lo)*h_lo(z) + S_up(r_up)*h_up(z), S_lo/S_up being the solves on those surfaces.
  assert(up->basis.poly_order == 1);

  // Field continuous along z whose value at the lower z surface of each cell is S_lo(r_lo).
  gkyl_deflated_fem_poisson_advance(up, rhs, phibc, phi);

  // Field continuous along z whose value at the upper z surface of each cell is S_up(r_up)
  // (the lowest plane keeps the solution above, which is not used).
  int zidx_lo = up->global_sub_range.lower[up->cdim - 1];
  int zidx_up = up->global_sub_range.upper[up->cdim - 1];
  for (int zidx = zidx_lo; zidx <= zidx_up; zidx++) {
    deflated_fem_poisson_solve_plane(up, up->deflator_up, zidx, zidx - zidx_lo + 1, rhs, phibc);
  }
  gkyl_nodal_ops_n2m(
    up->n2m, up->basis_on_dev, &up->grid, &up->nrange, &up->local, 1, up->nodal_fld, buffer, false
  );

  // Keep the lower surface value of the first field and take the upper surface value of the
  // second one. A pair of coefficients (even, odd) of the basis functions that are constant and
  // linear in z has the surface values even -/+ sqrt(3)*odd (up to a factor), so with
  // diff = second - first: even += diff_even/2 + sqrt(3)*diff_odd/2, and
  // odd += diff_even/(2*sqrt(3)) + diff_odd/2.
  gkyl_array_accumulate_range(buffer, -1.0, phi, &up->local);
  int num_basis = up->basis.num_basis;
  // Coefficients of the basis functions linear in z, and of their counterparts constant in z.
  static const int odd_2x[] = {2, 3}, even_2x[] = {0, 1};
  static const int odd_3x[] = {3, 5, 6, 7}, even_3x[] = {0, 1, 2, 4};
  const int *odd = up->cdim == 2 ? odd_2x : odd_3x;
  const int *even = up->cdim == 2 ? even_2x : even_3x;
  for (int i = 0; i < num_basis / 2; i++) {
    gkyl_array_set_offset_ncoeff_range(up->diff_even, 1.0, buffer, 1, even[i], 0, &up->local);
    gkyl_array_set_offset_ncoeff_range(up->diff_odd, 1.0, buffer, 1, odd[i], 0, &up->local);
    gkyl_array_accumulate_offset_range(phi, 0.5, up->diff_even, even[i], &up->local);
    gkyl_array_accumulate_offset_range(phi, 0.5 * sqrt(3.0), up->diff_odd, even[i], &up->local);
    gkyl_array_accumulate_offset_range(phi, 0.5 / sqrt(3.0), up->diff_even, odd[i], &up->local);
    gkyl_array_accumulate_offset_range(phi, 0.5, up->diff_odd, odd[i], &up->local);
  }
}

void
gkyl_deflated_fem_poisson_release(struct gkyl_deflated_fem_poisson *up)
{
  gkyl_array_release(up->nodal_fld);
  gkyl_array_release(up->diff_even);
  gkyl_array_release(up->diff_odd);
  gkyl_nodal_ops_release(up->n2m);
  gkyl_nodal_ops_release(up->n2m_deflated);
  gkyl_deflate_zsurf_release(up->deflator_lo);
  gkyl_deflate_zsurf_release(up->deflator_up);
  int ctr = 0;
  for (int zidx = up->local.lower[up->cdim - 1]; zidx <= up->local.upper[up->cdim - 1] + 1;
       zidx++) {
    gkyl_array_release(up->d_fem_data[ctr].deflated_rhs);
    gkyl_array_release(up->d_fem_data[ctr].deflated_phi);
    gkyl_array_release(up->d_fem_data[ctr].deflated_epsilon);
    if (up->isdirichletvar) {
      gkyl_array_release(up->d_fem_data[ctr].deflated_phibc);
    }
    if (up->ishelmholtz) {
      gkyl_array_release(up->d_fem_data[ctr].deflated_kSq);
    }
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
