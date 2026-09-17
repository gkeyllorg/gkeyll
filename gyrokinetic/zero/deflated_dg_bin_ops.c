#include <assert.h>

#include <gkyl_deflated_dg_bin_ops.h>
#include <gkyl_deflated_dg_bin_ops_priv.h>

#include <gkyl_array_rio.h>

struct gkyl_deflated_dg_bin_ops* 
gkyl_deflated_dg_bin_ops_new(struct gkyl_rect_grid grid, 
  struct gkyl_basis *basis_on_dev, struct gkyl_basis basis, 
  struct gkyl_range local, bool use_gpu)
{
  struct gkyl_deflated_dg_bin_ops *up = gkyl_malloc(sizeof(*up));
  up->use_gpu = use_gpu;
  up->grid = grid;
  up->basis = basis;
  up->basis_on_dev = basis_on_dev;

  // Assumes domain decomposition *only* in z.
  // local range and global sub_ranges can be different if the operation is being called on global arrays
  up->local = local;

  up->cdim = grid.ndim;
  up->num_solves_z = up->local.upper[up->cdim-1] - up->local.lower[up->cdim-1] + 2;
  up->d_bop_data = gkyl_malloc(sizeof(struct deflated_dg_bin_ops_data[up->num_solves_z]));

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

  gkyl_rect_grid_init(&up->deflated_grid, up->cdim-1, deflated_lower, 
    deflated_upper, deflated_cells);
  int deflated_nghost[GKYL_MAX_CDIM] = { 1 };
  gkyl_create_grid_ranges(&up->deflated_grid, deflated_nghost, 
    &up->deflated_local_ext, &up->deflated_local);
  gkyl_cart_modal_serendip(&up->deflated_basis, up->cdim-1, poly_order);

  if (up->use_gpu) {
    up->deflated_basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
    gkyl_cart_modal_serendip_cu_dev(up->deflated_basis_on_dev, up->cdim-1, poly_order);
    up->nodal_fld = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->cdim, up->nrange.volume);
  }
  else {
    up->deflated_basis_on_dev = &up->deflated_basis;
    up->nodal_fld = gkyl_array_new(GKYL_DOUBLE, up->cdim, up->nrange.volume);
  }

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

  // Allocate memory for division
  if (up->use_gpu)
    up->mem = gkyl_dg_bin_op_mem_cu_dev_new(up->deflated_local.volume, up->deflated_basis.num_basis);
  else
    up->mem = gkyl_dg_bin_op_mem_new(up->deflated_local.volume, up->deflated_basis.num_basis);

  // Allocate necessary fields and solvers for each z slice
  int ctr = 0;
  for (int zidx = up->local.lower[up->cdim-1]; zidx <= up->local.upper[up->cdim-1]+1; zidx++) {
    if (use_gpu) {
      up->d_bop_data[ctr].deflated_lop = gkyl_array_cu_dev_new(GKYL_DOUBLE, 
        up->deflated_basis.num_basis, up->deflated_local_ext.volume);
      up->d_bop_data[ctr].deflated_rop = gkyl_array_cu_dev_new(GKYL_DOUBLE, 
        up->deflated_basis.num_basis, up->deflated_local_ext.volume);
      up->d_bop_data[ctr].deflated_out = gkyl_array_cu_dev_new(GKYL_DOUBLE, 
        up->deflated_basis.num_basis, up->deflated_local_ext.volume);
    }
    else {
      up->d_bop_data[ctr].deflated_lop = gkyl_array_new(GKYL_DOUBLE, 
        up->deflated_basis.num_basis, up->deflated_local_ext.volume);
      up->d_bop_data[ctr].deflated_rop = gkyl_array_new(GKYL_DOUBLE, 
        up->deflated_basis.num_basis, up->deflated_local_ext.volume);
      up->d_bop_data[ctr].deflated_out = gkyl_array_new(GKYL_DOUBLE, 
        up->deflated_basis.num_basis, up->deflated_local_ext.volume);
    }
    ctr += 1;
  }

  return up;
}

void 
deflated_dg_bin_ops_advance(enum deflated_dg_bin_ops_type op_type, struct gkyl_deflated_dg_bin_ops* up,
  int c_oop, struct gkyl_array *out, int c_lop, struct gkyl_array *lop, int c_rop, struct gkyl_array* rop)
{

  int lop_ncomp = lop->ncomp/up->basis.num_basis;
  int rop_ncomp = lop->ncomp/up->basis.num_basis;
  int ctr = 0;
  for (int zidx = up->local.lower[up->cdim-1]; zidx <= up->local.upper[up->cdim-1]; zidx++) {
    // Deflate lop and rop
    gkyl_deflate_zsurf_advance(up->deflator_lo, zidx, 
      &up->local, &up->deflated_local, lop, up->d_bop_data[ctr].deflated_lop, lop_ncomp);
    gkyl_deflate_zsurf_advance(up->deflator_lo, zidx, 
      &up->local, &up->deflated_local, rop, up->d_bop_data[ctr].deflated_rop, rop_ncomp);

    // Divide or Multiply
    if (op_type == GKYL_DEFLATED_DIV)
      gkyl_dg_div_op_range(up->mem, &up->deflated_basis, c_oop, up->d_bop_data[ctr].deflated_out,
        c_lop, up->d_bop_data[ctr].deflated_lop, c_rop, up->d_bop_data[ctr].deflated_rop, &up->deflated_local);
    else if (op_type == GKYL_DEFLATED_MUL)
      gkyl_dg_mul_op_range(&up->deflated_basis, c_oop, up->d_bop_data[ctr].deflated_out, c_lop,
        up->d_bop_data[ctr].deflated_lop, c_rop, up->d_bop_data[ctr].deflated_rop, &up->deflated_local);

    // Modal to Nodal in 1d -> Store the result in the 2d nodal field
    gkyl_nodal_ops_m2n_deflated(up->n2m_deflated, up->deflated_basis_on_dev, 
      &up->deflated_grid, &up->nrange, &up->deflated_nrange, &up->deflated_local, 1, 
      up->nodal_fld, up->d_bop_data[ctr].deflated_out, ctr);
    ctr += 1;
    if (zidx == up->local.upper[up->cdim-1]) {
      // Deflate lop and rop
      gkyl_deflate_zsurf_advance(up->deflator_up, zidx, 
        &up->local, &up->deflated_local, lop, up->d_bop_data[ctr].deflated_lop, lop_ncomp);
      gkyl_deflate_zsurf_advance(up->deflator_up, zidx, 
        &up->local, &up->deflated_local, rop, up->d_bop_data[ctr].deflated_rop, rop_ncomp);


      // Divide or Multiply
      if (op_type == GKYL_DEFLATED_DIV)
        gkyl_dg_div_op_range(up->mem, &up->deflated_basis, c_oop, up->d_bop_data[ctr].deflated_out,
          c_lop, up->d_bop_data[ctr].deflated_lop, c_rop, up->d_bop_data[ctr].deflated_rop, &up->deflated_local);
      else if (op_type == GKYL_DEFLATED_MUL)
        gkyl_dg_mul_op_range(&up->deflated_basis, c_oop, up->d_bop_data[ctr].deflated_out,
          c_lop, up->d_bop_data[ctr].deflated_lop, c_rop, up->d_bop_data[ctr].deflated_rop, &up->deflated_local);

      // Modal to Nodal in 1d -> Store the result in the 2d nodal field
      gkyl_nodal_ops_m2n_deflated(up->n2m_deflated, up->deflated_basis_on_dev, 
        &up->deflated_grid, &up->nrange, &up->deflated_nrange, &up->deflated_local, 1, 
        up->nodal_fld, up->d_bop_data[ctr].deflated_out, ctr);
    }
  }

  gkyl_nodal_ops_n2m(up->n2m, up->basis_on_dev, &up->grid, &up->nrange, &up->local, 1, up->nodal_fld, out, false);

}


void 
gkyl_deflated_dg_bin_ops_mul(struct gkyl_deflated_dg_bin_ops* up, int c_oop, struct gkyl_array *out,
  int c_lop, struct gkyl_array *lop, int c_rop, struct gkyl_array* rop)
{
  deflated_dg_bin_ops_advance(GKYL_DEFLATED_MUL, up, c_oop, out, c_lop, lop, c_rop, rop);
}

void 
gkyl_deflated_dg_bin_ops_div(struct gkyl_deflated_dg_bin_ops* up, int c_oop, struct gkyl_array *out,
  int c_lop, struct gkyl_array *lop, int c_rop, struct gkyl_array* rop)
{
  deflated_dg_bin_ops_advance(GKYL_DEFLATED_DIV, up, c_oop, out, c_lop, lop, c_rop, rop);
}

void
gkyl_deflated_dg_bin_ops_release(struct gkyl_deflated_dg_bin_ops* up)
{
  gkyl_array_release(up->nodal_fld);
  gkyl_nodal_ops_release(up->n2m);
  gkyl_nodal_ops_release(up->n2m_deflated);
  gkyl_deflate_zsurf_release(up->deflator_lo);
  gkyl_deflate_zsurf_release(up->deflator_up);
  gkyl_dg_bin_op_mem_release(up->mem);
  int ctr = 0;
  for (int zidx = up->local.lower[up->cdim-1]; zidx <= up->local.upper[up->cdim-1] + 1; zidx++) {
    gkyl_array_release(up->d_bop_data[ctr].deflated_lop);
    gkyl_array_release(up->d_bop_data[ctr].deflated_rop);
    gkyl_array_release(up->d_bop_data[ctr].deflated_out);
    ctr += 1;
  }
  if (up->use_gpu) {
    gkyl_cu_free(up->deflated_basis_on_dev);
  }

  gkyl_free(up->d_bop_data);
  gkyl_free(up);
}


