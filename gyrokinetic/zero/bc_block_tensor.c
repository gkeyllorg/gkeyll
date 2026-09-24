#include <gkyl_bc_block_tensor.h>
#include <gkyl_bc_block_tensor_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_util.h>

static inline double dot_product(const double *v1, const double *v2)
{
  double out = 0.0;
  for(int i = 0; i < 3; i++)
    out += v1[i]*v2[i];
  return out;
}




struct bc_block_tensor*
gkyl_bc_block_tensor_new(const struct gkyl_rect_grid* grid, const struct gkyl_range *range, const struct gkyl_range* range_ext, 
  const struct gkyl_basis* basis, bool use_gpu)
{
  struct bc_block_tensor *up = gkyl_malloc(sizeof(* up));
  up->basis = *basis;
  up->range = *range;
  up->range_ext = *range_ext;
  up->grid = *grid;
  up->cdim = grid->ndim;
  up->poly_order = basis->poly_order;

  struct gkyl_basis surf_basis;
  gkyl_cart_modal_tensor(&surf_basis, up->cdim-1, up->poly_order);
  up->num_surf_nodes = surf_basis.num_basis;
  up->tensor = gkyl_array_new(GKYL_DOUBLE, up->cdim*up->cdim*up->num_surf_nodes, up->range_ext.volume);

  return up;
}


void calc_tensor(struct bc_block_tensor *up, int dir, int edge1, int edge2, const double *ej, const double *e_i, double *tj_i)
{
  // First evaluate at all the quadrature nodes
  double ej_surf[up->num_surf_nodes][9];
  double e_i_surf[up->num_surf_nodes][9];
  for(int n = 0; n < up->num_surf_nodes; n++) {
    for(int i = 0; i < 9; i++) {
      e_i_surf[n][i] = bc_block_tensor_choose_kernel(up->cdim, edge1, dir, n) (&e_i[i*up->basis.num_basis]);
      ej_surf[n][i] = bc_block_tensor_choose_kernel(up->cdim, edge2, dir, n) (&ej[i*up->basis.num_basis]);
    }
  }

  // Now take the dot prduct at the quadrature nodes and fill the tensor
  // Only Need T11,13,31,33 in 2d
  int jctr=0;
  for (int j = 0; j < 3; j++){
    if(up->cdim==2 &&  j==1)
      continue;
    if(jctr!=dir){ // We only want to fill elements needed at this interface.
                // For example if we are at a z edge then we only need T^3'_1 and T^3'_3
                // At the corner cell, T^1'_1 and T^1'_3 will be filled using another blocks tan vecs
      jctr+=1;
      continue;
    }
    for(int n = 0; n < up->num_surf_nodes; n++) {
      int ictr=0;
      for (int i = 0; i < 3; i++){
        if(up->cdim==2 &&  i==1)
          continue;
        tj_i[up->cdim*up->num_surf_nodes*jctr + up->cdim*n + ictr] = dot_product(&ej_surf[n][3*j], &e_i_surf[n][3*i]);
        //printf("\n\nj,i = %d, %d\n", j,i);
        //printf("t = %g\n", tj_i[up->cdim*up->num_surf_nodes*jctr + up->cdim*n + ictr]);
        //printf("index = %d\n", up->cdim*up->num_surf_nodes*jctr + up->cdim*n + ictr);
        //printf("jctr,ictr,n = %d, %d, %d\n\n", jctr,ictr,n);
        ictr+=1;
      }
    }
    jctr+=1;
  }

}

void gkyl_bc_block_tensor_advance(struct bc_block_tensor* up, int dir, int edge1, int edge2,
    struct gkyl_array* dxdz1, struct gkyl_array* dzdx2, struct gkyl_range *range1, struct gkyl_range *range2)
{
  // Need to loop along only the directions != dir
  // For block 1, the index in dir will be min/max based on edge1
  // For block 2, index in dir is based on edge2

  int idx1[GKYL_MAX_DIM] = { 0};
  int idx2[GKYL_MAX_DIM] = { 0};
  idx1[dir] = edge1 == 0 ? range1->lower[dir] : range1->upper[dir];
  idx2[dir] = edge2 == 0 ? range2->lower[dir] : range2->upper[dir];
  
  struct gkyl_range range_def;
  int remdir[GKYL_MAX_DIM] = {0};
  remdir[dir] = 1;
  int locdir[GKYL_MAX_DIM] = {0};
  locdir[dir] = idx2[dir];
  gkyl_range_deflate(&range_def, range2, remdir, locdir);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range_def);
  while (gkyl_range_iter_next(&iter)) {
    // Fill the indices and fetch loc in each block
    int ictr = 0;
    for (int i = 0; i < up->cdim; i++) {
      if (i!=dir) {
        idx1[i] = iter.idx[ictr];
        idx2[i] = iter.idx[ictr];
        ictr+=1;
      }
    }
    long loc1 = gkyl_range_idx(range1, idx1);
    long loc2 = gkyl_range_idx(range1, idx2);

    const double* ej = gkyl_array_cfetch(dzdx2, loc2);
    const double* e_i = gkyl_array_cfetch(dxdz1, loc1);
    double* tj_i = gkyl_array_fetch(up->tensor, loc2);
    calc_tensor(up, dir, edge1, edge2, ej, e_i, tj_i);
  }
}
void gkyl_bc_block_tensor_release(struct bc_block_tensor* up)
{
  gkyl_array_release(up->tensor);
  gkyl_free(up);
}
