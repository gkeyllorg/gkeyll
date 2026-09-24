#include <gkyl_nodal_ops.h>
#include <assert.h>

struct gkyl_nodal_ops*
gkyl_nodal_ops_new(const struct gkyl_basis *cbasis, const struct gkyl_rect_grid *grid, bool use_gpu)
{
  // Allocate space for new updater.
  struct gkyl_nodal_ops *up = gkyl_malloc(sizeof(*up));

  up->poly_order = cbasis->poly_order;

  struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, grid->ndim, cbasis->num_basis);
  cbasis->node_list(gkyl_array_fetch(nodes, 0));

  if (use_gpu) 
    up->nodes = gkyl_array_cu_dev_new(GKYL_DOUBLE, grid->ndim, cbasis->num_basis);
  else 
    up->nodes = gkyl_array_new(GKYL_DOUBLE, grid->ndim, cbasis->num_basis);

  // Copy the nodal values to the pre-allocated array
  gkyl_array_copy(up->nodes, nodes);
  gkyl_array_release(nodes);

  return up;
}

void 
gkyl_nodal_ops_n2m_corners(const struct gkyl_nodal_ops *nodal_ops, 
  const struct gkyl_basis *cbasis, const struct gkyl_rect_grid *grid, 
  const struct gkyl_range *nrange, const struct gkyl_range *update_range, int num_comp, 
  const struct gkyl_array *nodal_fld, struct gkyl_array *modal_fld) 
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(modal_fld)) {
    return gkyl_nodal_ops_n2m_cu(nodal_ops, cbasis, grid, 
      nrange, update_range, num_comp, 
      nodal_fld, modal_fld);
  }
#endif  
  double xc[GKYL_MAX_DIM];
  int num_basis = cbasis->num_basis;
  int cpoly_order = cbasis->poly_order;
  double fnodal[num_basis]; // to store nodal function values

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);
  int nidx[3];
  long lin_nidx[num_basis];
  
  while (gkyl_range_iter_next(&iter)) {
     gkyl_rect_grid_cell_center(grid, iter.idx, xc);

    for (int i=0; i<num_basis; ++i) {
      const double *temp  = gkyl_array_cfetch(nodal_ops->nodes,i);
      for( int j = 0; j < grid->ndim; j++) {
        if (cpoly_order==1) {
          nidx[j] = (iter.idx[j]-update_range->lower[j]) + (temp[j]+1)/2 ;
        }
        if (cpoly_order==2) {
          nidx[j] = 2*(iter.idx[j]-update_range->lower[j]) + (temp[j]+1) ;
        }
      }
      lin_nidx[i] = gkyl_range_idx(nrange, nidx);
    }

    long lidx = gkyl_range_idx(update_range, iter.idx);
    double *arr_p = gkyl_array_fetch(modal_fld, lidx);
    double fao[num_basis*num_comp];
  
    for (int i=0; i<num_basis; ++i) {
      const double *temp = gkyl_array_cfetch(nodal_fld, lin_nidx[i]);
      for (int j=0; j<num_comp; ++j) {
        fao[i*num_comp + j] = temp[j];
      }
    }

    for (int i=0; i<num_comp; ++i) {
      // copy so nodal values for each return value are contiguous
      // (recall that function can have more than one return value)
      for (int k=0; k<num_basis; ++k)
        fnodal[k] = fao[num_comp*k+i];
      // transform to modal expansion
      cbasis->nodal_to_modal(fnodal, &arr_p[num_basis*i]);
    }
  }
}

void 
gkyl_nodal_ops_n2m_interior(const struct gkyl_nodal_ops *nodal_ops, 
  const struct gkyl_basis *cbasis, const struct gkyl_rect_grid *grid, 
  const struct gkyl_range *nrange, const struct gkyl_range *update_range, int num_comp, 
  const struct gkyl_array *nodal_fld, struct gkyl_array *modal_fld) 
{

  int num_basis = cbasis->num_basis;
  int cpoly_order = cbasis->poly_order;
  double fnodal[num_basis]; // to store nodal function values

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);
  int nidx[3];
  long lin_nidx[num_basis];
  
  while (gkyl_range_iter_next(&iter)) {

    for (int i=0; i<num_basis; ++i) {
      const double *temp  = gkyl_array_cfetch(nodal_ops->nodes,i);
      for( int j = 0; j < grid->ndim; j++) {
        int mod = (j==1 && grid->ndim==3) ? 2 * (i/4) : 0;
        if (j < grid->ndim-1)
          nidx[j] = (iter.idx[j]-update_range->lower[j])*2 + i/((int) pow(2,grid->ndim-1-j)) - mod;
        else
          nidx[j] = (iter.idx[j]-update_range->lower[j])*2 + i%2 ;
      }
      lin_nidx[i] = gkyl_range_idx(nrange, nidx);
    }

    long lidx = gkyl_range_idx(update_range, iter.idx);
    double *arr_p = gkyl_array_fetch(modal_fld, lidx);
    double fao[num_basis*num_comp];
  
    for (int i=0; i<num_basis; ++i) {
      const double *temp = gkyl_array_cfetch(nodal_fld, lin_nidx[i]);
      for (int j=0; j<num_comp; ++j) {
        fao[i*num_comp + j] = temp[j];
      }
    }

    for (int i=0; i<num_comp; ++i) {
      // copy so nodal values for each return value are contiguous
      // (recall that function can have more than one return value)
      for (int k=0; k<num_basis; ++k) {
        fnodal[k] = fao[num_comp*k+i];
      }
      // transform to modal expansion
      for (int k=0; k<num_basis; ++k)
        cbasis->quad_nodal_to_modal(fnodal, &arr_p[num_basis*i],k);
    }
  }
}

void 
gkyl_nodal_ops_n2m_surface_1x(const struct gkyl_nodal_ops *nodal_ops, 
  const struct gkyl_rect_grid *grid, 
  const struct gkyl_range *nrange, const struct gkyl_range *update_range, int num_comp, 
  const struct gkyl_array *nodal_fld, struct gkyl_array *modal_fld, int dir) 
{

  int num_basis = 1;
  double fnodal[num_basis]; // to store nodal function values

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);
  int nidx[3];
  long lin_nidx[num_basis];
  
  while (gkyl_range_iter_next(&iter)) {
    for (int i=0; i<num_basis; ++i) {
      const double *temp  = gkyl_array_cfetch(nodal_ops->nodes,i);
      for( int j = 0; j < grid->ndim; j++) {
        nidx[j] = (iter.idx[j]-update_range->lower[j]) ;
      }
      lin_nidx[i] = gkyl_range_idx(nrange, nidx);
    }

    long lidx = gkyl_range_idx(update_range, iter.idx);
    double *arr_p = gkyl_array_fetch(modal_fld, lidx);
    double fao[num_basis*num_comp];
  
    for (int i=0; i<num_basis; ++i) {
      const double *temp = gkyl_array_cfetch(nodal_fld, lin_nidx[i]);
      for (int j=0; j<num_comp; ++j) {
        fao[i*num_comp + j] = temp[j];
      }
    }

    for (int i=0; i<num_comp; ++i) {
      // copy so nodal values for each return value are contiguous
      // (recall that function can have more than one return value)
      for (int k=0; k<num_basis; ++k) {
        fnodal[k] = fao[num_comp*k+i];
      }
      // transform to modal expansion
      for (int k=0; k<num_basis; ++k)
        arr_p[num_basis*i] = fnodal[k];
    }
  }
}

void 
gkyl_nodal_ops_n2m_surface(const struct gkyl_nodal_ops *nodal_ops, 
  const struct gkyl_basis *cbasis, const struct gkyl_rect_grid *grid, 
  const struct gkyl_range *nrange, const struct gkyl_range *update_range, int num_comp, 
  const struct gkyl_array *nodal_fld, struct gkyl_array *modal_fld, int dir) 
{
  if(grid->ndim==1)
    return gkyl_nodal_ops_n2m_surface_1x(nodal_ops, grid, nrange, update_range, num_comp, nodal_fld, modal_fld, dir);

  int num_basis = cbasis->num_basis;
  int cpoly_order = cbasis->poly_order;
  double fnodal[num_basis]; // to store nodal function values

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);
  int nidx[3];
  long lin_nidx[num_basis];
  
  while (gkyl_range_iter_next(&iter)) {
    for (int i=0; i<num_basis; ++i) {
      const double *temp  = gkyl_array_cfetch(nodal_ops->nodes,i);
      for( int j = 0; j < grid->ndim; j++) {
          if (j !=dir) {
            if (j < grid->ndim-1) {
              if (dir == 0)
                nidx[j] = (iter.idx[j]-update_range->lower[j])*2 + i/((int) pow(2,grid->ndim-1-j));
              else if (dir == 1)
                nidx[j] = (iter.idx[j]-update_range->lower[j])*2 + i/((int) pow(2,grid->ndim-2+j));
              else if (dir == 2) {
                if(j==0)
                  nidx[j] = (iter.idx[j]-update_range->lower[j])*2 + i/((int) pow(2,grid->ndim-2+j));
                if(j==1)
                  nidx[j] = (iter.idx[j]-update_range->lower[j])*2 + i%2;
              }
            }
            else
              nidx[j] = (iter.idx[j]-update_range->lower[j])*2 + i%2 ;
          }
          else {
            nidx[j] = (iter.idx[j]-update_range->lower[j]) ;
          }
      }
      lin_nidx[i] = gkyl_range_idx(nrange, nidx);

    }

    long lidx = gkyl_range_idx(update_range, iter.idx);
    double *arr_p = gkyl_array_fetch(modal_fld, lidx);
    double fao[num_basis*num_comp];
  
    for (int i=0; i<num_basis; ++i) {
      const double *temp = gkyl_array_cfetch(nodal_fld, lin_nidx[i]);
      for (int j=0; j<num_comp; ++j) {
        fao[i*num_comp + j] = temp[j];
      }
    }

    for (int i=0; i<num_comp; ++i) {
      // copy so nodal values for each return value are contiguous
      // (recall that function can have more than one return value)
      for (int k=0; k<num_basis; ++k) {
        fnodal[k] = fao[num_comp*k+i];
      }
      // transform to modal expansion
      for (int k=0; k<num_basis; ++k)
        cbasis->quad_nodal_to_modal(fnodal, &arr_p[num_basis*i],k);
    }
  }
}

void 
gkyl_nodal_ops_m2n_surface_1x(const struct gkyl_nodal_ops *nodal_ops, 
  const struct gkyl_rect_grid *grid, 
  const struct gkyl_range *nrange, const struct gkyl_range *update_range, int num_comp, 
  struct gkyl_array *nodal_fld, const struct gkyl_array *modal_fld, int dir) 
{

  int num_basis = 1;
  double fnodal[num_basis]; // to store nodal function values

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);
  int nidx[3];
  long lin_nidx[num_basis];
  
  while (gkyl_range_iter_next(&iter)) {
    for (int i=0; i<num_basis; ++i) {
      const double *temp  = gkyl_array_cfetch(nodal_ops->nodes,i);
      for( int j = 0; j < grid->ndim; j++) {
        nidx[j] = (iter.idx[j]-update_range->lower[j]) ;
      }
      lin_nidx[i] = gkyl_range_idx(nrange, nidx);
    }

    long lidx = gkyl_range_idx(update_range, iter.idx);
    const double *arr_p = gkyl_array_cfetch(modal_fld, lidx);
    double fao[num_basis*num_comp];

    for (int i=0; i<num_comp; ++i) {
      // transform to modal expansion
      for (int k=0; k<num_basis; ++k){
        fnodal[k] = arr_p[num_basis*i];
      }
      // copy so nodal values for each return value are contiguous
      // (recall that function can have more than one return value)
      for (int k=0; k<num_basis; ++k)
        fao[num_comp*k+i] = fnodal[k];
    }

    for (int i=0; i<num_basis; ++i) {
      double *temp = gkyl_array_fetch(nodal_fld, lin_nidx[i]);
      for (int j=0; j<num_comp; ++j) {
         temp[j] = fao[i*num_comp + j];
      }
    }
   }
}

void 
gkyl_nodal_ops_m2n_surface(const struct gkyl_nodal_ops *nodal_ops, 
  const struct gkyl_basis *cbasis, const struct gkyl_rect_grid *grid, 
  const struct gkyl_range *nrange, const struct gkyl_range *update_range, int num_comp, 
  struct gkyl_array *nodal_fld, const struct gkyl_array *modal_fld, int dir) 
{
  if(grid->ndim==1)
    return gkyl_nodal_ops_m2n_surface_1x(nodal_ops, grid, nrange, update_range, num_comp, nodal_fld, modal_fld, dir);

  int num_basis = cbasis->num_basis;
  int cpoly_order = cbasis->poly_order;
  double fnodal[num_basis]; // to store nodal function values

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);
  int nidx[3];
  long lin_nidx[num_basis];
  
  while (gkyl_range_iter_next(&iter)) {
    for (int i=0; i<num_basis; ++i) {
      const double *temp  = gkyl_array_cfetch(nodal_ops->nodes,i);
      for( int j = 0; j < grid->ndim; j++) {
          if (j !=dir) {
            if (j < grid->ndim-1) {
              if (dir == 0)
                nidx[j] = (iter.idx[j]-update_range->lower[j])*2 + i/((int) pow(2,grid->ndim-1-j));
              else if (dir == 1)
                nidx[j] = (iter.idx[j]-update_range->lower[j])*2 + i/((int) pow(2,grid->ndim-2+j));
              else if (dir == 2) {
                if(j==0)
                  nidx[j] = (iter.idx[j]-update_range->lower[j])*2 + i/((int) pow(2,grid->ndim-2+j));
                if(j==1)
                  nidx[j] = (iter.idx[j]-update_range->lower[j])*2 + i%2;
              }
            }
            else
              nidx[j] = (iter.idx[j]-update_range->lower[j])*2 + i%2 ;
          }
          else {
            nidx[j] = (iter.idx[j]-update_range->lower[j]) ;
          }
      }
      lin_nidx[i] = gkyl_range_idx(nrange, nidx);

    }

    long lidx = gkyl_range_idx(update_range, iter.idx);
    const double *arr_p = gkyl_array_cfetch(modal_fld, lidx);
    double fao[num_basis*num_comp];

    for (int i=0; i<num_comp; ++i) {
      // transform to modal expansion
      for (int k=0; k<cbasis->num_quad; ++k){
        cbasis->modal_to_quad_nodal(&arr_p[num_basis*i], fnodal, k);
      }
      // copy so nodal values for each return value are contiguous
      // (recall that function can have more than one return value)
      for (int k=0; k<num_basis; ++k)
        fao[num_comp*k+i] = fnodal[k];
    }

    for (int i=0; i<num_basis; ++i) {
      double *temp = gkyl_array_fetch(nodal_fld, lin_nidx[i]);
      for (int j=0; j<num_comp; ++j) {
         temp[j] = fao[i*num_comp + j];
      }
    }
   }
}



void 
gkyl_nodal_ops_m2n_corners_p2(const struct gkyl_nodal_ops *nodal_ops, 
  const struct gkyl_basis *cbasis, const struct gkyl_rect_grid *grid, 
  const struct gkyl_range *nrange, const struct gkyl_range *update_range, int num_comp, 
  struct gkyl_array *nodal_fld, const struct gkyl_array *modal_fld) 
{
  int num_basis = cbasis->num_basis;
  int cpoly_order = cbasis->poly_order;
  struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, grid->ndim, cbasis->num_basis);
  cbasis->node_list(gkyl_array_fetch(nodes, 0));
  double fnodal[num_basis]; // to store nodal function values
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);
  int nidx[3];
  long lin_nidx[num_basis];
  
  while (gkyl_range_iter_next(&iter)) {
    for (int i=0; i<num_basis; ++i) {
      const double* temp  = gkyl_array_cfetch(nodes,i);
      for( int j = 0; j < grid->ndim; j++){
        nidx[j] = 2*(iter.idx[j]-update_range->lower[j]) + (temp[j] + 1) ;
      }
      lin_nidx[i] = gkyl_range_idx(nrange, nidx);
    }
    long lidx = gkyl_range_idx(update_range, iter.idx);
    const double *arr_p = gkyl_array_cfetch(modal_fld, lidx);
    double fao[num_basis*num_comp];
    // already fetched modal coeffs in arr_p
    // Now we are going to need to fill the nodal values
    // So in the loop of num_nodes/num_basis we will fetch the nodal value at lin_nidx[i]
    // at each place do an eval_basis at logical coords
  
    for (int i=0; i<num_basis; ++i) {
      double* temp = gkyl_array_fetch(nodal_fld, lin_nidx[i]);
      const double*  node_i  = gkyl_array_cfetch(nodes,i);
      for (int j=0; j<num_comp; ++j) {
        temp[j] = cbasis->eval_expand(node_i, &arr_p[j*num_basis]);
      }
    }
  }

  gkyl_array_release(nodes);
}

void 
gkyl_nodal_ops_m2n_corners(const struct gkyl_nodal_ops *nodal_ops, 
  const struct gkyl_basis *cbasis, const struct gkyl_rect_grid *grid, 
  const struct gkyl_range *nrange, const struct gkyl_range *update_range, int num_comp, 
  struct gkyl_array *nodal_fld, const struct gkyl_array *modal_fld) 
{
    
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(nodal_fld)) {
    if (nodal_ops->poly_order==2)
      assert(false);
    return gkyl_nodal_ops_m2n_cu(nodal_ops, cbasis, grid, 
      nrange, update_range, num_comp, 
      nodal_fld, modal_fld);
  }
#endif 

  if(cbasis->poly_order == 2)
    return gkyl_nodal_ops_m2n_corners_p2(nodal_ops, cbasis, grid, nrange, update_range, num_comp, nodal_fld, modal_fld);

  int num_basis = cbasis->num_basis;
  struct gkyl_range_iter iter;
  // do the nodal loop instead
  gkyl_range_iter_init(&iter, nrange);
  int idx[3];
  long lin_nidx;
  while (gkyl_range_iter_next(&iter)) { // iter.idx = nidx
    lin_nidx = gkyl_range_idx(nrange, idx);
    int node_idx = 0;
    for( int j = 0; j < grid->ndim; j++){
      int mod = j==0 ? 1 : 0;
      if (iter.idx[j] == nrange->upper[j]) {
        idx[j] = iter.idx[j] + update_range->lower[j]-1;
        node_idx += 2*j + mod;
      }
      else {
        idx[j] = iter.idx[j] + update_range->lower[j];
      }
    }
    lin_nidx = gkyl_range_idx(nrange, iter.idx);
    long lidx = gkyl_range_idx(update_range, idx);
    const double *arr_p = gkyl_array_cfetch(modal_fld, lidx);
    double *temp = gkyl_array_fetch(nodal_fld, lin_nidx);


    const double *node_i  = gkyl_array_cfetch(nodal_ops->nodes, node_idx);
    for (int j=0; j<num_comp; ++j) {
      temp[j] = cbasis->eval_expand(node_i, &arr_p[j*num_basis]);
    }


  }
}

void 
gkyl_nodal_ops_m2n_interior(const struct gkyl_nodal_ops *nodal_ops, 
  const struct gkyl_basis *cbasis, const struct gkyl_rect_grid *grid, 
  const struct gkyl_range *nrange, const struct gkyl_range *update_range, int num_comp, 
  struct gkyl_array *nodal_fld, const struct gkyl_array *modal_fld) 
{

  int num_basis = cbasis->num_basis;
  int cpoly_order = cbasis->poly_order;
  double fnodal[num_basis]; // to store nodal function values

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);
  int nidx[3];
  long lin_nidx[num_basis];
  
  while (gkyl_range_iter_next(&iter)) {

    for (int i=0; i<num_basis; ++i) {
      const double *temp  = gkyl_array_cfetch(nodal_ops->nodes,i);
      for( int j = 0; j < grid->ndim; j++) {
        int mod = (j==1 && grid->ndim==3) ? 2 * (i/4) : 0;
        if (j < grid->ndim-1)
          nidx[j] = (iter.idx[j]-update_range->lower[j])*2 + i/((int) pow(2,grid->ndim-1-j)) - mod;
        else
          nidx[j] = (iter.idx[j]-update_range->lower[j])*2 + i%2 ;
      }
      lin_nidx[i] = gkyl_range_idx(nrange, nidx);
    }

    long lidx = gkyl_range_idx(update_range, iter.idx);
    const double *arr_p = gkyl_array_cfetch(modal_fld, lidx);
    double fao[num_basis*num_comp];
  


    for (int i=0; i<num_comp; ++i) {
      // transform to modal expansion
      for (int k=0; k<cbasis->num_quad; ++k){
        cbasis->modal_to_quad_nodal(&arr_p[num_basis*i], fnodal, k);
      }
      // copy so nodal values for each return value are contiguous
      // (recall that function can have more than one return value)
      for (int k=0; k<num_basis; ++k)
        fao[num_comp*k+i] = fnodal[k];
    }

    for (int i=0; i<num_basis; ++i) {
      double *temp = gkyl_array_fetch(nodal_fld, lin_nidx[i]);
      for (int j=0; j<num_comp; ++j) {
         temp[j] = fao[i*num_comp + j];
      }
    }

  }

}




void 
gkyl_nodal_ops_m2n_deflated(const struct gkyl_nodal_ops *nodal_ops, 
  const struct gkyl_basis *deflated_cbasis, const struct gkyl_rect_grid *deflated_grid, 
  const struct gkyl_range *nrange, const struct gkyl_range *deflated_nrange, const struct gkyl_range *deflated_update_range, int num_comp, 
  struct gkyl_array *nodal_fld, const struct gkyl_array *deflated_modal_fld, int extra_idx) 
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(nodal_fld)) {
    return gkyl_nodal_ops_m2n_deflated_cu(nodal_ops, deflated_cbasis, deflated_grid, 
      nrange, deflated_nrange, deflated_update_range, num_comp, 
      nodal_fld, deflated_modal_fld, extra_idx);
  }
#endif 
  int num_basis = deflated_cbasis->num_basis;
  int idx[3];
  int midx[3];
  idx[deflated_grid->ndim] = extra_idx;
  //while (gkyl_range_iter_next(&iter)) { // iter.idx = nidx
  for(int linc1 = 0; linc1 < deflated_nrange->volume; linc1++){
    gkyl_sub_range_inv_idx(deflated_nrange, linc1, idx);
    long linc = gkyl_range_idx(nrange, idx);
    int node_idx = 0;
    for( int j = 0; j < deflated_grid->ndim; j++){
      int mod = j==0 ? 1 : 0;
      if (idx[j] == deflated_nrange->upper[j]) {
        midx[j] = idx[j] + deflated_update_range->lower[j]-1;
        node_idx += 2*j + mod;
      }
      else {
        midx[j] = idx[j] + deflated_update_range->lower[j];
      }
    }
    long lidx = gkyl_range_idx(deflated_update_range, midx);
    const double *arr_p = gkyl_array_cfetch(deflated_modal_fld, lidx);
    double *temp = gkyl_array_fetch(nodal_fld, linc);
    const double *node_i  = gkyl_array_cfetch(nodal_ops->nodes, node_idx);
    for (int j=0; j<num_comp; ++j) {
      temp[j] = deflated_cbasis->eval_expand(node_i, &arr_p[j*num_basis]);
    }
  }
}


void 
gkyl_nodal_ops_n2m(const struct gkyl_nodal_ops *nodal_ops, 
  const struct gkyl_basis *cbasis, const struct gkyl_rect_grid *grid, 
  const struct gkyl_range *nrange, const struct gkyl_range *update_range, int num_comp, 
  const struct gkyl_array *nodal_fld, struct gkyl_array *modal_fld, bool use_quad) 
{
  if (use_quad)
    return gkyl_nodal_ops_n2m_interior(nodal_ops, cbasis, grid, nrange, update_range, num_comp, nodal_fld, modal_fld);
  else
    return gkyl_nodal_ops_n2m_corners(nodal_ops, cbasis, grid, nrange, update_range, num_comp, nodal_fld, modal_fld);
}

void 
gkyl_nodal_ops_m2n(const struct gkyl_nodal_ops *nodal_ops, 
  const struct gkyl_basis *cbasis, const struct gkyl_rect_grid *grid, 
  const struct gkyl_range *nrange, const struct gkyl_range *update_range, int num_comp, 
  struct gkyl_array *nodal_fld, const struct gkyl_array *modal_fld, bool use_quad) 
{
  if (use_quad)
    return gkyl_nodal_ops_m2n_interior(nodal_ops, cbasis, grid, nrange, update_range, num_comp, nodal_fld, modal_fld);
  else
    return gkyl_nodal_ops_m2n_corners(nodal_ops, cbasis, grid, nrange, update_range, num_comp, nodal_fld, modal_fld);
}


void 
gkyl_nodal_ops_release(struct gkyl_nodal_ops *up)
{
  gkyl_array_release(up->nodes);
  gkyl_free(up);
}
