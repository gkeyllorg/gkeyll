/* -*- c++ -*- */

#include "gkyl_util.h"
#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_nodal_ops.h>
}

__global__ static void
gkyl_nodal_ops_n2m_cu_kernel(const struct gkyl_basis *cbasis, 
  struct gkyl_rect_grid grid, struct gkyl_range nrange, struct gkyl_range update_range, 
  const struct gkyl_array *nodes, int num_comp, 
  const struct gkyl_array *nodal_fld, struct gkyl_array *modal_fld)
{
  double xc[GKYL_MAX_DIM];
  int idx[GKYL_MAX_DIM];
  int num_basis = cbasis->num_basis;
  int cpoly_order = cbasis->poly_order;
  double fnodal[20]; // to store nodal function values

  int nidx[3];
  long lin_nidx[20];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < update_range.volume; linc1 += blockDim.x*gridDim.x) {
    // inverse index from linc1 to idxc
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idxc={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&update_range, linc1, idx);
    gkyl_rect_grid_cell_center(&grid, idx, xc);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long linc = gkyl_range_idx(&update_range, idx);

    for (int i=0; i<num_basis; ++i) {
      const double *temp  = (const double *) gkyl_array_cfetch(nodes, i);
      for( int j = 0; j < grid.ndim; j++){
        if (cpoly_order==1) {
          nidx[j] = (idx[j]-update_range.lower[j]) + (temp[j]+1)/2 ;
        }
        if (cpoly_order==2) {
          nidx[j] = 2*(idx[j]-update_range.lower[j]) + (temp[j]+1) ;
        }
      }
      lin_nidx[i] = gkyl_range_idx(&nrange, nidx);
    }

    double *arr_p = (double *) gkyl_array_fetch(modal_fld, linc);
    double fao[20*9];
  
    for (int i=0; i<num_basis; ++i) {
      const double *temp = (const double *) gkyl_array_cfetch(nodal_fld, lin_nidx[i]);
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
gkyl_nodal_ops_n2m_cu(const struct gkyl_nodal_ops *nodal_ops, 
  const struct gkyl_basis *cbasis, const struct gkyl_rect_grid *grid, 
  const struct gkyl_range *nrange, const struct gkyl_range *update_range, int num_comp, 
  const struct gkyl_array *nodal_fld, struct gkyl_array *modal_fld) 
{
  int nblocks = update_range->nblocks;
  int nthreads = update_range->nthreads;

  gkyl_nodal_ops_n2m_cu_kernel<<<nblocks, nthreads>>>(cbasis, *grid, 
    *nrange, *update_range, nodal_ops->nodes->on_dev, num_comp, 
    nodal_fld->on_dev, modal_fld->on_dev);
}

__global__ static void
gkyl_nodal_ops_m2n_cu_kernel(const struct gkyl_basis *cbasis, 
  struct gkyl_rect_grid grid, struct gkyl_range nrange, struct gkyl_range update_range, 
  const struct gkyl_array *nodes, int num_comp, 
  struct gkyl_array *nodal_fld, const struct gkyl_array *modal_fld)
{
  int idx[GKYL_MAX_DIM];
  int midx[GKYL_MAX_DIM];
  int num_basis = cbasis->num_basis;

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < nrange.volume; linc1 += blockDim.x*gridDim.x) {
    // inverse index from linc1 to idxc
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idxc={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&nrange, linc1, idx);
    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long linc = gkyl_range_idx(&nrange, idx);
    int node_idx = 0;
    for( int j = 0; j < grid.ndim; j++){
      int mod = j==0 ? 1 : 0;
      if (idx[j] == nrange.upper[j]) {
        midx[j] = idx[j] + update_range.lower[j]-1;
        node_idx += 2*j + mod;
      }
      else {
        midx[j] = idx[j] + update_range.lower[j];
      }
    }
    long lidx = gkyl_range_idx(&update_range, midx);
    const double *arr_p = (const double *) gkyl_array_cfetch(modal_fld, lidx);
    double *temp = (double *) gkyl_array_fetch(nodal_fld, linc);
    const double *node_i  = (const double *) gkyl_array_cfetch(nodes, node_idx);
    for (int j=0; j<num_comp; ++j) {
      temp[j] = cbasis->eval_expand(node_i, &arr_p[j*num_basis]);
    }
  }
}

void 
gkyl_nodal_ops_m2n_cu(const struct gkyl_nodal_ops *nodal_ops, 
  const struct gkyl_basis *cbasis, const struct gkyl_rect_grid *grid, 
  const struct gkyl_range *nrange, const struct gkyl_range *update_range, int num_comp, 
  struct gkyl_array *nodal_fld, const struct gkyl_array *modal_fld) 
{
  int nblocks = update_range->nblocks;
  int nthreads = update_range->nthreads;

  gkyl_nodal_ops_m2n_cu_kernel<<<nblocks, nthreads>>>(cbasis, *grid, 
    *nrange, *update_range, nodal_ops->nodes->on_dev, num_comp, 
    nodal_fld->on_dev, modal_fld->on_dev);
}

__global__ static void
gkyl_nodal_ops_m2n_deflated_cu_kernel(const struct gkyl_basis *deflated_cbasis, 
  struct gkyl_rect_grid deflated_grid, struct gkyl_range nrange, struct gkyl_range def_nrange, struct gkyl_range deflated_update_range, 
  const struct gkyl_array *nodes, int num_comp, 
  struct gkyl_array *nodal_fld, const struct gkyl_array *deflated_modal_fld, int extra_idx)
{
  int idx[GKYL_MAX_DIM];
  int midx[GKYL_MAX_DIM];
  int num_basis = deflated_cbasis->num_basis;
  idx[deflated_grid.ndim] = extra_idx;

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
    linc1 < def_nrange.volume; linc1 += blockDim.x*gridDim.x) {
    // inverse index from linc1 to idxc
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idxc={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&def_nrange, linc1, idx);
    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long linc = gkyl_range_idx(&nrange, idx);
    int node_idx = 0;
    for( int j = 0; j < deflated_grid.ndim; j++){
      int mod = j==0 ? 1 : 0;
      if (idx[j] == def_nrange.upper[j]) {
        midx[j] = idx[j] + deflated_update_range.lower[j]-1;
        node_idx += 2*j + mod;
      }
      else {
        midx[j] = idx[j] + deflated_update_range.lower[j];
      }
    }
    long lidx = gkyl_range_idx(&deflated_update_range, midx);
    const double *arr_p = (const double *) gkyl_array_cfetch(deflated_modal_fld, lidx);
    double *temp = (double *) gkyl_array_fetch(nodal_fld, linc);
    const double *node_i  = (const double *) gkyl_array_cfetch(nodes, node_idx);
    for (int j=0; j<num_comp; ++j) {
      temp[j] = deflated_cbasis->eval_expand(node_i, &arr_p[j*num_basis]);
    }
  }
}

void 
gkyl_nodal_ops_m2n_deflated_cu(const struct gkyl_nodal_ops *nodal_ops, 
  const struct gkyl_basis *deflated_cbasis, const struct gkyl_rect_grid *deflated_grid, 
  const struct gkyl_range *nrange, const struct gkyl_range *deflated_nrange,  const struct gkyl_range *deflated_update_range, int num_comp, 
  struct gkyl_array *nodal_fld, const struct gkyl_array *deflated_modal_fld, int extra_idx) 
{
  int nblocks = deflated_update_range->nblocks;
  int nthreads = deflated_update_range->nthreads;

  gkyl_nodal_ops_m2n_deflated_cu_kernel<<<nblocks, nthreads>>>(deflated_cbasis, *deflated_grid, 
    *nrange, *deflated_nrange, *deflated_update_range, nodal_ops->nodes->on_dev, num_comp, 
    nodal_fld->on_dev, deflated_modal_fld->on_dev, extra_idx);
}
