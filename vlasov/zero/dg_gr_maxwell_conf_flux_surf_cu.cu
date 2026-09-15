/* -*- c++ -*- */

#include <math.h>
#include <time.h>

extern "C" {
#include <assert.h>
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_gr_maxwell_conf_flux_surf.h>
#include <gkyl_dg_gr_maxwell_conf_flux_surf_priv.h>
#include <gkyl_dg_gr_maxwell_surf_and_vol_nodes.h> 
#include <gkyl_util.h> 
}

#include <cassert>

__global__ void
gkyl_dg_gr_maxwell_conf_flux_surf_advance_cu_kernel(struct gkyl_dg_gr_maxwell_conf_flux_surf *up, 
  const struct gkyl_range conf_range, const struct gkyl_range conf_range_ext, 
  const struct gkyl_array *jacob_pos, 
  const struct gkyl_surf_and_vol_node_arrays *lapse, const struct gkyl_surf_and_vol_node_arrays *shift, 
  const struct gkyl_surf_and_vol_node_arrays *h_ij, const struct gkyl_surf_and_vol_node_arrays *h_ij_inv, const struct gkyl_surf_and_vol_node_arrays *det_h, 
  const struct gkyl_array *field_con, const struct gkyl_array *field_no_J_con, 
  struct gkyl_array *cflrate, struct gkyl_array *conf_flux_surf)
{

  int cdim = up->cdim;
  int idx[GKYL_MAX_DIM], idx_l[GKYL_MAX_DIM], idx_r[GKYL_MAX_DIM]; 

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);  
    long cidx = gkyl_range_idx(&conf_range, idx);

    // Grab the cell center location for NC bracket calculation 
    double xcC[GKYL_MAX_DIM], xcR[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(&up->conf_grid, idx, xcC);

    const double *field_no_J_con_c = (const double*) gkyl_array_cfetch(field_no_J_con, cidx); 
    const double *field_con_c = (const double*) gkyl_array_cfetch(field_con, cidx); 
    double *cflrate_d = (double*) gkyl_array_fetch(cflrate, cidx);
    double *flux = (double*) gkyl_array_fetch(conf_flux_surf, cidx); 

    // Each cell owns *lower* fluxes in each configuration-space direction. 
    // So we need the distribution function in our current cell, and the cell
    // one index lower in each direction. If we are at the lower configuration-space
    // edge, we call ghost cells
    for (int dir = 0; dir<cdim; ++dir) {

      // Select the geometry elements based on direction
      const double *lapse_d = 0;
      const double *shift_d = 0;
      const double *h_ij_d = 0;
      const double *h_ij_inv_d = 0;
      const double *det_h_d = 0;

      if (dir == 0) {
        lapse_d = (const double*) gkyl_array_cfetch(lapse->nodal_arr_surf_x, cidx);
        shift_d = (const double*) gkyl_array_cfetch(shift->nodal_arr_surf_x, cidx);
        h_ij_d = (const double*) gkyl_array_cfetch(h_ij->nodal_arr_surf_x, cidx);
        h_ij_inv_d = (const double*) gkyl_array_cfetch(h_ij_inv->nodal_arr_surf_x, cidx);
        det_h_d = (const double*) gkyl_array_cfetch(det_h->nodal_arr_surf_x, cidx);
      } else if (dir == 1) {
        lapse_d = (const double*) gkyl_array_cfetch(lapse->nodal_arr_surf_y, cidx);
        shift_d = (const double*) gkyl_array_cfetch(shift->nodal_arr_surf_y, cidx);
        h_ij_d = (const double*) gkyl_array_cfetch(h_ij->nodal_arr_surf_y, cidx);
        h_ij_inv_d = (const double*) gkyl_array_cfetch(h_ij_inv->nodal_arr_surf_y, cidx);
        det_h_d = (const double*) gkyl_array_cfetch(det_h->nodal_arr_surf_y, cidx);
      } else {
        lapse_d = (const double*) gkyl_array_cfetch(lapse->nodal_arr_surf_z, cidx);
        shift_d = (const double*) gkyl_array_cfetch(shift->nodal_arr_surf_z, cidx);
        h_ij_d = (const double*) gkyl_array_cfetch(h_ij->nodal_arr_surf_z, cidx);
        h_ij_inv_d = (const double*) gkyl_array_cfetch(h_ij_inv->nodal_arr_surf_z, cidx);
        det_h_d = (const double*) gkyl_array_cfetch(det_h->nodal_arr_surf_z, cidx);
      }

      // Create an index for the left cell (which may be a ghost cell) 
      gkyl_copy_int_arr(cdim, idx, idx_l);
      idx_l[dir] = idx_l[dir]-1;
      long cidx_l = gkyl_range_idx(&conf_range, idx_l); 
      const double *field_no_J_con_l = (const double*) gkyl_array_cfetch(field_no_J_con, cidx_l);
      const double *field_con_l = (const double*) gkyl_array_cfetch(field_con, cidx_l);
      const double *jacob_pos_l = (const double*) gkyl_array_cfetch(jacob_pos, cidx_l);
      const double *jacob_pos_c = (const double*) gkyl_array_cfetch(jacob_pos, cidx);

      // For Points not along the domain-edge in theta, compute the left hand surface 
      // conf-flux.
      int theta_pole = 0;

      // If at the left edge in the theta dir
      if(idx[dir] == conf_range.lower[dir] && up->theta_pole_lo[dir]) {
        theta_pole = 1;
      }

      cflrate_d[0] += up->conf_flux_surf(up, dir, xcC, up->conf_grid.dx, theta_pole,
        jacob_pos_l, jacob_pos_c,
        lapse_d, shift_d, h_ij_d, h_ij_inv_d, det_h_d, field_con_l, field_con_c,
         field_no_J_con_l, field_no_J_con_c, flux);     

      // If at the right boundary compute flux owned by the point in the ghost cell
      if (idx[dir] == conf_range.upper[dir]) {

        // Index the right cell (ghost cell)
        gkyl_copy_int_arr(cdim, idx, idx_r);
        idx_r[dir] = idx_r[dir]+1;
        long cidx_r = gkyl_range_idx(&conf_range_ext, idx_r); 

        const double *field_no_J_con_r = (const double*) gkyl_array_cfetch(field_no_J_con, cidx_r);
        const double *field_con_r = (const double*) gkyl_array_cfetch(field_con, cidx_r);
        double *flux_r = (double*) gkyl_array_fetch(conf_flux_surf, cidx_r); 
        double *cflrate_d_r = (double*) gkyl_array_fetch(cflrate, cidx_r);

        // Grab the geometry information at the right interface surface
        if (dir == 0) {
          lapse_d = (const double*) gkyl_array_cfetch(lapse->nodal_arr_surf_x, cidx_r);
          shift_d = (const double*) gkyl_array_cfetch(shift->nodal_arr_surf_x, cidx_r);
          h_ij_d = (const double*) gkyl_array_cfetch(h_ij->nodal_arr_surf_x, cidx_r);
          h_ij_inv_d = (const double*) gkyl_array_cfetch(h_ij_inv->nodal_arr_surf_x, cidx_r);
          det_h_d = (const double*) gkyl_array_cfetch(det_h->nodal_arr_surf_x, cidx_r);
        } else if (dir == 1) {
          lapse_d = (const double*) gkyl_array_cfetch(lapse->nodal_arr_surf_y, cidx_r);
          shift_d = (const double*) gkyl_array_cfetch(shift->nodal_arr_surf_y, cidx_r);
          h_ij_d = (const double*) gkyl_array_cfetch(h_ij->nodal_arr_surf_y, cidx_r);
          h_ij_inv_d = (const double*) gkyl_array_cfetch(h_ij_inv->nodal_arr_surf_y, cidx_r);
          det_h_d = (const double*) gkyl_array_cfetch(det_h->nodal_arr_surf_y, cidx_r);
        } else {
          lapse_d = (const double*) gkyl_array_cfetch(lapse->nodal_arr_surf_z, cidx_r);
          shift_d = (const double*) gkyl_array_cfetch(shift->nodal_arr_surf_z, cidx_r);
          h_ij_d = (const double*) gkyl_array_cfetch(h_ij->nodal_arr_surf_z, cidx_r);
          h_ij_inv_d = (const double*) gkyl_array_cfetch(h_ij_inv->nodal_arr_surf_z, cidx_r);
          det_h_d = (const double*) gkyl_array_cfetch(det_h->nodal_arr_surf_z, cidx_r);
        }

        // If at the upper surface (in the ghost cell, computing the left edge)
        if(idx[dir] == conf_range.upper[dir] && up->theta_pole_up[dir]) {
          theta_pole = 1;
        }

        gkyl_rect_grid_cell_center(&up->conf_grid, idx_r, xcR);
        const double *jacob_pos_r = (const double*) gkyl_array_cfetch(jacob_pos, cidx_r);
        cflrate_d_r[0] += up->conf_flux_surf(up, dir, xcR, up->conf_grid.dx, theta_pole,
        jacob_pos_c, jacob_pos_r,
        lapse_d, shift_d, h_ij_d, h_ij_inv_d, det_h_d, field_con_c, field_con_r,
        field_no_J_con_c, field_no_J_con_r, flux_r);
      }
    }
  }
}

void 
gkyl_dg_gr_maxwell_conf_flux_surf_advance_cu(struct gkyl_dg_gr_maxwell_conf_flux_surf *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *conf_range_ext, 
  const struct gkyl_surf_and_vol_node_arrays *lapse, const struct gkyl_surf_and_vol_node_arrays *shift, 
  const struct gkyl_surf_and_vol_node_arrays *h_ij, const struct gkyl_surf_and_vol_node_arrays *h_ij_inv, const struct gkyl_surf_and_vol_node_arrays *det_h, 
  const struct gkyl_array *field_con, const struct gkyl_array *field_no_J_con, 
  struct gkyl_array *cflrate, struct gkyl_array *conf_flux_surf)
{
  int nblocks = conf_range->nblocks;
  int nthreads = conf_range->nthreads;
  gkyl_dg_gr_maxwell_conf_flux_surf_advance_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, 
    *conf_range, *conf_range_ext, up->jacob_pos->on_dev, lapse->on_dev, shift->on_dev, h_ij->on_dev, h_ij_inv->on_dev, det_h->on_dev, 
    field_con->on_dev, field_no_J_con->on_dev, cflrate->on_dev, conf_flux_surf->on_dev);  
}


// CUDA kernel to set device pointers to kernel functions
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
gkyl_dg_gr_maxwell_conf_flux_surf_set_cu_dev_ptrs(struct gkyl_dg_gr_maxwell_conf_flux_surf *up, 
    enum gkyl_basis_type b_type, int cdim, int poly_order, enum gkyl_field_id field_id)
{

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      // Kernels to compute the lax flux using the maximum eigenvalue and the left and right fluxes.
      up->lax_flux[0] = ser_lax_flux_x_kernels[cdim-1].kernels[poly_order];
      up->lax_flux[1] = ser_lax_flux_y_kernels[cdim-1].kernels[poly_order];
      up->lax_flux[2] = ser_lax_flux_z_kernels[cdim-1].kernels[poly_order];

      // Kernels to compute the maximum of the eigenvalues and isolate the fluxes E^i, H^i. 
      up->dg_gr_maxwell_alpha_quad[0] = ser_dg_gr_maxwell_alpha_quad_x_kernels[cdim-1].kernels[poly_order];
      up->dg_gr_maxwell_alpha_quad[1] = ser_dg_gr_maxwell_alpha_quad_y_kernels[cdim-1].kernels[poly_order];
      up->dg_gr_maxwell_alpha_quad[2] = ser_dg_gr_maxwell_alpha_quad_z_kernels[cdim-1].kernels[poly_order];

      break;

    case GKYL_BASIS_MODAL_TENSOR:
    // Kernels to compute the lax flux using the maximum eigenvalue and the left and right fluxes.
      up->lax_flux[0] = ten_lax_flux_x_kernels[cdim-1].kernels[poly_order];
      up->lax_flux[1] = ten_lax_flux_y_kernels[cdim-1].kernels[poly_order];
      up->lax_flux[2] = ten_lax_flux_z_kernels[cdim-1].kernels[poly_order];

      // Kernels to compute the maximum of the eigenvalues and isolate the fluxes E^i, H^i. 
      up->dg_gr_maxwell_alpha_quad[0] = ten_dg_gr_maxwell_alpha_quad_x_kernels[cdim-1].kernels[poly_order];
      up->dg_gr_maxwell_alpha_quad[1] = ten_dg_gr_maxwell_alpha_quad_y_kernels[cdim-1].kernels[poly_order];
      up->dg_gr_maxwell_alpha_quad[2] = ten_dg_gr_maxwell_alpha_quad_z_kernels[cdim-1].kernels[poly_order];

      break;      

    default:
      assert(false);
      break;    
  } 
  // Set assembly functions for computing fluxes. 
  up->conf_flux_surf = conf_flux_surf_kernels[cdim-1].kernels[poly_order]; 
}


gkyl_dg_gr_maxwell_conf_flux_surf*
gkyl_dg_gr_maxwell_conf_flux_surf_cu_dev_inew(const struct gkyl_dg_gr_maxwell_conf_flux_surf_inp *inp)
{
  struct gkyl_dg_gr_maxwell_conf_flux_surf *up = (struct gkyl_dg_gr_maxwell_conf_flux_surf*) gkyl_malloc(sizeof(*up));

  int cdim = inp->conf_basis->ndim;
  int poly_order = inp->conf_basis->poly_order;

  up->conf_grid = *inp->conf_grid;
  up->cdim = cdim;
  // Position map: acquired for lifetime safety; device pointer unpacked at
  // launch from the borrowed host-side array object.
  assert(inp->pos_map);
  up->pos_map = gkyl_vlasov_position_map_acquire(inp->pos_map);
  up->jacob_pos = inp->pos_map->jacob_pos;
  up->gr_maxwell_data.chi = inp->chi;
  up->gr_maxwell_data.gamma = inp->gamma;
  up->gr_maxwell_data.K_phi = 0.0;
  up->gr_maxwell_data.K_psi = 0.0;
  up->use_gpu = true;

  // Set the theta direction pole bc flags
  for (int i=0; i<cdim; ++i) {
    up->theta_pole_lo[i] = inp->theta_pole_lo[i];
    up->theta_pole_up[i] = inp->theta_pole_up[i];
  }

  up->flags = 0;

  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_dg_gr_maxwell_conf_flux_surf *up_cu = (struct gkyl_dg_gr_maxwell_conf_flux_surf*) gkyl_cu_malloc(sizeof(*up_cu));
    gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_dg_gr_maxwell_conf_flux_surf), GKYL_CU_MEMCPY_H2D);

  gkyl_dg_gr_maxwell_conf_flux_surf_set_cu_dev_ptrs<<<1,1>>>(up_cu, inp->conf_basis->b_type, 
    cdim, poly_order, inp->field_id);  

  // set parent on_dev pointer
  up->on_dev = up_cu;
 
  return up;

}
