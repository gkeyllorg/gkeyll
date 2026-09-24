#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_calc_bmag.h>
#include <gkyl_calc_derived_geo.h>
#include <gkyl_calc_metric.h>
#include <gkyl_comm.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_priv.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_math.h>
#include <gkyl_nodal_ops.h>
#include <gkyl_dg_bin_ops.h>
#include <assert.h>

static
void gk_geometry_mapc2p_advance(struct gk_geometry* up, struct gkyl_range *nrange, double dzc[3], 
  evalf_t mapc2p_func, void* mapc2p_ctx, evalf_t bfield_func, void *bfield_ctx, struct gkyl_position_map *position_map)
{
  //Now project mapc2p and the FD array
  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates
                                
  double dtheta = up->grid.dx[TH_IDX],
    dpsi = up->grid.dx[PSI_IDX],
    dalpha = up->grid.dx[AL_IDX];

  double theta_lo = up->grid.lower[TH_IDX] + (up->local.lower[TH_IDX] - up->global.lower[TH_IDX])*up->grid.dx[TH_IDX],
    psi_lo = up->grid.lower[PSI_IDX] + (up->local.lower[PSI_IDX] - up->global.lower[PSI_IDX])*up->grid.dx[PSI_IDX],
    alpha_lo = up->grid.lower[AL_IDX] + (up->local.lower[AL_IDX] - up->global.lower[AL_IDX])*up->grid.dx[AL_IDX];

  double dx_fact = up->basis.poly_order == 1 ? 1 : 0.5;
  dtheta *= dx_fact; dpsi *= dx_fact; dalpha *= dx_fact;
                                
  int cidx[3] = { 0 };
  for(int ia=nrange->lower[AL_IDX]; ia<=nrange->upper[AL_IDX]; ++ia){
    cidx[AL_IDX] = ia;
    double alpha_curr = alpha_lo + ia*dalpha;
    position_map->maps[1](0.0, &alpha_curr, &alpha_curr, position_map->ctxs[1]);
    for (int ip=nrange->lower[PSI_IDX]; ip<=nrange->upper[PSI_IDX]; ++ip) {
      double psi_curr = psi_lo + ip*dpsi;
      position_map->maps[0](0.0, &psi_curr,   &psi_curr,   position_map->ctxs[0]);
      cidx[PSI_IDX] = ip;
      // set node coordinates
      for (int it=nrange->lower[TH_IDX]; it<=nrange->upper[TH_IDX]; ++it) {
        double theta_curr = theta_lo + it*dtheta;
        position_map->maps[2](0.0, &theta_curr, &theta_curr, position_map->ctxs[2]);
        cidx[TH_IDX] = it;

        double *mc2p_n = (double *) gkyl_array_fetch(up->geo_corn.mc2p_nodal, gkyl_range_idx(nrange, cidx));
        double *bmag_n = (double *) gkyl_array_fetch(up->geo_corn.bmag_nodal, gkyl_range_idx(nrange, cidx));
        double xyz[3] = {psi_curr, alpha_curr, theta_curr};
        double XYZ[3] = {0.};
        double B[3] = {0.};
        mapc2p_func(0.0, xyz, XYZ, mapc2p_ctx);
        bfield_func(0.0, xyz, B, bfield_ctx);

        mc2p_n[X_IDX] = XYZ[X_IDX];
        mc2p_n[Y_IDX] = XYZ[Y_IDX];
        mc2p_n[Z_IDX] = XYZ[Z_IDX];
        bmag_n[0] = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);

        double *mc2nu_n = gkyl_array_fetch(up->geo_corn.mc2nu_pos_nodal, gkyl_range_idx(nrange, cidx));
        mc2nu_n[X_IDX] = psi_curr;
        mc2nu_n[Y_IDX] = alpha_curr;
        mc2nu_n[Z_IDX] = theta_curr;
      }
    }
  }

  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&up->basis, &up->grid, false);
  gkyl_nodal_ops_n2m(n2m, &up->basis, &up->grid, nrange, &up->local, 3, up->geo_corn.mc2p_nodal, up->geo_corn.mc2p, false);
  gkyl_nodal_ops_n2m(n2m, &up->basis, &up->grid, nrange, &up->local, 3, up->geo_corn.mc2nu_pos_nodal, up->geo_corn.mc2nu_pos, false);
  gkyl_nodal_ops_n2m(n2m, &up->basis, &up->grid, nrange, &up->local, 1, up->geo_corn.bmag_nodal, up->geo_corn.bmag, false);
  gkyl_nodal_ops_release(n2m);

  // Need 1/B for LBO collisions, computed weakly.
  gkyl_dg_inv_op_range(&up->basis, 0, up->geo_corn.bmag_inv, 0, up->geo_corn.bmag, &up->local);

}

static
void gk_geometry_mapc2p_advance_interior(struct gk_geometry* up, struct gkyl_range *nrange, double dzc[3], 
  evalf_t mapc2p_func, void* mapc2p_ctx, evalf_t bfield_func, void *bfield_ctx, struct gkyl_position_map *position_map)
{

  //Now project mapc2p and the FD array
  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates
                                
  double dtheta = up->grid.dx[TH_IDX],
    dpsi = up->grid.dx[PSI_IDX],
    dalpha = up->grid.dx[AL_IDX];

  double theta_lo = up->grid.lower[TH_IDX] + (up->local.lower[TH_IDX] - up->global.lower[TH_IDX])*up->grid.dx[TH_IDX],
    psi_lo = up->grid.lower[PSI_IDX] + (up->local.lower[PSI_IDX] - up->global.lower[PSI_IDX])*up->grid.dx[PSI_IDX],
    alpha_lo = up->grid.lower[AL_IDX] + (up->local.lower[AL_IDX] - up->global.lower[AL_IDX])*up->grid.dx[AL_IDX];

  double dels[2] = {1.0/sqrt(3), 1.0-1.0/sqrt(3) };
  theta_lo = theta_lo + dels[1]*dtheta/2.0;
  psi_lo = psi_lo + dels[1]*dpsi/2.0;
  alpha_lo = alpha_lo + dels[1]*dalpha/2.0;

  double dx_fact = up->basis.poly_order == 1 ? 1 : 0.5;
  dtheta *= dx_fact; dpsi *= dx_fact; dalpha *= dx_fact;

  // used for finite differences 
  double delta_alpha = dalpha*1e-2;
  double delta_psi = dpsi*1e-2;
  double delta_theta = dtheta*1e-2;
  dzc[0] = delta_psi;
  dzc[1] = delta_alpha;
  dzc[2] = delta_theta;
  int modifiers[5] = {0, -1, 1, -2, 2};
                                
  int cidx[3] = { 0 };
  for(int ia=nrange->lower[AL_IDX]; ia<=nrange->upper[AL_IDX]; ++ia){
    cidx[AL_IDX] = ia;
    for(int ia_delta = 0; ia_delta < 3; ia_delta++){ // interior stencil
      double alpha_curr = calc_running_coord(alpha_lo, ia-nrange->lower[AL_IDX], dalpha) + modifiers[ia_delta]*delta_alpha;
            position_map->maps[1](0.0, &alpha_curr, &alpha_curr, position_map->ctxs[1]);

      for (int ip=nrange->lower[PSI_IDX]; ip<=nrange->upper[PSI_IDX]; ++ip) {
        int ip_delta_max = 3;// interior
        if(ia_delta != 0)
          ip_delta_max = 1;
        for(int ip_delta = 0; ip_delta < ip_delta_max; ip_delta++){
          double psi_curr = calc_running_coord(psi_lo, ip-nrange->lower[PSI_IDX], dpsi) + modifiers[ip_delta]*delta_psi;
          position_map->maps[0](0.0, &psi_curr,   &psi_curr,   position_map->ctxs[0]);
          cidx[PSI_IDX] = ip;
          // set node coordinates
          for (int it=nrange->lower[TH_IDX]; it<=nrange->upper[TH_IDX]; ++it) {
            int it_delta_max = 3; // interior
            if(ia_delta != 0 || ip_delta != 0 )
              it_delta_max = 1;
            for(int it_delta = 0; it_delta < it_delta_max; it_delta++){
              double theta_curr = calc_running_coord(theta_lo, it-nrange->lower[TH_IDX], dtheta) + modifiers[it_delta]*delta_theta;
              position_map->maps[2](0.0, &theta_curr, &theta_curr, position_map->ctxs[2]);

              cidx[TH_IDX] = it;
              int lidx = 0;
              if (ip_delta != 0)
                lidx = 3 + 3*(ip_delta-1);
              if (ia_delta != 0)
                lidx = 15 + 3*(ia_delta-1);
              if (it_delta != 0)
                lidx = 27 + 3*(it_delta-1);

              double *mc2p_quad_fd_n = (double *) gkyl_array_fetch(up->geo_int.mc2p_nodal_fd, gkyl_range_idx(nrange, cidx));
              double *mc2p_quad_n = (double *) gkyl_array_fetch(up->geo_int.mc2p_nodal, gkyl_range_idx(nrange, cidx));
              double *bmag_n = (double *) gkyl_array_fetch(up->geo_int.bmag_nodal, gkyl_range_idx(nrange, cidx));
              double *bhat_quad_fd_n = (double *) gkyl_array_fetch(up->geo_int.b_i_nodal_fd, gkyl_range_idx(nrange, cidx));

              double xyz[3] = {psi_curr, alpha_curr, theta_curr};
              double XYZ[3] = {0.};
              double B[3] = {0.};
              mapc2p_func(0.0, xyz, XYZ, mapc2p_ctx);
              bfield_func(0.0, xyz, B, bfield_ctx);

              mc2p_quad_fd_n[lidx+X_IDX] = XYZ[X_IDX];
              mc2p_quad_fd_n[lidx+Y_IDX] = XYZ[Y_IDX];
              mc2p_quad_fd_n[lidx+Z_IDX] = XYZ[Z_IDX];

              double bmag = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
              bhat_quad_fd_n[lidx+X_IDX] = B[X_IDX]/bmag;
              bhat_quad_fd_n[lidx+Y_IDX] = B[Y_IDX]/bmag;
              bhat_quad_fd_n[lidx+Z_IDX] = B[Z_IDX]/bmag;


              if(ip_delta==0 && ia_delta==0 && it_delta==0){
                mc2p_quad_n[X_IDX] = XYZ[X_IDX];
                mc2p_quad_n[Y_IDX] = XYZ[Y_IDX];
                mc2p_quad_n[Z_IDX] = XYZ[Z_IDX];
                bmag_n[0] = bmag;
              }
            }
          }
        }
      }
    }
  }

  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&up->basis, &up->grid, false);
  gkyl_nodal_ops_n2m(n2m, &up->basis, &up->grid, nrange, &up->local, 3, up->geo_int.mc2p_nodal, up->geo_int.mc2p, true);
  gkyl_nodal_ops_n2m(n2m, &up->basis, &up->grid, nrange, &up->local, 1, up->geo_int.bmag_nodal, up->geo_int.bmag, true);
  gkyl_nodal_ops_release(n2m);

  // Now calculate the metrics.
  struct gkyl_calc_metric* mcalc = gkyl_calc_metric_new(&up->basis, &up->grid,
    &up->global, &up->global_ext, &up->local, &up->local_ext, false, false);
  gkyl_calc_metric_advance_interior(mcalc, up);
  gkyl_array_copy(up->geo_int.g_ij_neut, up->geo_int.g_ij);
  
  // Calculate the derived geometric quantities.
  struct gkyl_calc_derived_geo *jcalculator = gkyl_calc_derived_geo_new(&up->basis, &up->grid, 1, false);
  gkyl_calc_derived_geo_advance(jcalculator, &up->local, up->geo_int.g_ij, up->geo_int.bmag, 
    up->geo_int.jacobgeo, up->geo_int.jacobgeo_inv, up->geo_int.gij, up->geo_int.b_i, up->geo_int.cmag, up->geo_int.jacobtot, up->geo_int.jacobtot_inv, 
    up->geo_int.gxxj, up->geo_int.gxyj, up->geo_int.gyyj, up->geo_int.gxzj, up->geo_int.eps2);
  gkyl_array_copy(up->geo_int.gij_neut, up->geo_int.gij);
  gkyl_calc_derived_geo_release(jcalculator);
  gkyl_calc_metric_advance_bcart(mcalc, nrange, up->geo_int.b_i, up->geo_int.dzdx, up->geo_int.bcart, &up->local);
  gkyl_calc_metric_release(mcalc);
}

void gk_geometry_mapc2p_advance_surface(struct gk_geometry* up, int dir, struct gkyl_range *nrange, double dzc[3], 
  evalf_t mapc2p_func, void* mapc2p_ctx, evalf_t bfield_func, void *bfield_ctx, 
  struct gkyl_position_map *position_map)
{

  //Now project mapc2p and the FD array
  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates
                                
  double dtheta = up->grid.dx[TH_IDX],
    dpsi = up->grid.dx[PSI_IDX],
    dalpha = up->grid.dx[AL_IDX];

  double theta_lo = up->grid.lower[TH_IDX] + (up->local.lower[TH_IDX] - up->global.lower[TH_IDX])*up->grid.dx[TH_IDX],
    psi_lo = up->grid.lower[PSI_IDX] + (up->local.lower[PSI_IDX] - up->global.lower[PSI_IDX])*up->grid.dx[PSI_IDX],
    alpha_lo = up->grid.lower[AL_IDX] + (up->local.lower[AL_IDX] - up->global.lower[AL_IDX])*up->grid.dx[AL_IDX];

  double dels[2] = {1.0/sqrt(3), 1.0-1.0/sqrt(3) };
  theta_lo += dir == 2 ? 0.0 : dels[1]*dtheta/2.0;
  psi_lo += dir == 0 ? 0.0 : dels[1]*dpsi/2.0;
  alpha_lo += dir == 1 ? 0. : dels[1]*dalpha/2.0;

  double dx_fact = up->basis.poly_order == 1 ? 1 : 0.5;
  dtheta *= dx_fact; dpsi *= dx_fact; dalpha *= dx_fact;

  // used for finite differences 
  double delta_alpha = dalpha*1e-2;
  double delta_psi = dpsi*1e-2;
  double delta_theta = dtheta*1e-2;
  dzc[0] = delta_psi;
  dzc[1] = delta_alpha;
  dzc[2] = delta_theta;
  int modifiers[5] = {0, -1, 1, -2, 2};
                                
  int cidx[3] = { 0 };
  for(int ia=nrange->lower[AL_IDX]; ia<=nrange->upper[AL_IDX]; ++ia){
    cidx[AL_IDX] = ia;
    for(int ia_delta = 0; ia_delta < 5; ia_delta++){ // should be <5
      if((ia == nrange->lower[AL_IDX]) && (up->local.lower[AL_IDX]== up->global.lower[AL_IDX]) && dir==1){
        if(ia_delta == 1 || ia_delta == 3)
          continue; // want to use one sided stencils at edge
      }
      else if((ia == nrange->upper[AL_IDX])  && (up->local.upper[AL_IDX]== up->global.upper[AL_IDX])&& dir==1){
          if(ia_delta == 2 || ia_delta == 4)
            continue; // want to use one sided stencils at edge
      }
      else{ //interior
        if( ia_delta == 3 || ia_delta == 4)
          continue; //dont do two away
      }
      double alpha_curr = dir==1 ? alpha_lo + ia*dalpha : calc_running_coord(alpha_lo, ia-nrange->lower[AL_IDX], dalpha);
      alpha_curr += modifiers[ia_delta]*delta_alpha;
      position_map->maps[1](0.0, &alpha_curr, &alpha_curr, position_map->ctxs[1]);

      for (int ip=nrange->lower[PSI_IDX]; ip<=nrange->upper[PSI_IDX]; ++ip) {
        int ip_delta_max = 5;// should be 5
        if(ia_delta != 0)
          ip_delta_max = 1;
        for(int ip_delta = 0; ip_delta < ip_delta_max; ip_delta++){
          if((ip == nrange->lower[PSI_IDX]) && (up->local.lower[PSI_IDX]== up->global.lower[PSI_IDX]) && dir==0){
            if(ip_delta == 1 || ip_delta == 3)
              continue; // want to use one sided stencils at edge
          }
          else if((ip == nrange->upper[PSI_IDX]) && (up->local.upper[PSI_IDX]== up->global.upper[PSI_IDX]) && dir==0){
            if(ip_delta == 2 || ip_delta == 4)
              continue; // want to use one sided stencils at edge
          }
          else{ // interior 
            if( ip_delta == 3 || ip_delta == 4)
              continue; //dont do two away
          }
          double psi_curr = dir == 0 ? psi_lo + ip*dpsi : calc_running_coord(psi_lo, ip-nrange->lower[PSI_IDX], dpsi) ;
          psi_curr += modifiers[ip_delta]*delta_psi;
          position_map->maps[0](0.0, &psi_curr,   &psi_curr,   position_map->ctxs[0]);
          cidx[PSI_IDX] = ip;
          // set node coordinates
          for (int it=nrange->lower[TH_IDX]; it<=nrange->upper[TH_IDX]; ++it) {
            int it_delta_max = 5; // should be 5
            if(ia_delta != 0 || ip_delta != 0 )
              it_delta_max = 1;
            for(int it_delta = 0; it_delta < it_delta_max; it_delta++){
              if((it == nrange->lower[TH_IDX]) && (up->local.lower[TH_IDX]== up->global.lower[TH_IDX]) && dir==2){
                if(it_delta == 1 || it_delta == 3)
                  continue; // want to use one sided stencils at edge
              }
              else if((it == nrange->upper[TH_IDX]) && (up->local.upper[TH_IDX]== up->global.upper[TH_IDX]) && dir==2){
                if(it_delta == 2 || it_delta == 4)
                  continue; // want to use one sided stencils at edge
              }
              else{
                if( it_delta == 3 || it_delta == 4)
                  continue; //dont do two away
              }
              double theta_curr = dir==2 ? theta_lo + it*dtheta: calc_running_coord(theta_lo, it-nrange->lower[TH_IDX], dtheta);
              theta_curr += modifiers[it_delta]*delta_theta;
              position_map->maps[2](0.0, &theta_curr, &theta_curr, position_map->ctxs[2]);

              cidx[TH_IDX] = it;
              int lidx = 0;
              if (ip_delta != 0)
                lidx = 3 + 3*(ip_delta-1);
              if (ia_delta != 0)
                lidx = 15 + 3*(ia_delta-1);
              if (it_delta != 0)
                lidx = 27 + 3*(it_delta-1);

              double *mc2p_fd_n = (double *) gkyl_array_fetch(up->geo_surf[dir].mc2p_nodal_fd, gkyl_range_idx(nrange, cidx));
              double *mc2p_n = (double *) gkyl_array_fetch(up->geo_surf[dir].mc2p_nodal, gkyl_range_idx(nrange, cidx));
              double *bmag_n = (double *) gkyl_array_fetch(up->geo_surf[dir].bmag_nodal, gkyl_range_idx(nrange, cidx));
              double *bhat_fd_n = (double *) gkyl_array_fetch(up->geo_surf[dir].b_i_nodal_fd, gkyl_range_idx(nrange, cidx));

              double xyz[3] = {psi_curr, alpha_curr, theta_curr};
              double XYZ[3] = {0.};
              double B[3] = {0.};
              mapc2p_func(0.0, xyz, XYZ, mapc2p_ctx);
              bfield_func(0.0, xyz, B, bfield_ctx);

              mc2p_fd_n[lidx+X_IDX] = XYZ[X_IDX];
              mc2p_fd_n[lidx+Y_IDX] = XYZ[Y_IDX];
              mc2p_fd_n[lidx+Z_IDX] = XYZ[Z_IDX];

              double bmag = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
              bhat_fd_n[lidx+X_IDX] = B[X_IDX]/bmag;
              bhat_fd_n[lidx+Y_IDX] = B[Y_IDX]/bmag;
              bhat_fd_n[lidx+Z_IDX] = B[Z_IDX]/bmag;

              if(ip_delta==0 && ia_delta==0 && it_delta==0){
                mc2p_n[X_IDX] = XYZ[X_IDX];
                mc2p_n[Y_IDX] = XYZ[Y_IDX];
                mc2p_n[Z_IDX] = XYZ[Z_IDX];
                bmag_n[0] = bmag;
              }
            }
          }
        }
      }
    }
  }

  // Now calculate the metrics.
  struct gkyl_calc_metric* mcalc = gkyl_calc_metric_new(&up->basis, &up->grid,
    &up->global, &up->global_ext, &up->local, &up->local_ext, false, false);
  gkyl_calc_metric_advance_surface(mcalc, dir, up);
  gkyl_calc_metric_release(mcalc);
  gk_geometry_surf_calc_expansions(up, dir, *nrange);
}

struct gk_geometry*
gk_geometry_mapc2p_init(struct gkyl_gk_geometry_inp *geometry_inp)
{

  struct gk_geometry *up = gkyl_malloc(sizeof(struct gk_geometry));
  up->geometry_id = geometry_inp->geometry_id;
  up->basis = geometry_inp->geo_basis;
  up->local = geometry_inp->geo_local;
  up->local_ext = geometry_inp->geo_local_ext;
  up->global = geometry_inp->geo_global;
  up->global_ext = geometry_inp->geo_global_ext;
  up->grid = geometry_inp->geo_grid;

  up->geqdsk_sign_convention = 0;
  up->half_domain = 0;

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
  up->parallel_lower_bc_shift_ctx  = geometry_inp->parallel_lower_bc_shift_ctx ;
  up->parallel_upper_bc_shift_ctx  = geometry_inp->parallel_upper_bc_shift_ctx ;

  gk_geometry_set_nodal_ranges(up) ;

  // Initialize surface basis abd allocate surface geo
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

  gkyl_position_map_optimize(geometry_inp->position_map, up->grid, up->global);

  // calculate mapc2p in cartesian coords at corner nodes for
  // getting cell coordinates (used only for plotting)
  gk_geometry_mapc2p_advance(up, &up->nrange_corn, up->dzc, geometry_inp->mapc2p, geometry_inp->c2p_ctx,
    geometry_inp->bfield_func, geometry_inp->bfield_ctx, geometry_inp->position_map);
  // calculate mapc2p in cartesian coords at interior nodes for
  // calculating geo quantity volume expansions 
  gk_geometry_mapc2p_advance_interior(up, &up->nrange_int, up->dzc, geometry_inp->mapc2p, geometry_inp->c2p_ctx,
    geometry_inp->bfield_func, geometry_inp->bfield_ctx, geometry_inp->position_map);
  // calculate mapc2p in cylindrical coords at surfaces
  for (int dir = 0; dir <up->grid.ndim; dir++) {
    gk_geometry_mapc2p_advance_surface(up, dir, &up->nrange_surf[dir], up->dzc, geometry_inp->mapc2p, 
      geometry_inp->c2p_ctx, geometry_inp->bfield_func, geometry_inp->bfield_ctx, geometry_inp->position_map);
  }

  // Store metadata for I/O.
  struct gkyl_msgpack_map_elem io_meta_basic[] = {
    { .key = "geometry_type", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = up->geometry_id },
  };
  up->io_meta_basic_len = sizeof(io_meta_basic)/sizeof(io_meta_basic[0]);
  up->io_meta_basic = gkyl_msgpack_map_elem_clone(up->io_meta_basic_len, io_meta_basic);

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->ref_count = gkyl_ref_count_init(gkyl_gk_geometry_free);
  up->on_dev = up; // CPU eqn obj points to itself

  return up;
}

struct gk_geometry*
gkyl_gk_geometry_mapc2p_new(struct gkyl_gk_geometry_inp *geometry_inp)
{
  struct gk_geometry* gk_geom_3d;
  struct gk_geometry* gk_geom;

  // First construct the uniform 3d geometry
  gk_geom_3d = gk_geometry_mapc2p_init(geometry_inp);
  if (geometry_inp->position_map->id == GKYL_PMAP_CONSTANT_DB_POLYNOMIAL || \
      geometry_inp->position_map->id == GKYL_PMAP_CONSTANT_DB_NUMERIC) {
    // The array mc2nu is computed using the uniform geometry, so we need to deflate it
    // Must deflate the 3D uniform geometry in order for the allgather to work
    if (geometry_inp->grid.ndim < 3)
      gk_geom = gkyl_gk_geometry_deflate(gk_geom_3d, geometry_inp);
    else
      gk_geom = gkyl_gk_geometry_acquire(gk_geom_3d);

    gkyl_position_map_set_bmag(geometry_inp->position_map, geometry_inp->comm, \
      gk_geom->geo_int.bmag);

    gkyl_gk_geometry_release(gk_geom_3d); // release temporary 3d geometry
    gkyl_gk_geometry_release(gk_geom); // release 3d geometry

    // Construct the non-uniform grid
    gk_geom_3d = gk_geometry_mapc2p_init(geometry_inp);
  }
  return gk_geom_3d;
}
