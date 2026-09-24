#include <gkyl_calc_metric_mirror.h>
#include <gkyl_calc_metric_mirror_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_nodal_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_mirror_grid_gen.h>
#include <gkyl_math.h>
#include <gkyl_dg_bin_ops.h>

gkyl_calc_metric_mirror*
gkyl_calc_metric_mirror_new(const struct gkyl_basis *cbasis, const struct gkyl_rect_grid *grid,
  const struct gkyl_range *local, const struct gkyl_range *local_ext, bool use_gpu)
{
  gkyl_calc_metric_mirror *up = gkyl_malloc(sizeof(gkyl_calc_metric_mirror));
  up->cbasis = cbasis;
  up->cdim = cbasis->ndim;
  up->cnum_basis = cbasis->num_basis;
  up->poly_order = cbasis->poly_order;
  up->grid = grid;
  up->use_gpu = use_gpu;
  up->num_cells = up->grid->cells;
  up->n2m = gkyl_nodal_ops_new(up->cbasis, up->grid, up->use_gpu);

  up->local = *local;
  up->local_ext = *local_ext;
  return up;
}

static inline double calc_metric(double dxdz[3][3], int i, int j) 
{
  double sum = 0;
  for (int k = 0; k < 3; ++k)
    sum += dxdz[k][i - 1] * dxdz[k][j - 1];
  return sum;
} 

// Calculates e^1 = e_2 x e_3 /J
static inline void
calc_dual(double J, const double e_2[3], const double e_3[3], double e1[3])
{
  e1[0] = (e_2[1]*e_3[2] - e_2[2]*e_3[1] )/J;
  e1[1] = -(e_2[0]*e_3[2] - e_2[2]*e_3[0] )/J;
  e1[2] = (e_2[0]*e_3[1] - e_2[1]*e_3[0] )/J;
}

void gkyl_calc_metric_mirror_advance(
  gkyl_calc_metric_mirror *up, struct gk_geometry *gk_geom, struct gkyl_mirror_grid_gen *mirror_grid)
{
  enum { PSI_IDX, AL_IDX, TH_IDX }; // Arrangement of computational coordinates.
  enum { R_IDX, Z_IDX, PHI_IDX }; // Arrangement of cartesian coordinates.
  int cidx[3];
  for(int ia=gk_geom->nrange_corn.lower[AL_IDX]; ia<=gk_geom->nrange_corn.upper[AL_IDX]; ++ia){
    for (int ip=gk_geom->nrange_corn.lower[PSI_IDX]; ip<=gk_geom->nrange_corn.upper[PSI_IDX]; ++ip) {
      for (int it=gk_geom->nrange_corn.lower[TH_IDX]; it<=gk_geom->nrange_corn.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;

        // First fetch the mirror stuff at this location
        const double *mirror_rza_n = gkyl_array_cfetch(mirror_grid->nodes_rza, gkyl_range_idx(&gk_geom->nrange_corn, cidx));
        const double *mirror_psi_n = gkyl_array_cfetch(mirror_grid->nodes_psi, gkyl_range_idx(&gk_geom->nrange_corn, cidx));
        const struct gkyl_mirror_grid_gen_geom *mirror_geo_n = gkyl_array_cfetch(mirror_grid->nodes_geom, gkyl_range_idx(&gk_geom->nrange_corn, cidx));

        // Next fetch the gk_geometry nodal values at this location
        double *mc2p_n = gkyl_array_fetch(gk_geom->geo_corn.mc2p_nodal, gkyl_range_idx(&gk_geom->nrange_corn, cidx));
        double *mc2nu_pos_n = gkyl_array_fetch(gk_geom->geo_corn.mc2nu_pos_nodal, gkyl_range_idx(&gk_geom->nrange_corn, cidx));
        double *bmag_n = gkyl_array_fetch(gk_geom->geo_corn.bmag_nodal, gkyl_range_idx(&gk_geom->nrange_corn, cidx));

        // Set mapc2p
        mc2p_n[0] = mirror_rza_n[0]; // R
        mc2p_n[1] = mirror_rza_n[1]; // Z
        mc2p_n[2] = mirror_rza_n[2]; // PHI
        // Set mapc2nu
        mc2nu_pos_n[0] = mirror_psi_n[0];
        mc2nu_pos_n[1] = mirror_rza_n[2];
        mc2nu_pos_n[2] = mirror_rza_n[1];
        // Set bmag
        struct gkyl_vec3 B_cart = gkyl_vec3_polar_con_to_cart(mirror_rza_n[0], mirror_rza_n[2], mirror_geo_n->B);
        bmag_n[0] = gkyl_vec3_len(B_cart);
      }
    }
  }
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_corn, &gk_geom->local, 3,
    gk_geom->geo_corn.mc2p_nodal, gk_geom->geo_corn.mc2p, false);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_corn, &gk_geom->local, 3,
    gk_geom->geo_corn.mc2nu_pos_nodal, gk_geom->geo_corn.mc2nu_pos, false);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_corn, &gk_geom->local, 1,
    gk_geom->geo_corn.bmag_nodal, gk_geom->geo_corn.bmag, false);

  // Need 1/B for LBO collisions, computed weakly.
  gkyl_dg_inv_op_range(up->cbasis, 0, gk_geom->geo_corn.bmag_inv, 0, gk_geom->geo_corn.bmag, &gk_geom->local);

}

void 
gkyl_calc_metric_mirror_advance_interior( gkyl_calc_metric_mirror *up, struct gk_geometry *gk_geom, struct gkyl_mirror_grid_gen *mirror_grid)
{
  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { R_IDX, Z_IDX, PHI_IDX }; // arrangement of cartesian coordinates
  int cidx[3];
  for(int ia=gk_geom->nrange_int.lower[AL_IDX]; ia<=gk_geom->nrange_int.upper[AL_IDX]; ++ia){
    for (int ip=gk_geom->nrange_int.lower[PSI_IDX]; ip<=gk_geom->nrange_int.upper[PSI_IDX]; ++ip) {
      for (int it=gk_geom->nrange_int.lower[TH_IDX]; it<=gk_geom->nrange_int.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;

        // First fetch the mirror stuff at this location
        const double *mirror_rza_n = gkyl_array_cfetch(mirror_grid->nodes_rza, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        const double *mirror_psi_n = gkyl_array_cfetch(mirror_grid->nodes_psi, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        const struct gkyl_mirror_grid_gen_geom *mirror_geo_n = gkyl_array_cfetch(mirror_grid->nodes_geom, gkyl_range_idx(&gk_geom->nrange_int, cidx));

        // Next fetch the gk_geometry nodal values at this location
        double *mc2p_n = gkyl_array_fetch(gk_geom->geo_int.mc2p_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        double *jFld_n= gkyl_array_fetch(gk_geom->geo_int.jacobgeo_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        double *bmag_n = gkyl_array_fetch(gk_geom->geo_int.bmag_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        double *gFld_n= gkyl_array_fetch(gk_geom->geo_int.g_ij_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        double *bcartFld_n= gkyl_array_fetch(gk_geom->geo_int.bcart_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        double *tanvecFld_n= gkyl_array_fetch(gk_geom->geo_int.dxdz_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        double *dualFld_n= gkyl_array_fetch(gk_geom->geo_int.dzdx_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        double *dualmagFld_n = gkyl_array_fetch(gk_geom->geo_int.dualmag_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        double *normFld_n = gkyl_array_fetch(gk_geom->geo_int.normals_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        double *curlbhat_n = gkyl_array_fetch(gk_geom->geo_int.curlbhat_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        double *dualcurlbhat_n = gkyl_array_fetch(gk_geom->geo_int.dualcurlbhat_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        double *B3_n = gkyl_array_fetch(gk_geom->geo_int.B3_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        double *dualcurlbhatoverB_n = gkyl_array_fetch(gk_geom->geo_int.dualcurlbhatoverB_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        double *rtg33inv_n = gkyl_array_fetch(gk_geom->geo_int.rtg33inv_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));
        double *bioverJB_n = gkyl_array_fetch(gk_geom->geo_int.bioverJB_nodal, gkyl_range_idx(&gk_geom->nrange_int, cidx));

        // Set mapc2p
        mc2p_n[0] = mirror_rza_n[0]; // R
        mc2p_n[1] = mirror_rza_n[1]; // Z
        mc2p_n[2] = mirror_rza_n[2]; // PHI
        // Set bmag
        struct gkyl_vec3 B_cart = gkyl_vec3_polar_con_to_cart(mirror_rza_n[0], mirror_rza_n[2], mirror_geo_n->B);
        bmag_n[0] = gkyl_vec3_len(B_cart);
        // Set Jacobian
        jFld_n[0] = mirror_geo_n->Jc;
        // Set curl(bhat)
        struct gkyl_vec3 curlbhat_cart = gkyl_vec3_polar_con_to_cart(mirror_rza_n[0], mirror_rza_n[2], mirror_geo_n->curlbhat);
        curlbhat_n[0] = curlbhat_cart.x[0];
        curlbhat_n[1] = curlbhat_cart.x[1];
        curlbhat_n[2] = curlbhat_cart.x[2];

        // Set tangents and duals
        for(int ivec=0; ivec<3; ivec++) {
          struct gkyl_vec3 tan = gkyl_vec3_polar_con_to_cart(mirror_rza_n[0], mirror_rza_n[2], mirror_geo_n->tang[ivec]);
          struct gkyl_vec3 dual = gkyl_vec3_polar_con_to_cart(mirror_rza_n[0], mirror_rza_n[2], mirror_geo_n->dual[ivec]);
          tanvecFld_n[ivec*3 + 0] = tan.x[0];
          tanvecFld_n[ivec*3 + 1] = tan.x[1];
          tanvecFld_n[ivec*3 + 2] = tan.x[2];
          dualFld_n[ivec*3 + 0] = dual.x[0];
          dualFld_n[ivec*3 + 1] = dual.x[1];
          dualFld_n[ivec*3 + 2] = dual.x[2];
        }
        // Set metric tensor
        int count = 0;
        for (int i=0; i<3; ++i) {
          for (int j=0; j<3; ++j) {
            if (i > j)
              continue;
            gFld_n[count] = tanvecFld_n[3*i+0]*tanvecFld_n[3*j+0] + tanvecFld_n[3*i+1]*tanvecFld_n[3*j+1] + tanvecFld_n[3*i+2]*tanvecFld_n[3*j+2];
            count++;
          }
        }


        // Now set derived geometric quantities
        double R = mc2p_n[R_IDX];
        double Z = mc2p_n[Z_IDX];
        double phi = mc2p_n[PHI_IDX];
        double J = jFld_n[0];

        // Calculate cartesian components of bhat
        double bi[3];
        bi[0] = gFld_n[2]/sqrt(gFld_n[5]);
        bi[1] = gFld_n[4]/sqrt(gFld_n[5]);
        bi[2] = gFld_n[5]/sqrt(gFld_n[5]);
        bcartFld_n[0] = dualFld_n[0]*bi[0] + dualFld_n[3]*bi[1] + dualFld_n[6]*bi[2];
        bcartFld_n[1] = dualFld_n[1]*bi[0] + dualFld_n[4]*bi[1] + dualFld_n[7]*bi[2];
        bcartFld_n[2] = dualFld_n[2]*bi[0] + dualFld_n[5]*bi[1] + dualFld_n[8]*bi[2];

        double norm1 = sqrt(dualFld_n[0]*dualFld_n[0] + dualFld_n[1]*dualFld_n[1] + dualFld_n[2]*dualFld_n[2]);
        double norm2 = sqrt(dualFld_n[3]*dualFld_n[3] + dualFld_n[4]*dualFld_n[4] + dualFld_n[5]*dualFld_n[5]);
        double norm3 = sqrt(dualFld_n[6]*dualFld_n[6] + dualFld_n[7]*dualFld_n[7] + dualFld_n[8]*dualFld_n[8]);

        dualmagFld_n[0] = norm1;
        dualmagFld_n[1] = norm2;
        dualmagFld_n[2] = norm3;
        
        // Set normal vectors
        normFld_n[0] = dualFld_n[0]/norm1;
        normFld_n[1] = dualFld_n[1]/norm1;
        normFld_n[2] = dualFld_n[2]/norm1;

        normFld_n[3] = dualFld_n[3]/norm2;
        normFld_n[4] = dualFld_n[4]/norm2;
        normFld_n[5] = dualFld_n[5]/norm2;

        normFld_n[6] = dualFld_n[6]/norm3;
        normFld_n[7] = dualFld_n[7]/norm3;
        normFld_n[8] = dualFld_n[8]/norm3;

        // Set e^m \dot curl(bhat) 
        dualcurlbhat_n[0] = dualFld_n[0]*curlbhat_n[0] +  dualFld_n[1]*curlbhat_n[1] + dualFld_n[2]*curlbhat_n[2];
        dualcurlbhat_n[1] = dualFld_n[3]*curlbhat_n[0] +  dualFld_n[4]*curlbhat_n[1] + dualFld_n[5]*curlbhat_n[2];
        dualcurlbhat_n[2] = dualFld_n[6]*curlbhat_n[0] +  dualFld_n[7]*curlbhat_n[1] + dualFld_n[8]*curlbhat_n[2];

        // Set e^3 \dot B 
        B3_n[0] = bmag_n[0]/sqrt(gFld_n[5]);

        // set e^3 \dot B /|B|
        dualcurlbhatoverB_n[0] = dualcurlbhat_n[0]/bmag_n[0];
        dualcurlbhatoverB_n[1] = dualcurlbhat_n[1]/bmag_n[0];
        dualcurlbhatoverB_n[2] = dualcurlbhat_n[2]/bmag_n[0];

        // set B^3/B = 1/sqrt(g_33)
        rtg33inv_n[0] = 1.0/sqrt(gFld_n[5]);

        // set b_i/JB
        bioverJB_n[0] = gFld_n[2]/sqrt(gFld_n[5])/J/bmag_n[0];
        bioverJB_n[1] = gFld_n[4]/sqrt(gFld_n[5])/J/bmag_n[0];
        bioverJB_n[2] = gFld_n[5]/sqrt(gFld_n[5])/J/bmag_n[0];
      }
    }
  }

  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &gk_geom->local, 3, gk_geom->geo_int.mc2p_nodal, gk_geom->geo_int.mc2p, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &gk_geom->local, 1, gk_geom->geo_int.bmag_nodal, gk_geom->geo_int.bmag, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &gk_geom->local, 6, gk_geom->geo_int.g_ij_nodal, gk_geom->geo_int.g_ij, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &gk_geom->local, 1, gk_geom->geo_int.jacobgeo_nodal, gk_geom->geo_int.jacobgeo, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &gk_geom->local, 3, gk_geom->geo_int.bcart_nodal, gk_geom->geo_int.bcart, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &gk_geom->local, 9, gk_geom->geo_int.dxdz_nodal, gk_geom->geo_int.dxdz, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &gk_geom->local, 9, gk_geom->geo_int.dzdx_nodal, gk_geom->geo_int.dzdx, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &gk_geom->local, 3, gk_geom->geo_int.dualmag_nodal, gk_geom->geo_int.dualmag, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &gk_geom->local, 9, gk_geom->geo_int.normals_nodal, gk_geom->geo_int.normals, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &gk_geom->local, 3, gk_geom->geo_int.dualcurlbhat_nodal, gk_geom->geo_int.dualcurlbhat, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &gk_geom->local, 3, gk_geom->geo_int.dualcurlbhatoverB_nodal, gk_geom->geo_int.dualcurlbhatoverB, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &gk_geom->local, 1, gk_geom->geo_int.rtg33inv_nodal, gk_geom->geo_int.rtg33inv, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &gk_geom->local, 3, gk_geom->geo_int.bioverJB_nodal, gk_geom->geo_int.bioverJB, true);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, &gk_geom->nrange_int, &gk_geom->local, 1, gk_geom->geo_int.B3_nodal, gk_geom->geo_int.B3, true);

}

void gkyl_calc_metric_mirror_advance_surface( gkyl_calc_metric_mirror *up, int dir, struct gk_geometry *gk_geom, struct gkyl_mirror_grid_gen *mirror_grid)
{
  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { R_IDX, Z_IDX, PHI_IDX }; // arrangement of cartesian coordinates
  int cidx[3];
  for(int ia=gk_geom->nrange_surf[dir].lower[AL_IDX]; ia<=gk_geom->nrange_surf[dir].upper[AL_IDX]; ++ia){
    for (int ip=gk_geom->nrange_surf[dir].lower[PSI_IDX]; ip<=gk_geom->nrange_surf[dir].upper[PSI_IDX]; ++ip) {
      for (int it=gk_geom->nrange_surf[dir].lower[TH_IDX]; it<=gk_geom->nrange_surf[dir].upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;

        // First fetch the mirror stuff at this location
        const double *mirror_rza_n = gkyl_array_cfetch(mirror_grid->nodes_rza, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        const double *mirror_psi_n = gkyl_array_cfetch(mirror_grid->nodes_psi, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        const struct gkyl_mirror_grid_gen_geom *mirror_geo_n = gkyl_array_cfetch(mirror_grid->nodes_geom, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));

        // Next fetch the gk_geometry nodal values at this location
        double *mc2p_n = gkyl_array_fetch(gk_geom->geo_surf[dir].mc2p_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        double *jFld_n= gkyl_array_fetch(gk_geom->geo_surf[dir].jacobgeo_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        double *bmag_n = gkyl_array_fetch(gk_geom->geo_surf[dir].bmag_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        double *gFld_n= gkyl_array_fetch(gk_geom->geo_surf[dir].g_ij_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        double *biFld_n= gkyl_array_fetch(gk_geom->geo_surf[dir].b_i_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        double *cmagFld_n= gkyl_array_fetch(gk_geom->geo_surf[dir].cmag_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        double *bcartFld_n= gkyl_array_fetch(gk_geom->geo_surf[dir].bcart_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        double *tanvecFld_n= gkyl_array_fetch(gk_geom->geo_surf[dir].dxdz_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        double *dualFld_n= gkyl_array_fetch(gk_geom->geo_surf[dir].dzdx_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        double *dualmagFld_n = gkyl_array_fetch(gk_geom->geo_surf[dir].dualmag_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        double *normFld_n = gkyl_array_fetch(gk_geom->geo_surf[dir].normals_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        double *lenr_n = gkyl_array_fetch(gk_geom->geo_surf[dir].lenr_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        double *B3_n = gkyl_array_fetch(gk_geom->geo_surf[dir].B3_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        double *curlbhat_n = gkyl_array_fetch(gk_geom->geo_surf[dir].curlbhat_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        double *normcurlbhat_n = gkyl_array_fetch(gk_geom->geo_surf[dir].normcurlbhat_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));

        // Set mapc2p
        mc2p_n[0] = mirror_rza_n[0]; // R
        mc2p_n[1] = mirror_rza_n[1]; // Z
        mc2p_n[2] = mirror_rza_n[2]; // PHI
        // Set bmag
        struct gkyl_vec3 B_cart = gkyl_vec3_polar_con_to_cart(mirror_rza_n[0], mirror_rza_n[2], mirror_geo_n->B);
        bmag_n[0] = gkyl_vec3_len(B_cart);
        // Set Jacobian
        jFld_n[0] = mirror_geo_n->Jc;
        // Set curl(bhat)
        struct gkyl_vec3 curlbhat_cart = gkyl_vec3_polar_con_to_cart(mirror_rza_n[0], mirror_rza_n[2], mirror_geo_n->curlbhat);
        curlbhat_n[0] = curlbhat_cart.x[0];
        curlbhat_n[1] = curlbhat_cart.x[1];
        curlbhat_n[2] = curlbhat_cart.x[2];
        // Set tangents and duals
        for(int ivec=0; ivec<3; ivec++) {
          struct gkyl_vec3 tan = gkyl_vec3_polar_con_to_cart(mirror_rza_n[0], mirror_rza_n[2], mirror_geo_n->tang[ivec]);
          struct gkyl_vec3 dual = gkyl_vec3_polar_con_to_cart(mirror_rza_n[0], mirror_rza_n[2], mirror_geo_n->dual[ivec]);
          tanvecFld_n[ivec*3 + 0] = tan.x[0];
          tanvecFld_n[ivec*3 + 1] = tan.x[1];
          tanvecFld_n[ivec*3 + 2] = tan.x[2];
          dualFld_n[ivec*3 + 0] = dual.x[0];
          dualFld_n[ivec*3 + 1] = dual.x[1];
          dualFld_n[ivec*3 + 2] = dual.x[2];
        }
        // Set metric tensor
        int count = 0;
        for (int i=0; i<3; ++i) {
          for (int j=0; j<3; ++j) {
            if (i > j)
              continue;
            gFld_n[count] = tanvecFld_n[3*i+0]*tanvecFld_n[3*j+0] + tanvecFld_n[3*i+1]*tanvecFld_n[3*j+1] + tanvecFld_n[3*i+2]*tanvecFld_n[3*j+2];
            count++;
          }
        }


        // Set derived quantities
        double R = mc2p_n[R_IDX];
        double Z = mc2p_n[Z_IDX];
        double J = jFld_n[0];
        double phi = mc2p_n[PHI_IDX];

        // Calculate cmag, bi, and jtot_inv
        biFld_n[0] = gFld_n[2]/sqrt(gFld_n[5]);
        biFld_n[1] = gFld_n[4]/sqrt(gFld_n[5]);
        biFld_n[2] = gFld_n[5]/sqrt(gFld_n[5]);

        cmagFld_n[0] = jFld_n[0]*bmag_n[0]/sqrt(gFld_n[5]);
        double *jtotinvFld_n= gkyl_array_fetch(gk_geom->geo_surf[dir].jacobtot_inv_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx));
        jtotinvFld_n[0] = 1.0/(jFld_n[0]*bmag_n[0]);

        double norm1 = sqrt(dualFld_n[0]*dualFld_n[0] + dualFld_n[1]*dualFld_n[1] + dualFld_n[2]*dualFld_n[2]);
        double norm2 = sqrt(dualFld_n[3]*dualFld_n[3] + dualFld_n[4]*dualFld_n[4] + dualFld_n[5]*dualFld_n[5]);
        double norm3 = sqrt(dualFld_n[6]*dualFld_n[6] + dualFld_n[7]*dualFld_n[7] + dualFld_n[8]*dualFld_n[8]);

        dualmagFld_n[0] = norm1;
        dualmagFld_n[1] = norm2;
        dualmagFld_n[2] = norm3;
        
        // Set normal vectors
        normFld_n[0] = dualFld_n[0]/norm1;
        normFld_n[1] = dualFld_n[1]/norm1;
        normFld_n[2] = dualFld_n[2]/norm1;

        normFld_n[3] = dualFld_n[3]/norm2;
        normFld_n[4] = dualFld_n[4]/norm2;
        normFld_n[5] = dualFld_n[5]/norm2;

        normFld_n[6] = dualFld_n[6]/norm3;
        normFld_n[7] = dualFld_n[7]/norm3;
        normFld_n[8] = dualFld_n[8]/norm3;

        // Set lenr
        lenr_n[0] = J*dualmagFld_n[dir];

        // Set n^3 \dot B 
        B3_n[0] = bmag_n[0]/sqrt(gFld_n[5])/norm3;

        // Set n^m \dot curl(bhat) 
        normcurlbhat_n[0] = normFld_n[3*dir+0]*curlbhat_n[0] +  normFld_n[3*dir+1]*curlbhat_n[1] + normFld_n[3*dir+2]*curlbhat_n[2];

      }
    }
  }
}

void
gkyl_calc_metric_mirror_release(gkyl_calc_metric_mirror* up)
{
  gkyl_nodal_ops_release(up->n2m);
  gkyl_free(up);
}
