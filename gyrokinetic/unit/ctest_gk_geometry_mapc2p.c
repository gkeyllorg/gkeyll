#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>


#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_ops.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <gkyl_basis.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_priv.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_nodal_ops.h>
#include <gkyl_position_map.h>

#include <gkyl_comm.h>

void
write_geometry(gk_geometry *up, struct gkyl_rect_grid grid, struct gkyl_range local, const char *name)
{
}


// Functions for this test
void mapc2p(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double r = xn[0], theta = xn[1], z = xn[2];
  fout[0] = r*cos(theta); 
  fout[1] = r*sin(theta); 
  fout[2] = z;
}

void exact_gij(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double r = xn[0], theta = xn[1], phi = xn[2];
  fout[0] = 1.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = r*r;
  fout[4] = 0.0;
  fout[5] = 1.0;
}

void bfield_func(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx){
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0398;
}

void
test_3x_p1()
{
  struct gkyl_basis basis;
  int poly_order = 1;
  gkyl_cart_modal_serendip(&basis, 3, poly_order);
  
  double Lz = 1.8049e+01;
  double Lx = 1.2534e+00;
  double Rmax = Lx/2;
  double Rmin = 0.5*Rmax;
  int Nz = 10;

  double lower[3] = {Rmin, -0.1, -Lz/2};
  double upper[3] = {Rmax,  0.1,  Lz/2};
  int cells[3] = { 8, 1, 8 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 3, lower, upper, cells);
  
  struct gkyl_range ext_range, range;
  int nghost[3] = { 1,1,1};
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  struct gkyl_position_map *pmap = gkyl_position_map_null_new();

  // Initialize geometry
  struct gkyl_gk_geometry_inp geometry_input = {
    .geometry_id = GKYL_GEOMETRY_MAPC2P,
    .mapc2p = mapc2p, // mapping of computational to physical space
    .c2p_ctx = 0,
    .bfield_func = bfield_func, // magnetic field magnitude
    .bfield_ctx =0 ,
    .position_map = pmap,
    .grid = grid,
    .local = range,
    .local_ext = ext_range,
    .global = range,
    .global_ext = ext_range,
    .basis = basis,
    .geo_grid = grid,
    .geo_local = range,
    .geo_local_ext = ext_range,
    .geo_global = range,
    .geo_global_ext = ext_range,
    .geo_basis = basis,
  };

  struct gk_geometry *gk_geom = gkyl_gk_geometry_mapc2p_new(&geometry_input);
  //write_geometry(gk_geom, grid, range, "geomapc2p");

  // Define nodal operations
  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  int cidx[3];
  int nodes[] = { 1, 1, 1 };
  for (int d=0; d<grid.ndim; ++d)
    nodes[d] = grid.cells[d] + 1;
  struct gkyl_range nrange;
  gkyl_range_init_from_shape(&nrange, grid.ndim, nodes);

  int nodes_quad_interior[] = { 1, 1, 1 };
  int num_quad_points=poly_order+1;
  for (int d=0; d<grid.ndim; ++d)
    nodes_quad_interior[d] = grid.cells[d]*num_quad_points;
  struct gkyl_range nrange_quad_interior;
  gkyl_range_init_from_shape(&nrange_quad_interior, grid.ndim, nodes_quad_interior);

  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&basis, &grid, false);

  struct gkyl_array* mapc2p_nodal = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange, &range, 3, mapc2p_nodal, gk_geom->geo_corn.mc2p, false);

  struct gkyl_array* mapc2p_nodal_interior = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange_quad_interior, &range, 3, mapc2p_nodal_interior, gk_geom->geo_int.mc2p, true);

  // Check that |bhat|=1 at nodes
  struct gkyl_array* bhat_nodal = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange_quad_interior, &range, 3, bhat_nodal, gk_geom->geo_int.bcart, true);
  for (int ia=nrange_quad_interior.lower[AL_IDX]; ia<=nrange_quad_interior.upper[AL_IDX]; ++ia){
    for (int ip=nrange_quad_interior.lower[PSI_IDX]; ip<=nrange_quad_interior.upper[PSI_IDX]; ++ip) {
      for (int it=nrange_quad_interior.lower[TH_IDX]; it<=nrange_quad_interior.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *bhat_n = gkyl_array_fetch(bhat_nodal, gkyl_range_idx(&nrange_quad_interior, cidx));
        double bhat_mag = sqrt(bhat_n[0]*bhat_n[0] + bhat_n[1]*bhat_n[1] + bhat_n[2]*bhat_n[2]);
        TEST_CHECK( gkyl_compare( bhat_mag, 1.0, 1e-12) );
      }
    }
  }

  // Check that the duals are what they should be
  struct gkyl_array* dualmag_nodal = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange_quad_interior, &range, 3, dualmag_nodal, gk_geom->geo_int.dualmag, true);
  for (int ia=nrange_quad_interior.lower[AL_IDX]; ia<=nrange_quad_interior.upper[AL_IDX]; ++ia){
    for (int ip=nrange_quad_interior.lower[PSI_IDX]; ip<=nrange_quad_interior.upper[PSI_IDX]; ++ip) {
      for (int it=nrange_quad_interior.lower[TH_IDX]; it<=nrange_quad_interior.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *dualmag_n = gkyl_array_fetch(dualmag_nodal, gkyl_range_idx(&nrange_quad_interior, cidx));
        double *mapc2p_n = gkyl_array_fetch(mapc2p_nodal_interior, gkyl_range_idx(&nrange_quad_interior, cidx));
        double e2mag = sqrt(1/(mapc2p_n[0]*mapc2p_n[0] +  mapc2p_n[1]*mapc2p_n[1])); // 1/R
        TEST_CHECK( gkyl_compare( dualmag_n[0], 1.0, 1e-8) );
        TEST_CHECK( gkyl_compare( dualmag_n[1], e2mag, 1e-6) );
        TEST_CHECK( gkyl_compare( dualmag_n[2], 1.0, 1e-8) );
      }
    }
  }


  // Check that Jacobgeo is what it should be. J = R in cylindrical coordinates
  struct gkyl_array* jacobgeo_nodal = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange_quad_interior, &range, 1, jacobgeo_nodal, gk_geom->geo_int.jacobgeo, true);
  for (int ia=nrange_quad_interior.lower[AL_IDX]; ia<=nrange_quad_interior.upper[AL_IDX]; ++ia){
    for (int ip=nrange_quad_interior.lower[PSI_IDX]; ip<=nrange_quad_interior.upper[PSI_IDX]; ++ip) {
      for (int it=nrange_quad_interior.lower[TH_IDX]; it<=nrange_quad_interior.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *jacobgeo_n = gkyl_array_fetch(jacobgeo_nodal, gkyl_range_idx(&nrange_quad_interior, cidx));
        // mapc2p_n[0] = x, mapc2p_n[1] = y, mapc2p_n[2] = z
        double *mapc2p_n = gkyl_array_fetch(mapc2p_nodal_interior, gkyl_range_idx(&nrange_quad_interior, cidx));
        double radius = sqrt(mapc2p_n[0]*mapc2p_n[0] + mapc2p_n[1]*mapc2p_n[1]);
        TEST_CHECK( gkyl_compare( jacobgeo_n[0], radius, 1e-6) );
      }
    }
  }

  // Check bmag is what it should be
  struct gkyl_array* bmag_nodal = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange_quad_interior, &range, 1, bmag_nodal, gk_geom->geo_int.bmag, true);
  for (int ia=nrange_quad_interior.lower[AL_IDX]; ia<=nrange_quad_interior.upper[AL_IDX]; ++ia){
    for (int ip=nrange_quad_interior.lower[PSI_IDX]; ip<=nrange_quad_interior.upper[PSI_IDX]; ++ip) {
      for (int it=nrange_quad_interior.lower[TH_IDX]; it<=nrange_quad_interior.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double psi = grid.lower[PSI_IDX] + ip*(grid.upper[PSI_IDX]-grid.lower[PSI_IDX])/grid.cells[PSI_IDX];
        double alpha = grid.lower[AL_IDX] + ia*(grid.upper[AL_IDX]-grid.lower[AL_IDX])/grid.cells[AL_IDX];
        double theta = grid.lower[TH_IDX] + it*(grid.upper[TH_IDX]-grid.lower[TH_IDX])/grid.cells[TH_IDX];
        double xn[3] = {psi, alpha, theta};
        double *bmag_n = gkyl_array_fetch(bmag_nodal, gkyl_range_idx(&nrange_quad_interior, cidx));
        double bmag_anal[3];
        bfield_func(0, xn, bmag_anal, 0);
        TEST_CHECK( gkyl_compare( bmag_n[0], bmag_anal[2], 1e-8) );
      }
    }
  }

  // Check gij
  struct gkyl_array* gij_nodal = gkyl_array_new(GKYL_DOUBLE, 6, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange_quad_interior, &range, 6, gij_nodal, gk_geom->geo_int.g_ij, true);
  for (int ia=nrange_quad_interior.lower[AL_IDX]; ia<=nrange_quad_interior.upper[AL_IDX]; ++ia){
    for (int ip=nrange_quad_interior.lower[PSI_IDX]; ip<=nrange_quad_interior.upper[PSI_IDX]; ++ip) {
      for (int it=nrange_quad_interior.lower[TH_IDX]; it<=nrange_quad_interior.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *gij_n = gkyl_array_fetch(gij_nodal, gkyl_range_idx(&nrange_quad_interior, cidx));
        double *mapc2p_n = gkyl_array_fetch(mapc2p_nodal_interior, gkyl_range_idx(&nrange_quad_interior, cidx));
        double r = sqrt(mapc2p_n[0]*mapc2p_n[0] + mapc2p_n[1]*mapc2p_n[1]);
        double xn[3] = {r, 0.0, 0.0};
        double fout[6];
        exact_gij(0.0, xn, fout, 0);
        for (int i=0; i<6; ++i)
          TEST_CHECK( gkyl_compare( gij_n[i], fout[i], 1e-6) );
      }
    }
  }

  // Check mapc2p
  for (int ia=nrange.lower[AL_IDX]; ia<=nrange.upper[AL_IDX]; ++ia){
    for (int ip=nrange.lower[PSI_IDX]; ip<=nrange.upper[PSI_IDX]; ++ip) {
      for (int it=nrange.lower[TH_IDX]; it<=nrange.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double psi = grid.lower[PSI_IDX] + ip*(grid.upper[PSI_IDX]-grid.lower[PSI_IDX])/grid.cells[PSI_IDX];
        double alpha = grid.lower[AL_IDX] + ia*(grid.upper[AL_IDX]-grid.lower[AL_IDX])/grid.cells[AL_IDX];
        double theta = grid.lower[TH_IDX] + it*(grid.upper[TH_IDX]-grid.lower[TH_IDX])/grid.cells[TH_IDX]; 
        // mapc2p_n[0] = x, mapc2p_n[1] = y, mapc2p_n[2] = z
        double *mapc2p_n = gkyl_array_fetch(mapc2p_nodal, gkyl_range_idx(&nrange, cidx));
        double xn[3] = {psi, alpha, theta};
        double fout[3];
        mapc2p(0.0, xn, fout, 0);
        for (int i=0; i<3; ++i)
          TEST_CHECK( gkyl_compare( mapc2p_n[i], fout[i], 1e-6) );
      }
    }
  }
  
  // Release memory

  gkyl_array_release(bhat_nodal);
  gkyl_array_release(dualmag_nodal);
  gkyl_array_release(mapc2p_nodal);
  gkyl_array_release(mapc2p_nodal_interior);
  gkyl_array_release(jacobgeo_nodal);
  gkyl_array_release(bmag_nodal);
  gkyl_array_release(gij_nodal);
  gkyl_nodal_ops_release(n2m);
  gkyl_position_map_release(pmap);
  gkyl_gk_geometry_release(gk_geom);
}

TEST_LIST = {
  { "test_3x_p1", test_3x_p1},
  { NULL, NULL },
};
