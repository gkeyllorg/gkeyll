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
  const char *fmt = "%s-%s.gkyl";
  int sz = gkyl_calc_strlen(fmt, name, "jacobtot_inv");
  char fileNm[sz+1]; // ensure no buffer overflow

  sprintf(fileNm, fmt, name, "mapc2p");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->geo_corn.mc2p, fileNm);
  sprintf(fileNm, fmt, name, "mapc2nu");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->geo_corn.mc2nu_pos, fileNm);
  sprintf(fileNm, fmt, name, "bmag_corn");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->geo_corn.bmag, fileNm);
  sprintf(fileNm, fmt, name, "bmag");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->geo_int.bmag, fileNm);
  sprintf(fileNm, fmt, name, "g_ij");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->geo_int.g_ij, fileNm);
  sprintf(fileNm, fmt, name, "dxdz");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->geo_int.dxdz, fileNm);
  sprintf(fileNm, fmt, name, "dzdx");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->geo_int.dzdx, fileNm);
  sprintf(fileNm, fmt, name, "normals");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->geo_int.normals, fileNm);
  sprintf(fileNm, fmt, name, "jacobgeo");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->geo_int.jacobgeo, fileNm);
  sprintf(fileNm, fmt, name, "jacobgeo_inv");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->geo_int.jacobgeo_inv, fileNm);
  sprintf(fileNm, fmt, name, "gij");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->geo_int.gij, fileNm);
  sprintf(fileNm, fmt, name, "b_i");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->geo_int.b_i, fileNm);
  sprintf(fileNm, fmt, name, "bcart");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->geo_int.bcart, fileNm);
  sprintf(fileNm, fmt, name, "cmag");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->geo_int.cmag, fileNm);
  sprintf(fileNm, fmt, name, "jacobtot");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->geo_int.jacobtot, fileNm);
  sprintf(fileNm, fmt, name, "jacobtot_inv");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->geo_int.jacobtot_inv, fileNm);
  sprintf(fileNm, fmt, name, "gxxj");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->geo_int.gxxj, fileNm);
  sprintf(fileNm, fmt, name, "gxyj");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->geo_int.gxyj, fileNm);
  sprintf(fileNm, fmt, name, "gyyj");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->geo_int.gyyj, fileNm);
  sprintf(fileNm, fmt, name, "gxzj");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->geo_int.gxzj, fileNm);
  sprintf(fileNm, fmt, name, "eps2");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->geo_int.eps2, fileNm);


  // Create Nodal Range and Grid and Write Nodal Coordinates
  struct gkyl_range nrange;
  gkyl_gk_geometry_init_nodal_range(&nrange, &local, up->basis.poly_order);
  struct gkyl_array* mc2p_nodal = gkyl_array_new(GKYL_DOUBLE, 3, nrange.volume);
  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&up->basis, &grid, false);
  gkyl_nodal_ops_m2n(n2m, &up->basis, &grid, &nrange, &local, 3, mc2p_nodal, up->geo_corn.mc2p, false);
  gkyl_nodal_ops_release(n2m);
  struct gkyl_rect_grid ngrid;
  gkyl_gk_geometry_init_nodal_grid(&ngrid, &grid, &nrange);
  sprintf(fileNm, fmt, name, "nodes");
  gkyl_grid_sub_array_write(&ngrid, &nrange, 0,  mc2p_nodal, fileNm);
  gkyl_array_release(mc2p_nodal);
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

void
mapz(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double Lz = 1.8049e+01;
  double a = -Lz/2;
  fout[0] = -1/(2*a) * pow(a - xn[0], 2) + a;
}

void
dmapz_dz(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double Lz = 1.8049e+01;
  double a = -Lz/2;
  fout[0] = 1 - xn[0]/a;
}

void
test_3x_p1_pmap()
{  
  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  struct gkyl_basis basis;
  int poly_order = 1;
  int cdim = 3;
  gkyl_cart_modal_serendip(&basis, cdim, poly_order);
  
  double Lz = 1.8049e+01;
  double Lx = 1.2534e+00;
  double Rmax = Lx/2;
  double Rmin = 0.05*Rmax;
  int Nz = 10;

  double lower[3] = {Rmin, -M_PI, -Lz/2};
  double upper[3] = {Rmax,  M_PI,  Lz/2};
  int cells[3] = { 18, 18, Nz };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);
  
  struct gkyl_range ext_range, range;
  int nghost[3] = { 1,1,1};
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  struct gkyl_position_map_inp pos_map_inp = {
    .maps = {0, 0, mapz},
    .ctxs = {0, 0, 0},
  };

  // Configuration space geometry initialization
  struct gkyl_position_map *pos_map = gkyl_position_map_new(pos_map_inp, grid, range, 
    ext_range, range, ext_range, basis);

  // Initialize geometry
  struct gkyl_gk_geometry_inp geometry_input = {
    .geometry_id = GKYL_GEOMETRY_MAPC2P,
    .mapc2p = mapc2p, // mapping of computational to physical space
    .c2p_ctx = 0,
    .bfield_func = bfield_func, // magnetic field magnitude
    .bfield_ctx =0 ,
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
    .position_map = pos_map,
  };

  struct gk_geometry *gk_geom = gkyl_gk_geometry_mapc2p_new(&geometry_input);
  
  gkyl_position_map_set_mc2nu(pos_map, gk_geom->geo_corn.mc2nu_pos);

  // Define the nodes for the script to calculate values at
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


  double dels[2] = {1.0/sqrt(3), 1.0-1.0/sqrt(3) };
  double theta_lo = grid.lower[TH_IDX] + dels[1]*grid.dx[TH_IDX]/2.0;
  double psi_lo = grid.lower[PSI_IDX] + dels[1]*grid.dx[PSI_IDX]/2.0;
  double alpha_lo = grid.lower[AL_IDX] + dels[1]*grid.dx[AL_IDX]/2.0;

  struct gkyl_array* mapc2p_nodal = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange, &range, 3, mapc2p_nodal, gk_geom->geo_corn.mc2p, false);

  struct gkyl_array* mapc2p_nodal_interior = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange_quad_interior, &range, 3, mapc2p_nodal_interior, gk_geom->geo_int.mc2p, true);

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
        mapz(0.0, &theta, &theta, 0);
        double xn[3] = {psi, alpha, theta};
        double fout[3];
        mapc2p(0.0, xn, fout, 0);
        for (int i=0; i<3; ++i){
          TEST_CHECK( gkyl_compare( mapc2p_n[i], fout[i], 1e-8) );
          TEST_MSG( "i %d, idx %d %d %d, mapc2p_n %g, fout %g\n", i, ip, ia, it, mapc2p_n[i], fout[i] );
        }
      }
    }
  }

  // Check that Jacobgeo is what it should be. J = R in cylindrical coordinates
  // We have a contribution from the position map too, given as dZ/dz
  struct gkyl_array* jacobgeo_nodal = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange_quad_interior, &range, 1, jacobgeo_nodal, gk_geom->geo_int.jacobgeo, true);
  for (int ia=nrange_quad_interior.lower[AL_IDX]; ia<=nrange_quad_interior.upper[AL_IDX]; ++ia){
    for (int ip=nrange_quad_interior.lower[PSI_IDX]; ip<=nrange_quad_interior.upper[PSI_IDX]; ++ip) {
      for (int it=nrange_quad_interior.lower[TH_IDX]; it<=nrange_quad_interior.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double psi = calc_running_coord(psi_lo, ip-nrange_quad_interior.lower[PSI_IDX], grid.dx[PSI_IDX]);
        double alpha = calc_running_coord(alpha_lo, ia-nrange_quad_interior.lower[AL_IDX], grid.dx[AL_IDX]);
        double theta = calc_running_coord(theta_lo, it-nrange_quad_interior.lower[TH_IDX], grid.dx[TH_IDX]);
        double *jacobgeo_n = gkyl_array_fetch(jacobgeo_nodal, gkyl_range_idx(&nrange_quad_interior, cidx));
        // mapc2p_n[0] = x, mapc2p_n[1] = y, mapc2p_n[2] = z
        double *mapc2p_n = gkyl_array_fetch(mapc2p_nodal_interior, gkyl_range_idx(&nrange_quad_interior, cidx));
        double radius = sqrt(mapc2p_n[0]*mapc2p_n[0] + mapc2p_n[1]*mapc2p_n[1]);
        double fout[1];
        dmapz_dz(0.0, &theta, fout, 0);
        double jacob_anal = radius * fout[0];
        TEST_CHECK( gkyl_compare( jacobgeo_n[0], jacob_anal, 1e-5) );
      }
    }
  }

  // Check bmag is what it should be
  struct gkyl_array* bmag_nodal = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange, &range, 1, bmag_nodal, gk_geom->geo_corn.bmag, false);
  for (int ia=nrange.lower[AL_IDX]; ia<=nrange.upper[AL_IDX]; ++ia){
    for (int ip=nrange.lower[PSI_IDX]; ip<=nrange.upper[PSI_IDX]; ++ip) {
      for (int it=nrange.lower[TH_IDX]; it<=nrange.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double psi = grid.lower[PSI_IDX] + ip*(grid.upper[PSI_IDX]-grid.lower[PSI_IDX])/grid.cells[PSI_IDX];
        double alpha = grid.lower[AL_IDX] + ia*(grid.upper[AL_IDX]-grid.lower[AL_IDX])/grid.cells[AL_IDX];
        double theta = grid.lower[TH_IDX] + it*(grid.upper[TH_IDX]-grid.lower[TH_IDX])/grid.cells[TH_IDX];
        mapz(0.0, &theta, &theta, 0);
        double xn[3] = {psi, alpha, theta};
        double *bmag_n = gkyl_array_fetch(bmag_nodal, gkyl_range_idx(&nrange, cidx));
        double bmag_anal[3];
        bfield_func(0, xn, bmag_anal, 0);
        TEST_CHECK( gkyl_compare( bmag_n[0], bmag_anal[2], 1e-8) );
      }
    }
  }


  // Release memory
  gkyl_array_release(jacobgeo_nodal);
  gkyl_array_release(mapc2p_nodal);
  gkyl_array_release(mapc2p_nodal_interior);
  gkyl_array_release(bmag_nodal);
  gkyl_nodal_ops_release(n2m);
  gkyl_position_map_release(pos_map);
  gkyl_gk_geometry_release(gk_geom);
}

TEST_LIST = {
  { "test_3x_p1", test_3x_p1},
  { "test_3x_p1_pmap", test_3x_p1_pmap},
  { NULL, NULL },
};
