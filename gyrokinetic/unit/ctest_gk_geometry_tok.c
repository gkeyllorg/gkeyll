#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_calc_bmag.h>
#include <gkyl_calc_derived_geo.h>
#include <gkyl_calc_metric.h>
#include <gkyl_efit.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_priv.h>
#include <gkyl_gk_geometry_tok.h>
#include <gkyl_nodal_ops.h>
#include <gkyl_position_map.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_tok_geo.h>
#include <gkyl_util.h>

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
  sprintf(fileNm, fmt, name, "bmag_inv_corn");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->geo_corn.bmag_inv, fileNm);
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
  sprintf(fileNm, fmt, name, "qprofile");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->geo_int.qprofile, fileNm);

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

void
test_tok_elliptical_ho()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  struct gkyl_efit_inp efit_inp = {
      // psiRZ and related inputs
      .filepath = "gyrokinetic/data/eqdsk/elliptical.geqdsk",
      .rz_poly_order = 2,
      .flux_poly_order = 1,
      .reflect = true,
    };

  double psisep = -4.0;
  double clower[] = { -5.0, -0.01, -M_PI+1e-14 };
  double cupper[] = {psisep, 0.01, M_PI-1e-14 };

  int ccells[] = { 2, 1, 2 };

  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  int cpoly_order = 1;
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);
  struct gkyl_position_map *pmap = gkyl_position_map_null_new();

  struct gkyl_tok_geo_grid_inp ginp = {
    .rmin = 0.0,
    .rmax = 5.0,
    .ftype = GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT,
    .rclose = 6.0,
    .rright = 6.0,
    .rleft = 0.0,
    .zmin = -3.0,
    .zmax = 3.0,
  }; 

  struct gkyl_gk_geometry_inp geometry_inp = {
    .geometry_id  = GKYL_GEOMETRY_TOKAMAK,
    .efit_info = efit_inp,
    .tok_grid_info = ginp,
    .position_map = pmap,
    .grid = cgrid,
    .local = clocal,
    .local_ext = clocal_ext,
    .global = clocal,
    .global_ext = clocal_ext,
    .basis = cbasis,
    .geo_grid = cgrid,
    .geo_local = clocal,
    .geo_local_ext = clocal_ext,
    .geo_global = clocal,
    .geo_global_ext = clocal_ext,
    .geo_basis = cbasis,
  };

  struct gk_geometry* up = gkyl_gk_geometry_tok_new(&geometry_inp); 
  gkyl_gk_geometry_release(up);
  gkyl_position_map_release(pmap);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
}


// Functions for test_3x_straight_cylinder
void mapc2p(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double psi = xn[0], alpha = xn[1], theta = xn[2];
  fout[0] = sqrt(psi * 4 ); // Function fed is psi = 0.5/2 * R^2 from the efit file
  fout[1] = theta * 1.0 / M_PI; // Note that this does not have pi-1e-2 in it because the coordinate zeta is always defined -pi to pi
  fout[2] = alpha; // There is a minus due to conventions
}

void exact_gij(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double r = xn[0], theta = xn[1], phi = xn[2];
  double psi = r*r/4;
  fout[0] = 1/psi; // g_11
  fout[1] = 0.0; // g_12
  fout[2] = 0.0; // g_13
  fout[3] = r*r; // g_22
  fout[4] = 0.0; // g_23
  fout[5] = 1/(M_PI*M_PI); // g_33
}

void exact_g_contra_ij(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double r = xn[0], theta = xn[1], phi = xn[2];
  double psi = r*r/4;
  fout[0] = r*r/4; // g_11
  fout[1] = 0.0; // g_12
  fout[2] = 0.0; // g_13
  fout[3] = 1/psi/4; // g_22
  fout[4] = 0.0; // g_23
  fout[5] = (M_PI*M_PI); // g_33
}

void exact_dual_magnitude(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double r = xn[0], theta = xn[1], phi = xn[2];
  double psi = r*r/4;
  fout[0]  = r/2;
  fout[1] = 1/(2*sqrt(psi));
  fout[2] = M_PI;
}

void exact_normals(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double psi = xn[0], alpha = xn[1], theta = xn[2];
  // Remember cylindrical angle = - alpha
  fout[0] = cos(alpha);
  fout[1] = sin(alpha);
  fout[2] = 0.0;
  fout[3] = -sin(alpha);
  fout[4] = cos(alpha);
  fout[5] = 0.0;
  fout[6] = 0.0;
  fout[7] = 0.0;
  fout[8] = 1.0;
}

void bmag_func(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx){
  fout[0] = 0.5;
}

void
test_tok_3x_p1_straight_cylinder_ho()
{
  // Very similar to the unit test in ctest_gk_geometry.c
  // The geometry is created to extend from Z = -1 to 1, R = (0.001, 1) in units meters
  // Magnetic field = 0.5 uniform everywhere
  // Psi = 0.5/2 * R^2
  // Check write_efit_straight_cylinder.py for how the efit file is written
  // There are some important differences between how the numerical geometry is calculated compared to mapc2p geometries

  struct gkyl_basis basis;
  int poly_order = 1;
  int cdim = 3;
  gkyl_cart_modal_serendip(&basis, cdim, poly_order);
  
  double psiMax = 0.2;
  double psiMin = 0.1;

  double lower[3] = {psiMin, -1.0, -M_PI+1e-14};
  double upper[3] = {psiMax,  1.0,  M_PI-1e-14};
  int cells[3] = { 8, 2, 8};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);
  
  struct gkyl_range ext_range, range;
  int nghost[3] = { 1,1,1};
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  struct gkyl_position_map *pmap = gkyl_position_map_null_new();

  struct gkyl_efit_inp inp = {
    // psiRZ and related inputs
    .filepath = "gyrokinetic/data/eqdsk/straight_cylinder.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect = true,
  };
  struct gkyl_tok_geo_grid_inp ginp = {
    .rclose = 0.5,
    .zmin = -1.,
    .zmax =  1.,
    .rleft = 0.001,
    .rmax = 1.0,
    .rright = 1.0,
  };
  // Initialize geometry
  struct gkyl_gk_geometry_inp geometry_input = {
    .geometry_id = GKYL_GEOMETRY_TOKAMAK,
    .efit_info = inp,
    .tok_grid_info = ginp,
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

  struct gk_geometry *gk_geom = gkyl_gk_geometry_tok_new(&geometry_input);

  //write_geometry(gk_geom, grid, range, "straight_cylinder");

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
  // There are errors at low psi, so we shift away from the axis by a few cells
  struct gkyl_array* dualmag_nodal = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange_quad_interior, &range, 3, dualmag_nodal, gk_geom->geo_int.dualmag, true);
  struct gkyl_array* mapc2p_nodal_interior = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange_quad_interior, &range, 3, mapc2p_nodal_interior, gk_geom->geo_int.mc2p, true);
  for (int ia=nrange_quad_interior.lower[AL_IDX]; ia<=nrange_quad_interior.upper[AL_IDX]; ++ia){
    for (int ip=nrange_quad_interior.lower[PSI_IDX]; ip<=nrange_quad_interior.upper[PSI_IDX]; ++ip) {
      for (int it=nrange_quad_interior.lower[TH_IDX]; it<=nrange_quad_interior.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *dualmag_n = gkyl_array_fetch(dualmag_nodal, gkyl_range_idx(&nrange_quad_interior, cidx));
        double *mapc2p_n = gkyl_array_fetch(mapc2p_nodal_interior, gkyl_range_idx(&nrange_quad_interior, cidx));
        double xn[3] = {mapc2p_n[0], mapc2p_n[1], mapc2p_n[2]};
        double dualmag_ana[3];
        exact_dual_magnitude(0, xn, dualmag_ana, 0);
        TEST_CHECK( gkyl_compare( dualmag_n[0], dualmag_ana[0], 1e-6) );
        TEST_CHECK( gkyl_compare( dualmag_n[1], dualmag_ana[1], 1e-6) );
        TEST_CHECK( gkyl_compare( dualmag_n[2], dualmag_ana[2], 1e-6) );
      }
    }
  }

  // Check bmag is what it should be
  struct gkyl_array* bmag_nodal = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange, &range, 1, bmag_nodal, gk_geom->geo_corn.bmag, false);
  struct gkyl_array* mapc2p_nodal = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange, &range, 3, mapc2p_nodal, gk_geom->geo_corn.mc2p, false);
  for (int ia=nrange.lower[AL_IDX]; ia<=nrange.upper[AL_IDX]; ++ia){
    for (int ip=nrange.lower[PSI_IDX]; ip<=nrange.upper[PSI_IDX]; ++ip) {
      for (int it=nrange.lower[TH_IDX]; it<=nrange.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double psi = grid.lower[PSI_IDX] + ip*(grid.upper[PSI_IDX]-grid.lower[PSI_IDX])/grid.cells[PSI_IDX];
        double alpha = grid.lower[AL_IDX] + ia*(grid.upper[AL_IDX]-grid.lower[AL_IDX])/grid.cells[AL_IDX];
        double theta = grid.lower[TH_IDX] + it*(grid.upper[TH_IDX]-grid.lower[TH_IDX])/grid.cells[TH_IDX];
        double xn[3] = {psi, alpha, theta};
        double *bmag_n = gkyl_array_fetch(bmag_nodal, gkyl_range_idx(&nrange, cidx));
        double bmag_ana[1];
        bmag_func(0, xn, bmag_ana, 0);
        TEST_CHECK( gkyl_compare( bmag_n[0], bmag_ana[0], 1e-8) );
      }
    }
  }

  // Check cmag = 1
  struct gkyl_array* cmag_nodal = gkyl_array_new(GKYL_DOUBLE, 1, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange_quad_interior, &range, 1, cmag_nodal, gk_geom->geo_int.cmag, true);
  for (int ia=nrange_quad_interior.lower[AL_IDX]; ia<=nrange_quad_interior.upper[AL_IDX]; ++ia){
    for (int ip=nrange_quad_interior.lower[PSI_IDX]; ip<=nrange_quad_interior.upper[PSI_IDX]; ++ip) {
      for (int it=nrange_quad_interior.lower[TH_IDX]; it<=nrange_quad_interior.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *cmag_n = gkyl_array_fetch(cmag_nodal, gkyl_range_idx(&nrange_quad_interior, cidx));
        TEST_CHECK( gkyl_compare( cmag_n[0], 1.0, 1e-6) );
      }
    }
  }

  // Check g_ij
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
        double r = mapc2p_n[0];
        double xn[3] = {r, 0.0, 0.0};
        double fout[6];
        exact_gij(0.0, xn, fout, 0);
        double tol;
        for (int i=0; i<6; ++i)
        {
          if (i == 4)
            tol = 1e-3;
          else
            tol = 1e-6;
          TEST_CHECK( gkyl_compare( gij_n[i], fout[i], tol) );
        }
      }
    }
  }

  // Check g^ij
  struct gkyl_array* gij_contra_nodal = gkyl_array_new(GKYL_DOUBLE, 6, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange_quad_interior, &range, 6, gij_contra_nodal, gk_geom->geo_int.gij, true);
  for (int ia=nrange_quad_interior.lower[AL_IDX]; ia<=nrange_quad_interior.upper[AL_IDX]; ++ia){
    for (int ip=nrange_quad_interior.lower[PSI_IDX]; ip<=nrange_quad_interior.upper[PSI_IDX]; ++ip) {
      for (int it=nrange_quad_interior.lower[TH_IDX]; it<=nrange_quad_interior.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *gij_contra_n = gkyl_array_fetch(gij_contra_nodal, gkyl_range_idx(&nrange_quad_interior, cidx));
        double *mapc2p_n = gkyl_array_fetch(mapc2p_nodal_interior, gkyl_range_idx(&nrange_quad_interior, cidx));
        double r = mapc2p_n[0];
        double xn[3] = {r, 0.0, 0.0};
        double fout[6];
        exact_g_contra_ij(0.0, xn, fout, 0);
        double tol;
        for (int i=0; i<6; ++i)
        {
          if (i == 4)
            tol = 1e-2;
          else
            tol = 1e-6;
          TEST_CHECK( gkyl_compare( gij_contra_n[i], fout[i], tol) );
        }
      }
    }
  }

  // Check that Jacobgeo is what it should be. This is the Jacobian for the problem
  struct gkyl_array* jacobgeo_nodal = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange_quad_interior, &range, 1, jacobgeo_nodal, gk_geom->geo_int.jacobgeo, true);
  for (int ia=nrange_quad_interior.lower[AL_IDX]; ia<=nrange_quad_interior.upper[AL_IDX]; ++ia){
    for (int ip=nrange_quad_interior.lower[PSI_IDX]; ip<=nrange_quad_interior.upper[PSI_IDX]; ++ip) {
      for (int it=nrange_quad_interior.lower[TH_IDX]; it<=nrange_quad_interior.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *jacobgeo_n = gkyl_array_fetch(jacobgeo_nodal, gkyl_range_idx(&nrange_quad_interior, cidx));
        double *mapc2p_n = gkyl_array_fetch(mapc2p_nodal_interior, gkyl_range_idx(&nrange_quad_interior, cidx));
        double jacobian_anaytic = 2/M_PI;
        TEST_CHECK( gkyl_compare( jacobgeo_n[0], jacobian_anaytic, 1e-6) );
      }
    }
  }

  // Check jacobgeo_inv
  struct gkyl_array* jacobgeo_inv_nodal = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange_quad_interior, &range, 1, jacobgeo_inv_nodal, gk_geom->geo_int.jacobgeo_inv, true);
  for (int ia=nrange_quad_interior.lower[AL_IDX]; ia<=nrange_quad_interior.upper[AL_IDX]; ++ia){
    for (int ip=nrange_quad_interior.lower[PSI_IDX]; ip<=nrange_quad_interior.upper[PSI_IDX]; ++ip) {
      for (int it=nrange_quad_interior.lower[TH_IDX]; it<=nrange_quad_interior.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *jacobgeo_inv_n = gkyl_array_fetch(jacobgeo_inv_nodal, gkyl_range_idx(&nrange_quad_interior, cidx));
        double *mapc2p_n = gkyl_array_fetch(mapc2p_nodal_interior, gkyl_range_idx(&nrange_quad_interior, cidx));
        double jacobian_anaytic = 2/M_PI;
        TEST_CHECK( gkyl_compare( jacobgeo_inv_n[0], 1/jacobian_anaytic, 1e-6) );
      }
    }
  }

  // Check jacobtot
  struct gkyl_array* jacobtot_nodal = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange_quad_interior, &range, 1, jacobtot_nodal, gk_geom->geo_int.jacobtot, true);
  for (int ia=nrange_quad_interior.lower[AL_IDX]; ia<=nrange_quad_interior.upper[AL_IDX]; ++ia){
    for (int ip=nrange_quad_interior.lower[PSI_IDX]; ip<=nrange_quad_interior.upper[PSI_IDX]; ++ip) {
      for (int it=nrange_quad_interior.lower[TH_IDX]; it<=nrange_quad_interior.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *jacobtot_n = gkyl_array_fetch(jacobtot_nodal, gkyl_range_idx(&nrange_quad_interior, cidx));
        // mapc2p_n[0] = x, mapc2p_n[1] = y, mapc2p_n[2] = z
        double *mapc2p_n = gkyl_array_fetch(mapc2p_nodal_interior, gkyl_range_idx(&nrange_quad_interior, cidx));
        double jacobian_anaytic = 2/M_PI;
        double magnetic_field = 0.5;
        double jacobtot_anaytic = jacobian_anaytic * magnetic_field;
        TEST_CHECK( gkyl_compare( jacobtot_n[0], jacobtot_anaytic, 1e-6) );
      }
    }
  }

  // Check jacobtot_inv
  struct gkyl_array* jacobtot_inv_nodal = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange_quad_interior, &range, 1, jacobtot_inv_nodal, gk_geom->geo_int.jacobtot_inv, true);
  for (int ia=nrange_quad_interior.lower[AL_IDX]; ia<=nrange_quad_interior.upper[AL_IDX]; ++ia){
    for (int ip=nrange_quad_interior.lower[PSI_IDX]; ip<=nrange_quad_interior.upper[PSI_IDX]; ++ip) {
      for (int it=nrange_quad_interior.lower[TH_IDX]; it<=nrange_quad_interior.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *jacobtot_inv_n = gkyl_array_fetch(jacobtot_inv_nodal, gkyl_range_idx(&nrange_quad_interior, cidx));
        // mapc2p_n[0] = x, mapc2p_n[1] = y, mapc2p_n[2] = z
        double *mapc2p_n = gkyl_array_fetch(mapc2p_nodal_interior, gkyl_range_idx(&nrange_quad_interior, cidx));
        double jacobian_anaytic = 2/M_PI;
        double magnetic_field = 0.5;
        double jacobtot_anaytic = jacobian_anaytic * magnetic_field;
        TEST_CHECK( gkyl_compare( jacobtot_inv_n[0], 1/jacobtot_anaytic, 1e-6) );
      }
    }
  }

  // Check mapc2p
  for (int ia=nrange.lower[AL_IDX]; ia<=nrange.upper[AL_IDX]; ++ia){
    for (int ip=nrange.lower[PSI_IDX]; ip<=nrange.upper[PSI_IDX]; ++ip) {
      for (int it=nrange.lower[TH_IDX]+1; it<=nrange.upper[TH_IDX]-1; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double psi = grid.lower[PSI_IDX] + ip*(grid.upper[PSI_IDX]-grid.lower[PSI_IDX])/grid.cells[PSI_IDX];
        double alpha = grid.lower[AL_IDX] + ia*(grid.upper[AL_IDX]-grid.lower[AL_IDX])/grid.cells[AL_IDX];
        double theta = grid.lower[TH_IDX] + it*(grid.upper[TH_IDX]-grid.lower[TH_IDX])/grid.cells[TH_IDX]; 
        // mapc2p_n[0] = R, mapc2p_n[1] = Theta, mapc2p_n[2] = Z_cylindrical
        double *mapc2p_n = gkyl_array_fetch(mapc2p_nodal, gkyl_range_idx(&nrange, cidx));
        double xn[3] = {psi, alpha, theta};
        double fout[3];
        mapc2p(0.0, xn, fout, 0);
        for (int i=0; i<3; ++i)
          TEST_CHECK( gkyl_compare( mapc2p_n[i], fout[i], 1e-6) );
      }
    }
  }

  // Check mc2nu_pos
  struct gkyl_array* mc2nu_pos_nodal = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange, &range, 3, mc2nu_pos_nodal, gk_geom->geo_corn.mc2nu_pos, false);
  for (int ia=nrange.lower[AL_IDX]; ia<=nrange.upper[AL_IDX]; ++ia){
    for (int ip=nrange.lower[PSI_IDX]; ip<=nrange.upper[PSI_IDX]; ++ip) {
      for (int it=nrange.lower[TH_IDX]; it<=nrange.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double psi = grid.lower[PSI_IDX] + ip*(grid.upper[PSI_IDX]-grid.lower[PSI_IDX])/grid.cells[PSI_IDX];
        double alpha = grid.lower[AL_IDX] + ia*(grid.upper[AL_IDX]-grid.lower[AL_IDX])/grid.cells[AL_IDX];
        double theta = grid.lower[TH_IDX] + it*(grid.upper[TH_IDX]-grid.lower[TH_IDX])/grid.cells[TH_IDX];
        double xn[3] = {psi, alpha, theta};
        double *mc2nu_pos_n = gkyl_array_fetch(mc2nu_pos_nodal, gkyl_range_idx(&nrange, cidx));
        for (int i=0; i<3; ++i) {
          TEST_CHECK( gkyl_compare( mc2nu_pos_n[i], xn[i], 1e-8) );
        }
      }
    }
  }

  // Check normals
  // Plus 3 away from axis to avoid errors
  double dels[2] = {1.0/sqrt(3), 1.0-1.0/sqrt(3) };
  double theta_lo = grid.lower[TH_IDX] + dels[1]*grid.dx[TH_IDX]/2.0;
  double psi_lo = grid.lower[PSI_IDX] + dels[1]*grid.dx[PSI_IDX]/2.0;
  double alpha_lo = grid.lower[AL_IDX] + dels[1]*grid.dx[AL_IDX]/2.0;
  struct gkyl_array* normals_nodal = gkyl_array_new(GKYL_DOUBLE, 9*basis.num_basis, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange_quad_interior, &range,  9*basis.num_basis, normals_nodal, gk_geom->geo_int.normals, true);
  for (int ia=nrange_quad_interior.lower[AL_IDX]; ia<=nrange_quad_interior.upper[AL_IDX]; ++ia){
    for (int ip=nrange_quad_interior.lower[PSI_IDX]; ip<=nrange_quad_interior.upper[PSI_IDX]; ++ip) {
      for (int it=nrange_quad_interior.lower[TH_IDX]; it<=nrange_quad_interior.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *normals_n = gkyl_array_fetch(normals_nodal, gkyl_range_idx(&nrange_quad_interior, cidx));
        double psi= calc_running_coord(psi_lo, ia-nrange_quad_interior.lower[PSI_IDX], grid.dx[PSI_IDX]);
        double alpha= calc_running_coord(alpha_lo, ia-nrange_quad_interior.lower[AL_IDX], grid.dx[AL_IDX]);
        double theta= calc_running_coord(theta_lo, ia-nrange_quad_interior.lower[TH_IDX], grid.dx[TH_IDX]);
        double xn[3] = {psi, alpha, theta};
        double fout[9];
        exact_normals(0.0, xn, fout, 0);
        double tol;
        for (int i=0; i<9; ++i) {
          TEST_CHECK( gkyl_compare( normals_n[i], fout[i], 1e-3) );
        }
      }
    }
  }

  gkyl_array_release(bhat_nodal);
  gkyl_array_release(dualmag_nodal);
  gkyl_array_release(bmag_nodal);
  gkyl_array_release(cmag_nodal);
  gkyl_array_release(gij_nodal);
  gkyl_array_release(gij_contra_nodal);
  gkyl_array_release(jacobgeo_nodal);
  gkyl_array_release(jacobgeo_inv_nodal);
  gkyl_array_release(jacobtot_nodal);
  gkyl_array_release(jacobtot_inv_nodal);
  gkyl_array_release(mapc2p_nodal);
  gkyl_array_release(mapc2p_nodal_interior);
  gkyl_array_release(mc2nu_pos_nodal);
  gkyl_array_release(normals_nodal);
  gkyl_nodal_ops_release(n2m);
  gkyl_position_map_release(pmap);
  gkyl_gk_geometry_release(gk_geom);
}

void
test_tok_asdex_qprofile_core_ho()
{
  double clower[] = { -0.09, -0.01, -M_PI+1e-14 };
  double cupper[] = {0.14975, 0.01, M_PI-1e-14 };
  int ccells[] = { 4, 1, 4 };

  int cpoly_order = 1;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  struct gkyl_rect_grid cgrid;
  struct gkyl_range clocal, clocal_ext;
  struct gkyl_basis cbasis;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_efit_inp efit_inp = {
    // psiRZ and related inputs
    .filepath = "gyrokinetic/data/eqdsk/asdex.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
  };
  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_GEOMETRY_TOKAMAK_CORE,
    .rmin = 0.0,
    .rmax = 5.0,
    .rclose = 2.5,
    .rright = 2.5,
    .rleft = 0.7,
    .zmin = -1.3,
    .zmax = 1.0,
    .zmin_left = -1.2,
    .zmin_right = -1.0,
  };
  // Initialize geometry
  struct gkyl_position_map *pmap = gkyl_position_map_null_new();

  struct gkyl_gk_geometry_inp geometry_input = {
    .geometry_id  = GKYL_GEOMETRY_TOKAMAK,
    .efit_info = efit_inp,
    .tok_grid_info = ginp,
    .position_map = pmap,
    .grid = cgrid,
    .local = clocal,
    .local_ext = clocal_ext,
    .global = clocal,
    .global_ext = clocal_ext,
    .basis = cbasis,
    .geo_grid = cgrid,
    .geo_local = clocal,
    .geo_local_ext = clocal_ext,
    .geo_global = clocal,
    .geo_global_ext = clocal_ext,
    .geo_basis = cbasis,
  };

  struct gk_geometry *gk_geom = gkyl_gk_geometry_tok_new(&geometry_input);
  write_geometry(gk_geom, cgrid, clocal, "asdex_core");

//  // Create Nodal Range and Grid and Write Nodal Coordinates
//  struct gkyl_range nrange;
//  gkyl_gk_geometry_init_nodal_range(&nrange, &clocal, cpoly_order);
//  struct gkyl_array* mc2p_nodal = gkyl_array_new(GKYL_DOUBLE, 3, nrange.volume);
//  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&cbasis, &cgrid, false);
//  gkyl_nodal_ops_m2n(n2m, &cbasis, &cgrid, &nrange, &clocal, 3, mc2p_nodal, gk_geom->geo_int.mc2p, true);
//  gkyl_nodal_ops_release(n2m);
//  struct gkyl_rect_grid ngrid;
//  gkyl_gk_geometry_init_nodal_grid(&ngrid, &cgrid, &nrange);
//
//  gkyl_grid_sub_array_write(&ngrid, &nrange, 0,  mc2p_nodal, "asdex_core_nodes.gkyl");
//  gkyl_array_release(mc2p_nodal);

  gkyl_gk_geometry_release(gk_geom);
  gkyl_position_map_release(pmap);
}

void
test_tok_asdex_qprofile_sol_ho()
{
  double clower[] = { 0.16, -0.01, -M_PI+1e-14 };
  double cupper[] = {0.17501, 0.01, M_PI-1e-14 };
  int ccells[] = { 4, 1, 4 };

  int cpoly_order = 1;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  struct gkyl_rect_grid cgrid;
  struct gkyl_range clocal, clocal_ext;
  struct gkyl_basis cbasis;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_efit_inp efit_inp = {
    // psiRZ and related inputs
    .filepath = "gyrokinetic/data/eqdsk/asdex.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
  };

  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_GEOMETRY_TOKAMAK_LSN_SOL,
    .rmin = 0.0,
    .rmax = 5.0,
    .rclose = 2.5,
    .rright = 2.5,
    .rleft = 0.7,
    .zmin = -1.3,
    .zmax = 1.0,
    .zmin_left = -1.2,
    .zmin_right = -1.0,
  };

  struct gkyl_position_map *pmap = gkyl_position_map_null_new();

  struct gkyl_gk_geometry_inp geometry_inp = {
    .geometry_id  = GKYL_GEOMETRY_TOKAMAK,
    .efit_info = efit_inp,
    .tok_grid_info = ginp,
    .position_map = pmap,
    .grid = cgrid,
    .local = clocal,
    .local_ext = clocal_ext,
    .global = clocal,
    .global_ext = clocal_ext,
    .basis = cbasis,
    .geo_grid = cgrid,
    .geo_local = clocal,
    .geo_local_ext = clocal_ext,
    .geo_global = clocal,
    .geo_global_ext = clocal_ext,
    .geo_basis = cbasis,
  };

  struct gk_geometry* gk_geom = gkyl_gk_geometry_tok_new(&geometry_inp);
  write_geometry(gk_geom, cgrid, clocal, "asdex_sol");
  gkyl_gk_geometry_release(gk_geom);
  gkyl_position_map_release(pmap);
}

TEST_LIST = {
  { "test_tok_elliptical_ho", test_tok_elliptical_ho},
  { "test_tok_3x_p1_straight_cylinder_ho", test_tok_3x_p1_straight_cylinder_ho},
  { "test_tok_asdex_qprofile_core_ho", test_tok_asdex_qprofile_core_ho},
  { "test_tok_asdex_qprofile_sol_ho", test_tok_asdex_qprofile_sol_ho},
  { NULL, NULL },
};
