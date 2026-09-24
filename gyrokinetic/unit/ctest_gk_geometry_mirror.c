#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_priv.h>
#include <gkyl_gk_geometry_mirror.h>
#include <gkyl_nodal_ops.h>
#include <gkyl_position_map.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

void
write_geometry(
  gk_geometry *up, struct gkyl_rect_grid grid, struct gkyl_range local, const char *name
)
{
  const char *fmt = "%s-%s.gkyl";
  int sz = gkyl_calc_strlen(fmt, name, "jacobtot_inv");
  char fileNm[sz + 1]; // ensure no buffer overflow

  sprintf(fileNm, fmt, name, "mapc2p");
  gkyl_grid_sub_array_write(&grid, &local, 0, up->geo_corn.mc2p, fileNm);
  sprintf(fileNm, fmt, name, "mapc2nu");
  gkyl_grid_sub_array_write(&grid, &local, 0, up->geo_corn.mc2nu_pos, fileNm);
  sprintf(fileNm, fmt, name, "bmag_corn");
  gkyl_grid_sub_array_write(&grid, &local, 0, up->geo_corn.bmag, fileNm);
  sprintf(fileNm, fmt, name, "bmag");
  gkyl_grid_sub_array_write(&grid, &local, 0, up->geo_int.bmag, fileNm);
  sprintf(fileNm, fmt, name, "g_ij");
  gkyl_grid_sub_array_write(&grid, &local, 0, up->geo_int.g_ij, fileNm);
  sprintf(fileNm, fmt, name, "dxdz");
  gkyl_grid_sub_array_write(&grid, &local, 0, up->geo_int.dxdz, fileNm);
  sprintf(fileNm, fmt, name, "dzdx");
  gkyl_grid_sub_array_write(&grid, &local, 0, up->geo_int.dzdx, fileNm);
  sprintf(fileNm, fmt, name, "normals");
  gkyl_grid_sub_array_write(&grid, &local, 0, up->geo_int.normals, fileNm);
  sprintf(fileNm, fmt, name, "jacobgeo");
  gkyl_grid_sub_array_write(&grid, &local, 0, up->geo_int.jacobgeo, fileNm);
  sprintf(fileNm, fmt, name, "jacobgeo_inv");
  gkyl_grid_sub_array_write(&grid, &local, 0, up->geo_int.jacobgeo_inv, fileNm);
  sprintf(fileNm, fmt, name, "gij");
  gkyl_grid_sub_array_write(&grid, &local, 0, up->geo_int.gij, fileNm);
  sprintf(fileNm, fmt, name, "b_i");
  gkyl_grid_sub_array_write(&grid, &local, 0, up->geo_int.b_i, fileNm);
  sprintf(fileNm, fmt, name, "bcart");
  gkyl_grid_sub_array_write(&grid, &local, 0, up->geo_int.bcart, fileNm);
  sprintf(fileNm, fmt, name, "cmag");
  gkyl_grid_sub_array_write(&grid, &local, 0, up->geo_int.cmag, fileNm);
  sprintf(fileNm, fmt, name, "jacobtot");
  gkyl_grid_sub_array_write(&grid, &local, 0, up->geo_int.jacobtot, fileNm);
  sprintf(fileNm, fmt, name, "jacobtot_inv");
  gkyl_grid_sub_array_write(&grid, &local, 0, up->geo_int.jacobtot_inv, fileNm);
  sprintf(fileNm, fmt, name, "gxxj");
  gkyl_grid_sub_array_write(&grid, &local, 0, up->geo_int.gxxj, fileNm);
  sprintf(fileNm, fmt, name, "gxyj");
  gkyl_grid_sub_array_write(&grid, &local, 0, up->geo_int.gxyj, fileNm);
  sprintf(fileNm, fmt, name, "gyyj");
  gkyl_grid_sub_array_write(&grid, &local, 0, up->geo_int.gyyj, fileNm);
  sprintf(fileNm, fmt, name, "gxzj");
  gkyl_grid_sub_array_write(&grid, &local, 0, up->geo_int.gxzj, fileNm);
  sprintf(fileNm, fmt, name, "eps2");
  gkyl_grid_sub_array_write(&grid, &local, 0, up->geo_int.eps2, fileNm);

  // Create Nodal Range and Grid and Write Nodal Coordinates
  struct gkyl_range nrange;
  gkyl_gk_geometry_init_nodal_range(&nrange, &local, up->basis.poly_order);
  struct gkyl_array *mc2p_nodal = gkyl_array_new(GKYL_DOUBLE, 3, nrange.volume);
  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&up->basis, &grid, false);
  gkyl_nodal_ops_m2n(
    n2m, &up->basis, &grid, &nrange, &local, 3, mc2p_nodal, up->geo_corn.mc2p, false
  );
  gkyl_nodal_ops_release(n2m);
  struct gkyl_rect_grid ngrid;
  gkyl_gk_geometry_init_nodal_grid(&ngrid, &grid, &nrange);
  sprintf(fileNm, fmt, name, "nodes");
  gkyl_grid_sub_array_write(&ngrid, &nrange, 0, mc2p_nodal, fileNm);
  gkyl_array_release(mc2p_nodal);
}

void
test_mirror_load_geometry_ho(void)
{
  struct gkyl_efit_inp inp = {
    // psiRZ and related inputs
    .filepath = "gyrokinetic/data/eqdsk/wham.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect = true,
  };

  clock_t start, end;
  double cpu_time_used;
  start = clock();

  double clower[] = {1e-10, -0.01, -2.0};
  double cupper[] = {3e-3, 0.01, 2.0};

  int ccells[] = {4, 4, 8};

  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);

  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = {1, 1, 1};
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);

  int cpoly_order = 1;
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_mirror_geo_grid_inp ginp = {
    .filename_psi = "gyrokinetic/data/unit/wham_hires.geqdsk_psi.gkyl", // psi file to use
    .rclose = 0.2, // closest R to region of interest
    .zmin = -2.0, // Z of lower boundary
    .zmax = 2.0, // Z of upper boundary
    .include_axis = false, // Include R=0 axis in grid
    .fl_coord = GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z, // coordinate system for psi grid
  };

  struct gkyl_position_map *pmap = gkyl_position_map_null_new();

  struct gkyl_gk_geometry_inp geometry_inp = {
    .geometry_id = GKYL_GEOMETRY_MIRROR,
    .mirror_grid_info = ginp,
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

  struct gk_geometry *up = gkyl_gk_geometry_mirror_new(&geometry_inp);
  //write_geometry(up, cgrid, clocal, "whamlores");
  gkyl_gk_geometry_release(up);
  gkyl_position_map_release(pmap);

  end = clock();
  cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
}

// def psi_f(R, Z):
//     Bmag = 0.5
//     return Bmag/2 * R**2

// Functions for test_3x_straight_cylinder
void
mapc2p(double t, const double *xn, double *GKYL_RESTRICT fout, void *ctx)
{
  double psi = xn[0], alpha = xn[1], zeta = xn[2];
  fout[0] = sqrt(psi * 4); // Function fed is psi = 0.25 * R^2 from the efit file
  fout[1] =
    zeta; // Note that this does not have pi-1e-2 in it because the coordinate zeta is always defined -pi to pi
  fout[2] = -alpha; // There is a minus due to conventions
}

void
exact_gij(double t, const double *xn, double *GKYL_RESTRICT fout, void *ctx)
{
  double r = xn[0], theta = xn[1], phi = xn[2];
  double psi = r * r / 4;
  fout[0] = 1 / psi; // g_11
  fout[1] = 0.0; // g_12
  fout[2] = 0.0; // g_13
  fout[3] = r * r; // g_22
  fout[4] = 0.0; // g_23
  fout[5] = 1, 0; // g_33
}

void
exact_g_contra_ij(double t, const double *xn, double *GKYL_RESTRICT fout, void *ctx)
{
  double r = xn[0], theta = xn[1], phi = xn[2];
  double psi = r * r / 4;
  fout[0] = r * r / 4; // g_11
  fout[1] = 0.0; // g_12
  fout[2] = 0.0; // g_13
  fout[3] = 1 / psi / 4; // g_22
  fout[4] = 0.0; // g_23
  fout[5] = 1.0; // g_33
}

void
exact_dual_magnitude(double t, const double *xn, double *GKYL_RESTRICT fout, void *ctx)
{
  double r = xn[0], theta = xn[1], phi = xn[2];
  double psi = r * r / 4;
  fout[0] = r / 2;
  fout[1] = 1 / (2 * sqrt(psi));
  fout[2] = 1.0;
}

void
exact_normals(double t, const double *xn, double *GKYL_RESTRICT fout, void *ctx)
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

void
exact_jacobian(double t, const double *xn, double *GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 2.0;
}

void
bmag_func(double t, const double *xn, double *GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.5;
}

void
test_mirror_3x_p1_straight_cylinder_ho(void)
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
  int Nz = 10;

  double lower[3] = {psiMin, -1.0, -1.0};
  double upper[3] = {psiMax, 1.0, 1.0};
  // int cells[3] = { 18, 18, Nz };
  int cells[3] = {2, 2, 2};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);

  struct gkyl_range ext_range, range;
  int nghost[3] = {1, 1, 1};
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  struct gkyl_mirror_geo_grid_inp ginp = {
    .filename_psi = "gyrokinetic/data/unit/straight_cylinder.geqdsk_psi.gkyl", // psi file to use
    .rclose = 0.5, // closest R to region of interest
    .zmin = -1.0, // Z of lower boundary
    .zmax = 1.0, // Z of upper boundary
    .include_axis = false, // Include R=0 axis in grid
    .fl_coord = GKYL_GEOMETRY_MIRROR_GRID_GEN_PSI_CART_Z, // coordinate system for psi grid
  };

  struct gkyl_position_map *pos_map = gkyl_position_map_null_new();

  // Initialize geometry
  struct gkyl_gk_geometry_inp geometry_input = {
    .geometry_id = GKYL_GEOMETRY_MIRROR,
    .mirror_grid_info = ginp,
    .position_map = pos_map,
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

  struct gk_geometry *gk_geom = gkyl_gk_geometry_mirror_new(&geometry_input);

  // write_geometry(gk_geom, grid, range, "straight_cylinder");

  // Define nodal operations
  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  int cidx[3];
  int nodes[] = {1, 1, 1};
  for (int d = 0; d < grid.ndim; ++d) {
    nodes[d] = grid.cells[d] + 1;
  }
  struct gkyl_range nrange;
  gkyl_range_init_from_shape(&nrange, grid.ndim, nodes);

  int nodes_quad_interior[] = {1, 1, 1};
  int num_quad_points = poly_order + 1;
  for (int d = 0; d < grid.ndim; ++d) {
    nodes_quad_interior[d] = grid.cells[d] * num_quad_points;
  }
  struct gkyl_range nrange_quad_interior;
  gkyl_range_init_from_shape(&nrange_quad_interior, grid.ndim, nodes_quad_interior);

  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&basis, &grid, false);

  // First, let's validate all objects in gk_geom_corn

  // Check mapc2p
  struct gkyl_array *mapc2p_nodal = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
  gkyl_nodal_ops_m2n(
    n2m, &basis, &grid, &nrange, &range, 3, mapc2p_nodal, gk_geom->geo_corn.mc2p, false
  );
  for (int ia = nrange.lower[AL_IDX]; ia <= nrange.upper[AL_IDX]; ++ia) {
    for (int ip = nrange.lower[PSI_IDX]; ip <= nrange.upper[PSI_IDX]; ++ip) {
      for (int it = nrange.lower[TH_IDX]; it <= nrange.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double psi = grid.lower[PSI_IDX] +
                     ip * (grid.upper[PSI_IDX] - grid.lower[PSI_IDX]) / grid.cells[PSI_IDX];
        double alpha =
          grid.lower[AL_IDX] + ia * (grid.upper[AL_IDX] - grid.lower[AL_IDX]) / grid.cells[AL_IDX];
        double theta =
          grid.lower[TH_IDX] + it * (grid.upper[TH_IDX] - grid.lower[TH_IDX]) / grid.cells[TH_IDX];

        // geo_corn.mapc2p_nodal
        double *mapc2p_n =
          gkyl_array_fetch(gk_geom->geo_corn.mc2p_nodal, gkyl_range_idx(&nrange, cidx));
        double xn[3] = {psi, -alpha, theta};
        double fout[3];
        mapc2p(0.0, xn, fout, 0);
        for (int i = 0; i < 3; ++i) {
          TEST_CHECK(gkyl_compare(mapc2p_n[i], fout[i], 1e-8));
        }

        // geo_corn.mc2p
        double *mapc2p = gkyl_array_fetch(mapc2p_nodal, gkyl_range_idx(&nrange, cidx));
        for (int i = 0; i < 3; ++i) {
          TEST_CHECK(gkyl_compare(mapc2p[i], fout[i], 1e-8));
        }

        // geo_corn.mc2p_deflated not tested because it's populated elsewhere
      }
    }
  }

  // Check mc2nu_pos
  struct gkyl_array *mc2nu_pos_nodal = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
  gkyl_nodal_ops_m2n(
    n2m, &basis, &grid, &nrange, &range, 3, mc2nu_pos_nodal, gk_geom->geo_corn.mc2nu_pos, false
  );
  for (int ia = nrange.lower[AL_IDX]; ia <= nrange.upper[AL_IDX]; ++ia) {
    for (int ip = nrange.lower[PSI_IDX]; ip <= nrange.upper[PSI_IDX]; ++ip) {
      for (int it = nrange.lower[TH_IDX]; it <= nrange.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double psi = grid.lower[PSI_IDX] +
                     ip * (grid.upper[PSI_IDX] - grid.lower[PSI_IDX]) / grid.cells[PSI_IDX];
        double alpha =
          grid.lower[AL_IDX] + ia * (grid.upper[AL_IDX] - grid.lower[AL_IDX]) / grid.cells[AL_IDX];
        double theta =
          grid.lower[TH_IDX] + it * (grid.upper[TH_IDX] - grid.lower[TH_IDX]) / grid.cells[TH_IDX];
        double xn[3] = {psi, alpha, theta};

        // geo_corn.mc2nu_pos
        double *mc2nu_pos_n = gkyl_array_fetch(mc2nu_pos_nodal, gkyl_range_idx(&nrange, cidx));
        for (int i = 0; i < 3; ++i) {
          TEST_CHECK(gkyl_compare(mc2nu_pos_n[i], xn[i], 1e-8));
        }

        // geo_corn.mc2nu_pos_nodal
        double *mc2nu_pos =
          gkyl_array_fetch(gk_geom->geo_corn.mc2nu_pos_nodal, gkyl_range_idx(&nrange, cidx));
        for (int i = 0; i < 3; ++i) {
          TEST_CHECK(gkyl_compare(mc2nu_pos[i], xn[i], 1e-8));
        }

        // geo_corn.mc2nu_pos_deflated not tested because it's populated elsewhere
      }
    }
  }

  // Check bmag is what it should be
  struct gkyl_array *bmag_nodal = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
  gkyl_nodal_ops_m2n(
    n2m, &basis, &grid, &nrange, &range, 1, bmag_nodal, gk_geom->geo_corn.bmag, false
  );
  for (int ia = nrange.lower[AL_IDX]; ia <= nrange.upper[AL_IDX]; ++ia) {
    for (int ip = nrange.lower[PSI_IDX]; ip <= nrange.upper[PSI_IDX]; ++ip) {
      for (int it = nrange.lower[TH_IDX]; it <= nrange.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double psi = grid.lower[PSI_IDX] +
                     ip * (grid.upper[PSI_IDX] - grid.lower[PSI_IDX]) / grid.cells[PSI_IDX];
        double alpha =
          grid.lower[AL_IDX] + ia * (grid.upper[AL_IDX] - grid.lower[AL_IDX]) / grid.cells[AL_IDX];
        double theta =
          grid.lower[TH_IDX] + it * (grid.upper[TH_IDX] - grid.lower[TH_IDX]) / grid.cells[TH_IDX];
        double xn[3] = {psi, -alpha, theta};

        // geo_corn.bmag
        double *bmag_n = gkyl_array_fetch(bmag_nodal, gkyl_range_idx(&nrange, cidx));
        double bmag_anal[1];
        bmag_func(0, xn, bmag_anal, 0);
        TEST_CHECK(gkyl_compare(bmag_n[0], bmag_anal[0], 1e-7));

        // geo_corn.bmag_nodal
        double *bmag =
          gkyl_array_fetch(gk_geom->geo_corn.bmag_nodal, gkyl_range_idx(&nrange, cidx));
        TEST_CHECK(gkyl_compare(bmag[0], bmag_anal[0], 1e-7));
      }
    }
  }

  // Second let's validate geo_int. Many coponents still need to be added

  // Plus 3 away from axis to avoid errors
  double dels[2] = {1.0 / sqrt(3), 1.0 - 1.0 / sqrt(3)};
  double theta_lo = grid.lower[TH_IDX] + dels[1] * grid.dx[TH_IDX] / 2.0;
  double psi_lo = grid.lower[PSI_IDX] + dels[1] * grid.dx[PSI_IDX] / 2.0;
  double alpha_lo = grid.lower[AL_IDX] + dels[1] * grid.dx[AL_IDX] / 2.0;

  // geo_int.mc2p
  struct gkyl_array *mc2p_nodal_interior =
    gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(
    n2m, &basis, &grid, &nrange_quad_interior, &range, 3, mc2p_nodal_interior,
    gk_geom->geo_int.mc2p, true
  );
  // Check that the mapc2p is what it should be
  for (int ia = nrange_quad_interior.lower[AL_IDX]; ia <= nrange_quad_interior.upper[AL_IDX];
       ++ia) {
    for (int ip = nrange_quad_interior.lower[PSI_IDX]; ip <= nrange_quad_interior.upper[PSI_IDX];
         ++ip) {
      for (int it = nrange_quad_interior.lower[TH_IDX]; it <= nrange_quad_interior.upper[TH_IDX];
           ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *mc2p_n =
          gkyl_array_fetch(mc2p_nodal_interior, gkyl_range_idx(&nrange_quad_interior, cidx));

        double psi =
          calc_running_coord(psi_lo, ip - nrange_quad_interior.lower[PSI_IDX], grid.dx[PSI_IDX]);
        double alpha =
          calc_running_coord(alpha_lo, ia - nrange_quad_interior.lower[AL_IDX], grid.dx[AL_IDX]);
        double theta =
          calc_running_coord(theta_lo, it - nrange_quad_interior.lower[TH_IDX], grid.dx[TH_IDX]);
        double xn[3] = {psi, -alpha, theta};
        double fout[3];
        mapc2p(0.0, xn, fout, 0);
        for (int i = 0; i < 3; ++i) {
          TEST_CHECK(gkyl_compare(mc2p_n[i], fout[i], 1e-8));
        }

        // geo_int.mc2p_nodal
        double *mc2p = gkyl_array_fetch(
          gk_geom->geo_int.mc2p_nodal, gkyl_range_idx(&nrange_quad_interior, cidx)
        );
        for (int i = 0; i < 3; ++i) {
          TEST_CHECK(gkyl_compare(mc2p[i], fout[i], 1e-8));
        }
      }
    }
  }

  // geo_int.bmag
  struct gkyl_array *bmag_nodal_interior =
    gkyl_array_new(GKYL_DOUBLE, 1, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(
    n2m, &basis, &grid, &nrange_quad_interior, &range, 1, bmag_nodal_interior,
    gk_geom->geo_int.bmag, true
  );
  for (int ia = nrange_quad_interior.lower[AL_IDX]; ia <= nrange_quad_interior.upper[AL_IDX];
       ++ia) {
    for (int ip = nrange_quad_interior.lower[PSI_IDX]; ip <= nrange_quad_interior.upper[PSI_IDX];
         ++ip) {
      for (int it = nrange_quad_interior.lower[TH_IDX]; it <= nrange_quad_interior.upper[TH_IDX];
           ++it) {
        double *bmag_n =
          gkyl_array_fetch(bmag_nodal_interior, gkyl_range_idx(&nrange_quad_interior, cidx));
        double fout[1];
        bmag_func(0.0, NULL, fout, 0);
        TEST_CHECK(gkyl_compare(bmag_n[0], fout[0], 1e-8));

        // geo_int.bmag_nodal
        double *bmag = gkyl_array_fetch(
          gk_geom->geo_int.bmag_nodal, gkyl_range_idx(&nrange_quad_interior, cidx)
        );
        TEST_CHECK(gkyl_compare(bmag[0], fout[0], 1e-8));
      }
    }
  }

  // geo_int.g_ij
  struct gkyl_array *gij_nodal = gkyl_array_new(GKYL_DOUBLE, 6, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(
    n2m, &basis, &grid, &nrange_quad_interior, &range, 6, gij_nodal, gk_geom->geo_int.g_ij, true
  );
  for (int ia = nrange_quad_interior.lower[AL_IDX]; ia <= nrange_quad_interior.upper[AL_IDX];
       ++ia) {
    for (int ip = nrange_quad_interior.lower[PSI_IDX]; ip <= nrange_quad_interior.upper[PSI_IDX];
         ++ip) {
      for (int it = nrange_quad_interior.lower[TH_IDX]; it <= nrange_quad_interior.upper[TH_IDX];
           ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *gij_n = gkyl_array_fetch(gij_nodal, gkyl_range_idx(&nrange_quad_interior, cidx));
        double *mapc2p_n =
          gkyl_array_fetch(mc2p_nodal_interior, gkyl_range_idx(&nrange_quad_interior, cidx));
        double r = mapc2p_n[0];
        double xn[3] = {r, 0.0, 0.0};
        double fout[6];
        exact_gij(0.0, xn, fout, 0);
        for (int i = 0; i < 6; ++i) {
          TEST_CHECK(gkyl_compare(gij_n[i], fout[i], 1e-6));
        }

        // geo_int.g_ij_nodal
        double *gij = gkyl_array_fetch(
          gk_geom->geo_int.g_ij_nodal, gkyl_range_idx(&nrange_quad_interior, cidx)
        );
        for (int i = 0; i < 6; ++i) {
          TEST_CHECK(gkyl_compare(gij[i], fout[i], 1e-6));
        }
      }
    }
  }

  // geo_int.gij
  struct gkyl_array *gij_contra_nodal = gkyl_array_new(GKYL_DOUBLE, 6, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(
    n2m, &basis, &grid, &nrange_quad_interior, &range, 6, gij_contra_nodal, gk_geom->geo_int.gij,
    true
  );
  for (int ia = nrange_quad_interior.lower[AL_IDX]; ia <= nrange_quad_interior.upper[AL_IDX];
       ++ia) {
    for (int ip = nrange_quad_interior.lower[PSI_IDX]; ip <= nrange_quad_interior.upper[PSI_IDX];
         ++ip) {
      for (int it = nrange_quad_interior.lower[TH_IDX]; it <= nrange_quad_interior.upper[TH_IDX];
           ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *gij_contra_n =
          gkyl_array_fetch(gij_contra_nodal, gkyl_range_idx(&nrange_quad_interior, cidx));
        double *mapc2p_n =
          gkyl_array_fetch(mc2p_nodal_interior, gkyl_range_idx(&nrange_quad_interior, cidx));
        double r = mapc2p_n[0];
        double xn[3] = {r, 0.0, 0.0};
        double fout[6];
        exact_g_contra_ij(0.0, xn, fout, 0);
        for (int i = 0; i < 6; ++i) {
          TEST_CHECK(gkyl_compare(gij_contra_n[i], fout[i], 1e-6));
        }

        // geo_int.gij_nodal does not exist
      }
    }
  }

  // Check that |bhat|=1 at nodes
  struct gkyl_array *bhat_nodal =
    gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(
    n2m, &basis, &grid, &nrange_quad_interior, &range, 3, bhat_nodal, gk_geom->geo_int.bcart, true
  );
  for (int ia = nrange_quad_interior.lower[AL_IDX]; ia <= nrange_quad_interior.upper[AL_IDX];
       ++ia) {
    for (int ip = nrange_quad_interior.lower[PSI_IDX]; ip <= nrange_quad_interior.upper[PSI_IDX];
         ++ip) {
      for (int it = nrange_quad_interior.lower[TH_IDX]; it <= nrange_quad_interior.upper[TH_IDX];
           ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *bhat_n = gkyl_array_fetch(bhat_nodal, gkyl_range_idx(&nrange_quad_interior, cidx));
        double bhat_mag =
          sqrt(bhat_n[0] * bhat_n[0] + bhat_n[1] * bhat_n[1] + bhat_n[2] * bhat_n[2]);
        TEST_CHECK(gkyl_compare(bhat_mag, 1.0, 1e-12));

        // geo_int.bcart_nodal
        double *bhat = gkyl_array_fetch(
          gk_geom->geo_int.bcart_nodal, gkyl_range_idx(&nrange_quad_interior, cidx)
        );
        double bhat_mag_nodal = sqrt(bhat[0] * bhat[0] + bhat[1] * bhat[1] + bhat[2] * bhat[2]);
        TEST_CHECK(gkyl_compare(bhat_mag_nodal, 1.0, 1e-12));
      }
    }
  }

  // Check that the duals are what they should be
  // There are errors at low psi, so we shift away from the axis by a few cells
  struct gkyl_array *dualmag_nodal =
    gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(
    n2m, &basis, &grid, &nrange_quad_interior, &range, 3, dualmag_nodal, gk_geom->geo_int.dualmag,
    true
  );
  for (int ia = nrange_quad_interior.lower[AL_IDX]; ia <= nrange_quad_interior.upper[AL_IDX];
       ++ia) {
    for (int ip = nrange_quad_interior.lower[PSI_IDX] + 3;
         ip <= nrange_quad_interior.upper[PSI_IDX]; ++ip) {
      for (int it = nrange_quad_interior.lower[TH_IDX]; it <= nrange_quad_interior.upper[TH_IDX];
           ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *dualmag_n =
          gkyl_array_fetch(dualmag_nodal, gkyl_range_idx(&nrange_quad_interior, cidx));
        double *mapc2p_n =
          gkyl_array_fetch(mc2p_nodal_interior, gkyl_range_idx(&nrange_quad_interior, cidx));
        double xn[3] = {mapc2p_n[0], mapc2p_n[1], mapc2p_n[2]};
        double dualmag_anal[3];
        exact_dual_magnitude(0, xn, dualmag_anal, 0);
        TEST_CHECK(gkyl_compare(dualmag_n[0], dualmag_anal[0], 1e-6));
        TEST_CHECK(gkyl_compare(dualmag_n[1], dualmag_anal[1], 1e-6));
        TEST_CHECK(gkyl_compare(dualmag_n[2], dualmag_anal[2], 1e-6));

        // geo_int.dualmag_nodal
        double *dualmag = gkyl_array_fetch(
          gk_geom->geo_int.dualmag_nodal, gkyl_range_idx(&nrange_quad_interior, cidx)
        );
        TEST_CHECK(gkyl_compare(dualmag[0], dualmag_anal[0], 1e-6));
        TEST_CHECK(gkyl_compare(dualmag[1], dualmag_anal[1], 1e-6));
        TEST_CHECK(gkyl_compare(dualmag[2], dualmag_anal[2], 1e-6));
      }
    }
  }

  // Check cmag = 1
  struct gkyl_array *cmag_nodal = gkyl_array_new(GKYL_DOUBLE, 1, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(
    n2m, &basis, &grid, &nrange_quad_interior, &range, 1, cmag_nodal, gk_geom->geo_int.cmag, true
  );
  for (int ia = nrange_quad_interior.lower[AL_IDX]; ia <= nrange_quad_interior.upper[AL_IDX];
       ++ia) {
    for (int ip = nrange_quad_interior.lower[PSI_IDX]; ip <= nrange_quad_interior.upper[PSI_IDX];
         ++ip) {
      for (int it = nrange_quad_interior.lower[TH_IDX]; it <= nrange_quad_interior.upper[TH_IDX];
           ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *cmag_n = gkyl_array_fetch(cmag_nodal, gkyl_range_idx(&nrange, cidx));
        TEST_CHECK(gkyl_compare(cmag_n[0], 1.0, 1e-8));
      }
    }
  }

  // Check that Jacobgeo is what it should be. This is the Jacobian for the problem
  struct gkyl_array *jacobgeo_nodal =
    gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(
    n2m, &basis, &grid, &nrange_quad_interior, &range, 1, jacobgeo_nodal, gk_geom->geo_int.jacobgeo,
    true
  );
  for (int ia = nrange_quad_interior.lower[AL_IDX]; ia <= nrange_quad_interior.upper[AL_IDX];
       ++ia) {
    for (int ip = nrange_quad_interior.lower[PSI_IDX]; ip <= nrange_quad_interior.upper[PSI_IDX];
         ++ip) {
      for (int it = nrange_quad_interior.lower[TH_IDX]; it <= nrange_quad_interior.upper[TH_IDX];
           ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *jacobgeo_n = gkyl_array_fetch(jacobgeo_nodal, gkyl_range_idx(&nrange, cidx));
        double fout[1];
        exact_jacobian(0.0, NULL, fout, 0);
        TEST_CHECK(gkyl_compare(jacobgeo_n[0], fout[0], 1e-6));

        // geo_int.jacobgeo_nodal
        double *jacobgeo =
          gkyl_array_fetch(gk_geom->geo_int.jacobgeo_nodal, gkyl_range_idx(&nrange, cidx));
        TEST_CHECK(gkyl_compare(jacobgeo[0], fout[0], 1e-6));
      }
    }
  }

  // Check jacobgeo_inv
  struct gkyl_array *jacobgeo_inv_nodal =
    gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(
    n2m, &basis, &grid, &nrange_quad_interior, &range, 1, jacobgeo_inv_nodal,
    gk_geom->geo_int.jacobgeo_inv, true
  );
  for (int ia = nrange_quad_interior.lower[AL_IDX]; ia <= nrange_quad_interior.upper[AL_IDX];
       ++ia) {
    for (int ip = nrange_quad_interior.lower[PSI_IDX]; ip <= nrange_quad_interior.upper[PSI_IDX];
         ++ip) {
      for (int it = nrange_quad_interior.lower[TH_IDX]; it <= nrange_quad_interior.upper[TH_IDX];
           ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *jacobgeo_inv_n =
          gkyl_array_fetch(jacobgeo_inv_nodal, gkyl_range_idx(&nrange, cidx));
        double fout[1];
        exact_jacobian(0.0, NULL, fout, 0);
        TEST_CHECK(gkyl_compare(jacobgeo_inv_n[0], 1 / fout[0], 1e-6));
      }
    }
  }

  // Check jacobtot
  struct gkyl_array *jacobtot_nodal =
    gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(
    n2m, &basis, &grid, &nrange_quad_interior, &range, 1, jacobtot_nodal, gk_geom->geo_int.jacobtot,
    true
  );
  for (int ia = nrange_quad_interior.lower[AL_IDX]; ia <= nrange_quad_interior.upper[AL_IDX];
       ++ia) {
    for (int ip = nrange_quad_interior.lower[PSI_IDX]; ip <= nrange_quad_interior.upper[PSI_IDX];
         ++ip) {
      for (int it = nrange_quad_interior.lower[TH_IDX]; it <= nrange_quad_interior.upper[TH_IDX];
           ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *jacobtot_n =
          gkyl_array_fetch(jacobtot_nodal, gkyl_range_idx(&nrange_quad_interior, cidx));
        // mapc2p_n[0] = x, mapc2p_n[1] = y, mapc2p_n[2] = z
        double *mapc2p_n = gkyl_array_fetch(mapc2p_nodal, gkyl_range_idx(&nrange, cidx));
        double foutJ[1];
        exact_jacobian(0.0, NULL, foutJ, 0);
        double foutB[1];
        bmag_func(0.0, NULL, foutB, 0);
        TEST_CHECK(gkyl_compare(jacobtot_n[0], foutJ[0] * foutB[0], 1e-6));
      }
    }
  }

  // Check jacobtot_inv
  struct gkyl_array *jacobtot_inv_nodal =
    gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(
    n2m, &basis, &grid, &nrange_quad_interior, &range, 1, jacobtot_inv_nodal,
    gk_geom->geo_int.jacobtot_inv, true
  );
  for (int ia = nrange_quad_interior.lower[AL_IDX]; ia <= nrange_quad_interior.upper[AL_IDX];
       ++ia) {
    for (int ip = nrange_quad_interior.lower[PSI_IDX]; ip <= nrange_quad_interior.upper[PSI_IDX];
         ++ip) {
      for (int it = nrange_quad_interior.lower[TH_IDX]; it <= nrange_quad_interior.upper[TH_IDX];
           ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *jacobtot_inv_n =
          gkyl_array_fetch(jacobtot_inv_nodal, gkyl_range_idx(&nrange_quad_interior, cidx));
        // mapc2p_n[0] = x, mapc2p_n[1] = y, mapc2p_n[2] = z
        double *mapc2p_n = gkyl_array_fetch(mapc2p_nodal, gkyl_range_idx(&nrange, cidx));
        double foutJ[1];
        exact_jacobian(0.0, NULL, foutJ, 0);
        double foutB[1];
        bmag_func(0.0, NULL, foutB, 0);
        TEST_CHECK(gkyl_compare(jacobtot_inv_n[0], 1 / foutJ[0] / foutB[0], 1e-6));
      }
    }
  }

  // Check normals
  // Plus 3 away from axis to avoid errors
  struct gkyl_array *normals_nodal =
    gkyl_array_new(GKYL_DOUBLE, 9 * basis.num_basis, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(
    n2m, &basis, &grid, &nrange_quad_interior, &range, 9 * basis.num_basis, normals_nodal,
    gk_geom->geo_int.normals, true
  );
  for (int ia = nrange_quad_interior.lower[AL_IDX]; ia <= nrange_quad_interior.upper[AL_IDX];
       ++ia) {
    for (int ip = nrange_quad_interior.lower[PSI_IDX] + 3;
         ip <= nrange_quad_interior.upper[PSI_IDX]; ++ip) {
      for (int it = nrange_quad_interior.lower[TH_IDX]; it <= nrange_quad_interior.upper[TH_IDX];
           ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *normals_n =
          gkyl_array_fetch(normals_nodal, gkyl_range_idx(&nrange_quad_interior, cidx));
        double psi =
          calc_running_coord(psi_lo, ia - nrange_quad_interior.lower[PSI_IDX], grid.dx[PSI_IDX]);
        double alpha =
          calc_running_coord(alpha_lo, ia - nrange_quad_interior.lower[AL_IDX], grid.dx[AL_IDX]);
        double theta =
          calc_running_coord(theta_lo, ia - nrange_quad_interior.lower[TH_IDX], grid.dx[TH_IDX]);
        double xn[3] = {psi, alpha, theta};
        double fout[9];
        exact_normals(0.0, xn, fout, 0);
        for (int i = 0; i < 9; ++i) {
          TEST_CHECK(gkyl_compare(normals_n[i], fout[i], 1e-3));
        }

        // geo_int.normals_nodal
        double *normals = gkyl_array_fetch(
          gk_geom->geo_int.normals_nodal, gkyl_range_idx(&nrange_quad_interior, cidx)
        );
        for (int i = 0; i < 9; ++i) {
          TEST_CHECK(gkyl_compare(normals[i], fout[i], 1e-3));
        }
      }
    }
  }

  // Finally, surface geometry must be tested

  // It's much easier to generate and find the nodes using the grid generators
  struct gkyl_rect_grid psi_grid;
  struct gkyl_array *psi = gkyl_grid_array_new_from_file(&psi_grid, ginp.filename_psi);
  // create mirror geometry for surfaces
  struct gkyl_mirror_grid_gen *mirror_grid_surf[3];
  for (int dir = 0; dir < cdim; dir++) {
    mirror_grid_surf[dir] = gkyl_mirror_grid_gen_surf_inew(&(struct gkyl_mirror_grid_gen_inp){
      .comp_grid = &grid,
      .nrange = gk_geom->nrange_surf[dir],
      .local = gk_geom->local,
      .global = gk_geom->global,
      .dir = dir,
      .position_map = pos_map,

      .R = {psi_grid.lower[0], psi_grid.upper[0]},
      .Z = {psi_grid.lower[1], psi_grid.upper[1]},

      // psi(R,Z) grid size
      .nrcells = psi_grid.cells[0] - 1, // cells and not nodes
      .nzcells = psi_grid.cells[1] - 1, // cells and not nodes

      .psiRZ = psi,
      .fl_coord = ginp.fl_coord,
      .include_axis = ginp.include_axis,
      .write_psi_cubic = false,
    });
  }

  for (int dir = 0; dir < cdim; dir++) {
    enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
    enum { R_IDX, Z_IDX, PHI_IDX }; // arrangement of cartesian coordinates
    int cidx[3];
    for (int ia = gk_geom->nrange_surf[dir].lower[AL_IDX];
         ia <= gk_geom->nrange_surf[dir].upper[AL_IDX]; ++ia) {
      for (int ip = gk_geom->nrange_surf[dir].lower[PSI_IDX];
           ip <= gk_geom->nrange_surf[dir].upper[PSI_IDX]; ++ip) {
        for (int it = gk_geom->nrange_surf[dir].lower[TH_IDX];
             it <= gk_geom->nrange_surf[dir].upper[TH_IDX]; ++it) {
          cidx[PSI_IDX] = ip;
          cidx[AL_IDX] = ia;
          cidx[TH_IDX] = it;

          // First fetch the mirror stuff at this location
          const double *mirror_rza_n = gkyl_array_cfetch(
            mirror_grid_surf[dir]->nodes_rza, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx)
          );
          const double *mirror_psi_n = gkyl_array_cfetch(
            mirror_grid_surf[dir]->nodes_psi, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx)
          );
          const struct gkyl_mirror_grid_gen_geom *mirror_geo_n = gkyl_array_cfetch(
            mirror_grid_surf[dir]->nodes_geom, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx)
          );

          const double psi = mirror_psi_n[0];
          const double alpha = mirror_rza_n[2];
          const double theta = mirror_rza_n[1];

          // Next fetch the gk_geometry nodal values at this location
          double *mc2p_n = gkyl_array_fetch(
            gk_geom->geo_surf[dir].mc2p_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx)
          );
          double *jFld_n = gkyl_array_fetch(
            gk_geom->geo_surf[dir].jacobgeo_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx)
          );
          double *bmag_n = gkyl_array_fetch(
            gk_geom->geo_surf[dir].bmag_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx)
          );
          double *gFld_n = gkyl_array_fetch(
            gk_geom->geo_surf[dir].g_ij_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx)
          );
          double *biFld_n = gkyl_array_fetch(
            gk_geom->geo_surf[dir].b_i_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx)
          );
          double *cmagFld_n = gkyl_array_fetch(
            gk_geom->geo_surf[dir].cmag_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx)
          );
          double *bcartFld_n = gkyl_array_fetch(
            gk_geom->geo_surf[dir].bcart_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx)
          );
          double *tanvecFld_n = gkyl_array_fetch(
            gk_geom->geo_surf[dir].dxdz_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx)
          );
          double *dualFld_n = gkyl_array_fetch(
            gk_geom->geo_surf[dir].dzdx_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx)
          );
          double *dualmagFld_n = gkyl_array_fetch(
            gk_geom->geo_surf[dir].dualmag_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx)
          );
          double *normFld_n = gkyl_array_fetch(
            gk_geom->geo_surf[dir].normals_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx)
          );
          double *lenr_n = gkyl_array_fetch(
            gk_geom->geo_surf[dir].lenr_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx)
          );
          double *B3_n = gkyl_array_fetch(
            gk_geom->geo_surf[dir].B3_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx)
          );
          double *curlbhat_n = gkyl_array_fetch(
            gk_geom->geo_surf[dir].curlbhat_nodal, gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx)
          );
          double *normcurlbhat_n = gkyl_array_fetch(
            gk_geom->geo_surf[dir].normcurlbhat_nodal,
            gkyl_range_idx(&gk_geom->nrange_surf[dir], cidx)
          );

          // Check mapc2p
          double fout[8];
          mapc2p(0.0, (double[]){psi, -alpha, theta}, fout, 0);
          for (int i = 0; i < 3; ++i) {
            TEST_CHECK(gkyl_compare(fout[i], mc2p_n[i], 1e-8));
          }

          // Check Jacobian
          double jgeo;
          exact_jacobian(0, (double[]){psi, -alpha, theta}, &jgeo, 0);
          TEST_CHECK(gkyl_compare(jgeo, jFld_n[0], 1e-6));

          // Check bmag
          double bmag_anal[1];
          bmag_func(0, (double[]){psi, -alpha, theta}, bmag_anal, 0);
          TEST_CHECK(gkyl_compare(bmag_n[0], bmag_anal[0], 1e-7));

          // Check b_i (magnetic unit vector should have magnitude 1)
          double bi_mag =
            sqrt(biFld_n[0] * biFld_n[0] + biFld_n[1] * biFld_n[1] + biFld_n[2] * biFld_n[2]);
          TEST_CHECK(gkyl_compare(bi_mag, 1.0, 1e-8));

          // Check cmag (curvature magnitude should be 1 for straight cylinder)
          TEST_CHECK(gkyl_compare(cmagFld_n[0], 1.0, 1e-8));

          // Check normals (surface normal vectors)
          double norm_mag = sqrt(
            normFld_n[0] * normFld_n[0] + normFld_n[1] * normFld_n[1] + normFld_n[2] * normFld_n[2]
          );
          TEST_CHECK(gkyl_compare(norm_mag, 1.0, 1e-8));

          // Check that dual magnitude is consistent
          double dual_mag = sqrt(
            dualFld_n[0] * dualFld_n[0] + dualFld_n[1] * dualFld_n[1] + dualFld_n[2] * dualFld_n[2]
          );
          TEST_CHECK(gkyl_compare(dual_mag, dualmagFld_n[0], 1e-8));
        }
      }
    }
  }
  gkyl_array_release(bhat_nodal);
  gkyl_array_release(dualmag_nodal);
  gkyl_array_release(bmag_nodal);
  gkyl_array_release(bmag_nodal_interior);
  gkyl_array_release(cmag_nodal);
  gkyl_array_release(gij_nodal);
  gkyl_array_release(gij_contra_nodal);
  gkyl_array_release(jacobgeo_nodal);
  gkyl_array_release(jacobgeo_inv_nodal);
  gkyl_array_release(jacobtot_nodal);
  gkyl_array_release(jacobtot_inv_nodal);
  gkyl_array_release(mapc2p_nodal);
  gkyl_array_release(mc2p_nodal_interior);
  gkyl_array_release(mc2nu_pos_nodal);
  gkyl_array_release(normals_nodal);
  gkyl_array_release(psi);
  gkyl_nodal_ops_release(n2m);
  gkyl_position_map_release(pos_map);
  for (int dir = 0; dir < cdim; dir++) {
    gkyl_mirror_grid_gen_release(mirror_grid_surf[dir]);
  }
  gkyl_gk_geometry_release(gk_geom);
}

void
mapz(double t, const double *xn, double *GKYL_RESTRICT fout, void *ctx)
{
  double a = 2.0;
  double s = 0.2;
  fout[0] = (-1 / (2 * a) * pow(a - xn[0], 2) + a) * (1 - s) + s * xn[0];
}

void
dmapz_dz(double t, const double *xn, double *GKYL_RESTRICT fout, void *ctx)
{
  double a = 2.0;
  double s = 0.2;
  fout[0] = (1 / a * pow(a - xn[0], 1)) * (1 - s) + s;
}

void
exact_gij_pmap(double t, const double *xn, double *GKYL_RESTRICT fout, void *ctx)
{
  double psi = xn[0], alpha = xn[1], theta = xn[2];
  double r = 2 * sqrt(psi);
  double dThetadtheta;
  dmapz_dz(0, &theta, &dThetadtheta, 0);
  fout[0] = 1 / psi; // g_11
  fout[1] = 0.0; // g_12
  fout[2] = 0.0; // g_13
  fout[3] = r * r; // g_22
  fout[4] = 0.0; // g_23
  fout[5] = 1.0 * pow(dThetadtheta, 2); // g_33
}

void
exact_g_contra_ij_pmap(double t, const double *xn, double *GKYL_RESTRICT fout, void *ctx)
{
  double psi = xn[0], alpha = xn[1], theta = xn[2];
  double r = 2 * sqrt(psi);
  double dThetadtheta;
  dmapz_dz(0, &theta, &dThetadtheta, 0);
  fout[0] = r * r / 4; // g_11
  fout[1] = 0.0; // g_12
  fout[2] = 0.0; // g_13
  fout[3] = 1 / psi / 4; // g_22
  fout[4] = 0.0; // g_23
  fout[5] = 1.0 / pow(dThetadtheta, 2); // g_33
}

void
exact_dual_magnitude_pmap(double t, const double *xn, double *GKYL_RESTRICT fout, void *ctx)
{
  double psi = xn[0], alpha = xn[1], theta = xn[2];
  double r = 2 * sqrt(psi);
  double dThetadtheta;
  dmapz_dz(0, &theta, &dThetadtheta, 0);
  fout[0] = r / 2;
  fout[1] = 1 / (2 * sqrt(psi));
  fout[2] = 1.0 / dThetadtheta;
}

void
exact_normals_pmap(double t, const double *xn, double *GKYL_RESTRICT fout, void *ctx)
{
  double psi = xn[0], alpha = xn[1], theta = xn[2];
  double r = 2 * sqrt(psi);
  double dThetadtheta;
  dmapz_dz(0, &theta, &dThetadtheta, 0);
  // Remember cylindrical angle = - alpha
  fout[0] = -cos(-alpha);
  fout[1] = -sin(-alpha);
  fout[2] = 0.0;
  fout[3] = -sin(-alpha);
  fout[4] = cos(-alpha);
  fout[5] = 0.0;
  fout[6] = 0.0;
  fout[7] = 0.0;
  fout[8] = -1.0;
}

void
test_mirror_3x_p1_pmap_straight_cylinder_ho(void)
{
  // Same as the above test, but using a quadratic position map
  struct gkyl_basis basis;
  int poly_order = 1;
  int cdim = 3;
  gkyl_cart_modal_serendip(&basis, cdim, poly_order);

  double psiMax = 0.2;
  double psiMin = 0.1;
  int Nz = 10;

  double lower[3] = {psiMin, -M_PI, -2.0};
  double upper[3] = {psiMax, M_PI, 2.0};

  // int cells[3] = { 18, 18, Nz };
  int cells[3] = {8, 1, 8};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);

  struct gkyl_range ext_range, range;
  int nghost[3] = {1, 1, 1};
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  struct gkyl_position_map_inp pos_map_inp = {.maps = {0, 0, mapz}, .ctxs = {0, 0, 0}};

  // Configuration space geometry initialization
  struct gkyl_position_map *pos_map =
    gkyl_position_map_new(pos_map_inp, grid, range, ext_range, range, ext_range, basis);

  struct gkyl_mirror_geo_grid_inp ginp = {
    .filename_psi = "gyrokinetic/data/unit/wham_hires.geqdsk_psi.gkyl", // psi file to use
    .rclose = 0.2, // closest R to region of interest
    .zmin = -2.0, // Z of lower boundary
    .zmax = 2.0, // Z of upper boundary
    .include_axis = false, // Include R=0 axis in grid
    .fl_coord = GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z, // coordinate system for psi grid
  };

  // Initialize geometry
  struct gkyl_gk_geometry_inp geometry_input = {
    .geometry_id = GKYL_GEOMETRY_MIRROR,
    .mirror_grid_info = ginp,
    .position_map = pos_map,
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

  struct gk_geometry *gk_geom = gkyl_gk_geometry_mirror_new(&geometry_input);

  gkyl_position_map_set_mc2nu(pos_map, gk_geom->geo_corn.mc2nu_pos);

  // write_geometry(gk_geom, grid, range, "straight_cylinder");

  int theta_shift = 1; // Because of the forward/backward difference
  // used to calculate derivatives of the map, the edge values of the map
  // have a linear order error so the derivatives are not accurate and the
  // Jacobian will be a little wrong there.

  // Define nodal operations
  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  int cidx[3];
  int nodes[] = {1, 1, 1};
  for (int d = 0; d < grid.ndim; ++d) {
    nodes[d] = grid.cells[d] + 1;
  }
  struct gkyl_range nrange;
  gkyl_range_init_from_shape(&nrange, grid.ndim, nodes);

  int nodes_quad_interior[] = {1, 1, 1};
  int num_quad_points = poly_order + 1;
  for (int d = 0; d < grid.ndim; ++d) {
    nodes_quad_interior[d] = grid.cells[d] * num_quad_points;
  }
  struct gkyl_range nrange_quad_interior;
  gkyl_range_init_from_shape(&nrange_quad_interior, grid.ndim, nodes_quad_interior);

  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&basis, &grid, false);

  double dels[2] = {1.0 / sqrt(3), 1.0 - 1.0 / sqrt(3)};
  double theta_lo = grid.lower[TH_IDX] + dels[1] * grid.dx[TH_IDX] / 2.0;
  double psi_lo = grid.lower[PSI_IDX] + dels[1] * grid.dx[PSI_IDX] / 2.0;
  double alpha_lo = grid.lower[AL_IDX] + dels[1] * grid.dx[AL_IDX] / 2.0;

  struct gkyl_array *mapc2p_nodal = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
  gkyl_nodal_ops_m2n(
    n2m, &basis, &grid, &nrange, &range, 3, mapc2p_nodal, gk_geom->geo_corn.mc2p, false
  );

  // Check that |bhat|=1 at nodes
  struct gkyl_array *bhat_nodal =
    gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(
    n2m, &basis, &grid, &nrange_quad_interior, &range, 3, bhat_nodal, gk_geom->geo_int.bcart, true
  );
  for (int ia = nrange_quad_interior.lower[AL_IDX]; ia <= nrange_quad_interior.upper[AL_IDX];
       ++ia) {
    for (int ip = nrange_quad_interior.lower[PSI_IDX]; ip <= nrange_quad_interior.upper[PSI_IDX];
         ++ip) {
      for (int it = nrange_quad_interior.lower[TH_IDX]; it <= nrange_quad_interior.upper[TH_IDX];
           ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *bhat_n = gkyl_array_fetch(bhat_nodal, gkyl_range_idx(&nrange_quad_interior, cidx));
        double bhat_mag =
          sqrt(bhat_n[0] * bhat_n[0] + bhat_n[1] * bhat_n[1] + bhat_n[2] * bhat_n[2]);
        TEST_CHECK(gkyl_compare(bhat_mag, 1.0, 1e-12));
      }
    }
  }

  // Check that the duals are what they should be
  // There are errors at low psi, so we shift away from the axis by a few cells
  struct gkyl_array *dualmag_nodal =
    gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(
    n2m, &basis, &grid, &nrange_quad_interior, &range, 3, dualmag_nodal, gk_geom->geo_int.dualmag,
    true
  );
  struct gkyl_array *mapc2p_nodal_interior =
    gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(
    n2m, &basis, &grid, &nrange_quad_interior, &range, 3, mapc2p_nodal_interior,
    gk_geom->geo_int.mc2p, true
  );
  for (int ia = nrange_quad_interior.lower[AL_IDX]; ia <= nrange_quad_interior.upper[AL_IDX];
       ++ia) {
    for (int ip = nrange_quad_interior.lower[PSI_IDX] + 3;
         ip <= nrange_quad_interior.upper[PSI_IDX]; ++ip) {
      for (int it = nrange_quad_interior.lower[TH_IDX] + theta_shift;
           it <= nrange_quad_interior.upper[TH_IDX] - theta_shift; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double psi =
          calc_running_coord(psi_lo, ia - nrange_quad_interior.lower[PSI_IDX], grid.dx[PSI_IDX]);
        double alpha =
          calc_running_coord(alpha_lo, ia - nrange_quad_interior.lower[AL_IDX], grid.dx[AL_IDX]);
        double theta =
          calc_running_coord(theta_lo, ia - nrange_quad_interior.lower[TH_IDX], grid.dx[TH_IDX]);
        double xn[3] = {psi, alpha, theta};
        double *dualmag_n =
          gkyl_array_fetch(dualmag_nodal, gkyl_range_idx(&nrange_quad_interior, cidx));
        double dualmag_anal[3];
        exact_dual_magnitude_pmap(0, xn, dualmag_anal, 0);
        TEST_CHECK(gkyl_compare(dualmag_n[0], dualmag_anal[0], 1e-6));
        TEST_CHECK(gkyl_compare(dualmag_n[1], dualmag_anal[1], 1e-6));
        TEST_CHECK(gkyl_compare(dualmag_n[2], dualmag_anal[2], 1e-6));
      }
    }
  }

  // Check bmag is what it should be
  struct gkyl_array *bmag_nodal = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
  gkyl_nodal_ops_m2n(
    n2m, &basis, &grid, &nrange, &range, 1, bmag_nodal, gk_geom->geo_corn.bmag, false
  );
  for (int ia = nrange.lower[AL_IDX]; ia <= nrange.upper[AL_IDX]; ++ia) {
    for (int ip = nrange.lower[PSI_IDX]; ip <= nrange.upper[PSI_IDX]; ++ip) {
      for (int it = nrange.lower[TH_IDX]; it <= nrange.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double psi = grid.lower[PSI_IDX] +
                     ip * (grid.upper[PSI_IDX] - grid.lower[PSI_IDX]) / grid.cells[PSI_IDX];
        double alpha =
          grid.lower[AL_IDX] + ia * (grid.upper[AL_IDX] - grid.lower[AL_IDX]) / grid.cells[AL_IDX];
        double theta =
          grid.lower[TH_IDX] + it * (grid.upper[TH_IDX] - grid.lower[TH_IDX]) / grid.cells[TH_IDX];
        double xn[3] = {psi, alpha, theta};
        double *bmag_n = gkyl_array_fetch(bmag_nodal, gkyl_range_idx(&nrange, cidx));
        double bmag_anal[1];
        bmag_func(0, xn, bmag_anal, 0);
        TEST_CHECK(gkyl_compare(bmag_n[0], bmag_anal[0], 1e-8));
      }
    }
  }

  // Check cmag = 1
  struct gkyl_array *cmag_nodal = gkyl_array_new(GKYL_DOUBLE, 1, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(
    n2m, &basis, &grid, &nrange_quad_interior, &range, 1, cmag_nodal, gk_geom->geo_int.cmag, true
  );
  for (int ia = nrange_quad_interior.lower[AL_IDX]; ia <= nrange_quad_interior.upper[AL_IDX];
       ++ia) {
    for (int ip = nrange_quad_interior.lower[PSI_IDX]; ip <= nrange_quad_interior.upper[PSI_IDX];
         ++ip) {
      for (int it = nrange_quad_interior.lower[TH_IDX]; it <= nrange_quad_interior.upper[TH_IDX];
           ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *cmag_n = gkyl_array_fetch(cmag_nodal, gkyl_range_idx(&nrange_quad_interior, cidx));
        TEST_CHECK(gkyl_compare(cmag_n[0], 1.0, 1e-6));
      }
    }
  }

  // Check g_ij
  struct gkyl_array *gij_nodal = gkyl_array_new(GKYL_DOUBLE, 6, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(
    n2m, &basis, &grid, &nrange_quad_interior, &range, 6, gij_nodal, gk_geom->geo_int.g_ij, true
  );
  for (int ia = nrange_quad_interior.lower[AL_IDX]; ia <= nrange_quad_interior.upper[AL_IDX];
       ++ia) {
    for (int ip = nrange_quad_interior.lower[PSI_IDX]; ip <= nrange_quad_interior.upper[PSI_IDX];
         ++ip) {
      for (int it = nrange_quad_interior.lower[TH_IDX] + theta_shift;
           it <= nrange_quad_interior.upper[TH_IDX] - theta_shift; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double psi =
          calc_running_coord(psi_lo, ia - nrange_quad_interior.lower[PSI_IDX], grid.dx[PSI_IDX]);
        double alpha =
          calc_running_coord(alpha_lo, ia - nrange_quad_interior.lower[AL_IDX], grid.dx[AL_IDX]);
        double theta =
          calc_running_coord(theta_lo, ia - nrange_quad_interior.lower[TH_IDX], grid.dx[TH_IDX]);
        double *gij_n = gkyl_array_fetch(gij_nodal, gkyl_range_idx(&nrange_quad_interior, cidx));
        double xn[3] = {psi, alpha, theta};
        double fout[6];
        exact_gij_pmap(0.0, xn, fout, 0);
        for (int i = 0; i < 6; ++i) {
          TEST_CHECK(gkyl_compare(gij_n[i], fout[i], 1e-6));
        }
      }
    }
  }

  // Check g^ij
  struct gkyl_array *gij_contra_nodal = gkyl_array_new(GKYL_DOUBLE, 6, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(
    n2m, &basis, &grid, &nrange_quad_interior, &range, 6, gij_contra_nodal, gk_geom->geo_int.gij,
    true
  );
  for (int ia = nrange_quad_interior.lower[AL_IDX]; ia <= nrange_quad_interior.upper[AL_IDX];
       ++ia) {
    for (int ip = nrange_quad_interior.lower[PSI_IDX]; ip <= nrange_quad_interior.upper[PSI_IDX];
         ++ip) {
      for (int it = nrange_quad_interior.lower[TH_IDX] + theta_shift;
           it <= nrange_quad_interior.upper[TH_IDX] - theta_shift; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *gij_contra_n =
          gkyl_array_fetch(gij_contra_nodal, gkyl_range_idx(&nrange_quad_interior, cidx));
        double psi =
          calc_running_coord(psi_lo, ia - nrange_quad_interior.lower[PSI_IDX], grid.dx[PSI_IDX]);
        double alpha =
          calc_running_coord(alpha_lo, ia - nrange_quad_interior.lower[AL_IDX], grid.dx[AL_IDX]);
        double theta =
          calc_running_coord(theta_lo, ia - nrange_quad_interior.lower[TH_IDX], grid.dx[TH_IDX]);
        double xn[3] = {psi, alpha, theta};
        double fout[6];
        exact_g_contra_ij_pmap(0.0, xn, fout, 0);
        for (int i = 0; i < 6; ++i) {
          TEST_CHECK(gkyl_compare(gij_contra_n[i], fout[i], 1e-6));
        }
      }
    }
  }

  // Check that Jacobgeo is what it should be. This is the Jacobian for the problem
  struct gkyl_array *jacobgeo_nodal =
    gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(
    n2m, &basis, &grid, &nrange_quad_interior, &range, 1, jacobgeo_nodal, gk_geom->geo_int.jacobgeo,
    true
  );
  for (int ia = nrange_quad_interior.lower[AL_IDX]; ia <= nrange_quad_interior.upper[AL_IDX];
       ++ia) {
    for (int ip = nrange_quad_interior.lower[PSI_IDX]; ip <= nrange_quad_interior.upper[PSI_IDX];
         ++ip) {
      for (int it = nrange_quad_interior.lower[TH_IDX] + theta_shift;
           it <= nrange_quad_interior.upper[TH_IDX] - theta_shift; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *jacobgeo_n =
          gkyl_array_fetch(jacobgeo_nodal, gkyl_range_idx(&nrange_quad_interior, cidx));
        double theta =
          grid.lower[TH_IDX] + it * (grid.upper[TH_IDX] - grid.lower[TH_IDX]) / grid.cells[TH_IDX];
        double xn[1] = {theta};
        double fout[1];
        dmapz_dz(0.0, xn, fout, 0);
        double jacobian_analytic = 2 / M_PI * fout[0];
        TEST_CHECK(gkyl_compare(jacobgeo_n[0], jacobian_analytic, 1e-6));
      }
    }
  }

  // Check jacobgeo_inv
  struct gkyl_array *jacobgeo_inv_nodal =
    gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(
    n2m, &basis, &grid, &nrange_quad_interior, &range, 1, jacobgeo_inv_nodal,
    gk_geom->geo_int.jacobgeo_inv, true
  );
  for (int ia = nrange_quad_interior.lower[AL_IDX]; ia <= nrange_quad_interior.upper[AL_IDX];
       ++ia) {
    for (int ip = nrange_quad_interior.lower[PSI_IDX]; ip <= nrange_quad_interior.upper[PSI_IDX];
         ++ip) {
      for (int it = nrange_quad_interior.lower[TH_IDX] + theta_shift;
           it <= nrange_quad_interior.upper[TH_IDX] - theta_shift; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *jacobgeo_inv_n =
          gkyl_array_fetch(jacobgeo_inv_nodal, gkyl_range_idx(&nrange_quad_interior, cidx));
        double theta =
          grid.lower[TH_IDX] + it * (grid.upper[TH_IDX] - grid.lower[TH_IDX]) / grid.cells[TH_IDX];
        double xn[1] = {theta};
        double fout[1];
        dmapz_dz(0.0, xn, fout, 0);
        double jacobian_analytic = 2 / M_PI * fout[0];
        TEST_CHECK(gkyl_compare(jacobgeo_inv_n[0], 1 / jacobian_analytic, 1e-6));
      }
    }
  }

  // Check jacobtot
  struct gkyl_array *jacobtot_nodal =
    gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(
    n2m, &basis, &grid, &nrange_quad_interior, &range, 1, jacobtot_nodal, gk_geom->geo_int.jacobtot,
    true
  );
  for (int ia = nrange_quad_interior.lower[AL_IDX]; ia <= nrange_quad_interior.upper[AL_IDX];
       ++ia) {
    for (int ip = nrange_quad_interior.lower[PSI_IDX]; ip <= nrange_quad_interior.upper[PSI_IDX];
         ++ip) {
      for (int it = nrange_quad_interior.lower[TH_IDX] + theta_shift;
           it <= nrange_quad_interior.upper[TH_IDX] - theta_shift; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *jacobtot_n =
          gkyl_array_fetch(jacobtot_nodal, gkyl_range_idx(&nrange_quad_interior, cidx));
        double theta =
          grid.lower[TH_IDX] + it * (grid.upper[TH_IDX] - grid.lower[TH_IDX]) / grid.cells[TH_IDX];
        double xn[1] = {theta};
        double fout[1];
        dmapz_dz(0.0, xn, fout, 0);
        double jacobian_analytic = 2 / M_PI * fout[0];
        double magnetic_field = 0.5;
        double jacobtot_analytic = jacobian_analytic * magnetic_field;
        TEST_CHECK(gkyl_compare(jacobtot_n[0], jacobtot_analytic, 1e-6));
      }
    }
  }

  // Check jacobtot_inv
  struct gkyl_array *jacobtot_inv_nodal =
    gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(
    n2m, &basis, &grid, &nrange_quad_interior, &range, 1, jacobtot_inv_nodal,
    gk_geom->geo_int.jacobtot_inv, true
  );
  for (int ia = nrange_quad_interior.lower[AL_IDX]; ia <= nrange_quad_interior.upper[AL_IDX];
       ++ia) {
    for (int ip = nrange_quad_interior.lower[PSI_IDX]; ip <= nrange_quad_interior.upper[PSI_IDX];
         ++ip) {
      for (int it = nrange_quad_interior.lower[TH_IDX] + theta_shift;
           it <= nrange_quad_interior.upper[TH_IDX] - theta_shift; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *jacobtot_inv_n =
          gkyl_array_fetch(jacobtot_inv_nodal, gkyl_range_idx(&nrange_quad_interior, cidx));
        double theta =
          grid.lower[TH_IDX] + it * (grid.upper[TH_IDX] - grid.lower[TH_IDX]) / grid.cells[TH_IDX];
        double xn[1] = {theta};
        double fout[1];
        dmapz_dz(0.0, xn, fout, 0);
        double jacobian_analytic = 2 / M_PI * fout[0];
        double magnetic_field = 0.5;
        double jacobtot_analytic = jacobian_analytic * magnetic_field;
        TEST_CHECK(gkyl_compare(jacobtot_inv_n[0], 1 / jacobtot_analytic, 1e-6));
      }
    }
  }

  // Check mapc2p
  for (int ia = nrange.lower[AL_IDX]; ia <= nrange.upper[AL_IDX]; ++ia) {
    for (int ip = nrange.lower[PSI_IDX]; ip <= nrange.upper[PSI_IDX]; ++ip) {
      for (int it = nrange.lower[TH_IDX]; it <= nrange.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double psi =
          calc_running_coord(psi_lo, ia - nrange_quad_interior.lower[PSI_IDX], grid.dx[PSI_IDX]);
        double alpha =
          calc_running_coord(alpha_lo, ia - nrange_quad_interior.lower[AL_IDX], grid.dx[AL_IDX]);
        double theta =
          calc_running_coord(theta_lo, ia - nrange_quad_interior.lower[TH_IDX], grid.dx[TH_IDX]);
        // mapc2p_n[0] = R, mapc2p_n[1] = Theta, mapc2p_n[2] = Z_cylindrical
        double *mapc2p_n = gkyl_array_fetch(mapc2p_nodal, gkyl_range_idx(&nrange, cidx));
        mapz(0.0, &theta, &theta, 0);
        double xn[3] = {psi, alpha, theta};
        double fout[3];
        mapc2p(0.0, xn, fout, 0);
        for (int i = 0; i < 3; ++i) {
          TEST_CHECK(gkyl_compare(mapc2p_n[i], fout[i], 1e-8));
        }
      }
    }
  }

  // Check mc2nu_pos
  struct gkyl_array *mc2nu_pos_nodal = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
  gkyl_nodal_ops_m2n(
    n2m, &basis, &grid, &nrange, &range, 3, mc2nu_pos_nodal, gk_geom->geo_corn.mc2nu_pos, false
  );
  for (int ia = nrange.lower[AL_IDX]; ia <= nrange.upper[AL_IDX]; ++ia) {
    for (int ip = nrange.lower[PSI_IDX]; ip <= nrange.upper[PSI_IDX]; ++ip) {
      for (int it = nrange.lower[TH_IDX]; it <= nrange.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double psi =
          calc_running_coord(psi_lo, ia - nrange_quad_interior.lower[PSI_IDX], grid.dx[PSI_IDX]);
        double alpha =
          calc_running_coord(alpha_lo, ia - nrange_quad_interior.lower[AL_IDX], grid.dx[AL_IDX]);
        double theta =
          calc_running_coord(theta_lo, ia - nrange_quad_interior.lower[TH_IDX], grid.dx[TH_IDX]);
        mapz(0.0, &theta, &theta, 0);
        double xn[3] = {psi, alpha, theta};
        double *mc2nu_pos_n = gkyl_array_fetch(mc2nu_pos_nodal, gkyl_range_idx(&nrange, cidx));
        for (int i = 0; i < 3; ++i) {
          TEST_CHECK(gkyl_compare(mc2nu_pos_n[i], xn[i], 1e-8));
        }
      }
    }
  }

  // Check normals
  // Plus 3 away from axis to avoid errors
  struct gkyl_array *normals_nodal =
    gkyl_array_new(GKYL_DOUBLE, 9 * basis.num_basis, nrange_quad_interior.volume);
  gkyl_nodal_ops_m2n(
    n2m, &basis, &grid, &nrange_quad_interior, &range, 9 * basis.num_basis, normals_nodal,
    gk_geom->geo_int.normals, true
  );
  for (int ia = nrange_quad_interior.lower[AL_IDX]; ia <= nrange_quad_interior.upper[AL_IDX];
       ++ia) {
    for (int ip = nrange_quad_interior.lower[PSI_IDX] + 3;
         ip <= nrange_quad_interior.upper[PSI_IDX]; ++ip) {
      for (int it = nrange_quad_interior.lower[TH_IDX]; it <= nrange_quad_interior.upper[TH_IDX];
           ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *normals_n =
          gkyl_array_fetch(normals_nodal, gkyl_range_idx(&nrange_quad_interior, cidx));
        double psi =
          calc_running_coord(psi_lo, ia - nrange_quad_interior.lower[PSI_IDX], grid.dx[PSI_IDX]);
        double alpha =
          calc_running_coord(alpha_lo, ia - nrange_quad_interior.lower[AL_IDX], grid.dx[AL_IDX]);
        double theta =
          calc_running_coord(theta_lo, ia - nrange_quad_interior.lower[TH_IDX], grid.dx[TH_IDX]);
        double xn[3] = {psi, alpha, theta};
        double fout[9];
        exact_normals_pmap(0.0, xn, fout, 0);
        for (int i = 0; i < 9; ++i) {
          TEST_CHECK(gkyl_compare(normals_n[i], fout[i], 1e-6));
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
  gkyl_position_map_release(pos_map);
  gkyl_gk_geometry_release(gk_geom);
}

struct mirror_direction_map {
  double shift, scale, curvature;
};

static void
mirror_direction_map(double t, const double *xn, double *out, void *ctx)
{
  const struct mirror_direction_map *map = ctx;
  out[0] = map->shift + xn[0] * (map->scale + map->curvature * xn[0]);
}

static void
mirror_direction_map_deriv(double t, const double *xn, double *out, void *ctx)
{
  const struct mirror_direction_map *map = ctx;
  out[0] = map->scale + 2.0 * map->curvature * xn[0];
}

// Check the assembled geometry, including the modal fields after the generic
// derived-geometry updater has run. That updater assumes bhat is along +e_3.
static void
check_mirror_field_direction(double curvature)
{
  const double strength = 2.7;
  struct gkyl_rect_grid psi_grid;
  gkyl_rect_grid_init(&psi_grid, 2, (double[]){0.0, -1.0}, (double[]){1.0, 1.0}, (int[]){9, 17});
  struct gkyl_range psi_range;
  gkyl_range_init_from_shape(&psi_range, 2, psi_grid.cells);
  struct gkyl_array *psi = gkyl_array_new(GKYL_DOUBLE, 1, psi_range.volume);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &psi_range);
  while (gkyl_range_iter_next(&iter)) {
    double radius = iter.idx[0] / 8.0, height = -1.0 + iter.idx[1] / 8.0;
    double *value = gkyl_array_fetch(psi, gkyl_range_idx(&psi_range, iter.idx));
    value[0] = 0.5 * strength * radius * radius * (1.0 + curvature * height * height);
  }
  char filename[] = "/tmp/ctest_gk_geometry_mirror_XXXXXX";
  int fd = mkstemp(filename);
  TEST_ASSERT(fd >= 0);
  close(fd);
  TEST_ASSERT(
    gkyl_grid_sub_array_write(&psi_grid, &psi_range, 0, psi, filename) == GKYL_ARRAY_RIO_SUCCESS
  );
  gkyl_array_release(psi);

  for (int sqrt_psi = 0; sqrt_psi < 2; ++sqrt_psi) {
    for (int mapped = 0; mapped < 2; ++mapped) {
      struct gkyl_rect_grid grid;
      gkyl_rect_grid_init(
        &grid, 3, (double[]){0.04, -0.4, -0.6}, (double[]){0.22, 0.4, 0.6}, (int[]){2, 2, 3}
      );
      struct gkyl_range local, ext;
      gkyl_create_grid_ranges(&grid, (int[]){1, 1, 1}, &ext, &local);
      struct gkyl_basis basis;
      gkyl_cart_modal_serendip(&basis, 3, 1);
      struct mirror_direction_map maps[3] = {
        {0.01, 0.8, 0.07}, {0.08, 1.1, 0.08}, {-0.03, 0.9, 0.13}
      };
      struct gkyl_position_map *pmap = gkyl_position_map_null_new();
      for (int dim = 0; dim < 3; ++dim) {
        if (!mapped) {
          maps[dim] = (struct mirror_direction_map){0.0, 1.0, 0.0};
        }
        pmap->maps[dim] = mirror_direction_map;
        pmap->map_derivs[dim] = mirror_direction_map_deriv;
        pmap->ctxs[dim] = &maps[dim];
      }
      pmap->use_map_derivs = true;
      struct gkyl_gk_geometry_inp inp = {
        .geometry_id = GKYL_GEOMETRY_MIRROR,
        .mirror_grid_info =
          {
            .fl_coord = sqrt_psi ? GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z :
                                   GKYL_GEOMETRY_MIRROR_GRID_GEN_PSI_CART_Z,
          },
        .position_map = pmap,
        .grid = grid,
        .local = local,
        .local_ext = ext,
        .global = local,
        .global_ext = ext,
        .basis = basis,
        .geo_grid = grid,
        .geo_local = local,
        .geo_local_ext = ext,
        .geo_global = local,
        .geo_global_ext = ext,
        .geo_basis = basis,
      };
      snprintf(
        inp.mirror_grid_info.filename_psi, sizeof inp.mirror_grid_info.filename_psi, "%s", filename
      );
      struct gk_geometry *geom = gkyl_gk_geometry_mirror_new(&inp);
      TEST_ASSERT(geom != NULL);
      struct gkyl_nodal_ops *ops = gkyl_nodal_ops_new(&basis, &grid, false);

      for (int kind = -1; kind < 3; ++kind) {
        bool surface = kind >= 0;
        struct gk_geom_int *vol = &geom->geo_int;
        struct gk_geom_surf *face = surface ? &geom->geo_surf[kind] : NULL;
        const struct gkyl_range *nodes = surface ? &geom->nrange_surf[kind] : &geom->nrange_int;
        enum { BMAG, BCART, BI, B3, CURL, CURL_PROJ, CURL_OVER_B, BPAR, BI_OVER_JB, CMAG, NFIELDS };
        struct {
          const char *name;
          int ncomp;
          const struct gkyl_array *nodal, *modal;
          struct gkyl_array *recovered;
        } fields[NFIELDS] = {
          {"bmag", 1, surface ? face->bmag_nodal : vol->bmag_nodal, surface ? face->bmag : vol->bmag
          },
          {"bcart", 3, surface ? face->bcart_nodal : vol->bcart_nodal, surface ? NULL : vol->bcart},
          {"b_i", 3, surface ? face->b_i_nodal : vol->b_i_nodal, surface ? face->b_i : vol->b_i},
          {"B3", 1, surface ? face->B3_nodal : vol->B3_nodal, surface ? face->B3 : vol->B3},
          {"curlbhat", 3, surface ? face->curlbhat_nodal : vol->curlbhat_nodal, NULL},
          {"projected curl", surface ? 1 : 3,
           surface ? face->normcurlbhat_nodal : vol->dualcurlbhat_nodal,
           surface ? face->normcurlbhat : vol->dualcurlbhat},
          {"dualcurlbhatoverB", 3, surface ? NULL : vol->dualcurlbhatoverB_nodal,
           surface ? NULL : vol->dualcurlbhatoverB},
          {"rtg33inv", 1, surface ? NULL : vol->rtg33inv_nodal, surface ? NULL : vol->rtg33inv},
          {"bioverJB", 3, surface ? NULL : vol->bioverJB_nodal, surface ? NULL : vol->bioverJB},
          {"cmag", 1, surface ? face->cmag_nodal : NULL, surface ? face->cmag : vol->cmag}
        };
        struct gkyl_range update = local;
        if (surface) {
          int upper[3];
          gkyl_copy_int_arr(3, local.upper, upper);
          upper[kind]++;
          gkyl_sub_range_init(&update, &ext, local.lower, upper);
        }
        for (int field = 0; field < NFIELDS; ++field) {
          if (fields[field].modal) {
            fields[field].recovered =
              gkyl_array_new(GKYL_DOUBLE, fields[field].ncomp, nodes->volume);
            if (surface) {
              gkyl_nodal_ops_m2n_surface(
                ops, &geom->surf_basis, &grid, nodes, &update, fields[field].ncomp,
                fields[field].recovered, fields[field].modal, kind
              );
            } else {
              gkyl_nodal_ops_m2n(
                ops, &basis, &grid, nodes, &update, fields[field].ncomp, fields[field].recovered,
                fields[field].modal, true
              );
            }
          }
        }

        gkyl_range_iter_init(&iter, nodes);
        while (gkyl_range_iter_next(&iter)) {
          double mapped_x[3], slope[3];
          for (int dim = 0; dim < 3; ++dim) {
            double lower = grid.lower[dim], upper = grid.upper[dim];
            if (dim == 0 && sqrt_psi) {
              lower = sqrt(lower);
              upper = sqrt(upper);
            }
            int node = iter.idx[dim] - nodes->lower[dim];
            double offset =
              dim == kind ? node : node / 2 + 0.5 * (1.0 + (node % 2 ? 1.0 : -1.0) / sqrt(3.0));
            double xc = lower + (upper - lower) * offset / grid.cells[dim];
            mirror_direction_map(0.0, &xc, &mapped_x[dim], &maps[dim]);
            mirror_direction_map_deriv(0.0, &xc, &slope[dim], &maps[dim]);
          }
          double flux = sqrt_psi ? mapped_x[0] * mapped_x[0] : mapped_x[0];
          double height = mapped_x[2], phi = mapped_x[1];
          double factor = 1.0 + curvature * height * height;
          double radius = sqrt(2.0 * flux / (strength * factor));
          double br = -curvature * radius * height, norm = sqrt(br * br + factor * factor);
          double radial_scale = slope[0] * (sqrt_psi ? 2.0 * mapped_x[0] : 1.0);
          double jac = radial_scale * slope[1] * slope[2] / (strength * factor);
          double curl_phi = -curvature * radius * factor *
                            (1.0 - 2.0 * curvature * height * height) / (norm * norm * norm);
          double expected[NFIELDS][3] = {0};
          expected[BMAG][0] = strength * norm;
          expected[BCART][0] = br * cos(phi) / norm;
          expected[BCART][1] = br * sin(phi) / norm;
          expected[BCART][2] = factor / norm;
          expected[BI][0] = -radial_scale * curvature * height / (strength * factor * norm);
          expected[BI][2] = slope[2] * norm / factor;
          expected[B3][0] = strength * factor / (surface ? 1.0 : slope[2]);
          expected[CURL][0] = -curl_phi * sin(phi);
          expected[CURL][1] = curl_phi * cos(phi);
          if (surface) {
            expected[CURL_PROJ][0] = kind == 1 ? curl_phi : 0.0;
          } else {
            expected[CURL_PROJ][1] = curl_phi / (radius * slope[1]);
            expected[CURL_OVER_B][1] = expected[CURL_PROJ][1] / expected[BMAG][0];
            expected[BPAR][0] = factor / (slope[2] * norm);
            for (int dim = 0; dim < 3; ++dim) {
              expected[BI_OVER_JB][dim] = expected[BI][dim] / (jac * expected[BMAG][0]);
            }
          }
          expected[CMAG][0] = radial_scale * slope[1];
          long loc = gkyl_range_idx(nodes, iter.idx);
          for (int field = 0; field < NFIELDS; ++field) {
            const struct gkyl_array *representations[] = {
              fields[field].nodal, fields[field].recovered
            };
            for (int rep = 0; rep < 2; ++rep) {
              if (!representations[rep]) {
                continue;
              }
              const double *actual = gkyl_array_cfetch(representations[rep], loc);
              for (int component = 0; component < fields[field].ncomp; ++component) {
                double reference = expected[field][component];
                TEST_CHECK(fabs(actual[component] - reference) < 2e-9 * (1.0 + fabs(reference)));
                TEST_MSG(
                  "%s[%d], kind=%d mapped=%d sqrt_psi=%d rep=%d: %.17g != %.17g",
                  fields[field].name, component, kind, mapped, sqrt_psi, rep, actual[component],
                  reference
                );
              }
            }
          }
        }
        for (int field = 0; field < NFIELDS; ++field) {
          if (fields[field].recovered) {
            gkyl_array_release(fields[field].recovered);
          }
        }
      }
      gkyl_nodal_ops_release(ops);
      gkyl_gk_geometry_release(geom);
      gkyl_position_map_release(pmap);
    }
  }
  TEST_CHECK(remove(filename) == 0);
}

static void
test_mirror_straight_field_direction_ho(void)
{
  check_mirror_field_direction(0.0);
}

static void
test_mirror_curved_field_direction_ho(void)
{
  check_mirror_field_direction(0.7);
}

TEST_LIST = {
  {"test_mirror_load_geometry_ho", test_mirror_load_geometry_ho},
  {"test_mirror_3x_p1_straight_cylinder_ho", test_mirror_3x_p1_straight_cylinder_ho},
  {"test_mirror_straight_field_direction_ho", test_mirror_straight_field_direction_ho},
  {"test_mirror_curved_field_direction_ho", test_mirror_curved_field_direction_ho},
  // { "test_mirror_3x_p1_pmap_straight_cylinder_ho", test_mirror_3x_p1_pmap_straight_cylinder_ho },
  {NULL, NULL}
};
