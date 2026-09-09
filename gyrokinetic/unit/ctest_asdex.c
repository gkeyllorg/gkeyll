#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>


#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_ops.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_position_map.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <gkyl_basis.h>

#include <gkyl_efit.h>
#include <gkyl_calc_bmag.h>
#include <gkyl_tok_geo.h>

#include <gkyl_calc_metric.h>
#include <gkyl_calc_derived_geo.h>

#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_tok.h>
#include <gkyl_nodal_ops.h>





// Z is constant at 8.429
// R goes from 4.9 to 5.9
void horizontal_pfunc_upper(double s, double* RZ){
  RZ[0] = 4.9 + s;
  RZ[1] = 8.429;
}

// Z is constant at -8.429
// R goes from 4.9 to 5.3
void horizontal_pfunc_lower(double s, double* RZ){
  RZ[0] = 4.9 + s;
  RZ[1] = -8.429;
}

// R is constant at 4.9
// R goes from 8.0 to 8.5
void vertical_pfunc_upper(double s, double* RZ){
  RZ[0] = 4.9;
  RZ[1] = 8.0 + s/2;
}

void vertical_pfunc_lower(double s, double* RZ){
  RZ[0] = 4.9;
  RZ[1] = -8.0 - s/2;
}

// Actual Plate info:
// p1 = [5.7413,8.5471]
// p2 = [5.8544, 8.5229]
// p3 = [5.8549, 8.4258]
// Try a Different (slanted) plate instead
// p1 [5.151,8.516]
// p2 [5.852, 8.434]
void shaped_pfunc_upper(double s, double* RZ){
    RZ[0] = 0.8 + (0.916 - 0.8)*s;
    RZ[1] = -1.2 + (-1.329 + 1.2)*s;
}

void shaped_pfunc_lower(double s, double* RZ){
    RZ[0] = 1.6 + (1.8 - 1.6)*s;
    RZ[1] = -1.26 + (-1.1 + 1.26)*s;
}

void
write_geometry(gk_geometry *up, struct gkyl_rect_grid grid, struct gkyl_basis basis, struct gkyl_range local, const char *name)
{
  const char *fmt = "%s-%s.gkyl";
  int sz = gkyl_calc_strlen(fmt, name, "jacobtot_inv");
  char fileNm[sz+1]; // ensure no buffer overflow

  sprintf(fileNm, fmt, name, "mapc2p");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->geo_corn.mc2p, fileNm);
  sprintf(fileNm, fmt, name, "bmag");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->geo_corn.bmag, fileNm);
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

  // Write Nodal Coordinates
  struct gkyl_range nrange;
  gkyl_gk_geometry_init_nodal_range(&nrange, &local, 1);
  struct gkyl_array* mc2p_nodal = gkyl_array_new(GKYL_DOUBLE, 3, nrange.volume);
  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&basis, &grid, false);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange, &local, 3, mc2p_nodal, up->geo_corn.mc2p, false);
  gkyl_nodal_ops_release(n2m);
  struct gkyl_rect_grid ngrid;
  gkyl_gk_geometry_init_nodal_grid(&ngrid, &grid, &nrange);
  sprintf(fileNm, fmt, name, "nodes");
  gkyl_grid_sub_array_write(&ngrid, &nrange, 0,  mc2p_nodal, fileNm);
  gkyl_array_release(mc2p_nodal);
}




void
test_fixed_z()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();



  struct gkyl_efit_inp efit_inp = {
      // psiRZ and related inputs
      .filepath = "gyrokinetic/data/eqdsk/asdex.geqdsk",
      .rz_poly_order = 2,
      .flux_poly_order = 1,
    };

  double clower[] = { 0.16, -0.01, -M_PI+1e-14 };
  double cupper[] = {0.17501, 0.01, M_PI-1e-14 };

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
    .ftype = GKYL_GEOMETRY_TOKAMAK_LSN_SOL,
    .rclose = 2.5,
    .rright = 2.5,
    .rleft = 0.7,
    .zmin = -1.3,
    .zmax = 1.0,
    .zmin_left = -1.3,
    .zmin_right = -1.3,
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

void
test_shaped_plate()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();



  struct gkyl_efit_inp efit_inp = {
      // psiRZ and related inputs
      .filepath = "gyrokinetic/data/eqdsk/asdex.geqdsk",
      .rz_poly_order = 2,
      .flux_poly_order = 1,
    };

  double clower[] = { 0.16, -0.01, -M_PI+1e-14 };
  double cupper[] = {0.17501, 0.01, M_PI-1e-14 };

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
  //write_geometry(up, cgrid, cbasis, clocal, "asdex");
  gkyl_gk_geometry_release(up);
  gkyl_position_map_release(pmap);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
}

void
test_lower()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();



  struct gkyl_efit_inp efit_inp = {
      // psiRZ and related inputs
      .filepath = "gyrokinetic/data/eqdsk/asdex.geqdsk",
      .rz_poly_order = 2,
      .flux_poly_order = 1,
    };

  double clower[] = { 0.16, -0.01, -M_PI+1e-14 };
  double cupper[] = {0.17501, 0.01, M_PI-1e-14 };
  int ccells[] = { 2, 1, 2 };

  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO,
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
  gkyl_gk_geometry_tok_set_grid_extents(efit_inp, ginp, &clower[2], &cupper[2]);

  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  int cpoly_order = 1;
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

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

  struct gk_geometry* up = gkyl_gk_geometry_tok_new(&geometry_inp); 
  //write_geometry(up, cgrid, cbasis, clocal, "asdexlo");
  gkyl_gk_geometry_release(up);
  gkyl_position_map_release(pmap);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
}

void
test_middle()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();



  struct gkyl_efit_inp efit_inp = {
      // psiRZ and related inputs
      .filepath = "gyrokinetic/data/eqdsk/asdex.geqdsk",
      .rz_poly_order = 2,
      .flux_poly_order = 1,
    };

  double clower[] = { 0.16, -0.01, -M_PI+1e-14 };
  double cupper[] = {0.17501, 0.01, M_PI-1e-14 };
  int ccells[] = { 2, 1, 2 };

  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID,
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
  gkyl_gk_geometry_tok_set_grid_extents(efit_inp, ginp, &clower[2], &cupper[2]);

  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  int cpoly_order = 1;
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

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

  struct gk_geometry* up = gkyl_gk_geometry_tok_new(&geometry_inp); 
  //write_geometry(up, cgrid, cbasis, clocal, "asdexmid");
  gkyl_gk_geometry_release(up);
  gkyl_position_map_release(pmap);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
}

void
test_upper()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();



  struct gkyl_efit_inp efit_inp = {
      // psiRZ and related inputs
      .filepath = "gyrokinetic/data/eqdsk/asdex.geqdsk",
      .rz_poly_order = 2,
      .flux_poly_order = 1,
    };

  double clower[] = { 0.16, -0.01, -M_PI+1e-14 };
  double cupper[] = {0.17501, 0.01, M_PI-1e-14 };
  int ccells[] = { 2, 1, 2 };

  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_GEOMETRY_TOKAMAK_LSN_SOL_UP,
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
  gkyl_gk_geometry_tok_set_grid_extents(efit_inp, ginp, &clower[2], &cupper[2]);

  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  int cpoly_order = 1;
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

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

  struct gk_geometry* up = gkyl_gk_geometry_tok_new(&geometry_inp); 
  //write_geometry(up, cgrid, cbasis, clocal, "asdexup");
  gkyl_gk_geometry_release(up);
  gkyl_position_map_release(pmap);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
}



TEST_LIST = {
  //{ "test_fixed_z", test_fixed_z},
  { "test_shaped_plate", test_shaped_plate},
  { "test_lower", test_lower},
  { "test_middle", test_middle},
  { "test_upper", test_upper},
  { NULL, NULL },
};
