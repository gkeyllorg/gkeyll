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
#include <gkyl_tok_geo.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_tok.h>
#include <gkyl_nodal_ops.h>

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

struct gkyl_efit_inp inp = {
  // psiRZ and related inputs
  .filepath = "gyrokinetic/data/eqdsk/ltx_miller.geqdsk",
  .rz_poly_order = 2,
  .flux_poly_order = 1,
  .reflect = true,
};

int cpoly_order = 1;
struct gkyl_basis cbasis;
int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
struct gkyl_rect_grid cgrid;
struct gkyl_range clocal, clocal_ext;

void
test_ltx_miller()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  double clower[] = { 0.0018, -0.01, -M_PI+1e-14 };
  double cupper[] = {0.0024, 0.01, M_PI-1e-14 };
  int ccells[] = { 2,1,2 };

  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_GEOMETRY_TOKAMAK_IWL,
    .rclose = 0.4,
    .rleft= 0.2,
    .rright= 0.45,
    .rmin=0.1,
    .rmax=0.65,
    .zmin = -0.3,
    .zmax = 0.3,
  }; 

  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_position_map *pmap = gkyl_position_map_null_new();

  struct gkyl_gk_geometry_inp geometry_inp = {
    .geometry_id  = GKYL_GEOMETRY_TOKAMAK,
    .efit_info = inp,
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
  write_geometry(up, cgrid, cbasis, clocal, "ltx_miller");
  gkyl_gk_geometry_release(up);
  gkyl_position_map_release(pmap);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

  gkyl_position_map_release(pmap);
}



TEST_LIST = {
  { "test_ltx_miller", test_ltx_miller},
  { NULL, NULL },
};
