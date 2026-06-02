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
