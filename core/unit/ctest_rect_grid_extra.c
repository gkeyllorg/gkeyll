// Additional unit tests for gkyl_rect_grid (cell geometry helpers)
#include <acutest.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <math.h>

void
test_grid_init_1d()
{
  double lower[] = { 0.0 }, upper[] = { 1.0 };
  int cells[] = { 10 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);

  TEST_CHECK( grid.ndim == 1 );
  TEST_CHECK( grid.cells[0] == 10 );
  TEST_CHECK( gkyl_compare_double(grid.lower[0], 0.0, 1e-15) );
  TEST_CHECK( gkyl_compare_double(grid.upper[0], 1.0, 1e-15) );
  TEST_CHECK( gkyl_compare_double(grid.dx[0], 0.1, 1e-15) );
  TEST_CHECK( gkyl_compare_double(grid.cellVolume, 0.1, 1e-15) );
}

void
test_grid_init_3d()
{
  double lower[] = { -1.0, 0.0, 2.0 }, upper[] = { 1.0, 4.0, 6.0 };
  int cells[] = { 4, 8, 2 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 3, lower, upper, cells);

  TEST_CHECK( grid.ndim == 3 );
  TEST_CHECK( gkyl_compare_double(grid.dx[0], 0.5, 1e-15) );
  TEST_CHECK( gkyl_compare_double(grid.dx[1], 0.5, 1e-15) );
  TEST_CHECK( gkyl_compare_double(grid.dx[2], 2.0, 1e-15) );
  // cellVolume = product of dx
  TEST_CHECK( gkyl_compare_double(grid.cellVolume, 0.5*0.5*2.0, 1e-15) );
}

void
test_grid_cell_center()
{
  double lower[] = { 0.0, 0.0 }, upper[] = { 2.0, 2.0 };
  int cells[] = { 2, 2 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  // dx = 1.0 in both directions; cell 1 center should be at 0.5
  int idx[] = { 1, 1 };
  double xc[2];
  gkyl_rect_grid_cell_center(&grid, idx, xc);
  TEST_CHECK( gkyl_compare_double(xc[0], 0.5, 1e-15) );
  TEST_CHECK( gkyl_compare_double(xc[1], 0.5, 1e-15) );

  int idx2[] = { 2, 2 };
  gkyl_rect_grid_cell_center(&grid, idx2, xc);
  TEST_CHECK( gkyl_compare_double(xc[0], 1.5, 1e-15) );
  TEST_CHECK( gkyl_compare_double(xc[1], 1.5, 1e-15) );
}

void
test_grid_ll_node()
{
  double lower[] = { 0.0 }, upper[] = { 1.0 };
  int cells[] = { 4 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);

  // dx = 0.25; lower-left node of cell i is lower + (i-1)*dx
  int idx[] = { 1 };
  double xc[1];
  gkyl_rect_grid_ll_node(&grid, idx, xc);
  TEST_CHECK( gkyl_compare_double(xc[0], 0.0, 1e-15) );

  idx[0] = 3;
  gkyl_rect_grid_ll_node(&grid, idx, xc);
  TEST_CHECK( gkyl_compare_double(xc[0], 0.5, 1e-15) );

  // cell center is ll_node + dx/2
  double cc[1];
  gkyl_rect_grid_cell_center(&grid, idx, cc);
  TEST_CHECK( gkyl_compare_double(cc[0], xc[0] + 0.5*grid.dx[0], 1e-15) );
}

void
test_grid_extents()
{
  double lower[] = { 0.0, 0.0 }, upper[] = { 1.0, 1.0 };
  int cells[] = { 5, 7 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  int ext[2];
  gkyl_rect_grid_extents(&grid, 0, ext);
  TEST_CHECK( ext[0] == 1 );
  TEST_CHECK( ext[1] == 5 );

  gkyl_rect_grid_extents(&grid, 1, ext);
  TEST_CHECK( ext[0] == 1 );
  TEST_CHECK( ext[1] == 7 );
}

void
test_grid_coord_idx()
{
  double lower[] = { 0.0 }, upper[] = { 10.0 };
  int cells[] = { 10 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);

  // dx = 1.0. x=0.5 -> cell 1, x=3.5 -> cell 4, x=9.9 -> cell 10
  double xn[1];
  int idx[1];

  xn[0] = 0.5; gkyl_rect_grid_coord_idx(&grid, xn, idx);
  TEST_CHECK( idx[0] == 1 );

  xn[0] = 3.5; gkyl_rect_grid_coord_idx(&grid, xn, idx);
  TEST_CHECK( idx[0] == 4 );

  xn[0] = 9.9; gkyl_rect_grid_coord_idx(&grid, xn, idx);
  TEST_CHECK( idx[0] == 10 );
}

void
test_grid_cmp()
{
  double lower[] = { 0.0, 0.0 }, upper[] = { 1.0, 1.0 };
  int cells[] = { 4, 4 };
  struct gkyl_rect_grid g1, g2, g3;
  gkyl_rect_grid_init(&g1, 2, lower, upper, cells);
  gkyl_rect_grid_init(&g2, 2, lower, upper, cells);
  TEST_CHECK( gkyl_rect_grid_cmp(&g1, &g2) );

  int cells2[] = { 4, 8 };
  gkyl_rect_grid_init(&g3, 2, lower, upper, cells2);
  TEST_CHECK( !gkyl_rect_grid_cmp(&g1, &g3) );
}

TEST_LIST = {
  { "grid_init_1d", test_grid_init_1d },
  { "grid_init_3d", test_grid_init_3d },
  { "grid_cell_center", test_grid_cell_center },
  { "grid_ll_node", test_grid_ll_node },
  { "grid_extents", test_grid_extents },
  { "grid_coord_idx", test_grid_coord_idx },
  { "grid_cmp", test_grid_cmp },
  { NULL, NULL },
};
