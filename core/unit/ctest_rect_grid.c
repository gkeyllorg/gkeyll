#include <acutest.h>

#include <gkyl_alloc.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

void test_grid_2d()
{
  double lower[] = {1.0, 1.0}, upper[] = {2.5, 5.0};
  int cells[] = {20, 20};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  TEST_CHECK( grid.ndim == 2 );
  for (int i=0; i<grid.ndim; ++i) {
    TEST_CHECK( grid.lower[i] == lower[i] );
    TEST_CHECK( grid.upper[i] == upper[i] );
    TEST_CHECK( grid.cells[i] == cells[i] );
    TEST_CHECK( grid.dx[i] == (upper[i]-lower[i])/cells[i] );
  }
  TEST_CHECK( grid.cellVolume == 0.075*0.2 );

  double xc[2];  
  int x_ext[2], y_ext[2];
  
  gkyl_rect_grid_extents(&grid, 0, x_ext);
  gkyl_rect_grid_extents(&grid, 1, y_ext);  

  for (int i=x_ext[0]; i<=x_ext[1]; ++i)
    for (int j=y_ext[0]; j<=y_ext[1]; ++j) {
      gkyl_rect_grid_cell_center(&grid, (int[2]) { i, j }, xc);

      TEST_CHECK( xc[0] == 1.0 + (i-0.5)*grid.dx[0] );
      TEST_CHECK( xc[1] == 1.0 + (j-0.5)*grid.dx[1] );
    }

  double xn[2] = { 0.0, 0.0 };
  int idx[2];

  xn[0] = grid.lower[0] + 0.5*grid.dx[0];
  xn[1] = grid.lower[1] + 0.5*grid.dx[1];
  gkyl_rect_grid_coord_idx(&grid, xn, idx);
  TEST_CHECK( (idx[0] == 1) && (idx[1] == 1) );

  xn[0] = grid.lower[0] + 1.5*grid.dx[0];
  xn[1] = grid.lower[1] + 1.5*grid.dx[1];
  gkyl_rect_grid_coord_idx(&grid, xn, idx);
  TEST_CHECK( (idx[0] == 2) && (idx[1] == 2) );

  TEST_CHECK( gkyl_rect_grid_cmp(&grid, &grid) == true );

  double lower2[] = {1.0, 0.5}, upper2[] = {2.5, 5.0};
  int cells2[] = {20, 19};
  struct gkyl_rect_grid grid2;
  gkyl_rect_grid_init(&grid2, 2, lower2, upper2, cells2);

  TEST_CHECK( gkyl_rect_grid_cmp(&grid, &grid2) == false );
}

/* Test rect_grid find cell in 1D
 * First test: Generic
 * Second test: Point nearly on cell boundary
 * Third test: Point on cell boundary, choose lower cell
 * Fourth test: One known index given (simply tests that it is the correct index)
 */
void test_find_cell_1d(){
  double lower[] = {0.0}, upper[] = {5.0};
  int cells[] = {5};
  double point[] = {2.5};
  bool pick_lower[1] = {false};
  const int known_index[1] = {-1};
  int cell_index[] = {0};
  struct gkyl_rect_grid grid;

  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);
  gkyl_rect_grid_find_cell(&grid, point, pick_lower, known_index, cell_index);
  int correct_index[1] = {3};
  for (int i=0; i<1; i++)
    TEST_CHECK( cell_index[i] == correct_index[i] );

  point[0]=2.0+1e-15;
  gkyl_rect_grid_find_cell(&grid, point, pick_lower, known_index, cell_index);
  for (int i=0; i<1; i++)
    TEST_CHECK( cell_index[i] == correct_index[i] );

  pick_lower[0] = true;
  correct_index[0] = 2;
  gkyl_rect_grid_find_cell(&grid, point, pick_lower, known_index, cell_index);
  for (int i=0; i<1; i++)
    TEST_CHECK( cell_index[i] == correct_index[i]);

  int idx = 1;
  correct_index[0] = idx;
  point[0] = 0.2;
  const int known_index2[1] = {idx};
  gkyl_rect_grid_find_cell(&grid, point, pick_lower, known_index2, cell_index);
  for (int i=0; i<1; i++)
    TEST_CHECK( cell_index[i] == correct_index[i]);
}

/* Test rect_grid find cell in 2D
 * First test: Generic 2D
 * Second test: Point on cell corner
 * Third test: Point on cell corner (same location as second), but pick lower index
 * Fourth test: One index is known 
 */
void test_find_cell_2d(){
  double lower[] = {0.0, -10.0}, upper[] = {5.0, 10.0};
  int cells[] = {5, 20};
  double point[] = {2.5, 1.3};
  bool pick_lower[2] = {false, false};
  const int known_index[2] = {-1, -1};
  int cell_index[] = {0, 0};
  struct gkyl_rect_grid grid;

  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);
  gkyl_rect_grid_find_cell(&grid, point, pick_lower, known_index, cell_index);
  int correct_index[2] = {3, 12};
  for (int i=0; i<2; i++)
    TEST_CHECK( cell_index[i] == correct_index[i]);

  point[0] = 2.0;
  point[1] = 1.0;
  gkyl_rect_grid_find_cell(&grid, point, pick_lower, known_index, cell_index);
  for (int i=0; i<2; i++)
    TEST_CHECK( cell_index[i] == correct_index[i]);

  pick_lower[0] = pick_lower[1] = true;
  correct_index[0] = 2;
  correct_index[1] = 11;
  gkyl_rect_grid_find_cell(&grid, point, pick_lower, known_index, cell_index);
  for (int i=0; i<2; i++)
    TEST_CHECK( cell_index[i] == correct_index[i]);

  pick_lower[0] = pick_lower[1] = false;
  int idx = 18;
  correct_index[0] = 3;
  correct_index[1] = idx;
  point[1] = 7.5;
  const int known_index2[2] = {-1, idx};
  gkyl_rect_grid_find_cell(&grid, point, pick_lower, known_index2, cell_index);
  for (int i=0; i<2; i++)
    TEST_CHECK( cell_index[i] == correct_index[i]);
}

/* Test rect_grid find cell in 3D
 * First test: Generic 3D
 * Second test: Point on just to one side of cell boundary
 * Third test: Point on just to one side of cell boundary, but pick lower index
 * Fourth test: 2 indecies are known
 */
void test_find_cell_3d(){
  double lower[] = {0.0, -10.0, 1.3}, upper[] = {5.0, 10.0, 2.5};
  int cells[] = {5, 20, 100};
  double point[] = {2.5, 1.3, 1.4};
  bool pick_lower[3] = {false, false, false};
  const int known_index[] = {-1, -1, -1};
  int cell_index[] = {0, 0, 0};
  struct gkyl_rect_grid grid;

  gkyl_rect_grid_init(&grid, 3, lower, upper, cells);
  gkyl_rect_grid_find_cell(&grid, point, pick_lower, known_index, cell_index);
  int correct_index[] = {3, 12, 9};
  for (int i=0; i<3; i++)
    TEST_CHECK( cell_index[i] == correct_index[i]);

  point[0] = 2.0;
  point[1] = 1.0+1e-15;
  gkyl_rect_grid_find_cell(&grid, point, pick_lower, known_index, cell_index);
  for (int i=0; i<3; i++)
    TEST_CHECK( cell_index[i] == correct_index[i]);

  pick_lower[0] = pick_lower[1] = pick_lower[2] = true;
  correct_index[0] = 2;
  correct_index[1] = 11;
  gkyl_rect_grid_find_cell(&grid, point, pick_lower, known_index, cell_index);
  for (int i=0; i<3; i++)
    TEST_CHECK( cell_index[i] == correct_index[i]);

  pick_lower[0] = pick_lower[1] = pick_lower[2] = false;
  correct_index[0] = 3;
  int idx = 18, idx2 = 35;
  correct_index[1] = idx;
  correct_index[2] = idx2;
  point[1] = 7.5;
  point[2] = 1.71;
  const int known_index2[] = {-1, idx, idx2};
  gkyl_rect_grid_find_cell(&grid, point, pick_lower, known_index2, cell_index);
  for (int i=0; i<3; i++)
    TEST_CHECK( cell_index[i] == correct_index[i]);
}

void test_grid_io()
{
  double lower[] = {1.0, 1.0}, upper[] = {2.5, 5.0};
  int cells[] = {20, 20};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  FILE *fp = 0;
  with_file (fp, "ctest_rect_grid.dat", "w")
    gkyl_rect_grid_write(&grid, 0, fp);

  struct gkyl_rect_grid grid2;
  with_file (fp, "ctest_rect_grid.dat", "r")
    gkyl_rect_grid_read(&grid2, fp);

  TEST_CHECK( grid.ndim == grid2.ndim );
  for (int d=0; d<grid.ndim; ++d) {
    TEST_CHECK( grid.lower[d] == grid2.lower[d] );
    TEST_CHECK( grid.upper[d] == grid2.upper[d] );
    TEST_CHECK( grid.cells[d] == grid2.cells[d] );
    TEST_CHECK( grid.dx[d] == grid2.dx[d] );
  }
  TEST_CHECK( grid.cellVolume == grid2.cellVolume );
}

// CUDA specific tests
#ifdef GKYL_HAVE_CUDA

int cu_rect_grid_test(const struct gkyl_rect_grid grid);

void test_cu_grid_2d()
{
  double lower[] = {1.0, 1.0}, upper[] = {2.5, 5.0};
  int cells[] = {20, 20};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  int nfail = cu_rect_grid_test(grid);
  TEST_CHECK( nfail == 0 );
}

#endif

TEST_LIST = {
  { "grid_2d", test_grid_2d },
  { "grid_find_cell_1d", test_find_cell_1d },
  { "grid_find_cell_2d", test_find_cell_2d },
  { "grid_find_cell_3d", test_find_cell_3d },
  { "grid_io", test_grid_io },
#ifdef GKYL_HAVE_CUDA
  { "cu_grid_2d", test_cu_grid_2d },
#endif  
  { NULL, NULL },
};
