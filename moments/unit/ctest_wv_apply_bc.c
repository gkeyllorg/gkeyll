#include <acutest.h>
#include <math.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wv_apply_bc.h>
#include <gkyl_wv_burgers.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_wv_euler.h>

static void
nomapc2p(double t, const double *xc, double *xp, void *ctx)
{
  int *ndim = ctx;
  for (int i=0; i<(*ndim); ++i) xp[i] = xc[i];
}

static void
rtheta_map(double t, const double *xc, double *xp, void *ctx)
{
  double r = xc[0], th  = xc[1];
  xp[0] = r*cos(th); xp[1] = r*sin(th);
}

static void
bc_copy(const struct gkyl_wv_eqn* eqn, double t, int nc, const double *skin, double *restrict ghost, void *ctx)
{
  for (int c=0; c<nc; ++c) ghost[c] = skin[c];
}

void
test_1()
{
  int ndim = 1;
  double lower[] = {-1.0}, upper[] = {1.0};
  int cells[] = {16};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  int nghost[GKYL_MAX_DIM] = { 2 };
  struct gkyl_range range, ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);
  
  struct gkyl_wave_geom *wg = gkyl_wave_geom_new(&grid, &ext_range, nomapc2p, &ndim, false);
  struct gkyl_wv_eqn *eqn = gkyl_wv_burgers_new(false);

  gkyl_wv_apply_bc *lbc = gkyl_wv_apply_bc_new(&grid, eqn, wg,
    0, GKYL_LOWER_EDGE, nghost, bc_copy, NULL);
  gkyl_wv_apply_bc *rbc = gkyl_wv_apply_bc_new(&grid, eqn, wg,
    0, GKYL_UPPER_EDGE, nghost, bc_copy, NULL);

  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, 1, ext_range.volume);

  gkyl_array_clear_range(distf, 1.0, &range);

  // check if ghost-cells on left/right edges of domain are 0.0
  double *data = distf->data;
  TEST_CHECK( 0.0 == data[0] );
  TEST_CHECK( 0.0 == data[1] );

  TEST_CHECK( 0.0 == data[18] );
  TEST_CHECK( 0.0 == data[19] );
  
  // apply BC
  gkyl_wv_apply_bc_advance(lbc, 0.0, &range, distf);
  gkyl_wv_apply_bc_advance(rbc, 0.0, &range, distf);

  // check if BCs applied correctly
  TEST_CHECK( 1.0 == data[0] );
  TEST_CHECK( 1.0 == data[1] );

  TEST_CHECK( 1.0 == data[18] );
  TEST_CHECK( 1.0 == data[19] );

  gkyl_wv_apply_bc_release(lbc);
  gkyl_wv_apply_bc_release(rbc);
  gkyl_wv_eqn_release(eqn);
  gkyl_wave_geom_release(wg);
  gkyl_array_release(distf);
}

void
test_2()
{
  int ndim = 2;
  double lower[] = {-1.0, -1.0}, upper[] = {1.0, 1.0};
  int cells[] = {16, 8};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  int nghost[] = { 2, 2 };

  struct gkyl_range ext_range, range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  struct gkyl_wave_geom *wg = gkyl_wave_geom_new(&grid, &ext_range, nomapc2p, &ndim, false);
  struct gkyl_wv_eqn *eqn = gkyl_wv_burgers_new(false);  

  gkyl_wv_apply_bc *lbc = gkyl_wv_apply_bc_new(&grid, eqn, wg,
    0, GKYL_LOWER_EDGE, nghost, bc_copy, NULL);
  gkyl_wv_apply_bc *rbc = gkyl_wv_apply_bc_new(&grid, eqn, wg,
    0, GKYL_UPPER_EDGE, nghost, bc_copy, NULL);

  gkyl_wv_apply_bc *bbc = gkyl_wv_apply_bc_new(&grid, eqn, wg,
    1, GKYL_LOWER_EDGE, nghost, bc_copy, NULL);
  gkyl_wv_apply_bc *tbc = gkyl_wv_apply_bc_new(&grid, eqn, wg,
    1, GKYL_UPPER_EDGE, nghost, bc_copy, NULL);

  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, 1, ext_range.volume);

  // clear interior of array
  gkyl_array_clear_range(distf, 1.0, &range);

  // check if only interior is cleared
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &ext_range);

  double vol = 0.0;
  while (gkyl_range_iter_next(&iter)) {
    double *f = gkyl_array_fetch(distf, gkyl_range_idx(&ext_range, iter.idx));
    vol += f[0];
  }
  TEST_CHECK( vol == range.volume );

  // apply various BCs

  gkyl_wv_apply_bc_advance(lbc, 0.0, &range, distf);
  gkyl_wv_apply_bc_advance(rbc, 0.0, &range, distf);

  gkyl_wv_apply_bc_advance(bbc, 0.0, &range, distf);
  gkyl_wv_apply_bc_advance(tbc, 0.0, &range, distf);

  gkyl_range_iter_init(&iter, &ext_range);

  vol = 0.0;
  while (gkyl_range_iter_next(&iter)) {
    double *f = gkyl_array_fetch(distf, gkyl_range_idx(&ext_range, iter.idx));
    vol += f[0];
  }

  // volume should be volume of ext_range but as corners are not
  // touched by BC updater we need to subtract the volume of the 4 corners
  TEST_CHECK( vol == ext_range.volume-4*4 );

  gkyl_wv_apply_bc_release(lbc);
  gkyl_wv_apply_bc_release(rbc);
  gkyl_wv_apply_bc_release(bbc);
  gkyl_wv_apply_bc_release(tbc);
  gkyl_wv_eqn_release(eqn);
  gkyl_wave_geom_release(wg);
  gkyl_array_release(distf);
}

void
test_3()
{
  int ndim = 2;
  double lower[] = {-1.0, -1.0}, upper[] = {1.0, 1.0};
  int cells[] = {16, 8};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  int nghost[] = { 2, 1 };

  struct gkyl_range ext_range, range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  struct gkyl_wave_geom *wg = gkyl_wave_geom_new(&grid, &ext_range, nomapc2p, &ndim, false);
  struct gkyl_wv_eqn *eqn = gkyl_wv_burgers_new(false);

  gkyl_wv_apply_bc *lbc = gkyl_wv_apply_bc_new(&grid, eqn, wg,
    0, GKYL_LOWER_EDGE, nghost, bc_copy, NULL);
  gkyl_wv_apply_bc *rbc = gkyl_wv_apply_bc_new(&grid, eqn, wg,
    0, GKYL_UPPER_EDGE, nghost, bc_copy, NULL);

  gkyl_wv_apply_bc *bbc = gkyl_wv_apply_bc_new(&grid, eqn, wg,
    1, GKYL_LOWER_EDGE, nghost, bc_copy, NULL);
  gkyl_wv_apply_bc *tbc = gkyl_wv_apply_bc_new(&grid, eqn, wg,
    1, GKYL_UPPER_EDGE, nghost, bc_copy, NULL);

  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, 1, ext_range.volume);

  // clear interior of array
  gkyl_array_clear(distf, 0.0);
  gkyl_array_clear_range(distf, 1.0, &range);

  // check if only interior is cleared
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &ext_range);

  double vol = 0.0;
  while (gkyl_range_iter_next(&iter)) {
    double *f = gkyl_array_fetch(distf, gkyl_range_idx(&ext_range, iter.idx));
    vol += f[0];
  }
  TEST_CHECK( vol == range.volume );

  /** apply BCs on restricted range: test 1 */
  struct gkyl_range sub_range;
  gkyl_sub_range_init(&sub_range, &ext_range, (int []) { 1, 2 }, (int []) { 10, 6 });

  gkyl_wv_apply_bc_advance(lbc, 0.0, &sub_range, distf);
  gkyl_wv_apply_bc_advance(rbc, 0.0, &sub_range, distf);

  gkyl_wv_apply_bc_advance(bbc, 0.0, &sub_range, distf);
  gkyl_wv_apply_bc_advance(tbc, 0.0, &sub_range, distf);

  gkyl_range_iter_init(&iter, &ext_range);

  vol = 0.0;
  while (gkyl_range_iter_next(&iter)) {
    double *f = gkyl_array_fetch(distf, gkyl_range_idx(&ext_range, iter.idx));
    vol += f[0];
  }

  // volume should be volume of range + nghost*5 (as only small
  // portion of boundary is updated)
  TEST_CHECK( vol == range.volume+2*5 );

  /** apply BCs on restricted range: test 2 */

  gkyl_array_clear(distf, 0.0);
  gkyl_array_clear_range(distf, 1.0, &range);
  
  gkyl_sub_range_init(&sub_range, &ext_range, (int []) { 2, 1 }, (int []) { 8, 4 });

  gkyl_wv_apply_bc_advance(lbc, 0.0, &sub_range, distf);
  gkyl_wv_apply_bc_advance(rbc, 0.0, &sub_range, distf);

  gkyl_wv_apply_bc_advance(bbc, 0.0, &sub_range, distf);
  gkyl_wv_apply_bc_advance(tbc, 0.0, &sub_range, distf);

  gkyl_range_iter_init(&iter, &ext_range);

  vol = 0.0;
  while (gkyl_range_iter_next(&iter)) {
    double *f = gkyl_array_fetch(distf, gkyl_range_idx(&ext_range, iter.idx));
    vol += f[0];
  }
  
  // volume should be volume of range + nghost*7 (as only small
  // portion of boundary is updated)
  TEST_CHECK( vol == range.volume+1*7 );

  /** apply BCs on restricted range: test 3 */

  gkyl_array_clear(distf, 0.0);
  gkyl_array_clear_range(distf, 1.0, &range);
  
  gkyl_sub_range_init(&sub_range, &ext_range, (int []) { 2, 2 }, (int []) { 8, 4 });

  gkyl_wv_apply_bc_advance(lbc, 0.0, &sub_range, distf);
  gkyl_wv_apply_bc_advance(rbc, 0.0, &sub_range, distf);

  gkyl_wv_apply_bc_advance(bbc, 0.0, &sub_range, distf);
  gkyl_wv_apply_bc_advance(tbc, 0.0, &sub_range, distf);

  gkyl_range_iter_init(&iter, &ext_range);

  vol = 0.0;
  while (gkyl_range_iter_next(&iter)) {
    double *f = gkyl_array_fetch(distf, gkyl_range_idx(&ext_range, iter.idx));
    vol += f[0];
  }

  // range does not touch boundaries
  TEST_CHECK( vol == range.volume );

  gkyl_wv_apply_bc_release(lbc);
  gkyl_wv_apply_bc_release(rbc);
  gkyl_wv_apply_bc_release(bbc);
  gkyl_wv_apply_bc_release(tbc);
  gkyl_wv_eqn_release(eqn);
  gkyl_wave_geom_release(wg);
  gkyl_array_release(distf);
}

// stuff to help in wedge BCs
struct skin_ghost_ranges {
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];

  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
};

// Create ghost and skin sub-ranges given a parent range
static void
skin_ghost_ranges_init(struct skin_ghost_ranges *sgr,
  const struct gkyl_range *parent, const int *ghost)
{
  int ndim = parent->ndim;
  
  for (int d=0; d<ndim; ++d) {
    gkyl_skin_ghost_ranges(&sgr->lower_skin[d], &sgr->lower_ghost[d],
      d, GKYL_LOWER_EDGE, parent, ghost);
    gkyl_skin_ghost_ranges(&sgr->upper_skin[d], &sgr->upper_ghost[d],
      d, GKYL_UPPER_EDGE, parent, ghost);
  }
}

void
test_bc_buff_rtheta()
{
  int ndim = 2;
  double lower[] = {0.25, 0.0}, upper[] = {1.25, 2*M_PI/4};
  int cells[] = {16, 8};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  struct gkyl_wv_eqn *eqn = gkyl_wv_euler_new(1.4, false);  

  int nghost[] = { 2, 2 };

  struct gkyl_range ext_range, range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  struct skin_ghost_ranges skin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &ext_range, nghost);

  long buff_sz = 0;
  // compute buffer size needed
  for (int d=0; d<grid.ndim; ++d) {
    long vol = skin_ghost.lower_skin[d].volume;
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  struct gkyl_array *bc_buffer = gkyl_array_new(GKYL_DOUBLE,
    eqn->num_equations, buff_sz);

  struct gkyl_wave_geom *wg = gkyl_wave_geom_new(&grid, &ext_range,
    rtheta_map, &ndim, false);

  gkyl_wv_apply_bc *bbc = gkyl_wv_apply_bc_new(&grid, eqn, wg,
    1, GKYL_LOWER_EDGE, nghost, bc_copy, NULL);
  gkyl_wv_apply_bc *tbc = gkyl_wv_apply_bc_new(&grid, eqn, wg,
    1, GKYL_UPPER_EDGE, nghost, bc_copy, NULL);

  struct gkyl_array *fluid = gkyl_array_new(GKYL_DOUBLE,
    eqn->num_equations, ext_range.volume);

  // set interior of array
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);

  while (gkyl_range_iter_next(&iter)) {
    long sloc = gkyl_range_idx(&range, iter.idx);
    double *q = gkyl_array_fetch(fluid, sloc);

    double xc[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(&grid, iter.idx, xc);

    double r = xc[0], th = xc[1];

    q[0] = r*th*0.1;
    q[1] = r*th*0.5;
    q[2] = r*th*1.5;
    q[3] = r*th*2.5;
    q[4] = r*th*10.5;
  }

  // apply various BCs, copying output to a buffer

  // bottom
  gkyl_array_clear(bc_buffer, 0.0);
  gkyl_wv_apply_bc_to_buff(bbc, 0.0, &range, fluid, bc_buffer->data);

  // check if data copied properly into buffer
  long count = 0;    
  gkyl_range_iter_init(&iter, &skin_ghost.lower_skin[1]);
  while (gkyl_range_iter_next(&iter)) {
    double xc[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(&grid, iter.idx, xc);

    double r = xc[0], th = xc[1];

    const double *q = gkyl_array_cfetch(bc_buffer, count++);

    TEST_CHECK( q[0] == r*th*0.1 );

    // rotate counter-clockwise by pi/4
    TEST_CHECK( gkyl_compare(q[1], -r*th*1.5, 1e-14) );
    TEST_CHECK( gkyl_compare(q[2], r*th*0.5, 1e-14) );
    
    TEST_CHECK( gkyl_compare(q[3], r*th*2.5, 1e-14) );
    TEST_CHECK( q[4] == r*th*10.5 );
  }

  // right
  gkyl_array_clear(bc_buffer, 0.0);
  gkyl_wv_apply_bc_to_buff(tbc, 0.0, &range, fluid, bc_buffer->data);

  // check if data copied properly into buffer
  count = 0;  
  gkyl_range_iter_init(&iter, &skin_ghost.upper_skin[1]);
  while (gkyl_range_iter_next(&iter)) {
    double xc[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(&grid, iter.idx, xc);

    double r = xc[0], th = xc[1];

    const double *q = gkyl_array_cfetch(bc_buffer, count++);

    TEST_CHECK( q[0] == r*th*0.1 );

    // rotate clockwise by pi/4
    TEST_CHECK( gkyl_compare(q[1], r*th*1.5, 1e-14) );
    TEST_CHECK( gkyl_compare(q[2], -r*th*0.5, 1e-14) );
    
    TEST_CHECK( gkyl_compare(q[3], r*th*2.5, 1e-14) );
    TEST_CHECK( q[4] == r*th*10.5 );
  }  

  gkyl_wv_apply_bc_release(bbc);
  gkyl_wv_apply_bc_release(tbc);
  gkyl_wv_eqn_release(eqn);
  gkyl_wave_geom_release(wg);
  gkyl_array_release(fluid);
  gkyl_array_release(bc_buffer);
}

TEST_LIST = {
  { "test_1", test_1 },
  { "test_2", test_2 },
  { "test_3", test_3 },
  { "test_bc_buff_rtheta", test_bc_buff_rtheta },
  { NULL, NULL },
};
