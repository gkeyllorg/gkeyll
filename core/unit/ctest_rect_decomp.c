#include <acutest.h>
#include <gkyl_rect_decomp.h>

#include <string.h>

void test_ranges_1d()
{
  double lower[] = { 1.0 }, upper[] = {2.5 };
  int cells[] = { 20 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);

  int nghost[] = { 1 };
  struct gkyl_range ext_range, range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  TEST_CHECK( ext_range.ndim == 1 );
  TEST_CHECK( ext_range.volume ==  22 );
  TEST_CHECK( ext_range.lower[0] == 0 );
  TEST_CHECK( ext_range.upper[0] == 21 );

  TEST_CHECK( range.ndim == 1 );
  TEST_CHECK( range.volume ==  20 );
  TEST_CHECK( range.lower[0] == 1 );
  TEST_CHECK( range.upper[0] == 20 );

  TEST_CHECK( gkyl_range_is_sub_range(&range) == 1 );
}

void test_ranges_2d()
{
  double lower[] = { 1.0, 1.0 }, upper[] = { 2.5, 5.0 };
  int cells[] = { 20, 40 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  int nghost[] = { 1, 0 };
  struct gkyl_range ext_range, range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  TEST_CHECK( ext_range.ndim == 2 );
  TEST_CHECK( ext_range.volume ==  22*40 );
  TEST_CHECK( ext_range.lower[0] == 0 );
  TEST_CHECK( ext_range.upper[0] == 21 );
  TEST_CHECK( ext_range.lower[1] == 1 );
  TEST_CHECK( ext_range.upper[1] == 40 );  

  TEST_CHECK( range.ndim == 2 );
  TEST_CHECK( range.volume ==  20*40 );
  TEST_CHECK( range.lower[0] == 1 );
  TEST_CHECK( range.upper[0] == 20 );
  TEST_CHECK( range.lower[1] == 1 );
  TEST_CHECK( range.upper[1] == 40 );  

  TEST_CHECK( gkyl_range_is_sub_range(&range) == 1 );  
}

void test_ranges_3d()
{
  double lower[] = { 1.0, 1.0, 1.0 }, upper[] = { 2.5, 5.0, 2.0 };
  int cells[] = { 20, 40, 10 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 3, lower, upper, cells);

  int nghost[] = { 1, 0, 2 };
  struct gkyl_range ext_range, range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  TEST_CHECK( ext_range.ndim == 3 );
  TEST_CHECK( ext_range.volume ==  22*40*14 );
  TEST_CHECK( ext_range.lower[0] == 0 );
  TEST_CHECK( ext_range.upper[0] == 21 );
  TEST_CHECK( ext_range.lower[1] == 1 );
  TEST_CHECK( ext_range.upper[1] == 40 );
  TEST_CHECK( ext_range.lower[2] == -1 );
  TEST_CHECK( ext_range.upper[2] == 12 ); 

  TEST_CHECK( range.ndim == 3 );
  TEST_CHECK( range.volume ==  20*40*10 );
  TEST_CHECK( range.lower[0] == 1 );
  TEST_CHECK( range.upper[0] == 20 );
  TEST_CHECK( range.lower[1] == 1 );
  TEST_CHECK( range.upper[1] == 40 );
  TEST_CHECK( range.lower[2] == 1 );
  TEST_CHECK( range.upper[2] == 10 );  

  TEST_CHECK( gkyl_range_is_sub_range(&range) == 1 );

}

static void
test_ranges_from_range_2d(void)
{
  struct gkyl_range inlocal;
  gkyl_range_init(&inlocal, 2, (int[]) { 1, 2 }, (int[]) { 10, 20 });

  struct gkyl_range local, local_ext;
  gkyl_create_ranges(&inlocal, (int[]) { 2, 1 }, &local_ext, &local);

  TEST_CHECK( local.ndim == inlocal.ndim );
  TEST_CHECK( local_ext.ndim == inlocal.ndim );
  
  for (int i=0; i<inlocal.ndim; ++i) {
    TEST_CHECK( local.lower[i] == inlocal.lower[i] );
    TEST_CHECK( local.upper[i] == inlocal.upper[i] );
  }

  TEST_CHECK( gkyl_range_is_sub_range(&local) == 1 );
}

static void
test_ranges_from_range_3d(void)
{
  struct gkyl_range inlocal;
  gkyl_range_init(&inlocal, 3, (int[]) { 1, 2, 3 }, (int[]) { 10, 20, 30 });

  struct gkyl_range local, local_ext;
  gkyl_create_ranges(&inlocal, (int[]) { 2, 1, 0 }, &local_ext, &local);

  TEST_CHECK( local.ndim == inlocal.ndim );
  TEST_CHECK( local_ext.ndim == inlocal.ndim );
  
  for (int i=0; i<inlocal.ndim; ++i) {
    TEST_CHECK( local.lower[i] == inlocal.lower[i] );
    TEST_CHECK( local.upper[i] == inlocal.upper[i] );
  }

  TEST_CHECK( gkyl_range_is_sub_range(&local) == 1 );
}

// some helper functions
static bool
is_on_corner(int ndim, const int *idx, const int *shape)
{
  bool isc = true;
  for (int i=0; i<ndim; ++i) {
    if ( !((idx[i] == 0) || (idx[i] == shape[i]-1)) )
      isc = false;
  }
  return isc;
}
static bool
is_on_dir_edge(int ndim, int dir, const int *idx, const int *shape)
{
  if ((idx[dir] == 0) || (idx[dir] == shape[dir]-1)) { // on a face
    for (int d=0; d<ndim; ++d) {
      if (d != dir)
        if ((idx[d] == 0) || (idx[d] == shape[d]-1))
          return true;
    }
  }
  return false;
}
static bool
is_on_edge(int ndim, const int *idx, const int *shape)
{
  for (int i=0; i<ndim; ++i)
  {
    if ((idx[i] == 0) || (idx[i] == shape[i]-1)) { // on a face
      for (int d=0; d<ndim; ++d) {
        if (d != i)
          if ((idx[d] == 0) || (idx[d] == shape[d]-1))
            return true;
      }
    }
  }
  return false;
}
static bool
is_on_face(int ndim, const int *idx, const int *shape)
{
  for (int i=0; i<ndim; ++i) {
    if ((idx[i] == 0) || (idx[i] == shape[i]-1))
      return true;
  }
  return false;
}

static void
test_rect_decomp_2d(void)
{
  struct gkyl_range range;
  gkyl_range_init(&range, 2, (int[]) { 1, 2 }, (int[]) { 100, 100 });
  
  int cuts[] = { 5, 6 };
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(2, cuts, &range);

  TEST_CHECK( decomp->ndim == 2 );
  TEST_CHECK( decomp->ndecomp == cuts[0]*cuts[1] );

  TEST_CHECK( gkyl_range_compare(&range, &decomp->parent_range) );

  long vol = 0;
  for (int i=0; i<decomp->ndecomp; ++i)
    vol += decomp->ranges[i].volume;

  TEST_CHECK( vol == range.volume );
  TEST_CHECK( gkyl_rect_decomp_check_covering(decomp) );

  long offs = 0;
  for (int i=0; i<decomp->ndecomp; ++i) {
    TEST_CHECK( offs == gkyl_rect_decomp_calc_offset(decomp, i) );
    offs += decomp->ranges[i].volume;
  }

  struct gkyl_range crange;
  gkyl_range_init_from_shape(&crange, 2, cuts);
  struct gkyl_range_iter iter;

  // check decomposition without corner neighbors
  gkyl_range_iter_init(&iter, &crange);
  while ( gkyl_range_iter_next(&iter) ) {

    struct gkyl_rect_decomp_neigh *neigh =
      gkyl_rect_decomp_calc_neigh(decomp, false, gkyl_range_idx(&crange, iter.idx));
    
    if (is_on_corner(2, iter.idx, cuts))
      TEST_CHECK( neigh->num_neigh == 2 );
    else if (is_on_face(2, iter.idx, cuts) )
      TEST_CHECK( neigh->num_neigh == 3 );
    else
      TEST_CHECK( neigh->num_neigh == 4 );

    gkyl_rect_decomp_neigh_release(neigh);
  }

  // check decomposition with corner neighbors
  gkyl_range_iter_init(&iter, &crange);
  while ( gkyl_range_iter_next(&iter) ) {

    struct gkyl_rect_decomp_neigh *neigh =
      gkyl_rect_decomp_calc_neigh(decomp, true, gkyl_range_idx(&crange, iter.idx));

    if (is_on_corner(2, iter.idx, cuts))
      TEST_CHECK( neigh->num_neigh == 3 );
    else if (is_on_face(2, iter.idx, cuts) )
      TEST_CHECK( neigh->num_neigh == 5 );
    else
      TEST_CHECK( neigh->num_neigh == 8 );

    gkyl_rect_decomp_neigh_release(neigh);
  }

  // check dir and edge without corners
  for (int i=0; i<decomp->ndecomp; ++i) {
    struct gkyl_rect_decomp_neigh *neigh =
      gkyl_rect_decomp_calc_neigh(decomp, false, i);

    for (int n=0; n<neigh->num_neigh; ++n) {
      struct gkyl_range_dir_edge dir_ed =
        gkyl_range_edge_match(&decomp->ranges[i],
          &decomp->ranges[neigh->neigh[n]]);
      
      TEST_CHECK( dir_ed.dir == neigh->dir[n] );
      TEST_CHECK( dir_ed.eloc == neigh->edge[n] );
    }

    gkyl_rect_decomp_neigh_release(neigh);
  }

  // Check the method to get cuts.
  int cuts_new[decomp->ndim];
  gkyl_rect_decomp_get_cuts(decomp, cuts_new);
  for (int d=0; d<decomp->ndim; d++)
    TEST_CHECK( cuts_new[d] == cuts[d] );

  gkyl_rect_decomp_release(decomp);
}

static void
test_rect_decomp_3d(void)
{
  struct gkyl_range range;
  gkyl_range_init(&range, 3, (int[]) { 1, 2, 3 }, (int[]) { 100, 200, 300 });
  
  int cuts[] = { 5, 6, 7 };
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(3, cuts, &range);

  TEST_CHECK( decomp->ndim == 3 );
  TEST_CHECK( decomp->ndecomp == cuts[0]*cuts[1]*cuts[2] );

  TEST_CHECK( gkyl_range_compare(&range, &decomp->parent_range) );

  long vol = 0;
  for (int i=0; i<decomp->ndecomp; ++i)
    vol += decomp->ranges[i].volume;

  TEST_CHECK( vol == range.volume );
  TEST_CHECK( gkyl_rect_decomp_check_covering(decomp) );

  struct gkyl_range crange;
  gkyl_range_init_from_shape(&crange, 3, cuts);
  struct gkyl_range_iter iter;

  // check decomposition without corner neighbors
  gkyl_range_iter_init(&iter, &crange);
  while ( gkyl_range_iter_next(&iter) ) {

    struct gkyl_rect_decomp_neigh *neigh =
      gkyl_rect_decomp_calc_neigh(decomp, false, gkyl_range_idx(&crange, iter.idx));
    
    if (is_on_corner(3, iter.idx, cuts))
      TEST_CHECK( neigh->num_neigh == 3 );
    else if (is_on_edge(3, iter.idx, cuts) )
      TEST_CHECK( neigh->num_neigh == 4 );
    else if (is_on_face(3, iter.idx, cuts) )
      TEST_CHECK( neigh->num_neigh == 5 );
    else
      TEST_CHECK( neigh->num_neigh == 6 );

    gkyl_rect_decomp_neigh_release(neigh);
  }

  // check decomposition with corner neighbors
  gkyl_range_iter_init(&iter, &crange);
  while ( gkyl_range_iter_next(&iter) ) {

    struct gkyl_rect_decomp_neigh *neigh =
      gkyl_rect_decomp_calc_neigh(decomp, true, gkyl_range_idx(&crange, iter.idx));
    
    if (is_on_corner(3, iter.idx, cuts))
      TEST_CHECK( neigh->num_neigh == 7 );
    else if (is_on_edge(3, iter.idx, cuts) )
      TEST_CHECK( neigh->num_neigh == 11 );
    else if (is_on_face(3, iter.idx, cuts) )
      TEST_CHECK( neigh->num_neigh == 17 );
    else
      TEST_CHECK( neigh->num_neigh == 26 );

    gkyl_rect_decomp_neigh_release(neigh);
  }  

  // check dir and edge without corners
  for (int i=0; i<decomp->ndecomp; ++i) {
    struct gkyl_rect_decomp_neigh *neigh =
      gkyl_rect_decomp_calc_neigh(decomp, false, i);

    for (int n=0; n<neigh->num_neigh; ++n) {
      struct gkyl_range_dir_edge dir_ed =
        gkyl_range_edge_match(&decomp->ranges[i],
          &decomp->ranges[neigh->neigh[n]]);
      
      TEST_CHECK( dir_ed.dir == neigh->dir[n] );
      TEST_CHECK( dir_ed.eloc == neigh->edge[n] );
    }

    gkyl_rect_decomp_neigh_release(neigh);
  }  
  
  // Check the method to get cuts.
  int cuts_new[decomp->ndim];
  gkyl_rect_decomp_get_cuts(decomp, cuts_new);
  for (int d=0; d<decomp->ndim; d++)
    TEST_CHECK( cuts_new[d] == cuts[d] );

  gkyl_rect_decomp_release(decomp);
}

static void
test_rect_decomp_4d(void)
{
  struct gkyl_range range;
  gkyl_range_init(&range, 4, (int[]) { 1, 2, 3, 4 }, (int[]) { 10, 20, 30, 40 });
  
  int cuts[] = { 2, 1, 3, 4 };
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(4, cuts, &range);

  TEST_CHECK( decomp->ndim == 4 );
  TEST_CHECK( decomp->ndecomp == cuts[0]*cuts[1]*cuts[2]*cuts[3] );

  TEST_CHECK( gkyl_range_compare(&range, &decomp->parent_range) );

  long vol = 0;
  for (int i=0; i<decomp->ndecomp; ++i)
    vol += decomp->ranges[i].volume;

  TEST_CHECK( vol == range.volume );
  TEST_CHECK( gkyl_rect_decomp_check_covering(decomp) );

  // Check the method to get cuts.
  int cuts_new[decomp->ndim];
  gkyl_rect_decomp_get_cuts(decomp, cuts_new);
  for (int d=0; d<decomp->ndim; d++)
    TEST_CHECK( cuts_new[d] == cuts[d] );

  gkyl_rect_decomp_release(decomp);
}

static void
test_rect_decomp_per_2d(void)
{
  struct gkyl_range range;
  gkyl_range_init(&range, 2, (int[]) { 1, 2 }, (int[]) { 100, 100 });
  
  int cuts[GKYL_MAX_DIM] = { 3, 3 };
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(2, cuts, &range);

  struct gkyl_range crange;
  gkyl_range_init_from_shape(&crange, 2, cuts);  

  struct gkyl_range_iter iter;

  gkyl_range_iter_init(&iter, &crange);
  while ( gkyl_range_iter_next(&iter) ) {
    
    for (int d=0; d<range.ndim; ++d) {
      struct gkyl_rect_decomp_neigh *neigh = gkyl_rect_decomp_calc_periodic_neigh(decomp,
        d, false, gkyl_range_idx(&crange, iter.idx));

      if (is_on_dir_edge(range.ndim, d, iter.idx, cuts)) {
        TEST_CHECK( neigh->num_neigh == 1 );
      }
      
      gkyl_rect_decomp_neigh_release(neigh);
    }
  }

  // Check the method to get cuts.
  int cuts_new[decomp->ndim];
  gkyl_rect_decomp_get_cuts(decomp, cuts_new);
  for (int d=0; d<decomp->ndim; d++)
    TEST_CHECK( cuts_new[d] == cuts[d] );

  gkyl_rect_decomp_release(decomp);
}

static void
test_rect_decomp_per_2d_2(void)
{
  struct gkyl_range range;
  gkyl_range_init(&range, 2, (int[]) { 1, 2 }, (int[]) { 100, 100 });
  
  int cuts[] = { 1, 2 };
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(2, cuts, &range);

  struct gkyl_range crange;
  gkyl_range_init_from_shape(&crange, 2, cuts);  

  struct gkyl_range_iter iter;

  gkyl_range_iter_init(&iter, &crange);
  while ( gkyl_range_iter_next(&iter) ) {
    
    for (int d=0; d<range.ndim; ++d) {
      struct gkyl_rect_decomp_neigh *neigh = gkyl_rect_decomp_calc_periodic_neigh(decomp,
        d, false, gkyl_range_idx(&crange, iter.idx));

      if (is_on_dir_edge(range.ndim, d, iter.idx, cuts)) {
        TEST_CHECK( neigh->num_neigh == 1 );
      }
      
      gkyl_rect_decomp_neigh_release(neigh);
    }
  }

  // Check the method to get cuts.
  int cuts_new[decomp->ndim];
  gkyl_rect_decomp_get_cuts(decomp, cuts_new);
  for (int d=0; d<decomp->ndim; d++)
    TEST_CHECK( cuts_new[d] == cuts[d] );

  gkyl_rect_decomp_release(decomp);
}

static void
test_rect_decomp_per_2d_corner(void)
{
  struct gkyl_range range;
  gkyl_range_init(&range, 2, (int[]) { 1, 2 }, (int[]) { 100, 100 });
  
  int cuts[] = { 2, 2 };
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(2, cuts, &range);

  struct gkyl_range crange;
  gkyl_range_init_from_shape(&crange, 2, cuts);  

  struct gkyl_range_iter iter;

  gkyl_range_iter_init(&iter, &crange);
  while ( gkyl_range_iter_next(&iter) ) {
    
    for (int d=0; d<range.ndim; ++d) {
      struct gkyl_rect_decomp_neigh *neigh = gkyl_rect_decomp_calc_periodic_neigh(decomp,
        d, true, gkyl_range_idx(&crange, iter.idx));

      if (is_on_dir_edge(range.ndim, d, iter.idx, cuts)) {
        TEST_CHECK( neigh->num_neigh == 2 ); // each domain has 2 neighbors
      }
      
      gkyl_rect_decomp_neigh_release(neigh);
    }
  }

  // Check the method to get cuts.
  int cuts_new[decomp->ndim];
  gkyl_rect_decomp_get_cuts(decomp, cuts_new);
  for (int d=0; d<decomp->ndim; d++)
    TEST_CHECK( cuts_new[d] == cuts[d] );

  gkyl_rect_decomp_release(decomp);
}

static void
test_rect_decomp_per_3d(void)
{
  struct gkyl_range range;
  gkyl_range_init(&range, 3, (int[]) { 1, 1, 1 }, (int[]) { 100, 100, 100 });
  
  int cuts[] = { 3, 3, 3 };
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(range.ndim, cuts, &range);

  struct gkyl_range crange;
  gkyl_range_init_from_shape(&crange, range.ndim, cuts);

  struct gkyl_range_iter iter;

  gkyl_range_iter_init(&iter, &crange);
  while ( gkyl_range_iter_next(&iter) ) {
    
    for (int d=0; d<range.ndim; ++d) {
      struct gkyl_rect_decomp_neigh *neigh = gkyl_rect_decomp_calc_periodic_neigh(decomp,
        d, false, gkyl_range_idx(&crange, iter.idx));

      if (is_on_dir_edge(range.ndim, d, iter.idx, cuts)) {
        TEST_CHECK( neigh->num_neigh == 1 );
      }
      
      gkyl_rect_decomp_neigh_release(neigh);
    }
  }

  // Check the method to get cuts.
  int cuts_new[decomp->ndim];
  gkyl_rect_decomp_get_cuts(decomp, cuts_new);
  for (int d=0; d<decomp->ndim; d++)
    TEST_CHECK( cuts_new[d] == cuts[d] );

  gkyl_rect_decomp_release(decomp);  
}

static void
test_rect_decomp_2d_2v(void)
{
  struct gkyl_range range;
  gkyl_range_init(&range, 2, (int[]) { 1, 2 }, (int[]) { 100, 100 });
  
  int cuts[] = { 5, 6 };
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(2, cuts, &range);

  TEST_CHECK( decomp->ndim == 2 );
  TEST_CHECK( decomp->ndecomp == cuts[0]*cuts[1] );

  struct gkyl_range vrange;
  gkyl_range_init(&vrange, 2, (int[]) { 1, 2 }, (int[]) { 16, 16 });

  struct gkyl_rect_decomp *ext_decomp = gkyl_rect_decomp_extended_new(
    &vrange, decomp);

  TEST_CHECK( ext_decomp->ndim == 4 );
  TEST_CHECK( ext_decomp->ndecomp == cuts[0]*cuts[1] );

  TEST_CHECK( ext_decomp->parent_range.volume = decomp->parent_range.volume*vrange.volume );
  for (int i=0; i<decomp->ndecomp; ++i)
    TEST_CHECK( ext_decomp->ranges[i].volume == decomp->ranges[i].volume*vrange.volume );

  for (int i=0; i<decomp->ndecomp; ++i) {
    
    for (int d=0; d<vrange.ndim; ++d) {
      TEST_CHECK( ext_decomp->ranges[d].lower[range.ndim+d] == vrange.lower[d] );
      TEST_CHECK( ext_decomp->ranges[d].upper[range.ndim+d] == vrange.upper[d] );
    }
  }

  // Check the method to get cuts.
  int cuts_new[decomp->ndim];
  gkyl_rect_decomp_get_cuts(decomp, cuts_new);
  for (int d=0; d<decomp->ndim; d++)
    TEST_CHECK( cuts_new[d] == cuts[d] );

  gkyl_rect_decomp_release(decomp);
  gkyl_rect_decomp_release(ext_decomp);
}

static void
test_rect_decomp_from_cuts_and_cells(void)
{
  int cuts[] = { 5, 6, 7 };
  int cells[] = { 100, 200, 300 };
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts_and_cells(3, cuts, cells);

  TEST_CHECK( decomp->ndim == 3 );
  TEST_CHECK( decomp->ndecomp == cuts[0]*cuts[1]*cuts[2] );

  TEST_CHECK( gkyl_rect_decomp_check_covering(decomp) );

  // Check the method to get cuts.
  int cuts_new[decomp->ndim];
  gkyl_rect_decomp_get_cuts(decomp, cuts_new);
  for (int d=0; d<decomp->ndim; d++)
    TEST_CHECK( cuts_new[d] == cuts[d] );

  gkyl_rect_decomp_release(decomp);
}

TEST_LIST = {
  { "ranges_1d", test_ranges_1d },
  { "ranges_2d", test_ranges_2d },
  { "ranges_3d", test_ranges_3d },

  { "ranges_from_range_2d", test_ranges_from_range_2d },
  { "ranges_from_range_3d", test_ranges_from_range_3d },  

  { "rect_decomp_2d", test_rect_decomp_2d },
  { "rect_decomp_3d", test_rect_decomp_3d },
  { "rect_decomp_4d", test_rect_decomp_4d },

  { "rect_decomp_per_2d", test_rect_decomp_per_2d },
  { "rect_decomp_per_2d_2", test_rect_decomp_per_2d_2 },
  { "rect_decomp_per_3d", test_rect_decomp_per_3d },

  { "rect_decomp_per_2d_corner", test_rect_decomp_per_2d_corner },

  { "rect_decomp_2d_2v", test_rect_decomp_2d_2v },

  { "rect_decomp_from_cuts_and_cells", test_rect_decomp_from_cuts_and_cells },
  
  { NULL, NULL },
};
