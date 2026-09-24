#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>

#include <string.h>

// define container to store vector of ints
#define i_val int
#include <stc/cvec.h>

// Private struct to manage the neighbor data struct
struct rect_decomp_neigh_cont {
  struct gkyl_rect_decomp_neigh neigh;
  cvec_int l_neigh;
  cvec_int l_dir;
  cvec_int l_edge;
};

static void
rect_decomp_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_rect_decomp *decomp = container_of(ref, struct gkyl_rect_decomp, ref_count);
  gkyl_free(decomp->ranges);
  gkyl_free(decomp);
}    

struct gkyl_rect_decomp*
gkyl_rect_decomp_new_from_cuts(int ndim, const int cuts[], const struct gkyl_range *range)
{
  struct gkyl_rect_decomp *decomp = gkyl_malloc(sizeof(*decomp));

  int ndecomp = 1;
  decomp->ndim = ndim;  

  for (int d=0; d<ndim; ++d) ndecomp *= cuts[d];  
  decomp->ndecomp = ndecomp;
  decomp->ranges = gkyl_malloc(sizeof(struct gkyl_range[ndecomp]));

  memcpy(&decomp->parent_range, range, sizeof(struct gkyl_range));

  div_t qr[GKYL_MAX_DIM];
  for (int d=0; d<ndim; ++d)
    qr[d] = div(gkyl_range_shape(range, d), cuts[d]);

  int *sidx[GKYL_MAX_DIM], *eidx[GKYL_MAX_DIM];
  for (int d=0; d<ndim; ++d) {
    sidx[d] = gkyl_malloc(sizeof(int[cuts[d]]));
    eidx[d] = gkyl_malloc(sizeof(int[cuts[d]]));

    int *shape = gkyl_malloc(sizeof(int[cuts[d]]));

    // compute shape in direction 'd'
    for (int i=0; i<cuts[d]; ++i)
      shape[i] = i<qr[d].rem ? qr[d].quot+1 : qr[d].quot;

    sidx[d][0] = range->lower[d];
    eidx[d][0] = sidx[d][0]+shape[0]-1;
    for (int i=1; i<cuts[d]; ++i) {
      sidx[d][i] = eidx[d][i-1]+1;
      eidx[d][i] = sidx[d][i]+shape[i]-1;
    }

    gkyl_free(shape);
  }

  struct gkyl_range rcuts;
  gkyl_range_init_from_shape(&rcuts, ndim, cuts);
  struct gkyl_range_iter citer;
  gkyl_range_iter_init(&citer, &rcuts);

  int dnum = 0;
  // loop over cuts range, constructing each of the sub-ranges in the
  // decomposition
  while( gkyl_range_iter_next(&citer) ) {
    int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];

    for (int d=0; d<ndim; ++d) {
      lower[d] = sidx[d][citer.idx[d]];
      upper[d] = eidx[d][citer.idx[d]];
    }

    gkyl_range_init(&decomp->ranges[dnum++], range->ndim, lower, upper);
  }

  for (int d=0; d<ndim; ++d) {
    gkyl_free(sidx[d]);
    gkyl_free(eidx[d]);
  }

  decomp->ref_count = gkyl_ref_count_init(rect_decomp_free);
  
  return decomp;
}

struct gkyl_rect_decomp*
gkyl_rect_decomp_new_from_cuts_and_cells(int ndim, const int cuts[], const int cells[])
{
  struct gkyl_range range;
  gkyl_create_global_range(ndim, cells, &range);
  return gkyl_rect_decomp_new_from_cuts(ndim, cuts, &range);
}

// ext_range = a X b 
static void
init_extend_range(struct gkyl_range *ext_range,
  const struct gkyl_range *a, const struct gkyl_range *b)
{
  int adim = a->ndim, bdim = b->ndim;
  int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];

  for (int d=0; d<adim; ++d) {
    lower[d] = a->lower[d];
    upper[d] = a->upper[d];
  }
  for (int d=0; d<bdim; ++d) {
    lower[adim+d] = b->lower[d];
    upper[adim+d] = b->upper[d];
  }

  gkyl_range_init(ext_range, adim+bdim, lower, upper);
}

struct gkyl_rect_decomp*
gkyl_rect_decomp_extended_new(const struct gkyl_range *arange,
  const struct gkyl_rect_decomp *decomp)
{
  struct gkyl_rect_decomp *extd = gkyl_malloc(sizeof(*extd));

  int ndecomp =  extd->ndecomp = decomp->ndecomp;  
  int ndim = extd->ndim = arange->ndim + decomp->ndim;
  extd->ranges = gkyl_malloc(sizeof(struct gkyl_range[ndecomp]));

  gkyl_range_ten_prod(&extd->parent_range, &decomp->parent_range, arange);
  for (int n=0; n<ndecomp; ++n)
    gkyl_range_ten_prod(&extd->ranges[n], &decomp->ranges[n], arange);

  extd->ref_count = gkyl_ref_count_init(rect_decomp_free);
  
  return extd;
}

struct gkyl_rect_decomp*
gkyl_rect_decomp_acquire(const struct gkyl_rect_decomp *decomp)
{
  gkyl_ref_count_inc(&decomp->ref_count);
  return (struct gkyl_rect_decomp*) decomp;
}

bool
gkyl_rect_decomp_check_covering(const struct gkyl_rect_decomp *decomp)
{
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 1, decomp->parent_range.volume);
  gkyl_array_clear(arr, 0.0);

  // following loops over each sub-range and increments the region it
  // indexes in 'arr'. Each index should be visited exactly once.
  for (int i=0; i<decomp->ndecomp; ++i) {
    // construct a sub-range so indexing into global array works fine
    struct gkyl_range lrange;
    gkyl_sub_range_intersect(&lrange, &decomp->parent_range, &decomp->ranges[i]);
    
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &lrange);

    while (gkyl_range_iter_next(&iter)) {
      double *d = gkyl_array_fetch(arr, gkyl_range_idx(&decomp->parent_range, iter.idx));
      d[0] += 1.0;
    }
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &decomp->parent_range);
  while (gkyl_range_iter_next(&iter)) {
    const double *d = gkyl_array_cfetch(arr, gkyl_range_idx(&decomp->parent_range, iter.idx));
    if (d[0] != 1.0)
      return false;
  }

  gkyl_array_release(arr);
  
  return true;
}

// compute neighbors accounting for corner neighbors
static struct gkyl_rect_decomp_neigh*
calc_neigh_with_corners(const struct gkyl_rect_decomp *decomp, int nidx)
{
  struct rect_decomp_neigh_cont *cont = gkyl_malloc(sizeof(*cont));
  cont->l_neigh = cvec_int_init();
  cont->l_dir = cvec_int_init();
  cont->l_edge = cvec_int_init();
  
  int elo[GKYL_MAX_DIM], eup[GKYL_MAX_DIM];
  for (int i=0; i<decomp->ndim; ++i)
    elo[i] = eup[i] = 1;
  
  struct gkyl_range erng;
  gkyl_range_extend(&erng, &decomp->ranges[nidx], elo, eup);

  for (int i=0; i<decomp->ndecomp; ++i)
    if (i != nidx) {
      struct gkyl_range irng;
      int is_inter = gkyl_range_intersect(&irng, &erng,
        &decomp->ranges[i]);
      if (is_inter) {
        cvec_int_push_back(&cont->l_neigh, i);
        
        struct gkyl_range_dir_edge dir_ed =
          gkyl_range_edge_match(&decomp->ranges[nidx], &decomp->ranges[i]);

        cvec_int_push_back(&cont->l_dir, dir_ed.dir);
        cvec_int_push_back(&cont->l_edge, dir_ed.eloc);
      }
    }

  cont->neigh.num_neigh = cvec_int_size(cont->l_neigh);
  cont->neigh.neigh = cvec_int_front(&cont->l_neigh);
  cont->neigh.dir = cvec_int_front(&cont->l_dir);
  cont->neigh.edge = cvec_int_front(&cont->l_edge);
  
  return &cont->neigh;
}

// compute neighbors leaving out corner neighbors: only face neighbors
// are included
static struct gkyl_rect_decomp_neigh*
calc_neigh_no_corners(const struct gkyl_rect_decomp *decomp, int nidx)
{
  struct rect_decomp_neigh_cont *cont = gkyl_malloc(sizeof(*cont));
  cont->l_neigh = cvec_int_init();
  cont->l_dir = cvec_int_init();
  cont->l_edge = cvec_int_init();  
  
  struct gkyl_range erng;

  for (int n=0; n<decomp->ndim; ++n) {
    
    int elo[GKYL_MAX_DIM] = { 0 }, eup[GKYL_MAX_DIM] = { 0 };
    elo[n] = eup[n] = 1; // only extend in 1 direction
    gkyl_range_extend(&erng, &decomp->ranges[nidx], elo, eup);

    for (int i=0; i<decomp->ndecomp; ++i)
      if (i != nidx) {
        struct gkyl_range irng;
        int is_inter = gkyl_range_intersect(&irng, &erng,
          &decomp->ranges[i]);
        if (is_inter) {
          cvec_int_push_back(&cont->l_neigh, i);

          struct gkyl_range_dir_edge dir_ed =
            gkyl_range_edge_match(&decomp->ranges[nidx], &decomp->ranges[i]);
          
          cvec_int_push_back(&cont->l_dir, dir_ed.dir);
          cvec_int_push_back(&cont->l_edge, dir_ed.eloc);
        }
      }
  }
  
  cont->neigh.num_neigh = cvec_int_size(cont->l_neigh);
  cont->neigh.neigh = cvec_int_front(&cont->l_neigh);
  cont->neigh.dir = cvec_int_front(&cont->l_dir);
  cont->neigh.edge = cvec_int_front(&cont->l_edge);  
  
  return &cont->neigh;
}

struct gkyl_rect_decomp_neigh*
gkyl_rect_decomp_calc_neigh(const struct gkyl_rect_decomp *decomp,
  bool inc_corners, int nidx)
{
  if (inc_corners)
    return calc_neigh_with_corners(decomp, nidx);
  return calc_neigh_no_corners(decomp, nidx);
}

struct gkyl_rect_decomp_neigh*
gkyl_rect_decomp_calc_periodic_neigh(const struct gkyl_rect_decomp *decomp,
  int dir, bool inc_corners, int nidx)
{
  struct rect_decomp_neigh_cont *cont = gkyl_malloc(sizeof(*cont));
  cont->l_neigh = cvec_int_init();
  cont->l_dir = cvec_int_init();
  cont->l_edge = cvec_int_init();  

  const struct gkyl_range *curr = &decomp->ranges[nidx];

  int elo[GKYL_MAX_DIM] = { 0 }, eup[GKYL_MAX_DIM] = { 0 };
  if (inc_corners)
    for (int i=0; i<decomp->ndim; ++i)
      elo[i] = eup[i] = 1;
  else
    elo[dir] = eup[dir] = 1;
  
  if (gkyl_range_is_on_lower_edge(dir, curr, &decomp->parent_range)) {
    int delta[GKYL_MAX_DIM] = { 0 };
    delta[dir] = gkyl_range_shape(&decomp->parent_range, dir);
    
    struct gkyl_range curr_shift;
    gkyl_range_shift(&curr_shift, curr, delta);
      
    struct gkyl_range shift_erng;
    gkyl_range_extend(&shift_erng, &curr_shift, elo, eup);

    for (int i=0; i<decomp->ndecomp; ++i)
      if (gkyl_range_is_on_upper_edge(dir, &decomp->ranges[i], &decomp->parent_range)) {
        struct gkyl_range irng;
        int is_inter = gkyl_range_intersect(&irng, &shift_erng,
          &decomp->ranges[i]);
        if (is_inter) {
          cvec_int_push_back(&cont->l_neigh, i);
          cvec_int_push_back(&cont->l_dir, dir);
          // this is not exactly correct: corner neighbors should not
          // be on any edge
          cvec_int_push_back(&cont->l_edge, GKYL_LOWER_EDGE);
        }
      }
  }
  else if (gkyl_range_is_on_upper_edge(dir, curr, &decomp->parent_range)) {
    int delta[GKYL_MAX_DIM] = { 0 };
    delta[dir] = -gkyl_range_shape(&decomp->parent_range, dir);
    
    struct gkyl_range curr_shift;
    gkyl_range_shift(&curr_shift, curr, delta);
      
    struct gkyl_range shift_erng;
    gkyl_range_extend(&shift_erng, &curr_shift, elo, eup);

    for (int i=0; i<decomp->ndecomp; ++i)
      if (gkyl_range_is_on_lower_edge(dir, &decomp->ranges[i], &decomp->parent_range)) {
        struct gkyl_range irng;
        int is_inter = gkyl_range_intersect(&irng, &shift_erng,
          &decomp->ranges[i]);
        if (is_inter) {
          cvec_int_push_back(&cont->l_neigh, i);
          cvec_int_push_back(&cont->l_dir, dir);
          // this is not exactly correct: corner neighbors should not
          // be on any edge          
          cvec_int_push_back(&cont->l_edge, GKYL_UPPER_EDGE);
        }
      }
  }

  cont->neigh.num_neigh = cvec_int_size(cont->l_neigh);  
  cont->neigh.neigh = cvec_int_front(&cont->l_neigh);
  cont->neigh.dir = cvec_int_front(&cont->l_dir);
  cont->neigh.edge = cvec_int_front(&cont->l_edge);  
  
  return &cont->neigh;
}

void
gkyl_rect_decomp_neigh_release(struct gkyl_rect_decomp_neigh *ng)
{
  struct rect_decomp_neigh_cont *cont = container_of(ng,
    struct rect_decomp_neigh_cont, neigh);
  cvec_int_drop(&cont->l_neigh);
  cvec_int_drop(&cont->l_dir);
  cvec_int_drop(&cont->l_edge);
  gkyl_free(cont);
}

long
gkyl_rect_decomp_calc_offset(const struct gkyl_rect_decomp *decomp, int nidx)
{
  long offset = 0;
  for (int i=0; i<nidx; ++i)
    offset += decomp->ranges[i].volume;
  return offset;
}

void      
gkyl_rect_decomp_release(struct gkyl_rect_decomp *decomp)
{
  gkyl_ref_count_dec(&decomp->ref_count);
}

// Utility functions

void
gkyl_create_global_range(int ndim, const int *cells, struct gkyl_range *range)
{
  int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  for (int i=0; i<ndim; ++i) {
    // this needs to be consistent with gkyl_create_grid_ranges below
    lower[i] = 1;
    upper[i] = cells[i];
  }
  gkyl_range_init(range, ndim, lower, upper);
}

void
gkyl_create_grid_ranges(const struct gkyl_rect_grid *grid,
  const int *nghost, struct gkyl_range *ext_range, struct gkyl_range *range)
{
  int lower_ext[GKYL_MAX_DIM], upper_ext[GKYL_MAX_DIM];
  int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  
  for (int i=0; i<grid->ndim; ++i) {
    lower_ext[i] = 1-nghost[i];
    upper_ext[i] = grid->cells[i]+nghost[i];

    // this needs to be consistent with gkyl_create_global_range above
    lower[i] = 1;
    upper[i] = grid->cells[i];
  }
  gkyl_range_init(ext_range, grid->ndim, lower_ext, upper_ext);
  gkyl_sub_range_init(range, ext_range, lower, upper);
}

void
gkyl_create_ranges(const struct gkyl_range *inrange,
  const int *nghost, struct gkyl_range *ext_range, struct gkyl_range *range)
{
  int lower_ext[GKYL_MAX_DIM], upper_ext[GKYL_MAX_DIM];
  int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  
  for (int i=0; i<inrange->ndim; ++i) {
    lower_ext[i] = inrange->lower[i]-nghost[i];
    upper_ext[i] = inrange->upper[i]+nghost[i];

    lower[i] = inrange->lower[i];
    upper[i] = inrange->upper[i];
  }
  gkyl_range_init(ext_range, inrange->ndim, lower_ext, upper_ext);
  gkyl_sub_range_init(range, ext_range, lower, upper);  
}

void
gkyl_create_vertex_ranges(const struct gkyl_range *inrange,
  const int *nghost, struct gkyl_range *ext_range, struct gkyl_range *range)
{
  int lower_ext[GKYL_MAX_DIM], upper_ext[GKYL_MAX_DIM];
  int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];

  for (int i=0; i<inrange->ndim; ++i) {
    lower_ext[i] = inrange->lower[i] - nghost[i];
    upper_ext[i] = inrange->upper[i] + 1 + nghost[i];

    lower[i] = inrange->lower[i];
    upper[i] = inrange->upper[i] + 1;
  }
  gkyl_range_init(ext_range, inrange->ndim, lower_ext, upper_ext);
  gkyl_sub_range_init(range, ext_range, lower, upper);  
}

void
gkyl_rect_decomp_get_cuts(struct gkyl_rect_decomp* decomp, int* cuts)
{
  int ndim = decomp->ndim;

  for (int d=0; d<ndim; d++) {
    int other_dim_lo[GKYL_MAX_DIM] = {0}, other_dim_up[GKYL_MAX_DIM] = {0};
    for (int i=0; i<ndim; i++) {
      if (i != d) {
        other_dim_lo[i] = decomp->ranges[0].lower[i];
        other_dim_up[i] = decomp->ranges[0].upper[i];
      }
    }

    int cuts_curr = 0, range_idx = 0;

    bool not_reached_upper = true;
    while (not_reached_upper) {
      struct gkyl_range range_curr = decomp->ranges[range_idx];
      bool same_other_lims = true;
      for (int i=0; i<ndim; i++) {
        if (i != d) {
          same_other_lims = same_other_lims &&
            ((range_curr.lower[i] == other_dim_lo[i]) && (range_curr.upper[i] == other_dim_up[i]));
        }
      }

      if (same_other_lims) {
        cuts_curr++;
        if (range_curr.upper[d] == decomp->parent_range.upper[d])
          not_reached_upper = false;
      }

      range_idx++;
    }
    cuts[d] = cuts_curr;
  }
}
