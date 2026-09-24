#include <gkyl_alloc.h>
#include <gkyl_mat_triples.h>
#include <gkyl_range.h>

// Index into matrix
struct mat_idx { size_t row, col; };

static inline int
cmp_long(long a, long b)
{
  return a==b ? 0 : a<b ? -1 : 1;
}

static int
(*mat_idx_cmp)(const struct mat_idx *mia, const struct mat_idx *mib);

static inline int
mat_idx_cmp_row(const struct mat_idx *mia, const struct mat_idx *mib)
{
  return cmp_long(mia->row, mib->row) == 0 ? cmp_long(mia->col, mib->col) :
    cmp_long(mia->row, mib->row);
}

static inline int
mat_idx_cmp_col(const struct mat_idx *mia, const struct mat_idx *mib)
{
  return cmp_long(mia->col, mib->col) == 0 ? cmp_long(mia->row, mib->row) :
    cmp_long(mia->col, mib->col);
}

// define map of mat_idx -> gkyl_mtriple: this is a sorted map, in
// which indices are sorted in column order (see mat_idx_cmp method)
#define i_key struct mat_idx
#define i_val struct gkyl_mtriple
#define i_cmp mat_idx_cmp
#define i_tag triple
#include <stc/csmap.h>
// complete definition of map

struct gkyl_mat_triples {
  struct gkyl_range range; // range representing matrix shape
  csmap_triple triples; // non-zero values in matrix
  int ordering;
};

struct gkyl_mat_triples_iter {
  bool is_first; // first time?
  int nrem; /// number of addresses remaining
  const csmap_triple *parent; // pointer to parent triples
  csmap_triple_iter it_curr, it_end; // iterator to start, end of map
};

gkyl_mat_triples*
gkyl_mat_triples_new(size_t nr, size_t nc)
{
  struct gkyl_mat_triples *tri = gkyl_malloc(sizeof(struct gkyl_mat_triples));

  gkyl_range_init_from_shape(&tri->range, 2, (const int[]) { nr, nc} );

  // set column-major order by default
  tri->ordering = COLMAJOR;

  tri->triples = csmap_triple_init();
  return tri;
}

void gkyl_mat_triples_set_rowmaj_order(gkyl_mat_triples *tri)
{
  tri->ordering = ROWMAJOR;
}

void gkyl_mat_triples_set_colmaj_order(gkyl_mat_triples *tri)
{
  tri->ordering = COLMAJOR;
}

bool gkyl_mat_triples_is_rowmaj(gkyl_mat_triples *tri) {
  return tri->ordering == ROWMAJOR;
}

bool gkyl_mat_triples_is_colmaj(gkyl_mat_triples *tri) {
  return tri->ordering == COLMAJOR;
}

GKYL_CU_DH double
gkyl_mat_triples_insert(gkyl_mat_triples *tri, size_t i, size_t j, double val)
{
  assert(i<gkyl_range_shape(&tri->range, 0) && j<gkyl_range_shape(&tri->range, 1));
  if(tri->ordering == COLMAJOR) 
    mat_idx_cmp = mat_idx_cmp_col;
  else
    mat_idx_cmp = mat_idx_cmp_row;
  
  long loc = gkyl_ridx(tri->range, i, j);
  csmap_triple_put(&tri->triples, (struct mat_idx) { .row = i, .col = j },
    (struct gkyl_mtriple) { .row = i, .col = j, .val = val }
  );
  return val;
}

GKYL_CU_DH double
gkyl_mat_triples_accum(gkyl_mat_triples *tri, size_t i, size_t j, double val)
{
  assert(i<gkyl_range_shape(&tri->range, 0) && j<gkyl_range_shape(&tri->range, 1));
  if(tri->ordering == COLMAJOR) 
    mat_idx_cmp = mat_idx_cmp_col;
  else
    mat_idx_cmp = mat_idx_cmp_row;
  
  long loc = gkyl_ridx(tri->range, i, j);
  struct csmap_triple_value *mt = csmap_triple_get_mut(&tri->triples,
    (struct mat_idx) { .row = i, .col = j }
  );
  double tot_val = val;
  if (mt) {
    // element exists, add to its current value
    tot_val = (mt->second.val += val);
  }
  else {
    csmap_triple_put(&tri->triples, (struct mat_idx) { .row = i, .col = j },
      (struct gkyl_mtriple) { .row = i, .col = j, .val = val }
    );
  }
  
  return tot_val;
}

double
gkyl_mat_triples_get(const gkyl_mat_triples *tri, size_t i, size_t j)
{
  long loc = gkyl_ridx(tri->range, i, j);
  const struct csmap_triple_value *mt = csmap_triple_get(&tri->triples,
    (struct mat_idx) { .row = i, .col = j }
  );
  return mt ? mt->second.val : 0.0;
}

size_t
gkyl_mat_triples_size(const gkyl_mat_triples *tri)
{
  return csmap_triple_size(tri->triples);
}

gkyl_mat_triples_iter*
gkyl_mat_triples_iter_new(const gkyl_mat_triples *tri)
{
  struct gkyl_mat_triples_iter *iter = gkyl_malloc(sizeof(*iter));

  iter->parent = &tri->triples;
  
  iter->is_first = true;
  iter->nrem = csmap_triple_size(tri->triples);
  iter->it_curr = csmap_triple_begin(&tri->triples);
  iter->it_end = csmap_triple_end(&tri->triples);
  
  return iter;
}

void
gkyl_mat_triples_iter_init(struct gkyl_mat_triples_iter *iter, const gkyl_mat_triples *tri)
{
  iter->is_first = true;
  iter->nrem = csmap_triple_size(tri->triples);
  iter->it_curr = csmap_triple_begin(&tri->triples);
}

bool
gkyl_mat_triples_iter_next(gkyl_mat_triples_iter *iter)
{
  if (iter->is_first) {
    iter->is_first = false;
    return (iter->nrem-- > 0);
  }
  if (iter->nrem-- > 0) {
    csmap_triple_next(&iter->it_curr);
    return true;
  }
  return false;
}

struct gkyl_mtriple
gkyl_mat_triples_iter_at(const gkyl_mat_triples_iter *iter)
{
  return iter->it_curr.ref->second;
}

void
gkyl_mat_triples_clear(struct gkyl_mat_triples *tri, double val)
{
  gkyl_mat_triples_iter *iter = gkyl_mat_triples_iter_new(tri);
  while (gkyl_mat_triples_iter_next(iter)) {
    struct gkyl_mtriple mt = gkyl_mat_triples_iter_at(iter);

    struct csmap_triple_value *mtm = csmap_triple_get_mut(&tri->triples,
      (struct mat_idx) { .row = mt.row, .col = mt.col });
    mtm->second.val = val;
  }
  gkyl_mat_triples_iter_release(iter);
}

void
gkyl_mat_triples_iter_release(gkyl_mat_triples_iter *iter)
{
  gkyl_free(iter);
}

void
gkyl_mat_triples_release(gkyl_mat_triples *tri)
{
  csmap_triple_drop(&tri->triples);
  gkyl_free(tri);
}
