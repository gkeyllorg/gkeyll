#pragma once

#include <gkyl_util.h>
#include <gkyl_vargm.h>

#include <stdint.h>
#include <stdio.h>

/**
 * Series of indexing "functions" to compute linear index into range
 */
#define gkyl_ridx1(r, i1)                       \
    ((r).ac[0]+(i1)*(r).ac[1])
#define gkyl_ridx2(r, i1, i2)                   \
    ((r).ac[0]+((i1)*(r).ac[1]+(i2)*(r).ac[2]))
#define gkyl_ridx3(r, i1, i2, i3)                                       \
    (((r).ac[0]+(i1)*(r).ac[1])+((i2)*(r).ac[2]+(i3)*(r).ac[3]))
#define gkyl_ridx4(r, i1, i2, i3, i4)                                   \
    (((r).ac[0]+(i1)*(r).ac[1])+((i2)*(r).ac[2]+(i3)*(r).ac[3]+(i4)*(r).ac[4]))
#define gkyl_ridx5(r, i1, i2, i3, i4, i5)                               \
    (((r).ac[0]+(i1)*(r).ac[1])+((i2)*(r).ac[2]+(i3)*(r).ac[3])+((i4)*(r).ac[4]+(i5)*(r).ac[5]))
#define gkyl_ridx6(r, i1, i2, i3, i4, i5, i6)                           \
    (((r).ac[0]+(i1)*(r).ac[1])+((i2)*(r).ac[2]+(i3)*(r).ac[3])+((i4)*(r).ac[4]+(i5)*(r).ac[5]+(i6)*(r).ac[6]))
#define gkyl_ridx7(r, i1, i2, i3, i4, i5, i6, i7)                       \
    (((r).ac[0]+(i1)*(r).ac[1])+((i2)*(r).ac[2]+(i3)*(r).ac[3])+((i4)*(r).ac[4]+(i5)*(r).ac[5]+(i6)*(r).ac[6]) + (i7)*(r).ac[7])

/** Generic indexing: works for 1D-7D (VFUNC1 is defined-ed in
 * gkyl_vargm.h) */
#define gkyl_ridx(r, ...) VFUNC1(gkyl_ridx, r, __VA_ARGS__)

/** Indexing macro taking index defined as array of int */
#define gkyl_ridxn(r, idx) gkyl_range_idx(&(r), idx)

// Constants to represent lower/upper edges
enum gkyl_edge_loc { GKYL_LOWER_EDGE = 0, GKYL_UPPER_EDGE, GKYL_NO_EDGE };

// Direction and location of range
struct gkyl_range_dir_edge {
  int dir;
  enum gkyl_edge_loc eloc;
};  

/**
 * Range object, representing an N-dimensional integer index
 * set. Lower and upper limits are inclusive.
 */
struct gkyl_range {
  int ndim; // number of dimension
  int lower[GKYL_MAX_DIM]; // lower bound
  int upper[GKYL_MAX_DIM]; // upper bound (inclusive)
  long volume; // total volume of range
    
  // do not access directly
  uint32_t flags; // Flags for internal use
  int ilo[GKYL_MAX_DIM]; // for use in inverse indexer
  long ac[GKYL_MAX_DIM+1]; // coefficients for indexing
  long iac[GKYL_MAX_DIM+1]; // for use in sub-range inverse indexer
  long linIdxZero; // linear index of {0,0,...}
  int nsplit, tid; // number of splits, split ID

  // FOR CUDA ONLY
  int nthreads, nblocks; // CUDA kernel launch specifiers for range-based ops
};

/**
 * Iterator object into the range. You can read the 'idx' pointer but
 * must not modify it or any other members of this struct.
 */
struct gkyl_range_iter {
  int idx[GKYL_MAX_DIM]; // current index (do not modify)

  // do not access
  int is_first, ndim;
  long bumps_left;
  int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
};

/**
 * Skip-list based iterator object.
 */
struct gkyl_range_skip_iter {
  long delta; // number of contiguous elements
  struct gkyl_range range; // outer range for iteration
};

/**
 * Initialize a range object.
 *
 * @param rng Range object to initialize
 * @param ndim Dimension of range to create.
 * @param lower Lower indices of range
 * @param upper Upper indices of range
 */
void gkyl_range_init(struct gkyl_range *rng, int ndim,
  const int *lower, const int *upper);

/**
 * Create and initialize a new range object. Free with gkyl_range_release.
 *
 * @param ndim Dimension of range to create.
 * @param lower Lower indices of range
 * @param upper Upper indices of range
 * @return Pointer to new range object.
 */
struct gkyl_range* gkyl_range_new(int ndim, const int *lower, const int *upper);

/**
 * Create new range object from specified shape. This sets the lower
 * indices to [0,...] and upper indices to [shape[0]-1, ...].
 *
 * @param rng Range object to initialize
 * @param ndim Dimensiom of range to create.
 * @param shape Shape of region
 */

void gkyl_range_init_from_shape(struct gkyl_range *rng, int ndim,
  const int *shape);

/**
 * Create new range object from specified shape. This sets the lower
 * indices to [1,...] and upper indices to [shape[0], ...].
 *
 * @param rng Range object to initialize
 * @param ndim Dimensiom of range to create.
 * @param shape Shape of region
 */

void gkyl_range_init_from_shape1(struct gkyl_range *rng, int ndim,
  const int *shape);

/**
 * Create a new range which is a tensor product of @a a and @a b input
 * ranges.
 *
 * @param rng On output, rng = a X b
 * @param a First operand of tensor-product
 * @param b Second operand of tensor-product
 */
void gkyl_range_ten_prod(struct gkyl_range *rng, const struct gkyl_range *a,
  const struct gkyl_range *b);

/**
 * Create a new range that is the same shape as inp range, but the
 * indices are shifted in each direction by delta[dir]
 *
 * @param rng On output new shifted range
 * @param inp Input range to shift
 * @param delta Range indices are shifted by delta[dir] in each direction
 */
void gkyl_range_shift(struct gkyl_range *rng, const struct gkyl_range *inp,
  const int *delta);

/**
 * Create a new range that is the same shape as inp range, but lower
 * indices are reset to the specified ones.
 *
 * @param rng On output new reset range
 * @param inp Input range to reset
 * @param lower New lower indices
 */
void gkyl_range_reset_lower(struct gkyl_range *rng, const struct gkyl_range *inp,
  const int *lower);

/**
 * Shape in direction dir
 *
 * @param rng Range object
 * @param dir Direction to compute shape
 * @return Shape in direction dit
 */
GKYL_CU_DH
static inline int gkyl_range_shape(const struct gkyl_range *rng, int dir)
{
  return rng->upper[dir]-rng->lower[dir]+1;  
}

/**
 * Return 1 if range is a sub-range.
 *
 * @param rng Range object
 * @return 1 if true, 0 otherwise
 */
int gkyl_range_is_sub_range(const struct gkyl_range *rng);

/**
 * Return 1 if idx is inside the range.
 *
 * @param rng Range obkect
 * @return 1 if true, 0 otherwise
 */
GKYL_CU_DH
static inline int gkyl_range_contains_idx(const struct gkyl_range *rng, const int *idx)
{
  for (int i=0; i<rng->ndim; ++i) {
    if ( (idx[i] < rng->lower[i]) || (idx[i] > rng->upper[i]) )
      return 0;
  }
  return 1;
}

/**
 * Create a sub-range from a given range. The sub-range must be fully
 * contained in the parent range or else it will be truncated. The
 * sub-range and the parent range will returns the same linear index
 * for a given index.
 *
 * @param rng New range object to initialize
 * @param bigrng Parent range object 
 * @param sublower Lower indices of sub-range
 * @param subupper Upper indices of sub-range
 */
void gkyl_sub_range_init(struct gkyl_range *rng,
  const struct gkyl_range *bigrng, const int *sublower, const int *subupper);

/**
 * Creates a new range that is a split of the given range. The only
 * place split matters is for iterators. Iterators for split-ranges
 * only walk over the set of indices owned by that split.
 *
 * @param rng Range object to split
 * @param nsplits Number of splits
 * @param tid Split ID [0, nsplits)
 * @return Split range
 */
struct gkyl_range gkyl_range_split(struct gkyl_range *rng, int nsplits, int tid);

/**
 * Return the number of elements looped over by iterator for this
 * range.
 *
 * @param rng Range object
 * @return number of elements looped over by iterator
 */
long gkyl_range_split_len(const struct gkyl_range *rng);

/**
 * Return range which has some directions removed by setting the index
 * in those directions to fixed values. The "deflated" range has lower
 * dimension than the parent 'rng' object. The indexing into the
 * returned lower dimensional range gives the same index as the
 * corresponding location in the parent range (with the missing
 * indices set to 'locDir[dir]'). 
 *
 * @param srng Deflated range.
 * @param rng Range object to deflate
 * @param remDir 'ndim' Array of flags: 0 to keep direction, 1 to remove
 * @param loc Index to set removed direction.
 */
void gkyl_range_deflate(struct gkyl_range* srng,
  const struct gkyl_range* rng, const int *remDir, const int *locDir);

/**
 * Return range which has 'dir' direction shortened to length
 * 'len', reducing the upper limit of 'range' in that direction.
 * The shortened range has the same dimensions and the same
 * start index in 'dir'.
 *
 * @param rng Shortened range.
 * @param range Range object to shorten
 * @param dir Direction to shorten
 * @param len Length of shortened direction
 */
void gkyl_range_shorten_from_above(struct gkyl_range *rng,
  const struct gkyl_range* range, int dir, int len);

/**
 * Return range which has 'dir' direction shortened to length
 * 'len', increasing the lower limit of 'range' in that direction.
 * The shortened range has the same dimensions and the same
 * start index in 'dir'.
 *
 * @param rng Shortened range.
 * @param range Range object to shorten
 * @param dir Direction to shorten
 * @param len Length of shortened direction
 */
void gkyl_range_shorten_from_below(struct gkyl_range *rng,
  const struct gkyl_range* range, int dir, int len);

/**
 * Return a new range that is an extension of the input range. The
 * lower index in dir is reduced by elo[dir] and upper index increased
 * by eup[dir].
 *
 * @param erng Extended range
 * @param rng Range to extend
 * @param elo Lower in dir is reduced by elo[dir]
 * @param eup Upper in dir is increased by eup[dir]
 */
void gkyl_range_extend(struct gkyl_range *erng, const struct gkyl_range *rng,
  const int *elo, const int *eup);

/**
 * Return a new range that is an extension of the input range. The
 * lower index in dir is reduced by elo[dir] and upper index increased
 * by eup[dir]. This method only extends the range in the directions
 * other than the input @a dir.
 *
 * @param erng Extended range
 * @param dir Direction to skip extension
 * @param rng Range to extend
 * @param elo Lower in dir is reduced by elo[dir]
 * @param eup Upper in dir is increased by eup[dir]
 */
void gkyl_range_perp_extend(struct gkyl_range *erng, int dir,
  const struct gkyl_range* rng, const int *elo, const int *eup);

/**
 * Return range in direction 'dir' which corresponds to the "lower
 * skin" cells.  Lower skin cells refer to the second inner-most layer
 * of cells on the lower end of the range.
 *
 * @param srng Skin range
 * @param range Range object to find lower skin cells
 * @param dir Direction to find lower skin cells in
 * @param nskin Number of skin cells
 */
void gkyl_range_lower_skin(struct gkyl_range* srng,
  const struct gkyl_range* range, int dir, int nskin);

/**
 * Return range in direction 'dir' which corresponds to the "upper
 * skin" cells.  Upper skin cells refer to the second inner-most layer
 * of cells on the upper end of the range.
 *
 * @param srng Skin range
 * @param range Range object to find upper skin cells
 * @param dir Direction to find upper skin cells in
 * @param nskin Number of skin cells
 */
void gkyl_range_upper_skin(struct gkyl_range* srng,
  const struct gkyl_range* range, int dir, int nskin);

/**
 * Create ghost and skin sub-ranges given parent *extended* range. The
 * skin and ghost ranges are sub-ranges of the parent range and DO NOT
 * include corners. For 2D, dir=1 and nghost = { 1, 1} skin and ghost
 * are the cells marked below ("S"kin, "G"ghost)
 *
 * Lower-edge:
 * +--+--+--+
 * |  |  |  |
 * +--+--+--+
 * |G |S |  |
 * +--+--+--+
 * |G |S |  |
 * +--+--+--+
 * |G |S |  |
 * +--+--+--+
 * |  |  |  |
 * +--+--+--+
 *
 * Upper-edge:
 * +--+--+--+
 * |  |  |  |
 * +--+--+--+
 * |  |S |G |
 * +--+--+--+
 * |  |S |G |
 * +--+--+--+
 * |  |S |G |
 * +--+--+--+
 * |  |  |  |
 * +--+--+--+
 *
 * @param skin On output, skin range
 * @param ghost On outout, ghost range
 * @param dir Direction in which skin/ghost are computed
 * @param edge Edge on which skin/ghost are computed
 * @param parent Range for which skin/ghost are computed
 * @param nghost Number of ghost cells in 'dir' are nghost[dir]
 */
void gkyl_skin_ghost_ranges(struct gkyl_range *skin, struct gkyl_range *ghost,
  int dir, enum gkyl_edge_loc edge, const struct gkyl_range *parent, const int *nghost);

/**
 * Create ghost and skin sub-ranges given parent *extended* range. The
 * skin and ghost ranges are sub-ranges of the parent range. The
 * ranges include the corners.  For 2D, dir=1 and nghost = { 1, 1}
 * skin and ghost are the cells marked below ("S"kin, "G"ghost)
 *
 * Lower-edge:
 * +--+--+--+
 * |G |S |  |
 * +--+--+--+
 * |G |S |  |
 * +--+--+--+
 * |G |S |  |
 * +--+--+--+
 * |G |S |  |
 * +--+--+--+
 * |G |S |  |
 * +--+--+--+
 *
 * Upper-edge:
 * +--+--+--+
 * |  |S |G |
 * +--+--+--+
 * |  |S |G |
 * +--+--+--+
 * |  |S |G |
 * +--+--+--+
 * |  |S |G |
 * +--+--+--+
 * |  |S |G |
 * +--+--+--+
 *
 * @param skin On output, skin range
 * @param ghost On outout, ghost range
 * @param dir Direction in which skin/ghost are computed
 * @param edge Edge on which skin/ghost are computed
 * @param parent Range for which skin/ghost are computed
 * @param nghost Number of ghost cells in 'dir' are nghost[dir]
 */
void gkyl_skin_ghost_with_corners_ranges(struct gkyl_range *skin, struct gkyl_range *ghost,
  int dir, enum gkyl_edge_loc edge, const struct gkyl_range *parent, const int *nghost);

/**
 * Compute intersection of two ranges. No sub-range information is
 * propagated to the new range object.
 * 
 * @param irng Intersection of r1 and r2
 * @param r1 Range to intersect
 * @param r2 Range to intersect
 * @return 1 if intersection is not-empty, 0 otherwise
 */
int gkyl_range_intersect(struct gkyl_range *irng, const struct gkyl_range *r1,
  const struct gkyl_range *r2);

/**
 * Compute intersection of two ranges. The intersection is a sub-range
 * of @a r1.
 * 
 * @param irng Intersection of r1 and r2. 
 * @param r1 Range to intersect. irng is sub-range of r1
 * @param r2 Range to intersect
 * @return 1 if intersection is not-empty, 0 otherwise
 */
int gkyl_sub_range_intersect(struct gkyl_range* irng,
  const struct gkyl_range *r1, const struct gkyl_range *r2);

/**
 * Check if range touches the lower edge of parent range in direction
 * dir.
 *
 * @param dir Direction to check
 * @param range Inner range
 * @param parent Parent range
 * @return true if range is on lower edge, false otherwise
 */
bool gkyl_range_is_on_lower_edge(int dir, const struct gkyl_range *range,
  const struct gkyl_range *parent);

/**
 * Check if range touches the upper edge of parent range in direction
 * dir.
 *
 * @param dir Direction to check
 * @param range Inner range
 * @param parent Parent range
 * @return true if range is on upper edge, false otherwise
 */
bool gkyl_range_is_on_upper_edge(int dir, const struct gkyl_range *range,
  const struct gkyl_range *parent);

/**
 * Check if @a targ range shares an edge with the @a base range. The
 * edges do not be fully shared but any edge overlap will be
 * checked.
 *
 * @param base Base range wrt which edge overlap is checked
 * @param targ Target range to check
 * @return direction and edge. Returned struct eloc is set
 *   to GKYL_NO_EDGE if ranges dont match.
 */
struct gkyl_range_dir_edge gkyl_range_edge_match(const struct gkyl_range *base,
  const struct gkyl_range *targ);
                                                     
/**
 * General indexing function. Returns linear index into the index
 * range mapped by 'range'.
 *
 * @param range Range object to index
 * @param idx Index for which to compute linear index
 */
GKYL_CU_DH
static inline long
gkyl_range_idx(const struct gkyl_range* range, const int *idx)
{
#define RI(...) gkyl_ridx(*range, __VA_ARGS__)
  switch (range->ndim) {
    case 0:
      return range->ac[0];
      break;    
    case 1:
      return RI(idx[0]); 
      break;
    case 2:
      return RI(idx[0], idx[1]);
      break;
    case 3:
      return RI(idx[0], idx[1], idx[2]);
      break;
    case 4:
      return RI(idx[0], idx[1], idx[2], idx[3]);
      break;
    case 5:
      return RI(idx[0], idx[1], idx[2], idx[3], idx[4]);
      break;
    case 6:
      return RI(idx[0], idx[1], idx[2], idx[3], idx[4], idx[5]);
      break;
    case 7:
      return RI(idx[0], idx[1], idx[2], idx[3], idx[4], idx[5], idx[6]);
      break;
  }
  return 0;
#undef RI
}

/**
 * Compute offset given relative index. So for a 2D range, idx[2] = {1,
 * 0} will compute the relative offset from {i, j} to {i+1, j} in the
 * mapping from indices to a linear integer space.
 *
 * @param range Range to find offset in.
 * @param idx Relative index for offset calculation
 * @return Relatice offset to idx.
 */
GKYL_CU_DH
static inline long
gkyl_range_offset(const struct gkyl_range* range, const int *idx)
{
  return gkyl_range_idx(range, idx) - range->linIdxZero;
}

/**
 * Inverse indexer, mapping a linear index to an N-dimension index
 * into 'range' object.
 *
 * @param range Range object to map into
 * @param loc Linear index in [0, range->volume)
 * @param idx On output, the N-dimensional index into 'range'
 */
GKYL_CU_DH
static inline void
gkyl_range_inv_idx(const struct gkyl_range *range, long loc, int *idx)
{
  long n = loc;
  for (int i=1; i<=range->ndim; ++i) {
    long quot = n/range->ac[i];
    long rem = n % range->ac[i];
    idx[i-1] = quot + range->ilo[i-1];
    n = rem;
  }
}

/**
 * Inverse indexer for use with a sub_range, mapping a linear index to
 * an N-dimension index into 'range' object.  Behavior is such that
 * loc = 0 gives idx = {0, 0, ...}.
 *
 * @param range Range object to map into
 * @param loc Linear index in [0, range->volume)
 * @param idx On output, the N-dimensional index into 'range'
 */
GKYL_CU_DH
static inline void
gkyl_sub_range_inv_idx(const struct gkyl_range *range, long loc, int *idx)
{
  long n = loc;
  for (int i=1; i<=range->ndim; ++i) {
    long quot = n/range->iac[i];
    long rem = n % range->iac[i];
    idx[i-1] = quot + range->lower[i-1];
    n = rem;
  }
}

/**
 * Create iterator. The returned iterator can be used in a 'while'
 * loop using gkyl_range_iter_next method to loop over the index set
 * spanned by the region.
 *
 * @param range Range object.
 * @return New iterator object for 'range'
 */
void gkyl_range_iter_init(struct gkyl_range_iter *iter,
  const struct gkyl_range* range);

/**
 * Create iterator, ignoring split information in range.
 *
 * @param range Range object.
 * @return New iterator object for 'range'
 */
void gkyl_range_iter_no_split_init(struct gkyl_range_iter *iter,
  const struct gkyl_range* range);

/**
 * Get next index into range. The iter->idx array holds the next
 * index. This should not be modified by the user!
 *
 * @param iter Iterator object. On exit, iter->idx has the next index
 * @return 1 if there are more indices remaining, 0 if done.
 */
int gkyl_range_iter_next(struct gkyl_range_iter *iter);

/**
 * Create skip-iterator. The returned iterator can be used in a nested
 * loop structure: a for loop inside a 'while' loop.
 *
 * @param range Range object.
 * @return New iterator object for 'range'
 */
void gkyl_range_skip_iter_init(struct gkyl_range_skip_iter *iter,
  const struct gkyl_range* range);

/**
 * Print range information to file object.
 *
 * @param range Range object to print
 * @param nm Name of range
 * @param fp File object to print range information
 */
void gkyl_print_range(const struct gkyl_range* range, const char *nm, FILE *fp);

/**
 * Compares two ranges: ranges are the same if they have the same
 * dimensions and lower and upper indices.
 *
 * @param r1 Range 1 to compare
 * @param r2 Range 2 to compare
 * @return true if ranges are same, false otherwise
 */
bool gkyl_range_compare(const struct gkyl_range* r1, const struct gkyl_range* r2);

/**
 * Get the number of dimensions from range (for opaque pointer).
 *
 * @return Number of dimensions.
 */
int gkyl_range_get_ndim(const struct gkyl_range* range);

/**
 * Get the lower extents of a range (for opaque pointer).
 *
 * @param range Range object.
 * @param lower Lower extents of the range (pre-allocated).
 */
void gkyl_range_get_lower(const struct gkyl_range* range, int *lower);

/**
 * Get the upper extents of a range (for opaque pointer).
 *
 * @param range Range object.
 * @param upper Upper extents of the range (pre-allocated).
 */
void gkyl_range_get_upper(const struct gkyl_range* range, int *upper);

/**
 * Get the volume of a range (for opaque pointer).
 *
 * @param range Range object.
 * @return Volume of the range.
 */
long gkyl_range_get_volume(const struct gkyl_range* range);

/**
 * Free dynamically allocated range object (with gkyl_range_new).
 *
 * @param range Range object.
 */
void gkyl_range_release(struct gkyl_range* range);

