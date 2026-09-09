// Additional unit tests for gkyl_range (indexing, shape, iteration)
#include <acutest.h>
#include <gkyl_range.h>

void
test_range_init_shape()
{
  int lower[] = { 1, 1 }, upper[] = { 4, 5 };
  struct gkyl_range rng;
  gkyl_range_init(&rng, 2, lower, upper);

  TEST_CHECK( rng.ndim == 2 );
  TEST_CHECK( rng.volume == 4*5 );
  TEST_CHECK( gkyl_range_shape(&rng, 0) == 4 );
  TEST_CHECK( gkyl_range_shape(&rng, 1) == 5 );
  TEST_CHECK( !gkyl_range_is_sub_range(&rng) );
}

void
test_range_init_from_shape()
{
  int shape[] = { 3, 4, 5 };
  struct gkyl_range rng;
  gkyl_range_init_from_shape(&rng, 3, shape);

  TEST_CHECK( rng.ndim == 3 );
  TEST_CHECK( rng.volume == 60 );
  TEST_CHECK( gkyl_range_shape(&rng, 0) == 3 );
  TEST_CHECK( gkyl_range_shape(&rng, 1) == 4 );
  TEST_CHECK( gkyl_range_shape(&rng, 2) == 5 );
  // default lower index is 0
  TEST_CHECK( rng.lower[0] == 0 );
}

void
test_range_index_roundtrip()
{
  int lower[] = { 1, 1 }, upper[] = { 3, 4 };
  struct gkyl_range rng;
  gkyl_range_init(&rng, 2, lower, upper);

  // Every linear index in [0, volume) maps back to a valid multi-index
  // whose forward index recovers the same linear index.
  for (long loc=0; loc<rng.volume; ++loc) {
    int idx[2];
    gkyl_range_inv_idx(&rng, loc, idx);
    TEST_CHECK( idx[0] >= rng.lower[0] && idx[0] <= rng.upper[0] );
    TEST_CHECK( idx[1] >= rng.lower[1] && idx[1] <= rng.upper[1] );
    long lin = gkyl_range_idx(&rng, idx);
    TEST_CHECK( lin == loc );
  }
}

void
test_range_index_unique()
{
  int lower[] = { 0, 0 }, upper[] = { 2, 2 };
  struct gkyl_range rng;
  gkyl_range_init(&rng, 2, lower, upper);

  // All linear indices over the range are distinct and cover [0, volume).
  int seen[9] = { 0 };
  for (int i=lower[0]; i<=upper[0]; ++i) {
    for (int j=lower[1]; j<=upper[1]; ++j) {
      int idx[] = { i, j };
      long lin = gkyl_range_idx(&rng, idx);
      TEST_CHECK( lin >= 0 && lin < rng.volume );
      seen[lin]++;
    }
  }
  for (int k=0; k<9; ++k)
    TEST_CHECK( seen[k] == 1 );
}

void
test_range_iter()
{
  int lower[] = { 1, 1 }, upper[] = { 3, 3 };
  struct gkyl_range rng;
  gkyl_range_init(&rng, 2, lower, upper);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &rng);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    TEST_CHECK( iter.idx[0] >= 1 && iter.idx[0] <= 3 );
    TEST_CHECK( iter.idx[1] >= 1 && iter.idx[1] <= 3 );
    count++;
  }
  TEST_CHECK( count == rng.volume );
  TEST_CHECK( count == 9 );
}

void
test_range_shorten()
{
  int lower[] = { 1, 1 }, upper[] = { 4, 6 };
  struct gkyl_range rng;
  gkyl_range_init(&rng, 2, lower, upper);

  struct gkyl_range srng;
  gkyl_range_shorten_from_above(&srng, &rng, 1, 1);
  // Shortened to 1 cell in direction 1.
  TEST_CHECK( gkyl_range_shape(&srng, 1) == 1 );
  TEST_CHECK( gkyl_range_shape(&srng, 0) == 4 );
  TEST_CHECK( srng.volume == 4 );
}

void
test_sub_range()
{
  int lower[] = { 1, 1 }, upper[] = { 10, 10 };
  struct gkyl_range rng;
  gkyl_range_init(&rng, 2, lower, upper);

  int sublo[] = { 2, 2 }, subup[] = { 5, 5 };
  struct gkyl_range sub;
  gkyl_sub_range_init(&sub, &rng, sublo, subup);

  TEST_CHECK( gkyl_range_is_sub_range(&sub) );
  TEST_CHECK( sub.volume == 16 );
  TEST_CHECK( gkyl_range_shape(&sub, 0) == 4 );
  TEST_CHECK( gkyl_range_shape(&sub, 1) == 4 );
}

TEST_LIST = {
  { "range_init_shape", test_range_init_shape },
  { "range_init_from_shape", test_range_init_from_shape },
  { "range_index_roundtrip", test_range_index_roundtrip },
  { "range_index_unique", test_range_index_unique },
  { "range_iter", test_range_iter },
  { "range_shorten", test_range_shorten },
  { "sub_range", test_sub_range },
  { NULL, NULL },
};
