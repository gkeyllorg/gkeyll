#include <acutest.h>
#include <gkyl_ref_count.h>
#include <gkyl_alloc.h>


// global to indicate if free was called
static int free_called = 0;

struct range {
  int value;
  struct gkyl_ref_count ref_count;
};

void
range_free(const struct gkyl_ref_count* rc)
{
  struct range *on_dev = container_of(rc, struct range, ref_count);
  free_called = 1;
  gkyl_free(on_dev);
}

struct range*
range_new(int value)
{
  struct range *rng = gkyl_malloc(sizeof(*rng));
  rng->value = value;
  rng->ref_count = (struct gkyl_ref_count) { range_free, 1 };
  return rng;
}

struct range*
range_acquire(const struct range *rng)
{
  gkyl_ref_count_inc(&rng->ref_count);
  return (struct range*) rng;
}

void
range_release(const struct range *rng)
{
  gkyl_ref_count_dec(&rng->ref_count);
}

void
test_ref_count()
{
  struct range *rng = range_new(10);
  TEST_CHECK( rng->ref_count.count == 1 );

  struct range *rngp = range_acquire(rng);
  TEST_CHECK( rng->ref_count.count == 2 );

  range_release(rngp);
  TEST_CHECK( rng->ref_count.count == 1 );

  range_release(rngp);
  TEST_CHECK( free_called == 1 );
}

TEST_LIST = {
  { "ref_count", test_ref_count },
  { NULL, NULL },  
};
