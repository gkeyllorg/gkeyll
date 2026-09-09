// Additional unit tests for gkyl_array element-wise operations (host)
#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_util.h>

void
test_array_new_meta()
{
  struct gkyl_array *a = gkyl_array_new(GKYL_DOUBLE, 3, 10);
  TEST_CHECK( a->type == GKYL_DOUBLE );
  TEST_CHECK( a->ncomp == 3 );
  TEST_CHECK( a->size == 10 );
  TEST_CHECK( gkyl_array_is_cu_dev(a) == false );
  gkyl_array_release(a);
}

void
test_array_clear()
{
  struct gkyl_array *a = gkyl_array_new(GKYL_DOUBLE, 1, 50);
  gkyl_array_clear(a, 3.5);
  double *d = a->data;
  for (unsigned i=0; i<a->size*a->ncomp; ++i)
    TEST_CHECK( d[i] == 3.5 );
  gkyl_array_release(a);
}

void
test_array_scale()
{
  struct gkyl_array *a = gkyl_array_new(GKYL_DOUBLE, 1, 20);
  gkyl_array_clear(a, 2.0);
  gkyl_array_scale(a, 4.0);
  double *d = a->data;
  for (unsigned i=0; i<a->size; ++i)
    TEST_CHECK( gkyl_compare_double(d[i], 8.0, 1e-14) );
  gkyl_array_release(a);
}

void
test_array_accumulate()
{
  // out = out + a*inp
  struct gkyl_array *out = gkyl_array_new(GKYL_DOUBLE, 1, 20);
  struct gkyl_array *inp = gkyl_array_new(GKYL_DOUBLE, 1, 20);
  gkyl_array_clear(out, 1.0);
  gkyl_array_clear(inp, 2.0);
  gkyl_array_accumulate(out, 3.0, inp);
  double *d = out->data;
  for (unsigned i=0; i<out->size; ++i)
    TEST_CHECK( gkyl_compare_double(d[i], 1.0 + 3.0*2.0, 1e-14) );
  gkyl_array_release(out); gkyl_array_release(inp);
}

void
test_array_set()
{
  // out = a*inp
  struct gkyl_array *out = gkyl_array_new(GKYL_DOUBLE, 1, 20);
  struct gkyl_array *inp = gkyl_array_new(GKYL_DOUBLE, 1, 20);
  gkyl_array_clear(out, 99.0);
  gkyl_array_clear(inp, 5.0);
  gkyl_array_set(out, 2.0, inp);
  double *d = out->data;
  for (unsigned i=0; i<out->size; ++i)
    TEST_CHECK( gkyl_compare_double(d[i], 10.0, 1e-14) );
  gkyl_array_release(out); gkyl_array_release(inp);
}

void
test_array_shiftc()
{
  // shift component k by a (out[k] += a)
  struct gkyl_array *a = gkyl_array_new(GKYL_DOUBLE, 3, 10);
  gkyl_array_clear(a, 0.0);
  gkyl_array_shiftc(a, 7.0, 1); // shift component 1
  for (unsigned i=0; i<a->size; ++i) {
    double *row = gkyl_array_fetch(a, i);
    TEST_CHECK( row[0] == 0.0 );
    TEST_CHECK( gkyl_compare_double(row[1], 7.0, 1e-14) );
    TEST_CHECK( row[2] == 0.0 );
  }
  gkyl_array_release(a);
}

void
test_array_copy()
{
  struct gkyl_array *a = gkyl_array_new(GKYL_DOUBLE, 2, 15);
  struct gkyl_array *b = gkyl_array_new(GKYL_DOUBLE, 2, 15);
  double *ad = a->data;
  for (unsigned i=0; i<a->size*a->ncomp; ++i) ad[i] = 0.5*i;
  gkyl_array_copy(b, a);
  double *bd = b->data;
  for (unsigned i=0; i<b->size*b->ncomp; ++i)
    TEST_CHECK( bd[i] == 0.5*i );
  gkyl_array_release(a); gkyl_array_release(b);
}

void
test_array_fetch_multicomp()
{
  // Write per-component values via fetch and read back.
  struct gkyl_array *a = gkyl_array_new(GKYL_DOUBLE, 4, 8);
  for (unsigned i=0; i<a->size; ++i) {
    double *row = gkyl_array_fetch(a, i);
    for (unsigned c=0; c<a->ncomp; ++c)
      row[c] = 100.0*i + c;
  }
  for (unsigned i=0; i<a->size; ++i) {
    const double *row = gkyl_array_cfetch(a, i);
    for (unsigned c=0; c<a->ncomp; ++c)
      TEST_CHECK( row[c] == 100.0*i + c );
  }
  gkyl_array_release(a);
}

TEST_LIST = {
  { "array_new_meta", test_array_new_meta },
  { "array_clear", test_array_clear },
  { "array_scale", test_array_scale },
  { "array_accumulate", test_array_accumulate },
  { "array_set", test_array_set },
  { "array_shiftc", test_array_shiftc },
  { "array_copy", test_array_copy },
  { "array_fetch_multicomp", test_array_fetch_multicomp },
  { NULL, NULL },
};
