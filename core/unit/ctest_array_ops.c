#include <acutest.h>
#include <mpack.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_rio_format_desc.h>
#include <gkyl_array_rio_priv.h>
#include <gkyl_elem_type_priv.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

static struct gkyl_array*
mkarr(bool use_gpu, long nc, long size)
{
  // Allocate array (filled with zeros)
  struct gkyl_array* a = use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
                                : gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

void test_array_clear()
{
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 1, 10);

  gkyl_array_clear(a1, 0.5);
  double *a1_d  = a1->data;  

  for (unsigned i=0; i<a1->size; ++i)
    TEST_CHECK( gkyl_compare(a1_d[i], 0.5, 1e-14) );

  gkyl_array_release(a1);
}

void test_array_clear_range()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 1, range.volume);
  gkyl_array_clear_range(a1, 0.5, &range);

  double *a1_d = a1->data;
  for (unsigned i=0; i<a1->size; ++i)
    TEST_CHECK( gkyl_compare(a1_d[i], 0.5, 1e-14) );

  gkyl_array_release(a1);
}

void test_array_accumulate()
{
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 1, 10);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 1, 10);

  double *a1_d  = a1->data, *a2_d = a2->data;
  for (unsigned i=0; i<a1->size; ++i) {
    a1_d[i] = i*1.0;
    a2_d[i] = i*0.1;
  }

  gkyl_array_accumulate(a1, 0.5, a2);

  for (unsigned i=0; i<a1->size; ++i)
    TEST_CHECK( gkyl_compare(a1_d[i], i*1.0+0.5*i*0.1, 1e-14) );

  gkyl_array_release(a1);
  gkyl_array_release(a2);
}

void test_array_accumulate_range()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 8, range.volume);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 3, range.volume);

  // test a1 = a1 + 0.5*a2
  gkyl_array_clear(a1, 0.5);
  gkyl_array_clear(a2, 1.5);

  gkyl_array_accumulate_range(a1, 0.5, a2, &range);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *a1d = gkyl_array_fetch(a1, loc);
    for (int i=0; i<3; ++i)
      TEST_CHECK( a1d[i] == 0.5 + 0.5*1.5 );
    for (int i=3; i<8; ++i)
      TEST_CHECK( a1d[i] == 0.5);
  }

  // test a2 = a2 + 0.5*a
  gkyl_array_clear(a1, 0.5);
  gkyl_array_clear(a2, 1.5);

  gkyl_array_accumulate_range(a2, 0.5, a1, &range);

  gkyl_range_iter_init(&iter, &range);
  
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *a2d = gkyl_array_fetch(a2, loc);
    for (int i=0; i<3; ++i)
      TEST_CHECK( a2d[i] == 1.5 + 0.5*0.5 );
  }  

  gkyl_array_release(a1);
  gkyl_array_release(a2);
}

void test_array_accumulate_offset()
{
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 2, 10);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 3*a1->ncomp, 10);

  double *a1_d  = a1->data, *a2_d = a2->data;

  // test a1 = 0.1*a2[a1->ncomp]
  for (unsigned i=0; i<a1->size; ++i)
    for (unsigned j=0; j<a1->ncomp; ++j)
      a1_d[i*a1->ncomp+j] = i*1.0+j;

  for (unsigned i=0; i<a2->size; ++i)
    for (unsigned j=0; j<a2->ncomp/a1->ncomp; ++j)
      for (unsigned k=0; k<a1->ncomp; ++k)
        a2_d[i*a2->ncomp+j*a1->ncomp+k] = i*0.1+k;

  gkyl_array_accumulate_offset(a1, 0.5, a2, 1*a1->ncomp);

  for (unsigned i=0; i<a1->size; ++i)
    for (unsigned j=0; j<a1->ncomp; ++j)
      TEST_CHECK( gkyl_compare(a1_d[i*a1->ncomp+j], i*1.0+j+0.5*(i*0.1+j), 1e-14) );

  // test a2[a1->ncomp] = 0.1*a1
  for (unsigned i=0; i<a1->size; ++i)
    for (unsigned j=0; j<a1->ncomp; ++j)
      a1_d[i*a1->ncomp+j] = i*1.0+j;

  gkyl_array_accumulate_offset(a2, 0.5, a1, 1*a1->ncomp);

  for (unsigned i=0; i<a1->size; ++i) {
    for (unsigned j=0; j<a1->ncomp; ++j) {
      TEST_CHECK( gkyl_compare(a2_d[i*a2->ncomp+0*a1->ncomp+j], i*0.1+j, 1e-14) );
      TEST_CHECK( gkyl_compare(a2_d[i*a2->ncomp+1*a1->ncomp+j], i*0.1+j+0.5*(i*1.0+j), 1e-14) );
      TEST_CHECK( gkyl_compare(a2_d[i*a2->ncomp+2*a1->ncomp+j], i*0.1+j, 1e-14) );
    }
  }

  gkyl_array_release(a1);
  gkyl_array_release(a2);
}

void test_array_accumulate_offset_range()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 2, range.volume);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 3*a1->ncomp, range.volume);

  // test a1 = a1+0.5*a2[a1->ncomp]
  gkyl_array_clear(a1, 0.5);
  gkyl_array_clear(a2, 1.5);

  gkyl_array_accumulate_offset_range(a1, 0.5, a2, 1*a1->ncomp, &range);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *a1d = gkyl_array_fetch(a1, loc);
    for (int i=0; i<a1->ncomp; ++i)
      TEST_CHECK( a1d[i] == 0.5+0.5*1.5 );
  }

  // test a2[a1->ncomp] = a2[a1->ncomp]+0.5*a1
  gkyl_array_clear(a1, 0.5);
  gkyl_array_clear(a2, 1.5);

  gkyl_array_accumulate_offset_range(a2, 0.5, a1, 1*a1->ncomp, &range);

  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *a2d = gkyl_array_fetch(a2, loc);
    for (int i=0; i<a1->ncomp; ++i) {
      TEST_CHECK( a2d[i+0*a1->ncomp] == 1.5 );
      TEST_CHECK( a2d[i+1*a1->ncomp] == 1.5+0.5*0.5 );
      TEST_CHECK( a2d[i+2*a1->ncomp] == 1.5 );
    }
  }  

  gkyl_array_release(a1);
  gkyl_array_release(a2);
}

void test_array_combine()
{
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 1, 10);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 1, 10);
  struct gkyl_array *b = gkyl_array_new(GKYL_DOUBLE, 1, 10);

  double *a1_d  = a1->data, *a2_d = a2->data, *b_d = b->data;
  for (unsigned i=0; i<a1->size; ++i) {
    a1_d[i] = i*1.0;
    a2_d[i] = i*0.1;
    b_d[i] = 10.5;
  }

  // b = 0.5*a1 + 2.5*a2
  gkyl_array_accumulate(gkyl_array_set(b, 0.5, a1), 2.5, a2);

  for (unsigned i=0; i<a1->size; ++i)
    TEST_CHECK( gkyl_compare(b_d[i], 0.5*i*1.0+2.5*i*0.1, 1e-14) );

  gkyl_array_release(a1);
  gkyl_array_release(a2);
  gkyl_array_release(b);
}

void test_array_set()
{
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 1, 10);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 1, 10);

  double *a1_d  = a1->data, *a2_d = a2->data;
  for (unsigned i=0; i<a1->size; ++i) {
    a1_d[i] = i*1.0;
    a2_d[i] = i*0.1;
  }

  gkyl_array_set(a1, 0.5, a2);

  for (unsigned i=0; i<a1->size; ++i)
    TEST_CHECK( gkyl_compare(a1_d[i], 0.5*i*0.1, 1e-14) );

  gkyl_array_release(a1);
  gkyl_array_release(a2);
}

void test_array_set_range()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 8, range.volume);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 3, range.volume);

  // test a1 = 0.5*a2
  gkyl_array_clear_range(a1, 0.5, &range);
  gkyl_array_clear_range(a2, 1.5, &range);

  gkyl_array_set_range(a1, 0.5, a2, &range);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *a1d = gkyl_array_fetch(a1, loc);
    for (int i=0; i<3; ++i)
      TEST_CHECK( a1d[i] == 0.5*1.5 );
    for (int i=3; i<8; ++i)
      TEST_CHECK( a1d[i] == 0.5);
  }

  // test a2 = 0.5*a1
  gkyl_array_clear(a1, 0.5);
  gkyl_array_clear(a2, 1.5);

  gkyl_array_set_range(a2, 0.5, a1, &range);

  gkyl_range_iter_init(&iter, &range);
  
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *a2d = gkyl_array_fetch(a2, loc);
    for (int i=0; i<3; ++i)
      TEST_CHECK( a2d[i] == 0.5*0.5 );
  }  

  gkyl_array_release(a1);
  gkyl_array_release(a2);
}

void test_array_set_offset()
{
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 2, 10);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 3*a1->ncomp, a1->size);

  double *a1_d = a1->data, *a2_d = a2->data;

  // Assign a component of the vector to the scalar.
  for (unsigned i=0; i<a1->size; ++i)
    for (unsigned j=0; j<a1->ncomp; ++j)
      a1_d[i*a1->ncomp+j] = i*0.2+2*j;

  for (unsigned i=0; i<a2->size; ++i)
    for (unsigned j=0; j<a2->ncomp/a1->ncomp; ++j)
      for (unsigned k=0; k<a1->ncomp; ++k)
        a2_d[i*a2->ncomp+j*a1->ncomp+k] = i*0.1+k;

  gkyl_array_set_offset(a1, 0.5, a2, 1*a1->ncomp);

  for (unsigned i=0; i<a1->size; ++i)
    for (unsigned j=0; j<a1->ncomp; ++j)
      TEST_CHECK( gkyl_compare(a1_d[i*a1->ncomp+j], 0.5*(i*0.1+j), 1e-14) );

  // Assign the scalar to a component of the vector.
  for (unsigned i=0; i<a1->size; ++i)
    for (unsigned j=0; j<a1->ncomp; ++j)
      a1_d[i*a1->ncomp+j] = i*0.2+2*j;

  gkyl_array_set_offset(a2, 2., a1, 1*a1->ncomp);

  for (unsigned i=0; i<a1->size; ++i) {
    for (unsigned j=0; j<a1->ncomp; ++j) {
      TEST_CHECK( gkyl_compare(a2_d[i*a2->ncomp+0*a1->ncomp+j], i*0.1+j, 1e-14) );
      TEST_CHECK( gkyl_compare(a2_d[i*a2->ncomp+1*a1->ncomp+j], 2.*(i*0.2+2*j), 1e-14) );
      TEST_CHECK( gkyl_compare(a2_d[i*a2->ncomp+2*a1->ncomp+j], i*0.1+j, 1e-14) );
    }
  }

  gkyl_array_release(a1);
  gkyl_array_release(a2);
}

void test_array_set_offset_range()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 2, range.volume);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 3*a1->ncomp, range.volume);

  // test a1 = 0.1*a2[a1->ncomp]
  gkyl_array_clear_range(a1, 0.5, &range);
  gkyl_array_clear_range(a2, 1.5, &range);

  gkyl_array_set_offset_range(a1, 0.1, a2, 1*a1->ncomp, &range);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *a1d = gkyl_array_fetch(a1, loc);
    for (int i=0; i<a1->ncomp; ++i)
      TEST_CHECK( a1d[i] == 0.1*1.5 );
  }

  // test a2[a1->ncomp] = 0.1*a1
  gkyl_array_clear(a1, 0.5);
  gkyl_array_clear(a2, 1.5);

  gkyl_array_set_offset_range(a2, 0.1, a1, 1*a1->ncomp, &range);

  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *a2d = gkyl_array_fetch(a2, loc);
    for (int i=0; i<a1->ncomp; ++i) {
      TEST_CHECK( a2d[i+0*a1->ncomp] == 1.5 );
      TEST_CHECK( a2d[i+1*a1->ncomp] == 0.1*0.5 );
      TEST_CHECK( a2d[i+2*a1->ncomp] == 1.5 );
    }
  }  

  gkyl_array_release(a1);
  gkyl_array_release(a2);
}

void test_array_scale()
{
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 1, 10);

  double *a1_d  = a1->data;
  for (unsigned i=0; i<a1->size; ++i) {
    a1_d[i] = i*1.0;
  }

  gkyl_array_scale(a1, 0.25);

  for (unsigned i=0; i<a1->size; ++i)
    TEST_CHECK( gkyl_compare(a1_d[i], i*0.25, 1e-14) );

  gkyl_array_release(a1);
}

void test_array_scale_by_cell()
{
  struct gkyl_array *a1 = mkarr(false, 3, 10);
  struct gkyl_array *s = mkarr(false, 1, 10);

  double *a1_d  = a1->data;
  for (unsigned i=0; i<a1->size; ++i) {
    a1_d[i] = i*1.0;
  }
  double *s_d  = s->data;
  for (unsigned i=0; i<s->size; ++i) {
    s_d[i] = i*1.0;
  }

  gkyl_array_scale_by_cell(a1, s);

  for (unsigned i=0; i<a1->size; ++i) {
    int fact = (i/a1->ncomp);
    TEST_CHECK( gkyl_compare(a1_d[i], i*1.0*fact, 1e-14) );
  }

  gkyl_array_release(a1);
  gkyl_array_release(s);
}

void test_array_divide_by_cell()
{
  struct gkyl_array *a1 = mkarr(false, 3, 10);
  struct gkyl_array *s = mkarr(false, 1, 10);

  double *a1_d  = a1->data;
  for (unsigned i=0; i<a1->size; ++i) {
    a1_d[i] = i*1.0;
  }
  double *s_d  = s->data;
  for (unsigned i=0; i<s->size; ++i) {
    s_d[i] = i*1.0;
  }

  gkyl_array_divide_by_cell(a1, s);

  for (unsigned i=0; i<a1->size; ++i) {
    int fact = (i/a1->ncomp);
    TEST_CHECK( gkyl_compare(a1_d[i], i*1.0/fact, 1e-14) );
  }

  gkyl_array_release(a1);
  gkyl_array_release(s);
}

void test_array_shiftc()
{
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 3, 10);

  double s = -0.5;
  double *a1_d = a1->data;
  for (unsigned i=0; i<a1->size; ++i) {
    for (size_t k=0; k<a1->ncomp; ++k) a1_d[i*a1->ncomp+k] = i*2.0+k;
  }

  gkyl_array_shiftc(a1, s, 0);

  TEST_CHECK( gkyl_compare(a1_d[0], 0*1.0+0+s, 1e-14) );
  for (unsigned i=0; i<a1->size; ++i) {
    for (size_t k=1; k<a1->ncomp; ++k) 
      TEST_CHECK( gkyl_compare(a1_d[i*a1->ncomp+k], i*2.0+k, 1e-14) );
  }

  gkyl_array_release(a1);

  // Repeat the test but shifting another coefficient as well.
  int shiftks[] = {0, 2};
  int nks = sizeof(shiftks)/sizeof(shiftks[0]);

  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 4, 8);
  double *a2_d = a2->data;
  for (unsigned i=0; i<a2->size; ++i) {
    for (size_t k=0; k<a2->ncomp; ++k) a2_d[i*a2->ncomp+k] = i*2.0+k;
  }

  for (size_t l=0; l<nks; l++)
    gkyl_array_shiftc(a2, s, shiftks[l]);

  for (unsigned i=0; i<a2->size; ++i) {
    for (size_t k=0; k<a2->ncomp; ++k) {
      bool isshifted = false;
      for (size_t l=0; l<nks; l++) {
        if (shiftks[l]==k) {isshifted = true; break;}
      }
      if (isshifted)
        TEST_CHECK( gkyl_compare(a2_d[i*a2->ncomp+k], i*2.0+k+s, 1e-14) );
      else
        TEST_CHECK( gkyl_compare(a2_d[i*a2->ncomp+k], i*2.0+k, 1e-14) );
    }
  }

  gkyl_array_release(a2);
}

void test_array_invert_by_cell()
{
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 3, 8);
  double *a1_d = a1->data;

  // Set values: 1.0, 2.0, 4.0, 5.0, 10.0, 0.5, 0.25, 0.1
  double test_vals[] = {1.0, 2.0, 4.0, 5.0, 10.0, 0.5, 0.25, 0.1};
  double expected_inv[] = {1.0, 0.5, 0.25, 0.2, 0.1, 2.0, 4.0, 10.0};
  
  for (unsigned i=0; i<a1->size; ++i)
    for (size_t k=0; k<a1->ncomp; ++k)
      a1_d[i*a1->ncomp+k] = test_vals[i];

  gkyl_array_invert_by_cell(a1, a1);

  // Check inverted values
  for (unsigned i=0; i<a1->size; ++i)
    for (size_t k=0; k<a1->ncomp; ++k)
      TEST_CHECK( gkyl_compare(a1_d[i*a1->ncomp+k], expected_inv[i], 1e-14) );

  gkyl_array_release(a1);
}

void test_array_shiftc_range(bool on_gpu)
{
  int lower[] = {1}, upper[] = {10};
  struct gkyl_range range;
  gkyl_range_init(&range, 1, lower, upper);

  struct gkyl_array *a1_ho = gkyl_array_new(GKYL_DOUBLE, 3, range.volume);
  struct gkyl_array *a1 = a1_ho;
  if (on_gpu) a1 = gkyl_array_cu_dev_new(GKYL_DOUBLE, a1_ho->ncomp, a1_ho->size);

  double s = -0.5;
  double *a1_ho_d = a1_ho->data;
  for (unsigned i=0; i<a1_ho->size; ++i)
    for (size_t k=0; k<a1_ho->ncomp; ++k) a1_ho_d[i*a1_ho->ncomp+k] = i*2.0+k;

  if (on_gpu) gkyl_array_copy(a1, a1_ho);

  int lowerSub[] = {2}, upperSub[] = {6};
  struct gkyl_range subrange;
  gkyl_sub_range_init(&subrange, &range, lowerSub, upperSub);

  gkyl_array_shiftc_range(a1, s, 0, &subrange);

  if (on_gpu) gkyl_array_copy(a1_ho, a1);

  for (size_t k=1; k<a1_ho->ncomp; ++k) {
    for (unsigned i=0; i<1; ++i)
      TEST_CHECK( gkyl_compare(a1_ho_d[i*a1_ho->ncomp+k], i*2.0+k, 1e-14) );
    for (unsigned i=6; i<10; ++i)
      TEST_CHECK( gkyl_compare(a1_ho_d[i*a1_ho->ncomp+k], i*2.0+k, 1e-14) );
  }
  for (unsigned i=1; i<6; ++i)
    TEST_CHECK( gkyl_compare(a1_ho_d[i*a1_ho->ncomp+0], i*2.0+0+s, 1e-14) );

  gkyl_array_release(a1_ho);
  if (on_gpu) gkyl_array_release(a1);

  // Repeat the test but shifting another coefficient as well.
  int shiftks[] = {0, 2};
  int nks = sizeof(shiftks)/sizeof(shiftks[0]);

  lower[0] = 1;  upper[0] = 8;
  struct gkyl_range range2;
  gkyl_range_init(&range2, 1, lower, upper);

  struct gkyl_array *a2_ho = gkyl_array_new(GKYL_DOUBLE, 4, range2.volume);
  struct gkyl_array *a2 = a2_ho;
  if (on_gpu) a2 = gkyl_array_cu_dev_new(GKYL_DOUBLE, a2_ho->ncomp, a2_ho->size);

  double *a2_ho_d = a2_ho->data;
  for (unsigned i=0; i<a2_ho->size; ++i) {
    for (size_t k=0; k<a2_ho->ncomp; ++k) a2_ho_d[i*a2_ho->ncomp+k] = i*2.0+k;
  }

  if (on_gpu) gkyl_array_copy(a2, a2_ho);

  for (size_t l=0; l<nks; l++)
    gkyl_array_shiftc_range(a2, s, shiftks[l], &subrange);

  if (on_gpu) gkyl_array_copy(a2_ho, a2);

  for (size_t k=0; k<a2_ho->ncomp; ++k) {
    for (unsigned i=0; i<1; ++i)
      TEST_CHECK( gkyl_compare(a2_ho_d[i*a2_ho->ncomp+k], i*2.0+k, 1e-14) );
    for (unsigned i=6; i<8; ++i)
      TEST_CHECK( gkyl_compare(a2_ho_d[i*a2_ho->ncomp+k], i*2.0+k, 1e-14) );
  }
  for (unsigned i=1; i<6; ++i) {
    for (size_t k=0; k<a2_ho->ncomp; ++k) {
      bool isshifted = false;
      for (size_t l=0; l<nks; l++) {
        if (shiftks[l]==k) {isshifted = true; break;}
      }
      if (isshifted)
        TEST_CHECK( gkyl_compare(a2_ho_d[i*a2_ho->ncomp+k], i*2.0+k+s, 1e-14) );
      else
        TEST_CHECK( gkyl_compare(a2_ho_d[i*a2_ho->ncomp+k], i*2.0+k, 1e-14) );
    }
  }

  gkyl_array_release(a2);
  if (on_gpu) gkyl_array_release(a2);
}

void test_array_min_by_cell(bool on_gpu)
{
  struct gkyl_array *a1_ho = gkyl_array_new(GKYL_DOUBLE, 3, 10);
  struct gkyl_array *a1 = a1_ho;
  if (on_gpu) a1 = gkyl_array_cu_dev_new(GKYL_DOUBLE, a1_ho->ncomp, a1_ho->size);

  double *a1_ho_d = a1_ho->data;
  for (unsigned i=0; i<a1_ho->size; ++i) {
    for (size_t k=0; k<a1_ho->ncomp; ++k) 
      a1_ho_d[i*a1_ho->ncomp+k] = i*10.0+k;
  }

  if (on_gpu) gkyl_array_copy(a1, a1_ho);

  // Apply min with threshold 15.0
  gkyl_array_min_by_cell(a1, a1, 15.0);

  if (on_gpu) gkyl_array_copy(a1_ho, a1);

  // Check that values > 15.0 are clamped to 15.0
  for (unsigned i=0; i<a1_ho->size; ++i) {
    for (size_t k=0; k<a1_ho->ncomp; ++k) {
      double expected = (i*10.0+k < 15.0) ? i*10.0+k : 15.0;
      TEST_CHECK( gkyl_compare(a1_ho_d[i*a1_ho->ncomp+k], expected, 1e-14) );
    }
  }

  gkyl_array_release(a1_ho);
  if (on_gpu) gkyl_array_release(a1);

  // Test with negative threshold
  struct gkyl_array *a2_ho = gkyl_array_new(GKYL_DOUBLE, 4, 8);
  struct gkyl_array *a2 = a2_ho;
  if (on_gpu) a2 = gkyl_array_cu_dev_new(GKYL_DOUBLE, a2_ho->ncomp, a2_ho->size);

  double *a2_ho_d = a2_ho->data;
  for (unsigned i=0; i<a2_ho->size; ++i) {
    for (size_t k=0; k<a2_ho->ncomp; ++k) 
      a2_ho_d[i*a2_ho->ncomp+k] = (double)i - 5.0 + k*0.1;
  }

  if (on_gpu) gkyl_array_copy(a2, a2_ho);

  gkyl_array_min_by_cell(a2, a2, -2.5);

  if (on_gpu) gkyl_array_copy(a2_ho, a2);

  for (unsigned i=0; i<a2_ho->size; ++i) {
    for (size_t k=0; k<a2_ho->ncomp; ++k) {
      double orig = (double)i - 5.0 + k*0.1;
      double expected = (orig < -2.5) ? orig : -2.5;
      TEST_CHECK( gkyl_compare(a2_ho_d[i*a2_ho->ncomp+k], expected, 1e-14) );
    }
  }

  gkyl_array_release(a2_ho);
  if (on_gpu) gkyl_array_release(a2);
}

void test_array_min_range(bool on_gpu)
{
  int lower[] = {1}, upper[] = {10};
  struct gkyl_range range;
  gkyl_range_init(&range, 1, lower, upper);

  struct gkyl_array *a1_ho = gkyl_array_new(GKYL_DOUBLE, 3, range.volume);
  struct gkyl_array *a1 = a1_ho;
  if (on_gpu) a1 = gkyl_array_cu_dev_new(GKYL_DOUBLE, a1_ho->ncomp, a1_ho->size);

  double *a1_ho_d = a1_ho->data;
  for (unsigned i=0; i<a1_ho->size; ++i)
    for (size_t k=0; k<a1_ho->ncomp; ++k) 
      a1_ho_d[i*a1_ho->ncomp+k] = i*10.0+k;

  if (on_gpu) gkyl_array_copy(a1, a1_ho);

  // Apply min only to subrange [2, 6]
  int lowerSub[] = {2}, upperSub[] = {6};
  struct gkyl_range subrange;
  gkyl_sub_range_init(&subrange, &range, lowerSub, upperSub);

  gkyl_array_min_by_cell_range(a1, a1, 12.0, &subrange);

  if (on_gpu) gkyl_array_copy(a1_ho, a1);

  // Check cells outside subrange are unchanged
  for (size_t k=0; k<a1_ho->ncomp; ++k) {
    for (unsigned i=0; i<1; ++i)
      TEST_CHECK( gkyl_compare(a1_ho_d[i*a1_ho->ncomp+k], i*10.0+k, 1e-14) );
    for (unsigned i=6; i<10; ++i)
      TEST_CHECK( gkyl_compare(a1_ho_d[i*a1_ho->ncomp+k], i*10.0+k, 1e-14) );
  }
  
  // Check cells in subrange have min applied
  for (unsigned i=1; i<6; ++i) {
    for (size_t k=0; k<a1_ho->ncomp; ++k) {
      double orig = i*10.0+k;
      double expected = (orig < 12.0) ? orig : 12.0;
      TEST_CHECK( gkyl_compare(a1_ho_d[i*a1_ho->ncomp+k], expected, 1e-14) );
    }
  }

  gkyl_array_release(a1_ho);
  if (on_gpu) gkyl_array_release(a1);

  // Test with different range and threshold
  lower[0] = 1;  upper[0] = 8;
  struct gkyl_range range2;
  gkyl_range_init(&range2, 1, lower, upper);

  struct gkyl_array *a2_ho = gkyl_array_new(GKYL_DOUBLE, 4, range2.volume);
  struct gkyl_array *a2 = a2_ho;
  if (on_gpu) a2 = gkyl_array_cu_dev_new(GKYL_DOUBLE, a2_ho->ncomp, a2_ho->size);

  double *a2_ho_d = a2_ho->data;
  for (unsigned i=0; i<a2_ho->size; ++i) {
    for (size_t k=0; k<a2_ho->ncomp; ++k) 
      a2_ho_d[i*a2_ho->ncomp+k] = i*5.0+k;
  }

  if (on_gpu) gkyl_array_copy(a2, a2_ho);

  gkyl_array_min_by_cell_range(a2, a2, 8.0, &subrange);

  if (on_gpu) gkyl_array_copy(a2_ho, a2);

  for (size_t k=0; k<a2_ho->ncomp; ++k) {
    for (unsigned i=0; i<1; ++i)
      TEST_CHECK( gkyl_compare(a2_ho_d[i*a2_ho->ncomp+k], i*5.0+k, 1e-14) );
    for (unsigned i=6; i<8; ++i)
      TEST_CHECK( gkyl_compare(a2_ho_d[i*a2_ho->ncomp+k], i*5.0+k, 1e-14) );
  }
  
  for (unsigned i=1; i<6; ++i) {
    for (size_t k=0; k<a2_ho->ncomp; ++k) {
      double orig = i*5.0+k;
      double expected = (orig < 8.0) ? orig : 8.0;
      TEST_CHECK( gkyl_compare(a2_ho_d[i*a2_ho->ncomp+k], expected, 1e-14) );
    }
  }

  gkyl_array_release(a2_ho);
  if (on_gpu) gkyl_array_release(a2);
}

void test_array_opcombine()
{
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 1, 10);  
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 1, 10);

  double *a1_d  = a1->data, *a2_d = a2->data;
  for (unsigned i=0; i<a1->size; ++i) {
    a1_d[i] = i*1.0;
    a2_d[i] = i*0.1;
  }

  // a1 <- 0.25*(a1 + 0.5*a2)
  gkyl_array_scale(gkyl_array_set(a1, 0.5, a2), 0.25);

  for (unsigned i=0; i<a1->size; ++i)
    TEST_CHECK( gkyl_compare(a1_d[i], 0.25*0.5*i*0.1, 1e-14) );

  gkyl_array_release(a1);
  gkyl_array_release(a2);
}

void test_array_ops_comp() // more than 1 "component" in array
{
  int nc = 5; // number of "components"
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, nc, 10);

  for (unsigned i=0; i<arr->size; ++i) {
    double *d = gkyl_array_fetch(arr, i);
    for (int k=0; k<nc; ++k)
      d[k] = i*1.0;
  }

  gkyl_array_clear(arr, 1.5f);

  for (unsigned i=0; i<arr->size; ++i) {
    const double *d = gkyl_array_fetch(arr, i);
    for (int k=0; k<nc; ++k)
      TEST_CHECK( gkyl_compare(d[k], 1.5f, 1e-10) );
  }

  gkyl_array_release(arr);
}

void test_array_copy_buffer()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 1, range.volume);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(arr, gkyl_range_idx(&range, iter.idx));
    d[0] = iter.idx[0] + 10.5*iter.idx[1];
  }

  int lower[] = {1, 1}, upper[] = {5, 10};
  struct gkyl_range sub_range;
  gkyl_sub_range_init(&sub_range, &range, lower, upper);

  double *buff = gkyl_malloc(sizeof(double)*sub_range.volume);
  gkyl_array_copy_to_buffer(buff, arr, &sub_range);

  long count = 0;
  gkyl_range_iter_init(&iter, &sub_range);
  while (gkyl_range_iter_next(&iter))
    TEST_CHECK( buff[count++] == iter.idx[0] + 10.5*iter.idx[1] );

  gkyl_array_clear(arr, 0.0);
  // copy back from buffer
  gkyl_array_copy_from_buffer(arr, buff, &sub_range);

  gkyl_range_iter_init(&iter, &sub_range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(arr, gkyl_range_idx(&sub_range, iter.idx));
    TEST_CHECK( d[0]  == iter.idx[0] + 10.5*iter.idx[1] );
  }

  gkyl_array_release(arr);
  gkyl_free(buff);
}

// function for use in the buffer_fn method
GKYL_CU_DH static void
buffer_fn(size_t nc, double *out, const double *inp, void *ctx)
{
  for (size_t i=0; i<nc; ++i)
    out[i] = 2*inp[i];
}
    
void test_array_copy_buffer_fn()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 1, range.volume);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(arr, gkyl_range_idx(&range, iter.idx));
    d[0] = iter.idx[0] + 10.5*iter.idx[1];
  }

  int lower[] = {1, 1}, upper[] = {5, 10};
  struct gkyl_range sub_range;
  gkyl_sub_range_init(&sub_range, &range, lower, upper);

  double *buff = gkyl_malloc(sizeof(double)*sub_range.volume);
  gkyl_array_copy_to_buffer_fn(buff, arr, &sub_range,
    &(struct gkyl_array_copy_func) { .func = buffer_fn, .ctx = 0 }
  );

  long count = 0;
  gkyl_range_iter_init(&iter, &sub_range);
  while (gkyl_range_iter_next(&iter))
    TEST_CHECK( buff[count++] == 2*(iter.idx[0] + 10.5*iter.idx[1]) );

  gkyl_array_clear(arr, 0.0);
  // copy back from buffer
  gkyl_array_copy_from_buffer(arr, buff, &sub_range);

  gkyl_range_iter_init(&iter, &sub_range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(arr, gkyl_range_idx(&sub_range, iter.idx));
    TEST_CHECK( d[0]  == 2*(iter.idx[0] + 10.5*iter.idx[1]) );
  }

  gkyl_array_release(arr);
  gkyl_free(buff);
}

void test_array_flip_copy_buffer_fn()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 1, range.volume);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(arr, gkyl_range_idx(&range, iter.idx));
    d[0] = iter.idx[0] + 10.5*iter.idx[1];
  }

  int lower[] = {1, 1}, upper[] = {5, 10};
  struct gkyl_range sub_range;
  gkyl_sub_range_init(&sub_range, &range, lower, upper);

  double *buff = gkyl_malloc(sizeof(double)*sub_range.volume);

  gkyl_array_flip_copy_to_buffer_fn(buff, arr, 0, &sub_range,
    &(struct gkyl_array_copy_func) { .func = buffer_fn, .ctx = 0 }
  );
  long count = 0;
  gkyl_range_iter_init(&iter, &sub_range);
  while (gkyl_range_iter_next(&iter))
    TEST_CHECK( buff[count++] == 2*((5+1)-iter.idx[0] + 10.5*iter.idx[1]) );

  gkyl_array_flip_copy_to_buffer_fn(buff, arr, 1, &sub_range,     
    &(struct gkyl_array_copy_func) { .func = buffer_fn, .ctx = 0 }
  );
  count = 0;
  gkyl_range_iter_init(&iter, &sub_range);
  while (gkyl_range_iter_next(&iter))
    TEST_CHECK( buff[count++] == 2*(iter.idx[0] + 10.5*((10+1)-iter.idx[1])) );

  gkyl_array_clear(arr, 0.0);
  // copy back from buffer
  gkyl_array_copy_from_buffer(arr, buff, &sub_range);

  gkyl_range_iter_init(&iter, &sub_range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(arr, gkyl_range_idx(&sub_range, iter.idx));
    //printf("%lg %lg\n", d[0], 2*(iter.idx[0] + 10.5*((10+1)-iter.idx[1])) );
    TEST_CHECK( d[0]  == 2*(iter.idx[0] + 10.5*((10+1)-iter.idx[1])) );
  }

  gkyl_array_release(arr);
  gkyl_free(buff);
}

void test_array_copy_range()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);

  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 10, range.volume);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 10, range.volume);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(a1, gkyl_range_idx(&range, iter.idx));
    d[0] = iter.idx[0] + 10.5*iter.idx[1];
  }

  // copy the contents of a1 into a2 over the specified range
  gkyl_array_copy_range(a2, a1, &range);
  
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(a2, gkyl_range_idx(&range, iter.idx));
    TEST_CHECK( d[0]  == iter.idx[0] + 10.5*iter.idx[1] );
  }

  // clear array for second test
  gkyl_array_clear(a2, 0.0);
  // initialize left sub-range
  int lower_l[] = {range.lower[0], range.lower[1]}, upper_l[] = {range.lower[0], shape[1]/2-1};
  struct gkyl_range sub_range_l;
  gkyl_sub_range_init(&sub_range_l, &range, lower_l, upper_l);
  
  // initialize right sub-range
  int lower_r[] = {range.upper[0], shape[1]/2}, upper_r[] = {range.upper[0], range.upper[1]};
  struct gkyl_range sub_range_r;
  gkyl_sub_range_init(&sub_range_r, &range, lower_r, upper_r);

  gkyl_array_copy_range_to_range(a2, a1, &sub_range_r, &sub_range_l);

  gkyl_range_iter_init(&iter, &sub_range_r);
  while (gkyl_range_iter_next(&iter)) {
    int idx_l[GKYL_MAX_DIM];
    idx_l[0] = range.lower[0], idx_l[1] = iter.idx[1]-shape[1]/2;
    double *d = gkyl_array_fetch(a2, gkyl_range_idx(&sub_range_r, iter.idx));
    TEST_CHECK( d[0]  == idx_l[0] + 10.5*idx_l[1] );
    TEST_MSG("Expected: %.13e in cell (%d,%d)", iter.idx[0] + 10.5*iter.idx[1], iter.idx[0], iter.idx[1]);
    TEST_MSG("Produced: %.13e", d[0]);
  }

  gkyl_array_release(a1);
  gkyl_array_release(a2);
}

void test_array_copy_split()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 1, range.volume);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(arr, gkyl_range_idx(&range, iter.idx));
    d[0] = iter.idx[0] + 10.5*iter.idx[1];
  }

  int lower[] = {1, 1}, upper[] = {5, 10};
  struct gkyl_range sub_range;
  gkyl_sub_range_init(&sub_range, &range, lower, upper);

  double *buff = gkyl_malloc(sizeof(double[sub_range.volume]));

  int num_split = 3;
  
  long loc = 0;
  for (int tid=0; tid<num_split; ++tid) {
    struct gkyl_range sr = gkyl_range_split(&sub_range, num_split, tid);
    gkyl_array_copy_to_buffer(buff+loc, arr, &sr);
    
    loc += gkyl_range_split_len(&sr);
  }

  long count = 0;
  gkyl_range_iter_init(&iter, &sub_range);
  while (gkyl_range_iter_next(&iter))
    TEST_CHECK( buff[count++] == iter.idx[0] + 10.5*iter.idx[1] );

  gkyl_array_clear(arr, 0.0);

  loc = 0;
  for (int tid=0; tid<num_split; ++tid) {
    // copy back from buffer
    struct gkyl_range sr = gkyl_range_split(&sub_range, num_split, tid);
    gkyl_array_copy_from_buffer(arr, buff+loc, &sr);
    loc += gkyl_range_split_len(&sr);
  }

  gkyl_range_iter_init(&iter, &sub_range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(arr, gkyl_range_idx(&sub_range, iter.idx));
    TEST_CHECK( d[0]  == iter.idx[0] + 10.5*iter.idx[1] );
  }

  gkyl_array_release(arr);
  gkyl_free(buff);
}

void test_array_shiftc_range_ho() {
  test_array_shiftc_range(false);
}

void test_array_copy_range_to_range_diff_range_dim(bool use_gpu)
{
  // Test array_copy_range_to_range with ranges of different dimensionality.
  // MF 2025/01/22: at the moment this only works for deflating the last
  // dimension.
  int shape_do[] = {10, 20};
  int ncomp = 6;

  int deflate_idx = 3;

  int ndim_do = sizeof(shape_do)/sizeof(shape_do[0]);
  struct gkyl_range range_do;
  gkyl_range_init_from_shape(&range_do, ndim_do, shape_do);

  struct gkyl_array *a_do = mkarr(use_gpu, ncomp, range_do.volume);
  struct gkyl_array *a_do_ho = use_gpu? mkarr(false, ncomp, range_do.volume)
                                      : gkyl_array_acquire(a_do);
  
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range_do);
  while (gkyl_range_iter_next(&iter)) {
    double *a_do_c = gkyl_array_fetch(a_do_ho, gkyl_range_idx(&range_do, iter.idx));
    for (int k=0; k<ncomp; k++) {
      a_do_c[k] = 0.2*k;
      for (int d=0; d<ndim_do; d++)
        a_do_c[k] += d*iter.idx[d];
    }
  }
  gkyl_array_copy(a_do, a_do_ho);

  for (int deflate_dir=1; deflate_dir<ndim_do; deflate_dir++) {

    int remaining_dir = deflate_dir==0? 1 : 0;
    int shape_tar[] = {shape_do[remaining_dir]};

    int ndim_tar = sizeof(shape_tar)/sizeof(shape_tar[0]);

    struct gkyl_range range_tar;
    gkyl_range_init_from_shape(&range_tar, ndim_tar, shape_tar);

    struct gkyl_array *a_tar = mkarr(use_gpu, ncomp, range_tar.volume);
    struct gkyl_array *a_tar_ho = use_gpu? mkarr(false, ncomp, range_tar.volume)
                                         : gkyl_array_acquire(a_tar);

    struct gkyl_range range_do_defl; // Deflated range as a subrange.
    int remove_dir[GKYL_MAX_DIM] = {0};
    int loc_in_dir[GKYL_MAX_DIM] = {0};
    remove_dir[deflate_dir] = 1;
    loc_in_dir[deflate_dir] = deflate_idx;
    gkyl_range_deflate(&range_do_defl, &range_do, remove_dir, loc_in_dir);

    gkyl_array_copy_range_to_range(a_tar, a_do, &range_tar, &range_do_defl);
    
    gkyl_array_copy(a_tar_ho, a_tar);
    gkyl_range_iter_init(&iter, &range_tar);
    while (gkyl_range_iter_next(&iter)) {
      double *a_tar_c = gkyl_array_fetch(a_tar_ho, gkyl_range_idx(&range_tar, iter.idx));
      double ref_val;
      for (int k=0; k<ncomp; k++) {
        ref_val = 0.2*k;
        for (int d=0; d<ndim_do; d++) {
          if (d != deflate_dir)
            ref_val += d*iter.idx[d];
          else
            ref_val += d*deflate_idx;
        }

        TEST_CHECK( a_tar_c[k] == ref_val );
        TEST_MSG("k:%2d Expected: %.13e in cell (%d) | Produced: %.13e", k, a_tar_c[k], iter.idx[remaining_dir], ref_val);
      }
    }

    gkyl_array_release(a_tar);
    gkyl_array_release(a_tar_ho);
  }

  gkyl_array_release(a_do);
  gkyl_array_release(a_do_ho);
}

void test_array_copy_range_to_range_diff_range_dim_ho() {
  test_array_copy_range_to_range_diff_range_dim(false);
}

void test_array_min_by_cell_ho() {
  test_array_min_by_cell(false);
}

void test_array_min_range_ho() {
  test_array_min_range(false);
}

// Cuda specific tests
#ifdef GKYL_HAVE_CUDA

void test_cu_array_clear()
{
  // create host array 
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 1, 10);
  // make device copy
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, 10);
  gkyl_array_copy(a1_cu, a1);

  // clear array on device
  gkyl_array_clear(a1_cu, 0.5);
  // copy from device and check if things are ok
  gkyl_array_copy(a1, a1_cu);
  double *a1_d  = a1->data;  
  for (unsigned i=0; i<a1->size; ++i)
    TEST_CHECK( gkyl_compare(a1_d[i], 0.5, 1e-14) );

  gkyl_array_release(a1);
  gkyl_array_release(a1_cu);
}

void test_cu_array_clear_range()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);

  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 1, range.volume);
  // make device copy of array
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, range.volume);
  gkyl_array_copy(a1_cu, a1);
  gkyl_array_clear_range(a1_cu, 0.5, &range);

  // copy from device and check if things are ok
  gkyl_array_copy(a1, a1_cu);
  double *a1_d  = a1->data;  
  for (unsigned i=0; i<a1->size; ++i)
    TEST_CHECK( gkyl_compare(a1_d[i], 0.5, 1e-14) );

  gkyl_array_release(a1);
  gkyl_array_release(a1_cu);
}

void test_cu_array_accumulate()
{
  // create host arrays 
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 1, 10);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 1, 10);
  // make device copies
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, 10);
  struct gkyl_array *a2_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, 10);

  // initialize data
  double *a1_d  = a1->data, *a2_d = a2->data;
  for (unsigned i=0; i<a1->size; ++i) {
    a1_d[i] = i*1.0;
    a2_d[i] = i*0.1;
  }

  // copy initialized arrays to device
  gkyl_array_copy(a1_cu, a1);
  gkyl_array_copy(a2_cu, a2);

  gkyl_array_accumulate(a1_cu, 0.5, a2_cu);

  // copy from device and check if things are ok
  gkyl_array_copy(a1, a1_cu);
  for (unsigned i=0; i<a1->size; ++i)
    TEST_CHECK( gkyl_compare(a1_d[i], i*1.0+0.5*i*0.1, 1e-14) );

  gkyl_array_release(a1);
  gkyl_array_release(a2);
  gkyl_array_release(a1_cu);
  gkyl_array_release(a2_cu);
}

void test_cu_array_accumulate_range()
{
  int shape[] = {20, 10};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 8, range.volume);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 3, range.volume);

  // make device copies of arrays
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 8, range.volume);
  struct gkyl_array *a2_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3, range.volume);

  // initialize data
  double *a1_d  = a1->data, *a2_d = a2->data;
  for (unsigned i=0; i<a1->size; ++i) {
    for(unsigned c=0; c<a1->ncomp; ++c) 
      a1_d[c+a1->ncomp*i] = i*1.0 + .01*c;
    for(unsigned c=0; c<a2->ncomp; ++c) 
      a2_d[c+a2->ncomp*i] = i*0.1 + .01*c;
  }

  // copy initialized arrays to device
  gkyl_array_copy(a1_cu, a1);
  gkyl_array_copy(a2_cu, a2);

  // a1 = a1 + 0.5*a2
  gkyl_array_accumulate_range(a1_cu, 0.5, a2_cu, &range);

 // copy from device and check if things are ok
  gkyl_array_copy(a1, a1_cu);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *a1d = gkyl_array_fetch(a1, loc);
    for (int c=0; c<a2->ncomp; ++c) {
      TEST_CHECK( gkyl_compare( a1d[c], (.01+.5*.01)*c + loc*1.0 + .5*loc*.1, 1e-14 ) );
    }
    for (int c=a2->ncomp; c<a1->ncomp; ++c)
      TEST_CHECK( a1d[c] == (.01)*c + loc*1.0 );
  }

  // test a2 = a2 + 0.5*a
  gkyl_array_clear(a1_cu, 0.5);
  gkyl_array_clear(a2_cu, 1.5);

  gkyl_array_accumulate_range(a2_cu, 0.5, a1_cu, &range);

 // copy from device and check if things are ok
  gkyl_array_copy(a2, a2_cu);
  gkyl_range_iter_init(&iter, &range);

  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *a2d = gkyl_array_fetch(a2, loc);
    for (int i=0; i<a2->ncomp; ++i)
      TEST_CHECK( a2d[i] == 1.5 + 0.5*0.5 );
  }  

  gkyl_array_release(a1);
  gkyl_array_release(a2);
  gkyl_array_release(a1_cu);
  gkyl_array_release(a2_cu);
}

void test_cu_array_accumulate_offset()
{
  // create host arrays 
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 2, 10);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 3*a1->ncomp, a1->size);
  // make device copies
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, a1->ncomp, a1->size);
  struct gkyl_array *a2_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, a2->ncomp, a2->size);

  // initialize data
  double *a1_d  = a1->data, *a2_d = a2->data;
  for (unsigned i=0; i<a1->size; ++i)
    for (unsigned j=0; j<a1->ncomp; ++j)
      a1_d[i*a1->ncomp+j] = i*1.0+j;

  for (unsigned i=0; i<a2->size; ++i)
    for (unsigned j=0; j<a2->ncomp/a1->ncomp; ++j)
      for (unsigned k=0; k<a1->ncomp; ++k)
        a2_d[i*a2->ncomp+j*a1->ncomp+k] = i*0.1+k;

  // copy initialized arrays to device
  gkyl_array_copy(a1_cu, a1);
  gkyl_array_copy(a2_cu, a2);

  // test a1 = 0.1*a2[a1->ncomp]
  gkyl_array_accumulate_offset(a1_cu, 0.5, a2_cu, 1*a1->ncomp);

  gkyl_array_copy(a1, a1_cu);
  for (unsigned i=0; i<a1->size; ++i)
    for (unsigned j=0; j<a1->ncomp; ++j)
      TEST_CHECK( gkyl_compare(a1_d[i*a1->ncomp+j], i*1.0+j+0.5*(i*0.1+j), 1e-14) );

  // test a2[a1->ncomp] = 0.1*a1
  for (unsigned i=0; i<a1->size; ++i)
    for (unsigned j=0; j<a1->ncomp; ++j)
      a1_d[i*a1->ncomp+j] = i*1.0+j;
  gkyl_array_copy(a1_cu, a1);

  gkyl_array_accumulate_offset(a2_cu, 0.5, a1_cu, 1*a1->ncomp);

  gkyl_array_copy(a2, a2_cu);
  for (unsigned i=0; i<a1->size; ++i) {
    for (unsigned j=0; j<a1->ncomp; ++j) {
      TEST_CHECK( gkyl_compare(a2_d[i*a2->ncomp+0*a1->ncomp+j], i*0.1+j, 1e-14) );
      TEST_CHECK( gkyl_compare(a2_d[i*a2->ncomp+1*a1->ncomp+j], i*0.1+j+0.5*(i*1.0+j), 1e-14) );
      TEST_CHECK( gkyl_compare(a2_d[i*a2->ncomp+2*a1->ncomp+j], i*0.1+j, 1e-14) );
    }
  }

  gkyl_array_release(a1);
  gkyl_array_release(a2);
  gkyl_array_release(a1_cu);
  gkyl_array_release(a2_cu);
}

void test_cu_array_accumulate_offset_range()
{
  int shape[] = {20, 10};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 2, range.volume);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 3*a1->ncomp, range.volume);

  // make device copies of arrays
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, a1->ncomp, range.volume);
  struct gkyl_array *a2_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, a2->ncomp, range.volume);

  // initialize data
  double *a1_d  = a1->data, *a2_d = a2->data;
  for (unsigned i=0; i<a1->size; ++i) {
    for(unsigned c=0; c<a1->ncomp; ++c) 
      a1_d[i*a1->ncomp+c] = i*1.0 + .01*c;
    for(unsigned j=0; j<a2->ncomp/a1->ncomp; ++j)
      for(unsigned c=0; c<a1->ncomp; ++c) 
        a2_d[i*a2->ncomp+j*a1->ncomp+c] = i*0.1 + .01*c;
  }

  // copy initialized arrays to device
  gkyl_array_copy(a1_cu, a1);
  gkyl_array_copy(a2_cu, a2);

  // a1 = a1 + 0.5*a2[1*a1->ncomp]
  gkyl_array_accumulate_offset_range(a1_cu, 0.5, a2_cu, 1*a1->ncomp, &range);

  gkyl_array_copy(a1, a1_cu);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *a1d = gkyl_array_fetch(a1, loc);
    for (int c=0; c<a1->ncomp; ++c)
      TEST_CHECK( gkyl_compare( a1d[c], loc*1.0+.01*c + 0.5*(loc*0.1+0.01*c), 1e-14 ) );
  }

  // test a2[1*a1->ncomp] = a2[1*a1->ncomp] + 0.5*a1
  gkyl_array_clear(a1_cu, 0.5);
  gkyl_array_clear(a2_cu, 1.5);

  gkyl_array_accumulate_offset_range(a2_cu, 0.5, a1_cu, 1*a1->ncomp, &range);

  gkyl_array_copy(a2, a2_cu);
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *a2d = gkyl_array_fetch(a2, loc);
    for (int i=0; i<a1->ncomp; ++i) {
      TEST_CHECK( a2d[i+0*a1->ncomp] == 1.5 );
      TEST_CHECK( a2d[i+1*a1->ncomp] == 1.5+0.5*0.5 );
      TEST_CHECK( a2d[i+2*a1->ncomp] == 1.5 );
    }
  }  

  gkyl_array_release(a1);
  gkyl_array_release(a2);
  gkyl_array_release(a1_cu);
  gkyl_array_release(a2_cu);
}

void test_cu_array_accumulate_range_4d()
{
  int lower[] = { 1, 1, 1, 1 };
  int upper[] = { 46, 46, 32, 32};
  struct gkyl_range range;
  gkyl_range_init(&range, 4, lower, upper);
  
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 8, range.volume);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 3, range.volume);

  // make device copies of arrays
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 8, range.volume);
  struct gkyl_array *a2_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3, range.volume);

  // initialize data
  gkyl_array_clear(a1, 0.5);
  gkyl_array_clear(a2, 1.5);
  
  // copy initialized arrays to device
  gkyl_array_copy(a1_cu, a1);
  gkyl_array_copy(a2_cu, a2);

  int slower[] = { 2, 2, 1, 1 };
  int supper[] = { 45, 45, 32, 32 };
  struct gkyl_range sub_range;
  gkyl_sub_range_init(&sub_range, &range, slower, supper);

  // a1 = a1 + 0.5*a2 (only first 3 components of a1 are modified)
  gkyl_array_accumulate_range(a1_cu, 0.5, a2_cu, &sub_range);

 // copy from device and check if things are ok
  gkyl_array_clear(a1, 0.0); 
  gkyl_array_copy(a1, a1_cu);
  
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &sub_range);
  
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *a1d = gkyl_array_fetch(a1, loc);
    for (int c=0; c<3; ++c)
      TEST_CHECK( a1d[c] == 0.5 + 0.5*1.5 );
    for (int c=3; c<8; ++c)
      TEST_CHECK( a1d[c] == 0.5 );
  }

  gkyl_array_release(a1);
  gkyl_array_release(a2);
  gkyl_array_release(a1_cu);
  gkyl_array_release(a2_cu);
}

void test_cu_array_combine()
{
  // create host arrays 
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 1, 10);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 1, 10);
  struct gkyl_array *b = gkyl_array_new(GKYL_DOUBLE, 1, 10);

  // make device copies
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, 10);
  struct gkyl_array *a2_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, 10);
  struct gkyl_array *b_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, 10);

  double *a1_d  = a1->data, *a2_d = a2->data, *b_d = b->data;
  for (unsigned i=0; i<a1->size; ++i) {
    a1_d[i] = i*1.0;
    a2_d[i] = i*0.1;
    b_d[i] = 10.5;
  }

  // copy initialized arrays to device
  gkyl_array_copy(a1_cu, a1);
  gkyl_array_copy(a2_cu, a2);
  gkyl_array_copy(b_cu, b);

  // b = 0.5*a1 + 2.5*a2
  gkyl_array_set(b_cu, 0.5, a1_cu);
  gkyl_array_accumulate(b_cu, 2.5, a2_cu);

  // copy from device and check if things are ok
  gkyl_array_copy(b, b_cu);
  for (unsigned i=0; i<a1->size; ++i)
    TEST_CHECK( gkyl_compare(b_d[i], 0.5*i*1.0+2.5*i*0.1, 1e-14) );

  gkyl_array_release(a1);
  gkyl_array_release(a2);
  gkyl_array_release(b);
  gkyl_array_release(a1_cu);
  gkyl_array_release(a2_cu);
  gkyl_array_release(b_cu);
}

void test_cu_array_set()
{
  // create host arrays 
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 1, 10);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 1, 10);
  // make device copies
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, 10);
  struct gkyl_array *a2_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, 10);

  // initialize data
  double *a1_d  = a1->data, *a2_d = a2->data;
  for (unsigned i=0; i<a1->size; ++i) {
    a1_d[i] = i*1.0;
    a2_d[i] = i*0.1;
  }

  // copy initialized arrays to device
  gkyl_array_copy(a1_cu, a1);
  gkyl_array_copy(a2_cu, a2);

  gkyl_array_set(a1_cu, 0.5, a2_cu);

  // copy from device and check if things are ok
  gkyl_array_copy(a1, a1_cu);
  for (unsigned i=0; i<a1->size; ++i)
    TEST_CHECK( gkyl_compare(a1_d[i], 0.5*i*0.1, 1e-14) );

  gkyl_array_release(a1);
  gkyl_array_release(a2);
  gkyl_array_release(a1_cu);
  gkyl_array_release(a2_cu);
}

void test_cu_array_set_range()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 8, range.volume);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 3, range.volume);

  // make device copies of arrays
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 8, range.volume);
  struct gkyl_array *a2_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3, range.volume);

  // copy host arrays to device
  gkyl_array_copy(a1_cu, a1);
  gkyl_array_copy(a2_cu, a2);

  // test a1 = 0.5*a2
  gkyl_array_clear(a1_cu, 0.5);
  gkyl_array_clear(a2_cu, 1.5);

  gkyl_array_set_range(a1_cu, 0.5, a2_cu, &range);

 // copy from device and check if things are ok
  gkyl_array_copy(a1, a1_cu);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *a1d = gkyl_array_fetch(a1, loc);
    for (int i=0; i<3; ++i)
      TEST_CHECK( a1d[i] == 0.5*1.5 );
    for (int i=3; i<8; ++i)
      TEST_CHECK( a1d[i] == 0.5);
  }

  // test a2 = 0.5*a1
  gkyl_array_clear(a1_cu, 0.5);
  gkyl_array_clear(a2_cu, 1.5);

  gkyl_array_set_range(a2_cu, 0.5, a1_cu, &range);

  // copy from device and check if things are ok
  gkyl_array_copy(a2, a2_cu);
  gkyl_range_iter_init(&iter, &range);
  
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *a2d = gkyl_array_fetch(a2, loc);
    for (int i=0; i<3; ++i)
      TEST_CHECK( a2d[i] == 0.5*0.5 );
  }  

  gkyl_array_release(a1);
  gkyl_array_release(a2);
  gkyl_array_release(a1_cu);
  gkyl_array_release(a2_cu);
}

void test_cu_array_set_offset()
{
  // create host arrays 
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 2, 10);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 3*a1->ncomp, a1->size);
  // make device copies
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, a1->ncomp, a1->size);
  struct gkyl_array *a2_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, a2->ncomp, a2->size);

  // initialize data
  double *a1_d  = a1->data, *a2_d = a2->data;
  for (unsigned i=0; i<a1->size; ++i)
    for (unsigned j=0; j<a1->ncomp; ++j)
      a1_d[i*a1->ncomp+j] = i*1.0+j;

  for (unsigned i=0; i<a2->size; ++i)
    for (unsigned j=0; j<a2->ncomp/a1->ncomp; ++j)
      for (unsigned k=0; k<a1->ncomp; ++k)
        a2_d[i*a2->ncomp+j*a1->ncomp+k] = i*0.1+k;

  // copy initialized arrays to device
  gkyl_array_copy(a1_cu, a1);
  gkyl_array_copy(a2_cu, a2);

  // test a1 = 0.5*a2[a1->ncomp]
  gkyl_array_set_offset(a1_cu, 0.5, a2_cu, 1*a1->ncomp);

  gkyl_array_copy(a1, a1_cu);
  for (unsigned i=0; i<a1->size; ++i)
    for (unsigned j=0; j<a1->ncomp; ++j)
      TEST_CHECK( gkyl_compare(a1_d[i*a1->ncomp+j], 0.5*(i*0.1+j), 1e-14) );

  // test a2[a1->ncomp] = 0.1*a1
  for (unsigned i=0; i<a1->size; ++i)
    for (unsigned j=0; j<a1->ncomp; ++j)
      a1_d[i*a1->ncomp+j] = i*1.0+j;
  gkyl_array_copy(a1_cu, a1);

  gkyl_array_set_offset(a2_cu, 0.5, a1_cu, 1*a1->ncomp);

  gkyl_array_copy(a2, a2_cu);
  for (unsigned i=0; i<a1->size; ++i) {
    for (unsigned j=0; j<a1->ncomp; ++j) {
      TEST_CHECK( gkyl_compare(a2_d[i*a2->ncomp+0*a1->ncomp+j], i*0.1+j, 1e-14) );
      TEST_CHECK( gkyl_compare(a2_d[i*a2->ncomp+1*a1->ncomp+j], 0.5*(i*1.0+j), 1e-14) );
      TEST_CHECK( gkyl_compare(a2_d[i*a2->ncomp+2*a1->ncomp+j], i*0.1+j, 1e-14) );
    }
  }

  gkyl_array_release(a1);
  gkyl_array_release(a2);
  gkyl_array_release(a1_cu);
  gkyl_array_release(a2_cu);
}

void test_cu_array_set_offset_range()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 2, range.volume);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 3*a1->ncomp, range.volume);

  // make device copies of arrays
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, a1->ncomp, range.volume);
  struct gkyl_array *a2_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, a2->ncomp, range.volume);

  // copy host arrays to device
  gkyl_array_copy(a1_cu, a1);
  gkyl_array_copy(a2_cu, a2);

  // test a1 = 0.5*a2[a1->ncomp]
  gkyl_array_clear(a1_cu, 0.5);
  gkyl_array_clear(a2_cu, 1.5);

  gkyl_array_set_offset_range(a1_cu, 0.5, a2_cu, 1*a1->ncomp, &range);

  gkyl_array_copy(a1, a1_cu);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *a1d = gkyl_array_fetch(a1, loc);
    for (int i=0; i<a1->ncomp; ++i)
      TEST_CHECK( a1d[i] == 0.5*1.5 );
  }

  // test a2[a1->ncomp] = 0.5*a1
  gkyl_array_clear(a1_cu, 0.5);
  gkyl_array_clear(a2_cu, 1.5);

  gkyl_array_set_offset_range(a2_cu, 0.5, a1_cu, 1*a1->ncomp, &range);

  gkyl_array_copy(a2, a2_cu);
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *a2d = gkyl_array_fetch(a2, loc);
    for (int i=0; i<a1->ncomp; ++i) {
      TEST_CHECK( a2d[i+0*a1->ncomp] == 1.5 );
      TEST_CHECK( a2d[i+1*a1->ncomp] == 0.5*0.5 );
      TEST_CHECK( a2d[i+2*a1->ncomp] == 1.5 );
    }
  }  

  gkyl_array_release(a1);
  gkyl_array_release(a2);
  gkyl_array_release(a1_cu);
  gkyl_array_release(a2_cu);
}

void test_cu_array_scale()
{
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 1, 10);
  // make device copies
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, 10);

  // initialize data
  double *a1_d  = a1->data;
  for (unsigned i=0; i<a1->size; ++i) {
    a1_d[i] = i*1.0;
  }

  // copy host arrays to device
  gkyl_array_copy(a1_cu, a1);

  gkyl_array_scale(a1_cu, 0.25);

 // copy from device and check if things are ok
  gkyl_array_copy(a1, a1_cu);
  for (unsigned i=0; i<a1->size; ++i)
    TEST_CHECK( gkyl_compare(a1_d[i], i*0.25, 1e-14) );

  gkyl_array_release(a1);
  gkyl_array_release(a1_cu);
}

void test_cu_array_scale_by_cell()
{
  struct gkyl_array *a1 = mkarr(true, 3, 10);
  struct gkyl_array *s = mkarr(true, 1, 10);
  // Make host copies.
  struct gkyl_array *a1_ho = mkarr(false, a1->ncomp, a1->size);
  struct gkyl_array *s_ho = mkarr(false, s->ncomp, s->size);

  // initialize data
  double *a1_d  = a1_ho->data;
  for (unsigned i=0; i<a1_ho->size; ++i) {
    a1_d[i] = i*1.0;
  }
  double *s_d  = s_ho->data;
  for (unsigned i=0; i<s_ho->size; ++i) {
    s_d[i] = i*1.0;
  }

  // Copy host arrays to device.
  gkyl_array_copy(a1, a1_ho);
  gkyl_array_copy(s, s_ho);

  gkyl_array_scale_by_cell(a1, s);

  // Copy from device and check if things are ok.
  gkyl_array_copy(a1_ho, a1);
  for (unsigned i=0; i<a1_ho->size; ++i) {
    int fact = (i/a1_ho->ncomp);    
    TEST_CHECK( gkyl_compare(a1_d[i], i*1.0*fact, 1e-14) );
  }

  gkyl_array_release(a1);
  gkyl_array_release(s);
  gkyl_array_release(a1_ho);
  gkyl_array_release(s_ho);
}

void test_cu_array_divide_by_cell()
{
  struct gkyl_array *a1 = mkarr(true, 3, 10);
  struct gkyl_array *s = mkarr(true, 1, 10);
  // Make host copies.
  struct gkyl_array *a1_ho = mkarr(false, a1->ncomp, a1->size);
  struct gkyl_array *s_ho = mkarr(false, s->ncomp, s->size);

  // initialize data
  double *a1_d  = a1_ho->data;
  for (unsigned i=0; i<a1_ho->size; ++i) {
    a1_d[i] = i*1.0;
  }
  double *s_d  = s_ho->data;
  for (unsigned i=0; i<s_ho->size; ++i) {
    s_d[i] = i*1.0;
  }

  // Copy host arrays to device.
  gkyl_array_copy(a1, a1_ho);
  gkyl_array_copy(s, s_ho);

  gkyl_array_divide_by_cell(a1, s);

  // Copy from device and check if things are ok.
  gkyl_array_copy(a1_ho, a1);
  for (unsigned i=0; i<a1_ho->size; ++i) {
    int fact = (i/a1_ho->ncomp);    
    TEST_CHECK( gkyl_compare(a1_d[i], i*1.0/fact, 1e-14) );
  }

  gkyl_array_release(a1);
  gkyl_array_release(s);
  gkyl_array_release(a1_ho);
  gkyl_array_release(s_ho);
}

void test_cu_array_shiftc()
{
  double s = -0.5;

  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 3, 10);
  // make device copies
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3, 10);

  // initialize data
  double *a1_d  = a1->data;
  for (unsigned i=0; i<a1->size; ++i) {
    for (size_t k=0; k<a1->ncomp; ++k) a1_d[i*a1->ncomp+k] = i*2.0+k;
  }

  // copy host arrays to device
  gkyl_array_copy(a1_cu, a1);

  gkyl_array_shiftc(a1_cu, s, 0);

  // copy from device and check if things are ok
  gkyl_array_copy(a1, a1_cu);
  TEST_CHECK( gkyl_compare(a1_d[0], 0*1.0+0-0.5, 1e-14) );
  for (unsigned i=0; i<a1->size; ++i) {
    for (size_t k=1; k<a1->ncomp; ++k)
      TEST_CHECK( gkyl_compare(a1_d[i*a1->ncomp+k], i*2.0+k, 1e-14) );
  }

  gkyl_array_release(a1);
  gkyl_array_release(a1_cu);

  // Repeat the test but shifting another coefficient as well.
  int shiftks[] = {0, 2};
  int nks = sizeof(shiftks)/sizeof(shiftks[0]);

  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 4, 8);
  struct gkyl_array *a2_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 4, 8);
  double *a2_d = a2->data;
  for (unsigned i=0; i<a2->size; ++i) {
    for (size_t k=0; k<a2->ncomp; ++k) a2_d[i*a2->ncomp+k] = i*2.0+k;
  }
  gkyl_array_copy(a2_cu, a2);

  for (size_t l=0; l<nks; l++)
    gkyl_array_shiftc(a2_cu, s, shiftks[l]);

  gkyl_array_copy(a2, a2_cu);
  for (unsigned i=0; i<a2->size; ++i) {
    for (size_t k=0; k<a2->ncomp; ++k) {
      bool isshifted = false;
      for (size_t l=0; l<nks; l++) {
        if (shiftks[l]==k) {isshifted = true; break;}
      }
      if (isshifted)
        TEST_CHECK( gkyl_compare(a2_d[i*a2->ncomp+k], i*2.0+k+s, 1e-14) );
      else
        TEST_CHECK( gkyl_compare(a2_d[i*a2->ncomp+k], i*2.0+k, 1e-14) );
    }
  }

  gkyl_array_release(a2);
  gkyl_array_release(a2_cu);
}

void test_cu_array_invert_by_cell()
{
  struct gkyl_array *a1_ho = gkyl_array_new(GKYL_DOUBLE, 3, 8);
  struct gkyl_array *a1 = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3, 8);
  double *a1_ho_d = a1_ho->data;

  // Set values: 1.0, 2.0, 4.0, 5.0, 10.0, 0.5, 0.25, 0.1
  double test_vals[] = {1.0, 2.0, 4.0, 5.0, 10.0, 0.5, 0.25, 0.1};
  double expected_inv[] = {1.0, 0.5, 0.25, 0.2, 0.1, 2.0, 4.0, 10.0};
  
  for (unsigned i=0; i<a1_ho->size; ++i)
    for (size_t k=0; k<a1_ho->ncomp; ++k)
      a1_ho_d[i*a1_ho->ncomp+k] = test_vals[i];

  // Copy to device
  gkyl_array_copy(a1, a1_ho);

  gkyl_array_invert_by_cell(a1, a1);

  // Copy back from device
  gkyl_array_copy(a1_ho, a1);

  // Check inverted values
  for (unsigned i=0; i<a1_ho->size; ++i)
    for (size_t k=0; k<a1_ho->ncomp; ++k)
      TEST_CHECK( gkyl_compare(a1_ho_d[i*a1_ho->ncomp+k], expected_inv[i], 1e-14) );

  gkyl_array_release(a1_ho);
  gkyl_array_release(a1);
}

void test_cu_array_copy_buffer()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 1, range.volume);
  // make device copies of arrays
  struct gkyl_array *arr_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, range.volume);

  // initialize array that will be copied
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(arr, gkyl_range_idx(&range, iter.idx));
    d[0] = iter.idx[0] + 10.5*iter.idx[1];
  }

  // copy host array to device
  gkyl_array_copy(arr_cu, arr);

  int lower[] = {1, 1}, upper[] = {5, 10};
  struct gkyl_range sub_range;
  gkyl_sub_range_init(&sub_range, &range, lower, upper);

  double *buff_cu = gkyl_cu_malloc(sizeof(double)*sub_range.volume);
  gkyl_array_copy_to_buffer(buff_cu, arr_cu, &sub_range);

  gkyl_array_clear(arr, 0.0);
  // copy back from buffer
  gkyl_array_copy_from_buffer(arr_cu, buff_cu, &sub_range);

  // copy from device and check if things are ok
  gkyl_array_copy(arr, arr_cu);
  gkyl_range_iter_init(&iter, &sub_range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(arr, gkyl_range_idx(&sub_range, iter.idx));
    TEST_CHECK( d[0]  == iter.idx[0] + 10.5*iter.idx[1] );
  }

  gkyl_array_release(arr);
  gkyl_array_release(arr_cu);
  gkyl_cu_free(buff_cu);
}

// declare so we can use below
void set_array_copy_fn(struct gkyl_array_copy_func *fn);

void test_cu_array_copy_buffer_fn()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 1, range.volume);

  // make device copy of array
  struct gkyl_array *arr_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, range.volume);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(arr, gkyl_range_idx(&range, iter.idx));
    d[0] = iter.idx[0] + 10.5*iter.idx[1];
  }

  // copy host array to device
  gkyl_array_copy(arr_cu, arr);

  int lower[] = {1, 1}, upper[] = {5, 10};
  struct gkyl_range sub_range;
  gkyl_sub_range_init(&sub_range, &range, lower, upper);

  // create function pointer on device
  struct gkyl_array_copy_func *fn = gkyl_cu_malloc(sizeof(*fn));
  set_array_copy_fn(fn);
  
  double *buff_cu = gkyl_cu_malloc(sizeof(double)*sub_range.volume);
  gkyl_array_copy_to_buffer_fn_cu(buff_cu, arr_cu, &sub_range, fn);
  // copy back from buffer
  gkyl_array_copy_from_buffer(arr_cu, buff_cu, &sub_range);

  // copy from device and check if things are ok
  gkyl_array_copy(arr, arr_cu);
  gkyl_range_iter_init(&iter, &sub_range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(arr, gkyl_range_idx(&sub_range, iter.idx));
    TEST_CHECK( d[0]  == 2*(iter.idx[0] + 10.5*iter.idx[1]) );
  }

  gkyl_array_release(arr);
  gkyl_array_release(arr_cu);
  gkyl_cu_free(buff_cu);
  gkyl_cu_free(fn);
}

void test_cu_array_flip_copy_buffer_fn()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 1, range.volume);

  // make device copy of array
  struct gkyl_array *arr_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, range.volume);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(arr, gkyl_range_idx(&range, iter.idx));
    d[0] = iter.idx[0] + 10.5*iter.idx[1];
  }
  // copy host array to device
  gkyl_array_copy(arr_cu, arr);

  int lower[] = {1, 1}, upper[] = {5, 10};
  struct gkyl_range sub_range;
  gkyl_sub_range_init(&sub_range, &range, lower, upper);

  // create function pointer on device
  struct gkyl_array_copy_func *fn = gkyl_cu_malloc(sizeof(*fn));
  set_array_copy_fn(fn);
  
  double *buff_cu = gkyl_cu_malloc(sizeof(double)*sub_range.volume);

  // test flip copy on first direction of 2D array
  gkyl_array_flip_copy_to_buffer_fn_cu(buff_cu, arr_cu, 0, &sub_range, fn);

  gkyl_array_clear(arr_cu, 0.0);
  // copy back from buffer
  gkyl_array_copy_from_buffer(arr_cu, buff_cu, &sub_range);

  gkyl_array_clear(arr, 0.0);  
  // copy from device and check if things are ok
  gkyl_array_copy(arr, arr_cu);
  gkyl_range_iter_init(&iter, &sub_range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(arr, gkyl_range_idx(&sub_range, iter.idx));
    TEST_CHECK( d[0]  == 2*((5+1)-iter.idx[0] + 10.5*iter.idx[1]) );
  }

  // re-initialize the array
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(arr, gkyl_range_idx(&range, iter.idx));
    d[0] = iter.idx[0] + 10.5*iter.idx[1];
  }
  // copy host array to device
  gkyl_array_copy(arr_cu, arr);
  
  // test flip copy on second direction of 2D array
  gkyl_array_flip_copy_to_buffer_fn(buff_cu, arr_cu, 1, &sub_range, fn);

  gkyl_array_clear(arr_cu, 0.0);
  // copy back from buffer
  gkyl_array_copy_from_buffer(arr_cu, buff_cu, &sub_range);

  gkyl_array_clear(arr, 0.0);
  // copy from device and check if things are ok
  gkyl_array_copy(arr, arr_cu);
  gkyl_range_iter_init(&iter, &sub_range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(arr, gkyl_range_idx(&sub_range, iter.idx));
    //printf("%lg %lg\n", d[0], 2*(iter.idx[0] + 10.5*((10+1)-iter.idx[1])) );
    TEST_CHECK( d[0]  == 2*(iter.idx[0] + 10.5*((10+1)-iter.idx[1])) );
  }

  gkyl_array_release(arr);
  gkyl_array_release(arr_cu);
  gkyl_cu_free(buff_cu);
  gkyl_cu_free(fn);
}

void test_cu_array_copy_range()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);

  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 10, range.volume);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 10, range.volume);

  // make device copies of arrays
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 10, range.volume);
  struct gkyl_array *a2_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 10, range.volume);  
  
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(a1, gkyl_range_idx(&range, iter.idx));
    d[0] = iter.idx[0] + 10.5*iter.idx[1];
  }

  // copy host array to device
  gkyl_array_copy(a1_cu, a1);
  gkyl_array_copy(a2_cu, a2);

  // copy the contents of a1 into a2 over the specified range
  gkyl_array_copy_range_cu(a2_cu, a1_cu, &range);

  /// copy back to host to check contents
  gkyl_array_copy(a2, a2_cu);
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(a2, gkyl_range_idx(&range, iter.idx));
    TEST_CHECK( d[0]  == iter.idx[0] + 10.5*iter.idx[1] );
  }

  // clear array for second test
  gkyl_array_clear_cu(a2_cu, 0.0);
  // initialize left sub-range
  int lower_l[] = {range.lower[0], range.lower[1]}, upper_l[] = {range.lower[0], shape[1]/2-1};
  struct gkyl_range sub_range_l;
  gkyl_sub_range_init(&sub_range_l, &range, lower_l, upper_l);
  
  // initialize right sub-range
  int lower_r[] = {range.upper[0], shape[1]/2}, upper_r[] = {range.upper[0], range.upper[1]};
  struct gkyl_range sub_range_r;
  gkyl_sub_range_init(&sub_range_r, &range, lower_r, upper_r);

  gkyl_array_copy_range_to_range_cu(a2_cu, a1_cu, &sub_range_r, &sub_range_l);

  // copy back to host to check contents
  gkyl_array_copy(a2, a2_cu);
  gkyl_range_iter_init(&iter, &sub_range_r);
  while (gkyl_range_iter_next(&iter)) {
    int idx_l[GKYL_MAX_DIM];
    idx_l[0] = range.lower[0], idx_l[1] = iter.idx[1]-shape[1]/2;
    double *d = gkyl_array_fetch(a2, gkyl_range_idx(&sub_range_r, iter.idx));
    TEST_CHECK( d[0]  == idx_l[0] + 10.5*idx_l[1] );
    TEST_MSG("Expected: %.13e in cell (%d,%d)", iter.idx[0] + 10.5*iter.idx[1], iter.idx[0], iter.idx[1]);
    TEST_MSG("Produced: %.13e", d[0]);
  }

  gkyl_array_release(a1);
  gkyl_array_release(a2);
  gkyl_array_release(a1_cu);
  gkyl_array_release(a2_cu);
}

void test_array_shiftc_range_dev() {
  test_array_shiftc_range(true);
}

void test_array_copy_range_to_range_diff_range_dim_dev() {
  test_array_copy_range_to_range_diff_range_dim(true);
}

void test_array_min_by_cell_dev() {
  test_array_min_by_cell(true);
}

void test_array_min_range_dev() {
  test_array_min_range(true);
}

#endif

TEST_LIST = {
  { "array_clear", test_array_clear },
  { "array_clear_range", test_array_clear_range },
  { "array_accumulate", test_array_accumulate },
  { "array_accumulate_range", test_array_accumulate_range },
  { "array_accumulate_offset", test_array_accumulate_offset },
  { "array_accumulate_offset_range", test_array_accumulate_offset_range },
  { "array_combine", test_array_combine },
  { "array_set", test_array_set },
  { "array_set_range", test_array_set_range },
  { "array_set_offset", test_array_set_offset },
  { "array_set_offset_range", test_array_set_offset_range },
  { "array_scale", test_array_scale },
  { "array_scale_by_cell", test_array_scale_by_cell },
  { "array_invert_by_cell", test_array_invert_by_cell },
  { "array_shiftc", test_array_shiftc },
  { "array_shiftc_range", test_array_shiftc_range_ho },
  { "array_min_by_cell", test_array_min_by_cell_ho },
  { "array_min_range", test_array_min_range_ho },
  { "array_opcombine", test_array_opcombine },
  { "array_ops_comp", test_array_ops_comp },
  { "array_copy_buffer", test_array_copy_buffer },
  { "array_copy_buffer_fn", test_array_copy_buffer_fn },
  { "array_flip_copy_buffer_fn", test_array_flip_copy_buffer_fn },
  { "array_copy_range", test_array_copy_range},
  { "array_copy_split", test_array_copy_split },
  { "array_copy_range_to_range_diff_range_dim", test_array_copy_range_to_range_diff_range_dim_ho},
#ifdef GKYL_HAVE_CUDA
  { "cu_array_clear", test_cu_array_clear},
  { "cu_array_clear_range", test_cu_array_clear_range},
  { "cu_array_accumulate", test_cu_array_accumulate},
  { "cu_array_accumulate_range", test_cu_array_accumulate_range},
  { "cu_array_accumulate_offset", test_cu_array_accumulate_offset},
  { "cu_array_accumulate_offset_range", test_cu_array_accumulate_offset_range},
  { "cu_array_accumulate_range_4d", test_cu_array_accumulate_range_4d  },
  { "cu_array_combine", test_cu_array_combine},
  { "cu_array_set", test_cu_array_set },
  { "cu_array_set_range", test_cu_array_set_range },
  { "cu_array_set_offset", test_cu_array_set_offset },
  { "cu_array_set_offset_range", test_cu_array_set_offset_range },
  { "cu_array_scale", test_cu_array_scale },
  { "cu_array_scale_by_cell", test_cu_array_scale_by_cell },
  { "cu_array_invert_by_cell", test_cu_array_invert_by_cell },
  { "cu_array_shiftc", test_cu_array_shiftc },
  { "cu_array_shiftc_range", test_array_shiftc_range_dev },
  { "cu_array_min_by_cell", test_array_min_by_cell_dev },
  { "cu_array_min_by_cell_range", test_array_min_range_dev },
  { "cu_array_copy_buffer", test_cu_array_copy_buffer },
  { "cu_array_copy_buffer_fn", test_cu_array_copy_buffer_fn },
  { "cu_array_flip_copy_buffer_fn", test_cu_array_flip_copy_buffer_fn },
  { "cu_array_copy_range", test_cu_array_copy_range },
  { "cu_array_copy_range_to_range_diff_range_dim", test_array_copy_range_to_range_diff_range_dim_dev},
#endif
  { NULL, NULL },
};
