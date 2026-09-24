#include <acutest.h>
#include <gkyl_array_rio_format_desc.h>
#include <gkyl_dynvec.h>

#include <math.h>

void
test_1()
{
  gkyl_dynvec dv = gkyl_dynvec_new(GKYL_DOUBLE, 3);

  TEST_CHECK( gkyl_dynvec_elem_type(dv) == GKYL_DOUBLE );
  TEST_CHECK( gkyl_dynvec_ncomp(dv) == 3 );

  double out[3];
  TEST_CHECK( gkyl_dynvec_size(dv) == 0 );
  TEST_CHECK( gkyl_dynvec_getlast(dv, out) == false );

  // add some data
  for (int i=0; i<2000; ++i)
    gkyl_dynvec_append(dv, 0.1*i, (double[3]) { i, i+1, i+2 } );
  TEST_CHECK( gkyl_dynvec_size(dv) == 2000 );

  TEST_CHECK( gkyl_dynvec_getlast_tm(dv) == 1999*0.1 );
  TEST_CHECK( gkyl_dynvec_getlast(dv, out) == true );
  TEST_CHECK( out[0] == 1999 );
  TEST_CHECK( out[1] == 2000 );
  TEST_CHECK( out[2] == 2001 );

  for (int i=0; i<2000; ++i) {
    double d[3];
    TEST_CHECK( gkyl_dynvec_get(dv, i, d) == true );
    TEST_CHECK( d[0] == i );
    TEST_CHECK( d[1] == i+1 );
    TEST_CHECK( d[2] == i+2 );

    TEST_CHECK( gkyl_dynvec_get_tm(dv, i) == i*0.1 );
  }

  TEST_CHECK( gkyl_dynvec_capacity(dv) > 2000 );
  TEST_CHECK( gkyl_dynvec_capacity(dv) < 5000 );

  gkyl_dynvec_reserve_more(dv, 5000);
  TEST_CHECK( gkyl_dynvec_capacity(dv) > 2000+5000 );
  
  for (int i=0; i<2000; ++i) {
    double d[3];
    TEST_CHECK( gkyl_dynvec_get(dv, i, d) == true );
    TEST_CHECK( d[0] == i );
    TEST_CHECK( d[1] == i+1 );
    TEST_CHECK( d[2] == i+2 );

    TEST_CHECK( gkyl_dynvec_get_tm(dv, i) == i*0.1 );
  }  

  gkyl_dynvec_clear(dv);
  TEST_CHECK( gkyl_dynvec_size(dv) == 0 );
  TEST_CHECK( gkyl_dynvec_getlast(dv, out) == false );

  gkyl_dynvec dv2 = gkyl_dynvec_acquire(dv);
  
  gkyl_dynvec_release(dv);
  gkyl_dynvec_release(dv2);
}

void
test_2()
{
  // store user-defined struct
  struct euler { double rho, rhou, rhov; };
  gkyl_dynvec dv = gkyl_dynvec_new(GKYL_USER, sizeof(struct euler));

  // add some data
  for (int i=0; i<2000; ++i)
    gkyl_dynvec_append(dv, 0.1*i, &(struct euler) { i, i+1, i+2 });
  
  TEST_CHECK( gkyl_dynvec_size(dv) == 2000 );

  for (int i=0; i<2000; ++i) {
    struct euler eu;
    TEST_CHECK( gkyl_dynvec_get(dv, i, &eu) == true );
    TEST_CHECK( eu.rho == i );
    TEST_CHECK( eu.rhou == i+1 );
    TEST_CHECK( eu.rhov == i+2 );

    TEST_CHECK( gkyl_dynvec_get_tm(dv, i) == i*0.1 );
  }  

  gkyl_dynvec_release(dv);
}

void
test_3()
{
  gkyl_dynvec dv = gkyl_dynvec_new(GKYL_DOUBLE, 3);

  double out[3];
  TEST_CHECK( gkyl_dynvec_size(dv) == 0 );
  TEST_CHECK( gkyl_dynvec_getlast(dv, out) == false );

  // add some data
  for (int i=0; i<2000; ++i)
    gkyl_dynvec_append(dv, 0.1*i, (double[3]) { i, i+1, i+2 } );
  
  TEST_CHECK( gkyl_dynvec_size(dv) == 2000 );

  TEST_CHECK( gkyl_dynvec_getlast_tm(dv) == 1999*0.1 );
  TEST_CHECK( gkyl_dynvec_getlast(dv, out) == true );
  TEST_CHECK( out[0] == 1999 );
  TEST_CHECK( out[1] == 2000 );
  TEST_CHECK( out[2] == 2001 );

  for (int i=0; i<2000; ++i) {
    double d[3];
    TEST_CHECK( gkyl_dynvec_get(dv, i, d) == true );
    TEST_CHECK( d[0] == i );
    TEST_CHECK( d[1] == i+1 );
    TEST_CHECK( d[2] == i+2 );

    TEST_CHECK( gkyl_dynvec_get_tm(dv, i) == i*0.1 );
  }

  gkyl_dynvec_clear_all_but(dv, 3);
  TEST_CHECK( gkyl_dynvec_size(dv) == 3 );

  TEST_CHECK( gkyl_dynvec_getlast(dv, out) == true );
  TEST_CHECK( out[0] == 1999 );
  TEST_CHECK( out[1] == 2000 );
  TEST_CHECK( out[2] == 2001 );
  TEST_CHECK( gkyl_dynvec_getlast_tm(dv) == 1999*0.1 );

  TEST_CHECK( gkyl_dynvec_get(dv, 0, out) == true );  
  TEST_CHECK( out[0] == 2000-3 );
  TEST_CHECK( out[1] == 2000-3+1 );
  TEST_CHECK( out[2] == 2000-3+2 );
  TEST_CHECK( gkyl_dynvec_get_tm(dv, 0) == (2000-3)*0.1 );
  
  TEST_CHECK( gkyl_dynvec_get(dv, 1, out) == true );
  TEST_CHECK( out[0] == 2000-2 );
  TEST_CHECK( out[1] == 2000-2+1 );
  TEST_CHECK( out[2] == 2000-2+2 );
  TEST_CHECK( gkyl_dynvec_get_tm(dv, 1) == (2000-2)*0.1 );
  
  TEST_CHECK( gkyl_dynvec_get(dv, 2, out) == true );
  TEST_CHECK( out[0] == 2000-1 );
  TEST_CHECK( out[1] == 2000-1+1 );
  TEST_CHECK( out[2] == 2000-1+2 );
  TEST_CHECK( gkyl_dynvec_get_tm(dv, 2) == (2000-1)*0.1 );
  
  gkyl_dynvec_release(dv);
}

void
test_4()
{
  gkyl_dynvec dv = gkyl_dynvec_new(GKYL_DOUBLE, 3);

  double out[3];
  TEST_CHECK( gkyl_dynvec_size(dv) == 0 );
  TEST_CHECK( gkyl_dynvec_getlast(dv, out) == false );

  gkyl_dynvec_clear_all_but(dv, 3);

  TEST_CHECK( gkyl_dynvec_size(dv) == 0 );
  TEST_CHECK( gkyl_dynvec_getlast(dv, out) == false );  

  // add some data
  for (int i=0; i<20; ++i)
    gkyl_dynvec_append(dv, 0.1*i, (double[3]) { i, i+1, i+2 } );
  
  TEST_CHECK( gkyl_dynvec_size(dv) == 20 );

  gkyl_dynvec_clear_all_but(dv, 30); // clearing more elements

  TEST_CHECK( gkyl_dynvec_size(dv) == 20 );

  TEST_CHECK( gkyl_dynvec_getlast_tm(dv) == 19*0.1 );
  TEST_CHECK( gkyl_dynvec_getlast(dv, out) == true );
  TEST_CHECK( out[0] == 19 );
  TEST_CHECK( out[1] == 20 );
  TEST_CHECK( out[2] == 21 );

  for (int i=0; i<20; ++i) {
    double d[3];
    TEST_CHECK( gkyl_dynvec_get(dv, i, d) == true );
    TEST_CHECK( d[0] == i );
    TEST_CHECK( d[1] == i+1 );
    TEST_CHECK( d[2] == i+2 );

    TEST_CHECK( gkyl_dynvec_get_tm(dv, i) == i*0.1 );
  }
  
  gkyl_dynvec_release(dv);
}

void
test_io()
{
  gkyl_dynvec dv = gkyl_dynvec_new(GKYL_DOUBLE, 3);

  double out[3];

  for (int i=0; i<1000; ++i) {
    out[0] = cos(0.1*i);
    out[1] = sin(0.1*i);
    out[2] = cos(0.1*i)*sin(0.1*i);
    gkyl_dynvec_append(dv, i*0.1, out);
  }

  // write and clear
  gkyl_dynvec_write(dv, "ctest_dynvec_test_io_1.gkyl");
  gkyl_dynvec_clear(dv);

  struct gkyl_dynvec_etype_ncomp enc = { };
  enc = gkyl_dynvec_read_ncomp("ctest_dynvec_test_io_1.gkyl");
  TEST_CHECK( enc.ncomp == 3 );
  TEST_CHECK( enc.type == GKYL_DOUBLE );

  int file_type = gkyl_get_gkyl_file_type("ctest_dynvec_test_io_1.gkyl");
  TEST_CHECK( 2 == file_type );  
  
  bool res = gkyl_dynvec_read(dv, "ctest_dynvec_test_io_1.gkyl");
  TEST_CHECK( res );

  TEST_CHECK( gkyl_dynvec_size(dv) == 1000 );

  for (int i=0; i<1000; ++i) {
    gkyl_dynvec_get(dv, i, out);
    TEST_CHECK( out[0] == cos(0.1*i) );
    TEST_CHECK( out[1] == sin(0.1*i) );
    TEST_CHECK( out[2] == cos(0.1*i)*sin(0.1*i) );

    TEST_CHECK( gkyl_dynvec_get_tm(dv, i) == i*0.1 );
  }

  gkyl_dynvec_clear(dv);

  // add some elements before reading
  for (int i=0; i<10; ++i) {
    out[0] = 0.1*i;
    out[1] = 0.2*i;
    out[2] = 0.3*i;
    gkyl_dynvec_append(dv, i*0.1, out);
  }

  TEST_CHECK( gkyl_dynvec_size(dv) == 10 );

  res = gkyl_dynvec_read(dv, "ctest_dynvec_test_io_1.gkyl");
  TEST_CHECK( res );

  for (int i=0; i<10; ++i) {
    gkyl_dynvec_get(dv, i, out);
    TEST_CHECK( out[0] == 0.1*i );
    TEST_CHECK( out[1] == 0.2*i );
    TEST_CHECK( out[2] == 0.3*i );

    TEST_CHECK( gkyl_dynvec_get_tm(dv, i) == i*0.1 );
  }

  TEST_CHECK( gkyl_dynvec_size(dv) == 10+1000 );

  for (int i=10; i<1000+10; ++i) {
    gkyl_dynvec_get(dv, i, out);
    TEST_CHECK( out[0] == cos(0.1*(i-10)) );
    TEST_CHECK( out[1] == sin(0.1*(i-10)) );
    TEST_CHECK( out[2] == cos(0.1*(i-10))*sin(0.1*(i-10)) );

    TEST_CHECK( gkyl_dynvec_get_tm(dv, i) == (i-10)*0.1 );
  }  

  gkyl_dynvec_release(dv);
}

void
test_io_2()
{
  gkyl_dynvec dv = gkyl_dynvec_new(GKYL_DOUBLE, 3);

  double out[3];

  for (int i=0; i<1000; ++i) {
    out[0] = cos(0.1*i);
    out[1] = sin(0.1*i);
    out[2] = cos(0.1*i)*sin(0.1*i);
    gkyl_dynvec_append(dv, i*0.1, out);
  }

  // write twice to the same file
  gkyl_dynvec_write(dv, "ctest_dynvec_test_io_2.gkyl");
  gkyl_dynvec_awrite(dv, "ctest_dynvec_test_io_2.gkyl");
  
  gkyl_dynvec_clear(dv);

  bool res = gkyl_dynvec_read(dv, "ctest_dynvec_test_io_2.gkyl");
  TEST_CHECK( res );

  TEST_CHECK( gkyl_dynvec_size(dv) == 2000 );

  for (int i=0; i<1000; ++i) {
    gkyl_dynvec_get(dv, i, out);
    TEST_CHECK( out[0] == cos(0.1*i) );
    TEST_CHECK( out[1] == sin(0.1*i) );
    TEST_CHECK( out[2] == cos(0.1*i)*sin(0.1*i) );

    TEST_CHECK( gkyl_dynvec_get_tm(dv, i) == i*0.1 );
  }
  for (int j=1000; j<2000; ++j) {
    int i = j-1000;
    gkyl_dynvec_get(dv, i, out);
    TEST_CHECK( out[0] == cos(0.1*i) );
    TEST_CHECK( out[1] == sin(0.1*i) );
    TEST_CHECK( out[2] == cos(0.1*i)*sin(0.1*i) );

    TEST_CHECK( gkyl_dynvec_get_tm(dv, i) == i*0.1 );
  }  

  gkyl_dynvec_release(dv);
}

void
test_to_array()
{
  gkyl_dynvec dv = gkyl_dynvec_new(GKYL_DOUBLE, 3);

  double out[3] = { 0.0 };
  for (int i=0; i<10; ++i) {
    out[0] = 0.1*i;
    out[1] = 0.2*i;
    out[2] = 0.3*i;
    gkyl_dynvec_append(dv, i*0.1, out);
  }

  struct gkyl_array *tm_mesh = gkyl_array_new(GKYL_DOUBLE, 1, gkyl_dynvec_size(dv));
  struct gkyl_array *dyn_data = gkyl_array_new(GKYL_DOUBLE, 3, gkyl_dynvec_size(dv));

  gkyl_dynvec_to_array(dv, tm_mesh, dyn_data);

  for (int i=0; i<10; ++i) {
    const double *tm_i = gkyl_array_cfetch(tm_mesh, i);
    const double *dd_i = gkyl_array_cfetch(dyn_data, i);

    TEST_CHECK( tm_i[0] == i*0.1 );

    TEST_CHECK( dd_i[0] == i*0.1 );
    TEST_CHECK( dd_i[1] == i*0.2 );
    TEST_CHECK( dd_i[2] == i*0.3 );
  }
  
  gkyl_dynvec_release(dv);
  gkyl_array_release(tm_mesh);
  gkyl_array_release(dyn_data);
}

TEST_LIST = {
  { "test_1", test_1 },
  { "test_2", test_2 },
  { "test_3", test_3 },
  { "test_4", test_4 },
  { "test_io", test_io },
  { "test_io_2", test_io_2 },
  { "test_to_array", test_to_array },
  { NULL, NULL },
};
