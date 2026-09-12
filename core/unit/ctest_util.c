// Unit tests for utility helpers in gkyl_util.h / util.c and constants in gkyl_const.h
#include <acutest.h>
#include <gkyl_util.h>
#include <gkyl_const.h>
#include <math.h>

void
test_compare_double()
{
  TEST_CHECK( gkyl_compare_double(1.0, 1.0, 1e-15) );
  TEST_CHECK( gkyl_compare_double(0.0, 0.0, 1e-15) );
  TEST_CHECK( gkyl_compare_double(1.0, 1.0 + 1e-16, 1e-12) );
  TEST_CHECK( !gkyl_compare_double(1.0, 2.0, 1e-12) );
  TEST_CHECK( !gkyl_compare_double(1.0, 1.1, 1e-6) );
  // near zero
  TEST_CHECK( gkyl_compare_double(0.0, 1e-20, 1e-12) );
  // symmetric
  TEST_CHECK( gkyl_compare_double(3.0, 3.0+1e-13, 1e-9) == gkyl_compare_double(3.0+1e-13, 3.0, 1e-9) );
}

void
test_compare_float()
{
  TEST_CHECK( gkyl_compare_float(1.0f, 1.0f, 1e-6f) );
  TEST_CHECK( !gkyl_compare_float(1.0f, 2.0f, 1e-6f) );
  TEST_CHECK( gkyl_compare_float(0.0f, 0.0f, 1e-6f) );
  TEST_CHECK( gkyl_compare_float(100.0f, 100.0f, 1e-6f) );
}

void
test_copy_int_arr()
{
  int in[5] = { 1, 2, 3, 4, 5 };
  int out[5] = { 0 };
  gkyl_copy_int_arr(5, in, out);
  for (int i=0; i<5; ++i)
    TEST_CHECK( out[i] == in[i] );

  // partial copy
  int out2[5] = { 9, 9, 9, 9, 9 };
  gkyl_copy_int_arr(3, in, out2);
  TEST_CHECK( out2[0]==1 && out2[1]==2 && out2[2]==3 );
  TEST_CHECK( out2[3]==9 && out2[4]==9 );
}

void
test_copy_long_arr()
{
  long in[4] = { 10L, 20L, 30L, 40L };
  long out[4] = { 0 };
  gkyl_copy_long_arr(4, in, out);
  for (int i=0; i<4; ++i)
    TEST_CHECK( out[i] == in[i] );
}

void
test_copy_double_arr()
{
  double in[4] = { 1.5, -2.5, 3.25, 0.0 };
  double out[4] = { 0 };
  gkyl_copy_double_arr(4, in, out);
  for (int i=0; i<4; ++i)
    TEST_CHECK( out[i] == in[i] );
}

void
test_int_div_up()
{
  TEST_CHECK( gkyl_int_div_up(10, 5) == 2 );
  TEST_CHECK( gkyl_int_div_up(11, 5) == 3 );
  TEST_CHECK( gkyl_int_div_up(9, 5) == 2 );
  TEST_CHECK( gkyl_int_div_up(0, 5) == 0 );
  TEST_CHECK( gkyl_int_div_up(5, 5) == 1 );
  TEST_CHECK( gkyl_int_div_up(1, 5) == 1 );
  TEST_CHECK( gkyl_int_div_up(100, 7) == 15 );
}

void
test_minmod_util()
{
  // all same sign positive -> min
  TEST_CHECK( gkyl_minmod(2.0, 3.0, 4.0) == 2.0 );
  // all negative -> max (closest to zero)
  TEST_CHECK( gkyl_minmod(-2.0, -3.0, -4.0) == -2.0 );
  // mixed signs -> 0
  TEST_CHECK( gkyl_minmod(2.0, -3.0, 4.0) == 0.0 );
  TEST_CHECK( gkyl_minmod(-1.0, 2.0, 3.0) == 0.0 );
}

void
test_sgn_macro()
{
  TEST_CHECK( GKYL_SGN(5.0) == 1.0 );
  TEST_CHECK( GKYL_SGN(-5.0) == -1.0 );
  TEST_CHECK( GKYL_SGN(0.0) == 1.0 );
}

void
test_minmax_macros()
{
  TEST_CHECK( GKYL_MIN2(3, 5) == 3 );
  TEST_CHECK( GKYL_MIN2(5, 3) == 3 );
  TEST_CHECK( GKYL_MAX2(3, 5) == 5 );
  TEST_CHECK( GKYL_MAX2(5, 3) == 5 );
  TEST_CHECK( GKYL_MIN2(-2.5, 1.0) == -2.5 );
  TEST_CHECK( GKYL_MAX2(-2.5, 1.0) == 1.0 );
}

void
test_tm_trigger()
{
  // Trigger every 1.0 time units.
  struct gkyl_tm_trigger tmt = { .dt = 1.0, .tcurr = 0.0, .curr = 0 };

  TEST_CHECK( gkyl_tm_trigger_check_and_bump(&tmt, 0.0) == 1 );
  TEST_CHECK( tmt.curr == 1 );

  // Below next threshold -> no trigger.
  TEST_CHECK( gkyl_tm_trigger_check_and_bump(&tmt, 0.5) == 0 );
  TEST_CHECK( tmt.curr == 1 );

  TEST_CHECK( gkyl_tm_trigger_check_and_bump(&tmt, 1.0) == 1 );
  TEST_CHECK( tmt.curr == 2 );

  TEST_CHECK( gkyl_tm_trigger_check_and_bump(&tmt, 1.9) == 0 );
  TEST_CHECK( gkyl_tm_trigger_check_and_bump(&tmt, 2.3) == 1 );
  TEST_CHECK( tmt.curr == 3 );
}

void
test_search_str_int_pair()
{
  struct gkyl_str_int_pair pairs[] = {
    { "one", 1 },
    { "two", 2 },
    { "three", 3 },
    { 0, 0 }
  };

  TEST_CHECK( gkyl_search_str_int_pair_by_str(pairs, "two", -1) == 2 );
  TEST_CHECK( gkyl_search_str_int_pair_by_str(pairs, "one", -1) == 1 );
  TEST_CHECK( gkyl_search_str_int_pair_by_str(pairs, "missing", -1) == -1 );

  const char *s = gkyl_search_str_int_pair_by_int(pairs, 3, "none");
  TEST_CHECK( strcmp(s, "three") == 0 );
  const char *s2 = gkyl_search_str_int_pair_by_int(pairs, 99, "none");
  TEST_CHECK( strcmp(s2, "none") == 0 );
}

void
test_constants()
{
  // c = 1/sqrt(mu0 * eps0)
  double c = 1.0/sqrt(GKYL_MU0*GKYL_EPSILON0);
  TEST_CHECK( gkyl_compare_double(c, GKYL_SPEED_OF_LIGHT, 1e-6) );

  // eV->Kelvin conversion factor consistency
  double ev2k = GKYL_ELEMENTARY_CHARGE/GKYL_BOLTZMANN_CONSTANT;
  TEST_CHECK( gkyl_compare_double(ev2k, GKYL_EV2KELVIN, 1e-9) );

  // proton mass much larger than electron mass
  TEST_CHECK( GKYL_PROTON_MASS > 1000.0*GKYL_ELECTRON_MASS );

  // pi and e values
  TEST_CHECK( gkyl_compare_double(GKYL_PI, M_PI, 1e-14) );
  TEST_CHECK( gkyl_compare_double(GKYL_E, M_E, 1e-14) );
}

TEST_LIST = {
  { "compare_double", test_compare_double },
  { "compare_float", test_compare_float },
  { "copy_int_arr", test_copy_int_arr },
  { "copy_long_arr", test_copy_long_arr },
  { "copy_double_arr", test_copy_double_arr },
  { "int_div_up", test_int_div_up },
  { "minmod_util", test_minmod_util },
  { "sgn_macro", test_sgn_macro },
  { "minmax_macros", test_minmax_macros },
  { "tm_trigger", test_tm_trigger },
  { "search_str_int_pair", test_search_str_int_pair },
  { "constants", test_constants },
  { NULL, NULL },
};
