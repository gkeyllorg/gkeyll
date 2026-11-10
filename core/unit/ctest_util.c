// Unit tests for utility helpers in gkyl_util.h / util.c and constants in gkyl_const.h
#include <acutest.h>
#include <gkyl_util.h>
#include <gkyl_const.h>
#include <math.h>
#include <sys/wait.h>
#include <unistd.h>

void
test_compare_double()
{
  TEST_CHECK(gkyl_compare_double(1.0, 1.0, 1e-15));
  TEST_CHECK(gkyl_compare_double(0.0, 0.0, 1e-15));
  TEST_CHECK(gkyl_compare_double(1.0, 1.0 + 1e-16, 1e-12));
  TEST_CHECK(!gkyl_compare_double(1.0, 2.0, 1e-12));
  TEST_CHECK(!gkyl_compare_double(1.0, 1.1, 1e-6));
  // near zero
  TEST_CHECK(gkyl_compare_double(0.0, 1e-20, 1e-12));
  // symmetric
  TEST_CHECK(
    gkyl_compare_double(3.0, 3.0 + 1e-13, 1e-9) == gkyl_compare_double(3.0 + 1e-13, 3.0, 1e-9)
  );
}

void
test_round_to_e16_9()
{
  TEST_CHECK(gkyl_round_to_e16_9(0.03483576483233173) == 0.03483576483);
  TEST_CHECK(gkyl_round_to_e16_9(-1.23456789056) == -1.234567891);
  TEST_CHECK(gkyl_round_to_e16_9(0.0) == 0.0);
}

void
test_compare_float()
{
  TEST_CHECK(gkyl_compare_float(1.0f, 1.0f, 1e-6f));
  TEST_CHECK(!gkyl_compare_float(1.0f, 2.0f, 1e-6f));
  TEST_CHECK(gkyl_compare_float(0.0f, 0.0f, 1e-6f));
  TEST_CHECK(gkyl_compare_float(100.0f, 100.0f, 1e-6f));
}

void
test_copy_int_arr()
{
  int in[5] = {1, 2, 3, 4, 5};
  int out[5] = {0};
  gkyl_copy_int_arr(5, in, out);
  for (int i = 0; i < 5; ++i) {
    TEST_CHECK(out[i] == in[i]);
  }

  // partial copy
  int out2[5] = {9, 9, 9, 9, 9};
  gkyl_copy_int_arr(3, in, out2);
  TEST_CHECK(out2[0] == 1 && out2[1] == 2 && out2[2] == 3);
  TEST_CHECK(out2[3] == 9 && out2[4] == 9);
}

void
test_copy_long_arr()
{
  long in[4] = {10L, 20L, 30L, 40L};
  long out[4] = {0};
  gkyl_copy_long_arr(4, in, out);
  for (int i = 0; i < 4; ++i) {
    TEST_CHECK(out[i] == in[i]);
  }
}

void
test_copy_double_arr()
{
  double in[4] = {1.5, -2.5, 3.25, 0.0};
  double out[4] = {0};
  gkyl_copy_double_arr(4, in, out);
  for (int i = 0; i < 4; ++i) {
    TEST_CHECK(out[i] == in[i]);
  }
}

void
test_int_div_up()
{
  TEST_CHECK(gkyl_int_div_up(10, 5) == 2);
  TEST_CHECK(gkyl_int_div_up(11, 5) == 3);
  TEST_CHECK(gkyl_int_div_up(9, 5) == 2);
  TEST_CHECK(gkyl_int_div_up(0, 5) == 0);
  TEST_CHECK(gkyl_int_div_up(5, 5) == 1);
  TEST_CHECK(gkyl_int_div_up(1, 5) == 1);
  TEST_CHECK(gkyl_int_div_up(100, 7) == 15);
}

void
test_minmod_util()
{
  // all same sign positive -> min
  TEST_CHECK(gkyl_minmod(2.0, 3.0, 4.0) == 2.0);
  // all negative -> max (closest to zero)
  TEST_CHECK(gkyl_minmod(-2.0, -3.0, -4.0) == -2.0);
  // mixed signs -> 0
  TEST_CHECK(gkyl_minmod(2.0, -3.0, 4.0) == 0.0);
  TEST_CHECK(gkyl_minmod(-1.0, 2.0, 3.0) == 0.0);
}

void
test_sgn_macro()
{
  TEST_CHECK(GKYL_SGN(5.0) == 1.0);
  TEST_CHECK(GKYL_SGN(-5.0) == -1.0);
  TEST_CHECK(GKYL_SGN(0.0) == 1.0);
}

void
test_minmax_macros()
{
  TEST_CHECK(GKYL_MIN2(3, 5) == 3);
  TEST_CHECK(GKYL_MIN2(5, 3) == 3);
  TEST_CHECK(GKYL_MAX2(3, 5) == 5);
  TEST_CHECK(GKYL_MAX2(5, 3) == 5);
  TEST_CHECK(GKYL_MIN2(-2.5, 1.0) == -2.5);
  TEST_CHECK(GKYL_MAX2(-2.5, 1.0) == 1.0);
}

void
test_tm_trigger()
{
  // Trigger every 1.0 time units.
  struct gkyl_tm_trigger tmt = {.dt = 1.0, .tcurr = 0.0, .curr = 0};

  TEST_CHECK(gkyl_tm_trigger_check_and_bump(&tmt, 0.0) == 1);
  TEST_CHECK(tmt.curr == 1);

  // Below next threshold -> no trigger.
  TEST_CHECK(gkyl_tm_trigger_check_and_bump(&tmt, 0.5) == 0);
  TEST_CHECK(tmt.curr == 1);

  TEST_CHECK(gkyl_tm_trigger_check_and_bump(&tmt, 1.0) == 1);
  TEST_CHECK(tmt.curr == 2);

  TEST_CHECK(gkyl_tm_trigger_check_and_bump(&tmt, 1.9) == 0);
  TEST_CHECK(gkyl_tm_trigger_check_and_bump(&tmt, 2.3) == 1);
  TEST_CHECK(tmt.curr == 3);
}

void
test_search_str_int_pair()
{
  struct gkyl_str_int_pair pairs[] = {{"one", 1}, {"two", 2}, {"three", 3}, {0, 0}};

  TEST_CHECK(gkyl_search_str_int_pair_by_str(pairs, "two", -1) == 2);
  TEST_CHECK(gkyl_search_str_int_pair_by_str(pairs, "one", -1) == 1);
  TEST_CHECK(gkyl_search_str_int_pair_by_str(pairs, "missing", -1) == -1);

  const char *s = gkyl_search_str_int_pair_by_int(pairs, 3, "none");
  TEST_CHECK(strcmp(s, "three") == 0);
  const char *s2 = gkyl_search_str_int_pair_by_int(pairs, 99, "none");
  TEST_CHECK(strcmp(s2, "none") == 0);
}

void
test_constants()
{
  // c = 1/sqrt(mu0 * eps0)
  double c = 1.0 / sqrt(GKYL_MU0 * GKYL_EPSILON0);
  TEST_CHECK(gkyl_compare_double(c, GKYL_SPEED_OF_LIGHT, 1e-6));

  // eV->Kelvin conversion factor consistency
  double ev2k = GKYL_ELEMENTARY_CHARGE / GKYL_BOLTZMANN_CONSTANT;
  TEST_CHECK(gkyl_compare_double(ev2k, GKYL_EV2KELVIN, 1e-9));

  // proton mass much larger than electron mass
  TEST_CHECK(GKYL_PROTON_MASS > 1000.0 * GKYL_ELECTRON_MASS);

  // pi and e values
  TEST_CHECK(gkyl_compare_double(GKYL_PI, M_PI, 1e-14));
  TEST_CHECK(gkyl_compare_double(GKYL_E, M_E, 1e-14));
}

void
test_msgpack_string_array()
{
  char first[] = "Cy";
  const char *keys[] = {first, "", "rho_s"};
  struct gkyl_msgpack_map_elem input = {
    .key = "keys",
    .elem_type = GKYL_MP_STRING_ARRAY,
    .sarray = {.size = 3, .data = keys},
  };

  struct gkyl_msgpack_map_elem *copy = gkyl_msgpack_map_elem_clone(1, &input);
  first[0] = 'X';
  keys[1] = "changed";
  TEST_CHECK(strcmp(copy[0].sarray.data[0], "Cy") == 0);
  TEST_CHECK(strcmp(copy[0].sarray.data[1], "") == 0);

  int lengths[] = {1, 0}, count;
  const struct gkyl_msgpack_map_elem *lists[] = {copy, NULL};
  struct gkyl_msgpack_map_elem *merged = gkyl_msgpack_map_elem_union(2, lengths, lists, &count);
  TEST_CHECK(count == 1);
  gkyl_msgpack_map_elem_release(1, copy);
  TEST_CHECK(merged[0].sarray.size == 3);
  TEST_CHECK(strcmp(merged[0].sarray.data[0], "Cy") == 0);
  TEST_CHECK(strcmp(merged[0].sarray.data[2], "rho_s") == 0);

  // Independent MessagePack encoding: {"keys": ["Cy", "", "rho_s"]}.
  const unsigned char expected[] = {0x81, 0xa4, 'k',  'e', 'y', 's', 0x93, 0xa2, 'C',
                                    'y',  0xa0, 0xa5, 'r', 'h', 'o', '_',  's'};
  struct gkyl_msgpack_data *packed = gkyl_msgpack_create(count, merged);
  TEST_ASSERT(packed != NULL);
  TEST_CHECK(packed->meta_sz == sizeof expected);
  TEST_CHECK(memcmp(packed->meta, expected, sizeof expected) == 0);
  gkyl_msgpack_data_release(packed);

  lists[0] = merged;
  packed = gkyl_msgpack_create_union(2, lengths, lists);
  TEST_ASSERT(packed != NULL);
  TEST_CHECK(packed->meta_sz == sizeof expected);
  TEST_CHECK(memcmp(packed->meta, expected, sizeof expected) == 0);
  gkyl_msgpack_data_release(packed);
  gkyl_msgpack_map_elem_release(count, merged);

  struct gkyl_msgpack_data reference = {.meta = (char *)expected, .meta_sz = sizeof expected};
  struct gkyl_msgpack_map_elem decoded = {.key = "keys", .elem_type = GKYL_MP_STRING_ARRAY};
  gkyl_msgpack_to_map_elem_list(&reference, 1, &decoded);
  TEST_CHECK(decoded.sarray.size == 3);
  TEST_CHECK(strcmp(decoded.sarray.data[0], "Cy") == 0);
  TEST_CHECK(strcmp(decoded.sarray.data[1], "") == 0);
  TEST_CHECK(strcmp(decoded.sarray.data[2], "rho_s") == 0);
  gkyl_msgpack_map_elem_release_value(&decoded);
  TEST_CHECK(decoded.sarray.size == 0 && decoded.sarray.data == NULL);

  const unsigned char empty[] = {0x81, 0xa4, 'k', 'e', 'y', 's', 0x90};
  input.sarray.size = 0;
  input.sarray.data = NULL;
  copy = gkyl_msgpack_map_elem_clone(1, &input);
  packed = gkyl_msgpack_create(1, copy);
  TEST_ASSERT(packed != NULL);
  TEST_CHECK(packed->meta_sz == sizeof empty);
  TEST_CHECK(memcmp(packed->meta, empty, sizeof empty) == 0);

  gkyl_msgpack_to_map_elem_list(packed, 1, &decoded);
  TEST_CHECK(decoded.sarray.size == 0);
  gkyl_msgpack_map_elem_release_value(&decoded);
  gkyl_msgpack_data_release(packed);
  gkyl_msgpack_map_elem_release(1, copy);
}

void
test_msgpack_scalar_copy()
{
  const struct gkyl_msgpack_map_elem input[] = {
    {.key = "bool", .elem_type = GKYL_MP_BOOL, .bval = true},
    {.key = "int", .elem_type = GKYL_MP_INT, .ival = -3},
    {.key = "uint", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 7},
    {.key = "float", .elem_type = GKYL_MP_FLOAT, .fval = 1.25f},
    {.key = "double", .elem_type = GKYL_MP_DOUBLE, .dval = -2.5},
    {.key = "string", .elem_type = GKYL_MP_STRING, .cval = "input"}
  };
  int count = sizeof input / sizeof input[0];
  TEST_CHECK(strcmp(gkyl_msgpack_map_elem_get_string(count, input, "string"), "input") == 0);
  TEST_CHECK(gkyl_msgpack_map_elem_get_string(count, input, "missing") == NULL);
  struct gkyl_msgpack_map_elem *copy = gkyl_msgpack_map_elem_clone(count, input);

  int lengths[] = {count}, merged_count;
  const struct gkyl_msgpack_map_elem *lists[] = {copy};
  struct gkyl_msgpack_map_elem *merged =
    gkyl_msgpack_map_elem_union(1, lengths, lists, &merged_count);
  gkyl_msgpack_map_elem_release(count, copy);

  TEST_CHECK(merged_count == count);
  TEST_CHECK(merged[0].bval == true);
  TEST_CHECK(merged[1].ival == -3);
  TEST_CHECK(merged[2].uval == 7);
  TEST_CHECK(merged[3].fval == 1.25f);
  TEST_CHECK(merged[4].dval == -2.5);
  TEST_CHECK(strcmp(merged[5].cval, "input") == 0);
  for (int i = 0; i < count; ++i) {
    TEST_CHECK(strcmp(merged[i].key, input[i].key) == 0);
    TEST_CHECK(merged[i].elem_type == input[i].elem_type);
  }

  // Decoded values own their storage, while this stack list borrows its key.
  struct gkyl_msgpack_data *packed = gkyl_msgpack_create(merged_count, merged);
  TEST_ASSERT(packed != NULL);
  struct gkyl_msgpack_map_elem decoded = {.key = "string", .elem_type = GKYL_MP_STRING};
  gkyl_msgpack_to_map_elem_list(packed, 1, &decoded);
  TEST_CHECK(strcmp(gkyl_msgpack_map_elem_get_string(1, &decoded, "string"), "input") == 0);
  gkyl_msgpack_map_elem_release_string(1, &decoded, "string");
  TEST_CHECK(decoded.cval == NULL);
  TEST_CHECK(strcmp(decoded.key, "string") == 0);
  gkyl_msgpack_map_elem_release_value(&decoded);
  gkyl_msgpack_data_release(packed);

  // Releasing a cloned string first must also leave full-list release safe.
  gkyl_msgpack_map_elem_release_string(merged_count, merged, "string");
  TEST_CHECK(gkyl_msgpack_map_elem_get_string(merged_count, merged, "string") == NULL);
  gkyl_msgpack_map_elem_release(merged_count, merged);
}

// An error must terminate the child, not the test runner, and identify the key.
static void
check_duplicate_metadata(int method, bool within_group, bool different_value)
{
  int errors[2];
  TEST_ASSERT(pipe(errors) == 0);
  fflush(NULL);
  pid_t child = fork();
  TEST_ASSERT(child >= 0);
  if (child == 0) {
    close(errors[0]);
    if (dup2(errors[1], STDERR_FILENO) < 0) {
      _exit(2);
    }
    close(errors[1]);
    struct gkyl_msgpack_map_elem entries[] = {
      {.key = "duplicate", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 1},
      {.key = "duplicate", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 1}
    };
    if (different_value) {
      entries[1].elem_type = GKYL_MP_STRING;
      entries[1].cval = "different";
    }
    int lengths[] = {within_group ? 2 : 1, within_group ? 0 : 1}, count;
    const struct gkyl_msgpack_map_elem *lists[] = {entries, entries + 1};
    if (method == 0) {
      gkyl_msgpack_map_elem_union(2, lengths, lists, &count);
    } else if (method == 1) {
      gkyl_msgpack_create_union(2, lengths, lists);
    } else {
      gkyl_msgpack_create(2, entries);
    }
    _exit(0);
  }
  close(errors[1]);
  FILE *stream = fdopen(errors[0], "r");
  TEST_ASSERT(stream != NULL);
  char message[128] = {0};
  TEST_CHECK(fgets(message, sizeof message, stream) != NULL);
  fclose(stream);
  int status;
  TEST_ASSERT(waitpid(child, &status, 0) == child);
  TEST_CHECK(WIFEXITED(status) && WEXITSTATUS(status) == EXIT_FAILURE);
  TEST_CHECK(strcmp(message, "Error: Duplicate metadata key \"duplicate\".\n") == 0);
}

void
test_msgpack_unique_keys()
{
  for (int method = 0; method < 3; ++method) {
    for (int different = 0; different < 2; ++different) {
      check_duplicate_metadata(method, true, different);
      if (method < 2) {
        check_duplicate_metadata(method, false, different);
      }
    }
  }

  // Keys are case-sensitive. Empty groups do not disturb order or encoding.
  struct gkyl_msgpack_map_elem entries[] = {
    {.key = "key", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 1},
    {.key = "Key", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 2}
  };
  int lengths[] = {1, 0, 1}, count;
  const struct gkyl_msgpack_map_elem *lists[] = {entries, NULL, entries + 1};
  struct gkyl_msgpack_map_elem *merged = gkyl_msgpack_map_elem_union(3, lengths, lists, &count);
  TEST_CHECK(count == 2);
  TEST_CHECK(strcmp(merged[0].key, "key") == 0);
  TEST_CHECK(strcmp(merged[1].key, "Key") == 0);
  const unsigned char expected[] = {0x82, 0xa3, 'k', 'e', 'y', 1, 0xa3, 'K', 'e', 'y', 2};
  struct gkyl_msgpack_data *packed[] = {
    gkyl_msgpack_create(count, merged), gkyl_msgpack_create_union(3, lengths, lists)
  };
  for (int i = 0; i < 2; ++i) {
    TEST_ASSERT(packed[i] != NULL);
    TEST_CHECK(packed[i]->meta_sz == sizeof expected);
    TEST_CHECK(memcmp(packed[i]->meta, expected, sizeof expected) == 0);
    gkyl_msgpack_data_release(packed[i]);
  }
  gkyl_msgpack_map_elem_release(count, merged);
  struct gkyl_msgpack_data *empty = gkyl_msgpack_create(0, NULL);
  TEST_ASSERT(empty != NULL);
  TEST_CHECK(empty->meta_sz == 1 && (unsigned char)empty->meta[0] == 0x80);
  gkyl_msgpack_data_release(empty);
}

void
test_msgpack_set_string()
{
  const struct gkyl_msgpack_map_elem input = {
    .key = "basis_type",
    .elem_type = GKYL_MP_STRING,
    .cval = "old basis",
  };
  struct gkyl_msgpack_map_elem *owned = gkyl_msgpack_map_elem_clone(1, &input);
  gkyl_msgpack_map_elem_set_string(1, owned, "basis_type", owned[0].cval + 4);
  TEST_CHECK(strcmp(owned[0].cval, "basis") == 0);
  TEST_CHECK(strcmp(input.cval, "old basis") == 0);
  gkyl_msgpack_map_elem_set_string(1, owned, "basis_type", "tensor");
  TEST_CHECK(strcmp(owned[0].cval, "tensor") == 0);
  TEST_CHECK(strcmp(owned[0].key, "basis_type") == 0);
  gkyl_msgpack_map_elem_set_string(1, owned, "missing", "ignored");
  TEST_CHECK(strcmp(owned[0].cval, "tensor") == 0);
  gkyl_msgpack_map_elem_release(1, owned);
}

TEST_LIST = {
  {"msgpack_unique_keys", test_msgpack_unique_keys},
  {"msgpack_set_string", test_msgpack_set_string},
  {"msgpack_string_array", test_msgpack_string_array},
  {"msgpack_scalar_copy", test_msgpack_scalar_copy},
  {"compare_double", test_compare_double},
  {"round_to_e16_9", test_round_to_e16_9},
  {"compare_float", test_compare_float},
  {"copy_int_arr", test_copy_int_arr},
  {"copy_long_arr", test_copy_long_arr},
  {"copy_double_arr", test_copy_double_arr},
  {"int_div_up", test_int_div_up},
  {"minmod_util", test_minmod_util},
  {"sgn_macro", test_sgn_macro},
  {"minmax_macros", test_minmax_macros},
  {"tm_trigger", test_tm_trigger},
  {"search_str_int_pair", test_search_str_int_pair},
  {"constants", test_constants},
  {NULL, NULL}
};
