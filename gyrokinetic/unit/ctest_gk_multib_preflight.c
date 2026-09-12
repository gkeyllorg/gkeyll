#include <acutest.h>
#include <gkyl_gyrokinetic_multib.h>

#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static struct gkyl_gk_block_geom *
make_pair(bool extended0, bool extended1)
{
  struct gkyl_gk_block_geom *bgeom = gkyl_gk_block_geom_new(2, 2);
  for (int b=0; b<2; ++b) {
    struct gkyl_gk_block_geom_info info = {
      .lower = { 0.0, -1.0 }, .upper = { 1.0, 1.0 }, .cells = { 2, 2 },
      .geometry = {
        .geometry_id = GKYL_GEOMETRY_TOKAMAK,
        .tok_grid_info = {
          .ftype = GKYL_GEOMETRY_TOKAMAK_CORE_L,
          .straight_xpt_ray = b ? extended1 : extended0,
        },
      },
    };
    for (int d=0; d<2; ++d)
      for (int e=0; e<2; ++e)
        info.connections[d][e] = (struct gkyl_target_edge) {
          .dir = d, .edge = GKYL_PHYSICAL,
        };
    info.connections[0][b ? 0 : 1] = (struct gkyl_target_edge) {
      .bid = 1-b, .dir = 0,
      .edge = b ? GKYL_UPPER_POSITIVE : GKYL_LOWER_POSITIVE,
    };
    gkyl_gk_block_geom_set_block(bgeom, b, &info);
  }
  return bgeom;
}

// Assert the machine-readable evidence emitted by the actual C entry point.
static void
check_preflight_report(const struct gkyl_gyrokinetic_multib *inp,
  int expected, const char *summary)
{
  FILE *capture = tmpfile();
  TEST_ASSERT(capture != NULL);
  fflush(stderr);
  int saved_stderr = dup(STDERR_FILENO);
  TEST_ASSERT(saved_stderr >= 0);
  TEST_ASSERT(dup2(fileno(capture), STDERR_FILENO) >= 0);
  int actual = gkyl_gyrokinetic_multib_app_geometry_preflight(inp);
  fflush(stderr);
  TEST_ASSERT(dup2(saved_stderr, STDERR_FILENO) >= 0);
  close(saved_stderr);

  char report[16384];
  rewind(capture);
  size_t len = fread(report, 1, sizeof(report)-1, capture);
  report[len] = '\0';
  fclose(capture);
  TEST_CHECK(actual == expected);
  TEST_MSG("preflight returned %d; expected %d\n%s", actual, expected, report);
  TEST_CHECK(strstr(report, summary) != NULL);
  TEST_MSG("missing summary: %s\nactual report: %s", summary, report);
}

static void
check_constructors_reject(struct gkyl_gk_block_geom *bgeom)
{
  // No communicator is supplied: both constructors must reject BEFORE any
  // communicator access, reference acquisition, or geometry construction.
  struct gkyl_gyrokinetic_multib inp = { .cdim = 2, .gk_block_geom = bgeom };
  TEST_CHECK(gkyl_gyrokinetic_multib_app_new_geom(&inp) == NULL);
  TEST_CHECK(gkyl_gyrokinetic_multib_app_new(&inp) == NULL);
}

static void
test_uniform(void)
{
  setenv("GKYL_TOK_STRICT_SEAM_PARTICIPATION", "1", 1);
  for (int extended=0; extended<2; ++extended) {
    struct gkyl_gk_block_geom *bgeom = make_pair(extended, extended);
    struct gkyl_gyrokinetic_multib inp = { .cdim = 2, .gk_block_geom = bgeom };
    check_preflight_report(&inp, 1,
      "GKYL_GEOMETRY_PREFLIGHT status=PASS scope=declaration num_blocks=2 "
      "strict=1 interfaces_examined=1 mixed=0");
    gkyl_gk_block_geom_release(bgeom);
  }
  unsetenv("GKYL_TOK_STRICT_SEAM_PARTICIPATION");
}

static void
test_mixed_default_and_strict(void)
{
  unsetenv("GKYL_TOK_STRICT_SEAM_PARTICIPATION");
  struct gkyl_gk_block_geom *bgeom = make_pair(true, false);
  struct gkyl_gyrokinetic_multib inp = { .cdim = 2, .gk_block_geom = bgeom };
  check_preflight_report(&inp, 1,
    "GKYL_GEOMETRY_PREFLIGHT status=PASS scope=declaration num_blocks=2 "
    "strict=0 interfaces_examined=1 mixed=1");
  setenv("GKYL_TOK_STRICT_SEAM_PARTICIPATION", "1", 1);
  check_preflight_report(&inp, 0,
    "GKYL_GEOMETRY_PREFLIGHT status=FAIL scope=declaration num_blocks=2 "
    "strict=1 interfaces_examined=1 mixed=1");
  check_constructors_reject(bgeom);
  // Strict rejection must leave the caller's declaration usable and owned.
  unsetenv("GKYL_TOK_STRICT_SEAM_PARTICIPATION");
  TEST_CHECK(gkyl_gk_block_geom_check_consistency(bgeom) == 1);
  gkyl_gk_block_geom_release(bgeom);
}

static void
test_malformed(void)
{
  unsetenv("GKYL_TOK_STRICT_SEAM_PARTICIPATION");
  const struct gkyl_target_edge malformed[] = {
    { .bid = 1, .dir = 0, .edge = 0 }, // unspecified
    { .bid = 2, .dir = 0, .edge = GKYL_LOWER_POSITIVE }, // invalid block
    { .bid = -1, .dir = 0, .edge = GKYL_LOWER_POSITIVE },
    { .bid = 1, .dir = 2, .edge = GKYL_LOWER_POSITIVE }, // invalid direction
    { .bid = 1, .dir = -1, .edge = GKYL_LOWER_POSITIVE },
    { .bid = 1, .dir = 1, .edge = GKYL_LOWER_POSITIVE }, // wrong valid direction
    { .bid = 1, .dir = 0, .edge = GKYL_PHYSICAL+1 }, // invalid enum
    { .bid = 1, .dir = 0, .edge = GKYL_LOWER_NEGATIVE }, // nonreciprocal
  };
  for (size_t i=0; i<sizeof(malformed)/sizeof(malformed[0]); ++i) {
    struct gkyl_gk_block_geom *bgeom = make_pair(false, false);
    struct gkyl_gk_block_geom_info info = *gkyl_gk_block_geom_get_block(bgeom, 0);
    info.connections[0][1] = malformed[i];
    gkyl_gk_block_geom_set_block(bgeom, 0, &info);
    check_constructors_reject(bgeom);
    gkyl_gk_block_geom_release(bgeom);
  }
}

static void
test_orientation_and_rotated_restriction(void)
{
  setenv("GKYL_TOK_STRICT_SEAM_PARTICIPATION", "1", 1);
  struct gkyl_gk_block_geom *bgeom = make_pair(true, true);
  struct gkyl_gk_block_geom_info a = *gkyl_gk_block_geom_get_block(bgeom, 0);
  struct gkyl_gk_block_geom_info b = *gkyl_gk_block_geom_get_block(bgeom, 1);
  a.connections[0][1].edge = GKYL_LOWER_NEGATIVE;
  b.connections[0][0].edge = GKYL_UPPER_NEGATIVE;
  gkyl_gk_block_geom_set_block(bgeom, 0, &a);
  gkyl_gk_block_geom_set_block(bgeom, 1, &b);
  struct gkyl_gyrokinetic_multib inp = { .cdim = 2, .gk_block_geom = bgeom };
  check_preflight_report(&inp, 1,
    "GKYL_GEOMETRY_PREFLIGHT status=PASS scope=declaration num_blocks=2 "
    "strict=1 interfaces_examined=1 mixed=0");

  // This connection is reciprocal but rotates logical directions. Preserve
  // the core topology restriction: gyrokinetic ranges assume aligned axes.
  a.connections[0][1] = (struct gkyl_target_edge) {
    .bid = 1, .dir = 1, .edge = GKYL_LOWER_POSITIVE,
  };
  b.connections[0][0].edge = GKYL_PHYSICAL;
  b.connections[1][0] = (struct gkyl_target_edge) {
    .bid = 0, .dir = 0, .edge = GKYL_UPPER_POSITIVE,
  };
  gkyl_gk_block_geom_set_block(bgeom, 0, &a);
  gkyl_gk_block_geom_set_block(bgeom, 1, &b);
  check_constructors_reject(bgeom);
  gkyl_gk_block_geom_release(bgeom);
  unsetenv("GKYL_TOK_STRICT_SEAM_PARTICIPATION");
}

static void
test_same_side_radial_edges(void)
{
  setenv("GKYL_TOK_STRICT_SEAM_PARTICIPATION", "1", 1);
  for (int side=0; side<2; ++side) {
    struct gkyl_gk_block_geom *bgeom = make_pair(true, true);
    for (int bid=0; bid<2; ++bid) {
      struct gkyl_gk_block_geom_info info = *gkyl_gk_block_geom_get_block(bgeom, bid);
      info.connections[0][0].edge = GKYL_PHYSICAL;
      info.connections[0][1].edge = GKYL_PHYSICAL;
      info.connections[0][side] = (struct gkyl_target_edge) {
        .bid = 1-bid, .dir = 0,
        .edge = side ? GKYL_UPPER_POSITIVE : GKYL_LOWER_POSITIVE,
      };
      gkyl_gk_block_geom_set_block(bgeom, bid, &info);
    }
    struct gkyl_gyrokinetic_multib inp = { .cdim = 2, .gk_block_geom = bgeom };
    check_preflight_report(&inp, 1,
      "GKYL_GEOMETRY_PREFLIGHT status=PASS scope=declaration num_blocks=2 "
      "strict=1 interfaces_examined=1 mixed=0");
    gkyl_gk_block_geom_release(bgeom);
  }
  unsetenv("GKYL_TOK_STRICT_SEAM_PARTICIPATION");
}

static void
test_missing_or_dimension_mismatch(void)
{
  unsetenv("GKYL_TOK_STRICT_SEAM_PARTICIPATION");
  check_preflight_report(NULL, 0,
    "GKYL_GEOMETRY_PREFLIGHT status=FAIL scope=declaration num_blocks=0 "
    "strict=0 interfaces_examined=0 mixed=0");
  TEST_CHECK(gkyl_gyrokinetic_multib_app_new_geom(NULL) == NULL);
  TEST_CHECK(gkyl_gyrokinetic_multib_app_new(NULL) == NULL);
  check_constructors_reject(NULL);

  struct gkyl_gk_block_geom *bgeom = make_pair(false, false);
  struct gkyl_gyrokinetic_multib inp = { .cdim = 3, .gk_block_geom = bgeom };
  check_preflight_report(&inp, 0,
    "GKYL_GEOMETRY_PREFLIGHT status=FAIL scope=declaration num_blocks=2");
  TEST_CHECK(gkyl_gyrokinetic_multib_app_new_geom(&inp) == NULL);
  TEST_CHECK(gkyl_gyrokinetic_multib_app_new(&inp) == NULL);
  gkyl_gk_block_geom_release(bgeom);
}

TEST_LIST = {
  { "uniform", test_uniform },
  { "mixed_default_and_strict", test_mixed_default_and_strict },
  { "malformed", test_malformed },
  { "orientation_and_rotated_restriction", test_orientation_and_rotated_restriction },
  { "same_side_radial_edges", test_same_side_radial_edges },
  { "missing_or_dimension_mismatch", test_missing_or_dimension_mismatch },
  { NULL, NULL },
};
