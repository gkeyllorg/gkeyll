// Defect 3.2: material preflight must check mapped psi, not computational psi.
// This exercises the constructors: the public geometry_preflight entry point
// only checks declarations. A PASS here certifies preflight, not a built grid.
#include <acutest.h>
#include <gkyl_comm_priv.h>
#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_tok_geo.h>

#include <errno.h>
#include <signal.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>

static void
plate_lower(double s, double *rz)
{ rz[0] = 1.5966+(1.6888-1.5966)*s; rz[1] = -1.1421+(-0.8781+1.1421)*s; }

static void
plate_upper(double s, double *rz)
{ rz[0] = 1.2686+(1.1886-1.2686)*s; rz[1] = -1.0520+(-0.7294+1.0520)*s; }

static void
shift_map(double t, const double *xn, double *out, void *ctx)
{ out[0] = xn[0] + *(const double *) ctx; }

// Terminate a passing child at the first communicator access. This bounds the
// fixture without relying on a NULL dereference or masking a sanitizer error.
static int
stop_at_rank(struct gkyl_comm *comm, int *rank)
{
  fflush(stderr);
  _exit(83);
}

static struct gkyl_gk_block_geom *
make_block(double lo, double hi, enum gkyl_position_map_id id, double *shift)
{
  struct gkyl_gk_block_geom *bg = gkyl_gk_block_geom_new(2, 1);
  struct gkyl_gk_block_geom_info info = {
    .lower = { lo, -0.01 }, .upper = { hi, 0.01 }, .cells = { 4, 4 },
    .geometry = {
      .geometry_id = GKYL_GEOMETRY_TOKAMAK,
      .efit_info = {
        .filepath = "gyrokinetic/data/eqdsk/asdex.geqdsk",
        .rz_poly_order = 2, .flux_poly_order = 1,
      },
      .tok_grid_info = {
        .ftype = GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO,
        .rmin = 0.0, .rmax = 5.0, .rclose = 2.5, .rright = 2.5, .rleft = 0.7,
        .zmin = -1.3, .zmax = 1.0, .zmin_left = -1.0, .zmin_right = -0.9,
        // Without plate_spec, coverage succeeds without inspecting a plate.
        .plate_spec = true,
        .plate_func_lower = plate_lower, .plate_func_upper = plate_upper,
      },
      .position_map_info = { .id = id },
    },
  };
  if (shift) {
    info.geometry.position_map_info.maps[0] = shift_map;
    info.geometry.position_map_info.ctxs[0] = shift;
  }
  for (int d=0; d<2; ++d)
    for (int e=0; e<2; ++e)
      info.connections[d][e] = (struct gkyl_target_edge) {
        .dir = d, .edge = GKYL_PHYSICAL,
      };
  gkyl_gk_block_geom_set_block(bg, 0, &info);
  return bg;
}

static void
check_preflight(double lo, double hi, enum gkyl_position_map_id id,
  double *shift, bool expected_pass, bool advisory)
{
  // Both constructors check material before accessing the communicator. The
  // sentinel communicator stops a passing child at that boundary. Require
  // both the preflight verdict and the expected clean termination.
  for (int full=0; full<2; ++full) {
    FILE *capture = tmpfile();
    TEST_ASSERT(capture != NULL);
    fflush(NULL);
    pid_t pid = fork();
    TEST_ASSERT(pid >= 0);
    if (pid == 0) {
      if (dup2(fileno(capture), STDERR_FILENO) < 0)
        _exit(120);
      signal(SIGSEGV, SIG_DFL);
      alarm(60);
      unsetenv("ADJUST_IF_EXCEEDING_WALL");
      unsetenv("EXTEND_TO_LIMITER");
      unsetenv("GKYL_TOK_STRICT_SEAM_PARTICIPATION");
      struct gkyl_gk_block_geom *bg = make_block(lo, hi, id, shift);
      struct gkyl_comm_priv sentinel = { .get_rank = stop_at_rank };
      struct gkyl_gyrokinetic_multib inp = {
        .cdim = 2, .gk_block_geom = bg, .comm = &sentinel.pub_comm,
      };
      gkyl_gyrokinetic_multib_app *app = full ?
        gkyl_gyrokinetic_multib_app_new(&inp) : gkyl_gyrokinetic_multib_app_new_geom(&inp);
      gkyl_gk_block_geom_release(bg);
      _exit(app ? 121 : 0);
    }
    int status = 0;
    pid_t waited;
    do { waited = waitpid(pid, &status, 0); } while (waited < 0 && errno == EINTR);
    TEST_ASSERT(waited == pid);
    char report[16384];
    rewind(capture);
    size_t len = fread(report, 1, sizeof(report)-1, capture);
    report[len] = '\0';
    TEST_CHECK(!ferror(capture) && len < sizeof(report)-1);
    fclose(capture);
    bool passed = strstr(report, "GKYL_MATERIAL_PREFLIGHT status=PASS\n") != NULL;
    bool failed = strstr(report, "GKYL_MATERIAL_PREFLIGHT status=FAIL block=0\n") != NULL;
    TEST_CHECK(passed == expected_pass && failed == !expected_pass);
    TEST_MSG("constructor=%s map=%d bounds=[%.6f,%.6f] expected=%s\n%s",
      full ? "full" : "geometry", id, lo, hi, expected_pass ? "PASS" : "FAIL", report);
    TEST_CHECK((strstr(report, "reason=position_map_not_resolved_yet") != NULL) == advisory);
    TEST_MSG("advisory=%d\n%s", advisory, report);
    TEST_CHECK(WIFEXITED(status) && WEXITSTATUS(status) == (expected_pass ? 83 : 0));
    TEST_MSG("PASS must reach the rank sentinel; rejection must return NULL first, status=%d\n%s", status, report);
  }
}

static void
test_plate_fixture(void)
{
  struct gkyl_gk_block_geom *bg = make_block(0.155, 0.165, GKYL_PMAP_USER_INPUT, NULL);
  const struct gkyl_gk_block_geom_info *bi = gkyl_gk_block_geom_get_block(bg, 0);
  struct gkyl_tok_geo *geo = gkyl_tok_geo_new(&bi->geometry.efit_info, &bi->geometry.tok_grid_info);
  TEST_ASSERT(geo != NULL);
  // Explicit native coverage anchors, independent of the map/preflight helper.
  const double psi[] = { 0.155, 0.157, 0.163, 0.165, 0.167, 0.173, 0.175 };
  for (size_t i=0; i<sizeof(psi)/sizeof(psi[0]); ++i) {
    TEST_CHECK(gkyl_tok_geo_check_plate_coverage(geo, &bi->geometry.tok_grid_info, psi[i]) == (i < 5));
    TEST_MSG("ASDEX native plate coverage at psi=%.6f", psi[i]);
  }
  gkyl_tok_geo_release(geo);
  gkyl_gk_block_geom_release(bg);
}

static void
test_mapped_rejection(void)
{
  double shift = 0.008;
  for (int id=GKYL_PMAP_USER_INPUT; id<=GKYL_PMAP_USER_INPUT_W_DERIVATIVE; ++id)
    check_preflight(0.155, 0.165, id, &shift, false, false);
}

static void
test_mapped_acceptance(void)
{
  double shift = -0.008;
  for (int id=GKYL_PMAP_USER_INPUT; id<=GKYL_PMAP_USER_INPUT_W_DERIVATIVE; ++id)
    check_preflight(0.165, 0.175, id, &shift, true, false);
}

static void
test_identity_controls(void)
{
  double shift = 0.0;
  check_preflight(0.155, 0.165, GKYL_PMAP_USER_INPUT, &shift, true, false);
  check_preflight(0.155, 0.165, GKYL_PMAP_USER_INPUT, NULL, true, false);
  check_preflight(0.165, 0.175, GKYL_PMAP_USER_INPUT, NULL, false, false);
}

static void
test_unresolved_maps_advisory(void)
{
  const enum gkyl_position_map_id ids[] = {
    GKYL_PMAP_CONSTANT_DB_POLYNOMIAL, GKYL_PMAP_CONSTANT_DB_NUMERIC, GKYL_PMAP_XPT_COMPRESSION,
  };
  // No map exists yet. Raw plate failure cannot decide mapped coverage.
  for (size_t i=0; i<sizeof(ids)/sizeof(ids[0]); ++i)
    check_preflight(0.165, 0.175, ids[i], NULL, true, true);
}

TEST_LIST = {
  { "plate_fixture", test_plate_fixture },
  { "mapped_rejection", test_mapped_rejection },
  { "mapped_acceptance", test_mapped_acceptance },
  { "identity_controls", test_identity_controls },
  { "unresolved_maps_advisory", test_unresolved_maps_advisory },
  { NULL, NULL },
};
