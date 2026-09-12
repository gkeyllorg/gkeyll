// Gate 4 regressions: attributing a wall refusal.
//
// Written BEFORE the fix, which is the campaign's rule.
//
// What this gate first thought, and why it was wrong
// --------------------------------------------------
// The initial reading was that `tok_wall_segment_inside` tests a straight chord
// where the real edge is curved, and so falsely refuses coarse theta grids.
// Both halves were measured and BOTH ARE FALSE at p1:
//
//   * gate4/probe_edge_linearity.c: 3D serendipity p1 is trilinear, so along an
//     edge it reduces to a linear function. Deviation of the mid node from the
//     chord midpoint is 2.16e-16 -- roundoff. At p1 the chord IS the represented
//     edge, and tok_wall_curve_inside is identical to tok_wall_segment_inside.
//
//   * The ASDEX contained band refused at ntheta=4 is a REAL excursion. Checked
//     independently against the 31-vertex outline: both endpoints inside, by
//     60.20 mm and 21.91 mm -- and the chord between them leaves the vessel by
//     1.00 mm. At p1 that chord is the cell edge, so the discrete geometry
//     genuinely puts an edge outside the machine. Refusing it is correct.
//
// So containment is inherently resolution-dependent, and that is right: a
// coarser representation of a curved flux surface can genuinely exit the vessel.
//
// What IS missing
// ---------------
// Two very different situations print the same line today:
//
//   declaration error      refused at EVERY resolution   (leg ends outside)
//   discretisation excursion  refused only when coarse   (1.00 mm out at nth=4)
//
// A user cannot tell which they have without running a refinement sweep by hand,
// and the excursion depth -- the one number that separates "move your cuts" from
// "refine theta" -- is never reported. That is this gate's first target.

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

#include <acutest.h>
#include <gkyl_basis.h>
#include <gkyl_efit.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_tok.h>
#include <gkyl_position_map.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_tok_geo.h>
// The containment predicates are private: this test is checking the library's
// internal machinery deliberately, not its public surface.
#include <gkyl_tok_geo_priv.h>

struct sol_case {
  const char *name;
  double psi_lo, psi_hi;
  double zmin_left, zmin_right;
};

// The two gate-2 ASDEX SOL declarations, verbatim.
static const struct sol_case contained  = { "contained",  0.1500, 0.1520, -1.0, -0.9 };
static const struct sol_case historical = { "historical", 0.1600, 0.17501, -1.2, -1.0 };

// Build one geometry in a child process: a wall violation calls abort(), so a
// trial cannot share an address space with the test runner.
// Returns true if the domain was ACCEPTED.
static bool
sol_accepted(const struct sol_case *c, int ntheta, char *report, size_t report_size)
{
  FILE *capture = tmpfile();
  TEST_ASSERT(capture != NULL);
  int errfd = fileno(capture);

  fflush(NULL);
  pid_t pid = fork();
  TEST_ASSERT(pid >= 0);
  if (pid == 0) {
    dup2(errfd, STDERR_FILENO);

    double clower[] = { c->psi_lo, -0.01, -M_PI+1e-14 };
    double cupper[] = { c->psi_hi,  0.01,  M_PI-1e-14 };
    int ccells[] = { 4, 1, ntheta };

    struct gkyl_rect_grid cgrid;
    struct gkyl_range clocal, clocal_ext;
    struct gkyl_basis cbasis;
    int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
    gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
    gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
    gkyl_cart_modal_serendip(&cbasis, 3, 1);

    struct gkyl_efit_inp efit_inp = {
      .filepath = "gyrokinetic/data/eqdsk/asdex.geqdsk",
      .rz_poly_order = 2,
      .flux_poly_order = 1,
    };
    struct gkyl_tok_geo_grid_inp ginp = {
      .ftype = GKYL_GEOMETRY_TOKAMAK_LSN_SOL,
      .rmin = 0.0, .rmax = 5.0,
      .rclose = 2.5, .rright = 2.5, .rleft = 0.7,
      .zmin = -1.3, .zmax = 1.0,
      .zmin_left = c->zmin_left,
      .zmin_right = c->zmin_right,
    };
    struct gkyl_position_map *pmap = gkyl_position_map_null_new();
    struct gkyl_gk_geometry_inp geometry_inp = {
      .geometry_id = GKYL_GEOMETRY_TOKAMAK,
      .efit_info = efit_inp,
      .tok_grid_info = ginp,
      .position_map = pmap,
      .grid = cgrid, .local = clocal, .local_ext = clocal_ext,
      .global = clocal, .global_ext = clocal_ext, .basis = cbasis,
      .geo_grid = cgrid, .geo_local = clocal, .geo_local_ext = clocal_ext,
      .geo_global = clocal, .geo_global_ext = clocal_ext, .geo_basis = cbasis,
    };

    struct gk_geometry *geom = gkyl_gk_geometry_tok_new(&geometry_inp);
    _exit(geom == NULL ? 1 : 0);
  }

  int status = 0;
  TEST_ASSERT(waitpid(pid, &status, 0) == pid);

  rewind(capture);
  size_t n = fread(report, 1, report_size-1, capture);
  report[n] = '\0';
  fclose(capture);

  return WIFEXITED(status) && WEXITSTATUS(status) == 0;
}

// Pull the diagnostic's own account of which test failed, so a refusal can be
// attributed rather than guessed.
static const char *
failing_scope(const char *report)
{
  if (strstr(report, "scope=segment_theta")) return "segment_theta (straight chord along theta)";
  if (strstr(report, "scope=segment_psi"))   return "segment_psi (straight chord along psi)";
  if (strstr(report, "scope=corner_node"))   return "corner_node (a node is genuinely outside)";
  if (strstr(report, "scope=radial_boundary_curve")) return "radial_boundary_curve (exact edge)";
  return "(no scope reported)";
}

static const int resolutions[] = { 4, 8, 16, 32 };

// THE regression. When a domain is refused because its DISCRETE edge leaves the
// vessel -- rather than because the declaration itself is outside -- the
// diagnostic must say so and must quantify the excursion. Today it prints the
// same `scope=segment_theta` line in both cases and no depth at all, so the two
// are indistinguishable without a hand-run refinement sweep.
//
// Fails today: no excursion figure is emitted.
static void
test_refusal_reports_excursion_depth(void)
{
  char report[16384];
  bool ok = sol_accepted(&contained, 4, report, sizeof(report));

  TEST_CHECK(!ok);
  TEST_MSG("the contained band at ntheta=4 should still be refused: its cell "
    "edge is genuinely 1.00 mm outside the vessel");

  // It must name the excursion depth, so "refine theta" is distinguishable
  // from "your declaration is outside the machine".
  TEST_CHECK(strstr(report, "excursion_m=") != NULL);
  TEST_MSG("refusal does not report how far outside the edge goes. Without it, "
    "a 1 mm discretisation excursion and a 257 mm declaration error read "
    "identically.\nscope was: %s\n%s", failing_scope(report), report);

  // And it must say that both endpoints were inside -- the fact that makes it a
  // segment excursion rather than a node placed outside.
  TEST_CHECK(strstr(report, "endpoints=inside") != NULL);
  TEST_MSG("refusal does not distinguish 'both nodes inside, edge bulges out' "
    "from 'a node is outside'.\n%s", report);
}

// The contained band's refusal must VANISH under refinement, and the historical
// one must not. That contrast is the evidence that the first is discretisation
// and the second is declaration -- it is the property a user would otherwise
// have to discover by hand.
static void
test_refinement_separates_discretisation_from_declaration(void)
{
  char report[16384];
  bool contained_ok[4], historical_ok[4];
  for (size_t k=0; k<sizeof(resolutions)/sizeof(*resolutions); ++k) {
    contained_ok[k]  = sol_accepted(&contained,  resolutions[k], report, sizeof(report));
    historical_ok[k] = sol_accepted(&historical, resolutions[k], report, sizeof(report));
  }

  // Contained: refused coarse, accepted once resolvable.
  TEST_CHECK(!contained_ok[0]);
  TEST_MSG("contained band accepted at ntheta=4, but its edge is 1.00 mm outside");
  for (size_t k=1; k<sizeof(resolutions)/sizeof(*resolutions); ++k) {
    TEST_CHECK(contained_ok[k]);
    TEST_MSG("contained band refused at ntheta=%d; it is representable there",
      resolutions[k]);
  }

  // Historical: refused everywhere. If a change ever makes this pass, the
  // containment check has been weakened, not sharpened.
  for (size_t k=0; k<sizeof(resolutions)/sizeof(*resolutions); ++k) {
    TEST_CHECK(!historical_ok[k]);
    TEST_MSG("historical declaration ACCEPTED at ntheta=%d. Its leg ends are "
      "outside the vessel; accepting it means the check was weakened.",
      resolutions[k]);
  }
}

// Pin the measurement that killed the original plan, so it is not re-proposed.
// At p1 the represented edge is straight, so the exact curve test and the chord
// test cannot disagree. Anyone "fixing" the chord test at p1 is changing nothing.
static void
test_curve_and_chord_are_identical_at_p1(void)
{
  struct gkyl_basis b;
  gkyl_cart_modal_serendip(&b, 3, 1);
  double c[8];
  srand(4321);
  double worst = 0.0;
  for (int t=0; t<200; ++t) {
    for (int i=0; i<8; ++i) c[i] = 2.0*((double) rand()/RAND_MAX) - 1.0;
    double e0[3]={-1.0,-1.0,-1.0}, em[3]={-1.0,-1.0,0.0}, e1[3]={-1.0,-1.0,1.0};
    double v0=b.eval_expand(e0,c), vm=b.eval_expand(em,c), v1=b.eval_expand(e1,c);
    double scale = fmax(1e-300, fabs(v0)+fabs(v1));
    worst = fmax(worst, fabs(vm - 0.5*(v0+v1))/scale);
  }
  TEST_CHECK(worst < 1.0e-14);
  TEST_MSG("p1 edge deviates from its chord by %.3e relative. If this ever "
    "becomes nonzero, p1 edges are curved and the chord test really would be "
    "an approximation -- but at 2e-16 it is not.", worst);
}

// The two predicates DO differ in general -- just not at p1, where the edge is
// straight. This pins that difference on a synthetic outline with no
// equilibrium in it, so the curve test's value is documented for the day p2
// geometry can be built. Until then tok_wall_curve_inside is correct, tested,
// and dormant: it cannot disagree with the chord test in any configuration
// this library will construct.
//
// The outline is a unit square with a V notch cut into its top edge:
//
//     (0,2)                       (2,2)
//       +-------+         +---------+
//       |        \       /          |
//       |         \     /           |     the V descends to (1.0,0.6);
//       |          \   /            |     everything inside the V is OUTSIDE
//       |           \ /             |     the polygon
//       |          (1,0.6)          |
//       +---------------------------+
//     (0,0)                       (2,0)
//
// a=(0.9,1.2) sits in the left prong, b=(1.1,1.2) in the right prong: both
// inside. At y=1.2 the notch spans x in (0.9143, 1.0857), so the straight chord
// from a to b crosses it. The quadratic through a, (1.0,0.4), b dives under the
// notch apex and stays in the body of the square.
static void
test_chord_and_curve_disagree_on_a_concave_notch(void)
{
  static const double R[] = { 0.0, 2.0, 2.0, 1.2, 1.0, 0.8, 0.0 };
  static const double Z[] = { 0.0, 0.0, 2.0, 2.0, 0.6, 2.0, 2.0 };

  struct gkyl_efit e = { 0 };
  e.limiter_status = 1;                       // usable outline
  e.limiter_n = sizeof(R)/sizeof(*R);
  e.limiter_R = (double *) R;
  e.limiter_Z = (double *) Z;
  e.rdim = 2.0; e.zdim = 2.0;                 // sets the on-edge tolerance
  TEST_ASSERT(gkyl_tok_wall_usable(&e));

  const double a[2]  = { 0.9, 1.2 };
  const double b[2]  = { 1.1, 1.2 };
  const double mid[2] = { 1.0, 0.4 };         // the represented edge's t=1/2 node

  // Both endpoints are inside; that is the whole difficulty -- an endpoint-only
  // test cannot tell these two cases apart.
  TEST_CHECK(tok_wall_point_inside(&e, a));
  TEST_MSG("endpoint a=(%g,%g) should be inside the left prong", a[0], a[1]);
  TEST_CHECK(tok_wall_point_inside(&e, b));
  TEST_MSG("endpoint b=(%g,%g) should be inside the right prong", b[0], b[1]);
  TEST_CHECK(tok_wall_point_inside(&e, mid));
  TEST_MSG("the curve's middle node (%g,%g) should be inside the body", mid[0], mid[1]);

  // The chord crosses the notch, so it must be refused.
  TEST_CHECK(!tok_wall_segment_inside(&e, a, b));
  TEST_MSG("the straight chord a->b crosses the notch and must NOT be called "
    "contained; if this passes, the fixture is wrong, not the library");

  // The actual represented edge does not, so it must be accepted.
  TEST_CHECK(tok_wall_curve_inside(&e, a, mid, b));
  TEST_MSG("the represented edge a->(1.0,0.4)->b stays inside the square and "
    "must be called contained");

  // Stated as the thing the gate is about: on this geometry the two predicates
  // give OPPOSITE answers, so which one the containment check calls decides
  // whether a contained domain is accepted.
  bool chord = tok_wall_segment_inside(&e, a, b);
  bool curve = tok_wall_curve_inside(&e, a, mid, b);
  TEST_CHECK(chord != curve);
  TEST_MSG("chord=%d curve=%d -- these must differ for the false refusal to be "
    "explicable by the chord approximation", chord, curve);
}

TEST_LIST = {
  { "refusal_reports_excursion_depth", test_refusal_reports_excursion_depth },
  { "refinement_separates_discretisation_from_declaration", test_refinement_separates_discretisation_from_declaration },
  { "curve_and_chord_are_identical_at_p1", test_curve_and_chord_are_identical_at_p1 },
  { "chord_and_curve_disagree_on_a_concave_notch", test_chord_and_curve_disagree_on_a_concave_notch },
  { NULL, NULL },
};
