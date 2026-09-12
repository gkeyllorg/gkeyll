#include <acutest.h>
#include <gkyl_efit.h>
#include <gkyl_tok_geo_priv.h>

#include <math.h>
#include <stdbool.h>
#include <stddef.h>

// These fixtures borrow stack arrays. They contain no equilibrium fields and
// must not be passed to gkyl_efit_release: wall predicates need only the outline.
static struct gkyl_efit
wall(int n, double *R, double *Z)
{
  return (struct gkyl_efit) {
    .limiter_status = 1, .limiter_n = n, .limiter_R = R, .limiter_Z = Z,
    .rdim = 4.0, .zdim = 4.0,
  };
}

static void
check_point(const struct gkyl_efit *efit, double R, double Z, bool expected)
{
  const double p[2] = { R, Z };
  TEST_CHECK(tok_wall_point_inside(efit, p) == expected);
  TEST_MSG("point (%g, %g): expected inside=%d", R, Z, expected);
}

static void
check_segment(const struct gkyl_efit *efit, const double a[2],
  const double b[2], bool expected)
{
  TEST_CHECK(tok_wall_segment_inside(efit, a, b) == expected);
  TEST_MSG("segment (%g, %g) -> (%g, %g): expected inside=%d",
    a[0], a[1], b[0], b[1], expected);
  TEST_CHECK(tok_wall_segment_inside(efit, b, a) == expected);
  TEST_MSG("reversed segment (%g, %g) -> (%g, %g): expected inside=%d",
    b[0], b[1], a[0], a[1], expected);
}

static void
test_square_points_and_winding(void)
{
  // No repeated closing vertex: the last-to-first edge must be included.
  double R[] = { 1.0, 3.0, 3.0, 1.0 };
  double Z[] = { -1.0, -1.0, 1.0, 1.0 };
  for (int reversed=0; reversed<2; ++reversed) {
    double r[4], z[4];
    for (int i=0; i<4; ++i) {
      int j = reversed ? 3-i : i;
      r[i] = R[j]; z[i] = Z[j];
    }
    struct gkyl_efit efit = wall(4, r, z);
    check_point(&efit, 2.0, 0.0, true);
    check_point(&efit, 1.0, 0.0, true);
    check_point(&efit, 3.0, 0.0, true);
    check_point(&efit, 2.0, -1.0, true);
    check_point(&efit, 2.0, 1.0, true);
    for (int i=0; i<4; ++i)
      check_point(&efit, R[i], Z[i], true);
    check_point(&efit, 0.9, 0.0, false);
    check_point(&efit, 3.1, 0.0, false);
    check_point(&efit, 2.0, -1.1, false);
    check_point(&efit, 2.0, 1.1, false);
  }
}

static void
test_square_segments(void)
{
  double R[] = { 1.0, 3.0, 3.0, 1.0 };
  double Z[] = { -1.0, -1.0, 1.0, 1.0 };
  struct gkyl_efit efit = wall(4, R, Z);
  const double center[] = { 2.0, 0.0 }, interior[] = { 2.5, 0.5 };
  const double lo[] = { 1.0, -1.0 }, hi[] = { 3.0, 1.0 };
  const double left_lo[] = { 1.0, -0.75 }, left_hi[] = { 1.0, 0.75 };
  const double bottom_lo[] = { 1.0, -1.0 }, bottom_hi[] = { 3.0, -1.0 };
  const double outside[] = { 3.5, 0.0 }, other_outside[] = { 0.5, 0.0 };
  const double extended_bottom[] = { 3.5, -1.0 };
  check_segment(&efit, center, interior, true);
  check_segment(&efit, lo, hi, true);
  check_segment(&efit, left_lo, left_hi, true); // closing boundary edge
  check_segment(&efit, bottom_lo, bottom_hi, true);
  check_segment(&efit, center, outside, false);
  check_segment(&efit, other_outside, outside, false); // crosses vessel
  check_segment(&efit, bottom_lo, extended_bottom, false);

  // Touching a polygon vertex alone cannot make an exterior segment valid.
  const double tangent_a[] = { 0.5, -0.5 }, tangent_b[] = { 1.5, -1.5 };
  check_segment(&efit, tangent_a, tangent_b, false);
}

static void
test_concave_segments(void)
{
  // A narrow notch removes 2<R<2.2, 1<Z<=3 from a rectangular vessel.
  double R[] = { 1.0, 4.0, 4.0, 2.2, 2.2, 2.0, 2.0, 1.0 };
  double Z[] = { 0.0, 0.0, 3.0, 3.0, 1.0, 1.0, 3.0, 3.0 };
  struct gkyl_efit efit = wall(8, R, Z);
  const double a[] = { 1.75, 2.0 }, b[] = { 2.45, 2.0 };
  check_point(&efit, a[0], a[1], true);
  check_point(&efit, b[0], b[1], true);
  check_point(&efit, 2.1, 2.0, false);
  check_segment(&efit, a, b, false);

  // Here even the midpoint is inside: checking only endpoints and midpoint
  // misses the short exterior interval near the first third of the segment.
  const double long_a[] = { 1.25, 2.0 }, long_b[] = { 3.75, 2.0 };
  check_point(&efit, long_a[0], long_a[1], true);
  check_point(&efit, long_b[0], long_b[1], true);
  check_point(&efit, 2.5, 2.0, true);
  check_segment(&efit, long_a, long_b, false);

  // Both endpoints lie on the notch boundary; the open segment is outside.
  const double notch_a[] = { 2.0, 2.0 }, notch_b[] = { 2.2, 2.0 };
  check_point(&efit, notch_a[0], notch_a[1], true);
  check_point(&efit, notch_b[0], notch_b[1], true);
  check_segment(&efit, notch_a, notch_b, false);

  // The notch floor is admissible boundary, and merely touching the reflex
  // vertex (2,1) while remaining inside on both sides is also admissible.
  const double floor_a[] = { 1.5, 1.0 }, floor_b[] = { 3.5, 1.0 };
  const double touch_a[] = { 1.5, 1.5 }, touch_b[] = { 2.5, 0.5 };
  check_segment(&efit, floor_a, floor_b, true);
  check_segment(&efit, touch_a, touch_b, true);

  // Reversing polygon winding must preserve a concave segment rejection.
  double reverse_R[8], reverse_Z[8];
  for (int i=0; i<8; ++i) {
    reverse_R[i] = R[7-i]; reverse_Z[i] = Z[7-i];
  }
  struct gkyl_efit reverse = wall(8, reverse_R, reverse_Z);
  check_segment(&reverse, long_a, long_b, false);
  check_segment(&reverse, touch_a, touch_b, true);
}

static void
test_collinear_boundary_vertices(void)
{
  // Extra collinear vertices and an explicitly repeated closing vertex must
  // not create a hole or make a boundary-contained segment leave the vessel.
  double R[] = { 1.0, 2.0, 3.0, 3.0, 3.0, 1.0, 1.0 };
  double Z[] = { -1.0, -1.0, -1.0, 0.0, 1.0, 1.0, -1.0 };
  struct gkyl_efit efit = wall(7, R, Z);
  const double a[] = { 1.25, -1.0 }, b[] = { 2.75, -1.0 };
  const double right_a[] = { 3.0, -0.75 }, right_b[] = { 3.0, 0.75 };
  check_point(&efit, 2.0, 0.0, true);
  check_point(&efit, 2.0, -1.0, true);
  check_point(&efit, 2.0, -1.1, false);
  check_segment(&efit, a, b, true);
  check_segment(&efit, right_a, right_b, true);
}

static void
test_zero_length_segments(void)
{
  double R[] = { 1.0, 3.0, 3.0, 1.0 };
  double Z[] = { -1.0, -1.0, 1.0, 1.0 };
  struct gkyl_efit efit = wall(4, R, Z);
  const double inside[] = { 2.0, 0.0 }, edge[] = { 1.0, 0.0 };
  const double vertex[] = { 1.0, -1.0 }, outside[] = { 0.5, 0.0 };
  check_segment(&efit, inside, inside, true);
  check_segment(&efit, edge, edge, true);
  check_segment(&efit, vertex, vertex, true);
  check_segment(&efit, outside, outside, false);
}

static void
test_unavailable_and_malformed_walls(void)
{
  double R[] = { 1.0, 3.0, 3.0, 1.0 };
  double Z[] = { -1.0, -1.0, 1.0, 1.0 };
  const double a[] = { 2.0, 0.0 }, b[] = { 2.5, 0.5 };
  for (int invalid=0; invalid<7; ++invalid) {
    struct gkyl_efit efit = wall(4, R, Z);
    switch (invalid) {
      case 0: efit.limiter_status = 0; break;
      case 1: efit.limiter_status = -1; break;
      case 2: efit.limiter_n = 0; break;
      case 3: efit.limiter_n = 1; break;
      case 4: efit.limiter_n = 2; break;
      case 5: efit.limiter_R = NULL; break;
      case 6: efit.limiter_Z = NULL; break;
    }
    check_point(&efit, a[0], a[1], false);
    check_segment(&efit, a, b, false);
    check_segment(&efit, a, a, false);
  }

  // A claimed-valid status must not allow nonfinite wall coordinates to
  // produce an apparent interior via partial ray crossings.
  struct gkyl_efit efit = wall(4, R, Z);
  R[0] = NAN;
  check_point(&efit, a[0], a[1], false);
  check_segment(&efit, a, b, false);
  R[0] = 1.0; Z[2] = INFINITY;
  check_point(&efit, a[0], a[1], false);
  check_segment(&efit, a, b, false);
}

static void
test_nonfinite_query_points(void)
{
  double R[] = { 1.0, 3.0, 3.0, 1.0 };
  double Z[] = { -1.0, -1.0, 1.0, 1.0 };
  struct gkyl_efit efit = wall(4, R, Z);
  const double good[] = { 2.0, 0.0 };
  const double bad[][2] = {
    { NAN, 0.0 }, { 2.0, NAN }, { INFINITY, 0.0 },
    { -INFINITY, 0.0 }, { 2.0, INFINITY }, { 2.0, -INFINITY },
  };
  for (size_t i=0; i<sizeof(bad)/sizeof(bad[0]); ++i) {
    check_point(&efit, bad[i][0], bad[i][1], false);
    check_segment(&efit, good, bad[i], false);
    check_segment(&efit, bad[i], bad[i], false);
  }
}

static void
test_quadratic_excursion_between_nodes(void)
{
  double R[]={1,4,4,1},Z[]={0,0,1,1};
  struct gkyl_efit e=wall(4,R,Z);
  // z(t)=1.0625-(t-.25)^2: all three nodes inside, but z(.25)>1.
  double a[]={1.25,1.0},m[]={2.5,1.0},b[]={3.75,.5};
  TEST_CHECK(tok_wall_point_inside(&e,a));
  TEST_CHECK(tok_wall_point_inside(&e,m));
  TEST_CHECK(tok_wall_point_inside(&e,b));
  TEST_CHECK(!tok_wall_curve_inside(&e,a,m,b));
  TEST_CHECK(!tok_wall_curve_inside(&e,b,m,a));
  // Move down to exact tangency, then further into the vessel.
  a[1]-=.0625;m[1]-=.0625;b[1]-=.0625;
  TEST_CHECK(tok_wall_curve_inside(&e,a,m,b));
  a[1]-=.1;m[1]-=.1;b[1]-=.1;
  TEST_CHECK(tok_wall_curve_inside(&e,a,m,b));
  // Straight edge is a degenerate quadratic and may lie on the wall.
  a[1]=m[1]=b[1]=0.0;
  TEST_CHECK(tok_wall_curve_inside(&e,a,m,b));
  // Collinear nodes in the vessel can still overshoot past a wall vertex.
  double c[]={4.0,.5},n[]={4.0,.5},d[]={3.5,.5};
  TEST_CHECK(!tok_wall_curve_inside(&e,c,n,d));
}

// The vessel-outline declaration is a two-way contract, so test both directions:
// it must permit a genuinely absent outline and must refuse every other case,
// including a declaration that contradicts a wall that is actually present.
static void
test_vessel_outline_declaration(void)
{
  double R[]={0.0,1.0,1.0,0.0}, Z[]={0.0,0.0,1.0,1.0};
  struct gkyl_efit usable=wall(4,R,Z);            // status 1, 4 vertices
  struct gkyl_efit absent=wall(4,R,Z); absent.limiter_status=0; absent.limiter_n=0;
  struct gkyl_efit degenerate=wall(4,R,Z); degenerate.limiter_status=2; degenerate.limiter_n=1;
  struct gkyl_efit malformed=wall(4,R,Z); malformed.limiter_status=-1; malformed.limiter_n=0;

  struct gkyl_tok_geo_grid_inp silent={0};
  struct gkyl_tok_geo_grid_inp declared={0}; declared.no_vessel_outline=true;

  // A usable outline is enforced, declaration or not -- and declaring absence
  // when a wall exists is a contradiction, never a way to switch enforcement off.
  TEST_CHECK(gkyl_tok_wall_policy_for(&silent,&usable)==GKYL_TOK_WALL_ENFORCE);
  TEST_CHECK(gkyl_tok_wall_policy_for(&declared,&usable)==GKYL_TOK_WALL_REJECT_CONTRADICTED);

  // Silence never disables the constraint.
  TEST_CHECK(gkyl_tok_wall_policy_for(&silent,&absent)==GKYL_TOK_WALL_REJECT_UNDECLARED);
  TEST_CHECK(gkyl_tok_wall_policy_for(&silent,&degenerate)==GKYL_TOK_WALL_REJECT_UNDECLARED);

  // Absent and degenerate outlines are benign and declarable.
  TEST_CHECK(gkyl_tok_wall_policy_for(&declared,&absent)==GKYL_TOK_WALL_NOT_ENFORCED);
  TEST_CHECK(gkyl_tok_wall_policy_for(&declared,&degenerate)==GKYL_TOK_WALL_NOT_ENFORCED);

  // Unreadable data is an input error and stays un-declarable either way.
  TEST_CHECK(gkyl_tok_wall_policy_for(&silent,&malformed)==GKYL_TOK_WALL_REJECT_UNDECLARED);
  TEST_CHECK(gkyl_tok_wall_policy_for(&declared,&malformed)==GKYL_TOK_WALL_REJECT_MALFORMED);

  // The usability predicate itself, which every site now shares.
  TEST_CHECK(gkyl_tok_wall_usable(&usable));
  TEST_CHECK(!gkyl_tok_wall_usable(&absent));
  TEST_CHECK(!gkyl_tok_wall_usable(&degenerate));
  TEST_CHECK(!gkyl_tok_wall_usable(&malformed));
  TEST_CHECK(!gkyl_tok_wall_usable(NULL));

  // A 2-vertex outline cannot bound a region even though it parses.
  struct gkyl_efit two=wall(4,R,Z); two.limiter_n=2;
  TEST_CHECK(!gkyl_tok_wall_usable(&two));
}

TEST_LIST = {
  { "square_points_and_winding", test_square_points_and_winding },
  { "square_segments", test_square_segments },
  { "concave_segments", test_concave_segments },
  { "collinear_boundary_vertices", test_collinear_boundary_vertices },
  { "zero_length_segments", test_zero_length_segments },
  { "unavailable_and_malformed_walls", test_unavailable_and_malformed_walls },
  { "nonfinite_query_points", test_nonfinite_query_points },
  { "quadratic_excursion_between_nodes", test_quadratic_excursion_between_nodes },
  { "vessel_outline_declaration", test_vessel_outline_declaration },
  { NULL, NULL },
};
