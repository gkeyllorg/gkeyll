// Gate 4: a disconnected plate must not manufacture a strike point.
//
// `plate_func` is not required to be continuous, and two separate divertor
// tiles are an ordinary thing to describe -- "a disconnected plate" is one of
// the fixtures this gate asks for. Where two disjoint pieces sit on opposite
// sides of the requested flux surface, f(s) = psi(plate(s)) - psi0 changes sign
// across the gap with no root anywhere between them. A bracket-and-bisect
// search cannot tell that from a real crossing: every value it evaluates while
// narrowing stays on one piece, so it converges onto the gap and never notices.
//
// Measured before the fix: neither piece touched the surface on its own, yet
// the pair reported success, for gaps from 2% to 40% of the plate. At the
// widest the accepted strike point sat about 55 mm from the surface. That is a
// FALSE ACCEPTANCE -- every other failure mode in this routine refuses, which
// is why this one mattered.
//
// The fix confirms the located point is on the requested surface before
// returning it. This test pins both halves: the gap is refused, and a plate
// that genuinely does reach the surface is still accepted.

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include <acutest.h>
#include <gkyl_tok_geo.h>

static void
base_plate(double t, double *RZ)
{
  RZ[0] = 1.5966 + (1.6888 - 1.5966)*t;
  RZ[1] = -1.1421 + (-0.8781 - (-1.1421))*t;
}

static double g_trunc = 1.0;
static void plate_truncated(double s, double *RZ) { base_plate(g_trunc*s, RZ); }

static double g_lo_end = 0.0, g_hi_start = 1.0;
static void plate_piece_a(double s, double *RZ) { base_plate(g_lo_end*s, RZ); }
static void plate_piece_b(double s, double *RZ)
{ base_plate(g_hi_start + (1.0-g_hi_start)*s, RZ); }

// Discontinuous at s = 0.5, which is sample 256 of 512, so consecutive samples
// straddle the jump and the spurious sign change is unmissable.
static void
plate_disconnected(double s, double *RZ)
{
  if (s < 0.5) base_plate(g_lo_end*(2.0*s), RZ);
  else         base_plate(g_hi_start + (1.0-g_hi_start)*(2.0*s-1.0), RZ);
}

static const double psi0 = 0.1600;   // inside ASDEX's covered band

static bool
covered(plate_func p)
{
  struct gkyl_efit_inp einp = { .rz_poly_order=2, .flux_poly_order=1 };
  snprintf(einp.filepath, sizeof(einp.filepath), "%s",
    "gyrokinetic/data/eqdsk/asdex.geqdsk");
  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype=GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO,
    .rmin=0.0, .rmax=5.0, .rclose=2.5, .rright=2.5, .rleft=0.7,
    .zmin=-1.3, .zmax=1.0, .zmin_left=-1.0, .zmin_right=-0.9,
    .plate_spec=true, .plate_func_lower=p, .plate_func_upper=p };
  struct gkyl_tok_geo *geo = gkyl_tok_geo_new(&einp,&ginp);
  if (!geo) return false;
  struct gkyl_tok_geo_grid_inp q = { .ftype=GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO };
  bool ok = gkyl_tok_geo_check_plate_coverage(geo,&q,psi0);
  gkyl_tok_geo_release(geo);
  return ok;
}

// Locate the crossing using the library as its own oracle: coverage of a plate
// truncated at L holds exactly while L >= t*.
static double
find_crossing_parameter(void)
{
  g_trunc = 1.0;
  TEST_ASSERT(covered(plate_truncated));
  double lo = 0.0, hi = 1.0;
  for (int k=0; k<60; ++k) {
    double mid = 0.5*(lo+hi);
    g_trunc = mid;
    if (covered(plate_truncated)) hi = mid; else lo = mid;
  }
  g_trunc = 1.0;
  return hi;
}

static void
test_gap_is_not_mistaken_for_a_crossing(void)
{
  const double tstar = find_crossing_parameter();
  const double deltas[] = { 0.01, 0.05, 0.10, 0.20 };

  for (size_t i=0; i<sizeof(deltas)/sizeof(*deltas); ++i) {
    double d = deltas[i];
    g_lo_end = tstar - d;
    g_hi_start = tstar + d;
    if (g_lo_end <= 0.0 || g_hi_start >= 1.0) continue;

    // Control: neither piece reaches the surface by itself. Without this, the
    // pair being refused would prove nothing.
    TEST_ASSERT(!covered(plate_piece_a));
    TEST_ASSERT(!covered(plate_piece_b));

    TEST_CHECK(!covered(plate_disconnected));
    TEST_MSG("a plate split by a gap of %.0f%% of its length (t in [%.4f,%.4f], "
      "crossing at %.4f) was accepted. Neither piece touches the surface, so "
      "the strike point returned does not exist.",
      200.0*d, g_lo_end, g_hi_start, tstar);
  }
  g_lo_end = 0.0; g_hi_start = 1.0;
}

// The other half: the fix must not have turned into "refuse everything". A
// continuous plate that genuinely reaches the surface is still accepted, and
// one truncated short of it is still refused.
static void
test_a_real_crossing_is_still_accepted(void)
{
  const double tstar = find_crossing_parameter();

  g_trunc = 1.0;
  TEST_CHECK(covered(plate_truncated));
  TEST_MSG("the full plate no longer reaches the surface -- the residual check "
    "is rejecting genuine roots");

  g_trunc = tstar*1.01;
  TEST_CHECK(covered(plate_truncated));
  TEST_MSG("a plate reaching just past the surface was refused");

  g_trunc = tstar*0.99;
  TEST_CHECK(!covered(plate_truncated));
  TEST_MSG("a plate stopping short of the surface was accepted");

  g_trunc = 1.0;
}

TEST_LIST = {
  { "gap_is_not_mistaken_for_a_crossing", test_gap_is_not_mistaken_for_a_crossing },
  { "a_real_crossing_is_still_accepted", test_a_real_crossing_is_still_accepted },
  { NULL, NULL },
};
