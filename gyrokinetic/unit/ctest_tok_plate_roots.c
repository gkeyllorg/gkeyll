// Gate 4: the resolution floor of plate strike-point finding.
//
// `tok_plate_intersection` samples the plate at 512 uniform points in s and
// brackets sign changes of f(s) = psi(plate(s)) - psi0. That is a fixed grid,
// so it has a feature-size floor, and this pins where the floor is.
//
// Measured (gate4/probe_plate_roots.c), with the sample spacing 1/512:
//
//   root pair half-separation > 0.5/512  -> found      (2.0e-3 … 1.0e-1 tested)
//   root pair half-separation < 0.5/512  -> MISSED     (5.0e-4 … 5.0e-5 tested)
//   tangent root, no sign change         -> MISSED
//
// The transition sits exactly at the predicted 9.77e-4 boundary.
//
// Which way the failure goes, and why that matters
// ------------------------------------------------
// A missed root makes `tok_plate_intersection` return false, which
// `tok_plate_coverage_status` turns into TOK_GEO_ACTIVE_PLATE_FAILED and a
// rejection. So the failure is CONSERVATIVE: the library refuses to build
// rather than silently accepting geometry whose strike point it could not
// locate. This is a robustness limit, not a correctness hole -- and it is why
// these tests assert the floor rather than treating the misses as defects.
//
// What is asserted here is the half of the contract that must not regress:
// features COARSER than the sampling are found. A change that reduced the
// sample count, or perturbed the bracketing, would break these. The misses
// below the floor are recorded in gate4/RECON.md as a documented limitation,
// deliberately not frozen into a test -- pinning them would make the floor
// harder to improve later.

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include <acutest.h>
#include <gkyl_tok_geo.h>

// The straight ASDEX outboard plate, from the commented-out form in
// creg/rt_gk_multib_asdex_2x2v_p1.c.
static void
base_plate(double t, double *RZ)
{
  RZ[0] = 1.5966 + (1.6888 - 1.5966)*t;
  RZ[1] = -1.1421 + (-0.8781 - (-1.1421))*t;
}

static double g_trunc = 1.0;
static void
plate_truncated(double s, double *RZ) { base_plate(g_trunc*s, RZ); }

// There-and-back, turning around at s0. s0 is placed MID-CELL so that no
// sample lands on the apex: an earlier version of this fixture put the apex at
// s = 0.5, which is sample 256 of 512, and the explicit |f| <= flux_tol
// endpoint test then "found" everything. The fixture, not the library, was
// being measured.
static double g_apex = 1.0;
static double g_s0 = 0.5;
static void
plate_there_and_back(double s, double *RZ)
{
  double M = g_s0 > 0.5 ? g_s0 : 1.0 - g_s0;
  double u = 1.0 - fabs(s - g_s0)/M;
  if (u < 0.0) u = 0.0;
  base_plate(g_apex*u, RZ);
}

// A plate that STARTS exactly on the surface: t runs from tstar to 1.
static double g_start = 0.0;
static void
plate_from_surface(double s, double *RZ)
{ base_plate(g_start + (1.0 - g_start)*s, RZ); }

static const double psi0 = 0.1600;   // inside ASDEX's covered band [0.1498, 0.16973]
static const int nsamp = 512;        // the sampling in tok_plate_intersection

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

// Find where the straight plate crosses the psi0 surface, using the library as
// its own oracle: truncate the plate at L, and coverage holds exactly while
// L >= t*. This is the gate plan's "shorten the callback" operation, and it
// avoids reimplementing psi(R,Z) -- which is static, and which a test copy of
// could silently drift from.
static double
find_crossing_parameter(void)
{
  g_trunc = 1.0;
  TEST_ASSERT(covered(plate_truncated));   // the full plate must reach it
  double lo = 0.0, hi = 1.0;
  for (int k=0; k<60; ++k) {
    double mid = 0.5*(lo+hi);
    g_trunc = mid;
    if (covered(plate_truncated)) hi = mid; else lo = mid;
  }
  g_trunc = 1.0;
  return hi;
}

// A plate that reaches the surface must be reported as reaching it, and one
// truncated just short of it must not. Without this the tests below could pass
// on a library that simply answered "no" to everything.
static void
test_truncation_boundary_is_sharp(void)
{
  double tstar = find_crossing_parameter();
  TEST_CHECK(tstar > 0.0 && tstar < 1.0);
  TEST_MSG("crossing parameter %.17g is not strictly inside the plate", tstar);

  g_trunc = tstar*1.001;
  TEST_CHECK(covered(plate_truncated));
  TEST_MSG("a plate reaching just PAST the surface was reported unreachable");

  g_trunc = tstar*0.999;
  TEST_CHECK(!covered(plate_truncated));
  TEST_MSG("a plate stopping just SHORT of the surface was reported reachable");
  g_trunc = 1.0;
}

// The contract that must not regress: a plate whose crossings are separated by
// more than the sample spacing is found. These are the separations a real
// shaped plate produces.
static void
test_features_above_the_sampling_floor_are_found(void)
{
  double tstar = find_crossing_parameter();
  g_s0 = (nsamp/2 + 0.5)/(double) nsamp;          // mid-cell: no sample on the apex
  double M = g_s0 > 0.5 ? g_s0 : 1.0 - g_s0;
  const double floor_half_sep = 0.5/(double) nsamp;

  const double ds[] = { 0.2, 0.05, 0.01, 0.004, 0.002 };
  for (size_t i=0; i<sizeof(ds)/sizeof(*ds); ++i) {
    double d = ds[i];
    g_apex = tstar/(1.0 - d);
    if (g_apex > 1.0) continue;                   // apex would run off the plate
    double half_sep = d*M;
    TEST_ASSERT(half_sep > floor_half_sep);       // this case must be above the floor
    TEST_CHECK(covered(plate_there_and_back));
    TEST_MSG("two crossings separated by %.3e (floor %.3e) were not found; the "
      "plate strike-point search has lost resolution it previously had",
      2*half_sep, 2*floor_half_sep);
  }
  g_apex = 1.0; g_s0 = 0.5;
}

// The floor is a property of the sampling, so state it where a plate author
// can act on it: features finer than this are not resolved, and the library
// responds by REFUSING rather than by guessing.
static void
test_the_floor_is_where_the_sampling_puts_it(void)
{
  double tstar = find_crossing_parameter();
  g_s0 = (nsamp/2 + 0.5)/(double) nsamp;
  double M = g_s0 > 0.5 ? g_s0 : 1.0 - g_s0;
  const double floor_half_sep = 0.5/(double) nsamp;

  // Comfortably above the floor: found.
  double d_above = 4.0*floor_half_sep/M;
  g_apex = tstar/(1.0 - d_above);
  bool above = g_apex <= 1.0 ? covered(plate_there_and_back) : true;
  TEST_CHECK(above);
  TEST_MSG("a pair at 4x the sampling floor was missed; the floor has moved");

  // Comfortably below it: not found, and that shows up as a REFUSAL, which is
  // the safe direction. Recorded so the asymmetry is visible in the suite.
  double d_below = 0.1*floor_half_sep/M;
  g_apex = tstar/(1.0 - d_below);
  bool below = g_apex <= 1.0 ? covered(plate_there_and_back) : false;
  TEST_CHECK(!below);
  TEST_MSG("a pair at a tenth of the sampling floor was found -- if the search "
    "has genuinely improved, raise the floor in the plate_func documentation "
    "and in gate4/RECON.md rather than leaving them disagreeing");

  g_apex = 1.0; g_s0 = 0.5;
}

// A root exactly at s = 0. The code handles this before the sampling loop
// begins, with an explicit |f| <= flux_tol test on plate(0). Worth asserting
// directly: it is a separate branch from everything the loop does, and the
// residual check added for disconnected plates now also runs on this path.
static void
test_root_exactly_at_s_equals_zero(void)
{
  const double tstar = find_crossing_parameter();
  g_start = tstar;
  TEST_CHECK(covered(plate_from_surface));
  TEST_MSG("a plate whose s=0 endpoint lies exactly on the surface was reported "
    "unreachable; the pre-loop endpoint test no longer fires");
  g_start = 0.0;
}

// A root exactly at s = 1, caught by the loop's endpoint test on the final
// sample. The companion to the case above, and the exact boundary that
// test_truncation_boundary_is_sharp only brackets from either side.
static void
test_root_exactly_at_s_equals_one(void)
{
  const double tstar = find_crossing_parameter();
  g_trunc = tstar;
  TEST_CHECK(covered(plate_truncated));
  TEST_MSG("a plate whose s=1 endpoint lies exactly on the surface was reported "
    "unreachable; the loop's final-sample endpoint test no longer fires");
  g_trunc = 1.0;
}

TEST_LIST = {
  { "root_exactly_at_s_equals_zero", test_root_exactly_at_s_equals_zero },
  { "root_exactly_at_s_equals_one", test_root_exactly_at_s_equals_one },
  { "truncation_boundary_is_sharp", test_truncation_boundary_is_sharp },
  { "features_above_the_sampling_floor_are_found", test_features_above_the_sampling_floor_are_found },
  { "the_floor_is_where_the_sampling_puts_it", test_the_floor_is_where_the_sampling_puts_it },
  { NULL, NULL },
};
