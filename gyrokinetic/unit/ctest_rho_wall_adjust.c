#include <acutest.h>
#include <gkyl_rho_wall_adjust_priv.h>

#include <float.h>
#include <math.h>
#include <stddef.h>

static void
reject_unchanged(double requested, double other, double axis, double sep,
  int family, int step)
{
  double rho = 123.25, psi = -456.75;
  TEST_CHECK(!gkyl_rho_wall_next(requested, other, axis, sep,
    family, step, &rho, &psi));
  TEST_MSG("expected rejection: requested=%g other=%g axis=%g sep=%g family=%d step=%d",
    requested, other, axis, sep, family, step);
  TEST_CHECK(rho == 123.25 && psi == -456.75);
  TEST_MSG("failure modified outputs: rho=%g psi=%g", rho, psi);
}

static void
test_type_ownership(void)
{
  // The enum's complete declared order is checked, including excluded CORE,
  // IWL and unimplemented USN. A new enum requires an explicit policy decision.
  const int families[] = {
    1, 1, 1, 0, 0,  // standalone DN/LSN, USN, CORE
    1, 1, 1, 1, 1, 1,  // six DN sections
    1, 1, 1,  // three LSN sections
    2, 2, 2, 2,  // upper/lower PF halves
    0, 0, 0,  // split cores and IWL
  };
  TEST_CHECK(sizeof families/sizeof families[0] == GKYL_GEOMETRY_TOKAMAK_IWL+1);
  for (size_t i=0; i<sizeof families/sizeof families[0]; ++i) {
    TEST_CHECK(gkyl_rho_wall_family((enum gkyl_tok_geo_type) i) == families[i]);
    TEST_MSG("wrong adjustment family for enum %zu", i);
  }
  TEST_CHECK(gkyl_rho_wall_family((enum gkyl_tok_geo_type) -1) == 0);
  TEST_CHECK(gkyl_rho_wall_family((enum gkyl_tok_geo_type) 999) == 0);
  reject_unchanged(1.05, 1.0, 0.0, 1.0,
    gkyl_rho_wall_family(GKYL_GEOMETRY_TOKAMAK_CORE), 1);
}

static void
test_sol_direction_and_flux_sign(void)
{
  double rho, psi;
  // Expectations are derived from GKYL_RHO_WALL_STEP, never written as a
  // literal: a test that hardcodes the increment fails for the right behaviour
  // the moment the increment is retuned, which is exactly what happened here.
  const double S = GKYL_RHO_WALL_STEP;
  TEST_ASSERT(gkyl_rho_wall_next(1.051, 1.0, -2.0, 2.0, 1, 1, &rho, &psi));
  TEST_CHECK(fabs(rho-(1.051-S)) < 8*DBL_EPSILON);
  TEST_CHECK(fabs(psi-((1.051-S)*(1.051-S)*4.0-2.0)) < 32*DBL_EPSILON);
  double positive = psi;
  TEST_ASSERT(gkyl_rho_wall_next(1.051, 1.0, 2.0, -2.0, 1, 1, &rho, &psi));
  TEST_CHECK(fabs(psi+positive) < 32*DBL_EPSILON);
  TEST_CHECK(psi < -2.0); // still on the SOL side of native separatrix

  // A nonseparatrix fixed boundary also limits shrinkage.
  const int collapse = (int) lround((1.25-1.125)/GKYL_RHO_WALL_STEP);
  TEST_ASSERT(gkyl_rho_wall_next(1.25, 1.125, 0.0, 1.0, 1, collapse-1, &rho, &psi));
  TEST_CHECK(rho > 1.125 && rho < 1.25);
  reject_unchanged(1.25, 1.125, 0.0, 1.0, 1, collapse);   // exact collapse
  reject_unchanged(1.25, 1.125, 0.0, 1.0, 1, collapse+1); // crossing
  reject_unchanged(1.125, 0.9, 0.0, 1.0, 1,
    (int) lround((1.125-1.0)/GKYL_RHO_WALL_STEP));        // separatrix stop
  reject_unchanged(1.0, 0.9, 0.0, 1.0, 1, 1);
  reject_unchanged(0.9, 0.8, 0.0, 1.0, 1, 1);
}

static void
test_pf_direction_and_flux_sign(void)
{
  double rho, psi;
  const double S = GKYL_RHO_WALL_STEP;
  TEST_ASSERT(gkyl_rho_wall_next(0.899, 1.0, -2.0, 2.0, 2, 1, &rho, &psi));
  TEST_CHECK(fabs(rho-(0.899+S)) < 8*DBL_EPSILON);
  TEST_CHECK(fabs(psi-((0.899+S)*(0.899+S)*4.0-2.0)) < 32*DBL_EPSILON);
  double positive = psi;
  TEST_ASSERT(gkyl_rho_wall_next(0.899, 1.0, 2.0, -2.0, 2, 1, &rho, &psi));
  TEST_CHECK(fabs(psi+positive) < 32*DBL_EPSILON);
  TEST_CHECK(psi > -2.0); // PF remains below rho=1 for either flux sign
  TEST_ASSERT(gkyl_rho_wall_next(0.0, 1.0, 0.0, 1.0, 2, 1, &rho, &psi));
  TEST_CHECK(rho == GKYL_RHO_WALL_STEP);

  const int pf_collapse = (int) lround((0.875-0.75)/GKYL_RHO_WALL_STEP);
  TEST_ASSERT(gkyl_rho_wall_next(0.75, 0.875, 0.0, 1.0, 2, pf_collapse-1, &rho, &psi));
  TEST_CHECK(rho < 0.875 && rho > 0.75);
  reject_unchanged(0.75, 0.875, 0.0, 1.0, 2, pf_collapse);
  reject_unchanged(0.75, 0.875, 0.0, 1.0, 2, pf_collapse+1);
  reject_unchanged(0.875, 1.2, 0.0, 1.0, 2,
    (int) lround((1.0-0.875)/GKYL_RHO_WALL_STEP)); // separatrix stop
  reject_unchanged(-0.1, 1.0, 0.0, 1.0, 2, 1);
  reject_unchanged(1.0, 1.1, 0.0, 1.0, 2, 1);
  reject_unchanged(1.05, 1.1, 0.0, 1.0, 2, 1);
}

static void
test_original_request_steps(void)
{
  // Preserve the sub-increment part of the request and avoid accumulated
  // subtraction error. Query order is deliberately nonmonotone.
  const double requested = 3.14159265;
  const int steps[] = { 2000, 1, 1700, 17, 2000 };
  double rho, psi, first = 0.0;
  for (size_t i=0; i<sizeof steps/sizeof steps[0]; ++i) {
    TEST_ASSERT(gkyl_rho_wall_next(requested, 1.0, 0.0, 1.0,
      1, steps[i], &rho, &psi));
    TEST_CHECK(fabs(rho-(requested-steps[i]*GKYL_RHO_WALL_STEP)) < 4*DBL_EPSILON);
    TEST_CHECK(fabs(sqrt(psi)-rho) < 8*DBL_EPSILON);
    if (i == 0) first = rho;
    if (i == 4) TEST_CHECK(rho == first);
  }
  TEST_ASSERT(gkyl_rho_wall_next(0.10025, 1.0, 0.0, 1.0,
    2, 899, &rho, &psi));
  TEST_CHECK(fabs(rho-(0.10025+899*GKYL_RHO_WALL_STEP)) < 4*DBL_EPSILON);

  // Full reach: MAX_STEPS*STEP, which must still land above the rho=1 stop.
  const double reach = GKYL_RHO_WALL_MAX_STEPS*GKYL_RHO_WALL_STEP;
  TEST_ASSERT(gkyl_rho_wall_next(12.0, 1.0, 0.0, 1.0,
    1, GKYL_RHO_WALL_MAX_STEPS, &rho, &psi));
  TEST_CHECK(fabs(rho-(12.0-reach)) < 8*DBL_EPSILON);
  TEST_CHECK(fabs(psi-(12.0-reach)*(12.0-reach)) < 64*DBL_EPSILON);
  reject_unchanged(12.0, 1.0, 0.0, 1.0, 1, GKYL_RHO_WALL_MAX_STEPS+1);
  reject_unchanged(1.05, 1.0, 0.0, 1.0, 1, 0);
  reject_unchanged(1.05, 1.0, 0.0, 1.0, 1, -1);
}

static void
test_invalid_and_nonfinite_inputs(void)
{
  const double bad[] = { NAN, INFINITY, -INFINITY };
  for (size_t i=0; i<sizeof bad/sizeof bad[0]; ++i) {
    TEST_CHECK(!gkyl_rho_wall_finite(bad[i]));
    reject_unchanged(bad[i], 1.0, 0.0, 1.0, 1, 1);
    reject_unchanged(1.05, bad[i], 0.0, 1.0, 1, 1);
    reject_unchanged(1.05, 1.0, bad[i], 1.0, 1, 1);
    reject_unchanged(1.05, 1.0, 0.0, bad[i], 1, 1);
    reject_unchanged(bad[i], 1.0, 0.0, 1.0, 2, 1);
  }
  TEST_CHECK(gkyl_rho_wall_finite(DBL_MAX));
  TEST_CHECK(gkyl_rho_wall_finite(-DBL_MAX));
  TEST_CHECK(gkyl_rho_wall_finite(0.0));
  reject_unchanged(1.05, 1.0, 2.0, 2.0, 1, 1);
  reject_unchanged(1.05, -0.5, 0.0, 1.0, 1, 1);
  reject_unchanged(1.05, 1.0, 0.0, 1.0, 3, 1);
  reject_unchanged(1.05, 1.0, 0.0, 1.0, -1, 1);
  reject_unchanged(2.0, 1.0, 0.0, DBL_MAX, 1, 1); // output flux overflow
  reject_unchanged(1.05, 1.0, DBL_MAX, -DBL_MAX, 1, 1); // span overflow
  reject_unchanged(DBL_MAX, 1.0, 0.0, 1.0, 1, 1); // no representable progress

  double value = 123.25;
  TEST_CHECK(!gkyl_rho_wall_next(1.05, 1.0, 0.0, 1.0, 1, 1, NULL, &value));
  TEST_CHECK(value == 123.25);
  TEST_CHECK(!gkyl_rho_wall_next(1.05, 1.0, 0.0, 1.0, 1, 1, &value, NULL));
  TEST_CHECK(value == 123.25);
  TEST_CHECK(!gkyl_rho_wall_next(1.05, 1.0, 0.0, 1.0, 1, 1, &value, &value));
  TEST_CHECK(value == 123.25);
}

// The two-phase search: a coarse approach that is only ever an accelerator,
// and a fine polish that must reproduce the single-phase answer exactly.
static void
test_phased_steps(void)
{
  double rho, psi, rho_ref, psi_ref;
  const double C = GKYL_RHO_WALL_COARSE_STEP, S = GKYL_RHO_WALL_STEP;

  // The coarse step must be an EXACT multiple of the fine one, or a polished
  // answer could land between representable fine points and the precision the
  // fine step buys would be lost to the lattice mismatch.
  TEST_CHECK(lround(C/S) == 10);
  TEST_CHECK(fabs(C - 10.0*S) < 1e-18);

  // One coarse step and ten fine steps are the SAME point, reached two ways.
  TEST_ASSERT(gkyl_rho_wall_next_phased(1.25, 1.0, 0.0, 1.0, 1, 1, 0, &rho, &psi));
  TEST_ASSERT(gkyl_rho_wall_next_phased(1.25, 1.0, 0.0, 1.0, 1, 0, 10, &rho_ref, &psi_ref));
  TEST_CHECK(fabs(rho-rho_ref) < 1e-15);
  TEST_CHECK(fabs(psi-psi_ref) < 1e-15);

  // Composition: the delta is coarse*C + fine*S, on both families.
  TEST_ASSERT(gkyl_rho_wall_next_phased(1.25, 1.0, 0.0, 1.0, 1, 3, 7, &rho, &psi));
  TEST_CHECK(fabs(rho - (1.25 - (3*C + 7*S))) < 1e-15);
  TEST_ASSERT(gkyl_rho_wall_next_phased(0.75, 1.0, 0.0, 1.0, 2, 3, 7, &rho, &psi));
  TEST_CHECK(fabs(rho - (0.75 + (3*C + 7*S))) < 1e-15);

  // The single-phase wrapper is the zero-coarse case, exactly.
  TEST_ASSERT(gkyl_rho_wall_next(1.25, 1.0, 0.0, 1.0, 1, 42, &rho_ref, &psi_ref));
  TEST_ASSERT(gkyl_rho_wall_next_phased(1.25, 1.0, 0.0, 1.0, 1, 0, 42, &rho, &psi));
  TEST_CHECK(rho == rho_ref && psi == psi_ref);

  // Rewinding a coarse step and polishing must never cross the bracket it came
  // from: coarse k-1 plus ten fine steps is exactly coarse k.
  TEST_ASSERT(gkyl_rho_wall_next_phased(1.25, 1.0, 0.0, 1.0, 1, 4, 0, &rho_ref, &psi_ref));
  TEST_ASSERT(gkyl_rho_wall_next_phased(1.25, 1.0, 0.0, 1.0, 1, 3, 10, &rho, &psi));
  TEST_CHECK(fabs(rho-rho_ref) < 1e-15);

  // Malformed counts are refused, and the outputs are untouched on refusal.
  rho = psi = -12345.0;
  TEST_CHECK(!gkyl_rho_wall_next_phased(1.25, 1.0, 0.0, 1.0, 1, -1, 5, &rho, &psi));
  TEST_CHECK(!gkyl_rho_wall_next_phased(1.25, 1.0, 0.0, 1.0, 1, 5, -1, &rho, &psi));
  TEST_CHECK(!gkyl_rho_wall_next_phased(1.25, 1.0, 0.0, 1.0, 1, 0, 0, &rho, &psi));
  TEST_CHECK(rho == -12345.0 && psi == -12345.0);

  // The total reach is still bounded in FINE units, so a coarse count cannot
  // smuggle the search past MAX_STEPS.
  TEST_CHECK(!gkyl_rho_wall_next_phased(2.0, 1.0, 0.0, 1.0, 1,
    GKYL_RHO_WALL_MAX_STEPS/10 + 1, 0, &rho, &psi));
}

TEST_LIST = {
  { "phased_steps", test_phased_steps },
  { "type_ownership", test_type_ownership },
  { "sol_direction_and_flux_sign", test_sol_direction_and_flux_sign },
  { "pf_direction_and_flux_sign", test_pf_direction_and_flux_sign },
  { "original_request_steps", test_original_request_steps },
  { "invalid_and_nonfinite_inputs", test_invalid_and_nonfinite_inputs },
  { NULL, NULL },
};
